#' Build LD/X_ref arguments for colocboost based on data type.
#'
#' When LD matrices are genotype X (non-square, rows=samples, cols=variants),
#' passes them as X_ref to colocboost. Otherwise passes as LD (correlation).
#'
#' @param ld_list A list of matrices (correlation R or genotype X).
#' @param subset Optional index vector to subset ld_list (e.g., from dict_sumstatLD).
#' @return A named list with either `LD = ...` or `X_ref = ...`.
#' @noRd
build_ld_args <- function(ld_list, subset = NULL) {
  if (!is.null(subset)) ld_list <- ld_list[subset]
  # Detect: if any matrix is non-square, it's genotype X (samples x variants)
  is_geno <- any(sapply(ld_list, function(m) nrow(m) != ncol(m)))
  if (is_geno) list(X_ref = ld_list) else list(LD = ld_list)
}

# Run colocboost with tryCatch and timing.
# @noRd
.run_colocboost <- function(label, ...) {
  t1 <- Sys.time()
  res <- tryCatch(
    colocboost(...),
    error = function(e) {
      message(label, " failed: ", conditionMessage(e))
      NULL
    }
  )
  list(result = res, time = Sys.time() - t1)
}

#' Multi-trait colocalization analysis pipeline
#'
#' This function perform a multi-trait colocalization using ColocBoost
#'
#' @param region_data A region data loaded from \code{load_regional_data}.
#' @param focal_trait Name of trait if perform focaled ColocBoost
#' @param event_filters A list of pattern for filtering events based on context names. Example: for sQTL, list(type_pattern = ".*clu_(\\d+_[+-?]).*",valid_pattern = "clu_(\\d+_[+-?]):PR:",exclude_pattern = "clu_(\\d+_[+-?]):IN:")
#' @param maf_cutoff A scalar to remove variants with maf < maf_cutoff, dafault is 0.005.
#' @param pip_cutoff_to_skip_ind A vector of cutoff values for skipping analysis based on PIP values for each context. Default is 0.
#' @param pip_cutoff_to_skip_sumstat A vector of cutoff values for skipping analysis based on PIP values for each sumstat Default is 0.
#' @param qc_method Quality control method to use. Options are "slalom" or "dentist" (default: "slalom").
#' @param impute Logical; if TRUE, performs imputation for outliers identified in the analysis (default: TRUE).
#' @param impute_opts A list of imputation options including rcond, R2_threshold, and minimum_ld (default: list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5)).
#'
#'
#' @return A list containing the individual_data and sumstat_data after QC:
#' individual_data contains the following components if exist
#' \itemize{
#'   \item Y: A list of residualized phenotype values for all tasks.
#'   \item X: A list of residualized genotype matrices all tasks.
#' }
#' sumstat_data contains the following components if exist
#' \itemize{
#'   \item sumstats: A list of summary statistics f or the matched LD_info, each sublist contains sumstats, n, var_y from \code{load_rss_data}.
#'   \item LD_info: A list of LD information, each sublist contains LD_variants, LD_matrix, ref_panel  \code{load_LD_matrix}.
#' }
#'
#' @export
colocboost_analysis_pipeline <- function(region_data,
                                         focal_trait = NULL,
                                         event_filters = NULL,
                                         # - analysis
                                         xqtl_coloc = TRUE,
                                         joint_gwas = FALSE,
                                         separate_gwas = FALSE,
                                         # - individual QC
                                         maf_cutoff = 0.0005,
                                         pip_cutoff_to_skip_ind = 0,
                                         # - sumstat QC
                                         keep_indel = TRUE,
                                         pip_cutoff_to_skip_sumstat = 0,
                                         qc_method = c("slalom", "dentist", "none"),
                                         impute = TRUE,
                                         impute_opts = list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5, lamb = 0.01),
                                         ...) {
  # - internal function by filtering events based on event_filters
  filter_events <- function(events, filters, condition) {
    # filters is a list of filter specifications
    # Each filter spec must have:
    #   type_pattern: pattern to identify event type
    #   And at least ONE of:
    #   valid_pattern: pattern that must exist in group
    #   exclude_pattern: pattern to exclude

    filtered_events <- events
    for (filter in filters) {
      if (is.null(filter$type_pattern) ||
        (is.null(filter$valid_pattern) && is.null(filter$exclude_pattern))) {
        stop("Each filter must specify type_pattern and at least one of valid_pattern or exclude_pattern")
      }
      # Get events of this type
      type_events <- filtered_events[grepl(filter$type_pattern, filtered_events)]

      if (length(type_events) == 0) next
      # Apply valid pattern if specified
      if (!is.null(filter$valid_pattern)) {
        valid_groups <- unique(gsub(
          filter$type_pattern, "\\1",
          type_events[grepl(filter$valid_pattern, type_events)]
        ))
        if (length(valid_groups) > 0) {
          type_events <- events[grepl(paste(valid_groups, collapse = "|"), type_events)]
        } else {
          type_events <- character(0)
        }
      }
      # Apply exclusions if specified
      if (!is.null(filter$exclude_pattern)) {
        type_events <- type_events[!grepl(filter$exclude_pattern, type_events)]
      }
      if (length(type_events) == length(events)) {
        message(paste("All events matching", filter$type_pattern, "in", condition, "included in following analysis."))
      } else if (length(type_events) == 0) {
        message(paste("No events matching", filter$type_pattern, "in", condition, "pass the filtering."))
        return(NULL)
      } else {
        exclude_events <- paste0(setdiff(events, type_events), collapse = ";")
        message(paste("Some events,", exclude_events, "in", condition, "are removed."))
      }
      # Update events list
      filtered_events <- unique(c(
        filtered_events[!grepl(filter$type_pattern, filtered_events)],
        type_events
      ))
    }

    return(filtered_events)
  }

  # - extract contexts and studies from region data, handling both pre- and post-QC
  extract_contexts_studies <- function(region_data, phenotypes_init = NULL) {
    individual_data <- region_data$individual_data
    sumstat_data <- region_data$sumstat_data
    phenotypes <- list("individual_contexts" = NULL, "sumstat_studies" = NULL)

    # Extract individual contexts
    if (!is.null(individual_data)) {
      if (is.null(phenotypes_init)) {
        phenotypes$individual_contexts <- names(individual_data$residual_Y)
      } else {
        null_Y <- which(sapply(individual_data$Y, is.null))
        if (length(null_Y) == 0) {
          message("All individual data pass QC steps.")
          phenotypes$individual_contexts <- names(individual_data$Y)
        } else if (length(null_Y) < length(individual_data$Y)) {
          message(paste(
            "Skipping follow-up analysis for individual traits",
            paste(names(individual_data$Y)[null_Y], collapse = ";"), "after QC."
          ))
          phenotypes$individual_contexts <- names(individual_data$Y)[-null_Y]
        } else {
          message("No individual data pass QC.")
        }
      }
    } else {
      message(if (is.null(phenotypes_init)) "No individual data in this region!" else "No individual data pass QC.")
    }

    # Extract sumstat studies
    if (!is.null(sumstat_data)) {
      if (is.null(phenotypes_init)) {
        phenotypes$sumstat_studies <- unlist(sapply(sumstat_data$sumstats, names))
      } else {
        phenotypes$sumstat_studies <- names(sumstat_data$sumstats)
        if (length(phenotypes_init$sumstat_studies) == length(phenotypes$sumstat_studies)) {
          message("All sumstat studies pass QC steps.")
        } else {
          message(paste(
            "Skipping follow-up analysis for sumstat studies",
            paste(setdiff(phenotypes_init$sumstat_studies, phenotypes$sumstat_studies), collapse = ";"), "after QC."
          ))
        }
      }
    } else {
      message(if (is.null(phenotypes_init)) "No sumstat data in this region!" else "No sumstat data pass QC.")
    }

    return(phenotypes)
  }

  ####### ========= resolve defaults ======== #######
  qc_method <- match.arg(qc_method)

  ####### ========= initial output results before QC ======== #######
  analysis_results <- list("xqtl_coloc" = NULL, "joint_gwas" = NULL, "separate_gwas" = NULL)
  analysis_results$computing_time <- list("QC" = NULL, "Analysis" = list("xqtl_coloc" = NULL, "joint_gwas" = NULL, "separate_gwas" = NULL))
  if (!xqtl_coloc & !joint_gwas & !separate_gwas) {
    message("No colocalization has been performed!")
    return(analysis_results)
  }
  phenotypes_init <- extract_contexts_studies(region_data)
  if (is.null(phenotypes_init$individual_contexts) & is.null(phenotypes_init$sumstat_studies)) {
    return(analysis_results)
  }
  if (!is.null(phenotypes_init$individual_contexts)) {
    analysis_results$xqtl_coloc <- list(NULL)
  }
  if (!is.null(phenotypes_init$sumstat_studies)) {
    analysis_results$joint_gwas <- list(NULL)
    if (length(phenotypes_init$sumstat_studies) > 1) {
      analysis_results$separate_gwas <- vector("list", length(phenotypes_init$sumstat_studies)) %>% setNames(phenotypes_init$sumstat_studies)
    } else {
      analysis_results$separate_gwas[[1]] <- list(NULL)
      names(analysis_results$separate_gwas) <- phenotypes_init$sumstat_studies
    }
  }

  ####### ========= Filtering events before QC =========== #########
  if (!is.null(event_filters) & !is.null(region_data$individual_data)) {
    Y <- region_data$individual_data$residual_Y
    Y <- lapply(seq_along(Y), function(i) {
      y <- Y[[i]]
      events <- colnames(y)
      condition <- names(Y)[i]
      filtered_events <- filter_events(events, event_filters, condition)
      if (is.null(filtered_events)) {
        return(NULL)
      }
      y[, filtered_events, drop = FALSE]
    }) %>% setNames(names(region_data$individual_data$residual_Y))
    region_data$individual_data$residual_Y <- Y
  }

  ####### ========= QC for the region_data ======== ########
  t01 <- Sys.time()
  region_data <- qc_regional_data(region_data,
    maf_cutoff = maf_cutoff,
    pip_cutoff_to_skip_ind = pip_cutoff_to_skip_ind,
    keep_indel = keep_indel,
    pip_cutoff_to_skip_sumstat = pip_cutoff_to_skip_sumstat,
    qc_method = qc_method,
    impute = impute,
    impute_opts = impute_opts
  )
  phenotypes_QC <- extract_contexts_studies(region_data, phenotypes_init = phenotypes_init)
  t02 <- Sys.time()
  analysis_results$computing_time$QC <- t02 - t01

  ####### ========= organize individual level data ======== ########
  individual_data <- region_data$individual_data
  if (!is.null(individual_data)) {
    X <- individual_data$X
    Y <- individual_data$Y
    null_Y <- which(sapply(Y, is.null))
    if (length(null_Y) != 0 & length(null_Y) != length(Y)) {
      X <- X[-null_Y]
      Y <- Y[-null_Y]
    } else if (length(null_Y) == length(Y)) {
      X <- NULL
      Y <- NULL
    }
    if (!is.null(Y)) {
      Y_split <- purrr::imap(Y, function(y, i) {
        purrr::map(seq_len(ncol(y)), function(j) setNames(y[, j, drop = FALSE], colnames(y)[j]))
      })
      dict_YX <- cbind(seq_along(unlist(Y_split, recursive = FALSE)), rep(seq_along(Y_split), purrr::map_int(Y_split, length)))
      Y <- unlist(Y_split, recursive = FALSE)
      names(Y) <- sapply(Y, colnames)
    }
  } else {
    X <- Y <- dict_YX <- NULL
  }

  ####### ========= organize summary statistics ======== ########
  sumstat_data <- region_data$sumstat_data
  if (!is.null(sumstat_data$sumstats)) {
    sumstats <- lapply(sumstat_data$sumstats, function(ss) {
      z <- ss$sumstats$z
      # Normalize variant IDs to canonical format (with chr prefix)
      variant <- normalize_variant_id(as.character(ss$sumstats$variant_id))
      n <- ss$n

      # Filter out NA values from z-scores and corresponding variants
      na_mask <- !is.na(z)
      if (sum(na_mask) == 0) {
        message("Warning: All z-scores are NA for this summary statistic dataset")
        return(data.frame("z" = numeric(0), "n" = numeric(0), "variant" = character(0)))
      }

      data.frame("z" = z[na_mask], "n" = n, "variant" = variant[na_mask])
    })
    names(sumstats) <- names(sumstat_data$sumstats)
    LD_mat <- lapply(sumstat_data$LD_mat, function(ld) {
      # Normalize LD dimnames to canonical format (with chr prefix)
      if (!is.null(colnames(ld))) {
        colnames(ld) <- normalize_variant_id(as.character(colnames(ld)))
      }
      if (!is.null(rownames(ld))) {
        rownames(ld) <- normalize_variant_id(as.character(rownames(ld)))
      }
      return(ld)
    })
    LD_match <- sumstat_data$LD_match
    dict_sumstatLD <- cbind(seq_along(sumstats), match(LD_match, names(sumstat_data$LD_mat)))
    # Validate and filter sumstat entries before analysis
    filtered <- filter_valid_sumstats(sumstats, LD_mat, LD_match)
    if (is.null(filtered)) {
      sumstats <- LD_mat <- dict_sumstatLD <- NULL
    } else {
      sumstats <- filtered$sumstats
      LD_mat <- filtered$LD_mat
      LD_match <- filtered$LD_match
      dict_sumstatLD <- filtered$dict_sumstatLD
    }
  } else {
    sumstats <- LD_mat <- dict_sumstatLD <- NULL
  }


  ####### ========= streamline three types of analyses ======== ########
  if (is.null(X) & is.null(sumstats)) {
    message("No data pass QC and will not perform analyses.")
    return(analysis_results)
  }
  # - run xQTL-only version of ColocBoost
  if (xqtl_coloc & !is.null(X)) {
    message(paste("====== Performing xQTL-only ColocBoost on", length(Y), "contexts. ====="))
    traits <- names(Y)
    focal_outcome_idx <- if (!is.null(focal_trait) && focal_trait %in% traits) which(traits == focal_trait) else NULL
    cb_res <- .run_colocboost("xQTL-only ColocBoost",
      X = X, Y = Y, dict_YX = dict_YX,
      outcome_names = traits, focal_outcome_idx = focal_outcome_idx,
      output_level = 2, ...
    )
    analysis_results$xqtl_coloc <- cb_res$result
    analysis_results$computing_time$Analysis$xqtl_coloc <- cb_res$time
  }
  # - run joint GWAS no focaled version of ColocBoost
  if (joint_gwas & !is.null(sumstats)) {
    message(paste("====== Performing non-focaled version GWAS-xQTL ColocBoost on", length(Y), "contexts and", length(sumstats), "GWAS. ====="))
    traits <- c(names(Y), names(sumstats))
    ld_args <- build_ld_args(LD_mat)
    cb_res <- do.call(.run_colocboost, c(list("Joint GWAS ColocBoost",
      X = X, Y = Y, sumstat = sumstats,
      dict_YX = dict_YX, dict_sumstatLD = dict_sumstatLD,
      outcome_names = traits, focal_outcome_idx = NULL,
      output_level = 2), ld_args, list(...)))
    analysis_results$joint_gwas <- cb_res$result
    analysis_results$computing_time$Analysis$joint_gwas <- cb_res$time
  }
  # - run focaled version of ColocBoost for each GWAS
  if (separate_gwas & !is.null(sumstats)) {
    t31 <- Sys.time()
    res_gwas_separate <- analysis_results$separate_gwas
    for (i_gwas in 1:nrow(dict_sumstatLD)) {
      current_study <- names(sumstats)[i_gwas]
      message(paste("====== Performing focaled version GWAS-xQTL ColocBoost on", length(Y), "contexts and ", current_study, "GWAS. ====="))
      dict <- dict_sumstatLD[i_gwas, ]
      traits <- c(names(Y), current_study)
      ld_args_sep <- build_ld_args(LD_mat, subset = dict[2])
      cb_res <- do.call(.run_colocboost, c(
        list(paste("Separate GWAS ColocBoost for", current_study),
          X = X, Y = Y, sumstat = sumstats[dict[1]],
          dict_YX = dict_YX,
          outcome_names = traits, focal_outcome_idx = length(traits),
          output_level = 2), ld_args_sep, list(...)))
      res_gwas_separate[[current_study]] <- cb_res$result
    }
    t32 <- Sys.time()
    analysis_results$separate_gwas <- res_gwas_separate
    analysis_results$computing_time$Analysis$separate_gwas <- list("total" = t32 - t31, "n_studies" = nrow(dict_sumstatLD), "average" = (t32 - t31) / nrow(dict_sumstatLD))
  }

  return(analysis_results)
}


#' Validate a single summary statistics entry
#'
#' Checks whether a sumstat data frame (with columns z, n, variant) has
#' sufficient data for colocboost analysis.
#'
#' @param ss_df A data.frame with columns "z", "n", and "variant" as produced
#'   by the sumstat processing block in \code{colocboost_analysis_pipeline}.
#' @param min_variants Minimum number of non-NA z-score variants required.
#'   Default is 2.
#' @return TRUE if the entry is valid; FALSE otherwise.
#' @noRd
is_valid_sumstat_entry <- function(ss_df, min_variants = 2) {
  if (is.null(ss_df) || !is.data.frame(ss_df)) return(FALSE)
  if (nrow(ss_df) < min_variants) return(FALSE)
  if (all(is.na(ss_df$z))) return(FALSE)
  if (all(ss_df$n <= 0 | is.na(ss_df$n))) return(FALSE)
  return(TRUE)
}


#' Filter summary statistics to retain only valid entries
#'
#' Applies \code{is_valid_sumstat_entry} to each element of the sumstats list
#' and removes invalid entries. Also updates the LD match mapping and
#' dict_sumstatLD accordingly.
#'
#' @param sumstats Named list of summary statistic data frames.
#' @param LD_mat List of LD matrices.
#' @param LD_match Character vector mapping each sumstat to an LD matrix name.
#' @param min_variants Minimum variant count passed to
#'   \code{is_valid_sumstat_entry}.
#' @return A list with filtered \code{sumstats}, \code{LD_mat}, \code{LD_match},
#'   and \code{dict_sumstatLD}, or NULL if no valid entries remain.
#' @noRd
filter_valid_sumstats <- function(sumstats, LD_mat, LD_match, min_variants = 2) {
  valid_idx <- vapply(sumstats, is_valid_sumstat_entry, logical(1),
                       min_variants = min_variants)
  if (!any(valid_idx)) {
    message("No valid summary statistic studies remain after validation.")
    return(NULL)
  }
  removed <- names(sumstats)[!valid_idx]
  if (length(removed) > 0) {
    message("Removed invalid sumstat studies: ", paste(removed, collapse = ", "))
  }
  sumstats <- sumstats[valid_idx]
  LD_match <- LD_match[valid_idx]
  dict_sumstatLD <- cbind(seq_along(sumstats), match(LD_match, names(LD_mat)))
  list(sumstats = sumstats, LD_mat = LD_mat, LD_match = LD_match,
       dict_sumstatLD = dict_sumstatLD)
}


#' Initial QC for the region data loading from \code{load_regional_data}
#'
#' This function do the initial QC including: check PIP; check maf for individual_data; check QC and impute for sumstat_data
#'
#' @section Loading individual level data from multiple corhorts
#' @param region_data A region data loaded from \code{load_regional_data}.
#' @param maf_cutoff A scalar to remove variants with maf < maf_cutoff, dafault is 0.005.
#' @param pip_cutoff_to_skip_ind A vector of cutoff values for skipping analysis based on PIP values for each context. Default is 0.
#' @param pip_cutoff_to_skip_sumstat A vector of cutoff values for skipping analysis based on PIP values for each sumstat Default is 0.
#' @param qc_method Quality control method to use. Options are "slalom" or "dentist" (default: "slalom").
#' @param impute Logical; if TRUE, performs imputation for outliers identified in the analysis (default: TRUE).
#' @param impute_opts A list of imputation options including rcond, R2_threshold, and minimum_ld (default: list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5)).
#'
#'
#' @return A list containing the individual_data and sumstat_data after QC:
#' individual_data contains the following components if exist
#' \itemize{
#'   \item Y: A list of residualized phenotype values for all tasks.
#'   \item X: A list of residualized genotype matrices all tasks.
#' }
#' sumstat_data contains the following components if exist
#' \itemize{
#'   \item sumstats: A list of summary statistics for the matched LD_info, each sublist contains sumstats, n, var_y from \code{load_rss_data}.
#'   \item LD_info: A list of LD information, each sublist contains LD_variants, LD_matrix, ref_panel  \code{load_LD_matrix}.
#' }
#'
#' @noRd
qc_regional_data <- function(region_data,
                             # - individual
                             maf_cutoff = 0.0005,
                             pip_cutoff_to_skip_ind = 0,
                             # - sumstat
                             keep_indel = TRUE,
                             pip_cutoff_to_skip_sumstat = 0,
                             qc_method = c("slalom", "dentist", "none"),
                             impute = TRUE,
                             impute_opts = list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5, lamb = 0.01)) {
  qc_method <- match.arg(qc_method)

  # Validate and recycle pip_cutoff_to_skip_ind: scalar -> recycled to n_contexts
  if (!is.null(region_data$individual_data)) {
    n_ind_contexts <- length(region_data$individual_data$residual_Y)
    if (length(pip_cutoff_to_skip_ind) == 1) {
      pip_cutoff_to_skip_ind <- rep(pip_cutoff_to_skip_ind, n_ind_contexts)
    } else if (length(pip_cutoff_to_skip_ind) != n_ind_contexts) {
      stop("pip_cutoff_to_skip_ind must be a scalar or a vector with length equal to the number of individual contexts (", n_ind_contexts, ").")
    }
  }

  # Validate pip_cutoff_to_skip_sumstat: scalar -> named vector for all studies
  if (!is.null(region_data$sumstat_data)) {
    all_study_names <- unlist(lapply(region_data$sumstat_data$sumstats, names))
    if (length(pip_cutoff_to_skip_sumstat) == 1 && is.null(names(pip_cutoff_to_skip_sumstat))) {
      pip_cutoff_to_skip_sumstat <- setNames(rep(pip_cutoff_to_skip_sumstat, length(all_study_names)), all_study_names)
    } else if (!is.null(names(pip_cutoff_to_skip_sumstat))) {
      # Named vector: fill missing studies with 0
      missing <- setdiff(all_study_names, names(pip_cutoff_to_skip_sumstat))
      if (length(missing) > 0) {
        pip_cutoff_to_skip_sumstat <- c(pip_cutoff_to_skip_sumstat, setNames(rep(0, length(missing)), missing))
      }
    }
  }

  #### related internal functions
  # Add context names to colname of Y if missing
  add_context_to_Y <- function(res_Y) {
    res <- lapply(seq_along(res_Y), function(iy) {
      y <- res_Y[[iy]]
      if (is.null(y)) {
        return(NULL)
      }
      if (is.null(colnames(y))) {
        colnames(y) <- names(res_Y)[iy]
      } else {
        colnames(y) <- paste0(names(res_Y)[iy], "_", colnames(y))
      }
      return(y)
    })
    names(res) <- names(res_Y)
    return(res)
  }

  # Initial PIP check for individual level data
  filter_resY_pip <- function(res_X, res_Y, pip_cutoff = 0, context = NULL) {
    # Initial PIP check
    if (pip_cutoff != 0) {
      if (pip_cutoff < 0) {
        # automatically determine the cutoff to use
        pip_cutoff <- 3 * 1 / ncol(res_X)
      }
      top_model_pip <- lapply(1:ncol(res_Y), function(i) susieR::susie(res_X, res_Y[, i], L = 1)$pip)
      check_model_pip <- sapply(top_model_pip, function(pip) any(pip > pip_cutoff))
      include_idx <- which(check_model_pip)
      if (length(include_idx) == 0) {
        message(paste(
          "Skipping follow-up analysis for individual-context", context,
          ". No signals above PIP threshold", pip_cutoff, "in initial model screening."
        ))
        return(NULL)
      } else if (length(include_idx) == ncol(res_Y)) {
        message(paste("Keep all individual-phenotypes in context", context, "."))
      } else {
        exclude_idx <- setdiff(1:ncol(res_Y), include_idx)
        exclude_pheno <- paste(colnames(res_Y)[exclude_idx], collapse = ";")
        message(paste(
          "Skipping follow-up analysis for individual-phenotypes", exclude_pheno, "in context", context,
          ". No signals above PIP threshold", pip_cutoff, "in initial model screening."
        ))
        res_Y <- res_Y[, include_idx, drop = FALSE] %>% .[, order(colnames(.)), drop = FALSE]
      }
    }
    return(res_Y)
  }

  # Initial check for all contexts with individual-level data
  data_initial_screen_individual <- function(X,
                                             Y,
                                             MAF,
                                             maf_cutoff = 0.0005,
                                             pip_cutoff_to_skip_ind = 0) {
    # - add context to colname of Y
    Y <- add_context_to_Y(Y)
    results <- purrr::imap(X, function(resX, context) {
      resY <- Y[[context]]
      maf <- MAF[[context]]
      i_context <- match(context, names(X))
      if (is.null(resY)) return(NULL)
      resX <- filter_X(resX, missing_rate_thresh = NULL, maf_thresh = maf_cutoff, maf = maf)
      resY <- filter_resY_pip(resX, resY, pip_cutoff = pip_cutoff_to_skip_ind[i_context], context = context)
      if (is.null(resY)) return(NULL)
      list(X = resX, Y = resY)
    }) %>% purrr::compact()

    if (length(results) == 0) {
      message("Skipping follow-up analysis for all contexts.")
      return(NULL)
    }
    keep_contexts <- names(results)
    message(paste("Region includes the following contexts after inital screening:", paste(keep_contexts, collapse = ";"), "."))
    list(
      X = purrr::map(results, "X"),
      Y = purrr::map(results, "Y")
    )
  }

  # - individual level data QC
  individual_data <- region_data$individual_data
  if (!is.null(individual_data)) {
    X <- individual_data$residual_X
    Y <- individual_data$residual_Y
    MAF <- individual_data$maf
    # 1. remove maf < maf_cutoff
    # 2. initial check PIP
    individual_data <- data_initial_screen_individual(
      X = X, Y = Y, MAF = MAF,
      maf_cutoff = maf_cutoff,
      pip_cutoff_to_skip_ind = pip_cutoff_to_skip_ind
    )
  }


  # sumstat_data QC, imputation using raiss
  # return A list of sumstat_data after initial checking and QC:
  # \itemize{
  #   \item sumstats: A list of summary statistics and ready to do analysis.
  #   \item LD_mat: A list of LD matrix and ready to do analysis.
  #   \item LD_match: A vector of strings to indicating sumstats and LD matching (save space since multiple sumstats may link to the same LD matrix).
  # }
  summary_stats_qc_multitask <- function(sumstat_data,
                                         keep_indel = TRUE,
                                         pip_cutoff_to_skip_sumstat = 0,
                                         qc_method = c("slalom", "dentist", "none"),
                                         impute = TRUE,
                                         impute_opts = list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5, lamb = 0.01)) {
    n_LD <- length(sumstat_data$LD_info)
    # Collect results into lists and flatten at the end
    collected_sumstats <- list()
    collected_LD <- list()
    collected_LD_match <- character()
    # Track LD matrices by variant signature for O(1) deduplication
    ld_variant_index <- list()

    for (i in 1:n_LD) {
      LD_data <- sumstat_data$LD_info[[i]]
      sumstats <- sumstat_data$sumstats[[i]]
      has_genotype <- isTRUE(LD_data$is_genotype)

      # When source is genotype X, derive R only where needed (QC, imputation).
      # Keep X as primary for colocboost X_ref.
      LD_data_for_qc <- LD_data
      if (has_genotype) {
        LD_data_for_qc$LD_matrix <- compute_LD(LD_data$LD_matrix, method = "sample")
        LD_data_for_qc$is_genotype <- FALSE
      }

      # Pre-compute LD partition once per block (shared across all GWAS studies)
      if (impute) {
        LD_matrix_partitioned <- partition_LD_matrix(LD_data_for_qc)
      }

      for (ii in seq_along(sumstats)) {
        sumstat <- sumstats[[ii]]
        if (nrow(sumstat$sumstats) == 0) next
        n <- sumstat$n
        conditions_sumstat <- names(sumstats)[ii]
        pip_cutoff_to_skip_ld <- if (conditions_sumstat %in% names(pip_cutoff_to_skip_sumstat)) {
          as.numeric(pip_cutoff_to_skip_sumstat[conditions_sumstat])
        } else {
          0
        }

        # Preprocess: allele QC + variant subsetting (needs R for [variants, variants] indexing)
        preprocess_results <- rss_basic_qc(sumstat$sumstats, LD_data_for_qc, keep_indel = keep_indel)
        sumstat$sumstats <- preprocess_results$sumstats
        R_mat <- preprocess_results$LD_mat

        # Initial PIP checking (uses X when available, R otherwise)
        if (pip_cutoff_to_skip_ld != 0) {
          pip_vars <- sumstat$sumstats$variant_id
          if (has_genotype) {
            pip <- susie_rss_wrapper(z = sumstat$sumstats$z,
              X = LD_data$LD_matrix[, pip_vars, drop = FALSE], L = 1, n = n)$pip
          } else {
            pip <- susie_rss_wrapper(z = sumstat$sumstats$z, R = R_mat, L = 1, n = n)$pip
          }
          if (pip_cutoff_to_skip_ld < 0) {
            pip_cutoff_to_skip_ld <- 3 * 1 / length(pip_vars)
          }
          if (!any(pip > pip_cutoff_to_skip_ld)) {
            message(paste(
              "Skipping follow-up analysis for sumstat study", conditions_sumstat,
              ". No signals above PIP threshold", pip_cutoff_to_skip_ld, "in initial model screening."
            ))
            next
          } else {
            message(paste("Keep summary study", conditions_sumstat, "."))
          }
        }

        # Quality control — remove outlier variants (needs R)
        if (!is.null(qc_method) && qc_method != "none") {
          qc_results <- summary_stats_qc(sumstat$sumstats, LD_data_for_qc, n = n, method = qc_method)
          sumstat$sumstats <- qc_results$sumstats
          R_mat <- qc_results$LD_mat
        }
        # Imputation (needs R via partitioned LD)
        if (impute) {
          impute_results <- raiss(LD_data_for_qc$ref_panel, sumstat$sumstats, LD_matrix_partitioned,
            rcond = impute_opts$rcond,
            R2_threshold = impute_opts$R2_threshold, minimum_ld = impute_opts$minimum_ld, lamb = impute_opts$lamb
          )
          sumstat$sumstats <- impute_results$result_filter
          R_mat <- impute_results$LD_mat
        }

        # Store: X subset if genotype source, R otherwise
        final_vars <- sumstat$sumstats$variant_id
        if (has_genotype) {
          missing <- setdiff(final_vars, colnames(LD_data$LD_matrix))
          if (length(missing) > 0) {
            stop("BUG: ", length(missing), " QC'd variants not found in genotype matrix X. ",
                 "First few: ", paste(head(missing, 3), collapse = ", "))
          }
          mat_to_store <- LD_data$LD_matrix[, final_vars, drop = FALSE]
        } else {
          mat_to_store <- R_mat
        }

        # Collect sumstat
        collected_sumstats[[conditions_sumstat]] <- sumstat

        # Deduplicate LD using variant signature hash
        variant_key <- paste(colnames(mat_to_store), collapse = ",")
        if (variant_key %in% names(ld_variant_index)) {
          collected_LD_match <- c(collected_LD_match, ld_variant_index[[variant_key]])
        } else {
          collected_LD[[conditions_sumstat]] <- mat_to_store
          ld_variant_index[[variant_key]] <- conditions_sumstat
          collected_LD_match <- c(collected_LD_match, conditions_sumstat)
        }
      }
    }
    return(list(sumstats = collected_sumstats, LD_mat = collected_LD, LD_match = collected_LD_match))
  }


  # - summary statistics QC
  sumstat_data <- region_data$sumstat_data
  if (!is.null(sumstat_data)) {
    # - initial check PIP, qc or impute
    sumstat_data <- summary_stats_qc_multitask(sumstat_data,
      keep_indel = keep_indel,
      pip_cutoff_to_skip_sumstat = pip_cutoff_to_skip_sumstat,
      qc_method = qc_method,
      impute = impute, impute_opts = impute_opts
    )
  }
  return(list(
    individual_data = individual_data,
    sumstat_data = sumstat_data
  ))
}
