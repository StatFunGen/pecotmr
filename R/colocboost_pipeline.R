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

#' Run colocboost with tryCatch and timing.
#' @importFrom colocboost colocboost
#' @noRd
.run_colocboost <- function(label, ...) {
  t1 <- Sys.time()
  res <- tryCatch(
    colocboost_analysis(...),
    error = function(e) {
      message(label, " failed: ", conditionMessage(e))
      NULL
    }
  )
  list(result = res, time = Sys.time() - t1)
}

.cb_call_colocboost <- function(args, dots) {
  if (!requireNamespace("colocboost", quietly = TRUE)) {
    stop("The colocboost package is required for colocboost_analysis().")
  }
  do.call(colocboost, c(args, dots))
}

#' Convert loaded regional data to ColocBoost inputs
#'
#' @param region_data A list returned by \code{load_multitask_regional_data()}.
#' @return A structured list containing \code{colocboost_input},
#'   \code{qc_input}, and \code{source_info}.
#' @export
region_data_to_colocboost_input <- function(region_data) {
  ind_records_from_input <- function(input) {
    X <- .cb_as_named_list(input$X, "individual")
    Y <- .cb_as_named_list(input$Y, "individual")
    contexts <- intersect(names(X), names(Y))
    records <- list()
    for (context in contexts) {
      if (is.null(X[[context]]) || .cb_y_ncol(Y[[context]]) == 0) next
      records[[context]] <- list(
        X = X[[context]],
        Y = Y[[context]],
        maf = .cb_list_value(input$maf, context),
        X_variance = .cb_list_value(input$X_variance, context)
      )
    }
    records
  }

  ind_input <- region_data_to_ind_input(region_data)
  rss_input <- region_data_to_rss_input(region_data)

  ind_records <- ind_records_from_input(ind_input)
  ind_args <- .cb_format_individual(ind_records)

  sumstat_records <- lapply(names(rss_input$rss_input), function(study) {
    ld_data <- .normalize_ld_data_for_qc(rss_input$LD_data[[study]])
    list(rss_input = rss_input$rss_input[[study]],
         LD_matrix = ld_data$LD_matrix)
  })
  names(sumstat_records) <- names(rss_input$rss_input)
  sumstat_args <- .cb_format_sumstat(sumstat_records)

  outcome_names <- c(ind_args$outcome_names, names(sumstat_args$sumstat))
  ind_args$outcome_names <- NULL
  colocboost_input <- .cb_merge_args(ind_args, sumstat_args)
  if (length(outcome_names) > 0) colocboost_input$outcome_names <- outcome_names

  list(
    colocboost_input = Filter(Negate(is.null), colocboost_input),
    qc_input = list(
      individual = ind_input[c("X", "Y", "maf", "X_variance")],
      sumstat = rss_input[c("rss_input", "LD_data")]
    ),
    source_info = list(individual = ind_input$source_info,
                       sumstat = rss_input$source_info)
  )
}

#' ColocBoost analysis with optional pipeline QC
#'
#' This wrapper keeps the direct \code{colocboost()} argument surface. All
#' ColocBoost inputs and model parameters are supplied through \code{...}. When
#' no QC options are requested, the call is passed directly to
#' \code{colocboost()}. When QC options are requested, the wrapper
#' inspects named \code{X}/\code{Y} and/or \code{sumstat}/\code{LD}/\code{X_ref}
#' arguments in \code{...}, runs the relevant reusable QC step, and then calls
#' ColocBoost on the cleaned inputs. If the required named inputs are not
#' available, QC is skipped with a warning and the original ColocBoost call is
#' used.
#'
#' @details
#' Use \code{colocboost_analysis()} the same way you would use
#' \code{colocboost()}: pass the native ColocBoost arguments by
#' name, for example \code{X}, \code{Y}, \code{sumstat}, \code{LD},
#' \code{X_ref}, \code{dict_YX}, \code{dict_sumstatLD},
#' \code{outcome_names}, \code{focal_outcome_idx}, \code{effect_est},
#' \code{effect_se}, \code{effect_n}, \code{M}, and other ColocBoost model or
#' post-processing options. These arguments are forwarded unchanged unless one
#' or more QC controls are requested.
#'
#' Individual-level QC is only attempted when at least one individual QC control
#' is non-\code{NULL} and named \code{X} and \code{Y} inputs are available in
#' \code{...}. Summary-statistic QC is only attempted when \code{qc_method},
#' \code{pip_cutoff_to_skip_sumstat}, \code{impute = TRUE}, or
#' \code{LD_reference_info} is supplied and named \code{sumstat} plus either
#' \code{LD}, \code{X_ref}, or \code{LD_reference_info} are available.
#' \code{qc_method = "none"} means run basic allele/variant harmonization
#' only; it does not run SLALOM/DENTIST
#' LD-mismatch QC. RAISS imputation is controlled separately by
#' \code{impute = TRUE}.
#'
#' If no QC controls are supplied, this function is a thin direct call to
#' \code{colocboost(...)}.
#' When QC removes outcomes, \code{outcome_names} and \code{focal_outcome_idx}
#' are updated to match the post-QC outcome order. If the requested focal outcome
#' is removed by QC, \code{focal_outcome_idx} is set to \code{NULL} with a
#' warning.
#'
#' @param ... Arguments passed to \code{colocboost()}, including
#'   data inputs such as \code{X}, \code{Y}, \code{sumstat}, \code{LD},
#'   \code{X_ref}, \code{dict_YX}, \code{dict_sumstatLD},
#'   \code{outcome_names}, and all ColocBoost model/post-processing options.
#'   QC can only inspect inputs that are supplied by name.
#' @param missing_rate_thresh,maf_cutoff,xvar_cutoff,ld_reference_meta_file,pip_cutoff_to_skip_ind
#'   Individual-level QC controls. If all are \code{NULL}, individual-level QC
#'   is not run.
#' @param keep_indel,pip_cutoff_to_skip_sumstat,qc_method,impute,impute_opts
#'   Summary-statistic QC controls. \code{qc_method = "none"} runs
#'   basic allele harmonization without
#'   LD-mismatch outlier detection. Imputation is only run when
#'   \code{impute = TRUE}.
#' @param LD_reference_info Optional LD reference information for
#'   summary-statistic QC. This is only needed when the native \code{LD} matrix
#'   row/column names or \code{X_ref} column names are missing or are not
#'   parseable genomic variant IDs. It can be a .bim/.pvar/.pvar.zst file path,
#'   a data.frame with variant metadata, or a \code{load_LD_matrix()} result.
#'   This is a QC-only argument and is not passed to
#'   \code{colocboost()}.
#' @param variant_convention Allele order used by native ColocBoost-style
#'   \code{sumstat$variant} and LD/X_ref names when deriving QC inputs:
#'   \code{"A2_A1"} for pecotmr canonical \code{chr:pos:A2:A1}, or
#'   \code{"A1_A2"} for \code{chr:pos:A1:A2}.
#' @return The object returned by \code{colocboost()}.
#' @examples
#' \dontrun{
#' # Direct ColocBoost call without QC.
#' fit <- colocboost_analysis(X = X, Y = Y, M = 500)
#'
#' # Summary-statistic input with basic allele/variant harmonization only.
#' fit <- colocboost_analysis(sumstat = sumstat, LD = LD,
#'                            qc_method = "none", M = 500)
#'
#' # Summary-statistic input with LD-mismatch QC and RAISS imputation.
#' fit <- colocboost_analysis(sumstat = sumstat, LD = LD,
#'                            qc_method = "slalom", impute = TRUE)
#'
#' # Use richer LD metadata from load_LD_matrix() for QC, while still passing
#' # ColocBoost's native LD input.
#' ld_data <- load_LD_matrix(ld_meta_file, region)
#' fit <- colocboost_analysis(sumstat = sumstat, LD = ld_data$LD_matrix,
#'                            LD_reference_info = ld_data, qc_method = "none")
#'
#' # Individual-level input with explicit genotype QC thresholds.
#' fit <- colocboost_analysis(X = X, Y = Y,
#'                            missing_rate_thresh = 0.1,
#'                            maf_cutoff = 0.0005)
#' }
#' @export
colocboost_analysis <- function(...,
                                # individual QC
                                missing_rate_thresh = NULL,
                                maf_cutoff = NULL,
                                xvar_cutoff = NULL,
                                ld_reference_meta_file = NULL,
                                pip_cutoff_to_skip_ind = NULL,
                                # sumstat QC
                                keep_indel = TRUE,
                                pip_cutoff_to_skip_sumstat = NULL,
                                qc_method = NULL,
                                impute = FALSE,
                                impute_opts = list(rcond = 0.01, R2_threshold = 0.6,
                                                   minimum_ld = 5, lamb = 0.01),
                                LD_reference_info = NULL,
                                variant_convention = c("A2_A1", "A1_A2")) {
  variant_convention <- match.arg(variant_convention)
  direct_args <- list(...)
  pre_qc_data_outcomes <- .cb_colocboost_outcome_names(direct_args, prefer_supplied = FALSE)
  pre_qc_display_outcomes <- .cb_colocboost_outcome_names(direct_args, prefer_supplied = TRUE)
  if (!is.null(qc_method)) qc_method <- .resolve_summary_qc_method(qc_method)

  individual_qc_requested <- !is.null(missing_rate_thresh) ||
    !is.null(maf_cutoff) || !is.null(xvar_cutoff) ||
    !is.null(ld_reference_meta_file) || !is.null(pip_cutoff_to_skip_ind)
  sumstat_qc_requested <- !is.null(qc_method) || isTRUE(impute) ||
    !is.null(pip_cutoff_to_skip_sumstat) || !is.null(LD_reference_info)
  qc_requested <- individual_qc_requested || sumstat_qc_requested
  if (!qc_requested) {
    return(.cb_call_colocboost(direct_args, list()))
  }

  X <- direct_args$X
  Y <- direct_args$Y
  sumstat <- direct_args$sumstat
  LD <- direct_args$LD
  X_ref <- direct_args$X_ref
  dict_YX <- direct_args$dict_YX
  dict_sumstatLD <- direct_args$dict_sumstatLD

  qc_skip_messages <- character()
  individual_qc_input <- NULL
  if (individual_qc_requested) {
    if (!is.null(X) && !is.null(Y)) {
      individual_qc_input <- .cb_individual_qc_input_from_colocboost(X, Y, dict_YX)
    } else {
      qc_skip_messages <- c(qc_skip_messages,
                            "Individual-level QC requested but named X and Y were not both supplied.")
    }
  }

  sumstat_qc_input <- NULL
  if (sumstat_qc_requested) {
    if (!is.null(sumstat) && (!is.null(LD) || !is.null(X_ref) || !is.null(LD_reference_info))) {
      sumstat_qc_input <- tryCatch(
        .cb_sumstat_qc_input_from_colocboost(
          sumstat, LD, X_ref, dict_sumstatLD,
          LD_reference_info = LD_reference_info,
          variant_convention = variant_convention
        ),
        error = function(e) {
          qc_skip_messages <<- c(
            qc_skip_messages,
            paste("Summary-statistic QC input could not be prepared:", conditionMessage(e))
          )
          NULL
        }
      )
      if (!is.null(sumstat_qc_input)) {
        if (length(sumstat_qc_input$skip_reasons) > 0) {
          qc_skip_messages <- c(qc_skip_messages, sumstat_qc_input$skip_reasons)
        }
        if (length(sumstat_qc_input$rss_input) == 0) {
          sumstat_qc_input <- NULL
        }
      }
    } else {
      qc_skip_messages <- c(qc_skip_messages,
                            "Summary-statistic QC requested but named sumstat plus LD, X_ref, or LD_reference_info were not supplied.")
    }
  }
  if (is.null(individual_qc_input) && is.null(sumstat_qc_input)) {
    warning("QC requested but required QC inputs are unavailable. Calling colocboost() directly. ",
            paste(qc_skip_messages, collapse = " "))
    return(.cb_call_colocboost(direct_args, list()))
  }
  if (length(qc_skip_messages) > 0) {
    warning(paste(qc_skip_messages, collapse = " "), " Skipping unavailable QC branch.")
  }

  qc_args <- tryCatch({
    args <- list()
    if (!is.null(individual_qc_input)) {
      message("QC track: processing individual-level inputs before ColocBoost.")
      ind <- qc_individual_data(
        X = individual_qc_input$X,
        Y = individual_qc_input$Y,
        maf = individual_qc_input$maf,
        X_variance = individual_qc_input$X_variance,
        missing_rate_thresh = missing_rate_thresh,
        maf_cutoff = maf_cutoff,
        xvar_cutoff = .cb_default(xvar_cutoff, 0),
        ld_reference_meta_file = ld_reference_meta_file,
        keep_indel = keep_indel,
        pip_cutoff_to_skip = .cb_default(pip_cutoff_to_skip_ind, 0)
      )
      args <- .cb_merge_args(args, .cb_format_individual(ind))
    }
    if (!is.null(sumstat_qc_input) && length(sumstat_qc_input$rss_input) > 0) {
      message("QC track: processing summary-statistic inputs before ColocBoost.")
      sumstat_qc <- summary_stats_qc(
        rss_input = sumstat_qc_input$rss_input,
        LD_data = sumstat_qc_input$LD_data,
        keep_indel = keep_indel,
        pip_cutoff_to_skip = .cb_default(pip_cutoff_to_skip_sumstat, 0),
        qc_method = if (is.null(qc_method)) "none" else qc_method,
        impute = impute,
        impute_opts = impute_opts
      )
      args <- .cb_merge_args(args, .cb_format_sumstat(sumstat_qc))
    }
    args
  }, error = function(e) {
    warning("QC requested but skipped: ", conditionMessage(e),
            ". Calling colocboost() directly.")
    NULL
  })

  if (is.null(qc_args) || length(qc_args) == 0) {
    return(.cb_call_colocboost(direct_args, list()))
  }
  merged_args <- .cb_merge_args(direct_args, qc_args)
  if (!is.null(qc_args$LD)) merged_args$X_ref <- NULL
  if (!is.null(qc_args$X_ref)) merged_args$LD <- NULL
  post_qc_data_outcomes <- .cb_colocboost_outcome_names(merged_args, prefer_supplied = FALSE)
  if (length(post_qc_data_outcomes) > 0) {
    merged_args$outcome_names <- .cb_resolve_qc_outcome_names(
      pre_qc_data_outcomes,
      pre_qc_display_outcomes,
      post_qc_data_outcomes
    )
    merged_args$focal_outcome_idx <- .cb_remap_focal_outcome_idx(
      focal_outcome_idx = direct_args$focal_outcome_idx,
      pre_qc_data_outcomes = pre_qc_data_outcomes,
      pre_qc_display_outcomes = pre_qc_display_outcomes,
      post_qc_data_outcomes = post_qc_data_outcomes,
      post_qc_display_outcomes = merged_args$outcome_names
    )
  }
  .cb_call_colocboost(Filter(Negate(is.null), merged_args), list())
}

#' Multi-trait colocalization analysis protocol pipeline
#'
#' This function performs protocol-level multi-trait colocalization using
#' ColocBoost. It accepts loaded regional data, performs QC once, then runs the
#' requested xQTL-only, joint GWAS, and separate GWAS analyses.
#'
#' @param region_data A region data loaded from \code{load_regional_data}.
#' @param focal_trait Name of trait if perform focaled ColocBoost
#' @param event_filters A list of pattern for filtering events based on context names. Example: for sQTL, list(type_pattern = ".*clu_(\\d+_[+-?]).*",valid_pattern = "clu_(\\d+_[+-?]):PR:",exclude_pattern = "clu_(\\d+_[+-?]):IN:")
#' @param maf_cutoff A scalar to remove variants with maf < maf_cutoff, dafault is 0.005.
#' @param pip_cutoff_to_skip_ind A vector of cutoff values for skipping analysis based on PIP values for each context. Default is 0.
#' @param pip_cutoff_to_skip_sumstat A vector of cutoff values for skipping analysis based on PIP values for each sumstat Default is 0.
#' @param qc_method Quality control method to use. Options are "none",
#'   "slalom", or "dentist". \code{NULL} is treated as \code{"none"} for
#'   basic-only summary-stat preprocessing.
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
#' @importFrom susieR susie_rss
#' @importFrom purrr imap map_int
#' @export
colocboost_pipeline <- function(
  region_data,
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
  qc_method = NULL,
  impute = TRUE,
  impute_opts = list(
    rcond = 0.01, R2_threshold = 0.6,
    minimum_ld = 5, lamb = 0.01
  ),
  ...
) {
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
      if (is.null(phenotypes_init)) {
        message("No sumstat data in this region!")
      } else if (!is.null(phenotypes_init$sumstat_studies)) {
        message(paste(
          "Skipping follow-up analysis for sumstat studies",
          paste(phenotypes_init$sumstat_studies, collapse = ";"), "after QC."
        ))
      } else {
        message("No sumstat data pass QC.")
      }
    }

    return(phenotypes)
  }

  ####### ========= resolve defaults ======== #######
  qc_method <- .resolve_summary_qc_method(qc_method)

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
  if (!is.null(phenotypes_init$sumstat_studies) &&
      is.null(phenotypes_QC$sumstat_studies)) {
    message("No valid summary statistic studies remain after validation.")
  }
  t02 <- Sys.time()
  analysis_results$computing_time$QC <- t02 - t01

  ####### ========= convert QC'd regional data to ColocBoost input ======== ########
  colocboost_input <- region_data_to_colocboost_input(region_data)$colocboost_input
  X <- colocboost_input$X
  Y <- colocboost_input$Y
  dict_YX <- colocboost_input$dict_YX
  sumstats <- colocboost_input$sumstat
  dict_sumstatLD <- colocboost_input$dict_sumstatLD
  LD_mat <- colocboost_input$LD
  if (is.null(LD_mat)) LD_mat <- colocboost_input$X_ref


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
    analysis_results["xqtl_coloc"] <- list(cb_res$result)
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
    analysis_results["joint_gwas"] <- list(cb_res$result)
    analysis_results$computing_time$Analysis$joint_gwas <- cb_res$time
  }
  # - run focaled version of ColocBoost for each GWAS
  if (separate_gwas & !is.null(sumstats)) {
    t31 <- Sys.time()
    res_gwas_separate <- analysis_results$separate_gwas
    for (i_gwas in 1:nrow(dict_sumstatLD)) {
      dict <- dict_sumstatLD[i_gwas, ]
      current_study <- names(sumstats)[dict[1]]
      message(paste("====== Performing focaled version GWAS-xQTL ColocBoost on", length(Y), "contexts and ", current_study, "GWAS. ====="))
      traits <- c(names(Y), current_study)
      ld_args_sep <- build_ld_args(LD_mat, subset = dict[2])
      cb_res <- do.call(.run_colocboost, c(
        list(paste("Separate GWAS ColocBoost for", current_study),
          X = X, Y = Y, sumstat = sumstats[dict[1]],
          dict_YX = dict_YX,
          outcome_names = traits, focal_outcome_idx = length(traits),
          output_level = 2), ld_args_sep, list(...)))
      res_gwas_separate[current_study] <- list(cb_res$result)
    }
    t32 <- Sys.time()
    analysis_results$separate_gwas <- res_gwas_separate
    analysis_results$computing_time$Analysis$separate_gwas <- list("total" = t32 - t31, "n_studies" = nrow(dict_sumstatLD), "average" = (t32 - t31) / nrow(dict_sumstatLD))
  }

  return(analysis_results)
}

#' Initial QC for the region data loaded from \code{load_regional_data}
#'
#' This compatibility wrapper converts loaded regional data to reusable individual
#' and RSS inputs, runs the shared QC helpers once, and returns the historical
#' post-QC structure consumed by \code{colocboost_pipeline()}.
#'
#' @param region_data A region data loaded from \code{load_regional_data}.
#' @param maf_cutoff A scalar to remove variants with maf < maf_cutoff.
#' @param pip_cutoff_to_skip_ind A vector of cutoff values for skipping individual contexts.
#' @param pip_cutoff_to_skip_sumstat A vector of cutoff values for skipping summary-stat studies.
#' @param qc_method Quality control method to use. Options are "none",
#'   "slalom", or "dentist". \code{NULL} is treated as \code{"none"} for
#'   basic-only summary-stat preprocessing.
#' @param impute Logical; if TRUE, performs imputation when required metadata are available.
#' @param impute_opts A list of imputation options.
#' @return A list containing post-QC \code{individual_data} and \code{sumstat_data}.
#' @noRd
qc_regional_data <- function(region_data,
                             # - individual
                             maf_cutoff = 0.0005,
                             pip_cutoff_to_skip_ind = 0,
                             # - sumstat
                             keep_indel = TRUE,
                             pip_cutoff_to_skip_sumstat = 0,
                             qc_method = NULL,
                             impute = TRUE,
                             impute_opts = list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5, lamb = 0.01)) {
  qc_method <- .resolve_summary_qc_method(qc_method)
  qced_individual_to_region_data <- function(ind_qc) {
    if (is.null(ind_qc) || length(ind_qc) == 0) return(NULL)
    list(
      X = lapply(ind_qc, `[[`, "X"),
      Y = lapply(ind_qc, `[[`, "Y")
    )
  }
  qced_sumstat_to_region_data <- function(sumstat_qc) {
    if (is.null(sumstat_qc) || length(sumstat_qc) == 0) return(NULL)
    if (!is.null(sumstat_qc$rss_input) && !is.null(sumstat_qc$LD_matrix)) {
      sumstat_qc <- list(study1 = sumstat_qc)
    }
    sumstats <- lapply(sumstat_qc, `[[`, "rss_input")
    LD_mat <- list()
    LD_match <- character()
    ld_variant_index <- list()
    for (study in names(sumstat_qc)) {
      ld <- sumstat_qc[[study]]$LD_matrix
      variant_key <- paste(colnames(ld), collapse = ",")
      if (variant_key %in% names(ld_variant_index)) {
        LD_match <- c(LD_match, ld_variant_index[[variant_key]])
      } else {
        LD_mat[[study]] <- ld
        ld_variant_index[[variant_key]] <- study
        LD_match <- c(LD_match, study)
      }
    }
    list(sumstats = sumstats, LD_mat = LD_mat, LD_match = LD_match)
  }

  individual_data <- NULL
  ind_input <- region_data_to_ind_input(region_data)
  if (isTRUE(ind_input$source_info$has_individual)) {
    ind_qc <- qc_individual_data(
      X = ind_input$X,
      Y = ind_input$Y,
      maf = ind_input$maf,
      X_variance = ind_input$X_variance,
      maf_cutoff = maf_cutoff,
      pip_cutoff_to_skip = pip_cutoff_to_skip_ind
    )
    individual_data <- qced_individual_to_region_data(ind_qc)
  }

  sumstat_data <- NULL
  rss_input <- region_data_to_rss_input(region_data)
  if (isTRUE(rss_input$source_info$has_sumstat)) {
    sumstat_qc <- summary_stats_qc(
      rss_input = rss_input$rss_input,
      LD_data = rss_input$LD_data,
      keep_indel = keep_indel,
      pip_cutoff_to_skip = pip_cutoff_to_skip_sumstat,
      qc_method = qc_method,
      impute = impute,
      impute_opts = impute_opts
    )
    sumstat_data <- qced_sumstat_to_region_data(sumstat_qc)
  }

  list(individual_data = individual_data, sumstat_data = sumstat_data)
}

#' Run reusable individual-level QC
#'
#' @param X Genotype matrix or named list of genotype matrices.
#' @param Y Phenotype vector/matrix or named list of phenotype matrices.
#' @param maf Optional MAF vector or named list.
#' @param X_variance Optional variant variance vector or named list.
#' @param missing_rate_thresh Maximum missing genotype rate.
#' @param maf_cutoff Minimum MAF cutoff.
#' @param xvar_cutoff Minimum genotype variance cutoff.
#' @param ld_reference_meta_file Optional LD reference metadata file.
#' @param keep_indel Whether indel variants are kept during LD-reference
#'   filtering.
#' @param pip_cutoff_to_skip Optional single-effect PIP cutoff.
#' @return A named list of cleaned context-level \code{X}/\code{Y} records, or
#'   one cleaned record for matrix inputs.
#' @export
qc_individual_data <- function(X, Y, maf = NULL, X_variance = NULL,
                               missing_rate_thresh = NULL,
                               maf_cutoff = 0.0005,
                               xvar_cutoff = 0,
                               ld_reference_meta_file = NULL,
                               keep_indel = TRUE,
                               pip_cutoff_to_skip = 0,
                               context = NULL) {
  qc_one <- function(X, Y, maf = NULL, X_variance = NULL, context = NULL,
                     pip_cutoff_to_skip = 0) {
    if (is.null(X) || is.null(Y)) return(NULL)
    if (is.null(colnames(X))) stop("X must have variant colnames for individual-level QC.")
    if (is.vector(Y)) Y <- matrix(Y, ncol = 1)
    if (is.null(colnames(Y))) colnames(Y) <- .cb_default(context, paste0("outcome", seq_len(ncol(Y))))
    if (!is.null(context)) colnames(Y) <- paste0(context, "_", colnames(Y))

    message("QC track: starting individual-level QC for ", .cb_default(context, "individual data"), ".")
    original_variants <- colnames(X)
    if (!is.null(names(maf)) || length(maf) == length(original_variants)) {
      if (is.null(names(maf))) names(maf) <- original_variants
    }
    if (!is.null(names(X_variance)) || length(X_variance) == length(original_variants)) {
      if (is.null(names(X_variance))) names(X_variance) <- original_variants
    }
    if (!is.null(ld_reference_meta_file)) {
      reference_filter <- filter_variants_by_ld_reference(original_variants, ld_reference_meta_file,
                                                          keep_indel = keep_indel)
      X <- X[, reference_filter$data, drop = FALSE]
      if (!is.null(names(maf))) maf <- maf[colnames(X)]
      if (!is.null(names(X_variance))) X_variance <- X_variance[colnames(X)]
    }
    X <- filter_X(X, missing_rate_thresh = missing_rate_thresh,
                  maf_thresh = maf_cutoff, var_thresh = xvar_cutoff,
                  maf = maf, X_variance = X_variance)
    if (!is.null(names(maf))) maf <- maf[colnames(X)]
    if (!is.null(names(X_variance))) X_variance <- X_variance[colnames(X)]

    if (!is.null(pip_cutoff_to_skip) && pip_cutoff_to_skip != 0) {
      cutoff <- pip_cutoff_to_skip
      if (cutoff < 0) cutoff <- 3 / ncol(X)
      keep_y <- logical(ncol(Y))
      for (j in seq_len(ncol(Y))) {
        observed <- !is.na(Y[, j])
        if (sum(observed) < 2) next
        pip <- susieR::susie(X[observed, , drop = FALSE], Y[observed, j],
                             L = 1, max_iter = 100)$pip
        keep_y[j] <- any(pip > cutoff)
      }
      if (!any(keep_y)) {
        message("QC track: skipping individual context ", context,
                ". No outcomes passed PIP threshold ", cutoff, ".")
        return(NULL)
      }
      Y <- Y[, keep_y, drop = FALSE]
    }
    message("QC track: retained ", ncol(X), " variants and ", ncol(Y),
            " outcome(s) for individual context ", .cb_default(context, "input"), ".")
    list(X = X, Y = Y, maf = maf, X_variance = X_variance)
  }

  if (is.list(X) && !is.matrix(X) && !is.data.frame(X)) {
    X <- .cb_as_named_list(X, "individual")
    Y <- .cb_as_named_list(Y, "individual")
    contexts <- intersect(names(X), names(Y))
    if (length(contexts) == 0) stop("No matched X/Y contexts for individual-level QC.")
    cutoffs <- .cb_named_cutoff(pip_cutoff_to_skip, contexts, "pip_cutoff_to_skip_ind")
    out <- list()
    for (context in contexts) {
      res <- qc_one(
        X[[context]], Y[[context]],
        maf = .cb_list_value(maf, context),
        X_variance = .cb_list_value(X_variance, context),
        context = context,
        pip_cutoff_to_skip = cutoffs[[context]]
      )
      if (!is.null(res)) out[[context]] <- res
    }
    return(out)
  }
  qc_one(
    X = X, Y = Y, maf = maf, X_variance = X_variance,
    pip_cutoff_to_skip = pip_cutoff_to_skip,
    context = context
  )
}


##### Generic ColocBoost helper functions #####

.cb_default <- function(x, y) if (is.null(x)) y else x

.cb_merge_args <- function(x, y) {
  for (nm in names(y)) {
    if (!is.null(y[[nm]])) x[[nm]] <- y[[nm]]
  }
  x
}

.cb_as_named_list <- function(x, default_name) {
  if (is.null(x)) return(NULL)
  if (is.list(x) && !is.matrix(x) && !is.data.frame(x)) return(x)
  stats::setNames(list(x), default_name)
}

.cb_list_value <- function(x, name, default = NULL) {
  if (is.null(x)) return(default)
  if (is.list(x) && !is.matrix(x) && !is.data.frame(x)) {
    if (name %in% names(x)) return(x[[name]])
    return(default)
  }
  x
}

.cb_named_cutoff <- function(x, names_to_fill, arg_name) {
  if (length(x) == 1 && is.null(names(x))) {
    return(stats::setNames(rep(x, length(names_to_fill)), names_to_fill))
  }
  if (!is.null(names(x))) {
    missing <- setdiff(names_to_fill, names(x))
    if (length(missing) > 0) x <- c(x, stats::setNames(rep(0, length(missing)), missing))
    return(x[names_to_fill])
  }
  if (length(x) == length(names_to_fill)) {
    return(stats::setNames(x, names_to_fill))
  }
  stop(arg_name, " must be a scalar, named vector, or match the number of inputs.")
}

.cb_y_ncol <- function(y) {
  if (is.null(y)) return(0L)
  if (is.null(dim(y))) return(as.integer(length(y) > 0))
  ncol(y)
}

##### Outcome helper functions #####

.cb_colocboost_outcome_names <- function(args, prefer_supplied = TRUE) {
  if (isTRUE(prefer_supplied) && !is.null(args$outcome_names)) {
    return(as.character(args$outcome_names))
  }
  Y <- args$Y
  y_outcomes <- character()
  if (!is.null(Y)) {
    if (is.data.frame(Y)) Y <- as.matrix(Y)
    if (is.matrix(Y)) {
      y_outcomes <- colnames(Y)
      if (is.null(y_outcomes)) y_outcomes <- paste0("Y", seq_len(ncol(Y)))
    } else if (is.atomic(Y) && !is.list(Y)) {
      y_outcomes <- "Y1"
    } else {
      y_names <- names(Y)
      y_cols <- unlist(Map(function(y, nm) {
        if (is.null(y)) return(character())
        if (is.data.frame(y)) y <- as.matrix(y)
        if (!is.null(dim(y))) {
          cn <- colnames(y)
          if (is.null(cn) || any(is.na(cn) | cn == "")) {
            if (!is.null(nm) && !is.na(nm) && nm != "" && ncol(y) == 1) return(nm)
            return(paste0(.cb_default(nm, "Y"), "_", seq_len(ncol(y))))
          }
          return(as.character(cn))
        }
        if (!is.null(nm) && !is.na(nm) && nm != "") return(nm)
        "Y"
      }, Y, .cb_default(y_names, rep("", length(Y)))), use.names = FALSE)
      if (length(y_cols) != length(Y) || anyDuplicated(y_cols) == 0) {
        y_outcomes <- y_cols
      } else {
        if (is.null(y_names) || any(is.na(y_names) | y_names == "")) {
          y_names <- paste0("Y", seq_along(Y))
        }
        y_outcomes <- y_names
      }
    }
  }

  sumstat <- args$sumstat
  effect_est <- args$effect_est
  sumstat_outcomes <- character()
  if (!is.null(sumstat)) {
    if (is.data.frame(sumstat)) {
      sumstat_outcomes <- "sumstat1"
    } else {
      ss_names <- names(sumstat)
      if (is.null(ss_names) || any(is.na(ss_names) | ss_names == "")) {
        ss_names <- paste0("sumstat", seq_along(sumstat))
      }
      sumstat_outcomes <- as.character(ss_names)
    }
  } else if (!is.null(effect_est)) {
    effect_est <- as.matrix(effect_est)
    effect_names <- colnames(effect_est)
    if (is.null(effect_names)) effect_names <- paste0("sumstat", seq_len(ncol(effect_est)))
    sumstat_outcomes <- as.character(effect_names)
  }
  c(as.character(y_outcomes), sumstat_outcomes)
}

.cb_resolve_qc_outcome_names <- function(pre_qc_data_outcomes,
                                         pre_qc_display_outcomes,
                                         post_qc_data_outcomes) {
  if (length(post_qc_data_outcomes) == 0) return(character())
  if (length(pre_qc_data_outcomes) == length(pre_qc_display_outcomes) &&
      length(pre_qc_data_outcomes) > 0 &&
      !anyDuplicated(pre_qc_data_outcomes)) {
    idx <- match(post_qc_data_outcomes, pre_qc_data_outcomes)
    if (all(!is.na(idx))) {
      return(pre_qc_display_outcomes[idx])
    }
  }
  post_qc_data_outcomes
}

.cb_remap_focal_outcome_idx <- function(focal_outcome_idx,
                                        pre_qc_data_outcomes,
                                        pre_qc_display_outcomes,
                                        post_qc_data_outcomes,
                                        post_qc_display_outcomes) {
  if (is.null(focal_outcome_idx)) return(NULL)
  if (length(focal_outcome_idx) != 1 || is.na(focal_outcome_idx)) {
    warning("focal_outcome_idx must be a single non-missing index. Passing it through unchanged.")
    return(focal_outcome_idx)
  }
  focal_outcome_idx <- as.integer(focal_outcome_idx)
  if (length(post_qc_data_outcomes) == 0) return(NULL)
  if (focal_outcome_idx < 1 ||
      focal_outcome_idx > max(length(pre_qc_display_outcomes), length(pre_qc_data_outcomes))) {
    warning("focal_outcome_idx is outside the pre-QC outcome range. Passing it through unchanged.")
    return(focal_outcome_idx)
  }

  focal_display <- if (focal_outcome_idx <= length(pre_qc_display_outcomes)) {
    pre_qc_display_outcomes[[focal_outcome_idx]]
  } else {
    NULL
  }
  focal_data <- if (focal_outcome_idx <= length(pre_qc_data_outcomes)) {
    pre_qc_data_outcomes[[focal_outcome_idx]]
  } else {
    NULL
  }

  candidates <- integer()
  for (needle in unique(c(focal_display, focal_data))) {
    if (is.null(needle) || is.na(needle) || needle == "") next
    candidates <- c(
      candidates,
      which(post_qc_display_outcomes == needle),
      which(post_qc_data_outcomes == needle)
    )
    suffix_pattern <- paste0("_", needle)
    candidates <- c(
      candidates,
      which(endsWith(post_qc_display_outcomes, suffix_pattern)),
      which(endsWith(post_qc_data_outcomes, suffix_pattern))
    )
  }
  candidates <- unique(candidates)
  if (length(candidates) > 0) {
    return(candidates[[1]])
  }
  if (length(pre_qc_data_outcomes) == length(post_qc_data_outcomes) ||
      length(pre_qc_display_outcomes) == length(post_qc_display_outcomes)) {
    return(focal_outcome_idx)
  }

  focal_label <- .cb_default(focal_display, .cb_default(focal_data, focal_outcome_idx))
  warning("focal_outcome_idx refers to outcome ", focal_label,
          ", which is not present after QC. Setting focal_outcome_idx to NULL.")
  NULL
}

##### Individual-level ColocBoost helper functions #####

.cb_individual_qc_input_from_colocboost <- function(X, Y, dict_YX = NULL) {
  bind_y_for_qc <- function(Y_list, y_names = NULL) {
    if (is.null(y_names)) y_names <- names(Y_list)
    if (is.null(y_names)) y_names <- rep("", length(Y_list))
    mats <- Map(function(y, nm) {
      if (is.null(dim(y))) y <- matrix(y, ncol = 1)
      if (is.null(colnames(y))) {
        colnames(y) <- if (!is.null(nm) && !is.na(nm) && nm != "") nm else paste0("Y", seq_len(ncol(y)))
      }
      y
    }, Y_list, y_names)
    do.call(cbind, mats)
  }

  X_list <- .cb_as_named_list(X, "X1")
  Y_list <- .cb_as_named_list(Y, "Y1")
  matched <- intersect(names(X_list), names(Y_list))
  if (length(matched) > 0 && is.null(dict_YX)) {
    return(list(X = X_list[matched], Y = Y_list[matched]))
  }

  if (!is.null(dict_YX)) {
    dict <- as.matrix(dict_YX)
    if (ncol(dict) < 2) stop("dict_YX must have at least two columns.")
    X_qc <- list()
    Y_qc <- list()
    for (x_idx in unique(dict[, 2])) {
      if (is.na(x_idx) || x_idx < 1 || x_idx > length(X_list)) next
      y_idx <- dict[dict[, 2] == x_idx, 1]
      y_idx <- y_idx[!is.na(y_idx) & y_idx >= 1 & y_idx <= length(Y_list)]
      if (length(y_idx) == 0) next
      context <- names(X_list)[x_idx]
      if (is.null(context) || is.na(context) || context == "") context <- paste0("X", x_idx)
      X_qc[[context]] <- X_list[[x_idx]]
      Y_qc[[context]] <- bind_y_for_qc(Y_list[y_idx], names(Y_list)[y_idx])
    }
    if (length(X_qc) > 0) return(list(X = X_qc, Y = Y_qc))
  }

  if (length(X_list) == 1 && length(Y_list) > 0) {
    context <- names(X_list)[1]
    if (is.null(context) || is.na(context) || context == "") context <- "X1"
    return(list(
      X = stats::setNames(list(X_list[[1]]), context),
      Y = stats::setNames(list(bind_y_for_qc(Y_list, names(Y_list))), context)
    ))
  }

  if (length(X_list) == length(Y_list)) {
    if (is.null(names(X_list))) names(X_list) <- paste0("X", seq_along(X_list))
    names(Y_list) <- names(X_list)
    return(list(X = X_list, Y = Y_list))
  }

  list(X = X_list, Y = Y_list)
}

.cb_format_individual <- function(ind) {
  if (length(ind) == 0) return(list())
  if (!is.null(ind$X) && !is.null(ind$Y)) {
    ind <- list(individual = ind)
  }
  ind <- Filter(function(x) {
    !is.null(x$X) && .cb_y_ncol(x$Y) > 0
  }, ind)
  if (length(ind) == 0) return(list())
  X <- lapply(ind, `[[`, "X")
  Y <- lapply(ind, `[[`, "Y")
  unique_X <- list()
  X_match <- integer(length(X))
  for (i in seq_along(X)) {
    matched <- names(unique_X)[vapply(unique_X, identical, logical(1), X[[i]])]
    if (length(matched) > 0) {
      X_match[[i]] <- match(matched[[1]], names(unique_X))
    } else {
      unique_X[[names(ind)[i]]] <- X[[i]]
      X_match[[i]] <- length(unique_X)
    }
  }
  Y_split <- list()
  dict <- matrix(nrow = 0, ncol = 2)
  y_index <- 0L
  all_outcome_names <- unlist(lapply(Y, colnames), use.names = FALSE)
  duplicated_outcome_names <- unique(all_outcome_names[
    duplicated(all_outcome_names) | duplicated(all_outcome_names, fromLast = TRUE)
  ])
  for (i in seq_along(Y)) {
    y <- Y[[i]]
    if (is.null(dim(y))) y <- matrix(y, ncol = 1)
    for (j in seq_len(ncol(y))) {
      y_index <- y_index + 1L
      outcome <- colnames(y)[j]
      if (is.null(outcome) || is.na(outcome) || outcome == "") {
        outcome <- paste0("outcome", y_index)
      }
      context <- names(ind)[i]
      if (outcome %in% duplicated_outcome_names &&
          !is.null(context) && !is.na(context) && context != "") {
        outcome <- paste0(context, "_", outcome)
      }
      if (outcome %in% names(Y_split)) {
        outcome <- make.unique(c(names(Y_split), outcome))[length(Y_split) + 1]
      }
      Y_split[[outcome]] <- y[, j, drop = FALSE]
      dict <- rbind(dict, c(y_index, X_match[[i]]))
    }
  }
  colnames(dict) <- c("Y", "X")
  list(X = unique_X, Y = Y_split, dict_YX = dict, outcome_names = names(Y_split))
}

##### Summary-statistic ColocBoost helper functions #####

.cb_sumstat_qc_input_from_colocboost <- function(sumstat, LD, X_ref, dict_sumstatLD,
                                                LD_reference_info = NULL,
                                                variant_convention = c("A2_A1", "A1_A2")) {
  is_ld_data <- function(x) {
    methods::is(x, "LDData") || (is.list(x) && !is.null(x$LD_matrix))
  }
  as_reference_info_list <- function(x) {
    if (is.null(x)) return(NULL)
    if (is_ld_data(x) || is.data.frame(x) || is.character(x)) {
      return(list(LD = x))
    }
    if (is.list(x) && length(x) > 0) return(x)
    stop("LD_reference_info must be a .bim/.pvar path, data.frame, load_LD_matrix() result, or a list of these.")
  }
  reference_info_for_index <- function(reference_info, index) {
    if (is.null(reference_info)) return(NULL)
    if (length(reference_info) == 1) return(reference_info[[1]])
    reference_info[[min(index, length(reference_info))]]
  }
  validate_sumstat_for_qc <- function(sumstat, study) {
    if (is.null(sumstat) || !is.data.frame(sumstat)) {
      return(paste0("Summary-statistic QC for study ", study,
                    " requires each sumstat input to be a data.frame."))
    }
    missing_cols <- setdiff(c("variant", "z"), colnames(sumstat))
    if (length(missing_cols) > 0) {
      return(paste0("Summary-statistic QC for study ", study,
                    " requires sumstat columns: ", paste(missing_cols, collapse = ", "), "."))
    }
    NULL
  }
  ld_variant_names <- function(ld) {
    if (!is.matrix(ld)) return(NULL)
    if (nrow(ld) == ncol(ld)) rownames(ld) else colnames(ld)
  }

  variant_convention <- match.arg(variant_convention)
  sumstat <- .cb_as_named_list(sumstat, "sumstat")
  using_x_ref <- is.null(LD) && !is.null(X_ref)
  ld_source <- .cb_as_named_list(if (!is.null(LD)) LD else X_ref, "LD")
  reference_info <- as_reference_info_list(LD_reference_info)
  if (is.null(dict_sumstatLD)) {
    dict_sumstatLD <- cbind(seq_along(sumstat), pmin(seq_along(sumstat), length(ld_source)))
  }
  rss_input <- list()
  LD_data <- list()
  skip_reasons <- character()
  for (i in seq_along(sumstat)) {
    study <- names(sumstat)[i]
    ld_index <- dict_sumstatLD[i, 2]
    validation_reason <- validate_sumstat_for_qc(sumstat[[i]], study)
    if (!is.null(validation_reason)) {
      skip_reasons <- c(skip_reasons, validation_reason)
      next
    }
    ld_mat <- ld_source[[ld_index]]
    ref_info <- reference_info_for_index(reference_info, ld_index)
    if (is.null(ref_info)) {
      message("QC track: checking LD/X_ref variant names for summary-stat study ", study, ".")
      ref_ids <- ld_variant_names(ld_mat)
      if (is.null(ref_ids)) {
        skip_reasons <- c(
          skip_reasons,
          paste0("Summary-statistic QC for study ", study,
                 " requires LD row/column names or X_ref column names parseable as genomic variant IDs; ",
                 "provide LD_reference_info for QC. Skipping summary-statistic QC for this study.")
        )
        next
      }
      ld_data <- .cb_make_ld_data(
        ld_mat,
        is_genotype = using_x_ref,
        variant_convention = variant_convention
      )
      if (is.null(ld_data)) {
        skip_reasons <- c(
          skip_reasons,
          paste0("Summary-statistic QC for study ", study,
                 " could not parse LD/X_ref names as genomic variant IDs; ",
                 "provide LD_reference_info for QC. Skipping summary-statistic QC for this study.")
        )
        next
      }
      message("QC track: LD/X_ref names are parseable for summary-stat study ", study, ".")
    } else if (is_ld_data(ref_info)) {
      message("QC track: using supplied LD_reference_info LD data for summary-stat study ", study, ".")
      ld_data <- .normalize_ld_data_for_qc(ref_info)
    } else {
      message("QC track: using supplied LD_reference_info variant metadata for summary-stat study ", study, ".")
      ld_data <- .cb_make_ld_data(
        ld_mat,
        is_genotype = using_x_ref,
        reference_info = ref_info,
        variant_convention = variant_convention
      )
    }
    parsed <- tryCatch(
      .cb_parse_variants_for_qc(sumstat[[i]]$variant, variant_convention),
      error = function(e) {
        skip_reasons <<- c(
          skip_reasons,
          paste0("Summary-statistic QC for study ", study,
                 " could not parse sumstat$variant as genomic variant IDs: ",
                 conditionMessage(e), " Skipping summary-statistic QC for this study.")
        )
        NULL
      }
    )
    if (is.null(parsed)) next
    variant_id <- format_variant_id(parsed$chrom, parsed$pos, parsed$A2, parsed$A1)
    ss <- data.frame(parsed, z = sumstat[[i]]$z,
                     variant_id = variant_id,
                     stringsAsFactors = FALSE)
    n <- if ("n" %in% colnames(sumstat[[i]])) unique(sumstat[[i]]$n)[1] else NULL
    var_y <- if ("var_y" %in% colnames(sumstat[[i]])) unique(sumstat[[i]]$var_y)[1] else 1
    rss_input[[study]] <- list(sumstats = ss, n = n, var_y = var_y)
    LD_data[[study]] <- ld_data
  }
  list(rss_input = rss_input, LD_data = LD_data, skip_reasons = skip_reasons)
}

.cb_format_sumstat <- function(sumstat_qc) {
  valid_sumstat_entry <- function(ss_df, min_variants = 2) {
    if (is.null(ss_df) || !is.data.frame(ss_df)) return(FALSE)
    if (nrow(ss_df) < min_variants) return(FALSE)
    if (all(is.na(ss_df$z))) return(FALSE)
    if (all(ss_df$n <= 0 | is.na(ss_df$n))) return(FALSE)
    TRUE
  }
  filter_valid_sumstats <- function(sumstats, LD_mat, min_variants = 2) {
    dedupe_ld <- function(LD_mat, studies) {
      unique_ld <- list()
      LD_match <- character()
      for (study in studies) {
        ld <- LD_mat[[study]]
        matched <- names(unique_ld)[vapply(unique_ld, identical, logical(1), ld)]
        if (length(matched) > 0) {
          LD_match <- c(LD_match, matched[[1]])
        } else {
          unique_ld[[study]] <- ld
          LD_match <- c(LD_match, study)
        }
      }
      list(LD_mat = unique_ld, LD_match = LD_match)
    }

    valid_idx <- vapply(sumstats, valid_sumstat_entry, logical(1),
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
    LD_mat <- LD_mat[valid_idx]
    deduped <- dedupe_ld(LD_mat, names(sumstats))
    LD_mat <- deduped$LD_mat
    LD_match <- deduped$LD_match
    dict_sumstatLD <- cbind(seq_along(sumstats), match(LD_match, names(LD_mat)))
    list(sumstats = sumstats, LD_mat = LD_mat, LD_match = LD_match,
         dict_sumstatLD = dict_sumstatLD)
  }
  if (length(sumstat_qc) == 0) return(list())
  if (!is.null(sumstat_qc$rss_input) && !is.null(sumstat_qc$LD_matrix)) {
    sumstat_qc <- list(sumstat = sumstat_qc)
  }
  sumstat <- lapply(sumstat_qc, function(x) {
    ss <- x$rss_input$sumstats
    variant_id <- if ("variant_id" %in% colnames(ss)) {
      ss$variant_id
    } else {
      format_variant_id(ss$chrom, ss$pos, ss$A2, ss$A1)
    }
    data.frame(z = ss$z, n = x$rss_input$n,
               variant = normalize_variant_id(variant_id),
               stringsAsFactors = FALSE)
  })
  LD_mat <- lapply(sumstat_qc, `[[`, "LD_matrix")
  filtered <- filter_valid_sumstats(sumstat, LD_mat)
  if (is.null(filtered)) return(list())
  c(
    list(
      sumstat = filtered$sumstats,
      dict_sumstatLD = filtered$dict_sumstatLD
    ),
    build_ld_args(filtered$LD_mat)
  )
}


##### LD/reference QC helper functions #####

.cb_make_ld_data <- function(ld, is_genotype = NULL,
                             reference_info = NULL,
                             variant_convention = c("A2_A1", "A1_A2")) {
  reference_info_to_ref_panel <- function(reference_info) {
    if (is.character(reference_info) && length(reference_info) == 1) {
      if (!file.exists(reference_info)) {
        stop("LD_reference_info file does not exist: ", reference_info)
      }
      reference_info <- read_variant_metadata(reference_info)
    }
    if (!is.data.frame(reference_info)) {
      stop("LD_reference_info must be a .bim/.pvar path or data.frame when it is not a load_LD_matrix() result.")
    }
    info <- as.data.frame(reference_info, stringsAsFactors = FALSE)
    names(info) <- sub("^#", "", names(info))
    names(info) <- sub("^chr$", "chrom", names(info), ignore.case = TRUE)
    names(info) <- sub("^CHROM$", "chrom", names(info))
    names(info) <- sub("^POS$", "pos", names(info))
    names(info) <- sub("^ID$", "id", names(info))
    names(info) <- sub("^REF$", "A2", names(info))
    names(info) <- sub("^ALT$", "A1", names(info))
    names(info) <- sub("^a0$", "A2", names(info))
    names(info) <- sub("^a1$", "A1", names(info))
    names(info) <- sub("^rsid$", "id", names(info), ignore.case = TRUE)

    if ("variants" %in% names(info) && !"variant_id" %in% names(info)) {
      info$variant_id <- normalize_variant_id(as.character(info$variants))
    }
    if (!"variant_id" %in% names(info)) {
      missing_cols <- setdiff(c("chrom", "pos", "A2", "A1"), names(info))
      if (length(missing_cols) > 0) {
        stop("LD_reference_info must contain variant_id or columns chrom, pos, A2, A1. Missing: ",
             paste(missing_cols, collapse = ", "), ".")
      }
      info$variant_id <- normalize_variant_id(format_variant_id(info$chrom, info$pos, info$A2, info$A1))
    } else {
      info$variant_id <- normalize_variant_id(as.character(info$variant_id))
    }
    if (!all(c("chrom", "pos", "A2", "A1") %in% names(info))) {
      parsed <- parse_variant_id(info$variant_id)
      info$chrom <- parsed$chrom
      info$pos <- parsed$pos
      info$A2 <- parsed$A2
      info$A1 <- parsed$A1
    }
    keep_cols <- intersect(c("chrom", "id", "pos", "A2", "A1", "variant_id",
                             "allele_freq", "variance", "n_nomiss"), names(info))
    info[, keep_cols, drop = FALSE]
  }

  align_ref_panel_to_ld <- function(ref_panel, ld_names, n_variants) {
    if (!is.null(ld_names)) {
      match_idx <- rep(NA_integer_, length(ld_names))
      if ("id" %in% names(ref_panel)) {
        match_idx <- match(ld_names, ref_panel$id)
      }
      if (any(is.na(match_idx))) {
        variant_match <- match(normalize_variant_id(ld_names), ref_panel$variant_id)
        match_idx[is.na(match_idx)] <- variant_match[is.na(match_idx)]
      }
      if (all(!is.na(match_idx))) {
        ref_panel <- ref_panel[match_idx, , drop = FALSE]
        rownames(ref_panel) <- NULL
        return(ref_panel)
      }
      if (length(ld_names) == nrow(ref_panel)) {
        message("QC track: LD_reference_info could not be matched by LD names; using LD_reference_info row order.")
        return(ref_panel[seq_along(ld_names), , drop = FALSE])
      }
      stop("LD_reference_info could not be matched to LD/X_ref names. ",
           "Provide an id/variant_id column matching LD names, or provide rows in LD matrix order.")
    }
    if (nrow(ref_panel) != n_variants) {
      stop("LD_reference_info has ", nrow(ref_panel), " variants but LD/X_ref has ",
           n_variants, " columns. Provide LD_reference_info in LD matrix order or with matching LD names.")
    }
    message("QC track: LD/X_ref names are missing; using LD_reference_info row order.")
    ref_panel
  }

  variant_convention <- match.arg(variant_convention)
  if (is.null(is_genotype)) is_genotype <- is.matrix(ld) && nrow(ld) != ncol(ld)
  ref_ids <- if (is.matrix(ld) && nrow(ld) == ncol(ld)) rownames(ld) else colnames(ld)

  if (!is.null(reference_info)) {
    ref_panel <- reference_info_to_ref_panel(reference_info)
    ref_panel <- align_ref_panel_to_ld(ref_panel, ref_ids, ncol(ld))
    variant_ids <- ref_panel$variant_id
    if (is.matrix(ld) && nrow(ld) == ncol(ld)) {
      rownames(ld) <- colnames(ld) <- variant_ids
    } else if (is.matrix(ld)) {
      colnames(ld) <- variant_ids
    }
    return(list(
      LD_matrix = ld,
      LD_variants = variant_ids,
      ref_panel = ref_panel,
      block_metadata = if (!isTRUE(is_genotype)) .infer_single_ld_block_metadata(ref_panel) else NULL,
      is_genotype = isTRUE(is_genotype)
    ))
  }

  if (is.null(ref_ids)) {
    return(NULL)
  }

  parsed <- if (!is.null(ref_ids)) {
    tryCatch(
      .cb_parse_variants_for_qc(ref_ids, variant_convention),
      error = function(e) NULL
    )
  } else {
    NULL
  }
  if (is.null(parsed)) return(NULL)
  variant_ids <- if (!is.null(parsed)) format_variant_id(parsed$chrom, parsed$pos, parsed$A2, parsed$A1) else NULL
  if (!is.null(variant_ids)) {
    if (is.matrix(ld) && nrow(ld) == ncol(ld)) {
      rownames(ld) <- colnames(ld) <- variant_ids
    } else if (is.matrix(ld)) {
      colnames(ld) <- variant_ids
    }
    parsed$variant_id <- variant_ids
  }
  list(
    LD_matrix = ld,
    LD_variants = variant_ids,
    ref_panel = parsed,
    block_metadata = if (!isTRUE(is_genotype) && !is.null(parsed)) .infer_single_ld_block_metadata(parsed) else NULL,
    is_genotype = isTRUE(is_genotype)
  )
}

.cb_parse_variants_for_qc <- function(ids, variant_convention = c("A2_A1", "A1_A2")) {
  variant_convention <- match.arg(variant_convention)
  parsed <- parse_variant_id(ids)
  if (any(is.na(parsed$chrom)) || any(is.na(parsed$pos)) ||
      any(is.na(parsed$A1)) || any(is.na(parsed$A2))) {
    stop("QC requires variant IDs parseable as chr:pos:A2:A1 by default. ",
         "If the input uses chr:pos:A1:A2, set variant_convention = 'A1_A2'.")
  }
  if (identical(variant_convention, "A1_A2")) {
    parsed <- data.frame(chrom = parsed$chrom, pos = parsed$pos,
                         A2 = parsed$A1, A1 = parsed$A2,
                         stringsAsFactors = FALSE)
  }
  parsed
}
