#' This function loads a mixture data sets for a specific region, including individual-level data (genotype, phenotype, covariate data)
#' or summary statistics (sumstats, LD). Run \code{load_regional_univariate_data} and \code{load_rss_data} multiple times for different datasets
#'
#' @section Loading individual level data from multiple corhorts
#' @param region A string of chr:start-end for the phenotype region.
#' @param genotype_list a vector of PLINK bed file containing genotype data.
#' @param phenotype_list A vector of phenotype file names.
#' @param covariate_list A vector of covariate file names corresponding to the phenotype file vector.
#' @param conditions_list_individual A vector of strings representing different conditions or groups.
#' @param match_geno_pheno A vector of index of phentoypes matched to genotype if mulitple genotype PLINK files
#' @param maf_cutoff Minimum minor allele frequency (MAF) cutoff. Default is 0.
#' @param mac_cutoff Minimum minor allele count (MAC) cutoff. Default is 0.
#' @param xvar_cutoff Minimum variance cutoff. Default is 0.
#' @param imiss_cutoff Maximum individual missingness cutoff. Default is 0.
#' @param association_window A string of chr:start-end for the association analysis window (cis or trans). If not provided, all genotype data will be loaded.
#' @param extract_region_name A list of vectors of strings (e.g., gene ID ENSG00000269699) to subset the information when there are multiple regions available. Default is NULL.
#' @param region_name_col Column name containing the region name. Default is NULL.
#' @param keep_indel Logical indicating whether to keep insertions/deletions (INDELs). Default is TRUE.
#' @param keep_samples A vector of sample names to keep. Default is NULL.
#' @param phenotype_header Number of rows to skip at the beginning of the transposed phenotype file (default is 4 for chr, start, end, and ID).
#' @param scale_residuals Logical indicating whether to scale residuals. Default is FALSE.
#' @param tabix_header Logical indicating whether the tabix file has a header. Default is TRUE.
#'
#' @section Loading summary statistics from multiple corhorts or data set
#' @param sumstat_path_list A vector of file path to the summary statistics.
#' @param column_file_path_list A vector of file path to the column file for mapping.
#' @param LD_meta_file_path_list A vector of path of LD_metadata, LD_metadata is a data frame specifying LD blocks with columns "chrom", "start", "end", and "path". "start" and "end" denote the positions of LD blocks. "path" is the path of each LD block, optionally including bim file paths.
#' @param match_LD_sumstat A vector of index of sumstat matched to LD if mulitple sumstat files
#' @param conditions_list_sumstat A vector of strings representing different sumstats.
#' @param n_samples User-specified sample size. If unknown, set as 0 to retrieve from the sumstat file.
#' @param n_cases User-specified number of cases.
#' @param n_controls User-specified number of controls.
#' @param region The region where tabix use to subset the input dataset.
#' @param extract_sumstats_region_name User-specified gene/phenotype name used to further subset the phenotype data.
#' @param sumstats_region_name_col Filter this specific column for the extract_sumstats_region_name.
#' @param comment_string comment sign in the column_mapping file, default is #
#' @param extract_coordinates Optional data frame with columns "chrom" and "pos" for specific coordinates extraction.
#'
#' @return A list containing the individual_data and sumstat_data:
#' individual_data contains the following components if exist
#' \itemize{
#'   \item residual_Y: A list of residualized phenotype values (either a vector or a matrix).
#'   \item residual_X: A list of residualized genotype matrices for each condition.
#'   \item residual_Y_scalar: Scaling factor for residualized phenotype values.
#'   \item residual_X_scalar: Scaling factor for residualized genotype values.
#'   \item dropped_sample: A list of dropped samples for X, Y, and covariates.
#'   \item covar: Covariate data.
#'   \item Y: Original phenotype data.
#'   \item X_data: Original genotype data.
#'   \item X: Filtered genotype matrix.
#'   \item maf: Minor allele frequency (MAF) for each variant.
#'   \item chrom: Chromosome of the region.
#'   \item grange: Genomic range of the region (start and end positions).
#'   \item Y_coordinates: Phenotype coordinates if a region is specified.
#' }
#' sumstat_data contains the following components if exist
#' \itemize{
#'   \item sumstats: A list of summary statistics for the matched LD_info, each sublist contains sumstats, n, var_y from \code{load_rss_data}.
#'   \item LD_info: A list of LD information, each sublist contains combined_LD_variants, combined_LD_matrix, ref_panel  \code{load_LD_matrix}.
#' }
#'
#' @export
load_multitask_regional_data <- function(region, # a string of chr:start-end for phenotype region
                                         genotype_list = NULL, # PLINK file
                                         phenotype_list = NULL, # a vector of phenotype file names
                                         covariate_list = NULL, # a vector of covariate file names corresponding to the phenotype file vector
                                         conditions_list_individual = NULL, # a vector of strings
                                         match_geno_pheno = NULL, # a vector of index of phentoypes matched to genotype if mulitple genotype files
                                         maf_cutoff = 0,
                                         mac_cutoff = 0,
                                         xvar_cutoff = 0,
                                         imiss_cutoff = 0,
                                         association_window = NULL,
                                         extract_region_name = NULL,
                                         region_name_col = NULL,
                                         keep_indel = TRUE,
                                         keep_samples = NULL,
                                         keep_variants = NULL,
                                         phenotype_header = 4, # skip first 4 rows of transposed phenotype for chr, start, end and ID
                                         scale_residuals = FALSE,
                                         tabix_header = TRUE,
                                         # sumstat if need
                                         sumstat_path_list = NULL,
                                         column_file_path_list = NULL,
                                         LD_meta_file_path_list = NULL,
                                         match_LD_sumstat = NULL, # a vector of index of sumstat matched to LD if mulitple sumstat files
                                         conditions_list_sumstat = NULL,
                                         n_samples = 0,
                                         n_cases = 0,
                                         n_controls = 0,
                                         extract_sumstats_region_name = NULL,
                                         sumstats_region_name_col = NULL,
                                         comment_string = "#",
                                         extract_coordinates = NULL) {
  if (is.null(genotype_list) & is.null(sumstat_path_list)) {
    stop("Data load error. Please make sure at least one data set (sumstat_path_list or genotype_list) exists.")
  }

  # - if exist individual level data
  individual_data <- NULL
  if (!is.null(genotype_list)) {
    #### FIXME: later if we have mulitple genotype list
    if (length(genotype_list) != 1 & is.null(match_geno_pheno)) {
      stop("Data load error. Please make sure 'match_geno_pheno' exists if you load data from multiple individual-level data.")
    } else if (length(genotype_list) == 1 & is.null(match_geno_pheno)) {
      match_geno_pheno <- rep(1, length(phenotype_list))
    }

    # - load individual data from multiple datasets
    n_dataset <- unique(match_geno_pheno) ### FIXME
    for (i_data in 1:n_dataset) {
      # extract genotype file name
      genotype <- genotype_list[i_data]
      # extract phenotype and covariate file names
      pos <- which(match_geno_pheno == i_data)
      phenotype <- phenotype_list[pos]
      covariate <- covariate_list[pos]
      conditions <- conditions_list_individual[pos]
      dat <- load_regional_univariate_data(
        genotype = genotype, phenotype = phenotype,
        covariate = covariate, region = region,
        association_window = association_window,
        conditions = conditions, xvar_cutoff = xvar_cutoff,
        maf_cutoff = maf_cutoff, mac_cutoff = mac_cutoff,
        imiss_cutoff = imiss_cutoff, keep_indel = keep_indel,
        keep_samples = keep_samples, keep_variants = keep_variants,
        extract_region_name = extract_region_name,
        phenotype_header = phenotype_header,
        region_name_col = region_name_col,
        scale_residuals = scale_residuals
      )
      if (is.null(individual_data)) {
        individual_data <- dat
      } else {
        individual_data <- lapply(names(dat), function(k) {
          c(individual_data[[k]], dat[[k]])
        })
        individual_data$chrom <- dat$chrom
        individual_data$grange <- dat$grange
      }
    }
  }

  # - if exist summstat data
  sumstat_data <- NULL
  if (!is.null(sumstat_path_list)) {
    if (length(match_LD_sumstat) == 0) {
      match_LD_sumstat[[1]] <- conditions_list_sumstat
    }
    if (length(match_LD_sumstat) != length(LD_meta_file_path_list)) {
      stop("Please make sure 'match_LD_sumstat' matched 'LD_meta_file_path_list' if you load data from multiple sumstats.")
    }
    # - load sumstat data from multiple datasets
    n_LD <- length(match_LD_sumstat)
    for (i_ld in 1:n_LD) {
      # extract LD meta file path name
      LD_meta_file_path <- LD_meta_file_path_list[i_ld]
      LD_info <- load_LD_matrix(LD_meta_file_path,
        region = association_window,
        extract_coordinates = extract_coordinates
      )
      # extract sumstat information
      conditions <- match_LD_sumstat[[i_ld]]
      pos <- match(conditions, conditions_list_sumstat)
      sumstats <- lapply(pos, function(ii) {
        sumstat_path <- sumstat_path_list[ii]
        column_file_path <- column_file_path_list[ii]
        # FIXME later: when consider multiple LD reference
        tmp <- load_rss_data(
          sumstat_path = sumstat_path, column_file_path = column_file_path,
          n_sample = n_samples[ii], n_case = n_cases[ii], n_control = n_controls[ii],
          region = association_window, extract_region_name = extract_sumstats_region_name,
          region_name_col = sumstats_region_name_col, comment_string = comment_string
        )
        if (!("variant_id" %in% colnames(tmp$sumstats))) {
          tmp$sumstats <- tmp$sumstats %>%
            rowwise() %>%
            mutate(variant_id = paste0(c(chrom, pos, A1, A2), collapse = ":"))
        }
        return(tmp)
      })
      names(sumstats) <- conditions
      sumstat_data$sumstats <- c(sumstat_data$sumstats, list(sumstats))
      sumstat_data$LD_info <- c(sumstat_data$LD_info, list(LD_info))
    }
    names(sumstat_data$sumstats) <- names(sumstat_data$LD_info) <- names(match_LD_sumstat)
  }

  return(list(
    individual_data = individual_data,
    sumstat_data = sumstat_data
  ))
}

#' Multi-trait colocalization analysis pipeline
#'
#' This function perform a multi-trait colocalization using ColocBoost
#'
#' @param region_data A region data loaded from \code{load_regional_data}.
#' @param target_trait Name of trait if perform targeted ColocBoost
#' @param event_filters A list of pattern for filtering events based on context names. Example: for sQTL, list(type_pattern = ".*clu_(\\d+_[+-?]).*",valid_pattern = "clu_(\\d+_[+-?]):PR:",exclude_pattern = "clu_(\\d+_[+-?]):IN:")
#' @param maf_cutoff A scalar to remove variants with maf < maf_cutoff, dafault is 0.005.
#' @param pip_cutoff_to_skip_ind A vector of cutoff values for skipping analysis based on PIP values for each context. Default is 0.
#' @param pip_cutoff_to_skip_sumstat A vector of cutoff values for skipping analysis based on PIP values for each sumstat Default is 0.
#' @param qc_method Quality control method to use. Options are "rss_qc", "dentist", or "slalom" (default: "rss_qc").
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
#'   \item LD_info: A list of LD information, each sublist contains combined_LD_variants, combined_LD_matrix, ref_panel  \code{load_LD_matrix}.
#' }
#'
#' @export
colocboost_analysis_pipeline <- function(region_data,
                                         target_trait = NULL,
                                         event_filters = NULL,
                                         # - analysis
                                         xqtl_coloc = TRUE,
                                         joint_gwas = FALSE,
                                         separate_gwas = FALSE,
                                         # - individual QC
                                         maf_cutoff = 0.0005,
                                         pip_cutoff_to_skip_ind = 0,
                                         # - sumstat QC
                                         remove_indels = FALSE,
                                         pip_cutoff_to_skip_sumstat = 0,
                                         qc_method = c("rss_qc", "dentist", "slalom"),
                                         impute = TRUE,
                                         impute_opts = list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5, lamb = 0.01),
                                         ...) {
  # - internel function by filtering events based on event_filters
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

  # - extract contexts and studies from region data
  extract_contexts_studies <- function(region_data, phenotypes_init = NULL) {
    individual_data <- region_data$individual_data
    sumstat_data <- region_data$sumstat_data

    if (is.null(phenotypes_init)) {
      # - inital setup
      phenotypes <- list("individual_contexts" = NULL, "sumstat_studies" = NULL)
      if (!is.null(individual_data)) {
        phenotypes$individual_contexts <- names(individual_data$residual_Y)
      } else {
        message("No individual data in this region!")
      }
      if (!is.null(sumstat_data)) {
        phenotypes$sumstat_studies <- sapply(sumstat_data$sumstats, function(ss) names(ss)) %>%
          unlist() %>%
          as.vector()
      } else {
        message("No sumstat data in this region!")
      }
    } else {
      # - after QC
      phenotypes <- list("individual_contexts" = NULL, "sumstat_studies" = NULL)
      if (!is.null(individual_data)) {
        null_Y <- which(sapply(individual_data$Y, is.null))
        if (length(null_Y) == 0) {
          message("All individual data pass QC steps.")
          phenotypes$individual_contexts <- names(individual_data$Y)
        } else if (length(null_Y) != length(individual_data$Y)) {
          message(paste(
            "Skipping follow-up analysis for individual traits",
            paste(names(individual_data$Y)[null_Y], collapse = ";"), "after QC."
          ))
          phenotypes$individual_contexts <- names(individual_data$Y)[-null_Y]
        } else if (length(null_Y) == length(individual_data$Y)) {
          message("No individual data pass QC.")
        }
      } else {
        message("No individual data pass QC.")
      }
      if (!is.null(sumstat_data)) {
        phenotypes$sumstat_studies <- names(sumstat_data$sumstats)
        sumstat_studies_init <- phenotypes_init$sumstat_studies
        if (length(sumstat_studies_init) == length(phenotypes$sumstat_studies)) {
          message("All sumstat studies pass QC steps.")
        } else {
          message(paste(
            "Skipping follow-up analysis for sumstat studies",
            paste(setdiff(sumstat_studies_init, phenotypes$sumstat_studies), collapse = ";"), "after QC."
          ))
        }
      } else {
        message("No sumstat data pass QC.")
      }
    }
    return(phenotypes)
  }

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
    Y <- lapply(1:length(Y), function(i) {
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
    remove_indels = remove_indels,
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
      Y <- lapply(1:length(Y), function(i) {
        y <- Y[[i]]
        lapply(seq_len(ncol(y)), function(j) y[, j, drop = FALSE] %>% setNames(colnames(y)[j]))
      })
      dict_YX <- cbind(seq_along(Reduce("c", Y)), rep(seq_along(Y), sapply(Y, length)))
      Y <- Reduce("c", Y)
      Y <- Y %>% setNames(sapply(Y, colnames))
    }
  } else {
    X <- Y <- dict_YX <- NULL
  }

  ####### ========= organize summary statistics ======== ########
  sumstat_data <- region_data$sumstat_data
  if (!is.null(sumstat_data$sumstats)) {
    sumstats <- lapply(sumstat_data$sumstats, function(ss) {
      z <- ss$sumstats$z
      variant <- paste0("chr", ss$sumstats$variant_id)
      n <- ss$n
      data.frame("z" = z, "n" = n, "variant" = variant)
    })
    names(sumstats) <- names(sumstat_data$sumstats)
    LD_mat <- lapply(sumstat_data$LD_mat, function(ld) {
      colnames(ld) <- rownames(ld) <- paste0("chr", colnames(ld))
      return(ld)
    })
    LD_match <- sumstat_data$LD_match
    dict_sumstatLD <- cbind(seq_along(sumstats), match(LD_match, names(sumstat_data$LD_mat)))
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
    t11 <- Sys.time()
    traits <- names(Y)
    target_idx <- NULL
    if (!is.null(target_trait)) {
      if (target_trait %in% traits) {
        target_idx <- which(traits == target_trait)
      }
    }
    res_xqtl <- colocboost(
      X = X, Y = Y, dict_YX = dict_YX,
      outcome_names = traits, target_idx = target_idx, ...
    )
    t12 <- Sys.time()
    analysis_results$xqtl_coloc <- res_xqtl
    analysis_results$computing_time$Analysis$xqtl_coloc <- t12 - t11
  }
  # - run joint GWAS no targeted version of ColocBoost
  if (joint_gwas & !is.null(sumstats)) {
    message(paste("====== Performing non-targeted version GWAS-xQTL ColocBoost on", length(Y), "contexts and", length(sumstats), "GWAS. ====="))
    t21 <- Sys.time()
    traits <- c(names(Y), names(sumstats))
    res_gwas <- colocboost(
      X = X, Y = Y, sumstat = sumstats, LD = LD_mat,
      dict_YX = dict_YX, dict_sumstatLD = dict_sumstatLD,
      outcome_names = traits, target_idx = NULL, ...
    )
    t22 <- Sys.time()
    analysis_results$joint_gwas <- res_gwas
    analysis_results$computing_time$Analysis$joint_gwas <- t22 - t21
  }
  # - run targeted version of ColocBoost for each GWAS
  if (separate_gwas & !is.null(sumstats)) {
    t31 <- Sys.time()
    res_gwas_separate <- analysis_results$separate_gwas
    for (i_gwas in 1:nrow(dict_sumstatLD)) {
      current_study <- names(sumstats)[i_gwas]
      message(paste("====== Performing targeted version GWAS-xQTL ColocBoost on", length(Y), "contexts and ", current_study, "GWAS. ====="))
      dict <- dict_sumstatLD[i_gwas, ]
      traits <- c(names(Y), current_study)
      res_gwas_separate[[current_study]] <- colocboost(
        X = X, Y = Y, sumstat = sumstats[dict[1]],
        LD = LD_mat[dict[2]], dict_YX = dict_YX,
        outcome_names = traits, target_idx = length(traits), ...
      )
    }
    t32 <- Sys.time()
    analysis_results$separate_gwas <- res_gwas_separate
    analysis_results$computing_time$Analysis$separate_gwas <- list("total" = t32 - t31, "n_studies" = nrow(dict_sumstatLD), "average" = (t32 - t31) / nrow(dict_sumstatLD))
  }

  return(analysis_results)
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
#' @param qc_method Quality control method to use. Options are "rss_qc", "dentist", or "slalom" (default: "rss_qc").
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
#'   \item LD_info: A list of LD information, each sublist contains combined_LD_variants, combined_LD_matrix, ref_panel  \code{load_LD_matrix}.
#' }
#'
#' @noRd
qc_regional_data <- function(region_data,
                             # - individual
                             maf_cutoff = 0.0005,
                             pip_cutoff_to_skip_ind = 0,
                             # - sumstat
                             remove_indels = FALSE,
                             pip_cutoff_to_skip_sumstat = 0,
                             qc_method = c("rss_qc", "dentist", "slalom"),
                             impute = TRUE,
                             impute_opts = list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5, lamb = 0.01)) {
  #### related internal functions
  # Add context names to colname of Y if missing
  add_context_to_Y <- function(res_Y) {
    res <- lapply(1:length(res_Y), function(iy) {
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
    n_context <- length(X)
    residual_X <- residual_Y <- list()
    keep_contexts <- c()
    for (i_context in 1:n_context) {
      resX <- X[[i_context]]
      resY <- Y[[i_context]]
      maf <- MAF[[i_context]]
      context <- names(Y)[i_context]
      if (is.null(resY)) next
      # - remove variants with maf < maf_cutoff
      # tmp <- filter_resX_maf(resX, maf, maf_cutoff = maf_cutoff)
      resX <- filter_X(resX, missing_rate_thresh = NULL, maf_thresh = maf_cutoff, maf = maf)
      # Initial PIP check
      resY <- filter_resY_pip(resX, resY, pip_cutoff = pip_cutoff_to_skip_ind[i_context], context = context)
      if (!is.null(resY)) {
        residual_X <- c(residual_X, list(resX))
        residual_Y <- c(residual_Y, list(resY))
        keep_contexts <- c(keep_contexts, context)
      }
    }
    if (length(keep_contexts) == 0) {
      message(paste("Skipping follow-up analysis for all contexts."))
      return(NULL)
    } else {
      message(paste("Region includes the following contexts after inital screening:", paste(keep_contexts, collapse = ";"), "."))
      names(residual_X) <- names(residual_Y) <- keep_contexts
      return(list(
        X = residual_X,
        Y = residual_Y
      ))
    }
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
                                         remove_indels = FALSE,
                                         pip_cutoff_to_skip_sumstat = 0,
                                         qc_method = c("rss_qc", "dentist", "slalom"),
                                         impute = TRUE,
                                         impute_opts = list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5, lamb = 0.01)) {
    n_LD <- length(sumstat_data$LD_info)
    final_sumstats <- final_LD <- NULL
    LD_match <- c()
    for (i in 1:n_LD) {
      LD_data <- sumstat_data$LD_info[[i]]
      sumstats <- sumstat_data$sumstats[[i]]
      for (ii in 1:length(sumstats)) {
        sumstat <- sumstats[[ii]]
        n <- sumstat$n
        var_y <- sumstat$var_y
        conditions_sumstat <- names(sumstats)[ii]
        pip_cutoff_to_skip_ld <- pip_cutoff_to_skip_sumstat[conditions_sumstat] %>% as.numeric()

        # Preprocess the input data
        preprocess_results <- rss_basic_qc(sumstat$sumstats, LD_data, remove_indels = remove_indels)
        sumstat$sumstats <- preprocess_results$sumstats
        LD_mat <- preprocess_results$LD_mat

        # initial PIP checking
        if (pip_cutoff_to_skip_ld != 0) {
          pip <- susie_rss_wrapper(z = sumstat$sumstats$z, R = LD_mat, L = 1, n = n, var_y = var_y)$pip
          if (pip_cutoff_to_skip_ld < 0) {
            # automatically determine the cutoff to use
            pip_cutoff_to_skip_ld <- 3 * 1 / nrow(LD_mat)
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

        # Perform quality control - remove
        if (!is.null(qc_method)) {
          qc_results <- summary_stats_qc(sumstat$sumstats, LD_data, n = n, var_y = var_y, method = qc_method)
          sumstat$sumstats <- qc_results$sumstats
          LD_mat <- qc_results$LD_mat
        }
        # Perform imputation
        if (impute) {
          LD_matrix <- partition_LD_matrix(LD_data)
          impute_results <- raiss(LD_data$ref_panel, sumstat$sumstats, LD_matrix,
            rcond = impute_opts$rcond,
            R2_threshold = impute_opts$R2_threshold, minimum_ld = impute_opts$minimum_ld, lamb = impute_opts$lamb
          )
          sumstat$sumstats <- impute_results$result_filter
          LD_mat <- impute_results$LD_mat
        }

        # - check if LD exist
        if (length(final_LD) == 0) {
          final_LD <- c(final_LD, list(LD_mat) %>% setNames(conditions_sumstat))
          final_sumstats <- c(final_sumstats, list(sumstat) %>% setNames(conditions_sumstat))
          LD_match <- c(LD_match, conditions_sumstat)
        } else {
          variants <- colnames(LD_mat)
          final_sumstats <- c(final_sumstats, list(sumstat) %>% setNames(conditions_sumstat))
          exist_variants <- lapply(final_LD, colnames)
          if_exist <- sapply(exist_variants, function(v) all(variants == v))
          pos <- which(if_exist)
          if (length(pos) == 0) {
            final_LD <- c(final_LD, list(LD_mat) %>% setNames(conditions_sumstat))
            LD_match <- c(LD_match, conditions_sumstat)
          } else {
            LD_match <- c(LD_match, LD_match[pos[1]])
          }
        }
      }
    }
    return(list(sumstats = final_sumstats, LD_mat = final_LD, LD_match = LD_match))
  }


  # - summary statistics QC
  sumstat_data <- region_data$sumstat_data
  if (!is.null(sumstat_data)) {
    # - initial check PIP, qc or impute
    sumstat_data <- summary_stats_qc_multitask(sumstat_data,
      remove_indels = remove_indels,
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
