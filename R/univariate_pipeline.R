#' Univariate Analysis Pipeline
#'
#' This function performs univariate analysis for fine-mapping and Transcriptome-Wide Association Study (TWAS)
#' with optional cross-validation.
#'
#' @param X A matrix of genotype data where rows represent samples and columns represent genetic variants.
#' @param Y A vector of phenotype measurements.
#' @param X_scalar A scalar or vector to rescale X to its original scale.
#' @param Y_scalar A scalar to rescale Y to its original scale.
#' @param maf A vector of minor allele frequencies for each variant in X.
#' @param X_variance Optional variance of X. Default is NULL.
#' @param other_quantities A list of other quantities to be passed to susie_post_processor. Default is an empty list.
#' @param imiss_cutoff Individual missingness cutoff. Default is 1.0.
#' @param maf_cutoff Minor allele frequency cutoff. Default is NULL.
#' @param xvar_cutoff Variance cutoff for X. Default is 0.05.
#' @param ld_reference_meta_file An optional path to a file containing linkage disequilibrium reference data. Default is NULL.
#' @param pip_cutoff_to_skip Cutoff value for skipping analysis based on PIP values. Default is 0.
#' @param init_L Initial number of components for SuSiE model optimization. Default is 5.
#' @param max_L The maximum number of components in SuSiE. Default is 20.
#' @param l_step Step size for increasing the number of components during SuSiE optimization. Default is 5.
#' @param signal_cutoff Cutoff value for signal identification in PIP values. Default is 0.025.
#' @param coverage A vector of coverage probabilities for credible sets. Default is c(0.95, 0.7, 0.5).
#' @param min_abs_corr Minimum absolute correlation for credible set purity filtering. Default is 0.8,
#'   which is stricter than the susieR default of 0.5.
#' @param twas_weights Whether to compute TWAS weights. Default is TRUE.
#' @param sample_partition Sample partition for cross-validation. Default is NULL.
#' @param max_cv_variants The maximum number of variants to be included in cross-validation. Default is -1 (no limit).
#' @param cv_folds The number of folds to use for cross-validation. Default is 5.
#' @param cv_threads The number of threads to use for parallel computation in cross-validation. Default is 1.
#' @param verbose Verbosity level. Default is 0.
#'
#' @return A list containing the univariate analysis results.
#' @importFrom susieR susie
#' @export
univariate_analysis_pipeline <- function(
    # input data
    X,
    Y,
    maf,
    X_scalar = 1,
    Y_scalar = 1,
    X_variance = NULL,
    other_quantities = list(),
    # filters
    imiss_cutoff = 1.0,
    maf_cutoff = NULL,
    xvar_cutoff = 0,
    ld_reference_meta_file = NULL,
    pip_cutoff_to_skip = 0,
    # methods parameter configuration
    init_L = 5,
    max_L = 20,
    l_step = 5,
    # fine-mapping results summary
    signal_cutoff = 0.025,
    coverage = c(0.95, 0.7, 0.5),
    min_abs_corr = 0.8,
    finemapping_extra_opts = list(refine = TRUE),
    # TWAS weights and CV for TWAS weights
    twas_weights = TRUE,
    sample_partition = NULL,
    max_cv_variants = -1,
    cv_folds = 5,
    cv_threads = 1,
    verbose = 0) {
  # Input validation
  if (!is.matrix(X) || !is.numeric(X)) stop("X must be a numeric matrix")
  if (!is.vector(Y) && !(is.matrix(Y) && ncol(Y) == 1) || !is.numeric(Y)) stop("Y must be a numeric vector or a single column matrix")
  if (nrow(X) != length(Y)) stop("X and Y must have the same number of rows/length")
  if (!is.numeric(maf) || length(maf) != ncol(X)) stop("maf must be a numeric vector with length equal to the number of columns in X")
  if (any(maf < 0 | maf > 1)) stop("maf values must be between 0 and 1")
  if (!is.numeric(X_scalar) || (length(X_scalar) != 1 && length(X_scalar) != ncol(X))) stop("X_scalar must be a numeric scalar or vector with length equal to the number of columns in X")
  if (!is.numeric(Y_scalar) || length(Y_scalar) != 1) stop("Y_scalar must be a numeric scalar")
  if (!is.numeric(init_L) || init_L <= 0) stop("init_L must be a positive integer")
  if (!is.numeric(max_L) || max_L <= 0) stop("max_L must be a positive integer")
  if (!is.numeric(l_step) || l_step <= 0) stop("l_step must be a positive integer")

  # Initial PIP check
  if (pip_cutoff_to_skip != 0) {
    if (pip_cutoff_to_skip < 0) {
      # automatically determine the cutoff to use
      pip_cutoff_to_skip <- 3 * 1 / ncol(X)
    }
    top_model_pip <- susie(X, Y, L = 1)$pip
    if (!any(top_model_pip > pip_cutoff_to_skip)) {
      message(paste("Skipping follow-up analysis: No signals above PIP threshold", pip_cutoff_to_skip, "in initial model screening."))
      return(list())
    } else {
      message(paste("Follow-up on region because signals above PIP threshold", pip_cutoff_to_skip, "were detected in initial model screening."))
    }
  }

  # Filter variants if LD reference is provided
  if (!is.null(ld_reference_meta_file)) {
    variants_kept <- filter_variants_by_ld_reference(colnames(X), ld_reference_meta_file)
    X <- X[, variants_kept$data, drop = FALSE]
    maf <- maf[variants_kept$idx]
    if (length(X_scalar) > 1) X_scalar <- X_scalar[variants_kept$idx]
  }

  # Filter X based on missingness, MAF, and variance
  if (!is.null(imiss_cutoff) || !is.null(maf_cutoff)) {
    X_filtered <- filter_X(X, imiss_cutoff, maf_cutoff, var_thresh = xvar_cutoff, maf = maf, X_variance = X_variance)
    kept_indices <- match(colnames(X_filtered), colnames(X))
    maf <- maf[kept_indices]
    if (length(X_scalar) > 1) X_scalar <- X_scalar[kept_indices]
    X <- X_filtered
  }

  # Main analysis
  st <- proc.time()
  res <- list()

  # SuSiE analysis with optimization
  message("Fitting SuSiE model on input data with L optimization...")
  base_susie_args <- list(X = X, y = Y, init_L = init_L, max_L = max_L, l_step = l_step, coverage = coverage[1])
  susie_args <- modifyList(finemapping_extra_opts, base_susie_args) # modifyList(A, B): B overrides A
  res$susie_fitted <- do.call(susie_wrapper, susie_args)

  # Process SuSiE results
  susie_result_trimmed <- susie_post_processor(
    res$susie_fitted, X, Y, X_scalar, Y_scalar, maf,
    secondary_coverage = if (length(coverage) > 1) coverage[-1] else NULL,
    signal_cutoff = signal_cutoff,
    min_abs_corr = min_abs_corr,
    other_quantities = other_quantities
  )
  res <- c(res, susie_result_trimmed)
  res$total_time_elapsed <- proc.time() - st

  # TWAS weights and cross-validation
  if (twas_weights) {
    res$twas_weights_result <- twas_weights_pipeline(
      X, Y, res$susie_fitted,
      cv_folds = cv_folds,
      max_cv_variants = max_cv_variants,
      cv_threads = cv_threads,
      sample_partition = sample_partition
    )
    if ("top_loci" %in% names(res)) res$twas_weights_result$susie_weights_intermediate$top_loci <- res$top_loci
  }

  return(res)
}

#' Load LD for a study, supporting single or mixture panels.
#'
#' @param ld_path A single LD metadata TSV path, or comma-separated paths for
#'   mixture panels (e.g., "ld_EUR.tsv,ld_AFR.tsv").
#' @param region Region string "chr:start-end".
#' @return An LD_data list from load_LD_matrix. For single panels, returns as-is.
#'   For mixture panels, LD_matrix is a list of X matrices (one per panel).
#' @export
load_study_LD <- function(ld_path, region) {
  paths <- strsplit(ld_path, ",")[[1]]
  if (length(paths) == 1) {
    return(load_LD_matrix(paths, region, return_genotype = "auto"))
  }
  # Mixture: load each panel as genotype X
  base <- load_LD_matrix(paths[1], region, return_genotype = TRUE)
  X_list <- c(
    list(base$LD_matrix),
    lapply(paths[-1], function(p) load_LD_matrix(p, region, return_genotype = TRUE)$LD_matrix)
  )
  base$LD_matrix <- X_list
  base
}

#' RSS Analysis Pipeline
#'
#' End-to-end pipeline for summary statistics fine-mapping via SuSiE RSS.
#' Supports both z+R (correlation matrix) and z+X (genotype matrix) interfaces.
#'
#' @param sumstat_path File path to the summary statistics.
#' @param column_file_path File path to the column mapping file.
#' @param LD_data A list from load_LD_matrix containing LD_matrix, LD_variants,
#'   ref_panel, block_metadata, and is_genotype flag. When is_genotype=TRUE
#'   (from return_genotype=TRUE), LD_matrix contains genotype X and susie_rss
#'   uses the z+X interface. R is computed internally for QC/imputation.
#' @param n_sample Sample size. If 0, retrieved from the sumstat file.
#' @param n_case Number of cases (for case-control studies).
#' @param n_control Number of controls (for case-control studies).
#' @param region Region string "chr:start-end" for tabix subsetting.
#' @param skip_region Character vector of regions to skip (format "chrom:start-end").
#' @param extract_region_name Gene/phenotype name to subset.
#' @param region_name_col Column to filter for extract_region_name.
#' @param qc_method QC method: "slalom" or "dentist".
#' @param finemapping_method One of "susie_rss", "single_effect", "bayesian_conditional_regression".
#' @param finemapping_opts List of fine-mapping options (init_L, max_L, l_step, coverage,
#'   signal_cutoff, min_abs_corr).
#' @param impute Whether to impute missing variants via RAISS (default TRUE).
#' @param impute_opts List of imputation options (rcond, R2_threshold, minimum_ld, lamb).
#' @param pip_cutoff_to_skip PIP threshold for early stopping (default 0, no skip).
#' @param stochastic_ld_sample Passed to susie_rss. NULL (default), TRUE, or integer.
#' @param keep_indel Whether to keep indel variants (default TRUE).
#' @param comment_string Comment character for sumstat file (default "#").
#' @param diagnostics Whether to include diagnostic info (default FALSE).
#'
#' @return A list with fine-mapping results and analyzed summary statistics.
#' @importFrom magrittr %>%
#' @export
rss_analysis_pipeline <- function(
    sumstat_path, column_file_path, LD_data,
    n_sample = 0, n_case = 0, n_control = 0, region = NULL, skip_region = NULL,
    extract_region_name = NULL, region_name_col = NULL,
    qc_method = c("slalom", "dentist"),
    finemapping_method = c("susie_rss", "single_effect", "bayesian_conditional_regression"),
    finemapping_opts = list(
      init_L = 5, max_L = 20, l_step = 5,
      coverage = c(0.95, 0.7, 0.5), signal_cutoff = 0.025,
      min_abs_corr = 0.8
    ),
    impute = TRUE, impute_opts = list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5, lamb = 0.01),
    pip_cutoff_to_skip = 0, stochastic_ld_sample = NULL,
    keep_indel = TRUE, comment_string = "#", diagnostics = FALSE) {
  # Detect genotype input: single X matrix or list of X matrices (mixture panel).
  # susie_rss accepts X=list(X1, X2, ...) for multi-panel mixture.
  is_X_list <- is.list(LD_data$LD_matrix) && !is.matrix(LD_data$LD_matrix)
  use_X <- isTRUE(LD_data$is_genotype) || is_X_list
  if (use_X) {
    X_data <- LD_data$LD_matrix
    # Compute R from first panel (or single panel) for QC/imputation
    X_for_R <- if (is_X_list) X_data[[1]] else X_data
    LD_data$LD_matrix <- compute_LD(X_for_R, method = "sample")
    LD_data$is_genotype <- FALSE
  }
  res <- list()
  rss_input <- load_rss_data(
    sumstat_path = sumstat_path, column_file_path = column_file_path,
    n_sample = n_sample, n_case = n_case, n_control = n_control,
    extract_region_name = extract_region_name, region = region,
    region_name_col = region_name_col, comment_string = comment_string
  )

  sumstats <- rss_input$sumstats
  n <- rss_input$n
  var_y <- rss_input$var_y

  if (nrow(sumstats) == 0) {
    return(list(rss_data_analyzed = sumstats))
  }

  # Preprocess: QC and imputation require LD_data with correlation matrix R.
  # When using X path, compute R from X for QC/imputation, then pass X to susie_rss.
  preprocess_results <- rss_basic_qc(sumstats, LD_data, skip_region = skip_region, keep_indel = keep_indel)
  sumstats <- preprocess_results$sumstats
  LD_mat <- preprocess_results$LD_mat

  if (nrow(sumstats) == 0) {
    message("No variants left after preprocessing. Returning empty results.")
    return(list(rss_data_analyzed = sumstats))
  }

  # PIP screening (always uses R)
  if (pip_cutoff_to_skip != 0) {
    if (pip_cutoff_to_skip < 0) pip_cutoff_to_skip <- 3 / nrow(sumstats)
    top_model_pip <- susie_rss_wrapper(z = sumstats$z, R = LD_mat, L = 1, n = n, var_y = var_y)$pip
    if (!any(top_model_pip > pip_cutoff_to_skip)) {
      message("Skipping follow-up analysis: No signals above PIP threshold ", pip_cutoff_to_skip)
      return(list(rss_data_analyzed = sumstats))
    }
    message("Follow-up on region: signals above PIP threshold ", pip_cutoff_to_skip, " detected.")
  }

  # Quality control (always uses R)
  if (!is.null(qc_method)) {
    qc_results <- summary_stats_qc(sumstats, LD_data, n = n, method = qc_method)
    sumstats <- qc_results$sumstats
    LD_mat <- qc_results$LD_mat
  }

  # Imputation (always uses R)
  if (impute) {
    LD_matrix <- partition_LD_matrix(LD_data)
    impute_results <- raiss(LD_data$ref_panel, sumstats, LD_matrix,
                            rcond = impute_opts$rcond, R2_threshold = impute_opts$R2_threshold,
                            minimum_ld = impute_opts$minimum_ld, lamb = impute_opts$lamb)
    sumstats <- impute_results$result_filter
    LD_mat <- impute_results$LD_mat
  }

  # Fine-mapping: use X_mat if available, otherwise R
  if (!is.null(finemapping_method)) {
    pri_coverage <- finemapping_opts$coverage[1]
    sec_coverage <- if (length(finemapping_opts$coverage) > 1) finemapping_opts$coverage[-1] else NULL

    # When using X path, subset X to QCed/imputed variants.
    # For mixture panels (list), subset each panel; susie_rss accepts X=list().
    if (use_X) {
      if (is_X_list) {
        X_mat_sub <- lapply(X_data, function(Xk) Xk[, sumstats$variant_id, drop = FALSE])
      } else {
        X_mat_sub <- X_data[, sumstats$variant_id, drop = FALSE]
      }
    } else {
      X_mat_sub <- NULL
    }

    res <- susie_rss_pipeline(sumstats,
      LD_mat = if (use_X) NULL else LD_mat,
      X_mat = X_mat_sub,
      n = n, var_y = var_y,
      L = finemapping_opts$init_L, max_L = finemapping_opts$max_L, l_step = finemapping_opts$l_step,
      analysis_method = finemapping_method,
      coverage = pri_coverage,
      secondary_coverage = sec_coverage,
      signal_cutoff = finemapping_opts$signal_cutoff,
      min_abs_corr = finemapping_opts$min_abs_corr,
      stochastic_ld_sample = stochastic_ld_sample
    )
    if (!is.null(qc_method)) {
      res$outlier_number <- qc_results$outlier_number
    }
  }
  if (impute & !is.null(qc_method)) {
    method_name <- paste0(finemapping_method, "_", toupper(qc_method), "_RAISS_imputed")
  } else if (!impute & !is.null(qc_method)) {
    method_name <- paste0(finemapping_method, "_", toupper(qc_method))
  } else {
    method_name <- paste0(finemapping_method, "_", "NO_QC")
  }
  result_list <- list()
  result_list[[method_name]] <- res
  result_list[["rss_data_analyzed"]] <- sumstats
  
  block_cs_metrics <- list()
  if (diagnostics) {
    if (length(res) > 0) { 
        bvsr_res = get_susie_result(res)
        bvsr_cs_num = if(!is.null(bvsr_res)) length(bvsr_res$sets$cs) else NULL
        if (isTRUE(bvsr_cs_num > 0)) { # have CS
            # Get the names of the credible sets
            cs_names_bvsr = names(bvsr_res$sets$cs)
            block_cs_metrics = extract_cs_info(con_data = res, cs_names = cs_names_bvsr, top_loci_table = res$top_loci)
        } else { # no CS
            if (sum(bvsr_res$pip > finemapping_opts$signal_cutoff) > 0) {
                block_cs_metrics = extract_top_pip_info(res)
            }
        }
    }
    # sensitive check for additional analyses
    if (!is.null(block_cs_metrics) && length(block_cs_metrics) > 0) {
      block_cs_metrics = parse_cs_corr(block_cs_metrics)
      cs_row = block_cs_metrics %>% filter(!is.na(block_cs_metrics$variants_per_cs))
      if (nrow(cs_row)>1) {# CS > 1
        block_cs_metrics <- block_cs_metrics %>%
          mutate(max_cs_corr_study_block = if(all(is.na(cs_corr_max))) {
            NA_real_
          } else {
            max(cs_corr_max, na.rm = TRUE)
          })
        if (any(block_cs_metrics$p_value > 1e-4 | block_cs_metrics$max_cs_corr_study_block > 0.5)) {
          finemapping_method <- "bayesian_conditional_regression"
          pri_coverage <- finemapping_opts$coverage[1]
          sec_coverage <- if (length(finemapping_opts$coverage) > 1) finemapping_opts$coverage[-1] else NULL
          bcr <- susie_rss_pipeline(sumstats, LD_mat,
            n = n, var_y = var_y,
            L = finemapping_opts$init_L, max_L = finemapping_opts$max_L, l_step = finemapping_opts$l_step,
            analysis_method = finemapping_method,
            coverage = pri_coverage,
            secondary_coverage = sec_coverage,
            signal_cutoff = finemapping_opts$signal_cutoff,
            min_abs_corr = finemapping_opts$min_abs_corr
          )
          if (!is.null(qc_method)) {
            bcr$outlier_number <- qc_results$outlier_number
          }
          if (impute & !is.null(qc_method)) {
              method_name <- paste0(finemapping_method, "_", toupper(qc_method), "_RAISS_imputed")
          } else if (!impute & !is.null(qc_method)) {
              method_name <- paste0(finemapping_method, "_", toupper(qc_method))
          } else {
              method_name <- paste0(finemapping_method, "_", "NO_QC")
          }
          result_list[[method_name]] <- bcr
          finemapping_method <- "single_effect"
          sumstats <- preprocess_results$sumstats
          LD_mat <- preprocess_results$LD_mat
          ser <- susie_rss_pipeline(sumstats, LD_mat,
            n = n, var_y = var_y,
            L = finemapping_opts$init_L, max_L = finemapping_opts$max_L, l_step = finemapping_opts$l_step,
            analysis_method = finemapping_method,
            coverage = pri_coverage,
            secondary_coverage = sec_coverage,
            signal_cutoff = finemapping_opts$signal_cutoff,
            min_abs_corr = finemapping_opts$min_abs_corr
          )
          result_list[[paste0(finemapping_method, "_NO_QC")]] <- ser
        }
      } else { # CS = 1 or NA
        finemapping_method <- "single_effect"
        sumstats <- preprocess_results$sumstats
        LD_mat <- preprocess_results$LD_mat
        ser <- susie_rss_pipeline(sumstats, LD_mat,
          n = n, var_y = var_y,
          L = finemapping_opts$init_L, max_L = finemapping_opts$max_L, l_step = finemapping_opts$l_step,
          analysis_method = finemapping_method,
          coverage = pri_coverage,
          secondary_coverage = sec_coverage,
          signal_cutoff = finemapping_opts$signal_cutoff,
          min_abs_corr = finemapping_opts$min_abs_corr
        )
        result_list[[paste0(finemapping_method, "_NO_QC")]] <- ser
      }
    result_list[["diagnostics"]] <- block_cs_metrics
    }
  }
  return(result_list)
}
