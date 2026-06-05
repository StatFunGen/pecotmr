#' Univariate Analysis Pipeline
#'
#' This function performs univariate analysis for fine-mapping and Transcriptome-Wide Association Study (TWAS)
#' with optional cross-validation. By default, fine-mapping fits SuSiE-inf first
#' and then fits SuSiE initialized from the SuSiE-inf result.
#'
#' @param X A matrix of genotype data where rows represent samples and columns represent genetic variants.
#' @param Y A vector of phenotype measurements.
#' @param X_scalar A scalar or vector to rescale X to its original scale.
#' @param Y_scalar A scalar to rescale Y to its original scale.
#' @param maf A vector of minor allele frequencies for each variant in X.
#' @param X_variance Optional variance of X. Default is NULL.
#' @param other_quantities A list of other quantities to be carried into fine-mapping post-processing. Default is an empty list.
#' @param region Optional \code{"chr:start-end"} string for the analysis region. Default is NULL.
#' @param imiss_cutoff Individual missingness cutoff. Default is 1.0.
#' @param maf_cutoff Minor allele frequency cutoff. Default is NULL.
#' @param xvar_cutoff Variance cutoff for X. Default is 0.05.
#' @param ld_reference_meta_file An optional path to a file containing linkage disequilibrium reference data. Default is NULL.
#' @param pip_cutoff_to_skip Cutoff value for skipping analysis based on PIP values. Default is 0.
#' @param L Maximum number of components in SuSiE. Default is 20.
#' @param L_greedy Initial greedy number of components in SuSiE. Default is 5.
#' @param signal_cutoff Cutoff value for signal identification in PIP values. Default is 0.025.
#' @param coverage A vector of coverage probabilities for credible sets. Default is c(0.95, 0.7, 0.5).
#' @param min_abs_corr Minimum absolute correlation for credible set purity filtering. Default is 0.8,
#'   which is stricter than the susieR default of 0.5.
#' @param finemapping_extra_opts Additional options passed to \code{susieR::susie()}.
#'   SuSiE-inf is always fitted with \code{refine = FALSE}; the ordinary SuSiE
#'   fit keeps these options and is initialized with \code{model_init}.
#' @param estimate_residual_variance Passed to \code{susieR::susie()}. Default is TRUE.
#' @param methods Optional character vector selecting which SuSiE variants to
#'   fit. Any subset of \code{c("susie", "susie_inf", "susie_ash")}. Default
#'   \code{NULL} falls back to the legacy \code{add_susie_inf} behavior:
#'   \code{add_susie_inf = TRUE} (default) maps to
#'   \code{methods = c("susie_inf", "susie")} with SuSiE-inf chained into the
#'   SuSiE fit as initialization; \code{add_susie_inf = FALSE} maps to
#'   \code{methods = "susie"} (plain SuSiE alone). When \code{methods} is
#'   passed explicitly, each requested method is fitted; if
#'   \code{"susie_inf"} is paired with \code{"susie"} or \code{"susie_ash"}
#'   (or both) and \code{add_susie_inf = TRUE}, the SuSiE-inf fit
#'   initialises each chained downstream method. This gives five distinct
#'   fitting modes: SuSiE alone, SuSiE with SuSiE-inf init, SuSiE-inf alone,
#'   SuSiE-ash alone, and SuSiE-ash with SuSiE-inf init.
#' @param add_susie_inf When \code{methods} is \code{NULL}, controls whether
#'   SuSiE-inf is fitted and chained into SuSiE. When \code{methods} is set
#'   explicitly, controls whether the chained-init shortcut is applied to
#'   any \code{"susie"} or \code{"susie_ash"} method present alongside
#'   \code{"susie_inf"}. Default \code{TRUE}.
#' @param twas_weights Whether to compute TWAS weights. Default is TRUE.
#' @param sample_partition Optional data frame with Sample and Fold columns for cross-validation. Default is NULL.
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
    region = NULL,
    # filters
    imiss_cutoff = 1.0,
    maf_cutoff = NULL,
    xvar_cutoff = 0,
    ld_reference_meta_file = NULL,
    pip_cutoff_to_skip = 0,
    # methods parameter configuration
    L = 20,
    L_greedy = 5,
    # fine-mapping results summary
    signal_cutoff = 0.025,
    coverage = c(0.95, 0.7, 0.5),
    min_abs_corr = 0.8,
    finemapping_extra_opts = list(refine = TRUE),
    estimate_residual_variance = TRUE,
    methods = NULL,
    add_susie_inf = TRUE,
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
  if (!is.numeric(L) || L <= 0) stop("L must be a positive integer")
  if (!is.null(L_greedy) && (!is.numeric(L_greedy) || L_greedy <= 0)) stop("L_greedy must be NULL or a positive integer")
  if (!is.logical(add_susie_inf) || length(add_susie_inf) != 1 || is.na(add_susie_inf)) {
    stop("add_susie_inf must be TRUE or FALSE")
  }

  # Resolve effective methods. NULL => backward-compat via add_susie_inf.
  valid_methods <- c("susie", "susie_inf", "susie_ash")
  if (is.null(methods)) {
    methods <- if (isTRUE(add_susie_inf)) c("susie_inf", "susie") else "susie"
  } else {
    if (!is.character(methods) || length(methods) == 0L) {
      stop("methods must be a non-empty character vector of method names.")
    }
    bad <- setdiff(methods, valid_methods)
    if (length(bad) > 0) {
      stop("Unknown method(s): ", paste(bad, collapse = ", "),
           ". Valid options: ", paste(valid_methods, collapse = ", "))
    }
    methods <- unique(methods)
  }
  # SuSiE-inf initialisation chains into SuSiE and/or SuSiE-ash whenever
  # either of them is requested alongside SuSiE-inf and add_susie_inf is TRUE.
  chain_inf_to_susie     <- isTRUE(add_susie_inf) &&
    all(c("susie_inf", "susie") %in% methods)
  chain_inf_to_susie_ash <- isTRUE(add_susie_inf) &&
    all(c("susie_inf", "susie_ash") %in% methods)
  any_chained_init <- chain_inf_to_susie || chain_inf_to_susie_ash
  if (isTRUE(twas_weights) && !("susie" %in% methods)) {
    stop("twas_weights = TRUE requires \"susie\" to be in methods.")
  }
  if (isTRUE(twas_weights) && !chain_inf_to_susie) {
    stop("twas_weights = TRUE requires SuSiE to be initialised from SuSiE-inf; ",
         "set methods = c(\"susie_inf\", \"susie\") and add_susie_inf = TRUE.")
  }

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

  susie_args <- modifyList(
    finemapping_extra_opts,
    list(L = L, L_greedy = L_greedy, coverage = coverage[1],
         estimate_residual_variance = estimate_residual_variance)
  )
  fitted_models <- list()

  if ("susie_inf" %in% methods || any_chained_init) {
    message("Fitting SuSiE-inf model on input data ...")
    inf_args <- modifyList(susie_args, list(
      X = X, y = Y,
      unmappable_effects = "inf",
      convergence_method = "pip",
      refine = FALSE, model_init = NULL
    ))
    inf_fit <- do.call(susie, inf_args)
    fitted_models[["susie_inf"]] <- .set_finemapping_fit_class(inf_fit, "susie_inf")
  }

  if ("susie" %in% methods) {
    if (chain_inf_to_susie) {
      message("Fitting SuSiE model initialized by SuSiE-inf ...")
      su_args <- prepare_susie_from_inf_args(susie_args,
                                             fitted_models[["susie_inf"]],
                                             refine_default = TRUE,
                                             unmappable_effects = "none")
      su_fit <- do.call(susie, c(list(X = X, y = Y), su_args))
    } else {
      message("Fitting SuSiE model on input data ...")
      su_fit <- do.call(susie, c(list(X = X, y = Y), susie_args))
    }
    fitted_models[["susie"]] <- .set_finemapping_fit_class(su_fit, "susie")
  }

  if ("susie_ash" %in% methods) {
    if (chain_inf_to_susie_ash) {
      message("Fitting SuSiE-ash model initialized by SuSiE-inf ...")
      ash_args <- prepare_susie_from_inf_args(susie_args,
                                              fitted_models[["susie_inf"]],
                                              refine_default = NULL,
                                              unmappable_effects = "ash")
      ash_fit <- do.call(susie, c(list(X = X, y = Y), ash_args))
    } else {
      message("Fitting SuSiE-ash model on input data ...")
      ash_args <- modifyList(susie_args, list(
        X = X, y = Y,
        unmappable_effects = "ash",
        convergence_method = "pip"
      ))
      ash_fit <- do.call(susie, ash_args)
    }
    fitted_models[["susie_ash"]] <- .set_finemapping_fit_class(ash_fit, "susie_ash")
  }

  # Drop susie_inf from post-processing if it was only fit to provide init for
  # SuSiE / SuSiE-ash (i.e. caller did not request "susie_inf" in methods).
  if (any_chained_init && !("susie_inf" %in% methods)) {
    fitted_models[["susie_inf"]] <- NULL
  }

  # Back-compat slots for the most common methods
  res$susie_inf_fitted <- fitted_models[["susie_inf"]]
  res$susie_fitted     <- fitted_models[["susie"]]
  res$susie_ash_fitted <- fitted_models[["susie_ash"]]

  # Process SuSiE results
  susie_post <- postprocess_finemapping_fits(
    fits = fitted_models,
    data_x = X,
    data_y = Y,
    X_scalar = X_scalar,
    y_scalar = Y_scalar,
    maf = maf,
    coverage = coverage[1],
    secondary_coverage = if (length(coverage) > 1) coverage[-1] else NULL,
    signal_cutoff = signal_cutoff,
    min_abs_corr = min_abs_corr,
    other_quantities = other_quantities,
    region = region
  )
  # Primary method drives root-level finemapping_result / sumstats / etc.
  # Preference order favors "susie" for backward compatibility, then
  # falls back to the first requested method actually fitted.
  primary_method <- if ("susie" %in% names(fitted_models)) "susie" else names(fitted_models)[1]
  res <- c(res, format_finemapping_output(susie_post, primary_method = primary_method))
  susie_inf_fm <- susie_post$finemapping_results$susie_inf$finemapping_result
  res$susie_inf_result_trimmed <- if (!is.null(susie_inf_fm)) getTrimmedFit(susie_inf_fm) else NULL
  susie_ash_fm <- susie_post$finemapping_results$susie_ash$finemapping_result
  res$susie_ash_result_trimmed <- if (!is.null(susie_ash_fm)) getTrimmedFit(susie_ash_fm) else NULL
  res$total_time_elapsed <- proc.time() - st

  # TWAS weights and cross-validation
  if (twas_weights) {
    res$twas_weights_result <- twas_weights_pipeline(
      X, Y, fitted_models = fitted_models,
      cv_folds = cv_folds,
      max_cv_variants = max_cv_variants,
      cv_threads = cv_threads,
      sample_partition = sample_partition
    )
    if ("top_loci" %in% names(res) && !is.null(res$twas_weights_result$susie_weights_intermediate)) {
      res$twas_weights_result$susie_weights_intermediate$top_loci <- res$top_loci
    }
  }

  return(res)
}

#' Load LD for a study, supporting single or mixture panels.
#'
#' @param ld_path A single LD metadata TSV path, or comma-separated paths for
#'   mixture panels (e.g., "ld_EUR.tsv,ld_AFR.tsv").
#' @param region Region string "chr:start-end".
#' @return An \code{LDData} S4 object. For single panels, returns the result of
#'   \code{load_LD_matrix()} unchanged. For mixture panels, \code{genotype_handle}
#'   is a list of per-panel genotype handles sharing the first panel's variants.
#' @export
load_study_LD <- function(ld_path, region) {
  paths <- strsplit(ld_path, ",")[[1]]
  if (length(paths) == 1) {
    return(load_LD_matrix(paths, region, return_genotype = "auto"))
  }
  # Mixture: load each panel; combine handles into a list
  base <- load_LD_matrix(paths[1], region, return_genotype = TRUE)
  other_handles <- lapply(paths[-1], function(p) {
    ld <- load_LD_matrix(p, region, return_genotype = TRUE)
    ld@genotype_handle
  })
  all_handles <- c(list(base@genotype_handle), other_handles)
  LDData(
    correlation = NULL,
    genotype_handle = all_handles,
    snp_idx = base@snp_idx,
    variants = base@variants,
    block_metadata = base@block_metadata,
    n_ref = base@n_ref
  )
}

#' RSS Analysis Pipeline
#'
#' End-to-end pipeline for summary statistics fine-mapping via SuSiE RSS.
#' Supports both z+R (correlation matrix) and z+X (genotype matrix) interfaces.
#'
#' @param sumstat_path File path to the summary statistics.
#' @param column_file_path File path to the column mapping file.
#' @param LD_data An \code{LDData} S4 object from \code{load_LD_matrix()}. When
#'   \code{hasGenotypes(LD_data)} is TRUE (from \code{return_genotype=TRUE}),
#'   susie_rss uses the z+X interface via \code{getGenotypes()}. Local R is
#'   computed only for QC stages that require a correlation matrix.
#' @param n_sample Sample size. If 0, retrieved from the sumstat file.
#' @param n_case Number of cases (for case-control studies).
#' @param n_control Number of controls (for case-control studies).
#' @param region Region string "chr:start-end" for tabix subsetting.
#' @param skip_region Character vector of regions to skip (format "chrom:start-end").
#' @param extract_region_name Gene/phenotype name to subset.
#' @param region_name_col Column to filter for extract_region_name.
#' @param qc_method Summary-statistic QC method. \code{"slalom"} and
#'   \code{"dentist"} run basic allele harmonization plus LD-mismatch QC;
#'   \code{"none"} runs basic allele harmonization only.
#' @param finemapping_method Iteration mode for the SuSiE-RSS fit (when
#'   \code{"susie_rss"} is among \code{methods}). One of \code{"susie_rss"}
#'   (default normal IBSS), \code{"single_effect"} (L=1, single iteration),
#'   or \code{"bayesian_conditional_regression"} (full L, single iteration).
#' @param methods Optional character vector selecting which SuSiE-RSS
#'   variants to fit. Any subset of \code{c("susie_rss", "susie_inf_rss",
#'   "susie_ash_rss")}. Default \code{NULL} preserves legacy single-method
#'   behavior via \code{finemapping_method}. When set explicitly, every
#'   requested method contributes rows to the unified \code{top_loci}; when
#'   \code{"susie_inf_rss"} is paired with \code{"susie_rss"} or
#'   \code{"susie_ash_rss"} (or both) and \code{add_susie_inf = TRUE}, the
#'   SuSiE-inf-RSS fit initialises the chained downstream method(s).
#' @param add_susie_inf Logical controlling chained init when
#'   \code{"susie_inf_rss"} is in \code{methods} alongside
#'   \code{"susie_rss"} and/or \code{"susie_ash_rss"}. Default \code{TRUE}.
#' @param finemapping_opts List of fine-mapping options (L, L_greedy, coverage,
#'   signal_cutoff, min_abs_corr).
#' @param impute Whether to impute missing variants via RAISS (default TRUE).
#' @param impute_opts List of imputation options (rcond, R2_threshold, minimum_ld, lamb).
#' @param pip_cutoff_to_skip PIP threshold for early stopping (default 0, no skip).
#' @param R_finite Controls variance inflation to account for finite reference LD.
#'   Passed to \code{susieR::susie_rss()}.
#' @param R_mismatch LD mismatch correction method passed to \code{susieR::susie_rss()}.
#'   Default NULL disables mismatch correction.
#' @param keep_indel Whether to keep indel variants (default TRUE).
#' @param comment_string Comment character for sumstat file (default "#").
#' @param diagnostics Whether to include diagnostic info (default FALSE).
#'
#' @return A list with fine-mapping results and analyzed summary statistics.
#' @importFrom magrittr %>%
#' @importFrom susieR susie_rss
#' @export
rss_analysis_pipeline <- function(
    sumstat_path, column_file_path, LD_data,
    n_sample = 0, n_case = 0, n_control = 0, region = NULL, skip_region = NULL,
    extract_region_name = NULL, region_name_col = NULL,
    qc_method = c("slalom", "dentist", "none"),
    finemapping_method = c("susie_rss", "single_effect", "bayesian_conditional_regression"),
    methods = NULL,
    add_susie_inf = TRUE,
    finemapping_opts = list(
      L = 20, L_greedy = 5,
      coverage = c(0.95, 0.7, 0.5), signal_cutoff = 0.025,
      min_abs_corr = 0.8
    ),
    impute = TRUE, impute_opts = list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5, lamb = 0.01),
    pip_cutoff_to_skip = 0, R_finite = NULL, R_mismatch = NULL,
    keep_indel = TRUE, comment_string = "#", diagnostics = FALSE) {
  if (!is(LD_data, "LDData")) {
    stop("LD_data must be an LDData object")
  }
  use_X <- hasGenotypes(LD_data)
  X_data <- if (use_X) getGenotypes(LD_data) else NULL
  is_X_list <- use_X && is.list(X_data)
  subset_X_data <- function(variants) {
    if (!use_X) return(NULL)
    if (is_X_list) {
      lapply(X_data, function(Xk) Xk[, variants, drop = FALSE])
    } else {
      X_data[, variants, drop = FALSE]
    }
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

  qc_method_arg <- if (is.null(qc_method)) NULL else match.arg(qc_method)
  qc_method <- qc_method_arg
  qc_record <- summary_stats_qc(
    rss_input = rss_input,
    LD_data = LD_data,
    keep_indel = keep_indel,
    skip_region = skip_region,
    pip_cutoff_to_skip = pip_cutoff_to_skip,
    qc_method = if (is.null(qc_method_arg)) "none" else qc_method_arg,
    impute = impute,
    impute_opts = impute_opts,
    return_on_skip = "preprocess",
    R_finite = R_finite,
    R_mismatch = R_mismatch
  )
  if (!is(qc_record, "QCResult")) {
    stop("summary_stats_qc must return a QCResult object.")
  }
  rss_record <- getRSSInput(qc_record)
  sumstats <- rss_record$sumstats
  qc_ld <- getLDData(qc_record)
  LD_mat <- if (is.null(qc_ld)) NULL else if (hasGenotypes(qc_ld)) getGenotypes(qc_ld) else getCorrelation(qc_ld)
  preprocess_snapshot <- getPreprocess(qc_record)
  preprocess_ld <- preprocess_snapshot$ld_data
  preprocess_results <- list(
    sumstats = preprocess_snapshot$sumstats,
    LD_mat = if (is.null(preprocess_ld)) NULL else if (hasGenotypes(preprocess_ld)) getGenotypes(preprocess_ld) else getCorrelation(preprocess_ld)
  )
  qc_results <- list(outlier_number = getOutlierNumber(qc_record))

  if (nrow(sumstats) == 0) {
    message("No variants left after preprocessing. Returning empty results.")
    return(list(rss_data_analyzed = sumstats))
  }
  if (isSkipped(qc_record)) {
    return(list(rss_data_analyzed = sumstats))
  }

  # Fine-mapping: use X_mat if available, otherwise R
  if (!is.null(finemapping_method)) {
    pri_coverage <- finemapping_opts$coverage[1]
    sec_coverage <- if (length(finemapping_opts$coverage) > 1) finemapping_opts$coverage[-1] else NULL

    X_mat_sub <- subset_X_data(sumstats$variant_id)

    res <- susie_rss_pipeline(sumstats,
      LD_mat = if (use_X) NULL else LD_mat,
      X_mat = X_mat_sub,
      n = n, var_y = var_y,
      L = finemapping_opts$L, L_greedy = finemapping_opts$L_greedy,
      analysis_method = finemapping_method,
      methods = methods,
      add_susie_inf = add_susie_inf,
      coverage = pri_coverage,
      secondary_coverage = sec_coverage,
      signal_cutoff = finemapping_opts$signal_cutoff,
      min_abs_corr = finemapping_opts$min_abs_corr,
      R_finite = R_finite,
      R_mismatch = R_mismatch
    )
    if (!is.null(qc_method)) {
      res$outlier_number <- qc_results$outlier_number
    }
  }
  .make_method_name <- function(method, qc_method, impute) {
    suffix <- if (!is.null(qc_method) && impute) {
      paste0(toupper(qc_method), "_RAISS_imputed")
    } else if (!is.null(qc_method)) {
      toupper(qc_method)
    } else {
      "NO_QC"
    }
    paste0(method, "_", suffix)
  }

  .run_reanalysis <- function(sumstats, LD_mat, method, finemapping_opts, pri_coverage, sec_coverage) {
    susie_rss_pipeline(sumstats,
      LD_mat = if (use_X) NULL else LD_mat,
      X_mat = subset_X_data(sumstats$variant_id),
      n = n, var_y = var_y,
      L = finemapping_opts$L, L_greedy = finemapping_opts$L_greedy,
      analysis_method = method,
      coverage = pri_coverage,
      secondary_coverage = sec_coverage,
      signal_cutoff = finemapping_opts$signal_cutoff,
      min_abs_corr = finemapping_opts$min_abs_corr,
      R_finite = R_finite,
      R_mismatch = R_mismatch
    )
  }

  method_name <- .make_method_name(finemapping_method, qc_method, impute)
  result_list <- list()
  result_list[[method_name]] <- res
  result_list[["rss_data_analyzed"]] <- sumstats

  block_cs_metrics <- list()
  if (diagnostics) {
    if (length(res) > 0) {
        bvsr_res = get_susie_result(res)
        bvsr_cs_num = if(!is.null(bvsr_res)) length(bvsr_res$sets$cs) else NULL
        if (isTRUE(bvsr_cs_num > 0)) { # have CS
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
          bcr <- .run_reanalysis(sumstats, LD_mat, "bayesian_conditional_regression",
            finemapping_opts, pri_coverage, sec_coverage)
          if (!is.null(qc_method)) {
            bcr$outlier_number <- qc_results$outlier_number
          }
          result_list[[.make_method_name("bayesian_conditional_regression", qc_method, impute)]] <- bcr
          ser <- .run_reanalysis(preprocess_results$sumstats, preprocess_results$LD_mat,
            "single_effect", finemapping_opts, pri_coverage, sec_coverage)
          result_list[["single_effect_NO_QC"]] <- ser
        }
      } else { # CS = 1 or NA
        ser <- .run_reanalysis(preprocess_results$sumstats, preprocess_results$LD_mat,
          "single_effect", finemapping_opts, pri_coverage, sec_coverage)
        result_list[["single_effect_NO_QC"]] <- ser
      }
    result_list[["diagnostics"]] <- block_cs_metrics
    }
  }
  return(result_list)
}
