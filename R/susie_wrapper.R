#' Convert Log Bayes Factors to Single Effects PIP
#'
#' This function converts log Bayes factors (LBF) to alpha values, optionally
#' using prior weights. It handles numerical stability by adjusting with the
#' maximum LBF value.
#'
#' @param lbf Numeric vector of log Bayes factors.
#' @param prior_weights Optional numeric vector of prior weights for each element in lbf.
#' @return A named numeric vector of alpha values corresponding to the input LBF.
#' @examples
#' lbf <- c(-0.5, 1.2, 0.3)
#' alpha <- lbf_to_alpha_vector(lbf)
#' print(alpha)
#' @noRd
lbf_to_alpha_vector <- function(lbf, prior_weights = NULL) {
  if (is.null(prior_weights)) prior_weights <- rep(1 / length(lbf), length(lbf))
  maxlbf <- max(lbf)

  # If maxlbf is 0, return a vector of zeros
  if (maxlbf == 0) {
    return(setNames(rep(0, length(lbf)), names(lbf)))
  }

  # w is proportional to BF, subtract max for numerical stability
  w <- exp(lbf - maxlbf)

  # Posterior prob for each SNP
  w_weighted <- w * prior_weights
  weighted_sum_w <- sum(w_weighted)
  alpha <- w_weighted / weighted_sum_w

  return(alpha)
}

#' Applies the 'lbf_to_alpha_vector' function row-wise to a matrix of log Bayes factors
#' to convert them to Single Effect PIP values.
#'
#' @param lbf Matrix of log Bayes factors.
#' @return A matrix of alpha values with the same dimensions as the input LBF matrix.
#' @examples
#' lbf_matrix <- matrix(c(-0.5, 1.2, 0.3, 0.7, -1.1, 0.4), nrow = 2)
#' alpha_matrix <- lbf_to_alpha(lbf_matrix)
#' print(alpha_matrix)
#' @export
lbf_to_alpha <- function(lbf) {
  alpha_matrix <- t(apply(as.matrix(lbf), 1, lbf_to_alpha_vector))
  if (ncol(lbf) == 1) alpha_matrix <- matrix(alpha_matrix, ncol = 1, dimnames = list(NULL, colnames(lbf)))
  return(alpha_matrix)
}

#' Adjust SuSiE Weights
#'
#' Adjusts SuSiE TWAS weights by subsetting to intersected variants and
#' optionally running allele QC against LD reference variants.
#'
#' @param twas_weights_results A list containing TWAS weight data (nested structure).
#' @param keep_variants Vector of variant names to keep.
#' @param run_allele_qc Whether to run allele_qc to align alleles. Default TRUE.
#' @param variable_name_obj Path to variant names in the nested list.
#' @param susie_obj Path to susie result in the nested list.
#' @param twas_weights_table Path to weights table in the nested list.
#' @param LD_variants Vector of LD reference variant IDs for allele QC.
#' @param match_min_prop Minimum proportion of matched variants. Default 0.2.
#' @return A list with adjusted_susie_weights and remained_variants_ids.
#' @export
adjust_susie_weights <- function(twas_weights_results, keep_variants, run_allele_qc = TRUE,
                                 variable_name_obj = c("susie_results", context, "variant_names"),
                                 susie_obj = c("susie_results", context, "susie_result_trimmed"),
                                 twas_weights_table = c("weights", context), LD_variants, match_min_prop = 0.2) {
  # Intersect the rownames of weights with keep_variants
  twas_weights_variants <- get_nested_element(twas_weights_results, variable_name_obj)
  # Normalize to canonical format (with chr prefix)
  twas_weights_variants <- normalize_variant_id(twas_weights_variants)
  # allele flip twas weights matrix variants name
  if (run_allele_qc) {
    weights_matrix <- get_nested_element(twas_weights_results, twas_weights_table)
    if (!all(c("chrom", "pos", "A2", "A1") %in% colnames(weights_matrix))) {
      weights_matrix <- cbind(parse_variant_id(twas_weights_variants), weights_matrix)
    }
    weights_matrix_qced <- allele_qc(weights_matrix, LD_variants, colnames(weights_matrix)[!colnames(weights_matrix) %in% c(
      "chrom",
      "pos", "A2", "A1"
    )], match_min_prop = match_min_prop)
    # allele_qc now outputs canonical variant_ids (with chr prefix)
    original_idx <- match(weights_matrix_qced$qc_summary$variants_id_original, twas_weights_variants)
    intersected_indices <- original_idx[weights_matrix_qced$qc_summary$keep == TRUE]
  } else {
    # Normalize keep_variants to canonical format for matching
    keep_variants_normalized <- normalize_variant_id(keep_variants)
    intersected_variants <- intersect(twas_weights_variants, keep_variants_normalized)
    intersected_indices <- match(intersected_variants, twas_weights_variants)
  }
  if (length(intersected_indices) == 0) {
    stop("Error: No intersected variants found. Please check 'twas_weights' and 'keep_variants' inputs to make sure there are variants left to use.")
  }
  # Subset lbf_matrix, mu, and x_column_scale_factors
  lbf_matrix <- get_nested_element(twas_weights_results, c(susie_obj, "lbf_variable"))
  mu <- get_nested_element(twas_weights_results, c(susie_obj, "mu"))
  x_column_scal_factors <- get_nested_element(twas_weights_results, c(susie_obj, "X_column_scale_factors"))

  lbf_matrix_subset <- lbf_matrix[, intersected_indices, drop = FALSE]
  mu_subset <- mu[, intersected_indices, drop = FALSE]
  x_column_scal_factors_subset <- x_column_scal_factors[intersected_indices]

  # Convert lbf_matrix to alpha and calculate adjusted xQTL coefficients
  adjusted_xqtl_alpha <- lbf_to_alpha(lbf_matrix_subset)
  adjusted_xqtl_coef <- colSums(adjusted_xqtl_alpha * mu_subset) / x_column_scal_factors_subset
  # allele_qc now outputs canonical variant_ids (with chr prefix) -- no need to add chr
  remained_variants_ids <- if (run_allele_qc) {
    weights_matrix_qced$target_data_qced$variant_id
  } else {
    intersected_variants
  }
  return(list(adjusted_susie_weights = adjusted_xqtl_coef, remained_variants_ids = remained_variants_ids))
}


#' @importFrom susieR susie
#' @export
susie_wrapper <- function(X, y, init_L = 5, max_L = 30, l_step = 5, ...) {
  if (init_L == max_L) {
    return(susie(X, y, L = init_L, ...))
  }
  L <- init_L
  # Perform SuSiE by dynamically increasing L
  gst <- proc.time()
  while (TRUE) {
    st <- proc.time()
    res <- susie(X, y, L = L, ...)
    res$time_elapsed <- proc.time() - st
    if (!is.null(res$sets$cs)) {
      if (length(res$sets$cs) >= L && L <= max_L) {
        L <- L + l_step
      } else {
        break
      }
    } else {
      break
    }
  }
  message(paste("Total time elapsed for susie_wrapper:", (proc.time() - gst)[3]))
  return(res)
}

#' Wrapper Function for SuSiE RSS with Dynamic L Adjustment
#'
#' Performs SuSiE RSS analysis with dynamic L adjustment. Supports both
#' z+R (correlation matrix) and z+X (genotype matrix) interfaces.
#'
#' @param z Z score vector.
#' @param R LD correlation matrix. Mutually exclusive with X.
#' @param X Genotype matrix (samples x variants). Mutually exclusive with R.
#'   When provided, susie_rss uses the low-rank X interface.
#' @param n Sample size.
#' @param L Initial number of causal configurations.
#' @param max_L Maximum number of causal configurations.
#' @param l_step Step size for increasing L when the limit is reached.
#' @param stochastic_ld_sample Stochastic LD parameter passed to susie_rss.
#'   NULL (default): no variance inflation. TRUE: infer sketch size from X
#'   (requires X, not R). Integer: explicit sketch size (for R only).
#' @param ... Extra parameters passed to susie_rss (e.g., var_y, coverage).
#' @return SuSiE RSS fit object after dynamic L adjustment
#' @importFrom susieR susie_rss
#' @export
susie_rss_wrapper <- function(z, R = NULL, X = NULL, n = NULL,
                              L = 10, max_L = 30, l_step = 5,
                              coverage = 0.95,
                              stochastic_ld_sample = NULL, ...) {
  # Validate: exactly one of R or X
  if (is.null(R) && is.null(X)) stop("Either R or X must be provided.")
  if (!is.null(R) && !is.null(X)) stop("Only one of R or X should be provided, not both.")

  # Build argument list for susie_rss
  base_args <- list(z = z, n = n, L = L, coverage = coverage,
                    stochastic_ld_sample = stochastic_ld_sample, ...)
  if (!is.null(X)) base_args$X <- X else base_args$R <- R

  run_susie <- function(args) do.call(susie_rss, args)

  if (L == 1) {
    base_args$max_iter <- 1
    result <- run_susie(base_args)
  } else if (L == max_L) {
    result <- run_susie(base_args)
  } else {
    while (TRUE) {
      st <- proc.time()
      base_args$L <- L
      result <- run_susie(base_args)
      result$time_elapsed <- proc.time() - st
      if (!is.null(result$sets$cs) && length(result$sets$cs) >= L && L <= max_L) {
        L <- L + l_step
      } else {
        break
      }
    }
  }

  result
}

#' Run the SuSiE RSS pipeline
#'
#' Runs SuSiE RSS analysis with the specified method. Supports both z+R
#' (correlation matrix) and z+X (genotype matrix) interfaces.
#'
#' @param sumstats Data frame with 'z' or ('beta' and 'se') columns.
#' @param LD_mat LD correlation matrix. Mutually exclusive with X_mat.
#' @param X_mat Genotype matrix (samples x variants). Mutually exclusive with LD_mat.
#' @param n Sample size.
#' @param L Initial number of causal configurations (default: 5).
#' @param max_L Maximum number of causal configurations (default: 30).
#' @param l_step Step size for dynamic L adjustment (default: 5).
#' @param analysis_method One of "susie_rss", "single_effect", "bayesian_conditional_regression".
#' @param coverage Coverage level (default: 0.95).
#' @param secondary_coverage Secondary coverage levels (default: c(0.7, 0.5)).
#' @param signal_cutoff PIP cutoff for susie_post_processor (default: 0.1).
#' @param min_abs_corr Minimum absolute correlation for CS purity (default: 0.8).
#' @param stochastic_ld_sample Passed to susie_rss. NULL, TRUE, or integer.
#' @param ... Additional parameters passed to susie_rss (e.g., var_y).
#' @return A list with post-processed SuSiE RSS results.
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange select
#' @export
susie_rss_pipeline <- function(sumstats, LD_mat = NULL, X_mat = NULL, n = NULL,
                               L = 5, max_L = 30, l_step = 5,
                               analysis_method = c("susie_rss", "single_effect", "bayesian_conditional_regression"),
                               coverage = 0.95,
                               secondary_coverage = c(0.7, 0.5),
                               signal_cutoff = 0.1,
                               min_abs_corr = 0.8,
                               stochastic_ld_sample = NULL, ...) {
  analysis_method <- match.arg(analysis_method)

  if (!is.null(sumstats$z)) {
    z <- sumstats$z
  } else if (!is.null(sumstats$beta) && !is.null(sumstats$se)) {
    z <- sumstats$beta / sumstats$se
  } else {
    stop("sumstats must have 'z' or ('beta' and 'se') columns.")
  }

  # Common args for susie_rss_wrapper
  common <- list(z = z, n = n, coverage = coverage,
                 stochastic_ld_sample = stochastic_ld_sample, ...)
  if (!is.null(X_mat)) common$X <- X_mat else common$R <- LD_mat

  if (analysis_method == "single_effect") {
    res <- do.call(susie_rss_wrapper, c(common, list(L = 1)))
  } else if (analysis_method == "bayesian_conditional_regression") {
    res <- do.call(susie_rss_wrapper, c(common, list(L = L, max_L = max_L, l_step = l_step, max_iter = 1)))
  } else {
    res <- do.call(susie_rss_wrapper, c(common, list(L = L, max_L = max_L, l_step = l_step)))
  }

  # For post-processing, need a square matrix (R or computed from X).
  # For mixture panels (list of X), use the first panel to compute R.
  if (!is.null(LD_mat)) {
    data_x <- LD_mat
  } else if (is.list(X_mat) && !is.matrix(X_mat)) {
    data_x <- compute_LD(X_mat[[1]][, seq_along(z), drop = FALSE], method = "sample")
  } else {
    data_x <- compute_LD(X_mat[, seq_along(z), drop = FALSE], method = "sample")
  }

  res <- susie_post_processor(res,
    data_x = data_x, data_y = list(z = z),
    signal_cutoff = signal_cutoff, secondary_coverage = secondary_coverage,
    min_abs_corr = min_abs_corr, mode = "susie_rss"
  )
  res
}

#' @noRd
get_cs_index <- function(snps_idx, susie_cs) {
  # Return ALL CS indices that contain this variant (not just one)
  idx <- which(vapply(susie_cs, function(x) snps_idx %in% x, logical(1)))
  if (length(idx) == 0) return(NA_integer_)
  return(idx)
}
#' @noRd
get_top_variants_idx <- function(susie_output, signal_cutoff) {
  c(which(susie_output$pip >= signal_cutoff), unlist(susie_output$sets$cs)) %>%
    unique() %>%
    sort()
}
# Returns a data.frame(variant_idx, cs_idx) with one row per (variant, CS) pair.
# Variants in multiple CSs get multiple rows.
#' @noRd
get_cs_info <- function(susie_output_sets_cs, top_variants_idx) {
  cs_names <- names(susie_output_sets_cs)
  rows <- lapply(top_variants_idx, function(vi) {
    idx <- get_cs_index(vi, susie_output_sets_cs)
    if (length(idx) == 1 && is.na(idx)) {
      data.frame(variant_idx = vi, cs_idx = 0L, stringsAsFactors = FALSE)
    } else {
      cs_nums <- as.integer(str_replace(cs_names[idx], "L", ""))
      data.frame(variant_idx = rep(vi, length(cs_nums)), cs_idx = cs_nums, stringsAsFactors = FALSE)
    }
  })
  do.call(rbind, rows)
}
#' @noRd
get_cs_and_corr <- function(susie_output, coverage, data_x, mode = c("susie", "susie_rss", "mvsusie"), min_abs_corr = NULL) {
  if (mode %in% c("susie", "mvsusie")) {
    susie_output_secondary <- list(sets = susie_get_cs(susie_output, X = data_x, coverage = coverage, min_abs_corr = min_abs_corr), pip = susie_output$pip)
    susie_output_secondary$cs_corr <- get_cs_correlation(susie_output_secondary, X = data_x)
  } else {
    susie_output_secondary <- list(sets = susie_get_cs(susie_output, Xcorr = data_x, coverage = coverage, min_abs_corr = min_abs_corr), pip = susie_output$pip)
    susie_output_secondary$cs_corr <- get_cs_correlation(susie_output_secondary, Xcorr = data_x)
  }
  susie_output_secondary
}

#' Post-process SuSiE Analysis Results
#'
#' This function processes the results from SuSiE (Sum of Single Effects) genetic analysis.
#' It extracts and processes various statistics and indices based on the provided SuSiE object and other parameters.
#' The function can operate in 3 modes: 'susie', 'susie_rss', 'mvsusie', based on the method used for the SuSiE analysis.
#'
#' @param susie_output Output from running susieR::susie() or susieR::susie_rss() or mvsusieR::mvsusie()
#' @param data_x Genotype data matrix for 'susie' or Xcorr matrix for 'susie_rss'.
#' @param data_y Phenotype data vector for 'susie' or summary stats object for 'susie_rss' (a list contain attribute betahat and sebetahat AND/OR z). i.e. data_y = list(betahat = ..., sebetahat = ...), or NULL for mvsusie
#' @param X_scalar Scalar for the genotype data, used in residual scaling.
#' @param y_scalar Scalar for the phenotype data, used in residual scaling.
#' @param maf Minor Allele Frequencies vector.
#' @param secondary_coverage Vector of coverage thresholds for secondary conditional analysis.
#' @param signal_cutoff Cutoff value for signal identification in PIP values.
#' @param other_quantities A list of other quantities to be added to the final object.
#' @param prior_eff_tol Prior effective tolerance.
#' @param min_abs_corr Minimum absolute correlation for credible set purity filtering.
#'   Default is 0.8, which is stricter than the susieR default of 0.5. Credible sets
#'   with purity below this threshold are excluded from the results.
#' @param mode Specify the analysis mode: 'susie' or 'susie_rss'.
#' @return A list containing modified SuSiE object along with additional post-processing information.
#' @examples
#' # Example usage for SuSiE
#' # result <- susie_post_processor(susie_output, X_data, y_data, maf, mode = "susie")
#' # Example usage for SuSiE RSS
#' # result <- susie_post_processor(susie_output, Xcorr, z, maf, mode = "susie_rss")
#' @importFrom dplyr full_join
#' @importFrom purrr map_int pmap
#' @importFrom susieR get_cs_correlation susie_get_cs
#' @importFrom stringr str_replace
#' @export
susie_post_processor <- function(susie_output, data_x, data_y, X_scalar, y_scalar, maf = NULL,
                                 secondary_coverage = c(0.5, 0.7), signal_cutoff = 0.1,
                                 other_quantities = NULL, prior_eff_tol = 1e-9, min_abs_corr = 0.8,
                                 mode = c("susie", "susie_rss", "mvsusie")) {
  mode <- match.arg(mode)
  # Initialize result list
  res <- list(
    variant_names = normalize_variant_id(names(susie_output$pip))
  )
  analysis_script <- load_script()
  if (analysis_script != "") res$analysis_script <- analysis_script
  if (!is.null(other_quantities)) res$other_quantities <- other_quantities
  if (mode == "mvsusie") {
    res$context_names <- susie_output$condition_names
  }
  if (!is.null(data_y)) {
    # Mode-specific processing
    if (mode == "susie") {
      # Processing specific to susie_post_processor
      res$sumstats <- univariate_regression(data_x, data_y)
      y_scalar <- if (is.null(y_scalar) || all(y_scalar == 1)) 1 else y_scalar
      X_scalar <- if (is.null(X_scalar) || all(X_scalar == 1)) 1 else X_scalar
      res$sumstats$betahat <- res$sumstats$betahat * y_scalar / X_scalar
      res$sumstats$sebetahat <- res$sumstats$sebetahat * y_scalar / X_scalar
      res$sample_names <- rownames(data_y)
    } else if (mode == "susie_rss") {
      # Processing specific to susie_rss_post_processor
      res$sumstats <- data_y
    }
  }
  max_L <- nrow(susie_output$alpha)
  if (!is.null(susie_output$V)) {
    # for fSuSiE there is no V for now
    eff_idx <- which(susie_output$V > prior_eff_tol)
  } else {
    eff_idx <- 1:max_L
  }

  # Re-filter primary CS purity (susieR default is 0.5, pecotmr default is 0.8)
  if (mode %in% c("susie", "mvsusie")) {
    susie_output$sets <- susie_get_cs(susie_output, X = data_x, coverage = susie_output$sets$requested_coverage, min_abs_corr = min_abs_corr)
  } else {
    susie_output$sets <- susie_get_cs(susie_output, Xcorr = data_x, coverage = susie_output$sets$requested_coverage, min_abs_corr = min_abs_corr)
  }

  if (length(eff_idx) > 0) {
    # Prepare for top loci table
    top_variants_idx_pri <- get_top_variants_idx(susie_output, signal_cutoff)
    # get_cs_info returns data.frame(variant_idx, cs_idx) with one row per (variant, CS) pair
    top_loci_pri <- get_cs_info(susie_output$sets$cs, top_variants_idx_pri)
    if (is.null(top_loci_pri)) top_loci_pri <- data.frame(variant_idx = integer(0), cs_idx = integer(0))
    susie_output$cs_corr <- if (mode %in% c("susie", "mvsusie")) get_cs_correlation(susie_output, X = data_x) else get_cs_correlation(susie_output, Xcorr = data_x)
    top_loci_list <- list("coverage_0.95" = top_loci_pri)

    ## Loop over each secondary coverage value independently
    sets_secondary <- list()
    if (!is.null(secondary_coverage) && length(secondary_coverage)) {
      for (sec_cov in secondary_coverage) {
        sets_secondary[[paste0("coverage_", sec_cov)]] <- get_cs_and_corr(susie_output, sec_cov, data_x, mode, min_abs_corr)
        top_variants_idx_sec <- get_top_variants_idx(sets_secondary[[paste0("coverage_", sec_cov)]], signal_cutoff)
        top_loci_sec <- get_cs_info(sets_secondary[[paste0("coverage_", sec_cov)]]$sets$cs, top_variants_idx_sec)
        if (is.null(top_loci_sec)) top_loci_sec <- data.frame(variant_idx = integer(0), cs_idx = integer(0))
        top_loci_list[[paste0("coverage_", sec_cov)]] <- top_loci_sec
      }
    }

    # Merge coverage tables via full_join
    names(top_loci_list[[1]])[2] <- paste0("cs_", names(top_loci_list)[1])
    top_loci <- top_loci_list[[1]]
    if (length(top_loci_list) > 1) {
      for (i in 2:length(top_loci_list)) {
        names(top_loci_list[[i]])[2] <- paste0("cs_", names(top_loci_list)[i])
        top_loci <- dplyr::full_join(top_loci, top_loci_list[[i]], by = "variant_idx")
      }
    }

    if (nrow(top_loci) > 0) {
      top_loci[is.na(top_loci)] <- 0
      variants <- res$variant_names[top_loci$variant_idx]
      pip <- susie_output$pip[top_loci$variant_idx]
      top_loci_cols <- c("variant_id", if (!is.null(res$sumstats$betahat)) "betahat", if (!is.null(res$sumstats$sebetahat)) "sebetahat", if (!is.null(res$sumstats$z)) "z", if (!is.null(maf)) "maf", "pip", colnames(top_loci)[-1])
      res$top_loci <- data.frame(variants, stringsAsFactors = FALSE)
      res$top_loci$betahat <- if (!is.null(res$sumstats$betahat)) res$sumstats$betahat[top_loci$variant_idx] else NULL
      res$top_loci$sebetahat <- if (!is.null(res$sumstats$sebetahat)) res$sumstats$sebetahat[top_loci$variant_idx] else NULL
      res$top_loci$z <- if (!is.null(res$sumstats$z)) res$sumstats$z[top_loci$variant_idx] else NULL
      res$top_loci$maf <- if (!is.null(maf)) maf[top_loci$variant_idx] else NULL
      res$top_loci$pip <- pip
      res$top_loci <- cbind(res$top_loci, top_loci[, -1])
      colnames(res$top_loci) <- top_loci_cols
      rownames(res$top_loci) <- NULL
    }
    names(susie_output$pip) <- NULL
    res$susie_result_trimmed <- list(
      pip = susie_output$pip,
      sets = susie_output$sets,
      cs_corr = susie_output$cs_corr,
      sets_secondary = if (length(sets_secondary)) lapply(sets_secondary, function(x) x[names(x) != "pip"]) else NULL,
      alpha = susie_output$alpha[eff_idx, , drop = FALSE],
      lbf_variable = susie_output$lbf_variable[eff_idx, , drop = FALSE],
      V = if (!is.null(susie_output$V)) susie_output$V[eff_idx] else NULL,
      niter = susie_output$niter,
      max_L = max_L
    )
    if (mode == "susie") {
      res$susie_result_trimmed$X_column_scale_factors <- susie_output$X_column_scale_factors
      res$susie_result_trimmed$mu <- susie_output$mu[eff_idx, , drop = FALSE]
      res$susie_result_trimmed$mu2 <- susie_output$mu2[eff_idx, , drop = FALSE]
    }
    if (mode == "mvsusie") {
      # res$susie_result_trimmed$b1 = susie_output$b1[eff_idx, , , drop = FALSE]
      # res$susie_result_trimmed$b2 = susie_output$b2[eff_idx, , , drop = FALSE]
      res$susie_result_trimmed$b1_rescaled <- susie_output$b1_rescaled[eff_idx, , , drop = FALSE]
      res$susie_result_trimmed$coef <- susie_output$coef
      res$susie_result_trimmed$clfsr <- susie_output$conditional_lfsr[eff_idx, , , drop = FALSE]
      # other lfsr can be computed:
      # se_lfsr <- mvsusie_single_effect_lfsr(clfsr, alpha)
      # lfsr <- mvsusie_get_lfsr(clfsr, alpha)
    }
    class(res$susie_result_trimmed) <- "susie"
  }
  return(res)
}
