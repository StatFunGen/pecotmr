#' Multivariate Analysis Pipeline
#'
#' This function performs weights computation for Transcriptome-Wide Association Study (TWAS) with fitting
#' models using mvSuSiE and mr.mash with the option of using a limited number of variants selected from
#' mvSuSiE fine-mapping for computing TWAS weights with cross-validation.
#'
#' @param X A matrix of genotype data where rows represent samples and columns represent genetic variants.
#' @param Y A matrix of phenotype measurements, representing samples and columns represent conditions.
#' @param maf A list of vectors for minor allele frequencies for each variant in X.
#' @param max_L The maximum number of components in mvSuSiE. Default is 30.
#' @param ld_reference_meta_file An optional path to a file containing linkage disequilibrium reference data. If provided, variants in X are filtered based on this reference.
#' @param pip_cutoff_to_skip Cutoff value for skipping conditions based on PIP values. Default is 0.
#' @param signal_cutoff Cutoff value for signal identification in PIP values for susie_post_processor. Default is 0.025.
#' @param coverage A vector of coverage probabilities, with the first element being the primary coverage and the rest being secondary coverage probabilities for credible set refinement. Defaults to c(0.95, 0.7, 0.5).
#' @param min_abs_corr Minimum absolute correlation for credible set purity filtering. Default is 0.8,
#'   which is stricter than the susieR default of 0.5.
#' @param data_driven_prior_matrices A list of data-driven covariance matrices for mr.mash weights.
#' @param data_driven_prior_matrices_cv A list of data-driven covariance matrices for mr.mash weights in cross-validation.
#' @param canonical_prior_matrices If set to TRUE, will compute canonical covariance matrices and add them into the prior covariance matrix list in mrmash_wrapper. Default is TRUE.
#' @param sample_partition Sample partition for cross-validation.
#' @param mrmash_max_iter The maximum number of iterations for mr.mash. Default is 5000.
#' @param mvsusie_max_iter The maximum number of iterations for mvSuSiE. Default is 200.
#' @param min_cv_maf The minimum minor allele frequency for variants to be included in cross-validation. Default is 0.05.
#' @param max_cv_variants The maximum number of variants to be included in cross-validation. Defaults to -1 which means no limit.
#' @param cv_folds The number of folds to use for cross-validation. Set to 0 to skip cross-validation. Default is 5.
#' @param cv_threads The number of threads to use for parallel computation in cross-validation. Defaults to 1.
#' @param data_driven_prior_weights_cutoff The minimum weight for prior covariance matrices. Default is 1e-4.
#' @param verbose Verbosity level. Default is 0.
#'
#' @return A list containing the multivariate analysis results.
#' @examples
#' library(pecotmr)
#'
#' data(multitrait_data)
#' attach(multitrait_data)
#'
#' data_driven_prior_matrices <- list(
#'   U = prior_matrices,
#'   w = rep(1 / length(prior_matrices), length(prior_matrices))
#' )
#'
#' data_driven_prior_matrices_cv <- lapply(prior_matrices_cv, function(x) {
#'   list(U = x, w = rep(1 / length(x), length(x)))
#' })
#'
#' result <- multivariate_analysis_pipeline(
#'   X = multitrait_data$X,
#'   Y = multitrait_data$Y,
#'   maf = colMeans(multitrait_data$X),
#'   X_variance = multitrait_data$X_variance,
#'   max_L = 10,
#'   ld_reference_meta_file = NULL,
#'   max_cv_variants = -1,
#'   pip_cutoff_to_skip = 0,
#'   signal_cutoff = 0.025,
#'   data_driven_prior_matrices = data_driven_prior_matrices,
#'   data_driven_prior_matrices_cv = data_driven_prior_matrices_cv,
#'   canonical_prior_matrices = TRUE,
#'   sample_partition = NULL,
#'   cv_folds = 5,
#'   cv_threads = 2,
#'   data_driven_prior_weights_cutoff = 1e-4
#' )
#' @export
multivariate_analysis_pipeline <- function(
    # input data
    X,
    Y,
    maf,
    X_variance = NULL,
    other_quantities = list(),
    # filters
    imiss_cutoff = 1.0,
    maf_cutoff = 0.01,
    xvar_cutoff = 0.01,
    ld_reference_meta_file = NULL,
    pip_cutoff_to_skip = 0,
    # methods parameter configuration
    max_L = -1,
    data_driven_prior_matrices = NULL,
    data_driven_prior_matrices_cv = NULL,
    data_driven_prior_weights_cutoff = 1e-4,
    canonical_prior_matrices = TRUE,
    mrmash_max_iter = 5000,
    mvsusie_max_iter = 200,
    # fine-mapping results summary
    signal_cutoff = 0.025,
    coverage = c(0.95, 0.7, 0.5),
    min_abs_corr = 0.8,
    # TWAS weights and CV for TWAS weights
    twas_weights = TRUE,
    sample_partition = NULL,
    max_cv_variants = -1,
    cv_folds = 5,
    cv_threads = 1,
    verbose = 0) {
  # Make sure mvsusieR is installed
  if (!requireNamespace("mvsusieR", quietly = TRUE)) {
    stop("To use this function, please install mvsusieR: https://github.com/stephenslab/mvsusieR")
  }
  # Skip conditions based on univariate PIP values
  skip_conditions <- function(X, Y, pip_cutoff_to_skip) {
    if (length(pip_cutoff_to_skip) == 1 && is.numeric(pip_cutoff_to_skip)) {
      pip_cutoff_to_skip <- rep(pip_cutoff_to_skip, ncol(Y))
    } else if (length(pip_cutoff_to_skip) != ncol(Y)) {
      stop("pip_cutoff_to_skip must be a single number or a vector of the same length as ncol(Y).")
    }
    cols_to_keep <- logical(ncol(Y))
    for (r in 1:ncol(Y)) {
      if (pip_cutoff_to_skip[r] != 0) {
        non_missing_indices <- which(!is.na(Y[, r]))
        X_non_missing <- X[match(names(Y[, r])[non_missing_indices], rownames(X)), ]
        Y_non_missing <- Y[non_missing_indices, r]
        if (pip_cutoff_to_skip[r] < 0) {
          # automatically determine the cutoff to use
          pip_cutoff_to_skip[r] <- 3 * 1 / ncol(X_non_missing)
        }
        top_model_pip <- susie(X_non_missing, Y_non_missing, L = 1)$pip

        if (any(top_model_pip > pip_cutoff_to_skip[r])) {
          cols_to_keep[r] <- TRUE
        } else {
          message(paste0(
            "Skipping condition ", colnames(Y)[r], ", because all top_model_pip < pip_cutoff_to_skip = ",
            pip_cutoff_to_skip[r], ". Top loci model does not show any potentially significant variants."
          ))
        }
      } else {
        cols_to_keep[r] <- TRUE
      }
    }

    Y_filtered <- Y[, cols_to_keep, drop = FALSE]

    if (ncol(Y_filtered) <= 1) {
      warning("After filtering by potential association signals, Y has ", ncol(Y_filtered), " context left. Returning NULL.")
      return(NULL)
    } else {
      message("After filtering by potential association signals, Y has ", ncol(Y_filtered), " contexts left.")
      return(Y_filtered)
    }
  }

  initialize_mvsusie_prior <- function(condition_names, data_driven_prior_matrices,
                                       data_driven_prior_matrices_cv, cv_folds, prior_weights, data_driven_prior_weights_cutoff) {
    if (!is.null(data_driven_prior_matrices)) {
      # update w based on mrmash prior weights
      message("Updating prior weights based on mrmash_fitted. ")
      data_driven_prior_matrices$w <- prior_weights
      data_driven_prior_matrices$U <- data_driven_prior_matrices$U[names(prior_weights)]
      data_driven_prior_matrices <- list(matrices = data_driven_prior_matrices$U, weights = data_driven_prior_matrices$w)
      data_driven_prior_matrices <- mvsusieR::create_mixture_prior(mixture_prior = data_driven_prior_matrices, weights_tol = data_driven_prior_weights_cutoff, include_indices = condition_names)
    } else {
      data_driven_prior_matrices <- mvsusieR::create_mixture_prior(R = length(condition_names), include_indices = condition_names)
    }

    if (!is.null(data_driven_prior_matrices_cv)) {
      data_driven_prior_matrices_cv <- lapply(
        data_driven_prior_matrices_cv,
        function(x) {
          x$U <- x$U[names(prior_weights)]
          x <- list(matrices = x$U, weights = prior_weights)
          mvsusieR::create_mixture_prior(mixture_prior = x, weights_tol = data_driven_prior_weights_cutoff, include_indices = condition_names)
        }
      )
    } else {
      if (!is.null(data_driven_prior_matrices)) {
        data_driven_prior_matrices_cv <- lapply(1:cv_folds, function(x) {
          return(data_driven_prior_matrices)
        })
      }
    }
    return(list(
      data_driven_prior_matrices = data_driven_prior_matrices, data_driven_prior_matrices_cv = data_driven_prior_matrices_cv
    ))
  }

  # filter X and Y missing, specific to multivariate analysis where some conditions are skipped we have to updated X matrix
  filter_X_Y_missing <- function(X, Y) {
    Y_rows_with_missing <- apply(Y, 1, function(row) all(is.na(row)))
    if (any(Y_rows_with_missing)) {
      Y_filtered <- Y[-which(Y_rows_with_missing), , drop = FALSE]
    } else {
      Y_filtered <- Y
    }
    X_filtered <- X[match(rownames(Y_filtered), rownames(X)), ]
    X_columns_with_missing <- apply(X_filtered, 2, function(column) all(is.na(column)))
    if (any(X_columns_with_missing)) {
      columns_to_remove <- which(X_columns_with_missing)
      X_filtered <- X_filtered[, -columns_to_remove, drop = FALSE]
    }
    return(list(X_filtered = X_filtered, Y_filtered = Y_filtered))
  }

  # Input validation
  if (!is.matrix(X) || !is.numeric(X)) stop("X must be a numeric matrix")
  if (!is.matrix(Y) || !is.numeric(Y)) stop("Y must be a numeric matrix")
  if (nrow(X) != nrow(Y)) stop("X and Y must have the same number of rows")
  if (!is.numeric(maf) || length(maf) != ncol(X)) stop("maf must be a numeric vector with length equal to the number of columns in X")
  if (any(maf < 0 | maf > 1)) stop("maf values must be between 0 and 1")

  # main analysis codes
  Y <- skip_conditions(X, Y, pip_cutoff_to_skip)
  if (is.null(Y)) {
    return(list())
  }

  # filter X and Y missing data
  X_Y_filtered <- filter_X_Y_missing(X, Y)
  X <- X_Y_filtered$X_filtered
  Y <- X_Y_filtered$Y_filtered
  if (nrow(Y) == 0 || is.null(Y)) {
    return(list())
  }

  # filter variants by ld reference panel
  if (!is.null(ld_reference_meta_file)) {
    variants_kept <- filter_variants_by_ld_reference(colnames(X), ld_reference_meta_file)
    X <- X[, variants_kept$data, drop = FALSE]
    maf <- maf[variants_kept$idx]
  }

  # filter X based on Y subjects
  if (!is.null(imiss_cutoff) || !is.null(maf_cutoff)) {
    X <- filter_X_with_Y(X, Y, imiss_cutoff, maf_cutoff, var_thresh = xvar_cutoff, maf = maf, X_variance = X_variance)
    maf <- maf[colnames(X)]
  }

  # filter data driven prior matrices
  if (!is.null(data_driven_prior_matrices)) {
    data_driven_prior_matrices <- filter_mixture_components(
      colnames(Y),
      data_driven_prior_matrices$U, data_driven_prior_matrices$w,
      data_driven_prior_weights_cutoff
    )
  }

  st <- proc.time()
  res <- list()
  message("Fitting mr.mash model on input data ...")
  res$mrmash_fitted <- mrmash_wrapper(
    X = X, Y = Y, data_driven_prior_matrices = data_driven_prior_matrices,
    canonical_prior_matrices = canonical_prior_matrices, max_iter = mrmash_max_iter
  )

  # For input into mvSuSiE
  resid_Y <- res$mrmash_fitted$V
  w0_updated <- rescale_cov_w0(res$mrmash_fitted$w0)
  if (length(w0_updated) == 0) {
    return(list())
  }
  w0_updated <- w0_updated[names(w0_updated) %in% names(data_driven_prior_matrices$U)]
  data_driven_prior_matrices$U <- data_driven_prior_matrices$U[names(w0_updated)]
  data_driven_prior_matrices$w <- data_driven_prior_matrices$w[names(w0_updated)]

  if (!is.null(data_driven_prior_matrices_cv)) {
    for (fold in seq_along(data_driven_prior_matrices_cv)) {
      data_driven_prior_matrices_cv[[fold]] <- filter_mixture_components(
        colnames(Y), data_driven_prior_matrices_cv[[fold]]$U,
        data_driven_prior_matrices_cv[[fold]]$w, data_driven_prior_weights_cutoff
      )
      data_driven_prior_matrices_cv[[fold]]$w <- data_driven_prior_matrices_cv[[fold]]$w[names(data_driven_prior_matrices_cv[[fold]]$w) %in% names(w0_updated)]
      data_driven_prior_matrices_cv[[fold]]$w <- w0_updated[names(data_driven_prior_matrices_cv[[fold]]$w)]
      data_driven_prior_matrices_cv[[fold]]$U <- data_driven_prior_matrices_cv[[fold]]$U[names(data_driven_prior_matrices_cv[[fold]]$U) %in% names(w0_updated)]
    }
  } else if (is.null(data_driven_prior_matrices_cv) && !is.null(data_driven_prior_matrices)) {
    data_driven_prior_matrices_cv <- lapply(1:cv_folds, function(fold) data_driven_prior_matrices)
    names(data_driven_prior_matrices_cv) <- paste0("fold_", 1:cv_folds)
  }

  if (max_L < 0) {
    # This is based on mr.mash fit
    # which can be a huge overestimate
    # so we bound it between 5 and 20
    max_L <- min(20, max(5, sum(1 - res$mrmash_fitted$w1[, 1])))
  }

  mvsusie_reweighted_mixture_prior <- initialize_mvsusie_prior(
    colnames(Y), data_driven_prior_matrices,
    data_driven_prior_matrices_cv, cv_folds, w0_updated, data_driven_prior_weights_cutoff
  )
  res$reweighted_mixture_prior <- mvsusie_reweighted_mixture_prior$data_driven_prior_matrices
  res$reweighted_mixture_prior_cv <- mvsusie_reweighted_mixture_prior$data_driven_prior_matrices_cv

  # Fit mvSuSiE
  message("Fitting mvSuSiE model on input data ...")
  res$mvsusie_fitted <- mvsusieR::mvsusie(X, Y,
    L = max_L, prior_variance = mvsusie_reweighted_mixture_prior$data_driven_prior_matrices,
    residual_variance = resid_Y, precompute_covariances = FALSE, compute_objective = TRUE,
    estimate_residual_variance = FALSE, estimate_prior_variance = TRUE, estimate_prior_method = "EM",
    max_iter = mvsusie_max_iter, n_thread = 1, approximate = FALSE, verbosity = verbose, coverage = coverage[1]
  )

  # Process mvSuSiE results
  sec_coverage <- if (length(coverage) > 1) coverage[-1] else NULL
  mvsusie_result_trimmed <- susie_post_processor(
    res$mvsusie_fitted, X, NULL, 1, 1,
    maf = maf, secondary_coverage = sec_coverage, signal_cutoff = signal_cutoff,
    min_abs_corr = min_abs_corr, mode = "mvsusie", other_quantities = other_quantities
  )
  res <- c(res, mvsusie_result_trimmed)
  res$total_time_elapsed <- proc.time() - st

  # Run TWAS weights and optionally CV
  if (twas_weights) {
    res$twas_weights_result <- twas_multivariate_weights_pipeline(X, Y, res,
      cv_folds = cv_folds, sample_partition = sample_partition,
      max_cv_variants = max_cv_variants,
      mvsusie_max_iter = mvsusie_max_iter, mrmash_max_iter = mrmash_max_iter,
      canonical_prior_matrices = canonical_prior_matrices, data_driven_prior_matrices = data_driven_prior_matrices,
      data_driven_prior_matrices_cv = data_driven_prior_matrices_cv,
      cv_threads = cv_threads, verbose = verbose
    )
  }
  return(res)
}
