# =============================================================================
# Molecular trait imputation
#
# User-facing imputation methods for filling missing values in molecular trait
# matrices (expression, splicing, protein abundance) before fine-mapping and
# TWAS-weight training. Three method wrappers plus a unified pipeline; see the
# molecular-imputation vignette for when to use which method.
#
# The genotype-side mean_impute() helper lives in misc.R because it is an
# internal building block of compute_LD() and genotype loading, not a
# user-facing imputer.
# =============================================================================

#' XGBoost-based iterative imputation of missing values
#'
#' Imputes missing values in a numeric matrix by iteratively training
#' per-column XGBoost models on observed entries and predicting missing ones.
#' Columns that are entirely missing are removed. Initial imputation uses
#' column means.
#'
#' @param data Numeric matrix with missing values (NA).
#' @param maxiter Maximum number of imputation iterations (default 10).
#' @param max_depth Maximum tree depth for XGBoost (default 2).
#' @param nrounds Number of boosting rounds per variable (default 50).
#' @param decreasing Logical. If TRUE, impute variables with most missing
#'   values first. Default FALSE (fewest missing first).
#' @param num_workers Number of parallel workers for BiocParallel. Default 1
#'   (sequential).
#' @param verbose Logical, print progress (default TRUE).
#' @return The imputed matrix with the same dimensions as the input (minus
#'   any all-NA columns).
#' @importFrom BiocParallel MulticoreParam SerialParam bplapply
#' @export
xgboost_imputation <- function(data, maxiter = 10L, max_depth = 2L,
                                nrounds = 50L, decreasing = FALSE,
                                num_workers = 1L, verbose = TRUE) {
  if (!requireNamespace("xgboost", quietly = TRUE))
    stop("Package 'xgboost' is required for xgboost_imputation")

  xmis <- as.matrix(data)
  n <- nrow(xmis)
  p <- ncol(xmis)

  # Remove completely missing columns
  all_na <- colSums(is.na(xmis)) == n
  if (any(all_na)) {
    if (verbose)
      message("Removed ", sum(all_na), " column(s) with all entries missing.")
    xmis <- xmis[, !all_na, drop = FALSE]
    p <- ncol(xmis)
  }

  # Initial mean imputation
  ximp <- xmis
  col_means <- colMeans(xmis, na.rm = TRUE)
  for (j in seq_len(p)) {
    ximp[is.na(xmis[, j]), j] <- col_means[j]
  }

  # Missing value locations
  NAloc <- is.na(xmis)
  noNAvar <- colSums(NAloc)
  sort_j <- order(noNAvar, decreasing = decreasing)
  nzsort_j <- sort_j[noNAvar[sort_j] > 0]

  if (length(nzsort_j) == 0) {
    if (verbose) message("No missing values to impute.")
    return(ximp)
  }

  # Set up BiocParallel
  if (num_workers > 1L) {
    BPPARAM <- MulticoreParam(workers = num_workers)
  } else {
    BPPARAM <- SerialParam()
  }

  iter <- 0L
  conv_new <- 0
  conv_old <- Inf
  ximp_history <- vector("list", maxiter)

  while (conv_new < conv_old && iter < maxiter) {
    if (iter > 0) conv_old <- conv_new
    if (verbose) message("  XGBoost iteration ", iter + 1L, " in progress...")

    ximp_old <- ximp

    # Impute each variable with missing values
    impute_one <- function(var_idx) {
      obsi <- !NAloc[, var_idx]
      misi <- NAloc[, var_idx]
      obsY <- ximp[obsi, var_idx]
      obsX <- ximp[obsi, -var_idx, drop = FALSE]
      misX <- ximp[misi, -var_idx, drop = FALSE]

      xgb_train <- xgboost::xgb.DMatrix(data = obsX, label = obsY)
      xgb_pred <- xgboost::xgb.DMatrix(data = misX)
      model <- xgboost::xgb.train(
        params = list(max_depth = max_depth, verbosity = 0),
        data = xgb_train, nrounds = nrounds)
      list(var_idx = var_idx, predicted = predict(model, xgb_pred))
    }

    results <- bplapply(nzsort_j, impute_one, BPPARAM = BPPARAM)

    for (res in results) {
      misi <- NAloc[, res$var_idx]
      ximp[misi, res$var_idx] <- res$predicted
    }

    iter <- iter + 1L
    ximp_history[[iter]] <- ximp

    # Convergence: relative change in imputed values
    conv_new <- sum((ximp - ximp_old)^2) / sum(ximp^2)
  }

  # Return last improving iteration
  if (iter == maxiter) ximp_history[[iter]] else ximp_history[[max(iter - 1L, 1L)]]
}

#' softImpute-based matrix imputation
#'
#' Imputes missing values in a numeric matrix via the soft-thresholded SVD
#' algorithm of \code{softImpute::softImpute}. Uses the \code{type = "als"}
#' alternating-least-squares solver by default.
#'
#' \code{lambda} drives the bias-variance trade-off. By default the
#' function tunes it via held-out cross-validation: a random fraction of
#' observed cells is masked, softImpute is fit across a geometric grid of
#' lambdas with warm-start, and the lambda minimising RMSE on the masked
#' cells is selected. The final imputation then uses that lambda on the
#' full input. Pass a scalar \code{lambda} to skip tuning.
#'
#' Columns with all entries missing are removed before fitting and a
#' one-line message records the count.
#'
#' @param data Numeric matrix (samples x features) with missing values
#'   (NA).
#' @param rank.max Maximum rank of the SVD solution. Default
#'   \code{min(dim(data) - 1)}.
#' @param lambda Nuclear-norm penalty. \code{NULL} (default) tunes over
#'   a default geometric grid via held-out CV. A scalar uses that lambda
#'   directly with no tuning. A vector defines the tuning grid.
#' @param tune_fraction Fraction of observed cells held out for CV.
#'   Default 0.05. Ignored when \code{lambda} is a scalar.
#' @param tune_seed Seed for the held-out mask. Default 42L.
#' @param type Optimisation type passed to \code{softImpute::softImpute};
#'   one of \code{"als"} (default, faster) or \code{"svd"}.
#' @param thresh Convergence threshold. Default 1e-5.
#' @param maxit Maximum iterations. Default 100.
#' @param verbose Logical, print progress (including per-lambda CV RMSE).
#'   Default FALSE.
#' @return The imputed matrix with the same dimensions as the input
#'   (minus any all-NA columns). Carries attributes
#'   \code{"softimpute_lambda"} (the lambda used for the final fit) and,
#'   when tuning was performed, \code{"softimpute_lambda_cv"} (a
#'   data.frame of grid-search RMSE values).
#' @export
softimpute_imputation <- function(data, rank.max = NULL, lambda = NULL,
                                   tune_fraction = 0.05,
                                   tune_seed = 42L,
                                   type = c("als", "svd"),
                                   thresh = 1e-5, maxit = 100,
                                   verbose = FALSE) {
  if (!requireNamespace("softImpute", quietly = TRUE)) {
    stop("Package 'softImpute' is required for softimpute_imputation.")
  }
  type <- match.arg(type)
  X <- as.matrix(data)
  n <- nrow(X); p <- ncol(X)
  all_na <- colSums(is.na(X)) == n
  if (any(all_na)) {
    if (verbose) message("Removed ", sum(all_na), " column(s) with all entries missing.")
    X <- X[, !all_na, drop = FALSE]
    p <- ncol(X)
  }
  if (is.null(rank.max)) rank.max <- max(1L, min(n, p) - 1L)
  lambda0_val <- softImpute::lambda0(X)

  # Decide whether to tune. lambda = scalar => no tuning; NULL/vector => tune.
  tune <- is.null(lambda) || length(lambda) > 1L
  if (tune) {
    if (is.null(lambda)) {
      # Default geometric grid from lambda0 down to 1e-3 * lambda0.
      lambda_grid <- exp(seq(log(lambda0_val),
                              log(max(1e-3 * lambda0_val, .Machine$double.eps)),
                              length.out = 10L))
    } else {
      lambda_grid <- sort(as.numeric(lambda), decreasing = TRUE)
    }
    cv_res <- .softimpute_cv_lambda(X, lambda_grid = lambda_grid,
                                     rank.max = rank.max,
                                     tune_fraction = tune_fraction,
                                     tune_seed = tune_seed,
                                     type = type, thresh = thresh, maxit = maxit,
                                     verbose = verbose)
    lambda <- cv_res$best_lambda
    if (verbose) {
      message(sprintf("softimpute_imputation: tuned lambda = %.4g (CV RMSE %.4f over %d cells)",
                      lambda, cv_res$best_rmse, cv_res$n_held_out))
    }
  }

  fit <- softImpute::softImpute(X, rank.max = rank.max, lambda = lambda,
                                 type = type, thresh = thresh, maxit = maxit,
                                 trace.it = isTRUE(verbose))
  Ximp <- softImpute::complete(X, fit)
  if (is.null(rownames(Ximp))) rownames(Ximp) <- rownames(X)
  if (is.null(colnames(Ximp))) colnames(Ximp) <- colnames(X)
  attr(Ximp, "softimpute_lambda") <- lambda
  if (tune) attr(Ximp, "softimpute_lambda_cv") <- cv_res$grid
  Ximp
}

# Internal: CV-tune softImpute's lambda over a grid using a held-out mask.
# Uses softImpute's warm.start feature to chain fits along the lambda path
# (largest lambda first). Returns the chosen lambda, its CV RMSE, and the
# full grid scores.
.softimpute_cv_lambda <- function(X, lambda_grid, rank.max,
                                   tune_fraction, tune_seed,
                                   type, thresh, maxit, verbose) {
  na_mask <- is.na(X)
  obs_idx <- which(!na_mask)
  n_hold <- max(1L, floor(tune_fraction * length(obs_idx)))
  if (n_hold >= length(obs_idx)) {
    stop("tune_fraction too large; would mask all observed entries.")
  }
  set.seed(tune_seed)
  hold_idx <- sample(obs_idx, size = n_hold, replace = FALSE)
  truth <- X[hold_idx]
  Xmasked <- X
  Xmasked[hold_idx] <- NA_real_

  lambdas <- sort(as.numeric(lambda_grid), decreasing = TRUE)
  rmse_vec <- numeric(length(lambdas))
  warm <- NULL
  for (i in seq_along(lambdas)) {
    fit_args <- list(x = Xmasked, rank.max = rank.max, lambda = lambdas[i],
                     type = type, thresh = thresh, maxit = maxit,
                     trace.it = FALSE)
    if (!is.null(warm)) fit_args$warm.start <- warm
    fit <- tryCatch(do.call(softImpute::softImpute, fit_args),
                    error = function(e) NULL)
    if (is.null(fit)) {
      rmse_vec[i] <- NA_real_
      next
    }
    Ximp <- softImpute::complete(Xmasked, fit)
    rmse_vec[i] <- sqrt(mean((Ximp[hold_idx] - truth)^2, na.rm = TRUE))
    warm <- fit
    if (verbose) {
      message(sprintf("  softimpute CV lambda = %.4g -> RMSE %.4f",
                      lambdas[i], rmse_vec[i]))
    }
  }
  grid <- data.frame(lambda = lambdas, rmse = rmse_vec,
                     stringsAsFactors = FALSE)
  ok <- which(is.finite(rmse_vec))
  if (length(ok) == 0L) {
    stop("All softImpute CV fits failed; supply lambda explicitly.")
  }
  best <- ok[which.min(rmse_vec[ok])]
  list(best_lambda = lambdas[best], best_rmse = rmse_vec[best],
       n_held_out = n_hold, grid = grid)
}

#' flashier-based matrix imputation
#'
#' Imputes missing values in a numeric matrix via empirical Bayes
#' matrix factorisation (\code{flashier::flash}). \code{flash} natively
#' supports NA entries and reconstructs them from the posterior fitted
#' values. The default prior families match those used by
#' \code{\link{compute_cov_flash}}: \code{ebnm::ebnm_normal} on the row
#' factor and \code{ebnm::ebnm_normal_scale_mixture} on the column
#' factor.
#'
#' @param data Numeric matrix (samples x features) with missing values
#'   (NA).
#' @param var_type Variance structure for the residuals, passed to
#'   \code{flashier::flash}. Default 2 (per-column).
#' @param ebnm_fn Prior family. Default
#'   \code{c(ebnm::ebnm_normal, ebnm::ebnm_normal_scale_mixture)}.
#' @param greedy_Kmax Maximum number of greedy factors. Default
#'   \code{min(20L, min(dim(data)) - 1L)}.
#' @param backfit Logical; run a final backfit pass. Default TRUE.
#' @param verbose Numeric verbosity (0 silent, default).
#' @param ... Additional arguments forwarded to \code{flashier::flash}.
#' @return The imputed matrix with the same dimensions as the input
#'   (minus any all-NA columns).
#' @export
flashier_imputation <- function(data, var_type = 2,
                                 ebnm_fn = c(ebnm::ebnm_normal,
                                             ebnm::ebnm_normal_scale_mixture),
                                 greedy_Kmax = NULL,
                                 backfit = TRUE, verbose = 0, ...) {
  if (!requireNamespace("flashier", quietly = TRUE)) {
    stop("Package 'flashier' is required for flashier_imputation.")
  }
  if (!requireNamespace("ebnm", quietly = TRUE)) {
    stop("Package 'ebnm' is required for flashier_imputation.")
  }
  X <- as.matrix(data)
  n <- nrow(X); p <- ncol(X)
  all_na <- colSums(is.na(X)) == n
  if (any(all_na)) {
    if (isTRUE(verbose) || (is.numeric(verbose) && verbose >= 1)) {
      message("Removed ", sum(all_na), " column(s) with all entries missing.")
    }
    X <- X[, !all_na, drop = FALSE]
    p <- ncol(X)
  }
  if (is.null(greedy_Kmax)) greedy_Kmax <- max(1L, min(20L, min(n, p) - 1L))
  fl <- flashier::flash(X, var_type = var_type, ebnm_fn = ebnm_fn,
                         greedy_Kmax = greedy_Kmax, backfit = backfit,
                         verbose = verbose, ...)
  # Posterior fitted value: L_pm %*% t(F_pm). When no factors are found,
  # fall back to column means.
  if (is.null(fl$n_factors) || fl$n_factors == 0) {
    fitted <- matrix(colMeans(X, na.rm = TRUE),
                     nrow = n, ncol = p, byrow = TRUE)
  } else {
    fitted <- tcrossprod(fl$L_pm, fl$F_pm)
  }
  Ximp <- X
  na_idx <- is.na(X)
  Ximp[na_idx] <- fitted[na_idx]
  if (is.null(rownames(Ximp))) rownames(Ximp) <- rownames(X)
  if (is.null(colnames(Ximp))) colnames(Ximp) <- colnames(X)
  Ximp
}

#' Pipeline for molecular trait imputation
#'
#' Unified entry point for imputing missing values in a molecular trait
#' matrix (e.g. expression, splicing, protein abundance). Dispatches to
#' one or more underlying methods and optionally evaluates imputation
#' quality on a held-out random mask of observed entries.
#'
#' Available methods:
#' \itemize{
#'   \item \code{"xgboost"} - \code{\link{xgboost_imputation}}, iterative
#'     XGBoost per-column models. No structural assumption on the matrix.
#'   \item \code{"softimpute"} - \code{\link{softimpute_imputation}},
#'     soft-thresholded SVD via \code{softImpute::softImpute}. Assumes a
#'     low-rank smooth signal.
#'   \item \code{"flashier"} - \code{\link{flashier_imputation}},
#'     empirical Bayes matrix factorisation via \code{flashier::flash}
#'     with sparse + scale-mixture priors.
#' }
#'
#' @param data Numeric matrix (samples x features) with missing values.
#' @param methods Character vector of methods to fit. Any subset of
#'   \code{c("xgboost", "softimpute", "flashier")}. Default fits all
#'   three.
#' @param method_args Named list of per-method argument lists. List
#'   names must match \code{methods}. Default empty.
#' @param evaluate Logical. If TRUE, mask a random fraction of observed
#'   entries, run each method, and report root-mean-squared error on the
#'   masked cells. Default FALSE.
#' @param eval_fraction Fraction of observed cells to hold out for
#'   evaluation. Default 0.05.
#' @param eval_seed Random seed for the held-out mask. Default 42.
#' @param verbose Logical, print progress. Default FALSE.
#'
#' @return A list with:
#' \describe{
#'   \item{imputed}{Named list of imputed matrices, one per requested
#'     method.}
#'   \item{missingness}{Per-column NA count and overall NA fraction.}
#'   \item{evaluation}{When \code{evaluate = TRUE}, a data.frame with
#'     one row per method: \code{method}, \code{rmse}, \code{n_held_out},
#'     \code{time_seconds}.}
#'   \item{method_args}{The resolved per-method arguments used.}
#' }
#' @export
molecular_trait_imputation_pipeline <- function(
    data,
    methods = c("xgboost", "softimpute", "flashier"),
    method_args = list(),
    evaluate = FALSE,
    eval_fraction = 0.05,
    eval_seed = 42L,
    verbose = FALSE) {
  valid <- c("xgboost", "softimpute", "flashier")
  if (!is.character(methods) || length(methods) == 0L) {
    stop("methods must be a non-empty character vector.")
  }
  bad <- setdiff(methods, valid)
  if (length(bad) > 0) {
    stop("Unknown imputation method(s): ", paste(bad, collapse = ", "),
         ". Valid options: ", paste(valid, collapse = ", "))
  }
  methods <- unique(methods)

  Xfull <- as.matrix(data)
  n <- nrow(Xfull); p <- ncol(Xfull)
  na_mask_full <- is.na(Xfull)
  na_per_col <- colSums(na_mask_full)
  missingness <- list(
    per_column = na_per_col,
    n_missing = sum(na_mask_full),
    fraction = sum(na_mask_full) / (n * p)
  )
  if (verbose) {
    message(sprintf("molecular_trait_imputation_pipeline: %d x %d matrix, %.2f%% missing.",
                    n, p, 100 * missingness$fraction))
  }

  dispatch <- list(xgboost = xgboost_imputation,
                   softimpute = softimpute_imputation,
                   flashier = flashier_imputation)

  run_one <- function(m, X_in) {
    args <- method_args[[m]] %||% list()
    do.call(dispatch[[m]], c(list(X_in), args))
  }

  # If evaluating, hold out a random subset of observed cells, run each
  # method against the masked matrix, then score RMSE on the held-out
  # cells. Final imputation is fit against the original matrix so the
  # returned matrices use every observed cell.
  eval_df <- NULL
  if (isTRUE(evaluate)) {
    obs_idx <- which(!na_mask_full)
    n_hold <- max(1L, floor(eval_fraction * length(obs_idx)))
    if (n_hold >= length(obs_idx)) {
      stop("eval_fraction too large; would mask all observed entries.")
    }
    set.seed(eval_seed)
    hold_idx <- sample(obs_idx, size = n_hold, replace = FALSE)
    Xmasked <- Xfull
    held_truth <- Xfull[hold_idx]
    Xmasked[hold_idx] <- NA_real_
    eval_records <- vector("list", length(methods))
    for (i in seq_along(methods)) {
      m <- methods[i]
      if (verbose) message("Evaluating ", m, " on held-out mask ...")
      t0 <- proc.time()
      Ximp_m <- run_one(m, Xmasked)
      t1 <- proc.time()
      pred <- Ximp_m[hold_idx]
      rmse <- sqrt(mean((pred - held_truth)^2, na.rm = TRUE))
      eval_records[[i]] <- data.frame(
        method = m, rmse = rmse, n_held_out = n_hold,
        time_seconds = as.numeric((t1 - t0)["elapsed"]),
        stringsAsFactors = FALSE
      )
    }
    eval_df <- do.call(rbind, eval_records)
    rownames(eval_df) <- NULL
  }

  # Final imputation: run each method on the full matrix
  imputed <- list()
  for (m in methods) {
    if (verbose) message("Imputing with ", m, " ...")
    imputed[[m]] <- run_one(m, Xfull)
  }

  list(
    imputed = imputed,
    missingness = missingness,
    evaluation = eval_df,
    method_args = method_args
  )
}
