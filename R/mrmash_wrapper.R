#' @title Mr.Mash Wrapper
#'
#' @description Compute weights with mr.mash using a precomputed prior grid and mixture prior.
#'
#' @param X An n x p matrix of genotype data, where n is the total number of individuals and p is the number of SNPs.
#' @param Y An n x r matrix of residual expression data, where n is the total number of individuals and r is the total number of conditions (tissue/cell-types).
#' @param data_driven_prior_matrices A list of data-driven covariance matrices. Default is NULL.
#' @param prior_grid A vector of scaling factors to be used in fitting the mr.mash model. Default is NULL.
#' @param nthreads The number of threads to use for parallel computation. Default is 2.
#' @param canonical_prior_matrices A logical indicating whether to use canonical matrices as priors. Default is FALSE.
#' @param standardize A logical indicating whether to standardize the input data. Default is FALSE.
#' @param update_w0 A logical indicating whether to update the prior mixture weights. Default is TRUE.
#' @param w0_threshold The threshold for updating prior mixture weights. Default is 1e-8.
#' @param update_V A logical indicating whether to update the residual covariance matrix. Default is TRUE.
#' @param update_V_method The method for updating the residual covariance matrix. Default is "full".
#' @param B_init_method The method for initializing the coefficient matrix. Default is "enet".
#' @param max_iter The maximum number of iterations. Default is 5000.
#' @param tol The tolerance for convergence. Default is 0.01.
#' @param verbose A logical indicating whether to print verbose output. Default is FALSE.
#' @param ... Additional arguments to be passed to mr.mash.
#'
#' @return A mr.mash fit, stored as a list with some or all of the following elements:
#' \item{mu1}{A p x r matrix of posterior means for the regression coefficients.}
#' \item{S1}{An r x r x p array of posterior covariances for the regression coefficients.}
#' \item{w1}{A p x K matrix of posterior assignment probabilities to the mixture components.}
#' \item{V}{An r x r residual covariance matrix.}
#' \item{w0}{A K-vector with (updated, if \code{update_w0=TRUE}) prior mixture weights, each associated with the respective covariance matrix in \code{S0}.}
#' \item{S0}{An r x r x K array of prior covariance matrices on the regression coefficients.}
#' \item{intercept}{An r-vector containing the posterior mean estimate of the intercept.}
#' \item{fitted}{An n x r matrix of fitted values.}
#' \item{G}{An r x r covariance matrix of fitted values.}
#' \item{pve}{An r-vector of proportion of variance explained by the covariates.}
#' \item{ELBO}{The Evidence Lower Bound (ELBO) at the last iteration.}
#' \item{progress}{A data frame including information regarding convergence criteria at each iteration.}
#' \item{converged}{A logical indicating whether the optimization algorithm converged to a solution within the chosen tolerance level.}
#' \item{elapsed_time}{The computation runtime for fitting mr.mash.}
#' \item{Y}{An n x r matrix of responses at the last iteration (only relevant when missing values are present in the input Y).}
#'
#' @examples
#' set.seed(123)
#' prior_grid <- runif(17, 0.00005, 0.05)
#'
#' sample_id <- paste0("P000", str_pad(1:400, 3, pad = "0"))
#' X <- matrix(sample(0:2, size = n * p, replace = TRUE, prob = c(0.65, 0.30, 0.05)), nrow = n)
#' rownames(X) <- sample_id
#' colnames(X) <- paste0("rs", sample(10000:100000, p))
#'
#' tissues <- c(
#'   "Adipose Tissue", "Muscle Tissue", "Brain Tissue", "Liver Tissue",
#'   "Kidney Tissue", "Heart Tissue", "Lung Tissue"
#' )
#' Y <- matrix(runif(n * r, -2, 2), nrow = n)
#' Y <- scale(Y)
#' colnames(Y) <- tissues
#' rownames(Y) <- sample_id
#'
#' set.seed(Sys.time())
#' components <- c(
#'   "XtX", "tFLASH_default", "FLASH_default", "tFLASH_nonneg",
#'   "FLASH_nonneg", "PCA"
#' )
#'
#' data_driven_prior_matrices <- list()
#' for (i in components) {
#'   A <- matrix(runif(r^2) * 2 - 1, ncol = r)
#'   cov <- t(A) %*% A
#'   colnames(cov) <- tissues
#'   rownames(cov) <- tissues
#'   data_driven_prior_matrices[[i]] <- cov
#' }
#'
#' res <- mrmash_wrapper(
#'   X = X, Y = Y,
#'   data_driven_prior_matrices = data_driven_prior_matrices,
#'   prior_grid = prior_grid
#' )
#'
#' @importFrom doFuture registerDoFuture
#' @importFrom future plan multicore
#' @export
mrmash_wrapper <- function(X,
                           Y,
                           V = NULL,
                           sumstats = NULL,
                           data_driven_prior_matrices = NULL,
                           prior_grid = NULL,
                           nthreads = 1,
                           canonical_prior_matrices = FALSE,
                           standardize = FALSE,
                           update_w0 = TRUE,
                           w0_threshold = 1e-8,
                           update_V = TRUE,
                           update_V_method = "full",
                           B_init_method = "enet",
                           max_iter = 5000,
                           tol = 0.01,
                           verbose = FALSE, ...) {
  # Make sure glmnet is installed
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("To use this function, please install glmnet: https://cran.r-project.org/web/packages/glmnet/index.html")
  }
  # Make sure mr.mashr is installed
  if (!requireNamespace("mr.mashr", quietly = TRUE)) {
    stop("To use this function, please install mr.mashr: https://github.com/stephenslab/mr.mashr")
  }
  # Check input data
  if (!exists(".Random.seed")) {
    message("! No seed has been set. Please set seed for reproducable result. ")
  }

  if (!is.matrix(X) || !is.matrix(Y)) {
    stop("X and Y must be matrices.")
  }

  if (nrow(X) != nrow(Y)) {
    stop("X and Y must have the same number of rows.")
  }
  if (!is.null(prior_grid) && !is.vector(prior_grid)) {
    stop("prior_grid must be a vector.")
  }
  if (is.null(data_driven_prior_matrices) && !isTRUE(canonical_prior_matrices)) {
    stop("Please provide data_driven_prior_matrices or set canonical_prior_matrices = TRUE.")
  }

  Y_has_missing <- any(is.na(Y))

  if (Y_has_missing && B_init_method == "glasso") {
    warning("B_init_method = 'glasso' can only be used without missing values in Y. Setting it to 'enet' instead")
    B_init_method <- "enet"
  }

  # Compute summary statistics and prior_grids
  if (is.null(sumstats)) {
    sumstats <- mr.mashr::compute_univariate_sumstats(X, Y,
      standardize = standardize,
      standardize.response = FALSE, mc.cores = nthreads
    )
  }

  prior_grid <- compute_grid(bhat = sumstats$Bhat, sbhat = sumstats$Shat)

  # Compute canonical matrices, if requested
  if (isTRUE(canonical_prior_matrices)) {
    canonical_prior_matrices <- mr.mashr::compute_canonical_covs(ncol(Y),
      singletons = TRUE,
      hetgrid = c(0, 0.25, 0.5, 0.75, 1)
    )
    if (!is.null(data_driven_prior_matrices)) {
      S0_raw <- c(canonical_prior_matrices, data_driven_prior_matrices$U)
    } else {
      S0_raw <- canonical_prior_matrices
    }
  } else {
    S0_raw <- data_driven_prior_matrices$U
  }

  # Compute prior covariance
  S0 <- mr.mashr::expand_covs(S0_raw, prior_grid, zeromat = TRUE)
  time1 <- proc.time()

  if (B_init_method == "glasso") {
    out <- compute_coefficients_glasso(X, Y,
      standardize = standardize,
      nthreads = nthreads, Xnew = NULL
    )
  } else {
    out <- compute_coefficients_univ_glmnet(X, Y,
      alpha = 0.5, standardize = standardize,
      nthreads = nthreads, Xnew = NULL
    )
  }

  B_init <- as.matrix(out$Bhat)
  w0 <- compute_w0(B_init, length(S0))

  # Robust initialization of V
  if (is.null(V)) {
    if (!Y_has_missing) {
      V <- mr.mashr:::compute_V_init(X, Y, matrix(0, nrow = ncol(X), ncol = ncol(Y)), rep(0, ncol(Y)), method = "cov")
    } else {
      muy <- colMeans(Y, na.rm = TRUE)
      V <- mr.mashr:::compute_V_init(X, Y, matrix(0, nrow = ncol(X), ncol = ncol(Y)), muy, method = "flash")
    }
    if (update_V_method == "diagonal") {
      V <- diag(diag(V))
    } else {
      if (any(eigen(V)$values < 1e-8)) {
        V <- V + diag(1e-8, nrow(V))
        update_V <- FALSE
      }
    }
  }

  # Fit mr.mash
  fit_mrmash <- mr.mashr::mr.mash(
    X = X, Y = Y, V = V, S0 = S0, w0 = w0, update_w0 = update_w0, tol = tol,
    max_iter = max_iter, convergence_criterion = "ELBO", compute_ELBO = TRUE,
    standardize = standardize, verbose = verbose, update_V = update_V,
    update_V_method = update_V_method, w0_threshold = w0_threshold,
    nthreads = nthreads, mu1_init = B_init
  )

  time2 <- proc.time()
  elapsed_time <- time2["elapsed"] - time1["elapsed"]
  fit_mrmash$analysis_time <- elapsed_time

  return(fit_mrmash)
}

### Function to compute initial estimates of the coefficients from group-lasso
compute_coefficients_glasso <- function(X, Y, standardize, nthreads, Xnew = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  r <- ncol(Y)
  condition_names <- colnames(Y)

  # Fit group-lasso
  cvfit_glmnet <- glmnet::cv.glmnet(
    x = X, y = Y, family = "mgaussian", alpha = 1,
    standardize = standardize, parallel = FALSE
  )
  coeff_glmnet <- coef(cvfit_glmnet, s = "lambda.min")

  # Build matrix of initial estimates for mr.mash
  B <- matrix(as.numeric(NA), nrow = p, ncol = r)

  for (i in seq_along(coeff_glmnet)) {
    B[, i] <- as.vector(coeff_glmnet[[i]])[-1]
  }

  # Make predictions if requested.
  if (!is.null(Xnew)) {
    Yhat_glmnet <- drop(predict(cvfit_glmnet, newx = Xnew, s = "lambda.min"))
    colnames(Yhat_glmnet) <- condition_names
    res <- list(Bhat = B, Ytrain = Y, Yhat_new = Yhat_glmnet)
  } else {
    res <- list(Bhat = B, Ytrain = Y)
  }
  return(res)
}

### Function to compute coefficients for univariate glmnet
compute_coefficients_univ_glmnet <- function(X, Y, alpha, standardize, nthreads, Xnew = NULL) {
  r <- ncol(Y)

  linreg <- function(i, X, Y, alpha, standardize, nthreads, Xnew) {
    samples_kept <- which(!is.na(Y[, i]))
    Ynomiss <- Y[samples_kept, i, drop = FALSE]
    Xnomiss <- X[samples_kept, , drop = FALSE]

    cvfit <- glmnet::cv.glmnet(
      x = Xnomiss, y = Ynomiss, family = "gaussian", alpha = alpha,
      standardize = standardize, parallel = FALSE
    )
    coeffic <- as.vector(coef(cvfit, s = "lambda.min"))
    lambda_seq <- cvfit$lambda

    # Make predictions if requested
    if (!is.null(Xnew)) {
      yhat_glmnet <- drop(predict(cvfit, newx = Xnew, s = "lambda.min"))
      res <- list(bhat = coeffic, lambda_seq = lambda_seq, yhat_new = yhat_glmnet)
    } else {
      res <- list(bhat = coeffic, lambda_seq = lambda_seq)
    }

    return(res)
  }

  out <- lapply(1:r, linreg, X, Y, alpha, standardize, nthreads, Xnew)

  Bhat <- sapply(out, "[[", "bhat")

  if (!is.null(Xnew)) {
    Yhat_new <- sapply(out, "[[", "yhat_new")
    colnames(Yhat_new) <- colnames(Y)
    results <- list(Bhat = Bhat[-1, ], intercept = Bhat[1, ], Yhat_new = Yhat_new)
  } else {
    results <- list(Bhat = Bhat[-1, ], intercept = Bhat[1, ])
  }
  return(results)
}

### Compute prior weights from coefficients estimates
compute_w0 <- function(Bhat, ncomps) {
  prop_nonzero <- sum(rowSums(abs(Bhat)) > 0) / nrow(Bhat)

  if (ncomps > 1) {
    w0 <- c((1 - prop_nonzero), rep(prop_nonzero / (ncomps - 1), (ncomps - 1)))
  } else {
    w0 <- 1
  }

  if (sum(w0 != 0) < 2) {
    w0 <- rep(1 / ncomps, ncomps)
  }

  return(w0)
}

#' Re-normalize mrmash weight w0 to have total weight sum to 1
#' @param w0 is the weight of mr.mash prior matrices that was generated from mr.mash() function.
rescale_cov_w0 <- function(w0) {
  # remove null component
  w0 <- w0[names(w0) != "null"]

  # split by prior group
  groups <- sub("_[^_]+$", "", names(w0))
  group_list <- split(w0, groups)

  # get per group sum
  group_weight <- lapply(group_list, sum)

  # Renormalize values within each group
  weights_list <- unlist(group_weight)
  sum_weights <- sum(weights_list)
  if (sum_weights > 0) {
    weights_list <- weights_list / sum_weights
  } else {
    # Use equal weights if all non null weights are zeros
    weights_list <- setNames(rep(1 / length(weights_list), length(weights_list)), names(weights_list))
  }
  # vector to store updated group w0
  updated_w0 <- rep(NA, length(unique(groups)))
  names(updated_w0) <- unique(groups)

  # replace with updated values
  updated_w0[names(weights_list)] <- weights_list
  return(updated_w0)
}


### Function to compute grids
compute_grid <- function(bhat, sbhat) {
  grid_mins <- c()
  grid_maxs <- c()

  include <- !(sbhat == 0 | !is.finite(sbhat) | is.na(sbhat) | is.na(bhat))
  gmax <- grid_max(bhat[include], sbhat[include])
  gmin <- grid_min(bhat[include], sbhat[include])
  grid_mins <- c(grid_mins, gmin)
  grid_maxs <- c(grid_maxs, gmax)

  gmin_tot <- min(grid_mins)
  gmax_tot <- max(grid_maxs)
  grid <- autoselect_mixsd(gmin_tot, gmax_tot, mult = sqrt(2))^2

  return(grid)
}

### Compute the minimum value for the grid
grid_min <- function(bhat, sbhat) {
  min(sbhat)
}

### Compute the maximum value for the grid
grid_max <- function(bhat, sbhat) {
  if (all(bhat^2 <= sbhat^2)) {
    8 * grid_min(bhat, sbhat) # the unusual case where we don't need much grid
  } else {
    2 * sqrt(max(bhat^2 - sbhat^2))
  }
}

### Function to compute the grid
autoselect_mixsd <- function(gmin, gmax, mult = 2) {
  if (mult == 0) {
    return(c(0, gmax / 2))
  } else {
    npoint <- ceiling(log2(gmax / gmin) / log2(mult))
    return(mult^((-npoint):0) * gmax)
  }
}
