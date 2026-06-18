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
#' sampleId <- paste0("P000", str_pad(1:400, 3, pad = "0"))
#' X <- matrix(sample(0:2, size = n * p, replace = TRUE, prob = c(0.65, 0.30, 0.05)), nrow = n)
#' rownames(X) <- sampleId
#' colnames(X) <- paste0("rs", sample(10000:100000, p))
#'
#' tissues <- c(
#'   "Adipose Tissue", "Muscle Tissue", "Brain Tissue", "Liver Tissue",
#'   "Kidney Tissue", "Heart Tissue", "Lung Tissue"
#' )
#' Y <- matrix(runif(n * r, -2, 2), nrow = n)
#' Y <- scale(Y)
#' colnames(Y) <- tissues
#' rownames(Y) <- sampleId
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
#' res <- mrmashWrapper(
#'   X = X, Y = Y,
#'   dataDrivenPriorMatrices = dataDrivenPriorMatrices,
#'   priorGrid = priorGrid
#' )
#'
#' @export
mrmashWrapper <- function(X,
                          Y,
                          V = NULL,
                          sumstats = NULL,
                          dataDrivenPriorMatrices = NULL,
                          priorGrid = NULL,
                          nthreads = 1,
                          canonicalPriorMatrices = FALSE,
                          standardize = FALSE,
                          updateW0 = TRUE,
                          w0Threshold = 1e-8,
                          updateV = TRUE,
                          updateVMethod = "full",
                          bInitMethod = "enet",
                          maxIter = 5000,
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
  if (!is.null(priorGrid) && !is.vector(priorGrid)) {
    stop("priorGrid must be a vector.")
  }
  if (is.null(dataDrivenPriorMatrices) && !isTRUE(canonicalPriorMatrices)) {
    stop("Please provide dataDrivenPriorMatrices or set canonicalPriorMatrices = TRUE.")
  }

  yHasMissing <- any(is.na(Y))

  if (yHasMissing && bInitMethod == "glasso") {
    warning("bInitMethod = 'glasso' can only be used without missing values in Y. Setting it to 'enet' instead")
    bInitMethod <- "enet"
  }

  # Compute summary statistics and prior_grids
  if (is.null(sumstats)) {
    sumstats <- mr.mashr::compute_univariate_sumstats(X, Y,
      standardize = standardize,
      standardize.response = FALSE, mc.cores = nthreads
    )
  }

  # Build prior covariance via shared helper (also used by mrmashRssWeights)
  priorBuilt <- buildMrmashPriorMatrices(
    Bhat = sumstats$Bhat, Shat = sumstats$Shat,
    K = ncol(Y),
    dataDrivenPriorMatrices = dataDrivenPriorMatrices,
    canonicalPriorMatrices = canonicalPriorMatrices,
    priorGrid = priorGrid
  )
  S0 <- priorBuilt$S0
  priorGrid <- priorBuilt$priorGrid
  time1 <- proc.time()

  if (bInitMethod == "glasso") {
    out <- computeCoefficientsGlasso(X, Y,
      standardize = standardize,
      nthreads = nthreads, Xnew = NULL
    )
  } else {
    out <- computeCoefficientsUnivGlmnet(X, Y,
      alpha = 0.5, standardize = standardize,
      nthreads = nthreads, Xnew = NULL
    )
  }

  B_init <- as.matrix(out$Bhat)
  w0 <- computeW0(B_init, length(S0))

  # Robust initialization of V
  if (is.null(V)) {
    if (!yHasMissing) {
      V <- mr.mashr:::compute_V_init(X, Y, matrix(0, nrow = ncol(X), ncol = ncol(Y)), rep(0, ncol(Y)), method = "cov")
    } else {
      muy <- colMeans(Y, na.rm = TRUE)
      V <- mr.mashr:::compute_V_init(X, Y, matrix(0, nrow = ncol(X), ncol = ncol(Y)), muy, method = "flash")
    }
    if (updateVMethod == "diagonal") {
      V <- diag(diag(V))
    } else {
      if (any(eigen(V)$values < 1e-8)) {
        V <- V + diag(1e-8, nrow(V))
        updateV <- FALSE
      }
    }
  }

  # Fit mr.mash
  fitMrmash <- mr.mashr::mr.mash(
    X = X, Y = Y, V = V, S0 = S0, w0 = w0, update_w0 = updateW0, tol = tol,
    max_iter = maxIter, convergence_criterion = "ELBO", compute_ELBO = TRUE,
    standardize = standardize, verbose = verbose, update_V = updateV,
    update_V_method = updateVMethod, w0_threshold = w0Threshold,
    nthreads = nthreads, mu1_init = B_init
  )

  time2 <- proc.time()
  elapsedTime <- time2["elapsed"] - time1["elapsed"]
  fitMrmash$analysis_time <- elapsedTime

  return(fitMrmash)
}

#' @export

### Function to compute initial estimates of the coefficients from group-lasso
computeCoefficientsGlasso <- function(X, Y, standardize, nthreads, Xnew = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  r <- ncol(Y)
  conditionNames <- colnames(Y)

  # Fit group-lasso
  cvfitGlmnet <- glmnet::cv.glmnet(
    x = X, y = Y, family = "mgaussian", alpha = 1,
    standardize = standardize, parallel = FALSE
  )
  coeffGlmnet <- coef(cvfitGlmnet, s = "lambda.min")

  # Build matrix of initial estimates for mr.mash
  B <- matrix(as.numeric(NA), nrow = p, ncol = r)

  for (i in seq_along(coeffGlmnet)) {
    B[, i] <- as.vector(coeffGlmnet[[i]])[-1]
  }

  # Make predictions if requested.
  if (!is.null(Xnew)) {
    YhatGlmnet <- drop(predict(cvfitGlmnet, newx = Xnew, s = "lambda.min"))
    colnames(YhatGlmnet) <- conditionNames
    res <- list(Bhat = B, Ytrain = Y, Yhat_new = YhatGlmnet)
  } else {
    res <- list(Bhat = B, Ytrain = Y)
  }
  return(res)
}


### Function to compute coefficients for univariate glmnet
computeCoefficientsUnivGlmnet <- function(X, Y, alpha, standardize, nthreads, Xnew = NULL) {
  r <- ncol(Y)

  linreg <- function(i, X, Y, alpha, standardize, nthreads, Xnew) {
    samplesKept <- which(!is.na(Y[, i]))
    Ynomiss <- Y[samplesKept, i, drop = FALSE]
    Xnomiss <- X[samplesKept, , drop = FALSE]

    cvfit <- glmnet::cv.glmnet(
      x = Xnomiss, y = Ynomiss, family = "gaussian", alpha = alpha,
      standardize = standardize, parallel = FALSE
    )
    coeffic <- as.vector(coef(cvfit, s = "lambda.min"))
    lambdaSeq <- cvfit$lambda

    # Make predictions if requested
    if (!is.null(Xnew)) {
      yhatGlmnet <- drop(predict(cvfit, newx = Xnew, s = "lambda.min"))
      res <- list(bhat = coeffic, lambda_seq = lambdaSeq, yhat_new = yhatGlmnet)
    } else {
      res <- list(bhat = coeffic, lambda_seq = lambdaSeq)
    }

    return(res)
  }

  out <- lapply(1:r, linreg, X, Y, alpha, standardize, nthreads, Xnew)

  Bhat <- sapply(out, "[[", "bhat")

  if (!is.null(Xnew)) {
    YhatNew <- sapply(out, "[[", "yhat_new")
    colnames(YhatNew) <- colnames(Y)
    results <- list(Bhat = Bhat[-1, ], intercept = Bhat[1, ], Yhat_new = YhatNew)
  } else {
    results <- list(Bhat = Bhat[-1, ], intercept = Bhat[1, ])
  }
  return(results)
}


### Compute prior weights from coefficients estimates
computeW0 <- function(Bhat, ncomps) {
  propNonzero <- sum(rowSums(abs(Bhat)) > 0) / nrow(Bhat)

  if (ncomps > 1) {
    w0 <- c((1 - propNonzero), rep(propNonzero / (ncomps - 1), (ncomps - 1)))
  } else {
    w0 <- 1
  }

  if (sum(w0 != 0) < 2) {
    w0 <- rep(1 / ncomps, ncomps)
  }

  return(w0)
}

compute_w0 <- computeW0

#' Re-normalize mrmash weight w0 to have total weight sum to 1
#' @param w0 is the weight of mr.mash prior matrices that was generated from mr.mash() function.
rescaleCovW0 <- function(w0) {
  # remove null component
  w0 <- w0[names(w0) != "null"]

  # split by prior group
  groups <- sub("_[^_]+$", "", names(w0))
  groupList <- split(w0, groups)

  # get per group sum
  groupWeight <- lapply(groupList, sum)

  # Renormalize values within each group
  weightsList <- unlist(groupWeight)
  sumWeights <- sum(weightsList)
  if (sumWeights > 0) {
    weightsList <- weightsList / sumWeights
  } else {
    # Use equal weights if all non null weights are zeros
    weightsList <- setNames(rep(1 / length(weightsList), length(weightsList)), names(weightsList))
  }
  # vector to store updated group w0
  updatedW0 <- rep(NA, length(unique(groups)))
  names(updatedW0) <- unique(groups)

  # replace with updated values
  updatedW0[names(weightsList)] <- weightsList
  return(updatedW0)
}

rescale_cov_w0 <- rescaleCovW0

### Function to compute grids
computeGrid <- function(bhat, sbhat) {
  gridMins <- c()
  gridMaxs <- c()

  include <- !(sbhat == 0 | !is.finite(sbhat) | is.na(sbhat) | is.na(bhat))
  gmax <- gridMax(bhat[include], sbhat[include])
  gmin <- gridMin(bhat[include], sbhat[include])
  gridMins <- c(gridMins, gmin)
  gridMaxs <- c(gridMaxs, gmax)

  gminTot <- min(gridMins)
  gmaxTot <- max(gridMaxs)
  grid <- autoselectMixsd(gminTot, gmaxTot, mult = sqrt(2))^2

  return(grid)
}


### Compute the minimum value for the grid
gridMin <- function(bhat, sbhat) {
  min(sbhat)
}


### Compute the maximum value for the grid
gridMax <- function(bhat, sbhat) {
  if (all(bhat^2 <= sbhat^2)) {
    8 * gridMin(bhat, sbhat) # the unusual case where we don't need much grid
  } else {
    2 * sqrt(max(bhat^2 - sbhat^2))
  }
}


### Function to compute the grid
autoselectMixsd <- function(gmin, gmax, mult = 2) {
  if (mult == 0) {
    return(c(0, gmax / 2))
  } else {
    npoint <- ceiling(log2(gmax / gmin) / log2(mult))
    return(mult^((-npoint):0) * gmax)
  }
}


#' Compute covariance matrix using FLASH
#'
#' Estimates a covariance matrix from a data matrix Y using empirical Bayes
#' matrix factorization (\code{flashier::flash}). When the FLASH fit finds
#' no shared factors, the returned covariance is diagonal with entries
#' \code{residuals_sd^2}; otherwise the factor contribution is added.
#' FLASH errors are not caught; callers should handle them explicitly or
#' supply a pre-computed prior covariance instead.
#'
#' @param Y Numeric matrix (samples x conditions).
#' @return A covariance matrix of dimension ncol(Y) x ncol(Y), rescaled by
#'   the column standard deviations of Y.
#' @export
computeCovFlash <- function(Y) {
  # flashier >= 1.0 API: var_type / ebnm_fn / verbose (renamed from var.type
  # / prior.family / verbose.lvl). Prior families now come from `ebnm`.
  fl <- flashier::flash(Y, var_type = 2,
    ebnm_fn = c(ebnm::ebnm_normal, ebnm::ebnm_normal_scale_mixture),
    backfit = TRUE, verbose = 0)
  if (fl$n_factors == 0) {
    covar <- diag(fl$residuals_sd^2)
  } else {
    # For each factor's right-side prior, marginal variance for a
    # mean-zero scale-mixture-of-normals is sum(pi * sd^2).
    fsd <- vapply(fl$F_ghat, function(g) sqrt(sum(g$pi * g$sd^2)), numeric(1))
    covar <- diag(fl$residuals_sd^2) + crossprod(t(fl$F_pm) * fsd)
  }
  if (nrow(covar) == 0) {
    stop("computeCovFlash: FLASH produced an empty covariance matrix.")
  }
  s <- apply(Y, 2, sd, na.rm = TRUE)
  if (length(s) > 1) s <- diag(s) else s <- matrix(s, 1, 1)
  s %*% cov2cor(covar) %*% s
}

#' Compute diagonal covariance matrix
#'
#' Returns a diagonal covariance matrix from the column-wise variances of Y.
#'
#' @param Y Numeric matrix (samples x conditions).
#' @return A diagonal covariance matrix of dimension ncol(Y) x ncol(Y).
#' @export
computeCovDiag <- function(Y) {
  diag(apply(Y, 2, var, na.rm = TRUE))
}

#' Build mr.mash prior covariance matrices
#'
#' Shared helper used by both \code{\link{mrmashWrapper}} (individual-level)
#' and \code{\link{mrmashRssWeights}} (summary statistics). Constructs the
#' \code{S0} list of prior covariance matrices via the canonical mixture
#' (\code{mr.mashr::compute_canonical_covs}) and optional data-driven
#' matrices, expanded over a scaling grid via
#' \code{mr.mashr::expand_covs}. The prior grid is derived from \code{Bhat}
#' and \code{Shat} via \code{\link{computeGrid}} when not supplied.
#'
#' @param Bhat Numeric matrix of effect-size estimates (variants x conditions).
#' @param Shat Numeric matrix of standard errors (variants x conditions).
#' @param K Number of conditions. When NULL, inferred from \code{ncol(Bhat)}.
#' @param dataDrivenPriorMatrices Optional list with element \code{U}
#'   (list of raw covariance matrices) computed e.g. by
#'   \code{\link{computeCovFlash}} / \code{\link{computeCovDiag}}.
#' @param canonicalPriorMatrices Logical. When TRUE (default for RSS),
#'   include the standard canonical mixture from
#'   \code{mr.mashr::compute_canonical_covs()}. When FALSE,
#'   \code{dataDrivenPriorMatrices} must be supplied.
#' @param priorGrid Optional pre-computed scaling grid (numeric vector).
#'   When NULL, derived from \code{Bhat}, \code{Shat} via
#'   \code{computeGrid()}.
#' @param hetgrid Heterogeneity grid passed to
#'   \code{mr.mashr::compute_canonical_covs()}. Default
#'   \code{c(0, 0.25, 0.5, 0.75, 1)}, matching the individual-level wrapper.
#' @param singletons Whether to include single-condition prior components.
#'   Default TRUE.
#' @return A list with components \code{S0} (the expanded list of prior
#'   covariance matrices) and \code{prior_grid} (the scaling grid that was
#'   used).
#' @export
buildMrmashPriorMatrices <- function(Bhat, Shat, K = NULL,
                                     dataDrivenPriorMatrices = NULL,
                                     canonicalPriorMatrices = TRUE,
                                     priorGrid = NULL,
                                     hetgrid = c(0, 0.25, 0.5, 0.75, 1),
                                     singletons = TRUE) {
  if (!requireNamespace("mr.mashr", quietly = TRUE)) {
    stop("Package 'mr.mashr' is required.")
  }
  if (is.null(dataDrivenPriorMatrices) && !isTRUE(canonicalPriorMatrices)) {
    stop("Supply dataDrivenPriorMatrices or set canonicalPriorMatrices = TRUE.")
  }
  if (is.null(K)) K <- ncol(Bhat)
  if (is.null(priorGrid)) priorGrid <- computeGrid(bhat = Bhat, sbhat = Shat)

  if (isTRUE(canonicalPriorMatrices)) {
    canonical <- mr.mashr::compute_canonical_covs(K, singletons = singletons, hetgrid = hetgrid)
    S0_raw <- if (!is.null(dataDrivenPriorMatrices)) {
      c(canonical, dataDrivenPriorMatrices$U)
    } else {
      canonical
    }
  } else {
    S0_raw <- dataDrivenPriorMatrices$U
  }

  S0 <- mr.mashr::expand_covs(S0_raw, priorGrid, zeromat = TRUE)
  list(S0 = S0, priorGrid = priorGrid)
}

