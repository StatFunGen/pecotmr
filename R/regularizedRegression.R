#' Extract weights from mr.ash.rss (susieR)
#' @return A numeric vector of the posterior mean of the coefficients.
#' @importFrom susieR mr.ash.rss
#' @export
mrAshRssWeights <- function(stat, LD, varY, sigma2E, s0, w0, z = numeric(0), ...) {
  model <- mr.ash.rss(
    bhat = stat$b, shat = stat$seb, z = z, R = LD,
    var_y = varY, n = median(stat$n), sigma2_e = sigma2E,
    s0 = s0, w0 = w0, ...
  )

  return(model$mu1)
}

#' PRS-CS: a polygenic prediction method that infers posterior SNP effect sizes under continuous shrinkage (CS) priors
#'
#' This function is a wrapper for the PRS-CS method implemented in C++. It takes marginal effect size estimates from regression and an external LD reference panel
#' and infers posterior SNP effect sizes using Bayesian regression with continuous shrinkage priors.
#'
#' @param bhat A vector of marginal effect sizes.
#' @param LD A list of LD blocks, where each element is a matrix representing an LD block.
#' @param n Sample size of the GWAS.
#' @param a Shape parameter for the prior distribution of psi. Default is 1.
#' @param b Scale parameter for the prior distribution of psi. Default is 0.5.
#' @param phi Global shrinkage parameter. If NULL, it will be estimated automatically. Default is NULL.
#' @param nIter Number of MCMC iterations. Default is 1000.
#' @param nBurnin Number of burn-in iterations. Default is 500.
#' @param thin Thinning factor for MCMC. Default is 5.
#' @param maf A vector of minor allele frequencies, if available, will standardize the effect sizes by MAF. Default is NULL.
#' @param verbose Whether to print verbose output. Default is FALSE.
#' @param seed Random seed for reproducibility. Default is NULL.
#'
#' @return A list containing the posterior estimates:
#'   - betaEst: Posterior estimates of SNP effect sizes.
#'   - psiEst: Posterior estimates of psi (shrinkage parameters).
#'   - sigmaEst: Posterior estimate of the residual variance.
#'   - phiEst: Posterior estimate of the global shrinkage parameter.
#' @examples
#' # Generate example data
#' set.seed(985115)
#' n <- 350
#' p <- 16
#' sigmasq_error <- 0.5
#' zeroes <- rbinom(p, 1, 0.6)
#' beta.true <- rnorm(p, 1, sd = 4)
#' beta.true[zeroes] <- 0
#'
#' X <- cbind(matrix(rnorm(n * p), nrow = n))
#' X <- scale(X, center = TRUE, scale = FALSE)
#' y <- X %*% matrix(beta.true, ncol = 1) + rnorm(n, 0, sqrt(sigmasq_error))
#' y <- scale(y, center = TRUE, scale = FALSE)
#'
#' # Calculate sufficient statistics
#' XtX <- t(X) %*% X
#' Xty <- t(X) %*% y
#' yty <- t(y) %*% y
#'
#' # Set the prior
#' K <- 9
#' sigma0 <- c(0.001, .1, .5, 1, 5, 10, 20, 30, .005)
#' omega0 <- rep(1 / K, K)
#'
#' # Calculate summary statistics
#' b.hat <- sapply(1:p, function(j) {
#'   summary(lm(y ~ X[, j]))$coefficients[-1, 1]
#' })
#' s.hat <- sapply(1:p, function(j) {
#'   summary(lm(y ~ X[, j]))$coefficients[-1, 2]
#' })
#' R.hat <- cor(X)
#' var_y <- var(y)
#' sigmasq_init <- 1.5
#'
#' # Run PRS CS
#' maf <- rep(0.5, length(b.hat)) # fake MAF
#' LD <- list(blk1 = R.hat)
#' out <- prsCs(b.hat, LD, n, maf = maf)
#' # In sample prediction correlations
#' cor(X %*% out$betaEst, y) # 0.9944553
#' @export
prsCs <- function(bhat, LD, n,
                  a = 1, b = 0.5, phi = NULL,
                  maf = NULL, nIter = 1000, nBurnin = 500,
                  thin = 5, verbose = FALSE, seed = NULL) {
  # Check input parameters
  if (missing(LD) || !is.list(LD)) {
    stop("Please provide a valid list of LD blocks using 'LD'.")
  }
  if (missing(n) || n <= 0) {
    stop("Please provide a valid sample size using 'n'.")
  }

  # Check if maf is provided and its length matches that of bhat
  if (!is.null(maf) && length(bhat) != length(maf)) {
    stop("The length of 'bhat' must be the same as 'maf'.")
  }

  # Check if the length of bhat matches the sum of the nrow of all elements in the LD list
  totalRowsInLd <- sum(sapply(LD, nrow))
  if (length(bhat) != totalRowsInLd) {
    stop("The length of 'bhat' must be the same as the sum of the number of rows of all elements in the 'LD' list.")
  }

  # Run PRS-CS
  # cpp11 requires exact integer types for int parameters
  result <- prsCsRcpp(
    a = a, b = b, phi = phi, bhat = bhat, maf = maf,
    n = as.integer(n), ldBlk = LD,
    nIter = as.integer(nIter), nBurnin = as.integer(nBurnin), thin = as.integer(thin),
    verbose = verbose, seed = seed
  )

  # Return the result as a list (camelCase to match the rest of the package API).
  list(
    betaEst = result$betaEst,
    psiEst = result$psiEst,
    sigmaEst = result$sigmaEst,
    phiEst = result$phiEst
  )
}

#' Extract weights from prsCs function
#' @return A numeric vector of the posterior SNP coefficients.
#' @export
prsCsWeights <- function(stat, LD, ...) {
  model <- prsCs(bhat = stat$b, LD = list(blk1 = LD), n = median(stat$n), ...)

  return(model$betaEst)
}

#' SDPR (Summary-Statistics-Based Dirichelt Process Regression for Polygenic Risk Prediction)
#'
#' This function is a wrapper for the SDPR C++ implementation, which performs Markov Chain Monte Carlo (MCMC)
#' for estimating effect sizes and heritability based on summary statistics and reference LD matrices.
#'
#' @param bhat A vector of marginal beta values for each SNP.
#' @param LD A list of LD matrices, where each matrix corresponds to a subset of SNPs.
#' @param n The total sample size of the GWAS.
#' @param perVariantSampleSize (Optional) A vector of sample sizes for each SNP. If NULL (default), it will be initialized
#'                    to a vector of length equal to `bhat`, with all values set to `n`.
#' @param array (Optional) A vector of genotyping array information for each SNP. If NULL (default), it will be
#'              initialized to a vector of 1's with length equal to `bhat`.
#' @param a Factor to shrink the reference LD matrix. Default is 0.1.
#' @param c Factor to correct for the deflation. Default is 1.
#' @param M Max number of variance components. Default is 1000.
#' @param a0k Hyperparameter for inverse gamma distribution. Default is 0.5.
#' @param b0k Hyperparameter for inverse gamma distribution. Default is 0.5.
#' @param iter Number of iterations for MCMC. Default is 1000.
#' @param burn Number of burn-in iterations for MCMC. Default is 200.
#' @param thin Thinning interval for MCMC. Default is 5.
#' @param nThreads Number of threads to use. Default is 1.
#' @param optLlk Which likelihood to evaluate. 1 for equation 6 (slightly shrink the correlation of SNPs)
#'                and 2 for equation 5 (SNPs genotyped on different arrays in a separate cohort).
#'                Default is 1.
#' @param verbose Whether to print verbose output. Default is true.
#'
#' @return A list containing the estimated effect sizes (beta) and heritability (h2).
#' @examples
#' # Generate example data
#' set.seed(985115)
#' n <- 350
#' p <- 16
#' sigmasq_error <- 0.5
#' zeroes <- rbinom(p, 1, 0.6)
#' beta.true <- rnorm(p, 1, sd = 4)
#' beta.true[zeroes] <- 0
#'
#' X <- cbind(matrix(rnorm(n * p), nrow = n))
#' X <- scale(X, center = TRUE, scale = FALSE)
#' y <- X %*% matrix(beta.true, ncol = 1) + rnorm(n, 0, sqrt(sigmasq_error))
#' y <- scale(y, center = TRUE, scale = FALSE)
#'
#' # Calculate sufficient statistics
#' XtX <- t(X) %*% X
#' Xty <- t(X) %*% y
#' yty <- t(y) %*% y
#'
#' # Set the prior
#' K <- 9
#' sigma0 <- c(0.001, .1, .5, 1, 5, 10, 20, 30, .005)
#' omega0 <- rep(1 / K, K)
#'
#' # Calculate summary statistics
#' b.hat <- sapply(1:p, function(j) {
#'   summary(lm(y ~ X[, j]))$coefficients[-1, 1]
#' })
#' s.hat <- sapply(1:p, function(j) {
#'   summary(lm(y ~ X[, j]))$coefficients[-1, 2]
#' })
#' R.hat <- cor(X)
#' var_y <- var(y)
#' sigmasq_init <- 1.5
#'
#' # Run SDPR
#' LD <- list(blk1 = R.hat)
#' out <- sdpr(b.hat, LD, n)
#' # In sample prediction correlations
#' cor(X %*% out$betaEst, y) #
#'
#' @note This function is a wrapper for the SDPR C++ implementation, which is a rewritten and adopted version
#'       of the SDPR package. The original SDPR documentation is available at
#'       https://htmlpreview.github.io/?https://github.com/eldronzhou/SDPR/blob/main/doc/Manual.html
#'
#' @export
sdpr <- function(bhat, LD, n, perVariantSampleSize = NULL, array = NULL, a = 0.1, c = 1.0, M = 1000,
                 a0k = 0.5, b0k = 0.5, iter = 1000, burn = 200, thin = 5, nThreads = 1,
                 optLlk = 1, verbose = TRUE, seed = NULL) {
  # Check if the sum of the rows in LD list is the same as length of bhat
  if (sum(sapply(LD, nrow)) != length(bhat)) {
    stop("The sum of the rows in LD list must be the same as the length of bhat.")
  }

  # Check if total sample size n is a positive integer
  if (missing(n) || n <= 0) {
    stop("The total sample size 'n' must be a positive integer.")
  }

  # M must be >= 4 (SDPR uses M-2 indexing in sample_V; M < 4 causes buffer overflow)
  if (M < 4) {
    stop("'M' must be at least 4.")
  }

  # Check if perVariantSampleSize vector contains only positive values (if provided)
  if (!is.null(perVariantSampleSize) && any(perVariantSampleSize <= 0)) {
    stop("The 'perVariantSampleSize' vector must contain only positive values.")
  }

  # Check if array vector contains only 0, 1, or 2 (if provided)
  if (!is.null(array) && any(!array %in% c(0, 1, 2))) {
    stop("The 'array' vector must contain only 0, 1, or 2.")
  }

  # cpp11 requires exact integer types for int parameters and sexp-wrapped vectors
  if (!is.null(array)) array <- as.integer(array)
  # Call the sdprRcpp function
  result <- sdprRcpp(
    bhatR = bhat, LD = LD, n = as.integer(n),
    perVariantSampleSize = perVariantSampleSize, array = array,
    a = a, c = c, M = as.integer(M),
    a0k = a0k, b0k = b0k,
    iter = as.integer(iter), burn = as.integer(burn), thin = as.integer(thin),
    nThreads = as.integer(nThreads), optLlk = as.integer(optLlk),
    verbose = verbose, seed = seed
  )

  return(result)
}

#' Extract weights from sdpr function
#' @return A numeric vector of the posterior SNP coefficients.
#' @export
sdprWeights <- function(stat, LD, ...) {
  model <- sdpr(bhat = stat$b, LD = list(blk1 = LD), n = median(stat$n), ...)

  return(model$betaEst)
}

# Shared helper for susie/susieAsh/susieInf weight extraction.
# @param fit A susie fit object (or NULL to fit from X, y).
# @param X Genotype matrix (optional).
# @param y Phenotype vector (optional).
# @param requiredFields Fields that must be present in the fit to extract weights.
# @param fitArgs Extra arguments passed to susieR::susie when fit is NULL.
# @param ... Additional arguments forwarded to susieR::susie.
#' @importFrom susieR coef.susie susie
#' @noRd
.susieExtractWeights <- function(fit, X, y, requiredFields, fitArgs = list(), retainFit = FALSE, ...) {
  if (is.null(fit)) {
    fit <- do.call(susie, c(list(X = X, y = y), fitArgs, list(...)))
  }
  if (!is.null(X) && length(fit$pip) != ncol(X)) {
    stop(paste0(
      "Dimension mismatch on number of variant in susie fit ", length(fit$pip),
      " and TWAS weights ", ncol(X), ". "
    ))
  }
  if (all(requiredFields %in% names(fit))) {
    fit$intercept <- 0
    weights <- coef.susie(fit)[-1]
  } else {
    weights <- rep(0, length(fit$pip))
  }
  if (retainFit) attr(weights, "fit") <- fit
  return(weights)
}

#' Compute SuSiE TWAS weights
#'
#' Extracts coefficients from an existing SuSiE fit or fits `susieR::susie()`
#' from `X` and `y` before extracting weights.
#'
#' @param X Genotype matrix. Required when `susieFit` is NULL.
#' @param y Phenotype vector. Required when `susieFit` is NULL.
#' @param susieFit Optional fitted SuSiE object.
#' @param retainFit If TRUE, stores the fitted object as an attribute on the returned weights.
#' @param ... Additional arguments passed to `susieR::susie()` when fitting.
#' @return Numeric vector of variant weights.
#' @export
susieWeights <- function(X = NULL, y = NULL, susieFit = NULL, retainFit = FALSE, ...) {
  .susieExtractWeights(susieFit, X, y,
    requiredFields = c("alpha", "mu", "X_column_scale_factors"),
    retainFit = retainFit, ...)
}

#' Compute SuSiE-ASH TWAS weights
#'
#' Extracts coefficients from an existing SuSiE-ASH fit or fits `susieR::susie()`
#' with `unmappable_effects = "ash"`.
#'
#' @param X Genotype matrix. Required when `susieAshFit` is NULL.
#' @param y Phenotype vector. Required when `susieAshFit` is NULL.
#' @param susieAshFit Optional fitted SuSiE-ASH object.
#' @param retainFit If TRUE, stores the fitted object as an attribute on the returned weights.
#' @param ... Additional arguments passed to `susieR::susie()` when fitting.
#' @return Numeric vector of variant weights.
#' @export
susieAshWeights <- function(X = NULL, y = NULL, susieAshFit = NULL, retainFit = FALSE, ...) {
  .susieExtractWeights(susieAshFit, X, y,
    requiredFields = c("alpha", "mu", "theta", "X_column_scale_factors"),
    fitArgs = list(unmappable_effects = "ash", convergence_method = "pip"),
    retainFit = retainFit, ...)
}

#' Compute SuSiE-inf TWAS weights
#'
#' Extracts coefficients from an existing SuSiE-inf fit or fits `susieR::susie()`
#' with `unmappable_effects = "inf"`.
#'
#' @section Non-zero weights with zero PIPs:
#' SuSiE-inf decomposes effects into a mappable component (driven by `alpha *
#' mu`, reported as per-variant PIPs) and an infinitesimal component (driven by
#' `theta`). When the fit converges with no mappable effects -- all `V` and `mu`
#' zero, so every `pip == 0` -- the returned weights are still non-zero because
#' `susieR::coef.susie` adds `theta / X_column_scale_factors` to the mappable
#' coefficient. This is intentional: it captures diffuse polygenic signal that
#' the mappable component could not localize to any credible set. Consumers
#' that interpret per-variant PIPs as a gate on whether to use the weights
#' should be aware that low or zero PIPs do not imply zero TWAS weights here.
#'
#' @param X Genotype matrix. Required when `susieInfFit` is NULL.
#' @param y Phenotype vector. Required when `susieInfFit` is NULL.
#' @param susieInfFit Optional fitted SuSiE-inf object.
#' @param retainFit If TRUE, stores the fitted object as an attribute on the returned weights.
#' @param ... Additional arguments passed to `susieR::susie()` when fitting.
#' @return Numeric vector of variant weights.
#' @export
susieInfWeights <- function(X = NULL, y = NULL, susieInfFit = NULL, retainFit = FALSE, ...) {
  .susieExtractWeights(susieInfFit, X, y,
    requiredFields = c("alpha", "mu", "theta", "X_column_scale_factors"),
    fitArgs = list(unmappable_effects = "inf", convergence_method = "pip"),
    retainFit = retainFit, ...)
}

# =============================================================================
# SuSiE-RSS weight functions
# =============================================================================

# Internal helper: extract weights from a susieRss fit.
# Mirrors .susie_extract_weights but uses the RSS interface.
#' @importFrom susieR coef.susie susieRss
#' @noRd
.susieRssExtractWeights <- function(fit, z, R, n,
                                    requiredFields, fitArgs = list(),
                                    retainFit = FALSE) {
  if (is.null(fit)) {
    fit <- do.call(susie_rss, c(list(z = z, R = R, n = n), fitArgs))
  }
  if (length(fit$pip) != nrow(R)) {
    stop(paste0(
      "Dimension mismatch: susieRss fit has ", length(fit$pip),
      " variants but R has ", nrow(R), " rows."))
  }
  if (all(requiredFields %in% names(fit))) {
    fit$intercept <- 0
    weights <- coef.susie(fit)[-1]
  } else {
    weights <- rep(0, length(fit$pip))
  }
  if (retainFit) attr(weights, "fit") <- fit
  return(weights)
}

#' Compute SuSiE-RSS TWAS weights
#'
#' Extracts coefficients from an existing SuSiE-RSS fit or fits
#' \code{susieR::susie_rss()} from summary statistics and LD.
#'
#' @param stat List with components \code{z} (z-scores), \code{n} (sample sizes).
#' @param LD LD correlation matrix.
#' @param susieRssFit Optional pre-fitted SuSiE-RSS object.
#' @param retainFit If TRUE, stores the fitted object as an attribute.
#' @param methodArgs Named list of additional arguments passed to
#'   \code{susieR::susie_rss()}. Use this instead of \code{...} to avoid
#'   partial matching of short argument names (e.g. \code{L}) to the
#'   \code{LD} parameter.
#' @return Numeric vector of variant weights.
#' @export
susieRssWeights <- function(stat, LD, susieRssFit = NULL, retainFit = TRUE,
                            methodArgs = list()) {
  .susieRssExtractWeights(fit = susieRssFit, z = stat$z, R = LD, n = median(stat$n),
    requiredFields = c("alpha", "mu", "X_column_scale_factors"),
    fitArgs = methodArgs,
    retainFit = retainFit)
}

#' Compute SuSiE-inf-RSS TWAS weights
#'
#' Extracts coefficients from an existing SuSiE-inf-RSS fit or fits
#' \code{susieR::susie_rss()} with \code{unmappable_effects = "inf"}.
#'
#' @inheritParams susieRssWeights
#' @param susieInfRssFit Optional pre-fitted SuSiE-inf-RSS object.
#' @return Numeric vector of variant weights.
#' @export
susieInfRssWeights <- function(stat, LD, susieInfRssFit = NULL, retainFit = TRUE,
                               methodArgs = list()) {
  .susieRssExtractWeights(fit = susieInfRssFit, z = stat$z, R = LD, n = median(stat$n),
    requiredFields = c("alpha", "mu", "theta", "X_column_scale_factors"),
    fitArgs = c(list(unmappable_effects = "inf", convergence_method = "pip"), methodArgs),
    retainFit = retainFit)
}

#' Compute SuSiE-ASH-RSS TWAS weights
#'
#' Extracts coefficients from an existing SuSiE-ASH-RSS fit or fits
#' \code{susieR::susie_rss()} with \code{unmappable_effects = "ash"}.
#'
#' @inheritParams susieRssWeights
#' @param susieAshRssFit Optional pre-fitted SuSiE-ASH-RSS object.
#' @return Numeric vector of variant weights.
#' @export
susieAshRssWeights <- function(stat, LD, susieAshRssFit = NULL, retainFit = TRUE,
                               methodArgs = list()) {
  .susieRssExtractWeights(fit = susieAshRssFit, z = stat$z, R = LD, n = median(stat$n),
    requiredFields = c("alpha", "mu", "theta", "X_column_scale_factors"),
    fitArgs = c(list(unmappable_effects = "ash", convergence_method = "pip"), methodArgs),
    retainFit = retainFit)
}

#' Compute mr.mash TWAS weights
#'
#' Extracts coefficients from an existing mr.mash fit or fits mr.mash from `X` and `Y`.
#'
#' @param mrmashFit Optional fitted mr.mash object.
#' @param X Genotype matrix. Required when `mrmashFit` is NULL.
#' @param Y Phenotype matrix. Required when `mrmashFit` is NULL.
#' @param ... Additional arguments passed to `mrmashWrapper()` when fitting.
#' @return Matrix of variant weights.
#' @export
mrmashWeights <- function(mrmashFit = NULL, X = NULL, Y = NULL, ...) {
  if (!requireNamespace("mr.mashr", quietly = TRUE)) {
    stop("Package 'mr.mashr' is required. Install with: devtools::install_github('stephenslab/mr.mashr')")
  }
  if (is.null(mrmashFit)) {
    message("mrmashFit is not provided; fitting mr.mash now ...")
    if (is.null(X) || is.null(Y)) {
      stop("Both X and Y must be provided if mrmashFit is NULL.")
    }
    mrmashFit <- mrmashWrapper(X, Y, ...)
  }
  return(mr.mashr::coef.mr.mash(mrmashFit)[-1, ])
}

#' Compute mvSuSiE TWAS weights
#'
#' Extracts coefficients from an existing mvSuSiE fit or fits `mvsusieR::mvsusie()`
#' from `X` and `Y`.
#'
#' @param mvsusieFit Optional fitted mvSuSiE object.
#' @param X Genotype matrix. Required when `mvsusieFit` is NULL.
#' @param Y Phenotype matrix. Required when `mvsusieFit` is NULL.
#' @param priorVariance Optional mvSuSiE prior variance list.
#' @param residualVariance Optional residual variance matrix.
#' @param L Maximum number of components.
#' @param LGreedy Initial greedy number of components.
#' @param verbose If TRUE, prints mvSuSiE fitting progress.
#' @param ... Additional arguments passed to `mvsusieR::mvsusie()` when fitting.
#' @return Matrix of variant weights.
#' @export
mvsusieWeights <- function(mvsusieFit = NULL, X = NULL, Y = NULL,
                           priorVariance = NULL, residualVariance = NULL,
                           L = 30, LGreedy = 5, verbose = FALSE, ...) {
  if (!requireNamespace("mvsusieR", quietly = TRUE)) {
    stop("Package 'mvsusieR' is required. Install with: devtools::install_github('stephenslab/mvsusieR')")
  }
  if (is.null(mvsusieFit)) {
    message("mvsusieFit is not provided; fitting mvSuSiE now ...")
    if (is.null(X) || is.null(Y)) {
      stop("Both X and Y must be provided if mvsusieFit is NULL.")
    }
    if (is.null(priorVariance)) priorVariance <- mvsusieR::create_mixture_prior(R = ncol(Y))
    if (!is.null(LGreedy)) LGreedy <- min(LGreedy, L)

    mvsusieFit <- mvsusieR::mvsusie(
      X = X, Y = Y, L = L, L_greedy = LGreedy, prior_variance = priorVariance,
      residual_variance = residualVariance,
      estimate_residual_variance = TRUE,
      verbose = verbose, ...
    )
  }
  return(mvsusieR::coef.mvsusie(mvsusieFit)[-1, ])
}

#' Compute mr.mash-RSS TWAS weights from summary statistics
#'
#' Multi-context summary-statistics analog of \code{\link{mrmashWeights}}:
#' extracts coefficients from an existing \code{mr.mashr::mr.mash.rss} fit,
#' or fits one from \code{stat} (variants x conditions) and \code{LD}.
#'
#' Follows the \code{*_rss_weights(stat, LD, ...)} contract. Expects
#' \code{stat$z} to be a numeric matrix (variants x conditions) and
#' \code{stat$n} a per-context numeric vector or scalar. \code{stat$Bhat}
#' and \code{stat$Shat} are used if present; otherwise derived from Z and
#' n.
#'
#' Prior construction reuses the same infrastructure as the individual-level
#' \code{\link{mrmashWrapper}}: \code{\link{computeGrid}} +
#' \code{mr.mashr::compute_canonical_covs()} +
#' \code{mr.mashr::expand_covs()} for \code{S0}, and
#' \code{\link{computeW0}} for the mixture weights. Supply
#' \code{dataDrivenPriorMatrices} (e.g. from
#' \code{\link{computeCovFlash}} / \code{\link{computeCovDiag}}) to add
#' data-driven covariance components alongside the canonical mixture.
#'
#' @param stat A list with \code{z} (variants x conditions matrix) and
#'   \code{n} (per-context numeric vector or scalar). May also include
#'   \code{Bhat}, \code{Shat} matrices.
#' @param LD LD correlation matrix.
#' @param mrmashRssFit Optional pre-fitted \code{mr.mash.rss} object;
#'   skips fitting when supplied.
#' @param dataDrivenPriorMatrices Optional list with element \code{U}
#'   (list of raw covariance matrices). Passed directly to
#'   \code{mr.mashr::expand_covs()} alongside the canonical mixture.
#' @param canonicalPriorMatrices Logical. When TRUE (default), include
#'   the standard canonical mixture from
#'   \code{mr.mashr::compute_canonical_covs()}. When FALSE,
#'   \code{dataDrivenPriorMatrices} must be supplied.
#' @param S0 Optional pre-built list of prior covariance matrices,
#'   bypassing the canonical / data-driven construction.
#' @param w0 Optional prior mixture weights; defaults to
#'   \code{\link{computeW0}(Bhat, length(S0))}.
#' @param V Optional residual covariance matrix (K x K). When NULL,
#'   defaults to the identity matrix of size K.
#' @param covY Optional response covariance matrix (K x K). When NULL,
#'   defaults to the identity matrix of size K.
#' @param retainFit If TRUE, attaches the fitted object as the
#'   \code{"fit"} attribute on the returned weights.
#' @param ... Additional arguments forwarded to
#'   \code{mr.mashr::mr.mash.rss}.
#'
#' @return A numeric matrix of per-variant per-context weights
#'   (variants x conditions).
#' @export
mrmashRssWeights <- function(stat, LD, mrmashRssFit = NULL,
                             dataDrivenPriorMatrices = NULL,
                             canonicalPriorMatrices = TRUE,
                             S0 = NULL, w0 = NULL, V = NULL, covY = NULL,
                             retainFit = FALSE, ...) {
  if (!requireNamespace("mr.mashr", quietly = TRUE)) {
    stop("Package 'mr.mashr' is required. ",
         "Install with: devtools::install_github('stephenslab/mr.mash.alpha')")
  }
  if (is.null(mrmashRssFit)) {
    Z <- if (is.matrix(stat$z)) stat$z else as.matrix(stat$z)
    if (ncol(Z) < 2) {
      stop("mrmashRssWeights expects stat$z to have >= 2 columns ",
           "(one per context). For single-context use mrAshRssWeights().")
    }
    K <- ncol(Z)
    nVec <- if (length(stat$n) > 1) stat$n else rep(stat$n, K)
    Bhat <- if (!is.null(stat$Bhat)) stat$Bhat else sweep(Z, 2, sqrt(nVec), "/")
    Shat <- if (!is.null(stat$Shat)) stat$Shat else matrix(1 / sqrt(rep(nVec, each = nrow(Z))),
                                                            nrow = nrow(Z), ncol = K)
    # Reuse the same prior-building helper as mrmashWrapper()
    if (is.null(S0)) {
      priorBuilt <- buildMrmashPriorMatrices(
        Bhat = Bhat, Shat = Shat, K = K,
        dataDrivenPriorMatrices = dataDrivenPriorMatrices,
        canonicalPriorMatrices = canonicalPriorMatrices
      )
      S0 <- priorBuilt$S0
    }
    if (is.null(w0)) {
      w0 <- computeW0(Bhat, length(S0))
    }
    if (is.null(V))    V    <- diag(K)
    if (is.null(covY)) covY <- diag(K)
    # mr.mash.rss expects either Z or (Bhat, Shat) but not both; prefer Bhat/Shat.
    # n must be a scalar (per the mr.mash.rss contract); use the median.
    nScalar <- as.numeric(stats::median(nVec))
    mrmashRssFit <- mr.mashr::mr.mash.rss(
      Bhat = Bhat, Shat = Shat, R = LD, n = nScalar,
      covY = covY, V = V, S0 = S0, w0 = w0, ...
    )
  }
  # coef.mr.mash.rss returns nrow(Bhat) rows (no intercept). Do not strip.
  weights <- mr.mashr::coef.mr.mash.rss(mrmashRssFit)
  if (retainFit) attr(weights, "fit") <- mrmashRssFit
  weights
}

#' Compute mvSuSiE-RSS TWAS weights from summary statistics
#'
#' Multi-context summary-statistics analog of \code{\link{mvsusieWeights}}:
#' extracts coefficients from an existing \code{mvsusieR::mvsusieRss} fit,
#' or fits one from \code{stat$z} (variants x conditions) and \code{LD}.
#'
#' Follows the \code{*_rss_weights(stat, LD, ...)} contract. Expects
#' \code{stat$z} to be a numeric matrix (variants x conditions) and
#' \code{stat$n} a per-context vector or scalar.
#'
#' @param stat A list with \code{z} (matrix variants x conditions) and
#'   \code{n} (numeric vector or scalar).
#' @param LD LD correlation matrix.
#' @param mvsusieRssFit Optional pre-fitted \code{mvsusieRss} object.
#' @param priorVariance Optional mvSuSiE prior variance specification.
#'   When NULL, \code{mvsusieR::create_mixture_prior()} is used with
#'   \code{R = ncol(stat$z)}.
#' @param residualVariance Optional residual covariance matrix.
#' @param L Maximum number of single effects (default 30).
#' @param LGreedy Initial greedy effect count (default 5).
#' @param retainFit If TRUE, attaches the fitted object as an attribute.
#' @param ... Additional arguments forwarded to \code{mvsusieR::mvsusieRss}.
#'
#' @return A numeric matrix of per-variant per-context weights
#'   (variants x conditions).
#' @export
mvsusieRssWeights <- function(stat, LD, mvsusieRssFit = NULL,
                              priorVariance = NULL,
                              residualVariance = NULL,
                              L = 30, LGreedy = 5,
                              retainFit = FALSE, ...) {
  if (!requireNamespace("mvsusieR", quietly = TRUE)) {
    stop("Package 'mvsusieR' is required. ",
         "Install with: devtools::install_github('stephenslab/mvsusieR')")
  }
  if (is.null(mvsusieRssFit)) {
    Z <- if (is.matrix(stat$z)) stat$z else as.matrix(stat$z)
    if (ncol(Z) < 2) {
      stop("mvsusieRssWeights expects stat$z to have >= 2 columns ",
           "(one per context). For single-context use susieRssWeights().")
    }
    # mvsusieR::mvsusieRss expects N to be a single scalar
    nScalar <- as.numeric(stats::median(stat$n))
    if (is.null(priorVariance)) {
      priorVariance <- mvsusieR::create_mixture_prior(R = ncol(Z))
    }
    if (!is.null(LGreedy)) LGreedy <- min(LGreedy, L)
    mvsusieRssFit <- mvsusieR::mvsusieRss(
      Z = Z, R = LD, N = nScalar,
      prior_variance = priorVariance,
      residual_variance = residualVariance, ...
    )
  }
  weights <- mvsusieR::coef.mvsusie(mvsusieRssFit)[-1, , drop = FALSE]
  if (retainFit) attr(weights, "fit") <- mvsusieRssFit
  weights
}

# Get a reasonable setting for the standard deviations of the mixture
# components in the mixture-of-normals prior based on the data (X, y).
# Input se is an estimate of the residual *variance*, and n is the
# number of standard deviations to return. This code is adapted from
# the autoselect.mixsd function in the ashr package.
#' @importFrom susieR univariate_regression
initPriorSd <- function(X, y, n = 30) {
  res <- univariate_regression(X, y)
  smax <- 3 * max(res$betahat)
  seq(0, smax, length.out = n)
}

# Identify zero-variance columns of X and warn the caller before they are
# dropped. Returns a logical vector of length ncol(X) where TRUE indicates a
# column to keep. The downstream solvers in pecotmr's regression wrappers
# (glmnet, ncvreg, L0Learn, qgg, BGLR, RcppDPR) all either error or behave
# poorly on constant columns, so wrappers should filter them out and zero-pad
# their results back to length p.
#' @importFrom matrixStats colSds
.dropZeroVariance <- function(X, fnName) {
  sds <- colSds(X)
  keep <- !is.na(sds) & sds != 0
  if (!all(keep)) {
    warning(sprintf(
      "%s: dropping %d zero-variance column(s) from X (indices: %s)",
      fnName, sum(!keep),
      paste(which(!keep), collapse = ", ")
    ), call. = FALSE)
  }
  keep
}

#' @importFrom stats coef
#' @export
glmnetWeights <- function(X, y, alpha) {
  # Check if glmnet is installed
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("To use this function, please install glmnet: https://cran.r-project.org/web/packages/glmnet/index.html")
  }
  eff.wgt <- matrix(0, ncol = 1, nrow = ncol(X))
  keep <- .dropZeroVariance(X, "glmnetWeights")
  enet <- glmnet::cv.glmnet(x = X[, keep, drop = FALSE], y = y, alpha = alpha, nfold = 5, intercept = TRUE, standardize = FALSE)
  eff.wgt[keep] <- coef(enet, s = "lambda.min")[2:(sum(keep) + 1)]
  return(eff.wgt)
}

#' @export
enetWeights <- function(X, y) glmnetWeights(X, y, 0.5)

#' @export
lassoWeights <- function(X, y) glmnetWeights(X, y, 1)

#' Compute Weights Using mr.ash Shrinkage
#'
#' This function fits the `mr.ash` model (adaptive shrinkage regression) to estimate weights
#' for a given set of predictors and response. It uses optional prior standard deviation initialization
#' and can accept custom initial beta values.
#'
#' @examples
#' wgt.mr.ash <- mrashWeights(eqtl$X, eqtl$y_res, beta.init = lassoWeights(X, y))
#' @importFrom susieR mr.ash
#' @importFrom stats predict
#' @export
mrashWeights <- function(X, y, initPriorSd = TRUE, retainFit = FALSE, ...) {
  eff.wgt <- rep(0, ncol(X))
  keep <- .dropZeroVariance(X, "mrashWeights")
  XKeep <- X[, keep, drop = FALSE]
  argsList <- list(...)
  if (!"beta.init" %in% names(argsList)) {
    argsList$beta.init <- lassoWeights(XKeep, y)
  } else if (length(argsList$beta.init) == ncol(X)) {
    argsList$beta.init <- argsList$beta.init[keep]
  }
  fit.mr.ash <- do.call(mr.ash, c(list(X = XKeep, y = y, sa2 = if (initPriorSd) initPriorSd(XKeep, y)^2 else NULL), argsList))
  eff.wgt[keep] <- predict(fit.mr.ash, type = "coefficients")[-1]
  if (retainFit) attr(eff.wgt, "fit") <- fit.mr.ash
  return(eff.wgt)
}
#' Extract Coefficients From Bayesian Linear Regression
#'
#' This function performs Bayesian linear regression using the `gbayes` function from
#' the `qgg` package. It then returns the estimated slopes.
#'
#' @param y A numeric vector of phenotypes.
#' @param X A numeric matrix of genotypes.
#' @param method A character string declaring the method/prior to be used. Options are
#' bayesN, bayesL, bayesA, bayesC, or bayesR.
#' @param Z An optional numeric matrix of covariates.
#' @return A vector containing the weights to be applied to each genotype in
#'   predicting the phenotype.
#' @details This function fits a Bayesian linear regression model with a range of priors.
#' @examples
#' X <- matrix(rnorm(100000), nrow = 1000)
#' Z <- matrix(round(runif(3000, 0, 0.8), 0), nrow = 1000)
#' set1 <- sample(1:ncol(X), 5)
#' set2 <- sample(1:ncol(X), 5)
#' sets <- list(set1, set2)
#' g <- rowSums(X[, c(set1, set2)])
#' e <- rnorm(nrow(X), mean = 0, sd = 1)
#' y <- g + e
#' bayesLWeights(y = y, X = X, Z = Z)
#' bayesRWeights(y = y, X = X, Z = Z)
#' @export
bayesAlphabetWeights <- function(X, y, method, Z = NULL, h2 = NULL, nit = 5000, nburn = 1000, nthin = 5, ...) {
  # Make sure qgg is installed
  if (!requireNamespace("qgg", quietly = TRUE)) {
    stop("To use this function, please install qgg: https://cran.r-project.org/web/packages/qgg/index.html")
  }
  # check for identical row lengths of response and genotype
  if (!(length(y) == nrow(X))) {
    stop("All objects must have the same number of rows")
  }
  # check for identical row lengths of genotype and covariates
  if (!is.null(Z)) {
    if (nrow(X) != nrow(Z)) {
      stop("Genotype and covariate matrices must have same number of rows")
    }
  }

  eff.wgt <- rep(0, ncol(X))
  keep <- .dropZeroVariance(X, "bayesAlphabetWeights")

  model <- qgg::gbayes(
    y = y,
    W = X[, keep, drop = FALSE],
    X = Z,
    method = method,
    h2 = h2,
    nit = nit,
    nburn = nburn,
    ...
  )

  eff.wgt[keep] <- model$bm
  return(eff.wgt)
}
#' Use Gaussian distribution as prior. Posterior means will be BLUP, equivalent to Ridge Regression.
#' @export
bayesNWeights <- function(X, y, Z = NULL, ...) {
  return(bayesAlphabetWeights(X, y, method = "bayesN", Z, ...))
}
#' Use laplace/double exponential distribution as prior. This is equivalent to Bayesian LASSO.
#' @export
bayesLWeights <- function(X, y, Z = NULL, ...) {
  return(bayesAlphabetWeights(X, y, method = "bayesL", Z, ...))
}
#' Use t-distribution as prior.
#' @export
bayesAWeights <- function(X, y, Z = NULL, ...) {
  return(bayesAlphabetWeights(X, y, method = "bayesA", Z, ...))
}
#' Use a rounded spike prior (low-variance Gaussian).
#' @export
bayesCWeights <- function(X, y, Z = NULL, pi = 0.1, ...) {
  return(bayesAlphabetWeights(X, y, method = "bayesC", Z, pi = pi, ...))
}
#' Use a hierarchical Bayesian mixture model with four Gaussian components. Variances are scaled
#' by 0, 0.0001 , 0.001 , and 0.01 .
#' @export
bayesRWeights <- function(X, y, Z = NULL, ...) {
  return(bayesAlphabetWeights(X, y, method = "bayesR", Z, ...))
}


# #' Bayesian linear regression using summary statistics
# #'
# #' @description
# #'
# #' This function is adapted from those written by Peter Sorensen in the qgg package.
# #' The following prior distributions are provided:
# #'
# #' Bayes N: Assigning a Gaussian prior to marker effects implies that the posterior means are the
# #' BLUP estimates (same as Ridge Regression).
# #'
# #' Bayes L: Assigning a double-exponential or Laplace prior is the density used in
# #' the Bayesian LASSO
# #'
# #' Bayes A: similar to ridge regression but t-distribution prior (rather than Gaussian)
# #' for the marker effects ; variance comes from an inverse-chi-square distribution instead of being fixed. Estimation
# #' via Gibbs sampling.
# #'
# #' Bayes C: uses a "rounded spike" (low-variance Gaussian) at origin many small
# #' effects can contribute to polygenic component, reduces the dimensionality of
# #' the model (makes Gibbs sampling feasible).
# #'
# #' Bayes R: Hierarchical Bayesian mixture model with 4 Gaussian components, with
# #' variances scaled by 0, 0.0001 , 0.001 , and 0.01 .
# #'
# #' @param sumstats dataframe with marker summary statistics. Required: beta coefficient (beta), standard
# #'        error of the beta coefficient (se), GWAS sample size (n). Optional: variant_id or rsid, alleles (A1
# #'        and A2), minor allele frequency (maf).
# #' @param LD is a the LD matrix corresponding to the same markers as in the stat dataframe
# #' @param variant_ids is an optional character vector of variant ids or rsids, provided outside of the rss dataframe
# #' @param nit is the number of iterations
# #' @param nburn is the number of burnin iterations
# #' @param nthin is the thinning parameter
# #' @param method specifies the methods used (method="bayesN","bayesA","bayesL","bayesC","bayesR")
# #' @param vg is a scalar or matrix of genetic (co)variances
# #' @param vb is a scalar or matrix of marker (co)variances
# #' @param ve is a scalar or matrix of residual (co)variances
# #' @param ssg_prior is a scalar or matrix of prior genetic (co)variances
# #' @param ssb_prior is a scalar or matrix of prior marker (co)variances
# #' @param sse_prior is a scalar or matrix of prior residual (co)variances
# #' @param lambda is a vector or matrix of lambda values
# #' @param h2 is the trait heritability
# #' @param pi is the proportion of markers in each marker variance class
# #' @param updateB is a logical for updating marker (co)variances
# #' @param updateG is a logical for updating genetic (co)variances
# #' @param updateE is a logical for updating residual (co)variances
# #' @param updatePi is a logical for updating pi
# #' @param adjustE is a logical for adjusting residual variance
# #' @param nug is a scalar or vector of prior degrees of freedom for prior genetic (co)variances
# #' @param nub is a scalar or vector of prior degrees of freedom for marker (co)variances
# #' @param nue is a scalar or vector of prior degrees of freedom for prior residual (co)variances
# #' @param mask is a vector or matrix of TRUE/FALSE specifying if marker should be ignored
# #' @param ve_prior is a scalar or matrix of prior residual (co)variances
# #' @param vg_prior is a scalar or matrix of prior genetic (co)variances
# #' @param algorithm is the algorithm to use. Should take on values ("mcmc", "em-mcmc")
# #' @param tol is tolerance, i.e. convergence criteria used in gbayes
# #' @param nit_local is the number of local iterations
# #' @param nit_global is the number of global iterations
# #'
# #' @return Returns a list structure including
# #' \item{bm}{vector of posterior means for marker effects}
# #' \item{dm}{vector of posterior means for marker inclusion probabilities}
# #' \item{vbs}{scalar or vector (t) of posterior means for marker variances}
# #' \item{vgs}{scalar or vector (t) of posterior means for genomic variances}
# #' \item{ves}{scalar or vector (t) of posterior means for residual variances}
# #' \item{pis}{vector of probabilites for each mcmc iteration}
# #' \item{pim}{posterior distribution probabilities}
# #' \item{r}{vector of residuals}
# #' \item{b}{vector of estimates from the final mcmc iteration}
# #' \item{param}{a list current parameters (same information as item listed above)
# #'              used for restart of the analysis}
# #' \item{stat}{matrix (mxt) of marker information and effects used for genomic risk scoring}
# #' \item{method}{the method used}
# #' \item{mask}{which loci were masked from analysis}
# #' \item{conv}{dataframe of convergence metrics}
# #' \item{post}{posterior parameter estimates}
# #' \item{ve}{mean residual variance}
# #' \item{vg}{mean genomic variance}
# #'
# #' @export
# gbayes_rss <- function(sumstats = NULL, LD = NULL, variant_ids = NULL, nit = 100, nburn = 0, nthin = 4, method = "bayesR",
#                        vg = NULL, vb = NULL, ve = NULL, ssg_prior = NULL, ssb_prior = NULL, sse_prior = NULL,
#                        lambda = NULL, h2 = NULL, pi = 0.001, updateB = TRUE, updateG = TRUE, updateE = TRUE,
#                        updatePi = TRUE, adjustE = TRUE, nug = 4, nub = 4, nue = 4, mask = NULL, ve_prior = NULL,
#                        vg_prior = NULL, algorithm = "mcmc", tol = 0.001, nit_local = NULL, nit_global = NULL) {
#   # Make sure qgg is installed
#   if (!requireNamespace("qgg", quietly = TRUE)) {
#     stop("To use this function, please install qgg: https://cran.r-project.org/web/packages/qgg/index.html")
#   }
#   # Check methods
#   methods <- c("bayesN", "bayesA", "bayesL", "bayesC", "bayesR")
#   method <- match(method, methods)
#   if (!sum(method %in% c(1:5)) == 1) stop("Method specified not valid")
#   if (method == 0) {
#     # BLUP and we do not estimate parameters
#     updateB <- FALSE
#     updateE <- FALSE
#   }
# 
#   # Set algorithm
#   if (algorithm == "em-mcmc") {
#     algo <- 2
#   } else {
#     algo <- 1
#   }
# 
#   # Check that LD matrix is provided and of same length as stats
#   if (is.null(LD)) stop("Must provide LD matrix")
#   if (nrow(sumstats) != nrow(LD)) stop("LD matrix must correspond to summary statistics")
# 
#   # Parameters from stat df
#   if (is.data.frame(sumstats)) {
#     if (!is.null(variant_ids)) {
#       variant_ids <- variant_ids
#     } else if (!is.null(sumstats$rsids)) {
#       variant_ids <- sumstats$rsids
#     } else if (!is.null(sumstats$variant_id)) {
#       variant_ids <- sumstats$variant_id
#     } else {
#       variant_ids <- paste0("snp", 1:nrow(sumstats))
#       sumstats$variant_id <- variant_ids
#     }
# 
#     m <- length(variant_ids)
#     b <- wy <- ww <- matrix(0, nrow = nrow(sumstats), ncol = 1)
#     mask <- matrix(FALSE, nrow = nrow(sumstats), ncol = 1)
#     rownames(b) <- rownames(wy) <- rownames(ww) <- rownames(mask) <- variant_ids
# 
#     if (is.null(sumstats$ww)) sumstats$ww <- 1 / (sumstats$se^2 + sumstats$beta^2 / sumstats$n)
#     if (is.null(sumstats$wy)) sumstats$wy <- sumstats$beta * sumstats$ww
#     if (!is.null(sumstats$n)) n <- as.integer(median(sumstats$n))
#     ww[, 1] <- sumstats$ww
#     wy[, 1] <- sumstats$wy
#     mask[, 1] <- FALSE
# 
#     if (any(is.na(wy))) stop("Missing values in wy")
#     if (any(is.na(ww))) stop("Missing values in ww")
# 
#     b2 <- sumstats$beta^2
#     seb2 <- sumstats$se^2
#     yy <- (b2 + (n - 2) * seb2) * sumstats$ww
#     yy <- median(yy)
# 
#     if (is.null(sumstats$A1)) sumstats$A1 <- rep("Unknown", length = nrow(sumstats))
#     if (is.null(sumstats$A2)) sumstats$A2 <- rep("Unknown", length = nrow(sumstats))
#     if (is.null(sumstats$maf)) {
#       sumstats$maf <- rep("Unknown", length = nrow(sumstats))
#       af_prov <- 0
#     } else {
#       af_prov <- 1
#     }
#   } else {
#     stop("Summary statistics must be provided in dataframe")
#   }
# 
# 
#   # prep LD for gbayes
#   LD_values <- lapply(1:nrow(LD), function(i) as.numeric(LD[i, ]))
#   names(LD_values) <- variant_ids
# 
#   LD_indices <- list(indices = vector("list", length = nrow(LD)))
#   for (i in 1:nrow(LD)) {
#     LD_indices[[i]] <- 1:nrow(LD) - 1
#   }
# 
#   bm <- dm <- fit <- res <- vector(length = 1, mode = "list")
#   names(bm) <- names(dm) <- names(fit) <- names(res) <- 1
# 
#   # Set parameters if not otherwise specified
#   if (is.null(m)) m <- length(LD_values)
#   vy <- yy / (n - 1)
#   if (is.null(pi)) pi <- 0.001
#   if (is.null(h2)) h2 <- 0.5
#   if (is.null(ve)) ve <- vy * (1 - h2)
#   if (is.null(vg)) vg <- vy * h2
#   if (method < 4 && is.null(vb)) vb <- vg / m
#   if (method >= 4 && is.null(vb)) vb <- vg / (m * pi)
#   if (is.null(lambda)) lambda <- rep(ve / vb, m)
#   if (method < 4 && is.null(ssb_prior)) ssb_prior <- ((nub - 2.0) / nub) * (vg / m)
#   if (method >= 4 && is.null(ssb_prior)) ssb_prior <- ((nub - 2.0) / nub) * (vg / (m * pi))
#   if (is.null(sse_prior)) sse_prior <- ((nue - 2.0) / nue) * ve
#   if (is.null(b)) b <- rep(0, m)
# 
#   pi <- c(1 - pi, pi)
#   gamma <- c(0, 1.0)
#   if (method == 5) pi <- c(0.95, 0.02, 0.02, 0.01)
#   if (method == 5) gamma <- c(0, 0.01, 0.1, 1.0)
# 
#   seed <- sample.int(.Machine$integer.max, 1)
# 
#   fit <- qgg:::sbayes_spa(
#     wy = wy,
#     ww = ww,
#     LDvalues = LD_values,
#     LDindices = LD_indices,
#     b = b,
#     lambda = lambda,
#     mask = mask,
#     yy = yy,
#     pi = pi,
#     gamma = gamma,
#     vg = vg,
#     vb = vb,
#     ve = ve,
#     ssb_prior = ssb_prior,
#     sse_prior = sse_prior,
#     nub = nub,
#     nue = nue,
#     updateB = updateB,
#     updateE = updateE,
#     updatePi = updatePi,
#     updateG = updateG,
#     adjustE = adjustE,
#     n = n,
#     nit = nit,
#     nburn = nburn,
#     nthin = nthin,
#     algo = algo,
#     method = as.integer(method),
#     seed = seed
#   )
# 
#   names(fit[[1]]) <- names(LD_values)
#   names(fit) <- c("bm", "dm", "coef", "vbs", "vgs", "ves", "pis", "pim", "r", "b", "param")
#   fit[3] <- NULL
# 
#   res <- data.frame(
#     variant_id = variant_ids, bm = fit$bm, dm = fit$dm,
#     pos = sumstats$pos, A1 = sumstats$A1,
#     A2 = sumstats$A2, maf = sumstats$maf,
#     stringsAsFactors = FALSE
#   )
#   rownames(res) <- variant_ids
# 
#   fit$sumstats <- res
#   if (af_prov == 1) {
#     fit$sumstats$vm <- 2 * (1 - fit$sumstats$maf) * fit$sumstats$maf * fit$sumstats$bm^2
#   }
#   fit$method <- methods[method]
#   fit$mask <- mask
# 
#   zve <- coda::geweke.diag(fit$ves[nburn:length(fit$ves)])$z
#   zvg <- coda::geweke.diag(fit$vgs[nburn:length(fit$vgs)])$z
#   zvb <- coda::geweke.diag(fit$vbs[nburn:length(fit$vbs)])$z
#   zpi <- coda::geweke.diag(fit$pis[nburn:length(fit$pis)])$z
# 
#   ve <- mean(fit$ves[nburn:length(fit$ves)])
#   vg <- mean(fit$vgs[nburn:length(fit$vgs)])
#   vb <- mean(fit$vbs[nburn:length(fit$vbs)])
#   pi <- 1 - fit$pim[1]
#   fit$conv <- data.frame(zve = zve, zvg = zvg, zvb = zvb, zpi = zpi)
#   fit$post <- data.frame(ve = ve, vg = vg, vb = vb, pi = pi)
#   fit$ve <- mean(ve)
#   fit$vg <- sum(vg)
# 
#   return(fit)
# }
# #' Extract weights from gbayes_rss function
# #' @return A numeric vector of the posterior mean of the coefficients.
# #' @export
# bayes_alphabet_rss_weights <- function(sumstats, LD, method, ...) {
#   model <- gbayes_rss(sumstats = sumstats, LD = LD, method = method, ...)
#   return(model$bm)
# }
# #' Use Gaussian distribution as prior. Posterior means will be BLUP, equivalent to Ridge Regression.
# #' @export
# bayes_n_rss_weights <- function(sumstats, LD, ...) {
#   return(bayes_alphabet_rss_weights(sumstats, LD, method = "bayesN", ...))
# }
# #' Use laplace/double exponential distribution as prior. This is equivalent to Bayesian LASSO.
# #' @export
# bayes_l_rss_weights <- function(sumstats, LD, ...) {
#   return(bayes_alphabet_rss_weights(sumstats, LD, method = "bayesL", ...))
# }
# #' Use t-distribution as prior.
# #' @export
# bayes_a_rss_weights <- function(sumstats, LD, ...) {
#   return(bayes_alphabet_rss_weights(sumstats, LD, method = "bayesA", ...))
# }
# #' Use a rounded spike prior (low-variance Gaussian).
# #' @export
# bayes_c_rss_weights <- function(sumstats, LD, ...) {
#   return(bayes_alphabet_rss_weights(sumstats, LD, method = "bayesC", ...))
# }
# #' Use a hierarchical Bayesian mixture model with four Gaussian components. Variances are scaled
# #' by 0, 0.0001 , 0.001 , and 0.01 .
# #' @export
# bayes_r_rss_weights <- function(sumstats, LD, ...) {
#   return(bayes_alphabet_rss_weights(sumstats, LD, method = "bayesR", ...))
# }

#' Lassosum RSS: LASSO on summary statistics with LD reference
#'
#' Coordinate descent to solve the penalized regression on summary statistics:
#' \deqn{f(\beta) = \beta' R \beta - 2\beta' r + 2\lambda ||\beta||_1}
#' where \eqn{R} is the LD matrix (pre-shrunk if desired) and \eqn{r = \hat\beta / \sqrt{n}}.
#'
#' Based on Mak et al (2017) "Polygenic scores via penalized regression on summary statistics",
#' Genetic Epidemiology 41(6):469-480.
#'
#' @param bhat A vector of marginal effect sizes.
#' @param LD A list of LD blocks, where each element is a matrix representing an LD block.
#'   If shrinkage is desired, apply it before passing (e.g., \code{(1-s)*R + s*I}).
#' @param n Sample size of the GWAS.
#' @param lambda A vector of L1 penalty values. Default: 20 values from 0.001 to 0.1 on log scale.
#' @param thr Convergence threshold. Default: 1e-4.
#' @param maxiter Maximum number of iterations. Default: 10000.
#'
#' @return A list containing:
#'   \item{betaEst}{Posterior estimates of SNP effect sizes at best lambda.}
#'   \item{beta}{Matrix of estimates (p x nlambda).}
#'   \item{lambda}{The lambda values used.}
#'   \item{conv}{Convergence indicators (1 = converged).}
#'   \item{loss}{Quadratic loss at each lambda.}
#'   \item{fbeta}{Full objective value at each lambda.}
#'   \item{nparams}{Number of non-zero coefficients at each lambda.}
#'
#' @examples
#' set.seed(42)
#' p <- 10
#' n <- 100
#' bhat <- rnorm(p, sd = 0.1)
#' R <- diag(p)
#' for (i in 1:(p - 1)) {
#'   R[i, i + 1] <- 0.3
#'   R[i + 1, i] <- 0.3
#' }
#' LD <- list(blk1 = R)
#' out <- lassosumRss(bhat, LD, n)
#' @export
lassosumRss <- function(bhat, LD, n,
                        lambda = exp(seq(log(0.0001), log(0.1), length.out = 20)),
                        thr = 1e-4, maxiter = 10000) {
  if (!is.list(LD)) {
    stop("Please provide a valid list of LD blocks using 'LD'.")
  }
  if (missing(n) || n <= 0) {
    stop("Please provide a valid sample size using 'n'.")
  }
  totalRowsInLd <- sum(sapply(LD, nrow))
  if (length(bhat) != totalRowsInLd) {
    stop("The length of 'bhat' must be the same as the sum of the number of rows of all elements in the 'LD' list.")
  }

  z <- bhat / sqrt(n)
  order <- order(lambda, decreasing = TRUE)
  # cpp11 requires exact integer types for int parameters
  result <- lassosumRssRcpp(zR = z, LD = LD, lambdaR = lambda[order],
                            thr = thr, maxiter = as.integer(maxiter))

  # Reorder back to original lambda order.
  # Must use inverse permutation to unsort: if order[i]=j, then
  # the result at position j in the sorted output goes to position i.
  invOrder <- order(order)
  result$beta  <- result$beta[, invOrder, drop = FALSE]
  result$conv  <- result$conv[invOrder]
  result$loss  <- result$loss[invOrder]
  result$fbeta <- result$fbeta[invOrder]
  result$lambda <- lambda
  result$nparams <- as.integer(colSums(result$beta != 0))
  result$betaEst <- as.numeric(result$beta[, which.min(result$fbeta)])
  result
}

.lassosumCorFromStat <- function(stat, n, p) {
  corInput <- if (!is.null(stat$cor)) {
    as.numeric(stat$cor)
  } else if (!is.null(stat$z)) {
    as.numeric(stat$z) / sqrt(n)
  } else if (!is.null(stat$b)) {
    as.numeric(stat$b)
  } else {
    stop("stat must contain one of 'cor', 'z', or 'b' for lassosum selection.")
  }
  if (length(corInput) != p) {
    stop("The length of lassosum input statistics (", length(corInput),
         ") must equal nrow(LD) (", p, ").")
  }
  corInput
}

.lassosumClampCor <- function(corInput) {
  maxAbsCor <- max(abs(corInput), na.rm = TRUE)
  if (is.finite(maxAbsCor) && maxAbsCor >= 1) {
    corInput <- corInput / (maxAbsCor / 0.9999)
  }
  corInput
}

.lassosumFirstMax <- function(x) {
  which(x == max(x, na.rm = TRUE))[1]
}

.lassosumSelectMinFbeta <- function(candidateBeta, candidateMeta) {
  idx <- which.min(candidateMeta$fbeta)
  list(
    beta = candidateBeta[, idx],
    index = idx,
    mode = "min_fbeta"
  )
}

.lassosumSelectLdQuadratic <- function(candidateBeta, corInput, LD) {
  ldBeta <- LD %*% candidateBeta
  bxy <- as.numeric(crossprod(corInput, candidateBeta))
  bxxb <- colSums(candidateBeta * ldBeta)
  scores <- rep(-Inf, length(bxy))
  positive <- is.finite(bxxb) & bxxb > 0
  scores[positive] <- bxy[positive] / sqrt(bxxb[positive])
  idx <- .lassosumFirstMax(scores)
  list(
    beta = candidateBeta[, idx],
    index = idx,
    mode = "ld_quadratic"
  )
}

#' Extract weights from lassosumRss with shrinkage grid search
#'
#' Searches over a grid of shrinkage parameters \code{s} (default:
#' \code{c(0.2, 0.5, 0.9, 1.0)}, matching the original lassosum and OTTERS).
#' For each \code{s}, the LD matrix is shrunk as \code{(1-s)*R + s*I}, then
#' \code{lassosumRss()} is called across the lambda path. Candidate selection
#' defaults to the LD-only quadratic pseudovalidation score
#' \deqn{\frac{c^T \beta}{\sqrt{\beta^T R \beta}}}
#' evaluated on the supplied LD matrix \code{R}. This uses the same candidate
#' beta path as \code{lassosumRss()}, but scores each candidate directly from
#' summary-statistics correlation \code{c} and LD, without requiring genotype.
#'
#' @details
#' The original lassosum pseudovalidation can be written as an LD quadratic
#' score after centering and standardizing the reference matrix columns by the
#' same per-variant scale:
#' \deqn{\mathrm{score}(\beta) = \frac{c^T \beta}{\sqrt{\beta^T R \beta}}.}
#' This implementation therefore uses the supplied LD matrix directly for
#' selection. \code{min(fbeta)} is retained only as an explicit debug option.
#'
#' @param stat A list with \code{$b} (effect sizes) and \code{$n} (per-variant sample sizes).
#' @param LD LD correlation matrix R (single matrix, NOT pre-shrunk).
#' @param s Numeric vector of shrinkage parameters to search over. Default:
#'   \code{c(0.2, 0.5, 0.9, 1.0)} following Mak et al (2017) and OTTERS.
#' @param selection Selection strategy. Default \code{"ld_quadratic"} uses
#'   \eqn{c^T \beta / \sqrt{\beta^T R \beta}} on the supplied LD matrix.
#'   \code{"min_fbeta"} is retained as an explicit alternative for debugging.
#' @param ... Additional arguments passed to \code{lassosumRss()}.
#'
#' @return A numeric vector of the posterior SNP coefficients at the best (s, lambda).
#' @export
lassosumRssWeights <- function(stat, LD, s = c(0.2, 0.5, 0.9, 1.0),
                               selection = c("ld_quadratic", "min_fbeta"),
                               ...) {
  selection <- match.arg(selection)
  n <- median(stat$n)
  p <- nrow(LD)
  corInput <- .lassosumClampCor(.lassosumCorFromStat(stat, n = n, p = p))
  solverInput <- corInput * sqrt(n)
  candidateBeta <- NULL
  candidateMeta <- list()

  for (sVal in s) {
    LDs <- (1 - sVal) * LD + sVal * diag(p)
    model <- lassosumRss(bhat = solverInput, LD = list(blk1 = LDs), n = n, ...)
    candidateBeta <- cbind(candidateBeta, model$beta)
    candidateMeta[[length(candidateMeta) + 1L]] <- data.frame(
      s = rep(sVal, length(model$lambda)),
      lambda = model$lambda,
      fbeta = model$fbeta,
      stringsAsFactors = FALSE
    )
  }
  candidateMeta <- do.call(rbind, candidateMeta)

  selectorResult <- if (selection == "ld_quadratic") {
    .lassosumSelectLdQuadratic(candidateBeta, corInput, LD)
  } else {
    .lassosumSelectMinFbeta(candidateBeta, candidateMeta)
  }

  bestBeta <- as.numeric(selectorResult$beta)
  attr(bestBeta, "lassosum_selection") <- c(
    mode = selectorResult$mode,
    index = selectorResult$index,
    s = candidateMeta$s[selectorResult$index],
    lambda = candidateMeta$lambda[selectorResult$index]
  )
  bestBeta
}

#' Penalized Regression on RSS (Summary Statistics) Objective
#'
#' Generalizes \code{lassosumRss()} to support LASSO, MCP, SCAD, L0, L0L1,
#' and L0L2 penalties.  Uses coordinate descent on the objective
#' \deqn{\beta^T R \beta - 2 \beta^T z + \mathrm{penalty}(\beta)}
#' where \eqn{R} is a (possibly pre-shrunk) LD matrix and \eqn{z = \hat\beta / \sqrt{n}}.
#'
#' @param bhat Numeric vector of marginal effect estimates (length p).
#' @param LD A list of LD correlation matrices (one per block), as in
#'   \code{lassosumRss()}.
#' @param n GWAS sample size (positive scalar).
#' @param penalty Penalty type: \code{"lasso"}, \code{"MCP"}, \code{"SCAD"},
#'   \code{"L0"}, \code{"L0L1"}, or \code{"L0L2"}.
#' @param lambda Numeric vector of regularization parameter values along which
#'   to trace a solution path (warm-started, largest-first).  For LASSO/MCP/SCAD
#'   this is the primary penalty strength; for L0 variants it controls the L1
#'   component.
#' @param gamma Concavity parameter for MCP (default 3) or SCAD (default 3.7).
#'   Ignored for LASSO and L0 variants.
#' @param alpha Elastic-net mixing for MCP/SCAD: \eqn{l_1 = \lambda \alpha},
#'   \eqn{l_2 = \lambda (1-\alpha)}.  Default 1 (pure L1, no ridge).
#' @param lambda0 L0 penalty weight (number of non-zeros).  Required for L0
#'   variants; ignored otherwise.  Default 0.
#' @param lambda2 L2 penalty weight for L0L2 variant.  Default 0.
#' @param thr Convergence threshold.  Default 1e-4.
#' @param maxiter Maximum coordinate descent iterations per lambda.  Default 10000.
#' @param max_swaps Maximum swap rounds for L0 variants.  Default 100.
#'   Set to 0 to disable swaps.
#'
#' @return A list with components:
#' \describe{
#'   \item{beta}{p x length(lambda) matrix of coefficient estimates.}
#'   \item{lambda}{The lambda values used.}
#'   \item{conv}{Convergence indicators (1 = converged).}
#'   \item{loss}{Quadratic loss at each lambda.}
#'   \item{fbeta}{Full penalized objective at each lambda.}
#'   \item{nparams}{Number of non-zero coefficients at each lambda.}
#'   \item{betaEst}{Coefficient vector at the lambda minimizing fbeta.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' p <- 10; n <- 100
#' bhat <- rnorm(p, sd = 0.1)
#' R <- diag(p)
#' # MCP
#' penalizedRss(bhat, list(blk1 = R), n, penalty = "MCP")
#' # SCAD
#' penalizedRss(bhat, list(blk1 = R), n, penalty = "SCAD")
#' # L0
#' penalizedRss(bhat, list(blk1 = R), n, penalty = "L0", lambda0 = 0.01,
#'               lambda = c(0))
#' }
#' @export
penalizedRss <- function(bhat, LD, n,
                         penalty = c("lasso", "MCP", "SCAD", "L0", "L0L1", "L0L2"),
                         lambda = exp(seq(log(0.0001), log(0.1), length.out = 20)),
                         gamma = NULL, alpha = 1.0,
                         lambda0 = 0, lambda2 = 0,
                         thr = 1e-4, maxiter = 10000, maxSwaps = 100) {
  penalty <- match.arg(penalty)
  if (!is.list(LD)) {
    stop("Please provide a valid list of LD blocks using 'LD'.")
  }
  if (missing(n) || n <= 0) {
    stop("Please provide a valid sample size using 'n'.")
  }
  totalRowsInLd <- sum(sapply(LD, nrow))
  if (length(bhat) != totalRowsInLd) {
    stop("The length of 'bhat' must be the same as the sum of the number of rows of all elements in the 'LD' list.")
  }

  # Default gamma per penalty
  if (is.null(gamma)) {
    gamma <- switch(penalty,
                    SCAD = 3.7,
                    MCP  = 3.0,
                    0.0)
  }

  z <- bhat / sqrt(n)
  order <- order(lambda, decreasing = TRUE)

  result <- penalizedRssRcpp(zR = z, LD = LD, lambdaR = lambda[order],
                             penaltyStr = penalty,
                             gamma = gamma, alpha = alpha,
                             lambda0 = lambda0, lambda2 = lambda2,
                             thr = thr, maxiter = as.integer(maxiter),
                             maxSwaps = as.integer(maxSwaps))

  # Reorder back to original lambda order
  invOrder <- order(order)
  result$beta   <- result$beta[, invOrder, drop = FALSE]
  result$conv   <- result$conv[invOrder]
  result$loss   <- result$loss[invOrder]
  result$fbeta  <- result$fbeta[invOrder]
  result$lambda <- lambda
  result$nparams <- as.integer(colSums(result$beta != 0))
  result$betaEst <- as.numeric(result$beta[, which.min(result$fbeta)])
  result
}

#' RSS Weights Helper for Penalized Methods
#'
#' Shared implementation for \code{scadRssWeights()}, \code{mcpRssWeights()},
#' and \code{l0learnRssWeights()}.  Searches over a shrinkage grid \code{s}
#' (LD matrix shrinkage \code{(1-s)R + sI}) and selects the best candidate via
#' LD-quadratic pseudovalidation or minimum penalized objective.
#'
#' @param stat,LD,s,selection,penalty,gamma,alpha,lambda0,lambda2,...
#'   See the public wrappers for details.
#' @return Numeric weight vector of length \code{nrow(LD)}.
#' @keywords internal
.penalizedRssWeights <- function(stat, LD, penalty,
                                 s = c(0.2, 0.5, 0.9, 1.0),
                                 gamma = NULL, alpha = 1.0,
                                 lambda0 = 0, lambda2 = 0,
                                 selection = c("ld_quadratic", "min_fbeta"),
                                 ...) {
  selection <- match.arg(selection)
  n <- median(stat$n)
  p <- nrow(LD)
  corInput <- .lassosumClampCor(.lassosumCorFromStat(stat, n = n, p = p))
  solverInput <- corInput * sqrt(n)
  candidateBeta <- NULL
  candidateMeta <- list()

  for (sVal in s) {
    LDs <- (1 - sVal) * LD + sVal * diag(p)
    model <- penalizedRss(bhat = solverInput, LD = list(blk1 = LDs), n = n,
                          penalty = penalty, gamma = gamma, alpha = alpha,
                          lambda0 = lambda0, lambda2 = lambda2, ...)
    candidateBeta <- cbind(candidateBeta, model$beta)
    candidateMeta[[length(candidateMeta) + 1L]] <- data.frame(
      s = rep(sVal, length(model$lambda)),
      lambda = model$lambda,
      fbeta = model$fbeta,
      stringsAsFactors = FALSE
    )
  }
  candidateMeta <- do.call(rbind, candidateMeta)

  selectorResult <- if (selection == "ld_quadratic") {
    .lassosumSelectLdQuadratic(candidateBeta, corInput, LD)
  } else {
    .lassosumSelectMinFbeta(candidateBeta, candidateMeta)
  }

  bestBeta <- as.numeric(selectorResult$beta)
  attr(bestBeta, "penalized_rss_selection") <- c(
    mode = selectorResult$mode,
    index = selectorResult$index,
    penalty = penalty,
    s = candidateMeta$s[selectorResult$index],
    lambda = candidateMeta$lambda[selectorResult$index]
  )
  bestBeta
}

#' Compute SCAD-Penalized Weights from Summary Statistics
#'
#' Fits SCAD-penalized regression on the RSS objective, searching over a
#' shrinkage grid \code{s} and lambda path.  Model selection uses LD-quadratic
#' pseudovalidation by default.
#'
#' @param stat A list with \code{$b} (effect sizes) and \code{$n} (per-variant sample sizes).
#' @param LD LD correlation matrix R (single matrix, NOT pre-shrunk).
#' @param s Numeric vector of LD shrinkage parameters.  Default:
#'   \code{c(0.2, 0.5, 0.9, 1.0)}.
#' @param gamma SCAD concavity parameter.  Default 3.7.
#' @param alpha Elastic-net mixing (1 = pure L1).  Default 1.
#' @param selection Selection strategy: \code{"ld_quadratic"} (default) or
#'   \code{"min_fbeta"}.
#' @param ... Additional arguments passed to \code{penalizedRss()}.
#' @return A numeric vector of SNP coefficient weights.
#' @export
scadRssWeights <- function(stat, LD, s = c(0.2, 0.5, 0.9, 1.0),
                           gamma = 3.7, alpha = 1.0,
                           selection = c("ld_quadratic", "min_fbeta"), ...) {
  .penalizedRssWeights(stat = stat, LD = LD, penalty = "SCAD",
                       s = s, gamma = gamma, alpha = alpha,
                       selection = selection, ...)
}

#' Compute MCP-Penalized Weights from Summary Statistics
#'
#' Fits MCP-penalized regression on the RSS objective, searching over a
#' shrinkage grid \code{s} and lambda path.  Model selection uses LD-quadratic
#' pseudovalidation by default.
#'
#' @param stat A list with \code{$b} (effect sizes) and \code{$n} (per-variant sample sizes).
#' @param LD LD correlation matrix R (single matrix, NOT pre-shrunk).
#' @param s Numeric vector of LD shrinkage parameters.  Default:
#'   \code{c(0.2, 0.5, 0.9, 1.0)}.
#' @param gamma MCP concavity parameter.  Default 3.
#' @param alpha Elastic-net mixing (1 = pure L1).  Default 1.
#' @param selection Selection strategy: \code{"ld_quadratic"} (default) or
#'   \code{"min_fbeta"}.
#' @param ... Additional arguments passed to \code{penalizedRss()}.
#' @return A numeric vector of SNP coefficient weights.
#' @export
mcpRssWeights <- function(stat, LD, s = c(0.2, 0.5, 0.9, 1.0),
                          gamma = 3.0, alpha = 1.0,
                          selection = c("ld_quadratic", "min_fbeta"), ...) {
  .penalizedRssWeights(stat = stat, LD = LD, penalty = "MCP",
                       s = s, gamma = gamma, alpha = alpha,
                       selection = selection, ...)
}

#' Compute L0-Penalized Weights from Summary Statistics
#'
#' Fits L0-penalized regression (with optional L1/L2 components) on the RSS
#' objective, searching over a shrinkage grid \code{s} and lambda0 path.
#' Model selection uses LD-quadratic pseudovalidation by default.
#'
#' The swap optimization from L0Learn is included: after coordinate descent
#' converges, non-zero coefficients are tested for swaps with zero ones to
#' escape local optima.
#'
#' @param stat A list with \code{$b} (effect sizes) and \code{$n} (per-variant sample sizes).
#' @param LD LD correlation matrix R (single matrix, NOT pre-shrunk).
#' @param penalty L0 variant: \code{"L0"}, \code{"L0L1"}, or \code{"L0L2"}.
#'   Default \code{"L0"}.
#' @param s Numeric vector of LD shrinkage parameters.  Default:
#'   \code{c(0.2, 0.5, 0.9, 1.0)}.
#' @param lambda0 Numeric vector of L0 penalty values to search over.
#'   Default: \code{exp(seq(log(0.001), log(1), length.out = 10))}.
#' @param lambda Numeric vector of L1 penalty values (for L0L1).  Default:
#'   \code{c(0)} (no L1 unless L0L1 is used).
#' @param lambda2 L2 penalty weight (for L0L2).  Default 0.
#' @param selection Selection strategy: \code{"ld_quadratic"} (default) or
#'   \code{"min_fbeta"}.
#' @param maxSwaps Maximum swap rounds per lambda.  Default 100.
#' @param ... Additional arguments passed to \code{penalizedRss()}.
#' @return A numeric vector of SNP coefficient weights.
#' @export
l0learnRssWeights <- function(stat, LD,
                              penalty = c("L0", "L0L1", "L0L2"),
                              s = c(0.2, 0.5, 0.9, 1.0),
                              lambda0 = exp(seq(log(0.001), log(1), length.out = 10)),
                              lambda = NULL, lambda2 = 0,
                              selection = c("ld_quadratic", "min_fbeta"),
                              maxSwaps = 100, ...) {
  penalty <- match.arg(penalty)
  selection <- match.arg(selection)

  # Default lambda (L1 component) depends on variant
  if (is.null(lambda)) {
    lambda <- if (penalty == "L0L1") {
      exp(seq(log(0.0001), log(0.1), length.out = 10))
    } else {
      c(0)
    }
  }

  n <- median(stat$n)
  p <- nrow(LD)
  corInput <- .lassosumClampCor(.lassosumCorFromStat(stat, n = n, p = p))
  solverInput <- corInput * sqrt(n)
  candidateBeta <- NULL
  candidateMeta <- list()

  # Grid search over s and lambda0
  for (sVal in s) {
    LDs <- (1 - sVal) * LD + sVal * diag(p)
    for (l0Val in lambda0) {
      model <- penalizedRss(bhat = solverInput, LD = list(blk1 = LDs), n = n,
                            penalty = penalty, lambda = lambda,
                            lambda0 = l0Val, lambda2 = lambda2,
                            maxSwaps = maxSwaps, ...)
      candidateBeta <- cbind(candidateBeta, model$beta)
      candidateMeta[[length(candidateMeta) + 1L]] <- data.frame(
        s = rep(sVal, length(model$lambda)),
        lambda0 = rep(l0Val, length(model$lambda)),
        lambda = model$lambda,
        fbeta = model$fbeta,
        stringsAsFactors = FALSE
      )
    }
  }
  candidateMeta <- do.call(rbind, candidateMeta)

  selectorResult <- if (selection == "ld_quadratic") {
    .lassosumSelectLdQuadratic(candidateBeta, corInput, LD)
  } else {
    .lassosumSelectMinFbeta(candidateBeta, candidateMeta)
  }

  bestBeta <- as.numeric(selectorResult$beta)
  attr(bestBeta, "penalized_rss_selection") <- c(
    mode = selectorResult$mode,
    index = selectorResult$index,
    penalty = penalty,
    s = candidateMeta$s[selectorResult$index],
    lambda0 = candidateMeta$lambda0[selectorResult$index],
    lambda = candidateMeta$lambda[selectorResult$index]
  )
  bestBeta
}

#' Compute Weights Using ncvreg with SCAD or MCP Penalty
#'
#' Internal helper that fits an `ncvreg` model with the specified non-convex
#' penalty using k-fold cross-validation, then returns the coefficients at
#' `lambda.min`. Following the convention of `glmnetWeights`, columns of `X`
#' with zero (or `NA`) standard deviation are dropped before fitting and their
#' weights are set to zero.
#'
#' @param X A numeric matrix of predictors (no intercept column; `ncvreg`
#'   standardizes internally and adds its own intercept).
#' @param y A numeric response vector.
#' @param penalty Either "SCAD" or "MCP".
#' @param nfolds Number of cross-validation folds. Default is 5.
#' @param ... Additional arguments passed through to `ncvreg::cv.ncvreg`.
#' @return A numeric vector of length `ncol(X)` of variant weights.
#' @importFrom stats coef
#' @keywords internal
ncvregWeights <- function(X, y, penalty, nfolds = 5, ...) {
  if (!requireNamespace("ncvreg", quietly = TRUE)) {
    stop("To use this function, please install ncvreg: https://cran.r-project.org/package=ncvreg")
  }
  eff.wgt <- matrix(0, ncol = 1, nrow = ncol(X))
  keep <- .dropZeroVariance(X, "ncvregWeights")
  fit <- ncvreg::cv.ncvreg(X = X[, keep, drop = FALSE], y = y, penalty = penalty, nfolds = nfolds, ...)
  eff.wgt[keep] <- coef(fit, lambda = fit$lambda.min)[-1]
  return(eff.wgt)
}

#' Compute Weights Using SCAD-Penalized Regression
#'
#' Fits a SCAD-penalized linear regression model via `ncvreg::cv.ncvreg` and
#' returns the coefficient vector at `lambda.min`.
#'
#' @param X A numeric matrix of predictors.
#' @param y A numeric response vector.
#' @param nfolds Number of cross-validation folds. Default is 5.
#' @param ... Additional arguments passed through to `ncvreg::cv.ncvreg`.
#' @return A numeric vector of length `ncol(X)` of variant weights.
#' @export
scadWeights <- function(X, y, nfolds = 5, ...) {
  ncvregWeights(X, y, penalty = "SCAD", nfolds = nfolds, ...)
}

#' Compute Weights Using MCP-Penalized Regression
#'
#' Fits an MCP-penalized linear regression model via `ncvreg::cv.ncvreg` and
#' returns the coefficient vector at `lambda.min`.
#'
#' @param X A numeric matrix of predictors.
#' @param y A numeric response vector.
#' @param nfolds Number of cross-validation folds. Default is 5.
#' @param ... Additional arguments passed through to `ncvreg::cv.ncvreg`.
#' @return A numeric vector of length `ncol(X)` of variant weights.
#' @export
mcpWeights <- function(X, y, nfolds = 5, ...) {
  ncvregWeights(X, y, penalty = "MCP", nfolds = nfolds, ...)
}

#' Compute Weights Using L0Learn
#'
#' Fits an L0-regularized linear regression model via `L0Learn::L0Learn.cvfit`
#' and returns the coefficient vector at the (lambda, gamma) pair minimizing
#' the cross-validation error. Default penalty is "L0"; the user can switch to
#' "L0L1" or "L0L2" (and tune the corresponding gamma grid) by passing the
#' relevant arguments through `...`.
#'
#' @param X A numeric matrix of predictors.
#' @param y A numeric response vector.
#' @param penalty Type of regularization: "L0", "L0L1", or "L0L2". Default is "L0".
#' @param nFolds Number of cross-validation folds. Default is 5.
#' @param ... Additional arguments passed through to `L0Learn::L0Learn.cvfit`
#'   (e.g. `nGamma`, `gammaMin`, `gammaMax`, `algorithm`, `maxSuppSize`).
#' @return A numeric vector of length `ncol(X)` of variant weights.
#' @export
l0learnWeights <- function(X, y, penalty = "L0", nFolds = 5, ...) {
  if (!requireNamespace("L0Learn", quietly = TRUE)) {
    stop("To use this function, please install L0Learn: https://cran.r-project.org/package=L0Learn")
  }
  eff.wgt <- matrix(0, ncol = 1, nrow = ncol(X))
  keep <- .dropZeroVariance(X, "l0learnWeights")
  fit <- L0Learn::L0Learn.cvfit(
    x = X[, keep, drop = FALSE], y = y, penalty = penalty, nFolds = nFolds, ...
  )
  # Find (gamma, lambda) minimizing CV error across the entire path.
  cvMins <- vapply(fit$cvMeans, function(v) min(as.numeric(v)), numeric(1))
  gammaIdx <- which.min(cvMins)
  lambdaIdx <- which.min(as.numeric(fit$cvMeans[[gammaIdx]]))
  bestGamma <- fit$fit$gamma[gammaIdx]
  bestLambda <- fit$fit$lambda[[gammaIdx]][lambdaIdx]
  coefs <- as.numeric(coef(fit, lambda = bestLambda, gamma = bestGamma))
  # If intercept was included, drop it (first row).
  if (length(coefs) == sum(keep) + 1L) {
    coefs <- coefs[-1L]
  }
  eff.wgt[keep] <- coefs
  return(eff.wgt)
}

#' Compute Weights Using a BGLR Linear Regression Model
#'
#' Internal helper that fits a `BGLR::BGLR` linear regression with a single
#' linear term whose `model` is one of BGLR's marker-effect priors (e.g.
#' "BayesB", "BL"), then returns the posterior mean of the marker effects.
#' BGLR writes per-call temporary files to disk; this helper sandboxes them in
#' a fresh `tempdir()` that is cleaned up on exit.
#'
#' @param X A numeric matrix of predictors.
#' @param y A numeric response vector.
#' @param model A BGLR marker-effect model name (e.g. "BayesB" or "BL").
#' @param nIter Number of MCMC iterations.
#' @param burnIn Number of burn-in iterations.
#' @param thin Thinning interval.
#' @param etaArgs Optional named list of additional arguments included in the
#'   `ETA` linear-term specification (e.g. `list(probIn = 0.05)` for BayesB).
#' @param ... Additional arguments passed through to `BGLR::BGLR`.
#' @return A numeric vector of length `ncol(X)` of variant weights.
#' @keywords internal
bglrWeights <- function(X, y, model, nIter, burnIn, thin, etaArgs = list(), ...) {
  if (!requireNamespace("BGLR", quietly = TRUE)) {
    stop("To use this function, please install BGLR: https://cran.r-project.org/package=BGLR")
  }
  eff.wgt <- rep(0, ncol(X))
  keep <- .dropZeroVariance(X, "bglrWeights")

  tmpdir <- tempfile("bglr_")
  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
  saveAt <- paste0(tmpdir, .Platform$file.sep)

  eta <- list(c(list(X = X[, keep, drop = FALSE], model = model), etaArgs))
  fit <- BGLR::BGLR(
    y = y, ETA = eta,
    nIter = nIter, burnIn = burnIn, thin = thin,
    saveAt = saveAt, verbose = FALSE, ...
  )
  eff.wgt[keep] <- as.numeric(fit$ETA[[1]]$b)
  return(eff.wgt)
}

#' Compute Weights Using BayesB
#'
#' Fits a BayesB linear regression model via `BGLR::BGLR` and returns the
#' posterior mean of the marker effects. BayesB places a "spike-and-slab"
#' mixture prior on each marker effect, with a scaled-t slab.
#'
#' Defaults for `nIter`, `burnIn`, and `thin` are larger than BGLR's package
#' defaults to better accommodate the high LD typical of cis-eQTL windows; see
#' Kim et al. (2022) which observed that the BGLR defaults can be inadequate
#' under correlated predictors. Override these arguments to recover the
#' package defaults if desired.
#'
#' @param X A numeric matrix of predictors.
#' @param y A numeric response vector.
#' @param nIter Number of MCMC iterations. Default is 10000.
#' @param burnIn Number of burn-in iterations. Default is 2000.
#' @param thin Thinning interval. Default is 5.
#' @param probIn Prior inclusion probability for each marker. Default is 0.2.
#' @param ... Additional arguments passed through to `BGLR::BGLR`.
#' @return A numeric vector of length `ncol(X)` of variant weights.
#' @export
bayesBWeights <- function(X, y, nIter = 10000, burnIn = 2000, thin = 5, probIn = 0.2, ...) {
  bglrWeights(
    X, y,
    model = "BayesB", nIter = nIter, burnIn = burnIn, thin = thin,
    etaArgs = list(probIn = probIn), ...
  )
}

#' Compute Weights Using the Bayesian LASSO (BGLR)
#'
#' Fits a Bayesian LASSO linear regression model via `BGLR::BGLR` (the "BL"
#' model, Park & Casella 2008) and returns the posterior mean of the marker
#' effects. This is the same "B-Lasso" implementation benchmarked in Kim et
#' al. (2022). Note that this is distinct from `bayesLWeights`, which uses a
#' different Bayesian LASSO implementation backed by `qgg`.
#'
#' Defaults for `nIter`, `burnIn`, and `thin` are larger than BGLR's package
#' defaults to better accommodate high-LD cis-eQTL windows; override to
#' recover the package defaults.
#'
#' @param X A numeric matrix of predictors.
#' @param y A numeric response vector.
#' @param nIter Number of MCMC iterations. Default is 10000.
#' @param burnIn Number of burn-in iterations. Default is 2000.
#' @param thin Thinning interval. Default is 5.
#' @param ... Additional arguments passed through to `BGLR::BGLR`.
#' @return A numeric vector of length `ncol(X)` of variant weights.
#' @export
bLassoWeights <- function(X, y, nIter = 10000, burnIn = 2000, thin = 5, ...) {
  bglrWeights(
    X, y,
    model = "BL", nIter = nIter, burnIn = burnIn, thin = thin, ...
  )
}

#' Compute Weights Using Dirichlet Process Regression (RcppDPR)
#'
#' Fits a Dirichlet Process Regression model via `RcppDPR::fit_model` and
#' returns the per-variant weights, computed as `beta + alpha` (matching
#' `RcppDPR:::predict.DPR_Model`, which uses `(beta + alpha) %*% x_new + pheno_mean`).
#'
#' By default the variational Bayes (`VB`) fitting method is used, which is
#' fast and deterministic. The user may switch to `Gibbs` or `Adaptive_Gibbs`
#' for full Bayesian MCMC inference. `rotate_variables` is held to `FALSE`
#' under the assumption that any covariates have already been regressed out
#' upstream; an intercept-only covariate matrix is supplied to `fit_model`.
#'
#' @param X A numeric matrix of predictors.
#' @param y A numeric response vector.
#' @param fittingMethod One of "VB", "Gibbs", or "Adaptive_Gibbs". Default is "VB".
#' @param ... Additional arguments passed through to `RcppDPR::fit_model`.
#' @return A numeric vector of length `ncol(X)` of variant weights.
#' @export
dprWeights <- function(X, y, fittingMethod = "VB", retainFit = FALSE, ...) {
  if (!requireNamespace("RcppDPR", quietly = TRUE)) {
    stop("To use this function, please install RcppDPR: https://cran.r-project.org/package=RcppDPR")
  }
  eff.wgt <- rep(0, ncol(X))
  keep <- .dropZeroVariance(X, "dprWeights")
  w <- matrix(1, nrow = nrow(X), ncol = 1)
  fit <- RcppDPR::fit_model(
    y = y, w = w, x = X[, keep, drop = FALSE],
    rotate_variables = FALSE, fitting_method = fittingMethod, ...
  )
  eff.wgt[keep] <- as.numeric(fit$beta + fit$alpha)
  if (retainFit) attr(eff.wgt, "fit") <- fit
  return(eff.wgt)
}

#' @rdname dprWeights
#' @export
dprVbWeights <- function(X, y, nK = 8, retainFit = FALSE, ...) dprWeights(X, y, fittingMethod = "VB", n_k = nK, retainFit = retainFit, ...)

#' @rdname dprWeights
#' @export
dprGibbsWeights <- function(X, y, sStep = 5000, retainFit = FALSE, ...) dprWeights(X, y, fittingMethod = "Gibbs", s_step = sStep, retainFit = retainFit, ...)

#' @rdname dprWeights
#' @export
dprAdaptiveGibbsWeights <- function(X, y, retainFit = FALSE, ...) dprWeights(X, y, fittingMethod = "Adaptive_Gibbs", retainFit = retainFit, ...)
