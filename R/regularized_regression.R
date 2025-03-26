#' Bayesian Multiple Regression with Mixture-of-Normals Prior
#'
#' This function performs Bayesian multiple regression with a mixture-of-normals prior using the `rcpp_mr_ash_rss` function from the C++ implementation.
#'
#' @param bhat Numeric vector of observed effect sizes (standardized).
#' @param shat Numeric vector of standard errors of effect sizes.
#' @param z Numeric vector of Z-scores.
#' @param R Numeric matrix of the correlation matrix.
#' @param var_y Numeric value of the variance of the outcome.
#' @param n Integer value of the sample size.
#' @param sigma2_e Numeric value of the error variance.
#' @param s0 Numeric vector of prior variances for the mixture components.
#' @param w0 Numeric vector of prior weights for the mixture components.
#' @param mu1_init Numeric vector of initial values for the posterior mean of the coefficients.
#' @param tol Numeric value of the convergence tolerance. Default is 1e-8.
#' @param max_iter Integer value of the maximum number of iterations. Default is 1e5.
#' @param update_w0 Logical value indicating whether to update the mixture weights. Default is TRUE.
#' @param update_sigma Logical value indicating whether to update the error variance. Default is TRUE.
#' @param compute_ELBO Logical value indicating whether to compute the Evidence Lower Bound (ELBO). Default is TRUE.
#' @param standardize Logical value indicating whether to standardize the input data. Default is FALSE.
#' @param ncpu An integer specifying the number of CPU cores to use for parallel computation. Default is 1.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{mu1}{Numeric vector of the posterior mean of the coefficients.}
#'   \item{sigma2_1}{Numeric vector of the posterior variance of the coefficients.}
#'   \item{w1}{Numeric matrix of the posterior assignment probabilities.}
#'   \item{sigma2_e}{Numeric value of the error variance.}
#'   \item{w0}{Numeric vector of the mixture weights.}
#'   \item{ELBO}{Numeric value of the Evidence Lower Bound (if `compute_ELBO = TRUE`).}
#' }
#'
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
#' # Run mr_ash_rss
#' out <- mr_ash_rss(b.hat, s.hat,
#'   R = R.hat, var_y = var_y, n = n,
#'   sigma2_e = sigmasq_init, s0 = sigma0, w0 = omega0,
#'   mu1_init = rep(0, ncol(X)), tol = 1e-8, max_iter = 1e5,
#'   update_w0 = TRUE, update_sigma = TRUE, compute_ELBO = TRUE,
#'   standardize = FALSE
#' )
#' # In sample prediction correlations
#' cor(X %*% out1$mu1, y) # 0.9984064
#' @export
mr_ash_rss <- function(bhat, shat, R, var_y, n,
                       sigma2_e, s0, w0, mu1_init = numeric(0),
                       tol = 1e-8, max_iter = 1e5, z = numeric(0),
                       update_w0 = TRUE, update_sigma = TRUE,
                       compute_ELBO = TRUE, standardize = FALSE, ncpu = 1L) {
  # Check if ncpu is greater than 0 and is an integer
  if (ncpu <= 0 || !is.integer(ncpu)) {
    stop("ncpu must be a positive integer.")
  }

  if (is.null(var_y)) var_y <- Inf
  if (identical(z, numeric(0))) z <- bhat / shat # rcpp_mr_ass_rss throws error at line 269 of mr_ash.h
  result <- rcpp_mr_ash_rss(
    bhat = bhat, shat = shat, z = z, R = R,
    var_y = var_y, n = n, sigma2_e = sigma2_e,
    s0 = s0, w0 = w0, mu1_init = mu1_init,
    tol = tol, max_iter = max_iter,
    update_w0 = update_w0, update_sigma = update_sigma,
    compute_ELBO = compute_ELBO, standardize = standardize,
    ncpus = ncpu
  )

  return(result)
}

#' Extract weights from mr_ash_rss function
#' @return A numeric vector of the posterior mean of the coefficients.
#' @export
mr_ash_rss_weights <- function(stat, LD, var_y, sigma2_e, s0, w0, z = numeric(0), ...) {
  model <- mr_ash_rss(
    bhat = stat$b, shat = stat$seb, z = z, R = LD,
    var_y = var_y, n = median(stat$n), sigma2_e = sigma2_e,
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
#' @param n_iter Number of MCMC iterations. Default is 1000.
#' @param n_burnin Number of burn-in iterations. Default is 500.
#' @param thin Thinning factor for MCMC. Default is 5.
#' @param maf A vector of minor allele frequencies, if available, will standardize the effect sizes by MAF. Default is NULL.
#' @param verbose Whether to print verbose output. Default is FALSE.
#' @param seed Random seed for reproducibility. Default is NULL.
#'
#' @return A list containing the posterior estimates:
#'   - beta_est: Posterior estimates of SNP effect sizes.
#'   - psi_est: Posterior estimates of psi (shrinkage parameters).
#'   - sigma_est: Posterior estimate of the residual variance.
#'   - phi_est: Posterior estimate of the global shrinkage parameter.
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
#' out <- prs_cs(b.hat, LD, n, maf = maf)
#' # In sample prediction correlations
#' cor(X %*% out$beta_est, y) # 0.9944553
#' @export
prs_cs <- function(bhat, LD, n,
                   a = 1, b = 0.5, phi = NULL,
                   maf = NULL, n_iter = 1000, n_burnin = 500,
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
  total_rows_in_LD <- sum(sapply(LD, nrow))
  if (length(bhat) != total_rows_in_LD) {
    stop("The length of 'bhat' must be the same as the sum of the number of rows of all elements in the 'LD' list.")
  }

  # Run PRS-CS
  result <- prs_cs_rcpp(
    a = a, b = b, phi = phi, bhat, maf,
    n = n, ld_blk = LD,
    n_iter = n_iter, n_burnin = n_burnin, thin = thin,
    verbose = verbose, seed = seed
  )

  # Return the result as a list
  list(
    beta_est = result$beta_est,
    psi_est = result$psi_est,
    sigma_est = result$sigma_est,
    phi_est = result$phi_est
  )
}

#' Extract weights from prs_cs function
#' @return A numeric vector of the posterior SNP coefficients.
#' @export
prs_cs_weights <- function(stat, LD, ...) {
  model <- prs_cs(bhat = stat$b, LD = list(blk1 = LD), n = median(stat$n), ...)

  return(model$beta_est)
}

#' SDPR (Summary-Statistics-Based Dirichelt Process Regression for Polygenic Risk Prediction)
#'
#' This function is a wrapper for the SDPR C++ implementation, which performs Markov Chain Monte Carlo (MCMC)
#' for estimating effect sizes and heritability based on summary statistics and reference LD matrices.
#'
#' @param bhat A vector of marginal beta values for each SNP.
#' @param LD A list of LD matrices, where each matrix corresponds to a subset of SNPs.
#' @param n The total sample size of the GWAS.
#' @param per_variant_sample_size (Optional) A vector of sample sizes for each SNP. If NULL (default), it will be initialized
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
#' @param n_threads Number of threads to use. Default is 1.
#' @param opt_llk Which likelihood to evaluate. 1 for equation 6 (slightly shrink the correlation of SNPs)
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
#' cor(X %*% out$beta_est, y) #
#'
#' @note This function is a wrapper for the SDPR C++ implementation, which is a rewritten and adopted version
#'       of the SDPR package. The original SDPR documentation is available at
#'       https://htmlpreview.github.io/?https://github.com/eldronzhou/SDPR/blob/main/doc/Manual.html
#'
#' @export
sdpr <- function(bhat, LD, n, per_variant_sample_size = NULL, array = NULL, a = 0.1, c = 1.0, M = 1000,
                 a0k = 0.5, b0k = 0.5, iter = 1000, burn = 200, thin = 5, n_threads = 1,
                 opt_llk = 1, verbose = TRUE) {
  # Check if the sum of the rows in LD list is the same as length of bhat
  if (sum(sapply(LD, nrow)) != length(bhat)) {
    stop("The sum of the rows in LD list must be the same as the length of bhat.")
  }

  # Check if total sample size n is a positive integer
  if (missing(n) || n <= 0) {
    stop("The total sample size 'n' must be a positive integer.")
  }

  # Check if per_variant_sample_size vector contains only positive values (if provided)
  if (!is.null(per_variant_sample_size) && any(per_variant_sample_size <= 0)) {
    stop("The 'per_variant_sample_size' vector must contain only positive values.")
  }

  # Check if array vector contains only 0, 1, or 2 (if provided)
  if (!is.null(array) && any(!array %in% c(0, 1, 2))) {
    stop("The 'array' vector must contain only 0, 1, or 2.")
  }

  # Call the sdpr_rcpp function
  result <- sdpr_rcpp(
    bhat, LD, n, per_variant_sample_size, array, a, c, M, a0k, b0k, iter, burn, thin,
    n_threads, opt_llk, verbose
  )

  return(result)
}

#' Extract weights from sdpr function
#' @return A numeric vector of the posterior SNP coefficients.
#' @export
sdpr_weights <- function(stat, LD, ...) {
  model <- sdpr(bhat = stat$b, LD = list(blk1 = LD), n = median(stat$n), ...)

  return(model$beta_est)
}

#' @importFrom susieR coef.susie
#' @export
susie_weights <- function(X = NULL, y = NULL, susie_fit = NULL, ...) {
  if (is.null(susie_fit)) {
    # get susie_fit object
    susie_fit <- susie_wrapper(X, y, ...)
  }
  if (!is.null(X)) {
    if (length(susie_fit$pip) != ncol(X)) {
      stop(paste0(
        "Dimension mismatch on number of variant in susie_fit ", length(susie_fit$pip),
        " and TWAS weights ", ncol(X), ". "
      ))
    }
  }
  if ("alpha" %in% names(susie_fit) && "mu" %in% names(susie_fit) && "X_column_scale_factors" %in% names(susie_fit)) {
    # This is designed to cope with output from pecotmr::susie_post_processor()
    # We set intercept to 0 and later trim it off anyways
    susie_fit$intercept <- 0
    return(coef.susie(susie_fit)[-1])
  } else {
    return(rep(0, length(susie_fit$pip)))
  }
}

#' @importFrom mr.mash.alpha coef.mr.mash
#' @export
mrmash_weights <- function(mrmash_fit = NULL, X = NULL, Y = NULL, ...) {
  if (is.null(mrmash_fit)) {
    message("mrmash_fit is not provided; fitting mr.mash now ...")
    if (is.null(X) || is.null(Y)) {
      stop("Both X and Y must be provided if mrmash_fit is NULL.")
    }
    mrmash_fit <- mrmash_wrapper(X, Y, ...)
  }
  return(coef.mr.mash(mrmash_fit)[-1, ])
}

#' @importFrom mvsusieR mvsusie coef.mvsusie create_mixture_prior
#' @export
mvsusie_weights <- function(mvsusie_fit = NULL, X = NULL, Y = NULL, prior_variance = NULL, residual_variance = NULL, L = 30, ...) {
  if (is.null(mvsusie_fit)) {
    message("mvsusie_fit is not provided; fitting mvSuSiE now ...")
    if (is.null(X) || is.null(Y)) {
      stop("Both X and Y must be provided if mvsusie_fit is NULL.")
    }
    if (is.null(prior_variance)) prior_variance <- create_mixture_prior(R = ncol(Y))
    if (is.null(residual_variance)) residual_variance <- mr.mash.alpha:::compute_cov_flash(Y)

    mvsusie_fit <- mvsusie(
      X = X, Y = Y, L = L, prior_variance = prior_variance,
      residual_variance = residual_variance, precompute_covariances = F,
      compute_objective = T, estimate_residual_variance = F, estimate_prior_variance = T,
      estimate_prior_method = "EM", approximate = F, ...
    )
  }
  return(coef.mvsusie(mvsusie_fit)[-1, ])
}

# Get a reasonable setting for the standard deviations of the mixture
# components in the mixture-of-normals prior based on the data (X, y).
# Input se is an estimate of the residual *variance*, and n is the
# number of standard deviations to return. This code is adapted from
# the autoselect.mixsd function in the ashr package.
#' @importFrom susieR univariate_regression
init_prior_sd <- function(X, y, n = 30) {
  res <- univariate_regression(X, y)
  smax <- 3 * max(res$betahat)
  seq(0, smax, length.out = n)
}

#' @importFrom stats coef
#' @export
glmnet_weights <- function(X, y, alpha) {
  # Check if glmnet is installed
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("To use this function, please install glmnet: https://cran.r-project.org/web/packages/glmnet/index.html")
  }
  eff.wgt <- matrix(0, ncol = 1, nrow = ncol(X))
  sds <- apply(X, 2, sd)
  keep <- sds != 0 & !is.na(sds)
  enet <- glmnet::cv.glmnet(x = X[, keep], y = y, alpha = alpha, nfold = 5, intercept = T, standardize = F)
  eff.wgt[keep] <- coef(enet, s = "lambda.min")[2:(sum(keep) + 1)]
  return(eff.wgt)
}

#' @export
enet_weights <- function(X, y) glmnet_weights(X, y, 0.5)

#' @export
lasso_weights <- function(X, y) glmnet_weights(X, y, 1)

#' @examples
#' wgt.mr.ash <- mrash_weights(eqtl$X, eqtl$y_res, beta.init = lasso_weights(X, y))
#' @importFrom stats predict
#' @export
mrash_weights <- function(X, y, init_prior_sd = TRUE, ...) {
  # Make sure mr.ash is installed
  if (!requireNamespace("mr.ash.alpha", quietly = TRUE)) {
    stop("To use this function, please install mr.ash: https://github.com/stephenslab/mr.ash.alpha")
  }
  args_list <- list(...)
  if (!"beta.init" %in% names(args_list)) {
    args_list$beta.init <- lasso_weights(X, y)
  }
  fit.mr.ash <- do.call(getFromNamespace("mr.ash", "mr.ash.alpha"), c(list(X = X, y = y, sa2 = if (init_prior_sd) init_prior_sd(X, y)^2 else NULL), args_list))
  predict(fit.mr.ash, type = "coefficients")[-1]
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
#' bayes_l_weights(y = y, X = X, Z = Z)
#' bayes_r_weights(y = y, X = X, Z = Z)
#' @export
bayes_alphabet_weights <- function(X, y, method, Z = NULL, nit = 5000, nburn = 1000, nthin = 5, ...) {
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

  model <- qgg::gbayes(
    y = y,
    W = X,
    X = Z,
    method = method,
    nit = nit,
    nburn = nburn,
    ...
  )

  return(model$bm)
}
#' Use Gaussian distribution as prior. Posterior means will be BLUP, equivalent to Ridge Regression.
#' @export
bayes_n_weights <- function(X, y, Z = NULL, ...) {
  return(bayes_alphabet_weights(X, y, method = "bayesN", Z, ...))
}
#' Use laplace/double exponential distribution as prior. This is equivalent to Bayesian LASSO.
#' @export
bayes_l_weights <- function(X, y, Z = NULL, ...) {
  return(bayes_alphabet_weights(X, y, method = "bayesL", Z, ...))
}
#' Use t-distribution as prior.
#' @export
bayes_a_weights <- function(X, y, Z = NULL, ...) {
  return(bayes_alphabet_weights(X, y, method = "bayesA", Z, ...))
}
#' Use a rounded spike prior (low-variance Gaussian).
#' @export
bayes_c_weights <- function(X, y, Z = NULL, ...) {
  return(bayes_alphabet_weights(X, y, method = "bayesC", Z, ...))
}
#' Use a hierarchical Bayesian mixture model with four Gaussian components. Variances are scaled
#' by 0, 0.0001 , 0.001 , and 0.01 .
#' @export
bayes_r_weights <- function(X, y, Z = NULL, ...) {
  return(bayes_alphabet_weights(X, y, method = "bayesR", Z, ...))
}


#' Bayesian linear regression using summary statistics
#'
#' @description
#'
#' This function is adapted from those written by Peter Sørensen in the qgg package.
#' The following prior distributions are provided:
#'
#' Bayes N: Assigning a Gaussian prior to marker effects implies that the posterior means are the
#' BLUP estimates (same as Ridge Regression).
#'
#' Bayes L: Assigning a double-exponential or Laplace prior is the density used in
#' the Bayesian LASSO
#'
#' Bayes A: similar to ridge regression but t-distribution prior (rather than Gaussian)
#' for the marker effects ; variance comes from an inverse-chi-square distribution instead of being fixed. Estimation
#' via Gibbs sampling.
#'
#' Bayes C: uses a “rounded spike” (low-variance Gaussian) at origin many small
#' effects can contribute to polygenic component, reduces the dimensionality of
#' the model (makes Gibbs sampling feasible).
#'
#' Bayes R: Hierarchical Bayesian mixture model with 4 Gaussian components, with
#' variances scaled by 0, 0.0001 , 0.001 , and 0.01 .
#'
#' @param sumstats dataframe with marker summary statistics. Required: beta coefficient (beta), standard
#'        error of the beta coefficient (se), GWAS sample size (n). Optional: variant_id or rsid, alleles (A1
#'        and A2), minor allele frequency (maf).
#' @param LD is a the LD matrix corresponding to the same markers as in the stat dataframe
#' @param variant_ids is an optional character vector of variant ids or rsids, provided outside of the rss dataframe
#' @param nit is the number of iterations
#' @param nburn is the number of burnin iterations
#' @param nthin is the thinning parameter
#' @param method specifies the methods used (method="bayesN","bayesA","bayesL","bayesC","bayesR")
#' @param vg is a scalar or matrix of genetic (co)variances
#' @param vb is a scalar or matrix of marker (co)variances
#' @param ve is a scalar or matrix of residual (co)variances
#' @param ssg_prior is a scalar or matrix of prior genetic (co)variances
#' @param ssb_prior is a scalar or matrix of prior marker (co)variances
#' @param sse_prior is a scalar or matrix of prior residual (co)variances
#' @param lambda is a vector or matrix of lambda values
#' @param h2 is the trait heritability
#' @param pi is the proportion of markers in each marker variance class
#' @param updateB is a logical for updating marker (co)variances
#' @param updateG is a logical for updating genetic (co)variances
#' @param updateE is a logical for updating residual (co)variances
#' @param updatePi is a logical for updating pi
#' @param adjustE is a logical for adjusting residual variance
#' @param nug is a scalar or vector of prior degrees of freedom for prior genetic (co)variances
#' @param nub is a scalar or vector of prior degrees of freedom for marker (co)variances
#' @param nue is a scalar or vector of prior degrees of freedom for prior residual (co)variances
#' @param mask is a vector or matrix of TRUE/FALSE specifying if marker should be ignored
#' @param ve_prior is a scalar or matrix of prior residual (co)variances
#' @param vg_prior is a scalar or matrix of prior genetic (co)variances
#' @param algorithm is the algorithm to use. Should take on values ("mcmc", "em-mcmc")
#' @param tol is tolerance, i.e. convergence criteria used in gbayes
#' @param nit_local is the number of local iterations
#' @param nit_global is the number of global iterations
#'
#' @return Returns a list structure including
#' \item{bm}{vector of posterior means for marker effects}
#' \item{dm}{vector of posterior means for marker inclusion probabilities}
#' \item{vbs}{scalar or vector (t) of posterior means for marker variances}
#' \item{vgs}{scalar or vector (t) of posterior means for genomic variances}
#' \item{ves}{scalar or vector (t) of posterior means for residual variances}
#' \item{pis}{vector of probabilites for each mcmc iteration}
#' \item{pim}{posterior distribution probabilities}
#' \item{r}{vector of residuals}
#' \item{b}{vector of estimates from the final mcmc iteration}
#' \item{param}{a list current parameters (same information as item listed above)
#'              used for restart of the analysis}
#' \item{stat}{matrix (mxt) of marker information and effects used for genomic risk scoring}
#' \item{method}{the method used}
#' \item{mask}{which loci were masked from analysis}
#' \item{conv}{dataframe of convergence metrics}
#' \item{post}{posterior parameter estimates}
#' \item{ve}{mean residual variance}
#' \item{vg}{mean genomic variance}
#'
#' @export
gbayes_rss <- function(sumstats = NULL, LD = NULL, variant_ids = NULL, nit = 100, nburn = 0, nthin = 4, method = "bayesR",
                       vg = NULL, vb = NULL, ve = NULL, ssg_prior = NULL, ssb_prior = NULL, sse_prior = NULL,
                       lambda = NULL, h2 = NULL, pi = 0.001, updateB = TRUE, updateG = TRUE, updateE = TRUE,
                       updatePi = TRUE, adjustE = TRUE, nug = 4, nub = 4, nue = 4, mask = NULL, ve_prior = NULL,
                       vg_prior = NULL, algorithm = "mcmc", tol = 0.001, nit_local = NULL, nit_global = NULL) {
  # Make sure qgg is installed
  if (!requireNamespace("qgg", quietly = TRUE)) {
    stop("To use this function, please install qgg: https://cran.r-project.org/web/packages/qgg/index.html")
  }
  # Check methods
  methods <- c("bayesN", "bayesA", "bayesL", "bayesC", "bayesR")
  method <- match(method, methods)
  if (!sum(method %in% c(1:5)) == 1) stop("Method specified not valid")
  if (method == 0) {
    # BLUP and we do not estimate parameters
    updateB <- FALSE
    updateE <- FALSE
  }

  # Set algorithm
  if (algorithm == "em-mcmc") {
    algo <- 2
  } else {
    algo <- 1
  }

  # Check that LD matrix is provided and of same length as stats
  if (is.null(LD)) stop("Must provide LD matrix")
  if (nrow(sumstats) != nrow(LD)) stop("LD matrix must correspond to summary statistics")

  # Parameters from stat df
  if (is.data.frame(sumstats)) {
    if (!is.null(variant_ids)) {
      variant_ids <- variant_ids
    } else if (!is.null(sumstats$rsids)) {
      variant_ids <- sumstats$rsids
    } else if (!is.null(sumstats$variant_id)) {
      variant_ids <- sumstats$variant_id
    } else {
      variant_ids <- paste0("snp", 1:nrow(sumstats))
      sumstats$variant_id <- variant_ids
    }

    m <- length(variant_ids)
    b <- wy <- ww <- matrix(0, nrow = nrow(sumstats), ncol = 1)
    mask <- matrix(FALSE, nrow = nrow(sumstats), ncol = 1)
    rownames(b) <- rownames(wy) <- rownames(ww) <- rownames(mask) <- variant_ids

    if (is.null(sumstats$ww)) sumstats$ww <- 1 / (sumstats$se^2 + sumstats$beta^2 / sumstats$n)
    if (is.null(sumstats$wy)) sumstats$wy <- sumstats$beta * sumstats$ww
    if (!is.null(sumstats$n)) n <- as.integer(median(sumstats$n))
    ww[, 1] <- sumstats$ww
    wy[, 1] <- sumstats$wy
    mask[, 1] <- FALSE

    if (any(is.na(wy))) stop("Missing values in wy")
    if (any(is.na(ww))) stop("Missing values in ww")

    b2 <- sumstats$beta^2
    seb2 <- sumstats$se^2
    yy <- (b2 + (n - 2) * seb2) * sumstats$ww
    yy <- median(yy)

    if (is.null(sumstats$A1)) sumstats$A1 <- rep("Unknown", length = nrow(sumstats))
    if (is.null(sumstats$A2)) sumstats$A2 <- rep("Unknown", length = nrow(sumstats))
    if (is.null(sumstats$maf)) {
      sumstats$maf <- rep("Unknown", length = nrow(sumstats))
      af_prov <- 0
    } else {
      af_prov <- 1
    }
  } else {
    stop("Summary statistics must be provided in dataframe")
  }


  # prep LD for gbayes
  LD_values <- lapply(1:nrow(LD), function(i) as.numeric(LD[i, ]))
  names(LD_values) <- variant_ids

  LD_indices <- list(indices = vector("list", length = nrow(LD)))
  for (i in 1:nrow(LD)) {
    LD_indices[[i]] <- 1:nrow(LD) - 1
  }

  bm <- dm <- fit <- res <- vector(length = 1, mode = "list")
  names(bm) <- names(dm) <- names(fit) <- names(res) <- 1

  # Set parameters if not otherwise specified
  if (is.null(m)) m <- length(LD_values)
  vy <- yy / (n - 1)
  if (is.null(pi)) pi <- 0.001
  if (is.null(h2)) h2 <- 0.5
  if (is.null(ve)) ve <- vy * (1 - h2)
  if (is.null(vg)) vg <- vy * h2
  if (method < 4 && is.null(vb)) vb <- vg / m
  if (method >= 4 && is.null(vb)) vb <- vg / (m * pi)
  if (is.null(lambda)) lambda <- rep(ve / vb, m)
  if (method < 4 && is.null(ssb_prior)) ssb_prior <- ((nub - 2.0) / nub) * (vg / m)
  if (method >= 4 && is.null(ssb_prior)) ssb_prior <- ((nub - 2.0) / nub) * (vg / (m * pi))
  if (is.null(sse_prior)) sse_prior <- ((nue - 2.0) / nue) * ve
  if (is.null(b)) b <- rep(0, m)

  pi <- c(1 - pi, pi)
  gamma <- c(0, 1.0)
  if (method == 5) pi <- c(0.95, 0.02, 0.02, 0.01)
  if (method == 5) gamma <- c(0, 0.01, 0.1, 1.0)

  seed <- sample.int(.Machine$integer.max, 1)

  fit <- qgg:::sbayes_spa(
    wy = wy,
    ww = ww,
    LDvalues = LD_values,
    LDindices = LD_indices,
    b = b,
    lambda = lambda,
    mask = mask,
    yy = yy,
    pi = pi,
    gamma = gamma,
    vg = vg,
    vb = vb,
    ve = ve,
    ssb_prior = ssb_prior,
    sse_prior = sse_prior,
    nub = nub,
    nue = nue,
    updateB = updateB,
    updateE = updateE,
    updatePi = updatePi,
    updateG = updateG,
    adjustE = adjustE,
    n = n,
    nit = nit,
    nburn = nburn,
    nthin = nthin,
    algo = algo,
    method = as.integer(method),
    seed = seed
  )

  names(fit[[1]]) <- names(LD_values)
  names(fit) <- c("bm", "dm", "coef", "vbs", "vgs", "ves", "pis", "pim", "r", "b", "param")
  fit[3] <- NULL

  res <- data.frame(
    variant_id = variant_ids, bm = fit$bm, dm = fit$dm,
    pos = sumstats$pos, A1 = sumstats$A1,
    A2 = sumstats$A2, maf = sumstats$maf,
    stringsAsFactors = FALSE
  )
  rownames(res) <- variant_ids

  fit$sumstats <- res
  if (af_prov == 1) {
    fit$sumstats$vm <- 2 * (1 - fit$sumstats$maf) * fit$sumstats$maf * fit$sumstats$bm^2
  }
  fit$method <- methods[method]
  fit$mask <- mask

  zve <- coda::geweke.diag(fit$ves[nburn:length(fit$ves)])$z
  zvg <- coda::geweke.diag(fit$vgs[nburn:length(fit$vgs)])$z
  zvb <- coda::geweke.diag(fit$vbs[nburn:length(fit$vbs)])$z
  zpi <- coda::geweke.diag(fit$pis[nburn:length(fit$pis)])$z

  ve <- mean(fit$ves[nburn:length(fit$ves)])
  vg <- mean(fit$vgs[nburn:length(fit$vgs)])
  vb <- mean(fit$vbs[nburn:length(fit$vbs)])
  pi <- 1 - fit$pim[1]
  fit$conv <- data.frame(zve = zve, zvg = zvg, zvb = zvb, zpi = zpi)
  fit$post <- data.frame(ve = ve, vg = vg, vb = vb, pi = pi)
  fit$ve <- mean(ve)
  fit$vg <- sum(vg)

  return(fit)
}
#' Extract weights from gbayes_rss function
#' @return A numeric vector of the posterior mean of the coefficients.
#' @export
bayes_alphabet_rss_weights <- function(sumstats, LD, method, ...) {
  model <- gbayes_rss(sumstats = sumstats, LD = LD, method = method, ...)
  return(model$bm)
}
#' Use Gaussian distribution as prior. Posterior means will be BLUP, equivalent to Ridge Regression.
#' @export
bayes_n_rss_weights <- function(sumstats, LD, ...) {
  return(bayes_alphabet_rss_weights(sumstats, LD, method = "bayesN", ...))
}
#' Use laplace/double exponential distribution as prior. This is equivalent to Bayesian LASSO.
#' @export
bayes_l_rss_weights <- function(sumstats, LD, ...) {
  return(bayes_alphabet_rss_weights(sumstats, LD, method = "bayesL", ...))
}
#' Use t-distribution as prior.
#' @export
bayes_a_rss_weights <- function(sumstats, LD, ...) {
  return(bayes_alphabet_rss_weights(sumstats, LD, method = "bayesA", ...))
}
#' Use a rounded spike prior (low-variance Gaussian).
#' @export
bayes_c_rss_weights <- function(sumstats, LD, ...) {
  return(bayes_alphabet_rss_weights(sumstats, LD, method = "bayesC", ...))
}
#' Use a hierarchical Bayesian mixture model with four Gaussian components. Variances are scaled
#' by 0, 0.0001 , 0.001 , and 0.01 .
#' @export
bayes_r_rss_weights <- function(sumstats, LD, ...) {
  return(bayes_alphabet_rss_weights(sumstats, LD, method = "bayesR", ...))
}

#' @export
susie_ash_weights <- function(susie_ash_fit, X = NULL, y = NULL, ...) {
  # If the fit object is missing or NULL, try to recover it from the parent frame.
  if (missing(susie_ash_fit) || is.null(susie_ash_fit)) {
    susie_ash_fit <- get0("susie_ash_fit", envir = parent.frame())
    if (is.null(susie_ash_fit)) {
      stop("A susie_ash_fit object is required.")
    }
  }
  if (!is.null(X)) {
    if (length(susie_ash_fit$marginal_PIP) != ncol(X)) {
      stop(paste0("Dimension mismatch on number of variants in susie_ash_fit ",
                  length(susie_ash_fit$marginal_PIP), " and TWAS weights ", ncol(X), ". "))
    }
  }
  # Calculate coefficients as per the provided formula.
  weights <- rowSums(susie_ash_fit$mu * susie_ash_fit$PIP) + susie_ash_fit$theta
  return(weights)
}

#' @export
susie_inf_weights <- function(susie_inf_fit, X = NULL, y = NULL, ...) {
  # If the fit object is missing or NULL, try to recover it from the parent frame.
  if (missing(susie_inf_fit) || is.null(susie_inf_fit)) {
    susie_inf_fit <- get0("susie_inf_fit", envir = parent.frame())
    if (is.null(susie_inf_fit)) {
      stop("A susie_inf_fit object is required.")
    }
  }
  if (!is.null(X)) {
    if (length(susie_inf_fit$marginal_PIP) != ncol(X)) {
      stop(paste0("Dimension mismatch on number of variants in susie_inf_fit ",
                  length(susie_inf_fit$marginal_PIP), " and TWAS weights ", ncol(X), ". "))
    }
  }
  # Calculate coefficients as per the provided formula.
  weights <- rowSums(susie_inf_fit$mu * susie_inf_fit$PIP) + susie_inf_fit$alpha
  return(weights)
}
