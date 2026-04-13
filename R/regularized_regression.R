#' Extract weights from mr.ash.rss (susieR)
#' @return A numeric vector of the posterior mean of the coefficients.
#' @importFrom susieR mr.ash.rss
#' @export
mr_ash_rss_weights <- function(stat, LD, var_y, sigma2_e, s0, w0, z = numeric(0), ...) {
  model <- mr.ash.rss(
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
                 opt_llk = 1, verbose = TRUE, seed = NULL) {
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
    n_threads, opt_llk, verbose, seed
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

#' @importFrom susieR coef.susie
#' @export
susie_ash_weights <- function(X = NULL, y = NULL, susie_ash_fit = NULL, ...) {
  if (is.null(susie_ash_fit)) {
    # get susie_ash_fit object
    susie_ash_fit <- susie_wrapper(X, y, unmappable_effects = "ash", convergence_method = "pip", ...)
  }
  if (!is.null(X)) {
    if (length(susie_ash_fit$pip) != ncol(X)) {
      stop(paste0(
        "Dimension mismatch on number of variant in susie_ash_fit ", length(susie_ash_fit$pip),
        " and TWAS weights ", ncol(X), ". "
      ))
    }
  }
  if ("alpha" %in% names(susie_ash_fit) && "mu" %in% names(susie_ash_fit) && "theta" %in% names(susie_ash_fit) && "X_column_scale_factors" %in% names(susie_ash_fit)) {
    # This is designed to cope with output from pecotmr::susie_post_processor()
    # We set intercept to 0 and later trim it off anyways
    susie_ash_fit$intercept <- 0
    return(coef.susie(susie_ash_fit)[-1])
  } else {
    return(rep(0, length(susie_ash_fit$pip)))
  }
}

#' @importFrom susieR coef.susie
#' @export
susie_inf_weights <- function(X = NULL, y = NULL, susie_inf_fit = NULL, ...) {
  if (is.null(susie_inf_fit)) {
    # get susie_inf_fit object
    susie_inf_fit <- susie_wrapper(X, y, unmappable_effects = "inf", convergence_method = "pip", ...)
  }
  if (!is.null(X)) {
    if (length(susie_inf_fit$pip) != ncol(X)) {
      stop(paste0(
        "Dimension mismatch on number of variant in susie_inf_fit ", length(susie_inf_fit$pip),
        " and TWAS weights ", ncol(X), ". "
      ))
    }
  }
  if ("alpha" %in% names(susie_inf_fit) && "mu" %in% names(susie_inf_fit) && "theta" %in% names(susie_inf_fit) && "X_column_scale_factors" %in% names(susie_inf_fit)) {
    # This is designed to cope with output from pecotmr::susie_post_processor()
    # We set intercept to 0 and later trim it off anyways
    susie_inf_fit$intercept <- 0
    return(coef.susie(susie_inf_fit)[-1])
  } else {
    return(rep(0, length(susie_inf_fit$pip)))
  }
}

#' @export
mrmash_weights <- function(mrmash_fit = NULL, X = NULL, Y = NULL, ...) {
  if (!requireNamespace("mr.mashr", quietly = TRUE)) {
    stop("Package 'mr.mashr' is required. Install with: devtools::install_github('stephenslab/mr.mashr')")
  }
  if (is.null(mrmash_fit)) {
    message("mrmash_fit is not provided; fitting mr.mash now ...")
    if (is.null(X) || is.null(Y)) {
      stop("Both X and Y must be provided if mrmash_fit is NULL.")
    }
    mrmash_fit <- mrmash_wrapper(X, Y, ...)
  }
  return(mr.mashr::coef.mr.mash(mrmash_fit)[-1, ])
}

#' @export
mvsusie_weights <- function(mvsusie_fit = NULL, X = NULL, Y = NULL, prior_variance = NULL, residual_variance = NULL, L = 30, ...) {
  if (!requireNamespace("mvsusieR", quietly = TRUE)) {
    stop("Package 'mvsusieR' is required. Install with: devtools::install_github('stephenslab/mvsusieR')")
  }
  if (is.null(mvsusie_fit)) {
    message("mvsusie_fit is not provided; fitting mvSuSiE now ...")
    if (is.null(X) || is.null(Y)) {
      stop("Both X and Y must be provided if mvsusie_fit is NULL.")
    }
    if (is.null(prior_variance)) prior_variance <- mvsusieR::create_mixture_prior(R = ncol(Y))
    if (is.null(residual_variance)) {
      if (!requireNamespace("mr.mashr", quietly = TRUE)) {
        stop("Package 'mr.mashr' is required for residual variance estimation. Install with: devtools::install_github('stephenslab/mr.mashr')")
      }
      residual_variance <- mr.mashr:::compute_cov_flash(Y)
    }

    mvsusie_fit <- mvsusieR::mvsusie(
      X = X, Y = Y, L = L, prior_variance = prior_variance,
      residual_variance = residual_variance, precompute_covariances = FALSE,
      compute_objective = TRUE, estimate_residual_variance = FALSE, estimate_prior_variance = TRUE,
      estimate_prior_method = "EM", approximate = FALSE, ...
    )
  }
  return(mvsusieR::coef.mvsusie(mvsusie_fit)[-1, ])
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

# Identify zero-variance columns of X and warn the caller before they are
# dropped. Returns a logical vector of length ncol(X) where TRUE indicates a
# column to keep. The downstream solvers in pecotmr's regression wrappers
# (glmnet, ncvreg, L0Learn, qgg, BGLR, RcppDPR) all either error or behave
# poorly on constant columns, so wrappers should filter them out and zero-pad
# their results back to length p.
.drop_zero_variance <- function(X, fn_name) {
  sds <- apply(X, 2, sd)
  keep <- !is.na(sds) & sds != 0
  if (!all(keep)) {
    warning(sprintf(
      "%s: dropping %d zero-variance column(s) from X (indices: %s)",
      fn_name, sum(!keep),
      paste(which(!keep), collapse = ", ")
    ), call. = FALSE)
  }
  keep
}

#' @importFrom stats coef
#' @export
glmnet_weights <- function(X, y, alpha) {
  # Check if glmnet is installed
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("To use this function, please install glmnet: https://cran.r-project.org/web/packages/glmnet/index.html")
  }
  eff.wgt <- matrix(0, ncol = 1, nrow = ncol(X))
  keep <- .drop_zero_variance(X, "glmnet_weights")
  enet <- glmnet::cv.glmnet(x = X[, keep, drop = FALSE], y = y, alpha = alpha, nfold = 5, intercept = TRUE, standardize = FALSE)
  eff.wgt[keep] <- coef(enet, s = "lambda.min")[2:(sum(keep) + 1)]
  return(eff.wgt)
}

#' @export
enet_weights <- function(X, y) glmnet_weights(X, y, 0.5)

#' @export
lasso_weights <- function(X, y) glmnet_weights(X, y, 1)

#' Compute Weights Using mr.ash Shrinkage
#'
#' This function fits the `mr.ash` model (adaptive shrinkage regression) to estimate weights
#' for a given set of predictors and response. It uses optional prior standard deviation initialization
#' and can accept custom initial beta values.
#'
#' @examples
#' wgt.mr.ash <- mrash_weights(eqtl$X, eqtl$y_res, beta.init = lasso_weights(X, y))
#' @importFrom susieR mr.ash
#' @importFrom stats predict
#' @export
mrash_weights <- function(X, y, init_prior_sd = TRUE, ...) {
  eff.wgt <- rep(0, ncol(X))
  keep <- .drop_zero_variance(X, "mrash_weights")
  X_keep <- X[, keep, drop = FALSE]
  args_list <- list(...)
  if (!"beta.init" %in% names(args_list)) {
    args_list$beta.init <- lasso_weights(X_keep, y)
  } else if (length(args_list$beta.init) == ncol(X)) {
    args_list$beta.init <- args_list$beta.init[keep]
  }
  fit.mr.ash <- do.call(mr.ash, c(list(X = X_keep, y = y, sa2 = if (init_prior_sd) init_prior_sd(X_keep, y)^2 else NULL), args_list))
  eff.wgt[keep] <- predict(fit.mr.ash, type = "coefficients")[-1]
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

  eff.wgt <- rep(0, ncol(X))
  keep <- .drop_zero_variance(X, "bayes_alphabet_weights")

  model <- qgg::gbayes(
    y = y,
    W = X[, keep, drop = FALSE],
    X = Z,
    method = method,
    nit = nit,
    nburn = nburn,
    ...
  )

  eff.wgt[keep] <- model$bm
  return(eff.wgt)
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
#'   \item{beta_est}{Posterior estimates of SNP effect sizes at best lambda.}
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
#' out <- lassosum_rss(bhat, LD, n)
#' @export
lassosum_rss <- function(bhat, LD, n,
                         lambda = exp(seq(log(0.0001), log(0.1), length.out = 20)),
                         thr = 1e-4, maxiter = 10000) {
  if (!is.list(LD)) {
    stop("Please provide a valid list of LD blocks using 'LD'.")
  }
  if (missing(n) || n <= 0) {
    stop("Please provide a valid sample size using 'n'.")
  }
  total_rows_in_LD <- sum(sapply(LD, nrow))
  if (length(bhat) != total_rows_in_LD) {
    stop("The length of 'bhat' must be the same as the sum of the number of rows of all elements in the 'LD' list.")
  }

  z <- bhat / sqrt(n)
  order <- order(lambda, decreasing = TRUE)
  result <- lassosum_rss_rcpp(z, LD, lambda[order], thr, maxiter)

  # Reorder back to original lambda order.
  # Must use inverse permutation to unsort: if order[i]=j, then
  # the result at position j in the sorted output goes to position i.
  inv_order <- order(order)
  result$beta  <- result$beta[, inv_order, drop = FALSE]
  result$conv  <- result$conv[inv_order]
  result$loss  <- result$loss[inv_order]
  result$fbeta <- result$fbeta[inv_order]
  result$lambda <- lambda
  result$nparams <- as.integer(colSums(result$beta != 0))
  result$beta_est <- as.numeric(result$beta[, which.min(result$fbeta)])
  result
}

#' Extract weights from lassosum_rss with shrinkage grid search
#'
#' Searches over a grid of shrinkage parameters \code{s} (default:
#' \code{c(0.2, 0.5, 0.9, 1.0)}, matching the original lassosum and OTTERS).
#' For each \code{s}, the LD matrix is shrunk as \code{(1-s)*R + s*I}, then
#' \code{lassosum_rss()} is called across the lambda path. The best
#' \code{(s, lambda)} combination is selected by the lowest objective value.
#'
#' @details
#' Model selection uses \code{min(fbeta)} (penalized objective) rather than
#' the pseudovalidation approach from the original lassosum R package. Empirical
#' comparison over 20 random trials (n=300, p=50, 3 causal) shows no systematic
#' advantage for either method: pseudovalidation won 4/20, min(fbeta) won 6/20,
#' tied 10/20. The shrinkage grid over \code{s} provides the primary regularization;
#' lambda selection within each \code{s} has minimal impact.
#'
#' @param stat A list with \code{$b} (effect sizes) and \code{$n} (per-variant sample sizes).
#' @param LD LD correlation matrix R (single matrix, NOT pre-shrunk).
#' @param s Numeric vector of shrinkage parameters to search over. Default:
#'   \code{c(0.2, 0.5, 0.9, 1.0)} following Mak et al (2017) and OTTERS.
#' @param ... Additional arguments passed to \code{lassosum_rss()}.
#'
#' @return A numeric vector of the posterior SNP coefficients at the best (s, lambda).
#' @export
lassosum_rss_weights <- function(stat, LD, s = c(0.2, 0.5, 0.9, 1.0), ...) {
  n <- median(stat$n)
  p <- nrow(LD)
  best_fbeta <- Inf
  best_beta  <- rep(0, p)

  # Clamp marginal correlations to (-1, 1) as required by lassosum.
  # This is lassosum-specific — other methods (PRS-CS, SDPR) handle
  # their own regularization and should not be globally rescaled.
  # Matches OTTERS shrink_factor logic (PRSmodels/lassosum.R lines 71-77).
  bhat <- stat$b
  max_abs_b <- max(abs(bhat))
  if (max_abs_b >= 1) {
    bhat <- bhat / (max_abs_b / 0.9999)
  }

  for (s_val in s) {
    # Shrink LD: R_s = (1 - s) * R + s * I
    LD_s <- (1 - s_val) * LD + s_val * diag(p)
    model <- lassosum_rss(bhat = bhat, LD = list(blk1 = LD_s), n = n, ...)
    min_fbeta <- min(model$fbeta)
    if (min_fbeta < best_fbeta) {
      best_fbeta <- min_fbeta
      best_beta  <- model$beta_est
    }
  }

  return(best_beta)
}

#' Compute Weights Using ncvreg with SCAD or MCP Penalty
#'
#' Internal helper that fits an `ncvreg` model with the specified non-convex
#' penalty using k-fold cross-validation, then returns the coefficients at
#' `lambda.min`. Following the convention of `glmnet_weights`, columns of `X`
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
ncvreg_weights <- function(X, y, penalty, nfolds = 5, ...) {
  if (!requireNamespace("ncvreg", quietly = TRUE)) {
    stop("To use this function, please install ncvreg: https://cran.r-project.org/package=ncvreg")
  }
  eff.wgt <- matrix(0, ncol = 1, nrow = ncol(X))
  keep <- .drop_zero_variance(X, "ncvreg_weights")
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
scad_weights <- function(X, y, nfolds = 5, ...) {
  ncvreg_weights(X, y, penalty = "SCAD", nfolds = nfolds, ...)
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
mcp_weights <- function(X, y, nfolds = 5, ...) {
  ncvreg_weights(X, y, penalty = "MCP", nfolds = nfolds, ...)
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
l0learn_weights <- function(X, y, penalty = "L0", nFolds = 5, ...) {
  if (!requireNamespace("L0Learn", quietly = TRUE)) {
    stop("To use this function, please install L0Learn: https://cran.r-project.org/package=L0Learn")
  }
  eff.wgt <- matrix(0, ncol = 1, nrow = ncol(X))
  keep <- .drop_zero_variance(X, "l0learn_weights")
  fit <- L0Learn::L0Learn.cvfit(
    x = X[, keep, drop = FALSE], y = y, penalty = penalty, nFolds = nFolds, ...
  )
  # Find (gamma, lambda) minimizing CV error across the entire path.
  cv_mins <- vapply(fit$cvMeans, function(v) min(as.numeric(v)), numeric(1))
  gamma_idx <- which.min(cv_mins)
  lambda_idx <- which.min(as.numeric(fit$cvMeans[[gamma_idx]]))
  best_gamma <- fit$fit$gamma[gamma_idx]
  best_lambda <- fit$fit$lambda[[gamma_idx]][lambda_idx]
  coefs <- as.numeric(coef(fit, lambda = best_lambda, gamma = best_gamma))
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
#' @param eta_args Optional named list of additional arguments included in the
#'   `ETA` linear-term specification (e.g. `list(probIn = 0.05)` for BayesB).
#' @param ... Additional arguments passed through to `BGLR::BGLR`.
#' @return A numeric vector of length `ncol(X)` of variant weights.
#' @keywords internal
bglr_weights <- function(X, y, model, nIter, burnIn, thin, eta_args = list(), ...) {
  if (!requireNamespace("BGLR", quietly = TRUE)) {
    stop("To use this function, please install BGLR: https://cran.r-project.org/package=BGLR")
  }
  eff.wgt <- rep(0, ncol(X))
  keep <- .drop_zero_variance(X, "bglr_weights")

  tmpdir <- tempfile("bglr_")
  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
  saveAt <- paste0(tmpdir, .Platform$file.sep)

  eta <- list(c(list(X = X[, keep, drop = FALSE], model = model), eta_args))
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
#' @param probIn Prior inclusion probability for each marker. Default is 0.05.
#' @param ... Additional arguments passed through to `BGLR::BGLR`.
#' @return A numeric vector of length `ncol(X)` of variant weights.
#' @export
bayes_b_weights <- function(X, y, nIter = 10000, burnIn = 2000, thin = 5, probIn = 0.05, ...) {
  bglr_weights(
    X, y,
    model = "BayesB", nIter = nIter, burnIn = burnIn, thin = thin,
    eta_args = list(probIn = probIn), ...
  )
}

#' Compute Weights Using the Bayesian LASSO (BGLR)
#'
#' Fits a Bayesian LASSO linear regression model via `BGLR::BGLR` (the "BL"
#' model, Park & Casella 2008) and returns the posterior mean of the marker
#' effects. This is the same "B-Lasso" implementation benchmarked in Kim et
#' al. (2022). Note that this is distinct from `bayes_l_weights`, which uses a
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
b_lasso_weights <- function(X, y, nIter = 10000, burnIn = 2000, thin = 5, ...) {
  bglr_weights(
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
#' @param fitting_method One of "VB", "Gibbs", or "Adaptive_Gibbs". Default is "VB".
#' @param ... Additional arguments passed through to `RcppDPR::fit_model`.
#' @return A numeric vector of length `ncol(X)` of variant weights.
#' @export
dpr_weights <- function(X, y, fitting_method = "VB", ...) {
  if (!requireNamespace("RcppDPR", quietly = TRUE)) {
    stop("To use this function, please install RcppDPR: https://cran.r-project.org/package=RcppDPR")
  }
  eff.wgt <- rep(0, ncol(X))
  keep <- .drop_zero_variance(X, "dpr_weights")
  w <- matrix(1, nrow = nrow(X), ncol = 1)
  fit <- RcppDPR::fit_model(
    y = y, w = w, x = X[, keep, drop = FALSE],
    rotate_variables = FALSE, fitting_method = fitting_method, ...
  )
  eff.wgt[keep] <- as.numeric(fit$beta + fit$alpha)
  return(eff.wgt)
}
