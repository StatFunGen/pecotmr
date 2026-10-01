#' Extract weights from mr.ash.rss (susieR)
#' @param stat A list of summary statistics with elements \code{b} (effect
#'   sizes), \code{seb} (standard errors) and \code{n} (per-variant sample
#'   sizes).
#' @param LD Numeric LD (correlation) matrix aligned to the variants in
#'   \code{stat}.
#' @param varY Numeric. Variance of the phenotype.
#' @param sigma2E Numeric. Residual error variance estimate.
#' @param s0 Numeric vector of prior mixture standard deviations (the sigma0
#'   grid).
#' @param w0 Numeric vector of prior mixture weights (summing to 1).
#' @param z Optional numeric vector of z-scores; defaults to \code{numeric(0)}
#'   (derived from \code{stat}).
#' @param methodArgs Options forwarded to \code{susieR::mr.ash.rss}, built
#'   with \code{\link{mrashConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of the posterior mean of the coefficients.
#' @importFrom susieR mr.ash.rss
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ss <- lapply(seq_len(ncol(X)), function(j) {
#'   coef(summary(lm(y ~ X[, j])))[2, 1:2]
#' })
#' stat <- list(b = vapply(ss, `[`, numeric(1), 1L),
#'   seb = vapply(ss, `[`, numeric(1), 2L), n = rep(nrow(X), ncol(X)))
#' mrashRssWeights(stat, cor(X), varY = var(y), sigma2E = var(y),
#'   s0 = c(0, 0.1, 0.5), w0 = c(0.8, 0.1, 0.1))
#' @importFrom checkmate assertList assertNumeric
#' @export
mrashRssWeights <- function(
    stat,
    LD,
    varY,
    sigma2E,
    s0,
    w0,
    z = numeric(0),
    methodArgs = mrashConfig()
) {
    .assertMethodConfig(methodArgs, "mrashConfig", "methodArgs")
    assertList(stat)
    assertNumeric(z)
    callArgs <- list_modify(
        list(
            bhat = stat$b,
            shat = stat$seb,
            z = z,
            R = LD,
            var_y = varY,
            n = median(stat$n),
            sigma2_e = sigma2E,
            s0 = s0,
            w0 = w0
        ),
        !!!methodArgs
    )
    model <- exec(mr.ash.rss, !!!callArgs)

    return(model$mu1)
}

#' PRS-CS: a polygenic prediction method that infers posterior SNP effect sizes
#' under continuous shrinkage (CS) priors
#'
#' This function is a wrapper for the PRS-CS method implemented in C++. It takes
#' marginal effect size estimates from regression and an external LD reference
#' panel and infers posterior SNP effect sizes using Bayesian regression with
#' continuous shrinkage priors.
#'
#' @param bhat A vector of marginal effect sizes.
#' @param R The LD correlation matrix (a single matrix over the analysed
#'   window), as in \code{susieR::susie_rss()}.
#' @param n Sample size of the GWAS.
#' @param a Shape parameter for the prior distribution of psi. Default is 1.
#' @param b Scale parameter for the prior distribution of psi. Default is 0.5.
#' @param phi Global shrinkage parameter. If NULL, it will be estimated
#'   automatically. Default is NULL.
#' @param nIter Number of MCMC iterations. Default is 1000.
#' @param nBurnin Number of burn-in iterations. Default is 500.
#' @param thin Thinning factor for MCMC. Default is 5.
#' @param maf A vector of minor allele frequencies, if available, will
#'   standardize the effect sizes by MAF. Default is NULL.
#' @param verbose Whether to print verbose output. Default is FALSE.
#' @param seed Random seed for reproducibility. Default is NULL.
#'
#' @return A list containing the posterior estimates: - betaEst: Posterior
#'   estimates of SNP effect sizes. - psiEst: Posterior estimates of psi
#'   (shrinkage parameters). - sigmaEst: Posterior estimate of the residual
#'   variance. - phiEst: Posterior estimate of the global shrinkage parameter.
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
#' out <- prsCs(b.hat, R.hat, n, maf = maf)
#' # In sample prediction correlations
#' cor(X %*% out$betaEst, y) # 0.9944553
#' @export
prsCs <- function(
    bhat,
    R,
    n,
    a = 1,
    b = 0.5,
    phi = NULL,
    maf = NULL,
    nIter = 1000,
    nBurnin = 500,
    thin = 5,
    verbose = FALSE,
    seed = NULL
) {
    # Shared LD-matrix validation, then the prsCs-specific maf length check.
    .rssValidateInputs(bhat, R, n)
    if (!is.null(maf) && length(bhat) != length(maf)) {
        abort("The length of 'bhat' must be the same as 'maf'.")
    }

    # Run PRS-CS
    # cpp11 requires exact integer types for int parameters
    result <- prsCsRcpp(
        a = a,
        b = b,
        phi = phi,
        bhat = bhat,
        maf = maf,
        n = as.integer(n),
        ldBlk = list(blk1 = R),
        nIter = as.integer(nIter),
        nBurnin = as.integer(nBurnin),
        thin = as.integer(thin),
        verbose = verbose,
        seed = seed
    )

    # Return the result as a list (camelCase to match the rest of the package
    # API).
    list(
        betaEst = result$betaEst,
        psiEst = result$psiEst,
        sigmaEst = result$sigmaEst,
        phiEst = result$phiEst
    )
}

#' Extract weights from prsCs function
#' @param stat A list of summary statistics with elements \code{b} (effect
#'   sizes) and \code{n} (per-variant sample sizes).
#' @param LD Numeric LD (correlation) matrix aligned to the variants in
#'   \code{stat}.
#' @param methodArgs Options forwarded to \code{prsCs}, built
#'   with \code{\link{prsCsConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of the posterior SNP coefficients.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ss <- lapply(
#'   seq_len(ncol(X)), function(j) coef(summary(lm(y ~ X[, j])))[2, 1:2])
#' stat <- list(
#'   b = vapply(ss, `[`, numeric(1), 1L),
#'   seb = vapply(ss, `[`, numeric(1), 2L),
#'   n = rep(nrow(X), ncol(X))
#' )
#' LD <- cor(X)
#' prsCsWeights(stat, LD, methodArgs = prsCsConfig(maf = rep(0.3, ncol(X))))
#' @export
prsCsWeights <- function(stat, LD, methodArgs = prsCsConfig()) {
    .assertMethodConfig(methodArgs, "prsCsConfig", "methodArgs")
    callArgs <- list_modify(
        list(bhat = stat$b, R = LD, n = median(stat$n)),
        !!!methodArgs
    )
    model <- exec(prsCs, !!!callArgs)

    return(model$betaEst)
}

#' SDPR (Summary-Statistics-Based Dirichelt Process Regression for Polygenic
#' Risk Prediction)
#'
#' This function is a wrapper for the SDPR C++ implementation, which performs
#' Markov Chain Monte Carlo (MCMC) for estimating effect sizes and heritability
#' based on summary statistics and reference LD matrices.
#'
#' @param bhat A vector of marginal beta values for each SNP.
#' @param R The LD correlation matrix (a single matrix over the analysed
#'   window).
#' @param n The total sample size of the GWAS.
#' @param perVariantSampleSize (Optional) A vector of sample sizes for each SNP.
#'   If NULL (default), it will be initialized to a vector of length equal to
#'   `bhat`, with all values set to `n`.
#' @param array (Optional) A vector of genotyping array information for each
#'   SNP. If NULL (default), it will be initialized to a vector of 1's with
#'   length equal to `bhat`.
#' @param a Factor to shrink the reference LD matrix. Default is 0.1.
#' @param c Factor to correct for the deflation. Default is 1.
#' @param M Max number of variance components. Default is 1000.
#' @param a0k Hyperparameter for inverse gamma distribution. Default is 0.5.
#' @param b0k Hyperparameter for inverse gamma distribution. Default is 0.5.
#' @param iter Number of iterations for MCMC. Default is 1000.
#' @param burn Number of burn-in iterations for MCMC. Default is 200.
#' @param thin Thinning interval for MCMC. Default is 5.
#' @param numThreads Number of threads to use. Default is 1.
#' @param optLlk Which likelihood to evaluate. 1 for equation 6 (slightly shrink
#'   the correlation of SNPs) and 2 for equation 5 (SNPs genotyped on different
#'   arrays in a separate cohort). Default is 1.
#' @param verbose Whether to print verbose output. Default is true.
#'
#' @param seed Integer or \code{NULL}. Random seed for the Gibbs sampler;
#'   \code{NULL} leaves the RNG state unchanged.
#' @return A list with \code{betaEst} (the estimated effect sizes) and
#'   \code{h2} (the estimated heritability). \code{sdprWeights()} returns
#'   only \code{betaEst}; call this directly when you also want \code{h2}.
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
#' out <- sdpr(b.hat, R.hat, n)
#' # In sample prediction correlations
#' cor(X %*% out$betaEst, y) #
#'
#' @note This function wraps a rewritten and adapted version of the SDPR C++
#'   implementation. SDPR is described in Zhou G, Zhao H (2021), "A fast and
#'   robust Bayesian nonparametric method for prediction of complex traits
#'   using summary statistics", PLoS Genetics 17(7): e1009697.
#'   \doi{10.1371/journal.pgen.1009697}
#'
#' @export
sdpr <- function(
    bhat,
    R,
    n,
    perVariantSampleSize = NULL,
    array = NULL,
    a = 0.1,
    c = 1.0,
    M = 1000,
    a0k = 0.5,
    b0k = 0.5,
    iter = 1000,
    burn = 200,
    thin = 5,
    numThreads = 1,
    optLlk = 1,
    verbose = TRUE,
    seed = NULL
) {
    .sdprValidate(bhat, R, n, M, perVariantSampleSize, array)
    # The C++ backend takes a block list; wrap the single-window matrix R as one
    # block. cpp11 requires exact integer types for int params + sexp vectors.
    sdprRcpp(
        bhatR = bhat,
        LD = list(blk1 = R),
        n = as.integer(n),
        perVariantSampleSize = perVariantSampleSize,
        array = if (!is.null(array)) as.integer(array) else NULL,
        a = a,
        c = c,
        M = as.integer(M),
        a0k = a0k,
        b0k = b0k,
        iter = as.integer(iter),
        burn = as.integer(burn),
        thin = as.integer(thin),
        numThreads = as.integer(numThreads),
        optLlk = as.integer(optLlk),
        verbose = verbose,
        seed = seed
    )
}

# Validate the SDPR inputs: RSS contract, M >= 4 (SDPR uses M-2 indexing in
# sample_V; M < 4 overflows), positive per-variant N, and array in {0, 1, 2}.
# @noRd
.sdprValidate <- function(bhat, R, n, M, perVariantSampleSize, array) {
    .rssValidateInputs(bhat, R, n)
    if (M < 4) {
        abort("'M' must be at least 4.")
    }
    if (!is.null(perVariantSampleSize) && any(perVariantSampleSize <= 0)) {
        msg <- glue(
            "The 'perVariantSampleSize' vector must contain only ",
            "positive values."
        )
        abort(msg)
    }
    if (!is.null(array) && any(!is_in(array, c(0, 1, 2)))) {
        abort("The 'array' vector must contain only 0, 1, or 2.")
    }
    invisible(NULL)
}

#' Extract weights from sdpr function
#' @param stat A list of summary statistics with elements \code{b} (effect
#'   sizes) and \code{n} (per-variant sample sizes).
#' @param LD Numeric LD (correlation) matrix aligned to the variants in
#'   \code{stat}.
#' @param methodArgs Options forwarded to \code{sdpr}, built
#'   with \code{\link{sdprConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of the posterior SNP coefficients.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ss <- lapply(
#'   seq_len(ncol(X)), function(j) coef(summary(lm(y ~ X[, j])))[2, 1:2])
#' stat <- list(
#'   b = vapply(ss, `[`, numeric(1), 1L),
#'   seb = vapply(ss, `[`, numeric(1), 2L),
#'   n = rep(nrow(X), ncol(X))
#' )
#' LD <- cor(X)
#' sdprWeights(stat, LD)
#' @export
sdprWeights <- function(stat, LD, methodArgs = sdprConfig()) {
    .assertMethodConfig(methodArgs, "sdprConfig", "methodArgs")
    callArgs <- list_modify(
        list(bhat = stat$b, R = LD, n = median(stat$n)),
        !!!methodArgs
    )
    model <- exec(sdpr, !!!callArgs)

    return(model$betaEst)
}


#' Compute mr.mash TWAS weights
#'
#' Extracts coefficients from an existing mr.mash fit or fits mr.mash from `X`
#' and `Y`.
#'
#' @param mrmashFit Optional fitted mr.mash object.
#' @param X Genotype matrix. Required when `mrmashFit` is NULL.
#' @param Y Phenotype matrix. Required when `mrmashFit` is NULL.
#' @param fitRetention How much of the mr.mash fit is kept as the
#'   \code{"fit"} attribute of the returned weights. \code{"none"} (default)
#'   attaches nothing. \code{"slim"} keeps what \code{fineMappingPipeline}
#'   needs to rebuild the mvSuSiE reweighted mixture prior and residual
#'   variance --- the original \code{dataDrivenPriorMatrices}, the fitted
#'   \code{w0} and \code{V}; the heavy coefficient matrix \code{mu1} is
#'   already returned as the weights, so it is not duplicated. \code{"full"}
#'   additionally retains the complete fit under \code{$fit}, at the cost of
#'   a larger payload.
#' @param methodArgs Options forwarded to \code{mr.mashr::mr.mash}, built
#'   with \code{\link{mrmashConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return Matrix of variant weights.
#' @examples
#' data(multiTraitData)
#' X <- multiTraitData$X[, 1:60]
#' Y <- multiTraitData$Y
#' ddpm <- multiTraitData$priorMatrices
#' fit <- mrmashWrapper(X = X, Y = Y, dataDrivenPriorMatrices = ddpm,
#'   prior = mrmashPriorConfig(canonicalPriorMatrices = TRUE))
#' mrmashWeights(mrmashFit = fit, X = X, Y = Y)
#' @param dataDrivenPriorMatrices Optional list of data-driven prior
#'   covariance matrices; forwarded to \code{mrmashWrapper} when it has to fit,
#'   and retained in the payload for mvSuSiE prior reconstruction.
#' @export
mrmashWeights <- function(
    mrmashFit = NULL,
    X = NULL,
    Y = NULL,
    fitRetention = c("none", "slim", "full"),
    dataDrivenPriorMatrices = NULL,
    methodArgs = mrmashConfig()
) {
    fitRetention <- arg_match(fitRetention)
    .assertMethodConfig(methodArgs, "mrmashConfig", "methodArgs")
    if (!requireNamespace("mr.mashr", quietly = TRUE)) {
        abort("Package 'mr.mashr' is required.")
    }
    if (is.null(mrmashFit)) {
        inform("mrmashFit is not provided; fitting mr.mash now ...")
        if (is.null(X) || is.null(Y)) {
            abort("Both X and Y must be provided if mrmashFit is NULL.")
        }
        # methodArgs is one flat bag; the prior-construction names belong to
        # mrmashWrapper's `prior` group, so they are routed there rather
        # than passed as top-level arguments it no longer has.
        split <- .mrmashSplitPriorArgs(as.list(methodArgs))
        mrmashFit <- exec(
            mrmashWrapper,
            X,
            Y,
            dataDrivenPriorMatrices = dataDrivenPriorMatrices,
            prior = exec(mrmashPriorConfig, !!!split$prior),
            !!!split$rest
        )
    }
    out <- mr.mashr::coef.mr.mash(mrmashFit)[-1, ]
    # mu1 (= out) is already the returned weights; the payload carries only the
    # mvSuSiE data-driven-prior reconstruction inputs (w0, V, and the raw prior
    # matrices), plus the whole fit at fitRetention = "full".
    .mrmashAttachFit(
        out,
        mrmashFit,
        dataDrivenPriorMatrices,
        fitRetention
    )
}


# Attach the mvSuSiE-prior reconstruction payload as the "fit" attribute,
# shared by mrmashWeights (individual) and mrmashRssWeights (summary stats):
# the data-driven prior matrices + fitted w0 + V, and -- at
# fitRetention = "full" -- the whole fit. The coefficients are already the
# returned `weights`, so mu1 is not duplicated. Returns `weights` unchanged
# at fitRetention = "none".
.mrmashAttachFit <- function(
    weights,
    fit,
    dataDrivenPriorMatrices,
    fitRetention
) {
    if (identical(fitRetention, "none")) {
        return(weights)
    }
    fitList <- c(
        list(
            dataDrivenPriorMatrices = dataDrivenPriorMatrices,
            w0 = fit$w0,
            V = fit$V
        ),
        compact(list(fit = if (fitRetention == "full") fit))
    )
    `attr<-`(weights, "fit", fitList)
}

# The `beta.init` override for mr.ash: lasso weights when the caller supplied
# none, or the caller's own initialisation restricted to the retained columns.
# @noRd
.mrashBetaInit <- function(methodArgs, XKeep, y, X, keep) {
    if (!is_in("beta.init", names(methodArgs))) {
        return(list(beta.init = lassoWeights(XKeep, y)))
    }
    if (length(methodArgs$beta.init) != ncol(X)) {
        return(list())
    }
    list(beta.init = methodArgs$beta.init[keep])
}

#' Compute mr.mash-RSS TWAS weights from summary statistics
#'
#' Multi-context summary-statistics analog of \code{\link{mrmashWeights}}:
#' extracts coefficients from an existing \code{mr.mashr::mr.mash.rss} fit, or
#' fits one from \code{stat} (variants x conditions) and \code{LD}.
#'
#' Follows the \code{*_rss_weights(stat, LD, ...)} contract. Expects
#' \code{stat$z} to be a numeric matrix (variants x conditions) and
#' \code{stat$n} a per-context numeric vector or scalar. \code{stat$Bhat} and
#' \code{stat$Shat} are used if present; otherwise derived from Z and n.
#'
#' Prior construction reuses the same infrastructure as the individual-level
#' \code{\link{mrmashWrapper}}: \code{computeGrid} +
#' \code{mr.mashr::compute_canonical_covs()} + \code{mr.mashr::expand_covs()}
#' for \code{S0}, and \code{computeW0} for the mixture weights. Supply
#' \code{dataDrivenPriorMatrices} (e.g. from
#' \code{\link{computeCovDiag}}, or \code{mashr::cov_flash}) to add
#' data-driven covariance components alongside the canonical mixture.
#'
#' @param stat A list with \code{z} (variants x conditions matrix) and \code{n}
#'   (per-context numeric vector or scalar). May also include \code{Bhat},
#'   \code{Shat} matrices.
#' @param LD LD correlation matrix.
#' @param mrmashRssFit Optional pre-fitted \code{mr.mash.rss} object; skips
#'   fitting when supplied.
#' @param dataDrivenPriorMatrices Optional list with element \code{U} (list of
#'   raw covariance matrices). Passed directly to \code{mr.mashr::expand_covs()}
#'   alongside the canonical mixture.
#' @param canonicalPriorMatrices Logical. When TRUE (default), include the
#'   standard canonical mixture from \code{mr.mashr::compute_canonical_covs()}.
#'   When FALSE, \code{dataDrivenPriorMatrices} must be supplied.
#' @param S0 Optional pre-built list of prior covariance matrices, bypassing the
#'   canonical / data-driven construction.
#' @param w0 Optional prior mixture weights; defaults to \code{computeW0(Bhat,
#'   length(S0))}.
#' @param V Optional residual covariance matrix (K x K). When NULL, defaults to
#'   the identity matrix of size K.
#' @param covY Optional response covariance matrix (K x K). When NULL, defaults
#'   to the identity matrix of size K.
#' @param fitRetention How much of the mr.mash fit is kept as the
#'   \code{"fit"} attribute of the returned weights. \code{"none"} (default)
#'   attaches nothing. \code{"slim"} keeps what \code{fineMappingPipeline}
#'   needs to rebuild the mvSuSiE reweighted mixture prior and residual
#'   variance --- the original \code{dataDrivenPriorMatrices}, the fitted
#'   \code{w0} and \code{V}; the heavy coefficient matrix \code{mu1} is
#'   already returned as the weights, so it is not duplicated. \code{"full"}
#'   additionally retains the complete fit under \code{$fit}, at the cost of
#'   a larger payload.
#' @param methodArgs Options forwarded to \code{mr.mashr::mr.mash.rss}, built
#'   with \code{\link{mrmashConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric matrix of per-variant per-context weights (variants x
#'   conditions).
#' @examples
#' data(eqtlRegionExample)
#' data(multiTraitData)
#' X <- eqtlRegionExample$X[, 1:30]
#' cond <- c("brain", "blood", "muscle")
#' zmat <- matrix(rnorm(30 * 3), 30, 3,
#'   dimnames = list(colnames(X), cond))
#' stat <- list(z = zmat, n = matrix(nrow(X), 30, 3))
#' mrmashRssWeights(stat, cor(X),
#'   dataDrivenPriorMatrices = multiTraitData$priorMatrices,
#'   canonicalPriorMatrices = TRUE)
#' @export
#' @importFrom checkmate assertFlag assertList
mrmashRssWeights <- function(
    stat,
    LD,
    mrmashRssFit = NULL,
    dataDrivenPriorMatrices = NULL,
    canonicalPriorMatrices = TRUE,
    S0 = NULL,
    w0 = NULL,
    V = NULL,
    covY = NULL,
    fitRetention = c("none", "slim", "full"),
    methodArgs = mrmashConfig()
) {
    fitRetention <- arg_match(fitRetention)
    .assertMethodConfig(methodArgs, "mrmashConfig", "methodArgs")
    assertList(stat)
    assertFlag(canonicalPriorMatrices)
    .mrmashRssRequirePackage()
    if (is.null(mrmashRssFit)) {
        mrmashRssFit <- .mrmashRssComputeFit(
            stat,
            LD,
            dataDrivenPriorMatrices,
            canonicalPriorMatrices,
            S0,
            w0,
            V,
            covY,
            methodArgs
        )
    }
    # coef.mr.mash.rss returns nrow(Bhat) rows (no intercept). Do not strip.
    weights <- mr.mashr::coef.mr.mash.rss(mrmashRssFit)
    .mrmashAttachFit(
        weights,
        mrmashRssFit,
        dataDrivenPriorMatrices,
        fitRetention
    )
}

# Require mr.mashr for the RSS mr.mash path.
# @noRd
.mrmashRssRequirePackage <- function() {
    if (!requireNamespace("mr.mashr", quietly = TRUE)) {
        msg <- glue(
            "Package 'mr.mashr' is required. ",
            "is required."
        )
        abort(msg)
    }
}

# Derive the (Z, K, nVec, Bhat, Shat) inputs from a multi-context stat list.
# @noRd
.mrmashRssStats <- function(stat) {
    Z <- if (is.matrix(stat$z)) stat$z else as.matrix(stat$z)
    if (ncol(Z) < 2) {
        msg <- glue(
            "mrmashRssWeights expects stat$z to have >= 2 columns ",
            "(one per context). For single-context use mrashRssWeights()."
        )
        abort(msg)
    }
    K <- ncol(Z)
    nVec <- if (length(stat$n) > 1) stat$n else rep(stat$n, K)
    Bhat <- if (!is.null(stat$Bhat)) stat$Bhat else sweep(Z, 2, sqrt(nVec), "/")
    Shat <- if (!is.null(stat$Shat)) {
        stat$Shat
    } else {
        matrix(1 / sqrt(rep(nVec, each = nrow(Z))), nrow = nrow(Z), ncol = K)
    }
    list(Z = Z, K = K, nVec = nVec, Bhat = Bhat, Shat = Shat)
}

# Fit mr.mash.rss, defaulting the prior / w0 / V / covY when not supplied.
# @noRd
.mrmashRssComputeFit <- function(
    stat,
    LD,
    dataDrivenPriorMatrices,
    canonicalPriorMatrices,
    S0,
    w0,
    V,
    covY,
    dots
) {
    ss <- .mrmashRssStats(stat)
    if (is.null(S0)) {
        S0 <- buildMrmashPriorMatrices(
            Bhat = ss$Bhat,
            Shat = ss$Shat,
            K = ss$K,
            dataDrivenPriorMatrices = dataDrivenPriorMatrices,
            canonicalPriorMatrices = canonicalPriorMatrices
        )$S0
    }
    if (is.null(w0)) {
        w0 <- computeW0(ss$Bhat, length(S0))
    }
    if (is.null(V)) {
        V <- diag(ss$K)
    }
    if (is.null(covY)) {
        covY <- diag(ss$K)
    }
    # mr.mash.rss expects either Z or (Bhat, Shat) but not both; prefer
    # Bhat/Shat. n must be a scalar (mr.mash.rss contract); use the median.
    rssConfig <- c(
        list(
            Bhat = ss$Bhat,
            Shat = ss$Shat,
            R = LD,
            n = as.numeric(stats::median(ss$nVec)),
            covY = covY,
            V = V,
            S0 = S0,
            w0 = w0
        ),
        # `dots` is the caller's mrmashConfig() record; c() would append it as
        # one opaque element instead of splicing its entries.
        as.list(dots)
    )
    exec(mr.mashr::mr.mash.rss, !!!rssConfig)
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
        nDrop <- sum(!keep)
        idxStr <- str_flatten(which(!keep), ", ")
        msg <- glue(
            "{fnName}: dropping {nDrop} zero-variance column(s) from X ",
            "(indices: {idxStr})"
        )
        warn(msg)
    }
    keep
}

# --- TWAS weight-method argument constructors -------------------------------
#
# One constructor per fitting engine, not per method token: several tokens
# share an engine (lasso and enet both fit with glmnet; bayesA/C/N/R all use
# qgg), and it is the engine that defines what the options are.

#' @title Arguments For A TWAS Weight-Fitting Engine
#' @description Options forwarded to the engine behind a TWAS weight method,
#'   checked against that engine's live formals where it enumerates them. Pass
#'   the result as an element of \code{\link{twasWeightsMethodsConfig}}, or
#'   directly to the matching weights function.
#'
#'   A method's individual-level and summary-statistic paths often run
#'   DIFFERENT engines, and then each takes its own constructor:
#'   \code{lassoWeights} fits with glmnet (\code{glmnetConfig}) while
#'   \code{lassosumRssWeights} runs pecotmr's lassosum solver
#'   (\code{lassosumConfig}); \code{scadWeights} / \code{mcpWeights} /
#'   \code{l0learnWeights} use ncvreg and L0Learn while their \code{*Rss}
#'   counterparts use \code{penalizedRssConfig}; \code{dprGibbsWeights} uses
#'   RcppDPR and \code{sdprWeights} uses \code{sdprConfig}. Splitting them is
#'   what lets each be checked exactly, rather than against a union in which
#'   an option meant for the other path would pass and then be dropped.
#'
#'   \code{glmnetConfig}, \code{ncvregConfig} and \code{dprConfig} cannot be
#'   checked: \code{cv.glmnet}, \code{cv.ncvreg} and \code{fit_model} each
#'   take \code{...}, so they accept any name and there is nothing to check
#'   against. Every other constructor here rejects an unknown option. See
#'   \code{\link{MethodConfig}}.
#' @param ... Arguments for the engine, under its own names.
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' glmnetConfig(nfold = 10)
#' qggConfig(nit = 1000)
#' lassosumConfig(thr = 1e-4)
#' @name twasEngineConfig
NULL

#' @rdname twasEngineConfig
#' @export
glmnetConfig <- function(...) {
    .newMethodConfig(
        "glmnet::cv.glmnet",
        defaults = list(),
        extra = list(...),
        label = "glmnetConfig",
        engine = "glmnet"
    )
}

#' @rdname twasEngineConfig
#' @export
qggConfig <- function(...) {
    .newMethodConfig(
        "qgg::gbayes",
        defaults = list(),
        extra = list(...),
        label = "qggConfig",
        engine = "qgg"
    )
}

#' @rdname twasEngineConfig
#' @export
bglrConfig <- function(...) {
    .newMethodConfig(
        "BGLR::BGLR",
        defaults = list(),
        extra = list(...),
        label = "bglrConfig",
        engine = "bglr"
    )
}

#' @rdname twasEngineConfig
#' @export
dprConfig <- function(...) {
    .newMethodConfig(
        "RcppDPR::fit_model",
        defaults = list(),
        extra = list(...),
        label = "dprConfig",
        engine = "dpr"
    )
}

#' @rdname twasEngineConfig
#' @export
ncvregConfig <- function(...) {
    .newMethodConfig(
        "ncvreg::cv.ncvreg",
        defaults = list(),
        extra = list(...),
        label = "ncvregConfig",
        engine = "ncvreg"
    )
}

#' @rdname twasEngineConfig
#' @export
l0learnConfig <- function(...) {
    .newMethodConfig(
        "L0Learn::L0Learn.cvfit",
        defaults = list(),
        extra = list(...),
        label = "l0learnConfig",
        engine = "l0learn"
    )
}

#' @rdname twasEngineConfig
#' @export
mrashConfig <- function(...) {
    .newMethodConfig(
        c("susieR::mr.ash", "susieR::mr.ash.rss"),
        defaults = list(),
        extra = list(...),
        label = "mrashConfig",
        engine = "mrash"
    )
}

#' @rdname twasEngineConfig
#' @export
prsCsConfig <- function(...) {
    .newMethodConfig(
        "pecotmr::prsCs",
        defaults = list(),
        extra = list(...),
        label = "prsCsConfig",
        engine = "prsCs"
    )
}

#' @rdname twasEngineConfig
#' @export
lassosumConfig <- function(...) {
    .newMethodConfig(
        "pecotmr::lassosumRss",
        defaults = list(),
        extra = list(...),
        label = "lassosumConfig",
        engine = "lassosum"
    )
}

#' @rdname twasEngineConfig
#' @export
penalizedRssConfig <- function(...) {
    .newMethodConfig(
        "pecotmr::penalizedRss",
        defaults = list(),
        extra = list(...),
        label = "penalizedRssConfig",
        engine = "penalizedRss"
    )
}

#' @rdname twasEngineConfig
#' @export
sdprConfig <- function(...) {
    .newMethodConfig(
        "pecotmr::sdpr",
        defaults = list(),
        extra = list(...),
        label = "sdprConfig",
        engine = "sdpr"
    )
}

#' @rdname twasEngineConfig
#' @export
mrmashConfig <- function(...) {
    .newMethodConfig(
        c("mr.mashr::mr.mash", "mr.mashr::mr.mash.rss"),
        defaults = list(),
        extra = list(...),
        label = "mrmashConfig",
        engine = "mrmash",
        # The bag is split before the engine sees it: .mrmashSplitPriorArgs
        # routes the prior-covariance names to mrmashPriorConfig() and the rest
        # to mr.mash. So the accepted set is both, read from the live
        # formals of each rather than transcribed.
        accepted = unique(c(
            .engineAcceptedNames(
                c("mr.mashr::mr.mash", "mr.mashr::mr.mash.rss")
            ),
            names(formals(mrmashPriorConfig))
        ))
    )
}

#' Compute TWAS weights via penalized regression (glmnet)
#'
#' Fit a cross-validated elastic-net / lasso model with \code{glmnet::cv.glmnet}
#' and return the per-variant coefficient vector at \code{lambda.min}.
#'
#' @param X Numeric genotype / design matrix (samples x variants).
#' @param y Numeric response (phenotype) vector of length \code{nrow(X)}.
#' @param alpha Elastic-net mixing parameter in \code{[0, 1]}: \code{1} = lasso,
#'   \code{0} = ridge, \code{0.5} = elastic net.
#' @param methodArgs Additional arguments forwarded to
#'   \code{glmnet::cv.glmnet}, built with \code{\link{glmnetConfig}}. A bare
#'   list is refused, since it cannot be checked against the engine's formals.
#' @return Numeric vector of length \code{ncol(X)} of per-variant weights.
#' @importFrom stats coef
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' glmnetWeights(X, y, alpha = 0.5)
#' @export
glmnetWeights <- function(X, y, alpha, methodArgs = glmnetConfig()) {
    # Check if glmnet is installed
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        abort("Package 'glmnet' is required for this function.")
    }
    .assertMethodConfig(methodArgs, "glmnetConfig", "methodArgs")
    eff.wgt <- matrix(0, ncol = 1, nrow = ncol(X))
    keep <- .dropZeroVariance(X, "glmnetWeights")
    # `alpha` is the wrapper's own dial -- it is what distinguishes lasso from
    # enet from ridge -- so it is supplied here rather than settable in
    # methodArgs. Everything else pecotmr fixed is now overridable.
    userArgs <- as.list(methodArgs)
    if (is_in("alpha", names(userArgs))) {
        abort(glue(
            "glmnetWeights: `alpha` selects the penalty family and is set by ",
            "the calling wrapper (lassoWeights = 1, enetWeights = 0.5), so ",
            "it cannot also be given in `methodArgs`."
        ))
    }
    callArgs <- list_modify(
        list(
            x = X[, keep, drop = FALSE],
            y = y,
            alpha = alpha,
            nfold = 5,
            intercept = TRUE,
            standardize = FALSE
        ),
        !!!userArgs
    )
    enet <- exec(glmnet::cv.glmnet, !!!callArgs)
    replace(eff.wgt, keep, coef(enet, s = "lambda.min")[2:(sum(keep) + 1)])
}

#' Compute TWAS weights via elastic net (glmnet, alpha = 0.5)
#'
#' Convenience wrapper for \code{\link{glmnetWeights}} with \code{alpha = 0.5}.
#' @inheritParams glmnetWeights
#' @return Numeric vector of length \code{ncol(X)} of per-variant weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' enetWeights(X, y)
#' @export
enetWeights <- function(X, y, methodArgs = glmnetConfig()) {
    .assertMethodConfig(methodArgs, "glmnetConfig", "methodArgs")
    glmnetWeights(X, y, 0.5, methodArgs = methodArgs)
}

#' Compute TWAS weights via lasso (glmnet, alpha = 1)
#'
#' Convenience wrapper for \code{\link{glmnetWeights}} with \code{alpha = 1}.
#' @inheritParams glmnetWeights
#' @return Numeric vector of length \code{ncol(X)} of per-variant weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' lassoWeights(X, y)
#' @export
lassoWeights <- function(X, y, methodArgs = glmnetConfig()) {
    .assertMethodConfig(methodArgs, "glmnetConfig", "methodArgs")
    glmnetWeights(X, y, 1, methodArgs = methodArgs)
}

#' Compute Weights Using mr.ash Shrinkage
#'
#' This function fits the `mr.ash` model (adaptive shrinkage regression) to
#' estimate weights for a given set of predictors and response. It uses optional
#' prior standard deviation initialization and can accept custom initial beta
#' values.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(200 * 10), nrow = 200)
#' y <- as.numeric(X %*% rnorm(10) + rnorm(200))
#' weights <- mrashWeights(X, y)
#' @importFrom susieR mr.ash
#' @importFrom stats predict
#' @param X Numeric genotype / design matrix (samples x variants).
#' @param y Numeric response (phenotype) vector of length \code{nrow(X)}.
#' @param initPriorSd Logical. Initialize the prior standard-deviation grid
#' from the data. Default \code{TRUE}.
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights: \code{"none"} (default) attaches
#'   nothing, \code{"slim"} and \code{"full"} attach the fitted model. This
#'   engine keeps no trimmable intermediates, so the two retaining levels
#'   behave alike here.
#' Default \code{FALSE}.
#' @param methodArgs Options forwarded to \code{susieR::mr.ash}, built
#'   with \code{\link{mrashConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of weights, one per variant (column of \code{X});
#'   zero-variance columns receive weight 0. Unless
#'   \code{fitRetention = "none"} the fitted \code{mr.ash} object is
#'   attached as attribute \code{"fit"}.
#' @importFrom checkmate assertFlag
#' @export
mrashWeights <- function(
    X,
    y,
    initPriorSd = TRUE,
    fitRetention = c("none", "slim", "full"),
    methodArgs = mrashConfig()
) {
    .assertMethodConfig(methodArgs, "mrashConfig", "methodArgs")
    assertFlag(initPriorSd)
    fitRetention <- arg_match(fitRetention)
    keep <- .dropZeroVariance(X, "mrashWeights")
    XKeep <- X[, keep, drop = FALSE]
    # as.list(): list_assign() rejects a MethodConfig outright as `.x`, the
    # same S4-SimpleList seam that makes c() append it instead of splicing.
    argsList <- list_assign(
        as.list(methodArgs),
        !!!.mrashBetaInit(methodArgs, XKeep, y, X, keep)
    )
    mrashConfig <- c(
        list(
            X = XKeep,
            y = y,
            sa2 = if (initPriorSd) initPriorSd(XKeep, y)^2 else NULL
        ),
        argsList
    )
    fit.mr.ash <- exec(mr.ash, !!!mrashConfig)
    # Zero-variance columns were never fitted and keep a zero weight.
    eff.wgt <- replace(
        rep(0, ncol(X)),
        keep,
        predict(fit.mr.ash, type = "coefficients")[-1]
    )
    if (identical(fitRetention, "none")) {
        return(eff.wgt)
    }
    `attr<-`(eff.wgt, "fit", fit.mr.ash)
}
#' Extract Coefficients From Bayesian Linear Regression
#'
#' This function performs Bayesian linear regression using the `gbayes` function
#' from the `qgg` package. It then returns the estimated slopes.
#'
#' @param y A numeric vector of phenotypes.
#' @param X A numeric matrix of genotypes.
#' @param method A character string declaring the method/prior to be used.
#'   Options are bayesN, bayesL, bayesA, bayesC, or bayesR.
#' @param Z An optional numeric matrix of covariates.
#' @param h2 Numeric or \code{NULL}. Prior heritability for the sampler;
#'   \code{NULL} lets \code{qgg} estimate it.
#' @param nit Integer. Total number of MCMC iterations. Default \code{5000}.
#' @param nburn Integer. Number of burn-in iterations discarded. Default
#'   \code{1000}.
#' @param nthin Integer. Thinning interval for retained MCMC samples. Default
#'   \code{5}.
#' @param methodArgs Options forwarded to \code{qgg::gbayes}, built
#'   with \code{\link{qggConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A vector containing the weights to be applied to each genotype in
#'   predicting the phenotype.
#' @details This function fits a Bayesian linear regression model with a range
#'   of priors.
#' @examples
#' X <- matrix(rnorm(100000), nrow = 1000)
#' Z <- matrix(round(runif(3000, 0, 0.8), 0), nrow = 1000)
#' set1 <- sample(seq_len(ncol(X)), 5)
#' set2 <- sample(seq_len(ncol(X)), 5)
#' sets <- list(set1, set2)
#' g <- rowSums(X[, c(set1, set2)])
#' e <- rnorm(nrow(X), mean = 0, sd = 1)
#' y <- g + e
#' bayesLWeights(y = y, X = X, Z = Z)
#' bayesRWeights(y = y, X = X, Z = Z)
#' @export
bayesAlphabetWeights <- function(
    X,
    y,
    method,
    Z = NULL,
    h2 = NULL,
    nit = 5000,
    nburn = 1000,
    nthin = 5,
    methodArgs = qggConfig()
) {
    .assertMethodConfig(methodArgs, "qggConfig", "methodArgs")
    .bayesAlphabetValidate(X, y, Z)

    eff.wgt <- rep(0, ncol(X))
    keep <- .dropZeroVariance(X, "bayesAlphabetWeights")

    callArgs <- list_modify(
        list(
            y = y,
            W = X[, keep, drop = FALSE],
            X = Z,
            method = method,
            h2 = h2,
            nit = nit,
            nburn = nburn,
            nthin = nthin
        ),
        !!!methodArgs
    )
    model <- exec(qgg::gbayes, !!!callArgs)
    replace(eff.wgt, keep, model$bm)
}

# Shared input validation for the gbayes-backed weight fitters: qgg present,
# and matching row counts for response / genotype / covariates.
# @noRd
#' @importFrom checkmate assertMatrix assertVector
.bayesAlphabetValidate <- function(X, y, Z) {
    if (!requireNamespace("qgg", quietly = TRUE)) {
        abort("Package 'qgg' is required for this function.")
    }
    assertVector(y, len = nrow(X))
    assertMatrix(Z, nrows = nrow(X), null.ok = TRUE)
}
#' @title BayesN TWAS weights (Gaussian prior, ridge-equivalent)
#' @description Use Gaussian distribution as prior. Posterior means will be
#'   BLUP, equivalent to Ridge Regression.
#' @param X Numeric genotype / design matrix (samples x variants).
#' @param y Numeric response (phenotype) vector of length \code{nrow(X)}.
#' @param Z Optional numeric matrix of fixed-effect covariates, or \code{NULL}.
#' @param methodArgs Options forwarded to \code{qgg::gbayes}, built
#'   with \code{\link{qggConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of effect-size weights, one per variant (column of
#'   \code{X}); columns dropped for zero variance receive weight 0.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' bayesNWeights(X, y)
#' @export
bayesNWeights <- function(X, y, Z = NULL, methodArgs = qggConfig()) {
    .assertMethodConfig(methodArgs, "qggConfig", "methodArgs")
    bayesAlphabetWeights(X, y, method = "bayesN", Z, methodArgs = methodArgs)
}
#' @title BayesL TWAS weights (Laplace prior, LASSO-equivalent)
#' @description Use laplace/double exponential distribution as prior. This is
#'   equivalent to Bayesian LASSO.
#' @param X Numeric genotype / design matrix (samples x variants).
#' @param y Numeric response (phenotype) vector of length \code{nrow(X)}.
#' @param Z Optional numeric matrix of fixed-effect covariates, or \code{NULL}.
#' @param methodArgs Options forwarded to \code{qgg::gbayes}, built
#'   with \code{\link{qggConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of effect-size weights, one per variant (column of
#'   \code{X}); columns dropped for zero variance receive weight 0.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' bayesLWeights(X, y)
#' @export
bayesLWeights <- function(X, y, Z = NULL, methodArgs = qggConfig()) {
    .assertMethodConfig(methodArgs, "qggConfig", "methodArgs")
    bayesAlphabetWeights(X, y, method = "bayesL", Z, methodArgs = methodArgs)
}
#' @title BayesA TWAS weights (t-distribution prior)
#' @description Use t-distribution as prior.
#' @param X Numeric genotype / design matrix (samples x variants).
#' @param y Numeric response (phenotype) vector of length \code{nrow(X)}.
#' @param Z Optional numeric matrix of fixed-effect covariates, or \code{NULL}.
#' @param methodArgs Options forwarded to \code{qgg::gbayes}, built
#'   with \code{\link{qggConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of effect-size weights, one per variant (column of
#'   \code{X}); columns dropped for zero variance receive weight 0.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' bayesAWeights(X, y)
#' @export
bayesAWeights <- function(X, y, Z = NULL, methodArgs = qggConfig()) {
    .assertMethodConfig(methodArgs, "qggConfig", "methodArgs")
    bayesAlphabetWeights(X, y, method = "bayesA", Z, methodArgs = methodArgs)
}
#' @title BayesC TWAS weights (rounded-spike prior)
#' @description Use a rounded spike prior (low-variance Gaussian).
#' @param X Numeric genotype / design matrix (samples x variants).
#' @param y Numeric response (phenotype) vector of length \code{nrow(X)}.
#' @param Z Optional numeric matrix of fixed-effect covariates, or \code{NULL}.
#' @param pi Numeric in (0, 1). Prior proportion of non-null effects for the
#'   BayesC mixture. Default \code{0.1}.
#' @param methodArgs Options forwarded to \code{qgg::gbayes}, built
#'   with \code{\link{qggConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of effect-size weights, one per variant (column of
#'   \code{X}); columns dropped for zero variance receive weight 0.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' bayesCWeights(X, y)
#' @export
bayesCWeights <- function(X, y, Z = NULL, pi = 0.1, methodArgs = qggConfig()) {
    .assertMethodConfig(methodArgs, "qggConfig", "methodArgs")
    # `pi` is a qgg option, not a formal of bayesAlphabetWeights, so it joins
    # the option list rather than the argument list.
    bayesAlphabetWeights(
        X,
        y,
        method = "bayesC",
        Z,
        methodArgs = .methodConfigWith(methodArgs, "qggConfig", pi = pi)
    )
}
#' @title BayesR TWAS weights (hierarchical mixture prior)
#' @description Use a hierarchical Bayesian mixture model with four Gaussian
#'   components. Variances are scaled by 0, 0.0001, 0.001, and 0.01.
#' @param X Numeric genotype / design matrix (samples x variants).
#' @param y Numeric response (phenotype) vector of length \code{nrow(X)}.
#' @param Z Optional numeric matrix of fixed-effect covariates, or \code{NULL}.
#' @param methodArgs Options forwarded to \code{qgg::gbayes}, built
#'   with \code{\link{qggConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of effect-size weights, one per variant (column of
#'   \code{X}); columns dropped for zero variance receive weight 0.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' bayesRWeights(X, y)
#' @export
bayesRWeights <- function(X, y, Z = NULL, methodArgs = qggConfig()) {
    .assertMethodConfig(methodArgs, "qggConfig", "methodArgs")
    bayesAlphabetWeights(X, y, method = "bayesR", Z, methodArgs = methodArgs)
}


# #' Bayesian linear regression using summary statistics
# #'
# #' @description
# #'
# #' This function is adapted from those written by Peter Sorensen in the qgg
# #' package.
# #' The following prior distributions are provided:
# #'
# #' Bayes N: Assigning a Gaussian prior to marker effects implies that the
# #' posterior means are the
# #' BLUP estimates (same as Ridge Regression).
# #'
# #' Bayes L: Assigning a double-exponential or Laplace prior is the density
# #' used in
# #' the Bayesian LASSO
# #'
# #' Bayes A: similar to ridge regression but t-distribution prior (rather than
# #' Gaussian)
# #' for the marker effects ; variance comes from an inverse-chi-square
# #' distribution instead of being fixed. Estimation
# #' via Gibbs sampling.
# #'
# #' Bayes C: uses a "rounded spike" (low-variance Gaussian) at origin many
# #' small
# #' effects can contribute to polygenic component, reduces the dimensionality
# #' of
# #' the model (makes Gibbs sampling feasible).
# #'
# #' Bayes R: Hierarchical Bayesian mixture model with 4 Gaussian components,
# #' with
# #' variances scaled by 0, 0.0001 , 0.001 , and 0.01 .
# #'
# #' @param sumstats dataframe with marker summary statistics. Required: beta
# #' coefficient (beta), standard
# #'        error of the beta coefficient (se), GWAS sample size (n). Optional:
# #'        variant_id or rsid, alleles (A1
# #'        and A2), minor allele frequency (maf).
# #' @param LD is a the LD matrix corresponding to the same markers as in the
# #' stat dataframe
# #' @param variant_ids is an optional character vector of variant ids or
# #' rsids, provided outside of the rss dataframe
# #' @param nit is the number of iterations
# #' @param nburn is the number of burnin iterations
# #' @param nthin is the thinning parameter
# #' @param method specifies the methods used
# #' (method="bayesN","bayesA","bayesL","bayesC","bayesR")
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
# #' @param nug is a scalar or vector of prior degrees of freedom for prior
# #' genetic (co)variances
# #' @param nub is a scalar or vector of prior degrees of freedom for marker
# #' (co)variances
# #' @param nue is a scalar or vector of prior degrees of freedom for prior
# #' residual (co)variances
# #' @param mask is a vector or matrix of TRUE/FALSE specifying if marker
# #' should be ignored
# #' @param ve_prior is a scalar or matrix of prior residual (co)variances
# #' @param vg_prior is a scalar or matrix of prior genetic (co)variances
# #' @param algorithm is the algorithm to use. Should take on values ("mcmc",
# #' "em-mcmc")
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
# #' \item{param}{a list current parameters (same information as item listed
# #' above)
# #'              used for restart of the analysis}
# #' \item{stat}{matrix (mxt) of marker information and effects used for
# #' genomic risk scoring}
# #' \item{method}{the method used}
# #' \item{mask}{which loci were masked from analysis}
# #' \item{conv}{dataframe of convergence metrics}
# #' \item{post}{posterior parameter estimates}
# #' \item{ve}{mean residual variance}
# #' \item{vg}{mean genomic variance}
# #'
# #' @export
# gbayes_rss <- function(
#     sumstats = NULL,
#     LD = NULL,
#     variant_ids = NULL,
#     nit = 100,
#     nburn = 0,
#     nthin = 4,
#     method = "bayesR",
#     vg = NULL,
#     vb = NULL,
#     ve = NULL,
#     ssg_prior = NULL,
#     ssb_prior = NULL,
#     sse_prior = NULL,
#     lambda = NULL,
#     h2 = NULL,
#     pi = 0.001,
#     updateB = TRUE,
#     updateG = TRUE,
#     updateE = TRUE,
#     updatePi = TRUE,
#     adjustE = TRUE,
#     nug = 4,
#     nub = 4,
#     nue = 4,
#     mask = NULL,
#     ve_prior = NULL,
#     vg_prior = NULL,
#     algorithm = "mcmc",
#     tol = 0.001,
#     nit_local = NULL,
#     nit_global = NULL
# ) {
#     # Make sure qgg is installed
#     if (!requireNamespace("qgg", quietly = TRUE)) {
#         stop(
#             paste0(
#                 "To use this function, please install qgg: ",
#
#             )
#         )
#     }
#     # Check methods
#     methods <- c("bayesN", "bayesA", "bayesL", "bayesC", "bayesR")
#     method <- match(method, methods)
#     if (!sum(method %in% c(1:5)) == 1) {
#         stop("Method specified not valid")
#     }
#     if (method == 0) {
#         # BLUP and we do not estimate parameters
#         updateB <- FALSE
#         updateE <- FALSE
#     }
#
#     # Set algorithm
#     if (algorithm == "em-mcmc") {
#         algo <- 2
#     } else {
#         algo <- 1
#     }
#
#     # Check that LD matrix is provided and of same length as stats
#     if (is.null(LD)) {
#         stop("Must provide LD matrix")
#     }
#     if (nrow(sumstats) != nrow(LD)) {
#         stop("LD matrix must correspond to summary statistics")
#     }
#
#     # Parameters from stat df
#     if (is.data.frame(sumstats)) {
#         if (!is.null(variant_ids)) {
#             variant_ids <- variant_ids
#         } else if (!is.null(sumstats$rsids)) {
#             variant_ids <- sumstats$rsids
#         } else if (!is.null(sumstats$variant_id)) {
#             variant_ids <- sumstats$variant_id
#         } else {
#             variant_ids <- paste0("snp", 1:nrow(sumstats))
#             sumstats$variant_id <- variant_ids
#         }
#
#         m <- length(variant_ids)
#         b <- wy <- ww <- matrix(0, nrow = nrow(sumstats), ncol = 1)
#         mask <- matrix(FALSE, nrow = nrow(sumstats), ncol = 1)
#         rownames(b) <- rownames(wy) <- rownames(ww) <- rownames(
#             mask
#         ) <- variant_ids
#
#         if (is.null(sumstats$ww)) {
#             sumstats$ww <- 1 / (sumstats$se^2 + sumstats$beta^2 / sumstats$n)
#         }
#         if (is.null(sumstats$wy)) {
#             sumstats$wy <- sumstats$beta * sumstats$ww
#         }
#         if (!is.null(sumstats$n)) {
#             n <- as.integer(median(sumstats$n))
#         }
#         ww[, 1] <- sumstats$ww
#         wy[, 1] <- sumstats$wy
#         mask[, 1] <- FALSE
#
#         if (any(is.na(wy))) {
#             stop("Missing values in wy")
#         }
#         if (any(is.na(ww))) {
#             stop("Missing values in ww")
#         }
#
#         b2 <- sumstats$beta^2
#         seb2 <- sumstats$se^2
#         yy <- (b2 + (n - 2) * seb2) * sumstats$ww
#         yy <- median(yy)
#
#         if (is.null(sumstats$A1)) {
#             sumstats$A1 <- rep("Unknown", length = nrow(sumstats))
#         }
#         if (is.null(sumstats$A2)) {
#             sumstats$A2 <- rep("Unknown", length = nrow(sumstats))
#         }
#         if (is.null(sumstats$maf)) {
#             sumstats$maf <- rep("Unknown", length = nrow(sumstats))
#             af_prov <- 0
#         } else {
#             af_prov <- 1
#         }
#     } else {
#         stop("Summary statistics must be provided in dataframe")
#     }
#
#     # prep LD for gbayes
#     LD_values <- lapply(1:nrow(LD), function(i) as.numeric(LD[i, ]))
#     names(LD_values) <- variant_ids
#
#     LD_indices <- list(indices = vector("list", length = nrow(LD)))
#     for (i in 1:nrow(LD)) {
#         LD_indices[[i]] <- 1:nrow(LD) - 1
#     }
#
#     bm <- dm <- fit <- res <- vector(length = 1, mode = "list")
#     names(bm) <- names(dm) <- names(fit) <- names(res) <- 1
#
#     # Set parameters if not otherwise specified
#     if (is.null(m)) {
#         m <- length(LD_values)
#     }
#     vy <- yy / (n - 1)
#     if (is.null(pi)) {
#         pi <- 0.001
#     }
#     if (is.null(h2)) {
#         h2 <- 0.5
#     }
#     if (is.null(ve)) {
#         ve <- vy * (1 - h2)
#     }
#     if (is.null(vg)) {
#         vg <- vy * h2
#     }
#     if (method < 4 && is.null(vb)) {
#         vb <- vg / m
#     }
#     if (method >= 4 && is.null(vb)) {
#         vb <- vg / (m * pi)
#     }
#     if (is.null(lambda)) {
#         lambda <- rep(ve / vb, m)
#     }
#     if (method < 4 && is.null(ssb_prior)) {
#         ssb_prior <- ((nub - 2.0) / nub) * (vg / m)
#     }
#     if (method >= 4 && is.null(ssb_prior)) {
#         ssb_prior <- ((nub - 2.0) / nub) * (vg / (m * pi))
#     }
#     if (is.null(sse_prior)) {
#         sse_prior <- ((nue - 2.0) / nue) * ve
#     }
#     if (is.null(b)) {
#         b <- rep(0, m)
#     }
#
#     pi <- c(1 - pi, pi)
#     gamma <- c(0, 1.0)
#     if (method == 5) {
#         pi <- c(0.95, 0.02, 0.02, 0.01)
#     }
#     if (method == 5) {
#         gamma <- c(0, 0.01, 0.1, 1.0)
#     }
#
#     seed <- sample.int(.Machine$integer.max, 1)
#
#     fit <- qgg:::sbayes_spa(
#         wy = wy,
#         ww = ww,
#         LDvalues = LD_values,
#         LDindices = LD_indices,
#         b = b,
#         lambda = lambda,
#         mask = mask,
#         yy = yy,
#         pi = pi,
#         gamma = gamma,
#         vg = vg,
#         vb = vb,
#         ve = ve,
#         ssb_prior = ssb_prior,
#         sse_prior = sse_prior,
#         nub = nub,
#         nue = nue,
#         updateB = updateB,
#         updateE = updateE,
#         updatePi = updatePi,
#         updateG = updateG,
#         adjustE = adjustE,
#         n = n,
#         nit = nit,
#         nburn = nburn,
#         nthin = nthin,
#         algo = algo,
#         method = as.integer(method),
#         seed = seed
#     )
#
#     names(fit[[1]]) <- names(LD_values)
#     names(fit) <- c(
#         "bm",
#         "dm",
#         "coef",
#         "vbs",
#         "vgs",
#         "ves",
#         "pis",
#         "pim",
#         "r",
#         "b",
#         "param"
#     )
#     fit[3] <- NULL
#
#     res <- data.frame(
#         variant_id = variant_ids,
#         bm = fit$bm,
#         dm = fit$dm,
#         pos = sumstats$pos,
#         A1 = sumstats$A1,
#         A2 = sumstats$A2,
#         maf = sumstats$maf,
#         stringsAsFactors = FALSE
#     )
#     rownames(res) <- variant_ids
#
#     fit$sumstats <- res
#     if (af_prov == 1) {
#         fit$sumstats$vm <- 2 *
#             (1 - fit$sumstats$maf) *
#             fit$sumstats$maf *
#             fit$sumstats$bm^2
#     }
#     fit$method <- methods[method]
#     fit$mask <- mask
#
#     zve <- coda::geweke.diag(fit$ves[nburn:length(fit$ves)])$z
#     zvg <- coda::geweke.diag(fit$vgs[nburn:length(fit$vgs)])$z
#     zvb <- coda::geweke.diag(fit$vbs[nburn:length(fit$vbs)])$z
#     zpi <- coda::geweke.diag(fit$pis[nburn:length(fit$pis)])$z
#
#     ve <- mean(fit$ves[nburn:length(fit$ves)])
#     vg <- mean(fit$vgs[nburn:length(fit$vgs)])
#     vb <- mean(fit$vbs[nburn:length(fit$vbs)])
#     pi <- 1 - fit$pim[1]
#     fit$conv <- data.frame(zve = zve, zvg = zvg, zvb = zvb, zpi = zpi)
#     fit$post <- data.frame(ve = ve, vg = vg, vb = vb, pi = pi)
#     fit$ve <- mean(ve)
#     fit$vg <- sum(vg)
#
#     return(fit)
# }
# #' Extract weights from gbayes_rss function
# #' @return A numeric vector of the posterior mean of the coefficients.
# #' @export
# bayes_alphabet_rss_weights <- function(sumstats, LD, method, ...) {
#     model <- gbayes_rss(sumstats = sumstats, LD = LD, method = method, ...)
#     return(model$bm)
# }
# #' Use Gaussian distribution as prior. Posterior means will be BLUP,
# #' equivalent to Ridge Regression.
# #' @export
# bayes_n_rss_weights <- function(sumstats, LD, methodArgs = list()) {
#     return(bayes_alphabet_rss_weights(sumstats, LD, method = "bayesN", ...))
# }
# #' Use laplace/double exponential distribution as prior. This is equivalent
# #' to Bayesian LASSO.
# #' @export
# bayes_l_rss_weights <- function(sumstats, LD, methodArgs = list()) {
#     return(bayes_alphabet_rss_weights(sumstats, LD, method = "bayesL", ...))
# }
# #' Use t-distribution as prior.
# #' @export
# bayes_a_rss_weights <- function(sumstats, LD, methodArgs = list()) {
#     return(bayes_alphabet_rss_weights(sumstats, LD, method = "bayesA", ...))
# }
# #' Use a rounded spike prior (low-variance Gaussian).
# #' @export
# bayes_c_rss_weights <- function(sumstats, LD, methodArgs = list()) {
#     return(bayes_alphabet_rss_weights(sumstats, LD, method = "bayesC", ...))
# }
# #' Use a hierarchical Bayesian mixture model with four Gaussian components.
# #' Variances are scaled
# #' by 0, 0.0001 , 0.001 , and 0.01 .
# #' @export
# bayes_r_rss_weights <- function(sumstats, LD, methodArgs = list()) {
#     return(bayes_alphabet_rss_weights(sumstats, LD, method = "bayesR", ...))
# }

#' Lassosum RSS: LASSO on summary statistics with LD reference
#'
#' Coordinate descent to solve the penalized regression on summary statistics:
#' \deqn{f(\beta) = \beta' R \beta - 2\beta' r + 2\lambda ||\beta||_1}
#' where \eqn{R} is the LD matrix (pre-shrunk if desired) and \eqn{r =
#' \hat\beta / \sqrt{n}}.
#'
#' Based on Mak et al (2017) "Polygenic scores via penalized regression on
#' summary statistics", Genetic Epidemiology 41(6):469-480.
#'
#' @param bhat A vector of marginal effect sizes.
#' @param R The LD correlation matrix (a single matrix over the analysed
#'   window), as in \code{susieR::susie_rss()}. If shrinkage is desired, apply
#'   it before passing (e.g., \code{(1-s)*R + s*I}).
#' @param n Sample size of the GWAS.
#' @param lambda A vector of L1 penalty values. Default: 20 values from 0.001 to
#'   0.1 on log scale.
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
#' out <- lassosumRss(bhat, R, n)
#' @export
#' @importFrom checkmate assertCount assertNumber assertNumeric
lassosumRss <- function(
    bhat,
    R,
    n,
    lambda = exp(seq(log(0.0001), log(0.1), length.out = 20)),
    thr = 1e-4,
    maxiter = 10000
) {
    assertNumeric(lambda, lower = 0, any.missing = FALSE)
    assertNumber(thr, lower = 0, finite = TRUE)
    assertCount(maxiter, positive = TRUE)
    # cpp11 requires exact integer types; the C++ backend takes a block list, so
    # the single-window matrix R is wrapped as one block here.
    .rssSolvePath(
        bhat,
        R,
        n,
        lambda,
        .rssLassosumSolve,
        solveArgs = list(thr = thr, maxiter = maxiter)
    )
}

.lassosumCorFromStat <- function(stat, n, p) {
    corInput <- if (!is.null(stat$cor)) {
        as.numeric(stat$cor)
    } else if (!is.null(stat$z)) {
        as.numeric(stat$z) / sqrt(n)
    } else if (!is.null(stat$b)) {
        as.numeric(stat$b)
    } else {
        msg <- glue(
            "stat must contain one of 'cor', 'z', or 'b' for lassosum ",
            "selection."
        )
        abort(msg)
    }
    if (length(corInput) != p) {
        nInput <- length(corInput)
        msg <- glue(
            "The length of lassosum input statistics ({nInput}) must ",
            "equal nrow(LD) ({p})."
        )
        abort(msg)
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
        mode = "minFbeta"
    )
}

.lassosumSelectLdQuadratic <- function(candidateBeta, corInput, LD) {
    ldBeta <- LD %*% candidateBeta
    bxy <- as.numeric(crossprod(corInput, candidateBeta))
    bxxb <- colSums(candidateBeta * ldBeta)
    positive <- is.finite(bxxb) & bxxb > 0
    scores <- replace(
        rep(-Inf, length(bxy)),
        positive,
        bxy[positive] / sqrt(bxxb[positive])
    )
    idx <- .lassosumFirstMax(scores)
    list(
        beta = candidateBeta[, idx],
        index = idx,
        mode = "ldQuadratic"
    )
}

# Validate the (bhat, R, n) inputs shared by the RSS solvers (lassosumRss /
# penalizedRss via .rssSolvePath, prsCs, and sdpr). R is a single LD correlation
# matrix over one cis-window, matching susieR::susie_rss; bhat must match
# nrow(R). missing() guards turn an omitted R / n into a clear message instead
# of an "argument is missing" error, and propagate through the public wrappers
# (verified two levels deep). Method-specific checks -- prsCs's maf length,
# sdpr's M / perVariantSampleSize / array -- stay in the caller.
#' @importFrom checkmate assertVector
.rssValidateInputs <- function(bhat, R, n) {
    if (missing(R) || !is.matrix(R)) {
        abort("Please provide the LD correlation matrix 'R' as a matrix.")
    }
    if (missing(n) || n <= 0) {
        abort("Please provide a valid sample size using 'n'.")
    }
    assertVector(bhat, len = nrow(R))
    invisible(NULL)
}

# Validate (bhat, LD, n), run a decreasing-lambda coordinate-descent sweep via
# `solveFn` (the penalty-specific Rcpp backend), then reorder back to the input
# lambda order via the inverse permutation and assemble the standard result
# list. Shared by lassosumRss and penalizedRss, which differ only in which Rcpp
# solver they pass as `solveFn` (and penalizedRss's per-penalty gamma default).
.rssSolvePath <- function(bhat, R, n, lambda, solveFn, solveArgs = list()) {
    .rssValidateInputs(bhat, R, n)
    z <- bhat / sqrt(n)
    order <- order(lambda, decreasing = TRUE)
    # `solveArgs` rather than `...`: the two solvers take different fixed
    # argument sets, and both callers know theirs statically, so an unknown
    # solver argument should be an error here rather than reaching the solver.
    solved <- exec(solveFn, z, lambda[order], R, !!!solveArgs)
    # Reorder back to original lambda order via the inverse permutation.
    invOrder <- order(order)
    beta <- solved$beta[, invOrder, drop = FALSE]
    fbeta <- solved$fbeta[invOrder]
    list_assign(
        solved,
        beta = beta,
        conv = solved$conv[invOrder],
        loss = solved$loss[invOrder],
        fbeta = fbeta,
        lambda = lambda,
        nparams = as.integer(colSums(beta != 0)),
        betaEst = as.numeric(beta[, which.min(fbeta)])
    )
}

# Per-`s` fit for one RSS method (`method` selects the solver + which `config`
# fields apply). Returns list(beta, meta); l0learn's inner lambda0 sweep is
# here.
# @noRd
.rssFitOne <- function(method, solverInput, LDs, n, sVal, config) {
    switch(
        method,
        lassosum = .rssFitLassosum(solverInput, LDs, n, sVal, config),
        penalized = .rssFitPenalized(solverInput, LDs, n, sVal, config),
        l0learn = .rssFitL0learn(solverInput, LDs, n, sVal, config)
    )
}

# Shared (s, lambda, fbeta) path-metadata frame for a single-path model.
# @noRd
.rssPathMeta <- function(sVal, model) {
    tibble(
        s = rep(sVal, length(model$lambda)),
        lambda = model$lambda,
        fbeta = model$fbeta
    )
}

# lassosum shrinkage path for one s.
# @noRd
.rssFitLassosum <- function(solverInput, LDs, n, sVal, config) {
    lsArgs <- c(
        list(bhat = solverInput, R = LDs, n = n),
        as.list(config$dotArgs)
    )
    model <- exec(lassosumRss, !!!lsArgs)
    list(beta = model$beta, meta = .rssPathMeta(sVal, model))
}

# ncvreg-penalized shrinkage path for one s.
# @noRd
.rssFitPenalized <- function(solverInput, LDs, n, sVal, config) {
    penArgs <- c(
        list(
            bhat = solverInput,
            R = LDs,
            n = n,
            penalty = config$penalty,
            gamma = config$gamma,
            alpha = config$alpha,
            lambda0 = config$lambda0,
            lambda2 = config$lambda2
        ),
        as.list(config$dotArgs)
    )
    model <- exec(penalizedRss, !!!penArgs)
    list(beta = model$beta, meta = .rssPathMeta(sVal, model))
}

# L0Learn path for one s: one penalizedRss fit per lambda0, column-bound.
# @noRd
.rssFitL0learn <- function(solverInput, LDs, n, sVal, config) {
    fits <- map(
        config$lambda0,
        .rssL0LambdaFit,
        solverInput = solverInput,
        LDs = LDs,
        n = n,
        sVal = sVal,
        config = config
    )
    betaList <- map(fits, "beta")
    metaList <- map(fits, "meta")
    list(
        beta = exec(cbind, !!!betaList),
        meta = bind_rows(metaList)
    )
}

# Attach the method-specific selection attribute to the chosen coefficient
# vector.
# @noRd
.rssFinalize <- function(method, bestBeta, sel, meta, config) {
    base <- c(mode = sel$mode, index = sel$index)
    selection <- switch(
        method,
        lassosum = c(
            base,
            s = meta$s[sel$index],
            lambda = meta$lambda[sel$index]
        ),
        penalized = c(
            base,
            penalty = config$penalty,
            s = meta$s[sel$index],
            lambda = meta$lambda[sel$index]
        ),
        l0learn = c(
            base,
            penalty = config$penalty,
            s = meta$s[sel$index],
            lambda0 = meta$lambda0[sel$index],
            lambda = meta$lambda[sel$index]
        )
    )
    attrName <- if (method == "lassosum") {
        "lassosum_selection"
    } else {
        "penalized_rss_selection"
    }
    `attr<-`(bestBeta, attrName, selection)
}

# Shared scaffold for the RSS shrinkage-grid weight functions
# (lassosumRssWeights / .penalizedRssWeights / l0learnRssWeights). Standardizes
# the stat -> solverInput conversion, the outer LD-shrinkage grid over `s`, the
# candidate accumulation, and the ldQuadratic / minFbeta selection. `method` +
# `config` pick the per-`s` solver (.rssFitOne) and the finalizer
# (.rssFinalize).
# One shrinkage level's fit against the correspondingly shrunk LD.
# @noRd
.rssFitAtS <- function(sVal, method, solverInput, LD, n, p, config) {
    .rssFitOne(
        method,
        solverInput,
        (1 - sVal) * LD + sVal * diag(p),
        n,
        sVal,
        config
    )
}

.rssShrinkGridWeights <- function(
    stat,
    LD,
    s,
    method,
    config,
    selection = c("ldQuadratic", "minFbeta")
) {
    selection <- arg_match(selection)
    n <- median(stat$n)
    p <- nrow(LD)
    corInput <- .lassosumClampCor(.lassosumCorFromStat(stat, n = n, p = p))
    solverInput <- corInput * sqrt(n)
    fits <- map(
        s,
        .rssFitAtS,
        method = method,
        solverInput = solverInput,
        LD = LD,
        n = n,
        p = p,
        config = config
    )
    candidateBeta <- exec(cbind, !!!map(fits, "beta"))
    candidateMeta <- bind_rows(map(fits, "meta"))
    selectorResult <- if (selection == "ldQuadratic") {
        .lassosumSelectLdQuadratic(candidateBeta, corInput, LD)
    } else {
        .lassosumSelectMinFbeta(candidateBeta, candidateMeta)
    }
    bestBeta <- as.numeric(selectorResult$beta)
    .rssFinalize(method, bestBeta, selectorResult, candidateMeta, config)
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
#' @param stat A list with \code{$b} (effect sizes) and \code{$n} (per-variant
#'   sample sizes).
#' @param LD LD correlation matrix R (single matrix, NOT pre-shrunk).
#' @param s Numeric vector of shrinkage parameters to search over. Default:
#'   \code{c(0.2, 0.5, 0.9, 1.0)} following Mak et al (2017) and OTTERS.
#' @param selection Selection strategy. Default \code{"ldQuadratic"} uses
#'   \eqn{c^T \beta / \sqrt{\beta^T R \beta}} on the supplied LD matrix.
#'   \code{"minFbeta"} is retained as an explicit alternative for debugging.
#' @param methodArgs Options forwarded to \code{lassosumRss}, built
#'   with \code{\link{lassosumConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of the posterior SNP coefficients at the best (s,
#'   lambda).
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ss <- lapply(
#'   seq_len(ncol(X)), function(j) coef(summary(lm(y ~ X[, j])))[2, 1:2])
#' stat <- list(
#'   b = vapply(ss, `[`, numeric(1), 1L),
#'   seb = vapply(ss, `[`, numeric(1), 2L),
#'   n = rep(nrow(X), ncol(X))
#' )
#' LD <- cor(X)
#' lassosumRssWeights(stat, LD)
#' @export
lassosumRssWeights <- function(
    stat,
    LD,
    s = c(0.2, 0.5, 0.9, 1.0),
    selection = c("ldQuadratic", "minFbeta"),
    methodArgs = lassosumConfig()
) {
    .assertMethodConfig(methodArgs, "lassosumConfig", "methodArgs")
    selection <- arg_match(selection)
    .rssShrinkGridWeights(
        stat,
        LD,
        s,
        "lassosum",
        list(dotArgs = methodArgs),
        selection
    )
}

#' Penalized Regression on RSS (Summary Statistics) Objective
#'
#' Generalizes \code{lassosumRss()} to support LASSO, MCP, SCAD, L0, L0L1,
#' and L0L2 penalties.  Uses coordinate descent on the objective
#' \deqn{\beta^T R \beta - 2 \beta^T z + \mathrm{penalty}(\beta)}
#' where \eqn{R} is a (possibly pre-shrunk) LD matrix and \eqn{z = \hat\beta /
#' \sqrt{n}}.
#'
#' @param bhat Numeric vector of marginal effect estimates (length p).
#' @param R The LD correlation matrix (a single matrix over the analysed
#'   window), as in \code{lassosumRss()}.
#' @param n GWAS sample size (positive scalar).
#' @param penalty Penalty type: \code{"lasso"}, \code{"MCP"}, \code{"SCAD"},
#'   \code{"L0"}, \code{"L0L1"}, or \code{"L0L2"}.
#' @param lambda Numeric vector of regularization parameter values along which
#'   to trace a solution path (warm-started, largest-first). For LASSO/MCP/SCAD
#'   this is the primary penalty strength; for L0 variants it controls the L1
#'   component.
#' @param gamma Concavity parameter for MCP (default 3) or SCAD (default 3.7).
#'   Ignored for LASSO and L0 variants.
#' @param alpha Elastic-net mixing for MCP/SCAD: \eqn{l_1 = \lambda \alpha},
#'   \eqn{l_2 = \lambda (1-\alpha)}. Default 1 (pure L1, no ridge).
#' @param lambda0 L0 penalty weight (number of non-zeros). Required for L0
#'   variants; ignored otherwise. Default 0.
#' @param lambda2 L2 penalty weight for L0L2 variant. Default 0.
#' @param thr Convergence threshold. Default 1e-4.
#' @param maxiter Maximum coordinate descent iterations per lambda. Default
#'   10000.
#' @param maxSwaps Maximum swap rounds for L0 variants. Default 100. Set to 0 to
#'   disable swaps.
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
#' set.seed(42)
#' p <- 10; n <- 100
#' bhat <- rnorm(p, sd = 0.1)
#' R <- diag(p)
#' # MCP
#' penalizedRss(bhat, R, n, penalty = "MCP")
#' # SCAD
#' penalizedRss(bhat, R, n, penalty = "SCAD")
#' # L0
#' penalizedRss(bhat, R, n, penalty = "L0", lambda0 = 0.01,
#'               lambda = c(0))
#' @export
penalizedRss <- function(
    bhat,
    R,
    n,
    penalty = c("lasso", "MCP", "SCAD", "L0", "L0L1", "L0L2"),
    lambda = exp(seq(log(0.0001), log(0.1), length.out = 20)),
    gamma = NULL,
    alpha = 1.0,
    lambda0 = 0,
    lambda2 = 0,
    thr = 1e-4,
    maxiter = 10000,
    maxSwaps = 100
) {
    penalty <- arg_match(penalty)
    # Default gamma per penalty
    if (is.null(gamma)) {
        gamma <- switch(penalty, SCAD = 3.7, MCP = 3.0, 0.0)
    }
    # C++ backend takes a block list; wrap the single-window matrix R as one
    # block.
    .rssSolvePath(
        bhat,
        R,
        n,
        lambda,
        .rssPenalizedSolve,
        solveArgs = list(
            penalty = penalty,
            gamma = gamma,
            alpha = alpha,
            lambda0 = lambda0,
            lambda2 = lambda2,
            thr = thr,
            maxiter = maxiter,
            maxSwaps = maxSwaps
        )
    )
}

#' RSS Weights Helper for Penalized Methods
#'
#' Shared implementation for \code{scadRssWeights()} and \code{mcpRssWeights()}.
#' (\code{l0learnRssWeights()} uses the lower-level \code{.rssShrinkGridWeights}
#' scaffold directly, since it sweeps an additional \code{lambda0} path.)
#' Searches over a shrinkage grid \code{s} (LD matrix shrinkage \code{(1-s)R +
#' sI}) and selects the best candidate via LD-quadratic pseudovalidation or
#' minimum penalized objective.
#'
#' @param stat,LD,s,selection,penalty,gamma,alpha,lambda0,lambda2 See the
#'   public wrappers for details.
#' @param methodArgs Options forwarded to \code{penalizedRss}, built with
#'   \code{\link{penalizedRssConfig}}.
#' @return Numeric weight vector of length \code{nrow(LD)}.
#' @keywords internal
.penalizedRssWeights <- function(
    stat,
    LD,
    penalty,
    s = c(0.2, 0.5, 0.9, 1.0),
    gamma = NULL,
    alpha = 1.0,
    lambda0 = 0,
    lambda2 = 0,
    selection = c("ldQuadratic", "minFbeta"),
    methodArgs = penalizedRssConfig()
) {
    .assertMethodConfig(methodArgs, "penalizedRssConfig", "methodArgs")
    selection <- arg_match(selection)
    .rssShrinkGridWeights(
        stat,
        LD,
        s,
        "penalized",
        list(
            penalty = penalty,
            gamma = gamma,
            alpha = alpha,
            lambda0 = lambda0,
            lambda2 = lambda2,
            dotArgs = methodArgs
        ),
        selection
    )
}

#' Compute SCAD-Penalized Weights from Summary Statistics
#'
#' Fits SCAD-penalized regression on the RSS objective, searching over a
#' shrinkage grid \code{s} and lambda path. Model selection uses LD-quadratic
#' pseudovalidation by default.
#'
#' @param stat A list with \code{$b} (effect sizes) and \code{$n} (per-variant
#'   sample sizes).
#' @param LD LD correlation matrix R (single matrix, NOT pre-shrunk).
#' @param s Numeric vector of LD shrinkage parameters. Default: \code{c(0.2,
#'   0.5, 0.9, 1.0)}.
#' @param gamma SCAD concavity parameter. Default 3.7.
#' @param alpha Elastic-net mixing (1 = pure L1). Default 1.
#' @param selection Selection strategy: \code{"ldQuadratic"} (default) or
#'   \code{"minFbeta"}.
#' @param methodArgs Options forwarded to \code{penalizedRss}, built
#'   with \code{\link{penalizedRssConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of SNP coefficient weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ss <- lapply(
#'   seq_len(ncol(X)), function(j) coef(summary(lm(y ~ X[, j])))[2, 1:2])
#' stat <- list(
#'   b = vapply(ss, `[`, numeric(1), 1L),
#'   seb = vapply(ss, `[`, numeric(1), 2L),
#'   n = rep(nrow(X), ncol(X))
#' )
#' LD <- cor(X)
#' scadRssWeights(stat, LD)
#' @importFrom checkmate assertList assertNumeric assertNumber
#' @export
scadRssWeights <- function(
    stat,
    LD,
    s = c(0.2, 0.5, 0.9, 1.0),
    gamma = 3.7,
    alpha = 1.0,
    selection = c("ldQuadratic", "minFbeta"),
    methodArgs = penalizedRssConfig()
) {
    .assertMethodConfig(methodArgs, "penalizedRssConfig", "methodArgs")
    assertList(stat)
    assertNumeric(s, lower = 0, any.missing = FALSE)
    assertNumber(gamma, finite = TRUE)
    assertNumber(alpha, finite = TRUE)
    .penalizedRssWeights(
        stat = stat,
        LD = LD,
        penalty = "SCAD",
        s = s,
        gamma = gamma,
        alpha = alpha,
        selection = selection,
        methodArgs = methodArgs
    )
}

#' Compute MCP-Penalized Weights from Summary Statistics
#'
#' Fits MCP-penalized regression on the RSS objective, searching over a
#' shrinkage grid \code{s} and lambda path. Model selection uses LD-quadratic
#' pseudovalidation by default.
#'
#' @param stat A list with \code{$b} (effect sizes) and \code{$n} (per-variant
#'   sample sizes).
#' @param LD LD correlation matrix R (single matrix, NOT pre-shrunk).
#' @param s Numeric vector of LD shrinkage parameters. Default: \code{c(0.2,
#'   0.5, 0.9, 1.0)}.
#' @param gamma MCP concavity parameter. Default 3.
#' @param alpha Elastic-net mixing (1 = pure L1). Default 1.
#' @param selection Selection strategy: \code{"ldQuadratic"} (default) or
#'   \code{"minFbeta"}.
#' @param methodArgs Options forwarded to \code{penalizedRss}, built
#'   with \code{\link{penalizedRssConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of SNP coefficient weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ss <- lapply(
#'   seq_len(ncol(X)), function(j) coef(summary(lm(y ~ X[, j])))[2, 1:2])
#' stat <- list(
#'   b = vapply(ss, `[`, numeric(1), 1L),
#'   seb = vapply(ss, `[`, numeric(1), 2L),
#'   n = rep(nrow(X), ncol(X))
#' )
#' LD <- cor(X)
#' mcpRssWeights(stat, LD)
#' @importFrom checkmate assertList assertNumeric assertNumber
#' @export
mcpRssWeights <- function(
    stat,
    LD,
    s = c(0.2, 0.5, 0.9, 1.0),
    gamma = 3.0,
    alpha = 1.0,
    selection = c("ldQuadratic", "minFbeta"),
    methodArgs = penalizedRssConfig()
) {
    .assertMethodConfig(methodArgs, "penalizedRssConfig", "methodArgs")
    assertList(stat)
    assertNumeric(s, lower = 0, any.missing = FALSE)
    assertNumber(gamma, finite = TRUE)
    assertNumber(alpha, finite = TRUE)
    .penalizedRssWeights(
        stat = stat,
        LD = LD,
        penalty = "MCP",
        s = s,
        gamma = gamma,
        alpha = alpha,
        selection = selection,
        methodArgs = methodArgs
    )
}

#' Compute L0-Penalized Weights from Summary Statistics
#'
#' Fits L0-penalized regression (with optional L1/L2 components) on the RSS
#' objective, searching over a shrinkage grid \code{s} and lambda0 path. Model
#' selection uses LD-quadratic pseudovalidation by default.
#'
#' The swap optimization from L0Learn is included: after coordinate descent
#' converges, non-zero coefficients are tested for swaps with zero ones to
#' escape local optima.
#'
#' @param stat A list with \code{$b} (effect sizes) and \code{$n} (per-variant
#'   sample sizes).
#' @param LD LD correlation matrix R (single matrix, NOT pre-shrunk).
#' @param penalty L0 variant: \code{"L0"}, \code{"L0L1"}, or \code{"L0L2"}.
#'   Default \code{"L0"}.
#' @param s Numeric vector of LD shrinkage parameters. Default: \code{c(0.2,
#'   0.5, 0.9, 1.0)}.
#' @param lambda0 Numeric vector of L0 penalty values to search over. Default:
#'   \code{exp(seq(log(0.001), log(1), length.out = 10))}.
#' @param lambda Numeric vector of L1 penalty values (for L0L1). Default:
#'   \code{c(0)} (no L1 unless L0L1 is used).
#' @param lambda2 L2 penalty weight (for L0L2). Default 0.
#' @param selection Selection strategy: \code{"ldQuadratic"} (default) or
#'   \code{"minFbeta"}.
#' @param maxSwaps Maximum swap rounds per lambda. Default 100.
#' @param methodArgs Options forwarded to \code{penalizedRss}, built
#'   with \code{\link{penalizedRssConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of SNP coefficient weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ss <- lapply(
#'   seq_len(ncol(X)), function(j) coef(summary(lm(y ~ X[, j])))[2, 1:2])
#' stat <- list(
#'   b = vapply(ss, `[`, numeric(1), 1L),
#'   seb = vapply(ss, `[`, numeric(1), 2L),
#'   n = rep(nrow(X), ncol(X))
#' )
#' LD <- cor(X)
#' l0learnRssWeights(stat, LD)
#' @export
l0learnRssWeights <- function(
    stat,
    LD,
    penalty = c("L0", "L0L1", "L0L2"),
    s = c(0.2, 0.5, 0.9, 1.0),
    lambda0 = exp(seq(log(0.001), log(1), length.out = 10)),
    lambda = NULL,
    lambda2 = 0,
    selection = c("ldQuadratic", "minFbeta"),
    maxSwaps = 100,
    methodArgs = penalizedRssConfig()
) {
    .assertMethodConfig(methodArgs, "penalizedRssConfig", "methodArgs")
    penalty <- arg_match(penalty)
    selection <- arg_match(selection)

    # Default lambda (L1 component) depends on variant
    if (is.null(lambda)) {
        lambda <- if (penalty == "L0L1") {
            exp(seq(log(0.0001), log(0.1), length.out = 10))
        } else {
            c(0)
        }
    }

    .rssShrinkGridWeights(
        stat,
        LD,
        s,
        "l0learn",
        list(
            penalty = penalty,
            lambda = lambda,
            lambda0 = lambda0,
            lambda2 = lambda2,
            maxSwaps = maxSwaps,
            dotArgs = methodArgs
        ),
        selection
    )
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
#' @param methodArgs Options forwarded to \code{ncvreg::cv.ncvreg}, built
#'   with \code{\link{ncvregConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of length `ncol(X)` of variant weights.
#' @importFrom stats coef
#' @keywords internal
#' @importFrom checkmate assertCount
ncvregWeights <- function(
    X,
    y,
    penalty,
    nfolds = 5,
    methodArgs = ncvregConfig()
) {
    .assertMethodConfig(methodArgs, "ncvregConfig", "methodArgs")
    assertCount(nfolds, positive = TRUE)
    if (!requireNamespace("ncvreg", quietly = TRUE)) {
        abort("Package 'ncvreg' is required for this function.")
    }
    eff.wgt <- matrix(0, ncol = 1, nrow = ncol(X))
    keep <- .dropZeroVariance(X, "ncvregWeights")
    callArgs <- list_modify(
        list(
            X = X[, keep, drop = FALSE],
            y = y,
            penalty = penalty,
            nfolds = nfolds
        ),
        !!!methodArgs
    )
    fit <- exec(ncvreg::cv.ncvreg, !!!callArgs)
    replace(eff.wgt, keep, coef(fit, lambda = fit$lambda.min)[-1])
}

#' Compute Weights Using SCAD-Penalized Regression
#'
#' Fits a SCAD-penalized linear regression model via `ncvreg::cv.ncvreg` and
#' returns the coefficient vector at `lambda.min`.
#'
#' @param X A numeric matrix of predictors.
#' @param y A numeric response vector.
#' @param nfolds Number of cross-validation folds. Default is 5.
#' @param methodArgs Options forwarded to \code{ncvreg::cv.ncvreg}, built
#'   with \code{\link{ncvregConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of length `ncol(X)` of variant weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' scadWeights(X, y)
#' @export
scadWeights <- function(X, y, nfolds = 5, methodArgs = ncvregConfig()) {
    .assertMethodConfig(methodArgs, "ncvregConfig", "methodArgs")
    ncvregWeights(
        X,
        y,
        penalty = "SCAD",
        nfolds = nfolds,
        methodArgs = methodArgs
    )
}

#' Compute Weights Using MCP-Penalized Regression
#'
#' Fits an MCP-penalized linear regression model via `ncvreg::cv.ncvreg` and
#' returns the coefficient vector at `lambda.min`.
#'
#' @param X A numeric matrix of predictors.
#' @param y A numeric response vector.
#' @param nfolds Number of cross-validation folds. Default is 5.
#' @param methodArgs Options forwarded to \code{ncvreg::cv.ncvreg}, built
#'   with \code{\link{ncvregConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of length `ncol(X)` of variant weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' mcpWeights(X, y)
#' @export
mcpWeights <- function(X, y, nfolds = 5, methodArgs = ncvregConfig()) {
    .assertMethodConfig(methodArgs, "ncvregConfig", "methodArgs")
    ncvregWeights(
        X,
        y,
        penalty = "MCP",
        nfolds = nfolds,
        methodArgs = methodArgs
    )
}

#' Compute Weights Using L0Learn
#'
#' Fits an L0-regularized linear regression model via `L0Learn::L0Learn.cvfit`
#' and returns the coefficient vector at the (lambda, gamma) pair minimizing the
#' cross-validation error. Default penalty is "L0"; the user can switch to
#' "L0L1" or "L0L2" (and tune the corresponding gamma grid) by passing the
#' relevant arguments through `...`.
#'
#' @param X A numeric matrix of predictors.
#' @param y A numeric response vector.
#' @param penalty Type of regularization: "L0", "L0L1", or "L0L2". Default is
#'   "L0".
#' @param nFolds Number of cross-validation folds. Default is 5.
#' @param methodArgs Options forwarded to \code{L0Learn::L0Learn.cvfit}, built
#'   with \code{\link{l0learnConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of length `ncol(X)` of variant weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' l0learnWeights(X, y)
#' @export
l0learnWeights <- function(
    X,
    y,
    penalty = "L0",
    nFolds = 5,
    methodArgs = l0learnConfig()
) {
    .assertMethodConfig(methodArgs, "l0learnConfig", "methodArgs")
    if (!requireNamespace("L0Learn", quietly = TRUE)) {
        abort("Package 'L0Learn' is required for this function.")
    }
    eff.wgt <- matrix(0, ncol = 1, nrow = ncol(X))
    keep <- .dropZeroVariance(X, "l0learnWeights")
    callArgs <- list_modify(
        list(
            x = X[, keep, drop = FALSE],
            y = y,
            penalty = penalty,
            nFolds = nFolds
        ),
        !!!methodArgs
    )
    fit <- exec(L0Learn::L0Learn.cvfit, !!!callArgs)
    # Find (gamma, lambda) minimizing CV error across the entire path.
    cvMins <- map_dbl(fit$cvMeans, .rssMinNumeric)
    gammaIdx <- which.min(cvMins)
    lambdaIdx <- which.min(as.numeric(fit$cvMeans[[gammaIdx]]))
    bestGamma <- fit$fit$gamma[gammaIdx]
    bestLambda <- fit$fit$lambda[[gammaIdx]][lambdaIdx]
    raw <- as.numeric(coef(fit, lambda = bestLambda, gamma = bestGamma))
    # If intercept was included, drop it (first row).
    coefs <- if (length(raw) == sum(keep) + 1L) raw[-1L] else raw
    replace(eff.wgt, keep, coefs)
}

#' Compute Weights Using a BGLR Linear Regression Model
#'
#' Internal helper that fits a `BGLR::BGLR` linear regression with a single
#' linear term whose `model` is one of BGLR's marker-effect priors (e.g.
#' "BayesB", "BL"), then returns the posterior mean of the marker effects. BGLR
#' writes per-call temporary files to disk; this helper sandboxes them in a
#' fresh `tempdir()` that is cleaned up on exit.
#'
#' @param X A numeric matrix of predictors.
#' @param y A numeric response vector.
#' @param model A BGLR marker-effect model name (e.g. "BayesB" or "BL").
#' @param nIter Number of MCMC iterations.
#' @param burnIn Number of burn-in iterations.
#' @param thin Thinning interval.
#' @param etaArgs Optional named list of additional arguments included in the
#'   `ETA` linear-term specification (e.g. `list(probIn = 0.05)` for BayesB).
#' @param methodArgs Options forwarded to \code{BGLR::BGLR}, built
#'   with \code{\link{bglrConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of length `ncol(X)` of variant weights.
#' @keywords internal
bglrWeights <- function(
    X,
    y,
    model,
    nIter,
    burnIn,
    thin,
    etaArgs = list(),
    methodArgs = bglrConfig()
) {
    .assertMethodConfig(methodArgs, "bglrConfig", "methodArgs")
    if (!requireNamespace("BGLR", quietly = TRUE)) {
        abort("Package 'BGLR' is required for this function.")
    }
    eff.wgt <- rep(0, ncol(X))
    keep <- .dropZeroVariance(X, "bglrWeights")

    tmpdir <- tempfile("bglr_")
    dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
    saveAt <- str_c(tmpdir, .Platform$file.sep)

    eta <- list(c(list(X = X[, keep, drop = FALSE], model = model), etaArgs))
    callArgs <- list_modify(
        list(
            y = y,
            ETA = eta,
            nIter = nIter,
            burnIn = burnIn,
            thin = thin,
            saveAt = saveAt,
            verbose = FALSE
        ),
        !!!methodArgs
    )
    fit <- exec(BGLR::BGLR, !!!callArgs)
    replace(eff.wgt, keep, as.numeric(fit$ETA[[1]]$b))
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
#' under correlated predictors. Override these arguments to recover the package
#' defaults if desired.
#'
#' @param X A numeric matrix of predictors.
#' @param y A numeric response vector.
#' @param nIter Number of MCMC iterations. Default is 10000.
#' @param burnIn Number of burn-in iterations. Default is 2000.
#' @param thin Thinning interval. Default is 5.
#' @param probIn Prior inclusion probability for each marker. Default is 0.2.
#' @param methodArgs Options forwarded to \code{BGLR::BGLR}, built
#'   with \code{\link{bglrConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of length `ncol(X)` of variant weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' bayesBWeights(X, y)
#' @importFrom checkmate assertCount assertNumber
#' @export
bayesBWeights <- function(
    X,
    y,
    nIter = 10000,
    burnIn = 2000,
    thin = 5,
    probIn = 0.2,
    methodArgs = bglrConfig()
) {
    .assertMethodConfig(methodArgs, "bglrConfig", "methodArgs")
    assertCount(nIter, positive = TRUE)
    assertCount(burnIn)
    assertCount(thin, positive = TRUE)
    assertNumber(probIn, lower = 0, upper = 1)
    bglrWeights(
        X,
        y,
        model = "BayesB",
        nIter = nIter,
        burnIn = burnIn,
        thin = thin,
        etaArgs = list(probIn = probIn),
        methodArgs = methodArgs
    )
}

#' Compute Weights Using the Bayesian LASSO (BGLR)
#'
#' Fits a Bayesian LASSO linear regression model via `BGLR::BGLR` (the "BL"
#' model, Park & Casella 2008) and returns the posterior mean of the marker
#' effects. This is the same "B-Lasso" implementation benchmarked in Kim et al.
#' (2022). Note that this is distinct from `bayesLWeights`, which uses a
#' different Bayesian LASSO implementation backed by `qgg`.
#'
#' Defaults for `nIter`, `burnIn`, and `thin` are larger than BGLR's package
#' defaults to better accommodate high-LD cis-eQTL windows; override to recover
#' the package defaults.
#'
#' @param X A numeric matrix of predictors.
#' @param y A numeric response vector.
#' @param nIter Number of MCMC iterations. Default is 10000.
#' @param burnIn Number of burn-in iterations. Default is 2000.
#' @param thin Thinning interval. Default is 5.
#' @param methodArgs Options forwarded to \code{BGLR::BGLR}, built
#'   with \code{\link{bglrConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of length `ncol(X)` of variant weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' bLassoWeights(X, y)
#' @importFrom checkmate assertCount
#' @export
bLassoWeights <- function(
    X,
    y,
    nIter = 10000,
    burnIn = 2000,
    thin = 5,
    methodArgs = bglrConfig()
) {
    .assertMethodConfig(methodArgs, "bglrConfig", "methodArgs")
    assertCount(nIter, positive = TRUE)
    assertCount(burnIn)
    assertCount(thin, positive = TRUE)
    bglrWeights(
        X,
        y,
        model = "BL",
        nIter = nIter,
        burnIn = burnIn,
        thin = thin,
        methodArgs = methodArgs
    )
}

#' Compute Weights Using Dirichlet Process Regression (RcppDPR)
#'
#' Fits a Dirichlet Process Regression model via `RcppDPR::fit_model` and
#' returns the per-variant weights, computed as `beta + alpha` (matching
#' RcppDPR's internal `predict.DPR_Model`, which uses `(beta + alpha) %*% x_new
#' + pheno_mean`).
#'
#' By default the variational Bayes (`VB`) fitting method is used, which is fast
#' and deterministic. The user may switch to `Gibbs` or `Adaptive_Gibbs` for
#' full Bayesian MCMC inference. `rotate_variables` is held to `FALSE` under the
#' assumption that any covariates have already been regressed out upstream; an
#' intercept-only covariate matrix is supplied to `fit_model`.
#'
#' @param X A numeric matrix of predictors.
#' @param y A numeric response vector.
#' @param fittingMethod One of "VB", "Gibbs", or "Adaptive_Gibbs". Default is
#'   "VB".
#' @param methodArgs Options forwarded to \code{RcppDPR::fit_model}, built
#'   with \code{\link{dprConfig}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights: \code{"none"} (default) attaches
#'   nothing, \code{"slim"} and \code{"full"} attach the fitted model. This
#'   engine keeps no trimmable intermediates, so the two retaining levels
#'   behave alike here.
#'   Default \code{FALSE}.
#' @param nK Integer. Number of variational mixture components for the VB fit
#'   (\code{dprVbWeights}). Default \code{8}.
#' @param sStep Integer. Number of Gibbs sampling steps for the Gibbs fit
#'   (\code{dprGibbsWeights}). Default \code{5000}.
#' @return A numeric vector of length `ncol(X)` of variant weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' dprWeights(X, y)
#' @export
dprWeights <- function(
    X,
    y,
    fittingMethod = "VB",
    fitRetention = c("none", "slim", "full"),
    methodArgs = dprConfig()
) {
    fitRetention <- arg_match(fitRetention)
    .assertMethodConfig(methodArgs, "dprConfig", "methodArgs")
    if (!requireNamespace("RcppDPR", quietly = TRUE)) {
        abort("Package 'RcppDPR' is required for this function.")
    }
    zeros <- rep(0, ncol(X))
    keep <- .dropZeroVariance(X, "dprWeights")
    w <- matrix(1, nrow = nrow(X), ncol = 1)
    callArgs <- list_modify(
        list(
            y = y,
            w = w,
            x = X[, keep, drop = FALSE],
            rotate_variables = FALSE,
            fitting_method = fittingMethod
        ),
        !!!methodArgs
    )
    fit <- exec(RcppDPR::fit_model, !!!callArgs)
    eff.wgt <- replace(zeros, keep, as.numeric(fit$beta + fit$alpha))
    if (identical(fitRetention, "none")) {
        return(eff.wgt)
    }
    `attr<-`(eff.wgt, "fit", fit)
}

#' @rdname dprWeights
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' dprVbWeights(X, y)
#' @export
dprVbWeights <- function(
    X,
    y,
    nK = 8,
    fitRetention = c("none", "slim", "full"),
    methodArgs = dprConfig()
) {
    fitRetention <- arg_match(fitRetention)
    .assertMethodConfig(methodArgs, "dprConfig", "methodArgs")
    dprWeights(
        X,
        y,
        fittingMethod = "VB",
        fitRetention = fitRetention,
        methodArgs = .methodConfigWith(methodArgs, "dprConfig", n_k = nK)
    )
}

#' @rdname dprWeights
#' @examples
#' set.seed(1)
#' n <- 50
#' p <- 8
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- sprintf("chr1:%d:A:G", 100L * (1:p))
#' y <- X[, 1] * 0.5 + rnorm(n)
#' dprGibbsWeights(X, y, sStep = 500)
#' @importFrom checkmate assertCount assertFlag
#' @export
dprGibbsWeights <- function(
    X,
    y,
    sStep = 5000,
    fitRetention = c("none", "slim", "full"),
    methodArgs = dprConfig()
) {
    .assertMethodConfig(methodArgs, "dprConfig", "methodArgs")
    assertCount(sStep, positive = TRUE)
    fitRetention <- arg_match(fitRetention)
    dprWeights(
        X,
        y,
        fittingMethod = "Gibbs",
        fitRetention = fitRetention,
        methodArgs = .methodConfigWith(
            methodArgs,
            "dprConfig",
            s_step = sStep
        )
    )
}

#' @rdname dprWeights
#' @examples
#' set.seed(1)
#' n <- 50
#' p <- 8
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- sprintf("chr1:%d:A:G", 100L * (1:p))
#' y <- X[, 1] * 0.5 + rnorm(n)
#' dprAdaptiveGibbsWeights(X, y)
#' @importFrom checkmate assertFlag
#' @export
dprAdaptiveGibbsWeights <- function(
    X,
    y,
    fitRetention = c("none", "slim", "full"),
    methodArgs = dprConfig()
) {
    .assertMethodConfig(methodArgs, "dprConfig", "methodArgs")
    fitRetention <- arg_match(fitRetention)
    dprWeights(
        X,
        y,
        fittingMethod = "Adaptive_Gibbs",
        fitRetention = fitRetention,
        methodArgs = methodArgs
    )
}
#' @title Mr.Mash Wrapper
#'
#' @description Compute weights with mr.mash using a precomputed prior grid and
#'   mixture prior.
#'
#' @param X An n x p matrix of genotype data, where n is the total number of
#'   individuals and p is the number of SNPs.
#' @param Y An n x r matrix of residual expression data, where n is the total
#'   number of individuals and r is the total number of conditions
#'   (tissue/cell-types).
#' @param dataDrivenPriorMatrices A list of data-driven covariance matrices.
#'   Default is NULL.
#' @param prior Prior-covariance construction options, built with
#'   \code{\link{mrmashPriorConfig}}.
#' @param bInitMethod The method for initializing the coefficient matrix.
#'   Default is "enet".
#' @param updateV A logical indicating whether to update the residual covariance
#'   matrix. Default is TRUE.
#' @param updateVMethod The method for updating the residual covariance matrix.
#'   Default is "full".
#' @param standardize A logical indicating whether to standardize the input
#'   data. Default is FALSE. Used by the summary-statistic and initialisation
#'   steps as well as the fit, so it is set here rather than in
#'   \code{methodArgs}.
#' @param numThreads The number of threads to use for parallel computation.
#'   Default is 1. Used by the same three steps as \code{standardize}.
#' @param methodArgs Options forwarded to \code{mr.mashr::mr.mash} under its
#'   own names, built with \code{\link{mrmashConfig}}
#'   (\code{update_w0}, \code{max_iter}, \code{tol}, \code{w0_threshold},
#'   \code{verbose}, ...). Arguments the wrapper derives --- the prior, the
#'   residual covariance, the coefficient initialisation --- cannot be set
#'   here.
#'
#' @param V Optional residual covariance matrix (conditions x conditions), or
#'   \code{NULL} to estimate it.
#' @param sumstats Optional list of summary statistics for the RSS variant of
#'   mr.mash, or \code{NULL} for the individual-data variant.
#' @return A mr.mash fit, stored as a list with some or all of the following
#' elements:
#' \item{mu1}{A p x r matrix of posterior means for the regression
#' coefficients.}
#' \item{S1}{An r x r x p array of posterior covariances for the regression
#' coefficients.}
#' \item{w1}{A p x K matrix of posterior assignment probabilities to the
#' mixture components.}
#' \item{V}{An r x r residual covariance matrix.}
#' \item{w0}{A K-vector with (updated, if \code{update_w0=TRUE}) prior mixture
#' weights, each associated with the respective covariance matrix in
#' \code{S0}.}
#' \item{S0}{An r x r x K array of prior covariance matrices on the regression
#' coefficients.}
#' \item{intercept}{An r-vector containing the posterior mean estimate of the
#' intercept.}
#' \item{fitted}{An n x r matrix of fitted values.}
#' \item{G}{An r x r covariance matrix of fitted values.}
#' \item{pve}{An r-vector of proportion of variance explained by the
#' covariates.}
#' \item{ELBO}{The Evidence Lower Bound (ELBO) at the last iteration.}
#' \item{progress}{A data frame including information regarding convergence
#' criteria at each iteration.}
#' \item{converged}{A logical indicating whether the optimization algorithm
#' converged to a solution within the chosen tolerance level.}
#' \item{elapsed_time}{The computation runtime for fitting mr.mash.}
#' \item{Y}{An n x r matrix of responses at the last iteration (only relevant
#' when missing values are present in the input Y).}
#'
#' @examples
#' data(multiTraitData)
#' res <- mrmashWrapper(
#'   X = multiTraitData$X[, 1:60], Y = multiTraitData$Y,
#'   dataDrivenPriorMatrices = multiTraitData$priorMatrices,
#'   prior = mrmashPriorConfig(canonicalPriorMatrices = TRUE)
#' )
#'
#' @export
mrmashWrapper <- function(
    X,
    Y,
    V = NULL,
    sumstats = NULL,
    dataDrivenPriorMatrices = NULL,
    prior = mrmashPriorConfig(),
    bInitMethod = "enet",
    updateV = TRUE,
    updateVMethod = "full",
    standardize = FALSE,
    numThreads = 1,
    methodArgs = mrmashConfig()
) {
    .mrmashRequirePackages()
    .assertMethodConfig(prior, "mrmashPriorConfig", "prior")
    .assertMethodConfig(methodArgs, "mrmashConfig", "methodArgs")
    priorOpts <- as.list(prior)
    .mrmashValidateWrapper(
        X = X,
        Y = Y,
        priorGrid = priorOpts$priorGrid,
        dataDrivenPriorMatrices = dataDrivenPriorMatrices,
        canonicalPriorMatrices = priorOpts$canonicalPriorMatrices %||% FALSE
    )
    bInitMethod <- .mrmashResolveBInit(Y, bInitMethod)
    if (is.null(sumstats)) {
        sumstats <- .mrmashComputeSumstats(X, Y, standardize, numThreads)
    }
    # Shared prior-covariance builder (also used by mrmashRssWeights).
    priorBuilt <- exec(
        buildMrmashPriorMatrices,
        Bhat = sumstats$Bhat,
        Shat = sumstats$Shat,
        K = ncol(Y),
        dataDrivenPriorMatrices = dataDrivenPriorMatrices,
        !!!priorOpts
    )
    time1 <- proc.time()
    bInit <- as.matrix(
        .mrmashInitCoefficients(
            X,
            Y,
            bInitMethod,
            standardize,
            numThreads
        )$Bhat
    )
    vInit <- .mrmashInitV(X, Y, V, updateV, updateVMethod)
    fitMrmash <- .mrmashFit(
        priorBuilt$S0,
        bInit,
        vInit,
        X = X,
        Y = Y,
        standardize = standardize,
        numThreads = numThreads,
        updateVMethod = updateVMethod,
        methodArgs = methodArgs
    )
    list_assign(
        fitMrmash,
        analysis_time = proc.time()["elapsed"] - time1["elapsed"]
    )
}

# Require glmnet + mr.mashr; also emit the no-seed reproducibility message.
# @noRd
.mrmashRequirePackages <- function() {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        abort("Package 'glmnet' is required for this function.")
    }
    if (!requireNamespace("mr.mashr", quietly = TRUE)) {
        abort("Package 'mr.mashr' is required for this function.")
    }
}

# Input validation for the individual-level mr.mash wrapper.
# @noRd
#' @importFrom checkmate assertMatrix
.mrmashValidateWrapper <- function(
    X,
    Y,
    priorGrid,
    dataDrivenPriorMatrices,
    canonicalPriorMatrices
) {
    if (!exists(".Random.seed")) {
        inform(
            "! No seed has been set. Please set seed for reproducable result. "
        )
    }
    if (!is.matrix(X) || !is.matrix(Y)) {
        abort("X and Y must be matrices.")
    }
    assertMatrix(Y, nrows = nrow(X))
    if (!is.null(priorGrid) && !is.vector(priorGrid)) {
        abort("priorGrid must be a vector.")
    }
    if (is.null(dataDrivenPriorMatrices) && !isTRUE(canonicalPriorMatrices)) {
        msg <- glue(
            "Please provide dataDrivenPriorMatrices or set ",
            "canonicalPriorMatrices = TRUE."
        )
        abort(msg)
    }
    invisible(NULL)
}

# glasso needs a complete Y; downgrade to enet with a warning when Y has NAs.
# @noRd
.mrmashResolveBInit <- function(Y, bInitMethod) {
    if (any(is.na(Y)) && bInitMethod == "glasso") {
        msg <- glue(
            "bInitMethod = 'glasso' can only be used without missing ",
            "values in Y. Setting it to 'enet' instead"
        )
        warn(msg)
        return("enet")
    }
    bInitMethod
}

# Univariate summary statistics (Bhat/Shat) for the prior + init.
# @noRd
.mrmashComputeSumstats <- function(X, Y, standardize, numThreads) {
    mr.mashr::compute_univariate_sumstats(
        X,
        Y,
        standardize = standardize,
        standardize.response = FALSE,
        mc.cores = numThreads
    )
}

# Initial coefficient matrix via graphical-lasso or univariate glmnet.
# @noRd
.mrmashInitCoefficients <- function(
    X,
    Y,
    bInitMethod,
    standardize,
    numThreads
) {
    if (bInitMethod == "glasso") {
        return(computeCoefficientsGlasso(
            X,
            Y,
            standardize = standardize,
            numThreads = numThreads,
            Xnew = NULL
        ))
    }
    computeCoefficientsUnivGlmnet(
        X,
        Y,
        alpha = 0.5,
        standardize = standardize,
        Xnew = NULL
    )
}

# Robust residual-covariance init. Returns list(V, updateV); a rank-deficient V
# is ridge-regularized and its update disabled.
# @noRd
.mrmashInitV <- function(X, Y, V, updateV, updateVMethod) {
    if (!is.null(V)) {
        return(list(V = V, updateV = updateV))
    }
    V <- .mrmashComputeVInit(X, Y, any(is.na(Y)))
    if (updateVMethod == "diagonal") {
        return(list(V = diag(diag(V)), updateV = updateV))
    }
    if (any(eigen(V)$values < 1e-8)) {
        return(list(V = V + diag(1e-8, nrow(V)), updateV = FALSE))
    }
    list(V = V, updateV = updateV)
}

# Compute V_init via mr.mashr (cov for complete Y, flash when Y has missing).
# @noRd
.mrmashComputeVInit <- function(X, Y, yHasMissing) {
    zeroB <- matrix(0, nrow = ncol(X), ncol = ncol(Y))
    if (!yHasMissing) {
        return(mr.mashr:::compute_V_init(
            X,
            Y,
            zeroB,
            rep(0, ncol(Y)),
            method = "cov"
        ))
    }
    mr.mashr:::compute_V_init(
        X,
        Y,
        zeroB,
        colMeans(Y, na.rm = TRUE),
        method = "flash"
    )
}

# Run mr.mash with the resolved prior / init / V.
# @noRd
# The mr.mash call. Everything pecotmr derives -- the prior, the residual
# covariance, the coefficient initialisation -- is supplied here and cannot be
# overridden; the rest of mr.mash's arguments come from `methodArgs` under
# mr.mash's own names, with pecotmr's defaults underneath.
# @noRd
.mrmashFit <- function(
    S0,
    bInit,
    vInit,
    X,
    Y,
    standardize,
    numThreads,
    updateVMethod,
    methodArgs
) {
    derived <- list(
        X = X,
        Y = Y,
        V = vInit$V,
        S0 = S0,
        w0 = computeW0(bInit, length(S0)),
        update_V = vInit$updateV,
        update_V_method = updateVMethod,
        mu1_init = bInit,
        standardize = standardize,
        nthreads = numThreads
    )
    user <- as.list(methodArgs)
    clash <- intersect(names(user), names(derived))
    if (length(clash) > 0L) {
        abort(glue(
            "mrmashWrapper: {str_flatten(clash, ', ')} ",
            "{if (length(clash) == 1L) 'is' else 'are'} derived by the ",
            "wrapper (from the data, the prior and the initialisation), so ",
            "{if (length(clash) == 1L) 'it' else 'they'} cannot be set in ",
            "`methodArgs`."
        ))
    }
    defaults <- list(
        update_w0 = TRUE,
        tol = 0.01,
        max_iter = 5000,
        convergence_criterion = "ELBO",
        compute_ELBO = TRUE,
        verbose = FALSE,
        w0_threshold = 1e-8
    )
    exec(mr.mashr::mr.mash, !!!c(derived, list_modify(defaults, !!!user)))
}

# @noRd
.rrDropIntercept <- function(coefs) {
    as.vector(coefs)[-1]
}

#' Compute initial mr.mash coefficients via group-lasso
#'
#' Fit a group-lasso (one group per response) to obtain initial estimates of the
#' \code{mr.mash} coefficient matrix, optionally predicting on a new design.
#'
#' @param X Numeric design matrix (samples x variants).
#' @param Y Numeric response matrix (samples x conditions).
#' @param standardize Logical; standardize the columns of \code{X} before
#'   fitting.
#' @param numThreads Integer; number of threads for the group-lasso fit.
#' @param Xnew Optional new design matrix for prediction; \code{NULL} returns
#'   in-sample coefficients.
#' @return A numeric matrix of initial coefficient estimates (variants x
#'   conditions).
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' Y <- matrix(rnorm(nrow(X) * 3), nrow(X), 3)
#' computeCoefficientsGlasso(X = X, Y = Y, standardize = TRUE,
#'   numThreads = 1L, Xnew = NULL)
#' @importFrom checkmate assertFlag assertInt
#' @export
computeCoefficientsGlasso <- function(
    X,
    Y,
    standardize,
    numThreads,
    Xnew = NULL
) {
    assertFlag(standardize)
    assertInt(numThreads)
    n <- nrow(X)
    p <- ncol(X)
    r <- ncol(Y)
    conditionNames <- colnames(Y)

    # Fit group-lasso
    cvfitGlmnet <- glmnet::cv.glmnet(
        x = X,
        y = Y,
        family = "mgaussian",
        alpha = 1,
        standardize = standardize,
        parallel = FALSE
    )
    coeffGlmnet <- coef(cvfitGlmnet, s = "lambda.min")

    # Build matrix of initial estimates for mr.mash: one column per outcome,
    # each the glmnet coefficients with the intercept dropped.
    B <- matrix(
        unname(list_c(map(coeffGlmnet, .rrDropIntercept))),
        nrow = p,
        ncol = r
    )

    # Make predictions if requested.
    if (!is.null(Xnew)) {
        YhatGlmnet <- `colnames<-`(
            drop(predict(cvfitGlmnet, newx = Xnew, s = "lambda.min")),
            conditionNames
        )
        res <- list(Bhat = B, Ytrain = Y, Yhat_new = YhatGlmnet)
    } else {
        res <- list(Bhat = B, Ytrain = Y)
    }
    return(res)
}


# Fit a cv.glmnet for outcome column `i` on its non-missing rows and return the
# lambda.min coefficients (plus predictions on `Xnew` when supplied).
# @noRd
.linreg <- function(i, X, Y, alpha, standardize, Xnew) {
    samplesKept <- which(!is.na(Y[, i]))
    Ynomiss <- Y[samplesKept, i, drop = FALSE]
    Xnomiss <- X[samplesKept, , drop = FALSE]

    cvfit <- glmnet::cv.glmnet(
        x = Xnomiss,
        y = Ynomiss,
        family = "gaussian",
        alpha = alpha,
        standardize = standardize,
        parallel = FALSE
    )
    coeffic <- as.vector(coef(cvfit, s = "lambda.min"))
    lambdaSeq <- cvfit$lambda

    # Make predictions if requested
    if (!is.null(Xnew)) {
        yhatGlmnet <- drop(predict(cvfit, newx = Xnew, s = "lambda.min"))
        res <- list(
            bhat = coeffic,
            lambda_seq = lambdaSeq,
            yhat_new = yhatGlmnet
        )
    } else {
        res <- list(bhat = coeffic, lambda_seq = lambdaSeq)
    }

    return(res)
}

### Function to compute coefficients for univariate glmnet
computeCoefficientsUnivGlmnet <- function(
    X,
    Y,
    alpha,
    standardize,
    Xnew = NULL
) {
    r <- ncol(Y)

    out <- map(seq_len(r), .linreg, X, Y, alpha, standardize, Xnew)

    bhatList <- map(out, "bhat")
    Bhat <- exec(cbind, !!!bhatList)

    if (!is.null(Xnew)) {
        yhatList <- map(out, "yhat_new")
        YhatNew <- `colnames<-`(exec(cbind, !!!yhatList), colnames(Y))
        results <- list(
            Bhat = Bhat[-1, ],
            intercept = Bhat[1, ],
            Yhat_new = YhatNew
        )
    } else {
        results <- list(Bhat = Bhat[-1, ], intercept = Bhat[1, ])
    }
    return(results)
}


### Compute prior weights from coefficients estimates
computeW0 <- function(Bhat, ncomps) {
    propNonzero <- sum(rowSums(abs(Bhat)) > 0) / nrow(Bhat)

    fromData <- if (ncomps > 1) {
        c((1 - propNonzero), rep(propNonzero / (ncomps - 1), (ncomps - 1)))
    } else {
        1
    }
    # Fewer than two non-zero components leaves nothing to mix: fall back to
    # a flat prior over all of them.
    if (sum(fromData != 0) < 2) {
        return(rep(1 / ncomps, ncomps))
    }
    fromData
}


#' Re-normalize mrmash weight w0 to have total weight sum to 1
#' @param w0 is the weight of mr.mash prior matrices that was generated from
#'   mr.mash() function.
#' @return A named numeric vector of prior-matrix weights with the \code{null}
#'   component removed and the remaining weights renormalized to sum to 1.
#' @keywords internal
rescaleCovW0 <- function(w0) {
    # remove null component
    nonNull <- w0[names(w0) != "null"]

    # split by prior group
    groups <- str_remove(names(nonNull), "_[^_]+$")
    groupList <- split(nonNull, groups)

    # get per group sum -- one scalar per group
    groupSums <- map_dbl(groupList, sum)
    sumWeights <- sum(groupSums)
    weightsList <- if (sumWeights > 0) {
        groupSums / sumWeights
    } else {
        # Use equal weights if all non null weights are zeros
        set_names(
            rep(1 / length(groupSums), length(groupSums)),
            names(groupSums)
        )
    }
    # One w0 slot per group, filled from the supplied weights.
    groupKeys <- unique(groups)
    replace(
        set_names(rep(NA, length(groupKeys)), groupKeys),
        names(weightsList),
        weightsList
    )
}


# The prior scaling grid, from mr.mash's own exported grid builder rather
# than a copy of it. The three helpers this replaces -- gridMin, gridMax,
# autoselectMixsd -- were ports, and the first had DRIFTED: upstream
# grid_min() is `min(Shat)/10`, the port was `min(Shat)`. Every grid built
# here therefore started ten times too high and, since the point count is
# ceiling(log2(gmax/gmin)/log2(mult)), carried fewer points than mr.mash
# would have used.
#
# `mult = sqrt(2)` (finer than upstream's default 2) and the square are
# pecotmr's own choices: this grid scales variances, not standard
# deviations. autoselect.mixsd() drops Shat that is zero or non-finite,
# which covers NA, so the port's extra is.na() term was redundant.
# @noRd
computeGrid <- function(bhat, sbhat) {
    mr.mashr::autoselect.mixsd(
        list(Bhat = bhat, Shat = sbhat),
        mult = sqrt(2)
    )^2
}


#' Compute diagonal covariance matrix
#'
#' Returns a diagonal covariance matrix from the column-wise variances of Y.
#'
#' @param Y Numeric matrix (samples x conditions).
#' @return A diagonal covariance matrix of dimension ncol(Y) x ncol(Y).
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' Y <- matrix(rnorm(nrow(X) * 3), nrow(X), 3)
#' computeCovDiag(Y = Y)
#' @export
computeCovDiag <- function(Y) {
    diag(apply(Y, 2, var, na.rm = TRUE))
}

# Split a flat mr.mash argument bag into the prior-construction names and
# everything else. mrmashWrapper takes the former as a `prior` record, so a
# caller who spells them flat -- as the pipeline's per-token kwargs do -- is
# routed rather than rejected.
# @noRd
.mrmashSplitPriorArgs <- function(args) {
    priorNames <- names(formals(mrmashPriorConfig))
    isPrior <- is_in(names(args), priorNames)
    list(prior = args[isPrior], rest = args[!isPrior])
}

#' @title Prior-Covariance Arguments For mr.mash
#' @description Options for \code{\link{buildMrmashPriorMatrices}}, the
#'   shared prior-covariance builder behind \code{\link{mrmashWrapper}} and
#'   \code{\link{mrmashRssWeights}}. Every field is a formal of that
#'   function, so an unknown name is rejected by R as an unused argument.
#'
#'   \code{dataDrivenPriorMatrices} is \emph{not} set here: it is an input
#'   computed upstream (typically by \code{\link{mashPriorCovariances}}),
#'   not a tuning choice, so it stays a parameter of the caller.
#' @param canonicalPriorMatrices Logical. Include the canonical mixture from
#'   \code{mr.mashr::compute_canonical_covs()}. Default \code{FALSE}.
#' @param priorGrid Optional pre-computed scaling grid. \code{NULL} (default)
#'   derives one from the summary statistics.
#' @param hetgrid Heterogeneity grid for the canonical mixture. \code{NULL}
#'   (default) leaves the builder's own grid in place.
#' @param singletons Logical or \code{NULL}. Include single-condition prior
#'   components; \code{NULL} (default) leaves the builder's own choice.
#' @param zeromat Logical or \code{NULL}. Include the all-zero (null)
#'   component in the expanded mixture; \code{NULL} (default) leaves the
#'   builder's own choice.
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' mrmashPriorConfig(canonicalPriorMatrices = TRUE)
#' @export
mrmashPriorConfig <- function(
    canonicalPriorMatrices = FALSE,
    priorGrid = NULL,
    hetgrid = NULL,
    singletons = NULL,
    zeromat = NULL
) {
    .newMethodConfig(
        NULL,
        defaults = list(
            canonicalPriorMatrices = canonicalPriorMatrices,
            priorGrid = priorGrid,
            hetgrid = hetgrid,
            singletons = singletons,
            zeromat = zeromat
        ),
        extra = list(),
        label = "mrmashPriorConfig",
        engine = "mrmashPrior"
    )
}

#' Build mr.mash prior covariance matrices
#'
#' Shared helper used by both \code{\link{mrmashWrapper}} (individual-level) and
#' \code{\link{mrmashRssWeights}} (summary statistics). Constructs the \code{S0}
#' list of prior covariance matrices via the canonical mixture
#' (\code{mr.mashr::compute_canonical_covs}) and optional data-driven matrices,
#' expanded over a scaling grid via \code{mr.mashr::expand_covs}. The prior grid
#' is derived from \code{Bhat} and \code{Shat} via \code{computeGrid} when not
#' supplied.
#'
#' @param Bhat Numeric matrix of effect-size estimates (variants x conditions).
#' @param Shat Numeric matrix of standard errors (variants x conditions).
#' @param K Number of conditions. When NULL, inferred from \code{ncol(Bhat)}.
#' @param dataDrivenPriorMatrices Optional list with element \code{U} (list of
#'   raw covariance matrices) computed e.g. by
#'   \code{\link{computeCovDiag}}, or by \code{mashr::cov_flash}, which
#'   returns such a list directly.
#' @param canonicalPriorMatrices Logical. When TRUE (default for RSS), include
#'   the standard canonical mixture from
#'   \code{mr.mashr::compute_canonical_covs()}. When FALSE,
#'   \code{dataDrivenPriorMatrices} must be supplied.
#' @param priorGrid Optional pre-computed scaling grid (numeric vector). When
#'   NULL, derived from \code{Bhat}, \code{Shat} via \code{computeGrid()}.
#' @param hetgrid Heterogeneity grid passed to
#'   \code{mr.mashr::compute_canonical_covs()}. Default \code{c(0, 0.25, 0.5,
#'   0.75, 1)}, matching the individual-level wrapper.
#' @param singletons Whether to include single-condition prior components.
#'   Default TRUE.
#' @param zeromat Logical. Passed to \code{mr.mashr::expand_covs()}: whether
#'   the expanded mixture includes the all-zero (null) component. Default
#'   \code{TRUE}, which is mr.mash's own convention --- a mixture with no
#'   null component cannot shrink an effect to zero.
#' @return A list with components \code{S0} (the expanded list of prior
#'   covariance matrices) and \code{prior_grid} (the scaling grid that was
#'   used).
#' @examples
#' Bhat <- matrix(rnorm(15), 5, 3)
#' Shat <- matrix(abs(rnorm(15)) + 0.1, 5, 3)
#' buildMrmashPriorMatrices(Bhat = Bhat, Shat = Shat)
#' @export
buildMrmashPriorMatrices <- function(
    Bhat,
    Shat,
    K = NULL,
    dataDrivenPriorMatrices = NULL,
    canonicalPriorMatrices = TRUE,
    priorGrid = NULL,
    hetgrid = c(0, 0.25, 0.5, 0.75, 1),
    singletons = TRUE,
    zeromat = TRUE
) {
    if (!requireNamespace("mr.mashr", quietly = TRUE)) {
        abort("Package 'mr.mashr' is required.")
    }
    assertFlag(zeromat)
    if (is.null(dataDrivenPriorMatrices) && !isTRUE(canonicalPriorMatrices)) {
        msg <- glue(
            "Supply dataDrivenPriorMatrices or set ",
            "canonicalPriorMatrices = TRUE."
        )
        abort(msg)
    }
    if (is.null(K)) {
        K <- ncol(Bhat)
    }
    if (is.null(priorGrid)) {
        priorGrid <- computeGrid(bhat = Bhat, sbhat = Shat)
    }

    if (isTRUE(canonicalPriorMatrices)) {
        canonical <- mr.mashr::compute_canonical_covs(
            K,
            singletons = singletons,
            hetgrid = hetgrid
        )
        S0_raw <- if (!is.null(dataDrivenPriorMatrices)) {
            c(canonical, dataDrivenPriorMatrices$U)
        } else {
            canonical
        }
    } else {
        S0_raw <- dataDrivenPriorMatrices$U
    }

    S0 <- mr.mashr::expand_covs(S0_raw, priorGrid, zeromat = zeromat)
    list(S0 = S0, priorGrid = priorGrid)
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# `.rssSolvePath` solver for lassosum: RSS lasso over the ordered lambda path.
# @noRd
.rssLassosumSolve <- function(z, lam, R, thr, maxiter) {
    lassosumRssRcpp(
        zR = z,
        LD = list(blk1 = R),
        lambdaR = lam,
        thr = thr,
        maxiter = as.integer(maxiter)
    )
}

# `.rssSolvePath` solver for the penalized (MCP/SCAD/L0...) RSS path.
# @noRd
.rssPenalizedSolve <- function(
    z,
    lam,
    R,
    penalty,
    gamma,
    alpha,
    lambda0,
    lambda2,
    thr,
    maxiter,
    maxSwaps
) {
    penalizedRssRcpp(
        zR = z,
        LD = list(blk1 = R),
        lambdaR = lam,
        penaltyStr = penalty,
        gamma = gamma,
        alpha = alpha,
        lambda0 = lambda0,
        lambda2 = lambda2,
        thr = thr,
        maxiter = as.integer(maxiter),
        maxSwaps = as.integer(maxSwaps)
    )
}

# One L0Learn fit at lambda0 = `l0Val`: penalizedRss over the lambda path.
# @noRd
.rssL0LambdaFit <- function(l0Val, solverInput, LDs, n, sVal, config) {
    penArgs <- c(
        list(
            bhat = solverInput,
            R = LDs,
            n = n,
            penalty = config$penalty,
            lambda = config$lambda,
            lambda0 = l0Val,
            lambda2 = config$lambda2,
            maxSwaps = config$maxSwaps
        ),
        as.list(config$dotArgs)
    )
    model <- exec(penalizedRss, !!!penArgs)
    list(
        beta = model$beta,
        meta = tibble(
            s = rep(sVal, length(model$lambda)),
            lambda0 = rep(l0Val, length(model$lambda)),
            lambda = model$lambda,
            fbeta = model$fbeta
        )
    )
}

# The minimum of one CV-mean path vector (coerced numeric).
# @noRd
.rssMinNumeric <- function(v) {
    min(as.numeric(v))
}
