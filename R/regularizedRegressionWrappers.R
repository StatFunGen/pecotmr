#' @include MethodParam.R
NULL

#' Extract weights from mr.ash.rss (susieR)
#' @param stat A list of summary statistics: \code{n} (per-variant sample
#'   sizes), \code{variantNames}, and the marginal effects as \code{z}, or as
#'   \code{b} with \code{seb}.
#' @param ldReference The LD reference the statistics' variants are matched
#'   to by \code{stat$variantNames}: a genotype panel (see
#'   \code{\link{readGenotypes}}), a \code{\link{GenotypeMixture}} or an
#'   \code{\link{LdMatrix}}. \code{mr.ash.rss} takes an LD matrix, which is
#'   formed from the reference over these variants.
#' @param varY Variance of the phenotype, or \code{NULL} (default) to fit on
#'   the standardized scale. Used only with \code{b} and \code{seb}.
#' @param sigma2E Initial residual variance, or \code{NULL} (default) for
#'   \code{mr.ash.rss}'s own.
#' @param s0 Prior variances of the mixture components, relative to the
#'   residual variance, or \code{NULL} (default) for the grid
#'   \code{susieR::mr.ash} uses by default, built for these statistics.
#' @param w0 Prior mixture weights (summing to 1), or \code{NULL} (default)
#'   for equal weights.
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights: \code{"none"} (default) attaches
#'   nothing, \code{"slim"} and \code{"full"} attach the fitted
#'   \code{mr.ash.rss} model. This engine keeps no trimmable intermediates,
#'   so the two retaining levels behave alike here.
#' @param methodArgs Options forwarded to \code{susieR::mr.ash.rss}, built
#'   with \code{\link{MrashRssOptions}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of the posterior mean of the coefficients. Unless
#'   \code{fitRetention = "none"}, the fitted model is attached as attribute
#'   \code{"fit"}.
#' @importFrom susieR mr.ash.rss
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ss <- lapply(seq_len(ncol(X)), function(j) {
#'   coef(summary(lm(y ~ X[, j])))[2, 1:2]
#' })
#' stat <- list(b = vapply(ss, `[`, numeric(1), 1L),
#'   seb = vapply(ss, `[`, numeric(1), 2L), n = rep(nrow(X), ncol(X)),
#'   variantNames = colnames(X))
#' mrashRssWeights(stat, LdMatrix(cor(X)))
#' @importFrom checkmate assertList
#' @export
mrashRssWeights <- function(
    stat,
    ldReference,
    varY = NULL,
    sigma2E = NULL,
    s0 = NULL,
    w0 = NULL,
    fitRetention = c("none", "slim", "full"),
    methodArgs = MrashRssOptions()
) {
    fitRetention <- arg_match(fitRetention)
    .assertMethodOptions(methodArgs, "MrashRssOptions", "methodArgs")
    assertList(stat)
    R <- .rssLdMatrix(stat, ldReference, "mrashRssWeights")
    n <- median(stat$n)
    shat <- as.numeric(stat$seb %||% numeric(0))
    z <- .mrashRssZ(stat)
    s0 <- s0 %||% .mrashRssDefaultGrid(n, z, shat, varY)
    callArgs <- list_modify(
        list(
            bhat = as.numeric(stat$b %||% numeric(0)),
            shat = shat,
            z = z,
            R = R,
            var_y = varY,
            n = n,
            sigma2_e = sigma2E,
            s0 = s0,
            # A copy: mr.ash.rss() writes its fitted weights over the vector
            # it is given, which would otherwise be the caller's.
            w0 = (w0 %||% rep(1 / length(s0), length(s0))) + 0
        ),
        !!!as.list(methodArgs)
    )
    model <- exec(mr.ash.rss, !!!callArgs)
    if (identical(fitRetention, "none")) {
        return(model$mu1)
    }
    `attr<-`(
        model$mu1,
        "fit",
        .withFitVariantIds(model, stat[["variantNames"]])
    )
}

# The z-scores mr.ash.rss() fits: the statistics' own, or b / seb.
# @noRd
.mrashRssZ <- function(stat) {
    hasBetaSe <- !is.null(stat$b) && !is.null(stat$seb)
    z <- stat$z %||% (if (hasBetaSe) stat$b / stat$seb)
    if (is.null(z)) {
        abort("mrashRssWeights: `stat` must carry `z`, or `b` with `seb`.")
    }
    as.numeric(z)
}

# susieR::mr.ash()'s own default prior grid, (2^((k - 1) / 25) - 1)^2 scaled
# by n / median(w) with w the columns' sums of squares, for the sufficient
# statistics mr.ash.rss() builds from these: w = n - 1 on the standardized
# scale it fits without both a finite varY and seb, and otherwise
# varY * adj / seb^2, with adj = (n - 1) / (z^2 + n - 2) the shrinkage it
# applies.
# @noRd
.mrashRssDefaultGrid <- function(n, z, shat, varY) {
    unstandardized <- length(shat) > 0L && !is.null(varY) && is.finite(varY)
    w <- if (unstandardized) {
        varY * (n - 1) / (z^2 + n - 2) / shat^2
    } else {
        rep(n - 1, length(z))
    }
    (2^((0:24) / 25) - 1)^2 / median(w) * n
}

#' PRS-CS: a polygenic prediction method that infers posterior SNP effect sizes
#' under continuous shrinkage (CS) priors
#'
#' This function is a wrapper for the PRS-CS method implemented in C++. It takes
#' marginal effect size estimates from regression and an external LD reference
#' panel and infers posterior SNP effect sizes using Bayesian regression with
#' continuous shrinkage priors.
#'
#' Internal: \code{\link{prsCsWeights}} is the public interface. It forwards
#' the arguments below other than \code{bhat}, \code{R} and \code{n} through
#' \code{\link{PrsCsOptions}}, and with \code{fitRetention} returns this
#' function's list as the \code{"fit"} attribute of its weights.
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
#' @param shrinkCorr The LD shrinkage \eqn{c}: the sampler reads
#'   \eqn{cR + (1 - c)I}. Default 1, none.
#' @param raoBlackwell Report the blocks' conditional posterior means,
#'   averaged over the kept iterations, rather than the draws' average.
#'   Default \code{FALSE}.
#' @param nChains The number of chains, from seeds \code{seed},
#'   \code{seed + 1}, ... (each drawn when \code{seed} is \code{NULL}).
#'   A single chain warns when its posterior mean effects' squares sum past
#'   twice the squared marginal effects, a sign that it diverged. With more
#'   than one, a chain stops once a kept draw's squared effects do, and is
#'   dropped, and the rest are filtered and averaged as LDpred2's vignette
#'   does (see \code{\link{ldpred2AutoWeights}}). Default 1.
#'
#' @return A list containing the posterior estimates: - betaEst: Posterior
#'   estimates of SNP effect sizes. - psiEst: Posterior estimates of psi
#'   (shrinkage parameters). - sigmaEst: Posterior estimate of the residual
#'   variance. - phiEst: Posterior estimate of the global shrinkage parameter.
#'   With several chains, these are the kept chains' means, and
#'   \code{chains} records each chain's seed and whether it diverged and was
#'   kept.
#' @keywords internal
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
    seed = NULL,
    shrinkCorr = 1,
    raoBlackwell = FALSE,
    nChains = 1
) {
    .rssValidateInputs(bhat, R, n)
    .prsCsRun(
        bhat,
        list(R = R),
        n,
        a,
        b,
        phi,
        maf,
        nIter,
        nBurnin,
        thin,
        verbose,
        seed,
        shrinkCorr,
        raoBlackwell,
        nChains
    )
}

# PRS-CS on an LD source (.ldSourceRows()): the sampler that factors R for a
# matrix, or for a factor the one that never forms R (FactorLdBlock in
# rss_mcmc.h), which draws differently from the same posterior.
# @noRd
.prsCsRun <- function(
    bhat,
    ld,
    n,
    a,
    b,
    phi,
    maf,
    nIter,
    nBurnin,
    thin,
    verbose,
    seed,
    shrinkCorr = 1,
    raoBlackwell = FALSE,
    nChains = 1
) {
    if (!is.null(maf) && length(bhat) != length(maf)) {
        abort("The length of 'bhat' must be the same as 'maf'.")
    }
    .mcmcAssertShared(shrinkCorr, raoBlackwell, nChains)
    # cpp11 requires exact integer types for int parameters
    settings <- list(
        a = a,
        b = b,
        phi = phi,
        bhat = bhat,
        maf = maf,
        n = as.integer(n),
        nIter = as.integer(nIter),
        nBurnin = as.integer(nBurnin),
        thin = as.integer(thin),
        verbose = verbose,
        shrinkCorr = shrinkCorr,
        raoBlackwell = raoBlackwell
    )
    .mcmcChains(
        ".prsCsChain",
        list(settings = settings, ld = ld),
        .mcmcSeeds(seed, nChains),
        c("betaEst", "psiEst", "sigmaEst", "phiEst"),
        "prsCs",
        bhat
    )
}

# One PRS-CS chain, on the LD matrix or a factor of it.
# @noRd
.prsCsChain <- function(settings, ld, seed, track) {
    if (is.null(ld$F)) {
        return(exec(
            prsCsCpp,
            !!!settings,
            ldBlk = list(blk1 = ld$R),
            seed = seed,
            track = track
        ))
    }
    exec(
        prsCsFactorCpp,
        !!!settings,
        factorR = ld$F,
        seed = seed,
        track = track
    )
}

#' Extract weights from prsCs function
#' @param stat A list of summary statistics: \code{n} (per-variant sample
#'   sizes) and the marginal effects as \code{cor} (standardized), \code{z},
#'   \code{b} with \code{seb}, or \code{b} alone, taken as already
#'   standardized. The engine receives the standardized effects z / sqrt(n),
#'   with n the median of \code{n}.
#' @param ldReference The LD reference the statistics' variants are matched
#'   to by \code{stat$variantNames}: a genotype panel (see
#'   \code{\link{readGenotypes}}), a \code{\link{GenotypeMixture}} or an
#'   \code{\link{LdMatrix}}. A genotype reference's LD is used through its
#'   standardized genotypes, never formed as a matrix.
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights: \code{"none"} (default) attaches
#'   nothing, \code{"slim"} and \code{"full"} attach the sampler's posterior
#'   estimates: \code{betaEst}, the local shrinkage parameters
#'   \code{psiEst}, the residual variance \code{sigmaEst} and the global
#'   shrinkage \code{phiEst}. This engine keeps no trimmable intermediates,
#'   so the two retaining levels behave alike here.
#' @param methodArgs Options forwarded to the sampler, \code{\link{prsCs}},
#'   built with \code{\link{PrsCsOptions}}. A bare list is refused: it cannot
#'   be checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of the posterior SNP coefficients. Unless
#'   \code{fitRetention = "none"}, the sampler's estimates are attached as
#'   attribute \code{"fit"}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ss <- lapply(
#'   seq_len(ncol(X)), function(j) coef(summary(lm(y ~ X[, j])))[2, 1:2])
#' stat <- list(
#'   b = vapply(ss, `[`, numeric(1), 1L),
#'   seb = vapply(ss, `[`, numeric(1), 2L),
#'   n = rep(nrow(X), ncol(X)),
#'   variantNames = colnames(X)
#' )
#' w <- prsCsWeights(
#'   stat,
#'   LdMatrix(cor(X)),
#'   fitRetention = "slim",
#'   methodArgs = PrsCsOptions(maf = rep(0.3, ncol(X)))
#' )
#' attr(w, "fit")$phiEst
#' @export
prsCsWeights <- function(
    stat,
    ldReference,
    fitRetention = c("none", "slim", "full"),
    methodArgs = PrsCsOptions()
) {
    fitRetention <- arg_match(fitRetention)
    .assertMethodOptions(methodArgs, "PrsCsOptions", "methodArgs")
    ld <- .rssLdSource(stat, ldReference, "prsCsWeights")
    n <- median(stat$n)
    bhat <- .rssMarginalCor(stat, n = n, p = .ldSourceSize(ld))
    .rssValidateSource(bhat, ld, n)
    settings <- list_modify(
        .formalDefaults(
            prsCs,
            c(
                "a",
                "b",
                "phi",
                "maf",
                "nIter",
                "nBurnin",
                "thin",
                "verbose",
                "seed",
                "shrinkCorr",
                "raoBlackwell",
                "nChains"
            )
        ),
        !!!as.list(methodArgs)
    )
    model <- exec(.prsCsRun, bhat = bhat, ld = ld, n = n, !!!settings)
    if (identical(fitRetention, "none")) {
        return(model$betaEst)
    }
    `attr<-`(
        model$betaEst,
        "fit",
        .withFitVariantIds(model, stat[["variantNames"]])
    )
}

#' SDPR (Summary-Statistics-Based Dirichelt Process Regression for Polygenic
#' Risk Prediction)
#'
#' This function is a wrapper for the SDPR C++ implementation, which performs
#' Markov Chain Monte Carlo (MCMC) for estimating effect sizes and heritability
#' based on summary statistics and reference LD matrices.
#'
#' Internal: \code{\link{sdprWeights}} is the public interface. It forwards
#' the arguments below other than \code{bhat}, \code{R} and \code{n} through
#' \code{\link{SdprOptions}}, and with \code{fitRetention} returns this
#' function's list as the \code{"fit"} attribute of its weights.
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
#' @param optLlk Which likelihood to evaluate. 1 for equation 6 (slightly shrink
#'   the correlation of SNPs) and 2 for equation 5 (SNPs genotyped on different
#'   arrays in a separate cohort). Default is 1.
#' @param verbose Whether to print verbose output. Default is true.
#'
#' @param seed Integer or \code{NULL}. Random seed for the Gibbs sampler;
#'   \code{NULL} leaves the RNG state unchanged.
#' @param raoBlackwell Report \eqn{\eta} times the effects' conditional
#'   posterior means, averaged over the kept iterations, rather than the
#'   draws' average. Default \code{FALSE}.
#' @param nChains The number of chains, from seeds \code{seed},
#'   \code{seed + 1}, ...; they are tested for divergence and, with more
#'   than one, combined as \code{\link{prsCs}}'s are. Default 1.
#' @return A list with \code{betaEst} (the estimated effect sizes) and
#'   \code{h2} (the estimated heritability), and with several chains
#'   \code{chains} as \code{\link{prsCs}} records it.
#' @note This function wraps a rewritten and adapted version of the SDPR C++
#'   implementation. SDPR is described in Zhou G, Zhao H (2021), "A fast and
#'   robust Bayesian nonparametric method for prediction of complex traits
#'   using summary statistics", PLoS Genetics 17(7): e1009697.
#'   \doi{10.1371/journal.pgen.1009697}
#'
#'   It departs from the original where a window has few causal variants.
#'   A lone causal variant is drawn jointly, as several are; the original's
#'   closed-form draw let the global scale of the effects run away. That
#'   scale has an N(1, 1) prior where the original's is flat, under which a
#'   window that drew no causal variant tended to stay null.
#'
#' @keywords internal
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
    optLlk = 1,
    verbose = TRUE,
    seed = NULL,
    raoBlackwell = FALSE,
    nChains = 1
) {
    .rssValidateInputs(bhat, R, n)
    .sdprRun(
        bhat,
        list(R = R),
        n,
        perVariantSampleSize,
        array,
        a,
        c,
        M,
        a0k,
        b0k,
        iter,
        burn,
        thin,
        optLlk,
        verbose,
        seed,
        raoBlackwell,
        nChains
    )
}

# SDPR on an LD source (.ldSourceRows()): the sampler on the matrix, or the
# same sampler on the low-rank form of SDPR's matrices from a factor of it
# (sdpr_factor_block() in rss_mcmc.h), which draws as the matrix route does
# without forming R.
# @noRd
.sdprRun <- function(
    bhat,
    ld,
    n,
    perVariantSampleSize,
    array,
    a,
    c,
    M,
    a0k,
    b0k,
    iter,
    burn,
    thin,
    optLlk,
    verbose,
    seed,
    raoBlackwell = FALSE,
    nChains = 1
) {
    .sdprValidate(M, perVariantSampleSize, array)
    .mcmcAssertShared(1, raoBlackwell, nChains)
    # cpp11 requires exact integer types for int params + sexp vectors.
    settings <- list(
        bhatR = bhat,
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
        optLlk = as.integer(optLlk),
        verbose = verbose,
        raoBlackwell = raoBlackwell
    )
    .mcmcChains(
        ".sdprChain",
        list(settings = settings, ld = ld),
        .mcmcSeeds(seed, nChains),
        c("betaEst", "h2"),
        "sdpr",
        bhat
    )
}

# One SDPR chain, on the LD matrix (a one-block list) or a factor of it.
# @noRd
.sdprChain <- function(settings, ld, seed, track) {
    if (is.null(ld$F)) {
        return(exec(
            sdprCpp,
            !!!settings,
            LD = list(blk1 = ld$R),
            seed = seed,
            track = track
        ))
    }
    exec(sdprFactorCpp, !!!settings, factorR = ld$F, seed = seed, track = track)
}

# Validate the SDPR settings: M >= 4 (SDPR uses M-2 indexing in sample_V;
# M < 4 overflows), positive per-variant N, and array in {0, 1, 2}. The RSS
# contract (bhat, LD, n) is checked by the callers.
# @noRd
.sdprValidate <- function(M, perVariantSampleSize, array) {
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
#' @param stat A list of summary statistics: \code{n} (per-variant sample
#'   sizes) and the marginal effects as \code{cor} (standardized), \code{z},
#'   \code{b} with \code{seb}, or \code{b} alone, taken as already
#'   standardized. The engine receives the standardized effects z / sqrt(n),
#'   with n the median of \code{n}.
#' @param ldReference The LD reference the statistics' variants are matched
#'   to by \code{stat$variantNames}: a genotype panel (see
#'   \code{\link{readGenotypes}}), a \code{\link{GenotypeMixture}} or an
#'   \code{\link{LdMatrix}}. A genotype reference's LD is used through its
#'   standardized genotypes and not formed as a matrix, except under
#'   \code{optLlk = 2} when the per-variant sample sizes and arrays are
#'   too varied for its low-rank form to pay.
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights: \code{"none"} (default) attaches
#'   nothing, \code{"slim"} and \code{"full"} attach the sampler's
#'   estimates, \code{betaEst} and the heritability \code{h2}. This engine
#'   keeps no trimmable intermediates, so the two retaining levels behave
#'   alike here.
#' @param methodArgs Options forwarded to the sampler, \code{\link{sdpr}},
#'   built with \code{\link{SdprOptions}}. A bare list is refused: it cannot
#'   be checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of the posterior SNP coefficients. Unless
#'   \code{fitRetention = "none"}, the sampler's estimates, the heritability
#'   among them, are attached as attribute \code{"fit"}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ss <- lapply(
#'   seq_len(ncol(X)), function(j) coef(summary(lm(y ~ X[, j])))[2, 1:2])
#' stat <- list(
#'   b = vapply(ss, `[`, numeric(1), 1L),
#'   seb = vapply(ss, `[`, numeric(1), 2L),
#'   n = rep(nrow(X), ncol(X)),
#'   variantNames = colnames(X)
#' )
#' w <- sdprWeights(stat, LdMatrix(cor(X)), fitRetention = "slim")
#' attr(w, "fit")$h2
#' @export
sdprWeights <- function(
    stat,
    ldReference,
    fitRetention = c("none", "slim", "full"),
    methodArgs = SdprOptions()
) {
    fitRetention <- arg_match(fitRetention)
    .assertMethodOptions(methodArgs, "SdprOptions", "methodArgs")
    ld <- .rssLdSource(stat, ldReference, "sdprWeights")
    n <- median(stat$n)
    bhat <- .rssMarginalCor(stat, n = n, p = .ldSourceSize(ld))
    .rssValidateSource(bhat, ld, n)
    settings <- list_modify(
        .formalDefaults(
            sdpr,
            c(
                "perVariantSampleSize",
                "array",
                "a",
                "c",
                "M",
                "a0k",
                "b0k",
                "iter",
                "burn",
                "thin",
                "optLlk",
                "verbose",
                "seed",
                "raoBlackwell",
                "nChains"
            )
        ),
        !!!as.list(methodArgs)
    )
    model <- exec(.sdprRun, bhat = bhat, ld = ld, n = n, !!!settings)
    if (identical(fitRetention, "none")) {
        return(model$betaEst)
    }
    `attr<-`(
        model$betaEst,
        "fit",
        .withFitVariantIds(model, stat[["variantNames"]])
    )
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
#'   with \code{\link{MrmashOptions}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return Matrix of variant weights.
#' @examples
#' data(multiTraitData)
#' X <- multiTraitData$X[, 1:60]
#' Y <- multiTraitData$Y
#' ddpm <- multiTraitData$priorMatrices
#' fit <- mrmashWrapper(X = X, Y = Y, dataDrivenPriorMatrices = ddpm,
#'   prior = MrmashPriorParam(canonicalPriorMatrices = TRUE))
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
    methodArgs = MrmashOptions()
) {
    fitRetention <- arg_match(fitRetention)
    .assertMethodOptions(methodArgs, "MrmashOptions", "methodArgs")
    if (!requireNamespace("mr.mashr", quietly = TRUE)) {
        abort("Package 'mr.mashr' is required.")
    }
    if (is.null(mrmashFit)) {
        inform("mrmashFit is not provided; fitting mr.mash now ...")
        if (is.null(X) || is.null(Y)) {
            abort("Both X and Y must be provided if mrmashFit is NULL.")
        }
        # methodArgs is one flat bag; the prior-construction names belong to
        # mrmashWrapper's `prior` group and the rest to its `methodArgs`, as
        # neither is a top-level argument of mrmashWrapper.
        split <- .mrmashSplitPriorArgs(as.list(methodArgs))
        mrmashFit <- exec(
            mrmashWrapper,
            X,
            Y,
            dataDrivenPriorMatrices = dataDrivenPriorMatrices,
            prior = exec(MrmashPriorParam, !!!split$prior),
            methodArgs = exec(MrmashOptions, !!!split$rest)
        )
    }
    out <- mr.mashr::coef.mr.mash(mrmashFit)[-1, , drop = FALSE]
    # mu1 (= out) is already the returned weights; the payload carries only the
    # mvSuSiE data-driven-prior reconstruction inputs (w0, V, and the raw prior
    # matrices), plus the whole fit at fitRetention = "full".
    .mrmashAttachFit(
        out,
        mrmashFit,
        dataDrivenPriorMatrices,
        fitRetention,
        colnames(X)
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
    fitRetention,
    variantIds = NULL
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
    `attr<-`(weights, "fit", .withFitVariantIds(fitList, variantIds))
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
#' @param stat A list with \code{z} (variants x conditions matrix),
#'   \code{n} (per-context numeric vector or scalar) and \code{variantNames}.
#'   May also include \code{Bhat}, \code{Shat} matrices.
#' @param ldReference The LD reference the statistics' variants are matched
#'   to by \code{stat$variantNames}: a genotype panel (see
#'   \code{\link{readGenotypes}}), a \code{\link{GenotypeMixture}} or an
#'   \code{\link{LdMatrix}}. \code{mr.mash.rss} takes an LD matrix, which is
#'   formed from the reference over these variants; it is not read when
#'   \code{mrmashRssFit} is supplied.
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
#'   with \code{\link{MrmashOptions}}. A bare list is refused: it cannot be
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
#' stat <- list(z = zmat, n = nrow(X), variantNames = colnames(X))
#' mrmashRssWeights(stat, LdMatrix(cor(X)),
#'   dataDrivenPriorMatrices = multiTraitData$priorMatrices,
#'   canonicalPriorMatrices = TRUE)
#' @export
#' @importFrom checkmate assertFlag assertList
mrmashRssWeights <- function(
    stat,
    ldReference,
    mrmashRssFit = NULL,
    dataDrivenPriorMatrices = NULL,
    canonicalPriorMatrices = TRUE,
    S0 = NULL,
    w0 = NULL,
    V = NULL,
    covY = NULL,
    fitRetention = c("none", "slim", "full"),
    methodArgs = MrmashRssOptions()
) {
    fitRetention <- arg_match(fitRetention)
    .assertMethodOptions(methodArgs, "MrmashRssOptions", "methodArgs")
    assertList(stat)
    assertFlag(canonicalPriorMatrices)
    .mrmashRssRequirePackage()
    if (is.null(mrmashRssFit)) {
        mrmashRssFit <- .mrmashRssComputeFit(
            stat,
            .rssLdMatrix(stat, ldReference, "mrmashRssWeights"),
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
        fitRetention,
        stat[["variantNames"]]
    )
}

# Require mr.mashr for the RSS mr.mash path.
# @noRd
.mrmashRssRequirePackage <- function() {
    if (!requireNamespace("mr.mashr", quietly = TRUE)) {
        abort("Package 'mr.mashr' is required for mrmashRssWeights().")
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
    prior <- .mrmashRssDefaults(
        ss,
        dataDrivenPriorMatrices,
        canonicalPriorMatrices,
        S0,
        w0,
        V,
        covY
    )
    # mr.mash.rss expects either Z or (Bhat, Shat) but not both; prefer
    # Bhat/Shat. n must be a scalar (mr.mash.rss contract); use the median.
    callArgs <- c(
        list(
            Bhat = ss$Bhat,
            Shat = ss$Shat,
            R = LD,
            n = as.numeric(stats::median(ss$nVec)),
            covY = prior$covY,
            V = prior$V,
            S0 = prior$S0,
            w0 = prior$w0,
            # One thread: see .mrmashUnavailable().
            nthreads = 1L
        ),
        # `dots` is the caller's MrmashOptions() record; c() would append it as
        # one opaque element instead of splicing its entries.
        as.list(dots)
    )
    exec(mr.mashr::mr.mash.rss, !!!callArgs)
}

# The prior (S0, w0) and the residual and Y covariances (V, covY) mr.mash.rss
# is given: each as supplied, else its default for the statistics `ss`
# (.mrmashRssStats()).
# @noRd
.mrmashRssDefaults <- function(
    ss,
    dataDrivenPriorMatrices,
    canonicalPriorMatrices,
    S0,
    w0,
    V,
    covY
) {
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
    list(S0 = S0, w0 = w0, V = V, covY = covY)
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
#'   the result as an element of \code{\link{TwasWeightsMethodsParam}}, or
#'   directly to the matching weights function.
#'
#'   A method's individual-level and summary-statistic paths often run
#'   DIFFERENT engines, and then each takes its own constructor:
#'   \code{lassoWeights} fits with glmnet (\code{GlmnetOptions}) while
#'   \code{lassosumRssWeights} runs pecotmr's lassosum
#'   (\code{LassosumOptions}), or with \code{backend = "bigsnpr"} bigsnpr's
#'   lassosum2 (\code{BigsnprOptions}, whose \code{delta} grid comes from
#'   the wrapper's \code{s}); \code{enetWeights} / \code{ridgeWeights} fit
#'   with glmnet too, while \code{enetRssWeights} / \code{ridgeRssWeights}
#'   use \code{PenalizedRssOptions}; \code{scadWeights} / \code{mcpWeights} /
#'   \code{l0learnWeights} use ncvreg and L0Learn while their \code{*Rss}
#'   counterparts use \code{PenalizedRssOptions}; \code{dprGibbsWeights} uses
#'   RcppDPR and \code{sdprWeights} uses \code{SdprOptions}. Splitting them is
#'   what lets each be checked exactly, rather than against a union in which
#'   an option meant for the other path would pass and then be dropped.
#'
#'   The Bayesian alphabet methods run on several engines each, chosen with
#'   their \code{backend}, and each engine takes its own constructor:
#'   \code{QggOptions} for qgg on either path, \code{BglrOptions} for BGLR,
#'   \code{HibayesOptions} for hibayes on individual-level data and
#'   \code{HibayesRssOptions} on summary statistics, and
#'   \code{BayesAlphabetRssOptions} for pecotmr's own summary-statistics
#'   sampler. Their wrappers set the iterations, burn-in, thinning and prior
#'   inclusion probability themselves (\code{nIter}, \code{nBurnin},
#'   \code{thin}, \code{pi}), so every constructor refuses its engine's
#'   spelling of them, and hibayes' \code{seed} (drawn from R's random
#'   number generator) and \code{threads} (one).
#'
#'   \code{GlmnetOptions}, \code{NcvregOptions} and \code{DprOptions} cannot
#'   check for an UNKNOWN name: \code{cv.glmnet}, \code{cv.ncvreg} and
#'   \code{fit_model} each take \code{...}, so they accept any name and there
#'   is nothing to check against. Every other constructor here rejects an
#'   unknown option. Separately, and in all of them, a name pecotmr supplies
#'   itself is refused --- \code{GlmnetOptions(alpha=)} and
#'   \code{NcvregOptions(penalty=)} select the penalty family, which is the
#'   calling wrapper's to set. See \code{\link{MethodOptions}}.
#'
#'   \code{GlmnetOptions} rejects two names by hand, which the \code{...}
#'   check above would otherwise let through in silence. \code{nfold} is not
#'   a \code{cv.glmnet} argument at all --- the fold count is \code{nfolds},
#'   and the shorter spelling only ever worked through R's partial matching.
#'   \code{parallel} is a real one, but pecotmr does not expose it: glmnet
#'   parallelises over its CV folds through a \code{foreach} backend that
#'   this package neither declares nor registers, and weight fitting already
#'   runs inside \code{BiocParallel} workers over the same folds, so
#'   enabling it would nest two frameworks and oversubscribe the machine.
#'   The glasso and univariate-glmnet \code{mr.mash} initialisations are
#'   single-threaded for the same reason.
#'
#'   \code{MrmashOptions} and \code{MrmashRssOptions} refuse
#'   \code{nthreads}: mr.mash itself runs on one thread, because on macOS
#'   its worker threads can later hang R inside the OpenMP runtime.
#' @param ... Arguments for the engine, under its own names.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' GlmnetOptions(nfolds = 10)
#' QggOptions(updateB = TRUE)
#' LassosumOptions(thr = 1e-4)
#' @name twasEngineOptions
NULL

#' @rdname twasEngineOptions
#' @export
GlmnetOptions <- function(...) {
    extra <- list(...)
    # cv.glmnet takes `...`, so neither of these would be reported by the
    # name check: the first is not an argument at all, the second is one
    # pecotmr will not expose. Both are called out by hand instead.
    .configRefuseUnavailable(
        extra,
        c(
            nfold = str_c(
                "it is not an argument of glmnet::cv.glmnet(), whose fold ",
                "count is `nfolds`. pecotmr passed this misspelling itself ",
                "until recently, so it reached the engine only through R's ",
                "partial matching"
            ),
            parallel = str_c(
                "glmnet parallelises over its CV folds through a foreach ",
                "backend, which pecotmr neither declares nor registers. ",
                "Weight fitting already runs inside BiocParallel workers ",
                "over the same folds, so a second framework underneath it ",
                "would oversubscribe the machine rather than speed it up"
            )
        ),
        "GlmnetOptions"
    )
    .configRefuseOwned(
        extra,
        c(
            x = "supplied from the data by the pipeline",
            y = "supplied from the data by the pipeline",
            # `alpha` selects the penalty family -- it is what distinguishes
            # lasso from enet from ridge -- so it belongs to the calling
            # wrapper. Refused here rather than in glmnetWeights() so the
            # failure lands where the value was written.
            alpha = str_c(
                "the calling wrapper: enetWeights' own `alpha`, or ",
                "lassoWeights / ridgeWeights, which fix it"
            )
        ),
        "GlmnetOptions"
    )
    .newMethodOptions(
        "glmnet::cv.glmnet",
        defaults = list(),
        extra = list(...),
        label = "GlmnetOptions",
        engine = "glmnet"
    )
}

#' @rdname twasEngineOptions
#' @export
QggOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            y = "supplied from the data by the pipeline",
            # qgg's spellings do not line up with the wrapper's: the
            # genotype matrix is `W` here, and matching only the names that
            # happen to agree on both sides is how `W` was missed.
            W = "supplied from the data by the pipeline",
            stat = "supplied from the data by the pipeline",
            Glist = "supplied from the data by the pipeline",
            LD = "supplied from the data by the pipeline",
            method = "the wrapper's `method` (the function called)",
            nit = "the wrapper's `nIter`",
            nburn = "the wrapper's `nBurnin`",
            nthin = "the wrapper's `thin`"
        ),
        "QggOptions"
    )
    .configRefuseUnavailable(
        extra,
        c(
            # qgg's `X` is its covariates.
            X = str_c(
                "covariates are regressed out of the phenotype before ",
                "weights are fitted, so none reach the engine"
            )
        ),
        "QggOptions"
    )
    .newMethodOptions(
        "qgg::gbayes",
        defaults = list(),
        extra = list(...),
        label = "QggOptions",
        engine = "qgg"
    )
}

#' @rdname twasEngineOptions
#' @export
BglrOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            y = "supplied from the data by the pipeline",
            ETA = "supplied from the data by the pipeline",
            nIter = "the wrapper's `nIter`",
            burnIn = "the wrapper's `nBurnin`",
            thin = "the wrapper's `thin`"
        ),
        "BglrOptions"
    )
    .newMethodOptions(
        "BGLR::BGLR",
        defaults = list(),
        extra = list(...),
        label = "BglrOptions",
        engine = "bglr"
    )
}

#' @rdname twasEngineOptions
#' @export
HibayesOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            formula = "supplied from the data by the pipeline",
            data = "supplied from the data by the pipeline",
            M = "supplied from the data by the pipeline",
            M.id = "supplied from the data by the pipeline",
            method = "the wrapper's `method` (the function called)",
            niter = "the wrapper's `nIter`",
            nburn = "the wrapper's `nBurnin`",
            thin = "the wrapper's `thin`"
        ),
        "HibayesOptions"
    )
    .configRefuseUnavailable(extra, .hibayesUnavailable(), "HibayesOptions")
    # The individual-level entry point. Its summary-statistics twin,
    # HibayesRssOptions(), is separate because hibayes::sbrm() takes
    # different arguments, so one record checked against both would accept
    # a name meant for the other path.
    .newMethodOptions(
        "hibayes::ibrm",
        defaults = list(),
        extra = list(...),
        label = "HibayesOptions",
        engine = "hibayes"
    )
}

#' @rdname twasEngineOptions
#' @export
HibayesRssOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            sumstat = "supplied from the data by the pipeline",
            ldm = "supplied from the data by the pipeline",
            method = "the wrapper's `method` (the function called)",
            niter = "the wrapper's `nIter`",
            nburn = "the wrapper's `nBurnin`",
            thin = "the wrapper's `thin`"
        ),
        "HibayesRssOptions"
    )
    .configRefuseUnavailable(extra, .hibayesUnavailable(), "HibayesRssOptions")
    # The summary-statistics twin of HibayesOptions(): hibayes::sbrm()
    # takes different arguments from hibayes::ibrm().
    .newMethodOptions(
        "hibayes::sbrm",
        defaults = list(),
        extra = list(...),
        label = "HibayesRssOptions",
        engine = "hibayesRss"
    )
}

#' @rdname twasEngineOptions
#' @export
BayesAlphabetRssOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            bhat = "supplied from the data by the pipeline",
            F = "supplied from the data by the pipeline",
            n = "supplied from the data by the pipeline",
            nPerVariant = "supplied from the data by the pipeline",
            method = "the wrapper's `method` (the function called)",
            nIter = "the wrapper's `nIter`",
            nBurnin = "the wrapper's `nBurnin`",
            thin = "the wrapper's `thin`",
            pi = "the wrapper's `pi`",
            shrinkCorr = "the wrapper's `shrinkCorr`"
        ),
        "BayesAlphabetRssOptions"
    )
    .newMethodOptions(
        "pecotmr::bayesAlphabetRss",
        defaults = list(),
        extra = list(...),
        label = "BayesAlphabetRssOptions",
        engine = "bayesAlphabetRss"
    )
}

# hibayes arguments pecotmr sets itself and does not expose, with why.
# @noRd
.hibayesUnavailable <- function() {
    c(
        seed = str_c(
            "it is drawn from R's random number generator, so set.seed() ",
            "sets it"
        ),
        threads = str_c(
            "hibayes runs on one thread, as weight fitting already runs ",
            "inside BiocParallel workers"
        )
    )
}

#' @rdname twasEngineOptions
#' @export
DprOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            y = "supplied from the data by the pipeline",
            x = "supplied from the data by the pipeline",
            w = "supplied from the data by the pipeline",
            fitting_method = "the wrapper's `fittingMethod`"
        ),
        "DprOptions"
    )
    .newMethodOptions(
        "RcppDPR::fit_model",
        defaults = list(),
        extra = list(...),
        label = "DprOptions",
        engine = "dpr"
    )
}

#' @rdname twasEngineOptions
#' @export
NcvregOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            X = "supplied from the data by the pipeline",
            y = "supplied from the data by the pipeline",
            nfolds = "the wrapper's `nfolds`",
            penalty = "the calling wrapper (scadWeights, mcpWeights)"
        ),
        "NcvregOptions"
    )
    .newMethodOptions(
        "ncvreg::cv.ncvreg",
        defaults = list(),
        extra = list(...),
        label = "NcvregOptions",
        engine = "ncvreg"
    )
}

#' @rdname twasEngineOptions
#' @export
L0learnOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            x = "supplied from the data by the pipeline",
            y = "supplied from the data by the pipeline",
            nFolds = "the wrapper's `nFolds`",
            penalty = "the wrapper's `penalty`"
        ),
        "L0learnOptions"
    )
    .newMethodOptions(
        "L0Learn::L0Learn.cvfit",
        defaults = list(),
        extra = list(...),
        label = "L0learnOptions",
        engine = "l0learn"
    )
}

#' @rdname twasEngineOptions
#' @export
MrashOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            X = "supplied from the data by the pipeline",
            y = "supplied from the data by the pipeline"
        ),
        "MrashOptions"
    )
    .newMethodOptions(
        "susieR::mr.ash",
        defaults = list(),
        extra = list(...),
        label = "MrashOptions",
        engine = "mrash"
    )
}

#' @rdname twasEngineOptions
#' @export
MrashRssOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            bhat = "supplied from the data by the pipeline",
            shat = "supplied from the data by the pipeline",
            z = "supplied from the data by the pipeline",
            R = "supplied from the data by the pipeline",
            n = "supplied from the data by the pipeline",
            s0 = "the wrapper's `s0`",
            w0 = "the wrapper's `w0`",
            # The engine writes these two in snake_case while the wrapper
            # spells them camelCase, so they read as unrelated names; both
            # are still fed straight from a formal.
            var_y = "the wrapper's `varY`",
            sigma2_e = "the wrapper's `sigma2E`"
        ),
        "MrashRssOptions"
    )
    # The summary-statistics twin of MrashOptions(). Separate because the two
    # engines take different arguments entirely -- mr.ash takes (X, y, ...),
    # mr.ash.rss takes (bhat, shat, R, ...) -- so one record validating
    # against the union of both would accept an individual-level name on a
    # summary-statistics run, and the reverse.
    .newMethodOptions(
        "susieR::mr.ash.rss",
        defaults = list(),
        extra = list(...),
        label = "MrashRssOptions",
        engine = "mrashRss"
    )
}

#' @rdname twasEngineOptions
#' @export
PrsCsOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            bhat = "supplied from the data by the pipeline",
            R = "supplied from the data by the pipeline",
            n = "supplied from the data by the pipeline"
        ),
        "PrsCsOptions"
    )
    .newMethodOptions(
        "pecotmr::prsCs",
        defaults = list(),
        extra = list(...),
        label = "PrsCsOptions",
        engine = "prsCs"
    )
}

#' @rdname twasEngineOptions
#' @export
LassosumOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            bhat = "supplied from the data by the pipeline",
            R = "supplied from the data by the pipeline",
            n = "supplied from the data by the pipeline"
        ),
        "LassosumOptions"
    )
    .newMethodOptions(
        "pecotmr::lassosumRss",
        defaults = list(),
        extra = list(...),
        label = "LassosumOptions",
        engine = "lassosum"
    )
}

#' @rdname twasEngineOptions
#' @export
BigsnprOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            corr = "supplied from the data by the pipeline",
            df_beta = "supplied from the data by the pipeline",
            ind.corr = "supplied from the data by the pipeline",
            delta = "lassosumRssWeights()'s `s`, as delta = s / (1 - s)",
            ncores = str_c(
                "one per window: BiocParallel already parallelizes over ",
                "windows"
            )
        ),
        "BigsnprOptions"
    )
    .newMethodOptions(
        "bigsnpr::snp_lassosum2",
        defaults = list(),
        extra = list(...),
        label = "BigsnprOptions",
        engine = "bigsnpr"
    )
}

#' @rdname twasEngineOptions
#' @export
PenalizedRssOptions <- function(...) {
    extra <- list(...)
    # The penalty and its shape belong to the calling wrapper: each one fixes
    # its penalty, and those that let it vary take `alpha` / `gamma` /
    # `lambda0` / `lambda2` / `maxSwaps` as their own arguments. `lambda`
    # stays settable here, the only route to it for the wrappers without one.
    .configRefuseOwned(
        extra,
        c(
            bhat = "supplied from the data by the pipeline",
            R = "supplied from the data by the pipeline",
            n = "supplied from the data by the pipeline",
            penalty = str_c(
                "the calling wrapper, which fixes it (scadRssWeights, ",
                "mcpRssWeights, ...)"
            ),
            alpha = str_c(
                "the calling wrapper's `alpha` (enetRssWeights, ",
                "scadRssWeights, mcpRssWeights); ridgeRssWeights fixes it at 0"
            ),
            gamma = str_c(
                "the calling wrapper's `gamma` (scadRssWeights, ",
                "mcpRssWeights)"
            ),
            lambda0 = "l0learnRssWeights' own `lambda0`",
            lambda2 = "l0learnRssWeights' own `lambda2`",
            maxSwaps = "l0learnRssWeights' own `maxSwaps`"
        ),
        "PenalizedRssOptions"
    )
    .newMethodOptions(
        "pecotmr::penalizedRss",
        defaults = list(),
        extra = list(...),
        label = "PenalizedRssOptions",
        engine = "penalizedRss"
    )
}

#' @rdname twasEngineOptions
#' @export
SdprOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            bhat = "supplied from the data by the pipeline",
            R = "supplied from the data by the pipeline",
            n = "supplied from the data by the pipeline"
        ),
        "SdprOptions"
    )
    .newMethodOptions(
        "pecotmr::sdpr",
        defaults = list(),
        extra = list(...),
        label = "SdprOptions",
        engine = "sdpr"
    )
}

#' @rdname twasEngineOptions
#' @export
MrmashOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            X = "supplied from the data by the pipeline",
            Y = "supplied from the data by the pipeline",
            S0 = "the pipeline's prior matrices",
            mu1_init = "the wrapper's `bInit`",
            standardize = "the wrapper's `standardize`",
            update_V_method = "the wrapper's `updateVMethod`"
        ),
        "MrmashOptions"
    )
    .configRefuseUnavailable(extra, .mrmashUnavailable(), "MrmashOptions")
    .newMethodOptions(
        "mr.mashr::mr.mash",
        defaults = list(),
        extra = list(...),
        label = "MrmashOptions",
        engine = "mrmash",
        # The bag is split before the engine sees it: .mrmashSplitPriorArgs
        # routes the prior-covariance names to MrmashPriorParam() and the rest
        # to mr.mash. So the accepted set is both, read from the live
        # formals of each rather than transcribed.
        accepted = unique(c(
            .engineAcceptedNames("mr.mashr::mr.mash"),
            names(formals(MrmashPriorParam))
        ))
    )
}

#' @rdname twasEngineOptions
#' @export
MrmashRssOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            Bhat = "supplied from the data by the pipeline",
            Shat = "supplied from the data by the pipeline",
            R = "supplied from the data by the pipeline",
            n = "supplied from the data by the pipeline",
            covY = "supplied from the data by the pipeline",
            V = "supplied from the data by the pipeline",
            S0 = "the pipeline's prior matrices",
            w0 = "the pipeline's prior weights"
        ),
        "MrmashRssOptions"
    )
    .configRefuseUnavailable(
        extra,
        .mrmashUnavailable(),
        "MrmashRssOptions"
    )
    # The summary-statistics twin of MrmashOptions(). mr.mash takes
    # (X, Y, S0, w0, V, mu1_init, ...) and mr.mash.rss takes
    # (Bhat, Shat, Z, R, covY, n, ...), so the two cannot share one
    # validated record without each accepting the other's names.
    .newMethodOptions(
        "mr.mashr::mr.mash.rss",
        defaults = list(),
        extra = list(...),
        label = "MrmashRssOptions",
        engine = "mrmashRss",
        # Same prior-name passthrough as the individual path:
        # .mrmashSplitPriorArgs routes those to MrmashPriorParam().
        accepted = unique(c(
            .engineAcceptedNames("mr.mashr::mr.mash.rss"),
            names(formals(MrmashPriorParam))
        ))
    )
}

# mr.mash arguments pecotmr sets itself and does not expose, with why. Its
# RcppParallel workers run OpenMP BLAS, which registers them with the OpenMP
# runtime. SparseArray pauses that runtime after every call, and on macOS
# (LLVM libomp) resuming it spins forever on such workers, so the next
# parallel region after any SparseArray call hangs R. One thread has no
# workers.
# @noRd
.mrmashUnavailable <- function() {
    c(
        nthreads = str_c(
            "mr.mash runs on one thread, as its worker threads can leave ",
            "OpenMP in a state that hangs R on macOS"
        )
    )
}

#' Fit a cross-validated glmnet model
#'
#' The single interface to \code{glmnet::cv.glmnet}, which every glmnet fit in
#' the package goes through. It returns the \emph{fit}, not weights:
#' \code{\link{glmnetWeights}} extracts the \code{lambda.min} coefficient
#' vector from it, while the \code{mr.mash} initialisation paths need the
#' coefficient matrix, the intercept row, the lambda sequence and
#' \code{predict()} on held-out rows --- the same fit-versus-extract split as
#' \code{\link{fitMvsusie}} and \code{\link{mvsusieWeights}}.
#'
#' The four knobs below are formals rather than \code{\link{GlmnetOptions}}
#' entries because they are what the package's glmnet call sites genuinely
#' differ on, and leaving them implicit is how three different conventions
#' came to coexist. Each defaults to \code{NULL}, meaning the argument is not
#' passed at all and \code{cv.glmnet}'s own default stands; a value supplied
#' through \code{methodArgs} still wins over the caller's, so none of them is
#' refused by the constructor.
#'
#' @param x Numeric design matrix (samples x features).
#' @param y Response: a vector, or a matrix for \code{family = "mgaussian"}.
#' @param alpha Elastic-net mixing parameter in \code{[0, 1]}: \code{1} =
#'   lasso, \code{0} = ridge, \code{0.5} = elastic net. The penalty family is
#'   the calling wrapper's to choose, which is why \code{GlmnetOptions()}
#'   refuses it.
#' @param family glmnet's response family, e.g. \code{"gaussian"} or
#'   \code{"mgaussian"} for a matrix response. \code{NULL} leaves glmnet's
#'   own default.
#' @param standardize Logical; standardize the columns of \code{x} before
#'   fitting. \code{NULL} leaves glmnet's own default (\code{TRUE}).
#' @param intercept Logical; fit an intercept. \code{NULL} leaves glmnet's
#'   own default (\code{TRUE}).
#' @param nfolds Number of cross-validation folds. \code{NULL} leaves
#'   \code{cv.glmnet}'s own default (10).
#' @param methodArgs Further arguments forwarded to
#'   \code{glmnet::cv.glmnet}, built with \code{\link{GlmnetOptions}}. A bare
#'   list is refused, since it cannot be checked against the engine.
#' @return The fit object returned by \code{glmnet::cv.glmnet}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' fit <- fitGlmnetCv(X, y, alpha = 0.5, nfolds = 5)
#' coef(fit, s = "lambda.min")[1:5]
#' @export
fitGlmnetCv <- function(
    x,
    y,
    alpha,
    family = NULL,
    standardize = NULL,
    intercept = NULL,
    nfolds = NULL,
    methodArgs = GlmnetOptions()
) {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        abort("Package 'glmnet' is required for this function.")
    }
    .assertMethodOptions(methodArgs, "GlmnetOptions", "methodArgs")
    # compact() is what makes a NULL formal mean "do not pass it": an
    # explicit NULL would reach cv.glmnet and override its default with
    # nothing. methodArgs is spliced last, so a user value wins over the
    # caller's -- these are literal defaults, not pecotmr-owned arguments.
    callArgs <- list_modify(
        c(
            list(x = x, y = y, alpha = alpha),
            compact(list(
                family = family,
                standardize = standardize,
                intercept = intercept,
                nfolds = nfolds
            ))
        ),
        !!!as.list(methodArgs)
    )
    exec(glmnet::cv.glmnet, !!!callArgs)
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
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights: \code{"none"} (default) attaches
#'   nothing; \code{"slim"} the selected fit, as a list of its coefficients
#'   (\code{coef}, over the variants fitted), \code{intercept},
#'   \code{lambda} and cross-validated error there (\code{cvm},
#'   \code{cvsd}); \code{"full"} the whole \code{cv.glmnet} object, with
#'   every \code{lambda}'s coefficients. Either names the variants it was
#'   fitted on.
#' @param methodArgs Additional arguments forwarded to
#'   \code{glmnet::cv.glmnet}, built with \code{\link{GlmnetOptions}}. A bare
#'   list is refused, since it cannot be checked against the engine's formals.
#' @return Numeric vector of length \code{ncol(X)} of per-variant weights,
#'   with the fit attached as attribute \code{"fit"} unless
#'   \code{fitRetention = "none"}.
#' @importFrom stats coef
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' glmnetWeights(X, y, alpha = 0.5)
#' @export
glmnetWeights <- function(
    X,
    y,
    alpha,
    fitRetention = c("none", "slim", "full"),
    methodArgs = GlmnetOptions()
) {
    fitRetention <- arg_match(fitRetention)
    eff.wgt <- matrix(0, ncol = 1, nrow = ncol(X))
    keep <- .dropZeroVariance(X, "glmnetWeights")
    # `alpha` is the wrapper's own dial and is refused by GlmnetOptions()
    # itself, so there is no second check here. The three settings below are
    # literal defaults and stay overridable through `methodArgs`; the glmnet
    # requirement check and the constructor check live in fitGlmnetCv().
    enet <- fitGlmnetCv(
        x = X[, keep, drop = FALSE],
        y = y,
        alpha = alpha,
        standardize = FALSE,
        intercept = TRUE,
        nfolds = 5,
        methodArgs = methodArgs
    )
    coefs <- as.numeric(coef(enet, s = "lambda.min"))
    weights <- replace(eff.wgt, keep, coefs[-1L])
    if (identical(fitRetention, "none")) {
        return(weights)
    }
    fit <- if (identical(fitRetention, "full")) {
        enet
    } else {
        .glmnetSlimFit(enet, coefs)
    }
    `attr<-`(weights, "fit", .withFitVariantIds(fit, colnames(X)[keep]))
}

# The fit cv.glmnet selected: its coefficients over the variants fitted, the
# intercept, its lambda and the cross-validated error there.
# @noRd
.glmnetSlimFit <- function(cvFit, coefs) {
    at <- match(cvFit$lambda.min, cvFit$lambda)
    list(
        coef = coefs[-1L],
        intercept = coefs[[1L]],
        lambda = cvFit$lambda.min,
        cvm = cvFit$cvm[[at]],
        cvsd = cvFit$cvsd[[at]]
    )
}

#' Compute TWAS weights via elastic net (glmnet)
#'
#' Convenience wrapper for \code{\link{glmnetWeights}} at a mixing between the
#' lasso and ridge, \code{alpha = 0.5} unless set.
#' @param alpha Elastic-net mixing in \code{[0, 1]}: the L1 share of the
#'   penalty. Default 0.5.
#' @inheritParams glmnetWeights
#' @return Numeric vector of length \code{ncol(X)} of per-variant weights,
#'   with the fit attached as attribute \code{"fit"} unless
#'   \code{fitRetention = "none"}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' enetWeights(X, y)
#' # A larger ridge share:
#' enetWeights(X, y, alpha = 0.2)
#' @importFrom checkmate assertNumber
#' @export
enetWeights <- function(
    X,
    y,
    alpha = 0.5,
    fitRetention = c("none", "slim", "full"),
    methodArgs = GlmnetOptions()
) {
    .assertMethodOptions(methodArgs, "GlmnetOptions", "methodArgs")
    assertNumber(alpha, lower = 0, upper = 1)
    glmnetWeights(
        X,
        y,
        alpha,
        fitRetention = arg_match(fitRetention),
        methodArgs = methodArgs
    )
}

#' Compute TWAS weights via ridge regression (glmnet, alpha = 0)
#'
#' Convenience wrapper for \code{\link{glmnetWeights}} with \code{alpha = 0}.
#' @inheritParams glmnetWeights
#' @return Numeric vector of length \code{ncol(X)} of per-variant weights,
#'   with the fit attached as attribute \code{"fit"} unless
#'   \code{fitRetention = "none"}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ridgeWeights(X, y)
#' @export
ridgeWeights <- function(
    X,
    y,
    fitRetention = c("none", "slim", "full"),
    methodArgs = GlmnetOptions()
) {
    .assertMethodOptions(methodArgs, "GlmnetOptions", "methodArgs")
    glmnetWeights(
        X,
        y,
        0,
        fitRetention = arg_match(fitRetention),
        methodArgs = methodArgs
    )
}

#' Compute TWAS weights via lasso (glmnet, alpha = 1)
#'
#' Convenience wrapper for \code{\link{glmnetWeights}} with \code{alpha = 1}.
#' @inheritParams glmnetWeights
#' @return Numeric vector of length \code{ncol(X)} of per-variant weights,
#'   with the fit attached as attribute \code{"fit"} unless
#'   \code{fitRetention = "none"}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' lassoWeights(X, y)
#' @export
lassoWeights <- function(
    X,
    y,
    fitRetention = c("none", "slim", "full"),
    methodArgs = GlmnetOptions()
) {
    .assertMethodOptions(methodArgs, "GlmnetOptions", "methodArgs")
    glmnetWeights(
        X,
        y,
        1,
        fitRetention = arg_match(fitRetention),
        methodArgs = methodArgs
    )
}

#' Compute Weights Using mr.ash Shrinkage
#'
#' This function fits the `mr.ash` model (adaptive shrinkage regression) to
#' estimate weights for a given set of predictors and response. The prior is
#' \code{susieR::mr.ash}'s own default mixture grid unless \code{sa2} is
#' supplied through \code{methodArgs}, and custom initial beta values are
#' accepted.
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
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights: \code{"none"} (default) attaches
#'   nothing, \code{"slim"} the fitted model without its copy of the data
#'   (\code{data}) or the order it updated the variants in each iteration
#'   (\code{update.order}), \code{"full"} the whole fitted model. Either
#'   names the variants it was fitted on.
#' @param methodArgs Options forwarded to \code{susieR::mr.ash}, built
#'   with \code{\link{MrashOptions}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of weights, one per variant (column of \code{X});
#'   zero-variance columns receive weight 0. Unless
#'   \code{fitRetention = "none"} the fitted \code{mr.ash} object is
#'   attached as attribute \code{"fit"}.
#' @export
mrashWeights <- function(
    X,
    y,
    fitRetention = c("none", "slim", "full"),
    methodArgs = MrashOptions()
) {
    .assertMethodOptions(methodArgs, "MrashOptions", "methodArgs")
    fitRetention <- arg_match(fitRetention)
    keep <- .dropZeroVariance(X, "mrashWeights")
    XKeep <- X[, keep, drop = FALSE]
    # as.list(): list_assign() rejects a MethodOptions outright as `.x`, the
    # same S4-SimpleList seam that makes c() append it instead of splicing.
    argsList <- list_assign(
        as.list(methodArgs),
        !!!.mrashBetaInit(methodArgs, XKeep, y, X, keep)
    )
    fit.mr.ash <- exec(mr.ash, X = XKeep, y = y, !!!argsList)
    # Zero-variance columns were never fitted and keep a zero weight.
    eff.wgt <- replace(
        rep(0, ncol(X)),
        keep,
        predict(fit.mr.ash, type = "coefficients")[-1]
    )
    if (identical(fitRetention, "none")) {
        return(eff.wgt)
    }
    kept <- if (identical(fitRetention, "full")) {
        fit.mr.ash
    } else {
        .withoutFields(fit.mr.ash, c("data", "update.order"))
    }
    `attr<-`(eff.wgt, "fit", .withFitVariantIds(kept, colnames(XKeep)))
}
# --- Bayesian alphabet ------------------------------------------------------
#
# Each alphabet method runs on several engines (backends), so its weight
# function takes `backend` and the settings every engine shares under one
# spelling -- nIter, nBurnin, thin and pi -- which each engine's adapter
# renames to its own. .alphabetBackends is the one record of which engines a
# method has on each data path and which Options constructor configures each.

# The engines each Bayesian alphabet method runs on, per data path, with the
# Options constructor that configures each. The first is the default.
# @noRd
.alphabetBackends <- list(
    individual = list(
        bayesA = c(qgg = "QggOptions", hibayes = "HibayesOptions"),
        bayesB = c(bglr = "BglrOptions", hibayes = "HibayesOptions"),
        bayesC = c(qgg = "QggOptions", hibayes = "HibayesOptions"),
        bayesL = c(
            qgg = "QggOptions",
            bglr = "BglrOptions",
            hibayes = "HibayesOptions"
        ),
        bayesN = c(qgg = "QggOptions", hibayes = "HibayesOptions"),
        bayesR = c(qgg = "QggOptions", hibayes = "HibayesOptions")
    ),
    # On summary statistics the first is the default for an LD matrix; a
    # genotype reference defaults to pecotmr (.rssBackend()).
    sumstat = list(
        bayesA = c(
            qgg = "QggOptions",
            hibayes = "HibayesRssOptions",
            pecotmr = "BayesAlphabetRssOptions"
        ),
        bayesB = c(
            hibayes = "HibayesRssOptions",
            pecotmr = "BayesAlphabetRssOptions"
        ),
        bayesC = c(
            qgg = "QggOptions",
            hibayes = "HibayesRssOptions",
            pecotmr = "BayesAlphabetRssOptions"
        ),
        bayesL = c(
            hibayes = "HibayesRssOptions",
            qgg = "QggOptions",
            pecotmr = "BayesAlphabetRssOptions"
        ),
        bayesN = c(
            qgg = "QggOptions",
            hibayes = "HibayesRssOptions",
            pecotmr = "BayesAlphabetRssOptions"
        ),
        bayesR = c(
            qgg = "QggOptions",
            hibayes = "HibayesRssOptions",
            pecotmr = "BayesAlphabetRssOptions"
        )
    )
)

# Every method with several engines on a data path: the Bayesian alphabet
# (.alphabetBackends), and on summary statistics the lasso and LDpred2. Each
# engine maps to its Options constructor, NA where it takes none. On summary
# statistics the first engine is an LD matrix's default (.rssBackend()).
# @noRd
.methodBackends <- list(
    individual = .alphabetBackends$individual,
    sumstat = c(
        .alphabetBackends$sumstat,
        list(
            lasso = c(
                pecotmr = "LassosumOptions",
                bigsnpr = "BigsnprOptions"
            ),
            ldpred2Auto = c(bigsnpr = NA_character_, pecotmr = NA_character_),
            ldpred2Grid = c(bigsnpr = NA_character_, pecotmr = NA_character_)
        )
    )
)

# The summary-statistics methods whose pecotmr engine reads only a genotype
# reference's factor of the LD: the alphabet's sampler and LDpred2's.
# @noRd
.rssFactorOnlyMethods <- function() {
    c(names(.alphabetBackends$sumstat), "ldpred2Auto", "ldpred2Grid")
}

# Each engine's name for each method. hibayes takes the variants that update
# pi, as qgg and BGLR do by default.
# @noRd
.alphabetEngineMethods <- list(
    qgg = c(
        bayesA = "bayesA",
        bayesC = "bayesC",
        bayesL = "bayesL",
        bayesN = "bayesN",
        bayesR = "bayesR"
    ),
    bglr = c(bayesB = "BayesB", bayesL = "BL"),
    hibayes = c(
        bayesA = "BayesA",
        bayesB = "BayesBpi",
        bayesC = "BayesCpi",
        bayesL = "BayesL",
        bayesN = "BayesRR",
        bayesR = "BayesR"
    )
)

# The adapter that runs each engine, per data path. An individual-level
# adapter returns the engine's fitted object, which .alphabetPosterior()
# reads.
# @noRd
.alphabetFitters <- list(
    individual = c(
        qgg = ".alphabetQggFit",
        bglr = ".alphabetBglrFit",
        hibayes = ".alphabetHibayesFit"
    ),
    sumstat = c(
        qgg = ".alphabetQggRssFit",
        hibayes = ".alphabetHibayesRssFit",
        pecotmr = ".alphabetPecotmrRssFit"
    )
)

# The prior proportion of non-null effects, for the methods that have one.
# @noRd
.alphabetPiDefaults <- c(bayesB = 0.2, bayesC = 0.1)

#' Bayesian Alphabet TWAS Weights
#'
#' Fits a Bayesian linear regression of \code{y} on \code{X} under one of
#' the "Bayesian alphabet" priors and returns the posterior mean effects.
#' Each method runs on several engines, chosen by \code{backend}:
#' \itemize{
#'   \item \code{bayesA} (a scaled-t prior), \code{bayesC} (spike and slab
#'     with one slab variance), \code{bayesN} (Gaussian, ridge-like) and
#'     \code{bayesR} (a mixture of normals): \code{"qgg"} (default,
#'     \code{qgg::gbayes}) or \code{"hibayes"} (\code{hibayes::ibrm}).
#'   \item \code{bayesB} (spike and slab with a variance per variant):
#'     \code{"bglr"} (default, \code{BGLR::BGLR}) or \code{"hibayes"}.
#'   \item \code{bayesL} (the Bayesian LASSO): \code{"qgg"} (default),
#'     \code{"bglr"} or \code{"hibayes"}, the last two in Park and
#'     Casella's form.
#' }
#' The same settings reach every engine: \code{nIter}, \code{nBurnin},
#' \code{thin} and \code{pi} are renamed to each engine's own spelling, and
#' everything else goes through \code{methodArgs}, built with the backend's
#' own constructor.
#'
#' @param X A numeric matrix of genotypes (samples x variants).
#' @param y A numeric vector of phenotypes, one per row of \code{X}.
#' @param method The prior: \code{"bayesA"}, \code{"bayesB"},
#'   \code{"bayesC"}, \code{"bayesL"}, \code{"bayesN"} or \code{"bayesR"}.
#' @param backend The engine, one the method has (see Description);
#'   \code{NULL} picks the method's default.
#' @param pi For \code{bayesB} and \code{bayesC}, the prior proportion of
#'   non-null effects, where each engine's estimate of it starts; \code{NULL}
#'   gives 0.2 for \code{bayesB} and 0.1 for \code{bayesC}. Refused for the
#'   other methods.
#' @param nIter Integer. Total number of MCMC iterations, burn-in included,
#'   on every engine. Default \code{5000}.
#' @param nBurnin Integer. Number of burn-in iterations discarded; less than
#'   \code{nIter}. Default \code{1000}.
#' @param thin Integer. Thinning interval for retained MCMC samples. Default
#'   \code{5}.
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights: \code{"none"} (default) attaches
#'   nothing; \code{"slim"} the posterior summaries, as a list of the
#'   effects' posterior means (\code{coef}, over the variants fitted) and,
#'   where the backend reports them, their posterior inclusion probabilities
#'   (\code{pip}); \code{"full"} the backend's whole fitted object. Either
#'   names the variants it was fitted on.
#' @param methodArgs Options for the backend's engine, built with its
#'   constructor: \code{\link{QggOptions}}, \code{\link{BglrOptions}} or
#'   \code{\link{HibayesOptions}}. \code{NULL} keeps the engine's defaults.
#'   A bare list is refused: it cannot be checked against the engine, so a
#'   misspelled option would be silently ignored.
#' @return A numeric vector of weights, one per variant (column of
#'   \code{X}); columns dropped for zero variance receive weight 0. Unless
#'   \code{fitRetention = "none"} the fit is attached as attribute
#'   \code{"fit"}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' bayesAlphabetWeights(X, y, "bayesR", nIter = 1000, nBurnin = 200)
#' @importFrom checkmate assertVector
#' @export
bayesAlphabetWeights <- function(
    X,
    y,
    method,
    backend = NULL,
    pi = NULL,
    nIter = 5000,
    nBurnin = 1000,
    thin = 5,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    method <- arg_match(method, names(.alphabetBackends$individual))
    fitRetention <- arg_match(fitRetention)
    backend <- .methodBackend(method, backend, "individual")
    .assertMcmcSchedule(nIter, nBurnin, thin)
    pi <- .alphabetPi(method, pi)
    methodArgs <- .backendOptions(
        methodArgs,
        .alphabetBackends$individual[[method]][[backend]],
        backend
    )
    assertVector(y, len = nrow(X))
    keep <- .dropZeroVariance(X, "bayesAlphabetWeights")
    fit <- exec(
        .alphabetFitters$individual[[backend]],
        X = X[, keep, drop = FALSE],
        y = y,
        method = method,
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        pi = pi,
        methodArgs = methodArgs
    )
    posterior <- .alphabetPosterior(fit, backend)
    weights <- replace(rep(0, ncol(X)), keep, posterior$coef)
    if (identical(fitRetention, "none")) {
        return(weights)
    }
    kept <- if (identical(fitRetention, "full")) fit else posterior
    `attr<-`(weights, "fit", .withFitVariantIds(kept, colnames(X)[keep]))
}

# A backend's posterior summaries over the variants it fitted: the effects'
# posterior means and, where it reports them, their posterior inclusion
# probabilities.
# @noRd
.alphabetPosterior <- function(fit, backend) {
    parts <- switch(
        backend,
        qgg = list(coef = fit$bm, pip = fit$dm),
        bglr = list(
            coef = pluck(fit, "ETA", 1L, "b"),
            pip = pluck(fit, "ETA", 1L, "d")
        ),
        hibayes = list(coef = fit$alpha, pip = fit$pip)
    )
    map(compact(parts), as.numeric)
}

# The backend `backend` names for `method` on one data path, or the method's
# default when NULL. Refuses an engine the method does not have there.
# @noRd
.methodBackend <- function(method, backend, path) {
    choices <- names(.methodBackends[[path]][[method]])
    if (is.null(backend)) {
        return(choices[[1L]])
    }
    if (length(backend) != 1L || !is_in(backend, choices)) {
        abort(glue(
            "`backend` for {method} must be one of ",
            "{str_flatten(choices, ', ')}."
        ))
    }
    backend
}

# The MCMC schedule every backend shares: nIter counts every iteration,
# burn-in included.
# @noRd
#' @importFrom checkmate assertCount
.assertMcmcSchedule <- function(nIter, nBurnin, thin) {
    assertCount(nIter, positive = TRUE)
    assertCount(nBurnin)
    assertCount(thin, positive = TRUE)
    if (nIter <= nBurnin) {
        abort("`nIter`, the total iterations, must exceed `nBurnin`.")
    }
    invisible(NULL)
}

# The starting prior proportion of non-null effects: `pi`, else the method's
# default, for the methods that have one; NULL for the others, which refuse
# one.
# @noRd
#' @importFrom checkmate assertNumber
.alphabetPi <- function(method, pi) {
    if (!is_in(method, names(.alphabetPiDefaults))) {
        if (!is.null(pi)) {
            abort(glue("`pi` applies to bayesB and bayesC, not {method}."))
        }
        return(NULL)
    }
    pi <- pi %||% .alphabetPiDefaults[[method]]
    assertNumber(pi, lower = 0, upper = 1)
    pi
}

# The record the backend's engine runs with: its constructor's defaults for
# NULL or an empty list, otherwise `methodArgs`, which must have been built
# by that constructor.
# @noRd
.backendOptions <- function(methodArgs, ctorName, backend) {
    ctor <- get(ctorName, envir = asNamespace("pecotmr"))
    if (is.null(methodArgs) || identical(methodArgs, list())) {
        return(ctor())
    }
    .assertMethodOptions(methodArgs, ctorName, "methodArgs")
    want <- metadata(ctor())$engine
    got <- metadata(methodArgs)$engine
    if (!identical(got, want)) {
        abort(glue(
            "`methodArgs` was built for engine '{got}', but backend ",
            "'{backend}' takes {ctorName}()."
        ))
    }
    methodArgs
}

# qgg's sampler on individual-level data. qgg's `X` is its covariates, which
# never reach it: they are regressed out of the phenotype beforehand. qgg
# runs nit + nburn iterations, so its `nit` counts the draws after the
# burn-in, where `nIter` counts them all.
# @noRd
.alphabetQggFit <- function(
    X,
    y,
    method,
    nIter,
    nBurnin,
    thin,
    pi,
    methodArgs
) {
    if (!requireNamespace("qgg", quietly = TRUE)) {
        abort("Package 'qgg' is required for backend = \"qgg\".")
    }
    if (!is.null(pi)) {
        methodArgs <- .methodOptionsWith(methodArgs, "QggOptions", pi = pi)
    }
    callArgs <- list_modify(
        list(
            y = y,
            W = X,
            method = .alphabetEngineMethods$qgg[[method]],
            nit = nIter - nBurnin,
            nburn = nBurnin,
            nthin = thin
        ),
        !!!as.list(methodArgs)
    )
    exec(qgg::gbayes, !!!callArgs)
}

# BGLR's sampler, through bglrWeights(), which keeps BGLR's files in a
# temporary directory; the fitted object is what it returns as the weights'
# fit. BayesB's prior inclusion probability belongs to the linear term rather
# than to BGLR() itself.
# @noRd
.alphabetBglrFit <- function(
    X,
    y,
    method,
    nIter,
    nBurnin,
    thin,
    pi,
    methodArgs
) {
    weights <- bglrWeights(
        X,
        y,
        model = .alphabetEngineMethods$bglr[[method]],
        nIter = nIter,
        burnIn = nBurnin,
        thin = thin,
        etaArgs = if (is.null(pi)) list() else list(probIn = pi),
        fitRetention = "full",
        methodArgs = methodArgs
    )
    attr(weights, "fit", exact = TRUE)
}

# hibayes' sampler on individual-level data. ibrm() matches the genotypes to
# the phenotypes by id and calls set.seed() itself, so it gets a seed drawn
# from R's generator and the caller's stream is put back afterwards.
# @noRd
.alphabetHibayesFit <- function(
    X,
    y,
    method,
    nIter,
    nBurnin,
    thin,
    pi,
    methodArgs
) {
    if (!requireNamespace("hibayes", quietly = TRUE)) {
        abort("Package 'hibayes' is required for backend = \"hibayes\".")
    }
    if (!is.null(pi)) {
        methodArgs <- .methodOptionsWith(
            methodArgs,
            "HibayesOptions",
            Pi = c(1 - pi, pi)
        )
    }
    ids <- as.character(seq_len(nrow(X)))
    callArgs <- list_modify(
        list(
            formula = y ~ 1,
            data = data.frame(id = ids, y = y),
            M = X,
            M.id = ids,
            method = .alphabetEngineMethods$hibayes[[method]],
            niter = nIter,
            nburn = nBurnin,
            thin = thin,
            seed = .sessionSeedDraw(),
            threads = 1L,
            verbose = FALSE
        ),
        !!!as.list(methodArgs)
    )
    withr::with_preserve_seed(exec(hibayes::ibrm, !!!callArgs))
}

#' @title BayesA TWAS weights (scaled-t prior)
#' @description Fits the BayesA prior, a normal on each effect whose
#'   variance is drawn per variant (a scaled-t), through
#'   \code{\link{bayesAlphabetWeights}}.
#' @inheritParams bayesAlphabetWeights
#' @param backend \code{"qgg"} (default) or \code{"hibayes"}.
#' @param methodArgs Options for the backend's engine:
#'   \code{\link{QggOptions}} for qgg, \code{\link{HibayesOptions}} for
#'   hibayes. \code{NULL} keeps the engine's defaults.
#' @return A numeric vector of weights, one per variant (column of
#'   \code{X}); columns dropped for zero variance receive weight 0.
#'   Unless \code{fitRetention = "none"} the fit is attached as
#'   attribute \code{"fit"}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' bayesAWeights(X, y)
#' @export
bayesAWeights <- function(
    X,
    y,
    backend = c("qgg", "hibayes"),
    nIter = 5000,
    nBurnin = 1000,
    thin = 5,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    bayesAlphabetWeights(
        X,
        y,
        "bayesA",
        backend = arg_match(backend),
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        fitRetention = arg_match(fitRetention),
        methodArgs = methodArgs
    )
}

#' @title BayesB TWAS weights (spike and slab, a variance per variant)
#' @description Fits the BayesB prior, a point mass at zero mixed with a
#'   scaled-t slab, through \code{\link{bayesAlphabetWeights}}.
#' @inheritParams bayesAlphabetWeights
#' @param backend \code{"bglr"} (default) or \code{"hibayes"}.
#' @param pi The prior proportion of non-null effects, where each engine's
#'   estimate of it starts. Default \code{0.2}.
#' @param methodArgs Options for the backend's engine:
#'   \code{\link{BglrOptions}} for BGLR, \code{\link{HibayesOptions}} for
#'   hibayes. \code{NULL} keeps the engine's defaults.
#' @return A numeric vector of weights, one per variant (column of
#'   \code{X}); columns dropped for zero variance receive weight 0.
#'   Unless \code{fitRetention = "none"} the fit is attached as
#'   attribute \code{"fit"}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' bayesBWeights(X, y)
#' @export
bayesBWeights <- function(
    X,
    y,
    backend = c("bglr", "hibayes"),
    pi = 0.2,
    nIter = 5000,
    nBurnin = 1000,
    thin = 5,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    bayesAlphabetWeights(
        X,
        y,
        "bayesB",
        backend = arg_match(backend),
        pi = pi,
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        fitRetention = arg_match(fitRetention),
        methodArgs = methodArgs
    )
}

#' @title BayesC TWAS weights (spike and slab, one slab variance)
#' @description Fits the BayesC prior, a point mass at zero mixed with a
#'   normal whose variance every variant shares, through
#'   \code{\link{bayesAlphabetWeights}}.
#' @inheritParams bayesAlphabetWeights
#' @param backend \code{"qgg"} (default) or \code{"hibayes"}.
#' @param pi The prior proportion of non-null effects, where each engine's
#'   estimate of it starts. Default \code{0.1}.
#' @param methodArgs Options for the backend's engine:
#'   \code{\link{QggOptions}} for qgg, \code{\link{HibayesOptions}} for
#'   hibayes. \code{NULL} keeps the engine's defaults.
#' @return A numeric vector of weights, one per variant (column of
#'   \code{X}); columns dropped for zero variance receive weight 0.
#'   Unless \code{fitRetention = "none"} the fit is attached as
#'   attribute \code{"fit"}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' bayesCWeights(X, y)
#' @export
bayesCWeights <- function(
    X,
    y,
    backend = c("qgg", "hibayes"),
    pi = 0.1,
    nIter = 5000,
    nBurnin = 1000,
    thin = 5,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    bayesAlphabetWeights(
        X,
        y,
        "bayesC",
        backend = arg_match(backend),
        pi = pi,
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        fitRetention = arg_match(fitRetention),
        methodArgs = methodArgs
    )
}

#' @title BayesL TWAS weights (the Bayesian LASSO)
#' @description Fits the Bayesian LASSO, a Laplace prior on each effect,
#'   through \code{\link{bayesAlphabetWeights}}. qgg's form is the default;
#'   BGLR's and hibayes' follow Park and Casella (2008).
#' @inheritParams bayesAlphabetWeights
#' @param backend \code{"qgg"} (default), \code{"bglr"} or
#'   \code{"hibayes"}.
#' @param methodArgs Options for the backend's engine:
#'   \code{\link{QggOptions}} for qgg, \code{\link{BglrOptions}} for BGLR,
#'   \code{\link{HibayesOptions}} for hibayes. \code{NULL} keeps the
#'   engine's defaults.
#' @return A numeric vector of weights, one per variant (column of
#'   \code{X}); columns dropped for zero variance receive weight 0.
#'   Unless \code{fitRetention = "none"} the fit is attached as
#'   attribute \code{"fit"}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' bayesLWeights(X, y)
#' @export
bayesLWeights <- function(
    X,
    y,
    backend = c("qgg", "bglr", "hibayes"),
    nIter = 5000,
    nBurnin = 1000,
    thin = 5,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    bayesAlphabetWeights(
        X,
        y,
        "bayesL",
        backend = arg_match(backend),
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        fitRetention = arg_match(fitRetention),
        methodArgs = methodArgs
    )
}

#' @title BayesN TWAS weights (Gaussian prior, ridge-like)
#' @description Fits a normal prior shared by every effect, whose posterior
#'   mean is the BLUP, a ridge regression, through
#'   \code{\link{bayesAlphabetWeights}}.
#' @inheritParams bayesAlphabetWeights
#' @param backend \code{"qgg"} (default) or \code{"hibayes"}.
#' @param methodArgs Options for the backend's engine:
#'   \code{\link{QggOptions}} for qgg, \code{\link{HibayesOptions}} for
#'   hibayes. \code{NULL} keeps the engine's defaults.
#' @return A numeric vector of weights, one per variant (column of
#'   \code{X}); columns dropped for zero variance receive weight 0.
#'   Unless \code{fitRetention = "none"} the fit is attached as
#'   attribute \code{"fit"}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' bayesNWeights(X, y)
#' @export
bayesNWeights <- function(
    X,
    y,
    backend = c("qgg", "hibayes"),
    nIter = 5000,
    nBurnin = 1000,
    thin = 5,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    bayesAlphabetWeights(
        X,
        y,
        "bayesN",
        backend = arg_match(backend),
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        fitRetention = arg_match(fitRetention),
        methodArgs = methodArgs
    )
}

#' @title BayesR TWAS weights (a mixture of normals)
#' @description Fits the BayesR prior, a mixture of a point mass at zero and
#'   normals of increasing variance, through
#'   \code{\link{bayesAlphabetWeights}}.
#' @inheritParams bayesAlphabetWeights
#' @param backend \code{"qgg"} (default) or \code{"hibayes"}.
#' @param methodArgs Options for the backend's engine:
#'   \code{\link{QggOptions}} for qgg, \code{\link{HibayesOptions}} for
#'   hibayes. \code{NULL} keeps the engine's defaults.
#' @return A numeric vector of weights, one per variant (column of
#'   \code{X}); columns dropped for zero variance receive weight 0.
#'   Unless \code{fitRetention = "none"} the fit is attached as
#'   attribute \code{"fit"}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' bayesRWeights(X, y)
#' @export
bayesRWeights <- function(
    X,
    y,
    backend = c("qgg", "hibayes"),
    nIter = 5000,
    nBurnin = 1000,
    thin = 5,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    bayesAlphabetWeights(
        X,
        y,
        "bayesR",
        backend = arg_match(backend),
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        fitRetention = arg_match(fitRetention),
        methodArgs = methodArgs
    )
}

# --- Bayesian alphabet on summary statistics --------------------------------
#
# qgg and hibayes read the LD as a matrix: given a genotype reference, the
# matrix is formed from the panel for them. pecotmr's own sampler
# (bayesAlphabetRss(), on rss_mcmc.h) reads a genotype reference's
# factor and never forms it, and runs on nothing else; it is the default
# there, while an LD matrix defaults to the first engine .alphabetBackends
# lists.

#' Bayesian Alphabet TWAS Weights From Summary Statistics
#'
#' Fits a Bayesian alphabet prior to marginal summary statistics and an LD
#' reference and returns the posterior mean effects. Each method runs on
#' several engines, chosen by \code{backend}:
#' \itemize{
#'   \item \code{"qgg"}: \code{qgg::gbayes} on the LD matrix, for
#'     \code{bayesA}, \code{bayesC}, \code{bayesL}, \code{bayesN} and
#'     \code{bayesR}.
#'   \item \code{"hibayes"}: \code{hibayes::sbrm} on the LD matrix, for
#'     every method.
#'   \item \code{"pecotmr"}: pecotmr's own sampler, on a genotype
#'     reference's factor of the LD without forming the matrix. Its priors
#'     follow qgg's for \code{bayesA}, \code{bayesC}, \code{bayesN} and
#'     \code{bayesR} and hibayes' for \code{bayesB}; \code{bayesL} takes
#'     Park and Casella's form (hibayes') or Legarra et al.'s (qgg's
#'     individual-level one), through
#'     \code{\link{BayesAlphabetRssOptions}}. By default it reads the
#'     statistics through SDPR's likelihood (Zhou and Zhao 2021), which
#'     drops the part of the marginal effects the reference LD cannot
#'     explain; \code{BayesAlphabetRssOptions(ldProjection = NULL)} turns
#'     that off (see \code{ldProjection} of the internal
#'     \code{bayesAlphabetRss()}).
#' }
#' \code{NULL} picks \code{"pecotmr"} for a genotype reference, and for an
#' LD matrix qgg (\code{bayesA}, \code{bayesC}, \code{bayesN},
#' \code{bayesR}) or hibayes (\code{bayesB}, \code{bayesL}). qgg and
#' hibayes run on a genotype reference through the LD matrix formed from
#' it; pecotmr's sampler does not run on an LD matrix.
#'
#' @param stat A list of summary statistics: \code{n} (per-variant sample
#'   sizes) and the marginal effects as \code{cor} (standardized), \code{z},
#'   \code{b} with \code{seb}, or \code{b} alone, taken as already
#'   standardized. The engines receive the standardized effects
#'   z / sqrt(n), with n the median of \code{n}.
#' @param ldReference The LD reference the statistics' variants are matched
#'   to by \code{stat$variantNames}: a genotype panel (see
#'   \code{\link{readGenotypes}}), a \code{\link{GenotypeMixture}} or an
#'   \code{\link{LdMatrix}}.
#' @param method The prior: \code{"bayesA"}, \code{"bayesB"},
#'   \code{"bayesC"}, \code{"bayesL"}, \code{"bayesN"} or \code{"bayesR"}.
#' @param backend The engine, one the method has (see Description);
#'   \code{NULL} picks the default for the LD reference.
#' @inheritParams bayesAlphabetWeights
#' @param shrinkCorr The LD shrinkage \eqn{c}: every engine reads
#'   \eqn{cR + (1 - c)I}, qgg and hibayes as the matrix they are given,
#'   pecotmr's sampler through the factor. Default 0.95, the LDpred2
#'   vignette's. A reference with fewer samples than the window has variants
#'   lacks directions the statistics have, and without shrinkage the
#'   samplers' effects can grow without bound along them (pecotmr's sampler
#'   also drops those directions by default, through SDPR's likelihood; see
#'   \code{\link{bayesAlphabetRssWeights}}). A small or
#'   mismatched reference does better with less, down to 0.4 as the LDpred2
#'   vignette advises; the study's own genotypes as the reference need none,
#'   up to 1.
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights: \code{"none"} (default) attaches
#'   nothing, \code{"slim"} and \code{"full"} attach the engine's estimates
#'   (\code{betaEst}, and where the engine gives them the inclusion
#'   probabilities \code{pip}, the non-null proportion \code{pi} and the
#'   residual variance), \code{backend} and the \code{shrinkCorr} used. No
#'   engine keeps trimmable intermediates, so the two retaining levels
#'   behave alike.
#' @param methodArgs Options for the backend's engine, built with its
#'   constructor: \code{\link{QggOptions}} for qgg,
#'   \code{\link{HibayesRssOptions}} for hibayes,
#'   \code{\link{BayesAlphabetRssOptions}} for pecotmr's sampler. \code{NULL}
#'   keeps the engine's defaults.
#' @return A numeric vector of the posterior mean effects on the
#'   standardized scale, one per variant. Unless \code{fitRetention =
#'   "none"}, the engine's estimates are attached as attribute \code{"fit"}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ss <- lapply(
#'   seq_len(ncol(X)), function(j) coef(summary(lm(y ~ X[, j])))[2, 1:2])
#' stat <- list(
#'   b = vapply(ss, `[`, numeric(1), 1L),
#'   seb = vapply(ss, `[`, numeric(1), 2L),
#'   n = rep(nrow(X), ncol(X)),
#'   variantNames = colnames(X)
#' )
#' w <- bayesAlphabetRssWeights(
#'   stat,
#'   LdMatrix(cor(X)),
#'   "bayesN",
#'   backend = "hibayes",
#'   nIter = 1000,
#'   nBurnin = 200
#' )
#' @export
bayesAlphabetRssWeights <- function(
    stat,
    ldReference,
    method,
    backend = NULL,
    pi = NULL,
    nIter = 5000,
    nBurnin = 1000,
    thin = 5,
    shrinkCorr = 0.95,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    fitRetention <- arg_match(fitRetention)
    method <- arg_match(method, names(.alphabetBackends$sumstat))
    backend <- .rssBackend(method, backend, ldReference)
    .assertMcmcSchedule(nIter, nBurnin, thin)
    .alphabetAssertSampleSize(stat$n, "bayesAlphabetRssWeights")
    .assertShrinkCorr(shrinkCorr)
    pi <- .alphabetPi(method, pi)
    methodArgs <- .backendOptions(
        methodArgs,
        .alphabetBackends$sumstat[[method]][[backend]],
        backend
    )
    fit <- exec(
        .alphabetFitters$sumstat[[backend]],
        stat = stat,
        ldReference = ldReference,
        method = method,
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        pi = pi,
        shrinkCorr = shrinkCorr,
        methodArgs = methodArgs
    )
    if (identical(fitRetention, "none")) {
        return(fit$betaEst)
    }
    `attr<-`(
        fit$betaEst,
        "fit",
        .withFitVariantIds(
            c(fit, list(shrinkCorr = shrinkCorr)),
            stat[["variantNames"]]
        )
    )
}

# The alphabet's samplers divide by n - 1 and the matrix engines' standard
# errors by n - 2, so the sample size (the median of `n`) must exceed 2.
# @noRd
.alphabetAssertSampleSize <- function(n, label) {
    assertNumeric(n, any.missing = FALSE, min.len = 1L)
    if (median(n) <= 2) {
        abort(glue("{label}: the sample size must exceed 2."))
    }
    invisible(NULL)
}

# The engine a summary-statistics method runs on: `backend`, or by default
# pecotmr where the method has it and a genotype reference supplies the LD,
# and otherwise the method's first engine. A pecotmr engine that reads only
# a factor of the LD (.rssFactorOnlyMethods()) refuses an LD matrix.
# @noRd
.rssBackend <- function(method, backend, ldReference) {
    choices <- names(.methodBackends$sumstat[[method]])
    genotypeRef <- .ldIsGenotypeRef(ldReference)
    if (is.null(backend)) {
        if (genotypeRef && is_in("pecotmr", choices)) {
            return("pecotmr")
        }
        return(choices[[1L]])
    }
    backend <- .methodBackend(method, backend, "sumstat")
    factorOnly <- is_in(method, .rssFactorOnlyMethods())
    if (identical(backend, "pecotmr") && factorOnly && !genotypeRef) {
        others <- setdiff(choices, "pecotmr")
        abort(glue(
            "backend = \"pecotmr\" runs on a genotype reference's factor of ",
            "the LD; with an LD matrix, use {str_flatten(others, ' or ')}."
        ))
    }
    backend
}

# The standardized statistics and the LD matrix the matrix engines read,
# formed from the panel for a genotype reference and shrunk to
# c R + (1 - c) I.
# @noRd
.alphabetRssDense <- function(stat, ldReference, shrinkCorr) {
    R <- .rssLdMatrix(stat, ldReference, "bayesAlphabetRssWeights")
    n <- median(stat$n)
    bhat <- .rssMarginalCor(stat, n = n, p = nrow(R))
    .rssValidateSource(bhat, list(R = R), n)
    list(
        R = .rssShrunkMatrix(R, shrinkCorr),
        bhat = bhat,
        seb = sqrt((1 - bhat^2) / (n - 2)),
        n = n,
        ids = stat$variantNames
    )
}

# qgg's summary-statistics sampler, on its in-memory LD: every variant's LD
# with every other, as the dense matrix gives it. qgg runs nit + nburn
# iterations, so nit is nIter net of the burn-in.
# @noRd
.alphabetQggRssFit <- function(
    stat,
    ldReference,
    method,
    nIter,
    nBurnin,
    thin,
    pi,
    shrinkCorr,
    methodArgs
) {
    if (!requireNamespace("qgg", quietly = TRUE)) {
        abort("Package 'qgg' is required for backend = \"qgg\".")
    }
    if (!is.null(pi)) {
        methodArgs <- .methodOptionsWith(methodArgs, "QggOptions", pi = pi)
    }
    inp <- .alphabetRssDense(stat, ldReference, shrinkCorr)
    qin <- .alphabetQggInputs(inp)
    callArgs <- list_modify(
        list(
            stat = qin$stat,
            Glist = qin$Glist,
            LD = qin$LD,
            method = .alphabetEngineMethods$qgg[[method]],
            nit = nIter - nBurnin,
            nburn = nBurnin,
            nthin = thin
        ),
        !!!as.list(methodArgs)
    )
    fit <- exec(qgg::gbayes, !!!callArgs)
    list(
        betaEst = unname(fit$stat[inp$ids, "bm"]),
        pip = unname(fit$stat[inp$ids, "dm"]),
        backend = "qgg"
    )
}

# qgg's summary statistics, map and LD for one window, built in memory: the
# statistics on the standardized scale (allele labels and frequencies are
# placeholders qgg carries but does not read), and each variant's LD with
# every other.
# @noRd
.alphabetQggInputs <- function(inp) {
    m <- length(inp$bhat)
    ids <- inp$ids
    pos <- seq_len(m)
    list(
        stat = data.frame(
            rsids = ids,
            chr = 1L,
            pos = pos,
            ea = "A",
            nea = "G",
            eaf = 0.5,
            b = inp$bhat,
            seb = inp$seb,
            p = 2 * pnorm(-abs(inp$bhat / inp$seb)),
            n = inp$n,
            row.names = ids
        ),
        Glist = list(
            nchr = 1L,
            rsidsLD = list(ids),
            rsids = list(ids),
            pos = list(pos),
            a1 = list(rep("A", m)),
            a2 = list(rep("G", m)),
            af = list(rep(0.5, m))
        ),
        LD = list(list(
            values = set_names(as.list(asplit(inp$R, 2L)), ids),
            indices = set_names(rep(list(seq_len(m) - 1L), m), ids)
        ))
    )
}

# hibayes' summary-statistics sampler on the LD matrix. Its genotype
# variances are the matrix's diagonal, so the correlation matrix and
# standardized statistics give standardized effects; it reads no allele
# frequencies, so the column carries a placeholder. sbrm() calls set.seed()
# itself, so it gets a seed drawn from R's generator and the caller's stream
# is put back afterwards.
# @noRd
.alphabetHibayesRssFit <- function(
    stat,
    ldReference,
    method,
    nIter,
    nBurnin,
    thin,
    pi,
    shrinkCorr,
    methodArgs
) {
    if (!requireNamespace("hibayes", quietly = TRUE)) {
        abort("Package 'hibayes' is required for backend = \"hibayes\".")
    }
    if (!is.null(pi)) {
        methodArgs <- .methodOptionsWith(
            methodArgs,
            "HibayesRssOptions",
            Pi = c(1 - pi, pi)
        )
    }
    inp <- .alphabetRssDense(stat, ldReference, shrinkCorr)
    callArgs <- list_modify(
        list(
            sumstat = data.frame(
                SNP = inp$ids,
                A1 = "A",
                A2 = "G",
                freq = 0.5,
                b = inp$bhat,
                se = inp$seb,
                p = 2 * pnorm(-abs(inp$bhat / inp$seb)),
                N = inp$n
            ),
            ldm = inp$R,
            method = .alphabetEngineMethods$hibayes[[method]],
            niter = nIter,
            nburn = nBurnin,
            thin = thin,
            seed = .sessionSeedDraw(),
            threads = 1L,
            verbose = FALSE
        ),
        !!!as.list(methodArgs)
    )
    fit <- withr::with_preserve_seed(exec(hibayes::sbrm, !!!callArgs))
    list(betaEst = as.numeric(fit$alpha), backend = "hibayes")
}

# pecotmr's sampler, on a genotype reference's factor of the LD.
# @noRd
.alphabetPecotmrRssFit <- function(
    stat,
    ldReference,
    method,
    nIter,
    nBurnin,
    thin,
    pi,
    shrinkCorr,
    methodArgs
) {
    ld <- .rssLdSource(stat, ldReference, "bayesAlphabetRssWeights")
    n <- median(stat$n)
    bhat <- .rssMarginalCor(stat, n = n, p = .ldSourceSize(ld))
    .rssValidateSource(bhat, ld, n)
    callArgs <- list_modify(
        list(
            bhat = bhat,
            F = ld$F,
            n = n,
            method = method,
            nIter = nIter,
            nBurnin = nBurnin,
            thin = thin,
            pi = pi,
            nPerVariant = stat$n,
            shrinkCorr = shrinkCorr
        ),
        !!!as.list(methodArgs)
    )
    c(exec(bayesAlphabetRss, !!!callArgs), list(backend = "pecotmr"))
}

#' pecotmr's Bayesian Alphabet Sampler on a Factor of the LD
#'
#' Internal: \code{\link{bayesAlphabetRssWeights}} is the public interface;
#' its \code{backend = "pecotmr"} runs this, and
#' \code{\link{BayesAlphabetRssOptions}} sets the arguments below other
#' than the data and the wrapper's own.
#'
#' The LD is a factor F of R = F'F / (k - 1), as a genotype reference gives
#' it, and is never formed: each effect's update costs O(k). The likelihood
#' is SBayes's, X'X = W^1/2 R W^1/2 and X'y = w * bhat for per-variant
#' weights w, read by default through SDPR's projection
#' (\code{ldProjection}). The priors are qgg's (CRAN 1.1.6) for
#' \code{bayesA}, \code{bayesC}, \code{bayesN} and \code{bayesR},
#' hibayes' for \code{bayesB} (as BayesBpi), and for \code{bayesL} either.
#'
#' @param bhat The standardized marginal effects.
#' @param F The factor of the LD (k x p).
#' @param n The sample size.
#' @param method The prior, as in \code{\link{bayesAlphabetRssWeights}}.
#' @param nIter,nBurnin,thin The chain's iterations (burn-in included),
#'   burn-in and thinning.
#' @param pi For \code{bayesB} and \code{bayesC}, where the estimate of the
#'   non-null proportion starts.
#' @param nPerVariant The per-variant sample sizes, read when
#'   \code{perVariantWeights = TRUE}.
#' @param bayesLForm For \code{bayesL} only: \code{"parkCasella"} (the
#'   default, hibayes' form, the Laplace prior scaled by the residual
#'   variance) or \code{"legarra"} (Legarra et al. 2011, qgg's
#'   individual-level form, unscaled).
#' @param perVariantWeights \code{FALSE} (default) weights every variant by
#'   n; \code{TRUE} weights each by qgg's 1 / (se^2 + b^2 / n_j), from its
#'   own sample size.
#' @param h2 The prior guess of the heritability that sets the starting
#'   variances. Default \code{0.5}, qgg's and hibayes'.
#' @param shrinkCorr The LD shrinkage \eqn{c}: the sampler reads
#'   \eqn{cR + (1 - c)I}. Default 1, none.
#' @param ldProjection The ridge \eqn{a} of SDPR's likelihood (Zhou and Zhao
#'   2021), under which the sampler reads the LD as
#'   \eqn{R (R + aI)^{-1} R} and the marginal effects as
#'   \eqn{R (R + aI)^{-1} \hat\beta}: each of R's eigen-directions keeps
#'   \eqn{\lambda / (\lambda + a)} of the effects, so their part outside
#'   the reference LD's span, which a reference smaller than the window or
#'   unlike the study leaves, never enters. Without it the samplers' effects
#'   can grow without bound along those directions. Default 0.1, SDPR's;
#'   \code{NULL} reads the LD and the effects as they are.
#' @param allowJumpSign Keep a drawn effect whose sign differs from the
#'   current one; \code{FALSE} sets it to zero, and out of the slab,
#'   instead. Default \code{TRUE}.
#' @param raoBlackwell Report each effect's conditional posterior mean,
#'   averaged over the kept iterations, rather than the draws' average.
#'   Default \code{FALSE}.
#' @param nChains The number of chains, each from a seed drawn from R's
#'   random number generator; they are tested for divergence and, with more
#'   than one, combined as \code{\link{prsCs}}'s are. Default 1.
#' @param verbose Whether to report the chain's progress.
#' @return A list: \code{betaEst} (the posterior mean effects), \code{pip}
#'   (for \code{bayesB}, \code{bayesC} and \code{bayesR}, each effect's
#'   posterior probability of being non-null), \code{pi} (the non-null
#'   proportion) and \code{residualVariance}.
#' @keywords internal
bayesAlphabetRss <- function(
    bhat,
    F,
    n,
    method,
    nIter = 5000,
    nBurnin = 1000,
    thin = 5,
    pi = NULL,
    nPerVariant = NULL,
    bayesLForm = NULL,
    perVariantWeights = FALSE,
    h2 = 0.5,
    shrinkCorr = 1,
    ldProjection = 0.1,
    allowJumpSign = TRUE,
    raoBlackwell = FALSE,
    nChains = 1,
    verbose = FALSE
) {
    method <- arg_match(method, names(.alphabetBackends$sumstat))
    .alphabetRssAssert(n, h2, ldProjection, allowJumpSign)
    .mcmcAssertShared(shrinkCorr, raoBlackwell, nChains)
    sampler <- list(
        method = method,
        legarra = .alphabetLegarra(method, bayesLForm),
        pi = pi %||% 0,
        h2 = h2,
        nIter = as.integer(nIter),
        nBurnin = as.integer(nBurnin),
        thin = as.integer(thin),
        verbose = verbose,
        shrinkCorr = shrinkCorr,
        allowJumpSign = allowJumpSign,
        raoBlackwell = raoBlackwell
    )
    data <- .alphabetData(
        bhat,
        F,
        n,
        nPerVariant,
        perVariantWeights,
        ldProjection
    )
    .alphabetRssChains(c(data, sampler), method, nChains, bhat)
}

# What bayesAlphabetRss() checks besides the settings every sampler shares.
# @noRd
.alphabetRssAssert <- function(n, h2, ldProjection, allowJumpSign) {
    .alphabetAssertSampleSize(n, "bayesAlphabetRss")
    assertNumber(h2, lower = 0, upper = 1)
    assertNumber(ldProjection, lower = 0, finite = TRUE, null.ok = TRUE)
    assertFlag(allowJumpSign)
    invisible(NULL)
}

# The data pecotmr's alphabet sampler reads: the marginal effects and the
# factor of the LD, through SDPR's likelihood unless `ldProjection` is NULL
# (.alphabetProjected()), the per-variant weights, taken from the marginal
# effects as given, and n.
# @noRd
.alphabetData <- function(
    bhat,
    F,
    n,
    nPerVariant,
    perVariantWeights,
    ldProjection
) {
    projected <- .alphabetProjected(bhat, F, ldProjection)
    list(
        bhat = projected$bhat,
        factorR = projected$F,
        weights = .alphabetRssWeightsFor(
            bhat,
            n,
            nPerVariant,
            perVariantWeights
        ),
        n = as.integer(n)
    )
}

# SDPR's likelihood (Zhou and Zhao 2021, its optLlk = 1) for the alphabet
# sampler: the RSS likelihood with the marginal effects' covariance
# (R + aI) / n, under which the sampler reads R (R + aI)^-1 R as the LD and
# R (R + aI)^-1 bhat as the effects. Each eigen-direction of R keeps
# lambda / (lambda + a) of the effects, so their part outside the reference
# LD's span -- left by a reference smaller than the window, or unlike the
# study -- never enters, nor does the runaway it drives through the variance
# updates. From the thin SVD of F (.ldFactorEigen()), over the r directions
# with a positive eigenvalue, the LD is the factor
# sqrt(r) diag(lambda / sqrt(lambda + a)) V' and a zero row, so that the
# sampler's normalization by its rows less one holds. NULL `a` leaves the
# effects and F as they are.
# @noRd
.alphabetProjected <- function(bhat, F, a) {
    if (is.null(a)) {
        return(list(bhat = bhat, F = F))
    }
    ev <- .ldFactorEigen(F)
    kept <- ev$values > sqrt(.Machine$double.eps) * max(ev$values)
    lambda <- ev$values[kept]
    V <- ev$vectors[, kept, drop = FALSE]
    list(
        bhat = as.numeric(V %*% (lambda / (lambda + a) * crossprod(V, bhat))),
        F = rbind(sqrt(length(lambda)) * (lambda / sqrt(lambda + a)) * t(V), 0)
    )
}

# The chains of pecotmr's alphabet sampler, combined (.mcmcChains()); its
# divergence warning tests the effects against the marginal effects as
# given. Only bayesB, bayesC and bayesR report inclusion probabilities.
# @noRd
.alphabetRssChains <- function(settings, method, nChains, bhat) {
    fit <- .mcmcChains(
        ".alphabetSamplerChain",
        list(settings = settings),
        as.list(sample.int(.Machine$integer.max, nChains)),
        c("betaEst", "pip", "pi", "residualVariance"),
        "bayesAlphabetRss",
        bhat
    )
    if (!is_in(method, c("bayesB", "bayesC", "bayesR"))) {
        fit$pip <- NULL
    }
    fit
}

# One chain of pecotmr's alphabet sampler.
# @noRd
.alphabetSamplerChain <- function(settings, seed, track) {
    exec(bayesAlphabetFactorCpp, !!!settings, seed = seed, track = track)
}

# Whether bayesL takes Legarra et al.'s form; `bayesLForm` is bayesL's
# alone.
# @noRd
.alphabetLegarra <- function(method, bayesLForm) {
    if (is.null(bayesLForm)) {
        return(FALSE)
    }
    if (!identical(method, "bayesL")) {
        abort(glue("`bayesLForm` applies to bayesL, not {method}."))
    }
    identical(arg_match0(bayesLForm, c("parkCasella", "legarra")), "legarra")
}

# The per-variant weights: n for every variant, or qgg's
# 1 / (se^2 + b^2 / n_j) from each one's own sample size.
# @noRd
.alphabetRssWeightsFor <- function(bhat, n, nPerVariant, perVariantWeights) {
    if (!isTRUE(perVariantWeights)) {
        return(rep(as.numeric(n), length(bhat)))
    }
    nj <- nPerVariant %||% rep(n, length(bhat))
    if (length(nj) != length(bhat)) {
        abort("`nPerVariant` must give one sample size per variant.")
    }
    seb2 <- (1 - bhat^2) / (nj - 2)
    1 / (seb2 + bhat^2 / nj)
}

#' @title BayesA TWAS weights from summary statistics
#' @description \code{\link{bayesAlphabetRssWeights}} with the BayesA prior.
#' @inheritParams bayesAlphabetRssWeights
#' @param backend \code{"qgg"}, \code{"hibayes"} or \code{"pecotmr"};
#'   \code{NULL} picks the default for the LD reference.
#' @return A numeric vector of the posterior mean effects, one per variant.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' stat <- list(
#'   z = as.vector(cor(X, eqtlRegionExample$yRes)) * sqrt(nrow(X)),
#'   n = rep(nrow(X), ncol(X)),
#'   variantNames = colnames(X)
#' )
#' bayesARssWeights(stat, LdMatrix(cor(X)), "hibayes", nIter = 1000,
#'   nBurnin = 200)
#' @export
bayesARssWeights <- function(
    stat,
    ldReference,
    backend = NULL,
    nIter = 5000,
    nBurnin = 1000,
    thin = 5,
    shrinkCorr = 0.95,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    bayesAlphabetRssWeights(
        stat,
        ldReference,
        "bayesA",
        backend = backend,
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        shrinkCorr = shrinkCorr,
        fitRetention = fitRetention,
        methodArgs = methodArgs
    )
}

#' @title BayesB TWAS weights from summary statistics
#' @description \code{\link{bayesAlphabetRssWeights}} with the BayesB prior.
#' @inheritParams bayesAlphabetRssWeights
#' @param backend \code{"hibayes"} or \code{"pecotmr"}; \code{NULL} picks
#'   the default for the LD reference.
#' @param pi The prior proportion of non-null effects, where each engine's
#'   estimate of it starts. Default \code{0.2}.
#' @return A numeric vector of the posterior mean effects, one per variant.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' stat <- list(
#'   z = as.vector(cor(X, eqtlRegionExample$yRes)) * sqrt(nrow(X)),
#'   n = rep(nrow(X), ncol(X)),
#'   variantNames = colnames(X)
#' )
#' bayesBRssWeights(stat, LdMatrix(cor(X)), nIter = 1000, nBurnin = 200)
#' @export
bayesBRssWeights <- function(
    stat,
    ldReference,
    backend = NULL,
    pi = 0.2,
    nIter = 5000,
    nBurnin = 1000,
    thin = 5,
    shrinkCorr = 0.95,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    bayesAlphabetRssWeights(
        stat,
        ldReference,
        "bayesB",
        backend = backend,
        pi = pi,
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        shrinkCorr = shrinkCorr,
        fitRetention = fitRetention,
        methodArgs = methodArgs
    )
}

#' @title BayesC TWAS weights from summary statistics
#' @description \code{\link{bayesAlphabetRssWeights}} with the BayesC prior.
#' @inheritParams bayesAlphabetRssWeights
#' @param backend \code{"qgg"}, \code{"hibayes"} or \code{"pecotmr"};
#'   \code{NULL} picks the default for the LD reference.
#' @param pi The prior proportion of non-null effects, where each engine's
#'   estimate of it starts. Default \code{0.1}.
#' @return A numeric vector of the posterior mean effects, one per variant.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' stat <- list(
#'   z = as.vector(cor(X, eqtlRegionExample$yRes)) * sqrt(nrow(X)),
#'   n = rep(nrow(X), ncol(X)),
#'   variantNames = colnames(X)
#' )
#' bayesCRssWeights(stat, LdMatrix(cor(X)), "hibayes", nIter = 1000,
#'   nBurnin = 200)
#' @export
bayesCRssWeights <- function(
    stat,
    ldReference,
    backend = NULL,
    pi = 0.1,
    nIter = 5000,
    nBurnin = 1000,
    thin = 5,
    shrinkCorr = 0.95,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    bayesAlphabetRssWeights(
        stat,
        ldReference,
        "bayesC",
        backend = backend,
        pi = pi,
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        shrinkCorr = shrinkCorr,
        fitRetention = fitRetention,
        methodArgs = methodArgs
    )
}

#' @title BayesL TWAS weights from summary statistics
#' @description \code{\link{bayesAlphabetRssWeights}} with the Bayesian
#'   LASSO.
#' @inheritParams bayesAlphabetRssWeights
#' @param backend \code{"hibayes"}, \code{"qgg"} or \code{"pecotmr"};
#'   \code{NULL} picks the default for the LD reference.
#' @return A numeric vector of the posterior mean effects, one per variant.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' stat <- list(
#'   z = as.vector(cor(X, eqtlRegionExample$yRes)) * sqrt(nrow(X)),
#'   n = rep(nrow(X), ncol(X)),
#'   variantNames = colnames(X)
#' )
#' bayesLRssWeights(stat, LdMatrix(cor(X)), nIter = 1000, nBurnin = 200)
#' @export
bayesLRssWeights <- function(
    stat,
    ldReference,
    backend = NULL,
    nIter = 5000,
    nBurnin = 1000,
    thin = 5,
    shrinkCorr = 0.95,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    bayesAlphabetRssWeights(
        stat,
        ldReference,
        "bayesL",
        backend = backend,
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        shrinkCorr = shrinkCorr,
        fitRetention = fitRetention,
        methodArgs = methodArgs
    )
}

#' @title BayesN TWAS weights from summary statistics
#' @description \code{\link{bayesAlphabetRssWeights}} with a normal prior
#'   shared by every effect.
#' @inheritParams bayesAlphabetRssWeights
#' @param backend \code{"qgg"}, \code{"hibayes"} or \code{"pecotmr"};
#'   \code{NULL} picks the default for the LD reference.
#' @return A numeric vector of the posterior mean effects, one per variant.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' stat <- list(
#'   z = as.vector(cor(X, eqtlRegionExample$yRes)) * sqrt(nrow(X)),
#'   n = rep(nrow(X), ncol(X)),
#'   variantNames = colnames(X)
#' )
#' bayesNRssWeights(stat, LdMatrix(cor(X)), "hibayes", nIter = 1000,
#'   nBurnin = 200)
#' @export
bayesNRssWeights <- function(
    stat,
    ldReference,
    backend = NULL,
    nIter = 5000,
    nBurnin = 1000,
    thin = 5,
    shrinkCorr = 0.95,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    bayesAlphabetRssWeights(
        stat,
        ldReference,
        "bayesN",
        backend = backend,
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        shrinkCorr = shrinkCorr,
        fitRetention = fitRetention,
        methodArgs = methodArgs
    )
}

#' @title BayesR TWAS weights from summary statistics
#' @description \code{\link{bayesAlphabetRssWeights}} with the BayesR prior.
#' @inheritParams bayesAlphabetRssWeights
#' @param backend \code{"qgg"}, \code{"hibayes"} or \code{"pecotmr"};
#'   \code{NULL} picks the default for the LD reference.
#' @return A numeric vector of the posterior mean effects, one per variant.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' stat <- list(
#'   z = as.vector(cor(X, eqtlRegionExample$yRes)) * sqrt(nrow(X)),
#'   n = rep(nrow(X), ncol(X)),
#'   variantNames = colnames(X)
#' )
#' bayesRRssWeights(stat, LdMatrix(cor(X)), "hibayes", nIter = 1000,
#'   nBurnin = 200)
#' @export
bayesRRssWeights <- function(
    stat,
    ldReference,
    backend = NULL,
    nIter = 5000,
    nBurnin = 1000,
    thin = 5,
    shrinkCorr = 0.95,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    bayesAlphabetRssWeights(
        stat,
        ldReference,
        "bayesR",
        backend = backend,
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        shrinkCorr = shrinkCorr,
        fitRetention = fitRetention,
        methodArgs = methodArgs
    )
}

# --- LDpred2 -----------------------------------------------------------------
#
# LDpred2-auto and LDpred2-grid run on bigsnpr, on the LD matrix written to a
# temporary file, or on pecotmr's port of its samplers (ldpred2FactorCpp()),
# on a genotype reference's factor of the LD. Both take the statistics in
# standardized units, z / sqrt(n_j) for each variant's own sample size, and
# return weights in the same units.

#' LDpred2-auto Weights from Summary Statistics
#'
#' LDpred2-auto (Prive et al. 2020, 2023) estimates the point-normal prior's
#' proportion of causal variants \eqn{p} and the heritability \eqn{h^2}
#' within its Gibbs sampler, so nothing is tuned. One chain runs from each
#' starting \eqn{p} in \code{pInit}. As in the LDpred2 vignette, a chain
#' whose fitted correlations span less than 0.95 times the 95th percentile of
#' the chains' spans is dropped, and the rest are averaged. A chain has
#' diverged, and is dropped too, once a kept draw's squared effects sum past
#' twice the squared marginal effects; bigsnpr's own sampler tests every
#' iteration, burn-in included.
#'
#' The defaults are the LDpred2 vignette's, without the MAF-dependent prior
#' variance.
#'
#' @details
#' \code{backend = "bigsnpr"} runs bigsnpr's \code{snp_ldpred2_auto()} on the
#' LD matrix, written to a temporary file for the call; for a genotype
#' reference the matrix is formed from its genotypes. \code{"pecotmr"} runs
#' pecotmr's port of the same sampler on a genotype reference's factor of the
#' LD, never forming the matrix. bigsnpr's functions call
#' \code{foreach::registerDoSEQ()}, which replaces any foreach parallel
#' backend registered in the session. Both draw their seeds from R's random
#' number generator, so \code{set.seed()} makes the weights reproducible.
#'
#' @param stat A list of summary statistics: \code{n} (per-variant sample
#'   sizes) and the marginal effects as \code{cor} (standardized), \code{z},
#'   \code{b} with \code{seb}, or \code{b} alone, taken as already
#'   standardized. The sampler reads z / sqrt(n) for each variant's own
#'   sample size.
#' @param ldReference The LD reference the statistics' variants are matched
#'   to by \code{stat$variantNames}: a genotype panel (see
#'   \code{\link{readGenotypes}}), a \code{\link{GenotypeMixture}} or an
#'   \code{\link{LdMatrix}}.
#' @param backend \code{"bigsnpr"} or \code{"pecotmr"}. \code{NULL}
#'   (default) is pecotmr for a genotype reference and bigsnpr for an
#'   \code{\link{LdMatrix}}; pecotmr's sampler reads only a genotype
#'   reference.
#' @param h2Init The starting heritability. \code{NULL} (default) estimates
#'   it for the window by LAVA, as \code{\link{estimateH2}} does for QTL
#'   summary statistics, which reads the reference's sample size, so an
#'   \code{LdMatrix} must record \code{nSamples}. Floored at 0.001, which is
#'   also the start where LAVA gives no estimate.
#' @param pInit The starting proportions of causal variants, one chain each.
#' @param nIter,nBurnin Each chain's iterations, burn-in included, and its
#'   burn-in.
#' @param shrinkCorr The LD shrinkage \eqn{c}: the sampler reads
#'   \eqn{cR + (1 - c)I}.
#' @param allowJumpSign Keep a drawn effect whose sign differs from the
#'   current one; \code{FALSE} sets it to zero instead.
#' @param pBounds The bounds on \eqn{p}.
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights. \code{"none"} (default) attaches
#'   nothing. \code{"slim"} keeps \code{chains}, a table with one row per
#'   chain (its starting \code{pInit}, its estimates \code{p} and \code{h2},
#'   whether it \code{diverged} and whether it was \code{kept}), and the kept
#'   chains' mean inclusion probabilities \code{postp}, \code{p} and
#'   \code{h2}. \code{"full"} also keeps \code{beta}, each chain's effects
#'   (variants x chains).
#' @return A numeric vector of weights, the kept chains' mean effects. Unless
#'   \code{fitRetention = "none"}, the fit is attached as attribute
#'   \code{"fit"}.
#' @examples
#' if (requireNamespace("bigsnpr", quietly = TRUE)) {
#'   data(eqtlRegionExample)
#'   X <- eqtlRegionExample$X[, 1:30]
#'   y <- eqtlRegionExample$yRes
#'   stat <- list(
#'     cor = as.vector(cor(X, y)),
#'     n = rep(nrow(X), ncol(X)),
#'     variantNames = colnames(X)
#'   )
#'   ld <- LdMatrix(cor(X), nSamples = nrow(X))
#'   set.seed(1)
#'   w <- ldpred2AutoWeights(stat, ld, h2Init = 0.1, pInit = c(0.01, 0.1),
#'     fitRetention = "slim")
#'   attr(w, "fit")$chains
#' }
#' @export
ldpred2AutoWeights <- function(
    stat,
    ldReference,
    backend = NULL,
    h2Init = NULL,
    pInit = exp(seq(log(1e-4), log(0.2), length.out = 30)),
    nIter = 700,
    nBurnin = 500,
    shrinkCorr = 0.95,
    allowJumpSign = FALSE,
    pBounds = c(1e-5, 1),
    fitRetention = c("none", "slim", "full")
) {
    fitRetention <- arg_match(fitRetention)
    label <- "ldpred2AutoWeights"
    backend <- .rssBackend("ldpred2Auto", backend, ldReference)
    .assertMcmcSchedule(nIter, nBurnin, 1L)
    .ldpred2AssertProportions(pInit, "pInit")
    .ldpred2AssertProportions(pBounds, "pBounds")
    .ldpred2AssertSampler(shrinkCorr, allowJumpSign)
    inp <- .rssPerVariantInput(stat, ldReference, label)
    config <- list(
        pInit = pInit,
        h2 = .ldpred2H2(h2Init, inp, label),
        nIter = nIter,
        nBurnin = nBurnin,
        shrinkCorr = shrinkCorr,
        allowJumpSign = allowJumpSign,
        pBounds = sort(pBounds)
    )
    engine <- .ldpred2Engine(inp, backend)
    on.exit(.bigsnprRelease(engine), add = TRUE)
    .withWeightsFitIds(
        .ldpred2AutoWeightsFrom(
            .ldpred2RunAuto(inp, engine, config),
            fitRetention,
            label
        ),
        stat[["variantNames"]]
    )
}

#' LDpred2-grid Weights from Summary Statistics
#'
#' LDpred2-grid (Prive et al. 2020) runs the point-normal prior's Gibbs
#' sampler at every combination of the proportion of causal variants
#' \code{p}, the heritability \code{h2} and \code{sparse}, and returns one of
#' the fits. Upstream picks it on a validation sample; here, as for the
#' penalized summary-statistics methods, it is picked on held-out summary
#' statistics by default (see Details of \code{\link{lassosumRssWeights}}). A
#' grid point has diverged, and is not picked, once a kept draw's squared
#' effects sum past twice the squared marginal effects; bigsnpr's own
#' sampler tests every iteration, burn-in included.
#'
#' @details
#' \code{backend = "bigsnpr"} runs bigsnpr's \code{snp_ldpred2_grid()} on the
#' LD matrix, written to a temporary file for the call; \code{"pecotmr"} runs
#' pecotmr's port of the same sampler on a genotype reference's factor of the
#' LD. bigsnpr's grid sampler has no LD shrinkage of its own, so for it the
#' matrix written is shrunk, \eqn{cR + (1 - c)I}, as the Bayesian alphabet's
#' matrix engines are given theirs; it always allows a drawn effect to
#' change sign, so it takes no \code{allowJumpSign}. pecotmr's sampler takes
#' both, as LDpred2-auto's does. See \code{\link{ldpred2AutoWeights}} for the
#' temporary file, the foreach backend and the seeds.
#'
#' @inheritParams ldpred2AutoWeights
#' @param p The proportions of causal variants in the grid.
#' @param h2 The heritabilities in the grid. \code{NULL} (default) is the
#'   window's LAVA estimate (see \code{h2Init} of
#'   \code{\link{ldpred2AutoWeights}}) times 0.3, 0.7, 1 and 1.4, as the
#'   LDpred2 vignette scales LD score regression's.
#' @param sparse Whether each fit sets an effect whose inclusion probability
#'   falls below \code{p} to zero; both by default.
#' @param shrinkCorr The LD shrinkage \eqn{c}: both engines read
#'   \eqn{cR + (1 - c)I}. Default 0.95, as for
#'   \code{\link{ldpred2AutoWeights}}; bigsnpr's grid itself has none.
#' @param allowJumpSign As for \code{\link{ldpred2AutoWeights}}; bigsnpr's
#'   grid sampler takes only \code{TRUE}, the default.
#' @param selection How the returned fit is chosen: \code{"heldOut"}
#'   (default) or \code{"ldQuadratic"}, as for
#'   \code{\link{lassosumRssWeights}}. \code{"heldOut"} needs the
#'   reference's sample size and refits the grid four times.
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights. \code{"none"} (default) attaches
#'   nothing. \code{"slim"} keeps \code{path}, the grid (\code{p},
#'   \code{h2}, \code{sparse}) with whether each point \code{diverged} and
#'   which was \code{selected}. \code{"full"} also keeps \code{beta}, every
#'   point's effects (variants x points).
#' @return A numeric vector of weights, the selected fit's effects, with the
#'   selection recorded in attribute \code{"ldpred2_selection"}. Unless
#'   \code{fitRetention = "none"}, the grid is attached as attribute
#'   \code{"fit"}.
#' @examples
#' if (requireNamespace("bigsnpr", quietly = TRUE)) {
#'   data(eqtlRegionExample)
#'   X <- eqtlRegionExample$X[, 1:30]
#'   y <- eqtlRegionExample$yRes
#'   stat <- list(
#'     cor = as.vector(cor(X, y)),
#'     n = rep(nrow(X), ncol(X)),
#'     variantNames = colnames(X)
#'   )
#'   ld <- LdMatrix(cor(X), nSamples = nrow(X))
#'   set.seed(1)
#'   w <- ldpred2GridWeights(stat, ld, p = c(0.01, 0.1, 1), h2 = 0.1)
#'   attr(w, "ldpred2_selection")
#' }
#' @export
ldpred2GridWeights <- function(
    stat,
    ldReference,
    backend = NULL,
    p = exp(seq(log(1e-5), 0, length.out = 21)),
    h2 = NULL,
    sparse = c(FALSE, TRUE),
    nIter = 150,
    nBurnin = 50,
    shrinkCorr = 0.95,
    allowJumpSign = TRUE,
    selection = c("heldOut", "ldQuadratic"),
    fitRetention = c("none", "slim", "full")
) {
    selection <- arg_match(selection)
    fitRetention <- arg_match(fitRetention)
    label <- "ldpred2GridWeights"
    backend <- .rssBackend("ldpred2Grid", backend, ldReference)
    .assertMcmcSchedule(nIter, nBurnin, 1L)
    .ldpred2AssertProportions(p, "p")
    .ldpred2AssertSampler(shrinkCorr, allowJumpSign)
    .ldpred2CheckGridBackend(backend, allowJumpSign, label)
    inp <- .rssPerVariantInput(stat, ldReference, label)
    if (identical(selection, "heldOut")) {
        .rssCheckHeldOutLd(inp$ld, label)
    }
    config <- list(
        grid = .ldpred2Grid(p, .ldpred2GridH2(h2, inp, label), sparse),
        nIter = nIter,
        nBurnin = nBurnin,
        shrinkCorr = shrinkCorr,
        allowJumpSign = allowJumpSign
    )
    engine <- .ldpred2Engine(inp, backend, shrinkCorr)
    on.exit(.bigsnprRelease(engine), add = TRUE)
    .withWeightsFitIds(
        .ldpred2GridWeightsFrom(
            .ldpred2RunGrid(inp, engine, config),
            list(inp = inp, engine = engine, config = config),
            selection,
            fitRetention,
            label
        ),
        stat[["variantNames"]]
    )
}

# Proportions strictly between 0 and 1, the upper end included.
# @noRd
.ldpred2AssertProportions <- function(x, arg) {
    assertNumeric(x, lower = 0, upper = 1, any.missing = FALSE, min.len = 1L)
    if (any(x <= 0)) {
        abort(glue("`{arg}` must be above 0."))
    }
    invisible(NULL)
}

# @noRd
.ldpred2AssertSampler <- function(shrinkCorr, allowJumpSign) {
    .assertShrinkCorr(shrinkCorr)
    assertFlag(allowJumpSign)
    invisible(NULL)
}

# bigsnpr's grid sampler has neither the LD shrinkage nor the sign rule.
# @noRd
.ldpred2CheckGridBackend <- function(backend, allowJumpSign, label) {
    if (identical(backend, "bigsnpr") && !allowJumpSign) {
        abort(c(
            glue(
                "{label}: bigsnpr's grid sampler always allows a drawn ",
                "effect to change sign, so it takes no `allowJumpSign`."
            ),
            i = "Use backend = \"pecotmr\" with a genotype reference."
        ))
    }
    invisible(NULL)
}

# The statistics bigsnpr's engines and pecotmr's LDpred2 read: the LD source
# (.rssLdSource()), the standardized effects z / sqrt(n_j) and each
# variant's sample size.
# @noRd
.rssPerVariantInput <- function(stat, ldReference, label) {
    ld <- .rssLdSource(stat, ldReference, label)
    p <- .ldSourceSize(ld)
    n <- rep_len(as.numeric(stat$n), p)
    bhat <- .lassosumClampCor(.rssMarginalCor(stat, n = n, p = p))
    .rssValidateSource(bhat, ld, median(n))
    list(ld = ld, bhat = bhat, n = n)
}

# The heritability LDpred2-auto starts from: the caller's, floored at 0.001
# as bigsnpr floors it, else the window's LAVA estimate.
# @noRd
.ldpred2H2 <- function(h2, inp, label) {
    if (is.null(h2)) {
        return(.ldpred2LavaH2(inp, label))
    }
    assertNumber(h2, lower = 0, upper = 1)
    max(h2, 1e-3)
}

# The heritabilities of the grid: the caller's, else the window's LAVA
# estimate times 0.3, 0.7, 1 and 1.4.
# @noRd
.ldpred2GridH2 <- function(h2, inp, label) {
    if (!is.null(h2)) {
        .ldpred2AssertProportions(h2, "h2")
        return(h2)
    }
    .ldpred2LavaH2(inp, label) * c(0.3, 0.7, 1, 1.4)
}

# The window's heritability by LAVA (.cisH2Lava(), estimateH2()'s default
# for QTL summary statistics), which reads the reference's sample size,
# floored at bigsnpr's 0.001: the start too where LAVA gives no estimate.
# @noRd
.ldpred2LavaH2 <- function(inp, label) {
    if (is.na(inp$ld$nSamples)) {
        abort(c(
            glue(
                "{label}: estimating the heritability by LAVA needs the ",
                "LD reference's sample size."
            ),
            i = "Build the LdMatrix with `nSamples`, or give the heritability."
        ))
    }
    h2 <- .cisH2Lava(inp$bhat * sqrt(inp$n), inp$n, inp$ld)$h2
    if (is.na(h2)) 1e-3 else max(h2, 1e-3)
}

# Every combination of p, h2 and sparse, large p first as bigsnpr runs them.
# @noRd
.ldpred2Grid <- function(p, h2, sparse) {
    assertLogical(sparse, any.missing = FALSE, min.len = 1L)
    grid <- expand.grid(p = p, h2 = h2, sparse = sparse)
    as_tibble(grid[order(-grid$p, grid$sparse, -grid$h2), ])
}

# What an LDpred2 run reads: pecotmr's sampler reads the factor in `inp`,
# bigsnpr its sparse matrix (.bigsnprEngine()), shrunk to `shrinkCorr` for a
# sampler of bigsnpr's with no shrinkage of its own (the grid's).
# @noRd
.ldpred2Engine <- function(inp, backend, shrinkCorr = 1) {
    if (identical(backend, "pecotmr")) {
        return(list(backend = backend))
    }
    .bigsnprEngine(inp, shrinkCorr)
}

# bigsnpr's LD: the matrix as bigsnpr's sparse matrix, written to a temporary
# file that .bigsnprRelease() deletes; for a genotype reference, formed from
# its factor; shrunk to `shrinkCorr` (.rssShrunkMatrix()).
# @noRd
.bigsnprEngine <- function(inp, shrinkCorr = 1) {
    if (!requireNamespace("bigsnpr", quietly = TRUE)) {
        abort("Package 'bigsnpr' is required for backend = \"bigsnpr\".")
    }
    R <- .rssShrunkMatrix(
        inp$ld$R %||% (crossprod(inp$ld$F) / (nrow(inp$ld$F) - 1)),
        shrinkCorr
    )
    corr <- bigsnpr::as_SFBM(
        methods::as(Matrix::Matrix(R, sparse = TRUE), "generalMatrix"),
        tempfile(),
        compact = TRUE
    )
    list(backend = "bigsnpr", corr = corr)
}

# The LD matrix shrunk toward the identity, c R + (1 - c) I, for an engine
# that reads the matrix and has no shrinkage of its own (qgg's and hibayes'
# samplers, bigsnpr's grid).
# @noRd
.rssShrunkMatrix <- function(R, shrinkCorr) {
    if (shrinkCorr == 1) {
        return(R)
    }
    R <- shrinkCorr * R
    diag(R) <- diag(R) + (1 - shrinkCorr)
    R
}

# @noRd
.bigsnprRelease <- function(engine) {
    if (!is.null(engine$corr)) {
        unlink(engine$corr$sbk)
    }
    invisible(NULL)
}

# bigsnpr's statistics: its rescaling sqrt(n se^2 + b^2) is then 1, so it
# reads the standardized effects as they are.
# @noRd
.bigsnprDfBeta <- function(inp) {
    data.frame(
        beta = inp$bhat,
        beta_se = sqrt((1 - inp$bhat^2) / inp$n),
        n_eff = inp$n
    )
}

# LDpred2-auto's chains, each a list of betaEst, corrEst, postp, h2, p,
# pInit and diverged. bigsnpr seeds from R's generator, so it runs under a
# seed drawn from it, leaving the caller's stream one draw on.
# @noRd
.ldpred2RunAuto <- function(inp, engine, config) {
    if (identical(engine$backend, "bigsnpr")) {
        return(withr::with_seed(
            .sessionSeedDraw(),
            .ldpred2BigsnprAuto(inp, engine$corr, config)
        ))
    }
    meanLd <- mean(.ldSquaredProduct(inp$ld, matrix(1, length(inp$bhat), 1L)))
    map(
        config$pInit,
        .ldpred2PecotmrAutoChain,
        inp = inp,
        config = config,
        meanLd = meanLd
    )
}

# @noRd
.ldpred2BigsnprAuto <- function(inp, corr, config) {
    runs <- bigsnpr::snp_ldpred2_auto(
        corr,
        .bigsnprDfBeta(inp),
        h2_init = config$h2,
        vec_p_init = config$pInit,
        burn_in = config$nBurnin,
        num_iter = config$nIter - config$nBurnin,
        allow_jump_sign = config$allowJumpSign,
        shrink_corr = config$shrinkCorr,
        use_MLE = FALSE,
        p_bounds = config$pBounds,
        ncores = 1
    )
    map(runs, .ldpred2BigsnprChain)
}

# @noRd
.ldpred2BigsnprChain <- function(run) {
    list(
        betaEst = run$beta_est,
        corrEst = run$corr_est,
        postp = run$postp_est,
        h2 = run$h2_est,
        p = run$p_est,
        pInit = run$p_init,
        diverged = anyNA(run$beta_est)
    )
}

# One pecotmr chain from `pInit`; p and h2 are their means over the kept
# sweeps, as bigsnpr reports them.
# @noRd
#' @importFrom utils tail
.ldpred2PecotmrAutoChain <- function(pInit, inp, config, meanLd) {
    run <- ldpred2FactorCpp(
        inp$bhat,
        inp$ld$F,
        inp$n,
        TRUE,
        pInit,
        config$h2,
        FALSE,
        config$shrinkCorr,
        config$allowJumpSign,
        config$pBounds[[1L]],
        config$pBounds[[2L]],
        meanLd,
        as.integer(config$nIter),
        as.integer(config$nBurnin),
        FALSE,
        .sessionSeedDraw()
    )
    kept <- config$nIter - config$nBurnin
    missing <- rep(NA_real_, length(inp$bhat))
    list(
        betaEst = if (run$diverged) missing else run$betaEst,
        corrEst = if (run$diverged) missing else run$corrEst,
        postp = if (run$diverged) missing else run$postp,
        h2 = if (run$diverged) NA_real_ else mean(tail(run$pathH2, kept)),
        p = if (run$diverged) NA_real_ else mean(tail(run$pathP, kept)),
        pInit = pInit,
        diverged = run$diverged
    )
}

# The grid's fits, variants x points, NA for a point that diverged.
# @noRd
.ldpred2RunGrid <- function(inp, engine, config) {
    if (identical(engine$backend, "bigsnpr")) {
        return(withr::with_seed(
            .sessionSeedDraw(),
            bigsnpr::snp_ldpred2_grid(
                engine$corr,
                .bigsnprDfBeta(inp),
                as.data.frame(config$grid),
                burn_in = config$nBurnin,
                num_iter = config$nIter - config$nBurnin,
                ncores = 1
            )
        ))
    }
    fits <- pmap(
        config$grid,
        .ldpred2PecotmrGridPoint,
        inp = inp,
        config = config
    )
    exec(cbind, !!!fits)
}

# @noRd
.ldpred2PecotmrGridPoint <- function(p, h2, sparse, inp, config) {
    run <- ldpred2FactorCpp(
        inp$bhat,
        inp$ld$F,
        inp$n,
        FALSE,
        p,
        h2,
        sparse,
        config$shrinkCorr,
        config$allowJumpSign,
        1e-5,
        1,
        1,
        as.integer(config$nIter),
        as.integer(config$nBurnin),
        FALSE,
        .sessionSeedDraw()
    )
    if (run$diverged) rep(NA_real_, length(inp$bhat)) else run$betaEst
}

# The chains LDpred2's vignette keeps: those whose fitted correlations span
# more than 0.95 times the 95th percentile of the spans, among the chains
# that did not diverge. Every sampler with several chains is combined so.
# @noRd
#' @importFrom stats quantile
.mcmcKeptChains <- function(chains, label) {
    span <- map_dbl(chains, .mcmcCorrSpan)
    kept <- !is.na(span) & span > 0.95 * quantile(span, 0.95, na.rm = TRUE)
    if (!any(kept)) {
        abort(glue("{label}: every chain diverged."))
    }
    kept
}

# @noRd
.mcmcCorrSpan <- function(chain) {
    if (chain$diverged) NA_real_ else diff(range(chain$corrEst))
}

# A sampler's estimates `fields` from one chain per seed in `seeds`. `run`
# names the function that runs one chain from `args`, a seed and whether to
# track it. One chain runs untracked and is returned as it is, with a
# warning when its effects outgrew the marginal effects `bhat`; several run
# tracked, so a diverging one stops, and are combined by .mcmcKeptChains(),
# their estimates averaged over the kept chains, with each chain's seed and
# fate in `chains`.
# @noRd
.mcmcChains <- function(run, args, seeds, fields, label, bhat) {
    if (length(seeds) == 1L) {
        one <- exec(run, !!!args, seed = seeds[[1L]], track = FALSE)
        .mcmcWarnDiverged(one$betaEst, bhat, label)
        return(one[fields])
    }
    chains <- map(seeds, .mcmcTrackedChain, run = run, args = args)
    kept <- .mcmcKeptChains(chains, label)
    c(
        map(set_names(fields), .mcmcMeanField, chains = chains[kept]),
        list(
            chains = tibble(
                seed = map_dbl(seeds, .mcmcSeedValue),
                diverged = map_lgl(chains, "diverged"),
                kept = kept
            )
        )
    )
}

# A single chain's posterior mean effects `beta` tested as a tracked chain's
# draws are: their squares summing past twice the marginal effects'. The
# draws themselves are not tested, since one noisy draw passes that bound
# on some null windows while the chain's mean stays near zero.
# @noRd
.mcmcWarnDiverged <- function(beta, bhat, label) {
    if (sum(beta^2) > 2 * sum(bhat^2)) {
        warn(glue(
            "{label}: the posterior mean effects' squares sum past twice ",
            "the squared marginal effects, so the chain diverged and its ",
            "effects are unreliable; several chains (`nChains`) drop the ",
            "chains that diverge."
        ))
    }
    invisible(NULL)
}

# @noRd
.mcmcTrackedChain <- function(seed, run, args) {
    exec(run, !!!args, seed = seed, track = TRUE)
}

# @noRd
.mcmcMeanField <- function(field, chains) {
    reduce(map(chains, field), `+`) / length(chains)
}

# @noRd
.mcmcSeedValue <- function(seed) {
    if (is.null(seed)) NA_real_ else as.numeric(seed)
}

# The chains' seeds: `seed`, `seed + 1`, ..., or NULL for each, which the
# sampler draws itself.
# @noRd
.mcmcSeeds <- function(seed, nChains) {
    if (is.null(seed)) {
        return(rep(list(NULL), nChains))
    }
    as.list(as.integer(seed) + seq_len(nChains) - 1L)
}

# The options every sampler shares.
# @noRd
.mcmcAssertShared <- function(shrinkCorr, raoBlackwell, nChains) {
    .assertShrinkCorr(shrinkCorr)
    assertFlag(raoBlackwell)
    assertCount(nChains, positive = TRUE)
    invisible(NULL)
}

# The LD shrinkage c of c R + (1 - c) I: above 0, at most 1.
# @noRd
.assertShrinkCorr <- function(shrinkCorr) {
    assertNumber(shrinkCorr, lower = 0, upper = 1)
    if (shrinkCorr <= 0) {
        abort("`shrinkCorr` must be above 0.")
    }
    invisible(NULL)
}

# The kept chains' mean effects, with the fit at the requested retention.
# @noRd
.ldpred2AutoWeightsFrom <- function(chains, fitRetention, label) {
    kept <- .mcmcKeptChains(chains, label)
    beta <- rowMeans(exec(cbind, !!!map(chains[kept], "betaEst")))
    if (identical(fitRetention, "none")) {
        return(beta)
    }
    fit <- list(
        chains = tibble(
            pInit = map_dbl(chains, "pInit"),
            p = map_dbl(chains, "p"),
            h2 = map_dbl(chains, "h2"),
            diverged = map_lgl(chains, "diverged"),
            kept = kept
        ),
        postp = rowMeans(exec(cbind, !!!map(chains[kept], "postp"))),
        p = mean(map_dbl(chains[kept], "p")),
        h2 = mean(map_dbl(chains[kept], "h2"))
    )
    if (identical(fitRetention, "full")) {
        fit$beta <- exec(cbind, !!!map(chains, "betaEst"))
    }
    `attr<-`(beta, "fit", fit)
}

# The grid's selected fit, by `selection` over the points that did not
# diverge, with the fit at the requested retention. `run` is the input, the
# engine and the configuration the grid ran with.
# @noRd
.ldpred2GridWeightsFrom <- function(
    candidates,
    run,
    selection,
    fitRetention,
    label
) {
    diverged <- colSums(is.na(candidates)) > 0
    if (all(diverged)) {
        abort(glue("{label}: every grid point diverged."))
    }
    usable <- replace(candidates, is.na(candidates), 0)
    sel <- .ldpred2SelectGrid(usable, run, selection)
    beta <- .ldpred2SelectionRecord(sel, run$config$grid)
    if (identical(fitRetention, "none")) {
        return(beta)
    }
    fit <- list(
        path = mutate(
            run$config$grid,
            diverged = diverged,
            selected = row_number() == sel$index
        )
    )
    if (identical(fitRetention, "full")) {
        fit$beta <- candidates
    }
    `attr<-`(beta, "fit", fit)
}

# @noRd
.ldpred2SelectGrid <- function(candidates, run, selection) {
    if (identical(selection, "ldQuadratic")) {
        return(.lassosumSelectLdQuadratic(candidates, run$inp$bhat, run$inp$ld))
    }
    .rssSelectHeldOut(
        candidates,
        run$inp$bhat,
        run$inp$ld,
        median(run$inp$n),
        list(fn = ".ldpred2RefitGrid", args = run)
    )
}

# The grid refit on a held-out split's training part; a point that diverges
# there is zero, so it scores -Inf.
# @noRd
.ldpred2RefitGrid <- function(trainCor, share, inp, engine, config) {
    train <- list_assign(inp, bhat = trainCor, n = share * inp$n)
    fits <- .ldpred2RunGrid(train, engine, config)
    replace(fits, is.na(fits), 0)
}

# The selected fit's effects, with its grid point recorded.
# @noRd
.ldpred2SelectionRecord <- function(sel, grid) {
    point <- grid[sel$index, ]
    `attr<-`(
        as.numeric(sel$beta),
        "ldpred2_selection",
        c(mode = sel$mode, p = point$p, h2 = point$h2, sparse = point$sparse)
    )
}

#' Lassosum RSS: LASSO on summary statistics with LD reference
#'
#' Coordinate descent to solve the penalized regression on summary statistics:
#' \deqn{f(\beta) = \beta' R \beta - 2\beta' r + 2\lambda ||\beta||_1}
#' where \eqn{R} is the LD matrix (pre-shrunk if desired) and \eqn{r =
#' \hat\beta / \sqrt{n}}.
#'
#' Based on Mak et al (2017) "Polygenic scores via penalized regression on
#' summary statistics", Genetic Epidemiology 41(6):469-480. It runs on the
#' solver of \code{\link{penalizedRss}}, whose LASSO at \code{alpha = 1} is
#' lassosum's coordinate descent: the same updates in the same order, so the
#' same coefficients bit for bit.
#'
#' Internal: \code{\link{lassosumRssWeights}} is the public interface. It
#' runs this over a grid of LD shrinkage values, forwards \code{lambda},
#' \code{thr} and \code{maxiter} through \code{\link{LassosumOptions}}, and
#' with \code{fitRetention} returns the paths as the \code{"fit"} attribute
#' of its weights.
#'
#' @param bhat A vector of marginal effect sizes.
#' @param R The LD correlation matrix (a single matrix over the analysed
#'   window), as in \code{susieR::susie_rss()}. If shrinkage is desired, apply
#'   it before passing (e.g., \code{(1-s)*R + s*I}).
#' @param n Sample size of the GWAS.
#' @param lambda A vector of L1 penalty values. Default: 20 values from 1e-4 to
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
#' @keywords internal
#' @importFrom checkmate assertCount assertNumber assertNumeric
lassosumRss <- function(
    bhat,
    R,
    n,
    lambda = exp(seq(log(0.0001), log(0.1), length.out = 20)),
    thr = 1e-4,
    maxiter = 10000
) {
    .penalizedRssPath(
        bhat,
        list(R = R),
        n,
        spec = .rssPenalty("lasso"),
        lambda = lambda,
        thr = thr,
        maxiter = maxiter
    )
}

# The values `fn` declares as the defaults of its formals `fmlNames`,
# evaluated as the function itself would evaluate them. The internal LD-source
# paths run lassosumRss(), penalizedRss(), prsCs() and sdpr() on a factor of
# the LD, and take those functions' defaults from here rather than restating
# them.
# @noRd
#' @importFrom rlang fn_fmls fn_env eval_bare
.formalDefaults <- function(fn, fmlNames) {
    map(fn_fmls(fn)[fmlNames], eval_bare, env = fn_env(fn))
}

# The standardized marginal effects, z / sqrt(n): the scale lassosum, PRS-CS
# and SDPR all fit. Read from `cor` as given, else from `z`, else from
# `b / seb`, which is a z-score whatever scale `b` is on, else from `b` alone,
# taken as already standardized.
.rssMarginalCor <- function(stat, n, p) {
    corInput <- if (!is.null(stat$cor)) {
        as.numeric(stat$cor)
    } else if (!is.null(stat$z)) {
        as.numeric(stat$z) / sqrt(n)
    } else if (!is.null(stat$b) && !is.null(stat$seb)) {
        as.numeric(stat$b) / as.numeric(stat$seb) / sqrt(n)
    } else if (!is.null(stat$b)) {
        as.numeric(stat$b)
    } else {
        abort("stat must contain one of 'cor', 'z', or 'b'.")
    }
    if (length(corInput) != p) {
        nInput <- length(corInput)
        msg <- glue(
            "The length of the input statistics ({nInput}) must equal the ",
            "number of LD variants ({p})."
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

# The first candidate whose score ties with the best (.rssScoreTies()).
# @noRd
.lassosumFirstMax <- function(x) {
    .rssScoreTies(x)[1L]
}

# The candidates whose scores tie with the best to rounding, within a
# relative sqrt(eps): candidates that differ only in scale score alike under
# a correlation, but not to the last bit. Where no score is finite, those
# equal to the best.
# @noRd
.rssScoreTies <- function(score) {
    best <- max(score, na.rm = TRUE)
    if (!is.finite(best)) {
        return(which(score == best))
    }
    which(score >= best - sqrt(.Machine$double.eps) * abs(best))
}

.lassosumSelectMinFbeta <- function(candidateBeta, candidateMeta) {
    idx <- which.min(candidateMeta$fbeta)
    list(
        beta = candidateBeta[, idx],
        index = idx,
        mode = "minFbeta"
    )
}

# lassosum's pseudovalidation: each candidate's in-sample score against the
# LD, with every marginal correlation first shrunk by its local false
# discovery rate, as upstream lassosum does.
# @noRd
.lassosumSelectLdQuadratic <- function(candidateBeta, corInput, ld) {
    shrunkCor <- corInput * (1 - .lassosumLfdr(corInput))
    idx <- .lassosumFirstMax(.rssShrunkLdScore(candidateBeta, shrunkCor, ld, 0))
    list(
        beta = candidateBeta[, idx],
        index = idx,
        mode = "ldQuadratic"
    )
}

# Each marginal correlation's local false discovery rate, from fdrtool's
# correlation null, as upstream lassosum's pseudovalidate() computes it.
# @noRd
.lassosumLfdr <- function(corInput) {
    if (!requireNamespace("fdrtool", quietly = TRUE)) {
        abort(str_c(
            "Package 'fdrtool' is required for selection = \"ldQuadratic\". ",
            "Install it, or use selection = \"heldOut\"."
        ))
    }
    fdrtool::fdrtool(
        corInput,
        statistic = "correlation",
        plot = FALSE,
        verbose = FALSE
    )$lfdr
}

# Each candidate's score against the LD shrunk toward the identity by `w`,
# c' b / sqrt((1 - w) b' R b + w b' b), and -Inf where that variance is not
# positive. At w = 0 it is lassosum's LD-quadratic score.
# @noRd
.rssShrunkLdScore <- function(candidateBeta, corInput, ld, w) {
    .rssCorrelationScore(.rssShrunkLdMoments(candidateBeta, corInput, ld, w))
}

# Each candidate's covariance with the statistics, c' b, and its variance
# against the LD shrunk toward the identity by `w`, (1 - w) b' R b + w b' b.
# @noRd
.rssShrunkLdMoments <- function(candidateBeta, corInput, ld, w) {
    list(
        bxy = as.numeric(crossprod(corInput, candidateBeta)),
        bxxb = (1 - w) *
            .ldQuadForms(ld, candidateBeta) +
            w * colSums(candidateBeta^2)
    )
}

# The correlation c' b / sqrt(b' Sigma b) from a candidate's moments, and -Inf
# where its variance is not positive.
# @noRd
.rssCorrelationScore <- function(moments) {
    positive <- is.finite(moments$bxxb) & moments$bxxb > 0
    replace(
        rep(-Inf, length(moments$bxy)),
        positive,
        moments$bxy[positive] / sqrt(moments$bxxb[positive])
    )
}

# How many times held-out selection splits the statistics, and the share of
# the samples each training part keeps.
# @noRd
.rssHeldOutSplits <- function() 4L

# @noRd
.rssHeldOutTrainShare <- function() 0.75

# Selection on held-out summary statistics, as PUMAS does it (Zhao et al.
# 2021). For the marginal correlations c of n samples, a training part of
# nTrain of them has c_t ~ N(c, (n - nTrain) / (n nTrain) Sigma), and its
# held-out remainder c_h = (n c - nTrain c_t) / (n - nTrain). Sigma is the
# LD shrunk toward the identity by the Schafer-Strimmer intensity for the
# reference's size: a reference smaller than the number of variants misses
# directions of the true LD, the least-penalized fits put weight exactly
# there, and draws and scores from R alone never see it. Each split refits
# the whole (s, lambda) grid on c_t and scores it on c_h; the full-data
# candidate with the best mean score is returned (.rssHeldOutPick()).
# @noRd
#' @importFrom stats rnorm
.rssSelectHeldOut <- function(candidateBeta, corInput, ld, n, refit) {
    w <- .ldShrinkIntensity(ld, ld$nSamples)
    share <- .rssHeldOutTrainShare()
    nSplits <- .rssHeldOutSplits()
    noise <- sqrt(1 - w) *
        .ldSourceDraws(ld, nSplits) +
        sqrt(w) * matrix(rnorm(length(corInput) * nSplits), ncol = nSplits)
    scores <- map(
        seq_len(nSplits),
        .rssHeldOutScore,
        noise = noise * sqrt((1 - share) / (share * n)),
        corInput = corInput,
        ld = ld,
        refit = refit,
        w = w
    )
    idx <- .rssHeldOutPick(
        rowMeans(exec(cbind, !!!map(scores, "cor"))),
        rowMeans(exec(cbind, !!!map(scores, "explained")))
    )
    list(beta = candidateBeta[, idx], index = idx, mode = "heldOut")
}

# Split `k`: the candidates refit on its training part and scored on its
# held-out part, by their correlation with it and by the variance of it they
# explain, 2 c_h' b - b' Sigma b, which unlike the correlation depends on b's
# scale. `refit` names the refitting function (`fn`) and its other
# arguments (`args`); it is given the training part's correlations and the
# share of the samples that part keeps, and returns the candidates in the
# full data's order.
# @noRd
.rssHeldOutScore <- function(k, noise, corInput, ld, refit, w) {
    share <- .rssHeldOutTrainShare()
    trainCor <- .lassosumClampCor(corInput + noise[, k])
    heldCor <- (corInput - share * trainCor) / (1 - share)
    beta <- exec(refit$fn, trainCor = trainCor, share = share, !!!refit$args)
    moments <- .rssShrunkLdMoments(beta, heldCor, ld, w)
    list(
        cor = .rssCorrelationScore(moments),
        explained = 2 * moments$bxy - moments$bxxb
    )
}

# The candidate with the best mean held-out correlation. Candidates that
# differ only in scale -- the ridge's at s = 1, where each fit is
# c / (1 + lambda) -- tie on it, and the one explaining the most held-out
# variance wins; the first in the grid where that ties too.
# @noRd
.rssHeldOutPick <- function(cor, explained) {
    ties <- .rssScoreTies(cor)
    ties[.lassosumFirstMax(explained[ties])]
}

# The shrinkage grid refit on a held-out split's training part.
# @noRd
.rssRefitShrinkGrid <- function(trainCor, share, s, ld, n, config) {
    fits <- map(
        s,
        .rssFitAtS,
        solverInput = trainCor * sqrt(share * n),
        ld = ld,
        n = share * n,
        config = config
    )
    exec(cbind, !!!map(fits, "beta"))
}

# Held-out selection draws the reference's sampling noise, so it needs the
# number of samples behind the LD.
# @noRd
.rssCheckHeldOutLd <- function(ld, label) {
    if (is.na(ld$nSamples)) {
        abort(c(
            glue(
                "{label}: selection = \"heldOut\" needs the LD reference's ",
                "sample size."
            ),
            i = str_c(
                "Build the LdMatrix with `nSamples`, or choose ",
                "selection = \"ldQuadratic\"."
            )
        ))
    }
    invisible(NULL)
}

# Validate the (bhat, R, n) inputs shared by the RSS solvers (lassosumRss /
# penalizedRss via .rssSolvePath, prsCs, and sdpr). R is a single LD correlation
# matrix over one cis-window, matching susieR::susie_rss; bhat must match
# nrow(R). missing() guards turn an omitted R / n into a clear message instead
# of an "argument is missing" error, and propagate through the solvers'
# callers (verified two levels deep). Method-specific checks -- prsCs's maf
# length, sdpr's M / perVariantSampleSize / array -- stay in the caller.
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

# .rssValidateInputs() for an LD source (.ldSourceRows()): the matrix
# checks for a matrix, and for a factor the same checks against its columns.
# @noRd
.rssValidateSource <- function(bhat, ld, n) {
    if (!is.null(ld$R)) {
        return(.rssValidateInputs(bhat, ld$R, n))
    }
    if (missing(n) || n <= 0) {
        abort("Please provide a valid sample size using 'n'.")
    }
    assertVector(bhat, len = ncol(ld$F))
    invisible(NULL)
}

# Validate (bhat, ld, n), run a decreasing-lambda coordinate-descent sweep via
# `solveFn` (the penalty-specific C++ solver), then reorder back to the input
# lambda order via the inverse permutation and assemble the standard result
# list. Shared by lassosumRss and penalizedRss, which differ only in which C++
# solver they pass as `solveFn` (and penalizedRss's per-penalty gamma default).
# `ld` is an LD source: the matrix, or a factor of it with the shrinkage the
# factor solver applies itself.
.rssSolvePath <- function(bhat, ld, n, lambda, solveFn, solveArgs = list()) {
    .rssValidateSource(bhat, ld, n)
    z <- bhat / sqrt(n)
    order <- order(lambda, decreasing = TRUE)
    # `solveArgs` rather than `...`: the two solvers take different fixed
    # argument sets, and both callers know theirs statically, so an unknown
    # solver argument should be an error here rather than reaching the solver.
    solved <- exec(solveFn, z, lambda[order], ld, !!!solveArgs)
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

# Per-`s` fit for one RSS method: one penalized path at the wrapper's
# penalty, or, for l0learn, one per value of its lambda0 grid. Returns
# list(beta, meta).
# @noRd
.rssFitOne <- function(solverInput, ld, n, sVal, config) {
    if (identical(config$method, "l0learn")) {
        return(.rssFitL0learn(solverInput, ld, n, sVal, config))
    }
    .rssFitPenalized(solverInput, ld, n, sVal, config)
}

# Shared per-candidate frame for a single-path model: the shrinkage s and
# each lambda, with the solver's objective, loss, convergence flag and number
# of non-zero coefficients there.
# @noRd
.rssPathMeta <- function(sVal, model) {
    tibble(
        s = rep(sVal, length(model$lambda)),
        lambda = model$lambda,
        fbeta = model$fbeta,
        loss = model$loss,
        conv = model$conv,
        nparams = model$nparams
    )
}

# The penalized path for one s, against `ld` shrunk to that s: the wrapper's
# penalty, `config$spec`, under the caller's path options.
# @noRd
.rssFitPenalized <- function(solverInput, ld, n, sVal, config) {
    model <- .penalizedRssOn(
        solverInput,
        ld,
        n,
        config$spec,
        config$options
    )
    list(beta = model$beta, meta = .rssPathMeta(sVal, model))
}

# L0Learn path for one s: one penalizedRss fit per lambda0, column-bound.
# @noRd
.rssFitL0learn <- function(solverInput, ld, n, sVal, config) {
    fits <- map(
        config$lambda0,
        .rssL0LambdaFit,
        solverInput = solverInput,
        ld = ld,
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
.rssFinalize <- function(bestBeta, sel, meta, config) {
    base <- c(mode = sel$mode, index = sel$index)
    selection <- switch(
        config$method,
        lassosum = c(
            base,
            s = meta$s[sel$index],
            lambda = meta$lambda[sel$index]
        ),
        penalized = c(
            base,
            penalty = config$spec$penalty,
            alpha = config$spec$alpha,
            s = meta$s[sel$index],
            lambda = meta$lambda[sel$index]
        ),
        l0learn = c(
            base,
            penalty = config$spec$penalty,
            s = meta$s[sel$index],
            lambda0 = meta$lambda0[sel$index],
            lambda = meta$lambda[sel$index]
        )
    )
    attrName <- if (config$method == "lassosum") {
        "lassosum_selection"
    } else {
        "penalized_rss_selection"
    }
    `attr<-`(bestBeta, attrName, selection)
}

# The LD the RSS weight methods read, matched to the statistics' variants
# (`stat$variantNames`): a factor of it for a genotype reference, which the
# factor solvers use without forming the p x p matrix, or the matrix an
# LdMatrix holds, unpacked once for the run: the solvers and the selection
# pass over it many times. Matching, allele flips, standardizing and padding
# are the LD reference layer's own (.ldFactorFromRef(), .ldFromRef()). The
# source also carries the number of samples behind the LD (.ldRefNoiseN(),
# NA for an LdMatrix that records none), which held-out selection needs.
# @noRd
.rssLdSource <- function(stat, ldReference, label) {
    ids <- .rssStatVariantIds(stat, label)
    ld <- if (.ldIsGenotypeRef(ldReference)) {
        list(F = .ldFactorFromRef(ldReference, ids, label))
    } else {
        list(R = .ldFromRef(ldReference, ids, label))
    }
    list_assign(ld, nSamples = .ldRefNoiseN(ldReference))
}

# The LD matrix itself, for the engines that take nothing else (mr.ash.rss,
# mr.mash.rss).
# @noRd
.rssLdMatrix <- function(stat, ldReference, label) {
    .ldFromRef(ldReference, .rssStatVariantIds(stat, label), label)
}

# The variant ids an LD reference is matched on.
# @noRd
.rssStatVariantIds <- function(stat, label) {
    ids <- stat$variantNames
    if (is.null(ids)) {
        abort(glue(
            "{label}: `stat` must carry `variantNames`, the variant ids ",
            "`ldReference` is matched on."
        ))
    }
    ids
}

# How many variants the statistics describe, for the weight extractors that
# check a fit against them: their variant ids, all of which must be in
# `ldReference` when one is given, or else their z-scores or effects (rows,
# for a multi-condition matrix). NA when the statistics carry none of these.
# @noRd
.rssVariantCount <- function(stat, ldReference, label) {
    if (!is.null(ldReference)) {
        ids <- .rssStatVariantIds(stat, label)
        .ldFromRefValidate(ldReference, label)
        .ldFromRefMatch(ldReference, ids, label, "error")
        return(length(ids))
    }
    described <- stat$variantNames %||% stat$z %||% stat$b
    if (is.null(described)) NA_integer_ else NROW(described)
}

# One shrinkage level's fit against R_s = (1 - s) R + s I.
# @noRd
.rssFitAtS <- function(sVal, solverInput, ld, n, config) {
    .rssFitOne(solverInput, .rssShrunkLd(ld, sVal), n, sVal, config)
}

# An LD source shrunk to R_s = (1 - s) R + s I: the source carries s to the
# solver, which applies it to the matrix or the factor without forming R_s,
# so the grid over s shares one LD matrix rather than copying it per s.
# @noRd
.rssShrunkLd <- function(ld, s) {
    list_assign(ld, shrink = s)
}

# Shared scaffold for the RSS shrinkage-grid weight functions
# (lassosumRssWeights / .penalizedRssWeights / l0learnRssWeights), run for the
# statistics `stat` against the LD source `ld` (see .rssLdSource()).
# Standardizes the stat -> solverInput conversion, the outer LD-shrinkage grid
# over `s`, the candidate accumulation, and the heldOut / ldQuadratic /
# minFbeta selection. `config` is the run: `method` (which selection record
# the weights carry; "l0learn" also sweeps `lambda0`), the penalty `spec`
# (.rssPenalty()), the caller's path `options` and its `label` for errors;
# `fitRetention` is how much of the search is kept.
# @noRd
.rssShrinkGridWeights <- function(
    stat,
    ld,
    s,
    config,
    selection = c("heldOut", "ldQuadratic", "minFbeta"),
    fitRetention = c("none", "slim", "full")
) {
    selection <- arg_match(selection)
    fitRetention <- arg_match(fitRetention)
    if (selection == "heldOut") {
        .rssCheckHeldOutLd(ld, config$label)
    }
    ld <- .rssSolverSource(ld, config$spec)
    n <- median(stat$n)
    p <- .ldSourceSize(ld)
    corInput <- .lassosumClampCor(.rssMarginalCor(stat, n = n, p = p))
    solverInput <- corInput * sqrt(n)
    fits <- map(
        s,
        .rssFitAtS,
        solverInput = solverInput,
        ld = ld,
        n = n,
        config = config
    )
    candidateBeta <- exec(cbind, !!!map(fits, "beta"))
    candidateMeta <- bind_rows(map(fits, "meta"))
    selectorResult <- .rssShrinkGridSelect(
        selection,
        candidateBeta,
        candidateMeta,
        corInput,
        ld,
        n,
        s,
        config
    )
    bestBeta <- as.numeric(selectorResult$beta)
    weights <- .rssFinalize(bestBeta, selectorResult, candidateMeta, config)
    .rssAttachGridFit(
        weights,
        candidateBeta,
        candidateMeta,
        selectorResult$index,
        fitRetention,
        stat[["variantNames"]]
    )
}

# The candidate `selection` picks from the grid: by held-out prediction,
# refitting the grid on a split of the statistics; by the LD quadratic form;
# or by the solver's smallest F-beta.
# @noRd
.rssShrinkGridSelect <- function(
    selection,
    candidateBeta,
    candidateMeta,
    corInput,
    ld,
    n,
    s,
    config
) {
    switch(
        selection,
        heldOut = .rssSelectHeldOut(
            candidateBeta,
            corInput,
            ld,
            n,
            list(
                fn = ".rssRefitShrinkGrid",
                args = list(s = s, ld = ld, n = n, config = config)
            )
        ),
        ldQuadratic = .lassosumSelectLdQuadratic(candidateBeta, corInput, ld),
        minFbeta = .lassosumSelectMinFbeta(candidateBeta, candidateMeta)
    )
}

# The search behind the grid weights as their "fit" attribute: `path`, one
# row per candidate fit with the solver's diagnostics and `selected` marking
# the one returned, and at fitRetention = "full" also `beta`, the candidates'
# coefficients (variants x candidates, in path's row order); it names the
# variants searched over, `variantIds`. Returns `weights` unchanged at
# fitRetention = "none".
# @noRd
.rssAttachGridFit <- function(
    weights,
    candidateBeta,
    candidateMeta,
    index,
    fitRetention,
    variantIds = NULL
) {
    if (identical(fitRetention, "none")) {
        return(weights)
    }
    fit <- c(
        list(path = mutate(candidateMeta, selected = row_number() == index)),
        compact(list(beta = if (fitRetention == "full") candidateBeta))
    )
    `attr<-`(weights, "fit", .withFitVariantIds(fit, variantIds))
}

#' Extract weights from lassosumRss with shrinkage grid search
#'
#' Searches over a grid of shrinkage parameters \code{s} (default:
#' \code{c(0.2, 0.5, 0.9, 1.0)}, matching the original lassosum and OTTERS).
#' For each \code{s}, the LD matrix is shrunk as \code{(1-s)*R + s*I}, then
#' \code{lassosumRss()} is called across the lambda path. One candidate
#' (\code{s}, lambda) is returned, by default the one that scores best on
#' held-out summary statistics.
#'
#' @details
#' \code{selection = "heldOut"} splits the statistics four times into a
#' training part, 75\% of the samples, and its held-out remainder, drawing
#' each split from the statistics and the LD as PUMAS does (Zhao et al.
#' 2021). The path is refit on each training part, every candidate is scored
#' on the held-out part by
#' \deqn{\frac{c_h^T \beta}{\sqrt{\beta^T \hat\Sigma \beta}}, \quad
#'   \hat\Sigma = (1 - w) R + w I,}
#' and the full-data fit at the (\code{s}, lambda) with the best mean score
#' is returned. \eqn{w} is the Schafer-Strimmer intensity for a correlation
#' matrix from the reference's sample size. A reference smaller than the
#' number of variants misses directions of the true LD, and the
#' least-penalized fits put weight exactly there, so scoring against
#' \eqn{R} itself would keep choosing them. The splits use R's random number
#' generator, so \code{set.seed()} makes the weights reproducible. The score
#' does not change with \eqn{\beta}'s scale, so candidates that differ only
#' in scale tie on it, as the ridge's do at \eqn{s = 1}, where each fit is
#' \eqn{c / (1 + \lambda)}. Of tied candidates, the one explaining the most
#' held-out variance, \eqn{2 c_h^T \beta - \beta^T \hat\Sigma \beta} averaged
#' over the splits, is returned, and the first in the grid where that ties
#' too.
#'
#' \code{selection = "ldQuadratic"} is lassosum's pseudovalidation (Mak et
#' al. 2017): the in-sample score
#' \eqn{\tilde c^T \beta / \sqrt{\beta^T R \beta}}, with each marginal
#' correlation shrunk by its local false discovery rate,
#' \eqn{\tilde c = c (1 - \mathrm{lfdr})}, from \code{fdrtool} as upstream
#' lassosum computes it. Of candidates tied on it, the first in the grid is
#' returned.
#'
#' @param stat A list of summary statistics: \code{n} (per-variant sample
#'   sizes) and the marginal effects as \code{cor} (standardized), \code{z},
#'   \code{b} with \code{seb}, or \code{b} alone, taken as already
#'   standardized. The engine receives the standardized effects z / sqrt(n),
#'   with n the median of \code{n}.
#' @param ldReference The LD reference the statistics' variants are matched
#'   to by \code{stat$variantNames}: a genotype panel (see
#'   \code{\link{readGenotypes}}), a \code{\link{GenotypeMixture}} or an
#'   \code{\link{LdMatrix}}. A genotype reference's LD is used through its
#'   standardized genotypes, never formed as a matrix.
#' @param backend \code{"pecotmr"} (default), pecotmr's lassosum, or
#'   \code{"bigsnpr"}, bigsnpr's lassosum2 (Prive et al. 2020) on the LD
#'   matrix, written to a temporary file for the call (formed from the
#'   genotypes for a genotype reference). lassosum2 fits the same objective
#'   with \eqn{\delta = s / (1 - s)}; its fits are rescaled by
#'   \eqn{1 + \delta} to lassosum's and recorded in lassosum's \code{s},
#'   and it traces its own lambda path down from the largest marginal
#'   effect. It reports no objective, so it has no
#'   \code{selection = "minFbeta"}, and its \code{"fit"} path counts each
#'   fit's \code{iterations} instead of \code{fbeta}, \code{loss} and
#'   \code{conv}. It is slow: \code{snp_lassosum2()} times each fit with
#'   \code{system.time()}, which first runs a full garbage collection of
#'   the R session, so a run costs one collection per grid point however
#'   few the variants, and \code{selection = "heldOut"} runs the grid five
#'   times. With the defaults that is 600 collections, minutes per call in
#'   a session with pecotmr loaded. \code{"pecotmr"} fits the same
#'   objective without that cost.
#' @param s Numeric vector of shrinkage parameters to search over.
#'   \code{NULL} (default) is \code{c(0.2, 0.5, 0.9, 1.0)} on pecotmr,
#'   following Mak et al (2017) and OTTERS, and lassosum2's own grid,
#'   \eqn{\delta = 0.001, 0.01, 0.1, 1}, on bigsnpr, where every \code{s}
#'   must lie strictly between 0 and 1.
#' @param selection How the returned candidate is chosen (see Details of
#'   \code{\link{lassosumRssWeights}}). \code{"heldOut"} (default) scores
#'   every candidate on held-out summary statistics and needs the
#'   reference's sample size, so an \code{\link{LdMatrix}} must record
#'   \code{nSamples}. \code{"ldQuadratic"} is lassosum's pseudovalidation and
#'   needs the \code{fdrtool} package. \code{"minFbeta"} picks the smallest
#'   objective, for debugging.
#' @param fitRetention How much of the search is kept as the \code{"fit"}
#'   attribute of the returned weights. \code{"none"} (default) attaches
#'   nothing. \code{"slim"} keeps \code{path}, a table with one row per
#'   candidate (\code{s}, \code{lambda}) giving the solver's objective
#'   \code{fbeta}, \code{loss}, convergence flag \code{conv} and number of
#'   non-zero coefficients \code{nparams}, with \code{selected} marking the
#'   candidate returned. \code{"full"} also keeps \code{beta}, the
#'   candidates' coefficients (variants x candidates, in the rows' order).
#' @param methodArgs Options forwarded to the solver: for pecotmr,
#'   \code{\link{lassosumRss}}'s, built with \code{\link{LassosumOptions}};
#'   for bigsnpr, \code{snp_lassosum2()}'s, built with
#'   \code{\link{BigsnprOptions}}. \code{NULL} (default) is the backend's
#'   constructor with its defaults. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be silently
#'   ignored.
#' @return A numeric vector of the posterior SNP coefficients at the best (s,
#'   lambda). Unless \code{fitRetention = "none"}, the search is attached as
#'   attribute \code{"fit"}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ss <- lapply(
#'   seq_len(ncol(X)), function(j) coef(summary(lm(y ~ X[, j])))[2, 1:2])
#' stat <- list(
#'   b = vapply(ss, `[`, numeric(1), 1L),
#'   seb = vapply(ss, `[`, numeric(1), 2L),
#'   n = rep(nrow(X), ncol(X)),
#'   variantNames = colnames(X)
#' )
#' ld <- LdMatrix(cor(X), nSamples = nrow(X))
#' set.seed(1)
#' w <- lassosumRssWeights(stat, ld, fitRetention = "slim")
#' head(attr(w, "fit")$path)
#' @export
lassosumRssWeights <- function(
    stat,
    ldReference,
    backend = NULL,
    s = NULL,
    selection = c("heldOut", "ldQuadratic", "minFbeta"),
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    selection <- arg_match(selection)
    fitRetention <- arg_match(fitRetention)
    backend <- .rssBackend("lasso", backend, ldReference)
    methodArgs <- .backendOptions(
        methodArgs,
        .methodBackends$sumstat$lasso[[backend]],
        backend
    )
    if (identical(backend, "bigsnpr")) {
        return(.lassosum2Weights(
            stat,
            ldReference,
            s,
            selection,
            fitRetention,
            methodArgs
        ))
    }
    .rssShrinkGridWeights(
        stat,
        .rssLdSource(stat, ldReference, "lassosumRssWeights"),
        s %||% c(0.2, 0.5, 0.9, 1.0),
        list(
            method = "lassosum",
            spec = .rssPenalty("lasso"),
            options = as.list(methodArgs),
            label = "lassosumRssWeights"
        ),
        selection,
        fitRetention
    )
}

# lassosum on bigsnpr's lassosum2, which minimizes
# b'(R + delta I)b - 2 b'bhat + 2 lambda |b|_1. That is lassosum's objective
# at s = delta / (1 + delta) divided by 1 - s, with bhat and lambda both
# scaled by 1 / (1 - s); a lasso solution scales with the two together, so at
# the same lambda lassosum's is lassosum2's times 1 + delta. Each fit is
# rescaled by that and recorded in lassosum's s.
# `s` gives delta = s / (1 - s), or lassosum2's own default when NULL; the
# lambda path is lassosum2's, from the largest marginal effect, so held-out
# selection matches the candidates by their place on it. lassosum2 reports no
# objective, so it has no "minFbeta".
# @noRd
.lassosum2Weights <- function(
    stat,
    ldReference,
    s,
    selection,
    fitRetention,
    methodArgs
) {
    label <- "lassosumRssWeights"
    if (identical(selection, "minFbeta")) {
        abort(glue(
            "{label}: lassosum2 reports no objective, so backend = ",
            "\"bigsnpr\" has no selection = \"minFbeta\"."
        ))
    }
    inp <- .rssPerVariantInput(stat, ldReference, label)
    if (identical(selection, "heldOut")) {
        .rssCheckHeldOutLd(inp$ld, label)
    }
    engine <- .bigsnprEngine(inp)
    on.exit(.bigsnprRelease(engine), add = TRUE)
    config <- list(delta = .lassosum2Delta(s), options = as.list(methodArgs))
    fits <- .lassosum2Run(inp, engine, config)
    .withWeightsFitIds(
        .lassosum2WeightsFrom(
            fits,
            inp,
            engine,
            config,
            selection,
            fitRetention
        ),
        stat[["variantNames"]]
    )
}

# lassosum2's delta from lassosum's s, or its own default grid.
# @noRd
.lassosum2Delta <- function(s) {
    if (is.null(s)) {
        return(.formalDefaults(bigsnpr::snp_lassosum2, "delta")[[1L]])
    }
    assertNumeric(s, lower = 0, upper = 1, any.missing = FALSE, min.len = 1L)
    if (any(s <= 0 | s >= 1)) {
        abort(str_c(
            "With backend = \"bigsnpr\", every `s` must lie strictly between ",
            "0 and 1: lassosum2's delta = s / (1 - s) must be positive and ",
            "finite."
        ))
    }
    s / (1 - s)
}

# lassosum2's fits, rescaled to lassosum's parametrization, with their grid
# in lassosum's s and lambda. A fit that diverged is NA.
# @noRd
.lassosum2Run <- function(inp, engine, config) {
    fits <- exec(
        bigsnpr::snp_lassosum2,
        engine$corr,
        .bigsnprDfBeta(inp),
        delta = config$delta,
        ncores = 1,
        !!!config$options
    )
    grid <- attr(fits, "grid_param")
    list(
        beta = sweep(
            fits[, seq_len(nrow(grid)), drop = FALSE],
            2L,
            1 + grid$delta,
            "*"
        ),
        path = tibble(
            s = grid$delta / (1 + grid$delta),
            lambda = grid$lambda,
            iterations = grid$num_iter
        )
    )
}

# The fit chosen by `selection`, recorded as lassosumRssWeights() records
# pecotmr's, with the search at the requested retention.
# @noRd
.lassosum2WeightsFrom <- function(
    fits,
    inp,
    engine,
    config,
    selection,
    fitRetention
) {
    usable <- replace(fits$beta, is.na(fits$beta), 0)
    sel <- if (identical(selection, "ldQuadratic")) {
        .lassosumSelectLdQuadratic(usable, inp$bhat, inp$ld)
    } else {
        .rssSelectHeldOut(
            usable,
            inp$bhat,
            inp$ld,
            median(inp$n),
            list(
                fn = ".lassosum2Refit",
                args = list(inp = inp, engine = engine, config = config)
            )
        )
    }
    weights <- `attr<-`(
        as.numeric(sel$beta),
        "lassosum_selection",
        c(
            mode = sel$mode,
            index = sel$index,
            s = fits$path$s[[sel$index]],
            lambda = fits$path$lambda[[sel$index]]
        )
    )
    path <- mutate(fits$path, nparams = as.integer(colSums(usable != 0)))
    .rssAttachGridFit(weights, fits$beta, path, sel$index, fitRetention)
}

# lassosum2 refit on a held-out split's training part; a fit that diverges
# there is zero, so it scores -Inf.
# @noRd
.lassosum2Refit <- function(trainCor, share, inp, engine, config) {
    train <- list_assign(inp, bhat = trainCor, n = share * inp$n)
    fits <- .lassosum2Run(train, engine, config)$beta
    replace(fits, is.na(fits), 0)
}

#' Penalized Regression on RSS (Summary Statistics) Objective
#'
#' The RSS penalized-regression solver: LASSO with its elastic-net and ridge
#' mixtures, MCP, SCAD, L0, L0L1 and L0L2. Uses coordinate descent on the
#' objective
#' \deqn{\beta^T R \beta - 2 \beta^T z + \mathrm{penalty}(\beta)}
#' where \eqn{R} is a (possibly pre-shrunk) LD matrix and \eqn{z = \hat\beta /
#' \sqrt{n}}. \code{lassosumRss()} is its LASSO at \code{alpha = 1}. The
#' ridge (LASSO at \code{alpha = 0}) has a closed form instead,
#' \eqn{(R + \lambda I)^{-1} z}, which is read for the whole path from one
#' eigendecomposition of \eqn{R}.
#'
#' Internal: \code{\link{enetRssWeights}}, \code{\link{ridgeRssWeights}},
#' \code{\link{scadRssWeights}}, \code{\link{mcpRssWeights}} and
#' \code{\link{l0learnRssWeights}} are the public interface, with
#' \code{\link{lassosumRssWeights}} for the lasso. They run this over a grid
#' of LD shrinkage values, forward its options through
#' \code{\link{PenalizedRssOptions}}, and with \code{fitRetention} return
#' the paths as the \code{"fit"} attribute of their weights.
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
#'   component. \code{NULL} (default) is 20 values from 1e-4 to 0.1 on a log
#'   scale, divided by \code{max(alpha, 0.001)} as glmnet sizes its path, so
#'   the L1 share spans the same range at any \code{alpha} and a ridge path
#'   runs from 0.1 to 100.
#' @param gamma Concavity parameter for MCP (default 3) or SCAD (default 3.7).
#'   Ignored for LASSO and L0 variants.
#' @param alpha Elastic-net mixing for LASSO, MCP and SCAD, in \code{[0, 1]}:
#'   \eqn{l_1 = \lambda \alpha} and \eqn{l_2 = \lambda (1-\alpha)}, so the
#'   penalty is \eqn{2 P(\beta; l_1) + l_2 \beta^T \beta}, glmnet's mixing.
#'   Default 1, no ridge term. For LASSO, 0.5 is the elastic net and 0 ridge;
#'   MCP and SCAD need \code{alpha > 0}. Ignored by the L0 variants.
#' @param lambda0 L0 penalty weight (number of non-zeros). Required for L0
#'   variants; ignored otherwise. Default 0.
#' @param lambda2 L2 penalty weight for L0L2 variant. Default 0.
#' @param thr Convergence threshold. Default 1e-4. The ridge, solved
#'   exactly, does not use it.
#' @param maxiter Maximum coordinate descent iterations per lambda. Default
#'   10000. The ridge, solved exactly, does not use it.
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
#' @keywords internal
penalizedRss <- function(
    bhat,
    R,
    n,
    penalty = c("lasso", "MCP", "SCAD", "L0", "L0L1", "L0L2"),
    lambda = NULL,
    gamma = NULL,
    alpha = 1.0,
    lambda0 = 0,
    lambda2 = 0,
    thr = 1e-4,
    maxiter = 10000,
    maxSwaps = 100
) {
    .penalizedRssPath(
        bhat,
        list(R = R),
        n,
        spec = .rssPenalty(
            arg_match(penalty),
            alpha = alpha,
            gamma = gamma,
            lambda0 = lambda0,
            lambda2 = lambda2,
            maxSwaps = maxSwaps
        ),
        lambda = lambda,
        thr = thr,
        maxiter = maxiter
    )
}

# The penalty a fit minimizes, as one checked record: the solver's penalty
# name, the elastic-net mixing `alpha`, the MCP / SCAD concavity `gamma`
# (their usual default when NULL) and the L0 settings. Each weight wrapper
# builds it from its own arguments before its shrinkage grid starts, so the
# penalty is the wrapper's alone and only the path (lambda, thr, maxiter) is
# left to the caller's options.
# @noRd
#' @importFrom checkmate assertCount assertNumber
.rssPenalty <- function(
    penalty,
    alpha = 1,
    gamma = NULL,
    lambda0 = 0,
    lambda2 = 0,
    maxSwaps = 0
) {
    .rssCheckAlpha(penalty, alpha)
    assertNumber(gamma, finite = TRUE, null.ok = TRUE)
    assertNumber(lambda0, lower = 0, finite = TRUE)
    assertNumber(lambda2, lower = 0, finite = TRUE)
    assertCount(maxSwaps)
    list(
        penalty = penalty,
        alpha = alpha,
        gamma = gamma %||% switch(penalty, SCAD = 3.7, MCP = 3.0, 0.0),
        lambda0 = lambda0,
        lambda2 = lambda2,
        maxSwaps = maxSwaps
    )
}

# `alpha` for `penalty`: a share in [0, 1], and above 0 for MCP and SCAD,
# where 0 leaves only the ridge term, which ridgeRssWeights() fits.
# @noRd
.rssCheckAlpha <- function(penalty, alpha) {
    assertNumber(alpha, lower = 0, upper = 1)
    if (alpha == 0 && is_in(penalty, c("MCP", "SCAD"))) {
        abort(glue(
            "`alpha` must be above 0 for {penalty}: at 0 only the ridge term ",
            "is left, which ridgeRssWeights() fits."
        ))
    }
    invisible(NULL)
}

# penalizedRss() on an LD source: the matrix, or a factor of it with its
# shrinkage, which the factor solver applies without forming R. `spec` is the
# penalty (.rssPenalty()); the rest is the path, its default sized by the L1
# share.
# @noRd
.penalizedRssPath <- function(bhat, ld, n, spec, lambda, thr, maxiter) {
    lambda <- lambda %||% .rssDefaultLambda(spec$alpha)
    assertNumeric(lambda, lower = 0, any.missing = FALSE)
    assertNumber(thr, lower = 0, finite = TRUE)
    assertCount(maxiter, positive = TRUE)
    .rssSolvePath(
        bhat,
        ld,
        n,
        lambda,
        .rssPenalizedSolve,
        solveArgs = c(spec, list(thr = thr, maxiter = maxiter))
    )
}

# The default lambda path: 20 values from 1e-4 to 0.1, sized for an L1
# penalty on correlation-scale effects, divided by max(alpha, 0.001) as glmnet
# sizes its own. The L1 share then spans the same range at any alpha, and a
# ridge path (alpha = 0) runs from 0.1 to 100; at alpha = 1 it is the
# lassosum default unchanged.
# @noRd
.rssDefaultLambda <- function(alpha) {
    exp(seq(log(1e-4), log(0.1), length.out = 20L)) / max(alpha, 1e-3)
}

# .penalizedRssPath() for the penalty `spec` under the caller's path
# `options` -- lambda, thr and maxiter, all an options bag may hold --
# penalizedRss()'s own defaults filling the rest. A matrix takes exactly the
# path penalizedRss() does.
# @noRd
.penalizedRssOn <- function(bhat, ld, n, spec, options) {
    path <- list_modify(
        .formalDefaults(penalizedRss, c("lambda", "thr", "maxiter")),
        !!!options
    )
    exec(
        .penalizedRssPath,
        bhat = bhat,
        ld = ld,
        n = n,
        spec = spec,
        !!!path
    )
}

#' RSS Weights Helper for Penalized Methods
#'
#' Shared implementation for \code{scadRssWeights()}, \code{mcpRssWeights()},
#' \code{enetRssWeights()} and \code{ridgeRssWeights()}.
#' (\code{l0learnRssWeights()} uses the lower-level \code{.rssShrinkGridWeights}
#' scaffold directly, since it sweeps an additional \code{lambda0} path.)
#' Searches over a shrinkage grid \code{s} (LD matrix shrinkage \code{(1-s)R +
#' sI}) and selects the best candidate on held-out summary statistics, by
#' lassosum's pseudovalidation, or by the minimum penalized objective.
#'
#' @param stat,ldReference,s,selection See the public wrappers for details.
#' @param spec The penalty, built by the calling wrapper from its own
#'   arguments (internal \code{.rssPenalty()}).
#' @param fitRetention How much of the search is kept as the \code{"fit"}
#'   attribute; see \code{\link{scadRssWeights}}.
#' @param methodArgs Options forwarded to \code{penalizedRss}, built with
#'   \code{\link{PenalizedRssOptions}}: the path only, since the penalty is
#'   the wrapper's.
#' @param label The calling wrapper's name, for its error messages.
#' @return Numeric weight vector, one weight per variant of \code{stat}.
#' @keywords internal
.penalizedRssWeights <- function(
    stat,
    ldReference,
    spec,
    s = c(0.2, 0.5, 0.9, 1.0),
    selection = c("heldOut", "ldQuadratic", "minFbeta"),
    fitRetention = c("none", "slim", "full"),
    methodArgs = PenalizedRssOptions(),
    label = str_c(str_to_lower(spec$penalty), "RssWeights")
) {
    .assertMethodOptions(methodArgs, "PenalizedRssOptions", "methodArgs")
    selection <- arg_match(selection)
    fitRetention <- arg_match(fitRetention)
    .rssShrinkGridWeights(
        stat,
        .rssLdSource(stat, ldReference, label),
        s,
        list(
            method = "penalized",
            spec = spec,
            options = as.list(methodArgs),
            label = label
        ),
        selection,
        fitRetention
    )
}

#' Compute SCAD-Penalized Weights from Summary Statistics
#'
#' Fits SCAD-penalized regression on the RSS objective, searching over a
#' shrinkage grid \code{s} and lambda path. Model selection uses held-out
#' summary statistics by default (see \code{\link{lassosumRssWeights}}).
#'
#' @param stat A list of summary statistics: \code{n} (per-variant sample
#'   sizes) and the marginal effects as \code{cor} (standardized), \code{z},
#'   \code{b} with \code{seb}, or \code{b} alone, taken as already
#'   standardized. The engine receives the standardized effects z / sqrt(n),
#'   with n the median of \code{n}.
#' @param ldReference The LD reference the statistics' variants are matched
#'   to by \code{stat$variantNames}: a genotype panel (see
#'   \code{\link{readGenotypes}}), a \code{\link{GenotypeMixture}} or an
#'   \code{\link{LdMatrix}}. A genotype reference's LD is used through its
#'   standardized genotypes, never formed as a matrix.
#' @param s Numeric vector of LD shrinkage parameters. Default: \code{c(0.2,
#'   0.5, 0.9, 1.0)}.
#' @param gamma SCAD concavity parameter. Default 3.7.
#' @param alpha Elastic-net mixing in \code{(0, 1]}: the SCAD share of
#'   the penalty, the rest a ridge term, as \code{\link{penalizedRss}}
#'   mixes them. Default 1, no ridge term.
#' @param methodArgs Options forwarded to the solver,
#'   \code{\link{penalizedRss}}, built with
#'   \code{\link{PenalizedRssOptions}}. A bare list is refused: it cannot
#'   be checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @inheritParams lassosumRssWeights
#' @return A numeric vector of SNP coefficient weights. Unless
#'   \code{fitRetention = "none"}, the search is attached as attribute
#'   \code{"fit"}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ss <- lapply(
#'   seq_len(ncol(X)), function(j) coef(summary(lm(y ~ X[, j])))[2, 1:2])
#' stat <- list(
#'   b = vapply(ss, `[`, numeric(1), 1L),
#'   seb = vapply(ss, `[`, numeric(1), 2L),
#'   n = rep(nrow(X), ncol(X)),
#'   variantNames = colnames(X)
#' )
#' scadRssWeights(stat, LdMatrix(cor(X), nSamples = nrow(X)))
#' @importFrom checkmate assertList assertNumeric assertNumber
#' @export
scadRssWeights <- function(
    stat,
    ldReference,
    s = c(0.2, 0.5, 0.9, 1.0),
    gamma = 3.7,
    alpha = 1.0,
    selection = c("heldOut", "ldQuadratic", "minFbeta"),
    fitRetention = c("none", "slim", "full"),
    methodArgs = PenalizedRssOptions()
) {
    .assertMethodOptions(methodArgs, "PenalizedRssOptions", "methodArgs")
    fitRetention <- arg_match(fitRetention)
    assertList(stat)
    assertNumeric(s, lower = 0, any.missing = FALSE)
    # Built before the grid starts, so a bad penalty fails here.
    spec <- .rssPenalty("SCAD", alpha = alpha, gamma = gamma)
    .penalizedRssWeights(
        stat = stat,
        ldReference = ldReference,
        spec = spec,
        s = s,
        selection = selection,
        fitRetention = fitRetention,
        methodArgs = methodArgs
    )
}

#' Compute MCP-Penalized Weights from Summary Statistics
#'
#' Fits MCP-penalized regression on the RSS objective, searching over a
#' shrinkage grid \code{s} and lambda path. Model selection uses held-out
#' summary statistics by default (see \code{\link{lassosumRssWeights}}).
#'
#' @param stat A list of summary statistics: \code{n} (per-variant sample
#'   sizes) and the marginal effects as \code{cor} (standardized), \code{z},
#'   \code{b} with \code{seb}, or \code{b} alone, taken as already
#'   standardized. The engine receives the standardized effects z / sqrt(n),
#'   with n the median of \code{n}.
#' @param ldReference The LD reference the statistics' variants are matched
#'   to by \code{stat$variantNames}: a genotype panel (see
#'   \code{\link{readGenotypes}}), a \code{\link{GenotypeMixture}} or an
#'   \code{\link{LdMatrix}}. A genotype reference's LD is used through its
#'   standardized genotypes, never formed as a matrix.
#' @param s Numeric vector of LD shrinkage parameters. Default: \code{c(0.2,
#'   0.5, 0.9, 1.0)}.
#' @param gamma MCP concavity parameter. Default 3.
#' @param alpha Elastic-net mixing in \code{(0, 1]}: the MCP share of
#'   the penalty, the rest a ridge term, as \code{\link{penalizedRss}}
#'   mixes them. Default 1, no ridge term.
#' @param methodArgs Options forwarded to the solver,
#'   \code{\link{penalizedRss}}, built with
#'   \code{\link{PenalizedRssOptions}}. A bare list is refused: it cannot
#'   be checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @inheritParams lassosumRssWeights
#' @return A numeric vector of SNP coefficient weights. Unless
#'   \code{fitRetention = "none"}, the search is attached as attribute
#'   \code{"fit"}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ss <- lapply(
#'   seq_len(ncol(X)), function(j) coef(summary(lm(y ~ X[, j])))[2, 1:2])
#' stat <- list(
#'   b = vapply(ss, `[`, numeric(1), 1L),
#'   seb = vapply(ss, `[`, numeric(1), 2L),
#'   n = rep(nrow(X), ncol(X)),
#'   variantNames = colnames(X)
#' )
#' mcpRssWeights(stat, LdMatrix(cor(X), nSamples = nrow(X)))
#' @importFrom checkmate assertList assertNumeric assertNumber
#' @export
mcpRssWeights <- function(
    stat,
    ldReference,
    s = c(0.2, 0.5, 0.9, 1.0),
    gamma = 3.0,
    alpha = 1.0,
    selection = c("heldOut", "ldQuadratic", "minFbeta"),
    fitRetention = c("none", "slim", "full"),
    methodArgs = PenalizedRssOptions()
) {
    .assertMethodOptions(methodArgs, "PenalizedRssOptions", "methodArgs")
    fitRetention <- arg_match(fitRetention)
    assertList(stat)
    assertNumeric(s, lower = 0, any.missing = FALSE)
    # Built before the grid starts, so a bad penalty fails here.
    spec <- .rssPenalty("MCP", alpha = alpha, gamma = gamma)
    .penalizedRssWeights(
        stat = stat,
        ldReference = ldReference,
        spec = spec,
        s = s,
        selection = selection,
        fitRetention = fitRetention,
        methodArgs = methodArgs
    )
}

#' Compute Elastic-Net Weights from Summary Statistics
#'
#' Fits the elastic net on the RSS objective -- the lasso penalty mixed with
#' a ridge term, as glmnet mixes them -- searching over a shrinkage grid
#' \code{s} and a lambda path. Model selection uses held-out summary
#' statistics by default (see \code{\link{lassosumRssWeights}}).
#'
#' \code{alpha} is the L1 share of the penalty: with \eqn{l_1 = \lambda
#' \alpha} and \eqn{l_2 = \lambda (1 - \alpha)} the penalty is
#' \eqn{2 l_1 \|\beta\|_1 + l_2 \beta^T \beta}, so \code{alpha = 1} is the
#' lasso and \code{alpha = 0} ridge (\code{\link{ridgeRssWeights}}). The
#' default lambda path is divided by \code{max(alpha, 0.001)}, as glmnet
#' sizes its own, so the L1 share spans the same range at any \code{alpha}.
#'
#' The LD shrinkage adds a ridge term of its own, since
#' \eqn{\beta^T R_s \beta = (1 - s) \beta^T R \beta + s \beta^T \beta}: every
#' RSS lasso is already partly a ridge, and \code{alpha} adds a second ridge
#' term that grows with \eqn{\lambda}. Expect the lasso and the elastic net
#' to differ less on summary statistics than glmnet's do on genotypes.
#'
#' @param alpha Elastic-net mixing in \code{[0, 1]}: the L1 share of the
#'   penalty. Default 0.5.
#' @inheritParams scadRssWeights
#' @inheritParams lassosumRssWeights
#' @return A numeric vector of SNP coefficient weights. Unless
#'   \code{fitRetention = "none"}, the search is attached as attribute
#'   \code{"fit"}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ss <- lapply(
#'   seq_len(ncol(X)), function(j) coef(summary(lm(y ~ X[, j])))[2, 1:2])
#' stat <- list(
#'   b = vapply(ss, `[`, numeric(1), 1L),
#'   seb = vapply(ss, `[`, numeric(1), 2L),
#'   n = rep(nrow(X), ncol(X)),
#'   variantNames = colnames(X)
#' )
#' enetRssWeights(stat, LdMatrix(cor(X), nSamples = nrow(X)))
#' # A larger ridge share:
#' enetRssWeights(stat, LdMatrix(cor(X), nSamples = nrow(X)), alpha = 0.2)
#' @importFrom checkmate assertList assertNumeric
#' @export
enetRssWeights <- function(
    stat,
    ldReference,
    alpha = 0.5,
    s = c(0.2, 0.5, 0.9, 1.0),
    selection = c("heldOut", "ldQuadratic", "minFbeta"),
    fitRetention = c("none", "slim", "full"),
    methodArgs = PenalizedRssOptions()
) {
    .assertMethodOptions(methodArgs, "PenalizedRssOptions", "methodArgs")
    fitRetention <- arg_match(fitRetention)
    assertList(stat)
    assertNumeric(s, lower = 0, any.missing = FALSE)
    # Built before the grid starts, so a bad mixing fails here.
    spec <- .rssPenalty("lasso", alpha = alpha)
    .penalizedRssWeights(
        stat = stat,
        ldReference = ldReference,
        spec = spec,
        s = s,
        selection = selection,
        fitRetention = fitRetention,
        methodArgs = methodArgs,
        label = "enetRssWeights"
    )
}

#' Compute Ridge Weights from Summary Statistics
#'
#' Fits ridge regression on the RSS objective, \code{\link{enetRssWeights}}
#' at \code{alpha = 0}: the penalty is \eqn{\lambda \beta^T \beta}, and the
#' default lambda path runs from 0.1 to 100. Searches over a shrinkage grid
#' \code{s} as well, whose own ridge term \eqn{s \beta^T \beta} adds to it;
#' model selection uses held-out summary statistics by default (see
#' \code{\link{lassosumRssWeights}}).
#'
#' Every fit is exact: the ridge has the closed form
#' \eqn{((1 - s) R + (s + \lambda) I)^{-1} z}, read from one
#' eigendecomposition of the window's LD that every shrinkage level, lambda
#' and held-out refit shares (for a genotype reference, from the singular
#' value decomposition of its standardized genotypes, never forming the LD).
#' The coordinate-descent settings \code{thr} and \code{maxiter} in
#' \code{methodArgs} therefore do not apply.
#'
#' @inheritParams scadRssWeights
#' @inheritParams lassosumRssWeights
#' @return A numeric vector of SNP coefficient weights. Unless
#'   \code{fitRetention = "none"}, the search is attached as attribute
#'   \code{"fit"}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ss <- lapply(
#'   seq_len(ncol(X)), function(j) coef(summary(lm(y ~ X[, j])))[2, 1:2])
#' stat <- list(
#'   b = vapply(ss, `[`, numeric(1), 1L),
#'   seb = vapply(ss, `[`, numeric(1), 2L),
#'   n = rep(nrow(X), ncol(X)),
#'   variantNames = colnames(X)
#' )
#' ridgeRssWeights(stat, LdMatrix(cor(X), nSamples = nrow(X)))
#' @importFrom checkmate assertList assertNumeric
#' @export
ridgeRssWeights <- function(
    stat,
    ldReference,
    s = c(0.2, 0.5, 0.9, 1.0),
    selection = c("heldOut", "ldQuadratic", "minFbeta"),
    fitRetention = c("none", "slim", "full"),
    methodArgs = PenalizedRssOptions()
) {
    .assertMethodOptions(methodArgs, "PenalizedRssOptions", "methodArgs")
    fitRetention <- arg_match(fitRetention)
    assertList(stat)
    assertNumeric(s, lower = 0, any.missing = FALSE)
    .penalizedRssWeights(
        stat = stat,
        ldReference = ldReference,
        spec = .rssPenalty("lasso", alpha = 0),
        s = s,
        selection = selection,
        fitRetention = fitRetention,
        methodArgs = methodArgs,
        label = "ridgeRssWeights"
    )
}

#' Compute L0-Penalized Weights from Summary Statistics
#'
#' Fits L0-penalized regression (with optional L1/L2 components) on the RSS
#' objective, searching over a shrinkage grid \code{s} and lambda0 path. Model
#' selection uses held-out summary statistics by default (see
#' \code{\link{lassosumRssWeights}}).
#'
#' The swap optimization from L0Learn is included: after coordinate descent
#' converges, non-zero coefficients are tested for swaps with zero ones to
#' escape local optima.
#'
#' @param stat A list of summary statistics: \code{n} (per-variant sample
#'   sizes) and the marginal effects as \code{cor} (standardized), \code{z},
#'   \code{b} with \code{seb}, or \code{b} alone, taken as already
#'   standardized. The engine receives the standardized effects z / sqrt(n),
#'   with n the median of \code{n}.
#' @param ldReference The LD reference the statistics' variants are matched
#'   to by \code{stat$variantNames}: a genotype panel (see
#'   \code{\link{readGenotypes}}), a \code{\link{GenotypeMixture}} or an
#'   \code{\link{LdMatrix}}. A genotype reference's LD is used through its
#'   standardized genotypes, never formed as a matrix.
#' @param penalty L0 variant: \code{"L0"}, \code{"L0L1"}, or \code{"L0L2"}.
#'   Default \code{"L0"}.
#' @param s Numeric vector of LD shrinkage parameters. Default: \code{c(0.2,
#'   0.5, 0.9, 1.0)}.
#' @param lambda0 Numeric vector of L0 penalty values to search over. Default:
#'   \code{exp(seq(log(0.001), log(1), length.out = 10))}.
#' @param lambda Numeric vector of L1 penalty values (for L0L1). Default:
#'   \code{c(0)} (no L1 unless L0L1 is used).
#' @param lambda2 L2 penalty weight (for L0L2). Default 0.
#' @param maxSwaps Maximum swap rounds per lambda. Default 100.
#' @param fitRetention How much of the search is kept as the \code{"fit"}
#'   attribute of the returned weights. \code{"none"} (default) attaches
#'   nothing. \code{"slim"} keeps \code{path}, a table with one row per
#'   candidate (\code{s}, \code{lambda0}, \code{lambda}) giving the solver's
#'   objective \code{fbeta}, \code{loss}, convergence flag \code{conv} and
#'   number of non-zero coefficients \code{nparams}, with \code{selected}
#'   marking the candidate returned. \code{"full"} also keeps \code{beta},
#'   the candidates' coefficients (variants x candidates, in the rows' order).
#' @param methodArgs Options forwarded to the solver,
#'   \code{\link{penalizedRss}}, built with
#'   \code{\link{PenalizedRssOptions}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be silently
#'   ignored.
#' @inheritParams lassosumRssWeights
#' @return A numeric vector of SNP coefficient weights. Unless
#'   \code{fitRetention = "none"}, the search is attached as attribute
#'   \code{"fit"}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ss <- lapply(
#'   seq_len(ncol(X)), function(j) coef(summary(lm(y ~ X[, j])))[2, 1:2])
#' stat <- list(
#'   b = vapply(ss, `[`, numeric(1), 1L),
#'   seb = vapply(ss, `[`, numeric(1), 2L),
#'   n = rep(nrow(X), ncol(X)),
#'   variantNames = colnames(X)
#' )
#' l0learnRssWeights(stat, LdMatrix(cor(X), nSamples = nrow(X)))
#' @export
l0learnRssWeights <- function(
    stat,
    ldReference,
    penalty = c("L0", "L0L1", "L0L2"),
    s = c(0.2, 0.5, 0.9, 1.0),
    lambda0 = exp(seq(log(0.001), log(1), length.out = 10)),
    lambda = NULL,
    lambda2 = 0,
    selection = c("heldOut", "ldQuadratic", "minFbeta"),
    maxSwaps = 100,
    fitRetention = c("none", "slim", "full"),
    methodArgs = PenalizedRssOptions()
) {
    .assertMethodOptions(methodArgs, "PenalizedRssOptions", "methodArgs")
    # Its own `lambda` is the path; one in the options would replace it.
    .configRefuseOwned(
        methodArgs,
        c(lambda = "l0learnRssWeights' own `lambda`"),
        "l0learnRssWeights"
    )
    penalty <- arg_match(penalty)
    selection <- arg_match(selection)
    fitRetention <- arg_match(fitRetention)
    # The L0 penalty, lambda0 left to the sweep over its grid.
    spec <- .rssPenalty(penalty, lambda2 = lambda2, maxSwaps = maxSwaps)
    assertNumeric(lambda0, lower = 0, finite = TRUE, any.missing = FALSE)
    lambda <- .l0learnRssLambda(lambda, penalty)
    .rssShrinkGridWeights(
        stat,
        .rssLdSource(stat, ldReference, "l0learnRssWeights"),
        s,
        list(
            method = "l0learn",
            spec = spec,
            lambda0 = lambda0,
            options = list_assign(as.list(methodArgs), lambda = lambda),
            label = "l0learnRssWeights"
        ),
        selection,
        fitRetention
    )
}

# The L1 component's grid: the caller's, else ten values on a log scale for
# L0L1 and none for the penalties without one.
# @noRd
.l0learnRssLambda <- function(lambda, penalty) {
    if (!is.null(lambda)) {
        return(lambda)
    }
    if (penalty == "L0L1") {
        exp(seq(log(0.0001), log(0.1), length.out = 10))
    } else {
        0
    }
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
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights: \code{"none"} (default) attaches
#'   nothing; \code{"slim"} the selected fit, as a list of its coefficients
#'   (\code{coef}, over the variants fitted), \code{intercept},
#'   \code{lambda} and cross-validated error there (\code{cve},
#'   \code{cvse}); \code{"full"} the whole \code{cv.ncvreg} object, with
#'   every \code{lambda}'s coefficients but without its copy of \code{X}.
#'   Either names the variants it was fitted on.
#' @param methodArgs Options forwarded to \code{ncvreg::cv.ncvreg}, built
#'   with \code{\link{NcvregOptions}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of length `ncol(X)` of variant weights, with the
#'   fit attached as attribute \code{"fit"} unless
#'   \code{fitRetention = "none"}.
#' @importFrom stats coef
#' @keywords internal
#' @importFrom checkmate assertCount
ncvregWeights <- function(
    X,
    y,
    penalty,
    nfolds = 5,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NcvregOptions()
) {
    .assertMethodOptions(methodArgs, "NcvregOptions", "methodArgs")
    assertCount(nfolds, positive = TRUE)
    fitRetention <- arg_match(fitRetention)
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
        !!!as.list(methodArgs)
    )
    fit <- exec(ncvreg::cv.ncvreg, !!!callArgs)
    coefs <- as.numeric(coef(fit, lambda = fit$lambda.min))
    weights <- replace(eff.wgt, keep, coefs[-1L])
    if (identical(fitRetention, "none")) {
        return(weights)
    }
    kept <- if (identical(fitRetention, "full")) {
        `[[<-`(fit, "fit", .withoutFields(fit$fit, "X"))
    } else {
        .ncvregSlimFit(fit, coefs)
    }
    `attr<-`(weights, "fit", .withFitVariantIds(kept, colnames(X)[keep]))
}

# The fit cv.ncvreg selected: its coefficients over the variants fitted, the
# intercept, its lambda and the cross-validated error there.
# @noRd
.ncvregSlimFit <- function(cvFit, coefs) {
    list(
        coef = coefs[-1L],
        intercept = coefs[[1L]],
        lambda = cvFit$lambda.min,
        cve = cvFit$cve[[cvFit$min]],
        cvse = cvFit$cvse[[cvFit$min]]
    )
}

# `fit` without its `fields`, keeping its class: an engine object's copies of
# the data it was given, which a stored fit has no use for.
# @noRd
.withoutFields <- function(fit, fields) {
    `class<-`(fit[setdiff(names(fit), fields)], class(fit))
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
#'   with \code{\link{NcvregOptions}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @inheritParams ncvregWeights
#' @return A numeric vector of length `ncol(X)` of variant weights, with the
#'   fit attached as attribute \code{"fit"} unless
#'   \code{fitRetention = "none"}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' scadWeights(X, y)
#' @export
scadWeights <- function(
    X,
    y,
    nfolds = 5,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NcvregOptions()
) {
    .assertMethodOptions(methodArgs, "NcvregOptions", "methodArgs")
    ncvregWeights(
        X,
        y,
        penalty = "SCAD",
        nfolds = nfolds,
        fitRetention = arg_match(fitRetention),
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
#'   with \code{\link{NcvregOptions}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @inheritParams ncvregWeights
#' @return A numeric vector of length `ncol(X)` of variant weights, with the
#'   fit attached as attribute \code{"fit"} unless
#'   \code{fitRetention = "none"}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' mcpWeights(X, y)
#' @export
mcpWeights <- function(
    X,
    y,
    nfolds = 5,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NcvregOptions()
) {
    .assertMethodOptions(methodArgs, "NcvregOptions", "methodArgs")
    ncvregWeights(
        X,
        y,
        penalty = "MCP",
        nfolds = nfolds,
        fitRetention = arg_match(fitRetention),
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
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights: \code{"none"} (default) attaches
#'   nothing; \code{"slim"} the selected fit, as a list of its coefficients
#'   (\code{coef}, over the variants fitted), \code{intercept},
#'   \code{gamma}, \code{lambda} and cross-validated error there
#'   (\code{cvMean}, \code{cvSd}); \code{"full"} the whole
#'   \code{L0Learn.cvfit} object, with every candidate's coefficients.
#'   Either names the variants it was fitted on, in \code{variantIds}.
#' @param methodArgs Options forwarded to \code{L0Learn::L0Learn.cvfit}, built
#'   with \code{\link{L0learnOptions}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of length `ncol(X)` of variant weights, with the
#'   fit attached as attribute \code{"fit"} unless
#'   \code{fitRetention = "none"}.
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
    fitRetention = c("none", "slim", "full"),
    methodArgs = L0learnOptions()
) {
    .assertMethodOptions(methodArgs, "L0learnOptions", "methodArgs")
    fitRetention <- arg_match(fitRetention)
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
        !!!as.list(methodArgs)
    )
    fit <- exec(L0Learn::L0Learn.cvfit, !!!callArgs)
    best <- .l0learnBest(fit)
    raw <- as.numeric(coef(fit, lambda = best$lambda, gamma = best$gamma))
    # If intercept was included, it is the first row.
    hasIntercept <- length(raw) == sum(keep) + 1L
    coefs <- if (hasIntercept) raw[-1L] else raw
    weights <- replace(eff.wgt, keep, coefs)
    if (identical(fitRetention, "none")) {
        return(weights)
    }
    kept <- if (identical(fitRetention, "full")) {
        fit
    } else {
        c(
            list(coef = coefs, intercept = if (hasIntercept) raw[[1L]] else 0),
            best
        )
    }
    `attr<-`(weights, "fit", .withFitVariantIds(kept, colnames(X)[keep]))
}

# The (gamma, lambda) pair minimizing the cross-validation error across the
# whole path, with that error and its standard deviation.
# @noRd
.l0learnBest <- function(fit) {
    gammaIdx <- which.min(map_dbl(fit$cvMeans, .rssMinNumeric))
    cvMeans <- as.numeric(fit$cvMeans[[gammaIdx]])
    lambdaIdx <- which.min(cvMeans)
    list(
        gamma = fit$fit$gamma[[gammaIdx]],
        lambda = fit$fit$lambda[[gammaIdx]][[lambdaIdx]],
        cvMean = cvMeans[[lambdaIdx]],
        cvSd = as.numeric(fit$cvSDs[[gammaIdx]])[[lambdaIdx]]
    )
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
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights: \code{"none"} (default) attaches
#'   nothing, \code{"slim"} and \code{"full"} the fitted \code{BGLR} object,
#'   which keeps no copy of the data, so the two retaining levels behave
#'   alike here.
#' @param methodArgs Options forwarded to \code{BGLR::BGLR}, built
#'   with \code{\link{BglrOptions}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of length `ncol(X)` of variant weights, with the
#'   fit attached as attribute \code{"fit"} unless
#'   \code{fitRetention = "none"}.
#' @keywords internal
bglrWeights <- function(
    X,
    y,
    model,
    nIter,
    burnIn,
    thin,
    etaArgs = list(),
    fitRetention = c("none", "slim", "full"),
    methodArgs = BglrOptions()
) {
    .assertMethodOptions(methodArgs, "BglrOptions", "methodArgs")
    fitRetention <- arg_match(fitRetention)
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
        !!!as.list(methodArgs)
    )
    fit <- exec(BGLR::BGLR, !!!callArgs)
    weights <- replace(eff.wgt, keep, as.numeric(fit$ETA[[1]]$b))
    if (identical(fitRetention, "none")) {
        return(weights)
    }
    `attr<-`(weights, "fit", .withFitVariantIds(fit, colnames(X)[keep]))
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
#'   with \code{\link{DprOptions}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights: \code{"none"} (default) attaches
#'   nothing, \code{"slim"} the fitted model without its copy of the
#'   genotypes (\code{x}), \code{"full"} the whole fitted model. Either
#'   names the variants it was fitted on.
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
    methodArgs = DprOptions()
) {
    fitRetention <- arg_match(fitRetention)
    .assertMethodOptions(methodArgs, "DprOptions", "methodArgs")
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
        !!!as.list(methodArgs)
    )
    fit <- exec(RcppDPR::fit_model, !!!callArgs)
    eff.wgt <- replace(zeros, keep, as.numeric(fit$beta + fit$alpha))
    if (identical(fitRetention, "none")) {
        return(eff.wgt)
    }
    kept <- if (identical(fitRetention, "full")) {
        fit
    } else {
        .withoutFields(fit, "x")
    }
    `attr<-`(eff.wgt, "fit", .withFitVariantIds(kept, colnames(X)[keep]))
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
    methodArgs = DprOptions()
) {
    fitRetention <- arg_match(fitRetention)
    .assertMethodOptions(methodArgs, "DprOptions", "methodArgs")
    dprWeights(
        X,
        y,
        fittingMethod = "VB",
        fitRetention = fitRetention,
        methodArgs = .methodOptionsWith(methodArgs, "DprOptions", n_k = nK)
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
    methodArgs = DprOptions()
) {
    .assertMethodOptions(methodArgs, "DprOptions", "methodArgs")
    assertCount(sStep, positive = TRUE)
    fitRetention <- arg_match(fitRetention)
    dprWeights(
        X,
        y,
        fittingMethod = "Gibbs",
        fitRetention = fitRetention,
        methodArgs = .methodOptionsWith(
            methodArgs,
            "DprOptions",
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
    methodArgs = DprOptions()
) {
    .assertMethodOptions(methodArgs, "DprOptions", "methodArgs")
    fitRetention <- arg_match(fitRetention)
    dprWeights(
        X,
        y,
        fittingMethod = "Adaptive_Gibbs",
        fitRetention = fitRetention,
        methodArgs = methodArgs
    )
}
# Starting coefficients as a plain matrix: .mrmashInitCoefficients()
# returns the fit it used, and only its Bhat is wanted here.
# @noRd
.mrmashInitB <- function(X, Y, bInitMethod, standardize, BPPARAM) {
    as.matrix(
        .mrmashInitCoefficients(
            X,
            Y,
            bInitMethod,
            standardize,
            BPPARAM
        )$Bhat
    )
}

# Everything mrmashWrapper needs in hand before the fit: the package and
# bundle checks, the resolved coefficient-initialisation method, the
# summary statistics the prior grid is built from (computed here when the
# caller did not supply them), and the prior covariance matrices.
# @noRd
.mrmashWrapperSetup <- function(
    X,
    Y,
    prior,
    methodArgs,
    dataDrivenPriorMatrices,
    bInitMethod,
    sumstats,
    standardize,
    threads
) {
    .mrmashRequirePackages()
    .assertMethodParam(prior, "MrmashPriorParam")
    .assertMethodOptions(methodArgs, "MrmashOptions", "methodArgs")
    priorOpts <- as.list(prior)
    .mrmashValidateWrapper(
        X = X,
        Y = Y,
        dataDrivenPriorMatrices = dataDrivenPriorMatrices,
        canonicalPriorMatrices = priorOpts$canonicalPriorMatrices %||% FALSE
    )
    if (is.null(sumstats)) {
        sumstats <- .mrmashComputeSumstats(X, Y, standardize, threads)
    }
    list(
        bInitMethod = .mrmashResolveBInit(Y, bInitMethod),
        # Shared prior-covariance builder (also used by mrmashRssWeights).
        priorBuilt = exec(
            buildMrmashPriorMatrices,
            Bhat = sumstats$Bhat,
            Shat = sumstats$Shat,
            K = ncol(Y),
            dataDrivenPriorMatrices = dataDrivenPriorMatrices,
            !!!priorOpts
        )
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
#'   \code{\link{MrmashPriorParam}}.
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
#' @param BPPARAM A \code{BiocParallelParam}. The univariate glmnet fits that
#'   initialise the coefficients run on its workers, and the summary
#'   statistics use as many cores as it has workers when they run on this
#'   machine (\code{SerialParam}, \code{MulticoreParam}); for other backends
#'   they use one. The mr.mash fit itself always runs on one thread (see
#'   \code{\link{MrmashOptions}}). Default \code{SerialParam()}.
#' @param methodArgs Options forwarded to \code{mr.mashr::mr.mash} under its
#'   own names, built with \code{\link{MrmashOptions}}
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
#'   prior = MrmashPriorParam(canonicalPriorMatrices = TRUE)
#' )
#'
#' @export
mrmashWrapper <- function(
    X,
    Y,
    V = NULL,
    sumstats = NULL,
    dataDrivenPriorMatrices = NULL,
    prior = MrmashPriorParam(),
    bInitMethod = "enet",
    updateV = TRUE,
    updateVMethod = "full",
    standardize = FALSE,
    BPPARAM = SerialParam(),
    methodArgs = MrmashOptions()
) {
    .assertBpparam(BPPARAM)
    threads <- .bpLocalThreads(BPPARAM, "mrmashWrapper")
    setup <- .mrmashWrapperSetup(
        X = X,
        Y = Y,
        prior = prior,
        methodArgs = methodArgs,
        dataDrivenPriorMatrices = dataDrivenPriorMatrices,
        bInitMethod = bInitMethod,
        sumstats = sumstats,
        standardize = standardize,
        threads = threads
    )
    bInitMethod <- setup$bInitMethod
    priorBuilt <- setup$priorBuilt
    time1 <- proc.time()
    bInit <- .mrmashInitB(X, Y, bInitMethod, standardize, BPPARAM)
    vInit <- .mrmashInitV(X, Y, V, updateV, updateVMethod)
    fitMrmash <- .mrmashFit(
        priorBuilt$S0,
        bInit,
        vInit,
        X = X,
        Y = Y,
        standardize = standardize,
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
.mrmashComputeSumstats <- function(X, Y, standardize, threads) {
    mr.mashr::compute_univariate_sumstats(
        X,
        Y,
        standardize = standardize,
        standardize.response = FALSE,
        mc.cores = threads
    )
}

# Initial coefficient matrix via graphical-lasso or univariate glmnet.
# @noRd
.mrmashInitCoefficients <- function(
    X,
    Y,
    bInitMethod,
    standardize,
    BPPARAM
) {
    if (bInitMethod == "glasso") {
        return(computeCoefficientsGlasso(
            X,
            Y,
            standardize = standardize,
            Xnew = NULL
        ))
    }
    computeCoefficientsUnivGlmnet(
        X,
        Y,
        alpha = 0.5,
        standardize = standardize,
        BPPARAM = BPPARAM,
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

# Run mr.mash with the resolved prior / init / V. Everything pecotmr derives --
# the prior, the residual covariance, the coefficient initialisation -- is
# supplied here and cannot be overridden; the rest of mr.mash's arguments come
# from `methodArgs` under mr.mash's own names, with pecotmr's defaults
# underneath.
# @noRd
.mrmashFit <- function(
    S0,
    bInit,
    vInit,
    X,
    Y,
    standardize,
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
        # One thread: see .mrmashUnavailable().
        nthreads = 1L
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
#' @param Xnew Optional new design matrix for prediction; \code{NULL} returns
#'   in-sample coefficients.
#' @return A numeric matrix of initial coefficient estimates (variants x
#'   conditions).
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' Y <- matrix(rnorm(nrow(X) * 3), nrow(X), 3)
#' computeCoefficientsGlasso(X = X, Y = Y, standardize = TRUE, Xnew = NULL)
#' @importFrom checkmate assertFlag
#' @export
computeCoefficientsGlasso <- function(X, Y, standardize, Xnew = NULL) {
    assertFlag(standardize)
    p <- ncol(X)
    r <- ncol(Y)
    conditionNames <- colnames(Y)

    # Fit group-lasso. `parallel = FALSE` used to be restated here; it is
    # cv.glmnet's own default, so transcribing it bought nothing.
    cvfitGlmnet <- fitGlmnetCv(
        x = X,
        y = Y,
        alpha = 1,
        family = "mgaussian",
        standardize = standardize
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

    cvfit <- fitGlmnetCv(
        x = Xnomiss,
        y = Ynomiss,
        alpha = alpha,
        family = "gaussian",
        standardize = standardize
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

# Initial mr.mash coefficients from one cv.glmnet per response. The fits are
# independent, so they run on BPPARAM's workers, each on its own seeded
# random-number stream (cv.glmnet draws its folds at random). The data go by
# position: a named `X` would be taken for bplapply()'s own X.
# @noRd
computeCoefficientsUnivGlmnet <- function(
    X,
    Y,
    alpha,
    standardize,
    BPPARAM = SerialParam(),
    Xnew = NULL
) {
    out <- .bpMap(
        seq_len(ncol(Y)),
        .linreg,
        X,
        Y,
        alpha,
        standardize,
        Xnew,
        BPPARAM = BPPARAM
    )

    bhatList <- map(out, "bhat")
    Bhat <- exec(cbind, !!!bhatList)

    if (!is.null(Xnew)) {
        yhatList <- map(out, "yhat_new")
        YhatNew <- `colnames<-`(exec(cbind, !!!yhatList), colnames(Y))
        results <- list(
            Bhat = Bhat[-1, , drop = FALSE],
            intercept = Bhat[1, ],
            Yhat_new = YhatNew
        )
    } else {
        results <- list(Bhat = Bhat[-1, , drop = FALSE], intercept = Bhat[1, ])
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
# `mult = sqrt(2)` and the square follow mr.mashr's own documented recipe
# (the mr.mash() and mr.mash.rss() examples build their grids exactly so):
# this grid scales variances, not standard deviations. autoselect.mixsd()
# drops Shat that is zero or non-finite, which covers NA, so the port's extra
# is.na() term was redundant.
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
    priorNames <- names(formals(MrmashPriorParam))
    isPrior <- is_in(names(args), priorNames)
    list(prior = args[isPrior], rest = args[!isPrior])
}

#' @title Arguments For mr.mashr's Canonical-Covariance Builder
#' @description Arguments for \code{mr.mashr::compute_canonical_covs()}, the
#'   function that builds the canonical prior mixture. \code{r} (the number
#'   of conditions) is pecotmr's to supply, so it is refused here; every
#'   other name is checked against that function's live formals.
#'
#'   No defaults: pecotmr used to carry \code{singletons = TRUE} and
#'   \code{hetgrid = c(0, 0.25, 0.5, 0.75, 1)}, which are character for
#'   character the engine's own defaults. Transcribing them bought nothing
#'   and would silently go stale if mr.mashr changed them.
#' @param ... Arguments for \code{mr.mashr::compute_canonical_covs()}, such
#'   as \code{singletons} and \code{hetgrid}.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' MrmashCanonicalCovsOptions(singletons = FALSE)
#' @export
MrmashCanonicalCovsOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(r = "the condition count, taken from the data"),
        "MrmashCanonicalCovsOptions"
    )
    .newMethodOptions(
        "mr.mashr::compute_canonical_covs",
        defaults = list(),
        extra = extra,
        label = "MrmashCanonicalCovsOptions",
        engine = "mrmashCanonicalCovs"
    )
}

#' @title Arguments For mr.mashr's Covariance Expansion
#' @description Arguments for \code{mr.mashr::expand_covs()}, which expands
#'   the prior covariances over the scaling grid. The covariance list and the
#'   grid are both pecotmr's to supply, so they are refused here; every other
#'   name is checked against that function's live formals.
#' @param ... Arguments for \code{mr.mashr::expand_covs()}, such as
#'   \code{zeromat}.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' MrmashExpandCovsOptions(zeromat = FALSE)
#' @export
MrmashExpandCovsOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            mats = "the covariances assembled from the canonical and ",
            grid = "`priorGrid` on MrmashPriorParam()"
        ),
        "MrmashExpandCovsOptions"
    )
    .newMethodOptions(
        "mr.mashr::expand_covs",
        defaults = list(),
        extra = extra,
        label = "MrmashExpandCovsOptions",
        engine = "mrmashExpandCovs"
    )
}

#' @rdname MrmashPriorParam
#' @aliases MrmashPriorParam-class
#' @exportClass MrmashPriorParam
setClass(
    "MrmashPriorParam",
    contains = "MethodParam",
    slots = c(
        canonicalPriorMatrices = "logical",
        priorGrid = "numeric_OR_NULL",
        canonicalCovs = "MethodOptions",
        expandCovs = "MethodOptions"
    )
)

#' @title Prior-Covariance Arguments For mr.mash
#' @description Settings for \code{\link{buildMrmashPriorMatrices}}, the
#'   shared prior-covariance builder behind \code{\link{mrmashWrapper}} and
#'   \code{\link{mrmashRssWeights}}.
#'
#'   The builder reaches two mr.mashr functions, so the arguments bound for
#'   each travel in their own nested bundle --- \code{canonicalCovs} for
#'   \code{compute_canonical_covs()} and \code{expandCovs} for
#'   \code{expand_covs()} --- and each is checked against that function's
#'   live formals. The two fields pecotmr itself reads stay flat.
#'
#'   \code{dataDrivenPriorMatrices} is \emph{not} set here: it is an input
#'   computed upstream (typically by \code{\link{mashPriorCovariances}}),
#'   not a tuning choice, so it stays a parameter of the caller.
#' @param canonicalPriorMatrices Logical. Include the canonical mixture from
#'   \code{mr.mashr::compute_canonical_covs()}. Default \code{FALSE}.
#' @param priorGrid Optional pre-computed scaling grid (numeric vector).
#'   \code{NULL} (default) derives one from the summary statistics via
#'   \code{computeGrid()}.
#' @param canonicalCovs Arguments for
#'   \code{mr.mashr::compute_canonical_covs()}, built with
#'   \code{\link{MrmashCanonicalCovsOptions}}. Only consulted when
#'   \code{canonicalPriorMatrices} is \code{TRUE}.
#' @param expandCovs Arguments for \code{mr.mashr::expand_covs()}, built
#'   with \code{\link{MrmashExpandCovsOptions}}.
#' @return A \code{MrmashPriorParam} object, a \code{\link{MethodParam}}.
#' @examples
#' MrmashPriorParam(canonicalPriorMatrices = TRUE)
#' MrmashPriorParam(
#'     canonicalCovs = MrmashCanonicalCovsOptions(singletons = FALSE)
#' )
#' @export
MrmashPriorParam <- function(
    canonicalPriorMatrices = FALSE,
    priorGrid = NULL,
    canonicalCovs = MrmashCanonicalCovsOptions(),
    expandCovs = MrmashExpandCovsOptions()
) {
    .assertMethodOptions(
        canonicalCovs,
        "MrmashCanonicalCovsOptions",
        "canonicalCovs"
    )
    .assertMethodOptions(expandCovs, "MrmashExpandCovsOptions", "expandCovs")
    new(
        "MrmashPriorParam",
        canonicalPriorMatrices = canonicalPriorMatrices,
        priorGrid = priorGrid,
        canonicalCovs = canonicalCovs,
        expandCovs = expandCovs
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
#' @param canonicalCovs Arguments for
#'   \code{mr.mashr::compute_canonical_covs()} --- \code{singletons},
#'   \code{hetgrid} --- built with
#'   \code{\link{MrmashCanonicalCovsOptions}}. Empty by default, which
#'   leaves that function's own defaults in place.
#' @param expandCovs Arguments for \code{mr.mashr::expand_covs()} ---
#'   \code{zeromat} --- built with \code{\link{MrmashExpandCovsOptions}}.
#'   Empty by default, which leaves that function's own defaults in place.
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
    canonicalCovs = MrmashCanonicalCovsOptions(),
    expandCovs = MrmashExpandCovsOptions()
) {
    if (!requireNamespace("mr.mashr", quietly = TRUE)) {
        abort("Package 'mr.mashr' is required.")
    }
    .assertMethodOptions(
        canonicalCovs,
        "MrmashCanonicalCovsOptions",
        "canonicalCovs"
    )
    .assertMethodOptions(expandCovs, "MrmashExpandCovsOptions", "expandCovs")
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
    S0_raw <- .mrmashRawPriorCovs(
        K,
        dataDrivenPriorMatrices,
        canonicalPriorMatrices,
        canonicalCovs
    )
    S0 <- exec(
        mr.mashr::expand_covs,
        S0_raw,
        priorGrid,
        !!!as.list(expandCovs)
    )
    list(S0 = S0, priorGrid = priorGrid)
}

# The unexpanded prior covariances: the canonical set (when requested) plus
# any data-driven matrices, or the data-driven matrices alone.
# @noRd
.mrmashRawPriorCovs <- function(
    K,
    dataDrivenPriorMatrices,
    canonicalPriorMatrices,
    canonicalCovs
) {
    if (!isTRUE(canonicalPriorMatrices)) {
        return(dataDrivenPriorMatrices$U)
    }
    canonical <- exec(
        mr.mashr::compute_canonical_covs,
        K,
        !!!as.list(canonicalCovs)
    )
    if (is.null(dataDrivenPriorMatrices)) {
        return(canonical)
    }
    c(canonical, dataDrivenPriorMatrices$U)
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# `.rssSolvePath` solver for the penalized RSS path (lasso and its elastic-net
# mixtures, MCP, SCAD, L0...) by coordinate descent, on the LD matrix or on a
# factor of it (never forming R); a ridge in its closed form
# (.rssRidgeSolve()).
# @noRd
.rssPenalizedSolve <- function(
    z,
    lam,
    ld,
    penalty,
    gamma,
    alpha,
    lambda0,
    lambda2,
    thr,
    maxiter,
    maxSwaps
) {
    if (.rssIsRidge(penalty, alpha)) {
        return(.rssRidgeSolve(z, lam, ld))
    }
    settings <- list(
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
    if (!is.null(ld$F)) {
        return(exec(
            penalizedRssFactorCpp,
            zR = z,
            factorR = ld$F,
            shrink = ld$shrink %||% 0,
            !!!settings
        ))
    }
    exec(
        penalizedRssCpp,
        zR = z,
        LD = list(blk1 = ld$R),
        shrink = ld$shrink %||% 0,
        !!!settings
    )
}

# Whether a penalty is a ridge: the lasso's at alpha = 0, its L2 term alone.
# @noRd
.rssIsRidge <- function(penalty, alpha) {
    identical(penalty, "lasso") && alpha == 0
}

# The LD source as the penalty's solver reads it. A ridge is solved from R's
# eigendecomposition (.rssRidgeSolve()), which the source then carries,
# computed here once for every shrinkage level, lambda and held-out refit of
# the grid; any other penalty reads the source as it is.
# @noRd
.rssSolverSource <- function(ld, spec) {
    if (!.rssIsRidge(spec$penalty, spec$alpha) || !is.null(ld$eigen)) {
        return(ld)
    }
    list_assign(ld, eigen = .ldSourceEigen(ld))
}

# The ridge's path in closed form. Its objective beta'R_s beta - 2 beta'z +
# lambda beta'beta, for R_s = (1 - s) R + s I, has the minimizer
# (R_s + lambda I)^-1 z, read for every lambda from R's eigendecomposition
# R = V D V' (.ldSourceEigen(); a factor's thin SVD leaves V narrower than
# R, which is zero off V's span): with t = s + lambda,
# beta = V (V'z / ((1 - s) D + t)) + (z - V V'z) / t. Where (1 - s) D + t
# vanishes (s = lambda = 0 on an LD of less than full rank) beta is the
# minimum-norm minimizer. Returns what penalizedRssCpp() returns, each fit
# exact.
# @noRd
.rssRidgeSolve <- function(z, lam, ld) {
    ev <- .ldSourceEigen(ld)
    s <- ld$shrink %||% 0
    proj <- as.numeric(crossprod(ev$vectors, z))
    rest <- z - as.numeric(ev$vectors %*% proj)
    top <- (1 - s) * max(ev$values, 0) + s + lam
    coef <- proj *
        .rssInverseAbove(outer((1 - s) * ev$values, s + lam, "+"), top)
    restCoef <- as.numeric(.rssInverseAbove(matrix(s + lam, 1L), top))
    restSq <- sum(rest^2)
    bb <- colSums(coef^2) + restSq * restCoef^2
    zb <- colSums(proj * coef) + restSq * restCoef
    loss <- (1 - s) * colSums(ev$values * coef^2) + s * bb - 2 * zb
    list(
        beta = ev$vectors %*% coef + outer(rest, restCoef),
        lambda = lam,
        conv = rep(1L, length(lam)),
        loss = loss,
        fbeta = loss + lam * bb
    )
}

# 1 / x for a spectrum `x`, one column per lambda, and 0 where an entry is
# zero to rounding: at most sqrt(machine epsilon) of its column's largest
# eigenvalue `top`, as MASS::ginv() cuts a pseudo-inverse.
# @noRd
.rssInverseAbove <- function(x, top) {
    kept <- x > rep(sqrt(.Machine$double.eps) * top, each = nrow(x))
    replace(array(0, dim(x)), kept, 1 / x[kept])
}

# One L0Learn fit at lambda0 = `l0Val`: penalizedRss over the lambda path.
# @noRd
.rssL0LambdaFit <- function(l0Val, solverInput, ld, n, sVal, config) {
    model <- .penalizedRssOn(
        solverInput,
        ld,
        n,
        list_assign(config$spec, lambda0 = l0Val),
        config$options
    )
    list(
        beta = model$beta,
        meta = mutate(.rssPathMeta(sVal, model), lambda0 = l0Val, .after = "s")
    )
}

# The minimum of one CV-mean path vector (coerced numeric).
# @noRd
.rssMinNumeric <- function(v) {
    min(as.numeric(v))
}
