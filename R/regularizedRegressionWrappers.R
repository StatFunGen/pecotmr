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
#' @param sigma2E The residual variance, or \code{NULL} (default) for
#'   \code{varY} when it is finite, else 1 (the standardized scale). It is
#'   held fixed, as the other summary-statistics fits hold theirs: against a
#'   reference LD, an estimate absorbs LD mismatch.
#'   \code{MrashRssOptions(update_sigma = TRUE)} estimates it from this
#'   start, as \code{mr.ash.rss} does by default.
#' @param s0 Prior variances of the mixture components, relative to the
#'   residual variance, or \code{NULL} (default) for the grid
#'   \code{susieR::mr.ash} uses by default, built for these statistics.
#' @param w0 Prior mixture weights (summing to 1), or \code{NULL} (default)
#'   for equal weights.
#' @param updateOrder The order the coordinate ascent updates the variants
#'   in: \code{"decreasingLogBf"} (default), by decreasing marginal log Bayes
#'   factor under the prior \code{mr.ash.rss} fits (\code{s0}, \code{w0}),
#'   or \code{"consecutive"}, their order in \code{stat}. The engine has no
#'   order of its own, so the variants are permuted before the fit, the LD
#'   formed in that order, and the results put back after. On the
#'   standardized scale (no \code{varY} with \code{seb}), or with
#'   \code{standardize = TRUE}, the Bayes factors order the variants exactly
#'   as \code{|z|} does. On real genotypes the default predicted slightly
#'   better (explained variance about +0.005 over the column order) at no
#'   cost in time.
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
    updateOrder = c("decreasingLogBf", "consecutive"),
    fitRetention = c("none", "slim", "full"),
    methodArgs = MrashRssOptions()
) {
    retention <- arg_match(fitRetention)
    ordering <- arg_match(updateOrder)
    .assertMethodOptions(methodArgs, "MrashRssOptions", "methodArgs")
    assertList(stat)
    .rssStatVariantIds(stat, "mrashRssWeights")
    n <- median(stat$n)
    inp <- list(
        b = as.numeric(stat$b %||% numeric(0)),
        shat = as.numeric(stat$seb %||% numeric(0)),
        z = .mrashRssZ(stat),
        n = n,
        varY = varY,
        sigma2E = sigma2E
    )
    prior <- .mrashRssPrior(s0, w0, inp)
    perm <- .mrashRssOrder(ordering, inp, prior, methodArgs)
    R <- .rssLdMatrix(stat, ldReference, "mrashRssWeights", perm)
    callArgs <- .mrashRssCallArgs(inp, R, prior, methodArgs, perm)
    model <- .mrashRssUnpermute(exec(mr.ash.rss, !!!callArgs), perm)
    if (identical(retention, "none")) {
        return(model$mu1)
    }
    `attr<-`(
        model$mu1,
        "fit",
        .withFitVariantIds(model, stat[["variantNames"]])
    )
}

# The prior mr.ash.rss() fits: the caller's grid `s0`, else susieR::mr.ash()'s
# default built for these statistics, and the caller's weights `w0`, else
# equal ones.
# @noRd
.mrashRssPrior <- function(s0, w0, inp) {
    grid <- s0 %||% .mrashRssDefaultGrid(inp$n, inp$z, inp$shat, inp$varY)
    list(s0 = grid, w0 = w0 %||% rep(1 / length(grid), length(grid)))
}

# The order the variants are given to mr.ash.rss() in, which it has none of
# its own: NULL for "consecutive", their own, else by decreasing marginal
# log Bayes factor under the prior it will fit (.mrashRssLogBf()), ties in
# their own order.
# @noRd
.mrashRssOrder <- function(ordering, inp, prior, methodArgs) {
    if (identical(ordering, "consecutive")) {
        return(NULL)
    }
    standardize <- isTRUE(as.list(methodArgs)$standardize)
    order(-.mrashRssLogBf(inp, prior, standardize))
}

# Each variant's marginal log Bayes factor under the prior mr.ash.rss()
# fits, from the sufficient statistics it builds (susieR src/mr_ash_rss.h):
# with adj = (n - 1) / (z^2 + n - 2), x'x = varY adj / seb^2 and x'y = z adj
# varY / seb for a finite varY and seb, else x'x = n - 1 and x'y = z
# sqrt(adj (n - 1)), the standardized scale; x'x = n - 1 again when it
# standardizes. The residual variance is sigma2E, else varY when finite,
# else 1. Where x'x is the same for every variant -- the standardized scale,
# or standardize = TRUE -- the order is that of |z|.
# @noRd
.mrashRssLogBf <- function(inp, prior, standardize) {
    adj <- (inp$n - 1) / (inp$z^2 + inp$n - 2)
    finiteVarY <- !is.null(inp$varY) && is.finite(inp$varY)
    unstandardized <- finiteVarY && length(inp$shat) > 0L
    xtx <- if (unstandardized) {
        inp$varY * adj / inp$shat^2
    } else {
        rep(inp$n - 1, length(inp$z))
    }
    xty <- if (unstandardized) {
        inp$z * adj * inp$varY / inp$shat
    } else {
        inp$z * sqrt(adj * (inp$n - 1))
    }
    sx <- if (standardize) sqrt(xtx / (inp$n - 1)) else 1
    sigma2 <- inp$sigma2E %||% (if (finiteVarY) inp$varY else 1)
    .mixtureLogBf(xtx / sx^2, xty / sx, sigma2, prior$s0, prior$w0)
}

# mr.ash.rss()'s arguments, in the variants' order `perm` when one is given:
# the statistics and a caller's mu1_init permuted (an empty bhat or seb left
# empty), R already formed in that order. The caller's options come last.
# @noRd
.mrashRssCallArgs <- function(inp, R, prior, methodArgs, perm) {
    idx <- perm %||% seq_along(inp$z)
    user <- .mrashRssPermuteUser(as.list(methodArgs), perm)
    list_modify(
        list(
            bhat = .rowsAt(inp$b, idx),
            shat = .rowsAt(inp$shat, idx),
            z = inp$z[idx],
            R = R,
            var_y = inp$varY,
            n = inp$n,
            sigma2_e = inp$sigma2E,
            # Held fixed, as the other RSS fits hold the residual variance:
            # against a reference LD, an estimate absorbs LD mismatch.
            update_sigma = FALSE,
            s0 = prior$s0,
            # A copy: mr.ash.rss() writes its fitted weights over the vector
            # it is given, which would otherwise be the caller's.
            w0 = prior$w0 + 0
        ),
        !!!user
    )
}

# The caller's options, with a mu1_init they give put in the order `perm`.
# @noRd
.mrashRssPermuteUser <- function(user, perm) {
    if (is.null(perm) || length(user$mu1_init) == 0L) {
        return(user)
    }
    list_assign(user, mu1_init = .rowsAt(user$mu1_init, perm))
}

# The fit's per-variant results put back in the statistics' order, from the
# order `perm` the variants were fitted in: mu1, beta and the rows of
# sigma2_1 and w1, as many as the fit has.
# @noRd
.mrashRssUnpermute <- function(model, perm) {
    if (is.null(perm)) {
        return(model)
    }
    fields <- intersect(c("mu1", "beta", "sigma2_1", "w1"), names(model))
    restored <- set_names(
        map(model[fields], .rowsAt, idx = order(perm)),
        fields
    )
    list_assign(model, !!!restored)
}

# `x` at `idx`, by row for a matrix; an empty `x` stays empty, where
# indexing it would give NAs.
# @noRd
.rowsAt <- function(x, idx) {
    if (length(x) == 0L) {
        return(x)
    }
    if (is.matrix(x)) x[idx, , drop = FALSE] else x[idx]
}

# The marginal log Bayes factor of each variant under mr.ash's mixture
# prior, for x'x = `xtx`, x'y = `xty` and residual variance `sigma2`: the
# effect is sum_k w_k N(0, sigma2 s_k), with s_k = `priorScale` and w_k =
# `priorWeights`, so component k's log Bayes factor is -log(1 + kappa) / 2
# + chi2 kappa / (2 (1 + kappa)), with kappa = s_k x'x and chi2 = (x'y)^2 /
# (x'x sigma2) (susieR's bayes_ridge_sufficient()), and the variant's is
# their log-sum-exp under the weights.
# @noRd
#' @importFrom matrixStats rowLogSumExps
.mixtureLogBf <- function(xtx, xty, sigma2, priorScale, priorWeights) {
    kappa <- outer(xtx, priorScale)
    chi2 <- xty^2 / (xtx * sigma2)
    logBf <- -0.5 * log1p(kappa) + 0.5 * chi2 * kappa / (1 + kappa)
    rowLogSumExps(sweep(logBf, 2, log(priorWeights), "+"))
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
#' @param rho Across conditions only (\code{\link{twasWeightsPipeline}}'s
#'   fits of a group on summary statistics, as JointPRS fits ancestries):
#'   the correlation of the conditions' effects. \code{NULL} (the default)
#'   learns each pair's on \eqn{[0, 0.99]} by JointPRS's Metropolis step; a
#'   number in \eqn{[0, 1)} fixes every pair at it; a \eqn{K \times K}
#'   positive-definite correlation matrix fixes each pair at its entry. One
#'   condition's fit has no pairs and ignores it.
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
    nChains = 1,
    rho = NULL
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
    if (is.null(ld$ldFactor)) {
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
        factorR = ld$ldFactor,
        seed = seed,
        track = track
    )
}

#' Extract weights from prsCs function
#' @param stat A list of summary statistics: \code{n} (per-variant sample
#'   sizes) and the marginal effects as \code{cor} (standardized), \code{z},
#'   \code{b} with \code{seb}, or \code{b} alone, taken as already
#'   standardized. The engine receives the correlation each z stands for,
#'   z / sqrt(z^2 + n - 2), with n the median of \code{n}.
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
    user <- as.list(methodArgs)
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
        !!!user
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
#' @param sharing Across studies only: the \code{dprGibbs} method fitted
#'   across the studies a \code{"study"} \code{jointSpecification} joins
#'   (\code{coupling} set in \code{\link{twasWeightsPipeline}}), as
#'   SDPRX fits ancestries. Which sets of the studies a non-null
#'   component's effects lie in: \code{"eachOrAll"} (the default), each
#'   study alone or all of them together; \code{"everySubset"}, every
#'   nonempty set, \eqn{2^K - 1} of them for \eqn{K} studies; \code{"all"},
#'   all of them together only, SDPRX's \code{force_shared}. The first two
#'   agree for two studies. A single study ignores it.
#' @param rho Across studies only: the correlation between a variant's
#'   effects in the studies a component spans. \code{NULL} (the default)
#'   learns each pair's on \eqn{[0, 0.99]}; a number in \eqn{[0, 1)} fixes
#'   every pair's, and a \eqn{K \times K} positive-definite correlation
#'   matrix each pair's. A single study ignores it.
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
#'   Across studies the same chain is SDPRX's, described in Zhou G, Chen T,
#'   Zhao H (2023), "SDPRX: A statistical method for cross-population
#'   prediction of complex traits", American Journal of Human Genetics
#'   110(1): 13-22. \doi{10.1016/j.ajhg.2022.11.007}. Each non-null
#'   component draws a variant's effects in the studies of its kind from
#'   \eqn{N(0, \sigma^2 \Sigma)}, \eqn{\Sigma} holding \code{rho}, for any
#'   number of studies. So that one study is SDPR draw for draw, the null
#'   stays the first stick of the process, the other kinds splitting the
#'   rest by a Dirichlet, and one concentration serves every kind; SDPRX
#'   puts the null among its kinds and gives each kind its own. A variant a
#'   study lacks can take any kind, its effect there integrated out, where
#'   SDPRX allows it only the kinds of the studies carrying it.
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
    nChains = 1,
    sharing = c("eachOrAll", "everySubset", "all"),
    rho = NULL
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
    if (is.null(ld$ldFactor)) {
        return(exec(
            sdprCpp,
            !!!settings,
            LD = list(blk1 = ld$R),
            seed = seed,
            track = track
        ))
    }
    exec(
        sdprFactorCpp,
        !!!settings,
        factorR = ld$ldFactor,
        seed = seed,
        track = track
    )
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
#'   standardized. The engine receives the correlation each z stands for,
#'   z / sqrt(z^2 + n - 2), with n the median of \code{n}.
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
    given <- as.list(methodArgs)
    # `sharing` and `rho` act across studies only, so one study drops them.
    user <- given[setdiff(names(given), .SDPR_ACROSS)]
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
        !!!user
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
#'   covariance matrices (\code{U}), and of the mash prior's residual
#'   correlation (\code{residualCorrelation}); forwarded to
#'   \code{mrmashWrapper} when it has to fit, and retained in the payload for
#'   mvSuSiE prior reconstruction.
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
        priorArgs <- split$prior
        optionArgs <- split$rest
        rows <- .mrmashRows(Y)
        mrmashFit <- exec(
            mrmashWrapper,
            X[rows, , drop = FALSE],
            Y[rows, , drop = FALSE],
            dataDrivenPriorMatrices = dataDrivenPriorMatrices,
            prior = exec(MrmashPriorParam, !!!priorArgs),
            methodArgs = exec(MrmashOptions, !!!optionArgs)
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


# The samples mr.mash is fitted on: every one, where the installed mr.mashr
# fits each condition on its own (`missing_y_method`), else those observing
# every condition.
# @noRd
.mrmashRows <- function(Y) {
    if (!anyNA(Y) || .mrmashTakesMissingY()) {
        return(seq_len(nrow(Y)))
    }
    .completeRows(
        Y,
        "mrmashWeights",
        "a mr.mashr with `missing_y_method`"
    )
}

# Whether the installed mr.mashr fits a condition on the samples it is
# observed in (its `missing_y_method`).
# @noRd
.mrmashTakesMissingY <- function() {
    is_in("missing_y_method", names(formals(mr.mashr::mr.mash)))
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

# mr.ash's `update.order` for the wrapper's `updateOrder`: none for
# "consecutive", which leaves the engine its column order, else the columns
# by decreasing marginal log Bayes factor under the prior it will fit
# (.mrashLogBf()), ties in column order.
# @noRd
.mrashUpdateOrder <- function(ordering, X, y, methodArgs) {
    if (identical(ordering, "consecutive")) {
        return(list())
    }
    list(update.order = order(-.mrashLogBf(X, y, as.list(methodArgs))))
}

# Each column's marginal log Bayes factor under the prior susieR::mr.ash()
# fits (R/mr.ash.R): X and y centred for its intercept and X scaled when it
# standardizes; the caller's grid `sa2`, else its default (2^((k - 1) / 25)
# - 1)^2 n / median(w) with w the columns' sums of squares; the caller's
# weights `pi`, else equal ones; the caller's residual variance `sigma2`,
# else the variance of y.
# @noRd
#' @importFrom stats var
.mrashLogBf <- function(X, y, settings) {
    intercept <- settings$intercept %||% TRUE
    Xc <- scale(X, center = intercept, scale = settings$standardize %||% FALSE)
    yc <- if (intercept) y - mean(y) else y
    w <- colSums(Xc^2)
    sa2 <- settings$sa2 %||% ((2^((0:24) / 25) - 1)^2 / median(w) * nrow(X))
    .mixtureLogBf(
        w,
        as.numeric(crossprod(Xc, yc)),
        settings$sigma2 %||% var(y),
        sa2,
        settings$pi %||% rep(1 / length(sa2), length(sa2))
    )
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
#'   May also include \code{Bhat}, \code{Shat} matrices, which with
#'   \code{covY} put the fit on their scale (without \code{covY} it is on
#'   the standardized scale of \code{z}: effects per standard deviation of
#'   the variant and of the condition), and
#'   \code{sampleGroups}, one label per column of \code{z} naming the
#'   samples it is measured on (such as its study): conditions with
#'   different labels are fitted as different samples even when they share
#'   an LD reference; and \code{h2}, one heritability per column of
#'   \code{z} (\code{NA} for none), which sets each condition's residual
#'   variance to \eqn{1 - h^2} of its variance, as SuSiE-RSS's residual
#'   variance, with the correlations of \code{covY} (a given \code{V}
#'   wins).
#' @param ldReference The LD reference the statistics' variants are matched
#'   to by \code{stat$variantNames}: a genotype panel (see
#'   \code{\link{readGenotypes}}), a \code{\link{GenotypeMixture}} or an
#'   \code{\link{LdMatrix}}. \code{mr.mash.rss} is given the LD over these
#'   variants: a factor of it for a genotype panel or a
#'   \code{GenotypeMixture} (its \code{X}), so the variants x variants
#'   matrix is never formed, and the matrix of an \code{LdMatrix} (its
#'   \code{R}). It is not read when
#'   \code{mrmashRssFit} is supplied. Conditions measured on different
#'   samples, such as studies of different ancestries, take a list of one
#'   reference per column of \code{stat$z}, the same reference for the
#'   conditions of one study. Each reference's LD is then taken over the
#'   variants its conditions carry (a non-missing z-score in any of them)
#'   that it holds, a variant they lack carrying no information in them,
#'   as factors when every reference has one and as matrices otherwise;
#'   the conditions on different references, or with different
#'   \code{stat$sampleGroups}, are fitted as different samples, with a
#'   residual covariance \code{V} and \code{covY} of zero between them, and
#'   the conditions of one sample share the median of their sample sizes.
#' @param mrmashRssFit Optional pre-fitted \code{mr.mash.rss} object; skips
#'   fitting when supplied.
#' @param dataDrivenPriorMatrices Optional list with element \code{U} (list of
#'   raw covariance matrices). Passed to \code{mr.mashr::expand_covs()}
#'   alongside the canonical mixture. It may also carry
#'   \code{residualCorrelation}, the residual correlation of the mash prior
#'   \code{U} comes from (conditions x conditions), the default \code{V} and
#'   \code{covY}. Matrices named by condition are taken on the columns of
#'   \code{stat$z} by name, so a prior learned on more conditions can be
#'   used (see \code{\link{buildMrmashPriorMatrices}}).
#' @param canonicalPriorMatrices Logical. When TRUE (default), include the
#'   standard canonical mixture from \code{mr.mashr::compute_canonical_covs()}.
#'   When FALSE, \code{dataDrivenPriorMatrices} must be supplied.
#' @param S0 Optional pre-built list of prior covariance matrices, bypassing the
#'   canonical / data-driven construction.
#' @param w0 Optional prior mixture weights; defaults to \code{computeW0(Bhat,
#'   length(S0))}.
#' @param V Optional residual covariance matrix (K x K). When NULL, the
#'   residual variances \code{stat$h2} implies, or without it
#'   \code{mr.mash.rss}'s default, \code{covY} (its correlations on the
#'   standardized scale).
#' @param covY Optional response covariance matrix (K x K). When NULL, the
#'   residual correlation of the mash prior, \code{dataDrivenPriorMatrices$
#'   residualCorrelation} (\code{\link{mashPipeline}}), zero between
#'   conditions on different references, which are fitted as different
#'   samples; without one, \code{mr.mash.rss} estimates it from the
#'   z-scores of the variants near-null (\eqn{|z| < 2}) in every condition
#'   on one reference, over the variants their reference holds, by mashr's
#'   \code{estimate_null_correlation_simple()} (fewer than
#'   \eqn{\max(2K_g, 20)} such variants, for the \eqn{K_g} conditions on a
#'   reference, give the identity there, with a warning).
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
#' # Three of the seven conditions the prior is over: it is taken on them by
#' # name.
#' cond <- colnames(multiTraitData$Y)[1:3]
#' zmat <- matrix(rnorm(30 * 3), 30, 3,
#'   dimnames = list(colnames(X), cond))
#' stat <- list(z = zmat, n = nrow(X), variantNames = colnames(X))
#' w <- mrmashRssWeights(stat, LdMatrix(cor(X)),
#'   dataDrivenPriorMatrices = multiTraitData$priorMatrices,
#'   canonicalPriorMatrices = TRUE,
#'   methodArgs = MrmashRssOptions(verbose = FALSE))
#' head(w)
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
            .rssConditionLd(stat, ldReference, "mrmashRssWeights"),
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

# Derive the (Z, K, nVec, Bhat, Shat) inputs from a multi-context stat list:
# Bhat and Shat, unless the statistics carry them, from z and the sample
# size of each variant in each condition (.mrmashRssSampleSizes()), and
# nVec, each condition's median sample size.
# @noRd
#' @importFrom matrixStats colMedians
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
    nMat <- .mrmashRssSampleSizes(stat$n, nrow(Z), K)
    list(
        Z = Z,
        K = K,
        nVec = colMedians(nMat),
        Bhat = stat$Bhat %||% (Z / sqrt(nMat)),
        Shat = stat$Shat %||% (1 / sqrt(nMat))
    )
}

# The statistics' sample sizes (`n`) as a variants x conditions matrix, from
# one for all, one per condition, one per variant (the same in every
# condition), or that matrix. A vector as long as both the conditions and
# the variants is read as one per condition.
# @noRd
.mrmashRssSampleSizes <- function(n, p, K) {
    if (is.matrix(n)) {
        assertMatrix(n, nrows = p, ncols = K, .var.name = "stat$n")
        return(n)
    }
    if (length(n) %in% c(1L, K)) {
        return(matrix(n, p, K, byrow = TRUE))
    }
    if (length(n) == p) {
        return(matrix(n, p, K))
    }
    abort(glue(
        "mrmashRssWeights: stat$n must be one sample size, one per ",
        "condition ({K}), one per variant ({p}), or a {p} x {K} matrix."
    ))
}

# Fit mr.mash.rss, defaulting the prior / w0 / V / covY when not supplied,
# with the data-driven prior taken on the statistics' conditions
# (.mrmashAlignPriorToContexts()).
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
    ss <- .mrmashRssStats(.rssCarriedStat(stat, LD))
    mashPrior <- .mrmashAlignPriorToContexts(
        dataDrivenPriorMatrices,
        colnames(ss$Z),
        ss$K
    )
    residual <- .mrmashRssResidual(stat, LD, mashPrior, V, covY)
    # The grid is in units of the residual standard deviations V starts from:
    # covY's on the original scale when V is not given, else 1 (a correlation)
    prior <- .mrmashRssDefaults(
        ss,
        mashPrior,
        canonicalPriorMatrices,
        S0,
        w0,
        gridSd = .mrmashGridSd(
            dots,
            residual$V %||% if (residual$originalScale) residual$covY
        )
    )
    callArgs <- .mrmashRssCallArgs(
        ss,
        .rssFitZ(stat, LD),
        LD,
        residual,
        prior,
        dots
    )
    .engineRequire(
        "mr.mashr::mr.mash.rss",
        names(callArgs),
        "mrmashRssWeights()"
    )
    .rssQuietMissingZ(exec(mr.mashr::mr.mash.rss, !!!callArgs))
}

# mr.mash.rss's arguments: the statistics (`ss`, .mrmashRssStats(), or the
# z-scores `z` the fit reads, .rssFitZ()) on the fit's scale, the
# covariances (`residual`, .mrmashRssResidual()), the prior (`prior`,
# .mrmashRssDefaults()) and the LD (`LD`), then the caller's options
# (`dots`).
# @noRd
.mrmashRssCallArgs <- function(ss, z, LD, residual, prior, dots) {
    # On one reference, n is one sample's: the median of the conditions'.
    # The LD is its matrix (R) or, unformed, its factor (X).
    c(
        .mrmashRssEffects(ss, z, residual$originalScale),
        list(
            n = LD$n %||% as.numeric(stats::median(ss$nVec)),
            S0 = prior$S0,
            w0 = prior$w0,
            # One thread: see .mrmashUnavailable().
            nthreads = 1L
        ),
        # The conditions' sample groups: mr.mash.rss would otherwise take
        # those on one reference with one n as one sample. The prior is
        # relative to the residual standard deviations unless `dots` says,
        # where mr.mashr can (.mrmashRelativeS0()).
        compact(list(
            covY = residual$covY,
            V = residual$V,
            R = LD$R,
            X = LD$X,
            sample_groups = LD$group,
            relative_S0 = if (
                is.null(as.list(dots)$relative_S0) && .mrmashRelativeS0()
            ) {
                TRUE
            }
        )),
        # `dots` is the caller's MrmashOptions() record; c() would append it as
        # one opaque element instead of splicing its entries.
        as.list(dots)
    )
}

# The prior (S0, w0) mr.mash.rss is given: each as supplied, else its default
# for the statistics `ss` (.mrmashRssStats()).
# @noRd
.mrmashRssDefaults <- function(
    ss,
    dataDrivenPriorMatrices,
    canonicalPriorMatrices,
    S0,
    w0,
    gridSd = 1
) {
    priorS0 <- S0 %||%
        buildMrmashPriorMatrices(
            Bhat = t(t(ss$Bhat) / gridSd),
            Shat = t(t(ss$Shat) / gridSd),
            K = ss$K,
            dataDrivenPriorMatrices = dataDrivenPriorMatrices,
            canonicalPriorMatrices = canonicalPriorMatrices
        )$S0
    list(S0 = priorS0, w0 = w0 %||% computeW0(ss$Bhat, length(priorS0)))
}

# The residual and Y covariances mr.mash.rss is given, and the scale of its
# fit: the original scale of the statistics' Bhat and Shat when covY is
# given, else the standardized scale of their z-scores. covY is as supplied,
# else the residual correlation the run's mash prior was fitted with, zero
# between conditions on different LD references (`LD$group`), which are
# fitted as different samples; without either, mr.mash.rss estimates it
# from the near-null z-scores of the conditions on each reference. V is as
# supplied, else each condition's residual variance under its heritability
# (`stat$h2`, .rssResidualVariances()), which mr.mash.rss gives covY's
# correlations; without either, mr.mash.rss starts from covY. A mr.mashr
# that predates fits across studies (.mrmashFitsStudies()) neither
# estimates covY nor takes V's variances alone: each is then the identity,
# as before, and a heritability cannot be used.
# @noRd
#' @importFrom checkmate assertMatrix
.mrmashRssResidual <- function(stat, LD, dataDrivenPriorMatrices, V, covY) {
    K <- ncol(as.matrix(stat$z))
    mash <- if (is.null(covY)) dataDrivenPriorMatrices$residualCorrelation
    if (!is.null(mash)) {
        assertMatrix(
            mash,
            nrows = K,
            ncols = K,
            .var.name = "dataDrivenPriorMatrices$residualCorrelation"
        )
    }
    originalScale <- !is.null(covY) &&
        !is.null(stat$Bhat) &&
        !is.null(stat$Shat)
    variance <- if (originalScale) diag(covY) else 1
    fitsStudies <- .mrmashFitsStudies()
    if (!fitsStudies && is.null(V) && !is.null(stat$h2)) {
        abort(glue(
            "mrmashRssWeights() with a heritability (`h2`) needs a newer ",
            "mr.mashr than the installed mr.mashr ",
            "{packageVersion('mr.mashr')}."
        ))
    }
    identity <- if (!fitsStudies) diag(K)
    list(
        covY = covY %||%
            .rssMaskBetweenRefs(unname(mash), LD$group) %||%
            identity,
        V = V %||% .rssResidualVariances(stat$h2, variance) %||% identity,
        originalScale = originalScale
    )
}

# The statistics mr.mash.rss fits: their Bhat and Shat (`ss`,
# .mrmashRssStats()) on the original scale, else their z-scores (`z`,
# .rssFitZ()), whose effects are per standard deviation of the variant and
# of the condition.
# @noRd
.mrmashRssEffects <- function(ss, z, originalScale) {
    if (originalScale) {
        return(list(Bhat = ss$Bhat, Shat = ss$Shat))
    }
    list(Z = z)
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
#'   \code{HibayesRssOptions} on summary statistics. Their wrappers set the
#'   iterations, burn-in, thinning and prior inclusion probability
#'   themselves (\code{nIter}, \code{nBurnin}, \code{thin}, \code{pi}), so
#'   every constructor refuses its engine's spelling of them, and hibayes'
#'   \code{seed} (drawn from R's random number generator) and
#'   \code{threads} (one). A wrapper's \code{h2} reaches its engine in that
#'   engine's spelling: \code{h2} for qgg, \code{R2} for BGLR, and for
#'   hibayes, which fixes its own at 0.5, the variances \code{vg} and
#'   \code{ve}. The constructors
#'   still take those, for a run that does not set \code{h2}; setting both
#'   is an error.
#'
#'   \code{BglrOptions} serves both of BGLR's functions: \code{BGLR()}, which
#'   fits one phenotype, and \code{Multitrait()}, which \code{bayesC} and
#'   \code{bayesN} run on several conditions at once (see
#'   \code{\link{bayesAlphabetWeights}}). It takes any option either one
#'   takes, and each call receives those it has, so a fit of one phenotype
#'   ignores \code{Multitrait()}'s \code{resCov} and \code{intercept}, and a
#'   joint fit ignores \code{BGLR()}'s \code{S0}, \code{df0} and the rest.
#'   It also takes \code{Cov}, the covariance of the effects across
#'   conditions, which goes into \code{Multitrait()}'s linear term as BGLR
#'   spells it there: \code{list(type = "UN")} (the default, with an
#'   inverse-Wishart prior), \code{"DIAG"}, \code{"FA"} or \code{"REC"}, with
#'   the prior settings that type takes (\code{df0}, \code{S0}, and the
#'   logical pattern \code{M} for \code{"FA"} and \code{"REC"}). BGLR draws
#'   the covariance afresh every iteration under that prior; it has no way
#'   to hold one fixed, so a large \code{df0} with \code{S0} set to the
#'   matrix it should stay near is the closest to supplying it. A matrix
#'   given there must have a row and a column per condition fitted.
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
#'   \code{GrpregOptions} and \code{AbessOptions} configure the
#'   \code{"grpreg"} and \code{"abess"} backends: grpreg's lasso, elastic
#'   net, MCP and SCAD for \code{lassoWeights}, \code{enetWeights},
#'   \code{mcpWeights} and \code{scadWeights}, and abess's best subset for
#'   \code{l0learnWeights}. Each is checked against its engine's named
#'   arguments (\code{cv.grpreg}'s and \code{grpreg}'s together, and
#'   \code{abess}'s), and each refuses the folds, the grouping and the
#'   family, which the wrappers set. \code{AbessOptions} also refuses
#'   \code{num.threads} (abess runs on one thread, for the reason
#'   \code{parallel} is refused above), \code{weight} (abess scales a sparse
#'   design's rows by the weights but its phenotype by their square roots)
#'   and \code{init.active.set} (abess replaces it with an index of its
#'   own).
#'
#'   \code{LassosumOptions} and \code{PenalizedRssOptions} hold the path:
#'   \code{lambda}, \code{thr}, \code{maxiter} and \code{strategy}, how the
#'   path is walked (\code{"strongRule"} or \code{"gradientRanked"}; see
#'   \code{\link{penalizedRss}}).
#'
#'   \code{MrashOptions} refuses \code{update.order}, which
#'   \code{mrashWeights()} sets from its own \code{updateOrder}.
#'
#'   \code{SparsevbOptions} takes \code{sparsevb::svb.fit}'s settings. The
#'   per-variant ones (\code{mu}, \code{sigma}, \code{gamma} and
#'   \code{update_order}) index the columns \code{sparsevbWeights()} fits,
#'   which leave out zero-variance ones, and \code{svb.fit} reads
#'   \code{update_order} counting from 0, though its documentation says 1.
#'
#'   \code{MrmashOptions} and \code{MrmashRssOptions} refuse
#'   \code{nthreads}: mr.mash itself runs on one thread, because on macOS
#'   its worker threads can later hang R inside the OpenMP runtime.
#'
#'   \code{MrmashOptions} also takes \code{vEigenFloor}, which is
#'   \code{\link{mrmashWrapper}}'s rather than mr.mash's: the smallest
#'   eigenvalue the estimated residual covariance may have, as a fraction of
#'   its spectral norm. Default \code{1e-3}. Below it, as near-duplicate
#'   contexts make it, the eigenvalues are raised to the floor, the covariance
#'   is not updated during the fit, and a warning says so.
#'
#'   It also takes \code{vInitMethod}, likewise the wrapper's: where the
#'   residual covariance starts. \code{"data"} (the default) starts it at the
#'   covariance of Y; \code{"mash"} at the residual correlation of the mash
#'   prior (\code{\link{mashPipeline}}), which \code{twasWeightsPipeline}
#'   passes in \code{dataDrivenPriorMatrices}, scaled by the standard
#'   deviations of Y's columns. Either start is zero for a pair of conditions
#'   never observed together, and is floored by \code{vEigenFloor}.
#'
#'   \code{SdprOptions} also takes \code{sharing} and \code{rho}, which act
#'   only when \code{dprGibbs} is fitted across studies as SDPRX fits
#'   ancestries (\code{\link{twasWeightsPipeline}}'s \code{coupling}):
#'   which sets of the studies a component's effects lie in, and their
#'   correlation (see \code{\link{sdpr}}). A single study ignores them.
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
        .BGLR_CALLEES,
        defaults = list(),
        extra = extra,
        label = "BglrOptions",
        engine = "bglr",
        accepted = .bglrAcceptedNames(),
        filtered = TRUE
    )
}

# BGLR's two functions: BGLR() on one phenotype, Multitrait() on several
# conditions at once. bglrWeights() hands each only the options it takes
# (.bglrOwnArgs()), so one record serves both.
# @noRd
.BGLR_CALLEES <- c("BGLR::BGLR", "BGLR::Multitrait")

# The names BglrOptions() takes: either function's arguments, and the
# settings of a Multitrait() linear term (.bglrTermSettings()), which
# bglrWeights() moves into that term. NULL, checking nothing, where BGLR's
# arguments cannot be read.
# @noRd
.bglrAcceptedNames <- function() {
    formalNames <- .engineAcceptedNames(.BGLR_CALLEES, filtered = TRUE)
    if (is.null(formalNames)) {
        return(NULL)
    }
    c(formalNames, fn_fmls_names(.bglrTermSettings))
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
GrpregOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            X = "supplied from the data by the pipeline",
            y = "supplied from the data by the pipeline",
            group = str_c(
                "the wrapper: one group per variant, or a variant's effects ",
                "across conditions"
            ),
            penalty = str_c(
                "the calling wrapper (lassoWeights, enetWeights, mcpWeights, ",
                "scadWeights)"
            ),
            family = "the wrapper, which fits quantitative phenotypes",
            nfolds = str_c(
                "the wrapper's `nfolds` (mcpWeights, scadWeights; ",
                "lassoWeights and enetWeights use 5)"
            ),
            fold = "the wrapper, which deals each sample into a fold",
            seed = "R's random number generator, which deals the folds",
            returnX = "the wrapper's `retainX`"
        ),
        "GrpregOptions"
    )
    .newMethodOptions(
        c("grpreg::cv.grpreg", "grpreg::grpreg"),
        defaults = list(),
        extra = extra,
        label = "GrpregOptions",
        engine = "grpreg",
        filtered = TRUE
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
AbessOptions <- function(...) {
    extra <- list(...)
    .configRefuseUnavailable(
        extra,
        c(
            num.threads = str_c(
                "abess runs on one thread, as weight fitting already runs ",
                "inside BiocParallel workers; abess's own default takes every ",
                "core"
            ),
            weight = str_c(
                "abess scales a sparse design's rows by the weights but its ",
                "phenotype by their square roots, and a fit across conditions ",
                "on their own samples is sparse"
            ),
            init.active.set = str_c(
                "abess replaces it with an index of its own before fitting"
            )
        ),
        "AbessOptions"
    )
    .configRefuseOwned(
        extra,
        c(
            x = "supplied from the data by the pipeline",
            y = "supplied from the data by the pipeline",
            family = "the wrapper, from the number of conditions",
            group.index = "the wrapper: a variant's effects across conditions",
            tune.type = str_c(
                "the wrapper, which picks the support size by ",
                "cross-validation"
            ),
            nfolds = "the wrapper's `nFolds` (l0learnWeights)",
            foldid = "the wrapper, which deals each sample into a fold",
            seed = "R's random number generator, which deals the folds",
            fit.intercept = "the wrapper, which centres each condition"
        ),
        "AbessOptions"
    )
    .newMethodOptions(
        "abess::abess.default",
        defaults = list(),
        extra = extra,
        label = "AbessOptions",
        engine = "abess",
        filtered = TRUE
    )
}

#' @rdname twasEngineOptions
#' @export
SparsevbOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            X = "supplied from the data by the pipeline",
            Y = "supplied from the data by the pipeline"
        ),
        "SparsevbOptions"
    )
    .newMethodOptions(
        "sparsevb::svb.fit",
        defaults = list(),
        extra = list(...),
        label = "SparsevbOptions",
        engine = "sparsevb"
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
            y = "supplied from the data by the pipeline",
            update.order = "mrashWeights()' own `updateOrder`"
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
    .mrmashCheckVEigenFloor(extra$vEigenFloor)
    .mrmashCheckVInitMethod(extra$vInitMethod)
    .newMethodOptions(
        "mr.mashr::mr.mash",
        defaults = list(),
        extra = list(...),
        label = "MrmashOptions",
        engine = "mrmash",
        # The bag is split before the engine sees it: .mrmashSplitPriorArgs
        # routes the prior-covariance names to MrmashPriorParam() and the rest
        # to mr.mash. So the accepted set is both, read from the live
        # formals of each rather than transcribed, plus the wrapper's own
        # vEigenFloor and vInitMethod, which never reach mr.mash.
        accepted = unique(c(
            .engineAcceptedNames("mr.mashr::mr.mash"),
            names(formals(MrmashPriorParam)),
            "vEigenFloor",
            "vInitMethod"
        ))
    )
}

# @noRd
.mrmashCheckVInitMethod <- function(vInitMethod) {
    if (is.null(vInitMethod)) {
        return(invisible(NULL))
    }
    isMethod <- is.character(vInitMethod) &&
        length(vInitMethod) == 1L &&
        is_in(vInitMethod, c("data", "mash"))
    if (!isMethod) {
        abort("MrmashOptions: `vInitMethod` must be \"data\" or \"mash\".")
    }
    invisible(NULL)
}

# @noRd
.mrmashCheckVEigenFloor <- function(vEigenFloor) {
    if (is.null(vEigenFloor)) {
        return(invisible(NULL))
    }
    isFraction <- is.numeric(vEigenFloor) &&
        length(vEigenFloor) == 1L &&
        !is.na(vEigenFloor) &&
        vEigenFloor >= 0 &&
        vEigenFloor < 1
    if (!isFraction) {
        abort("MrmashOptions: `vEigenFloor` must be one number in [0, 1).")
    }
    invisible(NULL)
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
            Z = "supplied from the data by the pipeline",
            R = "supplied from the data by the pipeline",
            X = "supplied from the data by the pipeline",
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
    user <- as.list(methodArgs)
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
        !!!user
    )
    exec(glmnet::cv.glmnet, !!!callArgs)
}
#' Compute TWAS weights via penalized regression (glmnet)
#'
#' Fit a cross-validated elastic-net / lasso model with \code{glmnet::cv.glmnet}
#' and return the per-variant coefficient vector at \code{lambda.min}.
#'
#' Given a matrix \code{y}, a column per condition measured on the same
#' samples, it fits the conditions at once with glmnet's \code{"mgaussian"}
#' family: each variant's effects across them form one group under the
#' penalty, so a variant is in the model for every condition or for none,
#' each condition with its own effect. Each phenotype is scaled to unit
#' variance for the fit and its weights are returned on its own scale. One
#' \code{lambda} serves every condition, the one with the lowest
#' cross-validated error summed over them. At \code{alpha = 0} the
#' conditions are separate ridge fits sharing that \code{lambda}. A missing
#' phenotype is an error: conditions measured on samples of their own take
#' \code{backend = "grpreg"} (see \code{\link{lassoWeights}}).
#'
#' @param X Numeric genotype / design matrix (samples x variants).
#' @param y Numeric response (phenotype) vector of length \code{nrow(X)}, or
#'   a matrix with a column per condition (see Details).
#' @param alpha Elastic-net mixing parameter in \code{[0, 1]}: \code{1} = lasso,
#'   \code{0} = ridge, \code{0.5} = elastic net.
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights: \code{"none"} (default) attaches
#'   nothing; \code{"slim"} the selected fit, as a list of its coefficients
#'   (\code{coef}, over the variants fitted, a column per condition for a
#'   matrix \code{y}), \code{intercept}, \code{lambda} and cross-validated
#'   error there (\code{cvm}, \code{cvsd}); \code{"full"} the whole
#'   \code{cv.glmnet} object, with every \code{lambda}'s coefficients but
#'   without the copy of \code{X} its stored calls would hold, fitted for a
#'   matrix \code{y} on the scaled phenotypes. Each names the variants it
#'   was fitted on.
#' @param retainX Logical. Keep that copy of \code{X} in a \code{"full"} fit
#'   (default \code{FALSE}). An error at any other \code{fitRetention}.
#' @param methodArgs Additional arguments forwarded to
#'   \code{glmnet::cv.glmnet}, built with \code{\link{GlmnetOptions}}. A bare
#'   list is refused, since it cannot be checked against the engine's formals.
#' @return Numeric vector of length \code{ncol(X)} of per-variant weights, or
#'   for a matrix \code{y} a matrix of them with a column per condition, with
#'   the fit attached as attribute \code{"fit"} unless
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
    retainX = FALSE,
    methodArgs = GlmnetOptions()
) {
    fitRetention <- arg_match(fitRetention)
    .assertRetainX(retainX, fitRetention)
    if (.isConditionMatrix(y)) {
        return(.glmnetGroupedWeights(
            X,
            y,
            alpha,
            fitRetention,
            retainX,
            methodArgs
        ))
    }
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
    fit <- switch(
        fitRetention,
        slim = .glmnetSlimFit(enet, coefs[-1L], coefs[[1L]]),
        full = if (retainX) enet else .glmnetWithoutData(enet)
    )
    `attr<-`(weights, "fit", .withFitVariantIds(fit, colnames(X)[keep]))
}

# cv.glmnet's whole object without the data it was given: called with values
# rather than names, its stored calls hold them inline.
# @noRd
.glmnetWithoutData <- function(cvFit) {
    inner <- cvFit$glmnet.fit
    `[[<-`(
        `[[<-`(cvFit, "call", value = .callWithoutData(cvFit$call)),
        "glmnet.fit",
        value = `[[<-`(inner, "call", value = .callWithoutData(inner$call))
    )
}

# A stored call naming its data rather than holding it.
# @noRd
.callWithoutData <- function(call) {
    held <- intersect(c("x", "y"), names(call))
    as.call(replace(as.list(call), held, map(held, as.name)))
}

# `retainX` keeps a fit's copy of X, which only a "full" fit can hold.
# @noRd
.assertRetainX <- function(retainX, fitRetention) {
    assertFlag(retainX)
    if (retainX && !identical(fitRetention, "full")) {
        abort(glue(
            "`retainX = TRUE` requires `fitRetention = \"full\"`, not ",
            "\"{fitRetention}\"."
        ))
    }
    invisible(NULL)
}

# The fit cv.glmnet selected: its coefficients over the variants fitted (a
# column per condition for several), the intercept, its lambda and the
# cross-validated error there.
# @noRd
.glmnetSlimFit <- function(cvFit, coef, intercept) {
    at <- match(cvFit$lambda.min, cvFit$lambda)
    list(
        coef = coef,
        intercept = intercept,
        lambda = cvFit$lambda.min,
        cvm = cvFit$cvm[[at]],
        cvsd = cvFit$cvsd[[at]]
    )
}

# glmnet's fit of several conditions (glmnetWeights()): its "mgaussian"
# family, which takes no missing phenotype, on the samples observing every
# condition (.completeRows()), each condition's phenotype scaled to unit
# variance and its weights put back on its own scale.
# @noRd
.glmnetGroupedWeights <- function(X, Y, alpha, retention, retainX, methodArgs) {
    rows <- .completeRows(Y, "glmnetWeights", "backend = \"grpreg\"")
    observed <- Y[rows, , drop = FALSE]
    keep <- .dropZeroVariance(X[rows, , drop = FALSE], "glmnetWeights")
    scales <- .groupedScales(observed)
    cvFit <- fitGlmnetCv(
        x = X[rows, keep, drop = FALSE],
        y = .groupedScaled(observed, scales),
        alpha = alpha,
        family = "mgaussian",
        standardize = FALSE,
        intercept = TRUE,
        nfolds = 5,
        methodArgs = methodArgs
    )
    coefs <- .glmnetGroupedCoef(cvFit, scales)
    fitted <- list(
        coef = coefs$coef,
        slim = .glmnetSlimFit(cvFit, coefs$coef, coefs$intercept),
        object = cvFit
    )
    .engineWeights(fitted, keep, X, Y, retention, .glmnetFullFit(retainX))
}

# An mgaussian fit's weights at lambda.min, a column per condition, and each
# condition's intercept, back on the scale of its phenotype.
# @noRd
.glmnetGroupedCoef <- function(cvFit, scales) {
    B <- do.call(cbind, map(coef(cvFit, s = "lambda.min"), as.numeric))
    list(
        coef = sweep(B[-1L, , drop = FALSE], 2L, scales$spread, "*"),
        intercept = scales$center + scales$spread * B[1L, ]
    )
}

# How a "full" glmnet fit is kept: whole, or without its copies of the data.
# @noRd
.glmnetFullFit <- function(retainX) {
    if (retainX) identity else .glmnetWithoutData
}

#' Compute TWAS weights via elastic net
#'
#' The elastic net at a mixing between the lasso and ridge, \code{alpha =
#' 0.5} unless set, on glmnet (\code{\link{glmnetWeights}}, the default) or
#' on grpreg. A matrix \code{y}, a column per condition, is fitted across
#' the conditions at once, as \code{\link{lassoWeights}} describes.
#' @param alpha Elastic-net mixing in \code{[0, 1]}: the L1 share of the
#'   penalty. Default 0.5. grpreg needs it above 0, and across conditions it
#'   scales the two shares of a variant's penalty differently from glmnet,
#'   so the same \code{alpha} is not the same penalty on the two backends.
#' @inheritParams lassoWeights
#' @return Numeric vector of length \code{ncol(X)} of per-variant weights, or
#'   for a matrix \code{y} a matrix of them with a column per condition, with
#'   the fit attached as attribute \code{"fit"} unless
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
    backend = NULL,
    fitRetention = c("none", "slim", "full"),
    retainX = FALSE,
    methodArgs = NULL
) {
    assertNumber(alpha, lower = 0, upper = 1)
    .enetFamilyWeights(
        X,
        y,
        alpha,
        "enet",
        backend,
        arg_match(fitRetention),
        retainX,
        methodArgs
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
    retainX = FALSE,
    methodArgs = GlmnetOptions()
) {
    .assertMethodOptions(methodArgs, "GlmnetOptions", "methodArgs")
    glmnetWeights(
        X,
        y,
        0,
        fitRetention = arg_match(fitRetention),
        retainX = retainX,
        methodArgs = methodArgs
    )
}

#' Compute TWAS weights via the lasso
#'
#' The lasso, on glmnet (\code{\link{glmnetWeights}} at \code{alpha = 1}, the
#' default) or on grpreg (\code{backend = "grpreg"}), whose group lasso with
#' each variant its own group is the lasso. grpreg always standardizes the
#' columns of \code{X}, where glmnet here fits them as given.
#'
#' Given a matrix \code{y}, a column per condition, it fits the conditions
#' at once: each variant's effects across them form one group under the
#' penalty, so a variant is in the model for every condition or for none,
#' each condition with its own effect, fitted on its own samples. The loss is
#' each condition's mean squared error on its phenotype scaled to unit
#' variance, summed over the conditions, and one \code{lambda} serves them
#' all: the one with the lowest cross-validated error summed over them, on
#' folds that hold each sample out of every condition it is in. Weights come
#' back on each condition's own scale. glmnet fits conditions measured on
#' the same samples (its \code{"mgaussian"} family). grpreg also fits
#' conditions measured on samples of their own, as several studies' are,
#' their phenotypes missing on the others' samples: each condition's rows
#' hold its own samples, centred within them, against its own copy of the
#' variants, and a variant's copies form its group. That design is dense,
#' \eqn{8 N C p} bytes for \eqn{N} rows in all, \eqn{C} conditions and
#' \eqn{p} variants.
#' @param backend The engine: \code{"glmnet"} (the default for \code{NULL})
#'   or \code{"grpreg"}.
#' @param methodArgs Options for the backend's engine, built with its
#'   constructor: \code{\link{GlmnetOptions}} for glmnet,
#'   \code{\link{GrpregOptions}} for grpreg. \code{NULL} (the default)
#'   keeps the engine's defaults. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be silently
#'   ignored.
#' @param retainX Logical. Keep the fit's copy of \code{X} in a
#'   \code{"full"} fit (default \code{FALSE}): glmnet's stored calls hold
#'   it, and on grpreg it is grpreg's standardized design (its
#'   \code{returnX}). An error at any other \code{fitRetention}.
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights: \code{"none"} (default) attaches
#'   nothing; \code{"slim"} the selected fit, as a list of its coefficients
#'   (\code{coef}, over the variants fitted, a column per condition for a
#'   matrix \code{y}), \code{intercept}, \code{lambda} and the
#'   cross-validated error there (\code{cvm} and \code{cvsd} on glmnet,
#'   \code{cve} and \code{cvse} on grpreg); \code{"full"} the whole
#'   \code{cv.glmnet} or \code{cv.grpreg} object, with every
#'   \code{lambda}'s coefficients but without its copies of the data. Each
#'   names the variants it was fitted on.
#' @inheritParams glmnetWeights
#' @return Numeric vector of length \code{ncol(X)} of per-variant weights, or
#'   for a matrix \code{y} a matrix of them with a column per condition, with
#'   the fit attached as attribute \code{"fit"} unless
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
    backend = NULL,
    fitRetention = c("none", "slim", "full"),
    retainX = FALSE,
    methodArgs = NULL
) {
    .enetFamilyWeights(
        X,
        y,
        1,
        "lasso",
        backend,
        arg_match(fitRetention),
        retainX,
        methodArgs
    )
}

# The lasso or elastic net (`method`) on its backend: glmnet
# (glmnetWeights()), or grpreg's group lasso (grpregWeights()), one group per
# variant for one phenotype.
# @noRd
.enetFamilyWeights <- function(
    X,
    y,
    alpha,
    method,
    backend,
    retention,
    retainX,
    methodArgs
) {
    engine <- .methodBackend(method, backend, "individual")
    options <- .backendOptions(
        methodArgs,
        .methodBackends$individual[[method]][[engine]],
        engine
    )
    if (identical(engine, "glmnet")) {
        return(glmnetWeights(X, y, alpha, retention, retainX, options))
    }
    .grpregAssertAlpha(options, alpha)
    grpregWeights(
        X,
        y,
        "grLasso",
        alpha = alpha,
        fitRetention = retention,
        retainX = retainX,
        methodArgs = options
    )
}

# The lasso and elastic net take their mixing from the wrapper on grpreg as
# on glmnet, and grpreg's penalties need some of the lasso in them.
# @noRd
.grpregAssertAlpha <- function(methodArgs, alpha) {
    if (is_in("alpha", names(methodArgs))) {
        abort(str_c(
            "GrpregOptions(alpha=) does not apply here: enetWeights takes its ",
            "own `alpha`, and lassoWeights fixes it at 1."
        ))
    }
    if (alpha <= 0) {
        abort("grpreg needs `alpha` > 0: it has no ridge fit.")
    }
    invisible(NULL)
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
#' @param updateOrder The order the coordinate ascent updates the variants
#'   in: \code{"consecutive"} (default), their column order, or
#'   \code{"decreasingLogBf"}, by decreasing marginal log Bayes factor under
#'   the prior \code{mr.ash} fits (its grid \code{sa2} and weights \code{pi},
#'   from \code{methodArgs} or its defaults, and residual variance
#'   \code{sigma2}, else the variance of \code{y}). It is passed as
#'   \code{susieR::mr.ash}'s \code{update.order}, which
#'   \code{\link{MrashOptions}} therefore refuses. Either way the fit starts
#'   from lasso weights unless \code{methodArgs} gives \code{beta.init}.
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights: \code{"none"} (default) attaches
#'   nothing, \code{"slim"} the fitted model without its copy of the data
#'   (\code{data}) or the order it updated the variants in each iteration
#'   (\code{update.order}), \code{"full"} the whole fitted model but for its
#'   copy of \code{X}. Each names the variants it was fitted on.
#' @param retainX Logical. Keep that copy of \code{X} in a \code{"full"} fit
#'   (default \code{FALSE}). An error at any other \code{fitRetention}.
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
    updateOrder = c("consecutive", "decreasingLogBf"),
    fitRetention = c("none", "slim", "full"),
    retainX = FALSE,
    methodArgs = MrashOptions()
) {
    .assertMethodOptions(methodArgs, "MrashOptions", "methodArgs")
    ordering <- arg_match(updateOrder)
    retention <- arg_match(fitRetention)
    .assertRetainX(retainX, retention)
    keep <- .dropZeroVariance(X, "mrashWeights")
    XKeep <- X[, keep, drop = FALSE]
    # as.list(): list_assign() rejects a MethodOptions outright as `.x`, the
    # same S4-SimpleList seam that makes c() append it instead of splicing.
    betaInit <- .mrashBetaInit(methodArgs, XKeep, y, X, keep)
    orderArg <- .mrashUpdateOrder(ordering, XKeep, y, methodArgs)
    argsList <- list_assign(as.list(methodArgs), !!!betaInit, !!!orderArg)
    fit.mr.ash <- exec(mr.ash, X = XKeep, y = y, !!!argsList)
    # Zero-variance columns were never fitted and keep a zero weight.
    eff.wgt <- replace(
        rep(0, ncol(X)),
        keep,
        predict(fit.mr.ash, type = "coefficients")[-1]
    )
    if (identical(retention, "none")) {
        return(eff.wgt)
    }
    kept <- switch(
        retention,
        slim = .withoutFields(fit.mr.ash, c("data", "update.order")),
        full = if (retainX) {
            fit.mr.ash
        } else {
            `[[<-`(
                fit.mr.ash,
                "data",
                value = .withoutFields(fit.mr.ash$data, "X")
            )
        }
    )
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
        bayesC = c(
            qgg = "QggOptions",
            hibayes = "HibayesOptions",
            bglr = "BglrOptions"
        ),
        bayesL = c(
            qgg = "QggOptions",
            bglr = "BglrOptions",
            hibayes = "HibayesOptions"
        ),
        bayesN = c(
            qgg = "QggOptions",
            hibayes = "HibayesOptions",
            bglr = "BglrOptions"
        ),
        bayesR = c(qgg = "QggOptions", hibayes = "HibayesOptions")
    ),
    # On summary statistics every engine reads the LD matrix, formed from the
    # panel for a genotype reference. bayesB has hibayes alone there.
    sumstat = list(
        bayesA = c(qgg = "QggOptions", hibayes = "HibayesRssOptions"),
        bayesB = c(hibayes = "HibayesRssOptions"),
        bayesC = c(qgg = "QggOptions", hibayes = "HibayesRssOptions"),
        bayesL = c(hibayes = "HibayesRssOptions", qgg = "QggOptions"),
        bayesN = c(qgg = "QggOptions", hibayes = "HibayesRssOptions"),
        bayesR = c(qgg = "QggOptions", hibayes = "HibayesRssOptions")
    )
)

# The engines of every method with a `backend` on a data path: the Bayesian
# alphabet (.alphabetBackends); on individual-level data the lasso, elastic
# net, MCP and SCAD, which also run on grpreg, and L0Learn's best subset,
# which also runs on abess; and on summary statistics the lasso. Each engine
# maps to its Options constructor; the first is the default. A method with
# one engine there (bayesB on summary statistics) leaves nothing to choose
# (.methodHasBackends()).
# @noRd
.methodBackends <- list(
    individual = c(
        .alphabetBackends$individual,
        list(
            lasso = c(glmnet = "GlmnetOptions", grpreg = "GrpregOptions"),
            enet = c(glmnet = "GlmnetOptions", grpreg = "GrpregOptions"),
            mcp = c(ncvreg = "NcvregOptions", grpreg = "GrpregOptions"),
            scad = c(ncvreg = "NcvregOptions", grpreg = "GrpregOptions"),
            l0learn = c(l0learn = "L0learnOptions", abess = "AbessOptions")
        )
    ),
    sumstat = c(
        .alphabetBackends$sumstat,
        list(
            lasso = c(
                pecotmr = "LassosumOptions",
                bigsnpr = "BigsnprOptions"
            )
        )
    )
)

# Whether `method` has engines to choose between on a data path.
# @noRd
.methodHasBackends <- function(method, path) {
    length(.methodBackends[[path]][[method]]) > 1L
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
    bglr = c(
        bayesB = "BayesB",
        bayesC = "BayesC",
        bayesL = "BL",
        bayesN = "BRR"
    ),
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
        hibayes = ".alphabetHibayesRssFit"
    )
)

# BGLR's models with a form in its Multitrait(), which fits several
# conditions' phenotypes at once: a spike and slab whose slab is correlated
# across conditions, and a Gaussian prior with a covariance across them.
# BayesB and BL have none.
# @noRd
.bglrMultitraitModels <- c(BayesC = "SpikeSlab", BRR = "BRR")

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
#' \code{thin}, \code{pi} and \code{h2} are renamed to each engine's own
#' spelling, and everything else goes through \code{methodArgs}, built with
#' the backend's own constructor.
#'
#' On \code{"bglr"}, \code{bayesC} and \code{bayesN} also fit several
#' conditions at once: given a matrix \code{y}, one column per condition,
#' they run BGLR's \code{Multitrait()}. \code{bayesC} becomes its spike and
#' slab, with each condition's own inclusion indicators and a slab
#' correlated across conditions, and \code{bayesN} its Gaussian prior with a
#' covariance across conditions. A missing phenotype is drawn within the
#' sampler, so each condition may be measured on different samples, and
#' \code{pi} starts every condition's inclusion probability. The covariance
#' of the effects across conditions and the residual covariance are learned,
#' under the structures and priors \code{BglrOptions(Cov = )} and
#' \code{BglrOptions(resCov = )} give them.
#'
#' @param X A numeric matrix of genotypes (samples x variants).
#' @param y A numeric vector of phenotypes, one per row of \code{X}; or, for
#'   \code{bayesC} and \code{bayesN} on \code{"bglr"}, a matrix with one
#'   column per condition, which may hold missing values (see Description).
#' @param method The prior: \code{"bayesA"}, \code{"bayesB"},
#'   \code{"bayesC"}, \code{"bayesL"}, \code{"bayesN"} or \code{"bayesR"}.
#' @param backend The engine, one the method has (see Description);
#'   \code{NULL} picks the method's default.
#' @param pi For \code{bayesB} and \code{bayesC}, the prior proportion of
#'   non-null effects, where each engine's estimate of it starts; \code{NULL}
#'   gives 0.2 for \code{bayesB} and 0.1 for \code{bayesC}. Refused for the
#'   other methods.
#' @param h2Estimate,h2EstimateParam The trait's heritability in the window,
#'   which sets every engine's starting variances and the scales of its
#'   priors: \code{h2Estimate}, a number in [0, 1) or a one-row
#'   \code{\linkS4class{H2Estimate}}; or \code{h2EstimateParam}, the
#'   \code{\link{H2EstimateParam}} settings to estimate it under first, on
#'   \code{X} and \code{y}, in-sample (as \code{\link{estimateLocalH2}}
#'   estimates a \code{QtlDataset}'s traits). Neither (the default) keeps
#'   each engine's own, 0.5; both is an error. It is used floored at 0.001
#'   and capped at 0.99, and an estimate that cannot be made leaves the
#'   engine's default, with a warning. For a matrix \code{y} the estimate is
#'   made on each column's observed samples and their mean is used, as
#'   \code{Multitrait()} takes one. qgg reads it as \code{h2} and BGLR as
#'   \code{R2}; hibayes, which fixes its own at 0.5, gets the genetic and
#'   residual variances (\code{vg}, \code{ve}) its defaults would give at
#'   it. Setting that engine's spelling through \code{methodArgs} as well is
#'   an error.
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
#'   \code{X}), or for a matrix \code{y} a matrix of them with one column
#'   per condition; columns dropped for zero variance receive weight 0.
#'   Unless \code{fitRetention = "none"} the fit is attached as attribute
#'   \code{"fit"}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' bayesAlphabetWeights(X, y, "bayesR", nIter = 1000, nBurnin = 200)
#' @export
bayesAlphabetWeights <- function(
    X,
    y,
    method,
    backend = NULL,
    pi = NULL,
    h2Estimate = NULL,
    h2EstimateParam = NULL,
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
    methodArgs <- .alphabetOptions(methodArgs, "individual", method, backend)
    response <- .alphabetResponse(y, X, method, backend)
    keep <- .dropZeroVariance(X, "bayesAlphabetWeights")
    XKeep <- X[, keep, drop = FALSE]
    fit <- exec(
        .alphabetFitters$individual[[backend]],
        X = XKeep,
        y = response,
        method = method,
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        pi = pi,
        h2 = .alphabetIndividualH2(
            h2Estimate,
            h2EstimateParam,
            XKeep,
            response
        ),
        methodArgs = methodArgs
    )
    .alphabetResult(fit, backend, fitRetention, keep, X, response)
}

# The engine options of a method on `path` ("individual" or "sumstat") and
# its `backend`, checked against the engine (.backendOptions()).
# @noRd
.alphabetOptions <- function(methodArgs, path, method, backend) {
    .backendOptions(
        methodArgs,
        .alphabetBackends[[path]][[method]][[backend]],
        backend
    )
}

# A fit's weights over every column of X (.weightsOnColumns()), with its
# posterior summaries, or the engine's whole fit, attached unless
# `fitRetention` is "none".
# @noRd
.alphabetResult <- function(fit, backend, fitRetention, keep, X, response) {
    posterior <- .alphabetPosterior(fit, backend)
    weights <- .weightsOnColumns(
        posterior$coef,
        keep,
        colnames(X),
        colnames(response)
    )
    if (identical(fitRetention, "none")) {
        return(weights)
    }
    kept <- if (identical(fitRetention, "full")) fit else posterior
    `attr<-`(weights, "fit", .withFitVariantIds(kept, colnames(X)[keep]))
}

# The phenotypes a fit reads: a vector, one per row of X, or, for a method
# whose engine fits several conditions at once (.alphabetMultitraitMethods()
# on BGLR), a matrix with a column per condition. A one-column matrix is that
# column.
# @noRd
#' @importFrom checkmate assertMatrix assertVector
.alphabetResponse <- function(y, X, method, backend) {
    if (!is.matrix(y)) {
        assertVector(y, len = nrow(X))
        return(y)
    }
    assertMatrix(y, mode = "numeric", nrows = nrow(X), min.cols = 1L)
    if (ncol(y) == 1L) {
        return(y[, 1L])
    }
    multitrait <- .alphabetMultitraitMethods()
    if (!identical(backend, "bglr") || !is_in(method, multitrait)) {
        abort(glue(
            "bayesAlphabetWeights: only backend = \"bglr\" fits a matrix ",
            "`y` (several conditions at once), for ",
            "{str_flatten(multitrait, ' and ')}; {method} on {backend} ",
            "takes a vector."
        ))
    }
    y
}

# The alphabet methods BGLR also fits on several conditions at once: those
# whose BGLR model has a Multitrait() form (.bglrMultitraitModels).
# @noRd
.alphabetMultitraitMethods <- function() {
    models <- .alphabetEngineMethods$bglr
    names(models)[is_in(models, names(.bglrMultitraitModels))]
}

# Weights over every column of X from a fit over the kept ones (`keep`), zero
# at a dropped column. A fit of several conditions gives a matrix, with a row
# per column of X (`ids`) and a column per condition.
# @noRd
.weightsOnColumns <- function(coef, keep, ids, conditions = NULL) {
    if (!is.matrix(coef)) {
        return(replace(rep(0, length(keep)), keep, coef))
    }
    rowOf <- replace(cumsum(keep), !keep, nrow(coef) + 1L)
    `dimnames<-`(
        rbind(coef, 0)[rowOf, , drop = FALSE],
        list(ids, conditions)
    )
}

# The heritability an individual-level fit starts from (.h2Resolve()), an
# estimate being made on the variants it fits. For several conditions the
# estimate is made on each one's observed samples and the estimates' mean is
# used, NULL when none could be made.
# @noRd
.alphabetIndividualH2 <- function(h2Estimate, h2EstimateParam, X, y) {
    .h2AssertOneSource(h2Estimate, h2EstimateParam, "bayesAlphabetWeights")
    if (!is.matrix(y) || is.null(h2EstimateParam)) {
        return(.alphabetConditionH2(h2Estimate, h2EstimateParam, X, y))
    }
    estimates <- compact(map(
        seq_len(ncol(y)),
        .alphabetColumnH2,
        h2EstimateParam = h2EstimateParam,
        X = X,
        y = y
    ))
    if (length(estimates) == 0L) NULL else mean(list_c(estimates))
}

# Column `k`'s heritability, on the samples it was measured on.
# @noRd
.alphabetColumnH2 <- function(k, h2EstimateParam, X, y) {
    seen <- !is.na(y[, k])
    .alphabetConditionH2(
        NULL,
        h2EstimateParam,
        X[seen, , drop = FALSE],
        y[seen, k]
    )
}

# One phenotype's heritability (.h2Resolve()), estimated in-sample.
# @noRd
.alphabetConditionH2 <- function(h2Estimate, h2EstimateParam, X, y) {
    .h2Resolve(
        h2Estimate,
        h2EstimateParam,
        "bayesAlphabetWeights",
        .h2IndividualFitParam,
        list(X = X, y = y),
        individual = TRUE
    )
}

# Each backend's spelling of the heritability its priors start from.
# hibayes fixes its own at 0.5, so it takes the variances that implies.
# @noRd
.alphabetH2Names <- list(
    qgg = "h2",
    bglr = "R2",
    hibayes = c("vg", "ve"),
    pecotmr = "h2"
)

# `methodArgs` with the heritability `h2` in its backend's spelling, added
# through the constructor `ctorName`, so a value the caller also set there
# is refused as given twice. hibayes gets the genetic and residual variances
# its own defaults give at 0.5 (src/Bayes.cpp, src/SBayesD.cpp), at `h2`:
# (dfvg - 2) / dfvg vary h2 and vary (1 - h2), `vary` being the phenotypic
# variance it computes from its input.
# @noRd
.alphabetWithH2 <- function(methodArgs, ctorName, backend, h2, vary = NULL) {
    if (is.null(h2)) {
        return(methodArgs)
    }
    values <- if (identical(backend, "hibayes")) {
        dfvg <- as.list(methodArgs)[["dfvg"]] %||% 4
        c((dfvg - 2) / dfvg * vary * h2, vary * (1 - h2))
    } else {
        h2
    }
    exec(
        .methodOptionsWith,
        methodArgs,
        ctorName,
        !!!set_names(as.list(values), .alphabetH2Names[[backend]])
    )
}

# A backend's posterior summaries over the variants it fitted: the effects'
# posterior means and, where it reports them, their posterior inclusion
# probabilities. A fit of several conditions keeps them as matrices, a
# column per condition.
# @noRd
.alphabetPosterior <- function(fit, backend) {
    parts <- switch(
        backend,
        qgg = list(coef = fit$bm, pip = fit$dm),
        bglr = .bglrPosterior(fit),
        hibayes = list(coef = fit$alpha, pip = fit$pip)
    )
    map(compact(parts), .posteriorValues)
}

# BGLR's posterior means for its one linear term. Multitrait() keeps a
# condition's effects as `beta`, the slab draw times its inclusion, where
# BGLR() keeps the effects as `b`.
# @noRd
.bglrPosterior <- function(fit) {
    term <- pluck(fit, "ETA", 1L)
    list(
        coef = if (is.matrix(term$beta)) term$beta else term$b,
        pip = term$d
    )
}

# Posterior summaries as plain numbers: a matrix of several conditions' kept
# whole, anything else (hibayes reports one phenotype's as a one-column
# matrix) as a vector.
# @noRd
.posteriorValues <- function(x) {
    if (is.matrix(x) && ncol(x) > 1L) unname(x) else as.numeric(x)
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
    h2,
    methodArgs
) {
    if (!requireNamespace("qgg", quietly = TRUE)) {
        abort("Package 'qgg' is required for backend = \"qgg\".")
    }
    user <- as.list(.alphabetWithH2(
        .alphabetQggPi(methodArgs, pi),
        "QggOptions",
        "qgg",
        h2
    ))
    callArgs <- list_modify(
        list(
            y = y,
            W = X,
            method = .alphabetEngineMethods$qgg[[method]],
            nit = nIter - nBurnin,
            nburn = nBurnin,
            nthin = thin
        ),
        !!!user
    )
    exec(qgg::gbayes, !!!callArgs)
}

# qgg's record with bayesC's prior proportion `pi`, for the methods that
# have one.
# @noRd
.alphabetQggPi <- function(methodArgs, pi) {
    if (is.null(pi)) {
        return(methodArgs)
    }
    .methodOptionsWith(methodArgs, "QggOptions", pi = pi)
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
    h2,
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
        methodArgs = .alphabetWithH2(methodArgs, "BglrOptions", "bglr", h2)
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
    h2,
    methodArgs
) {
    if (!requireNamespace("hibayes", quietly = TRUE)) {
        abort("Package 'hibayes' is required for backend = \"hibayes\".")
    }
    ids <- as.character(seq_len(nrow(X)))
    user <- as.list(.alphabetWithH2(
        .alphabetHibayesPi(methodArgs, pi, "HibayesOptions"),
        "HibayesOptions",
        "hibayes",
        h2,
        var(y)
    ))
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
        !!!user
    )
    withr::with_preserve_seed(exec(hibayes::ibrm, !!!callArgs))
}

# hibayes' record, built by `ctorName`, with the prior proportions `Pi` of
# the null and non-null effects, for the methods that have one.
# @noRd
.alphabetHibayesPi <- function(methodArgs, pi, ctorName) {
    if (is.null(pi)) {
        return(methodArgs)
    }
    .methodOptionsWith(methodArgs, ctorName, Pi = c(1 - pi, pi))
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
    h2Estimate = NULL,
    h2EstimateParam = NULL,
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
        h2Estimate = h2Estimate,
        h2EstimateParam = h2EstimateParam,
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
    h2Estimate = NULL,
    h2EstimateParam = NULL,
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
        h2Estimate = h2Estimate,
        h2EstimateParam = h2EstimateParam,
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
#' @param backend \code{"qgg"} (default), \code{"hibayes"} or
#'   \code{"bglr"}; only \code{"bglr"} takes a matrix \code{y}, fitting its
#'   conditions jointly (see \code{\link{bayesAlphabetWeights}}).
#' @param pi The prior proportion of non-null effects, where each engine's
#'   estimate of it starts. Default \code{0.1}.
#' @param methodArgs Options for the backend's engine:
#'   \code{\link{QggOptions}} for qgg, \code{\link{HibayesOptions}} for
#'   hibayes, \code{\link{BglrOptions}} for BGLR. \code{NULL} keeps the
#'   engine's defaults.
#' @return A numeric vector of weights, one per variant (column of
#'   \code{X}), or for a matrix \code{y} a matrix of them with one column
#'   per condition; columns dropped for zero variance receive weight 0.
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
    backend = c("qgg", "hibayes", "bglr"),
    pi = 0.1,
    h2Estimate = NULL,
    h2EstimateParam = NULL,
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
        h2Estimate = h2Estimate,
        h2EstimateParam = h2EstimateParam,
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
    h2Estimate = NULL,
    h2EstimateParam = NULL,
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
        h2Estimate = h2Estimate,
        h2EstimateParam = h2EstimateParam,
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
#' @param backend \code{"qgg"} (default), \code{"hibayes"} or
#'   \code{"bglr"}; only \code{"bglr"} takes a matrix \code{y}, fitting its
#'   conditions jointly (see \code{\link{bayesAlphabetWeights}}).
#' @param methodArgs Options for the backend's engine:
#'   \code{\link{QggOptions}} for qgg, \code{\link{HibayesOptions}} for
#'   hibayes, \code{\link{BglrOptions}} for BGLR. \code{NULL} keeps the
#'   engine's defaults.
#' @return A numeric vector of weights, one per variant (column of
#'   \code{X}), or for a matrix \code{y} a matrix of them with one column
#'   per condition; columns dropped for zero variance receive weight 0.
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
    backend = c("qgg", "hibayes", "bglr"),
    h2Estimate = NULL,
    h2EstimateParam = NULL,
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
        h2Estimate = h2Estimate,
        h2EstimateParam = h2EstimateParam,
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
    h2Estimate = NULL,
    h2EstimateParam = NULL,
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
        h2Estimate = h2Estimate,
        h2EstimateParam = h2EstimateParam,
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        fitRetention = arg_match(fitRetention),
        methodArgs = methodArgs
    )
}

# --- Bayesian alphabet on summary statistics --------------------------------
#
# qgg and hibayes read the LD as a matrix, formed from the panel for a
# genotype reference (.alphabetRssDense()). It is formed once, and every
# chain reads it.

#' Bayesian Alphabet TWAS Weights From Summary Statistics
#'
#' Fits a Bayesian alphabet prior to marginal summary statistics and an LD
#' reference and returns the posterior mean effects. Each method runs on
#' the LD matrix through one of these engines, chosen by \code{backend}:
#' \itemize{
#'   \item \code{"qgg"}: \code{qgg::gbayes}, for \code{bayesA},
#'     \code{bayesC}, \code{bayesL}, \code{bayesN} and \code{bayesR}.
#'   \item \code{"hibayes"}: \code{hibayes::sbrm}, for every method.
#' }
#' \code{NULL} picks qgg for \code{bayesA}, \code{bayesC}, \code{bayesN}
#' and \code{bayesR}, and hibayes for \code{bayesB} and \code{bayesL}. For a
#' genotype reference the matrix is formed from its genotypes.
#'
#' By default both engines read the statistics through SDPR's likelihood
#' (\code{ldProjection}), given to each in its own terms: qgg takes the
#' diagonal of the LD to be 1 and reads per-variant weights in its place,
#' so it gets the LD at unit diagonal with those weights; hibayes reads
#' the LD as it is.
#'
#' Several chains (\code{nChains}) run the engine once each, from seeds
#' drawn from R's random number generator, and are combined as
#' \code{\link{ldpred2AutoWeights}} combines LDpred2-auto's: a chain whose
#' posterior mean effects' squares sum past twice the squared marginal
#' effects has diverged and is dropped, and of the rest, those whose fitted
#' correlations span less than 0.95 times the 95th percentile of the
#' chains' spans are dropped too; the kept chains' estimates are averaged.
#' The engines keep their draws to themselves, so a chain is tested on its
#' posterior mean. A single chain that diverged is returned with a warning.
#'
#' @param stat A list of summary statistics: \code{n} (per-variant sample
#'   sizes) and the marginal effects as \code{cor} (standardized), \code{z},
#'   \code{b} with \code{seb}, or \code{b} alone, taken as already
#'   standardized. The engines receive the correlation each z stands for,
#'   z / sqrt(z^2 + n - 2), with n the median of \code{n}.
#' @param ldReference The LD reference the statistics' variants are matched
#'   to by \code{stat$variantNames}: a genotype panel (see
#'   \code{\link{readGenotypes}}), a \code{\link{GenotypeMixture}} or an
#'   \code{\link{LdMatrix}}.
#' @param method The prior: \code{"bayesA"}, \code{"bayesB"},
#'   \code{"bayesC"}, \code{"bayesL"}, \code{"bayesN"} or \code{"bayesR"}.
#' @param backend The engine, one the method has (see Description);
#'   \code{NULL} picks the method's default.
#' @param h2Estimate,h2EstimateParam The trait's heritability in the window,
#'   which sets every engine's starting variances and the scales of its
#'   priors: \code{h2Estimate}, a number in [0, 1) or a one-row
#'   \code{\linkS4class{H2Estimate}}; or \code{h2EstimateParam}, the
#'   \code{\link{H2EstimateParam}} settings to estimate it under first, on
#'   \code{stat} and \code{ldReference} (as \code{\link{estimateLocalH2}}
#'   estimates QTL summary statistics), taking the reference to share no
#'   samples with the study. Neither (the default) keeps each engine's own,
#'   0.5; both is an error. It is used floored at 0.001 and capped at 0.99,
#'   and an estimate that cannot be made leaves the engine's default, with a
#'   warning. qgg reads it as \code{h2}; hibayes, which fixes its own at 0.5,
#'   gets the genetic and residual variances (\code{vg}, \code{ve}) its
#'   defaults would give at it. Setting that engine's spelling through
#'   \code{methodArgs} as well is an error.
#' @inheritParams bayesAlphabetWeights
#' @param shrinkCorr The LD shrinkage \eqn{c}: both engines read
#'   \eqn{cR + (1 - c)I}. Default 0.95, the LDpred2 vignette's. A reference
#'   with fewer samples than the window has variants lacks directions the
#'   statistics have, and without shrinkage the samplers' effects can grow
#'   without bound along them. A small or mismatched reference does better
#'   with less, down to 0.4 as the LDpred2 vignette advises; the study's own
#'   genotypes as the reference need none, up to 1.
#' @param ldProjection The ridge \eqn{a} of SDPR's likelihood (Zhou and Zhao
#'   2021), under which both engines read the LD as
#'   \eqn{R (R + aI)^{-1} R} and the marginal effects as
#'   \eqn{R (R + aI)^{-1} \hat\beta}, before the shrinkage. Each of R's
#'   eigen-directions keeps \eqn{\lambda / (\lambda + a)} of the effects, so
#'   their part outside the reference LD's span, which a reference smaller
#'   than the window or unlike the study leaves, never enters; that part
#'   can make the statistics claim more of the trait's variance than any
#'   phenotype has, which sends the samplers' residual variance negative.
#'   Default 0.1, SDPR's; \code{NULL} reads the LD and the statistics as
#'   they are.
#' @param nChains The number of chains, each the engine run once from its
#'   own seed and combined as Description says. Default 1.
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights: \code{"none"} (default) attaches
#'   nothing, \code{"slim"} and \code{"full"} attach the engine's estimates
#'   (\code{betaEst}, and qgg's inclusion probabilities \code{pip}),
#'   \code{backend}, the \code{shrinkCorr} and \code{ldProjection} used and,
#'   with several chains, \code{chains}: each chain's seed, whether it
#'   diverged and whether it was kept. No engine keeps trimmable
#'   intermediates, so the two retaining levels behave alike.
#' @param methodArgs Options for the backend's engine, built with its
#'   constructor: \code{\link{QggOptions}} for qgg,
#'   \code{\link{HibayesRssOptions}} for hibayes. \code{NULL} keeps the
#'   engine's defaults.
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
    h2Estimate = NULL,
    h2EstimateParam = NULL,
    nIter = 5000,
    nBurnin = 1000,
    thin = 5,
    shrinkCorr = 0.95,
    ldProjection = 0.1,
    nChains = 1,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    retention <- arg_match(fitRetention)
    prior <- arg_match(method, names(.alphabetBackends$sumstat))
    engine <- .methodBackend(prior, backend, "sumstat")
    lds <- list(shrinkCorr = shrinkCorr, ldProjection = ldProjection)
    .alphabetRssAssert(stat, nIter, nBurnin, thin, lds, nChains)
    run <- list(
        fitter = .alphabetFitters$sumstat[[engine]],
        method = prior,
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        pi = .alphabetPi(prior, pi),
        methodArgs = .alphabetOptions(methodArgs, "sumstat", prior, engine),
        h2 = .alphabetRssH2(h2Estimate, h2EstimateParam, stat, ldReference),
        solverInput = .alphabetRssDense(stat, ldReference, lds, engine)
    )
    fit <- .mcmcChains(
        ".alphabetRssChain",
        run,
        as.list(sample.int(.Machine$integer.max, nChains)),
        .alphabetRssFields[[engine]],
        "bayesAlphabetRssWeights",
        run$solverInput$bhat
    )
    .alphabetRssResult(
        c(fit, list(backend = engine), lds),
        retention,
        stat[["variantNames"]]
    )
}

# What bayesAlphabetRssWeights() checks before it reads the LD; `lds` holds
# its shrinkCorr and ldProjection.
# @noRd
.alphabetRssAssert <- function(stat, nIter, nBurnin, thin, lds, nChains) {
    .assertMcmcSchedule(nIter, nBurnin, thin)
    .alphabetAssertSampleSize(stat$n, "bayesAlphabetRssWeights")
    .assertShrinkCorr(lds$shrinkCorr)
    assertNumber(
        lds$ldProjection,
        lower = 0,
        finite = TRUE,
        null.ok = TRUE,
        .var.name = "ldProjection"
    )
    assertCount(nChains, positive = TRUE)
    invisible(NULL)
}

# The heritability a summary-statistics fit starts from (.h2Resolve()), an
# estimate being made on its statistics and LD reference.
# @noRd
.alphabetRssH2 <- function(h2Estimate, h2EstimateParam, stat, ldReference) {
    label <- "bayesAlphabetRssWeights"
    .h2Resolve(
        h2Estimate,
        h2EstimateParam,
        label,
        .h2StatFit,
        list(stat = stat, ldReference = ldReference, label = label)
    )
}

# A summary-statistics fit's weights, with the fit (its estimates and the
# settings it carries) and the variants it was fitted on attached unless
# `fitRetention` is "none".
# @noRd
.alphabetRssResult <- function(fit, fitRetention, ids) {
    if (identical(fitRetention, "none")) {
        return(fit$betaEst)
    }
    `attr<-`(fit$betaEst, "fit", .withFitVariantIds(fit, ids))
}

# The matrix engines' standard errors divide by n - 2, so the sample size
# (the median of `n`) must exceed 2.
# @noRd
.alphabetAssertSampleSize <- function(n, label) {
    assertNumeric(n, any.missing = FALSE, min.len = 1L)
    if (median(n) <= 2) {
        abort(glue("{label}: the sample size must exceed 2."))
    }
    invisible(NULL)
}

# What each engine's chain reports, which several chains average.
# @noRd
.alphabetRssFields <- list(qgg = c("betaEst", "pip"), hibayes = "betaEst")

# One chain of a matrix engine: its adapter's fit, with R's random number
# generator seeded by `seed` for the call (qgg draws its own seed from it,
# and hibayes' is drawn from it). The engines keep their draws to
# themselves, so a tracked chain is tested on its posterior mean, against
# the marginal effects as given, and carries the fitted correlations that
# .mcmcKeptChains() filters on (.alphabetFittedCorr()).
# @noRd
.alphabetRssChain <- function(fitter, solverInput, ..., seed, track) {
    fit <- withr::with_seed(
        seed,
        exec(fitter, solverInput = solverInput, ...)
    )
    if (!track) {
        return(fit)
    }
    c(
        fit,
        list(
            corrEst = .alphabetFittedCorr(solverInput, fit$betaEst),
            diverged = .mcmcMeanDiverged(fit$betaEst, solverInput$bhat)
        )
    )
}

# The fitted correlations R b under the LD the likelihood reads: the
# engine's matrix, or for the unit-diagonal matrix qgg reads under the
# projection, that matrix taken back to its own diagonal through `scale`.
# @noRd
.alphabetFittedCorr <- function(solverInput, beta) {
    s <- solverInput$scale %||% 1
    as.vector(s * (solverInput$R %*% (s * beta)))
}

# What the matrix engines read: the LD matrix and the effects they fit, with
# the standard errors that set their y'y; the marginal effects as given
# (`bhat`), which a chain's divergence is tested against; n and the variant
# ids. As given, the LD is formed from the panel for a genotype reference,
# shrunk to c R + (1 - c) I, and both engines read it and the standardized
# effects alike. Under `lds$ldProjection`, they read SDPR's likelihood, each
# in its own terms (.alphabetRssProjected()).
# @noRd
.alphabetRssDense <- function(stat, ldReference, lds, engine) {
    if (!is.null(lds$ldProjection)) {
        return(.alphabetRssProjected(stat, ldReference, lds, engine))
    }
    R <- .rssLdMatrix(stat, ldReference, "bayesAlphabetRssWeights")
    n <- median(stat$n)
    bhat <- .rssMarginalCor(stat, n = n, p = nrow(R))
    .rssValidateSource(bhat, list(R = R), n)
    list(
        R = .rssShrunkMatrix(R, lds$shrinkCorr),
        b = bhat,
        seb = sqrt((1 - bhat^2) / (n - 2)),
        bhat = bhat,
        n = n,
        ids = stat$variantNames
    )
}

# SDPR's likelihood (Zhou and Zhao 2021, its optLlk = 1) for the matrix
# engines: the RSS likelihood with the marginal effects' covariance
# (R + aI) / n, under which a sampler reads R' = R (R + aI)^-1 R as the LD
# and b' = R (R + aI)^-1 bhat as the effects. R' and b' come from the LD
# source's eigendecomposition (.ldSourceEigen(): a genotype reference's
# factor, so its p x p matrix is never needed for it, or an LD matrix), and
# R' is then shrunk to c R' + (1 - c) I. Each engine gets them in its own
# terms (.alphabetProjectedInputs).
# @noRd
.alphabetRssProjected <- function(stat, ldReference, lds, engine) {
    ld <- .rssLdSource(stat, ldReference, "bayesAlphabetRssWeights")
    n <- median(stat$n)
    bhat <- .rssMarginalCor(stat, n = n, p = .ldSourceSize(ld))
    .rssValidateSource(bhat, ld, n)
    projected <- .alphabetProjection(
        .ldSourceEigen(ld),
        bhat,
        lds$ldProjection
    )
    exec(
        .alphabetProjectedInputs[[engine]],
        R = .rssShrunkMatrix(projected$R, lds$shrinkCorr),
        effects = projected$effects,
        bhat = bhat,
        n = n,
        ids = stat$variantNames
    )
}

# R' = V diag(lambda^2 / (lambda + a)) V' and
# b' = V diag(lambda / (lambda + a)) V' bhat from the LD's eigendecomposition
# `ev`, over the directions with a positive eigenvalue: each keeps
# lambda / (lambda + a) of the effects, so their part outside the LD's span
# is dropped.
# @noRd
.alphabetProjection <- function(ev, bhat, a) {
    kept <- ev$values > sqrt(.Machine$double.eps) * max(ev$values)
    lambda <- ev$values[kept]
    V <- ev$vectors[, kept, drop = FALSE]
    list(
        R = V %*% ((lambda^2 / (lambda + a)) * t(V)),
        effects = as.vector(V %*% (lambda / (lambda + a) * crossprod(V, bhat)))
    )
}

# The projected likelihood in each engine's terms.
# @noRd
.alphabetProjectedInputs <- c(
    qgg = ".alphabetProjectedQgg",
    hibayes = ".alphabetProjectedHibayes"
)

# qgg takes the diagonal of X'X to be its per-variant weight ww, reading
# the LD's diagonal as 1, and reads X'y as wy, both in `stat` when given.
# So it gets the LD at unit diagonal, D^-1/2 R D^-1/2 for R's diagonal d,
# with ww = n d and wy = n b': X'X = sqrt(ww ww') LD is then n R exactly,
# and X'y n b'. `scale` takes the unit-diagonal matrix back to R. Its y'y,
# median((b^2 + (n - 2) seb^2) ww), is n through the standard errors, with
# b = b' / d (.alphabetProjectedSe()).
# @noRd
.alphabetProjectedQgg <- function(R, effects, bhat, n, ids) {
    d <- diag(R)
    b <- effects / d
    list(
        R = R / sqrt(outer(d, d)),
        scale = sqrt(d),
        b = b,
        seb = .alphabetProjectedSe(b, d, n),
        ww = n * d,
        wy = n * effects,
        bhat = bhat,
        n = n,
        ids = ids
    )
}

# hibayes forms X'X = n ldm and X'y = n ldm_jj b, so given ldm = R and
# b = b' / d for R's diagonal d it reads n R and n b' exactly. Its y'y, the
# mean of n d (b^2 + (n - 2) se^2), is n through the standard errors
# (.alphabetProjectedSe()).
# @noRd
.alphabetProjectedHibayes <- function(R, effects, bhat, n, ids) {
    d <- diag(R)
    b <- effects / d
    list(
        R = R,
        b = b,
        seb = .alphabetProjectedSe(b, d, n),
        bhat = bhat,
        n = n,
        ids = ids
    )
}

# The standard errors under which the engines' y'y is n, as it is for the
# statistics as given: (b^2 + (n - 2) se^2) d = 1 for each variant. Only
# statistics that claim more of y'y than there is can have b^2 d > 1; such
# a variant gets se 0 and adds b^2 d to y'y's average instead of 1.
# @noRd
.alphabetProjectedSe <- function(b, d, n) {
    sqrt(pmax(1 / d - b^2, 0) / (n - 2))
}

# qgg's summary-statistics sampler, on its in-memory LD: every variant's LD
# with every other, as the dense matrix gives it. qgg runs nit + nburn
# iterations, so nit is nIter net of the burn-in.
# @noRd
.alphabetQggRssFit <- function(
    solverInput,
    method,
    nIter,
    nBurnin,
    thin,
    pi,
    h2,
    methodArgs
) {
    if (!requireNamespace("qgg", quietly = TRUE)) {
        abort("Package 'qgg' is required for backend = \"qgg\".")
    }
    qin <- .alphabetQggInputs(solverInput)
    user <- as.list(.alphabetWithH2(
        .alphabetQggPi(methodArgs, pi),
        "QggOptions",
        "qgg",
        h2
    ))
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
        !!!user
    )
    fit <- exec(qgg::gbayes, !!!callArgs)
    list(
        betaEst = unname(fit$stat[solverInput$ids, "bm"]),
        pip = unname(fit$stat[solverInput$ids, "dm"])
    )
}

# qgg's summary statistics, map and LD for one window, built in memory: the
# statistics on the standardized scale (allele labels and frequencies are
# placeholders qgg carries but does not read), with the per-variant weights
# `ww` and `wy` when the input gives them (.alphabetProjectedQgg()), and
# each variant's LD with every other.
# @noRd
.alphabetQggInputs <- function(inp) {
    m <- length(inp$b)
    ids <- inp$ids
    pos <- seq_len(m)
    list(
        stat = exec(
            data.frame,
            rsids = ids,
            chr = 1L,
            pos = pos,
            ea = "A",
            nea = "G",
            eaf = 0.5,
            b = inp$b,
            seb = inp$seb,
            p = 2 * pnorm(-abs(inp$b / inp$seb)),
            n = inp$n,
            !!!compact(list(ww = inp$ww, wy = inp$wy)),
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
    solverInput,
    method,
    nIter,
    nBurnin,
    thin,
    pi,
    h2,
    methodArgs
) {
    if (!requireNamespace("hibayes", quietly = TRUE)) {
        abort("Package 'hibayes' is required for backend = \"hibayes\".")
    }
    user <- as.list(.alphabetWithH2(
        .alphabetHibayesPi(methodArgs, pi, "HibayesRssOptions"),
        "HibayesRssOptions",
        "hibayes",
        h2,
        .hibayesRssVary(solverInput)
    ))
    callArgs <- list_modify(
        list(
            sumstat = data.frame(
                SNP = solverInput$ids,
                A1 = "A",
                A2 = "G",
                freq = 0.5,
                b = solverInput$b,
                se = solverInput$seb,
                p = 2 * pnorm(-abs(solverInput$b / solverInput$seb)),
                N = solverInput$n
            ),
            ldm = solverInput$R,
            method = .alphabetEngineMethods$hibayes[[method]],
            niter = nIter,
            nburn = nBurnin,
            thin = thin,
            seed = .sessionSeedDraw(),
            threads = 1L,
            verbose = FALSE
        ),
        !!!user
    )
    fit <- withr::with_preserve_seed(exec(hibayes::sbrm, !!!callArgs))
    list(betaEst = as.numeric(fit$alpha))
}

# The phenotypic variance hibayes' sbrm() computes from the statistics it is
# given (src/SBayesD.cpp): the mean over variants of n R_jj (b^2 + (N - 2)
# se^2), over n - 1, n being the mean sample size truncated to an integer.
# @noRd
.hibayesRssVary <- function(inp) {
    n <- as.integer(mean(inp$n))
    mean(diag(inp$R) * n * (inp$b^2 + (inp$n - 2) * inp$seb^2)) / (n - 1)
}

#' @title BayesA TWAS weights from summary statistics
#' @description \code{\link{bayesAlphabetRssWeights}} with the BayesA prior.
#' @inheritParams bayesAlphabetRssWeights
#' @param backend \code{"qgg"} or \code{"hibayes"}; \code{NULL} picks
#'   qgg.
#' @return A numeric vector of the posterior mean effects, one per variant.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' stat <- list(
#'   cor = as.vector(cor(X, eqtlRegionExample$yRes)),
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
    h2Estimate = NULL,
    h2EstimateParam = NULL,
    nIter = 5000,
    nBurnin = 1000,
    thin = 5,
    shrinkCorr = 0.95,
    ldProjection = 0.1,
    nChains = 1,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    bayesAlphabetRssWeights(
        stat,
        ldReference,
        "bayesA",
        backend = backend,
        h2Estimate = h2Estimate,
        h2EstimateParam = h2EstimateParam,
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        shrinkCorr = shrinkCorr,
        ldProjection = ldProjection,
        nChains = nChains,
        fitRetention = fitRetention,
        methodArgs = methodArgs
    )
}

#' @title BayesB TWAS weights from summary statistics
#' @description \code{\link{bayesAlphabetRssWeights}} with the BayesB prior,
#'   which hibayes fits.
#' @inheritParams bayesAlphabetRssWeights
#' @param pi The prior proportion of non-null effects, where hibayes'
#'   estimate of it starts. Default \code{0.2}.
#' @return A numeric vector of the posterior mean effects, one per variant.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' stat <- list(
#'   cor = as.vector(cor(X, eqtlRegionExample$yRes)),
#'   n = rep(nrow(X), ncol(X)),
#'   variantNames = colnames(X)
#' )
#' bayesBRssWeights(stat, LdMatrix(cor(X)), nIter = 1000, nBurnin = 200)
#' @export
bayesBRssWeights <- function(
    stat,
    ldReference,
    pi = 0.2,
    h2Estimate = NULL,
    h2EstimateParam = NULL,
    nIter = 5000,
    nBurnin = 1000,
    thin = 5,
    shrinkCorr = 0.95,
    ldProjection = 0.1,
    nChains = 1,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    bayesAlphabetRssWeights(
        stat,
        ldReference,
        "bayesB",
        pi = pi,
        h2Estimate = h2Estimate,
        h2EstimateParam = h2EstimateParam,
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        shrinkCorr = shrinkCorr,
        ldProjection = ldProjection,
        nChains = nChains,
        fitRetention = fitRetention,
        methodArgs = methodArgs
    )
}

#' @title BayesC TWAS weights from summary statistics
#' @description \code{\link{bayesAlphabetRssWeights}} with the BayesC prior.
#' @inheritParams bayesAlphabetRssWeights
#' @param backend \code{"qgg"} or \code{"hibayes"}; \code{NULL} picks
#'   qgg.
#' @param pi The prior proportion of non-null effects, where each engine's
#'   estimate of it starts. Default \code{0.1}.
#' @return A numeric vector of the posterior mean effects, one per variant.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' stat <- list(
#'   cor = as.vector(cor(X, eqtlRegionExample$yRes)),
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
    h2Estimate = NULL,
    h2EstimateParam = NULL,
    nIter = 5000,
    nBurnin = 1000,
    thin = 5,
    shrinkCorr = 0.95,
    ldProjection = 0.1,
    nChains = 1,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    bayesAlphabetRssWeights(
        stat,
        ldReference,
        "bayesC",
        backend = backend,
        pi = pi,
        h2Estimate = h2Estimate,
        h2EstimateParam = h2EstimateParam,
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        shrinkCorr = shrinkCorr,
        ldProjection = ldProjection,
        nChains = nChains,
        fitRetention = fitRetention,
        methodArgs = methodArgs
    )
}

#' @title BayesL TWAS weights from summary statistics
#' @description \code{\link{bayesAlphabetRssWeights}} with the Bayesian
#'   LASSO.
#' @inheritParams bayesAlphabetRssWeights
#' @param backend \code{"hibayes"} or \code{"qgg"}; \code{NULL} picks
#'   hibayes.
#' @return A numeric vector of the posterior mean effects, one per variant.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' stat <- list(
#'   cor = as.vector(cor(X, eqtlRegionExample$yRes)),
#'   n = rep(nrow(X), ncol(X)),
#'   variantNames = colnames(X)
#' )
#' bayesLRssWeights(stat, LdMatrix(cor(X)), nIter = 1000, nBurnin = 200)
#' @export
bayesLRssWeights <- function(
    stat,
    ldReference,
    backend = NULL,
    h2Estimate = NULL,
    h2EstimateParam = NULL,
    nIter = 5000,
    nBurnin = 1000,
    thin = 5,
    shrinkCorr = 0.95,
    ldProjection = 0.1,
    nChains = 1,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    bayesAlphabetRssWeights(
        stat,
        ldReference,
        "bayesL",
        backend = backend,
        h2Estimate = h2Estimate,
        h2EstimateParam = h2EstimateParam,
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        shrinkCorr = shrinkCorr,
        ldProjection = ldProjection,
        nChains = nChains,
        fitRetention = fitRetention,
        methodArgs = methodArgs
    )
}

#' @title BayesN TWAS weights from summary statistics
#' @description \code{\link{bayesAlphabetRssWeights}} with a normal prior
#'   shared by every effect.
#' @inheritParams bayesAlphabetRssWeights
#' @param backend \code{"qgg"} or \code{"hibayes"}; \code{NULL} picks
#'   qgg.
#' @return A numeric vector of the posterior mean effects, one per variant.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' stat <- list(
#'   cor = as.vector(cor(X, eqtlRegionExample$yRes)),
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
    h2Estimate = NULL,
    h2EstimateParam = NULL,
    nIter = 5000,
    nBurnin = 1000,
    thin = 5,
    shrinkCorr = 0.95,
    ldProjection = 0.1,
    nChains = 1,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    bayesAlphabetRssWeights(
        stat,
        ldReference,
        "bayesN",
        backend = backend,
        h2Estimate = h2Estimate,
        h2EstimateParam = h2EstimateParam,
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        shrinkCorr = shrinkCorr,
        ldProjection = ldProjection,
        nChains = nChains,
        fitRetention = fitRetention,
        methodArgs = methodArgs
    )
}

#' @title BayesR TWAS weights from summary statistics
#' @description \code{\link{bayesAlphabetRssWeights}} with the BayesR prior.
#' @inheritParams bayesAlphabetRssWeights
#' @param backend \code{"qgg"} or \code{"hibayes"}; \code{NULL} picks
#'   qgg.
#' @return A numeric vector of the posterior mean effects, one per variant.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' stat <- list(
#'   cor = as.vector(cor(X, eqtlRegionExample$yRes)),
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
    h2Estimate = NULL,
    h2EstimateParam = NULL,
    nIter = 5000,
    nBurnin = 1000,
    thin = 5,
    shrinkCorr = 0.95,
    ldProjection = 0.1,
    nChains = 1,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    bayesAlphabetRssWeights(
        stat,
        ldReference,
        "bayesR",
        backend = backend,
        h2Estimate = h2Estimate,
        h2EstimateParam = h2EstimateParam,
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        shrinkCorr = shrinkCorr,
        ldProjection = ldProjection,
        nChains = nChains,
        fitRetention = fitRetention,
        methodArgs = methodArgs
    )
}

# --- LDpred2 -----------------------------------------------------------------
#
# LDpred2-auto and LDpred2-grid run on bigsnpr, on the LD matrix written to a
# temporary file, formed from the panel for a genotype reference. Both take
# the statistics in standardized units, the correlations for each variant's
# own sample size, and return weights in the same units.

#' LDpred2-auto Weights from Summary Statistics
#'
#' LDpred2-auto (Prive et al. 2020, 2023) estimates the point-normal prior's
#' proportion of causal variants \eqn{p} and the heritability \eqn{h^2}
#' within its Gibbs sampler, so nothing is tuned. One chain runs from each
#' starting \eqn{p} in \code{pInit}. As in the LDpred2 vignette, a chain
#' whose fitted correlations span less than 0.95 times the 95th percentile of
#' the chains' spans is dropped, and the rest are averaged. A chain has
#' diverged, and is dropped too, once a draw's squared effects sum past
#' twice the squared marginal effects; bigsnpr tests every iteration,
#' burn-in included.
#'
#' The defaults are the LDpred2 vignette's, without the MAF-dependent prior
#' variance.
#'
#' @details
#' bigsnpr's \code{snp_ldpred2_auto()} runs on the LD matrix, written to a
#' temporary file for the call; for a genotype reference the matrix is
#' formed from its genotypes. bigsnpr's functions call
#' \code{foreach::registerDoSEQ()}, which replaces any foreach parallel
#' backend registered in the session. bigsnpr draws its seeds from R's
#' random number generator, so \code{set.seed()} makes the weights
#' reproducible.
#'
#' @param stat A list of summary statistics: \code{n} (per-variant sample
#'   sizes) and the marginal effects as \code{cor} (standardized), \code{z},
#'   \code{b} with \code{seb}, or \code{b} alone, taken as already
#'   standardized. The sampler reads the correlation each z stands for,
#'   z / sqrt(z^2 + n - 2), for each variant's own sample size.
#' @param ldReference The LD reference the statistics' variants are matched
#'   to by \code{stat$variantNames}: a genotype panel (see
#'   \code{\link{readGenotypes}}), a \code{\link{GenotypeMixture}} or an
#'   \code{\link{LdMatrix}}.
#' @param h2Estimate,h2EstimateParam The starting heritability:
#'   \code{h2Estimate}, a number in [0, 1) or a one-row
#'   \code{\linkS4class{H2Estimate}}; or \code{h2EstimateParam}, the
#'   \code{\link{H2EstimateParam}} settings to estimate it under for the
#'   window, as \code{\link{estimateLocalH2}} does for QTL summary
#'   statistics, the reference taken to share no samples with the study.
#'   Neither (the default) estimates it under the default settings, HDL-L's
#'   (\code{H2EstimateParam()}); both is an error. LAVA's and HDL-L's
#'   settings read the reference's sample size, so an \code{LdMatrix} must
#'   record \code{nSamples}. Floored at 0.001, which is also the start where
#'   no estimate can be made.
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
#'   w <- ldpred2AutoWeights(stat, ld, h2Estimate = 0.1,
#'     pInit = c(0.01, 0.1), fitRetention = "slim")
#'   attr(w, "fit")$chains
#' }
#' @export
ldpred2AutoWeights <- function(
    stat,
    ldReference,
    h2Estimate = NULL,
    h2EstimateParam = NULL,
    pInit = exp(seq(log(1e-4), log(0.2), length.out = 30)),
    nIter = 700,
    nBurnin = 500,
    shrinkCorr = 0.95,
    allowJumpSign = FALSE,
    pBounds = c(1e-5, 1),
    fitRetention = c("none", "slim", "full")
) {
    retention <- arg_match(fitRetention)
    label <- "ldpred2AutoWeights"
    .assertMcmcSchedule(nIter, nBurnin, 1L)
    .ldpred2AssertProportions(pInit, "pInit")
    .ldpred2AssertProportions(pBounds, "pBounds")
    .ldpred2AssertSampler(shrinkCorr, allowJumpSign)
    inp <- .rssPerVariantInput(stat, ldReference, label)
    config <- list(
        pInit = pInit,
        h2 = .ldpred2H2(h2Estimate, h2EstimateParam, inp, label),
        nIter = nIter,
        nBurnin = nBurnin,
        shrinkCorr = shrinkCorr,
        allowJumpSign = allowJumpSign,
        pBounds = sort(pBounds)
    )
    engine <- .bigsnprEngine(inp)
    on.exit(.bigsnprRelease(engine), add = TRUE)
    .withWeightsFitIds(
        .ldpred2AutoWeightsFrom(
            .ldpred2RunAuto(inp, engine, config),
            retention,
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
#' grid point has diverged, and is not picked, once a draw's squared effects
#' sum past twice the squared marginal effects; bigsnpr tests every
#' iteration, burn-in included.
#'
#' @details
#' bigsnpr's \code{snp_ldpred2_grid()} runs on the LD matrix, written to a
#' temporary file for the call. Its grid sampler has no LD shrinkage of its
#' own, so the matrix written is shrunk, \eqn{cR + (1 - c)I}, as the
#' Bayesian alphabet's engines are given theirs; it always lets a drawn
#' effect change sign. See \code{\link{ldpred2AutoWeights}} for the
#' temporary file, the foreach backend and the seeds.
#'
#' @inheritParams ldpred2AutoWeights
#' @param p The proportions of causal variants in the grid.
#' @param h2Grid The heritabilities in the grid; \code{NULL} (the default)
#'   makes the grid a centre times 0.3, 0.7, 1 and 1.4, at most 1, as the
#'   LDpred2 vignette scales LD score regression's estimate.
#' @param h2Estimate,h2EstimateParam The grid's centre, as
#'   \code{\link{ldpred2AutoWeights}} takes its start: a value, or the
#'   settings to estimate it under for the window; neither estimates it under
#'   the default settings, HDL-L's. Either with \code{h2Grid} is an error.
#' @param sparse Whether each fit sets an effect whose inclusion probability
#'   falls below \code{p} to zero; both by default.
#' @param shrinkCorr The LD shrinkage \eqn{c}: the matrix bigsnpr reads is
#'   \eqn{cR + (1 - c)I}. Default 0.95, as for
#'   \code{\link{ldpred2AutoWeights}}; bigsnpr's grid itself has none.
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
#'   w <- ldpred2GridWeights(stat, ld, p = c(0.01, 0.1, 1), h2Grid = 0.1)
#'   attr(w, "ldpred2_selection")
#' }
#' @export
ldpred2GridWeights <- function(
    stat,
    ldReference,
    p = exp(seq(log(1e-5), 0, length.out = 21)),
    h2Grid = NULL,
    h2Estimate = NULL,
    h2EstimateParam = NULL,
    sparse = c(FALSE, TRUE),
    nIter = 150,
    nBurnin = 50,
    shrinkCorr = 0.95,
    selection = c("heldOut", "ldQuadratic"),
    fitRetention = c("none", "slim", "full")
) {
    pick <- arg_match(selection)
    retention <- arg_match(fitRetention)
    label <- "ldpred2GridWeights"
    .assertMcmcSchedule(nIter, nBurnin, 1L)
    .ldpred2AssertProportions(p, "p")
    .assertShrinkCorr(shrinkCorr)
    inp <- .rssPerVariantInput(stat, ldReference, label)
    if (identical(pick, "heldOut")) {
        .rssCheckHeldOutLd(inp$ld, label)
    }
    h2 <- .ldpred2GridH2(h2Grid, h2Estimate, h2EstimateParam, inp, label)
    config <- list(
        grid = .ldpred2Grid(p, h2, sparse),
        nIter = nIter,
        nBurnin = nBurnin
    )
    engine <- .bigsnprEngine(inp, shrinkCorr)
    on.exit(.bigsnprRelease(engine), add = TRUE)
    .withWeightsFitIds(
        .ldpred2GridWeightsFrom(
            .ldpred2RunGrid(inp, engine, config),
            list(inp = inp, engine = engine, config = config),
            pick,
            retention,
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

# The statistics bigsnpr's engines read: the LD source (.rssLdSource()), the
# correlations (.rssMarginalCor()) and each variant's sample size.
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
# as bigsnpr floors it, else the window's estimate under the caller's
# settings, the default ones (HDL-L's) when none are given.
# @noRd
.ldpred2H2 <- function(h2Estimate, h2EstimateParam, inp, label) {
    .h2AssertOneSource(h2Estimate, h2EstimateParam, label)
    if (is.null(h2Estimate)) {
        settings <- h2EstimateParam %||% H2EstimateParam()
        return(.ldpred2EstimateH2(settings, inp, label))
    }
    max(.h2SuppliedValue(h2Estimate, label), 1e-3)
}

# The heritabilities of the grid: the caller's, else a centre (.ldpred2H2())
# times 0.3, 0.7, 1 and 1.4, at most 1.
# @noRd
.ldpred2GridH2 <- function(h2Grid, h2Estimate, h2EstimateParam, inp, label) {
    if (is.null(h2Grid)) {
        centre <- .ldpred2H2(h2Estimate, h2EstimateParam, inp, label)
        return(pmin(centre * c(0.3, 0.7, 1, 1.4), 1))
    }
    if (!is.null(h2Estimate) || !is.null(h2EstimateParam)) {
        abort(glue(
            "{label}: give the grid (`h2Grid`) or its centre (`h2Estimate` ",
            "or `h2EstimateParam`), not both."
        ))
    }
    .ldpred2AssertProportions(h2Grid, "h2Grid")
    h2Grid
}

# The window's heritability under the settings `param`, the reference taken
# to share no samples with the study, and floored at bigsnpr's 0.001: the
# start too where no estimate can be made.
# @noRd
.ldpred2EstimateH2 <- function(param, inp, label) {
    settings <- .h2EstimateSettings(param, FALSE, label)
    if (.h2ReadsReferenceSize(settings) && is.na(inp$ld$nSamples)) {
        abort(c(
            glue(
                "{label}: estimating the heritability under these settings ",
                "needs the LD reference's sample size."
            ),
            i = "Build the LdMatrix with `nSamples`, or give the heritability."
        ))
    }
    h2 <- .h2LocalFit(
        settings,
        z = .corToZ(inp$bhat, inp$n),
        nVariant = inp$n,
        ld = inp$ld
    )$h2
    if (is.na(h2)) 1e-3 else max(h2, 1e-3)
}

# Every combination of p, h2 and sparse, large p first as bigsnpr runs them.
# @noRd
.ldpred2Grid <- function(p, h2, sparse) {
    assertLogical(sparse, any.missing = FALSE, min.len = 1L)
    grid <- expand.grid(p = p, h2 = h2, sparse = sparse)
    as_tibble(grid[order(-grid$p, grid$sparse, -grid$h2), ])
}

# bigsnpr's LD: the matrix as bigsnpr's sparse matrix, written to a temporary
# file that .bigsnprRelease() deletes; for a genotype reference, formed from
# its factor; shrunk to `shrinkCorr` (.rssShrunkMatrix()) for a sampler of
# bigsnpr's with no shrinkage of its own (the grid's).
# @noRd
.bigsnprEngine <- function(inp, shrinkCorr = 1) {
    if (!requireNamespace("bigsnpr", quietly = TRUE)) {
        abort(str_c(
            "Package 'bigsnpr' is required for LDpred2 and for ",
            "backend = \"bigsnpr\"."
        ))
    }
    R <- .rssShrunkMatrix(
        inp$ld$R %||%
            (crossprod(inp$ld$ldFactor) / (nrow(inp$ld$ldFactor) - 1)),
        shrinkCorr
    )
    corr <- bigsnpr::as_SFBM(
        methods::as(Matrix::Matrix(R, sparse = TRUE), "generalMatrix"),
        tempfile(),
        compact = TRUE
    )
    list(corr = corr)
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
    withr::with_seed(
        .sessionSeedDraw(),
        .ldpred2BigsnprAuto(inp, engine$corr, config)
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

# The grid's fits, variants x points, NA for a point that diverged. bigsnpr
# seeds from R's generator, so it runs under a seed drawn from it.
# @noRd
.ldpred2RunGrid <- function(inp, engine, config) {
    withr::with_seed(
        .sessionSeedDraw(),
        bigsnpr::snp_ldpred2_grid(
            engine$corr,
            .bigsnprDfBeta(inp),
            as.data.frame(config$grid),
            burn_in = config$nBurnin,
            num_iter = config$nIter - config$nBurnin,
            ncores = 1
        )
    )
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
# draws are (.mcmcMeanDiverged()). The draws themselves are not tested,
# since one noisy draw passes that bound on some null windows while the
# chain's mean stays near zero.
# @noRd
.mcmcWarnDiverged <- function(beta, bhat, label) {
    if (.mcmcMeanDiverged(beta, bhat)) {
        warn(glue(
            "{label}: the posterior mean effects' squares sum past twice ",
            "the squared marginal effects, so the chain diverged and its ",
            "effects are unreliable; several chains (`nChains`) drop the ",
            "chains that diverge."
        ))
    }
    invisible(NULL)
}

# Whether effects `beta` have diverged: their squares summing past twice the
# marginal effects' `bhat`, LDpred2's bound.
# @noRd
.mcmcMeanDiverged <- function(beta, bhat) {
    sum(beta^2) > 2 * sum(bhat^2)
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
    keptBetas <- map(chains[kept], "betaEst")
    beta <- rowMeans(exec(cbind, !!!keptBetas))
    if (identical(fitRetention, "none")) {
        return(beta)
    }
    keptPostps <- map(chains[kept], "postp")
    fit <- list(
        chains = tibble(
            pInit = map_dbl(chains, "pInit"),
            p = map_dbl(chains, "p"),
            h2 = map_dbl(chains, "h2"),
            diverged = map_lgl(chains, "diverged"),
            kept = kept
        ),
        postp = rowMeans(exec(cbind, !!!keptPostps)),
        p = mean(map_dbl(chains[kept], "p")),
        h2 = mean(map_dbl(chains[kept], "h2"))
    )
    if (identical(fitRetention, "full")) {
        chainBetas <- map(chains, "betaEst")
        fit$beta <- exec(cbind, !!!chainBetas)
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
#' summary statistics", Genetic Epidemiology 41(6):469-480. It is
#' \code{\link{penalizedRss}}'s LASSO at \code{alpha = 1}: lassosum's
#' coordinate updates, walked down the path by ncvreg's strong-rule
#' strategy unless \code{strategy} says otherwise.
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
#' @param lambda A vector of L1 penalty values. \code{NULL} (default) is
#'   \code{\link{penalizedRss}}'s default path: 20 values on a log scale from
#'   the smallest at which every coefficient is zero down to ncvreg's
#'   \code{lambda.min}, 0.1\% of it with more samples than variants and 5\%
#'   otherwise.
#' @param thr Convergence threshold. Default: 1e-4.
#' @param maxiter Maximum number of iterations. Default: 10000.
#' @param strategy How the path is walked; see \code{\link{penalizedRss}}.
#'   \code{NULL} (default) is \code{"strongRule"}.
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
    lambda = NULL,
    thr = 1e-4,
    maxiter = 10000,
    strategy = NULL
) {
    .penalizedRssPath(
        bhat,
        list(R = R),
        n,
        spec = .rssPenalty("lasso"),
        lambda = lambda,
        thr = thr,
        maxiter = maxiter,
        strategy = strategy
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

# The correlation behind a t statistic z on n samples, the marginal effect of a
# standardized variant on a standardized trait. z / sqrt(n) is only close to it
# while z^2 is small against n, which a GWAS can assume and a strong QTL cannot.
# @noRd
.zToCor <- function(z, n) {
    z / sqrt(z^2 + n - 2)
}

# The standardized marginal effects, the scale lassosum, PRS-CS and SDPR all
# fit: the correlations (.zToCor()) rather than their usual stand-in z /
# sqrt(n), which past r^2 = 0.5 exceeds 1 and fails a strong QTL's fit. Read
# from `cor` as given, else from `z`, else from `b / seb`, which is a z-score
# whatever scale `b` is on, else from `b` alone, taken as already standardized.
.rssMarginalCor <- function(stat, n, p) {
    corInput <- if (!is.null(stat$cor)) {
        as.numeric(stat$cor)
    } else if (!is.null(stat$z)) {
        .zToCor(as.numeric(stat$z), n)
    } else if (!is.null(stat$b) && !is.null(stat$seb)) {
        .zToCor(as.numeric(stat$b) / as.numeric(stat$seb), n)
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
        bxxb = as.numeric(
            (1 - w) *
                .ldQuadForms(ld, candidateBeta) +
                w * colSums(candidateBeta^2)
        )
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

# Selection on held-out summary statistics, with PUMAS's splits (Zhao et al.
# 2021). For the marginal correlations c of n samples, a training part of
# nTrain of them has c_t ~ N(c, (n - nTrain) / (n nTrain) Sigma), and its
# held-out remainder c_h = (n c - nTrain c_t) / (n - nTrain). Sigma is the
# LD shrunk toward the identity by the Schafer-Strimmer intensity for the
# reference's size: a reference smaller than the number of variants misses
# directions of the true LD, the least-penalized fits put weight exactly
# there, and draws and scores from R alone never see it. Each split refits
# the whole (s, lambda) grid on c_t and scores it on c_h; the full-data
# candidate with the best mean score is returned (.rssHeldOutPick()). One
# whose full-data fit is empty has nothing to return, so it counts as having
# no correlation.
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
    corScores <- map(scores, "cor")
    explainedScores <- map(scores, "explained")
    idx <- .rssHeldOutPick(
        replace(
            rowMeans(exec(cbind, !!!corScores)),
            colSums(candidateBeta != 0, na.rm = TRUE) == 0,
            -Inf
        ),
        rowMeans(exec(cbind, !!!explainedScores))
    )
    list(beta = candidateBeta[, idx], index = idx, mode = "heldOut")
}

# Split `k`: the candidates refit on its training part and scored on its
# held-out part, by the variance of it they explain, 2 c_h' b - b' Sigma b,
# and by their correlation with it, which is -Inf for a fit empty there.
# `refit` names the refitting function (`fn`) and its other arguments
# (`args`); it is given the training part's correlations and the share of the
# samples that part keeps, and returns the candidates in the full data's
# order.
# @noRd
.rssHeldOutScore <- function(k, noise, corInput, ld, refit, w) {
    share <- .rssHeldOutTrainShare()
    trainCor <- .lassosumClampCor(corInput + noise[, k])
    heldCor <- (corInput - share * trainCor) / (1 - share)
    refitArgs <- refit$args
    beta <- exec(refit$fn, trainCor = trainCor, share = share, !!!refitArgs)
    moments <- .rssShrunkLdMoments(beta, heldCor, ld, w)
    list(
        cor = .rssCorrelationScore(moments),
        explained = 2 * moments$bxy - moments$bxxb
    )
}

# The candidate explaining the most held-out variance on average. PUMAS picks
# on the correlation instead, but that ignores scale: at s near 1 the fit nearly
# ignores the LD, every tag of a signal keeps its marginal weight, and it can
# win with a predicted variance thousands of times the trait's. A correlation
# of -Inf marks a fit empty in some split or on the full data. Such a candidate
# is out while any other is left, or its score of 0 would beat every fit that
# overshoots, so a fit is returned wherever one exists. A tie in the variance
# goes to the better correlation, then to the first in the grid.
# @noRd
.rssHeldOutPick <- function(cor, explained) {
    ties <- .rssScoreTies(replace(explained, !is.finite(cor), -Inf))
    ties[.lassosumFirstMax(cor[ties])]
}

# The shrinkage grid refit on a held-out split's training part, at the
# full-data fit's path values `grids`, so its candidates line up with the
# full data's.
# @noRd
.rssRefitShrinkGrid <- function(trainCor, share, s, ld, n, config, grids) {
    .rssFitShrinkGrid(
        s,
        trainCor * sqrt(share * n),
        ld,
        share * n,
        config,
        grids
    )$beta
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
    assertVector(bhat, len = ncol(ld$ldFactor))
    invisible(NULL)
}

# Validate (bhat, ld, n), run a decreasing-lambda coordinate-descent sweep via
# `solveFn` (the penalty-specific C++ solver), then reorder back to the input
# lambda order via the inverse permutation and assemble the standard result
# list. Shared by lassosumRss and penalizedRss, which differ only in which C++
# solver they pass as `solveFn` (and penalizedRss's per-penalty gamma default).
# `ld` is an LD source: the matrix, or a factor of it with the shrinkage the
# factor solver applies itself. An empty `lambda` leaves the path to the
# solver (the L0 family's default), returned in the order it solved it.
.rssSolvePath <- function(bhat, ld, n, lambda, solveFn, solveArgs = list()) {
    .rssValidateSource(bhat, ld, n)
    z <- bhat / sqrt(n)
    order <- order(lambda, decreasing = TRUE)
    # `solveArgs` rather than `...`: the two solvers take different fixed
    # argument sets, and both callers know theirs statically, so an unknown
    # solver argument should be an error here rather than reaching the solver.
    solved <- exec(solveFn, z, lambda[order], ld, !!!solveArgs)
    ownPath <- length(lambda) == 0L
    # Reorder back to original lambda order via the inverse permutation.
    invOrder <- if (ownPath) seq_along(solved$lambda) else order(order)
    beta <- solved$beta[, invOrder, drop = FALSE]
    fbeta <- solved$fbeta[invOrder]
    list_assign(
        solved,
        beta = beta,
        conv = solved$conv[invOrder],
        loss = solved$loss[invOrder],
        fbeta = fbeta,
        lambda = if (ownPath) as.numeric(solved$lambda) else lambda,
        nparams = as.integer(colSums(beta != 0)),
        betaEst = as.numeric(beta[, which.min(fbeta)])
    )
}

# Per-`s` fit for one RSS method: one penalized path at the wrapper's
# penalty, or, for l0learn, one lambda0 path per value of its L1 grid.
# `grid` is the path to solve at, a held-out refit's full-data values (NULL
# for the full-data fit itself). Returns list(beta, meta, grid), `grid` the
# path values solved at, for the refits.
# @noRd
.rssFitOne <- function(solverInput, ld, n, sVal, config, grid = NULL) {
    if (identical(config$method, "l0learn")) {
        return(.rssFitL0learn(solverInput, ld, n, sVal, config, grid))
    }
    .rssFitPenalized(solverInput, ld, n, sVal, config, grid)
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
# penalty, `config$spec`, under the caller's path options, at `grid` when
# given.
# @noRd
.rssFitPenalized <- function(solverInput, ld, n, sVal, config, grid = NULL) {
    model <- .penalizedRssOn(
        solverInput,
        ld,
        n,
        config$spec,
        .rssGridOptions(config$options, grid)
    )
    list(
        beta = model$beta,
        meta = .rssPathMeta(sVal, model),
        grid = model$lambda
    )
}

# The path options of one fit: the caller's, with `path` as the path when
# given.
# @noRd
.rssGridOptions <- function(options, path) {
    if (is.null(path)) {
        return(options)
    }
    list_assign(options, lambda = path)
}

# L0Learn path for one s: one warm-started lambda0 path per value of the L1
# grid, column-bound. `grid`, when given, holds each L1 value's lambda0 path.
# @noRd
.rssFitL0learn <- function(solverInput, ld, n, sVal, config, grid = NULL) {
    fits <- map2(
        config$lambda1,
        grid %||% vector("list", length(config$lambda1)),
        .rssL0PathFit,
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
        meta = bind_rows(metaList),
        grid = map(fits, "grid")
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
        list(ldFactor = .ldFactorFromRef(ldReference, ids, label))
    } else {
        list(R = .ldFromRef(ldReference, ids, label))
    }
    list_assign(ld, nSamples = .ldRefNoiseN(ldReference))
}

# The LD matrix itself, for the engines that take nothing else (mr.ash.rss,
# mr.mash.rss), formed over the variants in the order `perm` when given, so
# that a permuted fit needs no second p x p copy.
# @noRd
.rssLdMatrix <- function(stat, ldReference, label, perm = NULL) {
    ids <- .rssStatVariantIds(stat, label)
    .ldFromRef(ldReference, ids[perm %||% seq_along(ids)], label)
}

# The LD a fit taking one LD per condition (mr.mash.rss, mvsusie_rss) is
# given, as list(R, X, n, group): a factor of it (X) when every reference is a
# genotype reference, so the p x p matrix is never formed, or else the
# matrix (R). One reference's LD is over the statistics' variants (n and group
# NULL, left to the caller); for a list of one reference per condition, there
# is one LD per reference over the variants its conditions carry
# (.rssGroupLd()), each condition the median sample size of the conditions
# of its sample group, which the fit takes as one sample, and that group
# (.rssSampleGroups()): the conditions on one reference, split by
# `stat$sampleGroups` when given (one label per condition, such as its
# study), as studies that share a reference panel are on different samples.
# One reference shared by several sample groups is taken as one per
# condition.
# @noRd
#' @importFrom stats ave median
.rssConditionLd <- function(stat, ldReference, label) {
    Z <- as.matrix(stat$z)
    perCondition <- .jgIsPerConditionLd(ldReference)
    if (!perCondition && length(unique(stat$sampleGroups)) <= 1L) {
        return(.rssOneLd(stat, ldReference, label))
    }
    refs <- if (perCondition) {
        ldReference
    } else {
        rep(list(ldReference), ncol(Z))
    }
    if (length(refs) != ncol(Z)) {
        abort(glue(
            "{label}: a list `ldReference` holds one LD reference ",
            "per column of stat$z ({ncol(Z)}); got {length(refs)}."
        ))
    }
    refGroups <- .rssRefGroups(refs)
    useFactor <- all(map_lgl(refs, .ldIsGenotypeRef))
    ld <- map(
        seq_len(max(refGroups)),
        .rssGroupLd,
        groups = refGroups,
        refs = refs,
        Z = Z,
        ids = .rssStatVariantIds(stat, label),
        label = label,
        useFactor = useFactor
    )
    groups <- .rssSampleGroups(refGroups, stat$sampleGroups)
    n <- ave(rep_len(as.numeric(stat$n), ncol(Z)), groups, FUN = median)
    if (useFactor) {
        list(X = ld[refGroups], n = n, group = groups)
    } else {
        list(R = ld[refGroups], n = n, group = groups)
    }
}

# The groups of conditions a fit takes as measured on the same samples: the
# conditions on one LD reference (`refGroups`), split by `sampleGroups` (one
# label per condition) when given. Numbered in order of first appearance.
# @noRd
.rssSampleGroups <- function(refGroups, sampleGroups) {
    if (is.null(sampleGroups)) {
        return(refGroups)
    }
    key <- paste(refGroups, sampleGroups, sep = "\r")
    match(key, unique(key))
}

# The residual variance each condition's heritability implies (`h2`, one
# per condition, NA for none), as SuSiE-RSS's 1 - h2: that fraction of the
# condition's variance (`variance`, 1 on the standardized scale), which the
# multivariate fits give their residual correlations. NULL without `h2`.
# @noRd
.rssResidualVariances <- function(h2, variance = 1) {
    if (is.null(h2)) {
        return(NULL)
    }
    (1 - replace(h2, is.na(h2), 0)) * variance
}

# `V` with zero between conditions on different samples (`group`, as
# .rssConditionLd() gives it: on different LD references, or in different
# sample groups); as it is when they share their samples (`group` NULL).
# @noRd
.rssMaskBetweenRefs <- function(V, group) {
    if (is.null(V) || is.null(group)) {
        return(V)
    }
    replace(V, outer(group, group, "!="), 0)
}

# One reference's LD over the statistics' variants, as .rssConditionLd()
# gives it: its factor for a genotype reference, else its matrix.
# @noRd
.rssOneLd <- function(stat, ldReference, label) {
    if (!.ldIsGenotypeRef(ldReference)) {
        return(list(R = .rssLdMatrix(stat, ldReference, label), n = NULL))
    }
    ids <- .rssStatVariantIds(stat, label)
    list(X = .ldFactorFromRef(ldReference, ids, label), n = NULL)
}

# Each condition's reference group: conditions with one (identical) reference
# share a group, numbered in order of first appearance.
# @noRd
#' @importFrom purrr detect_index
.rssRefGroups <- function(refs) {
    first <- map_int(seq_along(refs), .rssFirstRef, refs = refs)
    match(first, unique(first))
}

# @noRd
.rssFirstRef <- function(i, refs) {
    detect_index(refs, identical, refs[[i]])
}

# The LD of the conditions in group `g`: their reference's over the variants
# they carry that it holds, laid on the statistics' variants with zero rows
# and columns (and diagonal) for the rest, which carry no information in
# them. With `useFactor`, its factor, with zero columns for the rest.
# @noRd
.rssGroupLd <- function(g, groups, refs, Z, ids, label, useFactor = FALSE) {
    cols <- which(groups == g)
    carried <- .rssGroupCarried(cols, Z, ids)
    fromRef <- if (useFactor) .ldFactorFromRef else .ldFromRef
    ld <- if (length(carried) > 0L) {
        fromRef(refs[[cols[[1L]]]], carried, label, onMissing = "drop")
    }
    if (useFactor) .rssEmbedFactor(ld, ids) else .rssEmbedLd(ld, ids)
}

# The variants the conditions `cols` carry: those with a statistic in at
# least one of them.
# @noRd
.rssGroupCarried <- function(cols, Z, ids) {
    ids[rowSums(!is.na(Z[, cols, drop = FALSE])) > 0L]
}

# The panels of the conditions in group `g` when their reference is a
# GenotypeMixture -- each panel's LD, or with `useFactor` its factor
# (.ldPanelsFromRef()) -- laid on the statistics' variants as .rssGroupLd()
# lays the reference's LD; NULL for any other reference.
# @noRd
.rssGroupPanels <- function(
    g,
    groups,
    refs,
    Z,
    ids,
    label,
    useFactor = FALSE
) {
    cols <- which(groups == g)
    ref <- refs[[cols[[1L]]]]
    if (!.ldIsMixtureRef(ref)) {
        return(NULL)
    }
    carried <- .rssGroupCarried(cols, Z, ids)
    panels <- if (length(carried) > 0L) {
        .ldPanelsFromRef(ref, carried, label, "drop", useFactor = useFactor)
    }
    map(
        panels %||% vector("list", length(.gmComponents(ref))),
        if (useFactor) .rssEmbedFactor else .rssEmbedLd,
        ids = ids
    )
}

# A factor `ldFactor` of LD laid on the variants `ids`: a zero column where
# it holds none (one row of zeros when it holds no variant at all), which
# gives that variant a zero row and column of the LD.
# @noRd
#' @importFrom dplyr coalesce
.rssEmbedFactor <- function(ldFactor, ids) {
    if (is.null(ldFactor)) {
        return(matrix(0, 1L, length(ids), dimnames = list(NULL, ids)))
    }
    padded <- cbind(ldFactor, 0)
    at <- coalesce(match(ids, colnames(ldFactor)), ncol(padded))
    `dimnames<-`(padded[, at, drop = FALSE], list(NULL, ids))
}

# `ld` laid on the variants `ids`: zero where it holds no row.
# @noRd
.rssEmbedLd <- function(ld, ids) {
    if (is.null(ld)) {
        return(matrix(0, length(ids), length(ids), dimnames = list(ids, ids)))
    }
    at <- match(ids, rownames(ld))
    embedded <- unname(as.matrix(ld)[at, at, drop = FALSE])
    `dimnames<-`(replace(embedded, is.na(embedded), 0), list(ids, ids))
}

# The LD a fit's credible sets' purity reads, as list(R, X): `ld`'s, or,
# when the fit estimated a mixture's weights (`omega`) from its panels
# (`ld$panels`), the panels mixed with them -- from their factors when
# there are, as for the block's own.
# @noRd
.rssPurityLd <- function(ld, omega) {
    if (is.null(ld$panels) || is.null(omega)) {
        return(list(R = ld$R, X = ld$X))
    }
    if (!is.null(ld$panels$X)) {
        return(list(X = .rssWeightedFactor(ld$panels$X, omega)))
    }
    list(R = .rssWeightedLd(ld$panels$R, omega))
}

# Whether a condition's LD is a list of panels, whose weights the fit
# estimates, rather than one matrix or factor.
# @noRd
.rssIsPanels <- function(ld) {
    is.list(ld) && !is.matrix(ld)
}

# A condition's LD mixed from its panels with the weights `omega`
# (.rssWeightedFactor(), or .rssWeightedLd() without `useFactor`); its LD as
# it is when it has no panels.
# @noRd
.rssMixedPanels <- function(ld, omega, useFactor) {
    if (!.rssIsPanels(ld)) {
        return(ld)
    }
    if (useFactor) .rssWeightedFactor(ld, omega) else .rssWeightedLd(ld, omega)
}

# The LD of panels mixed with the weights `omega`: the weighted sum of
# their matrices.
# @noRd
.rssWeightedLd <- function(panels, omega) {
    reduce(map2(panels, omega, `*`), `+`)
}

# A factor (.ldRefFactor()'s contract) of the panels, given by their
# factors, mixed with the weights `omega`: each panel's rows scaled so that
# its crossproduct is omega_k times its LD, stacked, and rescaled for the
# rows of the stack.
# @noRd
.rssWeightedFactor <- function(panels, omega) {
    stacked <- exec(rbind, !!!map2(panels, omega, .rssWeightedFactorBlock))
    stacked * sqrt(nrow(stacked) - 1)
}

# @noRd
.rssWeightedFactorBlock <- function(ldFactor, w) {
    ldFactor * sqrt(w / (nrow(ldFactor) - 1))
}

# The statistics as a fit on its LD (.rssConditionLd()) is given them: a
# variant a condition lacks -- no statistic there, or, with one LD per
# condition, no LD, its reference not holding it -- has z and Bhat 0 and
# Shat 1, carrying no information there.
# @noRd
.rssCarriedStat <- function(stat, LD) {
    lacked <- if (is.list(LD$R) || is.list(LD$X)) {
        .rssLackedVariants(LD)
    } else {
        FALSE
    }
    held <- intersect(c("z", "Bhat", "Shat"), names(compact(stat)))
    list_assign(
        stat,
        !!!set_names(
            map(held, .rssFillLacked, stat = stat, lacked = lacked),
            held
        )
    )
}

# The z-scores an RSS fit reads (`LD` as .rssConditionLd() gives it): 0 for
# a variant a condition's LD lacks, which carries no information there, and
# missing where the condition has no z-score for a variant its LD holds.
# The engines fit a missing z-score as 0 but leave it out of their estimate
# of the conditions' residual correlation, which a 0 would pull toward 0.
# @noRd
.rssFitZ <- function(stat, LD) {
    z <- as.matrix(stat$z)
    if (!is.list(LD$R) && !is.list(LD$X)) {
        return(z)
    }
    replace(z, .rssLackedVariants(LD), 0)
}

# `expr`, an RSS fit given missing z-scores (.rssFitZ()), without the
# warning its engine raises on fitting them as 0: they are missing by
# design, a condition with no statistic for a variant. Nothing else is
# muffled.
# @noRd
.rssQuietMissingZ <- function(expr) {
    withCallingHandlers(expr, warning = .rssMuffleMissingZ)
}

# @noRd
.rssMuffleMissingZ <- function(w) {
    if (str_detect(conditionMessage(w), "z-scores are replaced with 0")) {
        invokeRestart("muffleWarning")
    }
}

# The variants each condition's LD lacks (variants x conditions), from the
# per-condition LD of .rssConditionLd(): a zero diagonal of its matrix, or a
# zero column of its factor.
# @noRd
.rssLackedVariants <- function(LD) {
    diagonals <- if (is.null(LD$X)) {
        map(LD$R, diag)
    } else {
        map(LD$X, .rssFactorDiagonal)
    }
    matrix(list_c(diagonals), ncol = length(diagonals)) == 0
}

# The diagonal of the LD a factor gives, up to its scale.
# @noRd
.rssFactorDiagonal <- function(ldFactor) {
    colSums(ldFactor^2)
}

# @noRd
.rssFillLacked <- function(field, stat, lacked) {
    values <- as.matrix(stat[[field]])
    replace(
        values,
        is.na(values) | lacked,
        if (identical(field, "Shat")) 1 else 0
    )
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
# A list of one reference per condition holds one per column of `stat$z`,
# each a valid reference; the variants a condition's reference lacks are
# missing there, as the fit had them (.rssConditionLd()).
# @noRd
.rssVariantCount <- function(stat, ldReference, label) {
    if (.jgIsPerConditionLd(ldReference)) {
        nCond <- NCOL(stat$z)
        if (length(ldReference) != nCond) {
            abort(glue(
                "{label}: a list `ldReference` holds one LD reference ",
                "per column of stat$z ({nCond}); got {length(ldReference)}."
            ))
        }
        walk(ldReference, .ldFromRefValidate, label = label)
        return(length(.rssStatVariantIds(stat, label)))
    }
    if (!is.null(ldReference)) {
        ids <- .rssStatVariantIds(stat, label)
        .ldFromRefValidate(ldReference, label)
        .ldFromRefMatch(ldReference, ids, label, "error")
        return(length(ids))
    }
    described <- stat$variantNames %||% stat$z %||% stat$b
    if (is.null(described)) NA_integer_ else NROW(described)
}

# One shrinkage level's fit against R_s = (1 - s) R + s I, at `grid` when
# given.
# @noRd
.rssFitAtS <- function(sVal, grid, solverInput, ld, n, config) {
    .rssFitOne(solverInput, .rssShrunkLd(ld, sVal), n, sVal, config, grid)
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
    mode <- arg_match(selection)
    retention <- arg_match(fitRetention)
    .rssCheckStrategy(config)
    if (mode == "heldOut") {
        .rssCheckHeldOutLd(ld, config$label)
    }
    solverLd <- .rssSolverSource(ld, config$spec)
    n <- median(stat$n)
    corInput <- .lassosumClampCor(
        .rssMarginalCor(stat, n = n, p = .ldSourceSize(solverLd))
    )
    search <- .rssFitShrinkGrid(s, corInput * sqrt(n), solverLd, n, config)
    picked <- .rssShrinkGridSelect(
        mode,
        search,
        corInput,
        solverLd,
        n,
        s,
        config
    )
    weights <- .rssFinalize(
        as.numeric(picked$beta),
        picked,
        search$meta,
        config
    )
    .rssAttachGridFit(
        weights,
        search$beta,
        search$meta,
        picked$index,
        retention,
        stat[["variantNames"]]
    )
}

# A strategy the caller's options name is checked before the grid starts,
# so a misspelled one fails here rather than inside the first fit.
# @noRd
.rssCheckStrategy <- function(config) {
    .rssStrategy(config$options$strategy, config$spec$penalty)
    invisible(NULL)
}

# The grid over s: each s's path against the LD shrunk to it, at `grids`,
# the path values a held-out refit takes from the full-data fit (NULL for
# that fit, which solves at the caller's or the default path). Returns the
# candidates' coefficients (`beta`), their path table (`meta`) and the path
# values each s solved at (`grids`).
# @noRd
.rssFitShrinkGrid <- function(s, solverInput, ld, n, config, grids = NULL) {
    fits <- map2(
        s,
        grids %||% vector("list", length(s)),
        .rssFitAtS,
        solverInput = solverInput,
        ld = ld,
        n = n,
        config = config
    )
    fitBetas <- map(fits, "beta")
    list(
        beta = exec(cbind, !!!fitBetas),
        meta = bind_rows(map(fits, "meta")),
        grids = map(fits, "grid")
    )
}

# The candidate `selection` picks from the grid `search`
# (.rssFitShrinkGrid()): by held-out prediction, refitting the grid on a
# split of the statistics at the full-data path values; by the LD quadratic
# form; or by the solver's smallest F-beta.
# @noRd
.rssShrinkGridSelect <- function(
    selection,
    search,
    corInput,
    ld,
    n,
    s,
    config
) {
    switch(
        selection,
        heldOut = .rssSelectHeldOut(
            search$beta,
            corInput,
            ld,
            n,
            list(
                fn = ".rssRefitShrinkGrid",
                args = list(
                    s = s,
                    ld = ld,
                    n = n,
                    config = config,
                    grids = search$grids
                )
            )
        ),
        ldQuadratic = .lassosumSelectLdQuadratic(search$beta, corInput, ld),
        minFbeta = .lassosumSelectMinFbeta(search$beta, search$meta)
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
#' 2021). The path is refit on each training part, at the penalty values
#' the full data took, and every candidate is scored on the held-out part by
#' the variance of it that it explains,
#' \deqn{2 c_h^T \beta - \beta^T \hat\Sigma \beta, \quad
#'   \hat\Sigma = (1 - w) R + w I,}
#' and the full-data fit at the (\code{s}, lambda) with the best mean score
#' is returned. \eqn{w} is the Schafer-Strimmer intensity for a correlation
#' matrix from the reference's sample size. A reference smaller than the
#' number of variants misses directions of the true LD, and the
#' least-penalized fits put weight exactly there, so scoring against
#' \eqn{R} itself would keep choosing them. The splits use R's random number
#' generator, so \code{set.seed()} makes the weights reproducible.
#'
#' PUMAS scores by the correlation with the held-out part,
#' \eqn{c_h^T \beta / \sqrt{\beta^T \hat\Sigma \beta}}, which here only breaks
#' ties. It does not change with \eqn{\beta}'s scale, so it lets a fit win at
#' \eqn{s} near 1, where the LD is nearly ignored and every tag of a signal
#' keeps its marginal weight, with a predicted variance thousands of times the
#' trait's. A candidate whose fit is empty in a split or on the full data is
#' returned only if every candidate's is, and the first in the grid wins a tie
#' on both scores.
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
#'   standardized. The engine receives the correlation each z stands for,
#'   z / sqrt(z^2 + n - 2), with n the median of \code{n}.
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
    backend <- .methodBackend("lasso", backend, "sumstat")
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
    lassosumOptions <- config$options
    fits <- exec(
        bigsnpr::snp_lassosum2,
        engine$corr,
        .bigsnprDfBeta(inp),
        delta = config$delta,
        ncores = 1,
        !!!lassosumOptions
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
#' On standardized data the objective is exactly twice ncvreg's,
#' \eqn{\frac{1}{2n}\|y - X\beta\|^2 + P(\beta)}, and twice L0Learn's,
#' \eqn{\frac12\|r\|^2 + \lambda_0\|\beta\|_0 + \lambda_1\|\beta\|_1 +
#' \lambda_2\|\beta\|_2^2}, so a penalty value means here what it means
#' there.
#'
#' @section Path strategies:
#' The path is warm-started from its largest value down. Values at or above
#' the smallest at which every coefficient is zero are recorded without being
#' solved, and each later one is solved by one of two strategies:
#' \describe{
#'   \item{\code{"strongRule"}}{ncvreg's: the sequential strong rule screens
#'     each step; coordinate descent cycles, in index order, over every
#'     coordinate that has been active; then scans of the strong set, and of
#'     the rest, add any coordinate that leaves zero. The default for LASSO,
#'     MCP and SCAD.}
#'   \item{\code{"gradientRanked"}}{L0Learn's: each step cycles coordinates
#'     in decreasing size of their gradient at the previous solution,
#'     screened to the support and the 1000 next, sweeps only the support
#'     once it has held for two sweeps, and checks every zero coordinate at
#'     convergence. Unlike L0Learn, the support is always swept and the
#'     ranking always uses the current gradient. The default for the L0
#'     family.}
#' }
#' For LASSO, the elastic net and the ridge the strategies reach the same
#' solution; for MCP, SCAD and L0 they can reach different local optima.
#' Either stops a solve when no coefficient moved by \code{thr} or more in a
#' sweep, ncvreg's test; L0Learn's relative change in the objective is not
#' used.
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
#' @param lambda For LASSO, MCP and SCAD, the path: the penalty strengths to
#'   trace the solution along (warm-started, largest first). \code{NULL}
#'   (default) is ncvreg's path: 20 values on a log scale from
#'   \eqn{\lambda_{max} = \max_j |z_j| / \alpha}, the smallest at which every
#'   coefficient is zero, down to ncvreg's \code{lambda.min}:
#'   \eqn{10^{-3} \lambda_{max}} when the \code{n} samples outnumber the
#'   variants and \eqn{0.05 \lambda_{max}} otherwise, where the fits below it
#'   are near-unpenalized regressions on more variants than samples. A ridge,
#'   with no such value, is sized as glmnet sizes one, as if \code{alpha}
#'   were 0.001, and keeps the full three decades. For the L0 family, the
#'   fixed L1 weight \eqn{\lambda_1} (one value; \code{NULL} is 0).
#' @param gamma Concavity parameter for MCP (default 3) or SCAD (default 3.7).
#'   Ignored for LASSO and L0 variants.
#' @param alpha Elastic-net mixing for LASSO, MCP and SCAD, in \code{[0, 1]}:
#'   \eqn{l_1 = \lambda \alpha} and \eqn{l_2 = \lambda (1-\alpha)}, so the
#'   penalty is \eqn{2 P(\beta; l_1) + l_2 \beta^T \beta}, glmnet's mixing.
#'   Default 1, no ridge term. For LASSO, 0.5 is the elastic net and 0 ridge;
#'   MCP and SCAD need \code{alpha > 0}. Ignored by the L0 variants.
#' @param lambda0 For the L0 family, the path: the L0 weights to trace the
#'   solution along (warm-started, largest first). \code{NULL} (default) is
#'   the solver's own path: under \code{"gradientRanked"}, L0Learn's
#'   adaptive one, each step letting the next coordinate enter, until 100 are
#'   non-zero or 100 points are taken, solved at the geometric midpoints of
#'   its values, since a held-out refit at the values themselves, each just
#'   below a coordinate's entry threshold, often takes a different support;
#'   20 values on a log scale
#'   over three decades below the largest useful value under
#'   \code{"strongRule"}. Ignored otherwise.
#' @param lambda2 L2 penalty weight for L0L2 variant. Default 0.
#' @param thr Convergence threshold: a solve stops when no coefficient moved by
#'   \code{thr} or more in a sweep. Default 1e-4. The ridge, solved exactly,
#'   does not use it.
#' @param maxiter Maximum coordinate descent sweeps per solve. Default 10000.
#'   The scans for coordinates leaving zero do not count. The ridge, solved
#'   exactly, does not use it.
#' @param maxSwaps Maximum swap rounds for L0 variants. Default 100. Set to 0 to
#'   disable swaps.
#' @param strategy How the path is walked: \code{"strongRule"} or
#'   \code{"gradientRanked"} (see Path strategies). \code{NULL} (default) is
#'   \code{"strongRule"} for LASSO, MCP and SCAD and \code{"gradientRanked"}
#'   for the L0 family.
#'
#' @return A list with components:
#' \describe{
#'   \item{beta}{p x (path length) matrix of coefficient estimates.}
#'   \item{lambda}{The path values used: lambda, or lambda0 for the L0
#'     family.}
#'   \item{conv}{Convergence indicators: 1 when every solve at that point,
#'     including those after swaps, converged within \code{maxiter} sweeps.}
#'   \item{loss}{Quadratic loss at each point.}
#'   \item{fbeta}{Full penalized objective at each point.}
#'   \item{nparams}{Number of non-zero coefficients at each point.}
#'   \item{betaEst}{Coefficient vector at the point minimizing fbeta.}
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
    lambda0 = NULL,
    lambda2 = 0,
    thr = 1e-4,
    maxiter = 10000,
    maxSwaps = 100,
    strategy = NULL
) {
    pen <- arg_match(penalty)
    l0 <- .rssIsL0(pen)
    .penalizedRssPath(
        bhat,
        list(R = R),
        n,
        spec = .rssPenalty(
            pen,
            alpha = alpha,
            gamma = gamma,
            lambda1 = if (l0) lambda %||% 0 else 0,
            lambda2 = lambda2,
            maxSwaps = maxSwaps
        ),
        lambda = if (l0) lambda0 else lambda,
        thr = thr,
        maxiter = maxiter,
        strategy = strategy
    )
}

# The penalty a fit minimizes, as one checked record: the solver's penalty
# name, the elastic-net mixing `alpha`, the MCP / SCAD concavity `gamma`
# (their usual default when NULL) and the L0 family's fixed L1 and L2
# weights, `lambda1` and `lambda2`. Each weight wrapper builds it from its
# own arguments before its shrinkage grid starts, so the penalty is the
# wrapper's alone and only the path (lambda, thr, maxiter, strategy) is left
# to the caller's options.
# @noRd
#' @importFrom checkmate assertCount assertNumber
.rssPenalty <- function(
    penalty,
    alpha = 1,
    gamma = NULL,
    lambda1 = 0,
    lambda2 = 0,
    maxSwaps = 0
) {
    .rssCheckAlpha(penalty, alpha)
    assertNumber(gamma, finite = TRUE, null.ok = TRUE)
    assertNumber(lambda1, lower = 0, finite = TRUE)
    assertNumber(lambda2, lower = 0, finite = TRUE)
    assertCount(maxSwaps)
    list(
        penalty = penalty,
        alpha = alpha,
        gamma = gamma %||% switch(penalty, SCAD = 3.7, MCP = 3.0, 0.0),
        lambda1 = lambda1,
        lambda2 = lambda2,
        maxSwaps = maxSwaps
    )
}

# Whether `penalty` is one of the L0 family, whose path is over lambda0.
# @noRd
.rssIsL0 <- function(penalty) {
    is_in(penalty, c("L0", "L0L1", "L0L2"))
}

# The path strategy a penalty runs under: the caller's, else strong-rule
# screening for LASSO, MCP and SCAD and gradient ranking for the L0 family,
# each its reference implementation's (ncvreg's, L0Learn's).
# @noRd
.rssStrategy <- function(strategy, penalty) {
    if (is.null(strategy)) {
        return(if (.rssIsL0(penalty)) "gradientRanked" else "strongRule")
    }
    arg_match0(
        strategy,
        c("strongRule", "gradientRanked"),
        arg_nm = "strategy"
    )
}

# The constants the path solver works to, in one place for R and C++: a
# default path's `points` and its last value over its first, `ratio`, and
# for a lambda path on fewer samples than variants, `ratioFewSamples`
# (ncvreg's lambda.min, 0.001 and 0.05; .rssDefaultLambda()); for gradient
# ranking, the `screenSize` coordinates swept outside the support and the
# `stableSweeps` the support must hold before only it is swept; and for the
# adaptive lambda0 path, its `maxSupport` non-zeros and `maxPoints` points,
# L0Learn's defaults (screenSize, activeSetNum - 1, maxSuppSize, nLambda),
# and `offEdge`, whether the path is solved at the midpoints of its values
# rather than at the full data's entry thresholds, where held-out refits
# score a different support than the one returned.
# @noRd
.rssPathControl <- function() {
    list(
        points = 20L,
        ratio = 1e-3,
        ratioFewSamples = 0.05,
        screenSize = 1000L,
        stableSweeps = 2L,
        maxSupport = 100L,
        maxPoints = 100L,
        offEdge = TRUE
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
# penalty (.rssPenalty()); the rest is the path -- in lambda, or lambda0 for
# the L0 family -- and how it is walked (.rssStrategy()).
# @noRd
.penalizedRssPath <- function(
    bhat,
    ld,
    n,
    spec,
    lambda,
    thr,
    maxiter,
    strategy = NULL
) {
    assertNumeric(lambda, lower = 0, any.missing = FALSE, null.ok = TRUE)
    assertNumber(thr, lower = 0, finite = TRUE)
    assertCount(maxiter, positive = TRUE)
    walk <- .rssStrategy(strategy, spec$penalty)
    # The inputs are checked before the default path reads them.
    .rssValidateSource(bhat, ld, n)
    path <- lambda %||% .rssDefaultPath(bhat / sqrt(n), n, spec)
    .rssSolvePath(
        bhat,
        ld,
        n,
        path,
        .rssPenalizedSolve,
        solveArgs = c(
            spec,
            list(thr = thr, maxiter = maxiter, strategy = walk)
        )
    )
}

# The default path for the marginal correlations `z` of `n` samples: none
# for the L0 family, whose solver then takes its own (.rssSolvePath()), and
# the default lambda path otherwise.
# @noRd
.rssDefaultPath <- function(z, n, spec) {
    if (.rssIsL0(spec$penalty)) {
        return(numeric(0))
    }
    .rssDefaultLambda(z, spec$alpha, n)
}

# The default lambda path, ncvreg's: from lambda_max = max|z| / alpha, the
# smallest lambda at which every coefficient is zero, down to its
# lambda.min in 20 values evenly on a log scale (.rssPathControl()). That is
# 0.1% of lambda_max when the `n` samples outnumber the variants, and 5%
# otherwise, where the fits below it are near-unpenalized regressions on
# more variants than samples, which held-out selection does not pick and
# which take most of the solver's time. A ridge (alpha = 0) has no such
# lambda, so it is sized as glmnet sizes its own, as if alpha were 0.001,
# and keeps the full depth: none of its fits is sparse.
# @noRd
.rssDefaultLambda <- function(z, alpha, n) {
    control <- .rssPathControl()
    ratio <- if (alpha == 0 || n > length(z)) {
        control$ratio
    } else {
        control$ratioFewSamples
    }
    lambdaMax <- max(abs(z), 0) / max(alpha, 1e-3)
    lambdaMax * ratio^seq(0, 1, length.out = control$points)
}

# .penalizedRssPath() for the penalty `spec` under the caller's path
# `options` -- lambda, thr, maxiter and strategy, all an options bag may
# hold -- penalizedRss()'s own defaults filling the rest. A matrix takes
# exactly the path penalizedRss() does.
# @noRd
.penalizedRssOn <- function(bhat, ld, n, spec, options) {
    path <- list_modify(
        .formalDefaults(
            penalizedRss,
            c("lambda", "thr", "maxiter", "strategy")
        ),
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
#'   standardized. The engine receives the correlation each z stands for,
#'   z / sqrt(z^2 + n - 2), with n the median of \code{n}.
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
#'   standardized. The engine receives the correlation each z stands for,
#'   z / sqrt(z^2 + n - 2), with n the median of \code{n}.
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
#' default lambda path is ncvreg's: 20 values on a log scale from
#' \eqn{\max_j |z_j| / \alpha}, the smallest at which every weight is zero,
#' down to 0.1\% of it with more samples than variants and 5\% otherwise,
#' so the L1 share spans the same range at any \code{alpha}.
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
#' default lambda path is sized as glmnet sizes a ridge's, as if
#' \code{alpha} were 0.001: 20 values on a log scale from
#' \eqn{1000 \max_j |z_j|} down to \eqn{\max_j |z_j|}. Searches over a
#' shrinkage grid \code{s} as well, whose own ridge term
#' \eqn{s \beta^T \beta} adds to it; model selection uses held-out summary
#' statistics by default (see \code{\link{lassosumRssWeights}}).
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
#'   standardized. The engine receives the correlation each z stands for,
#'   z / sqrt(z^2 + n - 2), with n the median of \code{n}.
#' @param ldReference The LD reference the statistics' variants are matched
#'   to by \code{stat$variantNames}: a genotype panel (see
#'   \code{\link{readGenotypes}}), a \code{\link{GenotypeMixture}} or an
#'   \code{\link{LdMatrix}}. A genotype reference's LD is used through its
#'   standardized genotypes, never formed as a matrix.
#' @param penalty L0 variant: \code{"L0"}, \code{"L0L1"}, or \code{"L0L2"}.
#'   Default \code{"L0"}.
#' @param s Numeric vector of LD shrinkage parameters. Default: \code{c(0.2,
#'   0.5, 0.9, 1.0)}.
#' @param lambda0 Numeric vector of L0 penalty values to search over, solved
#'   as one warm-started path from the largest down. \code{NULL} (default) is
#'   the solver's own path for each shrinkage level and L1 value: L0Learn's
#'   adaptive one, from the smallest value at which every coefficient is zero,
#'   each step letting the next variant enter, until 100 are non-zero or 100
#'   values are taken, solved at the geometric midpoints of its values (see
#'   \code{\link{penalizedRss}}). Held-out refits solve at the values the
#'   full data took.
#' @param lambda Numeric vector of L1 penalty values (for L0L1), each with its
#'   own lambda0 path. Default: ten values from 1e-4 to 0.1 on a log scale for
#'   L0L1, and \code{0} (no L1) otherwise.
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
    lambda0 = NULL,
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
    pen <- arg_match(penalty)
    mode <- arg_match(selection)
    retention <- arg_match(fitRetention)
    # The L0 penalty; each L1 value of the grid gets its own lambda0 path.
    spec <- .rssPenalty(pen, lambda2 = lambda2, maxSwaps = maxSwaps)
    assertNumeric(
        lambda0,
        lower = 0,
        finite = TRUE,
        any.missing = FALSE,
        null.ok = TRUE
    )
    .rssShrinkGridWeights(
        stat,
        .rssLdSource(stat, ldReference, "l0learnRssWeights"),
        s,
        list(
            method = "l0learn",
            spec = spec,
            lambda0 = lambda0,
            lambda1 = .l0learnRssLambda(lambda, pen),
            options = as.list(methodArgs),
            label = "l0learnRssWeights"
        ),
        mode,
        retention
    )
}

# One s's lambda0 path at the fixed L1 weight `l1Val`, warm-started from the
# largest lambda0 down: at `path` when given (a held-out refit's full-data
# values), else the caller's lambda0 grid, else the solver's own path.
# @noRd
.rssL0PathFit <- function(l1Val, path, solverInput, ld, n, sVal, config) {
    model <- .penalizedRssOn(
        solverInput,
        ld,
        n,
        list_assign(config$spec, lambda1 = l1Val),
        .rssGridOptions(config$options, path %||% config$lambda0)
    )
    meta <- mutate(
        .rssPathMeta(sVal, model),
        lambda0 = .data$lambda,
        lambda = l1Val,
        .after = "s"
    )
    list(beta = model$beta, meta = meta, grid = model$lambda)
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
#'   Each names the variants it was fitted on.
#' @param retainX Logical. Keep that copy of \code{X} in a \code{"full"} fit
#'   (default \code{FALSE}). An error at any other \code{fitRetention}.
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
    retainX = FALSE,
    methodArgs = NcvregOptions()
) {
    .assertMethodOptions(methodArgs, "NcvregOptions", "methodArgs")
    assertCount(nfolds, positive = TRUE)
    fitRetention <- arg_match(fitRetention)
    .assertRetainX(retainX, fitRetention)
    if (!requireNamespace("ncvreg", quietly = TRUE)) {
        abort("Package 'ncvreg' is required for this function.")
    }
    eff.wgt <- matrix(0, ncol = 1, nrow = ncol(X))
    keep <- .dropZeroVariance(X, "ncvregWeights")
    user <- as.list(methodArgs)
    callArgs <- list_modify(
        list(
            X = X[, keep, drop = FALSE],
            y = y,
            penalty = penalty,
            nfolds = nfolds
        ),
        !!!user
    )
    fit <- exec(ncvreg::cv.ncvreg, !!!callArgs)
    coefs <- as.numeric(coef(fit, lambda = fit$lambda.min))
    weights <- replace(eff.wgt, keep, coefs[-1L])
    if (identical(fitRetention, "none")) {
        return(weights)
    }
    kept <- switch(
        fitRetention,
        slim = .ncvregSlimFit(fit, coefs),
        full = if (retainX) {
            fit
        } else {
            `[[<-`(fit, "fit", value = .withoutFields(fit$fit, "X"))
        }
    )
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
#' Fits a SCAD-penalized linear regression and returns the coefficient vector
#' at \code{lambda.min}: on ncvreg (\code{ncvreg::cv.ncvreg}, the default) or
#' on grpreg (\code{backend = "grpreg"}), whose group SCAD with each variant
#' its own group is SCAD. Both take \code{gamma = 3.7} unless
#' \code{methodArgs} sets it (grpreg's own default would be 4).
#'
#' On grpreg a matrix \code{y}, a column per condition, is fitted across the
#' conditions at once, each variant's effects across them one group under the
#' penalty, as \code{\link{lassoWeights}} describes; ncvreg fits one
#' phenotype.
#'
#' @param X A numeric matrix of predictors.
#' @param y A numeric response vector, or on grpreg a matrix with a column
#'   per condition.
#' @param nfolds Number of cross-validation folds. Default is 5.
#' @param backend The engine: \code{"ncvreg"} (the default for \code{NULL})
#'   or \code{"grpreg"}.
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights: \code{"none"} (default) attaches
#'   nothing; \code{"slim"} the selected fit, as a list of its coefficients
#'   (\code{coef}, over the variants fitted, a column per condition for a
#'   matrix \code{y}), \code{intercept}, \code{lambda} and cross-validated
#'   error there (\code{cve}, \code{cvse}); \code{"full"} the whole
#'   \code{cv.ncvreg} or \code{cv.grpreg} object, with every
#'   \code{lambda}'s coefficients but without its copies of the data. Each
#'   names the variants it was fitted on.
#' @param retainX Logical. Keep the fit's copy of \code{X} in a
#'   \code{"full"} fit (default \code{FALSE}): ncvreg's, or on grpreg its
#'   standardized design (its \code{returnX}). An error at any other
#'   \code{fitRetention}.
#' @param methodArgs Options for the backend's engine, built with its
#'   constructor: \code{\link{NcvregOptions}} for ncvreg,
#'   \code{\link{GrpregOptions}} for grpreg. \code{NULL} (the default) keeps
#'   the engine's defaults. A bare list is refused: it cannot be checked
#'   against the engine, so a misspelled option would be silently ignored.
#' @return A numeric vector of length `ncol(X)` of variant weights, or for a
#'   matrix \code{y} a matrix of them with a column per condition, with the
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
    backend = NULL,
    fitRetention = c("none", "slim", "full"),
    retainX = FALSE,
    methodArgs = NULL
) {
    .ncvFamilyWeights(
        X,
        y,
        "SCAD",
        "scad",
        nfolds,
        backend,
        arg_match(fitRetention),
        retainX,
        methodArgs
    )
}

#' Compute Weights Using MCP-Penalized Regression
#'
#' Fits an MCP-penalized linear regression and returns the coefficient vector
#' at \code{lambda.min}: on ncvreg (\code{ncvreg::cv.ncvreg}, the default) or
#' on grpreg (\code{backend = "grpreg"}), whose group MCP with each variant
#' its own group is MCP.
#'
#' On grpreg a matrix \code{y}, a column per condition, is fitted across the
#' conditions at once, each variant's effects across them one group under the
#' penalty, as \code{\link{lassoWeights}} describes; ncvreg fits one
#' phenotype.
#'
#' @inheritParams scadWeights
#' @return A numeric vector of length `ncol(X)` of variant weights, or for a
#'   matrix \code{y} a matrix of them with a column per condition, with the
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
    backend = NULL,
    fitRetention = c("none", "slim", "full"),
    retainX = FALSE,
    methodArgs = NULL
) {
    .ncvFamilyWeights(
        X,
        y,
        "MCP",
        "mcp",
        nfolds,
        backend,
        arg_match(fitRetention),
        retainX,
        methodArgs
    )
}

# MCP or SCAD (`penalty`, method `method`) on its backend: ncvreg
# (ncvregWeights()), or grpreg's group penalty (grpregWeights()), one group
# per variant for one phenotype.
# @noRd
.ncvFamilyWeights <- function(
    X,
    y,
    penalty,
    method,
    nfolds,
    backend,
    retention,
    retainX,
    methodArgs
) {
    engine <- .methodBackend(method, backend, "individual")
    options <- .backendOptions(
        methodArgs,
        .methodBackends$individual[[method]][[engine]],
        engine
    )
    if (identical(engine, "grpreg")) {
        return(grpregWeights(
            X,
            y,
            str_c("gr", penalty),
            nfolds = nfolds,
            fitRetention = retention,
            retainX = retainX,
            methodArgs = options
        ))
    }
    .assertOnePhenotype(y, "ncvreg", "grpreg")
    ncvregWeights(X, y, penalty, nfolds, retention, retainX, options)
}

# Refuse several conditions' phenotypes on an engine that fits one, naming
# the backend that fits them together.
# @noRd
.assertOnePhenotype <- function(y, engine, alternative) {
    if (.isConditionMatrix(y)) {
        abort(glue(
            "{engine} fits one phenotype at a time; a matrix `y`, a column ",
            "per condition, takes backend = \"{alternative}\"."
        ))
    }
    invisible(NULL)
}

#' Compute Weights Using an L0 Penalty
#'
#' Fits an L0-regularized linear regression, a best subset of the variants,
#' and returns the coefficient vector of the fit minimizing the
#' cross-validation error: on L0Learn (\code{L0Learn::L0Learn.cvfit}, the
#' default), at the (lambda, gamma) pair minimizing it, or on abess
#' (\code{backend = "abess"}), at the support size minimizing it. The default
#' penalty is "L0"; on L0Learn it can be "L0L1" or "L0L2" instead (with the
#' gamma grid tuned through \code{methodArgs}), and abess fits "L0" alone.
#'
#' On abess a matrix \code{y}, a column per condition, is fitted across the
#' conditions at once: a best subset of variants, each in the model for every
#' condition or for none, each condition with its own effects, on its
#' phenotype scaled to unit variance and fitted on its own samples. One
#' support size serves every condition, the one with the lowest
#' cross-validated error summed over them, on folds that hold each sample out
#' of every condition it is in. Conditions measured on the same samples run
#' abess's \code{"mgaussian"} family. Conditions measured on samples of their
#' own, a phenotype missing wherever a condition was not measured, run a best
#' subset of groups on a sparse design in which each condition's rows hold
#' its own samples against its own copy of the variants (see
#' \code{\link{lassoWeights}}). abess runs on one thread here, and R's random
#' number generator, which abess reseeds as it returns, is restored.
#'
#' @param X A numeric matrix of predictors.
#' @param y A numeric response vector, or on abess a matrix with a column per
#'   condition.
#' @param penalty Type of regularization: "L0", "L0L1", or "L0L2" on L0Learn;
#'   "L0" on abess. Default is "L0".
#' @param nFolds Number of cross-validation folds. Default is 5.
#' @param backend The engine: \code{"l0learn"} (the default for \code{NULL})
#'   or \code{"abess"}.
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights: \code{"none"} (default) attaches
#'   nothing; \code{"slim"} the selected fit, as a list of its coefficients
#'   (\code{coef}, over the variants fitted, a column per condition for a
#'   matrix \code{y}), \code{intercept}, and on L0Learn \code{gamma},
#'   \code{lambda} and the cross-validated error there (\code{cvMean},
#'   \code{cvSd}), on abess the support size (\code{supportSize}) and the
#'   cross-validated error there (\code{cvMean}); \code{"full"} the whole
#'   \code{L0Learn.cvfit} or \code{abess} object, with every candidate's
#'   coefficients. Either names the variants it was fitted on, in
#'   \code{variantIds}.
#' @param methodArgs Options for the backend's engine, built with its
#'   constructor: \code{\link{L0learnOptions}} for L0Learn,
#'   \code{\link{AbessOptions}} for abess. \code{NULL} (the default) keeps
#'   the engine's defaults. A bare list is refused: it cannot be checked
#'   against the engine, so a misspelled option would be silently ignored.
#' @return A numeric vector of length `ncol(X)` of variant weights, or for a
#'   matrix \code{y} a matrix of them with a column per condition, with the
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
    backend = NULL,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    retention <- arg_match(fitRetention)
    engine <- .methodBackend("l0learn", backend, "individual")
    options <- .backendOptions(
        methodArgs,
        .methodBackends$individual$l0learn[[engine]],
        engine
    )
    if (identical(engine, "abess")) {
        .abessRequireL0(penalty)
        return(abessWeights(X, y, nFolds, retention, options))
    }
    .assertOnePhenotype(y, "L0Learn", "abess")
    .l0learnFit(X, y, penalty, nFolds, retention, options)
}

# abess fits the L0 penalty alone.
# @noRd
.abessRequireL0 <- function(penalty) {
    if (!identical(penalty, "L0")) {
        abort(glue(
            "abess fits the L0 penalty alone, not {penalty}; L0L1 and L0L2 ",
            "take backend = \"l0learn\"."
        ))
    }
    invisible(NULL)
}

# L0Learn's cross-validated fit (l0learnWeights()): the coefficients at the
# (gamma, lambda) pair minimizing the cross-validation error.
# @noRd
.l0learnFit <- function(X, y, penalty, nFolds, retention, methodArgs) {
    if (!requireNamespace("L0Learn", quietly = TRUE)) {
        abort("Package 'L0Learn' is required for this function.")
    }
    eff.wgt <- matrix(0, ncol = 1, nrow = ncol(X))
    keep <- .dropZeroVariance(X, "l0learnWeights")
    user <- as.list(methodArgs)
    callArgs <- list_modify(
        list(
            x = X[, keep, drop = FALSE],
            y = y,
            penalty = penalty,
            nFolds = nFolds
        ),
        !!!user
    )
    fit <- exec(L0Learn::L0Learn.cvfit, !!!callArgs)
    best <- .l0learnBest(fit)
    raw <- as.numeric(coef(fit, lambda = best$lambda, gamma = best$gamma))
    # If intercept was included, it is the first row.
    hasIntercept <- length(raw) == sum(keep) + 1L
    coefs <- if (hasIntercept) raw[-1L] else raw
    weights <- replace(eff.wgt, keep, coefs)
    if (identical(retention, "none")) {
        return(weights)
    }
    kept <- if (identical(retention, "full")) {
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

# ---- Fits across conditions (a group penalty) ------------------------------
#
# A matrix `y` holds a column per condition, missing wherever a condition was
# not measured. Fitted across them, each variant's effects in the conditions
# are one group under the penalty, so a variant is in the model for every
# condition or for none, each condition with its own effect:
#
#   sum_c RSS_c / (2 n_c) + lambda sum_j P(b_j1, ..., b_jC)
#
# Each condition's loss is its own mean squared error, on its phenotype
# scaled to unit variance and fitted on its own samples. glmnet's mgaussian
# family fits conditions measured on the same samples
# (.glmnetGroupedWeights()); grpreg and abess fit any, through the group
# design (.groupedDesign()). One lambda, or for abess one support size, serves
# every condition: the one with the lowest cross-validated error summed over
# them, on folds that hold a sample out of every condition it is in.

#' Compute Weights Using grpreg
#'
#' The engine behind the \code{"grpreg"} backend of
#' \code{\link{lassoWeights}}, \code{\link{enetWeights}},
#' \code{\link{mcpWeights}} and \code{\link{scadWeights}}: fits
#' \code{grpreg::cv.grpreg} and returns the coefficients at
#' \code{lambda.min}. For a vector \code{y} each variant is its own group,
#' which makes grpreg's group lasso, MCP and SCAD the lasso, MCP and SCAD on
#' standardized columns. For a matrix \code{y}, a column per condition, each
#' variant's effects across the conditions form one group on the group
#' design, so a condition may be missing on some samples.
#'
#' @param X A numeric matrix of predictors.
#' @param y A numeric response vector, or a matrix with a column per
#'   condition, missing where a condition was not measured.
#' @param penalty grpreg's penalty: \code{"grLasso"}, \code{"grMCP"} or
#'   \code{"grSCAD"}.
#' @param alpha grpreg's mixing of its penalty with ridge, in (0, 1].
#'   Default 1; \code{methodArgs} may set it for MCP and SCAD.
#' @param nfolds Number of cross-validation folds. Default is 5.
#' @param methodArgs Options forwarded to \code{grpreg::cv.grpreg} and
#'   \code{grpreg::grpreg}, built with \code{\link{GrpregOptions}}.
#' @inheritParams scadWeights
#' @return A numeric matrix of variant weights, a row per column of \code{X}
#'   and a column per condition, with the fit attached as attribute
#'   \code{"fit"} unless \code{fitRetention = "none"}.
#' @importFrom checkmate assertCount
#' @keywords internal
grpregWeights <- function(
    X,
    y,
    penalty,
    alpha = 1,
    nfolds = 5,
    fitRetention = c("none", "slim", "full"),
    retainX = FALSE,
    methodArgs = GrpregOptions()
) {
    .assertMethodOptions(methodArgs, "GrpregOptions", "methodArgs")
    assertCount(nfolds, positive = TRUE)
    retention <- arg_match(fitRetention)
    .assertRetainX(retainX, retention)
    if (!requireNamespace("grpreg", quietly = TRUE)) {
        abort("Package 'grpreg' is required for this function.")
    }
    keep <- .dropZeroVariance(X, "grpregWeights")
    run <- list(
        penalty = penalty,
        alpha = alpha,
        nfolds = nfolds,
        retainX = retainX,
        user = as.list(methodArgs)
    )
    fitted <- if (.isConditionMatrix(y)) {
        .grpregGroupedFit(X[, keep, drop = FALSE], y, run)
    } else {
        .grpregFit(X[, keep, drop = FALSE], y, run)
    }
    .engineWeights(fitted, keep, X, y, retention, .grpregWithoutData)
}

# grpreg's fit of one phenotype, each variant its own group.
# @noRd
.grpregFit <- function(X, y, run) {
    cvFit <- .grpregCv(
        list(X = X, y = y, group = seq_len(ncol(X)), nfolds = run$nfolds),
        run
    )
    coefs <- as.numeric(coef(cvFit))
    list(
        coef = coefs[-1L],
        slim = .grpregSlimFit(cvFit, coefs[-1L], coefs[[1L]]),
        object = cvFit
    )
}

# grpreg's fit across conditions on the group design, a variant's copies one
# group (.groupedDesign()).
# @noRd
.grpregGroupedFit <- function(X, Y, run) {
    design <- .groupedDesign(X, Y, run$nfolds, sparse = FALSE)
    cvFit <- .grpregCv(
        list(
            X = design$X,
            y = design$y,
            group = design$group,
            fold = design$fold
        ),
        run
    )
    coefs <- .groupedCoef(as.numeric(coef(cvFit))[-1L], design)
    list(
        coef = coefs$coef,
        slim = .grpregSlimFit(cvFit, coefs$coef, coefs$intercept),
        object = cvFit
    )
}

# cv.grpreg on `data` (design, phenotype, groups and folds) at the run's
# penalty and mixing, keeping grpreg's standardized design only for
# `retainX`; the caller's options win over the mixing, which only MCP and
# SCAD leave them (.grpregAssertAlpha()), and over SCAD's concavity, which
# defaults to ncvreg's so a method's fits agree across its engines.
# @noRd
.grpregCv <- function(data, run) {
    callArgs <- list_modify(
        c(
            data,
            list(
                penalty = run$penalty,
                alpha = run$alpha,
                returnX = run$retainX
            ),
            .grpregDefaultGamma(run$penalty)
        ),
        !!!run$user
    )
    exec(grpreg::cv.grpreg, !!!callArgs)
}

# grpreg's concavity where it differs from ncvreg's: SCAD's 3.7, where
# grpreg's own is 4. MCP's is 3 on both.
# @noRd
.grpregDefaultGamma <- function(penalty) {
    if (identical(penalty, "grSCAD")) list(gamma = .ncvregGamma[["SCAD"]])
}

# The fit cv.grpreg selected: its coefficients over the variants fitted (a
# column per condition for several), the intercept, its lambda and the
# cross-validated error there.
# @noRd
.grpregSlimFit <- function(cvFit, coef, intercept) {
    list(
        coef = coef,
        intercept = intercept,
        lambda = cvFit$lambda.min,
        cve = cvFit$cve[[cvFit$min]],
        cvse = cvFit$cvse[[cvFit$min]]
    )
}

# cv.grpreg's whole object without its fit's copies of the data: the
# phenotype and the linear predictors at every lambda.
# @noRd
.grpregWithoutData <- function(cvFit) {
    `[[<-`(
        cvFit,
        "fit",
        value = .withoutFields(cvFit$fit, c("y", "linear.predictors"))
    )
}

#' Compute Weights Using abess
#'
#' The engine behind the \code{"abess"} backend of
#' \code{\link{l0learnWeights}}: fits \code{abess::abess}, a best subset of
#' the variants, picking the support size by cross-validation, and returns
#' the coefficients at that size. For a matrix \code{y}, a column per
#' condition, the subset is one of variants each in the model for every
#' condition or none: abess's \code{"mgaussian"} family when every condition
#' was measured on every sample, otherwise a best subset of groups on the
#' sparse group design.
#'
#' @param X A numeric matrix of predictors.
#' @param y A numeric response vector, or a matrix with a column per
#'   condition, missing where a condition was not measured.
#' @param nfolds Number of cross-validation folds. Default is 5.
#' @param methodArgs Options forwarded to \code{abess::abess}, built with
#'   \code{\link{AbessOptions}}.
#' @inheritParams l0learnWeights
#' @return A numeric matrix of variant weights, a row per column of \code{X}
#'   and a column per condition, with the fit attached as attribute
#'   \code{"fit"} unless \code{fitRetention = "none"}.
#' @importFrom checkmate assertCount
#' @keywords internal
abessWeights <- function(
    X,
    y,
    nfolds = 5,
    fitRetention = c("none", "slim", "full"),
    methodArgs = AbessOptions()
) {
    .assertMethodOptions(methodArgs, "AbessOptions", "methodArgs")
    assertCount(nfolds, positive = TRUE)
    retention <- arg_match(fitRetention)
    if (!requireNamespace("abess", quietly = TRUE)) {
        abort("Package 'abess' is required for this function.")
    }
    keep <- .dropZeroVariance(X, "abessWeights")
    fitted <- .abessFitFor(
        X[, keep, drop = FALSE],
        y,
        nfolds,
        as.list(methodArgs)
    )
    .engineWeights(fitted, keep, X, y, retention, .abessWithoutData)
}

# abess's fit of one phenotype, of conditions measured on the same samples
# (its mgaussian family), or of conditions on samples of their own (a best
# subset of groups on the sparse group design, .groupedDesign()).
# @noRd
.abessFitFor <- function(X, y, nfolds, user) {
    if (!.isConditionMatrix(y)) {
        return(.abessGaussianFit(X, y, nfolds, user))
    }
    if (!anyNA(y)) {
        return(.abessMgaussianFit(X, y, nfolds, user))
    }
    .abessGroupedFit(.groupedDesign(X, y, nfolds, sparse = TRUE), user)
}

# abess's best subset of one phenotype's variants.
# @noRd
.abessGaussianFit <- function(X, y, nfolds, user) {
    fit <- .abessCall(
        list(
            x = X,
            y = y,
            family = "gaussian",
            foldid = .coupledFoldIds(seq_len(nrow(X)), nfolds)
        ),
        user
    )
    coefs <- .abessBestCoef(fit)
    list(
        coef = coefs[-1L],
        slim = .abessSlimFit(fit, coefs[-1L], coefs[[1L]]),
        object = fit
    )
}

# abess's best subset shared by conditions measured on the same samples (its
# mgaussian family, which stops on a phenotype matrix without column names),
# each phenotype scaled to unit variance and its weights put back on its own
# scale.
# @noRd
.abessMgaussianFit <- function(X, Y, nfolds, user) {
    scales <- .groupedScales(Y)
    fit <- .abessCall(
        list(
            x = X,
            y = `colnames<-`(
                .groupedScaled(Y, scales),
                str_c("condition", seq_len(ncol(Y)))
            ),
            family = "mgaussian",
            foldid = .coupledFoldIds(seq_len(nrow(X)), nfolds)
        ),
        user
    )
    B <- unname(as.matrix(coef(fit, support.size = fit$best.size)[[1L]]))
    weights <- sweep(B[-1L, , drop = FALSE], 2L, scales$spread, "*")
    list(
        coef = weights,
        slim = .abessSlimFit(
            fit,
            weights,
            scales$center + scales$spread * B[1L, ]
        ),
        object = fit
    )
}

# abess's best subset of groups on the sparse group design, which it never
# normalizes and which is centred within each condition already.
# @noRd
.abessGroupedFit <- function(design, user) {
    fit <- .abessCall(
        list(
            x = design$X,
            y = design$y,
            family = "gaussian",
            group.index = design$group,
            foldid = design$fold,
            fit.intercept = FALSE
        ),
        user
    )
    coefs <- .groupedCoef(.abessBestCoef(fit)[-1L], design)
    list(
        coef = coefs$coef,
        slim = .abessSlimFit(fit, coefs$coef, coefs$intercept),
        object = fit
    )
}

# abess on `args`, picking the support size by cross-validation, on one
# thread, with the caller's options. abess reseeds R's random number
# generator as it returns, so the session's state is restored around it.
# @noRd
.abessCall <- function(args, user) {
    callArgs <- list_modify(
        c(args, list(tune.type = "cv", num.threads = 1L)),
        !!!user
    )
    withr::with_preserve_seed(exec(abess::abess, !!!callArgs))
}

# The intercept and coefficients of abess's fit at the support size its
# cross-validation picked.
# @noRd
.abessBestCoef <- function(fit) {
    as.numeric(coef(fit, support.size = fit$best.size))
}

# The fit abess selected: its coefficients over the variants fitted (a
# column per condition for several), the intercept, its support size and the
# cross-validated error there.
# @noRd
.abessSlimFit <- function(fit, coef, intercept) {
    list(
        coef = coef,
        intercept = intercept,
        supportSize = fit$best.size,
        cvMean = fit$tune.value[[match(fit$best.size, fit$support.size)]]
    )
}

# abess's whole object without the data its stored call holds.
# @noRd
.abessWithoutData <- function(fit) {
    `[[<-`(fit, "call", value = .callWithoutData(fit$call))
}

# Weights over every column of X from an engine's fit over the kept ones
# (`keep`), with the fit kept as `retention` asks. `fitted` holds the
# coefficients (`coef`, a column per condition for several), the selected
# fit (`slim`) and the engine's object (`object`), which `withoutData`
# strips of its copies of the data.
# @noRd
.engineWeights <- function(fitted, keep, X, y, retention, withoutData) {
    weights <- if (is.matrix(fitted$coef)) {
        .weightsOnColumns(fitted$coef, keep, colnames(X), colnames(y))
    } else {
        replace(matrix(0, nrow = ncol(X), ncol = 1L), keep, fitted$coef)
    }
    if (identical(retention, "none")) {
        return(weights)
    }
    kept <- if (identical(retention, "slim")) {
        fitted$slim
    } else {
        withoutData(fitted$object)
    }
    `attr<-`(weights, "fit", .withFitVariantIds(kept, colnames(X)[keep]))
}

# The rows of `Y` observing every condition, for an engine that takes no
# missing phenotype; an error, naming the engine (`label`) and one that does
# (`instead`), when fewer than two remain.
# @noRd
.completeRows <- function(Y, label, instead) {
    rows <- which(stats::complete.cases(Y))
    if (length(rows) < 2L) {
        abort(c(
            glue(
                "{label}: fewer than two samples observe every condition, ",
                "and this engine takes no missing phenotype."
            ),
            i = glue("Conditions on samples of their own take {instead}.")
        ))
    }
    rows
}

# Whether `y` holds several conditions' phenotypes: a matrix with a column
# per condition.
# @noRd
.isConditionMatrix <- function(y) {
    is.matrix(y) && ncol(y) > 1L
}

# The group design of a fit across the conditions of `Y`: condition c's rows
# hold its own samples, X centred within them, against its own copy of the
# variants, ordered so each variant's copies sit together as its group (as
# abess's group.index needs); its phenotype is scaled to unit variance; and
# its rows are weighted by sqrt(N / (C n_c)), for N rows in all and C
# conditions, so its loss is its own mean squared error. A copy that does not
# vary in its condition's samples is left out. Each sample is dealt into one
# fold (.groupedFolds()), which its rows take in every condition.
# @noRd
.groupedDesign <- function(X, Y, nfolds, sparse) {
    rows <- map(seq_len(ncol(Y)), .coupledStudyRows, Y = Y)
    scales <- .groupedScales(Y)
    blocks <- pmap(
        list(
            seq_len(ncol(Y)),
            rows,
            sqrt(sum(lengths(rows)) / (ncol(Y) * lengths(rows)))
        ),
        .groupedBlock,
        X = X,
        Y = Y,
        scales = scales
    )
    cells <- cbind(
        variant = list_c(map(blocks, "variants")),
        condition = rep(seq_along(blocks), lengths(map(blocks, "variants")))
    )
    byVariant <- order(cells[, "variant"], cells[, "condition"])
    design <- bdiag(map(blocks, "X"))[, byVariant, drop = FALSE]
    list(
        X = if (sparse) design else as.matrix(design),
        y = list_c(map(blocks, "y")),
        group = cells[byVariant, "variant"],
        cells = cells[byVariant, , drop = FALSE],
        fold = .groupedFolds(Y, nfolds)[list_c(rows)],
        means = do.call(cbind, map(blocks, "means")),
        scales = scales
    )
}

# Condition k's block of the group design: its samples' genotypes centred
# within them and weighted, over the variants that vary there, and its
# phenotype scaled and weighted.
# @noRd
#' @importFrom matrixStats colSds
.groupedBlock <- function(k, rows, weight, X, Y, scales) {
    Xk <- X[rows, , drop = FALSE]
    means <- colMeans(Xk)
    sds <- colSds(Xk)
    varies <- !is.na(sds) & sds > 0
    list(
        X = weight * sweep(Xk[, varies, drop = FALSE], 2L, means[varies]),
        y = weight * (Y[rows, k] - scales$center[[k]]) / scales$spread[[k]],
        variants = which(varies),
        means = means
    )
}

# Each condition's phenotype centre and spread, on the samples observing it.
# @noRd
#' @importFrom matrixStats colSds
.groupedScales <- function(Y) {
    spread <- colSds(Y, na.rm = TRUE)
    flat <- which(is.na(spread) | spread == 0)
    if (length(flat) > 0L) {
        msg <- glue(
            "A fit across conditions needs each one's phenotype to vary, ",
            "but condition(s) {str_flatten(flat, ', ')} do not."
        )
        abort(msg)
    }
    list(center = colMeans(Y, na.rm = TRUE), spread = spread)
}

# `Y` with each condition's phenotype centred and scaled (.groupedScales()).
# @noRd
.groupedScaled <- function(Y, scales) {
    sweep(sweep(Y, 2L, scales$center), 2L, scales$spread, "/")
}

# Each sample's fold, dealt at random among the samples observing the same
# conditions, so every fold holds out a share of each condition's samples.
# @noRd
.groupedFolds <- function(Y, nfolds) {
    strata <- split(
        seq_len(nrow(Y)),
        vec_group_id(as.data.frame(!is.na(Y)))
    )
    replace(
        integer(nrow(Y)),
        list_c(strata),
        list_c(map(strata, .coupledFoldIds, nfolds = nfolds))
    )
}

# A group fit's coefficients (`beta`, over the design's columns) as weights on
# each condition's own scale, a row per variant and a column per condition,
# and each condition's intercept there.
# @noRd
.groupedCoef <- function(beta, design) {
    B <- replace(
        matrix(0, nrow(design$means), ncol(design$means)),
        design$cells,
        beta
    )
    list(
        coef = sweep(B, 2L, design$scales$spread, "*"),
        intercept = design$scales$center -
            design$scales$spread * colSums(design$means * B)
    )
}

# ---- Fits across a joint group's conditions --------------------------------
#
# The implementations the pipeline names for a method fitted across a group's
# conditions at once (.jointTwasGrouped()). Each takes `Y`, which sends the
# pipeline's fits, and each cross-validation fold's, down the whole-Y path
# (.twasFitsWholeY()), and hands it to the method's own wrapper on the run's
# backend.

# @noRd
.lassoGroupedWeights <- function(
    X,
    Y,
    backend = NULL,
    fitRetention = c("none", "slim", "full"),
    retainX = FALSE,
    methodArgs = NULL
) {
    lassoWeights(
        X,
        Y,
        backend = backend,
        fitRetention = arg_match(fitRetention),
        retainX = retainX,
        methodArgs = methodArgs
    )
}

# @noRd
.enetGroupedWeights <- function(
    X,
    Y,
    alpha = 0.5,
    backend = NULL,
    fitRetention = c("none", "slim", "full"),
    retainX = FALSE,
    methodArgs = NULL
) {
    enetWeights(
        X,
        Y,
        alpha = alpha,
        backend = backend,
        fitRetention = arg_match(fitRetention),
        retainX = retainX,
        methodArgs = methodArgs
    )
}

# @noRd
.mcpGroupedWeights <- function(
    X,
    Y,
    nfolds = 5,
    backend = NULL,
    fitRetention = c("none", "slim", "full"),
    retainX = FALSE,
    methodArgs = NULL
) {
    mcpWeights(
        X,
        Y,
        nfolds = nfolds,
        backend = backend,
        fitRetention = arg_match(fitRetention),
        retainX = retainX,
        methodArgs = methodArgs
    )
}

# @noRd
.scadGroupedWeights <- function(
    X,
    Y,
    nfolds = 5,
    backend = NULL,
    fitRetention = c("none", "slim", "full"),
    retainX = FALSE,
    methodArgs = NULL
) {
    scadWeights(
        X,
        Y,
        nfolds = nfolds,
        backend = backend,
        fitRetention = arg_match(fitRetention),
        retainX = retainX,
        methodArgs = methodArgs
    )
}

# @noRd
.l0learnGroupedWeights <- function(
    X,
    Y,
    penalty = "L0",
    nFolds = 5,
    backend = NULL,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    l0learnWeights(
        X,
        Y,
        penalty = penalty,
        nFolds = nFolds,
        backend = backend,
        fitRetention = arg_match(fitRetention),
        methodArgs = methodArgs
    )
}

# ---- Fits across coupled conditions (difference penalty) -------------------
#
# A group of conditions -- studies, contexts or traits -- holds X over the
# union of their variants and Y one column per condition, missing outside the
# condition's own samples; the conditions of one study share its samples.
# lasso, elastic net, ridge, MCP, SCAD and L0 fit such a group with a penalty
# on the difference between coupled conditions' per-SD effects, the
# individual-level analogue of PROSPER:
#
#   sum_k [RSS_k / (2 n_k) + lambda P(b_k)]
#       + sum_j sum_{k < k'} (c_kk' / 2) (s_jk b_jk - s_jk' b_jk')^2
#
# Each condition's loss is its own mean squared error, on its phenotype
# scaled to unit variance, and s_jk is variant j's SD in condition k's
# samples. A pair differing along one coupled axis alone takes that axis's
# strength c, any other pair none (.coupledSpec()); several coupled axes give
# a grid of every combination of their strengths. The engines solve it as one
# augmented least-squares problem (.coupledDesign()), each condition on its
# own copy of the variants. An engine with a group penalty adds it, each
# variant's copies one group: grpreg and abess over every condition's copy,
# glmnet's mgaussian family over the conditions that share samples and are
# not coupled to one another (.coupledUnitsHold()). Each condition takes the
# grid point its own held-out samples score best, on folds that hold a sample
# out of every condition at once, the pseudo-rows always training
# (.coupledFit()). At c = 0 the fit is the conditions' separate fits (or
# their group fit); as c grows it becomes one pooled fit.
#
# ncvreg and grpreg always standardize their columns, and the pseudo-rows
# enter each column's scale. The weights on their columns put the lasso part
# of the penalty on the per-SD effects (.coupledNcvregPath()); MCP's
# concavity is then exact only where every column has the same coupling, and
# SCAD's, and grpreg's group norm, only at c = 0 (.coupledWarnApproximate()).

# lasso across a group's coupled conditions, on glmnet or grpreg
# (.penalizedCoupledWeights()). `coupling` is the group's coupling
# (.coupledSpec()), or a grid of strengths coupling every pair.
# @noRd
.lassoCoupledWeights <- function(
    X,
    Y,
    coupling,
    backend = NULL,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    .penalizedCoupledWeights(
        list(X = X, Y = Y, coupling = coupling),
        list(token = "lasso", alpha = 1, backend = backend),
        arg_match(fitRetention),
        methodArgs
    )
}

# The elastic net across a group's coupled conditions, at the mixing `alpha`
# enetWeights() takes.
# @noRd
.enetCoupledWeights <- function(
    X,
    Y,
    coupling,
    alpha = 0.5,
    backend = NULL,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    assertNumber(alpha, lower = 0, upper = 1)
    .penalizedCoupledWeights(
        list(X = X, Y = Y, coupling = coupling),
        list(token = "enet", alpha = alpha, backend = backend),
        arg_match(fitRetention),
        methodArgs
    )
}

# Ridge across a group's coupled conditions, on glmnet.
# @noRd
.ridgeCoupledWeights <- function(
    X,
    Y,
    coupling,
    fitRetention = c("none", "slim", "full"),
    methodArgs = GlmnetOptions()
) {
    .penalizedCoupledWeights(
        list(X = X, Y = Y, coupling = coupling),
        list(token = "ridge", alpha = 0),
        arg_match(fitRetention),
        methodArgs
    )
}

# Ridge across a group's conditions with no coupling: each condition's own
# ridge, all at the one lambda whose cross-validated error, averaged over the
# conditions, is least. Ridge has no group penalty, its penalty being separate
# for each condition; sharing lambda is what joins them.
# @noRd
.ridgeGroupedWeights <- function(
    X,
    Y,
    fitRetention = c("none", "slim", "full"),
    methodArgs = GlmnetOptions()
) {
    .penalizedCoupledWeights(
        list(X = X, Y = Y, coupling = .coupledNone(ncol(Y))),
        list(token = "ridge", alpha = 0),
        arg_match(fitRetention),
        methodArgs,
        choose = "pooled"
    )
}

# MCP across a group's coupled conditions, on ncvreg (the coupling alone) or
# grpreg (with the group penalty).
# @noRd
.mcpCoupledWeights <- function(
    X,
    Y,
    coupling,
    nfolds = 5,
    backend = NULL,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    assertCount(nfolds, positive = TRUE)
    .penalizedCoupledWeights(
        list(X = X, Y = Y, coupling = coupling),
        list(
            token = "mcp",
            penalty = "MCP",
            nfolds = nfolds,
            backend = backend
        ),
        arg_match(fitRetention),
        methodArgs
    )
}

# SCAD across a group's coupled conditions, as .mcpCoupledWeights().
# @noRd
.scadCoupledWeights <- function(
    X,
    Y,
    coupling,
    nfolds = 5,
    backend = NULL,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    assertCount(nfolds, positive = TRUE)
    .penalizedCoupledWeights(
        list(X = X, Y = Y, coupling = coupling),
        list(
            token = "scad",
            penalty = "SCAD",
            nfolds = nfolds,
            backend = backend
        ),
        arg_match(fitRetention),
        methodArgs
    )
}

# L0 across a group's coupled conditions, on L0Learn (the coupling alone) or
# abess (a best subset of groups as well). The L0 penalty alone: its lambda
# reads the same at every c, where the L1 and L2 terms' would not, as they
# see L0Learn's own rescaling of the columns, which the pseudo-rows change.
# @noRd
.l0learnCoupledWeights <- function(
    X,
    Y,
    coupling,
    penalty = "L0",
    nFolds = 5,
    backend = NULL,
    fitRetention = c("none", "slim", "full"),
    methodArgs = NULL
) {
    .coupledRequireL0(penalty)
    assertCount(nFolds, positive = TRUE)
    .penalizedCoupledWeights(
        list(X = X, Y = Y, coupling = coupling),
        list(token = "l0learn", nfolds = nFolds, backend = backend),
        arg_match(fitRetention),
        methodArgs
    )
}

# @noRd
.coupledRequireL0 <- function(penalty) {
    if (identical(penalty, "L0")) {
        return(invisible(NULL))
    }
    abort(c(
        glue(
            "l0learn across coupled conditions fits penalty = \"L0\" only, ",
            "not \"{penalty}\": the L1 and L2 terms see the engine's own ",
            "rescaling of the columns, which the coupling changes."
        ),
        i = str_c(
            "Use penalty = \"L0\", or pass coupling = NULL to fit each ",
            "condition on its own."
        )
    ))
}

# A penalized method's fit across the coupled conditions of `data` (X, Y and
# `coupling`), on the engine its backend names (`method`: the token, its
# backend and settings), each condition at its own grid point or, for
# `choose = "pooled"`, all at one.
# @noRd
.penalizedCoupledWeights <- function(
    data,
    method,
    retention,
    methodArgs,
    choose = "own"
) {
    backend <- .coupledBackend(method)
    options <- .backendOptions(
        methodArgs,
        .coupledCtorName(method$token, backend),
        backend
    )
    .coupledRequirePackage(backend)
    conds <- .coupledConditions(data$X, data$Y, data$coupling)
    engine <- .coupledEngine(backend, method, as.list(options), conds)
    .coupledWarnApproximate(backend, method, conds)
    .coupledFit(conds, engine, retention, choose)
}

# The engine a coupled fit runs on: ridge's is glmnet's, the others' the
# backend chosen for the method.
# @noRd
.coupledBackend <- function(method) {
    if (identical(method$token, "ridge")) {
        return("glmnet")
    }
    .methodBackend(method$token, method$backend, "individual")
}

# @noRd
.coupledCtorName <- function(token, backend) {
    if (identical(token, "ridge")) {
        return("GlmnetOptions")
    }
    .methodBackends$individual[[token]][[backend]]
}

# @noRd
.coupledRequirePackage <- function(backend) {
    package <- .coupledPackages[[backend]]
    if (!requireNamespace(package, quietly = TRUE)) {
        abort(glue("Package '{package}' is required for this function."))
    }
    invisible(NULL)
}

# The package each coupled engine is.
# @noRd
.coupledPackages <- c(
    glmnet = "glmnet",
    l0learn = "L0Learn",
    ncvreg = "ncvreg",
    grpreg = "grpreg",
    abess = "abess"
)

# The engine's path function, grid and settings for a coupled fit
# (`method`'s, with the caller's engine options `user`).
# @noRd
.coupledEngine <- function(backend, method, user, conds) {
    switch(
        backend,
        glmnet = .coupledGlmnetEngine(method, user, conds),
        l0learn = .coupledL0learnEngine(method, user, conds),
        ncvreg = .coupledNcvregEngine(method, user, conds),
        grpreg = .coupledGrpregEngine(method, user, conds),
        abess = .coupledAbessEngine(method, user, conds)
    )
}

# The coupling of a group of K conditions: which axis each pair is coupled
# along (`pairs`, K x K, NA for a pair not coupled), each coupled axis's grid
# of strengths (`grids`), and the conditions' units and columns, which
# glmnet's mgaussian family fits together (`unit`, `column`). A numeric
# `coupling` is a grid coupling every pair, each condition its own unit.
# @noRd
#' @importFrom checkmate assertNumeric assertList assertMatrix
.coupledSpec <- function(coupling, K) {
    if (is.numeric(coupling)) {
        .coupledAssertGrid(coupling)
        return(list(
            pairs = .coupledAllPairs(K),
            grids = list(strength = coupling),
            unit = seq_len(K),
            column = rep(1L, K)
        ))
    }
    assertList(coupling)
    assertMatrix(coupling$pairs, nrows = K, ncols = K)
    walk(coupling$grids, .coupledAssertGrid)
    list_modify(.coupledNone(K), !!!coupling)
}

# @noRd
.coupledAssertGrid <- function(strengths) {
    assertNumeric(
        strengths,
        lower = 0,
        finite = TRUE,
        any.missing = FALSE,
        min.len = 1L
    )
}

# Every pair of K conditions coupled along one axis, "strength".
# @noRd
.coupledAllPairs <- function(K) {
    `diag<-`(matrix("strength", K, K), NA_character_)
}

# A group of K conditions with no coupling.
# @noRd
.coupledNone <- function(K) {
    list(
        pairs = matrix(NA_character_, K, K),
        grids = list(),
        unit = seq_len(K),
        column = rep(1L, K)
    )
}

# The conditions of a coupled fit: X, Y, each condition's rows (those holding
# its phenotype), its coupling and the grid's cells, one row per combination
# of the coupled axes' strengths.
# @noRd
#' @importFrom tidyr expand_grid
.coupledConditions <- function(X, Y, coupling) {
    spec <- .coupledSpec(coupling, ncol(Y))
    list(
        X = X,
        Y = Y,
        rows = map(seq_len(ncol(Y)), .coupledStudyRows, Y = Y),
        spec = spec,
        cells = expand_grid(!!!spec$grids)
    )
}

# @noRd
.coupledStudyRows <- function(k, Y) {
    which(!is.na(Y[, k]))
}

# The strengths of every pair of conditions in cell `i` (K x K, zero for a
# pair not coupled).
# @noRd
.coupledStrength <- function(i, conds) {
    .coupledPairStrengths(c(numeric(0), unlist(conds$cells[i, ])), conds)
}

# The weakest coupling of the grid, each axis at its least strength: where a
# grid that serves every cell is set.
# @noRd
.coupledWeakest <- function(conds) {
    .coupledPairStrengths(map_dbl(conds$spec$grids, min), conds)
}

# The pairs' strengths (K x K) at the axes' strengths `values`, named by
# axis: a pair's is its axis's, zero for a pair not coupled.
# @noRd
.coupledPairStrengths <- function(values, conds) {
    axes <- conds$spec$pairs
    W <- matrix(unname(values[axes]), nrow(axes), ncol(axes))
    replace(W, is.na(W), 0)
}

# ---- engines of a coupled fit ----

# glmnet: the coupling alone, or for lasso and the elastic net with
# mgaussian's group penalty over the conditions of each unit, where the
# conditions fall into units (.coupledUnitsHold()).
# @noRd
.coupledGlmnetEngine <- function(method, user, conds) {
    opts <- .coupledGlmnetOptions(user)
    grouped <- !identical(method$token, "ridge") && .coupledUnitsHold(conds)
    engine <- list(
        path = .coupledGlmnetPath,
        alpha = method$alpha,
        args = opts$engine,
        layout = if (grouped) "units" else "conditions",
        nfolds = opts$nfolds,
        gridName = "lambda"
    )
    list_assign(
        engine,
        grid = opts$lambda %||% .coupledGlmnetLambda(conds, engine, opts)
    )
}

# The GlmnetOptions a coupled fit honours: its folds' count (`nfolds`), its
# lambda grid (`lambda`, `nlambda`, `lambda.min.ratio`) and glmnet's
# tolerances. The rest set up what this fit does not do -- cv.glmnet's own
# folds, or a per-variant layout the augmented design lacks.
# @noRd
.coupledGlmnetOptions <- function(user) {
    honoured <- c(
        "nfolds",
        "lambda",
        "nlambda",
        "lambda.min.ratio",
        "thresh",
        "maxit",
        "type.measure"
    )
    .coupledRefuseOptions(setdiff(names(user), honoured), "GlmnetOptions")
    .coupledRefuseMeasure(user$type.measure)
    list(
        nfolds = user$nfolds %||% 5L,
        lambda = user$lambda,
        nlambda = user$nlambda %||% 100L,
        ratio = user$lambda.min.ratio %||% 0.01,
        engine = user[intersect(names(user), c("thresh", "maxit"))]
    )
}

# A coupled fit picks its grid point by mean squared error, which is what
# cv.glmnet's "mse" and, for a gaussian response, "deviance" measure.
# @noRd
.coupledRefuseMeasure <- function(measure) {
    if (is.null(measure) || is_in(measure, c("default", "mse", "deviance"))) {
        return(invisible(NULL))
    }
    abort(glue(
        "GlmnetOptions(type.measure = \"{measure}\"): a fit across coupled ",
        "conditions chooses its grid point by mean squared error."
    ))
}

# Engine options that have no meaning for a coupled fit.
# @noRd
.coupledRefuseOptions <- function(unsupported, label) {
    if (length(unsupported) == 0L) {
        return(invisible(NULL))
    }
    abort(c(
        glue(
            "{label}({str_flatten(unsupported, ', ')}) does not apply to a ",
            "fit across coupled conditions, whose folds, grid and columns are ",
            "its own."
        ),
        i = str_c(
            "Leave it unset, or pass coupling = NULL to fit each condition ",
            "on its own."
        )
    ))
}

# glmnet's lambda grid: from the largest any unit's own fit starts at, where
# every effect is zero, down by `ratio`. The pseudo-rows add nothing to the
# gradient at zero, so one grid serves every c.
# @noRd
.coupledGlmnetLambda <- function(conds, engine, opts) {
    copies <- .coupledCopies(conds, conds$rows, engine$layout)
    top <- max(map_dbl(copies, .coupledLambdaMax, alpha = engine$alpha))
    exp(seq(log(top), log(top * opts$ratio), length.out = opts$nlambda))
}

# A copy's own largest lambda: glmnet's, the largest norm of a variant's
# correlations with the copy's phenotypes, which takes a ridge's from
# alpha = 0.001.
# @noRd
.coupledLambdaMax <- function(copy, alpha) {
    norms <- sqrt(rowSums(crossprod(copy$X, copy$Y)^2))
    max(norms) / (copy$n * max(alpha, 1e-3))
}

# One glmnet path on `design`, its coefficients over the grid (columns),
# padded with the last fit where glmnet stopped early, its deviance no longer
# moving: gaussian on one phenotype column, mgaussian on a unit's several,
# their coefficients stacked a column after another. The design's loss is
# the conditions' mean, so lambda is divided by K.
# @noRd
.coupledGlmnetPath <- function(design, engine) {
    several <- is.matrix(design$y)
    fit <- exec(
        glmnet::glmnet,
        x = design$X,
        y = design$y,
        family = if (several) "mgaussian" else "gaussian",
        alpha = engine$alpha,
        lambda = engine$grid / design$K,
        standardize = FALSE,
        intercept = FALSE,
        !!!c(if (several) list(standardize.response = FALSE), engine$args)
    )
    beta <- if (several) map(fit$beta, as.matrix) else list(fit$beta)
    .coupledPadPath(
        as.matrix(exec(rbind, !!!beta)),
        length(engine$grid)
    )
}

# L0Learn: the coupling alone, with the L0 penalty.
# @noRd
.coupledL0learnEngine <- function(method, user, conds) {
    engine <- list(
        path = .coupledL0learnPath,
        args = .coupledL0learnOptions(user),
        layout = "conditions",
        nfolds = method$nfolds,
        gridName = "lambda"
    )
    list_assign(engine, grid = .coupledL0learnLambda(conds, engine))
}

# The L0learnOptions a coupled fit honours: those of L0Learn's search and
# grid. `seed` sets L0Learn.cvfit()'s folds, which this fit does not use; its
# own folds follow the run's seed.
# @noRd
.coupledL0learnOptions <- function(user) {
    honoured <- c(
        "maxSuppSize",
        "nLambda",
        "scaleDownFactor",
        "algorithm",
        "maxSwaps",
        "maxIters",
        "rtol",
        "atol",
        "activeSet",
        "activeSetNum",
        "screenSize"
    )
    .coupledRefuseOptions(
        setdiff(names(user), c(honoured, "seed")),
        "L0learnOptions"
    )
    user[intersect(names(user), honoured)]
}

# L0Learn's own lambda grid for the full data at the weakest coupling. Its
# lambda is relative to the phenotype's norm, which every design fixes, so the
# grid reads the same at every c and in every fold.
# @noRd
.coupledL0learnLambda <- function(conds, engine) {
    design <- .coupledDesign(
        conds,
        conds$rows,
        .coupledWeakest(conds),
        engine$layout
    )
    fit <- exec(
        L0Learn::L0Learn.fit,
        x = design$X,
        y = design$y,
        penalty = "L0",
        intercept = FALSE,
        !!!engine$args
    )
    fit$lambda[[1L]]
}

# One L0Learn path on `design` over the grid, padded as glmnet's is.
# @noRd
.coupledL0learnPath <- function(design, engine) {
    gridArgs <- engine$args[
        setdiff(names(engine$args), c("nLambda", "scaleDownFactor"))
    ]
    fit <- exec(
        L0Learn::L0Learn.fit,
        x = design$X,
        y = design$y,
        penalty = "L0",
        intercept = FALSE,
        lambdaGrid = list(engine$grid),
        !!!gridArgs
    )
    .coupledPadPath(as.matrix(fit$beta[[1L]]), length(engine$grid))
}

# ncvreg: the coupling alone, MCP or SCAD on a dense design.
# @noRd
.coupledNcvregEngine <- function(method, user, conds) {
    opts <- .coupledPerSdOptions(user, "NcvregOptions", "convex")
    engine <- list(
        path = .coupledNcvregPath,
        penalty = method$penalty,
        gamma = opts$gamma %||% .ncvregGamma[[method$penalty]],
        args = opts$engine,
        layout = "conditions",
        nfolds = method$nfolds,
        gridName = "lambda"
    )
    list_assign(
        engine,
        grid = opts$lambda %||% .coupledPerSdLambda(conds, opts)
    )
}

# ncvreg's own concavity defaults.
# @noRd
.ncvregGamma <- c(MCP = 3, SCAD = 3.7)

# The concavity of a coupled fit: the caller's, else ncvreg's for MCP and
# SCAD, so a method's fits agree across its engines; none for the lasso.
# @noRd
.coupledGamma <- function(gamma, penalty) {
    if (is.null(penalty)) {
        return(gamma)
    }
    gamma %||% .ncvregGamma[[penalty]]
}

# The NcvregOptions or GrpregOptions (`label`) a coupled fit honours: the
# concavity, the grid and the tolerances, plus `extra`. The mixing with ridge
# (`alpha`) is refused: its ridge part sees the engine's rescaling of the
# columns, which the coupling changes.
# @noRd
.coupledPerSdOptions <- function(user, label, extra) {
    honoured <- c(
        "gamma",
        "nlambda",
        "lambda.min",
        "lambda",
        "eps",
        "max.iter",
        "dfmax",
        "warn",
        extra
    )
    .coupledRefuseOptions(setdiff(names(user), honoured), label)
    list(
        gamma = user$gamma,
        lambda = user$lambda,
        nlambda = user$nlambda %||% 100L,
        ratio = user$lambda.min,
        engine = user[intersect(names(user), setdiff(honoured, .gridOptions))]
    )
}

# The options setting the grid, which a coupled fit sets itself.
# @noRd
.gridOptions <- c("gamma", "nlambda", "lambda.min", "lambda")

# A lambda grid on the per-SD scale: from the largest any condition's own
# fit starts at, where every effect is zero, down by `lambda.min` (ncvreg's
# default: 0.001 with more samples than variants, else 0.05).
# @noRd
.coupledPerSdLambda <- function(conds, opts) {
    copies <- .coupledCopies(conds, conds$rows, "conditions")
    top <- max(map_dbl(copies, .coupledPerSdMax))
    few <- any(map_int(copies, "n") <= ncol(conds$X))
    ratio <- opts$ratio %||% (if (few) 0.05 else 0.001)
    exp(seq(log(top), log(top * ratio), length.out = opts$nlambda))
}

# @noRd
.coupledPerSdMax <- function(copy) {
    max(abs(crossprod(copy$X, copy$Y)) / copy$s) / copy$n
}

# One ncvreg path on `design`. ncvreg penalizes its standardized columns, the
# pseudo-rows among them, so column j's coefficient is r_j b_j, r_j^2 =
# s_j^2 (1 + C_j) / K, C_j the strengths summed over its coupled partners.
# The weight sqrt(K / (1 + C_j)) on each column puts the lasso part of the
# penalty on the per-SD effect s_j b_j at lambda / K, the design's loss being
# the conditions' mean; MCP's concavity there is gamma (1 + C_j), exact when
# C_j is the same for every column, and SCAD's is left at gamma.
# @noRd
.coupledNcvregPath <- function(design, engine) {
    scale <- 1 + design$partners
    gamma <- if (identical(engine$penalty, "MCP")) {
        engine$gamma * stats::median(scale)
    } else {
        engine$gamma
    }
    fit <- exec(
        ncvreg::ncvreg,
        X = as.matrix(design$X),
        y = design$y,
        penalty = engine$penalty,
        gamma = gamma,
        lambda = engine$grid / design$K,
        penalty.factor = sqrt(design$K / scale),
        !!!engine$args
    )
    .coupledPadPath(
        as.matrix(fit$beta)[-1L, , drop = FALSE],
        length(engine$grid)
    )
}

# grpreg: the coupling with the group penalty, each variant's copies one
# group, on a dense design.
# @noRd
.coupledGrpregEngine <- function(method, user, conds) {
    if (!is.null(method$alpha)) {
        .grpregAssertAlpha(user, method$alpha)
    }
    opts <- .coupledPerSdOptions(
        user,
        "GrpregOptions",
        c("tau", "gmax", "log.lambda")
    )
    engine <- list(
        path = .coupledGrpregPath,
        penalty = str_c("gr", method$penalty %||% "Lasso"),
        alpha = method$alpha %||% 1,
        args = c(
            opts$engine,
            compact(list(
                gamma = .coupledGamma(opts$gamma, method$penalty)
            ))
        ),
        layout = "groups",
        nfolds = method$nfolds %||% 5L,
        gridName = "lambda"
    )
    list_assign(
        engine,
        grid = opts$lambda %||% .coupledGrpregLambda(conds, engine, opts)
    )
}

# grpreg's own lambda grid for the full data at the weakest coupling, which
# then serves every c and fold.
# @noRd
.coupledGrpregLambda <- function(conds, engine, opts) {
    design <- .coupledDesign(
        conds,
        conds$rows,
        .coupledWeakest(conds),
        engine$layout
    )
    fit <- exec(
        grpreg::grpreg,
        X = as.matrix(design$X),
        y = design$y,
        group = design$group,
        penalty = engine$penalty,
        alpha = engine$alpha,
        nlambda = opts$nlambda,
        !!!c(compact(list(lambda.min = opts$ratio)), engine$args)
    )
    fit$lambda
}

# One grpreg path on `design` over the grid, padded as glmnet's is.
# @noRd
.coupledGrpregPath <- function(design, engine) {
    fit <- exec(
        grpreg::grpreg,
        X = as.matrix(design$X),
        y = design$y,
        group = design$group,
        penalty = engine$penalty,
        alpha = engine$alpha,
        lambda = engine$grid,
        !!!engine$args
    )
    .coupledPadPath(
        as.matrix(fit$beta)[-1L, , drop = FALSE],
        length(engine$grid)
    )
}

# abess: the coupling with a best subset of groups, each variant's copies one
# group, on a sparse design; its grid is the number of groups kept.
# @noRd
.coupledAbessEngine <- function(method, user, conds) {
    honoured <- c(
        "support.size",
        "c.max",
        "max.splicing.iter",
        "warm.start",
        "splicing.type",
        "important.search"
    )
    .coupledRefuseOptions(setdiff(names(user), honoured), "AbessOptions")
    engine <- list(
        path = .coupledAbessPath,
        args = user[intersect(names(user), setdiff(honoured, "support.size"))],
        layout = "groups",
        nfolds = method$nfolds,
        gridName = "supportSize"
    )
    list_assign(
        engine,
        grid = user$support.size %||% .coupledAbessSizes(conds, engine)
    )
}

# abess's own support sizes for the full data at the weakest coupling.
# @noRd
.coupledAbessSizes <- function(conds, engine) {
    design <- .coupledDesign(
        conds,
        conds$rows,
        .coupledWeakest(conds),
        engine$layout
    )
    .coupledAbessCall(design, engine, NULL)$support.size
}

# One abess path on `design` over the support sizes of the grid.
# @noRd
.coupledAbessPath <- function(design, engine) {
    fit <- .coupledAbessCall(design, engine, engine$grid)
    as.matrix(fit$beta)
}

# abess on `design`'s groups at `sizes` (its own when NULL), on one thread,
# the session's random number state restored around it (.abessCall()).
# @noRd
.coupledAbessCall <- function(design, engine, sizes) {
    callArgs <- c(
        list(
            x = design$X,
            y = design$y,
            family = "gaussian",
            group.index = design$group,
            tune.type = "gic",
            fit.intercept = FALSE,
            num.threads = 1L
        ),
        compact(list(support.size = sizes)),
        engine$args
    )
    withr::with_preserve_seed(exec(abess::abess, !!!callArgs))
}

# A warning, once a session for each, that ncvreg's or grpreg's
# standardization leaves the penalty approximate on this group
# (.coupledNcvregPath()): grpreg's group norm counts the coupling at any
# c > 0, SCAD's concavity shrinks with it, and MCP's is exact only when
# every column has the same coupling.
# @noRd
.coupledWarnApproximate <- function(backend, method, conds) {
    if (!is_in(backend, c("ncvreg", "grpreg"))) {
        return(invisible(NULL))
    }
    copies <- .coupledPlaced(
        .coupledCopies(conds, conds$rows, "conditions"),
        "conditions"
    )
    partners <- map(
        seq_len(nrow(conds$cells)),
        .coupledCellPartners,
        conds = conds,
        copies = copies
    )
    concaveOnly <- identical(backend, "ncvreg") &&
        identical(method$token, "mcp")
    exact <- all(map_lgl(
        partners,
        if (concaveOnly) .coupledEven else .coupledUncoupled
    ))
    if (!exact) {
        warn(
            glue(
                "{method$token} on {backend} across coupled conditions: ",
                "{backend} rescales its columns, the coupling's rows among ",
                "them, so the penalty is approximate on this group."
            ),
            .frequency = "once",
            .frequency_id = str_c("pecotmr_coupled_", backend, method$token)
        )
    }
    invisible(NULL)
}

# The design columns' couplings at cell `i` (.coupledPartners()).
# @noRd
.coupledCellPartners <- function(i, conds, copies) {
    .coupledPartners(
        .coupledPairs(copies, .coupledStrength(i, conds)),
        sum(map_int(copies, .coupledWidth))
    )
}

# @noRd
.coupledEven <- function(partners) {
    diff(range(partners)) <= 1e-12
}

# @noRd
.coupledUncoupled <- function(partners) {
    all(partners == 0)
}

# ---- the coupled fit ----

# A path's coefficients over `n` grid points: an engine that stopped early
# has reached its last fit, which then stands for the rest.
# @noRd
.coupledPadPath <- function(B, n) {
    if (ncol(B) >= n) {
        return(B[, seq_len(n), drop = FALSE])
    }
    cbind(B, B[, rep(ncol(B), n - ncol(B)), drop = FALSE])
}

# The coupled fit of `conds` on `engine` (its path function, grid and
# settings): every cell's full-data path, each condition's cross-validated
# error over the grid and cells, and the weights of each condition's own best
# point (`choose = "own"`) or of the one point best on average
# ("pooled"), with the fit kept as `retention` asks.
# @noRd
.coupledFit <- function(conds, engine, retention, choose) {
    folds <- .groupedFolds(conds$Y, engine$nfolds)
    full <- map(
        seq_len(nrow(conds$cells)),
        .coupledFullPath,
        conds = conds,
        engine = engine
    )
    cv <- .coupledCvError(conds, engine, folds)
    choice <- .coupledChoices(cv, choose)
    W <- .coupledWeightMatrix(conds, full, choice)
    if (identical(retention, "none")) {
        return(W)
    }
    fit <- .coupledFitRecord(conds, engine, full, cv, choice, retention)
    `attr<-`(W, "fit", .withFitVariantIds(fit, colnames(conds$X)))
}

# A sample's rows dealt into `nfolds` folds at random.
# @noRd
.coupledFoldIds <- function(rows, nfolds) {
    sample(rep_len(seq_len(nfolds), length(rows)))
}

# The full-data path at cell `i`, with the scales its coefficients are on.
# @noRd
.coupledFullPath <- function(i, conds, engine) {
    design <- .coupledDesign(
        conds,
        conds$rows,
        .coupledStrength(i, conds),
        engine$layout
    )
    list(blocks = design$blocks, beta = engine$path(design, engine))
}

# Each condition's cross-validated error over the grid and cells: every fold
# fitted at every cell on the other folds' samples -- the pseudo-rows always
# among them -- and scored on its own.
# @noRd
.coupledCvError <- function(conds, engine, folds) {
    perFold <- map(
        seq_len(engine$nfolds),
        .coupledFoldSse,
        conds = conds,
        engine = engine,
        folds = folds
    )
    map(
        seq_along(conds$rows),
        .coupledConditionCv,
        perFold = perFold,
        held = map(conds$rows, .coupledHeldCounts, folds, engine$nfolds)
    )
}

# @noRd
.coupledHeldCounts <- function(rows, folds, nfolds) {
    tabulate(folds[rows], nbins = nfolds)
}

# Fold `f`'s squared errors: per condition, its held-out rows' summed squared
# error at each grid point (rows: grid; columns: cells).
# @noRd
.coupledFoldSse <- function(f, conds, engine, folds) {
    train <- map(conds$rows, .coupledRowsOut, folds = folds, f = f)
    test <- map(conds$rows, .coupledRowsIn, folds = folds, f = f)
    byCell <- map(
        seq_len(nrow(conds$cells)),
        .coupledCellSse,
        conds = conds,
        train = train,
        test = test,
        engine = engine
    )
    map(seq_along(test), .coupledSseByCell, byCell = byCell)
}

# @noRd
.coupledRowsOut <- function(rows, folds, f) {
    rows[folds[rows] != f]
}

# @noRd
.coupledRowsIn <- function(rows, folds, f) {
    rows[folds[rows] == f]
}

# Each condition's held-out squared errors along the path fitted at cell `i`
# on the rows `train`.
# @noRd
.coupledCellSse <- function(i, conds, train, test, engine) {
    design <- .coupledDesign(
        conds,
        train,
        .coupledStrength(i, conds),
        engine$layout
    )
    beta <- engine$path(design, engine)
    map2(
        seq_along(test),
        test,
        .coupledHeldOutSse,
        design = design,
        beta = beta,
        conds = conds
    )
}

# @noRd
.coupledSseByCell <- function(k, byCell) {
    exec(cbind, !!!map(byCell, k))
}

# Condition `k`'s held-out rows' summed squared error at each grid point.
# @noRd
.coupledHeldOutSse <- function(k, rows, design, beta, conds) {
    if (length(rows) == 0L) {
        return(numeric(ncol(beta)))
    }
    pred <- .coupledPredict(
        design$blocks[[k]],
        beta,
        conds$X[rows, , drop = FALSE]
    )
    colSums((conds$Y[rows, k] - pred)^2)
}

# Phenotype predictions for the rows of `X` from a condition's coefficients
# along a path `beta`, on the centring and scale its fit was made on.
# @noRd
.coupledPredict <- function(block, beta, X) {
    Xc <- sweep(X[, block$carried, drop = FALSE], 2L, block$xMean)
    coefs <- block$sigma * beta[block$cols, , drop = FALSE]
    as.matrix(Xc %*% coefs) + block$yMean
}

# Condition `k`'s cross-validated mean squared error over the grid and cells,
# and its standard error, computed from the folds' errors as cv.glmnet
# computes them; `held` is each condition's count of rows in each fold.
# @noRd
.coupledConditionCv <- function(k, perFold, held) {
    sse <- map(perFold, k)
    counts <- held[[k]]
    cvm <- reduce(sse, `+`) / sum(counts)
    used <- which(counts > 0L)
    deviations <- map2(
        sse[used],
        counts[used],
        .coupledFoldDeviation,
        cvm = cvm
    )
    weighted <- map2(deviations, counts[used] / sum(counts), `*`)
    list(
        cvm = cvm,
        cvsd = sqrt(reduce(weighted, `+`) / max(length(used) - 1L, 1L))
    )
}

# @noRd
.coupledFoldDeviation <- function(sse, n, cvm) {
    (sse / n - cvm)^2
}

# Each condition's grid point: its own best, or the one whose error averaged
# over the conditions is least.
# @noRd
.coupledChoices <- function(cv, choose) {
    if (identical(choose, "own")) {
        return(map(cv, .coupledChoice))
    }
    pooled <- reduce(map(cv, "cvm"), `+`) / length(cv)
    rep(list(.coupledChoice(list(cvm = pooled))), length(cv))
}

# A grid point of least cross-validated error, the first of ties: its place
# in the grid and its cell.
# @noRd
.coupledChoice <- function(cv) {
    at <- arrayInd(which.min(cv$cvm), dim(cv$cvm))
    list(grid = at[[1L]], cell = at[[2L]])
}

# The weights (variants x conditions), each condition's from the full-data
# path at its grid point.
# @noRd
.coupledWeightMatrix <- function(conds, full, choice) {
    cols <- map2(
        seq_along(choice),
        choice,
        .coupledConditionWeights,
        full = full,
        nVariants = ncol(conds$X)
    )
    `dimnames<-`(
        exec(cbind, !!!cols),
        list(colnames(conds$X), colnames(conds$Y))
    )
}

# Condition `k`'s weights, on its phenotype's own scale, zero for each variant
# it does not carry.
# @noRd
.coupledConditionWeights <- function(k, pick, full, nVariants) {
    path <- full[[pick$cell]]
    block <- path$blocks[[k]]
    replace(
        numeric(nVariants),
        block$carried,
        block$sigma * path$beta[block$cols, pick$grid]
    )
}

# The fit as `retention` keeps it. "slim": each condition's grid point (its
# strengths, a column per coupled axis, and its lambda or support size) and
# the cross-validated error there. "full" adds the grid, its cells, every
# cell's full-data path (coefficients over each condition's carried variants,
# with their scales) and each condition's errors over the whole grid.
# @noRd
.coupledFitRecord <- function(conds, engine, full, cv, choice, retention) {
    named <- colnames(conds$Y)
    at <- map_int(choice, "grid")
    slim <- c(
        compact(list(
            strength = .coupledChosenStrengths(conds, choice, named)
        )),
        set_names(list(set_names(engine$grid[at], named)), engine$gridName),
        list(
            cvm = set_names(
                map2_dbl(cv, choice, .coupledAt, field = "cvm"),
                named
            ),
            cvsd = set_names(
                map2_dbl(cv, choice, .coupledAt, field = "cvsd"),
                named
            )
        )
    )
    if (identical(retention, "slim")) {
        return(slim)
    }
    c(
        slim,
        set_names(list(engine$grid), str_c(engine$gridName, "Path")),
        list(
            cells = conds$cells,
            paths = map(full, .coupledSparsePath),
            cvmPath = set_names(map(cv, "cvm"), named),
            cvsdPath = set_names(map(cv, "cvsd"), named)
        )
    )
}

# Each condition's chosen strengths, a row per condition and a column per
# coupled axis; NULL when nothing is coupled.
# @noRd
.coupledChosenStrengths <- function(conds, choice, named) {
    if (ncol(conds$cells) == 0L) {
        return(NULL)
    }
    picked <- conds$cells[map_int(choice, "cell"), , drop = FALSE]
    `dimnames<-`(as.matrix(picked), list(named, names(conds$cells)))
}

# @noRd
.coupledAt <- function(cv, pick, field) {
    cv[[field]][pick$grid, pick$cell]
}

# @noRd
#' @importFrom Matrix Matrix
.coupledSparsePath <- function(path) {
    list_assign(path, beta = Matrix(path$beta, sparse = TRUE))
}

# ---- the augmented design ----

# The augmented least-squares problem of the conditions' rows `train` (one
# index vector per condition) at the pairs' strengths `W`, laid out as
# `layout` asks: "conditions", each condition a copy of the variants with its
# own columns; "groups", the same with each variant's copies side by side as
# its group (as abess's group.index needs); "units", each unit's conditions
# sharing one copy, a phenotype column each (.coupledCopies()). Each copy's
# rows are centred, its phenotypes scaled to unit variance, and both scaled
# by sqrt(N / (K n)) for K conditions: each condition's loss is then its mean
# squared error. For each variant two coupled copies carry, the pseudo-rows
# +/- sqrt(N c / (2 K)) (s_j e_j - s'_j e'_j) add (c / (2 K)) (s_j b_j -
# s'_j b'_j)^2 for each phenotype column, the difference penalty on the
# per-SD effects; in +/- pairs they leave every column's mean at zero.
# @noRd
#' @importFrom Matrix bdiag
.coupledDesign <- function(conds, train, W, layout) {
    copies <- .coupledPlaced(.coupledCopies(conds, train, layout), layout)
    pairs <- .coupledPairs(copies, .coupledCopyStrength(W, copies))
    K <- length(conds$rows)
    nPseudo <- 2L * sum(map_int(pairs, .coupledPairSize))
    N <- sum(map_int(copies, "n")) + nPseudo
    a <- sqrt(N / (K * map_int(copies, "n")))
    columns <- .coupledColumnOrder(copies, layout)
    data <- bdiag(map2(copies, a, .coupledScaledX))[, columns, drop = FALSE]
    pseudo <- .coupledPseudo(pairs, copies, N / (2 * K), ncol(data))
    y <- exec(rbind, !!!map2(copies, a, .coupledScaledY))
    list(
        X = if (is.null(pseudo)) data else rbind(data, pseudo, -pseudo),
        y = .coupledResponse(
            rbind(y, matrix(0, nPseudo, ncol(y)))
        ),
        blocks = .coupledBlocks(copies),
        group = .coupledGroups(copies, ncol(data)),
        partners = .coupledPartners(pairs, ncol(data)),
        K = K
    )
}

# The phenotype column of a design with one, the matrix of a unit's several.
# @noRd
.coupledResponse <- function(y) {
    if (ncol(y) == 1L) y[, 1L] else y
}

# The copies of the variants a design fits: each condition's own, or under
# "units" one per unit, its conditions' phenotypes in the order of their
# columns.
# @noRd
.coupledCopies <- function(conds, train, layout) {
    members <- if (identical(layout, "units")) {
        .coupledUnitMembers(conds$spec)
    } else {
        as.list(seq_along(conds$rows))
    }
    map(members, .coupledCopy, train = train, conds = conds)
}

# Each unit's conditions, in the order of their columns.
# @noRd
.coupledUnitMembers <- function(spec) {
    units <- split(seq_along(spec$unit), spec$unit)
    map(unname(units), .coupledByColumn, column = spec$column)
}

# @noRd
.coupledByColumn <- function(members, column) {
    members[order(column[members])]
}

# The copy of the conditions `members` on their rows of `train`, the first's
# (a unit's conditions share theirs): the centred genotypes over the variants
# that vary there, their SDs, and the centred phenotypes in units of their SD
# (one that does not vary is left unscaled).
# @noRd
.coupledCopy <- function(members, train, conds) {
    rows <- train[[members[[1L]]]]
    Xk <- conds$X[rows, , drop = FALSE]
    xMean <- colMeans(Xk)
    centred <- sweep(Xk, 2L, xMean)
    s <- sqrt(colMeans(centred^2))
    carried <- which(s > sqrt(.Machine$double.eps))
    Yk <- conds$Y[rows, members, drop = FALSE]
    yMean <- colMeans(Yk)
    spread <- sqrt(colMeans(sweep(Yk, 2L, yMean)^2))
    sigma <- replace(spread, spread == 0, 1)
    list(
        X = centred[, carried, drop = FALSE],
        Y = sweep(sweep(Yk, 2L, yMean), 2L, sigma, "/"),
        carried = carried,
        xMean = xMean[carried],
        s = s[carried],
        yMean = unname(yMean),
        sigma = unname(sigma),
        members = members,
        n = length(rows)
    )
}

# The pairs' strengths between copies: a unit's are its first condition's,
# which every column of it shares.
# @noRd
.coupledCopyStrength <- function(W, copies) {
    firsts <- map_int(copies, .coupledFirstMember)
    W[firsts, firsts, drop = FALSE]
}

# @noRd
.coupledFirstMember <- function(copy) {
    copy$members[[1L]]
}

# The copies with their columns in the design (`cols`): copy after copy, or
# under "groups" each variant's copies side by side.
# @noRd
.coupledPlaced <- function(copies, layout) {
    widths <- map_int(copies, .coupledWidth)
    starts <- cumsum(c(0L, widths[-length(widths)]))
    stacked <- map2(starts, widths, .coupledStackedCols)
    if (!identical(layout, "groups")) {
        return(map2(copies, stacked, .coupledPlaceCopy))
    }
    byVariant <- .coupledVariantOrder(copies)
    map2(
        copies,
        map(stacked, .coupledRank, byVariant = byVariant),
        .coupledPlaceCopy
    )
}

# @noRd
.coupledWidth <- function(copy) {
    length(copy$carried)
}

# @noRd
.coupledStackedCols <- function(start, width) {
    start + seq_len(width)
}

# @noRd
.coupledPlaceCopy <- function(copy, cols) {
    list_assign(copy, cols = cols)
}

# @noRd
.coupledRank <- function(cols, byVariant) {
    match(cols, byVariant)
}

# The stacked columns in the order that puts each variant's copies side by
# side, by variant then copy.
# @noRd
.coupledVariantOrder <- function(copies) {
    variant <- list_c(map(copies, "carried"))
    copy <- rep(seq_along(copies), map_int(copies, .coupledWidth))
    order(variant, copy)
}

# The design's column order: the stacked one, or the variants' under
# "groups".
# @noRd
.coupledColumnOrder <- function(copies, layout) {
    stacked <- seq_len(sum(map_int(copies, .coupledWidth)))
    if (!identical(layout, "groups")) {
        return(stacked)
    }
    .coupledVariantOrder(copies)
}

# Each condition's block (.coupledPredict()): its copy's variants, columns
# and centring, its phenotype column's scale and centre, and under a unit
# the columns of its own phenotype.
# @noRd
.coupledBlocks <- function(copies) {
    width <- sum(map_int(copies, .coupledWidth))
    blocks <- list_flatten(map(copies, .coupledCopyBlocks, width = width))
    members <- list_c(map(copies, "members"))
    unname(blocks[order(members)])
}

# @noRd
.coupledCopyBlocks <- function(copy, width) {
    map(
        seq_along(copy$members),
        .coupledColumnBlock,
        copy = copy,
        width = width
    )
}

# @noRd
.coupledColumnBlock <- function(i, copy, width) {
    list(
        carried = copy$carried,
        cols = (i - 1L) * width + copy$cols,
        xMean = copy$xMean,
        yMean = copy$yMean[[i]],
        sigma = copy$sigma[[i]]
    )
}

# Each design column's variant, the group grpreg and abess take.
# @noRd
.coupledGroups <- function(copies, width) {
    replace(
        integer(width),
        list_c(map(copies, "cols")),
        list_c(map(copies, "carried"))
    )
}

# Each pair of coupled copies with the variants both carry, as positions in
# each one's copy, and the pair's strength; a pair sharing none is dropped.
# @noRd
.coupledPairs <- function(copies, W) {
    if (length(copies) < 2L) {
        return(list())
    }
    compact(map(
        utils::combn(length(copies), 2L, simplify = FALSE),
        .coupledPair,
        copies = copies,
        W = W
    ))
}

# @noRd
.coupledPair <- function(pair, copies, W) {
    strength <- W[pair[[1L]], pair[[2L]]]
    first <- copies[[pair[[1L]]]]$carried
    second <- copies[[pair[[2L]]]]$carried
    shared <- intersect(first, second)
    if (strength <= 0 || length(shared) == 0L) {
        return(NULL)
    }
    at1 <- match(shared, first)
    at2 <- match(shared, second)
    list(
        copies = pair,
        strength = strength,
        at1 = at1,
        at2 = at2,
        cols1 = copies[[pair[[1L]]]]$cols[at1],
        cols2 = copies[[pair[[2L]]]]$cols[at2]
    )
}

# @noRd
.coupledPairSize <- function(pair) {
    length(pair$at1)
}

# The pseudo-rows of every coupled pair, before their negated copies, each
# scaled by sqrt(`unit` c); NULL when nothing is coupled.
# @noRd
.coupledPseudo <- function(pairs, copies, unit, nCoef) {
    if (length(pairs) == 0L) {
        return(NULL)
    }
    exec(
        rbind,
        !!!map(
            pairs,
            .coupledPairRows,
            copies = copies,
            unit = unit,
            nCoef = nCoef
        )
    )
}

# One pair's pseudo-rows: one per shared variant, its scale times the
# variant's SD in the first copy and minus its SD in the second.
# @noRd
#' @importFrom Matrix sparseMatrix
.coupledPairRows <- function(pair, copies, unit, nCoef) {
    first <- copies[[pair$copies[[1L]]]]
    second <- copies[[pair$copies[[2L]]]]
    m <- length(pair$at1)
    scale <- sqrt(unit * pair$strength)
    sparseMatrix(
        i = rep(seq_len(m), 2L),
        j = c(first$cols[pair$at1], second$cols[pair$at2]),
        x = scale * c(first$s[pair$at1], -second$s[pair$at2]),
        dims = c(m, nCoef)
    )
}

# Each design column's coupling, its partners' strengths summed.
# @noRd
.coupledPartners <- function(pairs, width) {
    if (length(pairs) == 0L) {
        return(numeric(width))
    }
    columns <- list_c(map(pairs, .coupledPairCols))
    strengths <- list_c(map(pairs, .coupledPairRepeated))
    unname(map_dbl(
        split(strengths, factor(columns, levels = seq_len(width))),
        sum
    ))
}

# @noRd
.coupledPairCols <- function(pair) {
    c(pair$cols1, pair$cols2)
}

# @noRd
.coupledPairRepeated <- function(pair) {
    rep(pair$strength, 2L * length(pair$at1))
}

# @noRd
.coupledScaledX <- function(copy, a) {
    a * copy$X
}

# @noRd
.coupledScaledY <- function(copy, a) {
    a * copy$Y
}

# Whether glmnet's mgaussian family can fit the group penalty with the
# coupling: the conditions fall into units -- those alike along every coupled
# axis -- of two or more, every unit holding the same columns, those alike
# along the other axes, and each unit's conditions measured on the same
# samples.
# @noRd
.coupledUnitsHold <- function(conds) {
    spec <- conds$spec
    members <- .coupledUnitMembers(spec)
    if (max(lengths(members)) < 2L) {
        return(FALSE)
    }
    columns <- map(members, .coupledUnitColumns, column = spec$column)
    all(map_lgl(columns, identical, columns[[1L]])) &&
        all(map_lgl(members, .coupledSharesRows, rows = conds$rows))
}

# @noRd
.coupledUnitColumns <- function(members, column) {
    column[members]
}

# @noRd
.coupledSharesRows <- function(members, rows) {
    all(map_lgl(rows[members], identical, rows[[members[[1L]]]]))
}

# ---- Penalized RSS across coupled conditions (difference penalty) ----------
#
# A group of summary statistics holds its conditions' z-scores over the union
# of their variants, NA where a condition lacks one, each condition on its
# study's LD reference. The penalized RSS methods fit it under the penalty on
# the difference between coupled conditions' per-SD effects, PROSPER's
# coupling: penalizedRssCoupledCpp() solves
#
#   sum_k [b_k'R_k,s b_k - 2 b_k'z_k + penalty(b_k)]
#       + sum_j sum_{k < k'} c_kk' (b_jk - b_jk')^2,
#
# twice the individual-level objective of .penalizedCoupledWeights() on
# standardized data, so c means the same on either; a pair's c is that of
# the axis it differs along (.coupledSpec()). Each condition is fitted
# against its own LD over the grid of shrinkage s, the coupling's cells and
# the shared lambda path, and takes the grid point its own held-out
# statistics score best, by .rssSelectHeldOut()'s rule; each split refits the
# coupled grid on every condition's training part at once. The splits are
# drawn for each condition apart, as if no two shared samples.

# The weights (union variants `ids` x conditions) each condition's own best
# (s, c, lambda) gives it, zero at a variant the condition does not carry.
# `stats` holds one statistics record per condition (z or b with seb, n, and
# variantNames among `ids`, in their spelling), `ldReferences` each
# condition's reference, `spec` the penalty (.rssPenalty()), `coupling` the
# group's coupling (.coupledSpec()) and `options` the path options a
# PenalizedRssOptions holds. With `grouped`, each variant's effects across
# the conditions are one group under the penalty as well; with `choose =
# "pooled"` every condition takes the one grid point scoring best on
# average over them. `residual`, the conditions' residual correlation (K x
# K), correlates the held-out splits of conditions on one LD reference and
# variants (.rssCoupledCorrelate()); NULL draws every condition's apart.
# @noRd
.rssCoupledWeights <- function(
    stats,
    ldReferences,
    ids,
    spec,
    s,
    coupling,
    options,
    fitRetention,
    label,
    grouped = FALSE,
    choose = "own",
    residual = NULL
) {
    conds <- list(spec = .coupledSpec(coupling, length(stats)))
    lds <- .rssCoupledLd(map2(stats, ldReferences, .rssLdSource, label = label))
    walk(lds, .rssCheckHeldOutLd, label = label)
    ns <- map_dbl(stats, .rssCoupledN)
    run <- list(
        lds = lds,
        idx = map(stats, .rssCoupledIndex, ids = ids),
        nUnion = length(ids),
        spec = spec,
        options = .rssCoupledOptions(options, spec),
        cells = .rssCoupledCells(s, conds),
        grouped = grouped,
        choose = choose,
        residual = residual,
        sources = .rssSharedSources(ldReferences, stats)
    )
    cors <- pmap(list(stats, ns, lds), .rssCoupledCor)
    search <- .rssCoupledGrid(
        cors,
        run,
        run$options$lambda %||% .rssCoupledLambda(cors, spec, ns)
    )
    picks <- .rssCoupledSelect(search, cors, ns, run)
    .rssCoupledResult(search, picks, run, ids, fitRetention)
}

# The studies' LD sources in one form for the coupled solver: their factors
# when every study has one, else each one's matrix.
# @noRd
.rssCoupledLd <- function(lds) {
    if (all(map_lgl(lds, .rssHasFactor))) {
        return(lds)
    }
    map(lds, .rssDenseSource)
}

# @noRd
.rssHasFactor <- function(ld) {
    !is.null(ld$ldFactor)
}

# @noRd
.rssDenseSource <- function(ld) {
    list(R = .ldSourceMatrix(ld), nSamples = ld$nSamples)
}

# A study's sample size, the median of its per-variant ones.
# @noRd
.rssCoupledN <- function(stat) {
    stats::median(stat$n)
}

# A study's variants as 1-based positions in the union.
# @noRd
.rssCoupledIndex <- function(stat, ids) {
    match(stat$variantNames, ids)
}

# A study's marginal correlations, clamped as every RSS fit clamps them.
# @noRd
.rssCoupledCor <- function(stat, n, ld) {
    .lassosumClampCor(.rssMarginalCor(stat, n = n, p = .ldSourceSize(ld)))
}

# The path options, penalizedRss()'s own defaults filling what `options`
# leaves unset, with the strategy resolved for the penalty.
# @noRd
.rssCoupledOptions <- function(options, spec) {
    path <- list_modify(
        .formalDefaults(
            penalizedRss,
            c("lambda", "thr", "maxiter", "strategy")
        ),
        !!!options
    )
    list_assign(path, strategy = .rssStrategy(path$strategy, spec$penalty))
}

# The grid's cells, every shrinkage s with every cell of the coupling, in
# the order the candidates are laid out, each with its pairs' strengths
# (`W`, K x K).
# @noRd
.rssCoupledCells <- function(s, conds) {
    cells <- expand_grid(s = s, !!!conds$spec$grids)
    strengths <- map(
        seq_len(nrow(cells)),
        .rssCellStrength,
        cells = cells,
        conds = conds
    )
    mutate(cells, W = strengths)
}

# @noRd
.rssCellStrength <- function(i, cells, conds) {
    values <- unlist(cells[i, setdiff(names(cells), "s")])
    .coupledPairStrengths(c(numeric(0), values), conds)
}

# The shared lambda path, as .rssDefaultLambda() makes one study's: from the
# largest any study starts at, down by the ratio for the fewest samples any
# study has against its variants. None for the L0 family, whose solver takes
# a geometric one from the joint start.
# @noRd
.rssCoupledLambda <- function(cors, spec, ns) {
    if (.rssIsL0(spec$penalty)) {
        return(numeric(0))
    }
    control <- .rssPathControl()
    few <- spec$alpha != 0 && any(ns <= map_int(cors, length))
    ratio <- if (few) control$ratioFewSamples else control$ratio
    top <- max(map_dbl(cors, .rssAbsMax)) / max(spec$alpha, 1e-3)
    top * ratio^seq(0, 1, length.out = control$points)
}

# @noRd
.rssAbsMax <- function(x) {
    max(abs(x), 0)
}

# The coupled grid on the studies' correlations `cors`: one path per cell, at
# `path` or at each cell's own `grids`. Returns each study's candidates
# (`beta`, a matrix per study, cells' paths side by side), the cells' path
# values (`grids`) and the candidates' table (`meta`).
# @noRd
.rssCoupledGrid <- function(cors, run, path = NULL, grids = NULL) {
    fits <- pmap(
        list(
            run$cells$s,
            run$cells$W,
            grids %||% rep(list(path), nrow(run$cells))
        ),
        .rssCoupledFit,
        cors = cors,
        run = run
    )
    list(
        beta = map(seq_along(cors), .rssCoupledStudyBeta, fits = fits),
        grids = map(fits, "lambda"),
        meta = list_rbind(map2(
            seq_along(fits),
            fits,
            .rssCoupledCellMeta,
            cells = run$cells
        ))
    )
}

# One cell's coupled path, through the factor solver when the studies' LD
# sources are factors.
# @noRd
.rssCoupledFit <- function(sVal, W, path, cors, run) {
    spec <- run$spec
    settings <- list(
        zList = cors,
        unionIdxList = map(run$idx, as.integer),
        nUnion = as.integer(run$nUnion),
        strengthR = W,
        grouped = run$grouped,
        shrink = sVal,
        pathR = path,
        penaltyStr = spec$penalty,
        strategyStr = run$options$strategy,
        gamma = spec$gamma,
        alpha = spec$alpha,
        lambda1 = spec$lambda1,
        lambda2 = spec$lambda2,
        thr = run$options$thr,
        maxiter = as.integer(run$options$maxiter),
        maxSwaps = as.integer(spec$maxSwaps),
        control = .rssPathControl()
    )
    if (.rssHasFactor(run$lds[[1L]])) {
        return(exec(
            penalizedRssCoupledFactorCpp,
            factorList = map(run$lds, "ldFactor"),
            !!!settings
        ))
    }
    exec(penalizedRssCoupledCpp, ldList = map(run$lds, "R"), !!!settings)
}

# Study `k`'s candidates across the cells.
# @noRd
.rssCoupledStudyBeta <- function(k, fits) {
    exec(cbind, !!!map(fits, list("studies", k, "beta")))
}

# One cell's candidates: its s and coupling strengths, each lambda, the
# convergence flag and the coupling term there.
# @noRd
.rssCoupledCellMeta <- function(i, fit, cells) {
    lambda <- as.numeric(fit$lambda)
    bind_cols(
        cells[rep(i, length(lambda)), setdiff(names(cells), "W")],
        tibble(
            lambda = lambda,
            conv = as.integer(fit$studies[[1L]]$conv),
            coupling = as.numeric(fit$coupling)
        )
    )
}

# Each study's pick by held-out statistics (.rssSelectHeldOut()): its own
# splits, drawn from its own LD, with the coupled grid refit on every
# study's training part at once and each study scored on its own held-out
# part.
# @noRd
#' @importFrom stats rnorm
.rssCoupledSelect <- function(search, cors, ns, run) {
    nSplits <- .rssHeldOutSplits()
    ws <- map(run$lds, .rssCoupledIntensity)
    unit <- pmap(
        list(run$lds, ws, cors),
        .rssCoupledUnitNoise,
        nSplits = nSplits
    )
    noise <- map2(.rssCoupledCorrelate(unit, run), ns, .rssCoupledScaleNoise)
    scores <- map(
        seq_len(nSplits),
        .rssCoupledSplitScore,
        noise = noise,
        cors = cors,
        ws = ws,
        run = run,
        grids = search$grids
    )
    picked <- map(
        seq_along(cors),
        .rssCoupledScores,
        scores = scores,
        search = search
    )
    if (identical(run$choose, "pooled")) {
        return(rep(list(.rssPooledPick(picked)), length(cors)))
    }
    map(picked, .rssScoresPick)
}

# The one candidate every condition takes: the best on its held-out scores
# averaged over the conditions.
# @noRd
.rssPooledPick <- function(picked) {
    .rssScoresPick(list(
        cor = reduce(map(picked, "cor"), `+`) / length(picked),
        explained = reduce(map(picked, "explained"), `+`) / length(picked)
    ))
}

# @noRd
.rssScoresPick <- function(scores) {
    .rssHeldOutPick(scores$cor, scores$explained)
}

# @noRd
.rssCoupledIntensity <- function(ld) {
    .ldShrinkIntensity(ld, ld$nSamples)
}

# A condition's splits' noise in unit variance, as .rssSelectHeldOut() draws
# it, before its scale (.rssCoupledScaleNoise()).
# @noRd
.rssCoupledUnitNoise <- function(ld, w, cor, nSplits) {
    sqrt(1 - w) *
        .ldSourceDraws(ld, nSplits) +
        sqrt(w) * matrix(rnorm(length(cor) * nSplits), ncol = nSplits)
}

# @noRd
.rssCoupledScaleNoise <- function(draws, n) {
    share <- .rssHeldOutTrainShare()
    draws * sqrt((1 - share) / (share * n))
}

# The conditions' unit noise, those sharing an LD reference and variants
# (`run$sources`) correlated by their residual correlation (`run$residual`):
# each a combination of the source's independent draws by a square root S of
# the correlation (S S' = C), so its draws keep their variance and a pair's
# correlate at C's entry, as a split of samples shared by both would. Every
# condition's draws stay apart without a correlation.
# @noRd
.rssCoupledCorrelate <- function(unit, run) {
    if (is.null(run$residual)) {
        return(unit)
    }
    S <- .rssBlockRoot(run$residual, run$sources)
    map(seq_along(unit), .rssMixedNoise, unit = unit, S = S)
}

# A square root of the residual correlation within each shared source, the
# identity between sources.
# @noRd
.rssBlockRoot <- function(residual, sources) {
    within <- outer(sources, sources, `==`)
    C <- replace(residual, !within, 0)
    ev <- eigen(C, symmetric = TRUE)
    ev$vectors %*% (sqrt(pmax(ev$values, 0)) * t(ev$vectors))
}

# @noRd
.rssMixedNoise <- function(k, unit, S) {
    partners <- which(S[k, ] != 0)
    reduce(map2(unit[partners], S[k, partners], `*`), `+`)
}

# Each condition's source: the first condition on the same LD reference over
# the same variants, whose splits can be shared.
# @noRd
.rssSharedSources <- function(ldReferences, stats) {
    map_int(
        seq_along(stats),
        .rssFirstSharing,
        ldReferences = ldReferences,
        stats = stats
    )
}

# @noRd
.rssFirstSharing <- function(k, ldReferences, stats) {
    detect_index(
        seq_len(k),
        .rssSharesSource,
        k = k,
        ldReferences = ldReferences,
        stats = stats
    )
}

# @noRd
.rssSharesSource <- function(l, k, ldReferences, stats) {
    identical(ldReferences[[l]], ldReferences[[k]]) &&
        identical(stats[[l]]$variantNames, stats[[k]]$variantNames)
}

# Split `i`: the coupled grid refit on each study's training part and each
# study's candidates scored on its own held-out part.
# @noRd
.rssCoupledSplitScore <- function(i, noise, cors, ws, run, grids) {
    share <- .rssHeldOutTrainShare()
    train <- map2(cors, noise, .rssCoupledTrainCor, i = i)
    held <- map2(cors, train, .rssCoupledHeldCor, share = share)
    refit <- .rssCoupledGrid(train, run, grids = grids)
    pmap(list(refit$beta, held, run$lds, ws), .rssCoupledMoments)
}

# @noRd
.rssCoupledTrainCor <- function(cor, noise, i) {
    .lassosumClampCor(cor + noise[, i])
}

# @noRd
.rssCoupledHeldCor <- function(cor, train, share) {
    (cor - share * train) / (1 - share)
}

# A study's candidates' correlation with its held-out part and the variance
# of it they explain (.rssHeldOutScore()).
# @noRd
.rssCoupledMoments <- function(beta, held, ld, w) {
    moments <- .rssShrunkLdMoments(beta, held, ld, w)
    list(
        cor = .rssCorrelationScore(moments),
        explained = 2 * moments$bxy - moments$bxxb
    )
}

# Condition `k`'s candidates' held-out scores (.rssHeldOutPick() reads
# them), averaged over the splits: one empty on the full data has nothing to
# return.
# @noRd
.rssCoupledScores <- function(k, scores, search) {
    cor <- rowMeans(exec(cbind, !!!map(scores, list(k, "cor"))))
    explained <- rowMeans(exec(cbind, !!!map(scores, list(k, "explained"))))
    empty <- colSums(search$beta[[k]] != 0) == 0
    list(cor = replace(cor, empty, -Inf), explained = explained)
}

# The weights over the union, each study's at its pick, with the search kept
# as `fitRetention` asks: "slim" each study's grid point, "full" every
# candidate too.
# @noRd
.rssCoupledResult <- function(search, picks, run, ids, fitRetention) {
    cols <- pmap(
        list(search$beta, picks, run$idx),
        .rssCoupledStudyWeights,
        nUnion = length(ids)
    )
    W <- `rownames<-`(exec(cbind, !!!cols), ids)
    if (identical(fitRetention, "none")) {
        return(W)
    }
    chosen <- search$meta[unlist(picks), , drop = FALSE]
    fit <- if (identical(fitRetention, "slim")) {
        list(selection = chosen)
    } else {
        list(selection = chosen, beta = search$beta, meta = search$meta)
    }
    `attr<-`(W, "fit", .withFitVariantIds(fit, ids))
}

# @noRd
.rssCoupledStudyWeights <- function(beta, pick, idx, nUnion) {
    replace(numeric(nUnion), idx, beta[, pick])
}

# ---- Samplers across studies: PRS-CSx / JointPRS and SDPRX ------------------
#
# A study-axis group of summary statistics fitted by one sampler across its
# studies, each study against its own LD reference, over the union of the
# studies' variants. PRS-CS runs as JointPRS fits ancestries (PrsCsChain in
# rss_mcmc.h): the studies share each variant's local shrinkage and the
# global shrinkage, as PRS-CSx shares them, and their effects on a variant
# are correlated by rho, learned on [0, 0.99]. SDPR runs as SDPRX fits them
# (the SDPR chain in rss_mcmc.cpp): each non-null component draws a
# variant's effects in a set of the studies, correlated by rho.

# The studies' inputs to a sampler across studies: their LD sources in one
# form, sample sizes, variants' places in the union of size `nUnion`, and
# marginal correlations.
# @noRd
.mcmcCoupledInputs <- function(stats, ldReferences, ids, label) {
    lds <- .rssCoupledLd(map2(stats, ldReferences, .rssLdSource, label = label))
    ns <- map_dbl(stats, .rssCoupledN)
    list(
        lds = lds,
        ns = ns,
        idx = map(stats, .rssCoupledIndex, ids = ids),
        bhats = pmap(list(stats, ns, lds), .coupledMarginalCor),
        nUnion = length(ids)
    )
}

# A study's marginal correlations, as prsCsWeights() and sdprWeights() read
# them.
# @noRd
.coupledMarginalCor <- function(stat, n, ld) {
    .rssMarginalCor(stat, n = n, p = .ldSourceSize(ld))
}

# A study's values over the union, zero at a variant it does not carry.
# @noRd
.coupledUnionColumn <- function(values, idx, nUnion) {
    replace(numeric(nUnion), idx, values)
}

# The studies' marginal correlations over the union (union x studies), the
# bound a single chain's effects are checked against.
# @noRd
.coupledUnionMatrix <- function(run) {
    exec(
        cbind,
        !!!map2(run$bhats, run$idx, .coupledUnionColumn, nUnion = run$nUnion)
    )
}

# A sampler's weights across studies (union variants `ids` x studies), with
# its estimates as the "fit" attribute unless `fitRetention` is "none".
# @noRd
.mcmcCoupledResult <- function(fit, ids, fitRetention) {
    W <- `dimnames<-`(fit$betaEst, list(ids, NULL))
    if (identical(fitRetention, "none")) {
        return(W)
    }
    `attr<-`(W, "fit", .withFitVariantIds(fit, ids))
}

# PRS-CS's weights across studies (union variants `ids` x studies), zero at a
# variant a study does not carry. `stats` and `ldReferences` are as
# .rssCoupledWeights() takes them, and `user` the method's settings, as
# PrsCsOptions() holds them.
# @noRd
.prsCsCoupledWeights <- function(
    stats,
    ldReferences,
    ids,
    user,
    fitRetention,
    label,
    shared = list()
) {
    settings <- .prsCsCoupledSettings(user, label)
    fixed <- .sdprCoupledRho(settings$rho, length(stats), label)
    run <- c(
        list(
            settings = settings,
            rho = fixed %||% .samplerSeedRho(shared$seed),
            learnRho = is.null(fixed),
            noise = .samplerNoisePrecision(shared$residual, ldReferences, stats)
        ),
        .mcmcCoupledInputs(stats, ldReferences, ids, label)
    )
    fit <- .mcmcChains(
        ".jointPrsChain",
        list(run = run),
        .mcmcSeeds(settings$seed, settings$nChains),
        c("betaEst", "psiEst", "sigmaEst", "phiEst", "rhoEst", "residualEst"),
        label,
        .coupledUnionMatrix(run)
    )
    .mcmcCoupledResult(fit, ids, fitRetention)
}

# The sampler's settings, prsCs()'s defaults filling what `user` leaves unset.
# The weights across studies stay on the standardized scale, so a `maf` to put
# one study's on its alleles' is refused.
# @noRd
.prsCsCoupledSettings <- function(user, label) {
    if (!is.null(user$maf)) {
        abort(c(
            glue(
                "{label}: `maf` puts one study's weights on its alleles' ",
                "scale, and the weights across studies stay standardized."
            ),
            i = "Leave `maf` unset, or pass coupling = NULL."
        ))
    }
    honoured <- c(
        "a",
        "b",
        "phi",
        "nIter",
        "nBurnin",
        "thin",
        "verbose",
        "seed",
        "shrinkCorr",
        "raoBlackwell",
        "nChains",
        "rho"
    )
    settings <- list_modify(
        .formalDefaults(prsCs, honoured),
        !!!user[intersect(names(user), honoured)]
    )
    .mcmcAssertShared(
        settings$shrinkCorr,
        settings$raoBlackwell,
        settings$nChains
    )
    settings
}

# One chain of the sampler across conditions, on their LD matrices or
# factors: their correlations learned (a constraint of 1) from the identity or
# `run$rho`, or held at `run$rho` (`run$learnRho`), and the noise of those on
# shared samples correlated by `run$noise` (its precision, NULL for none).
# @noRd
.jointPrsChain <- function(run, seed, track) {
    s <- run$settings
    args <- list(
        a = s$a,
        b = s$b,
        phi = s$phi,
        constraint = rep(1, length(run$bhats)),
        rho = run$rho,
        learnRho = run$learnRho,
        noisePrecision = run$noise,
        bhatList = run$bhats,
        n = run$ns,
        unionIdxList = map(run$idx, as.integer),
        nUnion = as.integer(run$nUnion),
        nIter = as.integer(s$nIter),
        nBurnin = as.integer(s$nBurnin),
        thin = as.integer(s$thin),
        verbose = s$verbose,
        seed = seed,
        shrinkCorr = s$shrinkCorr,
        raoBlackwell = s$raoBlackwell,
        track = track
    )
    if (.rssHasFactor(run$lds[[1L]])) {
        return(exec(
            jointPrsFactorCpp,
            factorList = map(run$lds, "ldFactor"),
            !!!args
        ))
    }
    exec(jointPrsCpp, ldList = map(run$lds, "R"), !!!args)
}

# The settings of sdpr() that act across studies only.
# @noRd
.SDPR_ACROSS <- c("sharing", "rho")

# SDPR's weights across studies (union variants `ids` x studies), zero at a
# variant a study does not carry, as .prsCsCoupledWeights() takes its
# arguments, with `user` as SdprOptions() holds the settings.
# @noRd
.sdprCoupledWeights <- function(
    stats,
    ldReferences,
    ids,
    user,
    fitRetention,
    label,
    shared = list()
) {
    settings <- .sdprCoupledSettings(user, length(stats), label)
    run <- c(
        list(
            settings = settings,
            rho = settings$rho %||% .samplerSeedRho(shared$seed),
            learnRho = is.null(settings$rho),
            kindWeights = .sdprSeedKindWeights(shared$seed, settings$kinds),
            noise = .samplerNoisePrecision(shared$residual, ldReferences, stats)
        ),
        .mcmcCoupledInputs(stats, ldReferences, ids, label)
    )
    fit <- .mcmcChains(
        ".sdprxChain",
        list(run = run),
        .mcmcSeeds(settings$seed, settings$nChains),
        c("betaEst", "h2", "rhoEst"),
        label,
        .coupledUnionMatrix(run)
    )
    .mcmcCoupledResult(fit, ids, fitRetention)
}

# The sampler's settings, sdpr()'s defaults filling what `user` leaves unset,
# with the kinds `sharing` gives `K` studies and rho as a matrix (NULL to
# learn it). Per-variant sample sizes, arrays and optLlk = 2 describe one
# study's variants, so they are refused.
# @noRd
.sdprCoupledSettings <- function(user, K, label) {
    perVariant <- c("perVariantSampleSize", "array")
    given <- names(compact(user[intersect(names(user), perVariant)]))
    llk2 <- !is.null(user$optLlk) && user$optLlk != 1
    if (length(given) > 0L || llk2) {
        abort(c(
            glue(
                "{label}: `perVariantSampleSize`, `array` and `optLlk = 2` ",
                "describe one study's variants; a fit across studies reads ",
                "each study's sample size and the optLlk = 1 likelihood."
            ),
            i = "Leave them unset, or pass coupling = NULL."
        ))
    }
    honoured <- c(
        "a",
        "c",
        "M",
        "a0k",
        "b0k",
        "iter",
        "burn",
        "thin",
        "verbose",
        "seed",
        "raoBlackwell",
        "nChains",
        .SDPR_ACROSS
    )
    settings <- list_modify(
        .formalDefaults(sdpr, honoured),
        !!!user[intersect(names(user), honoured)]
    )
    .sdprValidate(settings$M, NULL, NULL)
    .mcmcAssertShared(1, settings$raoBlackwell, settings$nChains)
    sharing <- arg_match0(
        settings$sharing,
        .formalDefaults(sdpr, "sharing")[[1L]],
        "sharing"
    )
    c(
        settings[setdiff(names(settings), .SDPR_ACROSS)],
        list(
            kinds = .sdprKinds(sharing, K),
            rho = .sdprCoupledRho(settings$rho, K, label)
        )
    )
}

# The kinds of non-null component `sharing` gives K studies, each the
# studies its effects lie in: each alone and all together, every nonempty
# set by size, or all together only.
# @noRd
#' @importFrom utils combn
.sdprKinds <- function(sharing, K) {
    switch(
        sharing,
        eachOrAll = unique(c(as.list(seq_len(K)), list(seq_len(K)))),
        everySubset = list_flatten(
            map(seq_len(K), combn, x = K, simplify = FALSE)
        ),
        all = list(seq_len(K))
    )
}

# rho as the sampler takes it: NULL to learn it, else the K x K correlation
# matrix a number in [0, 1) gives for every pair, or the one given.
# @noRd
.sdprCoupledRho <- function(rho, K, label) {
    if (is.null(rho)) {
        return(NULL)
    }
    if (is.numeric(rho) && length(rho) == 1L && is.null(dim(rho))) {
        if (is.na(rho) || rho < 0 || rho >= 1) {
            abort(glue("{label}: a single `rho` must lie in [0, 1)."))
        }
        return(`diag<-`(matrix(rho, K, K), 1))
    }
    if (!.sdprIsCorrelation(rho, K)) {
        abort(glue(
            "{label}: `rho` must be NULL, a number in [0, 1) or a {K} x {K} ",
            "positive-definite correlation matrix."
        ))
    }
    unname((rho + t(rho)) / 2)
}

# Whether `rho` is a K x K positive-definite correlation matrix.
# @noRd
.sdprIsCorrelation <- function(rho, K) {
    is.matrix(rho) &&
        is.numeric(rho) &&
        identical(dim(rho), c(K, K)) &&
        !anyNA(rho) &&
        isSymmetric(unname(rho)) &&
        all(diag(rho) == 1) &&
        min(eigen(rho, symmetric = TRUE, only.values = TRUE)$values) > 0
}

# The precision of the conditions' noise the samplers read (NoiseCoupling in
# rss_mcmc.h): Q = E^-1 for E the residual correlation over conditions on
# one LD reference and the same variants (.rssSharedSources()), zero between
# others; NULL when no two share samples. solve() leaves the inverse
# asymmetric by rounding, and the samplers take only an exactly symmetric
# one, so it is symmetrized.
# @noRd
.samplerNoisePrecision <- function(residual, ldReferences, stats) {
    if (is.null(residual)) {
        return(NULL)
    }
    sources <- .rssSharedSources(ldReferences, stats)
    E <- replace(residual, !outer(sources, sources, `==`), 0)
    if (all(E[upper.tri(E)] == 0)) {
        return(NULL)
    }
    Q <- solve(E)
    unname((Q + t(Q)) / 2)
}

# A MashPrior's starting correlation for a sampler: that of its covariances
# averaged over their weights, sum_k w_k U_k, on the group's conditions
# (`seed`: U and w); NULL without a seed, or where it is not a positive-
# definite correlation.
# @noRd
.samplerSeedRho <- function(seed) {
    S <- .samplerSeedCovariance(seed)
    if (is.null(S) || any(diag(S) <= 0)) {
        return(NULL)
    }
    R <- stats::cov2cor(S)
    if (min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) <= 1e-8) {
        return(NULL)
    }
    unname(R)
}

# sum_k w_k U_k over a seed's covariances, equal weights where it has none.
# @noRd
.samplerSeedCovariance <- function(seed) {
    if (is.null(seed) || length(seed$U) == 0L) {
        return(NULL)
    }
    w <- .samplerSeedWeights(seed)
    reduce(map2(seed$U, w, `*`), `+`)
}

# @noRd
.samplerSeedWeights <- function(seed) {
    w <- seed$w[names(seed$U)] %||% seed$w
    if (is.null(w) || length(w) != length(seed$U) || anyNA(w)) {
        return(rep(1 / length(seed$U), length(seed$U)))
    }
    unname(w)
}

# SDPR's kinds' starting weights from a MashPrior seed: each covariance's
# weight to the kind holding the conditions it covers (its diagonal's
# non-zero entries), a covariance covering no kind's set adding to none; each
# kind at least a small share, so every one can be drawn. NULL without a
# seed.
# @noRd
.sdprSeedKindWeights <- function(seed, kinds) {
    if (is.null(seed) || length(seed$U) == 0L) {
        return(NULL)
    }
    covered <- map(seed$U, .samplerCovered)
    at <- map_int(covered, .sdprKindOf, kinds = kinds)
    w <- .samplerSeedWeights(seed)
    total <- map_dbl(seq_along(kinds), .sdprKindWeight, at = at, w = w)
    floor <- 1e-3 * max(sum(total), 1)
    pmax(total, floor)
}

# @noRd
.samplerCovered <- function(U) {
    which(diag(U) > 1e-8 * max(diag(U), 1e-300))
}

# @noRd
.sdprKindOf <- function(covered, kinds) {
    hit <- detect_index(kinds, setequal, covered)
    if (hit == 0L) NA_integer_ else hit
}

# @noRd
.sdprKindWeight <- function(q, at, w) {
    sum(w[!is.na(at) & at == q])
}

# One chain of SDPR across conditions, on their LD matrices or factors: rho
# fixed or learned from `run$rho`, the kinds starting at `run$kindWeights`,
# and the noise of conditions on shared samples correlated by `run$noise`.
# @noRd
.sdprxChain <- function(run, seed, track) {
    s <- run$settings
    args <- list(
        bhatList = run$bhats,
        n = as.integer(run$ns),
        unionIdxList = map(run$idx, as.integer),
        nUnion = as.integer(run$nUnion),
        kinds = map(s$kinds, as.integer),
        rho = run$rho,
        learnRho = run$learnRho,
        kindWeights = run$kindWeights,
        noisePrecision = run$noise,
        a = s$a,
        c = s$c,
        M = as.integer(s$M),
        a0k = s$a0k,
        b0k = s$b0k,
        iter = as.integer(s$iter),
        burn = as.integer(s$burn),
        thin = as.integer(s$thin),
        verbose = s$verbose,
        seed = seed,
        raoBlackwell = s$raoBlackwell,
        track = track
    )
    if (.rssHasFactor(run$lds[[1L]])) {
        return(exec(
            sdprxFactorCpp,
            factorList = map(run$lds, "ldFactor"),
            !!!args
        ))
    }
    exec(sdprxCpp, ldList = map(run$lds, "R"), !!!args)
}

#' Compute Weights Using a BGLR Linear Regression Model
#'
#' Internal helper that fits a `BGLR::BGLR` linear regression with a single
#' linear term whose `model` is one of BGLR's marker-effect priors (e.g.
#' "BayesB", "BL"), then returns the posterior mean of the marker effects.
#' Given a matrix `y`, one column per condition, it fits the conditions
#' jointly with `BGLR::Multitrait` instead, under the multivariate form of
#' `model`: "SpikeSlab" for "BayesC" and "BRR" for "BRR", the only two
#' Multitrait() has. A missing phenotype is drawn within its sampler. BGLR
#' writes per-call temporary files to disk; this helper sandboxes them in a
#' fresh `tempdir()` that is cleaned up on exit.
#'
#' @param X A numeric matrix of predictors.
#' @param y A numeric response vector, or a matrix with one column per
#'   condition (missing values allowed).
#' @param model A BGLR marker-effect model name (e.g. "BayesB" or "BL").
#' @param nIter Number of MCMC iterations.
#' @param burnIn Number of burn-in iterations.
#' @param thin Thinning interval.
#' @param etaArgs Optional named list of additional arguments included in the
#'   `ETA` linear-term specification (e.g. `list(probIn = 0.05)` for BayesB).
#'   For a matrix `y`, a `probIn` starts every condition's inclusion
#'   probability.
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights: \code{"none"} (default) attaches
#'   nothing, \code{"slim"} and \code{"full"} the fitted \code{BGLR} object,
#'   which keeps no copy of the data, so the two retaining levels behave
#'   alike here.
#' @param methodArgs Options for \code{BGLR::BGLR} or
#'   \code{BGLR::Multitrait}, built with \code{\link{BglrOptions}}; each call
#'   receives the options it takes, and \code{Cov} goes into Multitrait()'s
#'   linear term. A bare list is refused: it cannot be checked against the
#'   engine, so a misspelled option would be silently ignored.
#' @return A numeric vector of length `ncol(X)` of variant weights, or for a
#'   matrix `y` a matrix of them with one column per condition, with the
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
    keep <- .dropZeroVariance(X, "bglrWeights")
    tmpdir <- tempfile("bglr_")
    dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
    run <- list(
        nIter = nIter,
        burnIn = burnIn,
        thin = thin,
        saveAt = str_c(tmpdir, .Platform$file.sep),
        verbose = FALSE
    )
    XKeep <- X[, keep, drop = FALSE]
    multitrait <- is.matrix(y) && ncol(y) > 1L
    fit <- if (multitrait) {
        .bglrMultitraitFit(XKeep, y, model, etaArgs, methodArgs, run)
    } else {
        .bglrFit(XKeep, y, model, etaArgs, methodArgs, run)
    }
    coef <- if (multitrait) fit$ETA[[1L]]$beta else as.numeric(fit$ETA[[1L]]$b)
    weights <- .weightsOnColumns(coef, keep, colnames(X), colnames(y))
    if (identical(fitRetention, "none")) {
        return(weights)
    }
    `attr<-`(weights, "fit", .withFitVariantIds(fit, colnames(X)[keep]))
}

# bayesC or bayesN (`method`) on BGLR fitted across a joint group's
# conditions at once (Multitrait()), as the pipeline fits them on a group of
# two or more (.jointTwasMultitrait()): `Y` has a column per condition,
# missing where a condition's phenotype was not measured. Taking `Y` is what
# sends the pipeline's fits, and each cross-validation fold's, down the
# whole-Y path (.twasFitsWholeY()). `studies` names each condition's study:
# the residuals are correlated within a study and independent across
# (.bglrStudyResiduals()).
# @noRd
.bglrJointWeights <- function(
    X,
    Y,
    method,
    studies = NULL,
    backend = "bglr",
    pi = NULL,
    h2Estimate = NULL,
    nIter = 5000,
    nBurnin = 1000,
    thin = 5,
    fitRetention = c("none", "slim", "full"),
    methodArgs = BglrOptions()
) {
    bayesAlphabetWeights(
        X,
        Y,
        method,
        backend = backend,
        pi = pi,
        h2Estimate = h2Estimate,
        nIter = nIter,
        nBurnin = nBurnin,
        thin = thin,
        fitRetention = arg_match(fitRetention),
        methodArgs = .bglrStudyResiduals(methodArgs, studies)
    )
}

# The residual covariance of a joint fit over conditions in the studies
# `studies`: `resCov` as given within one study; across studies, whose
# samples are disjoint, independent. Each condition its own study makes it
# diagonal; several conditions in some study make it REC (Multitrait()'s
# recursive structure, R = (I - B)^-1 diag(psi) (I - B)^-T), its recursions
# B allowed only between conditions of one study, so R is free within each
# study and zero across. A diagonal `resCov` is kept as it is.
# @noRd
.bglrStudyResiduals <- function(methodArgs, studies) {
    settings <- as.list(methodArgs)
    if (n_distinct(studies) <= 1L || identical(settings$resCov$type, "DIAG")) {
        return(methodArgs)
    }
    resCov <- if (anyDuplicated(studies) == 0L) {
        list(type = "DIAG")
    } else {
        list(type = "REC", M = .bglrWithinStudy(studies))
    }
    exec(BglrOptions, !!!list_assign(settings, resCov = resCov))
}

# REC's recursions: TRUE where condition i (a row) may depend on an earlier
# condition j of the same study.
# @noRd
.bglrWithinStudy <- function(studies) {
    same <- outer(studies, studies, `==`)
    same & lower.tri(same)
}

# BGLR() on one phenotype, its one linear term `model` with `etaArgs`. Of
# the options only those BGLR() takes reach it: BglrOptions() also holds
# Multitrait()'s and the covariance of a joint fit's linear term.
# @noRd
.bglrFit <- function(X, y, model, etaArgs, methodArgs, run) {
    own <- .bglrOwnArgs(methodArgs, BGLR::BGLR)
    callArgs <- list_modify(
        c(
            list(
                y = if (is.matrix(y)) y[, 1L] else y,
                ETA = list(c(list(X = X, model = model), etaArgs))
            ),
            run
        ),
        !!!own
    )
    exec(BGLR::BGLR, !!!callArgs)
}

# BGLR's Multitrait() on a matrix of phenotypes, one column per condition,
# under the multivariate form of `model` (.bglrMultitraitModels), its linear
# term holding the covariance across conditions `methodArgs` gives
# (.bglrMultitraitTerm()). Only the options Multitrait() takes reach it. It
# prints its priors whatever `verbose` says, so the call is quieted.
# @noRd
.bglrMultitraitFit <- function(X, Y, model, etaArgs, methodArgs, run) {
    if (!is_in(model, names(.bglrMultitraitModels))) {
        abort(glue(
            "bglrWeights: BGLR's Multitrait() has no form of {model}; a ",
            "matrix `y` takes ",
            "{str_flatten(names(.bglrMultitraitModels), ' or ')}."
        ))
    }
    mtModel <- .bglrMultitraitModels[[model]]
    term <- c(
        list(X = X, model = mtModel),
        .bglrMultitraitTerm(mtModel, etaArgs, methodArgs, ncol(Y))
    )
    own <- .bglrOwnArgs(methodArgs, BGLR::Multitrait)
    callArgs <- list_modify(c(list(y = Y, ETA = list(term)), run), !!!own)
    .quietEval(exec(BGLR::Multitrait, !!!callArgs))
}

# The settings of a Multitrait() linear term: those `methodArgs` gives it
# (.bglrTermSettings()), and for a spike and slab a prior inclusion
# probability (`probIn`) starting every condition, with BGLR's default prior
# counts (2): Multitrait() fails when `probIn` is given without them.
# @noRd
.bglrMultitraitTerm <- function(mtModel, etaArgs, methodArgs, nCond) {
    termArgs <- .bglrOwnArgs(methodArgs, .bglrTermSettings)
    inclusion <- if (
        identical(mtModel, "SpikeSlab") && !is.null(etaArgs$probIn)
    ) {
        list(probIn = rep_len(etaArgs$probIn, nCond), counts = rep(2, nCond))
    }
    c(
        exec(.bglrTermSettings, !!!termArgs),
        compact(list(inclusionProb = inclusion))
    )
}

# The settings of a Multitrait() linear term BglrOptions() takes, as BGLR
# spells them there: `Cov`, the covariance of the effects across
# conditions. Its formals are what BglrOptions() and the pipeline's check of
# a method's options (.alphabetEngineCallees) read them from.
# @noRd
.bglrTermSettings <- function(Cov = NULL) {
    compact(list(Cov = Cov))
}

# The options in the record `methodArgs` that the BGLR function `fn` takes.
# @noRd
#' @importFrom rlang fn_fmls_names
.bglrOwnArgs <- function(methodArgs, fn) {
    user <- as.list(methodArgs)
    user[is_in(names(user), fn_fmls_names(fn))]
}

#' Compute Weights Using Dirichlet Process Regression (RcppDPR)
#'
#' Fits a Dirichlet Process Regression model via `RcppDPR::fit_model` and
#' returns its per-variant effects `beta` as the weights. Without a kinship
#' rotation, `fit_model`'s `alpha` is not a genetic effect (it maps the
#' in-sample residual back onto the genotypes), so it is left out.
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
#'   nothing, \code{"slim"} and \code{"full"} the fitted model without its
#'   copy of the genotypes (\code{x}). Each names the variants it was fitted
#'   on.
#' @param retainX Logical. Keep that copy of the genotypes in a \code{"full"}
#'   fit (default \code{FALSE}). An error at any other \code{fitRetention}.
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
    retainX = FALSE,
    methodArgs = DprOptions()
) {
    fitRetention <- arg_match(fitRetention)
    .assertRetainX(retainX, fitRetention)
    .assertMethodOptions(methodArgs, "DprOptions", "methodArgs")
    if (!requireNamespace("RcppDPR", quietly = TRUE)) {
        abort("Package 'RcppDPR' is required for this function.")
    }
    zeros <- rep(0, ncol(X))
    keep <- .dropZeroVariance(X, "dprWeights")
    w <- matrix(1, nrow = nrow(X), ncol = 1)
    user <- as.list(methodArgs)
    callArgs <- list_modify(
        list(
            y = y,
            w = w,
            x = X[, keep, drop = FALSE],
            rotate_variables = FALSE,
            fitting_method = fittingMethod
        ),
        !!!user
    )
    fit <- exec(RcppDPR::fit_model, !!!callArgs)
    eff.wgt <- replace(zeros, keep, as.numeric(fit$beta))
    if (identical(fitRetention, "none")) {
        return(eff.wgt)
    }
    kept <- if (retainX) fit else .withoutFields(fit, "x")
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
    retainX = FALSE,
    methodArgs = DprOptions()
) {
    fitRetention <- arg_match(fitRetention)
    .assertMethodOptions(methodArgs, "DprOptions", "methodArgs")
    dprWeights(
        X,
        y,
        fittingMethod = "VB",
        fitRetention = fitRetention,
        retainX = retainX,
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
#' # show_progress = FALSE turns off RcppDPR's progress bar.
#' dprGibbsWeights(X, y, sStep = 500,
#'   methodArgs = DprOptions(show_progress = FALSE))
#' @importFrom checkmate assertCount assertFlag
#' @export
dprGibbsWeights <- function(
    X,
    y,
    sStep = 5000,
    fitRetention = c("none", "slim", "full"),
    retainX = FALSE,
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
        retainX = retainX,
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
#' dprAdaptiveGibbsWeights(X, y,
#'   methodArgs = DprOptions(show_progress = FALSE))
#' @importFrom checkmate assertFlag
#' @export
dprAdaptiveGibbsWeights <- function(
    X,
    y,
    fitRetention = c("none", "slim", "full"),
    retainX = FALSE,
    methodArgs = DprOptions()
) {
    .assertMethodOptions(methodArgs, "DprOptions", "methodArgs")
    fitRetention <- arg_match(fitRetention)
    dprWeights(
        X,
        y,
        fittingMethod = "Adaptive_Gibbs",
        fitRetention = fitRetention,
        retainX = retainX,
        methodArgs = methodArgs
    )
}

#' Compute Weights Using Spike-and-Slab Variational Bayes
#'
#' Fits a linear model with a spike-and-slab prior by mean-field variational
#' Bayes, via \code{sparsevb::svb.fit} (Ray and Szabo 2022), and returns each
#' variant's posterior mean effect, \code{gamma * mu}: its inclusion
#' probability times its slab mean. The slab is Laplace unless
#' \code{methodArgs} sets \code{slab = "gaussian"}. sparsevb starts from
#' ridge and lasso fits, and by default updates the variants in decreasing
#' order of their ridge estimates.
#'
#' The noise SD is \code{methodArgs}' \code{noise_sd} when it sets one.
#' Otherwise it is estimated as sparsevb would, by
#' \code{selectiveInference::estimateSigma}, but on the distinct columns of
#' \code{X}: that estimator refuses identical columns, which variants in
#' perfect LD give it, and dropping exact copies of a column leaves its
#' lasso fit unchanged.
#'
#' Two behaviours of sparsevb 0.1.1 are left as they are. With
#' \code{intercept = TRUE}, the \code{intercept} it returns is the ridge
#' estimate it started from, not the fitted one; the weights do not use it.
#' With the Gaussian slab, the log-odds of inclusion carry an extra 0.5 that
#' the update it implements (Carbonetto and Stephens 2012) does not have.
#'
#' @param X A numeric matrix of predictors.
#' @param y A numeric response vector.
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights: \code{"none"} (default) attaches
#'   nothing, \code{"slim"} and \code{"full"} the fitted posterior
#'   (\code{mu}, \code{sigma}, \code{gamma} and \code{intercept}, over the
#'   variants fitted), which keeps no copy of the data, so the two retaining
#'   levels behave alike here. Either names the variants it was fitted on,
#'   in \code{variantIds}.
#' @param methodArgs Options forwarded to \code{sparsevb::svb.fit}, built
#'   with \code{\link{SparsevbOptions}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of length \code{ncol(X)} of variant weights;
#'   zero-variance columns receive weight 0. Unless
#'   \code{fitRetention = "none"} the fit is attached as attribute
#'   \code{"fit"}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' sparsevbWeights(X, y)
#' @export
sparsevbWeights <- function(
    X,
    y,
    fitRetention = c("none", "slim", "full"),
    methodArgs = SparsevbOptions()
) {
    .assertMethodOptions(methodArgs, "SparsevbOptions", "methodArgs")
    retention <- arg_match(fitRetention)
    if (!requireNamespace("sparsevb", quietly = TRUE)) {
        abort("Package 'sparsevb' is required for this function.")
    }
    keep <- .dropZeroVariance(X, "sparsevbWeights")
    XKeep <- X[, keep, drop = FALSE]
    user <- as.list(methodArgs)
    noise <- .sparsevbNoiseSd(XKeep, y, user)
    fit <- exec(sparsevb::svb.fit, X = XKeep, Y = y, !!!user, !!!noise)
    weights <- replace(rep(0, ncol(X)), keep, fit$mu * fit$gamma)
    if (identical(retention, "none")) {
        return(weights)
    }
    `attr<-`(weights, "fit", .withFitVariantIds(fit, colnames(XKeep)))
}

# svb.fit's `noise_sd` when the caller set none: its own estimator,
# selectiveInference::estimateSigma(), on the distinct columns. That
# estimator refuses identical columns, which variants in perfect LD give it,
# and dropping exact copies of a column leaves its lasso fit unchanged.
# @noRd
.sparsevbNoiseSd <- function(X, y, user) {
    if (is_in("noise_sd", names(user))) {
        return(list())
    }
    distinct <- X[, !duplicated(X, MARGIN = 2L), drop = FALSE]
    list(noise_sd = selectiveInference::estimateSigma(distinct, y)$sigmahat)
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
# caller did not supply them), the prior covariance matrices, and the mash
# prior's residual correlation, the prior taken on Y's conditions
# (.mrmashAlignPriorToContexts()).
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
    threads,
    gridSd = 1
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
    mashPrior <- .mrmashAlignPriorToContexts(
        dataDrivenPriorMatrices,
        colnames(Y),
        ncol(Y)
    )
    # The grid is in units of the residual standard deviations `gridSd`
    # when the prior is relative to them (.mrmashGridSd())
    list(
        bInitMethod = .mrmashResolveBInit(Y, bInitMethod),
        # Shared prior-covariance builder (also used by mrmashRssWeights).
        priorBuilt = exec(
            buildMrmashPriorMatrices,
            Bhat = t(t(sumstats$Bhat) / gridSd),
            Shat = t(t(sumstats$Shat) / gridSd),
            K = ncol(Y),
            dataDrivenPriorMatrices = mashPrior,
            !!!priorOpts
        ),
        residualCorrelation = mashPrior$residualCorrelation
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
#'   (tissue/cell-types). Missing values are allowed: each condition is
#'   fitted on the individuals it is observed in (mr.mash's
#'   \code{missing_y_method = "approximate"}), as when studies on different
#'   individuals are stacked.
#' @param dataDrivenPriorMatrices A list of data-driven covariance matrices
#'   (\code{U}), and optionally the residual correlation of the mash prior
#'   they come from (\code{residualCorrelation}, conditions x conditions),
#'   which \code{MrmashOptions(vInitMethod = "mash")} starts the residual
#'   covariance from. Matrices named by condition are taken on the columns
#'   of \code{Y} by name, so a prior learned on more conditions can be used
#'   (see \code{\link{buildMrmashPriorMatrices}}). Default is NULL.
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
#'   here. It also takes the wrapper's own \code{vEigenFloor}, the floor on
#'   the estimated residual covariance's eigenvalues, and \code{vInitMethod},
#'   where that covariance starts; see \code{\link{twasEngineOptions}}.
#'
#' @param V Optional residual covariance matrix (conditions x conditions), or
#'   \code{NULL} to estimate it, starting where \code{vInitMethod} says. A
#'   supplied \code{V} is its own start, so it cannot be given with
#'   \code{vInitMethod}.
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
#' \item{Y}{An n x r matrix of responses at the last iteration (only when
#' missing values in Y are imputed, \code{missing_y_method = "impute"}).}
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
        threads = threads,
        gridSd = .mrmashGridSd(methodArgs, V, Y)
    )
    time1 <- proc.time()
    bInit <- .mrmashInitB(X, Y, setup$bInitMethod, standardize, BPPARAM)
    vInit <- .mrmashInitV(
        Y,
        V,
        updateV,
        updateVMethod,
        as.list(methodArgs)$vEigenFloor %||% 1e-3,
        as.list(methodArgs)$vInitMethod,
        setup$residualCorrelation
    )
    fitMrmash <- .mrmashFit(
        setup$priorBuilt$S0,
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

# Robust residual-covariance init. Returns list(V, updateV): a supplied V as
# it is, else the start `vInitMethod` names (.mrmashVStart()), whose
# eigenvalues are floored (.mrmashFloorV()).
# @noRd
.mrmashInitV <- function(
    Y,
    V,
    updateV,
    updateVMethod,
    vEigenFloor,
    vInitMethod = NULL,
    residualCorrelation = NULL
) {
    if (!is.null(V)) {
        if (!is.null(vInitMethod)) {
            abort(glue(
                "mrmashWrapper: a supplied `V` is where the residual ",
                "covariance starts, so `vInitMethod` in MrmashOptions() ",
                "cannot also be given."
            ))
        }
        return(list(V = V, updateV = updateV))
    }
    start <- .mrmashVStart(Y, vInitMethod %||% "data", residualCorrelation)
    if (updateVMethod == "diagonal") {
        return(list(
            V = diag(diag(start), nrow = nrow(start)),
            updateV = updateV
        ))
    }
    .mrmashFloorV(start, updateV, vEigenFloor)
}

# A starting V whose smallest eigenvalue is under vEigenFloor times its
# spectral norm has its eigenvalues raised to that floor and its update
# disabled, since updating would take it back to singular. Near-duplicate
# contexts make such a V, and mr.mash's ELBO can turn NaN on one.
# @noRd
.mrmashFloorV <- function(V, updateV, vEigenFloor) {
    eigenvalues <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    floorValue <- vEigenFloor * max(eigenvalues)
    if (min(eigenvalues) >= floorValue) {
        return(list(V = V, updateV = updateV))
    }
    warn(glue(
        "mrmashWrapper: the residual covariance of Y is close to singular ",
        "(smallest / largest eigenvalue ",
        "{signif(min(eigenvalues) / max(eigenvalues), 3)}), as ",
        "near-duplicate contexts make it. Raising its eigenvalues to ",
        "vEigenFloor = {vEigenFloor} of the largest and not updating it ",
        "during the fit; `vEigenFloor` in MrmashOptions() sets the floor."
    ))
    list(
        V = V + diag(floorValue - min(eigenvalues), nrow(V)),
        updateV = FALSE
    )
}

# Where the residual covariance starts (MrmashOptions(vInitMethod)): the
# covariance of Y ("data", .mrmashComputeVInit()), or the mash prior's
# residual correlation scaled by the standard deviations of Y's columns
# ("mash"), zero for a pair of responses never observed together and only
# its diagonal when that is not positive definite, as the covariance of Y is.
# @noRd
#' @importFrom checkmate assertMatrix
#' @importFrom stats sd
.mrmashVStart <- function(Y, vInitMethod, residualCorrelation) {
    if (identical(vInitMethod, "data")) {
        return(.mrmashComputeVInit(Y))
    }
    if (is.null(residualCorrelation)) {
        abort(glue(
            "mrmashWrapper: vInitMethod = \"mash\" starts the residual ",
            "covariance at the mash prior's residual correlation, but ",
            "`dataDrivenPriorMatrices` carries no `residualCorrelation` (a ",
            "MashPrior from mashPipeline() gives one)."
        ))
    }
    assertMatrix(
        residualCorrelation,
        nrows = ncol(Y),
        ncols = ncol(Y),
        .var.name = "dataDrivenPriorMatrices$residualCorrelation"
    )
    sds <- unname(apply(Y, 2L, sd, na.rm = TRUE))
    together <- crossprod(!is.na(Y)) > 0L
    scaled <- replace(
        unname(residualCorrelation) * tcrossprod(sds),
        !together,
        0
    )
    .mrmashPdOrDiagonal(`dimnames<-`(scaled, list(colnames(Y), colnames(Y))))
}

# The residual covariance mr.mash starts from, its coefficients starting at
# zero: the covariance of Y. With missing values, each pair of responses
# over the rows where both are observed and zero for a pair never observed
# together, as mr.mash's approximate method for missing Y starts; only the
# diagonal when that is not positive definite.
# @noRd
#' @importFrom stats cov
.mrmashComputeVInit <- function(Y) {
    if (!anyNA(Y)) {
        return(cov(Y))
    }
    pairwise <- cov(Y, use = "pairwise.complete.obs")
    .mrmashPdOrDiagonal(replace(pairwise, is.na(pairwise), 0))
}

# A starting V as it is when positive definite, else its diagonal.
# @noRd
.mrmashPdOrDiagonal <- function(V) {
    if (min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) > 0) {
        return(V)
    }
    `dimnames<-`(diag(diag(V), nrow = nrow(V)), dimnames(V))
}

# Run mr.mash with the resolved prior / init / V. Everything pecotmr derives --
# the prior, the residual covariance, the coefficient initialisation -- is
# supplied here and cannot be overridden; the rest of mr.mash's arguments come
# from `methodArgs` under mr.mash's own names, with pecotmr's defaults
# underneath.
# @noRd
#' @importFrom purrr discard_at
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
    # vEigenFloor and vInitMethod are the wrapper's own, read by
    # .mrmashInitV().
    user <- discard_at(as.list(methodArgs), c("vEigenFloor", "vInitMethod"))
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
    mrMashArgs <- c(derived, list_modify(.mrmashFitDefaults(user), !!!user))
    .engineRequire("mr.mashr::mr.mash", names(mrMashArgs), "mrmashWrapper()")
    exec(mr.mashr::mr.mash, !!!mrMashArgs)
}

# pecotmr's defaults for the mr.mash options the caller (`user`) leaves:
# the prior relative to the residual standard deviations where mr.mashr can
# (.mrmashRelativeS0(); its grid then built in their units,
# .mrmashGridSd()); DAAREM cannot drop components.
# @noRd
.mrmashFitDefaults <- function(user) {
    c(
        list(
            update_w0 = TRUE,
            tol = 0.01,
            max_iter = 5000,
            convergence_criterion = "ELBO",
            compute_ELBO = TRUE,
            verbose = FALSE,
            w0_threshold = if (identical(user$acceleration, "daarem")) {
                0
            } else {
                1e-8
            }
        ),
        compact(list(relative_S0 = if (.mrmashRelativeS0()) TRUE))
    )
}

# Whether mr.mash and mr.mash.rss take a prior relative to the residual
# standard deviations (their relative_S0, which mr.mashr 0.3.44 lacks): the
# wrappers' default where they do; otherwise the prior grid is on the
# outcomes' scale, as it was before relative_S0.
# @noRd
.mrmashRelativeS0 <- function() {
    .engineTakes("mr.mashr::mr.mash", "relative_S0") &&
        .engineTakes("mr.mashr::mr.mash.rss", "relative_S0")
}

# Whether the installed mr.mashr fits summary statistics across studies:
# an LD reference per study (sample_groups), and mr.mash.rss estimating the
# residual and Y covariances it is not given (V, covY) or giving residual
# variances covY's correlations -- which the mr.mashr that has them marks
# by mr.mash.rss's sample_groups, the covariances' handling showing in no
# formal.
# @noRd
.mrmashFitsStudies <- function() {
    .engineTakes("mr.mashr::mr.mash.rss", "sample_groups")
}

# The residual standard deviations mr.mash's prior grid is built in units
# of when the prior is relative to them (mr.mash's relative_S0, the
# wrappers' default where mr.mashr takes it, .mrmashRelativeS0()): those of
# the residual covariance the fit starts from (.residualSd()); 1 when the
# prior is not relative.
# @noRd
.mrmashGridSd <- function(methodArgs, V, Y = NULL) {
    if (!isTRUE(as.list(methodArgs)$relative_S0 %||% .mrmashRelativeS0())) {
        return(1)
    }
    .residualSd(V, Y)
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
#' coefs <- computeCoefficientsGlasso(X = X, Y = Y, standardize = TRUE,
#'   Xnew = NULL)
#' head(coefs$Bhat)
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
#'   returns such a list directly. Matrices whose rows and columns are named
#'   by condition are taken on the conditions of \code{Bhat}'s columns, by
#'   name, so a prior learned on more conditions can be used; each must hold
#'   all of them. Unnamed matrices, and any when \code{Bhat}'s columns are
#'   unnamed, must be K x K, in \code{Bhat}'s column order.
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
#' priors <- buildMrmashPriorMatrices(Bhat = Bhat, Shat = Shat)
#' length(priors$S0)
#' head(names(priors$S0))
#' priors$priorGrid
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
        .mrmashAlignPriorToContexts(dataDrivenPriorMatrices, colnames(Bhat), K),
        canonicalPriorMatrices,
        canonicalCovs
    )
    user <- as.list(expandCovs)
    S0 <- exec(
        mr.mashr::expand_covs,
        S0_raw,
        priorGrid,
        !!!user
    )
    list(S0 = S0, priorGrid = priorGrid)
}

# A data-driven prior (`prior`, as dataDrivenPriorMatrices: its covariances
# `U`, and its residual correlation when it carries one) on the data's
# conditions `contexts` (NULL when they are unnamed), of which there are K.
# A matrix named by condition is taken on them by name (.mashSubsetMatrix()),
# so a prior learned on more conditions, such as a mash fit across cohorts,
# serves a trait observed in only some of them; an unnamed one must already
# be K x K, in the data's order. A prior with neither, such as a bare list
# of matrices, is refused: nothing would read it.
# @noRd
.mrmashAlignPriorToContexts <- function(prior, contexts, K) {
    if (is.null(prior)) {
        return(NULL)
    }
    U <- prior$U
    residual <- prior$residualCorrelation
    if (is.null(U) && is.null(residual)) {
        abort(glue(
            "dataDrivenPriorMatrices carries neither covariance matrices ",
            "(`U`) nor a residual correlation (`residualCorrelation`); a ",
            "list of covariance matrices goes under `U`: list(U = matrices)."
        ))
    }
    list_assign(
        prior,
        !!!compact(list(
            U = if (!is.null(U)) {
                map(U, .mrmashAlignPriorMatrix, contexts = contexts, K = K)
            },
            residualCorrelation = if (!is.null(residual)) {
                .mrmashAlignPriorMatrix(residual, contexts, K)
            }
        ))
    )
}

# One matrix of a data-driven prior on the data's conditions
# (.mrmashAlignPriorToContexts()): by name when both name them, else as it
# is, which must then be K x K.
# @noRd
.mrmashAlignPriorMatrix <- function(m, contexts, K) {
    mat <- as.matrix(m)
    named <- !is.null(rownames(mat)) && !is.null(colnames(mat))
    if (named && !is.null(contexts)) {
        return(.mashSubsetMatrix(mat, contexts))
    }
    if (nrow(mat) != K || ncol(mat) != K) {
        fix <- if (named) {
            "name the data's conditions (the columns of Y, Bhat or stat$z)"
        } else {
            "name its rows and columns by condition"
        }
        abort(glue(
            "A data-driven prior matrix is {nrow(mat)} x {ncol(mat)} but the ",
            "data hold {K} condition(s): {fix} to match them by name."
        ))
    }
    mat
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
    user <- as.list(canonicalCovs)
    canonical <- exec(
        mr.mashr::compute_canonical_covs,
        K,
        !!!user
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
    lambda1,
    lambda2,
    thr,
    maxiter,
    maxSwaps,
    strategy
) {
    if (.rssIsRidge(penalty, alpha)) {
        return(.rssRidgeSolve(z, lam, ld))
    }
    settings <- list(
        pathR = lam,
        penaltyStr = penalty,
        strategyStr = strategy,
        gamma = gamma,
        alpha = alpha,
        lambda1 = lambda1,
        lambda2 = lambda2,
        thr = thr,
        maxiter = as.integer(maxiter),
        maxSwaps = as.integer(maxSwaps),
        control = .rssPathControl()
    )
    if (!is.null(ld$ldFactor)) {
        return(exec(
            penalizedRssFactorCpp,
            zR = z,
            factorR = ld$ldFactor,
            shrink = ld$shrink %||% 0,
            !!!settings
        ))
    }
    exec(
        penalizedRssCpp,
        zR = z,
        ldR = ld$R,
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

# The minimum of one CV-mean path vector (coerced numeric).
# @noRd
.rssMinNumeric <- function(v) {
    min(as.numeric(v))
}
