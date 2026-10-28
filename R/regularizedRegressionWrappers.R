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
        !!!methodArgs
    )
    model <- exec(mr.ash.rss, !!!callArgs)
    if (identical(fitRetention, "none")) {
        return(model$mu1)
    }
    `attr<-`(model$mu1, "fit", model)
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
#'
#' @return A list containing the posterior estimates: - betaEst: Posterior
#'   estimates of SNP effect sizes. - psiEst: Posterior estimates of psi
#'   (shrinkage parameters). - sigmaEst: Posterior estimate of the residual
#'   variance. - phiEst: Posterior estimate of the global shrinkage parameter.
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
    seed = NULL
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
        seed
    )
}

# PRS-CS on an LD source (.ldSourceRows()): the sampler that factors R for a
# matrix, or for a factor the one that never forms R (prscs_mcmc.h), which
# draws differently from the same posterior.
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
    seed
) {
    if (!is.null(maf) && length(bhat) != length(maf)) {
        abort("The length of 'bhat' must be the same as 'maf'.")
    }
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
        seed = seed
    )
    result <- if (is.null(ld$F)) {
        exec(prsCsCpp, !!!settings, ldBlk = list(blk1 = ld$R))
    } else {
        exec(prsCsFactorCpp, !!!settings, factorR = ld$F)
    }
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
                "seed"
            )
        ),
        !!!as.list(methodArgs)
    )
    model <- exec(.prsCsRun, bhat = bhat, ld = ld, n = n, !!!settings)
    if (identical(fitRetention, "none")) {
        return(model$betaEst)
    }
    `attr<-`(model$betaEst, "fit", model)
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
#' @return A list with \code{betaEst} (the estimated effect sizes) and
#'   \code{h2} (the estimated heritability).
#' @note This function wraps a rewritten and adapted version of the SDPR C++
#'   implementation. SDPR is described in Zhou G, Zhao H (2021), "A fast and
#'   robust Bayesian nonparametric method for prediction of complex traits
#'   using summary statistics", PLoS Genetics 17(7): e1009697.
#'   \doi{10.1371/journal.pgen.1009697}
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
    seed = NULL
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
        seed
    )
}

# SDPR on an LD source (.ldSourceRows()): the sampler on the matrix, or the
# same sampler on the low-rank form of SDPR's matrices from a factor of it
# (sdpr_mcmc.h), which draws as the matrix route does without forming R.
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
    seed
) {
    .sdprValidate(M, perVariantSampleSize, array)
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
        # The kernel's thread pool splits its work across LD blocks, and a
        # window is always one block, so a second thread would sit idle.
        numThreads = 1L,
        optLlk = as.integer(optLlk),
        verbose = verbose,
        seed = seed
    )
    if (is.null(ld$F)) {
        # The C++ backend takes a block list; the window is one block.
        return(exec(sdprCpp, !!!settings, LD = list(blk1 = ld$R)))
    }
    exec(sdprFactorCpp, !!!settings, factorR = ld$F)
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
                "seed"
            )
        ),
        !!!as.list(methodArgs)
    )
    model <- exec(.sdprRun, bhat = bhat, ld = ld, n = n, !!!settings)
    if (identical(fitRetention, "none")) {
        return(model$betaEst)
    }
    `attr<-`(model$betaEst, "fit", model)
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
        # mrmashWrapper's `prior` group, so they are routed there rather
        # than passed as top-level arguments it no longer has.
        split <- .mrmashSplitPriorArgs(as.list(methodArgs))
        mrmashFit <- exec(
            mrmashWrapper,
            X,
            Y,
            dataDrivenPriorMatrices = dataDrivenPriorMatrices,
            prior = exec(MrmashPriorParam, !!!split$prior),
            !!!split$rest
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
        fitRetention
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
    callArgs <- c(
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
        # `dots` is the caller's MrmashOptions() record; c() would append it as
        # one opaque element instead of splicing its entries.
        as.list(dots)
    )
    exec(mr.mashr::mr.mash.rss, !!!callArgs)
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
#'   (\code{LassosumOptions}); \code{enetWeights} / \code{ridgeWeights} fit
#'   with glmnet too, while \code{enetRssWeights} / \code{ridgeRssWeights}
#'   use \code{PenalizedRssOptions}; \code{scadWeights} / \code{mcpWeights} /
#'   \code{l0learnWeights} use ncvreg and L0Learn while their \code{*Rss}
#'   counterparts use \code{PenalizedRssOptions}; \code{dprGibbsWeights} uses
#'   RcppDPR and \code{sdprWeights} uses \code{SdprOptions}. Splitting them is
#'   what lets each be checked exactly, rather than against a union in which
#'   an option meant for the other path would pass and then be dropped.
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
            # genotype matrix is `W` here and `X` is the covariates (the
            # wrapper's `Z`). Both are pecotmr's to supply, so both are
            # refused -- matching only the names that happen to agree on
            # both sides is how `W` was missed.
            W = "supplied from the data by the pipeline",
            X = "supplied from the data by the pipeline",
            h2 = "the wrapper's `h2`",
            method = "the wrapper's `method`",
            nit = "the wrapper's `nit`",
            nburn = "the wrapper's `nburn`",
            nthin = "the wrapper's `nthin`"
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
            burnIn = "the wrapper's `burnIn`",
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
            update_V_method = "the wrapper's `updateVMethod`",
            nthreads = "the wrapper's `BPPARAM`"
        ),
        "MrmashOptions"
    )
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
        !!!methodArgs
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
#' @param methodArgs Additional arguments forwarded to
#'   \code{glmnet::cv.glmnet}, built with \code{\link{GlmnetOptions}}. A bare
#'   list is refused, since it cannot be checked against the engine's formals.
#' @return Numeric vector of length \code{ncol(X)} of per-variant weights.
#' @importFrom stats coef
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' glmnetWeights(X, y, alpha = 0.5)
#' @export
glmnetWeights <- function(X, y, alpha, methodArgs = GlmnetOptions()) {
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
    replace(eff.wgt, keep, coef(enet, s = "lambda.min")[2:(sum(keep) + 1)])
}

#' Compute TWAS weights via elastic net (glmnet)
#'
#' Convenience wrapper for \code{\link{glmnetWeights}} at a mixing between the
#' lasso and ridge, \code{alpha = 0.5} unless set.
#' @param alpha Elastic-net mixing in \code{[0, 1]}: the L1 share of the
#'   penalty. Default 0.5.
#' @inheritParams glmnetWeights
#' @return Numeric vector of length \code{ncol(X)} of per-variant weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' enetWeights(X, y)
#' # A larger ridge share:
#' enetWeights(X, y, alpha = 0.2)
#' @importFrom checkmate assertNumber
#' @export
enetWeights <- function(X, y, alpha = 0.5, methodArgs = GlmnetOptions()) {
    .assertMethodOptions(methodArgs, "GlmnetOptions", "methodArgs")
    assertNumber(alpha, lower = 0, upper = 1)
    glmnetWeights(X, y, alpha, methodArgs = methodArgs)
}

#' Compute TWAS weights via ridge regression (glmnet, alpha = 0)
#'
#' Convenience wrapper for \code{\link{glmnetWeights}} with \code{alpha = 0}.
#' @inheritParams glmnetWeights
#' @return Numeric vector of length \code{ncol(X)} of per-variant weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ridgeWeights(X, y)
#' @export
ridgeWeights <- function(X, y, methodArgs = GlmnetOptions()) {
    .assertMethodOptions(methodArgs, "GlmnetOptions", "methodArgs")
    glmnetWeights(X, y, 0, methodArgs = methodArgs)
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
lassoWeights <- function(X, y, methodArgs = GlmnetOptions()) {
    .assertMethodOptions(methodArgs, "GlmnetOptions", "methodArgs")
    glmnetWeights(X, y, 1, methodArgs = methodArgs)
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
#'   nothing, \code{"slim"} and \code{"full"} attach the fitted model. This
#'   engine keeps no trimmable intermediates, so the two retaining levels
#'   behave alike here.
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
#'   with \code{\link{QggOptions}}. A bare list is refused: it cannot be
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
    methodArgs = QggOptions()
) {
    .assertMethodOptions(methodArgs, "QggOptions", "methodArgs")
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
#' @param h2 Numeric or \code{NULL}. Prior heritability; \code{NULL}
#'   leaves the engine's own default.
#' @param nit Integer. Total MCMC iterations. Default \code{5000}.
#' @param nburn Integer. Burn-in iterations. Default \code{1000}.
#' @param nthin Integer. Thinning interval for retained MCMC samples.
#'   Default \code{5}.
#' @param methodArgs Options forwarded to \code{qgg::gbayes}, built
#'   with \code{\link{QggOptions}}. A bare list is refused: it cannot be
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
bayesNWeights <- function(
    X,
    y,
    Z = NULL,
    h2 = NULL,
    nit = 5000,
    nburn = 1000,
    nthin = 5,
    methodArgs = QggOptions()
) {
    .assertMethodOptions(methodArgs, "QggOptions", "methodArgs")
    bayesAlphabetWeights(
        X,
        y,
        method = "bayesN",
        Z,
        h2 = h2,
        nit = nit,
        nburn = nburn,
        nthin = nthin,
        methodArgs = methodArgs
    )
}
#' @title BayesL TWAS weights (Laplace prior, LASSO-equivalent)
#' @description Use laplace/double exponential distribution as prior. This is
#'   equivalent to Bayesian LASSO.
#' @param X Numeric genotype / design matrix (samples x variants).
#' @param y Numeric response (phenotype) vector of length \code{nrow(X)}.
#' @param Z Optional numeric matrix of fixed-effect covariates, or \code{NULL}.
#' @param h2 Numeric or \code{NULL}. Prior heritability; \code{NULL}
#'   leaves the engine's own default.
#' @param nit Integer. Total MCMC iterations. Default \code{5000}.
#' @param nburn Integer. Burn-in iterations. Default \code{1000}.
#' @param nthin Integer. Thinning interval for retained MCMC samples.
#'   Default \code{5}.
#' @param methodArgs Options forwarded to \code{qgg::gbayes}, built
#'   with \code{\link{QggOptions}}. A bare list is refused: it cannot be
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
bayesLWeights <- function(
    X,
    y,
    Z = NULL,
    h2 = NULL,
    nit = 5000,
    nburn = 1000,
    nthin = 5,
    methodArgs = QggOptions()
) {
    .assertMethodOptions(methodArgs, "QggOptions", "methodArgs")
    bayesAlphabetWeights(
        X,
        y,
        method = "bayesL",
        Z,
        h2 = h2,
        nit = nit,
        nburn = nburn,
        nthin = nthin,
        methodArgs = methodArgs
    )
}
#' @title BayesA TWAS weights (t-distribution prior)
#' @description Use t-distribution as prior.
#' @param X Numeric genotype / design matrix (samples x variants).
#' @param y Numeric response (phenotype) vector of length \code{nrow(X)}.
#' @param Z Optional numeric matrix of fixed-effect covariates, or \code{NULL}.
#' @param h2 Numeric or \code{NULL}. Prior heritability; \code{NULL}
#'   leaves the engine's own default.
#' @param nit Integer. Total MCMC iterations. Default \code{5000}.
#' @param nburn Integer. Burn-in iterations. Default \code{1000}.
#' @param nthin Integer. Thinning interval for retained MCMC samples.
#'   Default \code{5}.
#' @param methodArgs Options forwarded to \code{qgg::gbayes}, built
#'   with \code{\link{QggOptions}}. A bare list is refused: it cannot be
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
bayesAWeights <- function(
    X,
    y,
    Z = NULL,
    h2 = NULL,
    nit = 5000,
    nburn = 1000,
    nthin = 5,
    methodArgs = QggOptions()
) {
    .assertMethodOptions(methodArgs, "QggOptions", "methodArgs")
    bayesAlphabetWeights(
        X,
        y,
        method = "bayesA",
        Z,
        h2 = h2,
        nit = nit,
        nburn = nburn,
        nthin = nthin,
        methodArgs = methodArgs
    )
}
#' @title BayesC TWAS weights (rounded-spike prior)
#' @description Use a rounded spike prior (low-variance Gaussian).
#' @param X Numeric genotype / design matrix (samples x variants).
#' @param y Numeric response (phenotype) vector of length \code{nrow(X)}.
#' @param Z Optional numeric matrix of fixed-effect covariates, or \code{NULL}.
#' @param pi Numeric in (0, 1). Prior proportion of non-null effects for the
#'   BayesC mixture. Default \code{0.1}.
#' @param h2 Numeric or \code{NULL}. Prior heritability; \code{NULL}
#'   leaves the engine's own default.
#' @param nit Integer. Total MCMC iterations. Default \code{5000}.
#' @param nburn Integer. Burn-in iterations. Default \code{1000}.
#' @param nthin Integer. Thinning interval for retained MCMC samples.
#'   Default \code{5}.
#' @param methodArgs Options forwarded to \code{qgg::gbayes}, built
#'   with \code{\link{QggOptions}}. A bare list is refused: it cannot be
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
bayesCWeights <- function(
    X,
    y,
    Z = NULL,
    pi = 0.1,
    h2 = NULL,
    nit = 5000,
    nburn = 1000,
    nthin = 5,
    methodArgs = QggOptions()
) {
    .assertMethodOptions(methodArgs, "QggOptions", "methodArgs")
    # `pi` is a qgg option, not a formal of bayesAlphabetWeights, so it joins
    # the option list rather than the argument list.
    bayesAlphabetWeights(
        X,
        y,
        method = "bayesC",
        Z,
        h2 = h2,
        nit = nit,
        nburn = nburn,
        nthin = nthin,
        methodArgs = .methodOptionsWith(methodArgs, "QggOptions", pi = pi)
    )
}
#' @title BayesR TWAS weights (hierarchical mixture prior)
#' @description Use a hierarchical Bayesian mixture model with four Gaussian
#'   components. Variances are scaled by 0, 0.0001, 0.001, and 0.01.
#' @param X Numeric genotype / design matrix (samples x variants).
#' @param y Numeric response (phenotype) vector of length \code{nrow(X)}.
#' @param Z Optional numeric matrix of fixed-effect covariates, or \code{NULL}.
#' @param h2 Numeric or \code{NULL}. Prior heritability; \code{NULL}
#'   leaves the engine's own default.
#' @param nit Integer. Total MCMC iterations. Default \code{5000}.
#' @param nburn Integer. Burn-in iterations. Default \code{1000}.
#' @param nthin Integer. Thinning interval for retained MCMC samples.
#'   Default \code{5}.
#' @param methodArgs Options forwarded to \code{qgg::gbayes}, built
#'   with \code{\link{QggOptions}}. A bare list is refused: it cannot be
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
bayesRWeights <- function(
    X,
    y,
    Z = NULL,
    h2 = NULL,
    nit = 5000,
    nburn = 1000,
    nthin = 5,
    methodArgs = QggOptions()
) {
    .assertMethodOptions(methodArgs, "QggOptions", "methodArgs")
    bayesAlphabetWeights(
        X,
        y,
        method = "bayesR",
        Z,
        h2 = h2,
        nit = nit,
        nburn = nburn,
        nthin = nthin,
        methodArgs = methodArgs
    )
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

.lassosumSelectLdQuadratic <- function(candidateBeta, corInput, ld) {
    bxy <- as.numeric(crossprod(corInput, candidateBeta))
    bxxb <- .ldQuadForms(ld, candidateBeta)
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
# LdMatrix holds. Matching, allele flips, standardizing and padding are the
# LD reference layer's own (.ldFactorFromRef(), .ldFromRef()).
# @noRd
.rssLdSource <- function(stat, ldReference, label) {
    ids <- .rssStatVariantIds(stat, label)
    if (.ldIsGenotypeRef(ldReference)) {
        return(list(F = .ldFactorFromRef(ldReference, ids, label)))
    }
    list(R = .ldFromRef(ldReference, ids, label))
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
# over `s`, the candidate accumulation, and the ldQuadratic / minFbeta
# selection. `config` is the run: `method` (which selection record the
# weights carry; "l0learn" also sweeps `lambda0`), the penalty `spec`
# (.rssPenalty()) and the caller's path `options`; `fitRetention` is how much
# of the search is kept.
# @noRd
.rssShrinkGridWeights <- function(
    stat,
    ld,
    s,
    config,
    selection = c("ldQuadratic", "minFbeta"),
    fitRetention = c("none", "slim", "full")
) {
    selection <- arg_match(selection)
    fitRetention <- arg_match(fitRetention)
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
    selectorResult <- if (selection == "ldQuadratic") {
        .lassosumSelectLdQuadratic(candidateBeta, corInput, ld)
    } else {
        .lassosumSelectMinFbeta(candidateBeta, candidateMeta)
    }
    bestBeta <- as.numeric(selectorResult$beta)
    weights <- .rssFinalize(bestBeta, selectorResult, candidateMeta, config)
    .rssAttachGridFit(
        weights,
        candidateBeta,
        candidateMeta,
        selectorResult$index,
        fitRetention
    )
}

# The search behind the grid weights as their "fit" attribute: `path`, one
# row per candidate fit with the solver's diagnostics and `selected` marking
# the one returned, and at fitRetention = "full" also `beta`, the candidates'
# coefficients (variants x candidates, in path's row order). Returns `weights`
# unchanged at fitRetention = "none".
# @noRd
.rssAttachGridFit <- function(
    weights,
    candidateBeta,
    candidateMeta,
    index,
    fitRetention
) {
    if (identical(fitRetention, "none")) {
        return(weights)
    }
    fit <- c(
        list(path = mutate(candidateMeta, selected = row_number() == index)),
        compact(list(beta = if (fitRetention == "full") candidateBeta))
    )
    `attr<-`(weights, "fit", fit)
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
#' @param s Numeric vector of shrinkage parameters to search over. Default:
#'   \code{c(0.2, 0.5, 0.9, 1.0)} following Mak et al (2017) and OTTERS.
#' @param selection Selection strategy. Default \code{"ldQuadratic"} uses
#'   \eqn{c^T \beta / \sqrt{\beta^T R \beta}} on the supplied LD matrix.
#'   \code{"minFbeta"} is retained as an explicit alternative for debugging.
#' @param fitRetention How much of the search is kept as the \code{"fit"}
#'   attribute of the returned weights. \code{"none"} (default) attaches
#'   nothing. \code{"slim"} keeps \code{path}, a table with one row per
#'   candidate (\code{s}, \code{lambda}) giving the solver's objective
#'   \code{fbeta}, \code{loss}, convergence flag \code{conv} and number of
#'   non-zero coefficients \code{nparams}, with \code{selected} marking the
#'   candidate returned. \code{"full"} also keeps \code{beta}, the
#'   candidates' coefficients (variants x candidates, in the rows' order).
#' @param methodArgs Options forwarded to the solver,
#'   \code{\link{lassosumRss}}, built with \code{\link{LassosumOptions}}. A
#'   bare list is refused: it cannot be checked against the engine, so a
#'   misspelled option would be silently ignored.
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
#' w <- lassosumRssWeights(stat, LdMatrix(cor(X)), fitRetention = "slim")
#' head(attr(w, "fit")$path)
#' @export
lassosumRssWeights <- function(
    stat,
    ldReference,
    s = c(0.2, 0.5, 0.9, 1.0),
    selection = c("ldQuadratic", "minFbeta"),
    fitRetention = c("none", "slim", "full"),
    methodArgs = LassosumOptions()
) {
    .assertMethodOptions(methodArgs, "LassosumOptions", "methodArgs")
    selection <- arg_match(selection)
    fitRetention <- arg_match(fitRetention)
    .rssShrinkGridWeights(
        stat,
        .rssLdSource(stat, ldReference, "lassosumRssWeights"),
        s,
        list(
            method = "lassosum",
            spec = .rssPenalty("lasso"),
            options = as.list(methodArgs)
        ),
        selection,
        fitRetention
    )
}

#' Penalized Regression on RSS (Summary Statistics) Objective
#'
#' The RSS penalized-regression solver: LASSO with its elastic-net and ridge
#' mixtures, MCP, SCAD, L0, L0L1 and L0L2. Uses coordinate descent on the
#' objective
#' \deqn{\beta^T R \beta - 2 \beta^T z + \mathrm{penalty}(\beta)}
#' where \eqn{R} is a (possibly pre-shrunk) LD matrix and \eqn{z = \hat\beta /
#' \sqrt{n}}. \code{lassosumRss()} is its LASSO at \code{alpha = 1}.
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
#' sI}) and selects the best candidate via LD-quadratic pseudovalidation or
#' minimum penalized objective.
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
    selection = c("ldQuadratic", "minFbeta"),
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
            options = as.list(methodArgs)
        ),
        selection,
        fitRetention
    )
}

#' Compute SCAD-Penalized Weights from Summary Statistics
#'
#' Fits SCAD-penalized regression on the RSS objective, searching over a
#' shrinkage grid \code{s} and lambda path. Model selection uses LD-quadratic
#' pseudovalidation by default.
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
#' @param selection Selection strategy: \code{"ldQuadratic"} (default) or
#'   \code{"minFbeta"}.
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
#' scadRssWeights(stat, LdMatrix(cor(X)))
#' @importFrom checkmate assertList assertNumeric assertNumber
#' @export
scadRssWeights <- function(
    stat,
    ldReference,
    s = c(0.2, 0.5, 0.9, 1.0),
    gamma = 3.7,
    alpha = 1.0,
    selection = c("ldQuadratic", "minFbeta"),
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
#' shrinkage grid \code{s} and lambda path. Model selection uses LD-quadratic
#' pseudovalidation by default.
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
#' @param selection Selection strategy: \code{"ldQuadratic"} (default) or
#'   \code{"minFbeta"}.
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
#' mcpRssWeights(stat, LdMatrix(cor(X)))
#' @importFrom checkmate assertList assertNumeric assertNumber
#' @export
mcpRssWeights <- function(
    stat,
    ldReference,
    s = c(0.2, 0.5, 0.9, 1.0),
    gamma = 3.0,
    alpha = 1.0,
    selection = c("ldQuadratic", "minFbeta"),
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
#' \code{s} and a lambda path. Model selection uses LD-quadratic
#' pseudovalidation by default.
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
#' enetRssWeights(stat, LdMatrix(cor(X)))
#' # A larger ridge share:
#' enetRssWeights(stat, LdMatrix(cor(X)), alpha = 0.2)
#' @importFrom checkmate assertList assertNumeric
#' @export
enetRssWeights <- function(
    stat,
    ldReference,
    alpha = 0.5,
    s = c(0.2, 0.5, 0.9, 1.0),
    selection = c("ldQuadratic", "minFbeta"),
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
#' model selection uses LD-quadratic pseudovalidation by default.
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
#' ridgeRssWeights(stat, LdMatrix(cor(X)))
#' @importFrom checkmate assertList assertNumeric
#' @export
ridgeRssWeights <- function(
    stat,
    ldReference,
    s = c(0.2, 0.5, 0.9, 1.0),
    selection = c("ldQuadratic", "minFbeta"),
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
#' selection uses LD-quadratic pseudovalidation by default.
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
#' @param selection Selection strategy: \code{"ldQuadratic"} (default) or
#'   \code{"minFbeta"}.
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
#' l0learnRssWeights(stat, LdMatrix(cor(X)))
#' @export
l0learnRssWeights <- function(
    stat,
    ldReference,
    penalty = c("L0", "L0L1", "L0L2"),
    s = c(0.2, 0.5, 0.9, 1.0),
    lambda0 = exp(seq(log(0.001), log(1), length.out = 10)),
    lambda = NULL,
    lambda2 = 0,
    selection = c("ldQuadratic", "minFbeta"),
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
        .rssLdSource(stat, ldReference, "l0learnRssWeights"),
        s,
        list(
            method = "l0learn",
            spec = spec,
            lambda0 = lambda0,
            options = list_assign(as.list(methodArgs), lambda = lambda)
        ),
        selection,
        fitRetention
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
#'   with \code{\link{NcvregOptions}}. A bare list is refused: it cannot be
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
    methodArgs = NcvregOptions()
) {
    .assertMethodOptions(methodArgs, "NcvregOptions", "methodArgs")
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
#'   with \code{\link{NcvregOptions}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of length `ncol(X)` of variant weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' scadWeights(X, y)
#' @export
scadWeights <- function(X, y, nfolds = 5, methodArgs = NcvregOptions()) {
    .assertMethodOptions(methodArgs, "NcvregOptions", "methodArgs")
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
#'   with \code{\link{NcvregOptions}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @return A numeric vector of length `ncol(X)` of variant weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' mcpWeights(X, y)
#' @export
mcpWeights <- function(X, y, nfolds = 5, methodArgs = NcvregOptions()) {
    .assertMethodOptions(methodArgs, "NcvregOptions", "methodArgs")
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
#'   with \code{\link{L0learnOptions}}. A bare list is refused: it cannot be
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
    methodArgs = L0learnOptions()
) {
    .assertMethodOptions(methodArgs, "L0learnOptions", "methodArgs")
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
#'   with \code{\link{BglrOptions}}. A bare list is refused: it cannot be
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
    methodArgs = BglrOptions()
) {
    .assertMethodOptions(methodArgs, "BglrOptions", "methodArgs")
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
#'   with \code{\link{BglrOptions}}. A bare list is refused: it cannot be
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
    methodArgs = BglrOptions()
) {
    .assertMethodOptions(methodArgs, "BglrOptions", "methodArgs")
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
#'   with \code{\link{BglrOptions}}. A bare list is refused: it cannot be
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
    methodArgs = BglrOptions()
) {
    .assertMethodOptions(methodArgs, "BglrOptions", "methodArgs")
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
#'   with \code{\link{DprOptions}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be
#'   silently ignored.
#' @param fitRetention How much of the fit is kept as the \code{"fit"}
#'   attribute of the returned weights: \code{"none"} (default) attaches
#'   nothing, \code{"slim"} and \code{"full"} attach the fitted model. This
#'   engine keeps no trimmable intermediates, so the two retaining levels
#'   behave alike here.
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
#'   initialise the coefficients run on its workers, and mr.mash's own
#'   threads (the summary statistics and the fit) number its workers when
#'   they run on this machine (\code{SerialParam}, \code{MulticoreParam});
#'   for other backends mr.mash runs on one thread. Default
#'   \code{SerialParam()}.
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
        threads = threads,
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
    threads,
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
        nthreads = threads
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
# and ridge mixtures, MCP, SCAD, L0...), on the LD matrix or on a factor of it
# (never forming R).
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
