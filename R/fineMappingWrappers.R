#' Convert Log Bayes Factors to Single Effects PIP
#'
#' This function converts log Bayes factors (LBF) to alpha values, optionally
#' using prior weights. It handles numerical stability by adjusting with the
#' maximum LBF value.
#'
#' @param lbf Numeric vector of log Bayes factors.
#' @param priorWeights Optional numeric vector of prior weights for each element
#'   in lbf.
#' @return A named numeric vector of alpha values corresponding to the input
#'   LBF.
#' @examples
#' lbf <- c(-0.5, 1.2, 0.3)
#' alpha <- lbfToAlphaVector(lbf)
#' print(alpha)
#' @noRd
lbfToAlphaVector <- function(lbf, priorWeights = NULL) {
    if (length(lbf) == 0L) {
        return(set_names(numeric(0), names(lbf)))
    }
    if (is.null(priorWeights)) {
        priorWeights <- rep(1 / length(lbf), length(lbf))
    }
    maxlbf <- max(lbf)

    # A non-finite maximum means every variant's Bayes factor underflowed, so
    # nothing updates the prior and the posterior IS the normalized prior.
    #
    # This deliberately does NOT trigger on `maxlbf == 0`, which the previous
    # guard did: an all-zero LBF vector is BF = 1 everywhere and falls out of
    # the softmax below as the normalized prior on its own, while a vector
    # whose best variant merely happens to sit at exactly 0 is perfectly
    # informative. Both used to be flattened to alpha = 0 (PIP 0).
    if (!is.finite(maxlbf)) {
        return(set_names(priorWeights / sum(priorWeights), names(lbf)))
    }

    # w is proportional to BF, subtract max for numerical stability
    w <- exp(lbf - maxlbf)

    # Posterior prob for each SNP
    wWeighted <- w * priorWeights
    weightedSumW <- sum(wWeighted)
    alpha <- wWeighted / weightedSumW

    set_names(alpha, names(lbf))
}

#' @title Convert a log-Bayes-factor matrix to Single Effect PIPs
#' @description Applies the 'lbfToAlphaVector' function row-wise to a matrix of
#'   log Bayes factors to convert them to Single Effect PIP values.
#'
#' @param lbf Matrix of log Bayes factors.
#' @return A matrix of alpha values with the same dimensions as the input LBF
#'   matrix.
#' @examples
#' lbfMatrix <- matrix(c(-0.5, 1.2, 0.3, 0.7, -1.1, 0.4), nrow = 2)
#' alphaMatrix <- lbfToAlpha(lbfMatrix)
#' print(alphaMatrix)
#' @export
lbfToAlpha <- function(lbf) {
    alphaMatrix <- t(apply(as.matrix(lbf), 1, lbfToAlphaVector))
    if (ncol(lbf) != 1) {
        return(alphaMatrix)
    }
    # t() turns a single-column lbf into a row vector; restore the shape.
    matrix(alphaMatrix, ncol = 1, dimnames = list(NULL, colnames(lbf)))
}

formatPipColumn <- function(method) {
    str_c("pip_", method)
}

resolvePipColumn <- function(topLoci, method = NULL) {
    if (is.null(topLoci) || nrow(topLoci) == 0) {
        return(NULL)
    }
    if (!is.null(method)) {
        pipCol <- formatPipColumn(method)
        if (is_in(pipCol, names(topLoci))) return(pipCol)
    }
    if (is_in("pip", names(topLoci))) {
        return("pip")
    }
    pipCols <- names(topLoci)[str_detect(names(topLoci), "^pip_")]
    if (length(pipCols) == 1) {
        return(pipCols)
    }
    NULL
}

formatCsColumn <- function(coverage, method) {
    pct <- as.numeric(coverage) * 100
    if (is.na(pct)) {
        abort("coverage must be numeric.")
    }
    label <- if (abs(pct - round(pct)) < 1e-8) {
        as.character(as.integer(round(pct)))
    } else {
        str_replace_all(
            format(pct, scientific = FALSE, trim = TRUE),
            "\\.",
            "_"
        )
    }
    str_c("CS_", label, "_", method)
}

.translateLegacyCsColumnName <- function(coverage) {
    if (is.null(coverage)) {
        return(NULL)
    }
    map_chr(coverage, .translateOneLegacyCsColumn)
}

# The legacy per-method `pip_susie` column is the plain `pip` column unless
# the table already carries one.
# @noRd
.translateLegacyTopLociNames <- function(nms) {
    translated <- .translateLegacyCsColumnName(nms)
    if (is_in("pip", translated)) {
        return(translated)
    }
    if_else(translated == "pip_susie", "pip", translated)
}

.translateLegacyTopLociCsColumns <- function(topLoci) {
    if (!is.data.frame(topLoci)) {
        return(topLoci)
    }
    `names<-`(topLoci, .translateLegacyTopLociNames(names(topLoci)))
}

# Translate a camelCase pecotmr method identifier (e.g. "susieInfRss") into the
# snake_case form (e.g. "susie_inf_rss") used in the documented top_loci schema.
# Single-word identifiers (e.g. "susie", "mvsusie", "fsusie") pass through.
.camelToSnakeMethod <- function(method) {
    if (is.null(method) || length(method) == 0L) {
        return(method)
    }
    lookup <- c(
        susieInf = "susie_inf",
        susieAsh = "susie_ash",
        susieRss = "susie_rss",
        susieInfRss = "susie_inf_rss",
        susieAshRss = "susie_ash_rss",
        singleEffect = "single_effect",
        bayesianConditionalRegression = "bayesian_conditional_regression"
    )
    map_chr(method, .camelToSnakeOne, lookup = lookup)
}

.setFinemappingFitClass <- function(fit, method) {
    if (is.null(fit)) {
        return(NULL)
    }
    methodClass <- switch(
        method,
        susie = "susie",
        susieInf = "susieInf",
        susieRss = "susieRss",
        singleEffect = "susieRss",
        bayesianConditionalRegression = "susieRss",
        fsusie = "susiF",
        mvsusie = "mvsusie",
        NULL
    )
    if (is.null(methodClass)) {
        return(fit)
    }
    `class<-`(fit, unique(c(methodClass, class(fit))))
}

# Build the argument list for a SuSiE / SuSiE-ash fit initialised from a
# prior SuSiE-inf fit. `unmappableEffects` controls which branch the
# downstream fit takes: "none" yields the standard SuSiE-inf-initialised
# SuSiE; "ash" yields SuSiE-ash with the SuSiE-inf warm start.
prepareSusieFromInfArgs <- function(
    args,
    susieInfFit,
    refineDefault = NULL,
    unmappableEffects = c("none", "ash")
) {
    unmappableEffects <- arg_match(unmappableEffects)
    L <- args[["L"]] %||% length(susieInfFit$V)
    list_assign(
        args,
        unmappable_effects = unmappableEffects,
        model_init = susieInfFit,
        !!!compact(list(
            refine = if (is.null(args[["refine"]])) refineDefault,
            convergence_method = if (unmappableEffects == "ash") {
                args[["convergence_method"]] %||% "pip"
            },
            # Clamped, not passed through raw: a caller's L_greedy above the
            # number of inf effects would ask susie for effects it cannot warm-
            # start.
            L_greedy = if (!is.null(args[["L_greedy"]])) {
                min(length(susieInfFit$V), L)
            }
        ))
    )
}

# Merge one stage's overrides onto the shared defaults of a two-stage SuSiE
# chain. Either may arrive as a constructor result or as the empty list that
# means "no options", so each is taken down to a plain list first: list_modify
# works on lists, not on the S4 SimpleList a constructor returns.
# @noRd
.fmChainStageArgs <- function(shared, stage) {
    list_modify(as.list(shared), !!!compact(as.list(stage)))
}

# The three argument bundles of a two-stage chain, each checked against the
# engine it reaches. `susie` and `susieInf` forward to the same two susieR
# entry points, so all three accept the same names and SusieOptions() names them
# in the error whichever bundle was at fault.
#
# The bundles default to list() rather than to SusieOptions(): a formal cannot
# default to a call on its own name -- `SusieOptions = SusieOptions()` is a
# recursive default argument reference and errors when forced. An empty list
# already means "no options" to .assertMethodOptions, so nothing is lost but the
# self-documenting signature.
# @noRd
.fmAssertChainArgs <- function(args, SusieInfOptions, SusieOptions) {
    .assertMethodOptions(args, "SusieOptions", "args")
    .assertMethodOptions(SusieInfOptions, "SusieInfOptions", "SusieInfOptions")
    .assertMethodOptions(SusieOptions, "SusieOptions", "SusieOptions")
    invisible(NULL)
}

#' @noRd
fitSusieInfThenSusie <- function(
    X,
    y,
    args = list(),
    SusieInfOptions = list(),
    SusieOptions = list(),
    fittedModels = NULL
) {
    .fmAssertChainArgs(args, SusieInfOptions, SusieOptions)
    # Two-stage chain built from the shared per-token fitter (.fmFitSusieIndiv),
    # so the susieInf fit arguments and the susieInf -> susie initialisation
    # live in one place rather than being duplicated here and in the pipeline.
    cached <- fittedModels %||% list()
    susieInfFit <- if (is.null(cached[["susieInf"]])) {
        .fmFitSusieIndiv(
            X,
            y,
            "susieInf",
            userArgs = .fmChainStageArgs(args, SusieInfOptions)
        )
    } else {
        .setFinemappingFitClass(cached[["susieInf"]], "susieInf")
    }
    susieFit <- if (is.null(cached[["susie"]])) {
        .fmFitSusieIndiv(
            X,
            y,
            "susie",
            chainFromInf = susieInfFit,
            userArgs = .fmChainStageArgs(args, SusieOptions)
        )
    } else {
        .setFinemappingFitClass(cached[["susie"]], "susie")
    }
    list(susie = susieFit, susieInf = susieInfFit)
}

#' Two-stage SuSiE-RSS Fine-mapping
#'
#' RSS analog of \code{fitSusieInfThenSusie}. Fits SuSiE-inf via \code{susieRss}
#' first, then initialises standard SuSiE-RSS from the SuSiE-inf result. The
#' single pair of fits can be used both for fine-mapping post-processing and
#' TWAS weight extraction.
#'
#' @param z Numeric vector of z-scores.
#' @param R LD correlation matrix.
#' @param n Sample size (scalar).
#' @param args Defaults forwarded to both fits, built with
#'   \code{\link{SusieOptions}}. A bare list is refused, since it cannot be
#'   checked; \code{list()} (the default) means no options.
#' @param SusieInfOptions SuSiE-inf-specific overrides, built with
#'   \code{\link{SusieInfOptions}}; they take precedence over \code{args}.
#' @param SusieOptions Standard SuSiE-RSS-specific overrides, built with
#'   \code{\link{SusieOptions}}; they take precedence over \code{args}.
#' @param fittedModels Optional list with pre-fitted \code{$susie} and/or
#'   \code{$susieInf} objects to skip re-fitting.
#' @return A list with \code{susie} and \code{susieInf} fit objects.
#' @importFrom susieR susie_rss
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ss <- lapply(seq_len(ncol(X)), function(j) {
#'   coef(summary(lm(y ~ X[, j])))[2, 1:2]
#' })
#' stat <- list(
#'   bhat = vapply(ss, `[`, numeric(1), 1L),
#'   shat = vapply(ss, `[`, numeric(1), 2L),
#'   z = vapply(ss, function(s) s[1] / s[2], numeric(1)),
#'   n = rep(nrow(X), ncol(X)))
#' LD <- cor(X)
#' fitSusieInfThenSusieRss(z = stat$z, R = LD, n = nrow(X))
#' @importFrom checkmate assertNumeric assertNumber assertList
#' @export
fitSusieInfThenSusieRss <- function(
    z,
    R,
    n,
    args = list(),
    SusieInfOptions = list(),
    SusieOptions = list(),
    fittedModels = NULL
) {
    assertNumeric(z)
    assertNumber(n, lower = 0, finite = TRUE)
    .fmAssertChainArgs(args, SusieInfOptions, SusieOptions)
    assertList(fittedModels, null.ok = TRUE)
    # RSS analog of fitSusieInfThenSusie, built from the shared per-token RSS
    # fitter (.fmFitSusieRss). .fmFitSusieRss tags every fit "susieRss", so the
    # inf fit is re-tagged "susieInf" to preserve this wrapper's contract.
    cached <- fittedModels %||% list()
    infRaw <- cached[["susieInf"]] %||%
        .fmFitSusieRss(
            z,
            R,
            n,
            "susieInf",
            userArgs = .fmChainStageArgs(args, SusieInfOptions)
        )
    susieInfFit <- .setFinemappingFitClass(infRaw, "susieInf")
    susieRaw <- cached[["susie"]] %||%
        .fmFitSusieRss(
            z,
            R,
            n,
            "susie",
            chainFromInf = susieInfFit,
            userArgs = .fmChainStageArgs(args, SusieOptions)
        )
    susieFit <- .setFinemappingFitClass(susieRaw, "susieRss")
    list(susie = susieFit, susieInf = susieInfFit)
}

#' Post-process Fine-mapping Fits
#'
#' Applies method-aware post-processing to one or more SuSiE-family fits and
#' builds both a method-specific result list and shared top-loci tables.
#'
#' @param fits Named list of fine-mapping fits. Names define method identity,
#'   for example \code{susie}, \code{susieInf}, \code{susieRss}, \code{mvsusie},
#'   or \code{fsusie}.
#' @param dataX Genotype matrix, LD/correlation matrix, or other method-specific
#'   input used for credible-set purity and correlations.
#' @param dataY Phenotype vector/matrix or summary statistics. Default NULL.
#' @param xScalar Scaling factor for genotype effects. Default 1.
#' @param yScalar Scaling factor for phenotype effects. Default 1.
#' @param af Effect-allele frequencies (exported as the \code{af} column; never
#'   MAF). Default NULL.
#' @param n Optional per-variant sample size, exported as the \code{N} column.
#'   Default NULL -> \code{N} falls back to the fit's own scalar sample size.
#' @param credibleSetArgs How credible sets are built and reported, built with
#'   \code{\link{CredibleSetParam}}: \code{coverage},
#'   \code{secondaryCoverage}, \code{signalCutoff} (the PIP cutoff for
#'   including non-credible-set variants in top loci), \code{minAbsCorr} and
#'   \code{medianAbsCorr} for purity, and \code{includeAllCs}.
#' @param fitRetention How much of each fit is kept: \code{"slim"} (default)
#'   trims the retained fit to what downstream needs, \code{"full"} keeps the
#'   whole fit object. The \code{topLoci} table's per-credible-set columns
#'   are governed by \code{credibleSet}'s \code{perCsColumns} instead.
#' @param otherQuantities Optional named list of extra per-method quantities to
#'   carry on the result untouched (returned as the \code{otherQuantities}
#'   element). Default \code{NULL}.
#' @param region Optional \code{"chr:start-end"} string naming the region the
#'   fits cover; recorded on the result. Default \code{NULL}.
#' @param priorEffTol Numeric (length 1). Effects whose prior variance
#'   \code{V} is at or below this tolerance are dropped as unconverged before
#'   the credible sets are read. Default \code{1e-9}.
#' @param csInput One of \code{"X"}, \code{"Xcorr"}, \code{"fsusie"}: how
#'   credible-set purity is computed from \code{dataX}. \code{NULL} (default)
#'   lets each method pick its own -- \code{"Xcorr"} for the RSS methods,
#'   \code{"fsusie"} for fsusie, \code{"X"} otherwise.
#' @param conditionIdx Integer or \code{NULL}. For a multi-condition
#'   (3-D) fit, the condition to slice out; \code{NULL} (default) keeps the
#'   unconditioned fit.
#' @return A list with \code{finemappingResults} (per-method post-processed
#'   objects, each carrying a trimmed fit and method-specific intermediates) and
#'   a single unified \code{top_loci} table in the fixed 22-column shape (see
#'   the internal \code{buildTopLoci}). Per-method contributions are row-bound
#'   into \code{top_loci} by an outer method for-loop.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:40]
#' y <- eqtlRegionExample$yRes
#' fit <- susieR::susie(X, y, L = 5)
#' postprocessFinemappingFits(fits = list(susie = fit), dataX = X, dataY = y)
#' @export
postprocessFinemappingFits <- function(
    fits,
    dataX,
    dataY = NULL,
    xScalar = 1,
    yScalar = 1,
    af = NULL,
    n = NULL,
    credibleSetArgs = CredibleSetParam(),
    fitRetention = "slim",
    otherQuantities = NULL,
    region = NULL,
    priorEffTol = 1e-9,
    csInput = NULL,
    conditionIdx = NULL
) {
    fits <- fits[!map_lgl(fits, is.null)]
    if (length(fits) == 0) {
        abort("At least one fine-mapping fit must be supplied.")
    }
    if (is.null(names(fits)) || any(names(fits) == "")) {
        abort("fits must be a named list; names define method identity.")
    }
    .ppFitsCombine(.ppFitsPerMethod(
        fits,
        dataX = dataX,
        dataY = dataY,
        xScalar = xScalar,
        yScalar = yScalar,
        af = af,
        n = n,
        credibleSetArgs = credibleSetArgs,
        fitRetention = fitRetention,
        otherQuantities = otherQuantities,
        region = region,
        priorEffTol = priorEffTol,
        csInput = csInput,
        conditionIdx = conditionIdx
    ))
}

# Post-process each method's fit once (buildTopLoci per fit); the per-method
# 22-column contributions are row-bound later into the single top_loci table.
.ppFitsPerMethod <- function(
    fits,
    dataX,
    dataY,
    xScalar,
    yScalar,
    af,
    n,
    credibleSetArgs,
    fitRetention,
    otherQuantities,
    region,
    priorEffTol,
    csInput,
    conditionIdx
) {
    posts <- map(
        names(fits),
        .ppOneFit,
        fits = fits,
        dataX = dataX,
        dataY = dataY,
        xScalar = xScalar,
        yScalar = yScalar,
        af = af,
        n = n,
        credibleSetArgs = credibleSetArgs,
        fitRetention = fitRetention,
        otherQuantities = otherQuantities,
        region = region,
        priorEffTol = priorEffTol,
        csInput = csInput,
        conditionIdx = conditionIdx
    )
    set_names(posts, names(fits))
}

# Row-bind the per-method top_loci tables and drop them from the per-method
# entries; returns the final finemappingResults + combined top_loci.
.ppFitsCombine <- function(posts) {
    perMethod <- compact(map(posts, "top_loci"))
    topLoci <- if (length(perMethod) == 0L) {
        .emptyTopLoci()
    } else {
        bind_rows(perMethod)
    }
    posts <- map(posts, .ppDropTopLoci)
    list(finemappingResults = posts, top_loci = topLoci)
}

postprocessFinemappingFit <- function(fit, ...) {
    UseMethod("postprocessFinemappingFit")
}

#' @exportS3Method
postprocessFinemappingFit.susie <- function(
    fit,
    method = "susie",
    csInput = NULL,
    ...
) {
    if (is.null(csInput)) {
        csInput <- "X"
    }
    .postprocessFinemappingFitCommon(
        fit,
        method = method,
        csInput = csInput,
        ...
    )
}

#' @exportS3Method
postprocessFinemappingFit.susieInf <- function(
    fit,
    method = "susieInf",
    csInput = NULL,
    ...
) {
    if (is.null(csInput)) {
        csInput <- "X"
    }
    .postprocessFinemappingFitCommon(
        fit,
        method = method,
        csInput = csInput,
        ...
    )
}

#' @exportS3Method
postprocessFinemappingFit.susieRss <- function(
    fit,
    method = "susieRss",
    csInput = NULL,
    ...
) {
    if (is.null(csInput)) {
        csInput <- "Xcorr"
    }
    .postprocessFinemappingFitCommon(
        fit,
        method = method,
        csInput = csInput,
        ...
    )
}

#' @exportS3Method
postprocessFinemappingFit.mvsusie <- function(
    fit,
    method = "mvsusie",
    csInput = NULL,
    ...
) {
    if (is.null(csInput)) {
        csInput <- "X"
    }
    .postprocessFinemappingFitCommon(
        fit,
        method = method,
        csInput = csInput,
        ...
    )
}

#' @exportS3Method
postprocessFinemappingFit.susiF <- function(
    fit,
    method = "fsusie",
    csInput = NULL,
    ...
) {
    if (is.null(csInput)) {
        csInput <- "fsusie"
    }
    .postprocessFinemappingFitCommon(
        fit,
        method = method,
        csInput = csInput,
        ...
    )
}

.postprocessFinemappingFitCommon <- function(
    fit,
    method,
    dataX,
    dataY = NULL,
    xScalar = 1,
    yScalar = 1,
    af = NULL,
    n = NULL,
    credibleSetArgs = CredibleSetParam(),
    fitRetention = "slim",
    otherQuantities = NULL,
    region = NULL,
    priorEffTol = 1e-9,
    conditionIdx = NULL,
    csInput = c("X", "Xcorr", "fsusie")
) {
    csInput <- arg_match(csInput)
    # "slim" keeps a trimmed view of the fit; "full" the whole susie()
    # return, so getSusieFit() and non-default-coverage getCs() can read
    # the full posterior matrices.
    trim <- identical(arg_match(fitRetention, c("slim", "full")), "slim")
    variantNames <- extractVariantNames(fit)
    sumstats <- extractSumstats(dataX, dataY, xScalar, yScalar, method)
    csTables <- .ppCsTables(
        csInput,
        fit = fit,
        dataX = dataX,
        credibleSetArgs = credibleSetArgs,
        method = method
    )
    .ppFinish(
        csTables,
        variantNames,
        sumstats,
        fit = fit,
        method = method,
        af = af,
        n = n,
        dataY = dataY,
        otherQuantities = otherQuantities,
        region = region,
        conditionIdx = conditionIdx,
        credibleSetArgs = credibleSetArgs,
        trim = trim,
        priorEffTol = priorEffTol
    )
}

# Build the canonical top-loci table and wrap it into the postprocess
# result. The table is always built unfiltered -- the FineMappingRow stores
# it as-is so accessors can filter by PIP at query time -- which is why
# `signalCutoff` is read here for the RESULT only, not for the table.
# @noRd
.ppFinish <- function(
    csTables,
    variantNames,
    sumstats,
    fit,
    method,
    af,
    n,
    dataY,
    otherQuantities,
    region,
    conditionIdx,
    credibleSetArgs,
    trim,
    priorEffTol
) {
    topLociFull <- .ppTopLoci(
        csTables,
        variantNames,
        sumstats,
        fit = fit,
        method = method,
        af = af,
        n = n,
        dataY = dataY,
        otherQuantities = otherQuantities,
        region = region,
        conditionIdx = conditionIdx,
        credibleSetArgs = credibleSetArgs
    )
    .ppEntryAndResult(
        topLociFull,
        csTables,
        variantNames,
        sumstats,
        fit = fit,
        method = method,
        dataY = dataY,
        otherQuantities = otherQuantities,
        signalCutoff = credibleSetArgs$signalCutoff,
        trim = trim,
        priorEffTol = priorEffTol
    )
}

# Wrap the finished tables into a FineMappingRow and the postprocess result
# around it. The stored fit is built here rather than earlier because `trim`
# decides how much of it survives, and nothing above this point reads it.
# @noRd
.ppEntryAndResult <- function(
    topLociFull,
    csTables,
    variantNames,
    sumstats,
    fit,
    method,
    dataY,
    otherQuantities,
    signalCutoff,
    trim,
    priorEffTol
) {
    fmEntry <- fineMappingRow(
        variantIds = variantNames,
        susieFit = .ppStoredFit(
            csTables,
            fit = fit,
            trim = trim,
            priorEffTol = priorEffTol,
            method = method
        ),
        topLoci = topLociFull
    )
    .ppAssembleRes(
        topLociFull,
        fmEntry,
        sumstats,
        fit = fit,
        method = method,
        dataY = dataY,
        otherQuantities = otherQuantities,
        signalCutoff = signalCutoff
    )
}

# The fit as stored: trim = TRUE keeps a minimal subset, FALSE the full
# untrimmed susie return (mu / mu2 / lbf_variable / V / ...).
# @noRd
.ppStoredFit <- function(
    csTables,
    fit,
    trim,
    priorEffTol,
    method
) {
    if (!isTRUE(trim)) {
        return(fit)
    }
    trimFinemappingFit(
        fit,
        selectEffects(fit, priorEffTol),
        method,
        csTables
    )
}

# Credible-set tables for the fit at the requested coverages.
.ppCsTables <- function(
    csInput,
    fit,
    dataX,
    credibleSetArgs,
    method
) {
    computeCsTables(
        fit,
        dataX = dataX,
        coverage = credibleSetArgs$coverage,
        secondaryCoverage = credibleSetArgs$secondaryCoverage,
        method = method,
        csInput = csInput,
        minAbsCorr = credibleSetArgs$minAbsCorr,
        medianAbsCorr = credibleSetArgs$medianAbsCorr
    )
}

# Canonical unfiltered top-loci table (signalCutoff = 0).
.ppTopLoci <- function(
    csTables,
    variantNames,
    sumstats,
    fit,
    method,
    af,
    n,
    dataY,
    otherQuantities,
    region,
    conditionIdx,
    credibleSetArgs
) {
    buildTopLoci(
        fit,
        csTables,
        variantNames = variantNames,
        sumstats = sumstats,
        af = af,
        n = n,
        method = method,
        signalCutoff = 0,
        dataY = dataY,
        otherQuantities = otherQuantities,
        region = region,
        conditionIdx = conditionIdx,
        credibleSetArgs = credibleSetArgs
    )
}

# Assemble the wrapper-facing result: PIP-filtered top_loci (legacy behaviour
# for non-S4 callers) + the entry + optional sumstats/sampleNames/context.
.ppAssembleRes <- function(
    topLociFull,
    fmEntry,
    sumstats,
    fit,
    method,
    dataY,
    otherQuantities,
    signalCutoff
) {
    filtering <- !is.null(signalCutoff) &&
        signalCutoff > 0 &&
        nrow(topLociFull) > 0L
    topLociWrapper <- if (!filtering) {
        topLociFull
    } else {
        topLociFull[
            !is.na(topLociFull$pip) & topLociFull$pip > signalCutoff,
            ,
            drop = FALSE
        ]
    }
    c(
        list(
            top_loci = topLociWrapper,
            finemappingEntry = fmEntry,
            method = method
        ),
        compact(list(
            sumstats = sumstats,
            sampleNames = .sampleNamesFromDataY(dataY),
            contextNames = if (method == "mvsusie") fit$outcome_names,
            otherQuantities = otherQuantities
        ))
    )
}

extractVariantNames <- function(fit) {
    variantNames <- names(fit$pip) %||%
        colnames(fit$alpha) %||%
        str_c("variant_", seq_along(fit$pip))
    try_fetch(
        normalizeVariantId(variantNames),
        error = function(cnd) {
            msg <- glue(
                "variant ids could not be normalised; using them as given. ",
                "Downstream joins that assume the canonical form may not ",
                "match."
            )
            warn(msg, parent = cnd)
            variantNames
        }
    )
}

extractSumstats <- function(
    dataX,
    dataY,
    xScalar = 1,
    yScalar = 1,
    method = "susie"
) {
    if (is.null(dataY)) {
        return(NULL)
    }
    if (method == "susieRss") {
        return(dataY)
    }
    if (
        is.list(dataY) &&
            !is.data.frame(dataY) &&
            any(is_in(c("betahat", "sebetahat", "z"), names(dataY)))
    ) {
        return(dataY)
    }
    if (is.null(dataX)) {
        return(NULL)
    }
    if (is.matrix(dataY) || is.data.frame(dataY)) {
        if (ncol(as.matrix(dataY)) != 1) return(NULL)
    }
    sumstats <- univariate_regression(dataX, dataY)
    yScalar <- if (is.null(yScalar) || all(yScalar == 1)) 1 else yScalar
    xScalar <- if (is.null(xScalar) || all(xScalar == 1)) 1 else xScalar
    scale <- yScalar / xScalar
    list_assign(
        sumstats,
        betahat = sumstats$betahat * scale,
        sebetahat = sumstats$sebetahat * scale
    )
}

.sampleNamesFromDataY <- function(dataY) {
    if (is.null(dataY) || is.list(dataY)) {
        return(NULL)
    }
    rownames(as.matrix(dataY))
}

selectEffects <- function(fit, priorEffTol = 1e-9) {
    alpha <- .asEffectMatrix(fit$alpha)
    nEffects <- nrow(alpha)
    if (nEffects == 0) {
        return(integer(0))
    }
    if (!is.null(fit$V)) {
        which(fit$V > priorEffTol)
    } else {
        seq_len(nEffects)
    }
}

.asEffectMatrix <- function(x) {
    if (is.null(x)) {
        return(matrix(numeric(0), nrow = 0))
    }
    if (is.list(x) && !is.data.frame(x)) {
        return(exec(rbind, !!!x))
    }
    as.matrix(x)
}

.asLbfMatrix <- function(fit) {
    if (!is.null(fit$lbf_variable)) {
        return(.asEffectMatrix(fit$lbf_variable))
    }
    if (!is.null(fit$lBF)) {
        return(.asEffectMatrix(fit$lBF))
    }
    NULL
}

#' Compute Credible-Set Tables From a Fine-Mapping Fit
#'
#' Build the per-coverage, purity-filtered credible-set tables from a SuSiE /
#' fSuSiE fit and its design matrix. These tables are the \code{csTables} input
#' to \code{\link{buildTopLoci}}.
#'
#' @param fit A SuSiE-family fit (e.g. from \code{susieR::susie}) carrying
#'   \code{sets}, \code{pip}, and (for fSuSiE) LBF fields.
#' @param dataX Numeric genotype / design matrix (samples x variants) the fit
#'   was computed on; used to assess credible-set purity.
#' @param coverage Numeric primary coverage level, or \code{NULL} to use the
#'   fit's requested coverage (falling back to \code{0.95}).
#' @param secondaryCoverage Numeric vector of additional coverage levels.
#' @param method Character method token (e.g. \code{"susie"}, \code{"fsusie"}).
#' @param csInput One of \code{"X"}, \code{"Xcorr"}, \code{"fsusie"}: how
#'   credible-set purity is computed.
#' @param minAbsCorr,medianAbsCorr Purity thresholds (minimum and median
#'   absolute correlation).
#' @return A named list of per-coverage credible-set tables.
#' @importFrom susieR get_cs_correlation
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:40]
#' y <- eqtlRegionExample$yRes
#' fit <- susieR::susie(X, y, L = 5)
#' computeCsTables(fit, dataX = X, method = "susie")
#' @export
computeCsTables <- function(
    fit,
    dataX,
    coverage = NULL,
    secondaryCoverage = c(0.7, 0.5),
    method = "susie",
    csInput = c("X", "Xcorr", "fsusie"),
    minAbsCorr = 0.8,
    medianAbsCorr = NULL
) {
    csInput <- arg_match(csInput)
    primaryCoverage <- coverage %||% fit$sets$requested_coverage %||% 0.95
    coverages <- discard(
        unique(c(primaryCoverage, secondaryCoverage)),
        is.na
    )

    tables <- map(
        coverages,
        .computeCsTableForCov,
        fit = fit,
        dataX = dataX,
        csInput = csInput,
        minAbsCorr = minAbsCorr,
        medianAbsCorr = medianAbsCorr
    )
    named <- set_names(
        tables,
        map_chr(coverages, formatCsColumn, method = method)
    )
    `attr<-`(named, "coverage", coverages)
}

computeCsTable <- function(
    fit,
    dataX,
    coverage,
    csInput = c("X", "Xcorr", "fsusie"),
    minAbsCorr = 0.8,
    medianAbsCorr = NULL
) {
    csInput <- arg_match(csInput)
    if (csInput == "fsusie") {
        return(.csTableFsusie(fit, dataX, coverage))
    }
    .csTableSusie(fit, dataX, coverage, csInput, minAbsCorr, medianAbsCorr)
}

# The credible sets with the within-CS purity attached, when fsusieR is
# available and returns one purity value per set.
# @noRd
.fsusieSetsWithPurity <- function(sets, dataX) {
    if (!requireNamespace("fsusieR", quietly = TRUE)) {
        return(sets)
    }
    purity <- try_fetch(
        # cal_purity returns one length-1 numeric per credible set, so
        # list_c() is exactly equivalent and refuses a non-numeric or
        # non-scalar element instead of silently producing a longer
        # vector that the length check below would then reject.
        as.numeric(list_c(fsusieR::cal_purity(sets$cs, dataX))),
        error = function(cnd) NULL
    )
    if (is.null(purity) || length(purity) != length(sets$cs)) {
        return(sets)
    }
    list_assign(sets, purity = tibble(min.abs.corr = purity))
}

# fSuSiE credible sets: purity is the min |correlation| WITHIN each CS
# (fsusieR::cal_purity), recorded as sets$purity$min.abs.corr for the canonical
# .csPurityVec() reader; cs_corr keeps the BETWEEN-CS correlation matrix.
.csTableFsusie <- function(fit, dataX, coverage) {
    sets <- try_fetch(
        fsusieGetCs(fit, dataX, requestedCoverage = coverage),
        error = function(cnd) list(cs = list(), requested_coverage = coverage)
    )
    if (
        is.null(sets$cs) ||
            length(sets$cs) == 0 ||
            all(map_lgl(sets$cs, is.null))
    ) {
        return(list(sets = list_assign(sets, cs = list()), pip = fit$pip))
    }
    list(sets = .fsusieSetsWithPurity(sets, dataX), pip = fit$pip)
}

# susieR credible sets from X (correlation computed on genotypes) or Xcorr
# (precomputed LD). min_abs_corr / median_abs_corr passed only when set.
.csTableSusie <- function(
    fit,
    dataX,
    coverage,
    csInput,
    minAbsCorr,
    medianAbsCorr
) {
    csArgs <- c(
        list(coverage = coverage),
        compact(list(
            min_abs_corr = minAbsCorr,
            median_abs_corr = medianAbsCorr
        ))
    )
    # X vs Xcorr only changes how susie_get_cs computes purity; the between-CS
    # correlation is derived on demand later by computeCsCorrelation(), so it is
    # no longer stored on the fit.
    ldArg <- if (csInput == "X") list(X = dataX) else list(Xcorr = dataX)
    sets <- exec(susie_get_cs, !!!c(list(fit), csArgs, ldArg))
    list(sets = sets, pip = fit$pip)
}

# --- computeCsCorrelation: between-CS correlation, derived on demand ----------
# The between-CS correlation is a view over the fit-time LD, which lives on the
# QtlDataset (genotypes) / SumStats (LD sketch) -- it is NEVER stored on the
# fit. get_cs_correlation() needs only the CS membership + PIP + the LD, with
# the LD columns/rows ALIGNED to the fit's variable order (getVariantIds).

# TRUE when the fit has fewer than two credible sets (no between-CS corr).
.csCountBelowTwo <- function(fit) {
    is.null(fit$sets) || is.null(fit$sets$cs) || length(fit$sets$cs) < 2L
}

# GRanges spanning the fit's variants (parsed from chrom:pos in the ids).
.csVariantRegion <- function(variantIds) {
    parts <- str_split(variantIds, ":", simplify = TRUE)
    chrom <- unique(parts[, 1L])
    if (length(chrom) != 1L) {
        abort(glue(
            "computeCsCorrelation(): the fit variants span multiple ",
            "chromosomes ({str_flatten(chrom, ', ')})."
        ))
    }
    pos <- as.integer(parts[, 2L])
    GenomicRanges::GRanges(
        seqnames = chrom,
        ranges = IRanges::IRanges(start = min(pos), end = max(pos))
    )
}

# QtlDataset genotypes for the fit's region, aligned to the fit's variable
# order; errors if any fit variant is absent (a missing one would misalign the
# 1..p credible-set indices with a shrunken genotype matrix).
.csGenotypesForFit <- function(qtlDataset, variantIds) {
    geno <- getGenotypes(qtlDataset, region = .csVariantRegion(variantIds))
    absent <- setdiff(variantIds, colnames(geno))
    if (length(absent) > 0L) {
        abort(glue(
            "computeCsCorrelation(): {length(absent)} fit variant(s) absent ",
            "from the QtlDataset genotypes."
        ))
    }
    geno[, variantIds, drop = FALSE]
}

#' @rdname computeCsCorrelation
setMethod(
    "computeCsCorrelation",
    signature(x = "FineMappingResultBase", ldSource = "SumStatsBase"),
    function(x, ldSource) {
        .rowCsCorrelationSumstats(.asFmRowPayload(x), ldSource)
    }
)

# @noRd
.rowCsCorrelationSumstats <- function(parts, ldSource) {
    {
        fit <- .fmrPartsSusieFit(parts)
        if (.csCountBelowTwo(fit)) {
            return(NULL)
        }
        ldSketch <- getLdSketch(ldSource)
        if (is.null(ldSketch)) {
            abort(glue(
                "computeCsCorrelation(): the summary-statistics ldSource ",
                "carries no LD sketch to derive the correlation from."
            ))
        }
        # onMissing = "error": every fit variant must be in the panel, else the
        # 1..p sets$cs indices would misalign with a shrunken LD matrix.
        xcorr <- .ldFromSketch(
            ldSketch,
            .fmrPartsVariantIds(parts),
            label = "computeCsCorrelation",
            onMissing = "error"
        )
        get_cs_correlation(list(sets = fit$sets, pip = fit$pip), Xcorr = xcorr)
    }
}

# Individual-level LD source: genotypes -> aligned X. susie fits derive the LD
# via get_cs_correlation(X = ); fSuSiE fits (class "susiF") via cal_cor_cs().
#' @rdname computeCsCorrelation
setMethod(
    "computeCsCorrelation",
    signature(x = "FineMappingResultBase", ldSource = "QtlDataset"),
    function(x, ldSource) {
        .rowCsCorrelationGeno(.asFmRowPayload(x), ldSource)
    }
)

# @noRd
.rowCsCorrelationGeno <- function(parts, ldSource) {
    {
        fit <- .fmrPartsSusieFit(parts)
        if (.csCountBelowTwo(fit)) {
            return(NULL)
        }
        geno <- .csGenotypesForFit(ldSource, .fmrPartsVariantIds(parts))
        if (inherits(fit, "susiF")) {
            if (!requireNamespace("fsusieR", quietly = TRUE)) {
                abort(glue(
                    "computeCsCorrelation(): the fit is an fSuSiE object but ",
                    "fsusieR is not installed."
                ))
            }
            fsusieR::cal_cor_cs(fit, geno)$cs_cor
        } else {
            get_cs_correlation(list(sets = fit$sets, pip = fit$pip), X = geno)
        }
    }
}

# Per-effect (per credible set) variant-level columns from the susie fit. Always
# returns `within_cs_pip` (the variant's alpha in the single effect of its
# assigned primary-coverage CS; NA for non-CS variants -- alpha is a
# probability,
# no scaling). With fullFit = TRUE it also widens the per-effect matrices, one
# column set per CS: `within_cs_pip_<lab>` (alpha) and -- unless
# perCsColumns = "full" -- `cs_logbf_<lab>` (lbf_variable),
# `cs_effect_<lab>` (mu /
# X_column_scale_factors) and `cs_effect_var_<lab>` ((mu2 - mu^2) / scale^2).
# includeAllCs = TRUE widens EVERY effect (label `L<k>`), else only effects that
# produced a passing CS (label `cs<pos>`, matching the cs_<cov> columns).
# alpha/mu/mu2/lbf are L x p per-effect matrices (mu/mu2 already
# condition-sliced upstream); missing on a trimmed / fSuSiE fit, in which case
# the values are NA.
# @noRd
.fullFitColumns <- function(
    alpha,
    mu,
    mu2,
    lbfMat,
    scale,
    primaryCsPos,
    effectOf,
    credibleSetArgs = CredibleSetParam(includeAllCs = FALSE)
) {
    perCs <- credibleSetArgs$perCsColumns %||% "none"
    nV <- if (is.null(alpha) || length(dim(alpha)) < 2L) {
        length(primaryCsPos)
    } else {
        ncol(alpha)
    }
    hasAlpha <- !is.null(alpha) && length(dim(alpha)) == 2L && nrow(alpha) > 0L
    withinPip <- .ffcWithinPip(alpha, primaryCsPos, effectOf, nV, hasAlpha)
    cols <- tibble(within_cs_pip = withinPip)
    if (identical(perCs, "none") || !hasAlpha) {
        return(cols)
    }
    .ffcWideColumns(
        cols,
        alpha,
        mu,
        mu2,
        lbfMat,
        scale,
        effectOf,
        nV,
        identical(perCs, "alpha"),
        credibleSetArgs$includeAllCs
    )
}

# Per-variant PIP within its primary-coverage credible set (NA outside any CS).
.ffcWithinPip <- function(alpha, primaryCsPos, effectOf, nV, hasAlpha) {
    if (!(hasAlpha && length(primaryCsPos) == nV && length(effectOf) > 0L)) {
        return(rep(NA_real_, nV))
    }
    # Indexing with NA yields NA, so an out-of-range credible-set position
    # drops out on its own rather than needing a per-variant branch.
    csOk <- !is.na(primaryCsPos) &
        primaryCsPos >= 1L &
        primaryCsPos <= length(effectOf)
    effect <- effectOf[if_else(csOk, as.integer(primaryCsPos), NA_integer_)]
    effectOk <- !is.na(effect) & effect >= 1L & effect <= nrow(alpha)
    if (!any(effectOk)) {
        return(rep(NA_real_, nV))
    }
    # Row 1 stands in wherever the effect is unusable; those entries are
    # masked back to NA immediately, and it keeps the index in bounds.
    picked <- alpha[cbind(if_else(effectOk, effect, 1L), seq_len(nV))]
    if_else(effectOk, picked, NA_real_)
}

# Wide per-effect columns: within_cs_pip_<lab> always, plus cs_logbf / cs_effect
# / cs_effect_var when perCsColumns is "full". Effects come from every effect
# (includeAllCs) or only the credible-set effects.
.ffcWideColumns <- function(
    cols,
    alpha,
    mu,
    mu2,
    lbfMat,
    scale,
    effectOf,
    nV,
    alphaOnly,
    includeAllCs
) {
    if (isTRUE(includeAllCs)) {
        effs <- seq_len(nrow(alpha))
        labs <- str_c("L", effs)
    } else {
        keep <- which(
            !is.na(effectOf) & effectOf >= 1L & effectOf <= nrow(alpha)
        )
        effs <- effectOf[keep]
        labs <- str_c("cs", keep)
    }
    if (is.null(scale) || length(scale) != nV) {
        scale <- rep(1, nV)
    }
    perEffect <- map(
        seq_along(effs),
        .ffcEffectColumns,
        effs = effs,
        labs = labs,
        alpha = alpha,
        mu = mu,
        mu2 = mu2,
        lbfMat = lbfMat,
        scale = scale,
        alphaOnly = alphaOnly
    )
    mutate(cols, !!!.fmwConcat(perEffect))
}

# Concatenate per-item lists, empty-safe.
# @noRd
.fmwConcat <- function(pieces) {
    if (length(pieces) == 0L) {
        return(list())
    }
    list_c(pieces)
}

# One effect's wide columns: the within-CS PIP always, plus the log-BF and
# effect columns when the fit carries them.
# @noRd
.ffcEffectColumns <- function(
    i,
    effs,
    labs,
    alpha,
    mu,
    mu2,
    lbfMat,
    scale,
    alphaOnly
) {
    L <- effs[[i]]
    lab <- labs[[i]]
    pip <- set_names(
        list(unname(alpha[L, ])),
        str_c("within_cs_pip_", lab)
    )
    if (isTRUE(alphaOnly)) {
        return(pip)
    }
    logbf <- if (!is.null(lbfMat) && L <= nrow(lbfMat)) {
        set_names(list(unname(lbfMat[L, ])), str_c("cs_logbf_", lab))
    } else {
        list()
    }
    effect <- if (!is.null(mu) && L <= nrow(mu)) {
        set_names(list(unname(mu[L, ] / scale)), str_c("cs_effect_", lab))
    } else {
        list()
    }
    effectVar <- if (!is.null(mu) && !is.null(mu2) && L <= nrow(mu2)) {
        set_names(
            list(unname((mu2[L, ] - mu[L, ]^2) / scale^2)),
            str_c("cs_effect_var_", lab)
        )
    } else {
        list()
    }
    c(pip, logbf, effect, effectVar)
}

# Slice a susie posterior array to the active condition (3-D fit) or coerce a
# 2-D array to matrix; NULL for a 3-D fit with no conditionIdx.
# @noRd
.fmSliceCond <- function(arr, conditionIdx) {
    if (is.null(arr)) {
        return(NULL)
    }
    if (length(dim(arr)) == 3L) {
        if (is.null(conditionIdx)) {
            return(NULL)
        }
        return(as.matrix(arr[,, conditionIdx]))
    }
    as.matrix(arr)
}

# Per-variant CS index at coverage `targetCov` (0 = not in any CS; on overlap
# the smallest cs_idx wins).
# @noRd
.fmCsIdxAtCoverage <- function(
    targetCov,
    coverageValues,
    csTables,
    nV,
    variantNames = NULL
) {
    empty <- integer(nV)
    hit <- which(abs(coverageValues - targetCov) < 1e-12)
    if (length(hit) == 0L) {
        return(empty)
    }
    sets <- csTables[[hit[1L]]]$sets$cs
    if (is.null(sets) || length(sets) == 0L) {
        return(empty)
    }
    # A variant in several sets goes to the SMALLEST containing set (ties ->
    # lowest position, so the answer is deterministic), not to whichever set
    # happened to come first in the list. Every membership is recorded so the
    # ambiguity can be reported rather than silently resolved.
    # Every (variant, set) membership, flattened. Ordering the memberships by
    # set size then set index makes the first one per variant the winner --
    # smallest set, ties to the lowest position -- which is exactly what the
    # running "is this smaller than the best so far" comparison decided.
    setSizes <- lengths(sets)
    memberships <- map(seq_along(sets), .fmCsMemberships, sets = sets, nV = nV)
    variantOf <- .fmwConcatInt(map(memberships, "variant"))
    csOf <- .fmwConcatInt(map(memberships, "cs"))
    if (length(variantOf) == 0L) {
        return(empty)
    }
    ord <- order(variantOf, setSizes[csOf], csOf)
    firstPerVariant <- ord[!duplicated(variantOf[ord])]
    winners <- csOf[firstPerVariant]
    winnerAt <- variantOf[firstPerVariant]
    out <- .fmScatter(empty, winnerAt, winners)
    bestSize <- .fmScatter(rep(Inf, nV), winnerAt, setSizes[winners])
    memb <- split(csOf, factor(variantOf, levels = seq_len(nV)))
    .fmWarnMultiCs(memb, out, bestSize, sets, variantNames)
    out
}

# Place `values` at positions `at` in `base`. The one write the scatter needs,
# named so it reads as a total operation rather than an accumulation.
# @noRd
.fmScatter <- function(base, at, values) {
    replace(base, at, values)
}

# @noRd
.fmwConcatInt <- function(pieces) {
    if (length(pieces) == 0L) {
        return(integer(0))
    }
    list_c(pieces)
}

# Set `csIdx`'s in-range variant memberships, as parallel (variant, cs) runs.
# @noRd
.fmCsMemberships <- function(csIdx, sets, nV) {
    raw <- as.integer(sets[[csIdx]])
    vi <- raw[raw >= 1L & raw <= nV]
    list(variant = vi, cs = rep(csIdx, length(vi)))
}

# Name the variants that fell in more than one credible set, and which set won.
# @noRd
.fmWarnMultiCs <- function(memb, out, bestSize, sets, variantNames) {
    effIdx <- .fmEffectIndices(sets)
    for (v in which(lengths(memb) > 1L)) {
        vname <- if (!is.null(variantNames) && v <= length(variantNames)) {
            variantNames[v]
        } else {
            str_c("#", v)
        }
        warn(glue(
            "Variant {vname} is in multiple credible sets: ",
            "{str_flatten(str_c('L', effIdx[memb[[v]]]), ', ')}. ",
            "Keeping the smallest: CS L{effIdx[out[v]]} ",
            "(size {bestSize[v]})."
        ))
    }
    invisible(NULL)
}

# The fit's true credible-set effect index (the k behind "L<k>") for each set
# POSITION in `sets$cs`; falls back to the position when the names are absent
# or unparseable. This is what the cs_<cov> label carries, so a gapped effect
# index (L1, L2, L4) survives instead of being renumbered 1..n.
# @noRd
.fmEffectIndices <- function(sets) {
    nm <- names(sets)
    if (is.null(nm)) {
        return(seq_along(sets))
    }
    e <- suppressWarnings(as.integer(str_remove(nm, "^L")))
    bad <- is.na(e)
    if (!any(bad)) {
        return(e)
    }
    replace(e, bad, seq_along(sets)[bad])
}

# Map a per-variant set-POSITION vector onto the fit's true effect indices;
# non-CS entries (0) stay 0.
# @noRd
.fmEffectIdxAtCoverage <- function(
    targetCov,
    posVec,
    coverageValues,
    csTables
) {
    hit <- which(abs(coverageValues - targetCov) < 1e-12)
    if (length(hit) == 0L) {
        return(posVec)
    }
    sets <- csTables[[hit[1L]]]$sets$cs
    if (is.null(sets) || length(sets) == 0L) {
        return(posVec)
    }
    effIdx <- .fmEffectIndices(sets)
    nz <- posVec > 0L
    replace(integer(length(posVec)), nz, effIdx[posVec[nz]])
}

# Per-variant CS purity (min.abs.corr) at coverage `targetCov`; 0 for non-CS
# variants.
# @noRd
.fmPurityAtCoverage <- function(targetCov, idxVec, coverageValues, csTables) {
    h <- which(abs(coverageValues - targetCov) < 1e-12)
    pv <- if (length(h) > 0L) .csPurityVec(csTables[[h[1L]]]) else numeric()
    map_dbl(idxVec, .csPurityAt, pv = pv)
}

#' Build the unified top-loci table for one fit and one method
#'
#' Returns the per-fit, per-method contribution to the unified \code{top_loci}
#' table in the fixed 22-column shape. \code{postprocessFinemappingFits()} calls
#' this once per method per fit and row-binds the results into the single
#' \code{top_loci} returned by \code{formatFinemappingOutput()}.
#'
#' Output columns, in order: \code{#chr}, \code{start}, \code{end}, \code{a1},
#' \code{a2}, \code{variant}, \code{gene}, \code{event}, \code{n}, \code{af},
#' \code{beta}, \code{se}, \code{pip}, \code{posterior_effect_mean},
#' \code{posterior_effect_se}, \code{cs_95}, \code{cs_70}, \code{cs_50},
#' \code{cs_95_purity}, \code{method}, \code{grange_start}, \code{grange_end}.
#'
#' \code{cs_95} / \code{cs_70} / \code{cs_50} are character strings of the form
#' \code{"<method>_<cs_index>"} where each method numbers credible sets
#' independently from 1. Variants retained by the PIP cutoff but not in any
#' credible set at a coverage carry \code{"<method>_0"}. \code{cs_95_purity} is
#' the 0.95-coverage purity for the row's \code{(method, cs_95)}; rows whose
#' \code{cs_95} is \code{"<method>_0"} carry \code{0}.
#'
#' Row uniqueness is \code{(variant, gene, cs_membership)} at the given
#' \code{method}; overlapping CS within the same method produces one row per CS.
#'
#' @param fit Fitted SuSiE-family object (must expose \code{alpha}, \code{mu},
#'   \code{mu2}, \code{pip}).
#' @param csTables List of CS tables (one per coverage) from
#'   \code{computeCsTables()}.
#' @param variantNames Character vector of variant IDs (\code{chr:pos:A2:A1}).
#' @param sumstats Optional marginal-association summary (\code{betahat},
#'   \code{sebetahat}) filling \code{beta} / \code{se}.
#' @param af Optional numeric vector of effect-allele frequencies (frequency of
#'   the final effect allele / \code{a1} after allele harmonization against the
#'   LD/reference variants). Exported directly as the \code{af} column. MAF is
#'   never exported; derive it from \code{af} at filter time. Default NULL ->
#'   \code{af = NA_real_}.
#' @param n Optional per-variant sample size, one entry per variant, exported as
#'   the \code{N} column. Used by RSS fine-mapping, where the effective sample
#'   size varies by variant. Default NULL -> \code{N} falls back to the fit's
#'   scalar sample size (the outcome-matrix \code{nrow} for individual-level
#'   fits, \code{NA} otherwise).
#' @param method Method name (e.g. \code{"susie"}, \code{"susieInf"}). Required.
#' @param signalCutoff PIP cutoff for retaining PIP-only (non-CS) variants.
#' @param dataY Optional regional phenotype vector or matrix;
#'   \code{nrow(as.matrix(dataY))} fills \code{n}, \code{colnames(dataY)[1]}
#'   fills \code{gene}. A \code{list} (the RSS \code{list(z = ...)} form)
#'   carries no sample count: \code{n} is left \code{NA} for the per-variant
#'   \code{n} argument to fill.
#' @param otherQuantities Optional list. Default is NULL.
#' @param region Optional \code{"chr:start-end"} string. Default is NULL.
#' @param conditionIdx Integer or \code{NULL}. Index of the conditioned effect
#'   (per-condition output); \code{NULL} for the unconditioned fit.
#' @param credibleSetArgs How credible sets are built and reported, built with
#'   \code{\link{CredibleSetParam}}. \code{perCsColumns} decides which
#'   per-credible-set variant-level columns this table carries, and
#'   \code{includeAllCs} their labels.
#' @return A data frame in the fixed 22-column shape for this fit and method, or
#'   an empty data frame if nothing is retained.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:40]
#' y <- eqtlRegionExample$yRes
#' fit <- susieR::susie(X, y, L = 5)
#' csTables <- computeCsTables(fit, dataX = X, method = "susie")
#' buildTopLoci(fit = fit, csTables = csTables, variantNames = colnames(X),
#'   method = "susie", dataY = y)
#' @export
buildTopLoci <- function(
    fit,
    csTables,
    variantNames,
    sumstats = NULL,
    af = NULL,
    n = NULL,
    method,
    signalCutoff = 0,
    dataY = NULL,
    otherQuantities = NULL,
    region = NULL,
    conditionIdx = NULL,
    credibleSetArgs = CredibleSetParam(includeAllCs = FALSE)
) {
    if (missing(method)) {
        method <- NULL
    }
    .btlValidateMethod(method)
    if (length(variantNames) == 0L) {
        return(.emptyTopLoci())
    }
    .btlBuild(
        fit = fit,
        csTables = csTables,
        variantNames = variantNames,
        sumstats = sumstats,
        af = af,
        n = n,
        method = method,
        signalCutoff = signalCutoff,
        dataY = dataY,
        otherQuantities = otherQuantities,
        region = region,
        conditionIdx = conditionIdx,
        credibleSetArgs = credibleSetArgs
    )
}

# Orchestrate the top-loci table from the buildTopLoci() arguments.
.btlBuild <- function(
    fit,
    csTables,
    variantNames,
    sumstats,
    af,
    n,
    method,
    signalCutoff,
    dataY,
    otherQuantities,
    region,
    conditionIdx,
    credibleSetArgs
) {
    nV <- length(variantNames)
    cov <- .btlCoverage(csTables)
    fc <- .btlFitConstants(dataY, otherQuantities, region)
    post <- .btlPosterior(fit, conditionIdx, nV)
    marg <- .btlMarginal(sumstats, nV)
    cs <- .btlCsMembership(cov, csTables, nV, variantNames)
    fullFitBlock <- .btlFullFitBlock(
        fit,
        post,
        cov,
        cs,
        csTables,
        nV,
        list(credibleSetArgs = credibleSetArgs)
    )
    .btlAssembleFinal(
        variantNames,
        fc = fc,
        marg = marg,
        post = post,
        fit = fit,
        af = af,
        n = n,
        method = method,
        cs = cs,
        cov = cov,
        csTables = csTables,
        fullFitBlock = fullFitBlock,
        nV = nV,
        conditionIdx = conditionIdx,
        signalCutoff = signalCutoff
    )
}

# Assemble the per-variant table, add the conditional columns, and apply
# the signal cutoff. Split from .btlBuild() so that function reads as the
# list of pieces it derives, and this one as what is done with them.
# @noRd
.btlAssembleFinal <- function(
    variantNames,
    fc,
    marg,
    post,
    fit,
    af,
    n,
    method,
    cs,
    cov,
    csTables,
    fullFitBlock,
    nV,
    conditionIdx,
    signalCutoff
) {
    out <- .btlAssemble(
        variantNames,
        .btlParseVariants(variantNames),
        fc,
        marg,
        post,
        fit,
        af,
        n,
        method,
        .btlCsBlock(method, cs, nV),
        fullFitBlock,
        nV
    )
    cond <- .btlConditional(
        fit,
        method,
        conditionIdx,
        cov,
        cs$covSorted,
        csTables,
        nV
    )
    .btlFinalize(out, cond, conditionIdx, signalCutoff)
}

# Attach per-condition columns (multi-condition fits) and apply the PIP cutoff.
.btlFinalize <- function(out, cond, conditionIdx, signalCutoff) {
    withCond <- if (is.null(conditionIdx)) {
        out
    } else {
        mutate(
            out,
            conditional_effect = cond$condEffect,
            lfsr = cond$condLfsr
        )
    }
    if (is.null(signalCutoff) || signalCutoff <= 0) {
        return(withCond)
    }
    filter(withCond, !is.na(.data$pip) & .data$pip > signalCutoff)
}

# buildTopLoci step helpers ---------------------------------------------------

# `method` is required and must be a single non-empty, non-NA string.
#' @importFrom checkmate checkString
.btlValidateMethod <- function(method) {
    res <- checkString(method, min.chars = 1L)
    if (!isTRUE(res)) {
        abort(
            "buildTopLoci: `method` is required (e.g. \"susie\", \"susieInf\")."
        )
    }
}

# Coverage levels attached to csTables (NA-filled when the attribute is absent).
.btlCoverage <- function(csTables) {
    cov <- attr(csTables, "coverage")
    if (is.null(cov)) rep(NA_real_, length(csTables)) else cov
}

# The per-variant N column. Numeric, so a fractional *effective* N (e.g.
# 4 / (1 / nCase + 1 / nControl)) matches getSumStatsDf(entry)$N exactly rather
# than being truncated. With no per-variant n this falls back to the fit's
# scalar N (integer nrow on the QTL path, NA otherwise); a scalar n recycles.
# @noRd
.btlNColumn <- function(n, fitN, nV) {
    if (is.null(n)) {
        return(rep(fitN, nV))
    }
    if (length(n) == 1L) {
        return(rep(as.numeric(n), nV))
    }
    as.numeric(n)
}

# Per-fit constants: sample size, gene (first phenotype column), event id, and
# the parsed genomic range.
.btlFitConstants <- function(dataY, otherQuantities, region) {
    # Only a genuine sample-by-outcome vector/matrix carries a meaningful nrow.
    # RSS passes dataY as a LIST (e.g. list(z = ...)); as.matrix() on a list
    # collapses it to a 1x1 cell, which used to make fitN (and every top_loci
    # N) 1. Excluding lists is the whole rule -- same test .sampleNamesFromDataY
    # applies. Requiring a matrix here is too strict: the individual-level path
    # passes a bare numeric vector of one value per sample, which as.matrix()
    # turns into the n x 1 column it already is, so demanding a matrix silently
    # dropped its N to NA. For the list case fitN stays NA and the per-variant
    # `n` (threaded from the RSS effective sample size) fills the N column.
    hasOutcomeMatrix <- !is.null(dataY) &&
        (is.matrix(dataY) ||
            is.data.frame(dataY) ||
            (is.atomic(dataY) && is.numeric(dataY)))
    dataYMat <- if (hasOutcomeMatrix) as.matrix(dataY) else NULL
    fitN <- if (is.null(dataYMat)) NA_integer_ else as.integer(nrow(dataYMat))
    fitGene <- if (!is.null(dataYMat) && !is.null(colnames(dataYMat))) {
        colnames(dataYMat)[1]
    } else {
        NA_character_
    }
    fitEvent <- if (
        !is.null(otherQuantities$condition_id) &&
            !is.na(fitGene) &&
            str_length(fitGene) > 0L
    ) {
        str_c(otherQuantities$condition_id, fitGene, sep = "_")
    } else {
        NA_character_
    }
    list(
        fitN = fitN,
        fitGene = fitGene,
        fitEvent = fitEvent,
        grange = .parseGrange(region)
    )
}

# Per-variant posterior mean/SD (from alpha, mu, mu2) and the strongest
# single-effect log Bayes factor. A conditionIdx slices 3-D mvsusie mu/mu2.
.btlPosterior <- function(fit, conditionIdx, nV) {
    alpha <- as.matrix(fit$alpha)
    mu <- .fmSliceCond(fit$mu, conditionIdx)
    mu2 <- .fmSliceCond(fit$mu2, conditionIdx)
    postMean <- if (!is.null(mu) && all(dim(alpha) == dim(mu))) {
        colSums(alpha * mu)
    } else {
        rep(NA_real_, nV)
    }
    postSd <- if (!is.null(mu2) && all(dim(alpha) == dim(mu2))) {
        sqrt(pmax(colSums(alpha * mu2) - postMean^2, 0))
    } else {
        rep(NA_real_, nV)
    }
    lbfMat <- .asLbfMatrix(fit)
    logBF <- if (!is.null(lbfMat) && ncol(lbfMat) == nV) {
        apply(lbfMat, 2, .finiteMax)
    } else {
        rep(NA_real_, nV)
    }
    list(
        alpha = alpha,
        mu = mu,
        mu2 = mu2,
        postMean = postMean,
        postSd = postSd,
        logBF = logBF
    )
}

# Parse variant IDs to chrom/pos/A1/A2; error on missing or invalid coordinates.
#' @importFrom rlang try_fetch
.btlParseVariants <- function(variantNames) {
    parsed <- try_fetch(
        suppressWarnings(parseVariantId(variantNames)),
        error = function(cnd) {
            abort("buildTopLoci: parseVariantId failed", parent = cnd)
        }
    )
    if (is.null(parsed) || nrow(parsed) != length(variantNames)) {
        abort(
            "buildTopLoci: parseVariantId did not return one row per variant."
        )
    }
    invalid <- is.na(parsed$chrom) |
        is.na(parsed$pos) |
        is.na(parsed$A1) |
        str_length(parsed$A1) == 0L |
        is.na(parsed$A2) |
        str_length(parsed$A2) == 0L
    if (any(invalid)) {
        badVar <- variantNames[which(invalid)[[1]]]
        msg <- glue(
            "buildTopLoci: parseVariantId produced invalid coordinates ",
            "for variant_id: {badVar}"
        )
        abort(msg)
    }
    parsed
}

# Marginal univariate effects (beta, se, z, p) from the sumstats list; z and p
# are derived when not supplied directly.
.btlMarginal <- function(sumstats, nV) {
    beta <- if (!is.null(sumstats$betahat)) {
        as.numeric(sumstats$betahat)
    } else {
        rep(NA_real_, nV)
    }
    se <- if (!is.null(sumstats$sebetahat)) {
        as.numeric(sumstats$sebetahat)
    } else {
        rep(NA_real_, nV)
    }
    z <- if (!is.null(sumstats$z)) {
        as.numeric(sumstats$z)
    } else if (any(!is.na(beta)) && any(!is.na(se))) {
        beta / se
    } else {
        rep(NA_real_, nV)
    }
    p <- if (!is.null(sumstats$p)) {
        as.numeric(sumstats$p)
    } else if (any(!is.na(z))) {
        2 * stats::pnorm(-abs(z))
    } else {
        rep(NA_real_, nV)
    }
    list(beta = beta, se = se, z = z, p = p)
}

# CS membership index + purity for every coverage present (high -> low), with
# the matching `cs_<coverage*100>` column names.
.btlCsMembership <- function(
    coverageValues,
    csTables,
    nV,
    variantNames = NULL
) {
    covSorted <- sort(
        unique(coverageValues[is.finite(coverageValues)]),
        decreasing = TRUE
    )
    csIdxByCov <- map(
        covSorted,
        .fmCsIdxAtCoverage,
        coverageValues,
        csTables,
        nV,
        variantNames
    )
    # Effect-index view of the same assignment: the cs_<cov> LABEL carries the
    # fit's true (possibly gapped) effect index, while purity and the full-fit
    # block keep indexing by set POSITION via csIdxByCov.
    csEffectByCov <- map2(
        covSorted,
        csIdxByCov,
        .fmEffectIdxAtCoverage,
        coverageValues = coverageValues,
        csTables = csTables
    )
    csPurityByCov <- map2(
        covSorted,
        csIdxByCov,
        .fmPurityAtCoverage,
        coverageValues = coverageValues,
        csTables = csTables
    )
    list(
        covSorted = covSorted,
        csIdxByCov = csIdxByCov,
        csEffectByCov = csEffectByCov,
        csPurityByCov = csPurityByCov,
        csColNames = str_c("cs_", covSorted * 100)
    )
}

# Per-condition conditional effect (coef / pip) for a multi-condition fit.
# Accepts a trimmed fit's `$coef` or a raw mvsusie fit via coef.mvsusie().
.btlCondEffect <- function(fit, method, conditionIdx, nV) {
    coefMat <- if (!is.null(fit$coef)) {
        as.matrix(fit$coef)
    } else if (
        identical(method, "mvsusie") &&
            requireNamespace("mvsusieR", quietly = TRUE)
    ) {
        cm <- try_fetch(mvsusieR::coef.mvsusie(fit), error = function(cnd) NULL)
        if (!is.null(cm)) as.matrix(cm)[-1L, , drop = FALSE] else NULL
    } else {
        NULL
    }
    if (
        is.null(coefMat) || nrow(coefMat) != nV || ncol(coefMat) < conditionIdx
    ) {
        return(rep(NA_real_, nV))
    }
    pipVec <- as.numeric(fit$pip)
    if_else(pipVec > 0, coefMat[, conditionIdx] / pipVec, NA_real_)
}

# Per-condition conditional lfsr: map each variant to its effect (L) via the
# primary (highest-coverage) credible set, then read that effect's lfsr.
.btlCondLfsr <- function(
    fit,
    conditionIdx,
    coverageValues,
    covSorted,
    csTables,
    nV
) {
    clf <- if (!is.null(fit$clfsr)) fit$clfsr else fit$conditional_lfsr
    condLfsr <- rep(NA_real_, nV)
    if (
        is.null(clf) ||
            length(dim(clf)) != 3L ||
            dim(clf)[3L] < conditionIdx ||
            length(covSorted) == 0L
    ) {
        return(condLfsr)
    }
    hPrim <- which(abs(coverageValues - covSorted[1L]) < 1e-12)
    if (length(hPrim) == 0L) {
        return(condLfsr)
    }
    setsPrim <- csTables[[hPrim[1L]]]$sets$cs
    if (is.null(setsPrim) || length(setsPrim) == 0L) {
        return(condLfsr)
    }
    effectOf <- suppressWarnings(as.integer(str_remove(names(setsPrim), "^L")))
    assignments <- compact(map(
        seq_along(setsPrim),
        .btlCondLfsrForSet,
        setsPrim = setsPrim,
        effectOf = effectOf,
        clf = clf,
        conditionIdx = conditionIdx,
        nV = nV
    ))
    # Later sets overwrite earlier ones at a shared variant, as the running
    # assignment did; one scatter replaces the per-set writes.
    .fmScatter(
        condLfsr,
        .fmwConcatInt(map(assignments, "variant")),
        .fmwConcatDbl(map(assignments, "value"))
    )
}

# @noRd
.fmwConcatDbl <- function(pieces) {
    if (length(pieces) == 0L) {
        return(numeric(0))
    }
    list_c(pieces)
}

# One credible set's conditional-lfsr values, or NULL when its effect is out
# of range or it covers no in-range variant.
# @noRd
.btlCondLfsrForSet <- function(
    csPos,
    setsPrim,
    effectOf,
    clf,
    conditionIdx,
    nV
) {
    L <- effectOf[csPos]
    if (is.na(L) || L < 1L || L > dim(clf)[1L]) {
        return(NULL)
    }
    raw <- as.integer(setsPrim[[csPos]])
    vi <- raw[raw >= 1L & raw <= nV]
    if (length(vi) == 0L) {
        return(NULL)
    }
    list(variant = vi, value = as.numeric(clf[L, vi, conditionIdx]))
}

# Per-condition posterior quantities (NA for univariate fits).
.btlConditional <- function(
    fit,
    method,
    conditionIdx,
    coverageValues,
    covSorted,
    csTables,
    nV
) {
    if (is.null(conditionIdx)) {
        return(list(
            condEffect = rep(NA_real_, nV),
            condLfsr = rep(NA_real_, nV)
        ))
    }
    list(
        condEffect = .btlCondEffect(fit, method, conditionIdx, nV),
        condLfsr = .btlCondLfsr(
            fit,
            conditionIdx,
            coverageValues,
            covSorted,
            csTables,
            nV
        )
    )
}

# Dynamic CS block: cs_<C> memberships then cs_<C>_purity, one pair per
# coverage.
.btlCsBlock <- function(method, cs, nV) {
    methodTag <- .camelToSnakeMethod(method)
    csList <- c(
        set_names(
            map(cs$csEffectByCov, .btlCsLabel, methodTag = methodTag),
            cs$csColNames
        ),
        set_names(cs$csPurityByCov, str_c(cs$csColNames, "_purity"))
    )
    if (length(csList) == 0L) {
        return(tibble(.rows = nV))
    }
    as_tibble(csList, .name_repair = "minimal")
}

# The effect indices behind the primary-coverage CS, read off the "L<k>" names
# of its sets$cs. integer(0) when there is no such CS or it names no effects.
# @noRd
.btlPrimaryEffects <- function(coverageValues, cs, csTables) {
    if (length(cs$covSorted) == 0L) {
        return(integer(0))
    }
    hP <- which(abs(coverageValues - cs$covSorted[1L]) < 1e-12)
    if (length(hP) == 0L) {
        return(integer(0))
    }
    spP <- csTables[[hP[1L]]]$sets$cs
    if (is.null(spP) || length(spP) == 0L) {
        return(integer(0))
    }
    suppressWarnings(as.integer(str_remove(names(spP), "^L")))
}

# within_cs_pip (+ optional fullFit-wide) columns, mapping each variant to its
# primary-coverage CS effect (position -> effect via the sets$cs "L<k>" names).
.btlFullFitBlock <- function(
    fit,
    post,
    coverageValues,
    cs,
    csTables,
    nV,
    opts
) {
    primaryCsPos <- if (length(cs$csIdxByCov) > 0L) {
        cs$csIdxByCov[[1L]]
    } else {
        integer(nV)
    }
    effectOfPrim <- .btlPrimaryEffects(coverageValues, cs, csTables)
    .fullFitColumns(
        post$alpha,
        post$mu,
        post$mu2,
        .asLbfMatrix(fit),
        fit$X_column_scale_factors,
        primaryCsPos,
        effectOfPrim,
        credibleSetArgs = opts$credibleSetArgs
    )
}

# Assemble the core per-variant table (identity + marginal + posterior) with the
# CS and fullFit blocks and per-fit metadata.
.btlAssemble <- function(
    variantNames,
    parsed,
    fc,
    marg,
    post,
    fit,
    af,
    n,
    method,
    csBlock,
    fullFitBlock,
    nV
) {
    core <- tibble(
        variant_id = as.character(variantNames),
        chrom = unname(parsed$chrom),
        pos = as.integer(parsed$pos),
        A1 = unname(parsed$A1),
        A2 = unname(parsed$A2),
        N = .btlNColumn(n, fc$fitN, nV),
        af = if (is.null(af)) rep(NA_real_, nV) else as.numeric(af),
        marginal_beta = unname(marg$beta),
        marginal_se = unname(marg$se),
        marginal_z = unname(marg$z),
        marginal_p = unname(marg$p),
        pip = as.numeric(fit$pip),
        posterior_mean = unname(post$postMean),
        posterior_sd = unname(post$postSd),
        logBF = unname(post$logBF)
    )
    meta <- tibble(
        method = rep(method, nV),
        gene = rep(fc$fitGene, nV),
        event = rep(fc$fitEvent, nV),
        grange_start = rep(fc$grange[["start"]], nV),
        grange_end = rep(fc$grange[["end"]], nV)
    )
    bind_cols(core, csBlock, fullFitBlock, meta)
}

# Per-CS purity from one cs_table: susieR's sets$purity$min.abs.corr, or NA
# when the fit carries no purity.
.csPurityVec <- function(ct) {
    sp <- ct$sets$purity
    if (!is.null(sp) && is_in("min.abs.corr", names(sp))) {
        return(as.numeric(sp$min.abs.corr))
    }
    rep(NA_real_, length(ct$sets$cs))
}

.emptyTopLoci <- function() {
    tibble(
        variant_id = character(),
        chrom = character(),
        pos = integer(),
        A1 = character(),
        A2 = character(),
        N = numeric(),
        af = numeric(),
        marginal_beta = numeric(),
        marginal_se = numeric(),
        marginal_z = numeric(),
        marginal_p = numeric(),
        pip = numeric(),
        posterior_mean = numeric(),
        posterior_sd = numeric(),
        logBF = numeric(),
        cs_95 = character(),
        cs_70 = character(),
        cs_50 = character(),
        cs_95_purity = numeric(),
        cs_70_purity = numeric(),
        cs_50_purity = numeric(),
        within_cs_pip = numeric(),
        method = character(),
        gene = character(),
        event = character(),
        grange_start = integer(),
        grange_end = integer()
    )
}

.parseGrange <- function(regionStr) {
    if (
        is.null(regionStr) ||
            length(regionStr) == 0L ||
            is.na(regionStr) ||
            str_length(as.character(regionStr)) == 0L
    ) {
        return(c(start = NA_integer_, end = NA_integer_))
    }
    pr <- try_fetch(
        parseRegion(as.character(regionStr)),
        error = function(cnd) NULL
    )
    if (is.null(pr) || !is.data.frame(pr)) {
        return(c(start = NA_integer_, end = NA_integer_))
    }
    c(start = as.integer(pr$start), end = as.integer(pr$end))
}

trimFinemappingFit <- function(fit, effectIdx, method, csTables) {
    common <- .trimBaseFit(fit, effectIdx, csTables) |>
        .trimAddCommon(fit, effectIdx)
    trimmed <- if (method == "mvsusie") {
        .trimAddMvsusie(common, fit, effectIdx)
    } else {
        common
    }
    # fSuSiE: keep the precomputed variants x features TWAS weight matrix
    # (fsusieWeights output, attached as $coef before trimming) so downstream
    # TWAS can read it without the dropped wavelet slots.
    withCoef <- if (method == "fsusie" && !is.null(fit$coef)) {
        list_assign(trimmed, coef = fit$coef)
    } else {
        trimmed
    }
    `class<-`(withCoef, unique(c(method, "susie")))
}

# The minimal always-kept subset of a susie fit (pip, credible sets, effect
# matrices for the selected effects).
.trimBaseFit <- function(fit, effectIdx, csTables) {
    alpha <- .asEffectMatrix(fit$alpha)
    lbfVariable <- .asLbfMatrix(fit)
    primary <- csTables[[1]]
    secondary <- if (length(csTables) > 1) {
        map(csTables[-1], .dropPipCol)
    } else {
        NULL
    }
    list(
        pip = as.numeric(fit$pip),
        sets = primary$sets,
        sets_secondary = secondary,
        alpha = alpha[effectIdx, , drop = FALSE],
        lbf_variable = if (!is.null(lbfVariable)) {
            lbfVariable[effectIdx, , drop = FALSE]
        } else {
            NULL
        },
        V = if (!is.null(fit$V)) fit$V[effectIdx] else NULL,
        niter = fit$niter,
        max_L = nrow(alpha),
        n_effects = nrow(alpha)
    )
}

# Optional slots common to susie/mvsusie: column scales, posterior mu/mu2
# (L x p, or L x p x R for multivariate), theta, omega_weights.
.trimAddCommon <- function(trimmed, fit, effectIdx) {
    withCommon <- list_assign(
        trimmed,
        !!!compact(list(
            X_column_scale_factors = fit$X_column_scale_factors,
            mu = .trimEffectSlice(fit$mu, effectIdx),
            mu2 = .trimEffectSlice(fit$mu2, effectIdx),
            theta = fit$theta,
            omega_weights = fit$omega_weights
        ))
    )
    .trimAddScalars(withCommon, fit)
}

# The selected effects of a posterior array, which is L x p for a univariate
# fit and L x p x R for a multivariate one.
# @noRd
.trimEffectSlice <- function(x, effectIdx) {
    if (is.null(x)) {
        return(NULL)
    }
    if (length(dim(x)) == 3) {
        return(x[effectIdx, , , drop = FALSE])
    }
    x[effectIdx, , drop = FALSE]
}

# Cheap fields that are not effect-indexed, kept so a trimmed fit is sufficient
# for variant-subset reconciliation:
#
#   sigma2      the moment recovery (post_var = mu2 - mu^2; pw = sigma2 *
#               (1/post_var - 1/V)) that makes an optional refit possible
#   null_index  detects a null_weight fit, whose alpha has p + 1 columns while
#               pip stays length p -- invisible to any pip-length check
#   converged   otherwise a non-converged fit is indistinguishable downstream
#   pi          the fit's true prior weights, for prior-aware refits
#   XtXr        quantifies what the dropped variants contributed
#
# Together ~2p + 3 numbers (~0.8 MB at p = 50k), against the ~10x blow-up of
# storing the fit untrimmed. Read with [[ so `pi` cannot prefix-match `pip`.
# @noRd
.trimAddScalars <- function(trimmed, fit) {
    reduce(
        c("sigma2", "null_index", "converged", "pi", "XtXr"),
        .trimCopyField,
        fit = fit,
        .init = trimmed
    )
}

# Copy one non-effect-indexed field from `fit` onto `trimmed`, when present.
# @noRd
.trimCopyField <- function(trimmed, nm, fit) {
    value <- fit[[nm]]
    if (is.null(value)) {
        return(trimmed)
    }
    list_assign(trimmed, !!!set_names(list(value), nm))
}

# mvsusie-specific slots: per-effect mu2_diag, the coefficient matrix, and the
# conditional lfsr array.
.trimAddMvsusie <- function(trimmed, fit, effectIdx) {
    list_assign(
        trimmed,
        !!!compact(list(
            mu2_diag = if (!is.null(fit$mu2_diag)) {
                fit$mu2_diag[effectIdx, , , drop = FALSE]
            },
            coef = if (requireNamespace("mvsusieR", quietly = TRUE)) {
                mvsusieR::coef.mvsusie(fit)[-1, , drop = FALSE]
            },
            clfsr = if (!is.null(fit$conditional_lfsr)) {
                fit$conditional_lfsr[effectIdx, , , drop = FALSE]
            }
        ))
    )
}

#' Format Fine-mapping Post-processing for Protocol Output
#'
#' Promotes the primary method's per-method post-processing payload to the root
#' level and attaches the unified \code{top_loci} table. The primary method's
#' bare \code{FineMappingRow} appears at \code{$finemappingEntry}; wrap it
#' into a \code{FineMappingResult} collection at the pipeline level once (study,
#' context, trait, method) identity tags are known.
#'
#' @param post Output from \code{\link{postprocessFinemappingFits}}.
#' @param primaryMethod Method whose result should populate root-level fields.
#' @return A list with root-level fields including \code{finemappingEntry} (a
#'   bare \code{FineMappingRow} S4 payload) and \code{top_loci}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:40]
#' y <- eqtlRegionExample$yRes
#' fit <- susieR::susie(X, y, L = 5)
#' post <- postprocessFinemappingFits(fits = list(susie = fit),
#'   dataX = X, dataY = y)
#' formatFinemappingOutput(post, primaryMethod = "susie")
#' @export
formatFinemappingOutput <- function(post, primaryMethod) {
    methodPost <- post$finemappingResults[[primaryMethod]]
    if (is.null(methodPost)) {
        msg <- glue(
            "primaryMethod was not found in finemappingResults: ",
            "{primaryMethod}"
        )
        abort(msg)
    }
    c(
        methodPost,
        list(
            top_loci = post$top_loci
        )
    )
}

#' @noRd
getCsIndex <- function(snpsIdx, susieCs) {
    # Return ALL CS indices that contain this variant (not just one)
    idx <- which(map_lgl(susieCs, .csContains, snpsIdx = snpsIdx))
    if (length(idx) == 0) {
        return(NA_integer_)
    }
    return(idx)
}
#' @noRd
getTopVariantsIdx <- function(susieOutput, signalCutoff) {
    # `sets$cs` is absent when no credible set was found; list_c() is strict
    # about NULL where unlist() silently returned it.
    cs <- list_c(susieOutput$sets$cs %||% list())
    c(which(susieOutput$pip >= signalCutoff), cs) |>
        unique() |>
        sort()
}
# Returns a data.frame(variant_idx, cs_idx) with one row per (variant, CS) pair.
# Variants in multiple CSs get multiple rows.
#' @importFrom stringr str_replace
#' @noRd
getCsInfo <- function(susieOutputSetsCs, topVariantsIdx) {
    csNames <- names(susieOutputSetsCs)
    rows <- map(
        topVariantsIdx,
        .csInfoRow,
        susieOutputSetsCs = susieOutputSetsCs,
        csNames = csNames
    )
    bind_rows(rows)
}
#' @title Calculate Purity Measures for Credible Sets
#'
#' @description As an extension of the internal cal_purity function. This
#'   function computes purity metrics (minimum, mean, and median absolute
#'   correlations) for each credible set in a list of credible set indices,
#'   based on the provided X matrix. The output Purity depends on the method
#'   specified: for the 'min' method, it returns a single value for
#'   single-element sets or the minimum absolute correlation for others. For
#'   other methods, it returns a vector of three values (min, mean, median) for
#'   each set.
#'
#' @param lCs A list of credible set indices, where each element is a vector of
#'   indices corresponding to variables in a credible set.
#' @param X The data matrix used to compute correlations between variables in
#'   each credible set.
#' @param method A character string specifying the method to use for calculating
#'   purity. Defaults to 'min'. Other methods return a vector of min, mean, and
#'   median absolute correlations for each credible set.
#' @return A list where each element corresponds to a credible set and contains
#'   either a single purity value (for 'min' method and single-element sets) or
#'   a vector of purity metrics (for other methods and multi-element sets).
#' @noRd

# Absolute off-diagonal correlations among one credible set's variants. The
# diagonal is blanked so a set is never judged pure by its self-correlation.
# @noRd
.cpOffDiagonalLd <- function(csIndices, X) {
    x <- abs(computeLd(X[, csIndices, drop = FALSE], method = "sample"))
    replace(x, col(x) == row(x), NA)
}

# One credible set's purity: the weakest pairwise correlation, or the
# (min, mean, median) triple when more than the minimum is asked for. A
# single-variant set is pure by definition.
# @noRd
.cpSetPurity <- function(csIndices, X, method) {
    if (method == "min") {
        if (length(csIndices) == 1) {
            return(1)
        }
        return(min(.cpOffDiagonalLd(csIndices, X), na.rm = TRUE))
    }
    if (length(csIndices) == 1) {
        return(c(1, 1, 1))
    }
    x <- .cpOffDiagonalLd(csIndices, X)
    c(
        min(x, na.rm = TRUE),
        mean(x, na.rm = TRUE),
        median(x, na.rm = TRUE)
    )
}

calPurity <- function(lCs, X, method = "min") {
    # Each `lCs[[k]]` is documented (and always supplied) as a plain index
    # vector; the unlist() that used to sit here was a no-op that quietly
    # tolerated a nested list, leaving the contract unsettled.
    map(lCs, .cpSetPurity, X = X, method = method)
}


#' @importFrom checkmate assertNumber
#' @title Create Sets Similar to SuSiE Output from fSuSiE Object
#'
#' @description This function constructs a list that mimics the structure of
#'   SuSiE output sets from a fSuSiE object. It includes credible sets (cs) with
#'   their names, a purity dataframe, coverage information, and the requested
#'   coverage level.
#'
#' @param fsusieObj A fSuSiE object containing the results from a fSuSiE
#'   analysis. expected to at least have 'cs' and 'alpha' components.
#' @param requestedCoverage A numeric value specifying the desired coverage
#'   level for the credible sets. This is purely for record purpose so should be
#'   manually ensured that it correctly reflect the actual coverage used.
#'   Defaults to 0.95.
#' @param X Numeric genotype matrix used to compute credible-set purity.
#' @return A list containing named credible sets (cs), a dataframe of purity
#'   metrics (a \code{cs} label column plus minAbsCorr, meanAbsCorr,
#'   medianAbsCorr), an index of credible sets
#'   (cs_index), coverage values for each set, and the requested coverage level.
#'   Similar to the SuSiE set output
#' @examples
#' data(fsusieFineMappingExample)
#' fit <- getSusieFit(fsusieFineMappingExample)
#' fsusieGetCs(fit)
#' @export
fsusieGetCs <- function(fsusieObj, X, requestedCoverage = 0.95) {
    assertNumber(requestedCoverage, lower = 0, upper = 1)
    # Create 'cs' set with names
    csNamed <- set_names(
        fsusieObj$cs,
        str_c("L", seq_along(fsusieObj$cs))
    )

    # Create 'purity' data frame
    purity <- `colnames<-`(
        bind_rows(
            map(calPurity(fsusieObj$cs, X = X, method = "susie"), .asDataFrameT)
        ),
        c("minAbsCorr", "meanAbsCorr", "medianAbsCorr")
    )
    # Credible-set label as a `cs` column (was rownames; tibbles carry none).
    purityDf <- bind_cols(tibble(cs = names(csNamed)), purity)

    # Create 'coverage' without
    coverageVector <- map_dbl(
        seq_along(fsusieObj$alpha),
        .fsusieSetCoverage,
        fsusieObj = fsusieObj
    )

    # Combine all elements into a list
    sets <- list(
        cs = csNamed,
        purity = purityDf,
        cs_index = seq_along(fsusieObj$cs),
        coverage = coverageVector,
        requested_coverage = requestedCoverage
    )

    return(sets)
}

#' @title Wrapper for fsusie Function with Automatic Post-Processing
#'
#' @description This function serves as a wrapper for the fsusie function,
#'   facilitating automatic post-processing such as removing dummy credible sets
#'   (cs) that don't meet the minimum purity threshold and calculating
#'   correlations for the remaining cs. The function parameters are identical to
#'   those of the fSuSiE function.
#'
#' @param X Residual genotype matrix.
#' @param Y Response phenotype matrix.
#' @param pos Genomics position of phenotypes, used for specifying the wavelet
#'   model.
#' @param L The maximum number of the credible set.
#' @param prior method to generate the prior.
#' @param maxSnpEm maximum number of SNP used for learning the prior.
#' @param covLev Coverage level for the credible sets.
#' @param maxScale numeric, define the maximum of wavelet coefficients used in
#'   the analysis (2^maxScale). Set 10 true by default.
#' @param minPurity Minimum purity threshold for credible sets to be retained.
#' @param methodArgs Options forwarded to \code{fsusieR::susiF}, built with
#'   \code{\link{FsusieOptions}}. A bare list is refused: it cannot be checked
#'   against the engine, so a misspelled option would be silently ignored.
#' @return A modified fsusie object with the susie sets list, correlations for
#'   cs, alpha as df like susie, and without the dummy cs that do not meet the
#'   minimum purity requirement.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:50]
#' n <- nrow(X)
#' nPos <- 16
#' base <- sin(seq(0, 2 * pi, length.out = nPos))
#' Y <- matrix(rep(base, each = n), n, nPos) +
#'   X[, 1] %o% (0.5 * cos(seq(0, pi, length.out = nPos)))
#' pos <- seq_len(nPos)
#' fsusieWrapper(X, Y, pos = pos, L = 2,
#'   prior = "mixture_normal_per_scale", maxSnpEm = 10,
#'   covLev = 0.95, minPurity = 0.5, maxScale = 6)
#' @export
fsusieWrapper <- function(
    X,
    Y,
    pos,
    L,
    prior,
    maxSnpEm,
    covLev,
    minPurity,
    maxScale,
    methodArgs = FsusieOptions()
) {
    .assertMethodOptions(methodArgs, "FsusieOptions", "methodArgs")
    if (!requireNamespace("fsusieR", quietly = TRUE)) {
        abort("Package 'fsusieR' is required for this function.")
    }
    callArgs <- list_modify(
        list(
            X = X,
            Y = Y,
            pos = pos,
            L = L,
            prior = prior,
            max_SNP_EM = maxSnpEm,
            cov_lev = covLev,
            min_purity = minPurity,
            max_scale = maxScale
        ),
        !!!methodArgs
    )
    fsusieObj <- exec(fsusieR::susiF, !!!callArgs)
    .fsusieWrapperPostprocess(fsusieObj, X, minPurity, covLev)
}

# Drop dummy credible sets below the purity threshold (else build sets + CS
# correlations), then reshape alpha (per-effect list) into a single data.frame.
.fsusieWrapperPostprocess <- function(fsusieObj, X, minPurity, covLev) {
    withSets <- if (all(abs(as.numeric(fsusieObj$purity)) < minPurity)) {
        list_assign(
            fsusieObj,
            cs = list(NULL),
            sets = list(cs = list(NULL), requested_coverage = covLev)
        )
    } else {
        list_assign(
            fsusieObj,
            sets = fsusieGetCs(fsusieObj, X, requestedCoverage = covLev)
        )
    }
    list_assign(
        withSets,
        alpha = bind_rows(map(fsusieObj$alpha, .asDataFrameT))
    )
}


# =============================================================================
# Uniform fit wrappers for mvSuSiE (individual + RSS)
# -----------------------------------------------------------------------------
# Thin wrappers around mvsusieR::mvsusie and mvsusieR::mvsusie_rss. Every
# inline call across the package routes through these so the indirection
# is testable in one place and so future changes to the underlying mvsusieR
# API only need updating here.
# =============================================================================

#' Fit mvSuSiE on individual-level (X, Y) data
#'
#' Wrapper around \code{mvsusieR::mvsusie} with the canonical argument set used
#' inside fine-mapping and TWAS-weight pipelines.
#'
#' @param X Numeric matrix of genotypes (samples x variants).
#' @param Y Numeric matrix of multi-trait / multi-context outcomes (samples x
#'   conditions).
#' @param prior_variance Prior variance matrix; pass the output of
#'   \code{mvsusieR::create_mixture_prior(R = ncol(Y))} unless you have a
#'   domain-specific prior.
#' @param coverage Credible set coverage (default 0.95).
#' @param methodArgs Options forwarded to \code{mvsusieR::mvsusie}, built
#'   with \code{\link{MvsusieOptions}}. A bare list is refused: it cannot be
#'   checked against the engine, so a misspelled option would be silently
#'   ignored.
#' @return The fit object returned by \code{mvsusieR::mvsusie}.
#' @examples
#' \donttest{
#' # mvsusieR 0.3.0 calls susieR::block_coordinate_ascent unqualified, so
#' # susieR must be attached until mvsusieR adds the NAMESPACE import.
#' library(susieR)
#' data(multiTraitData)
#' X <- multiTraitData$X[, 1:60]
#' Y <- multiTraitData$Y
#' fitMvsusie(X, Y,
#'   prior_variance = mvsusieR::create_mixture_prior(R = ncol(Y)))
#' }
#' @importFrom checkmate assertNumber
#' @export
fitMvsusie <- function(
    X,
    Y,
    prior_variance,
    coverage = 0.95,
    methodArgs = MvsusieOptions()
) {
    .assertMethodOptions(methodArgs, "MvsusieOptions", "methodArgs")
    assertNumber(coverage, lower = 0, upper = 1)
    callArgs <- list_modify(
        list(
            X = X,
            Y = Y,
            prior_variance = prior_variance,
            coverage = coverage
        ),
        !!!methodArgs
    )
    exec(mvsusieR::mvsusie, !!!callArgs)
}

#' Fit mvSuSiE-RSS on summary-statistic (Z, R, N) data
#'
#' Wrapper around \code{mvsusieR::mvsusie_rss}. The underlying function was
#' renamed from \code{mvsusieRss} to \code{mvsusie_rss} upstream; this wrapper
#' insulates pecotmr from that naming.
#'
#' @param Z Numeric matrix of Z-scores (variants x conditions).
#' @param R Variant-by-variant LD correlation matrix.
#' @param N Scalar sample size (median across conditions when N varies).
#' @param prior_variance Prior variance matrix.
#' @param coverage Credible set coverage (default 0.95).
#' @param methodArgs Options forwarded to \code{mvsusieR::mvsusie_rss},
#'   built with \code{\link{MvsusieRssOptions}} --- the RSS entry point has
#'   its own constructor, since it is a different function from
#'   \code{mvsusie} with different arguments.
#' @return The fit object returned by \code{mvsusieR::mvsusie_rss}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ss <- lapply(seq_len(ncol(X)), function(j) {
#'   coef(summary(lm(y ~ X[, j])))[2, 1:2]
#' })
#' stat <- list(
#'   bhat = vapply(ss, `[`, numeric(1), 1L),
#'   shat = vapply(ss, `[`, numeric(1), 2L),
#'   z = vapply(ss, function(s) s[1] / s[2], numeric(1)),
#'   n = rep(nrow(X), ncol(X)))
#' LD <- cor(X)
#' fitMvsusieRss(Z = stat$z, R = LD, N = nrow(X), prior_variance = 1)
#' @importFrom checkmate assertNumber
#' @export
fitMvsusieRss <- function(
    Z,
    R,
    N,
    prior_variance,
    coverage = 0.95,
    methodArgs = MvsusieRssOptions()
) {
    .assertMethodOptions(methodArgs, "MvsusieRssOptions", "methodArgs")
    assertNumber(N, lower = 0, finite = TRUE)
    assertNumber(coverage, lower = 0, upper = 1)
    callArgs <- list_modify(
        list(
            Z = Z,
            R = R,
            N = N,
            prior_variance = prior_variance,
            coverage = coverage
        ),
        !!!methodArgs
    )
    exec(mvsusieR::mvsusie_rss, !!!callArgs)
}

#' Fit fSuSiE on individual-level (X, Y, pos) data
#'
#' Thin wrapper around \code{fsusieR::susiF}.
#'
#' @param X Numeric matrix of genotypes (samples x variants).
#' @param Y Numeric matrix of multi-trait outcomes (samples x traits).
#' @param pos Numeric vector of trait positions (length \code{ncol(Y)}).
#' @param methodArgs Options forwarded to \code{fsusieR::susiF}, built with
#'   \code{\link{FsusieOptions}}. A bare list is refused: it cannot be checked
#'   against the engine, so a misspelled option would be silently ignored.
#' @return The fit object returned by \code{fsusieR::susiF}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:50]
#' n <- nrow(X)
#' nPos <- 16
#' base <- sin(seq(0, 2 * pi, length.out = nPos))
#' Y <- matrix(rep(base, each = n), n, nPos) +
#'   X[, 1] %o% (0.5 * cos(seq(0, pi, length.out = nPos)))
#' pos <- seq_len(nPos)
#' fitFsusie(X, Y, pos = pos, methodArgs = FsusieOptions(L = 2))
#' @export
fitFsusie <- function(X, Y, pos, methodArgs = FsusieOptions()) {
    .assertMethodOptions(methodArgs, "FsusieOptions", "methodArgs")
    callArgs <- list_modify(list(X = X, Y = Y, pos = pos), !!!methodArgs)
    exec(fsusieR::susiF, !!!callArgs)
}

# =============================================================================
# SuSiE / mvSuSiE / fSuSiE TWAS weight extractors
# (relocated from regularizedRegressionWrappers.R: the SuSiE-family weight
#  extractors live with the rest of the fine-mapping/SuSiE wrappers).
# =============================================================================

# Shared helper for susie/susieAsh/susieInf weight extraction.
# @param fit A susie fit object (or NULL to fit from X, y).
# @param X Genotype matrix (optional).
# @param y Phenotype vector (optional).
# Drop-intercept TWAS coefficient weights from a fitted SuSiE(-RSS) model
# (zero the intercept, then coef.susie without the intercept row).
# @noRd
.susieCoefWeights <- function(fit) {
    coef.susie(list_assign(fit, intercept = 0))[-1]
}

# @param requiredFields Fields that must be present in the fit to extract
# weights.
# @param token SuSiE-family token ("susie" / "susieInf" / "susieAsh") selecting
#   the unmappable_effects mode. The fit is delegated to .fmFitSusieIndiv so the
#   package keeps a single susie-invocation point.
# @param userArgs Extra arguments forwarded to susieR::susie via
# .fmFitSusieIndiv.
#' @importFrom susieR coef.susie susie
#' @noRd
.susieExtractWeights <- function(
    fit,
    X,
    requiredFields,
    token = "susie",
    fitRetention = c("none", "slim", "full")
) {
    fitRetention <- arg_match(fitRetention)
    if (is.null(fit)) {
        msg <- glue(
            "{token}Weights: no '{token}' fit supplied. These extract weights ",
            "from an existing fit and never run fine-mapping themselves; run ",
            "fineMappingPipeline() with method '{token}' first and pass the ",
            "fit in."
        )
        abort(msg)
    }
    if (!is.null(X) && length(fit$pip) != ncol(X)) {
        nPip <- length(fit$pip)
        nX <- ncol(X)
        msg <- glue(
            "Dimension mismatch on number of variant in susie fit ",
            "{nPip} and TWAS weights {nX}. "
        )
        abort(msg)
    }
    weights <- if (all(is_in(requiredFields, names(fit)))) {
        .susieCoefWeights(fit)
    } else {
        rep(0, length(fit$pip))
    }
    if (identical(fitRetention, "none")) {
        return(weights)
    }
    `attr<-`(weights, "fit", fit)
}

#' Compute single-effect (SER) TWAS weights
#'
#' Extracts coefficients from an existing single-effect fit, as produced by
#' \code{fineMappingPipeline(..., methods = "ser")} on summary statistics.
#' The SER model carries the same \code{alpha} / \code{mu} /
#' \code{X_column_scale_factors} structure as a SuSiE fit --- with a single
#' effect, so \code{alpha} has one row --- so extraction is identical.
#'
#' SER fits only exist on the summary-statistics path, because
#' \code{susieR::susie_ser} is the only implementation
#' \code{\link{fineMappingPipeline}} has for the \code{"ser"} method. That
#' does not restrict the weights: a fit is a fit, and weights read off an
#' RSS fit are as usable for TWAS as any from penalized regression.
#'
#' @param X Optional genotype matrix; when supplied it is only used to
#'   check that the fit covers the same number of variants.
#' @param y Unused; retained for signature compatibility.
#' @param serFit Optional fitted SER object.
#' @param fitRetention How much of the fit is kept on the returned weights:
#'   \code{"none"} attaches nothing, \code{"slim"} and \code{"full"} attach
#'   it as the \code{"fit"} attribute. These extractors hold no trimmable
#'   intermediates, so the two retaining levels behave alike here.
#' @return Numeric vector of variant weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' z <- apply(X, 2, function(g) summary(stats::lm(y ~ g))$coefficients[2, 3])
#' fit <- suppressWarnings(susieR::susie_ser(z = z, n = nrow(X)))
#' serWeights(serFit = fit)
#' @export
serWeights <- function(
    X = NULL,
    y = NULL,
    serFit = NULL,
    fitRetention = c("none", "slim", "full")
) {
    fitRetention <- arg_match(fitRetention)
    .susieExtractWeights(
        serFit,
        X,
        requiredFields = c("alpha", "mu", "X_column_scale_factors"),
        token = "ser",
        fitRetention = fitRetention
    )
}

#' Compute SuSiE TWAS weights
#'
#' Extracts coefficients from an existing SuSiE fit.
#' from `X` and `y` before extracting weights.
#'
#' @param X Optional genotype matrix; when supplied it is only used to
#'   check that the fit covers the same number of variants.
#' @param y Unused; retained for signature compatibility.
#' @param susieFit Optional fitted SuSiE object.
#' @param fitRetention How much of the fit is kept on the returned weights:
#'   \code{"none"} attaches nothing, \code{"slim"} and \code{"full"} attach
#'   it as the \code{"fit"} attribute. These extractors hold no trimmable
#'   intermediates, so the two retaining levels behave alike here.
#' @return Numeric vector of variant weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' fit <- susieR::susie(X, y, L = 5)
#' susieWeights(susieFit = fit)
#' @importFrom checkmate assertFlag
#' @export
susieWeights <- function(
    X = NULL,
    y = NULL,
    susieFit = NULL,
    fitRetention = c("none", "slim", "full")
) {
    fitRetention <- arg_match(fitRetention)
    .susieExtractWeights(
        susieFit,
        X,
        requiredFields = c("alpha", "mu", "X_column_scale_factors"),
        token = "susie",
        fitRetention = fitRetention
    )
}

#' Compute SuSiE-ASH TWAS weights
#'
#' Extracts coefficients from an existing SuSiE-ASH fit or fits
#' `susieR::susie()` with `unmappable_effects = "ash"`.
#'
#' @param X Optional genotype matrix; when supplied it is only used to
#'   check that the fit covers the same number of variants.
#' @param y Unused; retained for signature compatibility.
#' @param susieAshFit Optional fitted SuSiE-ASH object.
#' @param fitRetention How much of the fit is kept on the returned weights:
#'   \code{"none"} attaches nothing, \code{"slim"} and \code{"full"} attach
#'   it as the \code{"fit"} attribute. These extractors hold no trimmable
#'   intermediates, so the two retaining levels behave alike here.
#' @return Numeric vector of variant weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' fit <- susieR::susie(X, y, L = 5)
#' susieAshWeights(susieAshFit = fit)
#' @importFrom checkmate assertFlag
#' @export
susieAshWeights <- function(
    X = NULL,
    y = NULL,
    susieAshFit = NULL,
    fitRetention = c("none", "slim", "full")
) {
    fitRetention <- arg_match(fitRetention)
    .susieExtractWeights(
        susieAshFit,
        X,
        requiredFields = c("alpha", "mu", "theta", "X_column_scale_factors"),
        token = "susieAsh",
        fitRetention = fitRetention
    )
}

#' Compute SuSiE-inf TWAS weights
#'
#' Extracts coefficients from an existing SuSiE-inf fit or fits
#' `susieR::susie()` with `unmappable_effects = "inf"`.
#'
#' @section Non-zero weights with zero PIPs: SuSiE-inf decomposes effects into a
#'   mappable component (driven by `alpha * mu`, reported as per-variant PIPs)
#'   and an infinitesimal component (driven by `theta`). When the fit converges
#'   with no mappable effects -- all `V` and `mu` zero, so every `pip == 0` --
#'   the returned weights are still non-zero because `susieR::coef.susie` adds
#'   `theta / X_column_scale_factors` to the mappable coefficient. This is
#'   intentional: it captures diffuse polygenic signal that the mappable
#'   component could not localize to any credible set. Consumers that interpret
#'   per-variant PIPs as a gate on whether to use the weights should be aware
#'   that low or zero PIPs do not imply zero TWAS weights here.
#'
#' @param X Optional genotype matrix; when supplied it is only used to
#'   check that the fit covers the same number of variants.
#' @param y Unused; retained for signature compatibility.
#' @param susieInfFit Optional fitted SuSiE-inf object.
#' @param fitRetention How much of the fit is kept on the returned weights:
#'   \code{"none"} attaches nothing, \code{"slim"} and \code{"full"} attach
#'   it as the \code{"fit"} attribute. These extractors hold no trimmable
#'   intermediates, so the two retaining levels behave alike here.
#' @return Numeric vector of variant weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' fit <- susieR::susie(X, y, L = 5)
#' susieInfWeights(susieInfFit = fit)
#' @importFrom checkmate assertFlag
#' @export
susieInfWeights <- function(
    X = NULL,
    y = NULL,
    susieInfFit = NULL,
    fitRetention = c("none", "slim", "full")
) {
    fitRetention <- arg_match(fitRetention)
    .susieExtractWeights(
        susieInfFit,
        X,
        requiredFields = c("alpha", "mu", "theta", "X_column_scale_factors"),
        token = "susieInf",
        fitRetention = fitRetention
    )
}
# Internal helper: extract weights from a susieRss fit.
# Mirrors .susie_extract_weights but uses the RSS interface.
#' @importFrom susieR coef.susie susie_rss
#' @noRd
.susieRssExtractWeights <- function(
    fit,
    R,
    requiredFields,
    token = "susie",
    fitRetention = c("none", "slim", "full")
) {
    fitRetention <- arg_match(fitRetention)
    if (is.null(fit)) {
        msg <- glue(
            "{token}RssWeights: no '{token}' fit supplied. These extract ",
            "weights from an existing fit and never run fine-mapping ",
            "themselves; run fineMappingPipeline() with method '{token}' ",
            "first and pass the fit in."
        )
        abort(msg)
    }
    if (length(fit$pip) != nrow(R)) {
        nPip <- length(fit$pip)
        nR <- nrow(R)
        msg <- glue(
            "Dimension mismatch: susieRss fit has {nPip} variants but R ",
            "has {nR} rows."
        )
        abort(msg)
    }
    weights <- if (all(is_in(requiredFields, names(fit)))) {
        .susieCoefWeights(fit)
    } else {
        rep(0, length(fit$pip))
    }
    if (identical(fitRetention, "none")) {
        return(weights)
    }
    `attr<-`(weights, "fit", fit)
}

#' Compute SuSiE-RSS TWAS weights
#'
#' Extracts coefficients from an existing SuSiE-RSS fit or fits
#' \code{susieR::susie_rss()} from summary statistics and LD.
#'
#' @param stat List with components \code{z} (z-scores), \code{n} (sample
#'   sizes).
#' @param LD LD correlation matrix.
#' @param susieRssFit A fitted SuSiE-RSS object. Required: these wrappers
#'   extract weights and never run fine-mapping themselves.
#' @param fitRetention How much of the fit is kept on the returned weights:
#'   \code{"none"} attaches nothing, \code{"slim"} and \code{"full"} attach
#'   it as the \code{"fit"} attribute. These extractors hold no trimmable
#'   intermediates, so the two retaining levels behave alike here.
#' @return Numeric vector of variant weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ss <- lapply(seq_len(ncol(X)), function(j) {
#'   coef(summary(lm(y ~ X[, j])))[2, 1:2]
#' })
#' stat <- list(
#'   bhat = vapply(ss, `[`, numeric(1), 1L),
#'   shat = vapply(ss, `[`, numeric(1), 2L),
#'   z = vapply(ss, function(s) s[1] / s[2], numeric(1)),
#'   n = rep(nrow(X), ncol(X)))
#' LD <- cor(X)
#' fit <- susieR::susie_rss(z = stat$z, R = LD, n = nrow(X), L = 5)
#' susieRssWeights(stat, LD, susieRssFit = fit)
#' @importFrom checkmate assertList assertFlag
#' @export
susieRssWeights <- function(
    stat,
    LD,
    susieRssFit = NULL,
    fitRetention = c("slim", "none", "full")
) {
    assertList(stat)
    fitRetention <- arg_match(fitRetention)
    .susieRssExtractWeights(
        fit = susieRssFit,
        R = LD,
        requiredFields = c("alpha", "mu", "X_column_scale_factors"),
        token = "susie",
        fitRetention = fitRetention
    )
}

#' Compute SuSiE-inf-RSS TWAS weights
#'
#' Extracts coefficients from an existing SuSiE-inf-RSS fit or fits
#' \code{susieR::susie_rss()} with \code{unmappable_effects = "inf"}.
#'
#' @inheritParams susieRssWeights
#' @param susieInfRssFit Optional pre-fitted SuSiE-inf-RSS object.
#' @return Numeric vector of variant weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ss <- lapply(seq_len(ncol(X)), function(j) {
#'   coef(summary(lm(y ~ X[, j])))[2, 1:2]
#' })
#' stat <- list(
#'   bhat = vapply(ss, `[`, numeric(1), 1L),
#'   shat = vapply(ss, `[`, numeric(1), 2L),
#'   z = vapply(ss, function(s) s[1] / s[2], numeric(1)),
#'   n = rep(nrow(X), ncol(X)))
#' LD <- cor(X)
#' fit <- susieR::susie_rss(z = stat$z, R = LD, n = nrow(X), L = 5)
#' susieInfRssWeights(stat, LD, susieInfRssFit = fit)
#' @importFrom checkmate assertList assertFlag
#' @export
susieInfRssWeights <- function(
    stat,
    LD,
    susieInfRssFit = NULL,
    fitRetention = c("slim", "none", "full")
) {
    assertList(stat)
    fitRetention <- arg_match(fitRetention)
    .susieRssExtractWeights(
        fit = susieInfRssFit,
        R = LD,
        requiredFields = c("alpha", "mu", "theta", "X_column_scale_factors"),
        token = "susieInf",
        fitRetention = fitRetention
    )
}

#' Compute SuSiE-ASH-RSS TWAS weights
#'
#' Extracts coefficients from an existing SuSiE-ASH-RSS fit or fits
#' \code{susieR::susie_rss()} with \code{unmappable_effects = "ash"}.
#'
#' @inheritParams susieRssWeights
#' @param susieAshRssFit Optional pre-fitted SuSiE-ASH-RSS object.
#' @return Numeric vector of variant weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ss <- lapply(seq_len(ncol(X)), function(j) {
#'   coef(summary(lm(y ~ X[, j])))[2, 1:2]
#' })
#' stat <- list(
#'   bhat = vapply(ss, `[`, numeric(1), 1L),
#'   shat = vapply(ss, `[`, numeric(1), 2L),
#'   z = vapply(ss, function(s) s[1] / s[2], numeric(1)),
#'   n = rep(nrow(X), ncol(X)))
#' LD <- cor(X)
#' fit <- susieR::susie_rss(z = stat$z, R = LD, n = nrow(X), L = 5)
#' susieAshRssWeights(stat, LD, susieAshRssFit = fit)
#' @importFrom checkmate assertList assertFlag
#' @export
susieAshRssWeights <- function(
    stat,
    LD,
    susieAshRssFit = NULL,
    fitRetention = c("slim", "none", "full")
) {
    assertList(stat)
    fitRetention <- arg_match(fitRetention)
    .susieRssExtractWeights(
        fit = susieAshRssFit,
        R = LD,
        requiredFields = c("alpha", "mu", "theta", "X_column_scale_factors"),
        token = "susieAsh",
        fitRetention = fitRetention
    )
}
#' Compute mvSuSiE TWAS weights
#'
#' Extracts coefficients from an existing mvSuSiE fit. This never fits
#' mvSuSiE itself: fine-mapping belongs to \code{fineMappingPipeline()}, and a
#' missing fit is an error rather than an invitation to refit.
#'
#' @param mvsusieFit A fitted mvSuSiE object. Required.
#' @return Matrix of variant weights.
#' @examples
#' \donttest{
#' # Requires susieR attached (mvsusieR 0.3.0 packaging limitation).
#' library(susieR)
#' data(multiTraitData)
#' X <- multiTraitData$X[, 1:60]
#' Y <- multiTraitData$Y
#' fit <- fitMvsusie(X = X, Y = Y,
#'   prior_variance = mvsusieR::create_mixture_prior(R = ncol(Y)),
#'   methodArgs = MvsusieOptions(L = 5))
#' mvsusieWeights(mvsusieFit = fit)
#' }
#' @export
mvsusieWeights <- function(mvsusieFit = NULL) {
    if (!requireNamespace("mvsusieR", quietly = TRUE)) {
        abort("Package 'mvsusieR' is required.")
    }
    if (is.null(mvsusieFit)) {
        msg <- glue(
            "mvsusieWeights: `mvsusieFit` is required. This extracts weights ",
            "from an existing mvSuSiE fit and never runs fine-mapping ",
            "itself; fit it via fineMappingPipeline() and pass the result in."
        )
        abort(msg)
    }
    mvsusieR::coef.mvsusie(mvsusieFit)[-1, , drop = FALSE]
}

# One wavelet basis row: inverse-DWT (wr) of the unit coefficient vector e_k,
# using the fit's template DWT object.
# @noRd
.fmReconstructUnit <- function(k, nWac, scaleCols, template) {
    coeffRow <- replace(numeric(nWac), k, 1)
    scaling <- replace(
        template$C,
        length(template$C),
        sum(coeffRow[scaleCols])
    )
    temp <- list_assign(template, D = coeffRow[-scaleCols], C = scaling)
    as.numeric(wavethresh::wr(temp))
}

# Build the wavelet synthesis (inverse-DWT) matrix S (n_wac x nFeat) for the
# basis fSuSiE uses, by reconstructing each unit wavelet coefficient through the
# SAME $D / $C assignment as out_prep.susiF (detail columns -> $D, the coarsest
# scaling column -> last $C entry), then `wavethresh::wr`. A wavelet-coefficient
# row `c` then maps to the feature domain as `c %*% S`. `scaleCols` is the
# column index of the scaling coefficient(s) (per the prior family). fSuSiE's
# default basis (DaubLeAsymm, filter 10) matches `wavethresh::wd`'s default, the
# same one out_prep uses, so the plain `wd(rep(0, nWac))` template is
# consistent.
# @noRd
.fsusieSynthesisMatrix <- function(nWac, scaleCols) {
    template <- wavethresh::wd(rep(0, nWac))
    rows <- map(seq_len(nWac), .fmReconstructUnit, nWac, scaleCols, template)
    exec(rbind, !!!rows)
}

#' Compute fSuSiE feature-level TWAS weights
#'
#' Collapses a functional SuSiE (\code{fsusieR::susiF}) fit back to a
#' \code{variants x features} weight matrix usable for TWAS prediction of each
#' molecular feature. fSuSiE fits the regression in the wavelet domain, storing
#' per-SNP posterior-mean wavelet effects \code{fitted_wc[[l]]}
#' (\code{nSNP x n_wac}) and inclusion probabilities \code{alpha[[l]]}. Because
#' the inverse wavelet transform \code{wr()} is linear, the posterior-mean
#' prediction pushes through to a per-SNP, per-feature weight matrix:
#' \deqn{W[j, f] = \sum_l alpha[[l]][j] \cdot
#'   \mathrm{wr}\!\left(fitted\_wc[[l]][j, ] / csd\_X[j]\right)[f].}
#' This is the exact analog of \code{coef.susie} for scalar SuSiE (all SNPs,
#' alpha-weighted), which spreads weight across the credible set -- more robust
#' for out-of-sample TWAS than fSuSiE's in-sample lead-SNP summary
#' (\code{update_cal_indf}).
#'
#' The reconstruction uses the raw posterior wavelet coefficients
#' \code{fitted_wc}, so it is independent of the \code{post_processing} mode
#' (\code{"smash"}/\code{"TI"}/\code{"HMM"}/\code{"none"}) -- that smoothing
#' only denoises the alpha-collapsed display curve \code{fitted_func}, never the
#' per-SNP predictive coefficients. The \code{$D}/\code{$C} coefficient layout
#' and wavelet basis mirror \code{out_prep.susiF}, so the feature-domain output
#' matches fSuSiE's own conventions.
#'
#' @param fsusieFit A fitted \code{fsusieR::susiF} object. Must retain
#'   \code{fitted_wc}, \code{alpha}, \code{csd_X}, \code{n_wac}, and
#'   \code{outing_grid} (i.e. an untrimmed fit). Required.
#' @param X,Y Accepted for call-compatibility with the multivariate
#'   weight-method dispatch in \code{\link{learnTwasWeights}}, which invokes
#'   every method as \code{fn(X = ., Y = ., ...)}. fSuSiE is a functional method
#'   that cannot be refit from a bare \code{(X, Y)} pair (it needs feature
#'   positions and the wavelet model), so these are ignored: a fitted
#'   \code{fsusieFit} is always required.
#' @param variantIds Optional character vector of variant IDs (length = number
#'   of SNPs in the fit) for the matrix row names. Defaults to
#'   \code{names(fsusieFit$csd_X)} / \code{names(fsusieFit$pip)}.
#' @param featureNames Optional character vector of feature (outcome) names for
#'   the matrix column names. Defaults to the fit's \code{outing_grid}.
#' @param fitRetention How much of the fit is kept on the returned weights:
#'   \code{"none"} attaches nothing, \code{"slim"} and \code{"full"} attach
#'   it as the \code{"fit"} attribute. These extractors hold no trimmable
#'   intermediates, so the two retaining levels behave alike here.
#' @return A numeric matrix of variant (rows) by feature (columns) weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:50]
#' n <- nrow(X)
#' nPos <- 16
#' base <- sin(seq(0, 2 * pi, length.out = nPos))
#' Y <- matrix(rep(base, each = n), n, nPos) +
#'   X[, 1] %o% (0.5 * cos(seq(0, pi, length.out = nPos)))
#' pos <- seq_len(nPos)
#' fit <- fsusieWrapper(X, Y, pos = pos, L = 2,
#'   prior = "mixture_normal_per_scale", maxSnpEm = 10,
#'   covLev = 0.95, minPurity = 0.5, maxScale = 6)
#' fsusieWeights(fsusieFit = fit, X = X, Y = Y,
#'   variantIds = colnames(X))
#' @export
fsusieWeights <- function(
    fsusieFit = NULL,
    X = NULL,
    Y = NULL,
    variantIds = NULL,
    featureNames = NULL,
    fitRetention = c("none", "slim", "full")
) {
    fitRetention <- arg_match(fitRetention)
    if (is.null(fsusieFit)) {
        msg <- glue(
            "fsusieWeights: `fsusieFit` is required. fSuSiE is functional ",
            "and cannot be refit from a bare (X, Y); fit it via ",
            "fineMappingPipeline() and pass the fitted fsusieR::susiF ",
            "object."
        )
        abort(msg)
    }
    fast <- .fsusieWeightsFastPath(fsusieFit, variantIds, fitRetention)
    if (!is.null(fast)) {
        return(fast)
    }
    .fsusieWeightsRequire()
    fit <- fsusieFit
    .fsusieWeightsCheckSlots(fit)
    csdX <- as.numeric(fit$csd_X)
    alphaList <- .fsusieAlphaList(fit$alpha)
    S <- .fsusieSynthesisMatrix(fit$n_wac, .fsusieScaleCols(fit))
    W <- .fsusieComputeW(fit, alphaList, csdX, S) |>
        .fsusieWeightsNames(
            fit,
            variantIds,
            featureNames,
            length(csdX),
            ncol(S)
        )
    if (identical(fitRetention, "none")) {
        return(W)
    }
    `attr<-`(W, "fit", fit)
}

# Row names only mean the variants when there is one per row of the matrix.
# @noRd
.withRownamesIfSized <- function(W, variantIds) {
    if (is.null(variantIds) || length(variantIds) != nrow(W)) {
        return(W)
    }
    `rownames<-`(W, variantIds)
}

# Fast path: a trimmed fit carries the precomputed variants x features weight
# matrix in `$coef` (fineMappingPipeline computes it eagerly while the full fit
# is in hand, since trimming drops fitted_wc/csd_X/...). NULL if not applicable.
.fsusieWeightsFastPath <- function(fsusieFit, variantIds, fitRetention) {
    if (!(is.matrix(fsusieFit$coef) && is.null(fsusieFit$fitted_wc))) {
        return(NULL)
    }
    W <- .withRownamesIfSized(fsusieFit$coef, variantIds)
    if (!identical(fitRetention, "none")) {
        return(`attr<-`(W, "fit", fsusieFit))
    }
    W
}

# fSuSiE weight reconstruction needs fsusieR + wavethresh.
.fsusieWeightsRequire <- function() {
    if (!requireNamespace("fsusieR", quietly = TRUE)) {
        abort("Package 'fsusieR' is required for fsusieWeights().")
    }
    if (!requireNamespace("wavethresh", quietly = TRUE)) {
        abort("Package 'wavethresh' is required for fsusieWeights().")
    }
}

# A full (untrimmed) fit must retain the wavelet-reconstruction slots.
.fsusieWeightsCheckSlots <- function(fit) {
    missingSlots <- setdiff(
        c("fitted_wc", "alpha", "csd_X", "n_wac", "outing_grid"),
        names(fit)
    )
    if (length(missingSlots) > 0L) {
        slotStr <- str_flatten(missingSlots, ", ")
        msg <- glue(
            "fsusieWeights: the fSuSiE fit is missing required slot(s): ",
            "{slotStr}. Pass an untrimmed fit (these are dropped when ",
            "trimmed)."
        )
        abort(msg)
    }
}

# Normalize alpha to a list of per-effect vectors. fsusieR::susiF returns a
# list; fsusieWrapper reshaping yields an L x nSNP matrix/data.frame.
.fsusieAlphaList <- function(alpha) {
    if (is.list(alpha) && !is.data.frame(alpha)) {
        return(map(alpha, as.numeric))
    }
    am <- as.matrix(alpha)
    map(seq_len(nrow(am)), .amRow, am = am)
}

# Scaling-coefficient column(s): coarsest level for a per-scale prior, else the
# last column (mirrors the two branches of out_prep.susiF).
.fsusieScaleCols <- function(fit) {
    perScale <- is_in(
        "mixture_normal_per_scale",
        class(fsusieR::get_G_prior(fit))
    )
    indxLst <- fsusieR::gen_wavelet_indx(log2(length(fit$outing_grid)))
    if (perScale) {
        indxLst[[length(indxLst)]]
    } else {
        ncol(as.matrix(fit$fitted_wc[[1L]]))
    }
}

# W = sum_l (alpha_l/csd_X-scaled fitted_wc_l) %*% S, one wavelet inverse
# transform (S) applied to every SNP/effect via a matrix multiply.
# @noRd
.fsusieSetCoverage <- function(i, fsusieObj) {
    sum(fsusieObj$alpha[[i]][fsusieObj$cs[[i]]])
}

# One effect's contribution to the weight matrix.
# @noRd
.fsusieEffectW <- function(l, fit, alphaList, invCsd, S) {
    (alphaList[[l]] * invCsd * as.matrix(fit$fitted_wc[[l]])) %*% S
}

.fsusieComputeW <- function(fit, alphaList, csdX, S) {
    invCsd <- 1 / csdX
    # W is the sum of the per-effect contributions, so it is a fold rather
    # than a matrix added into repeatedly.
    reduce(
        map(
            seq_along(fit$fitted_wc),
            .fsusieEffectW,
            fit = fit,
            alphaList = alphaList,
            invCsd = invCsd,
            S = S
        ),
        `+`,
        .init = matrix(0, nrow = length(csdX), ncol = ncol(S))
    )
}

# Attach variant (row) and feature/grid (column) names to the weight matrix.
.fsusieWeightsNames <- function(W, fit, variantIds, featureNames, p, nFeat) {
    rn <- variantIds %||% names(fit$csd_X) %||% names(fit$pip)
    cn <- if (
        is.null(featureNames) &&
            !is.null(fit$outing_grid) &&
            length(fit$outing_grid) == nFeat
    ) {
        as.character(fit$outing_grid)
    } else {
        featureNames
    }
    `dimnames<-`(
        W,
        list(
            if (length(rn) == p) rn else rownames(W),
            if (length(cn) == nFeat) cn else colnames(W)
        )
    )
}
#' Compute mvSuSiE-RSS TWAS weights from summary statistics
#'
#' Multi-context summary-statistics analog of \code{\link{mvsusieWeights}}:
#' extracts coefficients from an existing \code{mvsusieR::mvsusie_rss} fit.
#' It never runs fine-mapping itself --- \code{mvsusieRssFit} is required.
#'
#' Follows the \code{*_rss_weights(stat, LD, ...)} contract: \code{stat} and
#' \code{LD} describe the block the fit is being extracted against, and the
#' fit is checked against them rather than refitted.
#'
#' @param stat A list with \code{z} (matrix variants x conditions) and \code{n}
#'   (numeric vector or scalar), describing the block the fit came from.
#' @param LD LD correlation matrix for that block. The fit must cover the
#'   same number of variants; a mismatch is an error, since weights from one
#'   block are not comparable against another's LD.
#' @param mvsusieRssFit A fitted \code{mvsusieRss} object. Required: this
#'   extracts weights and never runs fine-mapping itself.
#' @param fitRetention How much of the fit is kept on the returned weights:
#'   \code{"none"} attaches nothing, \code{"slim"} and \code{"full"} attach
#'   it as the \code{"fit"} attribute. These extractors hold no trimmable
#'   intermediates, so the two retaining levels behave alike here.
#'
#' @return A numeric matrix of per-variant per-context weights (variants x
#'   conditions).
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' ss <- lapply(seq_len(ncol(X)), function(j) {
#'   coef(summary(lm(y ~ X[, j])))[2, 1:2]
#' })
#' stat <- list(
#'   bhat = vapply(ss, `[`, numeric(1), 1L),
#'   shat = vapply(ss, `[`, numeric(1), 2L),
#'   z = vapply(ss, function(s) s[1] / s[2], numeric(1)),
#'   n = rep(nrow(X), ncol(X)))
#' LD <- cor(X)
#' fit <- fitMvsusieRss(Z = stat$z, R = LD, N = nrow(X),
#'   prior_variance = 1)
#' mvsusieRssWeights(stat, LD, mvsusieRssFit = fit)
#' @export
mvsusieRssWeights <- function(
    stat,
    LD,
    mvsusieRssFit = NULL,
    fitRetention = c("none", "slim", "full")
) {
    assertList(stat)
    fitRetention <- arg_match(fitRetention)
    if (!requireNamespace("mvsusieR", quietly = TRUE)) {
        abort("Package 'mvsusieR' is required.")
    }
    if (is.null(mvsusieRssFit)) {
        msg <- glue(
            "mvsusieRssWeights: `mvsusieRssFit` is required. This extracts ",
            "weights from an existing mvSuSiE-RSS fit and never runs ",
            "fine-mapping itself; fit it via fineMappingPipeline() and pass ",
            "the result in."
        )
        abort(msg)
    }
    weights <- mvsusieR::coef.mvsusie(mvsusieRssFit)[-1, , drop = FALSE]
    # `stat` / `LD` describe the block the fit is being combined with. The
    # susie *RssWeights wrappers check the fit against them; this one used
    # to take both and read neither, so a fit from a different block was
    # extracted without complaint.
    .mvsusieRssAssertBlock(weights, LD)
    if (identical(fitRetention, "none")) {
        return(weights)
    }
    `attr<-`(weights, "fit", mvsusieRssFit)
}

# The supplied fit must describe the LD block it is being extracted against;
# the pairing is what makes the weights comparable across a region.
# @noRd
.mvsusieRssAssertBlock <- function(weights, LD) {
    if (is.null(LD) || nrow(weights) == nrow(LD)) {
        return(invisible(NULL))
    }
    abort(glue(
        "Dimension mismatch: mvsusieRss fit has {nrow(weights)} variants ",
        "but LD has {nrow(LD)} rows."
    ))
}


# =============================================================================
# Cross-condition credible-set merging
# =============================================================================

# Identify variant IDs that are associated with more than one credible set.
# @noRd
# The sets a variant belongs to, but only when it belongs to more than one.
# @noRd
.ovlMultiSets <- function(entry) {
    sets <- entry[["sets"]]
    if (length(sets) > 1) sets else NULL
}

.identifyOverlapSets <- function(variantsSetsAndPipsList) {
    compact(map(variantsSetsAndPipsList, .ovlMultiSets))
}

# Union-find root of `x` following the `parent` map.
# @noRd
.ufFindRoot <- function(x, parent) {
    # Walking to the root recursively rather than reassigning `x`; union-find
    # trees are shallow, so the depth is not a concern.
    if (identical(parent[[x]], x)) {
        return(x)
    }
    .ufFindRoot(parent[[x]], parent)
}

# Union-find merge of `a` and `b` in `parent`; returns the updated parent map.
# @noRd
.ufUnion <- function(a, b, parent) {
    rootA <- .ufFindRoot(a, parent)
    rootB <- .ufFindRoot(b, parent)
    if (identical(rootA, rootB)) {
        return(parent)
    }
    `[[<-`(parent, rootB, value = rootA)
}

# Merge overlapping credible sets using connected components (union-find).
# @noRd
# @noRd
.ufPair <- function(s, first) {
    list(a = first, b = s)
}

# One overlap's (first set, other set) pairs.
# @noRd
.ufOverlapPairs <- function(sets) {
    if (length(sets) <= 1) {
        return(list())
    }
    map(sets[-1], .ufPair, first = sets[[1]])
}

# @noRd
.ufUnionPair <- function(parent, pair) {
    .ufUnion(pair$a, pair$b, parent)
}

# Every member of one component mapped to that component's joint label.
# @noRd
.ufComponentLabels <- function(members) {
    label <- str_flatten(sort(members), ",")
    set_names(as.list(rep(label, length(members))), members)
}

.mergeAndUpdateOverlapSets <- function(variantsSetsAndPipsList, overlapSets) {
    allSets <- unique(list_c(overlapSets))
    if (length(allSets) == 0) {
        return(list())
    }

    # Each overlap ties its sets to the first one; the merges are a fold over
    # those pairs, since every union sees the map the previous one produced.
    parent <- reduce(
        .fmwConcat(map(overlapSets, .ufOverlapPairs)),
        .ufUnionPair,
        .init = set_names(allSets, allSets)
    )
    components <- split(
        names(parent),
        map_chr(names(parent), .ufFindRoot, parent)
    )
    setNameMap <- .fmwConcat(map(components, .ufComponentLabels))

    # Update each variant's credible set names
    updatedCredibleSets <- map(
        set_names(
            names(variantsSetsAndPipsList),
            names(variantsSetsAndPipsList)
        ),
        .updateCredibleSet,
        variantsSetsAndPipsList = variantsSetsAndPipsList,
        setNameMap = setNameMap
    )
    return(updatedCredibleSets)
}

# Collapse the per-variant extracted-CS map into a top-loci data frame: merge
# overlapping credible sets, then one row per variant with its merged CS label,
# max PIP and median PIP.
# @noRd
.combineTopLoci <- function(extractedResult) {
    if (length(extractedResult) == 0) {
        return(NULL)
    }

    overlapSets <- .identifyOverlapSets(extractedResult)
    hasOverlaps <- length(overlapSets) != 0
    mergedSets <- if (hasOverlaps) {
        .mergeAndUpdateOverlapSets(extractedResult, overlapSets = overlapSets)
    } else {
        NULL
    }

    topLociDf <- bind_rows(map(
        names(extractedResult),
        .csMergedVariantRow,
        extractedResult = extractedResult,
        hasOverlaps = hasOverlaps,
        mergedSets = mergedSets
    ))
    return(topLociDf)
}

# Build the per-variant extracted-CS map from a fine-mapping result: for each
# entry, one record per (variant, credible set) labelled cs_<entry>_<set>,
# aggregated by variant preserving first-seen order.
# @noRd
.fmExtractTopLoci <- function(fineMappingResult, csCol) {
    entries <- .collectionEntries(fineMappingResult)
    rows <- map_dfr(
        seq_along(entries),
        .extractCsEntryRows,
        entries = entries,
        csCol = csCol
    )

    if (is.null(rows) || nrow(rows) == 0) {
        return(list())
    }

    # Aggregate by variant_id preserving first-seen order.
    seenOrder <- unique(rows$variant_id)
    splitRows <- split(rows, factor(rows$variant_id, levels = seenOrder))
    map(splitRows, .csSplitToList)
}

#' Merge SuSiE credible sets across conditions
#'
#' Reconciles per-condition (univariate) SuSiE fine-mapping into a single set of
#' merged credible sets. Each row of the supplied
#' \code{\link{QtlFineMappingResult}} is treated as one condition (its
#' \code{topLoci} carrying that condition's credible sets); credible sets that
#' share variants across conditions are unioned via connected components, and
#' every variant is reported with its merged credible-set label plus the maximum
#' and median PIP across the conditions it appears in. A typical use is
#' selecting a representative lead variant per merged credible set to assemble
#' the \code{"strong"} input for \code{\link{mashPipeline}}.
#'
#' @param fineMappingResult A \code{\link{QtlFineMappingResult}} (or any
#'   \code{FineMappingResult}) produced by per-condition SuSiE fine-mapping.
#'   Each entry's \code{topLoci} must carry a credible-set column
#'   (\code{cs_<coverage*100>}, e.g. \code{cs_95}, with values such as
#'   \code{"susie_1"} where the trailing integer is the set index and \code{_0}
#'   means "not in a credible set") and a PIP column.
#' @param coverage Credible-set coverage level selecting the \code{cs_*} column
#'   (default \code{0.95} -> \code{cs_95}).
#' @return A \code{data.frame} with one row per variant: \code{variant_id},
#'   \code{credibleSetNames} (the merged credible-set label), \code{maxPip} and
#'   \code{medianPip}; or \code{NULL} when no credible sets are present.
#' @seealso \code{\link{fineMappingPipeline}}, \code{\link{mashPipeline}}
#' @importFrom purrr map_dfr
#' @importFrom stats median
#' @examples
#' data(qtlFineMappingExample)
#' mergeSusieCs(fineMappingResult = qtlFineMappingExample)
#' @export
mergeSusieCs <- function(fineMappingResult, coverage = 0.95) {
    if (!is(fineMappingResult, "FineMappingResultBase")) {
        msg <- glue(
            "`fineMappingResult` must be a QtlFineMappingResult (or ",
            "FineMappingResult)."
        )
        abort(msg)
    }
    csCol <- str_c("cs_", as.integer(round(coverage * 100)))

    # Each row (entry) of the fine-mapping result is one condition. Build a flat
    # data frame of (variant_id, pip, set_name) across conditions, giving each
    # condition's credible sets a unique "cs_<conditionIdx>_<setIdx>" label.
    extractedTopLoci <- .fmExtractTopLoci(fineMappingResult, csCol)
    if (length(extractedTopLoci) == 0) {
        return(NULL)
    }
    combinedTopLociDf <- .combineTopLoci(extractedTopLoci)
    if (is.null(combinedTopLociDf) || nrow(combinedTopLociDf) == 0) {
        return(NULL)
    }
    combinedTopLociDf <- distinct(
        combinedTopLociDf,
        .data$variant_id,
        .keep_all = TRUE
    )
    return(combinedTopLociDf)
}


# =============================================================================
# Post-fine-mapping credible-set extraction
# -----------------------------------------------------------------------------
# Pull the trimmed SuSiE fit out of a pipeline result and reduce it to the
# per-credible-set / top-PIP diagnostic rows that summaryStatsQc() and the
# fine-mapping report consume. Relocated here from sumstatsQc.R: these read a
# fine-mapping result, they do not perform QC.
# =============================================================================

#' Extract the trimmed SuSiE fit from a finemapping pipeline result
#'
#' Returns the trimmed model fit underlying \code{con_data$finemappingEntry} (a
#' \code{FineMappingRow} S4 object), or NULL if no fine-mapping entry is
#' attached.
#'
#' @param conData List. The method-layer entry from a finemapping pipeline
#'   result, expected to carry \code{$finemappingEntry} as a
#'   \code{FineMappingRow} object.
#' @return The trimmed fit (a list with \code{pip}, \code{sets}, etc.) or NULL.
#' @examples
#' data(qtlSumStatsExample)
#' getSusieResult(qtlSumStatsExample)
#' @importFrom checkmate assertList
#' @export
getSusieResult <- function(conData) {
    # No type guard: this is duck-typed on `$` and `length()` and returns NULL
    # for anything without a `finemappingEntry`. Its own @example passes a
    # QtlSumStats, which assertList rejects.
    if (length(conData) == 0) {
        return(NULL)
    }
    fm <- conData$finemappingEntry
    if (is.null(fm) || !is(fm, "FineMappingResultBase")) {
        return(NULL)
    }
    trimmed <- .fmrPartsSusieFit(fm)
    if (length(trimmed) == 0) {
        return(NULL)
    }
    trimmed
}

#' Process Credible Sets (CS) from Finemapping Results
#'
#' This function extracts and processes information for each Credible Set (CS)
#' from finemapping results, typically obtained from a finemapping RDS file.
#'
#' @param fmRow A \code{\link{fineMappingRow}}, or a single-row
#'   fine-mapping collection as returned by
#'   \code{\link{getFineMappingResult}}, carrying the SuSiE fit and
#'   variant ids.
#' @param csNames Character vector. Names of the Credible Sets, usually in the
#'   format "L_<number>".
#' @param topLociTable Data frame. The top-loci table (e.g. from
#'   \code{\link{getTopLoci}}) carrying \code{variant_id}, \code{pip}, and
#'   \code{z} columns.
#' @param ldSource The LD source from which the between-credible-set correlation
#'   is derived on demand: a \code{QtlDataset} (individual-level) or a
#'   \code{QtlSumStats} / \code{GwasSumStats} (summary statistics). See
#'   \code{\link{computeCsCorrelation}}.
#'
#' @return A data frame with one row per CS, containing the following columns:
#'   \item{cs_name}{Name of the Credible Set}
#'   \item{variants_per_cs}{Number of variants in the CS}
#'   \item{top_variant}{ID of the variant with the highest PIP in the CS}
#'   \item{top_variant_index}{Global index of the top variant}
#'   \item{top_pip}{Highest Posterior Inclusion Probability (PIP) in the CS}
#'   \item{top_z}{Z-score of the top variant}
#'   \item{p_value}{P-value calculated from the top Z-score}
#'   \item{cs_corr_1, cs_corr_2, ...}{Each CS's pairwise correlation with every
#'     CS (its row of the between-CS matrix, self-correlation on the diagonal),
#'     computed on demand from \code{ldSource}. Absent when there are fewer than
#'     two credible sets.}
#'   \item{cs_corr_max}{Maximum absolute between-CS correlation (excluding the
#'     self == 1); \code{NA} for a single CS.}
#'   \item{cs_corr_min}{Minimum absolute between-CS correlation; \code{NA} for a
#'     single CS.}
#'
#' @details This function is designed to be used only when there is at least one
#'   Credible Set in the finemapping results usually for a given study and
#'   block. It processes each CS, extracting key information such as the top
#'   variant, its statistics, and correlation information between multiple CS if
#'   available.
#'
#' @importFrom purrr map map_dbl map_int
#' @importFrom dplyr bind_rows
#'
#' @examples
#' data(qtlSumStatsExample)
#' vids <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")
#' fit <- list(pip = c(0.1, 0.7, 0.2), sets = list(cs = list(L_1 = c(1, 2))))
#' tl <- data.frame(variant_id = vids, pip = c(0.1, 0.7, 0.2),
#'   z = c(1.0, 3.5, -0.5))
#' fe <- fineMappingRow(variantIds = vids, susieFit = fit, topLoci = tl)
#' # A single credible set has no between-CS correlation (cs_corr_* are NA), so
#' # the ldSource is not consulted here.
#' extractCsInfo(fe, csNames = "L_1", topLociTable = tl,
#'   ldSource = qtlSumStatsExample)
#'
#' @importFrom checkmate assertClass assertCharacter
#' @export
extractCsInfo <- function(fmRow, csNames, topLociTable, ldSource) {
    assertCharacter(csNames, any.missing = FALSE)
    fm <- fmRow
    trimmed <- .fmrPartsSusieFit(fm)
    variantNames <- .fmrPartsVariantIds(fm)
    csCorr <- .rowCsCorrelation(fm, ldSource)
    rows <- map(
        seq_along(csNames),
        .extractCsInfoRow,
        csNames = csNames,
        trimmed = trimmed,
        variantNames = variantNames,
        topLociTable = topLociTable
    )
    .csAppendCorrelationCols(bind_rows(rows), csCorr)
}

#' Extract Information for Top Variant from Finemapping Results
#'
#' This function extracts information about the variant with the highest
#' Posterior Inclusion Probability (PIP) from finemapping results, typically
#' used when no Credible Sets (CS) are identified in the analysis.
#'
#' @param fmRow A \code{\link{fineMappingRow}}, or a single-row
#'   fine-mapping collection as returned by
#'   \code{\link{getFineMappingResult}}, carrying the SuSiE fit and
#'   variant ids.
#' @param sumstats A list or data frame carrying a \code{z} element aligned to
#'   the fit's variants (\code{sumstats$z}).
#'
#' @return A data frame with one row containing the following columns:
#'   \item{cs_name}{NA (as no CS is identified)}
#'   \item{variants_per_cs}{NA (as no CS is identified)}
#'   \item{top_variant}{ID of the variant with the highest PIP}
#'   \item{top_variant_index}{Index of the top variant in the original data}
#'   \item{top_pip}{Highest Posterior Inclusion Probability (PIP)}
#'   \item{top_z}{Z-score of the top variant}
#'   \item{p_value}{P-value calculated from the top Z-score}
#'   \item{cs_corr_max}{NA (no between-CS correlation without a CS)}
#'   \item{cs_corr_min}{NA (no between-CS correlation without a CS)}
#'
#' @details This function is designed to be used when no Credible Sets are
#'   identified in the finemapping results, but information about the most
#'   significant variant is still desired. It identifies the variant with the
#'   highest PIP and extracts relevant statistical information.
#'
#' @note This function is particularly useful for capturing information about
#'   potentially important variants that might be included in Credible Sets
#'   under different analysis parameters or lower coverage. It maintains a
#'   structure similar to the output of `extract_cs_info()` for consistency in
#'   downstream analyses.
#'
#' @seealso \code{\link{extractCsInfo}} for processing when Credible Sets are
#'   present.
#'
#' @examples
#' vids <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")
#' fit <- list(pip = c(0.1, 0.7, 0.2))
#' tl <- data.frame(variant_id = vids, pip = c(0.1, 0.7, 0.2))
#' fe <- fineMappingRow(variantIds = vids, susieFit = fit, topLoci = tl)
#' extractTopPipInfo(fe, sumstats = list(z = c(1.0, 3.5, -0.5)))
#'
#' @importFrom checkmate assertClass
#' @export
extractTopPipInfo <- function(fmRow, sumstats) {
    fm <- fmRow
    trimmed <- .fmrPartsSusieFit(fm)
    variantNames <- .fmrPartsVariantIds(fm)
    # Find the variant with the highest PIP
    topPipIndex <- which.max(trimmed$pip)
    topPip <- trimmed$pip[topPipIndex]
    topVariant <- variantNames[topPipIndex]
    topZ <- sumstats$z[topPipIndex]
    pValue <- .zToPvalue(topZ)

    list(
        cs_name = NA,
        variants_per_cs = NA,
        top_variant = topVariant,
        top_variant_index = topPipIndex,
        top_pip = topPip,
        top_z = topZ,
        p_value = pValue,
        cs_corr_max = NA_real_,
        cs_corr_min = NA_real_
    )
}

# Reduce one credible set's correlation vector (a row of the between-CS matrix,
# the self-correlation == 1 on the diagonal) to its |corr| max/min, excluding
# every self / perfect correlation (== 1). An empty result yields NA.
# @noRd
.extractCorrelations <- function(x) {
    filtered <- abs(x[x != 1])
    if (length(filtered) == 0L) {
        return(list(max_corr = NA_real_, min_corr = NA_real_))
    }
    list(
        max_corr = max(filtered, na.rm = TRUE),
        min_corr = min(filtered, na.rm = TRUE)
    )
}

# Append the between-CS correlation columns to the per-CS summary `base` from
# the m x m matrix `csCorr` (whose rows are aligned to `base`): the expanded
# cs_corr_1..m (each CS's row, self-correlation on the diagonal) plus
# cs_corr_max / cs_corr_min (|corr| excluding the self == 1). A NULL matrix
# (fewer than two credible sets) yields NA max/min and no expanded columns.
# @noRd
.csAppendCorrelationCols <- function(base, csCorr) {
    if (is.null(csCorr)) {
        return(mutate(base, cs_corr_max = NA_real_, cs_corr_min = NA_real_))
    }
    perRow <- apply(csCorr, 1L, .extractCorrelations, simplify = FALSE)
    expanded <- `names<-`(
        as_tibble(csCorr, .name_repair = "minimal"),
        str_c("cs_corr_", seq_len(ncol(csCorr)))
    )
    # unname(): apply() names its result by the matrix rownames, which map_dbl
    # then carries into the column (tibbles preserve element names).
    base |>
        bind_cols(expanded) |>
        mutate(
            cs_corr_max = unname(map_dbl(perRow, "max_corr")),
            cs_corr_min = unname(map_dbl(perRow, "min_corr"))
        )
}

# One credible set's scalar summary row (top variant / PIP / z / p). The
# between-CS correlation columns are appended once by .csAppendCorrelationCols.
# @noRd
.extractCsInfoRow <- function(i, csNames, trimmed, variantNames, topLociTable) {
    csName <- csNames[i]
    indices <- trimmed$sets$cs[[csName]]
    csVariants <- variantNames[indices]
    csData <- filter(topLociTable, is_in(.data$variant_id, csVariants))
    topRow <- which.max(csData$pip)
    topVariant <- csData$variant_id[topRow]
    topZ <- csData$z[topRow]
    tibble(
        cs_name = csName,
        variants_per_cs = length(csVariants),
        top_variant = topVariant,
        top_variant_index = which(variantNames == topVariant),
        top_pip = csData$pip[topRow],
        top_z = topZ,
        p_value = .zToPvalue(topZ)
    )
}


# =============================================================================
# SuSiE-family fitters (single-fit wrappers + per-block dispatch)
# -----------------------------------------------------------------------------
# Relocated here from fineMappingPipeline.R so all method-fitting wrappers live
# in one file, alongside fitMvsusie / fitFsusie / fitSusieInfThenSusie.
# `.fmFitSusie{Indiv,Rss,Ser}` each invoke a single susieR entry point;
# `.fmFit{X,Rss}Block` fit every requested token on one (X, y) / (z, R, n)
# block. They call orchestration helpers that remain in fineMappingPipeline.R
# (.fmResolveSusieChain / .fmPostprocessOne / .fmMergeUserArgs /
# .fineMappingMethodCapabilities); all resolve within the package namespace.
# =============================================================================

# Fit one of the SuSiE-family individual-level methods on (X, y). When
# `chainFromInf` is non-NULL, the susieInf fit it points at is used as
# initialisation (with prepareSusieFromInfArgs); otherwise a plain fit
# with the requested `unmappable_effects` is performed. `userArgs` are
# spliced via .fmMergeUserArgs (user wins over chain/base/capability
# defaults), so the caller can override things like L, max_iter,
# estimate_residual_method, refine, etc.
# @noRd
.fmFitSusieIndiv <- function(
    X,
    y,
    token,
    chainFromInf = NULL,
    coverage = 0.95,
    userArgs = NULL
) {
    info <- .fineMappingMethodCapabilities[[token]]
    if (is.null(info) || identical(info$unmappableEffects, NA_character_)) {
        msg <- glue(
            ".fmFitSusieIndiv: token '{token}' is not a SuSiE-family method."
        )
        abort(msg)
    }
    baseArgs <- list(
        X = X,
        y = y,
        coverage = coverage,
        unmappable_effects = info$unmappableEffects
    )
    fitArgs <- if (!is.null(chainFromInf) && token != "susieInf") {
        # SuSiE(-ash) initialised from a SuSiE-inf fit. userArgs are folded
        # into the arg prep (not merged afterwards) so L_greedy is clamped to
        # min(#inf effects, L) rather than passed through raw.
        chainedArgs <- prepareSusieFromInfArgs(
            .fmMergeUserArgs(list(), token, userArgs),
            chainFromInf,
            refineDefault = if (token == "susie") TRUE else NULL,
            unmappableEffects = if (token == "susieAsh") "ash" else "none"
        )
        list_assign(
            list_modify(baseArgs, !!!compact(chainedArgs)),
            X = X,
            y = y,
            coverage = coverage
        )
    } else {
        .fmMergeUserArgs(
            list_assign(baseArgs, !!!.fmSusieTokenDefaults(token)),
            token,
            userArgs
        )
    }
    .setFinemappingFitClass(exec(susieR::susie, !!!fitArgs), token)
}


# Sumstat counterpart of .fmFitSusieIndiv. Calls susieR::susie_rss with
# the same unmappable_effects switch, chained init, and userArgs merge.
# @noRd
.fmFitSusieRss <- function(
    z,
    R,
    n,
    token,
    chainFromInf = NULL,
    coverage = 0.95,
    userArgs = NULL,
    rFinite = NULL,
    rMismatch = "none",
    rssControl = NULL
) {
    info <- .fmRssValidateToken(token)
    baseArgs <- .fmRssAddControl(
        c(
            list(
                z = z,
                R = R,
                n = n,
                coverage = coverage,
                unmappable_effects = info$unmappableEffects
            ),
            # rFinite = NULL omits the element -> susie_rss default; these sit
            # in baseArgs so they survive the chained modifyList / non-chained
            # userArgs merge, while user methodArgs (folded in after) still
            # override them.
            compact(list(R_finite = rFinite, R_mismatch = rMismatch))
        ),
        rssControl
    )
    fitArgs <- if (!is.null(chainFromInf) && token != "susieInf") {
        .fmRssChainedArgs(
            baseArgs,
            token,
            userArgs,
            chainFromInf,
            z,
            R,
            n,
            coverage
        )
    } else {
        .fmRssNonChainedArgs(baseArgs, token, userArgs)
    }
    # All susie_rss fits get the "susieRss" S3 class for post-processing (drives
    # the Xcorr cs-input mode); token distinction stays in the `method` column.
    .setFinemappingFitClass(exec(susieR::susie_rss, !!!fitArgs), "susieRss")
}

# Validate the method token and return its capability record.
.fmRssValidateToken <- function(token) {
    info <- .fineMappingMethodCapabilities[[token]]
    if (is.null(info) || identical(info$unmappableEffects, NA_character_)) {
        msg <- glue(
            ".fmFitSusieRss: token '{token}' is not a SuSiE-family method."
        )
        abort(msg)
    }
    info
}

# Optional susie_rss_control() settings, supplied as a named list and forwarded
# as susie_rss()'s `control` argument.
.fmRssAddControl <- function(baseArgs, rssControl) {
    if (is.null(rssControl)) {
        return(baseArgs)
    }
    if (
        !is.list(rssControl) ||
            is.null(names(rssControl)) ||
            any(str_length(names(rssControl)) == 0L)
    ) {
        msg <- glue(
            ".fmFitSusieRss: `rssControl` must be a named list of ",
            "susieR::susie_rss_control() settings."
        )
        abort(msg)
    }
    list_assign(
        baseArgs,
        control = exec(susieR::susie_rss_control, !!!rssControl)
    )
}

# SuSiE-RSS(-ash) initialised from a SuSiE-inf fit; userArgs folded into the arg
# prep so L_greedy is clamped rather than passed through raw.
.fmRssChainedArgs <- function(
    baseArgs,
    token,
    userArgs,
    chainFromInf,
    z,
    R,
    n,
    coverage
) {
    chainedArgs <- prepareSusieFromInfArgs(
        .fmMergeUserArgs(list(), token, userArgs),
        chainFromInf,
        refineDefault = if (token == "susie") TRUE else NULL,
        unmappableEffects = if (token == "susieAsh") "ash" else "none"
    )
    list_assign(
        list_modify(baseArgs, !!!compact(chainedArgs)),
        z = z,
        R = R,
        n = n,
        coverage = coverage
    )
}

# Non-chained fit: token-specific defaults then the user methodArgs merge.
.fmRssNonChainedArgs <- function(baseArgs, token, userArgs) {
    .fmMergeUserArgs(
        list_assign(baseArgs, !!!.fmSusieTokenDefaults(token)),
        token,
        userArgs
    )
}

# Token-specific susie defaults. `model_init` is deliberately not among them:
# baseArgs never carries one, so susieR's own default already applies.
# @noRd
.fmSusieTokenDefaults <- function(token) {
    if (token == "susieInf") {
        return(list(convergence_method = "pip", refine = FALSE))
    }
    if (token == "susieAsh") {
        return(list(convergence_method = "pip"))
    }
    list()
}

# Single-effect (SER) sumstat fit via susieR::susie_ser on z + n. LD-free (no R,
# no L, no unmappable_effects), so it cannot reuse .fmFitSusieRss. Tagged
# "susieRss" so the shared post-processing (credible sets + purity against the
# LD sketch) applies unchanged.
# @noRd
.fmFitSusieSer <- function(z, n, coverage = 0.95, userArgs = NULL) {
    baseArgs <- .fmMergeUserArgs(
        list(z = z, n = n, coverage = coverage),
        "ser",
        userArgs
    )
    .setFinemappingFitClass(exec(susieR::susie_ser, !!!baseArgs), "susieRss")
}

# Fit every requested univariate token on one residualized (X, y) block,
# returning a named list (token -> FineMappingRow). Extracted from the
# univariate dispatch so the same logic serves the cis path (one block), the
# jointRegions=TRUE path (one concatenated block) and the jointRegions=FALSE
# path (one block per region, merged afterwards via .fmMergeEntries).
# One token's fit plus its postprocessing, or NULL when the fit did not run.
# @noRd
.fmXFitAndPostprocess <- function(
    tk,
    chainLocal,
    infFit,
    X,
    y,
    credibleSetArgs,
    methodArgs,
    verbose,
    ctx,
    tid,
    af,
    fitRetention
) {
    fit <- .fmXFitOne(
        tk,
        chainLocal,
        infFit,
        X = X,
        y = y,
        coverage = credibleSetArgs$coverage,
        methodArgs = methodArgs,
        verbose = verbose,
        ctx = ctx,
        tid = tid
    )
    if (is.null(fit)) {
        return(NULL)
    }
    .fmXPostprocess(
        fit,
        tk,
        X = X,
        y = y,
        credibleSetArgs = credibleSetArgs,
        af = af,
        fitRetention = fitRetention
    )
}

# Fit every requested method on one individual-level block. The susie-inf
# pre-fit is built first because the chained methods take it as their
# starting point; `chainLocal` says which of them actually chain.
# @noRd
.fmXFitAll <- function(
    toRun,
    addSusieInf,
    X,
    y,
    credibleSetArgs,
    methodArgs,
    verbose,
    ctx,
    tid,
    af,
    fitRetention
) {
    chainLocal <- .fmResolveSusieChain(toRun, addSusieInf)
    infFit <- .fmXInfFit(
        chainLocal,
        X = X,
        y = y,
        coverage = credibleSetArgs$coverage,
        methodArgs = methodArgs,
        verbose = verbose,
        ctx = ctx,
        tid = tid
    )
    compact(set_names(
        map(
            toRun,
            .fmXFitAndPostprocess,
            chainLocal = chainLocal,
            infFit = infFit,
            X = X,
            y = y,
            credibleSetArgs = credibleSetArgs,
            methodArgs = methodArgs,
            verbose = verbose,
            ctx = ctx,
            tid = tid,
            af = af,
            fitRetention = fitRetention
        ),
        toRun
    ))
}

.fmFitXBlock <- function(
    X,
    y,
    toRun,
    addSusieInf,
    credibleSetArgs = CredibleSetParam(includeAllCs = FALSE),
    methodArgs,
    verbose,
    ctx,
    tid,
    cvFolds = 0,
    cvThreads = 1,
    samplePartition = NULL,
    af = NULL,
    fitRetention = "slim",
    seed = NULL
) {
    out <- .fmXFitAll(
        toRun,
        addSusieInf = addSusieInf,
        X = X,
        y = y,
        credibleSetArgs = credibleSetArgs,
        methodArgs = methodArgs,
        verbose = verbose,
        ctx = ctx,
        tid = tid,
        af = af,
        fitRetention = fitRetention
    )
    .fmXCrossValidate(
        out,
        X = X,
        y = y,
        coverage = credibleSetArgs$coverage,
        methodArgs = methodArgs,
        cvFolds = cvFolds,
        cvThreads = cvThreads,
        samplePartition = samplePartition,
        seed = seed,
        verbose = verbose,
        ctx = ctx,
        tid = tid
    )
}

# Fit the shared susieInf model once, if the requested chain needs it.
.fmXInfFit <- function(
    chainLocal,
    X,
    y,
    coverage,
    methodArgs,
    verbose,
    ctx,
    tid
) {
    if (!chainLocal$runInf) {
        return(NULL)
    }
    if (verbose >= 1) {
        msg <- glue(
            "Fitting susieInf for (context='{ctx}', trait='{tid}') ..."
        )
        inform(msg)
    }
    .fmFitSusieIndiv(
        X,
        y,
        "susieInf",
        coverage = coverage,
        userArgs = methodArgs[["susieInf"]]
    )
}

# Resolve the fit for one method token; NULL means "skip this token".
.fmXFitOne <- function(
    tk,
    chainLocal,
    infFit,
    X,
    y,
    coverage,
    methodArgs,
    verbose,
    ctx,
    tid
) {
    if (tk == "susieInf") {
        if (!chainLocal$keepInf) {
            return(NULL)
        }
        return(infFit)
    }
    chainFrom <- if (
        (tk == "susie" && chainLocal$chainSusie) ||
            (tk == "susieAsh" && chainLocal$chainAsh)
    ) {
        infFit
    } else {
        NULL
    }
    if (verbose >= 1) {
        msg <- glue(
            "Fitting {tk} for (context='{ctx}', trait='{tid}') ..."
        )
        inform(msg)
    }
    .fmFitSusieIndiv(
        X,
        y,
        tk,
        chainFromInf = chainFrom,
        coverage = coverage,
        userArgs = methodArgs[[tk]]
    )
}

# Post-process one individual-level fit into a finemapping entry.
.fmXPostprocess <- function(
    fit,
    tk,
    X,
    y,
    credibleSetArgs,
    af,
    fitRetention
) {
    .fmPostprocessOne(
        fit = fit,
        method = tk,
        dataX = X,
        dataY = y,
        af = af,
        csInput = "X",
        credibleSetArgs = credibleSetArgs,
        fitRetention = fitRetention
    )
}

# Per-fold cross-validation across the fitted methods; attach each method's
# out-of-fold predictions to its entry.
.fmXCrossValidate <- function(
    out,
    X,
    y,
    coverage,
    methodArgs,
    cvFolds,
    cvThreads,
    samplePartition,
    seed,
    verbose,
    ctx,
    tid
) {
    if (!(cvFolds > 1L && length(out) > 0L)) {
        return(out)
    }
    if (verbose >= 1) {
        msg <- glue(
            "Cross-validating ({cvFolds} folds) for ",
            "(context='{ctx}', trait='{tid}') ..."
        )
        inform(msg)
    }
    cv <- .fmWeightsCv(
        X,
        y,
        names(out),
        methodArgs,
        cvFolds,
        samplePartition = samplePartition,
        coverage = coverage,
        verbose = verbose,
        numThreads = cvThreads,
        seed = seed
    )
    set_names(
        map(names(out), .fmAttachCvAt, out = out, cv = cv),
        names(out)
    )
}

# @noRd
.fmAttachCvAt <- function(tk, out, cv) {
    .fmAttachCv(out[[tk]], .fmSliceCv(cv, tk))
}

# Fit every requested RSS token on one (z, R, n) sumstat block, returning a
# named list (token -> FineMappingRow). The sumstat analog of .fmFitXBlock:
# the QtlSumStats and GwasSumStats methods both call it and differ only in how
# they push the returned entries (tuple shape) and the progress `label`.
.fmFitRssBlock <- function(
    z,
    R,
    n,
    toRun,
    addSusieInf,
    credibleSetArgs = CredibleSetParam(includeAllCs = FALSE),
    methodArgs,
    verbose,
    label,
    af = NULL,
    nVar = NULL,
    fitRetention = "slim",
    rssArgs
) {
    chainLocal <- .fmResolveSusieChain(toRun, addSusieInf)
    infFit <- .fmRssInfFit(
        chainLocal,
        z = z,
        R = R,
        n = n,
        coverage = credibleSetArgs$coverage,
        methodArgs = methodArgs,
        rssArgs = rssArgs,
        verbose = verbose,
        label = label
    )
    compact(set_names(
        map(
            toRun,
            .fmRssFitAndPostprocess,
            chainLocal = chainLocal,
            infFit = infFit,
            z = z,
            R = R,
            n = n,
            credibleSetArgs = credibleSetArgs,
            methodArgs = methodArgs,
            rssArgs = rssArgs,
            verbose = verbose,
            label = label,
            af = af,
            nVar = nVar,
            fitRetention = fitRetention
        ),
        toRun
    ))
}

# One RSS token's fit, postprocessed and (when it fell back to the single
# effect model) labelled as such. NULL when the fit did not run.
# @noRd
.fmRssFitAndPostprocess <- function(
    tk,
    chainLocal,
    infFit,
    z,
    R,
    n,
    credibleSetArgs,
    methodArgs,
    rssArgs,
    verbose,
    label,
    af,
    nVar,
    fitRetention
) {
    f <- .fmRssFitOne(
        tk,
        chainLocal,
        infFit,
        z = z,
        R = R,
        n = n,
        coverage = credibleSetArgs$coverage,
        methodArgs = methodArgs,
        rssArgs = rssArgs,
        verbose = verbose,
        label = label
    )
    if (is.null(f)) {
        return(NULL)
    }
    entry <- .fmRssPostprocess(
        f$fit,
        R = R,
        z = z,
        credibleSetArgs = credibleSetArgs,
        af = af,
        nVar = nVar,
        fitRetention = fitRetention
    )
    if (f$isStd && isTRUE(rssArgs$serFallback)) {
        return(.fmRssRecordFallback(entry, f, rssArgs$keepFullFit))
    }
    entry
}

# Fit the shared susieInf (RSS) model once, if the requested chain needs it.
.fmRssInfFit <- function(
    chainLocal,
    z,
    R,
    n,
    coverage,
    methodArgs,
    rssArgs,
    verbose,
    label
) {
    if (!chainLocal$runInf) {
        return(NULL)
    }
    if (verbose >= 1) {
        msg <- glue("Fitting susieInf (RSS) for {label} ...")
        inform(msg)
    }
    .fmFitSusieRss(
        z,
        R,
        n,
        "susieInf",
        coverage = coverage,
        userArgs = methodArgs[["susieInf"]],
        rFinite = rssArgs$rFinite,
        rMismatch = rssArgs$rMismatch,
        rssControl = .rssControlList(rssArgs$control)
    )
}

# Standard multi-effect SuSiE-RSS fit (susie / susieAsh): the only branch that
# carries susieR's finite-sample R diagnostics and honours the SER fallback.
.fmRssFitStd <- function(
    tk,
    chainLocal,
    infFit,
    z,
    R,
    n,
    coverage,
    methodArgs,
    rssArgs,
    verbose,
    label
) {
    chainFrom <- if (
        (tk == "susie" && chainLocal$chainSusie) ||
            (tk == "susieAsh" && chainLocal$chainAsh)
    ) {
        infFit
    } else {
        NULL
    }
    if (verbose >= 1) {
        msg <- glue("Fitting {tk} (RSS) for {label} ...")
        inform(msg)
    }
    fit <- .fmFitSusieRss(
        z,
        R,
        n,
        tk,
        chainFromInf = chainFrom,
        coverage = coverage,
        userArgs = methodArgs[[tk]],
        rFinite = rssArgs$rFinite,
        rMismatch = rssArgs$rMismatch,
        rssControl = .rssControlList(rssArgs$control)
    )
    .fmRssSerFallback(fit, rssArgs)
}

# An unreliable LD matrix makes the multi-effect fit untrustworthy, so with
# `serFallback` the single-effect model susie_rss already produced is
# returned instead -- with the multi-effect fit kept alongside it, since
# the caller may still want to inspect what was rejected.
# @noRd
.fmRssSerFallback <- function(fit, rssArgs) {
    rfd <- fit$R_finite_diagnostics
    flag <- if (!is.null(rfd) && !is.null(rfd$R_reliability_flag)) {
        isTRUE(rfd$R_reliability_flag)
    } else {
        NA
    }
    if (
        !isTRUE(rssArgs$serFallback) || !isTRUE(flag) || is.null(rfd$ser_model)
    ) {
        return(list(fit = fit, flag = flag, multiFit = NULL))
    }
    list(
        fit = .setFinemappingFitClass(rfd$ser_model, "susieRss"),
        flag = flag,
        multiFit = fit
    )
}

# Resolve the fit for one method token; NULL means "skip this token".
.fmRssFitOne <- function(
    tk,
    chainLocal,
    infFit,
    z,
    R,
    n,
    coverage,
    methodArgs,
    rssArgs,
    verbose,
    label
) {
    if (tk == "susieInf") {
        if (!chainLocal$keepInf) {
            return(NULL)
        }
        return(list(fit = infFit, flag = NA, isStd = FALSE, multiFit = NULL))
    }
    if (tk == "ser") {
        if (verbose >= 1) {
            msg <- glue("Fitting ser (RSS single-effect) for {label} ...")
            inform(msg)
        }
        fit <- .fmFitSusieSer(
            z,
            n,
            coverage = coverage,
            userArgs = methodArgs[["ser"]]
        )
        return(list(fit = fit, flag = NA, isStd = FALSE, multiFit = NULL))
    }
    std <- .fmRssFitStd(
        tk,
        chainLocal,
        infFit,
        z = z,
        R = R,
        n = n,
        coverage = coverage,
        methodArgs = methodArgs,
        rssArgs = rssArgs,
        verbose = verbose,
        label = label
    )
    list(fit = std$fit, flag = std$flag, isStd = TRUE, multiFit = std$multiFit)
}

# Post-process one RSS fit into a finemapping entry.
.fmRssPostprocess <- function(
    fit,
    R,
    z,
    credibleSetArgs,
    af,
    nVar,
    fitRetention
) {
    .fmPostprocessOne(
        fit = fit,
        method = "susieRss",
        dataX = R,
        dataY = list(z = z),
        af = af,
        # Per-variant effective N (reporting-only, top_loci$N). NULL on any
        # path that has no per-variant N -> buildTopLoci leaves N as NA,
        # never 1. This is NOT `n` (the scalar median the RSS fit consumes).
        n = nVar,
        csInput = "Xcorr",
        credibleSetArgs = credibleSetArgs,
        fitRetention = fitRetention
    )
}

# Record the SER-fallback reliability decision (and retained multi-effect fit)
# on the entry's susieFit list. Gated on serFallback so the default path is
# byte-identical.
.fmRssRecordFallback <- function(ent, f, keepFullFit) {
    multiEffectFit <- if (
        !is.null(f$multiFit) && is_in(keepFullFit, c("fallback", "all"))
    ) {
        f$multiFit
    } else if (identical(keepFullFit, "all")) {
        f$fit
    }
    sf <- list_assign(
        .fmrPartsSusieFit(ent),
        !!!compact(list(
            R_reliability_flag = f$flag,
            serFallbackUsed = isTRUE(f$flag),
            multiEffectFit = multiEffectFit
        ))
    )
    # A row is immutable: rebuild it with the amended fit rather than
    # assigning into it. The variants and topLoci are unchanged, so this
    # round-trips them through the same builder the caller used.
    fineMappingRow(
        variantIds = getVariantIds(ent),
        susieFit = sf,
        topLoci = .fmrPartsTopLoci(ent),
        cvResult = getCvResult(ent)
    )
}

# =============================================================================
# Named helpers for map/apply call sites (no inline lambdas)
# =============================================================================

# Translate one legacy cs_coverage_<C> column name to the canonical form.
# @noRd
.translateOneLegacyCsColumn <- function(x) {
    x <- as.character(x)
    oldParts <- str_match(
        x,
        regex("^cs_coverage_([0-9.]+)$", ignore_case = TRUE)
    )[1L, ]
    if (!is.na(oldParts[[1L]])) {
        return(formatCsColumn(as.numeric(oldParts[[2L]]), "susie"))
    }
    x
}

# @noRd
.camelToSnakeOne <- function(m, lookup) {
    if (is_in(m, names(lookup))) lookup[[m]] else m
}

# Post-process one method's fit (buildTopLoci per fit).
# @noRd
.ppOneFit <- function(
    method,
    fits,
    dataX,
    dataY,
    xScalar,
    yScalar,
    af,
    n,
    credibleSetArgs,
    fitRetention,
    otherQuantities,
    region,
    priorEffTol,
    csInput,
    conditionIdx
) {
    fit <- .setFinemappingFitClass(fits[[method]], method)
    postprocessFinemappingFit(
        fit,
        method = method,
        dataX = dataX,
        dataY = dataY,
        xScalar = xScalar,
        yScalar = yScalar,
        af = af,
        n = n,
        credibleSetArgs = credibleSetArgs,
        fitRetention = fitRetention,
        otherQuantities = otherQuantities,
        region = region,
        priorEffTol = priorEffTol,
        csInput = csInput,
        conditionIdx = conditionIdx
    )
}

# @noRd
.ppDropTopLoci <- function(x) {
    list_modify(x, top_loci = zap())
}

# @noRd
.computeCsTableForCov <- function(
    cov,
    fit,
    dataX,
    csInput,
    minAbsCorr,
    medianAbsCorr
) {
    computeCsTable(
        fit,
        dataX,
        coverage = cov,
        csInput = csInput,
        minAbsCorr = minAbsCorr,
        medianAbsCorr = medianAbsCorr
    )
}

# Purity value for the i-th index (0 out of range / NA).
# @noRd
.csPurityAt <- function(i, pv) {
    if (i <= 0L || i > length(pv)) {
        return(0)
    }
    v <- pv[i]
    if (is.na(v)) 0 else as.numeric(v)
}

# Max of the finite entries of x (NA when none).
# @noRd
.finiteMax <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else max(x)
}

# @noRd
.btlCsLabel <- function(ix, methodTag) {
    str_c(methodTag, "_", ix)
}

# Trailing integer of a "<method>_<idx>" cs string (0 when empty / NA).
# @noRd
.cs95ToIndex <- function(s) {
    if (is.na(s) || str_length(s) == 0L) {
        return(0L)
    }
    suppressWarnings(as.integer(str_remove(s, "^.*_")))
}

# @noRd
.dropPipCol <- function(x) {
    x[names(x) != "pip"]
}

# @noRd
.csContains <- function(x, snpsIdx) {
    is_in(snpsIdx, x)
}

# One (variant, CS) block of rows for variant `vi`.
# @noRd
.csInfoRow <- function(vi, susieOutputSetsCs, csNames) {
    idx <- getCsIndex(vi, susieOutputSetsCs)
    if (length(idx) == 1 && is.na(idx)) {
        return(tibble(
            variant_idx = vi,
            cs_idx = 0L
        ))
    }
    csNums <- as.integer(str_replace(csNames[idx], "L", ""))
    tibble(
        variant_idx = rep(vi, length(csNums)),
        cs_idx = csNums
    )
}

# @noRd
.asDataFrameT <- function(x) {
    as.data.frame(t(x))
}

# @noRd
.amRow <- function(l, am) {
    as.numeric(am[l, ])
}

# Merged credible-set label for one variant (mapped set name, else joined set
# ids).
# @noRd
.updateCredibleSet <- function(variantId, variantsSetsAndPipsList, setNameMap) {
    currentSets <- variantsSetsAndPipsList[[variantId]][["sets"]]
    mapped <- intersect(currentSets, names(setNameMap))
    if (length(mapped) > 0) {
        setNameMap[[mapped[1]]]
    } else {
        str_flatten(sort(unique(currentSets)), ",")
    }
}

# One merged-CS summary row (variant id + set label + max/median PIP).
# @noRd
.csMergedVariantRow <- function(
    variantId,
    extractedResult,
    hasOverlaps,
    mergedSets
) {
    credibleSetNames <- if (hasOverlaps) {
        mergedSets[[variantId]]
    } else {
        str_flatten(
            sort(unique(extractedResult[[variantId]]$sets)),
            ","
        )
    }
    tibble(
        variant_id = variantId,
        credibleSetNames = credibleSetNames,
        maxPip = max(extractedResult[[variantId]]$pips),
        medianPip = median(extractedResult[[variantId]]$pips)
    )
}

# (variant, CS) rows for the i-th entry, labelled cs_<entry>_<set>.
# @noRd
.extractCsEntryRows <- function(i, entries, csCol) {
    topLoci <- .translateLegacyTopLociCsColumns(
        .fmrRowTopLoci(entries[[i]])
    )
    if (
        is.null(topLoci) ||
            nrow(topLoci) == 0 ||
            !is_in(csCol, names(topLoci))
    ) {
        return(NULL)
    }
    pipCol <- resolvePipColumn(topLoci)
    if (is.null(pipCol)) {
        return(NULL)
    }
    csIdx <- .fmCsIdx(topLoci[[csCol]])
    uniq <- unique(csIdx)
    setNum <- uniq[!is.na(uniq) & uniq != 0]
    if (length(setNum) == 0) {
        return(NULL)
    }
    map_dfr(
        setNum,
        .extractCsSetRows,
        csIdx = csIdx,
        topLoci = topLoci,
        pipCol = pipCol,
        i = i
    )
}

# @noRd
.extractCsSetRows <- function(sn, csIdx, topLoci, pipCol, i) {
    keep <- !is.na(csIdx) & csIdx == sn
    topLoci |>
        filter(keep) |>
        select("variant_id", pip = all_of(pipCol)) |>
        mutate(set_name = str_c("cs_", i, "_", sn))
}

# @noRd
.csSplitToList <- function(df) {
    list(sets = df$set_name, pips = df$pip)
}

# Between-credible-set correlation for one row, dispatched on the ldSource
# kind. The public computeCsCorrelation() generic wraps this; the pipeline
# calls it directly because it carries row payloads rather than collections.
# @noRd
.rowCsCorrelation <- function(parts, ldSource) {
    if (methods::is(ldSource, "SumStatsBase")) {
        return(.rowCsCorrelationSumstats(parts, ldSource))
    }
    if (methods::is(ldSource, "QtlDataset")) {
        return(.rowCsCorrelationGeno(parts, ldSource))
    }
    abort(glue(
        "computeCsCorrelation() requires a QtlDataset, QtlSumStats, or ",
        "GwasSumStats as `ldSource`: the between-credible-set correlation ",
        "is derived from that object's LD and is never stored on the fit."
    ))
}

# The canonical (non-reweighted) mvSuSiE mixture prior for residual variance
# `V`:
# create_mixture_prior(R) restricted to the group's conditions.
# @noRd
.fmCanonicalPrior <- function(V, conditionNames, R) {
    list(
        priorVariance = mvsusieR::create_mixture_prior(
            R = R,
            include_indices = conditionNames
        ),
        residualVariance = V
    )
}

# Rebuild the mvSuSiE data-driven *reweighted* mixture prior + residual variance
# from a stored mr.mash fit -- the lean payload
# (list(dataDrivenPriorMatrices, w0, V)) that mrmashWeights() attaches at
# fitRetention "slim" and twasWeightsPipeline keeps on the mrmash row. Shared
# by the fine-mapping mvsusie consumer and the twas mvsusie_weights consumer.
#
# Reproduces the deleted multivariate_pipeline.R reweighting bit-identically:
# rescaleCovW0(w0) collapses the expanded mr.mash weights onto the original
# data-driven covariance matrices ($U), filters to surviving components, and
# create_mixture_prior() wraps them, restricted to the fit's conditions
# (`conditionNames` = colnames(Y)). `V` becomes mvsusie's residual_variance.
# A NULL fit, NULL matrices, or no surviving component falls back to the
# canonical create_mixture_prior(R), matching the legacy `else` branch.
# Returns list(priorVariance, residualVariance) (residualVariance NULL only
# when no fit was supplied at all).
# @noRd
.buildMvsusieReweightedPrior <- function(
    fitParts,
    conditionNames,
    weightsTol = 1e-10,
    overrideU = NULL
) {
    R <- length(conditionNames)
    if (is.null(fitParts)) {
        return(.fmCanonicalPrior(NULL, conditionNames, R))
    }
    # `overrideU` (mode C / hybrid): reuse this fit's reweighted mixture weights
    # (w0) and residual variance (V) but swap in a different set of data-driven
    # covariance matrices -- the per-fold mash prior U. Components are matched
    # to w0 by name, so the override U must share component names with the fit.
    ddpm <- if (!is.null(overrideU)) {
        overrideU
    } else {
        fitParts$dataDrivenPriorMatrices
    }
    if (is.null(ddpm) || is.null(ddpm$U)) {
        return(.fmCanonicalPrior(fitParts$V, conditionNames, R))
    }
    rescaled <- rescaleCovW0(fitParts$w0)
    w0Updated <- rescaled[is_in(names(rescaled), names(ddpm$U))]
    if (length(w0Updated) == 0L) {
        return(.fmCanonicalPrior(fitParts$V, conditionNames, R))
    }
    mixture <- list(matrices = ddpm$U[names(w0Updated)], weights = w0Updated)
    list(
        priorVariance = mvsusieR::create_mixture_prior(
            mixture_prior = mixture,
            weights_tol = weightsTol,
            include_indices = conditionNames
        ),
        residualVariance = fitParts$V
    )
}

# Per-column marginal-association z-scores of y on each column of X (univariate
# regression z = betahat / sebetahat), used by the individual-level absZ screen.
# @noRd
.marginalZ <- function(X, y) {
    ur <- susieR::univariate_regression(X, y)
    ur$betahat / ur$sebetahat
}

# Single-effect (SER) pre-screen, individual-level. Reports whether a
# residualized (X, y) block shows a strong enough signal (by the chosen metric)
# to be worth a full fit. `screen` is a screen spec (see .asScreen): a legacy
# PIP cutoff (numeric scalar, 0 = off) OR a resolved list(metric, cutoff) for
# one of pip / absZ / bf / logBf. Ports the deleted multivariate_pipeline.R
# `skipConditions` / susie_twas `pip_cutoff_to_skip` logic (the individual-level
# analog of the sumstat-path `.applyEntryScreen`):
#   * no screen (NULL / 0 / non-scalar numeric) -> always keep.
#   * pip cutoff < 0 uses the adaptive 3 / nVariants threshold.
#   * absZ needs no susie fit; pip/bf/logBf fit susie L = 1 once (its
#     $lbf_variable gives the per-variant logBF for bf/logBf, $pip for pip).
#   * NA entries of `y` are dropped before fitting.
# The screen is advisory: too few samples/variants or a fit failure returns
# `fallback` -- TRUE (default) keeps the block rather than discard a potentially
# real signal (fineMapping / joint paths); colocboost passes FALSE to drop an
# outcome it cannot screen. This is the single L = 1 SuSiE pre-screen shared by
# .fmSerScreenColumns (joint) and .cbPipSkipOutcomes (colocboost).
# @noRd
.fmSerScreen <- function(X, y, screen, fallback = TRUE) {
    scr <- .asScreen(screen)
    if (is.null(scr)) {
        return(TRUE)
    }
    ok <- !is.na(y)
    if (sum(ok) < 2L || ncol(X) < 1L) {
        return(fallback)
    }
    raw <- X[ok, , drop = FALSE]
    # susieR needs a double X.
    Xs <- if (is.double(raw)) raw else `storage.mode<-`(raw, "double")
    ys <- y[ok]
    metric <- scr$metric
    cutoff <- scr$cutoff
    if (metric == "absZ") {
        z <- try_fetch(.marginalZ(Xs, ys), error = function(cnd) NULL)
        if (is.null(z)) {
            return(fallback)
        }
        return(any(abs(z) > cutoff, na.rm = TRUE))
    }
    fit <- try_fetch(
        suppressMessages(susieR::susie(Xs, ys, L = 1L)),
        error = function(cnd) NULL
    )
    if (is.null(fit)) {
        return(fallback)
    }
    if (metric == "pip") {
        thr <- if (cutoff < 0) 3 / ncol(Xs) else cutoff
        return(any(fit$pip > thr, na.rm = TRUE))
    }
    maxLbf <- suppressWarnings(max(as.numeric(fit$lbf_variable), na.rm = TRUE))
    if (!is.finite(maxLbf)) {
        return(fallback)
    }
    # bf: cutoff on the raw BF scale -> compare in log space; logBf: log scale.
    maxLbf > (if (metric == "bf") log(cutoff) else cutoff)
}

# Per-fold mvsusie weights. Reuses the data-driven reweighted prior + residual
# covariance from the full-data mr.mash fit on every fold -- the prior is over
# conditions, identical across folds (only samples are held out). NULL mvPrior
# -> canonical prior (unchanged behavior).
# @noRd
.fmFoldWeightsMv <- function(Xtr, Ytr, coverage, userArgs, mvPrior) {
    baseArgs <- list(
        X = Xtr,
        Y = Ytr,
        coverage = coverage,
        prior_variance = if (is.null(mvPrior)) {
            mvsusieR::create_mixture_prior(R = ncol(Ytr))
        } else {
            mvPrior$priorVariance
        }
    )
    withPrior <- list_assign(
        baseArgs,
        !!!compact(list(residual_variance = mvPrior$residualVariance))
    )
    mvArgs <- .fmMergeUserArgs(withPrior, "mvsusie", userArgs)
    fit <- exec(fitMvsusie, !!!.splitMethodArgs(fitMvsusie, mvArgs))
    raw <- as.matrix(mvsusieWeights(mvsusieFit = fit))
    W <- `rownames<-`(raw, rownames(raw) %||% colnames(Xtr))
    `attr<-`(W, "fit", .fmLeanFoldFit(fit, "mvsusie"))
}

#' @rdname fineMappingMethodOptions
#' @export
SusieOptions <- function(...) {
    .fmMethodOptions("susieR::susie", "SusieOptions", "susie", list(...))
}

#' @rdname fineMappingMethodOptions
#' @export
SusieRssOptions <- function(...) {
    .fmMethodOptions(
        "susieR::susie_rss",
        "SusieRssOptions",
        "susieRss",
        list(...)
    )
}

#' @rdname fineMappingMethodOptions
#' @export
SusieInfOptions <- function(...) {
    .fmMethodOptions("susieR::susie", "SusieInfOptions", "susieInf", list(...))
}

#' @rdname fineMappingMethodOptions
#' @export
SusieInfRssOptions <- function(...) {
    .fmMethodOptions(
        "susieR::susie_rss",
        "SusieInfRssOptions",
        "susieInfRss",
        list(...)
    )
}

#' @rdname fineMappingMethodOptions
#' @export
SusieAshOptions <- function(...) {
    .fmMethodOptions("susieR::susie", "SusieAshOptions", "susieAsh", list(...))
}

#' @rdname fineMappingMethodOptions
#' @export
SusieAshRssOptions <- function(...) {
    .fmMethodOptions(
        "susieR::susie_rss",
        "SusieAshRssOptions",
        "susieAshRss",
        list(...)
    )
}

#' @rdname fineMappingMethodOptions
#' @export
SerOptions <- function(...) {
    .fmMethodOptions("susieR::susie_ser", "SerOptions", "ser", list(...))
}

#' @rdname fineMappingMethodOptions
#' @export
MvsusieOptions <- function(...) {
    .fmMethodOptions(
        "mvsusieR::mvsusie",
        "MvsusieOptions",
        "mvsusie",
        list(...)
    )
}

#' @rdname fineMappingMethodOptions
#' @export
MvsusieRssOptions <- function(...) {
    .fmMethodOptions(
        "mvsusieR::mvsusie_rss",
        "MvsusieRssOptions",
        "mvsusieRss",
        list(...)
    )
}

#' @rdname fineMappingMethodOptions
#' @export
FsusieOptions <- function(...) {
    .fmMethodOptions("fsusieR::susiF", "FsusieOptions", "fsusie", list(...))
}

#' @title Arguments For susieR's RSS Control Block
#' @description Options forwarded to \code{susieR::susie_rss_control()} and
#'   from there to \code{susie_rss()}'s \code{control} argument. Names are
#'   checked against that function's live formals, so a misspelling fails at
#'   the call site rather than being dropped into an ignored list --- which is
#'   what a bare named list here used to do.
#' @param ... Arguments for \code{susieR::susie_rss_control()}, under its own
#'   names (\code{check_prior}, \code{mismatch_estimator}, ...).
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' SusieRssControlOptions(check_prior = TRUE)
#' @export
SusieRssControlOptions <- function(...) {
    .newMethodOptions(
        "susieR::susie_rss_control",
        defaults = list(),
        extra = list(...),
        label = "SusieRssControlOptions",
        engine = "susieRssControl"
    )
}

#' @title Options for the susieR Kriging RSS Diagnostic
#' @description Build a checked record of extra arguments for
#'   \code{susieR::kriging_rss()}, the engine behind
#'   \code{\link{krigingOutlierQc}}.
#' @param ... Arguments for \code{susieR::kriging_rss()}: \code{r_tol} (the
#'   eigenvalue tolerance below which an LD eigenvalue is treated as zero) and
#'   \code{s} (the estimated proportion of LD/sumstats mismatch, which
#'   defaults to \code{susieR::estimate_s_rss()}). \code{z}, \code{R} and
#'   \code{n} are supplied by pecotmr and refused.
#' @return A \code{MethodOptions} record for
#'   \code{krigingOutlierQc(methodArgs =)}.
#' @seealso \code{\link{krigingOutlierQc}}
#' @examples
#' KrigingOptions(r_tol = 1e-06)
#' @export
KrigingOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            z = "the caller's `zScore`",
            R = "the caller's `R`",
            n = "the caller's `n`"
        ),
        "KrigingOptions"
    )
    .newMethodOptions(
        "susieR::kriging_rss",
        defaults = list(),
        extra = extra,
        label = "KrigingOptions",
        engine = "kriging"
    )
}

#' Kriging-style LD-consistency outlier QC
#'
#' Flags variants whose observed z-score is inconsistent with the value
#' predicted from its LD neighbours, using susieR's kriging diagnostic.
#' \code{susieR::kriging_rss()} computes the leave-one-out conditional
#' distribution of each \code{z_i} given the rest (with the LD-mismatch scale
#' \code{s} defaulting to \code{susieR::estimate_s_rss()}) and a per-variant
#' \code{logLR} for the allele-switch hypothesis. This helper reuses susieR's
#' own allele-switch rule --- \code{logLR > logLRThreshold & abs(z) >
#' zThreshold} (the same \code{logLR > 2 & |z| > 2} used in
#' \code{susie_rss_utils}) --- to flag variants whose sign should be flipped.
#' RSS-only helper, opt-in via \code{alleleFlipKriging}; never wired into
#' \code{alleleQc()} / \code{matchRefPanel()}. Requires a susieR that provides
#' \code{kriging_rss()} and \code{estimate_s_rss()}.
#'
#' @param zScore Numeric vector of harmonized z-scores.
#' @param R Square LD correlation matrix aligned to \code{zScore}.
#' @param n Sample size, forwarded to \code{susieR::kriging_rss()} (whose
#'   default \code{s} is \code{susieR::estimate_s_rss()}).
#' @param variantIds Optional variant IDs for the diagnostics table.
#' @param zThreshold Absolute-z cutoff for the allele-switch rule (default
#'   \code{2}, matching susieR).
#' @param logLRThreshold Log-likelihood-ratio cutoff for the allele-switch rule
#'   (default \code{2}, matching susieR).
#' @param methodArgs Extra arguments for \code{susieR::kriging_rss()}, built
#'   with \code{\link{KrigingOptions}} -- \code{r_tol} and \code{s}.
#' @return A list with \code{flip} (logical vector; \code{TRUE} = allele switch,
#'   z-score should be sign-flipped) and \code{diagnostics} (data frame of
#'   per-variant \code{z}, \code{condmean}, \code{z_std_diff}, \code{logLR}, and
#'   the \code{flipped} flag).
#' @importFrom stats pnorm
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:20]
#' R <- cor(X)
#' krigingOutlierQc(
#'   zScore = rnorm(20), R = R, n = 415, variantIds = colnames(X))
#' @export
#' @importFrom checkmate assertMatrix
krigingOutlierQc <- function(
    zScore,
    R,
    n,
    variantIds = NULL,
    zThreshold = 2,
    logLRThreshold = 2,
    methodArgs = KrigingOptions()
) {
    .assertMethodOptions(methodArgs, "KrigingOptions", "methodArgs")
    zScore <- as.numeric(zScore)
    m <- length(zScore)
    assertMatrix(R, nrows = m, ncols = m, .var.name = "R (LD matrix)")
    if (missing(n) || length(n) != 1L || is.na(n) || !is.finite(n) || n <= 0) {
        abort("krigingOutlierQc requires a single positive sample size 'n'.")
    }
    .krigingCheckSusie()
    if (is.null(variantIds)) {
        variantIds <- rownames(R)
    }
    # susieR kriging RSS diagnostic (kriging_rss(z, R, n)): s defaults to
    # estimate_s_rss(); logLR matches susieR's allele-switch selection.
    cd <- exec(
        susieR::kriging_rss,
        z = zScore,
        R = R,
        n = n,
        !!!as.list(methodArgs)
    )$conditional_dist
    condMean <- as.numeric(cd$condmean)
    zStdDiff <- as.numeric(cd$z_std_diff)
    logLR <- as.numeric(cd$logLR)
    # susieR's allele-switch rule (susie_rss_utils.R): logLR > 2 & |z| > 2.
    flip <- !is.na(logLR) &
        !is.na(zScore) &
        logLR > logLRThreshold &
        abs(zScore) > zThreshold
    list(
        flip = flip,
        diagnostics = tibble(
            variant_id = if (is.null(variantIds)) seq_len(m) else variantIds,
            z = zScore,
            condmean = condMean,
            z_std_diff = zStdDiff,
            logLR = logLR,
            flipped = flip
        )
    )
}

# Decide whether the z-scores of one entry clear the chosen screen. Returns
# list(ok = logical, reason = character). susie_ser is fit at most once, and
# only for the metrics that need it (absZ stays model-free).
.entryScreenPass <- function(z, n, nVar, scr) {
    metric <- scr$metric
    cutoff <- scr$cutoff
    if (metric == "absZ") {
        m <- suppressWarnings(max(abs(as.numeric(z)), na.rm = TRUE))
        return(list(
            ok = is.finite(m) && m > cutoff,
            reason = sprintf(
                "no variant with |Z| above %g (max |Z| = %g)",
                cutoff,
                m
            )
        ))
    }
    ser <- susieR::susie_ser(z = z, n = n, coverage = NULL)
    if (metric == "pip") {
        eff <- if (cutoff < 0) 3 / nVar else cutoff
        return(list(
            ok = any(ser$pip > eff),
            reason = sprintf("no signals above PIP threshold %g", eff)
        ))
    }
    maxLbf <- suppressWarnings(max(as.numeric(ser$lbf_variable), na.rm = TRUE))
    if (metric == "logBf") {
        return(list(
            ok = is.finite(maxLbf) && maxLbf > cutoff,
            reason = sprintf(
                "no variant with logBF above %g (max logBF = %g)",
                cutoff,
                maxLbf
            )
        ))
    }
    # metric == "bf": compare in log space to avoid overflow of exp(maxLbf).
    list(
        ok = is.finite(maxLbf) && maxLbf > log(cutoff),
        reason = sprintf(
            "no variant with BF above %g (max BF = %g)",
            cutoff,
            exp(maxLbf)
        )
    )
}

# Populate `obj$cred_band` via fsusieR's wavethresh/GenW band computation. That
# function is registered as an S3 method but NOT exported, and the exported
# affected_reg() depends on cred_band already being populated, so this internal
# call is the only path that works across all post_processing modes. Guarded so
# an upstream fsusieR change surfaces as a clear error, not a silent NULL.
# @noRd
#' @importFrom rlang try_fetch
.fsusiePopulateCredibleBand <- function(fit) {
    fn <- try_fetch(
        get("update_cal_credible_band.susiF", envir = asNamespace("fsusieR")),
        error = function(cnd) NULL
    )
    if (is.null(fn)) {
        # Defensive guard against an upstream fsusieR rename; only reachable
        # if fsusieR drops this unexported S3 method.
        msg <- glue(
            "fsusieR's internal update_cal_credible_band.susiF not found; ",
            "cannot compute the fSuSiE credible band (upstream fsusieR API ",
            "changed)."
        )
        abort(msg)
    }
    indxLst <- fsusieR::gen_wavelet_indx(log2(length(fit$outing_grid)))
    fn(fit, indxLst)
}

# @noRd
.fsusieAffectedRegionsFit <- function(fit, topLoci = NULL) {
    if (!.isFsusieFit(fit)) {
        return(GenomicRanges::GRanges())
    }
    fit <- .fsusiePopulateCredibleBand(fit)
    raw <- try_fetch(fsusieR::affected_reg(fit), error = function(cnd) NULL)
    if (is.null(raw) || nrow(raw) == 0L) {
        return(GenomicRanges::GRanges())
    }
    reg <- as_tibble(raw)
    chrom <- .fsusieChrom(fit)
    grid <- as.numeric(fit$outing_grid)
    csMap <- .fsusieCsMapFromTopLoci(fit, topLoci)
    csKey <- as.character(reg$CS)
    # Effect direction over each region (sign of the fitted effect curve), which
    # upstream affected_reg() collapses away.
    direction <- map_chr(
        seq_len(nrow(reg)),
        .fsusieRegionDirection,
        grid = grid,
        reg = reg,
        fit = fit
    )
    GenomicRanges::GRanges(
        seqnames = str_c("chr", str_remove(chrom, "^chr")),
        ranges = IRanges::IRanges(
            start = as.integer(reg$Start),
            end = as.integer(reg$End)
        ),
        cs = unname(csMap$label[csKey]),
        purity = unname(csMap$purity[csKey]),
        direction = direction
    )
}

# Require a susieR that provides the kriging RSS diagnostic.
.krigingCheckSusie <- function() {
    if (
        !requireNamespace("susieR", quietly = TRUE) ||
            !all(
                is_in(
                    c("estimate_s_rss", "kriging_rss"),
                    getNamespaceExports("susieR")
                )
            )
    ) {
        msg <- glue(
            "krigingOutlierQc requires a susieR that provides ",
            "estimate_s_rss() and kriging_rss(); the installed susieR does ",
            "not. Install a susieR with the kriging RSS diagnostic, or ",
            "disable alleleFlipKriging."
        )
        abort(msg)
    }
}

# susieR's `control` argument is a plain named list, so the constructor
# result is flattened on the way out. NULL is preserved rather than becoming
# list(): to susie_rss() an absent control means "use susie_rss_control()'s
# own defaults", which an empty list does not.
# @noRd
.rssControlList <- function(control) {
    if (is.null(control) || length(control) == 0L) {
        return(NULL)
    }
    as.list(control)
}
