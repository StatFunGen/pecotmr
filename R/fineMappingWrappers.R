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
    if (is.null(priorWeights)) {
        priorWeights <- rep(1 / length(lbf), length(lbf))
    }
    maxlbf <- max(lbf)

    # If maxlbf is 0, return a vector of zeros
    if (maxlbf == 0) {
        return(set_names(rep(0, length(lbf)), names(lbf)))
    }

    # w is proportional to BF, subtract max for numerical stability
    w <- exp(lbf - maxlbf)

    # Posterior prob for each SNP
    wWeighted <- w * priorWeights
    weightedSumW <- sum(wWeighted)
    alpha <- wWeighted / weightedSumW

    return(alpha)
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
    if (ncol(lbf) == 1) {
        alphaMatrix <- matrix(
            alphaMatrix,
            ncol = 1,
            dimnames = list(NULL, colnames(lbf))
        )
    }
    return(alphaMatrix)
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
    pipCols <- grep("^pip_", names(topLoci), value = TRUE)
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

.translateLegacyTopLociCsColumns <- function(topLoci) {
    if (!is.data.frame(topLoci)) {
        return(topLoci)
    }
    names(topLoci) <- .translateLegacyCsColumnName(names(topLoci))
    if (is_in("pip_susie", names(topLoci)) && !is_in("pip", names(topLoci))) {
        names(topLoci)[names(topLoci) == "pip_susie"] <- "pip"
    }
    topLoci
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
    if (!is.null(methodClass)) {
        class(fit) <- unique(c(methodClass, class(fit)))
    }
    fit
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
    L <- args[["L"]]
    if (is.null(L)) {
        L <- length(susieInfFit$V)
    }
    if (is.null(args[["refine"]]) && !is.null(refineDefault)) {
        args[["refine"]] <- refineDefault
    }
    args[["unmappable_effects"]] <- unmappableEffects
    args[["model_init"]] <- susieInfFit
    if (unmappableEffects == "ash") {
        args[["convergence_method"]] <- args[["convergence_method"]] %||% "pip"
    }
    if (!is.null(args[["L_greedy"]])) {
        args[["L_greedy"]] <- min(length(susieInfFit$V), L)
    }
    args
}

#' @importFrom utils modifyList
#' @noRd
fitSusieInfThenSusie <- function(
    X,
    y,
    args = list(),
    susieInfArgs = list(),
    susieArgs = list(),
    fittedModels = NULL
) {
    # Two-stage chain built from the shared per-token fitter (.fmFitSusieIndiv),
    # so the susieInf fit arguments and the susieInf -> susie initialisation
    # live in one place rather than being duplicated here and in the pipeline.
    if (is.null(fittedModels)) {
        fittedModels <- list()
    }
    susieInfFit <- fittedModels[["susieInf"]]
    if (is.null(susieInfFit)) {
        susieInfFit <- .fmFitSusieIndiv(
            X,
            y,
            "susieInf",
            userArgs = modifyList(args, susieInfArgs)
        )
    } else {
        susieInfFit <- .setFinemappingFitClass(susieInfFit, "susieInf")
    }
    susieFit <- fittedModels[["susie"]]
    if (is.null(susieFit)) {
        susieFit <- .fmFitSusieIndiv(
            X,
            y,
            "susie",
            chainFromInf = susieInfFit,
            userArgs = modifyList(args, susieArgs)
        )
    } else {
        susieFit <- .setFinemappingFitClass(susieFit, "susie")
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
#' @param args Default arguments forwarded to both fits.
#' @param susieInfArgs SuSiE-inf-specific overrides.
#' @param susieArgs Standard SuSiE-RSS-specific overrides.
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
#' @export
fitSusieInfThenSusieRss <- function(
    z,
    R,
    n,
    args = list(),
    susieInfArgs = list(),
    susieArgs = list(),
    fittedModels = NULL
) {
    # RSS analog of fitSusieInfThenSusie, built from the shared per-token RSS
    # fitter (.fmFitSusieRss). .fmFitSusieRss tags every fit "susieRss", so the
    # inf fit is re-tagged "susieInf" to preserve this wrapper's contract.
    if (is.null(fittedModels)) {
        fittedModels <- list()
    }
    susieInfFit <- fittedModels[["susieInf"]]
    if (is.null(susieInfFit)) {
        susieInfFit <- .fmFitSusieRss(
            z,
            R,
            n,
            "susieInf",
            userArgs = modifyList(args, susieInfArgs)
        )
    }
    susieInfFit <- .setFinemappingFitClass(susieInfFit, "susieInf")
    susieFit <- fittedModels[["susie"]]
    if (is.null(susieFit)) {
        susieFit <- .fmFitSusieRss(
            z,
            R,
            n,
            "susie",
            chainFromInf = susieInfFit,
            userArgs = modifyList(args, susieArgs)
        )
    }
    susieFit <- .setFinemappingFitClass(susieFit, "susieRss")

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
#' @param coverage Primary credible-set coverage.
#' @param secondaryCoverage Additional credible-set coverages.
#' @param signalCutoff PIP cutoff for including non-CS variants in top loci.
#' @param otherQuantities Optional list carried into each method result.
#' @param priorEffTol Tolerance for retaining effects by prior variance.
#' @param minAbsCorr Minimum absolute correlation for credible-set purity.
#' @param region Optional genomic anchor (\code{"chr:start-end"} or
#'   \code{GRanges}) recorded on the result; \code{NULL} to omit.
#' @param medianAbsCorr Numeric or \code{NULL}. Median absolute within-CS
#'   correlation threshold for purity; \code{NULL} uses only \code{minAbsCorr}.
#' @param csInput Optional precomputed credible-set specification, or
#'   \code{NULL} to derive it from the fits.
#' @param conditionIdx Integer or \code{NULL}. Index of the conditioned effect
#'   (per-condition output); \code{NULL} for the unconditioned fit.
#' @param trim Logical. Trim the retained fit to the fields needed downstream.
#'   Default \code{TRUE}.
#' @param fullFit Logical. Retain the full fit object on each entry. Default
#'   \code{FALSE}.
#' @param fullFitAlphaOnly Logical. When retaining the full fit, keep only the
#'   per-effect alpha matrix. Default \code{TRUE}.
#' @param includeAllCs Logical. Include all credible sets rather than only the
#'   top one. Default \code{FALSE}.
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
    coverage = NULL,
    secondaryCoverage = c(0.7, 0.5),
    signalCutoff = 0.1,
    otherQuantities = NULL,
    region = NULL,
    priorEffTol = 1e-9,
    minAbsCorr = 0.8,
    medianAbsCorr = NULL,
    csInput = NULL,
    conditionIdx = NULL,
    trim = TRUE,
    fullFit = FALSE,
    fullFitAlphaOnly = TRUE,
    includeAllCs = FALSE
) {
    p <- as.list(environment())
    fits <- fits[!map_lgl(fits, is.null)]
    if (length(fits) == 0) {
        abort("At least one fine-mapping fit must be supplied.")
    }
    if (is.null(names(fits)) || any(names(fits) == "")) {
        abort("fits must be a named list; names define method identity.")
    }
    .ppFitsCombine(.ppFitsPerMethod(fits, p))
}

# Post-process each method's fit once (buildTopLoci per fit); the per-method
# 22-column contributions are row-bound later into the single top_loci table.
.ppFitsPerMethod <- function(fits, p) {
    posts <- map(names(fits), .ppOneFit, fits = fits, p = p)
    names(posts) <- names(fits)
    posts
}

# Row-bind the per-method top_loci tables and drop them from the per-method
# entries; returns the final finemappingResults + combined top_loci.
.ppFitsCombine <- function(posts) {
    perMethod <- map(posts, "top_loci")
    perMethod <- perMethod[!map_lgl(perMethod, is.null)]
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
    coverage = NULL,
    secondaryCoverage = c(0.7, 0.5),
    signalCutoff = 0.1,
    otherQuantities = NULL,
    region = NULL,
    priorEffTol = 1e-9,
    trim = TRUE,
    minAbsCorr = 0.8,
    medianAbsCorr = NULL,
    conditionIdx = NULL,
    fullFit = FALSE,
    fullFitAlphaOnly = TRUE,
    includeAllCs = FALSE,
    csInput = c("X", "Xcorr", "fsusie")
) {
    csInput <- arg_match(csInput)
    p <- as.list(environment())
    variantNames <- extractVariantNames(fit)
    sumstats <- extractSumstats(fit, dataX, dataY, xScalar, yScalar, method)
    csTables <- .ppCsTables(p, csInput)
    # Always build the canonical unfiltered table; the FineMappingEntry stores
    # it as-is so accessors can filter by PIP at query time.
    topLociFull <- .ppTopLoci(p, csTables, variantNames, sumstats)
    # trim = TRUE stores a minimal subset of the fit; FALSE keeps the full
    # untrimmed susie return (mu / mu2 / lbf_variable / V / ...).
    storedFit <- if (isTRUE(trim)) {
        trimFinemappingFit(
            fit,
            selectEffects(fit, priorEffTol),
            method,
            csTables
        )
    } else {
        fit
    }
    fmEntry <- FineMappingEntry(
        variantIds = variantNames,
        susieFit = storedFit,
        topLoci = topLociFull
    )
    .ppAssembleRes(p, topLociFull, fmEntry, sumstats)
}

# Credible-set tables for the fit at the requested coverages.
.ppCsTables <- function(p, csInput) {
    computeCsTables(
        p$fit,
        dataX = p$dataX,
        coverage = p$coverage,
        secondaryCoverage = p$secondaryCoverage,
        method = p$method,
        csInput = csInput,
        minAbsCorr = p$minAbsCorr,
        medianAbsCorr = p$medianAbsCorr
    )
}

# Canonical unfiltered top-loci table (signalCutoff = 0).
.ppTopLoci <- function(p, csTables, variantNames, sumstats) {
    buildTopLoci(
        p$fit,
        csTables,
        variantNames = variantNames,
        sumstats = sumstats,
        af = p$af,
        method = p$method,
        signalCutoff = 0,
        dataY = p$dataY,
        otherQuantities = p$otherQuantities,
        region = p$region,
        conditionIdx = p$conditionIdx,
        fullFit = p$fullFit,
        fullFitAlphaOnly = p$fullFitAlphaOnly,
        includeAllCs = p$includeAllCs
    )
}

# Assemble the wrapper-facing result: PIP-filtered top_loci (legacy behaviour
# for non-S4 callers) + the entry + optional sumstats/sampleNames/context.
.ppAssembleRes <- function(p, topLociFull, fmEntry, sumstats) {
    topLociWrapper <- topLociFull
    if (
        !is.null(p$signalCutoff) &&
            p$signalCutoff > 0 &&
            nrow(topLociWrapper) > 0L
    ) {
        keep <- !is.na(topLociWrapper$pip) &
            topLociWrapper$pip > p$signalCutoff
        topLociWrapper <- topLociWrapper[keep, , drop = FALSE]
    }
    res <- list(
        top_loci = topLociWrapper,
        finemappingEntry = fmEntry,
        method = p$method
    )
    if (!is.null(sumstats)) {
        res$sumstats <- sumstats
    }
    sampleNames <- .sampleNamesFromDataY(p$dataY)
    if (!is.null(sampleNames)) {
        res$sampleNames <- sampleNames
    }
    if (p$method == "mvsusie" && !is.null(p$fit$outcome_names)) {
        res$contextNames <- p$fit$outcome_names
    }
    if (!is.null(p$otherQuantities)) {
        res$otherQuantities <- p$otherQuantities
    }
    res
}

extractVariantNames <- function(fit) {
    variantNames <- names(fit$pip)
    if (is.null(variantNames)) {
        variantNames <- colnames(fit$alpha)
    }
    if (is.null(variantNames)) {
        variantNames <- str_c("variant_", seq_along(fit$pip))
    }
    tryCatch(normalizeVariantId(variantNames), error = function(e) variantNames)
}

extractSumstats <- function(
    fit,
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
    sumstats$betahat <- sumstats$betahat * yScalar / xScalar
    sumstats$sebetahat <- sumstats$sebetahat * yScalar / xScalar
    sumstats
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
    primaryCoverage <- coverage
    if (is.null(primaryCoverage)) {
        primaryCoverage <- fit$sets$requested_coverage
    }
    if (is.null(primaryCoverage)) {
        primaryCoverage <- 0.95
    }
    coverages <- unique(c(primaryCoverage, secondaryCoverage))
    coverages <- coverages[!is.na(coverages)]

    tables <- map(
        coverages,
        .computeCsTableForCov,
        fit = fit,
        dataX = dataX,
        csInput = csInput,
        minAbsCorr = minAbsCorr,
        medianAbsCorr = medianAbsCorr
    )
    names(tables) <- map_chr(
        coverages,
        formatCsColumn,
        method = method
    )
    attr(tables, "coverage") <- coverages
    tables
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

# fSuSiE credible sets: purity is the min |correlation| WITHIN each CS
# (fsusieR::cal_purity), stamped as sets$purity$min.abs.corr for the canonical
# .csPurityVec() reader; cs_corr keeps the BETWEEN-CS correlation matrix.
.csTableFsusie <- function(fit, dataX, coverage) {
    sets <- tryCatch(
        fsusieGetCs(fit, dataX, requestedCoverage = coverage),
        error = function(e) list(cs = list(), requested_coverage = coverage)
    )
    if (
        is.null(sets$cs) ||
            length(sets$cs) == 0 ||
            all(map_lgl(sets$cs, is.null))
    ) {
        sets$cs <- list()
        return(list(sets = sets, pip = fit$pip))
    }
    if (requireNamespace("fsusieR", quietly = TRUE)) {
        purity <- tryCatch(
            as.numeric(unlist(fsusieR::cal_purity(sets$cs, dataX))),
            error = function(e) NULL
        )
        if (!is.null(purity) && length(purity) == length(sets$cs)) {
            sets$purity <- tibble(min.abs.corr = purity)
        }
    }
    list(sets = sets, pip = fit$pip)
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
    csArgs <- list(coverage = coverage)
    if (!is.null(minAbsCorr)) {
        csArgs$min_abs_corr <- minAbsCorr
    }
    if (!is.null(medianAbsCorr)) {
        csArgs$median_abs_corr <- medianAbsCorr
    }
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
    signature(x = "FineMappingEntry", ldSource = "SumStatsBase"),
    function(x, ldSource, ...) {
        fit <- getSusieFit(x)
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
            getVariantIds(x),
            label = "computeCsCorrelation",
            onMissing = "error"
        )
        get_cs_correlation(list(sets = fit$sets, pip = fit$pip), Xcorr = xcorr)
    }
)

# Individual-level LD source: genotypes -> aligned X. susie fits derive the LD
# via get_cs_correlation(X = ); fSuSiE fits (class "susiF") via cal_cor_cs().
#' @rdname computeCsCorrelation
setMethod(
    "computeCsCorrelation",
    signature(x = "FineMappingEntry", ldSource = "QtlDataset"),
    function(x, ldSource, ...) {
        fit <- getSusieFit(x)
        if (.csCountBelowTwo(fit)) {
            return(NULL)
        }
        geno <- .csGenotypesForFit(ldSource, getVariantIds(x))
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
)

#' @rdname computeCsCorrelation
setMethod(
    "computeCsCorrelation",
    signature(x = "FineMappingEntry", ldSource = "ANY"),
    function(x, ldSource, ...) {
        abort(glue(
            "computeCsCorrelation() requires a QtlDataset, QtlSumStats, or ",
            "GwasSumStats as `ldSource`: the between-credible-set correlation ",
            "is derived from that object's LD and is never stored on the fit."
        ))
    }
)

# Per-effect (per credible set) variant-level columns from the susie fit. Always
# returns `within_cs_pip` (the variant's alpha in the single effect of its
# assigned primary-coverage CS; NA for non-CS variants -- alpha is a
# probability,
# no scaling). With fullFit = TRUE it also widens the per-effect matrices, one
# column set per CS: `within_cs_pip_<lab>` (alpha) and -- unless
# fullFitAlphaOnly -- `cs_logbf_<lab>` (lbf_variable), `cs_effect_<lab>` (mu /
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
    fullFit = FALSE,
    fullFitAlphaOnly = TRUE,
    includeAllCs = FALSE
) {
    nV <- if (is.null(alpha) || length(dim(alpha)) < 2L) {
        length(primaryCsPos)
    } else {
        ncol(alpha)
    }
    hasAlpha <- !is.null(alpha) && length(dim(alpha)) == 2L && nrow(alpha) > 0L
    withinPip <- .ffcWithinPip(alpha, primaryCsPos, effectOf, nV, hasAlpha)
    cols <- tibble(within_cs_pip = withinPip)
    if (!isTRUE(fullFit) || !hasAlpha) {
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
        fullFitAlphaOnly,
        includeAllCs
    )
}

# Per-variant PIP within its primary-coverage credible set (NA outside any CS).
.ffcWithinPip <- function(alpha, primaryCsPos, effectOf, nV, hasAlpha) {
    withinPip <- rep(NA_real_, nV)
    if (!(hasAlpha && length(primaryCsPos) == nV && length(effectOf) > 0L)) {
        return(withinPip)
    }
    for (v in seq_len(nV)) {
        cp <- primaryCsPos[[v]]
        if (!is.na(cp) && cp >= 1L && cp <= length(effectOf)) {
            L <- effectOf[[cp]]
            if (!is.na(L) && L >= 1L && L <= nrow(alpha)) {
                withinPip[[v]] <- alpha[L, v]
            }
        }
    }
    withinPip
}

# Wide per-effect columns: within_cs_pip_<lab> always, plus cs_logbf / cs_effect
# / cs_effect_var when fullFitAlphaOnly is FALSE. Effects come from every effect
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
    fullFitAlphaOnly,
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
    for (i in seq_along(effs)) {
        L <- effs[[i]]
        lab <- labs[[i]]
        cols[[str_c("within_cs_pip_", lab)]] <- unname(alpha[L, ])
        if (!isTRUE(fullFitAlphaOnly)) {
            if (!is.null(lbfMat) && L <= nrow(lbfMat)) {
                cols[[str_c("cs_logbf_", lab)]] <- unname(lbfMat[L, ])
            }
            if (!is.null(mu) && L <= nrow(mu)) {
                cols[[str_c("cs_effect_", lab)]] <- unname(mu[L, ] / scale)
            }
            if (!is.null(mu) && !is.null(mu2) && L <= nrow(mu2)) {
                cols[[str_c("cs_effect_var_", lab)]] <-
                    unname((mu2[L, ] - mu[L, ]^2) / scale^2)
            }
        }
    }
    cols
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
.fmCsIdxAtCoverage <- function(targetCov, coverageValues, csTables, nV) {
    out <- integer(nV)
    hit <- which(abs(coverageValues - targetCov) < 1e-12)
    if (length(hit) == 0L) {
        return(out)
    }
    sets <- csTables[[hit[1L]]]$sets$cs
    if (is.null(sets) || length(sets) == 0L) {
        return(out)
    }
    for (csIdx in seq_along(sets)) {
        vi <- as.integer(sets[[csIdx]])
        vi <- vi[vi >= 1L & vi <= nV & out[vi] == 0L]
        out[vi] <- csIdx
    }
    out
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
#' @param method Method name (e.g. \code{"susie"}, \code{"susieInf"}). Required.
#' @param signalCutoff PIP cutoff for retaining PIP-only (non-CS) variants.
#' @param dataY Optional regional phenotype matrix; \code{nrow(dataY)} fills
#'   \code{n}, \code{colnames(dataY)[1]} fills \code{gene}.
#' @param otherQuantities Optional list. Default is NULL.
#' @param region Optional \code{"chr:start-end"} string. Default is NULL.
#' @param conditionIdx Integer or \code{NULL}. Index of the conditioned effect
#'   (per-condition output); \code{NULL} for the unconditioned fit.
#' @param fullFit Logical. Retain the full fit object on each entry. Default
#'   \code{FALSE}.
#' @param fullFitAlphaOnly Logical. When retaining the full fit, keep only the
#'   per-effect alpha matrix. Default \code{TRUE}.
#' @param includeAllCs Logical. Include all credible sets rather than only the
#'   top one. Default \code{FALSE}.
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
    method,
    signalCutoff = 0,
    dataY = NULL,
    otherQuantities = NULL,
    region = NULL,
    conditionIdx = NULL,
    fullFit = FALSE,
    fullFitAlphaOnly = TRUE,
    includeAllCs = FALSE
) {
    if (missing(method)) {
        method <- NULL
    }
    p <- as.list(environment())
    .btlValidateMethod(method)
    if (length(variantNames) == 0L) {
        return(.emptyTopLoci())
    }
    .btlBuild(p)
}

# Orchestrate the top-loci table from the captured argument list `p`.
.btlBuild <- function(p) {
    nV <- length(p$variantNames)
    cov <- .btlCoverage(p$csTables)
    fc <- .btlFitConstants(p$dataY, p$otherQuantities, p$region)
    post <- .btlPosterior(p$fit, p$conditionIdx, nV)
    marg <- .btlMarginal(p$sumstats, nV)
    cs <- .btlCsMembership(cov, p$csTables, nV)
    fullFitBlock <- .btlFullFitBlock(
        p$fit,
        post,
        cov,
        cs,
        p$csTables,
        nV,
        p[c("fullFit", "fullFitAlphaOnly", "includeAllCs")]
    )
    out <- .btlAssemble(
        p$variantNames,
        .btlParseVariants(p$variantNames),
        fc,
        marg,
        post,
        p$fit,
        p$af,
        p$method,
        .btlCsBlock(p$method, cs, nV),
        fullFitBlock,
        nV
    )
    cond <- .btlConditional(
        p$fit,
        p$method,
        p$conditionIdx,
        cov,
        cs$covSorted,
        p$csTables,
        nV
    )
    .btlFinalize(out, cond, p$conditionIdx, p$signalCutoff)
}

# Attach per-condition columns (multi-condition fits) and apply the PIP cutoff.
.btlFinalize <- function(out, cond, conditionIdx, signalCutoff) {
    if (!is.null(conditionIdx)) {
        out$conditional_effect <- cond$condEffect
        out$lfsr <- cond$condLfsr
    }
    if (!is.null(signalCutoff) && signalCutoff > 0) {
        out <- filter(out, !is.na(.data$pip) & .data$pip > signalCutoff)
    }
    out
}

# buildTopLoci step helpers ---------------------------------------------------

# `method` is required and must be a single non-empty, non-NA string.
.btlValidateMethod <- function(method) {
    if (
        is.null(method) ||
            length(method) != 1L ||
            is.na(method) ||
            str_length(method) == 0L
    ) {
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

# Per-fit constants: sample size, gene (first phenotype column), event id, and
# the parsed genomic range.
.btlFitConstants <- function(dataY, otherQuantities, region) {
    dataYMat <- if (!is.null(dataY)) as.matrix(dataY) else NULL
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
.btlParseVariants <- function(variantNames) {
    parsed <- tryCatch(
        suppressWarnings(parseVariantId(variantNames)),
        error = function(e) {
            eMsg <- conditionMessage(e)
            msg <- glue("buildTopLoci: parseVariantId failed: {eMsg}")
            abort(msg)
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
.btlCsMembership <- function(coverageValues, csTables, nV) {
    covSorted <- sort(
        unique(coverageValues[is.finite(coverageValues)]),
        decreasing = TRUE
    )
    csIdxByCov <- map(
        covSorted,
        .fmCsIdxAtCoverage,
        coverageValues,
        csTables,
        nV
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
        cm <- tryCatch(mvsusieR::coef.mvsusie(fit), error = function(e) NULL)
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
    for (csPos in seq_along(setsPrim)) {
        L <- effectOf[csPos]
        if (is.na(L) || L < 1L || L > dim(clf)[1L]) {
            next
        }
        vi <- as.integer(setsPrim[[csPos]])
        vi <- vi[vi >= 1L & vi <= nV]
        if (length(vi) > 0L) {
            condLfsr[vi] <- as.numeric(clf[L, vi, conditionIdx])
        }
    }
    condLfsr
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
            map(cs$csIdxByCov, .btlCsLabel, methodTag = methodTag),
            cs$csColNames
        ),
        set_names(cs$csPurityByCov, str_c(cs$csColNames, "_purity"))
    )
    if (length(csList) == 0L) {
        return(tibble(.rows = nV))
    }
    as_tibble(csList, .name_repair = "minimal")
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
    effectOfPrim <- integer(0)
    if (length(cs$covSorted) > 0L) {
        hP <- which(abs(coverageValues - cs$covSorted[1L]) < 1e-12)
        if (length(hP) > 0L) {
            spP <- csTables[[hP[1L]]]$sets$cs
            if (!is.null(spP) && length(spP) > 0L) {
                effectOfPrim <- suppressWarnings(
                    as.integer(str_remove(names(spP), "^L"))
                )
            }
        }
    }
    .fullFitColumns(
        post$alpha,
        post$mu,
        post$mu2,
        .asLbfMatrix(fit),
        fit$X_column_scale_factors,
        primaryCsPos,
        effectOfPrim,
        fullFit = opts$fullFit,
        fullFitAlphaOnly = opts$fullFitAlphaOnly,
        includeAllCs = opts$includeAllCs
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
        N = rep(fc$fitN, nV),
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

# Translate susieR's snake-case `sets$purity` columns into pecotmr camelCase.
# Accepts a data.frame, matrix, or NULL; preserves type and column order.
.translateSusiePurity <- function(p) {
    if (is.null(p)) {
        return(p)
    }
    lookup <- c(
        "min.abs.corr" = "minAbsCorr",
        "mean.abs.corr" = "meanAbsCorr",
        "median.abs.corr" = "medianAbsCorr"
    )
    if (is.data.frame(p)) {
        nm <- names(p)
        names(p) <- if_else(is_in(nm, names(lookup)), unname(lookup[nm]), nm)
    } else if (is.matrix(p)) {
        cn <- colnames(p)
        if (!is.null(cn)) {
            colnames(p) <- if_else(
                is_in(cn, names(lookup)),
                unname(lookup[cn]),
                cn
            )
        }
    }
    p
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
    pr <- tryCatch(parseRegion(as.character(regionStr)), error = function(e) {
        NULL
    })
    if (is.null(pr) || !is.data.frame(pr)) {
        return(c(start = NA_integer_, end = NA_integer_))
    }
    c(start = as.integer(pr$start), end = as.integer(pr$end))
}

# Project the new 22-column `top_loci` into the legacy shape expected by the
# FineMappingResult S4 slot, vcf_writer, getPIP, and getCS. We add backward-
# compatible aliases without renaming any column in the wrapper-facing
# `top_loci`:
#
#   * `variant_id` -- copy of `variant`
#   * `cs`         -- integer credible-set index derived from `cs_95` strings of
#                    the form `<method>_<idx>` (PIP-only `<method>_0` -> 0L)
#
# This isolates the schema change to susieWrapper.R so allClasses.R,
# allMethods.R, and vcfWriter.R do not have to change.
.topLociForS4Slot <- function(topLoci) {
    if (is.null(topLoci) || nrow(topLoci) == 0) {
        return(tibble(
            variant_id = character(0),
            method = character(0)
        ))
    }
    out <- topLoci
    if (is_in("variant", names(out)) && !is_in("variant_id", names(out))) {
        out$variant_id <- out$variant
    }
    if (is_in("cs_95", names(out)) && !is_in("cs", names(out))) {
        out$cs <- map_int(out$cs_95, .cs95ToIndex)
        out$cs[is.na(out$cs)] <- 0L
    }
    out
}

trimFinemappingFit <- function(fit, effectIdx, method, csTables) {
    trimmed <- .trimBaseFit(fit, effectIdx, csTables)
    trimmed <- .trimAddCommon(trimmed, fit, effectIdx)
    if (method == "mvsusie") {
        trimmed <- .trimAddMvsusie(trimmed, fit, effectIdx)
    }
    # fSuSiE: keep the precomputed variants x features TWAS weight matrix
    # (fsusieWeights output, attached as $coef before trimming) so downstream
    # TWAS can read it without the dropped wavelet slots.
    if (method == "fsusie" && !is.null(fit$coef)) {
        trimmed$coef <- fit$coef
    }
    class(trimmed) <- unique(c(method, "susie"))
    trimmed
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
    if (!is.null(fit$X_column_scale_factors)) {
        trimmed$X_column_scale_factors <- fit$X_column_scale_factors
    }
    if (!is.null(fit$mu)) {
        trimmed$mu <- if (length(dim(fit$mu)) == 3) {
            fit$mu[effectIdx, , , drop = FALSE]
        } else {
            fit$mu[effectIdx, , drop = FALSE]
        }
    }
    if (!is.null(fit$mu2)) {
        trimmed$mu2 <- if (length(dim(fit$mu2)) == 3) {
            fit$mu2[effectIdx, , , drop = FALSE]
        } else {
            fit$mu2[effectIdx, , drop = FALSE]
        }
    }
    if (!is.null(fit$theta)) {
        trimmed$theta <- fit$theta
    }
    if (!is.null(fit$omega_weights)) {
        trimmed$omega_weights <- fit$omega_weights
    }
    trimmed
}

# mvsusie-specific slots: per-effect mu2_diag, the coefficient matrix, and the
# conditional lfsr array.
.trimAddMvsusie <- function(trimmed, fit, effectIdx) {
    if (!is.null(fit$mu2_diag)) {
        trimmed$mu2_diag <- fit$mu2_diag[effectIdx, , , drop = FALSE]
    }
    if (requireNamespace("mvsusieR", quietly = TRUE)) {
        trimmed$coef <- mvsusieR::coef.mvsusie(fit)[-1, , drop = FALSE]
    }
    if (!is.null(fit$conditional_lfsr)) {
        trimmed$clfsr <- fit$conditional_lfsr[effectIdx, , , drop = FALSE]
    }
    trimmed
}

#' Format Fine-mapping Post-processing for Protocol Output
#'
#' Promotes the primary method's per-method post-processing payload to the root
#' level and attaches the unified \code{top_loci} table. The primary method's
#' bare \code{FineMappingEntry} appears at \code{$finemappingEntry}; wrap it
#' into a \code{FineMappingResult} collection at the pipeline level once (study,
#' context, trait, method) identity tags are known.
#'
#' @param post Output from \code{\link{postprocessFinemappingFits}}.
#' @param primaryMethod Method whose result should populate root-level fields.
#' @return A list with root-level fields including \code{finemappingEntry} (a
#'   bare \code{FineMappingEntry} S4 payload) and \code{top_loci}.
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
    c(which(susieOutput$pip >= signalCutoff), unlist(susieOutput$sets$cs)) |>
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

calPurity <- function(lCs, X, method = "min") {
    tt <- list()

    for (k in seq_along(lCs)) {
        csIndices <- unlist(lCs[[k]])
        if (method == "min") {
            if (length(csIndices) == 1) {
                tt[[k]] <- 1
            } else {
                x <- abs(computeLd(
                    X[, csIndices, drop = FALSE],
                    method = "sample"
                ))
                x[col(x) == row(x)] <- NA
                tt[[k]] <- min(x, na.rm = TRUE)
            }
        } else {
            if (length(csIndices) == 1) {
                tt[[k]] <- c(1, 1, 1)
            } else {
                x <- abs(computeLd(
                    X[, csIndices, drop = FALSE],
                    method = "sample"
                ))
                x[col(x) == row(x)] <- NA
                tt[[k]] <- c(
                    min(x, na.rm = TRUE),
                    mean(x, na.rm = TRUE),
                    median(x, na.rm = TRUE)
                )
            }
        }
    }

    return(tt)
}


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
#' fit <- getSusieFit(fsusieFineMappingExample$entry[[1]])
#' fsusieGetCs(fit)
#' @export
fsusieGetCs <- function(fsusieObj, X, requestedCoverage = 0.95) {
    # Create 'cs' set with names
    csNamed <- set_names(
        fsusieObj$cs,
        str_c("L", seq_along(fsusieObj$cs))
    )

    # Create 'purity' data frame
    purityDf <- bind_rows(
        map(calPurity(fsusieObj$cs, X = X, method = "susie"), .asDataFrameT)
    )
    colnames(purityDf) <- c("minAbsCorr", "meanAbsCorr", "medianAbsCorr")
    # Credible-set label as a `cs` column (was rownames; tibbles carry none).
    purityDf <- bind_cols(tibble(cs = names(csNamed)), purityDf)

    # Create 'coverage' without
    coverageVector <- numeric(length(fsusieObj$alpha))
    for (i in seq_along(fsusieObj$alpha)) {
        alphaI <- fsusieObj$alpha[[i]]
        csI <- fsusieObj$cs[[i]]
        coverageVector[i] <- sum(alphaI[csI])
    }

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
#' @param ... Additional arguments passed to the fsusie function.
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
    ...
) {
    if (!requireNamespace("fsusieR", quietly = TRUE)) {
        # nocov start
        msg <- glue(
            "To use this function, please install fsusieR: ",
            "https://github.com/stephenslab/fsusieR"
        )
        abort(msg)
        # nocov end
    }
    fsusieObj <- fsusieR::susiF(
        X = X,
        Y = Y,
        pos = pos,
        L = L,
        prior = prior,
        max_SNP_EM = maxSnpEm,
        cov_lev = covLev,
        min_purity = minPurity,
        max_scale = maxScale,
        ...
    )
    .fsusieWrapperPostprocess(fsusieObj, X, minPurity, covLev)
}

# Drop dummy credible sets below the purity threshold (else build sets + CS
# correlations), then reshape alpha (per-effect list) into a single data.frame.
.fsusieWrapperPostprocess <- function(fsusieObj, X, minPurity, covLev) {
    if (all(abs(as.numeric(fsusieObj$purity)) < minPurity)) {
        fsusieObj$cs <- list(NULL)
        fsusieObj$sets <- list(cs = list(NULL), requested_coverage = covLev)
    } else {
        fsusieObj$sets <- fsusieGetCs(fsusieObj, X, requestedCoverage = covLev)
    }
    fsusieObj$alpha <- bind_rows(map(fsusieObj$alpha, .asDataFrameT))
    fsusieObj
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
#' @param ... Additional arguments forwarded to \code{mvsusieR::mvsusie}.
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
#' @export
fitMvsusie <- function(X, Y, prior_variance, coverage = 0.95, ...) {
    mvsusieR::mvsusie(
        X = X,
        Y = Y,
        prior_variance = prior_variance,
        coverage = coverage,
        ...
    )
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
#' @param ... Additional arguments forwarded to \code{mvsusieR::mvsusie_rss}.
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
#' @export
fitMvsusieRss <- function(Z, R, N, prior_variance, coverage = 0.95, ...) {
    mvsusieR::mvsusie_rss(
        Z = Z,
        R = R,
        N = N,
        prior_variance = prior_variance,
        coverage = coverage,
        ...
    )
}

#' Fit fSuSiE on individual-level (X, Y, pos) data
#'
#' Thin wrapper around \code{fsusieR::susiF}.
#'
#' @param X Numeric matrix of genotypes (samples x variants).
#' @param Y Numeric matrix of multi-trait outcomes (samples x traits).
#' @param pos Numeric vector of trait positions (length \code{ncol(Y)}).
#' @param ... Additional arguments forwarded to \code{fsusieR::susiF}.
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
#' fitFsusie(X, Y, pos = pos, L = 2)
#' @export
fitFsusie <- function(X, Y, pos, ...) {
    fsusieR::susiF(X = X, Y = Y, pos = pos, ...)
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
    fit$intercept <- 0
    coef.susie(fit)[-1]
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
    y,
    requiredFields,
    token = "susie",
    userArgs = list(),
    retainFit = FALSE
) {
    if (is.null(fit)) {
        fit <- .fmFitSusieIndiv(X, y, token, userArgs = userArgs)
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
    if (all(is_in(requiredFields, names(fit)))) {
        weights <- .susieCoefWeights(fit)
    } else {
        weights <- rep(0, length(fit$pip))
    }
    if (retainFit) {
        attr(weights, "fit") <- fit
    }
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
#' @param retainFit If TRUE, stores the fitted object as an attribute on the
#'   returned weights.
#' @param ... Additional arguments passed to `susieR::susie()` when fitting.
#' @return Numeric vector of variant weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' susieWeights(X, y)
#' @export
susieWeights <- function(
    X = NULL,
    y = NULL,
    susieFit = NULL,
    retainFit = FALSE,
    ...
) {
    .susieExtractWeights(
        susieFit,
        X,
        y,
        requiredFields = c("alpha", "mu", "X_column_scale_factors"),
        token = "susie",
        userArgs = list(...),
        retainFit = retainFit
    )
}

#' Compute SuSiE-ASH TWAS weights
#'
#' Extracts coefficients from an existing SuSiE-ASH fit or fits
#' `susieR::susie()` with `unmappable_effects = "ash"`.
#'
#' @param X Genotype matrix. Required when `susieAshFit` is NULL.
#' @param y Phenotype vector. Required when `susieAshFit` is NULL.
#' @param susieAshFit Optional fitted SuSiE-ASH object.
#' @param retainFit If TRUE, stores the fitted object as an attribute on the
#'   returned weights.
#' @param ... Additional arguments passed to `susieR::susie()` when fitting.
#' @return Numeric vector of variant weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' susieAshWeights(X, y)
#' @export
susieAshWeights <- function(
    X = NULL,
    y = NULL,
    susieAshFit = NULL,
    retainFit = FALSE,
    ...
) {
    .susieExtractWeights(
        susieAshFit,
        X,
        y,
        requiredFields = c("alpha", "mu", "theta", "X_column_scale_factors"),
        token = "susieAsh",
        userArgs = list(...),
        retainFit = retainFit
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
#' @param X Genotype matrix. Required when `susieInfFit` is NULL.
#' @param y Phenotype vector. Required when `susieInfFit` is NULL.
#' @param susieInfFit Optional fitted SuSiE-inf object.
#' @param retainFit If TRUE, stores the fitted object as an attribute on the
#'   returned weights.
#' @param ... Additional arguments passed to `susieR::susie()` when fitting.
#' @return Numeric vector of variant weights.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' y <- eqtlRegionExample$yRes
#' susieInfWeights(X, y)
#' @export
susieInfWeights <- function(
    X = NULL,
    y = NULL,
    susieInfFit = NULL,
    retainFit = FALSE,
    ...
) {
    .susieExtractWeights(
        susieInfFit,
        X,
        y,
        requiredFields = c("alpha", "mu", "theta", "X_column_scale_factors"),
        token = "susieInf",
        userArgs = list(...),
        retainFit = retainFit
    )
}
# Internal helper: extract weights from a susieRss fit.
# Mirrors .susie_extract_weights but uses the RSS interface.
#' @importFrom susieR coef.susie susie_rss
#' @noRd
.susieRssExtractWeights <- function(
    fit,
    z,
    R,
    n,
    requiredFields,
    token = "susie",
    userArgs = list(),
    retainFit = FALSE
) {
    if (is.null(fit)) {
        fit <- .fmFitSusieRss(z, R, n, token, userArgs = userArgs)
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
    if (all(is_in(requiredFields, names(fit)))) {
        weights <- .susieCoefWeights(fit)
    } else {
        weights <- rep(0, length(fit$pip))
    }
    if (retainFit) {
        attr(weights, "fit") <- fit
    }
    return(weights)
}

#' Compute SuSiE-RSS TWAS weights
#'
#' Extracts coefficients from an existing SuSiE-RSS fit or fits
#' \code{susieR::susie_rss()} from summary statistics and LD.
#'
#' @param stat List with components \code{z} (z-scores), \code{n} (sample
#'   sizes).
#' @param LD LD correlation matrix.
#' @param susieRssFit Optional pre-fitted SuSiE-RSS object.
#' @param retainFit If TRUE, stores the fitted object as an attribute.
#' @param methodArgs Named list of additional arguments passed to
#'   \code{susieR::susie_rss()}. Use this instead of \code{...} to avoid partial
#'   matching of short argument names (e.g. \code{L}) to the \code{LD}
#'   parameter.
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
#' susieRssWeights(stat, LD)
#' @export
susieRssWeights <- function(
    stat,
    LD,
    susieRssFit = NULL,
    retainFit = TRUE,
    methodArgs = list()
) {
    .susieRssExtractWeights(
        fit = susieRssFit,
        z = stat$z,
        R = LD,
        n = median(stat$n),
        requiredFields = c("alpha", "mu", "X_column_scale_factors"),
        token = "susie",
        userArgs = methodArgs,
        retainFit = retainFit
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
#' susieInfRssWeights(stat, LD)
#' @export
susieInfRssWeights <- function(
    stat,
    LD,
    susieInfRssFit = NULL,
    retainFit = TRUE,
    methodArgs = list()
) {
    .susieRssExtractWeights(
        fit = susieInfRssFit,
        z = stat$z,
        R = LD,
        n = median(stat$n),
        requiredFields = c("alpha", "mu", "theta", "X_column_scale_factors"),
        token = "susieInf",
        userArgs = methodArgs,
        retainFit = retainFit
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
#' susieAshRssWeights(stat, LD)
#' @export
susieAshRssWeights <- function(
    stat,
    LD,
    susieAshRssFit = NULL,
    retainFit = TRUE,
    methodArgs = list()
) {
    .susieRssExtractWeights(
        fit = susieAshRssFit,
        z = stat$z,
        R = LD,
        n = median(stat$n),
        requiredFields = c("alpha", "mu", "theta", "X_column_scale_factors"),
        token = "susieAsh",
        userArgs = methodArgs,
        retainFit = retainFit
    )
}
#' Compute mvSuSiE TWAS weights
#'
#' Extracts coefficients from an existing mvSuSiE fit or fits `fitMvsusie()`
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
#' @param ... Additional arguments passed to `fitMvsusie()` when fitting.
#' @return Matrix of variant weights.
#' @examples
#' \donttest{
#' # Requires susieR attached (mvsusieR 0.3.0 packaging limitation).
#' library(susieR)
#' data(multiTraitData)
#' X <- multiTraitData$X[, 1:60]
#' Y <- multiTraitData$Y
#' mvsusieWeights(X = X, Y = Y, L = 5, LGreedy = 2)
#' }
#' @export
mvsusieWeights <- function(
    mvsusieFit = NULL,
    X = NULL,
    Y = NULL,
    priorVariance = NULL,
    residualVariance = NULL,
    L = 30,
    LGreedy = 5,
    verbose = FALSE,
    ...
) {
    if (!requireNamespace("mvsusieR", quietly = TRUE)) {
        # nocov start
        msg <- glue(
            "Package 'mvsusieR' is required. Install with: ",
            "devtools::install_github('stephenslab/mvsusieR')"
        )
        abort(msg)
        # nocov end
    }
    if (is.null(mvsusieFit)) {
        inform("mvsusieFit is not provided; fitting mvSuSiE now ...")
        if (is.null(X) || is.null(Y)) {
            abort("Both X and Y must be provided if mvsusieFit is NULL.")
        }
        if (is.null(priorVariance)) {
            priorVariance <- mvsusieR::create_mixture_prior(R = ncol(Y))
        }
        if (!is.null(LGreedy)) {
            LGreedy <- min(LGreedy, L)
        }

        mvsusieFit <- fitMvsusie(
            X = X,
            Y = Y,
            L = L,
            L_greedy = LGreedy,
            prior_variance = priorVariance,
            residual_variance = residualVariance,
            estimate_residual_variance = TRUE,
            verbose = verbose,
            ...
        )
    }
    return(mvsusieR::coef.mvsusie(mvsusieFit)[-1, ])
}

# One wavelet basis row: inverse-DWT (wr) of the unit coefficient vector e_k,
# using the fit's template DWT object.
# @noRd
.fmReconstructUnit <- function(k, nWac, scaleCols, template) {
    coeffRow <- numeric(nWac)
    coeffRow[k] <- 1
    temp <- template
    temp$D <- coeffRow[-scaleCols]
    temp$C[length(temp$C)] <- sum(coeffRow[scaleCols])
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
#' @param retainFit If TRUE, stores the fit as an attribute on the result.
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
    retainFit = FALSE
) {
    if (is.null(fsusieFit)) {
        msg <- glue(
            "fsusieWeights: `fsusieFit` is required. fSuSiE is functional ",
            "and cannot be refit from a bare (X, Y); fit it via ",
            "fineMappingPipeline() and pass the fitted fsusieR::susiF ",
            "object."
        )
        abort(msg)
    }
    fast <- .fsusieWeightsFastPath(fsusieFit, variantIds, retainFit)
    if (!is.null(fast)) {
        return(fast)
    }
    .fsusieWeightsRequire()
    fit <- fsusieFit
    .fsusieWeightsCheckSlots(fit)
    csdX <- as.numeric(fit$csd_X)
    alphaList <- .fsusieAlphaList(fit$alpha)
    S <- .fsusieSynthesisMatrix(fit$n_wac, .fsusieScaleCols(fit))
    W <- .fsusieComputeW(fit, alphaList, csdX, S)
    W <- .fsusieWeightsNames(
        W,
        fit,
        variantIds,
        featureNames,
        length(csdX),
        ncol(S)
    )
    if (retainFit) {
        attr(W, "fit") <- fit
    }
    W
}

# Fast path: a trimmed fit carries the precomputed variants x features weight
# matrix in `$coef` (fineMappingPipeline computes it eagerly while the full fit
# is in hand, since trimming drops fitted_wc/csd_X/...). NULL if not applicable.
.fsusieWeightsFastPath <- function(fsusieFit, variantIds, retainFit) {
    if (!(is.matrix(fsusieFit$coef) && is.null(fsusieFit$fitted_wc))) {
        return(NULL)
    }
    W <- fsusieFit$coef
    if (!is.null(variantIds) && length(variantIds) == nrow(W)) {
        rownames(W) <- variantIds
    }
    if (retainFit) {
        attr(W, "fit") <- fsusieFit
    }
    W
}

# fSuSiE weight reconstruction needs fsusieR + wavethresh.
.fsusieWeightsRequire <- function() {
    # nocov start
    if (!requireNamespace("fsusieR", quietly = TRUE)) {
        abort("Package 'fsusieR' is required for fsusieWeights().")
    }
    if (!requireNamespace("wavethresh", quietly = TRUE)) {
        abort("Package 'wavethresh' is required for fsusieWeights().")
    }
    # nocov end
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
.fsusieComputeW <- function(fit, alphaList, csdX, S) {
    invCsd <- 1 / csdX
    W <- matrix(0, nrow = length(csdX), ncol = ncol(S))
    for (l in seq_along(fit$fitted_wc)) {
        wc <- as.matrix(fit$fitted_wc[[l]])
        W <- W + (alphaList[[l]] * invCsd * wc) %*% S
    }
    W
}

# Attach variant (row) and feature/grid (column) names to the weight matrix.
.fsusieWeightsNames <- function(W, fit, variantIds, featureNames, p, nFeat) {
    rn <- variantIds
    if (is.null(rn)) {
        rn <- names(fit$csd_X)
    }
    if (is.null(rn)) {
        rn <- names(fit$pip)
    }
    if (!is.null(rn) && length(rn) == p) {
        rownames(W) <- rn
    }
    cn <- featureNames
    if (
        is.null(cn) &&
            !is.null(fit$outing_grid) &&
            length(fit$outing_grid) == nFeat
    ) {
        cn <- as.character(fit$outing_grid)
    }
    if (!is.null(cn) && length(cn) == nFeat) {
        colnames(W) <- cn
    }
    W
}
#' Compute mvSuSiE-RSS TWAS weights from summary statistics
#'
#' Multi-context summary-statistics analog of \code{\link{mvsusieWeights}}:
#' extracts coefficients from an existing \code{mvsusieR::mvsusie_rss} fit, or
#' fits one from \code{stat$z} (variants x conditions) and \code{LD}.
#'
#' Follows the \code{*_rss_weights(stat, LD, ...)} contract. Expects
#' \code{stat$z} to be a numeric matrix (variants x conditions) and
#' \code{stat$n} a per-context vector or scalar.
#'
#' @param stat A list with \code{z} (matrix variants x conditions) and \code{n}
#'   (numeric vector or scalar).
#' @param LD LD correlation matrix.
#' @param mvsusieRssFit Optional pre-fitted \code{mvsusieRss} object.
#' @param priorVariance Optional mvSuSiE prior variance specification. When
#'   NULL, \code{mvsusieR::create_mixture_prior()} is used with \code{R =
#'   ncol(stat$z)}.
#' @param residualVariance Optional residual covariance matrix.
#' @param L Maximum number of single effects (default 30).
#' @param LGreedy Initial greedy effect count (default 5).
#' @param retainFit If TRUE, attaches the fitted object as an attribute.
#' @param ... Additional arguments forwarded to \code{mvsusieR::mvsusie_rss}.
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
    priorVariance = NULL,
    residualVariance = NULL,
    L = 30,
    LGreedy = 5,
    retainFit = FALSE,
    ...
) {
    if (!requireNamespace("mvsusieR", quietly = TRUE)) {
        # nocov start
        msg <- glue(
            "Package 'mvsusieR' is required. ",
            "Install with: devtools::install_github('stephenslab/mvsusieR')"
        )
        abort(msg)
        # nocov end
    }
    if (is.null(mvsusieRssFit)) {
        mvsusieRssFit <- .mvsusieRssBuildFit(
            stat,
            LD,
            priorVariance,
            residualVariance,
            L,
            LGreedy,
            ...
        )
    }
    weights <- mvsusieR::coef.mvsusie(mvsusieRssFit)[-1, , drop = FALSE]
    if (retainFit) {
        attr(weights, "fit") <- mvsusieRssFit
    }
    weights
}

# Build the mvsusie-RSS fit from summary stats when the caller supplied none.
# @noRd
.mvsusieRssBuildFit <- function(
    stat,
    LD,
    priorVariance,
    residualVariance,
    L,
    LGreedy,
    ...
) {
    Z <- if (is.matrix(stat$z)) stat$z else as.matrix(stat$z)
    if (ncol(Z) < 2) {
        msg <- glue(
            "mvsusieRssWeights expects stat$z to have >= 2 columns ",
            "(one per context). For single-context use ",
            "susieRssWeights()."
        )
        abort(msg)
    }
    # mvsusieR::mvsusie_rss expects N to be a single scalar
    nScalar <- as.numeric(stats::median(stat$n))
    if (is.null(priorVariance)) {
        priorVariance <- mvsusieR::create_mixture_prior(R = ncol(Z))
    }
    if (!is.null(LGreedy)) {
        LGreedy <- min(LGreedy, L)
    }
    fitMvsusieRss(
        Z = Z,
        R = LD,
        N = nScalar,
        prior_variance = priorVariance,
        residual_variance = residualVariance,
        ...
    )
}

# =============================================================================
# Cross-condition credible-set merging
# =============================================================================

# Identify variant IDs that are associated with more than one credible set.
# @noRd
.identifyOverlapSets <- function(variantsSetsAndPipsList) {
    overlapSets <- list()
    for (variantId in names(variantsSetsAndPipsList)) {
        sets <- variantsSetsAndPipsList[[variantId]][["sets"]]
        if (length(sets) > 1) {
            overlapSets[[variantId]] <- sets
        }
    }
    return(overlapSets)
}

# Union-find root of `x` following the `parent` map.
# @noRd
.ufFindRoot <- function(x, parent) {
    while (!identical(parent[[x]], x)) {
        x <- parent[[x]]
    }
    x
}

# Union-find merge of `a` and `b` in `parent`; returns the updated parent map.
# @noRd
.ufUnion <- function(a, b, parent) {
    rootA <- .ufFindRoot(a, parent)
    rootB <- .ufFindRoot(b, parent)
    if (!identical(rootA, rootB)) {
        parent[[rootB]] <- rootA
    }
    parent
}

# Merge overlapping credible sets using connected components (union-find).
# @noRd
.mergeAndUpdateOverlapSets <- function(variantsSetsAndPipsList, overlapSets) {
    allSets <- unique(unlist(overlapSets))
    if (length(allSets) == 0) {
        return(list())
    }

    parent <- set_names(allSets, allSets)
    for (sets in overlapSets) {
        if (length(sets) > 1) {
            for (s in sets[-1]) {
                parent <- .ufUnion(sets[[1]], s, parent)
            }
        }
    }

    components <- split(
        names(parent),
        map_chr(names(parent), .ufFindRoot, parent)
    )
    setNameMap <- list()
    for (members in components) {
        label <- str_flatten(sort(members), ",")
        for (s in members) {
            setNameMap[[s]] <- label
        }
    }

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
    entries <- fineMappingResult$entry
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
        variant_id,
        .keep_all = TRUE
    )
    return(combinedTopLociDf)
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
    if (!is.null(chainFromInf) && token != "susieInf") {
        # SuSiE(-ash) initialised from a SuSiE-inf fit. userArgs are folded
        # into the
        # arg prep (not merged afterwards) so L_greedy is clamped to
        # min(#inf effects, L) rather than passed through raw.
        chainedArgs <- prepareSusieFromInfArgs(
            .fmMergeUserArgs(list(), token, userArgs),
            chainFromInf,
            refineDefault = if (token == "susie") TRUE else NULL,
            unmappableEffects = if (token == "susieAsh") "ash" else "none"
        )
        baseArgs <- modifyList(baseArgs, chainedArgs)
        baseArgs$X <- X
        baseArgs$y <- y
        baseArgs$coverage <- coverage
    } else {
        if (token == "susieInf") {
            baseArgs$convergence_method <- "pip"
            baseArgs$refine <- FALSE
            baseArgs$model_init <- NULL
        } else if (token == "susieAsh") {
            baseArgs$convergence_method <- "pip"
        }
        baseArgs <- .fmMergeUserArgs(baseArgs, token, userArgs)
    }
    fit <- exec(susieR::susie, !!!baseArgs)
    .setFinemappingFitClass(fit, token)
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
    baseArgs <- list(
        z = z,
        R = R,
        n = n,
        coverage = coverage,
        unmappable_effects = info$unmappableEffects
    )
    # rFinite = NULL removes the element -> susie_rss default; these sit in
    # baseArgs so they survive the chained modifyList / non-chained userArgs
    # merge, while user methodArgs (folded in after) still override them.
    baseArgs$R_finite <- rFinite
    baseArgs$R_mismatch <- rMismatch
    baseArgs <- .fmRssAddControl(baseArgs, rssControl)
    baseArgs <- if (!is.null(chainFromInf) && token != "susieInf") {
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
    .setFinemappingFitClass(exec(susieR::susie_rss, !!!baseArgs), "susieRss")
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
    baseArgs$control <- exec(susieR::susie_rss_control, !!!rssControl)
    baseArgs
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
    baseArgs <- modifyList(baseArgs, chainedArgs)
    baseArgs$z <- z
    baseArgs$R <- R
    baseArgs$n <- n
    baseArgs$coverage <- coverage
    baseArgs
}

# Non-chained fit: token-specific defaults then the user methodArgs merge.
.fmRssNonChainedArgs <- function(baseArgs, token, userArgs) {
    if (token == "susieInf") {
        baseArgs$convergence_method <- "pip"
        baseArgs$refine <- FALSE
        baseArgs$model_init <- NULL
    } else if (token == "susieAsh") {
        baseArgs$convergence_method <- "pip"
    }
    .fmMergeUserArgs(baseArgs, token, userArgs)
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
# returning a named list (token -> FineMappingEntry). Extracted from the
# univariate dispatch so the same logic serves the cis path (one block), the
# jointRegions=TRUE path (one concatenated block) and the jointRegions=FALSE
# path (one block per region, merged afterwards via .fmMergeEntries).
.fmFitXBlock <- function(
    X,
    y,
    toRun,
    addSusieInf,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    methodArgs,
    verbose,
    ctx,
    tid,
    cvFolds = 0,
    cvThreads = 1,
    samplePartition = NULL,
    af = NULL,
    fullFit = FALSE,
    fullFitAlphaOnly = TRUE,
    includeAllCs = FALSE
) {
    p <- as.list(environment())
    chainLocal <- .fmResolveSusieChain(toRun, addSusieInf)
    infFit <- .fmXInfFit(p, chainLocal)
    out <- list()
    for (tk in toRun) {
        fit <- .fmXFitOne(tk, p, chainLocal, infFit)
        if (is.null(fit)) {
            next
        }
        out[[tk]] <- .fmXPostprocess(fit, tk, p)
    }
    .fmXCrossValidate(out, p)
}

# Fit the shared susieInf model once, if the requested chain needs it.
.fmXInfFit <- function(p, chainLocal) {
    if (!chainLocal$runInf) {
        return(NULL)
    }
    if (p$verbose >= 1) {
        msg <- glue(
            "Fitting susieInf for (context='{p$ctx}', trait='{p$tid}') ..."
        )
        inform(msg)
    }
    .fmFitSusieIndiv(
        p$X,
        p$y,
        "susieInf",
        coverage = p$coverage,
        userArgs = p$methodArgs[["susieInf"]]
    )
}

# Resolve the fit for one method token; NULL means "skip this token".
.fmXFitOne <- function(tk, p, chainLocal, infFit) {
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
    if (p$verbose >= 1) {
        msg <- glue(
            "Fitting {tk} for (context='{p$ctx}', trait='{p$tid}') ..."
        )
        inform(msg)
    }
    .fmFitSusieIndiv(
        p$X,
        p$y,
        tk,
        chainFromInf = chainFrom,
        coverage = p$coverage,
        userArgs = p$methodArgs[[tk]]
    )
}

# Post-process one individual-level fit into a finemapping entry.
.fmXPostprocess <- function(fit, tk, p) {
    .fmPostprocessOne(
        fit = fit,
        method = tk,
        dataX = p$X,
        dataY = p$y,
        coverage = p$coverage,
        secondaryCoverage = p$secondaryCoverage,
        signalCutoff = p$signalCutoff,
        minAbsCorr = p$minAbsCorr,
        af = p$af,
        csInput = "X",
        fullFit = p$fullFit,
        fullFitAlphaOnly = p$fullFitAlphaOnly,
        includeAllCs = p$includeAllCs
    )
}

# Per-fold cross-validation across the fitted methods; attach each method's
# out-of-fold predictions to its entry.
.fmXCrossValidate <- function(out, p) {
    if (!(p$cvFolds > 1L && length(out) > 0L)) {
        return(out)
    }
    if (p$verbose >= 1) {
        msg <- glue(
            "Cross-validating ({p$cvFolds} folds) for ",
            "(context='{p$ctx}', trait='{p$tid}') ..."
        )
        inform(msg)
    }
    cv <- .fmWeightsCv(
        p$X,
        p$y,
        names(out),
        p$methodArgs,
        p$cvFolds,
        samplePartition = p$samplePartition,
        coverage = p$coverage,
        verbose = p$verbose,
        numThreads = p$cvThreads
    )
    for (tk in names(out)) {
        out[[tk]] <- .fmAttachCv(out[[tk]], .fmSliceCv(cv, tk))
    }
    out
}

# Fit every requested RSS token on one (z, R, n) sumstat block, returning a
# named list (token -> FineMappingEntry). The sumstat analog of .fmFitXBlock:
# the QtlSumStats and GwasSumStats methods both call it and differ only in how
# they push the returned entries (tuple shape) and the progress `label`.
.fmFitRssBlock <- function(
    z,
    R,
    n,
    toRun,
    addSusieInf,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    methodArgs,
    verbose,
    label,
    af = NULL,
    fullFit = FALSE,
    fullFitAlphaOnly = TRUE,
    includeAllCs = FALSE,
    serFallback = FALSE,
    rFinite = NULL,
    rMismatch = "none",
    rssControl = NULL,
    keepFullFit = "fallback"
) {
    p <- as.list(environment())
    chainLocal <- .fmResolveSusieChain(toRun, addSusieInf)
    infFit <- .fmRssInfFit(p, chainLocal)
    out <- list()
    for (tk in toRun) {
        f <- .fmRssFitOne(tk, p, chainLocal, infFit)
        if (is.null(f)) {
            next
        }
        ent <- .fmRssPostprocess(f$fit, p)
        if (f$isStd && isTRUE(p$serFallback)) {
            ent <- .fmRssRecordFallback(ent, f, p$keepFullFit)
        }
        out[[tk]] <- ent
    }
    out
}

# Fit the shared susieInf (RSS) model once, if the requested chain needs it.
.fmRssInfFit <- function(p, chainLocal) {
    if (!chainLocal$runInf) {
        return(NULL)
    }
    if (p$verbose >= 1) {
        msg <- glue("Fitting susieInf (RSS) for {p$label} ...")
        inform(msg)
    }
    .fmFitSusieRss(
        p$z,
        p$R,
        p$n,
        "susieInf",
        coverage = p$coverage,
        userArgs = p$methodArgs[["susieInf"]],
        rFinite = p$rFinite,
        rMismatch = p$rMismatch,
        rssControl = p$rssControl
    )
}

# Standard multi-effect SuSiE-RSS fit (susie / susieAsh): the only branch that
# carries susieR's finite-sample R diagnostics and honours the SER fallback.
.fmRssFitStd <- function(tk, p, chainLocal, infFit) {
    chainFrom <- if (
        (tk == "susie" && chainLocal$chainSusie) ||
            (tk == "susieAsh" && chainLocal$chainAsh)
    ) {
        infFit
    } else {
        NULL
    }
    if (p$verbose >= 1) {
        msg <- glue("Fitting {tk} (RSS) for {p$label} ...")
        inform(msg)
    }
    fit <- .fmFitSusieRss(
        p$z,
        p$R,
        p$n,
        tk,
        chainFromInf = chainFrom,
        coverage = p$coverage,
        userArgs = p$methodArgs[[tk]],
        rFinite = p$rFinite,
        rMismatch = p$rMismatch,
        rssControl = p$rssControl
    )
    rfd <- fit$R_finite_diagnostics
    flag <- if (!is.null(rfd) && !is.null(rfd$R_reliability_flag)) {
        isTRUE(rfd$R_reliability_flag)
    } else {
        NA
    }
    multiFit <- NULL
    if (isTRUE(p$serFallback) && isTRUE(flag) && !is.null(rfd$ser_model)) {
        multiFit <- fit
        fit <- .setFinemappingFitClass(rfd$ser_model, "susieRss")
    }
    list(fit = fit, flag = flag, multiFit = multiFit)
}

# Resolve the fit for one method token; NULL means "skip this token".
.fmRssFitOne <- function(tk, p, chainLocal, infFit) {
    if (tk == "susieInf") {
        if (!chainLocal$keepInf) {
            return(NULL)
        }
        return(list(fit = infFit, flag = NA, isStd = FALSE, multiFit = NULL))
    }
    if (tk == "ser") {
        if (p$verbose >= 1) {
            msg <- glue("Fitting ser (RSS single-effect) for {p$label} ...")
            inform(msg)
        }
        fit <- .fmFitSusieSer(
            p$z,
            p$n,
            coverage = p$coverage,
            userArgs = p$methodArgs[["ser"]]
        )
        return(list(fit = fit, flag = NA, isStd = FALSE, multiFit = NULL))
    }
    std <- .fmRssFitStd(tk, p, chainLocal, infFit)
    list(fit = std$fit, flag = std$flag, isStd = TRUE, multiFit = std$multiFit)
}

# Post-process one RSS fit into a finemapping entry.
.fmRssPostprocess <- function(fit, p) {
    .fmPostprocessOne(
        fit = fit,
        method = "susieRss",
        dataX = p$R,
        dataY = list(z = p$z),
        coverage = p$coverage,
        secondaryCoverage = p$secondaryCoverage,
        signalCutoff = p$signalCutoff,
        minAbsCorr = p$minAbsCorr,
        af = p$af,
        csInput = "Xcorr",
        fullFit = p$fullFit,
        fullFitAlphaOnly = p$fullFitAlphaOnly,
        includeAllCs = p$includeAllCs
    )
}

# Record the SER-fallback reliability decision (and retained multi-effect fit)
# on the entry's susieFit list. Gated on serFallback so the default path is
# byte-identical.
.fmRssRecordFallback <- function(ent, f, keepFullFit) {
    sf <- ent@susieFit
    sf$R_reliability_flag <- f$flag
    sf$serFallbackUsed <- isTRUE(f$flag)
    if (!is.null(f$multiFit) && is_in(keepFullFit, c("fallback", "all"))) {
        sf$multiEffectFit <- f$multiFit
    } else if (identical(keepFullFit, "all")) {
        sf$multiEffectFit <- f$fit
    }
    ent@susieFit <- sf
    ent
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
.ppOneFit <- function(method, fits, p) {
    fit <- .setFinemappingFitClass(fits[[method]], method)
    postprocessFinemappingFit(
        fit,
        method = method,
        dataX = p$dataX,
        dataY = p$dataY,
        xScalar = p$xScalar,
        yScalar = p$yScalar,
        af = p$af,
        coverage = p$coverage,
        secondaryCoverage = p$secondaryCoverage,
        signalCutoff = p$signalCutoff,
        otherQuantities = p$otherQuantities,
        region = p$region,
        priorEffTol = p$priorEffTol,
        minAbsCorr = p$minAbsCorr,
        medianAbsCorr = p$medianAbsCorr,
        csInput = p$csInput,
        conditionIdx = p$conditionIdx,
        trim = p$trim,
        fullFit = p$fullFit,
        fullFitAlphaOnly = p$fullFitAlphaOnly,
        includeAllCs = p$includeAllCs
    )
}

# @noRd
.ppDropTopLoci <- function(x) {
    x$top_loci <- NULL
    x
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
            sort(unique(unlist(extractedResult[[variantId]]$sets))),
            ","
        )
    }
    tibble(
        variant_id = variantId,
        credibleSetNames = credibleSetNames,
        maxPip = max(unlist(extractedResult[[variantId]]$pips)),
        medianPip = median(unlist(extractedResult[[variantId]]$pips))
    )
}

# (variant, CS) rows for the i-th entry, labelled cs_<entry>_<set>.
# @noRd
.extractCsEntryRows <- function(i, entries, csCol) {
    topLoci <- .translateLegacyTopLociCsColumns(getTopLoci(entries[[i]]))
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
    setNum <- unique(csIdx)
    setNum <- setNum[!is.na(setNum) & setNum != 0]
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
        select(variant_id, pip = all_of(pipCol)) |>
        mutate(set_name = str_c("cs_", i, "_", sn))
}

# @noRd
.csSplitToList <- function(df) {
    list(sets = df$set_name, pips = df$pip)
}
