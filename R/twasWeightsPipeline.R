# =============================================================================
# Helpers + S4 dispatch surface for twasWeightsPipeline
# =============================================================================

# Concatenate two TwasWeights collections row-wise, carrying forward every
# column (delegates to the generic `.rbindCollections`, which unions columns
# and pads a side lacking an optional column such as joint* / region).
# @noRd
#' @importFrom checkmate assertClass
.rbindTwasWeights <- function(a, b, ldReference = NULL) {
    assertClass(a, "TwasWeights")
    assertClass(b, "TwasWeights")
    # Carry forward every column (joint*, region, ...) and reconcile the
    # collection-level slots via the shared combine.
    .combineTupleCollections(list(a, b), ldReference, ".rbindTwasWeights")
}

# Normalize combine() varargs: accept either N objects or a single list of
# them; drop NULLs; require at least one input of the expected class `cls`.
.asCombineList <- function(parts, cls, fn) {
    # A single list argument is the collection itself, not a one-element
    # variadic call.
    unwrapped <- if (
        length(parts) == 1L &&
            is.list(parts[[1L]]) &&
            !methods::is(parts[[1L]], cls)
    ) {
        parts[[1L]]
    } else {
        parts
    }
    present <- compact(unwrapped)
    if (length(present) == 0L) {
        msg <- glue("{fn}: nothing to combine (need at least one {cls}).")
        abort(msg)
    }
    if (!all(map_lgl(present, methods::is, cls))) {
        msg <- glue("{fn}: every input must be a {cls}.")
        abort(msg)
    }
    present
}

#' Combine TwasWeights collections
#'
#' Row-bind two or more \code{\link{TwasWeights}} collections into one -- e.g.
#' assembling per-gene weight sets into a single per-region collection for
#' cTWAS. Joint-specification metadata columns are carried through.
#'
#' @param ... Two or more \code{TwasWeights} objects, or a single \code{list} of
#'   them.
#' @param ldReference Optional LD reference -- a genotype panel (see
#'   \code{\link{readGenotypes}}), a \code{GenotypeMixture} or an
#'   \code{LdMatrix} -- that every row of the combined collection uses instead
#'   of its own. Default \code{NULL}: each row keeps the reference it had,
#'   and views of one source are stored once. Applied when combining two or
#'   more inputs; a single input is returned unchanged.
#' @return A single combined \code{TwasWeights}.
#' @seealso \code{\link{combineFineMappingResults}}
#' @examples
#' twe <- TwasWeightsRow(
#'   variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)), weights = rep(0.1, 4),
#'   cvResult = list(rsq = 0.5), weightStandardized = FALSE)
#' tw1 <- TwasWeights(studyName = "s1", context = "brain", trait = "g1",
#'   method = "susie", entry = list(twe))
#' tw2 <- TwasWeights(studyName = "s2", context = "brain", trait = "g1",
#'   method = "susie", entry = list(twe))
#' combineTwasWeights(tw1, tw2)
#' @importFrom stringr str_starts
#' @export
combineTwasWeights <- function(..., ldReference = NULL) {
    parts <- .asCombineList(list(...), "TwasWeights", "combineTwasWeights")
    .combineTupleCollections(parts, ldReference, "combineTwasWeights")
}

# --- Multi-region (jointRegions) helpers for the QtlDataset method ----------

# Label a region block for per-region reporting: the genomic coordinate of a
# single-range window, or "cis" for the trait-derived (region = NULL) block.
.twasRegionLabel <- function(rg) {
    if (is.null(rg)) {
        return("cis")
    }
    str_c(
        as.character(GenomicRanges::seqnames(rg))[[1L]],
        ":",
        GenomicRanges::start(rg)[[1L]],
        "-",
        GenomicRanges::end(rg)[[1L]]
    )
}

# Select the per-region fine-mapping fits for region block `i`. A
# jointRegions=FALSE multi-region fine-mapping stores its per-region SuSiE fits
# as a named list (region1, region2, ...); pick the matching element. With a
# single block the fits are returned unchanged; a non-region-list fit under
# multiple blocks cannot be aligned and is dropped (the method learns fresh).
.twasFitsForRegion <- function(fits, i, nBlocks) {
    if (length(fits) == 0L || nBlocks == 1L) {
        return(fits)
    }
    compact(map(fits, .twasRegionFitOf, i = i))
}

# Flat per-region cvResult reporting table: one row per region carrying the
# region label plus that region's CV metric columns. Per-sample predictions are
# intentionally omitted -- this is a summary-reporting structure.
.twasRegionCvDf <- function(entries, regionLabels) {
    rows <- compact(map2(entries, regionLabels, .twasCvRow))
    if (length(rows) == 0L) {
        return(NULL)
    }
    bind_rows(rows)
}

# Concatenate one method's per-region TwasWeightsRow payloads into a single
# entry. Variants/weights are stacked (regions are disjoint), the per-region
# fits are kept as a named list, and cvResult becomes the flat per-region
# reporting data.frame.
.twasMergeRegionEntries <- function(entries, regionLabels) {
    keep <- !map_lgl(entries, is.null)
    fitted <- entries[keep]
    fittedLabels <- regionLabels[keep]
    if (length(fitted) == 0L) {
        return(NULL)
    }
    if (length(fitted) == 1L) {
        return(fitted[[1L]])
    }
    payloads <- map(fitted, .asTwRowPayload)
    wList <- map(payloads, weights)
    stacked <- if (is.matrix(wList[[1L]])) {
        exec(rbind, !!!wList)
    } else {
        unname(list_c(wList))
    }
    TwasWeightsRow(
        variantIds = unname(list_c(map(payloads, .twrPartsVariantIds))),
        weights = stacked,
        methodFits = set_names(map(payloads, methodFits), fittedLabels),
        cvResult = .twasRegionCvDf(payloads, fittedLabels),
        weightStandardized = weightStandardized(payloads[[1L]]),
        weightsDataType = weightsDataType(payloads[[1L]])
    )
}

# Unpack a MashPrior input into the internal twasWeightsPipeline arguments:
# $fullPrior -> mr.mash full-data dataDrivenPriorMatrices
# $dataDrivenPriorMatricesCv -> per-fold priors for twasWeightsCv
# $samplePartition -> the CV folds (an explicit `samplePartition` arg wins;
# otherwise the partition the per-fold priors were computed on) NULL input
# returns all-NULL, preserving the supplied samplePartition.
# @noRd
#' @importFrom checkmate assertClass
.unpackMashPrior <- function(mashPrior, samplePartition = NULL) {
    if (is.null(mashPrior)) {
        return(list(
            fullPrior = NULL,
            dataDrivenPriorMatricesCv = NULL,
            samplePartition = samplePartition
        ))
    }
    assertClass(mashPrior, "MashPrior")
    cvFits <- cvFits(mashPrior)
    perFold <- if (!is.null(cvFits)) cvFits$perFoldFits else NULL
    sp <- samplePartition %||% cvFits$samplePartition
    list(
        fullPrior = fullFit(mashPrior),
        dataDrivenPriorMatricesCv = perFold,
        samplePartition = sp
    )
}

# User-facing TWAS method tokens are unified across input classes;
# auto-dispatch picks the individual-level vs sumstat implementation based
# on the QtlDataset / QtlSumStats input. Each entry records:
#   individualImpl  Function name to call on QtlDataset input (NULL = not
#                   supported on individual-level input).
#   sumstatImpl     Function name to call on QtlSumStats input (NULL = not
#                   supported on sumstat input).
#   multivariate    Whether the method requires multi-trait / multi-context
#                   structure (mvsusie / mrmash / mvsusieRss / mrmashRss).
#
# Per the design: BGLR / qgg "Bayes alphabet" methods (bayes_a/b/c/l/n/r,
# b_lasso) are individual-only until the qgg CRAN release adds qBayes
# sumstat support. dpr_gibbs has the SDPR sumstat counterpart;
# dpr_vb / dpr_adaptive_gibbs remain individual-only. lasso, enet and ridge
# are glmnet on individual data and the penalized RSS solver on summary
# statistics, at alpha = 1, 0.5 (or the caller's) and 0. prsCs has no
# individual-level counterpart (it is a sumstat-only Bayesian shrinkage
# method).
# @noRd
.twasMethodCapabilities <- list(
    # NOTE: fine-mapping methods (susie / susieInf / susieAsh / mvsusie /
    # fsusie) are NOT listed here. Their availability is governed by
    # .fineMappingMethodCapabilities (the same registry fineMappingPipeline
    # uses) and gated by .twasCheckFineMappingMethods, which delegates
    # input-class compatibility to .fmCheckMethodCapabilities.
    mrash = list(
        individualImpl = "mrashWeights",
        sumstatImpl = "mrashRssWeights",
        multivariate = FALSE
    ),
    lasso = list(
        individualImpl = "lassoWeights",
        sumstatImpl = "lassosumRssWeights",
        multivariate = FALSE
    ),
    scad = list(
        individualImpl = "scadWeights",
        sumstatImpl = "scadRssWeights",
        multivariate = FALSE
    ),
    mcp = list(
        individualImpl = "mcpWeights",
        sumstatImpl = "mcpRssWeights",
        multivariate = FALSE
    ),
    l0learn = list(
        individualImpl = "l0learnWeights",
        sumstatImpl = "l0learnRssWeights",
        multivariate = FALSE
    ),
    mrmash = list(
        individualImpl = "mrmashWeights",
        sumstatImpl = "mrmashRssWeights",
        multivariate = TRUE
    ),
    dprGibbs = list(
        individualImpl = "dprGibbsWeights",
        sumstatImpl = "sdprWeights",
        multivariate = FALSE
    ),
    enet = list(
        individualImpl = "enetWeights",
        sumstatImpl = "enetRssWeights",
        multivariate = FALSE
    ),
    ridge = list(
        individualImpl = "ridgeWeights",
        sumstatImpl = "ridgeRssWeights",
        multivariate = FALSE
    ),
    # Individual-only DPR variants (sumstat counterparts not implemented).
    dprVb = list(
        individualImpl = "dprVbWeights",
        sumstatImpl = NULL,
        multivariate = FALSE
    ),
    dprAdaptiveGibbs = list(
        individualImpl = "dprAdaptiveGibbsWeights",
        sumstatImpl = NULL,
        multivariate = FALSE
    ),
    # qgg Bayes alphabet -- individual-only until qgg CRAN release.
    bayesA = list(
        individualImpl = "bayesAWeights",
        sumstatImpl = NULL,
        multivariate = FALSE
    ),
    bayesB = list(
        individualImpl = "bayesBWeights",
        sumstatImpl = NULL,
        multivariate = FALSE
    ),
    bayesC = list(
        individualImpl = "bayesCWeights",
        sumstatImpl = NULL,
        multivariate = FALSE
    ),
    bayesL = list(
        individualImpl = "bLassoWeights",
        sumstatImpl = NULL,
        multivariate = FALSE
    ),
    bayesN = list(
        individualImpl = "bayesNWeights",
        sumstatImpl = NULL,
        multivariate = FALSE
    ),
    bayesR = list(
        individualImpl = "bayesRWeights",
        sumstatImpl = NULL,
        multivariate = FALSE
    ),
    bLasso = list(
        individualImpl = "bLassoWeights",
        sumstatImpl = NULL,
        multivariate = FALSE
    ),
    # Sumstat-only Bayesian shrinkage (no individual-level analogue).
    prsCs = list(
        individualImpl = NULL,
        sumstatImpl = "prsCsWeights",
        multivariate = FALSE
    )
)

# Normalize a user-supplied `methods` argument (character vector, preset
# string, or named list per `.twasMethodLookup`) into a (token, args) pair
# suitable for `.twasWeightsPipelineMatrix` / the sumstat sub-pipelines.
# Returns a list with `tokens` (canonical short names, used for capability
# lookup) and `methodList` (the `<token>_weights = args` list passed to
# `learnTwasWeights` / sumstat helpers).
# @noRd
.twasNormalizeMethods <- function(methods, inputKind = "QtlDataset") {
    if (is.null(methods)) {
        methodList <- .twasMethodLookup("default")
        return(list(
            tokens = .twasTokensFromMethodList(methodList),
            methodList = methodList
        ))
    }
    if (is.character(methods)) {
        return(.twasNormalizeCharMethods(methods))
    }
    if (is(methods, "MethodsSelectionParam") || is.list(methods)) {
        # Routable here, unlike in the constructor: this call knows its
        # input class, so a plain list of overrides goes into that path's
        # slot and is checked against the engine that receives it.
        resolved <- .methodsParamResolve(
            .methodsParamFor(
                methods,
                inputKind,
                "TwasWeightsMethodsParam",
                "twasWeightsPipeline"
            ),
            inputKind
        )
        return(.twasNormalizeListMethods(resolved$methodArgs))
    }
    msg <- glue(
        "`methods` must be a character vector, a preset string, a named ",
        "list of per-method options, or a TwasWeightsMethodsParam() record."
    )
    abort(msg)
}

# Character `methods`: resolve regular tokens via .twasMethodLookup and append
# empty stub entries for fine-mapping tokens with no learner counterpart (e.g.
# fsusie) so the downstream gate can produce a method-specific error rather than
# the generic "unknown method token(s)".
# @noRd
.twasNormalizeCharMethods <- function(methods) {
    fmExtra <- setdiff(
        intersect(methods, .twasFineMappingTokens()),
        .twasKnownMethodLookupNames()
    )
    regular <- setdiff(methods, fmExtra)
    # Fine-mapping tokens carry no learner arguments of their own, so each
    # gets an empty stub entry under its `<token>_weights` key.
    methodList <- c(
        if (length(regular) > 0L) .twasMethodLookup(regular) else list(),
        .twasEmptyMethodArgs(str_c(fmExtra, "_weights"))
    )
    # Tokens come from the user input (canonical camelCase) -- the snake keys in
    # methodList are an internal detail of learnTwasWeights.
    list(tokens = unique(methods), methodList = methodList)
}

# One `methods` entry re-keyed to its canonical name, with the caller's
# kwargs merged over the method's defaults. A token with no learner default
# is kept under its own name (the capability gate reports it downstream).
# @noRd
.twasNormalizeOneMethod <- function(tk, methods) {
    base <- try_fetch(.twasMethodLookup(tk), error = function(cnd) NULL)
    if (is.null(base)) {
        return(list(key = tk, args = methods[[tk]]))
    }
    snake <- names(base)[[1L]]
    merged <- list_modify(base[[snake]], !!!compact(methods[[tk]]))
    # `impl` tells .resolveMethodFunction which function backs the token;
    # without it the bare token is used, which is not a function.
    list(
        key = snake,
        args = `attr<-`(merged, "impl", attr(base[[snake]], "impl"))
    )
}

# Named-list `methods`: re-key each entry to its canonical `<token>_weights`
# name, merge the caller's kwargs over the method's defaults, and carry the
# token->impl map as an "impl" attribute (without this, downstream
# .resolveMethodFunction falls back to the bare token, which is not a function).
# @noRd
.twasNormalizeListMethods <- function(methods) {
    entries <- map(names(methods), .twasNormalizeOneMethod, methods = methods)
    # Later entries win on a repeated canonical key, as the keyed assignment
    # in the loop did.
    keyed <- set_names(map(entries, "args"), map_chr(entries, "key"))
    methodList <- keyed[!duplicated(names(keyed), fromLast = TRUE)]
    list(
        tokens = .twasTokensFromMethodList(methodList),
        methodList = methodList
    )
}

# Canonical (camelCase) tokens known to .twasMethodLookup, for use by
# .twasNormalizeMethods. Source of truth: the methodMap inside
# .twasMethodLookup.
# @noRd
.twasKnownMethodLookupNames <- function() {
    c(
        "susie",
        "susieAsh",
        "susieInf",
        "mrash",
        "enet",
        "ridge",
        "lasso",
        "bayesR",
        "bayesL",
        "bayesA",
        "bayesB",
        "bayesC",
        "bayesN",
        "bLasso",
        "dprVb",
        "dprGibbs",
        "dprAdaptiveGibbs",
        "scad",
        "mcp",
        "l0learn",
        "mvsusie",
        "mrmash"
    )
}

# Convert a methodList (snake_case keys like `susie_inf_weights`) back to
# canonical camelCase tokens (susieInf). Falls back to the snake form for
# unknown keys.
# @noRd
.twasTokensFromMethodList <- function(methodList) {
    snake <- str_remove(names(methodList), "(_weights|Weights)$")
    snakeToCanonical <- c(
        susie = "susie",
        susie_ash = "susieAsh",
        susie_inf = "susieInf",
        susie_ash_inf = "susieAsh",
        mrash = "mrash",
        enet = "enet",
        ridge = "ridge",
        lasso = "lasso",
        bayesR = "bayesR",
        bayesL = "bayesL",
        bayesA = "bayesA",
        bayesB = "bayesB",
        bayesC = "bayesC",
        bayesN = "bayesN",
        bLasso = "bLasso",
        dprVb = "dprVb",
        dprGibbs = "dprGibbs",
        dprAdaptiveGibbs = "dprAdaptiveGibbs",
        scad = "scad",
        mcp = "mcp",
        l0learn = "l0learn",
        mvsusie = "mvsusie",
        mrmash = "mrmash",
        prsCs = "prsCs",
        fsusie = "fsusie"
    )
    unname(map_chr(
        snake,
        .twasCanonicalMethod,
        snakeToCanonical = snakeToCanonical
    ))
}

# Enforce input-class / method compatibility against the TWAS
# capability table. Routes the input class through individual /
# sumstat branches; the twasWeightsPipeline has no GwasSumStats input
# path so that branch is omitted. Emits a single error listing every
# offending token.
# @noRd
.twasCheckMethodCapabilities <- function(tokens, inputKind) {
    if (length(tokens) == 0L) {
        return(invisible(NULL))
    }
    caps <- .twasMethodCapabilities
    # Fine-mapping tokens are governed by .twasCheckFineMappingMethods (and
    # delegate input-class compat to .fmCheckMethodCapabilities); skip them here
    # so they aren't reported as "unknown".
    tokens <- setdiff(tokens, intersect(tokens, .twasFineMappingTokens()))
    if (length(tokens) == 0L) {
        return(invisible(NULL))
    }
    .twasCheckUnknownTokens(tokens, caps)
    violations <- compact(map(
        tokens,
        .twasTokenViolation,
        caps = caps,
        inputKind = inputKind
    ))
    if (length(violations) > 0L) {
        .twasStopCapability(violations, inputKind)
    }
    invisible(NULL)
}

# Error out when any token is absent from the capability table.
# @noRd
.twasCheckUnknownTokens <- function(tokens, caps) {
    unknown <- setdiff(tokens, names(caps))
    if (length(unknown) == 0L) {
        return(invisible(NULL))
    }
    unknownStr <- str_flatten(unknown, ", ")
    knownStr <- str_flatten(c(names(caps), .twasFineMappingTokens()), ", ")
    msg <- glue(
        "twasWeightsPipeline: unknown method token(s): {unknownStr}. ",
        "Known tokens: {knownStr}."
    )
    abort(msg)
}

# Return a list(token, reason) record when `tk` is incompatible with the input
# class, else NULL.
# @noRd
.twasCapabilityViolation <- function(info, tk, inputKind) {
    individualKinds <- c("QtlDataset", "MultiStudyQtlDataset")
    if (is_in(inputKind, individualKinds) && is.null(info$individualImpl)) {
        return(list(
            token = tk,
            reason = "is sumstat-only (use a QtlSumStats input)"
        ))
    }
    if (inputKind == "QtlSumStats" && is.null(info$sumstatImpl)) {
        return(list(
            token = tk,
            reason = "is individual-only (use a QtlDataset input)"
        ))
    }
    # twasWeightsPipeline does not support GwasSumStats input.
    NULL
}

# Raise the aggregated incompatible-method error from a list of violations.
# @noRd
.twasStopCapability <- function(violations, inputKind) {
    bad <- map_chr(violations, "token")
    reason <- map_chr(violations, "reason")
    badStr <- str_flatten(unique(bad), ", ")
    detailStr <- str_flatten(glue("{bad} {reason}"), "; ")
    msg <- glue(
        "twasWeightsPipeline: the following method(s) are not available ",
        "for input class '{inputKind}': {badStr}. {detailStr}."
    )
    abort(msg)
}

# Adapter registry mapping each fine-mapping method (whose existence is
# governed by .fineMappingMethodCapabilities) to its TWAS-weight extractor
# wrapper. The wrapper names follow the *Weights / *RssWeights convention,
# and the *Fit argument receives the pre-fitted fine-mapping object.
# fSuSiE is multivariate (it collapses a functional fit to a variants x
# features weight matrix via fsusieWeights) and has no RSS counterpart.
# @noRd
.twasFineMappingMethodAdapters <- list(
    susie = list(
        weightFn = "susieWeights",
        rssWeightFn = "susieRssWeights",
        fitArg = "susieFit",
        rssFitArg = "susieRssFit",
        methodKey = "susie_weights"
    ),
    susieInf = list(
        weightFn = "susieInfWeights",
        rssWeightFn = "susieInfRssWeights",
        fitArg = "susieInfFit",
        rssFitArg = "susieInfRssFit",
        methodKey = "susie_inf_weights"
    ),
    susieAsh = list(
        weightFn = "susieAshWeights",
        rssWeightFn = "susieAshRssWeights",
        fitArg = "susieAshFit",
        rssFitArg = "susieAshRssFit",
        methodKey = "susie_ash_weights"
    ),
    mvsusie = list(
        weightFn = "mvsusieWeights",
        rssWeightFn = "mvsusieRssWeights",
        fitArg = "mvsusieFit",
        rssFitArg = "mvsusieRssFit",
        methodKey = "mvsusie_weights"
    ),
    fsusie = list(
        weightFn = "fsusieWeights",
        rssWeightFn = NULL,
        fitArg = "fsusieFit",
        rssFitArg = NULL,
        methodKey = "fsusie_weights"
    )
)

# Canonical list of fine-mapping tokens recognised by twasWeightsPipeline:
# fineMappingPipeline's registry, which now contains only fine-mapping methods
# (mr.mash is a TWAS method, kept out of that registry).
# @noRd
.twasFineMappingTokens <- function() {
    names(.fineMappingMethodCapabilities)
}

# Canonical fine-mapping tokens actually present as methods in a
# fine-mapping result, matched tolerantly across canonical / camelCase /
# snake_case spellings (mirrors the candidate logic in .twasFineMappingFits).
# @noRd
.twasFineMappingMethodsPresent <- function(fineMappingResult) {
    if (is.null(fineMappingResult)) {
        return(character(0))
    }
    methods <- str_to_lower(as.character(fineMappingResult$method))
    keep(.twasFineMappingTokens(), .twasTokenPresentIn, methods = methods)
}

# The spellings one canonical token may appear under: itself, its camelCase
# form, and its snake_case form -- all lowercased for comparison.
# @noRd
.twasSpellingCandidates <- function(canonical) {
    str_to_lower(c(
        canonical,
        str_c(
            str_to_lower(str_sub(canonical, 1L, 1L)),
            str_sub(canonical, 2L)
        ),
        str_replace_all(canonical, "([A-Z])", "_\\1")
    ))
}

# @noRd
.twasTokenPresentIn <- function(canonical, methods) {
    any(is_in(methods, .twasSpellingCandidates(canonical)))
}

# Reject fine-mapping methods (susie / susieInf / susieAsh / mvsusie /
# fsusie) when no fine-mapping result is supplied. twasWeightsPipeline is
# not allowed to re-fit fine-mapping models from scratch; users must run
# fineMappingPipeline() first and pass the result via `fineMappingResult`.
# Input-class compatibility (e.g. fsusie has no QtlSumStats path) is
# delegated to .fmCheckMethodCapabilities so the rule set stays in lock-
# step with fineMappingPipeline. Methods with no TWAS-weight extractor
# (fsusie) are rejected with a method-specific message.
# @noRd
.twasCheckFineMappingMethods <- function(
    tokens,
    fineMappingResult,
    inputKind,
    cvFolds = 0,
    fitFullData = TRUE
) {
    fmTokens <- intersect(tokens, .twasFineMappingTokens())
    if (length(fmTokens) == 0L) {
        return(invisible(NULL))
    }
    # Defer input-class compatibility to fineMappingPipeline. e.g. this rejects
    # fsusie on QtlSumStats (fsusie has no RSS impl).
    .fmCheckMethodCapabilities(fmTokens, inputKind)
    .twasCheckFmAdapters(fmTokens)
    .twasRequireFmResult(fineMappingResult, fmTokens)
    .twasCheckFmPresent(fmTokens, fineMappingResult)
    .twasCheckFmCvPresent(fmTokens, fineMappingResult, cvFolds)
    if (isTRUE(fitFullData)) {
        .twasCheckFmFullFits(fmTokens, fineMappingResult)
    }
    invisible(NULL)
}

# Whether row `i` of a fine-mapping result carries a cross-validation result.
# @noRd
.fmrRowHasCv <- function(i, fineMappingResult) {
    !is.null(cvResult(.fmrRowParts(fineMappingResult, i)))
}

# Whether a fine-mapping result carries any cross-validation result at all.
# @noRd
.fmrAnyCvResult <- function(fineMappingResult) {
    if (!is(fineMappingResult, "FineMappingResultBase")) {
        return(FALSE)
    }
    n <- nrow(fineMappingResult)
    if (is.null(n) || n == 0L) {
        return(FALSE)
    }
    any(map_lgl(
        seq_len(n),
        .fmrRowHasCv,
        fineMappingResult = fineMappingResult
    ))
}

# Cross-validating a fine-mapping method needs that method's per-fold fits,
# and this pipeline never fine-maps for itself. Only fineMappingPipeline() run
# with cvFolds > 1 produces them, so require them up front rather than
# discovering the gap once the folds are already being scored.
# @noRd
.twasCheckFmCvPresent <- function(fmTokens, fineMappingResult, cvFolds) {
    if (is.null(cvFolds) || cvFolds <= 1L) {
        return(invisible(NULL))
    }
    if (.fmrAnyCvResult(fineMappingResult)) {
        return(invisible(NULL))
    }
    fmStr <- str_flatten(fmTokens, ", ")
    msg <- glue(
        "twasWeightsPipeline: cross-validating method(s) {fmStr} needs each ",
        "fold's own fine-mapping fit, but the supplied fineMappingResult ",
        "carries no cross-validation. Run fineMappingPipeline() with ",
        "cvFolds > 1; its fold partition is then reused for every other ",
        "weight method so all of them are scored on the same folds."
    )
    abort(msg)
}

# Full-data weights for a fine-mapping-derived method are read off the
# fine-mapping fit, so a row from a cross-validation-only fine-mapping run
# (CrossValidationParam(fitFullData = FALSE)) cannot supply them. Say so up
# front rather than failing inside the extraction.
# @noRd
.twasCheckFmFullFits <- function(fmTokens, fineMappingResult) {
    rows <- seq_len(nrow(fineMappingResult))
    method <- as.character(fineMappingResult$method)
    noFit <- rows[
        is_in(method, fmTokens) &
            !map_lgl(rows, .fmrRowHasFit, fineMappingResult = fineMappingResult)
    ]
    if (length(noFit) == 0L) {
        return(invisible(NULL))
    }
    msg <- glue(
        "twasWeightsPipeline: {length(noFit)} fineMappingResult row(s) for ",
        "{str_flatten(unique(method[noFit]), ', ')} carry no full-data fit -- ",
        "they come from a cross-validation-only fine-mapping run. Rerun ",
        "fineMappingPipeline() with CrossValidationParam(fitFullData = TRUE), ",
        "or run twasWeightsPipeline() cross-validation-only as well."
    )
    abort(msg)
}

# Reject fine-mapping methods that have no TWAS-weight extractor (e.g. fsusie).
# @noRd
.twasCheckFmAdapters <- function(fmTokens) {
    noAdapter <- setdiff(fmTokens, names(.twasFineMappingMethodAdapters))
    if (length(noAdapter) == 0L) {
        return(invisible(NULL))
    }
    noAdapterStr <- str_flatten(noAdapter, ", ")
    msg <- glue(
        "twasWeightsPipeline: method(s) {noAdapterStr} have no ",
        "TWAS-weight extractor. For multi-trait fine-mapping use mvsusie ",
        "via fineMappingResult."
    )
    abort(msg)
}

# A supplied fineMappingResult is mandatory (fine-mapping is never re-fit) and
# must be a FineMappingResultBase.
# @noRd
#' @importFrom checkmate assertClass
.twasRequireFmResult <- function(fineMappingResult, fmTokens) {
    if (is.null(fineMappingResult)) {
        fmStr <- str_flatten(unique(fmTokens), ", ")
        msg <- glue(
            "twasWeightsPipeline: method(s) {fmStr} are fine-mapping ",
            "methods and may not be re-fit by twasWeightsPipeline. Run ",
            "fineMappingPipeline() first and pass the result via ",
            "`fineMappingResult = <QtlFineMappingResult>`."
        )
        abort(msg)
    }
    assertClass(fineMappingResult, "FineMappingResultBase")
    invisible(NULL)
}

# Every requested fine-mapping method must actually be present in the supplied
# fineMappingResult.
# @noRd
.twasCheckFmPresent <- function(fmTokens, fineMappingResult) {
    missingMethods <- setdiff(
        fmTokens,
        .twasFineMappingMethodsPresent(fineMappingResult)
    )
    if (length(missingMethods) == 0L) {
        return(invisible(NULL))
    }
    missingStr <- str_flatten(unique(missingMethods), ", ")
    msg <- glue(
        "twasWeightsPipeline: method(s) {missingStr} were requested but ",
        "the supplied fineMappingResult contains no such fine-mapping ",
        "fit. Run fineMappingPipeline() with method(s) {missingStr} first ",
        "and pass the result via `fineMappingResult = <QtlFineMappingResult>`."
    )
    abort(msg)
}

# Look up the multivariate flag for a token. Checks the TWAS-regression
# capability table first; if absent, falls back to the fine-mapping
# capability table (the source of truth for susie / mvsusie / fsusie /
# etc.). Returns FALSE for unknown tokens.
# @noRd
.twasIsMultivariateToken <- function(token) {
    info <- .twasMethodCapabilities[[token]]
    if (!is.null(info)) {
        return(isTRUE(info$multivariate))
    }
    fmInfo <- .fineMappingMethodCapabilities[[token]]
    if (!is.null(fmInfo)) {
        return(isTRUE(fmInfo$multivariate))
    }
    FALSE
}

# The multivariate rule across a MultiStudyQtlDataset's components.
#
# A multivariate method -- mvsusie, mr.mash, fsusie, and the PCA alternative
# to fsusie -- fits one model across a tuple's traits and contexts, so it
# needs a multivariate Y from EVERY component it will run on, not just from
# one. A collection mixing multivariate and univariate studies cannot answer
# that, so the request is refused here rather than failing partway through
# the per-component recursion with only some studies fitted.
#
# .twasCheckMultivariateY asks the same question of a single QtlDataset; this
# asks it of each component and names the ones that cannot comply.
# @noRd
.twasCheckMultivariateComponents <- function(tokens, data) {
    mv <- keep(tokens, .twasIsMultivariateToken)
    if (length(mv) == 0L) {
        return(invisible(NULL))
    }
    components <- as.list(data)
    if (length(components) == 0L) {
        return(invisible(NULL))
    }
    bad <- keep(
        names(components) %||% seq_along(components),
        .twasComponentIsUnivariate,
        components = components
    )
    if (length(bad) == 0L) {
        return(invisible(NULL))
    }
    abort(glue(
        "twasWeightsPipeline: method(s) {str_flatten(mv, ', ')} fit one ",
        "model across a tuple's traits and contexts, so every study in a ",
        "MultiStudyQtlDataset must be multivariate. These are not: ",
        "{str_flatten(as.character(bad), ', ')}."
    ))
}

# TRUE when one component offers neither multiple traits nor multiple
# contexts, which is what a multivariate fit needs.
# @noRd
.twasComponentIsUnivariate <- function(key, components) {
    d <- components[[key]]
    nTraits <- length(try_fetch(traitNames(d), error = function(cnd) {
        character(0)
    }))
    nCtx <- length(try_fetch(contexts(d), error = function(cnd) character(0)))
    nTraits < 2L && nCtx < 2L
}

# Enforce the multi-trait / multi-context rule for mvsusie / mr.mash
# methods (same family as the fine-mapping mvSuSiE rule in the design
# doc). Multivariate methods need at least 2 traits *or* at least 2
# contexts in the Y matrix passed to learnTwasWeights.
# @noRd
.twasCheckMultivariateY <- function(tokens, nTraits, nContexts) {
    multivariateTokens <- tokens[map_lgl(tokens, .twasIsMultivariateToken)]
    if (length(multivariateTokens) == 0L) {
        return(invisible(NULL))
    }
    if (nTraits < 2L && nContexts < 2L) {
        mvStr <- str_flatten(multivariateTokens, ", ")
        msg <- glue(
            "twasWeightsPipeline: method(s) {mvStr} require multi-trait or ",
            "multi-context input (got {nTraits} trait(s) x ",
            "{nContexts} context(s))."
        )
        abort(msg)
    }
}

# Reject SumStats inputs that have not been QC'd via summaryStatsQc(), which
# is mandatory before twasWeightsPipeline for SumStats input: it also drops the
# variants absent from the ldReference, so once this passes every entry's SNP
# set is a subset of the ldReference panel.
# @noRd
.twasAssertQcd <- function(sumstats) {
    if (length(qcOptions(sumstats)) == 0L) {
        cls <- class(sumstats)[[1L]]
        msg <- glue(
            "twasWeightsPipeline: the supplied {cls} has no QC record ",
            "(qcOptions() is empty). Call summaryStatsQc() first and pass the ",
            "QC-applied result."
        )
        abort(msg)
    }
}

# Optional resume-cache lookup for twasWeightsPipeline. Returns the
# matching TwasWeightsRow from `twasWeights` for the tuple (study,
# context, trait, method), or NULL when there is no hit. Returns NULL
# silently when twasWeights is NULL or not a TwasWeights collection.
# Mirrors .fmCacheLookup (R/fineMappingPipeline.R).
# @noRd
.twasCacheLookup <- function(twasWeights, studyName, context, trait, method) {
    if (is.null(twasWeights)) {
        return(NULL)
    }
    if (!is(twasWeights, "TwasWeights")) {
        return(NULL)
    }
    idx <- .matchTupleRows(
        twasWeights,
        list(
            study = studyName,
            context = context,
            trait = trait,
            method = method
        )
    )
    if (length(idx) == 0L) {
        return(NULL)
    }
    .twrRowParts(twasWeights, idx[[1L]])
}

# Convert a fine-mapping result (single-method susie/susie_inf row matched
# to the requested study/context/trait) into a `fittedModels` list
# suitable for `learnTwasWeights`. Pulls the trimmedFit from the matching
# entry. Returns a (possibly empty) list.
# @noRd
#' @importFrom checkmate assertClass
.twasFineMappingFits <- function(fineMappingResult, studyName, context, trait) {
    if (is.null(fineMappingResult)) {
        return(list())
    }
    assertClass(fineMappingResult, "FineMappingResultBase")
    tokens <- c("susie", "susieInf", "susieAsh", "mvsusie", "fsusie")
    found <- compact(set_names(
        map(
            tokens,
            .twasFitForToken,
            fineMappingResult = fineMappingResult,
            studyName = studyName,
            context = context,
            trait = trait
        ),
        tokens
    ))
    # compact() on an all-NULL named list leaves a zero-length names
    # attribute, which is not the bare list() the contract promises.
    if (length(found) == 0L) list() else found
}

# The fit this result holds for one token on one tuple, or NULL when it has
# none. The first matching row wins, as the keyed assignment it replaced did.
# @noRd
.twasFitForToken <- function(
    canonical,
    fineMappingResult,
    studyName,
    context,
    trait
) {
    idx <- which(
        is_in(
            str_to_lower(as.character(fineMappingResult$method)),
            .twasSpellingCandidates(canonical)
        ) &
            as.character(fineMappingResult$study) == studyName &
            as.character(fineMappingResult$context) == context &
            as.character(fineMappingResult$trait) == trait
    )
    if (length(idx) == 0L) {
        return(NULL)
    }
    susieFit(.fmrRowParts(fineMappingResult, idx[[1L]]))
}

# Locate a fine-mapping fit for one (study, context, trait, token) tuple.
# Used by the QtlSumStats sumstat dispatcher to pass the precomputed fit
# into susieRssWeights / susieInfRssWeights / susieAshRssWeights /
# mvsusieRssWeights via their respective *Fit arguments.
# @noRd
.twasFineMappingFitFor <- function(
    fineMappingResult,
    studyName,
    context,
    trait,
    token
) {
    if (is.null(fineMappingResult)) {
        return(NULL)
    }
    fits <- .twasFineMappingFits(
        fineMappingResult,
        studyName = studyName,
        context = context,
        trait = trait
    )
    fits[[token]]
}

# Concatenate per-row lists, empty-safe.
# @noRd
.twasCvConcat <- function(pieces) {
    if (length(pieces) == 0L) {
        return(list())
    }
    list_c(pieces)
}

# Row `i`'s cross-validation result, or NULL when it has none. Multi-region
# entries store cvResult as a named per-region list, so the first region that
# carries a partition stands for the row.
# @noRd
.twasRowCvResult <- function(i, fineMappingResult) {
    cv <- cvResult(.fmrRowParts(fineMappingResult, i))
    if (is.null(cv)) {
        return(NULL)
    }
    if (!is.null(cv$samplePartition)) {
        return(cv)
    }
    hit <- keep(cv, .twasCvHasPartition)
    if (length(hit) == 0L) {
        return(NULL)
    }
    hit[[1L]]
}

# Collect the cross-validation payload that fineMappingPipeline stored on the
# fine-mapping result for one (study, context, trait) tuple. It records one
# cvResult per (study, context, trait, method) entry (samplePartition +
# per-fold predictions/metrics, keyed by the TWAS snake method name); this
# merges them across the fine-mapping methods of the tuple into a single
# twasWeightsCv()-shaped list so twasWeightsPipeline can reuse the partition and
# feed those out-of-fold predictions into the SR-TWAS ensemble without re-
# fitting the fine-mapping models. A multi-region entry stores cvResult as a
# per-region list; the first region carrying CV is used. Returns NULL when no
# fine-mapping entry for the tuple recorded CV.
# @noRd
.twasCvResultFor <- function(fineMappingResult, studyName, context, trait) {
    if (is.null(fineMappingResult)) {
        return(NULL)
    }
    if (!is(fineMappingResult, "FineMappingResultBase")) {
        return(NULL)
    }
    idx <- which(
        as.character(fineMappingResult$study) == studyName &
            as.character(fineMappingResult$context) == context &
            as.character(fineMappingResult$trait) == trait
    )
    if (length(idx) == 0L) {
        return(NULL)
    }
    cvs <- compact(map(
        idx,
        .twasRowCvResult,
        fineMappingResult = fineMappingResult
    ))
    prediction <- .twasCvConcat(map(cvs, "prediction"))
    if (length(prediction) == 0L) {
        return(NULL)
    }
    list(
        # The first row that carries one defines the partition, as the
        # "only set it if still NULL" assignment did.
        samplePartition = cvs[[1L]]$samplePartition,
        prediction = prediction,
        performance = .twasCvConcat(map(cvs, "performance"))
    )
}

#' TWAS Weights Pipeline
#'
#' S4-dispatched per-region pipeline for learning TWAS prediction weights.
#' Accepts:
#' \itemize{
#'   \item a \code{\link{QtlDataset}} for individual-level cohort fits;
#'   \item a \code{\link{QtlSumStats}} for per-trait RSS fits;
#'   \item a \code{\link{MultiStudyQtlDataset}}, whose individual-level
#'         studies and optional summary-statistic studies are each fitted,
#'         then combined.
#' }
#'
#' Method-restriction rules (enforced):
#' \itemize{
#'   \item \code{mr.mash}, \code{mvsusie} follow the multi-trait /
#'         multi-context rules of the fine-mapping \code{mvsusie} family
#'         (require at least two traits OR at least two contexts).
#'   \item \code{prsCs} runs only on summary statistics, so it is rejected
#'         on \code{QtlDataset} input.
#'   \item \code{dprVb}, \code{dprAdaptiveGibbs} and the Bayesian alphabet
#'         (\code{bayesA}, \code{bayesB}, \code{bayesC}, \code{bayesL},
#'         \code{bayesN}, \code{bayesR}, \code{bLasso}) run only on
#'         individual-level data, so they are rejected on \code{QtlSumStats}
#'         input.
#'   \item \code{lasso}, \code{enet} and \code{ridge} are glmnet on
#'         individual-level data and the penalized RSS solver on summary
#'         statistics. \code{enet} mixes the two penalties at
#'         \code{alpha = 0.5} unless the method's options set it, on either
#'         input: \code{methods = list(enet = list(alpha = 0.3))}.
#' }
#'
#' A \code{QtlSumStats} input must have been QC'd via
#' \code{\link{summaryStatsQc}} first; otherwise an error is raised pointing
#' at that function.
#'
#' Each row of the returned \code{\link{TwasWeights}} collection records the
#' LD reference its weights were learned against: none for individual-level
#' fits, the reference of the input row(s) for RSS-derived fits.
#'
#' Optionally a \code{QtlFineMappingResult} may be supplied as a source of
#' pre-fit SuSiE / SuSiE-inf / SuSiE-ash objects; their stored fits are passed
#' through to \code{learnTwasWeights} / the RSS sub-pipelines via the
#' \code{fittedModels} slot, avoiding a re-fit.
#'
#' When the supplied result was produced with cross-validation (a
#' \code{\link{CrossValidationParam}} with \code{folds >= 2}), each matching
#' \code{(study, context, trait)} entry's \code{cvResult} is reused: its fold
#' partition becomes the CV partition (unless \code{samplePartition} is given
#' explicitly) and its per-fold out-of-fold predictions/metrics are fed
#' directly into the SR-TWAS ensemble in place of re-fitting those
#' fine-mapping methods here.
#' Non-fine-mapping methods (lasso, enet, ...) are still cross-validated on the
#' same shared partition.
#'
#' @param data A \code{QtlDataset}, \code{MultiStudyQtlDataset}, or
#'   \code{QtlSumStats}. The \code{MultiStudyQtlDataset} method iterates the
#'   embedded individual-level \code{QtlDataset} entries and the optional
#'   embedded \code{QtlSumStats}, then rbinds the results.
#' @param methods A character vector of short method names, a preset string
#'   (\code{"default"} or \code{"fastDefault"}), or a named list of
#'   \code{<method>_weights = args} entries. For QtlSumStats input the
#'   default is \code{c("lasso", "prsCs", "dprGibbs")}, which there run
#'   lassosum, PRS-CS and SDPR; the fine-mapping methods are not in it, and
#'   are requested by name together with a \code{fineMappingResult}.
#' @param contexts Optional context selection; \code{NULL} (the default)
#'   uses every context. A character vector applies to every study: a context
#'   no study has is an error, and one only some studies lack is skipped for
#'   those, with a warning. A list keyed by study name picks each study's
#'   contexts separately (\code{list(A = "liver")}); studies it does not name
#'   keep all of theirs.
#' @param traitId Optional trait selection; \code{NULL} (the default) uses
#'   every trait, or those overlapping \code{region} when one is given. A
#'   character vector applies to every study. A list keyed by study name
#'   picks each study's traits separately (\code{list(A = "ENSG1")}), and a
#'   study's entry may itself be a list keyed by context
#'   (\code{list(A = list(brain = "ENSG1"))}). Studies and contexts the list
#'   does not name keep all of their traits. A joint fit spanning a study's
#'   contexts uses every trait selected in any of them, minus the (context,
#'   trait) pairs a per-context entry leaves out. In a MultiStudyQtlDataset
#'   a study left with no trait selected sits out; the run is refused only
#'   when no study or summary-statistics row has one.
#' @param region Optional variant window for QtlDataset trait selection: a
#'   \code{GRanges}. Mutually exclusive with \code{traitId}.
#' @param cisWindow For QtlDataset: cis-window (bp) around each trait's genomic
#'   position when extracting variants. Required when \code{traitId} is
#'   supplied. Mutually exclusive with \code{region}. For QtlSumStats:
#'   optional; each selected tuple is analysed over its \code{traitPos}
#'   \eqn{\pm} \code{cisWindow}, within its \code{region}, which is then
#'   the window its result records. Every selected tuple needs a
#'   \code{traitPos}. A MultiStudyQtlDataset applies it to each of its
#'   studies, individual-level or summary-statistics, as above.
#' @section Panel filters on the RSS path: On \code{QtlSumStats} input there
#'   is no genotype matrix to filter, so \code{panelFilterParam}'s
#'   \code{mafCutoff} / \code{macCutoff} / \code{imissCutoff} are measured
#'   against the \strong{LD reference panel} instead: a variant whose panel
#'   genotypes fall below the cutoffs is dropped before the z-scores and LD
#'   matrix are built. The thresholds mean the same thing as
#'   \code{genotypeFilterParam}'s on the \code{QtlDataset} path (MAC is
#'   converted to a MAF equivalent and the stricter of the two applies), so
#'   one number carries across input types. The defaults filter nothing.
#'
#'   This discards \emph{observed} variants, unlike
#'   \code{summaryStatsQc(raissParam = ...)}, which only bounds which variants
#'   RAISS will impute.
#'
#' @param genotypeFilterParam For QtlDataset: which variants and samples the
#'   run keeps, built with \code{\link{GenotypeFilterParam}}: its
#'   \code{keepSamples} narrows the dataset once, and its other settings
#'   apply to each extracted block. A field left unset filters nothing.
#'   These are applied identically here and in
#'   \code{\link{fineMappingPipeline}} --- there is deliberately no
#'   TWAS-specific variant filter.
#' @param panelFilterParam For QtlSumStats: LD-reference-panel filters,
#'   built with
#'   \code{\link{PanelFilterParam}}. See \emph{Panel filters on the RSS path}
#'   above.
#' @param jointRegions For QtlDataset with a multi-range \code{region}:
#'   \code{FALSE} (default) learns weights for each range independently and
#'   concatenates them into one entry per (study, context, trait, method); the
#'   per-region fits are kept as a named list and per-region CV is recorded as a
#'   flat \code{cvResult} data frame (one row per region). \code{TRUE}
#'   concatenates the ranges' genotypes into one joint fit. Ignored for a
#'   single-range / cis request.
#' @param jointSpecification Optional joint-fit specification saying which
#'   axes \code{mrmash} fits across. When NULL (the default) it runs its
#'   implicit multi-trait / multi-context fit. Otherwise give a character
#'   vector of axes (one spec) or a list of specs. A spec is a character
#'   vector drawn from \code{"study"}, \code{"context"} and \code{"trait"},
#'   or \code{list(axes = <axes>, scope = <named list>)}, where \code{scope}
#'   (keys \code{study}, \code{context}, \code{trait}) limits the fit to the
#'   listed values, each of which must exist in \code{data}. A list-form
#'   spec always goes inside the outer list, even on its own. Other methods
#'   run per (study, context, trait) as usual. The \code{"study"} axis needs
#'   summary statistics: it is rejected on a \code{QtlDataset}, and on a
#'   \code{MultiStudyQtlDataset} only the summary-statistics studies take
#'   part.
#' @param fineMappingResult Optional \code{QtlFineMappingResult}. When supplied,
#'   its SuSiE / SuSiE-inf / SuSiE-ash trimmed fits for the matching (study,
#'   context, trait) tuples are injected into \code{learnTwasWeights} via
#'   \code{fittedModels} so SuSiE-family weight methods reuse the prior fit
#'   instead of refitting. Its fits and cross-validation are reused whole,
#'   so a row narrowed since its fit warns (condition class
#'   \code{pecotmrFitCoverage}).
#' @param twasWeights Optional \code{\link{TwasWeights}} resume cache. For each
#'   requested \code{(study, context, trait, method)} tuple already present in
#'   this collection, the cached \code{TwasWeightsRow} is copied through and
#'   the corresponding weight fit is skipped. Only the un-cached method subset
#'   is fit; the cached and fresh entries are concatenated in the returned
#'   collection. Per-tuple matching mirrors the \code{fineMappingResult} cache
#'   in \code{\link{fineMappingPipeline}}. Multivariate dispatch
#'   (\code{mvsusie}, \code{mr.mash}) is unaffected because those methods
#'   produce one fit jointly across multiple \code{(context, trait)} columns.
#' @param mashPrior Optional \code{\link{MashPrior}} (or coercible) supplying
#'   data-driven prior matrices for mr.mash weight methods; \code{NULL} to skip.
#' @param usePCA Logical (length 1). Individual-level data only: a
#'   \code{QtlDataset}, or the \code{QtlDataset} components of a
#'   \code{MultiStudyQtlDataset}. When
#'   \code{TRUE}, each multi-trait context is additionally PCA-reduced and
#'   weights are fitted for each top principal component \emph{as a trait}
#'   (\code{topPC1}, \code{topPC2}, ...), alongside the per-trait rows. The
#'   same reduction \code{\link{fineMappingPipeline}} applies, sharing its
#'   scores helper and producing groups the ordinary CV / ensemble /
#'   retention machinery handles unchanged. A single-trait context, or one
#'   with no usable component, contributes nothing.
#' @param nPCs Integer (length 1). Caps the number of top principal
#'   components fitted per context when \code{usePCA = TRUE} (default
#'   \code{10}). The effective count is \code{min(nPCs, usable traits)}.
#' @param fitRetention How much of each fit is kept on the entry:
#'   \code{"slim"} (default), \code{"none"} to keep no fit at all, or
#'   \code{"full"}. The same setting \code{\link{fineMappingPipeline}}
#'   takes, which has no \code{"none"} --- its \code{susieFit} slot is part
#'   of the returned object's contract.
#'
#'   \code{"full"} is honoured by the mr.mash engines and, on summary
#'   statistics, by \code{lasso}, \code{enet}, \code{ridge}, \code{scad},
#'   \code{mcp} and \code{l0learn}, whose full fit also keeps every candidate's
#'   coefficients; every other method keeps its fit whole or not at all.
#'   Requesting \code{"full"} for a run where no method can honour it is an
#'   error, and a mixed run warns rather than dropping it in silence.
#' @param phenotypeFilterParam The run's phenotype QC, built with
#'   \code{\link{PhenotypeFilterParam}}: how missing phenotype values are
#'   handled (unset drops the affected samples; \code{"impute"} fills each
#'   trait's mean; a fit cannot take \code{"keep"}) and whether outlier
#'   samples are dropped. Every fit of the run, joint ones included, reads
#'   its phenotypes under it.
#' @param crossValidationParam Cross-validation settings, built with
#'   \code{\link{CrossValidationParam}}: \code{folds} (default \code{0}, no
#'   CV), \code{threads}, \code{samplePartition}, \code{maxVariants} (cap
#'   on the CV design matrix), \code{weightMethods} (which methods to
#'   cross-validate; unset means every method with non-zero weights) and
#'   \code{fitFullData} (\code{FALSE} keeps only the cross-validation
#'   results; default \code{TRUE}).
#'
#'   \code{folds} defaulted to \code{5} before, while
#'   \code{\link{fineMappingPipeline}}'s defaulted to \code{0}. The two now
#'   agree on \code{0}, so an unspecified setting means one thing.
#'
#'   Cross-validation holds out \strong{samples}, so it is not possible on
#'   \code{QtlSumStats} input; that method \strong{errors} if it is set
#'   rather than returning results that were never cross-validated.
#' @param ensembleParam SR-TWAS ensemble settings, built with
#'   \code{\link{EnsembleParam}}: \code{enabled} (default \code{FALSE}),
#'   \code{r2Threshold}, \code{solver} and \code{alpha}. Stacking reads
#'   out-of-fold predictions, so \code{enabled = TRUE} requires
#'   \code{crossValidationParam} with \code{folds >= 2} and is an error
#'   otherwise --- it previously returned a result quietly missing its
#'   ensemble row.
#' @param estimatePi If TRUE, estimate spike-and-slab sparsity from mr.ash
#'   before BGLR / qgg spike-and-slab methods that consume it.
#' @param residualizationParam Covariate residualization settings, built with
#'   \code{\link{ResidualizationParam}} and forwarded to
#'   \code{\link{residualizedPhenotypes}} /
#'   \code{\link{residualizedGenotypes}}: \code{phenotypeCovariates} and
#'   \code{genotypeCovariates} name which covariates to regress out
#'   (\code{NULL}, the default, uses every available one), and
#'   \code{residualizePhenotype} / \code{residualizeGenotype} turn each side
#'   off.
#'
#'   \code{QtlDataset} / \code{MultiStudyQtlDataset} input only. A
#'   \code{QtlSumStats} run has no covariates to regress out, so the bundle
#'   is \strong{ignored} there rather than refused --- unlike
#'   \code{crossValidationParam}, which errors. Residualization is on by
#'   default,
#'   so refusing a non-default value would reject the default bundle; CV is
#'   off by default, so a non-default value there is an explicit request for
#'   something the input cannot do.
#' @param weightsDataType Optional data-type label recorded on every
#'   \code{TwasWeightsRow$dataType} (e.g. \code{"expression"}).
#' @section The genotype scale:
#'   Individual-level genotypes reach every method centred and scaled to unit
#'   variance, so all the weights are on one scale and the recorded
#'   \code{weight_standardized} says so. The residualization
#'   (\code{\link{residualizedGenotypes}}) applies it, governed by
#'   \code{ResidualizationParam(scaleResiduals = )}, which defaults to
#'   \code{TRUE}.
#'
#'   It is not cosmetic: measured through pecotmr's own entry points, every
#'   method returns coefficients on the scale of the \code{X} it was given
#'   --- exactly so for glmnet / ncvreg / L0Learn / mr.ash / the Bayesian
#'   alphabet and for \code{susieWeights()} --- while the DPR samplers are
#'   not scale-equivariant at all, so a uniform rescale of \code{X} changes
#'   their fit and not merely its units. The chosen scale is also the one the
#'   LD correlations and z-scores behind a TWAS Z already use. A
#'   summary-statistics fit needs nothing done.
#'
#'   \code{\link{learnTwasWeights}} still takes \code{weightStandardized}:
#'   there the caller supplies \code{X} themselves, so only they can say.
#' @param verbose Verbosity (0 silent, 1 default, 2 includes external package
#'   messages).
#' @param seed Integer or \code{NULL}. When supplied (\code{QtlDataset} path),
#'   seeds both the main-process RNG and the parallel method-fitting /
#'   cross-validation RNG (\code{BiocParallel} \code{RNGseed}) for
#'   reproducibility under multi-threading. The main-process seed is scoped to
#'   the call, so the session RNG is left as it was found. \code{NULL}
#'   (default) does not seed at all, so an outer \code{set.seed()} still
#'   governs the main-process draws.
#' @param ... Reserved for method-specific arguments.
#'
#' @return A \code{\link{TwasWeights}} collection keyed by \code{(study,
#'   context, trait, method)}. Each row's LD reference is none for
#'   individual-level fits and that of the input row(s) for RSS-derived fits.
#'   Its \code{seqinfo()} names the genome build the inputs name; inputs
#'   naming different builds are refused.
#' @examples
#' data(qtlDatasetExample)
#' twasWeightsPipeline(qtlDatasetExample, methods = "lasso", cisWindow = 1e6)
#' @importFrom purrr imap map map_chr map_int map_lgl compact list_flatten
#'   set_names
#' @export
setGeneric("twasWeightsPipeline", function(data, ...) {
    .twasAssertInputClass(data)
    standardGeneric("twasWeightsPipeline")
})

# @noRd
.twasAssertInputClass <- function(data) {
    ok <- c("QtlDataset", "MultiStudyQtlDataset", "QtlSumStats")
    if (any(map_lgl(ok, is, object = data))) {
        return(invisible(NULL))
    }
    abort(glue(
        "twasWeightsPipeline does not accept inputs of class ",
        "'{class(data)[[1L]]}'. Pass a QtlDataset, MultiStudyQtlDataset, ",
        "or QtlSumStats. (GwasSumStats inputs are not supported; GWAS-side ",
        "per-LD-block weights are produced inside the new ctwasPipeline / ",
        "qtlEnrichmentPipeline.)"
    ))
}

# Run joint specs through the engine over every region block of `ctx` and
# merge the per-region results by (study, context, trait, method). Shared by
# the per-tuple multivariate grid (a synthesized context x trait spec) and the
# explicit jointSpecification path, so both fit under the same pipeline
# `marker` and the same engine `args`. `ctx` is .twasJointSetup()'s record.
# @noRd
.twasRunJointSpecGrid <- function(spec, tokens, contexts, traits, marker, ctx) {
    allLabs <- map_chr(ctx$xRegions, .twasRegionLabel)
    allRegions <- map(
        seq_along(ctx$xRegions),
        .twasJointSpecRegion,
        spec = spec,
        tokens = tokens,
        contexts = contexts,
        traits = traits,
        marker = marker,
        ctx = ctx
    )
    keep <- !map_lgl(allRegions, is.null)
    perRegion <- allRegions[keep]
    labs <- allLabs[keep]
    if (length(perRegion) == 0L) {
        return(NULL)
    }
    if (length(perRegion) == 1L) {
        return(perRegion[[1L]])
    }
    .twasMergeResultsByKey(perRegion, labs)
}

# ---- QtlDataset pipeline worker + phase helpers ----------------------------

# The multivariate half of .twasQdsRunGrid: one fit across all traits,
# against the shared grid context.
# @noRd
.twasQdsRunMultivariate <- function(
    grid,
    marker,
    data,
    xRegions,
    norm,
    fineMappingResult,
    dataDrivenPriorMatricesCv,
    cisWindow,
    residualizationParam,
    genotypeFilter,
    verbose
) {
    # A per-context selection goes in whole: the joint scope then takes each
    # study's union of it and drops the (context, trait) pairs it leaves out.
    .twasRunJointSpecGrid(
        list(list(axes = c("context", "trait"), scope = NULL)),
        norm$tokens,
        grid$useCtx,
        grid$traitSelection %||% grid$allTraits,
        marker,
        list(
            xRegions = xRegions,
            data = data,
            norm = norm,
            fineMappingResult = fineMappingResult,
            dataDrivenPriorMatricesCv = dataDrivenPriorMatricesCv,
            cisWindow = cisWindow,
            residualizationParam = residualizationParam,
            genotypeFilter = genotypeFilter,
            verbose = verbose
        )
    )
}

# Run the resolved grid: the multivariate engine when the grid says so, the
# per-tuple univariate engine otherwise. `grid` carries the axes
# (.twasQdsResolveGrid) and `marker` the per-run fit settings.
# @noRd
.twasQdsRunGrid <- function(
    grid,
    marker,
    data,
    xRegions,
    norm,
    fineMappingResult,
    twasWeights,
    dataDrivenPriorMatricesCv,
    cisWindow,
    phenotypeFilter,
    residualizationParam,
    genotypeFilter,
    usePCA,
    nPCs,
    verbose
) {
    if (grid$multivariate) {
        return(.twasQdsRunMultivariate(
            grid,
            marker,
            data = data,
            xRegions = xRegions,
            norm = norm,
            fineMappingResult = fineMappingResult,
            dataDrivenPriorMatricesCv = dataDrivenPriorMatricesCv,
            cisWindow = cisWindow,
            residualizationParam = residualizationParam,
            genotypeFilter = genotypeFilter,
            verbose = verbose
        ))
    }
    .twasQdsUnivariateEngine(
        grid,
        marker = marker,
        data = data,
        xRegions = xRegions,
        norm = norm,
        fineMappingResult = fineMappingResult,
        twasWeights = twasWeights,
        dataDrivenPriorMatricesCv = dataDrivenPriorMatricesCv,
        cisWindow = cisWindow,
        phenotypeFilter = phenotypeFilter,
        residualizationParam = residualizationParam,
        genotypeFilter = genotypeFilter,
        usePCA = usePCA,
        nPCs = nPCs,
        verbose = verbose
    )
}

# The QtlDataset worker's argument checks, grouped so its entry point names
# the validation phase once instead of listing four calls inline.
# @noRd
.twasQdsAssertParams <- function(
    methods,
    fitRetention,
    crossValidationParam,
    ensembleParam
) {
    .twasWarnUnretainableDetail(fitRetention, methods, "QtlDataset")
    .assertMethodParam(
        crossValidationParam,
        "CrossValidationParam"
    )
    .assertMethodParam(ensembleParam, "EnsembleParam")
    .ensembleAssertCv(ensembleParam, crossValidationParam)
    invisible(NULL)
}

# The QtlDataset method of twasWeightsPipeline(). Each stage returns only the
# values it derives; nothing is grafted onto a captured environment. `cfg`
# (.twasJointCfg()) is the one record of the run's settings, and the pipeline
# marker is built from it BEFORE the joint phase, so explicit joint fits run
# under the same marker and engine args as the per-tuple loop. The joint phase
# consumes the mrmash token; `norm` comes back holding only the methods that
# still have to go through the per-tuple loop.
# @noRd
.twasPipelineQtlDataset <- function(
    data,
    methods = "default",
    contexts = NULL,
    traitId = NULL,
    region = NULL,
    cisWindow = NULL,
    genotypeFilterParam = GenotypeFilterParam(),
    jointRegions = FALSE,
    jointSpecification = NULL,
    fineMappingResult = NULL,
    twasWeights = NULL,
    mashPrior = NULL,
    crossValidationParam = CrossValidationParam(),
    ensembleParam = EnsembleParam(),
    estimatePi = TRUE,
    usePCA = FALSE,
    nPCs = 10L,
    fitRetention = c("slim", "none", "full"),
    residualizationParam = ResidualizationParam(),
    weightsDataType = NULL,
    phenotypeFilterParam = PhenotypeFilterParam(),
    verbose = 1,
    seed = NULL
) {
    phenotypeFilter <- .twasResolvePhenotypeFilter(phenotypeFilterParam)
    fitRetention <- arg_match(fitRetention)
    build <- .inputsBuild(
        list(data, region, fineMappingResult, twasWeights),
        "twasWeightsPipeline"
    )
    .fitCoverageCheckedRows(fineMappingResult, "twasWeightsPipeline()")
    contexts <- parseContexts(contexts, data)
    traitId <- parseTraitIds(traitId, data)
    .twasQdsAssertParams(
        methods,
        fitRetention,
        crossValidationParam,
        ensembleParam
    )
    resolved <- .twasQdsResolveInputs(
        data = data,
        region = region,
        cisWindow = cisWindow,
        jointRegions = jointRegions,
        genotypeFilterParam = genotypeFilterParam
    )
    data <- resolved$data
    cfg <- .twasJointCfg(
        genotypeFilter = resolved$genotypeFilter,
        phenotypeFilter = phenotypeFilter,
        contexts = contexts,
        traitId = traitId,
        cisWindow = cisWindow,
        xRegions = resolved$xRegions,
        fineMappingResult = fineMappingResult,
        mashPrior = mashPrior,
        crossValidationParam = crossValidationParam,
        ensembleParam = ensembleParam,
        estimatePi = estimatePi,
        fitRetention = fitRetention,
        residualizationParam = residualizationParam,
        weightsDataType = weightsDataType,
        verbose = verbose,
        seed = seed
    )
    resolvedTokens <- .twasQdsResolveTokens(
        jointSpecification,
        methods,
        fineMappingResult,
        .cvResolve(crossValidationParam),
        data
    )
    setup <- .twasJointSetup(data, resolvedTokens$norm, cfg)
    joint <- .twasQdsJointPhase(
        resolvedTokens$parsedJointSpec,
        setup$marker,
        setup$ctx,
        contexts,
        traitId
    )
    if (joint$done) {
        return(.withGenome(joint$result, build))
    }
    grid <- .twasQdsResolveGrid(
        data,
        contexts,
        traitId,
        region,
        joint$norm$tokens
    )
    tw <- .twasQdsRunGrid(
        grid,
        setup$marker,
        data = data,
        xRegions = resolved$xRegions,
        norm = joint$norm,
        fineMappingResult = fineMappingResult,
        twasWeights = twasWeights,
        dataDrivenPriorMatricesCv = setup$ctx$dataDrivenPriorMatricesCv,
        cisWindow = cisWindow,
        phenotypeFilter = phenotypeFilter,
        residualizationParam = residualizationParam,
        genotypeFilter = resolved$genotypeFilter,
        usePCA = usePCA,
        nPCs = nPCs,
        verbose = verbose
    )
    .withGenome(.twasQdsAssemble(tw, joint$result), build)
}

#' @rdname twasWeightsPipeline
#' @export
setMethod(
    "twasWeightsPipeline",
    "QtlDataset",
    .twasPipelineQtlDataset
)

# Unpack the MashPrior bundle: route the full-data prior into the mr.mash method
# args, the per-fold priors + fold partition into the CV machinery. Returns the
# updated parameter bundle.
# @noRd
.twasQdsUnpackMash <- function(
    mashPrior,
    samplePartition,
    norm
) {
    mp <- .unpackMashPrior(mashPrior, samplePartition)
    if (!is.null(mashPrior) && !is_in("mrmash", norm$tokens)) {
        msg <- glue(
            "`mashPrior` was supplied but 'mrmash' is not among `methods`; ",
            "the data-driven prior is ignored."
        )
        warn(msg)
    }
    withPrior <- if (
        !is.null(mp$fullPrior) &&
            is_in("mrmash_weights", names(norm$methodList))
    ) {
        .twasQdsSetMrmashPrior(norm, mp$fullPrior)
    } else {
        norm
    }
    list(
        samplePartition = mp$samplePartition,
        dataDrivenPriorMatricesCv = mp$dataDrivenPriorMatricesCv,
        norm = withPrior
    )
}

# The normalized argument list with mr.mash's full-data prior attached,
# leaving every other method's arguments untouched.
# @noRd
.twasQdsSetMrmashPrior <- function(norm, fullPrior) {
    updated <- list_assign(
        norm$methodList$mrmash_weights,
        dataDrivenPriorMatrices = fullPrior
    )
    list_assign(
        norm,
        methodList = list_assign(norm$methodList, mrmash_weights = updated)
    )
}

# Reject the region + cisWindow combination, derive the X windows, resolve the
# run's genotype filter and narrow the dataset to its samples. Mirrors
# .fmQdsResolveInputs; variant QC is applied identically to fine-mapping and
# TWAS, so there is no TWAS-specific variant filter.
# @noRd
.twasQdsResolveInputs <- function(
    data,
    region,
    cisWindow,
    jointRegions,
    genotypeFilterParam
) {
    .rejectRegionWithCisWindow(
        region,
        cisWindow,
        "twasWeightsPipeline(QtlDataset)"
    )
    genotypeFilter <- .qtlResolveFilter(genotypeFilterParam)
    list(
        data = .qtlRestrictSamples(data, genotypeFilter$keepSamples),
        xRegions = .makeXRegions(region, jointRegions),
        genotypeFilter = genotypeFilter
    )
}


# Explicit jointSpecification path: run the per-spec axis dispatcher for
# mr.mash. Returns list(done, result, norm): `done` requests an early return
# with `result`; otherwise `norm` is the mrmash-stripped normalization for the
# per-tuple loop below.
# @noRd
.twasQdsJointPhase <- function(
    parsedJointSpec,
    marker,
    ctx,
    contexts,
    traitId
) {
    norm <- ctx$norm
    if (length(parsedJointSpec) == 0L) {
        return(list(done = FALSE, result = NULL, norm = norm))
    }
    jointResult <- .twasDispatchJointSpecsQtlDataset(
        parsedJointSpec,
        norm$tokens,
        marker,
        ctx,
        contexts,
        traitId
    )
    keep <- setdiff(norm$tokens, intersect(norm$tokens, "mrmash"))
    if (length(keep) == 0L) {
        if (is.null(jointResult)) {
            msg <- glue(
                "twasWeightsPipeline(QtlDataset): no joint fits produced. ",
                "Check that the jointSpecification scope intersects the ",
                "available studies / contexts / traits."
            )
            abort(msg)
        }
        return(list(done = TRUE, result = jointResult))
    }
    keepKeys <- which(
        is_in(str_remove(names(norm$methodList), "(_weights|Weights)$"), keep)
    )
    list(
        done = FALSE,
        result = jointResult,
        norm = list_assign(
            norm,
            tokens = keep,
            methodList = norm$methodList[keepKeys]
        )
    )
}

# Parse the joint spec and normalize the methods, then run the capability,
# per-input argument and fine-mapping gates. The mash-prior unpacking and the
# joint phase follow once the pipeline marker exists (.twasJointSetup()), so
# that explicit joint fits run under the same marker as everything else.
# @noRd
.twasQdsResolveTokens <- function(
    jointSpecification,
    methods,
    fineMappingResult,
    cvCfg,
    data
) {
    parsedJointSpec <- parseJointSpecification(jointSpecification, data)
    norm <- .twasNormalizeMethods(methods)
    .twasCheckMethodCapabilities(norm$tokens, "QtlDataset")
    .twasCheckMethodArgsForInput(
        .twasMethodListArgs(norm$methodList),
        "QtlDataset"
    )
    .twasCheckFineMappingMethods(
        norm$tokens,
        fineMappingResult,
        "QtlDataset",
        cvFolds = cvCfg$folds,
        fitFullData = cvCfg$fitFullData
    )
    list(parsedJointSpec = parsedJointSpec, norm = norm)
}

# TWAS's phenotype QC: a weight fit cannot take a missing phenotype, so an
# unset naAction drops the affected samples and "keep" is refused.
# @noRd
.twasResolvePhenotypeFilter <- function(phenotypeFilterParam) {
    .qtlResolvePhenotypeFilter(
        phenotypeFilterParam,
        naDefault = "drop",
        naAllowed = c("drop", "impute"),
        caller = "twasWeightsPipeline"
    )
}

# The individual-level settings a TWAS run fits under, as one record: built
# once by the QtlDataset method and once by the MultiStudyQtlDataset method,
# with the same builder, so a setting added for one reaches the other.
# @noRd
.twasJointCfg <- function(
    genotypeFilter,
    phenotypeFilter,
    contexts,
    traitId,
    cisWindow,
    xRegions,
    fineMappingResult,
    mashPrior,
    crossValidationParam,
    ensembleParam,
    estimatePi,
    fitRetention,
    residualizationParam,
    weightsDataType,
    verbose,
    seed
) {
    list(
        genotypeFilter = genotypeFilter,
        phenotypeFilter = phenotypeFilter,
        contexts = contexts,
        traitId = traitId,
        cisWindow = cisWindow,
        xRegions = xRegions,
        fineMappingResult = fineMappingResult,
        mashPrior = mashPrior,
        crossValidationParam = crossValidationParam,
        ensembleParam = ensembleParam,
        estimatePi = estimatePi,
        fitRetention = fitRetention,
        residualizationParam = residualizationParam,
        weightsDataType = weightsDataType,
        verbose = verbose,
        seed = seed
    )
}

# Everything an individual-level TWAS fit on `data` needs from the record
# `cfg` (.twasJointCfg()): the mash prior unpacked into the methods, the
# pipeline marker, and the engine-args context. Used for the QtlDataset
# method's own run and for each MultiStudyQtlDataset component's joint fits.
# @noRd
.twasJointSetup <- function(data, norm, cfg) {
    cvCfg <- .cvResolve(cfg$crossValidationParam)
    mash <- .twasQdsUnpackMash(cfg$mashPrior, cvCfg$samplePartition, norm)
    marker <- .twasQdsMarker(
        crossValidationParam = list_assign(
            cvCfg,
            samplePartition = mash$samplePartition
        ),
        ensembleParam = .ensembleResolve(cfg$ensembleParam),
        weightsDataType = cfg$weightsDataType,
        standardized = .qtlScaleResiduals(cfg$residualizationParam),
        fitRetention = cfg$fitRetention,
        estimatePi = cfg$estimatePi,
        verbose = cfg$verbose,
        seed = cfg$seed
    )
    list(
        marker = marker,
        ctx = list(
            data = data,
            xRegions = cfg$xRegions,
            norm = mash$norm,
            fineMappingResult = cfg$fineMappingResult,
            dataDrivenPriorMatricesCv = mash$dataDrivenPriorMatricesCv,
            cisWindow = cfg$cisWindow,
            residualizationParam = cfg$residualizationParam,
            genotypeFilter = cfg$genotypeFilter,
            phenotypeFilter = cfg$phenotypeFilter,
            verbose = cfg$verbose
        )
    )
}

# Resolve the (context, trait) grid + multivariate flag and build the joint-
# pipeline marker + shared grid context. Returns the updated parameter bundle.
# @noRd
.twasQdsResolveGrid <- function(
    data,
    contexts,
    traitId,
    region,
    tokens
) {
    studyName <- studyName(data)
    useCtx <- .twasQdsResolveContexts(data, contexts)
    allTraits <- .twasQdsResolveTraits(data, useCtx, traitId, region)
    .twasCheckMultivariateY(tokens, length(allTraits), length(useCtx))
    list(
        study = studyName,
        useCtx = useCtx,
        allTraits = allTraits,
        # Set only when traits are chosen context by context: the grid spans
        # every (context, trait) pair of allTraits, and this drops the pairs
        # the selection leaves out.
        traitSelection = if (.hasNestedTraitSelection(traitId)) traitId,
        multivariate = any(map_lgl(tokens, .twasIsMultivariateToken))
    )
}

# Selected contexts: those parseContexts() kept for this study, or all.
# @noRd
.twasQdsResolveContexts <- function(data, contexts) {
    .selectedContexts(contexts, studyName(data)) %||% contexts(data)
}

# Traits to iterate: traitId when supplied, else per-context region overlap,
# else every trait in every selected context.
# @noRd
.twasQdsResolveTraits <- function(data, useCtx, traitId, region) {
    perCtxTraits <- map(
        useCtx,
        .contextSelectedTraits,
        qd = data,
        traitId = traitId,
        region = region
    )
    allTraits <- unique(list_c(perCtxTraits))
    if (length(allTraits) == 0L) {
        abort("twasWeightsPipeline(QtlDataset): no traits selected.")
    }
    allTraits
}

# Joint-pipeline marker carrying the CV / ensemble config for the engine.
# @noRd
.twasQdsMarker <- function(
    crossValidationParam,
    ensembleParam,
    weightsDataType,
    standardized,
    fitRetention,
    estimatePi,
    verbose,
    seed
) {
    new(
        "TwasJointPipeline",
        config = list(
            # The two groups travel as records, not as nine loose fields:
            # this function's only job is to store them, so unrolling them
            # into its signature would name each one twice for nothing.
            crossValidationParam = crossValidationParam,
            ensembleParam = ensembleParam,
            dataType = weightsDataType,
            fitRetention = fitRetention,
            # DERIVED from what the accessors actually do: genotypes reach
            # the fitter scaled to unit variance unless
            # `ResidualizationParam(scaleResiduals = FALSE)` says otherwise.
            # Resolved by the caller, which holds the bundle.
            standardized = standardized,
            estimatePi = estimatePi,
            verbose = verbose,
            seed = seed,
            ldReference = NULL
        )
    )
}

# Which joint cells run, and over which study/context/trait scope. usePCA
# adds top-PC rows ALONGSIDE the per-trait ones, matching
# fineMappingPipeline's c(univRows, pcaRows).
# @noRd
.twasQdsUnivPlan <- function(
    usePCA,
    studyName,
    useCtx,
    allTraits,
    traitSelection
) {
    univCell <- .lookupJointCell("univariate", "individual")
    list(
        cells = if (isTRUE(usePCA)) {
            list(univCell, .lookupJointCell("topPc", "individual"))
        } else {
            list(univCell)
        },
        scope = list(
            studies = studyName,
            contexts = set_names(list(useCtx), studyName),
            traits = set_names(list(allTraits), studyName),
            traitSelection = traitSelection
        )
    )
}

# Univariate methods ROUTED THROUGH THE ENGINE: one 1-condition group per
# (context, trait), per region -> the SAME per-method fitter (+ ensemble layer
# for >= 2 methods + resume cache) as the joint paths, merged across regions.
# `grid` is the dispatch row naming the study, contexts and traits.
# @noRd
.twasQdsUnivariateEngine <- function(
    grid,
    marker,
    data,
    xRegions,
    norm,
    fineMappingResult,
    twasWeights,
    dataDrivenPriorMatricesCv,
    cisWindow,
    phenotypeFilter,
    residualizationParam,
    genotypeFilter,
    usePCA,
    nPCs,
    verbose
) {
    plan <- .twasQdsUnivPlan(
        usePCA,
        grid$study,
        grid$useCtx,
        grid$allTraits,
        grid$traitSelection
    )
    cells <- plan$cells
    scope <- plan$scope
    labs <- map_chr(xRegions, .twasRegionLabel)
    perRegion <- map(
        seq_along(xRegions),
        .twasQdsUnivRegion,
        cells = cells,
        scope = scope,
        marker = marker,
        data = data,
        xRegions = xRegions,
        norm = norm,
        fineMappingResult = fineMappingResult,
        twasWeights = twasWeights,
        dataDrivenPriorMatricesCv = dataDrivenPriorMatricesCv,
        cisWindow = cisWindow,
        phenotypeFilter = phenotypeFilter,
        residualizationParam = residualizationParam,
        genotypeFilter = genotypeFilter,
        nPCs = nPCs,
        verbose = verbose
    )
    keep <- !map_lgl(perRegion, is.null)
    .twasMergeRegionResults(perRegion[keep], labs[keep])
}

# Per-region engine args for the univariate path.
# @noRd
.twasQdsUnivArgs <- function(
    bi,
    xRegions,
    norm,
    fineMappingResult,
    twasWeights,
    dataDrivenPriorMatricesCv,
    cisWindow,
    phenotypeFilter,
    residualizationParam,
    genotypeFilter,
    nPCs,
    verbose
) {
    list(
        methodList = norm$methodList,
        fineMappingResult = fineMappingResult,
        cache = twasWeights,
        dataDrivenPriorMatricesCv = dataDrivenPriorMatricesCv,
        cisWindow = cisWindow,
        region = xRegions[[bi]],
        regionIndex = bi,
        nRegions = length(xRegions),
        phenotypeFilter = phenotypeFilter,
        residualizationParam = residualizationParam,
        genotypeFilter = genotypeFilter,
        nPCs = nPCs,
        verbose = verbose
    )
}

# Merge per-region results (NULL when none; single passthrough; else by key).
# @noRd
.twasMergeRegionResults <- function(perRegion, labs) {
    if (length(perRegion) == 0L) {
        return(NULL)
    }
    if (length(perRegion) == 1L) {
        return(perRegion[[1L]])
    }
    .twasMergeResultsByKey(perRegion, labs)
}

# Combine the per-tuple result with any joint result (error if both empty).
# @noRd
.twasQdsAssemble <- function(tw, jointResult) {
    if (is.null(tw) && is.null(jointResult)) {
        msg <- glue(
            "twasWeightsPipeline(QtlDataset): no (context, trait) pair ",
            "produced any weights."
        )
        abort(msg)
    }
    if (is.null(tw)) {
        return(jointResult)
    }
    if (is.null(jointResult)) {
        return(tw)
    }
    .rbindTwasWeights(tw, jointResult, ldReference = NULL)
}

#' @rdname twasWeightsPipeline
#' @export
setMethod(
    "twasWeightsPipeline",
    "QtlSumStats",
    function(
        data,
        methods = NULL,
        contexts = NULL,
        traitId = NULL,
        cisWindow = NULL,
        jointSpecification = NULL,
        fineMappingResult = NULL,
        twasWeights = NULL,
        fitRetention = c("slim", "none", "full"),
        weightsDataType = NULL,
        verbose = 1L,
        panelFilterParam = PanelFilterParam(),
        crossValidationParam = CrossValidationParam()
    ) {
        fitRetention <- arg_match(fitRetention)
        .twasWarnUnretainableDetail(fitRetention, methods, "QtlSumStats")
        .cvRefuseOnSumstats(
            crossValidationParam,
            "twasWeightsPipeline",
            "QtlSumStats"
        )
        build <- .inputsBuild(
            list(data, fineMappingResult, twasWeights),
            "twasWeightsPipeline"
        )
        .fitCoverageCheckedRows(fineMappingResult, "twasWeightsPipeline()")
        res <- .twasPipelineQtlSumStats(
            data = data,
            methods = methods,
            contexts = contexts,
            traitId = traitId,
            cisWindow = cisWindow,
            jointSpecification = jointSpecification,
            fineMappingResult = fineMappingResult,
            twasWeights = twasWeights,
            fitRetention = fitRetention,
            weightsDataType = weightsDataType,
            verbose = verbose,
            panelFilterParam = panelFilterParam
        )
        .withGenome(res, build)
    }
)

# ---- QtlSumStats pipeline worker + phase helpers ---------------------------

# Univariate + multivariate per-tuple rows for the QtlSumStats path. `part`
# carries the selected rows and the token split from
# .twasQssSelectAndPartition(); `cfg` the settings both builders share.
# @noRd
.twasQssBuildRows <- function(part, methodArgs, cfg) {
    rowArgs <- list(
        part$selRows,
        studyCol = part$studyCol,
        contextCol = part$contextCol,
        traitCol = part$traitCol,
        methodArgs = methodArgs
    )
    c(
        exec(
            .twasQssUnivariateRows,
            !!!rowArgs,
            univariateTokens = part$univariateTokens,
            data = cfg$data,
            twasWeights = cfg$twasWeights,
            weightsDataType = cfg$dataType,
            fitRetention = cfg$fitRetention,
            fineMappingResult = cfg$fineMappingResult,
            panelFilterParam = cfg$panelFilterParam
        ),
        # `twasWeights` is the resume cache. The univariate builder looks
        # up one (study, context, trait, method) row; the multivariate one
        # asks for a whole (study, trait) group and reuses it only when
        # every context's row is present -- see .twasMvCacheHits.
        exec(
            .twasQssMultivariateRows,
            !!!rowArgs,
            multivariateTokens = part$multivariateTokens,
            twasWeights = cfg$twasWeights,
            data = cfg$data,
            weightsDataType = cfg$dataType,
            fitRetention = cfg$fitRetention,
            fineMappingResult = cfg$fineMappingResult,
            panelFilterParam = cfg$panelFilterParam
        )
    )
}

# A QtlSumStats TWAS-weights run's inputs: its data, checked for a QC
# record, and its selection parsed once, with `cisWindow` laid on the tuples
# selected.
# @noRd
.twasQssInputs <- function(data, contexts, traitId, cisWindow) {
    .twasAssertQcd(data)
    .qssSelection(data, contexts, traitId, cisWindow, "twasWeightsPipeline")
}

.twasPipelineQtlSumStats <- function(
    data,
    methods,
    contexts,
    traitId,
    cisWindow,
    jointSpecification,
    fineMappingResult,
    twasWeights,
    fitRetention,
    weightsDataType,
    verbose,
    panelFilterParam
) {
    inputs <- .twasQssInputs(data, contexts, traitId, cisWindow)
    # One record of the settings the joint phase and the shared dispatcher
    # both need, so each names them once.
    jointCfg <- list(
        data = inputs$data,
        contexts = inputs$contexts,
        traitId = inputs$traitId,
        weightsDataType = weightsDataType,
        verbose = verbose,
        fitRetention = fitRetention,
        panelFilterParam = panelFilterParam
    )
    joint <- .twasQssResolveTokens(
        jointSpecification,
        methods,
        fineMappingResult,
        jointCfg
    )
    if (joint$done) {
        return(joint$result)
    }
    .twasQssAfterJoint(
        inputs$data,
        joint = joint,
        contexts = inputs$contexts,
        traitId = inputs$traitId,
        twasWeights = twasWeights,
        weightsDataType = weightsDataType,
        fineMappingResult = fineMappingResult,
        fitRetention = fitRetention,
        panelFilterParam = panelFilterParam
    )
}

# The per-tuple loop, after the joint phase has consumed the tokens it
# owns (mrmash); `joint` carries the remainder, their options, and
# whatever the joint phase itself produced.
# @noRd
.twasQssAfterJoint <- function(
    data,
    joint,
    contexts,
    traitId,
    twasWeights,
    weightsDataType,
    fineMappingResult,
    fitRetention,
    panelFilterParam
) {
    part <- .twasQssSelectAndPartition(
        data,
        joint$tokens,
        contexts,
        traitId
    )
    rows <- .twasQssBuildRows(
        part,
        joint$methodArgs,
        list(
            data = data,
            twasWeights = twasWeights,
            dataType = weightsDataType,
            fineMappingResult = fineMappingResult,
            fitRetention = fitRetention,
            panelFilterParam = panelFilterParam
        )
    )
    .twasQssAssemble(rows, joint$result, data)
}

# Normalize the methods argument into (tokens, methodArgs). The default set
# excludes fine-mapping methods; those must be requested explicitly together
# with a fine-mapping result passed via `fineMappingResult`.
# @noRd
.twasSumStatsMethodTokens <- function(methods) {
    if (is.null(methods)) {
        tokens <- c("lasso", "prsCs", "dprGibbs")
        return(list(tokens = tokens, methodArgs = .twasEmptyMethodArgs(tokens)))
    }
    if (is.character(methods)) {
        return(list(
            tokens = methods,
            methodArgs = .twasEmptyMethodArgs(methods)
        ))
    }
    if (is(methods, "MethodsSelectionParam") || is.list(methods)) {
        # This resolver only ever serves the summary-statistics path, so
        # the path is not a parameter -- it is the function's identity.
        return(.methodsParamResolve(
            .methodsParamFor(
                methods,
                "QtlSumStats",
                "TwasWeightsMethodsParam",
                "twasWeightsPipeline"
            ),
            "QtlSumStats"
        ))
    }
    msg <- glue(
        "`methods` must be NULL, a character vector, a named list of ",
        "per-method options, or a TwasWeightsMethodsParam() record."
    )
    abort(msg)
}

# Empty-args map keyed by method token.
# @noRd
.twasEmptyMethodArgs <- function(tokens) {
    set_names(rep(list(list()), length(tokens)), tokens)
}


# Joint-specification dispatch for mrmash. Returns list(done, result, tokens,
# methodArgs): `done` requests an early return with `result`; otherwise the
# remaining (non-mrmash) tokens + args continue through the per-tuple loop.
# @noRd
.twasQssJointPhase <- function(parsedJointSpec, tokens, methodArgs, cfg) {
    if (length(parsedJointSpec) == 0L) {
        return(list(
            done = FALSE,
            result = NULL,
            tokens = tokens,
            methodArgs = methodArgs
        ))
    }
    jointResult <- .twasDispatchJointSpecsQtlSumStats(
        parsedJointSpec,
        tokens,
        methodArgs,
        cfg
    )
    keep <- setdiff(tokens, "mrmash")
    if (length(keep) == 0L) {
        if (is.null(jointResult)) {
            abort("twasWeightsPipeline(QtlSumStats): no joint fits produced.")
        }
        return(list(done = TRUE, result = jointResult))
    }
    list(
        done = FALSE,
        result = jointResult,
        tokens = keep,
        methodArgs = methodArgs[keep]
    )
}

# Normalize the method tokens, run the capability / fine-mapping gates, then
# hand off to the joint phase. Mirrors .fmQssResolveTokens.
# @noRd
.twasQssResolveTokens <- function(
    jointSpecification,
    methods,
    fineMappingResult,
    cfg
) {
    parsedJointSpec <- parseJointSpecification(jointSpecification, cfg$data)
    tm <- .twasSumStatsMethodTokens(methods)
    .twasCheckMethodCapabilities(tm$tokens, "QtlSumStats")
    .twasCheckMethodArgsForInput(tm$methodArgs, "QtlSumStats")
    .twasCheckFineMappingMethods(tm$tokens, fineMappingResult, "QtlSumStats")
    .twasQssJointPhase(parsedJointSpec, tm$tokens, tm$methodArgs, cfg)
}

# Resolve the selected rows, partition tokens into univariate vs multivariate,
# attach the LD sketch, and enforce the multivariate >=2-contexts rule. Returns
# the updated parameter bundle.
# @noRd
.twasQssSelectAndPartition <- function(data, tokens, contexts, traitId) {
    studyCol <- as.character(data$study)
    contextCol <- as.character(data$context)
    traitCol <- as.character(data$trait)
    selRows <- .twasQssSelectRows(data, contexts, traitId)
    isMv <- map_lgl(tokens, .twasIsMultivariateToken)
    multivariateTokens <- tokens[isMv]
    univariateTokens <- tokens[!isMv]
    .twasQssCheckMultivariate(
        multivariateTokens,
        selRows,
        studyCol,
        contextCol,
        traitCol
    )
    list(
        studyCol = studyCol,
        contextCol = contextCol,
        traitCol = traitCol,
        selRows = selRows,
        multivariateTokens = multivariateTokens,
        univariateTokens = univariateTokens
    )
}

# Row indices matching the contexts / traitId filters (error if none).
# @noRd
.twasQssSelectRows <- function(data, contexts, traitId) {
    selRows <- which(.sumStatsRowsSelected(data, contexts, traitId))
    if (length(selRows) == 0L) {
        msg <- glue(
            "twasWeightsPipeline(QtlSumStats): no entries matched the ",
            "supplied contexts / traitId filters."
        )
        abort(msg)
    }
    selRows
}

# Multivariate methods require at least two contexts within some (study, trait).
# @noRd
.twasQssCheckMultivariate <- function(
    multivariateTokens,
    selRows,
    studyCol,
    contextCol,
    traitCol
) {
    if (length(multivariateTokens) == 0L) {
        return(invisible(NULL))
    }
    groupKey <- str_c(
        studyCol[selRows],
        traitCol[selRows],
        sep = "||"
    )
    perGroupNCtx <- map_int(split(contextCol[selRows], groupKey), length)
    if (all(perGroupNCtx < 2L)) {
        mvStr <- str_flatten(multivariateTokens, ", ")
        msg <- glue(
            "twasWeightsPipeline(QtlSumStats): multivariate method(s) ",
            "{mvStr} require at least two contexts per (study, trait); the ",
            "supplied collection has only one context per trait."
        )
        abort(msg)
    }
    invisible(NULL)
}

# ---- Shared row-record helpers ---------------------------------------------

# A single TwasWeights row: (study, context, trait, method) + its entry.
# @noRd
.twasRowRecord <- function(studyName, context, trait, method, entry) {
    list(
        study = studyName,
        context = context,
        trait = trait,
        method = method,
        entry = entry
    )
}

# Assemble a TwasWeights collection from a flat list of row records, each row
# on the LD reference its record carries (none for an individual-level fit),
# and on the trait position and window `data` holds for its (study, context,
# trait) -- the window a cisWindow laid, which the result promises to record.
# @noRd
.twasRowsToWeights <- function(rows, data) {
    if (length(rows) == 0L) {
        return(NULL)
    }
    studies <- map_chr(rows, "study")
    contexts <- map_chr(rows, "context")
    traits <- map_chr(rows, "trait")
    TwasWeights(
        studyName = studies,
        context = contexts,
        trait = traits,
        method = map_chr(rows, "method"),
        entry = map(rows, "entry"),
        traitPos = .anchorVector(data, studies, contexts, traits, "traitPos"),
        region = .anchorVector(data, studies, contexts, traits, "region"),
        ldReference = map(rows, "ldReference")
    )
}

# Resolve the (weight function, fine-mapping adapter) for a method token.
# @noRd
.twasResolveWeightFn <- function(tk) {
    adapter <- .twasFineMappingMethodAdapters[[tk]]
    fn <- if (!is.null(adapter)) {
        adapter$rssWeightFn
    } else {
        .twasMethodCapabilities[[tk]]$sumstatImpl
    }
    list(fn = fn, adapter = adapter)
}

# User kwargs for a token (empty list when unset).
# @noRd
.twasUserArgs <- function(methodArgs, tk) {
    userArgs <- methodArgs[[tk]]
    if (is.null(userArgs)) list() else userArgs
}

# Retain-fit defaults for the mr.mash producer (only when tk == "mrmash" and it
# has no fine-mapping adapter); respects explicit caller overrides.
# @noRd
.twasMrmashRetainDefaults <- function(userArgs, adapter, tk, fitRetention) {
    if (!is.null(adapter) || tk != "mrmash") {
        return(userArgs)
    }
    # mr.mash is the producer of the mvSuSiE prior payload, so it retains
    # even when the run as a whole does not -- but at the level asked for.
    level <- if (identical(fitRetention, "none")) "slim" else fitRetention
    list_assign(
        userArgs,
        fitRetention = userArgs$fitRetention %||% level
    )
}

# The pipeline's fitRetention for a summary-statistics weight function that
# takes one, unless the caller set it for that method; a function without the
# argument keeps no fit. "none" is passed on too, unlike on the
# individual-level path: the fine-mapping adapters here keep a slim fit
# unless told otherwise, so leaving it unset would not drop it.
# @noRd
.twasQssRetainArgs <- function(args, fnName, fitRetention) {
    if (!is_in("fitRetention", names(formals(fnName)))) {
        return(args)
    }
    list_assign(args, fitRetention = args$fitRetention %||% fitRetention)
}

# Run a weight function, warning (with `errPrefix`) and returning NULL on error.
# @noRd
#' @importFrom rlang try_fetch
.twasTryWeights <- function(fn, stat, ldReference, userArgs, errPrefix) {
    try_fetch(
        {
            wfn <- get(fn, mode = "function")
            wArgs <- .twasWeightCallArgs(
                fn,
                list(stat = stat, ldReference = ldReference),
                userArgs
            )
            exec(wfn, !!!wArgs)
        },
        error = function(cnd) {
            warn(errPrefix, parent = cnd)
            NULL
        }
    )
}

# ---- Univariate dispatch: per (study, context, trait), per method ----------

.twasQssUnivariateRows <- function(
    selRows,
    univariateTokens,
    studyCol,
    contextCol,
    traitCol,
    data,
    twasWeights,
    weightsDataType,
    fitRetention,
    fineMappingResult,
    methodArgs,
    panelFilterParam
) {
    if (length(univariateTokens) == 0L) {
        return(list())
    }
    list_flatten(map(
        selRows,
        .twasQssUnivariateRowsForEntry,
        univariateTokens = univariateTokens,
        studyCol = studyCol,
        contextCol = contextCol,
        traitCol = traitCol,
        data = data,
        twasWeights = twasWeights,
        weightsDataType = weightsDataType,
        fitRetention = fitRetention,
        fineMappingResult = fineMappingResult,
        methodArgs = methodArgs,
        panelFilterParam = panelFilterParam
    ))
}

# Which of this tuple's tokens the resume cache already answers, as rows,
# and which still have to be fitted.
# @noRd
.twasQssCacheSplit <- function(twasWeights, st, ctx, tr, univariateTokens) {
    cacheHits <- .twasResolveCacheHits(
        twasWeights,
        st,
        ctx,
        tr,
        univariateTokens
    )
    list(
        rows = imap(
            cacheHits,
            .twasCachedRowRecord,
            st = st,
            ctx = ctx,
            tr = tr
        ),
        toFit = setdiff(univariateTokens, names(cacheHits))
    )
}

# Cached + freshly-fitted rows for one sumstats entry. Resume cache: pull cached
# entries up front and reduce the per-entry fit work to the un-cached tokens.
# @noRd
.twasQssUnivariateRowsForEntry <- function(
    i,
    univariateTokens,
    studyCol,
    contextCol,
    traitCol,
    data,
    twasWeights,
    weightsDataType,
    fitRetention,
    fineMappingResult,
    methodArgs,
    panelFilterParam
) {
    st <- studyCol[i]
    ctx <- contextCol[i]
    tr <- traitCol[i]
    cached <- .twasQssCacheSplit(twasWeights, st, ctx, tr, univariateTokens)
    ldReference <- .rowLdRef(data, i)
    cachedRows <- .recordsWithLdRef(unname(cached$rows), ldReference)
    toFit <- cached$toFit
    if (length(toFit) == 0L) {
        return(cachedRows)
    }
    fitCtx <- .twasQssUnivariateFitCtx(
        data,
        st,
        ctx,
        tr,
        ldReference,
        cutoffs = .panelCutoffs(panelFilterParam)
    )
    fitted <- compact(map(
        toFit,
        .twasQssUnivariateFitOne,
        st = st,
        ctx = ctx,
        tr = tr,
        fitCtx = fitCtx,
        methodArgs = methodArgs,
        fitRetention = fitRetention,
        fineMappingResult = fineMappingResult,
        weightsDataType = weightsDataType
    ))
    c(cachedRows, .recordsWithLdRef(fitted, ldReference))
}

# Cache hits for ONE multivariate group: a named list of the group's rows,
# or NULL when any of them is absent.
#
# A multivariate fit spans every context in its (study, trait) group and
# emits one row each, so it can only be resumed when the cache holds ALL of
# them. A partial hit is not usable: reusing some contexts' weights while
# refitting the others would mix two fits inside one group, and the whole
# point of a multivariate method is that the contexts were fitted together.
#
# This is why the resume cache used to skip multivariate rows entirely --
# `.twasCacheLookup` keys on a single trait/context, so there was no
# per-group question to ask. There is; it is just all-or-nothing.
# @noRd
.twasMvCacheHits <- function(twasWeights, st, tr, ctxNames, tk) {
    if (is.null(twasWeights)) {
        return(NULL)
    }
    hits <- map(
        ctxNames,
        .twasMvCacheHitOne,
        twasWeights = twasWeights,
        st = st,
        tr = tr,
        tk = tk
    )
    if (any(map_lgl(hits, is.null))) {
        return(NULL)
    }
    set_names(hits, ctxNames)
}

# @noRd
.twasMvCacheHitOne <- function(ctx, twasWeights, st, tr, tk) {
    .twasCacheLookup(twasWeights, st, ctx, tr, tk)
}

# Cache hits for a (study, context, trait): named list token -> cached entry.
# @noRd
.twasResolveCacheHits <- function(twasWeights, st, ctx, tr, tokens) {
    hits <- set_names(
        map(
            tokens,
            .twasCacheLookupTok,
            twasWeights = twasWeights,
            st = st,
            ctx = ctx,
            tr = tr
        ),
        tokens
    )
    compact(hits)
}

# Shared Z/N/varY/LD setup for one univariate entry.
# @noRd
.twasQssUnivariateFitCtx <- function(
    data,
    st,
    ctx,
    tr,
    ldReference,
    cutoffs = NULL
) {
    allDf <- .sumStatsDf(
        data,
        studyName = st,
        context = ctx,
        trait = tr,
        require = c("Z", "N"),
        derive = "zFromBetaSe"
    )
    # Narrowed here, where the ids are produced: z, the LD matrix and the
    # variant names all derive from this one vector, so they stay aligned by
    # construction rather than by three subsetting steps agreeing.
    label <- glue(
        "twasWeightsPipeline(QtlSumStats): study='{st}', ",
        "context='{ctx}', trait='{tr}'"
    )
    df <- allDf[
        .panelKeepMask(allDf$variant_id, ldReference, cutoffs, label),
        ,
        drop = FALSE
    ]
    variantIds <- df$variant_id
    varY <- varY(data, studyName = st, context = ctx, trait = tr) %||% 1
    stat <- list(
        z = df$z,
        n = stats::median(df$N, na.rm = TRUE),
        varY = varY,
        variantNames = variantIds
    )
    # The entry's block of the reference, read once and shared by every
    # method fitted to it.
    list(
        variantIds = variantIds,
        stat = stat,
        ldReference = .ldRefInMemory(
            ldReference,
            variantIds,
            label = "twasWeightsPipeline"
        )
    )
}

# The weight function's arguments for one tuple: the user's args and, when the
# token is a fine-mapping method, the precomputed fit threaded into its
# dedicated *Fit argument. NULL when such a token has no fit to thread.
# @noRd
.twasQssUnivariateArgs <- function(
    spec,
    methodArgs,
    tk,
    st,
    ctx,
    tr,
    fineMappingResult
) {
    baseArgs <- .twasUserArgs(methodArgs, tk)
    if (is.null(spec$adapter)) {
        return(baseArgs)
    }
    fit <- .twasFineMappingFitFor(
        fineMappingResult,
        studyName = st,
        context = ctx,
        trait = tr,
        token = tk
    )
    if (is.null(fit)) {
        .twasWarnNoFitUniv(tk, st, ctx, tr)
        return(NULL)
    }
    list_assign(baseArgs, !!!set_names(list(fit), spec$adapter$rssFitArg))
}

# Fit one univariate method for one entry -> a row record, or NULL on skip.
# @noRd
.twasQssUnivariateFitOne <- function(
    tk,
    st,
    ctx,
    tr,
    fitCtx,
    methodArgs,
    fitRetention,
    fineMappingResult,
    weightsDataType
) {
    spec <- .twasResolveWeightFn(tk)
    userArgs <- .twasQssUnivariateArgs(
        spec,
        methodArgs,
        tk,
        st,
        ctx,
        tr,
        fineMappingResult
    )
    # NULL means a fine-mapping token whose fit is missing; it has already
    # warned, and there is nothing to fit against.
    if (is.null(userArgs)) {
        return(NULL)
    }
    weights <- .twasTryWeights(
        spec$fn,
        fitCtx$stat,
        fitCtx$ldReference,
        .twasQssRetainArgs(userArgs, spec$fn, fitRetention),
        .twasFitErrUniv(tk, st, ctx, tr)
    )
    if (is.null(weights)) {
        return(NULL)
    }
    .twasQssRowFromWeights(weights, st, ctx, tr, tk, fitCtx, weightsDataType)
}

# The weights come back with the fit hung off them as an attribute; the row
# stores the two separately, so they are split apart here.
# @noRd
.twasQssRowFromWeights <- function(
    weights,
    st,
    ctx,
    tr,
    tk,
    fitCtx,
    weightsDataType
) {
    .twasRowRecord(
        st,
        ctx,
        tr,
        tk,
        TwasWeightsRow(
            variantIds = fitCtx$variantIds,
            weights = as.numeric(`attr<-`(weights, "fit", NULL)),
            methodFits = attr(weights, "fit"),
            cvResult = NULL,
            weightStandardized = TRUE,
            weightsDataType = weightsDataType
        )
    )
}

# Warning for a missing univariate fine-mapping fit.
# @noRd
.twasWarnNoFitUniv <- function(tk, st, ctx, tr) {
    msg <- glue(
        "twasWeightsPipeline: no '{tk}' fit found in fineMappingResult ",
        "for (study={st}, context={ctx}, trait={tr}); skipping."
    )
    warn(msg)
}

# Error-message prefix for a failed univariate weight fit.
# @noRd
.twasFitErrUniv <- function(tk, st, ctx, tr) {
    glue(
        "twasWeightsPipeline: method '{tk}' failed for (study={st}, ",
        "context={ctx}, trait={tr}): ",
        .trim = FALSE
    )
}

# ---- Multivariate dispatch: per (study, trait), all selected contexts ------

.twasQssMultivariateRows <- function(
    selRows,
    multivariateTokens,
    studyCol,
    contextCol,
    traitCol,
    data,
    methodArgs,
    fitRetention,
    fineMappingResult,
    weightsDataType,
    panelFilterParam,
    twasWeights = NULL
) {
    if (length(multivariateTokens) == 0L) {
        return(list())
    }
    groupKey <- str_c(
        studyCol[selRows],
        traitCol[selRows],
        sep = "||"
    )
    groups <- split(selRows, groupKey)
    list_flatten(map(
        groups,
        .twasQssMultivariateGroupRows,
        twasWeights = twasWeights,
        multivariateTokens = multivariateTokens,
        studyCol = studyCol,
        contextCol = contextCol,
        traitCol = traitCol,
        data = data,
        methodArgs = methodArgs,
        fitRetention = fitRetention,
        fineMappingResult = fineMappingResult,
        weightsDataType = weightsDataType,
        panelFilterParam = panelFilterParam
    ))
}

# Multivariate rows for one (study, trait) group across its contexts.
# @noRd
.twasQssMultivariateGroupRows <- function(
    gIdx,
    multivariateTokens,
    studyCol,
    contextCol,
    traitCol,
    data,
    methodArgs,
    fitRetention,
    fineMappingResult,
    weightsDataType,
    panelFilterParam,
    twasWeights = NULL
) {
    if (length(gIdx) < 2L) {
        return(list())
    }
    grp <- .twasQssMultivariateGroup(
        gIdx,
        data,
        studyCol,
        contextCol,
        traitCol,
        panelFilterParam
    )
    rows <- list_flatten(map(
        multivariateTokens,
        .twasQssMultivariateFitOne,
        twasWeights = twasWeights,
        st = grp$st,
        tr = grp$tr,
        ctxNames = grp$ctxNames,
        mvStat = grp$mvStat,
        blockReference = grp$blockReference,
        methodArgs = methodArgs,
        fitRetention = fitRetention,
        fineMappingResult = fineMappingResult,
        weightsDataType = weightsDataType
    ))
    .recordsWithLdRef(rows, grp$ldReference)
}

# What every method fitted to one (study, trait) group shares: the study,
# trait and contexts, their LD reference, the stacked statistics, and the
# group's block of the reference, read once. The contexts are fitted together
# against one LD matrix, so they must share an LD reference.
# @noRd
.twasQssMultivariateGroup <- function(
    gIdx,
    data,
    studyCol,
    contextCol,
    traitCol,
    panelFilterParam
) {
    ldReference <- .jointRowsLdRef(
        data,
        gIdx,
        "twasWeightsPipeline (QtlSumStats, multivariate)"
    )
    ctxNames <- contextCol[gIdx]
    mvStat <- .twasQssMultivariateStat(
        data,
        studyCol[gIdx[[1L]]],
        traitCol[gIdx[[1L]]],
        ctxNames,
        ldReference = ldReference,
        cutoffs = .panelCutoffs(panelFilterParam)
    )
    list(
        st = studyCol[gIdx[[1L]]],
        tr = traitCol[gIdx[[1L]]],
        ctxNames = ctxNames,
        ldReference = ldReference,
        mvStat = mvStat,
        blockReference = .ldRefInMemory(
            ldReference,
            mvStat$variantIds,
            label = "twasWeightsPipeline"
        )
    )
}

# Build the (variants x contexts) Z matrix + per-context N for a group. All
# entries in a (study, trait) group must share an identical variant order after
# summaryStatsQc().
# @noRd
.twasQssMultivariateStat <- function(
    data,
    st,
    tr,
    ctxNames,
    ldReference = NULL,
    cutoffs = NULL
) {
    allDf <- .sumStatsDf(
        data,
        studyName = st,
        context = ctxNames[[1L]],
        trait = tr,
        require = c("Z", "N"),
        derive = "zFromBetaSe"
    )
    # Every context shares one SNP order (asserted below), so filtering the
    # first context's ids filters the group: the Z matrix is built against
    # this vector and each context is checked against it.
    label <- glue(
        "twasWeightsPipeline(QtlSumStats, multivariate): study='{st}', ",
        "trait='{tr}'"
    )
    keep <- .panelKeepMask(allDf$variant_id, ldReference, cutoffs, label)
    firstDf <- allDf[keep, , drop = FALSE]
    variantIds <- firstDf$variant_id
    filled <- .twasQssFillContexts(data, st, tr, ctxNames, variantIds)
    list(
        variantIds = variantIds,
        stat = list(
            z = filled$z,
            n = filled$n,
            variantNames = variantIds
        )
    )
}

# One context's z column and median N, checked against the shared SNP order.
# @noRd
.twasQssContextStats <- function(ctx, data, st, tr, variantIds) {
    d <- .sumStatsDf(
        data,
        studyName = st,
        context = ctx,
        trait = tr,
        require = c("Z", "N"),
        derive = "zFromBetaSe"
    )
    kept <- d[is_in(d$variant_id, variantIds), , drop = FALSE]
    .twasQssCheckSnpOrder(kept$variant_id, variantIds, st, tr)
    list(z = kept$z, n = stats::median(kept$N, na.rm = TRUE))
}

# Read each context's z / N into the shared matrix. Every context must present
# the same SNPs in the same order -- the multivariate fitters index z by
# position, not by name, so a differing order would silently pair one context's
# variant with another's.
#
# Every context shares one SNP order (asserted per context), which is what
# lets the columns simply be laid side by side rather than filled into a
# preallocated matrix.
# @noRd
.twasQssFillContexts <- function(data, st, tr, ctxNames, variantIds) {
    stats <- map(
        ctxNames,
        .twasQssContextStats,
        data = data,
        st = st,
        tr = tr,
        variantIds = variantIds
    )
    list(
        z = matrix(
            unname(list_c(map(stats, "z"))),
            nrow = length(variantIds),
            ncol = length(ctxNames),
            dimnames = list(variantIds, ctxNames)
        ),
        n = set_names(map_dbl(stats, "n"), ctxNames)
    )
}

# @noRd
.twasQssCheckSnpOrder <- function(got, want, st, tr) {
    if (identical(got, want)) {
        return(invisible(NULL))
    }
    abort(glue(
        "twasWeightsPipeline(QtlSumStats, multivariate): every ",
        "entry for (study='{st}', trait='{tr}') must share an ",
        "identical SNP order after summaryStatsQc(). Use the same ",
        "ldReference on every entry."
    ))
}

# The engine arguments for one multivariate token. mvsusie is
# fine-mapping, so its pre-fit is threaded through (mr.mash is not);
# .twasMvThreadFit() answers NULL when that pre-fit is missing.
# @noRd
.twasQssMvArgs <- function(
    spec,
    methodArgs,
    tk,
    st,
    tr,
    ctxNames,
    fitRetention,
    fineMappingResult
) {
    baseArgs <- .twasQssRetainArgs(
        .twasMrmashRetainDefaults(
            .twasUserArgs(methodArgs, tk),
            spec$adapter,
            tk,
            fitRetention
        ),
        spec$fn,
        fitRetention
    )
    if (is.null(spec$adapter)) {
        return(baseArgs)
    }
    .twasMvThreadFit(
        spec,
        baseArgs,
        tk,
        st,
        tr,
        ctxNames,
        fineMappingResult
    )
}

# Fit one multivariate method for a group -> one row record per context (empty
# list on skip). `userArgs` is NULL when the pre-fit a fine-mapping token needs
# is missing.
# @noRd
.twasQssMultivariateFitOne <- function(
    tk,
    st,
    tr,
    ctxNames,
    mvStat,
    blockReference,
    methodArgs,
    fitRetention,
    fineMappingResult,
    weightsDataType,
    twasWeights = NULL
) {
    cached <- .twasMvCacheHits(twasWeights, st, tr, ctxNames, tk)
    if (!is.null(cached)) {
        return(.twasMvCachedRows(cached, st, tr, tk))
    }
    spec <- .twasResolveWeightFn(tk)
    userArgs <- .twasQssMvArgs(
        spec,
        methodArgs,
        tk = tk,
        st = st,
        tr = tr,
        ctxNames = ctxNames,
        fitRetention = fitRetention,
        fineMappingResult = fineMappingResult
    )
    if (is.null(userArgs)) {
        return(list())
    }
    weights <- .twasTryWeights(
        spec$fn,
        mvStat$stat,
        blockReference,
        userArgs,
        .twasFitErrMv(tk, st, tr)
    )
    .twasMvRowsFromWeights(
        weights,
        ctxNames,
        mvStat,
        st,
        tr,
        tk,
        weightsDataType
    )
}

# Per-context row records from one multivariate fit's weights (empty list when
# the fit failed and `weights` is NULL). The fit is hung off the weights as an
# attribute; the rows store the two separately.
# @noRd
.twasMvRowsFromWeights <- function(
    weights,
    ctxNames,
    mvStat,
    st,
    tr,
    tk,
    weightsDataType
) {
    if (is.null(weights)) {
        return(list())
    }
    wMatrix <- if (is.matrix(weights)) weights else as.matrix(weights)
    fitAttr <- attr(wMatrix, "fit")
    bare <- `attr<-`(wMatrix, "fit", NULL)
    .twasMvContextRows(
        bare,
        fitAttr,
        ctxNames,
        mvStat,
        st,
        tr,
        tk,
        weightsDataType
    )
}

# Thread the precomputed fine-mapping fit into a multivariate method's args;
# returns NULL (signalling skip) when the fit is absent.
# @noRd
.twasMvThreadFit <- function(
    spec,
    userArgs,
    tk,
    st,
    tr,
    ctxNames,
    fineMappingResult
) {
    fit <- .twasFineMappingFitFor(
        fineMappingResult,
        studyName = st,
        context = ctxNames[[1L]],
        trait = tr,
        token = tk
    )
    if (is.null(fit)) {
        .twasWarnNoFitMv(tk, st, tr)
        return(NULL)
    }
    list_assign(userArgs, !!!set_names(list(fit), spec$adapter$rssFitArg))
}

# Warning for a missing multivariate fine-mapping fit.
# @noRd
.twasWarnNoFitMv <- function(tk, st, tr) {
    msg <- glue(
        "twasWeightsPipeline: no '{tk}' fit found in fineMappingResult ",
        "for (study={st}, trait={tr}); skipping."
    )
    warn(msg)
}

# Error-message prefix for a failed multivariate weight fit.
# @noRd
.twasFitErrMv <- function(tk, st, tr) {
    glue(
        "twasWeightsPipeline: multivariate method '{tk}' failed for ",
        "(study={st}, trait={tr}): ",
        .trim = FALSE
    )
}

# One row record per context from a fitted multivariate weight matrix. The
# underlying joint fit is shared on the first row only; the remaining rows
# reference it by leaving fits NULL.
# @noRd
.twasMvContextRows <- function(
    weights,
    fitAttr,
    ctxNames,
    mvStat,
    st,
    tr,
    tk,
    weightsDataType
) {
    map(
        seq_along(ctxNames),
        .twasMvContextRow,
        weights = weights,
        fitAttr = fitAttr,
        ctxNames = ctxNames,
        mvStat = mvStat,
        st = st,
        tr = tr,
        tk = tk,
        weightsDataType = weightsDataType
    )
}

# Combine the per-tuple result with any joint result (error if both empty).
# @noRd
.twasQssAssemble <- function(rows, jointResult, data) {
    perTupleResult <- .twasRowsToWeights(rows, data)
    if (is.null(jointResult)) {
        if (is.null(perTupleResult)) {
            msg <- glue(
                "twasWeightsPipeline(QtlSumStats): no entries produced ",
                "weights."
            )
            abort(msg)
        }
        return(perTupleResult)
    }
    if (is.null(perTupleResult)) {
        return(jointResult)
    }
    .rbindTwasWeights(perTupleResult, jointResult)
}

# =============================================================================
# MultiStudyQtlDataset method
# =============================================================================
# Mirrors the fineMappingPipeline(MultiStudyQtlDataset) recursion: iterates
# the embedded individual-level QtlDataset entries, then processes the
# optional embedded QtlSumStats. The result rows from the two phases are
# rbind'd; the joint columns (when populated by either phase) are carried
# through .rbindTwasWeights.

# Per-embedded-study TWAS-weights worker for .multiStudyPipelineDriver: recurse
# twasWeightsPipeline on one QtlDataset. `cfg` bundles the parent's forwarded
# args.
#
# Every setting is named -- not the `...`/dotArgs channel this once spliced,
# which is what made multi-study TWAS inherit the per-study CV default and
# silently drop residualization / fitRetention. A name this entry point does
# not accept is an "unused argument" error; `cfg`'s panelFilterParam belongs
# to the summary-statistics path and is deliberately not forwarded.
# @noRd
.twasPerStudy <- function(qd, cfg) {
    # `cfg` carries the parsed selection for every study; this study gets its
    # own share, and a study left with no trait selected in its selected
    # contexts -- none of its contexts, none of the traits -- sits out.
    if (!.studyHasSelectedTraits(qd, cfg$contexts, cfg$traitId, cfg$region)) {
        return(NULL)
    }
    s <- studyName(qd)
    .muffleFitCoverage(twasWeightsPipeline(
        data = qd,
        jointSpecification = NULL,
        fitRetention = cfg$fitRetention,
        methods = cfg$methods,
        contexts = .selectedContexts(cfg$contexts, s),
        traitId = .studyTraitIdArg(cfg$traitId, s),
        region = cfg$region,
        cisWindow = cfg$cisWindow,
        jointRegions = cfg$jointRegions,
        fineMappingResult = cfg$fineMappingResult,
        twasWeights = cfg$twasWeights,
        phenotypeFilterParam = cfg$phenotypeFilterParam,
        verbose = cfg$verbose,
        crossValidationParam = cfg$crossValidationParam,
        ensembleParam = cfg$ensembleParam,
        seed = cfg$seed,
        genotypeFilterParam = cfg$genotypeFilterParam,
        mashPrior = cfg$mashPrior,
        estimatePi = cfg$estimatePi,
        usePCA = cfg$usePCA,
        nPCs = cfg$nPCs,
        weightsDataType = cfg$weightsDataType,
        residualizationParam = cfg$residualizationParam
    ))
}

# Embedded-sumstats TWAS-weights worker for .multiStudyPipelineDriver.
# @noRd
.twasSumStats <- function(ss, cfg) {
    # Summary statistics carry no genotypes: each tuple's window is in the
    # input, narrowed by `cisWindow` around its trait position as on the
    # QtlDataset side, and there is nothing to residualize, cross-validate,
    # stack or fit on full data. `cfg` carries those for the individual-level
    # sibling and they are not forwarded here.
    #
    # The parsed selection is applied here, to the rows, so the recursion runs
    # on what it keeps with nothing left to filter -- and rows it leaves out
    # are not an error when another component holds the selected tuples.
    keepRows <- .sumStatsRowsSelected(ss, cfg$contexts, cfg$traitId)
    if (!any(keepRows)) {
        return(NULL)
    }
    .muffleFitCoverage(twasWeightsPipeline(
        data = if (all(keepRows)) ss else ss[which(keepRows)],
        jointSpecification = NULL,
        fitRetention = cfg$fitRetention,
        methods = cfg$methods,
        contexts = NULL,
        traitId = NULL,
        cisWindow = cfg$cisWindow,
        fineMappingResult = cfg$fineMappingResult,
        twasWeights = cfg$twasWeights,
        verbose = cfg$verbose,
        panelFilterParam = cfg$panelFilterParam,
        weightsDataType = cfg$weightsDataType
    ))
}

# Few methods distinguish the retention levels (.twasFullDetailTokens());
# every other weight method keeps its fit whole or not at all, so "full" is
# the same as "slim" for them. Erroring when NO requested method can honour
# it, and warning when only some can, is the alternative to dropping it in
# silence. Which methods can depends on the input the run reads.
# @noRd
.twasWarnUnretainableDetail <- function(fitRetention, methods, inputKind) {
    if (!identical(fitRetention, "full")) {
        return(invisible(NULL))
    }
    tokens <- .twasMethodTokensFromArg(methods)
    if (length(tokens) == 0L) {
        return(invisible(NULL))
    }
    honouring <- .twasFullDetailTokens(inputKind)
    honours <- is_in(tokens, honouring)
    if (!any(honours)) {
        abort(glue(
            "twasWeightsPipeline: fitRetention = \"full\" is honoured on ",
            "{inputKind} input only by {str_flatten(honouring, ', ')}; ",
            "{str_flatten(tokens, ', ')} cannot vary the retained detail. ",
            "Use fitRetention = \"slim\"."
        ))
    }
    if (!all(honours)) {
        warn(glue(
            "twasWeightsPipeline: fitRetention = \"full\" applies to ",
            "{str_flatten(tokens[honours], ', ')}; ",
            "{str_flatten(tokens[!honours], ', ')} retain their usual fit."
        ))
    }
    invisible(NULL)
}

# The methods whose "full" fit keeps more than their "slim" one, on one input
# path: mr.mash's whole fit, and on summary statistics the shrinkage-grid
# methods' candidate coefficients (lassosumRssWeights() and the penalized
# RSS weights).
# @noRd
.twasFullDetailTokens <- function(inputKind) {
    if (identical(inputKind, "QtlSumStats")) {
        return(c("mrmash", "lasso", "enet", "ridge", "scad", "mcp", "l0learn"))
    }
    "mrmash"
}

# ---- MultiStudyQtlDataset pipeline worker + phase helpers ------------------

# Method tokens from a `methods` arg: a character vector as-is; a named list ->
# canonical bare tokens; otherwise empty.
# @noRd
.twasMethodTokensFromArg <- function(methods) {
    if (.twasMethodsIsParam(methods)) {
        return(.methodsParamTokens(methods))
    }
    if (is.character(methods)) {
        methods
    } else if (is.list(methods)) {
        str_remove(names(methods), "(_weights|Weights)$")
    } else {
        character(0)
    }
}

# A MethodsSelectionParam keeps its structure through these helpers: the
# token questions read across all three slots, and dropping a method drops
# it from each, so the per-component recursion still receives a Param that
# knows which path each entry is for.
# @noRd
.twasMethodsIsParam <- function(methods) {
    is(methods, "MethodsSelectionParam")
}

# Drop mrmash (handled by the joint dispatcher) from a `methods` arg.
# @noRd
.twasMsStripMrmash <- function(methods) {
    if (.twasMethodsIsParam(methods)) {
        return(.methodsParamDrop(methods, "mrmash"))
    }
    if (is.character(methods)) {
        setdiff(methods, "mrmash")
    } else if (is.list(methods)) {
        methods[.twasMethodTokensFromArg(methods) != "mrmash"]
    } else {
        methods
    }
}

# TRUE when a character/list `methods` arg has become empty.
# @noRd
.twasMethodsEmpty <- function(methods) {
    # A Param is neither a character vector nor a list, so it has to be
    # asked directly: with every slot empty it names no methods, which is
    # what the joint phase produces after stripping mrmash.
    if (.twasMethodsIsParam(methods)) {
        return(length(.methodsParamTokens(methods)) == 0L)
    }
    (is.character(methods) || is.list(methods)) && length(methods) == 0L
}

# The MultiStudyQtlDataset method of twasWeightsPipeline(). `jointCfg` is the
# QtlDataset method's own record (.twasJointCfg()) plus the collection and the
# two filter settings the components apply, so each component's joint fits
# run under the settings its own QtlDataset run would use.
# @noRd
.twasPipelineMultiStudy <- function(
    data,
    methods = "default",
    contexts = NULL,
    traitId = NULL,
    region = NULL,
    cisWindow = NULL,
    genotypeFilterParam = GenotypeFilterParam(),
    panelFilterParam = PanelFilterParam(),
    jointRegions = FALSE,
    jointSpecification = NULL,
    fineMappingResult = NULL,
    twasWeights = NULL,
    mashPrior = NULL,
    estimatePi = TRUE,
    usePCA = FALSE,
    nPCs = 10L,
    fitRetention = c("slim", "none", "full"),
    weightsDataType = NULL,
    phenotypeFilterParam = PhenotypeFilterParam(),
    verbose = 1,
    residualizationParam = ResidualizationParam(),
    crossValidationParam = CrossValidationParam(),
    ensembleParam = EnsembleParam(),
    seed = NULL
) {
    fitRetention <- arg_match(fitRetention)
    # Every collection holds individual-level studies, so their rule applies
    # here; a summary-statistic study is checked again by its own route.
    .twasWarnUnretainableDetail(fitRetention, methods, "QtlDataset")
    # Resolved once here so a bad choice fails before any study runs; the
    # record carries the bundle itself, as each component resolves its own.
    .twasResolvePhenotypeFilter(phenotypeFilterParam)
    build <- .inputsBuild(
        list(data, region, fineMappingResult, twasWeights),
        "twasWeightsPipeline(MultiStudyQtlDataset)"
    )
    # Once here; each study's run below reads the same rows and is muffled.
    .fitCoverageCheckedRows(fineMappingResult, "twasWeightsPipeline()")
    inputs <- .msResolveInputs(
        region,
        cisWindow,
        jointRegions,
        "twasWeightsPipeline(MultiStudyQtlDataset)"
    )
    cfg <- .twasMsCfg(
        methods = methods,
        contexts = parseContexts(contexts, data),
        traitId = parseTraitIds(traitId, data),
        region = region,
        cisWindow = cisWindow,
        jointRegions = jointRegions,
        fineMappingResult = fineMappingResult,
        twasWeights = twasWeights,
        phenotypeFilterParam = phenotypeFilterParam,
        verbose = verbose,
        crossValidationParam = crossValidationParam,
        ensembleParam = ensembleParam,
        seed = seed,
        genotypeFilterParam = genotypeFilterParam,
        panelFilterParam = panelFilterParam,
        mashPrior = mashPrior,
        estimatePi = estimatePi,
        usePCA = usePCA,
        nPCs = nPCs,
        fitRetention = fitRetention,
        weightsDataType = weightsDataType,
        residualizationParam = residualizationParam
    )
    .msAssertTraitsSelected(data, cfg, "twasWeightsPipeline")
    jointCfg <- .twasMsJointCfg(cfg, data, inputs$xRegions)
    joint <- .twasMsResolveTokens(
        jointSpecification,
        methods,
        fineMappingResult,
        jointCfg
    )
    if (joint$done) {
        return(.withGenome(joint$result, build))
    }
    .withGenome(
        .twasMsDriver(
            data,
            joint$result,
            list_assign(cfg, methods = joint$methods)
        ),
        build
    )
}

# The MultiStudyQtlDataset joint record, derived from the record the
# per-component recursion receives (`cfg`, .twasMsCfg()): the settings each
# component's own QtlDataset run fits under (.twasJointCfg()), plus the
# collection and the two filters the components apply. Deriving it, rather
# than re-listing the settings, keeps the per-study joint fits on exactly the
# settings of the per-study fits.
# @noRd
.twasMsJointCfg <- function(cfg, data, xRegions) {
    list_assign(
        .twasJointCfg(
            genotypeFilter = .qtlResolveFilter(cfg$genotypeFilterParam),
            phenotypeFilter = .twasResolvePhenotypeFilter(
                cfg$phenotypeFilterParam
            ),
            contexts = cfg$contexts,
            traitId = cfg$traitId,
            cisWindow = cfg$cisWindow,
            xRegions = xRegions,
            fineMappingResult = cfg$fineMappingResult,
            mashPrior = cfg$mashPrior,
            crossValidationParam = cfg$crossValidationParam,
            ensembleParam = cfg$ensembleParam,
            estimatePi = cfg$estimatePi,
            fitRetention = cfg$fitRetention,
            residualizationParam = cfg$residualizationParam,
            weightsDataType = cfg$weightsDataType,
            verbose = cfg$verbose,
            seed = cfg$seed
        ),
        data = data,
        genotypeFilterParam = cfg$genotypeFilterParam,
        panelFilterParam = cfg$panelFilterParam
    )
}

#' @rdname twasWeightsPipeline
#' @export
setMethod(
    "twasWeightsPipeline",
    "MultiStudyQtlDataset",
    .twasPipelineMultiStudy
)


# Joint-specification dispatch for mrmash. Returns list(done, result, methods)
# where `done` requests an early return with `result` and `methods` is the
# mrmash-stripped set for the per-component recursion.
# @noRd
.twasMsJointPhase <- function(parsedJointSpec, methods, cfg) {
    if (length(parsedJointSpec) == 0L) {
        # No joint spec: every token is fitted per component, so each must
        # find a multivariate Y there.
        .twasCheckMultivariateComponents(
            .twasMethodTokensFromArg(methods),
            cfg$data
        )
        return(list(done = FALSE, result = NULL, methods = methods))
    }
    jointResult <- .twasDispatchJointSpecsMultiStudy(
        parsedJointSpec,
        .twasMethodTokensFromArg(methods),
        methods,
        cfg
    )
    stripped <- .twasMsStripMrmash(methods)
    if (.twasMethodsEmpty(stripped)) {
        if (is.null(jointResult)) {
            msg <- glue(
                "twasWeightsPipeline(MultiStudyQtlDataset): no joint fits ",
                "produced."
            )
            abort(msg)
        }
        return(list(done = TRUE, result = jointResult))
    }
    # Whatever the joint phase did not handle falls to the per-component
    # recursion, and is held to the same requirement.
    .twasCheckMultivariateComponents(
        .twasMethodTokensFromArg(stripped),
        cfg$data
    )
    list(done = FALSE, result = jointResult, methods = stripped)
}

# Gate fine-mapping methods early so the recursion into the embedded
# QtlDataset / QtlSumStats components doesn't re-run fine-mapping, then hand
# off to the joint phase. Mirrors .fmMsResolveTokens.
# @noRd
.twasMsResolveTokens <- function(
    jointSpecification,
    methods,
    fineMappingResult,
    cfg
) {
    parsedJointSpec <- parseJointSpecification(jointSpecification, cfg$data)
    # Translated here, where the dataset is in hand, so the Param is what
    # travels into the per-component recursion; each component reads the
    # slot for its own input class.
    methods <- .methodsParamForMulti(
        methods,
        "TwasWeightsMethodsParam",
        "twasWeightsPipeline",
        !is.null(sumStats(cfg$data))
    )
    .twasCheckFineMappingMethods(
        .twasMethodTokensFromArg(methods),
        fineMappingResult,
        "MultiStudyQtlDataset",
        fitFullData = .cvResolve(cfg$crossValidationParam)$fitFullData
    )
    .twasMsJointPhase(parsedJointSpec, methods, cfg)
}

# Run the per-study / per-component recursion via the shared multi-study driver.
# `cfg` is the `.twasMsCfg()` record of the settings each component re-runs
# with.
# @noRd
.twasMsDriver <- function(data, jointResult, cfg) {
    .multiStudyPipelineDriver(
        data,
        jointResult,
        .twasPerStudy,
        .twasSumStats,
        cfg,
        .rbindTwasWeights,
        TwasWeights,
        "twasWeightsPipeline",
        noun = "weights"
    )
}

# The settings the per-study / per-component recursion re-runs each component
# with, as the record .multiStudyPipelineDriver() hands to .twasPerStudy /
# .twasSumStats.
# @noRd
.twasMsCfg <- function(
    methods,
    contexts,
    traitId,
    region,
    cisWindow,
    jointRegions,
    fineMappingResult,
    twasWeights,
    phenotypeFilterParam,
    verbose,
    crossValidationParam,
    ensembleParam,
    seed,
    genotypeFilterParam,
    panelFilterParam,
    mashPrior,
    estimatePi,
    usePCA,
    nPCs,
    fitRetention,
    weightsDataType,
    residualizationParam
) {
    list(
        methods = methods,
        contexts = contexts,
        traitId = traitId,
        region = region,
        cisWindow = cisWindow,
        jointRegions = jointRegions,
        fineMappingResult = fineMappingResult,
        twasWeights = twasWeights,
        phenotypeFilterParam = phenotypeFilterParam,
        verbose = verbose,
        crossValidationParam = crossValidationParam,
        ensembleParam = ensembleParam,
        seed = seed,
        genotypeFilterParam = genotypeFilterParam,
        panelFilterParam = panelFilterParam,
        mashPrior = mashPrior,
        estimatePi = estimatePi,
        usePCA = usePCA,
        nPCs = nPCs,
        fitRetention = fitRetention,
        weightsDataType = weightsDataType,
        residualizationParam = residualizationParam
    )
}

# =============================================================================
# SR-TWAS ensemble stacking solvers (used by ensembleWeights, the primitive the
# engine's .twasEnsembleLayer calls per context)
# =============================================================================

# Solve ensemble stacking via quadprog (constrained QP with sum-to-1 and
# non-negativity).
# @param Pvalid Matrix of CV predictions for valid methods (n x Kvalid).
# @param yObs Observed outcome vector (n).
# @param Kvalid Number of valid methods.
# @return Normalized coefficient vector of length Kvalid.
# @noRd
.solveEnsembleQuadprog <- function(Pvalid, yObs, Kvalid) {
    if (!requireNamespace("quadprog", quietly = TRUE)) {
        abort("Package 'quadprog' is required for solver='quadprog'.")
    }

    gram <- crossprod(Pvalid)
    dvec <- as.vector(crossprod(Pvalid, yObs))
    # Ridge term for numerical stability (small relative to trace)
    Dmat <- gram + 1e-8 * mean(diag(gram)) * diag(Kvalid)

    # Constraint matrix: first constraint is equality (sum = 1), then Kvalid
    # non-negativity constraints.
    Amat <- cbind(rep(1, Kvalid), diag(Kvalid))
    bvec <- c(1, rep(0, Kvalid))

    qpSol <- try_fetch(
        solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1),
        error = function(cnd) {
            msg <- glue(
                "QP solver failed. Falling back to equal weights among ",
                "valid methods."
            )
            warn(msg, parent = cnd)
            NULL
        }
    )

    if (is.null(qpSol)) {
        return(rep(1 / Kvalid, Kvalid))
    }

    # Numerical cleanup: clamp to non-negative and renormalize
    zetaValid <- pmax(qpSol$solution, 0)
    zetaSum <- sum(zetaValid)
    if (zetaSum <= 0) {
        warn("QP returned all-zero solution. Falling back to equal weights.")
        return(rep(1 / Kvalid, Kvalid))
    }
    zetaValid / zetaSum
}

# Solve ensemble stacking via NNLS (non-negative least squares, then normalize).
# This is the approach used by SuperLearner (Lawson-Hanson algorithm).
# @param Pvalid Matrix of CV predictions for valid methods (n x Kvalid).
# @param yObs Observed outcome vector (n).
# @param Kvalid Number of valid methods.
# @return Normalized coefficient vector of length Kvalid.
# @noRd
.solveEnsembleNnls <- function(Pvalid, yObs, Kvalid) {
    # Engine-call exception: nnls solves for the ensemble's stacking weights
    # over the per-method prediction matrix. It is an optimisation primitive
    # used by pipeline logic, not an interface to a TWAS weight-learning
    # method, so it does not belong in a *Wrapper.R file.
    if (!requireNamespace("nnls", quietly = TRUE)) {
        abort("Package 'nnls' is required for solver='nnls'.")
    }

    fit <- try_fetch(
        nnls::nnls(Pvalid, yObs),
        error = function(cnd) {
            msg <- "NNLS solver failed. Falling back to equal weights."
            warn(msg, parent = cnd)
            NULL
        }
    )

    if (is.null(fit)) {
        return(rep(1 / Kvalid, Kvalid))
    }

    zetaValid <- fit$x
    zetaSum <- sum(zetaValid)
    if (zetaSum <= 0) {
        warn(
            "NNLS returned all-zero solution. Falling back to equal weights."
        )
        return(rep(1 / Kvalid, Kvalid))
    }
    zetaValid / zetaSum
}

# Ensemble stacking objective: the sum of squared residuals of the stacked
# prediction `Pvalid %*% z`. `problem` is the record .solveEnsembleLbfgsb()
# builds once -- the predictions, the outcome and their cross-products -- and
# hands optim() as its single extra argument, so the objective and its
# gradient read one shared definition, each taking the part it needs.
# @noRd
.ensembleObj <- function(z, problem) {
    sum((problem$yObs - problem$Pvalid %*% z)^2)
}

# Gradient of the stacking objective, 2 (P'P z - P'y), from the
# cross-products in `problem`.
# @noRd
.ensembleGrad <- function(z, problem) {
    as.vector(2 * (problem$PtP %*% z - problem$Pty))
}

# Solve ensemble stacking via L-BFGS-B (box-constrained optimization, then
# normalize). Uses base R optim() with analytical gradient. No extra
# dependencies.
# @param Pvalid Matrix of CV predictions for valid methods (n x Kvalid).
# @param yObs Observed outcome vector (n).
# @param Kvalid Number of valid methods.
# @return Normalized coefficient vector of length Kvalid.
# @noRd
.solveEnsembleLbfgsb <- function(Pvalid, yObs, Kvalid) {
    problem <- list(
        Pvalid = Pvalid,
        yObs = yObs,
        PtP = crossprod(Pvalid),
        Pty = as.vector(crossprod(Pvalid, yObs))
    )
    fit <- try_fetch(
        optim(
            par = rep(1 / Kvalid, Kvalid),
            fn = .ensembleObj,
            gr = .ensembleGrad,
            problem = problem,
            method = "L-BFGS-B",
            lower = rep(0, Kvalid)
        ),
        error = function(cnd) {
            msg <- glue(
                "L-BFGS-B solver failed. Falling back to equal weights."
            )
            warn(msg, parent = cnd)
            NULL
        }
    )

    if (is.null(fit)) {
        return(rep(1 / Kvalid, Kvalid))
    }

    zetaValid <- pmax(fit$par, 0)
    zetaSum <- sum(zetaValid)
    if (zetaSum <= 0) {
        msg <- glue(
            "L-BFGS-B returned all-zero solution. Falling back to ",
            "equal weights."
        )
        warn(msg)
        return(rep(1 / Kvalid, Kvalid))
    }
    zetaValid / zetaSum
}

# Solve ensemble stacking via glmnet (penalized regression with non-negativity).
# Uses cv.glmnet for automatic lambda selection. The alpha parameter controls
# the elastic net mixing: alpha=1 is lasso (sparse), alpha=0 is ridge.
# @param Pvalid Matrix of CV predictions for valid methods (n x Kvalid).
# @param yObs Observed outcome vector (n).
# @param Kvalid Number of valid methods.
# @param alpha Elastic net mixing parameter (default 1 = lasso).
# @return Normalized coefficient vector of length Kvalid.
# @noRd
.solveEnsembleGlmnet <- function(Pvalid, yObs, Kvalid, alpha = 1) {
    # Here glmnet is a constrained-regression solver fitting the ensemble's
    # stacking weights -- a different role from weight learning, but the same
    # engine, so the call goes through fitGlmnetCv() like every other. Only
    # the stacking logic is pipeline-side. `lower.limits = 0` is what makes
    # the regression non-negative; it has no pecotmr formal, so it travels in
    # the engine's own options bag. `standardize` is left at glmnet's default
    # because the columns of `Pvalid` are already on one scale: they are all
    # predictions of the same outcome.
    fit <- try_fetch(
        fitGlmnetCv(
            x = Pvalid,
            y = yObs,
            alpha = alpha,
            intercept = FALSE,
            methodArgs = GlmnetOptions(lower.limits = 0)
        ),
        error = function(cnd) {
            msg <- glue(
                "glmnet solver failed. Falling back to equal weights."
            )
            warn(msg, parent = cnd)
            NULL
        }
    )

    if (is.null(fit)) {
        return(rep(1 / Kvalid, Kvalid))
    }

    # [-1] drops the intercept.
    zetaValid <- pmax(as.numeric(coef(fit, s = "lambda.min"))[-1], 0)
    zetaSum <- sum(zetaValid)
    if (zetaSum <= 0) {
        warn(
            "glmnet returned all-zero solution. Falling back to equal weights."
        )
        return(rep(1 / Kvalid, Kvalid))
    }
    zetaValid / zetaSum
}

#' Ensemble TWAS Weights via Stacked Regression
#'
#' Given cross-validated predictions from multiple TWAS weight methods, learns
#' non-negative combination coefficients (summing to 1) via constrained least
#' squares. Returns ensemble weights and per-method performance metrics.
#'
#' This implements the stacked regression approach of SR-TWAS (Dai et al.,
#' Nature Communications, 2024, \doi{10.1038/s41467-024-50983-w}). The ensemble
#' provides a principled way to combine predictions from many TWAS weight
#' methods without requiring the user to pick one method a priori or pay a
#' multiple-testing penalty for running several.
#'
#' For single-dataset usage, pass one \code{twasWeightsCv()} result directly.
#' For multi-dataset ensemble (e.g., combining cell types or reference panels
#' such as CUMC1 + MIT), pass a list of \code{twasWeightsCv()} results along
#' with a list of observed Y vectors - this learns a single joint set of
#' coefficients.
#'
#' @param cvResults Output of \code{\link{twasWeightsCv}}, with
#'   \code{$prediction} (named list of method -> out-of-fold prediction matrix,
#'   keys like \code{"susie_predicted"}). For multi-dataset: a list of such
#'   objects.
#' @param Y Observed outcome vector or matrix (samples x contexts). For
#'   multi-dataset: a list of vectors/matrices, one per dataset.
#' @param twasWeightList Optional named list of weight matrices from
#'   \code{\link{learnTwasWeights}}, with keys like \code{"susie_weights"}. Used
#'   to construct the final combined TWAS weight vector. For multi-dataset: a
#'   list of such lists (the first is used as the weight template).
#' @param contextIndex Integer indicating which column of Y to use when Y is a
#'   matrix. Default is 1 (univariate).
#' @param solver Character string specifying the optimization backend. One of
#'   \code{"quadprog"} (default), \code{"nnls"}, \code{"lbfgsb"}, or
#'   \code{"glmnet"}. \code{"quadprog"} solves a constrained QP with sum-to-1
#'   and non-negativity constraints. \code{"nnls"} uses non-negative least
#'   squares (Lawson-Hanson algorithm, as in SuperLearner) and normalizes
#'   post-hoc. \code{"lbfgsb"} uses \code{optim(method = "L-BFGS-B")} with
#'   non-negativity bounds and normalizes post-hoc. \code{"glmnet"} uses
#'   \code{cv.glmnet} with \code{lower.limits = 0} for penalized non-negative
#'   regression, providing automatic method selection via regularization. All
#'   solvers fall back to equal weights on failure.
#' @param alpha Elastic net mixing parameter, used only when \code{solver =
#'   "glmnet"}. \code{alpha = 1} (default) is lasso (sparse method selection),
#'   \code{alpha = 0} is ridge, and intermediate values give elastic net.
#'
#' @return A list with components:
#' \describe{
#'   \item{methodCoef}{Named numeric vector of combination coefficients
#'     (\eqn{\zeta_k}), non-negative and summing to 1. Names are method
#'     base names (e.g., \code{"susie"}, \code{"enet"}).}
#'   \item{ensembleTwasWeights}{Final combined weight vector
#'     \eqn{w = \sum_k \zeta_k w_k}, or NULL if \code{twasWeightList}
#'     is not provided. Returned as a vector for univariate Y, matrix
#'     otherwise.}
#'   \item{methodPerformance}{Named numeric vector of per-method R-squared
#'     computed from out-of-fold CV predictions. Preserved so users can still
#'     report individual method performance.}
#' }
#'
#' @details
#' The stacked regression solves:
#' \deqn{\min_{\zeta} \|y - P\zeta\|^2 \quad \text{s.t.} \quad
#'   \zeta_k \geq 0,\ \sum_k \zeta_k = 1}
#' where P is the \eqn{n \times K} matrix of out-of-fold predictions from K
#' methods. Four solver backends are available: \code{"quadprog"} enforces
#' both constraints during optimization; \code{"nnls"}, \code{"lbfgsb"}, and
#' \code{"glmnet"} enforce non-negativity only, then normalize coefficients
#' to sum to 1. The \code{"glmnet"} solver additionally applies
#' regularization, which can produce sparse solutions (method selection).
#' If any solver fails, the function falls back to equal weights with a
#' warning.
#'
#' Methods whose CV predictions have zero variance (e.g., when all weights are
#' zero) are excluded from the optimization and assigned \eqn{\zeta_k = 0}.
#'
#' Predictions and Y are aligned by sample names (rownames) when available,
#' rather than assuming positional order.
#'
#' @seealso \code{\link{twasWeightsCv}}, \code{\link{learnTwasWeights}},
#'   \code{\link{twasWeightsPipeline}}
#'
#' @examples
#' data(multiTraitData)
#' X <- multiTraitData$X[, 1:30]
#' y <- matrix(multiTraitData$Y[, 1], ncol = 1,
#'   dimnames = list(rownames(X), "outcome_1"))
#' # lasso/enet on this small toy panel only capture signal for some CV
#' # splits; a fixed seed keeps the example deterministic.
#' set.seed(42)
#' cv <- twasWeightsCv(X, y, fold = 3,
#'   weightMethods = list(lasso_weights = list(), enet_weights = list()))
#' ens <- ensembleWeights(cvResults = cv, Y = y)
#' ens$methodCoef # combination weights, sum to 1
#'
#' @importFrom stats optim coef complete.cases sd cor
#' @export
ensembleWeights <- function(
    cvResults,
    Y,
    twasWeightList = NULL,
    contextIndex = 1,
    solver = c("quadprog", "nnls", "lbfgsb", "glmnet"),
    alpha = 1
) {
    solver <- arg_match(solver)
    .ensembleValidateArgs(cvResults, Y, contextIndex)
    norm <- .ensembleNormalizeInput(cvResults, Y, twasWeightList)
    nm <- .ensembleMethodNames(norm$cvResults)
    stacked <- .ensembleStackPredictions(
        norm$cvResults,
        norm$Y,
        nm,
        contextIndex
    )
    methodSds <- apply(stacked$P, 2, sd)
    zeta <- .ensembleSolveZeta(
        stacked$P,
        stacked$yObs,
        methodSds,
        nm,
        solver,
        alpha
    )
    list(
        methodCoef = zeta,
        ensembleTwasWeights = .ensembleCombineWeights(
            norm$twasWeightList,
            nm$baseNames,
            zeta,
            nm$K
        ),
        methodPerformance = .ensembleMethodRsq(
            stacked$P,
            stacked$yObs,
            methodSds,
            nm$baseNames,
            nm$K
        )
    )
}

# Validate the required scalar / presence constraints on the raw inputs.
# @noRd
#' @importFrom checkmate assertCount
.ensembleValidateArgs <- function(cvResults, Y, contextIndex) {
    if (is.null(cvResults)) {
        abort("'cvResults' is required.")
    }
    if (is.null(Y)) {
        abort("'Y' is required.")
    }
    assertCount(contextIndex, positive = TRUE)
    invisible(NULL)
}

# Normalize single vs multi-dataset input to parallel lists. Single dataset:
# cvResults has $prediction directly (a twasWeightsCv() output). Multi-dataset:
# cvResults is a list of such outputs.
# @noRd
.ensembleNormalizeInput <- function(cvResults, Y, twasWeightList) {
    if (!is.null(cvResults$prediction)) {
        return(list(
            cvResults = list(cvResults),
            Y = list(Y),
            twasWeightList = if (is.null(twasWeightList)) {
                NULL
            } else {
                list(twasWeightList)
            }
        ))
    }
    .ensembleValidateMultiInput(cvResults, Y, twasWeightList)
    list(cvResults = cvResults, Y = Y, twasWeightList = twasWeightList)
}

# List-consistency checks for the multi-dataset ensemble path.
# @noRd
.ensembleValidateMultiInput <- function(cvResults, Y, twasWeightList) {
    if (!is.list(cvResults) || length(cvResults) == 0) {
        msg <- glue(
            "For multi-dataset ensemble, 'cvResults' must be a non-empty ",
            "list of twasWeightsCv() outputs."
        )
        abort(msg)
    }
    if (!is.list(Y) || length(Y) != length(cvResults)) {
        msg <- glue(
            "'Y' must be a list of the same length as 'cvResults' for ",
            "multi-dataset ensemble."
        )
        abort(msg)
    }
    if (
        !is.null(twasWeightList) &&
            (!is.list(twasWeightList) ||
                length(twasWeightList) != length(cvResults))
    ) {
        msg <- glue(
            "'twasWeightList' must be a list of the same length as ",
            "'cvResults'."
        )
        abort(msg)
    }
    for (d in seq_along(cvResults)) {
        if (is.null(cvResults[[d]]$prediction)) {
            msg <- glue(
                "cvResults[[{d}]] does not contain '$prediction'. ",
                "Expected a twasWeightsCv() output."
            )
            abort(msg)
        }
    }
    invisible(NULL)
}

# Extract + validate the method names shared across datasets. Returns
# list(predNames, baseNames, K).
# @noRd
.ensembleMethodNames <- function(cvResults) {
    predNames <- names(cvResults[[1]]$prediction)
    if (is.null(predNames) || any(predNames == "")) {
        msg <- glue(
            "cvResults$prediction must be a named list (output of ",
            "twasWeightsCv)."
        )
        abort(msg)
    }
    baseNames <- str_remove(predNames, "(_predicted|Predicted)$")
    K <- length(baseNames)
    if (K < 2) {
        msg <- glue(
            "Ensemble learning requires at least 2 methods. Found: {K}."
        )
        abort(msg)
    }
    for (d in seq_along(cvResults)) {
        if (!identical(names(cvResults[[d]]$prediction), predNames)) {
            pred1 <- str_flatten(predNames, ", ")
            predD <- str_flatten(names(cvResults[[d]]$prediction), ", ")
            msg <- glue(
                "All cvResults must have the same method names (in ",
                "$prediction) in the same order. Dataset 1 has: {pred1}; ",
                "dataset {d} has: {predD}"
            )
            abort(msg)
        }
    }
    list(predNames = predNames, baseNames = baseNames, K = K)
}

# Build the stacked prediction matrix P and observed y vector, dropping rows
# with any NA. Returns list(P, yObs).
# @noRd
.ensembleStackPredictions <- function(cvResults, Y, nm, contextIndex) {
    perDataset <- map(
        seq_along(cvResults),
        .ensembleDatasetRow,
        cvResults = cvResults,
        Y = Y,
        nm = nm,
        contextIndex = contextIndex
    )
    pMats <- map(perDataset, "P")
    P <- exec(rbind, !!!pMats)
    yObs <- list_c(map(perDataset, "y"))
    .ensembleDropIncomplete(P, yObs, nm$K)
}

# Per-dataset prediction matrix + aligned outcome. Returns list(P, y).
# @noRd
.ensembleDatasetMatrix <- function(predsD, yRaw, nm, contextIndex, d) {
    aln <- .ensembleAlignSamples(predsD, yRaw, nm$predNames, contextIndex, d)
    list(P = .ensembleBuildPd(predsD, nm, aln, contextIndex, d), y = aln$yD)
}

# Align prediction rows to the outcome, by sample name when available else
# positionally. Returns list(yD, predOrder, nD).
# @noRd
.ensembleAlignSamples <- function(predsD, yRaw, predNames, contextIndex, d) {
    predSamples <- rownames(predsD[[predNames[1]]])
    yNames <- if (is.matrix(yRaw) || is.data.frame(yRaw)) {
        rownames(yRaw)
    } else {
        names(yRaw)
    }
    if (!is.null(predSamples) && !is.null(yNames)) {
        return(.ensembleAlignByName(predSamples, yNames, yRaw, contextIndex, d))
    }
    .ensembleAlignPositional(yRaw, contextIndex, d)
}

# Name-based alignment over the intersection of sample names.
# @noRd
.ensembleAlignByName <- function(predSamples, yNames, yRaw, contextIndex, d) {
    common <- intersect(predSamples, yNames)
    if (length(common) == 0) {
        msg <- glue(
            "No common sample names between predictions and Y in ",
            "dataset {d}."
        )
        abort(msg)
    }
    if (
        length(common) < length(predSamples) ||
            length(common) < length(yNames)
    ) {
        nCommon <- length(common)
        nPred <- length(predSamples)
        nY <- length(yNames)
        msg <- glue(
            "Dataset {d}: using {nCommon} common samples ",
            "(predictions: {nPred}, Y: {nY})."
        )
        inform(msg)
    }
    yD <- if (is.matrix(yRaw) || is.data.frame(yRaw)) {
        .ensembleCheckContextIndex(contextIndex, ncol(yRaw), d)
        as.numeric(as.matrix(yRaw)[match(common, yNames), contextIndex])
    } else {
        as.numeric(yRaw[match(common, yNames)])
    }
    list(yD = yD, predOrder = match(common, predSamples), nD = length(common))
}

# Positional alignment fallback (no sample names on either side).
# @noRd
.ensembleAlignPositional <- function(yRaw, contextIndex, d) {
    yD <- if (is.matrix(yRaw) || is.data.frame(yRaw)) {
        .ensembleCheckContextIndex(contextIndex, ncol(yRaw), d)
        as.numeric(as.matrix(yRaw)[, contextIndex])
    } else {
        as.numeric(yRaw)
    }
    list(yD = yD, predOrder = seq_len(length(yD)), nD = length(yD))
}

# Guard: contextIndex must not exceed the outcome's column count.
# @noRd
.ensembleCheckContextIndex <- function(contextIndex, ncolY, d) {
    if (contextIndex > ncolY) {
        msg <- glue(
            "contextIndex ({contextIndex}) exceeds number of columns in ",
            "Y[[{d}]] ({ncolY})."
        )
        abort(msg)
    }
    invisible(NULL)
}

# One method's aligned prediction column for dataset `d`.
# @noRd
.ensemblePredColumn <- function(k, predsD, nm, aln, contextIndex, d) {
    methodName <- nm$predNames[k]
    predMat <- predsD[[methodName]]
    pCol <- if (is.matrix(predMat)) {
        predMat[aln$predOrder, contextIndex]
    } else {
        as.numeric(predMat)[aln$predOrder]
    }
    if (length(pCol) != aln$nD) {
        nCol <- length(pCol)
        nAligned <- aln$nD
        msg <- glue(
            "Prediction length for method '{methodName}' in dataset ",
            "{d} ({nCol}) does not match number of aligned samples ",
            "({nAligned})."
        )
        abort(msg)
    }
    pCol
}

# Assemble one dataset's (samples x methods) prediction matrix.
# @noRd
.ensembleBuildPd <- function(predsD, nm, aln, contextIndex, d) {
    cols <- map(
        seq_along(nm$predNames),
        .ensemblePredColumn,
        predsD = predsD,
        nm = nm,
        aln = aln,
        contextIndex = contextIndex,
        d = d
    )
    matrix(
        unname(list_c(cols)),
        nrow = aln$nD,
        ncol = nm$K,
        dimnames = list(NULL, nm$baseNames)
    )
}

# Drop rows with any NA prediction/outcome; error when too few remain.
# @noRd
.ensembleDropIncomplete <- function(P, yObs, K) {
    complete <- complete.cases(P, yObs)
    nDropped <- sum(!complete)
    if (nDropped > 0) {
        msg <- glue(
            "Dropping {nDropped} observation(s) with NA predictions or ",
            "outcomes."
        )
        inform(msg)
    }
    if (sum(complete) < K + 1) {
        nComplete <- sum(complete)
        nNeed <- K + 1
        msg <- glue(
            "Too few complete observations ({nComplete}) for {K} methods. ",
            "Need at least {nNeed}."
        )
        abort(msg)
    }
    list(P = P[complete, , drop = FALSE], yObs = yObs[complete])
}

# Solve for the method-combination coefficients zeta (length K), routing through
# the requested solver over the non-degenerate methods.
# @noRd
.ensembleSolveZeta <- function(P, yObs, methodSds, nm, solver, alpha) {
    validMethods <- methodSds > .Machine$double.eps
    nValid <- sum(validMethods)
    if (nValid < 1) {
        msg <- glue(
            "All methods have zero-variance predictions. Cannot compute ",
            "ensemble. This typically means all methods returned zero ",
            "weights - check that the input data has sufficient signal."
        )
        abort(msg)
    }
    if (nValid == 1) {
        return(.ensembleSingleMethodZeta(validMethods, nm$baseNames, nm$K))
    }
    zetaValid <- .ensembleSolveValid(
        P[, validMethods, drop = FALSE],
        yObs,
        solver,
        alpha
    )
    zeta <- replace(rep(0, nm$K), validMethods, zetaValid)
    set_names(zeta, nm$baseNames)
}

# Degenerate case: a single signal-bearing method takes full weight.
# @noRd
.ensembleSingleMethodZeta <- function(validMethods, baseNames, K) {
    zeta <- set_names(replace(rep(0, K), validMethods, 1), baseNames)
    methodName <- baseNames[validMethods]
    msg <- glue(
        "Only one method ('{methodName}') has non-zero variance ",
        "predictions. Assigning it full weight."
    )
    inform(msg)
    zeta
}

# Dispatch the valid-method coefficient solve to the chosen solver.
# @noRd
.ensembleSolveValid <- function(Pvalid, yObs, solver, alpha) {
    Kvalid <- ncol(Pvalid)
    switch(
        solver,
        quadprog = .solveEnsembleQuadprog(Pvalid, yObs, Kvalid),
        nnls = .solveEnsembleNnls(Pvalid, yObs, Kvalid),
        lbfgsb = .solveEnsembleLbfgsb(Pvalid, yObs, Kvalid),
        glmnet = .solveEnsembleGlmnet(Pvalid, yObs, Kvalid, alpha = alpha)
    )
}

# Per-method out-of-sample R^2 (NA for zero-variance methods).
# @noRd
.ensembleMethodRsq <- function(P, yObs, methodSds, baseNames, K) {
    set_names(
        map_dbl(
            seq_len(K),
            .ensembleMethodR2,
            methodSds = methodSds,
            yObs = yObs,
            P = P
        ),
        baseNames
    )
}

# Combine the per-method TWAS weight matrices (from the first dataset) using the
# fitted coefficients zeta. Returns NULL when no weights are supplied/matched.
# @noRd
.ensembleCombineWeights <- function(twasWeightList, baseNames, zeta, K) {
    if (is.null(twasWeightList)) {
        return(NULL)
    }
    wtList <- twasWeightList[[1]]
    if (!is.list(wtList) || length(wtList) == 0) {
        msg <- glue(
            "twasWeightList[[1]] is empty or not a list; skipping weight ",
            "combination."
        )
        warn(msg)
        return(NULL)
    }
    wtKeys <- str_c(baseNames, "_weights")
    matched <- is_in(wtKeys, names(wtList))
    if (!any(matched)) {
        keyExamples <- str_flatten(wtKeys[seq_len(min(3, K))], ", ")
        msg <- glue(
            "No matching weight keys found in twasWeightList. Expected keys ",
            "like: {keyExamples}"
        )
        warn(msg)
        return(NULL)
    }
    .ensembleAccumulateWeights(wtList, wtKeys, matched, zeta)
}

# Coerce a weight vector/matrix to a (variants x contexts) matrix.
# @noRd
.ensembleAsMatrix <- function(w) {
    if (is.matrix(w)) w else matrix(w, ncol = 1)
}

# One method's zeta-scaled contribution, or NULL when its weight matrix does
# not line up with the first one's shape.
# @noRd
.ensembleWeightTerm <- function(i, wtList, wtKeys, zeta, shape) {
    wMat <- .ensembleAsMatrix(wtList[[wtKeys[i]]])
    if (!identical(dim(wMat), shape)) {
        wtKey <- wtKeys[i]
        msg <- glue(
            "Weight matrix for '{wtKey}' has inconsistent dimensions; ",
            "skipping."
        )
        warn(msg)
        return(NULL)
    }
    zeta[i] * wMat
}

# Zeta-weighted sum of the matched weight matrices; univariate -> named vector.
# @noRd
.ensembleAccumulateWeights <- function(wtList, wtKeys, matched, zeta) {
    firstWt <- .ensembleAsMatrix(wtList[[wtKeys[which(matched)[1]]]])
    shape <- dim(firstWt)
    # The ensemble is the sum of the scaled contributions, so it is a fold
    # over them rather than a matrix added into repeatedly.
    ensembleTwasWt <- reduce(
        compact(map(
            which(matched),
            .ensembleWeightTerm,
            wtList = wtList,
            wtKeys = wtKeys,
            zeta = zeta,
            shape = shape
        )),
        `+`,
        .init = matrix(
            0,
            nrow = nrow(firstWt),
            ncol = ncol(firstWt),
            dimnames = dimnames(firstWt)
        )
    )
    # For the univariate case, return as a named vector.
    if (ncol(ensembleTwasWt) == 1) {
        return(set_names(as.numeric(ensembleTwasWt), rownames(ensembleTwasWt)))
    }
    ensembleTwasWt
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# One method's per-region fit for region `i`: pass region-list fits through by
# position, drop non-region-list fits (they can't be aligned across blocks).
# @noRd
.twasRegionFitOf <- function(f, i) {
    isRegionList <- is.list(f) &&
        !is.null(names(f)) &&
        length(f) > 0L &&
        all(str_length(names(f)) > 0L) &&
        all(str_starts(names(f), "region"))
    if (isRegionList && i <= length(f)) f[[i]] else NULL
}

# One per-region cvResult reporting row (region label + metric columns), or NULL
# when the entry carries no CV metrics.
# @noRd
.twasCvRow <- function(e, lab) {
    cv <- .rowCvResult(e)
    if (is.null(cv) || is.null(cv$metrics)) {
        return(NULL)
    }
    bind_cols(
        tibble(region = lab),
        as_tibble(as.list(cv$metrics), .name_repair = "minimal")
    )
}

# Canonical method name for a snake_case token (identity when not in the table).
# @noRd
.twasCanonicalMethod <- function(s, snakeToCanonical) {
    if (!is.na(snakeToCanonical[s])) snakeToCanonical[[s]] else s
}

# Capability violation (if any) for one method token against the input kind.
# @noRd
.twasTokenViolation <- function(tk, caps, inputKind) {
    .twasCapabilityViolation(caps[[tk]], tk, inputKind)
}

# One region block of .twasRunJointSpecGrid().
# @noRd
.twasJointSpecRegion <- function(
    bi,
    spec,
    tokens,
    contexts,
    traits,
    marker,
    ctx
) {
    .runJointSpecs(
        spec,
        ctx$data,
        dataForm = "individual",
        pipeline = marker,
        jointMethods = tokens,
        contexts = contexts,
        traitIds = traits,
        args = .twasJointArgs(ctx, bi)
    )
}

# The engine `args` of an individual-level TWAS joint run over region block
# `bi`. Residualization rides here: the group enumerators read it while
# extracting X and Y.
# @noRd
.twasJointArgs <- function(ctx, bi) {
    list(
        methodList = ctx$norm$methodList,
        fineMappingResult = ctx$fineMappingResult,
        dataDrivenPriorMatricesCv = ctx$dataDrivenPriorMatricesCv,
        cisWindow = ctx$cisWindow,
        region = ctx$xRegions[[bi]],
        regionIndex = bi,
        nRegions = length(ctx$xRegions),
        residualizationParam = ctx$residualizationParam,
        genotypeFilter = ctx$genotypeFilter,
        phenotypeFilter = ctx$phenotypeFilter,
        verbose = ctx$verbose
    )
}

# One region's univariate joint-cell fit (region `bi` of `xRegions`).
# @noRd
.twasQdsUnivRegion <- function(
    bi,
    cells,
    scope,
    marker,
    data,
    xRegions,
    norm,
    fineMappingResult,
    twasWeights,
    dataDrivenPriorMatricesCv,
    cisWindow,
    phenotypeFilter,
    residualizationParam,
    genotypeFilter,
    nPCs,
    verbose
) {
    # `cells` is the univariate cell, plus the top-PC cell when usePCA. Each
    # enumerates its own groups from the same data and returns its own rows.
    args <- .twasQdsUnivArgs(
        bi,
        xRegions = xRegions,
        norm = norm,
        fineMappingResult = fineMappingResult,
        twasWeights = twasWeights,
        dataDrivenPriorMatricesCv = dataDrivenPriorMatricesCv,
        cisWindow = cisWindow,
        phenotypeFilter = phenotypeFilter,
        residualizationParam = residualizationParam,
        genotypeFilter = genotypeFilter,
        nPCs = nPCs,
        verbose = verbose
    )
    out <- compact(map(
        cells,
        .runJointCell,
        pipeline = marker,
        data = data,
        scope = scope,
        tokens = norm$tokens,
        args = args
    ))
    .twasCombineCellRows(out)
}

# Combine the per-cell results of one region. rbind, NOT
# .twasMergeResultsByKey(): that one merges the SAME rows across regions,
# keeping results[[1]]'s row set. The top-PC rows are different rows, so
# merging would silently drop them.
# @noRd
.twasCombineCellRows <- function(out) {
    if (length(out) == 0L) {
        return(NULL)
    }
    if (length(out) == 1L) {
        return(out[[1L]])
    }
    reduce(out, .rbindTwasWeights)
}

# One cached row record from an imap over (entry, token) cache hits.
# @noRd
.twasCachedRowRecord <- function(entry, tk, st, ctx, tr) {
    .twasRowRecord(st, ctx, tr, tk, entry)
}

# A resumed multivariate group as row records. `.twasMvCacheHits` hands back
# bare TwasWeightsRow entries keyed by context, but the assembler reads the
# tuple axes off each row, so they are attached here exactly as the fitted
# path attaches them -- returning the bare entries leaves `study` absent.
# @noRd
.twasMvCachedRows <- function(cached, st, tr, tk) {
    unname(imap(cached, .twasMvCachedRow, st = st, tr = tr, tk = tk))
}

# One resumed row from an imap over (entry, context) cache hits.
# @noRd
.twasMvCachedRow <- function(entry, ctx, st, tr, tk) {
    .twasRowRecord(st, ctx, tr, tk, entry)
}

# Cache lookup for one token in a (study, context, trait).
# @noRd
.twasCacheLookupTok <- function(tk, twasWeights, st, ctx, tr) {
    .twasCacheLookup(twasWeights, st, ctx, tr, tk)
}

# One context row (`kk`) of a fitted multivariate weight matrix. The shared
# joint
# fit rides on the first row only; later rows leave fits NULL.
# @noRd
.twasMvContextRow <- function(
    kk,
    weights,
    fitAttr,
    ctxNames,
    mvStat,
    st,
    tr,
    tk,
    weightsDataType
) {
    .twasRowRecord(
        st,
        ctxNames[[kk]],
        tr,
        tk,
        TwasWeightsRow(
            variantIds = mvStat$variantIds,
            weights = as.numeric(weights[, kk]),
            methodFits = if (kk == 1L) fitAttr else NULL,
            cvResult = NULL,
            weightStandardized = TRUE,
            weightsDataType = weightsDataType
        )
    )
}

# The stacked prediction/observed matrix for ensemble dataset `d`.
# @noRd
.ensembleDatasetRow <- function(d, cvResults, Y, nm, contextIndex) {
    .ensembleDatasetMatrix(
        cvResults[[d]]$prediction,
        Y[[d]],
        nm,
        contextIndex,
        d
    )
}

# Out-of-sample R^2 for ensemble method `k` (NA for a zero-variance method).
# @noRd
.ensembleMethodR2 <- function(k, methodSds, yObs, P) {
    if (methodSds[k] > 0) cor(yObs, P[, k])^2 else NA_real_
}

# TRUE when a per-region cvResult element carries a sample partition.
# @noRd
.twasCvHasPartition <- function(z) {
    is.list(z) && !is.null(z$samplePartition)
}

#' @rdname TwasWeightsMethodsParam
#' @aliases TwasWeightsMethodsParam-class
#' @exportClass TwasWeightsMethodsParam
setClass("TwasWeightsMethodsParam", contains = "MethodsSelectionParam")

#' @title Which TWAS Weight Methods To Run, And How
#' @description Selects the methods \code{\link{twasWeightsPipeline}} runs and
#'   carries each one's engine arguments.
#'
#'   A method that runs on both individual-level and summary-statistic data
#'   reaches a different engine on each path, so for a
#'   \code{MultiStudyQtlDataset} carrying both there is no single set of
#'   arguments per method. Name such a method under
#'   \code{qtlDatasetMethods} or \code{qtlSumStatsMethods} to say which
#'   path its options are for; name it under \code{methods} when the path
#'   need not be stated --- a single-type run, or a method with nothing
#'   path-specific to configure.
#'
#'   Naming a method selects it. A method belongs in exactly one of the
#'   three slots.
#' @param methods Named list of per-method options whose input path need not
#'   be stated. Each entry is that method's \code{*Options()} record, or
#'   \code{list()} to run it with its defaults.
#' @param qtlDatasetMethods Named list of per-method options for the
#'   individual-level path.
#' @param qtlSumStatsMethods Named list of per-method options for the
#'   summary-statistics path.
#' @return A \code{TwasWeightsMethodsParam} object, a \code{\link{MethodParam}}.
#' @examples
#' TwasWeightsMethodsParam(methods = list(lasso = list(), susie = list()))
#' TwasWeightsMethodsParam(
#'     methods = list(susie = list()),
#'     qtlDatasetMethods = list(lasso = GlmnetOptions(nfolds = 5)),
#'     qtlSumStatsMethods = list(lasso = LassosumOptions())
#' )
#' @export
TwasWeightsMethodsParam <- function(
    methods = NULL,
    qtlDatasetMethods = NULL,
    qtlSumStatsMethods = NULL
) {
    new(
        "TwasWeightsMethodsParam",
        methods = .methodsNormalizeSlot(
            methods,
            "methods",
            "TwasWeightsMethodsParam"
        ),
        qtlDatasetMethods = .methodsNormalizeSlot(
            qtlDatasetMethods,
            "qtlDatasetMethods",
            "TwasWeightsMethodsParam"
        ),
        qtlSumStatsMethods = .methodsNormalizeSlot(
            qtlSumStatsMethods,
            "qtlSumStatsMethods",
            "TwasWeightsMethodsParam"
        )
    )
}

#' @rdname EnsembleParam
#' @aliases EnsembleParam-class
#' @exportClass EnsembleParam
setClass(
    "EnsembleParam",
    contains = "MethodParam",
    slots = c(
        enabled = "logical",
        r2Threshold = "numeric",
        solver = "character",
        alpha = "numeric"
    )
)

#' @title SR-TWAS Ensemble Settings
#' @description Whether and how \code{\link{twasWeightsPipeline}} stacks its
#'   per-method weights into an SR-TWAS ensemble.
#' @section Requires cross-validation:
#'   Stacking combines each method's \strong{out-of-fold} predictions, so it
#'   cannot run without cross-validation. \code{enabled = TRUE} together with
#'   \code{CrossValidationParam(folds < 2)} is an error. It used to be
#'   neither: the ensemble row was simply absent from the result, with
#'   \code{ensemble = TRUE} still reading as on.
#'
#'   \code{enabled} defaults to \code{FALSE} because
#'   \code{\link{CrossValidationParam}} defaults to no folds. To get an
#'   ensemble, ask for both.
#' @param enabled Logical. Compute SR-TWAS ensemble weights. Default
#'   \code{FALSE}.
#' @param r2Threshold Minimum cross-validated \eqn{R^2} for a method to enter
#'   the stack. Default \code{0.01}. Stacking needs at least two methods to
#'   clear it.
#' @param solver Stacking solver, \code{"quadprog"} (default) or
#'   \code{"glmnet"}.
#' @param alpha Elastic-net mixing parameter, used only when
#'   \code{solver = "glmnet"}. Default \code{1}.
#' @return A \code{EnsembleParam} object, a \code{\link{MethodParam}}.
#' @seealso \code{\link{CrossValidationParam}}
#' @examples
#' EnsembleParam(enabled = TRUE, r2Threshold = 0.05)
#' @export
EnsembleParam <- function(
    enabled = FALSE,
    r2Threshold = 0.01,
    solver = c("quadprog", "glmnet"),
    alpha = 1
) {
    solver <- arg_match(solver)
    new(
        "EnsembleParam",
        enabled = enabled,
        r2Threshold = r2Threshold,
        solver = solver,
        alpha = alpha
    )
}

# Refuse an ensemble that cannot be built. Stacking reads out-of-fold
# predictions, so without folds there is nothing to stack -- and the old
# behaviour was to return a result silently missing its ensemble row.
# @noRd
.ensembleAssertCv <- function(ensembleParam, crossValidationParam) {
    if (!isTRUE(ensembleParam$enabled) || .cvEnabled(crossValidationParam)) {
        return(invisible(NULL))
    }
    abort(glue(
        "twasWeightsPipeline: EnsembleParam(enabled = TRUE) needs ",
        "out-of-fold predictions to stack, so it requires ",
        "CrossValidationParam(folds >= 2); got folds = ",
        "{crossValidationParam$folds %||% 0}."
    ))
}

# The ensemble settings as a plain list, so a cfg record can carry the group
# as one named field instead of four loose ones.
# @noRd
.ensembleResolve <- function(ensembleParam) {
    list(
        enabled = isTRUE(ensembleParam$enabled),
        r2Threshold = ensembleParam$r2Threshold %||% 0.01,
        solver = ensembleParam$solver %||% "quadprog",
        alpha = ensembleParam$alpha %||% 1
    )
}
