# =============================================================================
# Helpers + S4 dispatch surface for twasWeightsPipeline
# =============================================================================

# Concatenate two TwasWeights collections row-wise, carrying forward every
# column (delegates to the generic `.rbindCollections`, which unions columns
# and pads a side lacking an optional column such as joint* / region).
# @noRd
#' @importFrom checkmate assertClass
.rbindTwasWeights <- function(a, b, ldSketch = NULL) {
    assertClass(a, "TwasWeights")
    assertClass(b, "TwasWeights")
    # Carry forward every column (joint*, region, ...) and reconcile the
    # collection-level slots via the shared combine.
    .combineTupleCollections(list(a, b), ldSketch, ".rbindTwasWeights")
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
#' @param ldSketch Optional genotype panel (see \code{\link{readGenotypes}}) to
#'   attach to the
#'   combined collection. Default \code{NULL}. Applied when combining two or
#'   more inputs; a single input is returned unchanged.
#' @return A single combined \code{TwasWeights}.
#' @seealso \code{\link{combineFineMappingResults}}
#' @examples
#' twe <- twasWeightsRow(
#'   variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)), weights = rep(0.1, 4),
#'   cvResult = list(rsq = 0.5), standardized = FALSE)
#' tw1 <- TwasWeights(study = "s1", context = "brain", trait = "g1",
#'   method = "susie", entry = list(twe))
#' tw2 <- TwasWeights(study = "s2", context = "brain", trait = "g1",
#'   method = "susie", entry = list(twe))
#' combineTwasWeights(tw1, tw2)
#' @importFrom stringr str_starts
#' @export
combineTwasWeights <- function(..., ldSketch = NULL) {
    parts <- .asCombineList(list(...), "TwasWeights", "combineTwasWeights")
    .combineTupleCollections(parts, ldSketch, "combineTwasWeights")
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
    wList <- map(payloads, getWeights)
    weights <- if (is.matrix(wList[[1L]])) {
        exec(rbind, !!!wList)
    } else {
        unname(list_c(wList))
    }
    twasWeightsRow(
        variantIds = unname(list_c(map(payloads, .twrPartsVariantIds))),
        weights = weights,
        fits = set_names(map(payloads, getFits), fittedLabels),
        cvResult = .twasRegionCvDf(payloads, fittedLabels),
        standardized = getStandardized(payloads[[1L]]),
        dataType = getDataType(payloads[[1L]])
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
    cvFits <- getCvFits(mashPrior)
    perFold <- if (!is.null(cvFits)) cvFits$perFoldFits else NULL
    sp <- samplePartition %||% cvFits$samplePartition
    list(
        fullPrior = getFullFit(mashPrior),
        dataDrivenPriorMatricesCv = perFold,
        samplePartition = sp
    )
}

# Mapping from short / canonical TWAS weight-method name to dispatch
# capability. Used to reject incompatible (input class, method) pairs.
#
# `allowsIndiv`  : may be invoked on a QtlDataset (individual-level X, Y).
# `allowsRss`    : may be invoked on a QtlSumStats / GwasSumStats (RSS).
# `multivariate` : requires a multi-trait / multi-context Y (mvsusie /
#                  mr.mash family).
#
# Rules from `dev/refactor-design.md` (`twasWeightsPipeline` row):
# - PRS-CS is RSS-only.
# - BGLR / CRAN-stable qgg methods (bayes_a/b/c/l/n/r, b_lasso, dpr_*)
#   are individual-level only.
# - mr.mash / mvsusie follow the multi-trait / multi-context rules of
#   the mvSuSiE fine-mapping family.
# @noRd
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
# dpr_vb / dpr_adaptive_gibbs remain individual-only. enet has no cpp11
# sumstat solver yet (lassosumRssRcpp is pure L1, no alpha mixing) and is
# documented as individual-only for now. prsCs has no individual-level
# counterpart (it is a sumstat-only Bayesian shrinkage method).
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
    # Individual-only -- no cpp11 sumstat solver yet.
    enet = list(
        individualImpl = "enetWeights",
        sumstatImpl = NULL,
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
    if (.isMethodOptions(methods)) {
        # The retired TwasWeightsMethodsParam() record.
        return(.twasNormalizeListMethods(map(as.list(methods), as.list)))
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

# Named-list `methods`: re-key each entry to its canonical `<token>_weights`
# name, merge the caller's kwargs over the method's defaults, and carry the
# token->impl map as an "impl" attribute (without this, downstream
# .resolveMethodFunction falls back to the bare token, which is not a function).
# @noRd
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
# FineMappingResult, matched tolerantly across canonical / camelCase /
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
# fsusie) when no FineMappingResult is supplied. twasWeightsPipeline is
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
    cvFolds = 0
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
    invisible(NULL)
}

# Whether row `i` of a FineMappingResult carries a cross-validation result.
# @noRd
.fmrRowHasCv <- function(i, fineMappingResult) {
    !is.null(getCvResult(.fmrRowParts(fineMappingResult, i)))
}

# Whether a FineMappingResult carries any cross-validation result at all.
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
# must be a FineMappingResult.
# @noRd
#' @importFrom checkmate assertClass
.twasRequireFmResult <- function(fineMappingResult, fmTokens) {
    if (is.null(fineMappingResult)) {
        fmStr <- str_flatten(unique(fmTokens), ", ")
        msg <- glue(
            "twasWeightsPipeline: method(s) {fmStr} are fine-mapping ",
            "methods and may not be re-fit by twasWeightsPipeline. Run ",
            "fineMappingPipeline() first and pass the result via ",
            "`fineMappingResult = <FineMappingResult>`."
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
        "and pass the result via `fineMappingResult = <FineMappingResult>`."
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

# Enforce the multi-trait / multi-context rule for mvsusie / mr.mash
# methods (same family as the fine-mapping mvSuSiE rule in the design
# doc). Multivariate methods need at least 2 traits *or* at least 2
# contexts in the Y matrix passed to learnTwasWeights.
# @noRd
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
    components <- getQtlDatasets(data)
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
    nTraits <- length(tryCatch(getTraits(d), error = function(cnd) {
        character(0)
    }))
    nCtx <- length(tryCatch(getContexts(d), error = function(cnd) character(0)))
    nTraits < 2L && nCtx < 2L
}

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

# Reject SumStats inputs that have not been QC'd via summaryStatsQc.
# @noRd
.twasAssertQcd <- function(sumstats) {
    if (length(getQcInfo(sumstats)) == 0L) {
        cls <- class(sumstats)[[1L]]
        msg <- glue(
            "twasWeightsPipeline: the supplied {cls} has no QC record ",
            "(qcInfo is empty). Call summaryStatsQc() first and pass the ",
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
.twasCacheLookup <- function(twasWeights, study, context, trait, method) {
    if (is.null(twasWeights)) {
        return(NULL)
    }
    if (!is(twasWeights, "TwasWeights")) {
        return(NULL)
    }
    idx <- .matchTupleRows(
        twasWeights,
        list(study = study, context = context, trait = trait, method = method)
    )
    if (length(idx) == 0L) {
        return(NULL)
    }
    .twrRowParts(twasWeights, idx[[1L]])
}

# Convert a FineMappingResult (single-method susie/susie_inf row matched
# to the requested study/context/trait) into a `fittedModels` list
# suitable for `learnTwasWeights`. Pulls the trimmedFit from the matching
# entry. Returns a (possibly empty) list.
# @noRd
#' @importFrom checkmate assertClass
.twasFineMappingFits <- function(fineMappingResult, study, context, trait) {
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
            study = study,
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
    study,
    context,
    trait
) {
    idx <- which(
        is_in(
            str_to_lower(as.character(fineMappingResult$method)),
            .twasSpellingCandidates(canonical)
        ) &
            as.character(fineMappingResult$study) == study &
            as.character(fineMappingResult$context) == context &
            as.character(fineMappingResult$trait) == trait
    )
    if (length(idx) == 0L) {
        return(NULL)
    }
    getSusieFit(.fmrRowParts(fineMappingResult, idx[[1L]]))
}

# Locate a fine-mapping fit for one (study, context, trait, token) tuple.
# Used by the QtlSumStats sumstat dispatcher to pass the precomputed fit
# into susieRssWeights / susieInfRssWeights / susieAshRssWeights /
# mvsusieRssWeights via their respective *Fit arguments.
# @noRd
.twasFineMappingFitFor <- function(
    fineMappingResult,
    study,
    context,
    trait,
    token
) {
    if (is.null(fineMappingResult)) {
        return(NULL)
    }
    fits <- .twasFineMappingFits(
        fineMappingResult,
        study = study,
        context = context,
        trait = trait
    )
    fits[[token]]
}

# Collect the cross-validation payload that fineMappingPipeline stored on the
# FineMappingResult for one (study, context, trait) tuple. fineMapping records
# one cvResult per (study, context, trait, method) entry (samplePartition +
# per-fold predictions/metrics, keyed by the TWAS snake method name); this
# merges them across the fine-mapping methods of the tuple into a single
# twasWeightsCv()-shaped list so twasWeightsPipeline can reuse the partition and
# feed those out-of-fold predictions into the SR-TWAS ensemble without re-
# fitting the fine-mapping models. A multi-region entry stores cvResult as a
# per-region list; the first region carrying CV is used. Returns NULL when no
# fine-mapping entry for the tuple recorded CV.
# @noRd
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
    cv <- getCvResult(.fmrRowParts(fineMappingResult, i))
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

.twasCvResultFor <- function(fineMappingResult, study, context, trait) {
    if (is.null(fineMappingResult)) {
        return(NULL)
    }
    if (!is(fineMappingResult, "FineMappingResultBase")) {
        return(NULL)
    }
    idx <- which(
        as.character(fineMappingResult$study) == study &
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
#'   \item a \code{\link{GwasSumStats}} for per-LD-block PRS-CS-style fits
#'         from GWAS summary statistics.
#' }
#'
#' Method-restriction rules (enforced):
#' \itemize{
#'   \item \code{mr.mash}, \code{mvsusie} follow the multi-trait /
#'         multi-context rules of the fine-mapping \code{mvsusie} family
#'         (require at least two traits OR at least two contexts).
#'   \item RSS-only methods (PRS-CS, \code{lassosumRss}, SDPR, all
#'         \code{*Rss} variants) are rejected on \code{QtlDataset}
#'         input.
#'   \item Individual-level-only methods (BGLR and CRAN-stable qgg:
#'         \code{bayes_a/b/c/l/n/r}, \code{b_lasso}, \code{dpr_*}) are
#'         rejected on \code{QtlSumStats} / \code{GwasSumStats} input.
#' }
#'
#' Both \code{QtlSumStats} and \code{GwasSumStats} inputs must have been QC'd
#' via \code{\link{summaryStatsQc}} first; otherwise an error is raised pointing
#' at that function.
#'
#' The returned \code{\link{TwasWeights}} collection's \code{ldSketch} slot is
#' set automatically: \code{NULL} for individual-level fits, the input's
#' \code{ldSketch} for RSS-derived fits.
#'
#' Optionally a \code{FineMappingResult} may be supplied as a source of pre-fit
#' SuSiE / SuSiE-inf / SuSiE-ash objects; their \code{trimmedFit} payloads are
#' passed through to \code{learnTwasWeights} / the RSS sub-pipelines via the
#' \code{fittedModels} slot, avoiding a re-fit.
#'
#' When the supplied \code{FineMappingResult} was produced with cross-validation
#' (\code{fineMappingPipeline(..., cvFolds > 1)}), each matching \code{(study,
#' context, trait)} entry's \code{cvResult} is reused: its fold partition
#' becomes the CV partition (unless \code{samplePartition} is given explicitly)
#' and its per-fold out-of-fold predictions/metrics are fed directly into the
#' SR-TWAS ensemble in place of re-fitting those fine-mapping methods here.
#' Non-fine-mapping methods (lasso, enet, ...) are still cross-validated on the
#' same shared partition.
#'
#' @param data A \code{QtlDataset}, \code{MultiStudyQtlDataset}, or
#'   \code{QtlSumStats}. The \code{MultiStudyQtlDataset} method iterates the
#'   embedded individual-level \code{QtlDataset} entries and the optional
#'   embedded \code{QtlSumStats}, then rbinds the results.
#' @param methods A character vector of short method names, a preset string
#'   (\code{"default"} or \code{"fastDefault"}), or a named list of
#'   \code{<method>_weights = args} entries. For QtlSumStats / GwasSumStats
#'   inputs the default switches to the RSS preset (\code{c("susieRss",
#'   "susieInfRss", "lassosumRss", "prsCs", "sdpr")}).
#' @param contexts Optional character vector of contexts to restrict processing
#'   to (QtlDataset / QtlSumStats inputs). Default \code{NULL} (use all
#'   contexts).
#' @param traitId Optional character vector of trait identifiers to restrict
#'   processing to (QtlDataset / QtlSumStats inputs). Default \code{NULL}.
#' @param region Optional variant window for QtlDataset trait selection: a
#'   \code{GRanges}, a \code{"chr:start-end"} string, or a one-row data.frame
#'   with \code{chrom}/\code{start}/\code{end}. Mutually exclusive with
#'   \code{traitId}.
#' @param cisWindow For QtlDataset: cis-window (bp) around each trait's genomic
#'   position when extracting variants. Required when \code{traitId} is
#'   supplied. Mutually exclusive with \code{region}.
#' @section Panel filters on the RSS path: On \code{QtlSumStats} input there
#'   is no genotype matrix to filter, so \code{panelFilter}'s
#'   \code{mafCutoff} / \code{macCutoff} / \code{imissCutoff} are measured
#'   against the \strong{LD reference panel} instead: a variant whose panel
#'   genotypes fall below the cutoffs is dropped before the z-scores and LD
#'   matrix are built. The thresholds mean the same thing as
#'   \code{genotypeFilter}'s on the \code{QtlDataset} path (MAC is converted
#'   to a MAF equivalent and the stricter of the two applies), so one number
#'   carries across input types. The defaults filter nothing.
#'
#'   This discards \emph{observed} variants, unlike
#'   \code{summaryStatsQc(imputeArgs = ...)}, which only bounds which variants
#'   RAISS will impute.
#'
#' @param genotypeFilterArgs For QtlDataset: per-call genotype-filter overrides,
#'   built with \code{\link{GenotypeFilterParam}}. Each field that is set
#'   replaces the corresponding construct-time \code{\link{QtlDataset}} slot
#'   for this call only (applied to a validated copy); a field left unset
#'   leaves the stored value in place. Variant QC is a property of the data,
#'   so these are applied identically here and in
#'   \code{\link{fineMappingPipeline}} --- there is deliberately no
#'   TWAS-specific variant filter.
#' @param panelFilterArgs For QtlSumStats: LD-reference-panel filters,
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
#' @param jointSpecification Optional joint-fit specification (NULL by default).
#'   When NULL, the pipeline runs the implicit multi-trait / multi-context
#'   mr.mash branches as before. When non-NULL, the argument is parsed and
#'   validated via the joint-spec grammar documented under
#'   \code{parseJointSpecification}; the per-spec axis dispatcher implementation
#'   is in progress and a non-NULL value currently errors with an informative
#'   message.
#' @param fineMappingResult Optional \code{FineMappingResult}. When supplied,
#'   its SuSiE / SuSiE-inf / SuSiE-ash trimmed fits for the matching (study,
#'   context, trait) tuples are injected into \code{learnTwasWeights} via
#'   \code{fittedModels} so SuSiE-family weight methods reuse the prior fit
#'   instead of refitting.
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
#' @param fitFullData Logical. If \code{TRUE}, fit final weights on the full
#'   data after cross-validation. Default \code{TRUE}.
#' @param usePCA Logical (length 1). \code{QtlDataset} only. When
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
#'   \code{"full"} is honoured only by the mr.mash engines, the only ones
#'   that distinguish the two retaining levels; every other method keeps its
#'   fit whole or not at all. Requesting \code{"full"} for a run where no
#'   method can honour it is an error, and a mixed run warns rather than
#'   dropping it in silence.
#' @param naAction Character. How to handle missing values in the extracted
#'   data.
#' @param crossValidationArgs Cross-validation settings, built with
#'   \code{\link{CrossValidationParam}}: \code{folds} (default \code{0}, no
#'   CV), \code{threads}, \code{samplePartition}, \code{maxVariants} (cap
#'   on the CV design matrix) and \code{weightMethods} (which methods to
#'   cross-validate; unset means every method with non-zero weights).
#'
#'   \code{folds} defaulted to \code{5} before, while
#'   \code{\link{fineMappingPipeline}}'s defaulted to \code{0}. The two now
#'   agree on \code{0}, so an unspecified setting means one thing.
#'
#'   Cross-validation holds out \strong{samples}, so it is not possible on
#'   \code{QtlSumStats} input; that method \strong{errors} if it is set
#'   rather than returning results that were never cross-validated.
#' @param ensembleArgs SR-TWAS ensemble settings, built with
#'   \code{\link{EnsembleParam}}: \code{enabled} (default \code{FALSE}),
#'   \code{r2Threshold}, \code{solver} and \code{alpha}. Stacking reads
#'   out-of-fold predictions, so \code{enabled = TRUE} requires
#'   \code{crossValidation} with \code{folds >= 2} and is an error
#'   otherwise --- it previously returned a result quietly missing its
#'   ensemble row.
#' @param estimatePi If TRUE, estimate spike-and-slab sparsity from mr.ash
#'   before BGLR / qgg spike-and-slab methods that consume it.
#' @param residualizationArgs Covariate residualization settings, built with
#'   \code{\link{ResidualizationParam}} and forwarded to
#'   \code{\link{getResidualizedPhenotypes}} /
#'   \code{\link{getResidualizedGenotypes}}: \code{phenotypeCovariates} and
#'   \code{genotypeCovariates} name which covariates to regress out
#'   (\code{NULL}, the default, uses every available one), and
#'   \code{residualizePhenotype} / \code{residualizeGenotype} turn each side
#'   off.
#'
#'   \code{QtlDataset} / \code{MultiStudyQtlDataset} input only. A
#'   \code{QtlSumStats} run has no covariates to regress out, so the bundle
#'   is \strong{ignored} there rather than refused --- unlike
#'   \code{crossValidation}, which errors. Residualization is on by default,
#'   so refusing a non-default value would reject the default bundle; CV is
#'   off by default, so a non-default value there is an explicit request for
#'   something the input cannot do.
#' @param dataType Optional data-type label recorded on every
#'   \code{TwasWeightsRow$dataType} (e.g. \code{"expression"}).
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
#'   context, trait, method)}. The \code{ldSketch} slot is \code{NULL} for
#'   individual-level fits and equals the input's \code{ldSketch} for
#'   RSS-derived fits.
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

# Run the multivariate joint TWAS-weight fit over the (context, trait) grid for
# `traits`: dispatch each cis-region through the joint engine and merge
# per-region results. `marker` = the TwasJointPipeline config; `ctx` bundles the
# shared state
# (xRegions, data, norm, useCtx, fineMappingResult, dataDrivenPriorMatricesCv,
# cisWindow, verbose).
# @noRd
.twasRunMultivariateGrid <- function(traits, marker, ctx) {
    synthSpec <- list(list(axes = c("context", "trait"), scope = NULL))
    allLabs <- map_chr(ctx$xRegions, .twasRegionLabel)
    allRegions <- map(
        seq_along(ctx$xRegions),
        .twasMvGridRegion,
        synthSpec = synthSpec,
        marker = marker,
        ctx = ctx,
        traits = traits
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
    residualizationArgs,
    verbose
) {
    .twasRunMultivariateGrid(
        grid$allTraits,
        marker,
        list(
            xRegions = xRegions,
            data = data,
            norm = norm,
            useCtx = grid$useCtx,
            fineMappingResult = fineMappingResult,
            dataDrivenPriorMatricesCv = dataDrivenPriorMatricesCv,
            cisWindow = cisWindow,
            residualizationArgs = residualizationArgs,
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
    naAction,
    residualizationArgs,
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
            residualizationArgs = residualizationArgs,
            verbose = verbose
        ))
    }
    .twasQdsUnivariateEngine(
        study = grid$study,
        useCtx = grid$useCtx,
        allTraits = grid$allTraits,
        marker = marker,
        data = data,
        xRegions = xRegions,
        norm = norm,
        fineMappingResult = fineMappingResult,
        twasWeights = twasWeights,
        dataDrivenPriorMatricesCv = dataDrivenPriorMatricesCv,
        cisWindow = cisWindow,
        naAction = naAction,
        residualizationArgs = residualizationArgs,
        usePCA = usePCA,
        nPCs = nPCs,
        verbose = verbose
    )
}

.twasPipelineQtlDataset <- function(
    data,
    methods = "default",
    contexts = NULL,
    traitId = NULL,
    region = NULL,
    cisWindow = NULL,
    genotypeFilterArgs = GenotypeFilterParam(),
    jointRegions = FALSE,
    jointSpecification = NULL,
    fineMappingResult = NULL,
    twasWeights = NULL,
    mashPrior = NULL,
    crossValidationArgs = CrossValidationParam(),
    fitFullData = TRUE,
    ensembleArgs = EnsembleParam(),
    estimatePi = TRUE,
    usePCA = FALSE,
    nPCs = 10L,
    fitRetention = c("slim", "none", "full"),
    residualizationArgs = ResidualizationParam(),
    dataType = NULL,
    naAction = c("drop", "impute"),
    verbose = 1,
    seed = NULL
) {
    naAction <- arg_match(naAction)
    fitRetention <- arg_match(fitRetention)
    .twasWarnUnretainableDetail(fitRetention, methods)
    .assertMethodParam(
        crossValidationArgs,
        "CrossValidationParam",
        "crossValidation"
    )
    .assertMethodParam(ensembleArgs, "EnsembleParam", "ensemble")
    .ensembleAssertCv(ensembleArgs, crossValidationArgs)
    cvCfg <- .cvResolve(crossValidationArgs)
    # Each stage returns only the values it derives; nothing is grafted onto a
    # captured environment.
    resolved <- .twasQdsResolveInputs(
        data = data,
        region = region,
        cisWindow = cisWindow,
        jointRegions = jointRegions,
        genotypeFilterArgs = genotypeFilterArgs
    )
    data <- resolved$data
    xRegions <- resolved$xRegions
    # One record of the settings the joint phase and the shared dispatcher
    # both need, so each names them once.
    jointCfg <- list(
        data = data,
        contexts = contexts,
        traitId = traitId,
        cisWindow = cisWindow,
        dataType = dataType,
        verbose = verbose,
        xRegions = xRegions,
        fitRetention = fitRetention,
        seed = seed
    )
    joint <- .twasQdsResolveTokens(
        jointSpecification,
        methods,
        fineMappingResult,
        cvCfg$folds,
        fitFullData,
        mashPrior,
        cvCfg$samplePartition,
        jointCfg
    )
    if (joint$done) {
        return(joint$result)
    }
    # The joint phase consumes the mrmash token; `norm` comes back holding
    # only the methods that still have to go through the per-tuple loop.
    norm <- joint$norm
    samplePartition <- joint$samplePartition
    dataDrivenPriorMatricesCv <- joint$dataDrivenPriorMatricesCv
    grid <- .twasQdsResolveGrid(data, contexts, traitId, region, norm$tokens)
    # The bundles are the user-facing form; below this line the marker and
    # the joint engine that reads it keep plain scalars.
    marker <- .twasQdsMarker(
        crossValidationArgs = list_assign(
            cvCfg,
            samplePartition = samplePartition
        ),
        ensembleArgs = .ensembleResolve(ensembleArgs),
        fitFullData = fitFullData,
        dataType = dataType,
        fitRetention = fitRetention,
        estimatePi = estimatePi,
        verbose = verbose,
        seed = seed
    )
    tw <- .twasQdsRunGrid(
        grid,
        marker,
        data = data,
        xRegions = xRegions,
        norm = norm,
        fineMappingResult = fineMappingResult,
        twasWeights = twasWeights,
        dataDrivenPriorMatricesCv = dataDrivenPriorMatricesCv,
        cisWindow = cisWindow,
        naAction = naAction,
        residualizationArgs = residualizationArgs,
        usePCA = usePCA,
        nPCs = nPCs,
        verbose = verbose
    )
    .twasQdsAssemble(tw, joint$result)
}

#' @rdname twasWeightsPipeline
#' @export
setMethod(
    "twasWeightsPipeline",
    "QtlDataset",
    .twasPipelineQtlDataset
)

# `cisWindow` expands a trait's own coordinates; `region` is literal. Supplying
# both signals a misunderstanding -> reject.
# @noRd
.twasQdsCheckRegionCisWindow <- function(
    region,
    cisWindow
) {
    if (!is.null(region) && !is.null(cisWindow)) {
        msg <- glue(
            "twasWeightsPipeline(QtlDataset): specify either `region` or ",
            "`cisWindow`, not both. `cisWindow` expands each trait's own ",
            "coordinates, whereas `region` is the literal variant window."
        )
        abort(msg)
    }
    invisible(NULL)
}

# fitFullData = FALSE (CV-only) is meaningful only with cross-validation.
# @noRd
.twasQdsCheckFitFull <- function(
    fitFullData,
    cvFolds
) {
    if (!isTRUE(fitFullData) && cvFolds <= 1L) {
        msg <- glue(
            "twasWeightsPipeline: fitFullData = FALSE requires ",
            "cross-validation (cvFolds > 1)."
        )
        abort(msg)
    }
    invisible(NULL)
}

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

# Reject the region + cisWindow combination, derive the X windows, and apply
# the per-call filter overrides to a validated copy of the dataset. Mirrors
# .fmQdsResolveInputs; variant QC is a data property applied identically to
# fine-mapping and TWAS, so there is no TWAS-specific variant filter.
# @noRd
.twasQdsResolveInputs <- function(
    data,
    region,
    cisWindow,
    jointRegions,
    genotypeFilterArgs
) {
    .twasQdsCheckRegionCisWindow(region, cisWindow)
    list(
        data = .qtlApplyFilterOverrides(data, genotypeFilterArgs),
        xRegions = .makeXRegions(region, jointRegions)
    )
}

# Run the joint engine for a QtlDataset with the given spec + token set.
# @noRd
.twasQdsJointDispatch <- function(jointSpec, tokens, cfg) {
    .twasDispatchJointSpecsQtlDataset(
        jointSpec,
        cfg$data,
        tokens,
        cfg$contexts,
        cfg$traitId,
        cfg$cisWindow,
        cfg$dataType,
        cfg$verbose,
        xRegions = cfg$xRegions,
        fitRetention = cfg$fitRetention,
        seed = cfg$seed
    )
}

# Explicit jointSpecification path: run the per-spec axis dispatcher for
# mr.mash. Returns list(done, result, norm): `done` requests an early return
# with `result`; otherwise `norm` is the mrmash-stripped normalization for the
# per-tuple loop below.
# @noRd
.twasQdsJointPhase <- function(parsedJointSpec, norm, cfg) {
    if (length(parsedJointSpec) == 0L) {
        return(list(done = FALSE, result = NULL, norm = norm))
    }
    jointResult <- .twasQdsJointDispatch(
        parsedJointSpec,
        intersect(norm$tokens, "mrmash"),
        cfg
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

# Normalize methods, run the capability / fine-mapping / fit-full gates and the
# mash-prior unpacking, then hand off to the joint phase. Returns the phase
# record extended with the mash-derived sample partition and CV prior matrices
# the per-tuple loop still needs. Mirrors .fmQdsResolveTokens.
# @noRd
.twasQdsResolveTokens <- function(
    jointSpecification,
    methods,
    fineMappingResult,
    cvFolds,
    fitFullData,
    mashPrior,
    samplePartition,
    cfg
) {
    parsedJointSpec <- parseJointSpecification(jointSpecification, cfg$data)
    rawNorm <- .twasNormalizeMethods(methods)
    .twasCheckMethodCapabilities(rawNorm$tokens, "QtlDataset")
    .twasCheckMethodArgsForInput(
        .twasMethodListArgs(rawNorm$methodList),
        "QtlDataset"
    )
    .twasCheckFineMappingMethods(
        rawNorm$tokens,
        fineMappingResult,
        "QtlDataset",
        cvFolds = cvFolds
    )
    .twasQdsCheckFitFull(fitFullData, cvFolds)
    mash <- .twasQdsUnpackMash(mashPrior, samplePartition, rawNorm)
    list_assign(
        .twasQdsJointPhase(parsedJointSpec, mash$norm, cfg),
        samplePartition = mash$samplePartition,
        dataDrivenPriorMatricesCv = mash$dataDrivenPriorMatricesCv
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
    study <- getStudy(data)
    useCtx <- .twasQdsResolveContexts(data, contexts)
    allTraits <- .twasQdsResolveTraits(data, useCtx, traitId, region)
    .twasCheckMultivariateY(tokens, length(allTraits), length(useCtx))
    list(
        study = study,
        useCtx = useCtx,
        allTraits = allTraits,
        multivariate = any(map_lgl(tokens, .twasIsMultivariateToken))
    )
}

# Selected contexts (all when NULL; else validated against the dataset).
# @noRd
.twasQdsResolveContexts <- function(data, contexts) {
    allCtx <- getContexts(data)
    if (is.null(contexts)) {
        return(allCtx)
    }
    bad <- setdiff(contexts, allCtx)
    if (length(bad) > 0L) {
        badStr <- str_flatten(bad, ", ")
        msg <- glue(
            "twasWeightsPipeline(QtlDataset): unknown context(s): {badStr}"
        )
        abort(msg)
    }
    contexts
}

# Traits to iterate: traitId when supplied, else per-context region overlap,
# else every trait in every selected context.
# @noRd
.twasQdsResolveTraits <- function(data, useCtx, traitId, region) {
    perCtxTraits <- map(
        useCtx,
        .twasQdsCtxTraits,
        data = data,
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
    crossValidationArgs,
    ensembleArgs,
    fitFullData,
    dataType,
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
            crossValidationArgs = crossValidationArgs,
            ensembleArgs = ensembleArgs,
            fitFullData = fitFullData,
            dataType = dataType,
            fitRetention = fitRetention,
            standardized = FALSE,
            estimatePi = estimatePi,
            verbose = verbose,
            seed = seed,
            ldSketch = NULL
        )
    )
}

# Which joint cells run, and over which study/context/trait scope. usePCA
# adds top-PC rows ALONGSIDE the per-trait ones, matching
# fineMappingPipeline's c(univRows, pcaRows).
# @noRd
.twasQdsUnivPlan <- function(usePCA, study, useCtx, allTraits) {
    univCell <- .lookupJointCell("univariate", "individual")
    list(
        cells = if (isTRUE(usePCA)) {
            list(univCell, .lookupJointCell("topPc", "individual"))
        } else {
            list(univCell)
        },
        scope = list(
            studies = study,
            contexts = set_names(list(useCtx), study),
            traits = set_names(list(allTraits), study)
        )
    )
}

# Univariate methods ROUTED THROUGH THE ENGINE: one 1-condition group per
# (context, trait), per region -> the SAME per-method fitter (+ ensemble layer
# for >= 2 methods + resume cache) as the joint paths, merged across regions.
# @noRd
.twasQdsUnivariateEngine <- function(
    study,
    useCtx,
    allTraits,
    marker,
    data,
    xRegions,
    norm,
    fineMappingResult,
    twasWeights,
    dataDrivenPriorMatricesCv,
    cisWindow,
    naAction,
    residualizationArgs,
    usePCA,
    nPCs,
    verbose
) {
    plan <- .twasQdsUnivPlan(usePCA, study, useCtx, allTraits)
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
        naAction = naAction,
        residualizationArgs = residualizationArgs,
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
    naAction,
    residualizationArgs,
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
        naAction = naAction,
        residualizationArgs = residualizationArgs,
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
    .rbindTwasWeights(tw, jointResult, ldSketch = NULL)
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
        jointSpecification = NULL,
        fineMappingResult = NULL,
        twasWeights = NULL,
        fitRetention = c("slim", "none", "full"),
        dataType = NULL,
        verbose = 1L,
        panelFilterArgs = PanelFilterParam(),
        crossValidationArgs = CrossValidationParam()
    ) {
        fitRetention <- arg_match(fitRetention)
        .twasWarnUnretainableDetail(fitRetention, methods)
        .cvRefuseOnSumstats(
            crossValidationArgs,
            "twasWeightsPipeline",
            "QtlSumStats"
        )
        .twasPipelineQtlSumStats(
            data = data,
            methods = methods,
            contexts = contexts,
            traitId = traitId,
            jointSpecification = jointSpecification,
            fineMappingResult = fineMappingResult,
            twasWeights = twasWeights,
            fitRetention = fitRetention,
            dataType = dataType,
            verbose = verbose,
            panelFilterArgs = panelFilterArgs
        )
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
        ldSketch = part$ldSketch,
        methodArgs = methodArgs
    )
    c(
        # The univariate builder keeps whatever the per-method fit returns,
        # so it takes no fitRetention.
        exec(
            .twasQssUnivariateRows,
            !!!rowArgs,
            univariateTokens = part$univariateTokens,
            data = cfg$data,
            twasWeights = cfg$twasWeights,
            dataType = cfg$dataType,
            fineMappingResult = cfg$fineMappingResult,
            panelFilterArgs = cfg$panelFilterArgs
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
            dataType = cfg$dataType,
            fitRetention = cfg$fitRetention,
            fineMappingResult = cfg$fineMappingResult,
            panelFilterArgs = cfg$panelFilterArgs
        )
    )
}

.twasPipelineQtlSumStats <- function(
    data,
    methods,
    contexts,
    traitId,
    jointSpecification,
    fineMappingResult,
    twasWeights,
    fitRetention,
    dataType,
    verbose,
    panelFilterArgs
) {
    # summaryStatsQc() is mandatory before twasWeightsPipeline for SumStats
    # input; it also drops variants not present in the ldSketch, so every
    # entry's SNP set is a subset of the ldSketch panel by the time we get here.
    .twasAssertQcd(data)
    # One record of the settings the joint phase and the shared dispatcher
    # both need, so each names them once.
    jointCfg <- list(
        data = data,
        contexts = contexts,
        traitId = traitId,
        dataType = dataType,
        verbose = verbose,
        fitRetention = fitRetention,
        panelFilterArgs = panelFilterArgs
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
        data,
        joint = joint,
        contexts = contexts,
        traitId = traitId,
        twasWeights = twasWeights,
        dataType = dataType,
        fineMappingResult = fineMappingResult,
        fitRetention = fitRetention,
        panelFilterArgs = panelFilterArgs
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
    dataType,
    fineMappingResult,
    fitRetention,
    panelFilterArgs
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
            dataType = dataType,
            fineMappingResult = fineMappingResult,
            fitRetention = fitRetention,
            panelFilterArgs = panelFilterArgs
        )
    )
    .twasQssAssemble(rows, joint$result, part$ldSketch)
}

# Normalize the methods argument into (tokens, methodArgs). The default set
# excludes fine-mapping methods; those must be requested explicitly together
# with a FineMappingResult passed via `fineMappingResult`.
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
    if (.isMethodOptions(methods)) {
        return(list(
            tokens = names(methods),
            methodArgs = map(as.list(methods), as.list)
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

# Run the joint engine for a QtlSumStats with the given spec + token set.
# @noRd
.twasQssJointDispatch <- function(jointSpec, tokens, cfg) {
    .twasDispatchJointSpecsQtlSumStats(
        jointSpec,
        cfg$data,
        tokens,
        cfg$contexts,
        cfg$traitId,
        cfg$dataType,
        cfg$verbose,
        fitRetention = cfg$fitRetention,
        panelFilterArgs = cfg$panelFilterArgs %||% PanelFilterParam()
    )
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
    jointResult <- .twasQssJointDispatch(
        parsedJointSpec,
        intersect(tokens, "mrmash"),
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
    selRows <- .twasQssSelectRows(
        data,
        contextCol,
        traitCol,
        contexts,
        traitId
    )
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
        univariateTokens = univariateTokens,
        ldSketch = getLdSketch(data)
    )
}

# Row indices matching the contexts / traitId filters (error if none).
# @noRd
.twasQssSelectRows <- function(
    data,
    contextCol,
    traitCol,
    contexts,
    traitId
) {
    byContext <- if (is.null(contexts)) {
        seq_len(nrow(data))
    } else {
        which(is_in(contextCol, contexts))
    }
    selRows <- if (is.null(traitId)) {
        byContext
    } else {
        byContext[is_in(traitCol[byContext], traitId)]
    }
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
.twasRowRecord <- function(study, context, trait, method, entry) {
    list(
        study = study,
        context = context,
        trait = trait,
        method = method,
        entry = entry
    )
}

# Assemble a TwasWeights collection from a flat list of row records.
# @noRd
.twasRowsToWeights <- function(rows, ldSketch) {
    if (length(rows) == 0L) {
        return(NULL)
    }
    TwasWeights(
        study = map_chr(rows, "study"),
        context = map_chr(rows, "context"),
        trait = map_chr(rows, "trait"),
        method = map_chr(rows, "method"),
        entry = map(rows, "entry"),
        ldSketch = ldSketch
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

# Run a weight function, warning (with `errPrefix`) and returning NULL on error.
# @noRd
#' @importFrom rlang try_fetch
.twasTryWeights <- function(fn, stat, ldMat, userArgs, errPrefix) {
    try_fetch(
        {
            wfn <- get(fn, mode = "function")
            wArgs <- .twasWeightCallArgs(
                fn,
                list(stat = stat, LD = ldMat),
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
    ldSketch,
    twasWeights,
    dataType,
    fineMappingResult,
    methodArgs,
    panelFilterArgs
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
        ldSketch = ldSketch,
        twasWeights = twasWeights,
        dataType = dataType,
        fineMappingResult = fineMappingResult,
        methodArgs = methodArgs,
        panelFilterArgs = panelFilterArgs
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
    ldSketch,
    twasWeights,
    dataType,
    fineMappingResult,
    methodArgs,
    panelFilterArgs
) {
    st <- studyCol[i]
    ctx <- contextCol[i]
    tr <- traitCol[i]
    cached <- .twasQssCacheSplit(twasWeights, st, ctx, tr, univariateTokens)
    cachedRows <- cached$rows
    toFit <- cached$toFit
    if (length(toFit) == 0L) {
        return(unname(cachedRows))
    }
    fitCtx <- .twasQssUnivariateFitCtx(
        data,
        st,
        ctx,
        tr,
        ldSketch,
        cutoffs = .panelCutoffs(panelFilterArgs)
    )
    fitted <- compact(map(
        toFit,
        .twasQssUnivariateFitOne,
        st = st,
        ctx = ctx,
        tr = tr,
        fitCtx = fitCtx,
        methodArgs = methodArgs,
        fineMappingResult = fineMappingResult,
        dataType = dataType
    ))
    c(unname(cachedRows), fitted)
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
    ldSketch,
    cutoffs = NULL
) {
    allDf <- getSumStatsDf(
        data,
        study = st,
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
        .panelKeepMask(allDf$variant_id, ldSketch, cutoffs, label),
        ,
        drop = FALSE
    ]
    variantIds <- df$variant_id
    varY <- getVarY(data, study = st, context = ctx, trait = tr) %||% 1
    stat <- list(
        z = df$z,
        n = stats::median(df$N, na.rm = TRUE),
        varY = varY,
        variantNames = variantIds
    )
    list(
        variantIds = variantIds,
        stat = stat,
        ldMat = .ldFromSketch(
            ldSketch,
            variantIds,
            label = "twasWeightsPipeline"
        )
    )
}

# Fit one univariate method for one entry -> a row record, or NULL on skip.
# @noRd
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
        study = st,
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

.twasQssUnivariateFitOne <- function(
    tk,
    st,
    ctx,
    tr,
    fitCtx,
    methodArgs,
    fineMappingResult,
    dataType
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
        fitCtx$ldMat,
        userArgs,
        .twasFitErrUniv(tk, st, ctx, tr)
    )
    if (is.null(weights)) {
        return(NULL)
    }
    .twasQssRowFromWeights(weights, st, ctx, tr, tk, fitCtx, dataType)
}

# The weights come back with the fit hung off them as an attribute; the row
# stores the two separately, so they are split apart here.
# @noRd
.twasQssRowFromWeights <- function(weights, st, ctx, tr, tk, fitCtx, dataType) {
    .twasRowRecord(
        st,
        ctx,
        tr,
        tk,
        twasWeightsRow(
            variantIds = fitCtx$variantIds,
            weights = as.numeric(`attr<-`(weights, "fit", NULL)),
            fits = attr(weights, "fit"),
            cvResult = NULL,
            standardized = TRUE,
            dataType = dataType
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
    ldSketch,
    methodArgs,
    fitRetention,
    fineMappingResult,
    dataType,
    panelFilterArgs,
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
        ldSketch = ldSketch,
        methodArgs = methodArgs,
        fitRetention = fitRetention,
        fineMappingResult = fineMappingResult,
        dataType = dataType,
        panelFilterArgs = panelFilterArgs
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
    ldSketch,
    methodArgs,
    fitRetention,
    fineMappingResult,
    dataType,
    panelFilterArgs,
    twasWeights = NULL
) {
    if (length(gIdx) < 2L) {
        return(list())
    }
    st <- studyCol[gIdx[[1L]]]
    tr <- traitCol[gIdx[[1L]]]
    ctxNames <- contextCol[gIdx]
    mvStat <- .twasQssMultivariateStat(
        data,
        st,
        tr,
        ctxNames,
        ldSketch = ldSketch,
        cutoffs = .panelCutoffs(panelFilterArgs)
    )
    ldMat <- .ldFromSketch(
        ldSketch,
        mvStat$variantIds,
        label = "twasWeightsPipeline"
    )
    list_flatten(map(
        multivariateTokens,
        .twasQssMultivariateFitOne,
        twasWeights = twasWeights,
        st = st,
        tr = tr,
        ctxNames = ctxNames,
        mvStat = mvStat,
        ldMat = ldMat,
        methodArgs = methodArgs,
        fitRetention = fitRetention,
        fineMappingResult = fineMappingResult,
        dataType = dataType
    ))
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
    ldSketch = NULL,
    cutoffs = NULL
) {
    allDf <- getSumStatsDf(
        data,
        study = st,
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
    keep <- .panelKeepMask(allDf$variant_id, ldSketch, cutoffs, label)
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

# Read each context's z / N into the shared matrix. Every context must present
# the same SNPs in the same order -- the multivariate fitters index z by
# position, not by name, so a differing order would silently pair one context's
# variant with another's.
# @noRd
# One context's z column and median N, checked against the shared SNP order.
# @noRd
.twasQssContextStats <- function(ctx, data, st, tr, variantIds) {
    d <- getSumStatsDf(
        data,
        study = st,
        context = ctx,
        trait = tr,
        require = c("Z", "N"),
        derive = "zFromBetaSe"
    )
    kept <- d[is_in(d$variant_id, variantIds), , drop = FALSE]
    .twasQssCheckSnpOrder(kept$variant_id, variantIds, st, tr)
    list(z = kept$z, n = stats::median(kept$N, na.rm = TRUE))
}

# Every context shares one SNP order (asserted per context), which is what
# lets the columns simply be laid side by side rather than filled into a
# preallocated matrix.
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
        "ldSketch on every entry."
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
    baseArgs <- .twasMrmashRetainDefaults(
        .twasUserArgs(methodArgs, tk),
        spec$adapter,
        tk,
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
# list on skip).
# @noRd
.twasQssMultivariateFitOne <- function(
    tk,
    st,
    tr,
    ctxNames,
    mvStat,
    ldMat,
    methodArgs,
    fitRetention,
    fineMappingResult,
    dataType,
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
    # NULL when the pre-fit a fine-mapping token needs is missing.
    if (is.null(userArgs)) {
        return(list())
    }
    weights <- .twasTryWeights(
        spec$fn,
        mvStat$stat,
        ldMat,
        userArgs,
        .twasFitErrMv(tk, st, tr)
    )
    if (is.null(weights)) {
        return(list())
    }
    # The fit is hung off the weights as an attribute; the rows store the
    # two separately.
    wMatrix <- if (is.matrix(weights)) weights else as.matrix(weights)
    fitAttr <- attr(wMatrix, "fit")
    bare <- `attr<-`(wMatrix, "fit", NULL)
    .twasMvContextRows(bare, fitAttr, ctxNames, mvStat, st, tr, tk, dataType)
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
        study = st,
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
    dataType
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
        dataType = dataType
    )
}

# Combine the per-tuple result with any joint result (error if both empty).
# @noRd
.twasQssAssemble <- function(
    rows,
    jointResult,
    ldSketch
) {
    perTupleResult <- .twasRowsToWeights(rows, ldSketch)
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
    .rbindTwasWeights(perTupleResult, jointResult, ldSketch = ldSketch)
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
# @noRd
.twasPerStudy <- function(qd, cfg) {
    # Checked forwarding, not the `...`/dotArgs channel this used to splice:
    # that is what made multi-study TWAS inherit the per-study CV default and
    # silently drop residualization / fitRetention.
    # Every setting is named -- NOT the `...`/dotArgs channel this once
    # spliced, which is what made multi-study TWAS inherit the per-study CV
    # default and silently drop residualization / fitRetention. A name this
    # entry point does not accept is an "unused argument" error; `cfg`'s
    # panelFilterArgs belongs to the summary-statistics path and is
    # deliberately not forwarded.
    twasWeightsPipeline(
        data = qd,
        jointSpecification = NULL,
        fitRetention = .twasRetentionEnum(cfg),
        methods = cfg$methods,
        contexts = cfg$contexts,
        traitId = cfg$traitId,
        region = cfg$region,
        cisWindow = cfg$cisWindow,
        jointRegions = cfg$jointRegions,
        fineMappingResult = cfg$fineMappingResult,
        twasWeights = cfg$twasWeights,
        naAction = cfg$naAction,
        verbose = cfg$verbose,
        crossValidationArgs = cfg$crossValidationArgs,
        ensembleArgs = cfg$ensembleArgs,
        seed = cfg$seed,
        genotypeFilterArgs = cfg$genotypeFilterArgs,
        mashPrior = cfg$mashPrior,
        fitFullData = cfg$fitFullData,
        estimatePi = cfg$estimatePi,
        dataType = cfg$dataType,
        residualizationArgs = cfg$residualizationArgs
    )
}

# Embedded-sumstats TWAS-weights worker for .multiStudyPipelineDriver.
# @noRd
.twasSumStats <- function(ss, cfg) {
    # Summary statistics carry no genotypes: nothing to select a region
    # from, residualize, cross-validate, stack or fit on full data. `cfg`
    # carries those for the individual-level sibling and they are not
    # forwarded here.
    twasWeightsPipeline(
        data = ss,
        jointSpecification = NULL,
        fitRetention = .twasRetentionEnum(cfg),
        methods = cfg$methods,
        contexts = cfg$contexts,
        traitId = cfg$traitId,
        fineMappingResult = cfg$fineMappingResult,
        twasWeights = cfg$twasWeights,
        verbose = cfg$verbose,
        panelFilterArgs = cfg$panelFilterArgs,
        dataType = cfg$dataType
    )
}

# The config carries the engine-facing pair; the per-study recursion re-enters
# the public method, which takes the single enum. One place converts back.
# @noRd
.twasRetentionEnum <- function(cfg) {
    if (!isTRUE(cfg$fitRetention)) {
        return("none")
    }
    cfg$fitRetention %||% "slim"
}

# Only the mr.mash engines distinguish the retention levels; every other
# weight method keeps its fit whole or not at all, so "full" is the same as
# "slim" for them. Erroring when NO requested method can honour it, and
# warning when only some can, is the alternative to dropping it in silence.
# @noRd
.twasWarnUnretainableDetail <- function(fitRetention, methods) {
    if (!identical(fitRetention, "full")) {
        return(invisible(NULL))
    }
    tokens <- .twasMethodTokensFromArg(methods)
    if (length(tokens) == 0L) {
        return(invisible(NULL))
    }
    honours <- str_detect(tokens, regex("mrmash", ignore_case = TRUE))
    if (!any(honours)) {
        abort(glue(
            "twasWeightsPipeline: fitRetention = \"full\" is only honoured ",
            "by the mr.mash methods; {str_flatten(tokens, ', ')} cannot vary ",
            "the retained detail. Use fitRetention = \"slim\"."
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

# ---- MultiStudyQtlDataset pipeline worker + phase helpers ------------------

# Method tokens from a `methods` arg: a character vector as-is; a named list ->
# canonical bare tokens; otherwise empty.
# @noRd
.twasMethodTokensFromArg <- function(methods) {
    if (.twasMethodsIsParam(methods)) {
        return(.methodsParamTokens(methods))
    }
    methods <- .twasMethodsAsList(methods)
    if (is.character(methods)) {
        methods
    } else if (is.list(methods)) {
        str_remove(names(methods), "(_weights|Weights)$")
    } else {
        character(0)
    }
}

# Drop mrmash (handled by the joint dispatcher) from a `methods` arg.
# @noRd
# A methods argument as an ordinary named list, whatever form it arrived in.
# A MethodOptions record is a SimpleList, so `is.list()` is FALSE for it and the
# helpers below would otherwise fall through to their "leave it alone" branch.
# @noRd
.twasMethodsAsList <- function(methods) {
    if (.isMethodOptions(methods)) as.list(methods) else methods
}

# A MethodsSelectionParam keeps its structure through these helpers: the
# token questions read across all three slots, and dropping a method drops
# it from each, so the per-component recursion still receives a Param that
# knows which path each entry is for.
# @noRd
.twasMethodsIsParam <- function(methods) {
    is(methods, "MethodsSelectionParam")
}

.twasMsStripMrmash <- function(methods) {
    if (.twasMethodsIsParam(methods)) {
        return(.methodsParamDrop(methods, "mrmash"))
    }
    methods <- .twasMethodsAsList(methods)
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
    methods <- .twasMethodsAsList(methods)
    (is.character(methods) || is.list(methods)) && length(methods) == 0L
}

.twasPipelineMultiStudy <- function(
    data,
    methods = "default",
    contexts = NULL,
    traitId = NULL,
    region = NULL,
    cisWindow = NULL,
    genotypeFilterArgs = GenotypeFilterParam(),
    panelFilterArgs = PanelFilterParam(),
    jointRegions = FALSE,
    jointSpecification = NULL,
    fineMappingResult = NULL,
    twasWeights = NULL,
    mashPrior = NULL,
    fitFullData = TRUE,
    estimatePi = TRUE,
    fitRetention = c("slim", "none", "full"),
    dataType = NULL,
    naAction = c("drop", "impute"),
    verbose = 1,
    residualizationArgs = ResidualizationParam(),
    crossValidationArgs = CrossValidationParam(),
    ensembleArgs = EnsembleParam(),
    seed = NULL
) {
    naAction <- arg_match(naAction)
    fitRetention <- arg_match(fitRetention)
    .twasWarnUnretainableDetail(fitRetention, methods)
    if (!is.null(region) && !is.null(cisWindow)) {
        msg <- glue(
            "twasWeightsPipeline(MultiStudyQtlDataset): specify either ",
            "`region` or `cisWindow`, not both."
        )
        abort(msg)
    }
    # One record of the settings the joint phase and the shared dispatcher
    # both need, so each names them once.
    jointCfg <- list(
        data = data,
        contexts = contexts,
        traitId = traitId,
        cisWindow = cisWindow,
        verbose = verbose,
        xRegions = .makeXRegions(region, jointRegions),
        fitRetention = fitRetention,
        seed = seed
    )
    joint <- .twasMsResolveTokens(
        jointSpecification,
        methods,
        fineMappingResult,
        jointCfg
    )
    if (joint$done) {
        return(joint$result)
    }
    .twasMsDriver(
        data = data,
        contexts = contexts,
        traitId = traitId,
        cisWindow = cisWindow,
        region = region,
        jointRegions = jointRegions,
        fineMappingResult = fineMappingResult,
        twasWeights = twasWeights,
        naAction = naAction,
        verbose = verbose,
        crossValidationArgs = crossValidationArgs,
        ensembleArgs = ensembleArgs,
        seed = seed,
        genotypeFilterArgs = genotypeFilterArgs,
        panelFilterArgs = panelFilterArgs,
        mashPrior = mashPrior,
        fitFullData = fitFullData,
        estimatePi = estimatePi,
        fitRetention = fitRetention,
        dataType = dataType,
        residualizationArgs = residualizationArgs,
        joint$result,
        joint$methods
    )
}

#' @rdname twasWeightsPipeline
#' @export
setMethod(
    "twasWeightsPipeline",
    "MultiStudyQtlDataset",
    .twasPipelineMultiStudy
)

# Run the joint engine for a MultiStudyQtlDataset with the given spec + token
# set. The NULL in the dataType slot is deliberate: the per-component
# recursion resolves each study's own data type.
# @noRd
.twasMsJointDispatch <- function(jointSpec, tokens, cfg) {
    .twasDispatchJointSpecsMultiStudy(
        jointSpec,
        cfg$data,
        tokens,
        cfg$contexts,
        cfg$traitId,
        cfg$cisWindow,
        NULL,
        cfg$verbose,
        xRegions = cfg$xRegions,
        fitRetention = cfg$fitRetention,
        seed = cfg$seed
    )
}

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
    jointResult <- .twasMsJointDispatch(
        parsedJointSpec,
        intersect(.twasMethodTokensFromArg(methods), "mrmash"),
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
        !is.null(getSumStats(cfg$data))
    )
    .twasCheckFineMappingMethods(
        .twasMethodTokensFromArg(methods),
        fineMappingResult,
        "MultiStudyQtlDataset"
    )
    .twasMsJointPhase(parsedJointSpec, methods, cfg)
}

# Run the per-study / per-component recursion via the shared multi-study driver.
# @noRd
.twasMsDriver <- function(
    data,
    contexts,
    traitId,
    cisWindow,
    region,
    jointRegions,
    fineMappingResult,
    twasWeights,
    naAction,
    verbose,
    crossValidationArgs,
    ensembleArgs,
    seed,
    genotypeFilterArgs,
    panelFilterArgs,
    mashPrior,
    fitFullData,
    estimatePi,
    fitRetention,
    dataType,
    residualizationArgs,
    jointResult,
    methods
) {
    cfg <- list(
        methods = methods,
        contexts = contexts,
        traitId = traitId,
        region = region,
        cisWindow = cisWindow,
        jointRegions = jointRegions,
        fineMappingResult = fineMappingResult,
        twasWeights = twasWeights,
        naAction = naAction,
        verbose = verbose,
        crossValidationArgs = crossValidationArgs,
        ensembleArgs = ensembleArgs,
        seed = seed,
        genotypeFilterArgs = genotypeFilterArgs,
        panelFilterArgs = panelFilterArgs,
        mashPrior = mashPrior,
        fitFullData = fitFullData,
        estimatePi = estimatePi,
        fitRetention = fitRetention,
        dataType = dataType,
        residualizationArgs = residualizationArgs
    )
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

# Ensemble stacking objective (sum of squared residuals). `...` absorbs the
# gradient's extra optim args (PtP, Pty).
# @noRd
.ensembleObj <- function(Pvalid, yObs) {
    # Captured by name, and forced here so the closure holds values rather
    # than promises into a frame that has already returned.
    force(Pvalid)
    force(yObs)
    function(z) sum((yObs - Pvalid %*% z)^2)
}

# Gradient of the ensemble stacking objective, same construction. Building
# both as closures of exactly what each needs is what lets optim() be called
# with no `...`: previously it forwarded the union of both callbacks'
# arguments to both, and each absorbed the other's half in a `...` tail.
# @noRd
.ensembleGrad <- function(PtP, Pty) {
    force(PtP)
    force(Pty)
    function(z) as.vector(2 * (PtP %*% z - Pty))
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
    PtP <- crossprod(Pvalid)
    Pty <- as.vector(crossprod(Pvalid, yObs))

    fit <- try_fetch(
        optim(
            par = rep(1 / Kvalid, Kvalid),
            fn = .ensembleObj(Pvalid, yObs),
            gr = .ensembleGrad(PtP, Pty),
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
    # Engine-call exception: glmnet is a constrained-regression solver here,
    # fitting the ensemble's stacking weights -- a different role from its use
    # as a weight-learning engine in regularizedRegressionWrappers.R. The
    # stacking step is pipeline logic, so the call stays with it.
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        abort("Package 'glmnet' is required for solver='glmnet'.")
    }

    fit <- try_fetch(
        glmnet::cv.glmnet(
            x = Pvalid,
            y = yObs,
            lower.limits = 0,
            alpha = alpha,
            intercept = FALSE
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

# Assemble one dataset's (samples x methods) prediction matrix.
# @noRd
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

# Zeta-weighted sum of the matched weight matrices; univariate -> named vector.
# @noRd
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

# The entry in one region's collection matching a (study, context, trait,
# method) key, or NULL when that region lacks it.
# @noRd
.twasEntryMatchingKey <- function(tw, key) {
    hit <- which(
        as.character(tw$study) == key[[1L]] &
            as.character(tw$context) == key[[2L]] &
            as.character(tw$trait) == key[[3L]] &
            as.character(tw$method) == key[[4L]]
    )
    if (length(hit)) .twrRowParts(tw, hit[[1L]]) else NULL
}

# The merged entry for base row `r`: gather that key's entry from every region
# and concatenate them.
# @noRd
.twasMergedEntryForRow <- function(r, base, twList, regionLabels) {
    key <- c(
        as.character(base$study[[r]]),
        as.character(base$context[[r]]),
        as.character(base$trait[[r]]),
        as.character(base$method[[r]])
    )
    perRegion <- map(twList, .twasEntryMatchingKey, key = key)
    .twasMergeRegionEntries(perRegion, regionLabels)
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

# One region's multivariate-grid joint fit (region `bi` of ctx$xRegions).
# @noRd
.twasMvGridRegion <- function(bi, synthSpec, marker, ctx, traits) {
    .runJointSpecs(
        synthSpec,
        ctx$data,
        dataForm = "individual",
        pipeline = marker,
        jointMethods = ctx$norm$tokens,
        contexts = ctx$useCtx,
        traitIds = traits,
        args = list(
            methodList = ctx$norm$methodList,
            fineMappingResult = ctx$fineMappingResult,
            dataDrivenPriorMatricesCv = ctx$dataDrivenPriorMatricesCv,
            cisWindow = ctx$cisWindow,
            region = ctx$xRegions[[bi]],
            regionIndex = bi,
            nRegions = length(ctx$xRegions),
            residualizationArgs = ctx$residualizationArgs,
            verbose = ctx$verbose
        )
    )
}

# The traits of one context: traitId when supplied, else region overlap, else
# every trait in the context.
# @noRd
.twasQdsCtxTraits <- function(ctx, data, traitId, region) {
    se <- getPhenotypes(data, contexts = ctx)
    ids <- rownames(se)
    if (!is.null(traitId)) {
        intersect(ids, traitId)
    } else if (!is.null(region)) {
        rr <- SummarizedExperiment::rowRanges(se)
        ids[IRanges::overlapsAny(rr, region)]
    } else {
        ids
    }
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
    naAction,
    residualizationArgs,
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
        naAction = naAction,
        residualizationArgs = residualizationArgs,
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
    dataType
) {
    .twasRowRecord(
        st,
        ctxNames[[kk]],
        tr,
        tk,
        twasWeightsRow(
            variantIds = mvStat$variantIds,
            weights = as.numeric(weights[, kk]),
            fits = if (kk == 1L) fitAttr else NULL,
            cvResult = NULL,
            standardized = TRUE,
            dataType = dataType
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
#'     qtlDatasetMethods = list(lasso = GlmnetOptions(alpha = 0.5)),
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
.ensembleAssertCv <- function(ensembleArgs, crossValidationArgs) {
    if (!isTRUE(ensembleArgs$enabled) || .cvEnabled(crossValidationArgs)) {
        return(invisible(NULL))
    }
    abort(glue(
        "twasWeightsPipeline: EnsembleParam(enabled = TRUE) needs ",
        "out-of-fold predictions to stack, so it requires ",
        "CrossValidationParam(folds >= 2); got folds = ",
        "{crossValidationArgs$folds %||% 0}."
    ))
}

# The ensemble settings as a plain list, so a cfg record can carry the group
# as one named field instead of four loose ones.
# @noRd
.ensembleResolve <- function(ensembleArgs) {
    list(
        enabled = isTRUE(ensembleArgs$enabled),
        r2Threshold = ensembleArgs$r2Threshold %||% 0.01,
        solver = ensembleArgs$solver %||% "quadprog",
        alpha = ensembleArgs$alpha %||% 1
    )
}
