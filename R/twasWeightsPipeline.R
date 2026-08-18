# =============================================================================
# Helpers + S4 dispatch surface for twasWeightsPipeline
# =============================================================================

# Concatenate two TwasWeights collections row-wise, carrying forward every
# column (delegates to the generic `.rbindCollections`, which unions columns
# and pads a side lacking an optional column such as joint* / region).
# @noRd
.rbindTwasWeights <- function(a, b, ldSketch = NULL) {
    if (!is(a, "TwasWeights") || !is(b, "TwasWeights")) {
        abort(".rbindTwasWeights expects two TwasWeights inputs.")
    }
    # Carry forward every column (joint*, region, ...) via the generic combine.
    .rbindCollections(list(a, b), ldSketch = ldSketch)
}

# Normalize combine() varargs: accept either N objects or a single list of
# them; drop NULLs; require at least one input of the expected class `cls`.
.asCombineList <- function(parts, cls, fn) {
    if (
        length(parts) == 1L &&
            is.list(parts[[1L]]) &&
            !methods::is(parts[[1L]], cls)
    ) {
        parts <- parts[[1L]]
    }
    parts <- compact(parts)
    if (length(parts) == 0L) {
        msg <- glue("{fn}: nothing to combine (need at least one {cls}).")
        abort(msg)
    }
    if (!all(map_lgl(parts, methods::is, cls))) {
        msg <- glue("{fn}: every input must be a {cls}.")
        abort(msg)
    }
    parts
}

#' Combine TwasWeights collections
#'
#' Row-bind two or more \code{\link{TwasWeights}} collections into one -- e.g.
#' assembling per-gene weight sets into a single per-region collection for
#' cTWAS. Joint-specification metadata columns are carried through.
#'
#' @param ... Two or more \code{TwasWeights} objects, or a single \code{list} of
#'   them.
#' @param ldSketch Optional \code{\link{GenotypeHandle}} to attach to the
#'   combined collection. Default \code{NULL}. Applied when combining two or
#'   more inputs; a single input is returned unchanged.
#' @return A single combined \code{TwasWeights}.
#' @seealso \code{\link{combineFineMappingResults}}
#' @examples
#' twe <- TwasWeightsEntry(variantIds = paste0("v", 1:4), weights = rep(0.1, 4),
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
    reduce(parts, .rbindTwasWeights, ldSketch = ldSketch)
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

# Concatenate one method's per-region TwasWeightsEntry payloads into a single
# entry. Variants/weights are stacked (regions are disjoint), the per-region
# fits are kept as a named list, and cvResult becomes the flat per-region
# reporting data.frame.
.twasMergeRegionEntries <- function(entries, regionLabels) {
    keep <- !map_lgl(entries, is.null)
    entries <- entries[keep]
    regionLabels <- regionLabels[keep]
    if (length(entries) == 0L) {
        return(NULL)
    }
    if (length(entries) == 1L) {
        return(entries[[1L]])
    }
    wList <- map(entries, getWeights)
    weights <- if (is.matrix(wList[[1L]])) {
        exec(rbind, !!!wList)
    } else {
        unlist(wList, use.names = FALSE)
    }
    TwasWeightsEntry(
        variantIds = unlist(map(entries, getVariantIds), use.names = FALSE),
        weights = weights,
        fits = set_names(map(entries, getFits), regionLabels),
        cvResult = .twasRegionCvDf(entries, regionLabels),
        standardized = getStandardized(entries[[1L]]),
        dataType = getDataType(entries[[1L]])
    )
}

# Merge per-region TwasWeights collections (same study/context/trait, same
# methods) into one collection by concatenating each method's entry.
.twasMergeRegions <- function(twList, regionLabels) {
    keep <- !map_lgl(twList, is.null)
    twList <- twList[keep]
    regionLabels <- regionLabels[keep]
    if (length(twList) == 0L) {
        return(NULL)
    }
    if (length(twList) == 1L) {
        return(twList[[1L]])
    }
    base <- twList[[1L]]
    mergedEntries <- map(
        seq_along(base$method),
        .twasMergedEntryForRow,
        base = base,
        twList = twList,
        regionLabels = regionLabels
    )
    TwasWeights(
        study = as.character(base$study),
        context = as.character(base$context),
        trait = as.character(base$trait),
        method = as.character(base$method),
        entry = mergedEntries
    )
}

# Unpack a MashPrior input into the internal twasWeightsPipeline arguments:
# $fullPrior -> mr.mash full-data dataDrivenPriorMatrices
# $dataDrivenPriorMatricesCv -> per-fold priors for twasWeightsCv
# $samplePartition -> the CV folds (an explicit `samplePartition` arg wins;
# otherwise the partition the per-fold priors were computed on) NULL input
# returns all-NULL, preserving the supplied samplePartition.
# @noRd
.unpackMashPrior <- function(mashPrior, samplePartition = NULL) {
    if (is.null(mashPrior)) {
        return(list(
            fullPrior = NULL,
            dataDrivenPriorMatricesCv = NULL,
            samplePartition = samplePartition
        ))
    }
    if (!is(mashPrior, "MashPrior")) {
        abort("`mashPrior` must be a MashPrior object (see ?MashPrior).")
    }
    cvFits <- getCvFits(mashPrior)
    perFold <- if (!is.null(cvFits)) cvFits$perFoldFits else NULL
    sp <- samplePartition
    if (is.null(sp) && !is.null(cvFits) && !is.null(cvFits$samplePartition)) {
        sp <- cvFits$samplePartition
    }
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
        sumstatImpl = "mrAshRssWeights",
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
    dpr_gibbs = list(
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
    dpr_vb = list(
        individualImpl = "dprVbWeights",
        sumstatImpl = NULL,
        multivariate = FALSE
    ),
    dpr_adaptive_gibbs = list(
        individualImpl = "dprAdaptiveGibbsWeights",
        sumstatImpl = NULL,
        multivariate = FALSE
    ),
    # qgg Bayes alphabet -- individual-only until qgg CRAN release.
    bayes_a = list(
        individualImpl = "bayesAWeights",
        sumstatImpl = NULL,
        multivariate = FALSE
    ),
    bayes_b = list(
        individualImpl = "bayesBWeights",
        sumstatImpl = NULL,
        multivariate = FALSE
    ),
    bayes_c = list(
        individualImpl = "bayesCWeights",
        sumstatImpl = NULL,
        multivariate = FALSE
    ),
    bayes_l = list(
        individualImpl = "bLassoWeights",
        sumstatImpl = NULL,
        multivariate = FALSE
    ),
    bayes_n = list(
        individualImpl = "bayesNWeights",
        sumstatImpl = NULL,
        multivariate = FALSE
    ),
    bayes_r = list(
        individualImpl = "bayesRWeights",
        sumstatImpl = NULL,
        multivariate = FALSE
    ),
    b_lasso = list(
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
.twasNormalizeMethods <- function(methods) {
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
    if (is.list(methods)) {
        return(.twasNormalizeListMethods(methods))
    }
    msg <- glue(
        "`methods` must be a character vector, preset string, or ",
        "named list."
    )
    abort(msg)
}

# Character `methods`: resolve regular tokens via .twasMethodLookup and append
# empty stub entries for fine-mapping tokens with no learner counterpart (e.g.
# fsusie) so the downstream gate can produce a method-specific error rather than
# "Unknown TWAS method".
# @noRd
.twasNormalizeCharMethods <- function(methods) {
    fmExtra <- setdiff(
        intersect(methods, .twasFineMappingTokens()),
        .twasKnownMethodLookupNames()
    )
    regular <- setdiff(methods, fmExtra)
    methodList <- if (length(regular) > 0L) {
        .twasMethodLookup(regular)
    } else {
        list()
    }
    for (tk in fmExtra) {
        methodList[[str_c(tk, "_weights")]] <- list()
    }
    # Tokens come from the user input (canonical camelCase) -- the snake keys in
    # methodList are an internal detail of learnTwasWeights.
    list(tokens = unique(methods), methodList = methodList)
}

# Named-list `methods`: re-key each entry to its canonical `<token>_weights`
# name, merge the caller's kwargs over the method's defaults, and carry the
# token->impl map as an "impl" attribute (without this, downstream
# .resolveMethodFunction falls back to the bare token, which is not a function).
# @noRd
.twasNormalizeListMethods <- function(methods) {
    methodList <- list()
    for (tk in names(methods)) {
        base <- tryCatch(.twasMethodLookup(tk), error = function(e) NULL)
        if (is.null(base)) {
            # Fine-mapping-only / unknown token with no learner default: keep
            # as-is (the downstream capability gate produces the message).
            methodList[[tk]] <- methods[[tk]]
            next
        }
        snake <- names(base)[[1L]]
        merged <- modifyList(base[[snake]], methods[[tk]])
        attr(merged, "impl") <- attr(base[[snake]], "impl")
        methodList[[snake]] <- merged
    }
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
        "bayes_r",
        "bayes_l",
        "bayes_a",
        "bayes_b",
        "bayes_c",
        "bayes_n",
        "b_lasso",
        "dpr_vb",
        "dpr_gibbs",
        "dpr_adaptive_gibbs",
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
        bayes_r = "bayes_r",
        bayes_l = "bayes_l",
        bayes_a = "bayes_a",
        bayes_b = "bayes_b",
        bayes_c = "bayes_c",
        bayes_n = "bayes_n",
        b_lasso = "b_lasso",
        dpr_vb = "dpr_vb",
        dpr_gibbs = "dpr_gibbs",
        dpr_adaptive_gibbs = "dpr_adaptive_gibbs",
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
    present <- character(0)
    for (canonical in .twasFineMappingTokens()) {
        candidates <- str_to_lower(c(
            canonical,
            str_c(
                str_to_lower(str_sub(canonical, 1L, 1L)),
                str_sub(canonical, 2L)
            ),
            str_replace_all(canonical, "([A-Z])", "_\\1")
        ))
        if (any(is_in(methods, candidates))) present <- c(present, canonical)
    }
    present
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
.twasCheckFineMappingMethods <- function(tokens, fineMappingResult, inputKind) {
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
    invisible(NULL)
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
    if (!is(fineMappingResult, "FineMappingResultBase")) {
        abort("`fineMappingResult` must be a FineMappingResult or NULL.")
    }
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

# Extract a correlation matrix from a GenotypeHandle (LD sketch) for the
# variant subset given by `variantIds`. Thin wrapper over the shared
# `.ldFromSketch` helper.
# @noRd
.twasLdFromSketch <- function(ldSketch, variantIds) {
    .ldFromSketch(ldSketch, variantIds, label = ".twasLdFromSketch")
}

# Optional resume-cache lookup for twasWeightsPipeline. Returns the
# matching TwasWeightsEntry from `twasWeights` for the tuple (study,
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
    twasWeights$entry[[idx[[1L]]]]
}

# Build a TwasWeights collection from a list of cached entries keyed by
# short-method-name (the value of the `method` column). Helper for the
# resume-cache short-circuit path in twasWeightsPipeline.
# @noRd
.twasBuildFromCachedRows <- function(
    cachedRows,
    study,
    context,
    trait,
    ldSketch = NULL
) {
    if (length(cachedRows) == 0L) {
        return(NULL)
    }
    TwasWeights(
        study = rep(study, length(cachedRows)),
        context = rep(context, length(cachedRows)),
        trait = rep(trait, length(cachedRows)),
        method = names(cachedRows),
        entry = unname(cachedRows),
        ldSketch = ldSketch
    )
}

# Convert a FineMappingResult (single-method susie/susie_inf row matched
# to the requested study/context/trait) into a `fittedModels` list
# suitable for `learnTwasWeights`. Pulls the trimmedFit from the matching
# entry. Returns a (possibly empty) list.
# @noRd
.twasFineMappingFits <- function(fineMappingResult, study, context, trait) {
    if (is.null(fineMappingResult)) {
        return(list())
    }
    if (!is(fineMappingResult, "FineMappingResultBase")) {
        abort("`fineMappingResult` must be a FineMappingResult or NULL.")
    }
    out <- list()
    methods <- as.character(fineMappingResult$method)
    for (canonical in c("susie", "susieInf", "susieAsh", "mvsusie", "fsusie")) {
        candidates <- c(
            canonical,
            str_c(
                str_to_lower(str_sub(canonical, 1L, 1L)),
                str_sub(canonical, 2L)
            ),
            str_replace_all(canonical, "([A-Z])", "_\\1")
        )
        candidates <- str_to_lower(candidates)
        idx <- which(
            is_in(str_to_lower(methods), candidates) &
                as.character(fineMappingResult$study) == study &
                as.character(fineMappingResult$context) == context &
                as.character(fineMappingResult$trait) == trait
        )
        if (length(idx) > 0L) {
            out[[canonical]] <- getSusieFit(fineMappingResult$entry[[idx[[
                1L
            ]]]])
        }
    }
    out
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
    samplePartition <- NULL
    prediction <- list()
    performance <- list()
    for (i in idx) {
        cv <- getCvResult(fineMappingResult$entry[[i]])
        if (is.null(cv)) {
            next
        }
        # Multi-region entries store cvResult as a named per-region list; pick
        # the first region that carries a partition.
        if (is.null(cv$samplePartition)) {
            hit <- keep(cv, .twasCvHasPartition)
            if (length(hit) == 0L) {
                next
            }
            cv <- hit[[1L]]
        }
        if (is.null(samplePartition)) {
            samplePartition <- cv$samplePartition
        }
        prediction <- c(prediction, cv$prediction)
        performance <- c(performance, cv$performance)
    }
    if (length(prediction) == 0L) {
        return(NULL)
    }
    list(
        samplePartition = samplePartition,
        prediction = prediction,
        performance = performance
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
#'   (\code{"default"} or \code{"fast_default"}), or a named list of
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
#' @param mafCutoff,macCutoff,xvarCutoff,imissCutoff For QtlDataset: optional
#'   per-call genotype-filter overrides. Each replaces the corresponding
#'   construct-time \code{\link{QtlDataset}} slot for this call only (applied to
#'   a validated copy); \code{NULL} (default) leaves the stored value in place.
#'   Variant QC is a property of the data, so these are applied identically here
#'   and in \code{\link{fineMappingPipeline}} -- there is deliberately no
#'   TWAS-specific variant filter.
#' @param keepIndel,keepSamples,keepVariants As above: per-call
#'   \code{\link{QtlDataset}} genotype-filter overrides.
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
#'   this collection, the cached \code{TwasWeightsEntry} is copied through and
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
#' @param retainFit Logical. Retain the fitted-model object on each entry.
#'   Default \code{FALSE}.
#' @param retainFitDetail Character. Level of fit detail to retain:
#'   \code{"slim"} (default) or \code{"full"}.
#' @param naAction Character. How to handle missing values in the extracted
#'   data.
#' @param cvFolds Integer. Cross-validation folds. Default 5. Set to 0 to skip
#'   CV (and ensemble).
#' @param samplePartition Optional pre-defined CV partition data.frame.
#' @param maxCvVariants Maximum number of variants for CV. Default -1 (no
#'   limit).
#' @param cvThreads Threads for CV parallelism. Default 1.
#' @param cvWeightMethods Optional override of methods used for CV.
#' @param ensemble Logical. Compute SR-TWAS ensemble weights. Default
#'   \code{TRUE}.
#' @param ensembleR2Threshold Minimum CV R-squared for ensemble inclusion.
#'   Default 0.01.
#' @param ensembleSolver Solver for ensemble stacking. Default
#'   \code{"quadprog"}.
#' @param ensembleAlpha Elastic-net mixing parameter (only when
#'   \code{ensembleSolver = "glmnet"}). Default 1.
#' @param estimatePi If TRUE, estimate spike-and-slab sparsity from mr.ash
#'   before BGLR / qgg spike-and-slab methods that consume it.
#' @param phenotypeCovariatesToResidualize,genotypeCovariatesToResidualize
#'   Character vector (or \code{NULL}) of covariate column names to residualize
#'   against. Forwarded to \code{\link{getResidualizedPhenotypes}} /
#'   \code{\link{getResidualizedGenotypes}} for \code{QtlDataset} /
#'   \code{MultiStudyQtlDataset} input. Default \code{NULL} (use all available
#'   covariates). Ignored for sumstat inputs.
#' @param residualizePhenotypeCovariates Logical (length 1). When \code{TRUE}
#'   (default) residualize against the phenotype-side covariates listed in
#'   \code{phenotypeCovariatesToResidualize}; set \code{FALSE} to disable.
#' @param residualizeGenotypeCovariates Logical (length 1). When \code{TRUE}
#'   (default) residualize against the genotype-side covariates listed in
#'   \code{genotypeCovariatesToResidualize}; set \code{FALSE} to disable.
#' @param dataType Optional data-type tag stamped into every
#'   \code{TwasWeightsEntry$dataType} (e.g. \code{"expression"}).
#' @param verbose Verbosity (0 silent, 1 default, 2 includes external package
#'   messages).
#' @param seed Integer or \code{NULL}. When supplied (\code{QtlDataset} path),
#'   seeds both the main-process RNG (\code{set.seed}) and the parallel
#'   method-fitting / cross-validation RNG (\code{BiocParallel} \code{RNGseed})
#'   for reproducibility under multi-threading. \code{NULL} (default) leaves the
#'   session RNG untouched, so an outer \code{set.seed()} still governs the
#'   main-process draws.
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
    standardGeneric("twasWeightsPipeline")
})

# Run the multivariate joint TWAS-weight fit over the (context, trait) grid for
# `traits`: dispatch each cis-region through the joint engine and merge
# per-region results. `marker` = the TwasJointPipeline config; `ctx` bundles the
# shared state
# (xRegions, data, norm, useCtx, fineMappingResult, dataDrivenPriorMatricesCv,
# cisWindow, verbose).
# @noRd
.twasRunMultivariateGrid <- function(traits, marker, ctx) {
    synthSpec <- list(list(axes = c("context", "trait"), scope = NULL))
    labs <- map_chr(ctx$xRegions, .twasRegionLabel)
    perRegion <- map(
        seq_along(ctx$xRegions),
        .twasMvGridRegion,
        synthSpec = synthSpec,
        marker = marker,
        ctx = ctx,
        traits = traits
    )
    keep <- !map_lgl(perRegion, is.null)
    perRegion <- perRegion[keep]
    labs <- labs[keep]
    if (length(perRegion) == 0L) {
        return(NULL)
    }
    if (length(perRegion) == 1L) {
        return(perRegion[[1L]])
    }
    .twasMergeResultsByKey(perRegion, labs)
}

#' @rdname twasWeightsPipeline
#' @export
setMethod(
    "twasWeightsPipeline",
    "QtlDataset",
    function(
        data,
        methods = "default",
        contexts = NULL,
        traitId = NULL,
        region = NULL,
        cisWindow = NULL,
        mafCutoff = NULL,
        macCutoff = NULL,
        xvarCutoff = NULL,
        imissCutoff = NULL,
        keepIndel = NULL,
        keepSamples = NULL,
        keepVariants = NULL,
        jointRegions = FALSE,
        jointSpecification = NULL,
        fineMappingResult = NULL,
        twasWeights = NULL,
        mashPrior = NULL,
        cvFolds = 5,
        samplePartition = NULL,
        fitFullData = TRUE,
        maxCvVariants = -1,
        cvThreads = 1,
        cvWeightMethods = NULL,
        ensemble = TRUE,
        ensembleR2Threshold = 0.01,
        ensembleSolver = "quadprog",
        ensembleAlpha = 1,
        estimatePi = TRUE,
        retainFit = TRUE,
        retainFitDetail = c("slim", "full"),
        phenotypeCovariatesToResidualize = NULL,
        genotypeCovariatesToResidualize = NULL,
        residualizePhenotypeCovariates = TRUE,
        residualizeGenotypeCovariates = TRUE,
        dataType = NULL,
        naAction = c("drop", "impute"),
        verbose = 1,
        seed = NULL,
        ...
    ) {
        naAction <- arg_match(naAction)
        retainFitDetail <- arg_match(retainFitDetail)
        p <- as.list(environment())
        p$dots <- list(...)
        .twasPipelineQtlDataset(p)
    }
)

# ---- QtlDataset pipeline worker + phase helpers ----------------------------

.twasPipelineQtlDataset <- function(p) {
    .twasQdsCheckRegionCisWindow(p)
    p$xRegions <- .makeXRegions(p$region, p$jointRegions)
    # Per-call filter overrides replace the construct-time slot values on a
    # validated copy. Variant QC is a data property applied identically to
    # fine-mapping and TWAS -- there is no TWAS-specific variant filter.
    p$data <- .qtlApplyFilterOverrides(
        p$data,
        p$mafCutoff,
        p$macCutoff,
        p$xvarCutoff,
        p$imissCutoff,
        p$keepIndel,
        p$keepSamples,
        p$keepVariants
    )
    p$parsedJointSpec <- parseJointSpecification(p$jointSpecification, p$data)
    p$norm <- .twasNormalizeMethods(p$methods)
    .twasCheckMethodCapabilities(p$norm$tokens, "QtlDataset")
    .twasCheckFineMappingMethods(
        p$norm$tokens,
        p$fineMappingResult,
        "QtlDataset"
    )
    .twasQdsCheckFitFull(p)
    p <- .twasQdsUnpackMash(p)
    joint <- .twasQdsJointPhase(p)
    if (joint$done) {
        return(joint$result)
    }
    p$norm <- joint$norm
    p$jointResult <- joint$result
    p <- .twasQdsResolveGrid(p)
    .twasQdsAssemble(.twasQdsDispatch(p), p$jointResult)
}

# `cisWindow` expands a trait's own coordinates; `region` is literal. Supplying
# both signals a misunderstanding -> reject.
# @noRd
.twasQdsCheckRegionCisWindow <- function(p) {
    if (!is.null(p$region) && !is.null(p$cisWindow)) {
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
.twasQdsCheckFitFull <- function(p) {
    if (!isTRUE(p$fitFullData) && p$cvFolds <= 1L) {
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
.twasQdsUnpackMash <- function(p) {
    mp <- .unpackMashPrior(p$mashPrior, p$samplePartition)
    p$samplePartition <- mp$samplePartition
    p$dataDrivenPriorMatricesCv <- mp$dataDrivenPriorMatricesCv
    if (!is.null(p$mashPrior) && !is_in("mrmash", p$norm$tokens)) {
        msg <- glue(
            "`mashPrior` was supplied but 'mrmash' is not among `methods`; ",
            "the data-driven prior is ignored."
        )
        warn(msg)
    }
    if (
        !is.null(mp$fullPrior) &&
            is_in("mrmash_weights", names(p$norm$methodList))
    ) {
        p$norm$methodList$mrmash_weights$dataDrivenPriorMatrices <- mp$fullPrior
    }
    p
}

# Explicit jointSpecification path: run the per-spec axis dispatcher for
# mr.mash. Returns list(done, result, norm): `done` requests an early return
# with `result`; otherwise `norm` is the mrmash-stripped normalization for the
# per-tuple loop below.
# @noRd
.twasQdsJointPhase <- function(p) {
    if (length(p$parsedJointSpec) == 0L) {
        return(list(done = FALSE, result = NULL, norm = p$norm))
    }
    jointResult <- .twasDispatchJointSpecsQtlDataset(
        p$parsedJointSpec,
        p$data,
        intersect(p$norm$tokens, "mrmash"),
        p$contexts,
        p$traitId,
        p$cisWindow,
        p$dataType,
        p$verbose,
        xRegions = p$xRegions,
        retainFit = p$retainFit,
        retainFitDetail = p$retainFitDetail,
        seed = p$seed
    )
    keep <- setdiff(p$norm$tokens, intersect(p$norm$tokens, "mrmash"))
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
    norm <- p$norm
    norm$tokens <- keep
    keepKeys <- which(
        is_in(str_remove(names(norm$methodList), "(_weights|Weights)$"), keep)
    )
    norm$methodList <- norm$methodList[keepKeys]
    list(done = FALSE, result = jointResult, norm = norm)
}

# Resolve the (context, trait) grid + multivariate flag and build the joint-
# pipeline marker + shared grid context. Returns the updated parameter bundle.
# @noRd
.twasQdsResolveGrid <- function(p) {
    p$study <- getStudy(p$data)
    p$useCtx <- .twasQdsResolveContexts(p$data, p$contexts)
    p$allTraits <- .twasQdsResolveTraits(
        p$data,
        p$useCtx,
        p$traitId,
        p$region
    )
    p$nCtx <- length(p$useCtx)
    .twasCheckMultivariateY(p$norm$tokens, length(p$allTraits), p$nCtx)
    p$multivariate <- any(map_lgl(p$norm$tokens, .twasIsMultivariateToken))
    p$marker <- .twasQdsMarker(p)
    p$twasGridCtx <- .twasQdsGridCtx(p)
    p
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
    allTraits <- unique(unlist(perCtxTraits))
    if (length(allTraits) == 0L) {
        abort("twasWeightsPipeline(QtlDataset): no traits selected.")
    }
    allTraits
}

# Joint-pipeline marker carrying the CV / ensemble config for the engine.
# @noRd
.twasQdsMarker <- function(p) {
    new(
        "TwasJointPipeline",
        config = list(
            cvFolds = p$cvFolds,
            samplePartition = p$samplePartition,
            fitFullData = p$fitFullData,
            dataType = p$dataType,
            retainFitDetail = p$retainFitDetail,
            standardized = FALSE,
            ensemble = p$ensemble,
            ensembleR2Threshold = p$ensembleR2Threshold,
            ensembleSolver = p$ensembleSolver,
            ensembleAlpha = p$ensembleAlpha,
            maxCvVariants = p$maxCvVariants,
            cvThreads = p$cvThreads,
            estimatePi = p$estimatePi,
            verbose = p$verbose,
            seed = p$seed,
            ldSketch = NULL
        )
    )
}

# Shared grid context consumed by both dispatch paths.
# @noRd
.twasQdsGridCtx <- function(p) {
    list(
        xRegions = p$xRegions,
        data = p$data,
        norm = p$norm,
        useCtx = p$useCtx,
        fineMappingResult = p$fineMappingResult,
        dataDrivenPriorMatricesCv = p$dataDrivenPriorMatricesCv,
        cisWindow = p$cisWindow,
        verbose = p$verbose
    )
}

# Top-level dispatch: multivariate joint grid vs univariate engine path.
# @noRd
.twasQdsDispatch <- function(p) {
    if (p$multivariate) {
        return(.twasRunMultivariateGrid(p$allTraits, p$marker, p$twasGridCtx))
    }
    .twasQdsUnivariateEngine(p)
}

# Univariate methods ROUTED THROUGH THE ENGINE: one 1-condition group per
# (context, trait), per region -> the SAME per-method fitter (+ ensemble layer
# for >= 2 methods + resume cache) as the joint paths, merged across regions.
# @noRd
.twasQdsUnivariateEngine <- function(p) {
    univCell <- .lookupJointCell("univariate", "individual")
    scope <- list(
        studies = p$study,
        contexts = set_names(list(p$useCtx), p$study),
        traits = set_names(list(p$allTraits), p$study)
    )
    labs <- map_chr(p$xRegions, .twasRegionLabel)
    perRegion <- map(
        seq_along(p$xRegions),
        .twasQdsUnivRegion,
        univCell = univCell,
        p = p,
        scope = scope
    )
    keep <- !map_lgl(perRegion, is.null)
    .twasMergeRegionResults(perRegion[keep], labs[keep])
}

# Per-region engine args for the univariate path.
# @noRd
.twasQdsUnivArgs <- function(p, bi) {
    list(
        methodList = p$norm$methodList,
        fineMappingResult = p$fineMappingResult,
        cache = p$twasWeights,
        dataDrivenPriorMatricesCv = p$dataDrivenPriorMatricesCv,
        cisWindow = p$cisWindow,
        region = p$xRegions[[bi]],
        regionIndex = bi,
        nRegions = length(p$xRegions),
        naAction = p$naAction,
        verbose = p$verbose
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
        retainFit = TRUE,
        retainFitDetail = c("slim", "full"),
        dataType = NULL,
        verbose = 1L,
        ...
    ) {
        retainFitDetail <- arg_match(retainFitDetail)
        p <- as.list(environment())
        p$dots <- list(...)
        .twasPipelineQtlSumStats(p)
    }
)

# ---- QtlSumStats pipeline worker + phase helpers ---------------------------

.twasPipelineQtlSumStats <- function(p) {
    # summaryStatsQc() is mandatory before twasWeightsPipeline for SumStats
    # input; it also drops variants not present in the ldSketch, so every
    # entry's SNP set is a subset of the ldSketch panel by the time we get here.
    .twasAssertQcd(p$data)
    p$parsedJointSpec <- parseJointSpecification(p$jointSpecification, p$data)
    tm <- .twasSumStatsMethodTokens(p$methods)
    p$tokens <- tm$tokens
    p$methodArgs <- tm$methodArgs
    .twasCheckMethodCapabilities(p$tokens, "QtlSumStats")
    .twasCheckFineMappingMethods(p$tokens, p$fineMappingResult, "QtlSumStats")
    joint <- .twasQssJointPhase(p)
    if (joint$done) {
        return(joint$result)
    }
    p$tokens <- joint$tokens
    p$methodArgs <- joint$methodArgs
    p <- .twasQssSelectAndPartition(p)
    rows <- c(.twasQssUnivariateRows(p), .twasQssMultivariateRows(p))
    .twasQssAssemble(rows, joint$result, p)
}

# Normalize the methods argument into (tokens, methodArgs). The default set
# excludes fine-mapping methods; those must be requested explicitly together
# with a FineMappingResult passed via `fineMappingResult`.
# @noRd
.twasSumStatsMethodTokens <- function(methods) {
    if (is.null(methods)) {
        tokens <- c("lasso", "prsCs", "dpr_gibbs")
        return(list(tokens = tokens, methodArgs = .twasEmptyMethodArgs(tokens)))
    }
    if (is.character(methods)) {
        return(list(
            tokens = methods,
            methodArgs = .twasEmptyMethodArgs(methods)
        ))
    }
    if (is.list(methods)) {
        return(list(tokens = names(methods), methodArgs = methods))
    }
    msg <- glue(
        "`methods` must be NULL, a character vector, or a named list ",
        "of <token> = <args> entries."
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
.twasQssJointPhase <- function(p) {
    if (length(p$parsedJointSpec) == 0L) {
        return(list(
            done = FALSE,
            result = NULL,
            tokens = p$tokens,
            methodArgs = p$methodArgs
        ))
    }
    jointResult <- .twasDispatchJointSpecsQtlSumStats(
        p$parsedJointSpec,
        p$data,
        intersect(p$tokens, "mrmash"),
        p$contexts,
        p$traitId,
        p$dataType,
        p$verbose,
        retainFit = p$retainFit,
        retainFitDetail = p$retainFitDetail
    )
    keep <- setdiff(p$tokens, "mrmash")
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
        methodArgs = p$methodArgs[keep]
    )
}

# Resolve the selected rows, partition tokens into univariate vs multivariate,
# attach the LD sketch, and enforce the multivariate >=2-contexts rule. Returns
# the updated parameter bundle.
# @noRd
.twasQssSelectAndPartition <- function(p) {
    p$studyCol <- as.character(p$data$study)
    p$contextCol <- as.character(p$data$context)
    p$traitCol <- as.character(p$data$trait)
    p$selRows <- .twasQssSelectRows(p)
    isMv <- map_lgl(p$tokens, .twasIsMultivariateToken)
    p$multivariateTokens <- p$tokens[isMv]
    p$univariateTokens <- p$tokens[!isMv]
    p$ldSketch <- getLdSketch(p$data)
    .twasQssCheckMultivariate(p)
    p
}

# Row indices matching the contexts / traitId filters (error if none).
# @noRd
.twasQssSelectRows <- function(p) {
    selRows <- seq_len(nrow(p$data))
    if (!is.null(p$contexts)) {
        selRows <- selRows[is_in(p$contextCol[selRows], p$contexts)]
    }
    if (!is.null(p$traitId)) {
        selRows <- selRows[is_in(p$traitCol[selRows], p$traitId)]
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
.twasQssCheckMultivariate <- function(p) {
    if (length(p$multivariateTokens) == 0L) {
        return(invisible(NULL))
    }
    groupKey <- str_c(
        p$studyCol[p$selRows],
        p$traitCol[p$selRows],
        sep = "||"
    )
    perGroupNCtx <- map_int(split(p$contextCol[p$selRows], groupKey), length)
    if (all(perGroupNCtx < 2L)) {
        mvStr <- str_flatten(p$multivariateTokens, ", ")
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
.twasMrmashRetainDefaults <- function(userArgs, adapter, tk, retainFitDetail) {
    if (!is.null(adapter) || tk != "mrmash") {
        return(userArgs)
    }
    if (is.null(userArgs$retainFit)) {
        userArgs$retainFit <- TRUE
    }
    if (is.null(userArgs$fitDetail)) {
        userArgs$fitDetail <- retainFitDetail
    }
    userArgs
}

# Run a weight function, warning (with `errPrefix`) and returning NULL on error.
# @noRd
.twasTryWeights <- function(fn, stat, ldMat, userArgs, errPrefix) {
    tryCatch(
        {
            wfn <- get(fn, mode = "function")
            wArgs <- c(list(stat = stat, LD = ldMat), userArgs)
            exec(wfn, !!!wArgs)
        },
        error = function(e) {
            eMsg <- conditionMessage(e)
            msg <- glue("{errPrefix}{eMsg}")
            warn(msg)
            NULL
        }
    )
}

# ---- Univariate dispatch: per (study, context, trait), per method ----------

.twasQssUnivariateRows <- function(p) {
    if (length(p$univariateTokens) == 0L) {
        return(list())
    }
    list_flatten(map(p$selRows, .twasQssUnivariateRowsForEntry, p = p))
}

# Cached + freshly-fitted rows for one sumstats entry. Resume cache: pull cached
# entries up front and reduce the per-entry fit work to the un-cached tokens.
# @noRd
.twasQssUnivariateRowsForEntry <- function(i, p) {
    st <- p$studyCol[i]
    ctx <- p$contextCol[i]
    tr <- p$traitCol[i]
    cacheHits <- .twasResolveCacheHits(
        p$twasWeights,
        st,
        ctx,
        tr,
        p$univariateTokens
    )
    cachedRows <- imap(
        cacheHits,
        .twasCachedRowRecord,
        st = st,
        ctx = ctx,
        tr = tr
    )
    toFit <- setdiff(p$univariateTokens, names(cacheHits))
    if (length(toFit) == 0L) {
        return(unname(cachedRows))
    }
    fitCtx <- .twasQssUnivariateFitCtx(p$data, st, ctx, tr, p$ldSketch)
    fitted <- compact(map(
        toFit,
        .twasQssUnivariateFitOne,
        st = st,
        ctx = ctx,
        tr = tr,
        fitCtx = fitCtx,
        p = p
    ))
    c(unname(cachedRows), fitted)
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
.twasQssUnivariateFitCtx <- function(data, st, ctx, tr, ldSketch) {
    df <- getSumstatDf(
        data,
        study = st,
        context = ctx,
        trait = tr,
        require = c("Z", "N"),
        derive = "zFromBetaSe"
    )
    variantIds <- df$variant_id
    varY <- getVarY(data, study = st, context = ctx, trait = tr)
    if (is.null(varY)) {
        varY <- 1
    }
    stat <- list(
        z = df$z,
        n = stats::median(df$N, na.rm = TRUE),
        varY = varY,
        variantNames = variantIds
    )
    list(
        variantIds = variantIds,
        stat = stat,
        ldMat = .twasLdFromSketch(ldSketch, variantIds)
    )
}

# Fit one univariate method for one entry -> a row record, or NULL on skip.
# @noRd
.twasQssUnivariateFitOne <- function(tk, st, ctx, tr, fitCtx, p) {
    spec <- .twasResolveWeightFn(tk)
    userArgs <- .twasUserArgs(p$methodArgs, tk)
    # When the token is a fine-mapping method, pass the precomputed fit into the
    # *Rss weight function via its dedicated *Fit arg. The gate above ensures
    # fineMappingResult is non-NULL here.
    if (!is.null(spec$adapter)) {
        fit <- .twasFineMappingFitFor(
            p$fineMappingResult,
            study = st,
            context = ctx,
            trait = tr,
            token = tk
        )
        if (is.null(fit)) {
            .twasWarnNoFitUniv(tk, st, ctx, tr)
            return(NULL)
        }
        userArgs[[spec$adapter$rssFitArg]] <- fit
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
    fitAttr <- attr(weights, "fit")
    attr(weights, "fit") <- NULL
    .twasRowRecord(
        st,
        ctx,
        tr,
        tk,
        TwasWeightsEntry(
            variantIds = fitCtx$variantIds,
            weights = as.numeric(weights),
            fits = fitAttr,
            cvResult = NULL,
            standardized = TRUE,
            dataType = p$dataType
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

.twasQssMultivariateRows <- function(p) {
    if (length(p$multivariateTokens) == 0L) {
        return(list())
    }
    groupKey <- str_c(
        p$studyCol[p$selRows],
        p$traitCol[p$selRows],
        sep = "||"
    )
    groups <- split(p$selRows, groupKey)
    list_flatten(map(groups, .twasQssMultivariateGroupRows, p = p))
}

# Multivariate rows for one (study, trait) group across its contexts.
# @noRd
.twasQssMultivariateGroupRows <- function(gIdx, p) {
    if (length(gIdx) < 2L) {
        return(list())
    }
    st <- p$studyCol[gIdx[[1L]]]
    tr <- p$traitCol[gIdx[[1L]]]
    ctxNames <- p$contextCol[gIdx]
    mvStat <- .twasQssMultivariateStat(p$data, st, tr, ctxNames)
    ldMat <- .twasLdFromSketch(p$ldSketch, mvStat$variantIds)
    list_flatten(map(
        p$multivariateTokens,
        .twasQssMultivariateFitOne,
        st = st,
        tr = tr,
        ctxNames = ctxNames,
        mvStat = mvStat,
        ldMat = ldMat,
        p = p
    ))
}

# Build the (variants x contexts) Z matrix + per-context N for a group. All
# entries in a (study, trait) group must share an identical variant order after
# summaryStatsQc().
# @noRd
.twasQssMultivariateStat <- function(data, st, tr, ctxNames) {
    firstDf <- getSumstatDf(
        data,
        study = st,
        context = ctxNames[[1L]],
        trait = tr,
        require = c("Z", "N"),
        derive = "zFromBetaSe"
    )
    variantIds <- firstDf$variant_id
    Z <- matrix(
        NA_real_,
        nrow = length(variantIds),
        ncol = length(ctxNames),
        dimnames = list(variantIds, ctxNames)
    )
    nVec <- numeric(length(ctxNames))
    for (kk in seq_along(ctxNames)) {
        d <- getSumstatDf(
            data,
            study = st,
            context = ctxNames[[kk]],
            trait = tr,
            require = c("Z", "N"),
            derive = "zFromBetaSe"
        )
        if (!identical(d$variant_id, variantIds)) {
            msg <- glue(
                "twasWeightsPipeline(QtlSumStats, multivariate): every ",
                "entry for (study='{st}', trait='{tr}') must share an ",
                "identical SNP order after summaryStatsQc(). Use the same ",
                "ldSketch on every entry."
            )
            abort(msg)
        }
        Z[, kk] <- d$z
        nVec[kk] <- stats::median(d$N, na.rm = TRUE)
    }
    names(nVec) <- ctxNames
    list(
        variantIds = variantIds,
        stat = list(z = Z, n = nVec, variantNames = variantIds)
    )
}

# Fit one multivariate method for a group -> one row record per context (empty
# list on skip).
# @noRd
.twasQssMultivariateFitOne <- function(tk, st, tr, ctxNames, mvStat, ldMat, p) {
    spec <- .twasResolveWeightFn(tk)
    userArgs <- .twasMrmashRetainDefaults(
        .twasUserArgs(p$methodArgs, tk),
        spec$adapter,
        tk,
        p$retainFitDetail
    )
    # mvsusie is fine-mapping; thread its pre-fit through (mr.mash is not).
    if (!is.null(spec$adapter)) {
        userArgs <- .twasMvThreadFit(spec, userArgs, tk, st, tr, ctxNames, p)
        if (is.null(userArgs)) {
            return(list())
        }
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
    if (!is.matrix(weights)) {
        weights <- as.matrix(weights)
    }
    fitAttr <- attr(weights, "fit")
    attr(weights, "fit") <- NULL
    .twasMvContextRows(
        weights,
        fitAttr,
        ctxNames,
        mvStat,
        st,
        tr,
        tk,
        p$dataType
    )
}

# Thread the precomputed fine-mapping fit into a multivariate method's args;
# returns NULL (signalling skip) when the fit is absent.
# @noRd
.twasMvThreadFit <- function(spec, userArgs, tk, st, tr, ctxNames, p) {
    fit <- .twasFineMappingFitFor(
        p$fineMappingResult,
        study = st,
        context = ctxNames[[1L]],
        trait = tr,
        token = tk
    )
    if (is.null(fit)) {
        .twasWarnNoFitMv(tk, st, tr)
        return(NULL)
    }
    userArgs[[spec$adapter$rssFitArg]] <- fit
    userArgs
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
.twasQssAssemble <- function(rows, jointResult, p) {
    perTupleResult <- .twasRowsToWeights(rows, p$ldSketch)
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
    .rbindTwasWeights(perTupleResult, jointResult, ldSketch = p$ldSketch)
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
    twArgs <- c(
        list(
            data = qd,
            methods = cfg$methods,
            contexts = cfg$contexts,
            traitId = cfg$traitId,
            region = cfg$region,
            cisWindow = cfg$cisWindow,
            jointRegions = cfg$jointRegions,
            jointSpecification = NULL,
            fineMappingResult = cfg$fineMappingResult,
            twasWeights = cfg$twasWeights,
            naAction = cfg$naAction,
            verbose = cfg$verbose,
            seed = cfg$seed
        ),
        cfg$dotArgs
    )
    exec(twasWeightsPipeline, !!!twArgs)
}

# Embedded-sumstats TWAS-weights worker for .multiStudyPipelineDriver.
# @noRd
.twasSumStats <- function(ss, cfg) {
    twArgs <- c(
        list(
            data = ss,
            methods = cfg$methods,
            contexts = cfg$contexts,
            traitId = cfg$traitId,
            jointSpecification = NULL,
            fineMappingResult = cfg$fineMappingResult,
            twasWeights = cfg$twasWeights,
            verbose = cfg$verbose
        ),
        cfg$dotArgs
    )
    exec(twasWeightsPipeline, !!!twArgs)
}

#' @rdname twasWeightsPipeline
#' @export
setMethod(
    "twasWeightsPipeline",
    "MultiStudyQtlDataset",
    function(
        data,
        methods = "default",
        contexts = NULL,
        traitId = NULL,
        region = NULL,
        cisWindow = NULL,
        jointRegions = FALSE,
        jointSpecification = NULL,
        fineMappingResult = NULL,
        twasWeights = NULL,
        retainFit = TRUE,
        retainFitDetail = c("slim", "full"),
        naAction = c("drop", "impute"),
        verbose = 1,
        phenotypeCovariatesToResidualize = NULL,
        genotypeCovariatesToResidualize = NULL,
        residualizePhenotypeCovariates = TRUE,
        residualizeGenotypeCovariates = TRUE,
        seed = NULL,
        ...
    ) {
        naAction <- arg_match(naAction)
        retainFitDetail <- arg_match(retainFitDetail)
        p <- as.list(environment())
        p$dots <- list(...)
        .twasPipelineMultiStudy(p)
    }
)

# ---- MultiStudyQtlDataset pipeline worker + phase helpers ------------------

# Method tokens from a `methods` arg: a character vector as-is; a named list ->
# canonical bare tokens; otherwise empty.
# @noRd
.twasMethodTokensFromArg <- function(methods) {
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
.twasMsStripMrmash <- function(methods) {
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
    (is.character(methods) || is.list(methods)) && length(methods) == 0L
}

.twasPipelineMultiStudy <- function(p) {
    if (!is.null(p$region) && !is.null(p$cisWindow)) {
        msg <- glue(
            "twasWeightsPipeline(MultiStudyQtlDataset): specify either ",
            "`region` or `cisWindow`, not both."
        )
        abort(msg)
    }
    xRegions <- .makeXRegions(p$region, p$jointRegions)
    parsedJointSpec <- parseJointSpecification(p$jointSpecification, p$data)
    # Gate fine-mapping methods early so the recursion into the embedded
    # QtlDataset / QtlSumStats components doesn't re-run fine-mapping.
    .twasCheckFineMappingMethods(
        .twasMethodTokensFromArg(p$methods),
        p$fineMappingResult,
        "MultiStudyQtlDataset"
    )
    joint <- .twasMsJointPhase(p, parsedJointSpec, xRegions)
    if (joint$done) {
        return(joint$result)
    }
    .twasMsDriver(p, joint$result, joint$methods)
}

# Joint-specification dispatch for mrmash. Returns list(done, result, methods)
# where `done` requests an early return with `result` and `methods` is the
# mrmash-stripped set for the per-component recursion.
# @noRd
.twasMsJointPhase <- function(p, parsedJointSpec, xRegions) {
    if (length(parsedJointSpec) == 0L) {
        return(list(done = FALSE, result = NULL, methods = p$methods))
    }
    jointMethods <- intersect(.twasMethodTokensFromArg(p$methods), "mrmash")
    jointResult <- .twasDispatchJointSpecsMultiStudy(
        parsedJointSpec,
        p$data,
        jointMethods,
        p$contexts,
        p$traitId,
        p$cisWindow,
        NULL,
        p$verbose,
        xRegions = xRegions,
        retainFit = p$retainFit,
        retainFitDetail = p$retainFitDetail,
        seed = p$seed
    )
    methods <- .twasMsStripMrmash(p$methods)
    if (.twasMethodsEmpty(methods)) {
        if (is.null(jointResult)) {
            msg <- glue(
                "twasWeightsPipeline(MultiStudyQtlDataset): no joint fits ",
                "produced."
            )
            abort(msg)
        }
        return(list(done = TRUE, result = jointResult))
    }
    list(done = FALSE, result = jointResult, methods = methods)
}

# Run the per-study / per-component recursion via the shared multi-study driver.
# @noRd
.twasMsDriver <- function(p, jointResult, methods) {
    cfg <- list(
        methods = methods,
        contexts = p$contexts,
        traitId = p$traitId,
        region = p$region,
        cisWindow = p$cisWindow,
        jointRegions = p$jointRegions,
        fineMappingResult = p$fineMappingResult,
        twasWeights = p$twasWeights,
        naAction = p$naAction,
        verbose = p$verbose,
        seed = p$seed,
        dotArgs = p$dots
    )
    .multiStudyPipelineDriver(
        p$data,
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


#' @rdname twasWeightsPipeline
#' @export
setMethod("twasWeightsPipeline", "ANY", function(data, ...) {
    cls <- class(data)[[1L]]
    msg <- glue(
        "twasWeightsPipeline does not accept inputs of class '{cls}'. ",
        "Pass a QtlDataset, MultiStudyQtlDataset, or QtlSumStats. ",
        "(GwasSumStats inputs are not supported; GWAS-side per-LD-block ",
        "weights are produced inside the new ctwasPipeline / ",
        "qtlEnrichmentPipeline.)"
    )
    abort(msg)
})

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
        # nocov start
        msg <- glue(
            "Package 'quadprog' is required for solver='quadprog'. ",
            "Install with: install.packages('quadprog')"
        )
        abort(msg)
        # nocov end
    }

    Dmat <- crossprod(Pvalid)
    dvec <- as.vector(crossprod(Pvalid, yObs))
    # Ridge term for numerical stability (small relative to trace)
    Dmat <- Dmat + 1e-8 * mean(diag(Dmat)) * diag(Kvalid)

    # Constraint matrix: first constraint is equality (sum = 1), then Kvalid
    # non-negativity constraints.
    Amat <- cbind(rep(1, Kvalid), diag(Kvalid))
    bvec <- c(1, rep(0, Kvalid))

    qpSol <- tryCatch(
        solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1),
        error = function(e) {
            eMsg <- conditionMessage(e)
            msg <- glue(
                "QP solver failed: {eMsg}. Falling back to equal weights ",
                "among valid methods."
            )
            warn(msg)
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
    if (!requireNamespace("nnls", quietly = TRUE)) {
        # nocov start
        msg <- glue(
            "Package 'nnls' is required for solver='nnls'. ",
            "Install with: install.packages('nnls')"
        )
        abort(msg)
        # nocov end
    }

    fit <- tryCatch(
        nnls::nnls(Pvalid, yObs),
        error = function(e) {
            eMsg <- conditionMessage(e)
            msg <- glue(
                "NNLS solver failed: {eMsg}. Falling back to equal weights."
            )
            warn(msg)
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
.ensembleObj <- function(z, Pvalid, yObs, ...) sum((yObs - Pvalid %*% z)^2)

# Gradient of the ensemble stacking objective. `...` absorbs the objective's
# extra optim args (Pvalid, yObs).
# @noRd
.ensembleGrad <- function(z, PtP, Pty, ...) as.vector(2 * (PtP %*% z - Pty))

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

    fit <- tryCatch(
        optim(
            par = rep(1 / Kvalid, Kvalid),
            fn = .ensembleObj,
            gr = .ensembleGrad,
            Pvalid = Pvalid,
            yObs = yObs,
            PtP = PtP,
            Pty = Pty,
            method = "L-BFGS-B",
            lower = rep(0, Kvalid)
        ),
        error = function(e) {
            eMsg <- conditionMessage(e)
            msg <- glue(
                "L-BFGS-B solver failed: {eMsg}. Falling back to equal ",
                "weights."
            )
            warn(msg)
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
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        # nocov start
        msg <- glue(
            "Package 'glmnet' is required for solver='glmnet'. ",
            "Install with: install.packages('glmnet')"
        )
        abort(msg)
        # nocov end
    }

    fit <- tryCatch(
        glmnet::cv.glmnet(
            x = Pvalid,
            y = yObs,
            lower.limits = 0,
            alpha = alpha,
            intercept = FALSE
        ),
        error = function(e) {
            eMsg <- conditionMessage(e)
            msg <- glue(
                "glmnet solver failed: {eMsg}. Falling back to equal ",
                "weights."
            )
            warn(msg)
            NULL
        }
    )

    if (is.null(fit)) {
        return(rep(1 / Kvalid, Kvalid))
    }

    zetaValid <- as.numeric(coef(fit, s = "lambda.min"))[-1] # drop intercept
    zetaValid <- pmax(zetaValid, 0)
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
.ensembleValidateArgs <- function(cvResults, Y, contextIndex) {
    if (is.null(cvResults)) {
        abort("'cvResults' is required.")
    }
    if (is.null(Y)) {
        abort("'Y' is required.")
    }
    if (
        !is.numeric(contextIndex) ||
            length(contextIndex) != 1 ||
            contextIndex < 1
    ) {
        abort("'contextIndex' must be a positive integer scalar.")
    }
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
.ensembleBuildPd <- function(predsD, nm, aln, contextIndex, d) {
    Pd <- matrix(NA_real_, nrow = aln$nD, ncol = nm$K)
    colnames(Pd) <- nm$baseNames
    for (k in seq_along(nm$predNames)) {
        predMat <- predsD[[nm$predNames[k]]]
        pCol <- if (is.matrix(predMat)) {
            predMat[aln$predOrder, contextIndex]
        } else {
            as.numeric(predMat)[aln$predOrder]
        }
        if (length(pCol) != aln$nD) {
            methodName <- nm$predNames[k]
            nCol <- length(pCol)
            nAligned <- aln$nD
            msg <- glue(
                "Prediction length for method '{methodName}' in dataset ",
                "{d} ({nCol}) does not match number of aligned samples ",
                "({nAligned})."
            )
            abort(msg)
        }
        Pd[, k] <- pCol
    }
    Pd
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
    zeta <- rep(0, nm$K)
    zeta[validMethods] <- zetaValid
    names(zeta) <- nm$baseNames
    zeta
}

# Degenerate case: a single signal-bearing method takes full weight.
# @noRd
.ensembleSingleMethodZeta <- function(validMethods, baseNames, K) {
    zeta <- rep(0, K)
    zeta[validMethods] <- 1
    names(zeta) <- baseNames
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
.ensembleAccumulateWeights <- function(wtList, wtKeys, matched, zeta) {
    firstWt <- .ensembleAsMatrix(wtList[[wtKeys[which(matched)[1]]]])
    ensembleTwasWt <- matrix(0, nrow = nrow(firstWt), ncol = ncol(firstWt))
    rownames(ensembleTwasWt) <- rownames(firstWt)
    colnames(ensembleTwasWt) <- colnames(firstWt)
    for (i in which(matched)) {
        wMat <- .ensembleAsMatrix(wtList[[wtKeys[i]]])
        if (!identical(dim(wMat), dim(ensembleTwasWt))) {
            wtKey <- wtKeys[i]
            msg <- glue(
                "Weight matrix for '{wtKey}' has inconsistent dimensions; ",
                "skipping."
            )
            warn(msg)
            next
        }
        ensembleTwasWt <- ensembleTwasWt + zeta[i] * wMat
    }
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
    cv <- getCvResult(e)
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
    if (length(hit)) tw$entry[[hit[[1L]]]] else NULL
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

# One region's univariate joint-cell fit (region `bi` of p$xRegions).
# @noRd
.twasQdsUnivRegion <- function(bi, univCell, p, scope) {
    .runJointCell(
        univCell,
        p$marker,
        p$data,
        scope,
        p$norm$tokens,
        args = .twasQdsUnivArgs(p, bi)
    )
}

# One cached row record from an imap over (entry, token) cache hits.
# @noRd
.twasCachedRowRecord <- function(entry, tk, st, ctx, tr) {
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
        TwasWeightsEntry(
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
