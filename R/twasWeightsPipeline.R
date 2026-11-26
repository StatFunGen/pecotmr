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
        weightStandardized = weightStandardized(payloads[[1L]])
    )
}

# Unpack a MashPrior input into the internal twasWeightsPipeline arguments:
# $priors -> list(full, folds): the full-data priors and each fold's, from
# which each mr.mash group takes its own (.jointTwasGroupPriors())
# $samplePartition -> the CV folds (an explicit `samplePartition` arg wins;
# otherwise the partition the per-fold priors were computed on). NULL input
# returns no priors, preserving the supplied samplePartition.
# @noRd
#' @importFrom checkmate assertClass
.unpackMashPrior <- function(mashPrior, samplePartition = NULL) {
    if (is.null(mashPrior)) {
        return(list(priors = NULL, samplePartition = samplePartition))
    }
    assertClass(mashPrior, "MashPrior")
    cvFits <- cvFits(mashPrior)
    list(
        priors = list(full = fullFit(mashPrior), folds = cvFits$perFoldFits),
        samplePartition = samplePartition %||% cvFits$samplePartition
    )
}

# The data-driven prior one mr.mash group takes from a set of priors (a
# MashPrior's full-data set, or one fold's): the prior over the axes of the
# group's joint spec, `axes`, whose conditions include each of the group's
# `cond` and whose scope, on the axes it pools, admits their values. Its
# covariances are taken on the group's conditions, in their order, named by
# `labels`, and put on the per-SD scale mr.mash fits on
# (.mashPriorScaleFactor()); its residual correlation, when it carries one,
# on the same conditions unscaled; and its mixture weights (`w`, by
# component), when it carries them. No such prior, or several, is an error.
# @noRd
.mashPriorForGroup <- function(priors, cond, axes, labels) {
    rows <- map(priors, .mashPriorRows, cond = cond, axes = axes)
    hit <- which(!map_lgl(rows, is.null))
    if (length(hit) != 1L) {
        .mashPriorMatchError(names(priors), length(hit), cond, axes)
    }
    prior <- priors[[hit]]
    residual <- prior$residualCorrelation
    compact(list(
        U = map(
            prior$U,
            .mashPriorOnGroup,
            rows = rows[[hit]],
            scale = .mashPriorScaleFactor(prior, rows[[hit]]),
            labels = labels
        ),
        residualCorrelation = if (!is.null(residual)) {
            .mashPriorOnGroup(residual, rows[[hit]], 1, labels)
        },
        w = prior$w
    ))
}

# The rows of a prior's covariances that are the group's conditions, in
# their order, or NULL when the prior is over other axes, leaves one of
# them out, or is scoped away from them.
# @noRd
#' @importFrom vctrs vec_match
.mashPriorRows <- function(prior, cond, axes) {
    if (
        !setequal(prior$axes, axes) ||
            !.mashPriorScopeAdmits(prior$scope, cond, axes)
    ) {
        return(NULL)
    }
    rows <- vec_match(
        .mashConditionKeys(cond, axes),
        .mashConditionKeys(prior$conditions, axes)
    )
    if (anyNA(rows)) NULL else rows
}

# The conditions' values on `axes`, as text, one row per condition.
# @noRd
.mashConditionKeys <- function(conditions, axes) {
    mutate(as_tibble(conditions[axes]), across(everything(), as.character))
}

# Whether a prior's scope admits the group: on each axis the prior pools, the
# scope lists every condition's value.
# @noRd
.mashPriorScopeAdmits <- function(scope, cond, axes) {
    pooled <- setdiff(intersect(names(scope), .spValidJointAxes), axes)
    all(map_lgl(pooled, .mashScopeAdmitsAxis, scope = scope, cond = cond))
}

# @noRd
.mashScopeAdmitsAxis <- function(axis, scope, cond) {
    all(is_in(as.character(cond[[axis]]), scope[[axis]]))
}

# What divides a prior's covariances to put them on the per-SD scale mr.mash
# fits on: (n_k n_l)^(e / 2), for the conditions' sample sizes n. mash fits
# Bhat / Shat^alpha (mash_set_data()), which is the per-SD effect times
# n^(alpha / 2) on the standardized scale and the z-score on the z scale,
# where Shat is 1; e is 1 there and `alpha` otherwise.
# @noRd
.mashPriorScaleFactor <- function(prior, rows) {
    e <- if (identical(prior$inputScale, "z")) 1 else prior$alpha
    if (e == 0) {
        return(1)
    }
    n <- prior$conditions$n[rows]
    if (anyNA(n)) {
        abort(glue(
            "mashPrior: the prior over {str_flatten(prior$axes, '+')} was ",
            "fitted on the scale of Bhat / Shat^{e}, which is put on ",
            "mr.mash's per-SD scale by each condition's sample size, and its ",
            "summary statistics carried no N for some of the group's ",
            "conditions."
        ))
    }
    tcrossprod(n^(e / 2))
}

# One covariance matrix on the group's conditions, rescaled and named.
# @noRd
.mashPriorOnGroup <- function(U, rows, scale, labels) {
    `dimnames<-`(U[rows, rows, drop = FALSE] / scale, list(labels, labels))
}

# @noRd
.mashPriorMatchError <- function(priorNames, nHit, cond, axes) {
    what <- if (nHit == 0L) "none" else "more than one"
    conditions <- exec(str_c, !!!unname(as.list(cond[axes])), sep = ":")
    abort(glue(
        "mashPrior: {what} of its priors ",
        "({str_flatten(priorNames, ', ')}) is over the axes ",
        "{str_flatten(axes, '+')} and the conditions ",
        "{str_flatten(conditions, ', ')}; mashPipeline() fits one per ",
        "joint spec (`jointSpecification`)."
    ))
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
#   coupled         How the method is fitted across a group's conditions
#                   under the penalty on the difference between their effects
#                   along the axes `coupling` names, by path: on
#                   individual-level data the weight function, on whichever
#                   engine the run's backend for the method is
#                   (.penalizedCoupledWeights()); on summary statistics the
#                   `penalty` the penalized RSS solver takes
#                   (.rssCoupledWeights()), or the `sampler` fitting the
#                   conditions at once. Absent where it has none.
#   grouped         The weight function fitting a group's conditions at once
#                   on individual-level data, each variant's effects across
#                   them one group under the penalty: glmnet's mgaussian
#                   family, grpreg's group penalties or abess's best subset,
#                   whichever the run's backend for the method is, or for
#                   ridge each condition's own fit at one lambda
#                   (.twasGroupedLayouts()). Absent where it has none.
#
# The Bayesian alphabet (bayesA/B/C/L/N/R) runs on several engines per data
# path, and so do the lasso and LDpred2 on summary statistics, and on
# individual-level data the lasso, elastic net, MCP and SCAD (glmnet or
# ncvreg, and grpreg) and L0 (L0Learn and abess); their
# implementations take `backend`, and .methodBackends records which engines
# and Options constructors each method has. dpr_gibbs has the
# SDPR sumstat counterpart;
# dpr_vb / dpr_adaptive_gibbs remain individual-only, as does sparsevb, whose
# engine has no summary-statistics form. lasso, enet and ridge
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
        multivariate = FALSE,
        coupled = list(
            individual = ".lassoCoupledWeights",
            sumstat = list(penalty = "lasso")
        ),
        grouped = ".lassoGroupedWeights"
    ),
    scad = list(
        individualImpl = "scadWeights",
        sumstatImpl = "scadRssWeights",
        multivariate = FALSE,
        coupled = list(
            individual = ".scadCoupledWeights",
            sumstat = list(penalty = "SCAD")
        ),
        grouped = ".scadGroupedWeights"
    ),
    mcp = list(
        individualImpl = "mcpWeights",
        sumstatImpl = "mcpRssWeights",
        multivariate = FALSE,
        coupled = list(
            individual = ".mcpCoupledWeights",
            sumstat = list(penalty = "MCP")
        ),
        grouped = ".mcpGroupedWeights"
    ),
    l0learn = list(
        individualImpl = "l0learnWeights",
        sumstatImpl = "l0learnRssWeights",
        multivariate = FALSE,
        coupled = list(
            individual = ".l0learnCoupledWeights",
            sumstat = list(penalty = "L0")
        ),
        grouped = ".l0learnGroupedWeights"
    ),
    mrmash = list(
        individualImpl = "mrmashWeights",
        sumstatImpl = "mrmashRssWeights",
        multivariate = TRUE
    ),
    dprGibbs = list(
        individualImpl = "dprGibbsWeights",
        sumstatImpl = "sdprWeights",
        multivariate = FALSE,
        # SDPR fitted across studies as SDPRX fits them.
        coupled = list(sumstat = list(sampler = ".sdprCoupledWeights"))
    ),
    enet = list(
        individualImpl = "enetWeights",
        sumstatImpl = "enetRssWeights",
        multivariate = FALSE,
        coupled = list(
            individual = ".enetCoupledWeights",
            sumstat = list(penalty = "enet")
        ),
        grouped = ".enetGroupedWeights"
    ),
    ridge = list(
        individualImpl = "ridgeWeights",
        sumstatImpl = "ridgeRssWeights",
        multivariate = FALSE,
        coupled = list(
            individual = ".ridgeCoupledWeights",
            sumstat = list(penalty = "ridge")
        ),
        grouped = ".ridgeGroupedWeights"
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
    # Spike-and-slab variational Bayes; sparsevb fits individual-level data
    # only.
    sparsevb = list(
        individualImpl = "sparsevbWeights",
        sumstatImpl = NULL,
        multivariate = FALSE
    ),
    # The Bayesian alphabet: several engines per path (.alphabetBackends).
    bayesA = list(
        individualImpl = "bayesAWeights",
        sumstatImpl = "bayesARssWeights",
        multivariate = FALSE
    ),
    bayesB = list(
        individualImpl = "bayesBWeights",
        sumstatImpl = "bayesBRssWeights",
        multivariate = FALSE
    ),
    bayesC = list(
        individualImpl = "bayesCWeights",
        sumstatImpl = "bayesCRssWeights",
        multivariate = FALSE
    ),
    bayesL = list(
        individualImpl = "bayesLWeights",
        sumstatImpl = "bayesLRssWeights",
        multivariate = FALSE
    ),
    bayesN = list(
        individualImpl = "bayesNWeights",
        sumstatImpl = "bayesNRssWeights",
        multivariate = FALSE
    ),
    bayesR = list(
        individualImpl = "bayesRWeights",
        sumstatImpl = "bayesRRssWeights",
        multivariate = FALSE
    ),
    # Sumstat-only Bayesian shrinkage (no individual-level analogue).
    prsCs = list(
        individualImpl = NULL,
        sumstatImpl = "prsCsWeights",
        multivariate = FALSE,
        # PRS-CS fitted across studies as JointPRS fits ancestries.
        coupled = list(sumstat = list(sampler = ".prsCsCoupledWeights"))
    ),
    # LDpred2, on summary statistics only, on bigsnpr or pecotmr's sampler
    # (.methodBackends).
    ldpred2Auto = list(
        individualImpl = NULL,
        sumstatImpl = "ldpred2AutoWeights",
        multivariate = FALSE
    ),
    ldpred2Grid = list(
        individualImpl = NULL,
        sumstatImpl = "ldpred2GridWeights",
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
    expanded <- .twasExpandPresets(methods)
    fmExtra <- setdiff(
        intersect(expanded, .twasFineMappingTokens()),
        .twasKnownMethodLookupNames()
    )
    regular <- setdiff(expanded, fmExtra)
    # Fine-mapping tokens carry no learner arguments of their own, so each
    # gets an empty stub entry under its `<token>_weights` key.
    methodList <- c(
        if (length(regular) > 0L) .twasMethodLookup(regular) else list(),
        .twasEmptyMethodArgs(str_c(fmExtra, "_weights"))
    )
    # Tokens come from the user input (canonical camelCase) -- the snake keys in
    # methodList are an internal detail of learnTwasWeights.
    list(tokens = unique(expanded), methodList = methodList)
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
    givenSettings <- compact(methods[[tk]])
    merged <- list_modify(base[[snake]], !!!givenSettings)
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
        "dprVb",
        "dprGibbs",
        "dprAdaptiveGibbs",
        "scad",
        "mcp",
        "l0learn",
        "sparsevb",
        "mvsusie",
        "mvsusieInf",
        "mvsusieAsh",
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
        dprVb = "dprVb",
        dprGibbs = "dprGibbs",
        dprAdaptiveGibbs = "dprAdaptiveGibbs",
        scad = "scad",
        mcp = "mcp",
        l0learn = "l0learn",
        sparsevb = "sparsevb",
        mvsusie = "mvsusie",
        mvsusie_inf = "mvsusieInf",
        mvsusie_ash = "mvsusieAsh",
        mrmash = "mrmash",
        prsCs = "prsCs",
        ldpred2Auto = "ldpred2Auto",
        ldpred2Grid = "ldpred2Grid",
        fsusie = "fsusie"
    )
    # Every method's own list key too (bayes_r_weights -> bayesR), read off
    # .twasMethodMap and matched whole: the table above spells the
    # multi-word tokens only in camelCase, so the keys .twasMethodLookup()
    # itself writes for them came back as unknown tokens ("bayes_r"), while
    # a bare snake_case token stays unknown.
    canonical <- map_chr(
        snake,
        .twasCanonicalMethod,
        snakeToCanonical = snakeToCanonical
    )
    fromMap <- set_names(names(.twasMethodMap), map_chr(.twasMethodMap, "fn"))
    hit <- is_in(names(methodList), names(fromMap))
    canonical[hit] <- fromMap[names(methodList)[hit]]
    unname(canonical)
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
    mvsusieInf = list(
        weightFn = "mvsusieInfWeights",
        rssWeightFn = "mvsusieInfRssWeights",
        fitArg = "mvsusieInfFit",
        rssFitArg = "mvsusieInfRssFit",
        methodKey = "mvsusie_inf_weights"
    ),
    mvsusieAsh = list(
        weightFn = "mvsusieAshWeights",
        rssWeightFn = "mvsusieAshRssWeights",
        fitArg = "mvsusieAshFit",
        rssFitArg = "mvsusieAshRssFit",
        methodKey = "mvsusie_ash_weights"
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

# How `token` is fitted coupled across a group's studies on `path`
# ("individual" or "sumstat"), as its capabilities name it; NULL when it is
# not.
# @noRd
.twasCoupledImpl <- function(token, path) {
    .twasMethodCapabilities[[token]][["coupled"]][[path]]
}

# The joint methods among `tokens` (mr.mash, mvsusie, fsusie): those only a
# jointSpecification groups, the same groups for both pipelines.
# @noRd
.twasJointTokens <- function(tokens) {
    keep(tokens, .twasIsMultivariateToken)
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
# context, trait, method) fitted in joint group `joint` (.noJointKey() for a
# fit of the tuple alone), or NULL when there is no hit: a row of another
# group's fit is another fit. Returns NULL silently when twasWeights is NULL
# or not a TwasWeights collection. Mirrors .fmCacheLookup
# (R/fineMappingPipeline.R).
# @noRd
.twasCacheLookup <- function(
    twasWeights,
    studyName,
    context,
    trait,
    method,
    joint = .noJointKey()
) {
    if (is.null(twasWeights)) {
        return(NULL)
    }
    if (!is(twasWeights, "TwasWeights")) {
        return(NULL)
    }
    cells <- .twasCells(twasWeights)
    hit <- .twasCellIndex(cells, studyName, context, trait, method, joint)
    if (length(hit) == 0L) {
        return(NULL)
    }
    .twrCellParts(
        twasWeights,
        cells$row[[hit[[1L]]]],
        cells$method[[hit[[1L]]]]
    )
}

# Convert a fine-mapping result (single-method susie/susie_inf row matched
# to the requested study/context/trait) into a `fittedModels` list
# suitable for `learnTwasWeights`. Pulls the trimmedFit from the matching
# entry, the fit of joint group `joint` (.noJointKey(), the tuple's own fits,
# by default). Returns a (possibly empty) list.
# @noRd
#' @importFrom checkmate assertClass
.twasFineMappingFits <- function(
    fineMappingResult,
    studyName,
    context,
    trait,
    joint = .noJointKey()
) {
    if (is.null(fineMappingResult)) {
        return(list())
    }
    assertClass(fineMappingResult, "FineMappingResultBase")
    tokens <- .twasFmFitTokens()
    found <- compact(set_names(
        map(
            tokens,
            .twasFitForToken,
            fineMappingResult = fineMappingResult,
            studyName = studyName,
            context = context,
            trait = trait,
            joint = joint
        ),
        tokens
    ))
    # compact() on an all-NULL named list leaves a zero-length names
    # attribute, which is not the bare list() the contract promises.
    if (length(found) == 0L) list() else found
}

# The fine-mapping methods whose fits TWAS reads its weights off.
# @noRd
.twasFmFitTokens <- function() {
    c("susie", "susieInf", "susieAsh", .fmMvsusieTokens(), "fsusie")
}

# The variants each of `tokens`' fits of .twasFineMappingFits() was made on,
# in the fit's order, by token (.twasFitIdsForToken()); a token whose fit
# does not say is left out.
# @noRd
.twasFineMappingFitIds <- function(
    fineMappingResult,
    studyName,
    context,
    trait,
    joint = .noJointKey(),
    tokens = .twasFmFitTokens()
) {
    if (is.null(fineMappingResult) || length(tokens) == 0L) {
        return(list())
    }
    found <- compact(set_names(
        map(
            tokens,
            .twasFitIdsForToken,
            fineMappingResult = fineMappingResult,
            studyName = studyName,
            context = context,
            trait = trait,
            joint = joint
        ),
        tokens
    ))
    if (length(found) == 0L) list() else found
}

# The variants the fit .twasFitForToken() returns was made on, in its order:
# its row's, when the row holds all of the fit, else the names the fit gives
# them; NULL when it has no fit, or neither says. A fit of unknown width
# (no PIPs or effects) says nothing of its variants.
# @noRd
.twasFitIdsForToken <- function(
    canonical,
    fineMappingResult,
    studyName,
    context,
    trait,
    joint
) {
    idx <- .twasFmTokenRows(
        fineMappingResult,
        canonical,
        studyName,
        context,
        trait,
        joint
    )
    if (length(idx) == 0L) {
        return(NULL)
    }
    parts <- .fmrRowParts(fineMappingResult, idx[[1L]])
    fit <- susieFit(parts)
    width <- .fitWidth(fit)
    if (is.na(width)) {
        return(NULL)
    }
    if (length(variants(parts)) == width) {
        return(variantIds(parts))
    }
    .fitVariableNames(fit)
}

# The fit this result holds for one token on one tuple in joint group
# `joint`, or NULL when it has none. The first matching row wins, as the
# keyed assignment it replaced did.
# @noRd
.twasFitForToken <- function(
    canonical,
    fineMappingResult,
    studyName,
    context,
    trait,
    joint
) {
    idx <- .twasFmTokenRows(
        fineMappingResult,
        canonical,
        studyName,
        context,
        trait,
        joint
    )
    if (length(idx) == 0L) {
        return(NULL)
    }
    susieFit(.fmrRowParts(fineMappingResult, idx[[1L]]))
}

# The rows of fine-mapping result `x` holding `token` for one tuple in joint
# group `joint` (in any group when NULL).
# @noRd
.twasFmTokenRows <- function(x, token, studyName, context, trait, joint) {
    idx <- .matchTupleRows(
        x,
        c(list(study = studyName, context = context, trait = trait), joint)
    )
    idx[is_in(
        str_to_lower(as.character(x$method)[idx]),
        .twasSpellingCandidates(token)
    )]
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
# per-region list; the first region carrying CV is used. The rows read are
# those of joint group `joint` (.noJointKey(), the tuple's own fits, by
# default): another group's folds can hold other samples. Returns NULL when
# no such entry for the tuple recorded CV.
# @noRd
.twasCvResultFor <- function(
    fineMappingResult,
    studyName,
    context,
    trait,
    joint = .noJointKey()
) {
    if (is.null(fineMappingResult)) {
        return(NULL)
    }
    if (!is(fineMappingResult, "FineMappingResultBase")) {
        return(NULL)
    }
    idx <- .matchTupleRows(
        fineMappingResult,
        c(list(study = studyName, context = context, trait = trait), joint)
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
#'   \item \code{mr.mash}, \code{mvsusie}, \code{mvsusieInf} and
#'         \code{fsusie} are grouped across conditions as
#'         \code{\link{fineMappingPipeline}} groups them, by
#'         \code{jointSpecification}, which they need: without one the run
#'         stops.
#'   \item \code{prsCs} runs only on summary statistics, so it is rejected
#'         on \code{QtlDataset} input.
#'   \item \code{dprVb}, \code{dprAdaptiveGibbs} and \code{sparsevb} run
#'         only on individual-level data, so they are rejected on
#'         \code{QtlSumStats} input.
#'   \item The Bayesian alphabet (\code{bayesA}, \code{bayesB},
#'         \code{bayesC}, \code{bayesL}, \code{bayesN}, \code{bayesR})
#'         runs on both, each method on several engines chosen with
#'         \code{backend}: see \code{\link{bayesAlphabetWeights}} and
#'         \code{\link{bayesAlphabetRssWeights}}. On summary statistics
#'         a row with a genotype reference defaults to pecotmr's own
#'         sampler, which reads its factor of the LD.
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
#'   \code{<method>_weights = args} entries. On individual-level input the
#'   default is \code{"default"}, which is \code{"fastDefault"} (\code{susie},
#'   \code{susieInf}, \code{mrash}, \code{enet}, \code{lasso}, \code{mcp},
#'   \code{scad} and \code{l0learn}) with \code{bayesR} and \code{bayesC}.
#'   \code{susie} and \code{susieInf} read their fits from a
#'   \code{fineMappingResult}, so either preset skips them, with a warning,
#'   where no \code{fineMappingResult} holds them; naming them outside a preset
#'   without one is an error. For QtlSumStats input
#'   the default, which either preset also names, is
#'   \code{c("lasso", "prsCs", "dprGibbs")}, which there run lassosum, PRS-CS
#'   and SDPR; the fine-mapping methods are not in it, and are requested by
#'   name together with a \code{fineMappingResult}. A
#'   multivariate method (\code{mrmash}, \code{mvsusie}, \code{mvsusieInf},
#'   \code{fsusie}) fits the groups \code{jointSpecification} makes and any
#'   other method each (context, trait) alone, so the two kinds can be
#'   requested at once.
#' @param backend The engine a method with several runs on, as a named
#'   character vector keyed by method, e.g.
#'   \code{c(bayesR = "hibayes", bayesL = "bglr")}; \code{NULL} (the
#'   default) leaves every method on its default. The Bayesian alphabet
#'   methods have several (see \code{\link{bayesAlphabetWeights}} and
#'   \code{\link{bayesAlphabetRssWeights}}; \code{bayesB} has one on summary
#'   statistics), and on summary statistics so does \code{lasso}
#'   (\code{"pecotmr"} or \code{"bigsnpr"}, which is slow: see
#'   \code{\link{lassosumRssWeights}}). On individual-level data
#'   \code{lasso} and \code{enet} run on \code{"glmnet"} (the default) or
#'   \code{"grpreg"}, \code{mcp} and \code{scad} on \code{"ncvreg"} (the
#'   default) or \code{"grpreg"}, and \code{l0learn} on \code{"l0learn"}
#'   (the default) or \code{"abess"}; each also fits groups of conditions
#'   with a group penalty on glmnet, grpreg and abess, and \code{bayesC} and
#'   \code{bayesN} on \code{"bglr"} fit each group a
#'   \code{jointSpecification} makes jointly (see there). Naming a
#'   method that is not being run, one with a single implementation, or an
#'   engine the method lacks on the input is an error, and so is a
#'   \code{backend} inside a \code{methods} entry.
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
#'   axes the joint methods -- \code{mrmash}, \code{mvsusie},
#'   \code{mvsusieInf} and \code{fsusie} -- group conditions across. It
#'   means what it means to
#'   \code{\link{fineMappingPipeline}}, so a run groups conditions exactly as
#'   a fine-mapping run given the same specification and selection did. When
#'   NULL (the default) every method is fitted one (study, context, trait) at
#'   a time, and the joint methods stop the run, since they need this to say
#'   which conditions they fit together. Otherwise give a character vector
#'   of axes (one spec) or a list of specs. A spec is a character vector drawn
#'   from \code{"study"}, \code{"context"} and \code{"trait"}, or
#'   \code{list(axes = <axes>, scope = <named list>)}, where \code{scope}
#'   (keys \code{study}, \code{context}, \code{trait}) limits the fit to the
#'   listed values, each of which must exist in \code{data}. A list-form
#'   spec always goes inside the outer list, even on its own. Other methods
#'   run per (study, context, trait) as usual. \code{mrmash} is fitted on
#'   each group; \code{mvsusie}, \code{mvsusieInf} and \code{fsusie} are
#'   read off \code{fineMappingResult}'s fit over exactly the group's
#'   conditions, each
#'   condition's own column of it, and a group the result holds no such fit
#'   of stops the run: fine-map and build weights with the same
#'   \code{jointSpecification} and \code{contexts} / \code{traitId}
#'   selection. The \code{"study"} axis needs more than one study, so it is
#'   rejected on a \code{QtlDataset}. On a \code{MultiStudyQtlDataset} it
#'   stacks the individual-level studies, each study's conditions observed
#'   on its own samples only: the study axis alone makes one group per
#'   (context, trait) two or more of them carry, and a composed spec
#'   including \code{"study"} one group per value of the axes it leaves out
#'   (each trait for \code{c("study", "context")}). On those groups it reads
#'   the mvSuSiE methods off the fine-mapping fits (no genotype is read) and
#'   fits
#'   \code{mrmash}, each condition on its own samples (mr.mash's
#'   \code{missing_y_method = "approximate"}); the penalized methods fit
#'   them as described below. The summary-statistics studies are fitted as
#'   well. On summary statistics a group along a coupled axis also fits the
#'   penalized methods under the penalty on the difference between their
#'   effects (see \code{coupling}), each condition against its study's LD
#'   reference, and \code{mrmash} and the mvSuSiE methods
#'   take each study's own reference too (see \code{\link{mrmashRssWeights}}
#'   and \code{\link{mvsusieRssWeights}}).
#'   On individual-level data, \code{bayesC} and \code{bayesN} on
#'   \code{backend = "bglr"} also fit each group of two or more conditions
#'   at once, through BGLR's \code{Multitrait()} (see
#'   \code{\link{bayesAlphabetWeights}}), alongside their fits of one
#'   (study, context, trait) at a time: a study's contexts or traits, with
#'   an unstructured residual covariance unless
#'   \code{BglrOptions(resCov = )} sets another, and on a
#'   \code{MultiStudyQtlDataset} the individual-level studies on the study
#'   axis alone, stacked as for the coupled methods, where the residual
#'   covariance is diagonal: their samples are disjoint. On a composed spec
#'   including \code{"study"} the residual covariance is BGLR's recursive
#'   (\code{REC}) structure, free among a study's contexts or traits and zero
#'   across studies, with that structure's own prior. Without a
#'   specification they run one tuple at a time only.
#'   On individual-level data the penalized methods (\code{lasso},
#'   \code{enet}, \code{ridge}, \code{mcp}, \code{scad}, \code{l0learn})
#'   also fit each group of two or more conditions a specification makes at
#'   once, alongside their fits of one tuple at a time, on the engine
#'   \code{backend} puts them on; without a specification they run one tuple
#'   at a time only. A group along an axis \code{coupling} names is fitted under
#'   the penalty on the difference between the conditions' effects, with
#'   the engine's group penalty as well where it has one; any other group
#'   under the group penalty alone: each variant's effects across the
#'   conditions are one group, so a variant is in every condition's weights
#'   or in none (see \code{\link{lassoWeights}}). The group penalty is
#'   glmnet's \code{"mgaussian"} family (\code{lasso}, \code{enet}) for a
#'   study's contexts or traits on its samples, grpreg's group penalties
#'   (\code{backend = "grpreg"}: \code{lasso}, \code{enet}, \code{mcp},
#'   \code{scad}) and abess's best subset of groups (\code{"abess"}:
#'   \code{l0learn}) for any group, each condition on its own samples; ridge
#'   has none, and fits each condition at one lambda chosen across them.
#'   ncvreg and L0Learn, the default engines of \code{mcp}, \code{scad} and
#'   \code{l0learn}, have no group penalty, so with no \code{backend} named
#'   those methods take their group fits on grpreg and abess, and their
#'   coupled ones on ncvreg and L0Learn; named \code{"ncvreg"} or
#'   \code{"l0learn"}, they fit a group only when it is coupled. A method
#'   named in \code{methods} that cannot fit a specification's groups is an
#'   error. On summary statistics every penalized method fits each group
#'   under the group penalty on pecotmr's own solver (not
#'   \code{lasso} on \code{backend = "bigsnpr"}), every condition at one
#'   grid point; a coupled group takes it with the coupling. A trait group
#'   there is fitted over the union of its genes' variants where their
#'   statistics differ, each gene's weights reaching as far as its
#'   statistics do. Within a study the groups keep every sample some
#'   condition observes unless \code{completeSamples = TRUE}.
#' @param fineMappingResult Optional \code{QtlFineMappingResult}. When supplied,
#'   its SuSiE / SuSiE-inf / SuSiE-ash trimmed fits for the matching (study,
#'   context, trait) tuples are injected into \code{learnTwasWeights} via
#'   \code{fittedModels} so SuSiE-family weight methods reuse the prior fit
#'   instead of refitting. Its fits and cross-validation are reused whole,
#'   so a row narrowed since its fit warns (condition class
#'   \code{pecotmrFitCoverage}). Each lookup reads the fits of the group it
#'   means: a method fitted per tuple reads the tuple's own fit, a joint
#'   method its group's (see \code{jointSpecification}), never a fit of the
#'   tuple in another joint group. On a \code{MultiStudyQtlDataset} whose
#'   \code{jointSpecification} has the study axis, the rows of a fit across
#'   its studies (\code{fineMappingPipeline(jointSpecification = "study")})
#'   are read one by one: each study's weights are its own column of that
#'   fit, over all of the fit's variants (zero at a variant the study lacks),
#'   with its row's cross-validation, in rows that keep the fit's
#'   \code{jointStudies}.
#' @param twasWeights Optional \code{\link{TwasWeights}} resume cache. For each
#'   requested \code{(study, context, trait, method)} entry already present in
#'   this collection, the cached entry is copied through and the
#'   corresponding weight fit is skipped. Only the un-cached methods are fit;
#'   the cached and fresh entries are combined in the returned collection.
#'   Per-entry matching mirrors the \code{fineMappingResult} cache in
#'   \code{\link{fineMappingPipeline}}. A joint method (\code{mvsusie},
#'   \code{mr.mash}) fits a group of conditions together, so its group is
#'   resumed only when the cache holds every condition's entry from a fit of
#'   that same group (its \code{jointStudies} / \code{jointContexts} /
#'   \code{jointTraits}); otherwise the whole group is refit. A method fitted
#'   per tuple resumes only from the tuple's own fit.
#' @param mashPrior Optional \code{\link{MashPrior}} supplying data-driven
#'   priors for mr.mash, one per joint specification
#'   (\code{\link{mashPipeline}}); \code{NULL} to skip. Each mr.mash group
#'   takes the prior over its own joint specification's axes, on its own
#'   conditions, put on the per-SD scale mr.mash fits on (see
#'   \code{\link{MashPrior}}). A \code{QtlSumStats} is not cross-validated,
#'   so its groups take the full-data priors, and a \code{MashPrior} that
#'   holds per-fold priors alone is an error there.
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
#'   Every method keeps a fit, and the fit names the variants it was fitted
#'   on. \code{"full"} is honoured by the mr.mash engines; on summary
#'   statistics also by \code{lasso}, \code{enet}, \code{ridge}, \code{scad},
#'   \code{mcp} and \code{l0learn}, whose full fit keeps every candidate's
#'   coefficients; and on individual-level data by \code{lasso}, \code{enet},
#'   \code{ridge}, \code{scad}, \code{mcp}, \code{l0learn}, the Bayesian
#'   alphabet, \code{mrash} and the DPR methods, whose full fit is the
#'   engine's whole object where the slim one keeps the selected fit or the
#'   posterior summaries, without any copy of \code{X} unless \code{retainX}
#'   asks for it. Every other method keeps its fit whole or not at all, except
#'   the fine-mapping methods (\code{susie}, \code{susieInf},
#'   \code{susieAsh}, \code{mvsusie}, \code{mvsusieInf}, \code{fsusie}),
#'   whose rows keep none:
#'   their weights come from the \code{fineMappingResult}'s fit, which is where
#'   \code{ctwasPipeline()} reads it. Requesting \code{"full"} for a run where
#'   no method can honour it is an error, and a mixed run warns rather than
#'   dropping it in silence.
#' @param retainX Logical. Keep a fit's copy of \code{X}, which can be its
#'   largest part, at \code{fitRetention = "full"} (default \code{FALSE}). An
#'   error at any other \code{fitRetention}, and on a run where no method's fit
#'   has a copy: only \code{lasso}, \code{enet}, \code{ridge}, \code{scad},
#'   \code{mcp}, \code{mrash} and the DPR methods do. No fine-mapping fit holds
#'   one, so \code{\link{fineMappingPipeline}} has no such argument, and summary
#'   statistics have no \code{X}: on a \code{MultiStudyQtlDataset} it reaches
#'   the individual-level studies.
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
#' @param coupling The axes along which the penalized methods couple a
#'   group's conditions, each with its grid of strengths \eqn{c}: a named
#'   list with any of \code{study}, \code{context} and \code{trait}, e.g.
#'   \code{list(study = c(0.01, 0.1, 1, 10), context = c(0.1, 1))}.
#'   Coupled conditions are fitted under the penalty
#'   \eqn{(c/2) \sum (b_{jk} - b_{jk'})^2} on the difference between their
#'   per-SD effects, as PROSPER couples ancestries: a pair of a group's
#'   conditions differing along one coupled axis alone takes that axis's
#'   \eqn{c}, any other pair none, and several coupled axes make a grid of
#'   every combination of their strengths. Each condition's loss is its own
#'   mean squared error, and each condition takes the grid point its own
#'   held-out data score best: on individual-level data, \eqn{(\lambda, c)}
#'   on folds that hold a sample out of every condition at once; on summary
#'   statistics, \eqn{(s, \lambda, c)} on held-out statistics drawn from
#'   each condition's own LD, every condition fitted against its study's LD
#'   reference (the draws are made for each condition apart, as if no two
#'   shared samples). \eqn{c = 0} gives the conditions' separate fits and a
#'   large \eqn{c} one pooled fit. The default couples studies only, at
#'   \code{c(0.01, 0.1, 1, 10)}, which spans nearly separate to nearly
#'   pooled; \code{NULL} couples nothing. Only the groups a
#'   \code{jointSpecification} makes are coupled, so naming \code{coupling}
#'   in a run without one is an error. Coupled on individual-level data:
#'   \code{lasso}, \code{enet} and \code{ridge} on glmnet, with the group
#'   penalty of glmnet's \code{"mgaussian"} family as well on the conditions
#'   of a group that share samples and are not coupled to one another;
#'   \code{mcp} and \code{scad} on ncvreg; \code{l0learn} on L0Learn; and
#'   on grpreg and abess every method with its group penalty as well.
#'   ncvreg and grpreg rescale their columns, the coupling among them, so
#'   their penalty is approximate on a coupled group (MCP's exact when every
#'   variant has the same coupling), which a warning says once a session;
#'   their designs are dense. On summary statistics those methods are
#'   fitted by pecotmr's own solver whatever their \code{backend}, with the
#'   group penalty across the conditions as well (the L0 family swapping
#'   whole groups, as L0Learn swaps single variants).
#'   \code{l0learn} coupled takes \code{penalty = "L0"} only. The methods'
#'   fits of one condition at a time run too; the coupled ones are told
#'   apart by their joint columns (\code{\link{WeightSelectionParam}}'s
#'   \code{jointAxes} selects them). On summary statistics \code{prsCs}
#'   and \code{dprGibbs} are fitted across a coupled group too, as JointPRS
#'   and SDPRX fit ancestries; \eqn{c} does not apply to them, though
#'   \code{NULL} turns them off as well. Under \code{prsCs} the conditions
#'   share PRS-CS's continuous shrinkage prior, as PRS-CSx shares it, and
#'   their effects on a variant are correlated, the correlation learned on
#'   \eqn{[0, 0.99]} (or fixed by \code{PrsCsOptions(rho = )}); its weights
#'   stay standardized, so \code{maf} is refused there. A \code{mashPrior}
#'   gives both samplers their starting correlation, its covariances averaged
#'   over their weights, and \code{dprGibbs} its kinds' starting weights,
#'   each covariance's weight to the kind holding the conditions it covers.
#'   Under \code{dprGibbs} each non-null component of SDPR's
#'   Dirichlet process mixture carries a variant's effects in a set of the
#'   conditions, correlated: \code{\link{SdprOptions}}' \code{sharing} sets
#'   which sets, and \code{rho} the correlation (learned by default). Its
#'   \code{perVariantSampleSize}, \code{array} and \code{optLlk = 2}
#'   describe one condition, so they are refused there.
#' @param completeSamples Logical (length 1). Individual-level data only:
#'   whether a joint group within a study (a study's contexts or traits, or
#'   each study's share of a group across studies) is fitted only on the
#'   samples observed in every one of its conditions. \code{FALSE}, the
#'   default, keeps every sample observed in some condition, each condition
#'   fitted on its own samples; a variant is then kept only when it passes
#'   \code{genotypeFilterParam}'s MAF, MAC and X-variance cutoffs on every
#'   condition's own samples, as well as on all of them. An engine that
#'   takes no missing phenotype (glmnet's \code{"mgaussian"} group fit,
#'   fSuSiE, and mr.mash where the installed mr.mashr has no
#'   \code{missing_y_method}) still fits on the samples observing every
#'   condition itself.
#' @param residualCorrelation Summary statistics only: the correlation of
#'   the conditions' z-scores under the null, which carries their residual
#'   correlation and their samples' overlap (what
#'   \code{\link{mashResidualCorrelation}} estimates), a matrix whose rows
#'   and columns are named by condition, \code{"study:context:trait"}.
#'   Conditions on one LD reference over the same variants are read as
#'   sharing samples by it: the penalized methods draw their held-out splits
#'   correlated by it, and \code{prsCs} and \code{dprGibbs} fitted across
#'   them correlate their noise in their likelihood: \code{dprGibbs} by it,
#'   and \code{prsCs} by its conditions' residual covariance, learned with
#'   it standing for their phenotypes' cross-products (exact where they share
#'   every sample; the residual correlation is reported as
#'   \code{residualEst}). Without it, the residual
#'   correlation of a \code{mashPrior}'s prior over the group's axes is used;
#'   given with a \code{mashPrior}, it overrides that one, with a warning.
#'   \code{NULL} (the default) with no \code{mashPrior} reads every
#'   condition's noise as independent, which favours stronger sharing where
#'   the conditions share samples.
#' @param ensembleParam SR-TWAS ensemble settings, built with
#'   \code{\link{EnsembleParam}}: \code{enabled} (default \code{FALSE}),
#'   \code{r2Threshold}, \code{solver}, \code{alpha}, \code{acrossStudies}
#'   and \code{transferScale}. On individual-level data stacking reads
#'   out-of-fold predictions, so \code{enabled = TRUE} requires
#'   \code{crossValidationParam} with \code{folds >= 2} and is an error
#'   otherwise --- it previously returned a result quietly missing its
#'   ensemble row. On summary statistics it reads held-out splits of the
#'   statistics instead. With \code{acrossStudies = TRUE} and a
#'   \code{"study"} \code{jointSpecification}, each study's ensemble also
#'   stacks the other studies' fits, as MAGEPRO and LEOPARD combine
#'   ancestries; those rows carry \code{jointStudies}.
#' @param estimatePi If TRUE, estimate spike-and-slab sparsity from mr.ash
#'   before BGLR / qgg spike-and-slab methods that consume it.
#' @param h2Estimate,h2EstimateParam The heritability given, trait by trait,
#'   to the methods that can start from one: the Bayesian alphabet
#'   (\code{bayesA} to \code{bayesR}) and, on summary statistics, LDpred2
#'   (\code{ldpred2Auto}'s start, and the centre of \code{ldpred2Grid}'s
#'   grid) and the joint mr.mash fits, whose residual variances it sets to
#'   each condition's \eqn{1 - h^2}, as SuSiE-RSS's, with the correlations of
#'   their residual covariance (see \code{\link{mrmashRssWeights}}). Neither
#'   (the default) gives none, so each keeps its own default; both is an
#'   error. \code{h2Estimate} is one number for every trait, or an
#'   \code{\linkS4class{H2Estimate}}, as \code{\link{estimateLocalH2}}
#'   returns, which gives each selected (study, context, trait) its row's
#'   estimate and must hold exactly one row for each.
#'   \code{h2EstimateParam}, an \code{\link{H2EstimateParam}}, estimates
#'   each trait's first, once, on the variants, samples and LD its methods
#'   fit (on a \code{QtlDataset} in-sample, its phenotype statistic, REML,
#'   included); on summary statistics the reference is taken to share no
#'   samples with the study, and LDpred2 estimates under the settings itself.
#'   A value is used floored at 0.001 and capped at 0.99, and a trait with no
#'   estimate keeps each method's default, with a warning. A method's own
#'   setting in \code{methods} wins.
#' @param residualizationParam Covariate residualization settings, built with
#'   \code{\link{ResidualizationParam}} and forwarded to
#'   \code{\link{residualizedPhenotypes}} /
#'   \code{\link{residualizedGenotypes}}: \code{phenotypeCovariates} and
#'   \code{genotypeCovariates} name which covariates to regress out
#'   (\code{NULL}, the default, uses every available one), and
#'   \code{residualizePhenotype} / \code{residualizeGenotype} turn each side
#'   off. \code{covariateScope} sets whose covariates a joint fit takes
#'   where its outcomes span contexts: by default each outcome's phenotypes
#'   its own context's, and the genotypes they share the genotype covariates
#'   and the phenotype covariates every context shares.
#'
#'   \code{QtlDataset} / \code{MultiStudyQtlDataset} input only. A
#'   \code{QtlSumStats} run has no covariates to regress out, so the bundle
#'   is \strong{ignored} there rather than refused --- unlike
#'   \code{crossValidationParam}, which errors. Residualization is on by
#'   default,
#'   so refusing a non-default value would reject the default bundle; CV is
#'   off by default, so a non-default value there is an explicit request for
#'   something the input cannot do.
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
#' @return A \code{\link{TwasWeights}} collection with one row per gene --
#'   \code{(study, context, trait)} and joint group -- holding its methods.
#'   Each row's LD reference is none for individual-level fits and that of
#'   the input row(s) for RSS-derived fits.
#'   Its \code{seqinfo()} names the genome build the inputs name; inputs
#'   naming different builds are refused.
#' @examples
#' data(qtlDatasetExample)
#' if (requireNamespace("PlinkMatrix", quietly = TRUE)) {
#'     twasWeightsPipeline(
#'         qtlDatasetExample,
#'         methods = "lasso",
#'         cisWindow = 1e6
#'     )
#' }
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
    .qtlRunWork(
        .twasJointSpecGridWork(spec, tokens, contexts, traits, marker, ctx),
        ctx$data,
        ctx$genotypeFilter,
        ctx$residualizationParam
    )
}

# .twasRunJointSpecGrid()'s fits, planned, for a run with others to fit
# alongside.
# @noRd
.twasJointSpecGridWork <- function(
    spec,
    tokens,
    contexts,
    traits,
    marker,
    ctx
) {
    .qtlWorkOf(
        map(
            seq_along(ctx$xRegions),
            .twasJointSpecRegionWork,
            spec = spec,
            tokens = tokens,
            contexts = contexts,
            traits = traits,
            marker = marker,
            ctx = ctx
        ),
        .twasJointGridResult,
        list(allLabs = map_chr(ctx$xRegions, .twasRegionLabel))
    )
}

# The regions' joint results merged.
# @noRd
.twasJointGridResult <- function(allRegions, allLabs) {
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

# A jointSpecification's groups must be ones every penalized method the
# caller names in `methods` fits, coupled or with its engine's group penalty;
# a preset's methods that fit neither are left to their fits one tuple at a
# time, as without a specification.
# @noRd
.twasAssertSpecFits <- function(parsedJointSpec, methods, norm, marker) {
    if (length(parsedJointSpec) == 0L || .twasIsPreset(methods)) {
        return(invisible(NULL))
    }
    walk(
        parsedJointSpec,
        .twasAssertOneSpecFits,
        tokens = intersect(norm$tokens, names(.groupedLayouts)),
        cfg = .jpConfig(marker)
    )
}

# @noRd
.twasAssertOneSpecFits <- function(spec, tokens, cfg) {
    misfit <- discard(
        tokens,
        .twasJointAdmits,
        layout = .jointSpecLayout(spec$axes),
        axes = spec$axes,
        path = "individual",
        cfg = cfg
    )
    if (length(misfit) == 0L) {
        return(invisible(NULL))
    }
    abort(c(
        glue(
            "jointSpecification joins {str_flatten(spec$axes, ' x ')}, ",
            "which {str_flatten(misfit, ', ')} cannot fit on the engine the ",
            "run puts it on without coupling one of those axes."
        ),
        i = str_c(
            "Choose a backend with a group penalty there (`backend`), ",
            "couple one of the axes (`coupling`), or leave the method out."
        )
    ))
}

# `norm` narrowed to the methods among `keep`.
# @noRd
.twasNormKeep <- function(norm, keep) {
    held <- is_in(.twasTokensFromMethodList(norm$methodList), keep)
    list_assign(norm, tokens = keep, methodList = norm$methodList[held])
}

# The QtlDataset worker's argument checks, grouped so its entry point names
# the validation phase once instead of listing four calls inline.
# @noRd
.twasQdsAssertParams <- function(
    methods,
    fitRetention,
    retainX,
    crossValidationParam,
    ensembleParam
) {
    .twasWarnUnretainableDetail(fitRetention, methods, "QtlDataset")
    .twasCheckRetainX(retainX, fitRetention, methods)
    .assertMethodParam(
        crossValidationParam,
        "CrossValidationParam"
    )
    .assertMethodParam(ensembleParam, "EnsembleParam")
    .ensembleAssertCv(ensembleParam, crossValidationParam)
    .ensembleAssertAcross(ensembleParam, list(), 1L, "QtlDataset")
    invisible(NULL)
}

# The QtlDataset method of twasWeightsPipeline(). Each stage returns only the
# values it derives; nothing is grafted onto a captured environment. `cfg`
# (.twasJointCfg()) is the one record of the run's settings, and the pipeline
# marker is built from it BEFORE the joint phase, so explicit joint fits run
# under the same marker and engine args as the per-tuple loop. The joint phase
# consumes the joint methods (.twasJointTokens()); `norm` comes back holding
# only the methods that still have to go through the per-tuple loop.
# @noRd
.twasPipelineQtlDataset <- function(
    data,
    methods = "default",
    backend = NULL,
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
    h2Estimate = NULL,
    h2EstimateParam = NULL,
    usePCA = FALSE,
    nPCs = 10L,
    fitRetention = c("slim", "none", "full"),
    retainX = FALSE,
    residualizationParam = ResidualizationParam(),
    phenotypeFilterParam = PhenotypeFilterParam(),
    verbose = 1,
    coupling = list(study = c(0.01, 0.1, 1, 10)),
    completeSamples = FALSE,
    seed = NULL
) {
    phenotypeFilter <- .twasResolvePhenotypeFilter(phenotypeFilterParam)
    fitRetention <- arg_match(fitRetention)
    .twasAssertCoupling(coupling)
    .twasAssertCouplingGroups(coupling, !missing(coupling), jointSpecification)
    assertFlag(completeSamples)
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
        retainX,
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
        h2 = .twasQdsH2Spec(
            h2Estimate,
            h2EstimateParam,
            data,
            contexts,
            traitId,
            region
        ),
        fitRetention = fitRetention,
        retainX = retainX,
        residualizationParam = residualizationParam,
        verbose = verbose,
        seed = seed,
        twasWeights = twasWeights,
        backend = backend,
        coupling = coupling,
        completeSamples = completeSamples
    )
    resolvedTokens <- .twasQdsResolveTokens(
        jointSpecification,
        methods,
        backend,
        fineMappingResult,
        .cvResolve(crossValidationParam),
        data
    )
    setup <- .twasJointSetup(data, resolvedTokens$norm, cfg)
    .twasAssertSpecFits(
        resolvedTokens$parsedJointSpec,
        methods,
        resolvedTokens$norm,
        setup$marker
    )
    joint <- .twasQdsJointPhase(
        resolvedTokens$parsedJointSpec,
        setup$marker,
        setup$ctx,
        contexts,
        traitId
    )
    if (joint$done) {
        return(.withGenome(
            .twasRecordDataTypes(
                .twasQdsJointOnly(
                    joint$work,
                    data,
                    resolved$genotypeFilter,
                    residualizationParam
                ),
                data
            ),
            build
        ))
    }
    grid <- .twasQdsResolveGrid(data, contexts, traitId, region)
    gridWork <- .twasQdsUnivariateWork(
        grid,
        marker = setup$marker,
        data = data,
        xRegions = resolved$xRegions,
        norm = joint$norm,
        fineMappingResult = fineMappingResult,
        twasWeights = twasWeights,
        mashPriors = setup$ctx$mashPriors,
        cisWindow = cisWindow,
        phenotypeFilter = phenotypeFilter,
        residualizationParam = residualizationParam,
        genotypeFilter = resolved$genotypeFilter,
        usePCA = usePCA,
        nPCs = nPCs,
        verbose = verbose
    )
    # Every path planned, so one read of a window serves them all.
    .withGenome(
        .twasRecordDataTypes(
            .qtlRunWork(
                .qtlWorkOf(list(gridWork, joint$work), .twasQdsAssembleResults),
                data,
                resolved$genotypeFilter,
                residualizationParam
            ),
            data
        ),
        build
    )
}

# The result of a run's paths: the grid's fits, then the joint specs'.
# @noRd
.twasQdsAssembleResults <- function(results) {
    .twasQdsAssemble(results[[1L]], results[[2L]])
}

#' @rdname twasWeightsPipeline
#' @export
setMethod(
    "twasWeightsPipeline",
    "QtlDataset",
    .twasPipelineQtlDataset
)

# Unpack the MashPrior bundle: its priors go to the joint groups, which take
# each mr.mash group's own (.jointTwasGroupPriors()), and its fold partition
# to the CV machinery.
# @noRd
.twasQdsUnpackMash <- function(
    mashPrior,
    samplePartition,
    norm
) {
    mp <- .unpackMashPrior(mashPrior, samplePartition)
    .twasWarnUnusedMashPrior(mashPrior, norm$tokens)
    list(samplePartition = mp$samplePartition, mashPriors = mp$priors)
}

# The method arguments with mr.mash's data-driven prior set, when it has
# arguments and a prior; every other method's are left as they are.
# @noRd
.twasWithMrmashPrior <- function(methodList, prior) {
    if (is.null(prior) || !is_in("mrmash_weights", names(methodList))) {
        return(methodList)
    }
    list_assign(
        methodList,
        mrmash_weights = list_assign(
            methodList$mrmash_weights,
            dataDrivenPriorMatrices = prior
        )
    )
}

# A mashPrior no requested method takes is ignored, and said to be.
# @noRd
.twasWarnUnusedMashPrior <- function(mashPrior, tokens) {
    if (!is.null(mashPrior) && !is_in("mrmash", tokens)) {
        msg <- glue(
            "`mashPrior` was supplied but 'mrmash' is not among `methods`; ",
            "the data-driven prior is ignored."
        )
        warn(msg)
    }
}

# The priors a summary-statistics run's mr.mash groups take theirs from
# (.jointMrmashRssPrior()): a MashPrior's full-data priors, or NULL without
# one or without mr.mash. Summary statistics are not cross-validated, so a
# MashPrior holding per-fold priors alone gives them none.
# @noRd
.twasQssMashPriors <- function(mashPrior, tokens) {
    if (is.null(mashPrior) || !is_in("mrmash", tokens)) {
        return(NULL)
    }
    priors <- .unpackMashPrior(mashPrior)$priors$full
    if (is.null(priors)) {
        abort(glue(
            "twasWeightsPipeline(QtlSumStats): `mashPrior` holds per-fold ",
            "priors only. Summary statistics are not cross-validated, so ",
            "mr.mash needs its full-data priors (`fullFit`)."
        ))
    }
    priors
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

# Explicit jointSpecification path: plan the per-spec axis dispatcher for the
# joint methods (.twasJointTokens()), which need one. Returns
# list(done, work, norm): `done` means its fits are the whole answer;
# otherwise `norm` is the normalization left for the per-tuple loop below,
# the joint methods gone. Planned, not run: they share their reads with the
# per-tuple fits.
# @noRd
.twasQdsJointPhase <- function(
    parsedJointSpec,
    marker,
    ctx,
    contexts,
    traitId
) {
    norm <- ctx$norm
    .assertJointSpecForJointMethods(
        parsedJointSpec,
        .twasJointTokens(norm$tokens),
        "twasWeightsPipeline(QtlDataset)"
    )
    if (length(parsedJointSpec) == 0L) {
        return(list(done = FALSE, work = .qtlWorkValue(NULL), norm = norm))
    }
    work <- .twasJointSpecsWork(
        parsedJointSpec,
        norm$tokens,
        marker,
        ctx,
        contexts,
        traitId
    )
    keep <- setdiff(norm$tokens, .twasJointTokens(norm$tokens))
    if (length(keep) == 0L) {
        return(list(done = TRUE, work = work))
    }
    list(done = FALSE, work = work, norm = .twasNormKeep(norm, keep))
}

# Parse the joint spec and normalize the methods, then run the capability,
# per-input argument and fine-mapping gates. The mash-prior unpacking and the
# joint phase follow once the pipeline marker exists (.twasJointSetup()), so
# that explicit joint fits run under the same marker as everything else.
# @noRd
.twasQdsResolveTokens <- function(
    jointSpecification,
    methods,
    backend,
    fineMappingResult,
    cvCfg,
    data
) {
    parsedJointSpec <- parseJointSpecification(jointSpecification, data)
    norm <- .twasNormalizeMethods(.twasRunMethods(methods, fineMappingResult))
    .twasCheckMethodCapabilities(norm$tokens, "QtlDataset")
    .twasCheckMethodArgsForInput(
        .twasMethodListArgs(norm$methodList),
        "QtlDataset"
    )
    .twasCheckBackend(backend, norm$tokens, "individual")
    .twasCheckFineMappingMethods(
        norm$tokens,
        fineMappingResult,
        "QtlDataset",
        cvFolds = cvCfg$folds,
        fitFullData = cvCfg$fitFullData
    )
    norm$methodList <- .twasInjectBackends(norm$methodList, backend)
    list(parsedJointSpec = parsedJointSpec, norm = norm)
}

# Check the pipeline's `backend`: a named character vector keyed by the
# methods being run, each naming an engine that method has on `path`.
# @noRd
.twasCheckBackend <- function(backend, tokens, path) {
    if (is.null(backend)) {
        return(invisible(NULL))
    }
    if (!is.character(backend) || !is_named(backend)) {
        abort(str_c(
            "`backend` must be a named character vector keyed by method, ",
            "e.g. c(bayesR = \"hibayes\")."
        ))
    }
    stray <- setdiff(names(backend), tokens)
    if (length(stray) > 0L) {
        abort(glue(
            "`backend` names methods that are not being run: ",
            "{str_flatten(stray, ', ')}."
        ))
    }
    walk(names(backend), .twasCheckOneBackend, backend = backend, path = path)
    invisible(NULL)
}

# One `backend` entry: the method must have several engines on `path`, and
# the entry must name one of them.
# @noRd
.twasCheckOneBackend <- function(token, backend, path) {
    choices <- names(.methodBackends[[path]][[token]])
    if (!.methodHasBackends(token, path)) {
        abort(glue(
            "`backend`: {token} has one implementation on this input, so ",
            "there is no backend to choose."
        ))
    }
    if (!is_in(backend[[token]], choices)) {
        abort(glue(
            "`backend` for {token} must be one of ",
            "{str_flatten(choices, ', ')} on this input."
        ))
    }
    invisible(NULL)
}

# The pipeline's `backend` choices, added to the arguments of the methods
# they name. Done after the caller's own arguments are checked, which refuse
# a `backend` there.
# @noRd
.twasInjectBackends <- function(methodList, backend) {
    if (is.null(backend)) {
        return(methodList)
    }
    reduce(
        names(backend),
        .twasInjectOneBackend,
        backend = backend,
        .init = methodList
    )
}

# @noRd
.twasInjectOneBackend <- function(methodList, token, backend) {
    key <- .twasMethodMap[[token]]$fn
    methodList[[key]]$backend <- backend[[token]]
    methodList
}

# A MultiStudyQtlDataset's methods normalized for its individual-level fits,
# with the run's `backend` choices of an individual-level engine for a method
# among them in that method's arguments, as each study's own run puts them
# (.twasQdsResolveTokens()). A choice for the summary-statistics path alone
# is that path's to read.
# @noRd
.twasNormWithBackends <- function(methods, backend) {
    norm <- .twasNormalizeMethods(methods)
    chosen <- backend[is_in(names(backend), norm$tokens)]
    individual <- keep(
        names(chosen) %||% character(0),
        .twasIndividualBackend,
        backend = chosen
    )
    list_assign(
        norm,
        methodList = .twasInjectBackends(norm$methodList, chosen[individual])
    )
}

# Whether a token's chosen engine is one it has on individual-level data.
# @noRd
.twasIndividualBackend <- function(token, backend) {
    is_in(backend[[token]], names(.methodBackends$individual[[token]]))
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
    retainX,
    residualizationParam,
    verbose,
    seed,
    h2 = NULL,
    twasWeights = NULL,
    backend = NULL,
    coupling = NULL,
    completeSamples = FALSE
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
        h2 = h2,
        fitRetention = fitRetention,
        retainX = retainX,
        residualizationParam = residualizationParam,
        verbose = verbose,
        seed = seed,
        twasWeights = twasWeights,
        backend = backend,
        coupling = coupling,
        completeSamples = completeSamples
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
        # The methods a group across studies stacks, which
        # .twasStudyAxisFits() alone sets.
        ensembleParam = list_assign(
            .ensembleResolve(cfg$ensembleParam),
            stackTokens = cfg$stackTokens
        ),
        standardized = .qtlScaleResiduals(cfg$residualizationParam),
        fitRetention = cfg$fitRetention,
        retainX = cfg$retainX,
        estimatePi = cfg$estimatePi,
        h2 = cfg$h2,
        verbose = cfg$verbose,
        seed = cfg$seed,
        coupling = cfg$coupling,
        multitraitTokens = .twasMultitraitTokens(norm$tokens, cfg$backend),
        engines = .twasRunEngines(norm$tokens, cfg$backend)
    )
    list(
        marker = marker,
        ctx = list(
            data = data,
            xRegions = cfg$xRegions,
            norm = norm,
            fineMappingResult = cfg$fineMappingResult,
            mashPriors = mash$mashPriors,
            cisWindow = cfg$cisWindow,
            residualizationParam = cfg$residualizationParam,
            genotypeFilter = cfg$genotypeFilter,
            phenotypeFilter = cfg$phenotypeFilter,
            completeSamples = cfg$completeSamples,
            twasWeights = cfg$twasWeights,
            verbose = cfg$verbose
        )
    )
}

# Resolve the (context, trait) grid the per-tuple fits run over.
# @noRd
.twasQdsResolveGrid <- function(data, contexts, traitId, region) {
    studyName <- studyName(data)
    useCtx <- .twasQdsResolveContexts(data, contexts)
    allTraits <- .twasQdsResolveTraits(data, useCtx, traitId, region)
    list(
        study = studyName,
        useCtx = useCtx,
        allTraits = allTraits,
        # Set only when traits are chosen context by context: the grid spans
        # every (context, trait) pair of allTraits, and this drops the pairs
        # the selection leaves out.
        traitSelection = if (.hasNestedTraitSelection(traitId)) traitId
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
.twasQdsResolveTraits <- function(
    data,
    useCtx,
    traitId,
    region,
    label = "twasWeightsPipeline(QtlDataset)"
) {
    perCtxTraits <- map(
        useCtx,
        .contextSelectedTraits,
        qd = data,
        traitId = traitId,
        region = region
    )
    allTraits <- unique(list_c(perCtxTraits))
    if (length(allTraits) == 0L) {
        abort(glue("{label}: no traits selected."))
    }
    allTraits
}

# Joint-pipeline marker carrying the CV / ensemble config for the engine.
# @noRd
.twasQdsMarker <- function(
    crossValidationParam,
    ensembleParam,
    standardized,
    fitRetention,
    retainX,
    estimatePi,
    h2,
    verbose,
    seed,
    coupling = NULL,
    multitraitTokens = NULL,
    engines = NULL
) {
    new(
        "TwasJointPipeline",
        config = list(
            # The two groups travel as records, not as nine loose fields:
            # this function's only job is to store them, so unrolling them
            # into its signature would name each one twice for nothing.
            crossValidationParam = crossValidationParam,
            ensembleParam = ensembleParam,
            fitRetention = fitRetention,
            retainX = retainX,
            # DERIVED from what the accessors actually do: genotypes reach
            # the fitter scaled to unit variance unless
            # `ResidualizationParam(scaleResiduals = FALSE)` says otherwise.
            # Resolved by the caller, which holds the bundle.
            standardized = standardized,
            estimatePi = estimatePi,
            # The run's `h2`, settled at entry (.h2PipelineSpec()).
            h2 = h2,
            verbose = verbose,
            seed = seed,
            ldReference = NULL,
            # The strengths a group across studies is coupled at
            # (.jointTwasCoupled()); NULL outside a MultiStudyQtlDataset.
            coupling = coupling,
            # The methods fitted jointly on a group of two or more
            # conditions as well as one condition at a time
            # (.jointTwasMultitrait()).
            multitraitTokens = multitraitTokens,
            # The engine of each method with several on individual-level
            # data (.twasRunEngines()), which says whether it fits a group
            # across conditions and on which layouts (.twasGroupedLayouts()).
            engines = engines
        )
    )
}

# The methods a run also fits jointly on each group a jointSpecification
# makes of two or more conditions: bayesC and bayesN run on BGLR (the run's
# `backend`), whose Multitrait() fits a group's conditions at once
# (.alphabetMultitraitMethods()). Their fits one condition at a time run as
# well.
# @noRd
.twasMultitraitTokens <- function(tokens, backend) {
    onBglr <- names(backend)[backend == "bglr"]
    as.character(intersect(
        intersect(tokens, .alphabetMultitraitMethods()),
        onBglr
    ))
}

# The methods a run also fits across each group of two or more conditions
# with a group penalty: those whose engine for the fit has one
# (.twasGroupedLayouts()). Their fits one condition at a time run as well.
# @noRd
.twasGroupedTokens <- function(tokens, backend) {
    keep(
        tokens,
        .twasGroups,
        path = "individual",
        cfg = list(engines = .twasRunEngines(tokens, backend))
    )
}

# The engine the run's `backend` puts each of `tokens` on for individual-level
# data, for the methods with several there (.methodBackends): a list keyed by
# method, holding only those it names. A method left out fits on its default
# engines: one condition at a time on its first, a group with the group
# penalty on the first that has one (.twasGroupEngine()).
# @noRd
.twasRunEngines <- function(tokens, backend) {
    several <- keep(tokens, .methodHasBackends, path = "individual")
    compact(set_names(map(several, .twasRunEngine, backend = backend), several))
}

# A choice of `backend` made for the summary-statistics path alone leaves the
# method on its defaults here.
# @noRd
.twasRunEngine <- function(token, backend, path = "individual") {
    engines <- names(.methodBackends[[path]][[token]])
    pluck(intersect(backend[names(backend) == token], engines), 1L)
}

# The engine the run's `backend` names for each of `tokens` on summary
# statistics, for the methods with several there: a list keyed by method,
# holding only those it names.
# @noRd
.twasSumstatEngines <- function(tokens, backend) {
    several <- keep(tokens, .methodHasBackends, path = "sumstat")
    compact(set_names(
        map(several, .twasRunEngine, backend = backend, path = "sumstat"),
        several
    ))
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
        scope = .twasQdsScope(studyName, useCtx, allTraits, traitSelection)
    )
}

# The engine's scope over one study's selected contexts and traits, with a
# per-context trait selection when the traits were chosen that way.
# @noRd
.twasQdsScope <- function(studyName, useCtx, allTraits, traitSelection) {
    list(
        studies = studyName,
        contexts = set_names(list(useCtx), studyName),
        traits = set_names(list(allTraits), studyName),
        traitSelection = traitSelection
    )
}

# Univariate methods ROUTED THROUGH THE ENGINE: one 1-condition group per
# (context, trait), per region -> the SAME per-method fitter (+ ensemble layer
# for >= 2 methods + resume cache) as the joint paths, merged across regions.
# `grid` is the dispatch row naming the study, contexts and traits.
# @noRd
.twasQdsUnivariateWork <- function(
    grid,
    marker,
    data,
    xRegions,
    norm,
    fineMappingResult,
    twasWeights,
    mashPriors,
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
    .qtlWorkOf(
        map(
            seq_along(xRegions),
            .twasQdsUnivRegionWork,
            cells = plan$cells,
            scope = plan$scope,
            marker = marker,
            data = data,
            xRegions = xRegions,
            norm = norm,
            fineMappingResult = fineMappingResult,
            twasWeights = twasWeights,
            mashPriors = mashPriors,
            cisWindow = cisWindow,
            phenotypeFilter = phenotypeFilter,
            residualizationParam = residualizationParam,
            genotypeFilter = genotypeFilter,
            nPCs = nPCs,
            verbose = verbose
        ),
        .twasUnivRegionsResult,
        list(labs = map_chr(xRegions, .twasRegionLabel))
    )
}

# @noRd
.twasUnivRegionsResult <- function(perRegion, labs) {
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
    mashPriors,
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
        mashPriors = mashPriors,
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
        backend = NULL,
        contexts = NULL,
        traitId = NULL,
        cisWindow = NULL,
        jointSpecification = NULL,
        fineMappingResult = NULL,
        twasWeights = NULL,
        mashPrior = NULL,
        h2Estimate = NULL,
        h2EstimateParam = NULL,
        fitRetention = c("slim", "none", "full"),
        verbose = 1L,
        panelFilterParam = PanelFilterParam(),
        crossValidationParam = CrossValidationParam(),
        ensembleParam = EnsembleParam(),
        coupling = list(study = c(0.01, 0.1, 1, 10)),
        residualCorrelation = NULL
    ) {
        fitRetention <- arg_match(fitRetention)
        .twasAssertCoupling(coupling)
        .twasAssertCouplingGroups(
            coupling,
            !missing(coupling),
            jointSpecification
        )
        .twasAssertResidualCorrelation(residualCorrelation)
        .twasWarnResidualOverride(mashPrior, residualCorrelation)
        .assertMethodParam(ensembleParam, "EnsembleParam")
        .ensembleAssertSumstats(ensembleParam, "QtlSumStats")
        chosen <- .twasSumStatsPreset(methods)
        .twasWarnUnretainableDetail(fitRetention, chosen, "QtlSumStats")
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
            methods = chosen,
            backend = backend,
            contexts = contexts,
            traitId = traitId,
            cisWindow = cisWindow,
            jointSpecification = jointSpecification,
            fineMappingResult = fineMappingResult,
            twasWeights = twasWeights,
            mashPrior = mashPrior,
            fitRetention = fitRetention,
            verbose = verbose,
            panelFilterParam = panelFilterParam,
            h2Estimate = h2Estimate,
            h2EstimateParam = h2EstimateParam,
            coupling = coupling,
            residualCorrelation = residualCorrelation,
            ensembleParam = ensembleParam
        )
        .withGenome(.twasRecordDataTypes(res, data), build)
    }
)

# ---- QtlSumStats pipeline worker + phase helpers ---------------------------

# Per-tuple rows for the QtlSumStats path. `part` carries the selected rows
# and tokens from .twasQssPerTupleSelection(); `cfg` the settings the
# builder takes. `twasWeights` is the resume cache, looked up one (study,
# context, trait, method) row at a time; the joint methods were fitted in the
# joint phase, which resumes a group only when every condition's row is
# cached from it.
# @noRd
.twasQssBuildRows <- function(part, methodArgs, cfg) {
    .twasQssUnivariateRows(
        part$selRows,
        studyCol = part$studyCol,
        contextCol = part$contextCol,
        traitCol = part$traitCol,
        methodArgs = methodArgs,
        univariateTokens = part$univariateTokens,
        data = cfg$data,
        twasWeights = cfg$twasWeights,
        fitRetention = cfg$fitRetention,
        fineMappingResult = cfg$fineMappingResult,
        panelFilterParam = cfg$panelFilterParam,
        h2Spec = cfg$h2Spec,
        ensemble = cfg$ensemble
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
    backend,
    contexts,
    traitId,
    cisWindow,
    jointSpecification,
    fineMappingResult,
    twasWeights,
    mashPrior = NULL,
    fitRetention,
    verbose,
    panelFilterParam,
    h2Estimate = NULL,
    h2EstimateParam = NULL,
    coupling = NULL,
    residualCorrelation = NULL,
    ensembleParam = NULL
) {
    inputs <- .twasQssInputs(data, contexts, traitId, cisWindow)
    .twasAssertResidualCovers(residualCorrelation, inputs)
    h2Spec <- .twasQssH2Spec(h2Estimate, h2EstimateParam, inputs)
    # One record of the settings the joint phase, the shared dispatcher and
    # the per-tuple loop all need, so each names them once.
    jointCfg <- list(
        data = inputs$data,
        unwindowed = inputs$unwindowed,
        backend = backend,
        contexts = inputs$contexts,
        traitId = inputs$traitId,
        verbose = verbose,
        fitRetention = fitRetention,
        panelFilterParam = panelFilterParam,
        fineMappingResult = fineMappingResult,
        twasWeights = twasWeights,
        mashPrior = mashPrior,
        coupling = coupling,
        residualCorrelation = residualCorrelation,
        ensemble = .ensembleResolve(ensembleParam),
        h2Spec = h2Spec
    )
    .twasQssRun(jointSpecification, methods, fineMappingResult, jointCfg)
}

# A QtlSumStats run's joint phase, then, unless it was the whole answer, the
# per-tuple fits of the methods it left.
# @noRd
.twasQssRun <- function(jointSpecification, methods, fineMappingResult, cfg) {
    joint <- .twasQssResolveTokens(
        jointSpecification,
        methods,
        fineMappingResult,
        cfg
    )
    if (joint$done) {
        return(joint$result)
    }
    .twasQssAfterJoint(joint, cfg, cfg$h2Spec)
}

# The per-tuple loop, after the joint phase has consumed the tokens it
# owns (.twasJointTokens()); `joint` carries the remainder, their options,
# and whatever the joint phase itself produced, and `cfg` is the run's joint
# record (.twasPipelineQtlSumStats()), whose settings it shares.
# @noRd
.twasQssAfterJoint <- function(joint, cfg, h2Spec = NULL) {
    part <- .twasQssPerTupleSelection(
        cfg$data,
        joint$tokens,
        cfg$contexts,
        cfg$traitId
    )
    rows <- .twasQssBuildRows(
        part,
        joint$methodArgs,
        list_assign(cfg, h2Spec = h2Spec)
    )
    .twasQssAssemble(rows, joint$result, cfg$data)
}

# Normalize the methods argument into (tokens, methodArgs). The default set
# excludes fine-mapping methods; those must be requested explicitly together
# with a fine-mapping result passed via `fineMappingResult`.
# @noRd
.twasSumStatsMethodTokens <- function(methods) {
    if (is.null(methods)) {
        tokens <- .twasSumStatsDefaults()
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

# The methods summary statistics run by default.
# @noRd
.twasSumStatsDefaults <- function() {
    c("lasso", "prsCs", "dprGibbs")
}

# Whether `methods` is a preset string (.twasExpandPresets()).
# @noRd
.twasIsPreset <- function(methods) {
    is.character(methods) &&
        length(methods) == 1L &&
        is_in(methods, c("default", "fastDefault"))
}

# The methods an individual-level run fits: those it names, a preset (no
# `methods` being the "default" one) for what it can run.
# @noRd
.twasRunMethods <- function(methods, fineMappingResult) {
    chosen <- methods %||% "default"
    if (.twasIsPreset(chosen)) {
        return(.twasPresetMethods(chosen, fineMappingResult))
    }
    chosen
}

# A preset's methods. Its fine-mapping methods read their fits from a
# fine-mapping result, so where none holds them they are skipped, which asking
# for them by name is not: that stops (.twasRequireFmResult()).
# @noRd
.twasPresetMethods <- function(preset, fineMappingResult) {
    expanded <- .twasExpandPresets(preset)
    skipped <- setdiff(
        intersect(expanded, .twasFineMappingTokens()),
        .twasFineMappingMethodsPresent(fineMappingResult)
    )
    if (length(skipped) > 0L) {
        warn(c(
            glue(
                "twasWeightsPipeline: the \"{preset}\" preset is skipping ",
                "{str_flatten(skipped, ', ')}, because no ",
                "`fineMappingResult` holds their fits."
            ),
            i = str_c(
                "To include them, run fineMappingPipeline() with these ",
                "methods and pass its result as `fineMappingResult`."
            )
        ))
    }
    setdiff(expanded, skipped)
}

# Summary statistics have one list of defaults, the one
# .twasSumStatsMethodTokens() takes NULL for, so a preset names that list.
# @noRd
.twasSumStatsPreset <- function(methods) {
    if (.twasIsPreset(methods)) NULL else methods
}

# A preset over a collection: each kind of component runs its own list.
# @noRd
.twasMsPreset <- function(methods, hasSumStats, fineMappingResult) {
    if (!.twasIsPreset(methods)) {
        return(methods)
    }
    forDatasets <- .twasPresetMethods(methods, fineMappingResult)
    if (!hasSumStats) {
        return(forDatasets)
    }
    TwasWeightsMethodsParam(
        qtlDatasetMethods = forDatasets,
        qtlSumStatsMethods = .twasSumStatsDefaults()
    )
}

# Empty-args map keyed by method token.
# @noRd
.twasEmptyMethodArgs <- function(tokens) {
    set_names(rep(list(list()), length(tokens)), tokens)
}

# The joint phase: every method fitted on joint groups -- mr.mash fitted,
# mvsusie read off its fine-mapping fits, and across studies the coupled and
# stacked ones -- over the groups the run's jointSpecification makes, as
# fine-mapping groups them; the joint methods refuse a run without one.
# Returns list(done, result, tokens, methodArgs):
# `done` requests an early return with `result`; otherwise the tokens left,
# all but the joint methods (.twasJointTokens()), and their args continue
# through the per-tuple loop.
# @noRd
.twasQssJointPhase <- function(parsedJointSpec, tokens, methodArgs, cfg) {
    jointTokens <- .twasJointTokens(tokens)
    .assertJointSpecForJointMethods(
        parsedJointSpec,
        jointTokens,
        "twasWeightsPipeline(QtlSumStats)"
    )
    jointResult <- if (length(parsedJointSpec) > 0L) {
        .twasDispatchJointSpecsQtlSumStats(
            parsedJointSpec,
            tokens,
            methodArgs,
            cfg
        )
    }
    if (length(jointTokens) == 0L) {
        return(list(
            done = FALSE,
            result = jointResult,
            tokens = tokens,
            methodArgs = methodArgs
        ))
    }
    keep <- setdiff(tokens, jointTokens)
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
    .ensembleAssertAcross(
        cfg$ensemble,
        parsedJointSpec,
        n_distinct(studyName(cfg$data)),
        "QtlSumStats"
    )
    .ensembleAssertAlleleFrequencies(cfg$ensemble, cfg$data, "QtlSumStats")
    tm <- .twasSumStatsMethodTokens(methods)
    .twasCheckMethodCapabilities(tm$tokens, "QtlSumStats")
    .twasCheckMethodArgsForInput(tm$methodArgs, "QtlSumStats")
    .twasCheckBackend(cfg$backend, tm$tokens, "sumstat")
    .twasCheckFineMappingMethods(tm$tokens, fineMappingResult, "QtlSumStats")
    .twasWarnUnusedMashPrior(cfg$mashPrior, tm$tokens)
    methodArgs <- .twasInjectTokenBackends(tm$methodArgs, cfg$backend)
    .twasQssJointPhase(
        parsedJointSpec,
        tm$tokens,
        methodArgs,
        list_assign(
            cfg,
            mashPriors = .twasQssMashPriors(cfg$mashPrior, tm$tokens)
        )
    )
}

# The pipeline's `backend` choices added to the summary-statistics methods'
# arguments, which are keyed by method.
# @noRd
.twasInjectTokenBackends <- function(methodArgs, backend) {
    if (is.null(backend)) {
        return(methodArgs)
    }
    reduce(
        names(backend),
        .twasInjectTokenBackend,
        backend = backend,
        .init = methodArgs
    )
}

# @noRd
.twasInjectTokenBackend <- function(methodArgs, token, backend) {
    methodArgs[[token]] <- list_assign(
        as.list(methodArgs[[token]] %||% list()),
        backend = backend[[token]]
    )
    methodArgs
}

# The per-tuple loop's selection: the selected rows and their identity
# columns, with the tokens the joint phase left, all univariate.
# @noRd
.twasQssPerTupleSelection <- function(data, tokens, contexts, traitId) {
    list(
        studyCol = as.character(data$study),
        contextCol = as.character(data$context),
        traitCol = as.character(data$trait),
        selRows = .twasQssSelectRows(data, contexts, traitId),
        univariateTokens = tokens
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

# The pipeline's fitRetention for a summary-statistics weight function, unless
# the caller set it for that method; every weight method takes one. "none" is
# passed on too, unlike on the individual-level path: the fine-mapping
# adapters here keep a slim fit unless told otherwise, so leaving it unset
# would not drop it.
# @noRd
.twasQssRetainArgs <- function(args, fitRetention) {
    list_assign(args, fitRetention = args$fitRetention %||% fitRetention)
}

# The retention level of `token`'s TWAS rows. A fine-mapping method's weights
# are extracted from the fineMappingResult's fit, so its rows keep none
# rather than a second copy of it; cTWAS reads it from the fine-mapping
# result.
# @noRd
.twasRowRetention <- function(token, fitRetention) {
    if (isTRUE(is_in(token, names(.twasFineMappingMethodAdapters)))) {
        return("none")
    }
    fitRetention
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
    fitRetention,
    fineMappingResult,
    methodArgs,
    panelFilterParam,
    h2Spec = NULL,
    ensemble = NULL
) {
    if (length(univariateTokens) == 0L) {
        return(list())
    }
    list_flatten(map(
        selRows,
        .twasQssUnivariateRowsForEntry,
        cols = list(studyCol, contextCol, traitCol),
        settings = list(
            univariateTokens = univariateTokens,
            data = data,
            twasWeights = twasWeights,
            fitRetention = fitRetention,
            fineMappingResult = fineMappingResult,
            methodArgs = methodArgs,
            panelFilterParam = panelFilterParam,
            h2Spec = h2Spec,
            ensemble = ensemble
        )
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
# With `settings$ensemble` enabled the entry's methods are stacked as well
# (.twasQssEnsembleRows()), which refits even the cached ones. `cols` are the
# rows' study, context and trait columns, and `settings` the run's.
# @noRd
.twasQssUnivariateRowsForEntry <- function(i, cols, settings) {
    key <- c(cols[[1L]][i], cols[[2L]][i], cols[[3L]][i])
    tokens <- settings$univariateTokens
    cached <- .twasQssCacheSplit(
        settings$twasWeights,
        key[[1L]],
        key[[2L]],
        key[[3L]],
        tokens
    )
    ldReference <- .rowLdRef(settings$data, i)
    cachedRows <- .recordsWithLdRef(unname(cached$rows), ldReference)
    stacking <- isTRUE(settings$ensemble$enabled)
    if (length(cached$toFit) == 0L && !stacking) {
        return(cachedRows)
    }
    fitCtx <- .twasQssEntryFitCtx(
        settings$data,
        key,
        ldReference,
        settings$panelFilterParam,
        settings$h2Spec,
        if (stacking) tokens else cached$toFit
    )
    fitted <- .twasQssFitTokens(cached$toFit, key, fitCtx, settings)
    .twasQssWithEnsemble(
        c(cachedRows, .recordsWithLdRef(fitted, ldReference)),
        fitCtx,
        settings$methodArgs,
        if (stacking) settings$ensemble,
        ldReference
    )
}

# Each of `tokens` fitted to one entry (.twasQssUnivariateFitOne()), for the
# (study, context, trait) `key`, under `settings`' method arguments,
# retention and fine-mapping result; those that fit nothing left out.
# @noRd
.twasQssFitTokens <- function(tokens, key, fitCtx, settings) {
    compact(map(
        tokens,
        .twasQssUnivariateFitOne,
        st = key[[1L]],
        ctx = key[[2L]],
        tr = key[[3L]],
        fitCtx = fitCtx,
        methodArgs = settings$methodArgs,
        fitRetention = settings$fitRetention,
        fineMappingResult = settings$fineMappingResult
    ))
}

# An entry's rows, with its ensemble's (.twasQssEnsembleRows()) when
# `ensemble` is given.
# @noRd
.twasQssWithEnsemble <- function(rows, fitCtx, methodArgs, ensemble, ref) {
    if (is.null(ensemble)) {
        return(rows)
    }
    c(
        rows,
        .recordsWithLdRef(
            .twasQssEnsembleRows(rows, fitCtx, methodArgs, ensemble),
            ref
        )
    )
}

# One entry's fit context (.twasQssUnivariateFitCtx()) for the (study,
# context, trait) `key`, with the heritability its methods `tokens` start
# from.
# @noRd
.twasQssEntryFitCtx <- function(
    data,
    key,
    ldReference,
    panelFilterParam,
    h2Spec,
    tokens
) {
    rowCtx <- .twasQssUnivariateFitCtx(
        data,
        key[[1L]],
        key[[2L]],
        key[[3L]],
        ldReference,
        cutoffs = .panelCutoffs(panelFilterParam)
    )
    c(
        rowCtx,
        list(
            h2 = .twasQssRowH2(
                h2Spec,
                key[[1L]],
                key[[2L]],
                key[[3L]],
                rowCtx,
                tokens
            )
        )
    )
}

# The SR-TWAS ensemble of one entry's methods on its summary statistics
# (.sumstatEnsemble()), as a one-row list, or an empty one when fewer than
# two of its methods clear the threshold. `rows` are the entry's methods'
# rows; a fine-mapping method's weights come from its fit, which a training
# part cannot refit, so it stays out.
# @noRd
.twasQssEnsembleRows <- function(rows, fitCtx, methodArgs, ensemble) {
    stacked <- compact(set_names(
        map(
            rows,
            .twasQssStackMethod,
            fitCtx = fitCtx,
            methodArgs = methodArgs
        ),
        map_chr(rows, "method")
    ))
    if (length(stacked) < 2L) {
        return(list())
    }
    stat <- fitCtx$stat
    ld <- .rssLdSource(stat, fitCtx$ldReference, "twasWeightsPipeline")
    ens <- .sumstatEnsemble(
        stacked,
        NULL,
        .lassosumClampCor(
            .rssMarginalCor(stat, n = stat$n, p = .ldSourceSize(ld))
        ),
        stat$n,
        ld,
        ensemble$r2Threshold
    )
    if (is.null(ens)) {
        return(list())
    }
    list(.twasRowRecord(
        rows[[1L]]$study,
        rows[[1L]]$context,
        rows[[1L]]$trait,
        "ensemble",
        .twasSumstatEnsembleRow(ens, fitCtx$variantIds)
    ))
}

# One row's method as the stack takes it: its full-data weights over the
# entry's variants and how to refit it (.twasQssRefit()). NULL for a
# fine-mapping method's row and for the ensemble's own.
# @noRd
.twasQssStackMethod <- function(row, fitCtx, methodArgs) {
    tk <- row$method
    spec <- .twasResolveWeightFn(tk)
    if (!is.null(spec$adapter) || is.null(spec$fn)) {
        return(NULL)
    }
    own <- set_names(weights(row$entry), variantIds(row$entry))
    userArgs <- .twasWithH2(
        .twasUserArgs(methodArgs, tk),
        tk,
        fitCtx$h2,
        "sumstat"
    )
    list(
        weights = replace_na(unname(own[fitCtx$variantIds]), 0),
        refit = list(
            fn = ".twasQssRefit",
            args = list(
                weightFn = spec$fn,
                userArgs = list_assign(userArgs, fitRetention = "none"),
                stat = fitCtx$stat,
                ldReference = fitCtx$ldReference,
                label = glue("twasWeightsPipeline: ensemble refit of '{tk}': ")
            )
        )
    )
}

# A method refit on a held-out split's training part of an entry's
# statistics: their correlations `trainCor` on `share` of the samples, given
# as the z-scores those correlations stand for (the inverse of .zToCor()), so
# that every weight function reads the same part. The weights come back over
# the entry's variants; a refit that fails gives none.
# @noRd
.twasQssRefit <- function(
    trainCor,
    share,
    weightFn,
    userArgs,
    stat,
    ldReference,
    label
) {
    n <- share * stat$n
    weights <- .twasTryWeights(
        weightFn,
        list_assign(
            stat,
            z = trainCor * sqrt((n - 2) / (1 - trainCor^2)),
            n = n
        ),
        ldReference,
        userArgs,
        label
    )
    if (is.null(weights)) {
        return(numeric(length(trainCor)))
    }
    as.numeric(`attr<-`(weights, "fit", NULL))
}

# The ensemble's row on summary statistics: its combined weights, and in
# place of a cross-validation result its coefficients, each method's
# held-out R^2 and its own held-out metrics.
# @noRd
.twasSumstatEnsembleRow <- function(ens, variantIds) {
    TwasWeightsRow(
        variantIds = variantIds,
        weights = ens$weights,
        cvResult = list(
            metrics = ens$metrics,
            methodCoef = ens$methodCoef,
            methodPerformance = ens$methodPerformance
        ),
        weightStandardized = TRUE
    )
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
    fitArg <- set_names(list(fit), spec$adapter$rssFitArg)
    list_assign(baseArgs, !!!fitArg)
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
    fineMappingResult
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
    withH2 <- .twasWithH2(userArgs, tk, fitCtx$h2, "sumstat")
    weights <- .twasTryWeights(
        spec$fn,
        fitCtx$stat,
        fitCtx$ldReference,
        .twasQssRetainArgs(withH2, .twasRowRetention(tk, fitRetention)),
        .twasFitErrUniv(tk, st, ctx, tr)
    )
    if (is.null(weights)) {
        return(NULL)
    }
    .twasQssRowFromWeights(weights, st, ctx, tr, tk, fitCtx)
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
    fitCtx
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
            weightStandardized = TRUE
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
    specs <- .msStudySpecs(cfg$studySpecs, qd)
    .muffleFitCoverage(twasWeightsPipeline(
        data = qd,
        jointSpecification = if (length(specs) > 0L) specs,
        fitRetention = cfg$fitRetention,
        retainX = cfg$retainX,
        methods = if (length(specs) > 0L) cfg$studyMethods else cfg$methods,
        backend = cfg$backend,
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
        ensembleParam = .ensembleWithinStudy(cfg$ensembleParam),
        seed = cfg$seed,
        genotypeFilterParam = cfg$genotypeFilterParam,
        mashPrior = cfg$mashPrior,
        estimatePi = cfg$estimatePi,
        # `[[` matches exactly, where `$` would take h2EstimateParam.
        h2Estimate = cfg[["h2Estimate"]],
        h2EstimateParam = cfg$h2EstimateParam,
        usePCA = cfg$usePCA,
        nPCs = cfg$nPCs,
        residualizationParam = cfg$residualizationParam,
        # A study with no specification of its own has no groups to couple.
        coupling = if (length(specs) > 0L) cfg$coupling,
        completeSamples = cfg$completeSamples
    ))
}

# Embedded-sumstats TWAS-weights worker for .multiStudyPipelineDriver.
# @noRd
.twasSumStats <- function(ss, cfg) {
    # Summary statistics carry no genotypes: each tuple's window is in the
    # input, narrowed by `cisWindow` around its trait position as on the
    # QtlDataset side, and there is nothing to residualize, cross-validate or
    # fit on full data. `cfg` carries those for the individual-level sibling
    # and they are not forwarded here. The ensemble is, as summary statistics
    # stack on held-out draws; the study axis, and with it `acrossStudies`
    # and the coupling `coupling` sets, is the collection's own run's, and
    # this run, with no specification, has no groups to couple.
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
        backend = cfg$backend,
        contexts = NULL,
        traitId = NULL,
        cisWindow = cfg$cisWindow,
        fineMappingResult = cfg$fineMappingResult,
        twasWeights = cfg$twasWeights,
        mashPrior = cfg$mashPrior,
        # `[[` matches exactly, where `$` would take h2EstimateParam.
        h2Estimate = cfg[["h2Estimate"]],
        h2EstimateParam = cfg$h2EstimateParam,
        verbose = cfg$verbose,
        panelFilterParam = cfg$panelFilterParam,
        ensembleParam = .ensembleWithinStudy(cfg$ensembleParam),
        coupling = NULL,
        residualCorrelation = cfg$residualCorrelation
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

# `retainX` is rejected unless the fit is "full", and when no requested method's
# fit has a copy of X to keep, as a "full" no method can honour is. A mixed run
# says nothing: a method with no copy loses nothing.
# @noRd
.twasCheckRetainX <- function(retainX, fitRetention, methods) {
    .assertRetainX(retainX, fitRetention)
    if (!retainX) {
        return(invisible(NULL))
    }
    tokens <- .twasMethodTokensFromArg(methods)
    holding <- .twasTokensKeepingX()
    if (length(tokens) > 0L && !any(is_in(tokens, holding))) {
        abort(glue(
            "twasWeightsPipeline: `retainX = TRUE` applies only to ",
            "{str_flatten(holding, ', ')}; the fits of ",
            "{str_flatten(tokens, ', ')} hold no copy of X."
        ))
    }
    invisible(NULL)
}

# The methods whose fit can hold a copy of X, by their weight functions' own
# formals.
# @noRd
.twasTokensKeepingX <- function() {
    fns <- compact(map(.twasMethodCapabilities, "individualImpl"))
    names(fns)[map_lgl(fns, .twasKeepsX)]
}

# The methods whose "full" fit keeps more than their "slim" one, on one input
# path. On summary statistics: mr.mash's whole fit and the shrinkage-grid
# methods' candidate coefficients (lassosumRssWeights() and the penalized
# RSS weights). On individual-level data: the engines' whole fitted objects,
# where "slim" keeps the selected fit or the posterior summaries (glmnet,
# ncvreg, L0Learn, the Bayesian alphabet), or drops more of the data (mr.ash);
# and mr.mash's whole fit. DPR holds nothing to drop but its copy of X, which
# goes at either level, so "full" loses it nothing.
# @noRd
.twasFullDetailTokens <- function(inputKind) {
    if (identical(inputKind, "QtlSumStats")) {
        return(c("mrmash", "lasso", "enet", "ridge", "scad", "mcp", "l0learn"))
    }
    c(
        "mrmash",
        "lasso",
        "enet",
        "ridge",
        "scad",
        "mcp",
        "l0learn",
        "mrash",
        "dprGibbs",
        "dprVb",
        "dprAdaptiveGibbs",
        "bayesA",
        "bayesB",
        "bayesC",
        "bayesL",
        "bayesN",
        "bayesR"
    )
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
        .twasExpandPresets(methods)
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

# A `methods` arg without the methods `drop`: mrmash, which the joint
# dispatcher fits, and fine-mapping methods read only off fits across
# studies, which leave the per-component recursion nothing to read.
# @noRd
.twasMsDropMethods <- function(methods, drop) {
    if (.twasMethodsIsParam(methods)) {
        return(.methodsParamDrop(methods, drop))
    }
    if (is.character(methods)) {
        setdiff(methods, drop)
    } else if (is.list(methods)) {
        methods[!is_in(.twasMethodTokensFromArg(methods), drop)]
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
    backend = NULL,
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
    h2Estimate = NULL,
    h2EstimateParam = NULL,
    usePCA = FALSE,
    nPCs = 10L,
    fitRetention = c("slim", "none", "full"),
    retainX = FALSE,
    phenotypeFilterParam = PhenotypeFilterParam(),
    verbose = 1,
    residualizationParam = ResidualizationParam(),
    crossValidationParam = CrossValidationParam(),
    ensembleParam = EnsembleParam(),
    coupling = list(study = c(0.01, 0.1, 1, 10)),
    completeSamples = FALSE,
    residualCorrelation = NULL,
    seed = NULL
) {
    fitRetention <- arg_match(fitRetention)
    .twasAssertCoupling(coupling)
    .twasAssertCouplingGroups(coupling, !missing(coupling), jointSpecification)
    assertFlag(completeSamples)
    .twasAssertResidualCorrelation(residualCorrelation)
    .twasWarnResidualOverride(mashPrior, residualCorrelation)
    .ensembleAssertCollection(ensembleParam, data)
    # Every collection holds individual-level studies, so their rule applies
    # here; a summary-statistic study is checked again by its own route.
    .twasWarnUnretainableDetail(fitRetention, methods, "QtlDataset")
    .twasCheckRetainX(retainX, fitRetention, methods)
    .msAssertH2(
        h2EstimateParam,
        data,
        "twasWeightsPipeline(MultiStudyQtlDataset)"
    )
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
    cfg <- list(
        methods = methods,
        backend = backend,
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
        h2Estimate = h2Estimate,
        h2EstimateParam = h2EstimateParam,
        usePCA = usePCA,
        nPCs = nPCs,
        fitRetention = fitRetention,
        retainX = retainX,
        residualizationParam = residualizationParam,
        coupling = coupling,
        completeSamples = completeSamples,
        residualCorrelation = residualCorrelation
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
        return(.withGenome(.twasRecordDataTypes(joint$result, data), build))
    }
    .withGenome(
        .twasRecordDataTypes(
            .twasMsDriver(
                data,
                joint$result,
                list_assign(
                    cfg,
                    methods = joint$methods,
                    studySpecs = joint$studySpecs,
                    studyMethods = joint$studyMethods
                )
            ),
            data
        ),
        build
    )
}

# The run's coupling: NULL for none, or a list naming the axes it couples
# ("study", "context", "trait"), each with its grid of strengths, distinct
# and not negative.
# @noRd
.twasAssertCoupling <- function(coupling) {
    if (is.null(coupling)) {
        return(invisible(NULL))
    }
    assertList(coupling, min.len = 1L)
    assertNames(
        names(coupling),
        type = "unique",
        subset.of = c("study", "context", "trait")
    )
    walk(coupling, .twasAssertCouplingGrid)
}

# A `coupling` the caller names couples a group's conditions, and a run with
# no jointSpecification makes no groups, so it would go unused. `supplied` is
# whether the caller named it: the default couples studies and is left
# unused without complaint.
# @noRd
.twasAssertCouplingGroups <- function(coupling, supplied, jointSpecification) {
    if (!supplied || is.null(coupling) || !is.null(jointSpecification)) {
        return(invisible(NULL))
    }
    abort(c(
        str_c(
            "twasWeightsPipeline: `coupling` couples the conditions of the ",
            "groups a `jointSpecification` makes, and this run has none."
        ),
        i = "Pass a `jointSpecification`, or leave `coupling` unset."
    ))
}

# The run's residual correlation for the summary-statistics held-out splits:
# NULL, or a correlation matrix whose rows and columns are named by condition,
# "study:context:trait".
# @noRd
.twasAssertResidualCorrelation <- function(residualCorrelation) {
    if (is.null(residualCorrelation)) {
        return(invisible(NULL))
    }
    assertMatrix(
        residualCorrelation,
        mode = "numeric",
        any.missing = FALSE,
        nrows = ncol(residualCorrelation),
        row.names = "unique",
        col.names = "unique"
    )
    ok <- identical(
        rownames(residualCorrelation),
        colnames(residualCorrelation)
    ) &&
        isSymmetric(unname(residualCorrelation)) &&
        all(abs(diag(residualCorrelation) - 1) < 1e-8) &&
        all(abs(residualCorrelation) <= 1 + 1e-8)
    if (!ok) {
        abort(str_c(
            "`residualCorrelation` must be a correlation matrix: symmetric, ",
            "unit diagonal, its rows and columns named alike by condition ",
            "(\"study:context:trait\")."
        ))
    }
    invisible(NULL)
}

# The full-data priors of a MashPrior, whose residual correlations the
# summary-statistics held-out splits take when no `residualCorrelation` is
# given; NULL without one, or with a `residualCorrelation`, which overrides
# it (.twasWarnResidualOverride()).
# @noRd
.twasResidualPriors <- function(mashPrior, residualCorrelation) {
    if (is.null(mashPrior) || !is.null(residualCorrelation)) {
        return(NULL)
    }
    .unpackMashPrior(mashPrior)$priors$full
}

# `residualCorrelation` must name every selected condition of the run.
# @noRd
.twasAssertResidualCovers <- function(residualCorrelation, inputs) {
    if (is.null(residualCorrelation)) {
        return(invisible(NULL))
    }
    data <- inputs$data
    rows <- which(.sumStatsRowsSelected(data, inputs$contexts, inputs$traitId))
    keys <- str_c(data$study, data$context, data$trait, sep = ":")[rows]
    missing <- setdiff(keys, rownames(residualCorrelation))
    if (length(missing) > 0L) {
        abort(glue(
            "`residualCorrelation` has no row for ",
            "{str_flatten(missing, ', ')}."
        ))
    }
    invisible(NULL)
}

# A warning, at entry, that `residualCorrelation` overrides the MashPrior's.
# @noRd
.twasWarnResidualOverride <- function(mashPrior, residualCorrelation) {
    if (is.null(mashPrior) || is.null(residualCorrelation)) {
        return(invisible(NULL))
    }
    warn(str_c(
        "twasWeightsPipeline: `residualCorrelation` overrides the ",
        "`mashPrior`'s residual correlation for the held-out splits."
    ))
}

# @noRd
.twasAssertCouplingGrid <- function(strengths) {
    assertNumeric(
        strengths,
        lower = 0,
        finite = TRUE,
        any.missing = FALSE,
        min.len = 1L,
        unique = TRUE
    )
}

# The MultiStudyQtlDataset joint record, derived from the record the
# per-component recursion receives (`cfg`, which .twasPipelineMultiStudy()
# builds): the settings each
# component's own QtlDataset run fits under (.twasJointCfg()), plus the
# collection, the `region` its traits are selected by, the two filters the
# components apply, and the coupling strengths of the fits across studies.
# Deriving it, rather than re-listing the settings, keeps the per-study joint
# fits on exactly the settings of the per-study fits.
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
            h2 = .twasMsH2Spec(cfg, data),
            fitRetention = cfg$fitRetention,
            retainX = cfg$retainX,
            residualizationParam = cfg$residualizationParam,
            verbose = cfg$verbose,
            seed = cfg$seed,
            twasWeights = cfg$twasWeights,
            backend = cfg$backend,
            completeSamples = cfg$completeSamples
        ),
        data = data,
        region = cfg$region,
        genotypeFilterParam = cfg$genotypeFilterParam,
        panelFilterParam = cfg$panelFilterParam,
        coupling = cfg$coupling
    )
}

# The run's heritability for the collection's joint fits, settled at entry
# (.h2PipelineSpec()) against every study's selected (study, context, trait)
# as each study's own run settles its own (.twasQdsH2Spec()); settings
# estimate each trait's in-sample.
# @noRd
.twasMsH2Spec <- function(cfg, data) {
    studies <- keep(
        as.list(data),
        .studyHasSelectedTraits,
        contexts = cfg$contexts,
        traitId = cfg$traitId,
        region = cfg$region
    )
    .h2PipelineSpec(
        cfg[["h2Estimate"]],
        cfg$h2EstimateParam,
        list_rbind(map(studies, .twasMsStudyTuples, cfg = cfg)),
        "twasWeightsPipeline(MultiStudyQtlDataset)",
        individual = TRUE
    )
}

# One study's selected (study, context, trait), its share of the selection.
# @noRd
.twasMsStudyTuples <- function(qd, cfg) {
    .twasQdsTuples(
        qd,
        cfg$contexts,
        .studyTraitIdArg(cfg$traitId, studyName(qd)),
        cfg$region
    )
}

#' @rdname twasWeightsPipeline
#' @export
setMethod(
    "twasWeightsPipeline",
    "MultiStudyQtlDataset",
    .twasPipelineMultiStudy
)

# Joint-specification dispatch for the joint methods (mr.mash, mvsusie,
# fsusie), grouped as fineMappingPipeline() groups them. Returns
# list(done, result, methods) where `done` requests an early return with
# `result` and `methods` is the set left for the per-component recursion.
# The joint methods refuse a run without a spec, which fits every other
# method one condition at a time in each component.
# @noRd
.twasMsJointPhase <- function(parsedJointSpec, methods, cfg) {
    .assertJointSpecForJointMethods(
        parsedJointSpec,
        .twasJointTokens(.twasMethodTokensFromArg(methods)),
        "twasWeightsPipeline(MultiStudyQtlDataset)"
    )
    if (length(parsedJointSpec) == 0L) {
        return(list(done = FALSE, result = NULL, methods = methods))
    }
    # A component with per-tuple fits to come fits its joint specs along with
    # them, on one read of each window. A method fitted jointly as well as
    # one tuple at a time (.twasMultitraitTokens(), .twasGroupedTokens())
    # stays among those.
    tokens <- .twasMethodTokensFromArg(methods)
    jointTokens <- .twasJointTokens(tokens)
    stripped <- .twasMsDropMethods(methods, jointTokens)
    acrossConditions <- c(
        jointTokens,
        .twasMultitraitTokens(tokens, cfg$backend),
        .twasGroupedTokens(tokens, cfg$backend)
    )
    withComponents <- !.twasMethodsEmpty(stripped) &&
        length(acrossConditions) > 0L
    jointResult <- .twasDispatchJointSpecsMultiStudy(
        parsedJointSpec,
        tokens,
        methods,
        cfg,
        components = !withComponents
    )
    if (.twasMethodsEmpty(stripped)) {
        return(list(done = TRUE, result = .twasMsRequireJoint(jointResult)))
    }
    list(
        done = FALSE,
        result = jointResult,
        methods = stripped,
        studySpecs = if (withComponents) {
            .fmSplitStudyAxisSpecs(parsedJointSpec)$nonStudy
        },
        studyMethods = methods
    )
}

# The joint phase's result when it leaves the components nothing to run,
# which then has to hold something.
# @noRd
.twasMsRequireJoint <- function(jointResult) {
    if (is.null(jointResult)) {
        msg <- glue(
            "twasWeightsPipeline(MultiStudyQtlDataset): no joint fits ",
            "produced."
        )
        abort(msg)
    }
    jointResult
}

# Gate fine-mapping methods early so the recursion into the embedded
# QtlDataset / QtlSumStats components doesn't re-run fine-mapping, then hand
# off to the collection's own phase. Mirrors .fmMsResolveTokens.
# @noRd
.twasMsResolveTokens <- function(
    jointSpecification,
    methods,
    fineMappingResult,
    cfg
) {
    parsedJointSpec <- parseJointSpecification(jointSpecification, cfg$data)
    .ensembleAssertAcross(
        cfg$ensembleParam,
        parsedJointSpec,
        .msStudyCount(cfg$data),
        "MultiStudyQtlDataset"
    )
    # Translated here, where the dataset is in hand, so the Param is what
    # travels into the per-component recursion; each component reads the
    # slot for its own input class.
    hasSumStats <- !is.null(sumStats(cfg$data))
    methods <- .methodsParamForMulti(
        .twasMsPreset(methods, hasSumStats, fineMappingResult),
        "TwasWeightsMethodsParam",
        "twasWeightsPipeline",
        hasSumStats
    )
    .twasCheckFineMappingMethods(
        .twasMethodTokensFromArg(methods),
        fineMappingResult,
        "MultiStudyQtlDataset",
        fitFullData = .cvResolve(cfg$crossValidationParam)$fitFullData
    )
    .twasMsJointPhase(parsedJointSpec, methods, cfg)
}

# ---- weights read off fine-mapping fits across studies ----------------------
#
# fineMappingPipeline()'s study axis fits one model across a collection's
# individual-level studies and stores the whole fit on each study's row. A run
# whose jointSpecification has the study axis groups the studies as that run
# did (.twasStudyAxisGroups()) and, as nothing here fine-maps, reads each
# group's weights off the fit over exactly its studies: one TWAS row per
# fine-mapping row, holding that row's column of the fit's weights. No
# genotype is read.

# One TWAS row per row of each mvSuSiE method's (`tokens`) fit over each of
# the study-axis groups `groups` (condition tables), stopping on a group the
# fine-mapping result holds no fit over exactly (.twasRequireFmGroup()); NULL
# for no rows. The rows are built by the joint engine's record builder, as a
# joint fit's are.
# @noRd
.twasAcrossStudyWeights <- function(groups, cfg, tokens) {
    fm <- cfg$fineMappingResult
    rows <- list_c(
        map(tokens, .twasStudyTokenRows, groups = groups, fm = fm),
        ptype = integer()
    )
    read <- list(
        crossValidationParam = .cvResolve(cfg$crossValidationParam),
        residualizationParam = cfg$residualizationParam
    )
    records <- map(rows, .twasAcrossStudyRecord, fm = fm, read = read)
    .buildJointResult(TwasWeights, compact(records))
}

# The rows of `token`'s fits over the study-axis groups `groups`.
# @noRd
.twasStudyTokenRows <- function(token, groups, fm) {
    list_c(
        map(groups, .twasStudyGroupRows, fm = fm, token = token),
        ptype = integer()
    )
}

# The rows of `token`'s fit over exactly group `cond`'s studies.
# @noRd
.twasStudyGroupRows <- function(cond, fm, token) {
    joint <- .jointKeyOfConditions(cond)
    .twasRequireFmGroup(fm, cond, joint, token)
    .twasFmConditionRows(fm, token, cond, joint)
}

# The groups the study-axis specs `studySpecs` make of the individual-level
# studies `studies`, made as fineMappingPipeline() makes them: each spec's
# scope over each region block, by its cross-study enumerator, read only as
# far as its groups' conditions (their phenotypes; no genotype). One of each.
# @noRd
.twasStudyAxisGroups <- function(studySpecs, studies, cfg) {
    groups <- list_flatten(map(
        studySpecs,
        .twasStudySpecGroups,
        studies = studies,
        cfg = cfg
    ))
    unname(groups[!duplicated(map_chr(groups, .twasConditionsKey))])
}

# @noRd
.twasStudySpecGroups <- function(spec, studies, cfg) {
    list_flatten(map(
        cfg$xRegions,
        .twasStudySpecRegionGroups,
        spec = spec,
        studies = studies,
        cfg = cfg
    ))
}

# One spec's groups over one region block (.jointSpecWork()'s scope and
# enumeration, under the arguments the fine-mapping enumeration reads).
# @noRd
.twasStudySpecRegionGroups <- function(region, spec, studies, cfg) {
    scope <- .jointSpecScope(
        spec,
        studies,
        "individual",
        cfg$contexts,
        cfg$traitId,
        region
    )
    cell <- .lookupJointCell(.jointSpecPattern(spec$axes), "individual")
    units <- .jcEnumerate(cell)(
        studies,
        scope,
        list(
            axes = spec$axes,
            region = region,
            cisWindow = cfg$cisWindow,
            residualizationParam = cfg$residualizationParam,
            phenotypeFilter = cfg$phenotypeFilter,
            completeSamples = cfg$completeSamples,
            verbose = cfg$verbose
        )
    )
    map(units, .xsUnitConditions)
}

# A condition table's key, telling one group from another.
# @noRd
.twasConditionsKey <- function(cond) {
    str_flatten(str_c(cond$study, cond$context, cond$trait, sep = "\r"), "\n")
}

# The TWAS record of fine-mapping row `i`, in the joint engine's record shape
# (.jointEntryRecordAt()): the row's identity, joint columns and windows, and
# the entry read off its fit. NULL when the row holds neither a fit nor a
# cross-validation.
# @noRd
.twasAcrossStudyRecord <- function(i, fm, read) {
    md <- mcols(fm, use.names = FALSE)
    parts <- .fmrRowParts(fm, i)
    row <- list(
        token = .twasFmTokenFor(as.character(md$method)[[i]]),
        study = as.character(md$study)[[i]],
        column = .fmrJointIndexAt(fm, i),
        stdz = .twasFmRowStandardized(md, i, read)
    )
    entry <- .twasAcrossStudyEntry(susieFit(parts), cvResult(parts), row, read)
    if (is.null(entry)) {
        return(NULL)
    }
    list(
        study = row$study,
        context = as.character(md$context)[[i]],
        trait = as.character(md$trait)[[i]],
        method = row$token,
        entry = entry,
        jointStudies = as.character(md$jointStudies)[[i]],
        jointContexts = .twasFmRowValue(md, "jointContexts", i, NA_character_),
        jointTraits = .twasFmRowValue(md, "jointTraits", i, NA_character_),
        traitPos = .twasFmRowValue(md, "traitPos", i),
        region = .twasFmRowValue(md, "region", i),
        jointIndex = .fmrJointIndexAt(fm, i)
    )
}

# Row `i`'s value in column `col` of fine-mapping columns `md`, or `absent`
# when there is no such column.
# @noRd
.twasFmRowValue <- function(md, col, i, absent = NULL) {
    v <- md[[col]]
    if (is.null(v)) absent else v[i]
}

# Whether row `i`'s weights are on the standardized genotype scale: as its
# fine-mapping run recorded (.fmRecordProvenance()), else as this run's
# residualization scales the genotypes.
# @noRd
.twasFmRowStandardized <- function(md, i, read) {
    recorded <- .twasFmRowValue(md, "standardized", i, NA)
    if (is.na(recorded)) {
        return(.qtlScaleResiduals(read$residualizationParam))
    }
    recorded
}

# A row's TwasWeightsRow from its stored fit `fits` and its cross-validation
# `cvs`. A row fitted over several ranges keeps one of each per range, and
# its pieces are merged as the per-tuple path merges ranges
# (.twasMergeRegionEntries()). `row` names the method, the study and the
# weights' scale.
# @noRd
.twasAcrossStudyEntry <- function(fits, cvs, row, read) {
    if (!.twasIsRegionList(fits) && !.twasIsRegionList(cvs)) {
        return(.twasAcrossStudyPiece(fits, cvs, row, read))
    }
    labels <- names(fits) %||% names(cvs)
    pieces <- map(
        labels,
        .twasAcrossStudyPieceAt,
        fits = fits,
        cvs = cvs,
        row = row,
        read = read
    )
    .twasMergeRegionEntries(pieces, labels)
}

# @noRd
.twasAcrossStudyPieceAt <- function(label, fits, cvs, row, read) {
    .twasAcrossStudyPiece(fits[[label]], cvs[[label]], row, read)
}

# One range's TwasWeightsRow (.jointTwasEntry()): the row's column of the
# fit's weights, unless the run is cross-validation only, and the range's
# cross-validation when the run asks for one. Like every fine-mapping
# method's row it keeps no copy of the fit; cTWAS reads that off the
# fineMappingResult.
# @noRd
.twasAcrossStudyPiece <- function(fit, cv, row, read) {
    if (is.null(fit) && is.null(cv)) {
        return(NULL)
    }
    W <- if (.jointFitsFullData(read)) {
        .twasFmFitColumn(fit, row$token, row$column, row$study)
    } else {
        NULL
    }
    .jointTwasEntry(
        1L,
        1L,
        list(W = W, fitParts = NULL, vids = rownames(W)),
        .twasAcrossStudyCv(cv, row$token, W, read),
        row$stdz,
        read
    )
}

# Column `column` -- the row's position in the fit's condition axis -- of the
# weights the method's own extractor reads off `fit`, as a one-column matrix
# over the fit's variants. A row merged over ranges whose fits dropped
# different conditions records no position, and has no column to read.
# @noRd
.twasFmFitColumn <- function(fit, token, column, study) {
    if (is.na(column)) {
        msg <- glue(
            "twasWeightsPipeline: study '{study}''s row of a {token} fit ",
            "across studies records no position in the fit: its ranges' ",
            "fits dropped different studies, so its column is not known."
        )
        abort(msg)
    }
    adapter <- .twasFineMappingMethodAdapters[[token]]
    fitArg <- set_names(list(fit), adapter$fitArg)
    W <- exec(adapter$weightFn, !!!fitArg)
    W[, column, drop = FALSE]
}

# A range's cross-validation when the run cross-validates: the fine-mapping
# run's own (.twasFmHandoffCv()), on its folds, which a partition given here
# must match. Without one the run is refused, unless this method's
# cross-validation was not asked for (.jointTwasCvBlocked()).
# @noRd
.twasAcrossStudyCv <- function(cv, token, W, read) {
    cvParam <- read$crossValidationParam
    if ((cvParam$folds %||% 0L) <= 1L) {
        return(NULL)
    }
    handoff <- .twasFmHandoffCv(cv, token)
    if (is.null(handoff)) {
        .jointTwasCvBlocked(W, read, token)
        return(NULL)
    }
    .jointCheckPartitionAgreement(
        handoff$samplePartition,
        cvParam$samplePartition
    )
    handoff
}

# Run the per-study / per-component recursion via the shared multi-study driver.
# `cfg` is the record of the settings each component re-runs with, which
# .twasPipelineMultiStudy() builds and the driver hands to .twasPerStudy /
# .twasSumStats.
# @noRd
.twasMsDriver <- function(data, jointResult, cfg) {
    .multiStudyPipelineDriver(
        data,
        jointResult,
        .twasPerStudy,
        .twasSumStats,
        cfg,
        .rbindTwasWeights,
        "twasWeightsPipeline",
        noun = "weights"
    )
}

# ---- the heritability a run's methods start from ----------------------------

# The methods that can start from a trait's heritability: each takes it as a
# value (`h2Estimate`), and LDpred2 also as the settings to estimate it under
# on its own terms (`h2EstimateParam`).
# @noRd
.twasH2Methods <- c(
    "bayesA",
    "bayesB",
    "bayesC",
    "bayesL",
    "bayesN",
    "bayesR",
    "ldpred2Auto",
    "ldpred2Grid"
)

# A QtlDataset run's heritability, settled at entry (.h2PipelineSpec())
# against the (study, context, trait) it selects; settings estimate each
# trait's in-sample.
# @noRd
.twasQdsH2Spec <- function(
    h2Estimate,
    h2EstimateParam,
    data,
    contexts,
    traitId,
    region
) {
    .h2PipelineSpec(
        h2Estimate,
        h2EstimateParam,
        .twasQdsTuples(data, contexts, traitId, region),
        "twasWeightsPipeline(QtlDataset)",
        individual = TRUE
    )
}

# The (study, context, trait) a QtlDataset run selects.
# @noRd
.twasQdsTuples <- function(data, contexts, traitId, region) {
    list_rbind(map(
        .twasQdsResolveContexts(data, contexts),
        .twasQdsContextTuples,
        data = data,
        traitId = traitId,
        region = region
    ))
}

# @noRd
.twasQdsContextTuples <- function(cx, data, traitId, region) {
    tibble(
        study = studyName(data),
        context = cx,
        trait = .contextSelectedTraits(cx, data, traitId, region)
    )
}

# A QtlSumStats run's heritability, settled at entry against its selected
# rows; settings that read the LD reference's sample size need each row's
# reference to record one.
# @noRd
.twasQssH2Spec <- function(h2Estimate, h2EstimateParam, inputs) {
    if (is.null(h2Estimate) && is.null(h2EstimateParam)) {
        return(NULL)
    }
    data <- inputs$data
    rows <- which(.sumStatsRowsSelected(data, inputs$contexts, inputs$traitId))
    label <- "twasWeightsPipeline(QtlSumStats)"
    spec <- .h2PipelineSpec(
        h2Estimate,
        h2EstimateParam,
        .h2RowTuples(data, rows),
        label
    )
    walk(
        rows,
        .h2AssertRowReference,
        data = data,
        param = spec$param,
        label = label
    )
    spec
}

# What one QtlSumStats row's trait gets from the run's heritability
# (.twasRowH2()), an estimate being made on the row's statistics and LD
# reference.
# @noRd
.twasQssRowH2 <- function(spec, st, ctx, tr, rowCtx, tokens) {
    label <- glue(
        "twasWeightsPipeline(QtlSumStats): study='{st}', context='{ctx}', ",
        "trait='{tr}'"
    )
    .twasRowH2(
        spec,
        .h2TupleKey(st, ctx, tr),
        label,
        .h2StatFit,
        list(
            stat = rowCtx$stat,
            ldReference = rowCtx$ldReference,
            label = label
        ),
        tokens
    )
}

# What a trait's methods get from the run's heritability: NULL when it sets
# none, or when none of `tokens` starts from one; else its settings, for
# LDpred2 to estimate under on its own terms, and the trait's value
# (.h2ForTrait()), made only when a method among `tokens` takes a number.
# @noRd
.twasRowH2 <- function(spec, key, label, fit, fitArgs, tokens) {
    if (is.null(spec) || !any(is_in(tokens, .twasH2Methods))) {
        return(NULL)
    }
    takesValue <- is.null(spec$param) ||
        any(is_in(tokens, names(.alphabetBackends$individual)))
    list(
        param = spec$param,
        value = if (takesValue) .h2ForTrait(spec, key, label, fit, fitArgs)
    )
}

# A method's arguments with its trait's heritability from the run
# (`rowH2`, .twasRowH2()), unless they set their own.
# @noRd
.twasWithH2 <- function(userArgs, tk, rowH2, path) {
    if (is.null(rowH2) || !is_in(tk, .twasH2Methods)) {
        return(userArgs)
    }
    given <- .twasH2Given(tk, rowH2)
    if (is.null(given) || .twasSetsOwnH2(userArgs, tk, path)) {
        return(userArgs)
    }
    list_assign(userArgs, !!!given)
}

# The argument a method takes the run's heritability as: for LDpred2 the
# run's settings, where it has them; otherwise the trait's value (LDpred2-grid
# centring its grid on it); NULL without one.
# @noRd
.twasH2Given <- function(tk, rowH2) {
    ldpred2 <- is_in(tk, c("ldpred2Auto", "ldpred2Grid"))
    if (ldpred2 && !is.null(rowH2$param)) {
        return(list(h2EstimateParam = rowH2$param))
    }
    if (is.null(rowH2$value)) {
        return(NULL)
    }
    list(h2Estimate = rowH2$value)
}

# Whether a method's own arguments set its heritability: one of its own
# arguments for it, or for the alphabet its backend's spelling of it.
# @noRd
.twasSetsOwnH2 <- function(userArgs, tk, path) {
    given <- c(
        names(userArgs),
        names(as.list(userArgs$methodArgs %||% list()))
    )
    engineNames <- if (is_in(tk, names(.alphabetBackends[[path]]))) {
        .alphabetH2Names[[.methodBackend(tk, userArgs$backend, path)]]
    }
    own <- c("h2Estimate", "h2EstimateParam", "h2Grid", engineNames)
    any(is_in(own, given))
}

# REML reads a trait's genotypes and phenotype, so a collection holding
# summary statistics cannot estimate under it on all of its studies.
# @noRd
.msAssertH2 <- function(h2EstimateParam, data, label) {
    reml <- is(h2EstimateParam, "H2EstimateParam") &&
        identical(h2EstimateParam@statistic, "phenotype")
    if (!reml || is.null(sumStats(data))) {
        return(invisible(NULL))
    }
    abort(c(
        glue(
            "{label}: the phenotype statistic (REML) needs individual-level ",
            "data, and the collection holds summary statistics."
        ),
        i = str_c(
            "Use summary-statistic settings, or estimate every study's with ",
            "estimateLocalH2() under one setting and pass them combined with ",
            "c() as `h2Estimate`."
        )
    ))
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
    .solveEnsembleMoments(
        crossprod(Pvalid),
        as.vector(crossprod(Pvalid, yObs)),
        Kvalid
    )
}

# The quadprog stack from its moments: the predictions' Gram matrix `gram`
# (P'P) and their products with the outcome `dvec` (P'y), which summary
# statistics give without the predictions themselves.
# @noRd
.solveEnsembleMoments <- function(gram, dvec, Kvalid) {
    if (!requireNamespace("quadprog", quietly = TRUE)) {
        abort("Package 'quadprog' is required for solver='quadprog'.")
    }

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
#'   keys like \code{"susie_predicted"}) and, to cross-validate the ensemble
#'   itself, the \code{$samplePartition} the predictions were made over. For
#'   multi-dataset: a list of such objects.
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
#'   \item{predictions}{The ensemble's own out-of-fold predictions, cross-fit
#'     (see Details) and named by sample, or \code{NULL} unless every sample
#'     is named and placed in a fold of \code{cvResults$samplePartition}.}
#'   \item{metrics}{The CV metrics of those predictions, as
#'     \code{\link{twasWeightsCv}} reports each method's: \code{corr},
#'     \code{rsq}, \code{adj_rsq}, \code{pval}, \code{RMSE} and \code{MAE}.
#'     \code{NULL} with \code{predictions}.}
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
#' The ensemble is cross-validated over the same folds as its methods, by
#' cross-fitting: each fold's samples are predicted with coefficients fit on
#' the other folds' out-of-fold predictions only, so its metrics compare with
#' the methods' own. The coefficients returned in \code{methodCoef} are the
#' ones fit on every sample.
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
    crossFit <- .ensembleCrossFit(stacked, solver, alpha)
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
        ),
        predictions = crossFit$predictions,
        metrics = crossFit$metrics
    )
}

# The ensemble's own out-of-fold predictions, cross-fit -- each fold's
# samples predicted with coefficients fit on the other folds' out-of-fold
# predictions only -- named by sample, and their CV metrics. Both NULL unless
# every sample has a fold.
# @noRd
.ensembleCrossFit <- function(stacked, solver, alpha) {
    if (is.null(stacked$fold) || n_distinct(stacked$fold) < 2L) {
        return(list(predictions = NULL, metrics = NULL))
    }
    parts <- map(
        unique(stacked$fold),
        .ensembleFoldPrediction,
        stacked = stacked,
        solver = solver,
        alpha = alpha
    )
    rows <- list_c(map(parts, "rows"))
    pred <- set_names(
        list_c(map(parts, "pred"))[order(rows)],
        stacked$samples
    )
    list(predictions = pred, metrics = .cvMetricRow(pred, stacked$yObs))
}

# Fold `f`'s rows of the stacked predictions, combined with the coefficients
# the other folds' rows give.
# @noRd
.ensembleFoldPrediction <- function(f, stacked, solver, alpha) {
    held <- stacked$fold == f
    zeta <- .ensembleTrainZeta(
        stacked$P[!held, , drop = FALSE],
        stacked$yObs[!held],
        solver,
        alpha
    )
    list(
        rows = which(held),
        pred = drop(stacked$P[held, , drop = FALSE] %*% zeta)
    )
}

# Coefficients over every method from training rows `P`: solved over the
# methods whose predictions vary there, the rest left at 0. A lone such method
# takes full weight, as .ensembleSolveZeta() gives it, without its message.
# @noRd
.ensembleTrainZeta <- function(P, y, solver, alpha) {
    valid <- apply(P, 2, sd) > .Machine$double.eps
    if (sum(valid) <= 1L) {
        return(as.numeric(valid))
    }
    replace(
        rep(0, ncol(P)),
        valid,
        .ensembleSolveValid(P[, valid, drop = FALSE], y, solver, alpha)
    )
}

# ---- The ensemble on summary statistics -------------------------------------
#
# Summary statistics hold no samples to cross-validate on, so the stack reads
# PUMAS splits of the statistics instead, as held-out selection draws them
# (.rssSelectHeldOut()): each split draws a training part's marginal
# correlations from the LD, refits every method on it, and scores the refits
# on the held-out remainder. Over a split's predictors B, the held-out part
# gives the moments the individual-level stack reads off its predictions:
# B' Sigma B for P'P and B' c_h for P'y, Sigma the LD shrunk toward the
# identity as held-out selection shrinks it. Summed over the splits they make
# one quadprog stack, whose coefficients combine the full-data weights.

# The stack over the methods `own` -- each a list of its full-data `weights`
# over the statistics' variants and its `refit`, the function (`fn`) and
# other arguments (`args`) that refit it on a training part's correlations
# `trainCor` holding the sample share `share` -- and the fixed predictors
# `external`, a variants x predictors matrix (possibly with no column).
# `cor`, `n` and `ld` are the statistics' marginal correlations, sample size
# and LD source. NULL when fewer than two predictors explain at least `r2Cut`
# of the held-out variance; else the combined weights, every predictor's
# coefficient and held-out R^2, and the stack's own held-out R^2 and
# correlation, each split scored with the coefficients the others give.
# @noRd
#' @importFrom stats rnorm
.sumstatEnsemble <- function(own, external, cor, n, ld, r2Cut) {
    w <- .ldShrinkIntensity(ld, ld$nSamples)
    share <- .rssHeldOutTrainShare()
    nSplits <- .rssHeldOutSplits()
    noise <- sqrt(1 - w) *
        .ldSourceDraws(ld, nSplits) +
        sqrt(w) * matrix(rnorm(length(cor) * nSplits), ncol = nSplits)
    splits <- map(
        seq_len(nSplits),
        .sumstatEnsembleSplit,
        noise = noise * sqrt((1 - share) / (share * n)),
        cor = cor,
        ld = ld,
        own = own,
        external = external,
        w = w
    )
    full <- cbind(exec(cbind, !!!map(own, "weights")), external)
    r2 <- reduce(map(splits, .sumstatSplitExplained), `+`) / nSplits
    pass <- is.finite(r2) & r2 >= r2Cut
    if (sum(pass) < 2L) {
        return(NULL)
    }
    zeta <- .sumstatEnsembleZeta(splits, pass)
    list(
        weights = drop(full[, pass, drop = FALSE] %*% zeta),
        methodCoef = set_names(
            replace(rep(0, ncol(full)), pass, zeta),
            colnames(full)
        ),
        methodPerformance = set_names(r2, colnames(full)),
        metrics = .sumstatEnsembleCrossFit(splits, pass)
    )
}

# Split `k`: every method refit on its training part, beside the fixed
# predictors, and the moments of the held-out part over them.
# @noRd
.sumstatEnsembleSplit <- function(k, noise, cor, ld, own, external, w) {
    share <- .rssHeldOutTrainShare()
    trainCor <- .lassosumClampCor(cor + noise[, k])
    heldCor <- (cor - share * trainCor) / (1 - share)
    refits <- map(
        own,
        .sumstatEnsembleRefit,
        trainCor = trainCor,
        share = share
    )
    B <- cbind(exec(cbind, !!!refits), external)
    list(
        gram = (1 - w) * .ldCrossForms(ld, B) + w * crossprod(B),
        dvec = as.vector(crossprod(B, heldCor))
    )
}

# @noRd
.sumstatEnsembleRefit <- function(method, trainCor, share) {
    exec(
        method$refit$fn,
        trainCor = trainCor,
        share = share,
        !!!method$refit$args
    )
}

# Each predictor's share of a split's held-out variance it explains,
# 2 c_h' b - b' Sigma b.
# @noRd
.sumstatSplitExplained <- function(split) {
    2 * split$dvec - diag(split$gram)
}

# The stack's coefficients over the predictors `pass` from the splits'
# summed moments.
# @noRd
.sumstatEnsembleZeta <- function(splits, pass) {
    .solveEnsembleMoments(
        reduce(map(splits, .sumstatSplitGram, pass = pass), `+`),
        reduce(map(splits, .sumstatSplitDvec, pass = pass), `+`),
        sum(pass)
    )
}

# @noRd
.sumstatSplitGram <- function(split, pass) {
    split$gram[pass, pass, drop = FALSE]
}

# @noRd
.sumstatSplitDvec <- function(split, pass) {
    split$dvec[pass]
}

# The stack's own held-out R^2 and correlation, averaged over the splits,
# each scored with the coefficients the other splits give, as the
# individual-level stack is cross-fit over its folds.
# @noRd
.sumstatEnsembleCrossFit <- function(splits, pass) {
    scores <- map(
        seq_along(splits),
        .sumstatEnsembleFoldScore,
        splits = splits,
        pass = pass
    )
    c(
        corr = mean(map_dbl(scores, "corr"), na.rm = TRUE),
        rsq = mean(map_dbl(scores, "rsq"))
    )
}

# @noRd
.sumstatEnsembleFoldScore <- function(k, splits, pass) {
    zeta <- .sumstatEnsembleZeta(splits[-k], pass)
    variance <- drop(crossprod(
        zeta,
        .sumstatSplitGram(splits[[k]], pass) %*%
            zeta
    ))
    covariance <- sum(zeta * .sumstatSplitDvec(splits[[k]], pass))
    list(
        rsq = 2 * covariance - variance,
        corr = if (variance > 0) covariance / sqrt(variance) else NA_real_
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
# with any NA. Returns list(P, yObs, fold, samples): each row's fold and
# sample, both NULL unless every dataset gives them.
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
    stacked <- list(
        P = exec(rbind, !!!pMats),
        yObs = list_c(map(perDataset, "y")),
        fold = .ensembleAllOrNone(map(perDataset, "fold")),
        samples = .ensembleAllOrNone(map(perDataset, "samples"))
    )
    .ensembleDropIncomplete(stacked, nm$K)
}

# `parts` joined, or NULL when any is NULL.
# @noRd
.ensembleAllOrNone <- function(parts) {
    if (any(map_lgl(parts, is.null))) {
        return(NULL)
    }
    list_c(parts)
}

# Per-dataset prediction matrix + aligned outcome. Returns list(P, y,
# samples).
# @noRd
.ensembleDatasetMatrix <- function(predsD, yRaw, nm, contextIndex, d) {
    aln <- .ensembleAlignSamples(predsD, yRaw, nm$predNames, contextIndex, d)
    list(
        P = .ensembleBuildPd(predsD, nm, aln, contextIndex, d),
        y = aln$yD,
        samples = aln$samples
    )
}

# Each of `samples`' fold in dataset `d`'s partition, keyed by the dataset so
# the folds of different datasets stay apart; NULL unless every sample has
# one.
# @noRd
.ensembleDatasetFolds <- function(partition, samples, d) {
    if (is.null(partition) || is.null(samples)) {
        return(NULL)
    }
    fold <- partition$Fold[match(samples, partition$Sample)]
    if (anyNA(fold)) {
        return(NULL)
    }
    str_c(d, "\r", fold)
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
    .ensembleAlignPositional(yRaw, contextIndex, d, predSamples)
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
    list(
        yD = yD,
        predOrder = match(common, predSamples),
        nD = length(common),
        samples = common
    )
}

# Positional alignment fallback (no sample names on one side). The samples are
# the predictions' own names, when they have them.
# @noRd
.ensembleAlignPositional <- function(yRaw, contextIndex, d, predSamples) {
    yD <- if (is.matrix(yRaw) || is.data.frame(yRaw)) {
        .ensembleCheckContextIndex(contextIndex, ncol(yRaw), d)
        as.numeric(as.matrix(yRaw)[, contextIndex])
    } else {
        as.numeric(yRaw)
    }
    list(
        yD = yD,
        predOrder = seq_len(length(yD)),
        nD = length(yD),
        samples = predSamples[seq_len(length(yD))]
    )
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

# Drop the rows of `stacked` with any NA prediction/outcome; error when too
# few remain.
# @noRd
.ensembleDropIncomplete <- function(stacked, K) {
    complete <- complete.cases(stacked$P, stacked$yObs)
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
    list(
        P = stacked$P[complete, , drop = FALSE],
        yObs = stacked$yObs[complete],
        fold = stacked$fold[complete],
        samples = stacked$samples[complete]
    )
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

# Coerce a weight vector/matrix to a (variants x contexts) matrix, a vector's
# names becoming its rownames.
# @noRd
.ensembleAsMatrix <- function(w) {
    if (is.matrix(w)) {
        return(w)
    }
    matrix(w, ncol = 1, dimnames = list(names(w), NULL))
}

# One method's zeta-scaled contribution, in the first method's variant order
# (`ids`), or NULL when its weight matrix does not line up with the first
# one's shape.
# @noRd
.ensembleWeightTerm <- function(i, wtList, wtKeys, zeta, shape, ids) {
    wMat <- .ensembleAlignWeights(
        .ensembleAsMatrix(wtList[[wtKeys[i]]]),
        ids,
        wtKeys[i]
    )
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

# `wMat` in the variant order `ids`, when both name their variants: adding
# by position would add one variant's weight to another's. Weights on other
# variants are refused; unnamed ones can only be taken in the order given.
# @noRd
.ensembleAlignWeights <- function(wMat, ids, key) {
    own <- rownames(wMat)
    if (is.null(ids) || is.null(own)) {
        return(wMat)
    }
    if (anyDuplicated(own) > 0L || !setequal(own, ids)) {
        abort(glue(
            "ensembleWeights: the weights of '{key}' are on other variants ",
            "than the first method's; give every method the same variants."
        ))
    }
    wMat[ids, , drop = FALSE]
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
            shape = shape,
            ids = rownames(firstWt)
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
    if (.twasIsRegionList(f) && i <= length(f)) f[[i]] else NULL
}

# Whether `f` is a per-region list (region1, region2, ...): what a multi-range
# run with jointRegions = FALSE keeps in place of one fit or one
# cross-validation.
# @noRd
.twasIsRegionList <- function(f) {
    is.list(f) &&
        !is.null(names(f)) &&
        length(f) > 0L &&
        all(str_length(names(f)) > 0L) &&
        all(str_starts(names(f), "region"))
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

# One region block of .twasJointSpecGridWork().
# @noRd
.twasJointSpecRegionWork <- function(
    bi,
    spec,
    tokens,
    contexts,
    traits,
    marker,
    ctx
) {
    .jointSpecsWork(
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
        cache = ctx$twasWeights,
        mashPriors = ctx$mashPriors,
        cisWindow = ctx$cisWindow,
        region = ctx$xRegions[[bi]],
        regionIndex = bi,
        nRegions = length(ctx$xRegions),
        residualizationParam = ctx$residualizationParam,
        genotypeFilter = ctx$genotypeFilter,
        phenotypeFilter = ctx$phenotypeFilter,
        completeSamples = ctx$completeSamples,
        verbose = ctx$verbose
    )
}

# One region's univariate joint-cell fits (region `bi` of `xRegions`).
# @noRd
.twasQdsUnivRegionWork <- function(
    bi,
    cells,
    scope,
    marker,
    data,
    xRegions,
    norm,
    fineMappingResult,
    twasWeights,
    mashPriors,
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
        mashPriors = mashPriors,
        cisWindow = cisWindow,
        phenotypeFilter = phenotypeFilter,
        residualizationParam = residualizationParam,
        genotypeFilter = genotypeFilter,
        nPCs = nPCs,
        verbose = verbose
    )
    .qtlWorkOf(
        map(
            cells,
            .jointCellWork,
            pipeline = marker,
            data = data,
            scope = scope,
            tokens = norm$tokens,
            args = args
        ),
        .twasCombineCellRows
    )
}

# Combine the per-cell results of one region. rbind, NOT
# .twasMergeResultsByKey(): that one merges the SAME rows across regions,
# keeping results[[1]]'s row set. The top-PC rows are different rows, so
# merging would silently drop them.
# @noRd
.twasCombineCellRows <- function(perCell) {
    out <- compact(perCell)
    if (length(out) == 0L) {
        return(NULL)
    }
    if (length(out) == 1L) {
        return(out[[1L]])
    }
    reduce(out, .rbindTwasWeights)
}

# A run whose joint fits are the whole answer.
# @noRd
.twasQdsJointOnly <- function(
    work,
    data,
    genotypeFilter,
    residualizationParam
) {
    jointResult <- .qtlRunWork(
        work,
        data,
        genotypeFilter,
        residualizationParam
    )
    if (is.null(jointResult)) {
        msg <- glue(
            "twasWeightsPipeline(QtlDataset): no joint fits produced. ",
            "Check that the jointSpecification scope intersects the ",
            "available studies / contexts / traits."
        )
        abort(msg)
    }
    jointResult
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

# The stacked prediction/observed matrix for ensemble dataset `d`.
# @noRd
.ensembleDatasetRow <- function(d, cvResults, Y, nm, contextIndex) {
    m <- .ensembleDatasetMatrix(
        cvResults[[d]]$prediction,
        Y[[d]],
        nm,
        contextIndex,
        d
    )
    list_assign(
        m,
        fold = .ensembleDatasetFolds(
            cvResults[[d]]$samplePartition,
            m$samples,
            d
        )
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
        alpha = "numeric",
        acrossStudies = "logical",
        transferScale = "character"
    )
)

#' @title SR-TWAS Ensemble Settings
#' @description Whether and how \code{\link{twasWeightsPipeline}} stacks its
#'   per-method weights into an SR-TWAS ensemble.
#' @section Out-of-fold predictions:
#'   Stacking combines each method's \strong{out-of-fold} predictions. On
#'   individual-level data they come from cross-validation, so
#'   \code{enabled = TRUE} together with \code{CrossValidationParam(folds <
#'   2)} is an error. It used to be neither: the ensemble row was simply
#'   absent from the result, with \code{ensemble = TRUE} still reading as on.
#'
#'   \code{enabled} defaults to \code{FALSE} because
#'   \code{\link{CrossValidationParam}} defaults to no folds. To get an
#'   ensemble on individual-level data, ask for both.
#'
#'   Summary statistics have no samples to hold out, so there each method is
#'   refit on training parts of the statistics drawn from the LD, and scored
#'   on the held-out remainder, as PUMAS (Zhao et al. 2021) splits them: four
#'   splits, each keeping three quarters of the samples. The stack is fit to
#'   the held-out parts, and its coefficients combine the methods' full-data
#'   weights. Methods whose weights come from a fine-mapping fit cannot be
#'   refit, so they stay out of it. The stack is solved by quadprog there.
#' @param enabled Logical. Compute SR-TWAS ensemble weights. Default
#'   \code{FALSE}.
#' @param r2Threshold Minimum cross-validated \eqn{R^2} for a method to enter
#'   the stack, the held-out \eqn{R^2} on summary statistics. Default
#'   \code{0.01}. Stacking needs at least two methods to clear it.
#' @param solver Stacking solver, \code{"quadprog"} (default) or
#'   \code{"glmnet"}, which needs individual-level data.
#' @param alpha Elastic-net mixing parameter, used only when
#'   \code{solver = "glmnet"}. Default \code{1}.
#' @param acrossStudies Logical. On several studies joined by a
#'   \code{"study"} \code{jointSpecification}, also stack the other
#'   studies' fits into each study's ensemble, as MAGEPRO (individual-level
#'   data) and LEOPARD (summary statistics) combine ancestries. Each study's
#'   stack then holds its own methods' out-of-fold predictions (the fits
#'   across studies under \code{coupling} among them) and, as fixed
#'   predictors, every other study's own fit of each method, applied to its
#'   samples or statistics. Another study's fit across studies saw this
#'   study's data, so it is left out. The rows carry \code{jointStudies}.
#'   Default \code{FALSE}.
#' @param transferScale The scale another study's weights are carried to a
#'   study on, with \code{acrossStudies}: \code{"perSd"} (the default), as
#'   fitted, per standard deviation of each study's genotypes;
#'   \code{"perAllele"}, per allele, each weight rescaled by the ratio of the
#'   two studies' genotype standard deviations, which needs every study's
#'   allele frequencies.
#' @return A \code{EnsembleParam} object, a \code{\link{MethodParam}}.
#' @seealso \code{\link{CrossValidationParam}}
#' @examples
#' EnsembleParam(enabled = TRUE, r2Threshold = 0.05)
#' EnsembleParam(enabled = TRUE, acrossStudies = TRUE)
#' @export
EnsembleParam <- function(
    enabled = FALSE,
    r2Threshold = 0.01,
    solver = c("quadprog", "glmnet"),
    alpha = 1,
    acrossStudies = FALSE,
    transferScale = c("perSd", "perAllele")
) {
    solver <- arg_match(solver)
    transferScale <- arg_match(transferScale)
    assertFlag(acrossStudies)
    new(
        "EnsembleParam",
        enabled = enabled,
        r2Threshold = r2Threshold,
        solver = solver,
        alpha = alpha,
        acrossStudies = acrossStudies,
        transferScale = transferScale
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

# An ensemble on summary statistics is stacked from held-out moments, which
# quadprog solves; glmnet's stack needs individual-level predictions.
# @noRd
.ensembleAssertSumstats <- function(ensembleParam, inputKind) {
    if (
        !isTRUE(ensembleParam$enabled) ||
            !identical(ensembleParam$solver, "glmnet")
    ) {
        return(invisible(NULL))
    }
    abort(c(
        glue(
            "twasWeightsPipeline({inputKind}): EnsembleParam(solver = ",
            "\"glmnet\") stacks individual-level predictions, and summary ",
            "statistics are stacked from held-out moments."
        ),
        i = "Use solver = \"quadprog\" on summary statistics."
    ))
}

# `acrossStudies` stacks the studies a "study" jointSpecification joins, so it
# needs that axis alone among the specs and two or more studies to join.
# @noRd
.ensembleAssertAcross <- function(
    ensembleParam,
    parsedJointSpec,
    nStudies,
    inputKind
) {
    if (!isTRUE(ensembleParam$acrossStudies)) {
        return(invisible(NULL))
    }
    studyAxis <- keep(
        .fmSplitStudyAxisSpecs(parsedJointSpec)$study,
        .jsIsStudyAxisOnly
    )
    if (nStudies >= 2L && length(studyAxis) > 0L) {
        return(invisible(NULL))
    }
    abort(c(
        glue(
            "twasWeightsPipeline({inputKind}): EnsembleParam(acrossStudies = ",
            "TRUE) stacks the studies a \"study\" jointSpecification ",
            "joins, two or more of them."
        ),
        i = str_c(
            "Add \"study\" to jointSpecification on several studies, or ",
            "leave acrossStudies = FALSE."
        )
    ))
}

# `transferScale = "perAllele"` carries summary statistics' weights between
# studies by their genotype variances, which an LD matrix records only with
# its variants' allele frequencies (.ldRefVariance()): refused before any fit
# when one records none.
# @noRd
.ensembleAssertAlleleFrequencies <- function(ensemble, data, inputKind) {
    if (
        !isTRUE(ensemble$acrossStudies) ||
            !identical(ensemble$transferScale, "perAllele")
    ) {
        return(invisible(NULL))
    }
    lacking <- keep(ldReference(data), .ldRefLacksAf)
    if (length(lacking) == 0L) {
        return(invisible(NULL))
    }
    abort(c(
        glue(
            "twasWeightsPipeline({inputKind}): EnsembleParam(transferScale = ",
            "\"perAllele\") needs each study's allele frequencies, and ",
            "{length(lacking)} of its LD references record none."
        ),
        i = "Give the LdMatrix's variants an `AF` column, or use \"perSd\"."
    ))
}

# @noRd
.ldRefLacksAf <- function(ref) {
    (.ldIsMatrixRef(ref) || .ldRefIsProjection(ref)) && is.null(.ldRefAf(ref))
}

# A collection's summary-statistic studies hold its ensemble to their rules:
# quadprog's stack, and allele frequencies for a per-allele one.
# @noRd
.ensembleAssertCollection <- function(ensembleParam, data) {
    if (is.null(sumStats(data))) {
        return(invisible(NULL))
    }
    .ensembleAssertSumstats(ensembleParam, "MultiStudyQtlDataset")
    .ensembleAssertAlleleFrequencies(
        .ensembleResolve(ensembleParam),
        sumStats(data),
        "MultiStudyQtlDataset"
    )
}

# How many studies a collection holds: its individual-level ones and those
# its summary statistics name.
# @noRd
.msStudyCount <- function(data) {
    ss <- sumStats(data)
    length(data) + if (is.null(ss)) 0L else n_distinct(studyName(ss))
}

# The ensemble settings a collection's component runs on its own: the same,
# less `acrossStudies`, which the collection's study axis carries out.
# @noRd
.ensembleWithinStudy <- function(ensembleParam) {
    if (is.null(ensembleParam)) {
        return(NULL)
    }
    EnsembleParam(
        enabled = isTRUE(ensembleParam$enabled),
        r2Threshold = ensembleParam$r2Threshold,
        solver = ensembleParam$solver,
        alpha = ensembleParam$alpha,
        transferScale = ensembleParam$transferScale
    )
}

# The ensemble settings as a plain list, so a cfg record can carry the group
# as one named field instead of four loose ones.
# @noRd
.ensembleResolve <- function(ensembleParam) {
    list(
        enabled = isTRUE(ensembleParam$enabled),
        r2Threshold = ensembleParam$r2Threshold %||% 0.01,
        solver = ensembleParam$solver %||% "quadprog",
        alpha = ensembleParam$alpha %||% 1,
        acrossStudies = isTRUE(ensembleParam$acrossStudies),
        transferScale = ensembleParam$transferScale %||% "perSd"
    )
}
