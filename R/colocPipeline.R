# fastenloc refuses to do both: given explicit p1/p2/p12 it prints
# "Applying user-specified colocalization priors, skipping enrichment
# analysis" and skips the enrichment step. pecotmr errors instead of choosing
# silently -- the silence is what made the two parameterisations easy to
# confuse in the first place.
# @noRd
#' @include MethodParam.R
#' @include ld.R
#' @include fineMappingPipeline.R
NULL

.colocRefusePriorsWithEnrichment <- function(priors, enrichment, wasDefault) {
    if (is.null(enrichment) || wasDefault) {
        return(invisible(NULL))
    }
    abort(glue(
        "colocPipeline: `priors` and `enrichment` cannot both be given. ",
        "With an enrichment table the priors are enloc's, derived from the ",
        "genome-wide QTL results by qtlEnrichmentPipeline(); supplying them ",
        "here as well would mean scoring with one set and enriching with ",
        "another. Drop `priors`, or drop `enrichment` to use your own."
    ))
}

# The three argument groups, checked in one place so the pipeline body does
# not open with three near-identical lines.
# @noRd
.colocAssertGroups <- function(priors, lbfFilterArgs, methodArgs) {
    .assertMethodParam(priors, "ColocPriorParam", "priors")
    .assertMethodParam(lbfFilterArgs, "ColocLbfFilterParam", "lbfFilter")
    .assertMethodOptions(methodArgs, "ColocOptions", "methodArgs")
}

# The coloc.bf_bf option list: the priors pecotmr owns, injected under
# coloc's names, plus whatever else the caller asked for. The priors are also
# read by the enrichment adjustment, so they must agree -- setting one in
# `methodArgs` too is refused rather than silently resolved.
# @noRd
.colocEngineArgs <- function(methodArgs, priors) {
    user <- as.list(methodArgs)
    clash <- intersect(names(user), c("p1", "p2", "p12"))
    if (length(clash) > 0L) {
        abort(glue(
            "colocPipeline: {str_flatten(clash, ', ')} ",
            "{if (length(clash) == 1L) 'is' else 'are'} set through ",
            "`priors`, which the enrichment adjustment reads as well, so ",
            "{if (length(clash) == 1L) 'it' else 'they'} cannot also be ",
            "given in `methodArgs`."
        ))
    }
    user
}

#' @rdname ColocPriorParam
#' @aliases ColocPriorParam-class
#' @exportClass ColocPriorParam
setClass(
    "ColocPriorParam",
    contains = "MethodParam",
    slots = c(
        p1 = "numeric",
        p2 = "numeric",
        p12 = "numeric"
    )
)

#' @title Prior Probabilities For Colocalisation
#' @description The per-variant prior probabilities coloc scores with, all
#'   three forwarded to \code{coloc::coloc.bf_bf}.
#'
#'   These describe the baseline analysis only. When \code{enrichment} is
#'   supplied to \code{\link{colocPipeline}} the priors are \emph{derived}
#'   from the enrichment estimate instead, the way enloc derives them, and
#'   supplying both is an error rather than a silent precedence rule.
#' @param p1 Prior probability a variant is causal for the QTL trait.
#'   Default \code{1e-4}.
#' @param p2 Prior probability a variant is causal for the GWAS trait.
#'   Default \code{1e-4}.
#' @param p12 Prior probability a variant is causal for both. Default
#'   \code{5e-6}.
#' @return A \code{ColocPriorParam} object, a \code{\link{MethodParam}}.
#' @examples
#' ColocPriorParam(p12 = 1e-5)
#' @export
ColocPriorParam <- function(
    p1 = 1e-4,
    p2 = 1e-4,
    p12 = 5e-6
) {
    new(
        "ColocPriorParam",
        p1 = p1,
        p2 = p2,
        p12 = p12
    )
}

#' @rdname ColocLbfFilterParam
#' @aliases ColocLbfFilterParam-class
#' @exportClass ColocLbfFilterParam
setClass(
    "ColocLbfFilterParam",
    contains = "MethodParam",
    slots = c(
        filterLbfCs = "logical",
        secondary = "numeric_OR_NULL",
        concentration = "numeric",
        priorTol = "numeric"
    )
)

#' @title Credible-Set Filtering For Colocalisation
#' @description Whether and how to restrict each fine-mapping result to its
#'   credible sets before scoring, using the log Bayes factors.
#' @param filterLbfCs Logical. Restrict to credible-set variants. Default
#'   \code{FALSE}.
#' @param secondary Optional secondary coverage levels used when selecting
#'   credible sets. \code{NULL} (default) uses the primary sets only.
#' @param concentration Minimum share of the credible set's posterior mass a
#'   variant must carry to be kept. Default \code{0.5}.
#' @param priorTol Prior-variance cutoff for the default filter: effects with
#'   \code{V <= priorTol} are dropped. Default \code{1e-9}. This is the
#'   default filter's own parameter, which is why it lives here rather than
#'   beside it --- \code{filterLbfCs} and \code{secondary} select a
#'   different filter, and then it does not apply.
#' @return A \code{ColocLbfFilterParam} object, a \code{\link{MethodParam}}.
#' @examples
#' ColocLbfFilterParam(filterLbfCs = TRUE)
#' @export
ColocLbfFilterParam <- function(
    filterLbfCs = FALSE,
    secondary = NULL,
    concentration = 0.5,
    priorTol = 1e-9
) {
    new(
        "ColocLbfFilterParam",
        filterLbfCs = filterLbfCs,
        secondary = secondary,
        concentration = concentration,
        priorTol = priorTol
    )
}

# `methods` takes either an engine name or a FineMappingMethodsParam()
# record carrying that engine's own options, exactly as fineMappingPipeline()
# does -- this bundle is forwarded there whole. Now exact: it was
# `MethodOptions` only while the aggregator was one, so any engine's
# argument bag was accepted here.
setClassUnion(
    "MethodSelection",
    c("character", "FineMappingMethodsParam")
)

#' @rdname GwasFineMappingParam
#' @aliases GwasFineMappingParam-class
#' @exportClass GwasFineMappingParam
setClass(
    "GwasFineMappingParam",
    contains = "MethodParam",
    slots = c(
        methods = "MethodSelection",
        credibleSetParam = "CredibleSetParam",
        susieRssParam = "SusieRssParam",
        panelFilterParam = "PanelFilterParam",
        initializeWithSusieInf = "logical",
        fitRetention = "character"
    )
)

#' @title Inline GWAS Fine-Mapping Settings
#' @description How \code{\link{colocPipeline}} fine-maps \code{gwasInput}
#'   when it is summary statistics rather than an existing fine-mapping
#'   result. Every field is forwarded to
#'   \code{\link{fineMappingPipeline}}; the settings that pipeline exposes
#'   but that cannot apply here are deliberately absent:
#'   \code{fineMappingResult} (the pipeline is being asked to produce one),
#'   \code{crossValidation} (refused on summary statistics) and
#'   \code{residualization} (nothing to regress out).
#'
#'   The whole bundle is inert when \code{gwasInput} is already a
#'   fine-mapping result, since no fit is run.
#' @param methods Fine-mapping methods, as a character vector or a
#'   \code{\link{FineMappingMethodsParam}} record. Default \code{"susie"}.
#' @param credibleSetParam How credible sets are built, built with
#'   \code{\link{CredibleSetParam}}. \code{coverage} and \code{L} matter
#'   most here: they decide the sets whose log Bayes factors coloc scores, so
#'   leaving them at the defaults while the QTL side used something else
#'   compares two differently-built sets.
#' @param susieRssParam The summary-statistics solver, built with
#'   \code{\link{SusieRssParam}} --- \code{serFallback}, \code{rMismatch},
#'   \code{rFinite} and the \code{susie_rss} control list. What a GWAS block
#'   with an imperfect LD panel needs.
#' @param panelFilterParam LD-reference-panel filters, built with
#'   \code{\link{PanelFilterParam}}.
#' @param initializeWithSusieInf Logical. Chain a SuSiE-inf fit, when
#'   \code{methods} asks for \code{susieInf} alongside \code{susie}.
#'   Default \code{TRUE}.
#' @param fitRetention How much of each fit is kept: \code{"slim"} (default)
#'   or \code{"full"}. Only observable when
#'   \code{returnGwasFineMapping = TRUE}, which is when the fine-mapping
#'   result is handed back for other uses.
#' @return \code{GwasFineMappingParam} returns a
#'   \code{GwasFineMappingParam} object, a \code{\link{MethodParam}}.
#'   Each accessor returns that setting's value --- for the three nested
#'   bundles, the \code{MethodParam} itself --- and each replacement form
#'   returns a modified copy.
#' @examples
#' GwasFineMappingParam(
#'     credibleSetParam = CredibleSetParam(coverage = 0.9)
#' )
#' @export
GwasFineMappingParam <- function(
    methods = "susie",
    credibleSetParam = CredibleSetParam(),
    susieRssParam = SusieRssParam(),
    panelFilterParam = PanelFilterParam(),
    initializeWithSusieInf = TRUE,
    fitRetention = c("slim", "full")
) {
    fitRetention <- arg_match(fitRetention)
    new(
        "GwasFineMappingParam",
        methods = methods,
        credibleSetParam = credibleSetParam,
        susieRssParam = susieRssParam,
        panelFilterParam = panelFilterParam,
        initializeWithSusieInf = initializeWithSusieInf,
        fitRetention = fitRetention
    )
}

#' @title Colocalization Pipeline (coloc.bf_bf over paired LBF matrices)
#' @description Per-region pipeline that pairs two fine-mapping result
#'   collections and runs \code{coloc::coloc.bf_bf} per (first-side tuple,
#'   second-side tuple) pair to produce per-pair colocalization posterior
#'   probabilities PP.H0-PP.H4.
#'
#'   Either side may be a \code{\link{QtlFineMappingResult}} or a
#'   \code{\link{GwasFineMappingResult}}: QTL-GWAS, QTL-QTL (two molecular
#'   phenotypes) and GWAS-GWAS (two diseases) all run through one code path,
#'   since nothing below the identity tuple depends on which flavour a side
#'   is. The second side may also be handed in as summary statistics
#'   (\code{\link{QtlSumStats}} or \code{\link{GwasSumStats}}), which are
#'   fine-mapped inline.
#'
#'   The argument names keep the QTL / GWAS wording of the common case. What
#'   they mean generally is: \code{qtlFineMappingResult} is the side whose
#'   identity is reported in the unprefixed \code{study} / \code{context} /
#'   \code{trait} / \code{method} columns, and \code{gwasInput} is the side
#'   reported in the \code{gwas}-prefixed ones. A GWAS side has no context or
#'   trait axis, so those two columns are \code{NA} for it.
#'
#' @section Why \code{coloc.bf_bf} and not \code{coloc.susie}:
#' The prior \code{colocWrapper} (now stubbed) used
#' \code{coloc::coloc.bf_bf} on the SuSiE \code{lbf_variable} matrices
#' directly. That choice carries three behaviours that
#' \code{coloc::coloc.susie} does not expose:
#' \itemize{
#'   \item \strong{fSuSiE support}: the LBF matrix lives at a different
#'     slot for fSuSiE fits (\code{fsusie_result$lBF}) and gets stacked
#'     into a single matrix.
#'   \item \strong{Effect filtering}: \code{lbfFilter}'s
#'     \code{filterLbfCs} keeps only effects that produced a credible
#'     set, its \code{secondary} keeps effects at a secondary coverage;
#'     otherwise the default filter drops effects whose prior variance is
#'     below \code{priorTol}.
#'   \item \strong{Per-tuple LBF reuse}: each second-side tuple's LBF
#'     matrix is extracted once and scored against every first-side
#'     tuple, so the filtering above is applied once per tuple rather
#'     than once per pair.
#' }
#' This pipeline preserves all three.
#'
#'   Second-side input dispatch:
#'   \itemize{
#'     \item \code{gwasInput} is a \code{\link{QtlSumStats}} or a
#'           \code{\link{GwasSumStats}}: it is fine-mapped inline by
#'           \code{\link{fineMappingPipeline}} with the supplied
#'           \code{GwasFineMappingParam()} (methods default \code{"susie"}).
#'     \item \code{gwasInput} is a fine-mapping result: used directly; no
#'           inline fine-mapping.
#'   }
#'
#' @section LD-sketch compatibility check: If
#'   \code{ldSketch(qtlFineMappingResult)} is non-\code{NULL}, it must come
#'   from the same reference panel as the LD sketch on \code{gwasInput}: the
#'   same samples, and the same allele orientation on the variants the two
#'   carry in common. The two need NOT carry the same variants --- running
#'   \code{\link{summaryStatsQc}} on the two sides separately normally leaves
#'   each sketch trimmed to its own surviving variants, and LD is looked up per
#'   variant, so a partial overlap only warns (once per session). No shared
#'   variant at all, a different sample set, or a swapped A1/A2 on a shared
#'   variant is a hard error. When the first side's \code{ldSketch} is
#'   \code{NULL} (individual-level fit), the validation is skipped on that side
#'   and the second side's panel is what the result carries forward.
#'
#' @param qtlFineMappingResult The first side: a
#'   \code{\link{QtlFineMappingResult}} or a
#'   \code{\link{GwasFineMappingResult}} (required).
#' @param gwasInput The second side: a \code{\link{QtlFineMappingResult}}, a
#'   \code{\link{GwasFineMappingResult}}, or the summary statistics to
#'   fine-map inline (\code{\link{QtlSumStats}} or
#'   \code{\link{GwasSumStats}}).
#' @param lbfFilterArgs Credible-set filtering applied before scoring,
#'   built with
#'   \code{\link{ColocLbfFilterParam}}. \code{filterLbfCs} keeps only effects
#'   that produced a credible set; supplying \code{secondary} instead runs a
#'   concentration filter at that coverage, where an effect is kept only if
#'   its credible set spans fewer than
#'   \code{nVariants * secondary * concentration} variants. Diffuse effects
#'   are dropped before the LBF matrix reaches
#'   \code{coloc::coloc.bf_bf}. \code{priorTol} is the default filter's own
#'   cutoff (effects with \code{V <= priorTol} are dropped) and does not
#'   apply once \code{filterLbfCs} or \code{secondary} selects another.
#' @param priors Per-variant prior probabilities, built with
#'   \code{\link{ColocPriorParam}}: \code{p1} (QTL signal), \code{p2} (GWAS
#'   signal) and \code{p12} (shared). Mutually exclusive with
#'   \code{enrichment}, which derives the same three.
#' @param gwasFineMappingArgs How \code{gwasInput} is fine-mapped when it is
#'   summary statistics rather than a fine-mapping result, built with
#'   \code{\link{GwasFineMappingParam}}: \code{methods},
#'   \code{credibleSetParam}, \code{susieRssParam},
#'   \code{panelFilterParam}, \code{initializeWithSusieInf} and
#'   \code{fitRetention}. Inert when \code{gwasInput} is already
#'   fine-mapped. \code{credibleSetParam} is the one to check: its
#'   \code{coverage} / \code{maxNumSingleEffects} decide the sets whose
#'   LBFs are
#'   scored, so defaults here against a differently-built QTL side compare
#'   two different things.
#' @param returnGwasFineMapping Logical. When \code{TRUE}, attach the
#'   fine-mapping result computed from \code{gwasInput} on the returned object
#'   as attribute \code{"gwasFineMapping"}. Default \code{FALSE}.
#' @param enrichment Optional data.frame of per-pair enrichment factors with
#'   columns \code{gwasStudy}, \code{qtlStudy}, \code{qtlContext},
#'   \code{enrichment}, and optionally \code{gwasContext} / \code{gwasTrait}
#'   (which the join uses when present, and which
#'   \code{\link{qtlEnrichmentPipeline}} emits for a QTL outcome side). Output
#'   of \code{\link{qtlEnrichmentPipeline}}. When non-\code{NULL} the pair is
#'   scored by enloc rather than coloc: the priors come from that table's
#'   \code{colocP1} / \code{colocP2} / \code{colocP12} columns, and the
#'   report gains \code{RCP} and \code{LCP}. Pairs without a matching
#'   enrichment row fall back to the baseline priors with a warning. Default
#'   \code{NULL} (baseline coloc).
#' @param adjustPips Logical, default \code{TRUE}. When TRUE, before any
#'   per-pair inference the QTL and GWAS fine-mapping result collections are
#'   reconciled with \code{\link{intersectVariants}}, so both sides are scored
#'   on the variants they share and their PIPs are renormalized to that shared
#'   set. This matters in two scenarios: (1) the user declined to impute
#'   missing variants in the GWAS \code{SumStats} and the QTL fine-mapping
#'   input has additional variants; (2) the GWAS fine-mapping result contains
#'   variants not present in the QTL fine-mapping result. Pass \code{FALSE} to
#'   use the FMRs as supplied.
#'
#'   How much of each effect's posterior survived that restriction is reported
#'   per result row as \code{qtlRetainedMass} / \code{gwasRetainedMass}; a low
#'   value means the effect was largely built on variants the other side does
#'   not carry, so its colocalization evidence rests on little retained signal.
#'   Both are \code{NA} when \code{adjustPips = FALSE}, since nothing was
#'   measured.
#' @param alleleFlip Logical, default \code{TRUE}. When TRUE, align LBF columns
#'   between the QTL and GWAS by (chrom, pos) with ref/alt swaps recognized (LBF
#'   is coding-invariant, so no sign change is needed); when FALSE, match on
#'   exact alleles only, so a ref/alt swap is treated as a distinct variant.
#' @param methodArgs Additional arguments forwarded to
#'   \code{coloc::coloc.bf_bf}, built with \code{\link{ColocOptions}}. The
#'   prior probabilities are set through \code{priors} instead, since the
#'   enrichment adjustment reads them too.
#' @return A \code{\linkS4class{ColocResult}}: one element per tested
#'   (first-side credible set, second-side credible set, block) pair, holding
#'   that pair's aligned variants with their \code{SNP.PP.H4}. Pair-level
#'   metadata carries the identity columns (\code{study}, \code{context},
#'   \code{trait}, \code{method}, \code{gwasStudy}, \code{gwasContext},
#'   \code{gwasTrait}, \code{gwasMethod}), the block and stable
#'   credible-set ids (\code{blockId}, \code{qtlCs}, \code{gwasCs}), the
#'   standard coloc fields (\code{idx1}, \code{idx2}, \code{nSnps},
#'   \code{hit1}, \code{hit2}, \code{PP.H0.abf} \ldots \code{PP.H4.abf}) and
#'   the reconciliation diagnostics \code{qtlRetainedMass} /
#'   \code{gwasRetainedMass}. When \code{enrichment} is supplied the pair is
#'   scored by enloc instead, and four further columns appear:
#'   \code{enrichment} and \code{p12Used} report the per-pair factor and the
#'   shared-signal prior enloc derived, while \code{RCP} and \code{LCP} are
#'   enloc's regional and locus colocalisation probabilities. Those two are
#'   the same numbers as \code{PP.H4.abf} and
#'   \code{PP.H3.abf + PP.H4.abf}, reported under enloc's names because that
#'   is what readers of enloc output expect to find.
#'
#'   Project it with \code{\link{colocPairs}} (the flat table this
#'   pipeline used to return, also available as \code{as.data.frame}),
#'   \code{\link{colocVariants}}, \code{\link{colocCredibleSets}} or
#'   \code{\link{colocGenes}}.
#' @examples
#' data(qtlFineMappingLbfExample)
#' data(gwasFineMappingLbfExample)
#' colocPipeline(qtlFineMappingLbfExample,
#'   gwasInput = gwasFineMappingLbfExample)
#' # Either side may be a QTL result. Pairing this collection against itself
#' # colocalizes its two contexts, and reports the second side's context and
#' # trait in gwasContext / gwasTrait.
#' res <- colocPipeline(qtlFineMappingLbfExample,
#'   gwasInput = qtlFineMappingLbfExample)
#' unique(colocPairs(res)[, c("context", "gwasContext")])
#' @export
colocPipeline <- function(
    qtlFineMappingResult,
    gwasInput,
    priors = ColocPriorParam(),
    lbfFilterArgs = ColocLbfFilterParam(),
    gwasFineMappingArgs = GwasFineMappingParam(),
    returnGwasFineMapping = FALSE,
    enrichment = NULL,
    adjustPips = TRUE,
    alleleFlip = TRUE,
    methodArgs = ColocOptions()
) {
    prepared <- .colocPrepare(
        qtlFineMappingResult = qtlFineMappingResult,
        gwasInput = gwasInput,
        priors = priors,
        gwasFineMappingArgs = gwasFineMappingArgs,
        enrichment = enrichment,
        adjustPips = adjustPips,
        methodArgs = methodArgs,
        lbfFilterArgs = lbfFilterArgs,
        priorsMissing = missing(priors)
    )
    # Pre-extract per-GWAS-tuple LBF matrices: group the GWAS FMR by study,
    # stack each study's LBF rows, and store per-(study, method) batched
    # matrices (reproduces the legacy row-wise combine per xQTL).
    gwasLbfByPair <- .colocPreextractGwasLbf(prepared$gwasFmr, lbfFilterArgs)
    # Both exits -- no scorable GWAS effect, and the scored result -- shape
    # their output from the same five values, so they are named once.
    shape <- list(
        useEnrichment = prepared$useEnrichment,
        qtlFineMappingResult = prepared$qtlFineMappingResult,
        gwasFmr = prepared$gwasFmr,
        gwasInput = gwasInput,
        returnGwasFineMapping = returnGwasFineMapping
    )
    if (length(gwasLbfByPair) == 0L) {
        return(exec(.colocEarlyReturn, !!!shape))
    }
    results <- .colocScoreAll(
        prepared,
        gwasLbfByPair,
        lbfFilterArgs = lbfFilterArgs,
        enrichment = enrichment,
        priors = priors,
        alleleFlip = alleleFlip
    )
    exec(.colocFinalize, results, !!!shape)
}

# Score every QTL tuple against the pre-extracted GWAS LBF matrices.
# @noRd
.colocScoreAll <- function(
    prepared,
    gwasLbfByPair,
    lbfFilterArgs,
    enrichment,
    priors,
    alleleFlip
) {
    list_flatten(map(
        seq_len(nrow(prepared$qtlFineMappingResult)),
        .colocScoreQtlTuple,
        qtlFineMappingResult = prepared$qtlFineMappingResult,
        gwasLbfByPair = gwasLbfByPair,
        lbfFilterArgs = lbfFilterArgs,
        useEnrichment = prepared$useEnrichment,
        enrichment = enrichment,
        priors = priors,
        alleleFlip = alleleFlip,
        ColocOptions = prepared$ColocOptions
    ))
}

# Validate the inputs, resolve the GWAS side to a fine-mapping result, check
# the two sides share an LD sketch, and apply the PIP adjustment. Everything
# that must happen before any tuple is scored, in one step because each part
# depends on the one before it.
# @noRd
.colocPrepare <- function(
    qtlFineMappingResult,
    gwasInput,
    priors,
    gwasFineMappingArgs,
    enrichment,
    adjustPips,
    methodArgs,
    lbfFilterArgs,
    priorsMissing
) {
    .colocAssertGroups(priors, lbfFilterArgs, methodArgs)
    .colocRefusePriorsWithEnrichment(priors, enrichment, priorsMissing)
    useEnrichment <- !is.null(enrichment)
    .colocValidateInputs(
        gwasInput = gwasInput,
        qtlFineMappingResult = qtlFineMappingResult,
        enrichment = enrichment,
        useEnrichment = useEnrichment
    )
    rawGwasFmr <- .colocResolveGwasFmr(gwasInput, gwasFineMappingArgs)
    .colocRequireMatchingLdSketches(
        ldSketch(qtlFineMappingResult),
        ldSketch(rawGwasFmr)
    )
    adjusted <- .colocMaybeAdjustPips(
        adjustPips = adjustPips,
        qtlFineMappingResult = qtlFineMappingResult,
        gwasFmr = rawGwasFmr
    )
    list(
        ColocOptions = .colocEngineArgs(methodArgs, priors),
        useEnrichment = useEnrichment,
        qtlFineMappingResult = adjusted$qtlFineMappingResult,
        gwasFmr = adjusted$gwasFmr
    )
}

# Validate the enrichment table (when supplied), the coloc package, and the
# input object classes.
# @noRd
.colocValidateInputs <- function(
    gwasInput,
    qtlFineMappingResult,
    enrichment,
    useEnrichment
) {
    if (useEnrichment) {
        .colocValidateEnrichment(enrichment)
    }
    if (!requireNamespace("coloc", quietly = TRUE)) {
        abort("Package 'coloc' is required for colocPipeline.")
    }
    if (!methods::is(qtlFineMappingResult, "FineMappingResultBase")) {
        msg <- glue(
            "`qtlFineMappingResult` must be a QtlFineMappingResult or a ",
            "GwasFineMappingResult ",
            "(got class '{class(qtlFineMappingResult)[[1L]]}')."
        )
        abort(msg)
    }
    if (
        !methods::is(gwasInput, "SumStatsBase") &&
            !methods::is(gwasInput, "FineMappingResultBase")
    ) {
        msg <- glue(
            "`gwasInput` must be a fine-mapping result ",
            "(QtlFineMappingResult / GwasFineMappingResult) or summary ",
            "statistics (QtlSumStats / GwasSumStats) ",
            "(got class '{class(gwasInput)[[1L]]}')."
        )
        abort(msg)
    }
    invisible(NULL)
}

# The enrichment table must be a data.frame carrying the required id + value
# columns.
# @noRd
#' @importFrom checkmate assertDataFrame assertNames
.colocValidateEnrichment <- function(enrichment) {
    # The producer is named in .var.name so a caller passing the wrong object
    # is pointed at what produces the right one.
    label <- "enrichment (output of qtlEnrichmentPipeline)"
    assertDataFrame(enrichment, .var.name = label)
    assertNames(
        colnames(enrichment),
        must.include = c("gwasStudy", "qtlStudy", "qtlContext", "enrichment"),
        what = "colnames",
        .var.name = label
    )
    .colocValidateEnrichmentKeys(enrichment)
    invisible(NULL)
}

# One factor per pair: rows that repeat the identity the lookup joins on would
# make the applied enrichment depend on row order, so they are refused here
# rather than resolved by taking the first.
# @noRd
.colocValidateEnrichmentKeys <- function(enrichment) {
    idCols <- intersect(
        c("gwasStudy", "gwasContext", "gwasTrait", "qtlStudy", "qtlContext"),
        colnames(enrichment)
    )
    ids <- select(as_tibble(enrichment), all_of(idCols))
    if (nrow(distinct(ids)) == nrow(ids)) {
        return(invisible(NULL))
    }
    msg <- glue(
        "`enrichment` has repeated ({str_flatten(idCols, ', ')}) rows; ",
        "each pair needs exactly one enrichment factor."
    )
    abort(msg)
}

# Resolve the second side to a fine-mapping collection (fine-map QC'd sumstats
# when summary statistics are passed, whichever flavour they are).
# @noRd
.colocResolveGwasFmr <- function(gwasInput, gwasFineMappingArgs) {
    if (methods::is(gwasInput, "FineMappingResultBase")) {
        # Already fine-mapped: nothing in the bundle applies, so a setting
        # given here would be dropped rather than honoured.
        .colocAssertGwasFmrUnset(gwasFineMappingArgs, class(gwasInput)[[1L]])
        return(gwasInput)
    }
    .assertMethodParam(
        gwasFineMappingArgs,
        "GwasFineMappingParam",
        "gwasFineMapping"
    )
    if (length(qcInfo(gwasInput)) == 0L) {
        msg <- glue(
            "colocPipeline: gwasInput ({class(gwasInput)[[1L]]}) has no QC ",
            "record. Call summaryStatsQc() first."
        )
        abort(msg)
    }
    # The bundle's terminal: fineMappingPipeline takes these as its own
    # arguments, so it is unrolled once, here.
    fineMappingPipeline(
        gwasInput,
        methods = gwasFineMappingArgs$methods %||% "susie",
        credibleSetParam = gwasFineMappingArgs$credibleSetParam %||%
            CredibleSetParam(),
        susieRssParam = gwasFineMappingArgs$susieRssParam %||% SusieRssParam(),
        panelFilterParam = gwasFineMappingArgs$panelFilterParam %||%
            PanelFilterParam(),
        initializeWithSusieInf = gwasFineMappingArgs$initializeWithSusieInf %||%
            TRUE,
        fitRetention = gwasFineMappingArgs$fitRetention %||% "slim"
    )
}

# `gwasFineMapping` configures the fine-mapping run colocPipeline does on a
# raw GWAS input. A caller who has already fine-mapped the GWAS has no such
# run, so a non-default bundle is an instruction that cannot be carried out.
# @noRd
.colocAssertGwasFmrUnset <- function(gwasFineMappingArgs, inputClass) {
    if (!.isMethodParam(gwasFineMappingArgs)) {
        return(invisible(NULL))
    }
    defaults <- GwasFineMappingParam()
    set <- keep(
        names(defaults),
        .colocGwasFmrFieldSet,
        given = gwasFineMappingArgs,
        defaults = defaults
    )
    if (length(set) == 0L) {
        return(invisible(NULL))
    }
    abort(glue(
        "colocPipeline: `gwasFineMapping` ",
        "({str_flatten(str_c(set, ' ='), ', ')}) configures the GWAS ",
        "fine-mapping run, but `gwasInput` is already a ",
        "{inputClass} -- there is no run to configure. ",
        "Drop `gwasFineMapping`, or pass the raw GwasSumStats instead."
    ))
}

# @noRd
.colocGwasFmrFieldSet <- function(field, given, defaults) {
    !isTRUE(all.equal(given[[field]], defaults[[field]]))
}

# Optional PIP renormalization, via the symmetric reconciliation verb.
#
# Coloc is the symmetric case: both sides have to end up scored on the SAME
# variant set, or a pair's two halves are not comparable. That is what
# intersectVariants() gives, and it replaces the older hand-rolled pass which
# adjusted each side to the UNION of the other's variants -- a union is not an
# intersection, so the two sides could still end up on different sets.
# @noRd
.colocMaybeAdjustPips <- function(
    adjustPips,
    qtlFineMappingResult,
    gwasFmr
) {
    unchanged <- list(
        qtlFineMappingResult = qtlFineMappingResult,
        gwasFmr = gwasFmr
    )
    if (!isTRUE(adjustPips)) {
        return(unchanged)
    }
    if (nrow(qtlFineMappingResult) == 0L || nrow(gwasFmr) == 0L) {
        return(unchanged)
    }
    both <- intersectVariants(qtlFineMappingResult, gwasFmr)
    list(qtlFineMappingResult = both$x, gwasFmr = both$y)
}

# The LD reference the result carries forward, so colocCredibleSets() can
# recompute purity (section 3.7) without being handed a sketch separately. The
# two sides' sketches are already required to match by
# .colocRequireMatchingLdSketches, so either one identifies the panel -- but a
# first side fit on individual-level data carries none, and then the second
# side's panel is the only one there is.
# @noRd
.colocLdSketch <- function(qtlFineMappingResult, gwasFmr) {
    ldSketch(qtlFineMappingResult) %||% ldSketch(gwasFmr)
}

# Empty-result early return (attaching the GWAS fine-mapping when requested).
# @noRd
.colocEarlyReturn <- function(
    useEnrichment,
    qtlFineMappingResult,
    gwasFmr,
    gwasInput,
    returnGwasFineMapping
) {
    out <- .colocEmptyResult(
        enriched = useEnrichment,
        ldSketch = .colocLdSketch(qtlFineMappingResult, gwasFmr)
    )
    .colocAttachGwasFm(out, returnGwasFineMapping, gwasInput, gwasFmr)
}

# Score one QTL tuple against every pre-extracted GWAS pair -> summary rows.
# @noRd
.colocScoreQtlTuple <- function(
    qi,
    qtlFineMappingResult,
    gwasLbfByPair,
    lbfFilterArgs,
    useEnrichment,
    enrichment,
    priors,
    alleleFlip,
    ColocOptions
) {
    q <- .colocQtlTupleInfo(qi, qtlFineMappingResult)
    qLbfInfo <- .colocExtractLbfFromEntry(
        q$parts,
        lbfFilterArgs,
        label = q$label
    )
    if (is.null(qLbfInfo)) {
        return(list())
    }
    scored <- list_assign(
        q,
        retainedMass = qLbfInfo$retainedMass,
        effect = qLbfInfo$effect
    )
    compact(map(
        gwasLbfByPair,
        .colocScorePairAt,
        qLbfInfo = qLbfInfo,
        q = scored,
        useEnrichment = useEnrichment,
        enrichment = enrichment,
        priors = priors,
        alleleFlip = alleleFlip,
        ColocOptions = ColocOptions
    ))
}

# Identity + row payload + log label for one first-side tuple.
# @noRd
.colocQtlTupleInfo <- function(qi, qtlFineMappingResult) {
    fmr <- qtlFineMappingResult
    ident <- .colocTupleIdentity(fmr, qi)
    c(
        ident,
        list(
            parts = .fmrRowParts(fmr, qi),
            label = .fmrTupleLabel(.fmrSideName(fmr), ident)
        )
    )
}

# The identity tuple of one fine-mapping row, whichever flavour the collection
# is (see .fmrIdentityAt for the absent-axis rule).
# @noRd
.colocTupleIdentity <- function(fmr, ri) {
    list(
        study = .fmrIdentityAt(fmr, "study", ri),
        context = .fmrIdentityAt(fmr, "context", ri),
        trait = .fmrIdentityAt(fmr, "trait", ri),
        method = .fmrIdentityAt(fmr, "method", ri)
    )
}

# With an enrichment table this is an enloc run: the prior depends on
# whether the GWAS hit is itself the causal eQTL, which coloc.bf_bf cannot
# express, so the scoring is ours. Without one it is plain coloc.
# @noRd
.colocRunPairEither <- function(
    aligned,
    p12Info,
    q,
    gInfo,
    useEnrichment,
    ColocOptions
) {
    if (useEnrichment) {
        return(.colocRunPairEnloc(aligned, p12Info))
    }
    .colocRunPair(
        aligned,
        p12Info$p12Used,
        q,
        gInfo,
        p1 = p12Info$p1,
        p2 = p12Info$p2,
        ColocOptions = ColocOptions
    )
}

# Score one (QTL, GWAS) pair via coloc.bf_bf -> a summary row, or NULL when the
# variants don't align or coloc fails / returns no summary.
# @noRd
.colocScorePair <- function(
    qLbf,
    gInfo,
    q,
    useEnrichment,
    enrichment,
    priors,
    alleleFlip,
    ColocOptions
) {
    # Align variants between the QTL and GWAS LBF matrices by (chrom, pos,
    # allele) tuple via matchVariants (see .colocAlignLbf).
    aligned <- .colocAlignLbf(qLbf, gInfo$lbf, alleleFlip = alleleFlip)
    if (is.null(aligned)) {
        return(NULL)
    }
    p12Info <- .colocResolvePriors(
        gInfo,
        q,
        useEnrichment = useEnrichment,
        enrichment = enrichment,
        priors = priors
    )
    pairRes <- .colocRunPairEither(
        aligned,
        p12Info,
        q,
        gInfo,
        useEnrichment = useEnrichment,
        ColocOptions = ColocOptions
    )
    if (is.null(pairRes) || is.null(pairRes$summary)) {
        return(NULL)
    }
    rows <- .colocSummaryRow(pairRes, q, gInfo, useEnrichment, p12Info)
    # $results is the per-variant layer that process_coloc_results() used to
    # consume and that this pipeline silently dropped. It is pivoted here, the
    # only place that knows which results column belongs to which summary row.
    list(
        rows = rows,
        variants = .crPivotColocResults(pairRes$results, nrow(rows))
    )
}

# The priors for one pair. Without an enrichment table these are the caller's
# own; with one they are enloc's, taken from the table rather than derived
# from a default -- which is the whole point of running in enrichment mode.
# @noRd
.colocResolvePriors <- function(gInfo, q, useEnrichment, enrichment, priors) {
    if (!useEnrichment) {
        return(list(
            enRow = NA_real_,
            p1 = priors$p1,
            p2 = priors$p2,
            p12Used = priors$p12
        ))
    }
    idx <- .colocLookupEnrichment(enrichment, gInfo, q)
    if (is.na(idx)) {
        msg <- glue(
            "colocPipeline: no enrichment entry for ",
            "(gwasStudy='{gInfo$study}', qtlStudy='{q$study}', ",
            "qtlContext='{q$context}'); this pair is scored with no ",
            "enrichment, i.e. the unconditional priors."
        )
        warn(msg)
        return(list(
            enRow = NA_real_,
            p1 = priors$p1,
            p2 = priors$p2,
            p12Used = priors$p12
        ))
    }
    en <- .colocEnrichmentPriors(enrichment, idx)
    list(
        enRow = en$enrichment,
        p1 = en$p1,
        p2 = en$p2,
        p12Used = en$p12
    )
}

# Score every (QTL effect, GWAS effect) pair with the enloc configuration
# probabilities, in the shape coloc.bf_bf returns so the two modes' rows line
# up. RCP and LCP are reported alongside the hypothesis probabilities they
# are: RCP is PP.H4 and LCP is PP.H3 + PP.H4.
# @noRd
.colocRunPairEnloc <- function(aligned, p12Info) {
    priors <- list(
        p1 = p12Info$p1,
        p2 = p12Info$p2,
        p12 = p12Info$p12Used
    )
    qtlRows <- seq_len(nrow(aligned$qtl))
    gwasRows <- seq_len(nrow(aligned$gwas))
    grid <- expand.grid(idx1 = qtlRows, idx2 = gwasRows)
    scored <- map2(
        grid$idx1,
        grid$idx2,
        .colocEnlocOneEffectPair,
        aligned = aligned,
        priors = priors
    )
    summary <- data.frame(
        idx1 = grid$idx1,
        idx2 = grid$idx2,
        nSnps = ncol(aligned$qtl),
        PP.H0.abf = map_dbl(scored, "ppH0"),
        PP.H1.abf = map_dbl(scored, "ppH1"),
        PP.H2.abf = map_dbl(scored, "ppH2"),
        PP.H3.abf = map_dbl(scored, "ppH3"),
        PP.H4.abf = map_dbl(scored, "ppH4"),
        RCP = map_dbl(scored, "rcp"),
        LCP = map_dbl(scored, "lcp"),
        stringsAsFactors = FALSE
    )
    list(summary = summary, results = .colocEnlocResults(scored, aligned))
}

# One (QTL effect, GWAS effect) pair, which is enloc's signal cluster.
# @noRd
.colocEnlocOneEffectPair <- function(i, j, aligned, priors) {
    .enlocScoreEffectPair(
        qtlLog10Bf = aligned$qtl[i, ],
        gwasLog10Bf = aligned$gwas[j, ],
        priors = priors
    )
}

# The per-variant layer, matching what coloc.bf_bf's `results` carries: one
# column per scored pair, holding that pair's SNP-level colocalisation
# probabilities (enloc's SCP, which sums to RCP).
# @noRd
.colocEnlocResults <- function(scored, aligned) {
    scp <- map(scored, "scp")
    # A data.frame rather than a tibble: this stands in for what
    # coloc.bf_bf's `results` carries, and .crPivotColocResults consumes both
    # paths' output through the same code.
    out <- as.data.frame(
        exec(cbind, !!!scp),
        stringsAsFactors = FALSE
    )
    colnames(out) <- sprintf("SNP.PP.H4.row%d", seq_along(scp))
    mutate(out, snp = colnames(aligned$qtl))
}

# The enrichment columns, present only when the pipeline used an enrichment
# table. A NULL `p12Info` carries the zero-row schema an empty result needs.
# @noRd
.colocEnrichmentCols <- function(useEnrichment, p12Info) {
    if (!useEnrichment) {
        return(list())
    }
    if (is.null(p12Info)) {
        return(list(enrichment = numeric(0), p12Used = numeric(0)))
    }
    list(enrichment = p12Info$enRow, p12Used = p12Info$p12Used)
}

# Which QTL tuple and which GWAS tuple a summary row belongs to.
# @noRd
.colocRowIdentity <- function(q, gInfo) {
    list(
        study = q$study,
        context = q$context,
        trait = q$trait,
        method = q$method,
        gwasStudy = gInfo$study,
        gwasContext = gInfo$context,
        gwasTrait = gInfo$trait,
        gwasMethod = gInfo$method
    )
}

# Build a coloc summary row carrying the QTL / GWAS identity + enrichment.
# @noRd
.colocSummaryRow <- function(
    pairRes,
    q,
    gInfo,
    useEnrichment,
    p12Info
) {
    sm <- .colocRenameNsnps(
        as.data.frame(pairRes$summary, stringsAsFactors = FALSE)
    )
    mutate(
        sm,
        !!!.colocRowIdentity(q, gInfo),
        # idx1 / idx2 index the LBF rows handed to coloc.bf_bf, which is
        # exactly what retainedMass runs parallel to -- so the mass reported
        # here is the mass of the two effects this row actually scores, not a
        # per-entry average.
        qtlRetainedMass = .colocPickAt(
            q$retainedMass,
            sm[["idx1"]],
            nrow(sm)
        ),
        gwasRetainedMass = .colocPickAt(
            gInfo$retainedMass,
            sm[["idx2"]],
            nrow(sm)
        ),
        # coloc's idx1 / idx2 number the rows of THIS call, so they are not
        # comparable across blocks. The fit's own effect indices are, and they
        # are what the credible-set and gene views group on.
        qtlCs = .colocPickAt(
            q$effect,
            sm[["idx1"]],
            nrow(sm),
            fill = NA_integer_
        ),
        gwasCs = .colocPickAt(
            gInfo$effect,
            sm[["idx2"]],
            nrow(sm),
            fill = NA_integer_
        ),
        blockId = gInfo$blockId %||% NA_character_,
        !!!.colocEnrichmentCols(useEnrichment, p12Info)
    )
}

# Assemble the result table + attach the GWAS fine-mapping when requested.
# @noRd
.colocFinalize <- function(
    results,
    useEnrichment,
    qtlFineMappingResult,
    gwasFmr,
    gwasInput,
    returnGwasFineMapping
) {
    out <- .colocAssemble(
        results,
        useEnrichment,
        .colocLdSketch(qtlFineMappingResult, gwasFmr)
    )
    .colocAttachGwasFm(out, returnGwasFineMapping, gwasInput, gwasFmr)
}

# The GWAS fine-mapping the pipeline produced, carried back on the result when
# the caller asked for it and the GWAS side was fine-mapped here.
# @noRd
.colocAttachGwasFm <- function(out, returnGwasFineMapping, gwasInput, gwasFmr) {
    if (!returnGwasFineMapping || !methods::is(gwasInput, "SumStatsBase")) {
        return(out)
    }
    `attr<-`(out, "gwasFineMapping", gwasFmr)
}

# Row-bind + column-order the per-pair summary rows (empty result when none).
# @noRd
.colocAssemble <- function(results, useEnrichment, ldSketch = NULL) {
    if (length(results) == 0L) {
        return(.colocEmptyResult(enriched = useEnrichment, ldSketch = ldSketch))
    }
    rows <- bind_rows(map(map(results, "rows"), .colocStandardiseRow))
    variants <- list_c(map(results, "variants"))
    ColocResult(.colocOrderColumns(rows, useEnrichment), variants, ldSketch)
}

# Identity columns first, then everything coloc produced.
# @noRd
.colocOrderColumns <- function(out, useEnrichment) {
    idCols <- c(
        "study",
        "context",
        "trait",
        "method",
        "gwasStudy",
        "gwasContext",
        "gwasTrait",
        "gwasMethod",
        "blockId",
        "qtlCs",
        "gwasCs",
        if (useEnrichment) c("enrichment", "p12Used")
    )
    select(out, all_of(idCols), everything())
}

# =============================================================================
# Internal helpers
# =============================================================================

# LD-sketch compatibility check. Thin wrapper over the shared
# `.requireMatchingLdSketches` helper (R/ld.R). Shared with
# qtlEnrichmentPipeline.
# @noRd
.colocRequireMatchingLdSketches <- function(qtlLd, gwasLd) {
    .requireMatchingLdSketches(qtlLd, gwasLd, pipelineName = "colocPipeline")
}

# SuSiE credible-set concentration filter. Given a trimmed SuSiE fit
# and a coverage level, return the L-effect indices whose credible set
# is "narrow enough" to be informative: |CS| < nVariants * coverage *
# concentration. With concentration = 0.5 a 50% CS is kept only if it
# spans fewer than 25% of the locus variants -- this prunes diffuse
# signals before they reach coloc.bf_bf.
#
# Returns an integer vector of kept effect indices (empty when nothing
# survives), or errors when susieR is unavailable.
# @noRd
#' @importFrom susieR susie_get_cs
#' @importFrom purrr map_lgl
.colocFilterCsByConcentration <- function(
    fit,
    coverage = 0.5,
    concentration = 0.5
) {
    # V zapped to disable V-based filtering inside susie_get_cs
    unfiltered <- list_modify(fit, V = zap())
    csList <- susie_get_cs(unfiltered, coverage = coverage, dedup = FALSE)
    totalVariants <- ncol(fit$alpha)
    maxSize <- totalVariants * coverage * concentration
    keep <- map_lgl(csList$cs, .colocCsUnderMax, maxSize = maxSize)
    as.numeric(str_remove_all(names(which(keep)), "L"))
}

# Extract an LBF matrix (effects x variants) from a FineMappingRow,
# applying the same filtering knobs as the legacy .extractLbfMatrix, now
# carried as one `lbfFilter` bundle rather than three loose arguments:
#   - filterLbfCs (CS-only)
#   - secondary (secondary coverage CS, with a concentration cutoff)
#   - priorTol drop on V (default)
# Handles the fSuSiE shape (where the LBF lives at a different slot).
# Returns list(lbf = <matrix>, variantIds = <character>) or NULL when
# the entry has no usable LBF matrix.
# @noRd
.colocExtractLbfFromEntry <- function(
    parts,
    lbfFilterArgs,
    label = "entry"
) {
    fit <- susieFit(parts)
    if (is.null(fit)) {
        msg <- glue("colocPipeline: {label} has no trimmedFit; skipping.")
        warn(msg)
        return(NULL)
    }
    allRows <- .colocLbfMatrix(fit, label)
    if (is.null(allRows)) {
        return(NULL)
    }
    allMass <- .colocEffectRetainedMass(fit, nrow(allRows))
    keep <- .colocSelectLbfRows(allRows, fit, lbfFilterArgs)
    kept <- allRows[keep, , drop = FALSE]
    mass <- allMass[keep]
    if (nrow(kept) == 0L) {
        return(NULL)
    }
    lbfMatrix <- .colocAssignLbfColnames(kept, parts)
    if (ncol(lbfMatrix) == 0L) {
        return(NULL)
    }
    list(
        lbf = lbfMatrix,
        variantIds = colnames(lbfMatrix),
        retainedMass = mass,
        # The fit's own effect indices, which -- unlike coloc's idx1 / idx2 --
        # are stable across calls and blocks, so they are what identifies a
        # credible set in the result.
        effect = keep
    )
}

# Extract the (effects x variants) LBF matrix from a fit: susie lbf_variable, or
# a stacked fSuSiE lBF (direct or nested), or NULL (with a warning) when absent
# / empty.
# @noRd
.colocLbfMatrix <- function(fit, label) {
    lbfMatrix <- if (!is.null(fit$lbf_variable)) {
        as.matrix(fit$lbf_variable)
    } else if (
        !is.null(fit$fsusie_result) &&
            is.list(fit$fsusie_result$lBF)
    ) {
        # fSuSiE path: stack per-trait lBF lists into a single matrix.
        lbfList <- fit$fsusie_result$lBF
        exec(rbind, !!!lbfList)
    } else if (
        is.list(fit) &&
            length(fit) >= 1L &&
            !is.null(fit[[1L]]$fsusie_result$lBF)
    ) {
        lbfList <- fit[[1L]]$fsusie_result$lBF
        exec(rbind, !!!lbfList)
    } else {
        msg <- glue(
            "colocPipeline: {label} trimmedFit has no lbf_variable / fsusie ",
            "lBF; skipping."
        )
        warn(msg)
        return(NULL)
    }
    if (is.null(lbfMatrix) || nrow(lbfMatrix) == 0L) {
        msg <- glue("colocPipeline: {label} LBF matrix is empty.")
        warn(msg)
        return(NULL)
    }
    lbfMatrix
}

# Row (effect) selection, in the original priority order: primary CS index,
# else secondary-CS-by-concentration, else prior-variance threshold.
#
# Returns the surviving row INDICES rather than the subsetted matrix. The
# caller has a second vector running parallel to those rows (the per-effect
# retained mass), and handing back a selector keeps the two narrowed by one
# decision instead of two hand-synchronised subsets.
#
# `fit` is read with `[[`, never `$`: `$` on a list falls back to prefix
# matching when the exact name is absent, so a fit carrying `sets_secondary`
# but no `sets` would silently filter on the wrong element.
# @noRd
.colocSelectLbfRows <- function(lbfMatrix, fit, lbfFilterArgs) {
    # Read off the bundle like the other three filter fields, rather than
    # taken as a second copy of a value the caller already has in hand.
    priorTol <- lbfFilterArgs$priorTol %||% 1e-9
    # The only place the filter's fields are read. Every function between
    # here and colocPipeline carries the bundle whole, so the names cannot
    # drift apart on the way down.
    allRows <- seq_len(nrow(lbfMatrix))
    if (isTRUE(lbfFilterArgs$filterLbfCs) && is.null(lbfFilterArgs$secondary)) {
        csIdx <- fit[["sets"]][["cs_index"]]
        if (!is.null(csIdx) && length(csIdx) > 0L) {
            return(csIdx)
        }
    } else if (!is.null(lbfFilterArgs$secondary)) {
        secIdx <- try_fetch(
            .colocFilterCsByConcentration(
                fit,
                coverage = lbfFilterArgs$secondary,
                concentration = lbfFilterArgs$concentration
            ),
            error = function(cnd) NULL
        )
        if (!is.null(secIdx) && length(secIdx) > 0L) {
            return(secIdx)
        }
    } else if (!is.null(fit[["V"]])) {
        return(allRows[fit[["V"]] > priorTol])
    }
    allRows
}

# Per-effect retained posterior mass, parallel to the LBF rows.
#
# `retained_mass` is written by adjustPips() at reconciliation time and records
# how much of each effect's alpha survived the shared-variant restriction. It
# is absent when no reconciliation ran (adjustPips = FALSE), which is reported
# as NA -- "not measured", distinct from a measured 0.
# @noRd
.colocEffectRetainedMass <- function(fit, nEffects) {
    mass <- fit[["retained_mass"]]
    if (is.null(mass) || length(mass) != nEffects) {
        return(rep(NA_real_, nEffects))
    }
    as.numeric(mass)
}

# The entry's variant ids, when the matrix does not already name its columns.
# @noRd
.colocLbfNamesFrom <- function(lbfMatrix, parts) {
    if (!is.null(colnames(lbfMatrix)) && !any(is.na(colnames(lbfMatrix)))) {
        return(lbfMatrix)
    }
    vids <- .fmrPartsVariantIds(parts)
    if (length(vids) != ncol(lbfMatrix)) {
        return(lbfMatrix)
    }
    `colnames<-`(lbfMatrix, vids)
}

# Assign variant-id column names (fit-provided, else the row's rendered
# variant ids) and drop columns with an NA id.
# @noRd
.colocAssignLbfColnames <- function(lbfMatrix, parts) {
    named <- .colocLbfNamesFrom(lbfMatrix, parts)
    named[, !is.na(colnames(named)), drop = FALSE]
}

# The second side's LBF matrices, one record per row of the collection, each
# carrying the row's identity so the pair it scores can be named.
#
# One record per ROW rather than per identity key: a QTL second side has many
# rows sharing (study, method, block) and differing only on context / trait, so
# keying on the GWAS 2-tuple would let one trait's matrix silently replace
# another's. Nothing downstream indexes this list by name, so positional
# records make the collision impossible instead of merely unlikely.
# @noRd
.colocPreextractGwasLbf <- function(gwasFmr, lbfFilterArgs) {
    if (nrow(gwasFmr) == 0L) {
        return(list())
    }
    compact(map(
        seq_len(nrow(gwasFmr)),
        .colocGwasLbfAt,
        gwasFmr = gwasFmr,
        blockIds = .colocGwasBlockIds(gwasFmr),
        side = .fmrSideName(gwasFmr),
        lbfFilterArgs = lbfFilterArgs
    ))
}

# One second-side row's LBF matrix plus its identity, or NULL when the row has
# no usable LBF.
# @noRd
.colocGwasLbfAt <- function(
    ri,
    gwasFmr,
    blockIds,
    side,
    lbfFilterArgs
) {
    ident <- .colocTupleIdentity(gwasFmr, ri)
    blockId <- blockIds[[ri]]
    label <- .fmrTupleLabel(side, ident, block = blockId)
    info <- .colocExtractLbfFromEntry(
        .fmrRowParts(gwasFmr, ri),
        lbfFilterArgs,
        label = label
    )
    if (is.null(info)) {
        return(NULL)
    }
    c(
        ident,
        list(
            lbf = info$lbf,
            retainedMass = info$retainedMass,
            effect = info$effect,
            blockId = blockId,
            label = label
        )
    )
}

# The LD block each GWAS fine-mapping row was computed on. The element's own
# range is the identity; an explicit `blockId` (which keys the external block
# manifest, and so carries the true block BOUNDARIES rather than the span of
# the variants that survived) is preferred when present.
# @noRd
.colocGwasBlockIds <- function(gwasFmr) {
    md <- mcols(gwasFmr, use.names = FALSE)
    if (is_in("blockId", colnames(md))) {
        return(as.character(md$blockId))
    }
    .rtlRangeKeys(gwasFmr)
}

# Align column names between a QTL and a GWAS LBF matrix by
# (chrom, pos, allele) tuple (via matchVariants), returning both restricted to
# the common variants under one shared id so coloc.bf_bf can align them.
# @noRd
.colocAlignLbf <- function(qtlLbf, gwasLbf, alleleFlip = TRUE) {
    qids <- colnames(qtlLbf)
    gids <- colnames(gwasLbf)
    # Match LBF columns by (chrom, pos, allele) tuple rather than by raw id
    # string, so chr-prefix / separator / allele-order differences resolve. LBF
    # is allele-coding-invariant, so this is pure identity alignment (no sign);
    # with alleleFlip = FALSE, ref/alt-swapped columns are not treated as
    # shared.
    m <- matchVariants(qids, gids, allowFlip = alleleFlip)
    if (length(m$idxA) == 0L) {
        return(NULL)
    }
    # Relabel both matrices to one shared id so coloc.bf_bf sees identical
    # names.
    sharedIds <- gids[m$idxB]
    qSub <- `colnames<-`(qtlLbf[, m$idxA, drop = FALSE], sharedIds)
    gSub <- `colnames<-`(gwasLbf[, m$idxB, drop = FALSE], sharedIds)
    list(qtl = qSub, gwas = gSub)
}

# A blank result data frame for the no-pair case so callers downstream
# do not have to special-case a NULL return. When `enriched = TRUE` the
# enrichment + p12Used columns are appended (the enloc-mode schema).
# @noRd
.colocEmptyResult <- function(enriched = FALSE, ldSketch = NULL) {
    base <- tibble(
        study = character(0),
        context = character(0),
        trait = character(0),
        method = character(0),
        gwasStudy = character(0),
        gwasContext = character(0),
        gwasTrait = character(0),
        gwasMethod = character(0),
        blockId = character(0),
        qtlCs = integer(0),
        gwasCs = integer(0),
        idx1 = integer(0),
        idx2 = integer(0),
        nSnps = integer(0),
        PP.H0.abf = numeric(0),
        PP.H1.abf = numeric(0),
        PP.H2.abf = numeric(0),
        PP.H3.abf = numeric(0),
        PP.H4.abf = numeric(0),
        qtlRetainedMass = numeric(0),
        gwasRetainedMass = numeric(0)
    )
    ColocResult(
        mutate(base, !!!.colocEnrichmentCols(enriched, NULL)),
        list(),
        ldSketch = ldSketch
    )
}

# Look up this pair's enrichment factor in the user-supplied enrichment table.
#
# The join uses whichever identity columns the table carries: the
# (gwasStudy, qtlStudy, qtlContext) triple always, plus gwasContext /
# gwasTrait when qtlEnrichmentPipeline ran with a QTL outcome side. Without
# those two, one study's molecular traits are indistinguishable in the table
# and every one of them would take the same row. Returns NA when the pair is
# not present; the caller falls back to the baseline p12 and warns.
# @noRd
.colocLookupEnrichment <- function(enrichment, gInfo, q) {
    wanted <- .colocEnrichmentKey(enrichment, gInfo, q)
    hits <- map(
        names(wanted),
        .colocEnrichmentColumnMatches,
        enrichment = enrichment,
        wanted = wanted
    )
    idx <- which(reduce(hits, `&`))
    if (length(idx) == 0L) {
        return(NA_integer_)
    }
    idx[[1L]]
}

# The enloc priors for one pair, read from the enrichment table rather than
# derived here. qtlEnrichmentPipeline already computes them the way fastenloc
# does -- p1 = (1 - P_eqtl) * expit(a0), p2 = P_eqtl / (1 + exp(a0 + a1)),
# p12 = P_eqtl * expit(a0 + a1) -- from the genome-wide QTL results, so using
# them is what makes this an enloc run rather than an approximation of one.
# @noRd
.colocEnrichmentPriors <- function(enrichment, idx) {
    needed <- c("colocP1", "colocP2", "colocP12")
    missing <- setdiff(needed, colnames(enrichment))
    if (length(missing) > 0L) {
        abort(glue(
            "colocPipeline: `enrichment` is missing ",
            "{str_flatten(missing, ', ')}. The enloc priors come from these ",
            "columns, which qtlEnrichmentPipeline() emits; a table carrying ",
            "only an `enrichment` factor is not enough to form them."
        ))
    }
    list(
        p1 = as.numeric(enrichment$colocP1[[idx]]),
        p2 = as.numeric(enrichment$colocP2[[idx]]),
        p12 = as.numeric(enrichment$colocP12[[idx]]),
        enrichment = as.numeric(enrichment$enrichment[[idx]])
    )
}

# The identity a pair is looked up by, narrowed to the columns the table has.
# @noRd
.colocEnrichmentKey <- function(enrichment, gInfo, q) {
    wanted <- list(
        gwasStudy = gInfo$study,
        gwasContext = gInfo$context,
        gwasTrait = gInfo$trait,
        qtlStudy = q$study,
        qtlContext = q$context
    )
    wanted[is_in(names(wanted), colnames(enrichment))]
}

# One key column's row match. An axis neither side has is NA on both, and
# matches -- `==` would evaluate to NA there and drop every row.
# @noRd
.colocEnrichmentColumnMatches <- function(column, enrichment, wanted) {
    values <- as.character(enrichment[[column]])
    if (is.na(wanted[[column]])) {
        return(is.na(values))
    }
    !is.na(values) & values == wanted[[column]]
}

# Ensure each row data.frame from coloc.bf_bf carries the standard PP
# columns even when the underlying call produced a slightly different
# shape.
# @noRd
.colocStandardiseRow <- function(sm) {
    missing <- setdiff(
        c(
            "idx1",
            "idx2",
            "nSnps",
            "PP.H0.abf",
            "PP.H1.abf",
            "PP.H2.abf",
            "PP.H3.abf",
            "PP.H4.abf"
        ),
        colnames(sm)
    )
    mutate(sm, !!!set_names(map(missing, .colocMissingColumn), missing))
}

# @noRd
.colocMissingColumn <- function(col) {
    NA
}

# coloc.bf_bf names the aligned-variant count `nsnps`; this package publishes
# it as `nSnps`. Without the rename .colocStandardiseRow() invents an all-NA
# `nSnps` and the real count is left sitting in a second column beside it.
# @noRd
.colocRenameNsnps <- function(sm) {
    if (is_in("nsnps", colnames(sm)) && !is_in("nSnps", colnames(sm))) {
        sm <- rename(sm, nSnps = "nsnps")
    }
    sm
}

# Index a per-effect vector (retained mass, effect id) by coloc's effect index,
# tolerating an index coloc did not report (NA) or one past the end.
# @noRd
.colocPickAt <- function(values, idx, n, fill = NA_real_) {
    if (is.null(idx) || is.null(values)) {
        return(rep(fill, n))
    }
    idx <- as.integer(idx)
    ok <- !is.na(idx) & idx >= 1L & idx <= length(values)
    replace(rep(fill, length(idx)), ok, values[idx[ok]])
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# The variant ids of one fine-mapping entry (S4 slot; not pluckable by name).
# @noRd

# Score the first side's LBF against one second-side record -> a summary row
# (or NULL).
# @noRd
.colocScorePairAt <- function(
    gInfo,
    qLbfInfo,
    q,
    useEnrichment,
    enrichment,
    priors,
    alleleFlip,
    ColocOptions
) {
    .colocScorePair(
        qLbfInfo$lbf,
        gInfo,
        q,
        useEnrichment = useEnrichment,
        enrichment = enrichment,
        priors = priors,
        alleleFlip = alleleFlip,
        ColocOptions = ColocOptions
    )
}

# TRUE when a credible set has fewer than `maxSize` variants.
# @noRd
.colocCsUnderMax <- function(x, maxSize) {
    length(x) < maxSize
}

# --- GwasFineMappingParam accessors ----------------------------------------
# The nested getters return the Param itself, so a caller edits it with its
# own replacement methods and hands it back:
#   cs <- credibleSetParam(g); csCoverage(cs) <- 0.9; credibleSetParam(g) <- cs

#' @rdname GwasFineMappingParam
setMethod("fineMappingMethods", "GwasFineMappingParam", function(x) {
    x@methods
})

#' @rdname GwasFineMappingParam
setReplaceMethod(
    "fineMappingMethods",
    "GwasFineMappingParam",
    function(x, value) {
        x@methods <- value
        validObject(x)
        x
    }
)

#' @rdname GwasFineMappingParam
setMethod("credibleSetParam", "GwasFineMappingParam", function(x) {
    x@credibleSetParam
})

#' @rdname GwasFineMappingParam
setReplaceMethod(
    "credibleSetParam",
    "GwasFineMappingParam",
    function(x, value) {
        x@credibleSetParam <- value
        validObject(x)
        x
    }
)

#' @rdname GwasFineMappingParam
setMethod("susieRssParam", "GwasFineMappingParam", function(x) x@susieRssParam)

#' @rdname GwasFineMappingParam
setReplaceMethod("susieRssParam", "GwasFineMappingParam", function(x, value) {
    x@susieRssParam <- value
    validObject(x)
    x
})

#' @rdname GwasFineMappingParam
setMethod("panelFilterParam", "GwasFineMappingParam", function(x) {
    x@panelFilterParam
})

#' @rdname GwasFineMappingParam
setReplaceMethod(
    "panelFilterParam",
    "GwasFineMappingParam",
    function(x, value) {
        x@panelFilterParam <- value
        validObject(x)
        x
    }
)

#' @rdname GwasFineMappingParam
setMethod("initializeWithSusieInf", "GwasFineMappingParam", function(x) {
    x@initializeWithSusieInf
})

#' @rdname GwasFineMappingParam
setReplaceMethod(
    "initializeWithSusieInf",
    "GwasFineMappingParam",
    function(x, value) {
        x@initializeWithSusieInf <- value
        validObject(x)
        x
    }
)

#' @rdname GwasFineMappingParam
setMethod("fitRetention", "GwasFineMappingParam", function(x) x@fitRetention)

#' @rdname GwasFineMappingParam
setReplaceMethod("fitRetention", "GwasFineMappingParam", function(x, value) {
    x@fitRetention <- value
    validObject(x)
    x
})

#' @title Arguments For coloc's Bayes-Factor Scorer
#' @description Options forwarded to \code{coloc::coloc.bf_bf}. The prior
#'   probabilities are \emph{not} settable here --- pecotmr derives them from
#'   \code{\link{ColocPriorParam}}, which its enrichment adjustment also
#'   reads --- and the two Bayes-factor matrices are supplied by the pipeline.
#' @param overlap.min Minimum overlap between the two variant sets.
#' @param trim_by_posterior Logical; trim by posterior before scoring.
#' @param ... Any other \code{coloc::coloc.bf_bf} argument.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' ColocOptions(trim_by_posterior = FALSE)
#' @export
ColocOptions <- function(overlap.min = NULL, trim_by_posterior = NULL, ...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            bf1 = "supplied from the data by the pipeline",
            bf2 = "supplied from the data by the pipeline",
            p1 = "`colocP1` on the pipeline",
            p2 = "`colocP2` on the pipeline",
            p12 = "`colocP12` on the pipeline"
        ),
        "ColocOptions"
    )
    .newMethodOptions(
        "coloc::coloc.bf_bf",
        defaults = list(
            overlap.min = overlap.min,
            trim_by_posterior = trim_by_posterior
        ),
        extra = list(...),
        label = "ColocOptions",
        engine = "coloc"
    )
}

# Run coloc.bf_bf for an aligned pair, warning + NULL on failure.
# @noRd
#' @importFrom rlang try_fetch
.colocRunPair <- function(aligned, p12Used, q, gInfo, p1, p2, ColocOptions) {
    # Engine-call exception: coloc's whole interface is two blocks -- this and
    # ColocOptions -- so a colocWrapper.R would hold ~50 lines against this
    # file's 1600. The engine boundary is stated here instead of split out.
    callArgs <- c(
        list(
            aligned$qtl,
            aligned$gwas,
            p1 = p1,
            p2 = p2,
            p12 = p12Used
        ),
        ColocOptions
    )
    try_fetch(
        exec(coloc::coloc.bf_bf, !!!callArgs),
        error = function(cnd) {
            msg <- glue(
                "colocPipeline: coloc.bf_bf failed for ",
                "{q$label} x {gInfo$label}"
            )
            warn(msg, parent = cnd)
            NULL
        }
    )
}
