#' @include MethodParam.R
NULL

#' @rdname CtwasPriorParam
#' @aliases CtwasPriorParam-class
#' @exportClass CtwasPriorParam
setClass(
    "CtwasPriorParam",
    contains = "MethodParam",
    slots = c(
        thin = "numeric",
        niterPrefit = "numeric",
        niter = "numeric",
        varStructure = "character",
        fallbackToPrefit = "logical"
    )
)

#' @title cTWAS Prior-Estimation Settings
#' @description The EM stage that estimates cTWAS's group priors
#'   (\code{group_prior} / \code{group_prior_var}). These are exactly the
#'   settings \code{\link{estCtwasGroupPriors}} takes, so the bundle travels
#'   from \code{\link{ctwasPipeline}} to that function whole.
#' @param thin Numeric. Proportion of SNPs retained when assembling region
#'   data for the EM. Default \code{0.1}. This is stochastic thinning for
#'   speed, unrelated to \code{\link{VariantPruningParam}}, which decides
#'   which variants enter a gene's weight matrix.
#' @param niterPrefit Integer. Iterations for the prefit EM. Default
#'   \code{3}.
#' @param niter Integer. Iterations for the accurate EM. Default \code{30}.
#' @param varStructure How the prior variance is shared across groups, a group
#'   being a study, context and molecular type (see
#'   \code{\link{ctwasPipeline}}). \code{"shared_nonSNP"} (default) pools it
#'   across the gene groups, so a group with few causal genes borrows from the
#'   others, and keeps the variants' apart. \code{"shared_all"} pools the
#'   variants in as well, \code{"shared_type"} pools the groups of each
#'   molecular type, \code{"shared_context"} those of each study and context,
#'   and \code{"independent"} gives each group its own. With one type per study
#'   and context, \code{"shared_context"} is \code{"independent"}.
#' @param fallbackToPrefit Logical. When the accurate EM returns NaN on
#'   underpowered data, fall back to the prefit estimates rather than failing.
#'   Default \code{FALSE}. A pecotmr recovery step, not a cTWAS argument.
#' @return A \code{CtwasPriorParam} object, a \code{\link{MethodParam}}.
#' @examples
#' CtwasPriorParam(thin = 1, niter = 10, fallbackToPrefit = TRUE)
#' @export
CtwasPriorParam <- function(
    thin = 0.1,
    niterPrefit = 3L,
    niter = 30L,
    varStructure = c(
        "shared_nonSNP",
        "shared_all",
        "shared_type",
        "shared_context",
        "independent"
    ),
    fallbackToPrefit = FALSE
) {
    varStructure <- arg_match(varStructure)
    new(
        "CtwasPriorParam",
        thin = thin,
        niterPrefit = niterPrefit,
        niter = niter,
        varStructure = varStructure,
        fallbackToPrefit = fallbackToPrefit
    )
}

#' @rdname VariantPruningParam
#' @aliases VariantPruningParam-class
#' @exportClass VariantPruningParam
setClass(
    "VariantPruningParam",
    contains = "MethodParam",
    slots = c(
        weightCutoff = "numeric",
        csMinCor = "numeric",
        pipRescue = "numeric",
        maxVariants = "numeric"
    )
)

#' @title Per-Gene Variant-Pruning Settings
#' @description Which variants survive into each gene's weight matrix before
#'   cTWAS sees it. \code{maxVariants} caps the count; \code{csMinCor} and
#'   \code{pipRescue} mark variants as must-keep so they survive the cap.
#' @param weightCutoff Numeric. Drop variants with
#'   \code{|weight| < weightCutoff}. Default \code{0} (no filter).
#' @param csMinCor Numeric. Variants in a 95\% credible set whose purity
#'   (\code{min_abs_corr}) is at least this are must-keep. Default
#'   \code{0.8}. Requires a \code{fineMappingResult}; ignored without one.
#' @param pipRescue Numeric. Variants with PIP above this are must-keep.
#'   Default \code{0} (no PIP rescue). Requires a \code{fineMappingResult}.
#' @param maxVariants Numeric. Cap on the per-gene variant count. Above it,
#'   all must-keep variants are retained and the remaining slots filled by
#'   descending PIP (when available) or descending \code{|weight|}. Default
#'   \code{Inf} (no cap).
#' @return A \code{VariantPruningParam} object, a \code{\link{MethodParam}}.
#' @examples
#' VariantPruningParam(weightCutoff = 1e-4, maxVariants = 5000)
#' @export
VariantPruningParam <- function(
    weightCutoff = 0,
    csMinCor = 0.8,
    pipRescue = 0,
    maxVariants = Inf
) {
    new(
        "VariantPruningParam",
        weightCutoff = weightCutoff,
        csMinCor = csMinCor,
        pipRescue = pipRescue,
        maxVariants = maxVariants
    )
}

#' @rdname BoundaryMergeParam
#' @aliases BoundaryMergeParam-class
#' @exportClass BoundaryMergeParam
setClass(
    "BoundaryMergeParam",
    contains = "MethodParam",
    slots = c(
        enabled = "logical",
        pipThreshold = "numeric",
        filterCs = "logical",
        maxSnp = "numeric"
    )
)

#' @title Boundary-Region Merge Settings
#' @description Whether a high-PIP gene whose weights cross an LD-block
#'   boundary has its regions merged and re-fine-mapped after each run
#'   (\code{\link{mergeCtwasBoundaryRegions}}). A first pass fine-maps such a
#'   gene with the LD of one block only, whatever its z-score drew from the
#'   other, and can call a gene whose real signal is a variant in that other
#'   block; merging is what removes those calls.
#' @param enabled Logical. Run the merge step. Default \code{TRUE}; the
#'   other three fields are read only when this is \code{TRUE}.
#' @param pipThreshold Numeric. PIP threshold for selecting which boundary
#'   genes to merge. Default \code{0.5}.
#' @param filterCs Logical. Require the boundary gene to be in a credible set
#'   to be selected. Default \code{FALSE}.
#' @param maxSnp Numeric. Per-merged-region SNP cap. Default \code{Inf}.
#' @return A \code{BoundaryMergeParam} object, a \code{\link{MethodParam}}.
#' @examples
#' BoundaryMergeParam(pipThreshold = 0.8)
#' @export
BoundaryMergeParam <- function(
    enabled = TRUE,
    pipThreshold = 0.5,
    filterCs = FALSE,
    maxSnp = Inf
) {
    new(
        "BoundaryMergeParam",
        enabled = enabled,
        pipThreshold = pipThreshold,
        filterCs = filterCs,
        maxSnp = maxSnp
    )
}

# Every bundle ctwasPipeline() takes, checked together. The argument names
# in the messages are the pipeline's own, which drop the `Args` suffix.
# @noRd
.ctwasAssertBundles <- function(
    ctwasPriorParam,
    variantPruningParam,
    boundaryMergeParam,
    methodArgs,
    weightSelectionParam
) {
    .assertMethodParam(ctwasPriorParam, "CtwasPriorParam")
    .assertMethodParam(
        variantPruningParam,
        "VariantPruningParam"
    )
    .assertMethodParam(
        boundaryMergeParam,
        "BoundaryMergeParam"
    )
    .assertMethodOptions(methodArgs, "CtwasOptions", "methodArgs")
    .assertMethodParam(weightSelectionParam, "WeightSelectionParam")
}

#' @title Causal TWAS Pipeline (cTWAS, multi LD block)
#' @description Pipeline that hands a per-block set of
#'   \code{\link{GwasSumStats}} of GWAS Z-scores together with the matching
#'   per-block per-gene TWAS weights and LD references to
#'   \code{ctwas::ctwas_sumstats}, producing per-gene posterior inclusion
#'   probabilities for causal genes. Optionally accepts a precomputed TWAS-Z
#'   \code{GRanges} from \code{\link{causalInferencePipeline}} as the
#'   \code{z_gene} input so the per-gene Z is not recomputed inside ctwas.
#'
#' @section LD block convention: \code{gwasSumStats} is ONE
#'   \code{\link{GwasSumStats}} whose elements are LD blocks, keyed by its
#'   \code{blockId} column -- build it with
#'   \code{loadGwasSumStatsFromManifest(..., region = <blocks>)}. Per-block
#'   \code{region_info}, \code{LD_map}, and \code{snp_map} entries are built
#'   automatically from the LD reference and concatenated before the call to
#'   \code{ctwas::ctwas_sumstats}. A single-block input is rejected: cTWAS's EM
#'   cannot converge on a single region, so callers must supply at least two
#'   blocks.
#'
#' @section Boundary genes: A gene is modelled as cTWAS models it. Its z-score
#'   reads all its weights, whichever blocks they fall in, its weight variance
#'   treating the blocks as independent; it is fine-mapped in the block holding
#'   the largest share of its absolute weights, with the part of its weights in
#'   that block. Such a first pass can call a gene whose real signal is a
#'   variant across the boundary, so by default a high-PIP gene crossing one has
#'   its regions merged and fine-mapped again (\code{boundaryMergeParam}).
#'
#' @section Groups and molecular types: cTWAS estimates one prior per group, and
#'   a group here is a study, context and molecular type, so two studies with
#'   one context name stay apart. The type is the \code{dataType} the weights
#'   record (\code{"QTL"} where they record none, all such contexts being one
#'   type), and it matters to the fit under \code{varStructure = "shared_type"},
#'   or where one study and context holds two types. Results name a group
#'   \code{{study}_{context}_{type}}, with the three also as columns of their
#'   own, and a gene \code{{trait}|{study}_{context}_{type}_{method}}, the form
#'   cTWAS reads the trait from. Groups whose names would coincide are refused,
#'   as is a name holding \code{"|"}, on which cTWAS's EM splits a group.
#'
#' @section cTWAS's log: cTWAS narrates each of its steps through the
#'   \code{logging} package, which writes to the console with \code{cat()},
#'   out of reach of \code{suppressMessages()}. The pipeline keeps that
#'   narration off the console, cTWAS's warnings excepted. To keep the full
#'   log, give a step a log file through its options, e.g.
#'   \code{CtwasEstParamOptions(logfile = "est_param.log")}.
#'
#' @section LD reference compatibility check: Per block: the block's GWAS
#'   rows must share one LD reference (a block's LD is one matrix), and every
#'   weight row placed in the block that carries a reference must be the same
#'   LD reference as theirs. Blocks may sit on different references.
#'   The two must be the same kind of reference. Two genotype panels must
#'   share their samples and allele orientation; two \code{LdMatrix} objects
#'   must agree on the LD they share; two \code{GenotypeMixture} objects must
#'   average the same panels with the same weights. The two need NOT carry the
#'   same variants; a partial overlap, which is what QC-ing the two sides
#'   separately produces, only warns. A mismatch of any of those, or no shared
#'   variant at all, is a hard error.
#'
#' @param gwasSumStats A \code{\link{GwasSumStats}} whose elements are LD
#'   blocks (at least two), keyed by its \code{blockId} column, with
#'   \code{qcOptions()} non-empty. Build it with
#'   \code{loadGwasSumStatsFromManifest(..., region = <blocks>)} and pass it
#'   through \code{\link{summaryStatsQc}}.
#' @param twasWeights The per-gene weight source. Either (a) a FLAT
#'   \code{\link{TwasWeights}} / \code{QtlFineMappingResult} (or a homogeneous
#'   list of them), each gene read with the LD block holding most of its
#'   weight; or (b) a pre-bucketed NAMED LIST keyed by \code{region_id} (keys a
#'   SUBSET of \code{gwasSumStats}'s), used as-is. A gene's \code{traitPos}, if
#'   it has one, gives its coordinates in the results, which are otherwise the
#'   span of its weight variants.
#'   Blocks without any TWAS weights do not feed the group priors by default:
#'   cTWAS skips regions with fewer genes than its \code{min_gene}, 1. To let
#'   their SNP-level signal in, set \code{min_gene = 0} in
#'   \code{\link{CtwasEstParamOptions}}. A
#'   \code{TwasWeights} row of the run's method that carries a SuSiE fit and
#'   was narrowed since that fit warns (condition class
#'   \code{pecotmrFitCoverage}): the fit no longer lines up with the row's
#'   variants, so its weights are not renormalized over the variants kept.
#' @param twasZ Optional \code{GRanges} of TWAS Z-scores (output of
#'   \code{\link{causalInferencePipeline}}). When supplied, each weighted
#'   gene's Z -- the TWAS row of its \code{(qtlStudy, context, trait,
#'   method)} for the GWAS study of \code{gwasSumStats} -- is used as the
#'   \code{z_gene} input to \code{ctwas_sumstats} so it is not recomputed.
#'   Every weighted gene needs exactly one, so \code{twasZ} must hold the
#'   method each gene is run with: chosen with the same
#'   \code{weightSelectionParam}, or a superset such as \code{methods =
#'   "all"}. Rows for other genes, methods or GWAS studies, and MR rows, are
#'   ignored.
#' @param fineMappingResult Optional \code{QtlFineMappingResult}: the genes'
#'   own fine-mapping, whose per-variant PIP and credible-set membership the
#'   CS / PIP rescue filters (\code{csMinCor} and \code{minPipCutoff}) read.
#'   Each gene is looked up by its weights' \code{(study, context, trait)} in
#'   the fits of \code{fineMappingMethod}; a gene with no fit there gets no
#'   rescue. When \code{NULL} (default) the smart filters are no-ops; only the
#'   magnitude filter (\code{twasWeightCutoff}) and the per-gene cap
#'   (\code{maxNumVariants}, ordered by \code{|weight|}) apply.
#' @param fineMappingMethod Character (length 1): the method of
#'   \code{fineMappingResult} whose fits the rescues read, e.g.
#'   \code{"susie"}. Required with a \code{fineMappingResult}, since the
#'   weights' own method -- an \code{"ensemble"} or \code{"lasso"} weight --
#'   names no fine-mapping fit. Default \code{NULL}.
#' @param weightSelectionParam Which method's weights each gene (each
#'   \code{(study, context, trait)}) is modelled with, built with
#'   \code{\link{WeightSelectionParam}}, as in
#'   \code{\link{causalInferencePipeline}}. A cTWAS run models one weight per
#'   gene. By default (\code{methods = "best"}) one run takes each gene's
#'   best method by cross-validated r-squared, the genes' methods mixed, and
#'   its results are labelled \code{method = "best"}; the gene ids keep each
#'   gene's own method. A gene whose methods carry no r-squared to pick one
#'   by is then refused, unless it has just one. \code{methods = "all"}, or
#'   several named methods, runs cTWAS once per method, one after another,
#'   the results labelled by method in one \code{CtwasResult}; a single
#'   named method is one run.
#' @param maxNumSingleEffects Integer. Max number of single effects in cTWAS
#'   fine-mapping. Default \code{5}. Pass-through to
#'   \code{ctwas::ctwas_sumstats}.
#' @param BPPARAM A \code{BiocParallelParam} whose workers number the cores
#'   ctwas's own \code{ncore} parallelism uses when they run on this machine
#'   (\code{SerialParam}, \code{MulticoreParam}); any other backend gives one
#'   core. Default \code{SerialParam()}.
#' @param ctwasPriorParam The EM stage that estimates cTWAS's group priors,
#'   built with \code{\link{CtwasPriorParam}}: \code{thin},
#'   \code{niterPrefit}, \code{niter}, \code{varStructure} and
#'   \code{fallbackToPrefit}. Handed whole to
#'   \code{\link{estCtwasGroupPriors}}.
#' @param variantPruningParam Which variants survive into each gene's weight
#'   matrix, built with \code{\link{VariantPruningParam}}:
#'   \code{weightCutoff}, \code{csMinCor}, \code{pipRescue} and
#'   \code{maxVariants}. Handed whole to
#'   \code{\link{assembleCtwasInputs}}. The \code{csMinCor} /
#'   \code{pipRescue} rescues need a \code{fineMappingResult}.
#' @param boundaryMergeParam Whether a high-PIP gene straddling an LD-block
#'   boundary has its regions merged and re-fine-mapped, built with
#'   \code{\link{BoundaryMergeParam}}: \code{enabled} (default
#'   \code{TRUE}), \code{pipThreshold}, \code{filterCs} and
#'   \code{maxSnp}.
#' @param keepSnps Logical (length 1). When \code{TRUE}, retain the
#'   context-agnostic SNP background of each run as extra \code{CtwasResult}
#'   rows, one per region (\code{study = context = "SNP"}, mirroring cTWAS's
#'   own \code{"SNP"} group) so the full ctwas output is reconstructable from
#'   \code{\link{ctwasPosteriors}} / \code{susieAlpha}. Default \code{FALSE} --
#'   the SNP rows are the null background and are dropped from the structured
#'   gene-level result.
#' @param methodArgs Additional arguments forwarded to ctwas, built
#'   with \code{\link{CtwasOptions}}. Names are checked against what
#'   the ctwas steps accept between them.
#' @return A \code{\link{CtwasResult}} collection: one row per
#'   \code{(gwasStudy, study, context, method)} run and cTWAS region. A
#'   single-context run gives each method its own rows; a multi-context
#'   (joint) run emits per-context rows sharing the same \code{jointContexts}
#'   set and the jointly-estimated group priors, and a multi-study run tags
#'   \code{jointStudies} the same way. The contexts of a joint run
#'   need not hold the same genes: each context's prior comes from its own
#'   genes, and a context too small for cTWAS's \code{min_group_size} is left
#'   out of the fit with a warning. Each row holds that region's
#'   per-gene fine-mapping posteriors as ranges, with the region's per-effect
#'   susie alpha table (\code{susieAlpha}), its \code{regionInfo} and the
#'   run's \code{groupPriors} in \code{mcols()}. Its \code{seqinfo()} names
#'   the genome build the inputs name; inputs naming different builds are
#'   refused. For the raw \code{ctwas::finemap_regions} list (e.g. to
#'   feed \code{\link{mergeCtwasBoundaryRegions}}), call the granular
#'   \code{\link{assembleCtwasInputs}} \eqn{\to}
#'   \code{\link{estCtwasGroupPriors}} \eqn{\to}
#'   \code{\link{screenCtwasRegions}} \eqn{\to}
#'   \code{\link{finemapCtwasRegions}} path instead.
#' @examples
#' data(ctwasWeightsExample)
#' ldDir <- system.file("extdata", "ld_reference", "chr22",
#'   package = "pecotmr")
#' ldStem <- file.path(ldDir, "protocol_example.LD.chr22")
#' gwasTsv <- system.file("extdata", "manifests",
#'   "protocol_example.twas.gwas_sumstats.chr22.tsv.gz", package = "pecotmr")
#' mani <- data.frame(study = "gwas1", sumStatsPath = gwasTsv)
#' blocks <- GenomicRanges::GRanges("chr22",
#'   IRanges::IRanges(c(10000000, 15000001), c(15000000, 19000000)),
#'   blockId = c("chr22_1", "chr22_2"))
#' gss <- loadGwasSumStatsFromManifest(manifest = mani, genome = "hg38",
#'   ldReference = ldStem, which = blocks)
#' gwasByRegion <- summaryStatsQc(gss,
#'   panelFilterParam = PanelFilterParam(mafCutoff = 0.0025))
#' ctwasPipeline(gwasSumStats = gwasByRegion,
#'   twasWeights = list(ctwasWeightsExample),
#'   weightSelectionParam = WeightSelectionParam(methods = "susie"),
#'   ctwasPriorParam = CtwasPriorParam(thin = 1, niterPrefit = 3, niter = 10,
#'     fallbackToPrefit = TRUE),
#'   methodArgs = CtwasOptions(
#'     regionData = CtwasRegionDataOptions(min_group_size = 1),
#'     estParam = CtwasEstParamOptions(
#'       min_group_size = 1, min_p_single_effect = 0
#'     )
#'   ))
#' @export
ctwasPipeline <- function(
    gwasSumStats,
    twasWeights,
    twasZ = NULL,
    fineMappingResult = NULL,
    fineMappingMethod = NULL,
    weightSelectionParam = WeightSelectionParam(),
    maxNumSingleEffects = 5L,
    BPPARAM = SerialParam(),
    ctwasPriorParam = CtwasPriorParam(),
    variantPruningParam = VariantPruningParam(),
    boundaryMergeParam = BoundaryMergeParam(),
    keepSnps = FALSE,
    methodArgs = CtwasOptions()
) {
    .ctwasAssertBundles(
        ctwasPriorParam,
        variantPruningParam,
        boundaryMergeParam,
        methodArgs,
        weightSelectionParam
    )
    .assertBpparam(BPPARAM)
    .ctwasRequireNamedLists(gwasSumStats, twasWeights)
    plan <- .ctwasPlanRuns(
        gwasSumStats,
        twasWeights,
        twasZ,
        fineMappingResult,
        fineMappingMethod,
        weightSelectionParam
    )
    rows <- list_c(map(
        plan$runs,
        .ctwasRunMethod,
        plan = plan,
        variantPruningParam = variantPruningParam,
        ctwasPriorParam = ctwasPriorParam,
        BPPARAM = BPPARAM,
        maxNumSingleEffects = maxNumSingleEffects,
        boundaryMergeParam = boundaryMergeParam,
        keepSnps = keepSnps,
        methodArgs = methodArgs
    ))
    .ctwasRequireRows(rows, map_chr(plan$runs, "label"))
    .withGenome(.ctwasRowsToResult(rows), plan$build)
}

# What a pipeline's runs share: the inputs' genome build, the resolved inputs
# `inp`, the runs the weight choice makes (.ctwasRunPlan()), their first pass
# `fp`, and the supplied `twasZ` with the GWAS study it is read for.
# @noRd
.ctwasPlanRuns <- function(
    gwasSumStats,
    twasWeights,
    twasZ,
    fineMappingResult,
    fineMappingMethod,
    weightSelectionParam
) {
    build <- .inputsBuild(
        list(gwasSumStats, twasWeights, twasZ, fineMappingResult),
        "ctwasPipeline"
    )
    inp <- .ctwasResolveInputs(
        gwasSumStats,
        twasWeights,
        twasZ,
        fineMappingResult,
        fineMappingMethod
    )
    list(
        build = build,
        inp = inp,
        runs = .ctwasRunPlan(
            inp$twasWeights,
            weightSelectionParam,
            "ctwasPipeline()"
        ),
        fp = .ctwasFirstPassOf(inp),
        twasZ = twasZ,
        gwasStudy = .ctwasGwasStudy(gwasSumStats)
    )
}

# A run that modeled no gene has no result to return.
# @noRd
.ctwasRequireRows <- function(rows, labels) {
    if (length(rows) > 0L) {
        return(invisible(NULL))
    }
    abort(glue(
        "ctwasPipeline: no genes were modeled (the weight sources ",
        "produced no usable gene weights for method(s) ",
        "{str_flatten(labels, ', ')})."
    ))
}

# The cTWAS run proper: estimate the group priors, screen regions against
# them, fine-map what survived, and re-fine-map the merged boundary regions
# unless told not to. Each stage takes the ctwas option list as one named
# argument, so it is passed rather than spliced in as loose top-level args.
# @noRd
.ctwasRunStages <- function(
    inputs,
    ctwasPriorParam,
    boundaryMergeParam,
    maxNumSingleEffects,
    BPPARAM,
    methodArgs
) {
    est <- estCtwasGroupPriors(
        inputs,
        ctwasPriorParam = ctwasPriorParam,
        BPPARAM = BPPARAM,
        methodArgs = methodArgs
    )
    screened <- screenCtwasRegions(
        est,
        BPPARAM = BPPARAM,
        methodArgs = methodArgs
    )
    finemap <- finemapCtwasRegions(
        screened,
        maxNumSingleEffects = maxNumSingleEffects,
        BPPARAM = BPPARAM,
        methodArgs = methodArgs
    )
    if (!isTRUE(boundaryMergeParam$enabled)) {
        return(finemap)
    }
    .ctwasMaybeMerge(
        finemap,
        boundaryMergeParam = boundaryMergeParam,
        maxNumSingleEffects = maxNumSingleEffects,
        BPPARAM = BPPARAM,
        methodArgs = methodArgs
    )
}

# One cTWAS run of `run` over what the pipeline's runs share, `plan`
# (.ctwasPlanRuns()): assemble inputs -> estimate params -> screen ->
# fine-map -> boundary-merge -> per-context row-specs, labelled by the run.
# @noRd
.ctwasRunMethod <- function(
    run,
    plan,
    variantPruningParam,
    ctwasPriorParam,
    BPPARAM,
    maxNumSingleEffects,
    boundaryMergeParam,
    keepSnps,
    methodArgs
) {
    # The run plan has already warned about the rows whose fits this reads.
    inputs <- .muffleFitCoverage(.ctwasAssembleRun(
        plan$inp,
        plan$fp,
        run,
        variantPruningParam,
        plan$twasZ
    ))
    merged <- .ctwasRunStages(
        inputs,
        ctwasPriorParam = ctwasPriorParam,
        boundaryMergeParam = boundaryMergeParam,
        maxNumSingleEffects = maxNumSingleEffects,
        BPPARAM = BPPARAM,
        methodArgs = methodArgs
    )
    .ctwasRunToRows(
        merged,
        gwasStudy = plan$gwasStudy,
        method = run$label,
        keepSnps = keepSnps
    )
}

# Boundary-gene region merging: merge a high-PIP boundary gene's adjacent
# regions and re-fine-map. Merge-transparent downstream (keyed by gene id).
# @noRd
.ctwasMaybeMerge <- function(
    posteriors,
    boundaryMergeParam,
    maxNumSingleEffects,
    BPPARAM,
    methodArgs
) {
    # The bundle's terminal; mergeCtwasBoundaryRegions() takes scalars.
    mergePipThresh <- boundaryMergeParam$pipThreshold
    mergeFilterCs <- boundaryMergeParam$filterCs
    mergeMaxSNP <- boundaryMergeParam$maxSnp
    # methodArgs stays a record all the way down: mergeCtwasBoundaryRegions
    # takes one, and c()-ing it into an argument list would append the S4
    # object as a single element rather than splice its entries.
    mergeCtwasBoundaryRegions(
        posteriors,
        pipThresh = mergePipThresh,
        filterCs = mergeFilterCs,
        maxSNP = mergeMaxSNP,
        maxNumSingleEffects = maxNumSingleEffects,
        BPPARAM = BPPARAM,
        methodArgs = methodArgs
    )
}

#' Assemble cTWAS inputs from S4 GwasSumStats / TwasWeights
#'
#' @description Builds the per-block ctwas-shape input set (\code{z_snp},
#'   \code{weights}, \code{region_info}, \code{snp_map}, \code{LD_map}, the LD-
#'   and SNP-info loader closures, plus optional \code{z_gene}) that the
#'   downstream ctwas steps consume. This is step 1 of the three-step
#'   \code{\link{ctwasPipeline}} split.
#'
#' @details The returned list is the SHARED STATE threaded through
#'   \code{\link{estCtwasGroupPriors}} -> \code{\link{screenCtwasRegions}} ->
#'   \code{\link{finemapCtwasRegions}}. Callers can short-circuit at any step
#'   (e.g. override the estimated priors before fine-mapping) or call
#'   \code{ctwasPipeline()} for the one-shot path.
#'
#' @inheritParams ctwasPipeline
#' @param weightSelectionParam Which method's weights each gene is modelled
#'   with, built with \code{\link{WeightSelectionParam}}, as in
#'   \code{\link{ctwasPipeline}}. These are one cTWAS run's inputs, so the
#'   choice must make one run: \code{"best"}, or a single method;
#'   \code{ctwasPipeline} runs one per method for several.
#' @return A list with elements \code{z_snp}, \code{z_gene} (NULL when no
#'   \code{twasZ}), \code{weights}, \code{region_info}, \code{snp_map},
#'   \code{LD_map}, \code{LD_loader_fun}, \code{snpinfo_loader_fun},
#'   \code{resolvedMethod}, the run's label (\code{"best"}, or its one
#'   method), and \code{gwas_n}, the GWAS's median sample size (\code{NA}
#'   when it records none).
#' @examples
#' data(ctwasWeightsExample)
#' ldDir <- system.file("extdata", "ld_reference", "chr22",
#'   package = "pecotmr")
#' ldStem <- file.path(ldDir, "protocol_example.LD.chr22")
#' gwasTsv <- system.file("extdata", "manifests",
#'   "protocol_example.twas.gwas_sumstats.chr22.tsv.gz", package = "pecotmr")
#' mani <- data.frame(study = "gwas1", sumStatsPath = gwasTsv)
#' blocks <- GenomicRanges::GRanges("chr22",
#'   IRanges::IRanges(c(10000000, 15000001), c(15000000, 19000000)),
#'   blockId = c("chr22_1", "chr22_2"))
#' gss <- loadGwasSumStatsFromManifest(manifest = mani, genome = "hg38",
#'   ldReference = ldStem, which = blocks)
#' gwasByRegion <- summaryStatsQc(gss,
#'   panelFilterParam = PanelFilterParam(mafCutoff = 0.0025))
#' inputs <- assembleCtwasInputs(gwasSumStats = gwasByRegion,
#'   twasWeights = list(ctwasWeightsExample),
#'   weightSelectionParam = WeightSelectionParam(methods = "susie"))
#' names(inputs)
#' inputs$z_gene
#' @export
assembleCtwasInputs <- function(
    gwasSumStats,
    twasWeights,
    twasZ = NULL,
    fineMappingResult = NULL,
    fineMappingMethod = NULL,
    weightSelectionParam = WeightSelectionParam(),
    variantPruningParam = VariantPruningParam()
) {
    .assertMethodParam(
        variantPruningParam,
        "VariantPruningParam"
    )
    .assertMethodParam(weightSelectionParam, "WeightSelectionParam")
    inp <- .ctwasResolveInputs(
        gwasSumStats,
        twasWeights,
        twasZ,
        fineMappingResult,
        fineMappingMethod
    )
    runs <- .ctwasRunPlan(
        inp$twasWeights,
        weightSelectionParam,
        "assembleCtwasInputs()"
    )
    .ctwasAssertOneRun(runs)
    .ctwasAssembleRun(
        inp,
        .ctwasFirstPassOf(inp),
        runs[[1L]],
        variantPruningParam,
        twasZ
    )
}

# One cTWAS run's inputs from the resolved inputs `inp` and their first pass
# `fp`: the genes' weights cut to the run's chosen weight sets, labelled by
# the run (.ctwasRunPlan()).
# @noRd
.ctwasAssembleRun <- function(inp, fp, run, variantPruningParam, twasZ) {
    regionIds <- names(inp$gwasSumStats)
    weightsList <- .ctwasSecondPass(
        regionIds,
        inp$twasWeights,
        run$cells,
        fp,
        inp$fineMappingResult,
        .ctwasPruningCutoffs(variantPruningParam),
        unique(list_c(map(fp$zSnpPieces, "id")))
    )
    .ctwasAssembleResult(
        regionIds,
        fp,
        weightsList,
        list(twasZ = twasZ, gwasStudy = inp$gwasStudy),
        run$label
    )
}

# The first pass over the resolved inputs `inp`, which no weight choice
# changes, so a pipeline's runs share it.
# @noRd
.ctwasFirstPassOf <- function(inp) {
    .ctwasFirstPass(
        names(inp$gwasSumStats),
        inp$gwasSumStats,
        inp$twasWeights,
        inp$globalPanelInfo
    )
}

# An assembled input set is one cTWAS run.
# @noRd
.ctwasAssertOneRun <- function(runs) {
    if (length(runs) == 1L) {
        return(invisible(NULL))
    }
    abort(c(
        glue(
            "assembleCtwasInputs: the weight choice makes {length(runs)} ",
            "cTWAS runs, one per method ",
            "({str_flatten(map_chr(runs, 'label'), ', ')}), and these are one ",
            "run's inputs."
        ),
        i = str_c(
            "Name one method, use methods = \"best\", or let ctwasPipeline() ",
            "run one per method."
        )
    ))
}

# The inputs of assembleCtwasInputs(), checked and resolved in the order a bad
# one should fail: the GWAS split by block with its panel summary, the weights
# placed on those blocks, and the fine-mapping source.
# @noRd
.ctwasResolveInputs <- function(
    gwasSumStats,
    twasWeights,
    twasZ,
    fineMappingResult,
    fineMappingMethod
) {
    .ctwasValidateGwasList(gwasSumStats)
    globalPanelInfo <- .ctwasGlobalPanelInfo(gwasSumStats)
    # Only a supplied twasZ needs the run's one GWAS study, to pick its rows.
    gwasStudy <- if (is.null(twasZ)) {
        NA_character_
    } else {
        .ctwasGwasStudy(gwasSumStats)
    }
    byBlock <- .ctwasGwasByBlock(gwasSumStats)
    placed <- .ctwasResolveAndValidateWeights(twasWeights, byBlock)
    .ctwasValidateOptional(twasZ, fineMappingResult)
    list(
        gwasSumStats = byBlock,
        twasWeights = placed,
        globalPanelInfo = globalPanelInfo,
        gwasStudy = gwasStudy,
        fineMappingResult = .ctwasRescueSource(
            fineMappingResult,
            fineMappingMethod
        )
    )
}

# The bundle's terminal: .ctwasGeneContext() and below read scalars, under the
# names they have always had there.
# @noRd
.ctwasPruningCutoffs <- function(variantPruningParam) {
    list(
        twasWeightCutoff = variantPruningParam$weightCutoff,
        csMinCor = variantPruningParam$csMinCor,
        minPipCutoff = variantPruningParam$pipRescue,
        maxNumVariants = variantPruningParam$maxVariants
    )
}

# Validate the gwasSumStats input: ctwas available, a QC'd GwasSumStats whose
# elements are per-LD-block (>= 2 of them).
#
# This used to demand a named list of GwasSumStats keyed by region_id, which
# existed only because the class could not hold multiple region-scoped blocks.
# It can now (§4.2): one collection, one element per block, keyed by `blockId`.
# @noRd
.ctwasValidateGwasList <- function(gwasSumStats) {
    if (!requireNamespace("ctwas", quietly = TRUE)) {
        abort("Package 'ctwas' is required for the cTWAS pipeline.")
    }
    if (missing(gwasSumStats) || !methods::is(gwasSumStats, "GwasSumStats")) {
        msg <- glue(
            "`gwasSumStats` must be a GwasSumStats whose elements are LD ",
            "blocks (got {class(gwasSumStats)[[1L]]}). Build one with ",
            "`loadGwasSumStatsFromManifest(..., region = <blocks>)`."
        )
        abort(msg)
    }
    if (!is_in("blockId", colnames(gwasSumStats))) {
        msg <- glue(
            "`gwasSumStats` has no `blockId` column, so its elements cannot ",
            "be keyed by LD block. Rebuild it with a current constructor."
        )
        abort(msg)
    }
    if (length(gwasSumStats) < 2L) {
        msg <- glue(
            "assembleCtwasInputs: at least two LD blocks are required ",
            "(got {length(gwasSumStats)}). cTWAS's EM cannot estimate the ",
            "SNP-group prior variance from a single region."
        )
        abort(msg)
    }
    .ctwasValidateGwasEntries(gwasSumStats)
}

# The collection must be QC'd, and its block keys must be usable as region ids.
# @noRd
.ctwasValidateGwasEntries <- function(gwasSumStats) {
    if (length(qcOptions(gwasSumStats)) == 0L) {
        msg <- glue(
            "assembleCtwasInputs: `gwasSumStats` has no QC record. ",
            "Call summaryStatsQc() first."
        )
        abort(msg)
    }
    keys <- as.character(gwasSumStats$blockId)
    if (any(is.na(keys)) || any(str_length(keys) == 0L)) {
        abort("`gwasSumStats` has empty or missing `blockId` value(s).")
    }
    if (anyDuplicated(keys) > 0L) {
        dup <- unique(keys[duplicated(keys)])
        msg <- glue(
            "`gwasSumStats` block ids must be unique; repeated: ",
            "{str_flatten(head(dup, 5L), ', ')}. Two elements sharing a ",
            "block id would silently overwrite one another as region ids."
        )
        abort(msg)
    }
}

# One single-element GwasSumStats per LD block, keyed by blockId -- the shape
# the two-pass assembly consumes. Each block keeps its own row's LD reference,
# so per-region LD lookups read the panel that block was harmonized against.
# @noRd
.ctwasGwasByBlock <- function(gwasSumStats) {
    set_names(
        map(seq_along(gwasSumStats), .ctwasPickBlock, x = gwasSumStats),
        as.character(gwasSumStats$blockId)
    )
}

# One block's slice of the collection. A row's LD reference can cover far more
# than its block -- combineGwasSumStats() unions the per-piece panels of one
# source -- so without narrowing here a block would reference the whole
# chromosome and its LD panel would cover all of it (a chr19-wide panel runs
# to ~155 GB). Narrowed to the block's own variants the same way
# summaryStatsQc narrows the sketch it retains.
# @noRd
.ctwasPickBlock <- function(i, x) {
    block <- x[i]
    refs <- ldReference(block)
    if (length(refs) == 0L) {
        return(block)
    }
    .withRowLdRefs(block, .subsetLdRefToIds(refs[[1L]], as.list(block)))
}

# Resolve a flat weight source into per-region buckets (.ctwasBucketWeights(); a
# pre-bucketed named list passes through), then validate it is a named list of
# TwasWeights / QtlFineMappingResult with no stray region keys.
# Returns the resolved twasWeights.
# @noRd
.ctwasResolveAndValidateWeights <- function(twasWeights, gwasSumStats) {
    if (missing(twasWeights) || is.null(twasWeights)) {
        msg <- glue(
            "`twasWeights` is required (a TwasWeights / QtlFineMappingResult ",
            "weight source, or a per-region named list keyed by region_id)."
        )
        abort(msg)
    }
    twasWeights <- .ctwasResolveWeightBuckets(twasWeights, gwasSumStats)
    if (
        is.null(names(twasWeights)) ||
            any(str_length(names(twasWeights)) == 0L)
    ) {
        msg <- glue(
            "`twasWeights` must resolve to a named list keyed by region_id ",
            "(got an unnamed or empty-named list)."
        )
        abort(msg)
    }
    extraKeys <- setdiff(names(twasWeights), names(gwasSumStats))
    if (length(extraKeys) > 0L) {
        msg <- glue(
            "`twasWeights` has region_id key(s) not present in ",
            "`gwasSumStats`: {str_flatten(extraKeys, ', ')}"
        )
        abort(msg)
    }
    .ctwasValidateWeightEntries(twasWeights)
    twasWeights
}

# Each twasWeights entry must be a TwasWeights or QtlFineMappingResult.
# @noRd
.ctwasValidateWeightEntries <- function(twasWeights) {
    for (rid in names(twasWeights)) {
        if (
            !methods::is(twasWeights[[rid]], "TwasWeights") &&
                !methods::is(twasWeights[[rid]], "QtlFineMappingResult")
        ) {
            msg <- glue(
                "twasWeights[['{rid}']] must be a TwasWeights or ",
                "QtlFineMappingResult (the per-gene weight source)."
            )
            abort(msg)
        }
    }
}

# Optional twasZ (GRanges) and fineMappingResult (FineMappingResultBase) types.
# @noRd
.ctwasValidateOptional <- function(twasZ, fineMappingResult) {
    if (!is.null(twasZ) && !methods::is(twasZ, "GRanges")) {
        msg <- glue(
            "`twasZ` must be a GRanges (output of causalInferencePipeline) ",
            "or NULL."
        )
        abort(msg)
    }
    # The rescues keep a gene's weight variants by that gene's own fine-mapping,
    # which is a QTL fit: the weights come from a TwasWeights or a
    # QtlFineMappingResult, and a GWAS result is keyed by a study and method
    # that no gene shares.
    if (
        !is.null(fineMappingResult) &&
            !methods::is(fineMappingResult, "QtlFineMappingResult")
    ) {
        msg <- glue(
            "`fineMappingResult` must be a QtlFineMappingResult -- the genes' ",
            "own fine-mapping -- or NULL."
        )
        abort(msg)
    }
}

# The fine-mapping the rescues read: `fineMappingResult` narrowed to the one
# method `fineMappingMethod` names. The weights' own method cannot pick it --
# an "ensemble", "lasso" or "mrash" weight has no fine-mapping fit of that
# name -- so the method is named separately, and is required with a result.
# @noRd
.ctwasRescueSource <- function(fineMappingResult, fineMappingMethod) {
    assertString(fineMappingMethod, null.ok = TRUE)
    if (is.null(fineMappingResult)) {
        if (!is.null(fineMappingMethod)) {
            abort(glue(
                "`fineMappingMethod` names a method of `fineMappingResult`, ",
                "which was not given."
            ))
        }
        return(NULL)
    }
    available <- unique(as.character(fineMappingResult$method))
    if (is.null(fineMappingMethod)) {
        abort(glue(
            "`fineMappingResult` is given, so `fineMappingMethod` must name ",
            "the method of it whose PIPs and credible sets the rescues read ",
            "(available: {str_flatten(available, ', ')})."
        ))
    }
    if (!is_in(fineMappingMethod, available)) {
        abort(glue(
            "`fineMappingMethod` '{fineMappingMethod}' is not a method of ",
            "`fineMappingResult` (available: {str_flatten(available, ', ')})."
        ))
    }
    fineMappingResult[fineMappingResult$method == fineMappingMethod]
}

# The ctwas-shaped snpInfo for the whole collection's panel, captured BEFORE
# .ctwasGwasByBlock() narrows each block's sketch to its own variants. A
# gene's cis SPAN has to cover every block it reaches or
# ctwas::get_boundary_genes cannot route it to merge_regions, so the span is
# measured against this table rather than against any one block's panel.
#
# NULL when the collection carries no sketch -- .ctwasRegionGwasLd() reports
# that case with its own message, so this must not pre-empt it. Rows on
# different references contribute each panel's variants, each once.
# @noRd
.ctwasGlobalPanelInfo <- function(gwasSumStats) {
    refs <- ldReference(gwasSumStats)
    if (length(refs) == 0L) {
        return(NULL)
    }
    info <- list_rbind(map(unname(refs), .ctwasSnpInfoForBlock))
    filter(info, !duplicated(.data$id))
}

# One region's GWAS LD sketch, checked against the region's weights. A
# region's LD is one matrix, so its GWAS rows must share one reference, and
# every weight row in the region is combined with them.
# @noRd
.ctwasRegionGwasLd <- function(rid, gwasSumStats, twasWeights) {
    gwas <- gwasSumStats[[rid]]
    refs <- ldReference(gwas)
    if (length(refs) == 0L) {
        msg <- glue(
            "ctwasPipeline: GwasSumStats for region '{rid}' carries no ",
            "ldReference (ldReference = NULL); cTWAS requires an LD reference."
        )
        abort(msg)
    }
    if (length(refs) > 1L) {
        msg <- glue(
            "ctwasPipeline: GwasSumStats for region '{rid}' has rows on ",
            "{length(refs)} different LD references; a region's LD is one ",
            "matrix."
        )
        abort(msg)
    }
    tw <- twasWeights[[rid]]
    if (!is.null(tw)) {
        .ctwasRequireMatchingLdRefs(tw, gwas)
    }
    refs[[1L]]
}

# @noRd
.ctwasZSnpAt <- function(rid, ldPanel, gwasSumStats) {
    .ctwasBuildZSnp(gwasSumStats[[rid]], ldPanel$snpInfo$id)
}

# @noRd
.ctwasSnpMapAt <- function(rid, ldPanel, gwasSumStats) {
    .ctwasSnpInfoForGwasBlock(gwasSumStats[[rid]], ldPanel$snpInfo)
}

# @noRd
.ctwasRegionInfoAt <- function(rid, gwasSumStats) {
    .ctwasBuildSingleRegionInfo(rid, gwasSumStats[[rid]])
}

# First pass: one LD panel per region, plus the z_snp / region_info / snp_map
# pieces built from them.
#
# Each region gets a token of its own. The panel cache used to be keyed on the
# genotype file backing the sketch, which every block on a chromosome shares,
# so only the first block's panel was ever computed and every other region
# silently read that one back out.
#
# `globalVariance` is assembled from the per-block variance vectors rather
# than re-read: a boundary gene needs the variance of its out-of-block weight
# variants too, and the blocks partition the panel, so their union is exactly
# the global vector. It costs nothing beyond the panels themselves.
# @noRd
.ctwasFirstPass <- function(
    regionIds,
    gwasSumStats,
    twasWeights,
    globalPanelInfo = NULL
) {
    sketches <- map(
        regionIds,
        .ctwasRegionGwasLd,
        gwasSumStats = gwasSumStats,
        twasWeights = twasWeights
    )
    ldTokens <- .ctwasRegionLdTokens(regionIds, sketches)
    ldPanelsByRegion <- set_names(
        map(sketches, .ctwasPanelFor),
        unname(ldTokens)
    )
    panels <- unname(ldPanelsByRegion)
    list(
        globalPanelInfo = globalPanelInfo,
        globalVariance = list_c(map(panels, "variance")),
        ldPanelsByRegion = ldPanelsByRegion,
        zSnpPieces = set_names(
            map2(regionIds, panels, .ctwasZSnpAt, gwasSumStats = gwasSumStats),
            regionIds
        ),
        regionInfoPieces = set_names(
            map(regionIds, .ctwasRegionInfoAt, gwasSumStats = gwasSumStats),
            regionIds
        ),
        snpMap = set_names(
            map2(
                regionIds,
                panels,
                .ctwasSnpMapAt,
                gwasSumStats = gwasSumStats
            ),
            regionIds
        ),
        ldFileByRegion = ldTokens
    )
}

# Second pass: the genes' weight entries, keyed by gene id, each gene's being
# the weight set of it `chosen` names (.ctwasRunPlan()). Each gene is prepared
# from the block its weights were placed in, whichever blocks they reach, and
# the variances then read each block's LD once, over all the genes.
# @noRd
.ctwasSecondPass <- function(
    regionIds,
    twasWeights,
    chosen,
    fp,
    fineMappingResult,
    cutoffs,
    globalGwasSnpIds
) {
    ctx <- .ctwasGeneContext(
        fp$globalPanelInfo,
        fp$globalVariance,
        globalGwasSnpIds,
        fineMappingResult,
        cutoffs
    )
    perRegion <- map(
        regionIds,
        .ctwasRegionGenes,
        twasWeights = twasWeights,
        chosen = chosen,
        ctx = ctx
    )
    .ctwasGeneEntries(
        .ctwasConcat(compact(perRegion)),
        unname(fp$ldPanelsByRegion)
    )
}

# Concatenate per-region lists, empty-safe.
# @noRd
.ctwasConcat <- function(pieces) {
    if (length(pieces) == 0L) {
        return(list())
    }
    list_c(pieces)
}

# The genes placed in one region, prepared up to their weight variances, or NULL
# when the region holds none of the chosen weight sets.
# @noRd
.ctwasRegionGenes <- function(rid, twasWeights, chosen, ctx) {
    tw <- twasWeights[[rid]]
    if (is.null(tw)) {
        return(NULL)
    }
    twChosen <- .ctwasKeepChosen(tw, chosen)
    if (is.null(twChosen)) {
        return(NULL)
    }
    .ctwasPreparedGenes(twChosen, ctx)
}

# Concatenate the per-region pieces into the ctwas-shape input list.
# @noRd
.ctwasAssembleResult <- function(
    regionIds,
    fp,
    weightsList,
    supplied,
    resolvedMethod
) {
    # The ctwas engine indexes its input frames POSITIONALLY (`df[, "col"]` must
    # yield a vector); a tibble returns a 1-col list and breaks min()/downstream
    # inside the ctwas run (verified). So every ctwas-input frame is coerced
    # to a
    # base data.frame HERE, at the external-ctwas boundary -- the builders stay
    # tidyverse-internal (tibbles), the coercion is confined to this assembly.
    zSnp <- as.data.frame(bind_rows(fp$zSnpPieces))
    regionInfo <- as.data.frame(bind_rows(fp$regionInfoPieces))
    ldMap <- data.frame(
        region_id = regionIds,
        LD_file = unname(fp$ldFileByRegion),
        SNP_file = unname(fp$ldFileByRegion),
        stringsAsFactors = FALSE
    )
    list(
        z_snp = zSnp,
        z_gene = if (!is.null(supplied$twasZ)) {
            .ctwasBuildZGene(supplied$twasZ, weightsList, supplied$gwasStudy)
        } else {
            .ctwasGeneZ(weightsList, zSnp)
        },
        weights = weightsList,
        region_info = regionInfo,
        snp_map = fp$snpMap,
        LD_map = ldMap,
        LD_loader_fun = .ctwasMultiBlockLdLoader(fp$ldPanelsByRegion),
        snpinfo_loader_fun = .ctwasMultiBlockSnpInfoLoader(fp$ldPanelsByRegion),
        resolvedMethod = resolvedMethod,
        gwas_n = .ctwasGwasN(zSnp)
    )
}

# The GWAS's sample size, one number as cTWAS's summaries read it: the median
# over its variants, NA when it records none.
# @noRd
.ctwasGwasN <- function(zSnp) {
    if (is.null(zSnp$N) || all(is.na(zSnp$N))) {
        return(NA_real_)
    }
    median(zSnp$N, na.rm = TRUE)
}

# estCtwasGroupPriors(): CtwasPriorParam is the bundle's terminal; the EM
# helpers take its fields as plain scalars. assemble_region_data does not echo
# z_gene back, so the precomputed z_gene passed in is propagated onto the
# result (inputs$z_gene is NULL when twasZ was not supplied) so $z_gene
# resolves to the right entry.
#' Estimate cTWAS group prior + prior variance
#'
#' @description Step 2 of the three-step \code{\link{ctwasPipeline}}: assembles
#'   \code{region_data} from the inputs and runs \code{ctwas::est_param} (prefit
#'   EM + accurate EM) to estimate the group prior probabilities and prior
#'   variances. Returns the input state plus \code{region_data},
#'   \code{boundary_genes}, \code{z_gene}, and \code{param}. A context with
#'   fewer genes than \code{min_group_size} gets no prior: cTWAS leaves it out
#'   of the fit, and a warning names it. To keep it, lower
#'   \code{min_group_size} in both \code{\link{CtwasRegionDataOptions}} and
#'   \code{\link{CtwasEstParamOptions}}.
#'
#' @param inputs A list returned by \code{\link{assembleCtwasInputs}}.
#' @param ctwasPriorParam The EM settings, built with
#'   \code{\link{CtwasPriorParam}}: \code{thin} (SNP thinning when
#'   assembling region data), \code{niterPrefit} and \code{niter} (prefit /
#'   accurate EM iterations), \code{varStructure}, and
#'   \code{fallbackToPrefit} -- when \code{TRUE}, an accurate EM that fails
#'   on a degenerate input is recovered by re-running only the prefit step
#'   and returning those (typically finite) priors.
#' @param BPPARAM A \code{BiocParallelParam} whose workers number the cores
#'   ctwas's own \code{ncore} parallelism uses when they run on this machine
#'   (\code{SerialParam}, \code{MulticoreParam}); any other backend gives one
#'   core. Default \code{SerialParam()}.
#' @param methodArgs Additional arguments forwarded to ctwas, built
#'   with \code{\link{CtwasOptions}}. Names are checked against what
#'   the ctwas steps accept between them.
#' @return The \code{inputs} list augmented with \code{region_data},
#'   \code{boundary_genes}, \code{z_gene}, and \code{param}.
#' @examples
#' data(ctwasInputsExample)
#' est <- estCtwasGroupPriors(ctwasInputsExample,
#'   ctwasPriorParam = CtwasPriorParam(thin = 1, niterPrefit = 3, niter = 10,
#'     fallbackToPrefit = TRUE),
#'   methodArgs = CtwasOptions(
#'     regionData = CtwasRegionDataOptions(min_group_size = 1),
#'     estParam = CtwasEstParamOptions(
#'       min_group_size = 1, min_p_single_effect = 0
#'     )
#'   ))
#' est$param[c("group_prior", "group_prior_var")]
#' @export
estCtwasGroupPriors <- function(
    inputs,
    ctwasPriorParam = CtwasPriorParam(),
    BPPARAM = SerialParam(),
    methodArgs = CtwasOptions()
) {
    .assertMethodParam(ctwasPriorParam, "CtwasPriorParam")
    .assertMethodOptions(methodArgs, "CtwasOptions", "methodArgs")
    .assertBpparam(BPPARAM)
    if (!requireNamespace("ctwas", quietly = TRUE)) {
        abort("Package 'ctwas' is required for estCtwasGroupPriors.")
    }
    threads <- .bpLocalThreads(BPPARAM, "estCtwasGroupPriors")
    inputs <- .ctwasResolveLdPaths(inputs)
    zGene <- .ctwasEnsureZGene(inputs, threads, methodArgs$geneZ)
    regionData <- .ctwasAssembleRegionData(
        inputs,
        zGene,
        ctwasPriorParam$thin,
        threads,
        methodArgs$regionData
    )
    .ctwasWarnDroppedGroups(zGene, regionData)
    boundaryGenes <- .ctwasBoundaryGenes(
        inputs,
        threads,
        methodArgs$boundaryGenes
    )
    paramRes <- .ctwasEstParamOrFallback(
        regionData,
        ctwasPriorParam$niterPrefit,
        ctwasPriorParam$niter,
        ctwasPriorParam$varStructure,
        threads,
        ctwasPriorParam$thin,
        ctwasPriorParam$fallbackToPrefit,
        methodArgs$estParam
    )
    c(
        list_assign(inputs, z_gene = zGene),
        list(
            region_data = regionData,
            boundary_genes = boundaryGenes,
            param = paramRes
        )
    )
}

# cTWAS leaves a group with too few genes (`min_group_size`) out of the fit and
# says so only in a log pecotmr quiets, so its context would come back with no
# priors and no posteriors, and nothing to say why. Contexts holding different
# genes make this reachable: one can fall short while another clears it. The
# advice names both options because `est_param` holds its own limit, and stops
# on a group under it where `assemble_region_data` drops one.
# @noRd
.ctwasWarnDroppedGroups <- function(zGene, regionData) {
    dropped <- !is_in(zGene$group, list_c(map(regionData, "groups")))
    if (!any(dropped)) {
        return(invisible())
    }
    counts <- lengths(split(zGene$id[dropped], zGene$context[dropped]))
    who <- str_c(names(counts), " (", counts, " gene(s))")
    warn(glue(
        "estCtwasGroupPriors: cTWAS left context(s) ",
        "{str_flatten(who, ', ')} out of the fit, so they have no priors or ",
        "posteriors: it drops a group with fewer genes than `min_group_size`. ",
        "To keep them, lower `min_group_size` in both ",
        "CtwasRegionDataOptions() and CtwasEstParamOptions(): the second ",
        "stops on a group under its own limit."
    ))
}

# est_param (accurate EM), falling back to a converged prefit EM on ANY accurate
# error when fallbackToPrefit is set. The accurate EM fails on degenerate inputs
# in several version-dependent ways (NAs / "No regions selected!" / NaN
# log-likelihood), so catch all rather than match brittle version messages. The
# prefit re-run runs to full `niter` (this is now the final prior).
# @noRd
#' @importFrom rlang try_fetch
.ctwasEstParamOrFallback <- function(
    regionData,
    niterPrefit,
    niter,
    groupPriorVarStructure,
    threads,
    thin,
    fallbackToPrefit,
    extra
) {
    # No fallback requested: run the accurate EM directly and let any error
    # propagate (no catch-then-rethrow).
    if (!fallbackToPrefit) {
        return(.ctwasEstParamAccurate(
            regionData,
            niterPrefit,
            niter,
            groupPriorVarStructure,
            threads,
            extra
        ))
    }
    try_fetch(
        .ctwasEstParamAccurate(
            regionData,
            niterPrefit,
            niter,
            groupPriorVarStructure,
            threads,
            extra
        ),
        error = function(cnd) {
            msg <- glue(
                "estCtwasGroupPriors: accurate EM unusable; falling back to ",
                "prefit estimates."
            )
            inform(msg, parent = cnd)
            .ctwasFitPrefitEm(
                regionData,
                niter = as.integer(niter),
                groupPriorVarStructure = groupPriorVarStructure,
                thin = thin,
                threads = threads,
                extra = extra
            )
        }
    )
}

#' Fine-map cTWAS regions
#'
#' @description Step 4 (final) of the three-step \code{\link{ctwasPipeline}}:
#'   runs \code{ctwas::finemap_regions} on the screened-region set from
#'   \code{\link{screenCtwasRegions}} and assembles the documented top-level
#'   ctwas output (\code{z_gene}, \code{param}, \code{finemap_res},
#'   \code{susie_alpha_res}, \code{region_data}, \code{boundary_genes},
#'   \code{screen_res}).
#'
#' @param screenResult A list returned by \code{\link{screenCtwasRegions}}.
#' @param maxNumSingleEffects Pass-through.
#' @param BPPARAM A \code{BiocParallelParam} whose workers number the cores
#'   ctwas's own \code{ncore} parallelism uses when they run on this machine
#'   (\code{SerialParam}, \code{MulticoreParam}); any other backend gives one
#'   core. Default \code{SerialParam()}.
#' @param methodArgs Additional arguments forwarded to ctwas, built
#'   with \code{\link{CtwasOptions}}. Names are checked against what
#'   the ctwas steps accept between them.
#' @return A list mirroring \code{ctwas::ctwas_sumstats}'s output:
#'   \code{z_gene}, \code{param}, \code{finemap_res}, \code{susie_alpha_res},
#'   \code{region_data}, \code{boundary_genes}, \code{screen_res}.
#' @examples
#' data(ctwasWeightsExample)
#' ldDir <- system.file("extdata", "ld_reference", "chr22",
#'   package = "pecotmr")
#' ldStem <- file.path(ldDir, "protocol_example.LD.chr22")
#' gwasTsv <- system.file("extdata", "manifests",
#'   "protocol_example.twas.gwas_sumstats.chr22.tsv.gz", package = "pecotmr")
#' mani <- data.frame(study = "gwas1", sumStatsPath = gwasTsv)
#' blocks <- GenomicRanges::GRanges("chr22",
#'   IRanges::IRanges(c(10000000, 15000001), c(15000000, 19000000)),
#'   blockId = c("chr22_1", "chr22_2"))
#' gss <- loadGwasSumStatsFromManifest(manifest = mani, genome = "hg38",
#'   ldReference = ldStem, which = blocks)
#' gwasByRegion <- summaryStatsQc(gss,
#'   panelFilterParam = PanelFilterParam(mafCutoff = 0.0025))
#' inp <- assembleCtwasInputs(gwasSumStats = gwasByRegion,
#'   twasWeights = list(ctwasWeightsExample),
#'   weightSelectionParam = WeightSelectionParam(methods = "susie"))
#' est <- estCtwasGroupPriors(inp,
#'   ctwasPriorParam = CtwasPriorParam(thin = 1, niterPrefit = 3, niter = 10,
#'     fallbackToPrefit = TRUE),
#'   methodArgs = CtwasOptions(
#'     regionData = CtwasRegionDataOptions(min_group_size = 1),
#'     estParam = CtwasEstParamOptions(
#'       min_group_size = 1, min_p_single_effect = 0
#'     )
#'   ))
#' screened <- screenCtwasRegions(est)
#' finemapped <- finemapCtwasRegions(screened, maxNumSingleEffects = 5L)
#' head(finemapped$finemap_res)
#' @export
finemapCtwasRegions <- function(
    screenResult,
    maxNumSingleEffects = 5L,
    BPPARAM = SerialParam(),
    methodArgs = CtwasOptions()
) {
    .assertMethodOptions(methodArgs, "CtwasOptions", "methodArgs")
    .assertBpparam(BPPARAM)
    if (!requireNamespace("ctwas", quietly = TRUE)) {
        abort("Package 'ctwas' is required for finemapCtwasRegions.")
    }
    threads <- .bpLocalThreads(BPPARAM, "finemapCtwasRegions")
    screenResult <- .ctwasResolveLdPaths(screenResult)
    fmRes <- .ctwasFinemapOrEmpty(
        screenResult,
        maxNumSingleEffects,
        threads,
        methodArgs
    )
    list(
        z_gene = screenResult$z_gene,
        param = screenResult$param,
        gwas_n = screenResult$gwas_n,
        finemap_res = .ctwasFixMolecularId(
            fmRes$finemap_res,
            screenResult$weights
        ),
        susie_alpha_res = .ctwasFixMolecularId(
            fmRes$susie_alpha_res,
            screenResult$weights
        ),
        region_data = screenResult$region_data,
        boundary_genes = screenResult$boundary_genes,
        screen_res = screenResult$screen_res,
        # Carried forward so mergeCtwasBoundaryRegions() can re-finemap the
        # merged boundary regions without re-deriving the assembled inputs.
        region_info = screenResult$region_info,
        z_snp = screenResult$z_snp,
        weights = screenResult$weights,
        snp_map = screenResult$snp_map,
        LD_map = screenResult$LD_map,
        LD_loader_fun = screenResult$LD_loader_fun,
        snpinfo_loader_fun = screenResult$snpinfo_loader_fun,
        # The run's label (assembleCtwasInputs()), which asCtwasResult()
        # labels the result with: a "best" run's genes mix methods.
        resolvedMethod = screenResult$resolvedMethod
    )
}

#' Merge boundary cTWAS regions and re-fine-map
#'
#' @description Step 4 of the cTWAS pipeline, on by default in
#'   \code{\link{ctwasPipeline}}. A gene whose weights cross an LD-block
#'   boundary (a \code{boundary_genes} member) is fine-mapped in one block in
#'   the first pass, with only that block's LD. This step selects the high-PIP
#'   boundary genes, merges each one's adjacent regions into a single region,
#'   re-runs fine-mapping on the merged regions, and splices the updated results
#'   back into the \code{\link{finemapCtwasRegions}} output. Thin wrapper over
#'   \code{ctwas::postprocess_region_merging()} (or
#'   \code{ctwas::postprocess_region_merging_noLD()} when the inputs carry no LD
#'   loaders).
#'
#' @param finemapResult A list returned by \code{\link{finemapCtwasRegions}}.
#'   Must carry \code{finemap_res}, \code{susie_alpha_res}, \code{region_data},
#'   \code{region_info}, \code{z_snp}, \code{z_gene}, \code{weights},
#'   \code{snp_map}, \code{param}, and -- on the LD path -- \code{LD_map} plus
#'   the \code{LD_loader_fun} / \code{snpinfo_loader_fun} closures (all retained
#'   by \code{finemapCtwasRegions}).
#' @param pipThresh Numeric (length 1). PIP threshold for selecting which
#'   boundary genes to merge (\code{select_boundary_genes} \code{pip_thresh}).
#'   Default \code{0.5}.
#' @param filterCs Logical (length 1). Require the gene to be in a credible set
#'   to be selected (\code{select_boundary_genes} \code{filter_cs}). Default
#'   \code{FALSE}.
#' @param maxSNP Numeric (length 1). Per-merged-region SNP cap. Default
#'   \code{Inf}.
#' @param maxNumSingleEffects Integer. Max number of single effects for the
#'   merged-region re-fine-mapping (LD path only). Default \code{5}.
#' @param BPPARAM A \code{BiocParallelParam} whose workers number the cores
#'   ctwas's own \code{ncore} parallelism uses when they run on this machine
#'   (\code{SerialParam}, \code{MulticoreParam}); any other backend gives one
#'   core. Default \code{SerialParam()}.
#' @param methodArgs Additional arguments forwarded to ctwas, built
#'   with \code{\link{CtwasOptions}}. Names are checked against what
#'   the ctwas steps accept between them.
#' @return The \code{finemapResult} list with \code{finemap_res},
#'   \code{susie_alpha_res}, \code{region_data}, \code{region_info},
#'   \code{LD_map}, and \code{snp_map} replaced by the post-merge ("updated")
#'   values, plus a \code{merge_res} element carrying the full ctwas postprocess
#'   output. When no boundary gene clears \code{pipThresh}, ctwas returns the
#'   inputs as the "updated" values, so the result is effectively unchanged.
#' @examples
#' data(ctwasFinemapExample)
#' merged <- mergeCtwasBoundaryRegions(ctwasFinemapExample)
#' head(merged$finemap_res)
#' @export
mergeCtwasBoundaryRegions <- function(
    finemapResult,
    pipThresh = 0.5,
    filterCs = FALSE,
    maxSNP = Inf,
    maxNumSingleEffects = 5L,
    BPPARAM = SerialParam(),
    methodArgs = CtwasOptions()
) {
    .assertMethodOptions(methodArgs, "CtwasOptions", "methodArgs")
    .assertBpparam(BPPARAM)
    if (!requireNamespace("ctwas", quietly = TRUE)) {
        abort("Package 'ctwas' is required for mergeCtwasBoundaryRegions.")
    }
    threads <- .bpLocalThreads(BPPARAM, "mergeCtwasBoundaryRegions")
    finemapResult <- .ctwasResolveLdPaths(finemapResult)
    fmRes <- finemapResult$finemap_res
    if (is.null(fmRes) || nrow(fmRes) == 0L) {
        msg <- glue(
            "mergeCtwasBoundaryRegions: no first-pass finemap result; ",
            "returning unchanged."
        )
        inform(msg)
        return(finemapResult)
    }
    common <- .ctwasMergeCommonArgs(
        finemapResult,
        pipThresh,
        filterCs,
        maxSNP,
        threads
    )
    fa <- .ctwasMergeDispatch(finemapResult, common, maxNumSingleEffects)
    # Flatten the record before merging: c() on a list and a MethodOptions (a
    # SimpleList) appends the S4 object as ONE unnamed element rather than
    # splicing its entries, and exec() would then hand that object to ctwas
    # as a positional argument.
    user <- as.list(methodArgs$merge)
    callArgs <- c(fa$args, user[setdiff(names(user), names(fa$args))])
    res <- .ctwasQuietly(fa$fn, callArgs)
    .ctwasApplyMergeResult(finemapResult, res)
}

# Shared postprocess_region_merging argument list built from a finemap result.
# @noRd
.ctwasMergeCommonArgs <- function(
    finemapResult,
    pipThresh,
    filterCs,
    maxSNP,
    threads
) {
    list(
        region_info = finemapResult$region_info,
        region_data = finemapResult$region_data,
        z_snp = finemapResult$z_snp,
        z_gene = finemapResult$z_gene,
        weights = finemapResult$weights,
        snp_map = finemapResult$snp_map,
        finemap_res = finemapResult$finemap_res,
        susie_alpha_res = finemapResult$susie_alpha_res,
        group_prior = finemapResult$param$group_prior,
        group_prior_var = finemapResult$param$group_prior_var,
        pip_thresh = pipThresh,
        filter_cs = filterCs,
        maxSNP = maxSNP,
        ncore = threads
    )
}

# Write region-merging outputs back onto the finemap result (only components
# ctwas actually returned). Rows fine-mapped again in a merged region have their
# molecular_id read off the id again, as the first pass's do.
# @noRd
.ctwasApplyMergeResult <- function(finemapResult, res) {
    regionUpdates <- compact(list(
        region_data = res$updated_region_data,
        region_info = res$updated_region_info,
        LD_map = res$updated_LD_map,
        snp_map = res$updated_snp_map
    ))
    list_assign(
        finemapResult,
        finemap_res = .ctwasFixMolecularId(
            res$updated_finemap_res,
            finemapResult$weights
        ),
        susie_alpha_res = .ctwasFixMolecularId(
            res$updated_susie_alpha_res,
            finemapResult$weights
        ),
        merge_res = res,
        !!!regionUpdates
    )
}

# Invoke a ctwas function with a fixed `args` list plus optional `extra`
# (typically the `...` collected by the wrapper). `extra` names that
# duplicate `args` names are silently dropped, so the wrapper's explicit
# arguments always win over caller-supplied `...`.
# @noRd
.ctwasInvoke <- function(fn, args, extra = list()) {
    # The bundle arrives as a MethodOptions record; everything below is ordinary
    # list work, so flatten it once at the boundary.
    extra <- as.list(extra)
    if (length(extra) > 0L) {
        deduped <- extra[setdiff(names(extra), names(args))]
        # `...` is forwarded uniformly to four different ctwas functions
        # (assemble_region_data / est_param / screen_regions /
        # finemap_regions). Restrict to fn's explicit formals so an arg
        # meant for a sibling step doesn't crash this one -- and so args
        # that fn would otherwise forward via its own `...` (e.g. into
        # susie_rss) don't bleed into incompatible downstream functions.
        accepted <- setdiff(names(formals(fn)), "...")
        args <- c(args, deduped[intersect(names(deduped), accepted)])
    }
    .ctwasQuietly(fn, args)
}

# Run ONLY ctwas's prefit EM step against `region_data` and return a
# param list shaped like ctwas::est_param normally produces. Used as
# the fallback path when est_param's accurate EM diverges to NaN on
# toy / underpowered data (matches the legacy ctwas_2 workaround).
# Calls ctwas's internal `fit_EM` (via getFromNamespace) for `niter`
# iterations -- run to convergence, because this fallback prior is the
# FINAL estimate (the accurate EM never ran), not a warm-up. `niter` here
# is the caller's full accurate-EM count, NOT niter_prefit; running only
# niter_prefit (a rough warm-up, e.g. 3) leaves the prior under-converged
# and depresses downstream gene PIPs. Then applies the same thin-adjustment
# to the SNP group_prior that est_param applies. p_single_effect is left as
# NA since the accurate EM never ran.
# @noRd
.ctwasFitPrefitEm <- function(
    region_data,
    niter,
    groupPriorVarStructure,
    thin,
    threads,
    extra = list()
) {
    fitEm <- utils::getFromNamespace("fit_EM", "ctwas")
    fitRegionData <- .ctwasPrefitRegionFilter(region_data, extra)
    fitArgs <- list(
        region_data = fitRegionData,
        niter = as.integer(niter),
        group_prior_var_structure = groupPriorVarStructure,
        ncore = threads
    )
    prefit <- .ctwasInvoke(fitEm, fitArgs, extra)
    adj <- .ctwasApplyThin(prefit$group_prior, prefit$group_size, thin)
    list(
        group_prior = adj$groupPrior,
        group_prior_var = prefit$group_prior_var,
        group_prior_iters = prefit$group_prior_iters,
        group_prior_var_iters = prefit$group_prior_var_iters,
        group_prior_var_structure = groupPriorVarStructure,
        group_size = adj$groupSize[names(adj$groupPrior)],
        p_single_effect = data.frame(
            region_id = names(region_data),
            p_single_effect = NA_real_,
            stringsAsFactors = FALSE
        )
    )
}

# Mirror ctwas::est_param's degenerate-region skip before the prefit fit_EM:
# drop regions with fewer than `min_var` total variables or fewer than
# `min_gene` genes (whose `sid` is unset), else ctwas::fit_EM errors inside
# extract_region_data ("regiondata$sid ... target is NULL") on a skipped
# region. Honors min_var / min_gene forwarded via `extra`.
# @noRd
.ctwasPrefitRegionFilter <- function(region_data, extra) {
    minVar <- if (!is.null(extra$min_var)) as.integer(extra$min_var) else 2L
    minGene <- if (!is.null(extra$min_gene)) as.integer(extra$min_gene) else 1L
    nGid <- lengths(map(region_data, "gid"))
    nSid <- lengths(map(region_data, "sid"))
    byVar <- if (minVar > 0L) {
        (nSid + nGid) >= minVar
    } else {
        rep(TRUE, length(region_data))
    }
    keep <- if (minGene > 0L) byVar & nGid >= minGene else byVar
    fitRegionData <- region_data[keep]
    if (length(fitRegionData) == 0L) {
        abort("No regions selected!")
    }
    fitRegionData
}

# Rescale the SNP group prior / size by `thin` (the SNP subsampling factor).
# @noRd
.ctwasApplyThin <- function(groupPrior, groupSize, thin) {
    if (thin == 1) {
        return(list(groupPrior = groupPrior, groupSize = groupSize))
    }
    list(
        groupPrior = .ctwasScaleSnpGroup(groupPrior, thin),
        groupSize = .ctwasScaleSnpGroup(groupSize, 1 / thin)
    )
}

# Scale the "SNP" entry of a per-group vector, leaving the molecular groups
# alone. fit_EM always reports a SNP group.
# @noRd
.ctwasScaleSnpGroup <- function(groups, factor) {
    replace(groups, "SNP", groups[["SNP"]] * factor)
}

# =============================================================================
# Internal helpers
# =============================================================================

# LD-sketch compatibility check for one region: every weight row in it is
# combined with every one of its GWAS rows. Thin wrapper over the shared
# `.requireMatchingRowLdRefs` helper (R/ld.R).
.ctwasRequireMatchingLdRefs <- function(tw, gwas) {
    nTw <- length(tw)
    nGwas <- length(gwas)
    .requireMatchingRowLdRefs(
        tw,
        gwas,
        qtlRows = rep(seq_len(nTw), each = nGwas),
        gwasRows = rep(seq_len(nGwas), times = nTw),
        pipelineName = "ctwasPipeline",
        label = "twasWeights"
    )
}

# The cTWAS runs `weightSelectionParam` asks of `twasWeights` (a flat weight
# source or a list of them), each modelling one weight per gene: under "best"
# one run, labelled "best", each gene on its own best method; under "all" or
# named methods one run per method, labelled by it. Each run holds its `label`
# and its chosen `cells` (.selectWeightCells()). The rows narrowed since the
# fits the choice and the runs read are warned about once.
# @noRd
.ctwasRunPlan <- function(twasWeights, weightSelectionParam, caller) {
    sources <- .ctwasWeightSources(twasWeights)
    chosen <- .ctwasChosenCells(sources, weightSelectionParam, caller)
    .ctwasCheckFitCoverage(sources, chosen, weightSelectionParam, caller)
    if (identical(weightSelectionParam$methods, "best")) {
        .ctwasAssertOnePerGene(chosen, caller)
        return(list(list(label = "best", cells = chosen)))
    }
    map(unique(chosen$method), .ctwasMethodRun, chosen = chosen)
}

# The weight sets of every source that `weightSelectionParam` keeps.
# @noRd
.ctwasChosenCells <- function(sources, weightSelectionParam, caller) {
    cells <- list_rbind(map(
        sources,
        .weightCells,
        weightSelectionParam = weightSelectionParam
    ))
    if (is.null(cells) || nrow(cells) == 0L) {
        abort(glue("{caller}: the weight sources carry no method entries."))
    }
    chosen <- .selectWeightCells(cells, weightSelectionParam, caller)
    if (nrow(chosen) == 0L) {
        abort(glue(
            "{caller}: no weights cleared WeightSelectionParam(cutoff = ",
            "{weightSelectionParam$cutoff}, pvalCutoff = ",
            "{weightSelectionParam$pvalCutoff})."
        ))
    }
    chosen
}

# The run of method `m`: its weight sets of `chosen`.
# @noRd
.ctwasMethodRun <- function(m, chosen) {
    list(label = m, cells = chosen[chosen$method == m, ])
}

# A "best" run models one weight per gene, so a gene left with several -- its
# methods carrying no CV metric to pick one by -- is refused.
# @noRd
.ctwasAssertOnePerGene <- function(chosen, caller) {
    gene <- str_c(chosen$study, chosen$context, chosen$trait, sep = "\r")
    nMethods <- map_int(split(chosen$method, gene), n_distinct)
    several <- names(nMethods)[nMethods > 1L]
    if (length(several) == 0L) {
        return(invisible(NULL))
    }
    first <- match(several[[1L]], gene)
    abort(c(
        glue(
            "{caller}: {length(several)} gene(s), such as ",
            "{chosen$trait[[first]]} ({chosen$study[[first]]}, ",
            "{chosen$context[[first]]}), hold several methods with no CV ",
            "metric to pick the best by, and a cTWAS run models one weight ",
            "per gene."
        ),
        i = str_c(
            "Name the method to use, WeightSelectionParam(methods = ",
            "\"<method>\"), or run cTWAS once per method with ",
            "methods = \"all\"."
        )
    ))
}

# The weight sources a `twasWeights` argument holds: itself when it is one
# flat source, else the non-NULL elements of its list.
# @noRd
.ctwasWeightSources <- function(twasWeights) {
    if (
        methods::is(twasWeights, "TwasWeights") ||
            methods::is(twasWeights, "QtlFineMappingResult")
    ) {
        return(list(twasWeights))
    }
    compact(as.list(twasWeights))
}

# Warn once about the narrowed rows whose fits the choice and the runs read:
# every row, when the choice reads the CV metrics, which describe the fits;
# else the chosen TwasWeights cells carrying a SuSiE fit, which renormalizes
# a weight whose variants shrank. A narrowed row's fit no longer lines up
# with its variants, so its weights go in without that step. Other weights
# are used as the row holds them, and a fine-mapping weight source supplies
# its topLoci effects, not its fit.
# @noRd
.ctwasCheckFitCoverage <- function(
    sources,
    chosen,
    weightSelectionParam,
    caller
) {
    if (.selectionReadsCv(weightSelectionParam)) {
        return(.fitCoverageCheckedAll(sources, caller))
    }
    .fitCoverageCheckedAll(
        map(
            keep(sources, methods::is, "TwasWeights"),
            .ctwasSusieFitRows,
            chosen = chosen
        ),
        caller
    )
}

# The chosen cells of TwasWeights `tw` whose stored fit is a SuSiE fit, or
# NULL when it has none.
# @noRd
.ctwasSusieFitRows <- function(tw, chosen) {
    cells <- .twasCells(tw)
    picked <- cells[is_in(.cellKey(cells), .cellKey(chosen)), ]
    susie <- map2_lgl(picked$row, picked$method, .ctwasCellHasSusieFit, x = tw)
    if (!any(susie)) NULL else .twasKeepCells(tw, picked[susie, ])
}

# @noRd
.ctwasCellHasSusieFit <- function(row, method, x) {
    .ctwasIsSusieFit(methodFits(.twrCellParts(x, row, method)))
}

# @noRd
.ctwasConcatChr <- function(pieces) {
    if (length(pieces) == 0L) {
        return(character(0))
    }
    as.character(list_c(pieces))
}

# Fail fast on the two cTWAS inputs. `gwasSumStats` defines the LD-block grid:
# one GwasSumStats whose elements are blocks. `twasWeights` may be a FLAT
# weight source (a single TwasWeights / QtlFineMappingResult, or a list of
# them) -- placed into blocks internally by `assembleCtwasInputs` -- or a
# pre-bucketed per-region named list.
# @noRd
.ctwasRequireNamedLists <- function(gwasSumStats, twasWeights) {
    if (!methods::is(gwasSumStats, "GwasSumStats")) {
        msg <- glue(
            "`gwasSumStats` must be a GwasSumStats whose elements are LD ",
            "blocks (got {class(gwasSumStats)[[1L]]}). Build one with ",
            "`loadGwasSumStatsFromManifest(..., region = <blocks>)`."
        )
        abort(msg)
    }
    okTw <- methods::is(twasWeights, "TwasWeights") ||
        methods::is(twasWeights, "QtlFineMappingResult") ||
        is.list(twasWeights)
    if (!okTw) {
        msg <- glue(
            "`twasWeights` must be a TwasWeights / QtlFineMappingResult ",
            "weight source (placed into LD blocks by its weights), or ",
            "a per-region named list keyed by region_id ",
            "(got {class(twasWeights)[[1L]]})."
        )
        abort(msg)
    }
}

# TRUE when `tw` is ALREADY a per-region named list (the pre-bucketed contract):
# a NAMED plain list (not an S4 collection) of weight collections. Its keys are
# validated against the block grid downstream (the `extra_tw_keys` check). A
# flat collection or an UNNAMED list is instead treated as a flat weight source
# to place by its weights.
# @noRd
.ctwasIsPreBucketed <- function(tw) {
    is.list(tw) &&
        !methods::is(tw, "DFrame") &&
        length(tw) > 0L &&
        !is.null(names(tw)) &&
        all(str_length(names(tw)) > 0L) &&
        all(map_lgl(tw, .ctwasIsWeightEntry))
}

# Combine a flat weight source into ONE collection. Accepts a single
# TwasWeights / QtlFineMappingResult, or a homogeneous list of one kind.
# @noRd
.ctwasCombineWeightSources <- function(weights) {
    if (
        methods::is(weights, "TwasWeights") ||
            methods::is(weights, "QtlFineMappingResult")
    ) {
        return(weights)
    }
    if (is.list(weights)) {
        weights <- weights[!map_lgl(weights, is.null)]
        if (length(weights) == 0L) {
            abort(
                "assembleCtwasInputs: `twasWeights` is an empty weight source."
            )
        }
        if (all(map_lgl(weights, methods::is, "TwasWeights"))) {
            parts <- unname(weights)
            return(exec(combineTwasWeights, !!!parts))
        }
        if (all(map_lgl(weights, methods::is, "QtlFineMappingResult"))) {
            parts <- unname(weights)
            return(exec(combineFineMappingResults, !!!parts))
        }
    }
    msg <- glue(
        "assembleCtwasInputs: `twasWeights` must be a TwasWeights or ",
        "QtlFineMappingResult (or a homogeneous list of one kind), or a ",
        "per-region named list keyed by region_id."
    )
    abort(msg)
}

# Place each gene (row) of a flat weight source into the LD block its weights
# are read with, from those weights: the block holding the largest sum of
# |weight| over the row's methods, the rule cTWAS homes a gene crossing blocks
# by (adjust_boundary_gene_region_assignment). Returns a region_id per row, NA
# when no weight falls in a block. This only divides the work: the block a gene
# is fine-mapped in is cTWAS's choice, from its weights (assemble_region_data).
# @noRd
#' @importFrom dplyr slice_max
.ctwasPlaceByWeights <- function(combined, gwasSumStats) {
    ids <- names(gwasSumStats)
    windows <- .ctwasBlockWindows(gwasSumStats, ids)
    sizes <- list_rbind(
        map(.ctwasWeightUnits(combined), .ctwasUnitWeightSizes, x = combined),
        ptype = tibble(
            row = integer(),
            chrom = character(),
            pos = integer(),
            size = double()
        )
    )
    placed <- sizes |>
        mutate(block = .ctwasBlockAt(.data$chrom, .data$pos, windows)) |>
        filter(!is.na(.data$block)) |>
        summarise(size = sum(.data$size), .by = all_of(c("row", "block"))) |>
        slice_max(.data$size, n = 1L, with_ties = FALSE, by = all_of("row"))
    ids[placed$block[match(seq_len(nrow(combined)), placed$row)]]
}

# Where one unit's weights lie and how large they are, a variant of several
# columns' weights counting them all. NULL when it holds none.
# @noRd
.ctwasUnitWeightSizes <- function(unit, x) {
    wr <- .rowResolveWeights(.ctwasUnitParts(x, unit))
    if (length(wr$variantIds) == 0L) {
        return(NULL)
    }
    parsed <- parseVariantId(wr$variantIds)
    tibble(
        row = unit$row,
        chrom = parsed$chrom,
        pos = as.integer(parsed$pos),
        size = rowSums(abs(as.matrix(wr$weights)))
    )
}

# The index of the block window holding each position, NA where none does: a
# window's stop is the half-open one region_info hands cTWAS, one past its end.
# Both sides take every chromosome either has, since one only a side has warns.
# @noRd
.ctwasBlockAt <- function(chrom, pos, windows) {
    windowChrom <- as.character(GenomicRanges::seqnames(windows))
    chroms <- union(windowChrom, chrom)
    GenomicRanges::findOverlaps(
        GenomicRanges::GRanges(
            factor(chrom, levels = chroms),
            IRanges::IRanges(pos, width = 1L)
        ),
        GenomicRanges::GRanges(
            factor(windowChrom, levels = chroms),
            IRanges::ranges(windows)
        ),
        select = "first"
    )
}

# The window of each LD block (.ctwasBlockWindow()), one per block. Built from
# its parts rather than by c(): each window knows only its own chromosome,
# and combining blocks on different chromosomes would warn that their
# seqlevels share nothing.
# @noRd
.ctwasBlockWindows <- function(gwasSumStats, ids) {
    windows <- unname(map2(gwasSumStats, ids, .ctwasBlockWindow))
    GenomicRanges::GRanges(
        map_chr(windows, .ctwasWindowSeqname),
        IRanges::IRanges(
            start = map_int(windows, GenomicRanges::start),
            end = map_int(windows, GenomicRanges::end)
        )
    )
}

# @noRd
.ctwasWindowSeqname <- function(w) {
    as.character(GenomicRanges::seqnames(w))
}

# The window of LD block `gss` (that block's GwasSumStats rows): its stored
# `region` -- the LD block it was split by, or the span of its variants as
# built -- not the span of the variants it holds now, so a gene or a SNP
# between two blocks' variants still lands in a block. A block id is only a
# name, so it gives no window. An error when the rows lie on several
# chromosomes.
# @noRd
.ctwasBlockWindow <- function(gss, regionId) {
    window <- range(.rowRegions(gss)$region, ignore.strand = TRUE)
    if (length(window) != 1L) {
        chrs <- as.character(GenomicRanges::seqnames(window))
        abort(glue(
            "ctwasPipeline: GwasSumStats block '{regionId}' spans multiple ",
            "chromosomes ({str_flatten(chrs, ', ')})."
        ))
    }
    window
}

# The genes placed in region `rid`, every row on that region's LD sketch, or
# NULL when nothing landed there.
# @noRd
.ctwasBucketForRegion <- function(rid, combined, home, gwasSumStats) {
    idx <- which(home == rid)
    if (length(idx) == 0L) {
        return(NULL)
    }
    .withRowLdRefs(combined[idx, ], ldReference(gwasSumStats[[rid]]))
}

# Bucket a flat weight source into a per-region named list keyed to the block
# grid, each gene where its weights lie. Each per-block sub-collection carries
# that block's GWAS LD sketch (the panel its weights are harmonized against, and
# what the downstream match-check expects).
# @noRd
.ctwasBucketWeights <- function(weights, gwasSumStats) {
    combined <- .ctwasCombineWeightSources(weights)
    home <- .ctwasPlaceByWeights(combined, gwasSumStats)
    unplaced <- sum(is.na(home))
    if (unplaced > 0L) {
        msg <- glue(
            "assembleCtwasInputs: {unplaced} gene(s) with no weight in any LD ",
            "block were dropped."
        )
        warn(msg)
    }
    out <- compact(set_names(
        map(
            names(gwasSumStats),
            .ctwasBucketForRegion,
            combined = combined,
            home = home,
            gwasSumStats = gwasSumStats
        ),
        names(gwasSumStats)
    ))
    if (length(out) == 0L) {
        msg <- glue(
            "assembleCtwasInputs: no gene placed into any LD block. ",
            "Check that the genes' weight variants fall inside the GWAS ",
            "blocks' windows (their `region`), on the same genome build."
        )
        abort(msg)
    }
    out
}

# Resolve `twasWeights` to a per-region named list: pass a pre-bucketed list
# through, otherwise place a flat weight source by its weights.
# @noRd
.ctwasResolveWeightBuckets <- function(twasWeights, gwasSumStats) {
    if (.ctwasIsPreBucketed(twasWeights)) {
        return(twasWeights)
    }
    .ctwasBucketWeights(twasWeights, gwasSumStats)
}


# The single GWAS (disease) study a cTWAS run models. cTWAS solves one disease
# per run (the z_snp carries a single z per SNP), so the input blocks must all
# reference the same GWAS study.
# @noRd
.ctwasGwasStudy <- function(gwasSumStats) {
    # Read the collection's own `study` column. map() over a GwasSumStats
    # iterates its ELEMENTS (per-block GRanges), which carry no study, so it
    # would silently yield NA for every block.
    allStudies <- unique(.ctwasStudyChr(gwasSumStats))
    studies <- allStudies[
        !is.na(allStudies) & str_length(allStudies) > 0L
    ]
    if (length(studies) == 0L) {
        return(NA_character_)
    }
    if (length(studies) > 1L) {
        msg <- glue(
            "ctwasPipeline: the input blocks reference multiple GWAS ",
            "studies ({str_flatten(studies, ', ')}); cTWAS models one ",
            "disease per run."
        )
        abort(msg)
    }
    studies
}

# The genes of a run, from their weight records: what cTWAS's strings do not
# carry apart, or carry only as one label. No id is parsed for this.
# @noRd
.ctwasGeneTable <- function(weights) {
    records <- unname(weights)
    tibble(
        id = names(weights) %||% character(0),
        study = map_chr(records, "study"),
        context = map_chr(records, "qtl_context"),
        trait = map_chr(records, "gene_name"),
        method = map_chr(records, "method")
    )
}

# cTWAS derives `molecular_id` from a gene id's first "|" field, which cuts a
# trait whose own id holds a "|" short. Restored from the gene's record;
# SNP-background rows are left untouched.
# @noRd
.ctwasFixMolecularId <- function(df, weights) {
    if (
        is.null(df) ||
            !is.data.frame(df) ||
            nrow(df) == 0L ||
            !all(is_in(c("id", "molecular_id"), names(df)))
    ) {
        return(df)
    }
    ids <- as.character(df$id)
    isGene <- is_in(ids, names(weights))
    if (!any(isGene)) {
        return(df)
    }
    mutate(
        df,
        molecular_id = replace(
            .data$molecular_id,
            isGene,
            map_chr(unname(weights[ids[isGene]]), "molecular_id")
        )
    )
}

# Subset a ctwas result frame (finemap_res / susie_alpha_res) to the rows whose
# `id` is in `ids`. Returns NULL when the frame is absent or nothing matches.
# @noRd
.ctwasSubsetById <- function(df, ids) {
    if (is.null(df)) {
        return(NULL)
    }
    sub <- df[is_in(as.character(df$id), ids), , drop = FALSE]
    if (nrow(sub) == 0L) NULL else `rownames<-`(sub, NULL)
}

# Subset a ctwas result frame to its SNP rows (type == "SNP"; anno_susie tags
# the non-gene background this way). Returns NULL when absent or none present.
# @noRd
.ctwasSubsetSnp <- function(df) {
    if (is.null(df) || is.null(df$type)) {
        return(NULL)
    }
    sub <- df[as.character(df$type) == "SNP", , drop = FALSE]
    if (nrow(sub) == 0L) NULL else `rownames<-`(sub, NULL)
}

# =============================================================================
# Ranging the cTWAS result payloads (spec 4.5)
# -----------------------------------------------------------------------------
# ctwas's finemap_res / susie_alpha_res carry no coordinates -- that is what
# anno_finemap_res(add_position = TRUE) exists for, and pecotmr never calls it.
# The coordinates are recoverable in-package without ctwas's mapping_table:
#
#   SNP rows   the `id` IS a variant id, so it renders its own range (4.1a)
#   gene rows  the run's weights carry that gene's trait_pos under its id
#              (falling back to the weight-variant span chrom/p0/p1)
#
# Without this a CtwasResult cannot answer "which genes did cTWAS implicate in
# this window" without the caller redoing the join by hand.
# =============================================================================

# Coordinates for every gene in a run, keyed by gene id.
# @noRd
.ctwasGeneCoords <- function(weights) {
    if (is.null(weights) || length(weights) == 0L) {
        return(NULL)
    }
    map(weights, .ctwasOneGeneCoord)
}

# @noRd
.ctwasOneGeneCoord <- function(w) {
    tp <- w$trait_pos
    if (methods::is(tp, "GRanges") && length(tp) == 1L) {
        return(list(
            chrom = as.character(seqnames(tp)),
            start = as.integer(start(tp)),
            end = as.integer(GenomicRanges::end(tp))
        ))
    }
    if (is.null(w$chrom) || is.null(w$p0)) {
        return(NULL)
    }
    list(
        chrom = withChrPrefix(as.character(w$chrom)),
        start = as.integer(w$p0),
        end = as.integer(w$p1)
    )
}

# A finemap / susieAlpha table as a GRanges, its original columns preserved as
# mcols. Returns the table unchanged when no row can be placed, so a payload
# that genuinely has no coordinates is not silently fabricated onto chr1; when
# only some rows can be, the rest are left out of the ranges, and named.
# @noRd
.ctwasRangePayload <- function(df, geneCoords) {
    if (is.null(df) || nrow(df) == 0L) {
        return(df)
    }
    ids <- as.character(df$id)
    coord <- map(seq_along(ids), .ctwasRowCoord, ids = ids, gc = geneCoords)
    placed <- !map_lgl(coord, is.null)
    if (!any(placed)) {
        return(df)
    }
    .ctwasWarnUnplaced(ids[!placed])
    gr <- GenomicRanges::GRanges(
        map_chr(coord, .ctwasCoordField, field = "chrom"),
        IRanges::IRanges(
            start = map_int(
                coord,
                .ctwasCoordField2,
                field = "start"
            ),
            end = map_int(
                coord,
                .ctwasCoordField2,
                field = "end"
            )
        )
    )
    S4Vectors::`mcols<-`(
        gr,
        value = S4Vectors::DataFrame(df, check.names = FALSE)
    )[placed]
}

# Rows with neither a gene coordinate nor a chrom:pos id cannot be ranged.
# @noRd
.ctwasWarnUnplaced <- function(ids) {
    if (length(ids) == 0L) {
        return(invisible(NULL))
    }
    shown <- unique(ids)
    warn(c(
        glue(
            "cTWAS: {length(ids)} result row(s) could not be placed on the ",
            "genome and are left out of the ranged results: they have no ",
            "gene coordinate and their id is not a chrom:pos variant id."
        ),
        set_names(utils::head(shown, 5L), rep("*", min(5L, length(shown))))
    ))
}

# @noRd
.ctwasRowCoord <- function(i, ids, gc) {
    id <- ids[[i]]
    if (!is.null(gc) && is_in(id, names(gc))) {
        return(gc[[id]])
    }
    parsed <- try_fetch(parseVariantId(id), error = function(cnd) NULL)
    if (is.null(parsed) || is.na(parsed$chrom[[1L]])) {
        return(NULL)
    }
    list(
        chrom = parsed$chrom[[1L]],
        start = as.integer(parsed$pos[[1L]]),
        end = as.integer(parsed$pos[[1L]])
    )
}

# @noRd
.ctwasCoordField <- function(co, field) {
    if (is.null(co)) "chrUnplaced" else as.character(co[[field]])
}

# @noRd
.ctwasCoordField2 <- function(co, field) {
    if (is.null(co)) 1L else as.integer(co[[field]])
}

# Build a CtwasResultEntry from a finemap + susieAlpha slice, recording the
# run's param + region_info on it.
# @noRd
.ctwasMkEntry <- function(fm, sa, runResult) {
    geneCoords <- .ctwasGeneCoords(runResult$weights)
    CtwasResultEntry(
        posteriors = .ctwasRangePayload(.ctwasReportGroups(fm), geneCoords),
        susieAlpha = .ctwasRangePayload(.ctwasReportGroups(sa), geneCoords),
        groupPriors = runResult$param,
        regionInfo = .ctwasRangeRegionInfo(runResult$region_info)
    )
}

# A result table with its groups named as results name them
# (.ctwasGroupLabel()).
# @noRd
.ctwasReportGroups <- function(df) {
    if (is.null(df) || !is_in("group", names(df))) {
        return(df)
    }
    mutate(df, group = .ctwasGroupLabel(.data$group))
}

# cTWAS's name for a group, `context|type`, as results name it. The "|" is the
# one cTWAS's EM needs (.ctwasGeneEntry()), of no use to anyone reading a
# result.
# @noRd
.ctwasGroupLabel <- function(groups) {
    str_replace(groups, fixed("|"), "_")
}

# region_info as a GRanges. pecotmr builds this table itself, each block's
# window as cTWAS's half-open [start, stop), so the range ends one base before
# `stop` -- the window again (a merged region, the span of the windows it
# merged).
# @noRd
.ctwasRangeRegionInfo <- function(ri) {
    if (is.null(ri) || nrow(ri) == 0L) {
        return(ri)
    }
    cols <- names(ri)
    if (!all(is_in(c("chrom", "start", "stop"), cols))) {
        return(ri)
    }
    gr <- GenomicRanges::GRanges(
        withChrPrefix(as.character(ri$chrom)),
        IRanges::IRanges(
            start = as.integer(ri$start),
            end = as.integer(ri$stop) - 1L
        )
    )
    S4Vectors::`mcols<-`(
        gr,
        value = S4Vectors::DataFrame(ri, check.names = FALSE)
    )
}

# Decompose one cTWAS run (a `finemapCtwasRegions` output) into per-context
# row-specs for a CtwasResult. The row skeleton comes from the ASSEMBLED weights
# (so every modeled (study, context) appears even if no gene reached
# fine-mapping); each row's finemap / susieAlpha payloads are the subsets whose
# gene id belongs to that context. Multi-context runs are annotated with the
# shared `jointContexts` set and share the jointly-estimated `param`.
#
# `keepSnps` (default FALSE) additionally retains the context-agnostic SNP
# background as ONE extra row (study = context = "SNP"), mirroring cTWAS's own
# "SNP" group in `group_prior`. Kept off by default because the SNP rows are the
# null background and bloat the structured gene-level result; when on, the full
# ctwas run is reconstructable from `ctwasPosteriors()` / `susieAlpha()`.
# @noRd
.ctwasRunToRows <- function(runResult, gwasStudy, method, keepSnps = FALSE) {
    genes <- .ctwasGeneTable(runResult$weights)
    if (nrow(genes) == 0L) {
        return(list())
    }
    runs <- distinct(genes, .data$study, .data$context)
    joint <- list(
        studies = .ctwasJointLabel(runs$study),
        contexts = .ctwasJointLabel(runs$context)
    )
    fmDf <- .ctwasAsDf(runResult$finemap_res)
    saDf <- .ctwasAsDf(runResult$susie_alpha_res)
    scored <- list_assign(runResult, param = .ctwasScoredParam(runResult))
    contextRows <- map2(
        runs$study,
        runs$context,
        .ctwasContextRow,
        genes = genes,
        gwasStudy = gwasStudy,
        method = method,
        joint = joint,
        fmDf = fmDf,
        saDf = saDf,
        runResult = scored
    )
    snpRow <- if (!keepSnps) {
        NULL
    } else {
        .ctwasSnpRow(gwasStudy, method, joint, fmDf, saDf, scored)
    }
    c(contextRows, compact(list(snpRow)))
}

# The studies or contexts a joint run spans, NA for a run over one.
# @noRd
.ctwasJointLabel <- function(values) {
    members <- sort(unique(values))
    if (length(members) > 1L) str_flatten(members, ",") else NA_character_
}

# The run's priors with what results report beside them: each group's study,
# context and type, and cTWAS's summary of the priors.
# @noRd
.ctwasScoredParam <- function(runResult) {
    param <- runResult$param
    if (is.null(param$group_prior)) {
        return(param)
    }
    c(
        param,
        list(
            group_info = .ctwasGroupInfo(runResult$weights),
            group_summary = .ctwasGroupSummary(param, runResult$gwas_n)
        )
    )
}

# Each group's study, context and molecular type, under cTWAS's name for it.
# @noRd
.ctwasGroupInfo <- function(weights) {
    records <- unname(weights)
    genes <- tibble(
        group = map_chr(records, .ctwasEntryGroup),
        study = map_chr(records, "study"),
        context = map_chr(records, "qtl_context"),
        type = map_chr(records, "type")
    )
    bind_rows(
        distinct(genes),
        tibble(group = "SNP", study = "SNP", context = "SNP", type = "SNP")
    )
}

# as.data.frame(), passing NULL through.
# @noRd
.ctwasAsDf <- function(x) {
    if (is.null(x)) NULL else as.data.frame(x)
}

# One (gwasStudy, study, context, method) row-spec.
# @noRd
.ctwasContextRow <- function(
    study,
    context,
    genes,
    gwasStudy,
    method,
    joint,
    fmDf,
    saDf,
    runResult
) {
    ids <- genes$id[genes$study == study & genes$context == context]
    list(
        gwasStudy = gwasStudy,
        study = study,
        context = context,
        method = method,
        jointStudies = joint$studies,
        jointContexts = joint$contexts,
        entry = .ctwasMkEntry(
            .ctwasSubsetById(fmDf, ids),
            .ctwasSubsetById(saDf, ids),
            runResult
        )
    )
}

# The SNP-level row-spec (study = context = "SNP"), or NULL when no SNP-level
# finemap / susie-alpha rows exist.
# @noRd
.ctwasSnpRow <- function(gwasStudy, method, joint, fmDf, saDf, runResult) {
    snpFm <- .ctwasSubsetSnp(fmDf)
    snpSa <- .ctwasSubsetSnp(saDf)
    if (is.null(snpFm) && is.null(snpSa)) {
        return(NULL)
    }
    list(
        gwasStudy = gwasStudy,
        study = "SNP",
        context = "SNP",
        method = method,
        jointStudies = joint$studies,
        jointContexts = joint$contexts,
        entry = .ctwasMkEntry(snpFm, snpSa, runResult)
    )
}

# Assemble accumulated per-run row-specs (from .ctwasRunToRows) into a single
# CtwasResult; a joint column is omitted when no row is joint on it. Shared by
# ctwasPipeline (across methods) and asCtwasResult.
# @noRd
.ctwasRowsToResult <- function(rows) {
    if (length(rows) == 0L) {
        msg <- glue(
            "cTWAS: no genes were modeled (the weight source produced ",
            "no usable gene weights)."
        )
        abort(msg)
    }
    jointStudies <- map_chr(rows, "jointStudies")
    jointContexts <- map_chr(rows, "jointContexts")
    CtwasResult(
        gwasStudy = map_chr(rows, "gwasStudy"),
        studyName = map_chr(rows, "study"),
        context = map_chr(rows, "context"),
        method = map_chr(rows, "method"),
        entry = map(rows, "entry"),
        jointStudies = if (any(!is.na(jointStudies))) jointStudies else NULL,
        jointContexts = if (any(!is.na(jointContexts))) jointContexts else NULL
    )
}

# The single GWAS study a finemap result models, read from z_snp$study (which
# `.ctwasBuildZSnp` fills in per row). Errors on multiple; NA when absent.
# @noRd
.ctwasGwasStudyFromZSnp <- function(zSnp) {
    if (is.null(zSnp) || is.null(zSnp$study)) {
        return(NA_character_)
    }
    allStudies <- unique(as.character(zSnp$study))
    s <- allStudies[!is.na(allStudies) & str_length(allStudies) > 0L]
    if (length(s) == 0L) {
        return(NA_character_)
    }
    if (length(s) > 1L) {
        msg <- glue(
            "asCtwasResult: z_snp references multiple GWAS studies ",
            "({str_flatten(s, ', ')}); cTWAS models one disease per run."
        )
        abort(msg)
    }
    s
}

# The single weight method a finemap result was built for, read from the
# assembled weights. Errors on a mix (the granular path is one method/run).
# @noRd
.ctwasMethodFromWeights <- function(weights) {
    if (is.null(weights) || length(weights) == 0L) {
        msg <- glue(
            "asCtwasResult: the finemap result carries no weights to derive a ",
            "method from."
        )
        abort(msg)
    }
    m <- unique(.ctwasGeneTable(weights)$method)
    if (length(m) != 1L) {
        msg <- glue(
            "asCtwasResult: the finemap result mixes weight methods ",
            "({str_flatten(m, ', ')}); expected one per run."
        )
        abort(msg)
    }
    m
}

#' @importFrom checkmate assertFlag
#' @title Structure a granular cTWAS finemap result as a CtwasResult
#' @description Decompose the raw list returned by
#'   \code{\link{finemapCtwasRegions}} (optionally after
#'   \code{\link{mergeCtwasBoundaryRegions}}) into the structured, per-(study,
#'   context) \code{\link{CtwasResult}} -- the same decomposition
#'   \code{\link{ctwasPipeline}} applies to its one-shot output, exposed for the
#'   granular \code{assembleCtwasInputs} \eqn{\to} \code{estCtwasGroupPriors}
#'   \eqn{\to} \code{screenCtwasRegions} \eqn{\to} \code{finemapCtwasRegions}
#'   path. The GWAS study is read from \code{z_snp$study}, and the method
#'   label from the run's \code{resolvedMethod} (\code{"best"}, or its one
#'   method; for a result made before that was carried, the one method of its
#'   weights).
#' @param finemapResult A list from \code{\link{finemapCtwasRegions}} or
#'   \code{\link{mergeCtwasBoundaryRegions}}.
#' @param keepSnps Logical (length 1). Retain the context-agnostic SNP
#'   background as one extra \code{study = context = "SNP"} row. Default
#'   \code{FALSE}. See \code{\link{ctwasPipeline}}.
#' @return A \code{\link{CtwasResult}}.
#' @seealso \code{\link{ctwasPipeline}}, \code{\link{finemapCtwasRegions}}
#' @examples
#' data(ctwasFinemapExample)
#' asCtwasResult(ctwasFinemapExample)
#' @export
asCtwasResult <- function(finemapResult, keepSnps = FALSE) {
    assertFlag(keepSnps)
    gwasStudy <- .ctwasGwasStudyFromZSnp(finemapResult$z_snp)
    # The run's label, as assembleCtwasInputs() records it; a result made
    # before it was carried names its one method in its weights.
    method <- finemapResult$resolvedMethod %||%
        .ctwasMethodFromWeights(finemapResult$weights)
    rows <- .ctwasRunToRows(
        finemapResult,
        gwasStudy = gwasStudy,
        method = method,
        keepSnps = keepSnps
    )
    .ctwasRowsToResult(rows)
}

# Cut a weight source to the weight sets `chosen` names (.ctwasRunPlan()), or
# NULL when it holds none of them. Used to enforce the "one ctwas gene per
# (study, context, trait)" semantics -- the legacy pipeline fed a single
# best-CV-method weight per gene; the S4 TwasWeights may carry many methods,
# but a run should only see one per gene. A TwasWeights is cut to its chosen
# cells, a fine-mapping result to its chosen rows.
# @noRd
.ctwasKeepChosen <- function(tw, chosen) {
    if (methods::is(tw, "TwasWeights")) {
        cells <- .twasCells(tw)
        kept <- cells[is_in(.cellKey(cells), .cellKey(chosen)), ]
        if (nrow(kept) == 0L) {
            return(NULL)
        }
        return(.twasKeepCells(tw, kept))
    }
    keep <- which(is_in(
        .cellKey(.weightRowIds(tw)),
        .cellKey(chosen)
    ))
    if (length(keep) == 0L) {
        return(NULL)
    }
    # Row subset carries every column forward (joint* / region / ...); the old
    # hand-listed rebuild silently dropped them.
    methods::initialize(tw[keep, ], ldReference = ldReference(tw))
}

# Match variant ids against the LD panel's ids the way every other LD consumer
# in the package does -- by (chrom, pos, allele) with ref/alt swaps tolerated --
# rather than by exact string identity.
#
# ctwas already harmonizes TWAS weights this way (see
# `.ctwasHarmonizeWeights()` below). The GWAS side used plain string joins, so
# one variant could be kept for the weights and dropped for the z-scores.
# `removeStrandAmbiguous = FALSE` matches `.ldFromRefMatch()`, the shared
# entry point: the panel defines the frame, so an A/T variant is not discarded
# merely for being palindromic.
# @noRd
.ctwasMatchToPanel <- function(ids, panelIds) {
    matchVariants(ids, panelIds, removeStrandAmbiguous = FALSE)
}

# Put the GWAS z-scores in the LD panel's allele frame.
#
# The z-scores arrive in the GWAS's own frame, while `R`, `variance` and the
# harmonized weights are all in the panel's. Where the two spell a variant with
# its alleles swapped, an exact-string join drops it silently -- it is present
# on both sides, so nothing reports a missing variant and the gene just loses an
# instrument. Matching allele-aware keeps it, and negating its z (with A1/A2
# swapped to match) is what makes the retained z agree in sign with the LD it
# gets modelled against. Variants with no panel match are left alone; ctwas
# drops them itself when it intersects with `snp_map`.
# @noRd
.ctwasHarmonizeZToPanel <- function(zSnp, panelIds) {
    if (nrow(zSnp) == 0L || length(panelIds) == 0L) {
        return(zSnp)
    }
    m <- .ctwasMatchToPanel(zSnp$id, as.character(panelIds))
    if (length(m$idxA) == 0L) {
        return(zSnp)
    }
    relabelled <- mutate(
        zSnp,
        id = replace(.data$id, m$idxA, as.character(panelIds)[m$idxB])
    )
    flip <- m$idxA[m$sign < 0]
    if (length(flip) == 0L) {
        return(relabelled)
    }
    # A swapped variant counts the other allele: negate z and exchange the
    # allele columns. Both replacements read the ORIGINAL frame, so the pair
    # swaps rather than each taking the other's already-swapped value.
    mutate(
        relabelled,
        z = replace(.data$z, flip, -relabelled$z[flip]),
        A1 = replace(.data$A1, flip, relabelled$A2[flip]),
        A2 = replace(.data$A2, flip, relabelled$A1[flip])
    )
}

# @noRd
.ctwasZSnpPiece <- function(i, gwasSumStats) {
    df <- .entryToSumstatDf(gwasSumStats[[i]], keepChrPrefix = FALSE)
    tibble(
        id = df$variant_id,
        chrom = as.integer(df$chrom),
        pos = df$pos,
        A1 = df$A1,
        A2 = df$A2,
        z = df$z,
        study = as.character(gwasSumStats$study)[[i]],
        N = df$N
    )
}

# Build the per-variant Z data.frame ctwas expects from a GwasSumStats.
# Stacks each study row's GRanges via the shared `.entryToSumstatDf`
# helper (R/sumstatsQc.R), then projects to ctwas's column shape and
# bolts on the `study` column ctwas uses to disambiguate stacked rows.
# @noRd
.ctwasBuildZSnp <- function(gwasSumStats, panelIds) {
    pieces <- map(
        seq_len(nrow(gwasSumStats)),
        .ctwasZSnpPiece,
        gwasSumStats = gwasSumStats
    )
    .ctwasHarmonizeZToPanel(bind_rows(pieces), panelIds)
}

# The single-row region_info for block `regionId`: its chromosome and its
# window (.ctwasBlockWindow()), which cTWAS reads as the half-open
# [start, stop) -- so `stop` is one past the window's last base, and every
# variant of the block is assigned to it. A block with no variants, or one
# spanning several chromosomes, is an error.
# @noRd
.ctwasBuildSingleRegionInfo <- function(regionId, gss) {
    # The window is the block's own -- NOT the LD sketch's span. When many
    # blocks share one whole-chromosome LD payload (the common one-file-per-chr
    # layout), the LD reference's variants span the entire chromosome, so every
    # region would collapse to the same whole-chromosome [start, stop] and every
    # SNP would be assigned to every region (inflating SNP group_size N-fold and
    # diluting the gene prior to ~0).
    #
    # Emptiness is checked FIRST: a block with no variants has nothing for
    # cTWAS to fit, and the window of one empty from the start is only a
    # placeholder.
    if (sum(lengths(gss)) == 0L) {
        msg <- glue(
            "ctwasPipeline: GwasSumStats block '{regionId}' has no variants ",
            "for cTWAS to fit."
        )
        abort(msg)
    }
    window <- .ctwasBlockWindow(gss, regionId)
    tibble(
        region_id = regionId,
        chrom = as.integer(str_remove(
            as.character(GenomicRanges::seqnames(window)),
            regex("^chr", ignore_case = TRUE)
        )),
        start = GenomicRanges::start(window),
        stop = GenomicRanges::end(window) + 1L
    )
}

# Per-block SNP info table (chrom, id, pos, alt, ref). ctwas requires
# these exact column names (read_snp_info_files asserts them). `alt`
# maps to A1 (effect allele) and `ref` to A2. NOTE: cTWAS requires an integer
# chrom, so the as.integer() cast below is intentional (an output-boundary
# format requirement) and this path is autosomal-only -- X/Y/MT are not
# supported by the downstream cTWAS model.
# @noRd
.ctwasSnpInfoForBlock <- function(gwasLd) {
    gr <- .ldRefRanges(gwasLd)
    mc <- S4Vectors::mcols(gr)
    # `.ldRefMatchIds()`, not the raw SNP label: this id becomes the R
    # dimnames, the `variance` names, `panelSnps` and the weight-harmonization
    # reference id all at once, and every one of those is compared against a
    # harmonized (reference-allele) id from the GWAS or the weights. A panel
    # entry spelling a tag where an allele belongs would fail all four
    # comparisons and silently drop the variant from the analysis.
    #
    # Base data.frame (not tibble): ctwas indexes snp_map positionally
    # (df[, "pos"] -> vector for region-bound min/max); a tibble column is a
    # 1-col list and errors min().
    data.frame(
        chrom = as.integer(canonChrom(.ldRefChrom(gwasLd))),
        id = .ldRefMatchIds(gwasLd),
        pos = as.integer(GenomicRanges::start(gr)),
        alt = as.character(mc$A1),
        ref = as.character(mc$A2),
        stringsAsFactors = FALSE
    )
}

# Everything about a block's LD panel that is cheap, plus the sketch it came
# from so the expensive part can be built later. Returns a list with:
#   snpInfo  : ctwas-shaped per-block table (chrom, id, pos, alt, ref)
#              -- both the snp_map element and the snpinfo loader return.
#   variance : named numeric vector of per-variant genotype variance from
#              the LD reference (see .ldRefVariance()), used to scale
#              non-standardized TWAS weights to the correlation scale ctwas
#              expects.
#   sketch   : the block's narrowed LD sketch, retained because the n x n
#              correlation is deliberately NOT built here.
#
# ctwas pulls a region's LD only while fine-mapping, and only for the regions
# that survive screening, so building every block's matrix up front spent
# O(n^2) apiece on panels that are frequently never read at all.
# @noRd
.ctwasPanelFor <- function(gwasLd) {
    # Share the validator with `.ldFromRef()`, the entry point every other
    # pipeline uses: skipping the guard meant a NULL or non-panel sketch
    # surfaced as "unable to find an inherited method for 'snpInfo'"
    # instead of saying the LD reference was missing.
    .ldFromRefValidate(gwasLd, "ctwasPipeline")
    snpInfo <- .ctwasSnpInfoForBlock(gwasLd)
    list(
        snpInfo = snpInfo,
        variance = set_names(
            .ldRefVariance(gwasLd, seq_len(nrow(snpInfo))),
            snpInfo$id
        ),
        sketch = gwasLd
    )
}

# The block's full n x n correlation, built on request. A panel that already
# carries one -- a caller who computed it, or a payload it travelled with --
# is used as it stands rather than rebuilt.
# @noRd
.ctwasPanelLd <- function(panel) {
    if (!is.null(panel$R)) {
        return(panel$R)
    }
    .ctwasPanelLdFor(
        panel,
        seq_len(nrow(panel$snpInfo)),
        panel$snpInfo$id
    )
}

# @noRd
.ctwasPanelLdFor <- function(panel, idx, ids) {
    `dimnames<-`(.ldRefLdAt(panel$sketch, idx), list(ids, ids))
}

# Harmonize TWAS weight variants against the LD reference panel. Same
# allele-matching semantics as the GWAS-side `harmonizeAlleles` flow:
# match by (chrom, pos), accept exact A1/A2 frame, sign-flip the weight
# when alleles are swapped, drop unmatched / strand-ambiguous variants.
# Returns a data.frame with columns:
#   variant_id : canonical (panel-frame) variant ID
#   w          : sign-flipped weight aligned to the panel's A1 frame
#   origIdx    : index back into the entry's original variantIds vector
#                (used by SuSiE renormalization to slice mu / lbf)
# Returns NULL when the entry has no variants in common with the panel.
# @noRd
.ctwasHarmonizeWeights <- function(origVids, origW, refVariants) {
    parsed <- try_fetch(parseVariantId(origVids), error = function(cnd) NULL)
    if (is.null(parsed) || nrow(parsed) == 0L) {
        return(NULL)
    }
    targetDf <- tibble(
        chrom = parsed$chrom,
        pos = as.integer(parsed$pos),
        A2 = as.character(parsed$A2),
        A1 = as.character(parsed$A1),
        w = as.numeric(origW),
        origIdx = seq_along(origVids)
    )
    res <- try_fetch(
        harmonizeAlleles(
            targetData = targetDf,
            refVariants = refVariants,
            colToFlip = "w",
            matchMinProp = 0,
            removeUnmatched = TRUE,
            removeStrandAmbiguous = TRUE
        ),
        error = function(cnd) NULL
    )
    if (is.null(res)) {
        return(NULL)
    }
    res$harmonizedData
}

# Does the entry's `fits` slot carry a SuSiE-shape intermediate (lbf,
# mu, X_column_scale_factors)? Used to gate the renormalization branch.
# @noRd
.ctwasIsSusieFit <- function(fits) {
    if (is.null(fits)) {
        return(FALSE)
    }
    # `alpha`, not `lbf_variable`: renormalization restricts the stored
    # per-effect alpha, which already carries the fit's prior.
    # Not X_column_scale_factors: a stored fit drops them when all are 1.
    needed <- c("alpha", "mu")
    all(is_in(needed, names(fits)))
}

# Renormalize SuSiE TWAS weights over the kept variant set. When some
# variants got dropped by allele harmonization / panel intersection,
# the posterior `alpha` values from the original fit no longer sum to
# 1 over the kept variants. We restrict each effect's stored `alpha` to
# `keptIdx` and renormalize it (.ctwasRenormAlpha()), sign-flip the columns
# of `mu[, keptIdx]` to match the panel's allele frame (carrying over the
# per-variant sign flip already applied to `harmonizedW`), and recompute the
# per-variant weight as
# `colSums(alpha * mu_subset) / X_column_scale_factors_subset`.
# Returns the new weight vector (length = length(keptIdx)), or NULL if
# the fit's dimensions don't line up with the entry's variantIds.
# @noRd
.ctwasRenormalizeSusieWeights <- function(
    fits,
    origVids,
    origW,
    keptIdx,
    harmonizedW
) {
    if (!.ctwasRenormalizable(fits, length(origVids))) {
        return(NULL)
    }
    # Per-variant sign flip applied by allele harmonization. NaN signs
    # (origW == 0) default to +1.
    rawSign <- sign(harmonizedW / origW[keptIdx])
    signFlip <- replace(rawSign, !is.finite(rawSign), 1)
    newAlpha <- .ctwasRenormAlpha(as.matrix(fits[["alpha"]]), keptIdx)
    if (is.null(newAlpha)) {
        return(NULL)
    }
    muSub <- sweep(fits[["mu"]][, keptIdx, drop = FALSE], 2L, signFlip, `*`)
    rawScale <- .fitScaleFactors(fits, length(origVids))[keptIdx]
    # Guard against zero scale factors (shouldn't happen in practice).
    xColSub <- replace(rawScale, rawScale == 0, 1)
    as.numeric(colSums(newAlpha * muSub) / xColSub)
}

# Whether a fit's weights can be recomputed from its alpha and mu over a
# subset of its `p` variants. Fields are READ with `[[`: `$` on a list falls
# back to prefix matching, so `fits$mu` would silently return `mu2` on a fit
# that has no `mu`.
# @noRd
.ctwasRenormalizable <- function(fits, p) {
    alpha <- fits[["alpha"]]
    mu <- fits[["mu"]]
    if (is.null(alpha) || is.null(mu)) {
        return(FALSE)
    }
    # susieInf / susieAsh carry an infinitesimal term: coef.susie is
    # colSums(alpha * mu) / scale + theta / scale. Recomputing the weight from
    # alpha and mu alone would silently drop theta, changing what the weight
    # means. Leave those fits to the plain subset-w-and-R path, which is
    # unbiased and only loses power.
    if (!is.null(fits[["theta"]]) || !is.null(fits[["omega_weights"]])) {
        return(FALSE)
    }
    # An mvSuSiE mu is L x p x R. .ctwasRenormFit() cuts it to the row's own
    # condition when the row records its position in the fit; one that
    # records none takes that same path.
    if (length(dim(mu)) > 2L) {
        return(FALSE)
    }
    # A fit-vs-entry dimension mismatch -- including a null_weight fit, whose
    # alpha carries an extra column -- is skipped rather than mis-sliced.
    # Absent scale factors are 1s (.fitScaleFactors()), which always fit.
    xCol <- fits[["X_column_scale_factors"]]
    nScale <- if (is.null(xCol)) p else length(xCol)
    all(c(ncol(as.matrix(alpha)), ncol(as.matrix(mu)), nScale) == p)
}

# Restrict each single-effect posterior to `keptIdx` and renormalize. The
# stored alpha already carries the fit's prior (alpha is proportional to
# pi * BF), so restricting and renormalizing is exact for any prior weights,
# whereas rebuilding alpha from lbf_variable substitutes a uniform one. Log
# space keeps an effect whose retained mass has underflowed from collapsing to
# NaN; NULL when an effect has no retained mass at all, so the caller falls
# back to the un-renormalized weights.
# @noRd
.ctwasRenormAlpha <- function(alpha, keptIdx) {
    logSub <- log(alpha[, keptIdx, drop = FALSE])
    rowMax <- apply(logSub, 1L, max)
    if (any(!is.finite(rowMax))) {
        return(NULL)
    }
    weights <- exp(sweep(logSub, 1L, rowMax, `-`))
    weights / rowSums(weights)
}

# Build the weights list ctwas expects: keyed by per-tuple gene id,
# each element a list with wgt (variants x 1 matrix; rownames = SNP id),
# wgt_var (the weight variance w'Rw), and gene metadata. ctwas's
# compute_gene_z would compute z.gene = crossprod(wgt, z.s) /
# sqrt(t(wgt) %*% R_wgt %*% wgt) from a per-gene LD submatrix R_wgt; the
# pipeline computes the same z from wgt_var instead (.ctwasGeneZ()) and
# supplies it, so wgt must still be a numeric matrix (not a vector) but no
# weight LD matrix is formed.
#
# The panel variants a gene's cis span may reach.
#
# ctwas's compute_gene_z asserts every weight variant exists in the block's
# z_snp$id. An LD sketch covering more than the block (e.g. a whole-chrom
# PLINK2) leaks variants outside it, so intersect with the caller's GWAS
# sumstats variant set when provided.
#
# Allele-aware, not `intersect()`: a swapped spelling is the same variant,
# and dropping it here silently shortens the gene's cis span.
# @noRd
.ctwasPanelSnpsForGwas <- function(panelSnps, gwasSnpIds) {
    if (is.null(gwasSnpIds)) {
        return(panelSnps)
    }
    m <- .ctwasMatchToPanel(panelSnps, as.character(gwasSnpIds))
    panelSnps[sort(m$idxA)]
}

# What preparing a gene reads, built once for the run since it is the same for
# every gene. A gene is read whole, whichever blocks its weights reach, so
# anything about it comes off the whole panel rather than off one block's --
# `spanInfo` for the coordinates and allele harmonization, `spanVariance` for
# the scaling -- and `gwasSnpIds` is the GLOBAL set of GWAS variants. `cutoffs`
# are .ctwasPruningCutoffs()'s.
# @noRd
.ctwasGeneContext <- function(
    spanInfo,
    spanVariance,
    gwasSnpIds,
    fineMappingResult,
    cutoffs
) {
    list(
        spanInfo = spanInfo,
        inGwas = is_in(
            spanInfo$id,
            .ctwasPanelSnpsForGwas(spanInfo$id, gwasSnpIds)
        ),
        refVariants = .ctwasRefVariants(spanInfo),
        spanVariance = spanVariance,
        varianceRows = match(spanInfo$id, names(spanVariance)),
        positions = .ctwasPositionIndex(spanInfo$pos),
        fineMappingResult = fineMappingResult,
        cutoffs = cutoffs
    )
}

# The panel's rows by position, for finding a gene's by searching the sorted
# positions rather than passing over the whole panel.
# @noRd
.ctwasPositionIndex <- function(pos) {
    rows <- order(pos)
    list(rows = rows, pos = pos[rows])
}

# The panel rows at these positions, in panel order.
# @noRd
.ctwasRowsAtPositions <- function(index, pos) {
    wanted <- unique(pos[!is.na(pos)])
    first <- findInterval(wanted, index$pos, left.open = TRUE) + 1L
    last <- findInterval(wanted, index$pos)
    sort(index$rows[sequence(last - first + 1L, first)])
}

# The part of the panel a gene reads: its rows at the positions of the gene's
# variants, the only ones those can be. Nothing about preparing a gene then
# depends on the size of the panel.
# @noRd
.ctwasPanelNear <- function(parts, ctx) {
    fields <- .variantIdFields(.rowResolveWeights(parts)$variantIds)
    rows <- .ctwasRowsAtPositions(ctx$positions, as.integer(fields[, 3L]))
    list(
        spanInfo = ctx$spanInfo[rows, , drop = FALSE],
        panelSnps = ctx$spanInfo$id[rows][ctx$inGwas[rows]],
        refVariants = ctx$refVariants[rows, ],
        spanVariance = ctx$spanVariance[ctx$varianceRows[rows]]
    )
}

# The genes of one weight source, prepared up to their weight variances from
# `ctx` (.ctwasGeneContext()). Variants absent from the panel are dropped from
# that gene's row set.
#
# Genes are prepared here and given their weight variances afterwards, across
# all blocks (.ctwasGeneEntries()), so each block's LD is read once.
# @noRd
.ctwasPreparedGenes <- function(twasWeights, ctx) {
    compact(map(
        .ctwasWeightUnits(twasWeights),
        .ctwasGenePrepared,
        twasWeights = twasWeights,
        ctx = ctx
    ))
}

# The prepared genes' ctwas weight entries, keyed by gene id, with their weight
# variances from `panels`, the blocks' LD panels.
# @noRd
.ctwasGeneEntries <- function(prepared, panels) {
    keys <- map_chr(prepared, "key")
    .ctwasRequireUniqueGenes(keys)
    .ctwasRequireDistinctGroups(prepared)
    set_names(
        map2(
            prepared,
            .ctwasWeightVariances(panels, prepared),
            .ctwasGeneEntry
        ),
        keys
    )
}

# One unit per weight vector a source holds: each (row, method) entry of a
# TwasWeights, each row of a fine-mapping result.
# @noRd
.ctwasWeightUnits <- function(x) {
    if (methods::is(x, "TwasWeights")) {
        cells <- .twasCells(x)
        return(map2(cells$row, cells$method, .ctwasUnit))
    }
    map2(seq_len(nrow(x)), as.character(x$method), .ctwasUnit)
}

# @noRd
.ctwasUnit <- function(row, method) {
    list(row = row, method = method)
}

# The payload of one unit: a TwasWeights entry, or a fine-mapping row.
# @noRd
.ctwasUnitParts <- function(x, unit) {
    if (methods::is(x, "TwasWeights")) {
        return(.twrCellParts(x, unit$row, unit$method))
    }
    .collectionEntry(x, unit$row)
}

# Two weight vectors for one (study, context, type, trait, method) -- rows that
# differ only in their joint group -- would both be modeled under one gene
# id. Refuse, rather than keep whichever came last.
# @noRd
.ctwasRequireUniqueGenes <- function(keys) {
    dup <- unique(keys[duplicated(keys)])
    if (length(dup) == 0L) {
        return(invisible(NULL))
    }
    abort(glue(
        "ctwasPipeline: the weights hold more than one fit of ",
        "{str_flatten(dup, ', ')} (rows that differ only in their joint ",
        "group, or one gene placed in several blocks). cTWAS models one ",
        "weight vector per gene, so keep one."
    ))
}

# cTWAS tells groups apart by their names, so distinct (study, context, type)
# must not read the same: the underscores of a result's label cannot be split
# back. And its EM finds a group's context and type by splitting the name on the
# "|" between them (.ctwasEntryGroup()), so none of the three may hold one.
# @noRd
.ctwasRequireDistinctGroups <- function(prepared) {
    metas <- map(prepared, "meta")
    groups <- distinct(tibble(
        study = map_chr(metas, "study"),
        context = map_chr(metas, "context"),
        type = map_chr(metas, "type")
    ))
    labels <- str_c(groups$study, "_", groups$context, "_", groups$type)
    repeated <- unique(labels[duplicated(labels)])
    if (length(repeated) > 0L) {
        abort(glue(
            "ctwasPipeline: distinct studies, contexts and data types read ",
            "the same once joined: {str_flatten(repeated, ', ')}. cTWAS would ",
            "pool them or results could not tell them apart, so rename one."
        ))
    }
    piped <- labels[str_detect(labels, fixed("|"))]
    if (length(piped) > 0L) {
        abort(glue(
            "ctwasPipeline: a study, context or data type holds a '|': ",
            "{str_flatten(piped, ', ')}. cTWAS splits its group names on it, ",
            "so rename it."
        ))
    }
    invisible(NULL)
}

# Panel variant info in the (chrom/pos/A2/A1/variant_id) frame harmonizeAlleles
# expects (A2 = ref, A1 = alt).
# @noRd
.ctwasRefVariants <- function(panelInfo) {
    tibble(
        chrom = as.integer(panelInfo$chrom),
        pos = as.integer(panelInfo$pos),
        A2 = as.character(panelInfo$ref),
        A1 = as.character(panelInfo$alt),
        variant_id = as.character(panelInfo$id)
    )
}

# The molecular type of a gene whose weights record none: every such context is
# one type.
# @noRd
.CTWAS_NO_TYPE <- "QTL"

# One unit's identity (.ctwasWeightUnits()) and the names cTWAS knows it by: its
# context, keeping two studies' contexts of one name apart as cTWAS's own FAQ
# suggests, its weight_name, as cTWAS's default forms it but for the study and
# the method, and its id, whose first "|" field cTWAS reads as the molecular id.
# @noRd
.ctwasGeneMeta <- function(twasWeights, unit) {
    i <- unit$row
    studyName <- as.character(twasWeights$study)[[i]]
    context <- as.character(twasWeights$context)[[i]]
    trait <- as.character(twasWeights$trait)[[i]]
    method <- unit$method
    type <- .ctwasUnitType(twasWeights, unit)
    engineContext <- str_c(studyName, "_", context)
    weightName <- str_c(engineContext, "_", type, "_", method)
    list(
        study = studyName,
        context = context,
        trait = trait,
        method = method,
        type = type,
        engineContext = engineContext,
        weightName = weightName,
        traitPos = .ctwasTraitPosAt(twasWeights, i),
        key = str_c(trait, "|", weightName)
    )
}

# The molecular type one unit's weights record, or .CTWAS_NO_TYPE.
# @noRd
.ctwasUnitType <- function(x, unit) {
    recorded <- .ctwasRecordedType(x, unit)
    if (length(recorded) == 0L || is.na(recorded)) {
        return(.CTWAS_NO_TYPE)
    }
    as.character(recorded)
}

# @noRd
.ctwasRecordedType <- function(x, unit) {
    md <- mcols(x, use.names = FALSE)
    if (methods::is(x, "TwasWeights")) {
        return(md$dataType[[unit$row]][[unit$method]])
    }
    md$dataType[unit$row]
}

# Steps 1-2: allele-harmonize a gene's (variantIds, weights) against the LD
# panel (matching by chrom/pos with sign/strand-flip detection; weights
# sign-flipped for swapped frames), then restrict to the panel variant set.
# Returns list(vids, w, keptIdx, origVids, origW), or NULL when nothing
# survives.
# @noRd
.ctwasAlignGeneWeights <- function(parts, refVariants, panelSnps) {
    wr <- .rowResolveWeights(parts)
    single <- .singleWeightColumn(wr$weights)
    if (is.null(single)) {
        inform(
            str_c(
                "ctwasPipeline: weights with several columns (one per ",
                "feature or outcome) give no single weight per variant for ",
                "cTWAS to model, so those genes are skipped."
            ),
            .frequency = "once",
            .frequency_id = "pecotmrCtwasMatrixWeights"
        )
        return(NULL)
    }
    wr$weights <- single
    if (length(wr$variantIds) == 0L) {
        return(NULL)
    }
    harm <- .ctwasHarmonizeWeights(wr$variantIds, wr$weights, refVariants)
    if (is.null(harm) || nrow(harm) == 0L) {
        return(NULL)
    }
    vids <- as.character(harm$variant_id)
    w <- as.numeric(harm$w)
    keptIdx <- as.integer(harm$origIdx)
    keep <- is_in(vids, panelSnps)
    if (!any(keep)) {
        return(NULL)
    }
    list(
        vids = vids[keep],
        w = w[keep],
        keptIdx = keptIdx[keep],
        origVids = wr$variantIds,
        origW = wr$weights
    )
}

# The fit a gene's weights are renormalized with: the TWAS row's own, else --
# a fine-mapping method's rows keep none (.twasRowRetention()) -- that
# method's fit for the gene in `fmr`, the fine-mapping result the run was
# given, cut to the gene's own condition when it is one of a joint fit's
# (.ctwasConditionFit()). A fine-mapping row as the weight source has none
# by design (.rowFits()).
# @noRd
.ctwasRenormFit <- function(parts, fmr, meta) {
    own <- .rowFits(parts)
    token <- .twasFmTokenFor(meta$method)
    if (
        !is.null(own) || is.na(token) || !methods::is(parts, "TwasWeightsRow")
    ) {
        return(own)
    }
    if (is.null(fmr)) {
        inform(
            str_c(
                "ctwasPipeline: a fine-mapping method's weights are ",
                "renormalized over the variants kept from that method's fit ",
                "in `fineMappingResult`; none was given, so they are used ",
                "as they are."
            ),
            .frequency = "once",
            .frequency_id = "pecotmrCtwasRenormFit"
        )
        return(NULL)
    }
    idx <- .matchTupleRows(
        fmr,
        list(study = meta$study, context = meta$context, trait = meta$trait)
    )
    fmTokens <- map_chr(as.character(fmr$method)[idx], .twasFmTokenFor)
    hit <- idx[is_in(fmTokens, token)]
    if (length(hit) == 0L) {
        return(NULL)
    }
    .ctwasConditionFit(
        .fmrPartsSusieFit(.fmrRowParts(fmr, hit[[1L]])),
        .fmrJointIndexAt(fmr, hit[[1L]])
    )
}

# A joint fit cut to one of its conditions: an mvSuSiE mu is L x p x R, and
# the row's position `r` in the fit's condition axis picks its own L x p. A
# fit with no condition axis, or a row recording no position, is kept whole.
# @noRd
.ctwasConditionFit <- function(fit, r) {
    mu <- fit[["mu"]]
    if (length(dim(mu)) <= 2L || is.na(r)) {
        return(fit)
    }
    d <- dim(mu)
    own <- array(mu[,, r], dim = d[1:2], dimnames = dimnames(mu)[1:2])
    list_assign(fit, mu = own)
}

# Row `i`'s position in its joint fit's condition axis, NA for a row from no
# joint fit (or one that records none).
# @noRd
.fmrJointIndexAt <- function(x, i) {
    jointIndex <- mcols(x, use.names = FALSE)$jointIndex
    if (is.null(jointIndex)) NA_integer_ else jointIndex[[i]]
}

# Steps 3-4: SuSiE alpha renormalization (when the kept variant set shrank the
# fit) and variance scaling for non-standardized weights (w * sqrt(per-variant
# genotype variance from the LD panel)). `fits` is what .ctwasRenormFit()
# found. Returns the adjusted weight vector.
# @noRd
.ctwasAdjustGeneWeights <- function(parts, aligned, spanVariance, fits) {
    shrank <- length(aligned$keptIdx) < length(aligned$origVids)
    # A NULL renormalization means the fit could not be re-keyed onto the
    # kept variants, so the harmonized weights stand as they are.
    renorm <- if (.ctwasIsSusieFit(fits) && shrank) {
        .ctwasRenormalizeSusieWeights(
            fits,
            origVids = aligned$origVids,
            origW = aligned$origW,
            keptIdx = aligned$keptIdx,
            harmonizedW = aligned$w
        )
    }
    renormalized <- renorm %||% aligned$w
    if (.rowStandardized(parts)) {
        return(renormalized)
    }
    varLookup <- spanVariance[aligned$vids]
    if (anyNA(varLookup)) {
        msg <- glue(
            "ctwasPipeline: missing genotype variance for ",
            "{sum(is.na(varLookup))} variant(s) in the LD panel; an LdMatrix ",
            "needs an `AF` column to scale non-standardized weights."
        )
        abort(msg)
    }
    renormalized * sqrt(varLookup)
}

# One prepared gene's ctwas weight entry (weight matrix, weight variance w'Rw,
# chrom/BP span, plus identity metadata). ctwas would read the weights' LD
# submatrix (R_wgt) only to compute the gene z (compute_gene_z); the pipeline
# supplies that z itself (.ctwasGeneZ()), so the entry carries the scalar
# w'Rw instead. The weights are all the gene's, wherever they fall: ctwas cuts
# them to a region's variants itself when it fine-maps (compute_region_cor), and
# reads the span to find the genes that cross a boundary (get_boundary_genes).
# @noRd
.ctwasGeneEntry <- function(g, wgtVar) {
    meta <- g$meta
    list(
        wgt = matrix(g$w, ncol = 1L, dimnames = list(g$vids, "wgt")),
        wgt_var = wgtVar,
        type = meta$type,
        context = meta$engineContext,
        gene_name = meta$trait,
        study = meta$study,
        qtl_context = meta$context,
        method = meta$method,
        n_wgt = length(g$vids),
        chrom = g$span$chrom,
        p0 = g$span$p0,
        p1 = g$span$p1,
        molecular_id = meta$trait,
        weight_name = meta$weightName,
        # The GENE's own coordinate, carried alongside the weight-variant span
        # (chrom/p0/p1) so the result payloads can be ranged on the gene rather
        # than on wherever its weight variants happen to fall.
        trait_pos = meta$traitPos
    )
}

# The chromosome and span of a gene's weight variants, from the whole-panel
# table: ctwas fine-maps one region at a time, but a gene can reach past it.
# @noRd
.ctwasGeneSpan <- function(vids, spanInfo) {
    rowIdx <- match(vids, spanInfo$id)
    list(
        chrom = as.integer(spanInfo$chrom[[rowIdx[1L]]]),
        p0 = min(as.integer(spanInfo$pos[rowIdx])),
        p1 = max(as.integer(spanInfo$pos[rowIdx]))
    )
}

# cTWAS's name for a gene's group, from the labels its entry carries.
# @noRd
.ctwasEntryGroup <- function(entry) {
    str_c(entry$context, "|", entry$type)
}

# The gene's own position for row i, or NULL when the weight source carries no
# traitPos provenance.
# @noRd
.ctwasTraitPosAt <- function(twasWeights, i) {
    if (!is_in("traitPos", .tupleColumnNames(twasWeights))) {
        return(NULL)
    }
    tp <- traitPosition(twasWeights)
    if (!methods::is(tp, "GRanges") || length(tp) < i) {
        return(NULL)
    }
    tp[i]
}

# Build one gene's ctwas weight record: align -> adjust -> smart-filter
# (PIP/CS + magnitude + cap) -> entry. Returns list(key, entry), or NULL when
# the gene contributes no usable variants.
# @noRd
.ctwasGenePrepared <- function(unit, twasWeights, ctx) {
    # Polymorphic: the weight source may be a TwasWeights or a
    # fine-mapping result, so the payload comes from the class-aware bridge.
    parts <- .ctwasUnitParts(twasWeights, unit)
    near <- .ctwasPanelNear(parts, ctx)
    aligned <- .ctwasAlignGeneWeights(parts, near$refVariants, near$panelSnps)
    if (is.null(aligned)) {
        return(NULL)
    }
    meta <- .ctwasGeneMeta(twasWeights, unit)
    w <- .ctwasAdjustGeneWeights(
        parts,
        aligned,
        near$spanVariance,
        .ctwasRenormFit(parts, ctx$fineMappingResult, meta)
    )
    kept <- .ctwasGeneKept(aligned$vids, w, meta, ctx)
    if (length(kept) < 1L) {
        return(NULL)
    }
    # No LD yet: the weight variances come in a second sweep, from one read
    # of each block's LD over every gene's variants in it together.
    list(
        key = meta$key,
        vids = kept$vids,
        w = kept$w,
        span = .ctwasGeneSpan(kept$vids, near$spanInfo),
        meta = meta
    )
}

# The gene's variants `vids` that pass the smart filter -- fine-mapping PIP
# and credible-set membership, weight magnitude, the cap -- with their
# weights (.ctwasFilterVariants()).
# @noRd
.ctwasGeneKept <- function(vids, w, meta, ctx) {
    finemapAux <- .ctwasGetFinemapAux(
        ctx$fineMappingResult,
        meta$study,
        meta$context,
        meta$trait
    )
    .ctwasFilterVariants(
        vids = vids,
        w = w,
        finemapAux = finemapAux,
        twasWeightCutoff = ctx$cutoffs$twasWeightCutoff,
        csMinCor = ctx$cutoffs$csMinCor,
        minPipCutoff = ctx$cutoffs$minPipCutoff,
        maxNumVariants = ctx$cutoffs$maxNumVariants
    )
}

# Each prepared gene's weight variance w'Rw. As in ctwas's own weights' LD
# (compute_weight_LD_from_ref), the blocks a gene crosses count as independent:
# the variance is the sum of its pieces' w'Rw, each from its own block's LD,
# read once over that block's variants of all the genes. ctwas reads a gene's
# weight LD only for this quadratic form (compute_gene_z), so the pipeline
# computes the gene z itself (.ctwasGeneZ()) and no per-gene or block-wide LD
# matrix is ever formed for it.
# @noRd
.ctwasWeightVariances <- function(panels, prepared) {
    if (length(prepared) == 0L) {
        return(numeric(0))
    }
    pieces <- .ctwasWeightPieces(panels, prepared)
    byBlock <- split(pieces, pieces$block)
    variances <- list_rbind(map2(
        byBlock,
        panels[as.integer(names(byBlock))],
        .ctwasBlockVariances
    ))
    as.numeric(tapply(variances$variance, variances$gene, sum))
}

# Every weight of every gene, with the panel its variant is in.
# @noRd
.ctwasWeightPieces <- function(panels, prepared) {
    vids <- map(prepared, "vids")
    panelIds <- map(panels, .ctwasPanelVariantIds)
    id <- list_c(vids)
    at <- match(id, list_c(panelIds))
    if (anyNA(at)) {
        abort(glue(
            "ctwasPipeline: {sum(is.na(at))} weight variant(s) absent from ",
            "the blocks' LD panels."
        ))
    }
    tibble(
        gene = rep(seq_along(vids), lengths(vids)),
        id = id,
        w = list_c(map(prepared, "w")),
        block = rep(seq_along(panelIds), lengths(panelIds))[at]
    )
}

# @noRd
.ctwasPanelVariantIds <- function(panel) {
    panel$snpInfo$id
}

# Each gene's w'Rw over its variants in one block, from one read of the block's
# LD. `pieces` are the block's rows of .ctwasWeightPieces().
# @noRd
.ctwasBlockVariances <- function(pieces, panel) {
    ld <- .ctwasPanelWeightLd(panel, unique(pieces$id))
    genes <- split(pieces, pieces$gene)
    tibble(
        gene = as.integer(names(genes)),
        variance = map_dbl(genes, .ctwasPieceVariance, ld = ld)
    )
}

# @noRd
.ctwasPieceVariance <- function(piece, ld) {
    .ctwasWeightVariance(list(vids = piece$id, w = piece$w), ld)
}

# The block's LD over `vids`, for weight variances: list(ldFactor = ) holding a
# genotype reference's factor (.ldRefFactor()), or list(R = ) holding the
# matrix a panel already carries -- a caller who computed it, or a payload it
# travelled with -- or an LdMatrix's.
# @noRd
.ctwasPanelWeightLd <- function(ldPanel, vids) {
    idx <- match(vids, ldPanel$snpInfo$id)
    if (anyNA(idx)) {
        msg <- glue(
            "ctwasPipeline: {sum(is.na(idx))} weight variant(s) absent from ",
            "the block's LD panel."
        )
        abort(msg)
    }
    if (!is.null(ldPanel$R)) {
        return(list(R = ldPanel$R[vids, vids, drop = FALSE]))
    }
    if (.ldIsGenotypeRef(ldPanel$sketch)) {
        factor <- .ldRefFactor(ldPanel$sketch, idx)
        return(list(ldFactor = `colnames<-`(factor, vids)))
    }
    list(R = .ctwasPanelLdFor(ldPanel, idx, vids))
}

# One gene's weight variance w'Rw from .ctwasPanelWeightLd()'s LD.
# @noRd
.ctwasWeightVariance <- function(g, ld) {
    if (!is.null(ld$ldFactor)) {
        f <- ld$ldFactor[, g$vids, drop = FALSE]
        return(sum((f %*% g$w)^2) / (nrow(f) - 1))
    }
    r <- ld$R[g$vids, g$vids, drop = FALSE]
    as.numeric(crossprod(g$w, r %*% g$w))
}

# The gene z-scores ctwas::compute_gene_z() would compute,
# w'z / sqrt(w'Rw), from each entry's weights and the weight variance the
# pipeline already holds (`wgt_var`), with ctwas's columns and its rule of
# dropping a gene whose z is infinite.
# @noRd
.ctwasGeneZ <- function(weightsList, zSnp) {
    if (length(weightsList) == 0L) {
        return(NULL)
    }
    zGene <- list_rbind(pmap(
        list(
            weightsList,
            names(weightsList),
            .ctwasGeneSnpZ(weightsList, zSnp)
        ),
        .ctwasGeneZRow
    ))
    infinite <- is.infinite(zGene$z)
    if (any(infinite)) {
        warn(glue(
            "ctwasPipeline: dropped {sum(infinite)} gene(s) with an infinite ",
            "z-score."
        ))
    }
    as.data.frame(filter(zGene, !infinite))
}

# Each gene's weight variants' GWAS z-scores, from one match over all of them.
# @noRd
.ctwasGeneSnpZ <- function(weightsList, zSnp) {
    ids <- map(unname(weightsList), .ctwasWeightIds)
    gene <- rep(seq_along(ids), lengths(ids))
    unname(split(
        zSnp$z[match(list_c(ids), zSnp$id)],
        factor(gene, levels = seq_along(ids))
    ))
}

# @noRd
.ctwasWeightIds <- function(entry) {
    rownames(entry$wgt)
}

# One gene's z row. Every weight variant must carry a GWAS z, as
# compute_gene_z() requires.
# @noRd
.ctwasGeneZRow <- function(entry, id, zs) {
    wgt <- entry$wgt
    if (anyNA(zs)) {
        abort(glue(
            "ctwasPipeline: {sum(is.na(zs))} weight variant(s) of gene ",
            "'{id}' have no GWAS z-score."
        ))
    }
    z <- as.numeric(crossprod(wgt, zs)) / sqrt(entry$wgt_var)
    .ctwasZGeneRow(entry, id, z)
}

# A gene's z_gene row as ctwas reads it: its weight id, its z, and the type,
# context and group its entry carries.
# @noRd
.ctwasZGeneRow <- function(entry, id, z) {
    tibble(
        id = id,
        z = z,
        type = entry$type,
        context = entry$context,
        group = .ctwasEntryGroup(entry)
    )
}

# Look up one gene's PIP vector and 95% credible-set membership / purity in
# the rescues' fine-mapping result, which .ctwasRescueSource() narrowed to the
# one method they read, so the gene's (study, context, trait) picks its fit.
# Returns NULL when no result was passed or the gene has no fit in it. Output
# is a list with:
#   pip       : named numeric vector keyed by variant_id
#   csMembers : list of character vectors (one per CS at 95% coverage)
#   csPurity  : numeric vector aligned with csMembers
# @noRd
.ctwasGetFinemapAux <- function(fmr, studyName, context, trait) {
    if (is.null(fmr)) {
        return(NULL)
    }
    idx <- .matchTupleRows(
        fmr,
        list(study = studyName, context = context, trait = trait)
    )
    if (length(idx) == 0L) {
        return(NULL)
    }
    # The stored table verbatim: every variant, no posterior projection.
    tl <- .fmrPartsTopLoci(.fmrRowParts(fmr, idx[[1L]]))
    if (nrow(tl) == 0L) {
        return(NULL)
    }
    # `pip` is part of the topLoci schema -- FineMappingRow() requires it on
    # any non-empty table, and the zero-row case returned above -- so there is
    # no pip-less frame to guard against here.
    pip <- set_names(as.numeric(tl$pip), as.character(tl$variant_id))
    cs <- .ctwasCsMembership(tl)
    list(pip = pip, csMembers = cs$csMembers, csPurity = cs$csPurity)
}

# One credible set's member variants and its purity.
# @noRd
.ctwasCsEntry <- function(k, tl, csIdx, keepIdx) {
    inCs <- csIdx == k & keepIdx
    purity <- if (is_in("cs_95_purity", names(tl))) {
        as.numeric(tl$cs_95_purity[which(inCs)[1L]])
    } else {
        NA_real_
    }
    list(members = as.character(tl$variant_id)[inCs], purity = purity)
}

# Per-CS membership + purity at 95% coverage from a topLoci table. cs_95 stores
# `<method>_<idx>` where idx == 0 means "not in any CS"; cs_95_purity (when
# present) broadcasts one purity value across a CS's rows.
# @noRd
.ctwasCsMembership <- function(tl) {
    if (!is_in("cs_95", names(tl))) {
        return(list(csMembers = list(), csPurity = numeric(0)))
    }
    csIdx <- .csLabelIndex(tl$cs_95)
    keepIdx <- !is.na(csIdx) & csIdx > 0L
    entries <- map(
        sort(unique(csIdx[keepIdx])),
        .ctwasCsEntry,
        tl = tl,
        csIdx = csIdx,
        keepIdx = keepIdx
    )
    list(
        csMembers = map(entries, "members"),
        csPurity = if (length(entries) == 0L) {
            numeric(0)
        } else {
            map_dbl(entries, "purity")
        }
    )
}

# Apply the four trimCtwasVariants filters to one gene's (vids, w)
# pair. Returns a list(vids, w) with the retained subset, or NULL when
# no variants survive. Filter order:
#   1. Magnitude:   drop variants with |w| < twasWeightCutoff
#   2. CS rescue:   when fineMappingResult is provided, mark variants
#                   in any high-purity CS (purity >= csMinCor) as
#                   "must-keep"
#   3. PIP rescue:  mark variants with PIP > minPipCutoff as must-keep
#   4. Cap:         if surviving variants > maxNumVariants, keep all
#                   must-keep variants and fill remaining slots by
#                   descending PIP (or |w| when no PIP available)
# @noRd
.ctwasFilterVariants <- function(
    vids,
    w,
    finemapAux,
    twasWeightCutoff,
    csMinCor,
    minPipCutoff,
    maxNumVariants
) {
    if (length(vids) == 0L) {
        return(NULL)
    }
    # Step 1: magnitude.
    if (twasWeightCutoff > 0) {
        magKeep <- !is.na(w) & abs(w) >= twasWeightCutoff
        keptVids <- vids[magKeep]
        keptW <- w[magKeep]
        if (length(keptVids) == 0L) {
            return(NULL)
        }
    } else {
        keptVids <- vids
        keptW <- w
    }
    # Steps 2-3: PIP / CS rescue (only when fineMappingResult was passed).
    # Relabel first: the auxiliaries are keyed in the fine-mapping result's
    # frame, the vids in the panel's, and every join from here on is exact.
    finemapAux <- .ctwasRelabelFinemapAux(finemapAux, keptVids)
    mustKeep <- .ctwasMustKeep(keptVids, finemapAux, csMinCor, minPipCutoff)
    # Step 4: cap, keeping must-keep variants first.
    capping <- length(keptVids) > maxNumVariants && is.finite(maxNumVariants)
    capped <- if (!capping) {
        list(vids = keptVids, w = keptW)
    } else {
        .ctwasCapVariants(
            keptVids,
            keptW,
            mustKeep,
            finemapAux,
            maxNumVariants
        )
    }
    list(vids = capped$vids, w = capped$w)
}

# Move a set of variant ids into the frame `vids` uses, leaving ids with no
# match alone.
# @noRd
.ctwasRelabelIds <- function(ids, vids) {
    if (length(ids) == 0L) {
        return(ids)
    }
    m <- .ctwasMatchToPanel(ids, vids)
    if (length(m$idxA) == 0L) {
        return(ids)
    }
    replace(ids, m$idxA, vids[m$idxB])
}

# Relabel a gene's fine-mapping auxiliaries into the weight/panel frame.
#
# `pip` and `csMembers` are keyed by the fine-mapping result's own spelling of
# a variant, while `vids` arrive harmonized to the LD panel. Where the two
# spell one variant with its alleles swapped, the exact-string joins below --
# `intersect(csMembers, vids)` and `pip[vids]` -- miss it: the credible-set
# member silently loses its must-keep protection, and its PIP silently reads
# NA so the cap falls back to |w|. The net effect is dropping the variant the
# fine-mapping was most confident about, with nothing reported.
#
# Only the LABEL moves here. A PIP and a credible-set membership do not depend
# on which allele is counted, so unlike the z-scores there is no sign to apply.
# @noRd
.ctwasRelabelFinemapAux <- function(finemapAux, vids) {
    if (is.null(finemapAux)) {
        return(NULL)
    }
    relabelled <- if (is.null(finemapAux$pip)) {
        finemapAux
    } else {
        list_assign(
            finemapAux,
            pip = set_names(
                finemapAux$pip,
                .ctwasRelabelIds(names(finemapAux$pip), vids)
            )
        )
    }
    if (length(finemapAux$csMembers) == 0L) {
        return(relabelled)
    }
    list_assign(
        relabelled,
        csMembers = map(finemapAux$csMembers, .ctwasRelabelIds, vids = vids)
    )
}

# Variants rescued by credible set `k`, when that set is pure enough.
# @noRd
.ctwasCsRescued <- function(k, finemapAux, vids, csMinCor) {
    if (is.na(finemapAux$csPurity[k]) || finemapAux$csPurity[k] < csMinCor) {
        return(character(0))
    }
    intersect(finemapAux$csMembers[[k]], vids)
}

# Variants that must survive the cap: members of any high-purity (>= csMinCor)
# credible set, plus any with PIP > minPipCutoff. Empty when no finemapAux.
# @noRd
.ctwasMustKeep <- function(vids, finemapAux, csMinCor, minPipCutoff) {
    if (is.null(finemapAux)) {
        return(character(0))
    }
    fromCs <- if (length(finemapAux$csMembers) > 0L && csMinCor > 0) {
        map(
            seq_along(finemapAux$csMembers),
            .ctwasCsRescued,
            finemapAux = finemapAux,
            vids = vids,
            csMinCor = csMinCor
        )
    } else {
        list()
    }
    fromPip <- if (!is.null(finemapAux$pip) && minPipCutoff > 0) {
        hits <- names(finemapAux$pip)[finemapAux$pip > minPipCutoff]
        intersect(hits, vids)
    } else {
        character(0)
    }
    # union() folded over the pieces, so each variant appears once and in
    # first-rescued order -- the same thing the repeated unions produced.
    reduce(c(fromCs, list(fromPip)), union, .init = character(0))
}

# Cap to maxNumVariants: must-keep variants first, then fill by descending PIP
# (falling back to |w| for variants the PIP table doesn't cover).
# @noRd
.ctwasCapVariants <- function(vids, w, mustKeep, finemapAux, maxNumVariants) {
    fromPip <- if (!is.null(finemapAux) && !is.null(finemapAux$pip)) {
        unname(finemapAux$pip[vids])
    } else {
        NULL
    }
    # No usable PIPs at all -> rank on |weight|; otherwise fill only the gaps.
    priorities <- if (is.null(fromPip) || all(is.na(fromPip))) {
        abs(w)
    } else {
        replace(fromPip, is.na(fromPip), abs(w)[is.na(fromPip)])
    }
    isMust <- is_in(vids, mustKeep)
    ord <- order(!isMust, -priorities)
    keepIdx <- ord[seq_len(min(maxNumVariants, length(vids)))]
    list(vids = vids[keepIdx], w = w[keepIdx])
}

# z_gene from a supplied TWAS-Z GRanges (output of causalInferencePipeline):
# each weighted gene's z, from the TWAS row of its (study, context, trait,
# method) for the run's GWAS study, under the gene's weight id with the
# columns ctwas reads (.ctwasZGeneRow()). ctwas keeps only the genes whose
# weight id z_gene names, so a weighted gene with no z is refused rather than
# dropped. Rows for other genes, methods or GWAS studies, and MR rows, are
# not this run's.
# @noRd
.ctwasBuildZGene <- function(twasZ, weightsList, gwasStudy) {
    if (length(weightsList) == 0L) {
        return(NULL)
    }
    rows <- .ctwasTwasZForStudy(
        .ctwasTwasRowsOnly(as_tibble(as.data.frame(S4Vectors::mcols(twasZ)))),
        gwasStudy
    )
    keys <- str_c(
        rows$qtlStudy,
        "|",
        rows$context,
        "|",
        rows$trait,
        "|",
        rows$method
    )
    geneKeys <- map_chr(unname(weightsList), .ctwasEntryKey)
    .ctwasRequireOneZ(names(weightsList), geneKeys, keys)
    zs <- as.numeric(rows$twasZ)[match(geneKeys, keys)]
    # Base data.frame (not tibble): z_gene is indexed positionally by the ctwas
    # engine (df[, "col"] -> vector); a tibble breaks it. See
    # .ctwasSnpInfoForBlock.
    as.data.frame(list_rbind(pmap(
        list(weightsList, names(weightsList), zs),
        .ctwasZGeneRow
    )))
}

# A gene's (study, context, trait, method) key, as twasZ identifies it.
# @noRd
.ctwasEntryKey <- function(entry) {
    str_c(
        entry$study,
        "|",
        entry$qtl_context,
        "|",
        entry$gene_name,
        "|",
        entry$method
    )
}

# The TWAS rows of a twasZ table: its MR rows carry no Z, and may share a
# TWAS row's key. A table with no `analysis` column holds TWAS rows only.
# @noRd
.ctwasTwasRowsOnly <- function(rows) {
    if (!is_in("analysis", names(rows))) {
        return(rows)
    }
    filter(rows, .data$analysis == "twas")
}

# The twasZ rows for the run's GWAS study `gwasStudyName` -- named apart from
# the rows' gwasStudy column, which the filter reads. One with no gwasStudy
# column, or a run whose GWAS names no study but whose twasZ holds just one,
# is used whole.
# @noRd
.ctwasTwasZForStudy <- function(rows, gwasStudyName) {
    if (!is_in("gwasStudy", names(rows))) {
        return(rows)
    }
    studies <- unique(as.character(rows$gwasStudy))
    if (is.na(gwasStudyName) && length(studies) > 1L) {
        abort(glue(
            "assembleCtwasInputs: `twasZ` holds several GWAS studies ",
            "({str_flatten(studies, ', ')}) and `gwasSumStats` names none, ",
            "so which z-scores apply is ambiguous."
        ))
    }
    if (is.na(gwasStudyName)) {
        return(rows)
    }
    filter(rows, .data$gwasStudy == gwasStudyName)
}

# Every weighted gene must have exactly one twasZ z.
# @noRd
.ctwasRequireOneZ <- function(ids, geneKeys, keys) {
    missing <- ids[!is_in(geneKeys, keys)]
    if (length(missing) > 0L) {
        abort(glue(
            "assembleCtwasInputs: `twasZ` has no z-score for ",
            "{length(missing)} weighted gene(s), e.g. ",
            "{str_flatten(utils::head(missing, 3L), ', ')}. ctwas would drop ",
            "them; supply their TWAS z-scores or leave `twasZ` NULL."
        ))
    }
    repeated <- unique(keys[duplicated(keys) & is_in(keys, geneKeys)])
    if (length(repeated) > 0L) {
        abort(glue(
            "assembleCtwasInputs: `twasZ` has more than one z-score for ",
            "{str_flatten(utils::head(repeated, 3L), ', ')}."
        ))
    }
    invisible(NULL)
}

# A cTWAS payload's `LD_file` token does double duty: ctwas asserts the file
# exists, and pecotmr dispatches on the same string into the cached per-region
# LD panels. Each token is a sentinel in the session's tempdir
# (.ctwasMintLdTokens()), so a payload that was serialised (the bundled
# ctwas*Example objects) carries tokens a later session no longer has; they
# are re-minted on the way in.
#
# The panels themselves travel with the closures and are never re-read, so
# this only has to keep the token and the cache keys agreeing with each other.
# @noRd
.ctwasResolveLdPaths <- function(payload) {
    ldMap <- payload$LD_map
    # Anything that is not a table with an LD_file column is passed through:
    # callers hand the granular steps hand-built payloads whose LD_map may be
    # a stub, and there is nothing to re-point in one.
    if (!is.list(ldMap) || is.null(ldMap[["LD_file"]])) {
        return(payload)
    }
    stored <- as.character(ldMap[["LD_file"]])
    resolved <- .ctwasRemintMissingTokens(
        stored,
        .ctwasLdMapRegionIds(ldMap, length(stored))
    )
    if (identical(resolved, stored)) {
        return(payload)
    }
    keyMap <- set_names(resolved, stored)
    rekeyed <- list_assign(
        payload,
        LD_map = list_assign(
            payload$LD_map,
            LD_file = unname(resolved),
            SNP_file = unname(resolved)
        )
    )
    .ctwasRekeyLdLoaders(rekeyed, keyMap)
}

# The region ids an LD_map carries, falling back to positions when the table
# has no region_id column (the granular steps accept hand-built stubs).
# @noRd
.ctwasLdMapRegionIds <- function(ldMap, n) {
    rid <- ldMap[["region_id"]]
    if (is.null(rid) || length(rid) != n) {
        return(as.character(seq_len(n)))
    }
    as.character(rid)
}

# Rebuild the loader closures over a cache re-keyed to the resolved tokens,
# so `LD_loader_fun(LD_file)` still finds its panel after the rename.
# @noRd
.ctwasRekeyLdLoaders <- function(payload, keyMap) {
    panels <- .ctwasCachedPanels(payload)
    if (is.null(panels)) {
        return(payload)
    }
    hit <- is_in(names(panels), names(keyMap))
    rekeyed <- `names<-`(
        panels,
        replace(names(panels), hit, unname(keyMap[names(panels)[hit]]))
    )
    list_assign(
        payload,
        LD_loader_fun = .ctwasMultiBlockLdLoader(rekeyed),
        snpinfo_loader_fun = .ctwasMultiBlockSnpInfoLoader(rekeyed)
    )
}

# The panel cache a payload's loader carries, as a plain list, or NULL when
# the payload carries loaders this package did not create. A payload saved
# before the cache became a SimpleList holds a plain list, which as.list()
# returns as it is.
# @noRd
.ctwasCachedPanels <- function(payload) {
    loader <- payload$LD_loader_fun
    if (!is.function(loader)) {
        return(NULL)
    }
    panels <- attr(loader, "ldPanelsByRegion")
    if (is.null(panels)) {
        return(NULL)
    }
    as.list(panels)
}

# Multi-block LD loader for ctwas. ctwas invokes
# `LD_loader_fun(LD_file)` per region during region_data assembly and
# fine-mapping; we dispatch by `LD_file` (the same string set on
# `LD_map$LD_file`) into the cached per-sketch ldPanel.
#
# ctwas passes the loader that one argument and nothing else, so the panel
# cache travels on the function object itself, as its `ldPanelsByRegion`
# attribute: .ctwasLoadRegionLd() reads it back through sys.function(), and
# .ctwasCachedPanels() reads the same attribute, so the cache has one home.
# The cache is a SimpleList, which prints as its length and tokens: printed
# as a plain list, it put every region's LD on the screen whenever a payload
# was printed. The function goes without its source references: a payload
# is saved, and a function loaded with them would carry its whole source
# file along.
# @noRd
#' @importFrom utils removeSource
.ctwasMultiBlockLdLoader <- function(ldPanelsByRegion) {
    `attr<-`(
        removeSource(.ctwasLoadRegionLd),
        "ldPanelsByRegion",
        SimpleList(ldPanelsByRegion)
    )
}

# The LD loader ctwas calls: region `LD_file`'s full LD, from the panel cache
# carried on this function's own attribute (see .ctwasMultiBlockLdLoader()).
# Built here rather than up front: this is the only place a region's full LD
# is actually needed, and ctwas reaches it only for regions that survived
# screening.
# @noRd
.ctwasLoadRegionLd <- function(LD_file) {
    .ctwasPanelLd(.ctwasLoaderPanel(sys.function(), LD_file, "LD loader"))
}

# Multi-block SNP-info loader for ctwas. Mirrors the LD loader, but carries
# only each panel's snpInfo: a saved payload would otherwise hold every
# region's LD twice, once on each loader.
# @noRd
.ctwasMultiBlockSnpInfoLoader <- function(ldPanelsByRegion) {
    `attr<-`(
        removeSource(.ctwasLoadRegionSnpInfo),
        "ldPanelsByRegion",
        SimpleList(map(ldPanelsByRegion, .ctwasSnpInfoPanel))
    )
}

# The part of a cached panel the SNP-info loader reads.
# @noRd
.ctwasSnpInfoPanel <- function(panel) {
    list(snpInfo = panel$snpInfo)
}

# The snpInfo loader ctwas calls, reading the same attribute.
# @noRd
.ctwasLoadRegionSnpInfo <- function(LD_file) {
    .ctwasLoaderPanel(sys.function(), LD_file, "snpInfo loader")$snpInfo
}

# The cached panel `loader` carries for `LD_file`; an error naming the token
# when its cache has none.
# @noRd
.ctwasLoaderPanel <- function(loader, LD_file, what) {
    panel <- attr(loader, "ldPanelsByRegion")[[LD_file]]
    if (is.null(panel)) {
        msg <- glue(
            "ctwasPipeline {what}: no cached panel for ",
            "LD_file = '{LD_file}'"
        )
        abort(msg)
    }
    panel
}

# The on-disk file that backs an LD panel's data, checked to exist before the
# panel builder reads dosages out of it. An LdMatrix holds its LD in memory,
# and a panel held in memory its dosages, so neither has a payload to find.
# @noRd
.ctwasLdPanelKey <- function(sketch) {
    if (.ldIsMatrixRef(sketch)) {
        return(invisible(NULL))
    }
    if (.ldIsMixtureRef(sketch)) {
        walk(.gmComponents(sketch), .ctwasLdPanelKey)
        return(invisible(NULL))
    }
    fmt <- .ldRefFormat(sketch)
    if (is.null(fmt)) {
        return(invisible(NULL))
    }
    stem <- .ldRefReadPath(sketch)
    candidates <- switch(
        fmt,
        "plink2" = c(str_c(stem, ".pgen")),
        "plink1" = c(str_c(stem, ".bed")),
        "gds" = c(stem),
        "vcf" = c(stem),
        stem
    )
    hit <- candidates[file.exists(candidates)]
    if (length(hit) == 0L) {
        msg <- glue(
            "ctwasPipeline: could not find the file behind the LD panel ",
            "(format={fmt}, path={stem}). Looked for: ",
            "{str_flatten(candidates, ', ')}"
        )
        abort(msg)
    }
    hit[[1L]]
}

# ctwas asserts `file.exists()` on every `LD_map$LD_file` / `SNP_file` and
# then dispatches `LD_loader_fun(LD_file)` on that same string, so a token
# has to be BOTH an existing path and unique per region. The genotype file
# backing a sketch is shared by every block on a chromosome, which is why
# using it directly handed every block the first one's panel. A per-region
# empty sentinel satisfies both constraints, and nothing ever reads it --
# the panels travel inside the loader closures.
# @noRd
.ctwasRegionLdTokens <- function(regionIds, sketches) {
    # Fail here, with the sketch-specific message, if a block's panel has no
    # readable payload: the panel builder is about to read dosages
    # out of it.
    walk(sketches, .ctwasLdPanelKey)
    .ctwasMintLdTokens(regionIds)
}

# The sentinels themselves. Separate from the validation above so the payload
# resolver can re-mint them without a sketch in hand.
# @noRd
.ctwasMintLdTokens <- function(regionIds) {
    dir <- file.path(tempdir(), "pecotmr-ctwas-ld")
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    tokens <- file.path(dir, str_c(.ctwasTokenSlug(regionIds), ".ld"))
    walk(tokens, .ctwasTouchToken)
    set_names(tokens, regionIds)
}

# Sentinels live in the session tempdir, so a payload serialised in one
# session carries tokens that no longer exist in the next -- and ctwas
# asserts file.exists() on every one. Mint replacements; the caller re-keys
# the cached panels onto them.
# @noRd
.ctwasRemintMissingTokens <- function(tokens, regionIds) {
    gone <- !file.exists(tokens)
    if (!any(gone)) {
        return(tokens)
    }
    replace(tokens, gone, unname(.ctwasMintLdTokens(regionIds[gone])))
}

# A filesystem-safe stem per region. The index prefix keeps two region ids
# that sanitize to the same string apart.
# @noRd
.ctwasTokenSlug <- function(regionIds) {
    str_c(
        seq_along(regionIds),
        "_",
        str_replace_all(as.character(regionIds), "[^A-Za-z0-9._-]+", "_")
    )
}

# @noRd
.ctwasTouchToken <- function(path) {
    if (!file.exists(path)) {
        file.create(path)
    }
    invisible(NULL)
}

# One entry's variant ids, or none when it does not carry a SNP column.
# @noRd
.ctwasEntrySnpIds <- function(i, gwasSumStats) {
    mc <- S4Vectors::mcols(gwasSumStats[[i]])
    if (!is_in("SNP", colnames(mc))) {
        return(character(0))
    }
    as.character(mc$SNP)
}

# Build a per-block snpInfo table restricted to variants present in the
# GwasSumStats entry. Mirrors `.ctwasSnpInfoForBlock` but restricts to
# the block's GWAS variants (intersected against the cached panel) so
# snp_map[[region_id]] is sized to the block, not the whole panel.
# @noRd
.ctwasSnpInfoForGwasBlock <- function(gwasSumStats, panelSnpInfo) {
    blockIds <- unique(.ctwasConcatChr(map(
        seq_len(nrow(gwasSumStats)),
        .ctwasEntrySnpIds,
        gwasSumStats = gwasSumStats
    )))
    if (length(blockIds) == 0L) {
        return(panelSnpInfo[FALSE, , drop = FALSE])
    }
    # Allele-aware, not `is_in()`: see `.ctwasMatchToPanel()`. Rows stay in
    # panel order, which the ctwas engine indexes positionally.
    m <- .ctwasMatchToPanel(panelSnpInfo$id, blockIds)
    panelSnpInfo[sort(m$idxA), , drop = FALSE]
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# TRUE when a pre-bucketed weight-list element is a modelable weight entry.
# @noRd
.ctwasIsWeightEntry <- function(x) {
    methods::is(x, "TwasWeights") ||
        methods::is(x, "QtlFineMappingResult")
}

# The GWAS study id of one block's sumstats record (character, possibly empty).
# @noRd
.ctwasStudyChr <- function(g) {
    as.character(g$study)
}
