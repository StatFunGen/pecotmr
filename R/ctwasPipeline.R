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
#' @param varStructure How the prior variance is shared across groups; one of
#'   \code{"shared_type"} (default), \code{"shared_context"},
#'   \code{"shared_nonSNP"}, \code{"shared_all"} or \code{"independent"}.
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
        "shared_type",
        "shared_context",
        "shared_nonSNP",
        "shared_all",
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
#' @description Whether a high-PIP gene whose cis window straddles an
#'   LD-block boundary has its adjacent regions merged and re-fine-mapped
#'   after each run (\code{\link{mergeCtwasBoundaryRegions}}).
#' @param enabled Logical. Run the merge step. Default \code{FALSE}; the
#'   other three fields are read only when this is \code{TRUE}.
#' @param pipThreshold Numeric. PIP threshold for selecting which boundary
#'   genes to merge. Default \code{0.5}.
#' @param filterCs Logical. Require the boundary gene to be in a credible set
#'   to be selected. Default \code{FALSE}.
#' @param maxSnp Numeric. Per-merged-region SNP cap. Default \code{Inf}.
#' @return A \code{BoundaryMergeParam} object, a \code{\link{MethodParam}}.
#' @examples
#' BoundaryMergeParam(enabled = TRUE, pipThreshold = 0.8)
#' @export
BoundaryMergeParam <- function(
    enabled = FALSE,
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

# cTWAS arguments the pipeline owns as named settings of its own. They are
# real ctwas formals, so the union check would accept them here -- but
# .ctwasInvoke drops any name the pipeline already supplies, so a value set
# through CtwasOptions() was silently discarded. Refuse it and say where the
# setting lives instead.
# @noRd
.ctwasPipelineOwnedArgs <- function() {
    c(
        thin = "CtwasPriorParam(thin =)",
        niter = "CtwasPriorParam(niter =)",
        niter_prefit = "CtwasPriorParam(niterPrefit =)",
        group_prior_var_structure = "CtwasPriorParam(varStructure =)",
        L = "the pipeline's own `L`",
        ncore = "the pipeline's own `numThreads`",
        maxSNP = "BoundaryMergeParam(maxSnp =)"
    )
}

# @noRd
.ctwasRefusePipelineOwned <- function(extra) {
    owned <- .ctwasPipelineOwnedArgs()
    clash <- intersect(names(extra), names(owned))
    if (length(clash) == 0L) {
        return(invisible(NULL))
    }
    where <- str_flatten(
        sprintf("`%s` -> %s", clash, unname(owned[clash])),
        collapse = "; "
    )
    abort(glue(
        "CtwasOptions: {str_flatten(clash, ', ')} ",
        "{if (length(clash) == 1L) 'is' else 'are'} set by the pipeline, ",
        "not through methodArgs -- a value given here is dropped. Use: ",
        "{where}."
    ))
}

# Every bundle ctwasPipeline() takes, checked together. The argument names
# in the messages are the pipeline's own, which drop the `Args` suffix.
# @noRd
.ctwasAssertBundles <- function(
    ctwasPriorArgs,
    variantPruningArgs,
    boundaryMergeArgs,
    methodArgs
) {
    .assertMethodParam(ctwasPriorArgs, "CtwasPriorParam", "ctwasPrior")
    .assertMethodParam(
        variantPruningArgs,
        "VariantPruningParam",
        "variantPruning"
    )
    .assertMethodParam(
        boundaryMergeArgs,
        "BoundaryMergeParam",
        "boundaryMerge"
    )
    .assertMethodOptions(methodArgs, "CtwasOptions", "methodArgs")
}

#' @title Causal TWAS Pipeline (cTWAS, multi LD block)
#' @description Pipeline that hands a per-block set of
#'   \code{\link{GwasSumStats}} of GWAS Z-scores together with the matching
#'   per-block per-gene TWAS weights and LD sketches to
#'   \code{ctwas::ctwas_sumstats}, producing per-gene posterior inclusion
#'   probabilities for causal genes. Optionally accepts a precomputed TWAS-Z
#'   \code{GRanges} from \code{\link{causalInferencePipeline}} as the
#'   \code{z_gene} input so the per-gene Z is not recomputed inside ctwas.
#'
#' @section LD block convention: \code{gwasSumStats} is ONE
#'   \code{\link{GwasSumStats}} whose elements are LD blocks, keyed by its
#'   \code{blockId} column -- build it with
#'   \code{loadGwasSumStatsFromManifest(..., ldBlocks = <blocks>)}. Per-block
#'   \code{region_info}, \code{LD_map}, and \code{snp_map} entries are built
#'   automatically from the LD sketch and concatenated before the call to
#'   \code{ctwas::ctwas_sumstats}. A single-block input is rejected: cTWAS's EM
#'   cannot converge on a single region, so callers must supply at least two
#'   blocks.
#'
#' @section LD-sketch compatibility check: Per block:
#'   \code{getLdSketch(twasWeights)} (when non-NULL) must come from the same
#'   reference panel as \code{getLdSketch(gwasSumStats)} --- same samples, same
#'   allele orientation on the shared variants. The two need NOT carry the same
#'   variants; a partial overlap, which is what QC-ing the two sides separately
#'   produces, only warns. A different sample set, a swapped A1/A2 on a shared
#'   variant, or no shared variant at all is a hard error.
#'
#' @param gwasSumStats A \code{\link{GwasSumStats}} whose elements are LD
#'   blocks (at least two), keyed by its \code{blockId} column, with
#'   \code{getQcInfo()} non-empty. Build it with
#'   \code{loadGwasSumStatsFromManifest(..., ldBlocks = <blocks>)} and pass it
#'   through \code{\link{summaryStatsQc}}.
#' @param twasWeights The per-gene weight source. Either (a) a FLAT
#'   \code{\link{TwasWeights}} / \code{QtlFineMappingResult} (or a homogeneous
#'   list of them) carrying \code{region} provenance -- each gene is placed into
#'   its home LD block internally by \code{start(region)} (matching cTWAS's
#'   \code{p0} assignment rule); or (b) a pre-bucketed NAMED LIST keyed by
#'   \code{region_id} (keys a SUBSET of \code{gwasSumStats}'s), used as-is.
#'   Blocks without any TWAS weights still contribute their SNP-level signal to
#'   ctwas's joint group prior estimate (the legacy whole-chromosome pattern
#'   where only a few of many LD blocks carry gene weights). A gene whose cis
#'   span straddles a block boundary is homed by its single anchor; the
#'   cross-block signal is cTWAS's boundary-gene concern
#'   (\code{\link{mergeCtwasBoundaryRegions}}), not placement.
#' @param twasZ Optional \code{GRanges} of TWAS Z-scores (output of
#'   \code{\link{causalInferencePipeline}}). When supplied, the per-(trait,
#'   context) Z is used as the \code{z_gene} input to \code{ctwas_sumstats} so
#'   it is not recomputed.
#' @param fineMappingResult Optional \code{QtlFineMappingResult} or
#'   \code{GwasFineMappingResult} carrying the per-variant PIP and credible-set
#'   membership data used by the CS / PIP rescue filters (\code{csMinCor} and
#'   \code{minPipCutoff}). When \code{NULL} (default) the smart filters are
#'   no-ops; only the magnitude filter (\code{twasWeightCutoff}) and the
#'   per-gene cap (\code{maxNumVariants}, ordered by \code{|weight|}) apply.
#' @param method Optional character (length 1). Picks which TWAS method's
#'   weights to feed into ctwas for each (study, context, trait) gene. When
#'   \code{NULL} (default): use \code{"ensemble"} if that method is present
#'   across the weight sources; otherwise use the sole method when only one is
#'   present; otherwise run \strong{every} method as an independent cTWAS run
#'   (one \code{CtwasResult} row-set per method). Passing the name explicitly
#'   (e.g. \code{"mrash"}) restricts the run to that single method.
#' @param L Integer. Max number of single effects in cTWAS fine-mapping.
#'   Default \code{5}. Pass-through to \code{ctwas::ctwas_sumstats}.
#' @param numThreads Number of cores. Default \code{1}.
#' @param ctwasPriorArgs The EM stage that estimates cTWAS's group priors, built
#'   with \code{\link{CtwasPriorParam}}: \code{thin}, \code{niterPrefit},
#'   \code{niter}, \code{varStructure} and \code{fallbackToPrefit}. Handed
#'   whole to \code{\link{estCtwasGroupPriors}}.
#' @param variantPruningArgs Which variants survive into each gene's weight
#'   matrix, built with \code{\link{VariantPruningParam}}:
#'   \code{weightCutoff}, \code{csMinCor}, \code{pipRescue} and
#'   \code{maxVariants}. Handed whole to
#'   \code{\link{assembleCtwasInputs}}. The \code{csMinCor} /
#'   \code{pipRescue} rescues need a \code{fineMappingResult}.
#' @param boundaryMergeArgs Whether a high-PIP gene straddling an LD-block
#'   boundary has its regions merged and re-fine-mapped, built with
#'   \code{\link{BoundaryMergeParam}}: \code{enabled} (default
#'   \code{FALSE}), \code{pipThreshold}, \code{filterCs} and
#'   \code{maxSnp}.
#' @param keepSnps Logical (length 1). When \code{TRUE}, retain the
#'   context-agnostic SNP background of each run as one extra \code{CtwasResult}
#'   row (\code{study = context = "SNP"}, mirroring cTWAS's own \code{"SNP"}
#'   group) so the full ctwas output is reconstructable from
#'   \code{\link{getFinemap}} / \code{getSusieAlpha}. Default \code{FALSE} --
#'   the SNP rows are the null background and are dropped from the structured
#'   gene-level result.
#' @param methodArgs Additional arguments forwarded to ctwas, built
#'   with \code{\link{CtwasOptions}}. Names are checked against what
#'   the ctwas steps accept between them.
#' @return A \code{\link{CtwasResult}} collection: one row per \code{(gwasStudy,
#'   study, context, method)}. A single-context run is one row per method; a
#'   multi-context (joint) run emits per-context rows sharing the same
#'   \code{jointContexts} set and the jointly-estimated group priors. Each row's
#'   \code{\link{CtwasResultEntry}} payload carries that context's per-gene
#'   fine-mapping posteriors (\code{finemap}), the run's \code{param}, and its
#'   \code{regionInfo}. For the raw \code{ctwas::finemap_regions} list (e.g. to
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
#'   ldSketch = ldStem, region = "chr22:10000000-19000000", ldBlocks = blocks)
#' gwasByRegion <- summaryStatsQc(gss,
#'   panelFilter = PanelFilterParam(mafCutoff = 0.0025))
#' ctwasPipeline(gwasSumStats = gwasByRegion,
#'   twasWeights = list(ctwasWeightsExample),
#'   ctwasPrior = CtwasPriorParam(thin = 1, niterPrefit = 3, niter = 10,
#'     fallbackToPrefit = TRUE),
#'   methodArgs = CtwasOptions(min_group_size = 1, min_p_single_effect = 0))
#' @export
ctwasPipeline <- function(
    gwasSumStats,
    twasWeights,
    twasZ = NULL,
    fineMappingResult = NULL,
    method = NULL,
    L = 5L,
    numThreads = 1L,
    ctwasPriorArgs = CtwasPriorParam(),
    variantPruningArgs = VariantPruningParam(),
    boundaryMergeArgs = BoundaryMergeParam(),
    keepSnps = FALSE,
    methodArgs = CtwasOptions()
) {
    .ctwasAssertBundles(
        ctwasPriorArgs,
        variantPruningArgs,
        boundaryMergeArgs,
        methodArgs
    )
    .ctwasRequireNamedLists(gwasSumStats, twasWeights)
    methods <- .ctwasResolveMethods(twasWeights, method)
    gwasStudy <- .ctwasGwasStudy(gwasSumStats)
    rows <- list_flatten(map(
        methods,
        .ctwasRunMethod,
        gwasSumStats = gwasSumStats,
        twasWeights = twasWeights,
        twasZ = twasZ,
        fineMappingResult = fineMappingResult,
        variantPruningArgs = variantPruningArgs,
        ctwasPriorArgs = ctwasPriorArgs,
        numThreads = numThreads,
        L = L,
        boundaryMergeArgs = boundaryMergeArgs,
        gwasStudy = gwasStudy,
        keepSnps = keepSnps,
        methodArgs = methodArgs
    ))
    if (length(rows) == 0L) {
        msg <- glue(
            "ctwasPipeline: no genes were modeled (the weight sources ",
            "produced no usable gene weights for method(s): ",
            "{str_flatten(methods, ', ')})."
        )
        abort(msg)
    }
    .ctwasRowsToResult(rows)
}

# The cTWAS run proper: estimate the group priors, screen regions against
# them, fine-map what survived, and optionally re-fine-map the merged
# boundary regions. Each stage takes the ctwas option list as one named
# argument, so it is passed rather than spliced in as loose top-level args.
# @noRd
.ctwasRunStages <- function(
    inputs,
    ctwasPriorArgs,
    boundaryMergeArgs,
    L,
    numThreads,
    methodArgs
) {
    est <- estCtwasGroupPriors(
        inputs,
        ctwasPriorArgs = ctwasPriorArgs,
        numThreads = numThreads,
        methodArgs = methodArgs
    )
    screened <- screenCtwasRegions(
        est,
        numThreads = numThreads,
        methodArgs = methodArgs
    )
    finemap <- finemapCtwasRegions(
        screened,
        L = L,
        numThreads = numThreads,
        methodArgs = methodArgs
    )
    if (!isTRUE(boundaryMergeArgs$enabled)) {
        return(finemap)
    }
    .ctwasMaybeMerge(
        finemap,
        boundaryMergeArgs = boundaryMergeArgs,
        L = L,
        numThreads = numThreads,
        methodArgs = methodArgs
    )
}

# One cTWAS run for method `m`: assemble inputs -> estimate params -> screen ->
# fine-map (optionally boundary-merge) -> per-context row-specs. `cfg` bundles
# the ctwasPipeline arguments (incl. the `CtwasOptions` option list).
# @noRd
.ctwasRunMethod <- function(
    m,
    gwasSumStats,
    twasWeights,
    twasZ,
    fineMappingResult,
    variantPruningArgs,
    ctwasPriorArgs,
    numThreads,
    L,
    boundaryMergeArgs,
    gwasStudy,
    keepSnps,
    methodArgs
) {
    inputs <- assembleCtwasInputs(
        gwasSumStats = gwasSumStats,
        twasWeights = twasWeights,
        twasZ = twasZ,
        fineMappingResult = fineMappingResult,
        method = m,
        variantPruningArgs = variantPruningArgs
    )
    merged <- .ctwasRunStages(
        inputs,
        ctwasPriorArgs = ctwasPriorArgs,
        boundaryMergeArgs = boundaryMergeArgs,
        L = L,
        numThreads = numThreads,
        methodArgs = methodArgs
    )
    .ctwasRunToRows(
        merged,
        gwasStudy = gwasStudy,
        method = m,
        keepSnps = keepSnps
    )
}

# Boundary-gene region merging: split a high-PIP straddling gene's adjacent
# regions and re-fine-map. Merge-transparent downstream (keyed by gene id).
# @noRd
.ctwasMaybeMerge <- function(
    finemap,
    boundaryMergeArgs,
    L,
    numThreads,
    methodArgs
) {
    # The bundle's terminal; mergeCtwasBoundaryRegions() takes scalars.
    mergePipThresh <- boundaryMergeArgs$pipThreshold
    mergeFilterCs <- boundaryMergeArgs$filterCs
    mergeMaxSNP <- boundaryMergeArgs$maxSnp
    # methodArgs stays a record all the way down: mergeCtwasBoundaryRegions
    # takes one, and c()-ing it into an argument list would append the S4
    # object as a single element rather than splice its entries.
    mergeCtwasBoundaryRegions(
        finemap,
        pipThresh = mergePipThresh,
        filterCs = mergeFilterCs,
        maxSNP = mergeMaxSNP,
        L = L,
        numThreads = numThreads,
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
#' @return A list with elements \code{z_snp}, \code{z_gene} (NULL when no
#'   \code{twasZ}), \code{weights}, \code{region_info}, \code{snp_map},
#'   \code{LD_map}, \code{LD_loader_fun}, \code{snpinfo_loader_fun}, and
#'   \code{resolvedMethod}.
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
#'   ldSketch = ldStem, region = "chr22:10000000-19000000", ldBlocks = blocks)
#' gwasByRegion <- summaryStatsQc(gss,
#'   panelFilter = PanelFilterParam(mafCutoff = 0.0025))
#' assembleCtwasInputs(gwasSumStats = gwasByRegion,
#'   twasWeights = list(ctwasWeightsExample))
#' @export
assembleCtwasInputs <- function(
    gwasSumStats,
    twasWeights,
    twasZ = NULL,
    fineMappingResult = NULL,
    method = NULL,
    variantPruningArgs = VariantPruningParam()
) {
    .assertMethodParam(
        variantPruningArgs,
        "VariantPruningParam",
        "variantPruning"
    )
    # The bundle's terminal: .ctwasBuildWeights and below take scalars.
    twasWeightCutoff <- variantPruningArgs$weightCutoff
    csMinCor <- variantPruningArgs$csMinCor
    minPipCutoff <- variantPruningArgs$pipRescue
    maxNumVariants <- variantPruningArgs$maxVariants
    .ctwasValidateGwasList(gwasSumStats)
    globalPanelInfo <- .ctwasGlobalPanelInfo(gwasSumStats)
    gwasSumStats <- .ctwasGwasByBlock(gwasSumStats)
    twasWeights <- .ctwasResolveAndValidateWeights(twasWeights, gwasSumStats)
    .ctwasValidateOptional(twasZ, fineMappingResult)
    regionIds <- names(gwasSumStats)
    resolvedMethod <- .ctwasResolveMethod(twasWeights, method)
    fp <- .ctwasFirstPass(
        regionIds,
        gwasSumStats,
        twasWeights,
        globalPanelInfo
    )
    globalGwasSnpIds <- unique(list_c(map(fp$zSnpPieces, "id")))
    cutoffs <- list(
        twasWeightCutoff = twasWeightCutoff,
        csMinCor = csMinCor,
        minPipCutoff = minPipCutoff,
        maxNumVariants = maxNumVariants
    )
    weightsList <- .ctwasSecondPass(
        regionIds,
        twasWeights,
        resolvedMethod,
        fp,
        fineMappingResult,
        cutoffs,
        globalGwasSnpIds
    )
    .ctwasAssembleResult(regionIds, fp, weightsList, twasZ, resolvedMethod)
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
            "`loadGwasSumStatsFromManifest(..., ldBlocks = <blocks>)`."
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
    if (length(getQcInfo(gwasSumStats)) == 0L) {
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
# the two-pass assembly consumes. The collection's own ldSketch rides along on
# each subset, so per-region LD lookups are unchanged.
# @noRd
.ctwasGwasByBlock <- function(gwasSumStats) {
    set_names(
        map(seq_along(gwasSumStats), .ctwasPickBlock, x = gwasSumStats),
        as.character(gwasSumStats$blockId)
    )
}

# One block's slice of the collection. `[` carries the collection-level
# ldSketch through verbatim, and combineGwasSumStats() deliberately UNIONS
# the per-piece panels, so without narrowing here every block references the
# whole chromosome and its LD panel would cover all of it (a
# chr19-wide panel runs to ~155 GB). Narrowed to the block's own variants
# the same way summaryStatsQc narrows the sketch it retains.
# @noRd
.ctwasPickBlock <- function(i, x) {
    block <- x[i]
    sketch <- getLdSketch(block)
    if (is.null(sketch)) {
        return(block)
    }
    methods::initialize(
        block,
        ldSketch = .subsetSketchToIds(sketch, as.list(block))
    )
}

# Resolve a flat weight source into per-region buckets (cTWAS's p0 start-of-
# region rule; a pre-bucketed named list passes through), then validate it is a
# named list of TwasWeights / QtlFineMappingResult with no stray region keys.
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
    if (
        !is.null(fineMappingResult) &&
            !methods::is(fineMappingResult, "FineMappingResultBase")
    ) {
        msg <- glue(
            "`fineMappingResult` must be a FineMappingResultBase ",
            "(QtlFineMappingResult or GwasFineMappingResult) or NULL."
        )
        abort(msg)
    }
}

# First pass: cache LD panels + build z_snp / region_info / snp_map per region.
# We need the union of GWAS variant ids ACROSS all blocks before filtering each
# per-block TwasWeights (a gene's weight variants can straddle adjacent blocks).
# Returns list(ldPanelsByRegion, zSnpPieces, regionInfoPieces, snpMap,
# ldFileByRegion).
# @noRd
# One region's GWAS LD sketch, checked against the region's weights.
# The ctwas-shaped snpInfo for the whole collection's panel, captured BEFORE
# .ctwasGwasByBlock() narrows each block's sketch to its own variants. A
# gene's cis SPAN has to cover every block it reaches or
# ctwas::get_boundary_genes cannot route it to merge_regions, so the span is
# measured against this table rather than against any one block's panel.
#
# NULL when the collection carries no sketch -- .ctwasRegionGwasLd() reports
# that case with its own message, so this must not pre-empt it.
# @noRd
.ctwasGlobalPanelInfo <- function(gwasSumStats) {
    sketch <- getLdSketch(gwasSumStats)
    if (is.null(sketch)) {
        return(NULL)
    }
    .ctwasSnpInfoForBlock(sketch)
}

# @noRd
.ctwasRegionGwasLd <- function(rid, gwasSumStats, twasWeights) {
    gwasLd <- getLdSketch(gwasSumStats[[rid]])
    if (is.null(gwasLd)) {
        msg <- glue(
            "ctwasPipeline: GwasSumStats for region '{rid}' carries no ",
            "ldSketch (ldSketch = NULL); cTWAS requires an LD reference."
        )
        abort(msg)
    }
    tw <- twasWeights[[rid]]
    if (!is.null(tw)) {
        .ctwasRequireMatchingLdSketches(getLdSketch(tw), gwasLd)
    }
    gwasLd
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

# Second pass: build per-block weight lists. The GLOBAL gwasSnpIds bounds each
# gene's cis SPAN, so a gene whose cis-window straddles block boundaries is
# still recognised as one; the per-region snpMap bounds the weight vector that
# is actually fitted, so susie_rss never sees a variant the region has no LD
# for. Weight names are prefixed with the region id.
# @noRd
.ctwasSecondPass <- function(
    regionIds,
    twasWeights,
    resolvedMethod,
    fp,
    fineMappingResult,
    cutoffs,
    globalGwasSnpIds
) {
    perRegion <- map(
        regionIds,
        .ctwasRegionBlockWeights,
        twasWeights = twasWeights,
        fp = fp,
        resolvedMethod = resolvedMethod,
        fineMappingResult = fineMappingResult,
        cutoffs = cutoffs,
        globalGwasSnpIds = globalGwasSnpIds
    )
    .ctwasConcat(compact(perRegion))
}

# Concatenate per-region lists, empty-safe.
# @noRd
.ctwasConcat <- function(pieces) {
    if (length(pieces) == 0L) {
        return(list())
    }
    list_c(pieces)
}

# One region's block weights, region-qualified, or NULL when the region has no
# weights for the resolved method.
# @noRd
.ctwasRegionBlockWeights <- function(
    rid,
    twasWeights,
    fp,
    resolvedMethod,
    fineMappingResult,
    cutoffs,
    globalGwasSnpIds
) {
    tw <- twasWeights[[rid]]
    if (is.null(tw)) {
        return(NULL)
    }
    twMethod <- .ctwasFilterMethod(tw, resolvedMethod)
    if (is.null(twMethod)) {
        return(NULL)
    }
    blockWeights <- .ctwasBuildWeights(
        twMethod,
        fp$ldPanelsByRegion[[fp$ldFileByRegion[[rid]]]],
        fineMappingResult = fineMappingResult,
        twasWeightCutoff = cutoffs$twasWeightCutoff,
        csMinCor = cutoffs$csMinCor,
        minPipCutoff = cutoffs$minPipCutoff,
        maxNumVariants = cutoffs$maxNumVariants,
        gwasSnpIds = globalGwasSnpIds,
        regionSnpIds = fp$snpMap[[rid]]$id,
        globalPanelInfo = fp$globalPanelInfo,
        globalVariance = fp$globalVariance
    )
    if (length(blockWeights) == 0L) {
        return(NULL)
    }
    set_names(blockWeights, str_c(rid, "|", names(blockWeights)))
}

# Concatenate the per-region pieces into the ctwas-shape input list.
# @noRd
.ctwasAssembleResult <- function(
    regionIds,
    fp,
    weightsList,
    twasZ,
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
        z_gene = if (!is.null(twasZ)) .ctwasBuildZGene(twasZ) else NULL,
        weights = weightsList,
        region_info = regionInfo,
        snp_map = fp$snpMap,
        LD_map = ldMap,
        LD_loader_fun = .ctwasMultiBlockLdLoader(fp$ldPanelsByRegion),
        snpinfo_loader_fun = .ctwasMultiBlockSnpInfoLoader(fp$ldPanelsByRegion),
        resolvedMethod = resolvedMethod
    )
}

#' Estimate cTWAS group prior + prior variance
#'
#' @description Step 2 of the three-step \code{\link{ctwasPipeline}}: assembles
#'   \code{region_data} from the inputs and runs \code{ctwas::est_param} (prefit
#'   EM + accurate EM) to estimate the group prior probabilities and prior
#'   variances. Returns the input state plus \code{region_data},
#'   \code{boundary_genes}, \code{z_gene}, and \code{param}.
#'
#' @param inputs A list returned by \code{\link{assembleCtwasInputs}}.
#' @param ctwasPriorArgs The EM settings, built with
#'   \code{\link{CtwasPriorParam}}: \code{thin} (SNP thinning when
#'   assembling region data), \code{niterPrefit} and \code{niter} (prefit /
#'   accurate EM iterations), \code{varStructure}, and
#'   \code{fallbackToPrefit} -- when \code{TRUE}, an accurate EM that fails
#'   on a degenerate input is recovered by re-running only the prefit step
#'   and returning those (typically finite) priors.
#' @param numThreads Number of cores.
#' @param methodArgs Additional arguments forwarded to ctwas, built
#'   with \code{\link{CtwasOptions}}. Names are checked against what
#'   the ctwas steps accept between them.
#' @return The \code{inputs} list augmented with \code{region_data},
#'   \code{boundary_genes}, \code{z_gene}, and \code{param}.
#' @examples
#' data(ctwasInputsExample)
#' estCtwasGroupPriors(ctwasInputsExample,
#'   ctwasPrior = CtwasPriorParam(thin = 1, niterPrefit = 3, niter = 10,
#'     fallbackToPrefit = TRUE),
#'   methodArgs = CtwasOptions(min_group_size = 1, min_p_single_effect = 0))
#' @export
estCtwasGroupPriors <- function(
    inputs,
    ctwasPriorArgs = CtwasPriorParam(),
    numThreads = 1L,
    methodArgs = CtwasOptions()
) {
    .assertMethodParam(ctwasPriorArgs, "CtwasPriorParam", "ctwasPrior")
    # The bundle's terminal: the EM helpers below take plain scalars.
    thin <- ctwasPriorArgs$thin
    niterPrefit <- ctwasPriorArgs$niterPrefit
    niter <- ctwasPriorArgs$niter
    groupPriorVarStructure <- ctwasPriorArgs$varStructure
    fallbackToPrefit <- ctwasPriorArgs$fallbackToPrefit
    .assertMethodOptions(methodArgs, "CtwasOptions", "methodArgs")
    if (!requireNamespace("ctwas", quietly = TRUE)) {
        abort("Package 'ctwas' is required for estCtwasGroupPriors.")
    }
    numThreads <- as.integer(numThreads)
    inputs <- .ctwasResolveLdPaths(inputs)
    zGene <- .ctwasEnsureZGene(inputs, numThreads, methodArgs)
    regionData <- .ctwasAssembleRegionData(
        inputs,
        zGene,
        thin,
        numThreads,
        methodArgs
    )
    boundaryGenes <- .ctwasBoundaryGenes(inputs, numThreads, methodArgs)
    paramRes <- .ctwasEstParamOrFallback(
        regionData,
        niterPrefit,
        niter,
        groupPriorVarStructure,
        numThreads,
        thin,
        fallbackToPrefit,
        methodArgs
    )
    # assemble_region_data does not echo z_gene back, so propagate the
    # precomputed z_gene we passed in (inputs$z_gene is NULL when twasZ was not
    # supplied) so $z_gene resolves to the right entry.
    c(
        list_assign(inputs, z_gene = zGene),
        list(
            region_data = regionData,
            boundary_genes = boundaryGenes,
            param = paramRes
        )
    )
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
    numThreads,
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
            numThreads,
            extra
        ))
    }
    try_fetch(
        .ctwasEstParamAccurate(
            regionData,
            niterPrefit,
            niter,
            groupPriorVarStructure,
            numThreads,
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
                numThreads = numThreads,
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
#' @param L Pass-through.
#' @param numThreads Number of cores.
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
#'   ldSketch = ldStem, region = "chr22:10000000-19000000", ldBlocks = blocks)
#' gwasByRegion <- summaryStatsQc(gss,
#'   panelFilter = PanelFilterParam(mafCutoff = 0.0025))
#' inp <- assembleCtwasInputs(gwasSumStats = gwasByRegion,
#'   twasWeights = list(ctwasWeightsExample))
#' est <- estCtwasGroupPriors(inp,
#'   ctwasPrior = CtwasPriorParam(thin = 1, niterPrefit = 3, niter = 10,
#'     fallbackToPrefit = TRUE),
#'   methodArgs = CtwasOptions(min_group_size = 1, min_p_single_effect = 0))
#' screened <- screenCtwasRegions(est)
#' finemapCtwasRegions(screened, L = 5L)
#' @export
finemapCtwasRegions <- function(
    screenResult,
    L = 5L,
    numThreads = 1L,
    methodArgs = CtwasOptions()
) {
    .assertMethodOptions(methodArgs, "CtwasOptions", "methodArgs")
    if (!requireNamespace("ctwas", quietly = TRUE)) {
        abort("Package 'ctwas' is required for finemapCtwasRegions.")
    }
    screenResult <- .ctwasResolveLdPaths(screenResult)
    fmRes <- .ctwasFinemapOrEmpty(screenResult, L, numThreads, methodArgs)
    # Repair cTWAS's molecular_id mislabel (first-"|" split of our composite
    # id).
    list(
        z_gene = screenResult$z_gene,
        param = screenResult$param,
        finemap_res = .ctwasFixMolecularId(fmRes$finemap_res),
        susie_alpha_res = .ctwasFixMolecularId(fmRes$susie_alpha_res),
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
        snpinfo_loader_fun = screenResult$snpinfo_loader_fun
    )
}

#' Merge boundary cTWAS regions and re-fine-map
#'
#' @description Optional step 4 of the cTWAS pipeline (default-off region
#'   merging). A gene whose cis window straddles an LD-block boundary (a
#'   \code{boundary_genes} member) is split across two regions in the first-pass
#'   fine-mapping. This step selects the high-PIP boundary genes, merges each
#'   one's adjacent regions into a single region, re-runs fine-mapping on the
#'   merged regions, and splices the updated results back into the
#'   \code{\link{finemapCtwasRegions}} output. Thin wrapper over
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
#' @param L Integer. Max number of single effects for the merged-region
#'   re-fine-mapping (LD path only). Default \code{5}.
#' @param numThreads Number of cores. Default \code{1}.
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
#' mergeCtwasBoundaryRegions(ctwasFinemapExample)
#' @export
mergeCtwasBoundaryRegions <- function(
    finemapResult,
    pipThresh = 0.5,
    filterCs = FALSE,
    maxSNP = Inf,
    L = 5L,
    numThreads = 1L,
    methodArgs = CtwasOptions()
) {
    .assertMethodOptions(methodArgs, "CtwasOptions", "methodArgs")
    if (!requireNamespace("ctwas", quietly = TRUE)) {
        abort("Package 'ctwas' is required for mergeCtwasBoundaryRegions.")
    }
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
        numThreads
    )
    fa <- .ctwasMergeDispatch(finemapResult, common, L)
    # Flatten the record before merging: c() on a list and a MethodOptions (a
    # SimpleList) appends the S4 object as ONE unnamed element rather than
    # splicing its entries, and exec() would then hand that object to ctwas
    # as a positional argument.
    user <- as.list(methodArgs)
    callArgs <- c(fa$args, user[setdiff(names(user), names(fa$args))])
    res <- exec(fa$fn, !!!callArgs)
    .ctwasApplyMergeResult(finemapResult, res)
}

# Shared postprocess_region_merging argument list built from a finemap result.
# @noRd
.ctwasMergeCommonArgs <- function(
    finemapResult,
    pipThresh,
    filterCs,
    maxSNP,
    numThreads
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
        ncore = as.integer(numThreads)
    )
}

# Write region-merging outputs back onto the finemap result (only components
# ctwas actually returned).
# @noRd
.ctwasApplyMergeResult <- function(finemapResult, res) {
    list_assign(
        finemapResult,
        finemap_res = res$updated_finemap_res,
        susie_alpha_res = res$updated_susie_alpha_res,
        merge_res = res,
        !!!compact(list(
            region_data = res$updated_region_data,
            region_info = res$updated_region_info,
            LD_map = res$updated_LD_map,
            snp_map = res$updated_snp_map
        ))
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
        formalsFn <- try_fetch(names(formals(fn)), error = function(cnd) NULL)
        usable <- if (is.null(formalsFn)) {
            deduped
        } else {
            deduped[intersect(names(deduped), setdiff(formalsFn, "..."))]
        }
        args <- c(args, usable)
    }
    exec(fn, !!!args)
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
    numThreads,
    extra = list()
) {
    fitEm <- utils::getFromNamespace("fit_EM", "ctwas")
    fitRegionData <- .ctwasPrefitRegionFilter(region_data, extra)
    fitArgs <- list(
        region_data = fitRegionData,
        niter = as.integer(niter),
        group_prior_var_structure = groupPriorVarStructure,
        ncore = as.integer(numThreads)
    )
    prefit <- .ctwasInvoke(fitEm, fitArgs, extra)
    adj <- .ctwasApplyThin(prefit$group_prior, prefit$group_size, thin)
    groupSize <- if (length(adj$groupPrior) > 0L) {
        adj$groupSize[names(adj$groupPrior)]
    } else {
        adj$groupSize
    }
    list(
        group_prior = adj$groupPrior,
        group_prior_var = prefit$group_prior_var,
        group_prior_iters = prefit$group_prior_iters,
        group_prior_var_iters = prefit$group_prior_var_iters,
        group_prior_var_structure = groupPriorVarStructure,
        group_size = groupSize,
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
# alone. A vector without a SNP group passes through unchanged.
# @noRd
.ctwasScaleSnpGroup <- function(groups, factor) {
    if (!is_in("SNP", names(groups))) {
        return(groups)
    }
    replace(groups, "SNP", groups[["SNP"]] * factor)
}

# =============================================================================
# Internal helpers
# =============================================================================

# LD-sketch compatibility check. Thin wrapper over the shared
# `.requireMatchingLdSketches` helper (R/ld.R).
.ctwasRequireMatchingLdSketches <- function(twLd, gwasLd) {
    .requireMatchingLdSketches(twLd, gwasLd, pipelineName = "ctwasPipeline")
}

# Resolve which TWAS method's weights to feed into ctwas given a
# TwasWeights collection that may carry multiple methods per
# (study, context, trait). Rules:
#   - Caller-supplied method (non-NULL, non-empty) wins, provided that
#     method exists in the TwasWeights's `method` column.
#   - Otherwise prefer "ensemble" when present.
#   - Otherwise return the sole method when only one is present.
#   - Otherwise: error.
# @noRd
.ctwasResolveMethod <- function(twasWeightsList, method = NULL) {
    available <- unique(list_c(map(twasWeightsList, .ctwasMethodChr)))
    if (length(available) == 0L) {
        abort("ctwasPipeline: TwasWeights collections have no method entries.")
    }
    if (!is.null(method) && str_length(method) > 0L) {
        if (!is_in(method, available)) {
            msg <- glue(
                "ctwasPipeline: method '{method}' not present in TwasWeights ",
                "(available: {str_flatten(available, ', ')})."
            )
            abort(msg)
        }
        return(method)
    }
    if (is_in("ensemble", available)) {
        return("ensemble")
    }
    if (length(available) == 1L) {
        return(available[[1L]])
    }
    msg <- glue(
        "ctwasPipeline: TwasWeights carries multiple methods ",
        "({str_flatten(available, ', ')}) with no 'ensemble' entry. ",
        "Supply a `method` argument to pick one (e.g. method = \"mrash\")."
    )
    abort(msg)
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
            "`loadGwasSumStatsFromManifest(..., ldBlocks = <blocks>)`."
        )
        abort(msg)
    }
    okTw <- methods::is(twasWeights, "TwasWeights") ||
        methods::is(twasWeights, "QtlFineMappingResult") ||
        is.list(twasWeights)
    if (!okTw) {
        msg <- glue(
            "`twasWeights` must be a TwasWeights / QtlFineMappingResult ",
            "weight source (placed into LD blocks internally by region), or ",
            "a per-region named list keyed by region_id ",
            "(got {class(twasWeights)[[1L]]})."
        )
        abort(msg)
    }
}

# Strip a leading "chr" and case so region_id-derived seqnames ("chr22") and
# phenotype rowRanges seqnames ("22" / "chr22") compare equal for placement.
# @noRd
.ctwasChrKey <- function(x) {
    str_remove(str_to_lower(as.character(x)), "^chr")
}

# TRUE when `tw` is ALREADY a per-region named list (the pre-bucketed contract):
# a NAMED plain list (not an S4 collection) of weight collections. Its keys are
# validated against the block grid downstream (the `extra_tw_keys` check). A
# flat collection or an UNNAMED list is instead treated as a flat weight source
# to place internally by region.
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

# Parse cTWAS block-manifest keys ("chr1_1000_2000" / "chr1:1000-2000") into a
# per-key GRanges.
#
# This is cTWAS's OWN `region_id` concept, not pecotmr's: cTWAS keys its
# region_info, snp_map and per-region data by these strings, so the package has
# to be able to read them even though pecotmr's own collections now key on the
# element range instead (section 4.4). A key that carries no coordinates
# becomes a 0-width chrUn sentinel, which matches no anchor and so yields the
# documented "NA when the anchor falls in no block".
# @noRd
# One block id parsed into chrom/start/end, or the unplaced sentinel when it
# does not parse as a range.
# @noRd
.ctwasBlockCoords <- function(id) {
    g <- try_fetch(
        asGranges(str_replace(
            as.character(id),
            "_([0-9]+)_([0-9]+)$",
            ":\\1-\\2"
        )),
        error = function(cnd) NULL
    )
    if (is.null(g) || length(g) < 1L) {
        return(list(chrom = "chrUn", start = 1L, end = 0L))
    }
    list(
        chrom = as.character(GenomicRanges::seqnames(g))[[1L]],
        start = GenomicRanges::start(g)[[1L]],
        end = GenomicRanges::end(g)[[1L]]
    )
}

.ctwasBlockGrFromIds <- function(ids) {
    coords <- map(ids, .ctwasBlockCoords)
    GenomicRanges::GRanges(
        map_chr(coords, "chrom"),
        IRanges::IRanges(
            start = map_int(coords, "start"),
            end = map_int(coords, "end")
        )
    )
}

# Place each gene (row) of a flat weight source into its home LD block. The
# anchor is start(region) -- matching cTWAS's own `assign_region_data` rule,
# which homes a gene by its p0 (single point) into the block where
# p0 in [start, stop). Returns a region_id per row (NA when the anchor falls in
# no block). A gene whose cis SPAN straddles a boundary is still homed by its
# single anchor here; the cross-block signal is cTWAS's boundary-gene concern
# (get_boundary_genes / postprocess_region_merging), not placement.
# @noRd
.ctwasPlaceByAnchor <- function(region, gwasSumStats) {
    ids <- names(gwasSumStats)
    blockGr <- .ctwasBlockWindows(gwasSumStats, ids)
    aChr <- .ctwasChrKey(as.character(GenomicRanges::seqnames(region)))
    aPos <- GenomicRanges::start(region)
    bChr <- .ctwasChrKey(as.character(GenomicRanges::seqnames(blockGr)))
    bS <- GenomicRanges::start(blockGr)
    bE <- GenomicRanges::end(blockGr)
    map_chr(
        seq_along(aPos),
        .ctwasBlockIdForVariant,
        aPos = aPos,
        aChr = aChr,
        bChr = bChr,
        bS = bS,
        bE = bE,
        ids = ids
    )
}

# The window of each LD block, taken from that block's own GWAS variants.
#
# These used to be parsed out of the region-id strings, which forced ids to
# encode coordinates ("chr1_100_200"). Each block is now an element of a
# GwasSumStats and carries its range directly, so an opaque id ("blockA")
# places just as well. A block with no variants yields an empty window and
# hosts no gene, which is the right answer: there is no GWAS signal there for
# a gene to be tested against.
# @noRd
.ctwasBlockWindows <- function(gwasSumStats, ids) {
    spans <- map(gwasSumStats, .ctwasBlockSpan)
    empty <- map_lgl(spans, .ctwasSpanIsEmpty)
    if (all(empty)) {
        # Nothing to place against; fall back to whatever the ids encode so a
        # coordinate-keyed caller still behaves as before.
        return(.ctwasBlockGrFromIds(ids))
    }
    # unname(): c() on a NAMED list of GRanges builds a list rather than
    # concatenating, and the result then fails seqnames() further down.
    exec(c, !!!unname(spans))
}

# @noRd
.ctwasBlockSpan <- function(gss) {
    # S4 dispatch, not list flattening: `gss` is a GwasSumStats collection and
    # unlist() returns the GRanges that range() below needs.
    variants <- unlist(gss, use.names = FALSE)
    if (length(variants) == 0L) {
        return(GenomicRanges::GRanges())
    }
    # A stored element spans exactly one seqname, so ignoring strand leaves a
    # single range -- the block's window.
    range(variants, ignore.strand = TRUE)
}

# @noRd
.ctwasSpanIsEmpty <- function(g) {
    length(g) == 0L
}

# Bucket a flat weight source into a per-region named list keyed to the block
# grid, homing each gene by start(region). Each per-block sub-collection carries
# that block's GWAS LD sketch (the panel its weights are harmonized against, and
# what the downstream match-check expects).
# @noRd
# The genes homed into region `rid`, carrying that region's LD sketch, or
# NULL when nothing landed there. `slot<-` applied as a function returns a
# copy rather than writing into the subset in place.
# @noRd
#' @importFrom methods slot<-
.ctwasBucketForRegion <- function(rid, combined, home, gwasSumStats) {
    idx <- which(home == rid)
    if (length(idx) == 0L) {
        return(NULL)
    }
    `slot<-`(
        combined[idx, ],
        "ldSketch",
        value = getLdSketch(gwasSumStats[[rid]])
    )
}

.ctwasBucketWeights <- function(weights, gwasSumStats) {
    combined <- .ctwasCombineWeightSources(weights)
    # Placement anchors on the GENE's own position, not on a stored analysis
    # window and not on the span of its weight variants. Two genes at
    # different loci can legitimately share a weight variant set, so the
    # variant span cannot tell them apart; traitPos can (spec 4.4, which keeps
    # gene coordinates in mcols for exactly this).
    if (!is_in("traitPos", .tupleColumnNames(combined))) {
        msg <- glue(
            "assembleCtwasInputs: the weight source carries no `traitPos` ",
            "provenance, which internal LD-block placement requires ",
            "(produced by twasWeightsPipeline / fineMappingPipeline). Supply ",
            "a pre-bucketed per-region named list if placement was done ",
            "upstream."
        )
        abort(msg)
    }
    home <- .ctwasPlaceByAnchor(
        getTraitPosition(combined),
        gwasSumStats
    )
    unplaced <- sum(is.na(home))
    if (unplaced > 0L) {
        msg <- glue(
            "assembleCtwasInputs: {unplaced} gene(s) whose traitPos anchor ",
            "fell in no LD block were dropped."
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
            "Check that the weight `region`s and the gwasSumStats region_id ",
            "keys share a coordinate system (e.g. 'chr22_1_1000000')."
        )
        abort(msg)
    }
    out
}

# Resolve `twasWeights` to a per-region named list: pass a pre-bucketed list
# through, otherwise place a flat weight source by region.
# @noRd
.ctwasResolveWeightBuckets <- function(twasWeights, gwasSumStats) {
    if (.ctwasIsPreBucketed(twasWeights)) {
        return(twasWeights)
    }
    .ctwasBucketWeights(twasWeights, gwasSumStats)
}

# Extract the character vector of method names carried by a weight source
# (NULL-safe: a NULL source contributes no methods).
# @noRd
.ctwasMethodsOf <- function(tw) {
    if (is.null(tw)) NULL else as.character(tw$method)
}

# Resolve the LIST of TWAS methods a `ctwasPipeline` run should iterate over
# (one independent cTWAS run per method -- weights are homogeneous within a
# run). - explicit `method`: exactly that one (validated present). - NULL + an
# "ensemble" method present: just "ensemble" (the pre-combined weight, the
# historical default). - NULL + a single method present: that one. - NULL +
# MULTIPLE methods, no "ensemble": ALL of them (the singular
# `.ctwasResolveMethod` errors here; the pipeline instead fans out).
# @noRd
.ctwasResolveMethods <- function(twasWeightsList, method = NULL) {
    available <- unique(
        if (
            methods::is(twasWeightsList, "TwasWeights") ||
                methods::is(twasWeightsList, "QtlFineMappingResult")
        ) {
            # a flat weight source
            .ctwasMethodsOf(twasWeightsList)
        } else {
            list_c(compact(map(twasWeightsList, .ctwasMethodsOf)))
        }
    ) # a list of them
    if (length(available) == 0L) {
        abort("ctwasPipeline: weight sources carry no method entries.")
    }
    if (!is.null(method) && str_length(method) > 0L) {
        if (!is_in(method, available)) {
            msg <- glue(
                "ctwasPipeline: method '{method}' not present in the weight ",
                "sources (available: {str_flatten(available, ', ')})."
            )
            abort(msg)
        }
        return(method)
    }
    if (is_in("ensemble", available)) {
        return("ensemble")
    }
    available # single -> length-1 (one run); multiple -> iterate over all
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

# Extract field `i` (as character) from each `region|study|context|trait|method`
# split in `parts`.
# @noRd
.ctwasPickField <- function(i, parts) {
    map_chr(parts, i)
}

# Parse the cTWAS gene ids (`region|study|context|trait|method`) that name the
# assembled weights list into their identity components. `method` is the LAST
# field and `trait` everything between context and method, so a trait that
# itself contains "|" is preserved.
# @noRd
.ctwasParseGeneIds <- function(ids) {
    parts <- str_split(ids, "\\|")
    n <- lengths(parts)
    if (any(n < 5L)) {
        msg <- glue(
            "ctwasPipeline: malformed cTWAS gene id(s): ",
            "{str_flatten(ids[n < 5L], ', ')} ",
            "(expected 'region|study|context|trait|method')."
        )
        abort(msg)
    }
    tibble(
        id = ids,
        rid = .ctwasPickField(1L, parts),
        study = .ctwasPickField(2L, parts),
        context = .ctwasPickField(3L, parts),
        trait = map_chr(parts, .ctwasTraitField),
        method = map_chr(parts, .ctwasMethodField)
    )
}

# cTWAS's finemap_regions derives `molecular_id` by splitting the gene id on the
# FIRST "|", which mislabels our composite `region|study|context|trait|method`
# id (it takes the region as the molecular_id). Restore the true trait for gene
# rows; SNP-background rows (variant ids, no "|") are left untouched.
# @noRd
.ctwasFixMolecularId <- function(df) {
    if (
        is.null(df) ||
            !is.data.frame(df) ||
            nrow(df) == 0L ||
            !all(is_in(c("id", "molecular_id"), names(df)))
    ) {
        return(df)
    }
    isGene <- lengths(str_split(as.character(df$id), "\\|")) >= 5L
    if (any(isGene)) {
        return(mutate(
            df,
            molecular_id = replace(
                .data$molecular_id,
                isGene,
                .ctwasParseGeneIds(as.character(df$id)[isGene])$trait
            )
        ))
    }
    df
}

# Enforce the multi-context joint-model invariant: every context in a run must
# carry the SAME set of genes (traits). A cTWAS joint fit couples the contexts
# through shared group priors; contexts with disjoint gene sets make the joint
# model meaningless. No-op for a single-context run.
# @noRd
.ctwasAssertSharedGenes <- function(parsed) {
    byCtx <- split(parsed$trait, parsed$context)
    if (length(byCtx) < 2L) {
        return(invisible())
    }
    geneSets <- map(byCtx, .ctwasSortUnique)
    ref <- geneSets[[1L]]
    mismatch <- names(geneSets)[
        !map_lgl(geneSets, identical, ref)
    ]
    if (length(mismatch) > 0L) {
        msg <- glue(
            "ctwasPipeline: multi-context cTWAS requires the SAME gene set in ",
            "every context (the joint model is only meaningful when genes are ",
            "shared across contexts). Context(s) differing from the ",
            "reference '{names(geneSets)[[1L]]}' ({length(ref)} gene(s)): ",
            "{str_flatten(mismatch, ', ')}."
        )
        abort(msg)
    }
    invisible()
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
#   gene rows  the id is region|study|context|trait|method, and the run's
#              weights carry that gene's trait_pos (falling back to the
#              weight-variant span chrom/p0/p1)
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
    # A weights element that is not a list carries no gene metadata (some
    # callers pass a bare weight vector), so there is nothing to place it by.
    if (!is.list(w)) {
        return(NULL)
    }
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
# that genuinely has no coordinates is not silently fabricated onto chr1.
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
        chrom = withChrPrefix(as.character(parsed$chrom[[1L]])),
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
        finemap = .ctwasRangePayload(fm, geneCoords),
        susieAlpha = .ctwasRangePayload(sa, geneCoords),
        groupPriors = runResult$param,
        regionInfo = .ctwasRangeRegionInfo(runResult$region_info)
    )
}

# region_info as a GRanges. pecotmr builds this table itself with `chr` and a
# variant-derived [start, stop], so the ranges are already there -- this only
# gives them their proper type.
# @noRd
.ctwasRangeRegionInfo <- function(ri) {
    if (is.null(ri) || nrow(ri) == 0L) {
        return(ri)
    }
    cols <- names(ri)
    if (!all(c("chrom", "start", "stop") %in% cols)) {
        return(ri)
    }
    gr <- GenomicRanges::GRanges(
        withChrPrefix(as.character(ri$chrom)),
        IRanges::IRanges(
            start = as.integer(ri$start),
            end = as.integer(ri$stop)
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
# ctwas run is reconstructable from `getFinemap()` / `getSusieAlpha()`.
# @noRd
.ctwasRunToRows <- function(runResult, gwasStudy, method, keepSnps = FALSE) {
    geneIds <- names(runResult$weights)
    if (is.null(geneIds) || length(geneIds) == 0L) {
        return(list())
    }
    parsed <- .ctwasParseGeneIds(geneIds)
    contexts <- unique(parsed$context)
    .ctwasAssertSharedGenes(parsed)
    jointStr <- if (length(contexts) > 1L) {
        str_flatten(sort(unique(contexts)), ",")
    } else {
        NA_character_
    }
    fmDf <- .ctwasAsDf(runResult$finemap_res)
    saDf <- .ctwasAsDf(runResult$susie_alpha_res)
    contextRows <- map(
        contexts,
        .ctwasContextRow,
        parsed = parsed,
        gwasStudy = gwasStudy,
        method = method,
        jointStr = jointStr,
        fmDf = fmDf,
        saDf = saDf,
        runResult = runResult
    )
    snpRow <- if (!keepSnps) {
        NULL
    } else {
        .ctwasSnpRow(gwasStudy, method, jointStr, fmDf, saDf, runResult)
    }
    rows <- c(contextRows, compact(list(snpRow)))
    rows
}

# as.data.frame(), passing NULL through.
# @noRd
.ctwasAsDf <- function(x) {
    if (is.null(x)) NULL else as.data.frame(x)
}

# One (gwasStudy, study, context, method) row-spec for a context. Errors if a
# context mixes multiple QTL studies (one study per context).
# @noRd
.ctwasContextRow <- function(
    cx,
    parsed,
    gwasStudy,
    method,
    jointStr,
    fmDf,
    saDf,
    runResult
) {
    inCx <- parsed$context == cx
    studyCx <- unique(parsed$study[inCx])
    if (length(studyCx) != 1L) {
        msg <- glue(
            "ctwasPipeline: context '{cx}' mixes multiple QTL studies ",
            "({str_flatten(studyCx, ', ')}); one study per context."
        )
        abort(msg)
    }
    idsCx <- parsed$id[inCx]
    list(
        gwasStudy = gwasStudy,
        study = studyCx,
        context = cx,
        method = method,
        jointContexts = jointStr,
        entry = .ctwasMkEntry(
            .ctwasSubsetById(fmDf, idsCx),
            .ctwasSubsetById(saDf, idsCx),
            runResult
        )
    )
}

# The SNP-level row-spec (study = context = "SNP"), or NULL when no SNP-level
# finemap / susie-alpha rows exist.
# @noRd
.ctwasSnpRow <- function(gwasStudy, method, jointStr, fmDf, saDf, runResult) {
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
        jointContexts = jointStr,
        entry = .ctwasMkEntry(snpFm, snpSa, runResult)
    )
}

# Assemble accumulated per-run row-specs (from .ctwasRunToRows) into a single
# CtwasResult; the jointContexts column is omitted when every row is
# single-context. Shared by ctwasPipeline (across methods) and asCtwasResult.
# @noRd
.ctwasRowsToResult <- function(rows) {
    if (length(rows) == 0L) {
        msg <- glue(
            "cTWAS: no genes were modeled (the weight source produced ",
            "no usable gene weights)."
        )
        abort(msg)
    }
    jointContexts <- map_chr(rows, "jointContexts")
    CtwasResult(
        gwasStudy = map_chr(rows, "gwasStudy"),
        study = map_chr(rows, "study"),
        context = map_chr(rows, "context"),
        method = map_chr(rows, "method"),
        entry = map(rows, "entry"),
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
# assembled weight ids. Errors on a mix (the granular path is one method/run).
# @noRd
.ctwasMethodFromWeights <- function(weights) {
    if (is.null(weights) || length(weights) == 0L) {
        msg <- glue(
            "asCtwasResult: the finemap result carries no weights to derive a ",
            "method from."
        )
        abort(msg)
    }
    m <- unique(.ctwasParseGeneIds(names(weights))$method)
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
#'   path. The GWAS study is read from \code{z_snp$study} and the (single)
#'   weight method from the gene ids.
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
    method <- .ctwasMethodFromWeights(finemapResult$weights)
    rows <- .ctwasRunToRows(
        finemapResult,
        gwasStudy = gwasStudy,
        method = method,
        keepSnps = keepSnps
    )
    .ctwasRowsToResult(rows)
}

# Subset a TwasWeights collection to rows whose `method` matches the
# resolved method. Used to enforce the "one ctwas gene per (study,
# context, trait)" semantics -- the legacy pipeline fed a single
# best-CV-method weight per gene; the new S4 TwasWeights may carry
# many methods, but ctwas should only see one.
# @noRd
.ctwasFilterMethod <- function(tw, method) {
    keep <- which(as.character(tw$method) == method)
    if (length(keep) == 0L) {
        return(NULL)
    }
    # Row subset carries every column forward (joint* / region / ...); the old
    # hand-listed rebuild silently dropped them.
    methods::initialize(tw[keep, ], ldSketch = getLdSketch(tw))
}

# Build the per-variant Z data.frame ctwas expects from a GwasSumStats.
# Stacks each study row's GRanges via the shared `.entryToSumstatDf`
# helper (R/sumstatsQc.R), then projects to ctwas's column shape and
# bolts on the `study` column ctwas uses to disambiguate stacked rows.
# @noRd
# Match variant ids against the LD panel's ids the way every other LD consumer
# in the package does -- by (chrom, pos, allele) with ref/alt swaps tolerated --
# rather than by exact string identity.
#
# ctwas already harmonizes TWAS weights this way (see
# `.ctwasHarmonizeWeights()` below). The GWAS side used plain string joins, so
# one variant could be kept for the weights and dropped for the z-scores.
# `removeStrandAmbiguous = FALSE` matches `.ldFromSketchMatch()`, the shared
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
        study = as.character(gwasSumStats$study)[[i]]
    )
}

.ctwasBuildZSnp <- function(gwasSumStats, panelIds) {
    pieces <- map(
        seq_len(nrow(gwasSumStats)),
        .ctwasZSnpPiece,
        gwasSumStats = gwasSumStats
    )
    .ctwasHarmonizeZToPanel(bind_rows(pieces), panelIds)
}

# Derive the single-row region_info from the LD sketch's snpInfo
# (min/max BP per chromosome). The sketch is assumed to cover exactly
# one block.
# @noRd
# @noRd
.ctwasEntryPositions <- function(i, gss) {
    gr <- gss[[i]]
    list(
        pos = as.integer(GenomicRanges::start(gr)),
        chrs = as.character(GenomicRanges::seqnames(gr))
    )
}

# @noRd
.ctwasConcatInt <- function(pieces) {
    if (length(pieces) == 0L) {
        return(integer(0))
    }
    list_c(pieces)
}

# @noRd
.ctwasConcatChr <- function(pieces) {
    if (length(pieces) == 0L) {
        return(character(0))
    }
    list_c(pieces)
}

.ctwasBuildSingleRegionInfo <- function(regionId, gss) {
    # Derive the block's [start, stop] from the GWAS variants actually in this
    # block (the GwasSumStats entry GRanges) -- NOT the LD sketch. When many
    # blocks share one whole-chromosome LD payload (the common one-file-per-chr
    # layout), getSnpInfo(ldSketch) spans the entire chromosome, so every region
    # would collapse to the same whole-chromosome [start, stop] and every SNP
    # would be assigned to every region (inflating SNP group_size N-fold and
    # diluting the gene prior to ~0).
    entries <- map(seq_len(nrow(gss)), .ctwasEntryPositions, gss = gss)
    pos <- .ctwasConcatInt(map(entries, "pos"))
    chrs <- .ctwasConcatChr(map(entries, "chrs"))
    # Emptiness is checked FIRST: `pos` and `chrs` are filled from the same
    # GRanges in the same loop, so an empty block has zero chromosomes too,
    # and the chromosome check below would report it as "spans multiple
    # chromosomes ()" -- which is both wrong and unactionable.
    if (length(pos) == 0L) {
        msg <- glue(
            "ctwasPipeline: GwasSumStats block '{regionId}' has no variants ",
            "to define region bounds."
        )
        abort(msg)
    }
    chr <- unique(as.integer(
        str_remove(chrs, regex("^chr", ignore_case = TRUE))
    ))
    if (length(chr) != 1L) {
        msg <- glue(
            "ctwasPipeline: GwasSumStats block '{regionId}' spans multiple ",
            "chromosomes ({str_flatten(chr, ', ')})."
        )
        abort(msg)
    }
    tibble(
        region_id = regionId,
        chrom = chr,
        start = min(pos),
        stop = max(pos)
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
    gr <- .ldSketchRanges(gwasLd)
    mc <- S4Vectors::mcols(gr)
    # `.ldSketchMatchIds()`, not the raw SNP label: this id becomes the R
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
        chrom = as.integer(.ldSketchChrom(gwasLd)),
        id = .ldSketchMatchIds(gwasLd),
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
#   variance : named numeric vector of per-variant dosage variance from the
#              LD reference, used to scale non-standardized TWAS weights to
#              the correlation scale ctwas expects.
#   sketch   : the block's narrowed LD sketch, retained because the n x n
#              correlation is deliberately NOT built here.
#
# ctwas pulls a region's LD only while fine-mapping, and only for the regions
# that survive screening, so building every block's matrix up front spent
# O(n^2) apiece on panels that are frequently never read at all.
# @noRd
.ctwasPanelFor <- function(gwasLd) {
    # Share the validator with `.ldFromSketch()`, the entry point every other
    # pipeline uses: skipping the guard meant a NULL or non-panel sketch
    # surfaced as "unable to find an inherited method for 'getSnpInfo'"
    # instead of saying the LD reference was missing.
    .ldFromSketchValidate(gwasLd, "ctwasPipeline")
    snpInfo <- .ctwasSnpInfoForBlock(gwasLd)
    geno <- .ldSketchDosage(
        gwasLd,
        seq_len(nrow(snpInfo)),
        meanImpute = TRUE
    )
    list(
        snpInfo = snpInfo,
        variance = set_names(
            apply(geno, 2, stats::var, na.rm = TRUE),
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

# The correlation over `vids` alone. Identical to slicing the full matrix --
# mean imputation is per column, so correlating a subset of the dosages gives
# the same pairwise values -- at O(|vids|^2) instead of O(n^2).
# @noRd
.ctwasPanelLdSubset <- function(panel, vids) {
    idx <- match(vids, panel$snpInfo$id)
    if (anyNA(idx)) {
        msg <- glue(
            "ctwasPipeline: {sum(is.na(idx))} weight variant(s) absent from ",
            "the block's LD panel."
        )
        abort(msg)
    }
    if (!is.null(panel$R)) {
        return(panel$R[vids, vids, drop = FALSE])
    }
    .ctwasPanelLdFor(panel, idx, vids)
}

# @noRd
.ctwasPanelLdFor <- function(panel, idx, ids) {
    geno <- .ldSketchDosage(panel$sketch, idx, meanImpute = TRUE)
    `dimnames<-`(
        computeLd(geno, method = "sample"),
        list(ids, ids)
    )
}

# Whether to build the block's full correlation once and slice it per gene,
# or to compute one small correlation per gene. Per-gene submatrices are what
# keep a block ctwas never fine-maps from ever forming its n x n matrix, but
# g genes on |vids| variants each cost sum(|vids|^2) against n^2 for the
# whole block, so the full matrix wins once the weights are dense enough.
# Returns the full matrix in that case and NULL in the other.
# @noRd
.ctwasRWgtSource <- function(ldPanel, nVids) {
    # Already built: slicing it is free, so there is nothing to weigh up.
    if (!is.null(ldPanel$R)) {
        return(ldPanel$R)
    }
    n <- as.numeric(nrow(ldPanel$snpInfo))
    if (length(nVids) == 0L || sum(as.numeric(nVids)^2) < n^2) {
        return(NULL)
    }
    .ctwasPanelLd(ldPanel)
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
        chrom = as.integer(parsed$chrom),
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
    needed <- c("alpha", "mu", "X_column_scale_factors")
    all(is_in(needed, names(fits)))
}

# Renormalize SuSiE TWAS weights over the kept variant set. When some
# variants got dropped by allele harmonization / panel intersection,
# the posterior `alpha` values from the original fit no longer sum to
# 1 over the kept variants. We re-softmax `lbf_variable[, keptIdx]`
# into a renormalized alpha, sign-flip the rows of `mu[, keptIdx]` to
# match the panel's allele frame (carrying over the per-variant sign
# flip already applied to `harmonizedW`), and recompute the per-variant
# weight as `colSums(alpha * mu_subset) / X_column_scale_factors_subset`.
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
    # Fields are READ with `[[`: `$` on a list falls back to prefix matching,
    # so `fits$mu` would silently return `mu2` on a fit that has no `mu`.
    rawAlpha <- fits[["alpha"]]
    mu <- fits[["mu"]]
    xCol <- fits[["X_column_scale_factors"]]
    if (is.null(rawAlpha) || is.null(mu) || is.null(xCol)) {
        return(NULL)
    }
    # susieInf / susieAsh carry an infinitesimal term: coef.susie is
    # colSums(alpha * mu) / scale + theta / scale. Recomputing the weight from
    # alpha and mu alone would silently drop theta, changing what the weight
    # means. Leave those fits to the plain subset-w-and-R path, which is
    # unbiased and only loses power.
    if (!is.null(fits[["theta"]]) || !is.null(fits[["omega_weights"]])) {
        return(NULL)
    }
    alpha <- as.matrix(rawAlpha)
    if (
        ncol(alpha) != length(origVids) ||
            ncol(mu) != length(origVids) ||
            length(xCol) != length(origVids)
    ) {
        # Fit-vs-entry dimension mismatch -- including a null_weight fit, whose
        # alpha carries an extra column. Skip rather than mis-slice.
        return(NULL)
    }
    # Per-variant sign flip applied by allele harmonization. NaN signs
    # (origW == 0) default to +1.
    rawSign <- sign(harmonizedW / origW[keptIdx])
    signFlip <- replace(rawSign, !is.finite(rawSign), 1)
    newAlpha <- .ctwasRenormAlpha(alpha, keptIdx)
    if (is.null(newAlpha)) {
        return(NULL)
    }
    muSub <- sweep(mu[, keptIdx, drop = FALSE], 2L, signFlip, `*`)
    rawScale <- xCol[keptIdx]
    # Guard against zero scale factors (shouldn't happen in practice).
    xColSub <- replace(rawScale, rawScale == 0, 1)
    as.numeric(colSums(newAlpha * muSub) / xColSub)
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
# R_wgt (per-gene LD submatrix), and gene metadata. ctwas's compute_gene_z
# pulls rownames(wgt) for the SNP IDs and computes z.gene = crossprod(wgt,
# z.s) / sqrt(t(wgt) %*% R_wgt %*% wgt), so wgt must be a numeric matrix
# (not a vector) and R_wgt must be the LD submatrix over the same SNPs.
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

# The ctwas per-gene weight entries for one block.
#
# Two variant sets, because they answer different questions. `gwasSnpIds` is
# the GLOBAL set: it bounds each gene's cis SPAN, which has to cover every
# block the gene reaches for boundary detection to work. `regionSnpIds` is
# this block's own set: it bounds the weight vector actually FITTED, because
# ctwas fine-maps one region at a time.
#
# So anything feeding the SPAN is read off the whole panel rather than off
# `ldPanel`, whose sketch is narrowed to this block -- `globalPanelInfo` for
# the coordinates, `globalVariance` for the scaling, and the same table for
# allele harmonization, since an out-of-block weight variant has to survive
# matching to be counted in p0/p1 at all. Reading any of them off the block
# would clip a boundary gene back to its home region, which is exactly what
# p0/p1 exist to prevent. The allele columns agree with the block's own
# table wherever the two overlap, so nothing about the in-block result
# changes.
#
# Genes are prepared in one sweep and given their `R_wgt` in a second, so
# that choice is made with every gene's fitted size in hand -- which is what
# lets a block no gene needs densely skip its n x n correlation entirely.
# .ctwasRWgtSource() picks, and says why. Variants absent from the panel are
# dropped from that gene's row set.
# @noRd
.ctwasBuildWeights <- function(
    twasWeights,
    ldPanel,
    fineMappingResult = NULL,
    twasWeightCutoff = 0,
    csMinCor = 0.8,
    minPipCutoff = 0,
    maxNumVariants = Inf,
    gwasSnpIds = NULL,
    regionSnpIds = NULL,
    globalPanelInfo = NULL,
    globalVariance = NULL
) {
    spanInfo <- globalPanelInfo %||% ldPanel$snpInfo
    panelSnps <- .ctwasPanelSnpsForGwas(spanInfo$id, gwasSnpIds)
    ctx <- list(
        ldPanel = ldPanel,
        spanInfo = spanInfo,
        panelSnps = panelSnps,
        regionSnps = regionSnpIds,
        refVariants = .ctwasRefVariants(spanInfo),
        spanVariance = globalVariance %||% ldPanel$variance,
        fineMappingResult = fineMappingResult,
        cutoffs = list(
            twasWeightCutoff = twasWeightCutoff,
            csMinCor = csMinCor,
            minPipCutoff = minPipCutoff,
            maxNumVariants = maxNumVariants
        )
    )
    prepared <- compact(map(
        seq_len(nrow(twasWeights)),
        .ctwasGenePrepared,
        twasWeights = twasWeights,
        ctx = ctx
    ))
    fullR <- .ctwasRWgtSource(ldPanel, map_int(prepared, .ctwasPreparedSize))
    genes <- map(
        prepared,
        .ctwasGeneEntryFor,
        ldPanel = ldPanel,
        fullR = fullR,
        spanInfo = spanInfo
    )
    # Later genes overwrite an earlier one sharing a key, which is what the
    # `out[[g$key]] <- ...` loop did.
    keyed <- set_names(map(genes, "entry"), map_chr(genes, "key"))
    keyed[!duplicated(names(keyed), fromLast = TRUE)]
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

# study/context/trait/method identity for gene row `i` + its collection key.
# @noRd
.ctwasGeneMeta <- function(twasWeights, i) {
    study <- as.character(twasWeights$study)[[i]]
    context <- as.character(twasWeights$context)[[i]]
    trait <- as.character(twasWeights$trait)[[i]]
    method <- as.character(twasWeights$method)[[i]]
    list(
        study = study,
        context = context,
        trait = trait,
        method = method,
        traitPos = .ctwasTraitPosAt(twasWeights, i),
        key = as.character(glue("{study}|{context}|{trait}|{method}"))
    )
}

# Steps 1-2: allele-harmonize a gene's (variantIds, weights) against the LD
# panel (matching by chrom/pos with sign/strand-flip detection; weights
# sign-flipped for swapped frames), then restrict to the panel variant set.
# Returns list(vids, w, keptIdx, origVids, origW), or NULL when nothing
# survives.
# @noRd
.ctwasAlignGeneWeights <- function(parts, refVariants, panelSnps) {
    wr <- .rowResolveWeights(parts)
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

# Steps 3-4: SuSiE alpha renormalization (when the kept variant set shrank the
# fit) and variance scaling for non-standardized weights (w * sqrt(per-variant
# genotype variance from the LD panel)). Returns the adjusted weight vector.
# @noRd
.ctwasAdjustGeneWeights <- function(parts, aligned, spanVariance) {
    fits <- .rowFits(parts)
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
            ".ctwasBuildWeights: missing genotype variance for ",
            "{sum(is.na(varLookup))} variant(s) in the LD panel."
        )
        abort(msg)
    }
    renormalized * sqrt(varLookup)
}

# The ctwas per-gene weight entry (weight matrix, LD submatrix, chrom/BP span,
# plus identity metadata).
#
# `spanVids` is the variant set the chrom/p0/p1 span is measured over, and it is
# deliberately allowed to be WIDER than the fitted `vids`: ctwas fine-maps one
# region at a time, so `wgt` / `R_wgt` must stay inside that region, while
# `ctwas::get_boundary_genes` reads p0/p1 to decide which genes straddle a
# region boundary. Measuring the span over the fitted subset would clip a
# boundary gene down to its home region and hide it from merge_regions.
# @noRd
.ctwasGeneEntry <- function(
    vids,
    w,
    rWgt,
    meta,
    spanVids = vids,
    spanInfo
) {
    # Span coordinates come from the whole-panel table, while `rWgt` was
    # computed on the block's own panel: ctwas fine-maps one region at a
    # time, but a boundary gene's span has to reach past it.
    panelInfo <- spanInfo
    rowIdx <- match(spanVids, panelInfo$id)
    list(
        wgt = matrix(w, ncol = 1L, dimnames = list(vids, "wgt")),
        R_wgt = rWgt,
        type = meta$context,
        context = meta$context,
        gene_name = meta$trait,
        study = meta$study,
        method = meta$method,
        n_wgt = length(vids),
        chrom = as.integer(panelInfo$chrom[[rowIdx[1L]]]),
        p0 = min(as.integer(panelInfo$pos[rowIdx])),
        p1 = max(as.integer(panelInfo$pos[rowIdx])),
        molecular_id = meta$trait,
        weight_name = str_c(meta$context, meta$context, sep = "_"),
        # The GENE's own coordinate, carried alongside the weight-variant span
        # (chrom/p0/p1) so the result payloads can be ranged on the gene rather
        # than on wherever its weight variants happen to fall.
        trait_pos = meta$traitPos
    )
}

# The gene's own position for row i, or NULL when the weight source carries no
# traitPos provenance.
# @noRd
.ctwasTraitPosAt <- function(twasWeights, i) {
    if (!is_in("traitPos", .tupleColumnNames(twasWeights))) {
        return(NULL)
    }
    tp <- getTraitPosition(twasWeights)
    if (!methods::is(tp, "GRanges") || length(tp) < i) {
        return(NULL)
    }
    tp[i]
}

# Build one gene's ctwas weight record: align -> adjust -> smart-filter
# (PIP/CS + magnitude + cap) -> entry. Returns list(key, entry), or NULL when
# the gene contributes no usable variants.
# @noRd
.ctwasGenePrepared <- function(i, twasWeights, ctx) {
    # Polymorphic: the weight source may be a TwasWeights or a
    # FineMappingResult, so the payload comes from the class-aware bridge.
    parts <- .rowParts(twasWeights, i)
    aligned <- .ctwasAlignGeneWeights(parts, ctx$refVariants, ctx$panelSnps)
    if (is.null(aligned)) {
        return(NULL)
    }
    w <- .ctwasAdjustGeneWeights(parts, aligned, ctx$spanVariance)
    meta <- .ctwasGeneMeta(twasWeights, i)
    finemapAux <- .ctwasGetFinemapAux(
        ctx$fineMappingResult,
        meta$study,
        meta$context,
        meta$trait,
        meta$method
    )
    kept <- .ctwasFilterVariants(
        vids = aligned$vids,
        w = w,
        finemapAux = finemapAux,
        twasWeightCutoff = ctx$cutoffs$twasWeightCutoff,
        csMinCor = ctx$cutoffs$csMinCor,
        minPipCutoff = ctx$cutoffs$minPipCutoff,
        maxNumVariants = ctx$cutoffs$maxNumVariants
    )
    if (length(kept) < 1L) {
        return(NULL)
    }
    # The weight vector ctwas fits must live inside the region being
    # fine-mapped: susie_rss gets that region's LD, and a gene reaching
    # outside it yields a non-finite ELBO (and a non-finite fit_EM
    # log-likelihood upstream). The SPAN stays the gene's full cis extent --
    # see .ctwasGeneEntry -- so a boundary gene is still detected and can be
    # recovered by merge_regions.
    inRegion <- .ctwasInRegion(kept$vids, ctx$regionSnps)
    if (!any(inRegion)) {
        return(NULL)
    }
    # No LD yet. R_wgt is built in a second sweep, once every gene's fitted
    # size is known and the block can choose between forming one full
    # correlation and one small correlation per gene.
    list(
        key = meta$key,
        vids = kept$vids[inRegion],
        w = kept$w[inRegion],
        spanVids = kept$vids,
        meta = meta
    )
}

# @noRd
.ctwasPreparedSize <- function(g) {
    length(g$vids)
}

# One prepared gene's ctwas entry. `fullR` is the block's correlation when
# .ctwasRWgtSource() chose to form it, and NULL when each gene computes its
# own small one.
# @noRd
.ctwasGeneEntryFor <- function(g, ldPanel, fullR, spanInfo) {
    rWgt <- if (is.null(fullR)) {
        .ctwasPanelLdSubset(ldPanel, g$vids)
    } else {
        fullR[g$vids, g$vids, drop = FALSE]
    }
    list(
        key = g$key,
        entry = .ctwasGeneEntry(
            g$vids,
            g$w,
            rWgt,
            g$meta,
            spanVids = g$spanVids,
            spanInfo = spanInfo
        )
    )
}

# Which of a gene's weight variants lie in the region being fine-mapped. A NULL
# region set means the caller supplied none, so nothing is restricted.
# @noRd
.ctwasInRegion <- function(vids, regionSnps) {
    if (is.null(regionSnps)) {
        return(rep(TRUE, length(vids)))
    }
    is_in(vids, as.character(regionSnps))
}

# Look up the per-(study, context, trait, method) PIP vector and the
# 95% credible-set membership / purity for one gene from the supplied
# FineMappingResult. Returns NULL when no FineMappingResult was passed
# or no matching tuple exists. Output is a list with:
#   pip       : named numeric vector keyed by variant_id
#   csMembers : list of character vectors (one per CS at 95% coverage)
#   csPurity  : numeric vector aligned with csMembers
# @noRd
.ctwasGetFinemapAux <- function(
    fineMappingResult,
    study,
    context,
    trait,
    method
) {
    if (is.null(fineMappingResult)) {
        return(NULL)
    }
    cols <- .tupleColumnNames(fineMappingResult)
    selectors <- c(
        list(study = study, method = method),
        compact(list(
            context = if (is_in("context", cols)) context,
            trait = if (is_in("trait", cols)) trait
        ))
    )
    selArgs <- c(list(fineMappingResult), selectors)
    entry <- try_fetch(
        exec(getFineMappingResult, !!!selArgs),
        error = function(cnd) NULL
    )
    if (is.null(entry)) {
        return(NULL)
    }
    tl <- getTopLoci(entry, raw = TRUE)
    if (nrow(tl) == 0L) {
        return(NULL)
    }
    # `pip` is part of the topLoci schema -- fineMappingRow() requires it on
    # any non-empty table, and the zero-row case returned above -- so there is
    # no pip-less frame to guard against here.
    pip <- set_names(as.numeric(tl$pip), as.character(tl$variant_id))
    cs <- .ctwasCsMembership(tl)
    list(pip = pip, csMembers = cs$csMembers, csPurity = cs$csPurity)
}

# Per-CS membership + purity at 95% coverage from a topLoci table. cs_95 stores
# `<method>_<idx>` where idx == 0 means "not in any CS"; cs_95_purity (when
# present) broadcasts one purity value across a CS's rows.
# @noRd
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

.ctwasCsMembership <- function(tl) {
    if (!is_in("cs_95", names(tl))) {
        return(list(csMembers = list(), csPurity = numeric(0)))
    }
    csIdx <- suppressWarnings(as.integer(str_remove(tl$cs_95, "^.*_")))
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

# Variants that must survive the cap: members of any high-purity (>= csMinCor)
# credible set, plus any with PIP > minPipCutoff. Empty when no finemapAux.
# @noRd
# Variants rescued by credible set `k`, when that set is pure enough.
# @noRd
.ctwasCsRescued <- function(k, finemapAux, vids, csMinCor) {
    if (is.na(finemapAux$csPurity[k]) || finemapAux$csPurity[k] < csMinCor) {
        return(character(0))
    }
    intersect(finemapAux$csMembers[[k]], vids)
}

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

# Build z_gene data.frame from a TWAS-Z GRanges (output of
# causalInferencePipeline). One row per (qtlStudy, context, trait,
# method, gwasStudy) tuple.
# @noRd
.ctwasBuildZGene <- function(twasZ) {
    mc <- as.data.frame(S4Vectors::mcols(twasZ))
    # Base data.frame (not tibble): z_gene is indexed positionally by the ctwas
    # engine (df[, "col"] -> vector); a tibble breaks it. See
    # .ctwasSnpInfoForBlock.
    data.frame(
        id = as.character(glue(
            "{mc$qtlStudy}|{mc$context}|{mc$trait}|{mc$method}"
        )),
        z = as.numeric(mc$twasZ),
        type = as.character(mc$context),
        context = as.character(mc$context),
        gene_name = as.character(mc$trait),
        study = as.character(mc$qtlStudy),
        method = as.character(mc$method),
        stringsAsFactors = FALSE
    )
}

# A cTWAS payload's `LD_file` token does double duty: ctwas asserts the file
# exists, and pecotmr dispatches on the same string into the cached per-region
# LD panels. Both make the token machine-specific, so a payload that was
# serialised (the bundled ctwas*Example objects) carries a path from whatever
# machine built it. Re-point it at this installation on the way in, using the
# same "pecotmr://extdata/<path>" convention the genotype handles use.
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
    resolved <- map_chr(stored, .resolveCtwasLdToken)
    resolved <- .ctwasRemintMissingTokens(
        resolved,
        .ctwasLdMapRegionIds(ldMap, length(resolved))
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

# Only "pecotmr://" tokens move; an ordinary path is the caller's own and is
# left alone (an unresolvable bundled reference is an error worth hearing).
# @noRd
.resolveCtwasLdToken <- function(token) {
    if (!str_detect(token, "^pecotmr://")) {
        return(token)
    }
    .resolveGenotypeResourcePath(token)
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

# The panel cache the loader closures were built over, or NULL when the
# payload carries loaders this package did not create.
# @noRd
.ctwasCachedPanels <- function(payload) {
    loader <- payload$LD_loader_fun
    if (!is.function(loader)) {
        return(NULL)
    }
    attr(loader, "ldPanelsByRegion")
}

# Multi-block LD loader for ctwas. ctwas invokes
# `LD_loader_fun(LD_file)` per region during region_data assembly and
# fine-mapping; we dispatch by `LD_file` (the same string set on
# `LD_map$LD_file`) into the cached per-sketch ldPanel.
# @noRd
.ctwasMultiBlockLdLoader <- function(ldPanelsByRegion) {
    force(ldPanelsByRegion)
    fn <- function(LD_file) {
        panel <- ldPanelsByRegion[[LD_file]]
        if (is.null(panel)) {
            msg <- glue(
                "ctwasPipeline LD loader: no cached panel for ",
                "LD_file = '{LD_file}'"
            )
            abort(msg)
        }
        # Built here rather than up front: this is the only place a region's
        # full LD is actually needed, and ctwas reaches it only for regions
        # that survived screening.
        .ctwasPanelLd(panel)
    }
    # Published explicitly so .ctwasCachedPanels can recover the cache from a
    # loader we built, instead of looking the name up inside the closure's
    # environment.
    `attr<-`(fn, "ldPanelsByRegion", ldPanelsByRegion)
}

# Multi-block SNP-info loader for ctwas. Mirrors the LD loader.
# @noRd
.ctwasMultiBlockSnpInfoLoader <- function(ldPanelsByRegion) {
    force(ldPanelsByRegion)
    function(LD_file) {
        panel <- ldPanelsByRegion[[LD_file]]
        if (is.null(panel)) {
            msg <- glue(
                "ctwasPipeline snpInfo loader: no cached panel for ",
                "LD_file = '{LD_file}'"
            )
            abort(msg)
        }
        panel$snpInfo
    }
}

# Derive the LD_file token for ctwas from a GenotypeHandle. We point
# at the on-disk file that already backs the sketch's data, so the
# `file.exists(LD_map$LD_file)` assertion in ctwas::ctwas_sumstats
# passes WITHOUT pecotmr doing any new I/O. The token also serves as
# the dispatch key for the multi-block LD / snpInfo loaders, so two
# blocks sharing the same on-disk LD payload share one cached panel.
# @noRd
.ctwasLdPanelKey <- function(sketch) {
    handle <- .ldSketchHandle(sketch)
    fmt <- getFormat(handle)
    stem <- .genotypeReadPath(handle)
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
            "ctwasPipeline: could not derive an existing LD-file token for ",
            "the GenotypeHandle (format={fmt}, path={stem}). Looked for: ",
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
# the cached panels onto them, exactly as it does for a "pecotmr://" move.
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

# Build a per-block snpInfo table restricted to variants present in the
# GwasSumStats entry. Mirrors `.ctwasSnpInfoForBlock` but restricts to
# the block's GWAS variants (intersected against the cached panel) so
# snp_map[[region_id]] is sized to the block, not the whole panel.
# @noRd
# One entry's variant ids, or none when it does not carry a SNP column.
# @noRd
.ctwasEntrySnpIds <- function(i, gwasSumStats) {
    mc <- S4Vectors::mcols(gwasSumStats[[i]])
    if (!is_in("SNP", colnames(mc))) {
        return(character(0))
    }
    as.character(mc$SNP)
}

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

# The weight method of one assembled TwasWeights record (character).
# @noRd
.ctwasMethodChr <- function(tw) {
    as.character(tw$method)
}

# The block id whose window contains anchor variant `i`, or NA when unplaced.
# @noRd
.ctwasBlockIdForVariant <- function(i, aPos, aChr, bChr, bS, bE, ids) {
    if (is.na(aPos[[i]])) {
        return(NA_character_)
    }
    hit <- which(bChr == aChr[[i]] & aPos[[i]] >= bS & aPos[[i]] < bE)
    if (length(hit) > 0L) ids[[hit[[1L]]]] else NA_character_
}

# The trait field of a split gene id: everything between context and the final
# method field, rejoined on "|" so a "|"-bearing trait survives.
# @noRd
.ctwasTraitField <- function(p) {
    str_flatten(p[4:(length(p) - 1L)], "|")
}

# The method field (last component) of a split gene id.
# @noRd
.ctwasMethodField <- function(p) {
    p[[length(p)]]
}

# The sorted unique gene (trait) set of one context.
# @noRd
.ctwasSortUnique <- function(g) {
    sort(unique(g))
}
