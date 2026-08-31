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
#' @section LD-sketch identity check: Per block: \code{getLdSketch(twasWeights)}
#'   (when non-NULL) must match \code{getLdSketch(gwasSumStats)}. Mismatch is a
#'   hard error.
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
#' @param thin,niterPrefit,niter,L Pass-throughs to
#'   \code{ctwas::ctwas_sumstats}.
#' @param groupPriorVarStructure Pass-through (defaults \code{"shared_type"}).
#' @param ncore Number of cores. Default \code{1}.
#' @param twasWeightCutoff Numeric (length 1). Drop variants with \code{|weight|
#'   < twasWeightCutoff} from each gene's weight matrix before ctwas sees it.
#'   Default \code{0} (no filter).
#' @param csMinCor Numeric (length 1). When \code{fineMappingResult} is
#'   provided, variants belonging to any 95\% credible set with purity
#'   (\code{min_abs_corr}) \code{>= csMinCor} are marked as must-keep and
#'   survive the per-gene cap. Default \code{0.8}. Ignored without a
#'   \code{fineMappingResult}.
#' @param minPipCutoff Numeric (length 1). When \code{fineMappingResult} is
#'   provided, variants with PIP greater than \code{minPipCutoff} are marked as
#'   must-keep and survive the per-gene cap. Default \code{0} (no PIP rescue).
#'   Ignored without a \code{fineMappingResult}.
#' @param maxNumVariants Numeric (length 1). Cap on per-gene variant count. When
#'   the gene has more variants than this, keep all must-keep variants and fill
#'   remaining slots by descending PIP (when available) or descending
#'   \code{|weight|}. Default \code{Inf} (no cap).
#' @param fallbackToPrefit Logical (length 1). Forwarded to
#'   \code{\link{estCtwasParam}}. When \code{TRUE}, ctwas's accurate-EM NaN
#'   failure is recovered by falling back to the prefit estimates (mirrors the
#'   legacy ctwas_2 workaround on underpowered data). Default \code{FALSE}.
#' @param keepSnps Logical (length 1). When \code{TRUE}, retain the
#'   context-agnostic SNP background of each run as one extra \code{CtwasResult}
#'   row (\code{study = context = "SNP"}, mirroring cTWAS's own \code{"SNP"}
#'   group) so the full ctwas output is reconstructable from
#'   \code{\link{getFinemap}} / \code{getSusieAlpha}. Default \code{FALSE} --
#'   the SNP rows are the null background and are dropped from the structured
#'   gene-level result.
#' @param mergeBoundary Logical (length 1). When \code{TRUE}, run
#'   \code{\link{mergeCtwasBoundaryRegions}} after fine-mapping each run: a
#'   high-PIP gene whose cis window straddles an LD-block boundary has its
#'   adjacent regions merged and re-fine-mapped (the legacy default-off
#'   \code{ctwas_3} post-processing). Default \code{FALSE}.
#' @param mergePipThresh Numeric (length 1). PIP threshold for selecting which
#'   boundary genes to merge (\code{\link{mergeCtwasBoundaryRegions}}
#'   \code{pipThresh}). Default \code{0.5}. Ignored unless \code{mergeBoundary =
#'   TRUE}.
#' @param mergeFilterCs Logical (length 1). Require the boundary gene to be in a
#'   credible set to be selected. Default \code{FALSE}. Ignored unless
#'   \code{mergeBoundary = TRUE}.
#' @param mergeMaxSNP Numeric (length 1). Per-merged-region SNP cap. Default
#'   \code{Inf}. Ignored unless \code{mergeBoundary = TRUE}.
#' @param ... Additional arguments forwarded to \code{ctwas::ctwas_sumstats}.
#' @return A \code{\link{CtwasResult}} collection: one row per \code{(gwasStudy,
#'   study, context, method)}. A single-context run is one row per method; a
#'   multi-context (joint) run emits per-context rows sharing the same
#'   \code{jointContexts} set and the jointly-estimated group priors. Each row's
#'   \code{\link{CtwasResultEntry}} payload carries that context's per-gene
#'   fine-mapping posteriors (\code{finemap}), the run's \code{param}, and its
#'   \code{regionInfo}. For the raw \code{ctwas::finemap_regions} list (e.g. to
#'   feed \code{\link{mergeCtwasBoundaryRegions}}), call the granular
#'   \code{\link{assembleCtwasInputs}} \eqn{\to} \code{\link{estCtwasParam}}
#'   \eqn{\to} \code{\link{screenCtwasRegions}} \eqn{\to}
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
#' gwasByRegion <- summaryStatsQc(gss, mafCutoff = 0.0025)
#' ctwasPipeline(gwasSumStats = gwasByRegion,
#'   twasWeights = list(ctwasWeightsExample), thin = 1, niterPrefit = 3,
#'   niter = 10, min_group_size = 1, min_p_single_effect = 0,
#'   fallbackToPrefit = TRUE)
#' @export
ctwasPipeline <- function(
    gwasSumStats,
    twasWeights,
    twasZ = NULL,
    fineMappingResult = NULL,
    method = NULL,
    thin = 0.1,
    niterPrefit = 3L,
    niter = 30L,
    L = 5L,
    groupPriorVarStructure = c(
        "shared_type",
        "shared_context",
        "shared_nonSNP",
        "shared_all",
        "independent"
    ),
    ncore = 1L,
    twasWeightCutoff = 0,
    csMinCor = 0.8,
    minPipCutoff = 0,
    maxNumVariants = Inf,
    fallbackToPrefit = FALSE,
    keepSnps = FALSE,
    mergeBoundary = FALSE,
    mergePipThresh = 0.5,
    mergeFilterCs = FALSE,
    mergeMaxSNP = Inf,
    ...
) {
    groupPriorVarStructure <- arg_match(groupPriorVarStructure)
    .ctwasRequireNamedLists(gwasSumStats, twasWeights)
    methods <- .ctwasResolveMethods(twasWeights, method)
    gwasStudy <- .ctwasGwasStudy(gwasSumStats)
    cfg <- as.list(environment())
    cfg$dots <- list(...)
    rows <- list_flatten(map(methods, .ctwasRunMethod, cfg = cfg))
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

# One cTWAS run for method `m`: assemble inputs -> estimate params -> screen ->
# fine-map (optionally boundary-merge) -> per-context row-specs. `cfg` bundles
# the ctwasPipeline arguments (incl. `dots` = the forwarded `...`).
# @noRd
.ctwasRunMethod <- function(m, cfg) {
    inputs <- assembleCtwasInputs(
        gwasSumStats = cfg$gwasSumStats,
        twasWeights = cfg$twasWeights,
        twasZ = cfg$twasZ,
        fineMappingResult = cfg$fineMappingResult,
        method = m,
        twasWeightCutoff = cfg$twasWeightCutoff,
        csMinCor = cfg$csMinCor,
        minPipCutoff = cfg$minPipCutoff,
        maxNumVariants = cfg$maxNumVariants
    )
    estArgs <- c(
        list(
            inputs,
            thin = cfg$thin,
            niterPrefit = cfg$niterPrefit,
            niter = cfg$niter,
            groupPriorVarStructure = cfg$groupPriorVarStructure,
            ncore = cfg$ncore,
            fallbackToPrefit = cfg$fallbackToPrefit
        ),
        cfg$dots
    )
    est <- exec(estCtwasParam, !!!estArgs)
    screenArgs <- c(list(est, L = cfg$L, ncore = cfg$ncore), cfg$dots)
    screened <- exec(screenCtwasRegions, !!!screenArgs)
    finemapArgs <- c(list(screened, L = cfg$L, ncore = cfg$ncore), cfg$dots)
    finemap <- exec(finemapCtwasRegions, !!!finemapArgs)
    if (cfg$mergeBoundary) {
        finemap <- .ctwasMaybeMerge(finemap, cfg)
    }
    .ctwasRunToRows(
        finemap,
        gwasStudy = cfg$gwasStudy,
        method = m,
        keepSnps = cfg$keepSnps
    )
}

# Boundary-gene region merging: split a high-PIP straddling gene's adjacent
# regions and re-fine-map. Merge-transparent downstream (keyed by gene id).
# @noRd
.ctwasMaybeMerge <- function(finemap, cfg) {
    mergeArgs <- c(
        list(
            finemap,
            pipThresh = cfg$mergePipThresh,
            filterCs = cfg$mergeFilterCs,
            maxSNP = cfg$mergeMaxSNP,
            L = cfg$L,
            ncore = cfg$ncore
        ),
        cfg$dots
    )
    exec(mergeCtwasBoundaryRegions, !!!mergeArgs)
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
#'   \code{\link{estCtwasParam}} -> \code{\link{screenCtwasRegions}} ->
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
#' gwasByRegion <- summaryStatsQc(gss, mafCutoff = 0.0025)
#' assembleCtwasInputs(gwasSumStats = gwasByRegion,
#'   twasWeights = list(ctwasWeightsExample))
#' @export
assembleCtwasInputs <- function(
    gwasSumStats,
    twasWeights,
    twasZ = NULL,
    fineMappingResult = NULL,
    method = NULL,
    twasWeightCutoff = 0,
    csMinCor = 0.8,
    minPipCutoff = 0,
    maxNumVariants = Inf
) {
    .ctwasValidateGwasList(gwasSumStats)
    # One single-block GwasSumStats per element, keyed by blockId: the region
    # grid the rest of the assembly walks.
    gwasSumStats <- .ctwasGwasByBlock(gwasSumStats)
    twasWeights <- .ctwasResolveAndValidateWeights(twasWeights, gwasSumStats)
    .ctwasValidateOptional(twasZ, fineMappingResult)
    regionIds <- names(gwasSumStats)
    resolvedMethod <- .ctwasResolveMethod(twasWeights, method)
    fp <- .ctwasFirstPass(regionIds, gwasSumStats, twasWeights)
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
        # nocov start
        msg <- glue(
            "Package 'ctwas' is required for the cTWAS pipeline. ",
            "Install from https://github.com/xinhe-lab/ctwas ."
        )
        abort(msg)
        # nocov end
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

# @noRd
.ctwasPickBlock <- function(i, x) {
    x[i]
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
.ctwasFirstPass <- function(regionIds, gwasSumStats, twasWeights) {
    ldPanelsByRegion <- list()
    zSnpPieces <- list()
    regionInfoPieces <- list()
    snpMap <- list()
    ldFileByRegion <- set_names(character(length(regionIds)), regionIds)
    for (rid in regionIds) {
        gss <- gwasSumStats[[rid]]
        tw <- twasWeights[[rid]]
        gwasLd <- getLdSketch(gss)
        if (is.null(gwasLd)) {
            msg <- glue(
                "ctwasPipeline: GwasSumStats for region '{rid}' carries no ",
                "ldSketch (ldSketch = NULL); cTWAS requires an LD reference."
            )
            abort(msg)
        }
        if (!is.null(tw)) {
            .ctwasRequireMatchingLdSketches(getLdSketch(tw), gwasLd)
        }
        ldKey <- .ctwasLdPanelKey(gwasLd)
        if (is.null(ldPanelsByRegion[[ldKey]])) {
            ldPanelsByRegion[[ldKey]] <- .ctwasComputeFullPanelLd(gwasLd)
        }
        ldPanel <- ldPanelsByRegion[[ldKey]]
        ldFileByRegion[[rid]] <- ldKey
        zSnpPieces[[rid]] <- .ctwasBuildZSnp(gss)
        regionInfoPieces[[rid]] <- .ctwasBuildSingleRegionInfo(rid, gss)
        snpMap[[rid]] <- .ctwasSnpInfoForGwasBlock(gss, ldPanel$snpInfo)
    }
    list(
        ldPanelsByRegion = ldPanelsByRegion,
        zSnpPieces = zSnpPieces,
        regionInfoPieces = regionInfoPieces,
        snpMap = snpMap,
        ldFileByRegion = ldFileByRegion
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
    weightsList <- list()
    for (rid in regionIds) {
        tw <- twasWeights[[rid]]
        if (is.null(tw)) {
            next
        }
        twMethod <- .ctwasFilterMethod(tw, resolvedMethod)
        if (is.null(twMethod)) {
            next
        }
        ldPanel <- fp$ldPanelsByRegion[[fp$ldFileByRegion[[rid]]]]
        blockWeights <- .ctwasBuildWeights(
            twMethod,
            ldPanel,
            fineMappingResult = fineMappingResult,
            twasWeightCutoff = cutoffs$twasWeightCutoff,
            csMinCor = cutoffs$csMinCor,
            minPipCutoff = cutoffs$minPipCutoff,
            maxNumVariants = cutoffs$maxNumVariants,
            gwasSnpIds = globalGwasSnpIds,
            regionSnpIds = fp$snpMap[[rid]]$id
        )
        if (length(blockWeights) > 0L) {
            names(blockWeights) <- str_c(rid, "|", names(blockWeights))
            weightsList <- c(weightsList, blockWeights)
        }
    }
    weightsList
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
#' @param thin,niterPrefit,niter Pass-throughs to
#'   \code{ctwas::assemble_region_data} / \code{ctwas::est_param}.
#' @param groupPriorVarStructure Pass-through.
#' @param ncore Number of cores.
#' @param fallbackToPrefit Logical (length 1). When \code{TRUE} (default
#'   \code{FALSE}), if \code{ctwas::est_param}'s accurate EM fails for ANY
#'   reason on a degenerate input, re-run only the prefit step via ctwas's
#'   internal \code{fit_EM} and return those (typically finite) priors as the
#'   param. The accurate-EM failure mode is version-dependent (ctwas <= 0.4.x:
#'   \code{"contains NAs"}; ctwas >= 0.6.0: \code{"No regions selected!"} or a
#'   NaN-loglik \code{"missing value where TRUE/FALSE needed"}), so the catch is
#'   deliberately broad; a genuinely broken input still surfaces because the
#'   prefit re-run will itself error. Mirrors the legacy ctwas_2 workaround on
#'   toy data where the accurate EM cannot be estimated.
#' @param ... Additional arguments forwarded to \code{ctwas::est_param} (e.g.
#'   \code{min_p_single_effect}, \code{min_group_size}).
#' @return The \code{inputs} list augmented with \code{region_data},
#'   \code{boundary_genes}, \code{z_gene}, and \code{param}.
#' @examples
#' data(ctwasInputsExample)
#' estCtwasParam(ctwasInputsExample, thin = 1, niterPrefit = 3,
#'   niter = 10, min_group_size = 1, min_p_single_effect = 0,
#'   fallbackToPrefit = TRUE)
#' @export
estCtwasParam <- function(
    inputs,
    thin = 0.1,
    niterPrefit = 3L,
    niter = 30L,
    groupPriorVarStructure = c(
        "shared_type",
        "shared_context",
        "shared_nonSNP",
        "shared_all",
        "independent"
    ),
    ncore = 1L,
    fallbackToPrefit = FALSE,
    ...
) {
    if (!requireNamespace("ctwas", quietly = TRUE)) {
        # nocov start
        abort("Package 'ctwas' is required for estCtwasParam.")
        # nocov end
    }
    groupPriorVarStructure <- arg_match(groupPriorVarStructure)
    ncore <- as.integer(ncore)
    extra <- list(...)
    zGene <- .ctwasEnsureZGene(inputs, ncore)
    regionData <- .ctwasAssembleRegionData(inputs, zGene, thin, ncore, extra)
    boundaryGenes <- .ctwasBoundaryGenes(inputs, ncore, extra)
    paramRes <- .ctwasEstParamOrFallback(
        regionData,
        niterPrefit,
        niter,
        groupPriorVarStructure,
        ncore,
        thin,
        fallbackToPrefit,
        extra
    )
    # assemble_region_data does not echo z_gene back, so propagate the
    # precomputed z_gene we passed in (inputs$z_gene is NULL when twasZ was not
    # supplied) so $z_gene resolves to the right entry.
    inputs$z_gene <- zGene
    c(
        inputs,
        list(
            region_data = regionData,
            boundary_genes = boundaryGenes,
            param = paramRes
        )
    )
}

# z_gene for assemble_region_data (which requires it non-NULL): use the caller's
# inputs$z_gene, else compute via ctwas::compute_gene_z (mirrors
# ctwas_sumstats).
# @noRd
.ctwasEnsureZGene <- function(inputs, ncore) {
    if (!is.null(inputs$z_gene)) {
        return(inputs$z_gene)
    }
    ctwas::compute_gene_z(inputs$z_snp, inputs$weights, ncore = ncore)
}

# assemble_region_data -> the per-region list (keyed by region_id). It does NOT
# echo z_gene or return boundary genes; both are recovered separately.
# @noRd
.ctwasAssembleRegionData <- function(inputs, zGene, thin, ncore, extra) {
    .ctwasInvoke(
        ctwas::assemble_region_data,
        list(
            region_info = inputs$region_info,
            z_snp = inputs$z_snp,
            z_gene = zGene,
            weights = inputs$weights,
            snp_map = inputs$snp_map,
            thin = thin,
            ncore = ncore
        ),
        extra = extra
    )
}

# Boundary genes (computed internally by assemble_region_data for adjustment but
# never returned) recovered via ctwas::get_boundary_genes; NULL for one region.
# @noRd
.ctwasBoundaryGenes <- function(inputs, ncore, extra) {
    if (nrow(inputs$region_info) <= 1L) {
        return(NULL)
    }
    .ctwasInvoke(
        ctwas::get_boundary_genes,
        list(
            region_info = inputs$region_info,
            weights = inputs$weights,
            ncore = ncore
        ),
        extra = extra
    )
}

# The accurate EM (ctwas::est_param) invocation, split out so
# .ctwasEstParamOrFallback can run it directly (no fallback) or inside a
# tryCatch
# (with fallback) without duplicating the argument assembly.
# @noRd
.ctwasEstParamAccurate <- function(
    regionData,
    niterPrefit,
    niter,
    groupPriorVarStructure,
    ncore,
    extra
) {
    .ctwasInvoke(
        ctwas::est_param,
        list(
            region_data = regionData,
            niter_prefit = as.integer(niterPrefit),
            niter = as.integer(niter),
            group_prior_var_structure = groupPriorVarStructure,
            ncore = ncore
        ),
        extra = extra
    )
}

# est_param (accurate EM), falling back to a converged prefit EM on ANY accurate
# error when fallbackToPrefit is set. The accurate EM fails on degenerate inputs
# in several version-dependent ways (NAs / "No regions selected!" / NaN
# log-likelihood), so catch all rather than match brittle version messages. The
# prefit re-run runs to full `niter` (this is now the final prior).
# @noRd
.ctwasEstParamOrFallback <- function(
    regionData,
    niterPrefit,
    niter,
    groupPriorVarStructure,
    ncore,
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
            ncore,
            extra
        ))
    }
    tryCatch(
        .ctwasEstParamAccurate(
            regionData,
            niterPrefit,
            niter,
            groupPriorVarStructure,
            ncore,
            extra
        ),
        error = function(e) {
            msg <- glue(
                "estCtwasParam: accurate EM unusable ",
                "({conditionMessage(e)}); falling back to prefit estimates."
            )
            inform(msg)
            .ctwasFitPrefitEm(
                regionData,
                niter = as.integer(niter),
                groupPriorVarStructure = groupPriorVarStructure,
                thin = thin,
                ncore = ncore,
                extra = extra
            )
        }
    )
}

#' Screen cTWAS regions
#'
#' @description Step 3 of the three-step \code{\link{ctwasPipeline}}: runs
#'   \code{ctwas::screen_regions} on the \code{\link{estCtwasParam}} result and
#'   returns the screened-region set. Use this entry point to substitute
#'   hand-tuned priors for the ones estimated in step 2 (e.g. when the accurate
#'   EM diverges to NaN and you want to recover the prefit values).
#'
#' @param estResult A list returned by \code{\link{estCtwasParam}}.
#' @param L Unused. Retained for call-site compatibility with
#'   \code{\link{ctwasPipeline}}; ctwas's screening always uses the
#'   single-effect (SER) model and ignores L. \code{L} is applied by
#'   \code{\link{finemapCtwasRegions}} downstream.
#' @param ncore Number of cores.
#' @param ... Additional arguments forwarded to \code{ctwas::screen_regions}
#'   (e.g. \code{min_nonSNP_PIP}, \code{min_snp_pval}, \code{min_var},
#'   \code{min_gene}).
#' @return The \code{estResult} list augmented with \code{screen_res} (the full
#'   ctwas output) and \code{screened_region_data}.
#' @importFrom purrr map compact
#' @examples
#' data(ctwasEstExample)
#' screenCtwasRegions(ctwasEstExample, L = 5L)
#' @export
screenCtwasRegions <- function(estResult, L = 5L, ncore = 1L, ...) {
    if (!requireNamespace("ctwas", quietly = TRUE)) {
        # nocov start
        abort("Package 'ctwas' is required for screenCtwasRegions.")
        # nocov end
    }
    # ctwas::screen_regions requires thin = 1 region_data; expand the
    # thinned set first when assemble_region_data was called with thin < 1
    # (matches ctwas_sumstats's own expand-before-screen step).
    thinVals <- compact(map(estResult$region_data, "thin"))
    needsExpand <- length(thinVals) > 0L && min(unlist(thinVals)) < 1
    regionDataForScreen <- if (needsExpand) {
        .ctwasInvoke(
            ctwas::expand_region_data,
            list(
                region_data = estResult$region_data,
                snp_map = estResult$snp_map,
                z_snp = estResult$z_snp,
                ncore = as.integer(ncore)
            ),
            extra = list(...)
        )
    } else {
        estResult$region_data
    }
    screenRes <- .ctwasInvoke(
        ctwas::screen_regions,
        list(
            region_data = regionDataForScreen,
            group_prior = estResult$param$group_prior,
            group_prior_var = estResult$param$group_prior_var,
            ncore = as.integer(ncore)
        ),
        extra = list(...)
    )
    c(
        estResult,
        list(
            screen_res = screenRes,
            screened_region_data = screenRes$screened_region_data
        )
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
#' @param ncore Number of cores.
#' @param ... Additional arguments forwarded to \code{ctwas::finemap_regions}.
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
#' gwasByRegion <- summaryStatsQc(gss, mafCutoff = 0.0025)
#' inp <- assembleCtwasInputs(gwasSumStats = gwasByRegion,
#'   twasWeights = list(ctwasWeightsExample))
#' est <- estCtwasParam(inp, thin = 1, niterPrefit = 3, niter = 10,
#'   min_group_size = 1, min_p_single_effect = 0, fallbackToPrefit = TRUE)
#' screened <- screenCtwasRegions(est, L = 5L)
#' finemapCtwasRegions(screened, L = 5L)
#' @export
finemapCtwasRegions <- function(screenResult, L = 5L, ncore = 1L, ...) {
    if (!requireNamespace("ctwas", quietly = TRUE)) {
        # nocov start
        abort("Package 'ctwas' is required for finemapCtwasRegions.")
        # nocov end
    }
    rd <- screenResult$screened_region_data
    fmRes <- if (length(rd) == 0L) {
        list(finemap_res = NULL, susie_alpha_res = NULL)
    } else {
        .ctwasInvoke(
            ctwas::finemap_regions,
            list(
                region_data = rd,
                LD_map = screenResult$LD_map,
                weights = screenResult$weights,
                group_prior = screenResult$param$group_prior,
                group_prior_var = screenResult$param$group_prior_var,
                L = as.integer(L),
                LD_format = "custom",
                LD_loader_fun = screenResult$LD_loader_fun,
                snpinfo_loader_fun = screenResult$snpinfo_loader_fun,
                ncore = as.integer(ncore)
            ),
            extra = list(...)
        )
    }
    # Repair cTWAS's molecular_id mislabel (first-"|" split of our composite
    # id).
    fmRes$finemap_res <- .ctwasFixMolecularId(fmRes$finemap_res)
    fmRes$susie_alpha_res <- .ctwasFixMolecularId(fmRes$susie_alpha_res)
    list(
        z_gene = screenResult$z_gene,
        param = screenResult$param,
        finemap_res = fmRes$finemap_res,
        susie_alpha_res = fmRes$susie_alpha_res,
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
#' @param ncore Number of cores. Default \code{1}.
#' @param ... Forwarded to the underlying ctwas postprocess function.
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
    ncore = 1L,
    ...
) {
    # nocov start
    if (!requireNamespace("ctwas", quietly = TRUE)) {
        abort("Package 'ctwas' is required for mergeCtwasBoundaryRegions.")
    }
    # nocov end
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
        ncore
    )
    fa <- .ctwasMergeDispatch(finemapResult, common, L)
    userExtra <- list(...)
    userExtra <- userExtra[setdiff(names(userExtra), names(fa$args))]
    callArgs <- c(fa$args, userExtra)
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
    ncore
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
        ncore = as.integer(ncore)
    )
}

# Pick the LD vs no-LD region-merging fn + args. ctwas's postprocess_*()
# forward `...` into finemap_regions, so the LD loader closures must ride in the
# explicit arg list (not filtered through .ctwasInvoke).
# @noRd
.ctwasMergeDispatch <- function(finemapResult, common, L) {
    if (is.null(finemapResult$LD_loader_fun)) {
        return(list(
            fn = ctwas::postprocess_region_merging_noLD,
            args = common
        ))
    }
    args <- c(
        common,
        list(
            LD_map = finemapResult$LD_map,
            L = as.integer(L),
            LD_format = "custom",
            LD_loader_fun = finemapResult$LD_loader_fun,
            snpinfo_loader_fun = finemapResult$snpinfo_loader_fun
        )
    )
    list(fn = ctwas::postprocess_region_merging, args = args)
}

# Write region-merging outputs back onto the finemap result (only components
# ctwas actually returned).
# @noRd
.ctwasApplyMergeResult <- function(finemapResult, res) {
    finemapResult$finemap_res <- res$updated_finemap_res
    finemapResult$susie_alpha_res <- res$updated_susie_alpha_res
    if (!is.null(res$updated_region_data)) {
        finemapResult$region_data <- res$updated_region_data
    }
    if (!is.null(res$updated_region_info)) {
        finemapResult$region_info <- res$updated_region_info
    }
    if (!is.null(res$updated_LD_map)) {
        finemapResult$LD_map <- res$updated_LD_map
    }
    if (!is.null(res$updated_snp_map)) {
        finemapResult$snp_map <- res$updated_snp_map
    }
    finemapResult$merge_res <- res
    finemapResult
}

# Invoke a ctwas function with a fixed `args` list plus optional `extra`
# (typically the `...` collected by the wrapper). `extra` names that
# duplicate `args` names are silently dropped, so the wrapper's explicit
# arguments always win over caller-supplied `...`.
# @noRd
.ctwasInvoke <- function(fn, args, extra = list()) {
    if (length(extra) > 0L) {
        extra <- extra[setdiff(names(extra), names(args))]
        # `...` is forwarded uniformly to four different ctwas functions
        # (assemble_region_data / est_param / screen_regions /
        # finemap_regions). Restrict to fn's explicit formals so an arg
        # meant for a sibling step doesn't crash this one -- and so args
        # that fn would otherwise forward via its own `...` (e.g. into
        # susie_rss) don't bleed into incompatible downstream functions.
        formalsFn <- tryCatch(names(formals(fn)), error = function(e) NULL)
        if (!is.null(formalsFn)) {
            explicitFormals <- setdiff(formalsFn, "...")
            extra <- extra[intersect(names(extra), explicitFormals)]
        }
        args <- c(args, extra)
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
    ncore,
    extra = list()
) {
    fitEm <- utils::getFromNamespace("fit_EM", "ctwas")
    fitRegionData <- .ctwasPrefitRegionFilter(region_data, extra)
    fitArgs <- list(
        region_data = fitRegionData,
        niter = as.integer(niter),
        group_prior_var_structure = groupPriorVarStructure,
        ncore = as.integer(ncore)
    )
    prefit <- .ctwasInvoke(fitEm, fitArgs, extra)
    adj <- .ctwasApplyThin(prefit$group_prior, prefit$group_size, thin)
    groupSize <- adj$groupSize
    if (length(adj$groupPrior) > 0L) {
        groupSize <- groupSize[names(adj$groupPrior)]
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
    keep <- rep(TRUE, length(region_data))
    if (minVar > 0L) {
        keep <- keep & (nSid + nGid) >= minVar
    }
    if (minGene > 0L) {
        keep <- keep & nGid >= minGene
    }
    fitRegionData <- region_data[keep]
    if (length(fitRegionData) == 0L) {
        abort("No regions selected!")
    }
    fitRegionData
}

# Rescale the SNP group prior / size by `thin` (the SNP subsampling factor).
# @noRd
.ctwasApplyThin <- function(groupPrior, groupSize, thin) {
    if (thin != 1) {
        if (is_in("SNP", names(groupPrior))) {
            groupPrior["SNP"] <- groupPrior["SNP"] * thin
        }
        if (is_in("SNP", names(groupSize))) {
            groupSize["SNP"] <- groupSize["SNP"] / thin
        }
    }
    list(groupPrior = groupPrior, groupSize = groupSize)
}

# =============================================================================
# Internal helpers
# =============================================================================

# LD-sketch identity check. Thin wrapper over the shared
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
    available <- unique(unlist(map(twasWeightsList, .ctwasMethodChr)))
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
.ctwasIsPreBucketed <- function(tw, gwasSumStats) {
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
.ctwasBlockGrFromIds <- function(ids) {
    n <- length(ids)
    chrom <- character(n)
    start <- integer(n)
    end <- integer(n)
    for (i in seq_len(n)) {
        g <- tryCatch(
            asGranges(str_replace(
                as.character(ids[[i]]),
                "_([0-9]+)_([0-9]+)$",
                ":\\1-\\2"
            )),
            error = function(e) NULL
        )
        if (!is.null(g) && length(g) >= 1L) {
            chrom[[i]] <- as.character(GenomicRanges::seqnames(g))[[1L]]
            start[[i]] <- GenomicRanges::start(g)[[1L]]
            end[[i]] <- GenomicRanges::end(g)[[1L]]
        } else {
            chrom[[i]] <- "chrUn"
            start[[i]] <- 1L
            end[[i]] <- 0L
        }
    }
    GenomicRanges::GRanges(chrom, IRanges::IRanges(start = start, end = end))
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
    out <- list()
    for (rid in names(gwasSumStats)) {
        idx <- which(home == rid)
        if (length(idx) == 0L) {
            next
        }
        sub <- combined[idx, ]
        sub@ldSketch <- getLdSketch(gwasSumStats[[rid]])
        out[[rid]] <- sub
    }
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
    if (.ctwasIsPreBucketed(twasWeights, gwasSumStats)) {
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
            unlist(map(twasWeightsList, .ctwasMethodsOf))
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
    studies <- unique(.ctwasStudyChr(gwasSumStats))
    studies <- studies[!is.na(studies) & str_length(studies) > 0L]
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
        df$molecular_id[isGene] <-
            .ctwasParseGeneIds(as.character(df$id)[isGene])$trait
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
        map_chr(coord, .ctwasCoordField, field = "chrom", placed = placed),
        IRanges::IRanges(
            start = map_int(
                coord,
                .ctwasCoordField2,
                field = "start",
                placed = placed
            ),
            end = map_int(
                coord,
                .ctwasCoordField2,
                field = "end",
                placed = placed
            )
        )
    )
    mcols(gr) <- S4Vectors::DataFrame(df, check.names = FALSE)
    gr[placed]
}

# @noRd
.ctwasRowCoord <- function(i, ids, gc) {
    id <- ids[[i]]
    if (!is.null(gc) && is_in(id, names(gc))) {
        return(gc[[id]])
    }
    parsed <- tryCatch(parseVariantId(id), error = function(e) NULL)
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
.ctwasCoordField <- function(co, field, placed) {
    if (is.null(co)) "chrUnplaced" else as.character(co[[field]])
}

# @noRd
.ctwasCoordField2 <- function(co, field, placed) {
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
        param = runResult$param,
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
    mcols(gr) <- S4Vectors::DataFrame(ri, check.names = FALSE)
    gr
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
    rows <- map(
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
    if (keepSnps) {
        snpRow <- .ctwasSnpRow(
            gwasStudy,
            method,
            jointStr,
            fmDf,
            saDf,
            runResult
        )
        if (!is.null(snpRow)) {
            rows <- c(rows, list(snpRow))
        }
    }
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
    s <- unique(as.character(zSnp$study))
    s <- s[!is.na(s) & str_length(s) > 0L]
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

#' @title Structure a granular cTWAS finemap result as a CtwasResult
#' @description Decompose the raw list returned by
#'   \code{\link{finemapCtwasRegions}} (optionally after
#'   \code{\link{mergeCtwasBoundaryRegions}}) into the structured, per-(study,
#'   context) \code{\link{CtwasResult}} -- the same decomposition
#'   \code{\link{ctwasPipeline}} applies to its one-shot output, exposed for the
#'   granular \code{assembleCtwasInputs} \eqn{\to} \code{estCtwasParam}
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
    out <- tw[keep, ]
    out@ldSketch <- getLdSketch(tw)
    out
}

# Build the per-variant Z data.frame ctwas expects from a GwasSumStats.
# Stacks each study row's GRanges via the shared `.entryToSumstatDf`
# helper (R/sumstatsQc.R), then projects to ctwas's column shape and
# bolts on the `study` column ctwas uses to disambiguate stacked rows.
# @noRd
.ctwasBuildZSnp <- function(gwasSumStats) {
    pieces <- list()
    for (i in seq_len(nrow(gwasSumStats))) {
        df <- .entryToSumstatDf(gwasSumStats[[i]], keepChrPrefix = FALSE)
        pieces[[i]] <- tibble(
            id = df$variant_id,
            chrom = as.integer(df$chrom),
            pos = df$pos,
            A1 = df$A1,
            A2 = df$A2,
            z = df$z,
            study = as.character(gwasSumStats$study)[[i]]
        )
    }
    bind_rows(pieces)
}

# Derive the single-row region_info from the LD sketch's snpInfo
# (min/max BP per chromosome). The sketch is assumed to cover exactly
# one block.
# @noRd
.ctwasBuildSingleRegionInfo <- function(regionId, gss) {
    # Derive the block's [start, stop] from the GWAS variants actually in this
    # block (the GwasSumStats entry GRanges) -- NOT the LD sketch. When many
    # blocks share one whole-chromosome LD payload (the common one-file-per-chr
    # layout), getSnpInfo(ldSketch) spans the entire chromosome, so every region
    # would collapse to the same whole-chromosome [start, stop] and every SNP
    # would be assigned to every region (inflating SNP group_size N-fold and
    # diluting the gene prior to ~0).
    pos <- integer(0)
    chrs <- character(0)
    for (i in seq_len(nrow(gss))) {
        gr <- gss[[i]]
        pos <- c(pos, as.integer(GenomicRanges::start(gr)))
        chrs <- c(chrs, as.character(GenomicRanges::seqnames(gr)))
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
    if (length(pos) == 0L) {
        msg <- glue(
            "ctwasPipeline: GwasSumStats block '{regionId}' has no variants ",
            "to define region bounds."
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

# Compute the full-panel LD ONCE and return everything the rest of the
# pipeline needs to consume it. Returns a list with:
#   R        : full-panel correlation matrix (n_var x n_var, dimnames =
#              SNP IDs). Single source of truth for both the per-region
#              LD loader closure and the per-gene R_wgt submatrices.
#   snpInfo  : ctwas-shaped per-block table (chrom, id, pos, alt, ref)
#              -- both the snp_map element and the snpinfo loader return.
#   variance : named numeric vector of per-variant dosage variance from
#              the LD reference. Used to scale non-standardized TWAS
#              weights to the correlation scale that ctwas expects.
# @noRd
.ctwasComputeFullPanelLd <- function(gwasLd) {
    snpInfoCtwas <- .ctwasSnpInfoForBlock(gwasLd)
    geno <- .ldSketchDosage(
        gwasLd,
        seq_len(nrow(snpInfoCtwas)),
        meanImpute = TRUE
    )
    R <- computeLd(geno, method = "sample")
    snpIds <- snpInfoCtwas$id
    dimnames(R) <- list(snpIds, snpIds)
    variance <- set_names(apply(geno, 2, stats::var, na.rm = TRUE), snpIds)
    list(R = R, snpInfo = snpInfoCtwas, variance = variance)
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
    parsed <- tryCatch(parseVariantId(origVids), error = function(e) NULL)
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
    res <- tryCatch(
        harmonizeAlleles(
            targetData = targetDf,
            refVariants = refVariants,
            colToFlip = "w",
            matchMinProp = 0,
            removeUnmatched = TRUE,
            removeStrandAmbiguous = TRUE
        ),
        error = function(e) NULL
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
    alpha <- fits[["alpha"]]
    mu <- fits[["mu"]]
    xCol <- fits[["X_column_scale_factors"]]
    if (is.null(alpha) || is.null(mu) || is.null(xCol)) {
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
    alpha <- as.matrix(alpha)
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
    signFlip <- sign(harmonizedW / origW[keptIdx])
    signFlip[!is.finite(signFlip)] <- 1
    newAlpha <- .ctwasRenormAlpha(alpha, keptIdx)
    if (is.null(newAlpha)) {
        return(NULL)
    }
    muSub <- sweep(mu[, keptIdx, drop = FALSE], 2L, signFlip, `*`)
    xColSub <- xCol[keptIdx]
    # Guard against zero scale factors (shouldn't happen in practice).
    xColSub[xColSub == 0] <- 1
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
# R_wgt is sliced from the cached full-panel LD by SNP ID -- no
# per-gene genotype re-extraction. Variants absent from the panel
# are dropped from that gene's row set.
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
    regionSnpIds = NULL
) {
    panelSnps <- rownames(ldPanel$R)
    # ctwas's compute_gene_z asserts every weight variant exists in the block's
    # z_snp$id. An LD sketch covering more than the block (e.g. a whole-chrom
    # PLINK2) leaks variants outside it, so intersect with the caller's GWAS
    # sumstats variant set when provided.
    #
    # Two variant sets, because they answer different questions. `gwasSnpIds` is
    # the GLOBAL set: it bounds the gene's cis SPAN, which has to cover every
    # block the gene reaches for boundary detection to work. `regionSnpIds` is
    # this block's own set: it bounds the weight vector actually FITTED, because
    # ctwas fine-maps one region at a time.
    if (!is.null(gwasSnpIds)) {
        panelSnps <- intersect(panelSnps, as.character(gwasSnpIds))
    }
    ctx <- list(
        ldPanel = ldPanel,
        panelSnps = panelSnps,
        regionSnps = regionSnpIds,
        refVariants = .ctwasRefVariants(ldPanel$snpInfo),
        fineMappingResult = fineMappingResult,
        cutoffs = list(
            twasWeightCutoff = twasWeightCutoff,
            csMinCor = csMinCor,
            minPipCutoff = minPipCutoff,
            maxNumVariants = maxNumVariants
        )
    )
    genes <- compact(map(
        seq_len(nrow(twasWeights)),
        .ctwasGeneWeight,
        twasWeights = twasWeights,
        ctx = ctx
    ))
    out <- list()
    for (g in genes) {
        out[[g$key]] <- g$entry
    }
    out
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
.ctwasAdjustGeneWeights <- function(parts, aligned, ldPanel) {
    w <- aligned$w
    fits <- .rowFits(parts)
    shrank <- length(aligned$keptIdx) < length(aligned$origVids)
    if (.ctwasIsSusieFit(fits) && shrank) {
        renorm <- .ctwasRenormalizeSusieWeights(
            fits,
            origVids = aligned$origVids,
            origW = aligned$origW,
            keptIdx = aligned$keptIdx,
            harmonizedW = w
        )
        if (!is.null(renorm)) {
            w <- renorm
        }
    }
    if (!.rowStandardized(parts)) {
        varLookup <- ldPanel$variance[aligned$vids]
        if (anyNA(varLookup)) {
            msg <- glue(
                ".ctwasBuildWeights: missing genotype variance for ",
                "{sum(is.na(varLookup))} variant(s) in the LD panel."
            )
            abort(msg)
        }
        w <- w * sqrt(varLookup)
    }
    w
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
.ctwasGeneEntry <- function(vids, w, ldPanel, meta, spanVids = vids) {
    panelInfo <- ldPanel$snpInfo
    rowIdx <- match(spanVids, panelInfo$id)
    list(
        wgt = matrix(w, ncol = 1L, dimnames = list(vids, "wgt")),
        R_wgt = ldPanel$R[vids, vids, drop = FALSE],
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
.ctwasGeneWeight <- function(i, twasWeights, ctx) {
    # Polymorphic: the weight source may be a TwasWeights or a
    # FineMappingResult, so the payload comes from the class-aware bridge.
    parts <- .rowParts(twasWeights, i)
    aligned <- .ctwasAlignGeneWeights(parts, ctx$refVariants, ctx$panelSnps)
    if (is.null(aligned)) {
        return(NULL)
    }
    w <- .ctwasAdjustGeneWeights(parts, aligned, ctx$ldPanel)
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
    list(
        key = meta$key,
        entry = .ctwasGeneEntry(
            kept$vids[inRegion],
            kept$w[inRegion],
            ctx$ldPanel,
            meta,
            spanVids = kept$vids
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
    selectors <- list(study = study, method = method)
    if (is_in("context", .tupleColumnNames(fineMappingResult))) {
        selectors$context <- context
    }
    if (is_in("trait", .tupleColumnNames(fineMappingResult))) {
        selectors$trait <- trait
    }
    selArgs <- c(list(fineMappingResult), selectors)
    entry <- tryCatch(
        exec(getFineMappingResult, !!!selArgs),
        error = function(e) NULL
    )
    if (is.null(entry)) {
        return(NULL)
    }
    tl <- getTopLoci(entry, raw = TRUE)
    if (nrow(tl) == 0L) {
        return(NULL)
    }
    pip <- if (is_in("pip", names(tl))) {
        set_names(as.numeric(tl$pip), as.character(tl$variant_id))
    } else {
        NULL
    }
    cs <- .ctwasCsMembership(tl)
    list(pip = pip, csMembers = cs$csMembers, csPurity = cs$csPurity)
}

# Per-CS membership + purity at 95% coverage from a topLoci table. cs_95 stores
# `<method>_<idx>` where idx == 0 means "not in any CS"; cs_95_purity (when
# present) broadcasts one purity value across a CS's rows.
# @noRd
.ctwasCsMembership <- function(tl) {
    csMembers <- list()
    csPurity <- numeric(0)
    if (!is_in("cs_95", names(tl))) {
        return(list(csMembers = csMembers, csPurity = csPurity))
    }
    csIdx <- suppressWarnings(as.integer(str_remove(tl$cs_95, "^.*_")))
    keepIdx <- !is.na(csIdx) & csIdx > 0L
    for (k in sort(unique(csIdx[keepIdx]))) {
        inCs <- csIdx == k & keepIdx
        csMembers[[length(csMembers) + 1L]] <- as.character(tl$variant_id)[inCs]
        p <- if (is_in("cs_95_purity", names(tl))) {
            as.numeric(tl$cs_95_purity[which(inCs)[1L]])
        } else {
            NA_real_
        }
        csPurity <- c(csPurity, p)
    }
    list(csMembers = csMembers, csPurity = csPurity)
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
        vids <- vids[magKeep]
        w <- w[magKeep]
        if (length(vids) == 0L) {
            return(NULL)
        }
    }
    # Steps 2-3: PIP / CS rescue (only when fineMappingResult was passed).
    mustKeep <- .ctwasMustKeep(vids, finemapAux, csMinCor, minPipCutoff)
    # Step 4: cap, keeping must-keep variants first.
    if (length(vids) > maxNumVariants && is.finite(maxNumVariants)) {
        capped <- .ctwasCapVariants(
            vids,
            w,
            mustKeep,
            finemapAux,
            maxNumVariants
        )
        vids <- capped$vids
        w <- capped$w
    }
    list(vids = vids, w = w)
}

# Variants that must survive the cap: members of any high-purity (>= csMinCor)
# credible set, plus any with PIP > minPipCutoff. Empty when no finemapAux.
# @noRd
.ctwasMustKeep <- function(vids, finemapAux, csMinCor, minPipCutoff) {
    mustKeep <- character(0)
    if (is.null(finemapAux)) {
        return(mustKeep)
    }
    if (length(finemapAux$csMembers) > 0L && csMinCor > 0) {
        for (k in seq_along(finemapAux$csMembers)) {
            if (
                !is.na(finemapAux$csPurity[k]) &&
                    finemapAux$csPurity[k] >= csMinCor
            ) {
                mustKeep <- union(
                    mustKeep,
                    intersect(finemapAux$csMembers[[k]], vids)
                )
            }
        }
    }
    if (!is.null(finemapAux$pip) && minPipCutoff > 0) {
        hits <- names(finemapAux$pip)[finemapAux$pip > minPipCutoff]
        mustKeep <- union(mustKeep, intersect(hits, vids))
    }
    mustKeep
}

# Cap to maxNumVariants: must-keep variants first, then fill by descending PIP
# (falling back to |w| for variants the PIP table doesn't cover).
# @noRd
.ctwasCapVariants <- function(vids, w, mustKeep, finemapAux, maxNumVariants) {
    priorities <- if (!is.null(finemapAux) && !is.null(finemapAux$pip)) {
        unname(finemapAux$pip[vids])
    } else {
        NULL
    }
    if (is.null(priorities) || all(is.na(priorities))) {
        priorities <- abs(w)
    } else {
        priorities[is.na(priorities)] <- abs(w)[is.na(priorities)]
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

# Multi-block LD loader for ctwas. ctwas invokes
# `LD_loader_fun(LD_file)` per region during region_data assembly and
# fine-mapping; we dispatch by `LD_file` (the same string set on
# `LD_map$LD_file`) into the cached per-sketch ldPanel.
# @noRd
.ctwasMultiBlockLdLoader <- function(ldPanelsByRegion) {
    function(LD_file, ...) {
        panel <- ldPanelsByRegion[[LD_file]]
        if (is.null(panel)) {
            msg <- glue(
                "ctwasPipeline LD loader: no cached panel for ",
                "LD_file = '{LD_file}'"
            )
            abort(msg)
        }
        panel$R
    }
}

# Multi-block SNP-info loader for ctwas. Mirrors the LD loader.
# @noRd
.ctwasMultiBlockSnpInfoLoader <- function(ldPanelsByRegion) {
    function(LD_file, ...) {
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

# Build a per-block snpInfo table restricted to variants present in the
# GwasSumStats entry. Mirrors `.ctwasSnpInfoForBlock` but restricts to
# the block's GWAS variants (intersected against the cached panel) so
# snp_map[[region_id]] is sized to the block, not the whole panel.
# @noRd
.ctwasSnpInfoForGwasBlock <- function(gwasSumStats, panelSnpInfo) {
    blockIds <- character(0)
    for (i in seq_len(nrow(gwasSumStats))) {
        mc <- S4Vectors::mcols(gwasSumStats[[i]])
        if (is_in("SNP", colnames(mc))) {
            blockIds <- c(blockIds, as.character(mc$SNP))
        }
    }
    blockIds <- unique(blockIds)
    if (length(blockIds) == 0L) {
        return(panelSnpInfo[FALSE, , drop = FALSE])
    }
    keep <- is_in(panelSnpInfo$id, blockIds)
    panelSnpInfo[keep, , drop = FALSE]
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
