#' @title Causal TWAS Pipeline (cTWAS, multi LD block)
#' @description Pipeline that hands a per-block set of
#'   \code{\link{GwasSumStats}} of GWAS Z-scores together with the
#'   matching per-block per-gene TWAS weights and LD sketches to
#'   \code{ctwas::ctwas_sumstats}, producing per-gene posterior
#'   inclusion probabilities for causal genes. Optionally accepts a
#'   precomputed TWAS-Z \code{GRanges} from
#'   \code{\link{causalInferencePipeline}} as the \code{z_gene} input
#'   so the per-gene Z is not recomputed inside ctwas.
#'
#' @section LD block convention:
#' Inputs are NAMED LISTS keyed by \code{region_id}
#' (\code{list(block1 = gss1, block2 = gss2, ...)}). Per-block
#' \code{region_info}, \code{LD_map}, and \code{snp_map} entries are
#' built automatically from each block's LD sketch and concatenated
#' before the call to \code{ctwas::ctwas_sumstats}. A single-block
#' input is rejected: cTWAS's EM cannot converge on a single region,
#' so callers must supply at least two blocks.
#'
#' @section LD-sketch identity check:
#' Per block: \code{getLdSketch(twasWeights)} (when non-NULL) must
#' match \code{getLdSketch(gwasSumStats)}. Mismatch is a hard error.
#'
#' @param gwasSumStats NAMED LIST of \code{\link{GwasSumStats}} keyed
#'   by \code{region_id} (at least two entries). Each must have
#'   \code{getQcInfo()} non-empty.
#' @param twasWeights The per-gene weight source. Either (a) a FLAT
#'   \code{\link{TwasWeights}} / \code{QtlFineMappingResult} (or a
#'   homogeneous list of them) carrying \code{region} provenance — each
#'   gene is placed into its home LD block internally by
#'   \code{start(region)} (matching cTWAS's \code{p0} assignment rule);
#'   or (b) a pre-bucketed NAMED LIST keyed by \code{region_id} (keys a
#'   SUBSET of \code{gwasSumStats}'s), used as-is. Blocks without any
#'   TWAS weights still contribute their SNP-level signal to ctwas's
#'   joint group prior estimate (the legacy whole-chromosome pattern
#'   where only a few of many LD blocks carry gene weights). A gene whose
#'   cis span straddles a block boundary is homed by its single anchor;
#'   the cross-block signal is cTWAS's boundary-gene concern
#'   (\code{\link{mergeCtwasBoundaryRegions}}), not placement.
#' @param twasZ Optional \code{GRanges} of TWAS Z-scores (output of
#'   \code{\link{causalInferencePipeline}}). When supplied, the
#'   per-(trait, context) Z is used as the \code{z_gene} input to
#'   \code{ctwas_sumstats} so it is not recomputed.
#' @param fineMappingResult Optional \code{QtlFineMappingResult} or
#'   \code{GwasFineMappingResult} carrying the per-variant PIP and
#'   credible-set membership data used by the CS / PIP rescue filters
#'   (\code{csMinCor} and \code{minPipCutoff}). When \code{NULL}
#'   (default) the smart filters are no-ops; only the magnitude filter
#'   (\code{twasWeightCutoff}) and the per-gene cap
#'   (\code{maxNumVariants}, ordered by \code{|weight|}) apply.
#' @param method Optional character (length 1). Picks which TWAS
#'   method's weights to feed into ctwas for each (study, context,
#'   trait) gene. When \code{NULL} (default): use \code{"ensemble"} if
#'   that method is present across the weight sources; otherwise use the
#'   sole method when only one is present; otherwise run \strong{every}
#'   method as an independent cTWAS run (one \code{CtwasResult} row-set
#'   per method). Passing the name explicitly (e.g. \code{"mrash"})
#'   restricts the run to that single method.
#' @param thin,niterPrefit,niter,L Pass-throughs to
#'   \code{ctwas::ctwas_sumstats}.
#' @param groupPriorVarStructure Pass-through (defaults
#'   \code{"shared_type"}).
#' @param ncore Number of cores. Default \code{1}.
#' @param twasWeightCutoff Numeric (length 1). Drop variants with
#'   \code{|weight| < twasWeightCutoff} from each gene's weight matrix
#'   before ctwas sees it. Default \code{0} (no filter).
#' @param csMinCor Numeric (length 1). When \code{fineMappingResult} is
#'   provided, variants belonging to any 95\% credible set with purity
#'   (\code{min_abs_corr}) \code{>= csMinCor} are marked as must-keep
#'   and survive the per-gene cap. Default \code{0.8}. Ignored without
#'   a \code{fineMappingResult}.
#' @param minPipCutoff Numeric (length 1). When
#'   \code{fineMappingResult} is provided, variants with PIP greater
#'   than \code{minPipCutoff} are marked as must-keep and survive the
#'   per-gene cap. Default \code{0} (no PIP rescue). Ignored without a
#'   \code{fineMappingResult}.
#' @param maxNumVariants Numeric (length 1). Cap on per-gene variant
#'   count. When the gene has more variants than this, keep all
#'   must-keep variants and fill remaining slots by descending PIP
#'   (when available) or descending \code{|weight|}. Default
#'   \code{Inf} (no cap).
#' @param fallbackToPrefit Logical (length 1). Forwarded to
#'   \code{\link{estCtwasParam}}. When \code{TRUE}, ctwas's accurate-EM
#'   NaN failure is recovered by falling back to the prefit estimates
#'   (mirrors the legacy ctwas_2 workaround on underpowered data).
#'   Default \code{FALSE}.
#' @param keepSnps Logical (length 1). When \code{TRUE}, retain the
#'   context-agnostic SNP background of each run as one extra
#'   \code{CtwasResult} row (\code{study = context = "SNP"}, mirroring
#'   cTWAS's own \code{"SNP"} group) so the full ctwas output is
#'   reconstructable from \code{\link{getFinemap}} / \code{getSusieAlpha}.
#'   Default \code{FALSE} — the SNP rows are the null background and are
#'   dropped from the structured gene-level result.
#' @param mergeBoundary Logical (length 1). When \code{TRUE}, run
#'   \code{\link{mergeCtwasBoundaryRegions}} after fine-mapping each
#'   run: a high-PIP gene whose cis window straddles an LD-block boundary
#'   has its adjacent regions merged and re-fine-mapped (the legacy
#'   default-off \code{ctwas_3} post-processing). Default \code{FALSE}.
#' @param mergePipThresh Numeric (length 1). PIP threshold for selecting
#'   which boundary genes to merge (\code{\link{mergeCtwasBoundaryRegions}}
#'   \code{pipThresh}). Default \code{0.5}. Ignored unless
#'   \code{mergeBoundary = TRUE}.
#' @param mergeFilterCs Logical (length 1). Require the boundary gene to
#'   be in a credible set to be selected. Default \code{FALSE}. Ignored
#'   unless \code{mergeBoundary = TRUE}.
#' @param mergeMaxSNP Numeric (length 1). Per-merged-region SNP cap.
#'   Default \code{Inf}. Ignored unless \code{mergeBoundary = TRUE}.
#' @param ... Additional arguments forwarded to
#'   \code{ctwas::ctwas_sumstats}.
#' @return A \code{\link{CtwasResult}} collection: one row per
#'   \code{(gwasStudy, study, context, method)}. A single-context run is
#'   one row per method; a multi-context (joint) run emits per-context
#'   rows sharing the same \code{jointContexts} set and the
#'   jointly-estimated group priors. Each row's
#'   \code{\link{CtwasResultEntry}} payload carries that context's
#'   per-gene fine-mapping posteriors (\code{finemap}), the run's
#'   \code{param}, and its \code{regionInfo}. For the raw
#'   \code{ctwas::finemap_regions} list (e.g. to feed
#'   \code{\link{mergeCtwasBoundaryRegions}}), call the granular
#'   \code{\link{assembleCtwasInputs}} \eqn{\to} \code{\link{estCtwasParam}}
#'   \eqn{\to} \code{\link{screenCtwasRegions}} \eqn{\to}
#'   \code{\link{finemapCtwasRegions}} path instead.
#' @export
ctwasPipeline <- function(gwasSumStats,
                          twasWeights,
                          twasZ                   = NULL,
                          fineMappingResult       = NULL,
                          method                  = NULL,
                          thin                    = 0.1,
                          niterPrefit             = 3L,
                          niter                   = 30L,
                          L                       = 5L,
                          groupPriorVarStructure  = c("shared_type",
                                                      "shared_context",
                                                      "shared_nonSNP",
                                                      "shared_all",
                                                      "independent"),
                          ncore                   = 1L,
                          twasWeightCutoff        = 0,
                          csMinCor                = 0.8,
                          minPipCutoff            = 0,
                          maxNumVariants          = Inf,
                          fallbackToPrefit        = FALSE,
                          keepSnps                = FALSE,
                          mergeBoundary           = FALSE,
                          mergePipThresh          = 0.5,
                          mergeFilterCs           = FALSE,
                          mergeMaxSNP             = Inf,
                          ...) {
  groupPriorVarStructure <- match.arg(groupPriorVarStructure)
  .ctwasRequireNamedLists(gwasSumStats, twasWeights)
  methods   <- .ctwasResolveMethods(twasWeights, method)
  gwasStudy <- .ctwasGwasStudy(gwasSumStats)

  # One independent cTWAS run per method (weights are homogeneous within a run),
  # each decomposed into per-context CtwasResult row-specs.
  rows <- list()
  for (m in methods) {
    inputs <- assembleCtwasInputs(
      gwasSumStats       = gwasSumStats,
      twasWeights        = twasWeights,
      twasZ              = twasZ,
      fineMappingResult  = fineMappingResult,
      method             = m,
      twasWeightCutoff   = twasWeightCutoff,
      csMinCor           = csMinCor,
      minPipCutoff       = minPipCutoff,
      maxNumVariants     = maxNumVariants)
    est <- estCtwasParam(
      inputs,
      thin                    = thin,
      niterPrefit             = niterPrefit,
      niter                   = niter,
      groupPriorVarStructure  = groupPriorVarStructure,
      ncore                   = ncore,
      fallbackToPrefit        = fallbackToPrefit,
      ...)
    screened <- screenCtwasRegions(est, L = L, ncore = ncore, ...)
    finemap  <- finemapCtwasRegions(screened, L = L, ncore = ncore, ...)
    # Optional boundary-gene region merging: a gene whose cis window straddles an
    # LD-block boundary is split across two regions in the first-pass fine-map;
    # this merges each high-PIP straddler's adjacent regions and re-fine-maps.
    # The decomposition below is merge-transparent (it keys by gene id).
    if (mergeBoundary)
      finemap <- mergeCtwasBoundaryRegions(
        finemap, pipThresh = mergePipThresh, filterCs = mergeFilterCs,
        maxSNP = mergeMaxSNP, L = L, ncore = ncore, ...)
    rows <- c(rows, .ctwasRunToRows(finemap, gwasStudy = gwasStudy, method = m,
                                    keepSnps = keepSnps))
  }
  if (length(rows) == 0L)
    stop("ctwasPipeline: no genes were modeled (the weight sources produced no ",
         "usable gene weights for method(s): ",
         paste(methods, collapse = ", "), ").")
  .ctwasRowsToResult(rows)
}

#' Assemble cTWAS inputs from S4 GwasSumStats / TwasWeights
#'
#' @description Builds the per-block ctwas-shape input set
#'   (\code{z_snp}, \code{weights}, \code{region_info}, \code{snp_map},
#'   \code{LD_map}, the LD- and SNP-info loader closures, plus optional
#'   \code{z_gene}) that the downstream ctwas steps consume.
#'   This is step 1 of the three-step \code{\link{ctwasPipeline}} split.
#'
#' @details The returned list is the SHARED STATE threaded through
#'   \code{\link{estCtwasParam}} → \code{\link{screenCtwasRegions}} →
#'   \code{\link{finemapCtwasRegions}}. Callers can short-circuit at any
#'   step (e.g. override the estimated priors before fine-mapping) or
#'   call \code{ctwasPipeline()} for the one-shot path.
#'
#' @inheritParams ctwasPipeline
#' @return A list with elements \code{z_snp}, \code{z_gene} (NULL when
#'   no \code{twasZ}), \code{weights}, \code{region_info},
#'   \code{snp_map}, \code{LD_map}, \code{LD_loader_fun},
#'   \code{snpinfo_loader_fun}, and \code{resolvedMethod}.
#' @export
assembleCtwasInputs <- function(gwasSumStats, twasWeights,
                                twasZ              = NULL,
                                fineMappingResult  = NULL,
                                method             = NULL,
                                twasWeightCutoff   = 0,
                                csMinCor           = 0.8,
                                minPipCutoff       = 0,
                                maxNumVariants     = Inf) {
  if (!requireNamespace("ctwas", quietly = TRUE)) {
    # nocov start
    stop("Package 'ctwas' is required for the cTWAS pipeline. ",
         "Install from https://github.com/xinhe-lab/ctwas .")
    # nocov end
  }
  if (missing(gwasSumStats) || !is.list(gwasSumStats) ||
      methods::is(gwasSumStats, "GwasSumStats"))
    stop("`gwasSumStats` must be a NAMED LIST of GwasSumStats keyed by ",
         "region_id (got ", class(gwasSumStats)[[1L]], "). cTWAS's EM ",
         "requires multi-block context to converge; single-block calls ",
         "are no longer supported.")
  if (is.null(names(gwasSumStats)) || any(!nzchar(names(gwasSumStats))))
    stop("`gwasSumStats` must be a named list keyed by region_id (got an ",
         "unnamed or empty-named list).")
  if (missing(twasWeights) || is.null(twasWeights))
    stop("`twasWeights` is required (a TwasWeights / QtlFineMappingResult ",
         "weight source, or a per-region named list keyed by region_id).")
  # Place a FLAT weight source into its home LD blocks by start(region)
  # (cTWAS's p0 rule); a pre-bucketed per-region named list passes through.
  twasWeights <- .ctwasResolveWeightBuckets(twasWeights, gwasSumStats)
  if (is.null(names(twasWeights)) || any(!nzchar(names(twasWeights))))
    stop("`twasWeights` must resolve to a named list keyed by region_id (got an ",
         "unnamed or empty-named list).")
  extra_tw_keys <- setdiff(names(twasWeights), names(gwasSumStats))
  if (length(extra_tw_keys) > 0L)
    stop("`twasWeights` has region_id key(s) not present in ",
         "`gwasSumStats`: ", paste(extra_tw_keys, collapse = ", "))
  if (length(gwasSumStats) < 2L)
    stop("assembleCtwasInputs: at least two LD blocks are required (got ",
         length(gwasSumStats), "). cTWAS's EM cannot estimate the SNP-",
         "group prior variance from a single region.")
  for (rid in names(gwasSumStats)) {
    if (!methods::is(gwasSumStats[[rid]], "GwasSumStats"))
      stop("gwasSumStats[['", rid, "']] is not a GwasSumStats.")
    if (length(getQcInfo(gwasSumStats[[rid]])) == 0L)
      stop("assembleCtwasInputs: gwasSumStats[['", rid,
           "']] has no QC record. Call summaryStatsQc() first.")
  }
  for (rid in names(twasWeights)) {
    if (!methods::is(twasWeights[[rid]], "TwasWeights") &&
        !methods::is(twasWeights[[rid]], "QtlFineMappingResult"))
      stop("twasWeights[['", rid, "']] must be a TwasWeights or ",
           "QtlFineMappingResult (the per-gene weight source).")
  }
  if (!is.null(twasZ) && !methods::is(twasZ, "GRanges"))
    stop("`twasZ` must be a GRanges (output of causalInferencePipeline) ",
         "or NULL.")
  if (!is.null(fineMappingResult) &&
      !methods::is(fineMappingResult, "FineMappingResultBase"))
    stop("`fineMappingResult` must be a FineMappingResultBase ",
         "(QtlFineMappingResult or GwasFineMappingResult) or NULL.")

  regionIds      <- names(gwasSumStats)
  resolvedMethod <- .ctwasResolveMethod(twasWeights, method)

  ldPanelsByRegion <- list()
  weightsList      <- list()
  zSnpPieces       <- list()
  regionInfoPieces <- list()
  snpMap           <- list()
  ldFileByRegion   <- setNames(character(length(regionIds)), regionIds)

  # First pass: cache ld panels, build z_snp pieces, region_info,
  # snp_map per region. We need the union of GWAS variant IDs ACROSS
  # all blocks before we can correctly filter each per-block TwasWeights
  # — a gene's weight variants can straddle adjacent LD blocks, and the
  # per-block GWAS subset would drop the cross-boundary variants.
  for (rid in regionIds) {
    gss    <- gwasSumStats[[rid]]
    tw     <- twasWeights[[rid]]   # may be NULL for SNP-only blocks
    gwasLd <- getLdSketch(gss)
    if (is.null(gwasLd))
      stop("ctwasPipeline: GwasSumStats for region '", rid, "' carries no ",
           "ldSketch (ldSketch = NULL); cTWAS requires an LD reference.")
    if (!is.null(tw))
      .ctwasRequireMatchingLdSketches(getLdSketch(tw), gwasLd)

    ldKey <- .ctwasLdPanelKey(gwasLd)
    if (is.null(ldPanelsByRegion[[ldKey]])) {
      ldPanelsByRegion[[ldKey]] <- .ctwasComputeFullPanelLd(gwasLd)
    }
    ldPanel <- ldPanelsByRegion[[ldKey]]
    ldFileByRegion[[rid]] <- ldKey

    zSnpPieces[[rid]]       <- .ctwasBuildZSnp(gss)
    regionInfoPieces[[rid]] <- .ctwasBuildSingleRegionInfo(rid, gss)
    snpMap[[rid]] <- .ctwasSnpInfoForGwasBlock(gss, ldPanel$snpInfo)
  }

  # Global union of GWAS variant IDs across all blocks. Used to filter
  # per-block TwasWeights so cross-boundary weight variants survive
  # (ctwas's compute_gene_z consumes the concatenated z_snp + the gene's
  # full weight vector regardless of which home block the gene's TSS
  # falls in).
  globalGwasSnpIds <- unique(unlist(lapply(zSnpPieces, function(p) p$id)))

  # Second pass: build per-block weight lists with the GLOBAL gwasSnpIds
  # filter, so a gene whose cis-window straddles block boundaries
  # contributes its full weight vector to ctwas.
  for (rid in regionIds) {
    tw <- twasWeights[[rid]]
    if (is.null(tw)) next
    twMethod <- .ctwasFilterMethod(tw, resolvedMethod)
    if (is.null(twMethod)) next
    ldKey   <- ldFileByRegion[[rid]]
    ldPanel <- ldPanelsByRegion[[ldKey]]
    blockWeights <- .ctwasBuildWeights(
      twMethod, ldPanel,
      fineMappingResult = fineMappingResult,
      twasWeightCutoff  = twasWeightCutoff,
      csMinCor          = csMinCor,
      minPipCutoff      = minPipCutoff,
      maxNumVariants    = maxNumVariants,
      gwasSnpIds        = globalGwasSnpIds)
    if (length(blockWeights) > 0L) {
      names(blockWeights) <- paste0(rid, "|", names(blockWeights))
      weightsList <- c(weightsList, blockWeights)
    }
  }

  zSnp       <- do.call(rbind, zSnpPieces)
  rownames(zSnp) <- NULL
  regionInfo <- do.call(rbind, regionInfoPieces)
  rownames(regionInfo) <- NULL

  ldMap <- data.frame(
    region_id = regionIds,
    LD_file   = unname(ldFileByRegion),
    SNP_file  = unname(ldFileByRegion),
    stringsAsFactors = FALSE)

  list(
    z_snp              = zSnp,
    z_gene             = if (!is.null(twasZ)) .ctwasBuildZGene(twasZ) else NULL,
    weights            = weightsList,
    region_info        = regionInfo,
    snp_map            = snpMap,
    LD_map             = ldMap,
    LD_loader_fun      = .ctwasMultiBlockLdLoader(ldPanelsByRegion),
    snpinfo_loader_fun = .ctwasMultiBlockSnpInfoLoader(ldPanelsByRegion),
    resolvedMethod     = resolvedMethod)
}

#' Estimate cTWAS group prior + prior variance
#'
#' @description Step 2 of the three-step \code{\link{ctwasPipeline}}:
#'   assembles \code{region_data} from the inputs and runs
#'   \code{ctwas::est_param} (prefit EM + accurate EM) to estimate the
#'   group prior probabilities and prior variances. Returns the input
#'   state plus \code{region_data}, \code{boundary_genes},
#'   \code{z_gene}, and \code{param}.
#'
#' @param inputs A list returned by \code{\link{assembleCtwasInputs}}.
#' @param thin,niterPrefit,niter Pass-throughs to
#'   \code{ctwas::assemble_region_data} / \code{ctwas::est_param}.
#' @param groupPriorVarStructure Pass-through.
#' @param ncore Number of cores.
#' @param fallbackToPrefit Logical (length 1). When \code{TRUE} (default
#'   \code{FALSE}), if \code{ctwas::est_param}'s accurate EM fails for ANY
#'   reason on a degenerate input, re-run only the prefit step via
#'   \code{ctwas:::fit_EM} and return those (typically finite) priors as the
#'   param. The accurate-EM failure mode is version-dependent (ctwas <= 0.4.x:
#'   \code{"contains NAs"}; ctwas >= 0.6.0: \code{"No regions selected!"} or a
#'   NaN-loglik \code{"missing value where TRUE/FALSE needed"}), so the catch is
#'   deliberately broad; a genuinely broken input still surfaces because the
#'   prefit re-run will itself error. Mirrors the legacy ctwas_2 workaround on
#'   toy data where the accurate EM cannot be estimated.
#' @param ... Additional arguments forwarded to \code{ctwas::est_param}
#'   (e.g. \code{min_p_single_effect}, \code{min_group_size}).
#' @return The \code{inputs} list augmented with \code{region_data},
#'   \code{boundary_genes}, \code{z_gene}, and \code{param}.
#' @export
estCtwasParam <- function(inputs,
                          thin                    = 0.1,
                          niterPrefit             = 3L,
                          niter                   = 30L,
                          groupPriorVarStructure  = c("shared_type",
                                                      "shared_context",
                                                      "shared_nonSNP",
                                                      "shared_all",
                                                      "independent"),
                          ncore                   = 1L,
                          fallbackToPrefit        = FALSE,
                          ...) {
  if (!requireNamespace("ctwas", quietly = TRUE)) {
    # nocov start
    stop("Package 'ctwas' is required for estCtwasParam.")
    # nocov end
  }
  groupPriorVarStructure <- match.arg(groupPriorVarStructure)
  # ctwas::assemble_region_data assumes z_gene is non-NULL; when the
  # caller did not supply a precomputed twasZ, compute it now via
  # ctwas::compute_gene_z, mirroring ctwas_sumstats's own behaviour.
  zGene <- inputs$z_gene
  if (is.null(zGene)) {
    zGene <- ctwas::compute_gene_z(
      inputs$z_snp, inputs$weights, ncore = as.integer(ncore))
  }
  # ctwas::assemble_region_data returns the region_data list directly
  # (a per-region list, keyed by region_id) — NOT a wrapper carrying
  # $region_data / $boundary_genes. Boundary genes are computed
  # internally for adjustment but never returned, so we recover them
  # separately via the exported ctwas::get_boundary_genes for the
  # downstream finemap return shape.
  regionData <- .ctwasInvoke(ctwas::assemble_region_data, list(
    region_info     = inputs$region_info,
    z_snp           = inputs$z_snp,
    z_gene          = zGene,
    weights         = inputs$weights,
    snp_map         = inputs$snp_map,
    thin            = thin,
    ncore           = as.integer(ncore)), extra = list(...))
  boundaryGenes <- if (nrow(inputs$region_info) > 1L) {
    .ctwasInvoke(ctwas::get_boundary_genes, list(
      region_info = inputs$region_info,
      weights     = inputs$weights,
      ncore       = as.integer(ncore)), extra = list(...))
  } else {
    NULL
  }
  paramRes <- tryCatch(
    .ctwasInvoke(ctwas::est_param, list(
      region_data               = regionData,
      niter_prefit              = as.integer(niterPrefit),
      niter                     = as.integer(niter),
      group_prior_var_structure = groupPriorVarStructure,
      ncore                     = as.integer(ncore)), extra = list(...)),
    error = function(e) {
      # The accurate EM fails on degenerate (e.g. single-gene) inputs in
      # several version-dependent ways: ctwas <= 0.4.x throws "contains NAs";
      # ctwas >= 0.6.0 throws "No regions selected!" (zero regions clear the
      # accurate pass) or "missing value where TRUE/FALSE needed" (NaN
      # log-likelihood in the EM convergence test). Rather than enumerate
      # brittle, version-specific messages, fall back on ANY accurate-EM error
      # when fallbackToPrefit is set: re-run the prefit EM only, which scores
      # every region and skips the p(single effect) selection gate. A genuinely
      # broken input still surfaces, because the prefit re-run will itself error.
      if (fallbackToPrefit) {
        message("estCtwasParam: accurate EM unusable (",
                conditionMessage(e), "); falling back to prefit estimates.")
        .ctwasFitPrefitEm(regionData,
                          niterPrefit            = as.integer(niterPrefit),
                          groupPriorVarStructure = groupPriorVarStructure,
                          thin                   = thin,
                          ncore                  = as.integer(ncore),
                          extra                  = list(...))
      } else {
        stop(e)
      }
    })
  # ctwas::assemble_region_data does not echo z_gene back, so propagate
  # the precomputed (or freshly computed) z_gene we passed into it.
  # Replace inputs$z_gene (which is NULL when twasZ wasn't supplied) so
  # $z_gene resolves to the right entry.
  inputs$z_gene <- zGene
  c(inputs, list(
    region_data    = regionData,
    boundary_genes = boundaryGenes,
    param          = paramRes))
}

#' Screen cTWAS regions
#'
#' @description Step 3 of the three-step \code{\link{ctwasPipeline}}:
#'   runs \code{ctwas::screen_regions} on the
#'   \code{\link{estCtwasParam}} result and returns the screened-region
#'   set. Use this entry point to substitute hand-tuned priors for the
#'   ones estimated in step 2 (e.g. when the accurate EM diverges to
#'   NaN and you want to recover the prefit values).
#'
#' @param estResult A list returned by \code{\link{estCtwasParam}}.
#' @param L Unused. Retained for call-site compatibility with
#'   \code{\link{ctwasPipeline}}; ctwas's screening always uses the
#'   single-effect (SER) model and ignores L. \code{L} is applied by
#'   \code{\link{finemapCtwasRegions}} downstream.
#' @param ncore Number of cores.
#' @param ... Additional arguments forwarded to
#'   \code{ctwas::screen_regions} (e.g. \code{min_nonSNP_PIP},
#'   \code{min_snp_pval}, \code{min_var}, \code{min_gene}).
#' @return The \code{estResult} list augmented with
#'   \code{screen_res} (the full ctwas output) and
#'   \code{screened_region_data}.
#' @export
screenCtwasRegions <- function(estResult,
                               L     = 5L,
                               ncore = 1L,
                               ...) {
  if (!requireNamespace("ctwas", quietly = TRUE)) {
    # nocov start
    stop("Package 'ctwas' is required for screenCtwasRegions.")
    # nocov end
  }
  # ctwas::screen_regions requires thin = 1 region_data; expand the
  # thinned set first when assemble_region_data was called with thin < 1
  # (matches ctwas_sumstats's own expand-before-screen step).
  thinVals <- sapply(estResult$region_data, function(rd) rd$thin)
  thinVals <- thinVals[!sapply(thinVals, is.null)]
  needsExpand <- length(thinVals) > 0L && min(unlist(thinVals)) < 1
  regionDataForScreen <- if (needsExpand) {
    .ctwasInvoke(ctwas::expand_region_data, list(
      region_data = estResult$region_data,
      snp_map     = estResult$snp_map,
      z_snp       = estResult$z_snp,
      ncore       = as.integer(ncore)), extra = list(...))
  } else {
    estResult$region_data
  }
  screenRes <- .ctwasInvoke(ctwas::screen_regions, list(
    region_data        = regionDataForScreen,
    group_prior        = estResult$param$group_prior,
    group_prior_var    = estResult$param$group_prior_var,
    ncore              = as.integer(ncore)), extra = list(...))
  c(estResult, list(
    screen_res            = screenRes,
    screened_region_data  = screenRes$screened_region_data))
}

#' Fine-map cTWAS regions
#'
#' @description Step 4 (final) of the three-step
#'   \code{\link{ctwasPipeline}}: runs \code{ctwas::finemap_regions} on
#'   the screened-region set from \code{\link{screenCtwasRegions}} and
#'   assembles the documented top-level ctwas output (\code{z_gene},
#'   \code{param}, \code{finemap_res}, \code{susie_alpha_res},
#'   \code{region_data}, \code{boundary_genes}, \code{screen_res}).
#'
#' @param screenResult A list returned by
#'   \code{\link{screenCtwasRegions}}.
#' @param L Pass-through.
#' @param ncore Number of cores.
#' @param ... Additional arguments forwarded to
#'   \code{ctwas::finemap_regions}.
#' @return A list mirroring \code{ctwas::ctwas_sumstats}'s output:
#'   \code{z_gene}, \code{param}, \code{finemap_res},
#'   \code{susie_alpha_res}, \code{region_data}, \code{boundary_genes},
#'   \code{screen_res}.
#' @export
finemapCtwasRegions <- function(screenResult,
                                L     = 5L,
                                ncore = 1L,
                                ...) {
  if (!requireNamespace("ctwas", quietly = TRUE)) {
    # nocov start
    stop("Package 'ctwas' is required for finemapCtwasRegions.")
    # nocov end
  }
  rd <- screenResult$screened_region_data
  fmRes <- if (length(rd) == 0L) {
    list(finemap_res = NULL, susie_alpha_res = NULL)
  } else {
    .ctwasInvoke(ctwas::finemap_regions, list(
      region_data        = rd,
      LD_map             = screenResult$LD_map,
      weights            = screenResult$weights,
      group_prior        = screenResult$param$group_prior,
      group_prior_var    = screenResult$param$group_prior_var,
      L                  = as.integer(L),
      LD_format          = "custom",
      LD_loader_fun      = screenResult$LD_loader_fun,
      snpinfo_loader_fun = screenResult$snpinfo_loader_fun,
      ncore              = as.integer(ncore)), extra = list(...))
  }
  # Repair cTWAS's molecular_id mislabel (first-"|" split of our composite id).
  fmRes$finemap_res     <- .ctwasFixMolecularId(fmRes$finemap_res)
  fmRes$susie_alpha_res <- .ctwasFixMolecularId(fmRes$susie_alpha_res)
  list(
    z_gene          = screenResult$z_gene,
    param           = screenResult$param,
    finemap_res     = fmRes$finemap_res,
    susie_alpha_res = fmRes$susie_alpha_res,
    region_data     = screenResult$region_data,
    boundary_genes  = screenResult$boundary_genes,
    screen_res      = screenResult$screen_res,
    # Carried forward so mergeCtwasBoundaryRegions() can re-finemap the merged
    # boundary regions without re-deriving the assembled inputs.
    region_info        = screenResult$region_info,
    z_snp              = screenResult$z_snp,
    weights            = screenResult$weights,
    snp_map            = screenResult$snp_map,
    LD_map             = screenResult$LD_map,
    LD_loader_fun      = screenResult$LD_loader_fun,
    snpinfo_loader_fun = screenResult$snpinfo_loader_fun)
}

#' Merge boundary cTWAS regions and re-fine-map
#'
#' @description Optional step 4 of the cTWAS pipeline (default-off region
#'   merging). A gene whose cis window straddles an LD-block boundary
#'   (a \code{boundary_genes} member) is split across two regions in the
#'   first-pass fine-mapping. This step selects the high-PIP boundary genes,
#'   merges each one's adjacent regions into a single region, re-runs
#'   fine-mapping on the merged regions, and splices the updated results back
#'   into the \code{\link{finemapCtwasRegions}} output. Thin wrapper over
#'   \code{ctwas::postprocess_region_merging()} (or
#'   \code{ctwas::postprocess_region_merging_noLD()} when the inputs carry no
#'   LD loaders).
#'
#' @param finemapResult A list returned by \code{\link{finemapCtwasRegions}}.
#'   Must carry \code{finemap_res}, \code{susie_alpha_res},
#'   \code{region_data}, \code{region_info}, \code{z_snp}, \code{z_gene},
#'   \code{weights}, \code{snp_map}, \code{param}, and — on the LD path —
#'   \code{LD_map} plus the \code{LD_loader_fun} / \code{snpinfo_loader_fun}
#'   closures (all retained by \code{finemapCtwasRegions}).
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
#' @export
mergeCtwasBoundaryRegions <- function(finemapResult,
                                      pipThresh = 0.5,
                                      filterCs  = FALSE,
                                      maxSNP    = Inf,
                                      L         = 5L,
                                      ncore     = 1L,
                                      ...) {
  # nocov start
  if (!requireNamespace("ctwas", quietly = TRUE))
    stop("Package 'ctwas' is required for mergeCtwasBoundaryRegions.")
  # nocov end
  fmRes <- finemapResult$finemap_res
  if (is.null(fmRes) || nrow(fmRes) == 0L) {
    message("mergeCtwasBoundaryRegions: no first-pass finemap result; ",
            "returning unchanged.")
    return(finemapResult)
  }

  hasLd <- !is.null(finemapResult$LD_loader_fun)
  common <- list(
    region_info     = finemapResult$region_info,
    region_data     = finemapResult$region_data,
    z_snp           = finemapResult$z_snp,
    z_gene          = finemapResult$z_gene,
    weights         = finemapResult$weights,
    snp_map         = finemapResult$snp_map,
    finemap_res     = fmRes,
    susie_alpha_res = finemapResult$susie_alpha_res,
    group_prior     = finemapResult$param$group_prior,
    group_prior_var = finemapResult$param$group_prior_var,
    pip_thresh      = pipThresh,
    filter_cs       = filterCs,
    maxSNP          = maxSNP,
    ncore           = as.integer(ncore))

  # ctwas's postprocess_*() forward `...` into finemap_regions, so the LD
  # loader closures must ride in the explicit arg list (not through
  # .ctwasInvoke, which would filter them to postprocess's own formals).
  if (hasLd) {
    fn   <- ctwas::postprocess_region_merging
    args <- c(common, list(
      LD_map             = finemapResult$LD_map,
      L                  = as.integer(L),
      LD_format          = "custom",
      LD_loader_fun      = finemapResult$LD_loader_fun,
      snpinfo_loader_fun = finemapResult$snpinfo_loader_fun))
  } else {
    fn   <- ctwas::postprocess_region_merging_noLD
    args <- common
  }
  userExtra <- list(...)
  userExtra <- userExtra[setdiff(names(userExtra), names(args))]
  res <- do.call(fn, c(args, userExtra))

  finemapResult$finemap_res     <- res$updated_finemap_res
  finemapResult$susie_alpha_res <- res$updated_susie_alpha_res
  if (!is.null(res$updated_region_data)) finemapResult$region_data <- res$updated_region_data
  if (!is.null(res$updated_region_info)) finemapResult$region_info <- res$updated_region_info
  if (!is.null(res$updated_LD_map))      finemapResult$LD_map      <- res$updated_LD_map
  if (!is.null(res$updated_snp_map))     finemapResult$snp_map     <- res$updated_snp_map
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
  do.call(fn, args)
}

# Run ONLY ctwas's prefit EM step against `region_data` and return a
# param list shaped like ctwas::est_param normally produces. Used as
# the fallback path when est_param's accurate EM diverges to NaN on
# toy / underpowered data (matches the legacy ctwas_2 workaround).
# Calls ctwas's internal `fit_EM` (via ::: getFromNamespace) with
# niter = niter_prefit, then applies the same thin-adjustment to the
# SNP group_prior that est_param applies. p_single_effect is left as
# NA since the accurate EM never ran.
# @noRd
.ctwasFitPrefitEm <- function(region_data, niterPrefit,
                              groupPriorVarStructure, thin, ncore,
                              extra = list()) {
  fitEm <- getFromNamespace("fit_EM", "ctwas")
  fitArgs <- list(
    region_data               = region_data,
    niter                     = as.integer(niterPrefit),
    group_prior_var_structure = groupPriorVarStructure,
    ncore                     = as.integer(ncore))
  prefit <- .ctwasInvoke(fitEm, fitArgs, extra)
  groupPrior <- prefit$group_prior
  groupSize  <- prefit$group_size
  if (thin != 1) {
    if ("SNP" %in% names(groupPrior))
      groupPrior["SNP"] <- groupPrior["SNP"] * thin
    if ("SNP" %in% names(groupSize))
      groupSize["SNP"]  <- groupSize["SNP"] / thin
  }
  if (length(groupPrior) > 0L)
    groupSize <- groupSize[names(groupPrior)]
  list(
    group_prior               = groupPrior,
    group_prior_var           = prefit$group_prior_var,
    group_prior_iters         = prefit$group_prior_iters,
    group_prior_var_iters     = prefit$group_prior_var_iters,
    group_prior_var_structure = groupPriorVarStructure,
    group_size                = groupSize,
    p_single_effect           = data.frame(
      region_id        = names(region_data),
      p_single_effect  = NA_real_,
      stringsAsFactors = FALSE))
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
  available <- unique(unlist(lapply(twasWeightsList, function(tw)
    as.character(tw$method))))
  if (length(available) == 0L)
    stop("ctwasPipeline: TwasWeights collections have no method entries.")
  if (!is.null(method) && nzchar(method)) {
    if (!method %in% available)
      stop("ctwasPipeline: method '", method, "' not present in TwasWeights ",
           "(available: ", paste(available, collapse = ", "), ").")
    return(method)
  }
  if ("ensemble" %in% available) return("ensemble")
  if (length(available) == 1L) return(available[[1L]])
  stop("ctwasPipeline: TwasWeights carries multiple methods (",
       paste(available, collapse = ", "),
       ") with no 'ensemble' entry. Supply a `method` argument to ",
       "pick one (e.g. method = \"mrash\").")
}

# Fail fast on the two cTWAS inputs. `gwasSumStats` defines the LD-block grid, so
# it must be a NAMED LIST keyed by region_id (a bare S4 collection is the common
# mistake). `twasWeights` may be a FLAT weight source (a single TwasWeights /
# QtlFineMappingResult, or a list of them) — placed into blocks internally by
# `assembleCtwasInputs` — or a pre-bucketed per-region named list.
# @noRd
.ctwasRequireNamedLists <- function(gwasSumStats, twasWeights) {
  if (!is.list(gwasSumStats) || methods::is(gwasSumStats, "GwasSumStats"))
    stop("`gwasSumStats` must be a NAMED LIST of GwasSumStats keyed by ",
         "region_id (got ", class(gwasSumStats)[[1L]], "). cTWAS's EM ",
         "requires multi-block context to converge; single-block calls ",
         "are no longer supported.")
  okTw <- methods::is(twasWeights, "TwasWeights") ||
          methods::is(twasWeights, "QtlFineMappingResult") || is.list(twasWeights)
  if (!okTw)
    stop("`twasWeights` must be a TwasWeights / QtlFineMappingResult weight ",
         "source (placed into LD blocks internally by region), or a per-region ",
         "named list keyed by region_id (got ", class(twasWeights)[[1L]], ").")
}

# Strip a leading "chr" and case so region_id-derived seqnames ("chr22") and
# phenotype rowRanges seqnames ("22" / "chr22") compare equal for placement.
# @noRd
.ctwasChrKey <- function(x) sub("^chr", "", tolower(as.character(x)))

# TRUE when `tw` is ALREADY a per-region named list (the pre-bucketed contract):
# a NAMED plain list (not an S4 collection) of weight collections. Its keys are
# validated against the block grid downstream (the `extra_tw_keys` check). A
# flat collection or an UNNAMED list is instead treated as a flat weight source
# to place internally by region.
# @noRd
.ctwasIsPreBucketed <- function(tw, gwasSumStats) {
  is.list(tw) && !methods::is(tw, "DFrame") && length(tw) > 0L &&
    !is.null(names(tw)) && all(nzchar(names(tw))) &&
    all(vapply(tw, function(x) methods::is(x, "TwasWeights") ||
                                methods::is(x, "QtlFineMappingResult"), logical(1)))
}

# Combine a flat weight source into ONE collection. Accepts a single
# TwasWeights / QtlFineMappingResult, or a homogeneous list of one kind.
# @noRd
.ctwasCombineWeightSources <- function(weights) {
  if (methods::is(weights, "TwasWeights") ||
      methods::is(weights, "QtlFineMappingResult"))
    return(weights)
  if (is.list(weights)) {
    weights <- weights[!vapply(weights, is.null, logical(1))]
    if (length(weights) == 0L)
      stop("assembleCtwasInputs: `twasWeights` is an empty weight source.")
    if (all(vapply(weights, methods::is, logical(1), "TwasWeights")))
      return(do.call(combineTwasWeights, unname(weights)))
    if (all(vapply(weights, methods::is, logical(1), "QtlFineMappingResult")))
      return(do.call(combineFineMappingResults, unname(weights)))
  }
  stop("assembleCtwasInputs: `twasWeights` must be a TwasWeights or ",
       "QtlFineMappingResult (or a homogeneous list of one kind), or a ",
       "per-region named list keyed by region_id.")
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
  ids     <- names(gwasSumStats)
  blockGr <- .regionFromIds(ids)
  aChr <- .ctwasChrKey(as.character(GenomicRanges::seqnames(region)))
  aPos <- GenomicRanges::start(region)
  bChr <- .ctwasChrKey(as.character(GenomicRanges::seqnames(blockGr)))
  bS   <- GenomicRanges::start(blockGr)
  bE   <- GenomicRanges::end(blockGr)
  vapply(seq_along(aPos), function(i) {
    if (is.na(aPos[[i]])) return(NA_character_)
    hit <- which(bChr == aChr[[i]] & aPos[[i]] >= bS & aPos[[i]] < bE)
    if (length(hit) > 0L) ids[[hit[[1L]]]] else NA_character_
  }, character(1))
}

# Bucket a flat weight source into a per-region named list keyed to the block
# grid, homing each gene by start(region). Each per-block sub-collection carries
# that block's GWAS LD sketch (the panel its weights are harmonized against, and
# what the downstream match-check expects).
# @noRd
.ctwasBucketWeights <- function(weights, gwasSumStats) {
  combined <- .ctwasCombineWeightSources(weights)
  if (!("region" %in% names(combined)))
    stop("assembleCtwasInputs: the weight source carries no `region` ",
         "provenance, which internal LD-block placement requires (produced by ",
         "twasWeightsPipeline / fineMappingPipeline). Supply a pre-bucketed ",
         "per-region named list if placement was done upstream.")
  home     <- .ctwasPlaceByAnchor(getRegion(combined), gwasSumStats)
  unplaced <- sum(is.na(home))
  if (unplaced > 0L)
    warning(sprintf(paste0("assembleCtwasInputs: %d gene(s) whose region anchor ",
                           "fell in no LD block were dropped."), unplaced))
  out <- list()
  for (rid in names(gwasSumStats)) {
    idx <- which(home == rid)
    if (length(idx) == 0L) next
    sub <- combined[idx, ]
    sub@ldSketch <- getLdSketch(gwasSumStats[[rid]])
    out[[rid]] <- sub
  }
  if (length(out) == 0L)
    stop("assembleCtwasInputs: no gene placed into any LD block. Check that the ",
         "weight `region`s and the gwasSumStats region_id keys share a ",
         "coordinate system (e.g. 'chr22_1_1000000').")
  out
}

# Resolve `twasWeights` to a per-region named list: pass a pre-bucketed list
# through, otherwise place a flat weight source by region.
# @noRd
.ctwasResolveWeightBuckets <- function(twasWeights, gwasSumStats) {
  if (.ctwasIsPreBucketed(twasWeights, gwasSumStats)) return(twasWeights)
  .ctwasBucketWeights(twasWeights, gwasSumStats)
}

# Resolve the LIST of TWAS methods a `ctwasPipeline` run should iterate over
# (one independent cTWAS run per method — weights are homogeneous within a run).
#   - explicit `method`: exactly that one (validated present).
#   - NULL + an "ensemble" method present: just "ensemble" (the pre-combined
#     weight, the historical default).
#   - NULL + a single method present: that one.
#   - NULL + MULTIPLE methods, no "ensemble": ALL of them (the singular
#     `.ctwasResolveMethod` errors here; the pipeline instead fans out).
# @noRd
.ctwasResolveMethods <- function(twasWeightsList, method = NULL) {
  methodsOf <- function(tw) if (is.null(tw)) NULL else as.character(tw$method)
  available <- unique(
    if (methods::is(twasWeightsList, "TwasWeights") ||
        methods::is(twasWeightsList, "QtlFineMappingResult"))
      methodsOf(twasWeightsList)                      # a flat weight source
    else unlist(lapply(twasWeightsList, methodsOf)))  # a list of them
  if (length(available) == 0L)
    stop("ctwasPipeline: weight sources carry no method entries.")
  if (!is.null(method) && nzchar(method)) {
    if (!method %in% available)
      stop("ctwasPipeline: method '", method, "' not present in the weight ",
           "sources (available: ", paste(available, collapse = ", "), ").")
    return(method)
  }
  if ("ensemble" %in% available) return("ensemble")
  available   # single -> length-1 (one run); multiple -> iterate over all
}

# The single GWAS (disease) study a cTWAS run models. cTWAS solves one disease
# per run (the z_snp carries a single z per SNP), so the input blocks must all
# reference the same GWAS study.
# @noRd
.ctwasGwasStudy <- function(gwasSumStats) {
  studies <- unique(unlist(lapply(gwasSumStats, function(g)
    as.character(g$study))))
  studies <- studies[!is.na(studies) & nzchar(studies)]
  if (length(studies) == 0L) return(NA_character_)
  if (length(studies) > 1L)
    stop("ctwasPipeline: the input blocks reference multiple GWAS studies (",
         paste(studies, collapse = ", "), "); cTWAS models one disease per run.")
  studies
}

# Parse the cTWAS gene ids (`region|study|context|trait|method`) that name the
# assembled weights list into their identity components. `method` is the LAST
# field and `trait` everything between context and method, so a trait that
# itself contains "|" is preserved.
# @noRd
.ctwasParseGeneIds <- function(ids) {
  parts <- strsplit(ids, "|", fixed = TRUE)
  n <- lengths(parts)
  if (any(n < 5L))
    stop("ctwasPipeline: malformed cTWAS gene id(s): ",
         paste(ids[n < 5L], collapse = ", "),
         " (expected 'region|study|context|trait|method').")
  pick <- function(i) vapply(parts, function(p) p[[i]], character(1))
  data.frame(
    id      = ids,
    rid     = pick(1L),
    study   = pick(2L),
    context = pick(3L),
    trait   = mapply(function(p, k) paste(p[4:(k - 1L)], collapse = "|"),
                     parts, n, USE.NAMES = FALSE),
    method  = mapply(function(p, k) p[[k]], parts, n, USE.NAMES = FALSE),
    stringsAsFactors = FALSE)
}

# cTWAS's finemap_regions derives `molecular_id` by splitting the gene id on the
# FIRST "|", which mislabels our composite `region|study|context|trait|method`
# id (it takes the region as the molecular_id). Restore the true trait for gene
# rows; SNP-background rows (variant ids, no "|") are left untouched.
# @noRd
.ctwasFixMolecularId <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L ||
      !all(c("id", "molecular_id") %in% names(df))) return(df)
  isGene <- lengths(strsplit(as.character(df$id), "|", fixed = TRUE)) >= 5L
  if (any(isGene))
    df$molecular_id[isGene] <-
      .ctwasParseGeneIds(as.character(df$id)[isGene])$trait
  df
}

# Enforce the multi-context joint-model invariant: every context in a run must
# carry the SAME set of genes (traits). A cTWAS joint fit couples the contexts
# through shared group priors; contexts with disjoint gene sets make the joint
# model meaningless. No-op for a single-context run.
# @noRd
.ctwasAssertSharedGenes <- function(parsed) {
  byCtx <- split(parsed$trait, parsed$context)
  if (length(byCtx) < 2L) return(invisible())
  geneSets <- lapply(byCtx, function(g) sort(unique(g)))
  ref      <- geneSets[[1L]]
  mismatch <- names(geneSets)[
    !vapply(geneSets, function(s) identical(s, ref), logical(1))]
  if (length(mismatch) > 0L)
    stop("ctwasPipeline: multi-context cTWAS requires the SAME gene set in ",
         "every context (the joint model is only meaningful when genes are ",
         "shared across contexts). Context(s) differing from the reference '",
         names(geneSets)[[1L]], "' (", length(ref), " gene(s)): ",
         paste(mismatch, collapse = ", "), ".")
  invisible()
}

# Subset a ctwas result frame (finemap_res / susie_alpha_res) to the rows whose
# `id` is in `ids`. Returns NULL when the frame is absent or nothing matches.
# @noRd
.ctwasSubsetById <- function(df, ids) {
  if (is.null(df)) return(NULL)
  sub <- df[as.character(df$id) %in% ids, , drop = FALSE]
  if (nrow(sub) == 0L) NULL else `rownames<-`(sub, NULL)
}

# Subset a ctwas result frame to its SNP rows (type == "SNP"; anno_susie tags
# the non-gene background this way). Returns NULL when absent or none present.
# @noRd
.ctwasSubsetSnp <- function(df) {
  if (is.null(df) || is.null(df$type)) return(NULL)
  sub <- df[as.character(df$type) == "SNP", , drop = FALSE]
  if (nrow(sub) == 0L) NULL else `rownames<-`(sub, NULL)
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
  if (is.null(geneIds) || length(geneIds) == 0L) return(list())
  parsed   <- .ctwasParseGeneIds(geneIds)
  contexts <- unique(parsed$context)
  .ctwasAssertSharedGenes(parsed)
  jointStr <- if (length(contexts) > 1L)
    paste(sort(unique(contexts)), collapse = ",") else NA_character_
  fmDf <- if (is.null(runResult$finemap_res)) NULL else
    as.data.frame(runResult$finemap_res)
  saDf <- if (is.null(runResult$susie_alpha_res)) NULL else
    as.data.frame(runResult$susie_alpha_res)
  mkEntry <- function(fm, sa) CtwasResultEntry(
    finemap = fm, susieAlpha = sa, param = runResult$param,
    regionInfo = runResult$region_info)

  rows <- lapply(contexts, function(cx) {
    inCx    <- parsed$context == cx
    studyCx <- unique(parsed$study[inCx])
    if (length(studyCx) != 1L)
      stop("ctwasPipeline: context '", cx, "' mixes multiple QTL studies (",
           paste(studyCx, collapse = ", "), "); one study per context.")
    idsCx <- parsed$id[inCx]
    list(gwasStudy = gwasStudy, study = studyCx, context = cx, method = method,
         jointContexts = jointStr,
         entry = mkEntry(.ctwasSubsetById(fmDf, idsCx),
                         .ctwasSubsetById(saDf, idsCx)))
  })

  if (keepSnps) {
    snpFm <- .ctwasSubsetSnp(fmDf)
    snpSa <- .ctwasSubsetSnp(saDf)
    if (!is.null(snpFm) || !is.null(snpSa))
      rows <- c(rows, list(list(
        gwasStudy = gwasStudy, study = "SNP", context = "SNP", method = method,
        jointContexts = jointStr, entry = mkEntry(snpFm, snpSa))))
  }
  rows
}

# Assemble accumulated per-run row-specs (from .ctwasRunToRows) into a single
# CtwasResult; the jointContexts column is omitted when every row is
# single-context. Shared by ctwasPipeline (across methods) and asCtwasResult.
# @noRd
.ctwasRowsToResult <- function(rows) {
  if (length(rows) == 0L)
    stop("cTWAS: no genes were modeled (the weight source produced no usable ",
         "gene weights).")
  jointContexts <- vapply(rows, function(r) r$jointContexts, character(1))
  CtwasResult(
    gwasStudy     = vapply(rows, function(r) r$gwasStudy, character(1)),
    study         = vapply(rows, function(r) r$study,     character(1)),
    context       = vapply(rows, function(r) r$context,   character(1)),
    method        = vapply(rows, function(r) r$method,    character(1)),
    entry         = lapply(rows, function(r) r$entry),
    jointContexts = if (any(!is.na(jointContexts))) jointContexts else NULL)
}

# The single GWAS study a finemap result models, read from z_snp$study (which
# `.ctwasBuildZSnp` stamps per row). Errors on multiple; NA when absent.
# @noRd
.ctwasGwasStudyFromZSnp <- function(zSnp) {
  if (is.null(zSnp) || is.null(zSnp$study)) return(NA_character_)
  s <- unique(as.character(zSnp$study))
  s <- s[!is.na(s) & nzchar(s)]
  if (length(s) == 0L) return(NA_character_)
  if (length(s) > 1L)
    stop("asCtwasResult: z_snp references multiple GWAS studies (",
         paste(s, collapse = ", "), "); cTWAS models one disease per run.")
  s
}

# The single weight method a finemap result was built for, read from the
# assembled weight ids. Errors on a mix (the granular path is one method/run).
# @noRd
.ctwasMethodFromWeights <- function(weights) {
  if (is.null(weights) || length(weights) == 0L)
    stop("asCtwasResult: the finemap result carries no weights to derive a ",
         "method from.")
  m <- unique(.ctwasParseGeneIds(names(weights))$method)
  if (length(m) != 1L)
    stop("asCtwasResult: the finemap result mixes weight methods (",
         paste(m, collapse = ", "), "); expected one per run.")
  m
}

#' @title Structure a granular cTWAS finemap result as a CtwasResult
#' @description Decompose the raw list returned by
#'   \code{\link{finemapCtwasRegions}} (optionally after
#'   \code{\link{mergeCtwasBoundaryRegions}}) into the structured, per-(study,
#'   context) \code{\link{CtwasResult}} — the same decomposition
#'   \code{\link{ctwasPipeline}} applies to its one-shot output, exposed for the
#'   granular \code{assembleCtwasInputs} \eqn{\to} \code{estCtwasParam} \eqn{\to}
#'   \code{screenCtwasRegions} \eqn{\to} \code{finemapCtwasRegions} path. The
#'   GWAS study is read from \code{z_snp$study} and the (single) weight method
#'   from the gene ids.
#' @param finemapResult A list from \code{\link{finemapCtwasRegions}} or
#'   \code{\link{mergeCtwasBoundaryRegions}}.
#' @param keepSnps Logical (length 1). Retain the context-agnostic SNP
#'   background as one extra \code{study = context = "SNP"} row. Default
#'   \code{FALSE}. See \code{\link{ctwasPipeline}}.
#' @return A \code{\link{CtwasResult}}.
#' @seealso \code{\link{ctwasPipeline}}, \code{\link{finemapCtwasRegions}}
#' @export
asCtwasResult <- function(finemapResult, keepSnps = FALSE) {
  gwasStudy <- .ctwasGwasStudyFromZSnp(finemapResult$z_snp)
  method    <- .ctwasMethodFromWeights(finemapResult$weights)
  rows <- .ctwasRunToRows(finemapResult, gwasStudy = gwasStudy, method = method,
                          keepSnps = keepSnps)
  .ctwasRowsToResult(rows)
}

# Subset a TwasWeights collection to rows whose `method` matches the
# resolved method. Used to enforce the "one ctwas gene per (study,
# context, trait)" semantics — the legacy pipeline fed a single
# best-CV-method weight per gene; the new S4 TwasWeights may carry
# many methods, but ctwas should only see one.
# @noRd
.ctwasFilterMethod <- function(tw, method) {
  keep <- which(as.character(tw$method) == method)
  if (length(keep) == 0L) return(NULL)
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
    df <- .entryToSumstatDf(gwasSumStats$entry[[i]],
                             keepChrPrefix = FALSE)
    pieces[[i]] <- data.frame(
      id    = df$variant_id,
      chrom = as.integer(df$chrom),
      pos   = df$pos,
      A1    = df$A1,
      A2    = df$A2,
      z     = df$z,
      study = as.character(gwasSumStats$study)[[i]],
      stringsAsFactors = FALSE)
  }
  do.call(rbind, pieces)
}

# Derive the single-row region_info from the LD sketch's snpInfo
# (min/max BP per chromosome). The sketch is assumed to cover exactly
# one block.
# @noRd
.ctwasBuildSingleRegionInfo <- function(regionId, gss) {
  # Derive the block's [start, stop] from the GWAS variants actually in this
  # block (the GwasSumStats entry GRanges) — NOT the LD sketch. When many
  # blocks share one whole-chromosome LD payload (the common one-file-per-chr
  # layout), getSnpInfo(ldSketch) spans the entire chromosome, so every region
  # would collapse to the same whole-chromosome [start, stop] and every SNP
  # would be assigned to every region (inflating SNP group_size N-fold and
  # diluting the gene prior to ~0).
  pos <- integer(0); chrs <- character(0)
  for (i in seq_len(nrow(gss))) {
    gr   <- gss$entry[[i]]
    pos  <- c(pos, as.integer(GenomicRanges::start(gr)))
    chrs <- c(chrs, as.character(GenomicRanges::seqnames(gr)))
  }
  chr <- unique(as.integer(sub("^chr", "", chrs, ignore.case = TRUE)))
  if (length(chr) != 1L)
    stop("ctwasPipeline: GwasSumStats block '", regionId, "' spans multiple ",
         "chromosomes (", paste(chr, collapse = ", "), ").")
  if (length(pos) == 0L)
    stop("ctwasPipeline: GwasSumStats block '", regionId,
         "' has no variants to define region bounds.")
  data.frame(
    region_id = regionId,
    chrom     = chr,
    start     = min(pos),
    stop      = max(pos),
    stringsAsFactors = FALSE)
}

# Per-block SNP info table (chrom, id, pos, alt, ref). ctwas requires
# these exact column names (read_snp_info_files asserts them). `alt`
# maps to A1 (effect allele) and `ref` to A2. NOTE: cTWAS requires an integer
# chrom, so the as.integer() cast below is intentional (an output-boundary
# format requirement) and this path is autosomal-only -- X/Y/MT are not
# supported by the downstream cTWAS model.
# @noRd
.ctwasSnpInfoForBlock <- function(gwasLd) {
  snpInfo <- getSnpInfo(gwasLd)
  chr <- as.integer(sub("^chr", "", as.character(snpInfo$CHR),
                         ignore.case = TRUE))
  data.frame(
    chrom = chr,
    id    = as.character(snpInfo$SNP),
    pos   = as.integer(snpInfo$BP),
    alt   = as.character(snpInfo$A1),
    ref   = as.character(snpInfo$A2),
    stringsAsFactors = FALSE)
}

# Compute the full-panel LD ONCE and return everything the rest of the
# pipeline needs to consume it. Returns a list with:
#   R        : full-panel correlation matrix (n_var x n_var, dimnames =
#              SNP IDs). Single source of truth for both the per-region
#              LD loader closure and the per-gene R_wgt submatrices.
#   snpInfo  : ctwas-shaped per-block table (chrom, id, pos, alt, ref)
#              — both the snp_map element and the snpinfo loader return.
#   variance : named numeric vector of per-variant dosage variance from
#              the LD reference. Used to scale non-standardized TWAS
#              weights to the correlation scale that ctwas expects.
# @noRd
.ctwasComputeFullPanelLd <- function(gwasLd) {
  snpInfoCtwas <- .ctwasSnpInfoForBlock(gwasLd)
  geno  <- .dosageMatrix(gwasLd, seq_len(nrow(snpInfoCtwas)), meanImpute = TRUE)
  R <- computeLd(geno, method = "sample")
  snpIds <- snpInfoCtwas$id
  dimnames(R) <- list(snpIds, snpIds)
  variance <- setNames(apply(geno, 2, stats::var, na.rm = TRUE), snpIds)
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
  if (is.null(parsed) || nrow(parsed) == 0L) return(NULL)
  targetDf <- data.frame(
    chrom   = as.integer(parsed$chrom),
    pos     = as.integer(parsed$pos),
    A2      = as.character(parsed$A2),
    A1      = as.character(parsed$A1),
    w       = as.numeric(origW),
    origIdx = seq_along(origVids),
    stringsAsFactors = FALSE)
  res <- tryCatch(
    harmonizeAlleles(
      targetData            = targetDf,
      refVariants           = refVariants,
      colToFlip             = "w",
      matchMinProp          = 0,
      removeUnmatched       = TRUE,
      removeStrandAmbiguous = TRUE),
    error = function(e) NULL)
  if (is.null(res)) return(NULL)
  res$harmonizedData
}

# Does the entry's `fits` slot carry a SuSiE-shape intermediate (lbf,
# mu, X_column_scale_factors)? Used to gate the renormalization branch.
# @noRd
.ctwasIsSusieFit <- function(fits) {
  if (is.null(fits)) return(FALSE)
  needed <- c("lbf_variable", "mu", "X_column_scale_factors")
  all(needed %in% names(fits))
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
.ctwasRenormalizeSusieWeights <- function(fits, origVids, origW,
                                          keptIdx, harmonizedW) {
  lbf  <- fits$lbf_variable
  mu   <- fits$mu
  xCol <- fits$X_column_scale_factors
  if (is.null(lbf) || is.null(mu) || is.null(xCol)) return(NULL)
  if (ncol(lbf) != length(origVids) ||
      ncol(mu)  != length(origVids) ||
      length(xCol) != length(origVids)) {
    # Fit-vs-entry dimension mismatch; skip rather than mis-slice.
    return(NULL)
  }
  # Per-variant sign flip applied by allele harmonization. NaN signs
  # (origW == 0) default to +1.
  signFlip <- sign(harmonizedW / origW[keptIdx])
  signFlip[!is.finite(signFlip)] <- 1
  lbfSub  <- lbf[, keptIdx, drop = FALSE]
  muSub   <- sweep(mu[, keptIdx, drop = FALSE], 2L, signFlip, `*`)
  xColSub <- xCol[keptIdx]
  # Guard against zero scale factors (shouldn't happen in practice).
  xColSub[xColSub == 0] <- 1
  newAlpha <- lbfToAlpha(lbfSub)
  as.numeric(colSums(newAlpha * muSub) / xColSub)
}

# Build the weights list ctwas expects: keyed by per-tuple gene id,
# each element a list with wgt (variants x 1 matrix; rownames = SNP id),
# R_wgt (per-gene LD submatrix), and gene metadata. ctwas's compute_gene_z
# pulls rownames(wgt) for the SNP IDs and computes z.gene = crossprod(wgt,
# z.s) / sqrt(t(wgt) %*% R_wgt %*% wgt), so wgt must be a numeric matrix
# (not a vector) and R_wgt must be the LD submatrix over the same SNPs.
#
# R_wgt is sliced from the cached full-panel LD by SNP ID — no
# per-gene genotype re-extraction. Variants absent from the panel
# are dropped from that gene's row set.
# @noRd
.ctwasBuildWeights <- function(twasWeights, ldPanel,
                               fineMappingResult = NULL,
                               twasWeightCutoff  = 0,
                               csMinCor          = 0.8,
                               minPipCutoff      = 0,
                               maxNumVariants    = Inf,
                               gwasSnpIds        = NULL) {
  panelSnps <- rownames(ldPanel$R)
  # ctwas's compute_gene_z asserts that every weight variant exists in the
  # block's z_snp$id. When the LD sketch covers more than the block (e.g.
  # a whole-chromosome PLINK2 used for a single-block GwasSumStats),
  # panelSnps alone leaks variants outside the block. Intersect with the
  # caller-supplied GWAS sumstats variant set when provided.
  if (!is.null(gwasSnpIds)) {
    panelSnps <- intersect(panelSnps, as.character(gwasSnpIds))
  }
  panelInfo <- ldPanel$snpInfo
  # Reference frame for allele-harmonization: panel variant info with
  # the column shape `harmonizeAlleles` expects (chrom/pos/A2/A1).
  refVariants <- data.frame(
    chrom      = as.integer(panelInfo$chrom),
    pos        = as.integer(panelInfo$pos),
    A2         = as.character(panelInfo$ref),
    A1         = as.character(panelInfo$alt),
    variant_id = as.character(panelInfo$id),
    stringsAsFactors = FALSE)

  out <- list()
  for (i in seq_len(nrow(twasWeights))) {
    entry    <- twasWeights$entry[[i]]
    # Shared resolver: aligned (variantIds, weights) from either a
    # TwasWeightsEntry or a FineMappingEntry (topLoci posterior effect).
    wr       <- resolveWeights(entry)
    origVids <- wr$variantIds
    origW    <- wr$weights
    if (length(origVids) == 0L) next

    # --- Step 1: allele-harmonize against the LD panel -------------
    # Parses chr:pos:A2:A1 IDs into the data.frame `harmonizeAlleles`
    # expects, then matches by (chrom, pos) with exact / sign-flip /
    # strand-flip detection. Returned canonical variant IDs are in the
    # panel's A1/A2 frame; weights are sign-flipped for variants whose
    # input A1/A2 frame was swapped relative to the panel.
    harm <- .ctwasHarmonizeWeights(origVids, origW, refVariants)
    if (is.null(harm) || nrow(harm) == 0L) next
    vids    <- as.character(harm$variant_id)
    w       <- as.numeric(harm$w)
    keptIdx <- as.integer(harm$origIdx)  # back-reference into origVids/origW

    # --- Step 2: restrict to panel ∩ gwasSnpIds --------------------
    keep <- vids %in% panelSnps
    if (!any(keep)) next
    vids    <- vids[keep]
    w       <- w[keep]
    keptIdx <- keptIdx[keep]

    # --- Step 3: SuSiE alpha renormalization -----------------------
    # When the entry carries a SuSiE-style fit (lbf_variable + mu +
    # X_column_scale_factors) and the kept variant set is smaller than
    # the original fit, the posterior probabilities `alpha` no longer
    # sum to 1 over the kept variants. Renormalize via softmax of
    # lbf_variable over the kept columns and recompute the per-variant
    # weight as colSums(new_alpha * mu_subset) /
    # X_column_scale_factors_subset. mu is sign-flipped per the allele
    # harmonization so the recomputed weight stays in the panel's
    # allele frame. Mirrors the legacy adjustSusieWeights helper.
    fits <- getFits(entry)
    if (.ctwasIsSusieFit(fits) && length(keptIdx) < length(origVids)) {
      renorm <- .ctwasRenormalizeSusieWeights(
        fits, origVids = origVids, origW = origW,
        keptIdx = keptIdx, harmonizedW = w)
      if (!is.null(renorm)) w <- renorm
    }

    # --- Step 4: variance scaling for non-standardized weights -----
    # w_scaled = w_raw * sqrt(per-variant genotype variance from the
    # LD reference panel). Standardized entries (RSS-style, already on
    # the correlation scale) pass through unchanged.
    if (!isTRUE(getStandardized(entry))) {
      varLookup <- ldPanel$variance[vids]
      if (anyNA(varLookup))
        stop(".ctwasBuildWeights: missing genotype variance for ",
             sum(is.na(varLookup)), " variant(s) in the LD panel.")
      w <- w * sqrt(varLookup)
    }

    gStudy   <- as.character(twasWeights$study)[[i]]
    gContext <- as.character(twasWeights$context)[[i]]
    gTrait   <- as.character(twasWeights$trait)[[i]]
    gMethod  <- as.character(twasWeights$method)[[i]]
    key <- sprintf("%s|%s|%s|%s", gStudy, gContext, gTrait, gMethod)

    # PIP / credible-set context for the smart filters (csMinCor +
    # minPipCutoff). Only available when the caller passed the matching
    # FineMappingResult. NULL means we fall back to weight-magnitude
    # priority only.
    finemapAux <- .ctwasGetFinemapAux(fineMappingResult, gStudy, gContext,
                                       gTrait, gMethod)

    # Apply the four filters in order.
    kept <- .ctwasFilterVariants(
      vids = vids, w = w, finemapAux = finemapAux,
      twasWeightCutoff = twasWeightCutoff,
      csMinCor         = csMinCor,
      minPipCutoff     = minPipCutoff,
      maxNumVariants   = maxNumVariants)
    if (length(kept) < 1L) next
    vids <- kept$vids; w <- kept$w

    Rwgt   <- ldPanel$R[vids, vids, drop = FALSE]
    wgtMat <- matrix(w, ncol = 1L, dimnames = list(vids, "wgt"))

    # Per-gene chromosome + BP span derived from the cached snpInfo
    # AFTER filtering (so p0/p1 reflect the retained variants).
    rowIdx <- match(vids, panelInfo$id)
    gChrom <- as.integer(panelInfo$chrom[[rowIdx[1L]]])
    gP0 <- min(as.integer(panelInfo$pos[rowIdx]))
    gP1 <- max(as.integer(panelInfo$pos[rowIdx]))

    out[[key]] <- list(
      wgt          = wgtMat,
      R_wgt        = Rwgt,
      type         = gContext,
      context      = gContext,
      gene_name    = gTrait,
      study        = gStudy,
      method       = gMethod,
      n_wgt        = length(vids),
      chrom        = gChrom,
      p0           = gP0,
      p1           = gP1,
      molecular_id = gTrait,
      weight_name  = paste(gContext, gContext, sep = "_"))
  }
  out
}

# Look up the per-(study, context, trait, method) PIP vector and the
# 95% credible-set membership / purity for one gene from the supplied
# FineMappingResult. Returns NULL when no FineMappingResult was passed
# or no matching tuple exists. Output is a list with:
#   pip       : named numeric vector keyed by variant_id
#   csMembers : list of character vectors (one per CS at 95% coverage)
#   csPurity  : numeric vector aligned with csMembers
# @noRd
.ctwasGetFinemapAux <- function(fineMappingResult, study, context, trait,
                                method) {
  if (is.null(fineMappingResult)) return(NULL)
  selectors <- list(study = study, method = method)
  if ("context" %in% names(fineMappingResult)) selectors$context <- context
  if ("trait"   %in% names(fineMappingResult)) selectors$trait   <- trait
  entry <- tryCatch(
    do.call(getFineMappingResult,
            c(list(fineMappingResult), selectors)),
    error = function(e) NULL)
  if (is.null(entry)) return(NULL)
  tl <- entry@topLoci
  if (nrow(tl) == 0L) return(NULL)
  pip <- if ("pip" %in% names(tl))
            setNames(as.numeric(tl$pip), as.character(tl$variant_id))
         else NULL
  # Per-CS membership at 95% coverage. cs_95 stores `<method>_<idx>`
  # where idx == 0 means "not in any CS".
  csMembers <- list(); csPurity <- numeric(0)
  if ("cs_95" %in% names(tl)) {
    csIdx <- suppressWarnings(as.integer(sub("^.*_", "", tl$cs_95)))
    keepIdx <- !is.na(csIdx) & csIdx > 0L
    for (k in sort(unique(csIdx[keepIdx]))) {
      members <- as.character(tl$variant_id)[csIdx == k & keepIdx]
      csMembers[[length(csMembers) + 1L]] <- members
      # Pull the purity from cs_95_purity if present; same value
      # broadcast to every row in the CS, so any row will do.
      p <- if ("cs_95_purity" %in% names(tl))
              as.numeric(tl$cs_95_purity[which(csIdx == k & keepIdx)[1L]])
           else NA_real_
      csPurity <- c(csPurity, p)
    }
  }
  list(pip = pip, csMembers = csMembers, csPurity = csPurity)
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
.ctwasFilterVariants <- function(vids, w, finemapAux,
                                 twasWeightCutoff, csMinCor,
                                 minPipCutoff, maxNumVariants) {
  if (length(vids) == 0L) return(NULL)
  # Step 1: magnitude.
  if (twasWeightCutoff > 0) {
    magKeep <- !is.na(w) & abs(w) >= twasWeightCutoff
    vids <- vids[magKeep]; w <- w[magKeep]
    if (length(vids) == 0L) return(NULL)
  }
  # Steps 2-3: PIP / CS rescue (only when fineMappingResult was passed).
  mustKeep <- character(0)
  if (!is.null(finemapAux)) {
    if (length(finemapAux$csMembers) > 0L && csMinCor > 0) {
      for (k in seq_along(finemapAux$csMembers)) {
        if (!is.na(finemapAux$csPurity[k]) &&
            finemapAux$csPurity[k] >= csMinCor) {
          mustKeep <- union(mustKeep,
                            intersect(finemapAux$csMembers[[k]], vids))
        }
      }
    }
    if (!is.null(finemapAux$pip) && minPipCutoff > 0) {
      hits <- names(finemapAux$pip)[finemapAux$pip > minPipCutoff]
      mustKeep <- union(mustKeep, intersect(hits, vids))
    }
  }
  # Step 4: cap. Always keep must-keep variants; fill the rest by
  # descending PIP (when PIP available) or descending |w|.
  if (length(vids) > maxNumVariants && is.finite(maxNumVariants)) {
    priorities <- if (!is.null(finemapAux) && !is.null(finemapAux$pip)) {
      unname(finemapAux$pip[vids])  # NAs for variants without PIP
    } else NULL
    if (is.null(priorities) || all(is.na(priorities))) {
      priorities <- abs(w)
    } else {
      # Fall back to |w| for variants the PIP table doesn't know about.
      priorities[is.na(priorities)] <- abs(w)[is.na(priorities)]
    }
    # Order: must-keep first, then the rest by descending priority.
    isMust <- vids %in% mustKeep
    ord <- order(!isMust, -priorities)
    keepIdx <- ord[seq_len(min(maxNumVariants, length(vids)))]
    vids <- vids[keepIdx]; w <- w[keepIdx]
  }
  list(vids = vids, w = w)
}

# Build z_gene data.frame from a TWAS-Z GRanges (output of
# causalInferencePipeline). One row per (qtlStudy, context, trait,
# method, gwasStudy) tuple.
# @noRd
.ctwasBuildZGene <- function(twasZ) {
  mc <- as.data.frame(S4Vectors::mcols(twasZ))
  data.frame(
    id        = sprintf("%s|%s|%s|%s",
                        mc$qtlStudy, mc$context,
                        mc$trait, mc$method),
    z         = as.numeric(mc$twasZ),
    type      = as.character(mc$context),
    context   = as.character(mc$context),
    gene_name = as.character(mc$trait),
    study     = as.character(mc$qtlStudy),
    method    = as.character(mc$method),
    stringsAsFactors = FALSE)
}

# Multi-block LD loader for ctwas. ctwas invokes
# `LD_loader_fun(LD_file)` per region during region_data assembly and
# fine-mapping; we dispatch by `LD_file` (the same string set on
# `LD_map$LD_file`) into the cached per-sketch ldPanel.
# @noRd
.ctwasMultiBlockLdLoader <- function(ldPanelsByRegion) {
  function(LD_file, ...) {
    panel <- ldPanelsByRegion[[LD_file]]
    if (is.null(panel))
      stop("ctwasPipeline LD loader: no cached panel for LD_file = '",
           LD_file, "'")
    panel$R
  }
}

# Multi-block SNP-info loader for ctwas. Mirrors the LD loader.
# @noRd
.ctwasMultiBlockSnpInfoLoader <- function(ldPanelsByRegion) {
  function(LD_file, ...) {
    panel <- ldPanelsByRegion[[LD_file]]
    if (is.null(panel))
      stop("ctwasPipeline snpInfo loader: no cached panel for LD_file = '",
           LD_file, "'")
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
.ctwasLdPanelKey <- function(handle) {
  fmt <- getFormat(handle)
  stem <- .genotypeReadPath(handle)
  candidates <- switch(fmt,
    "plink2" = c(paste0(stem, ".pgen")),
    "plink1" = c(paste0(stem, ".bed")),
    "gds"    = c(stem),
    "vcf"    = c(stem),
    stem)
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L)
    stop("ctwasPipeline: could not derive an existing LD-file token for ",
         "the GenotypeHandle (format=", fmt, ", path=", stem,
         "). Looked for: ", paste(candidates, collapse = ", "))
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
    mc <- S4Vectors::mcols(gwasSumStats$entry[[i]])
    if ("SNP" %in% colnames(mc))
      blockIds <- c(blockIds, as.character(mc$SNP))
  }
  blockIds <- unique(blockIds)
  if (length(blockIds) == 0L) return(panelSnpInfo[FALSE, , drop = FALSE])
  keep <- panelSnpInfo$id %in% blockIds
  panelSnpInfo[keep, , drop = FALSE]
}
