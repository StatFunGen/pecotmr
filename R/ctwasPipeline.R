#' @title Causal TWAS Pipeline (cTWAS, single LD block)
#' @description Per-LD-block pipeline that hands a
#'   \code{\link{GwasSumStats}} of GWAS Z-scores together with per-gene
#'   TWAS weights and the shared LD sketch to
#'   \code{ctwas::ctwas_sumstats}, producing per-gene posterior
#'   inclusion probabilities for causal genes. Optionally accepts a
#'   precomputed TWAS-Z \code{GRanges} from
#'   \code{\link{causalInferencePipeline}} as the \code{z_gene} input
#'   so the per-gene Z is not recomputed inside ctwas.
#'
#' @section LD block convention:
#' Each call assumes the inputs cover exactly one LD block — the user
#' is responsible for constructing the \code{GwasSumStats} and
#' \code{TwasWeights} over the block of interest before calling this
#' pipeline (the same convention used by
#' \code{\link{fineMappingPipeline}} on \code{GwasSumStats}). The
#' single-region \code{region_info}, \code{LD_map}, and \code{snp_map}
#' that \code{ctwas::ctwas_sumstats} requires are derived
#' automatically from the LD sketch on \code{gwasSumStats}.
#'
#' @section LD-sketch identity check:
#' \code{getLdSketch(twasWeights)} (when non-NULL) must match
#' \code{getLdSketch(gwasSumStats)}. Mismatch is a hard error.
#'
#' @param gwasSumStats A \code{\link{GwasSumStats}} over one LD block.
#'   Must have \code{getQcInfo()} non-empty.
#' @param twasWeights A \code{\link{TwasWeights}} carrying per-(study,
#'   context, trait, method) weights over the same LD block.
#' @param twasZ Optional \code{GRanges} of TWAS Z-scores (output of
#'   \code{\link{causalInferencePipeline}}). When supplied, the
#'   per-(trait, context) Z is used as the \code{z_gene} input to
#'   \code{ctwas_sumstats} so it is not recomputed.
#' @param regionId Optional character (length 1) label for the LD
#'   block. Default \code{"block1"}.
#' @param thin,niterPrefit,niter,L Pass-throughs to
#'   \code{ctwas::ctwas_sumstats}.
#' @param groupPriorVarStructure Pass-through (defaults
#'   \code{"shared_type"}).
#' @param ncore Number of cores. Default \code{1}.
#' @param ... Additional arguments forwarded to
#'   \code{ctwas::ctwas_sumstats}.
#' @return Whatever \code{ctwas::ctwas_sumstats} returns (a list with
#'   \code{susie_alpha_res}, \code{param}, and other diagnostics).
#' @export
ctwasPipeline <- function(gwasSumStats,
                          twasWeights,
                          twasZ                   = NULL,
                          regionId                = "block1",
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
                          ...) {
  if (!requireNamespace("ctwas", quietly = TRUE)) {
    stop("Package 'ctwas' is required for ctwasPipeline. ",
         "Install from https://github.com/xinhe-lab/ctwas .")
  }
  if (!methods::is(gwasSumStats, "GwasSumStats"))
    stop("`gwasSumStats` must be a GwasSumStats object.")
  if (length(getQcInfo(gwasSumStats)) == 0L)
    stop("ctwasPipeline: gwasSumStats has no QC record. Call ",
         "summaryStatsQc() first.")
  if (missing(twasWeights) || !methods::is(twasWeights, "TwasWeights"))
    stop("`twasWeights` must be a TwasWeights object.")
  if (!is.null(twasZ) && !methods::is(twasZ, "GRanges"))
    stop("`twasZ` must be a GRanges (output of causalInferencePipeline) ",
         "or NULL.")
  if (length(regionId) != 1L || !nzchar(regionId))
    stop("`regionId` must be a single non-empty character string.")
  groupPriorVarStructure <- match.arg(groupPriorVarStructure)

  twLd   <- getLdSketch(twasWeights)
  gwasLd <- getLdSketch(gwasSumStats)
  .ctwasRequireMatchingLdSketches(twLd, gwasLd)

  # --- Build the single-region ctwas inputs ---------------------------
  zSnp        <- .ctwasBuildZSnp(gwasSumStats)
  regionInfo  <- .ctwasBuildSingleRegionInfo(regionId, gwasLd)
  ldMap       <- data.frame(region_id = regionId, LD_file = regionId,
                            stringsAsFactors = FALSE)
  snpMap      <- list()
  snpMap[[regionId]] <- .ctwasSnpInfoForBlock(gwasLd)
  weightsList <- .ctwasBuildWeights(twasWeights)
  zGene       <- if (!is.null(twasZ)) .ctwasBuildZGene(twasZ) else NULL

  # --- Call the ctwas engine ------------------------------------------
  ctwas::ctwas_sumstats(
    z_snp                      = zSnp,
    weights                    = weightsList,
    region_info                = regionInfo,
    LD_map                     = ldMap,
    snp_map                    = snpMap,
    z_gene                     = zGene,
    thin                       = thin,
    niter_prefit               = as.integer(niterPrefit),
    niter                      = as.integer(niter),
    L                          = as.integer(L),
    group_prior_var_structure  = groupPriorVarStructure,
    LD_format                  = "custom",
    LD_loader_fun              = .ctwasSingleBlockLdLoader(gwasLd),
    snpinfo_loader_fun         = .ctwasSingleBlockSnpInfoLoader(gwasLd),
    ncore                      = as.integer(ncore),
    ...)
}

# =============================================================================
# Internal helpers
# =============================================================================

# LD-sketch identity check. Thin wrapper over the shared
# `.requireMatchingLdSketches` helper (R/ld.R).
.ctwasRequireMatchingLdSketches <- function(twLd, gwasLd) {
  .requireMatchingLdSketches(twLd, gwasLd, pipelineName = "ctwasPipeline")
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
.ctwasBuildSingleRegionInfo <- function(regionId, gwasLd) {
  snpInfo <- getSnpInfo(gwasLd)
  chr <- unique(as.integer(sub("^chr", "", as.character(snpInfo$CHR),
                                ignore.case = TRUE)))
  if (length(chr) != 1L)
    stop("ctwasPipeline: gwasSumStats LD sketch spans multiple ",
         "chromosomes (", paste(chr, collapse = ", "),
         "). ctwasPipeline assumes a single LD block per call.")
  data.frame(
    region_id = regionId,
    chrom     = chr,
    start     = min(as.integer(snpInfo$BP)),
    stop      = max(as.integer(snpInfo$BP)),
    stringsAsFactors = FALSE)
}

# Per-block SNP info table (id, chrom, pos, A1, A2). Used both as the
# single snp_map element and as the loader return value.
# @noRd
.ctwasSnpInfoForBlock <- function(gwasLd) {
  snpInfo <- getSnpInfo(gwasLd)
  chr <- as.integer(sub("^chr", "", as.character(snpInfo$CHR),
                         ignore.case = TRUE))
  data.frame(
    id    = as.character(snpInfo$SNP),
    chrom = chr,
    pos   = as.integer(snpInfo$BP),
    A1    = as.character(snpInfo$A1),
    A2    = as.character(snpInfo$A2),
    stringsAsFactors = FALSE)
}

# Build the weights list ctwas expects: keyed by per-tuple gene id,
# each element a list with id (SNP id), wgt (weight vector), and gene
# metadata. Walks every TwasWeights row.
# @noRd
.ctwasBuildWeights <- function(twasWeights) {
  out <- list()
  for (i in seq_len(nrow(twasWeights))) {
    entry  <- twasWeights$entry[[i]]
    vids   <- getVariantIds(entry)
    w      <- as.numeric(getWeights(entry))
    if (length(vids) == 0L || length(vids) != length(w)) next
    gStudy   <- as.character(twasWeights$study)[[i]]
    gContext <- as.character(twasWeights$context)[[i]]
    gTrait   <- as.character(twasWeights$trait)[[i]]
    gMethod  <- as.character(twasWeights$method)[[i]]
    key <- sprintf("%s|%s|%s|%s", gStudy, gContext, gTrait, gMethod)
    out[[key]] <- list(
      id        = vids,
      wgt       = w,
      type      = gContext,
      context   = gContext,
      gene_name = gTrait,
      study     = gStudy,
      method    = gMethod)
  }
  out
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

# Single-block LD loader for ctwas: ignores the LD_file token and just
# materialises the dosage matrix from the gwasSumStats's ldSketch
# (which is, by the single-block convention, the only block this call
# concerns).
# @noRd
.ctwasSingleBlockLdLoader <- function(gwasLd) {
  function(LD_file, ...) {
    snpInfo <- getSnpInfo(gwasLd)
    idx <- seq_len(nrow(snpInfo))
    block <- extractBlockGenotypes(gwasLd, idx, meanImpute = TRUE)
    geno  <- t(SummarizedExperiment::assay(block, "dosage"))
    ld    <- computeLd(geno, method = "sample")
    dimnames(ld) <- list(as.character(snpInfo$SNP),
                         as.character(snpInfo$SNP))
    ld
  }
}

# Single-block SNP-info loader for ctwas.
# @noRd
.ctwasSingleBlockSnpInfoLoader <- function(gwasLd) {
  function(LD_file, ...) .ctwasSnpInfoForBlock(gwasLd)
}
