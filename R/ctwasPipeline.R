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
#' @param fineMappingResult Optional \code{QtlFineMappingResult} or
#'   \code{GwasFineMappingResult} carrying the per-variant PIP and
#'   credible-set membership data used by the CS / PIP rescue filters
#'   (\code{csMinCor} and \code{minPipCutoff}). When \code{NULL}
#'   (default) the smart filters are no-ops; only the magnitude filter
#'   (\code{twasWeightCutoff}) and the per-gene cap
#'   (\code{maxNumVariants}, ordered by \code{|weight|}) apply.
#' @param regionId Optional character (length 1) label for the LD
#'   block. Default \code{"block1"}.
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
#' @param ... Additional arguments forwarded to
#'   \code{ctwas::ctwas_sumstats}.
#' @return Whatever \code{ctwas::ctwas_sumstats} returns (a list with
#'   \code{susie_alpha_res}, \code{param}, and other diagnostics).
#' @export
ctwasPipeline <- function(gwasSumStats,
                          twasWeights,
                          twasZ                   = NULL,
                          fineMappingResult       = NULL,
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
                          twasWeightCutoff        = 0,
                          csMinCor                = 0.8,
                          minPipCutoff            = 0,
                          maxNumVariants          = Inf,
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
  if (!is.null(fineMappingResult) &&
      !methods::is(fineMappingResult, "FineMappingResultBase"))
    stop("`fineMappingResult` must be a FineMappingResultBase ",
         "(QtlFineMappingResult or GwasFineMappingResult) or NULL.")
  if (length(regionId) != 1L || !nzchar(regionId))
    stop("`regionId` must be a single non-empty character string.")
  groupPriorVarStructure <- match.arg(groupPriorVarStructure)

  twLd   <- getLdSketch(twasWeights)
  gwasLd <- getLdSketch(gwasSumStats)
  .ctwasRequireMatchingLdSketches(twLd, gwasLd)

  # --- Compute the full-panel LD ONCE -------------------------------
  # Single source of truth for both the LD-loader closure (which ctwas
  # invokes per region during assemble + fine-map stages) and the
  # per-gene R_wgt submatrices (sliced from this cache by SNP ID).
  ldPanel <- .ctwasComputeFullPanelLd(gwasLd)

  # --- Build the single-region ctwas inputs ---------------------------
  zSnp        <- .ctwasBuildZSnp(gwasSumStats)
  regionInfo  <- .ctwasBuildSingleRegionInfo(regionId, gwasLd)
  # ctwas::ctwas_sumstats top-level asserts `file.exists(LD_map$LD_file)`
  # and `file.exists(LD_map$SNP_file)` unconditionally — even when
  # LD_format = "custom" routes all data through our loaders and the file
  # paths are never read. The right fix is upstream (gate the assertion
  # on `LD_format != "custom"` or drop it entirely; see
  # https://github.com/xinhe-lab/ctwas — `ctwas_sumstats()` L33-34). Until
  # that lands, point both columns at `tempdir()`: always exists, no disk
  # writes, no cleanup. The loader closures ignore the file token.
  vestigialPath <- tempdir()
  ldMap <- data.frame(region_id = regionId,
                      LD_file   = vestigialPath,
                      SNP_file  = vestigialPath,
                      stringsAsFactors = FALSE)
  snpMap      <- list()
  snpMap[[regionId]] <- ldPanel$snpInfo
  weightsList <- .ctwasBuildWeights(
    twasWeights, ldPanel,
    fineMappingResult = fineMappingResult,
    twasWeightCutoff  = twasWeightCutoff,
    csMinCor          = csMinCor,
    minPipCutoff      = minPipCutoff,
    maxNumVariants    = maxNumVariants)
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
    LD_loader_fun              = .ctwasSingleBlockLdLoader(ldPanel$R),
    snpinfo_loader_fun         = .ctwasSingleBlockSnpInfoLoader(ldPanel$snpInfo),
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

# Per-block SNP info table (chrom, id, pos, alt, ref). ctwas requires
# these exact column names (read_snp_info_files asserts them). `alt`
# maps to A1 (effect allele) and `ref` to A2.
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
#   R       : full-panel correlation matrix (n_var x n_var, dimnames =
#             SNP IDs). Single source of truth for both the per-region
#             LD loader closure and the per-gene R_wgt submatrices.
#   snpInfo : ctwas-shaped per-block table (chrom, id, pos, alt, ref)
#             — both the snp_map element and the snpinfo loader return.
# @noRd
.ctwasComputeFullPanelLd <- function(gwasLd) {
  snpInfoCtwas <- .ctwasSnpInfoForBlock(gwasLd)
  block <- extractBlockGenotypes(gwasLd, seq_len(nrow(snpInfoCtwas)),
                                  meanImpute = TRUE)
  geno  <- t(SummarizedExperiment::assay(block, "dosage"))
  R <- computeLd(geno, method = "sample")
  snpIds <- snpInfoCtwas$id
  dimnames(R) <- list(snpIds, snpIds)
  list(R = R, snpInfo = snpInfoCtwas)
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
                               maxNumVariants    = Inf) {
  panelSnps <- rownames(ldPanel$R)
  panelInfo <- ldPanel$snpInfo
  out <- list()
  for (i in seq_len(nrow(twasWeights))) {
    entry  <- twasWeights$entry[[i]]
    vids   <- getVariantIds(entry)
    w      <- as.numeric(getWeights(entry))
    if (length(vids) == 0L || length(vids) != length(w)) next

    # Drop variants not in the LD sketch panel.
    keep <- vids %in% panelSnps
    if (!any(keep)) next
    vids <- vids[keep]; w <- w[keep]

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

# Single-block LD loader for ctwas: captures the precomputed full-panel
# correlation matrix and returns it on every loader call (ctwas invokes
# the loader multiple times per region across assemble + fine-map). The
# `LD_file` argument is a vestigial region token — ignored.
# @noRd
.ctwasSingleBlockLdLoader <- function(R) {
  function(LD_file, ...) R
}

# Single-block SNP-info loader for ctwas: captures the precomputed
# per-block snpInfo table and returns it on every loader call.
# @noRd
.ctwasSingleBlockSnpInfoLoader <- function(snpInfo) {
  function(LD_file, ...) snpInfo
}
