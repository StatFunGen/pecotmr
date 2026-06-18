#' @title Causal Inference Pipeline (TWAS-Z + Mendelian Randomization)
#' @description Per-region pipeline that pairs QTL-derived weight vectors
#'   (\code{\link{TwasWeights}} and/or a QTL
#'   \code{\link{QtlFineMappingResult}}) with one or more GWAS studies
#'   (\code{\link{GwasSumStats}}) to produce per-tuple TWAS Z-scores and,
#'   when fine-mapping is supplied, Wald-ratio Mendelian Randomization
#'   estimates over the QTL credible sets.
#'
#'   Input combinations:
#'   \itemize{
#'     \item \code{twasWeights} alone (no \code{fineMappingResult}):
#'           TWAS Z only, no MR.
#'     \item \code{fineMappingResult} alone (no \code{twasWeights}):
#'           TWAS Z derived from SuSiE-style coefficients carried on the
#'           \code{topLoci} slot of each FineMappingEntry; plus MR.
#'     \item both: TWAS Z computed from \code{twasWeights};
#'           MR computed from \code{fineMappingResult}.
#'   }
#'
#' @section LD-sketch identity check:
#' If a QTL input (TwasWeights or QtlFineMappingResult) carries a
#' non-\code{NULL} \code{ldSketch}, it must match the \code{ldSketch} on
#' \code{gwasSumStats}. Mismatch is a hard error. A QTL input with
#' \code{ldSketch = NULL} (the fit was learned from individual-level
#' data) skips the validation for that input.
#'
#' @section Output shape:
#' A long-format \code{GRanges} with one row per
#' \code{(qtlStudy, context, trait, method, gwasStudy)} tuple,
#' positioned at the variant span of the QTL weight set, with mcols:
#' \describe{
#'   \item{\code{qtlStudy}, \code{context}, \code{trait},
#'         \code{method}, \code{gwasStudy}}{Identity columns.}
#'   \item{\code{twasZ}, \code{twasPval}}{Per-tuple TWAS Z and p-value.}
#'   \item{\code{waldRatio}, \code{waldRatioSe}, \code{mrPval}}{Per-tuple
#'         IVW-aggregated Wald-ratio MR estimate, standard error, and
#'         p-value. Present only when MR was computed; \code{NA}
#'         otherwise.}
#'   \item{\code{nIV}}{Number of instrumental variables used in the
#'         MR aggregation.}
#' }
#'
#' @param gwasSumStats A \code{\link{GwasSumStats}} object. Must be
#'   QC'd (\code{length(getQcInfo(x)) > 0L}).
#' @param twasWeights Optional \code{\link{TwasWeights}} carrying
#'   per-(study, context, trait, method) weights. When supplied, drives
#'   the TWAS-Z computation.
#' @param fineMappingResult Optional \code{\link{QtlFineMappingResult}}.
#'   When supplied, drives the MR computation and (when
#'   \code{twasWeights = NULL}) the TWAS-Z weights via the SuSiE-style
#'   coefficients on each entry's \code{topLoci}.
#' @param mrPipCutoff Numeric (length 1). PIP threshold for an entry's
#'   \code{topLoci} variant to be used as an instrumental variable. Default
#'   \code{0.5}.
#' @param combineMethods Optional character vector forwarded to
#'   \code{\link{combinePValues}} for cross-method combination per
#'   \code{(qtlStudy, context, trait, gwasStudy)} group. \code{NULL}
#'   (default) skips combination.
#' @param ... Reserved.
#' @return A \code{GRanges} as described above.
#' @export
causalInferencePipeline <- function(gwasSumStats,
                                    twasWeights = NULL,
                                    fineMappingResult = NULL,
                                    mrPipCutoff = 0.5,
                                    combineMethods = NULL,
                                    ...) {
  # --- Input validation --------------------------------------------------
  if (!methods::is(gwasSumStats, "GwasSumStats")) {
    stop("`gwasSumStats` must be a GwasSumStats object.")
  }
  if (length(getQcInfo(gwasSumStats)) == 0L) {
    stop("causalInferencePipeline: gwasSumStats has no QC record ",
         "(getQcInfo() is empty). Call summaryStatsQc() first.")
  }
  if (is.null(twasWeights) && is.null(fineMappingResult)) {
    stop("causalInferencePipeline: at least one of `twasWeights` or ",
         "`fineMappingResult` must be supplied.")
  }
  if (!is.null(twasWeights) && !methods::is(twasWeights, "TwasWeights")) {
    stop("`twasWeights` must be a TwasWeights object or NULL.")
  }
  if (!is.null(fineMappingResult) &&
      !methods::is(fineMappingResult, "QtlFineMappingResult")) {
    stop("`fineMappingResult` must be a QtlFineMappingResult or NULL ",
         "(causalInferencePipeline does not accept GWAS-side fine ",
         "mapping for the QTL slot).")
  }

  gwasLd <- getLdSketch(gwasSumStats)
  if (!is.null(twasWeights)) {
    twLd <- getLdSketch(twasWeights)
    .cipRequireMatchingLdSketches(twLd, gwasLd, label = "twasWeights")
  }
  if (!is.null(fineMappingResult)) {
    fmrLd <- getLdSketch(fineMappingResult)
    .cipRequireMatchingLdSketches(fmrLd, gwasLd, label = "fineMappingResult")
  }

  # --- Build the (qtlStudy, context, trait, method) work list ----
  qtlRows <- .cipBuildQtlWorkList(twasWeights, fineMappingResult)
  if (nrow(qtlRows) == 0L) {
    stop("causalInferencePipeline: no QTL tuples to score (the supplied ",
         "twasWeights / fineMappingResult collections are empty).")
  }
  qtlRows$useFmrForWeights <- is.null(twasWeights)

  # --- Per-tuple loop: compute TWAS Z + (optional) MR --------------------
  outRows <- list()
  for (qi in seq_len(nrow(qtlRows))) {
    qStudy   <- qtlRows$qtlStudy[[qi]]
    qContext <- qtlRows$context[[qi]]
    qTrait   <- qtlRows$trait[[qi]]
    qMethod  <- qtlRows$method[[qi]]

    weightsInfo <- .cipExtractWeights(
      twasWeights        = twasWeights,
      fineMappingResult  = fineMappingResult,
      study              = qStudy,
      context            = qContext,
      trait              = qTrait,
      method             = qMethod,
      useFmr             = qtlRows$useFmrForWeights[[qi]])
    if (is.null(weightsInfo)) next
    wVariantIds <- weightsInfo$variantIds
    wVec        <- weightsInfo$weights

    fmrEntry <- NULL
    if (!is.null(fineMappingResult) && .cipFmrHasTuple(
          fineMappingResult, qStudy, qContext, qTrait, qMethod)) {
      fmrEntry <- getFineMappingResult(
        fineMappingResult, study = qStudy, context = qContext,
        trait = qTrait, method = qMethod)
    }

    for (gi in seq_len(nrow(gwasSumStats))) {
      gStudy <- as.character(gwasSumStats$study)[[gi]]
      gGr    <- gwasSumStats$entry[[gi]]
      twasOut <- .cipComputeTwasZ(
        weights = wVec, variantIds = wVariantIds,
        gwasGr  = gGr, gwasLd = gwasLd)
      if (is.null(twasOut)) next

      mrOut <- if (!is.null(fmrEntry)) {
        .cipComputeMr(fmrEntry = fmrEntry, gwasGr = gGr,
                      pipCutoff = mrPipCutoff)
      } else {
        list(waldRatio = NA_real_, waldRatioSe = NA_real_,
             mrPval = NA_real_, nIV = NA_integer_)
      }

      outRows[[length(outRows) + 1L]] <- list(
        qtlStudy    = qStudy,
        context  = qContext,
        trait    = qTrait,
        method   = qMethod,
        gwasStudy   = gStudy,
        twasZ       = twasOut$Z,
        twasPval    = twasOut$pval,
        waldRatio   = mrOut$waldRatio,
        waldRatioSe = mrOut$waldRatioSe,
        mrPval      = mrOut$mrPval,
        nIV         = mrOut$nIV,
        chrom       = twasOut$chrom,
        startPos    = twasOut$startPos,
        endPos      = twasOut$endPos)
    }
  }

  if (length(outRows) == 0L) {
    stop("causalInferencePipeline: no (qtl, gwas) tuples produced a result.")
  }

  out <- .cipRowsToGranges(outRows)

  if (!is.null(combineMethods)) {
    out <- .cipCombineAcrossMethods(out, methods = combineMethods)
  }
  out
}

# =============================================================================
# Internal helpers
# =============================================================================

# Compare two GenotypeHandles for LD-sketch identity. Mirrors the
# colocboost helper but lets a NULL qtl-side handle skip the check.
.cipRequireMatchingLdSketches <- function(qtlLd, gwasLd, label) {
  if (is.null(qtlLd)) return(invisible(NULL))
  if (!methods::is(qtlLd, "GenotypeHandle") ||
      !methods::is(gwasLd, "GenotypeHandle")) {
    stop("causalInferencePipeline: ldSketch on `", label,
         "` and gwasSumStats must both be GenotypeHandle objects for ",
         "the cross-pipeline LD reference check.")
  }
  qSnp <- getSnpInfo(qtlLd)
  gSnp <- getSnpInfo(gwasLd)
  if (nrow(qSnp) != nrow(gSnp))
    stop("causalInferencePipeline: ldSketch panels on `", label,
         "` (", nrow(qSnp), " variants) and gwasSumStats (",
         nrow(gSnp), " variants) differ in size; the two ldSketch ",
         "GenotypeHandles must match exactly.")
  for (col in c("SNP", "CHR", "BP", "A1", "A2")) {
    if (!identical(as.character(qSnp[[col]]),
                   as.character(gSnp[[col]])))
      stop("causalInferencePipeline: ldSketch panels on `", label,
           "` and gwasSumStats differ in column ", col,
           "; use the same ldSketch on both.")
  }
  if (!identical(getSampleIds(qtlLd), getSampleIds(gwasLd)))
    stop("causalInferencePipeline: ldSketch panels on `", label,
         "` and gwasSumStats have different sample sets; use the same ",
         "ldSketch on both.")
  invisible(NULL)
}

# Build the (qtlStudy, context, trait, method) work list from
# whichever input was supplied. When both are supplied, prefer the
# TwasWeights tuples and only retain those that also appear in the FMR
# (so MR has something to attach).
.cipBuildQtlWorkList <- function(twasWeights, fineMappingResult) {
  if (!is.null(twasWeights)) {
    df <- data.frame(
      qtlStudy   = as.character(twasWeights$study),
      context = as.character(twasWeights$context),
      trait   = as.character(twasWeights$trait),
      method  = as.character(twasWeights$method),
      stringsAsFactors = FALSE)
  } else {
    df <- data.frame(
      qtlStudy   = as.character(fineMappingResult$study),
      context = as.character(fineMappingResult$context),
      trait   = as.character(fineMappingResult$trait),
      method  = as.character(fineMappingResult$method),
      stringsAsFactors = FALSE)
  }
  df
}

.cipFmrHasTuple <- function(fmr, study, context, trait, method) {
  any(as.character(fmr$study)   == study &
      as.character(fmr$context) == context &
      as.character(fmr$trait)   == trait &
      as.character(fmr$method)  == method)
}

# Extract the per-tuple weights vector. From TwasWeights: read the
# TwasWeightsEntry. From FineMappingResult: extract the SuSiE-style
# coefficient (betahat) from topLoci.
.cipExtractWeights <- function(twasWeights, fineMappingResult,
                               study, context, trait, method, useFmr) {
  if (!useFmr) {
    if (!any(as.character(twasWeights$study)   == study &
             as.character(twasWeights$context) == context &
             as.character(twasWeights$trait)   == trait &
             as.character(twasWeights$method)  == method)) return(NULL)
    twEntry <- getTwasWeights(twasWeights, study = study, context = context,
                              trait = trait, method = method)
    vids <- getVariantIds(twEntry)
    w    <- as.numeric(getWeights(twEntry))
    if (length(vids) != length(w) || length(vids) == 0L) return(NULL)
    return(list(variantIds = vids, weights = w))
  }
  # FMR-based weights: pull from the entry's topLoci$betahat column.
  if (!.cipFmrHasTuple(fineMappingResult, study, context, trait, method))
    return(NULL)
  ent <- getFineMappingResult(fineMappingResult, study = study,
                              context = context, trait = trait,
                              method = method)
  tl <- getTopLoci(ent)
  if (is.null(tl) || nrow(tl) == 0L) return(NULL)
  betaCol <- intersect(c("betahat", "beta", "bhat_x"), colnames(tl))
  if (length(betaCol) == 0L) return(NULL)
  vids <- as.character(tl$variant_id)
  w    <- as.numeric(tl[[betaCol[[1L]]]])
  ok   <- !is.na(vids) & !is.na(w)
  if (sum(ok) == 0L) return(NULL)
  list(variantIds = vids[ok], weights = w[ok])
}

# Compute the per-tuple TWAS Z for one (weights vector, GWAS GRanges).
# Returns NULL when the overlap is too small.
.cipComputeTwasZ <- function(weights, variantIds, gwasGr, gwasLd) {
  gwasIds <- as.character(S4Vectors::mcols(gwasGr)$SNP)
  gwasZ   <- as.numeric(S4Vectors::mcols(gwasGr)$Z)
  common <- intersect(variantIds, gwasIds)
  if (length(common) < 2L) return(NULL)

  wSub <- weights[match(common, variantIds)]
  zSub <- gwasZ[match(common, gwasIds)]
  ldMat <- .cipLdFromSketch(gwasLd, common)

  res <- twasZ(weights = wSub, z = zSub, R = ldMat)
  zMat <- res$Z
  zVal <- as.numeric(zMat[1L, "Z"])
  pVal <- as.numeric(zMat[1L, "pval"])
  # Position the row at the variant span.
  idx <- match(common, gwasIds)
  chrom    <- as.character(GenomicRanges::seqnames(gwasGr))[[idx[[1L]]]]
  startPos <- min(GenomicRanges::start(gwasGr)[idx])
  endPos   <- max(GenomicRanges::end(gwasGr)[idx])
  list(Z = zVal, pval = pVal,
       chrom = chrom, startPos = startPos, endPos = endPos)
}

# Build an LD correlation matrix for a given variant subset from an
# LD sketch. Same pattern as .twasLdFromSketch / .fmLdFromSketch.
.cipLdFromSketch <- function(ldSketch, variantIds) {
  snpInfo <- getSnpInfo(ldSketch)
  idx <- match(variantIds, as.character(snpInfo$SNP))
  if (anyNA(idx)) {
    stop("causalInferencePipeline: ", sum(is.na(idx)),
         " variant id(s) are absent from the ldSketch panel; the ",
         "summaryStatsQc step should have removed them.")
  }
  block <- extractBlockGenotypes(ldSketch, idx, meanImpute = TRUE)
  geno  <- t(SummarizedExperiment::assay(block, "dosage"))
  colnames(geno) <- variantIds
  ld <- computeLd(geno, method = "sample")
  dimnames(ld) <- list(variantIds, variantIds)
  ld
}

# Compute the Wald-ratio IVW MR estimate for a single tuple. Uses the
# FineMappingEntry's topLoci as the instrumental variable source: each
# variant with PIP > pipCutoff contributes one ratio = beta_y / beta_x.
# Returns list(waldRatio, waldRatioSe, mrPval, nIV) with NA fields when
# no IVs survive.
.cipComputeMr <- function(fmrEntry, gwasGr, pipCutoff) {
  tl <- getTopLoci(fmrEntry)
  if (is.null(tl) || nrow(tl) == 0L)
    return(list(waldRatio = NA_real_, waldRatioSe = NA_real_,
                mrPval = NA_real_, nIV = 0L))
  pipCol  <- intersect(c("pip", "PIP"), colnames(tl))
  betaCol <- intersect(c("betahat", "beta", "bhat_x"), colnames(tl))
  seCol   <- intersect(c("sebetahat", "se", "sbhat_x"), colnames(tl))
  if (length(pipCol) == 0L || length(betaCol) == 0L || length(seCol) == 0L)
    return(list(waldRatio = NA_real_, waldRatioSe = NA_real_,
                mrPval = NA_real_, nIV = 0L))
  keep <- !is.na(tl[[pipCol[[1L]]]]) & tl[[pipCol[[1L]]]] > pipCutoff
  ivVars <- as.character(tl$variant_id)[keep]
  if (length(ivVars) == 0L)
    return(list(waldRatio = NA_real_, waldRatioSe = NA_real_,
                mrPval = NA_real_, nIV = 0L))
  betaX <- as.numeric(tl[[betaCol[[1L]]]])[keep]
  seX   <- as.numeric(tl[[seCol[[1L]]]])[keep]
  gIds  <- as.character(S4Vectors::mcols(gwasGr)$SNP)
  gIdx  <- match(ivVars, gIds)
  ok    <- !is.na(gIdx)
  if (sum(ok) == 0L)
    return(list(waldRatio = NA_real_, waldRatioSe = NA_real_,
                mrPval = NA_real_, nIV = 0L))
  betaX <- betaX[ok]; seX <- seX[ok]; gIdx <- gIdx[ok]
  gZ <- as.numeric(S4Vectors::mcols(gwasGr)$Z)[gIdx]
  gN <- if ("N" %in% colnames(S4Vectors::mcols(gwasGr)))
    as.numeric(S4Vectors::mcols(gwasGr)$N)[gIdx]
        else rep(NA_real_, length(gIdx))
  gMaf <- if ("MAF" %in% colnames(S4Vectors::mcols(gwasGr)))
    as.numeric(S4Vectors::mcols(gwasGr)$MAF)[gIdx]
          else rep(NA_real_, length(gIdx))
  betaY <- .cipZToBeta(gZ, gMaf, gN)
  seY   <- .cipZToSe(gZ, gMaf, gN)
  ratio <- betaY / betaX
  # Standard Wald-ratio SE via delta method.
  rSe   <- sqrt((seY / betaX)^2 + (betaY * seX / betaX^2)^2)
  # IVW pooling: weight by 1/rSe^2; pooled SE = 1/sqrt(sum(w)).
  w <- 1 / rSe^2
  validW <- is.finite(w) & w > 0
  if (sum(validW) == 0L)
    return(list(waldRatio = NA_real_, waldRatioSe = NA_real_,
                mrPval = NA_real_, nIV = 0L))
  ratio <- ratio[validW]; rSe <- rSe[validW]; w <- w[validW]
  meta  <- sum(w * ratio) / sum(w)
  metaSe <- 1 / sqrt(sum(w))
  pval  <- 2 * stats::pnorm(-abs(meta / metaSe))
  list(waldRatio = meta, waldRatioSe = metaSe,
       mrPval = pval, nIV = length(ratio))
}

# Derive beta from z using maf + n. beta = z * se. With se = 1/sqrt(2*n*p*q),
# beta = z / sqrt(2*n*p*q).
.cipZToBeta <- function(z, maf, n) {
  if (any(is.na(maf)) || any(is.na(n)))
    return(z)  # fall back to z as a beta surrogate when no maf/n
  z / sqrt(2 * n * maf * (1 - maf))
}
.cipZToSe <- function(z, maf, n) {
  if (any(is.na(maf)) || any(is.na(n)))
    return(rep(1, length(z)))
  1 / sqrt(2 * n * maf * (1 - maf))
}

# Convert the accumulated list of row records to a GRanges with mcols.
.cipRowsToGranges <- function(rows) {
  df <- do.call(rbind.data.frame, lapply(rows, as.data.frame,
                                        stringsAsFactors = FALSE))
  chr <- paste0("chr", sub("^chr", "", as.character(df$chrom),
                            ignore.case = TRUE))
  gr <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges   = IRanges::IRanges(start = as.integer(df$startPos),
                                end   = as.integer(df$endPos)))
  mcols <- df[, c("qtlStudy", "context", "trait", "method",
                  "gwasStudy", "twasZ", "twasPval",
                  "waldRatio", "waldRatioSe", "mrPval", "nIV")]
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(mcols)
  gr
}

# Combine TWAS p-values across method for each (qtlStudy, context,
# trait, gwasStudy) group. Appends one row per group with a
# combined p-value and methodName = "combined.<methodToken>". Uses
# combinePValues() with the cross-method correlation set to the identity
# (we have no cross-method covariance available downstream).
.cipCombineAcrossMethods <- function(gr, methods) {
  mc <- as.data.frame(S4Vectors::mcols(gr))
  key <- paste(mc$qtlStudy, mc$context, mc$trait,
               mc$gwasStudy, sep = "||")
  groups <- split(seq_len(nrow(mc)), key)
  extras <- list()
  for (gkey in names(groups)) {
    rows <- groups[[gkey]]
    if (length(rows) < 2L) next
    pvals <- as.numeric(mc$twasPval[rows])
    zvec  <- as.numeric(mc$twasZ[rows])
    cp <- combinePValues(pvals = pvals, zScores = zvec,
                         methods = methods)
    for (m in methods) {
      newRow <- mc[rows[[1L]], , drop = FALSE]
      newRow$method <- paste0("combined.", m)
      newRow$twasZ    <- NA_real_
      newRow$twasPval <- as.numeric(cp$results[[m]]$pval)
      newRow$waldRatio   <- NA_real_
      newRow$waldRatioSe <- NA_real_
      newRow$mrPval      <- NA_real_
      newRow$nIV         <- NA_integer_
      extras[[length(extras) + 1L]] <- newRow
    }
  }
  if (length(extras) == 0L) return(gr)
  newMcs <- do.call(rbind, extras)
  newGr <- gr[rep(1L, nrow(newMcs))]
  S4Vectors::mcols(newGr) <- S4Vectors::DataFrame(newMcs)
  c(gr, newGr)
}
