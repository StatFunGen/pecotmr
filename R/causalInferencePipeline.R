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
#'   \code{topLoci} variant to be used as an instrumental variable.
#'   Used only when \code{mrMethod = "ivwPerVariant"}. Default \code{0.5}.
#' @param mrMethod One of \code{"ivwPerVariant"} (default) or
#'   \code{"csAware"}. The IVW-per-variant method filters topLoci
#'   variants by \code{pip > mrPipCutoff} and IVW-pools Wald ratios
#'   across variants. The CS-aware method groups variants by credible
#'   set (column \code{cs} in topLoci), computes a PIP-weighted
#'   composite Wald ratio per CS using \code{mrCpipCutoff} on the
#'   per-CS cumulative PIP, then IVW-pools across CSs and reports
#'   Cochran's Q + I-squared in the output columns \code{Q}, \code{I2}.
#' @param mrCpipCutoff Numeric (length 1). Cumulative-PIP cutoff for
#'   retaining a credible set. Used only when
#'   \code{mrMethod = "csAware"}. Default \code{0.5}.
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
                                    mrMethod = c("ivwPerVariant", "csAware"),
                                    mrCpipCutoff = 0.5,
                                    combineMethods = NULL,
                                    ...) {
  mrMethod <- match.arg(mrMethod)
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
      gdf    <- getSumstatDf(gwasSumStats, study = gStudy,
                             require = c("SNP", "Z"))
      twasOut <- .cipComputeTwasZ(
        weights = wVec, variantIds = wVariantIds,
        gwasDf  = gdf, gwasLd = gwasLd)
      if (is.null(twasOut)) next

      mrOut <- if (!is.null(fmrEntry)) {
        if (mrMethod == "csAware") {
          .cipComputeMrCsAware(fmrEntry = fmrEntry, gwasDf = gdf,
                                cpipCutoff = mrCpipCutoff)
        } else {
          .cipComputeMr(fmrEntry = fmrEntry, gwasDf = gdf,
                        pipCutoff = mrPipCutoff)
        }
      } else {
        list(waldRatio = NA_real_, waldRatioSe = NA_real_,
             mrPval = NA_real_, nIV = NA_integer_,
             Q = NA_real_, I2 = NA_real_, nCs = NA_integer_)
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
        Q           = mrOut$Q   %||% NA_real_,
        I2          = mrOut$I2  %||% NA_real_,
        nCs         = mrOut$nCs %||% NA_integer_,
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

# Compare two GenotypeHandles for LD-sketch identity. Thin wrapper over
# the shared `.requireMatchingLdSketches` helper (R/ld.R).
.cipRequireMatchingLdSketches <- function(qtlLd, gwasLd, label) {
  .requireMatchingLdSketches(qtlLd, gwasLd,
                             pipelineName = "causalInferencePipeline",
                             label = label)
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
  length(.matchTupleRows(fmr,
                          list(study = study, context = context,
                               trait = trait, method = method))) > 0L
}

# Extract the per-tuple weights vector. From TwasWeights: read the
# TwasWeightsEntry. From FineMappingResult: extract the SuSiE-style
# coefficient (betahat) from topLoci.
.cipExtractWeights <- function(twasWeights, fineMappingResult,
                               study, context, trait, method, useFmr) {
  if (!useFmr) {
    if (length(.matchTupleRows(twasWeights,
                               list(study = study, context = context,
                                    trait = trait, method = method))) == 0L)
      return(NULL)
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

# Compute the per-tuple TWAS Z from a single GwasSumStats tuple's
# unpacked data.frame (produced by getSumstatDf upstream). Returns
# NULL when the overlap is too small.
.cipComputeTwasZ <- function(weights, variantIds, gwasDf, gwasLd) {
  common <- intersect(variantIds, gwasDf$variant_id)
  if (length(common) < 2L) return(NULL)

  wSub <- weights[match(common, variantIds)]
  zSub <- gwasDf$z[match(common, gwasDf$variant_id)]
  ldMat <- .cipLdFromSketch(gwasLd, common)

  res <- twasZ(weights = wSub, z = zSub, R = ldMat)
  zMat <- res$Z
  zVal <- as.numeric(zMat[1L, "Z"])
  pVal <- as.numeric(zMat[1L, "pval"])
  # Position the row at the variant span.
  idx <- match(common, gwasDf$variant_id)
  chrom    <- gwasDf$chrom[[idx[[1L]]]]
  startPos <- min(gwasDf$pos[idx])
  endPos   <- max(gwasDf$pos[idx])
  list(Z = zVal, pval = pVal,
       chrom = chrom, startPos = startPos, endPos = endPos)
}

# Build an LD correlation matrix for a given variant subset from an
# LD sketch. Thin wrapper over the shared `.ldFromSketch` helper
# (R/ld.R).
.cipLdFromSketch <- function(ldSketch, variantIds) {
  .ldFromSketch(ldSketch, variantIds, label = "causalInferencePipeline")
}

# Compute the Wald-ratio IVW MR estimate for a single tuple. Uses the
# FineMappingEntry's topLoci as the instrumental variable source: each
# variant with PIP > pipCutoff contributes one ratio = beta_y / beta_x.
# Returns list(waldRatio, waldRatioSe, mrPval, nIV) with NA fields when
# no IVs survive.
.cipComputeMr <- function(fmrEntry, gwasDf, pipCutoff) {
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
  gIdx  <- match(ivVars, gwasDf$variant_id)
  ok    <- !is.na(gIdx)
  if (sum(ok) == 0L)
    return(list(waldRatio = NA_real_, waldRatioSe = NA_real_,
                mrPval = NA_real_, nIV = 0L))
  betaX <- betaX[ok]; seX <- seX[ok]; gIdx <- gIdx[ok]
  gZ <- gwasDf$z[gIdx]
  gN <- if (!is.null(gwasDf$N)) gwasDf$N[gIdx]
        else rep(NA_real_, length(gIdx))
  gMaf <- if (!is.null(gwasDf$maf)) gwasDf$maf[gIdx]
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

# Cochran's Q-based I-squared heterogeneity statistic. Ported from
# the legacy mr.R::calcI2. Q = 0 or near-zero -> I2 = 0; clipped to
# [0, 1] to stay in the usual heterogeneity convention.
# @noRd
.cipCalcI2 <- function(Q, nGroups) {
  if (!is.finite(Q) || Q <= 1e-3 || nGroups <= 1L) return(0)
  i2 <- (Q - (nGroups - 1)) / Q
  max(0, min(1, i2))
}

# CS-aware Wald-ratio MR estimate for a single tuple. Adapted from the
# legacy mr.R::mrAnalysis but operating on a `FineMappingEntry`'s topLoci
# data.frame + a GWAS `GRanges` instead of the old data.frame-only API.
#
# Logic:
#   1. Group topLoci variants by credible set (column `cs` in topLoci, or
#      first column matching ^cs).
#   2. Per CS, compute cumulative PIP; drop CSs with cpip < cpipCutoff.
#   3. Per surviving CS, compute a PIP-weighted composite Wald ratio
#      `composite_bhat = sum((bhatY/bhatX * pip) / cpip)` with the
#      delta-method composite SE.
#   4. IVW-pool composite Walds across CSs with weights wv = 1/se^2.
#   5. Report (metaEff, metaSe, metaPval, Q, I2, nCs, nIv) where the
#      result list keys are renamed (waldRatio/waldRatioSe/mrPval/nIV) so
#      the calling pipeline emits consistent column names.
#
# Returns list(waldRatio, waldRatioSe, mrPval, nIV, Q, I2, nCs) with NA
# fields when no usable CS survives.
# @noRd
.cipComputeMrCsAware <- function(fmrEntry, gwasDf, cpipCutoff) {
  naResult <- list(waldRatio = NA_real_, waldRatioSe = NA_real_,
                   mrPval = NA_real_, nIV = 0L,
                   Q = NA_real_, I2 = NA_real_, nCs = 0L)

  tl <- getTopLoci(fmrEntry)
  if (is.null(tl) || nrow(tl) == 0L) return(naResult)

  pipCol  <- intersect(c("pip", "PIP"), colnames(tl))
  betaCol <- intersect(c("betahat", "beta", "bhat_x"), colnames(tl))
  seCol   <- intersect(c("sebetahat", "se", "sbhat_x"), colnames(tl))
  csCol   <- intersect(c("cs"), colnames(tl))
  if (length(csCol) == 0L) {
    # Look for any column starting with "cs" (e.g. cs_0.95, cs_susie).
    csCandidates <- grep("^cs", colnames(tl), value = TRUE)
    if (length(csCandidates) > 0L) csCol <- csCandidates[[1L]]
  }
  if (length(pipCol) == 0L || length(betaCol) == 0L ||
      length(seCol) == 0L || length(csCol) == 0L) return(naResult)

  cs    <- as.integer(tl[[csCol[[1L]]]])
  pip   <- as.numeric(tl[[pipCol[[1L]]]])
  bhatX <- as.numeric(tl[[betaCol[[1L]]]])
  sbhatX <- as.numeric(tl[[seCol[[1L]]]])
  vids  <- as.character(tl$variant_id)

  # Drop rows with NA cs / non-positive cs / NA pip / NA beta.
  ok <- !is.na(cs) & cs > 0L & !is.na(pip) & !is.na(bhatX) & !is.na(sbhatX)
  if (!any(ok)) return(naResult)
  cs   <- cs[ok];   pip  <- pip[ok]
  bhatX <- bhatX[ok]; sbhatX <- sbhatX[ok]; vids <- vids[ok]

  # Match GWAS-side beta/se via the existing helpers.
  gIdx <- match(vids, gwasDf$variant_id)
  ok   <- !is.na(gIdx)
  if (!any(ok)) return(naResult)
  cs <- cs[ok]; pip <- pip[ok]
  bhatX <- bhatX[ok]; sbhatX <- sbhatX[ok]
  gIdx <- gIdx[ok]; vids <- vids[ok]

  gZ <- gwasDf$z[gIdx]
  gN <- if (!is.null(gwasDf$N)) gwasDf$N[gIdx]
        else rep(NA_real_, length(gIdx))
  gMaf <- if (!is.null(gwasDf$maf)) gwasDf$maf[gIdx]
          else rep(NA_real_, length(gIdx))
  bhatY  <- .cipZToBeta(gZ, gMaf, gN)
  sbhatY <- .cipZToSe(gZ, gMaf, gN)

  # Standardize bhatX -> z; sbhatX -> 1. (Matches the legacy mrAnalysis
  # rescaling step that puts the exposure on a unit-SE footing.)
  bhatX  <- bhatX / sbhatX
  sbhatX <- rep(1, length(bhatX))

  # Per-CS composite Wald.
  csIds <- sort(unique(cs))
  composite_bhat  <- numeric(0)
  composite_sbhat <- numeric(0)
  for (cid in csIds) {
    inCs <- which(cs == cid)
    cpip <- sum(pip[inCs])
    if (!is.finite(cpip) || cpip < cpipCutoff) next
    pNorm <- pip[inCs] / cpip
    bx <- bhatX[inCs]; sx <- sbhatX[inCs]
    by <- bhatY[inCs]; sy <- sbhatY[inCs]
    # Composite point estimate.
    cBhat <- sum((by / bx) * pNorm)
    # Composite SE: sqrt(E[(by/bx)^2 + sy^2/bx^2 + by^2*sx^2/bx^4] - cBhat^2)
    second <- sum(((by / bx)^2 + (sy^2 / bx^2) +
                   ((by^2 * sx^2) / bx^4)) * pNorm)
    cSe <- sqrt(max(0, second - cBhat^2))
    if (!is.finite(cBhat) || !is.finite(cSe) || cSe <= 0) next
    composite_bhat  <- c(composite_bhat,  cBhat)
    composite_sbhat <- c(composite_sbhat, cSe)
  }
  if (length(composite_bhat) == 0L) return(naResult)

  # IVW meta across CSs.
  wv <- 1 / composite_sbhat^2
  metaEff <- sum(wv * composite_bhat) / sum(wv)
  metaSe  <- sqrt(1 / sum(wv))
  Q       <- sum(wv * (composite_bhat - metaEff)^2)
  I2      <- .cipCalcI2(Q, length(composite_bhat))
  metaPval <- 2 * stats::pnorm(-abs(metaEff / metaSe))

  list(waldRatio = metaEff,
       waldRatioSe = metaSe,
       mrPval = metaPval,
       nIV = sum(cs %in% csIds),
       Q = Q, I2 = I2,
       nCs = length(composite_bhat))
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
                  "waldRatio", "waldRatioSe", "mrPval", "nIV",
                  "Q", "I2", "nCs")]
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


# =============================================================================
# Unified TWAS Z-statistic
# =============================================================================

# Internal: build the K x K covariance Wᵀ R W. Uses the SVD path when the
# triplet (V, D, nSketch) is supplied, otherwise the R / X path. Aligns LD
# rows/cols to the rownames of W when both are named; falls back to
# positional alignment otherwise.
.twasZCovY <- function(weights, R = NULL, X = NULL,
                       V = NULL, D = NULL, nSketch = NULL) {
  rn <- rownames(weights)
  useSvd <- !is.null(V) && !is.null(D) && !is.null(nSketch)
  if (useSvd) {
    if (!is.null(rownames(V)) && !is.null(rn)) {
      idx <- match(rn, rownames(V))
      if (anyNA(idx))
        stop("twasZ: V is missing rows for ", sum(is.na(idx)),
             " variant(s) named in weights.")
      vSub <- V[idx, , drop = FALSE]
    } else {
      if (nrow(V) != nrow(weights))
        stop("twasZ: positional alignment requires nrow(V) == nrow(weights).")
      vSub <- V
    }
    Lambda <- D^2 / (nSketch - 1)
    VtW    <- crossprod(vSub, weights)              # r x K
    covY   <- crossprod(VtW * sqrt(Lambda))          # K x K
    return(list(covY = covY))
  }
  if (is.null(R)) {
    if (is.null(X))
      stop("twasZ: provide R, X, or the (V, D, nSketch) SVD triplet.")
    R <- computeLd(X)
  }
  if (!is.null(rownames(R)) && !is.null(rn)) {
    idx <- match(rn, rownames(R))
    if (anyNA(idx))
      stop("twasZ: R is missing rows for ", sum(is.na(idx)),
           " variant(s) named in weights.")
    rSub <- R[idx, idx, drop = FALSE]
  } else {
    if (nrow(R) != nrow(weights))
      stop("twasZ: positional alignment requires nrow(R) == nrow(weights).")
    rSub <- R
  }
  covY <- crossprod(weights, rSub) %*% weights       # K x K
  list(covY = covY)
}

#' Calculate TWAS Z-Statistics for One or More Methods / Contexts
#'
#' Unified TWAS Z-statistic: accepts a weight vector (single
#' method/context) or a (variants x K) weight matrix, computes the
#' per-tuple TWAS Z-score and two-sided p-value, and optionally
#' delegates cross-tuple p-value combination to
#' \code{\link{combinePValues}}.
#'
#' For each column k of \code{weights}:
#' \itemize{
#'   \item \code{stat_k = w_kᵀ z}
#'   \item \code{denom_k = w_kᵀ R w_k}
#'   \item \code{Z_k = stat_k / sqrt(denom_k)}, \code{p_k = 2 * (1 - Phi(|Z_k|))}
#' }
#' When \code{combineMethods} is non-NULL and K >= 2, the cross-tuple
#' correlation matrix \code{rho_{i,j} = covY_{i,j} / sqrt(covY_{i,i} *
#' covY_{j,j})} is constructed once and forwarded to
#' \code{combinePValues} as the \code{R} argument. When K == 1, the
#' combined p-value trivially equals the per-tuple p-value.
#'
#' The SVD path (\code{V}, \code{D}, \code{nSketch}) lets the caller
#' avoid materializing the full LD matrix: \code{covY = (VᵀW · sqrt(Lambda))ᵀ
#' (VᵀW · sqrt(Lambda))} with \code{Lambda_i = D_i^2 / (nSketch - 1)}.
#' Use the \code{R} path when an LD correlation matrix is already
#' available; use the \code{X} path to compute \code{R} from a genotype
#' matrix.
#'
#' @param weights Numeric vector of weights (single tuple) or a numeric
#'   matrix with one column per tuple (method / context). When a
#'   vector, the column name defaults to \code{"method1"}.
#' @param z Numeric vector of GWAS Z-scores aligned to the rows of
#'   \code{weights}.
#' @param R Optional LD correlation matrix.
#' @param X Optional genotype matrix used to compute \code{R} when
#'   \code{R} is missing.
#' @param V,D,nSketch SVD components of the LD sketch (right-singular
#'   vectors, singular values, panel sample size). Supplying all three
#'   selects the SVD path.
#' @param combineMethods Optional character vector of method names to
#'   forward to \code{\link{combinePValues}} for cross-tuple
#'   combination. \code{NULL} (default) skips combination.
#' @return A list with:
#' \describe{
#'   \item{Z}{A \code{K x 2} numeric matrix with columns
#'     \code{c("Z", "pval")}; rownames are the column names of
#'     \code{weights}.}
#'   \item{combined}{Output of \code{combinePValues} (or its trivial
#'     K=1 equivalent) when \code{combineMethods} is non-NULL,
#'     otherwise \code{NULL}.}
#' }
#' @seealso \code{\link{combinePValues}} for the combination method menu.
#' @importFrom stats pnorm
#' @export
twasZ <- function(weights, z, R = NULL, X = NULL,
                  V = NULL, D = NULL, nSketch = NULL,
                  combineMethods = NULL) {
  # Coerce a numeric vector to a one-column matrix.
  if (is.numeric(weights) && is.null(dim(weights))) {
    nm <- if (!is.null(names(weights))) names(weights) else NULL
    weights <- matrix(weights, ncol = 1L,
                      dimnames = list(nm, "method1"))
  }
  if (!is.matrix(weights))
    stop("`weights` must be a numeric vector or a matrix.")
  if (is.null(colnames(weights))) {
    colnames(weights) <- paste0("method", seq_len(ncol(weights)))
  }
  if (nrow(weights) != length(z))
    stop("nrow(weights) must equal length(z).")
  K <- ncol(weights)

  covInfo <- .twasZCovY(weights = weights, R = R, X = X,
                        V = V, D = D, nSketch = nSketch)
  covY <- covInfo$covY
  ySd  <- sqrt(diag(covY))
  stats <- as.numeric(crossprod(weights, as.numeric(z)))
  zVec <- stats / ySd
  pVec <- 2 * pnorm(-abs(zVec))

  zMatrix <- cbind(Z = zVec, pval = pVec)
  rownames(zMatrix) <- colnames(weights)

  combined <- NULL
  if (!is.null(combineMethods)) {
    combineMethods <- as.character(combineMethods)
    if (K == 1L) {
      perMethod <- lapply(combineMethods, function(m) {
        list(method = m, pval = as.numeric(pVec[[1L]]))
      })
      names(perMethod) <- combineMethods
      combined <- list(
        input = list(nPvalsIn = 1L, nZScoresIn = 1L, nValid = 1L,
                     Raligned = matrix(1.0, 1L, 1L,
                                       dimnames = list(rownames(zMatrix),
                                                       rownames(zMatrix)))),
        results = perMethod)
    } else {
      sig <- covY / tcrossprod(ySd, ySd)
      rownames(sig) <- colnames(sig) <- rownames(zMatrix)
      names(pVec) <- rownames(zMatrix)
      names(zVec) <- rownames(zMatrix)
      combined <- combinePValues(
        pvals    = pVec,
        zScores  = zVec,
        methods  = combineMethods,
        R        = sig)
    }
  }

  list(Z = zMatrix, combined = combined)
}

