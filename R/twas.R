# Internal: RAISS-impute GWAS z-scores for LD-sketch variants missing from
# the harmonized GWAS subset. Returns the (possibly widened) sumstats data
# frame. Imputed rows fill `z` from RAISS; `beta` becomes the imputed z and
# `se` becomes 1 when those columns are present in the input. Other columns
# are filled with NA. Imputed variants with R^2 below the threshold are
# dropped by RAISS's internal filter.
imputeMissingGwasForSketch <- function(gwasDataSumstats, sketchRefPanel,
                                       sketchX, imputeOpts, contextLabel = "") {
  missingIds <- setdiff(sketchRefPanel$variant_id, gwasDataSumstats$variant_id)
  if (length(missingIds) == 0) return(gwasDataSumstats)

  refCols <- c("chrom", "pos", "variant_id", "A1", "A2")
  if (!all(refCols %in% colnames(sketchRefPanel))) {
    warning("imputeMissingGwasForSketch: sketch refPanel missing required columns; skipping imputation.")
    return(gwasDataSumstats)
  }
  if (!all(refCols %in% colnames(gwasDataSumstats)) || !"z" %in% colnames(gwasDataSumstats)) {
    warning("imputeMissingGwasForSketch: gwas sumstats missing required columns; skipping imputation.")
    return(gwasDataSumstats)
  }

  # RAISS requires inputs sorted by position (within each chromosome)
  refSorted <- sketchRefPanel[order(sketchRefPanel$chrom, sketchRefPanel$pos), refCols, drop = FALSE]
  knownSorted <- gwasDataSumstats[order(gwasDataSumstats$chrom, gwasDataSumstats$pos), c(refCols, "z"), drop = FALSE]
  # Reorder genotype matrix columns to match the sorted refPanel
  vidOrder <- match(refSorted$variant_id, colnames(sketchX))
  vidOrder <- vidOrder[!is.na(vidOrder)]
  sketchXSorted <- sketchX[, vidOrder, drop = FALSE]
  raissArgs <- c(list(
    refPanel = refSorted,
    knownZscores = knownSorted,
    genotypeMatrix = sketchXSorted,
    verbose = FALSE
  ), imputeOpts)
  raissOut <- tryCatch(do.call(raiss, raissArgs),
                       error = function(e) {
                         warning(sprintf("RAISS missing-variant imputation failed (%s): %s",
                                         contextLabel, e$message))
                         NULL
                       })
  if (is.null(raissOut) || is.null(raissOut$resultFilter)) return(gwasDataSumstats)

  imputedDf <- raissOut$resultFilter
  newRows <- imputedDf[!imputedDf$variant_id %in% gwasDataSumstats$variant_id, , drop = FALSE]
  if (nrow(newRows) == 0) return(gwasDataSumstats)

  added <- newRows[, c("variant_id", "chrom", "pos", "A1", "A2", "z"), drop = FALSE]
  if ("beta" %in% colnames(gwasDataSumstats)) added$beta <- newRows$z
  if ("se"   %in% colnames(gwasDataSumstats)) added$se   <- 1
  for (col in setdiff(colnames(gwasDataSumstats), colnames(added))) {
    added[[col]] <- NA
  }
  added <- added[, colnames(gwasDataSumstats), drop = FALSE]
  message(sprintf("RAISS imputed %d missing GWAS variants (%s).", nrow(added), contextLabel))
  rbind(gwasDataSumstats, added)
}

#' (Deprecated) Harmonize TWAS Weights and GWAS Against LD
#'
#' \strong{Deprecated.} Allele harmonization is now part of
#' \code{\link{summaryStatsQc}} (run it on your
#' \code{\link{GwasSumStats}} before fine-mapping or
#' \code{\link{twasWeightsPipeline}}). The TWAS-side panel
#' harmonization is absorbed into \code{\link{causalInferencePipeline}}.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
harmonizeTwas <- function(...) {
  .Deprecated(new = "causalInferencePipeline", package = "pecotmr",
    msg = paste(
      "harmonizeTwas() has been removed. Allele harmonization is now",
      "part of summaryStatsQc() (for the GWAS side) and",
      "causalInferencePipeline() (for the TWAS-against-LD side)."))
  invisible(NULL)
}

#' (Deprecated) Harmonize GWAS Summary Statistics
#'
#' \strong{Deprecated.} Use \code{\link{summaryStatsQc}} on the
#' \code{\link{GwasSumStats}} input directly; harmonization (allele
#' alignment against the LD sketch via \code{.matchRefPanel}, optional
#' MungeSumstats filters, optional liftover, optional imputation) is
#' all part of that single QC pass.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
harmonizeGwas <- function(...) {
  .Deprecated(new = "summaryStatsQc", package = "pecotmr",
    msg = paste(
      "harmonizeGwas() has been removed. Use summaryStatsQc() on the",
      "GwasSumStats input directly."))
  invisible(NULL)
}

#' (Deprecated) TWAS Pipeline
#'
#' \strong{Deprecated.} Use \code{\link{ctwasPipeline}} for the cTWAS
#' (causal TWAS) variant and \code{\link{causalInferencePipeline}} for
#' the standard TWAS z + MR computation. The new pipelines accept
#' \code{TwasWeights}, \code{FineMappingResult}, and
#' \code{GwasSumStats} S4 inputs and replace the old multi-file
#' driver.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
twasPipeline <- function(...) {
  .Deprecated(new = "ctwasPipeline", package = "pecotmr",
    msg = paste(
      "twasPipeline() has been removed. Use ctwasPipeline() for the",
      "cTWAS variant or causalInferencePipeline() for the TWAS-Z + MR",
      "path."))
  invisible(NULL)
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

