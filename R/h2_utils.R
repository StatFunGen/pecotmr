#' @title Shared Utilities for Heritability Estimation
#' @description Internal helper functions for block operations, regression,
#'   jackknife SE, enrichment computation, and meta-analysis.
#' @name pecotmr-h2-utils
#' @keywords internal
#' @importFrom GenomicRanges GRanges
#' @importFrom BiocParallel bplapply bpparam
NULL

# =============================================================================
# Block-level operations
# =============================================================================

#' @title Get SNP Indices Per Block
#' @description For each LD block, find the SNP indices from a reference
#'   that fall within the block boundaries.
#' @param snpInfo A data.frame with columns CHR, BP.
#' @param ldBlocks An \code{LdBlocks} object.
#' @return A list of integer vectors, one per block.
#' @keywords internal
snpsPerBlock <- function(snpInfo, ldBlocks) {
  blocksGr <- ldBlocks@blocks
  snpGr <- GRanges(
    seqnames = snpInfo$CHR,
    ranges = IRanges(start = snpInfo$BP, width = 1L)
  )
  hits <- findOverlaps(snpGr, blocksGr)
  split(queryHits(hits), subjectHits(hits))
}

#' @title Apply Function Per Block with BiocParallel
#' @description Apply a function to each LD block in parallel.
#' @param blockIndices List of SNP index vectors per block.
#' @param FUN Function to apply to each block's indices.
#' @param BPPARAM BiocParallel parameter object.
#' @param ... Additional arguments passed to FUN.
#' @return A list of results, one per block.
#' @keywords internal
bplapplyBlocks <- function(blockIndices, FUN, BPPARAM = NULL, ...) {
  if (is.null(BPPARAM)) {
    BPPARAM <- bpparam()
  }
  bplapply(blockIndices, FUN, BPPARAM = BPPARAM, ...)
}

# =============================================================================
# Regression utilities
# =============================================================================

#' @title Weighted Least Squares
#' @description Compute WLS estimate with standard errors.
#' @param y Numeric vector, response.
#' @param X Numeric matrix, predictors.
#' @param w Numeric vector, weights (inverse variance).
#' @return A list with coefficients, SE, residuals, fitted values.
#' @keywords internal
weightedLs <- function(y, X, w) {
  if (is.null(dim(X))) X <- matrix(X, ncol = 1)
  W <- diag(sqrt(w))
  Xw <- W %*% X
  yw <- W %*% y
  XtX <- crossprod(Xw)
  Xty <- crossprod(Xw, yw)
  coef <- solve(XtX, Xty)
  fitted <- X %*% coef
  resid <- y - fitted
  # Heteroskedasticity-robust SE (HC0)
  meat <- crossprod(Xw * as.vector(resid))
  bread <- solve(XtX)
  vcov <- bread %*% meat %*% bread
  se <- sqrt(diag(vcov))
  list(coef = as.vector(coef), se = se, residuals = as.vector(resid),
       fitted = as.vector(fitted), vcov = vcov)
}

#' @title Jackknife Standard Errors by Block
#' @description Compute jackknife SE estimates using leave-one-block-out.
#' @param estimatesFull Numeric vector, full-sample parameter estimates.
#' @param estimatesLoo A matrix (nBlocks x nParams), leave-one-out estimates.
#' @return Numeric vector of jackknife SEs.
#' @keywords internal
jackknifeSe <- function(estimatesFull, estimatesLoo) {
  nBlocks <- nrow(estimatesLoo)
  pseudoVals <- nBlocks * matrix(estimatesFull, nrow = nBlocks,
                                 ncol = length(estimatesFull),
                                 byrow = TRUE) -
    (nBlocks - 1) * estimatesLoo
  jkVar <- apply(pseudoVals, 2, var) / nBlocks
  sqrt(jkVar)
}

# =============================================================================
# Ridge-regularized WLS
# =============================================================================

#' @title Ridge-Regularized Weighted Least Squares
#' @description WLS with optional L2 penalty on coefficients.
#' @param y Numeric vector, response.
#' @param X Numeric matrix, predictors.
#' @param w Numeric vector, weights (inverse variance).
#' @param lambda Numeric, ridge penalty. 0 = no penalty (delegates to
#'   \code{weightedLs}).
#' @param penalizeIntercept Logical. If FALSE (default), the last column
#'   of X (assumed to be the intercept) is not penalized.
#' @return Same structure as \code{weightedLs}: coef, se, residuals, fitted,
#'   vcov.
#' @keywords internal
weightedLsRidge <- function(y, X, w, lambda = 0,
                            penalizeIntercept = FALSE) {
  if (lambda == 0) return(weightedLs(y, X, w))
  if (is.null(dim(X))) X <- matrix(X, ncol = 1)
  p <- ncol(X)
  W <- diag(sqrt(w))
  Xw <- W %*% X
  yw <- W %*% y
  XtX <- crossprod(Xw)
  Xty <- crossprod(Xw, yw)
  # Ridge penalty matrix (don't penalize intercept by default)
  penalty <- diag(lambda, p)
  if (!penalizeIntercept && p > 1) penalty[p, p] <- 0
  coef <- solve(XtX + penalty, Xty)
  fitted <- X %*% coef
  resid <- y - fitted
  # Sandwich SE accounting for ridge shrinkage
  bread <- solve(XtX + penalty)
  meat <- crossprod(Xw * as.vector(resid))
  vcov <- bread %*% meat %*% bread
  se <- sqrt(pmax(diag(vcov), 0))
  list(coef = as.vector(coef), se = se, residuals = as.vector(resid),
       fitted = as.vector(fitted), vcov = vcov)
}

# =============================================================================
# Baseline enrichment computation
# =============================================================================

#' @title Compute Baseline Annotation Enrichment Quantities
#' @description Given tau coefficients and a baseline annotation matrix,
#'   compute the full set of enrichment quantities: propH2, propSnps,
#'   enrichment ratio, enrichment SE (from jackknife or delta method),
#'   and p-value.
#' @param tau Numeric vector of per-annotation regression coefficients.
#' @param tauSe Numeric vector of SE for tau.
#' @param tauBlocks Numeric matrix (nBlocks x nAnnotations) of jackknife
#'   block-level tau values, or NULL.
#' @param baselineMat Numeric matrix (nSnps x nAnnotations).
#' @param annotNames Character vector of annotation names.
#' @param h2 Numeric scalar, total estimated h2.
#' @return A data.frame with columns: annotation, tau, tauSe, enrichment,
#'   enrichmentSe, enrichmentP, propH2, propSnps.
#' @keywords internal
computeBaselineEnrichment <- function(tau, tauSe, tauBlocks,
                                      baselineMat, annotNames, h2) {
  M <- nrow(baselineMat)
  M_a <- colSums(baselineMat)
  propSnps <- M_a / M

  # Per-annotation h2 and proportion
  h2_a <- tau * M_a
  propH2 <- h2_a / h2

  # Enrichment ratio: (propH2 / propSnps) = tau * M / h2
  enrichment <- tau * M / h2

  # Enrichment SE from jackknife blocks (preferred) or delta method (fallback)
  if (!is.null(tauBlocks)) {
    nBlocks <- nrow(tauBlocks)
    # Per-block enrichment: enrichment_b = tau_b * M / h2_b
    h2Blocks <- as.vector(tauBlocks %*% M_a)
    # Avoid division by zero
    h2Blocks[h2Blocks == 0] <- NA
    enrichmentBlocks <- sweep(tauBlocks, 1, h2Blocks, FUN = "/") * M
    # Jackknife variance: Var = (B-1)/B * sum((x_b - x_bar)^2)
    enrichmentMean <- colMeans(enrichmentBlocks, na.rm = TRUE)
    enrichmentVar <- (nBlocks - 1) / nBlocks *
      colSums(sweep(enrichmentBlocks, 2, enrichmentMean)^2, na.rm = TRUE)
    enrichmentSe <- sqrt(enrichmentVar)
  } else {
    # Delta method fallback: d(enrichment)/d(tau) = M / h2
    enrichmentSe <- tauSe * M / abs(h2)
  }

  # P-value from z-score
  enrichmentZ <- enrichment / enrichmentSe
  enrichmentP <- 2 * pnorm(-abs(enrichmentZ))

  data.frame(
    annotation = annotNames,
    tau = tau,
    tauSe = tauSe,
    enrichment = enrichment,
    enrichmentSe = enrichmentSe,
    enrichmentP = enrichmentP,
    propH2 = propH2,
    propSnps = propSnps,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# LD shrinkage
# =============================================================================

#' @title Apply LD Shrinkage
#' @description Apply shrinkage to sample LD matrix to reduce noise from
#'   finite reference panel size, following Wen & Stephens (2010).
#' @param R Numeric matrix, sample LD correlation matrix.
#' @param nRef Integer, reference panel sample size.
#' @param shrinkageType Character, one of "wen_stephens", "constant".
#' @param geneticMap Numeric vector, genetic map positions for SNPs in R.
#' @return Shrunk LD correlation matrix.
#' @keywords internal
shrinkLd <- function(R, nRef, shrinkageType = "wen_stephens",
                     geneticMap = NULL) {
  if (shrinkageType == "wen_stephens" && !is.null(geneticMap)) {
    # Wen & Stephens (2010) shrinkage based on genetic distance
    p <- nrow(R)
    theta <- 2 * nRef / (22 * nRef + 16)  # effective recombination
    distCm <- abs(outer(geneticMap, geneticMap, "-"))
    shrinkFactor <- exp(-4 * nRef * distCm / (100 * (2 * nRef + 16)))
    RShrunk <- R * shrinkFactor
    diag(RShrunk) <- 1
  } else {
    # Simple constant shrinkage
    lambda <- 1 / sqrt(nRef)
    RShrunk <- (1 - lambda) * R + lambda * diag(nrow(R))
  }
  RShrunk
}

# =============================================================================
# Genome build utilities
# =============================================================================

#' @title Validate Genome Build Consistency
#' @description Check that genome builds match between objects.
#' @param ... Objects with a \code{genome} slot.
#' @return TRUE if all match, error otherwise.
#' @keywords internal
checkGenomeBuild <- function(...) {
  objects <- list(...)
  genomes <- vapply(objects, function(x) {
    if (is(x, "GwasSumStats")) x@genome
    else if (is(x, "LdStatistic")) x@genome
    else if (is(x, "AnnotationMatrix")) x@genome
    else if (is(x, "LdBlocks")) x@genome
    else stop("Unknown object type for genome build check")
  }, character(1))
  if (length(unique(genomes)) > 1) {
    stop("Genome build mismatch: ", paste(genomes, collapse = ", "))
  }
  invisible(TRUE)
}

# =============================================================================
# Gazal tau* standardization
# =============================================================================

#' @title Standardize Tau to Tau-Star (Gazal et al. 2017)
#' @description Compute the Gazal-standardized per-annotation effect
#'   \eqn{\tau^*_C = \tau_C \cdot sd_C \cdot M_{ref} / h^2_g}, with
#'   jackknife SE from block-level tau values.
#' @param tau Numeric vector of per-annotation regression coefficients.
#' @param tauBlocks Numeric matrix (nBlocks x nAnnotations) of
#'   block-level tau estimates from delete-one jackknife.
#' @param sdAnnot Numeric vector of per-annotation standard deviations,
#'   same length as \code{tau}.
#' @param MRef Scalar integer, total number of reference-panel SNPs.
#' @param h2g Numeric scalar, total estimated SNP heritability.
#' @return A list with:
#'   \describe{
#'     \item{tauStar}{Numeric vector of standardized tau values.}
#'     \item{tauStarSe}{Numeric vector of jackknife SE for tauStar.}
#'   }
#' @keywords internal
standardizeTauStar <- function(tau, tauBlocks, sdAnnot, MRef, h2g) {
  if (length(tau) != length(sdAnnot)) {
    stop("standardizeTauStar: tau and sdAnnot must have the same length.")
  }
  if (h2g == 0) {
    stop("standardizeTauStar: h2g must be non-zero.")
  }

  # Gazal standardization: tau* = tau * sdAnnot * MRef / h2g
  coef <- sdAnnot * MRef / h2g
  tauStar <- tau * coef

  # Jackknife SE from block-level tau
  tauStarBlocks <- sweep(tauBlocks, 2L, coef, FUN = "*")
  nBlocks <- nrow(tauStarBlocks)
  jkVar <- apply(tauStarBlocks, 2L, function(x) var(x, na.rm = TRUE))
  tauStarSe <- sqrt((nBlocks - 1)^2 / nBlocks * jkVar)

  list(tauStar = tauStar, tauStarSe = tauStarSe)
}

# =============================================================================
# DerSimonian-Laird random-effects meta-analysis
# =============================================================================

#' @title Random-Effects Meta-Analysis (DerSimonian-Laird)
#' @description Perform a DerSimonian-Laird random-effects meta-analysis
#'   from a set of study-level point estimates and standard errors.
#' @param means Numeric vector of study-level point estimates.
#' @param ses Numeric vector of study-level standard errors (must be
#'   positive and finite).
#' @return A list with:
#'   \describe{
#'     \item{mean}{Pooled meta-analytic mean.}
#'     \item{se}{Standard error of the pooled mean.}
#'     \item{tau2}{Estimated between-study variance.}
#'     \item{I2}{Higgins I-squared heterogeneity statistic (proportion of
#'       total variance due to between-study variance), in [0, 1].}
#'     \item{Q}{Cochran's Q statistic for heterogeneity.}
#'   }
#' @keywords internal
metaRandomEffects <- function(means, ses) {
  k <- length(means)
  if (k != length(ses)) {
    stop("metaRandomEffects: means and ses must have the same length.")
  }
  if (k == 0L) {
    return(list(mean = NA_real_, se = NA_real_, tau2 = NA_real_,
                I2 = NA_real_, Q = NA_real_))
  }
  if (k == 1L) {
    return(list(mean = means[1], se = ses[1], tau2 = 0,
                I2 = 0, Q = 0))
  }
  if (any(!is.finite(ses) | ses <= 0)) {
    stop("metaRandomEffects: all ses must be positive and finite.")
  }

  # Fixed-effect weights
  wFe <- 1 / ses^2

  # Fixed-effect pooled estimate
  muFe <- sum(wFe * means) / sum(wFe)

  # Cochran's Q
  Q <- sum(wFe * (means - muFe)^2)

  # DerSimonian-Laird tau-squared estimator
  cDl <- sum(wFe) - sum(wFe^2) / sum(wFe)
  tau2 <- max(0, (Q - (k - 1)) / cDl)

  # Random-effects weights

  wRe <- 1 / (ses^2 + tau2)

  # Pooled random-effects estimate
  muRe <- sum(wRe * means) / sum(wRe)
  seRe <- sqrt(1 / sum(wRe))

  # Higgins I-squared
  I2 <- max(0, (Q - (k - 1)) / Q)

  list(mean = muRe, se = seRe, tau2 = tau2, I2 = I2, Q = Q)
}
