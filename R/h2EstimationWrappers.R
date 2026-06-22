# =============================================================================
# Heritability Estimation Wrappers
# -----------------------------------------------------------------------------
# Consolidated entry points for univariate heritability estimation. Three
# methods are exposed via the `estimateH2` S4 generic:
#
#   * gLDSC (Generalized LD Score Regression, Xiong et al. 2024)
#   * HDL/sHDL (High-Definition Likelihood, Ning et al. 2020 + Zhao 2023)
#   * LDER (LD Eigenvalue Regression, Song et al. 2022)
#
# Plus shared utilities (block ops, WLS, jackknife SE, enrichment, LD
# shrinkage, genome-build checks, tauStar standardization, meta-analysis)
# and the converter that bridges H2Estimate results into the sLDSC
# postprocessing pipeline.
#
# The H2Estimate result class itself lives in R/h2Estimate.R.
# =============================================================================

#' @title Shared Utilities for Heritability Estimation
#' @description Internal helper functions for block operations, regression,
#'   jackknife SE, enrichment computation, and meta-analysis.
#' @name pecotmr-h2-utils
#' @keywords internal
#' @importFrom GenomicRanges GRanges
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
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
  blocksGr <- getBlocks(ldBlocks)
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
#' @description Check that genome builds match between objects. Each
#'   object contributes a single genome build (from its \code{genome}
#'   slot).
#' @param ... Objects with a \code{genome} slot.
#' @return TRUE if all match, error otherwise.
#' @keywords internal
checkGenomeBuild <- function(...) {
  objects <- list(...)
  genomes <- vapply(objects, function(x) {
    if (is(x, "GwasSumStats") || is(x, "QtlSumStats") ||
        is(x, "LdStatistic")  || is(x, "AnnotationMatrix") ||
        is(x, "LdBlocks")) getGenome(x)
    else stop("Unknown object type for genome build check")
  }, character(1))
  if (length(unique(genomes)) > 1L) {
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


#' @title LDER: LD Eigenvalue Regression
#' @description Estimate heritability using LD eigenvalue regression
#'   (Song et al. 2022). Supports univariate global and local estimation,
#'   with optional annotation stratification.
#' @name pecotmr-h2-lder
#' @keywords internal
#' @references
#'   Song S, Jiang W, Zhang Y, Hou L, Zhao H (2022). Leveraging LD
#'   eigenvalue regression to improve the estimation of SNP heritability
#'   and confounding inflation. Am J Hum Genet, 109(5):802-811.
NULL

# =============================================================================
# Univariate LDER
# =============================================================================

#' @title Univariate LDER
#' @description Estimate SNP heritability using LD eigenvalue regression.
#' @param z Numeric vector of z-scores.
#' @param n Numeric, GWAS sample size.
#' @param eigenRef An \code{LdEigen} object.
#' @param annotations An \code{AnnotationMatrix}, or NULL.
#' @param local Logical, return per-block estimates.
#' @param lambda Numeric, ridge penalty (default 0).
#' @return A list with h2, h2Se, intercept, interceptSe, local estimates,
#'   and enrichment estimates.
#' @keywords internal
lderUnivariate <- function(z, n, eigenRef, annotations = NULL,
                           local = FALSE, lambda = 0) {
  eigenList <- getEigenList(eigenRef)
  nBlocks <- length(eigenList)
  nRef <- getNRef(eigenRef)
  inSample <- getInSample(eigenRef)
  M <- nrow(getSnpInfo(eigenRef))

  # Extract baseline annotations if provided
  baselineMat <- NULL
  if (!is.null(annotations)) {
    baseline <- getBaseline(annotations)
    if (ncol(getAnnotations(baseline)) > 0) {
      baselineMat <- getAnnotations(baseline)
    }
  }

  # Collect per-block eigenvalue regression quantities
  blockData <- lapply(seq_len(nBlocks), function(b) {
    block <- eigenList[[b]]
    idx <- block$snpIdx
    d <- block$values        # eigenvalues
    V <- block$vectors       # eigenvectors
    zBlock <- z[idx]

    # Rotate z-scores into eigenbasis
    zRot <- as.vector(t(V) %*% zBlock)
    chi2Rot <- zRot^2

    # Annotation-stratified eigenvalue scores for baseline annotations
    # ldAnnot[i, a] = sum_j V[j,i]^2 * annot[j, a]
    if (!is.null(baselineMat)) {
      ldAnnot <- crossprod(V^2, baselineMat[idx, , drop = FALSE])
    } else {
      ldAnnot <- NULL
    }

    list(
      chi2Rot = chi2Rot,
      eigenvalues = d,
      ldAnnot = ldAnnot,
      n_snps = length(idx),
      snpIdx = idx
    )
  })

  # Assemble regression data
  allChi2 <- unlist(lapply(blockData, `[[`, "chi2Rot"))
  allD <- unlist(lapply(blockData, `[[`, "eigenvalues"))

  # Build design matrix
  # Stratified model: E[chi2_rot_i - 1] = n/M * sum_a(tau_a * d_i * ld_annot_{a,i}) + n*a
  # Unstratified model (no baseline annotations): same with single base column
  if (!is.null(baselineMat)) {
    allLdAnnot <- do.call(rbind, lapply(blockData, `[[`, "ldAnnot"))
    X <- cbind(n * allD * allLdAnnot / M, rep(n, length(allD)))
    nTau <- ncol(baselineMat)
  } else {
    X <- cbind(n * allD / M, rep(n, length(allD)))
    nTau <- 1L
  }

  y <- allChi2 - 1
  w <- 1 / (2 * pmax(allChi2, 1)^2)

  fit <- weightedLsRidge(y, X, w, lambda = lambda, penalizeIntercept = FALSE)
  tau <- fit$coef[seq_len(nTau)]
  a <- fit$coef[nTau + 1]

  # Compute h2 from tau
  if (!is.null(baselineMat)) {
    # h2 = sum_a tau_a * M_a where M_a = sum_j annot_{j,a}
    h2 <- sum(tau * colSums(baselineMat))
  } else {
    h2 <- tau[1]
  }

  # Jackknife SE by block
  blockAssign <- rep(seq_len(nBlocks),
    vapply(blockData, function(x) length(x$eigenvalues), integer(1)))

  looEstimates <- matrix(NA, nrow = nBlocks, ncol = nTau + 1)
  for (b in seq_len(nBlocks)) {
    keep <- blockAssign != b
    fitLoo <- weightedLsRidge(y[keep], X[keep, , drop = FALSE], w[keep],
                              lambda = lambda, penalizeIntercept = FALSE)
    looEstimates[b, ] <- fitLoo$coef
  }

  # Extract per-annotation tau jackknife blocks and SE
  tauBlocks <- looEstimates[, seq_len(nTau), drop = FALSE]
  tauSe <- jackknifeSe(tau, tauBlocks)

  # Compute h2 for each LOO iteration, then jackknife
  if (!is.null(baselineMat)) {
    M_a <- colSums(baselineMat)
    h2Loo <- as.vector(tauBlocks %*% M_a)
  } else {
    h2Loo <- looEstimates[, 1]
  }
  aLoo <- looEstimates[, nTau + 1]
  se <- jackknifeSe(c(h2, a), cbind(h2Loo, aLoo))

  # Baseline enrichment (if annotations provided)
  baselineEnrichmentDf <- NULL
  if (!is.null(baselineMat)) {
    annotNames <- if (!is.null(colnames(baselineMat))) {
      colnames(baselineMat)
    } else {
      paste0("annot_", seq_len(ncol(baselineMat)))
    }
    baselineEnrichmentDf <- computeBaselineEnrichment(
      tau, tauSe, tauBlocks, baselineMat, annotNames, h2
    )
  }

  # Local heritability (if requested)
  localDf <- NULL
  if (local) {
    localDf <- .lderLocalH2(blockData, n, M, tau, a, baselineMat)
  }

  # Score statistics for candidate annotations (if provided)
  scoreStats <- NULL
  if (!is.null(annotations)) {
    strat <- .lderStratified(z, n, eigenRef, annotations, tau, a,
                             baselineMat)
    scoreStats <- strat$scoreStats
  }

  list(
    h2 = h2,
    h2Se = se[1],
    intercept = a,
    interceptSe = se[2],
    tau = tau,
    tauSe = tauSe,
    tauBlocks = tauBlocks,
    local = localDf,
    enrichment = baselineEnrichmentDf,
    scoreStats = scoreStats
  )
}

# =============================================================================
# Internal helpers
# =============================================================================

#' @title LDER local heritability
#' @description Per-block heritability using the Hessian-based SE.
#' @param blockData List of per-block eigenvalue regression quantities.
#' @param n Numeric, GWAS sample size.
#' @param M Integer, total number of SNPs.
#' @param tau Numeric vector of annotation coefficients.
#' @param aGlobal Numeric, global intercept.
#' @param baselineMat Matrix of baseline annotations, or NULL.
#' @return A data.frame with blockId, h2Local, h2LocalSe.
#' @keywords internal
.lderLocalH2 <- function(blockData, n, M, tau, aGlobal,
                         baselineMat = NULL) {
  # Per-block heritability using the Hessian-based SE
  localResults <- lapply(seq_along(blockData), function(b) {
    bd <- blockData[[b]]
    pBlock <- bd$n_snps
    d <- bd$eigenvalues
    chi2 <- bd$chi2Rot

    # Compute fitted baseline contribution for this block
    if (!is.null(baselineMat)) {
      ldAnnot <- bd$ldAnnot  # nEigenvalues x nAnnotations
      fittedBaseline <- as.vector(n / M * d *
                                    (ldAnnot %*% tau))
    } else {
      fittedBaseline <- n * tau[1] * d / M
    }

    # Local regression: residual after removing global baseline + intercept
    y <- chi2 - 1 - n * aGlobal - fittedBaseline
    x <- n * d / M
    if (length(y) < 3) {
      return(data.frame(blockId = b, h2Local = NA, h2LocalSe = NA))
    }
    w <- 1 / (2 * pmax(chi2, 1)^2)
    h2Local <- sum(w * x * y) / sum(w * x^2)

    # Fisher information SE
    info <- sum(w * x^2)
    seLocal <- 1 / sqrt(info)

    data.frame(blockId = b, h2Local = h2Local, h2LocalSe = seLocal)
  })
  do.call(rbind, localResults)
}

#' @title LDER stratified score statistics
#' @description Score-based approach: fit baseline jointly, compute scores
#'   for candidate annotations.
#' @param z Numeric vector of z-scores.
#' @param n Numeric, GWAS sample size.
#' @param eigenRef An \code{LdEigen} object.
#' @param annotations An \code{AnnotationMatrix}, or NULL.
#' @param tau Numeric vector of annotation coefficients.
#' @param a Numeric, intercept.
#' @param baselineMat Matrix of baseline annotations, or NULL.
#' @return A list with enrichment data.frame and scoreStats list.
#' @keywords internal
.lderStratified <- function(z, n, eigenRef, annotations, tau, a,
                            baselineMat = NULL) {
  # Score-based approach: fit baseline jointly, compute scores for candidates
  candidateAnnot <- getCandidates(annotations)
  candMat <- getAnnotations(candidateAnnot)
  nCandidates <- ncol(candMat)

  if (nCandidates == 0) {
    return(list(enrichment = NULL, scoreStats = NULL))
  }

  eigenList <- getEigenList(eigenRef)
  nBlocks <- length(eigenList)
  M <- nrow(getSnpInfo(eigenRef))

  # Collect per-block partial scores into a matrix (nBlocks x nCandidates)
  partialsMat <- matrix(0, nrow = nBlocks, ncol = nCandidates)

  for (b in seq_len(nBlocks)) {
    block <- eigenList[[b]]
    idx <- block$snpIdx
    V <- block$vectors
    d <- block$values
    zBlock <- z[idx]
    zRot <- as.vector(t(V) %*% zBlock)
    chi2Rot <- zRot^2

    # Compute residual from stratified baseline fit
    if (!is.null(baselineMat)) {
      ldAnnotBase <- crossprod(V^2, baselineMat[idx, , drop = FALSE])
      fittedBaseline <- as.vector(n / M * d *
                                    (ldAnnotBase %*% tau))
    } else {
      fittedBaseline <- n * tau[1] * d / M
    }
    residual <- chi2Rot - 1 - fittedBaseline - n * a
    w <- 1 / (2 * pmax(chi2Rot, 1)^2)

    for (ai in seq_len(nCandidates)) {
      annotCol <- candMat[, ai]
      annotBlock <- annotCol[idx]
      ldAnnot <- as.vector(t(V^2) %*% annotBlock)

      partialsMat[b, ai] <- sum(w * residual * n * ldAnnot / M)
    }
  }

  # Compute scoreZ from block partials
  scoreZ <- colSums(partialsMat) /
    sqrt(colSums(partialsMat^2) - colSums(partialsMat)^2 / nBlocks)

  # Score correlation matrix via jackknife
  # For each LOO iteration, recompute scoreZ excluding one block
  looScoreZ <- matrix(0, nrow = nBlocks, ncol = nCandidates)
  for (b in seq_len(nBlocks)) {
    partialsLoo <- partialsMat[-b, , drop = FALSE]
    nLoo <- nBlocks - 1
    looScoreZ[b, ] <- colSums(partialsLoo) /
      sqrt(colSums(partialsLoo^2) - colSums(partialsLoo)^2 / nLoo)
  }
  if (nCandidates > 1) {
    R <- cor(looScoreZ)
  } else {
    R <- matrix(1, 1, 1)
  }

  candMeta <- getAnnotationMeta(candidateAnnot)
  enrichmentDf <- data.frame(
    annotation = candMeta$name,
    scoreZ = scoreZ,
    scoreP = 2 * pnorm(-abs(scoreZ)),
    stringsAsFactors = FALSE
  )

  scoreStatsList <- list(
    z = scoreZ,
    R = R,
    annotationNames = candMeta$name
  )

  list(enrichment = enrichmentDf, scoreStats = scoreStatsList)
}


#' @title g-LDSC: Generalized LD Score Regression
#' @description Estimate heritability and enrichment using feasible
#'   generalized least squares on LD scores (Xiong et al. 2024).
#' @name pecotmr-h2-gldsc
#' @keywords internal
#' @references
#'   Xiong Z, Thach TQ, Zhang YD, Sham PC (2024). Improved estimation
#'   of functional enrichment in SNP heritability using feasible
#'   generalized least squares. HGG Advances, 5(2):100272.
NULL

# =============================================================================
# Univariate g-LDSC
# =============================================================================

#' @title Univariate g-LDSC
#' @description Estimate h2 using g-LDSC with FGLS estimation.
#' @param z Numeric vector of z-scores.
#' @param n Numeric, GWAS sample size.
#' @param ldRef An \code{LdScore} object.
#' @param annotations An \code{AnnotationMatrix}, or NULL.
#' @param local Logical, return per-block estimates.
#' @param lambda Numeric, ridge penalty (default 0).
#' @return A list with h2, h2Se, intercept, enrichment, scoreStats.
#' @keywords internal
gldscUnivariate <- function(z, n, ldRef, annotations = NULL,
                            local = FALSE, lambda = 0) {
  chi2 <- z^2
  M <- length(z)
  ldScores <- getLdScores(ldRef)
  weights <- getLdScoreWeights(ldRef)

  # --- Step 1: Initial S-LDSC estimate (WLS) ---
  # Model: E[chi2_j] = N/M * sum_a(tau_a * l_{j,a}) + N*a + 1
  if (!is.null(annotations)) {
    baseline <- getBaseline(annotations)
    ldStrat <- computeLdScores(ldRef, baseline)
    X <- cbind(n * ldStrat / M, rep(n, M))
  } else {
    X <- cbind(n * ldScores[, 1] / M, rep(n, M))
  }
  y <- chi2

  # Initial WLS (standard S-LDSC)
  fitWls <- weightedLsRidge(y, X, weights, lambda = lambda,
                            penalizeIntercept = FALSE)
  fittedWls <- fitWls$fitted

  # --- Step 2: Estimate residual covariance (FGLS) ---
  # The residual covariance arises from LD between SNPs
  # Cov(resid_j, resid_k) depends on R^2_{jk}
  # Approximate the residual covariance using the LD matrix list
  ldMatrixList <- getLdMatrixList(ldRef)
  if (length(ldMatrixList) > 0) {
    # Use stored LD matrices for FGLS
    OmegaInv <- .computeFglsWeights(ldMatrixList, fittedWls)
  } else {
    # Approximate: use LD scores as proxy for diagonal residual variance
    # This gives an intermediate estimator between S-LDSC and full g-LDSC
    residVar <- 2 * pmax(fittedWls, 1)^2
    OmegaInv <- 1 / residVar
    message("Note: g-LDSC without full LD matrices uses approximate FGLS. ",
            "For full g-LDSC, compute LD reference with ldMatrixList.")
  }

  # --- Step 3: FGLS estimate ---
  if (is.numeric(OmegaInv)) {
    # Diagonal approximation
    fitFgls <- weightedLsRidge(y, X, OmegaInv, lambda = lambda,
                               penalizeIntercept = FALSE)
  } else {
    # Full FGLS with block-diagonal OmegaInv
    fitFgls <- .fglsSolve(y, X, OmegaInv)
  }

  nParams <- ncol(X)
  nTau <- nParams - 1L
  tauFull <- fitFgls$coef[seq_len(nTau)]
  intercept <- fitFgls$coef[nParams]

  # Compute h2 from tau: h2 = sum_a(tau_a * M_a)
  # When annotations are present, computeLdScores prepends a total L2 column
  # (base = all-ones annotation), so tauFull[1] = tauBase with MBase = M,
  # and tauFull[2:nTau] = per-annotation tau.
  if (!is.null(annotations)) {
    baselineMat <- getAnnotations(getBaseline(annotations))
    M_a <- colSums(baselineMat)
    M_a_full <- c(M, M_a)  # base (all-ones) + annotation-specific
    h2 <- sum(tauFull * M_a_full)
    # For enrichment, use only the annotation-specific coefficients
    tau <- tauFull[-1]
  } else {
    baselineMat <- NULL
    M_a_full <- NULL
    h2 <- tauFull[1]  # unstratified: coefficient IS h2
    tau <- tauFull
  }

  # Jackknife SE and LOO estimates
  blockAssign <- .assignSnpsToJackknifeBlocks(ldRef, nBlocks = 200)
  jk <- .gldscJackknife(y, X, OmegaInv, fitFgls$coef, blockAssign,
                        lambda = lambda)

  # Extract per-annotation tau jackknife blocks and SE
  tauBlocksFull <- jk$looEstimates[, seq_len(nTau), drop = FALSE]

  # h2 and intercept SE
  if (!is.null(baselineMat)) {
    h2Loo <- as.vector(tauBlocksFull %*% M_a_full)
    # Annotation-specific blocks (exclude base L2 column)
    tauBlocks <- tauBlocksFull[, -1, drop = FALSE]
  } else {
    h2Loo <- jk$looEstimates[, 1]
    tauBlocks <- tauBlocksFull
  }
  tauSe <- jackknifeSe(tau, tauBlocks)
  aLoo <- jk$looEstimates[, nParams]
  se <- jackknifeSe(c(h2, intercept), cbind(h2Loo, aLoo))

  # Baseline enrichment (annotation-specific only, not base L2)
  baselineEnrichmentDf <- NULL
  if (!is.null(baselineMat)) {
    annotNames <- if (!is.null(colnames(baselineMat))) {
      colnames(baselineMat)
    } else {
      paste0("annot_", seq_len(ncol(baselineMat)))
    }
    baselineEnrichmentDf <- computeBaselineEnrichment(
      tau, tauSe, tauBlocks, baselineMat, annotNames, h2
    )
  }

  # Local estimates
  localDf <- NULL
  if (local) {
    localDf <- .gldscLocal(z, n, ldRef, h2, intercept)
  }

  # Score statistics for candidate annotations
  scoreStats <- NULL
  if (!is.null(annotations)) {
    strat <- .gldscScoreStats(z, n, ldRef, annotations,
                              fitFgls$coef, weights)
    scoreStats <- strat$scoreStats
  }

  list(
    h2 = h2,
    h2Se = se[1],
    intercept = intercept,
    interceptSe = se[2],
    tau = tau,
    tauSe = tauSe,
    tauBlocks = tauBlocks,
    local = localDf,
    enrichment = baselineEnrichmentDf,
    scoreStats = scoreStats
  )
}

# =============================================================================
# Internal helpers
# =============================================================================

#' Compute block-diagonal residual covariance for FGLS
#'
#' @description Compute per-block precision matrices from LD matrices
#'   and fitted values for feasible GLS estimation.
#' @param ldMatrixList List of LD matrix blocks, each with R and snpIdx.
#' @param fittedValues Numeric vector of fitted values from initial WLS.
#' @return A list of per-block precision matrices (inverse Omega blocks).
#' @keywords internal
.computeFglsWeights <- function(ldMatrixList, fittedValues) {
  # Omega_{jk} = 2 * (fitted_j * fitted_k * R^2_{jk} + R^4_{jk})
  result <- vector("list", length(ldMatrixList))
  for (b in seq_along(ldMatrixList)) {
    block <- ldMatrixList[[b]]
    RBlock <- block$R
    snpIdx <- block$snpIdx
    fittedBlock <- fittedValues[snpIdx]
    mBlock <- length(snpIdx)

    R2 <- RBlock^2
    R4 <- R2^2
    # OmegaBlock[j,k] = 2 * (fitted_j * fitted_k * R^2_jk + R^4_jk)
    fittedOuter <- outer(fittedBlock, fittedBlock)
    OmegaBlock <- 2 * (fittedOuter * R2 + R4)

    # Regularise to ensure positive definiteness
    OmegaBlock <- OmegaBlock + diag(1e-6, mBlock)
    OmegaInvBlock <- solve(OmegaBlock)

    result[[b]] <- list(omegaInv = OmegaInvBlock, snpIdx = snpIdx)
  }
  result
}

#' Solve GLS with block-diagonal precision
#'
#' @description Solve the GLS problem beta = (X' Omega^{-1} X)^{-1} X' Omega^{-1} y
#'   for block-diagonal OmegaInv.
#' @param y Numeric vector, response.
#' @param X Numeric matrix, design matrix.
#' @param OmegaInv Block-diagonal precision: either a numeric vector (diagonal)
#'   or a list of per-block precision matrices.
#' @return A list with coef, fitted, and residuals.
#' @keywords internal
.fglsSolve <- function(y, X, OmegaInv) {
  if (is.null(dim(X))) X <- matrix(X, ncol = 1)
  p <- ncol(X)

  if (is.numeric(OmegaInv)) {
    # Diagonal case: standard WLS
    return(weightedLs(y, X, OmegaInv))
  }

  # Block-diagonal case: accumulate across blocks
  XtOiX <- matrix(0, nrow = p, ncol = p)
  XtOiy <- numeric(p)

  for (b in seq_along(OmegaInv)) {
    idx <- OmegaInv[[b]]$snpIdx
    Oi <- OmegaInv[[b]]$omegaInv
    X_b <- X[idx, , drop = FALSE]
    y_b <- y[idx]

    XtOiX <- XtOiX + crossprod(X_b, Oi %*% X_b)
    XtOiy <- XtOiy + crossprod(X_b, Oi %*% y_b)
  }

  beta <- as.vector(solve(XtOiX, XtOiy))
  fittedVals <- as.vector(X %*% beta)
  resid <- y - fittedVals

  list(coef = beta, fitted = fittedVals, residuals = resid)
}

#' Assign SNPs to jackknife blocks
#'
#' @description Divide SNPs into approximately equal-sized jackknife blocks.
#' @param ldRef An \code{LdScore} object.
#' @param nBlocks Integer, number of jackknife blocks (default 200).
#' @return Integer vector of block assignments.
#' @keywords internal
.assignSnpsToJackknifeBlocks <- function(ldRef, nBlocks = 200) {
  nSnps <- nrow(getSnpInfo(ldRef))
  blockSize <- ceiling(nSnps / nBlocks)
  rep(seq_len(nBlocks), each = blockSize, length.out = nSnps)
}

#' Jackknife SE for g-LDSC
#'
#' @description Compute leave-one-block-out jackknife standard errors
#'   for g-LDSC coefficient estimates.
#' @param y Numeric vector, response.
#' @param X Numeric matrix, design matrix.
#' @param w Weights: numeric vector (diagonal) or list of precision blocks.
#' @param coefFull Numeric vector, full-sample coefficients.
#' @param blockAssign Integer vector, block assignments from
#'   \code{.assignSnpsToJackknifeBlocks}.
#' @param lambda Numeric, ridge penalty (default 0).
#' @return A list with se and looEstimates.
#' @keywords internal
.gldscJackknife <- function(y, X, w, coefFull, blockAssign,
                            lambda = 0) {
  nBlocks <- max(blockAssign)
  nParams <- length(coefFull)
  looEstimates <- matrix(NA, nrow = nBlocks, ncol = nParams)

  for (b in seq_len(nBlocks)) {
    keep <- blockAssign != b
    if (is.numeric(w)) {
      fitLoo <- weightedLsRidge(y[keep], X[keep, , drop = FALSE], w[keep],
                                lambda = lambda, penalizeIntercept = FALSE)
    } else {
      fitLoo <- weightedLsRidge(y[keep], X[keep, , drop = FALSE],
                                rep(1, sum(keep)),
                                lambda = lambda, penalizeIntercept = FALSE)
    }
    looEstimates[b, ] <- fitLoo$coef
  }
  list(
    se = jackknifeSe(coefFull, looEstimates),
    looEstimates = looEstimates
  )
}

#' Per-block local h2 from g-LDSC
#'
#' @description Estimate per-block local heritability using the
#'   g-LDSC intercept and LD scores.
#' @param z Numeric vector of z-scores.
#' @param n Numeric, GWAS sample size.
#' @param ldRef An \code{LdScore} object.
#' @param h2 Numeric, global h2 estimate.
#' @param intercept Numeric, g-LDSC intercept estimate.
#' @return A data.frame with blockId, h2Local, h2LocalSe.
#' @keywords internal
.gldscLocal <- function(z, n, ldRef, h2, intercept) {
  ldBlocksObj <- getLdBlocks(ldRef)
  nLdBlocks <- length(getBlocks(ldBlocksObj))
  if (nLdBlocks <= 22) {
    stop("Local g-LDSC requires fine-grained LD blocks (e.g., Berisa & Pickrell loci), ",
         "but the LD reference only has ", nLdBlocks, " blocks (likely per-chromosome). ",
         "Provide explicit ldBlocks when reading the LD reference with readLdRef().")
  }
  chi2 <- z^2
  M <- length(z)
  blockIndices <- snpsPerBlock(getSnpInfo(ldRef), ldBlocksObj)
  nBlocks <- length(blockIndices)

  blockId <- integer(nBlocks)
  h2Local <- numeric(nBlocks)
  h2LocalSe <- numeric(nBlocks)

  for (b in seq_len(nBlocks)) {
    idx <- blockIndices[[b]]
    MBlock <- length(idx)
    if (MBlock == 0) {
      blockId[b] <- b
      h2Local[b] <- NA_real_
      h2LocalSe[b] <- NA_real_
      next
    }
    chi2Block <- chi2[idx]
    ldBlock <- getLdScores(ldRef)[idx, 1]

    sumLd <- sum(ldBlock)
    h2_b <- (sum(chi2Block) - MBlock * (n * intercept + 1)) /
            (n * sumLd / M)

    # SE from block-level Fisher information:
    # Var(h2Local) ~ M_block / (n * sum_ld / M)^2
    h2LocalSeB <- sqrt(2 * MBlock) / (n * sumLd / M)

    blockId[b] <- b
    h2Local[b] <- h2_b
    h2LocalSe[b] <- h2LocalSeB
  }

  data.frame(
    blockId = blockId,
    h2Local = h2Local,
    h2LocalSe = h2LocalSe,
    stringsAsFactors = FALSE
  )
}

#' Score statistics for candidate annotations
#'
#' @description Compute score z-statistics and their correlation matrix
#'   for candidate annotations not included in the baseline model.
#' @param z Numeric vector of z-scores.
#' @param n Numeric, GWAS sample size.
#' @param ldRef An \code{LdScore} object.
#' @param annotations An \code{AnnotationMatrix} object.
#' @param coef Numeric vector, fitted coefficients from the baseline model.
#' @param w Numeric vector, regression weights.
#' @return A list with enrichment (data.frame) and scoreStats (list with
#'   z, R, annotationNames).
#' @keywords internal
.gldscScoreStats <- function(z, n, ldRef, annotations, coef, w) {
  # Score statistics for candidate annotations
  candidate <- getCandidates(annotations)
  candMat <- getAnnotations(candidate)
  nCand <- ncol(candMat)
  if (nCand == 0) return(list(enrichment = NULL, scoreStats = NULL))

  # Compute score for each candidate annotation
  # S_a = dLL/d(tau_a) at tau_a = 0, baseline fitted
  chi2 <- z^2
  M <- length(z)
  ldScores <- getLdScores(ldRef)
  fitted <- chi2 - weightedLs(chi2,
    cbind(n * ldScores[, 1] / M, rep(n, M)), w)$residuals
  resid <- chi2 - fitted

  # Precompute weighted annotation LD scores
  annotLdMat <- matrix(0, nrow = M, ncol = nCand)
  for (a in seq_len(nCand)) {
    annotLdMat[, a] <- ldScores[, 1] * candMat[, a]
  }

  # Full-sample score z-statistics
  scoreZ <- numeric(nCand)
  for (a in seq_len(nCand)) {
    S_a <- sum(w * resid * n * annotLdMat[, a] / M)
    V_a <- sum((w * n * annotLdMat[, a] / M)^2)
    scoreZ[a] <- S_a / sqrt(V_a)
  }

  # Block jackknife for score correlation matrix R
  nBlocks <- 200
  blockAssign <- .assignSnpsToJackknifeBlocks(ldRef, nBlocks = nBlocks)
  nBlocksActual <- max(blockAssign)
  scoreZLoo <- matrix(NA, nrow = nBlocksActual, ncol = nCand)

  for (b in seq_len(nBlocksActual)) {
    keep <- blockAssign != b
    w_b <- w[keep]
    resid_b <- resid[keep]
    for (a in seq_len(nCand)) {
      ald_b <- annotLdMat[keep, a]
      S_a <- sum(w_b * resid_b * n * ald_b / M)
      V_a <- sum((w_b * n * ald_b / M)^2)
      scoreZLoo[b, a] <- S_a / sqrt(V_a)
    }
  }

  # Score correlation from jackknife LOO z-statistics
  if (nCand > 1) {
    R <- cor(scoreZLoo)
    # Ensure valid correlation matrix
    R[is.na(R)] <- 0
  } else {
    R <- matrix(1, 1, 1)
  }

  candMeta <- getAnnotationMeta(candidate)
  list(
    enrichment = data.frame(
      annotation = candMeta$name,
      scoreZ = scoreZ,
      scoreP = 2 * pnorm(-abs(scoreZ)),
      stringsAsFactors = FALSE
    ),
    scoreStats = list(
      z = scoreZ, R = R,
      annotationNames = candMeta$name
    )
  )
}


#' @title HDL/sHDL: High-Definition Likelihood
#' @description Estimate heritability using eigenvalue-based likelihood
#'   maximization (Ning et al. 2020) with stratified extensions (sHDL,
#'   Kim et al. 2023).
#' @name pecotmr-h2-hdl
#' @keywords internal
#' @references
#'   Ning Z, Pawitan Y, Shen X (2020). High-definition likelihood
#'   inference of genetic correlations across human complex traits.
#'   Nat Genet, 52:859-864.
#'
#'   Kim SS, Dey KK, Weissbrod O, et al. (2023). Leveraging LD
#'   eigenvalue regression to improve the estimation of SNP heritability
#'   and functional enrichment. medRxiv.
#' @importFrom stats optimize
NULL

# =============================================================================
# Univariate HDL
# =============================================================================

#' @title Univariate HDL
#' @description Estimate h2 via HDL likelihood.
#' @param z Numeric vector of z-scores.
#' @param n Numeric, GWAS sample size.
#' @param eigenRef An \code{LdEigen} object.
#' @param annotations An \code{AnnotationMatrix}, or NULL.
#' @param local Logical, return per-block estimates.
#' @param lambda Numeric, L2 penalty on tau (default 0).
#' @return List with h2, h2Se, local, enrichment, scoreStats.
#' @keywords internal
hdlUnivariate <- function(z, n, eigenRef, annotations = NULL,
                          local = FALSE, lambda = 0) {
  eigenList <- getEigenList(eigenRef)
  nBlocks <- length(eigenList)
  M <- nrow(getSnpInfo(eigenRef))

  # Extract baseline annotations if provided
  baselineMat <- NULL
  if (!is.null(annotations)) {
    baseline <- getBaseline(annotations)
    if (ncol(getAnnotations(baseline)) > 0) {
      baselineMat <- getAnnotations(baseline)
    }
  }

  # Precompute per-block quantities including annotation-stratified
  # eigenvalue scores
  blockData <- lapply(seq_len(nBlocks), function(b) {
    block <- eigenList[[b]]
    idx <- block$snpIdx
    V <- block$vectors
    d <- block$values
    zRot <- as.vector(t(V) %*% z[idx])

    if (!is.null(baselineMat)) {
      ldAnnot <- crossprod(V^2, baselineMat[idx, , drop = FALSE])
    } else {
      ldAnnot <- NULL
    }

    list(zRot = zRot, d = d, ldAnnot = ldAnnot,
         snpIdx = idx, p = length(idx))
  })

  # Compute per-eigenvalue variance: sigma2_i = n/M * sum_a(tau_a * d_i * ld_annot_{a,i}) + 1
  .computeSigma2 <- function(tau, bd) {
    if (!is.null(bd$ldAnnot)) {
      n / M * bd$d * as.vector(bd$ldAnnot %*% tau) + 1
    } else {
      n * tau[1] * bd$d / M + 1
    }
  }

  if (!is.null(baselineMat)) {
    nTau <- ncol(baselineMat)

    # Negative log-likelihood as function of tau vector (with optional L2 penalty)
    nll <- function(tau) {
      val <- 0
      for (bd in blockData) {
        sigma2 <- .computeSigma2(tau, bd)
        if (any(sigma2 <= 0)) return(1e10)
        val <- val + 0.5 * sum(log(sigma2) + bd$zRot^2 / sigma2)
      }
      if (lambda > 0) val <- val + lambda * sum(tau^2)
      val
    }

    # Gradient of negative log-likelihood
    nllGrad <- function(tau) {
      grad <- numeric(nTau)
      for (bd in blockData) {
        sigma2 <- .computeSigma2(tau, bd)
        if (any(sigma2 <= 0)) return(rep(0, nTau))
        # dsigma2/dtau_a = n/M * d_i * ld_annot_{a,i}
        dsig <- n / M * bd$d * bd$ldAnnot  # (nEigen x nTau)
        # dNLL/dsigma2_i = 0.5 * (1/sigma2_i - z_rot_i^2/sigma2_i^2)
        dNLL_dsig <- 0.5 * (1 / sigma2 - bd$zRot^2 / sigma2^2)
        grad <- grad + as.vector(crossprod(dsig, dNLL_dsig))
      }
      if (lambda > 0) grad <- grad + 2 * lambda * tau
      grad
    }

    # Initialize with uniform h2 across annotations
    tauInit <- rep(0.5 / nTau, nTau)

    opt <- optim(tauInit, nll, gr = nllGrad, method = "BFGS",
                        control = list(maxit = 200, reltol = 1e-8))
    tau <- opt$par

    # h2 = sum_a tau_a * M_a
    h2 <- sum(tau * colSums(baselineMat))
  } else {
    # Single-parameter case: use optimize (current behavior)
    nll <- function(h2) {
      val <- 0
      for (bd in blockData) {
        sigma2 <- n * h2 * bd$d / M + 1
        val <- val + 0.5 * sum(log(sigma2) + bd$zRot^2 / sigma2)
      }
      val
    }

    opt <- optimize(nll, interval = c(-0.5, 1.5), tol = 1e-8)
    h2 <- opt$minimum
    tau <- h2  # scalar, used by downstream functions
  }

  # SE from observed Fisher information (returns h2Se, tauSe, tauVcov)
  fisherResult <- .hdlSeFisherStratified(tau, blockData, n, M,
                                         baselineMat, lambda = lambda)
  h2Se <- fisherResult$h2Se
  tauSe <- fisherResult$tauSe

  # Jackknife tauBlocks via score approximation
  tauBlocks <- NULL
  if (!is.null(baselineMat)) {
    jk <- .hdlJackknifeTau(tau, blockData, n, M, baselineMat,
                           lambda = lambda)
    tauBlocks <- jk$looEstimates
    # Use jackknife SE if available (more robust than Fisher for enrichment)
    tauSe <- jk$tauSe
  }

  # Baseline enrichment
  baselineEnrichmentDf <- NULL
  if (!is.null(baselineMat)) {
    annotNames <- if (!is.null(colnames(baselineMat))) {
      colnames(baselineMat)
    } else {
      paste0("annot_", seq_len(ncol(baselineMat)))
    }
    baselineEnrichmentDf <- computeBaselineEnrichment(
      tau, tauSe, tauBlocks, baselineMat, annotNames, h2
    )
  }

  # Local estimation
  localDf <- NULL
  if (local) {
    localDf <- .hdlLocal(blockData, n, M, tau, baselineMat)
  }

  # Score statistics for candidate annotations (sHDL)
  scoreStats <- NULL
  if (!is.null(annotations)) {
    strat <- .shdlStratified(z, n, eigenRef, annotations, tau,
                             baselineMat)
    scoreStats <- strat$scoreStats
  }

  list(h2 = h2, h2Se = h2Se, intercept = NA_real_,
       interceptSe = NA_real_, tau = tau, tauSe = tauSe,
       tauBlocks = tauBlocks, local = localDf,
       enrichment = baselineEnrichmentDf, scoreStats = scoreStats)
}

# =============================================================================
# Internal helpers
# =============================================================================

#' @keywords internal
.hdlSeFisher <- function(h2, blockData, n, M) {
  # Fisher information = -E[d^2 LL / d h2^2]
  info <- 0
  for (bd in blockData) {
    sigma2 <- n * h2 * bd$d / M + 1
    info <- info + 0.5 * sum((n * bd$d / M)^2 / sigma2^2)
  }
  h2Se <- 1 / sqrt(max(info, 1e-10))
  list(h2Se = h2Se, tauSe = h2Se, tauVcov = NULL)
}

#' @keywords internal
.hdlSeFisherStratified <- function(tau, blockData, n, M,
                                   baselineMat = NULL, lambda = 0) {
  if (is.null(baselineMat)) {
    return(.hdlSeFisher(tau, blockData, n, M))
  }

  # Fisher information matrix for tau, then delta method for h2
  nTau <- length(tau)
  infoMat <- matrix(0, nrow = nTau, ncol = nTau)
  for (bd in blockData) {
    sigma2 <- n / M * bd$d * as.vector(bd$ldAnnot %*% tau) + 1
    # dsigma2/dtau_a = n/M * d_i * ld_annot_{a,i}
    dsig <- n / M * bd$d * bd$ldAnnot  # (nEigen x nTau)
    # Fisher info contribution: sum_i dsig_a * dsig_b / sigma2_i^2
    w <- 1 / sigma2^2
    infoMat <- infoMat + 0.5 * crossprod(dsig * sqrt(w))
  }
  # Add ridge penalty contribution to Hessian
  if (lambda > 0) {
    infoMat <- infoMat + 2 * lambda * diag(nTau)
  }

  # Variance of tau
  tauVcov <- tryCatch(solve(infoMat), error = function(e) {
    diag(1 / pmax(diag(infoMat), 1e-10))
  })
  tauSe <- sqrt(pmax(diag(tauVcov), 0))

  # h2 = sum_a tau_a * M_a, so dh2/dtau = M_a
  M_a <- colSums(baselineMat)
  h2Var <- as.numeric(t(M_a) %*% tauVcov %*% M_a)
  h2Se <- sqrt(max(h2Var, 0))

  list(h2Se = h2Se, tauSe = tauSe, tauVcov = tauVcov)
}

#' @title Block-level jackknife for HDL tau via score approximation
#' @description Approximate leave-one-block-out tau estimates using the
#'   linear approximation: tau_loo_b ~ tau - H^{-1} * grad_b, where H is
#'   the Fisher information (Hessian of NLL) and grad_b is block b's
#'   gradient contribution at the MLE.
#' @keywords internal
.hdlJackknifeTau <- function(tau, blockData, n, M, baselineMat,
                             lambda = 0) {
  nBlocks <- length(blockData)
  nTau <- length(tau)

  infoMat <- matrix(0, nTau, nTau)
  blockGrads <- matrix(0, nBlocks, nTau)

  for (b in seq_len(nBlocks)) {
    bd <- blockData[[b]]
    sigma2 <- n / M * bd$d * as.vector(bd$ldAnnot %*% tau) + 1
    dsig <- n / M * bd$d * bd$ldAnnot
    # Block gradient contribution
    dNLL_dsig <- 0.5 * (1 / sigma2 - bd$zRot^2 / sigma2^2)
    blockGrads[b, ] <- as.vector(crossprod(dsig, dNLL_dsig))
    # Block Fisher info contribution
    w <- 1 / sigma2^2
    infoMat <- infoMat + 0.5 * crossprod(dsig * sqrt(w))
  }
  # Add ridge penalty to Hessian
  if (lambda > 0) {
    infoMat <- infoMat + 2 * lambda * diag(nTau)
  }

  H_inv <- tryCatch(solve(infoMat), error = function(e) {
    diag(1 / pmax(diag(infoMat), 1e-10))
  })

  # LOO estimates via linear approximation at MLE
  looEstimates <- matrix(NA, nBlocks, nTau)
  for (b in seq_len(nBlocks)) {
    looEstimates[b, ] <- tau - as.vector(H_inv %*% blockGrads[b, ])
  }

  list(
    looEstimates = looEstimates,
    tauSe = jackknifeSe(tau, looEstimates)
  )
}

#' @keywords internal
.hdlLocal <- function(blockData, n, M, tau = NULL,
                      baselineMat = NULL) {
  localResults <- lapply(seq_along(blockData), function(b) {
    bd <- blockData[[b]]
    if (bd$p < 3) {
      return(data.frame(blockId = b, h2Local = NA, h2LocalSe = NA))
    }

    # Compute the global baseline sigma2 for this block
    if (!is.null(baselineMat) && !is.null(bd$ldAnnot)) {
      sigma2Baseline <- n / M * bd$d *
        as.vector(bd$ldAnnot %*% tau) + 1
    } else if (!is.null(tau)) {
      sigma2Baseline <- n * tau[1] * bd$d / M + 1
    } else {
      sigma2Baseline <- NULL
    }

    # Local likelihood: optimize a local deviation deltaH2
    # sigma2_i = sigma2_baseline_i + n * delta_h2 * d_i / M
    nllLocal <- function(deltaH2) {
      if (!is.null(sigma2Baseline)) {
        sigma2 <- sigma2Baseline + n * deltaH2 * bd$d / M
      } else {
        sigma2 <- n * deltaH2 * bd$d / M + 1
      }
      if (any(sigma2 <= 0)) return(1e10)
      0.5 * sum(log(sigma2) + bd$zRot^2 / sigma2)
    }
    opt <- optimize(nllLocal, interval = c(-0.5, 0.5), tol = 1e-8)
    deltaH2 <- opt$minimum

    # Total local h2 = global baseline contribution + local deviation
    if (!is.null(sigma2Baseline)) {
      # The baseline already explains some h2 in this block; delta is
      # the additional local deviation
      h2Local <- deltaH2
    } else {
      h2Local <- deltaH2
    }

    # SE from Fisher information at optimum
    sigma2Opt <- if (!is.null(sigma2Baseline)) {
      sigma2Baseline + n * deltaH2 * bd$d / M
    } else {
      n * deltaH2 * bd$d / M + 1
    }
    info <- 0.5 * sum((n * bd$d / M)^2 / sigma2Opt^2)
    seLocal <- 1 / sqrt(max(info, 1e-10))

    data.frame(blockId = b, h2Local = h2Local, h2LocalSe = seLocal)
  })
  do.call(rbind, localResults)
}

#' @keywords internal
.shdlStratified <- function(z, n, eigenRef, annotations, tau,
                            baselineMat = NULL) {
  # sHDL: stratified HDL with annotation-specific h2
  # Score-based approach: baseline fit uses jointly fitted tau
  candidate <- getCandidates(annotations)
  candMat <- getAnnotations(candidate)
  nCand <- ncol(candMat)
  if (nCand == 0) return(list(enrichment = NULL, scoreStats = NULL))

  M <- nrow(getSnpInfo(eigenRef))
  eigenList <- getEigenList(eigenRef)
  nBlocks <- length(eigenList)

  # Score for each candidate: derivative of log-likelihood w.r.t. tau_a
  # at tau_a = 0, with baseline tau at the MLE
  # Collect per-block gradient and information for jackknife
  blockGrad <- matrix(0, nrow = nBlocks, ncol = nCand)
  blockInfo <- matrix(0, nrow = nBlocks, ncol = nCand)

  for (a in seq_len(nCand)) {
    annotCol <- candMat[, a]

    for (b in seq_len(nBlocks)) {
      block <- eigenList[[b]]
      idx <- block$snpIdx
      V <- block$vectors
      d <- block$values
      zRot <- as.vector(t(V) %*% z[idx])

      # Annotation-stratified eigenvalue score for this candidate
      annotBlock <- annotCol[idx]
      dAnnot <- as.vector(t(V^2) %*% annotBlock)

      # Compute sigma2 from stratified baseline fit
      if (!is.null(baselineMat)) {
        ldAnnotBase <- crossprod(V^2, baselineMat[idx, , drop = FALSE])
        sigma2 <- n / M * d * as.vector(ldAnnotBase %*% tau) + 1
      } else {
        sigma2 <- n * tau[1] * d / M + 1
      }

      # Gradient: sum_i [ (z_rot_i^2 / sigma2_i - 1) * n * d_annot_i / M ] /
      #           (2 * sigma2_i)
      grad <- 0.5 * sum((zRot^2 / sigma2 - 1) * n * dAnnot /
                          (M * sigma2))
      info <- 0.5 * sum((n * dAnnot / M)^2 / sigma2^2)

      blockGrad[b, a] <- grad
      blockInfo[b, a] <- info
    }
  }

  # Total score statistics
  totalGrad <- colSums(blockGrad)
  totalInfo <- colSums(blockInfo)
  scoreZ <- totalGrad / sqrt(pmax(totalInfo, 1e-10))

  # Jackknife score correlation R: LOO by block
  looScoreZ <- matrix(0, nrow = nBlocks, ncol = nCand)
  for (b in seq_len(nBlocks)) {
    looGrad <- totalGrad - blockGrad[b, ]
    looInfo <- totalInfo - blockInfo[b, ]
    looScoreZ[b, ] <- looGrad / sqrt(pmax(looInfo, 1e-10))
  }
  if (nCand > 1 && nBlocks > 2) {
    R <- cor(looScoreZ)
  } else {
    R <- diag(nCand)
  }

  candMeta <- getAnnotationMeta(candidate)
  list(
    enrichment = data.frame(
      annotation = candMeta$name,
      scoreZ = scoreZ,
      scoreP = 2 * pnorm(-abs(scoreZ)),
      stringsAsFactors = FALSE
    ),
    scoreStats = list(z = scoreZ, R = R,
                       annotationNames = candMeta$name)
  )
}


#' @title Heritability Estimation Entry Points and Converters
#' @description Top-level entry point for heritability estimation,
#'   LD score computation methods, H2Estimate accessors, and a converter
#'   to bridge H2Estimate into the sldscWrapper.R postprocessing pipeline.
#' @name pecotmr-h2-wrappers
#' @keywords internal
#' @include AllGenerics.R
#' @importFrom stats median
NULL

# =============================================================================
# estimateH2 — main dispatch
# =============================================================================

#' @rdname estimateH2
#' @export
setMethod("estimateH2",
  signature(sumstats = "GwasSumStats", ldRef = "LdStatistic"),
  function(sumstats, ldRef, method = "lder", annotations = NULL,
           local = FALSE, study = NULL, ...) {

    method <- match.arg(method, c("lder", "gldsc", "hdl"))
    .validateMethodRef(method, ldRef)

    # Resolve which study's stats to operate on. Defaults to the single
    # entry when the collection has one row.
    if (is.null(study)) {
      if (nrow(sumstats) != 1L) {
        stop("`study` is required when the GwasSumStats has ",
             nrow(sumstats), " entries.")
      }
      study <- as.character(sumstats$study[[1L]])
    }

    z <- getZ(sumstats, study = study)
    n <- median(getN(sumstats, study = study))
    M <- nSnps(sumstats, study = study)

    # Apply the legacy heritability-wrapper correction. This is separate from
    # the SuSiE RSS binaryTraitModel handling in the fine-mapping pipeline.
    varY <- getVarY(sumstats, study = study)
    if (!is.null(varY)) {
      n <- n / varY
    }

    # Dispatch to method-specific function
    result <- switch(method,
      "lder" = lderUnivariate(z, n, ldRef, annotations, local, ...),
      "gldsc" = gldscUnivariate(z, n, ldRef, annotations, local, ...),
      "hdl" = hdlUnivariate(z, n, ldRef, annotations, local, ...)
    )

    # Wrap into H2Estimate S4 object
    new("H2Estimate",
      h2 = result$h2,
      h2Se = result$h2Se,
      intercept = result$intercept %||% NA_real_,
      interceptSe = result$interceptSe %||% NA_real_,
      local = result$local,
      enrichment = result$enrichment,
      tauBlocks = result$tauBlocks,
      scoreStats = result$scoreStats,
      method = method,
      nSnps = as.integer(M),
      traitName = study
    )
  }
)

#' @keywords internal
.validateMethodRef <- function(method, ldRef) {
  if (method %in% c("lder", "hdl") && !is(ldRef, "LdEigen")) {
    stop("Method '", method, "' requires an LdEigen object, ",
         "got ", class(ldRef))
  }
  if (method == "gldsc" && !is(ldRef, "LdScore")) {
    stop("Method 'gldsc' requires an LdScore object, ",
         "got ", class(ldRef))
  }
  invisible(TRUE)
}

# =============================================================================
# computeLdScores — LD score computation
# =============================================================================

#' @rdname computeLdScores
#' @export
setMethod("computeLdScores",
  signature(ldRef = "LdEigen"),
  function(ldRef, annotations = NULL, ...) {
    # Reconstruct LD scores from eigendecompositions
    # l2[j] = sum_k r^2_{jk} = sum_b sum_{eigenvalues in b} V[j,.]^2 * d
    nSnps <- nrow(getSnpInfo(ldRef))
    eigenList <- getEigenList(ldRef)

    if (is.null(annotations)) {
      # Base LD scores only
      l2 <- numeric(nSnps)
      for (b in seq_along(eigenList)) {
        block <- eigenList[[b]]
        idx <- block$snpIdx
        V <- block$vectors
        d <- block$values
        # LD score for SNP j = sum_i V[j,i]^2 * d[i]^2
        # (since R = V D V', R^2_{jk} = sum_i V[j,i]^2 * d[i]^2 * V[k,i]^2)
        # Simplified: l2[j] = sum_i (V[j,i] * d[i])^2
        Vd <- sweep(V, 2, d, "*")
        l2[idx] <- rowSums(Vd^2)
      }
      return(matrix(l2, ncol = 1,
                     dimnames = list(NULL, "base_l2")))
    }

    # Annotation-stratified LD scores
    annotMat <- getAnnotations(annotations)
    nAnnot <- ncol(annotMat)
    # Base + annotation-stratified columns
    l2Strat <- matrix(0, nrow = nSnps, ncol = 1 + nAnnot)

    for (b in seq_along(eigenList)) {
      block <- eigenList[[b]]
      idx <- block$snpIdx
      V <- block$vectors
      d <- block$values

      # Base LD score
      Vd <- sweep(V, 2, d, "*")
      l2Strat[idx, 1] <- rowSums(Vd^2)

      # Annotation-stratified: l2_a[j] = sum_k r^2_{jk} * annot[k,a]
      # Using eigendecomposition: l2_a[j] = sum_i V[j,i]^2 * d[i]^2 *
      #   (sum_k V[k,i]^2 * annot[k,a])
      for (a in seq_len(nAnnot)) {
        annotCol <- annotMat[idx, a]
        # For each eigenvalue: weight = sum_k V[k,i]^2 * annot[k,a]
        annotWeights <- as.vector(crossprod(V^2, annotCol))
        l2Strat[idx, 1 + a] <- as.vector(Vd^2 %*% annotWeights)
      }
    }

    colNames <- c("base_l2", getAnnotationMeta(annotations)$name)
    colnames(l2Strat) <- colNames
    l2Strat
  }
)

#' @rdname computeLdScores
#' @export
setMethod("computeLdScores",
  signature(ldRef = "LdScore"),
  function(ldRef, annotations = NULL, ...) {
    if (is.null(annotations)) {
      return(getLdScores(ldRef))
    }

    # Compute annotation-stratified LD scores using LD matrices
    ldMatrixList <- getLdMatrixList(ldRef)
    if (length(ldMatrixList) == 0) {
      stop("Annotation-stratified LD scores require ldMatrixList in LdScore. ",
           "Recompute the LD reference with full LD matrices.")
    }

    nSnps <- nrow(getSnpInfo(ldRef))
    annotMat <- getAnnotations(annotations)
    nAnnot <- ncol(annotMat)

    # Base L2 + annotation-stratified columns
    l2Strat <- matrix(0, nrow = nSnps, ncol = 1 + nAnnot)
    l2Strat[, 1] <- getLdScores(ldRef)[, 1]

    for (b in seq_along(ldMatrixList)) {
      block <- ldMatrixList[[b]]
      R <- block$R
      idx <- block$snpIdx
      R2 <- R^2
      for (a in seq_len(nAnnot)) {
        # l2_a[j] = sum_k R^2_{jk} * annot[k, a]
        l2Strat[idx, 1 + a] <- as.vector(R2 %*% annotMat[idx, a])
      }
    }

    colNames <- c("base_l2", getAnnotationMeta(annotations)$name)
    colnames(l2Strat) <- colNames
    l2Strat
  }
)

# H2Estimate accessor methods (getLocal/getEnrichment/getScoreStats) live
# in R/h2Estimate.R alongside the class definition.

# =============================================================================
# Converter: H2Estimate -> sldsc_wrapper list format
# =============================================================================

#' @title Convert H2Estimate to S-LDSC Trait Format
#' @description Convert an \code{H2Estimate} object into the list format
#'   expected by \code{\link{standardizeSldscTrait}} and
#'   \code{\link{metaSldscRandom}}. This bridges the h2 estimation
#'   methods (LDER, gLDSC, HDL) into the sldscWrapper.R postprocessing
#'   pipeline.
#' @param h2Est An \code{H2Estimate} object with enrichment and tauBlocks.
#' @return A named list matching the format of \code{\link{readSldscTrait}}:
#'   \describe{
#'     \item{categories}{Character vector of annotation names}
#'     \item{tau}{Named numeric vector of per-annotation coefficients}
#'     \item{tauSe}{Named numeric vector of tau standard errors}
#'     \item{enrichment}{Named numeric vector of enrichment ratios}
#'     \item{enrichmentSe}{Named numeric vector of enrichment SEs}
#'     \item{enrichmentP}{Named numeric vector of enrichment p-values}
#'     \item{propH2}{Named numeric vector of proportion of h2}
#'     \item{propSnps}{Named numeric vector of proportion of SNPs}
#'     \item{h2g}{Numeric scalar, global h2 estimate}
#'     \item{tauBlocks}{Matrix (nBlocks x nCategories) for jackknife}
#'     \item{nBlocks}{Integer, number of jackknife blocks}
#'   }
#' @export
h2EstimateToSldscTrait <- function(h2Est) {
  if (!is(h2Est, "H2Estimate")) {
    stop("h2Est must be an H2Estimate object")
  }

  enrichDf <- getEnrichment(h2Est)
  if (is.null(enrichDf)) {
    stop("H2Estimate has no enrichment results. ",
         "Run estimateH2 with annotations to get enrichment estimates.")
  }

  cats <- as.character(enrichDf$annotation)
  nCats <- length(cats)

  tauBlocks <- getTauBlocks(h2Est)
  if (is.null(tauBlocks)) {
    # Create a dummy single-block matrix from the point estimates
    tauBlocks <- matrix(enrichDf$tau, nrow = 1)
    colnames(tauBlocks) <- cats
    nBlocks <- 1L
  } else {
    nBlocks <- nrow(tauBlocks)
    if (is.null(colnames(tauBlocks))) {
      colnames(tauBlocks) <- cats
    }
  }

  list(
    categories    = cats,
    tau           = setNames(enrichDf$tau, cats),
    tauSe        = setNames(enrichDf$tauSe, cats),
    enrichment    = setNames(enrichDf$enrichment, cats),
    enrichmentSe = setNames(enrichDf$enrichmentSe, cats),
    enrichmentP  = setNames(enrichDf$enrichmentP, cats),
    propH2       = setNames(enrichDf$propH2, cats),
    propSnps     = setNames(enrichDf$propSnps, cats),
    h2g           = getH2(h2Est),
    tauBlocks    = tauBlocks,
    nBlocks      = nBlocks
  )
}
