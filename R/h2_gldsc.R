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
#' @return A list with h2, h2_se, intercept, enrichment, score_stats.
#' @keywords internal
gldscUnivariate <- function(z, n, ldRef, annotations = NULL,
                            local = FALSE, lambda = 0) {
  chi2 <- z^2
  M <- length(z)
  ldScores <- ldRef@ldScores
  weights <- ldRef@ldScoreWeights

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
  if (length(ldRef@ldMatrixList) > 0) {
    # Use stored LD matrices for FGLS
    OmegaInv <- .computeFglsWeights(ldRef@ldMatrixList, fittedWls)
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
    baselineMat <- getBaseline(annotations)@annotations
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
  tauBlocksFull <- jk$loo_estimates[, seq_len(nTau), drop = FALSE]

  # h2 and intercept SE
  if (!is.null(baselineMat)) {
    h2Loo <- as.vector(tauBlocksFull %*% M_a_full)
    # Annotation-specific blocks (exclude base L2 column)
    tauBlocks <- tauBlocksFull[, -1, drop = FALSE]
  } else {
    h2Loo <- jk$loo_estimates[, 1]
    tauBlocks <- tauBlocksFull
  }
  tauSe <- jackknifeSe(tau, tauBlocks)
  aLoo <- jk$loo_estimates[, nParams]
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
    scoreStats <- strat$score_stats
  }

  list(
    h2 = h2,
    h2_se = se[1],
    intercept = intercept,
    intercept_se = se[2],
    tau = tau,
    tau_se = tauSe,
    tau_blocks = tauBlocks,
    local = localDf,
    enrichment = baselineEnrichmentDf,
    score_stats = scoreStats
  )
}

# =============================================================================
# Internal helpers
# =============================================================================

#' Compute block-diagonal residual covariance for FGLS
#'
#' @description Compute per-block precision matrices from LD matrices
#'   and fitted values for feasible GLS estimation.
#' @param ldMatrixList List of LD matrix blocks, each with R and snp_idx.
#' @param fittedValues Numeric vector of fitted values from initial WLS.
#' @return A list of per-block precision matrices (inverse Omega blocks).
#' @keywords internal
.computeFglsWeights <- function(ldMatrixList, fittedValues) {
  # Omega_{jk} = 2 * (fitted_j * fitted_k * R^2_{jk} + R^4_{jk})
  result <- vector("list", length(ldMatrixList))
  for (b in seq_along(ldMatrixList)) {
    block <- ldMatrixList[[b]]
    RBlock <- block$R
    snpIdx <- block$snp_idx
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

    result[[b]] <- list(Omega_inv = OmegaInvBlock, snp_idx = snpIdx)
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
    idx <- OmegaInv[[b]]$snp_idx
    Oi <- OmegaInv[[b]]$Omega_inv
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
  nSnps <- nrow(ldRef@snpInfo)
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
#' @return A list with se and loo_estimates.
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
    loo_estimates = looEstimates
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
#' @return A data.frame with block_id, h2_local, h2_local_se.
#' @keywords internal
.gldscLocal <- function(z, n, ldRef, h2, intercept) {
  nLdBlocks <- length(ldRef@ldBlocks@blocks)
  if (nLdBlocks <= 22) {
    stop("Local g-LDSC requires fine-grained LD blocks (e.g., Berisa & Pickrell loci), ",
         "but the LD reference only has ", nLdBlocks, " blocks (likely per-chromosome). ",
         "Provide explicit ldBlocks when reading the LD reference with readLdRef().")
  }
  chi2 <- z^2
  M <- length(z)
  blockIndices <- snpsPerBlock(ldRef@snpInfo, ldRef@ldBlocks)
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
    ldBlock <- ldRef@ldScores[idx, 1]

    sumLd <- sum(ldBlock)
    h2_b <- (sum(chi2Block) - MBlock * (n * intercept + 1)) /
            (n * sumLd / M)

    # SE from block-level Fisher information:
    # Var(h2_local) ~ M_block / (n * sum_ld / M)^2
    h2LocalSeB <- sqrt(2 * MBlock) / (n * sumLd / M)

    blockId[b] <- b
    h2Local[b] <- h2_b
    h2LocalSe[b] <- h2LocalSeB
  }

  data.frame(
    block_id = blockId,
    h2_local = h2Local,
    h2_local_se = h2LocalSe,
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
#' @return A list with enrichment (data.frame) and score_stats (list with
#'   z, R, annotation_names).
#' @keywords internal
.gldscScoreStats <- function(z, n, ldRef, annotations, coef, w) {
  # Score statistics for candidate annotations
  candidate <- getCandidates(annotations)
  nCand <- ncol(candidate@annotations)
  if (nCand == 0) return(list(enrichment = NULL, score_stats = NULL))

  # Compute score for each candidate annotation
  # S_a = dLL/d(tau_a) at tau_a = 0, baseline fitted
  chi2 <- z^2
  M <- length(z)
  fitted <- chi2 - weightedLs(chi2,
    cbind(n * ldRef@ldScores[, 1] / M, rep(n, M)), w)$residuals
  resid <- chi2 - fitted

  # Precompute weighted annotation LD scores
  annotLdMat <- matrix(0, nrow = M, ncol = nCand)
  for (a in seq_len(nCand)) {
    annotLdMat[, a] <- ldRef@ldScores[, 1] * candidate@annotations[, a]
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

  list(
    enrichment = data.frame(
      annotation = candidate@annotationMeta$name,
      score_z = scoreZ,
      score_p = 2 * pnorm(-abs(scoreZ)),
      stringsAsFactors = FALSE
    ),
    score_stats = list(
      z = scoreZ, R = R,
      annotation_names = candidate@annotationMeta$name
    )
  )
}
