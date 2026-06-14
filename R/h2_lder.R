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
  nBlocks <- length(eigenRef@eigenList)
  nRef <- eigenRef@nRef
  inSample <- eigenRef@inSample
  M <- nrow(eigenRef@snpInfo)

  # Extract baseline annotations if provided
  baselineMat <- NULL
  if (!is.null(annotations)) {
    baseline <- getBaseline(annotations)
    if (ncol(baseline@annotations) > 0) {
      baselineMat <- baseline@annotations
    }
  }

  # Collect per-block eigenvalue regression quantities
  blockData <- lapply(seq_len(nBlocks), function(b) {
    block <- eigenRef@eigenList[[b]]
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
  nCandidates <- ncol(candidateAnnot@annotations)

  if (nCandidates == 0) {
    return(list(enrichment = NULL, scoreStats = NULL))
  }

  nBlocks <- length(eigenRef@eigenList)
  M <- nrow(eigenRef@snpInfo)

  # Collect per-block partial scores into a matrix (nBlocks x nCandidates)
  partialsMat <- matrix(0, nrow = nBlocks, ncol = nCandidates)

  for (b in seq_len(nBlocks)) {
    block <- eigenRef@eigenList[[b]]
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
      annotCol <- candidateAnnot@annotations[, ai]
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

  enrichmentDf <- data.frame(
    annotation = candidateAnnot@annotationMeta$name,
    scoreZ = scoreZ,
    scoreP = 2 * pnorm(-abs(scoreZ)),
    stringsAsFactors = FALSE
  )

  scoreStatsList <- list(
    z = scoreZ,
    R = R,
    annotationNames = candidateAnnot@annotationMeta$name
  )

  list(enrichment = enrichmentDf, scoreStats = scoreStatsList)
}
