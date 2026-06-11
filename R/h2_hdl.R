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
#' @return List with h2, h2_se, local, enrichment, score_stats.
#' @keywords internal
hdlUnivariate <- function(z, n, eigenRef, annotations = NULL,
                          local = FALSE, lambda = 0) {
  nBlocks <- length(eigenRef@eigenList)
  M <- nrow(eigenRef@snpInfo)

  # Extract baseline annotations if provided
  baselineMat <- NULL
  if (!is.null(annotations)) {
    baseline <- getBaseline(annotations)
    if (ncol(baseline@annotations) > 0) {
      baselineMat <- baseline@annotations
    }
  }

  # Precompute per-block quantities including annotation-stratified
  # eigenvalue scores
  blockData <- lapply(seq_len(nBlocks), function(b) {
    block <- eigenRef@eigenList[[b]]
    idx <- block$snp_idx
    V <- block$vectors
    d <- block$values
    zRot <- as.vector(t(V) %*% z[idx])

    if (!is.null(baselineMat)) {
      ldAnnot <- crossprod(V^2, baselineMat[idx, , drop = FALSE])
    } else {
      ldAnnot <- NULL
    }

    list(z_rot = zRot, d = d, ld_annot = ldAnnot,
         snp_idx = idx, p = length(idx))
  })

  # Compute per-eigenvalue variance: sigma2_i = n/M * sum_a(tau_a * d_i * ld_annot_{a,i}) + 1
  .computeSigma2 <- function(tau, bd) {
    if (!is.null(bd$ld_annot)) {
      n / M * bd$d * as.vector(bd$ld_annot %*% tau) + 1
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
        val <- val + 0.5 * sum(log(sigma2) + bd$z_rot^2 / sigma2)
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
        dsig <- n / M * bd$d * bd$ld_annot  # (nEigen x nTau)
        # dNLL/dsigma2_i = 0.5 * (1/sigma2_i - z_rot_i^2/sigma2_i^2)
        dNLL_dsig <- 0.5 * (1 / sigma2 - bd$z_rot^2 / sigma2^2)
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
        val <- val + 0.5 * sum(log(sigma2) + bd$z_rot^2 / sigma2)
      }
      val
    }

    opt <- optimize(nll, interval = c(-0.5, 1.5), tol = 1e-8)
    h2 <- opt$minimum
    tau <- h2  # scalar, used by downstream functions
  }

  # SE from observed Fisher information (returns h2_se, tau_se, tau_vcov)
  fisherResult <- .hdlSeFisherStratified(tau, blockData, n, M,
                                         baselineMat, lambda = lambda)
  h2Se <- fisherResult$h2_se
  tauSe <- fisherResult$tau_se

  # Jackknife tau_blocks via score approximation
  tauBlocks <- NULL
  if (!is.null(baselineMat)) {
    jk <- .hdlJackknifeTau(tau, blockData, n, M, baselineMat,
                           lambda = lambda)
    tauBlocks <- jk$loo_estimates
    # Use jackknife SE if available (more robust than Fisher for enrichment)
    tauSe <- jk$tau_se
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
    scoreStats <- strat$score_stats
  }

  list(h2 = h2, h2_se = h2Se, intercept = NA_real_,
       intercept_se = NA_real_, tau = tau, tau_se = tauSe,
       tau_blocks = tauBlocks, local = localDf,
       enrichment = baselineEnrichmentDf, score_stats = scoreStats)
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
  list(h2_se = h2Se, tau_se = h2Se, tau_vcov = NULL)
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
    sigma2 <- n / M * bd$d * as.vector(bd$ld_annot %*% tau) + 1
    # dsigma2/dtau_a = n/M * d_i * ld_annot_{a,i}
    dsig <- n / M * bd$d * bd$ld_annot  # (nEigen x nTau)
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

  list(h2_se = h2Se, tau_se = tauSe, tau_vcov = tauVcov)
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
    sigma2 <- n / M * bd$d * as.vector(bd$ld_annot %*% tau) + 1
    dsig <- n / M * bd$d * bd$ld_annot
    # Block gradient contribution
    dNLL_dsig <- 0.5 * (1 / sigma2 - bd$z_rot^2 / sigma2^2)
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
    loo_estimates = looEstimates,
    tau_se = jackknifeSe(tau, looEstimates)
  )
}

#' @keywords internal
.hdlLocal <- function(blockData, n, M, tau = NULL,
                      baselineMat = NULL) {
  localResults <- lapply(seq_along(blockData), function(b) {
    bd <- blockData[[b]]
    if (bd$p < 3) {
      return(data.frame(block_id = b, h2_local = NA, h2_local_se = NA))
    }

    # Compute the global baseline sigma2 for this block
    if (!is.null(baselineMat) && !is.null(bd$ld_annot)) {
      sigma2Baseline <- n / M * bd$d *
        as.vector(bd$ld_annot %*% tau) + 1
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
      0.5 * sum(log(sigma2) + bd$z_rot^2 / sigma2)
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

    data.frame(block_id = b, h2_local = h2Local, h2_local_se = seLocal)
  })
  do.call(rbind, localResults)
}

#' @keywords internal
.shdlStratified <- function(z, n, eigenRef, annotations, tau,
                            baselineMat = NULL) {
  # sHDL: stratified HDL with annotation-specific h2
  # Score-based approach: baseline fit uses jointly fitted tau
  candidate <- getCandidates(annotations)
  nCand <- ncol(candidate@annotations)
  if (nCand == 0) return(list(enrichment = NULL, score_stats = NULL))

  M <- nrow(eigenRef@snpInfo)
  nBlocks <- length(eigenRef@eigenList)

  # Score for each candidate: derivative of log-likelihood w.r.t. tau_a
  # at tau_a = 0, with baseline tau at the MLE
  # Collect per-block gradient and information for jackknife
  blockGrad <- matrix(0, nrow = nBlocks, ncol = nCand)
  blockInfo <- matrix(0, nrow = nBlocks, ncol = nCand)

  for (a in seq_len(nCand)) {
    annotCol <- candidate@annotations[, a]

    for (b in seq_len(nBlocks)) {
      block <- eigenRef@eigenList[[b]]
      idx <- block$snp_idx
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

  list(
    enrichment = data.frame(
      annotation = candidate@annotationMeta$name,
      score_z = scoreZ,
      score_p = 2 * pnorm(-abs(scoreZ)),
      stringsAsFactors = FALSE
    ),
    score_stats = list(z = scoreZ, R = R,
                       annotation_names = candidate@annotationMeta$name)
  )
}
