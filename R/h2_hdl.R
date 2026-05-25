#' @title HDL/sHDL: High-Definition Likelihood
#' @description Estimate heritability using eigenvalue-based likelihood
#'   maximization (Ning et al. 2020) with stratified extensions (sHDL,
#'   Kim et al. 2023).
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
#' @param eigen_ref An \code{LDEigen} object.
#' @param annotations An \code{AnnotationMatrix}, or NULL.
#' @param local Logical, return per-block estimates.
#' @param lambda Numeric, L2 penalty on tau (default 0).
#' @return List with h2, h2_se, local, enrichment, score_stats.
#' @keywords internal
hdl_univariate <- function(z, n, eigen_ref, annotations = NULL,
                           local = FALSE, lambda = 0) {
  n_blocks <- length(eigen_ref@eigen_list)
  M <- nrow(eigen_ref@snp_info)

  # Extract baseline annotations if provided
  baseline_mat <- NULL
  if (!is.null(annotations)) {
    baseline <- getBaseline(annotations)
    if (ncol(baseline@annotations) > 0) {
      baseline_mat <- baseline@annotations
    }
  }

  # Precompute per-block quantities including annotation-stratified
  # eigenvalue scores
  block_data <- lapply(seq_len(n_blocks), function(b) {
    block <- eigen_ref@eigen_list[[b]]
    idx <- block$snp_idx
    V <- block$vectors
    d <- block$values
    z_rot <- as.vector(t(V) %*% z[idx])

    if (!is.null(baseline_mat)) {
      ld_annot <- crossprod(V^2, baseline_mat[idx, , drop = FALSE])
    } else {
      ld_annot <- NULL
    }

    list(z_rot = z_rot, d = d, ld_annot = ld_annot,
         snp_idx = idx, p = length(idx))
  })

  # Compute per-eigenvalue variance: sigma2_i = n/M * sum_a(tau_a * d_i * ld_annot_{a,i}) + 1
  .compute_sigma2 <- function(tau, bd) {
    if (!is.null(bd$ld_annot)) {
      n / M * bd$d * as.vector(bd$ld_annot %*% tau) + 1
    } else {
      n * tau[1] * bd$d / M + 1
    }
  }

  if (!is.null(baseline_mat)) {
    n_tau <- ncol(baseline_mat)

    # Negative log-likelihood as function of tau vector (with optional L2 penalty)
    nll <- function(tau) {
      val <- 0
      for (bd in block_data) {
        sigma2 <- .compute_sigma2(tau, bd)
        if (any(sigma2 <= 0)) return(1e10)
        val <- val + 0.5 * sum(log(sigma2) + bd$z_rot^2 / sigma2)
      }
      if (lambda > 0) val <- val + lambda * sum(tau^2)
      val
    }

    # Gradient of negative log-likelihood
    nll_grad <- function(tau) {
      grad <- numeric(n_tau)
      for (bd in block_data) {
        sigma2 <- .compute_sigma2(tau, bd)
        if (any(sigma2 <= 0)) return(rep(0, n_tau))
        # dsigma2/dtau_a = n/M * d_i * ld_annot_{a,i}
        dsig <- n / M * bd$d * bd$ld_annot  # (n_eigen x n_tau)
        # dNLL/dsigma2_i = 0.5 * (1/sigma2_i - z_rot_i^2/sigma2_i^2)
        dNLL_dsig <- 0.5 * (1 / sigma2 - bd$z_rot^2 / sigma2^2)
        grad <- grad + as.vector(crossprod(dsig, dNLL_dsig))
      }
      if (lambda > 0) grad <- grad + 2 * lambda * tau
      grad
    }

    # Initialize with uniform h2 across annotations
    tau_init <- rep(0.5 / n_tau, n_tau)

    opt <- optim(tau_init, nll, gr = nll_grad, method = "BFGS",
                        control = list(maxit = 200, reltol = 1e-8))
    tau <- opt$par

    # h2 = sum_a tau_a * M_a
    h2 <- sum(tau * colSums(baseline_mat))
  } else {
    # Single-parameter case: use optimize (current behavior)
    nll <- function(h2) {
      val <- 0
      for (bd in block_data) {
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
  fisher_result <- .hdl_se_fisher_stratified(tau, block_data, n, M,
                                             baseline_mat, lambda = lambda)
  h2_se <- fisher_result$h2_se
  tau_se <- fisher_result$tau_se

  # Jackknife tau_blocks via score approximation
  tau_blocks <- NULL
  if (!is.null(baseline_mat)) {
    jk <- .hdl_jackknife_tau(tau, block_data, n, M, baseline_mat,
                             lambda = lambda)
    tau_blocks <- jk$loo_estimates
    # Use jackknife SE if available (more robust than Fisher for enrichment)
    tau_se <- jk$tau_se
  }

  # Baseline enrichment
  baseline_enrichment_df <- NULL
  if (!is.null(baseline_mat)) {
    annot_names <- if (!is.null(colnames(baseline_mat))) {
      colnames(baseline_mat)
    } else {
      paste0("annot_", seq_len(ncol(baseline_mat)))
    }
    baseline_enrichment_df <- computeBaselineEnrichment(
      tau, tau_se, tau_blocks, baseline_mat, annot_names, h2
    )
  }

  # Local estimation
  local_df <- NULL
  if (local) {
    local_df <- .hdl_local(block_data, n, M, tau, baseline_mat)
  }

  # Score statistics for candidate annotations (sHDL)
  score_stats <- NULL
  if (!is.null(annotations)) {
    strat <- .shdl_stratified(z, n, eigen_ref, annotations, tau,
                              baseline_mat)
    score_stats <- strat$score_stats
  }

  list(h2 = h2, h2_se = h2_se, intercept = NA_real_,
       intercept_se = NA_real_, tau = tau, tau_se = tau_se,
       tau_blocks = tau_blocks, local = local_df,
       enrichment = baseline_enrichment_df, score_stats = score_stats)
}

# =============================================================================
# Internal helpers
# =============================================================================

#' @keywords internal
.hdl_se_fisher <- function(h2, block_data, n, M) {
  # Fisher information = -E[d^2 LL / d h2^2]
  info <- 0
  for (bd in block_data) {
    sigma2 <- n * h2 * bd$d / M + 1
    info <- info + 0.5 * sum((n * bd$d / M)^2 / sigma2^2)
  }
  h2_se <- 1 / sqrt(max(info, 1e-10))
  list(h2_se = h2_se, tau_se = h2_se, tau_vcov = NULL)
}

#' @keywords internal
.hdl_se_fisher_stratified <- function(tau, block_data, n, M,
                                      baseline_mat = NULL, lambda = 0) {
  if (is.null(baseline_mat)) {
    return(.hdl_se_fisher(tau, block_data, n, M))
  }

  # Fisher information matrix for tau, then delta method for h2
  n_tau <- length(tau)
  info_mat <- matrix(0, nrow = n_tau, ncol = n_tau)
  for (bd in block_data) {
    sigma2 <- n / M * bd$d * as.vector(bd$ld_annot %*% tau) + 1
    # dsigma2/dtau_a = n/M * d_i * ld_annot_{a,i}
    dsig <- n / M * bd$d * bd$ld_annot  # (n_eigen x n_tau)
    # Fisher info contribution: sum_i dsig_a * dsig_b / sigma2_i^2
    w <- 1 / sigma2^2
    info_mat <- info_mat + 0.5 * crossprod(dsig * sqrt(w))
  }
  # Add ridge penalty contribution to Hessian
  if (lambda > 0) {
    info_mat <- info_mat + 2 * lambda * diag(n_tau)
  }

  # Variance of tau
  tau_vcov <- tryCatch(solve(info_mat), error = function(e) {
    diag(1 / pmax(diag(info_mat), 1e-10))
  })
  tau_se <- sqrt(pmax(diag(tau_vcov), 0))

  # h2 = sum_a tau_a * M_a, so dh2/dtau = M_a
  M_a <- colSums(baseline_mat)
  h2_var <- as.numeric(t(M_a) %*% tau_vcov %*% M_a)
  h2_se <- sqrt(max(h2_var, 0))

  list(h2_se = h2_se, tau_se = tau_se, tau_vcov = tau_vcov)
}

#' @title Block-level jackknife for HDL tau via score approximation
#' @description Approximate leave-one-block-out tau estimates using the
#'   linear approximation: tau_loo_b ~ tau - H^{-1} * grad_b, where H is
#'   the Fisher information (Hessian of NLL) and grad_b is block b's
#'   gradient contribution at the MLE.
#' @keywords internal
.hdl_jackknife_tau <- function(tau, block_data, n, M, baseline_mat,
                               lambda = 0) {
  n_blocks <- length(block_data)
  n_tau <- length(tau)

  info_mat <- matrix(0, n_tau, n_tau)
  block_grads <- matrix(0, n_blocks, n_tau)

  for (b in seq_len(n_blocks)) {
    bd <- block_data[[b]]
    sigma2 <- n / M * bd$d * as.vector(bd$ld_annot %*% tau) + 1
    dsig <- n / M * bd$d * bd$ld_annot
    # Block gradient contribution
    dNLL_dsig <- 0.5 * (1 / sigma2 - bd$z_rot^2 / sigma2^2)
    block_grads[b, ] <- as.vector(crossprod(dsig, dNLL_dsig))
    # Block Fisher info contribution
    w <- 1 / sigma2^2
    info_mat <- info_mat + 0.5 * crossprod(dsig * sqrt(w))
  }
  # Add ridge penalty to Hessian
  if (lambda > 0) {
    info_mat <- info_mat + 2 * lambda * diag(n_tau)
  }

  H_inv <- tryCatch(solve(info_mat), error = function(e) {
    diag(1 / pmax(diag(info_mat), 1e-10))
  })

  # LOO estimates via linear approximation at MLE
  loo_estimates <- matrix(NA, n_blocks, n_tau)
  for (b in seq_len(n_blocks)) {
    loo_estimates[b, ] <- tau - as.vector(H_inv %*% block_grads[b, ])
  }

  list(
    loo_estimates = loo_estimates,
    tau_se = jackknifeSe(tau, loo_estimates)
  )
}

#' @keywords internal
.hdl_local <- function(block_data, n, M, tau = NULL,
                       baseline_mat = NULL) {
  local_results <- lapply(seq_along(block_data), function(b) {
    bd <- block_data[[b]]
    if (bd$p < 3) {
      return(data.frame(block_id = b, h2_local = NA, h2_local_se = NA))
    }

    # Compute the global baseline sigma2 for this block
    if (!is.null(baseline_mat) && !is.null(bd$ld_annot)) {
      sigma2_baseline <- n / M * bd$d *
        as.vector(bd$ld_annot %*% tau) + 1
    } else if (!is.null(tau)) {
      sigma2_baseline <- n * tau[1] * bd$d / M + 1
    } else {
      sigma2_baseline <- NULL
    }

    # Local likelihood: optimize a local deviation delta_h2
    # sigma2_i = sigma2_baseline_i + n * delta_h2 * d_i / M
    nll_local <- function(delta_h2) {
      if (!is.null(sigma2_baseline)) {
        sigma2 <- sigma2_baseline + n * delta_h2 * bd$d / M
      } else {
        sigma2 <- n * delta_h2 * bd$d / M + 1
      }
      if (any(sigma2 <= 0)) return(1e10)
      0.5 * sum(log(sigma2) + bd$z_rot^2 / sigma2)
    }
    opt <- optimize(nll_local, interval = c(-0.5, 0.5), tol = 1e-8)
    delta_h2 <- opt$minimum

    # Total local h2 = global baseline contribution + local deviation
    if (!is.null(sigma2_baseline)) {
      # The baseline already explains some h2 in this block; delta is
      # the additional local deviation
      h2_local <- delta_h2
    } else {
      h2_local <- delta_h2
    }

    # SE from Fisher information at optimum
    sigma2_opt <- if (!is.null(sigma2_baseline)) {
      sigma2_baseline + n * delta_h2 * bd$d / M
    } else {
      n * delta_h2 * bd$d / M + 1
    }
    info <- 0.5 * sum((n * bd$d / M)^2 / sigma2_opt^2)
    se_local <- 1 / sqrt(max(info, 1e-10))

    data.frame(block_id = b, h2_local = h2_local, h2_local_se = se_local)
  })
  do.call(rbind, local_results)
}

#' @keywords internal
.shdl_stratified <- function(z, n, eigen_ref, annotations, tau,
                             baseline_mat = NULL) {
  # sHDL: stratified HDL with annotation-specific h2
  # Score-based approach: baseline fit uses jointly fitted tau
  candidate <- getCandidates(annotations)
  n_cand <- ncol(candidate@annotations)
  if (n_cand == 0) return(list(enrichment = NULL, score_stats = NULL))

  M <- nrow(eigen_ref@snp_info)
  n_blocks <- length(eigen_ref@eigen_list)

  # Score for each candidate: derivative of log-likelihood w.r.t. tau_a
  # at tau_a = 0, with baseline tau at the MLE
  # Collect per-block gradient and information for jackknife
  block_grad <- matrix(0, nrow = n_blocks, ncol = n_cand)
  block_info <- matrix(0, nrow = n_blocks, ncol = n_cand)

  for (a in seq_len(n_cand)) {
    annot_col <- candidate@annotations[, a]

    for (b in seq_len(n_blocks)) {
      block <- eigen_ref@eigen_list[[b]]
      idx <- block$snp_idx
      V <- block$vectors
      d <- block$values
      z_rot <- as.vector(t(V) %*% z[idx])

      # Annotation-stratified eigenvalue score for this candidate
      annot_block <- annot_col[idx]
      d_annot <- as.vector(t(V^2) %*% annot_block)

      # Compute sigma2 from stratified baseline fit
      if (!is.null(baseline_mat)) {
        ld_annot_base <- crossprod(V^2, baseline_mat[idx, , drop = FALSE])
        sigma2 <- n / M * d * as.vector(ld_annot_base %*% tau) + 1
      } else {
        sigma2 <- n * tau[1] * d / M + 1
      }

      # Gradient: sum_i [ (z_rot_i^2 / sigma2_i - 1) * n * d_annot_i / M ] /
      #           (2 * sigma2_i)
      grad <- 0.5 * sum((z_rot^2 / sigma2 - 1) * n * d_annot /
                          (M * sigma2))
      info <- 0.5 * sum((n * d_annot / M)^2 / sigma2^2)

      block_grad[b, a] <- grad
      block_info[b, a] <- info
    }
  }

  # Total score statistics
  total_grad <- colSums(block_grad)
  total_info <- colSums(block_info)
  score_z <- total_grad / sqrt(pmax(total_info, 1e-10))

  # Jackknife score correlation R: LOO by block
  loo_score_z <- matrix(0, nrow = n_blocks, ncol = n_cand)
  for (b in seq_len(n_blocks)) {
    loo_grad <- total_grad - block_grad[b, ]
    loo_info <- total_info - block_info[b, ]
    loo_score_z[b, ] <- loo_grad / sqrt(pmax(loo_info, 1e-10))
  }
  if (n_cand > 1 && n_blocks > 2) {
    R <- cor(loo_score_z)
  } else {
    R <- diag(n_cand)
  }

  list(
    enrichment = data.frame(
      annotation = candidate@annotation_meta$name,
      score_z = score_z,
      score_p = 2 * pnorm(-abs(score_z)),
      stringsAsFactors = FALSE
    ),
    score_stats = list(z = score_z, R = R,
                       annotation_names = candidate@annotation_meta$name)
  )
}
