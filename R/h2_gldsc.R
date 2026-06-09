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
#' @param ld_ref An \code{LDScore} object.
#' @param annotations An \code{AnnotationMatrix}, or NULL.
#' @param local Logical, return per-block estimates.
#' @param lambda Numeric, ridge penalty (default 0).
#' @return A list with h2, h2_se, intercept, enrichment, score_stats.
#' @keywords internal
gldsc_univariate <- function(z, n, ld_ref, annotations = NULL,
                             local = FALSE, lambda = 0) {
  chi2 <- z^2
  M <- length(z)
  ld_scores <- ld_ref@ld_scores
  weights <- ld_ref@ld_score_weights

  # --- Step 1: Initial S-LDSC estimate (WLS) ---
  # Model: E[chi2_j] = N/M * sum_a(tau_a * l_{j,a}) + N*a + 1
  if (!is.null(annotations)) {
    baseline <- getBaseline(annotations)
    ld_strat <- computeLdScores(ld_ref, baseline)
    X <- cbind(n * ld_strat / M, rep(n, M))
  } else {
    X <- cbind(n * ld_scores[, 1] / M, rep(n, M))
  }
  y <- chi2

  # Initial WLS (standard S-LDSC)
  fit_wls <- weightedLsRidge(y, X, weights, lambda = lambda,
                             penalize_intercept = FALSE)
  fitted_wls <- fit_wls$fitted

  # --- Step 2: Estimate residual covariance (FGLS) ---
  # The residual covariance arises from LD between SNPs
  # Cov(resid_j, resid_k) depends on R^2_{jk}
  # Approximate the residual covariance using the LD matrix list
  if (length(ld_ref@ld_matrix_list) > 0) {
    # Use stored LD matrices for FGLS
    Omega_inv <- .compute_fgls_weights(ld_ref@ld_matrix_list, fitted_wls)
  } else {
    # Approximate: use LD scores as proxy for diagonal residual variance
    # This gives an intermediate estimator between S-LDSC and full g-LDSC
    resid_var <- 2 * pmax(fitted_wls, 1)^2
    Omega_inv <- 1 / resid_var
    message("Note: g-LDSC without full LD matrices uses approximate FGLS. ",
            "For full g-LDSC, compute LD reference with ld_matrix_list.")
  }

  # --- Step 3: FGLS estimate ---
  if (is.numeric(Omega_inv)) {
    # Diagonal approximation
    fit_fgls <- weightedLsRidge(y, X, Omega_inv, lambda = lambda,
                                penalize_intercept = FALSE)
  } else {
    # Full FGLS with block-diagonal Omega_inv
    fit_fgls <- .fgls_solve(y, X, Omega_inv)
  }

  n_params <- ncol(X)
  n_tau <- n_params - 1L
  tau_full <- fit_fgls$coef[seq_len(n_tau)]
  intercept <- fit_fgls$coef[n_params]

  # Compute h2 from tau: h2 = sum_a(tau_a * M_a)
  # When annotations are present, computeLdScores prepends a total L2 column
  # (base = all-ones annotation), so tau_full[1] = tau_base with M_base = M,
  # and tau_full[2:n_tau] = per-annotation tau.
  if (!is.null(annotations)) {
    baseline_mat <- getBaseline(annotations)@annotations
    M_a <- colSums(baseline_mat)
    M_a_full <- c(M, M_a)  # base (all-ones) + annotation-specific
    h2 <- sum(tau_full * M_a_full)
    # For enrichment, use only the annotation-specific coefficients
    tau <- tau_full[-1]
  } else {
    baseline_mat <- NULL
    M_a_full <- NULL
    h2 <- tau_full[1]  # unstratified: coefficient IS h2
    tau <- tau_full
  }

  # Jackknife SE and LOO estimates
  block_assign <- .assign_snps_to_jackknife_blocks(ld_ref, n_blocks = 200)
  jk <- .gldsc_jackknife(y, X, Omega_inv, fit_fgls$coef, block_assign,
                          lambda = lambda)

  # Extract per-annotation tau jackknife blocks and SE
  tau_blocks_full <- jk$loo_estimates[, seq_len(n_tau), drop = FALSE]

  # h2 and intercept SE
  if (!is.null(baseline_mat)) {
    h2_loo <- as.vector(tau_blocks_full %*% M_a_full)
    # Annotation-specific blocks (exclude base L2 column)
    tau_blocks <- tau_blocks_full[, -1, drop = FALSE]
  } else {
    h2_loo <- jk$loo_estimates[, 1]
    tau_blocks <- tau_blocks_full
  }
  tau_se <- jackknifeSe(tau, tau_blocks)
  a_loo <- jk$loo_estimates[, n_params]
  se <- jackknifeSe(c(h2, intercept), cbind(h2_loo, a_loo))

  # Baseline enrichment (annotation-specific only, not base L2)
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

  # Local estimates
  local_df <- NULL
  if (local) {
    local_df <- .gldsc_local(z, n, ld_ref, h2, intercept)
  }

  # Score statistics for candidate annotations
  score_stats <- NULL
  if (!is.null(annotations)) {
    strat <- .gldsc_score_stats(z, n, ld_ref, annotations,
                                fit_fgls$coef, weights)
    score_stats <- strat$score_stats
  }

  list(
    h2 = h2,
    h2_se = se[1],
    intercept = intercept,
    intercept_se = se[2],
    tau = tau,
    tau_se = tau_se,
    tau_blocks = tau_blocks,
    local = local_df,
    enrichment = baseline_enrichment_df,
    score_stats = score_stats
  )
}

# =============================================================================
# Internal helpers
# =============================================================================

#' Compute block-diagonal residual covariance for FGLS
#'
#' @description Compute per-block precision matrices from LD matrices
#'   and fitted values for feasible GLS estimation.
#' @param ld_matrix_list List of LD matrix blocks, each with R and snp_idx.
#' @param fitted_values Numeric vector of fitted values from initial WLS.
#' @return A list of per-block precision matrices (inverse Omega blocks).
#' @keywords internal
.compute_fgls_weights <- function(ld_matrix_list, fitted_values) {
  # Omega_{jk} = 2 * (fitted_j * fitted_k * R^2_{jk} + R^4_{jk})
  result <- vector("list", length(ld_matrix_list))
  for (b in seq_along(ld_matrix_list)) {
    block <- ld_matrix_list[[b]]
    R_block <- block$R
    snp_idx <- block$snp_idx
    fitted_block <- fitted_values[snp_idx]
    m_block <- length(snp_idx)

    R2 <- R_block^2
    R4 <- R2^2
    # Omega_block[j,k] = 2 * (fitted_j * fitted_k * R^2_jk + R^4_jk)
    fitted_outer <- outer(fitted_block, fitted_block)
    Omega_block <- 2 * (fitted_outer * R2 + R4)

    # Regularise to ensure positive definiteness
    Omega_block <- Omega_block + diag(1e-6, m_block)
    Omega_inv_block <- solve(Omega_block)

    result[[b]] <- list(Omega_inv = Omega_inv_block, snp_idx = snp_idx)
  }
  result
}

#' Solve GLS with block-diagonal precision
#'
#' @description Solve the GLS problem beta = (X' Omega^{-1} X)^{-1} X' Omega^{-1} y
#'   for block-diagonal Omega_inv.
#' @param y Numeric vector, response.
#' @param X Numeric matrix, design matrix.
#' @param Omega_inv Block-diagonal precision: either a numeric vector (diagonal)
#'   or a list of per-block precision matrices.
#' @return A list with coef, fitted, and residuals.
#' @keywords internal
.fgls_solve <- function(y, X, Omega_inv) {
  if (is.null(dim(X))) X <- matrix(X, ncol = 1)
  p <- ncol(X)

  if (is.numeric(Omega_inv)) {
    # Diagonal case: standard WLS
    return(weightedLs(y, X, Omega_inv))
  }

  # Block-diagonal case: accumulate across blocks
  XtOiX <- matrix(0, nrow = p, ncol = p)
  XtOiy <- numeric(p)

  for (b in seq_along(Omega_inv)) {
    idx <- Omega_inv[[b]]$snp_idx
    Oi <- Omega_inv[[b]]$Omega_inv
    X_b <- X[idx, , drop = FALSE]
    y_b <- y[idx]

    XtOiX <- XtOiX + crossprod(X_b, Oi %*% X_b)
    XtOiy <- XtOiy + crossprod(X_b, Oi %*% y_b)
  }

  beta <- as.vector(solve(XtOiX, XtOiy))
  fitted_vals <- as.vector(X %*% beta)
  resid <- y - fitted_vals

  list(coef = beta, fitted = fitted_vals, residuals = resid)
}

#' Assign SNPs to jackknife blocks
#'
#' @description Divide SNPs into approximately equal-sized jackknife blocks.
#' @param ld_ref An \code{LDScore} object.
#' @param n_blocks Integer, number of jackknife blocks (default 200).
#' @return Integer vector of block assignments.
#' @keywords internal
.assign_snps_to_jackknife_blocks <- function(ld_ref, n_blocks = 200) {
  n_snps <- nrow(ld_ref@snp_info)
  block_size <- ceiling(n_snps / n_blocks)
  rep(seq_len(n_blocks), each = block_size, length.out = n_snps)
}

#' Jackknife SE for g-LDSC
#'
#' @description Compute leave-one-block-out jackknife standard errors
#'   for g-LDSC coefficient estimates.
#' @param y Numeric vector, response.
#' @param X Numeric matrix, design matrix.
#' @param w Weights: numeric vector (diagonal) or list of precision blocks.
#' @param coef_full Numeric vector, full-sample coefficients.
#' @param block_assign Integer vector, block assignments from
#'   \code{.assign_snps_to_jackknife_blocks}.
#' @param lambda Numeric, ridge penalty (default 0).
#' @return A list with se and loo_estimates.
#' @keywords internal
.gldsc_jackknife <- function(y, X, w, coef_full, block_assign,
                             lambda = 0) {
  n_blocks <- max(block_assign)
  n_params <- length(coef_full)
  loo_estimates <- matrix(NA, nrow = n_blocks, ncol = n_params)

  for (b in seq_len(n_blocks)) {
    keep <- block_assign != b
    if (is.numeric(w)) {
      fit_loo <- weightedLsRidge(y[keep], X[keep, , drop = FALSE], w[keep],
                                 lambda = lambda, penalize_intercept = FALSE)
    } else {
      fit_loo <- weightedLsRidge(y[keep], X[keep, , drop = FALSE],
                                 rep(1, sum(keep)),
                                 lambda = lambda, penalize_intercept = FALSE)
    }
    loo_estimates[b, ] <- fit_loo$coef
  }
  list(
    se = jackknifeSe(coef_full, loo_estimates),
    loo_estimates = loo_estimates
  )
}

#' Per-block local h2 from g-LDSC
#'
#' @description Estimate per-block local heritability using the
#'   g-LDSC intercept and LD scores.
#' @param z Numeric vector of z-scores.
#' @param n Numeric, GWAS sample size.
#' @param ld_ref An \code{LDScore} object.
#' @param h2 Numeric, global h2 estimate.
#' @param intercept Numeric, g-LDSC intercept estimate.
#' @return A data.frame with block_id, h2_local, h2_local_se.
#' @keywords internal
.gldsc_local <- function(z, n, ld_ref, h2, intercept) {
  n_ld_blocks <- length(ld_ref@ld_blocks@blocks)
  if (n_ld_blocks <= 22) {
    stop("Local g-LDSC requires fine-grained LD blocks (e.g., Berisa & Pickrell loci), ",
         "but the LD reference only has ", n_ld_blocks, " blocks (likely per-chromosome). ",
         "Provide explicit ld_blocks when reading the LD reference with readLdRef().")
  }
  chi2 <- z^2
  M <- length(z)
  block_indices <- snpsPerBlock(ld_ref@snp_info, ld_ref@ld_blocks)
  n_blocks <- length(block_indices)

  block_id <- integer(n_blocks)
  h2_local <- numeric(n_blocks)
  h2_local_se <- numeric(n_blocks)

  for (b in seq_len(n_blocks)) {
    idx <- block_indices[[b]]
    M_block <- length(idx)
    if (M_block == 0) {
      block_id[b] <- b
      h2_local[b] <- NA_real_
      h2_local_se[b] <- NA_real_
      next
    }
    chi2_block <- chi2[idx]
    ld_block <- ld_ref@ld_scores[idx, 1]

    sum_ld <- sum(ld_block)
    h2_b <- (sum(chi2_block) - M_block * (n * intercept + 1)) /
            (n * sum_ld / M)

    # SE from block-level Fisher information:
    # Var(h2_local) ~ M_block / (n * sum_ld / M)^2
    h2_local_se_b <- sqrt(2 * M_block) / (n * sum_ld / M)

    block_id[b] <- b
    h2_local[b] <- h2_b
    h2_local_se[b] <- h2_local_se_b
  }

  data.frame(
    block_id = block_id,
    h2_local = h2_local,
    h2_local_se = h2_local_se,
    stringsAsFactors = FALSE
  )
}

#' Score statistics for candidate annotations
#'
#' @description Compute score z-statistics and their correlation matrix
#'   for candidate annotations not included in the baseline model.
#' @param z Numeric vector of z-scores.
#' @param n Numeric, GWAS sample size.
#' @param ld_ref An \code{LDScore} object.
#' @param annotations An \code{AnnotationMatrix} object.
#' @param coef Numeric vector, fitted coefficients from the baseline model.
#' @param w Numeric vector, regression weights.
#' @return A list with enrichment (data.frame) and score_stats (list with
#'   z, R, annotation_names).
#' @keywords internal
.gldsc_score_stats <- function(z, n, ld_ref, annotations, coef, w) {
  # Score statistics for candidate annotations
  candidate <- getCandidates(annotations)
  n_cand <- ncol(candidate@annotations)
  if (n_cand == 0) return(list(enrichment = NULL, score_stats = NULL))

  # Compute score for each candidate annotation
  # S_a = dLL/d(tau_a) at tau_a = 0, baseline fitted
  chi2 <- z^2
  M <- length(z)
  fitted <- chi2 - weightedLs(chi2,
    cbind(n * ld_ref@ld_scores[, 1] / M, rep(n, M)), w)$residuals
  resid <- chi2 - fitted

  # Precompute weighted annotation LD scores
  annot_ld_mat <- matrix(0, nrow = M, ncol = n_cand)
  for (a in seq_len(n_cand)) {
    annot_ld_mat[, a] <- ld_ref@ld_scores[, 1] * candidate@annotations[, a]
  }

  # Full-sample score z-statistics
  score_z <- numeric(n_cand)
  for (a in seq_len(n_cand)) {
    S_a <- sum(w * resid * n * annot_ld_mat[, a] / M)
    V_a <- sum((w * n * annot_ld_mat[, a] / M)^2)
    score_z[a] <- S_a / sqrt(V_a)
  }

  # Block jackknife for score correlation matrix R
  n_blocks <- 200
  block_assign <- .assign_snps_to_jackknife_blocks(ld_ref, n_blocks = n_blocks)
  n_blocks_actual <- max(block_assign)
  score_z_loo <- matrix(NA, nrow = n_blocks_actual, ncol = n_cand)

  for (b in seq_len(n_blocks_actual)) {
    keep <- block_assign != b
    w_b <- w[keep]
    resid_b <- resid[keep]
    for (a in seq_len(n_cand)) {
      ald_b <- annot_ld_mat[keep, a]
      S_a <- sum(w_b * resid_b * n * ald_b / M)
      V_a <- sum((w_b * n * ald_b / M)^2)
      score_z_loo[b, a] <- S_a / sqrt(V_a)
    }
  }

  # Score correlation from jackknife LOO z-statistics
  if (n_cand > 1) {
    R <- cor(score_z_loo)
    # Ensure valid correlation matrix
    R[is.na(R)] <- 0
  } else {
    R <- matrix(1, 1, 1)
  }

  list(
    enrichment = data.frame(
      annotation = candidate@annotation_meta$name,
      score_z = score_z,
      score_p = 2 * pnorm(-abs(score_z)),
      stringsAsFactors = FALSE
    ),
    score_stats = list(
      z = score_z, R = R,
      annotation_names = candidate@annotation_meta$name
    )
  )
}
