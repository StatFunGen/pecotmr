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
#' @param eigen_ref An \code{LDEigen} object.
#' @param annotations An \code{AnnotationMatrix}, or NULL.
#' @param local Logical, return per-block estimates.
#' @param lambda Numeric, ridge penalty (default 0).
#' @return A list with h2, h2_se, intercept, intercept_se, local estimates,
#'   and enrichment estimates.
#' @keywords internal
lder_univariate <- function(z, n, eigen_ref, annotations = NULL,
                            local = FALSE, lambda = 0) {
  n_blocks <- length(eigen_ref@eigen_list)
  n_ref <- eigen_ref@n_ref
  in_sample <- eigen_ref@in_sample
  M <- nrow(eigen_ref@snp_info)

  # Extract baseline annotations if provided
  baseline_mat <- NULL
  if (!is.null(annotations)) {
    baseline <- getBaseline(annotations)
    if (ncol(baseline@annotations) > 0) {
      baseline_mat <- baseline@annotations
    }
  }

  # Collect per-block eigenvalue regression quantities
  block_data <- lapply(seq_len(n_blocks), function(b) {
    block <- eigen_ref@eigen_list[[b]]
    idx <- block$snp_idx
    d <- block$values        # eigenvalues
    V <- block$vectors       # eigenvectors
    z_block <- z[idx]

    # Rotate z-scores into eigenbasis
    z_rot <- as.vector(t(V) %*% z_block)
    chi2_rot <- z_rot^2

    # Annotation-stratified eigenvalue scores for baseline annotations
    # ld_annot[i, a] = sum_j V[j,i]^2 * annot[j, a]
    if (!is.null(baseline_mat)) {
      ld_annot <- crossprod(V^2, baseline_mat[idx, , drop = FALSE])
    } else {
      ld_annot <- NULL
    }

    list(
      chi2_rot = chi2_rot,
      eigenvalues = d,
      ld_annot = ld_annot,
      n_snps = length(idx),
      snp_idx = idx
    )
  })

  # Assemble regression data
  all_chi2 <- unlist(lapply(block_data, `[[`, "chi2_rot"))
  all_d <- unlist(lapply(block_data, `[[`, "eigenvalues"))

  # Build design matrix
  # Stratified model: E[chi2_rot_i - 1] = n/M * sum_a(tau_a * d_i * ld_annot_{a,i}) + n*a
  # Unstratified model (no baseline annotations): same with single base column
  if (!is.null(baseline_mat)) {
    all_ld_annot <- do.call(rbind, lapply(block_data, `[[`, "ld_annot"))
    X <- cbind(n * all_d * all_ld_annot / M, rep(n, length(all_d)))
    n_tau <- ncol(baseline_mat)
  } else {
    X <- cbind(n * all_d / M, rep(n, length(all_d)))
    n_tau <- 1L
  }

  y <- all_chi2 - 1
  w <- 1 / (2 * pmax(all_chi2, 1)^2)

  fit <- weightedLsRidge(y, X, w, lambda = lambda, penalize_intercept = FALSE)
  tau <- fit$coef[seq_len(n_tau)]
  a <- fit$coef[n_tau + 1]

  # Compute h2 from tau
  if (!is.null(baseline_mat)) {
    # h2 = sum_a tau_a * M_a where M_a = sum_j annot_{j,a}
    h2 <- sum(tau * colSums(baseline_mat))
  } else {
    h2 <- tau[1]
  }

  # Jackknife SE by block
  block_assign <- rep(seq_len(n_blocks),
    vapply(block_data, function(x) length(x$eigenvalues), integer(1)))

  loo_estimates <- matrix(NA, nrow = n_blocks, ncol = n_tau + 1)
  for (b in seq_len(n_blocks)) {
    keep <- block_assign != b
    fit_loo <- weightedLsRidge(y[keep], X[keep, , drop = FALSE], w[keep],
                               lambda = lambda, penalize_intercept = FALSE)
    loo_estimates[b, ] <- fit_loo$coef
  }

  # Extract per-annotation tau jackknife blocks and SE
  tau_blocks <- loo_estimates[, seq_len(n_tau), drop = FALSE]
  tau_se <- jackknifeSe(tau, tau_blocks)

  # Compute h2 for each LOO iteration, then jackknife
  if (!is.null(baseline_mat)) {
    M_a <- colSums(baseline_mat)
    h2_loo <- as.vector(tau_blocks %*% M_a)
  } else {
    h2_loo <- loo_estimates[, 1]
  }
  a_loo <- loo_estimates[, n_tau + 1]
  se <- jackknifeSe(c(h2, a), cbind(h2_loo, a_loo))

  # Baseline enrichment (if annotations provided)
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

  # Local heritability (if requested)
  local_df <- NULL
  if (local) {
    local_df <- .lder_local_h2(block_data, n, M, tau, a, baseline_mat)
  }

  # Score statistics for candidate annotations (if provided)
  score_stats <- NULL
  if (!is.null(annotations)) {
    strat <- .lder_stratified(z, n, eigen_ref, annotations, tau, a,
                              baseline_mat)
    score_stats <- strat$score_stats
  }

  list(
    h2 = h2,
    h2_se = se[1],
    intercept = a,
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

#' @title LDER local heritability
#' @description Per-block heritability using the Hessian-based SE.
#' @param block_data List of per-block eigenvalue regression quantities.
#' @param n Numeric, GWAS sample size.
#' @param M Integer, total number of SNPs.
#' @param tau Numeric vector of annotation coefficients.
#' @param a_global Numeric, global intercept.
#' @param baseline_mat Matrix of baseline annotations, or NULL.
#' @return A data.frame with block_id, h2_local, h2_local_se.
#' @keywords internal
.lder_local_h2 <- function(block_data, n, M, tau, a_global,
                           baseline_mat = NULL) {
  # Per-block heritability using the Hessian-based SE
  local_results <- lapply(seq_along(block_data), function(b) {
    bd <- block_data[[b]]
    p_block <- bd$n_snps
    d <- bd$eigenvalues
    chi2 <- bd$chi2_rot

    # Compute fitted baseline contribution for this block
    if (!is.null(baseline_mat)) {
      ld_annot <- bd$ld_annot  # n_eigenvalues x n_annotations
      fitted_baseline <- as.vector(n / M * d *
                                     (ld_annot %*% tau))
    } else {
      fitted_baseline <- n * tau[1] * d / M
    }

    # Local regression: residual after removing global baseline + intercept
    y <- chi2 - 1 - n * a_global - fitted_baseline
    x <- n * d / M
    if (length(y) < 3) {
      return(data.frame(block_id = b, h2_local = NA, h2_local_se = NA))
    }
    w <- 1 / (2 * pmax(chi2, 1)^2)
    h2_local <- sum(w * x * y) / sum(w * x^2)

    # Fisher information SE
    info <- sum(w * x^2)
    se_local <- 1 / sqrt(info)

    data.frame(block_id = b, h2_local = h2_local, h2_local_se = se_local)
  })
  do.call(rbind, local_results)
}

#' @title LDER stratified score statistics
#' @description Score-based approach: fit baseline jointly, compute scores
#'   for candidate annotations.
#' @param z Numeric vector of z-scores.
#' @param n Numeric, GWAS sample size.
#' @param eigen_ref An \code{LDEigen} object.
#' @param annotations An \code{AnnotationMatrix}, or NULL.
#' @param tau Numeric vector of annotation coefficients.
#' @param a Numeric, intercept.
#' @param baseline_mat Matrix of baseline annotations, or NULL.
#' @return A list with enrichment data.frame and score_stats list.
#' @keywords internal
.lder_stratified <- function(z, n, eigen_ref, annotations, tau, a,
                             baseline_mat = NULL) {
  # Score-based approach: fit baseline jointly, compute scores for candidates
  candidate_annot <- getCandidates(annotations)
  n_candidates <- ncol(candidate_annot@annotations)

  if (n_candidates == 0) {
    return(list(enrichment = NULL, score_stats = NULL))
  }

  n_blocks <- length(eigen_ref@eigen_list)
  M <- nrow(eigen_ref@snp_info)

  # Collect per-block partial scores into a matrix (n_blocks x n_candidates)
  partials_mat <- matrix(0, nrow = n_blocks, ncol = n_candidates)

  for (b in seq_len(n_blocks)) {
    block <- eigen_ref@eigen_list[[b]]
    idx <- block$snp_idx
    V <- block$vectors
    d <- block$values
    z_block <- z[idx]
    z_rot <- as.vector(t(V) %*% z_block)
    chi2_rot <- z_rot^2

    # Compute residual from stratified baseline fit
    if (!is.null(baseline_mat)) {
      ld_annot_base <- crossprod(V^2, baseline_mat[idx, , drop = FALSE])
      fitted_baseline <- as.vector(n / M * d *
                                     (ld_annot_base %*% tau))
    } else {
      fitted_baseline <- n * tau[1] * d / M
    }
    residual <- chi2_rot - 1 - fitted_baseline - n * a
    w <- 1 / (2 * pmax(chi2_rot, 1)^2)

    for (ai in seq_len(n_candidates)) {
      annot_col <- candidate_annot@annotations[, ai]
      annot_block <- annot_col[idx]
      ld_annot <- as.vector(t(V^2) %*% annot_block)

      partials_mat[b, ai] <- sum(w * residual * n * ld_annot / M)
    }
  }

  # Compute score_z from block partials
  score_z <- colSums(partials_mat) /
    sqrt(colSums(partials_mat^2) - colSums(partials_mat)^2 / n_blocks)

  # Score correlation matrix via jackknife
  # For each LOO iteration, recompute score_z excluding one block
  loo_score_z <- matrix(0, nrow = n_blocks, ncol = n_candidates)
  for (b in seq_len(n_blocks)) {
    partials_loo <- partials_mat[-b, , drop = FALSE]
    n_loo <- n_blocks - 1
    loo_score_z[b, ] <- colSums(partials_loo) /
      sqrt(colSums(partials_loo^2) - colSums(partials_loo)^2 / n_loo)
  }
  if (n_candidates > 1) {
    R <- cor(loo_score_z)
  } else {
    R <- matrix(1, 1, 1)
  }

  enrichment_df <- data.frame(
    annotation = candidate_annot@annotation_meta$name,
    score_z = score_z,
    score_p = 2 * pnorm(-abs(score_z)),
    stringsAsFactors = FALSE
  )

  score_stats_list <- list(
    z = score_z,
    R = R,
    annotation_names = candidate_annot@annotation_meta$name
  )

  list(enrichment = enrichment_df, score_stats = score_stats_list)
}
