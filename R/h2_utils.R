#' @title Shared Utilities for Heritability Estimation
#' @description Internal helper functions for block operations, regression,
#'   jackknife SE, enrichment computation, and meta-analysis.
#' @importFrom GenomicRanges GRanges
#' @importFrom BiocParallel bplapply bpparam
NULL

# =============================================================================
# Block-level operations
# =============================================================================

#' @title Get SNP Indices Per Block
#' @description For each LD block, find the SNP indices from a reference
#'   that fall within the block boundaries.
#' @param snp_info A data.frame with columns CHR, BP.
#' @param ld_blocks An \code{LDBlocks} object.
#' @return A list of integer vectors, one per block.
#' @keywords internal
snpsPerBlock <- function(snp_info, ld_blocks) {
  blocks_gr <- ld_blocks@blocks
  snp_gr <- GRanges(
    seqnames = snp_info$CHR,
    ranges = IRanges(start = snp_info$BP, width = 1L)
  )
  hits <- findOverlaps(snp_gr, blocks_gr)
  split(queryHits(hits), subjectHits(hits))
}

#' @title Apply Function Per Block with BiocParallel
#' @description Apply a function to each LD block in parallel.
#' @param block_indices List of SNP index vectors per block.
#' @param FUN Function to apply to each block's indices.
#' @param BPPARAM BiocParallel parameter object.
#' @param ... Additional arguments passed to FUN.
#' @return A list of results, one per block.
#' @keywords internal
bplapplyBlocks <- function(block_indices, FUN, BPPARAM = NULL, ...) {
  if (is.null(BPPARAM)) {
    BPPARAM <- bpparam()
  }
  bplapply(block_indices, FUN, BPPARAM = BPPARAM, ...)
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
#' @param estimates_full Numeric vector, full-sample parameter estimates.
#' @param estimates_loo A matrix (n_blocks x n_params), leave-one-out estimates.
#' @return Numeric vector of jackknife SEs.
#' @keywords internal
jackknifeSe <- function(estimates_full, estimates_loo) {
  n_blocks <- nrow(estimates_loo)
  pseudo_vals <- n_blocks * matrix(estimates_full, nrow = n_blocks,
                                   ncol = length(estimates_full),
                                   byrow = TRUE) -
    (n_blocks - 1) * estimates_loo
  jk_var <- apply(pseudo_vals, 2, var) / n_blocks
  sqrt(jk_var)
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
#' @param penalize_intercept Logical. If FALSE (default), the last column
#'   of X (assumed to be the intercept) is not penalized.
#' @return Same structure as \code{weightedLs}: coef, se, residuals, fitted,
#'   vcov.
#' @keywords internal
weightedLsRidge <- function(y, X, w, lambda = 0,
                            penalize_intercept = FALSE) {
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
  if (!penalize_intercept && p > 1) penalty[p, p] <- 0
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
#'   compute the full set of enrichment quantities: prop_h2, prop_snps,
#'   enrichment ratio, enrichment SE (from jackknife or delta method),
#'   and p-value.
#' @param tau Numeric vector of per-annotation regression coefficients.
#' @param tau_se Numeric vector of SE for tau.
#' @param tau_blocks Numeric matrix (n_blocks x n_annotations) of jackknife
#'   block-level tau values, or NULL.
#' @param baseline_mat Numeric matrix (n_snps x n_annotations).
#' @param annot_names Character vector of annotation names.
#' @param h2 Numeric scalar, total estimated h2.
#' @return A data.frame with columns: annotation, tau, tau_se, enrichment,
#'   enrichment_se, enrichment_p, prop_h2, prop_snps.
#' @keywords internal
computeBaselineEnrichment <- function(tau, tau_se, tau_blocks,
                                      baseline_mat, annot_names, h2) {
  M <- nrow(baseline_mat)
  M_a <- colSums(baseline_mat)
  prop_snps <- M_a / M

  # Per-annotation h2 and proportion
  h2_a <- tau * M_a
  prop_h2 <- h2_a / h2

  # Enrichment ratio: (prop_h2 / prop_snps) = tau * M / h2
  enrichment <- tau * M / h2

  # Enrichment SE from jackknife blocks (preferred) or delta method (fallback)
  if (!is.null(tau_blocks)) {
    n_blocks <- nrow(tau_blocks)
    # Per-block enrichment: enrichment_b = tau_b * M / h2_b
    h2_blocks <- as.vector(tau_blocks %*% M_a)
    # Avoid division by zero
    h2_blocks[h2_blocks == 0] <- NA
    enrichment_blocks <- sweep(tau_blocks, 1, h2_blocks, FUN = "/") * M
    # Jackknife variance: Var = (B-1)/B * sum((x_b - x_bar)^2)
    enrichment_mean <- colMeans(enrichment_blocks, na.rm = TRUE)
    enrichment_var <- (n_blocks - 1) / n_blocks *
      colSums(sweep(enrichment_blocks, 2, enrichment_mean)^2, na.rm = TRUE)
    enrichment_se <- sqrt(enrichment_var)
  } else {
    # Delta method fallback: d(enrichment)/d(tau) = M / h2
    enrichment_se <- tau_se * M / abs(h2)
  }

  # P-value from z-score
  enrichment_z <- enrichment / enrichment_se
  enrichment_p <- 2 * pnorm(-abs(enrichment_z))

  data.frame(
    annotation = annot_names,
    tau = tau,
    tau_se = tau_se,
    enrichment = enrichment,
    enrichment_se = enrichment_se,
    enrichment_p = enrichment_p,
    prop_h2 = prop_h2,
    prop_snps = prop_snps,
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
#' @param n_ref Integer, reference panel sample size.
#' @param shrinkage_type Character, one of "wen_stephens", "constant".
#' @param genetic_map Numeric vector, genetic map positions for SNPs in R.
#' @return Shrunk LD correlation matrix.
#' @keywords internal
shrinkLd <- function(R, n_ref, shrinkage_type = "wen_stephens",
                      genetic_map = NULL) {
  if (shrinkage_type == "wen_stephens" && !is.null(genetic_map)) {
    # Wen & Stephens (2010) shrinkage based on genetic distance
    p <- nrow(R)
    theta <- 2 * n_ref / (22 * n_ref + 16)  # effective recombination
    dist_cm <- abs(outer(genetic_map, genetic_map, "-"))
    shrink_factor <- exp(-4 * n_ref * dist_cm / (100 * (2 * n_ref + 16)))
    R_shrunk <- R * shrink_factor
    diag(R_shrunk) <- 1
  } else {
    # Simple constant shrinkage
    lambda <- 1 / sqrt(n_ref)
    R_shrunk <- (1 - lambda) * R + lambda * diag(nrow(R))
  }
  R_shrunk
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
    if (is(x, "GWASSumStats")) x@genome
    else if (is(x, "LDStatistic")) x@genome
    else if (is(x, "AnnotationMatrix")) x@genome
    else if (is(x, "LDBlocks")) x@genome
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
#' @param tau_blocks Numeric matrix (n_blocks x n_annotations) of
#'   block-level tau estimates from delete-one jackknife.
#' @param sd_annot Numeric vector of per-annotation standard deviations,
#'   same length as \code{tau}.
#' @param M_ref Scalar integer, total number of reference-panel SNPs.
#' @param h2g Numeric scalar, total estimated SNP heritability.
#' @return A list with:
#'   \describe{
#'     \item{tau_star}{Numeric vector of standardized tau values.}
#'     \item{tau_star_se}{Numeric vector of jackknife SE for tau_star.}
#'   }
#' @keywords internal
standardize_tau_star <- function(tau, tau_blocks, sd_annot, M_ref, h2g) {
  if (length(tau) != length(sd_annot)) {
    stop("standardize_tau_star: tau and sd_annot must have the same length.")
  }
  if (h2g == 0) {
    stop("standardize_tau_star: h2g must be non-zero.")
  }

  # Gazal standardization: tau* = tau * sd_annot * M_ref / h2g
  coef <- sd_annot * M_ref / h2g
  tau_star <- tau * coef

  # Jackknife SE from block-level tau
  tau_star_blocks <- sweep(tau_blocks, 2L, coef, FUN = "*")
  n_blocks <- nrow(tau_star_blocks)
  jk_var <- apply(tau_star_blocks, 2L, function(x) var(x, na.rm = TRUE))
  tau_star_se <- sqrt((n_blocks - 1)^2 / n_blocks * jk_var)

  list(tau_star = tau_star, tau_star_se = tau_star_se)
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
meta_random_effects <- function(means, ses) {
  k <- length(means)
  if (k != length(ses)) {
    stop("meta_random_effects: means and ses must have the same length.")
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
    stop("meta_random_effects: all ses must be positive and finite.")
  }

  # Fixed-effect weights
  w_fe <- 1 / ses^2

  # Fixed-effect pooled estimate
  mu_fe <- sum(w_fe * means) / sum(w_fe)

  # Cochran's Q
  Q <- sum(w_fe * (means - mu_fe)^2)

  # DerSimonian-Laird tau-squared estimator
  c_dl <- sum(w_fe) - sum(w_fe^2) / sum(w_fe)
  tau2 <- max(0, (Q - (k - 1)) / c_dl)

  # Random-effects weights

  w_re <- 1 / (ses^2 + tau2)

  # Pooled random-effects estimate
  mu_re <- sum(w_re * means) / sum(w_re)
  se_re <- sqrt(1 / sum(w_re))

  # Higgins I-squared
  I2 <- max(0, (Q - (k - 1)) / Q)

  list(mean = mu_re, se = se_re, tau2 = tau2, I2 = I2, Q = Q)
}
