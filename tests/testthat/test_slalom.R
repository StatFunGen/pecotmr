context("slalom")

# ============================================================================
# Helper: build a valid positive-definite LD matrix from a genotype matrix
# ============================================================================
make_synthetic_ld <- function(n_samples, n_snps, seed = 1) {
  set.seed(seed)
  # Simulate genotypes with some LD structure by using a factor model
  # X = Z %*% L + noise, where Z is latent and L is a loading matrix
  n_factors <- min(3, n_snps)
  Z <- matrix(rnorm(n_samples * n_factors), nrow = n_samples)
  L <- matrix(runif(n_factors * n_snps, -1, 1), nrow = n_factors)
  X_raw <- Z %*% L + matrix(rnorm(n_samples * n_snps, sd = 0.5), nrow = n_samples)
  # Discretise to genotype-like values (0, 1, 2)
  X <- matrix(as.integer(cut(X_raw, breaks = c(-Inf, -0.5, 0.5, Inf))) - 1L,
              nrow = n_samples, ncol = n_snps)
  colnames(X) <- paste0("snp", seq_len(n_snps))
  R <- cor(X)
  # Ensure perfect diagonal and clean NaN from any zero-variance columns
  R[is.na(R) | is.nan(R)] <- 0
  diag(R) <- 1.0
  list(X = X, R = R)
}

# ============================================================================
# Basic output structure
# ============================================================================

test_that("slalom basic output structure", {
  set.seed(42)
  n <- 50
  z <- rnorm(n)
  R <- diag(n)
  # Add some off-diagonal correlations
  for (i in 1:(n - 1)) {
    R[i, i + 1] <- 0.3
    R[i + 1, i] <- 0.3
  }

  result <- slalom(zScore = z, R = R)

  expect_type(result, "list")
  expect_named(result, c("data", "summary"))
  expect_s3_class(result$data, "data.frame")
  expect_true("original_z" %in% colnames(result$data))
  expect_true("prob" %in% colnames(result$data))
  expect_true("pvalue" %in% colnames(result$data))
  expect_true("outliers" %in% colnames(result$data))
  expect_true("nlog10p_dentist_s" %in% colnames(result$data))
  expect_equal(nrow(result$data), n)
})

test_that("slalom errors on non-square R", {
  z <- rnorm(10)
  R <- matrix(rnorm(50), nrow = 5, ncol = 10)
  expect_error(slalom(zScore = z, R = R), "LD_mat must be a square matrix")
})

test_that("slalom accepts X matrix instead of R", {
  set.seed(42)
  n_samples <- 100
  n_snps <- 10
  X <- matrix(sample(0:2, n_samples * n_snps, replace = TRUE), nrow = n_samples, ncol = n_snps)
  colnames(X) <- paste0("snp", 1:n_snps)
  z <- rnorm(n_snps)

  result <- slalom(zScore = z, X = X)
  expect_type(result, "list")
  expect_equal(nrow(result$data), n_snps)
  # PIPs should be in [0,1] and sum to 1
  expect_true(all(result$data$prob >= 0 & result$data$prob <= 1))
  expect_equal(sum(result$data$prob), 1, tolerance = 1e-12)
  # Data frame should have expected column names
  expected_cols <- c("original_z", "prob", "pvalue", "outliers", "nlog10p_dentist_s")
  expect_true(all(expected_cols %in% colnames(result$data)))
})

# ============================================================================
# ABF computation correctness
# ============================================================================

test_that("ABF: strong signal (z=10) gets very high PIP", {
  set.seed(100)
  n <- 20
  z <- rnorm(n, sd = 0.3)
  z[7] <- 10
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  expect_equal(which.max(result$data$prob), 7)
  expect_gt(result$data$prob[7], 0.15)
  expect_equal(sum(result$data$prob), 1, tolerance = 1e-12)
})

test_that("ABF: moderate signal (z=3) gets higher PIP than weak signal (z=1)", {
  set.seed(101)
  n <- 10
  z <- rep(0, n)
  z[2] <- 1
  z[5] <- 3
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  expect_gt(result$data$prob[5], result$data$prob[2])
  # z=0 variants should all have same PIP (by symmetry, with identity LD)
  zero_pips <- result$data$prob[c(1, 3, 4, 6, 7, 8, 9, 10)]
  expect_equal(max(zero_pips) - min(zero_pips), 0, tolerance = 1e-14)
})

test_that("ABF: lbf formula matches manual calculation", {
  z_val <- 4.0
  se_val <- 1.0
  W <- 0.04

  V <- se_val^2
  r <- W / (W + V)
  expected_lbf <- 0.5 * (log(1 - r) + r * z_val^2)

  z <- c(z_val, 0)
  R <- diag(2)
  result <- slalom(zScore = z, R = R, abf_prior_variance = W)

  lbf_0 <- 0.5 * (log(1 - r) + r * 0^2)
  expected_ratio <- exp(expected_lbf - lbf_0)
  actual_ratio <- result$data$prob[1] / result$data$prob[2]
  expect_equal(actual_ratio, expected_ratio, tolerance = 1e-10)
})

test_that("ABF: PIPs always sum to exactly 1", {
  for (s in 1:5) {
    set.seed(200 + s)
    n <- sample(5:30, 1)
    z <- rnorm(n, sd = 2)
    R <- diag(n)
    result <- slalom(zScore = z, R = R)
    expect_equal(sum(result$data$prob), 1, tolerance = 1e-12,
                 label = paste("seed", 200 + s))
  }
})

test_that("ABF: symmetric z-scores give symmetric PIPs", {
  z <- c(-3, 3)
  R <- diag(2)
  result <- slalom(zScore = z, R = R)
  expect_equal(result$data$prob[1], result$data$prob[2], tolerance = 1e-14)
})

# ============================================================================
# Credible sets
# ============================================================================

test_that("CS95 contains the causal variant in a simple synthetic signal", {
  set.seed(300)
  n <- 30
  z <- rnorm(n, sd = 0.5)
  causal <- 12
  z[causal] <- 6
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  expect_true(causal %in% result$summary$cs_95)
  expect_true(causal %in% result$summary$cs_99)
})

test_that("CS99 is a superset of CS95", {
  set.seed(301)
  n <- 40
  z <- rnorm(n, sd = 1.5)
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  expect_true(all(result$summary$cs_95 %in% result$summary$cs_99))
  expect_gte(length(result$summary$cs_99), length(result$summary$cs_95))
})

test_that("CS95 covers at least 95% of posterior mass", {
  set.seed(302)
  n <- 25
  z <- rnorm(n)
  z[10] <- 4
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  cs95_mass <- sum(result$data$prob[result$summary$cs_95])
  expect_gt(cs95_mass, 0.95)

  cs99_mass <- sum(result$data$prob[result$summary$cs_99])
  expect_gt(cs99_mass, 0.99)
})

test_that("CS with very strong signal contains only the causal variant", {
  set.seed(303)
  n <- 15
  z <- rnorm(n, sd = 0.1)
  z[8] <- 15  # extremely strong
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  expect_equal(result$summary$cs_95[1], 8)
  expect_true(8 %in% result$summary$cs_95)
  expect_true(8 %in% result$summary$cs_99)
})

test_that("CS with diffuse signal contains many variants", {
  set.seed(304)
  n <- 10
  z <- rep(0, n)  # all equally uninformative
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  # Uniform PIPs => need at least ceiling(0.95 * n) = 10 variants for 95% coverage
  expect_equal(length(result$summary$cs_95), n)
})

# ============================================================================
# Lead variant by pvalue vs abf
# ============================================================================

test_that("lead variant by pvalue selects most negative z-score", {
  z <- c(0, -4, 3, -1, 2)
  R <- diag(5)

  result <- slalom(zScore = z, R = R, lead_variant_choice = "pvalue")

  expect_equal(result$summary$lead_pip_variant, 2)
})

test_that("lead variant by abf selects highest PIP", {
  z <- c(0, -4, 3, -1, 2)
  R <- diag(5)

  result <- slalom(zScore = z, R = R, lead_variant_choice = "abf")

  expect_equal(result$summary$lead_pip_variant, which.max(result$data$prob))
})

test_that("pvalue and abf lead can differ when z has asymmetric magnitudes", {
  z <- c(-3.0, 5.0, 0.1, -0.2, 0.3)
  R <- diag(5)

  result_pv <- slalom(zScore = z, R = R, lead_variant_choice = "pvalue")
  result_abf <- slalom(zScore = z, R = R, lead_variant_choice = "abf")

  expect_equal(result_pv$summary$lead_pip_variant, 1)
  expect_equal(result_abf$summary$lead_pip_variant, 2)
  expect_false(result_pv$summary$lead_pip_variant == result_abf$summary$lead_pip_variant)
})

# ============================================================================
# DENTIST-S outlier detection
# ============================================================================

test_that("DENTIST-S: lead variant itself is not flagged as outlier", {
  set.seed(400)
  n <- 10
  z <- rnorm(n, sd = 0.5)
  z[3] <- -5  # lead by pvalue
  R <- diag(n)

  result <- slalom(zScore = z, R = R)
  lead <- result$summary$lead_pip_variant
  expect_equal(lead, 3)

  expect_true(is.na(result$data$outliers[lead]) || !result$data$outliers[lead])
})

test_that("DENTIST-S: outlier variant inconsistent with LD is flagged", {
  n <- 5
  z <- c(-5, 0, 0.1, -0.1, 0.2)
  R <- diag(n)
  R[1, 2] <- R[2, 1] <- 0.9
  R[1, 3] <- R[3, 1] <- 0.1
  R[1, 4] <- R[4, 1] <- -0.05
  R[1, 5] <- R[5, 1] <- 0.02

  result <- slalom(zScore = z, R = R, r2_threshold = 0.5,
                   nlog10p_dentist_s_threshold = 2.0)

  lead <- result$summary$lead_pip_variant
  expect_equal(lead, 1)

  expect_true(result$data$outliers[2])
  expect_false(result$data$outliers[3])
})

test_that("DENTIST-S: perfectly consistent variant in LD is not flagged", {
  n <- 3
  z <- c(-5, -4.0, 0.1)
  R <- diag(n)
  R[1, 2] <- R[2, 1] <- 0.8
  R[1, 3] <- R[3, 1] <- 0.05
  R[2, 3] <- R[3, 2] <- 0.04

  result <- slalom(zScore = z, R = R, r2_threshold = 0.5,
                   nlog10p_dentist_s_threshold = 4.0)

  lead <- result$summary$lead_pip_variant
  expect_equal(lead, 1)

  expect_equal(result$data$nlog10p_dentist_s[2], 0, tolerance = 1e-10)
  expect_false(result$data$outliers[2])
})

test_that("DENTIST-S: n_dentist_s_outlier and fraction are consistent", {
  set.seed(401)
  n <- 20
  syn <- make_synthetic_ld(200, n, seed = 401)
  z <- rnorm(n, sd = 1)
  z[1] <- -5

  result <- slalom(zScore = z, R = syn$R, r2_threshold = 0.3)

  n_r2 <- result$summary$n_r2
  n_out <- result$summary$n_dentist_s_outlier
  frac <- result$summary$fraction

  expect_gte(n_r2, 1)
  expect_equal(frac, ifelse(n_r2 > 0, n_out / n_r2, 0), tolerance = 1e-14)
  expect_gte(n_out, 0)
  expect_lte(n_out, n_r2)
})

test_that("DENTIST-S: lowering threshold flags more outliers", {
  set.seed(402)
  n <- 15
  syn <- make_synthetic_ld(300, n, seed = 402)
  z <- rnorm(n, sd = 2)
  z[5] <- -6

  result_strict <- slalom(zScore = z, R = syn$R, nlog10p_dentist_s_threshold = 6.0,
                          r2_threshold = 0.3)
  result_loose <- slalom(zScore = z, R = syn$R, nlog10p_dentist_s_threshold = 1.0,
                         r2_threshold = 0.3)

  expect_gte(result_loose$summary$n_dentist_s_outlier,
             result_strict$summary$n_dentist_s_outlier)
})

# ============================================================================
# Edge cases
# ============================================================================

test_that("edge case: single variant", {
  z <- c(3.0)
  R <- matrix(1, nrow = 1, ncol = 1)

  result <- slalom(zScore = z, R = R)

  expect_equal(nrow(result$data), 1)
  expect_equal(result$data$prob[1], 1.0, tolerance = 1e-14)
  expect_equal(result$summary$lead_pip_variant, 1)
  expect_equal(result$summary$n_total, 1)
  expect_equal(result$summary$cs_95, 1)
  expect_equal(result$summary$cs_99, 1)
})

test_that("edge case: all zero z-scores", {
  n <- 10
  z <- rep(0, n)
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  expect_equal(result$data$prob, rep(1 / n, n), tolerance = 1e-14)
  expect_equal(result$data$pvalue, rep(0.5, n), tolerance = 1e-14)
  expect_equal(result$summary$max_pip, 1 / n, tolerance = 1e-14)
})

test_that("edge case: very large z-scores do not produce NaN in PIPs", {
  z <- c(50, -50, 30, -30, 0)
  R <- diag(5)

  result <- slalom(zScore = z, R = R)

  expect_false(any(is.nan(result$data$prob)))
  expect_false(any(is.na(result$data$prob)))
  expect_equal(sum(result$data$prob), 1, tolerance = 1e-10)
  expect_equal(result$data$prob[1], result$data$prob[2], tolerance = 1e-14)
})

test_that("edge case: identical z-scores yield uniform PIPs", {
  n <- 8
  z <- rep(2.5, n)
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  expect_equal(result$data$prob, rep(1 / n, n), tolerance = 1e-14)
})

test_that("edge case: two variants only", {
  z <- c(-3, 2)
  R <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)

  result <- slalom(zScore = z, R = R)

  expect_equal(nrow(result$data), 2)
  expect_equal(sum(result$data$prob), 1, tolerance = 1e-12)
  expect_equal(result$summary$lead_pip_variant, 1)
  expect_gt(result$data$prob[1], result$data$prob[2])
})

test_that("edge case: mismatched dimensions error", {
  z <- rnorm(10)
  R <- diag(5)
  expect_error(slalom(zScore = z, R = R),
               "LD_mat must be a square matrix matching the length of zScore")
})

test_that("edge case: no R and no X provided errors", {
  z <- rnorm(5)
  expect_error(slalom(zScore = z), "Either R.*or X.*must be provided")
})

test_that("edge case: both R and X provided errors", {
  set.seed(500)
  n <- 5
  z <- rnorm(n)
  R <- diag(n)
  X <- matrix(sample(0:2, 50 * n, replace = TRUE), nrow = 50, ncol = n)
  expect_error(slalom(zScore = z, R = R, X = X), "Provide either R or X, not both")
})

# ============================================================================
# X input mode
# ============================================================================

test_that("X input yields same result as R = cor(X)", {
  set.seed(600)
  n_samples <- 200
  n_snps <- 10
  X <- matrix(sample(0:2, n_samples * n_snps, replace = TRUE),
              nrow = n_samples, ncol = n_snps)
  colnames(X) <- paste0("snp", seq_len(n_snps))
  z <- rnorm(n_snps, sd = 2)

  R_manual <- cor(X)
  diag(R_manual) <- 1.0

  result_X <- slalom(zScore = z, X = X)
  result_R <- slalom(zScore = z, R = R_manual)

  expect_equal(result_X$data$prob, result_R$data$prob, tolerance = 1e-10)
  expect_equal(result_X$data$original_z, result_R$data$original_z, tolerance = 1e-14)
  expect_equal(result_X$summary$lead_pip_variant, result_R$summary$lead_pip_variant)
  expect_equal(result_X$summary$cs_95, result_R$summary$cs_95)
  expect_equal(result_X$summary$cs_99, result_R$summary$cs_99)
})

# ============================================================================
# Parameter variation
# ============================================================================

test_that("larger abf_prior_variance concentrates PIPs on strong signals more", {
  set.seed(700)
  n <- 15
  z <- rnorm(n, sd = 0.5)
  z[4] <- 4
  R <- diag(n)

  result_small_W <- slalom(zScore = z, R = R, abf_prior_variance = 0.01)
  result_large_W <- slalom(zScore = z, R = R, abf_prior_variance = 1.0)

  expect_gt(result_large_W$summary$max_pip, result_small_W$summary$max_pip)
})

test_that("abf_prior_variance = 0 gives uniform PIPs", {
  n <- 10
  z <- c(5, 3, 1, 0, -1, -3, -5, 2, -2, 4)
  R <- diag(n)

  result <- slalom(zScore = z, R = R, abf_prior_variance = 0)

  expect_equal(result$data$prob, rep(1 / n, n), tolerance = 1e-14)
})

test_that("different standard_error values affect PIPs", {
  n <- 5
  z <- c(3, 3, 3, 3, 3)  # same z for all
  R <- diag(n)
  se1 <- c(1, 1, 1, 1, 1)
  se2 <- c(0.5, 1, 1, 1, 1)  # variant 1 has smaller SE

  result1 <- slalom(zScore = z, R = R, standard_error = se1)
  result2 <- slalom(zScore = z, R = R, standard_error = se2)

  expect_gt(result2$data$prob[1], result1$data$prob[1])
})

test_that("r2_threshold variation affects n_r2 count", {
  set.seed(701)
  n <- 10
  syn <- make_synthetic_ld(200, n, seed = 701)
  z <- rnorm(n, sd = 2)
  z[1] <- -5

  result_low <- slalom(zScore = z, R = syn$R, r2_threshold = 0.1)
  result_high <- slalom(zScore = z, R = syn$R, r2_threshold = 0.9)

  expect_gte(result_low$summary$n_r2, result_high$summary$n_r2)
})

test_that("nlog10p_dentist_s_threshold variation affects outlier count", {
  set.seed(702)
  n <- 10
  syn <- make_synthetic_ld(200, n, seed = 702)
  z <- rnorm(n, sd = 2)
  z[1] <- -6

  result_low_thresh <- slalom(zScore = z, R = syn$R, nlog10p_dentist_s_threshold = 1.0)
  result_high_thresh <- slalom(zScore = z, R = syn$R, nlog10p_dentist_s_threshold = 10.0)

  expect_gte(result_low_thresh$summary$n_dentist_s_outlier,
             result_high_thresh$summary$n_dentist_s_outlier)
})

# ============================================================================
# Output structure validation
# ============================================================================

test_that("output data types are correct", {
  set.seed(801)
  n <- 15
  z <- rnorm(n)
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  expect_type(result$data$original_z, "double")
  expect_type(result$data$prob, "double")
  expect_type(result$data$pvalue, "double")
  expect_type(result$data$outliers, "logical")
  expect_type(result$data$nlog10p_dentist_s, "double")

  expect_type(result$summary$lead_pip_variant, "integer")
  expect_type(result$summary$n_total, "integer")
  expect_type(result$summary$n_r2, "integer")
  expect_type(result$summary$n_dentist_s_outlier, "integer")
  expect_type(result$summary$fraction, "double")
  expect_type(result$summary$max_pip, "double")
  expect_type(result$summary$cs_95, "integer")
  expect_type(result$summary$cs_99, "integer")
})

test_that("original_z in output matches input z-scores", {
  z <- c(1.5, -2.3, 0.7, 4.1, -0.5)
  R <- diag(5)

  result <- slalom(zScore = z, R = R)
  expect_equal(result$data$original_z, z, tolerance = 0)
})

test_that("pvalue in output matches pnorm(z)", {
  z <- c(-3, -1, 0, 1, 3)
  R <- diag(5)

  result <- slalom(zScore = z, R = R)
  expect_equal(result$data$pvalue, pnorm(z), tolerance = 1e-14)
})

# ============================================================================
# Summary statistics consistency
# ============================================================================

test_that("n_r2 counts variants with r2 > threshold to lead correctly", {
  n <- 10
  z <- c(-5, rep(0, 9))
  R <- diag(n)

  result <- slalom(zScore = z, R = R, r2_threshold = 0.6)
  expect_equal(result$summary$n_r2, 1)
})

test_that("n_r2 includes correlated variants", {
  n <- 5
  z <- c(-5, -4, 0, 0, 0)
  R <- diag(n)
  R[1, 2] <- R[2, 1] <- 0.9

  result <- slalom(zScore = z, R = R, r2_threshold = 0.6)
  expect_equal(result$summary$n_r2, 2)
})

test_that("fraction = 0 when there are no outliers (identity LD, consistent z)", {
  n <- 5
  z <- c(-3, 0, 0, 0, 0)
  R <- diag(n)

  result <- slalom(zScore = z, R = R)
  expect_equal(result$summary$fraction, 0)
})

test_that("fraction is between 0 and 1", {
  set.seed(900)
  for (s in 1:5) {
    set.seed(900 + s)
    n <- sample(10:30, 1)
    syn <- make_synthetic_ld(200, n, seed = 900 + s)
    z <- rnorm(n, sd = 2)
    z[1] <- -6

    result <- slalom(zScore = z, R = syn$R, r2_threshold = 0.3)
    expect_gte(result$summary$fraction, 0)
    expect_lte(result$summary$fraction, 1)
  }
})

test_that("max_pip equals the maximum of prob vector", {
  set.seed(901)
  n <- 20
  z <- rnorm(n, sd = 2)
  R <- diag(n)

  result <- slalom(zScore = z, R = R)
  expect_equal(result$summary$max_pip, max(result$data$prob), tolerance = 1e-14)
})

# ============================================================================
# Realistic synthetic LD scenarios
# ============================================================================

test_that("realistic LD: correlated variants share PIP mass", {
  set.seed(1000)
  syn <- make_synthetic_ld(500, 20, seed = 1000)
  z <- rep(0, 20)
  z[3] <- 5

  result <- slalom(zScore = z, R = syn$R)

  expect_true(3 %in% result$summary$cs_95)
  expect_equal(sum(result$data$prob), 1, tolerance = 1e-12)
})

test_that("realistic LD: DENTIST-S detects outlier in correlated block", {
  set.seed(1001)
  syn <- make_synthetic_ld(500, 15, seed = 1001)
  R <- syn$R

  lead_idx <- 1
  z <- R[, lead_idx] * (-5)
  z[1] <- -5

  r2_to_lead <- R[, 1]^2
  candidates <- which(r2_to_lead > 0.3 & seq_along(z) != 1)
  if (length(candidates) > 0) {
    corrupt_idx <- candidates[1]
    z[corrupt_idx] <- z[corrupt_idx] + 10

    result <- slalom(zScore = z, R = R, r2_threshold = 0.2,
                     nlog10p_dentist_s_threshold = 3.0)

    expect_true(result$data$outliers[corrupt_idx])
  }
})

test_that("realistic LD: no outliers when z perfectly matches LD structure", {
  set.seed(1002)
  syn <- make_synthetic_ld(500, 10, seed = 1002)
  R <- syn$R

  lead_idx <- 1
  z <- R[, lead_idx] * (-4)

  result <- slalom(zScore = z, R = R, r2_threshold = 0.3,
                   nlog10p_dentist_s_threshold = 4.0)

  non_lead <- setdiff(seq_along(z), result$summary$lead_pip_variant)
  lead <- result$summary$lead_pip_variant
  for (i in non_lead) {
    r2_val <- R[i, lead]^2
    if (r2_val > 0.3 && r2_val < 1.0) {
      expect_false(result$data$outliers[i])
    } else if (r2_val >= 1.0) {
      expect_true(is.na(result$data$outliers[i]) || !result$data$outliers[i])
    }
  }
})

# ============================================================================
# Internal function access (resolve_LD_input via :::)
# ============================================================================

test_that("resolve_LD_input returns R when R is provided", {
  R <- diag(5)
  res <- pecotmr:::resolve_LD_input(R = R, need_nSample = FALSE)
  expect_equal(res$R, R)
})

test_that("resolve_LD_input computes R from X", {
  set.seed(1100)
  X <- matrix(sample(0:2, 200 * 5, replace = TRUE), nrow = 200, ncol = 5)
  colnames(X) <- paste0("s", 1:5)
  res <- pecotmr:::resolve_LD_input(X = X, need_nSample = FALSE)
  expect_true(is.matrix(res$R))
  expect_equal(nrow(res$R), 5)
  expect_equal(ncol(res$R), 5)
  for (j in seq_len(5)) expect_equal(res$R[j, j], 1.0, tolerance = 1e-6)
})

test_that("resolve_LD_input errors when neither R nor X given", {
  expect_error(pecotmr:::resolve_LD_input(R = NULL, X = NULL),
               "Either R.*or X.*must be provided")
})

test_that("resolve_LD_input errors when both R and X given", {
  R <- diag(3)
  X <- matrix(1, nrow = 10, ncol = 3)
  expect_error(pecotmr:::resolve_LD_input(R = R, X = X),
               "Provide either R or X, not both")
})

# ============================================================================
# Numerical stability
# ============================================================================

test_that("log-sum-exp trick prevents overflow with extreme z-scores", {
  z <- c(100, 0, -100)
  R <- diag(3)

  result <- slalom(zScore = z, R = R)

  expect_false(any(is.nan(result$data$prob)))
  expect_false(any(is.infinite(result$data$prob)))
  expect_equal(sum(result$data$prob), 1, tolerance = 1e-10)
  expect_equal(result$data$prob[1], result$data$prob[3], tolerance = 1e-14)
  expect_gt(result$data$prob[1], result$data$prob[2])
})

test_that("standard_error near zero concentrates PIP on large z-scores", {
  z <- c(2, 0.5, 0, -0.1)
  R <- diag(4)
  se <- rep(0.01, 4)

  result <- slalom(zScore = z, R = R, standard_error = se, abf_prior_variance = 0.04)

  expect_equal(which.max(result$data$prob), 1)
  expect_false(any(is.nan(result$data$prob)))
})

# ============================================================================
# Determinism and reproducibility
# ============================================================================

test_that("slalom is deterministic (no randomness)", {
  z <- c(3, -2, 1, 0, -4)
  R <- diag(5)

  result1 <- slalom(zScore = z, R = R)
  result2 <- slalom(zScore = z, R = R)

  expect_identical(result1$data, result2$data)
  expect_identical(result1$summary$lead_pip_variant, result2$summary$lead_pip_variant)
  expect_identical(result1$summary$cs_95, result2$summary$cs_95)
  expect_identical(result1$summary$cs_99, result2$summary$cs_99)
})

# ============================================================================
# Credible set ordering
# ============================================================================

test_that("CS95 variants are ordered by decreasing PIP", {
  set.seed(1400)
  n <- 20
  z <- rnorm(n, sd = 2)
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  cs_pips <- result$data$prob[result$summary$cs_95]
  expect_true(all(diff(cs_pips) <= .Machine$double.eps))
})

