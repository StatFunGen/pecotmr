context("h2_utils")

# =============================================================================
# weightedLs
# =============================================================================

test_that("weightedLs returns correct structure", {
  set.seed(42)
  n <- 50
  x <- rnorm(n)
  X <- cbind(x, 1)
  y <- 2 * x + 1 + rnorm(n, sd = 0.1)
  w <- rep(1, n)
  res <- pecotmr:::weightedLs(y, X, w)
  expect_true(is.list(res))
  expect_named(res, c("coef", "se", "residuals", "fitted", "vcov"))
  expect_length(res$coef, 2)
  expect_length(res$se, 2)
  expect_length(res$residuals, n)
  expect_length(res$fitted, n)
  expect_equal(dim(res$vcov), c(2, 2))
})

test_that("weightedLs recovers coefficients for simple linear model", {
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  X <- cbind(x, 1)
  y <- 2 * x + 1 + rnorm(n, sd = 0.1)
  w <- rep(1, n)
  res <- pecotmr:::weightedLs(y, X, w)
  expect_equal(res$coef[1], 2, tolerance = 0.1)
  expect_equal(res$coef[2], 1, tolerance = 0.1)
})

test_that("weightedLs converts vector X to matrix", {
  set.seed(42)
  n <- 30
  x <- rnorm(n)
  y <- 3 * x + rnorm(n, sd = 0.1)
  w <- rep(1, n)
  res <- pecotmr:::weightedLs(y, x, w)
  expect_length(res$coef, 1)
  expect_equal(res$coef[1], 3, tolerance = 0.2)
})

test_that("weightedLs fitted + residuals = y", {
  set.seed(42)
  n <- 30
  X <- cbind(rnorm(n), 1)
  y <- rnorm(n)
  w <- rep(1, n)
  res <- pecotmr:::weightedLs(y, X, w)
  expect_equal(res$fitted + res$residuals, y)
})

# =============================================================================
# jackknifeSe
# =============================================================================

test_that("jackknifeSe returns zero when all LOO estimates equal full", {
  n_blocks <- 5
  n_params <- 3
  estimates_full <- c(1.0, 2.0, 3.0)
  estimates_loo <- matrix(rep(estimates_full, each = n_blocks),
                          nrow = n_blocks, ncol = n_params)
  se <- pecotmr:::jackknifeSe(estimates_full, estimates_loo)
  expect_length(se, n_params)
  expect_equal(se, rep(0, n_params))
})

test_that("jackknifeSe returns positive SE with varying LOO estimates", {
  set.seed(42)
  n_blocks <- 10
  estimates_full <- c(5.0, 10.0)
  estimates_loo <- matrix(rnorm(n_blocks * 2, mean = rep(estimates_full, each = n_blocks), sd = 0.5),
                          nrow = n_blocks, ncol = 2)
  se <- pecotmr:::jackknifeSe(estimates_full, estimates_loo)
  expect_length(se, 2)
  expect_true(all(se > 0))
})

test_that("jackknifeSe computes known case correctly", {
  # If full estimate is the mean of 5 values and LOO removes each one,
  # we can verify the formula
  vals <- c(2, 4, 6, 8, 10)
  full_mean <- mean(vals)
  n <- length(vals)
  # Leave-one-out means: remove element i, compute mean of remaining
  loo_means <- vapply(seq_len(n), function(i) mean(vals[-i]), numeric(1))
  estimates_loo <- matrix(loo_means, ncol = 1)
  se <- pecotmr:::jackknifeSe(full_mean, estimates_loo)
  # Pseudo-values: n * full - (n-1) * loo
  pseudo <- n * full_mean - (n - 1) * loo_means
  expected_se <- sqrt(var(pseudo) / n)
  expect_equal(se, expected_se)
})

# =============================================================================
# weightedLsRidge
# =============================================================================

test_that("weightedLsRidge with lambda=0 matches weightedLs", {
  set.seed(42)
  n <- 50
  X <- cbind(rnorm(n), rnorm(n), 1)
  y <- X %*% c(1, 2, 0.5) + rnorm(n, sd = 0.1)
  w <- rep(1, n)
  res_ridge <- pecotmr:::weightedLsRidge(y, X, w, lambda = 0)
  res_wls <- pecotmr:::weightedLs(y, X, w)
  expect_equal(res_ridge$coef, res_wls$coef)
  expect_equal(res_ridge$se, res_wls$se)
  expect_equal(res_ridge$residuals, res_wls$residuals)
  expect_equal(res_ridge$fitted, res_wls$fitted)
})

test_that("weightedLsRidge with lambda>0 shrinks coefficients", {
  set.seed(42)
  n <- 50
  X <- cbind(rnorm(n), rnorm(n), 1)
  y <- X %*% c(3, 3, 1) + rnorm(n, sd = 0.5)
  w <- rep(1, n)
  res_no_ridge <- pecotmr:::weightedLsRidge(y, X, w, lambda = 0)
  res_ridge <- pecotmr:::weightedLsRidge(y, X, w, lambda = 10)
  # Non-intercept coefficients should be shrunk toward zero
  expect_true(abs(res_ridge$coef[1]) < abs(res_no_ridge$coef[1]))
  expect_true(abs(res_ridge$coef[2]) < abs(res_no_ridge$coef[2]))
})

test_that("weightedLsRidge penalize_intercept=FALSE leaves last column unpenalized", {
  set.seed(42)
  n <- 100
  X <- cbind(rnorm(n), 1)
  y <- X %*% c(5, 3) + rnorm(n, sd = 0.5)
  w <- rep(1, n)
  # With large lambda but no intercept penalty, intercept should still be reasonable
  res <- pecotmr:::weightedLsRidge(y, X, w, lambda = 100, penalize_intercept = FALSE)
  # Intercept (col 2) should not be shrunk as aggressively as slope (col 1)
  res_pen <- pecotmr:::weightedLsRidge(y, X, w, lambda = 100, penalize_intercept = TRUE)
  # The intercept should differ between the two
  expect_false(isTRUE(all.equal(res$coef[2], res_pen$coef[2])))
})

# =============================================================================
# computeBaselineEnrichment
# =============================================================================

test_that("computeBaselineEnrichment returns correct structure", {
  M <- 100
  baseline_mat <- matrix(0, nrow = M, ncol = 2)
  baseline_mat[1:50, 1] <- 1
  baseline_mat[, 2] <- 1
  tau <- c(0.01, 0.005)
  tau_se <- c(0.002, 0.001)
  h2 <- 0.5
  annot_names <- c("annot1", "base")
  res <- pecotmr:::computeBaselineEnrichment(tau, tau_se, NULL,
                                              baseline_mat, annot_names, h2)
  expect_s3_class(res, "data.frame")
  expected_cols <- c("annotation", "tau", "tau_se", "enrichment",
                     "enrichment_se", "enrichment_p", "prop_h2", "prop_snps")
  expect_named(res, expected_cols)
  expect_equal(nrow(res), 2)
})

test_that("computeBaselineEnrichment computes enrichment = tau * M / h2", {
  M <- 100
  baseline_mat <- matrix(0, nrow = M, ncol = 2)
  baseline_mat[1:50, 1] <- 1
  baseline_mat[, 2] <- 1
  tau <- c(0.01, 0.005)
  tau_se <- c(0.002, 0.001)
  h2 <- 0.5
  annot_names <- c("annot1", "base")
  res <- pecotmr:::computeBaselineEnrichment(tau, tau_se, NULL,
                                              baseline_mat, annot_names, h2)
  expect_equal(res$enrichment, tau * M / h2)
})

test_that("computeBaselineEnrichment computes prop_snps correctly", {
  M <- 100
  baseline_mat <- matrix(0, nrow = M, ncol = 2)
  baseline_mat[1:50, 1] <- 1
  baseline_mat[, 2] <- 1
  tau <- c(0.01, 0.005)
  tau_se <- c(0.002, 0.001)
  h2 <- 0.5
  annot_names <- c("annot1", "base")
  res <- pecotmr:::computeBaselineEnrichment(tau, tau_se, NULL,
                                              baseline_mat, annot_names, h2)
  expect_equal(res$prop_snps, c(0.5, 1.0))
})

test_that("computeBaselineEnrichment uses jackknife blocks when provided", {
  M <- 100
  n_blocks <- 5
  baseline_mat <- matrix(0, nrow = M, ncol = 2)
  baseline_mat[1:50, 1] <- 1
  baseline_mat[, 2] <- 1
  tau <- c(0.01, 0.005)
  tau_se <- c(0.002, 0.001)
  h2 <- 0.5
  annot_names <- c("annot1", "base")
  tau_blocks <- matrix(rep(tau, each = n_blocks), nrow = n_blocks, ncol = 2)
  res <- pecotmr:::computeBaselineEnrichment(tau, tau_se, tau_blocks,
                                              baseline_mat, annot_names, h2)
  # With constant tau_blocks, enrichment_se should be 0
  expect_equal(res$enrichment_se, c(0, 0))
})

# =============================================================================
# standardize_tau_star
# =============================================================================

test_that("standardize_tau_star computes tau_star = tau * sd_annot * M_ref / h2g", {
  tau <- c(0.01, 0.02)
  sd_annot <- c(0.5, 0.3)
  M_ref <- 1000L
  h2g <- 0.4
  n_blocks <- 5
  tau_blocks <- matrix(rep(tau, each = n_blocks), nrow = n_blocks, ncol = 2)
  res <- pecotmr:::standardize_tau_star(tau, tau_blocks, sd_annot, M_ref, h2g)
  expected <- tau * sd_annot * M_ref / h2g
  expect_equal(res$tau_star, expected)
  expect_length(res$tau_star_se, 2)
})

test_that("standardize_tau_star errors when h2g == 0", {
  expect_error(
    pecotmr:::standardize_tau_star(c(0.01), matrix(0.01, nrow = 3, ncol = 1),
                                    c(0.5), 1000L, 0),
    "h2g must be non-zero"
  )
})

test_that("standardize_tau_star errors when tau and sd_annot differ in length", {
  expect_error(
    pecotmr:::standardize_tau_star(c(0.01, 0.02), matrix(0, nrow = 3, ncol = 2),
                                    c(0.5), 1000L, 0.5),
    "tau and sd_annot must have the same length"
  )
})

test_that("standardize_tau_star returns list with tau_star and tau_star_se", {
  tau <- c(0.01)
  sd_annot <- c(0.5)
  M_ref <- 1000L
  h2g <- 0.4
  tau_blocks <- matrix(rnorm(5, mean = 0.01, sd = 0.001), nrow = 5, ncol = 1)
  res <- pecotmr:::standardize_tau_star(tau, tau_blocks, sd_annot, M_ref, h2g)
  expect_true(is.list(res))
  expect_named(res, c("tau_star", "tau_star_se"))
  expect_length(res$tau_star, 1)
  expect_length(res$tau_star_se, 1)
})

# =============================================================================
# meta_random_effects
# =============================================================================

test_that("meta_random_effects returns all NA with k=0", {
  res <- pecotmr:::meta_random_effects(numeric(0), numeric(0))
  expect_true(is.na(res$mean))
  expect_true(is.na(res$se))
  expect_true(is.na(res$tau2))
  expect_true(is.na(res$I2))
  expect_true(is.na(res$Q))
})

test_that("meta_random_effects with k=1 returns input values", {
  res <- pecotmr:::meta_random_effects(5.0, 1.0)
  expect_equal(res$mean, 5.0)
  expect_equal(res$se, 1.0)
  expect_equal(res$tau2, 0)
})

test_that("meta_random_effects with identical means gives tau2=0", {
  means <- rep(3.0, 5)
  ses <- rep(1.0, 5)
  res <- pecotmr:::meta_random_effects(means, ses)
  expect_equal(res$tau2, 0)
  expect_equal(res$mean, 3.0)
})

test_that("meta_random_effects returns correct structure", {
  set.seed(42)
  means <- rnorm(5, mean = 2, sd = 0.5)
  ses <- rep(0.5, 5)
  res <- pecotmr:::meta_random_effects(means, ses)
  expect_true(is.list(res))
  expect_named(res, c("mean", "se", "tau2", "I2", "Q"))
  expect_true(res$se > 0)
  expect_true(res$I2 >= 0 && res$I2 <= 1)
  expect_true(res$Q >= 0)
  expect_true(res$tau2 >= 0)
})

test_that("meta_random_effects errors with non-positive ses", {
  expect_error(
    pecotmr:::meta_random_effects(c(1, 2), c(1, 0)),
    "all ses must be positive and finite"
  )
  expect_error(
    pecotmr:::meta_random_effects(c(1, 2), c(1, -1)),
    "all ses must be positive and finite"
  )
})

test_that("meta_random_effects known DerSimonian-Laird example", {
  # Three studies with known values
  means <- c(0.5, 0.8, 0.3)
  ses <- c(0.2, 0.3, 0.15)
  res <- pecotmr:::meta_random_effects(means, ses)

  # Fixed-effect weights
  w_fe <- 1 / ses^2
  mu_fe <- sum(w_fe * means) / sum(w_fe)
  Q <- sum(w_fe * (means - mu_fe)^2)
  c_dl <- sum(w_fe) - sum(w_fe^2) / sum(w_fe)
  tau2 <- max(0, (Q - 2) / c_dl)
  w_re <- 1 / (ses^2 + tau2)
  mu_re <- sum(w_re * means) / sum(w_re)
  se_re <- sqrt(1 / sum(w_re))

  expect_equal(res$Q, Q)
  expect_equal(res$tau2, tau2)
  expect_equal(res$mean, mu_re)
  expect_equal(res$se, se_re)
})

# =============================================================================
# snpsPerBlock
# =============================================================================

test_that("snpsPerBlock assigns SNPs to correct blocks", {
  # 10 SNPs across 2 blocks
  snp_info <- data.frame(
    CHR = rep("chr1", 10),
    BP = c(100, 500, 1000, 3000, 4500, 5500, 6000, 7500, 9000, 9999)
  )
  blocks_gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(1, 5001), end = c(5000, 10000))
  )
  ld_blocks <- new("LDBlocks", blocks = blocks_gr, genome = "hg19")
  res <- pecotmr:::snpsPerBlock(snp_info, ld_blocks)
  # Block 1 covers 1-5000: SNPs at 100, 500, 1000, 3000, 4500
  expect_equal(sort(res[["1"]]), c(1L, 2L, 3L, 4L, 5L))
  # Block 2 covers 5001-10000: SNPs at 5500, 6000, 7500, 9000, 9999
  expect_equal(sort(res[["2"]]), c(6L, 7L, 8L, 9L, 10L))
})

test_that("snpsPerBlock returns empty for blocks with no SNPs", {
  snp_info <- data.frame(
    CHR = rep("chr1", 3),
    BP = c(100, 200, 300)
  )
  blocks_gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr2"),
    ranges = IRanges::IRanges(start = c(1, 1), end = c(5000, 5000))
  )
  ld_blocks <- new("LDBlocks", blocks = blocks_gr, genome = "hg19")
  res <- pecotmr:::snpsPerBlock(snp_info, ld_blocks)
  # All SNPs on chr1, so block 2 (chr2) should have no entries
  expect_true("1" %in% names(res))
  expect_false("2" %in% names(res))
  expect_equal(sort(res[["1"]]), c(1L, 2L, 3L))
})

# =============================================================================
# checkGenomeBuild
# =============================================================================

test_that("checkGenomeBuild returns TRUE when all objects match", {
  blocks_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1, end = 5000)
  )
  ld_blocks <- new("LDBlocks", blocks = blocks_gr, genome = "hg19")

  set.seed(42)
  ss_df <- data.frame(
    SNP = paste0("rs", 1:5), CHR = "1", BP = 1:5,
    A1 = "A", A2 = "G", Z = rnorm(5), N = 1000
  )
  ss <- GWASSumStats(ss_df, genome = "hg19")

  expect_true(pecotmr:::checkGenomeBuild(ld_blocks, ss))
})

test_that("checkGenomeBuild errors when genome builds mismatch", {
  blocks_gr_19 <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1, end = 5000)
  )
  ld_blocks_19 <- new("LDBlocks", blocks = blocks_gr_19, genome = "hg19")

  blocks_gr_38 <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1, end = 5000)
  )
  ld_blocks_38 <- new("LDBlocks", blocks = blocks_gr_38, genome = "hg38")

  expect_error(
    pecotmr:::checkGenomeBuild(ld_blocks_19, ld_blocks_38),
    "Genome build mismatch"
  )
})

test_that("checkGenomeBuild works with AnnotationMatrix", {
  snp_gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 3),
    ranges = IRanges::IRanges(start = 1:3, width = 1)
  )
  annot_mat <- matrix(c(1, 0, 1), nrow = 3, ncol = 1)
  annot_meta <- data.frame(name = "annot1", tier = "baseline",
                           type = "binary", stringsAsFactors = FALSE)
  am <- new("AnnotationMatrix",
    snp_ranges = snp_gr,
    annotations = annot_mat,
    annotation_meta = annot_meta,
    genome = "hg19"
  )

  blocks_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1, end = 5000)
  )
  ld_blocks <- new("LDBlocks", blocks = blocks_gr, genome = "hg19")

  expect_true(pecotmr:::checkGenomeBuild(am, ld_blocks))
})

# =============================================================================
# shrinkLd
# =============================================================================

test_that("shrinkLd constant shrinkage applies (1-lambda)*R + lambda*I", {
  set.seed(42)
  p <- 5
  # Create a valid correlation matrix
  A <- matrix(rnorm(p * p), p, p)
  R <- cov2cor(crossprod(A))
  n_ref <- 100

  res <- pecotmr:::shrinkLd(R, n_ref, shrinkage_type = "constant")
  lambda <- 1 / sqrt(n_ref)
  expected <- (1 - lambda) * R + lambda * diag(p)
  expect_equal(res, expected)
})

test_that("shrinkLd constant shrinkage preserves diagonal of 1", {
  set.seed(42)
  p <- 4
  A <- matrix(rnorm(p * p), p, p)
  R <- cov2cor(crossprod(A))
  n_ref <- 200
  res <- pecotmr:::shrinkLd(R, n_ref, shrinkage_type = "constant")
  expect_equal(diag(res), rep(1, p))
})

test_that("shrinkLd constant shrinkage result is symmetric", {
  set.seed(42)
  p <- 4
  A <- matrix(rnorm(p * p), p, p)
  R <- cov2cor(crossprod(A))
  n_ref <- 200
  res <- pecotmr:::shrinkLd(R, n_ref, shrinkage_type = "constant")
  expect_equal(res, t(res))
})

test_that("shrinkLd wen_stephens uses genetic map when provided", {
  set.seed(42)
  p <- 4
  A <- matrix(rnorm(p * p), p, p)
  R <- cov2cor(crossprod(A))
  n_ref <- 500
  genetic_map <- c(0.0, 0.1, 0.5, 1.0)
  res <- pecotmr:::shrinkLd(R, n_ref, shrinkage_type = "wen_stephens",
                             genetic_map = genetic_map)
  # Diagonal should be 1
  expect_equal(diag(res), rep(1, p))
  # Result should be symmetric
  expect_equal(res, t(res))
  # Off-diagonal elements should be shrunk (closer to zero than original)
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      expect_true(abs(res[i, j]) <= abs(R[i, j]) + 1e-10)
    }
  }
})
