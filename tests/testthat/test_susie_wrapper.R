context("susie_wrapper")

# =============================================================================
# lbf_to_alpha_vector (internal)
# =============================================================================

test_that("lbf_to_alpha_vector converts correctly", {
  lbf <- c(a = -0.5, b = 1.2, c = 0.3)
  alpha <- pecotmr:::lbf_to_alpha_vector(lbf)
  expect_length(alpha, 3)
  expect_equal(sum(alpha), 1, tolerance = 1e-10)
  expect_true(alpha["b"] > alpha["a"])
})

test_that("lbf_to_alpha_vector with prior weights", {
  lbf <- c(a = 1, b = 1, c = 1)  # Equal LBFs
  pw <- c(0.5, 0.25, 0.25)
  alpha <- pecotmr:::lbf_to_alpha_vector(lbf, prior_weights = pw)
  expect_true(alpha[1] > alpha[2])
  expect_equal(sum(alpha), 1, tolerance = 1e-10)
})

test_that("lbf_to_alpha_vector returns zeros for all-zero lbf", {
  lbf <- c(a = 0, b = 0, c = 0)
  alpha <- pecotmr:::lbf_to_alpha_vector(lbf)
  expect_true(all(alpha == 0))
})

test_that("lbf_to_alpha_vector handles single element", {
  lbf <- c(a = 2.0)
  alpha <- pecotmr:::lbf_to_alpha_vector(lbf)
  expect_length(alpha, 1)
  expect_equal(alpha[["a"]], 1.0)
})

test_that("lbf_to_alpha_vector handles very large LBFs without overflow", {
  lbf <- c(a = 500, b = 500.1, c = 499)
  alpha <- pecotmr:::lbf_to_alpha_vector(lbf)
  expect_true(all(is.finite(alpha)))
  expect_equal(sum(alpha), 1, tolerance = 1e-10)
  expect_true(alpha["b"] > alpha["a"])
})

test_that("lbf_to_alpha_vector handles very negative LBFs", {
  lbf <- c(a = -1000, b = -999, c = -1001)
  alpha <- pecotmr:::lbf_to_alpha_vector(lbf)
  expect_true(all(is.finite(alpha)))
  expect_equal(sum(alpha), 1, tolerance = 1e-10)
  expect_true(alpha["b"] > alpha["a"])
})

test_that("lbf_to_alpha_vector with unequal prior weights", {
  lbf <- c(a = 0.5, b = 0.5, c = 0.5)
  pw <- c(0.8, 0.1, 0.1)
  alpha <- pecotmr:::lbf_to_alpha_vector(lbf, prior_weights = pw)
  expect_equal(sum(alpha), 1, tolerance = 1e-10)
  expect_true(alpha[1] > 0.7)
})

# =============================================================================
# lbf_to_alpha (matrix version)
# =============================================================================

test_that("lbf_to_alpha converts log BFs to posteriors", {
  lbf <- matrix(c(0, 3, 2, 1, 4, 0), nrow = 2, ncol = 3)
  alpha <- pecotmr:::lbf_to_alpha(lbf)
  expect_equal(dim(alpha), c(2, 3))
  expect_equal(rowSums(alpha), c(1, 1), tolerance = 1e-10)
  expect_true(alpha[1, 3] > alpha[1, 1])
  expect_true(alpha[2, 1] > alpha[2, 3])
})

test_that("lbf_to_alpha handles uniform lbf", {
  lbf <- matrix(1, nrow = 1, ncol = 5)
  alpha <- pecotmr:::lbf_to_alpha(lbf)
  expect_equal(as.numeric(alpha), rep(0.2, 5), tolerance = 1e-10)
})

test_that("lbf_to_alpha handles single-row matrix", {
  lbf <- matrix(c(1.0, 2.0, 0.5), nrow = 1)
  colnames(lbf) <- paste0("v", 1:3)
  result <- lbf_to_alpha(lbf)
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 3)
  expect_equal(sum(result), 1, tolerance = 1e-10)
})

test_that("lbf_to_alpha handles large matrix", {
  set.seed(42)
  lbf <- matrix(rnorm(100), nrow = 10, ncol = 10)
  colnames(lbf) <- paste0("v", 1:10)
  result <- lbf_to_alpha(lbf)
  expect_equal(dim(result), c(10, 10))
  expect_equal(rowSums(result), rep(1, 10), tolerance = 1e-10)
})

test_that("lbf_to_alpha with mixed zero and nonzero rows", {
  lbf <- matrix(c(0, 0, 0, 1, 2, 3), nrow = 2, byrow = TRUE)
  colnames(lbf) <- paste0("v", 1:3)
  result <- lbf_to_alpha(lbf)
  expect_true(all(result[1, ] == 0))
  expect_equal(sum(result[2, ]), 1, tolerance = 1e-10)
})

# =============================================================================
# get_cs_index (internal)
# =============================================================================

test_that("get_cs_index finds variant in credible set", {
  susie_cs <- list(L1 = c(1, 2, 3), L2 = c(4, 5))
  idx <- pecotmr:::get_cs_index(2, susie_cs)
  expect_equal(unname(idx), 1)
})

test_that("get_cs_index returns NA for variant not in any CS", {
  susie_cs <- list(L1 = c(1, 2), L2 = c(4, 5))
  idx <- pecotmr:::get_cs_index(3, susie_cs)
  expect_true(is.na(idx))
})

test_that("get_cs_index warns on variant in multiple CS", {
  susie_cs <- list(L1 = c(1, 2, 3), L2 = c(2, 4, 5))
  expect_warning(pecotmr:::get_cs_index(2, susie_cs), "found in multiple CS")
})

test_that("get_cs_index returns smallest CS when variant in multiple", {
  susie_cs <- list(L1 = c(1, 2, 3, 4, 5), L2 = c(2, 3))
  result <- suppressWarnings(pecotmr:::get_cs_index(2, susie_cs))
  expect_equal(unname(result), 2)
})

test_that("get_cs_index handles empty CS list", {
  susie_cs <- list()
  result <- pecotmr:::get_cs_index(1, susie_cs)
  expect_true(is.na(result))
})

test_that("get_cs_index returns correct CS assignment with real susie fit", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 200
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(2, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  y <- X %*% beta + rnorm(n, sd = 0.5)
  fit <- susieR::susie(X, y, L = 5)
  if (!is.null(fit$sets$cs)) {
    idx <- pecotmr:::get_cs_index(fit$sets$cs, 1)
    expect_true(is.numeric(unname(idx)))
  }
})

# =============================================================================
# get_top_variants_idx (internal)
# =============================================================================

test_that("get_top_variants_idx returns combined PIP and CS variants", {
  susie_output <- list(
    pip = c(0.01, 0.15, 0.02, 0.5, 0.01),
    sets = list(cs = list(L1 = c(1, 2)))
  )
  result <- pecotmr:::get_top_variants_idx(susie_output, signal_cutoff = 0.1)
  expect_true(1 %in% result)
  expect_true(2 %in% result)
  expect_true(4 %in% result)
  expect_true(all(result == sort(result)))
})

test_that("get_top_variants_idx with no CS", {
  susie_output <- list(
    pip = c(0.01, 0.5, 0.02, 0.8, 0.01),
    sets = list(cs = NULL)
  )
  result <- pecotmr:::get_top_variants_idx(susie_output, signal_cutoff = 0.1)
  expect_equal(result, c(2, 4))
})

test_that("get_top_variants_idx with all low PIPs", {
  susie_output <- list(
    pip = c(0.01, 0.02, 0.03),
    sets = list(cs = list(L1 = c(1, 2)))
  )
  result <- pecotmr:::get_top_variants_idx(susie_output, signal_cutoff = 0.5)
  expect_equal(result, c(1, 2))
})

test_that("get_top_variants_idx with high cutoff and no CS", {
  susie_output <- list(
    pip = c(0.01, 0.02, 0.03),
    sets = list(cs = NULL)
  )
  result <- pecotmr:::get_top_variants_idx(susie_output, signal_cutoff = 0.5)
  expect_length(result, 0)
})

# =============================================================================
# get_cs_info (internal)
# =============================================================================

test_that("get_cs_info maps variants to CS numbers", {
  susie_cs <- list(L1 = c(1, 2), L3 = c(4, 5, 6))
  top_idx <- c(1, 3, 5)
  result <- pecotmr:::get_cs_info(susie_cs, top_idx)
  # Now returns data.frame(variant_idx, cs_idx) with one row per (variant, CS) pair
  expect_true(is.data.frame(result))
  expect_equal(result$variant_idx, c(1, 3, 5))
  expect_equal(result$cs_idx, c(1L, 0L, 3L))
})

test_that("get_cs_info handles all variants outside CS", {
  susie_cs <- list(L1 = c(1, 2))
  top_idx <- c(5, 6, 7)
  result <- pecotmr:::get_cs_info(susie_cs, top_idx)
  expect_true(is.data.frame(result))
  expect_true(all(result$cs_idx == 0))
})

test_that("get_cs_info reports variant in multiple CSs as multiple rows", {
  susie_cs <- list(L1 = c(1, 2, 3), L3 = c(2, 3, 4))
  top_idx <- c(1, 2, 4)
  result <- pecotmr:::get_cs_info(susie_cs, top_idx)
  expect_true(is.data.frame(result))
  # variant 2 is in both L1 and L3, so it gets two rows
  expect_equal(nrow(result), 4)
  expect_equal(sum(result$variant_idx == 2), 2)
  expect_equal(sort(result$cs_idx[result$variant_idx == 2]), c(1L, 3L))
})

# =============================================================================
# susie_rss_pipeline
# =============================================================================

test_that("susie_rss_pipeline errors on missing z and beta/se", {
  sumstats <- data.frame(x = 1)
  LD_mat <- matrix(1)
  expect_error(susie_rss_pipeline(sumstats, LD_mat), "should have 'z'")
})

test_that("susie_rss_pipeline errors on invalid method", {
  sumstats <- list(z = rnorm(5))
  LD_mat <- diag(5)
  expect_error(susie_rss_pipeline(sumstats, LD_mat, analysis_method = "invalid"))
})

test_that("susie_rss_pipeline runs with single_effect method", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 20
  z <- rnorm(n)
  names(z) <- paste0("chr1:", seq_len(n), ":A:G")
  R <- diag(n)
  colnames(R) <- rownames(R) <- names(z)
  sumstats <- list(z = z)

  result <- susie_rss_pipeline(sumstats, R, analysis_method = "single_effect")
  expect_true(is.list(result))
  expect_true("variant_names" %in% names(result))
  expect_true("susie_result_trimmed" %in% names(result))
  # PIPs should be numeric, in [0,1], and sum to at most 1 (L=1)
  pip <- result$susie_result_trimmed$pip
  expect_true(is.numeric(pip))
  expect_length(pip, n)
  expect_true(all(pip >= 0 & pip <= 1))
  expect_true(sum(pip) <= 1 + 1e-6)
  # Credible sets, if any, should contain valid indices
  cs_list <- result$susie_result_trimmed$sets$cs
  if (!is.null(cs_list)) {
    for (cs in cs_list) {
      expect_true(all(cs >= 1 & cs <= n))
    }
  }
})

test_that("susie_rss_pipeline runs with bayesian_conditional_regression", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 20
  z <- rnorm(n)
  names(z) <- paste0("chr1:", seq_len(n), ":A:G")
  R <- diag(n)
  colnames(R) <- rownames(R) <- names(z)
  sumstats <- list(z = z)

  result <- susie_rss_pipeline(sumstats, R,
    analysis_method = "bayesian_conditional_regression",
    L = 5, max_L = 5
  )
  expect_true(is.list(result))
  expect_true("susie_result_trimmed" %in% names(result))
  pip <- result$susie_result_trimmed$pip
  expect_true(is.numeric(pip))
  expect_length(pip, n)
  expect_true(all(pip >= 0 & pip <= 1))
  # With L=5, sum of PIPs can be up to L
  expect_true(sum(pip) <= 5 + 1e-6)
  cs_list <- result$susie_result_trimmed$sets$cs
  if (!is.null(cs_list)) {
    for (cs in cs_list) {
      expect_true(all(cs >= 1 & cs <= n))
    }
  }
})

test_that("susie_rss_pipeline uses beta/se when z not provided", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 15
  beta <- rnorm(n, sd = 0.1)
  se <- rep(0.1, n)
  names(beta) <- paste0("chr1:", seq_len(n), ":A:G")
  R <- diag(n)
  colnames(R) <- rownames(R) <- names(beta)
  sumstats <- list(beta = beta, se = se)

  result <- susie_rss_pipeline(sumstats, R,
    analysis_method = "susie_rss",
    L = 5, max_L = 5
  )
  expect_true(is.list(result))
  expect_true("susie_result_trimmed" %in% names(result))
  pip <- result$susie_result_trimmed$pip
  expect_true(is.numeric(pip))
  expect_length(pip, n)
  expect_true(all(pip >= 0 & pip <= 1))
  expect_true(sum(pip) <= 5 + 1e-6)
  cs_list <- result$susie_result_trimmed$sets$cs
  if (!is.null(cs_list)) {
    for (cs in cs_list) {
      expect_true(all(cs >= 1 & cs <= n))
    }
  }
})

# =============================================================================
# susie_rss_wrapper
# =============================================================================

test_that("susie_rss_wrapper with L=1 runs single effect", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 10
  R <- diag(p)
  z <- rnorm(p)
  result <- susie_rss_wrapper(z = z, R = R, L = 1)
  expect_true("pip" %in% names(result))
  expect_length(result$pip, p)
  expect_true(is.numeric(result$pip))
  expect_true(all(result$pip >= 0 & result$pip <= 1))
  # L=1 so PIPs sum to at most 1
  expect_true(sum(result$pip) <= 1 + 1e-6)
  if (!is.null(result$sets$cs)) {
    for (cs in result$sets$cs) {
      expect_true(all(cs >= 1 & cs <= p))
    }
  }
})

test_that("susie_rss_wrapper with L equal to max_L", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 15
  z <- rnorm(n)
  R <- diag(n)
  result <- susie_rss_wrapper(z = z, R = R, L = 5, max_L = 5)
  expect_true("pip" %in% names(result))
})

test_that("susie_rss_wrapper dynamic L with no CS found", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 10
  R <- diag(p)
  z <- rep(0.1, p)
  result <- susie_rss_wrapper(z = z, R = R, L = 2, max_L = 10, l_step = 2)
  expect_true("pip" %in% names(result))
})

# =============================================================================
# susie_wrapper
# =============================================================================

test_that("susie_wrapper runs with init_L equal to max_L", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 50
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] * 2 + rnorm(n)

  result <- susie_wrapper(X, y, init_L = 5, max_L = 5)
  expect_true("pip" %in% names(result))
  expect_length(result$pip, p)
})

test_that("susie_wrapper dynamically adjusts L", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 200
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] + rnorm(n)
  result <- susie_wrapper(X, y, init_L = 1, max_L = 10, l_step = 2)
  expect_true("pip" %in% names(result))
  expect_true(!is.null(result$sets))
})

# =============================================================================
# susie_weights
# =============================================================================

test_that("susie_weights returns zeros when fit lacks alpha/mu", {
  fake_fit <- list(pip = rep(0.01, 5))
  result <- susie_weights(susie_fit = fake_fit)
  expect_equal(result, rep(0, 5))
})

test_that("susie_weights checks dimension mismatch", {
  set.seed(42)
  X <- matrix(rnorm(100), 20, 5)
  fake_fit <- list(pip = rep(0.01, 10))
  expect_error(susie_weights(X = X, susie_fit = fake_fit), "Dimension mismatch")
})

# =============================================================================
# susie_ash_weights
# =============================================================================

test_that("susie_ash_weights returns zeros without proper fit structure", {
  fake_fit <- list(pip = rep(0.01, 5))
  result <- susie_ash_weights(susie_ash_fit = fake_fit)
  expect_equal(result, rep(0, 5))
})

# =============================================================================
# susie_inf_weights
# =============================================================================

test_that("susie_inf_weights returns zeros without proper fit structure", {
  fake_fit <- list(pip = rep(0.01, 5))
  result <- susie_inf_weights(susie_inf_fit = fake_fit)
  expect_equal(result, rep(0, 5))
})

# =============================================================================
# glmnet_weights
# =============================================================================

test_that("glmnet_weights produces non-zero weights for correlated data", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 100
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  beta_true <- c(3, -2, rep(0, p - 2))
  y <- X %*% beta_true + rnorm(n)

  w <- glmnet_weights(X, y, alpha = 0.5)
  expect_length(w, p)
  expect_true(any(w != 0))
})

test_that("glmnet_weights handles zero-variance columns", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 100
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  X[, 3] <- 1  # zero variance column
  y <- X[, 1] * 2 + rnorm(n)

  w <- glmnet_weights(X, y, alpha = 1)
  expect_length(w, p)
  expect_equal(w[3], 0)
})

# =============================================================================
# init_prior_sd
# =============================================================================

test_that("init_prior_sd returns n standard deviations", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n_samples <- 50
  p <- 10
  X <- matrix(rnorm(n_samples * p), n_samples, p)
  y <- X[, 1] * 2 + rnorm(n_samples)

  sds <- pecotmr:::init_prior_sd(X, y, n = 15)
  expect_length(sds, 15)
  expect_equal(sds[1], 0)
  expect_true(all(diff(sds) >= 0))
})
