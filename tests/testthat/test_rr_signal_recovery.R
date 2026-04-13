context("regularized_regression — signal recovery")

# ============================================================================
# Signal recovery tests for X-y wrappers
#
# These tests catch breaking regressions in the post-solver coefficient
# extraction code (which the dispatch/end-to-end mocks short-circuit). They
# simulate a small sparse linear model and assert that the returned weights
# correlate with truth (or that the top-K nonzero indices overlap truth).
# ============================================================================

# Helper: simulate a small sparse linear model with binomial genotypes.
.simulate_sparse_xy <- function(n = 200, p = 20,
                                signal_idx = c(3, 10, 15),
                                signal_beta = c(0.8, -0.6, 0.5),
                                seed = 2024) {
  set.seed(seed)
  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n)
  beta_true <- rep(0, p)
  beta_true[signal_idx] <- signal_beta
  y <- as.numeric(X %*% beta_true) + rnorm(n, sd = 0.5)
  list(X = X, y = y, beta_true = beta_true, signal_idx = signal_idx)
}

test_that("lasso_weights recovers signal direction on simulated data", {
  skip_if_not_installed("glmnet")
  sim <- .simulate_sparse_xy()
  w <- lasso_weights(sim$X, sim$y)
  expect_equal(length(w), ncol(sim$X))
  expect_true(all(is.finite(w)))
  expect_gt(cor(w, sim$beta_true), 0.5)
})

test_that("scad_weights recovers signal indices on simulated data", {
  skip_if_not_installed("ncvreg")
  sim <- .simulate_sparse_xy()
  w <- scad_weights(sim$X, sim$y, nfolds = 5)
  expect_equal(length(w), ncol(sim$X))
  expect_true(all(is.finite(w)))
  # SCAD is sparse: top-3 by absolute weight should heavily overlap truth.
  top3 <- order(abs(w), decreasing = TRUE)[1:3]
  expect_gte(length(intersect(top3, sim$signal_idx)), 2L)
})

test_that("l0learn_weights with default L0 penalty recovers signal indices", {
  skip_if_not_installed("L0Learn")
  sim <- .simulate_sparse_xy()
  w <- l0learn_weights(sim$X, sim$y, penalty = "L0", nFolds = 5)
  expect_equal(length(w), ncol(sim$X))
  expect_true(all(is.finite(w)))
  top3 <- order(abs(w), decreasing = TRUE)[1:3]
  expect_gte(length(intersect(top3, sim$signal_idx)), 2L)
})

test_that("l0learn_weights with L0L2 penalty recovers signal direction", {
  skip_if_not_installed("L0Learn")
  sim <- .simulate_sparse_xy()
  w <- l0learn_weights(sim$X, sim$y, penalty = "L0L2", nFolds = 5)
  expect_equal(length(w), ncol(sim$X))
  expect_true(all(is.finite(w)))
  # L0L2 selects (gamma, lambda) over a 2D grid; correlation is the safer bet.
  expect_gt(cor(w, sim$beta_true), 0.4)
})

test_that("dpr_weights with VB fitting recovers signal direction", {
  skip_if_not_installed("RcppDPR")
  sim <- .simulate_sparse_xy()
  w <- dpr_weights(sim$X, sim$y, fitting_method = "VB")
  expect_equal(length(w), ncol(sim$X))
  expect_true(all(is.finite(w)))
  expect_gt(cor(w, sim$beta_true), 0.4)
})

test_that("lassosum_rss_weights recovers signal direction on simulated data", {
  skip_if_not_installed("Rcpp")
  set.seed(2024)
  n <- 500
  p <- 20
  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n)
  beta_true <- rep(0, p)
  beta_true[c(3, 10, 15)] <- c(0.4, -0.3, 0.2)
  y <- as.numeric(X %*% beta_true) + rnorm(n)
  bhat <- as.vector(cor(y, X))
  R <- cor(X)
  # Heterogeneous n forces the median(stat$n) path; mean would be wrong.
  stat <- list(b = bhat, n = c(rep(n, p - 1), 10 * n))
  w <- lassosum_rss_weights(stat = stat, LD = R, s = c(0.2, 0.5, 0.9))
  expect_equal(length(w), p)
  expect_true(all(is.finite(w)))
  expect_gt(cor(w, beta_true), 0.5)
})
