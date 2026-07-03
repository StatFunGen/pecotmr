context("regularized_regression — penalizedRss")

# ---- penalizedRss (low-level solver) ----

test_that("penalizedRss errors on non-matrix R input", {
  expect_error(penalizedRss(bhat = rnorm(5), R = "not_a_matrix", n = 100,
                             penalty = "MCP"),
               "as a matrix")
})

test_that("penalizedRss errors on non-positive sample size", {
  expect_error(penalizedRss(bhat = rnorm(5), R = diag(5),
                             n = -1, penalty = "SCAD"),
               "valid sample size")
})

test_that("penalizedRss errors on mismatched bhat and R dimensions", {
  expect_error(
    penalizedRss(bhat = rnorm(10), R = diag(5), n = 100,
                  penalty = "MCP"),
    "number of rows of 'R'"
  )
})

test_that("penalizedRss with large lambda gives all-zero betas (MCP)", {
  set.seed(42)
  p <- 10; n <- 100
  bhat <- rnorm(p, sd = 0.1)
  result <- penalizedRss(bhat = bhat, R = diag(p), n = n,
                          penalty = "MCP", lambda = c(100))
  expect_true(all(result$betaEst == 0))
})

test_that("penalizedRss with large lambda gives all-zero betas (SCAD)", {
  set.seed(42)
  p <- 10; n <- 100
  bhat <- rnorm(p, sd = 0.1)
  result <- penalizedRss(bhat = bhat, R = diag(p), n = n,
                          penalty = "SCAD", lambda = c(100))
  expect_true(all(result$betaEst == 0))
})

test_that("penalizedRss with large lambda0 gives all-zero betas (L0)", {
  set.seed(42)
  p <- 10; n <- 100
  bhat <- rnorm(p, sd = 0.1)
  result <- penalizedRss(bhat = bhat, R = diag(p), n = n,
                          penalty = "L0", lambda = c(0), lambda0 = 1e6)
  expect_true(all(result$betaEst == 0))
})

test_that("penalizedRss runs with MCP and returns correct structure", {
  set.seed(42)
  p <- 10; n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  for (i in 1:(p - 1)) {
    R[i, i + 1] <- 0.3
    R[i + 1, i] <- 0.3
  }
  result <- penalizedRss(bhat = bhat, R = R, n = n,
                          penalty = "MCP")
  expect_type(result, "list")
  expect_true("betaEst" %in% names(result))
  expect_equal(length(result$betaEst), p)
  expect_true(all(is.finite(result$betaEst)))
  expect_equal(nrow(result$beta), p)
  expect_equal(ncol(result$beta), 20)
  expect_true(all(result$conv %in% c(0L, 1L)))
})

test_that("penalizedRss runs with SCAD and returns correct structure", {
  set.seed(42)
  p <- 10; n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  result <- penalizedRss(bhat = bhat, R = R, n = n,
                          penalty = "SCAD")
  expect_type(result, "list")
  expect_equal(length(result$betaEst), p)
  expect_true(all(is.finite(result$betaEst)))
})

test_that("penalizedRss runs with L0 and returns correct structure", {
  set.seed(42)
  p <- 10; n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  result <- penalizedRss(bhat = bhat, R = R, n = n,
                          penalty = "L0", lambda = c(0), lambda0 = 0.01)
  expect_type(result, "list")
  expect_equal(length(result$betaEst), p)
  expect_true(all(is.finite(result$betaEst)))
})

test_that("penalizedRss LASSO matches lassosumRss on identity LD", {
  set.seed(42)
  p <- 10; n <- 200
  bhat <- rnorm(p, sd = 0.2)
  R <- diag(p)
  lam <- exp(seq(log(0.001), log(0.1), length.out = 5))
  res_lasso <- penalizedRss(bhat = bhat, R = R, n = n,
                             penalty = "lasso", lambda = lam)
  res_lassosum <- lassosumRss(bhat = bhat, R = R, n = n,
                               lambda = lam)
  expect_equal(res_lasso$beta, res_lassosum$beta, tolerance = 1e-4)
})

# ---- Signal recovery on simulated summary statistics ----

test_that("penalizedRss MCP recovers signal direction", {
  set.seed(2024)
  n <- 500; p <- 20
  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n)
  beta_true <- rep(0, p)
  beta_true[c(3, 10, 15)] <- c(0.4, -0.3, 0.2)
  y <- as.numeric(X %*% beta_true) + rnorm(n)
  bhat <- as.vector(cor(y, X))
  R <- cor(X)
  stat <- list(b = bhat, n = rep(n, p))
  w <- mcpRssWeights(stat = stat, LD = R, s = c(0.5, 0.9))
  expect_equal(length(w), p)
  expect_true(all(is.finite(w)))
  expect_gt(cor(w, beta_true), 0.4)
})

test_that("penalizedRss SCAD recovers signal direction", {
  set.seed(2024)
  n <- 500; p <- 20
  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n)
  beta_true <- rep(0, p)
  beta_true[c(3, 10, 15)] <- c(0.4, -0.3, 0.2)
  y <- as.numeric(X %*% beta_true) + rnorm(n)
  bhat <- as.vector(cor(y, X))
  R <- cor(X)
  stat <- list(b = bhat, n = rep(n, p))
  w <- scadRssWeights(stat = stat, LD = R, s = c(0.5, 0.9))
  expect_equal(length(w), p)
  expect_true(all(is.finite(w)))
  expect_gt(cor(w, beta_true), 0.4)
})

test_that("penalizedRss L0 recovers signal direction", {
  set.seed(2024)
  n <- 500; p <- 20
  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n)
  beta_true <- rep(0, p)
  beta_true[c(3, 10, 15)] <- c(0.4, -0.3, 0.2)
  y <- as.numeric(X %*% beta_true) + rnorm(n)
  bhat <- as.vector(cor(y, X))
  R <- cor(X)
  stat <- list(b = bhat, n = rep(n, p))
  w <- l0learnRssWeights(stat = stat, LD = R, penalty = "L0",
                           s = c(0.5, 0.9),
                           lambda0 = exp(seq(log(0.001), log(0.5), length.out = 5)))
  expect_equal(length(w), p)
  expect_true(all(is.finite(w)))
  expect_gt(cor(w, beta_true), 0.3)
})

# ---- Weight wrapper structure ----

test_that("scadRssWeights returns correct-length vector with attributes", {
  set.seed(42)
  p <- 10; n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  stat <- list(b = bhat, n = rep(n, p))
  w <- scadRssWeights(stat = stat, LD = R, s = 0.5)
  expect_equal(length(w), p)
  expect_true(is.numeric(w))
  sel <- attr(w, "penalized_rss_selection")
  expect_true(!is.null(sel))
  expect_equal(unname(sel["penalty"]), "SCAD")
})

test_that("mcpRssWeights returns correct-length vector with attributes", {
  set.seed(42)
  p <- 10; n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  stat <- list(b = bhat, n = rep(n, p))
  w <- mcpRssWeights(stat = stat, LD = R, s = 0.5)
  expect_equal(length(w), p)
  expect_true(is.numeric(w))
  sel <- attr(w, "penalized_rss_selection")
  expect_true(!is.null(sel))
  expect_equal(unname(sel["penalty"]), "MCP")
})

test_that("l0learnRssWeights returns correct-length vector with attributes", {
  set.seed(42)
  p <- 10; n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  stat <- list(b = bhat, n = rep(n, p))
  w <- l0learnRssWeights(stat = stat, LD = R, penalty = "L0",
                           s = 0.5, lambda0 = c(0.01, 0.1))
  expect_equal(length(w), p)
  expect_true(is.numeric(w))
  sel <- attr(w, "penalized_rss_selection")
  expect_true(!is.null(sel))
  expect_equal(unname(sel["penalty"]), "L0")
})
