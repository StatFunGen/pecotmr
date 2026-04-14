context("regularized_regression — sdpr")

# ---- sdpr ----
test_that("sdpr errors on mismatched bhat and LD dimensions", {
  expect_error(
    sdpr(bhat = rnorm(10), LD = list(blk1 = diag(5)), n = 100),
    "same as the length of bhat"
  )
})

test_that("sdpr errors on non-positive sample size", {
  expect_error(
    sdpr(bhat = rnorm(5), LD = list(blk1 = diag(5)), n = -1),
    "positive integer"
  )
})

test_that("sdpr errors when M is less than 4", {
  expect_error(
    sdpr(bhat = rnorm(5), LD = list(blk1 = diag(5)), n = 100, M = 3),
    "'M' must be at least 4"
  )
})

test_that("sdpr errors on invalid per_variant_sample_size", {
  expect_error(
    sdpr(bhat = rnorm(5), LD = list(blk1 = diag(5)), n = 100,
         per_variant_sample_size = c(100, -1, 100, 100, 100)),
    "positive values"
  )
})

test_that("sdpr errors on invalid array values", {
  expect_error(
    sdpr(bhat = rnorm(5), LD = list(blk1 = diag(5)), n = 100,
         array = c(0, 1, 3, 1, 0)),
    "0, 1, or 2"
  )
})

test_that("sdpr runs successfully", {
  set.seed(42)
  p <- 10
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  result <- sdpr(bhat = bhat, LD = list(blk1 = R), n = 100,
                 iter = 50, burn = 10, thin = 2, verbose = FALSE)
  expect_type(result, "list")
  expect_true("beta_est" %in% names(result))
  expect_equal(length(result$beta_est), p)
  expect_true(all(is.finite(result$beta_est)))
})

test_that("sdpr with per_variant_sample_size", {
  set.seed(42)
  p <- 10
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  result <- sdpr(bhat = bhat, LD = list(blk1 = R), n = 100,
                 per_variant_sample_size = rep(100, p),
                 iter = 50, burn = 10, thin = 2, verbose = FALSE)
  expect_equal(length(result$beta_est), p)
  expect_true(all(is.finite(result$beta_est)))
})

test_that("sdpr with valid array parameter", {
  set.seed(42)
  p <- 10
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  result <- sdpr(bhat = bhat, LD = list(blk1 = R), n = 100,
                 array = rep(1, p),
                 iter = 50, burn = 10, thin = 2, verbose = FALSE)
  expect_equal(length(result$beta_est), p)
  expect_true(all(is.finite(result$beta_est)))
})

# ---- sdpr signal recovery ----
test_that("sdpr recovers signal direction on simulated genotype data", {
  set.seed(2024)
  n <- 500
  p <- 20
  # Realistic genotype matrix with non-trivial LD (binomial, MAF=0.3)
  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n)
  beta_true <- rep(0, p)
  beta_true[c(3, 10, 15)] <- c(0.4, -0.3, 0.2)
  y <- X %*% beta_true + rnorm(n)
  bhat <- as.vector(cor(y, X))
  R <- cor(X)
  result <- sdpr(bhat = bhat, LD = list(blk1 = R), n = n,
                 iter = 500, burn = 200, thin = 5, verbose = FALSE, seed = 42L)
  expect_true("beta_est" %in% names(result))
  expect_equal(length(result$beta_est), p)
  expect_true(all(is.finite(result$beta_est)))
  # Correlation with truth should be positive (signal recovery)
  expect_gt(cor(result$beta_est, beta_true), 0.3)
})

test_that("sdpr accepts multiple LD blocks with realistic genotype data", {
  set.seed(2024)
  n <- 500
  p <- 20
  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n)
  beta_true <- rep(0, p)
  beta_true[c(3, 15)] <- c(0.4, 0.2)
  y <- X %*% beta_true + rnorm(n)
  bhat <- as.vector(cor(y, X))
  R <- cor(X)
  R1 <- R[1:10, 1:10]
  R2 <- R[11:20, 11:20]
  result <- sdpr(bhat = bhat, LD = list(blk1 = R1, blk2 = R2), n = n,
                 iter = 500, burn = 200, thin = 5, verbose = FALSE, seed = 42L)
  expect_equal(length(result$beta_est), p)
  expect_true(all(is.finite(result$beta_est)))
})

# ---- sdpr verbose output ----
test_that("sdpr with verbose = TRUE produces output", {
  set.seed(42)
  p <- 10
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  # iter >= 100 triggers the verbose print inside the MCMC loop
  result <- sdpr(bhat = bhat, LD = list(blk1 = R), n = 100,
                 iter = 110, burn = 10, thin = 2, verbose = TRUE, seed = 42L)
  expect_type(result, "list")
  expect_equal(length(result$beta_est), p)
  expect_true(all(is.finite(result$beta_est)))
})

# ---- sdpr opt_llk = 2 (multi-array) ----
test_that("sdpr runs with opt_llk = 2 and mixed array values", {
  set.seed(42)
  p <- 10
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  # Mixed array: some variants from array 1, some from array 2
  arr <- c(rep(1, 5), rep(2, 5))
  result <- sdpr(bhat = bhat, LD = list(blk1 = R), n = 100,
                 per_variant_sample_size = rep(100, p),
                 array = arr, opt_llk = 2,
                 iter = 50, burn = 10, thin = 2, verbose = FALSE, seed = 42L)
  expect_type(result, "list")
  expect_equal(length(result$beta_est), p)
  expect_true(all(is.finite(result$beta_est)))
})

test_that("sdpr opt_llk = 2 with realistic genotype data", {
  set.seed(2024)
  n <- 500
  p <- 20
  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n)
  beta_true <- rep(0, p)
  beta_true[c(3, 15)] <- c(0.4, 0.2)
  y <- X %*% beta_true + rnorm(n)
  bhat <- as.vector(cor(y, X))
  R <- cor(X)
  # Mixed arrays with varying sample sizes
  arr <- rep(c(1, 2), length.out = p)
  per_n <- rep(c(400, 500), length.out = p)
  result <- sdpr(bhat = bhat, LD = list(blk1 = R), n = n,
                 per_variant_sample_size = per_n, array = arr, opt_llk = 2,
                 iter = 100, burn = 30, thin = 5, verbose = FALSE, seed = 42L)
  expect_equal(length(result$beta_est), p)
  expect_true(all(is.finite(result$beta_est)))
})

# ---- sdpr_weights (wrapper) ----
test_that("sdpr_weights calls sdpr and returns beta_est", {
  set.seed(42)
  p <- 10
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  stat <- list(b = bhat, n = rep(100, p))
  result <- sdpr_weights(stat = stat, LD = R,
                         iter = 50, burn = 10, thin = 2, verbose = FALSE)
  expect_equal(length(result), p)
  expect_true(is.numeric(result))
})
