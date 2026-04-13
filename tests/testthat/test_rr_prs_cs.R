context("regularized_regression — prs_cs")

# ---- prs_cs ----
test_that("prs_cs errors on invalid LD input", {
  expect_error(prs_cs(bhat = rnorm(5), LD = "not_a_list", n = 100),
               "valid list of LD blocks")
})

test_that("prs_cs errors on non-positive sample size", {
  expect_error(prs_cs(bhat = rnorm(5), LD = list(blk1 = diag(5)), n = -1),
               "valid sample size")
})

test_that("prs_cs errors on mismatched maf length", {
  expect_error(
    prs_cs(bhat = rnorm(5), LD = list(blk1 = diag(5)), n = 100, maf = rep(0.3, 3)),
    "same as 'maf'"
  )
})

test_that("prs_cs errors on mismatched bhat and LD dimensions", {
  expect_error(
    prs_cs(bhat = rnorm(10), LD = list(blk1 = diag(5)), n = 100),
    "same as the sum"
  )
})

test_that("prs_cs runs successfully with valid input", {
  set.seed(42)
  p <- 10
  n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  for (i in 1:(p - 1)) {
    R[i, i + 1] <- 0.3
    R[i + 1, i] <- 0.3
  }
  result <- prs_cs(bhat = bhat, LD = list(blk1 = R), n = n,
                   maf = rep(0.3, p), n_iter = 50, n_burnin = 10, thin = 2)
  expect_type(result, "list")
  expect_true("beta_est" %in% names(result))
  expect_equal(length(result$beta_est), p)
  expect_true(all(is.finite(result$beta_est)))
})

test_that("prs_cs accepts multiple LD blocks whose dimensions sum to length of bhat", {
  set.seed(42)
  p1 <- 5
  p2 <- 5
  p <- p1 + p2
  n <- 100

  bhat <- rnorm(p, sd = 0.1)
  R1 <- diag(p1)
  R2 <- diag(p2)

  result <- prs_cs(
    bhat = bhat,
    LD = list(blk1 = R1, blk2 = R2),
    n = n,
    maf = rep(0.3, p),
    n_iter = 50, n_burnin = 10, thin = 2
  )

  expect_type(result, "list")
  expect_true("beta_est" %in% names(result))
  expect_equal(length(result$beta_est), p)
  expect_true(all(is.finite(result$beta_est)))
})

test_that("prs_cs with phi = NULL estimates phi automatically", {
  set.seed(42)
  p <- 10; n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  result <- prs_cs(bhat = bhat, LD = list(blk1 = R), n = n,
                   phi = NULL, maf = rep(0.3, p),
                   n_iter = 50, n_burnin = 10, thin = 2)
  expect_true("phi_est" %in% names(result))
})

test_that("prs_cs with explicit phi value", {
  set.seed(42)
  p <- 10; n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  result <- prs_cs(bhat = bhat, LD = list(blk1 = R), n = n,
                   phi = 0.01, maf = rep(0.3, p),
                   n_iter = 50, n_burnin = 10, thin = 2)
  expect_true("phi_est" %in% names(result))
  expect_true("sigma_est" %in% names(result))
  expect_true("psi_est" %in% names(result))
})

test_that("prs_cs works without maf (maf = NULL)", {
  set.seed(42)
  p <- 10; n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  result <- prs_cs(bhat = bhat, LD = list(blk1 = R), n = n,
                   maf = NULL, n_iter = 50, n_burnin = 10, thin = 2)
  expect_equal(length(result$beta_est), p)
})

# ---- prs_cs signal recovery ----
test_that("prs_cs recovers signal direction on simulated genotype data", {
  set.seed(42)
  n <- 500
  p <- 20
  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n)
  beta_true <- rep(0, p)
  beta_true[c(3, 10, 15)] <- c(0.4, -0.3, 0.2)
  y <- X %*% beta_true + rnorm(n)
  bhat <- as.vector(cor(y, X))
  R <- cor(X)
  result <- prs_cs(bhat = bhat, LD = list(blk1 = R), n = n,
                   n_iter = 1000, n_burnin = 500, thin = 5, seed = 42)
  expect_true("beta_est" %in% names(result))
  expect_equal(length(result$beta_est), p)
  expect_true(all(is.finite(result$beta_est)))
  # Sigma should be reasonable (near 1 for standardized data)
  expect_true(result$sigma_est > 0.1 && result$sigma_est < 10)
  # Correlation with truth should be positive (signal recovery)
  expect_gt(cor(result$beta_est, beta_true), 0.5)
})

# ---- prs_cs_weights (wrapper) ----
test_that("prs_cs_weights calls prs_cs and returns beta_est", {
  set.seed(42)
  p <- 10
  n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  for (i in 1:(p - 1)) {
    R[i, i + 1] <- 0.3
    R[i + 1, i] <- 0.3
  }
  stat <- list(b = bhat, n = rep(n, p))
  result <- prs_cs_weights(stat = stat, LD = R,
                           maf = rep(0.3, p), n_iter = 50, n_burnin = 10, thin = 2)
  expect_equal(length(result), p)
  expect_true(is.numeric(result))
})
