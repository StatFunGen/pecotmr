context("regularized_regression — lassosum_rss")

# ---- lassosum_rss ----
test_that("lassosum_rss errors on invalid LD input", {
  expect_error(lassosum_rss(bhat = rnorm(5), LD = "not_a_list", n = 100),
               "valid list of LD blocks")
})

test_that("lassosum_rss errors on non-positive sample size", {
  expect_error(lassosum_rss(bhat = rnorm(5), LD = list(blk1 = diag(5)), n = -1),
               "valid sample size")
})

test_that("lassosum_rss errors on mismatched bhat and LD dimensions", {
  expect_error(
    lassosum_rss(bhat = rnorm(10), LD = list(blk1 = diag(5)), n = 100),
    "same as the sum"
  )
})

test_that("lassosum_rss runs successfully with valid input", {
  set.seed(42)
  p <- 10
  n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  for (i in 1:(p - 1)) {
    R[i, i + 1] <- 0.3
    R[i + 1, i] <- 0.3
  }
  result <- lassosum_rss(bhat = bhat, LD = list(blk1 = R), n = n)
  expect_type(result, "list")
  expect_true("beta_est" %in% names(result))
  expect_equal(length(result$beta_est), p)
  expect_true(all(is.finite(result$beta_est)))
  expect_equal(nrow(result$beta), p)
  expect_equal(ncol(result$beta), 20)
})

test_that("lassosum_rss accepts multiple LD blocks", {
  set.seed(42)
  p1 <- 5
  p2 <- 5
  p <- p1 + p2
  n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R1 <- diag(p1)
  R2 <- diag(p2)
  result <- lassosum_rss(bhat = bhat, LD = list(blk1 = R1, blk2 = R2), n = n)
  expect_type(result, "list")
  expect_true("beta_est" %in% names(result))
  expect_equal(length(result$beta_est), p)
  expect_true(all(is.finite(result$beta_est)))
})

test_that("lassosum_rss with large lambda gives all-zero betas", {
  set.seed(42)
  p <- 10
  n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  result <- lassosum_rss(bhat = bhat, LD = list(blk1 = R), n = n,
                         lambda = c(100))
  expect_true(all(result$beta_est == 0))
})

# ---- lassosum_rss_weights (wrapper) ----
test_that("lassosum_rss_weights calls lassosum_rss and returns beta_est", {
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
  result <- lassosum_rss_weights(stat = stat, LD = R)
  expect_equal(length(result), p)
  expect_true(is.numeric(result))
})

test_that("lassosum_rss_weights rescales bhat when max(abs(bhat)) >= 1", {
  set.seed(42)
  p <- 5
  n <- 100
  # Construct bhat with max abs >= 1 to trigger the clamp branch
  bhat <- c(1.5, -0.2, 0.1, 0.05, -0.3)
  R <- diag(p)
  stat <- list(b = bhat, n = rep(n, p))
  captured <- NULL
  local_mocked_bindings(
    lassosum_rss = function(bhat, LD, n, ...) {
      captured <<- bhat
      list(beta_est = rep(0, length(bhat)), fbeta = 1)
    }
  )
  result <- lassosum_rss_weights(stat = stat, LD = R, s = 0.5)
  # Captured bhat should have been rescaled so max abs is just under 1
  expect_true(max(abs(captured)) < 1)
  expect_equal(max(abs(captured)), 0.9999, tolerance = 1e-9)
  # Sign-preserving rescale
  expect_true(all(sign(captured) == sign(bhat)))
})
