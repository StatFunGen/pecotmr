context("regularized_regression - dpr")

# ---- dprWeights ----
test_that("dprWeights computes weights with VB fitting method", {
  skip_if_not_installed("RcppDPR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- dprWeights(X, y, fittingMethod = "VB")
  expect_equal(length(result), p)
  expect_true(is.numeric(result))
  expect_true(all(is.finite(result)))
})

test_that("dprWeights computes weights with Gibbs fitting method", {
  skip_if_not_installed("RcppDPR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- dprWeights(X, y, fittingMethod = "Gibbs")
  expect_equal(length(result), p)
  expect_true(is.numeric(result))
  expect_true(all(is.finite(result)))
})

test_that("dprWeights warns and zero-pads when X has zero-variance columns", {
  skip_if_not_installed("RcppDPR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  X[, 8] <- 6
  y <- X[, 1] * 0.5 + rnorm(n)
  expect_warning(result <- dprWeights(X, y),
                 "dprWeights: dropping 1 zero-variance column")
  expect_equal(length(result), p)
  expect_equal(result[8], 0)
})
