context("regularized_regression — dpr")

# ---- dpr_weights ----
test_that("dpr_weights computes weights with VB fitting method", {
  skip_if_not_installed("RcppDPR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- dpr_weights(X, y, fitting_method = "VB")
  expect_equal(length(result), p)
  expect_true(is.numeric(result))
  expect_true(all(is.finite(result)))
})

test_that("dpr_weights computes weights with Gibbs fitting method", {
  skip_if_not_installed("RcppDPR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- dpr_weights(X, y, fitting_method = "Gibbs")
  expect_equal(length(result), p)
  expect_true(is.numeric(result))
  expect_true(all(is.finite(result)))
})

test_that("dpr_weights warns and zero-pads when X has zero-variance columns", {
  skip_if_not_installed("RcppDPR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  X[, 8] <- 6
  y <- X[, 1] * 0.5 + rnorm(n)
  expect_warning(result <- dpr_weights(X, y),
                 "dpr_weights: dropping 1 zero-variance column")
  expect_equal(length(result), p)
  expect_equal(result[8], 0)
})
