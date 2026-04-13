context("regularized_regression — l0learn")

# ---- l0learn_weights ----
test_that("l0learn_weights computes weights with default L0 penalty", {
  skip_if_not_installed("L0Learn")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- l0learn_weights(X, y)
  expect_equal(nrow(result), p)
  expect_equal(ncol(result), 1)
  expect_true(all(is.finite(result)))
})

test_that("l0learn_weights computes weights with L0L2 penalty", {
  skip_if_not_installed("L0Learn")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- l0learn_weights(X, y, penalty = "L0L2")
  expect_equal(nrow(result), p)
  expect_equal(ncol(result), 1)
  expect_true(all(is.finite(result)))
})

test_that("l0learn_weights warns and zero-pads when X has zero-variance columns", {
  skip_if_not_installed("L0Learn")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  X[, 6] <- 4
  y <- X[, 1] * 0.5 + rnorm(n)
  expect_warning(result <- l0learn_weights(X, y),
                 "l0learn_weights: dropping 1 zero-variance column")
  expect_equal(nrow(result), p)
  expect_equal(result[6, 1], 0)
})

test_that("l0learn_weights passes nGamma through to L0Learn.cvfit", {
  skip_if_not_installed("L0Learn")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  # nGamma = 0 is invalid; the underlying function should reject it,
  # proving the argument reached L0Learn.cvfit.
  expect_error(l0learn_weights(X, y, penalty = "L0L2", nGamma = 0))
})
