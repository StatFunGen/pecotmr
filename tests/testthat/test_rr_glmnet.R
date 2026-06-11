context("regularized_regression - glmnet")

# ---- glmnetWeights ----
test_that("glmnetWeights computes LASSO weights", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- glmnetWeights(X, y, alpha = 1)
  expect_equal(nrow(result), p)
  expect_equal(ncol(result), 1)
})

test_that("enetWeights computes elastic net weights", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- enetWeights(X, y)
  expect_equal(nrow(result), p)
})

test_that("lassoWeights computes LASSO weights", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- lassoWeights(X, y)
  expect_equal(nrow(result), p)
})

test_that("glmnetWeights handles zero-variance columns", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  X[, 3] <- 5
  y <- X[, 1] * 0.5 + rnorm(n)
  expect_warning(result <- glmnetWeights(X, y, alpha = 1),
                 "glmnetWeights: dropping 1 zero-variance column")
  expect_equal(result[3, 1], 0)
})

test_that("glmnetWeights errors when all columns are constant", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50
  p <- 5
  X <- matrix(rep(1:p, each = n), nrow = n, ncol = p)
  y <- rnorm(n)

  expect_error(glmnetWeights(X, y, alpha = 1), "matrix with 2 or more columns")
})

test_that("glmnetWeights handles NA-producing columns gracefully", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50
  p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)

  X[1, 4] <- NA

  expect_warning(result <- glmnetWeights(X, y, alpha = 1),
                 "glmnetWeights: dropping 1 zero-variance column")
  expect_equal(nrow(result), p)
  expect_equal(ncol(result), 1)
  expect_equal(result[4, 1], 0)
})

test_that("glmnetWeights with alpha = 0 (ridge regression) works", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- glmnetWeights(X, y, alpha = 0)
  expect_equal(nrow(result), p)
})
