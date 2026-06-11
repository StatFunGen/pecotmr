context("regularized_regression - mrash")

# ---- mrashWeights ----
test_that("mrashWeights calls lassoWeights as default beta.init", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  local_mocked_bindings(
    lassoWeights = function(X, y) rep(0.01, ncol(X)),
    initPriorSd = function(X, y, n = 30) seq(0, 3, length.out = n)
  )
  result <- mrashWeights(X, y)
  expect_true(is.numeric(result))
  expect_equal(length(result), p)
})

test_that("mrashWeights with initPriorSd = FALSE passes NULL sa2", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  local_mocked_bindings(
    lassoWeights = function(X, y) rep(0.01, ncol(X))
  )
  result <- mrashWeights(X, y, initPriorSd = FALSE)
  expect_true(is.numeric(result))
  expect_equal(length(result), p)
})

test_that("mrashWeights subsets a user-supplied beta.init of length ncol(X) by keep", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  X[, 5] <- 7  # zero-variance column to force subset by keep
  y <- X[, 1] * 0.5 + rnorm(n)
  user_beta_init <- seq(0.01, by = 0.01, length.out = p)
  expect_warning(
    result <- mrashWeights(X, y, beta.init = user_beta_init),
    "mrashWeights: dropping 1 zero-variance column"
  )
  expect_equal(length(result), p)
  expect_equal(result[5], 0)
})

test_that("mrashWeights warns and zero-pads when X has zero-variance columns", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  X[, 3] <- 7
  y <- X[, 1] * 0.5 + rnorm(n)
  local_mocked_bindings(
    lassoWeights = function(X, y) rep(0.01, ncol(X))
  )
  expect_warning(result <- mrashWeights(X, y),
                 "mrashWeights: dropping 1 zero-variance column")
  expect_equal(length(result), p)
  expect_equal(result[3], 0)
})

gc()
