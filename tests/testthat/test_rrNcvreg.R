context("regularized_regression - ncvreg")

# ---- ncvregWeights / scadWeights / mcpWeights ----
test_that("ncvregWeights computes weights with SCAD penalty", {
  skip_if_not_installed("ncvreg")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- ncvregWeights(X, y, penalty = "SCAD")
  expect_equal(nrow(result), p)
  expect_equal(ncol(result), 1)
  expect_true(all(is.finite(result)))
})

test_that("ncvregWeights warns and zero-pads when X has zero-variance columns", {
  skip_if_not_installed("ncvreg")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  X[, 5] <- 3
  y <- X[, 1] * 0.5 + rnorm(n)
  expect_warning(result <- ncvregWeights(X, y, penalty = "SCAD"),
                 "ncvregWeights: dropping 1 zero-variance column")
  expect_equal(nrow(result), p)
  expect_equal(result[5, 1], 0)
})

test_that("scadWeights computes weights and dispatches to ncvregWeights", {
  skip_if_not_installed("ncvreg")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- scadWeights(X, y)
  expect_equal(nrow(result), p)
  expect_equal(ncol(result), 1)
  expect_true(all(is.finite(result)))
})

test_that("mcpWeights computes weights and dispatches to ncvregWeights", {
  skip_if_not_installed("ncvreg")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- mcpWeights(X, y)
  expect_equal(nrow(result), p)
  expect_equal(ncol(result), 1)
  expect_true(all(is.finite(result)))
})

test_that("scadWeights passes nfolds through to cv.ncvreg", {
  skip_if_not_installed("ncvreg")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    cv.ncvreg = function(X, y, penalty, nfolds = 5, ...) {
      captured$nfolds <- nfolds
      stop("STOP_AFTER_CAPTURE")
    },
    .package = "ncvreg"
  )
  expect_error(scadWeights(X, y, nfolds = 7), "STOP_AFTER_CAPTURE")
  expect_equal(captured$nfolds, 7)
})
