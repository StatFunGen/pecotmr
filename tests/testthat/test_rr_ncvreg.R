context("regularized_regression — ncvreg")

# ---- ncvreg_weights / scad_weights / mcp_weights ----
test_that("ncvreg_weights computes weights with SCAD penalty", {
  skip_if_not_installed("ncvreg")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- ncvreg_weights(X, y, penalty = "SCAD")
  expect_equal(nrow(result), p)
  expect_equal(ncol(result), 1)
  expect_true(all(is.finite(result)))
})

test_that("ncvreg_weights warns and zero-pads when X has zero-variance columns", {
  skip_if_not_installed("ncvreg")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  X[, 5] <- 3
  y <- X[, 1] * 0.5 + rnorm(n)
  expect_warning(result <- ncvreg_weights(X, y, penalty = "SCAD"),
                 "ncvreg_weights: dropping 1 zero-variance column")
  expect_equal(nrow(result), p)
  expect_equal(result[5, 1], 0)
})

test_that("scad_weights computes weights and dispatches to ncvreg_weights", {
  skip_if_not_installed("ncvreg")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- scad_weights(X, y)
  expect_equal(nrow(result), p)
  expect_equal(ncol(result), 1)
  expect_true(all(is.finite(result)))
})

test_that("mcp_weights computes weights and dispatches to ncvreg_weights", {
  skip_if_not_installed("ncvreg")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- mcp_weights(X, y)
  expect_equal(nrow(result), p)
  expect_equal(ncol(result), 1)
  expect_true(all(is.finite(result)))
})

# test_that("scad_weights passes nfolds through to cv.ncvreg", {
#   skip_if_not_installed("ncvreg")
#   set.seed(42)
#   n <- 50; p <- 10
#   X <- matrix(rnorm(n * p), nrow = n)
#   y <- X[, 1] * 0.5 + rnorm(n)
#   # cv.ncvreg requires nfolds >= 2; passing 1 should error from inside ncvreg,
#   # which proves nfolds reached it.
#   expect_error(scad_weights(X, y, nfolds = 1))
# })
