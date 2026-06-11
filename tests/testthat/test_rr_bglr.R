context("regularized_regression - bglr")

# ---- bglrWeights / bayesBWeights / bLassoWeights ----
test_that("bglrWeights computes weights with BayesB model", {
  skip_if_not_installed("BGLR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- bglrWeights(X, y, model = "BayesB", nIter = 100, burnIn = 20, thin = 2,
                         etaArgs = list(probIn = 0.05))
  expect_equal(length(result), p)
  expect_true(is.numeric(result))
  expect_true(all(is.finite(result)))
})

test_that("bayesBWeights computes weights and dispatches to bglrWeights", {
  skip_if_not_installed("BGLR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- bayesBWeights(X, y, nIter = 100, burnIn = 20, thin = 2)
  expect_equal(length(result), p)
  expect_true(all(is.finite(result)))
})

test_that("bayesBWeights passes probIn through to BGLR ETA", {
  skip_if_not_installed("BGLR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)

  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    BGLR = function(y, ETA, ...) {
      captured$eta <- ETA
      list(ETA = list(list(b = rep(0, ncol(ETA[[1]]$X)))))
    },
    .package = "BGLR"
  )

  result <- bayesBWeights(X, y, nIter = 100, burnIn = 20, thin = 2, probIn = 0.42)
  expect_equal(length(result), p)
  expect_equal(captured$eta[[1]]$model, "BayesB")
  expect_equal(captured$eta[[1]]$probIn, 0.42)
})

test_that("bLassoWeights computes weights with BL model", {
  skip_if_not_installed("BGLR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- bLassoWeights(X, y, nIter = 100, burnIn = 20, thin = 2)
  expect_equal(length(result), p)
  expect_true(all(is.finite(result)))
})

test_that("bglrWeights warns and zero-pads when X has zero-variance columns", {
  skip_if_not_installed("BGLR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  X[, 7] <- 9
  y <- X[, 1] * 0.5 + rnorm(n)
  expect_warning(
    result <- bayesBWeights(X, y, nIter = 100, burnIn = 20, thin = 2),
    "bglrWeights: dropping 1 zero-variance column"
  )
  expect_equal(length(result), p)
  expect_equal(result[7], 0)
})

test_that("bglrWeights cleans up its tempdir on exit", {
  skip_if_not_installed("BGLR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  before <- list.files(tempdir(), pattern = "^bglr_")
  bayesBWeights(X, y, nIter = 100, burnIn = 20, thin = 2)
  after <- list.files(tempdir(), pattern = "^bglr_")
  expect_setequal(before, after)
})
