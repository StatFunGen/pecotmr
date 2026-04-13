context("regularized_regression — bglr")

# ---- bglr_weights / bayes_b_weights / b_lasso_weights ----
test_that("bglr_weights computes weights with BayesB model", {
  skip_if_not_installed("BGLR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- bglr_weights(X, y, model = "BayesB", nIter = 100, burnIn = 20, thin = 2,
                         eta_args = list(probIn = 0.05))
  expect_equal(length(result), p)
  expect_true(is.numeric(result))
  expect_true(all(is.finite(result)))
})

test_that("bayes_b_weights computes weights and dispatches to bglr_weights", {
  skip_if_not_installed("BGLR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- bayes_b_weights(X, y, nIter = 100, burnIn = 20, thin = 2)
  expect_equal(length(result), p)
  expect_true(all(is.finite(result)))
})

test_that("bayes_b_weights passes probIn through to BGLR ETA", {
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

  result <- bayes_b_weights(X, y, nIter = 100, burnIn = 20, thin = 2, probIn = 0.42)
  expect_equal(length(result), p)
  expect_equal(captured$eta[[1]]$model, "BayesB")
  expect_equal(captured$eta[[1]]$probIn, 0.42)
})

test_that("b_lasso_weights computes weights with BL model", {
  skip_if_not_installed("BGLR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- b_lasso_weights(X, y, nIter = 100, burnIn = 20, thin = 2)
  expect_equal(length(result), p)
  expect_true(all(is.finite(result)))
})

test_that("bglr_weights warns and zero-pads when X has zero-variance columns", {
  skip_if_not_installed("BGLR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  X[, 7] <- 9
  y <- X[, 1] * 0.5 + rnorm(n)
  expect_warning(
    result <- bayes_b_weights(X, y, nIter = 100, burnIn = 20, thin = 2),
    "bglr_weights: dropping 1 zero-variance column"
  )
  expect_equal(length(result), p)
  expect_equal(result[7], 0)
})

test_that("bglr_weights cleans up its tempdir on exit", {
  skip_if_not_installed("BGLR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  before <- list.files(tempdir(), pattern = "^bglr_")
  bayes_b_weights(X, y, nIter = 100, burnIn = 20, thin = 2)
  after <- list.files(tempdir(), pattern = "^bglr_")
  expect_setequal(before, after)
})
