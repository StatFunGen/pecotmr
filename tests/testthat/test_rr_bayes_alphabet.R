context("regularized_regression — bayes_alphabet")

# ---- bayes_alphabet_weights ----
test_that("bayes_alphabet_weights errors on dimension mismatch", {
  skip_if_not_installed("qgg")
  X <- matrix(rnorm(100), nrow = 10)
  y <- rnorm(5)
  expect_error(bayes_alphabet_weights(X, y, method = "bayesN"), "same number of rows")
})

test_that("bayes_alphabet_weights errors on covariate dimension mismatch", {
  skip_if_not_installed("qgg")
  X <- matrix(rnorm(100), nrow = 10)
  y <- rnorm(10)
  Z <- matrix(rnorm(15), nrow = 5)
  expect_error(bayes_alphabet_weights(X, y, method = "bayesN", Z = Z),
               "same number of rows")
})

test_that("bayes_alphabet_weights does not error on valid Z dimensions", {
  skip_if_not_installed("qgg")
  set.seed(42)
  n <- 50
  p <- 10
  q <- 3
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  Z <- matrix(round(runif(n * q, 0, 0.8), 0), nrow = n)

  result <- bayes_alphabet_weights(X, y, method = "bayesN", Z = Z, nit = 50, nburn = 10)
  expect_equal(length(result), p)
  expect_true(is.numeric(result))
})

test_that("bayes_alphabet_weights with Z = NULL does not error on dimension check", {
  skip_if_not_installed("qgg")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- bayes_alphabet_weights(X, y, method = "bayesN", Z = NULL, nit = 50, nburn = 10)
  expect_equal(length(result), p)
  expect_true(is.numeric(result))
})

# ---- bayes_n/l/a/c/r_weights (wrapper dispatchers) ----
test_that("bayes_n_weights dispatches to bayes_alphabet_weights with bayesN", {
  skip_if_not_installed("qgg")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- bayes_n_weights(X, y, nit = 50, nburn = 10)
  expect_equal(length(result), p)
  expect_true(is.numeric(result))
  expect_true(all(is.finite(result)))
})

# ---- bayes_alphabet_weights zero-variance handling ----
test_that("bayes_alphabet_weights warns and zero-pads when X has zero-variance columns", {
  skip_if_not_installed("qgg")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  X[, 4] <- 2
  y <- X[, 1] * 0.5 + rnorm(n)
  expect_warning(result <- bayes_alphabet_weights(X, y, method = "bayesN", nit = 50, nburn = 10),
                 "bayes_alphabet_weights: dropping 1 zero-variance column")
  expect_equal(length(result), p)
  expect_equal(result[4], 0)
  expect_true(all(is.finite(result)))
})
