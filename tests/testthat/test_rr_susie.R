context("regularized_regression — susie")

# ---- init_prior_sd ----
test_that("init_prior_sd returns correct number of values with expected properties", {
  set.seed(123)
  n <- 50
  p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 2 + rnorm(n)

  result <- init_prior_sd(X, y)
  expect_length(result, 30)
  expect_equal(result[1], 0)
  expect_true(all(diff(result) > 0))

  result_15 <- init_prior_sd(X, y, n = 15)
  expect_length(result_15, 15)
  expect_equal(result_15[1], 0)
  expect_true(all(diff(result_15) > 0))
})

# ---- susie_weights ----
test_that("susie_weights returns zeros when no alpha/mu/X_column_scale_factors", {
  mock_fit <- list(pip = c(0.1, 0.2, 0.3))
  result <- susie_weights(susie_fit = mock_fit)
  expect_equal(result, rep(0, 3))
})

test_that("susie_weights errors on dimension mismatch", {
  mock_fit <- list(pip = c(0.1, 0.2))
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  expect_error(susie_weights(X = X, susie_fit = mock_fit), "Dimension mismatch")
})

test_that("susie_weights with alpha/mu/X_column_scale_factors calls coef.susie", {
  p <- 5
  L <- 2
  mock_fit <- list(
    alpha = matrix(c(
      0.9, 0.05, 0.02, 0.02, 0.01,
      0.1, 0.1, 0.6, 0.1, 0.1
    ), nrow = L, ncol = p, byrow = TRUE),
    mu = matrix(c(
      2.0, 0.1, 0.05, 0.03, 0.01,
      0.5, 0.2, 1.5, 0.1, 0.05
    ), nrow = L, ncol = p, byrow = TRUE),
    X_column_scale_factors = rep(1.0, p),
    pip = runif(p),
    intercept = 0
  )
  result <- susie_weights(susie_fit = mock_fit)
  expect_length(result, p)
  expect_true(is.numeric(result))
  expect_true(any(result != 0))
})

test_that("susie_weights calls susie_wrapper when susie_fit is NULL", {
  set.seed(42)
  p <- 5
  n <- 50
  X <- matrix(rnorm(n * p), nrow = n)
  y <- rnorm(n)
  local_mocked_bindings(
    susie_wrapper = function(...) {
      list(pip = rep(0.1, p))
    }
  )
  result <- susie_weights(X = X, y = y)
  expect_equal(result, rep(0, p))
})

# ---- susie_ash_weights ----
test_that("susie_ash_weights returns zeros when no expected fields", {
  mock_fit <- list(pip = c(0.1, 0.2, 0.3))
  result <- susie_ash_weights(susie_ash_fit = mock_fit)
  expect_equal(result, rep(0, 3))
})

test_that("susie_ash_weights errors on dimension mismatch", {
  mock_fit <- list(pip = c(0.1, 0.2))
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  expect_error(susie_ash_weights(X = X, susie_ash_fit = mock_fit), "Dimension mismatch")
})

test_that("susie_ash_weights with proper fields calls coef.susie", {
  p <- 4
  L <- 2
  mock_fit <- list(
    alpha = matrix(runif(L * p), nrow = L, ncol = p),
    mu = matrix(rnorm(L * p), nrow = L, ncol = p),
    theta = matrix(rnorm(L * p), nrow = L, ncol = p),
    X_column_scale_factors = rep(1.0, p),
    pip = runif(p),
    intercept = 0
  )
  result <- susie_ash_weights(susie_ash_fit = mock_fit)
  expect_true(is.numeric(result))
  expect_true(length(result) >= p)
})

test_that("susie_ash_weights calls susie_wrapper when fit is NULL", {
  set.seed(42)
  p <- 4
  n <- 30
  X <- matrix(rnorm(n * p), nrow = n)
  y <- rnorm(n)
  local_mocked_bindings(
    susie_wrapper = function(...) {
      list(pip = rep(0.1, p))
    }
  )
  result <- susie_ash_weights(X = X, y = y)
  expect_equal(result, rep(0, p))
})

# ---- susie_inf_weights ----
test_that("susie_inf_weights returns zeros when no expected fields", {
  mock_fit <- list(pip = c(0.4, 0.5))
  result <- susie_inf_weights(susie_inf_fit = mock_fit)
  expect_equal(result, rep(0, 2))
})

test_that("susie_inf_weights errors on dimension mismatch", {
  mock_fit <- list(pip = c(0.1, 0.2))
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  expect_error(susie_inf_weights(X = X, susie_inf_fit = mock_fit), "Dimension mismatch")
})

test_that("susie_inf_weights with proper fields calls coef.susie", {
  p <- 4
  L <- 2
  mock_fit <- list(
    alpha = matrix(runif(L * p), nrow = L, ncol = p),
    mu = matrix(rnorm(L * p), nrow = L, ncol = p),
    theta = matrix(rnorm(L * p), nrow = L, ncol = p),
    X_column_scale_factors = rep(1.0, p),
    pip = runif(p),
    intercept = 0
  )
  result <- susie_inf_weights(susie_inf_fit = mock_fit)
  expect_true(is.numeric(result))
  expect_true(length(result) >= p)
})

test_that("susie_inf_weights calls susie_wrapper when fit is NULL", {
  set.seed(42)
  p <- 4
  n <- 30
  X <- matrix(rnorm(n * p), nrow = n)
  y <- rnorm(n)
  local_mocked_bindings(
    susie_wrapper = function(...) {
      list(pip = rep(0.1, p))
    }
  )
  result <- susie_inf_weights(X = X, y = y)
  expect_equal(result, rep(0, p))
})
