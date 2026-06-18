context("regularized_regression - susie")

# ---- initPriorSd ----
test_that("initPriorSd returns correct number of values with expected properties", {
  set.seed(123)
  n <- 50
  p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 2 + rnorm(n)

  result <- initPriorSd(X, y)
  expect_length(result, 30)
  expect_equal(result[1], 0)
  expect_true(all(diff(result) > 0))

  result_15 <- initPriorSd(X, y, n = 15)
  expect_length(result_15, 15)
  expect_equal(result_15[1], 0)
  expect_true(all(diff(result_15) > 0))
})

# ---- susieWeights ----
test_that("susieWeights returns zeros when no alpha/mu/X_column_scale_factors", {
  mock_fit <- list(pip = c(0.1, 0.2, 0.3))
  result <- susieWeights(susieFit = mock_fit)
  expect_equal(result, rep(0, 3))
})

test_that("susieWeights errors on dimension mismatch", {
  mock_fit <- list(pip = c(0.1, 0.2))
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  expect_error(susieWeights(X = X, susieFit = mock_fit), "Dimension mismatch")
})

test_that("susieWeights with alpha/mu/X_column_scale_factors calls coef.susie", {
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
  result <- susieWeights(susieFit = mock_fit)
  expect_length(result, p)
  expect_true(is.numeric(result))
  expect_true(any(result != 0))
})

test_that("susieWeights calls susie when susie_fit is NULL", {
  set.seed(42)
  p <- 5
  n <- 50
  X <- matrix(rnorm(n * p), nrow = n)
  y <- rnorm(n)
  local_mocked_bindings(
    susie = function(...) {
      list(pip = rep(0.1, p))
    }
  )
  result <- susieWeights(X = X, y = y)
  expect_equal(result, rep(0, p))
})

# ---- susieAshWeights ----
test_that("susieAshWeights returns zeros when no expected fields", {
  mock_fit <- list(pip = c(0.1, 0.2, 0.3))
  result <- susieAshWeights(susieAshFit = mock_fit)
  expect_equal(result, rep(0, 3))
})

test_that("susieAshWeights errors on dimension mismatch", {
  mock_fit <- list(pip = c(0.1, 0.2))
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  expect_error(susieAshWeights(X = X, susieAshFit = mock_fit), "Dimension mismatch")
})

test_that("susieAshWeights with proper fields calls coef.susie", {
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
  result <- susieAshWeights(susieAshFit = mock_fit)
  expect_true(is.numeric(result))
  expect_true(length(result) >= p)
})

test_that("susieAshWeights calls susie when fit is NULL", {
  set.seed(42)
  p <- 4
  n <- 30
  X <- matrix(rnorm(n * p), nrow = n)
  y <- rnorm(n)
  local_mocked_bindings(
    susie = function(...) {
      list(pip = rep(0.1, p))
    }
  )
  result <- susieAshWeights(X = X, y = y)
  expect_equal(result, rep(0, p))
})

# ---- susieInfWeights ----
test_that("susieInfWeights returns zeros when no expected fields", {
  mock_fit <- list(pip = c(0.4, 0.5))
  result <- susieInfWeights(susieInfFit = mock_fit)
  expect_equal(result, rep(0, 2))
})

test_that("susieInfWeights errors on dimension mismatch", {
  mock_fit <- list(pip = c(0.1, 0.2))
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  expect_error(susieInfWeights(X = X, susieInfFit = mock_fit), "Dimension mismatch")
})

test_that("susieInfWeights with proper fields calls coef.susie", {
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
  result <- susieInfWeights(susieInfFit = mock_fit)
  expect_true(is.numeric(result))
  expect_true(length(result) >= p)
})

test_that("susieInfWeights calls susie when fit is NULL", {
  set.seed(42)
  p <- 4
  n <- 30
  X <- matrix(rnorm(n * p), nrow = n)
  y <- rnorm(n)
  local_mocked_bindings(
    susie = function(...) {
      list(pip = rep(0.1, p))
    }
  )
  result <- susieInfWeights(X = X, y = y)
  expect_equal(result, rep(0, p))
})
