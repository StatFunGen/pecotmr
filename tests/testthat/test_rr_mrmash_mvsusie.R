context("regularized_regression - mrmash / mvsusie")

# ---- mrmashWeights ----
test_that("mrmashWeights errors when mr.mashr package is not available", {
  skip_if(requireNamespace("mr.mashr", quietly = TRUE),
          "mr.mashr is installed; skipping missing-package test")

  expect_error(
    mrmashWeights(mrmashFit = NULL, X = matrix(1, 10, 5), Y = matrix(1, 10, 3)),
    "mr\\.mash\\.alpha"
  )
})

test_that("mrmashWeights errors when X and Y are NULL and fit is NULL", {
  skip_if_not(requireNamespace("mr.mashr", quietly = TRUE),
              "mr.mashr not installed")
  expect_error(mrmashWeights(mrmashFit = NULL, X = NULL, Y = NULL),
               "Both X and Y must be provided")
})

# ---- mvsusieWeights ----
test_that("mvsusieWeights errors when mvsusieR package is not available", {
  skip_if(requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR is installed; skipping missing-package test")

  expect_error(
    mvsusieWeights(mvsusieFit = NULL, X = matrix(1, 10, 5), Y = matrix(1, 10, 3)),
    "mvsusieR"
  )
})

test_that("mvsusieWeights errors when X and Y are NULL and fit is NULL", {
  skip_if_not(requireNamespace("mvsusieR", quietly = TRUE),
              "mvsusieR not installed")
  expect_error(mvsusieWeights(mvsusieFit = NULL, X = NULL, Y = NULL),
               "Both X and Y must be provided")
})

test_that("mvsusieWeights fits model and returns coefficients when fit is NULL", {
  skip_if_not(requireNamespace("mvsusieR", quietly = TRUE),
              "mvsusieR not installed")
  set.seed(42)
  n <- 30
  p <- 5
  R <- 3
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * R), n, R)
  fake_coef <- matrix(rnorm((p + 1) * R), nrow = p + 1, ncol = R)
  captured <- list()

  local_mocked_bindings(
    create_mixture_prior = function(...) list(),
    mvsusie = function(...) {
      captured <<- list(...)
      "mock_fit"
    },
    coef.mvsusie = function(...) fake_coef,
    .package = "mvsusieR"
  )

  result <- expect_message(
    mvsusieWeights(X = X, Y = Y, L = 12, LGreedy = 4),
    "mvsusieFit is not provided"
  )
  # Should return coef without intercept row
  expect_equal(dim(result), c(p, R))
  expect_equal(result, fake_coef[-1, ])
  expect_equal(captured$L, 12)
  expect_equal(captured$L_greedy, 4)
})

test_that("mvsusieWeights returns coefficients from provided fit", {
  skip_if_not(requireNamespace("mvsusieR", quietly = TRUE),
              "mvsusieR not installed")
  p <- 5
  R <- 3
  fake_coef <- matrix(rnorm((p + 1) * R), nrow = p + 1, ncol = R)

  local_mocked_bindings(
    coef.mvsusie = function(...) fake_coef,
    .package = "mvsusieR"
  )

  result <- mvsusieWeights(mvsusieFit = "precomputed_fit")
  expect_equal(dim(result), c(p, R))
  expect_equal(result, fake_coef[-1, ])
})
