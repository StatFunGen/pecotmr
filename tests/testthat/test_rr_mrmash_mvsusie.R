context("regularized_regression — mrmash / mvsusie")

# ---- mrmash_weights ----
test_that("mrmash_weights errors when mr.mashr package is not available", {
  skip_if(requireNamespace("mr.mashr", quietly = TRUE),
          "mr.mashr is installed; skipping missing-package test")

  expect_error(
    mrmash_weights(mrmash_fit = NULL, X = matrix(1, 10, 5), Y = matrix(1, 10, 3)),
    "mr\\.mash\\.alpha"
  )
})

test_that("mrmash_weights errors when X and Y are NULL and fit is NULL", {
  skip_if_not(requireNamespace("mr.mashr", quietly = TRUE),
              "mr.mashr not installed")
  expect_error(mrmash_weights(mrmash_fit = NULL, X = NULL, Y = NULL),
               "Both X and Y must be provided")
})

# ---- mvsusie_weights ----
test_that("mvsusie_weights errors when mvsusieR package is not available", {
  skip_if(requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR is installed; skipping missing-package test")

  expect_error(
    mvsusie_weights(mvsusie_fit = NULL, X = matrix(1, 10, 5), Y = matrix(1, 10, 3)),
    "mvsusieR"
  )
})

test_that("mvsusie_weights errors when X and Y are NULL and fit is NULL", {
  skip_if_not(requireNamespace("mvsusieR", quietly = TRUE),
              "mvsusieR not installed")
  expect_error(mvsusie_weights(mvsusie_fit = NULL, X = NULL, Y = NULL),
               "Both X and Y must be provided")
})

test_that("mvsusie_weights fits model and returns coefficients when fit is NULL", {
  skip_if_not(requireNamespace("mvsusieR", quietly = TRUE),
              "mvsusieR not installed")
  set.seed(42)
  n <- 30
  p <- 5
  R <- 3
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * R), n, R)
  fake_coef <- matrix(rnorm((p + 1) * R), nrow = p + 1, ncol = R)

  local_mocked_bindings(
    create_mixture_prior = function(...) list(),
    mvsusie = function(...) "mock_fit",
    coef.mvsusie = function(...) fake_coef,
    .package = "mvsusieR"
  )

  result <- expect_message(
    mvsusie_weights(X = X, Y = Y),
    "mvsusie_fit is not provided"
  )
  # Should return coef without intercept row
  expect_equal(dim(result), c(p, R))
  expect_equal(result, fake_coef[-1, ])
})

test_that("mvsusie_weights returns coefficients from provided fit", {
  skip_if_not(requireNamespace("mvsusieR", quietly = TRUE),
              "mvsusieR not installed")
  p <- 5
  R <- 3
  fake_coef <- matrix(rnorm((p + 1) * R), nrow = p + 1, ncol = R)

  local_mocked_bindings(
    coef.mvsusie = function(...) fake_coef,
    .package = "mvsusieR"
  )

  result <- mvsusie_weights(mvsusie_fit = "precomputed_fit")
  expect_equal(dim(result), c(p, R))
  expect_equal(result, fake_coef[-1, ])
})
