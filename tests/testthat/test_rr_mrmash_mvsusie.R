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
