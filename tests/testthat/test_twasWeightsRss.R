context("SS-TWAS: weights, pipeline, and omnibus combination")

# Previous TwasWeights-class tests used the legacy constructor
# `TwasWeights(weights = list(...), variantIds = ..., standardized = ...)`.
# The new `TwasWeights` is a DFrame collection class with (study,
# context, trait, method, entry) columns where each entry is a
# `TwasWeightsEntry` S4 object carrying weights / fits / cvPerformance.
# Class-shape tests for the new collection should live alongside the
# pipeline tests and assert via accessors (`getWeights`, `getStudy`,
# `getCvPerformance`, etc.) — not against legacy slot shapes.
#
# `twasAnalysis()` was collapsed into the unified `twasZ()` dispatcher
# (task #37); its tests are removed here.
# `twasWeightsSumstatPipeline()` was removed without replacement in the
# S4 refactor (twasWeightsPipeline now dispatches directly on
# `QtlSumStats` / `QtlDataset` / `MultiTaskQtlDataset`).
#
# What remains: tests of the internal SuSiE-RSS weight extractors that
# are still present in `R/twasWeights.R` (`.susieRssExtractWeights`,
# `susieRssWeights`, `susieInfRssWeights`, `fitSusieInfThenSusieRss`).

# =============================================================================
# SuSiE-RSS weight extraction
# =============================================================================

test_that(".susie_rss_extract_weights returns correct-length vector", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 20
  n <- 500
  R <- diag(p)
  z <- rnorm(p)
  w <- pecotmr:::.susieRssExtractWeights(
    fit = NULL, z = z, R = R, n = n,
    requiredFields = c("alpha", "mu", "X_column_scale_factors"),
    fitArgs = list(L = 5)
  )
  expect_equal(length(w), p)
  expect_true(all(is.finite(w)))
})

test_that("susieRssWeights follows (stat, LD) convention", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 20
  n <- 500
  R <- diag(p)
  z <- rnorm(p)
  stat <- list(b = z / sqrt(n), cor = z / sqrt(n), z = z, n = rep(n, p))
  w <- susieRssWeights(stat, R, methodArgs = list(L = 5))
  expect_equal(length(w), p)
  expect_true(all(is.finite(w)))
})

test_that("susieRssWeights retains fit when retainFit = TRUE", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 20
  n <- 500
  R <- diag(p)
  z <- rnorm(p)
  stat <- list(b = z / sqrt(n), cor = z / sqrt(n), z = z, n = rep(n, p))
  w <- susieRssWeights(stat, R, retainFit = TRUE, methodArgs = list(L = 5))
  expect_false(is.null(attr(w, "fit")))
})

test_that("susieInfRssWeights works", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 20
  n <- 500
  R <- diag(p)
  z <- rnorm(p)
  stat <- list(b = z / sqrt(n), cor = z / sqrt(n), z = z, n = rep(n, p))
  w <- susieInfRssWeights(stat, R, methodArgs = list(L = 5))
  expect_equal(length(w), p)
  expect_true(all(is.finite(w)))
})

# =============================================================================
# Two-stage SuSiE-RSS fitting
# =============================================================================

test_that("fitSusieInfThenSusieRss returns two fits", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 20
  n <- 500
  R <- diag(p)
  z <- rnorm(p)
  fits <- fitSusieInfThenSusieRss(z, R, n, args = list(L = 5))
  expect_true(is.list(fits))
  expect_true("susie" %in% names(fits))
  expect_true("susieInf" %in% names(fits))
  expect_true("susieInf" %in% class(fits$susieInf))
  expect_true("susieRss" %in% class(fits$susie))
})
