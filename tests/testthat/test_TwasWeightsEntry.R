# Tests for R/TwasWeightsEntry.R

# === Tests migrated from test_s4Constructors.R (TwasWeightsEntry) ===

test_that("TwasWeightsEntry: constructor and accessors round-trip", {
  e <- TwasWeightsEntry(
    variantIds    = c("v1", "v2", "v3"),
    weights       = c(0.1, -0.2, 0.05),
    fits          = list(model = "lasso"),
    cvPerformance = list(rsq = 0.4),
    standardized  = TRUE,
    dataType      = "expression")
  expect_s4_class(e, "TwasWeightsEntry")
  expect_equal(getVariantIds(e), c("v1", "v2", "v3"))
  expect_equal(getWeights(e), c(0.1, -0.2, 0.05))
  expect_equal(getFits(e), list(model = "lasso"))
  expect_equal(getCvPerformance(e), list(rsq = 0.4))
  expect_true(getStandardized(e))
  expect_equal(getDataType(e), "expression")
})


test_that("TwasWeightsEntry: standardized is coerced via isTRUE() semantics", {
  # isTRUE() only returns TRUE for a length-1 logical TRUE. Non-TRUE
  # input lands as FALSE (the safe default for the standardized flag).
  e_logical <- TwasWeightsEntry(
    variantIds = "v1", weights = 0.1, standardized = TRUE)
  expect_true(getStandardized(e_logical))

  e_default <- TwasWeightsEntry(
    variantIds = "v1", weights = 0.1, standardized = "yes-please")
  expect_false(getStandardized(e_default))
})


test_that("TwasWeightsEntry: validity rejects matrix weights with wrong nrow", {
  expect_error(
    TwasWeightsEntry(
      variantIds = c("v1", "v2"),
      weights    = matrix(0, nrow = 5, ncol = 1)),
    "nrow\\(weights\\) must equal length\\(variantIds\\)"
  )
})

# ===========================================================================
# QtlFineMappingResult
# ===========================================================================


test_that("TwasWeights: rejects non-TwasWeightsEntry rows", {
  expect_error(
    TwasWeights(
      study = "s1", context = "c1", trait = "t1", method = "lasso",
      entry = list("not_an_entry")),
    "every element of the `entry` column must be a TwasWeightsEntry"
  )
})



# === Tests migrated from test_showMethods.R (TwasWeightsEntry) ===

test_that("show.TwasWeightsEntry reports standardized flag and CV availability", {
  e <- .sh_makeTwEntry(p = 5, standardized = TRUE)
  out <- capture.output(show(e))
  expect_true(any(grepl("TwasWeightsEntry: 5 variants.*standardized=TRUE", out)))
  expect_true(any(grepl("CV performance: TRUE", out)))

  e_no_cv <- TwasWeightsEntry(variantIds = c("v1", "v2"),
                               weights = c(0.1, 0.2))
  out2 <- capture.output(show(e_no_cv))
  expect_true(any(grepl("CV performance: FALSE", out2)))
})


