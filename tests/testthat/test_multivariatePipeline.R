context("multivariate_pipeline")

# multivariateAnalysisPipeline and the legacy regionDataToMvsusieRssInput /
# regionDataToSusieRssInput entry points have been removed in the S4
# refactor. Joint multi-trait/multi-context analyses now run through
# fineMappingPipeline(qtlDataset, methods = "mvsusie") (individual-level)
# or fineMappingPipeline(..., methods = "mvsusieRSS") (RSS-based).
# Coverage for those pipelines lives alongside their implementations.

test_that("multivariateAnalysisPipeline is a deprecated no-op", {
  expect_warning(
    res <- multivariateAnalysisPipeline(),
    "has been removed",
    ignore.case = TRUE)
  expect_null(res)
})
