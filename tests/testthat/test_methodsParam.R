context("univariate_pipeline_methods_param")

# univariateAnalysisPipeline has been removed in the S4 refactor. The
# fine-mapping / TWAS-weights pipelines are now invoked via
# fineMappingPipeline(qtlDataset, methods = "susie") and
# twasWeightsPipeline(qtlDataset, ...). Method-selection coverage
# (susie / susieInf / susieAsh dispatch, chained init, etc.) belongs
# alongside those pipelines.

test_that("univariateAnalysisPipeline is a deprecated no-op", {
  expect_warning(
    res <- univariateAnalysisPipeline(),
    "has been removed",
    ignore.case = TRUE)
  expect_null(res)
})
