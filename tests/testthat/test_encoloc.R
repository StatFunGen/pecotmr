context("encoloc")

# The file-path colocalization / enrichment wrappers (colocWrapper,
# xqtlEnrichmentWrapper, colocPostProcessor) and their helpers
# (filterAndOrderColocResults, calculateCumsum, calculate_purity,
# processColocResults, extract_ld_for_variants) have been removed in
# favor of the S4 colocPipeline / qtlEnrichmentPipeline / enlocPipeline
# entry points. The previous tests against the file-path wrappers no
# longer apply; they relied on mocking rssAnalysisPipeline (also
# removed). New tests for colocPipeline / qtlEnrichmentPipeline /
# enlocPipeline live alongside their pipeline implementations.

test_that("xqtlEnrichmentWrapper is a deprecated no-op", {
  expect_warning(
    res <- xqtlEnrichmentWrapper(),
    "has been removed",
    ignore.case = TRUE)
  expect_null(res)
})

test_that("colocWrapper is a deprecated no-op", {
  expect_warning(
    res <- colocWrapper(),
    "has been removed",
    ignore.case = TRUE)
  expect_null(res)
})

test_that("colocPostProcessor is a deprecated no-op", {
  expect_warning(
    res <- colocPostProcessor(),
    "has been removed",
    ignore.case = TRUE)
  expect_null(res)
})
