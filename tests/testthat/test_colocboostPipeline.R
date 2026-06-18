context("colocboostPipeline")

# ============================================================================
# Legacy colocboostPipeline tests removed for the post-S4-refactor API.
# ============================================================================
#
# This file previously contained ~100 test_that blocks targeting the legacy
# colocboost pipeline API (rssData/ldData list-shape inputs, RegionalData
# fixtures, qcRegionalData driver, colocboostAnalysis adapters, etc.). Every
# block exercised functions and classes that no longer exist in pecotmr:
#
#   Functions removed:
#     - colocboostAnalysis()         (direct colocboost adapter)
#     - qcRegionalData()             (legacy QC driver)
#     - regionDataToIndInput()       (now a .Deprecated() no-op)
#     - regionDataToRssInput()       (now a .Deprecated() no-op)
#     - regionDataToColocboostInput()
#     - rssAnalysisPipeline()        (replaced by fineMappingPipeline)
#     - rssBasicQc()                 (folded into summaryStatsQc)
#     - loadRssData()                (replaced by SumStats constructors)
#     - getrssinput(), getlddata(), getoutliernumber()
#     - colocWrapper(), xqtlEnrichmentWrapper(), colocPostProcessor()
#       (now .Deprecated() no-ops returning NULL)
#     - .runColocboost()             (replaced by internal .cbRun(label, args))
#     - buildLdArgs()                (replaced by internal .cbBuildLdArgs())
#
#   Classes removed:
#     - RegionalData, MultivariateRegionalData
#     - QcResult, AlleleQcResult
#
# The replacement API has a fundamentally different contract:
#
#   colocboostPipeline() is now an S4 generic dispatching on the QTL input
#   class. Signatures live in R/colocboostPipeline.R:
#
#     setMethod("colocboostPipeline", "QtlDataset",         ...)
#     setMethod("colocboostPipeline", "QtlSumStats",        ...)
#     setMethod("colocboostPipeline", "MultiTaskQtlDataset", ...)
#
#   - QTL data is supplied as a QtlDataset / QtlSumStats / MultiTaskQtlDataset
#     (DFrame-based S4 objects with an ldSketch slot for sumstats inputs and
#     getResidualizedGenotypes() / getResidualizedPhenotypes() accessors for
#     individual-level inputs).
#   - GWAS is supplied separately via the gwasSumStats = GwasSumStats(...)
#     argument.
#   - Individual-level QC (MAF / X-variance / sample missingness / event
#     selection) lives on the QtlDataset constructor and is applied lazily
#     by its accessors. There is no separate qcRegionalData() pass.
#   - All summary-statistic QC lives in summaryStatsQc(). The pipeline
#     rejects QtlSumStats / GwasSumStats whose getQcInfo() is empty.
#
# Rewriting the legacy tests in place would require fabricating new
# QtlDataset / QtlSumStats / GwasSumStats / MultiTaskQtlDataset fixtures
# and asserting against a different result shape -- i.e. inventing new
# coverage rather than porting existing coverage. That is out of scope
# for this legacy-cleanup pass.
#
# New tests for the S4 colocboostPipeline() methods, summaryStatsQc(), and
# the QtlDataset / QtlSumStats / GwasSumStats / MultiTaskQtlDataset
# constructors should be added in dedicated files alongside this one
# (e.g. test_colocboost_pipeline_qtl_dataset.R,
# test_colocboost_pipeline_qtl_sumstats.R,
# test_colocboost_pipeline_multitask.R). See:
#
#   - R/colocboostPipeline.R    (the new S4 generic + methods)
#   - R/allClasses.R             (QtlDataset, QtlSumStats, GwasSumStats,
#                                 MultiTaskQtlDataset, FineMappingEntry)
#   - R/allMethods.R             (constructors)
#   - R/sumstatsQc.R            (summaryStatsQc, which is now the only
#                                 summary-statistic QC entry point)
#   - tests/testthat/test_sumstatsQc.R  (existing QC coverage)
# ============================================================================

# Sentinel test so the testthat context is non-empty.
test_that("colocboostPipeline is exported as an S4 generic", {
  expect_true("colocboostPipeline" %in% getNamespaceExports("pecotmr"))
  expect_true(methods::isGeneric("colocboostPipeline"))
})
