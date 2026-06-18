context("univariate_pipeline")

# ===========================================================================
# Post-S4-refactor cleanup
# ---------------------------------------------------------------------------
# The S4 refactor removed `univariateAnalysisPipeline()`,
# `rssAnalysisPipeline()`, `loadStudyLd()`, `loadRssData()`, and
# `rssBasicQc()` as functional pipelines. They are now `.Deprecated()`
# no-ops that return `NULL` invisibly. The `QcResult` /
# `FineMappingResult` classes and `regionDataToSusieRssInput()` helper
# were also removed. Their replacements live behind a different S4 API
# (`fineMappingPipeline()` dispatched on `GwasSumStats` / `QtlSumStats` /
# `QtlDataset`, `summaryStatsQc()`, `FineMappingEntry()`) with different
# signatures and contracts, so the legacy mocked pipeline tests cannot
# be ported in a meaningful way and have been removed.
#
# What survives in this file:
#   * Deprecation no-op checks for the removed public wrappers.
#   * Tests for `resolveLdInput()`, which is still a live internal helper
#     in `dentistQc.R` with an unchanged signature.
# ===========================================================================

# ===========================================================================
# resolveLdInput (still-live internal helper in dentistQc.R)
# ===========================================================================

test_that("resolveLdInput errors when both R and X are NULL", {
  expect_error(
    pecotmr:::resolveLdInput(R = NULL, X = NULL),
    "Either R .* or X .* must be provided"
  )
})

test_that("resolveLdInput errors when both R and X are provided", {
  R <- diag(5)
  X <- matrix(rnorm(50), 10, 5)
  expect_error(
    pecotmr:::resolveLdInput(R = R, X = X),
    "Provide either R or X, not both"
  )
})

test_that("resolveLdInput with R returns R unchanged", {
  R <- diag(5)
  result <- pecotmr:::resolveLdInput(R = R, nSample = 100)
  expect_equal(result$R, R)
  expect_equal(result$nSample, 100)
})

test_that("resolveLdInput with X computes LD and infers nSample", {
  set.seed(42)
  X <- matrix(rbinom(200, 2, 0.3), 20, 10)
  result <- pecotmr:::resolveLdInput(X = X)
  expect_equal(nrow(result$R), 10)
  expect_equal(ncol(result$R), 10)
  expect_equal(result$nSample, 20)
})

test_that("resolveLdInput errors when nSample required but missing", {
  R <- diag(5)
  expect_error(
    pecotmr:::resolveLdInput(R = R, needNSample = TRUE),
    "nSample is required"
  )
})

test_that("resolveLdInput does not error when nSample not needed", {
  R <- diag(5)
  result <- pecotmr:::resolveLdInput(R = R, needNSample = FALSE)
  expect_equal(result$R, R)
  expect_null(result$nSample)
})

# ===========================================================================
# Deprecation no-op checks for removed public wrappers
# ===========================================================================

test_that("univariateAnalysisPipeline is a deprecated no-op", {
  expect_warning(res <- univariateAnalysisPipeline(), "deprecated|removed")
  expect_null(res)
})

test_that("rssAnalysisPipeline is a deprecated no-op", {
  expect_warning(res <- rssAnalysisPipeline(), "deprecated|removed")
  expect_null(res)
})

test_that("loadStudyLd is a deprecated no-op", {
  expect_warning(res <- loadStudyLd(), "deprecated|removed")
  expect_null(res)
})

test_that("loadRssData is a deprecated no-op", {
  expect_warning(res <- loadRssData(), "deprecated|removed")
  expect_null(res)
})

# ===========================================================================
# Removed during the post-S4-refactor cleanup (for traceability)
# ---------------------------------------------------------------------------
# univariateAnalysisPipeline input-validation tests (12):
#   * X non-matrix / non-numeric / Y non-numeric / multi-column Y /
#     single-column Y / X-Y row mismatch / maf length mismatch /
#     maf out of bounds / invalid xScalar / wrong-length xScalar /
#     non-numeric yScalar / non-positive L / non-positive lGreedy.
# univariateAnalysisPipeline functional tests (4):
#   * minimal valid input; twasWeights = TRUE; respects xScalar;
#     cvFolds = 0.
# univariateAnalysisPipeline mocked-pipeline tests (16):
#   * pip_cutoff_to_skip branch (3); LD reference filtering (2);
#     filter_X branch (2); main susie + post-processing (4);
#     TWAS weights (4); coverage forwarding (2); combined LD +
#     filter_X (1); null imiss/maf cutoffs (1).
# univariateAnalysisPipeline af/maf single-source tests (6):
#   * af supplied (no maf); maf supplied (no af); both disagreeing;
#     neither with mafCutoff errors; neither without mafCutoff runs;
#     af-derived MAF parity with explicit maf.
# rssAnalysisPipeline tests (~21):
#   * input-file validation; empty sumstats branches; pip_cutoff_to_skip
#     branches; full QC + imputation + fine-mapping; method-name
#     conventions (NO_QC / SLALOM / DENTIST / RAISS); outlierNumber
#     plumbing; finemappingMethod = NULL skip; zMismatchQc = NULL;
#     impute = FALSE; diagnostics branches (BCR + SER re-analysis);
#     finemappingOpts passthrough; mafCutoff (af single source);
#     finemapping-defaults passthrough; X-as-LD genotype path;
#     mixture-LD passthrough.
# regionDataToSusieRssInput tests (3):
#   * correlation-backed input; genotype-backed input; rejects
#     default rownames as variant IDs.
# loadStudyLd tests (2):
#   * single path returns loadLdMatrix output unchanged; comma-separated
#     paths build a mixture LdData with a list of genotype handles.
#
# Replacement APIs (do NOT speculatively port — that is out of scope for
# this cleanup): `fineMappingPipeline()` dispatched on `GwasSumStats` /
# `QtlSumStats` / `QtlDataset`; `summaryStatsQc()` (returns a SumStats
# with `getQcInfo()` populated); `FineMappingEntry(variantIds, ...)`;
# `result$finemappingEntry` (was `result$finemappingResult`).
# ===========================================================================
