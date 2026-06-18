context("sumstats_qc")

# Previous tests against `rssBasicQc()` and the legacy
# `summaryStatsQc(rssInput, ldData)` signature have been removed because:
#   * `rssBasicQc()` was folded into `summaryStatsQc()`
#   * `summaryStatsQc()` now dispatches on `GwasSumStats` / `QtlSumStats`
#     S4 collections (DFrame subclasses), not on a (rssInput list, ldData)
#     pair
#   * `QcResult` and the accessors `getRssInput()`, `getLdData()`,
#     `getOutlierNumber()` have been removed; QC audit lives on
#     `getQcInfo(<SumStats>)`
# Integration coverage of `summaryStatsQc(<SumStats>)` lives in the
# pipeline test files (test_colocboostPipeline.R, etc.).
#
# What remains here: tests of internal QC helpers
# (`ldMismatchQc`, `.resolveZMismatchQc`, `krigingOutlierQc`) whose
# signatures and contracts are unchanged.

# ===========================================================================
# ldMismatchQc
# ===========================================================================

test_that("ldMismatchQc with dentist method returns data frame with outlier column", {
  set.seed(42)
  p <- 20
  R <- diag(p)
  z <- rnorm(p)
  result <- ldMismatchQc(z, R = R, nSample = 1000, method = "dentist")
  expect_true(is.data.frame(result) || is.list(result))
  expect_true("outlier" %in% names(result))
})

test_that("ldMismatchQc with slalom method returns data frame with outlier column", {
  set.seed(42)
  p <- 20
  R <- diag(p)
  z <- rnorm(p)
  result <- ldMismatchQc(z, R = R, method = "slalom")
  expect_true(is.data.frame(result) || is.list(result))
  expect_true("outlier" %in% names(result))
})

test_that("ldMismatchQc method argument is validated", {
  z <- rnorm(5)
  R <- diag(5)
  expect_error(ldMismatchQc(z, R = R, method = "invalid"))
})

# ===========================================================================
# zMismatchQc resolver
# ===========================================================================

test_that(".resolveZMismatchQc resolves none/slalom/dentist and defaults to none", {
  expect_equal(pecotmr:::.resolveZMismatchQc(NULL), "none")
  expect_equal(pecotmr:::.resolveZMismatchQc("none"), "none")
  expect_equal(pecotmr:::.resolveZMismatchQc("slalom"), "slalom")
  expect_equal(pecotmr:::.resolveZMismatchQc("dentist"), "dentist")
})

test_that(".resolveZMismatchQc rejects stale rss_qc and other invalid tokens", {
  expect_error(pecotmr:::.resolveZMismatchQc("rss_qc"), "should be one of")
  expect_error(pecotmr:::.resolveZMismatchQc("bad"), "should be one of")
})

# ===========================================================================
# krigingOutlierQc
# ===========================================================================

test_that("krigingOutlierQc flags an LD-inconsistent variant and spares the rest", {
  m <- 6
  rho <- 0.7
  R <- matrix(rho, m, m); diag(R) <- 1
  ids <- paste0("1:", seq_len(m) * 100, ":A:G")
  rownames(R) <- colnames(R) <- ids
  z <- rep(3, m)
  z[3] <- -8                       # strongly inconsistent with its neighbours
  kr <- krigingOutlierQc(z, R, variantIds = ids)
  expect_true(kr$outlier[3])
  expect_false(any(kr$outlier[-3]))
  expect_equal(nrow(kr$diagnostics), m)
  expect_true(all(c("predicted", "residual", "statistic", "p_value") %in%
                    colnames(kr$diagnostics)))
})
