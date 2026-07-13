context("mash_wrapper")

# Build a minimal FineMappingEntry for unit-testing find_nested /
# extractFlattenSumstatsFromNested. Note: the legacy FineMappingResult
# constructor (with `variantNames`/`method` args) has been removed; the new
# per-entry payload class is `FineMappingEntry`, and method identity now
# lives on the parent `FineMappingResult` collection row.
.testFineMappingEntry <- function(variantNames) {
    FineMappingEntry(
        variantIds = variantNames,
        susieFit = list(pip = rep(0.5, length(variantNames))),
        topLoci = data.frame(variant_id = character(0),
                              pip = numeric(0),
                              stringsAsFactors = FALSE)
    )
}

# ===========================================================================
# filterInvalidSummaryStat
# ===========================================================================

test_that("filterInvalidSummaryStat replaces NaN/Inf in bhat", {
  dat <- list(
    bhat = data.frame(a = c(1, NaN, 3), b = c(Inf, 2, -Inf)),
    sbhat = data.frame(a = c(0.1, 0.2, 0.3), b = c(0.1, NA, 0.3))
  )
  result <- filterInvalidSummaryStat(dat, bhat = "bhat", sbhat = "sbhat")
  expect_true(all(!is.nan(result$bhat)))
  expect_true(all(!is.infinite(result$bhat)))
  # NaN/Inf in bhat replaced with 0
  expect_equal(unname(result$bhat[1, 2]), 0)  # Inf -> 0
})

test_that("filterInvalidSummaryStat replaces NaN/Inf in sbhat", {
  dat <- list(
    bhat = data.frame(a = c(1, 2, 3)),
    sbhat = data.frame(a = c(0.1, NaN, Inf))
  )
  result <- filterInvalidSummaryStat(dat, bhat = "bhat", sbhat = "sbhat")
  # NaN/Inf in sbhat replaced with 1000
  expect_equal(unname(result$sbhat[1, "a"]), 0.1)
  expect_equal(unname(result$sbhat[2, "a"]), 1000)
  expect_equal(unname(result$sbhat[3, "a"]), 1000)
})

test_that("filterInvalidSummaryStat filters by missing_rate when null.b present", {
  dat <- list(
    bhat = data.frame(a = c(0, 0, 1, 2), b = c(0, 0, 0, 3)),
    sbhat = data.frame(a = c(1, 1, 1, 1), b = c(1, 1, 1, 1)),
    null.b = TRUE
  )
  result <- filterInvalidSummaryStat(dat, bhat = "bhat", sbhat = "sbhat",
                                        filterByMissingRate = 0.5)
  expect_equal(nrow(result$bhat), 2) # rows 3 and 4 survive
})

test_that("filterInvalidSummaryStat filters by missing_rate when random.b present", {
  dat <- list(
    bhat = data.frame(a = c(0, 1, 2), b = c(0, 1, 3)),
    sbhat = data.frame(a = c(1, 1, 1), b = c(1, 1, 1)),
    random.b = TRUE
  )
  result <- filterInvalidSummaryStat(dat, bhat = "bhat", sbhat = "sbhat",
                                        filterByMissingRate = 0.5)
  expect_true(nrow(result$bhat) < 3)
})

test_that("filterInvalidSummaryStat btoz with .b and .s pattern creates condition.z", {
  dat <- list(
    strong.b = data.frame(a = c(1, 2), b = c(3, 4)),
    strong.s = data.frame(a = c(0.1, 0.2), b = c(0.3, 0.4))
  )
  result <- filterInvalidSummaryStat(dat, bhat = "strong.b", sbhat = "strong.s",
                                        btoz = TRUE, sigPCutoff = NULL)
  expect_true("strong.z" %in% names(result))
  expect_true(is.matrix(result$strong.z))
  expect_equal(nrow(result$strong.z), 2)
  expect_equal(ncol(result$strong.z), 2)
  expect_equal(as.numeric(result$strong.z), c(10, 10, 10, 10), tolerance = 1e-10)
})

test_that("filterInvalidSummaryStat btoz when bhat/sbhat data is NULL creates NULL z", {
  dat <- list(
    strong.b = NULL,
    strong.s = data.frame(a = c(0.1))
  )
  result <- filterInvalidSummaryStat(dat, bhat = "strong.b", sbhat = "strong.s",
                                        btoz = TRUE, sigPCutoff = NULL)
  expect_true("strong.z" %in% names(result))
  expect_null(result$strong.z)
})

test_that("filterInvalidSummaryStat btoz without .b/.s pattern creates generic z", {
  dat <- list(
    bhat = data.frame(a = c(1, 2, 3)),
    sbhat = data.frame(a = c(0.5, 1, 0.5))
  )
  result <- filterInvalidSummaryStat(dat, bhat = "bhat", sbhat = "sbhat",
                                        btoz = TRUE)
  expect_true("z" %in% names(result))
  expect_equal(as.numeric(result$z[, 1]), c(2, 2, 6))
})

test_that("filterInvalidSummaryStat btoz creates NULL z when bhat is NULL (no .b suffix)", {
  dat <- list(
    bhat = NULL,
    sbhat = data.frame(a = c(0.1))
  )
  result <- filterInvalidSummaryStat(dat, bhat = "bhat", sbhat = "sbhat",
                                        btoz = TRUE)
  expect_true("z" %in% names(result))
  expect_null(result$z)
})

test_that("filterInvalidSummaryStat btoz filters strong.z by significance cutoff", {
  dat <- list(
    strong.b = data.frame(a = c(1, 0.01, 0.02, 2), b = c(0.01, 0.01, 0.01, 0.01)),
    strong.s = data.frame(a = c(0.1, 0.1, 0.1, 0.1), b = c(0.1, 0.1, 0.1, 0.1)),
    strong.z = NULL
  )
  result <- filterInvalidSummaryStat(dat, bhat = "strong.b", sbhat = "strong.s",
                                        btoz = TRUE, sigPCutoff = 1E-6)
  expect_true("strong.z" %in% names(result))
  expect_equal(nrow(result$strong.z), 2)
  expect_equal(nrow(result$strong.b), 2)
  expect_equal(nrow(result$strong.s), 2)
})

test_that("filterInvalidSummaryStat processes z directly with null component", {
  dat <- list(
    strong = list(z = data.frame(a = c(5, NaN, 0.1), b = c(1, 2, Inf))),
    random = list(z = data.frame(a = c(0.5, 0.2), b = c(0.3, 0.4))),
    null = list(z = data.frame(a = c(0.01, NaN), b = c(Inf, 0.02)))
  )
  result <- filterInvalidSummaryStat(dat, z = "z")
  expect_true(all(!is.nan(result$strong$z)))
  expect_true(all(!is.infinite(result$strong$z)))
  expect_true(all(!is.nan(result$null$z)))
  expect_true(all(!is.infinite(result$null$z)))
})

test_that("filterInvalidSummaryStat z path applies significance cutoff to strong.z", {
  dat <- list(
    strong = list(z = data.frame(a = c(10, 0.1, 0.2), b = c(0.1, 0.1, 0.1))),
    random = list(z = data.frame(a = c(0.5, 0.2, 0.3), b = c(0.3, 0.4, 0.2)))
  )
  result <- filterInvalidSummaryStat(dat, z = "z", sigPCutoff = 1E-6)
  expect_equal(nrow(result$strong$z), 1)
})

test_that("filterInvalidSummaryStat z path with filterByMissingRate", {
  dat <- list(
    random = list(z = data.frame(a = c(0, 0, 5), b = c(0, 3, 4)))
  )
  result <- filterInvalidSummaryStat(dat, z = "z", filterByMissingRate = 0.5)
  expect_equal(nrow(result$random$z), 2)
})

test_that("filterInvalidSummaryStat processes bhat/sbhat without filterByMissingRate when no null.b/random.b", {
  dat <- list(
    bhat = data.frame(a = c(1, NaN, 3), b = c(Inf, 2, -Inf)),
    sbhat = data.frame(a = c(0.1, NA, 0.3), b = c(0.1, 0.2, NaN))
  )
  result <- filterInvalidSummaryStat(dat, bhat = "bhat", sbhat = "sbhat",
                                        filterByMissingRate = 0.5)
  expect_equal(nrow(result$bhat), 3)
  expect_equal(unname(result$bhat[2, "a"]), 0)
  expect_equal(unname(result$sbhat[3, "b"]), 1000)
})

test_that("filterInvalidSummaryStat with filterByMissingRate=NULL keeps all rows even with null.b", {
  dat <- list(
    bhat = data.frame(a = c(0, 0, 1), b = c(0, 0, 1)),
    sbhat = data.frame(a = c(1, 1, 1), b = c(1, 1, 1)),
    null.b = TRUE
  )
  result <- filterInvalidSummaryStat(dat, bhat = "bhat", sbhat = "sbhat",
                                        filterByMissingRate = NULL)
  expect_equal(nrow(result$bhat), 3)
})

test_that("filterInvalidSummaryStat z path handles NULL strong component", {
  dat <- list(
    strong = NULL,
    random = list(z = data.frame(a = c(0.5, 0.2), b = c(0.3, 0.4)))
  )
  result <- filterInvalidSummaryStat(dat, z = "z")
  expect_null(result$strong)
  expect_true(!is.null(result$random$z))
})

# ===========================================================================
# filterMixtureComponents
# ===========================================================================

test_that("filterMixtureComponents removes zero matrices", {
  U <- list(
    comp1 = matrix(c(1, 0, 0, 2), 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
    comp2 = matrix(c(0, 0, 0, 0), 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
    A = matrix(c(3, 0, 0, 4), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  )
  w <- c(comp1 = 0.5, comp2 = 0.3, A = 0.2)
  result <- filterMixtureComponents(c("A", "B"), U, w)
  expect_false("comp2" %in% names(result$U))
})

test_that("filterMixtureComponents removes matrices below weight cutoff", {
  U <- list(
    comp1 = matrix(c(1, 0, 0, 2), 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
    comp2 = matrix(c(0.1, 0, 0, 0.1), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  )
  w <- c(comp1 = 0.9, comp2 = 0.00001)
  result <- filterMixtureComponents(c("A", "B"), U, w, wCutoff = 1e-4)
  expect_false("comp2" %in% names(result$U))
})

test_that("filterMixtureComponents errors on missing conditions", {
  U <- list(
    comp1 = matrix(c(1, 0, 0, 2), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  )
  expect_error(
    filterMixtureComponents(c("A", "C"), U),
    "not found in matrix"
  )
})

test_that("filterMixtureComponents removes components named as filtered conditions", {
  U <- list(
    comp1 = matrix(c(1,0,0, 0,2,0, 0,0,3), 3, 3, dimnames = list(c("A","B","C"), c("A","B","C"))),
    C = matrix(c(4,0,0, 0,5,0, 0,0,6), 3, 3, dimnames = list(c("A","B","C"), c("A","B","C")))
  )
  w <- c(comp1 = 0.6, C = 0.4)
  result <- filterMixtureComponents(c("A", "B"), U, w)
  expect_false("C" %in% names(result$U))
  expect_equal(nrow(result$U$comp1), 2)
  expect_equal(ncol(result$U$comp1), 2)
})

test_that("filterMixtureComponents renormalizes weights to preserve original sum", {
  U <- list(
    comp1 = matrix(c(1, 0, 0, 2), 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
    comp2 = matrix(c(3, 0, 0, 4), 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
    comp3 = matrix(c(0, 0, 0, 0), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  )
  w <- c(comp1 = 0.5, comp2 = 0.3, comp3 = 0.2)
  original_sum <- sum(w)

  result <- filterMixtureComponents(c("A", "B"), U, w)
  expect_false("comp3" %in% names(result$U))
  expect_equal(sum(result$w), original_sum, tolerance = 1e-10)
  expect_true(result$w["comp1"] > 0.5)
  expect_true(result$w["comp2"] > 0.3)
})

test_that("filterMixtureComponents subsets 3x3 matrices to 2x2 and removes filtered condition names", {
  U <- list(
    comp1 = matrix(c(1, 0.1, 0, 0.1, 2, 0, 0, 0, 3), 3, 3,
                   dimnames = list(c("A", "B", "C"), c("A", "B", "C"))),
    B = matrix(c(4, 0, 0, 0, 5, 0, 0, 0, 6), 3, 3,
               dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  )
  w <- c(comp1 = 0.7, B = 0.3)

  result <- filterMixtureComponents(c("A", "C"), U, w)
  expect_false("B" %in% names(result$U))
  expect_equal(nrow(result$U$comp1), 2)
  expect_equal(ncol(result$U$comp1), 2)
  expect_equal(rownames(result$U$comp1), c("A", "C"))
})

test_that("filterMixtureComponents handles NULL weights gracefully", {
  U <- list(
    comp1 = matrix(c(1, 0, 0, 2), 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
    comp2 = matrix(c(0, 0, 0, 0), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  )
  result <- filterMixtureComponents(c("A", "B"), U, w = NULL)
  expect_false("comp2" %in% names(result$U))
  expect_true("comp1" %in% names(result$U))
})

# ===========================================================================
# mergeMashData
# ===========================================================================

test_that("mergeMashData combines two datasets with identical columns", {
  d1 <- list(random = data.frame(a = 1:3, b = 4:6))
  d2 <- list(random = data.frame(a = 7:8, b = 9:10))
  result <- mergeMashData(d1, d2)
  expect_equal(nrow(result$random), 5)
  expect_equal(ncol(result$random), 2)
  expect_equal(colnames(result$random), c("a", "b"))
  expect_equal(result$random$a, c(1, 2, 3, 7, 8))
})

test_that("mergeMashData handles NULL input", {
  d1 <- NULL
  d2 <- list(random = data.frame(a = 1:3))
  result <- mergeMashData(d1, d2)
  expect_equal(nrow(result$random), 3)
})

test_that("mergeMashData aligns different column names correctly", {
  d1 <- list(random = data.frame(a = 1:2, b = 3:4))
  d2 <- list(random = data.frame(a = 5:6, c = 7:8))
  result <- mergeMashData(d1, d2)
  expect_equal(nrow(result$random), 4)
  expect_true(all(c("a", "b", "c") %in% colnames(result$random)))
  expect_true(is.nan(result$random[3, "b"]))
  expect_true(is.nan(result$random[1, "c"]))
})

test_that("mergeMashData preserves data when one side is empty data.frame", {
  d1 <- list(random = data.frame(a = 1:3))
  d2 <- list(random = data.frame())
  result <- mergeMashData(d1, d2)
  expect_equal(nrow(result$random), 3)
})

test_that("mergeMashData preserves data when one side is NULL element", {
  d1 <- list(random = NULL)
  d2 <- list(random = data.frame(a = 1:3))
  result <- mergeMashData(d1, d2)
  expect_equal(nrow(result$random), 3)
})

test_that("mergeMashData handles multiple named elements", {
  d1 <- list(
    random = data.frame(a = 1:2, b = 3:4),
    null = data.frame(x = 10:11)
  )
  d2 <- list(
    random = data.frame(a = 5:6, b = 7:8),
    null = data.frame(x = 12:13)
  )
  result <- mergeMashData(d1, d2)
  expect_equal(nrow(result$random), 4)
  expect_equal(nrow(result$null), 4)
})

test_that("mergeMashData uses one_data when res_data element has zero rows", {
  d1 <- list(random = data.frame(a = numeric(0), b = numeric(0)))
  d2 <- list(random = data.frame(a = 1:3, b = 4:6))
  result <- mergeMashData(d1, d2)
  expect_true(nrow(result$random) >= 3)
})

# ===========================================================================
# mashRandNullSample
# ===========================================================================

test_that("mashRandNullSample with z scores returns random and null", {
  set.seed(42)
  dat <- list(
    z = data.frame(
      cond1 = c(5, 0.1, 0.2, 0.3, 0.5, 6, 0.1, 0.2, 0.4, 0.3),
      cond2 = c(0.2, 0.3, 0.1, 0.5, 0.4, 0.1, 0.3, 0.2, 0.1, 0.5)
    )
  )
  result <- mashRandNullSample(dat, nRandom = 5, nNull = 3,
                                   excludeCondition = c(), seed = 123)
  expect_type(result, "list")
  expect_true("random" %in% names(result))
  expect_true("null" %in% names(result))
  expect_true("z" %in% names(result$random))
  expect_equal(nrow(result$random$z), 5)
})

test_that("mashRandNullSample with seed is reproducible", {
  dat <- list(
    z = data.frame(
      cond1 = c(0.1, 0.2, 0.3, 0.4, 0.5),
      cond2 = c(0.5, 0.4, 0.3, 0.2, 0.1)
    )
  )
  result1 <- mashRandNullSample(dat, nRandom = 3, nNull = 2,
                                    excludeCondition = c(), seed = 42)
  result2 <- mashRandNullSample(dat, nRandom = 3, nNull = 2,
                                    excludeCondition = c(), seed = 42)
  expect_equal(result1$random$z, result2$random$z)
})

test_that("mashRandNullSample NULL input returns NULL", {
  result <- mashRandNullSample(NULL, nRandom = 5, nNull = 3,
                                   excludeCondition = c())
  expect_null(result)
})

test_that("mashRandNullSample warns when no null variants found (all abs_z > 2)", {
  dat <- list(
    z = data.frame(
      cond1 = c(5, 6, 7, 8, 9),
      cond2 = c(5, 6, 7, 8, 9)
    )
  )
  expect_warning(
    result <- mashRandNullSample(dat, nRandom = 3, nNull = 2,
                                     excludeCondition = c(), seed = 42),
    "no variants are included in the null"
  )
  expect_equal(length(result$null), 0)
})

test_that("mashRandNullSample warns when not enough null data", {
  dat <- list(
    z = data.frame(
      cond1 = c(5, 6, 0.1),
      cond2 = c(5, 6, 0.1),
      cond3 = c(5, 6, 0.1)
    )
  )
  expect_warning(
    result <- mashRandNullSample(dat, nRandom = 2, nNull = 1,
                                     excludeCondition = c(), seed = 42),
    "not enough null data"
  )
  expect_equal(length(result$null), 0)
})

test_that("mashRandNullSample with bhat/sbhat processes random and null samples", {
  dat <- list(
    bhat = data.frame(
      cond1 = c(0.1, 0.05, 0.02, 0.01, 0.03),
      cond2 = c(0.02, 0.01, 0.03, 0.05, 0.04)
    ),
    sbhat = data.frame(
      cond1 = c(0.1, 0.1, 0.1, 0.1, 0.1),
      cond2 = c(0.1, 0.1, 0.1, 0.1, 0.1)
    )
  )
  result <- mashRandNullSample(dat, nRandom = 3, nNull = 3,
                                   excludeCondition = c(), seed = 42)
  expect_equal(ncol(result$random$bhat), 2)
  expect_equal(nrow(result$random$bhat), 3)
  expect_true(length(result$null) > 0)
})

test_that("mashRandNullSample errors when excludeCondition not found (z path)", {
  dat <- list(
    z = data.frame(cond1 = 1:5, cond2 = 1:5)
  )
  expect_error(
    mashRandNullSample(dat, nRandom = 3, nNull = 2,
                           excludeCondition = "nonexistent", seed = 42),
    "excludeCondition are not present"
  )
})

test_that("mashRandNullSample errors when excludeCondition not found (bhat path)", {
  dat <- list(
    bhat = data.frame(cond1 = 1:5, cond2 = 1:5),
    sbhat = data.frame(cond1 = rep(1, 5), cond2 = rep(1, 5))
  )
  expect_error(
    mashRandNullSample(dat, nRandom = 3, nNull = 2,
                           excludeCondition = "nonexistent", seed = 42),
    "excludeCondition are not present"
  )
})

test_that("mashRandNullSample drops excluded condition by column name", {
  dat <- list(
    z = data.frame(
      cond1 = c(0.1, 0.2, 0.3, 0.4, 0.5),
      cond2 = c(0.5, 0.4, 0.3, 0.2, 0.1),
      cond3 = c(0.3, 0.3, 0.3, 0.3, 0.3)
    )
  )
  result <- mashRandNullSample(dat, nRandom = 3, nNull = 2,
                                   excludeCondition = "cond3", seed = 42)
  expect_equal(colnames(result$random$z), c("cond1", "cond2"))
  expect_equal(colnames(result$null$z), c("cond1", "cond2"))
  expect_false("cond3" %in% colnames(result$random$z))
})

test_that("mashRandNullSample extracts null data with z scores when enough null variants exist", {
  dat <- list(
    z = data.frame(
      cond1 = c(0.1, 0.3, 0.2, 0.5, 0.4, 0.1, 0.3, 0.2, 0.5, 0.4),
      cond2 = c(0.2, 0.1, 0.4, 0.3, 0.5, 0.2, 0.1, 0.4, 0.3, 0.5)
    )
  )
  result <- mashRandNullSample(dat, nRandom = 5, nNull = 4,
                                   excludeCondition = c(), seed = 42)
  expect_true("null" %in% names(result))
  expect_true("z" %in% names(result$null))
  expect_equal(nrow(result$null$z), 4)
  expect_equal(ncol(result$null$z), 2)
  expect_equal(nrow(result$random$z), 5)
})

test_that("mashRandNullSample null data capped at available null variants", {
  dat <- list(
    z = data.frame(
      cond1 = c(0.1, 0.2, 0.3),
      cond2 = c(0.1, 0.2, 0.3)
    )
  )
  result <- mashRandNullSample(dat, nRandom = 2, nNull = 100,
                                   excludeCondition = c(), seed = 42)
  expect_true(length(result$null) > 0)
  expect_equal(nrow(result$null$z), 3)
})

test_that("mashRandNullSample extracts null data with bhat/sbhat when enough null variants", {
  dat <- list(
    bhat = data.frame(
      cond1 = c(0.01, 0.02, 0.01, 0.03, 0.02, 0.01),
      cond2 = c(0.02, 0.01, 0.03, 0.01, 0.02, 0.01)
    ),
    sbhat = data.frame(
      cond1 = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
      cond2 = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
    )
  )
  result <- mashRandNullSample(dat, nRandom = 3, nNull = 4,
                                   excludeCondition = c(), seed = 42)
  expect_true("null" %in% names(result))
  expect_true("bhat" %in% names(result$null))
  expect_true("sbhat" %in% names(result$null))
  expect_equal(nrow(result$null$bhat), 4)
  expect_equal(nrow(result$null$sbhat), 4)
})

test_that("mashRandNullSample excludeCondition with numeric index errors on z path", {
  dat <- list(
    z = data.frame(
      cond1 = c(0.1, 0.2, 0.3, 0.4, 0.5),
      cond2 = c(0.5, 0.4, 0.3, 0.2, 0.1),
      cond3 = c(0.3, 0.3, 0.3, 0.3, 0.3)
    )
  )
  expect_error(
    mashRandNullSample(dat, nRandom = 3, nNull = 2,
                           excludeCondition = 3, seed = 42),
    "excludeCondition are not present"
  )
})

test_that("mashRandNullSample excludeCondition with numeric index errors on bhat path", {
  dat <- list(
    bhat = data.frame(
      cond1 = c(0.01, 0.02, 0.01, 0.03, 0.02),
      cond2 = c(0.02, 0.01, 0.03, 0.01, 0.02),
      cond3 = c(0.01, 0.01, 0.01, 0.01, 0.01)
    ),
    sbhat = data.frame(
      cond1 = c(0.1, 0.1, 0.1, 0.1, 0.1),
      cond2 = c(0.1, 0.1, 0.1, 0.1, 0.1),
      cond3 = c(0.1, 0.1, 0.1, 0.1, 0.1)
    )
  )
  expect_error(
    mashRandNullSample(dat, nRandom = 3, nNull = 2,
                           excludeCondition = 3, seed = 42),
    "excludeCondition are not present"
  )
})

test_that("mashRandNullSample caps random sample at available rows", {
  dat <- list(
    z = data.frame(
      cond1 = c(0.1, 0.2, 0.3),
      cond2 = c(0.2, 0.1, 0.3)
    )
  )
  result <- mashRandNullSample(dat, nRandom = 100, nNull = 2,
                                   excludeCondition = c(), seed = 42)
  expect_equal(nrow(result$random$z), 3)
})

# mashPipeline integration tests were removed in the S4 refactor: the
# legacy matrix-list input (strong.b/strong.s/random.b/random.s/null.b/null.s)
# is no longer accepted. The new API takes a named list of QtlSumStats /
# GwasSumStats objects, which is non-trivial to mock without exercising the
# full SumStats QC pipeline. Cover mashPipeline behavior end-to-end via the
# pipeline-level integration tests instead.



# === Tests migrated from test_mrmashWrapper.R (filterMixtureComponents) ===

test_that("filterMixtureComponents filters zero matrices", {
  U <- list(
    mat1 = matrix(c(1, 0.5, 0.5, 1), 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
    mat2 = matrix(0, 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
    mat3 = matrix(c(0.8, 0.3, 0.3, 0.9), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  )
  w <- c(mat1 = 0.5, mat2 = 0.3, mat3 = 0.2)
  conditions_to_keep <- c("A", "B")

  result <- filterMixtureComponents(conditions_to_keep, U, w)

  # mat2 should be removed (all zeros)
  expect_true(!"mat2" %in% names(result$U))
  # weights should be rescaled to maintain sum
  expect_equal(sum(result$w), sum(w), tolerance = 1e-10)
})


test_that("filterMixtureComponents removes low weight components", {
  U <- list(
    mat1 = matrix(c(1, 0.5, 0.5, 1), 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
    mat2 = matrix(c(0.8, 0.3, 0.3, 0.9), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  )
  w <- c(mat1 = 0.999, mat2 = 0.00001)  # mat2 below default cutoff
  conditions_to_keep <- c("A", "B")

  result <- filterMixtureComponents(conditions_to_keep, U, w, wCutoff = 1e-04)
  expect_true(!"mat2" %in% names(result$U))
})


test_that("filterMixtureComponents errors on missing condition", {
  U <- list(
    mat1 = matrix(c(1, 0.5, 0.5, 1), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  )
  w <- c(mat1 = 1.0)
  expect_error(filterMixtureComponents(c("A", "C"), U, w), "not found in matrix")
})


test_that("filterMixtureComponents subsets conditions", {
  U <- list(
    mat1 = matrix(c(1, 0.5, 0.2, 0.5, 1, 0.3, 0.2, 0.3, 1), 3, 3,
                  dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  )
  w <- c(mat1 = 1.0)

  result <- filterMixtureComponents(c("A", "B"), U, w)
  expect_equal(nrow(result$U[[1]]), 2)
  expect_equal(ncol(result$U[[1]]), 2)
})

# ===========================================================================
# Tests from test_misc_round3.R (mrmashWrapper coverage boost)
# ===========================================================================

# =========================================================================
# mrmashWrapper.R: compute_w0 (lines 284-298)
# =========================================================================



# ===========================================================================
# .mashSumStatsToMatrices — inputScale resolution
# ===========================================================================

# Helper: tiny QtlSumStats with 2 contexts on one trait. The mcols
# layout is controlled so each test can vary which columns are present.
.mssm_makeQtlSumStats <- function(mcolsBuilder, contexts = c("brain", "liver"),
                                   nSnp = 5L) {
  set.seed(13L)
  gh <- new("GenotypeHandle",
    path = "/tmp/sketch.gds", format = "gds",
    snpInfo = data.frame(
      SNP = paste0("v", seq_len(nSnp)), CHR = "1",
      BP = seq(100L, by = 100L, length.out = nSnp),
      A1 = "A", A2 = "G", stringsAsFactors = FALSE),
    nSamples = 50L, sampleIds = paste0("s", seq_len(50L)), pgenPtr = NULL)
  ranges <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = seq(100L, by = 100L, length.out = nSnp),
                              width = 1L))
  entries <- lapply(seq_along(contexts), function(i) {
    gr <- ranges
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(mcolsBuilder(i, nSnp))
    gr
  })
  QtlSumStats(
    study   = rep("s1", length(contexts)),
    context = contexts,
    trait   = rep("g1", length(contexts)),
    entry   = entries,
    genome  = "hg19",
    ldSketch = gh,
    qcInfo  = list(prebuilt = "synthetic"))
}

test_that(".mashSumStatsToMatrices: auto picks BETA+SE when present", {
  ss <- .mssm_makeQtlSumStats(function(i, n)
    list(SNP = paste0("v", seq_len(n)), A1 = "A", A2 = "G",
         Z    = rnorm(n),
         BETA = rnorm(n, sd = 0.1),
         SE   = abs(rnorm(n, sd = 0.05)) + 0.01))
  out <- pecotmr:::.mashSumStatsToMatrices(ss, "strong", inputScale = "auto")
  expect_equal(ncol(out$b), 2L)        # 2 contexts
  # On BETA scale, Shat values should be the small SEs we generated.
  expect_true(all(out$s[out$s < 1000] < 1))
})

test_that(".mashSumStatsToMatrices: auto falls back to Z when no BETA/SE", {
  ss <- .mssm_makeQtlSumStats(function(i, n)
    list(SNP = paste0("v", seq_len(n)), A1 = "A", A2 = "G",
         Z = rnorm(n)))
  out <- pecotmr:::.mashSumStatsToMatrices(ss, "strong", inputScale = "auto")
  # Shat should be 1 on the Z scale.
  expect_true(all(out$s == 1 | out$s == 1000))
})

test_that(".mashSumStatsToMatrices: inputScale='beta' errors when BETA missing", {
  ss <- .mssm_makeQtlSumStats(function(i, n)
    list(SNP = paste0("v", seq_len(n)), A1 = "A", A2 = "G",
         Z = rnorm(n)))
  expect_error(
    pecotmr:::.mashSumStatsToMatrices(ss, "strong", inputScale = "beta"),
    "BETA and SE")
})

test_that(".mashSumStatsToMatrices: inputScale='z' forces Z+1 even when BETA present", {
  ss <- .mssm_makeQtlSumStats(function(i, n)
    list(SNP = paste0("v", seq_len(n)), A1 = "A", A2 = "G",
         Z    = rnorm(n),
         BETA = rnorm(n, sd = 0.1),
         SE   = abs(rnorm(n, sd = 0.05)) + 0.01))
  out <- pecotmr:::.mashSumStatsToMatrices(ss, "strong", inputScale = "z")
  # Forced Z scale: Shat must be 1 everywhere except the NA fill (1000).
  expect_true(all(out$s == 1 | out$s == 1000))
})

test_that(".mashSumStatsToMatrices: errors when no usable scale", {
  ss <- .mssm_makeQtlSumStats(function(i, n)
    list(SNP = paste0("v", seq_len(n)), A1 = "A", A2 = "G",
         N = rep(1000L, n)))   # only N — no Z, no BETA/SE
  expect_error(
    pecotmr:::.mashSumStatsToMatrices(ss, "strong", inputScale = "auto"),
    "no usable scale")
})

test_that(".mashSumStatsToMatrices: inputScale='z' errors when Z missing", {
  ss <- .mssm_makeQtlSumStats(function(i, n)
    list(SNP = paste0("v", seq_len(n)), A1 = "A", A2 = "G",
         BETA = rnorm(n, sd = 0.1),
         SE   = abs(rnorm(n, sd = 0.05)) + 0.01))   # BETA/SE but no Z
  expect_error(
    pecotmr:::.mashSumStatsToMatrices(ss, "strong", inputScale = "z"),
    "carry a Z mcol")
})

# ===========================================================================
# .mashObjectMatrices / .mashObjectPartitions
# ===========================================================================

test_that(".mashObjectMatrices errors when marginal effects lack the required columns", {
  res <- QtlFineMappingResult(study = "s", context = "brain", trait = "t",
                              method = "susie", entry = list(.sc_makeFineMappingEntry(3)))
  # Force a marginal-effects table with no `context` column (as an mv/f-SuSiE
  # result trimmed of marginal sumstats would yield): mash cannot pivot it.
  testthat::local_mocked_bindings(
    getMarginalEffects = function(x, ...) data.frame(variant_id = "v", beta = 1, se = 1),
    .package = "pecotmr")
  expect_error(
    pecotmr:::.mashObjectMatrices(res, inputScale = "auto", coverage = 0.95),
    ">= 2 contexts")
})

test_that(".mashObjectMatrices warns and pins the first method on a multi-method result", {
  res <- QtlFineMappingResult(
    study = c("s", "s"), context = c("brain", "liver"), trait = c("t", "t"),
    method = c("susie", "mvsusie"),
    entry = list(.sc_makeFineMappingEntry(3), .sc_makeFineMappingEntry(3)))
  expect_warning(
    pecotmr:::.mashObjectMatrices(res, inputScale = "auto", coverage = 0.95),
    "multiple methods")
})

test_that(".mashObjectPartitions errors when < 2 conditions remain after excludeCondition", {
  ss <- .mssm_makeQtlSumStats(function(i, n)
    list(SNP = paste0("v", seq_len(n)), A1 = "A", A2 = "G",
         BETA = rnorm(n, sd = 0.1), SE = abs(rnorm(n, sd = 0.05)) + 0.01))
  expect_error(
    pecotmr:::.mashObjectPartitions(ss, nRandom = 3, nNull = 3,
      excludeCondition = "liver", coverage = 0.95, inputScale = "auto", seed = 1),
    "fewer than 2 conditions")
})

test_that(".mashObjectPartitions warns when no variants match the independent-variant list", {
  ss <- .mssm_makeQtlSumStats(function(i, n)
    list(SNP = paste0("v", seq_len(n)), A1 = "A", A2 = "G",
         BETA = rnorm(n, sd = 0.1), SE = abs(rnorm(n, sd = 0.05)) + 0.01))
  # No independent variant matches the panel, so the random/null background is
  # empty for this object; the warning fires before the (empty) sampling.
  expect_warning(
    tryCatch(
      pecotmr:::.mashObjectPartitions(ss, nRandom = 2, nNull = 2,
        excludeCondition = character(0), coverage = 0.95, inputScale = "auto",
        seed = 1, independentVariants = c("chrX:999999:N:N")),
      error = function(e) invisible(NULL)),
    "no variants matched")
})

# ===========================================================================
# .mashSumStatsToMatrices — matrix assembly behaviour (NA fill,
#   rowname disambiguation, multi-context shape) on the bundled fixture
#   and on hand-built partial-coverage data.
# ===========================================================================

test_that(".mashSumStatsToMatrices on bundled multicontext fixture: shape and rowname format", {
  data(qtl_sumstats_multicontext_example)
  out <- pecotmr:::.mashSumStatsToMatrices(
    qtl_sumstats_multicontext_example, "strong", inputScale = "auto")
  expect_equal(ncol(out$b), 3L)
  expect_equal(colnames(out$b), c("brain", "blood", "muscle"))
  # One (study, trait) block, 200 variants -> 200 rows
  expect_equal(nrow(out$b), 200L)
  # Rownames are disambiguated by (study::trait::variant)
  expect_true(all(grepl("^study1::ENSG_example::", rownames(out$b))))
  # On the BETA scale, sbhat values are small; no NA fill needed since every
  # context has every variant
  expect_true(all(out$s < 1))
})

test_that(".mashSumStatsToMatrices fills missing variants with bhat=0 / sbhat=1000", {
  set.seed(42L)
  gh <- new("GenotypeHandle",
    path = "/tmp/sketch.gds", format = "gds",
    snpInfo = data.frame(SNP = paste0("v", 1:5), CHR = "1",
                         BP = seq(100L, by = 100L, length.out = 5L),
                         A1 = "A", A2 = "G", stringsAsFactors = FALSE),
    nSamples = 50L, sampleIds = paste0("s", seq_len(50L)), pgenPtr = NULL)
  mkGr <- function(snpIds) {
    gr <- GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges = IRanges::IRanges(
        start = seq(100L, by = 100L, length.out = length(snpIds)), width = 1L))
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
      SNP = snpIds, A1 = "A", A2 = "G",
      Z = rnorm(length(snpIds)),
      BETA = rnorm(length(snpIds), sd = 0.1),
      SE   = rep(0.05, length(snpIds)))
    gr
  }
  # ctx1 has all 5 variants; ctx2 has only the first 3
  ss <- QtlSumStats(
    study = c("s1", "s1"), context = c("ctx1", "ctx2"),
    trait = c("g1", "g1"),
    entry = list(mkGr(paste0("v", 1:5)), mkGr(paste0("v", 1:3))),
    genome = "hg19", ldSketch = gh,
    qcInfo = list(prebuilt = "synthetic"))
  out <- pecotmr:::.mashSumStatsToMatrices(ss, "strong", inputScale = "auto")
  expect_equal(dim(out$b), c(5L, 2L))
  expect_setequal(colnames(out$b), c("ctx1", "ctx2"))
  # In ctx2, the last 2 variants are missing -> bhat NA -> 0, shat NA -> 1000
  expect_equal(unname(out$b[4:5, "ctx2"]), c(0, 0))
  expect_equal(unname(out$s[4:5, "ctx2"]), c(1000, 1000))
  # ctx1 has them present
  expect_true(all(abs(out$b[, "ctx1"]) < 1))
  expect_true(all(out$s[, "ctx1"] < 1))
})

test_that(".mashSumStatsToMatrices disambiguates rownames across (study, trait) blocks", {
  set.seed(7L)
  gh <- new("GenotypeHandle",
    path = "/tmp/sketch.gds", format = "gds",
    snpInfo = data.frame(SNP = paste0("v", 1:3), CHR = "1",
                         BP = c(100L, 200L, 300L),
                         A1 = "A", A2 = "G", stringsAsFactors = FALSE),
    nSamples = 50L, sampleIds = paste0("s", seq_len(50L)), pgenPtr = NULL)
  mkGr <- function(snpIds) {
    gr <- GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges = IRanges::IRanges(
        start = seq(100L, by = 100L, length.out = length(snpIds)), width = 1L))
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
      SNP = snpIds, A1 = "A", A2 = "G",
      Z = rnorm(length(snpIds)),
      BETA = rnorm(length(snpIds), sd = 0.1),
      SE   = rep(0.05, length(snpIds)))
    gr
  }
  # Two (study, trait) blocks but they share SNP IDs v1, v2, v3 — without
  # the prefix the rbind would silently merge them.
  ss <- QtlSumStats(
    study   = c("s1", "s1"),
    context = c("ctx1", "ctx1"),
    trait   = c("g1", "g2"),
    entry   = list(mkGr(paste0("v", 1:3)), mkGr(paste0("v", 1:3))),
    genome  = "hg19", ldSketch = gh,
    qcInfo  = list(prebuilt = "synthetic"))
  out <- pecotmr:::.mashSumStatsToMatrices(ss, "strong", inputScale = "auto")
  # 3 variants per block × 2 blocks = 6 rows
  expect_equal(nrow(out$b), 6L)
  expect_setequal(rownames(out$b),
                   c(paste0("s1::g1::v", 1:3),
                     paste0("s1::g2::v", 1:3)))
})

test_that(".mashSumStatsToMatrices errors when entry lacks SNP mcol", {
  set.seed(8L)
  gh <- new("GenotypeHandle",
    path = "/tmp/sketch.gds", format = "gds",
    snpInfo = data.frame(SNP = paste0("v", 1:3), CHR = "1",
                         BP = c(100L, 200L, 300L),
                         A1 = "A", A2 = "G", stringsAsFactors = FALSE),
    nSamples = 50L, sampleIds = paste0("s", seq_len(50L)), pgenPtr = NULL)
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(100L, 200L, 300L), width = 1L))
  # NO SNP mcol — should trigger the variant-alignment error
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    A1 = "A", A2 = "G",
    Z    = rnorm(3),
    BETA = rnorm(3, sd = 0.1), SE = rep(0.05, 3))
  ss <- QtlSumStats(
    study = "s1", context = "ctx1", trait = "g1",
    entry = list(gr),
    genome = "hg19", ldSketch = gh,
    qcInfo = list(prebuilt = "synthetic"))
  expect_error(
    pecotmr:::.mashSumStatsToMatrices(ss, "strong", inputScale = "auto"),
    "SNP")
})

# ===========================================================================
# .mashSumStatsToMatrices — GwasSumStats path + input-validation errors
# ===========================================================================

test_that(".mashSumStatsToMatrices on GwasSumStats: studies become columns", {
  set.seed(11L)
  gh <- new("GenotypeHandle",
    path = "/tmp/sketch.gds", format = "gds",
    snpInfo = data.frame(SNP = paste0("v", 1:3), CHR = "1",
                         BP = c(100L, 200L, 300L),
                         A1 = "A", A2 = "G", stringsAsFactors = FALSE),
    nSamples = 50L, sampleIds = paste0("s", seq_len(50L)), pgenPtr = NULL)
  mkGr <- function(snpIds) {
    gr <- GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges = IRanges::IRanges(
        start = seq(100L, by = 100L, length.out = length(snpIds)), width = 1L))
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
      SNP = snpIds, A1 = "A", A2 = "G",
      Z = rnorm(length(snpIds)),
      BETA = rnorm(length(snpIds), sd = 0.1),
      SE   = rep(0.05, length(snpIds)))
    gr
  }
  # Each study is its own (study) block; columns of the mash matrix are the
  # studies, so the result is block-diagonal with NA-fill off the diagonal.
  ss <- GwasSumStats(
    study = c("studyA", "studyB"),
    entry = list(mkGr(paste0("v", 1:3)), mkGr(paste0("v", 1:3))),
    genome = "hg19", ldSketch = gh,
    qcInfo = list(prebuilt = "synthetic"))
  out <- pecotmr:::.mashSumStatsToMatrices(ss, "strong", inputScale = "auto")
  expect_equal(ncol(out$b), 2L)
  expect_equal(colnames(out$b), c("studyA", "studyB"))
  # 2 studies x 3 variants = 6 rows; rownames prefixed by the study block key.
  expect_equal(nrow(out$b), 6L)
  expect_setequal(rownames(out$b),
                  c(paste0("studyA::v", 1:3), paste0("studyB::v", 1:3)))
  # studyA's rows are absent from studyB's column -> bhat 0 / shat 1000 fill.
  studyArows <- grep("^studyA::", rownames(out$b))
  expect_equal(unname(out$b[studyArows, "studyB"]), rep(0, 3))
  expect_equal(unname(out$s[studyArows, "studyB"]), rep(1000, 3))
  # On the BETA scale, the present cells carry the small generated SEs.
  expect_true(all(out$s[studyArows, "studyA"] < 1))
})

test_that(".mashSumStatsToMatrices errors on a non-SumStats input", {
  expect_error(
    pecotmr:::.mashSumStatsToMatrices(list(a = 1), "strong"),
    "must be a QtlSumStats or GwasSumStats")
})

test_that(".mashSumStatsToMatrices errors when SumStats has empty QC info", {
  gh <- new("GenotypeHandle",
    path = "/tmp/sketch.gds", format = "gds",
    snpInfo = data.frame(SNP = paste0("v", 1:3), CHR = "1",
                         BP = c(100L, 200L, 300L),
                         A1 = "A", A2 = "G", stringsAsFactors = FALSE),
    nSamples = 50L, sampleIds = paste0("s", seq_len(50L)), pgenPtr = NULL)
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(100L, 200L, 300L), width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = paste0("v", 1:3), A1 = "A", A2 = "G",
    Z = rnorm(3), BETA = rnorm(3, sd = 0.1), SE = rep(0.05, 3))
  ss <- QtlSumStats(study = "s1", context = "c1", trait = "g1",
                    entry = list(gr), genome = "hg19", ldSketch = gh,
                    qcInfo = list())  # empty QC info
  expect_error(
    pecotmr:::.mashSumStatsToMatrices(ss, "strong"),
    "no QC info")
})

test_that(".mashSumStatsToMatrices errors when SumStats has zero entries", {
  gh <- new("GenotypeHandle",
    path = "/tmp/sketch.gds", format = "gds",
    snpInfo = data.frame(SNP = paste0("v", 1:3), CHR = "1",
                         BP = c(100L, 200L, 300L),
                         A1 = "A", A2 = "G", stringsAsFactors = FALSE),
    nSamples = 50L, sampleIds = paste0("s", seq_len(50L)), pgenPtr = NULL)
  ss <- QtlSumStats(study = character(0), context = character(0),
                    trait = character(0), entry = list(), genome = "hg19",
                    ldSketch = gh, varY = numeric(0),
                    qcInfo = list(prebuilt = "synthetic"))
  expect_error(
    pecotmr:::.mashSumStatsToMatrices(ss, "strong"),
    "no entries")
})

# ===========================================================================
# mergeMashData — oneData-empty branch (returns resData unchanged)
# ===========================================================================

test_that("mergeMashData returns resData when oneData is NULL", {
  d1 <- list(random = data.frame(a = 1:3, b = 4:6))
  expect_equal(mergeMashData(d1, NULL), d1)
})

test_that("mergeMashData returns resData when oneData is an empty list", {
  d1 <- list(random = data.frame(a = 1:3))
  expect_equal(mergeMashData(d1, list()), d1)
})

# ===========================================================================
# qtlSumStatsFromZMatrix
# ===========================================================================

.qszm_gh <- function() {
  new("GenotypeHandle",
    path = "/tmp/sketch.gds", format = "gds",
    snpInfo = data.frame(SNP = "v1", CHR = "1", BP = 100L,
                         A1 = "A", A2 = "G", stringsAsFactors = FALSE),
    nSamples = 10L, sampleIds = paste0("s", seq_len(10L)), pgenPtr = NULL)
}

test_that("qtlSumStatsFromZMatrix: one row per context, Z preserved verbatim", {
  z <- matrix(c(1.1, -2.2, 0.3, 0.4, -0.5, 0.6), nrow = 3,
              dimnames = list(c("chr1:100:A:G", "chr1:200:C:T", "chr2:50:A:T"),
                              c("brain", "liver")))
  qss <- qtlSumStatsFromZMatrix(z, study = "s1", ldSketch = .qszm_gh())
  expect_s4_class(qss, "QtlSumStats")
  expect_equal(nrow(qss), 2L)
  expect_equal(as.character(qss$context), c("brain", "liver"))
  expect_equal(unique(as.character(qss$study)), "s1")
  expect_equal(unique(as.character(qss$trait)), "mash")
  # Z column for each context matches the input matrix column (values only;
  # the mcols column is unnamed whereas z[, j] carries the row ids as names).
  expect_equal(S4Vectors::mcols(qss$entry[[1]])$Z, unname(z[, 1]))
  expect_equal(S4Vectors::mcols(qss$entry[[2]])$Z, unname(z[, 2]))
})

test_that("qtlSumStatsFromZMatrix: decodes chrom/pos from ids, synthesises when they don't parse", {
  z <- matrix(rnorm(2), ncol = 1,
              dimnames = list(c("chr1:250:A:G", "not_a_variant"), "ctx"))
  qss <- qtlSumStatsFromZMatrix(z, study = "s1", ldSketch = .qszm_gh())
  e <- qss$entry[[1]]
  expect_equal(GenomicRanges::start(e)[1], 250L)          # decoded
  expect_true(GenomicRanges::start(e)[2] >= 1L)           # synthetic fallback
  # un-parseable chrom falls back to chr1 (never NA)
  expect_false(any(is.na(as.character(GenomicRanges::seqnames(e)))))
  # SNP ids are carried through unchanged
  expect_equal(S4Vectors::mcols(e)$SNP, rownames(z))
})

test_that("qtlSumStatsFromZMatrix: NULL rownames get synthetic variant ids", {
  z <- matrix(rnorm(4), nrow = 2, dimnames = list(NULL, c("a", "b")))
  qss <- qtlSumStatsFromZMatrix(z, study = "s1", ldSketch = .qszm_gh())
  expect_equal(S4Vectors::mcols(qss$entry[[1]])$SNP, c("var1", "var2"))
})

test_that("qtlSumStatsFromZMatrix: placeholders and pass-through qcInfo are set", {
  z <- matrix(rnorm(6), nrow = 3, dimnames = list(NULL, c("x", "y")))
  qss <- qtlSumStatsFromZMatrix(z, study = "s1", ldSketch = .qszm_gh(),
                                n = 500L, a1 = "T", a2 = "C", role = "strong")
  mc <- S4Vectors::mcols(qss$entry[[1]])
  expect_equal(unique(mc$A1), "T")
  expect_equal(unique(mc$A2), "C")
  expect_equal(unique(mc$N), 500L)
  expect_equal(qss@qcInfo$role, "strong")
  expect_equal(length(qss@qcInfo$entryAudit), 2L)   # one slot per context
})

test_that("qtlSumStatsFromZMatrix: columns can map to traits or context x trait pairs", {
  z <- matrix(rnorm(6), nrow = 3, dimnames = list(NULL, c("geneA", "geneB")))
  # columns as traits: constant context, one trait per column
  qss <- qtlSumStatsFromZMatrix(z, study = "s1", ldSketch = .qszm_gh(),
                                context = "brain", trait = colnames(z))
  expect_equal(as.character(qss$context), c("brain", "brain"))
  expect_equal(as.character(qss$trait), c("geneA", "geneB"))
  # columns as (context, trait) pairs
  qss2 <- qtlSumStatsFromZMatrix(z, study = "s1", ldSketch = .qszm_gh(),
                                 context = c("brain", "liver"),
                                 trait   = c("geneA", "geneA"))
  expect_equal(as.character(qss2$context), c("brain", "liver"))
  expect_equal(as.character(qss2$trait), c("geneA", "geneA"))
})

test_that("qtlSumStatsFromZMatrix: a condition label of the wrong length errors", {
  z <- matrix(rnorm(6), nrow = 3, dimnames = list(NULL, c("a", "b")))
  expect_error(
    qtlSumStatsFromZMatrix(z, study = "s1", ldSketch = .qszm_gh(),
                           trait = c("t1", "t2", "t3")),
    "must be length 1 or ncol")
})

test_that("qtlSumStatsFromZMatrix: rejects non-matrix input and unlabelled conditions", {
  expect_error(qtlSumStatsFromZMatrix(1:5, study = "s1", ldSketch = .qszm_gh()),
               "variants x conditions matrix")
  # no colnames -> the default context = colnames(z) is NULL
  z <- matrix(rnorm(4), nrow = 2)
  expect_error(qtlSumStatsFromZMatrix(z, study = "s1", ldSketch = .qszm_gh()),
               "column names")
})

# ===========================================================================
# qtlSumStatsFromBetaMatrix (beta-scale sibling of qtlSumStatsFromZMatrix)
# ===========================================================================

.qszmBeta <- function() {
  rn <- c("chr1:100:A:G", "chr1:200:C:T", "chr2:50:A:T")
  cn <- c("brain", "liver")
  list(bhat = matrix(c(0.5, -0.3, 0.2, 0.1, 0.4, -0.6), nrow = 3,
                     dimnames = list(rn, cn)),
       shat = matrix(c(0.1, 0.2, 0.15, 0.12, 0.09, 0.2), nrow = 3,
                     dimnames = list(rn, cn)))
}

test_that("qtlSumStatsFromBetaMatrix: one entry per context, BETA/SE/Z mcols set", {
  d <- .qszmBeta()
  qss <- qtlSumStatsFromBetaMatrix(d$bhat, d$shat, study = "s1",
                                   ldSketch = .qszm_gh())
  expect_s4_class(qss, "QtlSumStats")
  expect_equal(nrow(qss), 2L)
  expect_equal(as.character(qss$context), c("brain", "liver"))
  mc <- S4Vectors::mcols(qss$entry[[1]])
  expect_true(all(c("BETA", "SE", "Z") %in% colnames(mc)))
  expect_equal(mc$BETA, unname(d$bhat[, 1]))
  expect_equal(mc$SE, unname(d$shat[, 1]))
  expect_equal(mc$Z, unname(d$bhat[, 1] / d$shat[, 1]))
})

test_that("qtlSumStatsFromBetaMatrix: feeds .mashSumStatsToMatrices on both scales", {
  d <- .qszmBeta()
  qss <- qtlSumStatsFromBetaMatrix(d$bhat, d$shat, study = "s1",
                                   ldSketch = .qszm_gh())
  mb <- .mashSumStatsToMatrices(qss, "strong", inputScale = "beta")
  expect_equal(unname(mb$b), unname(d$bhat))
  expect_equal(unname(mb$s), unname(d$shat))
  mz <- .mashSumStatsToMatrices(qss, "strong", inputScale = "z")
  expect_equal(unname(mz$b), unname(d$bhat / d$shat))
})

test_that("qtlSumStatsFromBetaMatrix: validates matrices and matching dimensions", {
  d <- .qszmBeta()
  expect_error(qtlSumStatsFromBetaMatrix(1:5, d$shat, "s1", .qszm_gh()),
               "`bhat` must be a numeric")
  expect_error(qtlSumStatsFromBetaMatrix(d$bhat, 1:5, "s1", .qszm_gh()),
               "`shat` must be a numeric")
  expect_error(qtlSumStatsFromBetaMatrix(d$bhat, d$shat[1:2, ], "s1", .qszm_gh()),
               "identical dimensions")
})

test_that("qtlSumStatsFromBetaMatrix: NULL rownames -> synthetic ids; placeholders set", {
  bhat <- matrix(rnorm(4), nrow = 2, dimnames = list(NULL, c("x", "y")))
  shat <- matrix(abs(rnorm(4)) + 0.1, nrow = 2, dimnames = list(NULL, c("x", "y")))
  qss <- qtlSumStatsFromBetaMatrix(bhat, shat, study = "s1", ldSketch = .qszm_gh(),
                                   n = 500L, a1 = "T", a2 = "C", role = "strong")
  mc <- S4Vectors::mcols(qss$entry[[1]])
  expect_equal(mc$SNP, c("var1", "var2"))
  expect_equal(unique(mc$A1), "T")
  expect_equal(unique(mc$N), 500L)
  expect_equal(qss@qcInfo$role, "strong")
})

# ===========================================================================
# mashInput  (unified strong/random/null assembly from S4 objects)
# ===========================================================================

# A multi-context QtlFineMappingResult fixture: two contexts sharing the same
# 6 variants, each with one credible set whose lead (max PIP) is a distinct
# variant, plus low-signal background for random/null sampling.
.mi_makeFmr <- function() {
  mkTL <- function(zvec, csIdx) {
    n <- length(zvec)
    data.frame(
      variant_id = paste0("chr1:", 100 * seq_len(n), ":A:G"),
      chrom = "1", pos = as.integer(100 * seq_len(n)), A1 = "G", A2 = "A",
      N = 1000, MAF = 0.2,
      marginal_beta = zvec * 0.05, marginal_se = 0.05, marginal_z = zvec,
      marginal_p = 2 * pnorm(-abs(zvec)),
      pip = { p <- rep(0.05, n); p[csIdx] <- seq(0.9, by = -0.2,
              length.out = length(csIdx)); p },
      posterior_mean = zvec * 0.05, posterior_sd = 0.02,
      cs_95 = { cc <- rep("susie_0", n); cc[csIdx] <- "susie_1"; cc },
      stringsAsFactors = FALSE)
  }
  vids <- paste0("chr1:", 100 * seq_len(6), ":A:G")
  e1 <- FineMappingEntry(variantIds = vids, susieFit = list(x = 1),
                         topLoci = mkTL(c(0.5, -1, 6.0, 0.2, 1.1, -0.3), c(3, 2)))
  e2 <- FineMappingEntry(variantIds = vids, susieFit = list(x = 1),
                         topLoci = mkTL(c(-0.4, 0.7, 0.1, 5.0, -1.2, 0.6), c(4, 5)))
  QtlFineMappingResult(study = c("s1", "s1"), context = c("brain", "blood"),
                       trait = c("t1", "t1"), method = c("susie", "susie"),
                       entry = list(e1, e2))
}

test_that("mashInput: QtlSumStats path returns the flat b/s/z + XtX contract", {
  data(qtl_sumstats_multicontext_example)
  ss <- qtl_sumstats_multicontext_example
  out <- mashInput(list(geneA = ss), nRandom = 5, nNull = 5, seed = 1)
  for (k in c("strong.b", "strong.s", "strong.z", "random.b", "random.s",
              "random.z", "null.b", "null.s", "null.z", "XtX")) {
    expect_true(k %in% names(out), info = k)
  }
  nCond <- ncol(out$strong.z)
  expect_equal(nrow(out$random.z), 5L)
  expect_equal(nrow(out$null.z), 5L)
  # XtX is conditions x conditions and symmetric
  expect_equal(dim(out$XtX), c(nCond, nCond))
  expect_equal(unname(out$XtX), unname(t(out$XtX)))
})

test_that("mashInput: FineMappingResult strong = CS lead (max PIP) per condition", {
  fmr <- .mi_makeFmr()
  out <- mashInput(list(geneA = fmr), nRandom = 4, nNull = 4, coverage = 0.95,
                   sigPCutoff = 0.5, seed = 7)
  expect_equal(colnames(out$strong.z), c("brain", "blood"))
  # brain CS lead = variant 3 (chr1:300), blood CS lead = variant 4 (chr1:400)
  expect_setequal(sub("_geneA$", "", rownames(out$strong.z)),
                  c("chr1:300:A:G", "chr1:400:A:G"))
  expect_equal(nrow(out$random.z), 4L)
})

test_that("mashInput: a FineMappingResult with no credible set yields no strong", {
  vids <- paste0("chr1:", 100 * seq_len(6), ":A:G")
  noCsTL <- function(zvec) data.frame(
    variant_id = vids, chrom = "1", pos = as.integer(100 * seq_len(6)),
    A1 = "G", A2 = "A", N = 1000, MAF = 0.2,
    marginal_beta = zvec * 0.05, marginal_se = 0.05, marginal_z = zvec,
    marginal_p = 2 * pnorm(-abs(zvec)),
    pip = rep(0.05, 6), posterior_mean = zvec * 0.05, posterior_sd = 0.02,
    cs_95 = rep("susie_0", 6), stringsAsFactors = FALSE)
  e1 <- FineMappingEntry(variantIds = vids, susieFit = list(x = 1),
                         topLoci = noCsTL(rnorm(6)))
  e2 <- FineMappingEntry(variantIds = vids, susieFit = list(x = 1),
                         topLoci = noCsTL(rnorm(6)))
  fmr <- QtlFineMappingResult(study = c("s1", "s1"),
                              context = c("brain", "blood"),
                              trait = c("t1", "t1"),
                              method = c("susie", "susie"), entry = list(e1, e2))
  out <- mashInput(list(g = fmr), nRandom = 3, nNull = 3, seed = 2)
  expect_null(out$strong.z)
  expect_equal(nrow(out$random.z), 3L)
})

test_that("mashInput: multiple regions accumulate rows and disambiguate names", {
  fmr <- .mi_makeFmr()
  one <- mashInput(list(a = fmr), nRandom = 4, nNull = 4, sigPCutoff = 0.5, seed = 7)
  two <- mashInput(list(a = fmr, b = fmr), nRandom = 4, nNull = 4,
                   sigPCutoff = 0.5, seed = 7)
  expect_equal(nrow(two$strong.z), 2L * nrow(one$strong.z))
  expect_equal(nrow(two$random.z), 2L * nrow(one$random.z))
  expect_false(any(duplicated(rownames(two$random.z))))
})

test_that("mashInput: zOnly = TRUE drops the .b/.s matrices", {
  fmr <- .mi_makeFmr()
  out <- mashInput(list(g = fmr), nRandom = 3, nNull = 3, zOnly = TRUE,
                   sigPCutoff = 0.5, seed = 7)
  expect_false(any(grepl("\\.(b|s)$", names(out))))
  expect_true(all(c("strong.z", "random.z", "null.z", "XtX") %in% names(out)))
})

test_that("mashInput: excludeCondition drops the condition column everywhere", {
  # Exclude one of three conditions (mash needs >= 2), leaving brain + blood.
  data(qtl_sumstats_multicontext_example)
  ss <- qtl_sumstats_multicontext_example
  out <- mashInput(list(g = ss), nRandom = 4, nNull = 4,
                   excludeCondition = "muscle", seed = 7)
  expect_false("muscle" %in% colnames(out$random.z))
  expect_setequal(colnames(out$random.z), c("brain", "blood"))
})

test_that("mashInput: rejects a non-SumStats/non-FineMapping element", {
  expect_error(mashInput(list(1L)),
               "must be a QtlSumStats, GwasSumStats, or FineMappingResult")
  expect_error(mashInput(list()), "non-empty list")
})

test_that("mashInput: independentVariants restricts random/null pool (strong untouched)", {
  fmr <- .mi_makeFmr()                       # 6 variants chr1:100..600, 2 contexts
  vids <- paste0("chr1:", 100 * seq_len(6), ":A:G")
  indep <- vids[1:3]                         # only the first three are "independent"
  out <- mashInput(list(g = fmr), nRandom = 3, nNull = 3,
                   independentVariants = indep, sigPCutoff = 0.5, seed = 7)
  # random / null variant ids (strip the "_g" region suffix) must lie in indep.
  rn_ids <- sub("_g$", "", rownames(out$random.z))
  expect_true(length(rn_ids) > 0L && all(rn_ids %in% indep))
  if (!is.null(out$null.z)) {
    expect_true(all(sub("_g$", "", rownames(out$null.z)) %in% indep))
  }
  # strong is NOT filtered -- its CS lead may lie outside the independent set.
  expect_true(nrow(out$strong.z) >= 1L)
})

test_that("mashInput: independentVariants matches across allele flip + chr prefix", {
  fmr <- .mi_makeFmr()
  # same positions as the first three variants, but ref/alt swapped and the
  # chr prefix dropped -- matchVariants must still match them (not a string cmp).
  indep_flipped <- c("1:100:G:A", "1:200:G:A", "1:300:G:A")
  out <- mashInput(list(g = fmr), nRandom = 3, nNull = 3,
                   independentVariants = indep_flipped, sigPCutoff = 0.5, seed = 7)
  rn_ids <- sub("_g$", "", rownames(out$random.z))
  expect_true(length(rn_ids) > 0L)
  expect_true(all(rn_ids %in% paste0("chr1:", c(100, 200, 300), ":A:G")))
})
