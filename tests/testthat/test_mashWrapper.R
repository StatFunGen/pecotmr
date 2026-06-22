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
# mergeSusieCs
# ===========================================================================

test_that("mergeSusieCs merges credible sets correctly", {
  # Test case 1: No overlapping credible sets
  susie_fit_1 <- list(
    list(
      condition_1 = list(
        top_loci = data.frame(
          variant_id = c("variant1", "variant2"),
          pip = c(0.8, 0.6),
          CS_95_susie = c(1, 1)
        )
      ),
      condition_2 = list(
        top_loci = data.frame(
          variant_id = c("variant3", "variant4"),
          pip = c(0.9, 0.7),
          CS_95_susie = c(1, 2)
        )
      )
    )
  )

  expected_output_1 <- data.frame(
    variant_id = c("variant1", "variant2", "variant3", "variant4"),
    credibleSetNames = c("cs_1_1", "cs_1_1", "cs_2_1", "cs_2_2"),
    maxPip = c(0.8, 0.6, 0.9, 0.7),
    medianPip = c(0.8, 0.6, 0.9, 0.7),
    stringsAsFactors = FALSE
  )

  expect_equal(mergeSusieCs(susie_fit_1), expected_output_1)

  # Test case 2: Overlapping credible sets
  susie_fit_2 <- list(
    list(
      condition_1 = list(
        top_loci = data.frame(
          variant_id = c("variant1", "variant2"),
          pip = c(0.8, 0.6),
          CS_95_susie = c(1, 1)
        )
      ),
      condition_2 = list(
        top_loci = data.frame(
          variant_id = c("variant2", "variant3"),
          pip = c(0.7, 0.9),
          CS_95_susie = c(2, 2)
        )
      )
    )
  )

  expected_output_2 <- data.frame(
    variant_id = c("variant1", "variant2", "variant3"),
    credibleSetNames = c("cs_1_1,cs_2_2", "cs_1_1,cs_2_2", "cs_1_1,cs_2_2"),
    maxPip = c(0.8, 0.7, 0.9),
    medianPip = c(0.8, 0.65, 0.9),
    stringsAsFactors = FALSE
  )

  expect_equal(mergeSusieCs(susie_fit_2), expected_output_2)

  # Test case 3: Empty input
  susie_fit_3 <- list(condition_1 = list(top_loci = data.frame(
    variant_id = character(),
    credibleSetNames = character(),
    maxPip = numeric(),
    medianPip = numeric(),
    stringsAsFactors = FALSE
  )))

  expected_output_3 <- NULL

  expect_equal(mergeSusieCs(susie_fit_3), expected_output_3)

  # Test case 4: Different coverage parameter
  susie_fit_5 <- list(
    list(
      condition_1 = list(
        top_loci = data.frame(
          variant_id = c("variant1", "variant2"),
          pip = c(0.8, 0.6),
          CS_90_susie = c(1, 1)
        )
      ),
      condition_2 = list(
        top_loci = data.frame(
          variant_id = c("variant3", "variant4"),
          pip = c(0.9, 0.7),
          CS_90_susie = c(2, 2)
        )
      )
    )
  )

  expected_output_5 <- data.frame(
    variant_id = c("variant1", "variant2", "variant3", "variant4"),
    credibleSetNames = c("cs_1_1", "cs_1_1", "cs_2_2", "cs_2_2"),
    maxPip = c(0.8, 0.6, 0.9, 0.7),
    medianPip = c(0.8, 0.6, 0.9, 0.7),
    stringsAsFactors = FALSE
  )

  expect_equal(mergeSusieCs(susie_fit_5, coverage = "CS_90_susie"), expected_output_5)

  # Test case 6: Multiple top_loci tables with mixed coverage indices
  susie_fit_6 <- list(
    list(
      condition_1 = list(
        top_loci = data.frame(
          variant_id = c("variant1", "variant2", "variant3"),
          pip = c(0.8, 0.6, 0.7),
          CS_95_susie = c(1, 1, 2)
        )
      ),
      condition_2 = list(
        top_loci = data.frame(
          variant_id = c("variant4", "variant5"),
          pip = c(0.9, 0.7),
          CS_95_susie = c(2, 3)
        )
      ),
      condition_3 = list(
        top_loci = data.frame(
          variant_id = c("variant6", "variant7", "variant8"),
          pip = c(0.85, 0.75, 0.8),
          CS_95_susie = c(1, 3, 2)
        )
      )
    )
  )

  expected_output_6 <- data.frame(
    variant_id = c("variant1", "variant2", "variant3", "variant4", "variant5", "variant6", "variant7", "variant8"),
    credibleSetNames = c("cs_1_1", "cs_1_1", "cs_1_2", "cs_2_2", "cs_2_3", "cs_3_1", "cs_3_3", "cs_3_2"),
    maxPip = c(0.8, 0.6, 0.7, 0.9, 0.7, 0.85, 0.75, 0.8),
    medianPip = c(0.8, 0.6, 0.7, 0.9, 0.7, 0.85, 0.75, 0.8),
    stringsAsFactors = FALSE
  )

  expect_equal(mergeSusieCs(susie_fit_6), expected_output_6)

  # Test case 7: Multiple top_loci tables with overlapping sets and mixed coverage indices
  susie_fit_7 <- list(
    list(
      condition_1 = list(
        top_loci = data.frame(
          variant_id = c("variant1", "variant2", "variant3"),
          pip = c(0.8, 0.6, 0.7),
          CS_95_susie = c(1, 1, 2)
        )
      ),
      condition_2 = list(
        top_loci = data.frame(
          variant_id = c("variant2", "variant3", "variant4"),
          pip = c(0.7, 0.9, 0.85),
          CS_95_susie = c(2, 2, 1)
        )
      ),
      condition_3 = list(
        top_loci = data.frame(
          variant_id = c("variant4", "variant5"),
          pip = c(0.75, 0.8),
          CS_95_susie = c(3, 2)
        )
      )
    )
  )

  expected_output_7 <- data.frame(
    variant_id = c("variant1", "variant2", "variant3", "variant4", "variant5"),
    credibleSetNames = c("cs_1_1,cs_1_2,cs_2_2", "cs_1_1,cs_1_2,cs_2_2", "cs_1_1,cs_1_2,cs_2_2","cs_2_1,cs_3_3", "cs_3_2"),
    maxPip = c(0.8, 0.7, 0.9, 0.85, 0.8),
    medianPip = c(0.8, 0.65, 0.8, 0.8, 0.8),
    stringsAsFactors = FALSE
  )

  expect_equal(mergeSusieCs(susie_fit_7), expected_output_7)

  # Test case 8: Multiple top_loci tables with different coverage indices and no overlapping sets
  susie_fit_8 <- list(
    list(
      condition_1 = list(
        top_loci = data.frame(
          variant_id = c("variant1", "variant2", "variant3"),
          pip = c(0.8, 0.6, 0.7),
          CS_95_susie = c(1, 2, 3)
        )
      ),
      condition_2 = list(
        top_loci = data.frame(
          variant_id = c("variant4", "variant5"),
          pip = c(0.9, 0.7),
          CS_95_susie = c(3, 1)
        )
      ),
      condition_3 = list(
        top_loci = data.frame(
          variant_id = c("variant6", "variant7", "variant8"),
          pip = c(0.85, 0.75, 0.8),
          CS_95_susie = c(2, 3, 1)
        )
      )
    )
  )

  expected_output_8 <- data.frame(
    variant_id = c("variant1", "variant2", "variant3", "variant4", "variant5", "variant6", "variant7", "variant8"),
    credibleSetNames = c("cs_1_1", "cs_1_2", "cs_1_3", "cs_2_3", "cs_2_1", "cs_3_2", "cs_3_3", "cs_3_1"),
    maxPip = c(0.8, 0.6, 0.7, 0.9, 0.7, 0.85, 0.75, 0.8),
    medianPip = c(0.8, 0.6, 0.7, 0.9, 0.7, 0.85, 0.75, 0.8),
    stringsAsFactors = FALSE
  )

  expect_equal(mergeSusieCs(susie_fit_8), expected_output_8)

  # Test case 9: Single top_loci table with mixed coverage indices
  susie_fit_9 <- list(
    list(
      condition_1 = list(
        top_loci = data.frame(
          variant_id = c("variant1", "variant2", "variant3", "variant4", "variant5"),
          pip = c(0.8, 0.6, 0.7, 0.9, 0.85),
          CS_95_susie = c(1, 1, 2, 3, 2)
        )
      )
    )
  )

  expected_output_9 <- data.frame(
    variant_id = c("variant1", "variant2", "variant3", "variant5", "variant4"),
    credibleSetNames = c("cs_1_1", "cs_1_1", "cs_1_2", "cs_1_2", "cs_1_3"),
    maxPip = c(0.8, 0.6, 0.7, 0.85, 0.9),
    medianPip = c(0.8, 0.6, 0.7, 0.85, 0.9),
    stringsAsFactors = FALSE
  )

  expect_equal(mergeSusieCs(susie_fit_9), expected_output_9)

  # Test case 10: Multiple top_loci tables with mixed coverage indices and overlapping sets
  susie_fit_10 <- list(
    list(
      condition_1 = list(
        top_loci = data.frame(
          variant_id = c("variant1", "variant2", "variant3"),
          pip = c(0.8, 0.6, 0.7),
          CS_95_susie = c(1, 2, 1)
        )
      ),
      condition_2 = list(
        top_loci = data.frame(
          variant_id = c("variant2", "variant4", "variant5"),
          pip = c(0.75, 0.9, 0.85),
          CS_95_susie = c(2, 1, 3)
        )
      ),
      condition_3 = list(
        top_loci = data.frame(
          variant_id = c("variant3", "variant5", "variant6"),
          pip = c(0.65, 0.8, 0.7),
          CS_95_susie = c(3, 2, 1)
        )
      )
    )
  )

  expected_output_10 <- data.frame(
    variant_id = c("variant1", "variant3", "variant2", "variant4", "variant5", "variant6"),
    credibleSetNames = c("cs_1_1,cs_3_3", "cs_1_1,cs_3_3", "cs_1_2,cs_2_2", "cs_2_1", "cs_2_3,cs_3_2", "cs_3_1"),
    maxPip = c(0.8, 0.7, 0.75, 0.9, 0.85, 0.7),
    medianPip = c(0.8, 0.675, 0.675, 0.9, 0.825, 0.7),
    stringsAsFactors = FALSE
  )

  expect_equal(mergeSusieCs(susie_fit_10), expected_output_10)
})

test_that("mergeSusieCs handles single condition with single CS", {
  susie_fit <- list(list(
    cond1 = list(
      top_loci = data.frame(
        variant_id = c("1:100:A:G", "1:200:C:T"),
        pip = c(0.9, 0.1),
        CS_95_susie = c(1, 1),
        stringsAsFactors = FALSE
      )
    )
  ))

  result <- pecotmr:::mergeSusieCs(susie_fit)
  expect_s3_class(result, "data.frame")
  expect_true("variant_id" %in% colnames(result))
  expect_true("maxPip" %in% colnames(result))
})

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

test_that("mashRandNullSample excludeCondition with column names errors due to numeric indexing", {
  dat <- list(
    z = data.frame(
      cond1 = c(0.1, 0.2, 0.3, 0.4, 0.5),
      cond2 = c(0.5, 0.4, 0.3, 0.2, 0.1),
      cond3 = c(0.3, 0.3, 0.3, 0.3, 0.3)
    )
  )
  expect_error(
    mashRandNullSample(dat, nRandom = 3, nNull = 2,
                           excludeCondition = "cond3", seed = 42),
    "invalid argument to unary operator"
  )
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
