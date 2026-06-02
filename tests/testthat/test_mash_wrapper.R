context("mash_wrapper")

# Build a minimal FineMappingResult for unit-testing find_nested / extract_flatten_sumstats_from_nested
.test_fm_result <- function(variant_names) {
    FineMappingResult(
        variant_names = variant_names,
        trimmed_fit = list(pip = rep(0.5, length(variant_names))),
        top_loci = data.frame(variant_id = character(0), method = character(0),
                              stringsAsFactors = FALSE),
        method = "susie"
    )
}

# ===========================================================================
# merge_susie_cs
# ===========================================================================

test_that("merge_susie_cs merges credible sets correctly", {
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
    credible_set_names = c("cs_1_1", "cs_1_1", "cs_2_1", "cs_2_2"),
    max_pip = c(0.8, 0.6, 0.9, 0.7),
    median_pip = c(0.8, 0.6, 0.9, 0.7),
    stringsAsFactors = FALSE
  )

  expect_equal(merge_susie_cs(susie_fit_1), expected_output_1)

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
    credible_set_names = c("cs_1_1,cs_2_2", "cs_1_1,cs_2_2", "cs_1_1,cs_2_2"),
    max_pip = c(0.8, 0.7, 0.9),
    median_pip = c(0.8, 0.65, 0.9),
    stringsAsFactors = FALSE
  )

  expect_equal(merge_susie_cs(susie_fit_2), expected_output_2)

  # Test case 3: Empty input
  susie_fit_3 <- list(condition_1 = list(top_loci = data.frame(
    variant_id = character(),
    credible_set_names = character(),
    max_pip = numeric(),
    median_pip = numeric(),
    stringsAsFactors = FALSE
  )))

  expected_output_3 <- NULL

  expect_equal(merge_susie_cs(susie_fit_3), expected_output_3)

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
    credible_set_names = c("cs_1_1", "cs_1_1", "cs_2_2", "cs_2_2"),
    max_pip = c(0.8, 0.6, 0.9, 0.7),
    median_pip = c(0.8, 0.6, 0.9, 0.7),
    stringsAsFactors = FALSE
  )

  expect_equal(merge_susie_cs(susie_fit_5, coverage = "CS_90_susie"), expected_output_5)

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
    credible_set_names = c("cs_1_1", "cs_1_1", "cs_1_2", "cs_2_2", "cs_2_3", "cs_3_1", "cs_3_3", "cs_3_2"),
    max_pip = c(0.8, 0.6, 0.7, 0.9, 0.7, 0.85, 0.75, 0.8),
    median_pip = c(0.8, 0.6, 0.7, 0.9, 0.7, 0.85, 0.75, 0.8),
    stringsAsFactors = FALSE
  )

  expect_equal(merge_susie_cs(susie_fit_6), expected_output_6)

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
    credible_set_names = c("cs_1_1,cs_1_2,cs_2_2", "cs_1_1,cs_1_2,cs_2_2", "cs_1_1,cs_1_2,cs_2_2","cs_2_1,cs_3_3", "cs_3_2"),
    max_pip = c(0.8, 0.7, 0.9, 0.85, 0.8),
    median_pip = c(0.8, 0.65, 0.8, 0.8, 0.8),
    stringsAsFactors = FALSE
  )

  expect_equal(merge_susie_cs(susie_fit_7), expected_output_7)

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
    credible_set_names = c("cs_1_1", "cs_1_2", "cs_1_3", "cs_2_3", "cs_2_1", "cs_3_2", "cs_3_3", "cs_3_1"),
    max_pip = c(0.8, 0.6, 0.7, 0.9, 0.7, 0.85, 0.75, 0.8),
    median_pip = c(0.8, 0.6, 0.7, 0.9, 0.7, 0.85, 0.75, 0.8),
    stringsAsFactors = FALSE
  )

  expect_equal(merge_susie_cs(susie_fit_8), expected_output_8)

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
    credible_set_names = c("cs_1_1", "cs_1_1", "cs_1_2", "cs_1_2", "cs_1_3"),
    max_pip = c(0.8, 0.6, 0.7, 0.85, 0.9),
    median_pip = c(0.8, 0.6, 0.7, 0.85, 0.9),
    stringsAsFactors = FALSE
  )

  expect_equal(merge_susie_cs(susie_fit_9), expected_output_9)

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
    credible_set_names = c("cs_1_1,cs_3_3", "cs_1_1,cs_3_3", "cs_1_2,cs_2_2", "cs_2_1", "cs_2_3,cs_3_2", "cs_3_1"),
    max_pip = c(0.8, 0.7, 0.75, 0.9, 0.85, 0.7),
    median_pip = c(0.8, 0.675, 0.675, 0.9, 0.825, 0.7),
    stringsAsFactors = FALSE
  )

  expect_equal(merge_susie_cs(susie_fit_10), expected_output_10)
})

test_that("merge_susie_cs handles single condition with single CS", {
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

  result <- pecotmr:::merge_susie_cs(susie_fit)
  expect_s3_class(result, "data.frame")
  expect_true("variant_id" %in% colnames(result))
  expect_true("max_pip" %in% colnames(result))
})

# ===========================================================================
# filter_invalid_summary_stat
# ===========================================================================

test_that("filter_invalid_summary_stat replaces NaN/Inf in bhat", {
  dat <- list(
    bhat = data.frame(a = c(1, NaN, 3), b = c(Inf, 2, -Inf)),
    sbhat = data.frame(a = c(0.1, 0.2, 0.3), b = c(0.1, NA, 0.3))
  )
  result <- filter_invalid_summary_stat(dat, bhat = "bhat", sbhat = "sbhat")
  expect_true(all(!is.nan(result$bhat)))
  expect_true(all(!is.infinite(result$bhat)))
  # NaN/Inf in bhat replaced with 0
  expect_equal(unname(result$bhat[1, 2]), 0)  # Inf -> 0
})

test_that("filter_invalid_summary_stat replaces NaN/Inf in sbhat", {
  dat <- list(
    bhat = data.frame(a = c(1, 2, 3)),
    sbhat = data.frame(a = c(0.1, NaN, Inf))
  )
  result <- filter_invalid_summary_stat(dat, bhat = "bhat", sbhat = "sbhat")
  # NaN/Inf in sbhat replaced with 1000
  expect_equal(unname(result$sbhat[1, "a"]), 0.1)
  expect_equal(unname(result$sbhat[2, "a"]), 1000)
  expect_equal(unname(result$sbhat[3, "a"]), 1000)
})

test_that("filter_invalid_summary_stat filters by missing_rate when null.b present", {
  dat <- list(
    bhat = data.frame(a = c(0, 0, 1, 2), b = c(0, 0, 0, 3)),
    sbhat = data.frame(a = c(1, 1, 1, 1), b = c(1, 1, 1, 1)),
    null.b = TRUE
  )
  result <- filter_invalid_summary_stat(dat, bhat = "bhat", sbhat = "sbhat",
                                        filter_by_missing_rate = 0.5)
  expect_equal(nrow(result$bhat), 2) # rows 3 and 4 survive
})

test_that("filter_invalid_summary_stat filters by missing_rate when random.b present", {
  dat <- list(
    bhat = data.frame(a = c(0, 1, 2), b = c(0, 1, 3)),
    sbhat = data.frame(a = c(1, 1, 1), b = c(1, 1, 1)),
    random.b = TRUE
  )
  result <- filter_invalid_summary_stat(dat, bhat = "bhat", sbhat = "sbhat",
                                        filter_by_missing_rate = 0.5)
  expect_true(nrow(result$bhat) < 3)
})

test_that("filter_invalid_summary_stat btoz with .b and .s pattern creates condition.z", {
  dat <- list(
    strong.b = data.frame(a = c(1, 2), b = c(3, 4)),
    strong.s = data.frame(a = c(0.1, 0.2), b = c(0.3, 0.4))
  )
  result <- filter_invalid_summary_stat(dat, bhat = "strong.b", sbhat = "strong.s",
                                        btoz = TRUE, sig_p_cutoff = NULL)
  expect_true("strong.z" %in% names(result))
  expect_true(is.matrix(result$strong.z))
  expect_equal(nrow(result$strong.z), 2)
  expect_equal(ncol(result$strong.z), 2)
  expect_equal(as.numeric(result$strong.z), c(10, 10, 10, 10), tolerance = 1e-10)
})

test_that("filter_invalid_summary_stat btoz when bhat/sbhat data is NULL creates NULL z", {
  dat <- list(
    strong.b = NULL,
    strong.s = data.frame(a = c(0.1))
  )
  result <- filter_invalid_summary_stat(dat, bhat = "strong.b", sbhat = "strong.s",
                                        btoz = TRUE, sig_p_cutoff = NULL)
  expect_true("strong.z" %in% names(result))
  expect_null(result$strong.z)
})

test_that("filter_invalid_summary_stat btoz without .b/.s pattern creates generic z", {
  dat <- list(
    bhat = data.frame(a = c(1, 2, 3)),
    sbhat = data.frame(a = c(0.5, 1, 0.5))
  )
  result <- filter_invalid_summary_stat(dat, bhat = "bhat", sbhat = "sbhat",
                                        btoz = TRUE)
  expect_true("z" %in% names(result))
  expect_equal(as.numeric(result$z[, 1]), c(2, 2, 6))
})

test_that("filter_invalid_summary_stat btoz creates NULL z when bhat is NULL (no .b suffix)", {
  dat <- list(
    bhat = NULL,
    sbhat = data.frame(a = c(0.1))
  )
  result <- filter_invalid_summary_stat(dat, bhat = "bhat", sbhat = "sbhat",
                                        btoz = TRUE)
  expect_true("z" %in% names(result))
  expect_null(result$z)
})

test_that("filter_invalid_summary_stat btoz filters strong.z by significance cutoff", {
  dat <- list(
    strong.b = data.frame(a = c(1, 0.01, 0.02, 2), b = c(0.01, 0.01, 0.01, 0.01)),
    strong.s = data.frame(a = c(0.1, 0.1, 0.1, 0.1), b = c(0.1, 0.1, 0.1, 0.1)),
    strong.z = NULL
  )
  result <- filter_invalid_summary_stat(dat, bhat = "strong.b", sbhat = "strong.s",
                                        btoz = TRUE, sig_p_cutoff = 1E-6)
  expect_true("strong.z" %in% names(result))
  expect_equal(nrow(result$strong.z), 2)
  expect_equal(nrow(result$strong.b), 2)
  expect_equal(nrow(result$strong.s), 2)
})

test_that("filter_invalid_summary_stat processes z directly with null component", {
  dat <- list(
    strong = list(z = data.frame(a = c(5, NaN, 0.1), b = c(1, 2, Inf))),
    random = list(z = data.frame(a = c(0.5, 0.2), b = c(0.3, 0.4))),
    null = list(z = data.frame(a = c(0.01, NaN), b = c(Inf, 0.02)))
  )
  result <- filter_invalid_summary_stat(dat, z = "z")
  expect_true(all(!is.nan(result$strong$z)))
  expect_true(all(!is.infinite(result$strong$z)))
  expect_true(all(!is.nan(result$null$z)))
  expect_true(all(!is.infinite(result$null$z)))
})

test_that("filter_invalid_summary_stat z path applies significance cutoff to strong.z", {
  dat <- list(
    strong = list(z = data.frame(a = c(10, 0.1, 0.2), b = c(0.1, 0.1, 0.1))),
    random = list(z = data.frame(a = c(0.5, 0.2, 0.3), b = c(0.3, 0.4, 0.2)))
  )
  result <- filter_invalid_summary_stat(dat, z = "z", sig_p_cutoff = 1E-6)
  expect_equal(nrow(result$strong$z), 1)
})

test_that("filter_invalid_summary_stat z path with filter_by_missing_rate", {
  dat <- list(
    random = list(z = data.frame(a = c(0, 0, 5), b = c(0, 3, 4)))
  )
  result <- filter_invalid_summary_stat(dat, z = "z", filter_by_missing_rate = 0.5)
  expect_equal(nrow(result$random$z), 2)
})

test_that("filter_invalid_summary_stat processes bhat/sbhat without filter_by_missing_rate when no null.b/random.b", {
  dat <- list(
    bhat = data.frame(a = c(1, NaN, 3), b = c(Inf, 2, -Inf)),
    sbhat = data.frame(a = c(0.1, NA, 0.3), b = c(0.1, 0.2, NaN))
  )
  result <- filter_invalid_summary_stat(dat, bhat = "bhat", sbhat = "sbhat",
                                        filter_by_missing_rate = 0.5)
  expect_equal(nrow(result$bhat), 3)
  expect_equal(unname(result$bhat[2, "a"]), 0)
  expect_equal(unname(result$sbhat[3, "b"]), 1000)
})

test_that("filter_invalid_summary_stat with filter_by_missing_rate=NULL keeps all rows even with null.b", {
  dat <- list(
    bhat = data.frame(a = c(0, 0, 1), b = c(0, 0, 1)),
    sbhat = data.frame(a = c(1, 1, 1), b = c(1, 1, 1)),
    null.b = TRUE
  )
  result <- filter_invalid_summary_stat(dat, bhat = "bhat", sbhat = "sbhat",
                                        filter_by_missing_rate = NULL)
  expect_equal(nrow(result$bhat), 3)
})

test_that("filter_invalid_summary_stat z path handles NULL strong component", {
  dat <- list(
    strong = NULL,
    random = list(z = data.frame(a = c(0.5, 0.2), b = c(0.3, 0.4)))
  )
  result <- filter_invalid_summary_stat(dat, z = "z")
  expect_null(result$strong)
  expect_true(!is.null(result$random$z))
})

# ===========================================================================
# filter_mixture_components
# ===========================================================================

test_that("filter_mixture_components removes zero matrices", {
  U <- list(
    comp1 = matrix(c(1, 0, 0, 2), 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
    comp2 = matrix(c(0, 0, 0, 0), 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
    A = matrix(c(3, 0, 0, 4), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  )
  w <- c(comp1 = 0.5, comp2 = 0.3, A = 0.2)
  result <- filter_mixture_components(c("A", "B"), U, w)
  expect_false("comp2" %in% names(result$U))
})

test_that("filter_mixture_components removes matrices below weight cutoff", {
  U <- list(
    comp1 = matrix(c(1, 0, 0, 2), 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
    comp2 = matrix(c(0.1, 0, 0, 0.1), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  )
  w <- c(comp1 = 0.9, comp2 = 0.00001)
  result <- filter_mixture_components(c("A", "B"), U, w, w_cutoff = 1e-4)
  expect_false("comp2" %in% names(result$U))
})

test_that("filter_mixture_components errors on missing conditions", {
  U <- list(
    comp1 = matrix(c(1, 0, 0, 2), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  )
  expect_error(
    filter_mixture_components(c("A", "C"), U),
    "not found in matrix"
  )
})

test_that("filter_mixture_components removes components named as filtered conditions", {
  U <- list(
    comp1 = matrix(c(1,0,0, 0,2,0, 0,0,3), 3, 3, dimnames = list(c("A","B","C"), c("A","B","C"))),
    C = matrix(c(4,0,0, 0,5,0, 0,0,6), 3, 3, dimnames = list(c("A","B","C"), c("A","B","C")))
  )
  w <- c(comp1 = 0.6, C = 0.4)
  result <- filter_mixture_components(c("A", "B"), U, w)
  expect_false("C" %in% names(result$U))
  expect_equal(nrow(result$U$comp1), 2)
  expect_equal(ncol(result$U$comp1), 2)
})

test_that("filter_mixture_components renormalizes weights to preserve original sum", {
  U <- list(
    comp1 = matrix(c(1, 0, 0, 2), 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
    comp2 = matrix(c(3, 0, 0, 4), 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
    comp3 = matrix(c(0, 0, 0, 0), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  )
  w <- c(comp1 = 0.5, comp2 = 0.3, comp3 = 0.2)
  original_sum <- sum(w)

  result <- filter_mixture_components(c("A", "B"), U, w)
  expect_false("comp3" %in% names(result$U))
  expect_equal(sum(result$w), original_sum, tolerance = 1e-10)
  expect_true(result$w["comp1"] > 0.5)
  expect_true(result$w["comp2"] > 0.3)
})

test_that("filter_mixture_components subsets 3x3 matrices to 2x2 and removes filtered condition names", {
  U <- list(
    comp1 = matrix(c(1, 0.1, 0, 0.1, 2, 0, 0, 0, 3), 3, 3,
                   dimnames = list(c("A", "B", "C"), c("A", "B", "C"))),
    B = matrix(c(4, 0, 0, 0, 5, 0, 0, 0, 6), 3, 3,
               dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  )
  w <- c(comp1 = 0.7, B = 0.3)

  result <- filter_mixture_components(c("A", "C"), U, w)
  expect_false("B" %in% names(result$U))
  expect_equal(nrow(result$U$comp1), 2)
  expect_equal(ncol(result$U$comp1), 2)
  expect_equal(rownames(result$U$comp1), c("A", "C"))
})

test_that("filter_mixture_components handles NULL weights gracefully", {
  U <- list(
    comp1 = matrix(c(1, 0, 0, 2), 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
    comp2 = matrix(c(0, 0, 0, 0), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  )
  result <- filter_mixture_components(c("A", "B"), U, w = NULL)
  expect_false("comp2" %in% names(result$U))
  expect_true("comp1" %in% names(result$U))
})

# ===========================================================================
# merge_mash_data
# ===========================================================================

test_that("merge_mash_data combines two datasets with identical columns", {
  d1 <- list(random = data.frame(a = 1:3, b = 4:6))
  d2 <- list(random = data.frame(a = 7:8, b = 9:10))
  result <- merge_mash_data(d1, d2)
  expect_equal(nrow(result$random), 5)
  expect_equal(ncol(result$random), 2)
  expect_equal(colnames(result$random), c("a", "b"))
  expect_equal(result$random$a, c(1, 2, 3, 7, 8))
})

test_that("merge_mash_data handles NULL input", {
  d1 <- NULL
  d2 <- list(random = data.frame(a = 1:3))
  result <- merge_mash_data(d1, d2)
  expect_equal(nrow(result$random), 3)
})

test_that("merge_mash_data aligns different column names correctly", {
  d1 <- list(random = data.frame(a = 1:2, b = 3:4))
  d2 <- list(random = data.frame(a = 5:6, c = 7:8))
  result <- merge_mash_data(d1, d2)
  expect_equal(nrow(result$random), 4)
  expect_true(all(c("a", "b", "c") %in% colnames(result$random)))
  expect_true(is.nan(result$random[3, "b"]))
  expect_true(is.nan(result$random[1, "c"]))
})

test_that("merge_mash_data preserves data when one side is empty data.frame", {
  d1 <- list(random = data.frame(a = 1:3))
  d2 <- list(random = data.frame())
  result <- merge_mash_data(d1, d2)
  expect_equal(nrow(result$random), 3)
})

test_that("merge_mash_data preserves data when one side is NULL element", {
  d1 <- list(random = NULL)
  d2 <- list(random = data.frame(a = 1:3))
  result <- merge_mash_data(d1, d2)
  expect_equal(nrow(result$random), 3)
})

test_that("merge_mash_data handles multiple named elements", {
  d1 <- list(
    random = data.frame(a = 1:2, b = 3:4),
    null = data.frame(x = 10:11)
  )
  d2 <- list(
    random = data.frame(a = 5:6, b = 7:8),
    null = data.frame(x = 12:13)
  )
  result <- merge_mash_data(d1, d2)
  expect_equal(nrow(result$random), 4)
  expect_equal(nrow(result$null), 4)
})

test_that("merge_mash_data uses one_data when res_data element has zero rows", {
  d1 <- list(random = data.frame(a = numeric(0), b = numeric(0)))
  d2 <- list(random = data.frame(a = 1:3, b = 4:6))
  result <- merge_mash_data(d1, d2)
  expect_true(nrow(result$random) >= 3)
})

# ===========================================================================
# mash_rand_null_sample
# ===========================================================================

test_that("mash_rand_null_sample with z scores returns random and null", {
  set.seed(42)
  dat <- list(
    z = data.frame(
      cond1 = c(5, 0.1, 0.2, 0.3, 0.5, 6, 0.1, 0.2, 0.4, 0.3),
      cond2 = c(0.2, 0.3, 0.1, 0.5, 0.4, 0.1, 0.3, 0.2, 0.1, 0.5)
    )
  )
  result <- mash_rand_null_sample(dat, n_random = 5, n_null = 3,
                                   exclude_condition = c(), seed = 123)
  expect_type(result, "list")
  expect_true("random" %in% names(result))
  expect_true("null" %in% names(result))
  expect_true("z" %in% names(result$random))
  expect_equal(nrow(result$random$z), 5)
})

test_that("mash_rand_null_sample with seed is reproducible", {
  dat <- list(
    z = data.frame(
      cond1 = c(0.1, 0.2, 0.3, 0.4, 0.5),
      cond2 = c(0.5, 0.4, 0.3, 0.2, 0.1)
    )
  )
  result1 <- mash_rand_null_sample(dat, n_random = 3, n_null = 2,
                                    exclude_condition = c(), seed = 42)
  result2 <- mash_rand_null_sample(dat, n_random = 3, n_null = 2,
                                    exclude_condition = c(), seed = 42)
  expect_equal(result1$random$z, result2$random$z)
})

test_that("mash_rand_null_sample NULL input returns NULL", {
  result <- mash_rand_null_sample(NULL, n_random = 5, n_null = 3,
                                   exclude_condition = c())
  expect_null(result)
})

test_that("mash_rand_null_sample warns when no null variants found (all abs_z > 2)", {
  dat <- list(
    z = data.frame(
      cond1 = c(5, 6, 7, 8, 9),
      cond2 = c(5, 6, 7, 8, 9)
    )
  )
  expect_warning(
    result <- mash_rand_null_sample(dat, n_random = 3, n_null = 2,
                                     exclude_condition = c(), seed = 42),
    "no variants are included in the null"
  )
  expect_equal(length(result$null), 0)
})

test_that("mash_rand_null_sample warns when not enough null data", {
  dat <- list(
    z = data.frame(
      cond1 = c(5, 6, 0.1),
      cond2 = c(5, 6, 0.1),
      cond3 = c(5, 6, 0.1)
    )
  )
  expect_warning(
    result <- mash_rand_null_sample(dat, n_random = 2, n_null = 1,
                                     exclude_condition = c(), seed = 42),
    "not enough null data"
  )
  expect_equal(length(result$null), 0)
})

test_that("mash_rand_null_sample with bhat/sbhat processes random and null samples", {
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
  result <- mash_rand_null_sample(dat, n_random = 3, n_null = 3,
                                   exclude_condition = c(), seed = 42)
  expect_equal(ncol(result$random$bhat), 2)
  expect_equal(nrow(result$random$bhat), 3)
  expect_true(length(result$null) > 0)
})

test_that("mash_rand_null_sample errors when exclude_condition not found (z path)", {
  dat <- list(
    z = data.frame(cond1 = 1:5, cond2 = 1:5)
  )
  expect_error(
    mash_rand_null_sample(dat, n_random = 3, n_null = 2,
                           exclude_condition = "nonexistent", seed = 42),
    "exclude_condition are not present"
  )
})

test_that("mash_rand_null_sample errors when exclude_condition not found (bhat path)", {
  dat <- list(
    bhat = data.frame(cond1 = 1:5, cond2 = 1:5),
    sbhat = data.frame(cond1 = rep(1, 5), cond2 = rep(1, 5))
  )
  expect_error(
    mash_rand_null_sample(dat, n_random = 3, n_null = 2,
                           exclude_condition = "nonexistent", seed = 42),
    "exclude_condition are not present"
  )
})

test_that("mash_rand_null_sample exclude_condition with column names errors due to numeric indexing", {
  dat <- list(
    z = data.frame(
      cond1 = c(0.1, 0.2, 0.3, 0.4, 0.5),
      cond2 = c(0.5, 0.4, 0.3, 0.2, 0.1),
      cond3 = c(0.3, 0.3, 0.3, 0.3, 0.3)
    )
  )
  expect_error(
    mash_rand_null_sample(dat, n_random = 3, n_null = 2,
                           exclude_condition = "cond3", seed = 42),
    "invalid argument to unary operator"
  )
})

test_that("mash_rand_null_sample extracts null data with z scores when enough null variants exist", {
  dat <- list(
    z = data.frame(
      cond1 = c(0.1, 0.3, 0.2, 0.5, 0.4, 0.1, 0.3, 0.2, 0.5, 0.4),
      cond2 = c(0.2, 0.1, 0.4, 0.3, 0.5, 0.2, 0.1, 0.4, 0.3, 0.5)
    )
  )
  result <- mash_rand_null_sample(dat, n_random = 5, n_null = 4,
                                   exclude_condition = c(), seed = 42)
  expect_true("null" %in% names(result))
  expect_true("z" %in% names(result$null))
  expect_equal(nrow(result$null$z), 4)
  expect_equal(ncol(result$null$z), 2)
  expect_equal(nrow(result$random$z), 5)
})

test_that("mash_rand_null_sample null data capped at available null variants", {
  dat <- list(
    z = data.frame(
      cond1 = c(0.1, 0.2, 0.3),
      cond2 = c(0.1, 0.2, 0.3)
    )
  )
  result <- mash_rand_null_sample(dat, n_random = 2, n_null = 100,
                                   exclude_condition = c(), seed = 42)
  expect_true(length(result$null) > 0)
  expect_equal(nrow(result$null$z), 3)
})

test_that("mash_rand_null_sample extracts null data with bhat/sbhat when enough null variants", {
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
  result <- mash_rand_null_sample(dat, n_random = 3, n_null = 4,
                                   exclude_condition = c(), seed = 42)
  expect_true("null" %in% names(result))
  expect_true("bhat" %in% names(result$null))
  expect_true("sbhat" %in% names(result$null))
  expect_equal(nrow(result$null$bhat), 4)
  expect_equal(nrow(result$null$sbhat), 4)
})

test_that("mash_rand_null_sample exclude_condition with numeric index errors on z path", {
  dat <- list(
    z = data.frame(
      cond1 = c(0.1, 0.2, 0.3, 0.4, 0.5),
      cond2 = c(0.5, 0.4, 0.3, 0.2, 0.1),
      cond3 = c(0.3, 0.3, 0.3, 0.3, 0.3)
    )
  )
  expect_error(
    mash_rand_null_sample(dat, n_random = 3, n_null = 2,
                           exclude_condition = 3, seed = 42),
    "exclude_condition are not present"
  )
})

test_that("mash_rand_null_sample exclude_condition with numeric index errors on bhat path", {
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
    mash_rand_null_sample(dat, n_random = 3, n_null = 2,
                           exclude_condition = 3, seed = 42),
    "exclude_condition are not present"
  )
})

test_that("mash_rand_null_sample caps random sample at available rows", {
  dat <- list(
    z = data.frame(
      cond1 = c(0.1, 0.2, 0.3),
      cond2 = c(0.2, 0.1, 0.3)
    )
  )
  result <- mash_rand_null_sample(dat, n_random = 100, n_null = 2,
                                   exclude_condition = c(), seed = 42)
  expect_equal(nrow(result$random$z), 3)
})

# ===========================================================================
# mash_pipeline
# ===========================================================================

test_that("mash_pipeline uses residual_correlation when null data is empty", {
  skip_if_not_installed("mashr")
  skip_if_not_installed("flashier")

  n_cond <- 3
  n_strong <- 15
  mock_input <- list(
    null.b = numeric(0),
    null.s = numeric(0),
    random.b = matrix(rnorm(n_strong * n_cond), n_strong, n_cond),
    random.s = matrix(abs(rnorm(n_strong * n_cond)) + 0.1, n_strong, n_cond),
    strong.b = matrix(rnorm(n_strong * n_cond), n_strong, n_cond),
    strong.s = matrix(abs(rnorm(n_strong * n_cond)) + 0.1, n_strong, n_cond)
  )
  custom_vhat <- matrix(c(1, 0.3, 0.1, 0.3, 1, 0.2, 0.1, 0.2, 1), 3, 3)

  result <- mash_pipeline(mock_input, alpha = 1, residual_correlation = custom_vhat)
  expect_true(is.list(result))
  expect_null(result$error)
})

test_that("mash_pipeline does not error about mashr when mashr is available", {
  skip_if_not_installed("mashr")
  skip_if_not_installed("flashier")

  n_cond <- 3
  n_strong <- 15
  mock_input <- list(
    null.b = numeric(0),
    null.s = numeric(0),
    random.b = matrix(rnorm(n_strong * n_cond), n_strong, n_cond),
    random.s = matrix(abs(rnorm(n_strong * n_cond)) + 0.1, n_strong, n_cond),
    strong.b = matrix(rnorm(n_strong * n_cond), n_strong, n_cond),
    strong.s = matrix(abs(rnorm(n_strong * n_cond)) + 0.1, n_strong, n_cond)
  )
  result <- mash_pipeline(mock_input, alpha = 1)
  expect_true(is.list(result))
})

test_that("mash_pipeline uses identity matrix when null data empty and no residual_correlation", {
  skip_if_not_installed("mashr")
  skip_if_not_installed("flashier")

  n_cond <- 3
  n_strong <- 15
  mock_input <- list(
    null.b = numeric(0),
    null.s = numeric(0),
    random.b = matrix(rnorm(n_strong * n_cond), n_strong, n_cond),
    random.s = matrix(abs(rnorm(n_strong * n_cond)) + 0.1, n_strong, n_cond),
    strong.b = matrix(rnorm(n_strong * n_cond), n_strong, n_cond),
    strong.s = matrix(abs(rnorm(n_strong * n_cond)) + 0.1, n_strong, n_cond)
  )

  result <- mash_pipeline(mock_input, alpha = 1, residual_correlation = NULL)
  expect_true(is.list(result))
  expect_null(result$error)
})

# ===========================================================================
# merge_sumstats_matrices
# ===========================================================================

test_that("merge_sumstats_matrices validates empty list", {
  expect_error(merge_sumstats_matrices(list(), "z"), "non-empty list")
})

test_that("merge_sumstats_matrices validates value_column is single string", {
  expect_error(merge_sumstats_matrices(list(data.frame(a = 1)), c("a", "b")), "single string")
})

test_that("merge_sumstats_matrices validates id_column is single string", {
  expect_error(merge_sumstats_matrices(list(data.frame(a = 1)), "a", id_column = 123), "single string")
})

test_that("merge_sumstats_matrices merges multiple data frames", {
  m1 <- data.frame(variants = c("chr1:100:A:G", "chr1:200:C:T"), z = c(1.5, -2.0))
  m2 <- data.frame(variants = c("chr1:100:A:G", "chr1:300:G:A"), z = c(2.0, 0.5))
  result <- merge_sumstats_matrices(list(m1, m2), value_column = "z")
  expect_true("variants" %in% colnames(result))
  expect_equal(nrow(result), 3)
})

test_that("merge_sumstats_matrices with remove_any_missing=TRUE drops incomplete rows", {
  m1 <- data.frame(variants = c("chr1:100:A:G", "chr1:200:C:T"), z = c(1.5, -2.0))
  m2 <- data.frame(variants = c("chr1:100:A:G"), z = c(2.0))
  result <- merge_sumstats_matrices(list(m1, m2), value_column = "z", remove_any_missing = TRUE)
  expect_equal(nrow(result), 1)
})

test_that("merge_sumstats_matrices with remove_any_missing=FALSE keeps incomplete rows", {
  m1 <- data.frame(variants = c("chr1:100:A:G", "chr1:200:C:T"), z = c(1.5, -2.0))
  m2 <- data.frame(variants = c("chr1:100:A:G"), z = c(2.0))
  result <- merge_sumstats_matrices(list(m1, m2), value_column = "z",
                                    remove_any_missing = FALSE)
  expect_equal(nrow(result), 2)
})

test_that("merge_sumstats_matrices handles error in processing gracefully", {
  m1 <- data.frame(variants = c("chr1:100:A:G"), z = c(1.5))
  m2 <- data.frame(wrong_col = c("chr1:200:C:T"), z = c(2.0))
  result <- merge_sumstats_matrices(list(m1, m2), value_column = "z")
  expect_true(!is.null(result))
})

test_that("merge_sumstats_matrices returns NULL when all datasets fail processing", {
  m1 <- data.frame(wrong1 = c("a"), wrong2 = c(1))
  m2 <- data.frame(wrong3 = c("b"), wrong4 = c(2))

  expect_message(
    result <- merge_sumstats_matrices(list(m1, m2), value_column = "z"),
    "Error processing dataset"
  )
  expect_null(result)
})

test_that("merge_sumstats_matrices handles mix of valid and invalid datasets", {
  m1 <- data.frame(variants = c("chr1:100:A:G"), z = c(1.5))
  m2 <- data.frame(wrong_col = c("chr1:200:C:T"), z = c(2.0))
  m3 <- data.frame(variants = c("chr1:100:A:G", "chr1:300:G:A"), z = c(3.0, 0.5))

  result <- merge_sumstats_matrices(list(m1, m2, m3), value_column = "z")
  expect_true(!is.null(result))
  expect_true("variants" %in% colnames(result))
  expect_equal(nrow(result), 2)
})

test_that("merge_sumstats_matrices with single valid dataset returns properly", {
  m1 <- data.frame(variants = c("chr1:100:A:G", "chr1:200:C:T"), z = c(1.5, -2.0))
  result <- merge_sumstats_matrices(list(m1), value_column = "z")
  expect_equal(nrow(result), 2)
  expect_true("z_1" %in% colnames(result))
})

# ===========================================================================
# extract_flatten_sumstats_from_nested
# ===========================================================================

test_that("extract_flatten_sumstats_from_nested computes z from betahat/sebetahat", {
  data <- list(
    finemapping_result = .test_fm_result(c("1:100:A:G", "1:200:C:T")),
    sumstats = list(
      betahat = c(0.5, -0.3),
      sebetahat = c(0.1, 0.15)
    )
  )
  result <- extract_flatten_sumstats_from_nested(data, extract_inf = "z")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_true("z" %in% colnames(result))
  expect_equal(result$z[1], 5.0)
  expect_equal(result$z[2], -2.0)
})

test_that("extract_flatten_sumstats_from_nested uses z directly when available", {
  data <- list(
    finemapping_result = .test_fm_result(c("1:100:A:G")),
    sumstats = list(z = c(3.5))
  )
  result <- extract_flatten_sumstats_from_nested(data, extract_inf = "z")
  expect_equal(result$z[1], 3.5)
})

test_that("extract_flatten_sumstats_from_nested extracts beta from direct sumstats", {
  data <- list(
    finemapping_result = .test_fm_result(c("chr1:100:A:G", "chr1:200:C:T")),
    sumstats = list(
      betahat = c(0.5, -0.3),
      sebetahat = c(0.1, 0.15)
    )
  )
  result <- extract_flatten_sumstats_from_nested(data, extract_inf = "beta")
  expect_equal(result$beta, c(0.5, -0.3))
})

test_that("extract_flatten_sumstats_from_nested extracts se from direct sumstats", {
  data <- list(
    finemapping_result = .test_fm_result(c("chr1:100:A:G", "chr1:200:C:T")),
    sumstats = list(
      betahat = c(0.5, -0.3),
      sebetahat = c(0.1, 0.15)
    )
  )
  result <- extract_flatten_sumstats_from_nested(data, extract_inf = "se")
  expect_equal(result$se, c(0.1, 0.15))
})

test_that("extract_flatten_sumstats_from_nested reaches max_depth and returns NULL", {
  data <- list(level1 = list(level2 = list(level3 = list(level4 = list(
    finemapping_result = .test_fm_result(c("1:100:A:G")),
    sumstats = list(z = c(2.0))
  )))))
  result <- extract_flatten_sumstats_from_nested(data, extract_inf = "z", max_depth = 2)
  expect_null(result)
})

test_that("extract_flatten_sumstats_from_nested handles missing betahat for z", {
  data <- list(
    finemapping_result = .test_fm_result(c("1:100:A:G")),
    sumstats = list(something_else = c(1.0))
  )
  result <- expect_message(
    extract_flatten_sumstats_from_nested(data, extract_inf = "z"),
    "Cannot compute z"
  )
  expect_null(result)
})

test_that("extract_flatten_sumstats_from_nested handles missing betahat for beta", {
  data <- list(
    finemapping_result = .test_fm_result(c("1:100:A:G")),
    sumstats = list(z = c(2.0))
  )
  result <- expect_message(
    extract_flatten_sumstats_from_nested(data, extract_inf = "beta"),
    "Missing 'betahat'"
  )
  expect_null(result)
})

test_that("extract_flatten_sumstats_from_nested handles missing sebetahat for se", {
  data <- list(
    finemapping_result = .test_fm_result(c("1:100:A:G")),
    sumstats = list(betahat = c(0.5))
  )
  result <- expect_message(
    extract_flatten_sumstats_from_nested(data, extract_inf = "se"),
    "Missing 'sebetahat'"
  )
  expect_null(result)
})

test_that("extract_flatten_sumstats_from_nested rejects invalid extract_inf values", {
  data <- list(
    finemapping_result = .test_fm_result(c("1:100:A:G")),
    sumstats = list(z = c(1.0))
  )
  expect_error(
    extract_flatten_sumstats_from_nested(data, extract_inf = "invalid_type"),
    "must be one of"
  )
})

test_that("extract_flatten_sumstats_from_nested normalizes variant IDs to chr prefix", {
  data <- list(
    finemapping_result = .test_fm_result(c("1:100:A:G", "2:200:C:T")),
    sumstats = list(z = c(1.0, 2.0))
  )
  result <- extract_flatten_sumstats_from_nested(data, extract_inf = "z")
  expect_true(all(grepl("^chr", result$variants)))
})

test_that("extract_flatten_sumstats_from_nested normalizes variant IDs from nested search", {
  data <- list(
    nested = list(
      finemapping_result = .test_fm_result(c("1:100:A:G")),
      sumstats = list(z = c(3.0))
    )
  )
  result <- extract_flatten_sumstats_from_nested(data, extract_inf = "z")
  expect_true(grepl("^chr", result$variants[1]))
})

test_that("extract_flatten_sumstats_from_nested recurses through multiple nesting levels", {
  data <- list(
    level1 = list(
      level2 = list(
        finemapping_result = .test_fm_result(c("chr1:100:A:G", "chr1:200:C:T")),
        sumstats = list(
          betahat = c(0.5, -0.3),
          sebetahat = c(0.1, 0.15)
        )
      )
    )
  )
  result <- extract_flatten_sumstats_from_nested(data, extract_inf = "z", max_depth = 3)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_equal(result$z[1], 5.0, tolerance = 1e-10)
})

test_that("extract_flatten_sumstats_from_nested returns NULL for deeply nested beyond max_depth", {
  data <- list(
    a = list(
      b = list(
        c = list(
          finemapping_result = .test_fm_result(c("1:100:A:G")),
          sumstats = list(z = c(2.0))
        )
      )
    )
  )
  expect_message(
    result <- extract_flatten_sumstats_from_nested(data, extract_inf = "z", max_depth = 2),
    "Maximum search depth reached"
  )
  expect_null(result)
})

test_that("extract_flatten_sumstats_from_nested returns NULL when not found", {
  data <- list(a = list(b = 42))
  result <- extract_flatten_sumstats_from_nested(data, extract_inf = "z")
  expect_null(result)
})

# ===========================================================================
# load_multitrait_tensorqtl_sumstat: input validation
# ===========================================================================

test_that("load_multitrait_tensorqtl_sumstat errors when sumstats_paths files do not exist", {
  expect_error(
    load_multitrait_tensorqtl_sumstat(
      sumstats_paths = c("/nonexistent/file1.txt", "/nonexistent/file2.txt"),
      region = "chr1:1000-2000",
      trait_names = c("trait1", "trait2")
    ),
    "sumstats_paths must be a vector of existing file paths"
  )
})

test_that("load_multitrait_tensorqtl_sumstat errors when sumstats_paths is a list (not vector)", {
  expect_error(
    load_multitrait_tensorqtl_sumstat(
      sumstats_paths = list("a", "b"),
      region = "chr1:1000-2000",
      trait_names = c("trait1", "trait2")
    )
  )
})

test_that("load_multitrait_tensorqtl_sumstat errors when region is not single character", {
  f1 <- tempfile()
  f2 <- tempfile()
  file.create(f1, f2)
  on.exit(unlink(c(f1, f2)))

  expect_error(
    load_multitrait_tensorqtl_sumstat(
      sumstats_paths = c(f1, f2),
      region = 12345,
      trait_names = c("trait1", "trait2")
    ),
    "region must be a single character string"
  )

  expect_error(
    load_multitrait_tensorqtl_sumstat(
      sumstats_paths = c(f1, f2),
      region = c("chr1:1-100", "chr2:1-100"),
      trait_names = c("trait1", "trait2")
    ),
    "region must be a single character string"
  )
})

test_that("load_multitrait_tensorqtl_sumstat errors when trait_names is not character", {
  f1 <- tempfile()
  f2 <- tempfile()
  file.create(f1, f2)
  on.exit(unlink(c(f1, f2)))

  expect_error(
    load_multitrait_tensorqtl_sumstat(
      sumstats_paths = c(f1, f2),
      region = "chr1:1000-2000",
      trait_names = c(1, 2)
    ),
    "trait_names must be a vector of character strings"
  )
})
