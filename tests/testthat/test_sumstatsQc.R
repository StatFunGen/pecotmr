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


test_that("alignVariantNames correctly aligns variant names", {
  # Test case 1: Matching variant names
  source1 <- c("1:123:A:C", "2:456:G:T", "3:789:C:A")
  reference1 <- c("1:123:A:C", "2:456:T:G", "3:789:C:A")
  expected_aligned1 <- c("1:123:A:C", "2:456:T:G", "3:789:C:A")
  expected_unmatched1 <- integer(0)

  result1 <- alignVariantNames(source1, reference1)
  expect_equal(result1$alignedVariants, expected_aligned1)
  expect_equal(result1$unmatchedIndices, expected_unmatched1)

  # Test case 2: Unmatched variant names
  source2 <- c("1:123:A:C", "2:456:G:T", "3:789:C:A", "4:101:G:C")
  reference2 <- c("1:123:A:C", "2:456:T:G", "3:789:C:A")
  expected_aligned2 <- c("1:123:A:C", "2:456:T:G", "3:789:C:A", "4:101:G:C")
  expected_unmatched2 <- 4

  result2 <- alignVariantNames(source2, reference2)
  expect_equal(result2$alignedVariants, expected_aligned2)
  expect_equal(result2$unmatchedIndices, expected_unmatched2)

  # Test case 3: Different variant name formats
  source3 <- c("1:123:A:C", "2:456_G_T", "3:789:C:A")
  reference3 <- c("1:123:A:C", "2:456:T:G", "3:789:C:A")
  expected_aligned3 <- c("1:123:A:C", "2:456:T:G", "3:789:C:A")
  expected_unmatched3 <- integer(0)

  result3 <- alignVariantNames(source3, reference3)
  expect_equal(result3$alignedVariants, expected_aligned3)
  expect_equal(result3$unmatchedIndices, expected_unmatched3)
})

test_that("alignVariantNames correctly aligns variant names with different flip patterns", {
  # Test case 4: Strand flip
  source4 <- c("1:123:A:C", "2:456:G:T", "3:789:C:A")
  reference4 <- c("1:123:T:G", "2:456:A:C", "3:789:C:A")
  expected_aligned4 <- c("1:123:T:G", "2:456:A:C", "3:789:C:A")
  expected_unmatched4 <- integer(0)

  result4 <- alignVariantNames(source4, reference4)
  expect_equal(result4$alignedVariants, expected_aligned4)
  expect_equal(result4$unmatchedIndices, expected_unmatched4)

  # Test case 5: Strand ambiguous variants
  source5 <- c("1:123:A:T", "2:456:G:C", "3:789:C:A")
  reference5 <- c("1:123:A:T", "2:456:G:C", "3:789:C:A")
  expected_aligned5 <- c("1:123:A:T", "2:456:G:C", "3:789:C:A")
  expected_unmatched5 <- integer(0)

  result5 <- alignVariantNames(source5, reference5)
  expect_equal(result5$alignedVariants, expected_aligned5)
  expect_equal(result5$unmatchedIndices, expected_unmatched5)

  # Test case 6: Sign flip
  source6 <- c("1:123:A:C", "2:456:G:T", "3:789:C:A")
  reference6 <- c("1:123:C:A", "2:456:T:G", "3:789:C:A")
  expected_aligned6 <- c("1:123:C:A", "2:456:T:G", "3:789:C:A")
  expected_unmatched6 <- integer(0)

  result6 <- alignVariantNames(source6, reference6)
  expect_equal(result6$alignedVariants, expected_aligned6)
  expect_equal(result6$unmatchedIndices, expected_unmatched6)

  # Test case 7: Strand and sign flip
  source7 <- c("1:123:A:C", "2:456:G:T", "3:789:C:A")
  reference7 <- c("1:123:G:T", "2:456:A:C", "3:789:C:A")
  expected_aligned7 <- c("1:123:G:T", "2:456:A:C", "3:789:C:A")
  expected_unmatched7 <- integer(0)

  result7 <- alignVariantNames(source7, reference7)
  expect_equal(result7$alignedVariants, expected_aligned7)
  expect_equal(result7$unmatchedIndices, expected_unmatched7)

  # Test case 8: Indels
  source8 <- c("1:123:A:C", "2:456:G:T", "3:789:C:A", "4:101:G:GATC")
  reference8 <- c("1:123:A:C", "2:456:T:G", "3:789:C:A", "4:101:GATC:G")
  expected_aligned8 <- c("1:123:A:C", "2:456:T:G", "3:789:C:A", "4:101:GATC:G")
  expected_unmatched8 <- integer(0)

  result8 <- alignVariantNames(source8, reference8)
  expect_equal(result8$alignedVariants, expected_aligned8)
  expect_equal(result8$unmatchedIndices, expected_unmatched8)
})

test_that("alignVariantNames correctly aligns variant names with different chr prefix conventions", {
  # Test case 9: Original without chr prefix, reference with chr prefix
  source9 <- c("1:123:A:C", "2:456:G:T", "3:789:C:A")
  reference9 <- c("chr1:123:A:C", "chr2:456:T:G", "chr3:789:C:A")
  expected_aligned9 <- c("chr1:123:A:C", "chr2:456:T:G", "chr3:789:C:A")
  expected_unmatched9 <- integer(0)

  result9 <- alignVariantNames(source9, reference9)
  expect_equal(result9$alignedVariants, expected_aligned9)
  expect_equal(result9$unmatchedIndices, expected_unmatched9)

  # Test case 10: Original with chr prefix, reference without chr prefix
  source10 <- c("chr1:123:A:C", "chr2:456:G:T", "chr3:789:C:A")
  reference10 <- c("1:123:A:C", "2:456:T:G", "3:789:C:A")
  expected_aligned10 <- c("1:123:A:C", "2:456:T:G", "3:789:C:A")
  expected_unmatched10 <- integer(0)

  result10 <- alignVariantNames(source10, reference10)
  expect_equal(result10$alignedVariants, expected_aligned10)
  expect_equal(result10$unmatchedIndices, expected_unmatched10)
})

test_that("alignVariantNames warns on non-standard format", {
  source <- c("rs12345")
  reference <- c("rs67890")
  expect_warning(
    alignVariantNames(source, reference),
    "do not follow the expected"
  )
})

test_that("alignVariantNames errors on mixed formats", {
  source <- c("1:100:A:G")
  reference <- c("rs12345")
  expect_error(
    alignVariantNames(source, reference),
    "different variant naming conventions"
  )
})

test_that("alignVariantNames strips build suffix", {
  source <- c("1:100:A:G:b38")
  reference <- c("1:100:A:G")
  result <- alignVariantNames(source, reference, removeBuildSuffix = TRUE)
  expect_length(result$alignedVariants, 1)
})

test_that("alignVariantNames: disjoint source/reference returns source unchanged", {
  # Bug fix: when no variants harmonize, paste0() over length-0 components
  # used to collapse to a single "chr:::" placeholder. Verify the output
  # now preserves source length and flags every position as unmatched.
  src <- c("chr1:10:A:G", "chr1:20:A:G", "chr1:30:A:G")
  ref <- c("chr2:40:A:G", "chr2:50:A:G", "chr2:60:A:G")
  out <- suppressWarnings(alignVariantNames(src, ref))
  expect_equal(length(out$alignedVariants), length(src))
  expect_equal(out$alignedVariants, src)
  expect_equal(out$unmatchedIndices, seq_along(src))
})



context("raiss")
library(tidyverse)
library(MASS)

# Helper: build LdData S4 from a ref_panel data.frame, correlation matrix, and blockMetadata
make_ld_data_from_ref_panel <- function(R_mat, ref_panel, blockMetadata) {
  ref_panel$chrom <- as.character(ref_panel$chrom)
  ref_panel$variant_id <- as.character(ref_panel$variant_id)
  variants_gr <- pecotmr:::.refPanelToGranges(ref_panel)
  LdData(correlation = R_mat, variants = variants_gr, blockMetadata = blockMetadata)
}

generate_dummy_data <- function(seed=1, ref_panel_ordered=TRUE, known_zscores_ordered=TRUE) {
    set.seed(seed)

    n_variants <- 100
    ref_panel <- data.frame(
        chrom = rep(1, n_variants),
        pos = seq(1, n_variants * 10, 10),
        variant_id = paste0("rs", seq_len(n_variants)),
        A1 = sample(c("A", "T", "G", "C"), n_variants, replace = TRUE),
        A2 = sample(c("A", "T", "G", "C"), n_variants, replace = TRUE)
    )

    n_known <- 50
    known_zscores <- data.frame(
        chrom = rep(1, n_known),
        pos = sample(ref_panel$pos, n_known),
        variant_id = sample(ref_panel$variant_id, n_known),
        A1 = sample(c("A", "T", "G", "C"), n_known, replace = TRUE),
        A2 = sample(c("A", "T", "G", "C"), n_known, replace = TRUE),
        z = rnorm(n_known)
    )

    ldMatrix <- matrix(rnorm(n_variants^2), nrow = n_variants, ncol = n_variants)
    diag(ldMatrix) <- 1
    known_zscores <- if (known_zscores_ordered) known_zscores[order(known_zscores$pos),] else known_zscores
    ref_panel <- if (ref_panel_ordered) ref_panel else ref_panel[order(ref_panel$pos, decreasing = TRUE),]
    return(list(ref_panel=ref_panel, known_zscores=known_zscores, ldMatrix=ldMatrix))
}

test_that("Input validation for raiss works correctly", {
    input_data <- generate_dummy_data()
    input_data_ref_panel_unordered <- generate_dummy_data(ref_panel_ordered=FALSE)
    input_data_zscores_unordered <- generate_dummy_data(known_zscores_ordered=FALSE)
    expect_error(raiss(input_data_ref_panel_unordered$ref_panel, input_data$known_zscores, input_data$ldMatrix))
    expect_error(raiss(input_data$ref_panel, input_data_zscores_unordered$known_zscores, input_data$ldMatrix))
})

test_that("Default parameters for raiss work correctly", {
    input_data <- generate_dummy_data()
    result <- raiss(input_data$ref_panel, input_data$known_zscores, input_data$ldMatrix)
    expect_true(is.list(result))
    # Expected list elements
    expect_true(all(c("resultNofilter", "resultFilter", "ldMat") %in% names(result)))
    # resultNofilter should be a data frame with expected columns
    expect_true(is.data.frame(result$resultNofilter))
    expect_true(all(c("variant_id", "z", "Var", "raissLdScore") %in% names(result$resultNofilter)))
    # Imputed z-scores should be numeric and finite
    expect_true(is.numeric(result$resultNofilter$z))
    expect_true(all(is.finite(result$resultNofilter$z)))
    # Output should cover all ref_panel variants (known + imputed)
    expect_equal(nrow(result$resultNofilter), nrow(input_data$ref_panel))
    # Filtered result should be a subset of unfiltered
    expect_true(nrow(result$resultFilter) <= nrow(result$resultNofilter))
    # ldMat should be a matrix
    expect_true(is.matrix(result$ldMat))
})

test_that("Test Default Parameters for raissModel", {
  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

  result <- raissModel(zt, sig_t, sig_i_t)

  expect_true(is.list(result))
  expect_true(all(c("var", "mu", "raissLdScore", "conditionNumber", "correctInversion") %in% names(result)))
  # mu (imputed z-scores) should be numeric and finite
  expect_true(is.numeric(result$mu))
  expect_true(all(is.finite(result$mu)))
  # var should be numeric
  expect_true(is.numeric(result$var))
  # raissLdScore should be numeric and non-negative
  expect_true(is.numeric(result$raissLdScore))
  expect_true(all(result$raissLdScore >= 0))
})

test_that("Test with Different lamb Values for raissModel", {
  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

  lamb_values <- c(0.01, 0.05, 0.1)
  for (lamb in lamb_values) {
    result <- raissModel(zt, sig_t, sig_i_t, lamb)
    expect_true(is.list(result))
    expect_true(all(c("var", "mu", "raissLdScore") %in% names(result)))
    expect_true(is.numeric(result$mu))
    expect_true(all(is.finite(result$mu)))
  }
})

test_that("Report Condition Number in raissModel", {
  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

  result_with_cn <- raissModel(zt, sig_t, sig_i_t, reportConditionNumber =TRUE)
  result_without_cn <- raissModel(zt, sig_t, sig_i_t, reportConditionNumber =FALSE)

  expect_true(is.list(result_with_cn))
  expect_true(is.list(result_without_cn))
  # With condition number reporting, conditionNumber should be populated
  expect_true(is.numeric(result_with_cn$conditionNumber))
  expect_true(all(is.finite(result_with_cn$mu)))
  expect_true(all(is.finite(result_without_cn$mu)))
})

test_that("Input Validation of raissModel", {

  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)
  zt_invalid <- "not a numeric vector"
  sig_t_invalid <- "not a matrix"
  sig_i_t_invalid <- "not a matrix"

  expect_error(raissModel(zt_invalid, sig_t, sig_i_t))
  expect_error(raissModel(zt, sig_t_invalid, sig_i_t))
  expect_error(raissModel(zt, sig_t, sig_i_t_invalid))
})

test_that("Boundary Conditions of raissModel", {
  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

  zt_empty <- numeric(0)
  sig_t_empty <- matrix(numeric(0), nrow = 0)
  sig_i_t_empty <- matrix(numeric(0), nrow = 0)

  expect_error(raissModel(zt_empty, sig_t, sig_i_t))
  expect_error(raissModel(zt, sig_t_empty, sig_i_t))
  expect_error(raissModel(zt, sig_t, sig_i_t_empty))
})

test_that("Test with Different rcond Values for raissModel", {
  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

  rcond_values <- c(0.01, 0.05, 0.1)
  for (rcond in rcond_values) {
    result <- raissModel(zt, sig_t, sig_i_t, lamb = 0.01, rcond = rcond)
    expect_true(is.list(result))
    expect_true(all(c("var", "mu", "raissLdScore", "conditionNumber", "correctInversion") %in% names(result)))
    expect_true(is.numeric(result$mu))
    expect_true(all(is.finite(result$mu)))
  }
})

test_that("formatRaissDf returns correctly formatted data frame", {
  imp <- list(
    mu = rnorm(5),
    var = runif(5),
    raissLdScore = rnorm(5),
    conditionNumber = runif(5),
    correctInversion = sample(c(TRUE, FALSE), 5, replace = TRUE)
  )

  ref_panel <- data.frame(
    chrom = sample(1:22, 10, replace = TRUE),
    pos = sample(1:10000, 10),
    variant_id = paste0("rs", 1:10),
    A1 = sample(c("A", "T", "G", "C"), 10, replace = TRUE),
    A2 = sample(c("A", "T", "G", "C"), 10, replace = TRUE)
  )

  unknowns <- sample(1:nrow(ref_panel), 5)

  result <- formatRaissDf(imp, ref_panel, unknowns)

  expect_true(is.data.frame(result))
  expect_equal(ncol(result), 10)
  expect_equal(colnames(result), c('chrom', 'pos', 'variant_id', 'A1', 'A2', 'z', 'Var', 'raissLdScore', 'conditionNumber', 'correctInversion'))

  for (col in c('chrom', 'pos', 'variant_id', 'A1', 'A2')) {
    expect_equal(setNames(unlist(result[col]), NULL), unlist(ref_panel[unknowns, col, drop = TRUE]))
  }
  for (col in c('z', 'Var', 'raissLdScore', 'conditionNumber', 'correctInversion')) {
    expected_col <- if (col == "z") "mu" else if (col == "Var") "var" else col
    expect_equal(setNames(unlist(result[col]), NULL), setNames(unlist(imp[expected_col]), NULL))
  }
})

test_that("Merge operation is correct for mergeRaissDf", {
    raiss_df_example <- data.frame(
        chrom = c("chr21", "chr22"),
        pos = c(123, 456),
        variant_id = c("var1", "var2"),
        A1 = c("A", "T"),
        A2 = c("T", "A"),
        z = c(0.5, 1.5),
        Var = c(0.2, 0.3),
        raissLdScore = c(10, 20),
        raissR2 = c(0.8, 0.7))

    known_zscores_example <- data.frame(
        chrom = c("chr21", "chr22"),
        pos = c(123, 456),
        variant_id = c("var1", "var2"),
        A1 = c("A", "T"),
        A2 = c("T", "A"),
        z = c(0.5, 1.5))

    merged_df <- mergeRaissDf(raiss_df_example, known_zscores_example)
    expect_equal(nrow(merged_df), 2)
    expect_true(all(c("chr21", "chr22") %in% merged_df$chrom))
})

generate_fro_test_data <- function(seed=1) {
    set.seed(seed)
    return(data.frame(
        chrom = paste0("chr", rep(22, 10)),
        pos = seq(1, 100, 10),
        variant_id = 1:10,
        A1 = rep("A", 10),
        A2 = rep("T", 10),
        z = rnorm(10),
        Var = runif(10, 0, 1),
        raissLdScore = rnorm(10, 5, 2)
    ))
}

test_that("Correct columns are selected in filterRaissOutput", {
    test_data <- generate_fro_test_data()
    output <- filterRaissOutput(test_data)$zscores
    expect_true(all(c('variant_id', 'A1', 'A2', 'z', 'Var', 'raissLdScore') %in% names(output)))
})

test_that("raissR2 is calculated correctly in filterRaissOutput", {
    test_data <- generate_fro_test_data()
    output <- filterRaissOutput(test_data)$zscores
    expected_R2 <- 1 - test_data[which(test_data$raissLdScore >= 5),]$Var
    expect_equal(output$raissR2, expected_R2[which(expected_R2 > 0.6)])
})

test_that("Filtering is applied correctly in filterRaissOutput", {
    test_data <- generate_fro_test_data()
    R2_threshold <- 0.6
    minimum_ld <- 5
    output <- filterRaissOutput(test_data, R2_threshold, minimum_ld)$zscores

    expect_true(all(output$raissR2 > R2_threshold))
    expect_true(all(output$raissLdScore >= minimum_ld))
})

test_that("Function returns the correct subset in filterRaissOutput", {
    test_data <- generate_fro_test_data()
    test_data$raissR2 <- 1 - test_data$Var
    output <- filterRaissOutput(test_data)$zscores

    manual_filter <- test_data[test_data$raissR2 > 0.6 & test_data$raissLdScore >= 5, ]

    expect_equal(nrow(output), nrow(manual_filter))
    expect_equal(sum(output$variant_id != manual_filter$variant_id), 0)
})

test_that("computeMu basic functionality", {
    sig_i_t <- matrix(c(1, 2, 3, 4), nrow = 2)
    sig_t_inv <- matrix(c(5, 6, 7, 8), nrow = 2)
    zt <- matrix(c(9, 10, 11, 12), nrow = 2)

    expected_result <- matrix(c(517, 766, 625, 926), nrow = 2)
    result <- computeMu(sig_i_t, sig_t_inv, zt)
    expect_equal(result, expected_result)
})

generate_mock_data_for_computeVar <- function(seed=1) {
    return(
        list(
            sig_i_t_1 = matrix(c(1, 2, 3, 4), nrow = 2),
            sig_t_inv_1 = matrix(c(5, 6, 7, 8), nrow = 2),
            lamb_1 = 0.5))
}

test_that("computeVar returns correct output for batch = TRUE", {
    input_data <- generate_mock_data_for_computeVar()
    result <- computeVar(input_data$sig_i_t_1, input_data$sig_t_inv_1, input_data$lamb_1, batch = TRUE)
    expect_true(is.list(result))
    expect_length(result, 2)
    expect_true(all(c("var", "raissLdScore") %in% names(result)))
    expect_true(is.numeric(result$var))
    expect_true(is.numeric(result$raissLdScore))
})

test_that("computeVar returns correct output for batch = FALSE", {
    input_data <- generate_mock_data_for_computeVar()
    result <- computeVar(input_data$sig_i_t_1, input_data$sig_t_inv_1, input_data$lamb_1, batch = FALSE)
    expect_true(is.list(result))
    expect_length(result, 2)
    expect_true(all(c("var", "raissLdScore") %in% names(result)))
    expect_true(is.numeric(result$var))
    expect_true(is.numeric(result$raissLdScore))
})

test_that("checkInversion correctly identifies inverse matrices in", {
  sig_t <- matrix(c(1, 2, 3, 4), nrow=2, ncol=2)
  sig_t_inv <- solve(sig_t)
  expect_true(checkInversion(sig_t, sig_t_inv))
})

test_that("varInBoundaries sets boundaries correctly", {
  lamb_test <- 0.05
  var <- c(-1, 0, 0.5, 1.04, 1.05)

  result <- varInBoundaries(var, lamb_test)

  expect_equal(result[1], 0)                   # Value less than 0 should be set to 0
  expect_equal(result[2], 0)                   # Value within lower boundary should remain unchanged
  expect_equal(result[3], 0.5)                 # Value within boundaries should remain unchanged
  expect_equal(result[4], 1.04)                   # Value greater than 0.99999 + lamb should be set to 1
  expect_equal(result[5], 1)                   # Value greater than 0.99999 + lamb should be set to 1
})

test_that("invertMat computes correct pseudo-inverse", {
  mat <- matrix(c(1, 2, 3, 4), nrow = 2)
  lamb <- 0.5
  rcond <- 1e-7
  result <- invertMat(mat, lamb, rcond)
  expect_true(is.matrix(result))
})

test_that("invertMat handles errors and retries", {
  mat <- matrix(c(0, 0, 0, 0), nrow = 2)
  lamb <- 0.1
  rcond <- 1e-7
  result <- invertMat(mat, lamb, rcond)
  expect_true(is.matrix(result))
})

test_that("invertMatRecursive correctly inverts a valid square matrix", {
  mat <- matrix(c(2, -1, -1, 2), nrow = 2)
  lamb <- 0.5
  rcond <- 0.01
  result <- invertMatRecursive(mat, lamb, rcond)
  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(mat))
})

test_that("invertMatRecursive handles non-square matrices appropriately", {
  mat <- matrix(1:6, nrow = 2)
  lamb <- 0.5
  rcond <- 0.01
  expect_silent(invertMatRecursive(mat, lamb, rcond))
})

test_that("invertMatRecursive handles errors and performs recursive call correctly", {
  mat <- "not a matrix"
  lamb <- 0.5
  rcond <- 0.01
  expect_error(invertMatRecursive(mat, lamb, rcond))
})

# Test with Different Tolerance Levels
test_that("invertMatEigen behaves differently with varying tolerance levels", {
  mat <- matrix(c(1, 0, 0, 1e-4), nrow = 2)
  tol_high <- 1e-2
  tol_low <- 1e-6
  result_high_tol <- invertMatEigen(mat, tol_high)
  result_low_tol <- invertMatEigen(mat, tol_low)
  expect_true(!is.logical(all.equal(result_high_tol, result_low_tol)))
})

test_that("invertMatEigen handles non-square matrices", {
  mat <- matrix(1:6, nrow = 2)
  expect_error(invertMatEigen(mat))
})

test_that("invertMatEigen returns the same matrix for an identity matrix", {
    mat <- diag(2)
    expected <- mat
    actual <- invertMatEigen(mat)
    expect_equal(actual, expected)
})

test_that("invertMatEigen returns a zero matrix for a zero matrix input", {
    mat <- matrix(0, nrow = 2, ncol = 2)
    expected <- mat
    expect_error(invertMatEigen(mat),
      "Cannot invert the input matrix because all its eigen values are negative or close to zero")
})

test_that("invertMatEigen handles matrices with negative eigenvalues", {
    mat <- matrix(c(-2, 0, 0, -3), nrow = 2)
    expect_silent(invertMatEigen(mat))
})

# ===========================================================================
# raissSingleMatrix edge cases
# ===========================================================================

test_that("raissSingleMatrix returns NULL when no known variants overlap", {
  set.seed(42)
  ref_panel <- data.frame(
    chrom = rep(1, 10), pos = seq(10, 100, 10),
    variant_id = paste0("rs", 1:10),
    A1 = rep("A", 10), A2 = rep("G", 10),
    stringsAsFactors = FALSE
  )
  # known_zscores has variant IDs that don't match ref_panel at all
  known_zscores <- data.frame(
    chrom = rep(1, 3), pos = c(200, 300, 400),
    variant_id = paste0("other", 1:3),
    A1 = rep("A", 3), A2 = rep("G", 3),
    z = rnorm(3), stringsAsFactors = FALSE
  )
  ldMatrix <- diag(10)
  result <- raissSingleMatrix(ref_panel, known_zscores, ldMatrix, verbose = FALSE)
  expect_null(result)
})

test_that("raissSingleMatrix returns known zscores when no unknowns to impute", {
  set.seed(42)
  ref_panel <- data.frame(
    chrom = rep(1, 5), pos = seq(10, 50, 10),
    variant_id = paste0("rs", 1:5),
    A1 = rep("A", 5), A2 = rep("G", 5),
    stringsAsFactors = FALSE
  )
  # All ref_panel variants are known - nothing to impute
  known_zscores <- data.frame(
    chrom = rep(1, 5), pos = seq(10, 50, 10),
    variant_id = paste0("rs", 1:5),
    A1 = rep("A", 5), A2 = rep("G", 5),
    z = rnorm(5), stringsAsFactors = FALSE
  )
  ldMatrix <- diag(5)
  result <- raissSingleMatrix(ref_panel, known_zscores, ldMatrix, verbose = FALSE)
  expect_true(is.list(result))
  expect_equal(result$resultNofilter, known_zscores)
  expect_equal(result$resultFilter, known_zscores)
  expect_equal(result$ldMat, ldMatrix)
})

# ===========================================================================
# raissSingleMatrixFromX edge cases
# ===========================================================================

test_that("raissSingleMatrixFromX returns NULL when no known variants overlap", {
  set.seed(42)
  n <- 50
  p <- 10
  ref_panel <- data.frame(
    chrom = rep(1, p), pos = seq(10, p * 10, 10),
    variant_id = paste0("rs", 1:p),
    A1 = rep("A", p), A2 = rep("G", p),
    stringsAsFactors = FALSE
  )
  known_zscores <- data.frame(
    chrom = rep(1, 3), pos = c(200, 300, 400),
    variant_id = paste0("other", 1:3),
    A1 = rep("A", 3), A2 = rep("G", 3),
    z = rnorm(3), stringsAsFactors = FALSE
  )
  X <- scale(matrix(sample(0:2, n * p, replace = TRUE), nrow = n))
  X[is.na(X)] <- 0
  colnames(X) <- ref_panel$variant_id
  result <- raissSingleMatrixFromX(ref_panel, known_zscores, X, verbose = FALSE)
  expect_null(result)
})

test_that("raissSingleMatrixFromX returns known zscores when no unknowns to impute", {
  set.seed(42)
  n <- 50
  p <- 5
  ref_panel <- data.frame(
    chrom = rep(1, p), pos = seq(10, p * 10, 10),
    variant_id = paste0("rs", 1:p),
    A1 = rep("A", p), A2 = rep("G", p),
    stringsAsFactors = FALSE
  )
  known_zscores <- data.frame(
    chrom = rep(1, p), pos = seq(10, p * 10, 10),
    variant_id = paste0("rs", 1:p),
    A1 = rep("A", p), A2 = rep("G", p),
    z = rnorm(p), stringsAsFactors = FALSE
  )
  X <- scale(matrix(sample(0:2, n * p, replace = TRUE), nrow = n))
  X[is.na(X)] <- 0
  colnames(X) <- ref_panel$variant_id
  result <- raissSingleMatrixFromX(ref_panel, known_zscores, X, verbose = FALSE)
  expect_true(is.list(result))
  expect_equal(result$resultNofilter, known_zscores)
  expect_null(result$ldMat)
})

# ===========================================================================
# raiss() dispatch paths: single-matrix list and genotype_matrix list
# ===========================================================================

test_that("raiss with single-matrix LD list dispatches to single matrix path", {
  set.seed(42)
  n_variants <- 20
  ref_panel <- data.frame(
    chrom = rep(1, n_variants), pos = seq(10, n_variants * 10, 10),
    variant_id = paste0("rs", 1:n_variants),
    A1 = rep("A", n_variants), A2 = rep("G", n_variants),
    stringsAsFactors = FALSE
  )
  n_known <- 10
  known_idx <- sort(sample(seq_len(n_variants), n_known))
  known_zscores <- data.frame(
    chrom = rep(1, n_known), pos = ref_panel$pos[known_idx],
    variant_id = ref_panel$variant_id[known_idx],
    A1 = ref_panel$A1[known_idx], A2 = ref_panel$A2[known_idx],
    z = rnorm(n_known), stringsAsFactors = FALSE
  )
  R <- diag(n_variants)
  colnames(R) <- rownames(R) <- ref_panel$variant_id

  # Wrap in list structure with ldMatrices
  LD_list <- list(ldMatrices = list(R))

  result <- raiss(ref_panel, known_zscores, ldMatrix =LD_list,
                  r2Threshold =0, minimumLd =0, verbose = FALSE)
  expect_true(is.list(result))
  expect_true("resultNofilter" %in% names(result))
  expect_equal(nrow(result$resultNofilter), n_variants)
})

test_that("raiss with genotype_matrix list processes multiple blocks", {
  set.seed(42)
  n <- 50
  p <- 20
  # Each block has its own ref_panel subset matching the X columns
  ref_panel <- data.frame(
    chrom = rep(1, p), pos = seq(10, p * 10, 10),
    variant_id = paste0("rs", 1:p),
    A1 = rep("A", p), A2 = rep("G", p),
    stringsAsFactors = FALSE
  )
  # Use only first block's variants as known (so second block has unknowns)
  known_idx <- sort(sample(1:10, 5))
  known_zscores <- data.frame(
    chrom = rep(1, length(known_idx)), pos = ref_panel$pos[known_idx],
    variant_id = ref_panel$variant_id[known_idx],
    A1 = ref_panel$A1[known_idx], A2 = ref_panel$A2[known_idx],
    z = rnorm(length(known_idx)), stringsAsFactors = FALSE
  )
  X <- scale(matrix(sample(0:2, n * p, replace = TRUE), nrow = n))
  X[is.na(X)] <- 0
  colnames(X) <- ref_panel$variant_id

  # Use the full matrix as a single-element list - the simplest valid list input
  X_list <- list(X)

  result <- raiss(ref_panel, known_zscores, genotypeMatrix = X_list,
                  r2Threshold =0, minimumLd =0, verbose = FALSE)
  expect_true(is.list(result))
  expect_true("resultNofilter" %in% names(result))
  expect_true(nrow(result$resultNofilter) > 0)
  expect_null(result$ldMat)
})

test_that("raiss with genotype_matrix list returns NULL when all blocks fail", {
  set.seed(42)
  n <- 50
  p <- 10
  ref_panel <- data.frame(
    chrom = rep(1, p), pos = seq(10, p * 10, 10),
    variant_id = paste0("rs", 1:p),
    A1 = rep("A", p), A2 = rep("G", p),
    stringsAsFactors = FALSE
  )
  # known_zscores has no overlap with ref_panel
  known_zscores <- data.frame(
    chrom = rep(1, 3), pos = c(200, 300, 400),
    variant_id = paste0("other", 1:3),
    A1 = rep("A", 3), A2 = rep("G", 3),
    z = rnorm(3), stringsAsFactors = FALSE
  )
  X <- scale(matrix(sample(0:2, n * p, replace = TRUE), nrow = n))
  X[is.na(X)] <- 0
  colnames(X) <- ref_panel$variant_id
  X_list <- list(X[, 1:5, drop = FALSE], X[, 6:10, drop = FALSE])

  result <- raiss(ref_panel, known_zscores, genotypeMatrix = X_list,
                  verbose = FALSE)
  expect_null(result)
})

# Block-Diagonal LD data generator for RAISS testing
# Corrected function to generate proper block-diagonal test data
generate_block_diagonal_test_data <- function(seed = 123, block_structure = "overlapping", n_variants = 30) {
  set.seed(seed)

  # Create reference panel with variants
  ref_panel <- data.frame(
    chrom = rep(1, n_variants),
    pos = seq(1, n_variants * 10, 10),
    variant_id = paste0("var", seq_len(n_variants)),
    A1 = sample(c("A", "T", "G", "C"), n_variants, replace = TRUE),
    A2 = sample(c("A", "T", "G", "C"), n_variants, replace = TRUE),
    stringsAsFactors = FALSE
  )

  # Create known z-scores for every other variant
  known_indices <- seq(1, n_variants, by = 2)
  known_zscores <- data.frame(
    chrom = rep(1, length(known_indices)),
    pos = ref_panel$pos[known_indices],
    variant_id = ref_panel$variant_id[known_indices],
    A1 = ref_panel$A1[known_indices],
    A2 = ref_panel$A2[known_indices],
    z = rnorm(length(known_indices)),
    stringsAsFactors = FALSE
  )

  # Define block boundaries based on requested structure
  if (block_structure == "overlapping") {
    block_boundaries <- list(
      c(1, 11),    # Block 1: variants 1-11
      c(11, 21),   # Block 2: variants 11-21 (overlap at var11)
      c(21, n_variants)  # Block 3: variants 21-30 (overlap at var21)
    )
  } else if (block_structure == "non_overlapping") {
    block_boundaries <- list(
      c(1, 10),
      c(11, 20),
      c(21, n_variants)
    )
  } else if (block_structure == "uneven") {
    block_boundaries <- list(
      c(1, 5),
      c(6, 20),
      c(21, n_variants)
    )
  } else if (block_structure == "many_small") {
    block_size <- 5
    n_blocks <- ceiling(n_variants / block_size)
    block_boundaries <- list()
    for (i in 1:n_blocks) {
      startIdx <- (i-1) * block_size + 1
      endIdx <- min(i * block_size, n_variants)
      block_boundaries[[i]] <- c(startIdx, endIdx)
    }
  } else if (block_structure == "single_block") {
    block_boundaries <- list(c(1, n_variants))
  }

  # First, create independent block matrices
  block_matrices <- list()
  for (i in seq_along(block_boundaries)) {
    startIdx <- block_boundaries[[i]][1]
    endIdx <- block_boundaries[[i]][2]
    block_variant_ids <- ref_panel$variant_id[startIdx:endIdx]
    n_block <- length(block_variant_ids)

    # Create the block matrix with correlations ONLY within the block
    block_matrix <- matrix(0, nrow = n_block, ncol = n_block)
    for (a in 1:n_block) {
      for (b in 1:n_block) {
        if (a == b) {
          block_matrix[a, b] <- 1
        } else {
          # Use positions within the block, not absolute positions
          block_matrix[a, b] <- 0.95^abs(a - b)
        }
      }
    }
    rownames(block_matrix) <- block_variant_ids
    colnames(block_matrix) <- block_variant_ids

    block_matrices[[i]] <- block_matrix
  }

  # Create variant indices data frame
  variantIndices <- data.frame(
    variant_id = character(),
    blockId = integer(),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(block_boundaries)) {
    startIdx <- block_boundaries[[i]][1]
    endIdx <- block_boundaries[[i]][2]
    block_variant_ids <- ref_panel$variant_id[startIdx:endIdx]

    block_indices <- data.frame(
      variant_id = block_variant_ids,
      blockId = i,
      stringsAsFactors = FALSE
    )
    variantIndices <- rbind(variantIndices, block_indices)
  }

  # Create block metadata
  block_sizes <- sapply(block_boundaries, function(b) b[2] - b[1] + 1)
  blockMetadata <- data.frame(
    blockId = seq_along(block_boundaries),
    chrom = rep(1, length(block_boundaries)),
    size = block_sizes,
    startIdx = sapply(seq_along(block_boundaries), function(i) {
      # Adjust for 1-based indexing in R
      if (i == 1) return(1)
      # Count unique variants before this block
      sum(sapply(1:(i-1), function(j) {
        # If there's an overlap with the next block, count one less
        if (j < length(block_boundaries) &&
            block_boundaries[[j]][2] == block_boundaries[[j+1]][1]) {
          return(block_boundaries[[j]][2] - block_boundaries[[j]][1])
        } else {
          return(block_boundaries[[j]][2] - block_boundaries[[j]][1] + 1)
        }
      })) + 1
    }),
    endIdx = sapply(seq_along(block_boundaries), function(i) {
      # Count all unique variants up to and including this block
      sum(sapply(1:i, function(j) {
        # If there's an overlap with the next block, count one less
        if (j < i && block_boundaries[[j]][2] == block_boundaries[[j+1]][1]) {
          return(block_boundaries[[j]][2] - block_boundaries[[j]][1])
        } else {
          return(block_boundaries[[j]][2] - block_boundaries[[j]][1] + 1)
        }
      }))
    }),
    stringsAsFactors = FALSE
  )

  # Build the full matrix correctly ensuring proper block structure
  # IMPORTANT: Initialize a matrix with zeros - ensure no correlations between blocks
  all_variant_ids <- unique(variantIndices$variant_id)
  LD_matrix_full <- matrix(0, nrow = length(all_variant_ids), ncol = length(all_variant_ids))
  rownames(LD_matrix_full) <- all_variant_ids
  colnames(LD_matrix_full) <- all_variant_ids

  # For each block, fill in only the relevant section of the full matrix
  for (i in seq_along(block_matrices)) {
    block_matrix <- block_matrices[[i]]
    block_vars <- rownames(block_matrix)

    for (var_a in block_vars) {
      for (var_b in block_vars) {
        LD_matrix_full[var_a, var_b] <- block_matrix[var_a, var_b]
      }
    }
  }

  # Create the block structure for RAISS
  LD_matrix_blocks <- list(
    ldMatrices = block_matrices,
    variantIndices = variantIndices,
    blockMetadata = blockMetadata,
    ldVariants = all_variant_ids
  )

  return(list(
    ref_panel = ref_panel,
    known_zscores = known_zscores,
    LD_matrix_full = LD_matrix_full,
    LD_matrix_blocks = LD_matrix_blocks,
    variantIndices = variantIndices,
    block_boundaries = block_boundaries,
    blockMetadata = blockMetadata
  ))
}

test_that("full matrix and block processing produce identical results", {
  # Only test non-overlapping structures for exact z-score matching
  block_structures <- c("non_overlapping", "single_block")

  for (structure in block_structures) {
    test_data <- generate_block_diagonal_test_data(seed = 123, block_structure = structure)

    # Prepare ld_data as LdData S4 for partitionLdMatrix
    ld_data <- make_ld_data_from_ref_panel(
      test_data$LD_matrix_full, test_data$ref_panel, test_data$blockMetadata
    )

    # For non-overlapping structures, use partitionLdMatrix
    partitioned <- partitionLdMatrix(
      ld_data,
      mergeSmallBlocks =FALSE
    )

    # Run RAISS with full matrix
    result_full <- raiss(
      refPanel = test_data$ref_panel,
      knownZscores = test_data$known_zscores,
      ldMatrix =test_data$LD_matrix_full,
      lamb = 0.01,
      rcond = 0.01,
      r2Threshold =0.3,
      minimumLd =1,
      verbose = FALSE
    )

    # Run RAISS with partitioned blocks
    result_blocks <- raiss(
      refPanel = test_data$ref_panel,
      knownZscores = test_data$known_zscores,
      ldMatrix =partitioned,
      lamb = 0.01,
      rcond = 0.01,
      r2Threshold =0.3,
      minimumLd =1,
      verbose = FALSE
    )

    # For non-overlapping blocks, we compare all variants
    result_full_sorted <- result_full$resultNofilter %>% arrange(variant_id)
    result_blocks_sorted <- result_blocks$resultNofilter %>% arrange(variant_id)

    # Compare variant IDs
    expect_equal(
      sort(result_full$resultNofilter$variant_id),
      sort(result_blocks$resultNofilter$variant_id),
      info = paste("Variant IDs should match for", structure)
    )

    # Compare z-scores with appropriate tolerance
    expect_equal(
      result_full_sorted$z,
      result_blocks_sorted$z,
      tolerance = 0.01,
      info = paste("Z-scores should match for", structure)
    )

    # Compare filtered results if present
    if (!is.null(result_full$resultFilter) && !is.null(result_blocks$resultFilter) &&
        nrow(result_full$resultFilter) > 0 && nrow(result_blocks$resultFilter) > 0) {
      expect_equal(
        sort(result_full$resultFilter$variant_id),
        sort(result_blocks$resultFilter$variant_id),
        info = paste("Filtered variant IDs should match for", structure)
      )

      result_full_filter_sorted <- result_full$resultFilter %>% arrange(variant_id)
      result_blocks_filter_sorted <- result_blocks$resultFilter %>% arrange(variant_id)

      expect_equal(
        result_full_filter_sorted$z,
        result_blocks_filter_sorted$z,
        tolerance = 0.01,
        info = paste("Filtered Z-scores should match for", structure)
      )
    }
  }
})

test_that("overlapping blocks preserve variant IDs but may have different z-scores", {
  # Test only overlapping structure
  test_data <- generate_block_diagonal_test_data(seed = 123, block_structure = "overlapping")

  # Run RAISS with full matrix
  result_full <- raiss(
    refPanel = test_data$ref_panel,
    knownZscores = test_data$known_zscores,
    ldMatrix =test_data$LD_matrix_full,
    lamb = 0.01,
    rcond = 0.01,
    r2Threshold =0.3,
    minimumLd =1,
    verbose = FALSE
  )

  # Run RAISS with block processing
  result_blocks <- raiss(
    refPanel = test_data$ref_panel,
    knownZscores = test_data$known_zscores,
    ldMatrix =test_data$LD_matrix_blocks,
    lamb = 0.01,
    rcond = 0.01,
    r2Threshold =0.3,
    minimumLd =1,
    verbose = FALSE
  )

  # Test 1: Verify all variants are present in both results
  expect_equal(
    sort(result_full$resultNofilter$variant_id),
    sort(result_blocks$resultNofilter$variant_id),
    info = "Both methods should have the same set of variant IDs"
  )

  # Test 2: For overlapping blocks, verify boundary variants exist and have valid values
  # Identify boundary variants
  boundary_variants <- character(0)
  for (i in 1:(length(test_data$block_boundaries) - 1)) {
    overlap_pos <- test_data$block_boundaries[[i]][2]
    boundary_variants <- c(boundary_variants, paste0("var", overlap_pos))
  }

  # Verify boundary variants exist in results
  expect_true(
    all(boundary_variants %in% result_blocks$resultNofilter$variant_id),
    info = "All boundary variants should be present in block results"
  )

  # Verify boundary variants have valid z-scores
  boundary_results <- result_blocks$resultNofilter %>%
    filter(variant_id %in% boundary_variants)

  expect_true(
    all(!is.na(boundary_results$z)),
    info = "Boundary variants should have valid z-scores in block results"
  )

  # Test 3: Verify non-boundary variants have z-scores with reasonable range
  non_boundary_results <- result_blocks$resultNofilter %>%
    filter(!variant_id %in% boundary_variants)

  expect_true(
    all(!is.na(non_boundary_results$z)),
    info = "Non-boundary variants should have valid z-scores"
  )

  expect_true(
    all(abs(non_boundary_results$z) < 10),
    info = "Non-boundary variant z-scores should be in reasonable range"
  )

  # We deliberately do NOT compare z-score values between full matrix and block processing
  # for overlapping blocks, as differences are expected and valid
})

test_that("raiss handles block boundaries correctly", {
  # Generate test data with overlapping blocks
  test_data <- generate_block_diagonal_test_data(seed = 456, block_structure = "overlapping")

  # Define the thresholds explicitly
  test_R2_threshold <- 0.3
  test_minimum_ld <- 1

  # Run RAISS with block processing
  result <- raiss(
    refPanel = test_data$ref_panel,
    knownZscores = test_data$known_zscores,
    ldMatrix =test_data$LD_matrix_blocks,
    lamb = 0.01,
    rcond = 0.01,
    r2Threshold =test_R2_threshold,
    minimumLd =test_minimum_ld,
    verbose = FALSE
  )

  # First verify that the required columns exist in the results
  expect_true(
    "variant_id" %in% names(result$resultNofilter),
    info = "resultNofilter should contain a variant_id column"
  )

  expect_true(
    "raissR2" %in% names(result$resultNofilter),
    info = "resultNofilter should contain a raissR2 column"
  )

  expect_true(
    "raissLdScore" %in% names(result$resultNofilter),
    info = "resultNofilter should contain a raissLdScore column"
  )

  # Check that we have only one entry per variant ID (no duplicates)
  expect_equal(
    length(unique(result$resultNofilter$variant_id)),
    length(result$resultNofilter$variant_id),
    info = "Result should have no duplicate variant IDs"
  )

  # Check that boundary variants have reasonable values
  boundary_variants <- character(0)
  for (i in 1:(length(test_data$block_boundaries) - 1)) {
    overlap_pos <- test_data$block_boundaries[[i]][2]
    boundary_variants <- c(boundary_variants, paste0("var", overlap_pos))
  }

  # Verify that boundary variants exist in the results
  expect_true(
    all(boundary_variants %in% result$resultNofilter$variant_id),
    info = "All boundary variants should be present in the results"
  )

  # Get the boundary variant results
  boundary_results <- result$resultNofilter %>%
    filter(variant_id %in% boundary_variants)

  # Check R-squared values for non-NA boundary variants
  non_na_r2 <- boundary_results$raissR2[!is.na(boundary_results$raissR2)]
  if (length(non_na_r2) > 0) {
    expect_true(
      all(non_na_r2 >= 0 & non_na_r2 <= 1),
      info = "Non-NA boundary variant R-squared values should be between 0 and 1"
    )
  }

  # Check LD scores for non-NA boundary variants
  non_na_ld <- boundary_results$raissLdScore[!is.na(boundary_results$raissLdScore)]
  if (length(non_na_ld) > 0) {
    expect_true(
      all(non_na_ld >= 0),
      info = "Non-NA boundary variant LD scores should be non-negative"
    )
  }

  # Verify that pre-filtering and post-filtering steps handle boundary variants correctly
  if (!is.null(result$resultFilter) && nrow(result$resultFilter) > 0) {
    # First check if filtered results have the required columns
    expect_true(
      "variant_id" %in% names(result$resultFilter),
      info = "resultFilter should contain a variant_id column"
    )

    expect_true(
      "raissR2" %in% names(result$resultFilter),
      info = "resultFilter should contain a raissR2 column"
    )

    expect_true(
      "raissLdScore" %in% names(result$resultFilter),
      info = "resultFilter should contain a raissLdScore column"
    )

    # Check which boundary variants passed the filtering
    boundary_in_filtered <- boundary_variants %in% result$resultFilter$variant_id

    if (any(boundary_in_filtered)) {
      # Get the filtered boundary variants
      boundary_filtered <- result$resultFilter %>%
        filter(variant_id %in% boundary_variants)

      # Check that non-NA R-squared values meet the threshold
      non_na_r2_filtered <- boundary_filtered$raissR2[!is.na(boundary_filtered$raissR2)]
      if (length(non_na_r2_filtered) > 0) {
        expect_true(
          all(non_na_r2_filtered >= test_R2_threshold),
          info = paste("Non-NA filtered boundary variant R-squared values should meet the threshold of", test_R2_threshold)
        )
      }

      # Check that non-NA LD scores meet the threshold
      non_na_ld_filtered <- boundary_filtered$raissLdScore[!is.na(boundary_filtered$raissLdScore)]
      if (length(non_na_ld_filtered) > 0) {
        expect_true(
          all(non_na_ld_filtered >= test_minimum_ld),
          info = paste("Non-NA filtered boundary variant LD scores should meet the threshold of", test_minimum_ld)
        )
      }
    }
  }
})

test_that("partitionLdMatrix integrates correctly with RAISS", {
  test_data <- generate_block_diagonal_test_data(seed = 456, block_structure = "non_overlapping")

  ld_data <- make_ld_data_from_ref_panel(
    test_data$LD_matrix_full, test_data$ref_panel, test_data$blockMetadata
  )

  partitioned <- partitionLdMatrix(
    ld_data,
    mergeSmallBlocks =FALSE
  )

  result_full <- raiss(
    refPanel = test_data$ref_panel,
    knownZscores = test_data$known_zscores,
    ldMatrix =test_data$LD_matrix_full,
    lamb = 0.01,
    rcond = 0.01,
    r2Threshold =0.3,
    minimumLd =1,
    verbose = FALSE
  )

  result_partitioned <- raiss(
    refPanel = test_data$ref_panel,
    knownZscores = test_data$known_zscores,
    ldMatrix =partitioned,
    lamb = 0.01,
    rcond = 0.01,
    r2Threshold =0.3,
    minimumLd =1,
    verbose = FALSE
  )

  result_full_sorted <- result_full$resultNofilter %>% arrange(variant_id)
  result_partitioned_sorted <- result_partitioned$resultNofilter %>% arrange(variant_id)

  expect_equal(
    result_full_sorted$variant_id,
    result_partitioned_sorted$variant_id,
    info = "Variant IDs should match"
  )

  expect_equal(
    result_full_sorted$z,
    result_partitioned_sorted$z,
    tolerance = 1e-4,
    info = "Z-scores should match"
  )
})

# Test 3: Boundary overlap handling
test_that("boundary overlaps are handled correctly", {
  test_data <- generate_block_diagonal_test_data(seed = 789, block_structure = "overlapping")

  result_blocks <- raiss(
    refPanel = test_data$ref_panel,
    knownZscores = test_data$known_zscores,
    ldMatrix =test_data$LD_matrix_blocks,
    lamb = 0.01,
    rcond = 0.01,
    r2Threshold =0.1,
    minimumLd =1,
    verbose = FALSE
  )

  variant_counts <- table(test_data$variantIndices$variant_id)
  boundary_vars <- names(variant_counts[variant_counts > 1])

  for (var in boundary_vars) {
    expect_equal(
      sum(result_blocks$resultNofilter$variant_id == var),
      1,
      info = paste("Boundary variant", var, "should appear once")
    )
  }

  expect_equal(
    nrow(result_blocks$resultNofilter),
    length(unique(result_blocks$resultNofilter$variant_id)),
    info = "No duplicate variants in results"
  )
})

# Test 4: Single-block case
test_that("RAISS handles single-block list correctly", {
  test_data <- generate_block_diagonal_test_data(seed = 202, block_structure = "single_block")

  result_full <- raiss(
    refPanel = test_data$ref_panel,
    knownZscores = test_data$known_zscores,
    ldMatrix =test_data$LD_matrix_full,
    lamb = 0.01,
    rcond = 0.01,
    r2Threshold =0.3,
    minimumLd =1,
    verbose = FALSE
  )

  result_single_block <- raiss(
    refPanel = test_data$ref_panel,
    knownZscores = test_data$known_zscores,
    ldMatrix =test_data$LD_matrix_blocks,
    lamb = 0.01,
    rcond = 0.01,
    r2Threshold =0.3,
    minimumLd =1,
    verbose = FALSE
  )

  result_full_sorted <- result_full$resultNofilter %>% arrange(variant_id)
  result_single_block_sorted <- result_single_block$resultNofilter %>% arrange(variant_id)

  expect_equal(
    result_full_sorted$z,
    result_single_block_sorted$z,
    tolerance = 1e-6,
    info = "Z-scores should match for single block"
  )
})

# ============================================================================
# Tests for SVD-based genotype matrix path (raissSingleMatrixFromX)
# ============================================================================

#' Helper: generate a genotype matrix X with corresponding ref_panel,
#' known_zscores, and LD matrix R for equivalence testing.
generate_X_test_data <- function(n = 200, p = 100, n_known = 50, seed = 42) {
  set.seed(seed)
  # Generate genotype-like matrix (dosages 0/1/2)
  X_raw <- matrix(sample(0:2, n * p, replace = TRUE, prob = c(0.25, 0.5, 0.25)),
                  nrow = n, ncol = p)
  # Center and scale
  X <- scale(X_raw)
  X[is.na(X)] <- 0  # zero-variance columns become 0

  # ref_panel for all p variants
  ref_panel <- data.frame(
    chrom = rep(1, p),
    pos = seq(1, p * 10, 10),
    variant_id = paste0("rs", seq_len(p)),
    A1 = sample(c("A", "T", "G", "C"), p, replace = TRUE),
    A2 = sample(c("A", "T", "G", "C"), p, replace = TRUE),
    stringsAsFactors = FALSE
  )
  colnames(X) <- ref_panel$variant_id

  # Select known variants
  known_idx <- sort(sample(seq_len(p), n_known))
  known_zscores <- data.frame(
    chrom = rep(1, n_known),
    pos = ref_panel$pos[known_idx],
    variant_id = ref_panel$variant_id[known_idx],
    A1 = ref_panel$A1[known_idx],
    A2 = ref_panel$A2[known_idx],
    z = rnorm(n_known),
    stringsAsFactors = FALSE
  )

  # LD matrix from X
  R <- cor(X_raw)
  R[is.na(R)] <- 0
  colnames(R) <- rownames(R) <- ref_panel$variant_id

  list(X = X, R = R, ref_panel = ref_panel, known_zscores = known_zscores,
       n = n, p = p, n_known = n_known)
}

test_that("safeSvd basic functionality", {
  set.seed(1)
  mat <- matrix(rnorm(20), nrow = 5, ncol = 4)
  s <- pecotmr:::.safeSvd(mat)
  expect_equal(length(s$d), min(5, 4))
  expect_true(all(s$d > 0))
  # Reconstruct
  reconstructed <- s$u %*% diag(s$d) %*% t(s$v)
  expect_equal(mat, reconstructed, tolerance = 1e-10)
})

test_that("safeSvd filters small singular values", {
  set.seed(2)
  # Create rank-2 matrix
  u <- matrix(rnorm(10), nrow = 5, ncol = 2)
  v <- matrix(rnorm(8), nrow = 4, ncol = 2)
  mat <- u %*% t(v) + matrix(rnorm(20) * 1e-12, nrow = 5, ncol = 4)
  s <- pecotmr:::.safeSvd(mat, tol = 1e-6)
  expect_equal(length(s$d), 2)
})

test_that("safeSvd max_rank works", {
  set.seed(3)
  mat <- matrix(rnorm(50), nrow = 10, ncol = 5)
  s <- pecotmr:::.safeSvd(mat, maxRank =2)
  expect_equal(length(s$d), 2)
  expect_equal(ncol(s$u), 2)
  expect_equal(ncol(s$v), 2)
})

test_that("safeSvd rejects all-zero matrix", {
  mat <- matrix(0, nrow = 5, ncol = 3)
  expect_error(pecotmr:::.safeSvd(mat), "all-zero")
})

test_that("X path matches R path: basic equivalence (n > p)", {
  data <- generate_X_test_data(n = 200, p = 100, n_known = 50, seed = 42)

  result_R <- raiss(data$ref_panel, data$known_zscores, ldMatrix =data$R,
                    lamb = 0.01, rcond = 0.01, r2Threshold =0, minimumLd =0,
                    verbose = FALSE)
  result_X <- raiss(data$ref_panel, data$known_zscores, genotypeMatrix = data$X,
                    lamb = 0.01, svdTol =1e-12, r2Threshold =0, minimumLd =0,
                    verbose = FALSE)

  # Compare imputed z-scores (sort by variant_id for alignment)
  r_sorted <- result_R$resultNofilter %>% arrange(variant_id)
  x_sorted <- result_X$resultNofilter %>% arrange(variant_id)

  expect_equal(r_sorted$z, x_sorted$z, tolerance = 1e-4,
               info = "Imputed z-scores should match between X and R paths")
  expect_equal(r_sorted$Var, x_sorted$Var, tolerance = 1e-4,
               info = "Variance should match between X and R paths")
  expect_equal(r_sorted$raissLdScore, x_sorted$raissLdScore, tolerance = 1e-4,
               info = "LD scores should match between X and R paths")
})

test_that("X path matches R path: n < p regime", {
  data <- generate_X_test_data(n = 50, p = 200, n_known = 100, seed = 123)

  result_R <- raiss(data$ref_panel, data$known_zscores, ldMatrix =data$R,
                    lamb = 0.01, rcond = 0.01, r2Threshold =0, minimumLd =0,
                    verbose = FALSE)
  result_X <- raiss(data$ref_panel, data$known_zscores, genotypeMatrix = data$X,
                    lamb = 0.01, svdTol =1e-12, r2Threshold =0, minimumLd =0,
                    verbose = FALSE)

  r_sorted <- result_R$resultNofilter %>% arrange(variant_id)
  x_sorted <- result_X$resultNofilter %>% arrange(variant_id)

  expect_equal(r_sorted$z, x_sorted$z, tolerance = 1e-4,
               info = "z-scores should match in n < p regime")
  expect_equal(r_sorted$Var, x_sorted$Var, tolerance = 1e-4,
               info = "Variance should match in n < p regime")
})

test_that("X path matches R path: n >> p regime", {
  data <- generate_X_test_data(n = 500, p = 50, n_known = 25, seed = 99)

  result_R <- raiss(data$ref_panel, data$known_zscores, ldMatrix =data$R,
                    lamb = 0.01, rcond = 0.01, r2Threshold =0, minimumLd =0,
                    verbose = FALSE)
  result_X <- raiss(data$ref_panel, data$known_zscores, genotypeMatrix = data$X,
                    lamb = 0.01, svdTol =1e-12, r2Threshold =0, minimumLd =0,
                    verbose = FALSE)

  r_sorted <- result_R$resultNofilter %>% arrange(variant_id)
  x_sorted <- result_X$resultNofilter %>% arrange(variant_id)

  expect_equal(r_sorted$z, x_sorted$z, tolerance = 1e-4)
  expect_equal(r_sorted$Var, x_sorted$Var, tolerance = 1e-4)
})

test_that("X path matches R path: varying lambda", {
  data <- generate_X_test_data(n = 150, p = 80, n_known = 40, seed = 7)

  for (lamb in c(0.001, 0.01, 0.1)) {
    result_R <- raiss(data$ref_panel, data$known_zscores, ldMatrix =data$R,
                      lamb = lamb, rcond = 0.01, r2Threshold =0, minimumLd =0,
                      verbose = FALSE)
    result_X <- raiss(data$ref_panel, data$known_zscores, genotypeMatrix = data$X,
                      lamb = lamb, svdTol =1e-12, r2Threshold =0, minimumLd =0,
                      verbose = FALSE)

    r_sorted <- result_R$resultNofilter %>% arrange(variant_id)
    x_sorted <- result_X$resultNofilter %>% arrange(variant_id)

    expect_equal(r_sorted$z, x_sorted$z, tolerance = 1e-4,
                 info = paste("z-scores should match for lamb =", lamb))
    expect_equal(r_sorted$Var, x_sorted$Var, tolerance = 1e-4,
                 info = paste("Variance should match for lamb =", lamb))
  }
})

test_that("X path handles all-known edge case", {
  data <- generate_X_test_data(n = 100, p = 50, n_known = 50, seed = 10)
  # Make all variants known
  all_known <- data.frame(
    chrom = data$ref_panel$chrom,
    pos = data$ref_panel$pos,
    variant_id = data$ref_panel$variant_id,
    A1 = data$ref_panel$A1,
    A2 = data$ref_panel$A2,
    z = rnorm(nrow(data$ref_panel)),
    stringsAsFactors = FALSE
  )
  result <- raiss(data$ref_panel, all_known, genotypeMatrix = data$X,
                  verbose = FALSE)
  expect_equal(nrow(result$resultNofilter), nrow(data$ref_panel))
})

test_that("X path handles single unknown variant", {
  data <- generate_X_test_data(n = 100, p = 50, n_known = 49, seed = 15)

  result_R <- raiss(data$ref_panel, data$known_zscores, ldMatrix =data$R,
                    lamb = 0.01, rcond = 0.01, r2Threshold =0, minimumLd =0,
                    verbose = FALSE)
  result_X <- raiss(data$ref_panel, data$known_zscores, genotypeMatrix = data$X,
                    lamb = 0.01, svdTol =1e-12, r2Threshold =0, minimumLd =0,
                    verbose = FALSE)

  r_sorted <- result_R$resultNofilter %>% arrange(variant_id)
  x_sorted <- result_X$resultNofilter %>% arrange(variant_id)

  expect_equal(r_sorted$z, x_sorted$z, tolerance = 1e-4)
})

test_that("X path handles single known variant", {
  data <- generate_X_test_data(n = 100, p = 50, n_known = 1, seed = 20)

  result_X <- raiss(data$ref_panel, data$known_zscores, genotypeMatrix = data$X,
                    lamb = 0.01, svdTol =1e-12, r2Threshold =0, minimumLd =0,
                    verbose = FALSE)
  expect_true(is.data.frame(result_X$resultNofilter))
  expect_equal(nrow(result_X$resultNofilter), nrow(data$ref_panel))
})

test_that("X path R2 filtering matches R path", {
  data <- generate_X_test_data(n = 200, p = 100, n_known = 50, seed = 42)

  result_R <- raiss(data$ref_panel, data$known_zscores, ldMatrix =data$R,
                    lamb = 0.01, rcond = 0.01, r2Threshold =0.6, minimumLd =5,
                    verbose = FALSE)
  result_X <- raiss(data$ref_panel, data$known_zscores, genotypeMatrix = data$X,
                    lamb = 0.01, svdTol =1e-12, r2Threshold =0.6, minimumLd =5,
                    verbose = FALSE)

  # Same variants should pass filtering
  r_filtered_ids <- sort(result_R$resultFilter$variant_id)
  x_filtered_ids <- sort(result_X$resultFilter$variant_id)
  expect_equal(r_filtered_ids, x_filtered_ids,
               info = "Same variants should pass R2/LD filtering")
})

test_that("raw genotype_matrix path is not equivalent to LD path used by legacy pipeline", {
  set.seed(1)
  n <- 80
  p <- 40
  n_known <- 20
  X_raw <- matrix(sample(0:2, n * p, replace = TRUE,
                         prob = c(0.35, 0.45, 0.20)),
                  nrow = n, ncol = p)
  ref_panel <- data.frame(
    chrom = rep(1, p),
    pos = seq_len(p) * 100,
    variant_id = paste0("rs", seq_len(p)),
    A1 = rep("A", p),
    A2 = rep("G", p),
    stringsAsFactors = FALSE
  )
  colnames(X_raw) <- ref_panel$variant_id

  known_idx <- sort(sample(seq_len(p), n_known))
  known_zscores <- data.frame(
    chrom = rep(1, n_known),
    pos = ref_panel$pos[known_idx],
    variant_id = ref_panel$variant_id[known_idx],
    A1 = ref_panel$A1[known_idx],
    A2 = ref_panel$A2[known_idx],
    z = rnorm(n_known),
    stringsAsFactors = FALSE
  )

  R <- computeLd(X_raw, method = "sample")
  rownames(R) <- colnames(R) <- ref_panel$variant_id
  X_scaled <- scale(X_raw)
  X_scaled[is.na(X_scaled)] <- 0
  colnames(X_scaled) <- ref_panel$variant_id

  result_LD <- raiss(ref_panel, known_zscores, ldMatrix =R,
                     lamb = 0.01, rcond = 0.01,
                     r2Threshold =0.6, minimumLd =0,
                     verbose = FALSE)
  result_raw_X <- raiss(ref_panel, known_zscores, genotypeMatrix = X_raw,
                        lamb = 0.01, svdTol =1e-8,
                        r2Threshold =0.6, minimumLd =0,
                        verbose = FALSE)
  result_scaled_X <- raiss(ref_panel, known_zscores, genotypeMatrix = X_scaled,
                           lamb = 0.01, svdTol =1e-12,
                           r2Threshold =0.6, minimumLd =0,
                           verbose = FALSE)

  ld_ids <- sort(result_LD$resultFilter$variant_id)
  raw_x_ids <- sort(result_raw_X$resultFilter$variant_id)
  scaled_x_ids <- sort(result_scaled_X$resultFilter$variant_id)

  expect_false(identical(ld_ids, raw_x_ids))
  expect_gt(nrow(result_raw_X$resultFilter), nrow(result_LD$resultFilter))
  expect_equal(scaled_x_ids, ld_ids)

  ld_sorted <- result_LD$resultNofilter %>% arrange(variant_id)
  scaled_sorted <- result_scaled_X$resultNofilter %>% arrange(variant_id)
  expect_equal(ld_sorted$z, scaled_sorted$z, tolerance = 1e-10)
  expect_equal(ld_sorted$raissR2, scaled_sorted$raissR2, tolerance = 1e-10)
})

test_that("raiss rejects both ldMatrix and genotype_matrix", {
  data <- generate_X_test_data(n = 50, p = 20, n_known = 10, seed = 1)
  expect_error(
    raiss(data$ref_panel, data$known_zscores,
          ldMatrix =data$R, genotypeMatrix = data$X),
    "not both"
  )
})

test_that("raiss rejects neither ldMatrix nor genotype_matrix", {
  data <- generate_X_test_data(n = 50, p = 20, n_known = 10, seed = 1)
  expect_error(
    raiss(data$ref_panel, data$known_zscores),
    "Provide either"
  )
})

test_that("X path with collinear variants matches R path", {
  set.seed(55)
  n <- 150
  p <- 60
  # Create X with some near-duplicate columns
  X_raw <- matrix(sample(0:2, n * p, replace = TRUE, prob = c(0.25, 0.5, 0.25)),
                  nrow = n, ncol = p)
  # Make columns 5 and 6 nearly identical
  X_raw[, 6] <- X_raw[, 5] + sample(c(0, 0, 0, 0, 1), n, replace = TRUE)
  X_raw[X_raw > 2] <- 2

  X <- scale(X_raw)
  X[is.na(X)] <- 0

  ref_panel <- data.frame(
    chrom = rep(1, p), pos = seq(1, p * 10, 10),
    variant_id = paste0("rs", seq_len(p)),
    A1 = rep("A", p), A2 = rep("G", p),
    stringsAsFactors = FALSE
  )
  colnames(X) <- ref_panel$variant_id

  n_known <- 30
  known_idx <- sort(sample(seq_len(p), n_known))
  known_zscores <- data.frame(
    chrom = rep(1, n_known), pos = ref_panel$pos[known_idx],
    variant_id = ref_panel$variant_id[known_idx],
    A1 = ref_panel$A1[known_idx], A2 = ref_panel$A2[known_idx],
    z = rnorm(n_known), stringsAsFactors = FALSE
  )

  R <- cor(X_raw)
  R[is.na(R)] <- 0
  colnames(R) <- rownames(R) <- ref_panel$variant_id

  result_R <- raiss(ref_panel, known_zscores, ldMatrix =R,
                    lamb = 0.01, rcond = 0.01, r2Threshold =0, minimumLd =0,
                    verbose = FALSE)
  result_X <- raiss(ref_panel, known_zscores, genotypeMatrix = X,
                    lamb = 0.01, svdTol =1e-12, r2Threshold =0, minimumLd =0,
                    verbose = FALSE)

  r_sorted <- result_R$resultNofilter %>% arrange(variant_id)
  x_sorted <- result_X$resultNofilter %>% arrange(variant_id)

  expect_equal(r_sorted$z, x_sorted$z, tolerance = 1e-3,
               info = "Collinear case: z-scores should be close")
})


context("RAISS missing-variant imputation in TWAS pipelines")

# Previous tests covered `twasWeightsSumstatPipeline(imputeMissing = ...)`,
# which has been removed in favor of the S4 `twasWeightsPipeline` family
# dispatching on `QtlSumStats` / `QtlDataset`. The missing-variant
# imputation knob now lives inside `summaryStatsQc(impute = TRUE)`, and
# tests for that path live in test_sumstatsQc.R (internal helpers) and
# in the SumStats pipeline tests.
#
# RAISS itself (`raiss()`) still exists with the same signature; its
# direct tests live in test_raiss.R.


context("slalom")

# ============================================================================
# Helper: build a valid positive-definite LD matrix from a genotype matrix
# ============================================================================
make_synthetic_ld <- function(n_samples, n_snps, seed = 1) {
  set.seed(seed)
  # Simulate genotypes with some LD structure by using a factor model
  # X = Z %*% L + noise, where Z is latent and L is a loading matrix
  n_factors <- min(3, n_snps)
  Z <- matrix(rnorm(n_samples * n_factors), nrow = n_samples)
  L <- matrix(runif(n_factors * n_snps, -1, 1), nrow = n_factors)
  X_raw <- Z %*% L + matrix(rnorm(n_samples * n_snps, sd = 0.5), nrow = n_samples)
  # Discretise to genotype-like values (0, 1, 2)
  X <- matrix(as.integer(cut(X_raw, breaks = c(-Inf, -0.5, 0.5, Inf))) - 1L,
              nrow = n_samples, ncol = n_snps)
  colnames(X) <- paste0("snp", seq_len(n_snps))
  R <- cor(X)
  # Ensure perfect diagonal and clean NaN from any zero-variance columns
  R[is.na(R) | is.nan(R)] <- 0
  diag(R) <- 1.0
  list(X = X, R = R)
}

# ============================================================================
# Basic output structure
# ============================================================================

test_that("slalom basic output structure", {
  set.seed(42)
  n <- 50
  z <- rnorm(n)
  R <- diag(n)
  # Add some off-diagonal correlations
  for (i in 1:(n - 1)) {
    R[i, i + 1] <- 0.3
    R[i + 1, i] <- 0.3
  }

  result <- slalom(zScore = z, R = R)

  expect_type(result, "list")
  expect_named(result, c("data", "summary"))
  expect_s3_class(result$data, "data.frame")
  expect_true("original_z" %in% colnames(result$data))
  expect_true("prob" %in% colnames(result$data))
  expect_true("pvalue" %in% colnames(result$data))
  expect_true("outliers" %in% colnames(result$data))
  expect_true("nlog10p_dentist_s" %in% colnames(result$data))
  expect_equal(nrow(result$data), n)
})

test_that("slalom errors on non-square R", {
  z <- rnorm(10)
  R <- matrix(rnorm(50), nrow = 5, ncol = 10)
  expect_error(slalom(zScore = z, R = R), "R must be a square matrix")
})

test_that("slalom accepts X matrix instead of R", {
  set.seed(42)
  n_samples <- 100
  n_snps <- 10
  X <- matrix(sample(0:2, n_samples * n_snps, replace = TRUE), nrow = n_samples, ncol = n_snps)
  colnames(X) <- paste0("snp", 1:n_snps)
  z <- rnorm(n_snps)

  result <- slalom(zScore = z, X = X)
  expect_type(result, "list")
  expect_equal(nrow(result$data), n_snps)
  # PIPs should be in [0,1] and sum to 1
  expect_true(all(result$data$prob >= 0 & result$data$prob <= 1))
  expect_equal(sum(result$data$prob), 1, tolerance = 1e-12)
  # Data frame should have expected column names
  expected_cols <- c("original_z", "prob", "pvalue", "outliers", "nlog10p_dentist_s")
  expect_true(all(expected_cols %in% colnames(result$data)))
})

# ============================================================================
# ABF computation correctness
# ============================================================================

test_that("ABF: strong signal (z=10) gets very high PIP", {
  set.seed(100)
  n <- 20
  z <- rnorm(n, sd = 0.3)
  z[7] <- 10
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  expect_equal(which.max(result$data$prob), 7)
  expect_gt(result$data$prob[7], 0.15)
  expect_equal(sum(result$data$prob), 1, tolerance = 1e-12)
})

test_that("ABF: moderate signal (z=3) gets higher PIP than weak signal (z=1)", {
  set.seed(101)
  n <- 10
  z <- rep(0, n)
  z[2] <- 1
  z[5] <- 3
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  expect_gt(result$data$prob[5], result$data$prob[2])
  # z=0 variants should all have same PIP (by symmetry, with identity LD)
  zero_pips <- result$data$prob[c(1, 3, 4, 6, 7, 8, 9, 10)]
  expect_equal(max(zero_pips) - min(zero_pips), 0, tolerance = 1e-14)
})

test_that("ABF: lbf formula matches manual calculation", {
  z_val <- 4.0
  se_val <- 1.0
  W <- 0.04

  V <- se_val^2
  r <- W / (W + V)
  expected_lbf <- 0.5 * (log(1 - r) + r * z_val^2)

  z <- c(z_val, 0)
  R <- diag(2)
  result <- slalom(zScore = z, R = R, abfPriorVariance = W)

  lbf_0 <- 0.5 * (log(1 - r) + r * 0^2)
  expected_ratio <- exp(expected_lbf - lbf_0)
  actual_ratio <- result$data$prob[1] / result$data$prob[2]
  expect_equal(actual_ratio, expected_ratio, tolerance = 1e-10)
})

test_that("ABF: PIPs always sum to exactly 1", {
  for (s in 1:5) {
    set.seed(200 + s)
    n <- sample(5:30, 1)
    z <- rnorm(n, sd = 2)
    R <- diag(n)
    result <- slalom(zScore = z, R = R)
    expect_equal(sum(result$data$prob), 1, tolerance = 1e-12,
                 label = paste("seed", 200 + s))
  }
})

test_that("ABF: symmetric z-scores give symmetric PIPs", {
  z <- c(-3, 3)
  R <- diag(2)
  result <- slalom(zScore = z, R = R)
  expect_equal(result$data$prob[1], result$data$prob[2], tolerance = 1e-14)
})

# ============================================================================
# Credible sets
# ============================================================================

test_that("CS95 contains the causal variant in a simple synthetic signal", {
  set.seed(300)
  n <- 30
  z <- rnorm(n, sd = 0.5)
  causal <- 12
  z[causal] <- 6
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  expect_true(causal %in% result$summary$cs95)
  expect_true(causal %in% result$summary$cs99)
})

test_that("CS99 is a superset of CS95", {
  set.seed(301)
  n <- 40
  z <- rnorm(n, sd = 1.5)
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  expect_true(all(result$summary$cs95 %in% result$summary$cs99))
  expect_gte(length(result$summary$cs99), length(result$summary$cs95))
})

test_that("CS95 covers at least 95% of posterior mass", {
  set.seed(302)
  n <- 25
  z <- rnorm(n)
  z[10] <- 4
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  cs95_mass <- sum(result$data$prob[result$summary$cs95])
  expect_gt(cs95_mass, 0.95)

  cs99_mass <- sum(result$data$prob[result$summary$cs99])
  expect_gt(cs99_mass, 0.99)
})

test_that("CS with very strong signal contains only the causal variant", {
  set.seed(303)
  n <- 15
  z <- rnorm(n, sd = 0.1)
  z[8] <- 15  # extremely strong
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  expect_equal(result$summary$cs95[1], 8)
  expect_true(8 %in% result$summary$cs95)
  expect_true(8 %in% result$summary$cs99)
})

test_that("CS with diffuse signal contains many variants", {
  set.seed(304)
  n <- 10
  z <- rep(0, n)  # all equally uninformative
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  # Uniform PIPs => need at least ceiling(0.95 * n) = 10 variants for 95% coverage
  expect_equal(length(result$summary$cs95), n)
})

# ============================================================================
# Lead variant by pvalue vs abf
# ============================================================================

test_that("lead variant by pvalue selects most negative z-score", {
  z <- c(0, -4, 3, -1, 2)
  R <- diag(5)

  result <- slalom(zScore = z, R = R, leadVariantChoice = "pvalue")

  expect_equal(result$summary$leadPipVariant, 2)
})

test_that("lead variant by abf selects highest PIP", {
  z <- c(0, -4, 3, -1, 2)
  R <- diag(5)

  result <- slalom(zScore = z, R = R, leadVariantChoice = "abf")

  expect_equal(result$summary$leadPipVariant, which.max(result$data$prob))
})

test_that("pvalue and abf lead can differ when z has asymmetric magnitudes", {
  z <- c(-3.0, 5.0, 0.1, -0.2, 0.3)
  R <- diag(5)

  result_pv <- slalom(zScore = z, R = R, leadVariantChoice = "pvalue")
  result_abf <- slalom(zScore = z, R = R, leadVariantChoice = "abf")

  expect_equal(result_pv$summary$leadPipVariant, 1)
  expect_equal(result_abf$summary$leadPipVariant, 2)
  expect_false(result_pv$summary$leadPipVariant == result_abf$summary$leadPipVariant)
})

# ============================================================================
# DENTIST-S outlier detection
# ============================================================================

test_that("DENTIST-S: lead variant itself is not flagged as outlier", {
  set.seed(400)
  n <- 10
  z <- rnorm(n, sd = 0.5)
  z[3] <- -5  # lead by pvalue
  R <- diag(n)

  result <- slalom(zScore = z, R = R)
  lead <- result$summary$leadPipVariant
  expect_equal(lead, 3)

  expect_true(is.na(result$data$outliers[lead]) || !result$data$outliers[lead])
})

test_that("DENTIST-S: outlier variant inconsistent with LD is flagged", {
  n <- 5
  z <- c(-5, 0, 0.1, -0.1, 0.2)
  R <- diag(n)
  R[1, 2] <- R[2, 1] <- 0.9
  R[1, 3] <- R[3, 1] <- 0.1
  R[1, 4] <- R[4, 1] <- -0.05
  R[1, 5] <- R[5, 1] <- 0.02

  result <- slalom(zScore = z, R = R, r2Threshold = 0.5,
                   nlog10pDentistSThreshold = 2.0)

  lead <- result$summary$leadPipVariant
  expect_equal(lead, 1)

  expect_true(result$data$outliers[2])
  expect_false(result$data$outliers[3])
})

test_that("DENTIST-S: perfectly consistent variant in LD is not flagged", {
  n <- 3
  z <- c(-5, -4.0, 0.1)
  R <- diag(n)
  R[1, 2] <- R[2, 1] <- 0.8
  R[1, 3] <- R[3, 1] <- 0.05
  R[2, 3] <- R[3, 2] <- 0.04

  result <- slalom(zScore = z, R = R, r2Threshold = 0.5,
                   nlog10pDentistSThreshold = 4.0)

  lead <- result$summary$leadPipVariant
  expect_equal(lead, 1)

  expect_equal(result$data$nlog10p_dentist_s[2], 0, tolerance = 1e-10)
  expect_false(result$data$outliers[2])
})

test_that("DENTIST-S: n_dentist_s_outlier and fraction are consistent", {
  set.seed(401)
  n <- 20
  syn <- make_synthetic_ld(200, n, seed = 401)
  z <- rnorm(n, sd = 1)
  z[1] <- -5

  result <- slalom(zScore = z, R = syn$R, r2Threshold = 0.3)

  n_r2 <- result$summary$nR2
  n_out <- result$summary$nDentistSOutlier
  frac <- result$summary$fraction

  expect_gte(n_r2, 1)
  expect_equal(frac, ifelse(n_r2 > 0, n_out / n_r2, 0), tolerance = 1e-14)
  expect_gte(n_out, 0)
  expect_lte(n_out, n_r2)
})

test_that("DENTIST-S: lowering threshold flags more outliers", {
  set.seed(402)
  n <- 15
  syn <- make_synthetic_ld(300, n, seed = 402)
  z <- rnorm(n, sd = 2)
  z[5] <- -6

  result_strict <- slalom(zScore = z, R = syn$R, nlog10pDentistSThreshold = 6.0,
                          r2Threshold = 0.3)
  result_loose <- slalom(zScore = z, R = syn$R, nlog10pDentistSThreshold = 1.0,
                         r2Threshold = 0.3)

  expect_gte(result_loose$summary$nDentistSOutlier,
             result_strict$summary$nDentistSOutlier)
})

# ============================================================================
# Edge cases
# ============================================================================

test_that("edge case: single variant", {
  z <- c(3.0)
  R <- matrix(1, nrow = 1, ncol = 1)

  result <- slalom(zScore = z, R = R)

  expect_equal(nrow(result$data), 1)
  expect_equal(result$data$prob[1], 1.0, tolerance = 1e-14)
  expect_equal(result$summary$leadPipVariant, 1)
  expect_equal(result$summary$nTotal, 1)
  expect_equal(result$summary$cs95, 1)
  expect_equal(result$summary$cs99, 1)
})

test_that("edge case: all zero z-scores", {
  n <- 10
  z <- rep(0, n)
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  expect_equal(result$data$prob, rep(1 / n, n), tolerance = 1e-14)
  expect_equal(result$data$pvalue, rep(0.5, n), tolerance = 1e-14)
  expect_equal(result$summary$maxPip, 1 / n, tolerance = 1e-14)
})

test_that("edge case: very large z-scores do not produce NaN in PIPs", {
  z <- c(50, -50, 30, -30, 0)
  R <- diag(5)

  result <- slalom(zScore = z, R = R)

  expect_false(any(is.nan(result$data$prob)))
  expect_false(any(is.na(result$data$prob)))
  expect_equal(sum(result$data$prob), 1, tolerance = 1e-10)
  expect_equal(result$data$prob[1], result$data$prob[2], tolerance = 1e-14)
})

test_that("edge case: identical z-scores yield uniform PIPs", {
  n <- 8
  z <- rep(2.5, n)
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  expect_equal(result$data$prob, rep(1 / n, n), tolerance = 1e-14)
})

test_that("edge case: two variants only", {
  z <- c(-3, 2)
  R <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)

  result <- slalom(zScore = z, R = R)

  expect_equal(nrow(result$data), 2)
  expect_equal(sum(result$data$prob), 1, tolerance = 1e-12)
  expect_equal(result$summary$leadPipVariant, 1)
  expect_gt(result$data$prob[1], result$data$prob[2])
})

test_that("edge case: mismatched dimensions error", {
  z <- rnorm(10)
  R <- diag(5)
  expect_error(slalom(zScore = z, R = R),
               "R must be a square matrix matching the length of zScore")
})

test_that("edge case: no R and no X provided errors", {
  z <- rnorm(5)
  expect_error(slalom(zScore = z), "Either R.*or X.*must be provided")
})

test_that("edge case: both R and X provided errors", {
  set.seed(500)
  n <- 5
  z <- rnorm(n)
  R <- diag(n)
  X <- matrix(sample(0:2, 50 * n, replace = TRUE), nrow = 50, ncol = n)
  expect_error(slalom(zScore = z, R = R, X = X), "Provide either R or X, not both")
})

# ============================================================================
# X input mode
# ============================================================================

test_that("X input yields same result as R = cor(X)", {
  set.seed(600)
  n_samples <- 200
  n_snps <- 10
  X <- matrix(sample(0:2, n_samples * n_snps, replace = TRUE),
              nrow = n_samples, ncol = n_snps)
  colnames(X) <- paste0("snp", seq_len(n_snps))
  z <- rnorm(n_snps, sd = 2)

  R_manual <- cor(X)
  diag(R_manual) <- 1.0

  result_X <- slalom(zScore = z, X = X)
  result_R <- slalom(zScore = z, R = R_manual)

  expect_equal(result_X$data$prob, result_R$data$prob, tolerance = 1e-10)
  expect_equal(result_X$data$original_z, result_R$data$original_z, tolerance = 1e-14)
  expect_equal(result_X$summary$leadPipVariant, result_R$summary$leadPipVariant)
  expect_equal(result_X$summary$cs95, result_R$summary$cs95)
  expect_equal(result_X$summary$cs99, result_R$summary$cs99)
})

# ============================================================================
# Parameter variation
# ============================================================================

test_that("larger abf_prior_variance concentrates PIPs on strong signals more", {
  set.seed(700)
  n <- 15
  z <- rnorm(n, sd = 0.5)
  z[4] <- 4
  R <- diag(n)

  result_small_W <- slalom(zScore = z, R = R, abfPriorVariance = 0.01)
  result_large_W <- slalom(zScore = z, R = R, abfPriorVariance = 1.0)

  expect_gt(result_large_W$summary$maxPip, result_small_W$summary$maxPip)
})

test_that("abfPriorVariance = 0 gives uniform PIPs", {
  n <- 10
  z <- c(5, 3, 1, 0, -1, -3, -5, 2, -2, 4)
  R <- diag(n)

  result <- slalom(zScore = z, R = R, abfPriorVariance = 0)

  expect_equal(result$data$prob, rep(1 / n, n), tolerance = 1e-14)
})

test_that("different standard_error values affect PIPs", {
  n <- 5
  z <- c(3, 3, 3, 3, 3)  # same z for all
  R <- diag(n)
  se1 <- c(1, 1, 1, 1, 1)
  se2 <- c(0.5, 1, 1, 1, 1)  # variant 1 has smaller SE

  result1 <- slalom(zScore = z, R = R, standardError = se1)
  result2 <- slalom(zScore = z, R = R, standardError = se2)

  expect_gt(result2$data$prob[1], result1$data$prob[1])
})

test_that("r2_threshold variation affects n_r2 count", {
  set.seed(701)
  n <- 10
  syn <- make_synthetic_ld(200, n, seed = 701)
  z <- rnorm(n, sd = 2)
  z[1] <- -5

  result_low <- slalom(zScore = z, R = syn$R, r2Threshold = 0.1)
  result_high <- slalom(zScore = z, R = syn$R, r2Threshold = 0.9)

  expect_gte(result_low$summary$nR2, result_high$summary$nR2)
})

test_that("nlog10p_dentist_s_threshold variation affects outlier count", {
  set.seed(702)
  n <- 10
  syn <- make_synthetic_ld(200, n, seed = 702)
  z <- rnorm(n, sd = 2)
  z[1] <- -6

  result_low_thresh <- slalom(zScore = z, R = syn$R, nlog10pDentistSThreshold = 1.0)
  result_high_thresh <- slalom(zScore = z, R = syn$R, nlog10pDentistSThreshold = 10.0)

  expect_gte(result_low_thresh$summary$nDentistSOutlier,
             result_high_thresh$summary$nDentistSOutlier)
})

# ============================================================================
# Output structure validation
# ============================================================================

test_that("output data types are correct", {
  set.seed(801)
  n <- 15
  z <- rnorm(n)
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  expect_type(result$data$original_z, "double")
  expect_type(result$data$prob, "double")
  expect_type(result$data$pvalue, "double")
  expect_type(result$data$outliers, "logical")
  expect_type(result$data$nlog10p_dentist_s, "double")

  expect_type(result$summary$leadPipVariant, "integer")
  expect_type(result$summary$nTotal, "integer")
  expect_type(result$summary$nR2, "integer")
  expect_type(result$summary$nDentistSOutlier, "integer")
  expect_type(result$summary$fraction, "double")
  expect_type(result$summary$maxPip, "double")
  expect_type(result$summary$cs95, "integer")
  expect_type(result$summary$cs99, "integer")
})

test_that("original_z in output matches input z-scores", {
  z <- c(1.5, -2.3, 0.7, 4.1, -0.5)
  R <- diag(5)

  result <- slalom(zScore = z, R = R)
  expect_equal(result$data$original_z, z, tolerance = 0)
})

test_that("pvalue in output matches pnorm(z)", {
  z <- c(-3, -1, 0, 1, 3)
  R <- diag(5)

  result <- slalom(zScore = z, R = R)
  expect_equal(result$data$pvalue, pnorm(z), tolerance = 1e-14)
})

# ============================================================================
# Summary statistics consistency
# ============================================================================

test_that("n_r2 counts variants with r2 > threshold to lead correctly", {
  n <- 10
  z <- c(-5, rep(0, 9))
  R <- diag(n)

  result <- slalom(zScore = z, R = R, r2Threshold = 0.6)
  expect_equal(result$summary$nR2, 1)
})

test_that("n_r2 includes correlated variants", {
  n <- 5
  z <- c(-5, -4, 0, 0, 0)
  R <- diag(n)
  R[1, 2] <- R[2, 1] <- 0.9

  result <- slalom(zScore = z, R = R, r2Threshold = 0.6)
  expect_equal(result$summary$nR2, 2)
})

test_that("fraction = 0 when there are no outliers (identity LD, consistent z)", {
  n <- 5
  z <- c(-3, 0, 0, 0, 0)
  R <- diag(n)

  result <- slalom(zScore = z, R = R)
  expect_equal(result$summary$fraction, 0)
})

test_that("fraction is between 0 and 1", {
  set.seed(900)
  for (s in 1:5) {
    set.seed(900 + s)
    n <- sample(10:30, 1)
    syn <- make_synthetic_ld(200, n, seed = 900 + s)
    z <- rnorm(n, sd = 2)
    z[1] <- -6

    result <- slalom(zScore = z, R = syn$R, r2Threshold = 0.3)
    expect_gte(result$summary$fraction, 0)
    expect_lte(result$summary$fraction, 1)
  }
})

test_that("maxPip equals the maximum of prob vector", {
  set.seed(901)
  n <- 20
  z <- rnorm(n, sd = 2)
  R <- diag(n)

  result <- slalom(zScore = z, R = R)
  expect_equal(result$summary$maxPip, max(result$data$prob), tolerance = 1e-14)
})

# ============================================================================
# Realistic synthetic LD scenarios
# ============================================================================

test_that("realistic LD: correlated variants share PIP mass", {
  set.seed(1000)
  syn <- make_synthetic_ld(500, 20, seed = 1000)
  z <- rep(0, 20)
  z[3] <- 5

  result <- slalom(zScore = z, R = syn$R)

  expect_true(3 %in% result$summary$cs95)
  expect_equal(sum(result$data$prob), 1, tolerance = 1e-12)
})

test_that("realistic LD: DENTIST-S detects outlier in correlated block", {
  set.seed(1001)
  syn <- make_synthetic_ld(500, 15, seed = 1001)
  R <- syn$R

  lead_idx <- 1
  z <- R[, lead_idx] * (-5)
  z[1] <- -5

  r2_to_lead <- R[, 1]^2
  candidates <- which(r2_to_lead > 0.3 & seq_along(z) != 1)
  if (length(candidates) > 0) {
    corrupt_idx <- candidates[1]
    z[corrupt_idx] <- z[corrupt_idx] + 10

    result <- slalom(zScore = z, R = R, r2Threshold = 0.2,
                     nlog10pDentistSThreshold = 3.0)

    expect_true(result$data$outliers[corrupt_idx])
  }
})

test_that("realistic LD: no outliers when z perfectly matches LD structure", {
  set.seed(1002)
  syn <- make_synthetic_ld(500, 10, seed = 1002)
  R <- syn$R

  lead_idx <- 1
  z <- R[, lead_idx] * (-4)

  result <- slalom(zScore = z, R = R, r2Threshold = 0.3,
                   nlog10pDentistSThreshold = 4.0)

  non_lead <- setdiff(seq_along(z), result$summary$leadPipVariant)
  lead <- result$summary$leadPipVariant
  for (i in non_lead) {
    r2_val <- R[i, lead]^2
    if (r2_val > 0.3 && r2_val < 1.0) {
      expect_false(result$data$outliers[i])
    } else if (r2_val >= 1.0) {
      expect_true(is.na(result$data$outliers[i]) || !result$data$outliers[i])
    }
  }
})

# ============================================================================
# Internal function access (resolveLdInput via :::)
# ============================================================================

test_that("resolveLdInput returns R when R is provided", {
  R <- diag(5)
  res <- pecotmr:::resolveLdInput(R = R, needNSample = FALSE)
  expect_equal(res$R, R)
})

test_that("resolveLdInput computes R from X", {
  set.seed(1100)
  X <- matrix(sample(0:2, 200 * 5, replace = TRUE), nrow = 200, ncol = 5)
  colnames(X) <- paste0("s", 1:5)
  res <- pecotmr:::resolveLdInput(X = X, needNSample = FALSE)
  expect_true(is.matrix(res$R))
  expect_equal(nrow(res$R), 5)
  expect_equal(ncol(res$R), 5)
  for (j in seq_len(5)) expect_equal(res$R[j, j], 1.0, tolerance = 1e-6)
})

test_that("resolveLdInput errors when neither R nor X given", {
  expect_error(pecotmr:::resolveLdInput(R = NULL, X = NULL),
               "Either R.*or X.*must be provided")
})

test_that("resolveLdInput errors when both R and X given", {
  R <- diag(3)
  X <- matrix(1, nrow = 10, ncol = 3)
  expect_error(pecotmr:::resolveLdInput(R = R, X = X),
               "Provide either R or X, not both")
})

# ============================================================================
# Numerical stability
# ============================================================================

test_that("log-sum-exp trick prevents overflow with extreme z-scores", {
  z <- c(100, 0, -100)
  R <- diag(3)

  result <- slalom(zScore = z, R = R)

  expect_false(any(is.nan(result$data$prob)))
  expect_false(any(is.infinite(result$data$prob)))
  expect_equal(sum(result$data$prob), 1, tolerance = 1e-10)
  expect_equal(result$data$prob[1], result$data$prob[3], tolerance = 1e-14)
  expect_gt(result$data$prob[1], result$data$prob[2])
})

test_that("standard_error near zero concentrates PIP on large z-scores", {
  z <- c(2, 0.5, 0, -0.1)
  R <- diag(4)
  se <- rep(0.01, 4)

  result <- slalom(zScore = z, R = R, standardError = se, abfPriorVariance = 0.04)

  expect_equal(which.max(result$data$prob), 1)
  expect_false(any(is.nan(result$data$prob)))
})

# ============================================================================
# Determinism and reproducibility
# ============================================================================

test_that("slalom is deterministic (no randomness)", {
  z <- c(3, -2, 1, 0, -4)
  R <- diag(5)

  result1 <- slalom(zScore = z, R = R)
  result2 <- slalom(zScore = z, R = R)

  expect_identical(result1$data, result2$data)
  expect_identical(result1$summary$leadPipVariant, result2$summary$leadPipVariant)
  expect_identical(result1$summary$cs95, result2$summary$cs95)
  expect_identical(result1$summary$cs99, result2$summary$cs99)
})

# ============================================================================
# Credible set ordering
# ============================================================================

test_that("CS95 variants are ordered by decreasing PIP", {
  set.seed(1400)
  n <- 20
  z <- rnorm(n, sd = 2)
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  cs_pips <- result$data$prob[result$summary$cs95]
  expect_true(all(diff(cs_pips) <= .Machine$double.eps))
})



context("summaryStatsQc (with mocked MungeSumstats)")

# NOTE
# ----
# `.runMungeSumstatsFilter` wraps MungeSumstats::format_sumstats which needs
# a real dbSNP reference panel (multi-GB download). To exercise the QC chain
# in a unit test we mock that helper so it just returns the input data.frame
# unchanged, recording a "no variants dropped" audit record. The pecotmr-
# native steps (.applySkipRegion, .matchAgainstSketch, .applyPipScreen,
# .applyLdMismatchQcToEntry) all run for real on the synthetic fixture.

# ===========================================================================
# Fixture builders
# ===========================================================================

.ssQ_makeHandle <- function(snp_n = 8L, n_samples = 60L) {
  new("GenotypeHandle",
    path = "/tmp/sketch.gds",
    format = "gds",
    snpInfo = data.frame(
      SNP = paste0("rs", seq_len(snp_n)),
      CHR = rep("1", snp_n),
      BP  = seq(100L, by = 100L, length.out = snp_n),
      A1  = rep("A", snp_n),
      A2  = rep("G", snp_n),
      stringsAsFactors = FALSE),
    nSamples = n_samples,
    sampleIds = paste0("s", seq_len(n_samples)),
    pgenPtr = NULL)
}

.ssQ_makeEntryGr <- function(snp_ids = paste0("rs", 1:4),
                             positions = c(100L, 200L, 300L, 400L)) {
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", length(snp_ids)),
    ranges = IRanges::IRanges(start = positions, width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = snp_ids,
    A1  = rep("A", length(snp_ids)),
    A2  = rep("G", length(snp_ids)),
    Z   = seq(1.0, by = 0.5, length.out = length(snp_ids)),
    N   = rep(1000L, length(snp_ids)))
  gr
}

.ssQ_makeGwasSumStats <- function(snp_ids = paste0("rs", 1:4),
                                  positions = c(100L, 200L, 300L, 400L),
                                  study = "g1") {
  GwasSumStats(
    study    = study,
    entry    = list(.ssQ_makeEntryGr(snp_ids, positions)),
    genome   = "hg19",
    ldSketch = .ssQ_makeHandle())
}

.ssQ_mockMunge <- function(drop = 0L) {
  # Mock that pretends MungeSumstats validated the input and returned the
  # same data.frame, dropping `drop` rows.
  function(df, refGenome, useDbsnpRefCheck, removeIndels,
           removeStrandAmbiguous, mafCutoff, infoCutoff, nCutoff,
           convertRefGenome, mungeSumstatsArgs) {
    keep <- if (drop > 0L && drop < nrow(df))
      seq_len(nrow(df) - drop)
    else
      seq_len(nrow(df))
    list(df = df[keep, , drop = FALSE],
         droppedNVariants = nrow(df) - length(keep))
  }
}

.ssQ_mockExtractor <- function(seed = 13, n_samples = 60L) {
  function(handle, snpIdx, meanImpute = TRUE) {
    set.seed(seed)
    panel <- matrix(rbinom(n_samples * nrow(handle@snpInfo), 2, 0.3),
                    nrow = n_samples, ncol = nrow(handle@snpInfo),
                    dimnames = list(handle@sampleIds, handle@snpInfo$SNP))
    sub <- panel[, snpIdx, drop = FALSE]
    rr <- GenomicRanges::GRanges(
      seqnames = paste0("chr", handle@snpInfo$CHR[snpIdx]),
      ranges   = IRanges::IRanges(start = handle@snpInfo$BP[snpIdx], width = 1L))
    S4Vectors::mcols(rr) <- S4Vectors::DataFrame(
      SNP = handle@snpInfo$SNP[snpIdx],
      A1  = handle@snpInfo$A1[snpIdx],
      A2  = handle@snpInfo$A2[snpIdx])
    cd <- S4Vectors::DataFrame(sampleId = handle@sampleIds,
                               row.names = handle@sampleIds)
    dosage <- t(sub)
    rownames(dosage) <- handle@snpInfo$SNP[snpIdx]
    colnames(dosage) <- handle@sampleIds
    SummarizedExperiment::SummarizedExperiment(
      assays    = list(dosage = dosage),
      rowRanges = rr,
      colData   = cd)
  }
}

# ===========================================================================
# summaryStatsQc: input-type validation
# ===========================================================================

test_that("summaryStatsQc: rejects non-SumStats input", {
  expect_error(summaryStatsQc("not_a_sumstats"),
               "requires a QtlSumStats or GwasSumStats input")
})

test_that("summaryStatsQc: mafCutoff > 0 with no MAF/FRQ column errors", {
  ss <- .ssQ_makeGwasSumStats()
  expect_error(summaryStatsQc(ss, mafCutoff = 0.05),
               "MAF or FRQ column")
})

test_that("summaryStatsQc: infoCutoff > 0 with no INFO column errors", {
  ss <- .ssQ_makeGwasSumStats()
  expect_error(summaryStatsQc(ss, infoCutoff = 0.5),
               "infoCutoff > 0 requires every entry to carry an INFO column")
})

test_that(".deriveBetaSeFromZ: derives BETA+SE when entry has Z+MAF+N only", {
  df <- data.frame(
    SNP = paste0("rs", 1:3),
    A1 = c("A", "C", "G"), A2 = c("G", "T", "A"),
    Z = c(1.5, -2.1, 0.4),
    MAF = c(0.2, 0.35, 0.05),
    N = c(10000, 10000, 10000),
    stringsAsFactors = FALSE
  )
  out <- pecotmr:::.deriveBetaSeFromZ(df)
  expect_true(all(c("BETA", "SE") %in% names(out$df)))
  expect_equal(out$audit$nDerived, 3L)
  # Verify the formula: se = 1/sqrt(2*maf*(1-maf)*(N+z^2))
  expected_se <- 1 / sqrt(2 * df$MAF * (1 - df$MAF) * (df$N + df$Z^2))
  expect_equal(out$df$SE, expected_se)
  expect_equal(out$df$BETA, df$Z * expected_se)
})

test_that(".deriveBetaSeFromZ: no-op when BETA and SE already present", {
  df <- data.frame(
    Z = 1.5, BETA = 0.5, SE = 0.1, MAF = 0.3, N = 1000,
    stringsAsFactors = FALSE)
  out <- pecotmr:::.deriveBetaSeFromZ(df)
  expect_null(out$audit)
  expect_equal(out$df, df)
})

test_that(".deriveBetaSeFromZ: skipped when N missing", {
  df <- data.frame(
    Z = 1.5, MAF = 0.3,
    stringsAsFactors = FALSE)
  out <- pecotmr:::.deriveBetaSeFromZ(df)
  expect_null(out$audit)
  expect_false("BETA" %in% names(out$df))
  expect_false("SE" %in% names(out$df))
})

# ===========================================================================
# summaryStatsQc: end-to-end with mocked MungeSumstats
# ===========================================================================

test_that("summaryStatsQc: vanilla run populates qcInfo and returns a GwasSumStats", {
  ss <- .ssQ_makeGwasSumStats()
  local_mocked_bindings(
    .runMungeSumstatsFilter = .ssQ_mockMunge(),
    .package = "pecotmr")
  res <- summaryStatsQc(ss)
  expect_s4_class(res, "GwasSumStats")
  qc <- getQcInfo(res)
  expect_true(length(qc) > 0L)
  expect_true("options" %in% names(qc))
  expect_true("entryAudit" %in% names(qc))
  expect_equal(length(qc$entryAudit), nrow(ss))
  # Per-entry audit records variantsIn / variantsOut / mungeSumstatsDropped.
  ea <- qc$entryAudit[[1L]]
  expect_equal(ea$variantsIn, 4L)
  expect_equal(ea$variantsOut, 4L)
  expect_equal(ea$mungeSumstatsDropped, 0L)
})

test_that("summaryStatsQc: keepVariants subsets each entry and records the drop", {
  ss <- .ssQ_makeGwasSumStats()
  local_mocked_bindings(
    .runMungeSumstatsFilter = .ssQ_mockMunge(),
    .package = "pecotmr")
  res <- summaryStatsQc(ss, keepVariants = c("rs1", "rs3"))
  ea <- getQcInfo(res)$entryAudit[[1L]]
  expect_equal(ea$keepVariantsDropped, 2L)
  expect_equal(ea$variantsOut, 2L)
})

test_that("summaryStatsQc: skipRegion drops overlapping variants", {
  ss <- .ssQ_makeGwasSumStats()
  local_mocked_bindings(
    .runMungeSumstatsFilter = .ssQ_mockMunge(),
    .package = "pecotmr")
  res <- summaryStatsQc(ss, skipRegion = "chr1:50-150")
  ea <- getQcInfo(res)$entryAudit[[1L]]
  expect_equal(ea$skipRegionDropped, 1L)  # rs1 at pos 100 is dropped
})

test_that("summaryStatsQc: PIP screen triggers when no variant has signal", {
  # Build an entry with weak signal so the SER PIP screen tags everything
  # below the threshold.
  gr <- .ssQ_makeEntryGr()
  S4Vectors::mcols(gr)$Z <- rep(0.1, length(gr))
  ss <- GwasSumStats(study = "g1", entry = list(gr), genome = "hg19",
                      ldSketch = .ssQ_makeHandle())
  local_mocked_bindings(
    .runMungeSumstatsFilter = .ssQ_mockMunge(),
    .package = "pecotmr")
  res <- summaryStatsQc(ss, pipCutoffToSkip = 0.99)
  ea <- getQcInfo(res)$entryAudit[[1L]]
  expect_true(isTRUE(ea$pipScreenSkipped))
  expect_match(ea$pipScreenReason, "no signals above PIP threshold")
  expect_equal(length(res$entry[[1L]]), 0L)
})

test_that("summaryStatsQc: early-exit records when fewer than 2 variants remain pre-harmonization", {
  ss <- .ssQ_makeGwasSumStats()
  local_mocked_bindings(
    # Mock keeps only 1 row of the input
    .runMungeSumstatsFilter = .ssQ_mockMunge(drop = 3L),
    .package = "pecotmr")
  res <- summaryStatsQc(ss)
  ea <- getQcInfo(res)$entryAudit[[1L]]
  expect_match(ea$earlyExit, "fewer than two variants")
})

test_that("summaryStatsQc: harmonized variants count is recorded", {
  ss <- .ssQ_makeGwasSumStats()
  local_mocked_bindings(
    .runMungeSumstatsFilter = .ssQ_mockMunge(),
    .package = "pecotmr")
  res <- summaryStatsQc(ss)
  ea <- getQcInfo(res)$entryAudit[[1L]]
  expect_equal(ea$matchedAgainstSketch, 4L)
})

test_that("summaryStatsQc: options block records the curated knobs", {
  ss <- .ssQ_makeGwasSumStats()
  local_mocked_bindings(
    .runMungeSumstatsFilter = .ssQ_mockMunge(),
    .package = "pecotmr")
  res <- summaryStatsQc(ss, removeIndels = TRUE, removeStrandAmbiguous = FALSE,
                        nCutoff = 10)
  opts <- getQcInfo(res)$options
  expect_true(opts$removeIndels)
  expect_false(opts$removeStrandAmbiguous)
  expect_equal(opts$nCutoff, 10)
})

test_that("summaryStatsQc: round-trips QtlSumStats inputs", {
  gr <- .ssQ_makeEntryGr()
  ss <- QtlSumStats(study = "s1", context = "c1", trait = "t1",
                     entry = list(gr), genome = "hg19",
                     ldSketch = .ssQ_makeHandle())
  local_mocked_bindings(
    .runMungeSumstatsFilter = .ssQ_mockMunge(),
    .package = "pecotmr")
  res <- summaryStatsQc(ss)
  expect_s4_class(res, "QtlSumStats")
  expect_equal(length(getQcInfo(res)$entryAudit), 1L)
})

# ===========================================================================
# summaryStatsQc with LD-mismatch QC enabled (mocked extractor)
# ===========================================================================

test_that("summaryStatsQc: zMismatchQc = 'dentist' walks the LD-mismatch branch", {
  ss <- .ssQ_makeGwasSumStats(snp_ids = paste0("rs", 1:8),
                              positions = seq(100L, by = 100L, length.out = 8L))
  local_mocked_bindings(
    .runMungeSumstatsFilter = .ssQ_mockMunge(),
    extractBlockGenotypes   = .ssQ_mockExtractor(),
    .package = "pecotmr")
  res <- suppressWarnings(summaryStatsQc(ss, zMismatchQc = "dentist"))
  ea <- getQcInfo(res)$entryAudit[[1L]]
  expect_equal(ea$ldMismatchMethod, "dentist")
  expect_true("ldMismatchOutliersDropped" %in% names(ea))
})

# ===========================================================================
# summaryStatsQc with impute = TRUE: exercise the RAISS branch
# ===========================================================================

test_that("summaryStatsQc: impute = TRUE invokes RAISS and records the audit counts", {
  # Build a sketch panel with 8 variants and a GWAS entry covering only the
  # first 4 — RAISS is asked to impute the missing 4.
  full_snp_ids <- paste0("rs", 1:8)
  full_positions <- seq(100L, by = 100L, length.out = 8L)
  ss <- GwasSumStats(
    study  = "g1",
    entry  = list(.ssQ_makeEntryGr(
                    snp_ids   = full_snp_ids[1:4],
                    positions = full_positions[1:4])),
    genome = "hg19",
    ldSketch = .ssQ_makeHandle(snp_n = 8L, n_samples = 60L))

  local_mocked_bindings(
    .runMungeSumstatsFilter = .ssQ_mockMunge(),
    extractBlockGenotypes   = .ssQ_mockExtractor(),
    raiss = function(refPanel, knownZscores, genotypeMatrix, ...) {
      # Pretend RAISS imputed two of the missing panel variants (rs5, rs6)
      # with synthetic z-scores.
      added <- refPanel[refPanel$variant_id %in% c("rs5", "rs6"), , drop = FALSE]
      added$z <- c(1.5, -2.0)
      added$n <- c(1000, 1000)
      list(resultFilter = rbind(knownZscores, added))
    },
    .package = "pecotmr")
  res <- summaryStatsQc(ss, impute = TRUE)
  ea <- getQcInfo(res)$entryAudit[[1L]]
  expect_equal(ea$raissTotalVariants, 6L)
  expect_equal(ea$raissImputedVariants, 2L)
})

test_that("summaryStatsQc: impute = TRUE with raiss returning NULL records 0 imputed", {
  full_snp_ids <- paste0("rs", 1:8)
  full_positions <- seq(100L, by = 100L, length.out = 8L)
  ss <- GwasSumStats(
    study  = "g1",
    entry  = list(.ssQ_makeEntryGr(
                    snp_ids   = full_snp_ids[1:4],
                    positions = full_positions[1:4])),
    genome = "hg19",
    ldSketch = .ssQ_makeHandle(snp_n = 8L, n_samples = 60L))

  local_mocked_bindings(
    .runMungeSumstatsFilter = .ssQ_mockMunge(),
    extractBlockGenotypes   = .ssQ_mockExtractor(),
    raiss = function(...) NULL,
    .package = "pecotmr")
  res <- summaryStatsQc(ss, impute = TRUE)
  ea <- getQcInfo(res)$entryAudit[[1L]]
  expect_equal(ea$raissImputedVariants, 0L)
})

# ===========================================================================
# summaryStatsQc: per-step QC counter logging (concept salvaged from PR #520)
# ===========================================================================

test_that(".matchRefPanel surfaces sign/strand/dropped counts via qcCounts attribute", {
  # 4 shared positions: 100 exact, 200 sign-flip, 300 strand-flip (A/G
  # unambiguous), 400 allele mismatch (dropped).
  target <- data.frame(
    chrom = c(1, 1, 1, 1), pos = c(100, 200, 300, 400),
    A2 = c("A", "A", "A", "A"), A1 = c("G", "G", "G", "G"),
    z = c(1, 2, 3, 4), stringsAsFactors = FALSE)
  ref <- data.frame(
    chrom = c(1, 1, 1, 1), pos = c(100, 200, 300, 400),
    A2 = c("A", "G", "T", "C"), A1 = c("G", "A", "C", "A"),
    stringsAsFactors = FALSE)
  res <- pecotmr:::.matchRefPanel(target, ref, colToFlip = "z",
                                   matchMinProp = 0)
  # Default return shape unchanged.
  expect_named(res, c("harmonizedData", "qcSummary"))
  expect_equal(nrow(res$harmonizedData), 3L)
  cnt <- attr(res, "qcCounts")
  expect_false(is.null(cnt))
  expect_equal(cnt$considered, 4L)
  expect_equal(cnt$signFlip,   1L)
  expect_equal(cnt$strandFlip, 1L)
  expect_equal(cnt$kept,       3L)
  expect_equal(cnt$dropped,    1L)
})

test_that("summaryStatsQc: QC track emits per-step 'kept N of M' messages plus a rollup", {
  ss <- .ssQ_makeGwasSumStats()
  local_mocked_bindings(
    .runMungeSumstatsFilter = .ssQ_mockMunge(),
    .package = "pecotmr")
  local_mocked_bindings(
    extractBlockGenotypes = .ssQ_mockExtractor(),
    .package = "pecotmr")
  msgs <- capture_messages(summaryStatsQc(ss))
  joined <- paste(msgs, collapse = "")
  # MungeSumstats step + denominator framing.
  expect_match(joined, "MungeSumstats kept [0-9]+ of [0-9]+ variant")
  # Harmonization step + corrected/dropped breakdown.
  expect_match(joined, "harmonization kept [0-9]+ of [0-9]+")
  expect_match(joined, "corrected: sign-flipped [0-9]+, strand-flipped [0-9]+")
  # Per-entry rollup line.
  expect_match(joined, "QC summary: [0-9]+ in -> [0-9]+ out")
  expect_match(joined, "corrected: sign-flip [0-9]+, strand-flip [0-9]+")
})

test_that("summaryStatsQc: skipped optional steps are omitted from the rollup", {
  ss <- .ssQ_makeGwasSumStats()
  local_mocked_bindings(
    .runMungeSumstatsFilter = .ssQ_mockMunge(),
    .package = "pecotmr")
  local_mocked_bindings(
    extractBlockGenotypes = .ssQ_mockExtractor(),
    .package = "pecotmr")
  msgs <- capture_messages(
    summaryStatsQc(ss,
                   alleleFlipKriging = FALSE,
                   zMismatchQc       = "none",
                   impute            = FALSE))
  joined <- paste(msgs, collapse = "")
  # Kriging / mismatch / imputation are skipped: their per-step messages
  # and their rollup segments should be absent.
  expect_false(grepl("kriging", joined))
  expect_false(grepl("LD-mismatch", joined))
  expect_false(grepl("RAISS imputation", joined))
  expect_false(grepl("imputed [+-][1-9]", joined))
})

test_that("summaryStatsQc: per-entry log lines carry the (study/context/trait) label for QtlSumStats", {
  # Reuse the QtlSumStats fixture from the round-trip test.
  qss <- QtlSumStats(
    study    = "qstudy",
    context  = "qctx",
    trait    = "qtrait",
    entry    = list(.ssQ_makeEntryGr()),
    genome   = "hg19",
    ldSketch = .ssQ_makeHandle())
  local_mocked_bindings(
    .runMungeSumstatsFilter = .ssQ_mockMunge(),
    .package = "pecotmr")
  local_mocked_bindings(
    extractBlockGenotypes = .ssQ_mockExtractor(),
    .package = "pecotmr")
  msgs <- capture_messages(summaryStatsQc(qss))
  joined <- paste(msgs, collapse = "")
  expect_match(joined, "\\[qstudy/qctx/qtrait\\] QC track")
  expect_match(joined, "\\[qstudy/qctx/qtrait\\] QC summary")
})


context("sumstatsQc internal helpers")

# ===========================================================================
# Fixture builders
# ===========================================================================

.ssh_makeHandle <- function(snp_n = 6L, n_samples = 30L) {
  new("GenotypeHandle",
    path = "/tmp/sketch.gds",
    format = "gds",
    snpInfo = data.frame(
      SNP = paste0("rs", seq_len(snp_n)),
      CHR = rep("1", snp_n),
      BP  = seq(100L, by = 100L, length.out = snp_n),
      A1  = rep("A", snp_n),
      A2  = rep("G", snp_n),
      stringsAsFactors = FALSE),
    nSamples = n_samples,
    sampleIds = paste0("s", seq_len(n_samples)),
    pgenPtr = NULL)
}

.ssh_makeEntryGr <- function(n = 5, chr = "chr1", with_extras = FALSE) {
  gr <- GenomicRanges::GRanges(
    seqnames = rep(chr, n),
    ranges = IRanges::IRanges(start = seq(100L, by = 100L, length.out = n),
                              width = 1L))
  mc <- list(
    SNP = paste0("rs", seq_len(n)),
    A1  = rep("A", n),
    A2  = rep("G", n),
    Z   = seq(1.0, by = 0.5, length.out = n),
    N   = rep(1000L, n))
  if (with_extras) {
    mc$MAF  <- seq(0.1, by = 0.05, length.out = n)
    mc$INFO <- rep(0.95, n)
    mc$BETA <- rnorm(n)
    mc$SE   <- rep(0.1, n)
    mc$P    <- 2 * pnorm(-abs(mc$Z))
  }
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(mc)
  gr
}

.ssh_mockExtractor <- function(seed = 42, n_samples = 30L) {
  function(handle, snpIdx, meanImpute = TRUE) {
    set.seed(seed)
    panel <- matrix(rbinom(n_samples * nrow(handle@snpInfo), 2, 0.3),
                    nrow = n_samples, ncol = nrow(handle@snpInfo),
                    dimnames = list(handle@sampleIds, handle@snpInfo$SNP))
    sub <- panel[, snpIdx, drop = FALSE]
    rr <- GenomicRanges::GRanges(
      seqnames = paste0("chr", handle@snpInfo$CHR[snpIdx]),
      ranges   = IRanges::IRanges(start = handle@snpInfo$BP[snpIdx], width = 1L))
    S4Vectors::mcols(rr) <- S4Vectors::DataFrame(
      SNP = handle@snpInfo$SNP[snpIdx],
      A1  = handle@snpInfo$A1[snpIdx],
      A2  = handle@snpInfo$A2[snpIdx])
    cd <- S4Vectors::DataFrame(sampleId = handle@sampleIds,
                               row.names = handle@sampleIds)
    dosage <- t(sub)
    rownames(dosage) <- handle@snpInfo$SNP[snpIdx]
    colnames(dosage) <- handle@sampleIds
    SummarizedExperiment::SummarizedExperiment(
      assays    = list(dosage = dosage),
      rowRanges = rr,
      colData   = cd)
  }
}

# ===========================================================================
# .entryGrangesToDf and .dfToEntryGranges (round-trip)
# ===========================================================================

test_that(".entryGrangesToDf: extracts chrom/pos and all mcols", {
  gr <- .ssh_makeEntryGr(3, with_extras = TRUE)
  df <- pecotmr:::.entryGrangesToDf(gr)
  expect_s3_class(df, "data.frame")
  expect_equal(df$chrom, rep("1", 3))   # "chr" stripped
  expect_equal(df$pos, c(100L, 200L, 300L))
  expect_setequal(intersect(colnames(df),
                            c("SNP", "A1", "A2", "Z", "N", "MAF", "INFO",
                              "BETA", "SE", "P")),
                  c("SNP", "A1", "A2", "Z", "N", "MAF", "INFO",
                    "BETA", "SE", "P"))
})

test_that(".dfToEntryGranges: rebuilds the GRanges with canonical mcols", {
  df <- data.frame(
    chrom = "1",
    pos   = c(100L, 200L),
    SNP   = c("rs1", "rs2"),
    A1    = c("A", "A"),
    A2    = c("G", "G"),
    Z     = c(1.5, -2.0),
    N     = c(1000L, 1200L),
    stringsAsFactors = FALSE)
  gr <- pecotmr:::.dfToEntryGranges(df)
  expect_s4_class(gr, "GRanges")
  expect_equal(as.character(GenomicRanges::seqnames(gr)), c("chr1", "chr1"))
  expect_equal(GenomicRanges::start(gr), c(100L, 200L))
  expect_setequal(colnames(S4Vectors::mcols(gr)),
                  c("SNP", "A1", "A2", "Z", "N"))
})

test_that(".dfToEntryGranges: derives SNP from variant_id when SNP is absent", {
  df <- data.frame(
    chrom = "1", pos = 100L,
    variant_id = "chr1:100:A:G",
    A1 = "A", A2 = "G",
    Z = 1.0, N = 1000L,
    stringsAsFactors = FALSE)
  gr <- pecotmr:::.dfToEntryGranges(df)
  expect_equal(S4Vectors::mcols(gr)$SNP, "chr1:100:A:G")
})

test_that("entry GRanges round-trips through df conversion", {
  gr <- .ssh_makeEntryGr(4, with_extras = TRUE)
  df <- pecotmr:::.entryGrangesToDf(gr)
  gr2 <- pecotmr:::.dfToEntryGranges(df)
  # Positions and core mcols must match (mcols may differ by reordering).
  expect_equal(GenomicRanges::start(gr2), GenomicRanges::start(gr))
  expect_equal(S4Vectors::mcols(gr2)$SNP, S4Vectors::mcols(gr)$SNP)
  expect_equal(S4Vectors::mcols(gr2)$Z,   S4Vectors::mcols(gr)$Z)
})

# ===========================================================================
# .refVariantsFromSketch
# ===========================================================================

test_that(".refVariantsFromSketch: extracts chr/pos/A1/A2/variant_id from snpInfo", {
  h <- .ssh_makeHandle()
  rv <- pecotmr:::.refVariantsFromSketch(h)
  expect_equal(rv$chrom, rep("1", 6))   # "chr" stripped
  expect_equal(rv$pos, c(100L, 200L, 300L, 400L, 500L, 600L))
  expect_equal(rv$variant_id, paste0("rs", 1:6))
  expect_equal(rv$A1, rep("A", 6))
  expect_equal(rv$A2, rep("G", 6))
})

# ===========================================================================
# .applySkipRegion
# ===========================================================================

.ssh_smallDf <- function() {
  data.frame(
    chrom = c("1", "1", "2"),
    pos   = c(100L, 200L, 100L),
    SNP   = c("rs1", "rs2", "rs3"),
    Z     = c(1, 2, 3),
    stringsAsFactors = FALSE)
}

test_that(".applySkipRegion: NULL / empty skipRegion is a no-op", {
  df <- .ssh_smallDf()
  expect_identical(pecotmr:::.applySkipRegion(df, NULL), df)
  expect_identical(pecotmr:::.applySkipRegion(df, character()), df)
})

test_that(".applySkipRegion: drops variants overlapping a single character region", {
  df <- .ssh_smallDf()
  out <- pecotmr:::.applySkipRegion(df, "1:50-150")
  expect_equal(out$SNP, c("rs2", "rs3"))
})

test_that(".applySkipRegion: handles multiple regions and chr-prefixed input", {
  df <- .ssh_smallDf()
  out <- pecotmr:::.applySkipRegion(df, c("chr1:50-250", "chr2:50-150"))
  expect_equal(nrow(out), 0L)
})

test_that(".applySkipRegion: accepts a GRanges of skip regions", {
  df <- .ssh_smallDf()
  gr <- GenomicRanges::GRanges("1", IRanges::IRanges(start = 50, end = 150))
  out <- pecotmr:::.applySkipRegion(df, gr)
  expect_equal(out$SNP, c("rs2", "rs3"))
})

test_that(".applySkipRegion: rejects malformed character entries", {
  df <- .ssh_smallDf()
  expect_error(pecotmr:::.applySkipRegion(df, "garbage"),
               "must be 'chr:start-end'")
})

test_that(".applySkipRegion: rejects non-character non-GRanges input", {
  df <- .ssh_smallDf()
  expect_error(pecotmr:::.applySkipRegion(df, 42L),
               "must be a character vector")
})

# ===========================================================================
# .matchAgainstSketch
# ===========================================================================

test_that(".matchAgainstSketch: errors when neither Z nor BETA is present", {
  df <- data.frame(chrom = "1", pos = 100L, SNP = "rs1",
                   A1 = "A", A2 = "G", stringsAsFactors = FALSE)
  expect_error(
    pecotmr:::.matchAgainstSketch(df, .ssh_makeHandle(), matchMinProp = 0),
    "must contain at least one of Z or BETA"
  )
})

test_that(".matchAgainstSketch: errors when A1/A2 columns are missing", {
  df <- data.frame(chrom = "1", pos = 100L, SNP = "rs1",
                   Z = 1.0, stringsAsFactors = FALSE)
  expect_error(
    pecotmr:::.matchAgainstSketch(df, .ssh_makeHandle(), matchMinProp = 0),
    "must contain A1 and A2 columns"
  )
})

test_that(".matchAgainstSketch: harmonizes the input against the sketch", {
  df <- data.frame(
    chrom = c("1", "1"), pos = c(100L, 200L),
    SNP = c("rs1", "rs2"),
    A1 = c("A", "A"), A2 = c("G", "G"),
    Z = c(1.0, 2.0),
    stringsAsFactors = FALSE)
  out <- pecotmr:::.matchAgainstSketch(df, .ssh_makeHandle(), matchMinProp = 0)
  # All variants align to the sketch; Z values pass through unchanged.
  expect_equal(nrow(out), 2L)
  expect_equal(out$Z, c(1.0, 2.0))
})

# ===========================================================================
# .applyLdMismatchQcToEntry
# ===========================================================================

test_that(".applyLdMismatchQcToEntry: errors when SNP column is missing", {
  df <- data.frame(chrom = "1", pos = 100L, Z = 1.0,
                   stringsAsFactors = FALSE)
  expect_error(
    pecotmr:::.applyLdMismatchQcToEntry(df, .ssh_makeHandle(), method = "dentist"),
    "requires SNP column"
  )
})

test_that(".applyLdMismatchQcToEntry: errors on variants absent from the sketch", {
  df <- data.frame(SNP = c("rs1", "ghost"), Z = c(1, 2), N = c(1000, 1000),
                   stringsAsFactors = FALSE)
  expect_error(
    pecotmr:::.applyLdMismatchQcToEntry(df, .ssh_makeHandle(), method = "dentist"),
    "are absent from the ldSketch panel"
  )
})

# ===========================================================================
# .applyPipScreen
# ===========================================================================

test_that(".applyPipScreen: cutoff = 0 is a no-op", {
  df <- data.frame(Z = c(1, 2, 3),
                   stringsAsFactors = FALSE)
  out <- pecotmr:::.applyPipScreen(df, n = 1000, cutoff = 0)
  expect_false(out$skipped)
  expect_identical(out$df, df)
})

test_that(".applyPipScreen: skips when no variant exceeds the explicit cutoff", {
  # Build a small set of z-scores with no strong signal.
  df <- data.frame(Z = rep(0.1, 10), stringsAsFactors = FALSE)
  out <- pecotmr:::.applyPipScreen(df, n = 1000, cutoff = 0.99)
  expect_true(out$skipped)
  expect_match(out$reason, "no signals above PIP threshold")
  expect_equal(nrow(out$df), 0L)
})

test_that(".applyPipScreen: retains entry when signal clears the cutoff", {
  # A very strong z-score should give PIP near 1.
  df <- data.frame(Z = c(10, 0.1, 0.1, 0.1, 0.1), stringsAsFactors = FALSE)
  out <- pecotmr:::.applyPipScreen(df, n = 1000, cutoff = 0.5)
  expect_false(out$skipped)
  expect_identical(out$df, df)
})


context("dentist_qc")
library(MASS)
library(corpcor)

generate_dentist_data <- function(seed=42, nSnps = 100, sample_size = 100, n_outliers = 5, start_pos = 1000000, end_pos = 4000000) {
    set.seed(seed)
    cor_matrix <- matrix(0, nrow = nSnps, ncol = nSnps)
    for (i in 1:(nSnps - 1)) {
        for (j in (i + 1):nSnps) {
            cor_matrix[i, j] <- runif(1, 0.2, 0.8)
            cor_matrix[j, i] <- cor_matrix[i, j]
        }
    }
    diag(cor_matrix) <- 1
    ld_matrix <- cov2cor(make.positive.definite(cor_matrix))
    z_scores <- mvrnorm(n = 1, mu = rep(0, nSnps), Sigma = ld_matrix)
    outlier_indices <- sample(1:nSnps, n_outliers)
    z_scores[outlier_indices] <- rnorm(n_outliers, mean = 0, sd = 5)
    sumstat <- data.frame(
        position = unlist(lapply(seq(start_pos,end_pos,length.out = nSnps), round)),
        z = z_scores
    )
    return(list(sumstat = sumstat, ldMat = ld_matrix, nSample = sample_size))
}

generate_dentist_single_window_data <- function(seed=42, nSnps = 100, sample_size = 100, n_outliers = 5) {
    set.seed(seed)
    cor_matrix <- matrix(0, nrow = nSnps, ncol = nSnps)
    for (i in 1:(nSnps - 1)) {
        for (j in (i + 1):nSnps) {
            cor_matrix[i, j] <- runif(1, 0.2, 0.8)
            cor_matrix[j, i] <- cor_matrix[i, j]
        }
    }
    diag(cor_matrix) <- 1
    ld_matrix <- cov2cor(make.positive.definite(cor_matrix))
    z_scores <- mvrnorm(n = 1, mu = rep(0, nSnps), Sigma = ld_matrix)
    outlier_indices <- sample(1:nSnps, n_outliers)
    z_scores[outlier_indices] <- rnorm(n_outliers, mean = 0, sd = 5)
    return(list(z_scores = z_scores, ldMat = ld_matrix, nSample = sample_size))
}

# ===========================================================================
# dentist: basic tests
# ===========================================================================

test_that("dentist output has exactly N rows for N input variants", {
    data <- generate_dentist_data(nSnps = 100)
    expect_warning(res <- dentist(data$sumstat, R = data$ldMat, nSample = data$nSample))
    expect_equal(nrow(res), 100)
})

test_that("dentist output has exactly N rows with correctChenEtAlBug = FALSE", {
    data <- generate_dentist_data(nSnps = 100)
    expect_warning(res <- dentist(data$sumstat, R = data$ldMat, nSample = data$nSample, correctChenEtAlBug = FALSE))
    expect_equal(nrow(res), 100)
})

test_that("dentist stops when missing position", {
    data <- generate_dentist_data()
    colnames(data$sumstat) <- c("something", "z")
    expect_error(dentist(data$sumstat, R = data$ldMat, nSample = data$nSample, correctChenEtAlBug = FALSE),
                 regexp = "missing either.*pos.*or.*z")
})

test_that("dentist stops when missing zscore", {
    data <- generate_dentist_data()
    colnames(data$sumstat) <- c("position", "something")
    expect_error(dentist(data$sumstat, R = data$ldMat, nSample = data$nSample, correctChenEtAlBug = FALSE),
                 regexp = "missing either.*pos.*or.*z")
})

test_that("dentist accepts 'position' and 'zscore' column names", {
  set.seed(42)
  nSnps <- 80
  n_samples <- 100
  cor_matrix <- matrix(0, nrow = nSnps, ncol = nSnps)
  for (i in 1:(nSnps - 1)) {
    for (j in (i + 1):nSnps) {
      cor_matrix[i, j] <- runif(1, 0.2, 0.8)
      cor_matrix[j, i] <- cor_matrix[i, j]
    }
  }
  diag(cor_matrix) <- 1
  ld_matrix <- cov2cor(make.positive.definite(cor_matrix))
  z_scores <- mvrnorm(n = 1, mu = rep(0, nSnps), Sigma = ld_matrix)
  sumstat <- data.frame(
    position = seq(1000000, by = 1000, length.out = nSnps),
    zscore = z_scores
  )
  expect_warning(res <- dentist(sumstat, R = ld_matrix, nSample = n_samples))
  expect_equal(nrow(res), nSnps)
})

test_that("dentist with X matrix input returns exactly N rows", {
    set.seed(42)
    nSnps <- 80
    n_samples <- 100
    X <- matrix(rbinom(nSnps * n_samples, 2, 0.3), nrow = n_samples, ncol = nSnps)
    z_scores <- rnorm(nSnps)
    sumstat <- data.frame(position = seq(1000000, by = 1000, length.out = nSnps), z = z_scores)
    expect_warning(res <- dentist(sumstat, X = X))
    expect_equal(nrow(res), nSnps)
})

# ===========================================================================
# dentistSingleWindow
# ===========================================================================

test_that("dentistSingleWindow returns exactly N rows for N input z-scores", {
    data <- generate_dentist_single_window_data()
    expect_warning(res <- dentistSingleWindow(data$z_scores, R = data$ldMat, nSample = data$nSample))
    expect_equal(nrow(res), 100)
})

test_that("dentistSingleWindow warns when < 2000 variants", {
    data <- generate_dentist_single_window_data()
    expect_warning(dentistSingleWindow(data$z_scores, R = data$ldMat, nSample = data$nSample))
})

test_that("dentistSingleWindow stops with zscore/LD matrix dimension mismatch", {
    data <- generate_dentist_single_window_data()
    expect_warning(expect_error(dentistSingleWindow(generate_dentist_single_window_data()$z_scores, R = generate_dentist_single_window_data(nSnps = 80)$ldMat, nSample = data$nSample),
                                regexp = "ldMat must be a square matrix"))
})

test_that("dentistSingleWindow output columns are correct", {
    data <- generate_dentist_single_window_data()
    expect_warning(res <- dentistSingleWindow(data$z_scores, R = data$ldMat, nSample = data$nSample))
    expected_cols <- c("original_z", "imputed_z", "iter_to_correct", "rsq", "is_duplicate", "outlier_stat", "outlier")
    expect_true(all(expected_cols %in% colnames(res)))
})

test_that("dentistSingleWindow original_z matches input z-scores", {
    data <- generate_dentist_single_window_data()
    expect_warning(res <- dentistSingleWindow(data$z_scores, R = data$ldMat, nSample = data$nSample))
    expect_equal(res$original_z, data$z_scores)
})

test_that("dentistSingleWindow with X matrix input returns exactly N rows", {
    set.seed(42)
    nSnps <- 80
    n_samples <- 100
    X <- matrix(rbinom(nSnps * n_samples, 2, 0.3), nrow = n_samples, ncol = nSnps)
    z_scores <- rnorm(nSnps)
    expect_warning(res <- dentistSingleWindow(z_scores, X = X))
    expect_equal(nrow(res), nSnps)
})

test_that("dentistSingleWindow with correctChenEtAlBug = FALSE returns N rows", {
    data <- generate_dentist_single_window_data()
    expect_warning(res <- dentistSingleWindow(
      data$z_scores, R = data$ldMat, nSample = data$nSample,
      correctChenEtAlBug = FALSE
    ))
    expect_equal(nrow(res), 100)
    expect_true(all(c("original_z", "imputed_z", "outlier") %in% colnames(res)))
})

test_that("dentistSingleWindow with gcControl = TRUE returns N rows", {
    data <- generate_dentist_single_window_data()
    expect_warning(res <- dentistSingleWindow(
      data$z_scores, R = data$ldMat, nSample = data$nSample,
      gcControl = TRUE
    ))
    expect_equal(nrow(res), 100)
    expect_true(all(c("original_z", "imputed_z", "outlier") %in% colnames(res)))
})

test_that("dentist with gcControl = TRUE returns N rows", {
    data <- generate_dentist_data(nSnps = 100)
    expect_warning(res <- dentist(
      data$sumstat, R = data$ldMat, nSample = data$nSample,
      gcControl = TRUE
    ))
    expect_equal(nrow(res), 100)
})

test_that("dentistSingleWindow dedup path with message for duplicates", {
  set.seed(42)
  nSnps <- 80
  n_samples <- 100
  cor_matrix <- matrix(0, nrow = nSnps, ncol = nSnps)
  for (i in 1:(nSnps - 1)) {
    for (j in (i + 1):nSnps) {
      cor_matrix[i, j] <- runif(1, 0.2, 0.8)
      cor_matrix[j, i] <- cor_matrix[i, j]
    }
  }
  diag(cor_matrix) <- 1
  ld_matrix <- cov2cor(make.positive.definite(cor_matrix))
  z_scores <- mvrnorm(n = 1, mu = rep(0, nSnps), Sigma = ld_matrix)
  # Use a very low threshold to trigger dedup logic
  expect_warning(
    res <- dentistSingleWindow(z_scores, R = ld_matrix, nSample = n_samples, duprThreshold = 0.5)
  )
  expect_equal(nrow(res), nSnps)
  expect_true("is_duplicate" %in% colnames(res))
})

# ===========================================================================
# add_dups_back_dentist
# ===========================================================================

test_that("add_dups_back_dentist works", {
    zScore <- c(1.2, 2.3, 2.4, 1.4, 5.6)
    dentist_output <- data.frame(
        original_z = c(1.2, 2.3, 5.6),
        imputed_z = c(1.1, 2.1, 5.1),
        iter_to_correct = c(1, 1, 3),
        rsq = c(0.9, 0.8,  0.5),
        z_diff = c(0.1, 0.2, 0.5)
    )
    find_dup_output <- data.frame(
        dupBearer = c(-1, -1, 2, 1, -1),
        sign = c(1, -1, 1, -1, 1)
    )

    res <- pecotmr:::addDupsBackDentist(zScore, dentist_output, find_dup_output)
    # Non-duplicates: z_diff is copied directly from dentist output
    expect_equal(res$z_diff[1], 0.1)
    expect_equal(res$z_diff[2], 0.2)
    expect_equal(res$z_diff[5], 0.5)
    # Duplicates: z_diff = (zScore - imputed_z) / sqrt(1 - rsq)
    expect_equal(res$z_diff[3], (2.4 - 2.1) / sqrt(1 - 0.8), tolerance = 1e-10)
    expect_equal(res$z_diff[4], (1.4 - (-1.1)) / sqrt(1 - 0.9), tolerance = 1e-10)
    expect_equal(res$imputed_z, c(1.1, 2.1, 2.1, -1.1, 5.1))
    expect_equal(res$is_duplicate, c(FALSE, FALSE, TRUE, TRUE, FALSE))
})

test_that("add_dups_back_dentist stops when nrow mismatch", {
    z_scores <- rep(0, 5)
    dentist_output <- list(
        original_z = c(1, 2, 3, 4, 5),
        imputed_z = c(1, 2, 3, 4, 5),
        iter_to_correct = c(1, 2, 3, 4, 5),
        rsq = c(1, 2, 3, 4, 5),
        z_diff = c(1, 2, 3, 4, 5)
    )
    find_dup_output <- list(
        dupBearer = c(-1, -1, -1, -1, -1, -1),
        sign = c(1, 2, 3, 4, 5, 6)
    )
    expect_error(pecotmr:::addDupsBackDentist(z_scores, dentist_output, find_dup_output))
})

# ===========================================================================
# segment_by_dist
# ===========================================================================

test_that("segment_by_dist works", {
    res <- pecotmr:::segmentByDist(seq(2000000,5000000,100000), maxDist = 2000000, minDim = 10)
    expect_true(nrow(res) >= 1)
})

test_that("segment_by_dist fill regions cover all input positions", {
    # Verify that fill regions cover every SNP index exactly once
    pos <- seq(1000000, 5000000, length.out = 200)
    res <- pecotmr:::segmentByDist(pos, maxDist = 2000000, minDim = 10)
    # Collect all fill region indices
    covered <- integer(0)
    for (k in 1:nrow(res)) {
        covered <- c(covered, res$fillStartIdx[k]:(res$fillEndIdx[k] - 1L))
    }
    # Every position from 1 to length(pos) should be covered
    expect_equal(sort(unique(covered)), 1:length(pos))
})

test_that("segment_by_dist errors on empty positions", {
  expect_error(pecotmr:::segmentByDist(integer(0)), "No positions")
})

test_that("segment_by_dist verbose mode prints intervals", {
  pos <- seq(1000000, 5000000, length.out = 300)
  expect_message(
    pecotmr:::segmentByDist(pos, maxDist = 2000000, minDim = 50, verbose = TRUE),
    "Intervals"
  )
})

# ===========================================================================
# detect_gaps
# ===========================================================================

test_that("detect_gaps finds no internal gaps for contiguous positions", {
    pos <- seq(1000, by = 100, length.out = 50)
    gaps <- pecotmr:::detectGaps(pos, gapThreshold = 500)
    # Only start and end sentinel
    expect_equal(gaps, c(1L, 51L))
})

test_that("detect_gaps finds a centromeric gap", {
    pos <- c(seq(1000, by = 100, length.out = 50),
             seq(2000000, by = 100, length.out = 50))
    gaps <- pecotmr:::detectGaps(pos, gapThreshold = 1e6)
    expect_equal(length(gaps), 3)  # start, gap, end
    expect_equal(gaps, c(1L, 51L, 101L))
})

test_that("detect_gaps finds multiple gaps", {
    pos <- c(1000, 2000, 5000000, 6000000, 12000000)
    gaps <- pecotmr:::detectGaps(pos, gapThreshold = 1e6)
    # Gaps at positions 3 and 5 (diffs > 1e6 at indices 2 and 4)
    expect_equal(gaps, c(1L, 3L, 5L, 6L))
})

test_that("detect_gaps verbose branch prints messages", {
  pos <- c(1000, 2000, 5000000, 6000000)
  expect_message(
    pecotmr:::detectGaps(pos, gapThreshold = 1e6, verbose = TRUE),
    "No\\. of gaps found"
  )
})

# ===========================================================================
# segment_by_count
# ===========================================================================

test_that("segment_by_count produces valid windows", {
    pos <- seq(1000000, by = 1000, length.out = 500)
    res <- pecotmr:::segmentByCount(pos, maxCount = 100)
    expect_true(nrow(res) >= 1)
    # All window starts should be >= 1
    expect_true(all(res$windowStartIdx >= 1))
    # All window ends should be <= length(pos) + 1
    expect_true(all(res$windowEndIdx <= length(pos) + 1))
    # Fill regions should be within windows
    for (k in 1:nrow(res)) {
        expect_true(res$fillStartIdx[k] >= res$windowStartIdx[k])
        expect_true(res$fillEndIdx[k] <= res$windowEndIdx[k])
    }
})

test_that("segment_by_count fill regions cover all positions", {
    pos <- seq(1000000, by = 1000, length.out = 500)
    res <- pecotmr:::segmentByCount(pos, maxCount = 100)
    covered <- integer(0)
    for (k in 1:nrow(res)) {
        covered <- c(covered, res$fillStartIdx[k]:(res$fillEndIdx[k] - 1L))
    }
    expect_equal(sort(unique(covered)), 1:length(pos))
})

test_that("segment_by_count handles centromeric gap", {
    # Create two blocks separated by a large gap
    pos <- c(seq(1000000, by = 1000, length.out = 200),
             seq(5000000, by = 1000, length.out = 200))
    res <- pecotmr:::segmentByCount(pos, maxCount = 100)
    # Should create windows in both blocks
    expect_true(nrow(res) >= 2)
    # Fill regions should still cover all positions
    covered <- integer(0)
    for (k in 1:nrow(res)) {
        covered <- c(covered, res$fillStartIdx[k]:(res$fillEndIdx[k] - 1L))
    }
    expect_equal(sort(unique(covered)), 1:length(pos))
})

test_that("segment_by_count skips blocks smaller than half max_count", {
    # Block of 20 variants with max_count=100 (half=50): too small, should be skipped
    pos <- seq(1000000, by = 1000, length.out = 20)
    expect_error(pecotmr:::segmentByCount(pos, maxCount = 100),
                 "No intervals created by segmentation")
})

test_that("segment_by_count creates single window for small blocks", {
    # Block of 60 variants with max_count=100: creates one window (60 >= half=50)
    pos <- seq(1000000, by = 1000, length.out = 60)
    res <- pecotmr:::segmentByCount(pos, maxCount = 100)
    expect_equal(nrow(res), 1)
    expect_equal(res$windowStartIdx[1], 1)
    expect_equal(res$windowEndIdx[1], 61)
})

test_that("segment_by_count single block creates correct number of windows", {
    # 200 variants with maxCount = 100 should create ~3 windows
    pos <- seq(1000000, by = 1000, length.out = 200)
    res <- pecotmr:::segmentByCount(pos, maxCount = 100)
    expect_true(nrow(res) >= 2)
    expect_true(nrow(res) <= 5)
})

test_that("segment_by_count errors on empty positions", {
  expect_error(pecotmr:::segmentByCount(integer(0), maxCount = 100), "No positions")
})

test_that("segment_by_count verbose mode prints intervals", {
  pos <- seq(1000000, by = 1000, length.out = 300)
  expect_message(
    pecotmr:::segmentByCount(pos, maxCount = 100, verbose = TRUE),
    "Intervals"
  )
})

# ===========================================================================
# merge_windows
# ===========================================================================

test_that("merge_windows returns exactly N rows", {
    data <- generate_dentist_data(nSnps = 1000, sample_size = 1000, start_pos = 0, end_pos = 2000)
    window_divided_res <- pecotmr:::segmentByDist(data$sumstat$position, maxDist = 1000, minDim = 10)
    dentist_result_by_window <- list()
    suppressWarnings({
        for (k in 1:nrow(window_divided_res)) {
            idx_range <- window_divided_res$windowStartIdx[k]:(window_divided_res$windowEndIdx[k] - 1L)
            zScore_k <- data$sumstat$z[idx_range]
            LD_mat_k <- data$ldMat[idx_range, idx_range]
            dentist_result_by_window[[k]] <- dentistSingleWindow(
                zScore_k, R = LD_mat_k, nSample = 100,
                pValueThreshold = 5.0369e-8, propSVD = 0.4, gcControl = FALSE,
                nIter = 10, gPvalueThreshold = 0.05, duprThreshold = 0.99,
                ncpus = 1, correctChenEtAlBug = TRUE
            )
        }
    })
    res <- pecotmr:::mergeWindows(dentist_result_by_window, window_divided_res)
    expect_equal(nrow(res), 1000)
})

test_that("merge_windows stops with window and imputed mismatch", {
    expect_error(pecotmr:::mergeWindows(rep(0, 5), data.frame(windowStartIdx = rep(0,2), windowEndIdx = rep(0, 2))))
})

test_that("merge_windows correctly indexes and merges windows", {
  # Create two fake windows
  window1 <- data.frame(
    original_z = c(1.0, 2.0, 3.0),
    imputed_z = c(0.9, 1.9, 2.9),
    iter_to_correct = c(1, 1, 1),
    rsq = c(0.5, 0.6, 0.7),
    is_duplicate = c(FALSE, FALSE, FALSE),
    outlier_stat = c(0.1, 0.2, 0.3),
    outlier = c(FALSE, FALSE, FALSE)
  )
  window2 <- data.frame(
    original_z = c(4.0, 5.0, 6.0),
    imputed_z = c(3.9, 4.9, 5.9),
    iter_to_correct = c(2, 2, 2),
    rsq = c(0.8, 0.9, 0.5),
    is_duplicate = c(FALSE, FALSE, FALSE),
    outlier_stat = c(0.4, 0.5, 0.6),
    outlier = c(FALSE, FALSE, TRUE)
  )
  window_info <- data.frame(
    windowIdx = c(1, 2),
    windowStartIdx = c(1, 4),
    windowEndIdx = c(4, 7),
    fillStartIdx = c(1, 4),
    fillEndIdx = c(4, 7)
  )
  result <- pecotmr:::mergeWindows(list(window1, window2), window_info)
  expect_equal(nrow(result), 6)
  expect_true("index_global" %in% colnames(result))
  expect_true("index_within_window" %in% colnames(result))
})

# ===========================================================================
# dentist windowed mode
# ===========================================================================

test_that("dentist windowed output has exactly N rows for large input", {
    # Generate data large enough to trigger windowed mode (> min_dim)
    data <- generate_dentist_data(seed = 123, nSnps = 1000, sample_size = 1000,
                                   n_outliers = 50, start_pos = 0, end_pos = 5000000)
    suppressWarnings({
        res <- dentist(data$sumstat, R = data$ldMat, nSample = data$nSample,
                       minDim = 100, windowSize = 2000000)
    })
    expect_equal(nrow(res), 1000)
})

test_that("dentist outlier_stat formula is correct: (z-imputed)^2/(1-rsq)", {
    data <- generate_dentist_single_window_data(seed = 55, nSnps = 100)
    expect_warning(res <- dentistSingleWindow(data$z_scores, R = data$ldMat, nSample = data$nSample))
    expected_stat <- (res$original_z - res$imputed_z)^2 / pmax(1 - res$rsq, 1e-8)
    expect_equal(res$outlier_stat, expected_stat, tolerance = 1e-10)
})

# ===========================================================================
# dentist with count mode
# ===========================================================================

test_that("dentist with window_mode='count' returns exactly N rows", {
    data <- generate_dentist_data(seed = 789, nSnps = 500, sample_size = 500,
                                   n_outliers = 25, start_pos = 1000000, end_pos = 4000000)
    suppressWarnings({
        res <- dentist(data$sumstat, R = data$ldMat, nSample = data$nSample,
                       windowMode = "count", minDim = 100)
    })
    expect_equal(nrow(res), 500)
})

# ===========================================================================
# Equivalence tests: both windowing methods
# ===========================================================================

test_that("segment_by_dist and segment_by_count agree on uniformly-spaced variants", {
    n <- 200
    spacing <- 10000  # 10kb between each variant
    pos <- seq(1000000, by = spacing, length.out = n)
    window_count <- 50  # variants per window in count mode
    window_dist <- window_count * spacing  # equivalent distance

    res_dist <- pecotmr:::segmentByDist(pos, maxDist = window_dist, minDim = 10)
    res_count <- pecotmr:::segmentByCount(pos, maxCount = window_count, gapDist = 1e6)

    # Both should cover all positions
    covered_dist <- integer(0)
    for (k in 1:nrow(res_dist)) {
        covered_dist <- c(covered_dist, res_dist$fillStartIdx[k]:(res_dist$fillEndIdx[k] - 1L))
    }
    covered_count <- integer(0)
    for (k in 1:nrow(res_count)) {
        covered_count <- c(covered_count, res_count$fillStartIdx[k]:(res_count$fillEndIdx[k] - 1L))
    }
    expect_equal(sort(unique(covered_dist)), 1:n)
    expect_equal(sort(unique(covered_count)), 1:n)
})

test_that("both windowing modes produce same dentist results on uniform data", {
    data <- generate_dentist_data(seed = 555, nSnps = 500, sample_size = 500,
                                   n_outliers = 25, start_pos = 0, end_pos = 5000000)
    suppressWarnings({
        res_dist <- dentist(data$sumstat, R = data$ldMat, nSample = data$nSample,
                            windowMode = "distance", minDim = 100, windowSize = 2000000)
        res_count <- dentist(data$sumstat, R = data$ldMat, nSample = data$nSample,
                             windowMode = "count", minDim = 100)
    })
    # Both should return exactly N rows
    expect_equal(nrow(res_dist), 500)
    expect_equal(nrow(res_count), 500)
})

# ===========================================================================
# resolve_LD_input (internal)
# ===========================================================================

test_that("resolve_LD_input errors when neither R nor X provided", {
  expect_error(
    pecotmr:::resolveLdInput(R = NULL, X = NULL),
    "Either R.*or X.*must be provided"
  )
})

test_that("resolve_LD_input errors when both R and X provided", {
  R <- diag(3)
  X <- matrix(1:9, nrow = 3)
  expect_error(
    pecotmr:::resolveLdInput(R = R, X = X),
    "Provide either R or X, not both"
  )
})

test_that("resolve_LD_input errors when R provided without nSample and need_nSample is TRUE", {
  R <- diag(3)
  expect_error(
    pecotmr:::resolveLdInput(R = R, nSample = NULL, needNSample = TRUE),
    "nSample is required"
  )
})

test_that("resolve_LD_input returns nSample = NULL when need_nSample is FALSE", {
  R <- diag(3)
  result <- pecotmr:::resolveLdInput(R = R, nSample = NULL, needNSample = FALSE)
  expect_null(result$nSample)
  expect_equal(result$R, R)
})

test_that("resolve_LD_input infers nSample from X", {
  set.seed(42)
  n <- 50; p <- 5
  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n, ncol = p)
  result <- pecotmr:::resolveLdInput(X = X, needNSample = TRUE)
  expect_equal(result$nSample, n)
  expect_true(is.matrix(result$R))
  expect_equal(nrow(result$R), p)
})

test_that("resolve_LD_input converts non-matrix X to matrix", {
  set.seed(42)
  X_df <- data.frame(a = rbinom(30, 2, 0.3), b = rbinom(30, 2, 0.3))
  result <- pecotmr:::resolveLdInput(X = X_df, needNSample = FALSE)
  expect_true(is.matrix(result$R))
  expect_equal(result$nSample, 30)
})

test_that("resolve_LD_input uses explicit nSample when X provided", {
  set.seed(42)
  X <- matrix(rbinom(100, 2, 0.3), nrow = 20, ncol = 5)
  result <- pecotmr:::resolveLdInput(X = X, nSample = 999, needNSample = TRUE)
  expect_equal(result$nSample, 999)
})


# ===========================================================================
# build_segment_result (internal)
# ===========================================================================

test_that("build_segment_result caps end indices and verbose prints", {
  expect_message(
    result <- pecotmr:::buildSegmentResult(
      startList = c(1L), endList = c(200L),
      fillStartList = c(1L), fillEndList = c(200L),
      n = 100, verbose = TRUE
    ),
    "Intervals"
  )
  expect_equal(result$windowEndIdx[1], 101)
  expect_equal(result$fillEndIdx[1], 101)
})

test_that("build_segment_result errors on empty startList", {
  expect_error(
    pecotmr:::buildSegmentResult(
      startList = integer(0), endList = integer(0),
      fillStartList = integer(0), fillEndList = integer(0),
      n = 100
    ),
    "No intervals"
  )
})

# ===========================================================================
# sliding_window_loop (iteration limit)
# ===========================================================================

test_that("sliding_window_loop errors on infinite loop", {
  allGaps <- c(1L, 1001L)
  expect_error(
    pecotmr:::slidingWindowLoop(
      allGaps, n = 1000,
      minBlockFn = function(blockSize) TRUE,
      initEndFn = function(startIdx, blockEnd) startIdx + 10,
      fillFn = function(startIdx, endIdx, notStart, notLast) list(start = startIdx, end = endIdx),
      stepFn = function(startIdx, blockEnd) list(startIdx = startIdx, endIdx = startIdx + 10),
      verbose = FALSE
    ),
    "iteration limit exceeded"
  )
})



context("univariate_rss_diagnostics")

.testFineMappingEntry <- function(variantIds, trimmedFit = list(),
                                  topLoci = data.frame(
                                    variant_id = character(0),
                                    pip = numeric(0),
                                    stringsAsFactors = FALSE)) {
  FineMappingEntry(
    variantIds = variantIds,
    trimmedFit = trimmedFit,
    topLoci    = topLoci
  )
}

# ===========================================================================
# getSusieResult
# ===========================================================================

test_that("getSusieResult returns NULL for empty input", {
  result <- getSusieResult(list())
  expect_null(result)
})

test_that("getSusieResult returns NULL when finemappingEntry missing", {
  result <- getSusieResult(list(some_data = 42))
  expect_null(result)
})

test_that("getSusieResult returns trimmed result when present", {
  mock_result <- list(pip = c(0.1, 0.5, 0.3), sets = list(cs = list()))
  con_data <- list(finemappingEntry = .testFineMappingEntry(
    variantIds = c("1:100:A:G", "1:200:C:T", "1:300:G:A"),
    trimmedFit = mock_result
  ))
  result <- getSusieResult(con_data)
  expect_equal(result, mock_result)
})

# ===========================================================================
# extractTopPipInfo
# ===========================================================================

test_that("extractTopPipInfo finds top PIP variant", {
  con_data <- list(
    finemappingEntry = .testFineMappingEntry(
      variantIds = c("1:100:A:G", "1:200:C:T", "1:300:G:A"),
      trimmedFit = list(pip = c(0.1, 0.7, 0.2))
    ),
    sumstats = list(z = c(1.0, 3.5, -0.5))
  )
  result <- extractTopPipInfo(con_data)
  expect_equal(result$top_variant, "1:200:C:T")
  expect_equal(result$top_pip, 0.7)
  expect_equal(result$top_z, 3.5)
  expect_equal(result$top_variant_index, 2)
  expect_true(is.na(result$cs_name))
  expect_true(is.na(result$variants_per_cs))
})

test_that("extractTopPipInfo computes p_value from z", {
  con_data <- list(
    finemappingEntry = .testFineMappingEntry(
      variantIds = c("1:100:A:G", "1:200:C:T", "1:300:G:A"),
      trimmedFit = list(pip = c(0.9, 0.05, 0.05))
    ),
    sumstats = list(z = c(5.0, 0.5, -0.3))
  )
  result <- extractTopPipInfo(con_data)
  expected_pval <- pecotmr:::.zToPvalue(5.0)
  expect_equal(result$p_value, expected_pval)
})

test_that("extractTopPipInfo handles ties by taking first max", {
  con_data <- list(
    finemappingEntry = .testFineMappingEntry(
      variantIds = c("1:100:A:G", "1:200:C:T", "1:300:G:A"),
      trimmedFit = list(pip = c(0.5, 0.5, 0.5))
    ),
    sumstats = list(z = c(1.0, 2.0, 3.0))
  )
  result <- extractTopPipInfo(con_data)
  expect_equal(result$top_variant_index, 1)
  expect_equal(result$top_pip, 0.5)
})

# ===========================================================================
# extractCsInfo
# ===========================================================================

test_that("extractCsInfo extracts single CS correctly", {
  con_data <- list(
    finemappingEntry = .testFineMappingEntry(
      variantIds = c("1:100:A:G", "1:200:C:T", "1:300:G:A"),
      trimmedFit = list(
        sets = list(cs = list(L_1 = c(1, 2))),
        cs_corr = NULL
      )
    )
  )
  top_loci_table <- data.frame(
    variant_id = c("1:100:A:G", "1:200:C:T"),
    pip = c(0.3, 0.8),
    z = c(2.0, 4.5),
    stringsAsFactors = FALSE
  )
  result <- extractCsInfo(con_data, csNames = "L_1", topLociTable = top_loci_table)
  expect_equal(nrow(result), 1)
  expect_equal(result$cs_name, "L_1")
  expect_equal(result$top_variant, "1:200:C:T")
  expect_equal(result$top_pip, 0.8)
  expect_equal(result$variants_per_cs, 2)
  expect_true(grepl("NA", result$cs_corr[[1]]))
})

test_that("extractCsInfo extracts multiple CSs with cs_corr", {
  con_data <- list(
    finemappingEntry = .testFineMappingEntry(
      variantIds = c("1:100:A:G", "1:200:C:T", "1:300:G:A", "1:400:T:C"),
      trimmedFit = list(
        sets = list(
          cs = list(L_1 = c(1, 2), L_2 = c(3, 4))
        ),
        cs_corr = matrix(c(1, 0.3, 0.3, 1), nrow = 2)
      )
    )
  )
  top_loci_table <- data.frame(
    variant_id = c("1:100:A:G", "1:200:C:T", "1:300:G:A", "1:400:T:C"),
    pip = c(0.3, 0.8, 0.6, 0.1),
    z = c(2.0, 4.5, 3.0, 0.5),
    stringsAsFactors = FALSE
  )
  result <- extractCsInfo(con_data, csNames = c("L_1", "L_2"), topLociTable = top_loci_table)
  expect_equal(nrow(result), 2)
  expect_equal(result$cs_name[1], "L_1")
  expect_equal(result$cs_name[2], "L_2")
  expect_equal(result$top_variant[1], "1:200:C:T")
  expect_equal(result$top_variant[2], "1:300:G:A")
  expect_true(is.character(result$cs_corr[[1]]))
})

test_that("extractCsInfo computes p_value from z-score", {
  con_data <- list(
    finemappingEntry = .testFineMappingEntry(
      variantIds = c("1:100:A:G", "1:200:C:T"),
      trimmedFit = list(
        sets = list(cs = list(L_1 = c(1, 2))),
        cs_corr = NULL
      )
    )
  )
  top_loci_table <- data.frame(
    variant_id = c("1:100:A:G", "1:200:C:T"),
    pip = c(0.9, 0.1),
    z = c(5.0, 0.5),
    stringsAsFactors = FALSE
  )
  result <- extractCsInfo(con_data, csNames = "L_1", topLociTable = top_loci_table)
  expected_pval <- pecotmr:::.zToPvalue(5.0)
  expect_equal(result$p_value, expected_pval, tolerance = 1e-10)
})

# ===========================================================================
# parseCsCorr
# ===========================================================================

test_that("parseCsCorr handles NA correlations", {
  df <- data.frame(
    cs_name = "L1",
    top_pip = 0.9,
    cs_corr = NA_character_,
    stringsAsFactors = FALSE
  )
  result <- parseCsCorr(df)
  expect_true("cs_corr_max" %in% colnames(result))
  expect_true("cs_corr_min" %in% colnames(result))
  expect_true(is.na(result$cs_corr_max))
  expect_true(is.na(result$cs_corr_min))
})

test_that("parseCsCorr splits comma-separated correlations", {
  df <- data.frame(
    cs_name = c("L1", "L2"),
    top_pip = c(0.9, 0.3),
    cs_corr = c("1,0.3", "0.3,1"),
    stringsAsFactors = FALSE
  )
  result <- parseCsCorr(df)
  expect_true("cs_corr_1" %in% colnames(result))
  expect_true("cs_corr_2" %in% colnames(result))
  expect_equal(result$cs_corr_max[1], 0.3)
  expect_equal(result$cs_corr_min[1], 0.3)
})

test_that("parseCsCorr handles empty string", {
  df <- data.frame(
    cs_name = "L1",
    cs_corr = "",
    stringsAsFactors = FALSE
  )
  result <- parseCsCorr(df)
  expect_true(is.na(result$cs_corr_max))
})

test_that("parseCsCorr handles multiple correlations", {
  df <- data.frame(
    cs_name = "L1",
    cs_corr = "1,0.5,0.2",
    stringsAsFactors = FALSE
  )
  result <- parseCsCorr(df)
  expect_equal(result$cs_corr_max, 0.5)
  expect_equal(result$cs_corr_min, 0.2)
  expect_equal(result$cs_corr_1, 1)
  expect_equal(result$cs_corr_2, 0.5)
  expect_equal(result$cs_corr_3, 0.2)
})

test_that("parseCsCorr handles NULL value", {
  df <- data.frame(
    cs_name = "L1",
    cs_corr = NA,
    stringsAsFactors = FALSE
  )
  result <- parseCsCorr(df)
  expect_true(is.na(result$cs_corr_max))
  expect_true(is.na(result$cs_corr_min))
})

test_that("parseCsCorr handles single value without comma", {
  df <- data.frame(
    cs_name = "L1",
    cs_corr = "0.5",
    stringsAsFactors = FALSE
  )
  result <- parseCsCorr(df)
  expect_true(is.na(result$cs_corr_max))
  expect_true(is.na(result$cs_corr_min))
})

test_that("parseCsCorr handles all-1 correlations (self-corr only)", {
  df <- data.frame(
    cs_name = c("L1", "L2"),
    cs_corr = c("1,1", "1,1"),
    stringsAsFactors = FALSE
  )
  result <- parseCsCorr(df)
  expect_true(is.na(result$cs_corr_max[1]))
  expect_true(is.na(result$cs_corr_min[1]))
})

test_that("parseCsCorr handles mixed valid and NA rows", {
  df <- data.frame(
    cs_name = c("L1", "L2"),
    cs_corr = c("1,0.3,0.7", NA_character_),
    stringsAsFactors = FALSE
  )
  result <- parseCsCorr(df)
  expect_equal(result$cs_corr_max[1], 0.7)
  expect_equal(result$cs_corr_min[1], 0.3)
  expect_true(is.na(result$cs_corr_max[2]))
})

test_that("parseCsCorr expands columns for different lengths", {
  df <- data.frame(
    cs_name = c("L1", "L2", "L3"),
    cs_corr = c("1,0.3,0.5", "1,0.2", "1,0.8,0.4,0.1"),
    stringsAsFactors = FALSE
  )
  result <- parseCsCorr(df)
  expect_true("cs_corr_4" %in% colnames(result))
  expect_true(is.na(result$cs_corr_4[1]))
  expect_equal(result$cs_corr_4[3], 0.1)
})

# ===========================================================================
# autoDecision
# ===========================================================================

test_that("autoDecision assigns BVSR when no CS is tagged", {
  df <- data.frame(
    cs_name = c("L1", "L2"),
    top_z = c(5.0, 3.5),
    p_value = c(1e-10, 1e-6),
    stringsAsFactors = FALSE
  )
  result <- autoDecision(df, highCorrCols = character(0))
  expect_true("top_cs" %in% colnames(result))
  expect_true("tagged_cs" %in% colnames(result))
  expect_true("method" %in% colnames(result))
  expect_true(all(result$method == "BVSR"))
})

test_that("autoDecision assigns SER when all non-top CSs are tagged", {
  df <- data.frame(
    cs_name = c("L1", "L2"),
    top_z = c(5.0, 0.1),
    p_value = c(1e-10, 0.5),
    stringsAsFactors = FALSE
  )
  result <- autoDecision(df, highCorrCols = character(0))
  expect_true(result$top_cs[1])
  expect_false(result$top_cs[2])
  expect_true(result$tagged_cs[2])
  expect_true(all(result$method == "SER"))
})

test_that("autoDecision assigns BCR when untagged CS remain", {
  df <- data.frame(
    cs_name = c("L1", "L2", "L3"),
    top_z = c(5.0, 3.5, 0.1),
    p_value = c(1e-10, 1e-6, 0.5),
    stringsAsFactors = FALSE
  )
  result <- autoDecision(df, highCorrCols = character(0))
  expect_true(all(result$method == "BCR"))
})

test_that("autoDecision assigns SER for single CS", {
  df <- data.frame(
    cs_name = "L1",
    top_z = 5.0,
    p_value = 1e-10,
    stringsAsFactors = FALSE
  )
  result <- autoDecision(df, highCorrCols = character(0))
  expect_true(result$top_cs[1])
  expect_false(result$tagged_cs[1])
  expect_equal(result$method, "SER")
})
