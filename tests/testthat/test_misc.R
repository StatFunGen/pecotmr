context("misc")
library(tidyverse)

# =============================================================================
# compute_maf
# =============================================================================

test_that("Test compute_maf freq 0.5",{
    expect_equal(compute_maf(rep(1, 20)), 0.5)
})

test_that("Test compute_maf freq 0.6",{
    expect_equal(compute_maf(rep(1.2, 20)), 0.4)
})

test_that("Test compute_maf freq 0.3",{
    expect_equal(compute_maf(rep(0.6, 20)), 0.3)
})

test_that("Test compute_maf with NA",{
    set.seed(1)
    generate_small_dataset <- function(sample_size = 20) {
        vals <- c(1.2, NA)
        return(sample(vals, sample_size, replace = TRUE))
    }
    expect_equal(compute_maf(generate_small_dataset()), 0.4)
})

test_that("compute_maf returns 0 for monomorphic (all 0)", {
  expect_equal(compute_maf(rep(0, 10)), 0)
})

# =============================================================================
# compute_missing
# =============================================================================

test_that("test compute_missing",{
    small_dataset <- c(rep(NA, 20), rep(1, 80))
    expect_equal(compute_missing(small_dataset), 0.2)
})

# =============================================================================
# compute_non_missing_y and compute_all_missing_y
# =============================================================================

test_that("Test compute_non_missing_y",{
    small_dataset <- c(rep(NA, 20), rep(1, 80))
    expect_equal(compute_non_missing_y(small_dataset), 80)
})

test_that("Test compute_all_missing_y",{
    small_dataset <- c(rep(NA, 20), rep(1, 80))
    expect_equal(compute_all_missing_y(small_dataset), F)
})

test_that("compute_all_missing_y returns TRUE for all-NA vector", {
  expect_true(pecotmr:::compute_all_missing_y(rep(NA, 5)))
})

test_that("compute_all_missing_y returns FALSE for partially NA vector", {
  expect_false(pecotmr:::compute_all_missing_y(c(NA, 1, NA)))
})

# =============================================================================
# mean_impute
# =============================================================================

test_that("Test mean_impute",{
    dummy_data <- matrix(c(1,2,NA,1,2,3), nrow=3, ncol=2)
    expect_equal(mean_impute(dummy_data)[3,1], 1.5)
})

test_that("mean_impute with all NAs in a column imputes NaN", {
  X <- matrix(c(NA, NA, NA, 1, 2, 3), nrow = 3, ncol = 2)
  result <- pecotmr:::mean_impute(X)
  expect_true(all(is.nan(result[, 1])))
  expect_equal(result[, 2], c(1, 2, 3))
})

# =============================================================================
# is_zero_variance
# =============================================================================

test_that("Test is_zero_variance",{
    dummy_data <- matrix(c(1,2,3,1,1,1), nrow=3, ncol=2)
    col <- which(apply(dummy_data, 2, is_zero_variance))
    expect_equal(col, 2)
})

test_that("is_zero_variance with NA values treats them as distinct", {
  expect_false(pecotmr:::is_zero_variance(c(1, NA, 1)))
})

# =============================================================================
# filter_X
# =============================================================================

test_that("Test filter_X",{
    dummy_data <- matrix(
        c(1,NA,NA,NA, 0,0,1,1, 2,2,2,2, 1,1,1,2, 2,2,0,1, 0,1,1,2),
        # Missing Rate, MAF thresh, Zero Var, Var Thresh, Regular values
        nrow=4, ncol=6)
    var_thres <- 0.3
    expect_equal(filter_X(dummy_data, 0.70, 0.3, var_thres = 0.3), matrix(c(2,2,0,1, 0,1,1,2), nrow=4, ncol=2))
})

test_that("filter_X drops most columns when nearly all are zero variance", {
  X <- matrix(c(
    1, 1, 1, 1, 1,
    2, 2, 2, 2, 2,
    0, 0, 0, 0, 0,
    0, 1, 2, 0, 1
  ), nrow = 5, ncol = 4)
  result <- pecotmr:::filter_X(X, missing_rate_thresh = 1.0, maf_thresh = 0)
  expect_equal(ncol(result), 1)
})

test_that("filter_X with external maf vector uses it for filtering", {
  set.seed(42)
  X <- matrix(sample(0:2, 40, replace = TRUE), nrow = 10, ncol = 4)
  external_maf <- c(0.01, 0.05, 0.3, 0.4)
  result <- pecotmr:::filter_X(X, missing_rate_thresh = 1.0, maf_thresh = 0.1, maf = external_maf)
  # Columns 1 and 2 have MAF <= 0.1, so they are dropped; columns 3 and 4 remain
  expect_equal(ncol(result), 2)
})

test_that("filter_X skips MAF filtering for non-0/1/2 genotypes without external MAF", {
  set.seed(42)
  X <- matrix(runif(40, 0, 2), nrow = 10, ncol = 4)
  expect_message(
    result <- pecotmr:::filter_X(X, missing_rate_thresh = 1.0, maf_thresh = 0.1),
    "Skipping MAF filtering"
  )
  expect_true(ncol(result) >= 1)
})

test_that("filter_X applies var_thresh with external X_variance", {
  set.seed(42)
  X <- matrix(sample(0:2, 40, replace = TRUE), nrow = 10, ncol = 4)
  external_var <- c(0.01, 0.5, 1.0, 0.02)
  result <- pecotmr:::filter_X(X, missing_rate_thresh = 1.0, maf_thresh = 0, var_thresh = 0.1, X_variance = external_var)
  # Columns 1 and 4 have variance < 0.1, so they are dropped; columns 2 and 3 remain
  expect_equal(ncol(result), 2)
})

test_that("filter_X with NULL thresholds does not filter", {
  set.seed(42)
  X <- matrix(sample(0:2, 40, replace = TRUE), nrow = 10, ncol = 4)
  result <- pecotmr:::filter_X(X, missing_rate_thresh = NULL, maf_thresh = NULL, var_thresh = 0)
  # No filtering applied: all 4 columns should remain (zero-variance check still runs but none are zero-variance)
  expect_equal(ncol(result), 4)
})

test_that("filter_X with missing_rate_thresh=0 drops columns with any NA", {
  X <- matrix(c(0, 1, 2, 1, 0, 1, NA, 2, 0, 1, 2, 1), nrow = 4, ncol = 3)
  result <- pecotmr:::filter_X(X, missing_rate_thresh = 0, maf_thresh = 0, var_thresh = 0)
  expect_true(ncol(result) <= 2)
})

test_that("filter_X with maf_thresh=0.5 removes nearly all columns", {
  X <- matrix(c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 2,
    0, 0, 0, 0, 0, 0, 0, 0, 2, 2,
    0, 0, 0, 0, 0, 2, 2, 2, 2, 2
  ), nrow = 10, ncol = 3)
  result <- pecotmr:::filter_X(X, missing_rate_thresh = 1.0, maf_thresh = 0.45, var_thresh = 0)
  expect_equal(ncol(result), 1)
})

# =============================================================================
# filter_Y
# =============================================================================

test_that("Test filter_Y non-matrix",{
    dummy_data <- matrix(c(1,NA,NA,NA, 1,1,2,NA), nrow=4, ncol=2)
    res <- filter_Y(as.data.frame(dummy_data), 3)
    expect_equal(length(res$Y), 3)
    expect_equal(res$rm_rows, NULL)
})

test_that("Test filter_Y is-matrix",{
    dummy_data <- matrix(c(1,NA,NA,NA, 1,1,2,NA, 2,1,2,NA), nrow=4, ncol=3)
    expect_equal(nrow(filter_Y(dummy_data, 3)$Y), 3)
    expect_equal(ncol(filter_Y(dummy_data, 3)$Y), 2)
    expect_equal(length(filter_Y(dummy_data, 3)$rm_rows), 1)
})

test_that("filter_Y removes all columns with insufficient observations", {
  Y <- matrix(c(NA, NA, NA, 1, NA, NA, NA, 2), nrow = 4, ncol = 2)
  result <- pecotmr:::filter_Y(Y, n_nonmiss = 3)
  expect_true(length(result$Y) == 0 || ncol(as.matrix(result$Y)) == 0)
})

test_that("filter_Y removes columns with too few non-missing values", {
  Y <- matrix(c(1, NA, NA, NA, 1, 2, 3, 4), nrow = 4)
  result <- pecotmr:::filter_Y(Y, n_nonmiss = 3)
  expect_true(length(result$Y) >= 3)
})

test_that("filter_Y removes all-NA rows from matrix", {
  Y <- matrix(c(NA, NA, 1, 2, NA, NA, 3, 4), nrow = 4)
  result <- pecotmr:::filter_Y(Y, n_nonmiss = 1)
  expect_true(nrow(result$Y) < 4)
})

# =============================================================================
# format_variant_id
# =============================================================================

test_that("Test format_variant_id",{
    expect_equal(format_variant_id(c(1, 1), c(123, 132), c("G", "A"), c("C", "T")), c("chr1:123:G:C", "chr1:132:A:T"))
})

test_that("format_variant_id uses convention parameter automatically", {
  expect_equal(pecotmr:::format_variant_id(1, 100, "A", "G"), "chr1:100:A:G")

  conv_mixed <- list(has_chr = TRUE, allele_sep = "_")
  expect_equal(pecotmr:::format_variant_id(1, 100, "A", "G", convention = conv_mixed), "chr1:100_A_G")

  conv_nochr <- list(has_chr = FALSE, allele_sep = "_")
  expect_equal(pecotmr:::format_variant_id(1, 100, "A", "G", convention = conv_nochr), "1:100_A_G")

  expect_equal(pecotmr:::format_variant_id(1, 100, "A", "G", chr_prefix = FALSE, convention = conv_mixed), "chr1:100_A_G")
})

test_that("format_variant_id constructs canonical IDs", {
  expect_equal(pecotmr:::format_variant_id(c(1, 2), c(100, 200), c("A", "C"), c("G", "T")),
               c("chr1:100:A:G", "chr2:200:C:T"))
  expect_equal(pecotmr:::format_variant_id(1, 100, "A", "G", chr_prefix = FALSE), "1:100:A:G")
  expect_equal(pecotmr:::format_variant_id("chr1", 100, "A", "G"), "chr1:100:A:G")
})

# =============================================================================
# find_duplicate_variants
# =============================================================================

z <- c(1, 2, 3, 4, 5)
LD <- matrix(c(1.0, 0.8, 0.2, 0.1, 0.3,
               0.8, 1.0, 0.4, 0.2, 0.5,
               0.2, 0.4, 1.0, 0.6, 0.1,
               0.1, 0.2, 0.6, 1.0, 0.3,
               0.3, 0.5, 0.1, 0.3, 1.0), nrow = 5, ncol = 5)

test_that("find_duplicate_variants returns the expected output", {
  rThreshold <- 0.5
  expected_output <- list(
    filteredZ = c(1, 3, 5),
    filteredLD= LD[c(1,3,5), c(1,3,5)],
    dupBearer = c(-1, 1, -1, 2, -1),
    corABS = c(0, 0.8, 0, 0.6, 0),
    sign = c(1, 1, 1, 1, 1),
    minValue = 0.1
  )

  result <- find_duplicate_variants(z, LD, rThreshold)
  expect_equal(result, expected_output)
})

test_that("find_duplicate_variants handles a high correlation threshold", {
  rThreshold <- 1.0
  expected_output <- list(
    filteredZ = c(1, 2, 3, 4, 5),
    filteredLD=LD,
    dupBearer = c(-1, -1, -1, -1, -1),
    corABS = c(0, 0, 0, 0, 0),
    sign = c(1, 1, 1, 1, 1),
    minValue = 0.1
  )

  result <- find_duplicate_variants(z, LD, rThreshold)
  expect_equal(result, expected_output)
})

test_that("find_duplicate_variants handles a low correlation threshold", {
  rThreshold <- 0.0
  expected_output <- list(
    filteredZ = c(1),
    filteredLD=LD[1,1,drop=F],
    dupBearer = c(-1, 1, 1, 1, 1),
    corABS = c(0, 0.8, 0.2, 0.1, 0.3),
    sign = c(1, 1, 1, 1, 1),
    minValue = 0.1
  )

  result <- find_duplicate_variants(z, LD, rThreshold)
  expect_equal(result, expected_output)
})

test_that("find_duplicate_variants handles negative correlations", {
  LD_negative <- LD
  LD_negative[1, 2] <- -0.8
  LD_negative[2, 1] <- -0.8
  rThreshold <- 0.5
  expected_output <- list(
    filteredZ = c(1, 3, 5),
    filteredLD= LD[c(1,3,5), c(1,3,5)],
    dupBearer = c(-1, 1, -1, 2, -1),
    corABS = c(0, 0.8, 0, 0.6, 0),
    sign = c(1, -1, 1, 1, 1),
    minValue = 0.1
  )

  result <- find_duplicate_variants(z, LD_negative, rThreshold)
  expect_equal(result, expected_output)
})

# =============================================================================
# pval_global
# =============================================================================

test_that("pval_global with ACAT method returns valid combined p-value", {
  pvals <- c(0.01, 0.05, 0.5, 0.8)
  result <- pecotmr:::pval_global(pvals, comb_method = "ACAT", naive = FALSE)
  expect_true(is.numeric(result))
  # ACAT statistic: T = mean(tan(pi*(0.5 - p_i))), p = P[Cauchy >= T]
  expected <- pcauchy(mean(tan(pi * (0.5 - pvals))), lower.tail = FALSE)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("pval_global with naive=TRUE returns Bonferroni-corrected p-value", {
  pvals <- c(0.01, 0.05, 0.1, 0.5)
  result <- pecotmr:::pval_global(pvals, comb_method = "HMP", naive = TRUE)
  n_unique <- length(unique(pvals))
  expected <- min(n_unique * min(pvals), 1.0)
  expect_equal(result, expected)
})

test_that("pval_global naive method caps at 1.0", {
  pvals <- seq(0.1, 0.9, by = 0.01)
  result <- pecotmr:::pval_global(pvals, comb_method = "ACAT", naive = TRUE)
  expect_true(result <= 1.0)
})

test_that("pval_global naive method with single p-value returns that p-value", {
  result <- pecotmr:::pval_global(0.03, comb_method = "ACAT", naive = TRUE)
  expect_equal(result, 0.03)
})

test_that("pval_global ACAT with identical p-values", {
  pvals <- rep(0.05, 5)
  result <- pecotmr:::pval_global(pvals, comb_method = "ACAT", naive = FALSE)
  expect_true(is.numeric(result))
  expect_true(result > 0 && result < 1)
})

test_that("pval_global ACAT with single p-value delegates correctly", {
  result <- pecotmr:::pval_global(0.05, comb_method = "ACAT", naive = FALSE)
  expect_equal(result, 0.05)
})

test_that("pval_global ACAT with very significant p-values", {
  pvals <- c(1e-8, 1e-6, 1e-4)
  result <- pecotmr:::pval_global(pvals, comb_method = "ACAT", naive = FALSE)
  expect_true(is.numeric(result))
  expect_true(result > 0 && result <= 1)
})

test_that("pval_global HMP method returns valid p-value when harmonicmeanp available", {
  skip_if_not_installed("harmonicmeanp")
  pvals <- c(0.01, 0.05, 0.2, 0.7)
  result <- pecotmr:::pval_global(pvals, comb_method = "HMP", naive = FALSE)
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("pval_global HMP errors when harmonicmeanp not installed", {
  skip_if(requireNamespace("harmonicmeanp", quietly = TRUE),
          "harmonicmeanp is installed, cannot test missing-package path")
  pvals <- c(0.01, 0.05)
  expect_error(pecotmr:::pval_global(pvals, comb_method = "HMP", naive = FALSE),
               "harmonicmeanp")
})

# =============================================================================
# pval_hmp
# =============================================================================

test_that("pval_hmp returns valid p-value", {
  skip_if_not_installed("harmonicmeanp")
  pvals <- c(0.01, 0.05, 0.1)
  result <- pecotmr:::pval_hmp(pvals)
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
  # The harmonic mean is L/sum(1/p) where L = length(unique(pvals))
  L <- length(unique(pvals))
  HMP <- L / sum(1 / pvals)
  # Result should be based on pLandau(1/HMP, ...) and be smaller than the arithmetic mean
  expect_true(result < mean(pvals))
  # Verify the result is less than the smallest individual p-value is not required,
  # but it should be in a reasonable range relative to the harmonic mean
  expect_true(result < 0.1)
})

test_that("pval_hmp uses unique p-values only", {
  skip_if_not_installed("harmonicmeanp")
  pvals <- c(0.01, 0.01, 0.05, 0.05, 0.3)
  result <- pecotmr:::pval_hmp(pvals)
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("pval_hmp errors when package not available", {
  skip_if(requireNamespace("harmonicmeanp", quietly = TRUE),
          "harmonicmeanp is installed, cannot test missing-package path")
  expect_error(pecotmr:::pval_hmp(c(0.01, 0.05)), "harmonicmeanp")
})

# =============================================================================
# pval_acat
# =============================================================================

test_that("pval_acat returns single p-value unchanged", {
  expect_equal(pecotmr:::pval_acat(0.05), 0.05)
})

test_that("pval_acat combines multiple p-values", {
  pvals <- c(0.001, 0.01, 0.1)
  combined <- pecotmr:::pval_acat(pvals)
  expect_true(combined > 0 && combined < 1)
})

test_that("pval_acat with very small p-values does not return NA", {
  result <- pecotmr:::pval_acat(c(1e-10, 1e-8, 1e-6))
  expect_true(is.numeric(result))
  expect_true(!is.na(result))
  expect_true(result > 0 && result <= 1)
})

test_that("pval_acat with all large p-values returns large combined p-value", {
  result <- pecotmr:::pval_acat(c(0.8, 0.9, 0.95))
  expect_true(is.numeric(result))
  expect_true(result > 0.5 && result <= 1)
})

# =============================================================================
# pval_cauchy
# =============================================================================

test_that("pval_cauchy combines p-values", {
  pvals <- c(0.01, 0.05, 0.5)
  combined <- pecotmr:::pval_cauchy(pvals)
  # Manual: CCT stat = mean(tan((0.5 - p) * pi)), result = 1 - pcauchy(stat)
  cct_stat <- mean(tan((0.5 - pvals) * pi))
  expected <- 1 - pcauchy(cct_stat)
  expect_equal(combined, expected, tolerance = 1e-10)
})

test_that("pval_cauchy handles NAs with na.rm", {
  pvals <- c(0.01, NA, 0.05)
  combined <- pecotmr:::pval_cauchy(pvals, na.rm = TRUE)
  expect_true(!is.na(combined))
})

test_that("pval_cauchy with very small p-values", {
  pvals <- c(1e-20, 1e-15)
  combined <- pecotmr:::pval_cauchy(pvals)
  expect_true(combined < 1e-10)
})

test_that("pval_cauchy with all NA and na.rm=TRUE returns NA", {
  result <- pecotmr:::pval_cauchy(c(NA, NA, NA), na.rm = TRUE)
  expect_true(is.na(result))
})

test_that("pval_cauchy with p-values near 1 caps them at 0.99", {
  result <- pecotmr:::pval_cauchy(c(0.999, 0.9999, 0.5))
  expect_true(is.numeric(result))
  expect_true(!is.na(result))
})

test_that("pval_cauchy with na.rm=FALSE and NA present still computes result", {
  result <- pecotmr:::pval_cauchy(c(0.01, NA, 0.05), na.rm = FALSE)
  expect_true(is.numeric(result))
})

# =============================================================================
# compute_qvalues
# =============================================================================

test_that("compute_qvalues returns NA vector when all pvalues are NA", {
  result <- pecotmr:::compute_qvalues(rep(NA_real_, 5))
  expect_true(all(is.na(result)))
  expect_length(result, 5)
})

test_that("compute_qvalues returns single p-value unchanged", {
  result <- pecotmr:::compute_qvalues(0.05)
  expect_equal(result, 0.05)
})

test_that("compute_qvalues returns empty input unchanged", {
  result <- pecotmr:::compute_qvalues(numeric(0))
  expect_length(result, 0)
})

test_that("compute_qvalues works with valid p-value vector", {
  skip_if_not_installed("qvalue")
  set.seed(42)
  pvals <- runif(100, 0, 1)
  result <- pecotmr:::compute_qvalues(pvals)
  expect_length(result, 100)
  expect_true(all(result >= 0 & result <= 1))
})

test_that("compute_qvalues falls back to BH when too few p-values", {
  skip_if_not_installed("qvalue")
  pvals <- c(0.001, 0.999)
  result <- pecotmr:::compute_qvalues(pvals)
  expect_length(result, 2)
  expect_true(all(result >= 0 & result <= 1))
})

test_that("compute_qvalues errors when qvalue not installed", {
  skip_if(requireNamespace("qvalue", quietly = TRUE),
          "qvalue is installed, cannot test missing-package path")
  expect_error(pecotmr:::compute_qvalues(c(0.01, 0.05)), "qvalue")
})

# =============================================================================
# filter_molecular_events
# =============================================================================

test_that("filter_molecular_events errors when filter lacks required fields", {
  events <- c("gene_A_splicing", "gene_B_expression")
  bad_filter <- list(list(type_pattern = "gene"))
  expect_error(
    pecotmr:::filter_molecular_events(events, bad_filter),
    "type_pattern and at least one of"
  )
})

test_that("filter_molecular_events errors when type_pattern is NULL", {
  events <- c("gene_A_splicing")
  bad_filter <- list(list(type_pattern = NULL, valid_pattern = "splicing"))
  expect_error(
    pecotmr:::filter_molecular_events(events, bad_filter),
    "type_pattern and at least one of"
  )
})

test_that("filter_molecular_events keeps events matching valid_pattern", {
  events <- c("gene_A_splicing", "gene_A_expression", "gene_B_splicing", "protein_X")
  filters <- list(list(
    type_pattern = "gene_",
    valid_pattern = "splicing"
  ))
  result <- pecotmr:::filter_molecular_events(events, filters)
  expect_true("protein_X" %in% result)
  gene_events <- result[grepl("gene_", result)]
  expect_true(all(grepl("splicing", gene_events)))
})

test_that("filter_molecular_events excludes events matching exclude_pattern", {
  events <- c("gene_A_splicing", "gene_A_expression", "gene_B_splicing")
  filters <- list(list(
    type_pattern = "gene_",
    exclude_pattern = "expression"
  ))
  result <- pecotmr:::filter_molecular_events(events, filters)
  expect_false("gene_A_expression" %in% result)
  expect_true("gene_A_splicing" %in% result)
  expect_true("gene_B_splicing" %in% result)
})

test_that("filter_molecular_events returns NULL when no events pass filtering", {
  events <- c("gene_A_expression", "gene_B_expression")
  filters <- list(list(
    type_pattern = "gene_",
    valid_pattern = "splicing"
  ))
  result <- pecotmr:::filter_molecular_events(events, filters)
  expect_null(result)
})

test_that("filter_molecular_events skips filter when no type events match", {
  events <- c("protein_X", "protein_Y")
  filters <- list(list(
    type_pattern = "gene_",
    valid_pattern = "splicing"
  ))
  result <- pecotmr:::filter_molecular_events(events, filters)
  expect_equal(sort(result), sort(events))
})

test_that("filter_molecular_events handles comma-separated valid_pattern", {
  events <- c("gene_A_splicing", "gene_A_expression", "gene_B_methylation")
  filters <- list(list(
    type_pattern = "gene_",
    valid_pattern = "splicing,expression"
  ))
  result <- pecotmr:::filter_molecular_events(events, filters)
  expect_true("gene_A_splicing" %in% result)
  expect_true("gene_A_expression" %in% result)
  expect_false("gene_B_methylation" %in% result)
})

test_that("filter_molecular_events handles comma-separated exclude_pattern", {
  events <- c("gene_A_splicing", "gene_A_expression", "gene_B_methylation")
  filters <- list(list(
    type_pattern = "gene_",
    exclude_pattern = "expression,methylation"
  ))
  result <- pecotmr:::filter_molecular_events(events, filters)
  expect_true("gene_A_splicing" %in% result)
  expect_false("gene_A_expression" %in% result)
  expect_false("gene_B_methylation" %in% result)
})

test_that("filter_molecular_events with both valid_pattern and exclude_pattern", {
  events <- c("gene_A_splicing_good", "gene_A_splicing_bad",
              "gene_B_expression", "protein_X")
  filters <- list(list(
    type_pattern = "gene_",
    valid_pattern = "splicing",
    exclude_pattern = "bad"
  ))
  result <- pecotmr:::filter_molecular_events(events, filters)
  expect_true("gene_A_splicing_good" %in% result)
  expect_false("gene_A_splicing_bad" %in% result)
  expect_true("protein_X" %in% result)
})

test_that("filter_molecular_events with condition parameter", {
  events <- c("gene_A_splicing", "gene_A_expression")
  filters <- list(list(
    type_pattern = "gene_",
    exclude_pattern = "expression"
  ))
  result <- pecotmr:::filter_molecular_events(events, filters, condition = "test_context")
  expect_true("gene_A_splicing" %in% result)
  expect_false("gene_A_expression" %in% result)
})

test_that("filter_molecular_events with multiple filters", {
  events <- c("gene_A_splicing", "gene_A_expression",
              "protein_X_high", "protein_X_low")
  filters <- list(
    list(type_pattern = "gene_", valid_pattern = "splicing"),
    list(type_pattern = "protein_", exclude_pattern = "low")
  )
  result <- pecotmr:::filter_molecular_events(events, filters)
  expect_true("gene_A_splicing" %in% result)
  expect_false("gene_A_expression" %in% result)
  expect_true("protein_X_high" %in% result)
  expect_false("protein_X_low" %in% result)
})

test_that("filter_molecular_events returns all when all events match", {
  events <- c("gene_A_splicing", "gene_B_splicing")
  filters <- list(list(
    type_pattern = "gene_",
    valid_pattern = "splicing"
  ))
  result <- pecotmr:::filter_molecular_events(events, filters)
  expect_equal(sort(result), sort(events))
})

# =============================================================================
# find_valid_file_path and find_valid_file_paths
# =============================================================================

test_that("find_valid_file_path returns target when it exists directly", {
  pkg_root <- normalizePath(file.path(test_path(), "..", ".."), mustWork = TRUE)
  target <- file.path(pkg_root, "DESCRIPTION")
  ref <- file.path(pkg_root, "NAMESPACE")
  skip_if_not(file.exists(target) && file.exists(ref), "Package root files not found")
  result <- pecotmr:::find_valid_file_path(ref, target)
  expect_equal(result, target)
})

test_that("find_valid_file_path constructs path from reference directory", {
  pkg_root <- normalizePath(file.path(test_path(), "..", ".."), mustWork = TRUE)
  ref <- file.path(pkg_root, "NAMESPACE")
  skip_if_not(file.exists(ref), "NAMESPACE not found")
  result <- pecotmr:::find_valid_file_path(ref, "DESCRIPTION")
  expect_true(file.exists(result))
  expect_true(grepl("DESCRIPTION$", result))
})

test_that("find_valid_file_path errors when both paths are invalid", {
  expect_error(
    pecotmr:::find_valid_file_path("/nonexistent/dir/ref.txt", "/nonexistent/target.txt"),
    "Both reference and target file paths do not work"
  )
})

test_that("find_valid_file_path returns reference when target is invalid but reference exists", {
  pkg_root <- normalizePath(file.path(test_path(), "..", ".."), mustWork = TRUE)
  ref <- file.path(pkg_root, "DESCRIPTION")
  skip_if_not(file.exists(ref), "DESCRIPTION not found")
  result <- pecotmr:::find_valid_file_path(ref, "/totally/bogus/path.txt")
  expect_equal(result, ref)
})

test_that("find_valid_file_paths resolves multiple targets", {
  pkg_root <- normalizePath(file.path(test_path(), "..", ".."), mustWork = TRUE)
  ref <- file.path(pkg_root, "NAMESPACE")
  skip_if_not(file.exists(ref), "NAMESPACE not found")
  targets <- c("DESCRIPTION", "NAMESPACE")
  result <- pecotmr:::find_valid_file_paths(ref, targets)
  expect_length(result, 2)
  expect_true(all(file.exists(result)))
})

test_that("find_valid_file_paths errors on all-invalid targets", {
  ref <- "/nonexistent/ref.txt"
  targets <- c("/bogus/a.txt", "/bogus/b.txt")
  expect_error(pecotmr:::find_valid_file_paths(ref, targets))
})

# =============================================================================
# compute_LD
# =============================================================================

test_that("compute_LD sample method produces valid correlation matrix", {
  set.seed(42)
  X <- matrix(sample(0:2, 200, replace = TRUE), nrow = 50)
  colnames(X) <- paste0("rs", 1:4)

  R <- compute_LD(X, method = "sample")
  expect_equal(nrow(R), 4)
  expect_equal(ncol(R), 4)
  expect_equal(unname(diag(R)), rep(1, 4))
  expect_true(isSymmetric(R))
  expect_true(all(R >= -1 & R <= 1))
})

test_that("compute_LD population method produces valid matrix", {
  set.seed(42)
  X <- matrix(sample(0:2, 200, replace = TRUE), nrow = 50)
  colnames(X) <- paste0("rs", 1:4)

  R <- compute_LD(X, method = "population")
  expect_equal(nrow(R), 4)
  expect_equal(unname(diag(R)), rep(1, 4))
  expect_true(isSymmetric(R))
})

test_that("compute_LD with a single SNP returns 1x1 identity matrix", {
  X <- matrix(c(0, 1, 2, 1, 0), ncol = 1)
  colnames(X) <- "rs1"
  R <- compute_LD(X, method = "sample")
  expect_equal(dim(R), c(1L, 1L))
  expect_equal(R[1, 1], 1.0)
  expect_equal(colnames(R), "rs1")
})

test_that("compute_LD handles column with all NA gracefully", {
  set.seed(123)
  X <- matrix(sample(0:2, 100, replace = TRUE), nrow = 20, ncol = 5)
  X[, 3] <- NA
  colnames(X) <- paste0("rs", 1:5)

  R <- compute_LD(X, method = "sample")
  expect_equal(dim(R), c(5L, 5L))
  expect_equal(unname(diag(R)), rep(1, 5))
  expect_equal(R[3, 1], 0)
  expect_equal(R[1, 3], 0)
})

test_that("compute_LD population method handles column with all NA", {
  set.seed(123)
  X <- matrix(sample(0:2, 100, replace = TRUE), nrow = 20, ncol = 5)
  X[, 2] <- NA
  colnames(X) <- paste0("rs", 1:5)

  R <- compute_LD(X, method = "population")
  expect_equal(dim(R), c(5L, 5L))
  expect_equal(unname(diag(R)), rep(1, 5))
  expect_equal(R[2, 4], 0)
})

test_that("compute_LD with larger matrix (100 SNPs) is fast and valid", {
  set.seed(99)
  X <- matrix(sample(0:2, 5000, replace = TRUE), nrow = 50, ncol = 100)
  colnames(X) <- paste0("rs", 1:100)

  R <- compute_LD(X, method = "sample")
  expect_equal(dim(R), c(100L, 100L))
  expect_equal(unname(diag(R)), rep(1, 100))
  expect_true(isSymmetric(R))
  expect_true(all(R >= -1 & R <= 1))
})

test_that("compute_LD population method with larger matrix is valid", {
  set.seed(99)
  X <- matrix(sample(0:2, 5000, replace = TRUE), nrow = 50, ncol = 100)
  colnames(X) <- paste0("rs", 1:100)

  R <- compute_LD(X, method = "population")
  expect_equal(dim(R), c(100L, 100L))
  expect_equal(unname(diag(R)), rep(1, 100))
  expect_true(isSymmetric(R))
})

test_that("compute_LD with perfectly correlated SNPs returns correlation of 1", {
  set.seed(42)
  col1 <- sample(0:2, 50, replace = TRUE)
  X <- matrix(c(col1, col1), ncol = 2)
  colnames(X) <- c("rs1", "rs2")

  R <- compute_LD(X, method = "sample")
  expect_equal(R[1, 2], 1.0, tolerance = 1e-10)
  expect_equal(R[2, 1], 1.0, tolerance = 1e-10)
})

test_that("compute_LD population method with trim_samples trims correctly", {
  set.seed(42)
  X <- matrix(sample(0:2, 33, replace = TRUE), nrow = 11, ncol = 3)
  colnames(X) <- paste0("rs", 1:3)

  R_trimmed <- compute_LD(X, method = "population", trim_samples = TRUE)
  expect_equal(dim(R_trimmed), c(3L, 3L))
  R_full <- compute_LD(X, method = "population", trim_samples = FALSE)
  expect_equal(dim(R_full), c(3L, 3L))
})

test_that("compute_LD with two monomorphic SNPs produces 0 off-diagonal", {
  X <- matrix(c(rep(1, 50), rep(2, 50)), nrow = 50, ncol = 2)
  colnames(X) <- c("mono1", "mono2")

  R <- compute_LD(X, method = "sample")
  expect_equal(R[1, 2], 0)
  expect_equal(R[2, 1], 0)
  expect_equal(unname(diag(R)), c(1, 1))
})

test_that("compute_LD preserves column names", {
  set.seed(42)
  X <- matrix(sample(0:2, 60, replace = TRUE), nrow = 20, ncol = 3)
  colnames(X) <- c("snp_alpha", "snp_beta", "snp_gamma")

  R <- compute_LD(X, method = "sample")
  expect_equal(colnames(R), c("snp_alpha", "snp_beta", "snp_gamma"))
  expect_equal(rownames(R), c("snp_alpha", "snp_beta", "snp_gamma"))
})

test_that("compute_LD with heavy missingness still produces valid matrix", {
  set.seed(42)
  X <- matrix(sample(0:2, 200, replace = TRUE), nrow = 40, ncol = 5)
  na_idx <- sample(length(X), size = floor(0.5 * length(X)))
  X[na_idx] <- NA
  colnames(X) <- paste0("rs", 1:5)

  R <- compute_LD(X, method = "sample")
  expect_true(all(!is.na(R)))
  expect_equal(unname(diag(R)), rep(1, 5))

  R_pop <- compute_LD(X, method = "population")
  expect_true(all(!is.na(R_pop)))
  expect_equal(unname(diag(R_pop)), rep(1, 5))
})

test_that("compute_LD with NA genotypes and sample method", {
  set.seed(42)
  X <- matrix(sample(0:2, 200, replace = TRUE), nrow = 50)
  X[1, 1] <- NA
  X[5, 3] <- NA
  colnames(X) <- paste0("rs", 1:4)

  R <- compute_LD(X, method = "sample")
  expect_true(all(!is.na(R)))
  expect_equal(unname(diag(R)), rep(1, 4))
})

test_that("compute_LD errors on NULL input", {
  expect_error(compute_LD(NULL), "X must be provided")
})

test_that("compute_LD sample vs population differ but are close", {
  set.seed(42)
  X <- matrix(sample(0:2, 500, replace = TRUE), nrow = 100)
  colnames(X) <- paste0("rs", 1:5)

  R_sample <- compute_LD(X, method = "sample")
  R_pop <- compute_LD(X, method = "population")

  expect_false(identical(R_sample, R_pop))
  expect_true(max(abs(R_sample - R_pop)) < 0.1)
})

# =============================================================================
# filter_X_with_Y
# =============================================================================

test_that("filter_X_with_Y preserves variants when Y has no missing data", {
  set.seed(42)
  X <- matrix(sample(0:2, 40, replace = TRUE), nrow = 10, ncol = 4)
  Y <- matrix(rnorm(20), nrow = 10, ncol = 2)
  colnames(Y) <- c("ctx1", "ctx2")
  rownames(X) <- rownames(Y) <- paste0("s", 1:10)

  result <- pecotmr:::filter_X_with_Y(X, Y, missing_rate_thresh = 1, maf_thresh = 0)
  expect_true(ncol(result) >= 1)
})

test_that("filter_X_with_Y drops variants monomorphic due to Y missingness", {
  X <- matrix(c(0, 0, 1, 2, 0, 1, 1, 2), nrow = 4, ncol = 2)
  Y <- matrix(c(1, NA, NA, 1, 1, 1, 1, NA), nrow = 4, ncol = 2)
  colnames(Y) <- c("ctx1", "ctx2")
  rownames(X) <- rownames(Y) <- paste0("s", 1:4)

  result <- pecotmr:::filter_X_with_Y(X, Y, missing_rate_thresh = 1, maf_thresh = 0)
  expect_true(is.matrix(result))
})

# =============================================================================
# matxMax
# =============================================================================

test_that("matxMax finds location of maximum", {
  mtx <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)
  result <- pecotmr:::matxMax(mtx)
  expect_equal(result[1], 2)
  expect_equal(result[2], 3)
})

test_that("matxMax finds max in row vector", {
  mtx <- matrix(c(1, 5, 3), nrow = 1)
  result <- pecotmr:::matxMax(mtx)
  expect_equal(result[1], 1)
  expect_equal(result[2], 2)
})

test_that("matxMax finds max in column vector", {
  mtx <- matrix(c(1, 5, 3), ncol = 1)
  result <- pecotmr:::matxMax(mtx)
  expect_equal(result[1], 2)
  expect_equal(result[2], 1)
})

test_that("matxMax with negative values finds the least negative", {
  mtx <- matrix(c(-10, -5, -20, -1), nrow = 2)
  result <- pecotmr:::matxMax(mtx)
  expect_equal(result[1], 2)
  expect_equal(result[2], 2)
})

# =============================================================================
# wald_test_pval
# =============================================================================

test_that("wald_test_pval computes correct p-values", {
  pval <- wald_test_pval(beta = 5, se = 1, n = 100)
  expect_true(pval < 0.001)

  pval_zero <- wald_test_pval(beta = 0, se = 1, n = 100)
  expect_equal(pval_zero, 1.0, tolerance = 1e-10)
})

test_that("wald_test_pval handles vector inputs", {
  betas <- c(0, 1, 2, 5)
  ses <- c(1, 1, 1, 1)
  pvals <- wald_test_pval(betas, ses, n = 100)
  expect_length(pvals, 4)
  expect_true(pvals[1] > pvals[4])
})

test_that("wald_test_pval is symmetric in beta sign", {
  pval_pos <- wald_test_pval(beta = 3, se = 1, n = 50)
  pval_neg <- wald_test_pval(beta = -3, se = 1, n = 50)
  expect_equal(pval_pos, pval_neg, tolerance = 1e-10)
})

test_that("wald_test_pval with very large beta gives p near 0", {
  pval <- wald_test_pval(beta = 100, se = 1, n = 1000)
  expect_true(pval < 1e-10)
})

test_that("wald_test_pval with very large se gives p near 1", {
  pval <- wald_test_pval(beta = 1, se = 1000, n = 100)
  expect_true(pval > 0.99)
})

# =============================================================================
# parse_region
# =============================================================================

test_that("parse_region parses valid region string", {
  result <- parse_region("chr1:100-200")
  expect_s3_class(result, "data.frame")
  expect_equal(result$chrom, "1")
  expect_equal(result$start, 100L)
  expect_equal(result$end, 200L)
})

test_that("parse_region handles X chromosome", {
  result <- parse_region("chrX:500-1000")
  expect_equal(result$chrom, "X")
  expect_equal(result$start, 500L)
  expect_equal(result$end, 1000L)
})

test_that("parse_region errors on invalid format", {
  expect_error(parse_region("1:100-200"), "format must be")
  expect_error(parse_region("chr1-100-200"), "format must be")
  expect_error(parse_region("chr1:abc-200"), "format must be")
})

test_that("parse_region returns non-string input unchanged", {
  df <- data.frame(chrom = 1, start = 100, end = 200)
  result <- parse_region(df)
  expect_identical(result, df)
})

test_that("parse_region returns non-single-string input unchanged", {
  input <- c("chr1:100-200", "chr2:300-400")
  result <- parse_region(input)
  expect_identical(result, input)
})

# =============================================================================
# parse_variant_id
# =============================================================================

test_that("parse_variant_id parses single variant with chr prefix", {
  result <- parse_variant_id("chr1:12345:A:G")
  expect_equal(result$chrom, 1L)
  expect_equal(result$pos, 12345L)
  expect_equal(result$A2, "A")
  expect_equal(result$A1, "G")
  conv <- attr(result, "convention")
  expect_true(conv$has_chr)
  expect_equal(conv$allele_sep, ":")
})

test_that("parse_variant_id parses single variant without chr prefix", {
  result <- parse_variant_id("5:12345:A:G")
  expect_equal(result$chrom, 5L)
  expect_equal(result$pos, 12345L)
  expect_equal(result$A2, "A")
  expect_equal(result$A1, "G")
  conv <- attr(result, "convention")
  expect_false(conv$has_chr)
})

test_that("parse_variant_id parses multiple variants", {
  ids <- c("chr1:100:A:G", "chr2:200:C:T", "chr3:300:G:A")
  result <- parse_variant_id(ids)
  expect_equal(nrow(result), 3)
  expect_equal(result$chrom, c(1L, 2L, 3L))
  expect_equal(result$pos, c(100L, 200L, 300L))
  expect_equal(result$A2, c("A", "C", "G"))
  expect_equal(result$A1, c("G", "T", "A"))
})

# =============================================================================
# detect_variant_convention
# =============================================================================

test_that("detect_variant_convention detects chr prefix and allele separators", {
  conv <- pecotmr:::detect_variant_convention(c("chr1:100:A:G", "chr2:200:C:T"))
  expect_true(conv$has_chr)
  expect_equal(conv$allele_sep, ":")
  expect_false(conv$has_build)

  conv2 <- pecotmr:::detect_variant_convention(c("1_100_A_G", "2_200_C_T"))
  expect_false(conv2$has_chr)
  expect_equal(conv2$allele_sep, "_")

  conv3 <- pecotmr:::detect_variant_convention(c("chr1:100:A:G:b38"))
  expect_true(conv3$has_build)

  conv4 <- pecotmr:::detect_variant_convention(c("chr1:100_A_G"))
  expect_true(conv4$has_chr)
  expect_equal(conv4$allele_sep, "_")

  conv5 <- pecotmr:::detect_variant_convention(c("1:100_A_G"))
  expect_false(conv5$has_chr)
  expect_equal(conv5$allele_sep, "_")
})

# =============================================================================
# normalize_variant_id
# =============================================================================

test_that("normalize_variant_id normalizes various formats", {
  expect_equal(normalize_variant_id("1_100_A_G"), "chr1:100:A:G")
  expect_equal(normalize_variant_id("chr1:100:A:G"), "chr1:100:A:G")
  expect_equal(normalize_variant_id("1:100:A:G"), "chr1:100:A:G")
  expect_equal(normalize_variant_id("chr1:100:A:G", chr_prefix = FALSE), "1:100:A:G")
  expect_equal(normalize_variant_id("chr1:100:A:G:b38"), "chr1:100:A:G")
  expect_equal(normalize_variant_id("chr1:100_A_G"), "chr1:100:A:G")
  conv <- pecotmr:::detect_variant_convention(c("chr1:100_A_G"))
  expect_equal(normalize_variant_id("1:200:C:T", convention = conv), "chr1:200_C_T")
})

# =============================================================================
# variant_id_to_df
# =============================================================================

test_that("variant_id_to_df handles colon-separated format", {
  ids <- c("1:100:A:G", "2:200:C:T")
  result <- pecotmr:::variant_id_to_df(ids)
  expect_equal(nrow(result), 2)
  expect_equal(result$chrom, c(1L, 2L))
  expect_equal(result$pos, c(100L, 200L))
  expect_equal(result$A2, c("A", "C"))
  expect_equal(result$A1, c("G", "T"))
})

test_that("variant_id_to_df handles underscore-separated format", {
  ids <- c("1:100_A_G", "2:200_C_T")
  result <- pecotmr:::variant_id_to_df(ids)
  expect_equal(nrow(result), 2)
  expect_equal(result$A2, c("A", "C"))
})

test_that("variant_id_to_df strips chr prefix", {
  ids <- c("chr1:100:A:G", "chr2:200:C:T")
  result <- pecotmr:::variant_id_to_df(ids)
  expect_equal(result$chrom, c(1L, 2L))
})

test_that("variant_id_to_df handles data.frame input with named columns", {
  df <- data.frame(chrom = c("chr1", "2"), pos = c(100, 200),
                   A2 = c("A", "C"), A1 = c("G", "T"))
  suppressWarnings(result <- pecotmr:::variant_id_to_df(df))
  expect_equal(result$chrom, c(1L, 2L))
  expect_equal(result$pos, c(100L, 200L))
})

test_that("variant_id_to_df handles 5-part IDs with build suffix", {
  ids <- c("chr1:100:A:G:b38", "chr2:200:T:C")
  result <- pecotmr:::variant_id_to_df(ids)
  expect_equal(ncol(result), 4)
  expect_equal(colnames(result), c("chrom", "pos", "A2", "A1"))
  expect_equal(result$chrom, c(1L, 2L))
  expect_equal(result$A2, c("A", "T"))
  expect_equal(result$A1, c("G", "C"))
})

test_that("variant_id_to_df handles mixed 4/5-part IDs", {
  ids <- c("1:100:A:G", "chr2:200:T:C:b38", "3:300:G:A:b37")
  suppressWarnings(result <- pecotmr:::variant_id_to_df(ids))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 4)
  expect_equal(result$A1, c("G", "C", "A"))
})

# =============================================================================
# get_nested_element
# =============================================================================

test_that("get_nested_element retrieves deeply nested values", {
  nested <- list(a = list(b = list(c = 42)))
  expect_equal(get_nested_element(nested, c("a", "b", "c")), 42)
})

test_that("get_nested_element returns NULL for NULL name_vector", {
  nested <- list(a = 1)
  expect_null(get_nested_element(nested, NULL))
})

test_that("get_nested_element errors on missing element", {
  nested <- list(a = list(b = 1))
  expect_error(get_nested_element(nested, c("a", "x")), "Element not found")
})

test_that("get_nested_element handles single level", {
  nested <- list(a = "hello")
  expect_equal(get_nested_element(nested, "a"), "hello")
})

test_that("get_nested_element skips empty strings", {
  nested <- list(a = list(b = 99))
  expect_equal(get_nested_element(nested, c("", "a", "b")), 99)
})

# =============================================================================
# region_to_df
# =============================================================================

test_that("region_to_df converts underscore-separated region IDs", {
  ids <- c("1_100_200", "2_300_400")
  result <- region_to_df(ids)
  expect_equal(nrow(result), 2)
  expect_equal(result$chrom, c(1L, 2L))
  expect_equal(result$start, c(100L, 300L))
  expect_equal(result$end, c(200L, 400L))
})

test_that("region_to_df handles chr prefix", {
  ids <- c("chr1_100_200")
  result <- region_to_df(ids)
  expect_equal(result$chrom, 1L)
})

test_that("region_to_df allows custom column names", {
  ids <- c("1_100_200")
  result <- region_to_df(ids, colnames = c("chr", "begin", "finish"))
  expect_true(all(c("chr", "begin", "finish") %in% colnames(result)))
})

# =============================================================================
# z_to_pvalue
# =============================================================================

test_that("z_to_pvalue returns correct p-values", {
  expect_equal(z_to_pvalue(0), 1.0, tolerance = 1e-10)
  expect_true(z_to_pvalue(1.96) < 0.05)
  expect_true(z_to_pvalue(-1.96) < 0.05)
  expect_equal(z_to_pvalue(1.96), z_to_pvalue(-1.96), tolerance = 1e-10)
})

test_that("z_to_pvalue handles vector input", {
  z <- c(0, 1, 2, 3)
  p <- z_to_pvalue(z)
  expect_length(p, 4)
  expect_true(all(diff(p) < 0))
})

test_that("z_to_pvalue handles extreme values", {
  expect_true(z_to_pvalue(10) < 1e-20)
  expect_true(z_to_pvalue(40) >= 0)
})

# =============================================================================
# z_to_beta_se
# =============================================================================

test_that("z_to_beta_se produces correct conversions", {
  z <- c(2.0, -1.0)
  maf <- c(0.3, 0.1)
  n <- 10000

  result <- pecotmr:::z_to_beta_se(z, maf, n)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("beta", "se", "maf") %in% names(result)))
  expect_equal(nrow(result), 2)
  expect_true(result$beta[1] > 0)
  expect_true(result$beta[2] < 0)
})

test_that("z_to_beta_se errors on mismatched lengths", {
  expect_error(pecotmr:::z_to_beta_se(c(1, 2), c(0.3), 1000), "same length")
})

test_that("z_to_beta_se adjusts MAF > 0.5", {
  result <- pecotmr:::z_to_beta_se(c(1.0), c(0.7), 1000)
  expect_equal(result$maf, 0.3)
})

# =============================================================================
# lbf_to_alpha
# =============================================================================

test_that("lbf_to_alpha converts matrix correctly", {
  lbf <- matrix(c(-0.5, 1.2, 0.3, 0.7, -1.1, 0.4), nrow = 2)
  colnames(lbf) <- paste0("v", 1:3)
  result <- lbf_to_alpha(lbf)

  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 3)
  expect_equal(rowSums(result), c(1, 1), tolerance = 1e-10)
})

test_that("lbf_to_alpha handles single column", {
  lbf <- matrix(c(0.5, 1.0), ncol = 1)
  colnames(lbf) <- "v1"
  result <- lbf_to_alpha(lbf)
  expect_equal(ncol(result), 1)
})

test_that("lbf_to_alpha handles all-zero row", {
  lbf <- matrix(c(0, 0, 0), nrow = 1)
  colnames(lbf) <- paste0("v", 1:3)
  result <- lbf_to_alpha(lbf)
  expect_true(all(result == 0))
})

# =============================================================================
# find_data
# =============================================================================

test_that("find_data retrieves from nested list at depth 2", {
  x <- list(a = list(val = 42), b = list(val = 99))
  result <- find_data(x, c(2, "val"))
  expect_equal(result, c(42, 99))
})

test_that("find_data returns list at depth 0", {
  x <- list(a = 1, b = 2)
  result <- find_data(x, c(0))
  expect_type(result, "list")
})

test_that("find_data with show_path=TRUE returns list structure", {
  x <- list(a = list(val = 10), b = list(val = 20))
  result <- find_data(x, c(2, "val"), show_path = TRUE)
  expect_type(result, "list")
  expect_true("a" %in% names(result))
  expect_true("b" %in% names(result))
  expect_equal(result$a, 10)
  expect_equal(result$b, 20)
})

test_that("find_data with rm_dup=TRUE removes duplicate values", {
  x <- list(
    a = list(val = 42),
    b = list(val = 42),
    c = list(val = 99)
  )
  result <- find_data(x, c(2, "val"), rm_dup = TRUE)
  expect_true("shared_list_names" %in% names(result))
})

test_that("find_data at depth 3 retrieves deeply nested values", {
  x <- list(
    level1_a = list(
      level2_a = list(target = "found_a"),
      level2_b = list(target = "found_b")
    ),
    level1_b = list(
      level2_c = list(target = "found_c")
    )
  )
  result <- find_data(x, c(3, "target"))
  expect_true("found_a" %in% result)
  expect_true("found_b" %in% result)
  expect_true("found_c" %in% result)
})

test_that("find_data with depth=1 and list_name returns element directly", {
  x <- list(a = 1, b = 2, c = 3)
  result <- find_data(x, c(1, "b"))
  expect_equal(result, 2)
})

test_that("find_data with rm_null=TRUE removes NULL results", {
  x2 <- list(
    a = list(val = 10),
    b = "not a list"
  )
  result <- find_data(x2, c(2, "val"), rm_null = TRUE)
  expect_equal(result, 10)
})

test_that("find_data with show_path=TRUE and rm_dup=TRUE at depth 2", {
  x <- list(
    a = list(val = 42),
    b = list(val = 42),
    c = list(val = 99)
  )
  result <- find_data(x, c(2, "val"), show_path = TRUE, rm_dup = TRUE)
  expect_type(result, "list")
  expect_true("shared_list_names" %in% names(result))
})

test_that("find_data with depth=1 and no list_name returns whole object", {
  x <- list(a = 1, b = 2)
  result <- find_data(x, c(1))
  expect_identical(result, x)
})

test_that("find_data with docall=list preserves list structure", {
  x <- list(
    a = list(val = c(1, 2)),
    b = list(val = c(3, 4))
  )
  result <- find_data(x, c(2, "val"), docall = list)
  expect_type(result, "list")
  expect_equal(result[[1]], c(1, 2))
  expect_equal(result[[2]], c(3, 4))
})
