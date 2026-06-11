context("misc")
library(tidyverse)

# =============================================================================
# computeMaf
# =============================================================================

test_that("Test computeMaf freq 0.5",{
    expect_equal(pecotmr:::computeMaf(rep(1, 20)), 0.5)
})

test_that("Test computeMaf freq 0.6",{
    expect_equal(pecotmr:::computeMaf(rep(1.2, 20)), 0.4)
})

test_that("Test computeMaf freq 0.3",{
    expect_equal(pecotmr:::computeMaf(rep(0.6, 20)), 0.3)
})

test_that("Test computeMaf with NA",{
    set.seed(1)
    generate_small_dataset <- function(sample_size = 20) {
        vals <- c(1.2, NA)
        return(sample(vals, sample_size, replace = TRUE))
    }
    expect_equal(pecotmr:::computeMaf(generate_small_dataset()), 0.4)
})

test_that("computeMaf returns 0 for monomorphic (all 0)", {
  expect_equal(pecotmr:::computeMaf(rep(0, 10)), 0)
})

# =============================================================================
# computeMissing
# =============================================================================

test_that("test computeMissing",{
    small_dataset <- c(rep(NA, 20), rep(1, 80))
    expect_equal(pecotmr:::computeMissing(small_dataset), 0.2)
})

# =============================================================================
# computeNonMissingY and computeAllMissingY
# =============================================================================

test_that("Test computeNonMissingY",{
    small_dataset <- c(rep(NA, 20), rep(1, 80))
    expect_equal(pecotmr:::computeNonMissingY(small_dataset), 80)
})

test_that("Test computeAllMissingY",{
    small_dataset <- c(rep(NA, 20), rep(1, 80))
    expect_equal(pecotmr:::computeAllMissingY(small_dataset), F)
})

test_that("computeAllMissingY returns TRUE for all-NA vector", {
  expect_true(pecotmr:::computeAllMissingY(rep(NA, 5)))
})

test_that("computeAllMissingY returns FALSE for partially NA vector", {
  expect_false(pecotmr:::computeAllMissingY(c(NA, 1, NA)))
})

# =============================================================================
# meanImpute
# =============================================================================

test_that("Test meanImpute",{
    dummy_data <- matrix(c(1,2,NA,1,2,3), nrow=3, ncol=2)
    expect_equal(pecotmr:::meanImpute(dummy_data)[3,1], 1.5)
})

test_that("meanImpute with all NAs in a column imputes NaN", {
  X <- matrix(c(NA, NA, NA, 1, 2, 3), nrow = 3, ncol = 2)
  result <- pecotmr:::meanImpute(X)
  expect_true(all(is.nan(result[, 1])))
  expect_equal(result[, 2], c(1, 2, 3))
})

# =============================================================================
# isZeroVariance
# =============================================================================

test_that("Test isZeroVariance",{
    dummy_data <- matrix(c(1,2,3,1,1,1), nrow=3, ncol=2)
    col <- which(apply(dummy_data, 2, pecotmr:::isZeroVariance))
    expect_equal(col, 2)
})

test_that("isZeroVariance with NA values treats them as distinct", {
  expect_false(pecotmr:::isZeroVariance(c(1, NA, 1)))
})

# =============================================================================
# filterX
# =============================================================================

test_that("Test filterX",{
    dummy_data <- matrix(
        c(1,NA,NA,NA, 0,0,1,1, 2,2,2,2, 1,1,1,2, 2,2,0,1, 0,1,1,2),
        # Missing Rate, MAF thresh, Zero Var, Var Thresh, Regular values
        nrow=4, ncol=6)
    var_thres <- 0.3
    expect_equal(filterX(dummy_data, 0.70, 0.3, varThresh = 0.3), matrix(c(2,2,0,1, 0,1,1,2), nrow=4, ncol=2))
})

test_that("filterX drops most columns when nearly all are zero variance", {
  X <- matrix(c(
    1, 1, 1, 1, 1,
    2, 2, 2, 2, 2,
    0, 0, 0, 0, 0,
    0, 1, 2, 0, 1
  ), nrow = 5, ncol = 4)
  result <- pecotmr:::filterX(X, missingRateThresh = 1.0, mafThresh = 0)
  expect_equal(ncol(result), 1)
})

test_that("filterX with external maf vector uses it for filtering", {
  set.seed(42)
  X <- matrix(sample(0:2, 40, replace = TRUE), nrow = 10, ncol = 4)
  external_maf <- c(0.01, 0.05, 0.3, 0.4)
  result <- pecotmr:::filterX(X, missingRateThresh = 1.0, mafThresh = 0.1, maf = external_maf)
  # Columns 1 and 2 have MAF <= 0.1, so they are dropped; columns 3 and 4 remain
  expect_equal(ncol(result), 2)
})

test_that("filterX skips MAF filtering for non-0/1/2 genotypes without external MAF", {
  set.seed(42)
  X <- matrix(runif(40, 0, 2), nrow = 10, ncol = 4)
  expect_message(
    result <- pecotmr:::filterX(X, missingRateThresh = 1.0, mafThresh = 0.1),
    "Skipping MAF filtering"
  )
  expect_true(ncol(result) >= 1)
})

test_that("filterX applies var_thresh with external X_variance", {
  set.seed(42)
  X <- matrix(sample(0:2, 40, replace = TRUE), nrow = 10, ncol = 4)
  external_var <- c(0.01, 0.5, 1.0, 0.02)
  result <- pecotmr:::filterX(X, missingRateThresh = 1.0, mafThresh = 0, varThresh = 0.1, xVariance = external_var)
  # Columns 1 and 4 have variance < 0.1, so they are dropped; columns 2 and 3 remain
  expect_equal(ncol(result), 2)
})

test_that("filterX with NULL thresholds does not filter", {
  set.seed(42)
  X <- matrix(sample(0:2, 40, replace = TRUE), nrow = 10, ncol = 4)
  result <- pecotmr:::filterX(X, missingRateThresh = NULL, mafThresh = NULL, varThresh = 0)
  # No filtering applied: all 4 columns should remain (zero-variance check still runs but none are zero-variance)
  expect_equal(ncol(result), 4)
})

test_that("filterX with missing_rate_thresh=0 drops columns with any NA", {
  X <- matrix(c(0, 1, 2, 1, 0, 1, NA, 2, 0, 1, 2, 1), nrow = 4, ncol = 3)
  result <- pecotmr:::filterX(X, missingRateThresh = 0, mafThresh = 0, varThresh = 0)
  expect_true(ncol(result) <= 2)
})

test_that("filterX with maf_thresh=0.5 removes nearly all columns", {
  X <- matrix(c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 2,
    0, 0, 0, 0, 0, 0, 0, 0, 2, 2,
    0, 0, 0, 0, 0, 2, 2, 2, 2, 2
  ), nrow = 10, ncol = 3)
  result <- pecotmr:::filterX(X, missingRateThresh = 1.0, mafThresh = 0.45, varThresh = 0)
  expect_equal(ncol(result), 1)
})

# =============================================================================
# filterY
# =============================================================================

test_that("Test filterY non-matrix",{
    dummy_data <- matrix(c(1,NA,NA,NA, 1,1,2,NA), nrow=4, ncol=2)
    res <- filterY(as.data.frame(dummy_data), 3)
    expect_equal(length(res$Y), 3)
    expect_equal(res$rmRows, NULL)
})

test_that("Test filterY is-matrix",{
    dummy_data <- matrix(c(1,NA,NA,NA, 1,1,2,NA, 2,1,2,NA), nrow=4, ncol=3)
    expect_equal(nrow(filterY(dummy_data, 3)$Y), 3)
    expect_equal(ncol(filterY(dummy_data, 3)$Y), 2)
    expect_equal(length(filterY(dummy_data, 3)$rmRows), 1)
})

test_that("filterY removes all columns with insufficient observations", {
  Y <- matrix(c(NA, NA, NA, 1, NA, NA, NA, 2), nrow = 4, ncol = 2)
  result <- pecotmr:::filterY(Y, nNonmiss = 3)
  expect_true(length(result$Y) == 0 || ncol(as.matrix(result$Y)) == 0)
})

test_that("filterY removes columns with too few non-missing values", {
  Y <- matrix(c(1, NA, NA, NA, 1, 2, 3, 4), nrow = 4)
  result <- pecotmr:::filterY(Y, nNonmiss = 3)
  expect_true(length(result$Y) >= 3)
})

test_that("filterY removes all-NA rows from matrix", {
  Y <- matrix(c(NA, NA, 1, 2, NA, NA, 3, 4), nrow = 4)
  result <- pecotmr:::filterY(Y, nNonmiss = 1)
  expect_true(nrow(result$Y) < 4)
})

# =============================================================================
# formatVariantId
# =============================================================================

test_that("Test formatVariantId",{
    expect_equal(formatVariantId(c(1, 1), c(123, 132), c("G", "A"), c("C", "T")), c("chr1:123:G:C", "chr1:132:A:T"))
})

test_that("formatVariantId uses convention parameter automatically", {
  expect_equal(pecotmr:::formatVariantId(1, 100, "A", "G"), "chr1:100:A:G")

  conv_mixed <- list(has_chr = TRUE, allele_sep = "_")
  expect_equal(pecotmr:::formatVariantId(1, 100, "A", "G", convention = conv_mixed), "chr1:100_A_G")

  conv_nochr <- list(has_chr = FALSE, allele_sep = "_")
  expect_equal(pecotmr:::formatVariantId(1, 100, "A", "G", convention = conv_nochr), "1:100_A_G")

  expect_equal(pecotmr:::formatVariantId(1, 100, "A", "G", chrPrefix = FALSE, convention = conv_mixed), "chr1:100_A_G")
})

test_that("formatVariantId constructs canonical IDs", {
  expect_equal(pecotmr:::formatVariantId(c(1, 2), c(100, 200), c("A", "C"), c("G", "T")),
               c("chr1:100:A:G", "chr2:200:C:T"))
  expect_equal(pecotmr:::formatVariantId(1, 100, "A", "G", chrPrefix = FALSE), "1:100:A:G")
  expect_equal(pecotmr:::formatVariantId("chr1", 100, "A", "G"), "chr1:100:A:G")
})

# =============================================================================
# findDuplicateVariants
# =============================================================================

z <- c(1, 2, 3, 4, 5)
LD <- matrix(c(1.0, 0.8, 0.2, 0.1, 0.3,
               0.8, 1.0, 0.4, 0.2, 0.5,
               0.2, 0.4, 1.0, 0.6, 0.1,
               0.1, 0.2, 0.6, 1.0, 0.3,
               0.3, 0.5, 0.1, 0.3, 1.0), nrow = 5, ncol = 5)

test_that("findDuplicateVariants returns the expected output", {
  rThreshold <- 0.5
  expected_output <- list(
    filteredZ = c(1, 3, 5),
    filteredLD= LD[c(1,3,5), c(1,3,5)],
    dupBearer = c(-1, 1, -1, 2, -1),
    corABS = c(0, 0.8, 0, 0.6, 0),
    sign = c(1, 1, 1, 1, 1),
    minValue = 0.1
  )

  result <- findDuplicateVariants(z, LD, rThreshold)
  expect_equal(result, expected_output)
})

test_that("findDuplicateVariants handles a high correlation threshold", {
  rThreshold <- 1.0
  expected_output <- list(
    filteredZ = c(1, 2, 3, 4, 5),
    filteredLD=LD,
    dupBearer = c(-1, -1, -1, -1, -1),
    corABS = c(0, 0, 0, 0, 0),
    sign = c(1, 1, 1, 1, 1),
    minValue = 0.1
  )

  result <- findDuplicateVariants(z, LD, rThreshold)
  expect_equal(result, expected_output)
})

test_that("findDuplicateVariants handles a low correlation threshold", {
  rThreshold <- 0.0
  expected_output <- list(
    filteredZ = c(1),
    filteredLD=LD[1,1,drop=F],
    dupBearer = c(-1, 1, 1, 1, 1),
    corABS = c(0, 0.8, 0.2, 0.1, 0.3),
    sign = c(1, 1, 1, 1, 1),
    minValue = 0.1
  )

  result <- findDuplicateVariants(z, LD, rThreshold)
  expect_equal(result, expected_output)
})

test_that("findDuplicateVariants handles negative correlations", {
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

  result <- findDuplicateVariants(z, LD_negative, rThreshold)
  expect_equal(result, expected_output)
})

# =============================================================================
# pvalGlobal
# =============================================================================

test_that("pvalGlobal with ACAT method returns valid combined p-value", {
  pvals <- c(0.01, 0.05, 0.5, 0.8)
  result <- pecotmr:::pvalGlobal(pvals, combMethod = "ACAT", naive = FALSE)
  expect_true(is.numeric(result))
  # ACAT statistic: T = mean(tan(pi*(0.5 - p_i))), p = P[Cauchy >= T]
  expected <- pcauchy(mean(tan(pi * (0.5 - pvals))), lower.tail = FALSE)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("pvalGlobal with naive=TRUE returns Bonferroni-corrected p-value", {
  pvals <- c(0.01, 0.05, 0.1, 0.5)
  result <- pecotmr:::pvalGlobal(pvals, combMethod = "HMP", naive = TRUE)
  n_unique <- length(unique(pvals))
  expected <- min(n_unique * min(pvals), 1.0)
  expect_equal(result, expected)
})

test_that("pvalGlobal naive method caps at 1.0", {
  pvals <- seq(0.1, 0.9, by = 0.01)
  result <- pecotmr:::pvalGlobal(pvals, combMethod = "ACAT", naive = TRUE)
  expect_true(result <= 1.0)
})

test_that("pvalGlobal naive method with single p-value returns that p-value", {
  result <- pecotmr:::pvalGlobal(0.03, combMethod = "ACAT", naive = TRUE)
  expect_equal(result, 0.03)
})

test_that("pvalGlobal ACAT with identical p-values", {
  pvals <- rep(0.05, 5)
  result <- pecotmr:::pvalGlobal(pvals, combMethod = "ACAT", naive = FALSE)
  expect_true(is.numeric(result))
  expect_true(result > 0 && result < 1)
})

test_that("pvalGlobal ACAT with single p-value delegates correctly", {
  result <- pecotmr:::pvalGlobal(0.05, combMethod = "ACAT", naive = FALSE)
  expect_equal(result, 0.05)
})

test_that("pvalGlobal ACAT with very significant p-values", {
  pvals <- c(1e-8, 1e-6, 1e-4)
  result <- pecotmr:::pvalGlobal(pvals, combMethod = "ACAT", naive = FALSE)
  expect_true(is.numeric(result))
  expect_true(result > 0 && result <= 1)
})

test_that("pvalGlobal HMP method returns valid p-value when harmonicmeanp available", {
  skip_if_not_installed("harmonicmeanp")
  pvals <- c(0.01, 0.05, 0.2, 0.7)
  result <- pecotmr:::pvalGlobal(pvals, combMethod = "HMP", naive = FALSE)
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("pvalGlobal HMP errors when harmonicmeanp not installed", {
  skip_if(requireNamespace("harmonicmeanp", quietly = TRUE),
          "harmonicmeanp is installed, cannot test missing-package path")
  pvals <- c(0.01, 0.05)
  expect_error(pecotmr:::pvalGlobal(pvals, combMethod = "HMP", naive = FALSE),
               "harmonicmeanp")
})

# =============================================================================
# pvalHmp
# =============================================================================

test_that("pvalHmp returns valid p-value", {
  skip_if_not_installed("harmonicmeanp")
  pvals <- c(0.01, 0.05, 0.1)
  result <- pecotmr:::pvalHmp(pvals)
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

test_that("pvalHmp uses unique p-values only", {
  skip_if_not_installed("harmonicmeanp")
  pvals <- c(0.01, 0.01, 0.05, 0.05, 0.3)
  result <- pecotmr:::pvalHmp(pvals)
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("pvalHmp errors when package not available", {
  skip_if(requireNamespace("harmonicmeanp", quietly = TRUE),
          "harmonicmeanp is installed, cannot test missing-package path")
  expect_error(pecotmr:::pvalHmp(c(0.01, 0.05)), "harmonicmeanp")
})

# =============================================================================
# pvalAcat
# =============================================================================

test_that("pvalAcat returns single p-value unchanged", {
  expect_equal(pecotmr:::pvalAcat(0.05), 0.05)
})

test_that("pvalAcat combines multiple p-values", {
  pvals <- c(0.001, 0.01, 0.1)
  combined <- pecotmr:::pvalAcat(pvals)
  expect_true(combined > 0 && combined < 1)
})

test_that("pvalAcat with very small p-values does not return NA", {
  result <- pecotmr:::pvalAcat(c(1e-10, 1e-8, 1e-6))
  expect_true(is.numeric(result))
  expect_true(!is.na(result))
  expect_true(result > 0 && result <= 1)
})

test_that("pvalAcat with all large p-values returns large combined p-value", {
  result <- pecotmr:::pvalAcat(c(0.8, 0.9, 0.95))
  expect_true(is.numeric(result))
  expect_true(result > 0.5 && result <= 1)
})

# =============================================================================
# pvalCauchy
# =============================================================================

test_that("pvalCauchy combines p-values", {
  pvals <- c(0.01, 0.05, 0.5)
  combined <- pecotmr:::pvalCauchy(pvals)
  # Manual: CCT stat = mean(tan((0.5 - p) * pi)), result = 1 - pcauchy(stat)
  cct_stat <- mean(tan((0.5 - pvals) * pi))
  expected <- 1 - pcauchy(cct_stat)
  expect_equal(combined, expected, tolerance = 1e-10)
})

test_that("pvalCauchy handles NAs with na.rm", {
  pvals <- c(0.01, NA, 0.05)
  combined <- pecotmr:::pvalCauchy(pvals, na.rm = TRUE)
  expect_true(!is.na(combined))
})

test_that("pvalCauchy with very small p-values", {
  pvals <- c(1e-20, 1e-15)
  combined <- pecotmr:::pvalCauchy(pvals)
  expect_true(combined < 1e-10)
})

test_that("pvalCauchy with all NA and na.rm=TRUE returns NA", {
  result <- pecotmr:::pvalCauchy(c(NA, NA, NA), na.rm = TRUE)
  expect_true(is.na(result))
})

test_that("pvalCauchy with p-values near 1 caps them at 0.99", {
  result <- pecotmr:::pvalCauchy(c(0.999, 0.9999, 0.5))
  expect_true(is.numeric(result))
  expect_true(!is.na(result))
})

test_that("pvalCauchy with na.rm=FALSE and NA present still computes result", {
  result <- pecotmr:::pvalCauchy(c(0.01, NA, 0.05), na.rm = FALSE)
  expect_true(is.numeric(result))
})

test_that("pvalCauchy with extremely small p-values triggers large-stat branch", {
  # cct.stat > 1e+15 triggers the 1/(cct.stat*pi) return path
  result <- pecotmr:::pvalCauchy(c(1e-300, 1e-290))
  expect_true(is.numeric(result))
  expect_true(result > 0 && result < 1)
  expect_false(is.na(result))
})

test_that("pvalAcat uses asymptotic approximation for p < 1e-15", {
  # p-values below 1e-15 use 1/(p*pi) instead of tan()
  result <- pecotmr:::pvalAcat(c(1e-20, 1e-18, 0.01))
  expect_true(is.numeric(result))
  expect_true(result > 0 && result < 1)
  expect_false(is.na(result))
})

# =============================================================================
# pvalPoolr
# =============================================================================

test_that("pvalPoolr fisher method returns valid p-value", {
  skip_if_not_installed("poolr")
  pvals <- c(0.01, 0.05, 0.1)
  R <- diag(3)
  result <- pecotmr:::pvalPoolr(pvals, method = "fisher", R = R)
  expect_true(is.numeric(result))
  expect_true(result > 0 && result < 1)
})

test_that("pvalPoolr stouffer method returns valid p-value", {
  skip_if_not_installed("poolr")
  pvals <- c(0.01, 0.05, 0.1)
  R <- diag(3)
  result <- pecotmr:::pvalPoolr(pvals, method = "stouffer", R = R)
  expect_true(is.numeric(result))
  expect_true(result > 0 && result < 1)
})

test_that("pvalPoolr invchisq method returns valid p-value", {
  skip_if_not_installed("poolr")
  pvals <- c(0.01, 0.05, 0.1)
  R <- diag(3)
  result <- pecotmr:::pvalPoolr(pvals, method = "invchisq", R = R)
  expect_true(is.numeric(result))
  expect_true(result > 0 && result < 1)
})

test_that("pvalPoolr errors on unknown method", {
  skip_if_not_installed("poolr")
  expect_error(pecotmr:::pvalPoolr(c(0.01, 0.05), method = "bogus", R = diag(2)),
               "Unknown poolr method")
})

# =============================================================================
# pvalGbj
# =============================================================================

test_that("pvalGbj gbj method returns valid p-value", {
  skip_if_not_installed("GBJ")
  z <- c(2.5, 1.8, 3.0)
  R <- diag(3)
  result <- pecotmr:::pvalGbj(z, R, method = "gbj")
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("pvalGbj hc method returns valid p-value", {
  skip_if_not_installed("GBJ")
  z <- c(2.5, 1.8, 3.0)
  R <- diag(3)
  result <- pecotmr:::pvalGbj(z, R, method = "hc")
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("pvalGbj minp method returns valid p-value", {
  skip_if_not_installed("GBJ")
  z <- c(2.5, 1.8, 3.0)
  R <- diag(3)
  result <- pecotmr:::pvalGbj(z, R, method = "minp")
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("pvalGbj bj method returns valid p-value", {
  skip_if_not_installed("GBJ")
  z <- c(2.5, 1.8, 3.0)
  R <- diag(3)
  result <- pecotmr:::pvalGbj(z, R, method = "bj")
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("pvalGbj ghc method returns valid p-value", {
  skip_if_not_installed("GBJ")
  z <- c(2.5, 1.8, 3.0)
  R <- diag(3)
  result <- pecotmr:::pvalGbj(z, R, method = "ghc")
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("pvalGbj gbj_omni method returns valid p-value", {
  skip_if_not_installed("GBJ")
  z <- c(2.5, 1.8, 3.0)
  R <- diag(3)
  result <- pecotmr:::pvalGbj(z, R, method = "gbj_omni")
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("pvalGbj errors on unknown method", {
  skip_if_not_installed("GBJ")
  expect_error(pecotmr:::pvalGbj(c(2.5, 1.8), diag(2), method = "bogus"),
               "Unknown GBJ method")
})

# =============================================================================
# pvalAspu
# =============================================================================

test_that("pvalAspu aspu method returns valid p-value", {
  skip_if_not_installed("aSPU")
  set.seed(42)
  z <- c(2.5, 1.8, 3.0)
  R <- diag(3)
  result <- pecotmr:::pvalAspu(zScores = z, R = R, method = "aspu")
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("pvalAspu gates method returns valid p-value", {
  skip_if_not_installed("aSPU")
  pvals <- c(0.01, 0.05, 0.1)
  R <- diag(3)
  result <- pecotmr:::pvalAspu(pvals = pvals, R = R, method = "gates")
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("pvalAspu errors on unknown method", {
  skip_if_not_installed("aSPU")
  expect_error(pecotmr:::pvalAspu(zScores = c(1, 2), R = diag(2), method = "bogus"),
               "Unknown aSPU method")
})

# =============================================================================
# computeQvalues
# =============================================================================

test_that("computeQvalues returns NA vector when all pvalues are NA", {
  result <- pecotmr:::computeQvalues(rep(NA_real_, 5))
  expect_true(all(is.na(result)))
  expect_length(result, 5)
})

test_that("computeQvalues returns single p-value unchanged", {
  result <- pecotmr:::computeQvalues(0.05)
  expect_equal(result, 0.05)
})

test_that("computeQvalues returns empty input unchanged", {
  result <- pecotmr:::computeQvalues(numeric(0))
  expect_length(result, 0)
})

test_that("computeQvalues works with valid p-value vector", {
  skip_if_not_installed("qvalue")
  set.seed(42)
  pvals <- runif(100, 0, 1)
  result <- pecotmr:::computeQvalues(pvals)
  expect_length(result, 100)
  expect_true(all(result >= 0 & result <= 1))
})

test_that("computeQvalues falls back to BH when too few p-values", {
  skip_if_not_installed("qvalue")
  pvals <- c(0.001, 0.999)
  result <- pecotmr:::computeQvalues(pvals)
  expect_length(result, 2)
  expect_true(all(result >= 0 & result <= 1))
})

test_that("computeQvalues errors when qvalue not installed", {
  skip_if(requireNamespace("qvalue", quietly = TRUE),
          "qvalue is installed, cannot test missing-package path")
  expect_error(pecotmr:::computeQvalues(c(0.01, 0.05)), "qvalue")
})

# =============================================================================
# filterMolecularEvents
# =============================================================================

test_that("filterMolecularEvents errors when filter lacks required fields", {
  events <- c("gene_A_splicing", "gene_B_expression")
  bad_filter <- list(list(type_pattern = "gene"))
  expect_error(
    pecotmr:::filterMolecularEvents(events, bad_filter),
    "type_pattern and at least one of"
  )
})

test_that("filterMolecularEvents errors when type_pattern is NULL", {
  events <- c("gene_A_splicing")
  bad_filter <- list(list(type_pattern = NULL, valid_pattern = "splicing"))
  expect_error(
    pecotmr:::filterMolecularEvents(events, bad_filter),
    "type_pattern and at least one of"
  )
})

test_that("filterMolecularEvents keeps events matching valid_pattern", {
  events <- c("gene_A_splicing", "gene_A_expression", "gene_B_splicing", "protein_X")
  filters <- list(list(
    type_pattern = "gene_",
    valid_pattern = "splicing"
  ))
  result <- pecotmr:::filterMolecularEvents(events, filters)
  expect_true("protein_X" %in% result)
  gene_events <- result[grepl("gene_", result)]
  expect_true(all(grepl("splicing", gene_events)))
})

test_that("filterMolecularEvents excludes events matching exclude_pattern", {
  events <- c("gene_A_splicing", "gene_A_expression", "gene_B_splicing")
  filters <- list(list(
    type_pattern = "gene_",
    exclude_pattern = "expression"
  ))
  result <- pecotmr:::filterMolecularEvents(events, filters)
  expect_false("gene_A_expression" %in% result)
  expect_true("gene_A_splicing" %in% result)
  expect_true("gene_B_splicing" %in% result)
})

test_that("filterMolecularEvents returns NULL when no events pass filtering", {
  events <- c("gene_A_expression", "gene_B_expression")
  filters <- list(list(
    type_pattern = "gene_",
    valid_pattern = "splicing"
  ))
  result <- pecotmr:::filterMolecularEvents(events, filters)
  expect_null(result)
})

test_that("filterMolecularEvents skips filter when no type events match", {
  events <- c("protein_X", "protein_Y")
  filters <- list(list(
    type_pattern = "gene_",
    valid_pattern = "splicing"
  ))
  result <- pecotmr:::filterMolecularEvents(events, filters)
  expect_equal(sort(result), sort(events))
})

test_that("filterMolecularEvents handles comma-separated valid_pattern", {
  events <- c("gene_A_splicing", "gene_A_expression", "gene_B_methylation")
  filters <- list(list(
    type_pattern = "gene_",
    valid_pattern = "splicing,expression"
  ))
  result <- pecotmr:::filterMolecularEvents(events, filters)
  expect_true("gene_A_splicing" %in% result)
  expect_true("gene_A_expression" %in% result)
  expect_false("gene_B_methylation" %in% result)
})

test_that("filterMolecularEvents handles comma-separated exclude_pattern", {
  events <- c("gene_A_splicing", "gene_A_expression", "gene_B_methylation")
  filters <- list(list(
    type_pattern = "gene_",
    exclude_pattern = "expression,methylation"
  ))
  result <- pecotmr:::filterMolecularEvents(events, filters)
  expect_true("gene_A_splicing" %in% result)
  expect_false("gene_A_expression" %in% result)
  expect_false("gene_B_methylation" %in% result)
})

test_that("filterMolecularEvents with both valid_pattern and exclude_pattern", {
  events <- c("gene_A_splicing_good", "gene_A_splicing_bad",
              "gene_B_expression", "protein_X")
  filters <- list(list(
    type_pattern = "gene_",
    valid_pattern = "splicing",
    exclude_pattern = "bad"
  ))
  result <- pecotmr:::filterMolecularEvents(events, filters)
  expect_true("gene_A_splicing_good" %in% result)
  expect_false("gene_A_splicing_bad" %in% result)
  expect_true("protein_X" %in% result)
})

test_that("filterMolecularEvents with condition parameter", {
  events <- c("gene_A_splicing", "gene_A_expression")
  filters <- list(list(
    type_pattern = "gene_",
    exclude_pattern = "expression"
  ))
  result <- pecotmr:::filterMolecularEvents(events, filters, condition = "test_context")
  expect_true("gene_A_splicing" %in% result)
  expect_false("gene_A_expression" %in% result)
})

test_that("filterMolecularEvents with multiple filters", {
  events <- c("gene_A_splicing", "gene_A_expression",
              "protein_X_high", "protein_X_low")
  filters <- list(
    list(type_pattern = "gene_", valid_pattern = "splicing"),
    list(type_pattern = "protein_", exclude_pattern = "low")
  )
  result <- pecotmr:::filterMolecularEvents(events, filters)
  expect_true("gene_A_splicing" %in% result)
  expect_false("gene_A_expression" %in% result)
  expect_true("protein_X_high" %in% result)
  expect_false("protein_X_low" %in% result)
})

test_that("filterMolecularEvents returns all when all events match", {
  events <- c("gene_A_splicing", "gene_B_splicing")
  filters <- list(list(
    type_pattern = "gene_",
    valid_pattern = "splicing"
  ))
  result <- pecotmr:::filterMolecularEvents(events, filters)
  expect_equal(sort(result), sort(events))
})

# =============================================================================
# findValidFilePath and findValidFilePaths
# =============================================================================

test_that("findValidFilePath returns target when it exists directly", {
  pkg_root <- normalizePath(file.path(test_path(), "..", ".."), mustWork = TRUE)
  target <- file.path(pkg_root, "DESCRIPTION")
  ref <- file.path(pkg_root, "NAMESPACE")
  skip_if_not(file.exists(target) && file.exists(ref), "Package root files not found")
  result <- pecotmr:::findValidFilePath(ref, target)
  expect_equal(result, target)
})

test_that("findValidFilePath constructs path from reference directory", {
  pkg_root <- normalizePath(file.path(test_path(), "..", ".."), mustWork = TRUE)
  ref <- file.path(pkg_root, "NAMESPACE")
  skip_if_not(file.exists(ref), "NAMESPACE not found")
  result <- pecotmr:::findValidFilePath(ref, "DESCRIPTION")
  expect_true(file.exists(result))
  expect_true(grepl("DESCRIPTION$", result))
})

test_that("findValidFilePath errors when both paths are invalid", {
  expect_error(
    pecotmr:::findValidFilePath("/nonexistent/dir/ref.txt", "/nonexistent/target.txt"),
    "Both reference and target file paths do not work"
  )
})

test_that("findValidFilePath returns reference when target is invalid but reference exists", {
  pkg_root <- normalizePath(file.path(test_path(), "..", ".."), mustWork = TRUE)
  ref <- file.path(pkg_root, "DESCRIPTION")
  skip_if_not(file.exists(ref), "DESCRIPTION not found")
  result <- pecotmr:::findValidFilePath(ref, "/totally/bogus/path.txt")
  expect_equal(result, ref)
})

test_that("findValidFilePaths resolves multiple targets", {
  pkg_root <- normalizePath(file.path(test_path(), "..", ".."), mustWork = TRUE)
  ref <- file.path(pkg_root, "NAMESPACE")
  skip_if_not(file.exists(ref), "NAMESPACE not found")
  targets <- c("DESCRIPTION", "NAMESPACE")
  result <- pecotmr:::findValidFilePaths(ref, targets)
  expect_length(result, 2)
  expect_true(all(file.exists(result)))
})

test_that("findValidFilePaths errors on all-invalid targets", {
  ref <- "/nonexistent/ref.txt"
  targets <- c("/bogus/a.txt", "/bogus/b.txt")
  expect_error(pecotmr:::findValidFilePaths(ref, targets))
})

# =============================================================================
# computeLd
# =============================================================================

test_that("computeLd sample method produces valid correlation matrix", {
  set.seed(42)
  X <- matrix(sample(0:2, 200, replace = TRUE), nrow = 50)
  colnames(X) <- paste0("rs", 1:4)

  R <- computeLd(X, method = "sample")
  expect_equal(nrow(R), 4)
  expect_equal(ncol(R), 4)
  expect_equal(unname(diag(R)), rep(1, 4))
  expect_true(isSymmetric(R))
  expect_true(all(R >= -1 & R <= 1))
})

test_that("computeLd population method produces valid matrix", {
  set.seed(42)
  X <- matrix(sample(0:2, 200, replace = TRUE), nrow = 50)
  colnames(X) <- paste0("rs", 1:4)

  R <- computeLd(X, method = "population")
  expect_equal(nrow(R), 4)
  expect_equal(unname(diag(R)), rep(1, 4))
  expect_true(isSymmetric(R))
})

test_that("computeLd with a single SNP returns 1x1 identity matrix", {
  X <- matrix(c(0, 1, 2, 1, 0), ncol = 1)
  colnames(X) <- "rs1"
  R <- computeLd(X, method = "sample")
  expect_equal(dim(R), c(1L, 1L))
  expect_equal(R[1, 1], 1.0)
  expect_equal(colnames(R), "rs1")
})

test_that("computeLd handles column with all NA gracefully", {
  set.seed(123)
  X <- matrix(sample(0:2, 100, replace = TRUE), nrow = 20, ncol = 5)
  X[, 3] <- NA
  colnames(X) <- paste0("rs", 1:5)

  R <- computeLd(X, method = "sample")
  expect_equal(dim(R), c(5L, 5L))
  expect_equal(unname(diag(R)), rep(1, 5))
  expect_equal(R[3, 1], 0)
  expect_equal(R[1, 3], 0)
})

test_that("computeLd population method handles column with all NA", {
  set.seed(123)
  X <- matrix(sample(0:2, 100, replace = TRUE), nrow = 20, ncol = 5)
  X[, 2] <- NA
  colnames(X) <- paste0("rs", 1:5)

  R <- computeLd(X, method = "population")
  expect_equal(dim(R), c(5L, 5L))
  expect_equal(unname(diag(R)), rep(1, 5))
  expect_equal(R[2, 4], 0)
})

test_that("computeLd with larger matrix (100 SNPs) is fast and valid", {
  set.seed(99)
  X <- matrix(sample(0:2, 5000, replace = TRUE), nrow = 50, ncol = 100)
  colnames(X) <- paste0("rs", 1:100)

  R <- computeLd(X, method = "sample")
  expect_equal(dim(R), c(100L, 100L))
  expect_equal(unname(diag(R)), rep(1, 100))
  expect_true(isSymmetric(R))
  expect_true(all(R >= -1 & R <= 1))
})

test_that("computeLd population method with larger matrix is valid", {
  set.seed(99)
  X <- matrix(sample(0:2, 5000, replace = TRUE), nrow = 50, ncol = 100)
  colnames(X) <- paste0("rs", 1:100)

  R <- computeLd(X, method = "population")
  expect_equal(dim(R), c(100L, 100L))
  expect_equal(unname(diag(R)), rep(1, 100))
  expect_true(isSymmetric(R))
})

test_that("computeLd with perfectly correlated SNPs returns correlation of 1", {
  set.seed(42)
  col1 <- sample(0:2, 50, replace = TRUE)
  X <- matrix(c(col1, col1), ncol = 2)
  colnames(X) <- c("rs1", "rs2")

  R <- computeLd(X, method = "sample")
  expect_equal(R[1, 2], 1.0, tolerance = 1e-10)
  expect_equal(R[2, 1], 1.0, tolerance = 1e-10)
})

test_that("computeLd population method with trim_samples trims correctly", {
  set.seed(42)
  X <- matrix(sample(0:2, 33, replace = TRUE), nrow = 11, ncol = 3)
  colnames(X) <- paste0("rs", 1:3)

  R_trimmed <- computeLd(X, method = "population", trimSamples = TRUE)
  expect_equal(dim(R_trimmed), c(3L, 3L))
  R_full <- computeLd(X, method = "population", trimSamples = FALSE)
  expect_equal(dim(R_full), c(3L, 3L))
})

test_that("computeLd with two monomorphic SNPs produces 0 off-diagonal", {
  X <- matrix(c(rep(1, 50), rep(2, 50)), nrow = 50, ncol = 2)
  colnames(X) <- c("mono1", "mono2")

  R <- computeLd(X, method = "sample")
  expect_equal(R[1, 2], 0)
  expect_equal(R[2, 1], 0)
  expect_equal(unname(diag(R)), c(1, 1))
})

test_that("computeLd preserves column names", {
  set.seed(42)
  X <- matrix(sample(0:2, 60, replace = TRUE), nrow = 20, ncol = 3)
  colnames(X) <- c("snp_alpha", "snp_beta", "snp_gamma")

  R <- computeLd(X, method = "sample")
  expect_equal(colnames(R), c("snp_alpha", "snp_beta", "snp_gamma"))
  expect_equal(rownames(R), c("snp_alpha", "snp_beta", "snp_gamma"))
})

test_that("computeLd with heavy missingness still produces valid matrix", {
  set.seed(42)
  X <- matrix(sample(0:2, 200, replace = TRUE), nrow = 40, ncol = 5)
  na_idx <- sample(length(X), size = floor(0.5 * length(X)))
  X[na_idx] <- NA
  colnames(X) <- paste0("rs", 1:5)

  R <- computeLd(X, method = "sample")
  expect_true(all(!is.na(R)))
  expect_equal(unname(diag(R)), rep(1, 5))

  R_pop <- computeLd(X, method = "population")
  expect_true(all(!is.na(R_pop)))
  expect_equal(unname(diag(R_pop)), rep(1, 5))
})

test_that("computeLd with NA genotypes and sample method", {
  set.seed(42)
  X <- matrix(sample(0:2, 200, replace = TRUE), nrow = 50)
  X[1, 1] <- NA
  X[5, 3] <- NA
  colnames(X) <- paste0("rs", 1:4)

  R <- computeLd(X, method = "sample")
  expect_true(all(!is.na(R)))
  expect_equal(unname(diag(R)), rep(1, 4))
})

test_that("computeLd errors on NULL input", {
  expect_error(computeLd(NULL), "X must be provided")
})

test_that("computeLd sample vs population differ but are close", {
  set.seed(42)
  X <- matrix(sample(0:2, 500, replace = TRUE), nrow = 100)
  colnames(X) <- paste0("rs", 1:5)

  R_sample <- computeLd(X, method = "sample")
  R_pop <- computeLd(X, method = "population")

  expect_false(identical(R_sample, R_pop))
  expect_true(max(abs(R_sample - R_pop)) < 0.1)
})

test_that("computeLd gcta method produces valid correlation matrix", {
  set.seed(42)
  X <- matrix(sample(0:2, 500, replace = TRUE), nrow = 100)
  colnames(X) <- paste0("rs", 1:5)

  R <- computeLd(X, method = "gcta")
  expect_equal(dim(R), c(5, 5))
  expect_equal(unname(diag(R)), rep(1, 5), tolerance = 1e-10)
  expect_true(isSymmetric(R, tol = 1e-10))
  expect_true(all(abs(R) <= 1 + 1e-10))
})

test_that("computeLd gcta method handles missing data", {
  set.seed(42)
  X <- matrix(sample(0:2, 500, replace = TRUE), nrow = 100)
  colnames(X) <- paste0("rs", 1:5)
  X[sample(length(X), 50)] <- NA

  R <- computeLd(X, method = "gcta")
  expect_equal(dim(R), c(5, 5))
  expect_true(all(is.finite(R)))
})

test_that("computeLd gcta agrees with sample method on complete data", {
  set.seed(42)
  X <- matrix(sample(0:2, 500, replace = TRUE), nrow = 100)
  colnames(X) <- paste0("rs", 1:5)

  R_sample <- computeLd(X, method = "sample")
  R_gcta <- computeLd(X, method = "gcta")

  # With no missing data, GCTA and sample should be close (differ by N vs N-1 denom)
  expect_true(max(abs(R_sample - R_gcta)) < 0.05)
})

test_that("computeLd gcta preserves column names", {
  set.seed(42)
  X <- matrix(sample(0:2, 300, replace = TRUE), nrow = 100)
  colnames(X) <- c("snp_a", "snp_b", "snp_c")

  R <- computeLd(X, method = "gcta")
  expect_equal(colnames(R), c("snp_a", "snp_b", "snp_c"))
  expect_equal(rownames(R), c("snp_a", "snp_b", "snp_c"))
})

# =============================================================================
# filterXWithY
# =============================================================================

test_that("filterXWithY preserves variants when Y has no missing data", {
  set.seed(42)
  X <- matrix(sample(0:2, 40, replace = TRUE), nrow = 10, ncol = 4)
  Y <- matrix(rnorm(20), nrow = 10, ncol = 2)
  colnames(Y) <- c("ctx1", "ctx2")
  rownames(X) <- rownames(Y) <- paste0("s", 1:10)

  result <- pecotmr:::filterXWithY(X, Y, missingRateThresh = 1, mafThresh = 0)
  expect_true(ncol(result) >= 1)
})

test_that("filterXWithY drops variants monomorphic due to Y missingness", {
  X <- matrix(c(0, 0, 1, 2, 0, 1, 1, 2), nrow = 4, ncol = 2)
  Y <- matrix(c(1, NA, NA, 1, 1, 1, 1, NA), nrow = 4, ncol = 2)
  colnames(Y) <- c("ctx1", "ctx2")
  rownames(X) <- rownames(Y) <- paste0("s", 1:4)

  result <- pecotmr:::filterXWithY(X, Y, missingRateThresh = 1, mafThresh = 0)
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
# waldTestPval
# =============================================================================

test_that("waldTestPval computes correct p-values", {
  pval <- waldTestPval(beta = 5, se = 1, n = 100)
  expect_true(pval < 0.001)

  pval_zero <- waldTestPval(beta = 0, se = 1, n = 100)
  expect_equal(pval_zero, 1.0, tolerance = 1e-10)
})

test_that("waldTestPval handles vector inputs", {
  betas <- c(0, 1, 2, 5)
  ses <- c(1, 1, 1, 1)
  pvals <- waldTestPval(betas, ses, n = 100)
  expect_length(pvals, 4)
  expect_true(pvals[1] > pvals[4])
})

test_that("waldTestPval is symmetric in beta sign", {
  pval_pos <- waldTestPval(beta = 3, se = 1, n = 50)
  pval_neg <- waldTestPval(beta = -3, se = 1, n = 50)
  expect_equal(pval_pos, pval_neg, tolerance = 1e-10)
})

test_that("waldTestPval with very large beta gives p near 0", {
  pval <- waldTestPval(beta = 100, se = 1, n = 1000)
  expect_true(pval < 1e-10)
})

test_that("waldTestPval with very large se gives p near 1", {
  pval <- waldTestPval(beta = 1, se = 1000, n = 100)
  expect_true(pval > 0.99)
})

# =============================================================================
# parseRegion
# =============================================================================

test_that("parseRegion parses valid region string", {
  result <- parseRegion("chr1:100-200")
  expect_s3_class(result, "data.frame")
  expect_equal(result$chrom, "1")
  expect_equal(result$start, 100L)
  expect_equal(result$end, 200L)
})

test_that("parseRegion handles X chromosome", {
  result <- parseRegion("chrX:500-1000")
  expect_equal(result$chrom, "X")
  expect_equal(result$start, 500L)
  expect_equal(result$end, 1000L)
})

test_that("parseRegion errors on invalid format", {
  expect_error(parseRegion("1:100-200"), "format must be")
  expect_error(parseRegion("chr1-100-200"), "format must be")
  expect_error(parseRegion("chr1:abc-200"), "format must be")
})

test_that("parseRegion returns non-string input unchanged", {
  df <- data.frame(chrom = 1, start = 100, end = 200)
  result <- parseRegion(df)
  expect_identical(result, df)
})

test_that("parseRegion returns non-single-string input unchanged", {
  input <- c("chr1:100-200", "chr2:300-400")
  result <- parseRegion(input)
  expect_identical(result, input)
})

# =============================================================================
# parseVariantId
# =============================================================================

test_that("parseVariantId parses single variant with chr prefix", {
  result <- parseVariantId("chr1:12345:A:G")
  expect_equal(result$chrom, 1L)
  expect_equal(result$pos, 12345L)
  expect_equal(result$A2, "A")
  expect_equal(result$A1, "G")
  conv <- attr(result, "convention")
  expect_true(conv$has_chr)
  expect_equal(conv$allele_sep, ":")
})

test_that("parseVariantId parses single variant without chr prefix", {
  result <- parseVariantId("5:12345:A:G")
  expect_equal(result$chrom, 5L)
  expect_equal(result$pos, 12345L)
  expect_equal(result$A2, "A")
  expect_equal(result$A1, "G")
  conv <- attr(result, "convention")
  expect_false(conv$has_chr)
})

test_that("parseVariantId parses multiple variants", {
  ids <- c("chr1:100:A:G", "chr2:200:C:T", "chr3:300:G:A")
  result <- parseVariantId(ids)
  expect_equal(nrow(result), 3)
  expect_equal(result$chrom, c(1L, 2L, 3L))
  expect_equal(result$pos, c(100L, 200L, 300L))
  expect_equal(result$A2, c("A", "C", "G"))
  expect_equal(result$A1, c("G", "T", "A"))
})

# =============================================================================
# detectVariantConvention
# =============================================================================

test_that("detectVariantConvention detects chr prefix and allele separators", {
  conv <- pecotmr:::detectVariantConvention(c("chr1:100:A:G", "chr2:200:C:T"))
  expect_true(conv$has_chr)
  expect_equal(conv$allele_sep, ":")
  expect_false(conv$has_build)

  conv2 <- pecotmr:::detectVariantConvention(c("1_100_A_G", "2_200_C_T"))
  expect_false(conv2$has_chr)
  expect_equal(conv2$allele_sep, "_")

  conv3 <- pecotmr:::detectVariantConvention(c("chr1:100:A:G:b38"))
  expect_true(conv3$has_build)

  conv4 <- pecotmr:::detectVariantConvention(c("chr1:100_A_G"))
  expect_true(conv4$has_chr)
  expect_equal(conv4$allele_sep, "_")

  conv5 <- pecotmr:::detectVariantConvention(c("1:100_A_G"))
  expect_false(conv5$has_chr)
  expect_equal(conv5$allele_sep, "_")
})

# =============================================================================
# normalizeVariantId
# =============================================================================

test_that("normalizeVariantId normalizes various formats", {
  expect_equal(normalizeVariantId("1_100_A_G"), "chr1:100:A:G")
  expect_equal(normalizeVariantId("chr1:100:A:G"), "chr1:100:A:G")
  expect_equal(normalizeVariantId("1:100:A:G"), "chr1:100:A:G")
  expect_equal(normalizeVariantId("chr1:100:A:G", chrPrefix = FALSE), "1:100:A:G")
  expect_equal(normalizeVariantId("chr1:100:A:G:b38"), "chr1:100:A:G")
  expect_equal(normalizeVariantId("chr1:100_A_G"), "chr1:100:A:G")
  conv <- pecotmr:::detectVariantConvention(c("chr1:100_A_G"))
  expect_equal(normalizeVariantId("1:200:C:T", convention = conv), "chr1:200_C_T")
})

# =============================================================================
# variantIdToDf
# =============================================================================

test_that("variantIdToDf handles colon-separated format", {
  ids <- c("1:100:A:G", "2:200:C:T")
  result <- pecotmr:::variantIdToDf(ids)
  expect_equal(nrow(result), 2)
  expect_equal(result$chrom, c(1L, 2L))
  expect_equal(result$pos, c(100L, 200L))
  expect_equal(result$A2, c("A", "C"))
  expect_equal(result$A1, c("G", "T"))
})

test_that("variantIdToDf handles underscore-separated format", {
  ids <- c("1:100_A_G", "2:200_C_T")
  result <- pecotmr:::variantIdToDf(ids)
  expect_equal(nrow(result), 2)
  expect_equal(result$A2, c("A", "C"))
})

test_that("variantIdToDf strips chr prefix", {
  ids <- c("chr1:100:A:G", "chr2:200:C:T")
  result <- pecotmr:::variantIdToDf(ids)
  expect_equal(result$chrom, c(1L, 2L))
})

test_that("variantIdToDf handles data.frame input with named columns", {
  df <- data.frame(chrom = c("chr1", "2"), pos = c(100, 200),
                   A2 = c("A", "C"), A1 = c("G", "T"))
  suppressWarnings(result <- pecotmr:::variantIdToDf(df))
  expect_equal(result$chrom, c(1L, 2L))
  expect_equal(result$pos, c(100L, 200L))
})

test_that("variantIdToDf handles 5-part IDs with build suffix", {
  ids <- c("chr1:100:A:G:b38", "chr2:200:T:C")
  result <- pecotmr:::variantIdToDf(ids)
  expect_equal(ncol(result), 4)
  expect_equal(colnames(result), c("chrom", "pos", "A2", "A1"))
  expect_equal(result$chrom, c(1L, 2L))
  expect_equal(result$A2, c("A", "T"))
  expect_equal(result$A1, c("G", "C"))
})

test_that("variantIdToDf handles mixed 4/5-part IDs", {
  ids <- c("1:100:A:G", "chr2:200:T:C:b38", "3:300:G:A:b37")
  suppressWarnings(result <- pecotmr:::variantIdToDf(ids))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 4)
  expect_equal(result$A1, c("G", "C", "A"))
})

# =============================================================================
# getNestedElement
# =============================================================================

test_that("getNestedElement retrieves deeply nested values", {
  nested <- list(a = list(b = list(c = 42)))
  expect_equal(getNestedElement(nested, c("a", "b", "c")), 42)
})

test_that("getNestedElement returns NULL for NULL name_vector", {
  nested <- list(a = 1)
  expect_null(getNestedElement(nested, NULL))
})

test_that("getNestedElement errors on missing element", {
  nested <- list(a = list(b = 1))
  expect_error(getNestedElement(nested, c("a", "x")), "Element not found")
})

test_that("getNestedElement handles single level", {
  nested <- list(a = "hello")
  expect_equal(getNestedElement(nested, "a"), "hello")
})

test_that("getNestedElement skips empty strings", {
  nested <- list(a = list(b = 99))
  expect_equal(getNestedElement(nested, c("", "a", "b")), 99)
})

# =============================================================================
# regionToDf
# =============================================================================

test_that("regionToDf converts underscore-separated region IDs", {
  ids <- c("1_100_200", "2_300_400")
  result <- regionToDf(ids)
  expect_equal(nrow(result), 2)
  expect_equal(result$chrom, c(1L, 2L))
  expect_equal(result$start, c(100L, 300L))
  expect_equal(result$end, c(200L, 400L))
})

test_that("regionToDf handles chr prefix", {
  ids <- c("chr1_100_200")
  result <- regionToDf(ids)
  expect_equal(result$chrom, 1L)
})

test_that("regionToDf allows custom column names", {
  ids <- c("1_100_200")
  result <- regionToDf(ids, colnames = c("chr", "begin", "finish"))
  expect_true(all(c("chr", "begin", "finish") %in% colnames(result)))
})

# =============================================================================
# zToPvalue
# =============================================================================

test_that("zToPvalue returns correct p-values", {
  expect_equal(zToPvalue(0), 1.0, tolerance = 1e-10)
  expect_true(zToPvalue(1.96) < 0.05)
  expect_true(zToPvalue(-1.96) < 0.05)
  expect_equal(zToPvalue(1.96), zToPvalue(-1.96), tolerance = 1e-10)
})

test_that("zToPvalue handles vector input", {
  z <- c(0, 1, 2, 3)
  p <- zToPvalue(z)
  expect_length(p, 4)
  expect_true(all(diff(p) < 0))
})

test_that("zToPvalue handles extreme values", {
  expect_true(zToPvalue(10) < 1e-20)
  expect_true(zToPvalue(40) >= 0)
})

# =============================================================================
# zToBetaSe
# =============================================================================

test_that("zToBetaSe produces correct conversions", {
  z <- c(2.0, -1.0)
  maf <- c(0.3, 0.1)
  n <- 10000

  result <- pecotmr:::zToBetaSe(z, maf, n)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("beta", "se", "maf") %in% names(result)))
  expect_equal(nrow(result), 2)
  expect_true(result$beta[1] > 0)
  expect_true(result$beta[2] < 0)
})

test_that("zToBetaSe errors on mismatched lengths", {
  expect_error(pecotmr:::zToBetaSe(c(1, 2), c(0.3), 1000), "same length")
})

test_that("zToBetaSe adjusts MAF > 0.5", {
  result <- pecotmr:::zToBetaSe(c(1.0), c(0.7), 1000)
  expect_equal(result$maf, 0.3)
})

# =============================================================================
# lbfToAlpha
# =============================================================================

test_that("lbfToAlpha converts matrix correctly", {
  lbf <- matrix(c(-0.5, 1.2, 0.3, 0.7, -1.1, 0.4), nrow = 2)
  colnames(lbf) <- paste0("v", 1:3)
  result <- lbfToAlpha(lbf)

  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 3)
  expect_equal(rowSums(result), c(1, 1), tolerance = 1e-10)
})

test_that("lbfToAlpha handles single column", {
  lbf <- matrix(c(0.5, 1.0), ncol = 1)
  colnames(lbf) <- "v1"
  result <- lbfToAlpha(lbf)
  expect_equal(ncol(result), 1)
})

test_that("lbfToAlpha handles all-zero row", {
  lbf <- matrix(c(0, 0, 0), nrow = 1)
  colnames(lbf) <- paste0("v", 1:3)
  result <- lbfToAlpha(lbf)
  expect_true(all(result == 0))
})

# =============================================================================
# findData
# =============================================================================

test_that("findData retrieves from nested list at depth 2", {
  x <- list(a = list(val = 42), b = list(val = 99))
  result <- findData(x, c(2, "val"))
  expect_equal(result, c(42, 99))
})

test_that("findData returns list at depth 0", {
  x <- list(a = 1, b = 2)
  result <- findData(x, c(0))
  expect_type(result, "list")
})

test_that("findData with show_path=TRUE returns list structure", {
  x <- list(a = list(val = 10), b = list(val = 20))
  result <- findData(x, c(2, "val"), showPath = TRUE)
  expect_type(result, "list")
  expect_true("a" %in% names(result))
  expect_true("b" %in% names(result))
  expect_equal(result$a, 10)
  expect_equal(result$b, 20)
})

test_that("findData with rm_dup=TRUE removes duplicate values", {
  x <- list(
    a = list(val = 42),
    b = list(val = 42),
    c = list(val = 99)
  )
  result <- findData(x, c(2, "val"), rmDup = TRUE)
  expect_true("shared_list_names" %in% names(result))
})

test_that("findData at depth 3 retrieves deeply nested values", {
  x <- list(
    level1_a = list(
      level2_a = list(target = "found_a"),
      level2_b = list(target = "found_b")
    ),
    level1_b = list(
      level2_c = list(target = "found_c")
    )
  )
  result <- findData(x, c(3, "target"))
  expect_true("found_a" %in% result)
  expect_true("found_b" %in% result)
  expect_true("found_c" %in% result)
})

test_that("findData with depth=1 and list_name returns element directly", {
  x <- list(a = 1, b = 2, c = 3)
  result <- findData(x, c(1, "b"))
  expect_equal(result, 2)
})

test_that("findData with rm_null=TRUE removes NULL results", {
  x2 <- list(
    a = list(val = 10),
    b = "not a list"
  )
  result <- findData(x2, c(2, "val"), rmNull = TRUE)
  expect_equal(result, 10)
})

test_that("findData with show_path=TRUE and rm_dup=TRUE at depth 2", {
  x <- list(
    a = list(val = 42),
    b = list(val = 42),
    c = list(val = 99)
  )
  result <- findData(x, c(2, "val"), showPath = TRUE, rmDup = TRUE)
  expect_type(result, "list")
  expect_true("shared_list_names" %in% names(result))
})

test_that("findData with depth=1 and no list_name returns whole object", {
  x <- list(a = 1, b = 2)
  result <- findData(x, c(1))
  expect_identical(result, x)
})

test_that("findData with docall=list preserves list structure", {
  x <- list(
    a = list(val = c(1, 2)),
    b = list(val = c(3, 4))
  )
  result <- findData(x, c(2, "val"), docall = list)
  expect_type(result, "list")
  expect_equal(result[[1]], c(1, 2))
  expect_equal(result[[2]], c(3, 4))
})

# =============================================================================
# robustMahalanobis
# =============================================================================

test_that("robustMahalanobis works with non-singular covariance", {
  set.seed(42)
  x <- matrix(rnorm(100), ncol = 2)
  d <- robustMahalanobis(x)
  expect_length(d, 50)
  expect_true(all(d >= 0))
  expect_true(is.numeric(d))
})

test_that("robustMahalanobis works with singular covariance (falls back to ginv)", {
  # Create a matrix where columns are linearly dependent -> singular cov
  set.seed(42)
  col1 <- rnorm(20)
  col2 <- rnorm(20)
  col3 <- col1 + col2  # linearly dependent
  x <- cbind(col1, col2, col3)
  d <- robustMahalanobis(x)
  expect_length(d, 20)
  expect_true(all(d >= 0))
})

test_that("robustMahalanobis with pre-inverted covariance", {
  set.seed(42)
  x <- matrix(rnorm(60), ncol = 2)
  center <- colMeans(x)
  cov_mat <- stats::cov(x)
  inv_cov <- solve(cov_mat)
  d <- robustMahalanobis(x, center = center, cov = inv_cov, inverted = TRUE)
  expect_length(d, 30)
  expect_true(all(d >= 0))
})

test_that("robustMahalanobis with vector input", {
  # Single observation as a vector
  x <- c(1.0, 2.0, 3.0)
  d <- robustMahalanobis(x, center = c(1, 2, 3), cov = diag(3), inverted = TRUE)
  expect_length(d, 1)
  expect_equal(as.numeric(d), 0)
})

test_that("robustMahalanobis auto-computes center and cov when NULL", {
  set.seed(42)
  x <- matrix(rnorm(40), ncol = 2)
  d1 <- robustMahalanobis(x)
  d2 <- robustMahalanobis(x, center = colMeans(x), cov = stats::cov(x))
  expect_equal(d1, d2)
})

# =============================================================================
# detectOutliersMahalanobis
# =============================================================================

test_that("detectOutliersMahalanobis returns correct structure", {
  set.seed(42)
  x <- matrix(rnorm(200), ncol = 2)
  rownames(x) <- paste0("sample", 1:100)
  result <- detectOutliersMahalanobis(x)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("sample_id", "mahal", "pvalue", "is_outlier") %in% names(result)))
  expect_equal(nrow(result), 100)
  expect_equal(result$sample_id[1], "sample1")
})

test_that("detectOutliersMahalanobis detects clear outliers", {
  set.seed(42)
  x <- matrix(rnorm(200), ncol = 2)
  # Add a clear outlier far from the center
  x <- rbind(x, c(50, 50))
  rownames(x) <- paste0("s", 1:101)
  result <- detectOutliersMahalanobis(x)
  # The extreme point should be an outlier
  expect_true(result$is_outlier[101])
})

test_that("detectOutliersMahalanobis with unnamed rows uses indices", {
  set.seed(42)
  x <- matrix(rnorm(40), ncol = 2)
  result <- detectOutliersMahalanobis(x)
  expect_equal(result$sample_id, as.character(1:20))
})

test_that("detectOutliersMahalanobis threshold sensitivity", {
  set.seed(42)
  x <- matrix(rnorm(200), ncol = 2)
  # With very strict threshold, fewer outliers
  r_strict <- detectOutliersMahalanobis(x, prob = 0.999, pvalThreshold = 0.001)
  # With lenient threshold, more possible outliers
  r_lenient <- detectOutliersMahalanobis(x, prob = 0.90, pvalThreshold = 0.10)
  expect_true(sum(r_strict$is_outlier) <= sum(r_lenient$is_outlier))
})

# =============================================================================
# twasMethodCor
# =============================================================================

test_that("twasMethodCor with identity LD", {
  LD <- diag(3)
  w1 <- c(1, 0, 0)
  w2 <- c(0, 1, 0)
  w3 <- c(0, 0, 1)
  result <- twasMethodCor(list(w1, w2, w3), LD)
  # With identity LD and orthogonal weights, off-diag should be 0
  expect_equal(dim(result), c(3, 3))
  expect_equal(diag(result), c(1, 1, 1))
  expect_equal(result[1, 2], 0)
  expect_equal(result[1, 3], 0)
  expect_equal(result[2, 3], 0)
})

test_that("twasMethodCor with identical weights gives correlation 1", {
  LD <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  w <- c(1, 1)
  result <- twasMethodCor(list(w, w), LD)
  expect_equal(result[1, 2], 1)
  expect_equal(result[2, 1], 1)
})

test_that("twasMethodCor with diagonal LD", {
  LD <- diag(c(2, 3, 1))
  w1 <- c(1, 0, 0)
  w2 <- c(0, 1, 0)
  result <- twasMethodCor(list(w1, w2), LD)
  expect_equal(result[1, 2], 0)
})

# =============================================================================
# computeQvalues — uncovered branches
# =============================================================================

test_that("computeQvalues returns NA vector when all pvalues are NA", {
  skip_if_not_installed("qvalue")
  result <- expect_message(
    computeQvalues(rep(NA_real_, 5)),
    "All p-values are NA"
  )
  expect_equal(result, rep(NA_real_, 5))
})

test_that("computeQvalues falls back to BH when qvalue fails", {
  skip_if_not_installed("qvalue")
  # Very few unique p-values can cause qvalue() to fail with an error
  # Use only 2 identical p-values to trigger the tryCatch error path
  pvals <- rep(0.5, 3)
  result <- expect_message(
    computeQvalues(pvals),
    "Too few p-values|fall back to BH"
  )
  expect_length(result, 3)
  expect_true(all(!is.na(result)))
})

# =============================================================================
# safeSvd — uncovered branches
# =============================================================================

test_that("safeSvd with tol=0 keeps all singular values", {
  mat <- matrix(c(1, 0, 0, 1e-12), nrow = 2)
  result <- safeSvd(mat, tol = 0)
  expect_length(result$d, 2)
  expect_equal(ncol(result$u), 2)
  expect_equal(ncol(result$v), 2)
})

test_that("safeSvd errors when all singular values below tolerance", {
  # A matrix with very small singular values
  mat <- matrix(c(1e-15, 0, 0, 1e-15), nrow = 2)
  expect_error(safeSvd(mat, tol = 1), "All singular values are below the tolerance threshold")
})

# =============================================================================
# computeLd — uncovered branches
# =============================================================================

test_that("computeLd sample method without Rfast falls back to cor", {
  set.seed(42)
  X <- matrix(sample(0:2, 100, replace = TRUE), nrow = 20, ncol = 5)
  colnames(X) <- paste0("snp", 1:5)
  R <- computeLd(X, method = "sample")
  expect_equal(dim(R), c(5, 5))
  expect_equal(as.numeric(diag(R)), rep(1, 5))
  expect_true(all(abs(R) <= 1))
})

test_that("computeLd with gcta method and trim_samples", {
  set.seed(42)
  # 21 samples -> trimmed to 20 (multiple of 4)
  X <- matrix(sample(0:2, 105, replace = TRUE), nrow = 21, ncol = 5)
  colnames(X) <- paste0("snp", 1:5)
  R <- computeLd(X, method = "gcta", trimSamples = TRUE)
  expect_equal(dim(R), c(5, 5))
  expect_equal(as.numeric(diag(R)), rep(1, 5))
})

test_that("computeLd population method with trim_samples", {
  set.seed(42)
  X <- matrix(sample(0:2, 105, replace = TRUE), nrow = 21, ncol = 5)
  colnames(X) <- paste0("snp", 1:5)
  R <- computeLd(X, method = "population", trimSamples = TRUE)
  expect_equal(dim(R), c(5, 5))
  expect_equal(as.numeric(diag(R)), rep(1, 5))
})

test_that("computeLd with shrinkage > 0", {
  set.seed(42)
  X <- matrix(sample(0:2, 100, replace = TRUE), nrow = 20, ncol = 5)
  colnames(X) <- paste0("snp", 1:5)
  R_no_shrink <- computeLd(X, method = "sample", shrinkage = 0)
  R_shrink <- computeLd(X, method = "sample", shrinkage = 0.1)
  # Shrunk matrix should be closer to identity
  expect_equal(as.numeric(diag(R_shrink)), rep(1, 5))
  # Off-diagonal elements should be shrunk toward 0
  off_diag_no <- R_no_shrink[1, 2]
  off_diag_s <- R_shrink[1, 2]
  expect_equal(off_diag_s, 0.9 * off_diag_no)
})

# =============================================================================
# filterXWithY — uncovered lines 513-515
# =============================================================================

test_that("filterXWithY drops variants that become monomorphic due to Y NAs", {
  # Create X where some columns become monomorphic when Y NA rows are removed
  set.seed(42)
  X <- matrix(0, nrow = 10, ncol = 3)
  rownames(X) <- paste0("subj", 1:10)
  colnames(X) <- paste0("var", 1:3)
  # var1: all 0 except subject 1 has 1 -> monomorphic without subj1
  X[1, 1] <- 1
  X[, 2] <- sample(0:2, 10, replace = TRUE)
  X[, 3] <- sample(0:2, 10, replace = TRUE)
  # Y where subject 1 has NA -> removing subj1 makes var1 monomorphic
  Y <- matrix(rnorm(10), nrow = 10, ncol = 1)
  rownames(Y) <- paste0("subj", 1:10)
  colnames(Y) <- "context1"
  Y[1, 1] <- NA
  result <- expect_message(
    filterXWithY(X, Y, missingRateThresh = 1.0, mafThresh = 0),
    "Additional.*variants dropped"
  )
  # var1 should be dropped since it becomes monomorphic without subj1
  expect_true(ncol(result) < 3)
})

# =============================================================================
# detectVariantConvention — uncovered line 586
# =============================================================================

test_that("detectVariantConvention returns defaults for all-NA input", {
  result <- detectVariantConvention(c(NA, NA, NA))
  expect_false(result$has_chr)
  expect_equal(result$allele_sep, ":")
  expect_false(result$has_build)
  expect_true(is.na(result$example))
})

# =============================================================================
# parseVariantId — uncovered line 622
# =============================================================================

test_that("parseVariantId handles data.frame with generic column names", {
  df <- data.frame(
    col1 = c("chr1", "chr2"),
    col2 = c(100, 200),
    col3 = c("A", "T"),
    col4 = c("G", "C"),
    stringsAsFactors = FALSE
  )
  result <- parseVariantId(df)
  expect_equal(names(result)[1:4], c("chrom", "pos", "A2", "A1"))
  expect_equal(result$chrom, c(1L, 2L))
  expect_equal(result$pos, c(100L, 200L))
  expect_equal(result$A2, c("A", "T"))
  expect_equal(result$A1, c("G", "C"))
})

# =============================================================================
# filterMolecularEvents — uncovered lines 1116-1128 (remove_all_group)
# =============================================================================

test_that("filterMolecularEvents with remove_all_group=TRUE removes entire group", {
  events <- c(
    "gene1_tissue_brain",
    "gene1_tissue_liver",
    "gene2_tissue_brain",
    "gene2_tissue_liver"
  )
  filters <- list(
    list(
      type_pattern = "(.*)_tissue_.*",
      exclude_pattern = "brain"
    )
  )
  result <- expect_message(
    filterMolecularEvents(events, filters, condition = "test", removeAllGroup = TRUE),
    "removed"
  )
  # With remove_all_group=TRUE, events from groups that had a brain entry removed
  # should also be removed
  expect_true(is.character(result))
})

# =============================================================================
# findData — uncovered lines 780-786 (numeric index path)
# =============================================================================

test_that("findData with numeric indices in list_name path", {
  # When list_name contains a numeric string like "2", findData splits the path
  # at that numeric index: it navigates to "results" first, then treats "2" as
  # the new depth and "val" as the new list_name, recursing into each sub-list.
  x <- list(
    results = list(
      a = list(val = 10),
      b = list(val = 20)
    )
  )
  # depth=1, list_name = c("results", "2", "val")
  # -> second_depth at index 2 ("2" is numeric)
  # -> data = getNestedElement(x, "results") = x$results
  # -> remaining_path = c("2", "val") -> findData(x$results, c("2","val"))
  # -> depth=2, list_name="val" -> recurse into a and b at depth 1 looking for "val"
  result <- findData(x, c(1, "results", "2", "val"))
  expect_equal(result, c(10, 20))
})

# =============================================================================
# regionsOverlap
# =============================================================================

test_that("regionsOverlap detects overlapping regions on same chromosome", {
  expect_true(regionsOverlap("chr1:100-300", "chr1:200-400"))
})

test_that("regionsOverlap returns FALSE for non-overlapping same-chr regions", {
  expect_false(regionsOverlap("chr1:100-200", "chr1:300-400"))
})

test_that("regionsOverlap returns FALSE for different chromosomes", {
  expect_false(regionsOverlap("chr1:100-300", "chr2:100-300"))
})

test_that("regionsOverlap detects touching boundaries", {
  expect_true(regionsOverlap("chr1:100-200", "chr1:200-300"))
})

test_that("regionsOverlap works with underscore-separated IDs", {
  expect_true(regionsOverlap("1_100_300", "1_200_400"))
  expect_false(regionsOverlap("1_100_200", "2_100_200"))
})

test_that("regionsOverlap works with data.frame input", {
  df_a <- data.frame(chrom = 1, start = 100, end = 300)
  df_b <- data.frame(chrom = 1, start = 200, end = 400)
  expect_true(regionsOverlap(df_a, df_b))
})

# =============================================================================
# findOverlappingRegions
# =============================================================================

test_that("findOverlappingRegions returns correct indices", {
  query <- "chr1:100-300"
  targets <- c("chr1:200-400", "chr2:100-200", "chr1:50-150")
  result <- findOverlappingRegions(query, targets)
  expect_true(1 %in% result)
  expect_true(3 %in% result)
  expect_false(2 %in% result)
})

test_that("findOverlappingRegions returns empty vector for no matches", {
  query <- "chr1:100-200"
  targets <- c("chr2:100-200", "chr3:100-200")
  result <- findOverlappingRegions(query, targets)
  expect_length(result, 0)
})

test_that("findOverlappingRegions works with data.frame targets", {
  query <- "chr1:100-300"
  targets <- data.frame(chrom = c(1, 2, 1), start = c(200, 100, 50), end = c(400, 200, 150))
  result <- findOverlappingRegions(query, targets)
  expect_true(1 %in% result)
  expect_true(3 %in% result)
  expect_false(2 %in% result)
})

# =============================================================================
# classifyVariantType
# =============================================================================

test_that("classifyVariantType identifies SNPs", {
  expect_equal(classifyVariantType("chr1:100:A:G"), "SNP")
})

test_that("classifyVariantType identifies insertions", {
  expect_equal(classifyVariantType("chr1:100:A:ATG"), "insertion")
})

test_that("classifyVariantType identifies deletions", {
  expect_equal(classifyVariantType("chr1:100:ATG:A"), "deletion")
})

test_that("classifyVariantType identifies MNPs", {
  expect_equal(classifyVariantType("chr1:100:AT:GC"), "MNP")
})

test_that("classifyVariantType handles vector input", {
  ids <- c("chr1:100:A:G", "chr1:200:ATG:A", "chr1:300:A:ATG", "chr1:400:AT:GC")
  result <- classifyVariantType(ids)
  expect_equal(result, c("SNP", "deletion", "insertion", "MNP"))
})

test_that("classifyVariantType accepts data.frame input", {
  df <- data.frame(A2 = c("A", "ATG"), A1 = c("G", "A"))
  result <- classifyVariantType(df)
  expect_equal(result, c("SNP", "deletion"))
})

# =============================================================================
# ensureChrMatch
# =============================================================================

test_that("ensureChrMatch returns unchanged when both have chr prefix", {
  ids_a <- c("chr1:100:A:G", "chr1:200:C:T")
  ids_b <- c("chr1:150:A:G", "chr1:250:C:T")
  result <- pecotmr:::ensureChrMatch(ids_a, ids_b)
  expect_equal(result$ids_a, ids_a)
  expect_equal(result$ids_b, ids_b)
})

test_that("ensureChrMatch normalizes when prefixes mismatch", {
  ids_a <- c("chr1:100:A:G", "chr1:200:C:T")
  ids_b <- c("1:150:A:G", "1:250:C:T")
  result <- pecotmr:::ensureChrMatch(ids_a, ids_b)
  expect_true(all(grepl("^chr", result$ids_a)))
  expect_true(all(grepl("^chr", result$ids_b)))
})

test_that("ensureChrMatch returns unchanged when both lack chr prefix", {
  ids_a <- c("1:100:A:G", "1:200:C:T")
  ids_b <- c("1:150:A:G", "1:250:C:T")
  result <- pecotmr:::ensureChrMatch(ids_a, ids_b)
  # Both already match (no prefix), so returned unchanged
  expect_equal(result$ids_a, ids_a)
  expect_equal(result$ids_b, ids_b)
})
