context("SS-TWAS: weights, pipeline, and omnibus combination")

# =============================================================================
# TwasWeights S4 class with standardized slot
# =============================================================================

test_that("TwasWeights accepts standardized = TRUE", {
  wt <- TwasWeights(
    weights = list(method1 = matrix(1:5, ncol = 1)),
    variantIds = paste0("v", 1:5),
    standardized = TRUE
  )
  expect_true(wt@standardized)
  expect_equal(length(wt@methods), 1L)
})

test_that("TwasWeights defaults to standardized = FALSE", {
  wt <- TwasWeights(
    weights = list(method1 = matrix(1:5, ncol = 1)),
    variantIds = paste0("v", 1:5)
  )
  expect_false(wt@standardized)
})

test_that("TwasWeights show method includes standardized", {
  wt <- TwasWeights(
    weights = list(m1 = matrix(0, nrow = 3, ncol = 1)),
    variantIds = paste0("v", 1:3),
    standardized = TRUE
  )
  out <- capture.output(show(wt))
  expect_true(any(grepl("Standardized: TRUE", out)))
})

# =============================================================================
# SuSiE-RSS weight extraction
# =============================================================================

test_that(".susie_rss_extract_weights returns correct-length vector", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 20
  n <- 500
  R <- diag(p)
  z <- rnorm(p)
  w <- pecotmr:::.susieRssExtractWeights(
    fit = NULL, z = z, R = R, n = n,
    requiredFields = c("alpha", "mu", "X_column_scale_factors"),
    fitArgs = list(L = 5)
  )
  expect_equal(length(w), p)
  expect_true(all(is.finite(w)))
})

test_that("susieRssWeights follows (stat, LD) convention", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 20
  n <- 500
  R <- diag(p)
  z <- rnorm(p)
  stat <- list(b = z / sqrt(n), cor = z / sqrt(n), z = z, n = rep(n, p))
  w <- susieRssWeights(stat, R, methodArgs = list(L = 5))
  expect_equal(length(w), p)
  expect_true(all(is.finite(w)))
})

test_that("susieRssWeights retains fit when retainFit = TRUE", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 20
  n <- 500
  R <- diag(p)
  z <- rnorm(p)
  stat <- list(b = z / sqrt(n), cor = z / sqrt(n), z = z, n = rep(n, p))
  w <- susieRssWeights(stat, R, retainFit = TRUE, methodArgs = list(L = 5))
  expect_false(is.null(attr(w, "fit")))
})

test_that("susieInfRssWeights works", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 20
  n <- 500
  R <- diag(p)
  z <- rnorm(p)
  stat <- list(b = z / sqrt(n), cor = z / sqrt(n), z = z, n = rep(n, p))
  w <- susieInfRssWeights(stat, R, methodArgs = list(L = 5))
  expect_equal(length(w), p)
  expect_true(all(is.finite(w)))
})

# =============================================================================
# Two-stage SuSiE-RSS fitting
# =============================================================================

test_that("fitSusieInfThenSusieRss returns two fits", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 20
  n <- 500
  R <- diag(p)
  z <- rnorm(p)
  fits <- fitSusieInfThenSusieRss(z, R, n, args = list(L = 5))
  expect_true(is.list(fits))
  expect_true("susie" %in% names(fits))
  expect_true("susie_inf" %in% names(fits))
  expect_true("susie_inf" %in% class(fits$susie_inf))
  expect_true("susie_rss" %in% class(fits$susie))
})

# =============================================================================
# twasAnalysis omnibus combination
# =============================================================================

test_that("twasAnalysis adds omnibus when combine_if_no_cv = TRUE", {
  set.seed(42)
  p <- 10
  R <- diag(p)
  rownames(R) <- colnames(R) <- paste0("v", 1:p)
  weights_matrix <- matrix(rnorm(p * 3), ncol = 3)
  rownames(weights_matrix) <- paste0("v", 1:p)
  colnames(weights_matrix) <- c("m1", "m2", "m3")
  gwas_db <- data.frame(
    variant_id = paste0("v", 1:p),
    z = rnorm(p)
  )

  result <- twasAnalysis(
    weights_matrix, gwas_db, ldMatrix = R,
    extractVariantsObjs = paste0("v", 1:p),
    combineIfNoCv = TRUE
  )
  expect_true("omnibus" %in% names(result))
  expect_true(!is.null(result$omnibus$pval))
})

test_that("twasAnalysis skips omnibus when combine_if_no_cv = FALSE", {
  set.seed(42)
  p <- 10
  R <- diag(p)
  rownames(R) <- colnames(R) <- paste0("v", 1:p)
  weights_matrix <- matrix(rnorm(p * 3), ncol = 3)
  rownames(weights_matrix) <- paste0("v", 1:p)
  colnames(weights_matrix) <- c("m1", "m2", "m3")
  gwas_db <- data.frame(
    variant_id = paste0("v", 1:p),
    z = rnorm(p)
  )

  result <- twasAnalysis(
    weights_matrix, gwas_db, ldMatrix = R,
    extractVariantsObjs = paste0("v", 1:p),
    combineIfNoCv = FALSE
  )
  expect_false("omnibus" %in% names(result))
})

# =============================================================================
# twasWeightsSumstatPipeline end-to-end
# =============================================================================

test_that("twasWeightsSumstatPipeline produces TwasWeights with standardized = TRUE", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 30
  n <- 1000
  R <- diag(p)
  rownames(R) <- colnames(R) <- paste0("1:", 1000 + seq_len(p), ":A:T")
  z <- rnorm(p, sd = 2)
  sumstats <- data.frame(
    variant_id = rownames(R),
    chrom = "1",
    pos = 1000 + seq_len(p),
    A1 = "T",
    A2 = "A",
    z = z
  )

  ref_panel <- data.frame(
    chrom = "1",
    pos = 1000 + seq_len(p),
    variant_id = rownames(R),
    A1 = "T",
    A2 = "A",
    stringsAsFactors = FALSE
  )
  variants_gr <- pecotmr:::.refPanelToGranges(ref_panel)
  block_metadata <- S4Vectors::DataFrame(
    region = paste0("chr1:", 1001, "-", 1000 + p),
    start = 1001L, end = as.integer(1000 + p), chrom = "chr1",
    start_idx = 1L, end_idx = as.integer(p), size = as.integer(p)
  )
  ld_data <- new("LdData",
    correlation = R,
    genotypeHandle = NULL,
    variants = variants_gr,
    snpIdx = seq_len(p),
    blockMetadata = block_metadata
  )

  result <- twasWeightsSumstatPipeline(
    sumstats = sumstats,
    ldData = ld_data,
    n = n,
    methods = list(susie_rss = list(L = 5)),
    pThresholds = c(0.05),
    checkLdMethod = NULL,
    verbose = 0
  )

  expect_false(is.null(result$twas_weights))
  expect_true(is(result$twas_weights, "TwasWeights"))
  expect_true(result$twas_weights@standardized)
  expect_true(length(result$twas_weights@variantIds) > 0)
  expect_false(result$qc_summary$skipped)
})
