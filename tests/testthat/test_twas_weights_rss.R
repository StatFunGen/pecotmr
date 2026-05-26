context("SS-TWAS: weights, pipeline, and omnibus combination")

# =============================================================================
# TWASWeights S4 class with standardized slot
# =============================================================================

test_that("TWASWeights accepts standardized = TRUE", {
  wt <- TWASWeights(
    weights = list(method1 = matrix(1:5, ncol = 1)),
    variant_ids = paste0("v", 1:5),
    standardized = TRUE
  )
  expect_true(wt@standardized)
  expect_equal(length(wt@methods), 1L)
})

test_that("TWASWeights defaults to standardized = FALSE", {
  wt <- TWASWeights(
    weights = list(method1 = matrix(1:5, ncol = 1)),
    variant_ids = paste0("v", 1:5)
  )
  expect_false(wt@standardized)
})

test_that("TWASWeights show method includes standardized", {
  wt <- TWASWeights(
    weights = list(m1 = matrix(0, nrow = 3, ncol = 1)),
    variant_ids = paste0("v", 1:3),
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
  w <- pecotmr:::.susie_rss_extract_weights(
    fit = NULL, z = z, R = R, n = n,
    required_fields = c("alpha", "mu", "X_column_scale_factors"),
    fit_args = list(L = 5)
  )
  expect_equal(length(w), p)
  expect_true(all(is.finite(w)))
})

test_that("susie_rss_weights follows (stat, LD) convention", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 20
  n <- 500
  R <- diag(p)
  z <- rnorm(p)
  stat <- list(b = z / sqrt(n), cor = z / sqrt(n), z = z, n = rep(n, p))
  w <- susie_rss_weights(stat, R, method_args = list(L = 5))
  expect_equal(length(w), p)
  expect_true(all(is.finite(w)))
})

test_that("susie_rss_weights retains fit when retain_fit = TRUE", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 20
  n <- 500
  R <- diag(p)
  z <- rnorm(p)
  stat <- list(b = z / sqrt(n), cor = z / sqrt(n), z = z, n = rep(n, p))
  w <- susie_rss_weights(stat, R, retain_fit = TRUE, method_args = list(L = 5))
  expect_false(is.null(attr(w, "fit")))
})

test_that("susie_inf_rss_weights works", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 20
  n <- 500
  R <- diag(p)
  z <- rnorm(p)
  stat <- list(b = z / sqrt(n), cor = z / sqrt(n), z = z, n = rep(n, p))
  w <- susie_inf_rss_weights(stat, R, method_args = list(L = 5))
  expect_equal(length(w), p)
  expect_true(all(is.finite(w)))
})

# =============================================================================
# Two-stage SuSiE-RSS fitting
# =============================================================================

test_that("fit_susie_inf_then_susie_rss returns two fits", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 20
  n <- 500
  R <- diag(p)
  z <- rnorm(p)
  fits <- fit_susie_inf_then_susie_rss(z, R, n, args = list(L = 5))
  expect_true(is.list(fits))
  expect_true("susie" %in% names(fits))
  expect_true("susie_inf" %in% names(fits))
  expect_true("susie_inf" %in% class(fits$susie_inf))
  expect_true("susie_rss" %in% class(fits$susie))
})

# =============================================================================
# twas_analysis omnibus combination
# =============================================================================

test_that("twas_analysis adds omnibus when combine_if_no_cv = TRUE", {
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

  result <- twas_analysis(
    weights_matrix, gwas_db, LD_matrix = R,
    extract_variants_objs = paste0("v", 1:p),
    combine_if_no_cv = TRUE
  )
  expect_true("omnibus" %in% names(result))
  expect_true(!is.null(result$omnibus$pval))
})

test_that("twas_analysis skips omnibus when combine_if_no_cv = FALSE", {
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

  result <- twas_analysis(
    weights_matrix, gwas_db, LD_matrix = R,
    extract_variants_objs = paste0("v", 1:p),
    combine_if_no_cv = FALSE
  )
  expect_false("omnibus" %in% names(result))
})

# =============================================================================
# twas_weights_sumstat_pipeline end-to-end
# =============================================================================

test_that("twas_weights_sumstat_pipeline produces TWASWeights with standardized = TRUE", {
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
  variants_gr <- pecotmr:::.ref_panel_to_granges(ref_panel)
  block_metadata <- S4Vectors::DataFrame(
    region = paste0("chr1:", 1001, "-", 1000 + p),
    start = 1001L, end = as.integer(1000 + p), chrom = "chr1",
    start_idx = 1L, end_idx = as.integer(p), size = as.integer(p)
  )
  ld_data <- new("LDData",
    correlation = R,
    genotype_handle = NULL,
    variants = variants_gr,
    snp_idx = seq_len(p),
    block_metadata = block_metadata
  )

  result <- twas_weights_sumstat_pipeline(
    sumstats = sumstats,
    LD_data = ld_data,
    n = n,
    methods = list(susie_rss = list(L = 5)),
    p_thresholds = c(0.05),
    check_ld_method = NULL,
    verbose = 0
  )

  expect_false(is.null(result$twas_weights))
  expect_true(is(result$twas_weights, "TWASWeights"))
  expect_true(result$twas_weights@standardized)
  expect_true(length(result$twas_weights@variant_ids) > 0)
  expect_false(result$qc_summary$skipped)
})
