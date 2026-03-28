context("sumstats_qc")

# ===========================================================================
# Helper: build matching sumstats and LD_data
# ===========================================================================
make_test_sumstats_ld <- function(n_variants = 5, chrom_val = 1, with_indels = FALSE) {
  set.seed(42)
  positions <- seq_len(n_variants) * 100

  if (with_indels) {
    a1 <- c(rep("G", n_variants - 1), "ACGT")
    a2 <- c(rep("A", n_variants - 1), "A")
  } else {
    a1 <- rep("G", n_variants)
    a2 <- rep("A", n_variants)
  }

  variant_ids <- paste0(chrom_val, ":", positions, ":", a2, ":", a1)

  sumstats <- data.frame(
    chrom      = rep(chrom_val, n_variants),
    pos        = positions,
    A1         = a1,
    A2         = a2,
    beta       = rnorm(n_variants, 0, 0.5),
    se         = runif(n_variants, 0.05, 0.2),
    z          = rnorm(n_variants, 0, 2),
    stringsAsFactors = FALSE
  )

  LD_mat <- diag(n_variants) + matrix(0.01, n_variants, n_variants)
  diag(LD_mat) <- 1
  rownames(LD_mat) <- colnames(LD_mat) <- variant_ids

  ref_panel <- data.frame(
    chrom = rep(chrom_val, n_variants),
    pos   = positions,
    A2    = a2,
    A1    = a1,
    stringsAsFactors = FALSE
  )

  LD_data <- list(
    LD_variants = ref_panel,
    LD_matrix   = LD_mat,
    ref_panel            = ref_panel
  )

  list(sumstats = sumstats, LD_data = LD_data, variant_ids = variant_ids)
}

# ===========================================================================
# rss_basic_qc
# ===========================================================================

test_that("rss_basic_qc requires correct columns", {
  sumstats <- data.frame(beta = 1, se = 0.5)
  LD_data <- list(LD_variants = data.frame())
  expect_error(rss_basic_qc(sumstats, LD_data), "Missing columns")
})

test_that("rss_basic_qc processes matching variants correctly", {
  variant_ids <- c("1:100:A:G", "1:200:C:T", "1:300:G:A")

  sumstats <- data.frame(
    chrom = c(1, 1, 1),
    pos = c(100, 200, 300),
    A1 = c("G", "T", "A"),
    A2 = c("A", "C", "G"),
    beta = c(0.5, -0.3, 0.1),
    se = c(0.1, 0.15, 0.2),
    z = c(5.0, -2.0, 0.5),
    stringsAsFactors = FALSE
  )

  LD_mat <- diag(3)
  rownames(LD_mat) <- colnames(LD_mat) <- variant_ids

  ref_panel <- data.frame(
    chrom = c(1, 1, 1),
    pos = c(100, 200, 300),
    A2 = c("A", "C", "G"),
    A1 = c("G", "T", "A"),
    stringsAsFactors = FALSE
  )

  LD_data <- list(
    LD_variants = ref_panel,
    LD_matrix = LD_mat
  )

  result <- rss_basic_qc(sumstats, LD_data)
  expect_type(result, "list")
  expect_true("sumstats" %in% names(result))
  expect_true("LD_mat" %in% names(result))
})

test_that("rss_basic_qc skips variants in specified region", {
  td <- make_test_sumstats_ld(n_variants = 5)

  result <- rss_basic_qc(td$sumstats, td$LD_data, skip_region = "1:150-350")

  expect_type(result, "list")
  expect_true("sumstats" %in% names(result))
  expect_true("LD_mat" %in% names(result))
  remaining_pos <- result$sumstats$pos
  expect_false(200 %in% remaining_pos)
  expect_false(300 %in% remaining_pos)
})

test_that("rss_basic_qc with skip_region preserves non-skipped variants", {
  td <- make_test_sumstats_ld(n_variants = 5)
  result <- rss_basic_qc(td$sumstats, td$LD_data, skip_region = "1:150-250")
  remaining_pos <- result$sumstats$pos
  expect_false(200 %in% remaining_pos)
  expect_true(100 %in% remaining_pos)
  expect_true(300 %in% remaining_pos)
})

test_that("rss_basic_qc with remove_indels=TRUE removes indel variants", {
  td <- make_test_sumstats_ld(n_variants = 5, with_indels = TRUE)
  result <- rss_basic_qc(td$sumstats, td$LD_data, remove_indels = TRUE)
  expect_type(result, "list")
  expect_lte(nrow(result$sumstats), nrow(td$sumstats))
})

test_that("rss_basic_qc errors when no variants overlap", {
  set.seed(99)
  sumstats <- data.frame(
    chrom = c(1, 1),
    pos   = c(10000, 20000),
    A1    = c("G", "T"),
    A2    = c("A", "C"),
    beta  = c(0.5, -0.3),
    se    = c(0.1, 0.15),
    z     = c(5.0, -2.0),
    stringsAsFactors = FALSE
  )

  ld_ids <- c("1:50000:A:G", "1:60000:C:T")
  LD_mat <- diag(2)
  rownames(LD_mat) <- colnames(LD_mat) <- ld_ids

  ref_panel <- data.frame(
    chrom = c(1, 1),
    pos   = c(50000, 60000),
    A2    = c("A", "C"),
    A1    = c("G", "T"),
    stringsAsFactors = FALSE
  )

  LD_data <- list(
    LD_variants = ref_panel,
    LD_matrix   = LD_mat
  )

  expect_error(rss_basic_qc(sumstats, LD_data), "No overlapping|No matching")
})

test_that("rss_basic_qc aligns variant IDs by stripping build suffix", {
  set.seed(55)
  sumstats <- data.frame(
    chrom = c(1, 1, 1),
    pos   = c(100, 200, 300),
    A1    = c("G", "T", "A"),
    A2    = c("A", "C", "G"),
    beta  = c(0.5, -0.3, 0.1),
    se    = c(0.1, 0.15, 0.2),
    z     = c(5.0, -2.0, 0.5),
    stringsAsFactors = FALSE
  )

  ld_ids <- c("1:100:A:G_b38", "1:200:C:T_b38", "1:300:G:A_b38")
  LD_mat <- diag(3) + 0.01
  diag(LD_mat) <- 1
  rownames(LD_mat) <- colnames(LD_mat) <- ld_ids

  ref_panel <- data.frame(
    chrom = c(1, 1, 1),
    pos   = c(100, 200, 300),
    A2    = c("A", "C", "G"),
    A1    = c("G", "T", "A"),
    stringsAsFactors = FALSE
  )

  LD_data <- list(
    LD_variants = ref_panel,
    LD_matrix   = LD_mat
  )

  result <- rss_basic_qc(sumstats, LD_data)
  expect_type(result, "list")
  expect_true(nrow(result$sumstats) > 0)
})

test_that("rss_basic_qc handles chr prefix differences during alignment", {
  set.seed(77)
  sumstats <- data.frame(
    chrom = c(1, 1),
    pos   = c(100, 200),
    A1    = c("G", "T"),
    A2    = c("A", "C"),
    beta  = c(0.5, -0.3),
    se    = c(0.1, 0.15),
    z     = c(5.0, -2.0),
    stringsAsFactors = FALSE
  )

  ld_ids <- c("chr1:100:A:G", "chr1:200:C:T")
  LD_mat <- diag(2)
  rownames(LD_mat) <- colnames(LD_mat) <- ld_ids

  ref_panel <- data.frame(
    chrom = c(1, 1),
    pos   = c(100, 200),
    A2    = c("A", "C"),
    A1    = c("G", "T"),
    stringsAsFactors = FALSE
  )

  LD_data <- list(
    LD_variants = ref_panel,
    LD_matrix   = LD_mat
  )

  result <- rss_basic_qc(sumstats, LD_data)
  expect_type(result, "list")
  expect_true(nrow(result$sumstats) > 0)
})

test_that("rss_basic_qc output LD_mat has same dimension as sumstats rows", {
  td <- make_test_sumstats_ld(n_variants = 6)
  result <- rss_basic_qc(td$sumstats, td$LD_data)
  expect_equal(nrow(result$LD_mat), nrow(result$sumstats))
  expect_equal(ncol(result$LD_mat), nrow(result$sumstats))
})

test_that("rss_basic_qc errors when LD matrix has NULL rownames", {
  td <- make_test_sumstats_ld(n_variants = 3)
  ld_mat <- td$LD_data$LD_matrix
  rownames(ld_mat) <- NULL
  colnames(ld_mat) <- NULL
  td$LD_data$LD_matrix <- ld_mat

  expect_error(rss_basic_qc(td$sumstats, td$LD_data), "rownames are NULL|cannot align")
})

test_that("rss_basic_qc handles multiple skip regions", {
  td <- make_test_sumstats_ld(n_variants = 10)
  result <- rss_basic_qc(td$sumstats, td$LD_data,
                          skip_region = c("1:099-250", "1:650-850"))
  remaining_pos <- result$sumstats$pos
  expect_false(100 %in% remaining_pos)
  expect_false(200 %in% remaining_pos)
  expect_false(700 %in% remaining_pos)
  expect_false(800 %in% remaining_pos)
  expect_true(500 %in% remaining_pos)
})

# ===========================================================================
# summary_stats_qc
# ===========================================================================

test_that("summary_stats_qc errors on invalid method", {
  sumstats <- data.frame(variant_id = "1:100:A:G", z = 2.0)
  LD_data <- list(LD_matrix = matrix(1, 1, 1, dimnames = list("1:100:A:G", "1:100:A:G")))
  expect_error(summary_stats_qc(sumstats, LD_data, method = "invalid"),
               "Invalid quality control method")
})

test_that("summary_stats_qc with slalom method returns correct structure", {
  td <- make_test_sumstats_ld(n_variants = 5)
  basic_result <- rss_basic_qc(td$sumstats, td$LD_data)

  local_mocked_bindings(
    slalom = function(zScore, R, ...) {
      n <- length(zScore)
      list(
        data = data.frame(
          zScore   = zScore,
          outliers = c(rep(FALSE, n - 1), TRUE)
        )
      )
    }
  )

  result <- summary_stats_qc(
    basic_result$sumstats, td$LD_data,
    n = 10000, var_y = 1, method = "slalom"
  )

  expect_type(result, "list")
  expect_true("sumstats" %in% names(result))
  expect_true("LD_mat" %in% names(result))
  expect_true("outlier_number" %in% names(result))
  expect_equal(result$outlier_number, 1)
  expect_equal(nrow(result$sumstats), nrow(basic_result$sumstats) - 1)
})

test_that("summary_stats_qc with slalom and no outliers keeps all variants", {
  td <- make_test_sumstats_ld(n_variants = 4)
  basic_result <- rss_basic_qc(td$sumstats, td$LD_data)

  local_mocked_bindings(
    slalom = function(zScore, R, ...) {
      n <- length(zScore)
      list(
        data = data.frame(
          zScore   = zScore,
          outliers = rep(FALSE, n)
        )
      )
    }
  )

  result <- summary_stats_qc(
    basic_result$sumstats, td$LD_data,
    n = 10000, var_y = 1, method = "slalom"
  )
  expect_equal(result$outlier_number, 0)
  expect_equal(nrow(result$sumstats), nrow(basic_result$sumstats))
})

test_that("summary_stats_qc with dentist method returns correct structure", {
  td <- make_test_sumstats_ld(n_variants = 5)
  basic_result <- rss_basic_qc(td$sumstats, td$LD_data)

  local_mocked_bindings(
    dentist_single_window = function(zScore, R, nSample, ...) {
      n <- length(zScore)
      data.frame(
        z_score = zScore,
        outlier = c(TRUE, rep(FALSE, n - 1))
      )
    }
  )

  result <- summary_stats_qc(
    basic_result$sumstats, td$LD_data,
    n = 10000, var_y = 1, method = "dentist"
  )

  expect_type(result, "list")
  expect_true("sumstats" %in% names(result))
  expect_true("LD_mat" %in% names(result))
  expect_true("outlier_number" %in% names(result))
  expect_equal(result$outlier_number, 1)
  expect_equal(nrow(result$sumstats), nrow(basic_result$sumstats) - 1)
})

test_that("summary_stats_qc with dentist and all outliers returns empty", {
  td <- make_test_sumstats_ld(n_variants = 3)
  basic_result <- rss_basic_qc(td$sumstats, td$LD_data)

  local_mocked_bindings(
    dentist_single_window = function(zScore, R, nSample, ...) {
      n <- length(zScore)
      data.frame(
        z_score = zScore,
        outlier = rep(TRUE, n)
      )
    }
  )

  result <- summary_stats_qc(
    basic_result$sumstats, td$LD_data,
    n = 10000, var_y = 1, method = "dentist"
  )
  expect_equal(nrow(result$sumstats), 0)
  expect_equal(result$outlier_number, nrow(basic_result$sumstats))
})

test_that("summary_stats_qc returns LD_mat matching filtered sumstats dimensions", {
  td <- make_test_sumstats_ld(n_variants = 6)
  basic_result <- rss_basic_qc(td$sumstats, td$LD_data)

  local_mocked_bindings(
    slalom = function(zScore, R, ...) {
      n <- length(zScore)
      outlier_flags <- rep(FALSE, n)
      outlier_flags[c(1, 3)] <- TRUE
      list(data = data.frame(zScore = zScore, outliers = outlier_flags))
    }
  )

  result <- summary_stats_qc(
    basic_result$sumstats, td$LD_data,
    n = 10000, var_y = 1, method = "slalom"
  )
  expect_equal(nrow(result$LD_mat), nrow(result$sumstats))
  expect_equal(ncol(result$LD_mat), nrow(result$sumstats))
})
