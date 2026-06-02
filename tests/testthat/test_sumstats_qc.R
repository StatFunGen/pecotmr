context("sumstats_qc")

# ===========================================================================
# Helper: build an LDData S4 object from a correlation matrix and variant info
# ===========================================================================
make_ld_data_s4 <- function(R_mat, variant_ids, chrom_val = 1, positions = NULL) {
  ref_panel <- parse_variant_id(variant_ids)
  ref_panel$variant_id <- variant_ids
  ref_panel$chrom <- as.character(ref_panel$chrom)
  if (!is.null(positions)) {
    ref_panel$pos <- positions
  }
  variants_gr <- pecotmr:::.ref_panel_to_granges(ref_panel)
  bm <- data.frame(
    block_id = 1L,
    chrom = ref_panel$chrom[1],
    block_start = min(ref_panel$pos),
    block_end = max(ref_panel$pos),
    size = length(variant_ids),
    start_idx = 1L,
    end_idx = length(variant_ids),
    stringsAsFactors = FALSE
  )
  LDData(correlation = R_mat, variants = variants_gr, block_metadata = bm)
}

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

  LD_data <- make_ld_data_s4(LD_mat, variant_ids, chrom_val = chrom_val,
                             positions = positions)

  list(sumstats = sumstats, LD_data = LD_data, variant_ids = variant_ids)
}

# ===========================================================================
# rss_basic_qc
# ===========================================================================

test_that("rss_basic_qc requires correct columns", {
  sumstats <- data.frame(beta = 1, se = 0.5)
  R_mat <- matrix(1, 1, 1, dimnames = list("1:100:A:G", "1:100:A:G"))
  LD_data <- make_ld_data_s4(R_mat, "1:100:A:G")
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

  LD_data <- make_ld_data_s4(LD_mat, variant_ids)

  result <- rss_basic_qc(sumstats, LD_data)
  expect_true(is(result, "QCResult"))
  expect_true(!is.null(getRSSInput(result)$sumstats))
  expect_true(!is.null(getLDData(result)))
})

test_that("rss_basic_qc skips variants in specified region", {
  td <- make_test_sumstats_ld(n_variants = 5)

  result <- rss_basic_qc(td$sumstats, td$LD_data, skip_region = "1:150-350")

  expect_true(is(result, "QCResult"))
  expect_true(!is.null(getRSSInput(result)$sumstats))
  expect_true(!is.null(getLDData(result)))
  remaining_pos <- getRSSInput(result)$sumstats$pos
  expect_false(200 %in% remaining_pos)
  expect_false(300 %in% remaining_pos)
})

test_that("rss_basic_qc with skip_region preserves non-skipped variants", {
  td <- make_test_sumstats_ld(n_variants = 5)
  result <- rss_basic_qc(td$sumstats, td$LD_data, skip_region = "1:150-250")
  remaining_pos <- getRSSInput(result)$sumstats$pos
  expect_false(200 %in% remaining_pos)
  expect_true(100 %in% remaining_pos)
  expect_true(300 %in% remaining_pos)
})

test_that("rss_basic_qc with keep_indel=FALSE removes indel variants", {
  td <- make_test_sumstats_ld(n_variants = 5, with_indels = TRUE)
  result <- rss_basic_qc(td$sumstats, td$LD_data, keep_indel = FALSE)
  expect_true(is(result, "QCResult"))
  expect_lte(nrow(getRSSInput(result)$sumstats), nrow(td$sumstats))
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

  LD_data <- make_ld_data_s4(LD_mat, ld_ids)

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

  # For the build-suffix variant IDs, construct the LDData with the
  # base IDs (without suffix) for variant metadata so parse_variant_id works,
  # while the correlation matrix retains the suffixed rownames.
  base_ids <- c("1:100:A:G", "1:200:C:T", "1:300:G:A")
  LD_data <- make_ld_data_s4(LD_mat, base_ids)

  result <- rss_basic_qc(sumstats, LD_data)
  expect_true(is(result, "QCResult"))
  expect_true(nrow(getRSSInput(result)$sumstats) > 0)
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

  # Use base IDs (without chr prefix) for variant metadata, while
  # the correlation matrix has chr-prefixed rownames.
  base_ids <- c("1:100:A:G", "1:200:C:T")
  LD_data <- make_ld_data_s4(LD_mat, base_ids)

  result <- rss_basic_qc(sumstats, LD_data)
  expect_true(is(result, "QCResult"))
  expect_true(nrow(getRSSInput(result)$sumstats) > 0)
})

test_that("rss_basic_qc output LD_mat has same dimension as sumstats rows", {
  td <- make_test_sumstats_ld(n_variants = 6)
  result <- rss_basic_qc(td$sumstats, td$LD_data)
  result_ld_mat <- getCorrelation(getLDData(result))
  result_sumstats <- getRSSInput(result)$sumstats
  expect_equal(nrow(result_ld_mat), nrow(result_sumstats))
  expect_equal(ncol(result_ld_mat), nrow(result_sumstats))
})

test_that("rss_basic_qc errors when LD matrix has NULL rownames", {
  td <- make_test_sumstats_ld(n_variants = 3)
  ld_mat <- getCorrelation(td$LD_data)
  rownames(ld_mat) <- NULL
  colnames(ld_mat) <- NULL
  # Rebuild the LDData with a NULL-rownames correlation matrix
  LD_data_bad <- LDData(
    correlation = ld_mat,
    variants = getVariantInfo(td$LD_data),
    block_metadata = getBlockMetadata(td$LD_data)
  )

  expect_error(rss_basic_qc(td$sumstats, LD_data_bad), "rownames are NULL|cannot align")
})

test_that("rss_basic_qc handles multiple skip regions", {
  td <- make_test_sumstats_ld(n_variants = 10)
  result <- rss_basic_qc(td$sumstats, td$LD_data,
                          skip_region = c("1:099-250", "1:650-850"))
  remaining_pos <- getRSSInput(result)$sumstats$pos
  expect_false(100 %in% remaining_pos)
  expect_false(200 %in% remaining_pos)
  expect_false(700 %in% remaining_pos)
  expect_false(800 %in% remaining_pos)
  expect_true(500 %in% remaining_pos)
})

test_that("rss_basic_qc can skip LD matrix subsetting for genotype references", {
  td <- make_test_sumstats_ld(n_variants = 3)
  X_ref <- matrix(rnorm(30), 10, 3)
  colnames(X_ref) <- td$variant_ids
  # Store X_ref as correlation; with return_LD_mat=FALSE the matrix is not subsetted
  LD_data_geno <- LDData(
    correlation = X_ref,
    variants = getVariantInfo(td$LD_data),
    block_metadata = getBlockMetadata(td$LD_data)
  )

  result <- rss_basic_qc(td$sumstats, LD_data_geno, return_LD_mat = FALSE)

  expect_true(nrow(getRSSInput(result)$sumstats) > 0)
  expect_null(getLDData(result))
})

# ===========================================================================
# summary_stats_qc
# ===========================================================================

test_that("summary_stats_qc errors on invalid method", {
  sumstats <- data.frame(variant_id = "1:100:A:G", z = 2.0)
  R_mat <- matrix(1, 1, 1, dimnames = list("1:100:A:G", "1:100:A:G"))
  LD_data <- make_ld_data_s4(R_mat, "1:100:A:G")
  expect_error(summary_stats_qc(sumstats, LD_data, method = "invalid"),
               "should be one of")
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
    getRSSInput(basic_result)$sumstats, td$LD_data,
    n = 10000, method = "slalom"
  )

  expect_true(is(result, "QCResult"))
  expect_true(!is.null(getRSSInput(result)$sumstats))
  expect_true(!is.null(getLDData(result)))
  expect_equal(getOutlierNumber(result), 1)
  expect_equal(nrow(getRSSInput(result)$sumstats),
               nrow(getRSSInput(basic_result)$sumstats) - 1)
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
    getRSSInput(basic_result)$sumstats, td$LD_data,
    n = 10000, method = "slalom"
  )
  expect_equal(getOutlierNumber(result), 0)
  expect_equal(nrow(getRSSInput(result)$sumstats),
               nrow(getRSSInput(basic_result)$sumstats))
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
    getRSSInput(basic_result)$sumstats, td$LD_data,
    n = 10000, method = "dentist"
  )

  expect_true(is(result, "QCResult"))
  expect_true(!is.null(getRSSInput(result)$sumstats))
  expect_true(!is.null(getLDData(result)))
  expect_equal(getOutlierNumber(result), 1)
  expect_equal(nrow(getRSSInput(result)$sumstats),
               nrow(getRSSInput(basic_result)$sumstats) - 1)
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
    getRSSInput(basic_result)$sumstats, td$LD_data,
    n = 10000, method = "dentist"
  )
  expect_equal(nrow(getRSSInput(result)$sumstats), 0)
  expect_equal(getOutlierNumber(result),
               nrow(getRSSInput(basic_result)$sumstats))
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
    getRSSInput(basic_result)$sumstats, td$LD_data,
    n = 10000, method = "slalom"
  )
  result_ld_mat <- getCorrelation(getLDData(result))
  result_sumstats <- getRSSInput(result)$sumstats
  expect_equal(nrow(result_ld_mat), nrow(result_sumstats))
  expect_equal(ncol(result_ld_mat), nrow(result_sumstats))
})

test_that("summary_stats_qc basic genotype-backed path does not compute LD", {
  td <- make_test_sumstats_ld(n_variants = 5)
  X_ref <- matrix(rnorm(50), 10, 5)
  colnames(X_ref) <- td$variant_ids
  LD_data_geno <- make_ld_data_s4(cor(X_ref), td$variant_ids)
  rss_input <- list(sumstats = td$sumstats, n = 1000, var_y = 1)

  local_mocked_bindings(
    compute_LD = function(...) stop("compute_LD should not be called"),
    hasGenotypes = function(x) TRUE,
    getGenotypes = function(x) X_ref,
    .package = "pecotmr"
  )

  expect_message(
    result <- summary_stats_qc(rss_input = rss_input, LD_data = LD_data_geno,
                               qc_method = "none", impute = FALSE),
    "basic harmonization retained"
  )
  result_ld <- getLDData(result)
  result_geno <- getGenotypes(result_ld)
  expect_equal(nrow(result_geno), nrow(X_ref))
  expect_equal(ncol(result_geno), nrow(getRSSInput(result)$sumstats))
})

test_that("summary_stats_qc accepts genotype-backed LDData", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  tmp <- tempfile("lddata_qc_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  prefix <- "test_variants"
  for (ext in c("pgen", "pvar", "psam", "afreq")) {
    file.copy(file.path(td, paste0(prefix, ".", ext)),
              file.path(tmp, paste0(prefix, ".", ext)))
  }
  meta_file <- file.path(tmp, "ld_meta.tsv")
  writeLines(c("chrom\tstart\tend\tpath", "21\t0\t0\ttest_variants"), meta_file)

  ld_data <- suppressWarnings(suppressMessages(load_LD_matrix(
    meta_file,
    region = "chr21:17513228-17550000",
    return_genotype = TRUE
  )))
  variant_info <- getVariantInfo(ld_data)
  ref_panel <- as.data.frame(S4Vectors::mcols(variant_info))
  ref_panel$chrom <- as.character(GenomicRanges::seqnames(variant_info))
  ref_panel$pos <- GenomicRanges::start(variant_info)
  is_snp <- nchar(ref_panel$A1) == 1 & nchar(ref_panel$A2) == 1
  allele_pair <- apply(cbind(ref_panel$A1, ref_panel$A2), 1, function(x) {
    paste(sort(x), collapse = "")
  })
  ref_panel <- ref_panel[is_snp & !allele_pair %in% c("AT", "CG"), , drop = FALSE]
  ref_panel <- utils::head(ref_panel, 5)

  sumstats <- data.frame(
    chrom = ref_panel$chrom,
    pos = ref_panel$pos,
    A1 = ref_panel$A1,
    A2 = ref_panel$A2,
    beta = seq_len(nrow(ref_panel)) / 10,
    se = rep(0.1, nrow(ref_panel)),
    z = seq_len(nrow(ref_panel)),
    stringsAsFactors = FALSE
  )
  rss_input <- list(sumstats = sumstats, n = 1000, var_y = 1)

  local_mocked_bindings(
    compute_LD = function(...) stop("compute_LD should not be called")
  )
  result <- suppressMessages(summary_stats_qc(
    rss_input = rss_input,
    LD_data = ld_data,
    qc_method = "none",
    impute = FALSE
  ))

  result_ld <- getLDData(result)
  result_geno <- getGenotypes(result_ld)
  expect_equal(nrow(result_geno), 100L)
  expect_equal(ncol(result_geno), nrow(getRSSInput(result)$sumstats))
})

test_that("summary_stats_qc PIP screening uses LD-independent SER", {
  td <- make_test_sumstats_ld(n_variants = 5)
  X_ref <- matrix(rnorm(50), 10, 5)
  colnames(X_ref) <- td$variant_ids
  LD_data_geno <- make_ld_data_s4(cor(X_ref), td$variant_ids)
  rss_input <- list(sumstats = td$sumstats, n = 1000, var_y = 1)

  local_mocked_bindings(
    compute_LD = function(...) stop("compute_LD should not be called"),
    susie_ser = function(z, n = NULL, coverage = 0.95, ...) {
      expect_equal(n, rss_input$n)
      expect_null(coverage)
      list(pip = rep(1, length(z)))
    },
    .package = "pecotmr"
  )

  result <- suppressMessages(summary_stats_qc(
    rss_input = rss_input,
    LD_data = LD_data_geno,
    qc_method = "none",
    pip_cutoff_to_skip = 0.1,
    impute = FALSE
  ))
  result_ld <- getLDData(result)
  result_R <- getCorrelation(result_ld)
  expect_equal(ncol(result_R), nrow(getRSSInput(result)$sumstats))
})

test_that("summary_stats_qc treats NULL qc_method as basic-only none", {
  td <- make_test_sumstats_ld(n_variants = 5)
  rss_input <- list(sumstats = td$sumstats, n = 1000, var_y = 1)

  local_mocked_bindings(
    ld_mismatch_qc = function(...) stop("ld_mismatch_qc should not be called")
  )

  expect_message(
    result <- summary_stats_qc(
      rss_input = rss_input,
      LD_data = td$LD_data,
      qc_method = NULL,
      impute = FALSE
    ),
    "basic harmonization retained"
  )
  expect_equal(nrow(getRSSInput(result)$sumstats), nrow(td$sumstats))
})

test_that("summary_stats_qc rejects invalid qc_method values", {
  td <- make_test_sumstats_ld(n_variants = 5)
  rss_input <- list(sumstats = td$sumstats, n = 1000, var_y = 1)

  expect_error(
    summary_stats_qc(
      rss_input = rss_input,
      LD_data = td$LD_data,
      qc_method = "bad_method"
    ),
    "should be one of"
  )
})

test_that("summary_stats_qc LD-mismatch QC computes only filtered local LD from X_ref", {
  td <- make_test_sumstats_ld(n_variants = 5)
  X_ref <- matrix(rnorm(50), 10, 5)
  colnames(X_ref) <- td$variant_ids
  LD_data_geno <- make_ld_data_s4(cor(X_ref), td$variant_ids)
  rss_input <- list(sumstats = td$sumstats, n = 1000, var_y = 1)
  compute_calls <- 0

  local_mocked_bindings(
    compute_LD = function(X, ...) {
      compute_calls <<- compute_calls + 1
      expect_equal(ncol(X), 3)
      R <- diag(ncol(X))
      rownames(R) <- colnames(R) <- colnames(X)
      R
    },
    hasGenotypes = function(x) TRUE,
    getGenotypes = function(x) X_ref,
    ld_mismatch_qc = function(zScore, R, nSample = NULL, method = NULL, ...) {
      expect_equal(nrow(R), length(zScore))
      expect_equal(ncol(R), length(zScore))
      data.frame(outlier = rep(FALSE, length(zScore)))
    },
    .package = "pecotmr"
  )

  result <- suppressMessages(summary_stats_qc(
    rss_input = rss_input,
    LD_data = LD_data_geno,
    qc_method = "slalom",
    skip_region = "1:150-350",
    impute = FALSE
  ))
  expect_equal(compute_calls, 2)
  result_ld <- getLDData(result)
  # getGenotypes is mocked above to always return full X_ref, so read the
  # subsetted handle stored in the LDData slot directly to verify subsetting.
  result_geno <- result_ld@genotype_handle
  expect_equal(ncol(result_geno), nrow(getRSSInput(result)$sumstats))
  expect_equal(ncol(result_geno), 3)
})

# ===========================================================================
# ld_mismatch_qc
# ===========================================================================

test_that("ld_mismatch_qc with dentist method returns data frame with outlier column", {
  set.seed(42)
  p <- 20
  R <- diag(p)
  z <- rnorm(p)
  result <- ld_mismatch_qc(z, R = R, nSample = 1000, method = "dentist")
  expect_true(is.data.frame(result) || is.list(result))
  expect_true("outlier" %in% names(result))
})

test_that("ld_mismatch_qc with slalom method returns data frame with outlier column", {
  set.seed(42)
  p <- 20
  R <- diag(p)
  z <- rnorm(p)
  result <- ld_mismatch_qc(z, R = R, method = "slalom")
  expect_true(is.data.frame(result) || is.list(result))
  expect_true("outlier" %in% names(result))
})

test_that("ld_mismatch_qc method argument is validated", {
  z <- rnorm(5)
  R <- diag(5)
  expect_error(ld_mismatch_qc(z, R = R, method = "invalid"))
})
