context("sumstats_qc")

# ===========================================================================
# Helper: build an LdData S4 object from a correlation matrix and variant info
# ===========================================================================
make_ld_data_s4 <- function(R_mat, variant_ids, chrom_val = 1, positions = NULL) {
  ref_panel <- parseVariantId(variant_ids)
  ref_panel$variant_id <- variant_ids
  ref_panel$chrom <- as.character(ref_panel$chrom)
  if (!is.null(positions)) {
    ref_panel$pos <- positions
  }
  variants_gr <- pecotmr:::.refPanelToGranges(ref_panel)
  bm <- data.frame(
    blockId = 1L,
    chrom = ref_panel$chrom[1],
    blockStart = min(ref_panel$pos),
    blockEnd = max(ref_panel$pos),
    size = length(variant_ids),
    startIdx = 1L,
    endIdx = length(variant_ids),
    stringsAsFactors = FALSE
  )
  LdData(correlation = R_mat, variants = variants_gr, blockMetadata = bm)
}

# ===========================================================================
# Helper: build matching sumstats and ldData
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

  ldMat <- diag(n_variants) + matrix(0.01, n_variants, n_variants)
  diag(ldMat) <- 1
  rownames(ldMat) <- colnames(ldMat) <- variant_ids

  ldData <- make_ld_data_s4(ldMat, variant_ids, chrom_val = chrom_val,
                             positions = positions)

  list(sumstats = sumstats, ldData = ldData,
       variantIds = variant_ids, variant_ids = variant_ids)
}

# ===========================================================================
# rssBasicQc
# ===========================================================================

test_that("rssBasicQc requires correct columns", {
  sumstats <- data.frame(beta = 1, se = 0.5)
  R_mat <- matrix(1, 1, 1, dimnames = list("1:100:A:G", "1:100:A:G"))
  ldData <- make_ld_data_s4(R_mat, "1:100:A:G")
  expect_error(rssBasicQc(sumstats, ldData), "Missing columns")
})

test_that("rssBasicQc processes matching variants correctly", {
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

  ldMat <- diag(3)
  rownames(ldMat) <- colnames(ldMat) <- variant_ids

  ldData <- make_ld_data_s4(ldMat, variant_ids)

  result <- rssBasicQc(sumstats, ldData)
  expect_true(is(result, "QcResult"))
  expect_true(!is.null(getRssInput(result)$sumstats))
  expect_true(!is.null(getLdData(result)))
})

test_that("rssBasicQc skips variants in specified region", {
  td <- make_test_sumstats_ld(n_variants = 5)

  result <- rssBasicQc(td$sumstats, td$ldData, skipRegion = "1:150-350")

  expect_true(is(result, "QcResult"))
  expect_true(!is.null(getRssInput(result)$sumstats))
  expect_true(!is.null(getLdData(result)))
  remaining_pos <- getRssInput(result)$sumstats$pos
  expect_false(200 %in% remaining_pos)
  expect_false(300 %in% remaining_pos)
})

test_that("rssBasicQc with skip_region preserves non-skipped variants", {
  td <- make_test_sumstats_ld(n_variants = 5)
  result <- rssBasicQc(td$sumstats, td$ldData, skipRegion = "1:150-250")
  remaining_pos <- getRssInput(result)$sumstats$pos
  expect_false(200 %in% remaining_pos)
  expect_true(100 %in% remaining_pos)
  expect_true(300 %in% remaining_pos)
})

test_that("rssBasicQc with keep_indel=FALSE removes indel variants", {
  td <- make_test_sumstats_ld(n_variants = 5, with_indels = TRUE)
  result <- rssBasicQc(td$sumstats, td$ldData, keepIndel = FALSE)
  expect_true(is(result, "QcResult"))
  expect_lte(nrow(getRssInput(result)$sumstats), nrow(td$sumstats))
})

test_that("rssBasicQc errors when no variants overlap", {
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
  ldMat <- diag(2)
  rownames(ldMat) <- colnames(ldMat) <- ld_ids

  ldData <- make_ld_data_s4(ldMat, ld_ids)

  expect_error(rssBasicQc(sumstats, ldData), "No overlapping|No matching")
})

test_that("rssBasicQc aligns variant IDs by stripping build suffix", {
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
  ldMat <- diag(3) + 0.01
  diag(ldMat) <- 1
  rownames(ldMat) <- colnames(ldMat) <- ld_ids

  # For the build-suffix variant IDs, construct the LdData with the
  # base IDs (without suffix) for variant metadata so parseVariantId works,
  # while the correlation matrix retains the suffixed rownames.
  base_ids <- c("1:100:A:G", "1:200:C:T", "1:300:G:A")
  ldData <- make_ld_data_s4(ldMat, base_ids)

  result <- rssBasicQc(sumstats, ldData)
  expect_true(is(result, "QcResult"))
  expect_true(nrow(getRssInput(result)$sumstats) > 0)
})

test_that("rssBasicQc handles chr prefix differences during alignment", {
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
  ldMat <- diag(2)
  rownames(ldMat) <- colnames(ldMat) <- ld_ids

  # Use base IDs (without chr prefix) for variant metadata, while
  # the correlation matrix has chr-prefixed rownames.
  base_ids <- c("1:100:A:G", "1:200:C:T")
  ldData <- make_ld_data_s4(ldMat, base_ids)

  result <- rssBasicQc(sumstats, ldData)
  expect_true(is(result, "QcResult"))
  expect_true(nrow(getRssInput(result)$sumstats) > 0)
})

test_that("rssBasicQc output ldMat has same dimension as sumstats rows", {
  td <- make_test_sumstats_ld(n_variants = 6)
  result <- rssBasicQc(td$sumstats, td$ldData)
  result_ld_mat <- getCorrelation(getLdData(result))
  result_sumstats <- getRssInput(result)$sumstats
  expect_equal(nrow(result_ld_mat), nrow(result_sumstats))
  expect_equal(ncol(result_ld_mat), nrow(result_sumstats))
})

test_that("rssBasicQc errors when LD matrix has NULL rownames", {
  td <- make_test_sumstats_ld(n_variants = 3)
  ld_mat <- getCorrelation(td$ldData)
  rownames(ld_mat) <- NULL
  colnames(ld_mat) <- NULL
  # Rebuild the LdData with a NULL-rownames correlation matrix
  LD_data_bad <- LdData(
    correlation = ld_mat,
    variants = getVariantInfo(td$ldData),
    blockMetadata = getBlockMetadata(td$ldData)
  )

  expect_error(rssBasicQc(td$sumstats, LD_data_bad), "rownames are NULL|cannot align")
})

test_that("rssBasicQc handles multiple skip regions", {
  td <- make_test_sumstats_ld(n_variants = 10)
  result <- rssBasicQc(td$sumstats, td$ldData,
                          skipRegion = c("1:099-250", "1:650-850"))
  remaining_pos <- getRssInput(result)$sumstats$pos
  expect_false(100 %in% remaining_pos)
  expect_false(200 %in% remaining_pos)
  expect_false(700 %in% remaining_pos)
  expect_false(800 %in% remaining_pos)
  expect_true(500 %in% remaining_pos)
})

test_that("rssBasicQc can skip LD matrix subsetting for genotype references", {
  td <- make_test_sumstats_ld(n_variants = 3)
  X_ref <- matrix(rnorm(30), 10, 3)
  colnames(X_ref) <- td$variant_ids
  # Store X_ref as correlation; with return_LD_mat=FALSE the matrix is not subsetted
  LD_data_geno <- LdData(
    correlation = X_ref,
    variants = getVariantInfo(td$ldData),
    blockMetadata = getBlockMetadata(td$ldData)
  )

  result <- rssBasicQc(td$sumstats, LD_data_geno, returnLdMat = FALSE)

  expect_true(nrow(getRssInput(result)$sumstats) > 0)
  expect_null(getLdData(result))
})

# ===========================================================================
# summaryStatsQc
# ===========================================================================

test_that("summaryStatsQc errors on invalid method", {
  sumstats <- data.frame(variant_id = "1:100:A:G", z = 2.0)
  R_mat <- matrix(1, 1, 1, dimnames = list("1:100:A:G", "1:100:A:G"))
  ldData <- make_ld_data_s4(R_mat, "1:100:A:G")
  expect_error(summaryStatsQc(sumstats, ldData, method = "invalid"),
               "should be one of")
})

test_that("summaryStatsQc with slalom method returns correct structure", {
  td <- make_test_sumstats_ld(n_variants = 5)
  basic_result <- rssBasicQc(td$sumstats, td$ldData)

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

  result <- summaryStatsQc(
    getRssInput(basic_result)$sumstats, td$ldData,
    n = 10000, method = "slalom"
  )

  expect_true(is(result, "QcResult"))
  expect_true(!is.null(getRssInput(result)$sumstats))
  expect_true(!is.null(getLdData(result)))
  expect_equal(getOutlierNumber(result), 1)
  expect_equal(nrow(getRssInput(result)$sumstats),
               nrow(getRssInput(basic_result)$sumstats) - 1)
})

test_that("summaryStatsQc with slalom and no outliers keeps all variants", {
  td <- make_test_sumstats_ld(n_variants = 4)
  basic_result <- rssBasicQc(td$sumstats, td$ldData)

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

  result <- summaryStatsQc(
    getRssInput(basic_result)$sumstats, td$ldData,
    n = 10000, method = "slalom"
  )
  expect_equal(getOutlierNumber(result), 0)
  expect_equal(nrow(getRssInput(result)$sumstats),
               nrow(getRssInput(basic_result)$sumstats))
})

test_that("summaryStatsQc with dentist method returns correct structure", {
  td <- make_test_sumstats_ld(n_variants = 5)
  basic_result <- rssBasicQc(td$sumstats, td$ldData)

  local_mocked_bindings(
    dentistSingleWindow = function(zScore, R, nSample, ...) {
      n <- length(zScore)
      data.frame(
        z_score = zScore,
        outlier = c(TRUE, rep(FALSE, n - 1))
      )
    }
  )

  result <- summaryStatsQc(
    getRssInput(basic_result)$sumstats, td$ldData,
    n = 10000, method = "dentist"
  )

  expect_true(is(result, "QcResult"))
  expect_true(!is.null(getRssInput(result)$sumstats))
  expect_true(!is.null(getLdData(result)))
  expect_equal(getOutlierNumber(result), 1)
  expect_equal(nrow(getRssInput(result)$sumstats),
               nrow(getRssInput(basic_result)$sumstats) - 1)
})

test_that("summaryStatsQc with dentist and all outliers returns empty", {
  td <- make_test_sumstats_ld(n_variants = 3)
  basic_result <- rssBasicQc(td$sumstats, td$ldData)

  local_mocked_bindings(
    dentistSingleWindow = function(zScore, R, nSample, ...) {
      n <- length(zScore)
      data.frame(
        z_score = zScore,
        outlier = rep(TRUE, n)
      )
    }
  )

  result <- summaryStatsQc(
    getRssInput(basic_result)$sumstats, td$ldData,
    n = 10000, method = "dentist"
  )
  expect_equal(nrow(getRssInput(result)$sumstats), 0)
  expect_equal(getOutlierNumber(result),
               nrow(getRssInput(basic_result)$sumstats))
})

test_that("summaryStatsQc returns ldMat matching filtered sumstats dimensions", {
  td <- make_test_sumstats_ld(n_variants = 6)
  basic_result <- rssBasicQc(td$sumstats, td$ldData)

  local_mocked_bindings(
    slalom = function(zScore, R, ...) {
      n <- length(zScore)
      outlier_flags <- rep(FALSE, n)
      outlier_flags[c(1, 3)] <- TRUE
      list(data = data.frame(zScore = zScore, outliers = outlier_flags))
    }
  )

  result <- summaryStatsQc(
    getRssInput(basic_result)$sumstats, td$ldData,
    n = 10000, method = "slalom"
  )
  result_ld_mat <- getCorrelation(getLdData(result))
  result_sumstats <- getRssInput(result)$sumstats
  expect_equal(nrow(result_ld_mat), nrow(result_sumstats))
  expect_equal(ncol(result_ld_mat), nrow(result_sumstats))
})

test_that("summaryStatsQc basic genotype-backed path does not compute LD", {
  td <- make_test_sumstats_ld(n_variants = 5)
  X_ref <- matrix(rnorm(50), 10, 5)
  colnames(X_ref) <- td$variant_ids
  LD_data_geno <- make_ld_data_s4(cor(X_ref), td$variant_ids)
  rss_input <- list(sumstats = td$sumstats, n = 1000, varY = 1)

  local_mocked_bindings(
    computeLd = function(...) stop("computeLd should not be called"),
    hasGenotypes = function(x) TRUE,
    getGenotypes = function(x) X_ref,
    .package = "pecotmr"
  )

  expect_message(
    result <- summaryStatsQc(rssInput = rss_input, ldData = LD_data_geno,
                               zMismatchQc = "none", impute = FALSE),
    "basic harmonization retained"
  )
  result_ld <- getLdData(result)
  result_geno <- getGenotypes(result_ld)
  expect_equal(nrow(result_geno), nrow(X_ref))
  expect_equal(ncol(result_geno), nrow(getRssInput(result)$sumstats))
})

test_that("summaryStatsQc accepts genotype-backed LdData", {
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

  ld_data <- suppressWarnings(suppressMessages(loadLdMatrix(
    meta_file,
    region = "chr21:17513228-17550000",
    returnGenotype = TRUE
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
  rss_input <- list(sumstats = sumstats, n = 1000, varY = 1)

  local_mocked_bindings(
    computeLd = function(...) stop("computeLd should not be called")
  )
  result <- suppressMessages(summaryStatsQc(
    rssInput = rss_input,
    ldData = ld_data,
    zMismatchQc = "none",
    impute = FALSE
  ))

  result_ld <- getLdData(result)
  result_geno <- getGenotypes(result_ld)
  expect_equal(nrow(result_geno), 100L)
  expect_equal(ncol(result_geno), nrow(getRssInput(result)$sumstats))
})

test_that("summaryStatsQc PIP screening uses LD-independent SER", {
  td <- make_test_sumstats_ld(n_variants = 5)
  X_ref <- matrix(rnorm(50), 10, 5)
  colnames(X_ref) <- td$variant_ids
  LD_data_geno <- make_ld_data_s4(cor(X_ref), td$variant_ids)
  rss_input <- list(sumstats = td$sumstats, n = 1000, varY = 1)

  local_mocked_bindings(
    computeLd = function(...) stop("computeLd should not be called"),
    susie_ser = function(z, n = NULL, coverage = 0.95, ...) {
      expect_equal(n, rss_input$n)
      expect_null(coverage)
      list(pip = rep(1, length(z)))
    },
    .package = "pecotmr"
  )

  result <- suppressMessages(summaryStatsQc(
    rssInput = rss_input,
    ldData = LD_data_geno,
    zMismatchQc = "none",
    pipCutoffToSkip = 0.1,
    impute = FALSE
  ))
  result_ld <- getLdData(result)
  result_R <- getCorrelation(result_ld)
  expect_equal(ncol(result_R), nrow(getRssInput(result)$sumstats))
})

test_that("summaryStatsQc treats NULL qc_method as basic-only none", {
  td <- make_test_sumstats_ld(n_variants = 5)
  rss_input <- list(sumstats = td$sumstats, n = 1000, varY = 1)

  local_mocked_bindings(
    ldMismatchQc = function(...) stop("ldMismatchQc should not be called")
  )

  expect_message(
    result <- summaryStatsQc(
      rssInput = rss_input,
      ldData = td$ldData,
      zMismatchQc = NULL,
      impute = FALSE
    ),
    "basic harmonization retained"
  )
  expect_equal(nrow(getRssInput(result)$sumstats), nrow(td$sumstats))
})

test_that("summaryStatsQc rejects invalid qc_method values", {
  td <- make_test_sumstats_ld(n_variants = 5)
  rss_input <- list(sumstats = td$sumstats, n = 1000, varY = 1)

  expect_error(
    summaryStatsQc(
      rssInput = rss_input,
      ldData = td$ldData,
      zMismatchQc = "bad_method"
    ),
    "should be one of"
  )
})

test_that("summaryStatsQc LD-mismatch QC computes only filtered local LD from X_ref", {
  td <- make_test_sumstats_ld(n_variants = 5)
  X_ref <- matrix(rnorm(50), 10, 5)
  colnames(X_ref) <- td$variant_ids
  LD_data_geno <- make_ld_data_s4(cor(X_ref), td$variant_ids)
  rss_input <- list(sumstats = td$sumstats, n = 1000, varY = 1)
  compute_calls <- 0

  local_mocked_bindings(
    computeLd = function(X, ...) {
      compute_calls <<- compute_calls + 1
      expect_equal(ncol(X), 3)
      R <- diag(ncol(X))
      rownames(R) <- colnames(R) <- colnames(X)
      R
    },
    hasGenotypes = function(x) TRUE,
    getGenotypes = function(x) X_ref,
    ldMismatchQc = function(zScore, R, nSample = NULL, method = NULL, ...) {
      expect_equal(nrow(R), length(zScore))
      expect_equal(ncol(R), length(zScore))
      data.frame(outlier = rep(FALSE, length(zScore)))
    },
    .package = "pecotmr"
  )

  result <- suppressMessages(summaryStatsQc(
    rssInput = rss_input,
    ldData = LD_data_geno,
    zMismatchQc = "slalom",
    skipRegion = "1:150-350",
    impute = FALSE
  ))
  expect_equal(compute_calls, 2)
  result_ld <- getLdData(result)
  # getGenotypes is mocked above to always return full X_ref, so read the
  # subsetted handle stored in the LdData slot directly to verify subsetting.
  result_geno <- result_ld@genotypeHandle
  expect_equal(ncol(result_geno), nrow(getRssInput(result)$sumstats))
  expect_equal(ncol(result_geno), 3)
})

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
# zMismatchQc resolver + hard rename (rss-qc-parity)
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

test_that("removed qcMethod argument is rejected as unknown (no alias)", {
  td <- make_test_sumstats_ld(n_variants = 5)
  rssInput <- list(sumstats = td$sumstats, n = 1000, varY = 1)
  expect_error(
    summaryStatsQc(rssInput = rssInput, ldData = td$ldData, qcMethod = "slalom"),
    "unused argument"
  )
})

# ===========================================================================
# krigingOutlierQc + alleleFlipKriging (rss-qc-parity)
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

test_that("kriging is not run by default (alleleFlipKriging = FALSE)", {
  td <- make_test_sumstats_ld(n_variants = 5)
  rssInput <- list(sumstats = td$sumstats, n = 1000, varY = 1)
  local_mocked_bindings(
    krigingOutlierQc = function(...) stop("krigingOutlierQc should not be called")
  )
  expect_message(
    result <- summaryStatsQc(
      rssInput = rssInput, ldData = td$ldData,
      zMismatchQc = "none", impute = FALSE
    ),
    "basic harmonization retained"
  )
  expect_equal(nrow(getRssInput(result)$sumstats), nrow(td$sumstats))
})

test_that("alleleFlipKriging drops a planted outlier and composes with zMismatchQc", {
  td <- make_test_sumstats_ld(n_variants = 6)
  rssInput <- list(sumstats = td$sumstats, n = 1000, varY = 1)
  seenN <- NULL
  local_mocked_bindings(
    krigingOutlierQc = function(zScore, R, variantIds = NULL, ...) {
      out <- rep(FALSE, length(zScore)); out[2] <- TRUE
      list(outlier = out, diagnostics = data.frame(outlier = out))
    },
    ldMismatchQc = function(zScore, R, ...) {
      seenN <<- length(zScore)
      data.frame(outlier = rep(FALSE, length(zScore)))
    },
    .package = "pecotmr"
  )
  result <- suppressMessages(summaryStatsQc(
    rssInput = rssInput, ldData = td$ldData,
    zMismatchQc = "dentist", alleleFlipKriging = TRUE, impute = FALSE
  ))
  expect_true(is(result, "QcResult"))
  expect_equal(seenN, nrow(td$sumstats) - 1)   # kriging removed one before dentist
  expect_equal(getOutlierNumber(result), 1)
})
