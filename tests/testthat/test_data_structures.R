# Tests for S4 data structure classes: LDData, RegionalData,
# FineMappingResult, TWASWeights

# =============================================================================
# LDData
# =============================================================================

test_that("LDData constructor works with correlation matrix", {
  R <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  rownames(R) <- colnames(R) <- c("chr1:100:A:G", "chr1:200:C:T")
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(100L, 200L), width = 1L)
  )
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    variant_id = c("chr1:100:A:G", "chr1:200:C:T"),
    A1 = c("G", "T"), A2 = c("A", "C")
  )
  bm <- data.frame(block_id = 1L, chrom = "1", block_start = 100L,
                    block_end = 200L, size = 2L, start_idx = 1L, end_idx = 2L)

  ld <- LDData(correlation = R, variants = gr, block_metadata = bm)
  expect_s4_class(ld, "LDData")
  expect_false(hasGenotypes(ld))
  expect_true(is.matrix(getCorrelation(ld)))
  expect_equal(getVariantIds(ld), c("chr1:100:A:G", "chr1:200:C:T"))
  expect_equal(nrow(getBlockMetadata(ld)), 1L)
  expect_null(getGenotypes(ld))
})

test_that("LDData validation rejects empty variants", {
  R <- diag(2)
  gr <- GenomicRanges::GRanges()
  expect_error(
    LDData(correlation = R, variants = gr,
           block_metadata = data.frame()),
    "must not be empty"
  )
})

test_that("LDData validation rejects NULL correlation and handle", {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 100L, width = 1L)
  )
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    variant_id = "chr1:100:A:G", A1 = "G", A2 = "A"
  )
  expect_error(
    LDData(correlation = NULL, genotype_handle = NULL,
           variants = gr, block_metadata = data.frame()),
    "At least one"
  )
})

test_that("LDData show method works", {
  R <- diag(3)
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 3),
    ranges = IRanges::IRanges(start = c(100L, 200L, 300L), width = 1L)
  )
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    variant_id = paste0("v", 1:3), A1 = rep("A", 3), A2 = rep("G", 3)
  )
  ld <- LDData(correlation = R, variants = gr,
               block_metadata = data.frame())
  expect_output(show(ld), "LDData: 3 variants")
})

test_that("LDData supports block-diagonal correlation", {
  R1 <- diag(2)
  R2 <- diag(3)
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 5),
    ranges = IRanges::IRanges(start = seq(100L, 500L, 100L), width = 1L)
  )
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    variant_id = paste0("v", 1:5), A1 = rep("A", 5), A2 = rep("G", 5)
  )
  ld <- LDData(correlation = list(R1, R2), variants = gr,
               block_metadata = data.frame())
  corr <- getCorrelation(ld)
  expect_true(is.list(corr))
  expect_equal(length(corr), 2)
})

test_that("ld_data_to_list converts LDData to legacy format", {
  R <- diag(2)
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(100L, 200L), width = 1L)
  )
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    variant_id = c("v1", "v2"), A1 = c("A", "C"), A2 = c("G", "T")
  )
  ld <- LDData(correlation = R, variants = gr,
               block_metadata = data.frame(block_id = 1L))
  legacy <- pecotmr:::ld_data_to_list(ld)
  expect_true(is.list(legacy))
  expect_true("LD_matrix" %in% names(legacy))
  expect_true("LD_variants" %in% names(legacy))
  expect_true("ref_panel" %in% names(legacy))
  expect_false(legacy$is_genotype)
  expect_equal(legacy$LD_variants, c("v1", "v2"))
})

test_that(".ref_panel_to_granges builds GRanges from data.frame", {
  rp <- data.frame(
    chrom = c("1", "1"), pos = c(100L, 200L),
    A1 = c("G", "T"), A2 = c("A", "C"),
    variant_id = c("v1", "v2"),
    allele_freq = c(0.3, 0.7),
    stringsAsFactors = FALSE
  )
  gr <- pecotmr:::.ref_panel_to_granges(rp)
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 2)
  expect_equal(S4Vectors::mcols(gr)$variant_id, c("v1", "v2"))
  expect_equal(S4Vectors::mcols(gr)$allele_freq, c(0.3, 0.7))
})

# =============================================================================
# RegionalData
# =============================================================================

test_that("RegionalData constructor and lazy residuals work", {
  set.seed(42)
  n <- 50; p <- 10; k <- 2
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("v", 1:p)
  rownames(X) <- paste0("s", 1:n)
  Y1 <- matrix(rnorm(n), n, 1, dimnames = list(paste0("s", 1:n), "trait1"))
  Y2 <- matrix(rnorm(n), n, 1, dimnames = list(paste0("s", 1:n), "trait2"))
  C <- matrix(rnorm(n * 3), n, 3)
  rownames(C) <- paste0("s", 1:n)

  rd <- RegionalData(
    genotype_matrix = X,
    phenotypes = list(cond1 = Y1, cond2 = Y2),
    covariates = list(cond1 = C, cond2 = C),
    scale_residuals = FALSE,
    maf = list(cond1 = runif(p, 0.05, 0.5), cond2 = runif(p, 0.05, 0.5))
  )
  expect_s4_class(rd, "RegionalData")

  # Lazy residuals
  rx <- getResidualX(rd, 1L)
  expect_true(is.matrix(rx))
  expect_equal(ncol(rx), p)
  expect_equal(nrow(rx), n)

  ry <- getResidualY(rd, 1L)
  expect_true(is.matrix(ry))
  expect_equal(nrow(ry), n)

  # Scalars when not scaling
  expect_equal(getResidualXScalar(rd, 1L), rep(1, p))
  expect_equal(getResidualYScalar(rd, 1L), 1)
})

test_that("RegionalData with scale_residuals=TRUE computes scaled residuals", {
  set.seed(1)
  n <- 30; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("v", 1:p)
  rownames(X) <- paste0("s", 1:n)
  Y <- matrix(rnorm(n), n, 1, dimnames = list(paste0("s", 1:n), "y"))
  C <- matrix(rnorm(n * 2), n, 2)
  rownames(C) <- paste0("s", 1:n)

  rd <- RegionalData(
    genotype_matrix = X,
    phenotypes = list(cond1 = Y),
    covariates = list(cond1 = C),
    scale_residuals = TRUE,
    maf = list(cond1 = runif(p))
  )

  rx <- getResidualX(rd, 1L)
  # Scaled columns should have unit variance (approximately)
  col_sds <- apply(rx, 2, sd, na.rm = TRUE)
  expect_true(all(abs(col_sds - 1) < 0.1 | is.na(col_sds)))
})

test_that("RegionalData validation rejects mismatched phenotypes/covariates", {
  X <- matrix(1, 5, 3)
  colnames(X) <- paste0("v", 1:3)
  rownames(X) <- paste0("s", 1:5)
  Y <- matrix(1, 5, 1, dimnames = list(paste0("s", 1:5), "y"))
  C <- matrix(1, 5, 2)
  rownames(C) <- paste0("s", 1:5)

  expect_error(
    RegionalData(
      genotype_matrix = X,
      phenotypes = list(cond1 = Y, cond2 = Y),
      covariates = list(cond1 = C)
    ),
    "same length"
  )
})

test_that("RegionalData show method works", {
  X <- matrix(1, 5, 3)
  colnames(X) <- paste0("v", 1:3)
  rownames(X) <- paste0("s", 1:5)
  Y <- matrix(1, 5, 1, dimnames = list(paste0("s", 1:5), "y"))
  C <- matrix(1, 5, 2)
  rownames(C) <- paste0("s", 1:5)

  rd <- RegionalData(
    genotype_matrix = X,
    phenotypes = list(cond1 = Y),
    covariates = list(cond1 = C)
  )
  expect_output(show(rd), "RegionalData: 1 conditions")
})

# =============================================================================
# FineMappingResult
# =============================================================================

test_that("FineMappingResult constructor and accessors work", {
  tl <- data.frame(
    variant_id = c("v1", "v2", "v3"),
    method = rep("susie_rss", 3),
    pip = c(0.9, 0.05, 0.8),
    cs = c(1L, 0L, 2L),
    stringsAsFactors = FALSE
  )
  fm <- FineMappingResult(
    variant_names = c("v1", "v2", "v3"),
    trimmed_fit = list(alpha = matrix(0, 2, 3)),
    top_loci = tl,
    method = "susie_rss"
  )
  expect_s4_class(fm, "FineMappingResult")

  pip <- getPIP(fm)
  expect_equal(length(pip), 3)
  expect_equal(unname(pip[1]), 0.9)

  cs_df <- getCS(fm)
  expect_equal(nrow(cs_df), 2)  # v1 and v3 have cs > 0
})

test_that("FineMappingResult validation rejects missing method", {
  expect_error(
    FineMappingResult(
      variant_names = "v1",
      trimmed_fit = list(),
      top_loci = data.frame(variant_id = "v1", method = "x"),
      method = character(0)
    ),
    "single character"
  )
})

test_that("FineMappingResult show method works", {
  tl <- data.frame(
    variant_id = c("v1", "v2"),
    method = rep("susie_rss", 2),
    pip = c(0.9, 0.1),
    cs = c(1L, 0L),
    stringsAsFactors = FALSE
  )
  fm <- FineMappingResult(
    variant_names = c("v1", "v2"),
    trimmed_fit = list(),
    top_loci = tl,
    method = "susie_rss"
  )
  expect_output(show(fm), "FineMappingResult.*susie_rss.*2 variants.*1 credible")
})

# =============================================================================
# TWASWeights
# =============================================================================

test_that("TWASWeights constructor and accessors work", {
  w1 <- matrix(rnorm(10), 5, 2)
  w2 <- matrix(rnorm(10), 5, 2)
  rownames(w1) <- rownames(w2) <- paste0("v", 1:5)

  tw <- TWASWeights(
    weights = list(lasso = w1, enet = w2),
    variant_ids = paste0("v", 1:5)
  )
  expect_s4_class(tw, "TWASWeights")
  expect_equal(tw@methods, c("lasso", "enet"))

  all_w <- getWeights(tw)
  expect_true(is.list(all_w))
  expect_equal(length(all_w), 2)

  w_lasso <- getWeights(tw, "lasso")
  expect_true(is.matrix(w_lasso))
  expect_equal(nrow(w_lasso), 5)
})

test_that("TWASWeights validation rejects dimension mismatch", {
  w1 <- matrix(0, 3, 1)
  expect_error(
    TWASWeights(
      weights = list(method1 = w1),
      variant_ids = paste0("v", 1:5)
    ),
    "3 rows"
  )
})

test_that("TWASWeights show method works", {
  tw <- TWASWeights(
    weights = list(a = matrix(0, 2, 1), b = matrix(0, 2, 1)),
    variant_ids = c("v1", "v2")
  )
  expect_output(show(tw), "TWASWeights: 2 methods, 2 variants")
})

test_that("getWeights errors on unknown method", {
  tw <- TWASWeights(
    weights = list(a = matrix(0, 2, 1)),
    variant_ids = c("v1", "v2")
  )
  expect_error(getWeights(tw, "nonexistent"), "not found")
})

# =============================================================================
# top_loci_to_granges
# =============================================================================

test_that("top_loci_to_granges converts data.frame to GRanges", {
  tl <- data.frame(
    variant_id = c("1:100:A:G", "1:200:C:T"),
    pip = c(0.9, 0.1),
    cs = c(1L, 0L),
    method = "susie",
    stringsAsFactors = FALSE
  )
  gr <- top_loci_to_granges(tl)
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 2)
  expect_equal(S4Vectors::mcols(gr)$pip, c(0.9, 0.1))
})

test_that("top_loci_to_granges handles empty input", {
  gr <- top_loci_to_granges(NULL)
  expect_equal(length(gr), 0)
  gr2 <- top_loci_to_granges(data.frame())
  expect_equal(length(gr2), 0)
})

# =============================================================================
# extractBlockGenotypes returns RSE
# =============================================================================

test_that("extractBlockGenotypes returns SummarizedExperiment", {
  skip_if_not_installed("pgenlibr")
  # Find test plink2 data
  test_dir <- system.file("testdata", package = "pecotmr")
  if (test_dir == "" || !dir.exists(test_dir)) skip("No testdata directory")
  pgen_files <- list.files(test_dir, pattern = "\\.pgen$", full.names = TRUE)
  if (length(pgen_files) == 0) skip("No pgen files in testdata")

  stem <- tools::file_path_sans_ext(pgen_files[1])
  handle <- readGenotypes(stem, format = "plink2")
  n_snps <- nrow(handle@snp_info)
  if (n_snps == 0) skip("No SNPs in handle")

  rse <- extractBlockGenotypes(handle, seq_len(min(5, n_snps)))
  expect_s4_class(rse, "SummarizedExperiment")
  expect_true("dosage" %in% SummarizedExperiment::assayNames(rse))
  dosage <- SummarizedExperiment::assay(rse, "dosage")
  # Bioc convention: variants x samples
  expect_equal(nrow(dosage), min(5, n_snps))
  expect_equal(ncol(dosage), handle@n_samples)
  # rowRanges should have variant info
  rr <- SummarizedExperiment::rowRanges(rse)
  expect_true("A1" %in% names(S4Vectors::mcols(rr)))
  expect_true("A2" %in% names(S4Vectors::mcols(rr)))
})

# =============================================================================
# getLBF on FineMappingResult
# =============================================================================

test_that("getLBF returns data.frame with variant_id and L columns", {
  # 2 effects x 4 variants
  lbf_mat <- matrix(c(1.1, 2.2, 3.3, 4.4,
                       5.5, 6.6, 7.7, 8.8), nrow = 2, byrow = TRUE)
  pip_vec <- c(v1 = 0.9, v2 = 0.1, v3 = 0.8, v4 = 0.05)

  fm <- FineMappingResult(
    variant_names = c("v1", "v2", "v3", "v4"),
    trimmed_fit = list(lbf_variable = lbf_mat, pip = pip_vec),
    top_loci = data.frame(variant_id = names(pip_vec),
                          method = rep("susie_rss", 4),
                          stringsAsFactors = FALSE),
    method = "susie_rss"
  )
  result <- getLBF(fm)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 4)
  expect_true("variant_id" %in% names(result))
  expect_true(all(c("L1", "L2") %in% names(result)))
  expect_equal(result$variant_id, c("v1", "v2", "v3", "v4"))
  # lbf_variable rows are effects, cols are variants; transposed so L1 column
  # corresponds to the first effect row

  expect_equal(result$L1, c(1.1, 2.2, 3.3, 4.4))
  expect_equal(result$L2, c(5.5, 6.6, 7.7, 8.8))
})

test_that("getLBF returns empty data.frame when trimmed_fit is NULL", {
  fm <- FineMappingResult(
    variant_names = c("v1", "v2"),
    trimmed_fit = NULL,
    top_loci = data.frame(variant_id = c("v1", "v2"),
                          method = rep("susie_rss", 2),
                          stringsAsFactors = FALSE),
    method = "susie_rss"
  )
  result <- getLBF(fm)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
  expect_equal(ncol(result), 0)
})

test_that("getLBF returns empty data.frame when lbf_variable is absent", {
  fm <- FineMappingResult(
    variant_names = c("v1", "v2"),
    trimmed_fit = list(pip = c(v1 = 0.5, v2 = 0.5)),
    top_loci = data.frame(variant_id = c("v1", "v2"),
                          method = rep("susie_rss", 2),
                          stringsAsFactors = FALSE),
    method = "susie_rss"
  )
  result <- getLBF(fm)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
  expect_equal(ncol(result), 0)
})

test_that("getLBF falls back to pip names when variant_names is empty", {
  lbf_mat <- matrix(c(1, 2, 3, 4), nrow = 1)
  fm <- FineMappingResult(
    variant_names = character(0),
    trimmed_fit = list(lbf_variable = lbf_mat,
                       pip = c(a = 0.9, b = 0.5, c = 0.3, d = 0.1)),
    top_loci = data.frame(variant_id = character(0),
                          method = character(0),
                          stringsAsFactors = FALSE),
    method = "susie_rss"
  )
  result <- getLBF(fm)
  expect_equal(result$variant_id, c("a", "b", "c", "d"))
  expect_equal(result$L1, c(1, 2, 3, 4))
})

# =============================================================================
# getEffects on FineMappingResult
# =============================================================================

test_that("getEffects returns per-effect summary with correct columns", {
  # 3 effects, 5 variants
  purity_mat <- matrix(
    c(0.8, 0.9, 0.85,   # min.abs.corr column
      0.85, 0.95, 0.90), # mean.abs.corr column
    nrow = 3, ncol = 2
  )
  rownames(purity_mat) <- c("L1", "L2", "L3")

  fm <- FineMappingResult(
    variant_names = c("v1", "v2", "v3", "v4", "v5"),
    trimmed_fit = list(
      V = c(0.01, 0.02, 0.03),
      lbf = c(10.5, 20.3, 5.1),
      alpha = matrix(0, nrow = 3, ncol = 5),
      sets = list(
        cs = list(L1 = c(1L, 3L), L2 = c(2L, 4L, 5L)),
        coverage = c(0.95, 0.99),
        purity = purity_mat
      )
    ),
    top_loci = data.frame(variant_id = paste0("v", 1:5),
                          method = rep("susie_rss", 5),
                          stringsAsFactors = FALSE),
    method = "susie_rss"
  )

  result <- getEffects(fm)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)
  expected_cols <- c("effect_id", "V", "cs_log10bf", "cs_min_r2",
                     "cs_avg_r2", "coverage", "cs")
  expect_true(all(expected_cols %in% names(result)))

  # effect_ids

  expect_equal(result$effect_id, c("L1", "L2", "L3"))
  # V
  expect_equal(result$V, c(0.01, 0.02, 0.03))
  # cs_log10bf
  expect_equal(result$cs_log10bf, c(10.5, 20.3, 5.1))
})

test_that("getEffects cs column has semicolon-separated variant names", {
  purity_mat <- matrix(
    c(0.8, 0.9, 0.85, 0.95), nrow = 2, ncol = 2
  )
  rownames(purity_mat) <- c("L1", "L2")

  fm <- FineMappingResult(
    variant_names = c("v1", "v2", "v3", "v4"),
    trimmed_fit = list(
      V = c(0.01, 0.02),
      lbf = c(10.0, 20.0),
      alpha = matrix(0, nrow = 2, ncol = 4),
      sets = list(
        cs = list(L1 = c(1L, 3L), L2 = c(2L, 4L)),
        coverage = c(0.95, 0.99),
        purity = purity_mat
      )
    ),
    top_loci = data.frame(variant_id = paste0("v", 1:4),
                          method = rep("susie_rss", 4),
                          stringsAsFactors = FALSE),
    method = "susie_rss"
  )

  result <- getEffects(fm)
  expect_equal(result$cs[1], "v1;v3")
  expect_equal(result$cs[2], "v2;v4")
  expect_equal(result$coverage[1], 0.95)
  expect_equal(result$coverage[2], 0.99)
  expect_equal(result$cs_min_r2[1], 0.8)
  expect_equal(result$cs_avg_r2[1], 0.85)
  expect_equal(result$cs_min_r2[2], 0.9)
  expect_equal(result$cs_avg_r2[2], 0.95)
})

test_that("getEffects reports None for effects without CS", {
  fm <- FineMappingResult(
    variant_names = c("v1", "v2", "v3"),
    trimmed_fit = list(
      V = c(0.01, 0.02),
      lbf = c(5.0, 3.0),
      alpha = matrix(0, nrow = 2, ncol = 3),
      sets = list(cs = NULL, coverage = NULL, purity = NULL)
    ),
    top_loci = data.frame(variant_id = paste0("v", 1:3),
                          method = rep("susie_rss", 3),
                          stringsAsFactors = FALSE),
    method = "susie_rss"
  )

  result <- getEffects(fm)
  expect_equal(nrow(result), 2)
  expect_true(all(result$cs == "None"))
  expect_true(all(result$coverage == 0))
  expect_true(all(result$cs_min_r2 == 0))
})

test_that("getEffects returns empty data.frame when trimmed_fit is NULL", {
  fm <- FineMappingResult(
    variant_names = c("v1", "v2"),
    trimmed_fit = NULL,
    top_loci = data.frame(variant_id = c("v1", "v2"),
                          method = rep("susie_rss", 2),
                          stringsAsFactors = FALSE),
    method = "susie_rss"
  )
  result <- getEffects(fm)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
  expect_equal(ncol(result), 0)
})

test_that("getEffects returns empty data.frame when no V or alpha", {
  fm <- FineMappingResult(
    variant_names = c("v1", "v2"),
    trimmed_fit = list(pip = c(v1 = 0.5, v2 = 0.5)),
    top_loci = data.frame(variant_id = c("v1", "v2"),
                          method = rep("susie_rss", 2),
                          stringsAsFactors = FALSE),
    method = "susie_rss"
  )
  result <- getEffects(fm)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
  expect_equal(ncol(result), 0)
})
