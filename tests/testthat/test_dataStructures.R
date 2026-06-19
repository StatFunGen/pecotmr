# Tests for S4 data structure classes: LdData, RegionalData,
# FineMappingResult, TwasWeights

# =============================================================================
# LdData
# =============================================================================

test_that("LdData constructor works with correlation matrix", {
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
  bm <- data.frame(blockId = 1L, chrom = "1", blockStart = 100L,
                    blockEnd = 200L, size = 2L, startIdx = 1L, endIdx = 2L)

  ld <- LdData(correlation = R, variants = gr, blockMetadata = bm)
  expect_s4_class(ld, "LdData")
  expect_false(hasGenotypes(ld))
  expect_true(is.matrix(getCorrelation(ld)))
  expect_equal(getVariantIds(ld), c("chr1:100:A:G", "chr1:200:C:T"))
  expect_equal(nrow(getBlockMetadata(ld)), 1L)
  expect_null(getGenotypes(ld))
})

test_that("LdData validation rejects empty variants", {
  R <- diag(2)
  gr <- GenomicRanges::GRanges()
  expect_error(
    LdData(correlation = R, variants = gr,
           blockMetadata = data.frame()),
    "must not be empty"
  )
})

test_that("LdData validation rejects NULL correlation and handle", {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 100L, width = 1L)
  )
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    variant_id = "chr1:100:A:G", A1 = "G", A2 = "A"
  )
  expect_error(
    LdData(correlation = NULL, genotypeHandle = NULL,
           variants = gr, blockMetadata = data.frame()),
    "At least one"
  )
})

test_that("LdData show method works", {
  R <- diag(3)
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 3),
    ranges = IRanges::IRanges(start = c(100L, 200L, 300L), width = 1L)
  )
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    variant_id = paste0("v", 1:3), A1 = rep("A", 3), A2 = rep("G", 3)
  )
  ld <- LdData(correlation = R, variants = gr,
               blockMetadata = data.frame())
  expect_output(show(ld), "LdData: 3 variants")
})

test_that("LdData supports block-diagonal correlation", {
  R1 <- diag(2)
  R2 <- diag(3)
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 5),
    ranges = IRanges::IRanges(start = seq(100L, 500L, 100L), width = 1L)
  )
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    variant_id = paste0("v", 1:5), A1 = rep("A", 5), A2 = rep("G", 5)
  )
  ld <- LdData(correlation = list(R1, R2), variants = gr,
               blockMetadata = data.frame())
  corr <- getCorrelation(ld)
  expect_true(is.list(corr))
  expect_equal(length(corr), 2)
})

test_that("LdData S4 accessors return correct data", {
  R <- diag(2)
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(100L, 200L), width = 1L)
  )
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    variant_id = c("v1", "v2"), A1 = c("A", "C"), A2 = c("G", "T")
  )
  ld <- LdData(correlation = R, variants = gr,
               blockMetadata = data.frame(blockId = 1L))
  expect_equal(getCorrelation(ld), R)
  expect_equal(getVariantIds(ld), c("v1", "v2"))
  expect_false(hasGenotypes(ld))
  rp <- getRefPanel(ld)
  expect_true(is.data.frame(rp))
  expect_true("variant_id" %in% names(rp))
  expect_equal(rp$variant_id, c("v1", "v2"))
})

test_that(".refPanelToGranges builds GRanges from data.frame", {
  rp <- data.frame(
    chrom = c("1", "1"), pos = c(100L, 200L),
    A1 = c("G", "T"), A2 = c("A", "C"),
    variant_id = c("v1", "v2"),
    allele_freq = c(0.3, 0.7),
    stringsAsFactors = FALSE
  )
  gr <- pecotmr:::.refPanelToGranges(rp)
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 2)
  expect_equal(S4Vectors::mcols(gr)$variant_id, c("v1", "v2"))
  expect_equal(S4Vectors::mcols(gr)$allele_freq, c(0.3, 0.7))
})
# =============================================================================
# topLociToGranges
# =============================================================================

test_that("topLociToGranges converts data.frame to GRanges", {
  tl <- data.frame(
    variant_id = c("1:100:A:G", "1:200:C:T"),
    pip = c(0.9, 0.1),
    cs = c(1L, 0L),
    method = "susie",
    stringsAsFactors = FALSE
  )
  gr <- topLociToGranges(tl)
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 2)
  expect_equal(S4Vectors::mcols(gr)$pip, c(0.9, 0.1))
})

test_that("topLociToGranges handles empty input", {
  gr <- topLociToGranges(NULL)
  expect_equal(length(gr), 0)
  gr2 <- topLociToGranges(data.frame())
  expect_equal(length(gr2), 0)
})

# =============================================================================
# extractBlockGenotypes returns RSE
# =============================================================================

test_that("extractBlockGenotypes returns SummarizedExperiment", {
  skip_if_not_installed("pgenlibr")
  stem <- test_path("test_data", "test_variants")
  handle <- readGenotypes(stem, format = "plink2")
  n_snps <- nrow(handle@snpInfo)
  skip_if(n_snps == 0, "No SNPs in handle")

  rse <- extractBlockGenotypes(handle, seq_len(min(5L, n_snps)))
  expect_s4_class(rse, "SummarizedExperiment")
  expect_true("dosage" %in% SummarizedExperiment::assayNames(rse))
  dosage <- SummarizedExperiment::assay(rse, "dosage")
  # Bioc convention: variants x samples
  expect_equal(nrow(dosage), min(5L, n_snps))
  expect_equal(ncol(dosage), handle@nSamples)
  # rowRanges should have variant info
  rr <- SummarizedExperiment::rowRanges(rse)
  expect_true("A1" %in% names(S4Vectors::mcols(rr)))
  expect_true("A2" %in% names(S4Vectors::mcols(rr)))
})
