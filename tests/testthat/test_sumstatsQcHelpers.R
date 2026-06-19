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
