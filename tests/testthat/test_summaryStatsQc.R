context("summaryStatsQc (with mocked MungeSumstats)")

# NOTE
# ----
# `.runMungeSumstatsFilter` wraps MungeSumstats::format_sumstats which needs
# a real dbSNP reference panel (multi-GB download). To exercise the QC chain
# in a unit test we mock that helper so it just returns the input data.frame
# unchanged, recording a "no variants dropped" audit record. The pecotmr-
# native steps (.applySkipRegion, .matchAgainstSketch, .applyPipScreen,
# .applyLdMismatchQcToEntry) all run for real on the synthetic fixture.

# ===========================================================================
# Fixture builders
# ===========================================================================

.ssQ_makeHandle <- function(snp_n = 8L, n_samples = 60L) {
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

.ssQ_makeEntryGr <- function(snp_ids = paste0("rs", 1:4),
                             positions = c(100L, 200L, 300L, 400L)) {
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", length(snp_ids)),
    ranges = IRanges::IRanges(start = positions, width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = snp_ids,
    A1  = rep("A", length(snp_ids)),
    A2  = rep("G", length(snp_ids)),
    Z   = seq(1.0, by = 0.5, length.out = length(snp_ids)),
    N   = rep(1000L, length(snp_ids)))
  gr
}

.ssQ_makeGwasSumStats <- function(snp_ids = paste0("rs", 1:4),
                                  positions = c(100L, 200L, 300L, 400L),
                                  study = "g1") {
  GwasSumStats(
    study    = study,
    entry    = list(.ssQ_makeEntryGr(snp_ids, positions)),
    genome   = "hg19",
    ldSketch = .ssQ_makeHandle())
}

.ssQ_mockMunge <- function(drop = 0L) {
  # Mock that pretends MungeSumstats validated the input and returned the
  # same data.frame, dropping `drop` rows.
  function(df, refGenome, useDbsnpRefCheck, removeIndels,
           removeStrandAmbiguous, mafCutoff, infoCutoff, nCutoff,
           convertRefGenome, mungeSumstatsArgs) {
    keep <- if (drop > 0L && drop < nrow(df))
      seq_len(nrow(df) - drop)
    else
      seq_len(nrow(df))
    list(df = df[keep, , drop = FALSE],
         droppedNVariants = nrow(df) - length(keep))
  }
}

.ssQ_mockExtractor <- function(seed = 13, n_samples = 60L) {
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
# summaryStatsQc: input-type validation
# ===========================================================================

test_that("summaryStatsQc: rejects non-SumStats input", {
  expect_error(summaryStatsQc("not_a_sumstats"),
               "requires a QtlSumStats or GwasSumStats input")
})

test_that("summaryStatsQc: mafCutoff > 0 with no MAF column errors", {
  ss <- .ssQ_makeGwasSumStats()
  expect_error(summaryStatsQc(ss, mafCutoff = 0.05),
               "mafCutoff > 0 requires every entry to carry a MAF column")
})

test_that("summaryStatsQc: infoCutoff > 0 with no INFO column errors", {
  ss <- .ssQ_makeGwasSumStats()
  expect_error(summaryStatsQc(ss, infoCutoff = 0.5),
               "infoCutoff > 0 requires every entry to carry an INFO column")
})

# ===========================================================================
# summaryStatsQc: end-to-end with mocked MungeSumstats
# ===========================================================================

test_that("summaryStatsQc: vanilla run populates qcInfo and returns a GwasSumStats", {
  ss <- .ssQ_makeGwasSumStats()
  local_mocked_bindings(
    .runMungeSumstatsFilter = .ssQ_mockMunge(),
    .package = "pecotmr")
  res <- summaryStatsQc(ss)
  expect_s4_class(res, "GwasSumStats")
  qc <- getQcInfo(res)
  expect_true(length(qc) > 0L)
  expect_true("options" %in% names(qc))
  expect_true("entryAudit" %in% names(qc))
  expect_equal(length(qc$entryAudit), nrow(ss))
  # Per-entry audit records variantsIn / variantsOut / mungeSumstatsDropped.
  ea <- qc$entryAudit[[1L]]
  expect_equal(ea$variantsIn, 4L)
  expect_equal(ea$variantsOut, 4L)
  expect_equal(ea$mungeSumstatsDropped, 0L)
})

test_that("summaryStatsQc: keepVariants subsets each entry and records the drop", {
  ss <- .ssQ_makeGwasSumStats()
  local_mocked_bindings(
    .runMungeSumstatsFilter = .ssQ_mockMunge(),
    .package = "pecotmr")
  res <- summaryStatsQc(ss, keepVariants = c("rs1", "rs3"))
  ea <- getQcInfo(res)$entryAudit[[1L]]
  expect_equal(ea$keepVariantsDropped, 2L)
  expect_equal(ea$variantsOut, 2L)
})

test_that("summaryStatsQc: skipRegion drops overlapping variants", {
  ss <- .ssQ_makeGwasSumStats()
  local_mocked_bindings(
    .runMungeSumstatsFilter = .ssQ_mockMunge(),
    .package = "pecotmr")
  res <- summaryStatsQc(ss, skipRegion = "chr1:50-150")
  ea <- getQcInfo(res)$entryAudit[[1L]]
  expect_equal(ea$skipRegionDropped, 1L)  # rs1 at pos 100 is dropped
})

test_that("summaryStatsQc: PIP screen triggers when no variant has signal", {
  # Build an entry with weak signal so the SER PIP screen tags everything
  # below the threshold.
  gr <- .ssQ_makeEntryGr()
  S4Vectors::mcols(gr)$Z <- rep(0.1, length(gr))
  ss <- GwasSumStats(study = "g1", entry = list(gr), genome = "hg19",
                      ldSketch = .ssQ_makeHandle())
  local_mocked_bindings(
    .runMungeSumstatsFilter = .ssQ_mockMunge(),
    .package = "pecotmr")
  res <- summaryStatsQc(ss, pipCutoffToSkip = 0.99)
  ea <- getQcInfo(res)$entryAudit[[1L]]
  expect_true(isTRUE(ea$pipScreenSkipped))
  expect_match(ea$pipScreenReason, "no signals above PIP threshold")
  expect_equal(length(res$entry[[1L]]), 0L)
})

test_that("summaryStatsQc: early-exit records when fewer than 2 variants remain pre-harmonization", {
  ss <- .ssQ_makeGwasSumStats()
  local_mocked_bindings(
    # Mock keeps only 1 row of the input
    .runMungeSumstatsFilter = .ssQ_mockMunge(drop = 3L),
    .package = "pecotmr")
  res <- summaryStatsQc(ss)
  ea <- getQcInfo(res)$entryAudit[[1L]]
  expect_match(ea$earlyExit, "fewer than two variants")
})

test_that("summaryStatsQc: harmonized variants count is recorded", {
  ss <- .ssQ_makeGwasSumStats()
  local_mocked_bindings(
    .runMungeSumstatsFilter = .ssQ_mockMunge(),
    .package = "pecotmr")
  res <- summaryStatsQc(ss)
  ea <- getQcInfo(res)$entryAudit[[1L]]
  expect_equal(ea$matchedAgainstSketch, 4L)
})

test_that("summaryStatsQc: options block records the curated knobs", {
  ss <- .ssQ_makeGwasSumStats()
  local_mocked_bindings(
    .runMungeSumstatsFilter = .ssQ_mockMunge(),
    .package = "pecotmr")
  res <- summaryStatsQc(ss, removeIndels = TRUE, removeStrandAmbiguous = FALSE,
                        nCutoff = 10)
  opts <- getQcInfo(res)$options
  expect_true(opts$removeIndels)
  expect_false(opts$removeStrandAmbiguous)
  expect_equal(opts$nCutoff, 10)
})

test_that("summaryStatsQc: round-trips QtlSumStats inputs", {
  gr <- .ssQ_makeEntryGr()
  ss <- QtlSumStats(study = "s1", context = "c1", trait = "t1",
                     entry = list(gr), genome = "hg19",
                     ldSketch = .ssQ_makeHandle())
  local_mocked_bindings(
    .runMungeSumstatsFilter = .ssQ_mockMunge(),
    .package = "pecotmr")
  res <- summaryStatsQc(ss)
  expect_s4_class(res, "QtlSumStats")
  expect_equal(length(getQcInfo(res)$entryAudit), 1L)
})

# ===========================================================================
# summaryStatsQc with LD-mismatch QC enabled (mocked extractor)
# ===========================================================================

test_that("summaryStatsQc: zMismatchQc = 'dentist' walks the LD-mismatch branch", {
  ss <- .ssQ_makeGwasSumStats(snp_ids = paste0("rs", 1:8),
                              positions = seq(100L, by = 100L, length.out = 8L))
  local_mocked_bindings(
    .runMungeSumstatsFilter = .ssQ_mockMunge(),
    extractBlockGenotypes   = .ssQ_mockExtractor(),
    .package = "pecotmr")
  res <- suppressWarnings(summaryStatsQc(ss, zMismatchQc = "dentist"))
  ea <- getQcInfo(res)$entryAudit[[1L]]
  expect_equal(ea$ldMismatchMethod, "dentist")
  expect_true("ldMismatchOutliersDropped" %in% names(ea))
})

# ===========================================================================
# summaryStatsQc with impute = TRUE: exercise the RAISS branch
# ===========================================================================

test_that("summaryStatsQc: impute = TRUE invokes RAISS and records the audit counts", {
  # Build a sketch panel with 8 variants and a GWAS entry covering only the
  # first 4 — RAISS is asked to impute the missing 4.
  full_snp_ids <- paste0("rs", 1:8)
  full_positions <- seq(100L, by = 100L, length.out = 8L)
  ss <- GwasSumStats(
    study  = "g1",
    entry  = list(.ssQ_makeEntryGr(
                    snp_ids   = full_snp_ids[1:4],
                    positions = full_positions[1:4])),
    genome = "hg19",
    ldSketch = .ssQ_makeHandle(snp_n = 8L, n_samples = 60L))

  local_mocked_bindings(
    .runMungeSumstatsFilter = .ssQ_mockMunge(),
    extractBlockGenotypes   = .ssQ_mockExtractor(),
    raiss = function(refPanel, knownZscores, genotypeMatrix, ...) {
      # Pretend RAISS imputed two of the missing panel variants (rs5, rs6)
      # with synthetic z-scores.
      added <- refPanel[refPanel$variant_id %in% c("rs5", "rs6"), , drop = FALSE]
      added$z <- c(1.5, -2.0)
      added$n <- c(1000, 1000)
      list(resultFilter = rbind(knownZscores, added))
    },
    .package = "pecotmr")
  res <- summaryStatsQc(ss, impute = TRUE)
  ea <- getQcInfo(res)$entryAudit[[1L]]
  expect_equal(ea$raissTotalVariants, 6L)
  expect_equal(ea$raissImputedVariants, 2L)
})

test_that("summaryStatsQc: impute = TRUE with raiss returning NULL records 0 imputed", {
  full_snp_ids <- paste0("rs", 1:8)
  full_positions <- seq(100L, by = 100L, length.out = 8L)
  ss <- GwasSumStats(
    study  = "g1",
    entry  = list(.ssQ_makeEntryGr(
                    snp_ids   = full_snp_ids[1:4],
                    positions = full_positions[1:4])),
    genome = "hg19",
    ldSketch = .ssQ_makeHandle(snp_n = 8L, n_samples = 60L))

  local_mocked_bindings(
    .runMungeSumstatsFilter = .ssQ_mockMunge(),
    extractBlockGenotypes   = .ssQ_mockExtractor(),
    raiss = function(...) NULL,
    .package = "pecotmr")
  res <- summaryStatsQc(ss, impute = TRUE)
  ea <- getQcInfo(res)$entryAudit[[1L]]
  expect_equal(ea$raissImputedVariants, 0L)
})
