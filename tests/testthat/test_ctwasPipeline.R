context("ctwasPipeline")

# ===========================================================================
# Strategy: ctwas::ctwas_sumstats does the heavy work. We mock it to a
# function that just returns its inputs back, so we can verify how the
# pipeline assembles z_snp / weights / region_info / LD loader inputs.
# ===========================================================================

.ctp_makeHandle <- function(snp_n = 6L, n_samples = 30L) {
  new("GenotypeHandle",
    path = "/tmp/sketch.gds",
    format = "gds",
    snpInfo = data.frame(
      SNP = paste0("v", seq_len(snp_n)),
      CHR = rep("1", snp_n),
      BP  = seq(100L, by = 100L, length.out = snp_n),
      A1  = rep("A", snp_n),
      A2  = rep("G", snp_n),
      stringsAsFactors = FALSE),
    nSamples = n_samples,
    sampleIds = paste0("s", seq_len(n_samples)),
    pgenPtr = NULL)
}

.ctp_mockExtractor <- function(seed = 5, n_samples = 30L) {
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

.ctp_makeGwasSumstats <- function(qc = TRUE) {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = seq(100L, by = 100L, length.out = 6L),
                              width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = paste0("v", 1:6),
    A1  = rep("A", 6), A2  = rep("G", 6),
    Z   = rnorm(6), N = rep(1000L, 6))
  GwasSumStats(
    study    = "G1",
    entry    = list(gr),
    genome   = "hg19",
    ldSketch = .ctp_makeHandle(),
    qcInfo   = if (qc) list(step1 = "ok") else list())
}

.ctp_makeTwasWeights <- function() {
  e <- TwasWeightsEntry(
    variantIds = paste0("v", 1:5),
    weights    = c(0.1, 0.05, -0.2, 0.3, 0.0))
  TwasWeights(
    study    = "Q1", context = "c1", trait = "t1", method = "susie",
    entry    = list(e),
    ldSketch = .ctp_makeHandle())
}

# ===========================================================================
# Input validation
# ----------------------------------------------------------------------------
# The top of ctwasPipeline requires the (non-CRAN) `ctwas` package; without
# it the function errors out before any input-validation branch fires. Skip
# entry-point tests when ctwas isn't installed, but exercise the input-
# building helpers directly (they don't gate on ctwas).
# ===========================================================================

test_that("ctwasPipeline: rejects non-GwasSumStats gwasSumStats", {
  skip_if_not_installed("ctwas")
  expect_error(
    ctwasPipeline(gwasSumStats = "no",
                  twasWeights  = .ctp_makeTwasWeights()),
    "must be a GwasSumStats"
  )
})

test_that("ctwasPipeline: rejects un-QCd GwasSumStats", {
  skip_if_not_installed("ctwas")
  expect_error(
    ctwasPipeline(gwasSumStats = .ctp_makeGwasSumstats(qc = FALSE),
                  twasWeights  = .ctp_makeTwasWeights()),
    "has no QC record"
  )
})

test_that("ctwasPipeline: rejects missing twasWeights", {
  skip_if_not_installed("ctwas")
  expect_error(
    ctwasPipeline(gwasSumStats = .ctp_makeGwasSumstats()),
    "must be a TwasWeights"
  )
})

test_that("ctwasPipeline: rejects non-GRanges twasZ", {
  skip_if_not_installed("ctwas")
  expect_error(
    ctwasPipeline(gwasSumStats = .ctp_makeGwasSumstats(),
                  twasWeights  = .ctp_makeTwasWeights(),
                  twasZ        = "not a GRanges"),
    "must be a GRanges"
  )
})

test_that("ctwasPipeline: rejects bad regionId", {
  skip_if_not_installed("ctwas")
  expect_error(
    ctwasPipeline(gwasSumStats = .ctp_makeGwasSumstats(),
                  twasWeights  = .ctp_makeTwasWeights(),
                  regionId     = ""),
    "non-empty character"
  )
})

test_that("ctwasPipeline: rejects unknown groupPriorVarStructure value", {
  skip_if_not_installed("ctwas")
  expect_error(
    ctwasPipeline(gwasSumStats = .ctp_makeGwasSumstats(),
                  twasWeights  = .ctp_makeTwasWeights(),
                  groupPriorVarStructure = "bogus"),
    "'arg'"
  )
})

# ===========================================================================
# .ctwasRequireMatchingLdSketches
# ===========================================================================

test_that(".ctwasRequireMatchingLdSketches: NULL twas-side handle is allowed", {
  twNoLd <- TwasWeights(
    study = "Q1", context = "c1", trait = "t1", method = "susie",
    entry = list(TwasWeightsEntry(variantIds = paste0("v", 1:5),
                                   weights = rep(0.1, 5))),
    ldSketch = NULL)
  expect_silent(pecotmr:::.ctwasRequireMatchingLdSketches(
    twLd = NULL, gwasLd = .ctp_makeHandle()))
})

test_that(".ctwasRequireMatchingLdSketches: panel-size mismatch errors", {
  twLd  <- .ctp_makeHandle(snp_n = 5L)
  gwasLd <- .ctp_makeHandle(snp_n = 6L)
  expect_error(
    pecotmr:::.ctwasRequireMatchingLdSketches(twLd, gwasLd),
    "ldSketch panels differ in size"
  )
})

# ===========================================================================
# Input-building helpers
# ===========================================================================

test_that(".ctwasBuildZSnp: produces a flat data.frame keyed by SNP/study", {
  ss <- .ctp_makeGwasSumstats()
  df <- pecotmr:::.ctwasBuildZSnp(ss)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 6L)
  expect_setequal(colnames(df),
                  c("id", "chrom", "pos", "A1", "A2", "z", "study"))
  expect_setequal(df$id, paste0("v", 1:6))
  expect_setequal(unique(df$study), "G1")
})

test_that(".ctwasBuildSingleRegionInfo: pulls chrom + bp span from the ldSketch", {
  ri <- pecotmr:::.ctwasBuildSingleRegionInfo("block1", .ctp_makeHandle())
  expect_equal(ri$region_id, "block1")
  expect_equal(ri$chrom, 1L)
  expect_equal(ri$start, 100L)
  expect_equal(ri$stop, 600L)
})

test_that(".ctwasBuildSingleRegionInfo: multi-chromosome sketch errors", {
  h <- .ctp_makeHandle()
  h@snpInfo$CHR[1:3] <- "2"
  expect_error(
    pecotmr:::.ctwasBuildSingleRegionInfo("block1", h),
    "spans multiple chromosomes"
  )
})

test_that(".ctwasSnpInfoForBlock: returns id/chrom/pos/A1/A2", {
  df <- pecotmr:::.ctwasSnpInfoForBlock(.ctp_makeHandle())
  expect_setequal(colnames(df),
                  c("id", "chrom", "pos", "A1", "A2"))
})

test_that(".ctwasBuildWeights: keys per-tuple weights and stamps gene metadata", {
  tw <- .ctp_makeTwasWeights()
  wl <- pecotmr:::.ctwasBuildWeights(tw)
  expect_equal(length(wl), 1L)
  expect_equal(names(wl), "Q1|c1|t1|susie")
  expect_equal(wl[[1L]]$study, "Q1")
  expect_equal(wl[[1L]]$context, "c1")
  expect_equal(wl[[1L]]$gene_name, "t1")
  expect_equal(length(wl[[1L]]$id), 5L)
})

test_that(".ctwasBuildZGene: builds z_gene from a TWAS-Z GRanges", {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 200))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    qtlStudy = "Q1", context = "c1", trait = "t1", method = "susie",
    twasZ = 1.5)
  df <- pecotmr:::.ctwasBuildZGene(gr)
  expect_equal(nrow(df), 1L)
  expect_setequal(colnames(df),
                  c("id", "z", "type", "context", "gene_name",
                    "study", "method"))
  expect_equal(df$id, "Q1|c1|t1|susie")
})

# ===========================================================================
# LD loader / SNP-info loader closures
# ===========================================================================

test_that(".ctwasSingleBlockLdLoader: returns a function that produces an LD matrix", {
  local_mocked_bindings(extractBlockGenotypes = .ctp_mockExtractor(),
                        .package = "pecotmr")
  loader <- pecotmr:::.ctwasSingleBlockLdLoader(.ctp_makeHandle())
  R <- loader("ignored.cor")
  expect_true(is.matrix(R))
  expect_equal(dim(R), c(6L, 6L))
  expect_equal(rownames(R), paste0("v", 1:6))
})

test_that(".ctwasSingleBlockSnpInfoLoader: returns a function that produces snp_info", {
  loader <- pecotmr:::.ctwasSingleBlockSnpInfoLoader(.ctp_makeHandle())
  df <- loader("ignored.cor")
  expect_setequal(colnames(df), c("id", "chrom", "pos", "A1", "A2"))
})

# ===========================================================================
# End-to-end with mocked ctwas::ctwas_sumstats
# ===========================================================================

test_that("ctwasPipeline: assembles the documented input shape for ctwas_sumstats", {
  skip_if_not_installed("ctwas")
  ss <- .ctp_makeGwasSumstats()
  tw <- .ctp_makeTwasWeights()
  capturedArgs <- NULL
  local_mocked_bindings(
    ctwas_sumstats = function(...) {
      capturedArgs <<- list(...)
      list(susie_alpha_res = "mocked")
    },
    .package = "ctwas")
  out <- ctwasPipeline(gwasSumStats = ss, twasWeights = tw,
                       regionId = "myBlock")
  expect_equal(out$susie_alpha_res, "mocked")
  expect_equal(capturedArgs$region_info$region_id, "myBlock")
  expect_equal(capturedArgs$z_snp$id, paste0("v", 1:6))
  expect_equal(length(capturedArgs$weights), 1L)
  expect_equal(capturedArgs$L, 5L)
})

test_that("ctwasPipeline: forwards a twasZ argument as z_gene", {
  skip_if_not_installed("ctwas")
  ss <- .ctp_makeGwasSumstats()
  tw <- .ctp_makeTwasWeights()
  twasZ <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 200))
  S4Vectors::mcols(twasZ) <- S4Vectors::DataFrame(
    qtlStudy = "Q1", context = "c1", trait = "t1", method = "susie",
    twasZ = 1.5)
  capturedArgs <- NULL
  local_mocked_bindings(
    ctwas_sumstats = function(...) {
      capturedArgs <<- list(...)
      list(ok = TRUE)
    },
    .package = "ctwas")
  ctwasPipeline(gwasSumStats = ss, twasWeights = tw, twasZ = twasZ,
                regionId = "block1")
  expect_equal(capturedArgs$z_gene$id, "Q1|c1|t1|susie")
})
