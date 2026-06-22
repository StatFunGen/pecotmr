context("ctwasPipeline")

# ===========================================================================
# Strategy: ctwas::ctwas_sumstats does the heavy work. We mock it to a
# function that just returns its inputs back, so we can verify how the
# pipeline assembles z_snp / weights / region_info / LD loader inputs.
# ===========================================================================

.ctp_makeHandle <- function(snp_n = 6L, n_samples = 30L) {
  # Use a per-process tempfile so .ctwasLdPanelKey's file.exists check
  # succeeds against the fixture handle (real LD-sketch payloads exist
  # by construction; mock fixtures need an equivalent on-disk anchor).
  gdsPath <- file.path(tempdir(), "ctp_sketch.gds")
  if (!file.exists(gdsPath)) file.create(gdsPath)
  positions <- seq(100L, by = 100L, length.out = snp_n)
  # SNP IDs follow the canonical chr:pos:A2:A1 layout so allele
  # harmonization inside .ctwasBuildWeights / .ctwasHarmonizeWeights can
  # parse them via parseVariantId().
  snpIds <- sprintf("chr1:%d:G:A", positions)
  new("GenotypeHandle",
    path = gdsPath,
    format = "gds",
    snpInfo = data.frame(
      SNP = snpIds,
      CHR = rep("1", snp_n),
      BP  = positions,
      A1  = rep("A", snp_n),
      A2  = rep("G", snp_n),
      stringsAsFactors = FALSE),
    nSamples = n_samples,
    sampleIds = paste0("s", seq_len(n_samples)),
    pgenPtr = NULL)
}

# Canonical SNP IDs the fixtures use. .ctp_makeHandle() emits these in
# chr:pos:A2:A1 form; tests reference them by index via .ctp_snpId(i).
.ctp_snpId <- function(i) sprintf("chr1:%d:G:A", 100L * i)

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
    SNP = vapply(1:6, .ctp_snpId, character(1)),
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
    variantIds = vapply(1:5, .ctp_snpId, character(1)),
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

# Helper: minimal two-block input set for the multi-block API tests.
.ctp_makeMultiBlockInputs <- function(qc = TRUE) {
  ss <- .ctp_makeGwasSumstats(qc = qc)
  tw <- .ctp_makeTwasWeights()
  list(
    gwasSumStats = list(block1 = ss, block2 = ss),
    twasWeights  = list(block1 = tw, block2 = tw))
}

test_that("ctwasPipeline: rejects a bare (non-list) GwasSumStats", {
  skip_if_not_installed("ctwas")
  expect_error(
    ctwasPipeline(gwasSumStats = .ctp_makeGwasSumstats(),
                  twasWeights  = .ctp_makeTwasWeights()),
    "NAMED LIST of GwasSumStats"
  )
})

test_that("ctwasPipeline: rejects a single-block named list", {
  skip_if_not_installed("ctwas")
  expect_error(
    ctwasPipeline(gwasSumStats = list(block1 = .ctp_makeGwasSumstats()),
                  twasWeights  = list(block1 = .ctp_makeTwasWeights())),
    "at least two LD blocks"
  )
})

test_that("ctwasPipeline: rejects un-QCd GwasSumStats in any region", {
  skip_if_not_installed("ctwas")
  ss_qc   <- .ctp_makeGwasSumstats(qc = TRUE)
  ss_noqc <- .ctp_makeGwasSumstats(qc = FALSE)
  tw      <- .ctp_makeTwasWeights()
  expect_error(
    ctwasPipeline(gwasSumStats = list(block1 = ss_qc, block2 = ss_noqc),
                  twasWeights  = list(block1 = tw,    block2 = tw)),
    "has no QC record"
  )
})

test_that("ctwasPipeline: rejects twasWeights keys not present in gwasSumStats", {
  skip_if_not_installed("ctwas")
  ss <- .ctp_makeGwasSumstats()
  tw <- .ctp_makeTwasWeights()
  expect_error(
    ctwasPipeline(gwasSumStats = list(blockA = ss, blockB = ss),
                  twasWeights  = list(blockA = tw, blockC = tw)),
    "key.*not present in.*gwasSumStats"
  )
})

test_that("assembleCtwasInputs: allows twasWeights keys to be a subset of gwasSumStats", {
  ss <- .ctp_makeGwasSumstats()
  tw <- .ctp_makeTwasWeights()
  local_mocked_bindings(extractBlockGenotypes = .ctp_mockExtractor(),
                        .package = "pecotmr")
  # Two blocks supply zSnp; only block1 supplies TwasWeights.
  inputs <- assembleCtwasInputs(
    gwasSumStats = list(block1 = ss, block2 = ss),
    twasWeights  = list(block1 = tw))
  # Both blocks appear in region_info / snp_map (SNP-only block2 contributes
  # its zSnp), but the weights list only has block1-keyed entries.
  expect_setequal(inputs$region_info$region_id, c("block1", "block2"))
  expect_setequal(names(inputs$snp_map), c("block1", "block2"))
  expect_true(all(grepl("^block1\\|", names(inputs$weights))))
})

test_that("assembleCtwasInputs: filters TwasWeights against UNION of all blocks' GWAS variants", {
  # Build two blocks with NON-OVERLAPPING GWAS variants. Block 1 covers
  # v1..v3, block 2 covers v4..v6. The gene's weight spans v2..v5 — i.e.
  # crosses the block boundary. With a per-block filter the gene would
  # lose v4..v5 (block-2 variants); with a global-union filter all four
  # weight variants survive.
  mkBlockGss <- function(study, snpIds, qc = TRUE) {
    gr <- GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges = IRanges::IRanges(start = as.integer(gsub(".*:([0-9]+):.*", "\\1", snpIds)),
                                width = 1L))
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
      SNP = snpIds,
      A1  = rep("A", length(snpIds)), A2 = rep("G", length(snpIds)),
      Z   = rnorm(length(snpIds)), N = rep(1000L, length(snpIds)))
    GwasSumStats(study = study, entry = list(gr), genome = "hg19",
                 ldSketch = .ctp_makeHandle(),
                 qcInfo   = if (qc) list(step1 = "ok") else list())
  }
  ss1 <- mkBlockGss("G1", vapply(1:3, .ctp_snpId, character(1)))
  ss2 <- mkBlockGss("G2", vapply(4:6, .ctp_snpId, character(1)))
  # Cross-boundary weights: v2..v5 (4 variants spanning both blocks).
  crossEntry <- TwasWeightsEntry(
    variantIds = vapply(2:5, .ctp_snpId, character(1)),
    weights    = c(0.1, 0.2, 0.3, 0.4))
  tw <- TwasWeights(
    study = "G1", context = "c1", trait = "t1", method = "susie",
    entry = list(crossEntry), ldSketch = .ctp_makeHandle())
  local_mocked_bindings(extractBlockGenotypes = .ctp_mockExtractor(),
                        .package = "pecotmr")
  inputs <- assembleCtwasInputs(
    gwasSumStats = list(block1 = ss1, block2 = ss2),
    twasWeights  = list(block1 = tw))
  # All four cross-boundary variants should appear in the weights list
  # (would be only 2 with a per-block filter).
  wgt <- inputs$weights[[1L]]$wgt
  expect_equal(nrow(wgt), 4L)
  expect_setequal(rownames(wgt), vapply(2:5, .ctp_snpId, character(1)))
})

test_that("ctwasPipeline: rejects bare (non-list) TwasWeights", {
  skip_if_not_installed("ctwas")
  expect_error(
    ctwasPipeline(gwasSumStats = list(block1 = .ctp_makeGwasSumstats(),
                                       block2 = .ctp_makeGwasSumstats()),
                  twasWeights  = .ctp_makeTwasWeights()),
    "NAMED LIST of TwasWeights"
  )
})

test_that("ctwasPipeline: rejects non-GRanges twasZ", {
  skip_if_not_installed("ctwas")
  inp <- .ctp_makeMultiBlockInputs()
  expect_error(
    ctwasPipeline(gwasSumStats = inp$gwasSumStats,
                  twasWeights  = inp$twasWeights,
                  twasZ        = "not a GRanges"),
    "must be a GRanges"
  )
})

test_that("ctwasPipeline: rejects unknown groupPriorVarStructure value", {
  skip_if_not_installed("ctwas")
  inp <- .ctp_makeMultiBlockInputs()
  expect_error(
    ctwasPipeline(gwasSumStats = inp$gwasSumStats,
                  twasWeights  = inp$twasWeights,
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
  expect_setequal(df$id, vapply(1:6, .ctp_snpId, character(1)))
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

test_that(".ctwasSnpInfoForBlock: returns ctwas-required columns chrom/id/pos/alt/ref", {
  df <- pecotmr:::.ctwasSnpInfoForBlock(.ctp_makeHandle())
  # ctwas's read_snp_info_files asserts these exact column names
  expect_setequal(colnames(df),
                  c("chrom", "id", "pos", "alt", "ref"))
})

test_that(".ctwasLdPanelKey: returns the on-disk path for an existing GDS sketch", {
  handle <- .ctp_makeHandle()
  key <- pecotmr:::.ctwasLdPanelKey(handle)
  expect_true(file.exists(key))
  expect_equal(key, getPath(handle))
})

test_that(".ctwasLdPanelKey: errors when no candidate file exists", {
  ghost <- new("GenotypeHandle",
    path = file.path(tempdir(), "does_not_exist_for_test.pgen_stem"),
    format = "plink2",
    snpInfo = data.frame(SNP = "v1", CHR = "1", BP = 100L, A1 = "A",
                          A2 = "G", stringsAsFactors = FALSE),
    nSamples = 1L, sampleIds = "s1", pgenPtr = NULL)
  expect_error(pecotmr:::.ctwasLdPanelKey(ghost),
                "could not derive an existing LD-file token")
})

test_that(".ctwasResolveMethod: caller-supplied method wins when present", {
  e <- TwasWeightsEntry(variantIds = paste0("v", 1:3),
                         weights = c(0.1, 0.2, 0.3))
  tw <- TwasWeights(
    study = "Q1", context = "c1", trait = "t1", method = "mrash",
    entry = list(e), ldSketch = .ctp_makeHandle())
  expect_equal(pecotmr:::.ctwasResolveMethod(list(r1 = tw), "mrash"),
                "mrash")
})

test_that(".ctwasResolveMethod: caller-supplied unknown method errors", {
  e <- TwasWeightsEntry(variantIds = paste0("v", 1:3),
                         weights = c(0.1, 0.2, 0.3))
  tw <- TwasWeights(
    study = "Q1", context = "c1", trait = "t1", method = "mrash",
    entry = list(e), ldSketch = .ctp_makeHandle())
  expect_error(pecotmr:::.ctwasResolveMethod(list(r1 = tw), "bogus"),
                "not present in TwasWeights")
})

test_that(".ctwasResolveMethod: defaults to ensemble when present among multiple", {
  mkTw <- function(m) {
    e <- TwasWeightsEntry(variantIds = paste0("v", 1:3),
                           weights = c(0.1, 0.2, 0.3))
    TwasWeights(study = "Q1", context = "c1", trait = "t1", method = m,
                 entry = list(e), ldSketch = .ctp_makeHandle())
  }
  # Build a multi-method TwasWeights by stitching two methods together.
  tw <- TwasWeights(
    study   = c("Q1", "Q1"), context = c("c1", "c1"),
    trait   = c("t1", "t1"), method  = c("mrash", "ensemble"),
    entry   = list(
      TwasWeightsEntry(variantIds = paste0("v", 1:3), weights = c(0.1, 0.2, 0.3)),
      TwasWeightsEntry(variantIds = paste0("v", 1:3), weights = c(0.4, 0.5, 0.6))),
    ldSketch = .ctp_makeHandle())
  expect_equal(pecotmr:::.ctwasResolveMethod(list(r1 = tw)), "ensemble")
})

test_that(".ctwasResolveMethod: single method auto-picked when only one available", {
  e <- TwasWeightsEntry(variantIds = paste0("v", 1:3),
                         weights = c(0.1, 0.2, 0.3))
  tw <- TwasWeights(
    study = "Q1", context = "c1", trait = "t1", method = "mrash",
    entry = list(e), ldSketch = .ctp_makeHandle())
  expect_equal(pecotmr:::.ctwasResolveMethod(list(r1 = tw)), "mrash")
})

test_that(".ctwasResolveMethod: multi-method + no ensemble + no caller method errors", {
  tw <- TwasWeights(
    study   = c("Q1", "Q1"), context = c("c1", "c1"),
    trait   = c("t1", "t1"), method  = c("mrash", "susie"),
    entry   = list(
      TwasWeightsEntry(variantIds = paste0("v", 1:3), weights = c(0.1, 0.2, 0.3)),
      TwasWeightsEntry(variantIds = paste0("v", 1:3), weights = c(0.4, 0.5, 0.6))),
    ldSketch = .ctp_makeHandle())
  expect_error(pecotmr:::.ctwasResolveMethod(list(r1 = tw)),
                "Supply a `method` argument")
})

test_that(".ctwasFilterMethod: subsets rows to the requested method", {
  tw <- TwasWeights(
    study   = c("Q1", "Q1"), context = c("c1", "c1"),
    trait   = c("t1", "t1"), method  = c("mrash", "susie"),
    entry   = list(
      TwasWeightsEntry(variantIds = paste0("v", 1:3), weights = c(0.1, 0.2, 0.3)),
      TwasWeightsEntry(variantIds = paste0("v", 1:3), weights = c(0.4, 0.5, 0.6))),
    ldSketch = .ctp_makeHandle())
  twSub <- pecotmr:::.ctwasFilterMethod(tw, "susie")
  expect_equal(nrow(twSub), 1L)
  expect_equal(as.character(twSub$method), "susie")
})

# Build an ldPanel fixture (matches .ctwasComputeFullPanelLd's return
# shape) for the 6-SNP toy panel from .ctp_makeHandle().
.ctp_makeLdPanel <- function(snp_n = 6L) {
  h <- .ctp_makeHandle(snp_n = snp_n)
  snpInfo <- pecotmr:::.ctwasSnpInfoForBlock(h)
  R <- diag(1, snp_n)
  dimnames(R) <- list(snpInfo$id, snpInfo$id)
  # Unit dosage variance — sqrt(1) = 1, so the variance-scaling step in
  # .ctwasBuildWeights is a no-op for this fixture.
  variance <- setNames(rep(1, snp_n), snpInfo$id)
  list(R = R, snpInfo = snpInfo, variance = variance)
}

test_that(".ctwasBuildWeights: scales non-standardized weights by sqrt(variance)", {
  panel <- .ctp_makeLdPanel()
  # Replace the default unit variance with non-trivial values; raw
  # weights should be multiplied by sqrt(variance) before reaching the
  # final wgt matrix.
  panel$variance <- setNames(c(0.5, 1, 2, 4, 8, 16), panel$snpInfo$id)
  ids5  <- vapply(1:5, .ctp_snpId, character(1))
  rawW  <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  tw <- TwasWeights(
    study = "Q1", context = "c1", trait = "t1", method = "susie",
    entry = list(TwasWeightsEntry(variantIds = ids5,
                                    weights    = rawW)),
    ldSketch = .ctp_makeHandle())
  wl <- pecotmr:::.ctwasBuildWeights(tw, panel)
  expected <- unname(rawW * sqrt(panel$variance[ids5]))
  expect_equal(as.numeric(wl[[1L]]$wgt), expected, tolerance = 1e-12)
})

test_that(".ctwasBuildWeights: standardized weights bypass variance scaling", {
  panel <- .ctp_makeLdPanel()
  panel$variance <- setNames(c(0.5, 1, 2, 4, 8, 16), panel$snpInfo$id)
  ids5 <- vapply(1:5, .ctp_snpId, character(1))
  rawW <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  tw <- TwasWeights(
    study = "Q1", context = "c1", trait = "t1", method = "susie",
    entry = list(TwasWeightsEntry(variantIds = ids5,
                                    weights      = rawW,
                                    standardized = TRUE)),
    ldSketch = .ctp_makeHandle())
  wl <- pecotmr:::.ctwasBuildWeights(tw, panel)
  expect_equal(as.numeric(wl[[1L]]$wgt), rawW, tolerance = 1e-12)
})

# Build an LD-panel fixture with realistic chr:pos:A2:A1 variant IDs so
# `.ctwasHarmonizeWeights` (which calls parseVariantId + .matchRefPanel)
# can do real allele matching.
.ctp_makeAllelePanel <- function() {
  ids <- c("1:100:C:T", "1:200:G:A", "1:300:A:G", "1:400:T:C")
  snpInfo <- data.frame(
    chrom = 1L,
    id    = ids,
    pos   = c(100L, 200L, 300L, 400L),
    alt   = c("T", "A", "G", "C"),   # A1 (effect)
    ref   = c("C", "G", "A", "T"),   # A2 (other)
    stringsAsFactors = FALSE)
  R <- diag(1, 4); dimnames(R) <- list(ids, ids)
  list(R = R, snpInfo = snpInfo,
       variance = setNames(rep(1, 4), ids))
}

test_that(".ctwasHarmonizeWeights: sign-flips weights for swapped A1/A2", {
  panel <- .ctp_makeAllelePanel()
  refVariants <- data.frame(
    chrom = panel$snpInfo$chrom, pos = panel$snpInfo$pos,
    A2 = panel$snpInfo$ref, A1 = panel$snpInfo$alt,
    variant_id = panel$snpInfo$id, stringsAsFactors = FALSE)
  # Variant 1: alleles match panel ("1:100:C:T" — same A2/A1 ordering)
  # Variant 2: A1/A2 swapped vs panel ("1:200:A:G" — flips relative to "1:200:G:A")
  # Output variant_id values are rebuilt via formatVariantId(), which
  # always emits a `chr` prefix.
  res <- pecotmr:::.ctwasHarmonizeWeights(
    origVids = c("1:100:C:T", "1:200:A:G"),
    origW    = c(0.5, 0.3),
    refVariants = refVariants)
  expect_equal(nrow(res), 2L)
  # Variant 1 keeps its sign; variant 2 should be sign-flipped to -0.3.
  matches <- match(c("chr1:100:C:T", "chr1:200:G:A"), res$variant_id)
  expect_equal(res$w[matches[[1L]]],  0.5, tolerance = 1e-12)
  expect_equal(res$w[matches[[2L]]], -0.3, tolerance = 1e-12)
})

test_that(".ctwasHarmonizeWeights: drops variants not present in the panel", {
  panel <- .ctp_makeAllelePanel()
  refVariants <- data.frame(
    chrom = panel$snpInfo$chrom, pos = panel$snpInfo$pos,
    A2 = panel$snpInfo$ref, A1 = panel$snpInfo$alt,
    variant_id = panel$snpInfo$id, stringsAsFactors = FALSE)
  res <- pecotmr:::.ctwasHarmonizeWeights(
    origVids = c("1:100:C:T", "1:999:A:T"),  # 1:999 not in panel
    origW    = c(0.5, 0.3),
    refVariants = refVariants)
  expect_equal(nrow(res), 1L)
  expect_equal(res$variant_id, "chr1:100:C:T")
})

test_that(".ctwasIsSusieFit: recognizes the susie intermediate shape", {
  fits <- list(lbf_variable = matrix(0, 2, 3),
                mu = matrix(0, 2, 3),
                X_column_scale_factors = rep(1, 3))
  expect_true(pecotmr:::.ctwasIsSusieFit(fits))
  expect_false(pecotmr:::.ctwasIsSusieFit(NULL))
  expect_false(pecotmr:::.ctwasIsSusieFit(list(lbf_variable = matrix(0, 2, 3))))
})

test_that(".ctwasRenormalizeSusieWeights: lbfToAlpha + colSums recomputation", {
  # Original fit covers 4 variants; drop variant 4, renormalize over {1,2,3}.
  origVids <- paste0("v", 1:4)
  origW    <- c(0.1, 0.2, 0.3, 0.4)
  # Toy lbf_variable: 2 effects, 4 variants. Constructed so the kept
  # subset yields easily-predictable softmax weights.
  lbf <- rbind(c(0, 0, 0, 100),    # effect 1: only v4 has signal
                c(10, 0, -10, 0))   # effect 2: v1 dominates
  mu  <- matrix(c(1, 2, 3, 4, 1, 2, 3, 4), nrow = 2, byrow = TRUE)
  xCol <- rep(1, 4)
  fits <- list(lbf_variable = lbf, mu = mu, X_column_scale_factors = xCol)
  out <- pecotmr:::.ctwasRenormalizeSusieWeights(
    fits,
    origVids = origVids,
    origW = origW,
    keptIdx = c(1L, 2L, 3L),
    harmonizedW = origW[1:3])
  # Effect 1 lbf over {v1,v2,v3} = c(0,0,0) -> uniform alpha = 1/3
  # Effect 2 lbf over {v1,v2,v3} = c(10,0,-10) -> v1 ≈ 1
  # weight[v1] = (1/3)*1 + ~1*1 = ~1.33
  expect_equal(length(out), 3L)
  expect_true(out[[1L]] > out[[2L]] && out[[1L]] > out[[3L]])
})

test_that(".ctwasRenormalizeSusieWeights: returns NULL on fit/entry dimension mismatch", {
  fits <- list(lbf_variable = matrix(0, 2, 5),
                mu = matrix(0, 2, 5),
                X_column_scale_factors = rep(1, 5))
  out <- pecotmr:::.ctwasRenormalizeSusieWeights(
    fits,
    origVids = paste0("v", 1:3),  # entry says 3, fit covers 5 -> mismatch
    origW = rep(0.1, 3),
    keptIdx = 1:3,
    harmonizedW = rep(0.1, 3))
  expect_null(out)
})

test_that(".ctwasRenormalizeSusieWeights: signFlip carries over to mu", {
  # All variants kept; harmonized weights have opposite sign on v2.
  origVids <- paste0("v", 1:3)
  origW    <- c(0.1, 0.2, 0.3)
  # Strong lbf concentrated on a single effect / single variant per row,
  # so the recomputed weight directly mirrors mu (alpha ≈ identity rows).
  lbf <- rbind(c(100, -100, -100),
                c(-100, 100, -100))
  mu  <- rbind(c(1, 2, 3),
                c(4, 5, 6))
  fits <- list(lbf_variable = lbf, mu = mu,
                X_column_scale_factors = rep(1, 3))
  harmW_noflip <- c(0.1, 0.2, 0.3)   # all positive
  harmW_v2flip <- c(0.1, -0.2, 0.3)  # v2 flipped
  outNoFlip <- pecotmr:::.ctwasRenormalizeSusieWeights(
    fits, origVids, origW, keptIdx = 1:3, harmonizedW = harmW_noflip)
  outV2flip <- pecotmr:::.ctwasRenormalizeSusieWeights(
    fits, origVids, origW, keptIdx = 1:3, harmonizedW = harmW_v2flip)
  # v1, v3 should be unchanged between the two; v2 should flip sign.
  expect_equal(outNoFlip[[1L]],  outV2flip[[1L]], tolerance = 1e-9)
  expect_equal(outNoFlip[[3L]],  outV2flip[[3L]], tolerance = 1e-9)
  expect_equal(outNoFlip[[2L]], -outV2flip[[2L]], tolerance = 1e-9)
})

test_that(".ctwasSnpInfoForGwasBlock: restricts panel snpInfo to block GWAS variants", {
  ss <- .ctp_makeGwasSumstats()
  panelInfo <- data.frame(
    chrom = 1L,
    id    = vapply(1:6, .ctp_snpId, character(1)),  # whole panel
    pos   = seq(100L, by = 100L, length.out = 6L),
    alt   = "A", ref = "G",
    stringsAsFactors = FALSE)
  blockInfo <- pecotmr:::.ctwasSnpInfoForGwasBlock(ss, panelInfo)
  # Restricted to the variant IDs on the GwasSumStats entry.
  expect_true(all(blockInfo$id %in%
                  as.character(S4Vectors::mcols(ss$entry[[1L]])$SNP)))
  expect_true(nrow(blockInfo) <= nrow(panelInfo))
})

test_that(".ctwasBuildWeights: keys per-tuple weights and stamps gene metadata", {
  tw <- .ctp_makeTwasWeights()
  panel <- .ctp_makeLdPanel()
  ids5  <- vapply(1:5, .ctp_snpId, character(1))
  wl <- pecotmr:::.ctwasBuildWeights(tw, panel)
  expect_equal(length(wl), 1L)
  expect_equal(names(wl), "Q1|c1|t1|susie")
  expect_equal(wl[[1L]]$study, "Q1")
  expect_equal(wl[[1L]]$context, "c1")
  expect_equal(wl[[1L]]$gene_name, "t1")
  # wgt is a variants x 1 matrix with rownames = SNP IDs
  expect_true(is.matrix(wl[[1L]]$wgt))
  expect_equal(dim(wl[[1L]]$wgt), c(5L, 1L))
  expect_equal(rownames(wl[[1L]]$wgt), ids5)
  # R_wgt is a 5x5 slice of the cached panel R
  expect_true(is.matrix(wl[[1L]]$R_wgt))
  expect_equal(dim(wl[[1L]]$R_wgt), c(5L, 5L))
  expect_equal(rownames(wl[[1L]]$R_wgt), ids5)
  expect_equal(wl[[1L]]$n_wgt, 5L)
  # And it is literally a slice of the panel R (no recompute path).
  expect_equal(wl[[1L]]$R_wgt, panel$R[ids5, ids5])
})

test_that(".ctwasBuildWeights: drops variants not present in the LD panel", {
  ids3   <- vapply(1:3, .ctp_snpId, character(1))
  missing <- c("chr1:99900:G:A", "chr1:99910:G:A")  # not in panel
  tw <- TwasWeights(
    study = "Q1", context = "c1", trait = "t1", method = "susie",
    entry = list(TwasWeightsEntry(
      variantIds = c(ids3, missing),
      weights    = c(0.1, 0.2, 0.3, 0.4, 0.5))),
    ldSketch = .ctp_makeHandle())
  panel <- .ctp_makeLdPanel()
  wl <- pecotmr:::.ctwasBuildWeights(tw, panel)
  expect_equal(nrow(wl[[1L]]$wgt), 3L)
  expect_equal(rownames(wl[[1L]]$wgt), ids3)
  expect_equal(wl[[1L]]$n_wgt, 3L)
})

test_that(".ctwasBuildWeights: intersects with gwasSnpIds when supplied", {
  # The LD panel covers ids 1..6; the per-block GWAS sumstats covers only
  # ids 1, 2, 4 (a subset). Weight variants that live in the panel but
  # outside the block (id 3 here) must be dropped, otherwise ctwas's
  # compute_gene_z asserts the weight variant is missing from z_snp.
  ids5     <- vapply(1:5, .ctp_snpId, character(1))
  blockIds <- vapply(c(1, 2, 4), .ctp_snpId, character(1))
  tw <- TwasWeights(
    study = "Q1", context = "c1", trait = "t1", method = "susie",
    entry = list(TwasWeightsEntry(
      variantIds = ids5,
      weights    = c(0.1, 0.2, 0.3, 0.4, 0.5))),
    ldSketch = .ctp_makeHandle())
  panel <- .ctp_makeLdPanel()
  wl <- pecotmr:::.ctwasBuildWeights(
    tw, panel, gwasSnpIds = blockIds)
  expect_equal(rownames(wl[[1L]]$wgt), blockIds)
  expect_equal(wl[[1L]]$n_wgt, 3L)
})

test_that(".ctwasComputeFullPanelLd: extracts once + returns cached R + snpInfo + variance", {
  local_mocked_bindings(extractBlockGenotypes = .ctp_mockExtractor(),
                        .package = "pecotmr")
  out <- pecotmr:::.ctwasComputeFullPanelLd(.ctp_makeHandle())
  ids6 <- vapply(1:6, .ctp_snpId, character(1))
  expect_named(out, c("R", "snpInfo", "variance"))
  expect_true(is.matrix(out$R))
  expect_equal(dim(out$R), c(6L, 6L))
  expect_equal(rownames(out$R), ids6)
  expect_setequal(colnames(out$snpInfo),
                   c("chrom", "id", "pos", "alt", "ref"))
  expect_named(out$variance, ids6)
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

test_that(".ctwasMultiBlockLdLoader: dispatches by LD_file token", {
  RA <- matrix(runif(4), 2, 2,
                dimnames = list(c("v1", "v2"), c("v1", "v2")))
  RB <- matrix(runif(4), 2, 2,
                dimnames = list(c("v3", "v4"), c("v3", "v4")))
  loader <- pecotmr:::.ctwasMultiBlockLdLoader(
    list(tokenA = list(R = RA), tokenB = list(R = RB)))
  expect_identical(loader("tokenA"), RA)
  expect_identical(loader("tokenB"), RB)
  expect_error(loader("unknown_token"), "no cached panel")
})

test_that(".ctwasMultiBlockSnpInfoLoader: dispatches by LD_file token", {
  infoA <- data.frame(chrom = 1L, id = c("v1", "v2"),
                       pos = c(100L, 200L),
                       alt = "A", ref = "G", stringsAsFactors = FALSE)
  infoB <- data.frame(chrom = 1L, id = c("v3", "v4"),
                       pos = c(300L, 400L),
                       alt = "C", ref = "T", stringsAsFactors = FALSE)
  loader <- pecotmr:::.ctwasMultiBlockSnpInfoLoader(
    list(tokenA = list(snpInfo = infoA),
         tokenB = list(snpInfo = infoB)))
  expect_identical(loader("tokenA"), infoA)
  expect_identical(loader("tokenB"), infoB)
  expect_error(loader("unknown_token"), "no cached panel")
})

# ===========================================================================
# Input-assembly shape checks via assembleCtwasInputs (no ctwas engine needed)
# ===========================================================================

test_that("assembleCtwasInputs: assembles the documented input shape for ctwas", {
  inp <- .ctp_makeMultiBlockInputs()
  local_mocked_bindings(extractBlockGenotypes = .ctp_mockExtractor(),
                        .package = "pecotmr")
  inputs <- assembleCtwasInputs(
    gwasSumStats = inp$gwasSumStats,
    twasWeights  = inp$twasWeights)
  # Two regions, two LD_map rows.
  expect_equal(nrow(inputs$region_info), 2L)
  expect_setequal(inputs$region_info$region_id, c("block1", "block2"))
  # zSnp is the concatenation of both blocks' Z columns.
  expect_equal(nrow(inputs$z_snp), 12L)
  # snp_map keyed by region_id.
  expect_setequal(names(inputs$snp_map), c("block1", "block2"))
  # Per-region weights keys are prefixed with the region_id.
  expect_true(all(grepl("^(block1|block2)\\|", names(inputs$weights))))
  # LD_map carries the same number of rows as regions.
  expect_equal(nrow(inputs$LD_map), 2L)
  # LD / snpInfo loader closures are present.
  expect_true(is.function(inputs$LD_loader_fun))
  expect_true(is.function(inputs$snpinfo_loader_fun))
})

test_that("assembleCtwasInputs: forwards a twasZ argument as z_gene", {
  inp <- .ctp_makeMultiBlockInputs()
  twasZ <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 200))
  S4Vectors::mcols(twasZ) <- S4Vectors::DataFrame(
    qtlStudy = "Q1", context = "c1", trait = "t1", method = "susie",
    twasZ = 1.5)
  local_mocked_bindings(extractBlockGenotypes = .ctp_mockExtractor(),
                        .package = "pecotmr")
  inputs <- assembleCtwasInputs(
    gwasSumStats = inp$gwasSumStats,
    twasWeights  = inp$twasWeights,
    twasZ        = twasZ)
  expect_equal(inputs$z_gene$id, "Q1|c1|t1|susie")
})

# ===========================================================================
# Step-wise dispatch: estCtwasParam → screenCtwasRegions → finemapCtwasRegions
# ===========================================================================

test_that("ctwasPipeline: dispatches assemble → est → screen → finemap and accumulates state", {
  skip_if_not_installed("ctwas")
  inp <- .ctp_makeMultiBlockInputs()
  capturedAssemble <- list(); capturedEst <- list()
  capturedScreen   <- list(); capturedFinemap <- list()
  local_mocked_bindings(
    assemble_region_data = function(...) {
      capturedAssemble <<- list(...)
      list(region_data    = list(block1 = list(stub = TRUE)),
           boundary_genes = list(stub = TRUE),
           z_gene         = data.frame(id = "t1", z = 1.5))
    },
    est_param = function(...) {
      capturedEst <<- list(...)
      list(group_prior = c(g = 0.1, SNP = 0.0001),
           group_prior_var = c(g = 5, SNP = 5))
    },
    screen_regions = function(...) {
      capturedScreen <<- list(...)
      list(screened_region_data = list(block1 = list(stub = TRUE)),
           screen_res_meta = "mocked")
    },
    finemap_regions = function(...) {
      capturedFinemap <<- list(...)
      list(finemap_res = data.frame(id = "t1"),
           susie_alpha_res = data.frame(susie_pip = 0.9))
    },
    .package = "ctwas")
  local_mocked_bindings(extractBlockGenotypes = .ctp_mockExtractor(),
                        .package = "pecotmr")
  out <- ctwasPipeline(
    gwasSumStats = inp$gwasSumStats,
    twasWeights  = inp$twasWeights)
  # Each ctwas step was invoked exactly once.
  expect_true(length(capturedAssemble) > 0L)
  expect_true(length(capturedEst) > 0L)
  expect_true(length(capturedScreen) > 0L)
  expect_true(length(capturedFinemap) > 0L)
  # est sees the assembled region_data.
  expect_named(capturedEst$region_data, "block1")
  # screen receives the param estimates as group_prior / group_prior_var.
  expect_equal(unname(capturedScreen$group_prior), c(0.1, 0.0001))
  # finemap consumes screen's screened_region_data.
  expect_named(capturedFinemap$region_data, "block1")
  # Output mirrors the ctwas_sumstats top-level shape.
  expect_setequal(
    names(out),
    c("z_gene", "param", "finemap_res", "susie_alpha_res",
      "region_data", "boundary_genes", "screen_res"))
})

test_that("estCtwasParam: fallbackToPrefit recovers from accurate-EM NaN divergence", {
  skip_if_not_installed("ctwas")
  inp <- .ctp_makeMultiBlockInputs()
  # Mock est_param to throw the documented NaN error, and fit_EM to
  # produce a stub prefit result. Verify estCtwasParam catches the
  # NaN error AND that the returned param is the prefit estimate.
  local_mocked_bindings(
    assemble_region_data = function(...) list(
      region_data    = list(block1 = list(stub = TRUE)),
      boundary_genes = list(),
      z_gene         = data.frame(id = "t1", z = 1.0)),
    compute_gene_z = function(...) data.frame(id = "t1", z = 1.0),
    est_param = function(...) stop("Estimated group_prior_var contains NAs!"),
    # ctwas:::fit_EM is internal — mock via the same `ctwas` namespace.
    fit_EM = function(region_data, ...) list(
      group_prior     = c(g = 0.05, SNP = 1e-4),
      group_prior_var = c(g = 4.0, SNP = 5.0),
      group_size      = c(g = 1, SNP = 100)),
    .package = "ctwas")
  local_mocked_bindings(extractBlockGenotypes = .ctp_mockExtractor(),
                        .package = "pecotmr")
  # Without fallback: the NaN error propagates.
  expect_error(
    estCtwasParam(assembleCtwasInputs(inp$gwasSumStats, inp$twasWeights),
                  fallbackToPrefit = FALSE),
    "contains NAs")
  # With fallback: prefit estimates are returned as the param.
  # .ctwasFitPrefitEm thin-scales the SNP group_prior (mirroring ctwas's
  # est_param), so the mocked SNP prior 1e-4 emerges as 1e-4 * thin
  # (default thin = 0.1) → 1e-5. The group_prior_var is not thinned.
  est <- estCtwasParam(
    assembleCtwasInputs(inp$gwasSumStats, inp$twasWeights),
    fallbackToPrefit = TRUE)
  expect_equal(unname(est$param$group_prior),     c(0.05, 1e-5))
  expect_equal(unname(est$param$group_prior_var), c(4.0, 5.0))
})

test_that("estCtwasParam / screenCtwasRegions / finemapCtwasRegions can be called independently", {
  skip_if_not_installed("ctwas")
  inp <- .ctp_makeMultiBlockInputs()
  local_mocked_bindings(
    assemble_region_data = function(...) list(
      region_data    = list(block1 = list(stub = TRUE)),
      boundary_genes = list(),
      z_gene         = data.frame(id = "t1", z = 1.0)),
    est_param = function(...) list(
      group_prior = c(g = 0.05, SNP = 1e-4),
      group_prior_var = c(g = 4, SNP = 5)),
    screen_regions = function(...) list(
      screened_region_data = list(block1 = list(stub = TRUE))),
    finemap_regions = function(...) list(
      finemap_res = data.frame(id = "t1"),
      susie_alpha_res = data.frame(pip = 0.5)),
    .package = "ctwas")
  local_mocked_bindings(extractBlockGenotypes = .ctp_mockExtractor(),
                        .package = "pecotmr")
  # Step 1
  inputs <- assembleCtwasInputs(inp$gwasSumStats, inp$twasWeights)
  expect_true("region_info" %in% names(inputs))
  expect_true("LD_loader_fun" %in% names(inputs))
  # Step 2
  est <- estCtwasParam(inputs)
  expect_true("region_data" %in% names(est))
  expect_true("param" %in% names(est))
  # User can OVERRIDE the estimated priors before screen/finemap — this is
  # the escape hatch for NaN-on-iter-2 EM divergence.
  est$param$group_prior     <- c(g = 0.2, SNP = 1e-4)
  est$param$group_prior_var <- c(g = 4.5, SNP = 5)
  # Step 3
  screened <- screenCtwasRegions(est)
  expect_true("screened_region_data" %in% names(screened))
  # Step 4
  final <- finemapCtwasRegions(screened)
  expect_setequal(
    names(final),
    c("z_gene", "param", "finemap_res", "susie_alpha_res",
      "region_data", "boundary_genes", "screen_res"))
})

# ===========================================================================
# Real-engine end-to-end: drives ctwas::ctwas_sumstats with the bundled
# example PLINK panel + synthetic TwasWeights. Exercises the LD-loader
# and snp-info-loader closure bodies as actually invoked by ctwas
# (mocked tests only construct the closures, never invoke them).
# ===========================================================================

test_that("ctwasPipeline: real-engine end-to-end on the bundled example panel", {
  skip_if_not_installed("ctwas")
  data(gwas_sumstats_s4_example)
  data(qtl_dataset_example)
  gss <- fixupExampleGenotypePaths(gwas_sumstats_s4_example)
  qd  <- fixupExampleGenotypePaths(qtl_dataset_example)
  gh  <- qd@genotypes

  # 5-variant synthetic gene from the bundled panel.
  vids <- gh@snpInfo$SNP[1:5]
  ent  <- TwasWeightsEntry(
    variantIds = vids,
    weights    = c(0.1, -0.2, 0.05, 0.0, 0.3))
  tw <- TwasWeights(
    study   = "study1", context = "brain",
    trait   = "ENSG_example", method = "susie",
    entry   = list(ent),
    ldSketch = gh)

  # ctwasPipeline now requires >= 2 blocks; replicate the same GWAS +
  # TWAS pair under two region_ids so the joint EM has something to
  # estimate against. The bundled toy panel is still too sparse for
  # the convergence checks, so we relax the gates.
  res <- suppressMessages(suppressWarnings(
    ctwasPipeline(
      gwasSumStats = list(blockA = gss, blockB = gss),
      twasWeights  = list(blockA = tw,  blockB = tw),
      niter        = 5L,
      niterPrefit  = 2L,
      # Toy panel: relax the production filters that gate out tiny inputs.
      min_group_size       = 1L,
      min_p_single_effect  = 0,
      filter_L             = FALSE)))

  # ctwas_sumstats returns these 7 elements on success.
  expect_named(res, c("z_gene", "param", "finemap_res", "susie_alpha_res",
                       "region_data", "boundary_genes", "screen_res"),
                ignore.order = TRUE)
  # The gene we passed in came through.
  expect_true(any(grepl("study1\\|brain\\|ENSG_example\\|susie", res$z_gene$id)))
  expect_true(all(c("susie_pip", "susie_alpha", "region_id")
                   %in% colnames(res$susie_alpha_res)))
})

# ===========================================================================
# .ctwasFilterVariants — ported from R/ctwasWrapper.R::trimCtwasVariants
# ===========================================================================
# The filter has four knobs:
#   1. twasWeightCutoff — drop |w| < cutoff
#   2. csMinCor         — high-purity CS rescue (must-keep)
#   3. minPipCutoff     — high-PIP rescue (must-keep)
#   4. maxNumVariants   — per-gene cap, prioritized by PIP then |w|

test_that(".ctwasFilterVariants: twasWeightCutoff drops low-magnitude variants", {
  vids <- paste0("v", 1:6)
  w    <- c(0.5, 0.001, 0.3, 0.0005, -0.4, 0)
  out <- pecotmr:::.ctwasFilterVariants(
    vids = vids, w = w, finemapAux = NULL,
    twasWeightCutoff = 0.01, csMinCor = 0.8,
    minPipCutoff = 0, maxNumVariants = Inf)
  # Survivors: v1 (0.5), v3 (0.3), v5 (-0.4) — three with |w| >= 0.01
  expect_setequal(out$vids, c("v1", "v3", "v5"))
})

test_that(".ctwasFilterVariants: maxNumVariants caps by |w| when no PIP", {
  vids <- paste0("v", 1:5)
  w    <- c(0.1, 0.5, 0.2, 0.4, 0.05)
  out <- pecotmr:::.ctwasFilterVariants(
    vids = vids, w = w, finemapAux = NULL,
    twasWeightCutoff = 0, csMinCor = 0.8,
    minPipCutoff = 0, maxNumVariants = 3)
  # Top 3 by |w|: v2 (0.5), v4 (0.4), v3 (0.2)
  expect_setequal(out$vids, c("v2", "v4", "v3"))
})

test_that(".ctwasFilterVariants: minPipCutoff rescues high-PIP variants from cap", {
  vids <- paste0("v", 1:5)
  w    <- c(0.5, 0.4, 0.3, 0.2, 0.1)
  finemapAux <- list(
    pip = setNames(c(0.01, 0.02, 0.8, 0.01, 0.95), vids),
    csMembers = list(),
    csPurity  = numeric(0))
  out <- pecotmr:::.ctwasFilterVariants(
    vids = vids, w = w, finemapAux = finemapAux,
    twasWeightCutoff = 0, csMinCor = 0.8,
    minPipCutoff = 0.5, maxNumVariants = 2)
  # Must-keep (PIP > 0.5): v3, v5. Cap is 2 → both kept.
  expect_setequal(out$vids, c("v3", "v5"))
})

test_that(".ctwasFilterVariants: csMinCor rescues high-purity CS members from cap", {
  vids <- paste0("v", 1:6)
  w    <- c(0.5, 0.4, 0.3, 0.2, 0.1, 0.05)
  finemapAux <- list(
    pip = setNames(rep(0, length(vids)), vids),
    csMembers = list(c("v3", "v6"), c("v2", "v4")),
    csPurity  = c(0.9, 0.5))  # CS 1 (v3, v6) is high-purity
  out <- pecotmr:::.ctwasFilterVariants(
    vids = vids, w = w, finemapAux = finemapAux,
    twasWeightCutoff = 0, csMinCor = 0.8,
    minPipCutoff = 0, maxNumVariants = 3)
  # Must-keep from high-purity CS: v3, v6. Remaining slot filled by
  # next-highest |w| that isn't must-keep: v1 (0.5).
  expect_setequal(out$vids, c("v3", "v6", "v1"))
})

test_that(".ctwasFilterVariants: returns NULL when no variants survive", {
  vids <- paste0("v", 1:3)
  w    <- c(0.001, 0.0005, 0.002)
  out <- pecotmr:::.ctwasFilterVariants(
    vids = vids, w = w, finemapAux = NULL,
    twasWeightCutoff = 0.5, csMinCor = 0.8,
    minPipCutoff = 0, maxNumVariants = Inf)
  expect_null(out)
})

test_that(".ctwasBuildWeights: maxNumVariants caps the per-gene weight matrix", {
  data(qtl_dataset_example)
  qd <- fixupExampleGenotypePaths(qtl_dataset_example)
  gh <- qd@genotypes
  vids <- gh@snpInfo$SNP[1:5]
  ent <- TwasWeightsEntry(
    variantIds = vids,
    weights    = c(0.1, -0.2, 0.05, 0.3, 0.15))
  tw <- TwasWeights(
    study = "study1", context = "brain",
    trait = "ENSG_example", method = "susie",
    entry = list(ent), ldSketch = gh)
  ldPanel <- pecotmr:::.ctwasComputeFullPanelLd(gh)
  wl <- pecotmr:::.ctwasBuildWeights(tw, ldPanel, maxNumVariants = 3L)
  expect_equal(wl[[1L]]$n_wgt, 3L)
  expect_equal(nrow(wl[[1L]]$wgt), 3L)
  # Top 3 by |w| from c(0.1, -0.2, 0.05, 0.3, 0.15): 0.3, -0.2, 0.15
  expect_setequal(rownames(wl[[1L]]$wgt), vids[c(4L, 2L, 5L)])
})

test_that(".ctwasBuildWeights: twasWeightCutoff drops low-magnitude variants", {
  data(qtl_dataset_example)
  qd <- fixupExampleGenotypePaths(qtl_dataset_example)
  gh <- qd@genotypes
  vids <- gh@snpInfo$SNP[1:5]
  ent <- TwasWeightsEntry(
    variantIds = vids,
    # v1 (0.005) and v3 (0.001) will be dropped at cutoff 0.01
    weights    = c(0.005, 0.2, 0.001, 0.3, 0.1))
  tw <- TwasWeights(
    study = "study1", context = "brain",
    trait = "ENSG_example", method = "susie",
    entry = list(ent), ldSketch = gh)
  ldPanel <- pecotmr:::.ctwasComputeFullPanelLd(gh)
  wl <- pecotmr:::.ctwasBuildWeights(tw, ldPanel, twasWeightCutoff = 0.01)
  expect_equal(wl[[1L]]$n_wgt, 3L)
  expect_setequal(rownames(wl[[1L]]$wgt), vids[c(2L, 4L, 5L)])
})
