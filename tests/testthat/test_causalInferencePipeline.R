context("causalInferencePipeline")

# ===========================================================================
# Strategy: mock extractBlockGenotypes so .cipLdFromSketch returns a real
# LD matrix on a small panel. Everything else (twasZ, MR, p-value combine)
# runs for real on the tiny fixture.
# ===========================================================================

.cip_makeHandle <- function(snp_n = 6L, n_samples = 30L,
                            sample_prefix = "s") {
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
    sampleIds = paste0(sample_prefix, seq_len(n_samples)),
    pgenPtr = NULL)
}

.cip_mockExtractor <- function(seed = 7, n_samples = 30L) {
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

.cip_makeGwasSumstats <- function(qc = TRUE) {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = seq(100L, by = 100L, length.out = 5L),
                              width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = paste0("v", 1:5),
    A1  = rep("A", 5), A2  = rep("G", 5),
    Z   = c(2.0, -1.5, 0.5, 1.2, -0.8),
    N   = rep(1000L, 5),
    MAF = rep(0.3, 5))
  GwasSumStats(
    study    = "G1",
    entry    = list(gr),
    genome   = "hg19",
    ldSketch = .cip_makeHandle(),
    qcInfo   = if (qc) list(step1 = "ok") else list())
}

.cip_makeTwasWeights <- function(method = "susie",
                                  variant_ids = paste0("v", 1:5)) {
  entry <- TwasWeightsEntry(
    variantIds = variant_ids,
    weights    = c(0.1, 0.05, -0.2, 0.3, 0.0))
  TwasWeights(
    study    = "Q1", context = "c1", trait = "t1", method = method,
    entry    = list(entry),
    ldSketch = .cip_makeHandle())
}

.cip_makeQtlFmr <- function(variant_ids = paste0("v", 1:5)) {
  tl <- data.frame(
    variant_id = variant_ids,
    pip        = c(0.9, 0.05, 0.5, 0.8, 0.01),
    betahat    = c(0.2, 0.05, -0.1, 0.3, 0.0),
    sebetahat  = rep(0.05, length(variant_ids)),
    stringsAsFactors = FALSE)
  e <- FineMappingEntry(variantIds = variant_ids,
                        trimmedFit = list(),
                        topLoci    = tl)
  QtlFineMappingResult(
    study   = "Q1", context = "c1", trait = "t1", method = "susie",
    entry   = list(e),
    ldSketch = .cip_makeHandle())
}

# ===========================================================================
# Input-type validation
# ===========================================================================

test_that("causalInferencePipeline: rejects non-GwasSumStats input", {
  expect_error(causalInferencePipeline(gwasSumStats = "no"),
               "must be a GwasSumStats")
})

test_that("causalInferencePipeline: rejects un-QCd GwasSumStats", {
  expect_error(
    causalInferencePipeline(gwasSumStats = .cip_makeGwasSumstats(qc = FALSE),
                            twasWeights  = .cip_makeTwasWeights()),
    "has no QC record"
  )
})

test_that("causalInferencePipeline: requires at least one of twasWeights/fineMappingResult", {
  expect_error(
    causalInferencePipeline(gwasSumStats = .cip_makeGwasSumstats()),
    "at least one of"
  )
})

test_that("causalInferencePipeline: rejects non-TwasWeights twasWeights arg", {
  expect_error(
    causalInferencePipeline(gwasSumStats = .cip_makeGwasSumstats(),
                            twasWeights  = "not a TwasWeights"),
    "must be a TwasWeights"
  )
})

test_that("causalInferencePipeline: rejects GwasFineMappingResult for the QTL slot", {
  e <- FineMappingEntry(variantIds = "v1", trimmedFit = list(),
                        topLoci = data.frame(variant_id = "v1", pip = 0.1,
                                              stringsAsFactors = FALSE))
  gfmr <- GwasFineMappingResult(
    study = "G1", method = "susie", entry = list(e))
  expect_error(
    causalInferencePipeline(gwasSumStats      = .cip_makeGwasSumstats(),
                            fineMappingResult = gfmr),
    "does not accept GWAS-side fine"
  )
})

# ===========================================================================
# .cipRequireMatchingLdSketches branches
# ===========================================================================

test_that(".cipRequireMatchingLdSketches: NULL twas-side ldSketch is allowed", {
  # Build a TwasWeights without an ldSketch.
  twNoLd <- TwasWeights(
    study = "Q1", context = "c1", trait = "t1", method = "lasso",
    entry = list(TwasWeightsEntry(variantIds = paste0("v", 1:5),
                                   weights = rep(0.1, 5))),
    ldSketch = NULL)
  local_mocked_bindings(extractBlockGenotypes = .cip_mockExtractor(),
                        .package = "pecotmr")
  out <- causalInferencePipeline(
    gwasSumStats = .cip_makeGwasSumstats(),
    twasWeights  = twNoLd)
  expect_s4_class(out, "GRanges")
})

test_that(".cipRequireMatchingLdSketches: panel size mismatch errors", {
  bigSketch <- .cip_makeHandle(snp_n = 7L)
  twBig <- TwasWeights(
    study = "Q1", context = "c1", trait = "t1", method = "lasso",
    entry = list(TwasWeightsEntry(variantIds = paste0("v", 1:5),
                                   weights = rep(0.1, 5))),
    ldSketch = bigSketch)
  expect_error(
    causalInferencePipeline(gwasSumStats = .cip_makeGwasSumstats(),
                            twasWeights  = twBig),
    "differ in size"
  )
})

# ===========================================================================
# Happy path: TwasWeights only
# ===========================================================================

test_that("causalInferencePipeline: returns GRanges with TWAS Z per tuple", {
  tw <- .cip_makeTwasWeights()
  local_mocked_bindings(extractBlockGenotypes = .cip_mockExtractor(),
                        .package = "pecotmr")
  out <- causalInferencePipeline(gwasSumStats = .cip_makeGwasSumstats(),
                                  twasWeights  = tw)
  expect_s4_class(out, "GRanges")
  expect_equal(length(out), 1L)
  mc <- S4Vectors::mcols(out)
  expect_equal(as.character(mc$qtlStudy[[1L]]), "Q1")
  expect_equal(as.character(mc$gwasStudy[[1L]]), "G1")
  expect_true(is.finite(mc$twasZ[[1L]]))
  expect_true(is.finite(mc$twasPval[[1L]]))
  # No FMR -> MR fields stay NA.
  expect_true(is.na(mc$waldRatio[[1L]]))
  expect_true(is.na(mc$mrPval[[1L]]))
})

# ===========================================================================
# Happy path: TwasWeights + matching FineMappingResult enables MR
# ===========================================================================

test_that("causalInferencePipeline: with both inputs, MR fields are populated", {
  tw <- .cip_makeTwasWeights()
  fmr <- .cip_makeQtlFmr()
  local_mocked_bindings(extractBlockGenotypes = .cip_mockExtractor(),
                        .package = "pecotmr")
  out <- causalInferencePipeline(
    gwasSumStats       = .cip_makeGwasSumstats(),
    twasWeights        = tw,
    fineMappingResult  = fmr,
    mrPipCutoff        = 0.5)
  mc <- S4Vectors::mcols(out)
  # MR uses PIP > 0.5 variants from the FMR (v1, v3, v4).
  expect_true(is.finite(mc$waldRatio[[1L]]))
  expect_true(is.finite(mc$mrPval[[1L]]))
  expect_gt(mc$nIV[[1L]], 0L)
})

# ===========================================================================
# FineMappingResult-only path: weights come from topLoci$betahat
# ===========================================================================

test_that("causalInferencePipeline: FMR-only path extracts weights from topLoci$betahat", {
  fmr <- .cip_makeQtlFmr()
  local_mocked_bindings(extractBlockGenotypes = .cip_mockExtractor(),
                        .package = "pecotmr")
  out <- causalInferencePipeline(
    gwasSumStats       = .cip_makeGwasSumstats(),
    fineMappingResult  = fmr)
  mc <- S4Vectors::mcols(out)
  expect_true(is.finite(mc$twasZ[[1L]]))
})

# ===========================================================================
# combineMethods integration
# ===========================================================================

test_that("causalInferencePipeline: combineMethods appends combined rows", {
  tw1 <- TwasWeights(
    study   = c("Q1", "Q1"), context = c("c1", "c1"),
    trait   = c("t1", "t1"), method = c("lasso", "enet"),
    entry = list(
      TwasWeightsEntry(variantIds = paste0("v", 1:5), weights = rep(0.1, 5)),
      TwasWeightsEntry(variantIds = paste0("v", 1:5), weights = rep(0.05, 5))),
    ldSketch = .cip_makeHandle())
  local_mocked_bindings(extractBlockGenotypes = .cip_mockExtractor(),
                        .package = "pecotmr")
  out <- causalInferencePipeline(
    gwasSumStats   = .cip_makeGwasSumstats(),
    twasWeights    = tw1,
    combineMethods = "acat")
  # 2 per-tuple rows + 1 combined row = 3.
  expect_equal(length(out), 3L)
  mc <- S4Vectors::mcols(out)
  expect_true(any(grepl("^combined\\.", as.character(mc$method))))
})

# ===========================================================================
# .cipZToBeta / .cipZToSe fallbacks
# ===========================================================================

test_that(".cipZToBeta: falls back to z when maf/n are NA", {
  res <- pecotmr:::.cipZToBeta(z = c(1, 2), maf = NA, n = NA)
  expect_equal(res, c(1, 2))
})

test_that(".cipZToSe: falls back to vector of 1 when maf/n are NA", {
  res <- pecotmr:::.cipZToSe(z = c(1, 2), maf = NA, n = NA)
  expect_equal(res, c(1, 1))
})
