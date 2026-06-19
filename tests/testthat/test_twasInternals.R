context("twasWeights internal helpers")

# ===========================================================================
# .twasNormalizeMethods
# ===========================================================================

test_that(".twasNormalizeMethods: NULL falls through to the 'default' preset", {
  res <- pecotmr:::.twasNormalizeMethods(NULL)
  default_names <- names(pecotmr:::.twasMethodLookup("default"))
  expect_equal(sort(names(res$methodList)), sort(default_names))
})

test_that(".twasNormalizeMethods: character preset string forwards to .twasMethodLookup", {
  res <- pecotmr:::.twasNormalizeMethods("fast_default")
  fast_names <- names(pecotmr:::.twasMethodLookup("fast_default"))
  expect_equal(sort(names(res$methodList)), sort(fast_names))
})

test_that(".twasNormalizeMethods: character vector of short names forwards to .twasMethodLookup", {
  res <- pecotmr:::.twasNormalizeMethods(c("lasso", "enet"))
  expect_equal(sort(names(res$methodList)),
               sort(c("lasso_weights", "enet_weights")))
})

test_that(".twasNormalizeMethods: named list passes through unchanged", {
  ml <- list(lassoWeights = list(), enetWeights = list(alpha = 0.5))
  res <- pecotmr:::.twasNormalizeMethods(ml)
  expect_identical(res$methodList, ml)
})

test_that(".twasNormalizeMethods: tokens strip both _weights and Weights suffixes", {
  res_snake <- pecotmr:::.twasNormalizeMethods(list(lasso_weights = list(),
                                                    enet_weights  = list()))
  res_camel <- pecotmr:::.twasNormalizeMethods(list(lassoWeights = list(),
                                                    enetWeights  = list()))
  expect_equal(res_snake$tokens, c("lasso", "enet"))
  expect_equal(res_camel$tokens, c("lasso", "enet"))
})

test_that(".twasNormalizeMethods: unrecognised input type errors", {
  expect_error(
    pecotmr:::.twasNormalizeMethods(42L),
    "must be a character vector, preset string, or named list"
  )
})

# ===========================================================================
# .twasCheckMethodCapabilities
# ===========================================================================

test_that(".twasCheckMethodCapabilities: empty token list is a no-op", {
  expect_silent(pecotmr:::.twasCheckMethodCapabilities(character(0), "QtlDataset"))
})

test_that(".twasCheckMethodCapabilities: individual-only token + QtlDataset is fine", {
  expect_silent(pecotmr:::.twasCheckMethodCapabilities(c("lasso", "enet"),
                                                       "QtlDataset"))
  expect_silent(pecotmr:::.twasCheckMethodCapabilities(c("lasso", "enet"),
                                                       "MultiTaskQtlDataset"))
})

test_that(".twasCheckMethodCapabilities: sumstat-only token + QtlDataset errors", {
  # prsCs has sumstatImpl only -- not legal for QtlDataset.
  expect_error(
    pecotmr:::.twasCheckMethodCapabilities("prsCs", "QtlDataset"),
    "is sumstat-only"
  )
})

test_that(".twasCheckMethodCapabilities: individual-only token + QtlSumStats errors", {
  # enet has individualImpl only -- not legal for QtlSumStats.
  expect_error(
    pecotmr:::.twasCheckMethodCapabilities("enet", "QtlSumStats"),
    "is individual-only"
  )
})

test_that(".twasCheckMethodCapabilities: unknown token errors with the full menu", {
  expect_error(
    pecotmr:::.twasCheckMethodCapabilities("bogus", "QtlDataset"),
    "unknown method token"
  )
})

# ===========================================================================
# .twasCheckMultivariateY
# ===========================================================================

test_that(".twasCheckMultivariateY: no multivariate tokens is a no-op", {
  expect_silent(pecotmr:::.twasCheckMultivariateY(c("lasso", "enet"),
                                                   nTraits = 1L, nContexts = 1L))
})

test_that(".twasCheckMultivariateY: mvsusie passes when contexts >= 2", {
  expect_silent(pecotmr:::.twasCheckMultivariateY("mvsusie",
                                                   nTraits = 1L, nContexts = 2L))
})

test_that(".twasCheckMultivariateY: mvsusie passes when traits >= 2", {
  expect_silent(pecotmr:::.twasCheckMultivariateY("mvsusie",
                                                   nTraits = 2L, nContexts = 1L))
})

test_that(".twasCheckMultivariateY: mvsusie/mrmash error with single trait, single context", {
  expect_error(
    pecotmr:::.twasCheckMultivariateY(c("mvsusie", "mrmash"),
                                       nTraits = 1L, nContexts = 1L),
    "require multi-trait or multi-context input"
  )
})

# ===========================================================================
# .twasAssertQcd
# ===========================================================================

.tw_makeSumStatsBare <- function() {
  gr <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(start = c(100L, 200L), width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = c("rs1", "rs2"), A1 = c("A", "A"), A2 = c("G", "G"),
    Z = c(1, 2), N = c(100L, 100L))
  gh <- new("GenotypeHandle",
            path = "/tmp/x.gds", format = "gds",
            snpInfo = data.frame(), nSamples = 0L,
            sampleIds = character(), pgenPtr = NULL)
  GwasSumStats(study = "g1", entry = list(gr),
                genome = "hg19", ldSketch = gh)
}

test_that(".twasAssertQcd: errors when qcInfo is empty", {
  ss <- .tw_makeSumStatsBare()
  expect_error(
    pecotmr:::.twasAssertQcd(ss),
    "has no QC record"
  )
})

test_that(".twasAssertQcd: passes when qcInfo is populated", {
  gr <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(start = c(100L, 200L), width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = c("rs1", "rs2"), A1 = c("A", "A"), A2 = c("G", "G"),
    Z = c(1, 2), N = c(100L, 100L))
  gh <- new("GenotypeHandle",
            path = "/tmp/x.gds", format = "gds",
            snpInfo = data.frame(), nSamples = 0L,
            sampleIds = character(), pgenPtr = NULL)
  ss <- GwasSumStats(study = "g1", entry = list(gr),
                      genome = "hg19", ldSketch = gh,
                      qcInfo = list(step1 = "ok"))
  expect_silent(pecotmr:::.twasAssertQcd(ss))
})

# ===========================================================================
# .twasSumstatsEntryToDf
# ===========================================================================

test_that(".twasSumstatsEntryToDf: returns canonical column layout", {
  gr <- GenomicRanges::GRanges(
    c("chr1", "chr1"),
    IRanges::IRanges(start = c(100L, 200L), width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = c("rs1", "rs2"), A1 = c("A", "A"), A2 = c("G", "G"),
    Z = c(1.0, -2.5), N = c(1000L, 1500L), MAF = c(0.1, 0.3))
  df <- pecotmr:::.twasSumstatsEntryToDf(gr)
  expect_s3_class(df, "data.frame")
  expect_equal(df$variant_id, c("rs1", "rs2"))
  expect_equal(df$chrom, c("chr1", "chr1"))
  expect_equal(df$pos, c(100L, 200L))
  expect_equal(df$z, c(1.0, -2.5))
  expect_equal(df$N, c(1000, 1500))
  expect_equal(df$maf, c(0.1, 0.3))
})

test_that(".twasSumstatsEntryToDf: derives z from beta/se when z is absent", {
  gr <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(start = 100L, width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = "rs1", A1 = "A", A2 = "G",
    BETA = 0.5, SE = 0.1, N = 1000L)
  df <- pecotmr:::.twasSumstatsEntryToDf(gr)
  expect_equal(df$beta, 0.5)
  expect_equal(df$se, 0.1)
  expect_equal(df$z, 0.5 / 0.1)
})

test_that(".twasSumstatsEntryToDf: omits optional columns when absent", {
  gr <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(start = 100L, width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = "rs1", A1 = "A", A2 = "G")
  df <- pecotmr:::.twasSumstatsEntryToDf(gr)
  expect_false("z"    %in% names(df))
  expect_false("beta" %in% names(df))
  expect_false("se"   %in% names(df))
  expect_false("N"    %in% names(df))
  expect_false("maf"  %in% names(df))
})

# ===========================================================================
# estimateSparsity
# ===========================================================================

test_that("estimateSparsity: legacy list input reads attr(.,'fit')$pi", {
  # Build a weight result that mimics learnTwasWeights(retainFits = TRUE):
  # element name carries the `_weights` suffix; attr 'fit' carries an mr.ash
  # object whose pi[1] is the spike weight.
  fake_w <- structure(c(0.1, 0, 0.3),
                      fit = list(pi = c(0.6, 0.2, 0.2)))
  weightResults <- list(mrash_weights = fake_w)
  expect_equal(estimateSparsity(weightResults), 1 - 0.6,
               tolerance = 1e-12)
})

test_that("estimateSparsity: TwasWeights collection input reads from the mrash entry", {
  entry <- TwasWeightsEntry(
    variantIds    = c("v1", "v2", "v3"),
    weights       = c(0.1, 0, 0.3),
    fits          = list(pi = c(0.7, 0.2, 0.1)),
    standardized  = FALSE)
  tw <- TwasWeights(
    study = "s1", context = "c1", trait = "t1", method = "mrash",
    entry = list(entry))
  expect_equal(estimateSparsity(tw), 1 - 0.7, tolerance = 1e-12)
})

test_that("estimateSparsity: TwasWeights without mrash entry errors", {
  entry <- TwasWeightsEntry(variantIds = c("v1"), weights = c(0.1))
  tw <- TwasWeights(
    study = "s1", context = "c1", trait = "t1", method = "lasso",
    entry = list(entry))
  expect_error(
    estimateSparsity(tw),
    "mr.ash entry not found in TwasWeights"
  )
})

test_that("estimateSparsity: TwasWeights with mrash entry but no fit$pi errors", {
  entry <- TwasWeightsEntry(variantIds = c("v1"), weights = c(0.1),
                            fits = list(other = 1))
  tw <- TwasWeights(
    study = "s1", context = "c1", trait = "t1", method = "mrash",
    entry = list(entry))
  expect_error(
    estimateSparsity(tw),
    "mr.ash fit object not found"
  )
})

test_that("estimateSparsity: legacy list input without mrash_weights errors", {
  expect_error(
    estimateSparsity(list(lasso_weights = c(0.1, 0.2))),
    "'mrash_weights'.*not found"
  )
})

test_that("estimateSparsity: legacy list input without fit attr errors", {
  expect_error(
    estimateSparsity(list(mrash_weights = c(0.1, 0.2))),
    "mr.ash fit object not found"
  )
})
