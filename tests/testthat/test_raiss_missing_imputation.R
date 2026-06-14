context("RAISS missing-variant imputation in TWAS pipelines")

# Helper: build an LdData S4 object from a genotype matrix
.build_ld_data_from_X <- function(X) {
  R <- computeLd(X, method = "sample")
  ref <- parseVariantId(colnames(X))
  ref$variant_id <- colnames(X)
  gr <- GenomicRanges::GRanges(
    seqnames = ref$chrom,
    ranges = IRanges::IRanges(start = ref$pos, width = 1L)
  )
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    variant_id = ref$variant_id, A1 = ref$A1, A2 = ref$A2
  )
  bm <- data.frame(
    blockId = 1L, chrom = ref$chrom[1],
    blockStart = min(ref$pos), blockEnd = max(ref$pos),
    size = ncol(X), startIdx = 1L, endIdx = ncol(X)
  )
  LdData(correlation = R, variants = gr, blockMetadata = bm,
         nRef = nrow(X))
}

# Helper: drop a random subset of sumstats rows to create missing variants
.drop_random <- function(sumstats, frac_keep, seed) {
  set.seed(seed)
  keep <- sort(sample(nrow(sumstats), round(frac_keep * nrow(sumstats))))
  sumstats[keep, , drop = FALSE]
}

# =============================================================================
# Case 2: twasWeightsSumstatPipeline impute_missing
# =============================================================================

test_that("twasWeightsSumstatPipeline: imputeMissing = FALSE leaves sumstats as-is", {
  skip_if_not_installed("susieR")
  data(gwas_sumstats_example)
  data(eqtl_region_example)
  ld_data <- .build_ld_data_from_X(eqtl_region_example$X)
  ss_partial <- .drop_random(
    gwas_sumstats_example[, c("variant_id", "chrom", "pos", "A1", "A2", "z")],
    frac_keep = 0.7, seed = 42
  )
  result <- twasWeightsSumstatPipeline(
    sumstats = ss_partial, ldData = ld_data, n = 1000,
    methods = list(susieRss = list(max_iter = 30)),
    pThresholds = NULL, checkLdMethod = NULL,
    impute = FALSE, imputeMissing = FALSE, verbose = 0
  )
  expect_equal(length(getVariantIds(result$twasWeights)), nrow(ss_partial))
})

test_that("twasWeightsSumstatPipeline: imputeMissing = TRUE widens variant panel", {
  skip_if_not_installed("susieR")
  data(gwas_sumstats_example)
  data(eqtl_region_example)
  ld_data <- .build_ld_data_from_X(eqtl_region_example$X)
  ss_partial <- .drop_random(
    gwas_sumstats_example[, c("variant_id", "chrom", "pos", "A1", "A2", "z")],
    frac_keep = 0.7, seed = 42
  )

  result_on <- twasWeightsSumstatPipeline(
    sumstats = ss_partial, ldData = ld_data, n = 1000,
    methods = list(susieRss = list(max_iter = 30)),
    pThresholds = NULL, checkLdMethod = NULL,
    impute = FALSE, imputeMissing = TRUE, verbose = 0
  )

  n_kept <- length(getVariantIds(result_on$twasWeights))
  expect_gt(n_kept, nrow(ss_partial))
  expect_lte(n_kept, ncol(eqtl_region_example$X))
})

test_that("twasWeightsSumstatPipeline: low-R^2 imputations are dropped", {
  skip_if_not_installed("susieR")
  data(gwas_sumstats_example)
  data(eqtl_region_example)
  ld_data <- .build_ld_data_from_X(eqtl_region_example$X)
  ss_partial <- .drop_random(
    gwas_sumstats_example[, c("variant_id", "chrom", "pos", "A1", "A2", "z")],
    frac_keep = 0.7, seed = 42
  )

  # High R^2 threshold should yield fewer imputed variants than a low one
  result_strict <- twasWeightsSumstatPipeline(
    sumstats = ss_partial, ldData = ld_data, n = 1000,
    methods = list(susieRss = list(max_iter = 30)),
    pThresholds = NULL, checkLdMethod = NULL,
    impute = FALSE, imputeMissing = TRUE,
    imputeOpts = list(rcond = 0.01, r2Threshold = 0.95,
                       minimumLd = 5, lamb = 0.01),
    verbose = 0
  )
  result_lenient <- twasWeightsSumstatPipeline(
    sumstats = ss_partial, ldData = ld_data, n = 1000,
    methods = list(susieRss = list(max_iter = 30)),
    pThresholds = NULL, checkLdMethod = NULL,
    impute = FALSE, imputeMissing = TRUE,
    imputeOpts = list(rcond = 0.01, r2Threshold = 0.2,
                       minimumLd = 5, lamb = 0.01),
    verbose = 0
  )
  expect_lt(length(getVariantIds(result_strict$twasWeights)),
            length(getVariantIds(result_lenient$twasWeights)))
})

# =============================================================================
# Case 1: harmonizeTwas helper imputeMissingGwasForSketch
# =============================================================================
# Unit-test the helper directly, since harmonizeTwas() itself requires
# PLINK-style reference files on disk.

test_that("imputeMissingGwasForSketch: returns unchanged when nothing is missing", {
  data(gwas_sumstats_example)
  data(eqtl_region_example)
  X <- eqtl_region_example$X
  ref <- parseVariantId(colnames(X))
  ref$variant_id <- colnames(X)
  gwas <- gwas_sumstats_example[, c("variant_id", "chrom", "pos", "A1", "A2", "z")]

  result <- pecotmr:::imputeMissingGwasForSketch(
    gwasDataSumstats = gwas,
    sketchRefPanel = ref,
    sketchX = X,
    imputeOpts = list(rcond = 0.01, r2Threshold = 0.6,
                       minimumLd = 5, lamb = 0.01),
    contextLabel = "test-noop"
  )
  expect_equal(nrow(result), nrow(gwas))
})

test_that("imputeMissingGwasForSketch: widens GWAS to include sketch variants", {
  data(gwas_sumstats_example)
  data(eqtl_region_example)
  X <- eqtl_region_example$X
  ref <- parseVariantId(colnames(X))
  ref$variant_id <- colnames(X)
  gwas_full <- gwas_sumstats_example[, c("variant_id", "chrom", "pos", "A1", "A2", "z")]
  gwas_partial <- .drop_random(gwas_full, frac_keep = 0.6, seed = 7)

  result <- pecotmr:::imputeMissingGwasForSketch(
    gwasDataSumstats = gwas_partial,
    sketchRefPanel = ref,
    sketchX = X,
    imputeOpts = list(rcond = 0.01, r2Threshold = 0.6,
                       minimumLd = 5, lamb = 0.01),
    contextLabel = "test-widen"
  )
  expect_gt(nrow(result), nrow(gwas_partial))
  expect_lte(nrow(result), nrow(ref))
  # Imputed rows must have non-NA z and originally-missing IDs
  added_ids <- setdiff(result$variant_id, gwas_partial$variant_id)
  expect_gt(length(added_ids), 0)
  expect_true(all(!is.na(result$z[result$variant_id %in% added_ids])))
})

test_that("imputeMissingGwasForSketch: imputed beta = z and se = 1 when columns present", {
  data(gwas_sumstats_example)
  data(eqtl_region_example)
  X <- eqtl_region_example$X
  ref <- parseVariantId(colnames(X))
  ref$variant_id <- colnames(X)
  gwas <- gwas_sumstats_example[, c("variant_id", "chrom", "pos", "A1", "A2",
                                     "beta", "se", "z")]
  gwas_partial <- .drop_random(gwas, frac_keep = 0.6, seed = 9)

  result <- pecotmr:::imputeMissingGwasForSketch(
    gwasDataSumstats = gwas_partial,
    sketchRefPanel = ref,
    sketchX = X,
    imputeOpts = list(rcond = 0.01, r2Threshold = 0.6,
                       minimumLd = 5, lamb = 0.01),
    contextLabel = "test-beta-se"
  )
  added_ids <- setdiff(result$variant_id, gwas_partial$variant_id)
  if (length(added_ids) > 0) {
    added <- result[result$variant_id %in% added_ids, , drop = FALSE]
    expect_equal(added$beta, added$z)
    expect_true(all(added$se == 1))
  }
})

test_that("imputeMissingGwasForSketch: gracefully skips when required cols missing", {
  data(gwas_sumstats_example)
  data(eqtl_region_example)
  X <- eqtl_region_example$X
  ref <- parseVariantId(colnames(X))
  ref$variant_id <- colnames(X)
  # Drop required column and drop variants so the early-return shortcut is bypassed
  bad_gwas <- .drop_random(
    gwas_sumstats_example[, c("variant_id", "chrom", "pos", "A1", "A2")],
    frac_keep = 0.6, seed = 11
  )
  expect_warning(
    out <- pecotmr:::imputeMissingGwasForSketch(
      gwasDataSumstats = bad_gwas,
      sketchRefPanel = ref,
      sketchX = X,
      imputeOpts = list(rcond = 0.01, r2Threshold = 0.6,
                         minimumLd = 5, lamb = 0.01),
      contextLabel = "test-skip"
    ),
    "missing required columns"
  )
  expect_equal(nrow(out), nrow(bad_gwas))
})
