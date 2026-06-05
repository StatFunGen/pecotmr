context("RAISS missing-variant imputation in TWAS pipelines")

# Helper: build an LDData S4 object from a genotype matrix
.build_ld_data_from_X <- function(X) {
  R <- compute_LD(X, method = "sample")
  ref <- parse_variant_id(colnames(X))
  ref$variant_id <- colnames(X)
  gr <- GenomicRanges::GRanges(
    seqnames = ref$chrom,
    ranges = IRanges::IRanges(start = ref$pos, width = 1L)
  )
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    variant_id = ref$variant_id, A1 = ref$A1, A2 = ref$A2
  )
  bm <- data.frame(
    block_id = 1L, chrom = ref$chrom[1],
    block_start = min(ref$pos), block_end = max(ref$pos),
    size = ncol(X), start_idx = 1L, end_idx = ncol(X)
  )
  LDData(correlation = R, variants = gr, block_metadata = bm,
         n_ref = nrow(X))
}

# Helper: drop a random subset of sumstats rows to create missing variants
.drop_random <- function(sumstats, frac_keep, seed) {
  set.seed(seed)
  keep <- sort(sample(nrow(sumstats), round(frac_keep * nrow(sumstats))))
  sumstats[keep, , drop = FALSE]
}

# =============================================================================
# Case 2: twas_weights_sumstat_pipeline impute_missing
# =============================================================================

test_that("twas_weights_sumstat_pipeline: impute_missing = FALSE leaves sumstats as-is", {
  skip_if_not_installed("susieR")
  data(gwas_sumstats_example)
  data(eqtl_region_example)
  ld_data <- .build_ld_data_from_X(eqtl_region_example$X)
  ss_partial <- .drop_random(
    gwas_sumstats_example[, c("variant_id", "chrom", "pos", "A1", "A2", "z")],
    frac_keep = 0.7, seed = 42
  )
  result <- twas_weights_sumstat_pipeline(
    sumstats = ss_partial, LD_data = ld_data, n = 1000,
    methods = list(susie_rss = list(max_iter = 30)),
    p_thresholds = NULL, check_ld_method = NULL,
    impute = FALSE, impute_missing = FALSE, verbose = 0
  )
  expect_equal(length(getVariantIds(result$twas_weights)), nrow(ss_partial))
})

test_that("twas_weights_sumstat_pipeline: impute_missing = TRUE widens variant panel", {
  skip_if_not_installed("susieR")
  data(gwas_sumstats_example)
  data(eqtl_region_example)
  ld_data <- .build_ld_data_from_X(eqtl_region_example$X)
  ss_partial <- .drop_random(
    gwas_sumstats_example[, c("variant_id", "chrom", "pos", "A1", "A2", "z")],
    frac_keep = 0.7, seed = 42
  )

  result_on <- twas_weights_sumstat_pipeline(
    sumstats = ss_partial, LD_data = ld_data, n = 1000,
    methods = list(susie_rss = list(max_iter = 30)),
    p_thresholds = NULL, check_ld_method = NULL,
    impute = FALSE, impute_missing = TRUE, verbose = 0
  )

  n_kept <- length(getVariantIds(result_on$twas_weights))
  expect_gt(n_kept, nrow(ss_partial))
  expect_lte(n_kept, ncol(eqtl_region_example$X))
})

test_that("twas_weights_sumstat_pipeline: low-R^2 imputations are dropped", {
  skip_if_not_installed("susieR")
  data(gwas_sumstats_example)
  data(eqtl_region_example)
  ld_data <- .build_ld_data_from_X(eqtl_region_example$X)
  ss_partial <- .drop_random(
    gwas_sumstats_example[, c("variant_id", "chrom", "pos", "A1", "A2", "z")],
    frac_keep = 0.7, seed = 42
  )

  # High R^2 threshold should yield fewer imputed variants than a low one
  result_strict <- twas_weights_sumstat_pipeline(
    sumstats = ss_partial, LD_data = ld_data, n = 1000,
    methods = list(susie_rss = list(max_iter = 30)),
    p_thresholds = NULL, check_ld_method = NULL,
    impute = FALSE, impute_missing = TRUE,
    impute_opts = list(rcond = 0.01, R2_threshold = 0.95,
                       minimum_ld = 5, lamb = 0.01),
    verbose = 0
  )
  result_lenient <- twas_weights_sumstat_pipeline(
    sumstats = ss_partial, LD_data = ld_data, n = 1000,
    methods = list(susie_rss = list(max_iter = 30)),
    p_thresholds = NULL, check_ld_method = NULL,
    impute = FALSE, impute_missing = TRUE,
    impute_opts = list(rcond = 0.01, R2_threshold = 0.2,
                       minimum_ld = 5, lamb = 0.01),
    verbose = 0
  )
  expect_lt(length(getVariantIds(result_strict$twas_weights)),
            length(getVariantIds(result_lenient$twas_weights)))
})

# =============================================================================
# Case 1: harmonize_twas helper impute_missing_gwas_for_sketch
# =============================================================================
# Unit-test the helper directly, since harmonize_twas() itself requires
# PLINK-style reference files on disk.

test_that("impute_missing_gwas_for_sketch: returns unchanged when nothing is missing", {
  data(gwas_sumstats_example)
  data(eqtl_region_example)
  X <- eqtl_region_example$X
  ref <- parse_variant_id(colnames(X))
  ref$variant_id <- colnames(X)
  gwas <- gwas_sumstats_example[, c("variant_id", "chrom", "pos", "A1", "A2", "z")]

  result <- pecotmr:::impute_missing_gwas_for_sketch(
    gwas_data_sumstats = gwas,
    sketch_ref_panel = ref,
    sketch_X = X,
    impute_opts = list(rcond = 0.01, R2_threshold = 0.6,
                       minimum_ld = 5, lamb = 0.01),
    context_label = "test-noop"
  )
  expect_equal(nrow(result), nrow(gwas))
})

test_that("impute_missing_gwas_for_sketch: widens GWAS to include sketch variants", {
  data(gwas_sumstats_example)
  data(eqtl_region_example)
  X <- eqtl_region_example$X
  ref <- parse_variant_id(colnames(X))
  ref$variant_id <- colnames(X)
  gwas_full <- gwas_sumstats_example[, c("variant_id", "chrom", "pos", "A1", "A2", "z")]
  gwas_partial <- .drop_random(gwas_full, frac_keep = 0.6, seed = 7)

  result <- pecotmr:::impute_missing_gwas_for_sketch(
    gwas_data_sumstats = gwas_partial,
    sketch_ref_panel = ref,
    sketch_X = X,
    impute_opts = list(rcond = 0.01, R2_threshold = 0.6,
                       minimum_ld = 5, lamb = 0.01),
    context_label = "test-widen"
  )
  expect_gt(nrow(result), nrow(gwas_partial))
  expect_lte(nrow(result), nrow(ref))
  # Imputed rows must have non-NA z and originally-missing IDs
  added_ids <- setdiff(result$variant_id, gwas_partial$variant_id)
  expect_gt(length(added_ids), 0)
  expect_true(all(!is.na(result$z[result$variant_id %in% added_ids])))
})

test_that("impute_missing_gwas_for_sketch: imputed beta = z and se = 1 when columns present", {
  data(gwas_sumstats_example)
  data(eqtl_region_example)
  X <- eqtl_region_example$X
  ref <- parse_variant_id(colnames(X))
  ref$variant_id <- colnames(X)
  gwas <- gwas_sumstats_example[, c("variant_id", "chrom", "pos", "A1", "A2",
                                     "beta", "se", "z")]
  gwas_partial <- .drop_random(gwas, frac_keep = 0.6, seed = 9)

  result <- pecotmr:::impute_missing_gwas_for_sketch(
    gwas_data_sumstats = gwas_partial,
    sketch_ref_panel = ref,
    sketch_X = X,
    impute_opts = list(rcond = 0.01, R2_threshold = 0.6,
                       minimum_ld = 5, lamb = 0.01),
    context_label = "test-beta-se"
  )
  added_ids <- setdiff(result$variant_id, gwas_partial$variant_id)
  if (length(added_ids) > 0) {
    added <- result[result$variant_id %in% added_ids, , drop = FALSE]
    expect_equal(added$beta, added$z)
    expect_true(all(added$se == 1))
  }
})

test_that("impute_missing_gwas_for_sketch: gracefully skips when required cols missing", {
  data(gwas_sumstats_example)
  data(eqtl_region_example)
  X <- eqtl_region_example$X
  ref <- parse_variant_id(colnames(X))
  ref$variant_id <- colnames(X)
  # Drop required column and drop variants so the early-return shortcut is bypassed
  bad_gwas <- .drop_random(
    gwas_sumstats_example[, c("variant_id", "chrom", "pos", "A1", "A2")],
    frac_keep = 0.6, seed = 11
  )
  expect_warning(
    out <- pecotmr:::impute_missing_gwas_for_sketch(
      gwas_data_sumstats = bad_gwas,
      sketch_ref_panel = ref,
      sketch_X = X,
      impute_opts = list(rcond = 0.01, R2_threshold = 0.6,
                         minimum_ld = 5, lamb = 0.01),
      context_label = "test-skip"
    ),
    "missing required columns"
  )
  expect_equal(nrow(out), nrow(bad_gwas))
})
