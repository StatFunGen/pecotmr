# Tests for SVD-based TWAS computation (LD sketch path)

# Phase 1: twas_z() and twas_analysis() SVD branch

test_that("twas_z: SVD path matches R path for full-rank genotype matrix", {
  set.seed(42)
  n <- 100 # samples (sketch size)
  p <- 20 # variants

  # Generate genotype-like matrix (dosages 0/1/2)
  X <- matrix(rbinom(n * p, 2, runif(p, 0.1, 0.9)), nrow = n, ncol = p)

  # HWE-based standardization
  af <- colMeans(X) / 2
  X_std <- sweep(X, 2, 2 * af)
  X_std <- sweep(X_std, 2, sqrt(2 * af * (1 - af)), "/")

  # Compute R from standardized X
  R <- crossprod(X_std) / (n - 1)

  # SVD of standardized X
  svd_result <- svd(X_std)

  # Random weights and z-scores
  weights <- rnorm(p)
  z <- rnorm(p)

  # R path
  result_R <- pecotmr:::twas_z(weights, z, R = R)

  # SVD path
  result_SVD <- pecotmr:::twas_z(weights, z, V = svd_result$v, D = svd_result$d, n_sketch = n)

  expect_equal(as.numeric(result_SVD$z), as.numeric(result_R$z), tolerance = 1e-10)
  expect_equal(as.numeric(result_SVD$pval), as.numeric(result_R$pval), tolerance = 1e-10)
})

test_that("twas_z: SVD path matches R path for rank-deficient matrix (n < p)", {
  set.seed(123)
  n <- 15 # fewer samples than variants

  p <- 30

  X <- matrix(rbinom(n * p, 2, runif(p, 0.15, 0.85)), nrow = n, ncol = p)

  af <- colMeans(X) / 2
  # Remove monomorphic columns
  keep <- af > 0 & af < 1
  X <- X[, keep, drop = FALSE]
  af <- af[keep]
  p <- ncol(X)

  X_std <- sweep(X, 2, 2 * af)
  X_std <- sweep(X_std, 2, sqrt(2 * af * (1 - af)), "/")

  R <- crossprod(X_std) / (n - 1)
  svd_result <- svd(X_std)

  weights <- rnorm(p)
  z <- rnorm(p)

  result_R <- pecotmr:::twas_z(weights, z, R = R)
  result_SVD <- pecotmr:::twas_z(weights, z, V = svd_result$v, D = svd_result$d, n_sketch = n)

  expect_equal(as.numeric(result_SVD$z), as.numeric(result_R$z), tolerance = 1e-10)
  expect_equal(as.numeric(result_SVD$pval), as.numeric(result_R$pval), tolerance = 1e-10)
})

test_that("twas_z: error when weights and z have different lengths", {
  expect_error(
    pecotmr:::twas_z(rnorm(5), rnorm(3), V = matrix(1, 5, 2), D = c(1, 1), n_sketch = 10),
    "Weights and z-scores must have the same length"
  )
})

test_that("twas_analysis: SVD path produces same results as R path", {
  set.seed(99)
  n <- 50
  p <- 10
  variant_ids <- paste0("chr1:", seq(1000, by = 100, length.out = p), ":A:G")

  X <- matrix(rbinom(n * p, 2, runif(p, 0.2, 0.8)), nrow = n, ncol = p)
  af <- colMeans(X) / 2
  X_std <- sweep(X, 2, 2 * af)
  X_std <- sweep(X_std, 2, sqrt(2 * af * (1 - af)), "/")

  R <- crossprod(X_std) / (n - 1)
  rownames(R) <- colnames(R) <- variant_ids
  svd_result <- svd(X_std)

  # Weights matrix (2 methods)
  weights_matrix <- matrix(rnorm(p * 2), nrow = p, ncol = 2)
  rownames(weights_matrix) <- variant_ids
  colnames(weights_matrix) <- c("lasso_weights", "enet_weights")

  # GWAS data
  gwas_df <- data.frame(variant_id = variant_ids, z = rnorm(p))

  # Use a subset of variants
  extract_variants <- variant_ids[3:8]

  result_R <- pecotmr:::twas_analysis(weights_matrix, gwas_df, LD_matrix = R,
                                       extract_variants_objs = extract_variants)
  result_SVD <- pecotmr:::twas_analysis(weights_matrix, gwas_df,
                                         extract_variants_objs = extract_variants,
                                         V = svd_result$v, D = svd_result$d,
                                         n_sketch = n, ld_variant_ids = variant_ids)

  expect_equal(as.numeric(result_SVD[[1]]$z), as.numeric(result_R[[1]]$z), tolerance = 1e-10)
  expect_equal(as.numeric(result_SVD[[2]]$z), as.numeric(result_R[[2]]$z), tolerance = 1e-10)
  expect_equal(as.numeric(result_SVD[[1]]$pval), as.numeric(result_R[[1]]$pval), tolerance = 1e-10)
  expect_equal(as.numeric(result_SVD[[2]]$pval), as.numeric(result_R[[2]]$pval), tolerance = 1e-10)
})

test_that("twas_analysis: SVD path handles partial variant overlap", {
  set.seed(77)
  n <- 40
  p <- 8
  variant_ids <- paste0("chr1:", seq(1000, by = 100, length.out = p), ":A:G")

  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n, ncol = p)
  af <- colMeans(X) / 2
  X_std <- sweep(X, 2, 2 * af)
  X_std <- sweep(X_std, 2, sqrt(2 * af * (1 - af)), "/")
  svd_result <- svd(X_std)

  weights_matrix <- matrix(rnorm(p), nrow = p, ncol = 1)
  rownames(weights_matrix) <- variant_ids
  colnames(weights_matrix) <- "lasso_weights"

  gwas_df <- data.frame(variant_id = variant_ids, z = rnorm(p))

  # Request variants where some are NOT in ld_variant_ids
  extra_variant <- "chr1:5000:A:G"
  extract_variants <- c(variant_ids[1:4], extra_variant)

  result <- pecotmr:::twas_analysis(weights_matrix, gwas_df,
                                     extract_variants_objs = extract_variants,
                                     V = svd_result$v, D = svd_result$d,
                                     n_sketch = n, ld_variant_ids = variant_ids)

  # Should succeed using only the 4 valid variants
  expect_false(is.null(result))
  expect_equal(length(result[[1]]$z), 1)
})

test_that("twas_analysis: SVD path returns NULL when no variants overlap", {
  variant_ids <- paste0("chr1:", 1:5, ":A:G")
  other_ids <- paste0("chr2:", 1:5, ":A:G")

  weights_matrix <- matrix(1, nrow = 5, ncol = 1)
  rownames(weights_matrix) <- variant_ids
  colnames(weights_matrix) <- "w"
  gwas_df <- data.frame(variant_id = variant_ids, z = rnorm(5))

  result <- suppressWarnings(pecotmr:::twas_analysis(
    weights_matrix, gwas_df,
    extract_variants_objs = variant_ids,
    V = matrix(1, 5, 2), D = c(1, 1),
    n_sketch = 10, ld_variant_ids = other_ids
  ))

  expect_null(result)
})

# Phase 2: load_ld_sketch() and standardize_genotype_hwe()

test_that("standardize_genotype_hwe: centers by 2p and scales by sqrt(2p(1-p))", {
  set.seed(42)
  n <- 30
  p <- 5
  af <- runif(p, 0.1, 0.9)
  X <- matrix(rbinom(n * p, 2, rep(af, each = n)), nrow = n, ncol = p)

  X_std <- pecotmr:::standardize_genotype_hwe(X, af)

  # Manual verification
  expected <- sweep(sweep(X, 2, 2 * af), 2, sqrt(2 * af * (1 - af)), "/")
  expect_equal(X_std, expected, tolerance = 1e-14)
})

test_that("load_ld_sketch: returns LDData with raw genotypes and metadata", {
  set.seed(55)
  n <- 30
  p <- 12
  variant_ids <- paste0("chr1:", seq(1000, by = 100, length.out = p), ":A:G")

  # Create a mock genotype matrix
  af_true <- runif(p, 0.1, 0.9)
  X <- matrix(rbinom(n * p, 2, rep(af_true, each = n)), nrow = n, ncol = p)

  # Build mock ref_panel
  ref_panel <- data.frame(
    chrom = 1L, pos = seq(1000, by = 100, length.out = p),
    A2 = "A", A1 = "G",
    variant_id = variant_ids,
    allele_freq = colMeans(X) / 2,
    stringsAsFactors = FALSE
  )

  variants_gr <- pecotmr:::.ref_panel_to_granges(ref_panel)
  block_metadata <- S4Vectors::DataFrame(
    region = "chr1:1000-2100", start = 1000L, end = 2100L, chrom = "chr1"
  )
  # Store genotype matrix directly in genotype_handle (matching load_ld_sketch output)
  mock_ld_data <- new("LDData",
    correlation = NULL,
    genotype_handle = X,
    variants = variants_gr,
    snp_idx = NULL,
    block_metadata = block_metadata
  )

  local_mocked_bindings(
    load_LD_matrix = function(ld_meta_file_path, region, return_genotype = FALSE, n_sample = NULL, ...) {
      mock_ld_data
    },
    .package = "pecotmr"
  )

  result <- pecotmr::load_ld_sketch("fake_path.tsv", "chr1:1000-2100")

  # Check structure -- returns an LDData S4 object
  expect_true(is(result, "LDData"))
  result_X <- getGenotypes(result)
  result_ref <- getRefPanel(result)
  result_ids <- getVariantIds(result)
  expect_equal(nrow(result_X), n)
  expect_equal(ncol(result_X), p)
  expect_equal(length(result_ids), p)

  # Raw genotype matrix is returned unchanged
  expect_equal(result_X, X)
})

test_that("load_ld_sketch: removes monomorphic variants", {
  set.seed(66)
  n <- 20
  p <- 5
  variant_ids <- paste0("chr1:", 1:p, ":A:G")

  # Make column 3 monomorphic (all 0)
  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n, ncol = p)
  X[, 3] <- 0  # monomorphic

  ref_panel <- data.frame(
    chrom = 1L, pos = 1:p,
    A2 = "A", A1 = "G",
    variant_id = variant_ids,
    allele_freq = colMeans(X) / 2,
    stringsAsFactors = FALSE
  )

  variants_gr <- pecotmr:::.ref_panel_to_granges(ref_panel)
  block_metadata <- S4Vectors::DataFrame(
    region = "chr1:1-5", start = 1L, end = 5L, chrom = "chr1"
  )
  # Store genotype matrix directly in genotype_handle
  mock_ld_data <- new("LDData",
    correlation = NULL,
    genotype_handle = X,
    variants = variants_gr,
    snp_idx = NULL,
    block_metadata = block_metadata
  )

  local_mocked_bindings(
    load_LD_matrix = function(ld_meta_file_path, region, return_genotype = FALSE, n_sample = NULL, ...) {
      mock_ld_data
    },
    .package = "pecotmr"
  )

  result <- pecotmr::load_ld_sketch("fake_path.tsv", "chr1:1-5")

  # Returns LDData with monomorphic variant removed
  expect_true(is(result, "LDData"))
  result_ids <- getVariantIds(result)
  result_ref <- getRefPanel(result)
  result_X <- getGenotypes(result)
  expect_equal(length(result_ids), p - 1)
  expect_false(variant_ids[3] %in% result_ids)
  expect_equal(nrow(result_ref), p - 1)
  expect_equal(ncol(result_X), p - 1)
})

test_that("SVD from raw sketch matches direct computation", {
  set.seed(77)
  n <- 25
  p <- 8
  af <- runif(p, 0.15, 0.85)
  X <- matrix(rbinom(n * p, 2, rep(af, each = n)), nrow = n, ncol = p)

  # Two-step process: standardize then SVD
  X_std <- pecotmr:::standardize_genotype_hwe(X, af)
  svd_result <- svd(X_std)

  # Verify this matches manual computation
  X_manual <- sweep(sweep(X, 2, 2 * af), 2, sqrt(2 * af * (1 - af)), "/")
  svd_manual <- svd(X_manual)

  expect_equal(svd_result$d, svd_manual$d, tolerance = 1e-14)
  expect_equal(abs(svd_result$v), abs(svd_manual$v), tolerance = 1e-14)
})
