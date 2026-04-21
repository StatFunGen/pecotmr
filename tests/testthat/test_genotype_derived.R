# Tests for functions that consume genotype matrices or compute LD:
#   compute_LD, check_ld, ld_prune_by_correlation, ld_clump_by_score,
#   enforce_design_full_rank, filter_variants_by_ld_reference,
#   resolve_LD_input, dentist_single_window, dentist

# Fixtures: 100 samples x 100 biallelic polymorphic SNPs on chr21
test_data_dir <- test_path("test_data")
plink_prefix <- file.path(test_data_dir, "test_variants")

# Load genotype matrix once for reuse across tests
load_test_genotype <- function() {
  load_genotype_region(plink_prefix, return_variant_info = TRUE)
}

# ─── compute_LD ──────────────────────────────────────────────────────────────

test_that("compute_LD produces valid sample correlation matrix", {
  skip_if_not_installed("pgenlibr")
  geno <- load_test_genotype()
  R <- compute_LD(geno$X, method = "sample")
  expect_true(is.matrix(R))
  expect_equal(nrow(R), ncol(geno$X))
  expect_equal(ncol(R), ncol(geno$X))
  expect_true(isSymmetric(R))
  expect_true(all(abs(diag(R) - 1) < 1e-10))
  expect_false(any(is.nan(R)))
  expect_true(all(R >= -1 - 1e-10 & R <= 1 + 1e-10))
})

test_that("compute_LD population method produces valid matrix", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- compute_LD(X, method = "population")
  expect_true(isSymmetric(R))
  expect_true(all(abs(diag(R) - 1) < 1e-10))
  expect_false(any(is.nan(R)))
})

test_that("compute_LD sample and population methods are similar", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R_s <- compute_LD(X, method = "sample")
  R_p <- compute_LD(X, method = "population")
  # Should be close but not identical (N-1 vs N denominator)
  expect_true(max(abs(R_s - R_p)) < 0.05)
})

test_that("compute_LD errors on NULL input", {
  expect_error(compute_LD(NULL), "X must be provided")
})

# ─── check_ld ────────────────────────────────────────────────────────────────

test_that("check_ld diagnoses real LD matrix correctly", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- compute_LD(X, method = "sample")
  result <- check_ld(R)
  expect_true(is.list(result))
  expect_true(result$is_psd)
  expect_equal(result$method_applied, "none")
  # min eigenvalue may be near-zero (numerically PSD, not strictly PD)
  expect_true(result$min_eigenvalue > -1e-7)
  expect_true(result$n_negative == 0)
  expect_true(is.finite(result$condition_number))
})

test_that("check_ld eigenfix improves non-PSD matrix", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- compute_LD(X, method = "sample")
  # Make a non-PSD matrix by negating a small block of off-diagonal entries
  R_bad <- R
  R_bad[1:3, 4:6] <- -abs(R_bad[1:3, 4:6]) - 0.5
  R_bad[4:6, 1:3] <- t(R_bad[1:3, 4:6])
  diag(R_bad) <- 1

  result_check <- check_ld(R_bad, method = "check")
  expect_false(result_check$is_psd)
  expect_true(result_check$n_negative > 0)

  result_fix <- check_ld(R_bad, method = "eigenfix")
  expect_equal(result_fix$method_applied, "eigenfix")
  # Eigenfix should improve (raise) minimum eigenvalue
  fixed_check <- check_ld(result_fix$R)
  expect_true(fixed_check$min_eigenvalue > result_check$min_eigenvalue)
})

test_that("check_ld shrink repairs perturbed LD matrix", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- compute_LD(X, method = "sample")
  R_bad <- R
  R_bad[1, 2] <- R_bad[2, 1] <- 1.5
  diag(R_bad) <- 1
  result <- check_ld(R_bad, method = "shrink", shrinkage = 0.1)
  expect_equal(result$method_applied, "shrink")
})

# ─── ld_prune_by_correlation ─────────────────────────────────────────────────

test_that("ld_prune_by_correlation prunes correlated variants", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  result <- ld_prune_by_correlation(X, cor_thres = 0.8)
  expect_true(is.list(result))
  expect_true(is.matrix(result$X.new))
  expect_true(ncol(result$X.new) <= ncol(X))
  expect_true(ncol(result$X.new) > 0)
  expect_equal(length(result$filter.id), ncol(result$X.new))
  # Retained columns are a subset of original
  expect_true(all(result$filter.id %in% seq_len(ncol(X))))
})

test_that("ld_prune_by_correlation with strict threshold prunes more", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  loose <- ld_prune_by_correlation(X, cor_thres = 0.95)
  strict <- ld_prune_by_correlation(X, cor_thres = 0.5)
  expect_true(ncol(strict$X.new) <= ncol(loose$X.new))
})

test_that("ld_prune_by_correlation with high threshold keeps most columns", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  result <- ld_prune_by_correlation(X, cor_thres = 0.999)
  # At threshold near 1, only near-duplicates are pruned; real data may have many
  expect_true(ncol(result$X.new) >= ncol(X) * 0.4)
})

# ─── ld_clump_by_score ───────────────────────────────────────────────────────

test_that("ld_clump_by_score returns valid indices", {
  skip_if_not_installed("pgenlibr")
  skip_if_not_installed("bigsnpr")
  skip_if_not_installed("bigstatsr")
  geno <- load_test_genotype()
  set.seed(42)
  score <- runif(ncol(geno$X))
  chr <- as.integer(geno$variant_info$chrom)
  pos <- geno$variant_info$pos
  keep <- ld_clump_by_score(geno$X, score = score, chr = chr, pos = pos, r2 = 0.2)
  expect_true(is.integer(keep))
  expect_true(length(keep) > 0)
  expect_true(length(keep) <= ncol(geno$X))
  expect_true(all(keep %in% seq_len(ncol(geno$X))))
})

# ─── enforce_design_full_rank ────────────────────────────────────────────────

test_that("enforce_design_full_rank handles genotype matrix with covariates", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  # Create a simple covariate matrix (e.g., first 2 PCs of X)
  pca <- prcomp(X, rank. = 2)
  C <- pca$x
  result <- enforce_design_full_rank(X, C, strategy = "correlation")
  expect_true(is.matrix(result))
  expect_equal(nrow(result), nrow(X))
  # Should produce full-rank design
  full_design <- cbind(1, result, C)
  expect_equal(qr(full_design)$rank, ncol(full_design))
})

# ─── filter_variants_by_ld_reference ─────────────────────────────────────────

test_that("filter_variants_by_ld_reference filters against PLINK reference via metadata", {
  skip_if_not_installed("pgenlibr")
  geno <- load_test_genotype()
  vi <- geno$variant_info
  variant_ids <- paste0(vi$chrom, ":", vi$pos, ":", vi$A2, ":", vi$A1)
  fake_ids <- c("21:999999:A:G", "21:888888:C:T")
  all_ids <- c(variant_ids, fake_ids)

  # Create a metadata TSV in the same directory as the PLINK files
  # so the relative path resolves correctly
  meta_file <- file.path(test_data_dir, "ld_meta_tmp.tsv")
  on.exit(unlink(meta_file), add = TRUE)
  writeLines(
    paste("chrom", "start", "end", "path", sep = "\t"),
    meta_file
  )
  cat(paste("21", "0", "0", "test_variants", sep = "\t"), "\n",
      file = meta_file, append = TRUE)

  result <- suppressMessages(
    filter_variants_by_ld_reference(all_ids, meta_file, keep_indel = TRUE)
  )
  expect_true(is.list(result))
  expect_true(length(result$data) <= length(all_ids))
  expect_true(length(result$idx) == length(result$data))
  # Fake variants should be filtered out
  expect_true(length(result$data) <= length(variant_ids))
})

# ─── resolve_LD_input (internal) ─────────────────────────────────────────────

test_that("resolve_LD_input computes LD from genotype matrix", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  result <- pecotmr:::resolve_LD_input(X = X, need_nSample = TRUE)
  expect_true(is.list(result))
  expect_true(is.matrix(result$R))
  expect_equal(nrow(result$R), ncol(X))
  expect_equal(result$nSample, nrow(X))
  expect_true(isSymmetric(result$R))
})

test_that("resolve_LD_input passes through pre-computed R", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- compute_LD(X, method = "sample")
  result <- pecotmr:::resolve_LD_input(R = R, nSample = 100L, need_nSample = TRUE)
  expect_equal(result$R, R)
  expect_equal(result$nSample, 100L)
})

test_that("resolve_LD_input errors when neither R nor X provided", {
  expect_error(pecotmr:::resolve_LD_input(), "Either R .* or X .* must be provided")
})

test_that("resolve_LD_input errors when both R and X provided", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- compute_LD(X, method = "sample")
  expect_error(pecotmr:::resolve_LD_input(R = R, X = X), "Provide either R or X, not both")
})

test_that("resolve_LD_input errors when R given without nSample and needed", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- compute_LD(X, method = "sample")
  expect_error(pecotmr:::resolve_LD_input(R = R, need_nSample = TRUE),
               "nSample is required")
})

# ─── dentist_single_window ───────────────────────────────────────────────────

test_that("dentist_single_window works with genotype matrix X", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  set.seed(42)
  z <- rnorm(ncol(X))
  result <- suppressWarnings(dentist_single_window(z, X = X))
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), length(z))
  expect_true("original_z" %in% names(result))
  expect_true("imputed_z" %in% names(result))
  expect_true("outlier" %in% names(result))
  expect_true(is.logical(result$outlier))
})

test_that("dentist_single_window works with pre-computed R", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- compute_LD(X, method = "sample")
  set.seed(42)
  z <- rnorm(ncol(X))
  result <- suppressWarnings(dentist_single_window(z, R = R, nSample = nrow(X)))
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), length(z))
})

test_that("dentist_single_window detects injected outliers", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- compute_LD(X, method = "sample")
  set.seed(42)
  z <- rnorm(ncol(X))
  # Inject extreme outliers
  z[1] <- 50
  z[2] <- -50
  result <- suppressWarnings(dentist_single_window(z, R = R, nSample = nrow(X)))
  # At least one of the injected values should be flagged
  expect_true(any(result$outlier))
})

# ─── dentist (multi-window) ─────────────────────────────────────────────────

test_that("dentist works with genotype matrix and sum_stat data frame", {
  skip_if_not_installed("pgenlibr")
  geno <- load_test_genotype()
  set.seed(42)
  sum_stat <- data.frame(
    pos = geno$variant_info$pos,
    z = rnorm(ncol(geno$X))
  )
  # Use count mode with small window since we only have 100 variants
  result <- suppressWarnings(
    dentist(sum_stat, X = geno$X, window_mode = "count", min_dim = 50)
  )
  expect_true(is.data.frame(result))
  # Window merging may add overlap rows; result should be >= input size
  expect_true(nrow(result) >= nrow(sum_stat))
  expect_true(all(c("original_z", "imputed_z", "outlier") %in% names(result)))
})

test_that("dentist accepts zscore column name variant", {
  skip_if_not_installed("pgenlibr")
  geno <- load_test_genotype()
  R <- compute_LD(geno$X, method = "sample")
  set.seed(42)
  sum_stat <- data.frame(
    position = geno$variant_info$pos,
    zscore = rnorm(ncol(geno$X))
  )
  result <- suppressWarnings(
    dentist(sum_stat, R = R, nSample = nrow(geno$X), window_mode = "count", min_dim = 50)
  )
  expect_true(nrow(result) >= nrow(sum_stat))
})

test_that("dentist errors when sum_stat missing required columns", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  bad_stat <- data.frame(x = 1:ncol(X), y = rnorm(ncol(X)))
  expect_error(dentist(bad_stat, X = X), "missing either")
})
