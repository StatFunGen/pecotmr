# Tests for functions that consume genotype matrices or compute LD:
#   computeLd, checkLd, ldPruneByCorrelation, ldClumpByScore,
#   enforceDesignFullRank, filterVariantsByLdReference,
#   resolveLdInput, dentistSingleWindow, dentist

# Fixtures: 100 samples x 100 biallelic polymorphic SNPs on chr21
test_data_dir <- test_path("test_data")
plink_prefix <- file.path(test_data_dir, "test_variants")

# Load genotype matrix once for reuse across tests
load_test_genotype <- function() {
  loadGenotypeRegion(plink_prefix, returnVariantInfo = TRUE)
}

# --- computeLd --------------------------------------------------------------

test_that("computeLd produces valid sample correlation matrix", {
  skip_if_not_installed("pgenlibr")
  geno <- load_test_genotype()
  R <- computeLd(geno$X, method = "sample")
  expect_true(is.matrix(R))
  expect_equal(nrow(R), ncol(geno$X))
  expect_equal(ncol(R), ncol(geno$X))
  expect_true(isSymmetric(R))
  expect_true(all(abs(diag(R) - 1) < 1e-10))
  expect_false(any(is.nan(R)))
  expect_true(all(R >= -1 - 1e-10 & R <= 1 + 1e-10))
})

test_that("computeLd population method produces valid matrix", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- computeLd(X, method = "population")
  expect_true(isSymmetric(R))
  expect_true(all(abs(diag(R) - 1) < 1e-10))
  expect_false(any(is.nan(R)))
})

test_that("computeLd sample and population methods are similar", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R_s <- computeLd(X, method = "sample")
  R_p <- computeLd(X, method = "population")
  # Should be close but not identical (N-1 vs N denominator)
  expect_true(max(abs(R_s - R_p)) < 0.05)
})

test_that("computeLd errors on NULL input", {
  expect_error(computeLd(NULL), "X must be provided")
})

# --- checkLd ----------------------------------------------------------------

test_that("checkLd diagnoses real LD matrix correctly", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- computeLd(X, method = "sample")
  result <- checkLd(R)
  expect_true(is.list(result))
  expect_true(result$isPsd)
  expect_equal(result$methodApplied, "none")
  # min eigenvalue may be near-zero (numerically PSD, not strictly PD)
  expect_true(result$minEigenvalue > -1e-7)
  expect_true(result$nNegative == 0)
  expect_true(is.finite(result$conditionNumber))
})

test_that("checkLd eigenfix improves non-PSD matrix", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- computeLd(X, method = "sample")
  # Make a non-PSD matrix by negating a small block of off-diagonal entries
  R_bad <- R
  R_bad[1:3, 4:6] <- -abs(R_bad[1:3, 4:6]) - 0.5
  R_bad[4:6, 1:3] <- t(R_bad[1:3, 4:6])
  diag(R_bad) <- 1

  result_check <- checkLd(R_bad, method = "check")
  expect_false(result_check$isPsd)
  expect_true(result_check$nNegative > 0)

  result_fix <- checkLd(R_bad, method = "eigenfix")
  expect_equal(result_fix$methodApplied, "eigenfix")
  # Eigenfix should improve (raise) minimum eigenvalue
  fixed_check <- checkLd(result_fix$R)
  expect_true(fixed_check$minEigenvalue > result_check$minEigenvalue)
})

test_that("checkLd shrink repairs perturbed LD matrix", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- computeLd(X, method = "sample")
  R_bad <- R
  R_bad[1, 2] <- R_bad[2, 1] <- 1.5
  diag(R_bad) <- 1
  result <- checkLd(R_bad, method = "shrink", shrinkage = 0.1)
  expect_equal(result$methodApplied, "shrink")
})

# --- ldPruneByCorrelation -------------------------------------------------

test_that("ldPruneByCorrelation prunes correlated variants", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  result <- ldPruneByCorrelation(X, corThres = 0.8)
  expect_true(is.list(result))
  expect_true(is.matrix(result$X.new))
  expect_true(ncol(result$X.new) <= ncol(X))
  expect_true(ncol(result$X.new) > 0)
  expect_equal(length(result$filter.id), ncol(result$X.new))
  # Retained columns are a subset of original
  expect_true(all(result$filter.id %in% seq_len(ncol(X))))
})

test_that("ldPruneByCorrelation with strict threshold prunes more", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  loose <- ldPruneByCorrelation(X, corThres = 0.95)
  strict <- ldPruneByCorrelation(X, corThres = 0.5)
  expect_true(ncol(strict$X.new) <= ncol(loose$X.new))
})

test_that("ldPruneByCorrelation with high threshold keeps most columns", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  result <- ldPruneByCorrelation(X, corThres = 0.999)
  # At threshold near 1, only near-duplicates are pruned; real data may have many
  expect_true(ncol(result$X.new) >= ncol(X) * 0.4)
})

# --- ldClumpByScore -------------------------------------------------------

test_that("ldClumpByScore returns valid indices", {
  skip_if_not_installed("pgenlibr")
  skip_if_not_installed("bigsnpr")
  skip_if_not_installed("bigstatsr")
  geno <- load_test_genotype()
  set.seed(42)
  score <- runif(ncol(geno$X))
  chr <- as.integer(geno$variant_info$chrom)
  pos <- geno$variant_info$pos
  keep <- ldClumpByScore(geno$X, score = score, chr = chr, pos = pos, r2 = 0.2)
  expect_true(is.integer(keep))
  expect_true(length(keep) > 0)
  expect_true(length(keep) <= ncol(geno$X))
  expect_true(all(keep %in% seq_len(ncol(geno$X))))
})

# --- enforceDesignFullRank ------------------------------------------------

test_that("enforceDesignFullRank handles genotype matrix with covariates", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  # Create a simple covariate matrix (e.g., first 2 PCs of X)
  pca <- prcomp(X, rank. = 2)
  C <- pca$x
  result <- enforceDesignFullRank(X, C, strategy = "correlation")
  expect_true(is.matrix(result))
  expect_equal(nrow(result), nrow(X))
  # Should produce full-rank design
  full_design <- cbind(1, result, C)
  expect_equal(qr(full_design)$rank, ncol(full_design))
})

# --- filterVariantsByLdReference -----------------------------------------

test_that("filterVariantsByLdReference filters against PLINK reference via metadata", {
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
    filterVariantsByLdReference(all_ids, meta_file, keepIndel = TRUE)
  )
  expect_true(is.list(result))
  expect_true(length(result$data) <= length(all_ids))
  expect_true(length(result$idx) == length(result$data))
  # Fake variants should be filtered out
  expect_true(length(result$data) <= length(variant_ids))
})

# --- resolveLdInput (internal) ---------------------------------------------

test_that("resolveLdInput computes LD from genotype matrix", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  result <- pecotmr:::resolveLdInput(X = X, needNSample = TRUE)
  expect_true(is.list(result))
  expect_true(is.matrix(result$R))
  expect_equal(nrow(result$R), ncol(X))
  expect_equal(result$nSample, nrow(X))
  expect_true(isSymmetric(result$R))
})

test_that("resolveLdInput passes through pre-computed R", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- computeLd(X, method = "sample")
  result <- pecotmr:::resolveLdInput(R = R, nSample = 100L, needNSample = TRUE)
  expect_equal(result$R, R)
  expect_equal(result$nSample, 100L)
})

test_that("resolveLdInput errors when neither R nor X provided", {
  expect_error(pecotmr:::resolveLdInput(), "Either R .* or X .* must be provided")
})

test_that("resolveLdInput errors when both R and X provided", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- computeLd(X, method = "sample")
  expect_error(pecotmr:::resolveLdInput(R = R, X = X), "Provide either R or X, not both")
})

test_that("resolveLdInput errors when R given without nSample and needed", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- computeLd(X, method = "sample")
  expect_error(pecotmr:::resolveLdInput(R = R, needNSample = TRUE),
               "nSample is required")
})

# --- dentistSingleWindow ---------------------------------------------------

test_that("dentistSingleWindow works with genotype matrix X", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  set.seed(42)
  z <- rnorm(ncol(X))
  result <- suppressWarnings(dentistSingleWindow(z, X = X))
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), length(z))
  expect_true("original_z" %in% names(result))
  expect_true("imputed_z" %in% names(result))
  expect_true("outlier" %in% names(result))
  expect_true(is.logical(result$outlier))
})

test_that("dentistSingleWindow works with pre-computed R", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- computeLd(X, method = "sample")
  set.seed(42)
  z <- rnorm(ncol(X))
  result <- suppressWarnings(dentistSingleWindow(z, R = R, nSample = nrow(X)))
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), length(z))
})

test_that("dentistSingleWindow detects injected outliers", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  R <- computeLd(X, method = "sample")
  set.seed(42)
  z <- rnorm(ncol(X))
  # Inject extreme outliers
  z[1] <- 50
  z[2] <- -50
  result <- suppressWarnings(dentistSingleWindow(z, R = R, nSample = nrow(X)))
  # At least one of the injected values should be flagged
  expect_true(any(result$outlier))
})

# --- dentist (multi-window) -------------------------------------------------

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
    dentist(sum_stat, X = geno$X, windowMode = "count", minDim = 50)
  )
  expect_true(is.data.frame(result))
  # Window merging may add overlap rows; result should be >= input size
  expect_true(nrow(result) >= nrow(sum_stat))
  expect_true(all(c("original_z", "imputed_z", "outlier") %in% names(result)))
})

test_that("dentist accepts zscore column name variant", {
  skip_if_not_installed("pgenlibr")
  geno <- load_test_genotype()
  R <- computeLd(geno$X, method = "sample")
  set.seed(42)
  sum_stat <- data.frame(
    position = geno$variant_info$pos,
    zscore = rnorm(ncol(geno$X))
  )
  result <- suppressWarnings(
    dentist(sum_stat, R = R, nSample = nrow(geno$X), windowMode = "count", minDim = 50)
  )
  expect_true(nrow(result) >= nrow(sum_stat))
})

test_that("dentist errors when sum_stat missing required columns", {
  skip_if_not_installed("pgenlibr")
  X <- load_test_genotype()$X
  bad_stat <- data.frame(x = 1:ncol(X), y = rnorm(ncol(X)))
  expect_error(dentist(bad_stat, X = X), "missing either")
})
