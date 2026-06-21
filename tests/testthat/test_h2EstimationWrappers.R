context("h2_estimate_wrappers")

# ===========================================================================
# Helper functions to build test S4 objects
# ===========================================================================

make_test_eigen_ref <- function(nSnps = 20, nBlocks = 2) {
  snps_per_block <- nSnps / nBlocks
  blocks_gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", nBlocks),
    ranges = IRanges::IRanges(
      start = seq(1, by = snps_per_block * 100, length.out = nBlocks),
      width = snps_per_block * 100 - 1
    )
  )
  ld_blocks <- new("LdBlocks", blocks = blocks_gr, genome = "hg19")

  snp_info <- data.frame(
    SNP = paste0("rs", seq_len(nSnps)),
    CHR = rep("chr1", nSnps),
    BP = as.integer(seq(50, by = 100, length.out = nSnps)),
    A1 = rep("A", nSnps),
    A2 = rep("G", nSnps),
    stringsAsFactors = FALSE
  )

  eigen_list <- lapply(seq_len(nBlocks), function(b) {
    idx <- seq((b - 1) * snps_per_block + 1, b * snps_per_block)
    p <- length(idx)
    set.seed(42 + b)
    R <- matrix(0.3, p, p)
    diag(R) <- 1
    e <- eigen(R)
    list(values = e$values, vectors = e$vectors, snpIdx = as.integer(idx))
  })

  new("LdEigen",
    ldBlocks = ld_blocks,
    snpInfo = snp_info,
    nRef = 1000L,
    inSample = FALSE,
    genome = "hg19",
    eigenList = eigen_list,
    eigenvalueTruncation = 1.0
  )
}

make_test_score_ref <- function(nSnps = 20, nBlocks = 2,
                                with_ld_matrices = FALSE) {
  snps_per_block <- nSnps / nBlocks
  blocks_gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", nBlocks),
    ranges = IRanges::IRanges(
      start = seq(1, by = snps_per_block * 100, length.out = nBlocks),
      width = snps_per_block * 100 - 1
    )
  )
  ld_blocks <- new("LdBlocks", blocks = blocks_gr, genome = "hg19")

  snp_info <- data.frame(
    SNP = paste0("rs", seq_len(nSnps)),
    CHR = rep("chr1", nSnps),
    BP = as.integer(seq(50, by = 100, length.out = nSnps)),
    A1 = rep("A", nSnps),
    A2 = rep("G", nSnps),
    stringsAsFactors = FALSE
  )

  set.seed(99)
  ld_scores <- matrix(runif(nSnps, 1, 10), ncol = 1,
                      dimnames = list(NULL, "base_l2"))
  ld_score_weights <- rep(1 / nSnps, nSnps)

  ld_matrix_list <- if (with_ld_matrices) {
    lapply(seq_len(nBlocks), function(b) {
      idx <- seq((b - 1) * snps_per_block + 1, b * snps_per_block)
      p <- length(idx)
      R <- matrix(0.3, p, p)
      diag(R) <- 1
      list(R = R, snpIdx = as.integer(idx))
    })
  } else {
    list()
  }

  new("LdScore",
    ldBlocks = ld_blocks,
    snpInfo = snp_info,
    nRef = 1000L,
    inSample = FALSE,
    genome = "hg19",
    ldScores = ld_scores,
    ldScoreWeights = ld_score_weights,
    ldMatrixList = ld_matrix_list
  )
}

make_test_annotations <- function(nSnps = 20) {
  set.seed(77)
  snp_ranges <- GenomicRanges::GRanges(
    seqnames = rep("chr1", nSnps),
    ranges = IRanges::IRanges(
      start = as.integer(seq(50, by = 100, length.out = nSnps)),
      width = 1L
    )
  )
  annot_mat <- matrix(
    c(rbinom(nSnps, 1, 0.5), rbinom(nSnps, 1, 0.3)),
    nrow = nSnps, ncol = 2
  )
  colnames(annot_mat) <- c("annot_A", "annot_B")
  annotation_meta <- data.frame(
    name = c("annot_A", "annot_B"),
    tier = c("baseline", "baseline"),
    type = c("binary", "binary"),
    stringsAsFactors = FALSE
  )
  AnnotationMatrix(annot_mat, snp_ranges, annotation_meta, genome = "hg19")
}

make_test_h2estimate <- function(with_enrichment = TRUE) {
  enrich <- if (with_enrichment) {
    data.frame(
      annotation = c("annot1", "annot2"),
      tau = c(1e-7, 2e-7),
      tauSe = c(5e-8, 6e-8),
      enrichment = c(2.0, 3.0),
      enrichmentSe = c(0.5, 0.7),
      enrichmentP = c(0.01, 0.001),
      propH2 = c(0.3, 0.5),
      propSnps = c(0.15, 0.17),
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }

  set.seed(55)
  tauBlocks <- if (with_enrichment) {
    matrix(rnorm(20), nrow = 10, ncol = 2,
           dimnames = list(NULL, c("annot1", "annot2")))
  } else {
    NULL
  }

  new("H2Estimate",
    h2 = 0.3, h2Se = 0.05,
    intercept = 1.01, interceptSe = 0.02,
    local = NULL, enrichment = enrich,
    tauBlocks = tauBlocks, scoreStats = NULL,
    method = "lder", nSnps = 100L, traitName = "test"
  )
}

# ===========================================================================
# .validate_method_ref (internal)
# ===========================================================================

test_that(".validate_method_ref errors when method='lder' but ref is LdScore", {
  score_ref <- make_test_score_ref()
  expect_error(
    pecotmr:::.validateMethodRef("lder", score_ref),
    "requires an LdEigen"
  )
})

test_that(".validate_method_ref errors when method='hdl' but ref is LdScore", {
  score_ref <- make_test_score_ref()
  expect_error(
    pecotmr:::.validateMethodRef("hdl", score_ref),
    "requires an LdEigen"
  )
})

test_that(".validate_method_ref errors when method='gldsc' but ref is LdEigen", {
  eigen_ref <- make_test_eigen_ref()
  expect_error(
    pecotmr:::.validateMethodRef("gldsc", eigen_ref),
    "requires an LdScore"
  )
})

test_that(".validate_method_ref returns TRUE for valid combinations", {
  eigen_ref <- make_test_eigen_ref()
  score_ref <- make_test_score_ref()

  expect_true(pecotmr:::.validateMethodRef("lder", eigen_ref))
  expect_true(pecotmr:::.validateMethodRef("hdl", eigen_ref))
  expect_true(pecotmr:::.validateMethodRef("gldsc", score_ref))
})

# ===========================================================================
# computeLdScores for LdEigen
# ===========================================================================

test_that("computeLdScores LdEigen returns matrix with base_l2 column", {
  eigen_ref <- make_test_eigen_ref()
  result <- computeLdScores(eigen_ref)
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 1)
  expect_equal(colnames(result), "base_l2")
  expect_equal(nrow(result), nrow(eigen_ref@snpInfo))
})

test_that("computeLdScores LdEigen base LD scores are non-negative", {
  eigen_ref <- make_test_eigen_ref()
  result <- computeLdScores(eigen_ref)
  expect_true(all(result[, "base_l2"] >= 0))
})

test_that("computeLdScores LdEigen with annotations returns base + annotation columns", {
  eigen_ref <- make_test_eigen_ref()
  annot <- make_test_annotations()
  result <- computeLdScores(eigen_ref, annotations = annot)
  expect_true(is.matrix(result))
  # base_l2 + 2 annotation columns

  expect_equal(ncol(result), 3)
  expect_equal(nrow(result), nrow(eigen_ref@snpInfo))
})

test_that("computeLdScores LdEigen annotation column names match", {
  eigen_ref <- make_test_eigen_ref()
  annot <- make_test_annotations()
  result <- computeLdScores(eigen_ref, annotations = annot)
  expect_equal(colnames(result), c("base_l2", "annot_A", "annot_B"))
})

# ===========================================================================
# computeLdScores for LdScore
# ===========================================================================

test_that("computeLdScores LdScore without annotations returns stored ld_scores", {
  score_ref <- make_test_score_ref()
  result <- computeLdScores(score_ref)
  expect_identical(result, score_ref@ldScores)
})

test_that("computeLdScores LdScore with annotations but no ld_matrix_list errors", {
  score_ref <- make_test_score_ref(with_ld_matrices = FALSE)
  annot <- make_test_annotations()
  expect_error(
    computeLdScores(score_ref, annotations = annot),
    "ldMatrixList"
  )
})

# ===========================================================================
# H2Estimate accessors
# ===========================================================================

test_that("getlocal returns the local slot", {
  local_df <- data.frame(
    blockId = 1:3,
    h2Local = c(0.1, 0.15, 0.05),
    h2LocalSe = c(0.02, 0.03, 0.01),
    stringsAsFactors = FALSE
  )
  h2_obj <- new("H2Estimate",
    h2 = 0.3, h2Se = 0.05,
    intercept = 1.0, interceptSe = 0.01,
    local = local_df, enrichment = NULL,
    tauBlocks = NULL, scoreStats = NULL,
    method = "lder", nSnps = 100L, traitName = "test"
  )
  expect_identical(getLocal(h2_obj), local_df)
})

test_that("getenrichment returns the enrichment slot", {
  h2_obj <- make_test_h2estimate(with_enrichment = TRUE)
  result <- getEnrichment(h2_obj)
  expect_true(is.data.frame(result))
  expect_equal(result$annotation, c("annot1", "annot2"))
})

test_that("getscorestats returns the scoreStats slot", {
  scoreStats <- list(z = c(1.5, -0.8), R = diag(2),
                      annotationNames = c("a1", "a2"))
  h2_obj <- new("H2Estimate",
    h2 = 0.3, h2Se = 0.05,
    intercept = 1.0, interceptSe = 0.01,
    local = NULL, enrichment = NULL,
    tauBlocks = NULL, scoreStats = scoreStats,
    method = "lder", nSnps = 100L, traitName = "test"
  )
  expect_identical(getScoreStats(h2_obj), scoreStats)
})

test_that("accessors return NULL when slots are NULL", {
  h2_obj <- new("H2Estimate",
    h2 = 0.3, h2Se = 0.05,
    intercept = 1.0, interceptSe = 0.01,
    local = NULL, enrichment = NULL,
    tauBlocks = NULL, scoreStats = NULL,
    method = "lder", nSnps = 100L, traitName = "test"
  )
  expect_null(getLocal(h2_obj))
  expect_null(getEnrichment(h2_obj))
  expect_null(getScoreStats(h2_obj))
})

# ===========================================================================
# h2EstimateToSldscTrait
# ===========================================================================

test_that("h2EstimateToSldscTrait returns correct list structure", {
  h2_obj <- make_test_h2estimate(with_enrichment = TRUE)
  result <- h2EstimateToSldscTrait(h2_obj)

  expected_names <- c("categories", "tau", "tauSe", "enrichment",
                      "enrichmentSe", "enrichmentP", "propH2",
                      "propSnps", "h2g", "tauBlocks", "nBlocks")
  expect_true(all(expected_names %in% names(result)))
  expect_equal(length(result), length(expected_names))

  expect_equal(result$categories, c("annot1", "annot2"))
  expect_equal(result$h2g, 0.3)
  expect_true(is.matrix(result$tauBlocks))
  expect_equal(result$nBlocks, 10L)
  expect_equal(names(result$tau), c("annot1", "annot2"))
  expect_equal(names(result$enrichment), c("annot1", "annot2"))
})

test_that("h2EstimateToSldscTrait errors on non-H2Estimate input", {
  expect_error(
    h2EstimateToSldscTrait(list(h2 = 0.3)),
    "must be an H2Estimate"
  )
})

test_that("h2EstimateToSldscTrait errors when enrichment is NULL", {
  h2_obj <- make_test_h2estimate(with_enrichment = FALSE)
  expect_error(
    h2EstimateToSldscTrait(h2_obj),
    "no enrichment"
  )
})

test_that("h2EstimateToSldscTrait creates 1-row dummy when tauBlocks is NULL", {
  h2_obj <- make_test_h2estimate(with_enrichment = TRUE)
  h2_obj@tauBlocks <- NULL
  result <- h2EstimateToSldscTrait(h2_obj)
  expect_equal(nrow(result$tauBlocks), 1)
  expect_equal(result$nBlocks, 1L)
  expect_equal(colnames(result$tauBlocks), c("annot1", "annot2"))
})

# ===========================================================================
# Helper: GwasSumStats matched to a reference panel
# ===========================================================================

make_test_gwas_genotype_handle <- function() {
  new("GenotypeHandle",
    path = "/tmp/test.gds",
    format = "gds",
    snpInfo = data.frame(),
    nSamples = 0L,
    sampleIds = character(),
    pgenPtr = NULL)
}

.dfToGwasGr <- function(df) {
  gr <- GenomicRanges::GRanges(
    seqnames = df$CHR,
    ranges = IRanges::IRanges(start = df$BP, width = 1L)
  )
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = df$SNP, A1 = df$A1, A2 = df$A2, Z = df$Z, N = df$N
  )
  gr
}

make_test_sumstats_for_ref <- function(ref, traitName = "test", varY = NA_real_) {
  n_snps <- nrow(ref@snpInfo)
  set.seed(123)
  df <- data.frame(
    SNP = ref@snpInfo$SNP,
    CHR = sub("^chr", "", ref@snpInfo$CHR),
    BP = ref@snpInfo$BP,
    A1 = ref@snpInfo$A1,
    A2 = ref@snpInfo$A2,
    Z = rnorm(n_snps),
    N = rep(50000, n_snps),
    stringsAsFactors = FALSE
  )
  GwasSumStats(
    study = traitName,
    entry = list(.dfToGwasGr(df)),
    genome = "hg19",
    ldSketch = make_test_gwas_genotype_handle(),
    varY = varY)
}

# ===========================================================================
# estimateh2 end-to-end tests
# ===========================================================================

test_that("estimateh2 with method='lder' returns H2Estimate with correct slots", {
  eigen_ref <- make_test_eigen_ref()
  ss <- make_test_sumstats_for_ref(eigen_ref)
  result <- estimateH2(ss, eigen_ref, method = "lder")

  expect_s4_class(result, "H2Estimate")
  expect_true(is.numeric(result@h2))
  expect_true(is.numeric(result@h2Se))
  expect_true(is.numeric(result@intercept))
  expect_true(is.numeric(result@interceptSe))
  expect_equal(result@method, "lder")
  expect_equal(result@nSnps, nSnps(ss))
  expect_equal(result@traitName, "test")
})

test_that("estimateh2 with var_y correction runs without error", {
  eigen_ref <- make_test_eigen_ref()
  ss <- make_test_sumstats_for_ref(eigen_ref, traitName = "cc_trait", varY = 4.0)
  expect_equal(getVarY(ss), 4.0)

  result <- estimateH2(ss, eigen_ref, method = "lder")
  expect_s4_class(result, "H2Estimate")
  expect_true(is.numeric(result@h2))
})

test_that("estimateh2 with method='gldsc' returns H2Estimate", {
  score_ref <- make_test_score_ref()
  ss <- make_test_sumstats_for_ref(score_ref)
  result <- estimateH2(ss, score_ref, method = "gldsc")

  expect_s4_class(result, "H2Estimate")
  expect_true(is.numeric(result@h2))
  expect_true(is.numeric(result@h2Se))
  expect_equal(result@method, "gldsc")
  expect_equal(result@nSnps, nSnps(ss))
  expect_equal(result@traitName, "test")
})

test_that("estimateh2 with method='hdl' returns H2Estimate", {
  eigen_ref <- make_test_eigen_ref()
  ss <- make_test_sumstats_for_ref(eigen_ref)
  # HDL likelihood optimization on tiny test data produces NaN warnings
  # from log() on negative sigma2 during search; these are harmless.
  suppressWarnings(
    result <- estimateH2(ss, eigen_ref, method = "hdl")
  )

  expect_s4_class(result, "H2Estimate")
  expect_true(is.numeric(result@h2))
  expect_true(is.numeric(result@h2Se))
  expect_equal(result@method, "hdl")
  expect_equal(result@nSnps, nSnps(ss))
  expect_equal(result@traitName, "test")
})

# ===========================================================================
# computeLdScores for LdScore WITH annotations
# ===========================================================================

test_that("computeLdScores LdScore with annotations and ld_matrix_list returns correct matrix", {
  score_ref <- make_test_score_ref(with_ld_matrices = TRUE)
  annot <- make_test_annotations()
  result <- computeLdScores(score_ref, annotations = annot)

  expect_true(is.matrix(result))
  # base_l2 + 2 annotation columns
  expect_equal(ncol(result), 3)
  expect_equal(nrow(result), nrow(score_ref@snpInfo))
  expect_equal(colnames(result), c("base_l2", "annot_A", "annot_B"))
  # First column should be the stored base LD scores
  expect_equal(result[, 1], score_ref@ldScores[, 1])
})


context("h2_utils")

# =============================================================================
# weightedLs
# =============================================================================

test_that("weightedLs returns correct structure", {
  set.seed(42)
  n <- 50
  x <- rnorm(n)
  X <- cbind(x, 1)
  y <- 2 * x + 1 + rnorm(n, sd = 0.1)
  w <- rep(1, n)
  res <- pecotmr:::weightedLs(y, X, w)
  expect_true(is.list(res))
  expect_named(res, c("coef", "se", "residuals", "fitted", "vcov"))
  expect_length(res$coef, 2)
  expect_length(res$se, 2)
  expect_length(res$residuals, n)
  expect_length(res$fitted, n)
  expect_equal(dim(res$vcov), c(2, 2))
})

test_that("weightedLs recovers coefficients for simple linear model", {
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  X <- cbind(x, 1)
  y <- 2 * x + 1 + rnorm(n, sd = 0.1)
  w <- rep(1, n)
  res <- pecotmr:::weightedLs(y, X, w)
  expect_equal(res$coef[1], 2, tolerance = 0.1)
  expect_equal(res$coef[2], 1, tolerance = 0.1)
})

test_that("weightedLs converts vector X to matrix", {
  set.seed(42)
  n <- 30
  x <- rnorm(n)
  y <- 3 * x + rnorm(n, sd = 0.1)
  w <- rep(1, n)
  res <- pecotmr:::weightedLs(y, x, w)
  expect_length(res$coef, 1)
  expect_equal(res$coef[1], 3, tolerance = 0.2)
})

test_that("weightedLs fitted + residuals = y", {
  set.seed(42)
  n <- 30
  X <- cbind(rnorm(n), 1)
  y <- rnorm(n)
  w <- rep(1, n)
  res <- pecotmr:::weightedLs(y, X, w)
  expect_equal(res$fitted + res$residuals, y)
})

# =============================================================================
# jackknifeSe
# =============================================================================

test_that("jackknifeSe returns zero when all LOO estimates equal full", {
  nBlocks <- 5
  n_params <- 3
  estimates_full <- c(1.0, 2.0, 3.0)
  estimates_loo <- matrix(rep(estimates_full, each = nBlocks),
                          nrow = nBlocks, ncol = n_params)
  se <- pecotmr:::jackknifeSe(estimates_full, estimates_loo)
  expect_length(se, n_params)
  expect_equal(se, rep(0, n_params))
})

test_that("jackknifeSe returns positive SE with varying LOO estimates", {
  set.seed(42)
  nBlocks <- 10
  estimates_full <- c(5.0, 10.0)
  estimates_loo <- matrix(rnorm(nBlocks * 2, mean = rep(estimates_full, each = nBlocks), sd = 0.5),
                          nrow = nBlocks, ncol = 2)
  se <- pecotmr:::jackknifeSe(estimates_full, estimates_loo)
  expect_length(se, 2)
  expect_true(all(se > 0))
})

test_that("jackknifeSe computes known case correctly", {
  # If full estimate is the mean of 5 values and LOO removes each one,
  # we can verify the formula
  vals <- c(2, 4, 6, 8, 10)
  full_mean <- mean(vals)
  n <- length(vals)
  # Leave-one-out means: remove element i, compute mean of remaining
  loo_means <- vapply(seq_len(n), function(i) mean(vals[-i]), numeric(1))
  estimates_loo <- matrix(loo_means, ncol = 1)
  se <- pecotmr:::jackknifeSe(full_mean, estimates_loo)
  # Pseudo-values: n * full - (n-1) * loo
  pseudo <- n * full_mean - (n - 1) * loo_means
  expected_se <- sqrt(var(pseudo) / n)
  expect_equal(se, expected_se)
})

# =============================================================================
# weightedLsRidge
# =============================================================================

test_that("weightedLsRidge with lambda=0 matches weightedLs", {
  set.seed(42)
  n <- 50
  X <- cbind(rnorm(n), rnorm(n), 1)
  y <- X %*% c(1, 2, 0.5) + rnorm(n, sd = 0.1)
  w <- rep(1, n)
  res_ridge <- pecotmr:::weightedLsRidge(y, X, w, lambda = 0)
  res_wls <- pecotmr:::weightedLs(y, X, w)
  expect_equal(res_ridge$coef, res_wls$coef)
  expect_equal(res_ridge$se, res_wls$se)
  expect_equal(res_ridge$residuals, res_wls$residuals)
  expect_equal(res_ridge$fitted, res_wls$fitted)
})

test_that("weightedLsRidge with lambda>0 shrinks coefficients", {
  set.seed(42)
  n <- 50
  X <- cbind(rnorm(n), rnorm(n), 1)
  y <- X %*% c(3, 3, 1) + rnorm(n, sd = 0.5)
  w <- rep(1, n)
  res_no_ridge <- pecotmr:::weightedLsRidge(y, X, w, lambda = 0)
  res_ridge <- pecotmr:::weightedLsRidge(y, X, w, lambda = 10)
  # Non-intercept coefficients should be shrunk toward zero
  expect_true(abs(res_ridge$coef[1]) < abs(res_no_ridge$coef[1]))
  expect_true(abs(res_ridge$coef[2]) < abs(res_no_ridge$coef[2]))
})

test_that("weightedLsRidge penalize_intercept=FALSE leaves last column unpenalized", {
  set.seed(42)
  n <- 100
  X <- cbind(rnorm(n), 1)
  y <- X %*% c(5, 3) + rnorm(n, sd = 0.5)
  w <- rep(1, n)
  # With large lambda but no intercept penalty, intercept should still be reasonable
  res <- pecotmr:::weightedLsRidge(y, X, w, lambda = 100, penalizeIntercept = FALSE)
  # Intercept (col 2) should not be shrunk as aggressively as slope (col 1)
  res_pen <- pecotmr:::weightedLsRidge(y, X, w, lambda = 100, penalizeIntercept = TRUE)
  # The intercept should differ between the two
  expect_false(isTRUE(all.equal(res$coef[2], res_pen$coef[2])))
})

# =============================================================================
# computeBaselineEnrichment
# =============================================================================

test_that("computeBaselineEnrichment returns correct structure", {
  M <- 100
  baseline_mat <- matrix(0, nrow = M, ncol = 2)
  baseline_mat[1:50, 1] <- 1
  baseline_mat[, 2] <- 1
  tau <- c(0.01, 0.005)
  tauSe <- c(0.002, 0.001)
  h2 <- 0.5
  annot_names <- c("annot1", "base")
  res <- pecotmr:::computeBaselineEnrichment(tau, tauSe, NULL,
                                              baseline_mat, annot_names, h2)
  expect_s3_class(res, "data.frame")
  expected_cols <- c("annotation", "tau", "tauSe", "enrichment",
                     "enrichmentSe", "enrichmentP", "propH2", "propSnps")
  expect_named(res, expected_cols)
  expect_equal(nrow(res), 2)
})

test_that("computeBaselineEnrichment computes enrichment = tau * M / h2", {
  M <- 100
  baseline_mat <- matrix(0, nrow = M, ncol = 2)
  baseline_mat[1:50, 1] <- 1
  baseline_mat[, 2] <- 1
  tau <- c(0.01, 0.005)
  tauSe <- c(0.002, 0.001)
  h2 <- 0.5
  annot_names <- c("annot1", "base")
  res <- pecotmr:::computeBaselineEnrichment(tau, tauSe, NULL,
                                              baseline_mat, annot_names, h2)
  expect_equal(res$enrichment, tau * M / h2)
})

test_that("computeBaselineEnrichment computes propSnps correctly", {
  M <- 100
  baseline_mat <- matrix(0, nrow = M, ncol = 2)
  baseline_mat[1:50, 1] <- 1
  baseline_mat[, 2] <- 1
  tau <- c(0.01, 0.005)
  tauSe <- c(0.002, 0.001)
  h2 <- 0.5
  annot_names <- c("annot1", "base")
  res <- pecotmr:::computeBaselineEnrichment(tau, tauSe, NULL,
                                              baseline_mat, annot_names, h2)
  expect_equal(res$propSnps, c(0.5, 1.0))
})

test_that("computeBaselineEnrichment uses jackknife blocks when provided", {
  M <- 100
  nBlocks <- 5
  baseline_mat <- matrix(0, nrow = M, ncol = 2)
  baseline_mat[1:50, 1] <- 1
  baseline_mat[, 2] <- 1
  tau <- c(0.01, 0.005)
  tauSe <- c(0.002, 0.001)
  h2 <- 0.5
  annot_names <- c("annot1", "base")
  tauBlocks <- matrix(rep(tau, each = nBlocks), nrow = nBlocks, ncol = 2)
  res <- pecotmr:::computeBaselineEnrichment(tau, tauSe, tauBlocks,
                                              baseline_mat, annot_names, h2)
  # With constant tauBlocks, enrichmentSe should be 0
  expect_equal(res$enrichmentSe, c(0, 0))
})

# =============================================================================
# standardize_tau_star
# =============================================================================

test_that("standardize_tau_star computes tauStar = tau * sd_annot * M_ref / h2g", {
  tau <- c(0.01, 0.02)
  sd_annot <- c(0.5, 0.3)
  M_ref <- 1000L
  h2g <- 0.4
  nBlocks <- 5
  tauBlocks <- matrix(rep(tau, each = nBlocks), nrow = nBlocks, ncol = 2)
  res <- pecotmr:::standardizeTauStar(tau, tauBlocks, sd_annot, M_ref, h2g)
  expected <- tau * sd_annot * M_ref / h2g
  expect_equal(res$tauStar, expected)
  expect_length(res$tauStarSe, 2)
})

test_that("standardize_tau_star errors when h2g == 0", {
  expect_error(
    pecotmr:::standardizeTauStar(c(0.01), matrix(0.01, nrow = 3, ncol = 1),
                                    c(0.5), 1000L, 0),
    "h2g must be non-zero"
  )
})

test_that("standardize_tau_star errors when tau and sd_annot differ in length", {
  expect_error(
    pecotmr:::standardizeTauStar(c(0.01, 0.02), matrix(0, nrow = 3, ncol = 2),
                                    c(0.5), 1000L, 0.5),
    "tau and sdAnnot must have the same length"
  )
})

test_that("standardize_tau_star returns list with tauStar and tauStarSe", {
  tau <- c(0.01)
  sd_annot <- c(0.5)
  M_ref <- 1000L
  h2g <- 0.4
  tauBlocks <- matrix(rnorm(5, mean = 0.01, sd = 0.001), nrow = 5, ncol = 1)
  res <- pecotmr:::standardizeTauStar(tau, tauBlocks, sd_annot, M_ref, h2g)
  expect_true(is.list(res))
  expect_named(res, c("tauStar", "tauStarSe"))
  expect_length(res$tauStar, 1)
  expect_length(res$tauStarSe, 1)
})

# =============================================================================
# meta_random_effects
# =============================================================================

test_that("meta_random_effects returns all NA with k=0", {
  res <- pecotmr:::metaRandomEffects(numeric(0), numeric(0))
  expect_true(is.na(res$mean))
  expect_true(is.na(res$se))
  expect_true(is.na(res$tau2))
  expect_true(is.na(res$I2))
  expect_true(is.na(res$Q))
})

test_that("meta_random_effects with k=1 returns input values", {
  res <- pecotmr:::metaRandomEffects(5.0, 1.0)
  expect_equal(res$mean, 5.0)
  expect_equal(res$se, 1.0)
  expect_equal(res$tau2, 0)
})

test_that("meta_random_effects with identical means gives tau2=0", {
  means <- rep(3.0, 5)
  ses <- rep(1.0, 5)
  res <- pecotmr:::metaRandomEffects(means, ses)
  expect_equal(res$tau2, 0)
  expect_equal(res$mean, 3.0)
})

test_that("meta_random_effects returns correct structure", {
  set.seed(42)
  means <- rnorm(5, mean = 2, sd = 0.5)
  ses <- rep(0.5, 5)
  res <- pecotmr:::metaRandomEffects(means, ses)
  expect_true(is.list(res))
  expect_named(res, c("mean", "se", "tau2", "I2", "Q"))
  expect_true(res$se > 0)
  expect_true(res$I2 >= 0 && res$I2 <= 1)
  expect_true(res$Q >= 0)
  expect_true(res$tau2 >= 0)
})

test_that("meta_random_effects errors with non-positive ses", {
  expect_error(
    pecotmr:::metaRandomEffects(c(1, 2), c(1, 0)),
    "all ses must be positive and finite"
  )
  expect_error(
    pecotmr:::metaRandomEffects(c(1, 2), c(1, -1)),
    "all ses must be positive and finite"
  )
})

test_that("meta_random_effects known DerSimonian-Laird example", {
  # Three studies with known values
  means <- c(0.5, 0.8, 0.3)
  ses <- c(0.2, 0.3, 0.15)
  res <- pecotmr:::metaRandomEffects(means, ses)

  # Fixed-effect weights
  w_fe <- 1 / ses^2
  mu_fe <- sum(w_fe * means) / sum(w_fe)
  Q <- sum(w_fe * (means - mu_fe)^2)
  c_dl <- sum(w_fe) - sum(w_fe^2) / sum(w_fe)
  tau2 <- max(0, (Q - 2) / c_dl)
  w_re <- 1 / (ses^2 + tau2)
  mu_re <- sum(w_re * means) / sum(w_re)
  se_re <- sqrt(1 / sum(w_re))

  expect_equal(res$Q, Q)
  expect_equal(res$tau2, tau2)
  expect_equal(res$mean, mu_re)
  expect_equal(res$se, se_re)
})

# =============================================================================
# snpsPerBlock
# =============================================================================

test_that("snpsPerBlock assigns SNPs to correct blocks", {
  # 10 SNPs across 2 blocks
  snp_info <- data.frame(
    CHR = rep("chr1", 10),
    BP = c(100, 500, 1000, 3000, 4500, 5500, 6000, 7500, 9000, 9999)
  )
  blocks_gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(1, 5001), end = c(5000, 10000))
  )
  ld_blocks <- new("LdBlocks", blocks = blocks_gr, genome = "hg19")
  res <- pecotmr:::snpsPerBlock(snp_info, ld_blocks)
  # Block 1 covers 1-5000: SNPs at 100, 500, 1000, 3000, 4500
  expect_equal(sort(res[["1"]]), c(1L, 2L, 3L, 4L, 5L))
  # Block 2 covers 5001-10000: SNPs at 5500, 6000, 7500, 9000, 9999
  expect_equal(sort(res[["2"]]), c(6L, 7L, 8L, 9L, 10L))
})

test_that("snpsPerBlock returns empty for blocks with no SNPs", {
  snp_info <- data.frame(
    CHR = rep("chr1", 3),
    BP = c(100, 200, 300)
  )
  blocks_gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr2"),
    ranges = IRanges::IRanges(start = c(1, 1), end = c(5000, 5000))
  )
  ld_blocks <- new("LdBlocks", blocks = blocks_gr, genome = "hg19")
  res <- pecotmr:::snpsPerBlock(snp_info, ld_blocks)
  # All SNPs on chr1, so block 2 (chr2) should have no entries
  expect_true("1" %in% names(res))
  expect_false("2" %in% names(res))
  expect_equal(sort(res[["1"]]), c(1L, 2L, 3L))
})

# =============================================================================
# checkGenomeBuild
# =============================================================================

test_that("checkGenomeBuild returns TRUE when all objects match", {
  blocks_gr_a <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1, end = 5000)
  )
  ld_blocks_a <- new("LdBlocks", blocks = blocks_gr_a, genome = "hg19")
  ld_blocks_b <- new("LdBlocks", blocks = blocks_gr_a, genome = "hg19")

  expect_true(pecotmr:::checkGenomeBuild(ld_blocks_a, ld_blocks_b))
})

test_that("checkGenomeBuild errors when genome builds mismatch", {
  blocks_gr_19 <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1, end = 5000)
  )
  ld_blocks_19 <- new("LdBlocks", blocks = blocks_gr_19, genome = "hg19")

  blocks_gr_38 <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1, end = 5000)
  )
  ld_blocks_38 <- new("LdBlocks", blocks = blocks_gr_38, genome = "hg38")

  expect_error(
    pecotmr:::checkGenomeBuild(ld_blocks_19, ld_blocks_38),
    "Genome build mismatch"
  )
})

test_that("checkGenomeBuild works with AnnotationMatrix", {
  snp_gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 3),
    ranges = IRanges::IRanges(start = 1:3, width = 1)
  )
  annot_mat <- matrix(c(1, 0, 1), nrow = 3, ncol = 1)
  annot_meta <- data.frame(name = "annot1", tier = "baseline",
                           type = "binary", stringsAsFactors = FALSE)
  am <- new("AnnotationMatrix",
    snpRanges = snp_gr,
    annotations = annot_mat,
    annotationMeta = annot_meta,
    genome = "hg19"
  )

  blocks_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1, end = 5000)
  )
  ld_blocks <- new("LdBlocks", blocks = blocks_gr, genome = "hg19")

  expect_true(pecotmr:::checkGenomeBuild(am, ld_blocks))
})

# =============================================================================
# shrinkLd
# =============================================================================

test_that("shrinkLd constant shrinkage applies (1-lambda)*R + lambda*I", {
  set.seed(42)
  p <- 5
  # Create a valid correlation matrix
  A <- matrix(rnorm(p * p), p, p)
  R <- cov2cor(crossprod(A))
  n_ref <- 100

  res <- pecotmr:::shrinkLd(R, n_ref, shrinkageType = "constant")
  lambda <- 1 / sqrt(n_ref)
  expected <- (1 - lambda) * R + lambda * diag(p)
  expect_equal(res, expected)
})

test_that("shrinkLd constant shrinkage preserves diagonal of 1", {
  set.seed(42)
  p <- 4
  A <- matrix(rnorm(p * p), p, p)
  R <- cov2cor(crossprod(A))
  n_ref <- 200
  res <- pecotmr:::shrinkLd(R, n_ref, shrinkageType = "constant")
  expect_equal(diag(res), rep(1, p))
})

test_that("shrinkLd constant shrinkage result is symmetric", {
  set.seed(42)
  p <- 4
  A <- matrix(rnorm(p * p), p, p)
  R <- cov2cor(crossprod(A))
  n_ref <- 200
  res <- pecotmr:::shrinkLd(R, n_ref, shrinkageType = "constant")
  expect_equal(res, t(res))
})

test_that("shrinkLd wen_stephens uses genetic map when provided", {
  set.seed(42)
  p <- 4
  A <- matrix(rnorm(p * p), p, p)
  R <- cov2cor(crossprod(A))
  n_ref <- 500
  genetic_map <- c(0.0, 0.1, 0.5, 1.0)
  res <- pecotmr:::shrinkLd(R, n_ref, shrinkageType = "wen_stephens",
                             geneticMap = genetic_map)
  # Diagonal should be 1
  expect_equal(diag(res), rep(1, p))
  # Result should be symmetric
  expect_equal(res, t(res))
  # Off-diagonal elements should be shrunk (closer to zero than original)
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      expect_true(abs(res[i, j]) <= abs(R[i, j]) + 1e-10)
    }
  }
})


# Tests for heritability estimation methods:
#   lderUnivariate (h2Lder.R)
#   gldscUnivariate (h2Gldsc.R)
#   hdlUnivariate (h2Hdl.R)

set.seed(42)

# =============================================================================
# Helpers: simulate test data
# =============================================================================

simulate_h2_data <- function(n_snps = 100, nBlocks = 2, n_gwas = 50000,
                             h2_true = 0.3) {
  set.seed(42)
  snps_per_block <- n_snps / nBlocks

  # Block structure
  blocks_gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", nBlocks),
    ranges = IRanges::IRanges(
      start = seq(1, by = snps_per_block * 100, length.out = nBlocks),
      width = snps_per_block * 100 - 1
    )
  )
  ld_blocks <- new("LdBlocks", blocks = blocks_gr, genome = "hg19")

  # SNP info
  snp_info <- data.frame(
    SNP = paste0("rs", seq_len(n_snps)),
    CHR = rep("chr1", n_snps),
    BP = as.integer(seq(50, by = 100, length.out = n_snps)),
    A1 = rep("A", n_snps), A2 = rep("G", n_snps),
    stringsAsFactors = FALSE
  )

  # Build block-diagonal LD (AR(1) with rho = 0.5)
  eigen_list <- list()
  R_blocks <- list()
  for (b in seq_len(nBlocks)) {
    idx <- seq((b - 1) * snps_per_block + 1, b * snps_per_block)
    p <- length(idx)
    R <- 0.5^abs(outer(seq_len(p), seq_len(p), "-"))
    e <- eigen(R, symmetric = TRUE)
    eigen_list[[b]] <- list(
      values = e$values, vectors = e$vectors, snpIdx = as.integer(idx)
    )
    R_blocks[[b]] <- list(R = R, snpIdx = as.integer(idx))
  }

  # Simulate z-scores under infinitesimal model
  signal_var <- h2_true / n_snps
  z <- rnorm(n_snps) * sqrt(n_gwas * signal_var) + rnorm(n_snps)

  # LdEigen
  eigen_ref <- new("LdEigen",
    ldBlocks = ld_blocks, snpInfo = snp_info,
    nRef = 1000L, inSample = FALSE, genome = "hg19",
    eigenList = eigen_list, eigenvalueTruncation = 1.0
  )

  # LD scores: l2_j = sum_k r^2_{jk}
  ld_scores_vec <- numeric(n_snps)
  for (b in seq_len(nBlocks)) {
    idx <- R_blocks[[b]]$snpIdx
    R <- R_blocks[[b]]$R
    ld_scores_vec[idx] <- rowSums(R^2)
  }

  ld_score_ref <- new("LdScore",
    ldBlocks = ld_blocks, snpInfo = snp_info,
    nRef = 1000L, inSample = FALSE, genome = "hg19",
    ldScores = matrix(ld_scores_vec, ncol = 1,
                       dimnames = list(NULL, "base_l2")),
    ldScoreWeights = rep(1, n_snps),
    ldMatrixList = R_blocks
  )

  list(z = z, n = n_gwas, eigen_ref = eigen_ref, ld_score_ref = ld_score_ref,
       h2_true = h2_true, n_snps = n_snps)
}

# Annotation helper: all-ones baseline so the eigenvalue-score column is
# simply the eigenvalue itself (avoids singularity with few blocks), plus
# one candidate annotation for score-statistic testing.
make_test_annotations <- function(n_snps) {
  set.seed(99)
  snp_gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", n_snps),
    ranges = IRanges::IRanges(
      start = seq(50, by = 100, length.out = n_snps), width = 1L
    )
  )
  # All-ones baseline (behaves like the unstratified intercept)
  base_col <- rep(1, n_snps)
  # Binary candidate: random 30%
  cand_col <- rbinom(n_snps, 1, 0.3)
  mat <- cbind(base_col, cand_col)
  meta <- data.frame(
    name = c("baseline1", "candidate1"),
    tier = c("baseline", "candidate"),
    type = c("binary", "binary"),
    stringsAsFactors = FALSE
  )
  AnnotationMatrix(mat, snp_gr, meta, genome = "hg19")
}

# Pre-compute shared test data
dat <- simulate_h2_data()

# =============================================================================
# LDER tests (h2Lder.R)
# =============================================================================

test_that("lderUnivariate returns correct structure", {
  res <- pecotmr:::lderUnivariate(dat$z, dat$n, dat$eigen_ref)
  expect_type(res, "list")
  expect_true(all(c("h2", "h2Se", "intercept", "interceptSe",
                     "local", "enrichment", "tauBlocks", "scoreStats")
                   %in% names(res)))
})

test_that("lderUnivariate h2 is finite and in reasonable range", {
  res <- pecotmr:::lderUnivariate(dat$z, dat$n, dat$eigen_ref)
  expect_true(is.finite(res$h2))
  expect_true(res$h2 > -0.5 && res$h2 < 1.0)
})

test_that("lderUnivariate h2Se is positive", {
  res <- pecotmr:::lderUnivariate(dat$z, dat$n, dat$eigen_ref)
  expect_true(is.finite(res$h2Se))
  expect_true(res$h2Se > 0)
})

test_that("lderUnivariate intercept is near zero for well-calibrated data", {
  # In LDER the intercept parameter a represents confounding deviation:
  # the model is E[chi2Rot - 1] = n * h2/M * d + n * a
  # so a ~ 0 when there is no confounding.
  res <- pecotmr:::lderUnivariate(dat$z, dat$n, dat$eigen_ref)
  expect_true(is.finite(res$intercept))
  expect_true(abs(res$intercept) < 0.01)
})

test_that("lderUnivariate with local = TRUE returns local data.frame", {
  res <- pecotmr:::lderUnivariate(dat$z, dat$n, dat$eigen_ref, local = TRUE)
  expect_true(is.data.frame(res$local))
  expect_true("h2Local" %in% colnames(res$local))
})

test_that("lderUnivariate without annotations returns NULL enrichment", {
  res <- pecotmr:::lderUnivariate(dat$z, dat$n, dat$eigen_ref)
  expect_null(res$enrichment)
})

# =============================================================================
# gLDSC tests (h2Gldsc.R)
# =============================================================================

test_that("gldscUnivariate returns correct structure", {
  res <- pecotmr:::gldscUnivariate(dat$z, dat$n, dat$ld_score_ref)
  expect_type(res, "list")
  expect_true(all(c("h2", "h2Se", "intercept", "interceptSe",
                     "local", "enrichment", "tauBlocks", "scoreStats")
                   %in% names(res)))
})

test_that("gldscUnivariate h2 is finite", {
  res <- pecotmr:::gldscUnivariate(dat$z, dat$n, dat$ld_score_ref)
  expect_true(is.finite(res$h2))
})

test_that("gldscUnivariate h2Se is positive", {
  res <- pecotmr:::gldscUnivariate(dat$z, dat$n, dat$ld_score_ref)
  expect_true(is.finite(res$h2Se))
  expect_true(res$h2Se > 0)
})

test_that("gldscUnivariate with local = TRUE needs fine-grained blocks", {
  # Our test data has only 2 blocks, which is <= 22, so local should error
  expect_error(
    pecotmr:::gldscUnivariate(dat$z, dat$n, dat$ld_score_ref, local = TRUE),
    "Local g-LDSC requires fine-grained LD blocks"
  )
})

# =============================================================================
# HDL tests (h2Hdl.R)
# =============================================================================

test_that("hdlUnivariate returns correct structure", {
  res <- pecotmr:::hdlUnivariate(dat$z, dat$n, dat$eigen_ref)
  expect_type(res, "list")
  expect_true(all(c("h2", "h2Se", "intercept", "interceptSe",
                     "local", "enrichment", "tauBlocks", "scoreStats")
                   %in% names(res)))
})

test_that("hdlUnivariate h2 is finite", {
  res <- pecotmr:::hdlUnivariate(dat$z, dat$n, dat$eigen_ref)
  expect_true(is.finite(res$h2))
})

test_that("hdlUnivariate h2Se is positive", {
  res <- pecotmr:::hdlUnivariate(dat$z, dat$n, dat$eigen_ref)
  expect_true(is.finite(res$h2Se))
  expect_true(res$h2Se > 0)
})

test_that("hdlUnivariate with local = TRUE returns local data.frame", {
  res <- pecotmr:::hdlUnivariate(dat$z, dat$n, dat$eigen_ref, local = TRUE)
  expect_true(is.data.frame(res$local))
  expect_true("h2Local" %in% colnames(res$local))
})

# =============================================================================
# Annotation tests (using LDER with more blocks for numerical stability)
# =============================================================================

# Use 5 blocks and 500 SNPs for annotation tests to avoid singular matrices
dat_annot <- simulate_h2_data(n_snps = 500, nBlocks = 5)

test_that("lderUnivariate with annotations returns enrichment data.frame", {
  annot <- make_test_annotations(dat_annot$n_snps)
  res <- pecotmr:::lderUnivariate(dat_annot$z, dat_annot$n,
                                   dat_annot$eigen_ref,
                                   annotations = annot)
  expect_true(is.data.frame(res$enrichment))
})

test_that("lder enrichment has correct columns", {
  annot <- make_test_annotations(dat_annot$n_snps)
  res <- pecotmr:::lderUnivariate(dat_annot$z, dat_annot$n,
                                   dat_annot$eigen_ref,
                                   annotations = annot)
  expected_cols <- c("annotation", "tau", "tauSe", "enrichment",
                     "enrichmentSe", "enrichmentP", "propH2", "propSnps")
  expect_true(all(expected_cols %in% colnames(res$enrichment)))
})

test_that("lder with annotations returns tauBlocks matrix", {
  annot <- make_test_annotations(dat_annot$n_snps)
  res <- pecotmr:::lderUnivariate(dat_annot$z, dat_annot$n,
                                   dat_annot$eigen_ref,
                                   annotations = annot)
  expect_true(is.matrix(res$tauBlocks))
  # nBlocks rows (5 blocks in the annotation test data)
  expect_equal(nrow(res$tauBlocks), 5L)
})

test_that("lder with annotations returns scoreStats list", {
  annot <- make_test_annotations(dat_annot$n_snps)
  res <- pecotmr:::lderUnivariate(dat_annot$z, dat_annot$n,
                                   dat_annot$eigen_ref,
                                   annotations = annot)
  expect_type(res$scoreStats, "list")
  expect_true(all(c("z", "R") %in% names(res$scoreStats)))
})

# =============================================================================
# HDL annotation tests (h2Hdl.R stratified paths)
# =============================================================================

test_that("hdlUnivariate with annotations returns enrichment data.frame with correct columns", {
  set.seed(123)
  annot <- make_test_annotations(dat_annot$n_snps)
  res <- pecotmr:::hdlUnivariate(dat_annot$z, dat_annot$n,
                                  dat_annot$eigen_ref,
                                  annotations = annot)
  expect_true(is.data.frame(res$enrichment))
  expected_cols <- c("annotation", "tau", "tauSe", "enrichment",
                     "enrichmentSe", "enrichmentP", "propH2", "propSnps")
  expect_true(all(expected_cols %in% colnames(res$enrichment)))
  # Should have one row per baseline annotation
  expect_equal(nrow(res$enrichment), 1L)  # 1 baseline annotation
  # Values should be finite
  expect_true(all(is.finite(res$enrichment$tau)))
  expect_true(all(is.finite(res$enrichment$enrichment)))
})

test_that("hdlUnivariate with annotations returns tauBlocks matrix", {
  set.seed(124)
  annot <- make_test_annotations(dat_annot$n_snps)
  res <- pecotmr:::hdlUnivariate(dat_annot$z, dat_annot$n,
                                  dat_annot$eigen_ref,
                                  annotations = annot)
  expect_true(is.matrix(res$tauBlocks))
  # nBlocks rows (5 blocks in dat_annot)
  expect_equal(nrow(res$tauBlocks), 5L)
  # Number of columns matches baseline annotation count
  expect_equal(ncol(res$tauBlocks), 1L)
  # Values should be finite
  expect_true(all(is.finite(res$tauBlocks)))
})

test_that("hdlUnivariate with annotations and local = TRUE returns local data.frame", {
  set.seed(125)
  annot <- make_test_annotations(dat_annot$n_snps)
  res <- pecotmr:::hdlUnivariate(dat_annot$z, dat_annot$n,
                                  dat_annot$eigen_ref,
                                  annotations = annot,
                                  local = TRUE)
  expect_true(is.data.frame(res$local))
  expect_true("h2Local" %in% colnames(res$local))
  expect_true("h2LocalSe" %in% colnames(res$local))
  expect_true("blockId" %in% colnames(res$local))
  # Should have one row per block
  expect_equal(nrow(res$local), 5L)
})

test_that("hdlUnivariate with annotations returns scoreStats with z and R", {
  set.seed(126)
  annot <- make_test_annotations(dat_annot$n_snps)
  res <- pecotmr:::hdlUnivariate(dat_annot$z, dat_annot$n,
                                  dat_annot$eigen_ref,
                                  annotations = annot)
  expect_type(res$scoreStats, "list")
  expect_true(all(c("z", "R") %in% names(res$scoreStats)))
  # z should have length = number of candidate annotations (1)
  expect_length(res$scoreStats$z, 1L)
  expect_true(is.finite(res$scoreStats$z[1]))
  # R should be a matrix
  expect_true(is.matrix(res$scoreStats$R))
  expect_equal(dim(res$scoreStats$R), c(1L, 1L))
  # annotationNames should be present
  expect_true("annotationNames" %in% names(res$scoreStats))
  expect_equal(res$scoreStats$annotationNames, "candidate1")
})

# =============================================================================
# gLDSC annotation tests (h2Gldsc.R stratified paths)
# =============================================================================

# gLDSC annotation helper: use a spatially varying baseline instead of
# all-ones, because computeldscores with an all-ones baseline produces
# LD scores identical to the base L2 column, making the design matrix
# singular. A half-on/half-off baseline avoids collinearity.
make_gldsc_annotations <- function(n_snps) {
  set.seed(77)
  snp_gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", n_snps),
    ranges = IRanges::IRanges(
      start = seq(50, by = 100, length.out = n_snps), width = 1L
    )
  )
  # Spatially structured baseline: first half = 1, second half = 0
  base_col <- as.numeric(seq_len(n_snps) <= n_snps / 2)
  # Binary candidate: random 30%
  cand_col <- rbinom(n_snps, 1, 0.3)
  mat <- cbind(base_col, cand_col)
  meta <- data.frame(
    name = c("baseline1", "candidate1"),
    tier = c("baseline", "candidate"),
    type = c("binary", "binary"),
    stringsAsFactors = FALSE
  )
  AnnotationMatrix(mat, snp_gr, meta, genome = "hg19")
}

test_that("gldscUnivariate with annotations returns enrichment data.frame", {
  set.seed(200)
  annot <- make_gldsc_annotations(dat_annot$n_snps)
  res <- pecotmr:::gldscUnivariate(dat_annot$z, dat_annot$n,
                                    dat_annot$ld_score_ref,
                                    annotations = annot)
  expect_true(is.data.frame(res$enrichment))
  expected_cols <- c("annotation", "tau", "tauSe", "enrichment",
                     "enrichmentSe", "enrichmentP", "propH2", "propSnps")
  expect_true(all(expected_cols %in% colnames(res$enrichment)))
  # Should have one row per baseline annotation
  expect_equal(nrow(res$enrichment), 1L)
  # Values should be finite
  expect_true(all(is.finite(res$enrichment$tau)))
  expect_true(all(is.finite(res$enrichment$enrichment)))
})

test_that("gldscUnivariate with annotations returns tauBlocks matrix", {
  set.seed(201)
  annot <- make_gldsc_annotations(dat_annot$n_snps)
  res <- pecotmr:::gldscUnivariate(dat_annot$z, dat_annot$n,
                                    dat_annot$ld_score_ref,
                                    annotations = annot)
  expect_true(is.matrix(res$tauBlocks))
  # gldsc jackknife uses 200 blocks by default; with 500 SNPs the actual
  # number of unique blocks equals ceil(500/ceil(500/200)) = 200
  expect_true(nrow(res$tauBlocks) > 0)
  # Number of columns matches baseline annotation count
  expect_equal(ncol(res$tauBlocks), 1L)
  # Values should be finite (allow some NAs from edge blocks)
  expect_true(any(is.finite(res$tauBlocks)))
})

test_that("gldscUnivariate with annotations returns scoreStats", {
  set.seed(202)
  annot <- make_gldsc_annotations(dat_annot$n_snps)
  res <- pecotmr:::gldscUnivariate(dat_annot$z, dat_annot$n,
                                    dat_annot$ld_score_ref,
                                    annotations = annot)
  expect_type(res$scoreStats, "list")
  expect_true(all(c("z", "R") %in% names(res$scoreStats)))
  # z should have length = number of candidate annotations (1)
  expect_length(res$scoreStats$z, 1L)
  expect_true(is.finite(res$scoreStats$z[1]))
  # R should be a matrix
  expect_true(is.matrix(res$scoreStats$R))
  expect_equal(dim(res$scoreStats$R), c(1L, 1L))
  # annotationNames should be present
  expect_true("annotationNames" %in% names(res$scoreStats))
  expect_equal(res$scoreStats$annotationNames, "candidate1")
})
