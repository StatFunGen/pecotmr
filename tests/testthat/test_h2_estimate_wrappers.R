context("h2_estimate_wrappers")

# ===========================================================================
# Helper functions to build test S4 objects
# ===========================================================================

make_test_eigen_ref <- function(nSnps = 20, n_blocks = 2) {
  snps_per_block <- nSnps / n_blocks
  blocks_gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", n_blocks),
    ranges = IRanges::IRanges(
      start = seq(1, by = snps_per_block * 100, length.out = n_blocks),
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

  eigen_list <- lapply(seq_len(n_blocks), function(b) {
    idx <- seq((b - 1) * snps_per_block + 1, b * snps_per_block)
    p <- length(idx)
    set.seed(42 + b)
    R <- matrix(0.3, p, p)
    diag(R) <- 1
    e <- eigen(R)
    list(values = e$values, vectors = e$vectors, snp_idx = as.integer(idx))
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

make_test_score_ref <- function(nSnps = 20, n_blocks = 2,
                                with_ld_matrices = FALSE) {
  snps_per_block <- nSnps / n_blocks
  blocks_gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", n_blocks),
    ranges = IRanges::IRanges(
      start = seq(1, by = snps_per_block * 100, length.out = n_blocks),
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
    lapply(seq_len(n_blocks), function(b) {
      idx <- seq((b - 1) * snps_per_block + 1, b * snps_per_block)
      p <- length(idx)
      R <- matrix(0.3, p, p)
      diag(R) <- 1
      list(R = R, snp_idx = as.integer(idx))
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
      tau_se = c(5e-8, 6e-8),
      enrichment = c(2.0, 3.0),
      enrichment_se = c(0.5, 0.7),
      enrichment_p = c(0.01, 0.001),
      prop_h2 = c(0.3, 0.5),
      prop_snps = c(0.15, 0.17),
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }

  set.seed(55)
  tau_blocks <- if (with_enrichment) {
    matrix(rnorm(20), nrow = 10, ncol = 2,
           dimnames = list(NULL, c("annot1", "annot2")))
  } else {
    NULL
  }

  new("H2Estimate",
    h2 = 0.3, h2Se = 0.05,
    intercept = 1.01, interceptSe = 0.02,
    local = NULL, enrichment = enrich,
    tauBlocks = tau_blocks, scoreStats = NULL,
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
    block_id = 1:3,
    h2_local = c(0.1, 0.15, 0.05),
    h2_local_se = c(0.02, 0.03, 0.01),
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

test_that("getscorestats returns the score_stats slot", {
  score_stats <- list(z = c(1.5, -0.8), R = diag(2),
                      annotation_names = c("a1", "a2"))
  h2_obj <- new("H2Estimate",
    h2 = 0.3, h2Se = 0.05,
    intercept = 1.0, interceptSe = 0.01,
    local = NULL, enrichment = NULL,
    tauBlocks = NULL, scoreStats = score_stats,
    method = "lder", nSnps = 100L, traitName = "test"
  )
  expect_identical(getScoreStats(h2_obj), score_stats)
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

  expected_names <- c("categories", "tau", "tau_se", "enrichment",
                      "enrichment_se", "enrichment_p", "prop_h2",
                      "prop_snps", "h2g", "tau_blocks", "n_blocks")
  expect_true(all(expected_names %in% names(result)))
  expect_equal(length(result), length(expected_names))

  expect_equal(result$categories, c("annot1", "annot2"))
  expect_equal(result$h2g, 0.3)
  expect_true(is.matrix(result$tau_blocks))
  expect_equal(result$n_blocks, 10L)
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

test_that("h2EstimateToSldscTrait creates 1-row dummy when tau_blocks is NULL", {
  h2_obj <- make_test_h2estimate(with_enrichment = TRUE)
  h2_obj@tauBlocks <- NULL
  result <- h2EstimateToSldscTrait(h2_obj)
  expect_equal(nrow(result$tau_blocks), 1)
  expect_equal(result$n_blocks, 1L)
  expect_equal(colnames(result$tau_blocks), c("annot1", "annot2"))
})

# ===========================================================================
# Helper: GwasSumStats matched to a reference panel
# ===========================================================================

make_test_sumstats_for_ref <- function(ref, traitName = "test") {
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
  GwasSumStats(df, traitName = traitName, genome = "hg19")
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
  n_snps <- nrow(eigen_ref@snpInfo)
  set.seed(123)
  df <- data.frame(
    SNP = eigen_ref@snpInfo$SNP,
    CHR = sub("^chr", "", eigen_ref@snpInfo$CHR),
    BP = eigen_ref@snpInfo$BP,
    A1 = eigen_ref@snpInfo$A1,
    A2 = eigen_ref@snpInfo$A2,
    Z = rnorm(n_snps),
    N = rep(50000, n_snps),
    stringsAsFactors = FALSE
  )
  ss <- GwasSumStats(df, traitName = "cc_trait", genome = "hg19", varY = 4.0)
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
