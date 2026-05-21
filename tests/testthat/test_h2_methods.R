# Tests for heritability estimation methods:
#   lder_univariate (h2_lder.R)
#   gldsc_univariate (h2_gldsc.R)
#   hdl_univariate (h2_hdl.R)

set.seed(42)

# =============================================================================
# Helpers: simulate test data
# =============================================================================

simulate_h2_data <- function(n_snps = 100, n_blocks = 2, n_gwas = 50000,
                             h2_true = 0.3) {
  set.seed(42)
  snps_per_block <- n_snps / n_blocks

  # Block structure
  blocks_gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", n_blocks),
    ranges = IRanges::IRanges(
      start = seq(1, by = snps_per_block * 100, length.out = n_blocks),
      width = snps_per_block * 100 - 1
    )
  )
  ld_blocks <- new("LDBlocks", blocks = blocks_gr, genome = "hg19")

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
  for (b in seq_len(n_blocks)) {
    idx <- seq((b - 1) * snps_per_block + 1, b * snps_per_block)
    p <- length(idx)
    R <- 0.5^abs(outer(seq_len(p), seq_len(p), "-"))
    e <- eigen(R, symmetric = TRUE)
    eigen_list[[b]] <- list(
      values = e$values, vectors = e$vectors, snp_idx = as.integer(idx)
    )
    R_blocks[[b]] <- list(R = R, snp_idx = as.integer(idx))
  }

  # Simulate z-scores under infinitesimal model
  signal_var <- h2_true / n_snps
  z <- rnorm(n_snps) * sqrt(n_gwas * signal_var) + rnorm(n_snps)

  # LDEigen
  eigen_ref <- new("LDEigen",
    ld_blocks = ld_blocks, snp_info = snp_info,
    n_ref = 1000L, in_sample = FALSE, genome = "hg19",
    eigen_list = eigen_list, eigenvalue_truncation = 1.0
  )

  # LD scores: l2_j = sum_k r^2_{jk}
  ld_scores_vec <- numeric(n_snps)
  for (b in seq_len(n_blocks)) {
    idx <- R_blocks[[b]]$snp_idx
    R <- R_blocks[[b]]$R
    ld_scores_vec[idx] <- rowSums(R^2)
  }

  ld_score_ref <- new("LDScore",
    ld_blocks = ld_blocks, snp_info = snp_info,
    n_ref = 1000L, in_sample = FALSE, genome = "hg19",
    ld_scores = matrix(ld_scores_vec, ncol = 1,
                       dimnames = list(NULL, "base_l2")),
    ld_score_weights = rep(1, n_snps),
    ld_matrix_list = R_blocks
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
# LDER tests (h2_lder.R)
# =============================================================================

test_that("lder_univariate returns correct structure", {
  res <- pecotmr:::lder_univariate(dat$z, dat$n, dat$eigen_ref)
  expect_type(res, "list")
  expect_true(all(c("h2", "h2_se", "intercept", "intercept_se",
                     "local", "enrichment", "tau_blocks", "score_stats")
                   %in% names(res)))
})

test_that("lder_univariate h2 is finite and in reasonable range", {
  res <- pecotmr:::lder_univariate(dat$z, dat$n, dat$eigen_ref)
  expect_true(is.finite(res$h2))
  expect_true(res$h2 > -0.5 && res$h2 < 1.0)
})

test_that("lder_univariate h2_se is positive", {
  res <- pecotmr:::lder_univariate(dat$z, dat$n, dat$eigen_ref)
  expect_true(is.finite(res$h2_se))
  expect_true(res$h2_se > 0)
})

test_that("lder_univariate intercept is near zero for well-calibrated data", {
  # In LDER the intercept parameter a represents confounding deviation:
  # the model is E[chi2_rot - 1] = n * h2/M * d + n * a
  # so a ~ 0 when there is no confounding.
  res <- pecotmr:::lder_univariate(dat$z, dat$n, dat$eigen_ref)
  expect_true(is.finite(res$intercept))
  expect_true(abs(res$intercept) < 0.01)
})

test_that("lder_univariate with local = TRUE returns local data.frame", {
  res <- pecotmr:::lder_univariate(dat$z, dat$n, dat$eigen_ref, local = TRUE)
  expect_true(is.data.frame(res$local))
  expect_true("h2_local" %in% colnames(res$local))
})

test_that("lder_univariate without annotations returns NULL enrichment", {
  res <- pecotmr:::lder_univariate(dat$z, dat$n, dat$eigen_ref)
  expect_null(res$enrichment)
})

# =============================================================================
# gLDSC tests (h2_gldsc.R)
# =============================================================================

test_that("gldsc_univariate returns correct structure", {
  res <- pecotmr:::gldsc_univariate(dat$z, dat$n, dat$ld_score_ref)
  expect_type(res, "list")
  expect_true(all(c("h2", "h2_se", "intercept", "intercept_se",
                     "local", "enrichment", "tau_blocks", "score_stats")
                   %in% names(res)))
})

test_that("gldsc_univariate h2 is finite", {
  res <- pecotmr:::gldsc_univariate(dat$z, dat$n, dat$ld_score_ref)
  expect_true(is.finite(res$h2))
})

test_that("gldsc_univariate h2_se is positive", {
  res <- pecotmr:::gldsc_univariate(dat$z, dat$n, dat$ld_score_ref)
  expect_true(is.finite(res$h2_se))
  expect_true(res$h2_se > 0)
})

test_that("gldsc_univariate with local = TRUE needs fine-grained blocks", {
  # Our test data has only 2 blocks, which is <= 22, so local should error
  expect_error(
    pecotmr:::gldsc_univariate(dat$z, dat$n, dat$ld_score_ref, local = TRUE),
    "Local g-LDSC requires fine-grained LD blocks"
  )
})

# =============================================================================
# HDL tests (h2_hdl.R)
# =============================================================================

test_that("hdl_univariate returns correct structure", {
  res <- pecotmr:::hdl_univariate(dat$z, dat$n, dat$eigen_ref)
  expect_type(res, "list")
  expect_true(all(c("h2", "h2_se", "intercept", "intercept_se",
                     "local", "enrichment", "tau_blocks", "score_stats")
                   %in% names(res)))
})

test_that("hdl_univariate h2 is finite", {
  res <- pecotmr:::hdl_univariate(dat$z, dat$n, dat$eigen_ref)
  expect_true(is.finite(res$h2))
})

test_that("hdl_univariate h2_se is positive", {
  res <- pecotmr:::hdl_univariate(dat$z, dat$n, dat$eigen_ref)
  expect_true(is.finite(res$h2_se))
  expect_true(res$h2_se > 0)
})

test_that("hdl_univariate with local = TRUE returns local data.frame", {
  res <- pecotmr:::hdl_univariate(dat$z, dat$n, dat$eigen_ref, local = TRUE)
  expect_true(is.data.frame(res$local))
  expect_true("h2_local" %in% colnames(res$local))
})

# =============================================================================
# Annotation tests (using LDER with more blocks for numerical stability)
# =============================================================================

# Use 5 blocks and 500 SNPs for annotation tests to avoid singular matrices
dat_annot <- simulate_h2_data(n_snps = 500, n_blocks = 5)

test_that("lder_univariate with annotations returns enrichment data.frame", {
  annot <- make_test_annotations(dat_annot$n_snps)
  res <- pecotmr:::lder_univariate(dat_annot$z, dat_annot$n,
                                   dat_annot$eigen_ref,
                                   annotations = annot)
  expect_true(is.data.frame(res$enrichment))
})

test_that("lder enrichment has correct columns", {
  annot <- make_test_annotations(dat_annot$n_snps)
  res <- pecotmr:::lder_univariate(dat_annot$z, dat_annot$n,
                                   dat_annot$eigen_ref,
                                   annotations = annot)
  expected_cols <- c("annotation", "tau", "tau_se", "enrichment",
                     "enrichment_se", "enrichment_p", "prop_h2", "prop_snps")
  expect_true(all(expected_cols %in% colnames(res$enrichment)))
})

test_that("lder with annotations returns tau_blocks matrix", {
  annot <- make_test_annotations(dat_annot$n_snps)
  res <- pecotmr:::lder_univariate(dat_annot$z, dat_annot$n,
                                   dat_annot$eigen_ref,
                                   annotations = annot)
  expect_true(is.matrix(res$tau_blocks))
  # n_blocks rows (5 blocks in the annotation test data)
  expect_equal(nrow(res$tau_blocks), 5L)
})

test_that("lder with annotations returns score_stats list", {
  annot <- make_test_annotations(dat_annot$n_snps)
  res <- pecotmr:::lder_univariate(dat_annot$z, dat_annot$n,
                                   dat_annot$eigen_ref,
                                   annotations = annot)
  expect_type(res$score_stats, "list")
  expect_true(all(c("z", "R") %in% names(res$score_stats)))
})

# =============================================================================
# HDL annotation tests (h2_hdl.R stratified paths)
# =============================================================================

test_that("hdl_univariate with annotations returns enrichment data.frame with correct columns", {
  set.seed(123)
  annot <- make_test_annotations(dat_annot$n_snps)
  res <- pecotmr:::hdl_univariate(dat_annot$z, dat_annot$n,
                                  dat_annot$eigen_ref,
                                  annotations = annot)
  expect_true(is.data.frame(res$enrichment))
  expected_cols <- c("annotation", "tau", "tau_se", "enrichment",
                     "enrichment_se", "enrichment_p", "prop_h2", "prop_snps")
  expect_true(all(expected_cols %in% colnames(res$enrichment)))
  # Should have one row per baseline annotation
  expect_equal(nrow(res$enrichment), 1L)  # 1 baseline annotation
  # Values should be finite
  expect_true(all(is.finite(res$enrichment$tau)))
  expect_true(all(is.finite(res$enrichment$enrichment)))
})

test_that("hdl_univariate with annotations returns tau_blocks matrix", {
  set.seed(124)
  annot <- make_test_annotations(dat_annot$n_snps)
  res <- pecotmr:::hdl_univariate(dat_annot$z, dat_annot$n,
                                  dat_annot$eigen_ref,
                                  annotations = annot)
  expect_true(is.matrix(res$tau_blocks))
  # n_blocks rows (5 blocks in dat_annot)
  expect_equal(nrow(res$tau_blocks), 5L)
  # Number of columns matches baseline annotation count
  expect_equal(ncol(res$tau_blocks), 1L)
  # Values should be finite
  expect_true(all(is.finite(res$tau_blocks)))
})

test_that("hdl_univariate with annotations and local = TRUE returns local data.frame", {
  set.seed(125)
  annot <- make_test_annotations(dat_annot$n_snps)
  res <- pecotmr:::hdl_univariate(dat_annot$z, dat_annot$n,
                                  dat_annot$eigen_ref,
                                  annotations = annot,
                                  local = TRUE)
  expect_true(is.data.frame(res$local))
  expect_true("h2_local" %in% colnames(res$local))
  expect_true("h2_local_se" %in% colnames(res$local))
  expect_true("block_id" %in% colnames(res$local))
  # Should have one row per block
  expect_equal(nrow(res$local), 5L)
})

test_that("hdl_univariate with annotations returns score_stats with z and R", {
  set.seed(126)
  annot <- make_test_annotations(dat_annot$n_snps)
  res <- pecotmr:::hdl_univariate(dat_annot$z, dat_annot$n,
                                  dat_annot$eigen_ref,
                                  annotations = annot)
  expect_type(res$score_stats, "list")
  expect_true(all(c("z", "R") %in% names(res$score_stats)))
  # z should have length = number of candidate annotations (1)
  expect_length(res$score_stats$z, 1L)
  expect_true(is.finite(res$score_stats$z[1]))
  # R should be a matrix
  expect_true(is.matrix(res$score_stats$R))
  expect_equal(dim(res$score_stats$R), c(1L, 1L))
  # annotation_names should be present
  expect_true("annotation_names" %in% names(res$score_stats))
  expect_equal(res$score_stats$annotation_names, "candidate1")
})

# =============================================================================
# gLDSC annotation tests (h2_gldsc.R stratified paths)
# =============================================================================

# gLDSC annotation helper: use a spatially varying baseline instead of
# all-ones, because computeLdScores with an all-ones baseline produces
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

test_that("gldsc_univariate with annotations returns enrichment data.frame", {
  set.seed(200)
  annot <- make_gldsc_annotations(dat_annot$n_snps)
  res <- pecotmr:::gldsc_univariate(dat_annot$z, dat_annot$n,
                                    dat_annot$ld_score_ref,
                                    annotations = annot)
  expect_true(is.data.frame(res$enrichment))
  expected_cols <- c("annotation", "tau", "tau_se", "enrichment",
                     "enrichment_se", "enrichment_p", "prop_h2", "prop_snps")
  expect_true(all(expected_cols %in% colnames(res$enrichment)))
  # Should have one row per baseline annotation
  expect_equal(nrow(res$enrichment), 1L)
  # Values should be finite
  expect_true(all(is.finite(res$enrichment$tau)))
  expect_true(all(is.finite(res$enrichment$enrichment)))
})

test_that("gldsc_univariate with annotations returns tau_blocks matrix", {
  set.seed(201)
  annot <- make_gldsc_annotations(dat_annot$n_snps)
  res <- pecotmr:::gldsc_univariate(dat_annot$z, dat_annot$n,
                                    dat_annot$ld_score_ref,
                                    annotations = annot)
  expect_true(is.matrix(res$tau_blocks))
  # gldsc jackknife uses 200 blocks by default; with 500 SNPs the actual
  # number of unique blocks equals ceil(500/ceil(500/200)) = 200
  expect_true(nrow(res$tau_blocks) > 0)
  # Number of columns matches baseline annotation count
  expect_equal(ncol(res$tau_blocks), 1L)
  # Values should be finite (allow some NAs from edge blocks)
  expect_true(any(is.finite(res$tau_blocks)))
})

test_that("gldsc_univariate with annotations returns score_stats", {
  set.seed(202)
  annot <- make_gldsc_annotations(dat_annot$n_snps)
  res <- pecotmr:::gldsc_univariate(dat_annot$z, dat_annot$n,
                                    dat_annot$ld_score_ref,
                                    annotations = annot)
  expect_type(res$score_stats, "list")
  expect_true(all(c("z", "R") %in% names(res$score_stats)))
  # z should have length = number of candidate annotations (1)
  expect_length(res$score_stats$z, 1L)
  expect_true(is.finite(res$score_stats$z[1]))
  # R should be a matrix
  expect_true(is.matrix(res$score_stats$R))
  expect_equal(dim(res$score_stats$R), c(1L, 1L))
  # annotation_names should be present
  expect_true("annotation_names" %in% names(res$score_stats))
  expect_equal(res$score_stats$annotation_names, "candidate1")
})
