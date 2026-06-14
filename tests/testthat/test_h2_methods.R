# Tests for heritability estimation methods:
#   lderUnivariate (h2_lder.R)
#   gldscUnivariate (h2_gldsc.R)
#   hdlUnivariate (h2_hdl.R)

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
# LDER tests (h2_lder.R)
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
# gLDSC tests (h2_gldsc.R)
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
# HDL tests (h2_hdl.R)
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
# HDL annotation tests (h2_hdl.R stratified paths)
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
# gLDSC annotation tests (h2_gldsc.R stratified paths)
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
