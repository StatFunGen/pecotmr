# =============================================================================
# Tests for utility helpers exported from R/mashPipeline.R that don't
# require mashr / flashier installs:
#   sanitizeMashData, makePairwiseContrastCol, sliceMashData,
#   metaAnalysisPerCell
# =============================================================================

# ---------------------------------------------------------------------------
# sanitizeMashData
# ---------------------------------------------------------------------------

test_that("sanitizeMashData replaces NaN in bhat with 0", {
  d <- list(
    bhat  = matrix(c(1, NaN, 3, 4), nrow = 2),
    sbhat = matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2))
  out <- sanitizeMashData(d)
  expect_equal(out$bhat[1, 2], 3)         # untouched
  expect_equal(out$bhat[2, 1], 0)         # NaN -> 0
  expect_equal(out$sbhat, d$sbhat)        # sbhat untouched
})

test_that("sanitizeMashData replaces NaN/Inf in sbhat with 1e3", {
  d <- list(
    bhat  = matrix(c(1, 2, 3, 4), nrow = 2),
    sbhat = matrix(c(0.1, NaN, Inf, 0.4), nrow = 2))
  out <- sanitizeMashData(d)
  expect_equal(out$sbhat[2, 1], 1e3)
  expect_equal(out$sbhat[1, 2], 1e3)
  expect_equal(out$sbhat[1, 1], 0.1)
  expect_equal(out$sbhat[2, 2], 0.4)
  expect_equal(out$bhat, d$bhat)          # bhat untouched (no NaN there)
})

test_that("sanitizeMashData is idempotent on already-clean data", {
  d <- list(
    bhat  = matrix(c(1, 2, 3, 4), nrow = 2),
    sbhat = matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2))
  expect_equal(sanitizeMashData(d), d)
  expect_equal(sanitizeMashData(sanitizeMashData(d)), d)
})

test_that("sanitizeMashData leaves -Inf in bhat alone (only NaN is replaced)", {
  d <- list(
    bhat  = matrix(c(-Inf, Inf, 3, 4), nrow = 2),
    sbhat = matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2))
  out <- sanitizeMashData(d)
  expect_true(is.infinite(out$bhat[1, 1]))
  expect_true(is.infinite(out$bhat[2, 1]))
})

# ---------------------------------------------------------------------------
# makePairwiseContrastCol
# ---------------------------------------------------------------------------

test_that("makePairwiseContrastCol sets +1/-1 at named pair positions", {
  tmpl <- setNames(rep(0, 4), c("a", "b", "c", "d"))
  out  <- makePairwiseContrastCol(c("b", "d"), tmpl)
  expect_equal(out[["a"]],  0)
  expect_equal(out[["b"]],  1)
  expect_equal(out[["c"]],  0)
  expect_equal(out[["d"]], -1)
})

test_that("makePairwiseContrastCol preserves template names", {
  tmpl <- setNames(rep(0, 3), c("x", "y", "z"))
  out  <- makePairwiseContrastCol(c("x", "z"), tmpl)
  expect_equal(names(out), c("x", "y", "z"))
})

test_that("makePairwiseContrastCol overwrites pre-existing values in template", {
  tmpl <- setNames(c(5, -3, 7), c("a", "b", "c"))
  out  <- makePairwiseContrastCol(c("a", "c"), tmpl)
  expect_equal(out[["a"]],  1)            # was 5, now 1
  expect_equal(out[["b"]], -3)            # untouched
  expect_equal(out[["c"]], -1)            # was 7, now -1
})

# ---------------------------------------------------------------------------
# sliceMashData
# ---------------------------------------------------------------------------

test_that("sliceMashData subsets bhat / sbhat / Z by SNP and sample", {
  snps    <- c("s1", "s2", "s3")
  samples <- c("ctxA", "ctxB", "ctxC")
  data <- list(
    bhat  = matrix(seq_len(9),  3, 3, dimnames = list(snps, samples)),
    sbhat = matrix(seq_len(9) / 10, 3, 3, dimnames = list(snps, samples)),
    Z     = matrix(seq_len(9) * 2, 3, 3, dimnames = list(snps, samples)),
    snp   = snps)
  vhat <- diag(1, 3, 3); dimnames(vhat) <- list(samples, samples)

  out <- sliceMashData(data, vhat,
                       snps = c("s1", "s3"),
                       samples = c("ctxA", "ctxC"))
  expect_equal(dim(out$data$bhat), c(2, 2))
  expect_equal(dim(out$vhat),      c(2, 2))
  expect_equal(colnames(out$data$bhat),  c("ctxA", "ctxC"))
  expect_equal(colnames(out$data$sbhat), c("ctxA", "ctxC"))
  expect_equal(colnames(out$data$Z),     c("ctxA", "ctxC"))
  expect_equal(colnames(out$vhat),       c("ctxA", "ctxC"))
  expect_equal(out$data$snp, c("s1", "s3"))
})

test_that("sliceMashData restricts data$snp to intersection of snps argument", {
  snps    <- c("s1", "s2", "s3", "s4")
  samples <- c("ctxA", "ctxB")
  data <- list(
    bhat  = matrix(1, 4, 2, dimnames = list(snps, samples)),
    sbhat = matrix(1, 4, 2, dimnames = list(snps, samples)),
    Z     = matrix(1, 4, 2, dimnames = list(snps, samples)),
    snp   = snps)
  vhat <- diag(1, 2, 2); dimnames(vhat) <- list(samples, samples)
  out  <- sliceMashData(data, vhat,
                        snps    = c("s2", "s4"),
                        samples = samples)
  expect_equal(out$data$snp, c("s2", "s4"))
})

# ---------------------------------------------------------------------------
# metaAnalysisPerCell
# ---------------------------------------------------------------------------

test_that("metaAnalysisPerCell returns single-effect p-value when only one feature passes filter", {
  feat <- "var1"
  cols <- c("mean_contrast_brain_vs_blood")
  es <- matrix(0.5, nrow = 1, ncol = 1, dimnames = list(feat, cols))
  se <- matrix(0.1, nrow = 1, ncol = 1, dimnames = list(feat, cols))
  out <- metaAnalysisPerCell(es, se)
  # 2 conditions: brain, blood -> 2 rows
  expect_equal(nrow(out), 2L)
  expect_true(all(c("cell", "condition", "meta_pvalue",
                    "meta_effect", "meta_se", "tau2", "I2") %in% names(out)))
  # With a single effect both rows return single-effect p-values (not NA)
  expect_false(any(is.na(out$meta_pvalue)))
  # tau2 / I2 only meaningful with >= 2 effects
  expect_true(all(is.na(out$tau2)))
  expect_true(all(is.na(out$I2)))
})

test_that("metaAnalysisPerCell returns NA pvalue when SE cutoff drops everything", {
  feat <- c("v1", "v2")
  cols <- "mean_contrast_brain_vs_blood"
  es <- matrix(c(0.1, 0.2), nrow = 2, ncol = 1, dimnames = list(feat, cols))
  se <- matrix(c(0.01, 0.02), nrow = 2, ncol = 1, dimnames = list(feat, cols))
  # seCutoff = 0.5 drops both rows
  out <- metaAnalysisPerCell(es, se, seCutoff = 0.5)
  expect_true(all(is.na(out$meta_pvalue)))
  expect_true(all(is.na(out$meta_effect)))
})

test_that("metaAnalysisPerCell runs DerSimonian-Laird when >=2 effects survive", {
  feat <- c("v1", "v2", "v3")
  cols <- "mean_contrast_brain_vs_blood"
  es <- matrix(c(0.3, 0.5, 0.4), nrow = 3, ncol = 1, dimnames = list(feat, cols))
  se <- matrix(c(0.1, 0.1, 0.1), nrow = 3, ncol = 1, dimnames = list(feat, cols))
  out <- metaAnalysisPerCell(es, se)
  expect_equal(nrow(out), 2L)
  # All three effects survive (SE = 0.1 > 0): meta_effect and meta_se populated
  expect_false(any(is.na(out$meta_effect)))
  expect_false(any(is.na(out$meta_se)))
  expect_false(any(is.na(out$tau2)))
  expect_false(any(is.na(out$I2)))
  # I2 in [0, 1]
  expect_true(all(out$I2 >= 0 & out$I2 <= 1))
})

test_that("metaAnalysisPerCell unique-cells extraction handles >2 conditions", {
  feat <- "v1"
  cols <- c(
    "mean_contrast_brain_vs_blood",
    "mean_contrast_brain_vs_muscle",
    "mean_contrast_blood_vs_muscle")
  es <- matrix(c(0.3, 0.4, 0.1), nrow = 1, dimnames = list(feat, cols))
  se <- matrix(c(0.1, 0.1, 0.1), nrow = 1, dimnames = list(feat, cols))
  out <- metaAnalysisPerCell(es, se)
  # 3 cells (brain, blood, muscle), each with 2 vs-comparisons -> 6 rows
  expect_setequal(unique(out$cell), c("brain", "blood", "muscle"))
  expect_equal(nrow(out), 6L)
})

# ---------------------------------------------------------------------------
# updateMashModelCov — works on a hand-built mock fitted_g (no mashr dep)
# ---------------------------------------------------------------------------

test_that("updateMashModelCov drops dropped conditions + resizes remaining cov matrices", {
  R <- 3L; samples <- c("brain", "blood", "muscle")
  U <- list(
    brain    = diag(c(1, 0, 0)),
    blood    = diag(c(0, 1, 0)),
    muscle   = diag(c(0, 0, 1)),
    identity = diag(1, R),
    PCA_1    = matrix(seq_len(R * R), R, R))   # no dimnames -> last branch
  pi <- setNames(rep(0.2, 5L),
                  c("brain.scale1", "blood.scale1",
                    "muscle.scale1", "identity.scale1", "PCA_1.scale1"))
  m <- list(fitted_g = list(Ulist = U, pi = pi))
  m2 <- updateMashModelCov(m,
                            allSamples = samples,
                            samples    = c("brain", "blood"))
  expect_false("muscle" %in% names(m2$fitted_g$Ulist))
  expect_true(all(vapply(m2$fitted_g$Ulist,
                          function(x) all(dim(x) == c(2L, 2L)),
                          logical(1))))
  expect_false(any(grepl("muscle", names(m2$fitted_g$pi))))
  # Brain matrix has a single 1 at the brain position (the first of the
  # retained `samples` ordering).
  expect_equal(m2$fitted_g$Ulist$brain[1, 1], 1)
  expect_equal(sum(m2$fitted_g$Ulist$brain), 1)
})

# ---------------------------------------------------------------------------
# fitMashContrast — fabricated posterior inputs (no mashr dep)
# ---------------------------------------------------------------------------

test_that("fitMashContrast returns NULL when fewer than 2 tested conditions", {
  origMean <- matrix(0, nrow = 1, ncol = 3,
                      dimnames = list("v1", c("a", "b", "c")))
  origMean[1, "b"] <- 1
  pm <- matrix(0, nrow = 1, ncol = 3,
                dimnames = list("v1", c("a", "b", "c")))
  pv <- array(diag(3), dim = c(3, 3, 1))
  dimnames(pv) <- list(c("a", "b", "c"), c("a", "b", "c"), NULL)
  expect_null(fitMashContrast(1L, origMean, pm, pv))
})

test_that("fitMashContrast: 2-tested-conditions fast path yields one pairwise contrast", {
  origMean <- matrix(c(0.5, 0.3, 0), nrow = 1,
                      dimnames = list("v1", c("a", "b", "c")))
  pm <- matrix(c(0.5, 0.3, 0), nrow = 1,
                dimnames = list("v1", c("a", "b", "c")))
  pv <- array(diag(3) * 0.1, dim = c(3, 3, 1))
  dimnames(pv) <- list(c("a", "b", "c"), c("a", "b", "c"), NULL)
  out <- fitMashContrast(1L, origMean, pm, pv)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1L)
  # 2 tested -> 1 pairwise contrast -> 3 columns (mean, se, p)
  expect_equal(ncol(out), 3L)
  expect_setequal(names(out),
                   c("mean_contrast_a_vs_b",
                     "se_contrast_a_vs_b",
                     "p_contrast_a_vs_b"))
  expect_equal(out[["mean_contrast_a_vs_b"]], 0.5 - 0.3)
})

test_that("fitMashContrast: 3-tested-conditions yields deviation + pairwise contrasts", {
  origMean <- matrix(c(0.5, 0.3, -0.2), nrow = 1,
                      dimnames = list("v1", c("a", "b", "c")))
  pm <- matrix(c(0.5, 0.3, -0.2), nrow = 1,
                dimnames = list("v1", c("a", "b", "c")))
  pv <- array(diag(3) * 0.1, dim = c(3, 3, 1))
  dimnames(pv) <- list(c("a", "b", "c"), c("a", "b", "c"), NULL)
  out <- fitMashContrast(1L, origMean, pm, pv)
  expect_s3_class(out, "data.frame")
  # 3 deviation + choose(3,2)=3 pairwise = 6 contrasts -> 18 cols
  expect_equal(ncol(out), 18L)
  contrastSuffix <- sub("^(mean|se|p)_contrast_", "", names(out))
  expect_true(any(grepl("_deviation$", contrastSuffix)))
  expect_true(any(grepl("_vs_", contrastSuffix)))
})

# ---------------------------------------------------------------------------
# mashPipeline — end-to-end on the bundled multi-context example
# ---------------------------------------------------------------------------

test_that("mashPipeline runs end-to-end on qtl_sumstats_multicontext_example", {
  skip_if_not_installed("mashr")
  skip_if_not_installed("flashier")
  data(qtl_sumstats_multicontext_example)
  ss <- qtl_sumstats_multicontext_example
  # Use the same fixture for strong/random; nPcs <= ncol - 1 (3 contexts).
  res <- suppressMessages(suppressWarnings(
    mashPipeline(
      sumStatsList = list(strong = ss, random = ss),
      alpha        = 0,
      nPcs         = 2L,
      setSeed      = 1L)))
  expect_named(res, c("U", "w"))
  expect_type(res$U, "list")
  expect_gt(length(res$U), 0L)
  # Every covariance matrix is 3x3 (one row/col per context)
  expect_true(all(vapply(res$U,
                          function(m) all(dim(m) == c(3L, 3L)),
                          logical(1))))
  expect_type(res$w, "double")
  expect_equal(sum(res$w), 1, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# fitMashContrast — condition grouping (>2 conditions, grouped replicates)
# ---------------------------------------------------------------------------

test_that("fitMashContrast applies deviation + pairwise group adjustments", {
  conds <- c("a", "b", "c", "d")
  origMean <- matrix(c(0.5, 0.3, -0.2, 0.4), nrow = 1,
                     dimnames = list("v1", conds))
  pm <- matrix(c(0.5, 0.3, -0.2, 0.4), nrow = 1,
               dimnames = list("v1", conds))
  pv <- array(0, dim = c(4, 4, 1), dimnames = list(conds, conds, NULL))
  pv[, , 1] <- diag(4) * 0.1
  # a,b share group 1 (replicates); c is its own group 2; d ungrouped (0).
  # Non-NULL grouping triggers `grouping <- grouping[tested]`, the >2-condition
  # deviation re-weighting loop, and the pairwise group-adjustment loop.
  grouping <- setNames(c(1L, 1L, 2L, 0L), conds)
  out <- fitMashContrast(1L, origMean, pm, pv, grouping = grouping)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1L)
  # 4 deviation + choose(4,2)=6 pairwise = 10 contrasts -> 30 columns.
  expect_equal(ncol(out), 30L)
  contrastSuffix <- sub("^(mean|se|p)_contrast_", "", names(out))
  expect_true(any(grepl("_deviation$", contrastSuffix)))
  expect_true(any(grepl("_vs_", contrastSuffix)))
  expect_true(all(is.finite(unlist(out))))
})

# ---------------------------------------------------------------------------
# updateMashModelCov — named data-driven cov matrices (the `[samples, samples]`
# else branch, distinct from the no-dimnames positional-slice branch)
# ---------------------------------------------------------------------------

test_that("updateMashModelCov slices named data-driven cov matrices by sample", {
  allSamples <- c("brain", "blood", "muscle")
  ddMat <- matrix(seq_len(9), 3, 3, dimnames = list(allSamples, allSamples))
  U  <- list(identity = diag(1, 3), dataDriven = ddMat)
  pi <- setNames(c(0.5, 0.5), c("identity.scale1", "dataDriven.scale1"))
  m  <- list(fitted_g = list(Ulist = U, pi = pi))
  m2 <- updateMashModelCov(m, allSamples = allSamples,
                           samples = c("brain", "muscle"))
  # The named data-driven matrix is sliced by name: cov[[d]][samples, samples].
  expect_equal(dim(m2$fitted_g$Ulist$dataDriven), c(2L, 2L))
  expect_equal(m2$fitted_g$Ulist$dataDriven,
               ddMat[c("brain", "muscle"), c("brain", "muscle")])
  # identity collapses to a single 1 in the top-left corner.
  expect_equal(m2$fitted_g$Ulist$identity[1, 1], 1)
  expect_equal(sum(m2$fitted_g$Ulist$identity), 1)
})

# ---------------------------------------------------------------------------
# metaAnalysisPerCell — cells matching no contrast column are skipped
# ---------------------------------------------------------------------------

test_that("metaAnalysisPerCell skips cells whose name matches no column", {
  # The condition "x$y" yields a derived cell name "x$y"; used as a grep()
  # pattern the embedded `$` anchor matches nothing, exercising the
  # `if (length(cellIdx) == 0) next` skip branch.
  cols <- "mean_contrast_x$y_vs_z"
  es <- matrix(c(0.3, 0.5), nrow = 2, dimnames = list(c("v1", "v2"), cols))
  se <- matrix(c(0.1, 0.1), nrow = 2, dimnames = list(c("v1", "v2"), cols))
  out <- metaAnalysisPerCell(es, se)
  # "x$y" is skipped; only the "z" cell survives.
  expect_false("x$y" %in% out$cell)
  expect_true("z" %in% out$cell)
  expect_equal(nrow(out), 1L)
})

# ---------------------------------------------------------------------------
# mashPipeline — input validation (errors fire before any mashr call).
# Guarded by skips because the requireNamespace() checks run first; without
# mashr/flashier the function stops with an install message instead.
# ---------------------------------------------------------------------------

test_that("mashPipeline rejects a sumStatsList that is not a named list", {
  skip_if_not_installed("mashr")
  skip_if_not_installed("flashier")
  expect_error(mashPipeline(list(1, 2), alpha = 0),
               "must be a named list")
})

test_that("mashPipeline errors when a required entry is missing", {
  skip_if_not_installed("mashr")
  skip_if_not_installed("flashier")
  expect_error(mashPipeline(list(strong = 1), alpha = 0),
               "missing required entr")
})

test_that("mashPipeline errors on unrecognised entries", {
  skip_if_not_installed("mashr")
  skip_if_not_installed("flashier")
  expect_error(
    mashPipeline(list(strong = 1, random = 1, bogus = 1), alpha = 0),
    "unrecognised entries")
})

test_that("mashPipeline coerces a SimpleList before validating its names", {
  skip_if_not_installed("mashr")
  skip_if_not_installed("flashier")
  skip_if_not_installed("S4Vectors")
  # A SimpleList is converted to a base list first (the as.list branch), then
  # validation runs; here it is missing both required entries.
  expect_error(
    mashPipeline(S4Vectors::SimpleList(bogus = 1), alpha = 0),
    "missing required entr")
})

# ---------------------------------------------------------------------------
# mashPipeline — priorCovariances validation + bypass path.
# Supplying residualCorrelation makes random/null optional and short-circuits
# null-correlation estimation, so the supplied Vhat branch is also exercised.
# ---------------------------------------------------------------------------

test_that("mashPipeline rejects priorCovariances not a non-empty named list", {
  skip_if_not_installed("mashr")
  skip_if_not_installed("flashier")
  data(qtl_sumstats_multicontext_example)
  ss <- qtl_sumstats_multicontext_example
  vhat <- diag(3)
  # Empty list.
  expect_error(suppressMessages(suppressWarnings(
    mashPipeline(list(strong = ss), alpha = 0,
                 residualCorrelation = vhat, priorCovariances = list()))),
    "non-empty named")
  # Unnamed list.
  expect_error(suppressMessages(suppressWarnings(
    mashPipeline(list(strong = ss), alpha = 0,
                 residualCorrelation = vhat,
                 priorCovariances = list(diag(3))))),
    "non-empty named")
})

test_that("mashPipeline rejects priorCovariances with wrong dimensions", {
  skip_if_not_installed("mashr")
  skip_if_not_installed("flashier")
  data(qtl_sumstats_multicontext_example)
  ss <- qtl_sumstats_multicontext_example
  vhat <- diag(3)
  expect_error(suppressMessages(suppressWarnings(
    mashPipeline(list(strong = ss), alpha = 0,
                 residualCorrelation = vhat,
                 priorCovariances = list(myU = diag(2))))),
    "3 x 3 matrix")
})

test_that("mashPipeline passes supplied residualCorrelation + priorCovariances through", {
  skip_if_not_installed("mashr")
  skip_if_not_installed("flashier")
  data(qtl_sumstats_multicontext_example)
  ss <- qtl_sumstats_multicontext_example
  vhat <- diag(3)
  U0 <- list(identity = diag(3), effectA = diag(c(1, 0, 0)))
  res <- suppressMessages(suppressWarnings(
    mashPipeline(list(strong = ss), alpha = 0,
                 residualCorrelation = vhat, priorCovariances = U0)))
  expect_named(res, c("U", "w"))
  # priorCovariances passed straight through as the covariance list (bypass of
  # the cov_canonical / cov_pca / cov_flash / cov_ed chain).
  expect_identical(res$U, U0)
  expect_type(res$w, "double")
  expect_equal(sum(res$w), 1, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# mashPipeline — null-based Vhat estimation + default nPcs
# ---------------------------------------------------------------------------

test_that("mashPipeline estimates Vhat from a null set and defaults nPcs", {
  skip_if_not_installed("mashr")
  skip_if_not_installed("flashier")
  data(qtl_sumstats_multicontext_example)
  ss <- qtl_sumstats_multicontext_example
  # Supplying `null` triggers estimate_null_correlation_simple; leaving nPcs
  # NULL exercises the `nPcs <- ncol(Bhat) - 1` default in the cov_* chain.
  res <- suppressMessages(suppressWarnings(
    mashPipeline(list(strong = ss, random = ss, null = ss),
                 alpha = 0, setSeed = 1L)))
  expect_named(res, c("U", "w"))
  expect_gt(length(res$U), 0L)
  expect_true(all(vapply(res$U,
                         function(m) all(dim(m) == c(3L, 3L)),
                         logical(1))))
  expect_equal(sum(res$w), 1, tolerance = 1e-6)
})
