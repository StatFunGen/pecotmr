context("Fine-mapping methods parameter (SuSiE / SuSiE-inf / SuSiE-ash)")

# =============================================================================
# univariateAnalysisPipeline: methods parameter
# =============================================================================

# Helper: a minimal input set drawn from shipped example data, subsetted for
# fast vignette-style fits.
.make_uvp_inputs <- function(n_var = 200, seed = 42) {
  data(eqtl_region_example, envir = environment())
  X <- eqtl_region_example$X
  y <- as.numeric(eqtl_region_example$y_res)
  maf <- apply(X, 2, function(g) {
    f <- mean(g, na.rm = TRUE) / 2
    min(f, 1 - f)
  })
  set.seed(seed)
  sub <- sort(sample(ncol(X), n_var))
  list(X = X[, sub], y = y, maf = maf[sub])
}

test_that("univariateAnalysisPipeline: default behavior (addSusieInf = TRUE) fits two-stage", {
  skip_if_not_installed("susieR")
  inp <- .make_uvp_inputs()
  r <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$y, maf = inp$maf,
    twasWeights = FALSE,
    finemappingExtraOpts = list(refine = FALSE),
    signalCutoff = 0, verbose = 0
  )
  expect_false(is.null(r$susie_inf_fitted))
  expect_false(is.null(r$susie_fitted))
  expect_true(is.null(r$susie_ash_fitted))
  expect_true(all(c("susie", "susie_inf") %in% unique(r$top_loci$method)) ||
              "susie" %in% unique(r$top_loci$method))
})

test_that("univariateAnalysisPipeline: addSusieInf = FALSE fits plain SuSiE alone", {
  skip_if_not_installed("susieR")
  inp <- .make_uvp_inputs()
  r <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$y, maf = inp$maf,
    twasWeights = FALSE, addSusieInf = FALSE,
    finemappingExtraOpts = list(refine = FALSE),
    signalCutoff = 0, verbose = 0
  )
  expect_true(is.null(r$susie_inf_fitted))
  expect_false(is.null(r$susie_fitted))
  expect_equal(unique(r$top_loci$method), "susie")
})

test_that("univariateAnalysisPipeline: methods = \"susie_inf\" fits SuSiE-inf alone", {
  skip_if_not_installed("susieR")
  inp <- .make_uvp_inputs()
  r <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$y, maf = inp$maf,
    methods = c("susie_inf"), twasWeights = FALSE,
    finemappingExtraOpts = list(refine = FALSE),
    signalCutoff = 0, verbose = 0
  )
  expect_false(is.null(r$susie_inf_fitted))
  expect_true(is.null(r$susie_fitted))
  expect_true(is.null(r$susie_ash_fitted))
})

test_that("univariateAnalysisPipeline: methods = \"susie_ash\" fits SuSiE-ash alone", {
  skip_if_not_installed("susieR")
  inp <- .make_uvp_inputs()
  r <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$y, maf = inp$maf,
    methods = c("susie_ash"), twasWeights = FALSE,
    finemappingExtraOpts = list(refine = FALSE),
    signalCutoff = 0, verbose = 0
  )
  expect_true(is.null(r$susie_inf_fitted))
  expect_true(is.null(r$susie_fitted))
  expect_false(is.null(r$susie_ash_fitted))
  expect_equal(unique(r$top_loci$method), "susie_ash")
})

test_that("univariateAnalysisPipeline: methods = c(susie_inf, susie, susie_ash) returns rows for each", {
  skip_if_not_installed("susieR")
  inp <- .make_uvp_inputs()
  r <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$y, maf = inp$maf,
    methods = c("susie_inf", "susie", "susie_ash"),
    twasWeights = FALSE,
    finemappingExtraOpts = list(refine = FALSE),
    signalCutoff = 0, verbose = 0
  )
  expect_false(is.null(r$susie_inf_fitted))
  expect_false(is.null(r$susie_fitted))
  expect_false(is.null(r$susie_ash_fitted))
  # SuSiE-inf may report zero PIPs and contribute zero rows; expect at least susie + susie_ash
  expect_true(all(c("susie", "susie_ash") %in% unique(r$top_loci$method)))
})

test_that("univariateAnalysisPipeline: unknown method rejected", {
  skip_if_not_installed("susieR")
  inp <- .make_uvp_inputs()
  expect_error(
    univariateAnalysisPipeline(
      X = inp$X, Y = inp$y, maf = inp$maf,
      methods = c("not_a_real_method"),
      twasWeights = FALSE, verbose = 0
    ),
    "Unknown method"
  )
})

# Fifth fitting mode: SuSiE-ash with SuSiE-inf init -----------------------

# Helper: signal-rich synthetic data where SuSiE-inf finds non-zero mappable
# effects, so model_init genuinely carries information into SuSiE / SuSiE-ash.
# (The shipped eqtl_region_example is too weak: SuSiE-inf converges with all V=0
# / mu=0, which susieR's extract_model_init_fields correctly returns NULL for,
# collapsing chained and independent runs to identical results.)
.make_strong_signal_inputs <- function(seed = 1) {
  set.seed(seed)
  n <- 400; p <- 200
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p) * 1000, ":A:T")
  rownames(X) <- paste0("S", seq_len(n))
  beta_sparse <- rep(0, p); beta_sparse[c(20, 70, 130)] <- 2
  beta_poly   <- rnorm(p, sd = 0.05)
  y <- as.numeric(X %*% (beta_sparse + beta_poly) + rnorm(n))
  maf <- rep(0.3, p)
  list(X = X, y = y, maf = maf)
}

test_that("univariateAnalysisPipeline: chained SuSiE-inf -> SuSiE-ash differs from independent susie_ash", {
  skip_if_not_installed("susieR")
  inp <- .make_strong_signal_inputs()
  chained <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$y, maf = inp$maf,
    methods = c("susie_inf", "susie_ash"), addSusieInf = TRUE,
    twasWeights = FALSE,
    finemappingExtraOpts = list(refine = FALSE),
    signalCutoff = 0, verbose = 0
  )
  indep <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$y, maf = inp$maf,
    methods = c("susie_inf", "susie_ash"), addSusieInf = FALSE,
    twasWeights = FALSE,
    finemappingExtraOpts = list(refine = FALSE),
    signalCutoff = 0, verbose = 0
  )
  expect_true(all(c("susie_inf", "susie_ash") %in% unique(chained$top_loci$method)))
  expect_true(all(c("susie_inf", "susie_ash") %in% unique(indep$top_loci$method)))
  expect_false(is.null(chained$susie_inf_fitted))
  expect_false(is.null(chained$susie_ash_fitted))
  # SuSiE-inf must have non-trivial mappable effects for model_init to matter
  expect_true(any(chained$susie_inf_fitted$V > 0))
  # With non-trivial init, chained SuSiE-ash PIPs differ from independent run
  pip_chained <- chained$susie_ash_fitted$pip
  pip_indep   <- indep$susie_ash_fitted$pip
  expect_equal(length(pip_chained), length(pip_indep))
  expect_true(max(abs(pip_chained - pip_indep)) > 1e-3)
})

test_that("univariateAnalysisPipeline: chain dispatch emits the chained-init message", {
  skip_if_not_installed("susieR")
  inp <- .make_strong_signal_inputs()
  expect_message(
    univariateAnalysisPipeline(
      X = inp$X, Y = inp$y, maf = inp$maf,
      methods = c("susie_inf", "susie_ash"), addSusieInf = TRUE,
      twasWeights = FALSE,
      finemappingExtraOpts = list(refine = FALSE),
      signalCutoff = 0, verbose = 0
    ),
    "SuSiE-ash model initialized by SuSiE-inf"
  )
})

test_that("univariateAnalysisPipeline: add_susie_inf=FALSE fits SuSiE-ash without chained init", {
  skip_if_not_installed("susieR")
  inp <- .make_strong_signal_inputs()
  expect_message(
    univariateAnalysisPipeline(
      X = inp$X, Y = inp$y, maf = inp$maf,
      methods = c("susie_inf", "susie_ash"), addSusieInf = FALSE,
      twasWeights = FALSE,
      finemappingExtraOpts = list(refine = FALSE),
      signalCutoff = 0, verbose = 0
    ),
    "Fitting SuSiE-ash model on input data"
  )
})

test_that("univariateAnalysisPipeline: all three methods + chain initialises both susie and susie_ash", {
  skip_if_not_installed("susieR")
  inp <- .make_strong_signal_inputs()
  r <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$y, maf = inp$maf,
    methods = c("susie_inf", "susie", "susie_ash"), addSusieInf = TRUE,
    twasWeights = FALSE,
    finemappingExtraOpts = list(refine = FALSE),
    signalCutoff = 0, verbose = 0
  )
  expect_false(is.null(r$susie_inf_fitted))
  expect_false(is.null(r$susie_fitted))
  expect_false(is.null(r$susie_ash_fitted))
  expect_true(all(c("susie_inf", "susie", "susie_ash") %in% unique(r$top_loci$method)))
})

test_that("univariateAnalysisPipeline: twasWeights = TRUE requires chained two-stage", {
  skip_if_not_installed("susieR")
  inp <- .make_uvp_inputs()
  # methods = "susie" alone disables chaining => twasWeights must error
  expect_error(
    univariateAnalysisPipeline(
      X = inp$X, Y = inp$y, maf = inp$maf,
      methods = c("susie"), twasWeights = TRUE, verbose = 0
    ),
    "SuSiE-inf"
  )
  expect_error(
    univariateAnalysisPipeline(
      X = inp$X, Y = inp$y, maf = inp$maf,
      methods = c("susie_ash"), twasWeights = TRUE, verbose = 0
    ),
    "susie"
  )
})

# =============================================================================
# susieRssPipeline: methods parameter
# =============================================================================

.make_rss_inputs <- function() {
  data(gwas_sumstats_example, envir = environment())
  data(eqtl_region_example, envir = environment())
  X_ref <- eqtl_region_example$X
  R <- computeLd(X_ref, method = "sample")
  list(
    sumstats = gwas_sumstats_example,
    ldMat = R[gwas_sumstats_example$variant_id, gwas_sumstats_example$variant_id],
    n = nrow(X_ref)
  )
}

test_that("susieRssPipeline: default (methods = NULL) is single-method via analysis_method", {
  skip_if_not_installed("susieR")
  rss <- .make_rss_inputs()
  r <- susieRssPipeline(
    sumstats = rss$sumstats, ldMat = rss$ldMat, n = rss$n,
    L = 5, signalCutoff = 0
  )
  expect_equal(unique(r$top_loci$method), "susie_rss")
})

test_that("susieRssPipeline: legacy analysis_method = single_effect respected when methods = NULL", {
  skip_if_not_installed("susieR")
  rss <- .make_rss_inputs()
  r <- susieRssPipeline(
    sumstats = rss$sumstats, ldMat = rss$ldMat, n = rss$n,
    L = 5, analysisMethod = "single_effect", signalCutoff = 0
  )
  expect_equal(unique(r$top_loci$method), "single_effect")
})

test_that("susieRssPipeline: methods = \"susie_inf_rss\" alone", {
  skip_if_not_installed("susieR")
  rss <- .make_rss_inputs()
  r <- susieRssPipeline(
    sumstats = rss$sumstats, ldMat = rss$ldMat, n = rss$n,
    L = 5, methods = c("susie_inf_rss"), signalCutoff = 0
  )
  expect_equal(unique(r$top_loci$method), "susie_inf_rss")
})

test_that("susieRssPipeline: methods = \"susie_ash_rss\" alone", {
  skip_if_not_installed("susieR")
  rss <- .make_rss_inputs()
  r <- susieRssPipeline(
    sumstats = rss$sumstats, ldMat = rss$ldMat, n = rss$n,
    L = 5, methods = c("susie_ash_rss"), signalCutoff = 0
  )
  expect_equal(unique(r$top_loci$method), "susie_ash_rss")
})

test_that("susieRssPipeline: chained init produces both rows", {
  skip_if_not_installed("susieR")
  rss <- .make_rss_inputs()
  r <- susieRssPipeline(
    sumstats = rss$sumstats, ldMat = rss$ldMat, n = rss$n,
    L = 5, methods = c("susie_inf_rss", "susie_rss"),
    addSusieInf = TRUE, signalCutoff = 0
  )
  expect_true(all(c("susie_inf_rss", "susie_rss") %in% unique(r$top_loci$method)))
})

test_that("susieRssPipeline: independent both (addSusieInf = FALSE)", {
  skip_if_not_installed("susieR")
  rss <- .make_rss_inputs()
  r <- susieRssPipeline(
    sumstats = rss$sumstats, ldMat = rss$ldMat, n = rss$n,
    L = 5, methods = c("susie_inf_rss", "susie_rss"),
    addSusieInf = FALSE, signalCutoff = 0
  )
  expect_true(all(c("susie_inf_rss", "susie_rss") %in% unique(r$top_loci$method)))
})

test_that("susieRssPipeline: unknown method rejected", {
  skip_if_not_installed("susieR")
  rss <- .make_rss_inputs()
  expect_error(
    susieRssPipeline(
      sumstats = rss$sumstats, ldMat = rss$ldMat, n = rss$n,
      L = 5, methods = c("not_a_real_rss_method")
    ),
    "Unknown RSS method"
  )
})
