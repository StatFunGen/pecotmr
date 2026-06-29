context("SS-TWAS: weights, pipeline, and omnibus combination")

# Previous TwasWeights-class tests used the legacy constructor
# `TwasWeights(weights = list(...), variantIds = ..., standardized = ...)`.
# The new `TwasWeights` is a DFrame collection class with (study,
# context, trait, method, entry) columns where each entry is a
# `TwasWeightsEntry` S4 object carrying weights / fits / cvResult.
# Class-shape tests for the new collection should live alongside the
# pipeline tests and assert via accessors (`getWeights`, `getStudy`,
# `getCvResult`, etc.) — not against legacy slot shapes.
#
# `twasAnalysis()` was collapsed into the unified `twasZ()` dispatcher
# (task #37); its tests are removed here.
# `twasWeightsSumstatPipeline()` was removed without replacement in the
# S4 refactor (twasWeightsPipeline now dispatches directly on
# `QtlSumStats` / `QtlDataset` / `MultiStudyQtlDataset`).
#
# What remains: tests of the internal SuSiE-RSS weight extractors that
# are still present in `R/twasWeights.R` (`.susieRssExtractWeights`,
# `susieRssWeights`, `susieInfRssWeights`, `fitSusieInfThenSusieRss`).

# =============================================================================
# SuSiE-RSS weight extraction
# =============================================================================




test_that("mrmashWeights fitDetail: slim default omits the full fit, full keeps it", {
  skip_if_not_installed("mr.mashr")
  fakeFit <- list(w0 = c(a_1 = 0.5, a_2 = 0.5), V = diag(2))
  ddpm <- list(U = list(a = diag(2)))
  # Mock coef extraction so we exercise only the retain payload logic, not a
  # real mr.mash fit. coef.mr.mash(fit)[-1, ] -> drop the intercept row.
  local_mocked_bindings(
    coef.mr.mash = function(object, ...) rbind(c(0, 0), c(0.1, 0.2)),
    .package = "mr.mashr")
  fitSlim <- attr(mrmashWeights(mrmashFit = fakeFit, retainFit = TRUE,
                                dataDrivenPriorMatrices = ddpm), "fit")
  expect_setequal(names(fitSlim), c("dataDrivenPriorMatrices", "w0", "V"))
  expect_null(fitSlim$fit)                       # slim: no full fit
  expect_identical(fitSlim$dataDrivenPriorMatrices, ddpm)
  expect_identical(fitSlim$w0, fakeFit$w0)

  fitFull <- attr(mrmashWeights(mrmashFit = fakeFit, retainFit = TRUE,
                                fitDetail = "full",
                                dataDrivenPriorMatrices = ddpm), "fit")
  expect_true("fit" %in% names(fitFull))         # full: the whole fit retained
  expect_identical(fitFull$fit, fakeFit)
  expect_identical(fitFull$w0, fakeFit$w0)       # slim fields still present
})


# =============================================================================
# Two-stage SuSiE-RSS fitting
# =============================================================================

test_that("fitSusieInfThenSusieRss returns two fits", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 20
  n <- 500
  R <- diag(p)
  z <- rnorm(p)
  fits <- fitSusieInfThenSusieRss(z, R, n, args = list(L = 5))
  expect_true(is.list(fits))
  expect_true("susie" %in% names(fits))
  expect_true("susieInf" %in% names(fits))
  expect_true("susieInf" %in% class(fits$susieInf))
  expect_true("susieRss" %in% class(fits$susie))
})

# === Tests migrated from test_mrmashWrapper.R (mr.mash + glasso/glmnet coef helpers) ===

test_that("compute_w0 returns uniform weights when ncomps == 1", {
  Bhat <- matrix(c(1, 0, 0, 2, 0, 0), nrow = 3, ncol = 2)
  result <- pecotmr:::compute_w0(Bhat, ncomps = 1)
  expect_equal(result, 1)
})


test_that("compute_w0 handles all-zero Bhat by returning uniform weights", {
  # When Bhat is all zero, prop_nonzero = 0
  # w0 = c(1, 0, ..., 0) => sum(w0 != 0) < 2 => fallback to uniform
  Bhat <- matrix(0, nrow = 5, ncol = 3)
  result <- pecotmr:::compute_w0(Bhat, ncomps = 4)
  expect_equal(result, rep(1 / 4, 4))
  expect_equal(sum(result), 1)
})


test_that("compute_w0 distributes weight based on nonzero rows when ncomps > 1", {
  # 2 out of 4 rows have nonzero entries
  Bhat <- matrix(0, nrow = 4, ncol = 2)
  Bhat[1, 1] <- 1
  Bhat[3, 2] <- 2
  result <- pecotmr:::compute_w0(Bhat, ncomps = 3)
  expect_equal(length(result), 3)
  expect_equal(sum(result), 1, tolerance = 1e-10)
  # First element should be (1 - prop_nonzero) = 0.5
  expect_equal(result[1], 0.5)
})

# =========================================================================
# mrmashWrapper.R: rescale_cov_w0 (lines 300-329)
# =========================================================================


test_that("rescale_cov_w0 removes null component and renormalizes", {
  w0 <- c(null = 0.3, XtX_1 = 0.2, XtX_2 = 0.1, FLASH_1 = 0.15, FLASH_2 = 0.25)
  result <- pecotmr:::rescale_cov_w0(w0)
  expect_false("null" %in% names(result))
  expect_equal(sum(result), 1, tolerance = 1e-10)
})


test_that("rescale_cov_w0 handles all-zero non-null weights", {
  w0 <- c(null = 1.0, XtX_1 = 0, XtX_2 = 0, FLASH_1 = 0)
  result <- pecotmr:::rescale_cov_w0(w0)
  # All non-null weights are zero -> equal weights
  expect_equal(sum(result), 1, tolerance = 1e-10)
  expect_true(all(result == result[1]))  # all equal
})


test_that("rescale_cov_w0 groups correctly by prior group prefix", {
  w0 <- c(null = 0.5, PCA_1 = 0.1, PCA_2 = 0.2, tFLASH_1 = 0.1, tFLASH_2 = 0.1)
  result <- pecotmr:::rescale_cov_w0(w0)
  expect_true("PCA" %in% names(result))
  expect_true("tFLASH" %in% names(result))
  expect_equal(sum(result), 1, tolerance = 1e-10)
})

# =========================================================================
# mrmashWrapper.R: compute_grid, grid_min, grid_max, autoselect_mixsd
# (lines 333-372)
# =========================================================================


test_that("grid_max returns scaled grid_min when bhat^2 <= sbhat^2", {
  bhat <- c(0.1, 0.2)
  sbhat <- c(1.0, 1.0)
  result <- pecotmr:::gridMax(bhat, sbhat)
  expect_equal(result, 8 * pecotmr:::gridMin(bhat, sbhat))
})


test_that("grid_max returns 2*sqrt(max(bhat^2 - sbhat^2)) otherwise", {
  bhat <- c(5, 1)
  sbhat <- c(0.5, 0.5)
  expected <- 2 * sqrt(max(bhat^2 - sbhat^2))
  result <- pecotmr:::gridMax(bhat, sbhat)
  expect_equal(result, expected)
})


test_that("autoselect_mixsd returns 2-element vector when mult == 0", {
  result <- pecotmr:::autoselectMixsd(0.01, 1.0, mult = 0)
  expect_equal(result, c(0, 0.5))
})


test_that("autoselect_mixsd returns valid grid with sqrt(2) mult", {
  result <- pecotmr:::autoselectMixsd(0.01, 1.0, mult = sqrt(2))
  expect_true(length(result) > 1)
  expect_equal(result[length(result)], 1.0)  # last element is gmax
  # All elements should be positive

  expect_true(all(result > 0))
})


test_that("compute_grid produces a valid grid from summary statistics", {
  set.seed(42)
  bhat <- matrix(rnorm(20, sd = 2), nrow = 10, ncol = 2)
  sbhat <- matrix(abs(rnorm(20, mean = 0.5, sd = 0.1)), nrow = 10, ncol = 2)
  result <- pecotmr:::computeGrid(bhat, sbhat)
  expect_true(is.numeric(result))
  expect_true(length(result) > 0)
  expect_true(all(result > 0))
})


test_that("compute_grid handles NA and zero sbhat values", {
  bhat <- c(1, 2, NA, 4, 5)
  sbhat <- c(0.5, 0, NA, 0.3, 0.8)
  result <- pecotmr:::computeGrid(bhat, sbhat)
  expect_true(is.numeric(result))
  expect_true(length(result) > 0)
})

# =========================================================================
# mrmashWrapper.R: mrmashWrapper input validation (lines 99-130)
# =========================================================================

# Note: Cannot mock requireNamespace via local_mocked_bindings because it is a
# base R function, not in pecotmr's namespace. We skip these tests and instead
# test the downstream validation that we CAN exercise.


test_that("mrmashWrapper errors when X and Y are not matrices", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("mr.mashr")
  expect_error(mrmashWrapper(data.frame(x = 1:3), matrix(1:6, nrow = 3, ncol = 2)),
               "matrices")
})


test_that("mrmashWrapper errors when X and Y row counts differ", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("mr.mashr")
  expect_error(mrmashWrapper(matrix(1:6, nrow = 3, ncol = 2), matrix(1:8, nrow = 4, ncol = 2)),
               "same number of rows")
})


test_that("mrmashWrapper errors when prior_grid is not a vector", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("mr.mashr")
  X <- matrix(rnorm(12), nrow = 3, ncol = 4)
  Y <- matrix(rnorm(6), nrow = 3, ncol = 2)
  expect_error(mrmashWrapper(X, Y, priorGrid = matrix(1:4, nrow = 2)),
               "priorGrid must be a vector")
})


test_that("mrmashWrapper errors when no prior matrices and canonical_prior_matrices is FALSE", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("mr.mashr")
  X <- matrix(rnorm(12), nrow = 3, ncol = 4)
  Y <- matrix(rnorm(6), nrow = 3, ncol = 2)
  expect_error(mrmashWrapper(X, Y, dataDrivenPriorMatrices = NULL,
                               canonicalPriorMatrices = FALSE),
               "dataDrivenPriorMatrices")
})


test_that("mrmashWrapper warns when Y has missing and B_init_method is glasso", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("mr.mashr")
  set.seed(42)
  n <- 20; p <- 5; r <- 2
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
  Y[1, 1] <- NA  # introduce missing values
  colnames(Y) <- c("cond1", "cond2")

  # Should produce warning about glasso and NAs, then likely fail on the
  # downstream mr.mashr call, but the warning is what we test
  expect_warning(
    tryCatch(
      mrmashWrapper(X, Y, bInitMethod = "glasso",
                     dataDrivenPriorMatrices = list(U = list(matrix(1, 2, 2))),
                     canonicalPriorMatrices = FALSE),
      error = function(e) NULL
    ),
    "glasso"
  )
})

# =========================================================================
# mrmashWrapper.R: computeCoefficientsGlasso (lines 211-240)
# =========================================================================


test_that("computeCoefficientsGlasso runs without Xnew", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 5; r <- 3
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
  colnames(Y) <- paste0("cond", 1:r)
  result <- pecotmr:::computeCoefficientsGlasso(X, Y, standardize = FALSE,
                                                   nthreads = 1, Xnew = NULL)
  expect_true("Bhat" %in% names(result))
  expect_true("Ytrain" %in% names(result))
  expect_equal(nrow(result$Bhat), p)
  expect_equal(ncol(result$Bhat), r)
  expect_null(result$Yhat_new)
})


test_that("computeCoefficientsGlasso runs with Xnew", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 5; r <- 3
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
  colnames(Y) <- paste0("cond", 1:r)
  Xnew <- matrix(rnorm(10 * p), nrow = 10, ncol = p)
  result <- pecotmr:::computeCoefficientsGlasso(X, Y, standardize = FALSE,
                                                   nthreads = 1, Xnew = Xnew)
  expect_true("Yhat_new" %in% names(result))
  expect_equal(nrow(result$Yhat_new), 10)
  expect_equal(ncol(result$Yhat_new), r)
  expect_equal(colnames(result$Yhat_new), colnames(Y))
})

# =========================================================================
# mrmashWrapper.R: computeCoefficientsUnivGlmnet (lines 243-281)
# =========================================================================


test_that("computeCoefficientsUnivGlmnet runs without Xnew", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 60; p <- 5; r <- 2
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
  colnames(Y) <- paste0("cond", 1:r)
  result <- pecotmr:::computeCoefficientsUnivGlmnet(X, Y, alpha = 0.5,
                                                        standardize = FALSE,
                                                        nthreads = 1, Xnew = NULL)
  expect_true("Bhat" %in% names(result))
  expect_true("intercept" %in% names(result))
  expect_equal(nrow(result$Bhat), p)
  expect_equal(ncol(result$Bhat), r)
  expect_null(result$Yhat_new)
})


test_that("computeCoefficientsUnivGlmnet runs with Xnew", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 60; p <- 5; r <- 2
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
  colnames(Y) <- paste0("cond", 1:r)
  Xnew <- matrix(rnorm(8 * p), nrow = 8, ncol = p)
  result <- pecotmr:::computeCoefficientsUnivGlmnet(X, Y, alpha = 0.5,
                                                        standardize = FALSE,
                                                        nthreads = 1, Xnew = Xnew)
  expect_true("Yhat_new" %in% names(result))
  expect_equal(nrow(result$Yhat_new), 8)
  expect_equal(ncol(result$Yhat_new), r)
  expect_equal(colnames(result$Yhat_new), colnames(Y))
})


test_that("computeCoefficientsUnivGlmnet handles NA in Y", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 60; p <- 5; r <- 2
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
  colnames(Y) <- paste0("cond", 1:r)
  Y[1:5, 1] <- NA  # introduce missing values in one condition
  result <- pecotmr:::computeCoefficientsUnivGlmnet(X, Y, alpha = 0.5,
                                                        standardize = FALSE,
                                                        nthreads = 1, Xnew = NULL)
  expect_true("Bhat" %in% names(result))
  expect_equal(nrow(result$Bhat), p)
})

# =========================================================================
# mrmashWrapper.R: mrmashWrapper seed warning (line 107-108)
# =========================================================================

# Note: Cannot mock base::exists() via local_mocked_bindings.
# The seed-check message on line 107-108 would require removing .Random.seed
# from the global environment, which is not safe to do in tests.

# =============================================================================
# Real-fit coverage for the solver wrappers in R/regularizedRegressionWrappers.R
# -----------------------------------------------------------------------------
# The fine-mapping / TWAS pipelines MOCK these wrappers, so their bodies are
# otherwise untested. Here we drive each wrapper on a SMALL real fixture and
# assert the return shape (weight length == #variants; matrix for multivariate;
# attr(.,"fit") when retainFit = TRUE). Bayesian/MCMC iterations are kept tiny.
# fsusieWeights and the mock-based mrmash/mvsusie payload tests live in
# test_rrMrmashMvsusie.R and are not duplicated here.
# =============================================================================

# (Shared .rrwXy / .rrwStatLd / .rrwMulti fixtures live in helper-rrwFixtures.R
#  so test_fineMappingWrappers.R can reuse them for the SuSiE weight extractors.)

# -------------------------------- individual --------------------------------

test_that("lassoWeights / enetWeights (glmnet) return length-p weights", {
  skip_if_not_installed("glmnet")
  f <- .rrwXy()
  expect_length(as.numeric(lassoWeights(f$X, f$y)), f$p)
  expect_length(as.numeric(enetWeights(f$X, f$y)), f$p)
})

test_that("scadWeights / mcpWeights (ncvreg) return length-p weights", {
  skip_if_not_installed("ncvreg")
  f <- .rrwXy()
  expect_length(as.numeric(scadWeights(f$X, f$y)), f$p)
  expect_length(as.numeric(mcpWeights(f$X, f$y)), f$p)
})

test_that("l0learnWeights returns length-p weights", {
  skip_if_not_installed("L0Learn")
  f <- .rrwXy()
  expect_length(as.numeric(l0learnWeights(f$X, f$y)), f$p)
})

test_that("mrashWeights returns length-p weights and can retain the fit", {
  skip_if_not_installed("susieR")
  skip_if_not_installed("glmnet")
  f <- .rrwXy()
  w <- mrashWeights(f$X, f$y, retainFit = TRUE)
  expect_length(w, f$p)
  expect_false(is.null(attr(w, "fit")))
})

test_that("qgg Bayes-alphabet weights (N/L/A/C/R) return length-p weights", {
  skip_if_not_installed("qgg")
  f <- .rrwXy()
  mc <- list(nit = 200, nburn = 20, nthin = 1)
  expect_length(do.call(bayesNWeights, c(list(f$X, f$y), mc)), f$p)
  expect_length(do.call(bayesLWeights, c(list(f$X, f$y), mc)), f$p)
  expect_length(do.call(bayesAWeights, c(list(f$X, f$y), mc)), f$p)
  expect_length(do.call(bayesCWeights, c(list(f$X, f$y), mc)), f$p)
  expect_length(do.call(bayesRWeights, c(list(f$X, f$y), mc)), f$p)
})

test_that("bayesAlphabetWeights validates matching row counts before fitting", {
  skip_if_not_installed("qgg")
  f <- .rrwXy()
  expect_error(bayesAlphabetWeights(f$X, f$y[-1], method = "bayesN"),
               "same number of rows")
  expect_error(
    bayesAlphabetWeights(f$X, f$y, method = "bayesN", Z = matrix(1, f$n - 1, 1)),
    "same number of rows")
})

test_that("bayesBWeights / bLassoWeights (BGLR) return length-p weights", {
  skip_if_not_installed("BGLR")
  f <- .rrwXy()
  expect_length(bayesBWeights(f$X, f$y, nIter = 200, burnIn = 20, thin = 1), f$p)
  expect_length(bLassoWeights(f$X, f$y, nIter = 200, burnIn = 20, thin = 1), f$p)
})

test_that("dprVbWeights returns length-p weights and retains the fit", {
  skip_if_not_installed("RcppDPR")
  f <- .rrwXy()
  w <- dprVbWeights(f$X, f$y, retainFit = TRUE)
  expect_length(w, f$p)
  expect_false(is.null(attr(w, "fit")))
})

test_that("dprGibbsWeights returns length-p weights", {
  skip_if_not_installed("RcppDPR")
  f <- .rrwXy()
  invisible(capture.output(w <- dprGibbsWeights(f$X, f$y, sStep = 200)))
  expect_length(w, f$p)
})

test_that("dprAdaptiveGibbsWeights returns length-p weights", {
  skip_if_not_installed("RcppDPR")
  f <- .rrwXy()
  invisible(capture.output(w <- dprAdaptiveGibbsWeights(f$X, f$y, s_step = 100)))
  expect_length(w, f$p)
})

test_that("mrmashWeights fits from (X, Y) and returns p x K weights", {
  skip_if_not_installed("mr.mashr")
  skip_if_not_installed("glmnet")
  set.seed(3)
  m <- .rrwMulti(n = 60, p = 6, K = 3)
  w <- suppressMessages(mrmashWeights(X = m$X, Y = m$Y, canonicalPriorMatrices = TRUE))
  expect_equal(dim(w), c(m$p, m$K))
  expect_true(all(is.finite(w)))
})


# ----------------------------- RSS solvers (C++) ----------------------------

test_that("lassosumRss returns a p x nlambda beta matrix", {
  f <- .rrwStatLd()
  out <- lassosumRss(f$stat$b, list(blk1 = f$LD), f$n)
  expect_equal(nrow(out$beta), f$p)
  expect_equal(ncol(out$beta), length(out$lambda))
  expect_length(out$conv, length(out$lambda))
})

test_that("penalizedRss traces a solution path for MCP / SCAD / L0", {
  f <- .rrwStatLd()
  for (pen in c("MCP", "SCAD")) {
    out <- penalizedRss(f$stat$b, list(blk1 = f$LD), f$n, penalty = pen)
    expect_equal(nrow(out$beta), f$p)
  }
  outL0 <- penalizedRss(f$stat$b, list(blk1 = f$LD), f$n,
                        penalty = "L0", lambda0 = 0.01, lambda = c(0))
  expect_equal(nrow(outL0$beta), f$p)
})

test_that("prsCs returns posterior betaEst of length p", {
  f <- .rrwStatLd()
  out <- prsCs(f$stat$b, list(blk1 = f$LD), f$n, nIter = 100, nBurnin = 20, thin = 1)
  expect_length(out$betaEst, f$p)
  expect_true(all(is.finite(out$betaEst)))
})

test_that("sdpr returns betaEst of length p", {
  f <- .rrwStatLd()
  out <- sdpr(f$stat$b, list(blk1 = f$LD), f$n,
              iter = 100, burn = 20, thin = 1, verbose = FALSE)
  expect_length(out$betaEst, f$p)
})

test_that("RSS solvers validate their LD-list / sample-size / length inputs", {
  f <- .rrwStatLd()
  expect_error(prsCs(f$stat$b, f$LD, f$n), "list of LD blocks")
  expect_error(prsCs(f$stat$b, list(blk1 = f$LD), -1), "sample size")
  expect_error(prsCs(f$stat$b[-1], list(blk1 = f$LD), f$n), "same as the sum")
  expect_error(sdpr(f$stat$b[-1], list(blk1 = f$LD), f$n), "same as the length")
  expect_error(sdpr(f$stat$b, list(blk1 = f$LD), f$n, M = 2), "at least 4")
  expect_error(lassosumRss(f$stat$b, f$LD, f$n), "list of LD blocks")
  expect_error(penalizedRss(f$stat$b, f$LD, f$n), "list of LD blocks")
})

# --------------------------- RSS weight wrappers ----------------------------

test_that("lassosumRssWeights returns length-p weights and records the selection", {
  f <- .rrwStatLd()
  w <- lassosumRssWeights(f$stat, f$LD)
  expect_length(w, f$p)
  expect_equal(unname(attr(w, "lassosum_selection")["mode"]), "ld_quadratic")
  expect_length(lassosumRssWeights(f$stat, f$LD, selection = "min_fbeta"), f$p)
})

test_that("scadRssWeights / mcpRssWeights / l0learnRssWeights return length-p weights", {
  f <- .rrwStatLd()
  expect_length(scadRssWeights(f$stat, f$LD), f$p)
  expect_length(mcpRssWeights(f$stat, f$LD), f$p)
  expect_length(l0learnRssWeights(f$stat, f$LD), f$p)
})

test_that("prsCsWeights and sdprWeights follow the (stat, LD) contract", {
  f <- .rrwStatLd()
  expect_length(prsCsWeights(f$stat, f$LD, nIter = 100, nBurnin = 20, thin = 1), f$p)
  expect_length(
    sdprWeights(f$stat, f$LD, iter = 100, burn = 20, thin = 1, verbose = FALSE), f$p)
})

test_that("mrAshRssWeights returns posterior-mean weights of length p", {
  skip_if_not_installed("susieR")
  f <- .rrwStatLd()
  w <- mrAshRssWeights(f$stat, f$LD, varY = 1, sigma2E = 1,
                       s0 = c(0, 0.01, 0.1, 0.5, 1), w0 = rep(1 / 5, 5))
  expect_length(w, f$p)
  expect_true(all(is.finite(w)))
})


test_that("mrmashRssWeights fits mr.mash.rss and returns p x K weights", {
  skip_if_not_installed("mr.mashr")
  m <- .rrwMulti(n = 60, p = 6, K = 3)
  w <- mrmashRssWeights(m$stat, m$LD)
  expect_equal(dim(w), c(m$p, m$K))
  expect_true(all(is.finite(w)))
})

test_that("mrmashRssWeights errors on single-context stat$z", {
  skip_if_not_installed("mr.mashr")
  f <- .rrwStatLd()
  oneCol <- list(z = matrix(f$stat$z, ncol = 1), n = f$n)
  expect_error(mrmashRssWeights(oneCol, f$LD), ">= 2 columns")
})



# ------------------------------- pure helpers -------------------------------

test_that(".lassosumCorFromStat reads cor / z / b and validates length", {
  f <- .rrwStatLd()
  expect_length(pecotmr:::.lassosumCorFromStat(list(cor = f$stat$cor), n = f$n, p = f$p), f$p)
  expect_equal(pecotmr:::.lassosumCorFromStat(list(z = f$stat$z), n = f$n, p = f$p),
               as.numeric(f$stat$z) / sqrt(f$n))
  expect_equal(pecotmr:::.lassosumCorFromStat(list(b = f$stat$b), n = f$n, p = f$p),
               as.numeric(f$stat$b))
  expect_error(pecotmr:::.lassosumCorFromStat(list(), n = f$n, p = f$p), "one of")
  expect_error(pecotmr:::.lassosumCorFromStat(list(z = f$stat$z[-1]), n = f$n, p = f$p),
               "must equal")
})

test_that(".lassosumClampCor scales values with |cor| >= 1 below 1", {
  expect_equal(pecotmr:::.lassosumClampCor(c(0.1, 0.5)), c(0.1, 0.5))
  expect_lt(max(abs(pecotmr:::.lassosumClampCor(c(0.5, 1.5, -2)))), 1)
})

test_that(".lassosumFirstMax returns the first index of the maximum", {
  expect_equal(pecotmr:::.lassosumFirstMax(c(1, 3, 3, 2)), 2L)
  expect_equal(pecotmr:::.lassosumFirstMax(c(5, 1, 2)), 1L)
})

test_that(".lassosumSelectMinFbeta picks the minimum-fbeta candidate", {
  set.seed(7)
  cb <- matrix(rnorm(6 * 4), 6, 4)
  r <- pecotmr:::.lassosumSelectMinFbeta(cb, data.frame(fbeta = c(3, 1, 2, 4)))
  expect_equal(r$index, 2)
  expect_equal(r$mode, "min_fbeta")
  expect_equal(r$beta, cb[, 2])
})

test_that(".lassosumSelectLdQuadratic scores candidates by c'b / sqrt(b'Rb)", {
  f <- .rrwStatLd()
  set.seed(8)
  cb <- matrix(rnorm(f$p * 4), f$p, 4)
  r <- pecotmr:::.lassosumSelectLdQuadratic(cb, f$stat$b, f$LD)
  expect_equal(r$mode, "ld_quadratic")
  expect_true(r$index %in% seq_len(4))
  expect_equal(r$beta, cb[, r$index])
})

test_that("computeCovDiag returns a diagonal condition covariance", {
  m <- .rrwMulti(n = 60, p = 6, K = 3)
  cv <- computeCovDiag(m$Y)
  expect_equal(dim(cv), c(m$K, m$K))
  expect_equal(cv[upper.tri(cv)], rep(0, sum(upper.tri(cv))))
  expect_equal(unname(diag(cv)), unname(apply(m$Y, 2, var)))
})

test_that("computeCovFlash returns a finite K x K covariance from FLASH", {
  skip_if_not_installed("flashier")
  skip_if_not_installed("ebnm")
  m <- .rrwMulti(n = 80, p = 6, K = 3)
  cv <- computeCovFlash(m$Y)
  expect_equal(dim(cv), c(m$K, m$K))
  expect_true(all(is.finite(cv)))
})

test_that("buildMrmashPriorMatrices builds an expanded S0 list and a prior grid", {
  skip_if_not_installed("mr.mashr")
  set.seed(9)
  res <- buildMrmashPriorMatrices(
    Bhat = matrix(rnorm(18), 6, 3), Shat = matrix(0.2, 6, 3), K = 3)
  expect_true(is.list(res$S0))
  expect_gt(length(res$S0), 1)
  expect_true(is.numeric(res$priorGrid))
  expect_true(all(vapply(res$S0, function(s) all(dim(s) == c(3, 3)), logical(1))))
})

test_that("buildMrmashPriorMatrices errors without canonical or data-driven priors", {
  skip_if_not_installed("mr.mashr")
  expect_error(
    buildMrmashPriorMatrices(Bhat = matrix(rnorm(6), 3, 2), Shat = matrix(0.2, 3, 2),
                             K = 2, canonicalPriorMatrices = FALSE),
    "dataDrivenPriorMatrices")
})


# ===========================================================================
# mr.mash weight tests (relocated from test_rrMrmashMvsusie.R)
# ===========================================================================

# ---- mrmashWeights ----
test_that("mrmashWeights errors when mr.mashr package is not available", {
  skip_if(requireNamespace("mr.mashr", quietly = TRUE),
          "mr.mashr is installed; skipping missing-package test")

  expect_error(
    mrmashWeights(mrmashFit = NULL, X = matrix(1, 10, 5), Y = matrix(1, 10, 3)),
    "mr\\.mash\\.alpha"
  )
})

test_that("mrmashWeights errors when X and Y are NULL and fit is NULL", {
  skip_if_not(requireNamespace("mr.mashr", quietly = TRUE),
              "mr.mashr not installed")
  expect_error(mrmashWeights(mrmashFit = NULL, X = NULL, Y = NULL),
               "Both X and Y must be provided")
})

test_that("mrmashWeights(retainFit=TRUE) attaches {dataDrivenPriorMatrices, w0, V}", {
  skip_if_not(requireNamespace("mr.mashr", quietly = TRUE),
              "mr.mashr not installed")
  # These are exactly the parts fineMappingPipeline needs to rebuild the
  # mvSuSiE reweighted mixture prior (w0 -> rescaleCovW0, original $U) and the
  # residual variance (V); the heavy mu1 coefficient matrix is not retained.
  ddpm    <- list(U = list(comp = diag(2)))
  fakeFit <- structure(
    list(w0 = c(null = 0.4, comp_grid1 = 0.6), V = diag(2) * 2),
    class = "mr.mash")
  fakeCoef <- matrix(0.1, nrow = 5, ncol = 2)
  local_mocked_bindings(coef.mr.mash = function(object, ...) fakeCoef,
                        .package = "mr.mashr")
  w <- mrmashWeights(mrmashFit = fakeFit,
                     dataDrivenPriorMatrices = ddpm, retainFit = TRUE)
  fit <- attr(w, "fit")
  expect_true(is.list(fit))
  expect_identical(fit$dataDrivenPriorMatrices, ddpm)
  expect_identical(fit$w0, fakeFit$w0)
  expect_identical(fit$V,  fakeFit$V)
  # Default (retainFit = FALSE) leaves the weights free of the fit attribute.
  expect_null(attr(
    mrmashWeights(mrmashFit = fakeFit, dataDrivenPriorMatrices = ddpm), "fit"))
})
