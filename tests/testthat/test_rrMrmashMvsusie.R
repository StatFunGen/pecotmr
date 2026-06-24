context("regularized_regression - mrmash / mvsusie / fsusie")

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

# ---- mvsusieWeights ----
test_that("mvsusieWeights errors when mvsusieR package is not available", {
  skip_if(requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR is installed; skipping missing-package test")

  expect_error(
    mvsusieWeights(mvsusieFit = NULL, X = matrix(1, 10, 5), Y = matrix(1, 10, 3)),
    "mvsusieR"
  )
})

test_that("mvsusieWeights errors when X and Y are NULL and fit is NULL", {
  skip_if_not(requireNamespace("mvsusieR", quietly = TRUE),
              "mvsusieR not installed")
  expect_error(mvsusieWeights(mvsusieFit = NULL, X = NULL, Y = NULL),
               "Both X and Y must be provided")
})

test_that("mvsusieWeights fits model and returns coefficients when fit is NULL", {
  skip_if_not(requireNamespace("mvsusieR", quietly = TRUE),
              "mvsusieR not installed")
  set.seed(42)
  n <- 30
  p <- 5
  R <- 3
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * R), n, R)
  fake_coef <- matrix(rnorm((p + 1) * R), nrow = p + 1, ncol = R)
  captured <- list()

  local_mocked_bindings(
    create_mixture_prior = function(...) list(),
    mvsusie = function(...) {
      captured <<- list(...)
      "mock_fit"
    },
    coef.mvsusie = function(...) fake_coef,
    .package = "mvsusieR"
  )

  result <- expect_message(
    mvsusieWeights(X = X, Y = Y, L = 12, LGreedy = 4),
    "mvsusieFit is not provided"
  )
  # Should return coef without intercept row
  expect_equal(dim(result), c(p, R))
  expect_equal(result, fake_coef[-1, ])
  expect_equal(captured$L, 12)
  expect_equal(captured$L_greedy, 4)
})

test_that("mvsusieWeights returns coefficients from provided fit", {
  skip_if_not(requireNamespace("mvsusieR", quietly = TRUE),
              "mvsusieR not installed")
  p <- 5
  R <- 3
  fake_coef <- matrix(rnorm((p + 1) * R), nrow = p + 1, ncol = R)

  local_mocked_bindings(
    coef.mvsusie = function(...) fake_coef,
    .package = "mvsusieR"
  )

  result <- mvsusieWeights(mvsusieFit = "precomputed_fit")
  expect_equal(dim(result), c(p, R))
  expect_equal(result, fake_coef[-1, ])
})

# ---- fsusieWeights ----
# Collapse a functional SuSiE fit to a variants x features TWAS weight matrix
# (the all-SNP wavelet posterior mean, the coef.susie analog).

.fw_makeFsusieFit <- function(seed = 1, n = 150L, p = 24L, J = 16L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)), paste0("v", seq_len(p))))
  b1 <- sin(seq(0, 2 * pi, length.out = J))
  b2 <- cos(seq(0, pi, length.out = J))
  Y <- X[, 3] %o% b1 + X[, 10] %o% b2 +
    matrix(rnorm(n * J, sd = 0.3), n, J)
  colnames(Y) <- paste0("f", seq_len(J))
  list(X = X, Y = Y,
       fit = suppressWarnings(fsusieR::susiF(
         X = X, Y = Y, pos = seq_len(J), L = 5,
         post_processing = "none", verbose = FALSE)))
}

test_that("fsusieWeights returns a variants x features matrix with variant rownames", {
  skip_if_not_installed("fsusieR")
  skip_if_not_installed("wavethresh")
  obj <- .fw_makeFsusieFit()
  W <- fsusieWeights(fsusieFit = obj$fit, variantIds = colnames(obj$X))
  expect_true(is.matrix(W))
  expect_equal(nrow(W), ncol(obj$X))
  expect_equal(ncol(W), ncol(obj$Y))
  expect_equal(rownames(W), colnames(obj$X))
})

test_that("fsusieWeights matches fsusieR's own out_prep reconstruction (post_processing='none')", {
  skip_if_not_installed("fsusieR")
  skip_if_not_installed("wavethresh")
  obj <- .fw_makeFsusieFit()
  fit <- obj$fit
  # The alpha-weighted sum over SNPs of the per-SNP feature-domain curves that
  # fsusieWeights reconstructs must equal fSuSiE's own fitted_func[[l]] (built
  # by out_prep.susiF) for every effect l.
  csdX <- as.numeric(fit$csd_X)
  perScale <- "mixture_normal_per_scale" %in% class(fsusieR::get_G_prior(fit))
  indxLst <- fsusieR::gen_wavelet_indx(log2(length(fit$outing_grid)))
  scaleCols <- if (perScale) indxLst[[length(indxLst)]]
               else ncol(as.matrix(fit$fitted_wc[[1L]]))
  S <- pecotmr:::.fsusieSynthesisMatrix(fit$n_wac, scaleCols)
  maxErr <- 0
  for (l in seq_along(fit$fitted_wc)) {
    al <- as.numeric(fit$alpha[[l]])
    contrib <- colSums((al * (1 / csdX) * as.matrix(fit$fitted_wc[[l]])) %*% S)
    maxErr <- max(maxErr, max(abs(contrib - as.numeric(fit$fitted_func[[l]]))))
  }
  expect_lt(maxErr, 1e-8)
})

test_that("fsusieWeights concentrates weight on the causal SNPs", {
  skip_if_not_installed("fsusieR")
  skip_if_not_installed("wavethresh")
  obj <- .fw_makeFsusieFit()
  W <- fsusieWeights(fsusieFit = obj$fit, variantIds = colnames(obj$X))
  rowNorm <- sqrt(rowSums(W^2))
  top2 <- names(sort(rowNorm, decreasing = TRUE))[1:2]
  expect_setequal(top2, c("v3", "v10"))
})

test_that("fsusieWeights fast path returns precomputed $coef for a trimmed fit", {
  # A trimmed fSuSiE fit drops fitted_wc but keeps the precomputed weight
  # matrix in $coef; fsusieWeights returns it without touching wavelet slots.
  W0 <- matrix(c(1, 0, 2, 0, 0, 3), nrow = 3,
               dimnames = list(c("v1", "v2", "v3"), c("f1", "f2")))
  trimmed <- list(coef = W0, pip = c(0.1, 0.2, 0.7))
  class(trimmed) <- c("fsusie", "susie")
  W <- fsusieWeights(fsusieFit = trimmed)
  expect_identical(W, W0)
})

test_that("fsusieWeights errors without a fit and on an unusable (trimmed, no coef) fit", {
  expect_error(fsusieWeights(fsusieFit = NULL), "is required")
  bad <- list(pip = c(0.1, 0.9))  # no coef, no fitted_wc
  class(bad) <- c("fsusie", "susie")
  expect_error(fsusieWeights(fsusieFit = bad), "missing required slot")
})
