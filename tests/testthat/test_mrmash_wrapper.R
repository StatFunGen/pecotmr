context("mrmash_wrapper")

test_that("filter_mixture_components filters zero matrices", {
  U <- list(
    mat1 = matrix(c(1, 0.5, 0.5, 1), 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
    mat2 = matrix(0, 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
    mat3 = matrix(c(0.8, 0.3, 0.3, 0.9), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  )
  w <- c(mat1 = 0.5, mat2 = 0.3, mat3 = 0.2)
  conditions_to_keep <- c("A", "B")

  result <- filter_mixture_components(conditions_to_keep, U, w)

  # mat2 should be removed (all zeros)
  expect_true(!"mat2" %in% names(result$U))
  # weights should be rescaled to maintain sum
  expect_equal(sum(result$w), sum(w), tolerance = 1e-10)
})

test_that("filter_mixture_components removes low weight components", {
  U <- list(
    mat1 = matrix(c(1, 0.5, 0.5, 1), 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
    mat2 = matrix(c(0.8, 0.3, 0.3, 0.9), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  )
  w <- c(mat1 = 0.999, mat2 = 0.00001)  # mat2 below default cutoff
  conditions_to_keep <- c("A", "B")

  result <- filter_mixture_components(conditions_to_keep, U, w, w_cutoff = 1e-04)
  expect_true(!"mat2" %in% names(result$U))
})

test_that("filter_mixture_components errors on missing condition", {
  U <- list(
    mat1 = matrix(c(1, 0.5, 0.5, 1), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  )
  w <- c(mat1 = 1.0)
  expect_error(filter_mixture_components(c("A", "C"), U, w), "not found in matrix")
})

test_that("filter_mixture_components subsets conditions", {
  U <- list(
    mat1 = matrix(c(1, 0.5, 0.2, 0.5, 1, 0.3, 0.2, 0.3, 1), 3, 3,
                  dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  )
  w <- c(mat1 = 1.0)

  result <- filter_mixture_components(c("A", "B"), U, w)
  expect_equal(nrow(result$U[[1]]), 2)
  expect_equal(ncol(result$U[[1]]), 2)
})

# ===========================================================================
# Tests from test_misc_round3.R (mrmash_wrapper coverage boost)
# ===========================================================================

# =========================================================================
# mrmash_wrapper.R: compute_w0 (lines 284-298)
# =========================================================================

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
# mrmash_wrapper.R: rescale_cov_w0 (lines 300-329)
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
# mrmash_wrapper.R: compute_grid, grid_min, grid_max, autoselect_mixsd
# (lines 333-372)
# =========================================================================

test_that("grid_max returns scaled grid_min when bhat^2 <= sbhat^2", {
  bhat <- c(0.1, 0.2)
  sbhat <- c(1.0, 1.0)
  result <- pecotmr:::grid_max(bhat, sbhat)
  expect_equal(result, 8 * pecotmr:::grid_min(bhat, sbhat))
})

test_that("grid_max returns 2*sqrt(max(bhat^2 - sbhat^2)) otherwise", {
  bhat <- c(5, 1)
  sbhat <- c(0.5, 0.5)
  expected <- 2 * sqrt(max(bhat^2 - sbhat^2))
  result <- pecotmr:::grid_max(bhat, sbhat)
  expect_equal(result, expected)
})

test_that("autoselect_mixsd returns 2-element vector when mult == 0", {
  result <- pecotmr:::autoselect_mixsd(0.01, 1.0, mult = 0)
  expect_equal(result, c(0, 0.5))
})

test_that("autoselect_mixsd returns valid grid with sqrt(2) mult", {
  result <- pecotmr:::autoselect_mixsd(0.01, 1.0, mult = sqrt(2))
  expect_true(length(result) > 1)
  expect_equal(result[length(result)], 1.0)  # last element is gmax
  # All elements should be positive

  expect_true(all(result > 0))
})

test_that("compute_grid produces a valid grid from summary statistics", {
  set.seed(42)
  bhat <- matrix(rnorm(20, sd = 2), nrow = 10, ncol = 2)
  sbhat <- matrix(abs(rnorm(20, mean = 0.5, sd = 0.1)), nrow = 10, ncol = 2)
  result <- pecotmr:::compute_grid(bhat, sbhat)
  expect_true(is.numeric(result))
  expect_true(length(result) > 0)
  expect_true(all(result > 0))
})

test_that("compute_grid handles NA and zero sbhat values", {
  bhat <- c(1, 2, NA, 4, 5)
  sbhat <- c(0.5, 0, NA, 0.3, 0.8)
  result <- pecotmr:::compute_grid(bhat, sbhat)
  expect_true(is.numeric(result))
  expect_true(length(result) > 0)
})

# =========================================================================
# mrmash_wrapper.R: mrmash_wrapper input validation (lines 99-130)
# =========================================================================

# Note: Cannot mock requireNamespace via local_mocked_bindings because it is a
# base R function, not in pecotmr's namespace. We skip these tests and instead
# test the downstream validation that we CAN exercise.

test_that("mrmash_wrapper errors when X and Y are not matrices", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("mr.mashr")
  expect_error(mrmash_wrapper(data.frame(x = 1:3), matrix(1:6, nrow = 3, ncol = 2)),
               "matrices")
})

test_that("mrmash_wrapper errors when X and Y row counts differ", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("mr.mashr")
  expect_error(mrmash_wrapper(matrix(1:6, nrow = 3, ncol = 2), matrix(1:8, nrow = 4, ncol = 2)),
               "same number of rows")
})

test_that("mrmash_wrapper errors when prior_grid is not a vector", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("mr.mashr")
  X <- matrix(rnorm(12), nrow = 3, ncol = 4)
  Y <- matrix(rnorm(6), nrow = 3, ncol = 2)
  expect_error(mrmash_wrapper(X, Y, prior_grid = matrix(1:4, nrow = 2)),
               "prior_grid must be a vector")
})

test_that("mrmash_wrapper errors when no prior matrices and canonical_prior_matrices is FALSE", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("mr.mashr")
  X <- matrix(rnorm(12), nrow = 3, ncol = 4)
  Y <- matrix(rnorm(6), nrow = 3, ncol = 2)
  expect_error(mrmash_wrapper(X, Y, data_driven_prior_matrices = NULL,
                               canonical_prior_matrices = FALSE),
               "data_driven_prior_matrices")
})

test_that("mrmash_wrapper warns when Y has missing and B_init_method is glasso", {
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
      mrmash_wrapper(X, Y, B_init_method = "glasso",
                     data_driven_prior_matrices = list(U = list(matrix(1, 2, 2))),
                     canonical_prior_matrices = FALSE),
      error = function(e) NULL
    ),
    "glasso"
  )
})

# =========================================================================
# mrmash_wrapper.R: compute_coefficients_glasso (lines 211-240)
# =========================================================================

test_that("compute_coefficients_glasso runs without Xnew", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 5; r <- 3
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
  colnames(Y) <- paste0("cond", 1:r)
  result <- pecotmr:::compute_coefficients_glasso(X, Y, standardize = FALSE,
                                                   nthreads = 1, Xnew = NULL)
  expect_true("Bhat" %in% names(result))
  expect_true("Ytrain" %in% names(result))
  expect_equal(nrow(result$Bhat), p)
  expect_equal(ncol(result$Bhat), r)
  expect_null(result$Yhat_new)
})

test_that("compute_coefficients_glasso runs with Xnew", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 5; r <- 3
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
  colnames(Y) <- paste0("cond", 1:r)
  Xnew <- matrix(rnorm(10 * p), nrow = 10, ncol = p)
  result <- pecotmr:::compute_coefficients_glasso(X, Y, standardize = FALSE,
                                                   nthreads = 1, Xnew = Xnew)
  expect_true("Yhat_new" %in% names(result))
  expect_equal(nrow(result$Yhat_new), 10)
  expect_equal(ncol(result$Yhat_new), r)
  expect_equal(colnames(result$Yhat_new), colnames(Y))
})

# =========================================================================
# mrmash_wrapper.R: compute_coefficients_univ_glmnet (lines 243-281)
# =========================================================================

test_that("compute_coefficients_univ_glmnet runs without Xnew", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 60; p <- 5; r <- 2
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
  colnames(Y) <- paste0("cond", 1:r)
  result <- pecotmr:::compute_coefficients_univ_glmnet(X, Y, alpha = 0.5,
                                                        standardize = FALSE,
                                                        nthreads = 1, Xnew = NULL)
  expect_true("Bhat" %in% names(result))
  expect_true("intercept" %in% names(result))
  expect_equal(nrow(result$Bhat), p)
  expect_equal(ncol(result$Bhat), r)
  expect_null(result$Yhat_new)
})

test_that("compute_coefficients_univ_glmnet runs with Xnew", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 60; p <- 5; r <- 2
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
  colnames(Y) <- paste0("cond", 1:r)
  Xnew <- matrix(rnorm(8 * p), nrow = 8, ncol = p)
  result <- pecotmr:::compute_coefficients_univ_glmnet(X, Y, alpha = 0.5,
                                                        standardize = FALSE,
                                                        nthreads = 1, Xnew = Xnew)
  expect_true("Yhat_new" %in% names(result))
  expect_equal(nrow(result$Yhat_new), 8)
  expect_equal(ncol(result$Yhat_new), r)
  expect_equal(colnames(result$Yhat_new), colnames(Y))
})

test_that("compute_coefficients_univ_glmnet handles NA in Y", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 60; p <- 5; r <- 2
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
  colnames(Y) <- paste0("cond", 1:r)
  Y[1:5, 1] <- NA  # introduce missing values in one condition
  result <- pecotmr:::compute_coefficients_univ_glmnet(X, Y, alpha = 0.5,
                                                        standardize = FALSE,
                                                        nthreads = 1, Xnew = NULL)
  expect_true("Bhat" %in% names(result))
  expect_equal(nrow(result$Bhat), p)
})

# =========================================================================
# mrmash_wrapper.R: mrmash_wrapper seed warning (line 107-108)
# =========================================================================

# Note: Cannot mock base::exists() via local_mocked_bindings.
# The seed-check message on line 107-108 would require removing .Random.seed
# from the global environment, which is not safe to do in tests.
