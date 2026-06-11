context("ensembleWeights")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
# Build a synthetic twasWeightsCv() output with K methods. Each method's
# prediction is a convex combination of the truth + noise, letting us control
# per-method accuracy. Returns a list shaped exactly like twasWeightsCv()'s
# output (with $prediction, $performance, $sample_partition).
make_cv_result <- function(n = 100, K = 4, seed = 1, method_quality = NULL) {
  set.seed(seed)
  y <- rnorm(n)
  sample_names <- paste0("sample_", seq_len(n))

  if (is.null(method_quality)) {
    # Methods with decreasing quality (noise amounts)
    method_quality <- seq(0.1, 0.9, length.out = K)
  }
  stopifnot(length(method_quality) == K)

  method_names <- paste0("method", seq_len(K))
  pred_names <- paste0(method_names, "_predicted")

  prediction <- setNames(lapply(seq_len(K), function(k) {
    noise_sd <- method_quality[k]
    pred <- y + rnorm(n, sd = noise_sd)
    mat <- matrix(pred, ncol = 1)
    rownames(mat) <- sample_names
    colnames(mat) <- "outcome_1"
    mat
  }), pred_names)

  # Dummy performance (not used by ensembleWeights)
  performance <- setNames(lapply(seq_len(K), function(k) {
    m <- matrix(NA, nrow = 1, ncol = 6)
    colnames(m) <- c("corr", "rsq", "adj_rsq", "pval", "RMSE", "MAE")
    m
  }), paste0(method_names, "_performance"))

  list(
    samplePartition = data.frame(Sample = sample_names,
                                   Fold = rep(1:5, length.out = n),
                                   stringsAsFactors = FALSE),
    prediction = prediction,
    performance = performance,
    time_elapsed = 0,
    .y = y,
    .method_names = method_names
  )
}

# Build synthetic twasWeights() output
make_weight_list <- function(p = 20, method_names, seed = 2) {
  set.seed(seed)
  setNames(lapply(method_names, function(m) {
    w <- matrix(rnorm(p), ncol = 1)
    rownames(w) <- paste0("var_", seq_len(p))
    colnames(w) <- "outcome_1"
    w
  }), paste0(method_names, "_weights"))
}

# ===========================================================================
#  Input validation
# ===========================================================================

test_that("ensembleWeights: NULL cv_results errors", {
  expect_error(ensembleWeights(NULL, Y = rnorm(10)), "cvResults")
})

test_that("ensembleWeights: NULL Y errors", {
  cv <- make_cv_result(n = 20, K = 3)
  expect_error(ensembleWeights(cv, Y = NULL), "'Y' is required")
})

test_that("ensembleWeights: single method errors (need >= 2 for ensemble)", {
  cv <- make_cv_result(n = 20, K = 1)
  expect_error(ensembleWeights(cv, Y = cv$.y),
               "at least 2 methods")
})

test_that("ensembleWeights: invalid context_index errors", {
  cv <- make_cv_result(n = 20, K = 3)
  expect_error(ensembleWeights(cv, Y = cv$.y, contextIndex = 0),
               "contextIndex")
  expect_error(ensembleWeights(cv, Y = cv$.y, contextIndex = "a"),
               "contextIndex")
})

test_that("ensembleWeights: context_index beyond Y columns errors", {
  cv <- make_cv_result(n = 20, K = 3)
  Y_mat <- matrix(cv$.y, ncol = 1)
  expect_error(ensembleWeights(cv, Y = Y_mat, contextIndex = 5),
               "contextIndex")
})

test_that("ensembleWeights: multi-dataset with mismatched lengths errors", {
  cv1 <- make_cv_result(n = 20, K = 3, seed = 1)
  cv2 <- make_cv_result(n = 20, K = 3, seed = 2)
  expect_error(ensembleWeights(list(cv1, cv2), Y = list(cv1$.y)),
               "same length")
})

test_that("ensembleWeights: multi-dataset with different methods errors", {
  cv1 <- make_cv_result(n = 20, K = 3, seed = 1)
  cv2 <- make_cv_result(n = 20, K = 4, seed = 2)
  expect_error(
    ensembleWeights(list(cv1, cv2), Y = list(cv1$.y, cv2$.y)),
    "same method names"
  )
})

# ===========================================================================
#  Core algorithm correctness
# ===========================================================================

test_that("ensembleWeights: coefficients are non-negative and sum to 1", {
  cv <- make_cv_result(n = 100, K = 4, seed = 42)
  res <- ensembleWeights(cv, Y = cv$.y)

  expect_true(all(res$method_coef >= 0))
  expect_equal(sum(res$method_coef), 1, tolerance = 1e-6)
})

test_that("ensembleWeights: best method receives the largest coefficient", {
  # Method 1 is best (lowest noise), method K is worst
  cv <- make_cv_result(n = 200, K = 4, seed = 7,
                        method_quality = c(0.1, 0.5, 0.8, 1.2))
  res <- ensembleWeights(cv, Y = cv$.y)

  expect_equal(names(which.max(res$method_coef)), "method1")
})

test_that("ensembleWeights: does not return ensemble_performance (in-sample R^2 omitted)", {
  cv <- make_cv_result(n = 300, K = 5, seed = 13)
  res <- ensembleWeights(cv, Y = cv$.y)

  expect_null(res$ensemble_performance)
  expect_false("ensemble_performance" %in% names(res))
})

test_that("ensembleWeights: per-method R^2 values are sensible (between 0 and 1)", {
  cv <- make_cv_result(n = 200, K = 4, seed = 21)
  res <- ensembleWeights(cv, Y = cv$.y)

  expect_true(all(res$method_performance >= 0, na.rm = TRUE))
  expect_true(all(res$method_performance <= 1, na.rm = TRUE))
  expect_equal(length(res$method_performance), 4)
})

test_that("ensembleWeights: method names are stripped of _predicted suffix", {
  cv <- make_cv_result(n = 50, K = 3, seed = 1)
  res <- ensembleWeights(cv, Y = cv$.y)

  expect_equal(names(res$method_coef),
               c("method1", "method2", "method3"))
  expect_equal(names(res$method_performance),
               c("method1", "method2", "method3"))
})

# ===========================================================================
#  Sample name alignment
# ===========================================================================

test_that("ensembleWeights: aligns Y and predictions by sample name", {
  cv <- make_cv_result(n = 50, K = 3, seed = 10)

  # Shuffle Y order relative to predictions
  shuffled_order <- sample(50)
  y_shuffled <- cv$.y[shuffled_order]
  names(y_shuffled) <- paste0("sample_", shuffled_order)

  res_aligned <- ensembleWeights(cv, Y = y_shuffled)
  res_original <- ensembleWeights(cv, Y = cv$.y)

  # Results should be identical regardless of Y order
  expect_equal(res_aligned$method_coef, res_original$method_coef, tolerance = 1e-10)
})

test_that("ensembleWeights: aligns Y matrix and predictions by sample name", {
  cv <- make_cv_result(n = 50, K = 3, seed = 10)

  # Create Y as a matrix with shuffled row order
  shuffled_order <- sample(50)
  Y_mat <- matrix(cv$.y[shuffled_order], ncol = 1)
  rownames(Y_mat) <- paste0("sample_", shuffled_order)

  res_aligned <- ensembleWeights(cv, Y = Y_mat)
  res_original <- ensembleWeights(cv, Y = cv$.y)

  expect_equal(res_aligned$method_coef, res_original$method_coef, tolerance = 1e-10)
})

test_that("ensembleWeights: errors when no common sample names", {
  cv <- make_cv_result(n = 20, K = 3, seed = 1)
  y_bad <- setNames(rnorm(20), paste0("other_", seq_len(20)))

  expect_error(ensembleWeights(cv, Y = y_bad), "No common sample names")
})

# ===========================================================================
#  Zero-variance / edge cases
# ===========================================================================

test_that("ensembleWeights: zero-variance method gets coefficient 0", {
  cv <- make_cv_result(n = 100, K = 3, seed = 5)
  # Force method 2 to have constant predictions
  cv$prediction$method2_predicted[, 1] <- 0.5
  res <- ensembleWeights(cv, Y = cv$.y)

  expect_equal(res$method_coef["method2"], c(method2 = 0))
  expect_equal(sum(res$method_coef), 1, tolerance = 1e-6)
})

test_that("ensembleWeights: NA predictions in some samples are dropped", {
  cv <- make_cv_result(n = 100, K = 3, seed = 5)
  cv$prediction$method1_predicted[1:5, 1] <- NA
  expect_message(
    res <- ensembleWeights(cv, Y = cv$.y),
    "Dropping"
  )
  expect_equal(sum(res$method_coef), 1, tolerance = 1e-6)
})

test_that("ensembleWeights: all zero-variance methods errors", {
  cv <- make_cv_result(n = 50, K = 2, seed = 5)
  cv$prediction$method1_predicted[, 1] <- 0
  cv$prediction$method2_predicted[, 1] <- 0
  expect_error(ensembleWeights(cv, Y = cv$.y),
               "zero-variance predictions")
})

# ===========================================================================
#  Weight combination
# ===========================================================================

test_that("ensembleWeights: ensemble_twas_weights is sum of zeta_k * w_k", {
  cv <- make_cv_result(n = 100, K = 3, seed = 42)
  wt <- make_weight_list(p = 10, method_names = cv$.method_names)

  res <- ensembleWeights(cv, Y = cv$.y, twasWeightList = wt)

  expect_false(is.null(res$ensemble_twas_weights))

  # Verify the combination is correct
  expected <- matrix(0, nrow = 10, ncol = 1)
  for (k in seq_along(cv$.method_names)) {
    m <- cv$.method_names[k]
    expected <- expected + res$method_coef[m] * wt[[paste0(m, "_weights")]]
  }
  expect_equal(as.numeric(res$ensemble_twas_weights),
               as.numeric(expected),
               tolerance = 1e-10)
})

test_that("ensembleWeights: NULL twas_weight_list returns NULL ensemble_twas_weights", {
  cv <- make_cv_result(n = 50, K = 3, seed = 1)
  res <- ensembleWeights(cv, Y = cv$.y, twasWeightList = NULL)
  expect_null(res$ensemble_twas_weights)
})

test_that("ensembleWeights: weights with no matching keys warns and skips", {
  cv <- make_cv_result(n = 50, K = 2, seed = 1)
  wt <- list(unknown_weights = matrix(1, nrow = 10, ncol = 1))

  expect_warning(
    res <- ensembleWeights(cv, Y = cv$.y, twasWeightList = wt),
    "No matching weight keys"
  )
  expect_null(res$ensemble_twas_weights)
})

# ===========================================================================
#  Multi-dataset ensemble
# ===========================================================================

test_that("ensembleWeights: multi-dataset combines predictions correctly", {
  cv1 <- make_cv_result(n = 80, K = 3, seed = 1)
  cv2 <- make_cv_result(n = 80, K = 3, seed = 2)

  res <- ensembleWeights(
    cvResults = list(cv1, cv2),
    Y = list(cv1$.y, cv2$.y)
  )

  expect_true(all(res$method_coef >= 0))
  expect_equal(sum(res$method_coef), 1, tolerance = 1e-6)
  expect_equal(length(res$method_performance), 3)
})

test_that("ensembleWeights: Y as matrix with context_index works", {
  cv <- make_cv_result(n = 50, K = 3, seed = 1)
  Y_mat <- matrix(cv$.y, ncol = 1)
  colnames(Y_mat) <- "ctx1"

  res <- ensembleWeights(cv, Y = Y_mat, contextIndex = 1)
  expect_equal(sum(res$method_coef), 1, tolerance = 1e-6)
})

# ===========================================================================
#  End-to-end with twasWeightsCv (integration)
# ===========================================================================

test_that("ensembleWeights: end-to-end with twasWeightsCv output", {
  skip_if_not_installed("glmnet")

  set.seed(42)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("var_", seq_len(p))
  rownames(X) <- paste0("sample_", seq_len(n))

  beta <- c(1.5, -1.0, 0.8, rep(0, p - 3))
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.5))

  cv <- suppressMessages(twasWeightsCv(
    X, y, fold = 3,
    weightMethods = list(
      lassoWeights = list(),
      enetWeights = list()
    )
  ))

  res <- ensembleWeights(cv, Y = y)

  expect_equal(sum(res$method_coef), 1, tolerance = 1e-6)
  expect_true(all(res$method_coef >= 0))
  expect_equal(names(res$method_coef), c("lasso", "enet"))
  expect_null(res$ensemble_performance)
})

# ===========================================================================
#  twasWeightsPipeline ensemble integration
# ===========================================================================

test_that("pipeline: ensemble=TRUE with only 1 method prints skip message", {
  skip_if_not_installed("glmnet")

  set.seed(42)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("var_", seq_len(p))
  rownames(X) <- paste0("sample_", seq_len(n))

  beta <- c(1.5, -1.0, 0.8, rep(0, p - 3))
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.5))

  msgs <- testthat::capture_messages(
    res <- twasWeightsPipeline(
      X, y, cvFolds = 3,
      weightMethods = list(lassoWeights = list()),
      ensemble = TRUE
    )
  )

  # Should see the skip message
  expect_true(any(grepl("Ensemble model skipped.*only 1 weight method provided", msgs)))

  # No ensemble result should be present
  expect_null(res$ensemble)
  expect_false("ensembleWeights" %in% getMethodNames(res$twas_weights))
})

test_that("pipeline: ensemble=TRUE skips when methods fail R^2 cutoff", {
  skip_if_not_installed("glmnet")

  set.seed(42)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("var_", seq_len(p))
  rownames(X) <- paste0("sample_", seq_len(n))

  # Use signal so methods produce non-zero weights, but set threshold very high
  beta <- c(1.5, -1.0, 0.8, rep(0, p - 3))
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.5))

  msgs <- testthat::capture_messages(
    res <- twasWeightsPipeline(
      X, y, cvFolds = 3,
      weightMethods = list(lassoWeights = list(), enetWeights = list()),
      ensemble = TRUE,
      ensembleR2Threshold = 0.99  # impossibly high threshold
    )
  )

  expect_true(any(grepl("Ensemble TWAS skipped", msgs)))
  expect_null(res$ensemble)
  expect_false("ensembleWeights" %in% getMethodNames(res$twas_weights))
})

test_that("pipeline: ensemble=TRUE succeeds and adds ensembleWeights", {
  skip_if_not_installed("glmnet")

  set.seed(42)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("var_", seq_len(p))
  rownames(X) <- paste0("sample_", seq_len(n))

  beta <- c(1.5, -1.0, 0.8, rep(0, p - 3))
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.5))

  msgs <- testthat::capture_messages(
    res <- twasWeightsPipeline(
      X, y, cvFolds = 3,
      weightMethods = list(lassoWeights = list(), enetWeights = list()),
      ensemble = TRUE
    )
  )

  expect_true(any(grepl("Computing ensemble TWAS weights", msgs)))

  # Ensemble weights added alongside individual methods
  expect_true("ensembleWeights" %in% getMethodNames(res$twas_weights))
  expect_true("lassoWeights" %in% getMethodNames(res$twas_weights))
  expect_true("enetWeights" %in% getMethodNames(res$twas_weights))

  # Ensemble predictions added
  expect_true("ensemble_predicted" %in% names(res$twas_predictions))

  # Ensemble result metadata present
  expect_false(is.null(res$ensemble))
  expect_true(all(res$ensemble$method_coef >= 0))
  expect_equal(sum(res$ensemble$method_coef), 1, tolerance = 1e-6)

  # Ensemble weights should have same length as individual weights
  expect_equal(length(getWeights(res$twas_weights,"ensembleWeights")),
               length(getWeights(res$twas_weights,"lassoWeights")))
})

test_that("pipeline: ensemble=FALSE does not run ensemble", {
  skip_if_not_installed("glmnet")

  set.seed(42)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("var_", seq_len(p))
  rownames(X) <- paste0("sample_", seq_len(n))

  beta <- c(1.5, -1.0, 0.8, rep(0, p - 3))
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.5))

  res <- suppressMessages(twasWeightsPipeline(
    X, y, cvFolds = 3,
    weightMethods = list(lassoWeights = list(), enetWeights = list()),
    ensemble = FALSE
  ))

  expect_null(res$ensemble)
  expect_false("ensembleWeights" %in% getMethodNames(res$twas_weights))
})

test_that("pipeline: ensemble_r2_threshold filters methods for ensemble", {
  skip_if_not_installed("glmnet")

  set.seed(42)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("var_", seq_len(p))
  rownames(X) <- paste0("sample_", seq_len(n))

  beta <- c(1.5, -1.0, 0.8, rep(0, p - 3))
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.5))

  # Run with very low threshold - both methods should pass
  msgs_low <- testthat::capture_messages(
    res_low <- twasWeightsPipeline(
      X, y, cvFolds = 3,
      weightMethods = list(lassoWeights = list(), enetWeights = list()),
      ensemble = TRUE,
      ensembleR2Threshold = 0.001
    )
  )
  expect_false(is.null(res_low$ensemble))

  # Run with very high threshold - neither should pass
  msgs_high <- testthat::capture_messages(
    res_high <- twasWeightsPipeline(
      X, y, cvFolds = 3,
      weightMethods = list(lassoWeights = list(), enetWeights = list()),
      ensemble = TRUE,
      ensembleR2Threshold = 0.99
    )
  )
  expect_true(any(grepl("Ensemble TWAS skipped", msgs_high)))
  expect_null(res_high$ensemble)
})

# ===========================================================================
#  Solver alternatives
# ===========================================================================

for (slv in c("quadprog", "nnls", "lbfgsb", "glmnet")) {
  test_that(paste0("ensembleWeights: solver='", slv, "' produces valid coefficients"), {
    if (slv == "quadprog") skip_if_not_installed("quadprog")
    if (slv == "nnls") skip_if_not_installed("nnls")
    if (slv == "glmnet") skip_if_not_installed("glmnet")

    cv <- make_cv_result(n = 100, K = 4, seed = 42)
    res <- ensembleWeights(cv, Y = cv$.y, solver = slv)

    expect_true(all(res$method_coef >= 0))
    expect_equal(sum(res$method_coef), 1, tolerance = 1e-6)
    expect_equal(length(res$method_coef), 4)
  })

  test_that(paste0("ensembleWeights: solver='", slv, "' assigns best method largest coef"), {
    if (slv == "quadprog") skip_if_not_installed("quadprog")
    if (slv == "nnls") skip_if_not_installed("nnls")
    if (slv == "glmnet") skip_if_not_installed("glmnet")

    cv <- make_cv_result(n = 200, K = 4, seed = 7,
                          method_quality = c(0.1, 0.5, 0.8, 1.2))
    res <- ensembleWeights(cv, Y = cv$.y, solver = slv)

    expect_equal(names(which.max(res$method_coef)), "method1")
  })

  test_that(paste0("ensembleWeights: solver='", slv, "' combines weights correctly"), {
    if (slv == "quadprog") skip_if_not_installed("quadprog")
    if (slv == "nnls") skip_if_not_installed("nnls")
    if (slv == "glmnet") skip_if_not_installed("glmnet")

    cv <- make_cv_result(n = 100, K = 3, seed = 42)
    wt <- make_weight_list(p = 10, method_names = cv$.method_names)
    res <- ensembleWeights(cv, Y = cv$.y, twasWeightList = wt, solver = slv)

    expect_false(is.null(res$ensemble_twas_weights))

    expected <- matrix(0, nrow = 10, ncol = 1)
    for (k in seq_along(cv$.method_names)) {
      m <- cv$.method_names[k]
      expected <- expected + res$method_coef[m] * wt[[paste0(m, "_weights")]]
    }
    expect_equal(as.numeric(res$ensemble_twas_weights),
                 as.numeric(expected),
                 tolerance = 1e-10)
  })
}

test_that("ensembleWeights: invalid solver errors", {
  cv <- make_cv_result(n = 50, K = 3, seed = 1)
  expect_error(ensembleWeights(cv, Y = cv$.y, solver = "bogus"),
               "arg")
})

test_that("pipeline: ensemble_solver='nnls' works end-to-end", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("nnls")

  set.seed(42)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("var_", seq_len(p))
  rownames(X) <- paste0("sample_", seq_len(n))

  beta <- c(1.5, -1.0, 0.8, rep(0, p - 3))
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.5))

  msgs <- testthat::capture_messages(
    res <- twasWeightsPipeline(
      X, y, cvFolds = 3,
      weightMethods = list(lassoWeights = list(), enetWeights = list()),
      ensemble = TRUE,
      ensembleSolver = "nnls"
    )
  )

  expect_true(any(grepl("Computing ensemble TWAS weights", msgs)))
  expect_true("ensembleWeights" %in% getMethodNames(res$twas_weights))
  expect_true(all(res$ensemble$method_coef >= 0))
  expect_equal(sum(res$ensemble$method_coef), 1, tolerance = 1e-6)
})

test_that("pipeline: ensemble_solver='lbfgsb' works end-to-end", {
  skip_if_not_installed("glmnet")

  set.seed(42)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("var_", seq_len(p))
  rownames(X) <- paste0("sample_", seq_len(n))

  beta <- c(1.5, -1.0, 0.8, rep(0, p - 3))
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.5))

  msgs <- testthat::capture_messages(
    res <- twasWeightsPipeline(
      X, y, cvFolds = 3,
      weightMethods = list(lassoWeights = list(), enetWeights = list()),
      ensemble = TRUE,
      ensembleSolver = "lbfgsb"
    )
  )

  expect_true(any(grepl("Computing ensemble TWAS weights", msgs)))
  expect_true("ensembleWeights" %in% getMethodNames(res$twas_weights))
  expect_true(all(res$ensemble$method_coef >= 0))
  expect_equal(sum(res$ensemble$method_coef), 1, tolerance = 1e-6)
})

test_that("pipeline: ensemble_solver='glmnet' works end-to-end", {
  skip_if_not_installed("glmnet")

  set.seed(42)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("var_", seq_len(p))
  rownames(X) <- paste0("sample_", seq_len(n))

  beta <- c(1.5, -1.0, 0.8, rep(0, p - 3))
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.5))

  msgs <- testthat::capture_messages(
    res <- twasWeightsPipeline(
      X, y, cvFolds = 3,
      weightMethods = list(lassoWeights = list(), enetWeights = list()),
      ensemble = TRUE,
      ensembleSolver = "glmnet"
    )
  )

  expect_true(any(grepl("Computing ensemble TWAS weights", msgs)))
  expect_true("ensembleWeights" %in% getMethodNames(res$twas_weights))
  expect_true(all(res$ensemble$method_coef >= 0))
  expect_equal(sum(res$ensemble$method_coef), 1, tolerance = 1e-6)
})

test_that("ensembleWeights: solver='glmnet' respects alpha parameter", {
  skip_if_not_installed("glmnet")

  cv <- make_cv_result(n = 200, K = 4, seed = 42)

  res_lasso <- ensembleWeights(cv, Y = cv$.y, solver = "glmnet", alpha = 1)
  res_ridge <- ensembleWeights(cv, Y = cv$.y, solver = "glmnet", alpha = 0)

  # Both should be valid
  expect_true(all(res_lasso$method_coef >= 0))
  expect_equal(sum(res_lasso$method_coef), 1, tolerance = 1e-6)
  expect_true(all(res_ridge$method_coef >= 0))
  expect_equal(sum(res_ridge$method_coef), 1, tolerance = 1e-6)

  # Lasso should be at least as sparse as ridge (fewer or equal non-zero coefs)
  n_nonzero_lasso <- sum(res_lasso$method_coef > 1e-8)
  n_nonzero_ridge <- sum(res_ridge$method_coef > 1e-8)
  expect_true(n_nonzero_lasso <= n_nonzero_ridge)
})
