context("twas_weights")

# ---------------------------------------------------------------------------
# Shared synthetic data generator
# ---------------------------------------------------------------------------
make_data <- function(n = 50, p = 10, seed = 42, add_zero_var_col = FALSE) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("var_", seq_len(p))
  rownames(X) <- paste0("sample_", seq_len(n))

  beta <- rep(0, p)
  beta[1:3] <- c(1.5, -0.8, 0.5)
  noise <- rnorm(n, sd = 0.5)
  Y <- X %*% beta + noise
  Y <- matrix(Y, ncol = 1)
  colnames(Y) <- "outcome_1"
  rownames(Y) <- rownames(X)

  if (add_zero_var_col) {
    # Append a constant column (zero variance)
    X <- cbind(X, zero_var = rep(7, n))
    colnames(X)[p + 1] <- "zero_var"
  }

  list(X = X, Y = Y, beta = beta)
}

# ===========================================================================
#
#  twas_predict
#
# ===========================================================================

test_that("twas_predict: basic matrix multiplication is correct", {
  d <- make_data(n = 20, p = 5)
  set.seed(99)
  w <- matrix(runif(5), ncol = 1)
  rownames(w) <- colnames(d$X)
  wl <- list(test_weights = w)
  res <- twas_predict(d$X, wl)

  expected <- d$X %*% w
  expect_equal(res[["test_predicted"]], expected)
})

test_that("twas_predict: multiple weight methods in list", {
  d <- make_data(n = 20, p = 5)
  w1 <- matrix(c(1, 0, 0, 0, 0), ncol = 1)
  w2 <- matrix(c(0, 0, 0, 0, 1), ncol = 1)
  wl <- list(method_a_weights = w1, method_b_weights = w2)

  res <- twas_predict(d$X, wl)

  expect_length(res, 2)
  expect_equal(res[["method_a_predicted"]], d$X %*% w1)
  expect_equal(res[["method_b_predicted"]], d$X %*% w2)
})

test_that("twas_predict: name transformation weights -> predicted", {
  set.seed(42)
  wl <- list(
    lasso_weights = matrix(1, nrow = 3, ncol = 1),
    enet_weights  = matrix(1, nrow = 3, ncol = 1),
    susie_weights = matrix(1, nrow = 3, ncol = 1)
  )
  X <- matrix(rnorm(9), nrow = 3, ncol = 3)
  res <- twas_predict(X, wl)

  expect_equal(names(res), c("lasso_predicted", "enet_predicted", "susie_predicted"))
})

test_that("twas_predict: names without _weights suffix are kept unchanged", {
  wl <- list(custom_method = matrix(1, nrow = 2, ncol = 1))
  X <- matrix(1:4, nrow = 2, ncol = 2)
  res <- twas_predict(X, wl)

  # gsub("_weights", "_predicted", "custom_method") == "custom_method"
  expect_equal(names(res), "custom_method")
})

test_that("twas_predict: single column Y dimension preserved", {
  d <- make_data(n = 10, p = 4)
  w <- matrix(rep(0.25, 4), ncol = 1)
  wl <- list(avg_weights = w)
  res <- twas_predict(d$X, wl)

  expect_true(is.matrix(res[["avg_predicted"]]))
  expect_equal(nrow(res[["avg_predicted"]]), 10)
  expect_equal(ncol(res[["avg_predicted"]]), 1)
})

test_that("twas_predict: multi-column weights produce multi-column predictions", {
  set.seed(42)
  n <- 15
  p <- 5
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  W <- matrix(rnorm(p * 3), nrow = p, ncol = 3)
  wl <- list(multi_weights = W)
  res <- twas_predict(X, wl)

  expect_equal(ncol(res[["multi_predicted"]]), 3)
  expect_equal(res[["multi_predicted"]], X %*% W)
})

test_that("twas_predict: zero weights give zero predictions", {
  X <- matrix(1:6, nrow = 2, ncol = 3)
  wl <- list(null_weights = matrix(0, nrow = 3, ncol = 1))
  res <- twas_predict(X, wl)
  expect_true(all(res[["null_predicted"]] == 0))
})

# ===========================================================================
#
#  twas_weights  (input validation and basic behavior)
#
# ===========================================================================

test_that("twas_weights: X must be a matrix", {
  d <- make_data()
  expect_error(
    twas_weights(as.data.frame(d$X), d$Y, weight_methods = list()),
    "X must be a matrix"
  )
})

test_that("twas_weights: Y must be a matrix or vector", {
  d <- make_data()
  # In R, is.vector(list(...)) returns TRUE, so a list passes the initial
  # type check and gets converted via matrix(). The resulting matrix has
  # 1 row which mismatches X's 50 rows, triggering the row count error.
  expect_error(
    twas_weights(d$X, list(d$Y), weight_methods = list()),
    "The number of rows in X and Y must be the same"
  )
})

test_that("twas_weights: Y as vector gets converted to matrix internally", {
  d <- make_data()
  y_vec <- as.numeric(d$Y)

  # Mock lasso_weights (an existing package function) to return trivial weights
  local_mocked_bindings(
    lasso_weights = function(X, y, ...) rep(0, ncol(X))
  )
  result <- twas_weights(d$X, y_vec, weight_methods = list(lasso_weights = list()))
  expect_true(is.list(result))
  expect_equal(length(result), 1)
  expect_equal(nrow(result[[1]]), ncol(d$X))
  # Weight vector length must equal number of predictors and be numeric/finite
  w <- result[[1]][, 1]
  expect_equal(length(w), ncol(d$X))
  expect_true(is.numeric(w))
  expect_true(all(is.finite(w)))
})

test_that("twas_weights: mismatched row counts error", {
  d <- make_data(n = 50, p = 10)
  Y_short <- d$Y[1:30, , drop = FALSE]
  expect_error(
    twas_weights(d$X, Y_short, weight_methods = list()),
    "The number of rows in X and Y must be the same"
  )
})

test_that("twas_weights: character weight_methods input is accepted", {
  d <- make_data()
  local_mocked_bindings(
    lasso_weights = function(X, y, ...) rep(0, ncol(X))
  )
  # Character vector should be converted internally to named list

  result <- twas_weights(d$X, d$Y, weight_methods = c("lasso_weights"))
  expect_true(is.list(result))
  expect_equal(names(result), "lasso_weights")
})

test_that("twas_weights: zero variance columns are filtered and padded back with zeros", {
  d <- make_data(n = 50, p = 10, add_zero_var_col = TRUE)
  p_with_extra <- ncol(d$X)  # 11 columns, last is zero-var

  local_mocked_bindings(
    lasso_weights = function(X, y, ...) {
      # After filtering, the zero-var column should be removed
      # So ncol(X) should be p (10), not p+1 (11)
      rep(1, ncol(X))
    }
  )
  result <- twas_weights(d$X, d$Y, weight_methods = list(lasso_weights = list()))

  # The returned weight matrix should have rows equal to total columns (including zero-var)
  expect_equal(nrow(result[["lasso_weights"]]), p_with_extra)
  # The zero-var column weight should be 0 (padded back)
  expect_equal(result[["lasso_weights"]]["zero_var", 1], 0)
})

test_that("twas_weights: rownames of result match colnames of X", {
  d <- make_data()
  local_mocked_bindings(
    enet_weights = function(X, y, ...) rep(0.1, ncol(X))
  )
  result <- twas_weights(d$X, d$Y, weight_methods = list(enet_weights = list()))
  expect_equal(rownames(result[["enet_weights"]]), colnames(d$X))
})

test_that("twas_weights: result dimensions match ncol(X) x ncol(Y)", {
  d <- make_data()
  local_mocked_bindings(
    enet_weights = function(X, y, ...) rep(0, ncol(X))
  )
  result <- twas_weights(d$X, d$Y, weight_methods = list(enet_weights = list()))
  expect_equal(dim(result[["enet_weights"]]), c(ncol(d$X), ncol(d$Y)))
})

test_that("twas_weights: multiple methods return named list with one entry per method", {
  d <- make_data()
  local_mocked_bindings(
    lasso_weights = function(X, y, ...) rep(0.1, ncol(X)),
    enet_weights  = function(X, y, ...) rep(0.2, ncol(X))
  )
  result <- twas_weights(
    d$X, d$Y,
    weight_methods = list(lasso_weights = list(), enet_weights = list())
  )
  expect_equal(length(result), 2)
  expect_true("lasso_weights" %in% names(result))
  expect_true("enet_weights" %in% names(result))
})

# ===========================================================================
#
#  twas_weights with actual glmnet (lasso/enet) -- skip if not available
#
# ===========================================================================

test_that("twas_weights: lasso_weights produces correct structure with real glmnet", {
  skip_if_not_installed("glmnet")
  d <- make_data(n = 50, p = 10)
  result <- twas_weights(d$X, d$Y, weight_methods = list(lasso_weights = list()))

  expect_true(is.list(result))
  expect_equal(names(result), "lasso_weights")
  expect_equal(nrow(result[["lasso_weights"]]), ncol(d$X))
  expect_equal(ncol(result[["lasso_weights"]]), 1)
  # At least some weights should be non-zero for this strong signal
  expect_true(any(result[["lasso_weights"]] != 0))
})

test_that("twas_weights: enet_weights produces correct structure with real glmnet", {
  skip_if_not_installed("glmnet")
  d <- make_data(n = 50, p = 10)
  result <- twas_weights(d$X, d$Y, weight_methods = list(enet_weights = list()))

  expect_true(is.list(result))
  expect_equal(names(result), "enet_weights")
  expect_equal(nrow(result[["enet_weights"]]), ncol(d$X))
})

# ===========================================================================
#
#  twas_weights_cv  (input validation)
#
# ===========================================================================

test_that("twas_weights_cv: fold must be positive integer", {
  d <- make_data()
  expect_error(
    twas_weights_cv(d$X, d$Y, fold = -1),
    "Invalid value for 'fold'"
  )
  expect_error(
    twas_weights_cv(d$X, d$Y, fold = 0),
    "Invalid value for 'fold'"
  )
})

test_that("twas_weights_cv: fold as string errors", {
  d <- make_data()
  expect_error(
    twas_weights_cv(d$X, d$Y, fold = "abc"),
    "Invalid value for 'fold'"
  )
})

test_that("twas_weights_cv: X must be a matrix", {
  d <- make_data()
  expect_error(
    twas_weights_cv(as.data.frame(d$X), d$Y, fold = 5),
    "X must be a matrix"
  )
})

test_that("twas_weights_cv: Y must be a matrix or vector", {
  d <- make_data()
  # In R, is.vector(list(...)) returns TRUE, so a list passes the initial
  # type check and gets converted via matrix(). The resulting 3-row matrix
  # mismatches X's 50 rows, triggering the row count error.
  expect_error(
    twas_weights_cv(d$X, list(1, 2, 3), fold = 5),
    "The number of rows in X and Y must be the same"
  )
})

test_that("twas_weights_cv: fold or sample_partitions must be provided", {
  d <- make_data()
  expect_error(
    twas_weights_cv(d$X, d$Y),
    "Either 'fold' or 'sample_partitions' must be provided"
  )
})

test_that("twas_weights_cv: row count mismatch between X and Y errors", {
  d <- make_data()
  Y_wrong <- d$Y[1:20, , drop = FALSE]
  expect_error(
    twas_weights_cv(d$X, Y_wrong, fold = 5),
    "The number of rows in X and Y must be the same"
  )
})

test_that("twas_weights_cv: Y as vector is accepted and converted", {
  d <- make_data()
  y_vec <- as.numeric(d$Y)

  # With NULL weight_methods, should return just sample_partition
  expect_message(
    result <- twas_weights_cv(d$X, y_vec, fold = 3, weight_methods = NULL),
    "Y converted to matrix"
  )
  expect_true(is.list(result))
  expect_true("sample_partition" %in% names(result))
})

test_that("twas_weights_cv: NULL weight_methods returns only sample_partition", {
  d <- make_data()
  result <- twas_weights_cv(d$X, d$Y, fold = 3, weight_methods = NULL)
  expect_equal(names(result), "sample_partition")
  expect_true(is.data.frame(result$sample_partition))
})

test_that("twas_weights_cv: sample_partition structure is correct", {
  d <- make_data()
  result <- twas_weights_cv(d$X, d$Y, fold = 5, weight_methods = NULL)
  sp <- result$sample_partition

  expect_true("Sample" %in% colnames(sp))
  expect_true("Fold" %in% colnames(sp))
  expect_equal(nrow(sp), nrow(d$X))
  expect_equal(length(unique(sp$Fold)), 5)
  # All sample names should appear
  expect_true(all(rownames(d$X) %in% sp$Sample))
})

test_that("twas_weights_cv: character weight_methods are accepted", {
  d <- make_data()
  local_mocked_bindings(
    lasso_weights = function(X, y, ...) rep(0, ncol(X))
  )
  set.seed(42)
  result <- twas_weights_cv(
    d$X, d$Y, fold = 2,
    weight_methods = c("lasso_weights")
  )
  expect_true(is.list(result))
  expect_true("prediction" %in% names(result))
})

test_that("twas_weights_cv: max_num_variants subsets X columns", {
  d <- make_data(n = 50, p = 20)
  local_mocked_bindings(
    lasso_weights = function(X, y, ...) rep(0, ncol(X))
  )
  set.seed(42)
  expect_message(
    result <- twas_weights_cv(
      d$X, d$Y, fold = 2,
      weight_methods = list(lasso_weights = list()),
      max_num_variants = 5
    ),
    "Randomly selecting 5 out of 20"
  )
  expect_true(is.list(result))
  # The result should have prediction and performance entries for the one method
  expect_true("prediction" %in% names(result))
  expect_true("performance" %in% names(result))
  expect_true("lasso_predicted" %in% names(result$prediction))
  # Weight method returned zero weights, so predictions should exist but all be zero
  pred <- result$prediction[["lasso_predicted"]]
  expect_equal(nrow(pred), nrow(d$X))
  expect_true(all(pred == 0))
})

test_that("twas_weights_cv: max_num_variants with variants_to_keep", {
  d <- make_data(n = 50, p = 20)
  keep_vars <- colnames(d$X)[1:3]
  local_mocked_bindings(
    lasso_weights = function(X, y, ...) rep(0, ncol(X))
  )
  set.seed(42)
  expect_message(
    result <- twas_weights_cv(
      d$X, d$Y, fold = 2,
      weight_methods = list(lasso_weights = list()),
      max_num_variants = 8,
      variants_to_keep = keep_vars
    ),
    "Including 3 specified variants"
  )
  expect_true(is.list(result))
})

test_that("twas_weights_cv: sample_partitions with mismatched samples errors", {
  d <- make_data()
  bad_partitions <- data.frame(
    Sample = c("nonexistent_1", "nonexistent_2"),
    Fold = c(1, 2),
    stringsAsFactors = FALSE
  )
  expect_error(
    twas_weights_cv(d$X, d$Y, sample_partitions = bad_partitions),
    "Some samples in 'sample_partitions' do not match"
  )
})

test_that("twas_weights_cv: provided sample_partitions are used", {
  d <- make_data(n = 20, p = 5)
  sp <- data.frame(
    Sample = rownames(d$X),
    Fold = rep(1:2, each = 10),
    stringsAsFactors = FALSE
  )
  result <- twas_weights_cv(d$X, d$Y, sample_partitions = sp, weight_methods = NULL)
  expect_equal(result$sample_partition, sp)
})

test_that("twas_weights_cv: rownames are auto-generated when missing", {
  set.seed(42)
  n <- 30
  p <- 5
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- matrix(rnorm(n), ncol = 1)
  # No row names on X or Y

  result <- twas_weights_cv(X, Y, fold = 2, weight_methods = NULL)
  sp <- result$sample_partition

  # Should have auto-generated sample names
  expect_true(all(grepl("^sample_", sp$Sample)))
})

test_that("twas_weights_cv: colnames are auto-generated when missing", {
  set.seed(42)
  n <- 30
  p <- 5
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  rownames(X) <- paste0("s", 1:n)
  Y <- matrix(rnorm(n), ncol = 1)
  rownames(Y) <- paste0("s", 1:n)
  # No col names on X or Y

  local_mocked_bindings(
    lasso_weights = function(X, y, ...) rep(0, ncol(X))
  )
  set.seed(42)
  result <- twas_weights_cv(X, Y, fold = 2, weight_methods = list(lasso_weights = list()))
  # Should not error; column names are auto-generated
  expect_true(!is.null(result$prediction))
})

test_that("twas_weights_cv: fold and sample_partitions mismatch prints message", {
  d <- make_data(n = 20, p = 5)
  sp <- data.frame(
    Sample = rownames(d$X),
    Fold = rep(1:4, each = 5),
    stringsAsFactors = FALSE
  )
  # fold=2 but sample_partitions has 4 folds
  expect_message(
    twas_weights_cv(d$X, d$Y, fold = 2, sample_partitions = sp, weight_methods = NULL),
    "fold number provided does not match"
  )
})

test_that("twas_weights_cv: zero-variance predictions yield NA metrics with message", {
  d <- make_data()
  local_mocked_bindings(
    lasso_weights = function(X, y, ...) rep(0, ncol(X))
  )
  set.seed(42)
  expect_message(
    result <- twas_weights_cv(
      d$X, d$Y, fold = 2,
      weight_methods = list(lasso_weights = list())
    ),
    "zero variance"
  )
  perf <- result$performance[["lasso_performance"]]
  expect_true(all(is.na(perf)))
})

test_that("twas_weights_cv: performance names use _performance suffix", {
  d <- make_data()
  local_mocked_bindings(
    lasso_weights = function(X, y, ...) rep(0, ncol(X)),
    enet_weights  = function(X, y, ...) rep(0, ncol(X))
  )
  set.seed(42)
  result <- twas_weights_cv(
    d$X, d$Y, fold = 2,
    weight_methods = list(lasso_weights = list(), enet_weights = list())
  )
  expect_equal(
    sort(names(result$performance)),
    sort(c("lasso_performance", "enet_performance"))
  )
})

test_that("twas_weights_cv: prediction names use _predicted suffix", {
  d <- make_data()
  local_mocked_bindings(
    lasso_weights = function(X, y, ...) rep(0, ncol(X)),
    enet_weights  = function(X, y, ...) rep(0, ncol(X))
  )
  set.seed(42)
  result <- twas_weights_cv(
    d$X, d$Y, fold = 2,
    weight_methods = list(lasso_weights = list(), enet_weights = list())
  )
  expect_equal(
    sort(names(result$prediction)),
    sort(c("lasso_predicted", "enet_predicted"))
  )
})

# ---------------------------------------------------------------------------
# CV with real lasso_weights (integration test)
# ---------------------------------------------------------------------------

test_that("twas_weights_cv: basic CV with lasso_weights produces correct metrics structure", {
  skip_if_not_installed("glmnet")
  d <- make_data(n = 50, p = 10)

  set.seed(123)
  result <- twas_weights_cv(
    d$X, d$Y,
    fold = 3,
    weight_methods = list(lasso_weights = list())
  )

  # Structure checks
  expect_true(is.list(result))
  expect_true("sample_partition" %in% names(result))
  expect_true("prediction" %in% names(result))
  expect_true("performance" %in% names(result))
  expect_true("time_elapsed" %in% names(result))

  # Prediction name transformation
  expect_equal(names(result$prediction), "lasso_predicted")

  # Prediction dimensions should match Y
  pred <- result$prediction[["lasso_predicted"]]
  expect_equal(dim(pred), dim(d$Y))

  # Performance table structure
  perf <- result$performance[["lasso_performance"]]
  expect_true(is.matrix(perf))
  expect_equal(colnames(perf), c("corr", "rsq", "adj_rsq", "pval", "RMSE", "MAE"))
  expect_equal(nrow(perf), ncol(d$Y))

  # With strong signal, correlation should be positive
  expect_true(perf[1, "corr"] > 0)
})

test_that("twas_weights_cv: metrics table has correct column names (mocked)", {
  d <- make_data()
  local_mocked_bindings(
    lasso_weights = function(X, y, ...) {
      # Return weights that produce non-zero-variance predictions
      rep(0.1, ncol(X))
    }
  )
  set.seed(42)
  result <- twas_weights_cv(
    d$X, d$Y, fold = 2,
    weight_methods = list(lasso_weights = list())
  )
  perf <- result$performance[["lasso_performance"]]
  expect_equal(colnames(perf), c("corr", "rsq", "adj_rsq", "pval", "RMSE", "MAE"))
})

test_that("twas_weights_cv: multiple real methods produce per-method metrics", {
  skip_if_not_installed("glmnet")
  d <- make_data(n = 50, p = 10)

  set.seed(99)
  result <- twas_weights_cv(
    d$X, d$Y,
    fold = 3,
    weight_methods = list(
      lasso_weights = list(),
      enet_weights  = list()
    )
  )

  expect_equal(length(result$prediction), 2)
  expect_equal(length(result$performance), 2)
  expect_true("lasso_predicted" %in% names(result$prediction))
  expect_true("enet_predicted" %in% names(result$prediction))
  expect_true("lasso_performance" %in% names(result$performance))
  expect_true("enet_performance" %in% names(result$performance))
})

test_that("twas_weights_cv: all samples appear exactly once in predictions", {
  skip_if_not_installed("glmnet")
  d <- make_data(n = 50, p = 10)

  set.seed(77)
  result <- twas_weights_cv(
    d$X, d$Y,
    fold = 5,
    weight_methods = list(lasso_weights = list())
  )

  pred <- result$prediction[["lasso_predicted"]]
  # No NAs -- every sample was predicted in exactly one fold
  expect_false(any(is.na(pred)))
  expect_equal(nrow(pred), nrow(d$X))
})

# ===========================================================================
#
#  twas_weights_cv: multivariate Y
#
# ===========================================================================

test_that("twas_weights_cv: multivariate Y with multiple columns", {
  d <- make_data(n = 50, p = 10)
  # Create multi-column Y
  set.seed(42)
  Y_multi <- cbind(d$Y, d$X %*% c(0, 0, 0, 0, 0, 1, -1, 0, 0, 0) + rnorm(50, sd = 0.5))
  colnames(Y_multi) <- c("outcome_1", "outcome_2")
  rownames(Y_multi) <- rownames(d$X)

  local_mocked_bindings(
    lasso_weights = function(X, y, ...) rep(0.1, ncol(X))
  )
  set.seed(42)
  result <- twas_weights_cv(
    d$X, Y_multi, fold = 2,
    weight_methods = list(lasso_weights = list())
  )

  pred <- result$prediction[["lasso_predicted"]]
  expect_equal(ncol(pred), 2)
  expect_equal(nrow(pred), 50)

  perf <- result$performance[["lasso_performance"]]
  expect_equal(nrow(perf), 2)
  expect_equal(rownames(perf), c("outcome_1", "outcome_2"))
})

# ===========================================================================
#
#  twas_weights_pipeline  (structure and input validation)
#
# ===========================================================================

test_that("twas_weights_pipeline: returns list with expected structure (mocked)", {
  d <- make_data(n = 50, p = 10)
  y_vec <- as.numeric(d$Y)

  local_mocked_bindings(
    enet_weights  = function(X, y, ...) rep(0.1, ncol(X)),
    lasso_weights = function(X, y, ...) rep(0.2, ncol(X)),
    bayes_r_weights = function(X, y, ...) rep(0, ncol(X)),
    bayes_l_weights = function(X, y, ...) rep(0, ncol(X)),
    mrash_weights = function(X, y, ...) rep(0, ncol(X)),
    susie_weights = function(X, y, ...) rep(0, ncol(X))
  )

  result <- twas_weights_pipeline(d$X, y_vec, susie_fit = NULL, cv_folds = 0)

  expect_true(is.list(result))
  expect_true("twas_weights" %in% names(result))
  expect_true("twas_predictions" %in% names(result))
  expect_true("total_time_elapsed" %in% names(result))
  # Verify that mock values appear in the weight matrices
  enet_w <- result$twas_weights[["enet_weights"]]
  expect_true(all(enet_w[, 1] == 0.1))
  lasso_w <- result$twas_weights[["lasso_weights"]]
  expect_true(all(lasso_w[, 1] == 0.2))
  # The number of weight methods should equal the 6 default methods
  expect_equal(length(result$twas_weights), 6)
})

test_that("twas_weights_pipeline: twas_weights contains all default methods", {
  d <- make_data(n = 50, p = 10)
  y_vec <- as.numeric(d$Y)

  local_mocked_bindings(
    enet_weights  = function(X, y, ...) rep(0.1, ncol(X)),
    lasso_weights = function(X, y, ...) rep(0.2, ncol(X)),
    bayes_r_weights = function(X, y, ...) rep(0, ncol(X)),
    bayes_l_weights = function(X, y, ...) rep(0, ncol(X)),
    mrash_weights = function(X, y, ...) rep(0, ncol(X)),
    susie_weights = function(X, y, ...) rep(0, ncol(X))
  )

  result <- twas_weights_pipeline(d$X, y_vec, susie_fit = NULL, cv_folds = 0)

  expected_methods <- c(
    "enet_weights", "lasso_weights", "bayes_r_weights",
    "bayes_l_weights", "mrash_weights", "susie_weights"
  )
  expect_true(all(expected_methods %in% names(result$twas_weights)))
})

test_that("twas_weights_pipeline: predictions have _predicted suffix", {
  d <- make_data(n = 50, p = 10)
  y_vec <- as.numeric(d$Y)

  local_mocked_bindings(
    enet_weights  = function(X, y, ...) rep(0, ncol(X)),
    lasso_weights = function(X, y, ...) rep(0, ncol(X)),
    bayes_r_weights = function(X, y, ...) rep(0, ncol(X)),
    bayes_l_weights = function(X, y, ...) rep(0, ncol(X)),
    mrash_weights = function(X, y, ...) rep(0, ncol(X)),
    susie_weights = function(X, y, ...) rep(0, ncol(X))
  )

  result <- twas_weights_pipeline(d$X, y_vec, susie_fit = NULL, cv_folds = 0)

  expected_pred_names <- c(
    "enet_predicted", "lasso_predicted", "bayes_r_predicted",
    "bayes_l_predicted", "mrash_predicted", "susie_predicted"
  )
  expect_true(all(expected_pred_names %in% names(result$twas_predictions)))
})

test_that("twas_weights_pipeline: cv_folds=0 skips cross-validation", {
  d <- make_data(n = 50, p = 10)
  y_vec <- as.numeric(d$Y)

  local_mocked_bindings(
    enet_weights  = function(X, y, ...) rep(0, ncol(X)),
    lasso_weights = function(X, y, ...) rep(0, ncol(X)),
    bayes_r_weights = function(X, y, ...) rep(0, ncol(X)),
    bayes_l_weights = function(X, y, ...) rep(0, ncol(X)),
    mrash_weights = function(X, y, ...) rep(0, ncol(X)),
    susie_weights = function(X, y, ...) rep(0, ncol(X))
  )

  result <- twas_weights_pipeline(d$X, y_vec, susie_fit = NULL, cv_folds = 0)

  expect_false("twas_cv_result" %in% names(result))
  # All mock weights were zero, so all predictions should be zero
  for (pred_name in names(result$twas_predictions)) {
    expect_true(all(result$twas_predictions[[pred_name]] == 0),
                info = paste("Non-zero prediction in", pred_name))
  }
  # Weight dimensions should match ncol(X)
  for (w_name in names(result$twas_weights)) {
    expect_equal(nrow(result$twas_weights[[w_name]]), ncol(d$X),
                 info = paste("Wrong nrow for", w_name))
  }
})

test_that("twas_weights_pipeline: custom weight_methods are respected", {
  d <- make_data(n = 50, p = 10)
  y_vec <- as.numeric(d$Y)

  local_mocked_bindings(
    lasso_weights = function(X, y, ...) rep(1, ncol(X)),
    enet_weights  = function(X, y, ...) rep(2, ncol(X))
  )

  result <- twas_weights_pipeline(
    d$X, y_vec, susie_fit = NULL, cv_folds = 0,
    weight_methods = list(lasso_weights = list(), enet_weights = list())
  )

  expect_equal(sort(names(result$twas_weights)), sort(c("lasso_weights", "enet_weights")))
})

test_that("twas_weights_pipeline: with susie_fit stores intermediate results", {
  d <- make_data(n = 50, p = 10)
  y_vec <- as.numeric(d$Y)

  # Build a minimal fake susie_fit object
  fake_susie <- list(
    mu = matrix(0, nrow = 5, ncol = 10),
    lbf_variable = matrix(0, nrow = 5, ncol = 10),
    X_column_scale_factors = rep(1, 10),
    pip = rep(0.1, 10),
    V = rep(0.5, 5),
    sets = list(cs = NULL, purity = NULL)
  )

  local_mocked_bindings(
    enet_weights  = function(X, y, ...) rep(0, ncol(X)),
    lasso_weights = function(X, y, ...) rep(0, ncol(X)),
    bayes_r_weights = function(X, y, ...) rep(0, ncol(X)),
    bayes_l_weights = function(X, y, ...) rep(0, ncol(X)),
    mrash_weights = function(X, y, ...) rep(0, ncol(X)),
    susie_weights = function(X, y, ...) rep(0, ncol(X))
  )

  result <- twas_weights_pipeline(d$X, y_vec, susie_fit = fake_susie, cv_folds = 0)

  expect_true("susie_weights_intermediate" %in% names(result))
  expect_true("mu" %in% names(result$susie_weights_intermediate))
  expect_true("pip" %in% names(result$susie_weights_intermediate))
  expect_true("lbf_variable" %in% names(result$susie_weights_intermediate))
  expect_true("X_column_scale_factors" %in% names(result$susie_weights_intermediate))
})

test_that("twas_weights_pipeline: with susie_fit, susie_weights gets susie_fit argument", {
  d <- make_data(n = 50, p = 10)
  y_vec <- as.numeric(d$Y)

  fake_susie <- list(
    mu = matrix(0, nrow = 5, ncol = 10),
    lbf_variable = matrix(0, nrow = 5, ncol = 10),
    X_column_scale_factors = rep(1, 10),
    pip = rep(0.1, 10),
    V = rep(0.5, 5),
    sets = list(cs = NULL, purity = NULL)
  )

  susie_received_fit <- FALSE
  local_mocked_bindings(
    enet_weights  = function(X, y, ...) rep(0, ncol(X)),
    lasso_weights = function(X, y, ...) rep(0, ncol(X)),
    bayes_r_weights = function(X, y, ...) rep(0, ncol(X)),
    bayes_l_weights = function(X, y, ...) rep(0, ncol(X)),
    mrash_weights = function(X, y, ...) rep(0, ncol(X)),
    susie_weights = function(X, y, ...) {
      args <- list(...)
      # The pipeline should pass susie_fit as an argument
      if ("susie_fit" %in% names(args)) {
        susie_received_fit <<- TRUE
      }
      rep(0, ncol(X))
    }
  )

  result <- twas_weights_pipeline(d$X, y_vec, susie_fit = fake_susie, cv_folds = 0)
  expect_true(susie_received_fit)
})

test_that("twas_weights_pipeline: weight dimensions match input", {
  d <- make_data(n = 50, p = 10)
  y_vec <- as.numeric(d$Y)

  local_mocked_bindings(
    lasso_weights = function(X, y, ...) rep(0.5, ncol(X)),
    enet_weights  = function(X, y, ...) rep(0.3, ncol(X))
  )

  result <- twas_weights_pipeline(
    d$X, y_vec, susie_fit = NULL, cv_folds = 0,
    weight_methods = list(lasso_weights = list(), enet_weights = list())
  )

  for (method_name in names(result$twas_weights)) {
    w <- result$twas_weights[[method_name]]
    expect_equal(nrow(w), ncol(d$X))
    expect_equal(ncol(w), 1)
  }
})

# ===========================================================================
# twas_weights_cv: extra split_data / sample-name / variant-selection branches
# ===========================================================================

test_that("twas_weights_cv: split_data errors when a fold leaves train or test empty", {
  d <- make_data(n = 10, p = 5)
  # All samples in fold 1 -> with fold = 1, every sample is a "test" row and
  # the train set has zero rows, hitting the split_data zero-row stop.
  sp <- data.frame(
    Sample = rownames(d$X),
    Fold = rep(1L, nrow(d$X)),
    stringsAsFactors = FALSE
  )
  local_mocked_bindings(
    lasso_weights = function(X, y, ...) rep(0, ncol(X))
  )
  expect_error(
    suppressMessages(twas_weights_cv(
      d$X, d$Y,
      sample_partitions = sp,
      weight_methods = list(lasso_weights = list())
    )),
    "One of the datasets \\(train or test\\) has zero rows"
  )
})

test_that("twas_weights_cv: rownames(X) get reassigned to rownames(Y) when they differ", {
  d <- make_data(n = 20, p = 5)
  rownames(d$X) <- paste0("xname_", seq_len(nrow(d$X)))  # differ from rownames(Y)
  set.seed(42)
  result <- twas_weights_cv(d$X, d$Y, fold = 2, weight_methods = NULL)
  # sample_partition$Sample should now use rownames(Y), not rownames(X)
  expect_true(all(result$sample_partition$Sample %in% rownames(d$Y)))
  expect_false(any(grepl("^xname_", result$sample_partition$Sample)))
})

test_that("twas_weights_cv: sample_names taken from Y when only Y has rownames", {
  set.seed(42)
  n <- 20; p <- 5
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  rownames(X) <- NULL
  Y <- matrix(rnorm(n), ncol = 1)
  rownames(Y) <- paste0("yonly_", seq_len(n))
  result <- twas_weights_cv(X, Y, fold = 2, weight_methods = NULL)
  expect_true(all(grepl("^yonly_", result$sample_partition$Sample)))
})

test_that("twas_weights_cv: variants_to_keep >= max_num_variants samples from variants_to_keep only", {
  d <- make_data(n = 50, p = 20)
  # 10 keep variants, max_num_variants = 5 => length(variants_to_keep) >= max_num_variants
  keep_vars <- colnames(d$X)[1:10]
  local_mocked_bindings(
    lasso_weights = function(X, y, ...) rep(0, ncol(X))
  )
  set.seed(42)
  expect_message(
    result <- twas_weights_cv(
      d$X, d$Y, fold = 2,
      weight_methods = list(lasso_weights = list()),
      max_num_variants = 5,
      variants_to_keep = keep_vars
    ),
    "Randomly selecting 5 out of 10 input variants"
  )
  expect_true("prediction" %in% names(result))
})

test_that("twas_weights_cv: NA values in Y trigger NA-removal branch in metrics", {
  set.seed(42)
  n <- 30; p <- 5
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("v", seq_len(p))
  rownames(X) <- paste0("s", seq_len(n))
  Y <- matrix(rnorm(n), ncol = 1)
  rownames(Y) <- rownames(X); colnames(Y) <- "outcome"
  Y[c(3, 11, 17), 1] <- NA  # introduce NAs

  # Mock to return non-zero (so prediction has nonzero variance and lm_fit runs)
  local_mocked_bindings(
    lasso_weights = function(X, y, ...) {
      w <- rep(0, ncol(X)); w[1] <- 0.5; w
    }
  )
  set.seed(42)
  result <- twas_weights_cv(
    X, Y, fold = 2,
    weight_methods = list(lasso_weights = list())
  )
  perf <- result$performance[["lasso_performance"]]
  # NA-removal branch ran; metrics should be finite (not all-NA)
  expect_true(is.finite(perf[1, "rsq"]))
})

test_that("twas_weights_cv: multivariate cv_args data_driven_prior_matrices_cv is plumbed through", {
  set.seed(42)
  n <- 20; p <- 4
  X <- matrix(rnorm(n * p), nrow = n)
  colnames(X) <- paste0("v", seq_len(p)); rownames(X) <- paste0("s", seq_len(n))
  Y <- matrix(rnorm(n * 2), nrow = n)
  colnames(Y) <- c("y1", "y2"); rownames(Y) <- rownames(X)

  captured_args <- list()
  local_mocked_bindings(
    mrmash_weights = function(X, Y, ...) {
      captured_args[[length(captured_args) + 1]] <<- list(...)
      matrix(0, nrow = ncol(X), ncol = ncol(Y),
             dimnames = list(colnames(X), colnames(Y)))
    }
  )
  prior_cv <- list(matrix(1, 2, 2), matrix(2, 2, 2))
  set.seed(42)
  result <- twas_weights_cv(
    X, Y, fold = 2,
    weight_methods = list(mrmash_weights = list()),
    data_driven_prior_matrices_cv = prior_cv
  )
  # mrmash_weights mock should have been called and received the per-fold prior matrix
  expect_true(length(captured_args) >= 1)
  expect_true(any(vapply(captured_args, function(a)
    "data_driven_prior_matrices" %in% names(a), logical(1))))
})

# ===========================================================================
# twas_weights_pipeline: removed_methods warning + max_cv_variants subsampling
# ===========================================================================

test_that("twas_weights_pipeline: warns when methods are removed because all weights are zero", {
  d <- make_data(n = 30, p = 6)
  local_mocked_bindings(
    lasso_weights = function(X, y, ...) rep(0, ncol(X)),
    enet_weights  = function(X, y, ...) rep(0, ncol(X))
  )
  set.seed(42)
  expect_warning(
    suppressMessages(twas_weights_pipeline(
      d$X, d$Y, susie_fit = NULL, cv_folds = 2,
      weight_methods = list(lasso_weights = list(), enet_weights = list())
    )),
    "are removed from CV because all their weights are zeros"
  )
})

test_that("twas_weights_pipeline: max_cv_variants subsamples colnames of X", {
  d <- make_data(n = 30, p = 20)
  captured_keep <- NULL
  local_mocked_bindings(
    lasso_weights = function(X, y, ...) {
      w <- rep(0, ncol(X)); w[1] <- 0.5; w
    },
    twas_weights_cv = function(X, y, fold, sample_partition, weight_methods,
                               max_num_variants, num_threads, variants_to_keep, ...) {
      captured_keep <<- variants_to_keep
      list(sample_partition = data.frame(Sample = rownames(X), Fold = 1),
           prediction = list(), performance = list(), time_elapsed = 0)
    }
  )
  set.seed(42)
  suppressMessages(suppressWarnings(twas_weights_pipeline(
    d$X, d$Y, susie_fit = NULL, cv_folds = 2,
    weight_methods = list(lasso_weights = list()),
    max_cv_variants = 5
  )))
  expect_equal(length(captured_keep), 5)
  expect_true(all(captured_keep %in% colnames(d$X)))
})

# ===========================================================================
# twas_weights: dim-fix branch when nrow(weights_matrix) != length(valid_columns)
# ===========================================================================

test_that("twas_weights: multivariate weights_matrix is reduced to valid_columns when row counts mismatch", {
  set.seed(42)
  n <- 20; p <- 5
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("v", seq_len(p))  # all columns valid (no zero variance)
  Y <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  colnames(Y) <- c("y1", "y2")

  local_mocked_bindings(
    mrmash_weights = function(X, Y, ...) {
      # Return more rows than valid_columns (length p) so the dim-fix branch
      # subsets the matrix back to names(valid_columns).
      extra_rows <- p + 2
      m <- matrix(seq_len(extra_rows * ncol(Y)), nrow = extra_rows, ncol = ncol(Y))
      rownames(m) <- c(paste0("v", seq_len(p)), "extra1", "extra2")
      colnames(m) <- colnames(Y)
      m
    }
  )
  result <- twas_weights(X, Y, weight_methods = list(mrmash_weights = list()))
  # After the dim-fix, the weights matrix is restricted to v1..v5 -> shape p x ncol(Y)
  expect_equal(nrow(result$mrmash_weights), p)
  expect_equal(ncol(result$mrmash_weights), 2)
  expect_equal(rownames(result$mrmash_weights), paste0("v", seq_len(p)))
})
