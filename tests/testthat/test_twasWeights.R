context("twasWeights")

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

make_fake_susie_fit <- function(p = 10, L = 3, inf = FALSE) {
  fit <- list(
    alpha = matrix(1 / p, nrow = L, ncol = p),
    mu = matrix(0, nrow = L, ncol = p),
    lbf_variable = matrix(0, nrow = L, ncol = p),
    X_column_scale_factors = rep(1, p),
    pip = rep(0.1, p),
    V = rep(0.5, L),
    sets = list(cs = NULL, purity = NULL)
  )
  if (inf) fit$theta <- rep(0, p)
  fit
}

mock_susie <- function(...) {
  args <- list(...)
  L <- if (is.null(args$L)) 3 else args$L
  make_fake_susie_fit(ncol(args$X), L = L, inf = identical(args$unmappable_effects, "inf"))
}

# Test helper: fetch the weights matrix for a given method token from a
# TwasWeights collection. Accepts either the short token ("lasso") or the
# legacy suffixed name ("lassoWeights" / "lasso_weights"). Single-outcome
# entries store a bare numeric vector internally (drop()'d from the
# learnTwasWeights matrix); promote those back to a 1-column matrix here
# so test assertions on nrow/ncol/rownames keep working.
.weightsByMethod <- function(tw, method) {
  shortName <- sub("_?[Ww]eights$", "", method)
  idx <- which(as.character(tw$method) == shortName)
  if (length(idx) == 0L) return(NULL)
  entry <- tw$entry[[idx[[1L]]]]
  w <- entry@weights
  if (is.numeric(w) && is.null(dim(w))) {
    nm <- names(w)
    if (is.null(nm) && length(entry@variantIds) == length(w))
      nm <- entry@variantIds
    w <- matrix(w, ncol = 1L, dimnames = list(nm, NULL))
  }
  w
}

# ===========================================================================
#
#  .twas_method_lookup
#
# ===========================================================================

test_that(".twas_method_lookup: 'default' preset returns 10 methods", {
  result <- pecotmr:::.twasMethodLookup("default")
  expected_names <- c(
    "susie_weights", "susie_inf_weights", "mrash_weights", "enet_weights",
    "lasso_weights", "mcp_weights", "scad_weights", "l0learn_weights",
    "bayes_r_weights", "bayes_c_weights"
  )
  expect_equal(sort(names(result)), sort(expected_names))
})

test_that(".twas_method_lookup: 'fast_default' preset returns 8 methods", {
  result <- pecotmr:::.twasMethodLookup("fast_default")
  expected_names <- c(
    "susie_weights", "susie_inf_weights", "mrash_weights", "enet_weights",
    "lasso_weights", "mcp_weights", "scad_weights", "l0learn_weights"
  )
  expect_equal(sort(names(result)), sort(expected_names))
})

test_that(".twas_method_lookup: custom vector of short names", {
  result <- pecotmr:::.twasMethodLookup(c("susie", "enet", "dpr_vb"))
  expect_equal(sort(names(result)), sort(c("susie_weights", "enet_weights", "dpr_vb_weights")))
})

test_that(".twas_method_lookup: unknown method produces error", {
  expect_error(
    pecotmr:::.twasMethodLookup(c("susie", "nonexistent_method")),
    "Unknown TWAS method"
  )
})

test_that(".twas_method_lookup: default args are set for susie and mrash", {
  result <- pecotmr:::.twasMethodLookup("fast_default")
  expect_equal(result$susie_weights$refine, FALSE)
  expect_equal(result$susie_weights$L, 20)
  expect_equal(result$susie_weights$L_greedy, 5)
  expect_equal(result$mrash_weights$initPriorSd, TRUE)
  expect_equal(result$mrash_weights$max.iter, 100)
})

test_that(".twas_method_lookup: methods with no special args get empty list", {
  result <- pecotmr:::.twasMethodLookup(c("enet", "lasso"))
  expect_equal(length(result$enet_weights), 0L)
  expect_equal(length(result$lasso_weights), 0L)
})

test_that(".twas_method_lookup: all DPR variants can coexist", {
  result <- pecotmr:::.twasMethodLookup(c("dpr_vb", "dpr_gibbs", "dpr_adaptive_gibbs"))
  expect_equal(
    sort(names(result)),
    sort(c("dpr_vb_weights", "dpr_gibbs_weights", "dpr_adaptive_gibbs_weights"))
  )
})

# ===========================================================================
#
#  twasPredict
#
# ===========================================================================

test_that("twasPredict: basic matrix multiplication is correct", {
  d <- make_data(n = 20, p = 5)
  set.seed(99)
  w <- matrix(runif(5), ncol = 1)
  rownames(w) <- colnames(d$X)
  wl <- list(test_weights = w)
  res <- twasPredict(d$X, wl)

  expected <- d$X %*% w
  expect_equal(res[["test_predicted"]], expected)
})

test_that("twasPredict: multiple weight methods in list", {
  d <- make_data(n = 20, p = 5)
  w1 <- matrix(c(1, 0, 0, 0, 0), ncol = 1)
  w2 <- matrix(c(0, 0, 0, 0, 1), ncol = 1)
  wl <- list(method_a_weights = w1, method_b_weights = w2)

  res <- twasPredict(d$X, wl)

  expect_length(res, 2)
  expect_equal(res[["method_a_predicted"]], d$X %*% w1)
  expect_equal(res[["method_b_predicted"]], d$X %*% w2)
})

test_that("twasPredict: name transformation weights -> predicted", {
  set.seed(42)
  wl <- list(
    lassoWeights = matrix(1, nrow = 3, ncol = 1),
    enetWeights  = matrix(1, nrow = 3, ncol = 1),
    susieWeights = matrix(1, nrow = 3, ncol = 1)
  )
  X <- matrix(rnorm(9), nrow = 3, ncol = 3)
  res <- twasPredict(X, wl)

  expect_equal(names(res), c("lassoPredicted", "enetPredicted", "susiePredicted"))
})

test_that("twasPredict: names without _weights suffix are kept unchanged", {
  wl <- list(custom_method = matrix(1, nrow = 2, ncol = 1))
  X <- matrix(1:4, nrow = 2, ncol = 2)
  res <- twasPredict(X, wl)

  # gsub("_weights", "_predicted", "custom_method") == "custom_method"
  expect_equal(names(res), "custom_method")
})

test_that("twasPredict: single column Y dimension preserved", {
  d <- make_data(n = 10, p = 4)
  w <- matrix(rep(0.25, 4), ncol = 1)
  wl <- list(avg_weights = w)
  res <- twasPredict(d$X, wl)

  expect_true(is.matrix(res[["avg_predicted"]]))
  expect_equal(nrow(res[["avg_predicted"]]), 10)
  expect_equal(ncol(res[["avg_predicted"]]), 1)
})

test_that("twasPredict: multi-column weights produce multi-column predictions", {
  set.seed(42)
  n <- 15
  p <- 5
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  W <- matrix(rnorm(p * 3), nrow = p, ncol = 3)
  wl <- list(multi_weights = W)
  res <- twasPredict(X, wl)

  expect_equal(ncol(res[["multi_predicted"]]), 3)
  expect_equal(res[["multi_predicted"]], X %*% W)
})

test_that("twasPredict: zero weights give zero predictions", {
  X <- matrix(1:6, nrow = 2, ncol = 3)
  wl <- list(null_weights = matrix(0, nrow = 3, ncol = 1))
  res <- twasPredict(X, wl)
  expect_true(all(res[["null_predicted"]] == 0))
})

# ===========================================================================
#
#  twasWeights  (input validation and basic behavior)
#
# ===========================================================================

test_that("twasWeights: X must be a matrix", {
  d <- make_data()
  expect_error(
    learnTwasWeights(as.data.frame(d$X), d$Y, weightMethods = list()),
    "X must be a matrix"
  )
})

test_that("twasWeights: Y must be a matrix or vector", {
  d <- make_data()
  # In R, is.vector(list(...)) returns TRUE, so a list passes the initial
  # type check and gets converted via matrix(). The resulting matrix has
  # 1 row which mismatches X's 50 rows, triggering the row count error.
  expect_error(
    learnTwasWeights(d$X, list(d$Y), weightMethods = list()),
    "The number of rows in X and Y must be the same"
  )
})

test_that("twasWeights: Y as vector gets converted to matrix internally", {
  d <- make_data()
  y_vec <- as.numeric(d$Y)

  # Mock lassoWeights (an existing package function) to return trivial weights
  local_mocked_bindings(
    lassoWeights = function(X, y, ...) rep(0, ncol(X))
  )
  result <- learnTwasWeights(d$X, y_vec, weightMethods = list(lassoWeights = list()))
  expect_true(is(result, "TwasWeights"))
  expect_equal(length(getMethodNames(result)), 1)
  expect_equal(nrow(.weightsByMethod(result, "lassoWeights")), ncol(d$X))
  # Weight vector length must equal number of predictors and be numeric/finite
  w <- .weightsByMethod(result, "lassoWeights")[, 1]
  expect_equal(length(w), ncol(d$X))
  expect_true(is.numeric(w))
  expect_true(all(is.finite(w)))
})

test_that("twasWeights: mismatched row counts error", {
  d <- make_data(n = 50, p = 10)
  Y_short <- d$Y[1:30, , drop = FALSE]
  expect_error(
    learnTwasWeights(d$X, Y_short, weightMethods = list()),
    "The number of rows in X and Y must be the same"
  )
})

test_that("twasWeights: character weight_methods input is accepted", {
  d <- make_data()
  local_mocked_bindings(
    lassoWeights = function(X, y, ...) rep(0, ncol(X))
  )
  # Short name should be resolved via .twas_method_lookup
  result <- learnTwasWeights(d$X, d$Y, weightMethods = c("lasso"))
  expect_true(is(result, "TwasWeights"))
  expect_equal(getMethodNames(result), "lasso")
})

test_that("twasWeights: zero variance columns are filtered and padded back with zeros", {
  d <- make_data(n = 50, p = 10, add_zero_var_col = TRUE)
  p_with_extra <- ncol(d$X)  # 11 columns, last is zero-var

  local_mocked_bindings(
    lassoWeights = function(X, y, ...) {
      # After filtering, the zero-var column should be removed
      # So ncol(X) should be p (10), not p+1 (11)
      rep(1, ncol(X))
    }
  )
  result <- learnTwasWeights(d$X, d$Y, weightMethods = list(lassoWeights = list()))

  # The returned weight matrix should have rows equal to total columns (including zero-var)
  expect_equal(nrow(.weightsByMethod(result, "lassoWeights")), p_with_extra)
  # The zero-var column weight should be 0 (padded back)
  expect_equal(unname(.weightsByMethod(result, "lassoWeights")["zero_var", 1]), 0)
})

test_that("twasWeights: rownames of result match colnames of X", {
  d <- make_data()
  local_mocked_bindings(
    enetWeights = function(X, y, ...) rep(0.1, ncol(X))
  )
  result <- learnTwasWeights(d$X, d$Y, weightMethods = list(enetWeights = list()))
  expect_equal(rownames(.weightsByMethod(result, "enetWeights")), colnames(d$X))
})

test_that("twasWeights: result dimensions match ncol(X) x ncol(Y)", {
  d <- make_data()
  local_mocked_bindings(
    enetWeights = function(X, y, ...) rep(0, ncol(X))
  )
  result <- learnTwasWeights(d$X, d$Y, weightMethods = list(enetWeights = list()))
  expect_equal(dim(.weightsByMethod(result, "enetWeights")), c(ncol(d$X), ncol(d$Y)))
})

test_that("twasWeights: multiple methods return named list with one entry per method", {
  d <- make_data()
  local_mocked_bindings(
    lassoWeights = function(X, y, ...) rep(0.1, ncol(X)),
    enetWeights  = function(X, y, ...) rep(0.2, ncol(X))
  )
  result <- learnTwasWeights(
    d$X, d$Y,
    weightMethods = list(lassoWeights = list(), enetWeights = list())
  )
  expect_equal(length(getMethodNames(result)), 2)
  expect_true("lasso" %in% getMethodNames(result))
  expect_true("enet" %in% getMethodNames(result))
})

# ===========================================================================
#
#  twasWeights with actual glmnet (lasso/enet) -- skip if not available
#
# ===========================================================================

test_that("twasWeights: lassoWeights produces correct structure with real glmnet", {
  skip_if_not_installed("glmnet")
  d <- make_data(n = 50, p = 10)
  result <- learnTwasWeights(d$X, d$Y, weightMethods = list(lassoWeights = list()))

  expect_true(is(result, "TwasWeights"))
  expect_equal(getMethodNames(result), "lasso")
  expect_equal(nrow(.weightsByMethod(result, "lassoWeights")), ncol(d$X))
  expect_equal(ncol(.weightsByMethod(result, "lassoWeights")), 1)
  # At least some weights should be non-zero for this strong signal
  expect_true(any(.weightsByMethod(result, "lassoWeights") != 0))
})

test_that("twasWeights: enetWeights produces correct structure with real glmnet", {
  skip_if_not_installed("glmnet")
  d <- make_data(n = 50, p = 10)
  result <- learnTwasWeights(d$X, d$Y, weightMethods = list(enetWeights = list()))

  expect_true(is(result, "TwasWeights"))
  expect_equal(getMethodNames(result), "enet")
  expect_equal(nrow(.weightsByMethod(result, "enetWeights")), ncol(d$X))
})

# ===========================================================================
#
#  twasWeightsCv  (input validation)
#
# ===========================================================================

test_that("twasWeightsCv: fold must be positive integer", {
  d <- make_data()
  expect_error(
    twasWeightsCv(d$X, d$Y, fold = -1),
    "Invalid value for 'fold'"
  )
  expect_error(
    twasWeightsCv(d$X, d$Y, fold = 0),
    "Invalid value for 'fold'"
  )
})

test_that("twasWeightsCv: fold as string errors", {
  d <- make_data()
  expect_error(
    twasWeightsCv(d$X, d$Y, fold = "abc"),
    "Invalid value for 'fold'"
  )
})

test_that("twasWeightsCv: X must be a matrix", {
  d <- make_data()
  expect_error(
    twasWeightsCv(as.data.frame(d$X), d$Y, fold = 5),
    "X must be a matrix"
  )
})

test_that("twasWeightsCv: Y must be a matrix or vector", {
  d <- make_data()
  # In R, is.vector(list(...)) returns TRUE, so a list passes the initial
  # type check and gets converted via matrix(). The resulting 3-row matrix
  # mismatches X's 50 rows, triggering the row count error.
  expect_error(
    twasWeightsCv(d$X, list(1, 2, 3), fold = 5),
    "The number of rows in X and Y must be the same"
  )
})

test_that("twasWeightsCv: fold or sample_partitions must be provided", {
  d <- make_data()
  expect_error(
    twasWeightsCv(d$X, d$Y),
    "Either 'fold' or 'samplePartitions' must be provided"
  )
})

test_that("twasWeightsCv: row count mismatch between X and Y errors", {
  d <- make_data()
  Y_wrong <- d$Y[1:20, , drop = FALSE]
  expect_error(
    twasWeightsCv(d$X, Y_wrong, fold = 5),
    "The number of rows in X and Y must be the same"
  )
})

test_that("twasWeightsCv: Y as vector is accepted and converted", {
  d <- make_data()
  y_vec <- as.numeric(d$Y)

  # With NULL weight_methods, should return just samplePartition
  expect_message(
    result <- twasWeightsCv(d$X, y_vec, fold = 3, weightMethods = NULL),
    "Y converted to matrix"
  )
  expect_true(is.list(result))
  expect_true("samplePartition" %in% names(result))
})

test_that("twasWeightsCv: NULL weight_methods returns only samplePartition", {
  d <- make_data()
  result <- twasWeightsCv(d$X, d$Y, fold = 3, weightMethods = NULL)
  expect_equal(names(result), "samplePartition")
  expect_true(is.data.frame(result$samplePartition))
})

test_that("twasWeightsCv: samplePartition structure is correct", {
  d <- make_data()
  result <- twasWeightsCv(d$X, d$Y, fold = 5, weightMethods = NULL)
  sp <- result$samplePartition

  expect_true("Sample" %in% colnames(sp))
  expect_true("Fold" %in% colnames(sp))
  expect_equal(nrow(sp), nrow(d$X))
  expect_equal(length(unique(sp$Fold)), 5)
  # All sample names should appear
  expect_true(all(rownames(d$X) %in% sp$Sample))
})

test_that("twasWeightsCv: character weight_methods are accepted", {
  d <- make_data()
  local_mocked_bindings(
    lassoWeights = function(X, y, ...) rep(0, ncol(X))
  )
  set.seed(42)
  result <- twasWeightsCv(
    d$X, d$Y, fold = 2,
    weightMethods = c("lasso")
  )
  expect_true(is.list(result))
  expect_true("prediction" %in% names(result))
})

test_that("twasWeightsCv: max_num_variants subsets X columns", {
  d <- make_data(n = 50, p = 20)
  local_mocked_bindings(
    lassoWeights = function(X, y, ...) rep(0, ncol(X))
  )
  set.seed(42)
  expect_message(
    result <- twasWeightsCv(
      d$X, d$Y, fold = 2,
      weightMethods = list(lassoWeights = list()),
      maxNumVariants = 5
    ),
    "Randomly selecting 5 out of 20"
  )
  expect_true(is.list(result))
  # The result should have prediction and performance entries for the one method
  expect_true("prediction" %in% names(result))
  expect_true("performance" %in% names(result))
  expect_true("lassoPredicted" %in% names(result$prediction))
  # Weight method returned zero weights, so predictions should exist but all be zero
  pred <- result$prediction[["lassoPredicted"]]
  expect_equal(nrow(pred), nrow(d$X))
  expect_true(all(pred == 0))
})

test_that("twasWeightsCv: max_num_variants with variants_to_keep", {
  d <- make_data(n = 50, p = 20)
  keep_vars <- colnames(d$X)[1:3]
  local_mocked_bindings(
    lassoWeights = function(X, y, ...) rep(0, ncol(X))
  )
  set.seed(42)
  expect_message(
    result <- twasWeightsCv(
      d$X, d$Y, fold = 2,
      weightMethods = list(lassoWeights = list()),
      maxNumVariants = 8,
      variantsToKeep = keep_vars
    ),
    "Including 3 specified variants"
  )
  expect_true(is.list(result))
})

test_that("twasWeightsCv: sample_partitions with mismatched samples errors", {
  d <- make_data()
  bad_partitions <- data.frame(
    Sample = c("nonexistent_1", "nonexistent_2"),
    Fold = c(1, 2),
    stringsAsFactors = FALSE
  )
  expect_error(
    twasWeightsCv(d$X, d$Y, samplePartitions = bad_partitions),
    "Some samples in 'samplePartitions' do not match"
  )
})

test_that("twasWeightsCv: provided sample_partitions are used", {
  d <- make_data(n = 20, p = 5)
  sp <- data.frame(
    Sample = rownames(d$X),
    Fold = rep(1:2, each = 10),
    stringsAsFactors = FALSE
  )
  result <- twasWeightsCv(d$X, d$Y, samplePartitions = sp, weightMethods = NULL)
  expect_equal(result$samplePartition, sp)
})

test_that("twasWeightsCv: rownames are auto-generated when missing", {
  set.seed(42)
  n <- 30
  p <- 5
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- matrix(rnorm(n), ncol = 1)
  # No row names on X or Y

  result <- twasWeightsCv(X, Y, fold = 2, weightMethods = NULL)
  sp <- result$samplePartition

  # Should have auto-generated sample names
  expect_true(all(grepl("^sample_", sp$Sample)))
})

test_that("twasWeightsCv: colnames are auto-generated when missing", {
  set.seed(42)
  n <- 30
  p <- 5
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  rownames(X) <- paste0("s", 1:n)
  Y <- matrix(rnorm(n), ncol = 1)
  rownames(Y) <- paste0("s", 1:n)
  # No col names on X or Y

  local_mocked_bindings(
    lassoWeights = function(X, y, ...) rep(0, ncol(X))
  )
  set.seed(42)
  result <- twasWeightsCv(X, Y, fold = 2, weightMethods = list(lassoWeights = list()))
  # Should not error; column names are auto-generated
  expect_true(!is.null(result$prediction))
})

test_that("twasWeightsCv: fold and sample_partitions mismatch prints message", {
  d <- make_data(n = 20, p = 5)
  sp <- data.frame(
    Sample = rownames(d$X),
    Fold = rep(1:4, each = 5),
    stringsAsFactors = FALSE
  )
  # fold=2 but sample_partitions has 4 folds
  expect_message(
    twasWeightsCv(d$X, d$Y, fold = 2, samplePartitions = sp, weightMethods = NULL),
    "fold number provided does not match"
  )
})

test_that("twasWeightsCv: zero-variance predictions yield NA metrics with message", {
  d <- make_data()
  local_mocked_bindings(
    lassoWeights = function(X, y, ...) rep(0, ncol(X))
  )
  set.seed(42)
  expect_message(
    result <- twasWeightsCv(
      d$X, d$Y, fold = 2,
      weightMethods = list(lassoWeights = list())
    ),
    "zero variance"
  )
  perf <- result$performance[["lassoPerformance"]]
  expect_true(all(is.na(perf)))
})

test_that("twasWeightsCv: performance names use _performance suffix", {
  d <- make_data()
  local_mocked_bindings(
    lassoWeights = function(X, y, ...) rep(0, ncol(X)),
    enetWeights  = function(X, y, ...) rep(0, ncol(X))
  )
  set.seed(42)
  result <- twasWeightsCv(
    d$X, d$Y, fold = 2,
    weightMethods = list(lassoWeights = list(), enetWeights = list())
  )
  expect_equal(
    sort(names(result$performance)),
    sort(c("lassoPerformance", "enetPerformance"))
  )
})

test_that("twasWeightsCv: prediction names use _predicted suffix", {
  d <- make_data()
  local_mocked_bindings(
    lassoWeights = function(X, y, ...) rep(0, ncol(X)),
    enetWeights  = function(X, y, ...) rep(0, ncol(X))
  )
  set.seed(42)
  result <- twasWeightsCv(
    d$X, d$Y, fold = 2,
    weightMethods = list(lassoWeights = list(), enetWeights = list())
  )
  expect_equal(
    sort(names(result$prediction)),
    sort(c("lassoPredicted", "enetPredicted"))
  )
})

# ---------------------------------------------------------------------------
# CV with real lassoWeights (integration test)
# ---------------------------------------------------------------------------

test_that("twasWeightsCv: basic CV with lassoWeights produces correct metrics structure", {
  skip_if_not_installed("glmnet")
  d <- make_data(n = 50, p = 10)

  set.seed(123)
  result <- twasWeightsCv(
    d$X, d$Y,
    fold = 3,
    weightMethods = list(lassoWeights = list())
  )

  # Structure checks
  expect_true(is.list(result))
  expect_true("samplePartition" %in% names(result))
  expect_true("prediction" %in% names(result))
  expect_true("performance" %in% names(result))
  expect_true("timeElapsed" %in% names(result))

  # Prediction name transformation
  expect_equal(names(result$prediction), "lassoPredicted")

  # Prediction dimensions should match Y
  pred <- result$prediction[["lassoPredicted"]]
  expect_equal(dim(pred), dim(d$Y))

  # Performance table structure
  perf <- result$performance[["lassoPerformance"]]
  expect_true(is.matrix(perf))
  expect_equal(colnames(perf), c("corr", "rsq", "adj_rsq", "pval", "RMSE", "MAE"))
  expect_equal(nrow(perf), ncol(d$Y))

  # With strong signal, correlation should be positive
  expect_true(perf[1, "corr"] > 0)
})

test_that("twasWeightsCv: metrics table has correct column names (mocked)", {
  d <- make_data()
  local_mocked_bindings(
    lassoWeights = function(X, y, ...) {
      # Return weights that produce non-zero-variance predictions
      rep(0.1, ncol(X))
    }
  )
  set.seed(42)
  result <- twasWeightsCv(
    d$X, d$Y, fold = 2,
    weightMethods = list(lassoWeights = list())
  )
  perf <- result$performance[["lassoPerformance"]]
  expect_equal(colnames(perf), c("corr", "rsq", "adj_rsq", "pval", "RMSE", "MAE"))
})

test_that("twasWeightsCv: multiple real methods produce per-method metrics", {
  skip_if_not_installed("glmnet")
  d <- make_data(n = 50, p = 10)

  set.seed(99)
  result <- twasWeightsCv(
    d$X, d$Y,
    fold = 3,
    weightMethods = list(
      lassoWeights = list(),
      enetWeights  = list()
    )
  )

  expect_equal(length(result$prediction), 2)
  expect_equal(length(result$performance), 2)
  expect_true("lassoPredicted" %in% names(result$prediction))
  expect_true("enetPredicted" %in% names(result$prediction))
  expect_true("lassoPerformance" %in% names(result$performance))
  expect_true("enetPerformance" %in% names(result$performance))
})

test_that("twasWeightsCv: all samples appear exactly once in predictions", {
  skip_if_not_installed("glmnet")
  d <- make_data(n = 50, p = 10)

  set.seed(77)
  result <- twasWeightsCv(
    d$X, d$Y,
    fold = 5,
    weightMethods = list(lassoWeights = list())
  )

  pred <- result$prediction[["lassoPredicted"]]
  # No NAs -- every sample was predicted in exactly one fold
  expect_false(any(is.na(pred)))
  expect_equal(nrow(pred), nrow(d$X))
})

# ===========================================================================
#
#  twasWeightsCv: multivariate Y
#
# ===========================================================================

test_that("twasWeightsCv: multivariate Y with multiple columns", {
  d <- make_data(n = 50, p = 10)
  # Create multi-column Y
  set.seed(42)
  Y_multi <- cbind(d$Y, d$X %*% c(0, 0, 0, 0, 0, 1, -1, 0, 0, 0) + rnorm(50, sd = 0.5))
  colnames(Y_multi) <- c("outcome_1", "outcome_2")
  rownames(Y_multi) <- rownames(d$X)

  local_mocked_bindings(
    lassoWeights = function(X, y, ...) rep(0.1, ncol(X))
  )
  set.seed(42)
  result <- twasWeightsCv(
    d$X, Y_multi, fold = 2,
    weightMethods = list(lassoWeights = list())
  )

  pred <- result$prediction[["lassoPredicted"]]
  expect_equal(ncol(pred), 2)
  expect_equal(nrow(pred), 50)

  perf <- result$performance[["lassoPerformance"]]
  expect_equal(nrow(perf), 2)
  expect_equal(rownames(perf), c("outcome_1", "outcome_2"))
})

# ===========================================================================
#
#  twasWeightsPipeline  (structure and input validation)
#
# ===========================================================================

test_that("twasWeightsPipeline: returns list with expected structure (mocked)", {
  d <- make_data(n = 50, p = 10)
  y_vec <- as.numeric(d$Y)

  local_mocked_bindings(
    susie = mock_susie,
    enetWeights  = function(X, y, ...) rep(0.1, ncol(X)),
    lassoWeights = function(X, y, ...) rep(0.2, ncol(X)),
    bayesRWeights = function(X, y, ...) rep(0, ncol(X)),
    bayesCWeights = function(X, y, ...) rep(0, ncol(X)),
    mrashWeights = function(X, y, ...) rep(0, ncol(X)),
    mcpWeights   = function(X, y, ...) rep(0, ncol(X)),
    scadWeights  = function(X, y, ...) rep(0, ncol(X)),
    l0learnWeights =function(X, y, ...) rep(0, ncol(X)),
    susieWeights = function(X, y, ...) rep(0, ncol(X)),
    susieInfWeights = function(X, y, ...) rep(0, ncol(X))
  )

  result <- pecotmr:::.twasWeightsPipelineMatrix(d$X, y_vec, susieFit = NULL, cvFolds = 0,
                                  estimatePi = FALSE)

  expect_true(is.list(result))
  expect_true("twasWeights" %in% names(result))
  expect_true("twasPredictions" %in% names(result))
  expect_true("totalTimeElapsed" %in% names(result))
  # Verify that mock values appear in the weight matrices
  enet_w <- .weightsByMethod(result$twasWeights, "enet")
  expect_true(all(enet_w[, 1] == 0.1))
  lasso_w <- .weightsByMethod(result$twasWeights, "lasso")
  expect_true(all(lasso_w[, 1] == 0.2))
  # The number of weight methods should equal the 10 default methods
  expect_equal(length(getMethodNames(result$twasWeights)), 10)
})

test_that("twasWeightsPipeline: twasWeights contains all default methods", {
  d <- make_data(n = 50, p = 10)
  y_vec <- as.numeric(d$Y)

  local_mocked_bindings(
    susie = mock_susie,
    enetWeights  = function(X, y, ...) rep(0.1, ncol(X)),
    lassoWeights = function(X, y, ...) rep(0.2, ncol(X)),
    bayesRWeights = function(X, y, ...) rep(0, ncol(X)),
    bayesCWeights = function(X, y, ...) rep(0, ncol(X)),
    mrashWeights = function(X, y, ...) rep(0, ncol(X)),
    mcpWeights   = function(X, y, ...) rep(0, ncol(X)),
    scadWeights  = function(X, y, ...) rep(0, ncol(X)),
    l0learnWeights =function(X, y, ...) rep(0, ncol(X)),
    susieWeights = function(X, y, ...) rep(0, ncol(X)),
    susieInfWeights = function(X, y, ...) rep(0, ncol(X))
  )

  result <- pecotmr:::.twasWeightsPipelineMatrix(d$X, y_vec, susieFit = NULL, cvFolds = 0,
                                  estimatePi = FALSE)

  expected_methods <- c(
    "enet", "lasso", "bayes_r",
    "bayes_c", "mrash", "mcp",
    "scad", "l0learn", "susie",
    "susie_inf"
  )
  expect_true(all(expected_methods %in% getMethodNames(result$twasWeights)))
})

test_that("twasWeightsPipeline: stores ensemble weights when ensemble is fitted", {
  d <- make_data(n = 50, p = 10)
  y_vec <- as.numeric(d$Y)

  cv_perf <- matrix(NA_real_, nrow = 1, ncol = 6)
  colnames(cv_perf) <- c("corr", "rsq", "adj_rsq", "pval", "RMSE", "MAE")
  rownames(cv_perf) <- "outcome_1"
  cv_perf[1, "rsq"] <- 0.5

  local_mocked_bindings(
    enetWeights = function(X, y, ...) rep(0.1, ncol(X)),
    lassoWeights = function(X, y, ...) rep(0.2, ncol(X)),
    twasWeightsCv = function(X, Y, ...) {
      list(
        prediction = list(
          enetPredicted = matrix(as.numeric(Y), ncol = 1, dimnames = list(rownames(X), "outcome_1")),
          lassoPredicted = matrix(as.numeric(Y), ncol = 1, dimnames = list(rownames(X), "outcome_1"))
        ),
        performance = list(
          enetPerformance = cv_perf,
          lassoPerformance = cv_perf
        )
      )
    },
    ensembleWeights = function(cvResults, Y, twasWeightList, ...) {
      list(
        methodCoef = c(enet = 0.5, lasso = 0.5),
        ensembleTwasWeights = (twasWeightList$enet_weights + twasWeightList$lasso_weights) / 2
      )
    }
  )

  result <- pecotmr:::.twasWeightsPipelineMatrix(
    d$X, y_vec,
    weightMethods = list(enetWeights = list(), lassoWeights = list()),
    cvFolds = 2,
    ensemble = TRUE,
    ensembleR2Threshold = 0,
    estimatePi = FALSE
  )

  expect_true("ensemble" %in% getMethodNames(result$twasWeights))
  expect_true("ensemble_predicted" %in% names(result$twasPredictions))
  expect_true("ensemble" %in% names(result))
})

test_that("twasWeightsPipeline: predictions have _predicted suffix", {
  d <- make_data(n = 50, p = 10)
  y_vec <- as.numeric(d$Y)

  local_mocked_bindings(
    susie = mock_susie,
    enetWeights  = function(X, y, ...) rep(0, ncol(X)),
    lassoWeights = function(X, y, ...) rep(0, ncol(X)),
    bayesRWeights = function(X, y, ...) rep(0, ncol(X)),
    bayesCWeights = function(X, y, ...) rep(0, ncol(X)),
    mrashWeights = function(X, y, ...) rep(0, ncol(X)),
    mcpWeights   = function(X, y, ...) rep(0, ncol(X)),
    scadWeights  = function(X, y, ...) rep(0, ncol(X)),
    l0learnWeights =function(X, y, ...) rep(0, ncol(X)),
    susieWeights = function(X, y, ...) rep(0, ncol(X)),
    susieInfWeights = function(X, y, ...) rep(0, ncol(X))
  )

  result <- pecotmr:::.twasWeightsPipelineMatrix(d$X, y_vec, susieFit = NULL, cvFolds = 0,
                                  estimatePi = FALSE)

  expected_pred_names <- c(
    "enet_predicted", "lasso_predicted", "bayes_r_predicted",
    "bayes_c_predicted", "mrash_predicted", "mcp_predicted",
    "scad_predicted", "l0learn_predicted", "susie_predicted",
    "susie_inf_predicted"
  )
  expect_true(all(expected_pred_names %in% names(result$twasPredictions)))
})

test_that("twasWeightsPipeline: cv_folds=0 skips cross-validation", {
  d <- make_data(n = 50, p = 10)
  y_vec <- as.numeric(d$Y)

  local_mocked_bindings(
    susie = mock_susie,
    enetWeights  = function(X, y, ...) rep(0, ncol(X)),
    lassoWeights = function(X, y, ...) rep(0, ncol(X)),
    bayesRWeights = function(X, y, ...) rep(0, ncol(X)),
    bayesCWeights = function(X, y, ...) rep(0, ncol(X)),
    mrashWeights = function(X, y, ...) rep(0, ncol(X)),
    mcpWeights   = function(X, y, ...) rep(0, ncol(X)),
    scadWeights  = function(X, y, ...) rep(0, ncol(X)),
    l0learnWeights =function(X, y, ...) rep(0, ncol(X)),
    susieWeights = function(X, y, ...) rep(0, ncol(X)),
    susieInfWeights = function(X, y, ...) rep(0, ncol(X))
  )

  result <- pecotmr:::.twasWeightsPipelineMatrix(d$X, y_vec, susieFit = NULL, cvFolds = 0,
                                  estimatePi = FALSE)

  expect_false("twasCvResult" %in% names(result))
  # All mock weights were zero, so all predictions should be zero
  for (pred_name in names(result$twasPredictions)) {
    expect_true(all(result$twasPredictions[[pred_name]] == 0),
                info = paste("Non-zero prediction in", pred_name))
  }
  # Weight dimensions should match ncol(X)
  for (w_name in getMethodNames(result$twasWeights)) {
    expect_equal(nrow(.weightsByMethod(result$twasWeights, w_name)), ncol(d$X),
                 info = paste("Wrong nrow for", w_name))
  }
})

test_that("twasWeightsPipeline: custom weight_methods are respected", {
  d <- make_data(n = 50, p = 10)
  y_vec <- as.numeric(d$Y)

  local_mocked_bindings(
    lassoWeights = function(X, y, ...) rep(1, ncol(X)),
    enetWeights  = function(X, y, ...) rep(2, ncol(X))
  )

  result <- pecotmr:::.twasWeightsPipelineMatrix(
    d$X, y_vec, susieFit = NULL, cvFolds = 0,
    weightMethods = list(lassoWeights = list(), enetWeights = list())
  )

  expect_equal(sort(getMethodNames(result$twasWeights)), sort(c("lasso", "enet")))
})

test_that("twasWeightsPipeline: accepts 'fast_default' preset string", {
  d <- make_data(n = 50, p = 10)
  y_vec <- as.numeric(d$Y)

  local_mocked_bindings(
    susie = mock_susie,
    enetWeights  = function(X, y, ...) rep(0, ncol(X)),
    lassoWeights = function(X, y, ...) rep(0, ncol(X)),
    mrashWeights = function(X, y, ...) rep(0, ncol(X)),
    mcpWeights   = function(X, y, ...) rep(0, ncol(X)),
    scadWeights  = function(X, y, ...) rep(0, ncol(X)),
    l0learnWeights =function(X, y, ...) rep(0, ncol(X)),
    susieWeights = function(X, y, ...) rep(0, ncol(X)),
    susieInfWeights = function(X, y, ...) rep(0, ncol(X))
  )

  result <- pecotmr:::.twasWeightsPipelineMatrix(
    d$X, y_vec, susieFit = NULL, cvFolds = 0,
    weightMethods = "fast_default"
  )

  expected_methods <- c("susie", "susie_inf", "mrash",
                        "enet", "lasso", "mcp",
                        "scad", "l0learn")
  expect_equal(sort(getMethodNames(result$twasWeights)), sort(expected_methods))
})

test_that("twasWeightsPipeline: accepts custom short-name vector", {
  d <- make_data(n = 50, p = 10)
  y_vec <- as.numeric(d$Y)

  local_mocked_bindings(
    lassoWeights = function(X, y, ...) rep(1, ncol(X)),
    enetWeights  = function(X, y, ...) rep(2, ncol(X))
  )

  result <- pecotmr:::.twasWeightsPipelineMatrix(
    d$X, y_vec, susieFit = NULL, cvFolds = 0,
    weightMethods = c("lasso", "enet")
  )

  expect_equal(sort(getMethodNames(result$twasWeights)), sort(c("lasso", "enet")))
})

test_that("twasWeightsPipeline: with fitted_models stores SuSiE intermediates", {
  d <- make_data(n = 50, p = 10)
  y_vec <- as.numeric(d$Y)

  fake_susie <- make_fake_susie_fit(p = 10, L = 5)
  local_mocked_bindings(
    enetWeights  = function(X, y, ...) rep(0, ncol(X)),
    lassoWeights = function(X, y, ...) rep(0, ncol(X)),
    bayesRWeights = function(X, y, ...) rep(0, ncol(X)),
    bayesCWeights = function(X, y, ...) rep(0, ncol(X)),
    mrashWeights = function(X, y, ...) rep(0, ncol(X)),
    mcpWeights   = function(X, y, ...) rep(0, ncol(X)),
    scadWeights  = function(X, y, ...) rep(0, ncol(X)),
    l0learnWeights =function(X, y, ...) rep(0, ncol(X)),
    susieWeights = function(X, y, ...) rep(0, ncol(X)),
    susieInfWeights = function(X, y, ...) rep(0, ncol(X))
  )

  result <- pecotmr:::.twasWeightsPipelineMatrix(
    d$X, y_vec,
    fittedModels = list(susie = fake_susie),
    cvFolds = 0,
    estimatePi = FALSE
  )

  expect_true("susieWeightsIntermediate" %in% names(result))
  expect_true("mu" %in% names(result$susieWeightsIntermediate))
})

test_that("twasWeightsPipeline: fitted_models are injected into SuSiE-family weights", {
  d <- make_data(n = 50, p = 10)
  y_vec <- as.numeric(d$Y)

  fake_susie <- make_fake_susie_fit(p = 10, L = 5)
  fake_susie_inf <- make_fake_susie_fit(p = 10, L = 5, inf = TRUE)

  susie_received_fit <- FALSE
  susie_inf_received_fit <- FALSE
  local_mocked_bindings(
    enetWeights  = function(X, y, ...) rep(0, ncol(X)),
    lassoWeights = function(X, y, ...) rep(0, ncol(X)),
    bayesRWeights = function(X, y, ...) rep(0, ncol(X)),
    bayesCWeights = function(X, y, ...) rep(0, ncol(X)),
    mrashWeights = function(X, y, ...) rep(0, ncol(X)),
    mcpWeights   = function(X, y, ...) rep(0, ncol(X)),
    scadWeights  = function(X, y, ...) rep(0, ncol(X)),
    l0learnWeights =function(X, y, ...) rep(0, ncol(X)),
    susieInfWeights = function(X, y, ...) {
      args <- list(...)
      if (!is.null(args$susieInfFit) && "susieInf" %in% class(args$susieInfFit)) {
        susie_inf_received_fit <<- TRUE
      }
      rep(0, ncol(X))
    },
    susieWeights = function(X, y, ...) {
      args <- list(...)
      if (!is.null(args$susieFit) && "susie" %in% class(args$susieFit)) {
        susie_received_fit <<- TRUE
      }
      rep(0, ncol(X))
    }
  )

  result <- pecotmr:::.twasWeightsPipelineMatrix(
    d$X, y_vec,
    fittedModels = list(susie = fake_susie, susieInf = fake_susie_inf),
    cvFolds = 0,
    estimatePi = FALSE
  )
  expect_true(susie_received_fit)
  expect_true(susie_inf_received_fit)
})

test_that("twasWeights: SuSiE-inf is fitted before and initializes ordinary SuSiE", {
  d <- make_data(n = 50, p = 10)
  y_vec <- as.numeric(d$Y)
  susie_calls <- list()

  local_mocked_bindings(
    susie = function(...) {
      args <- list(...)
      susie_calls[[length(susie_calls) + 1]] <<- args
      make_fake_susie_fit(
        p = ncol(args$X),
        L = if (identical(args$unmappable_effects, "inf")) 7 else args$L,
        inf = identical(args$unmappable_effects, "inf")
      )
    },
    susieInfWeights = function(X, y, ...) rep(0, ncol(X)),
    susieWeights = function(X, y, ...) {
      rep(0, ncol(X))
    }
  )

  result <- learnTwasWeights(
    d$X,
    y_vec,
    weightMethods = list(
      susie_weights = list(L = 5, L_greedy = 3),
      susie_inf_weights = list()
    )
  )

  expect_equal(getMethodNames(result), c("susie", "susie_inf"))
  expect_length(susie_calls, 2)
  expect_equal(susie_calls[[1]]$unmappable_effects, "inf")
  expect_equal(susie_calls[[1]]$convergence_method, "pip")
  expect_equal(susie_calls[[2]]$unmappable_effects, "none")
  expect_true("susieInf" %in% class(susie_calls[[2]]$model_init))
  expect_equal(susie_calls[[2]]$L_greedy, 5)
})

test_that("twasWeightsPipeline: weight dimensions match input", {
  d <- make_data(n = 50, p = 10)
  y_vec <- as.numeric(d$Y)

  local_mocked_bindings(
    lassoWeights = function(X, y, ...) rep(0.5, ncol(X)),
    enetWeights  = function(X, y, ...) rep(0.3, ncol(X))
  )

  result <- pecotmr:::.twasWeightsPipelineMatrix(
    d$X, y_vec, susieFit = NULL, cvFolds = 0,
    weightMethods = list(lassoWeights = list(), enetWeights = list())
  )

  for (method_name in getMethodNames(result$twasWeights)) {
    w <- .weightsByMethod(result$twasWeights, method_name)
    expect_equal(nrow(w), ncol(d$X))
    expect_equal(ncol(w), 1)
  }
})

# ===========================================================================
# twasWeightsCv: extra split_data / sample-name / variant-selection branches
# ===========================================================================

test_that("twasWeightsCv: split_data errors when a fold leaves train or test empty", {
  d <- make_data(n = 10, p = 5)
  # All samples in fold 1 -> with fold = 1, every sample is a "test" row and
  # the train set has zero rows, hitting the split_data zero-row stop.
  sp <- data.frame(
    Sample = rownames(d$X),
    Fold = rep(1L, nrow(d$X)),
    stringsAsFactors = FALSE
  )
  local_mocked_bindings(
    lassoWeights = function(X, y, ...) rep(0, ncol(X))
  )
  expect_error(
    suppressMessages(twasWeightsCv(
      d$X, d$Y,
      samplePartitions = sp,
      weightMethods = list(lassoWeights = list())
    )),
    "One of the datasets \\(train or test\\) has zero rows"
  )
})

test_that("twasWeightsCv: rownames(X) get reassigned to rownames(Y) when they differ", {
  d <- make_data(n = 20, p = 5)
  rownames(d$X) <- paste0("xname_", seq_len(nrow(d$X)))  # differ from rownames(Y)
  set.seed(42)
  result <- twasWeightsCv(d$X, d$Y, fold = 2, weightMethods = NULL)
  # samplePartition$Sample should now use rownames(Y), not rownames(X)
  expect_true(all(result$samplePartition$Sample %in% rownames(d$Y)))
  expect_false(any(grepl("^xname_", result$samplePartition$Sample)))
})

test_that("twasWeightsCv: sampleNames taken from Y when only Y has rownames", {
  set.seed(42)
  n <- 20; p <- 5
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  rownames(X) <- NULL
  Y <- matrix(rnorm(n), ncol = 1)
  rownames(Y) <- paste0("yonly_", seq_len(n))
  result <- twasWeightsCv(X, Y, fold = 2, weightMethods = NULL)
  expect_true(all(grepl("^yonly_", result$samplePartition$Sample)))
})

test_that("twasWeightsCv: variants_to_keep >= max_num_variants samples from variants_to_keep only", {
  d <- make_data(n = 50, p = 20)
  # 10 keep variants, maxNumVariants = 5 => length(variants_to_keep) >= max_num_variants
  keep_vars <- colnames(d$X)[1:10]
  local_mocked_bindings(
    lassoWeights = function(X, y, ...) rep(0, ncol(X))
  )
  set.seed(42)
  expect_message(
    result <- twasWeightsCv(
      d$X, d$Y, fold = 2,
      weightMethods = list(lassoWeights = list()),
      maxNumVariants = 5,
      variantsToKeep = keep_vars
    ),
    "Randomly selecting 5 out of 10 input variants"
  )
  expect_true("prediction" %in% names(result))
})

test_that("twasWeightsCv: NA values in Y trigger NA-removal branch in metrics", {
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
    lassoWeights = function(X, y, ...) {
      w <- rep(0, ncol(X)); w[1] <- 0.5; w
    }
  )
  set.seed(42)
  result <- twasWeightsCv(
    X, Y, fold = 2,
    weightMethods = list(lassoWeights = list())
  )
  perf <- result$performance[["lassoPerformance"]]
  # NA-removal branch ran; metrics should be finite (not all-NA)
  expect_true(is.finite(perf[1, "rsq"]))
})

test_that("twasWeightsCv: multivariate cv_args data_driven_prior_matrices_cv is plumbed through", {
  set.seed(42)
  n <- 20; p <- 4
  X <- matrix(rnorm(n * p), nrow = n)
  colnames(X) <- paste0("v", seq_len(p)); rownames(X) <- paste0("s", seq_len(n))
  Y <- matrix(rnorm(n * 2), nrow = n)
  colnames(Y) <- c("y1", "y2"); rownames(Y) <- rownames(X)

  captured_args <- list()
  local_mocked_bindings(
    mrmashWeights = function(X, Y, ...) {
      captured_args[[length(captured_args) + 1]] <<- list(...)
      matrix(0, nrow = ncol(X), ncol = ncol(Y),
             dimnames = list(colnames(X), colnames(Y)))
    }
  )
  prior_cv <- list(matrix(1, 2, 2), matrix(2, 2, 2))
  set.seed(42)
  result <- twasWeightsCv(
    X, Y, fold = 2,
    weightMethods = list(mrmashWeights = list()),
    data_driven_prior_matrices_cv = prior_cv
  )
  # mrmashWeights mock should have been called and received the per-fold prior matrix
  expect_true(length(captured_args) >= 1)
  expect_true(any(vapply(captured_args, function(a)
    "data_driven_prior_matrices" %in% names(a), logical(1))))
})

# ===========================================================================
# twasWeightsPipeline: removed_methods warning + max_cv_variants subsampling
# ===========================================================================

test_that("twasWeightsPipeline: warns when methods are removed because all weights are zero", {
  d <- make_data(n = 30, p = 6)
  local_mocked_bindings(
    lassoWeights = function(X, y, ...) rep(0, ncol(X)),
    enetWeights  = function(X, y, ...) rep(0, ncol(X))
  )
  set.seed(42)
  expect_warning(
    suppressMessages(pecotmr:::.twasWeightsPipelineMatrix(
      d$X, d$Y, susieFit = NULL, cvFolds = 2,
      weightMethods = list(lassoWeights = list(), enetWeights = list())
    )),
    "are removed from CV because all their weights are zeros"
  )
})

test_that("twasWeightsPipeline: max_cv_variants subsamples colnames of X", {
  d <- make_data(n = 30, p = 20)
  captured_keep <- NULL
  local_mocked_bindings(
    lassoWeights = function(X, y, ...) {
      w <- rep(0, ncol(X)); w[1] <- 0.5; w
    },
    twasWeightsCv = function(X, Y, fold, samplePartitions, weightMethods,
                               maxNumVariants, numThreads, variantsToKeep, ...) {
      captured_keep <<- variantsToKeep
      list(samplePartition = data.frame(Sample = rownames(X), Fold = 1),
           prediction = list(), performance = list(), timeElapsed = 0)
    }
  )
  set.seed(42)
  suppressMessages(suppressWarnings(pecotmr:::.twasWeightsPipelineMatrix(
    d$X, d$Y, susieFit = NULL, cvFolds = 2,
    weightMethods = list(lassoWeights = list()),
    maxCvVariants = 5
  )))
  expect_equal(length(captured_keep), 5)
  expect_true(all(captured_keep %in% colnames(d$X)))
})

# ===========================================================================
# twasWeights: dim-fix branch when nrow(weights_matrix) != length(valid_columns)
# ===========================================================================

test_that("twasWeights: multivariate weights_matrix is reduced to valid_columns when row counts mismatch", {
  set.seed(42)
  n <- 20; p <- 5
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("v", seq_len(p))  # all columns valid (no zero variance)
  Y <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  colnames(Y) <- c("y1", "y2")

  local_mocked_bindings(
    mrmashWeights = function(X, Y, ...) {
      # Return more rows than valid_columns (length p) so the dim-fix branch
      # subsets the matrix back to names(valid_columns).
      extra_rows <- p + 2
      m <- matrix(seq_len(extra_rows * ncol(Y)), nrow = extra_rows, ncol = ncol(Y))
      rownames(m) <- c(paste0("v", seq_len(p)), "extra1", "extra2")
      colnames(m) <- colnames(Y)
      m
    }
  )
  result <- learnTwasWeights(X, Y, weightMethods = list(mrmashWeights = list()))
  # After the dim-fix, the weights matrix is restricted to v1..v5 -> shape p x ncol(Y)
  expect_equal(nrow(.weightsByMethod(result, "mrmashWeights")), p)
  expect_equal(ncol(.weightsByMethod(result, "mrmashWeights")), 2)
  expect_equal(rownames(.weightsByMethod(result, "mrmashWeights")), paste0("v", seq_len(p)))
})

# ===========================================================================
# Tests migrated from test_twas.R (twasWeightsCv, learnTwasWeights, twasPredict)
# ===========================================================================

test_that("twasWeightsCv is reproducible with seed", {
    sim <- generate_X_Y(seed=1)
    X <- sim$X
    y = sim$Y
    local_mocked_bindings(
        susieWeights = function(X, y, ...) rnorm(ncol(X)),
        glmnetWeights = function(X, y, ...) runif(ncol(X))
    )
    weight_methods_test <- list(susieWeights = list(), glmnetWeights = list())
    set.seed(1)
    result_seed1 <- twasWeightsCv(X, y, fold = 2, weightMethods = weight_methods_test)
    set.seed(1)
    result_seed2 <- twasWeightsCv(X, y, fold = 2, weightMethods = weight_methods_test)
    expect_equal(result_seed1$samplePartition, result_seed2$samplePartition)
})


test_that("twasWeightsCv handles errors appropriately", {
    sim <- generate_X_Y(seed=1)
    X <- sim$X
    y = sim$Y
    local_mocked_bindings(
        susieWeights = function(X, y, ...) rnorm(ncol(X)),
        glmnetWeights = function(X, y, ...) runif(ncol(X))
    )
    weight_methods_test <- list(susieWeights = list(), glmnetWeights = list())
    expect_error(twasWeightsCv(X, y, fold = NULL), "fold.*samplePartitions")
    expect_error(twasWeightsCv(X, y, fold = "invalid"), "positive integer")
    expect_error(twasWeightsCv(X, y, fold = -1), "positive integer")
    expect_error(twasWeightsCv(2, y, fold = 2), "must be a matrix")
    expect_error(twasWeightsCv(X, 2, fold = 2), "number of rows")
    expect_error(twasWeightsCv(matrix(rnorm(4, nrow=2)), matrix(rnorm(2, nrow=1)), fold = 2), "unused argument")
    expect_error(twasWeightsCv(X, y), "fold.*samplePartitions")
})


test_that("learnTwasWeights handles errors appropriately", {
    sim <- generate_X_Y(seed=1)
    X <- sim$X
    y = sim$Y
    local_mocked_bindings(
        susieWeights = function(X, y, ...) rnorm(ncol(X)),
        glmnetWeights = function(X, y, ...) runif(ncol(X))
    )
    weight_methods_test <- list(susieWeights = list(), glmnetWeights = list())
    expect_error(learnTwasWeights(matrix(rnorm(4, nrow=2)), matrix(rnorm(2, nrow=1))), "unused argument")
    expect_error(learnTwasWeights(X, y), "weightMethods")
})

# ===========================================================================
# twasZ: mathematical correctness (single-method / vector path)
# ===========================================================================


test_that("twasPredict multiplies X by weights", {
  X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
  weights_list <- list(method1_weights = c(0.5, -0.5))
  result <- twasPredict(X, weights_list)
  expect_length(result, 1)
  expect_equal(names(result), "method1_predicted")
  expected <- X %*% c(0.5, -0.5)
  expect_equal(result[[1]], expected)
})


test_that("twasPredict handles multiple weight methods", {
  set.seed(42)
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  weights_list <- list(
    lassoWeights = c(1, 0, -1),
    enetWeights = c(0.5, 0.3, 0.2),
    susieWeights = c(0, 0, 1)
  )
  result <- twasPredict(X, weights_list)
  expect_length(result, 3)
  expect_equal(names(result), c("lassoPredicted", "enetPredicted", "susiePredicted"))
  # Verify computation for one method
  expect_equal(result$lassoPredicted, X %*% c(1, 0, -1))
})


test_that("twasPredict with zero weights gives zero predictions", {
  X <- matrix(1:6, nrow = 2, ncol = 3)
  weights_list <- list(null_weights = c(0, 0, 0))
  result <- twasPredict(X, weights_list)
  expect_true(all(result$null_predicted == 0))
})


test_that("twasPredict with single variant", {
  X <- matrix(c(1, 2, 3), nrow = 3, ncol = 1)
  weights_list <- list(single_weights = 2.0)
  result <- twasPredict(X, weights_list)
  expect_equal(as.numeric(result$single_predicted), c(2, 4, 6))
})


# === Tests migrated from test_s4Constructors.R (TwasWeights) ===

test_that("TwasWeights: builds a collection keyed by 4-tuple", {
  e1 <- .sc_makeTwasWeightsEntry()
  e2 <- .sc_makeTwasWeightsEntry()
  tw <- TwasWeights(
    study   = c("s1", "s1"),
    context = c("c1", "c1"),
    trait   = c("t1", "t1"),
    method  = c("lasso", "enet"),
    entry   = list(e1, e2))
  expect_s4_class(tw, "TwasWeights")
  expect_equal(nrow(tw), 2L)
  expect_setequal(getMethodNames(tw), c("lasso", "enet"))
})


test_that("TwasWeights: getStudy / getContexts / getTraits / getMethodNames", {
  e <- .sc_makeTwasWeightsEntry()
  tw <- TwasWeights(
    study   = c("s1", "s2"),
    context = c("c1", "c2"),
    trait   = c("t1", "t1"),
    method  = c("lasso", "lasso"),
    entry   = list(e, e))
  expect_setequal(getContexts(tw), c("c1", "c2"))
  expect_equal(getTraits(tw), "t1")
  expect_equal(getMethodNames(tw), "lasso")
})


test_that("TwasWeights: rejects duplicate 4-tuples", {
  e <- .sc_makeTwasWeightsEntry()
  expect_error(
    TwasWeights(
      study   = c("s1", "s1"),
      context = c("c1", "c1"),
      trait   = c("t1", "t1"),
      method  = c("lasso", "lasso"),
      entry   = list(e, e)),
    "uniqueness violated"
  )
})


test_that("TwasWeights: joint columns work the same as on the FMR class", {
  e <- .sc_makeTwasWeightsEntry()
  tw <- TwasWeights(
    study   = c("s1", "s1"),
    context = c("c1", "joint"),
    trait   = c("t1", "t1"),
    method  = c("lasso", "mrmash"),
    entry   = list(e, e),
    jointContexts = c(NA_character_, "c1;c2"))
  expect_true("jointContexts" %in% names(tw))
  expect_identical(tw$jointContexts, c(NA_character_, "c1;c2"))
  # uniqueness: distinct jointContexts -> distinct rows
  tw2 <- TwasWeights(
    study   = c("s1", "s1"),
    context = c("joint", "joint"),
    trait   = c("t1", "t1"),
    method  = c("mrmash", "mrmash"),
    entry   = list(e, e),
    jointContexts = c("c1;c2", "c1;c3"))
  expect_equal(nrow(tw2), 2L)
})


test_that("TwasWeights: getTwasWeights extracts the entry for a tuple", {
  e1 <- .sc_makeTwasWeightsEntry()
  e2 <- .sc_makeTwasWeightsEntry()
  tw <- TwasWeights(
    study   = c("s1", "s1"),
    context = c("c1", "c1"),
    trait   = c("t1", "t1"),
    method  = c("lasso", "enet"),
    entry   = list(e1, e2))
  expect_identical(
    getTwasWeights(tw, study = "s1", context = "c1",
                   trait = "t1", method = "enet"),
    e2)
})

# ===========================================================================
# LdData
# ===========================================================================



# === Tests migrated from test_showMethods.R (TwasWeights) ===

test_that("show.TwasWeights prints entry/study/context/trait/method counts", {
  e <- .sh_makeTwEntry()
  tw <- TwasWeights(
    study   = c("s1", "s1"),
    context = c("c1", "c2"),
    trait   = c("t1", "t1"),
    method  = c("lasso", "enet"),
    entry   = list(e, e))
  out <- capture.output(show(tw))
  expect_true(any(grepl("TwasWeights: 2 entries", out)))
  expect_true(any(grepl("1 studies.*2 contexts.*1 traits.*2 methods", out)))
})


test_that("show.TwasWeights reports ldSketch when present", {
  e <- .sh_makeTwEntry()
  tw <- TwasWeights(
    study = "s1", context = "c1", trait = "t1", method = "lasso",
    entry = list(e),
    ldSketch = .sh_makeGenotypeHandle())
  out <- capture.output(show(tw))
  expect_true(any(grepl("LD sketch: gds @ /tmp/test.gds", out)))
})


