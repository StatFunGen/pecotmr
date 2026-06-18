context("twas")
library(testthat)

generate_X_Y <- function(seed=1, num_samples=10, num_features=10, X_rownames=TRUE, y_rownames=TRUE) {
  set.seed(seed)
  X <- scale(
    matrix(rnorm(num_samples * num_features), nrow = num_samples),
    center = TRUE, scale = TRUE)

  if (X_rownames) {
    rownames(X) <- paste0("sample", 1:num_samples)
  } else {
    rownames(X) <- NULL
  }

  beta = rep(0, num_features)
  beta[1:4] = 1
  y <- X %*% beta + rnorm(num_samples)
  y <- matrix(y, nrow = num_samples, ncol = 1)
  if (y_rownames) {
    rownames(y) <- paste0("sample", 1:num_samples)
  } else {
    rownames(y) <- NULL
  }
  colnames(y) <- c("Outcome")

  return(list(X=X, Y=y))
}

# ===========================================================================
# twasZ (post-S4-refactor): unified Z-statistic; output shape:
#   list(Z = K x 2 matrix with columns c("Z", "pval"), combined = NULL|list)
# ===========================================================================

test_that("twasZ errors when weights and z lengths differ", {
  expect_error(twasZ(c(1, 2), c(1, 2, 3)), "must equal")
})

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

test_that("twasZ: single weight and single z-score returns valid result", {
  weights <- 0.7
  z <- 2.5
  R <- matrix(1, nrow = 1, ncol = 1)
  result <- twasZ(weights, z, R = R)
  expect_true(is.list(result))
  expect_equal(colnames(result$Z), c("Z", "pval"))
  # With single variant: stat = 0.7 * 2.5 = 1.75, denom = 0.7 * 1 * 0.7 = 0.49
  # zscore = 1.75 / sqrt(0.49) = 1.75 / 0.7 = 2.5
  expect_equal(as.numeric(result$Z[, "Z"]), 2.5, tolerance = 1e-10)
  expect_true(result$Z[, "pval"] > 0 && result$Z[, "pval"] < 1)
})

test_that("twasZ: all-zero weights produce NaN z-score", {
  weights <- c(0, 0, 0)
  z <- c(1.5, -0.5, 2.0)
  R <- diag(3)
  result <- twasZ(weights, z, R = R)
  # stat = 0, denom = 0, so zscore = 0/0 = NaN
  expect_true(is.nan(as.numeric(result$Z[, "Z"])))
})

test_that("twasZ: very large z-scores still produce finite results", {
  set.seed(42)
  p <- 5
  weights <- rnorm(p)
  z <- rep(1e6, p)
  R <- diag(p)
  result <- twasZ(weights, z, R = R)
  expect_true(is.finite(as.numeric(result$Z[, "Z"])))
  # p-value should be extremely small for large z
  expect_true(result$Z[, "pval"] < 1e-10 || result$Z[, "pval"] == 0)
})

test_that("twasZ: identical z-scores with equal weights gives proportional result", {
  p <- 5
  weights <- rep(1, p)
  z <- rep(3.0, p)
  R <- diag(p)
  result <- twasZ(weights, z, R = R)
  # stat = sum(weights * z) = 5 * 3 = 15
  # denom = t(w) %*% I %*% w = 5
  # zscore = 15 / sqrt(5) = 6.7082...
  expect_equal(as.numeric(result$Z[, "Z"]), 15 / sqrt(5), tolerance = 1e-10)
})

test_that("twasZ: negative weights flip the sign of the z-score", {
  weights_pos <- c(0.5, 0.3)
  weights_neg <- c(-0.5, -0.3)
  z <- c(2.0, 1.0)
  R <- diag(2)
  result_pos <- twasZ(weights_pos, z, R = R)
  result_neg <- twasZ(weights_neg, z, R = R)
  # z-score should have opposite sign but same p-value
  expect_equal(as.numeric(result_pos$Z[, "Z"]),
               -as.numeric(result_neg$Z[, "Z"]),
               tolerance = 1e-10)
  expect_equal(as.numeric(result_pos$Z[, "pval"]),
               as.numeric(result_neg$Z[, "pval"]),
               tolerance = 1e-10)
})

test_that("twasZ: off-diagonal correlation in R changes the result", {
  weights <- c(0.5, 0.5)
  z <- c(2.0, 2.0)
  R_identity <- diag(2)
  R_correlated <- matrix(c(1, 0.8, 0.8, 1), nrow = 2)
  result_identity <- twasZ(weights, z, R = R_identity)
  result_correlated <- twasZ(weights, z, R = R_correlated)
  # Same stat but different denominators, so different z-scores
  expect_false(isTRUE(all.equal(
    as.numeric(result_identity$Z[, "Z"]),
    as.numeric(result_correlated$Z[, "Z"]))))
})

test_that("twasZ: computing R from X matches providing R directly", {
  set.seed(123)
  n <- 20
  p <- 5
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("SNP", 1:p)
  R <- cor(X)
  weights <- rnorm(p)
  z <- rnorm(p)
  result_with_R <- twasZ(weights, z, R = R)
  result_with_X <- twasZ(weights, z, X = X)
  expect_equal(as.numeric(result_with_R$Z[, "Z"]),
               as.numeric(result_with_X$Z[, "Z"]),
               tolerance = 1e-6)
  expect_equal(as.numeric(result_with_R$Z[, "pval"]),
               as.numeric(result_with_X$Z[, "pval"]),
               tolerance = 1e-6)
})

test_that("twasZ: p-value is always between 0 and 1 for random inputs", {
  set.seed(999)
  for (i in 1:5) {
    p <- sample(2:10, 1)
    weights <- rnorm(p)
    z <- rnorm(p)
    R <- diag(p) # use identity to avoid singularity
    result <- twasZ(weights, z, R = R)
    pval <- as.numeric(result$Z[, "pval"])
    expect_true(pval >= 0 && pval <= 1,
      info = paste("Iteration", i, "p-value out of range:", pval))
  }
})

test_that("twasZ: single very large weight with tiny z gives moderate result", {
  weights <- c(1e6, 0, 0)
  z <- c(1e-6, 5.0, 5.0)
  R <- diag(3)
  result <- twasZ(weights, z, R = R)
  # stat = 1e6 * 1e-6 + 0 + 0 = 1.0
  # denom = (1e6)^2 * 1 = 1e12
  # zscore = 1 / sqrt(1e12) = 1e-6
  expect_equal(as.numeric(result$Z[, "Z"]), 1e-6, tolerance = 1e-10)
})

# ===========================================================================
# twasZ: more mathematical edge cases
# ===========================================================================

test_that("twasZ: perfectly correlated R matrix (all ones off-diagonal)", {
  p <- 4
  R <- matrix(1, nrow = p, ncol = p) # perfectly correlated
  weights <- c(0.25, 0.25, 0.25, 0.25)
  z <- c(2, 2, 2, 2)
  result <- twasZ(weights, z, R = R)
  # stat = sum(0.25 * 2) = 2
  # denom = t(w) %*% ones_matrix %*% w = (sum(w))^2 = 1
  # zscore = 2 / sqrt(1) = 2
  expect_equal(as.numeric(result$Z[, "Z"]), 2.0, tolerance = 1e-10)
})

test_that("twasZ: sparse weights (only one non-zero) extracts single SNP signal", {
  p <- 5
  weights <- c(0, 0, 1, 0, 0)
  z <- c(1, 2, 3, 4, 5)
  R <- diag(p)
  result <- twasZ(weights, z, R = R)
  # With identity R and sparse weight: zscore = w3 * z3 / sqrt(w3^2) = 3
  expect_equal(as.numeric(result$Z[, "Z"]), 3.0, tolerance = 1e-10)
})

test_that("twasZ: z-scores of zero give zero TWAS z-score", {
  weights <- c(0.5, 0.3, 0.2)
  z <- c(0, 0, 0)
  R <- diag(3)
  result <- twasZ(weights, z, R = R)
  expect_equal(as.numeric(result$Z[, "Z"]), 0.0, tolerance = 1e-10)
  # p-value for z=0 should be 1
  expect_equal(as.numeric(result$Z[, "pval"]), 1.0, tolerance = 1e-10)
})

# ===========================================================================
# twasZ: edge case with near-singular R matrix
# ===========================================================================

test_that("twasZ: near-singular R matrix still produces a result", {
  set.seed(777)
  p <- 3
  # Create a nearly singular R by making two rows almost identical
  R <- matrix(c(1, 0.999, 0.999, 0.999, 1, 0.999, 0.999, 0.999, 1), nrow = 3)
  weights <- c(0.3, 0.4, 0.3)
  z <- c(2.0, 2.5, 1.8)
  result <- twasZ(weights, z, R = R)
  expect_true(is.list(result))
  expect_true(is.finite(as.numeric(result$Z[, "Z"])))
})

test_that("twasZ: length-one input vectors produce correct scalar output", {
  result <- twasZ(1.0, 3.0, R = matrix(1, 1, 1))
  expect_equal(as.numeric(result$Z[, "Z"]), 3.0, tolerance = 1e-10)
  expect_true(result$Z[, "pval"] > 0 && result$Z[, "pval"] < 1)
})

test_that("twasZ: R=NULL and X=NULL still errors consistently", {
  # When neither R nor X is provided, twasZ should error
  expect_error(twasZ(c(1, 2), c(3, 4)))
})

# ===========================================================================
# twasZ: matrix weights path -- multiple methods/conditions
# ===========================================================================

test_that("twasZ: matrix weights produce one Z row per column", {
  set.seed(10)
  p <- 5
  k <- 3
  weights <- matrix(rnorm(p * k), nrow = p, ncol = k)
  rownames(weights) <- paste0("SNP", 1:p)
  colnames(weights) <- paste0("Cond", 1:k)
  z <- rnorm(p)
  R <- diag(p)
  rownames(R) <- colnames(R) <- paste0("SNP", 1:p)
  result <- twasZ(weights, z, R = R)
  expect_equal(nrow(result$Z), k)
  expect_equal(rownames(result$Z), paste0("Cond", 1:k))
  expect_equal(colnames(result$Z), c("Z", "pval"))
  # combineMethods omitted -> combined is NULL
  expect_null(result$combined)
})

test_that("twasZ: combineMethods returns combined p-value summary", {
  skip_if_not_installed("ACAT")
  set.seed(11)
  p <- 4
  k <- 2
  weights <- matrix(rnorm(p * k), nrow = p, ncol = k,
                    dimnames = list(paste0("SNP", 1:p), paste0("Cond", 1:k)))
  z <- rnorm(p)
  R <- diag(p)
  rownames(R) <- colnames(R) <- paste0("SNP", 1:p)
  result <- twasZ(weights, z, R = R, combineMethods = "ACAT")
  expect_false(is.null(result$combined))
})

# ===========================================================================
# Tests from test_twas_predict.R
# ===========================================================================

# ---- twasPredict ----
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

# ---- cleanContextNames ----
test_that("cleanContextNames removes gene suffix", {
  context <- c("brain_ENSG00000123", "liver_ENSG00000123")
  gene <- "ENSG00000123"
  result <- cleanContextNames(context, gene)
  expect_equal(result, c("brain", "liver"))
})

test_that("cleanContextNames handles multiple genes", {
  context <- c("brain_GENE1", "liver_GENE2", "heart_GENE1")
  gene <- c("GENE1", "GENE2")
  result <- cleanContextNames(context, gene)
  expect_equal(result, c("brain", "liver", "heart"))
})

test_that("cleanContextNames no match leaves context unchanged", {
  context <- c("brain_ABC", "liver_XYZ")
  gene <- "NONEXISTENT"
  result <- cleanContextNames(context, gene)
  expect_equal(result, c("brain_ABC", "liver_XYZ"))
})

test_that("cleanContextNames longer gene removed first", {
  context <- c("tissue_ENSG00000123456")
  gene <- c("ENSG00000123", "ENSG00000123456")
  result <- cleanContextNames(context, gene)
  # Longer gene should be processed first
  expect_equal(result, "tissue")
})
