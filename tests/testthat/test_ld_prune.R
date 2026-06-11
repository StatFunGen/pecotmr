library(testthat)

# ===========================================================================
# ldPruneByCorrelation
# ===========================================================================

test_that("ldPruneByCorrelation removes highly correlated columns", {
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  colnames(X) <- paste0("v", 1:p)
  X[, 2] <- X[, 1] + rnorm(n, sd = 0.01)
  result <- ldPruneByCorrelation(X, corThres =0.9)
  expect_true(ncol(result$X.new) < p)
  expect_equal(length(result$filter.id), ncol(result$X.new))
})

test_that("ldPruneByCorrelation keeps all columns when uncorrelated", {
  set.seed(42)
  n <- 100; p <- 5
  X <- matrix(rnorm(n * p), nrow = n)
  colnames(X) <- paste0("v", 1:p)
  result <- ldPruneByCorrelation(X, corThres =0.99)
  expect_equal(ncol(result$X.new), p)
  expect_equal(result$filter.id, 1:p)
})

test_that("ldPruneByCorrelation preserves colnames for single remaining column", {
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n * 3), nrow = n)
  colnames(X) <- c("a", "b", "c")
  X[, 2] <- X[, 1] + rnorm(n, sd = 0.001)
  X[, 3] <- X[, 1] + rnorm(n, sd = 0.001)
  result <- ldPruneByCorrelation(X, corThres =0.5)
  expect_true(ncol(result$X.new) >= 1)
  expect_true(!is.null(colnames(result$X.new)))
})

test_that("ldPruneByCorrelation errors on single-column input", {
  set.seed(42)
  n <- 30
  X <- matrix(rnorm(n), nrow = n, ncol = 1)
  colnames(X) <- "v1"
  expect_error(ldPruneByCorrelation(X, corThres =0.8))
})

test_that("ldPruneByCorrelation strict threshold removes at least as many as lenient", {
  set.seed(42)
  n <- 100; p <- 5
  X <- matrix(rnorm(n * p), nrow = n)
  colnames(X) <- paste0("v", 1:p)
  X[, 2] <- X[, 1] + rnorm(n, sd = 0.1)
  X[, 3] <- X[, 1] + rnorm(n, sd = 0.1)
  X[, 5] <- X[, 4] + rnorm(n, sd = 0.1)
  result_strict <- ldPruneByCorrelation(X, corThres =0.3)
  result_lenient <- ldPruneByCorrelation(X, corThres =0.99)
  expect_true(ncol(result_strict$X.new) <= ncol(result_lenient$X.new))
})

test_that("ldPruneByCorrelation preserves colnames when no columns deleted", {
  set.seed(42)
  n <- 100; p <- 3
  X <- matrix(rnorm(n * p), nrow = n)
  colnames(X) <- c("snp_a", "snp_b", "snp_c")
  result <- ldPruneByCorrelation(X, corThres =0.999)
  expect_equal(colnames(result$X.new), colnames(X))
})

test_that("ldPruneByCorrelation is silent by default, chatty with verbose", {
  set.seed(1)
  X <- matrix(rnorm(100), 20, 5)
  colnames(X) <- paste0("v", 1:5)
  X[, 2] <- X[, 1] + rnorm(20, sd = 1e-3)
  expect_silent(ldPruneByCorrelation(X, corThres =0.9))
  expect_message(ldPruneByCorrelation(X, corThres =0.9, verbose = TRUE),
                 "ldPruneByCorrelation")
})

# ===========================================================================
# dropCollinearColumns
# ===========================================================================

# dropCollinearColumns and enforceDesignFullRank are unexported helpers;
# access them via pecotmr::: in these tests.

test_that("dropCollinearColumns returns X unchanged when problematicCols is empty", {
  X <- matrix(rnorm(100), nrow = 20, ncol = 5)
  colnames(X) <- paste0("v", 1:5)
  result <- pecotmr:::dropCollinearColumns(X, problematicCols = character(0), strategy = "correlation")
  expect_equal(ncol(result), 5)
})

test_that("dropCollinearColumns removes single problematic column", {
  X <- matrix(rnorm(100), nrow = 20, ncol = 5)
  colnames(X) <- paste0("v", 1:5)
  result <- pecotmr:::dropCollinearColumns(X, problematicCols = "v3", strategy = "correlation")
  expect_equal(ncol(result), 4)
  expect_false("v3" %in% colnames(result))
})

test_that("dropCollinearColumns variance strategy removes lowest variance column", {
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  colnames(X) <- c("low_var", "mid_var", "high_var")
  X[, 1] <- X[, 1] * 0.01
  X[, 3] <- X[, 3] * 10
  result <- pecotmr:::dropCollinearColumns(X, problematicCols = c("low_var", "mid_var", "high_var"),
                                             strategy = "variance")
  expect_equal(ncol(result), 2)
  expect_false("low_var" %in% colnames(result))
})

test_that("dropCollinearColumns correlation strategy with two columns removes one", {
  set.seed(42)
  X <- matrix(rnorm(100), nrow = 20, ncol = 5)
  colnames(X) <- paste0("v", 1:5)
  result <- pecotmr:::dropCollinearColumns(X, problematicCols = c("v1", "v2"), strategy = "correlation")
  expect_equal(ncol(result), 4)
})

test_that("dropCollinearColumns correlation strategy with 3+ cols removes highest sum", {
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
  colnames(X) <- paste0("v", 1:4)
  X[, 2] <- X[, 1] + rnorm(n, sd = 0.01)
  X[, 3] <- X[, 1] + rnorm(n, sd = 0.01)
  result <- pecotmr:::dropCollinearColumns(X, problematicCols = c("v1", "v2", "v3"),
                                             strategy = "correlation")
  expect_equal(ncol(result), 3)
})

test_that("dropCollinearColumns response_correlation strategy removes lowest |cor| with response", {
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  colnames(X) <- c("v1", "v2", "v3")
  response <- X[, 1] * 2 + rnorm(n, sd = 0.1)
  result <- pecotmr:::dropCollinearColumns(X, problematicCols = c("v1", "v2", "v3"),
                                             strategy = "response_correlation",
                                             response = response)
  expect_equal(ncol(result), 2)
  expect_true("v1" %in% colnames(result))
})

test_that("dropCollinearColumns errors on response_correlation without response", {
  X <- matrix(rnorm(60), 20, 3)
  colnames(X) <- c("v1", "v2", "v3")
  expect_error(
    pecotmr:::dropCollinearColumns(X, problematicCols = c("v1", "v2"),
                                     strategy = "response_correlation"),
    "response"
  )
})

test_that("dropCollinearColumns errors on invalid strategy", {
  X <- matrix(rnorm(100), nrow = 20, ncol = 5)
  colnames(X) <- paste0("v", 1:5)
  expect_error(
    pecotmr:::dropCollinearColumns(X, problematicCols = c("v1", "v2"),
                                     strategy = "invalid_strategy"),
    "arg"
  )
})

test_that("dropCollinearColumns preserves column name when single column remains", {
  set.seed(42)
  X <- matrix(rnorm(40), nrow = 20, ncol = 2)
  colnames(X) <- c("keeper", "removed")
  result <- pecotmr:::dropCollinearColumns(X, problematicCols = "removed", strategy = "correlation")
  expect_equal(ncol(result), 1)
  expect_equal(colnames(result), "keeper")
})

test_that("dropCollinearColumns is silent by default", {
  X <- matrix(rnorm(40), 20, 2)
  colnames(X) <- c("a", "b")
  expect_silent(pecotmr:::dropCollinearColumns(X, problematicCols = "b", strategy = "correlation"))
})

# ===========================================================================
# enforceDesignFullRank (unexported)
# ===========================================================================

test_that("enforceDesignFullRank returns full-rank matrix when already full rank", {
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  colnames(X) <- paste0("v", 1:3)
  C <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  result <- enforceDesignFullRank(X = X, C = C, strategy = "correlation")
  expect_true(is.matrix(result))
  expect_true(ncol(result) >= 1)
})

test_that("enforceDesignFullRank handles rank-deficient design via correlation fallback", {
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
  colnames(X) <- paste0("v", 1:4)
  X[, 4] <- X[, 1] + X[, 2]
  result <- enforceDesignFullRank(X = X, C = NULL, strategy = "correlation")
  design <- cbind(1, result)
  expect_equal(qr(design)$rank, ncol(design))
})

test_that("enforceDesignFullRank preserves colname for single-column input", {
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n), nrow = n, ncol = 1)
  colnames(X) <- "only_snp"
  result <- enforceDesignFullRank(X = X, C = NULL, strategy = "correlation")
  expect_equal(colnames(result), "only_snp")
})

test_that("enforceDesignFullRank is silent by default", {
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n * 3), n, 3)
  colnames(X) <- paste0("v", 1:3)
  expect_silent(enforceDesignFullRank(X = X, C = NULL, strategy = "correlation"))
})

# ===========================================================================
# ldClumpByScore
# ===========================================================================

test_that("ldClumpByScore skips clumping on single-column input", {
  skip_if_not_installed("bigsnpr")
  skip_if_not_installed("bigstatsr")
  set.seed(1)
  X <- matrix(rbinom(100, 2, 0.3), 100, 1)
  colnames(X) <- "chr1:100:A:G"
  keep <- ldClumpByScore(X, score = 1.0, chr = 1L, pos = 100L, r2 = 0.2)
  expect_equal(keep, 1L)
})

test_that("ldClumpByScore validates input lengths", {
  skip_if_not_installed("bigsnpr")
  skip_if_not_installed("bigstatsr")
  set.seed(1)
  X <- matrix(rbinom(100 * 3, 2, 0.3), 100, 3)
  colnames(X) <- paste0("chr1:", seq_len(3) * 1000, ":A:G")
  expect_error(
    ldClumpByScore(X, score = 1:2, chr = rep(1L, 3), pos = seq_len(3) * 1000L),
    "score"
  )
  expect_error(
    ldClumpByScore(X, score = 1:3, chr = rep(1L, 2), pos = seq_len(3) * 1000L),
    "chr and pos"
  )
})

test_that("ldClumpByScore returns indices on real data", {
  skip_if_not_installed("bigsnpr")
  skip_if_not_installed("bigstatsr")
  set.seed(1)
  n <- 500; p <- 20
  X <- matrix(rbinom(n * p, 2, 0.3), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p) * 1000, ":A:G")
  # Introduce perfect LD between variants 1 and 2
  X[, 2] <- X[, 1]
  score <- c(2, 1, runif(p - 2))
  chr <- rep(1L, p)
  pos <- seq_len(p) * 1000L
  keep <- ldClumpByScore(X, score = score, chr = chr, pos = pos, r2 = 0.2)
  expect_true(1L %in% keep)
  expect_false(2L %in% keep)   # pruned: same as variant 1 but lower score
  expect_true(length(keep) < p)
})
