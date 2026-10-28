context("regularized_regression - glmnet")

# ---- glmnetWeights ----
test_that("glmnetWeights computes LASSO weights", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    y <- X[, 1] * 0.5 + rnorm(n)
    result <- glmnetWeights(X, y, alpha = 1)
    expect_equal(nrow(result), p)
    expect_equal(ncol(result), 1)
})

test_that("enetWeights computes elastic net weights", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    y <- X[, 1] * 0.5 + rnorm(n)
    result <- enetWeights(X, y)
    expect_equal(nrow(result), p)
})

test_that("glmnetWeights handles zero-variance columns", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    X[, 3] <- 5
    y <- X[, 1] * 0.5 + rnorm(n)
    expect_warning(
        result <- glmnetWeights(X, y, alpha = 1),
        "glmnetWeights: dropping 1 zero-variance column"
    )
    expect_equal(nrow(result), p)
    expect_equal(result[3, 1], 0)
})

test_that("glmnetWeights errors when all columns are constant", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 50
    p <- 5
    X <- matrix(rep(seq_len(p), each = n), nrow = n, ncol = p)
    y <- rnorm(n)

    expect_error(
        glmnetWeights(X, y, alpha = 1),
        "matrix with 2 or more columns"
    )
})

test_that("glmnetWeights handles NA-producing columns gracefully", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    y <- X[, 1] * 0.5 + rnorm(n)

    X[1, 4] <- NA

    expect_warning(
        result <- glmnetWeights(X, y, alpha = 1),
        "glmnetWeights: dropping 1 zero-variance column"
    )
    expect_equal(nrow(result), p)
    expect_equal(ncol(result), 1)
    expect_equal(result[4, 1], 0)
})

test_that("enetWeights takes its mixing, and ridgeWeights fixes it at 0", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 60
    p <- 8
    X <- matrix(rnorm(n * p), nrow = n)
    y <- X[, 1] * 0.5 + rnorm(n)
    # cv.glmnet draws its folds, so each pair is fitted from the same seed.
    set.seed(11)
    mixed <- enetWeights(X, y, alpha = 0.3)
    set.seed(11)
    expect_equal(mixed, glmnetWeights(X, y, 0.3))
    set.seed(11)
    ridge <- ridgeWeights(X, y)
    set.seed(11)
    expect_equal(ridge, glmnetWeights(X, y, 0))
    expect_equal(nrow(ridge), p)
    # A ridge keeps every variant.
    expect_true(all(ridge != 0))
    expect_error(enetWeights(X, y, alpha = 2), "alpha.*<= 1")
})

test_that("glmnetWeights produces non-zero weights for correlated data", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 100
    p <- 10
    X <- matrix(rnorm(n * p), n, p)
    beta_true <- c(3, -2, rep(0, p - 2))
    y <- X %*% beta_true + rnorm(n)

    w <- glmnetWeights(X, y, alpha = 0.5)
    expect_length(w, p)
    expect_true(any(w != 0))
})
