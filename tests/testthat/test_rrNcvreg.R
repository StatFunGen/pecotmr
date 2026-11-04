context("regularized_regression - ncvreg")

# ---- ncvregWeights / scadWeights / mcpWeights ----
test_that("ncvregWeights computes weights with SCAD penalty", {
    skip_if_not_installed("ncvreg")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    y <- X[, 1] * 0.5 + rnorm(n)
    result <- ncvregWeights(X, y, penalty = "SCAD")
    expect_equal(nrow(result), p)
    expect_equal(ncol(result), 1)
    expect_true(all(is.finite(result)))
})

test_that("ncvregWeights warns, zero-pads X's zero-variance columns", {
    skip_if_not_installed("ncvreg")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    X[, 5] <- 3
    y <- X[, 1] * 0.5 + rnorm(n)
    expect_warning(
        result <- ncvregWeights(X, y, penalty = "SCAD"),
        "ncvregWeights: dropping 1 zero-variance column"
    )
    expect_equal(nrow(result), p)
    expect_equal(result[5, 1], 0)
})

test_that("mcpWeights computes weights and dispatches to ncvregWeights", {
    skip_if_not_installed("ncvreg")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    y <- X[, 1] * 0.5 + rnorm(n)
    result <- mcpWeights(X, y)
    expect_equal(nrow(result), p)
    expect_equal(ncol(result), 1)
    expect_true(all(is.finite(result)))
})

test_that("ncvregWeights keeps the selected fit or the path without X", {
    skip_if_not_installed("ncvreg")
    set.seed(42)
    X <- matrix(rnorm(500), nrow = 50)
    colnames(X) <- sprintf("chr1:%d:A:G", seq_len(10))
    y <- X[, 1] * 0.5 + rnorm(50)
    set.seed(1)
    none <- scadWeights(X, y)
    set.seed(1)
    slim <- scadWeights(X, y, fitRetention = "slim")
    expect_equal(as.numeric(slim), as.numeric(none))
    fit <- attr(slim, "fit")
    expect_named(
        fit,
        c("coef", "intercept", "lambda", "cve", "cvse", "variantIds")
    )
    expect_equal(fit$variantIds, colnames(X))
    full <- attr(mcpWeights(X, y, fitRetention = "full"), "fit")
    expect_s3_class(full, "cv.ncvreg")
    # Its standardized copy of the genotypes is not kept.
    expect_null(full$fit$X)
    expect_null(full$variantIds)
    expect_equal(pecotmr:::.twasFitVariantIds(full), colnames(X))
    withX <- attr(
        mcpWeights(X, y, fitRetention = "full", retainX = TRUE),
        "fit"
    )
    expect_true(is.matrix(withX$fit$X))
    expect_equal(pecotmr:::.twasFitVariantIds(withX), colnames(X))
})
