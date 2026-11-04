context("regularized_regression - dpr")

# ---- dprWeights ----
test_that("dprWeights computes weights with Gibbs fitting method", {
    skip_if_not_installed("RcppDPR")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    y <- X[, 1] * 0.5 + rnorm(n)
    result <- dprWeights(X, y, fittingMethod = "Gibbs")
    expect_equal(length(result), p)
    expect_true(is.numeric(result))
    expect_true(all(is.finite(result)))
})

test_that("dprWeights warns and zero-pads when X has zero-variance columns", {
    skip_if_not_installed("RcppDPR")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    X[, 8] <- 6
    y <- X[, 1] * 0.5 + rnorm(n)
    expect_warning(
        result <- dprWeights(X, y),
        "dprWeights: dropping 1 zero-variance column"
    )
    expect_equal(length(result), p)
    expect_equal(result[8], 0)
})

test_that("dprWeights' slim fit drops the copy of the genotypes", {
    skip_if_not_installed("RcppDPR")
    set.seed(42)
    X <- matrix(rnorm(500), nrow = 50)
    colnames(X) <- sprintf("chr1:%d:A:G", seq_len(10))
    y <- X[, 1] * 0.5 + rnorm(50)
    slim <- attr(
        dprWeights(X, y, fittingMethod = "VB", fitRetention = "slim"),
        "fit"
    )
    expect_null(slim$x)
    expect_equal(slim$variantIds, colnames(X))
    full <- attr(
        dprWeights(X, y, fittingMethod = "VB", fitRetention = "full"),
        "fit"
    )
    expect_null(full$x)
    expect_equal(pecotmr:::.twasFitVariantIds(full), colnames(X))
    withX <- attr(
        dprWeights(
            X,
            y,
            fittingMethod = "VB",
            fitRetention = "full",
            retainX = TRUE
        ),
        "fit"
    )
    expect_true(is.matrix(withX$x))
    expect_equal(pecotmr:::.twasFitVariantIds(withX), colnames(X))
})
