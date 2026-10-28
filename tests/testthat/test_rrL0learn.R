context("regularized_regression - l0learn")

# ---- l0learnWeights ----
test_that("l0learnWeights warns and zero-pads when X has zero-variance columns", {
    skip_if_not_installed("L0Learn")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    X[, 6] <- 4
    y <- X[, 1] * 0.5 + rnorm(n)
    expect_warning(
        result <- l0learnWeights(X, y),
        "l0learnWeights: dropping 1 zero-variance column"
    )
    expect_equal(nrow(result), p)
    expect_equal(result[6, 1], 0)
})

test_that("l0learnWeights passes nGamma through to L0Learn.cvfit", {
    skip_if_not_installed("L0Learn")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    y <- X[, 1] * 0.5 + rnorm(n)
    # nGamma is not a formal of l0learnWeights; it reaches the engine through
    # methodArgs.
    captured <- new.env(parent = emptyenv())
    local_mocked_bindings(
        L0Learn.cvfit = function(...) {
            captured$args <- list(...)
            stop("STOP_AFTER_CAPTURE")
        },
        .package = "L0Learn"
    )
    expect_error(
        l0learnWeights(
            X,
            y,
            penalty = "L0L2",
            methodArgs = L0learnOptions(nGamma = 3)
        ),
        "STOP_AFTER_CAPTURE"
    )
    expect_equal(captured$args$nGamma, 3)
    expect_equal(captured$args$penalty, "L0L2")
})
