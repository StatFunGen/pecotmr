context("regularized_regression - mrash")

# ---- mrashWeights ----
test_that("mrashWeights subsets a full-length beta.init by keep", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    X[, 5] <- 7 # zero-variance column to force subset by keep
    y <- X[, 1] * 0.5 + rnorm(n)
    user_beta_init <- seq(0.01, by = 0.01, length.out = p)
    expect_warning(
        result <- mrashWeights(
            X,
            y,
            methodArgs = MrashOptions(beta.init = user_beta_init)
        ),
        "mrashWeights: dropping 1 zero-variance column"
    )
    expect_equal(length(result), p)
    expect_equal(result[5], 0)
})

test_that("mrashWeights warns and zero-pads when X has zero-variance columns", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    X[, 3] <- 7
    y <- X[, 1] * 0.5 + rnorm(n)
    local_mocked_bindings(
        lassoWeights = function(X, y) rep(0.01, ncol(X))
    )
    expect_warning(
        result <- mrashWeights(X, y),
        "mrashWeights: dropping 1 zero-variance column"
    )
    expect_equal(length(result), p)
    expect_equal(result[3], 0)
})

gc()

test_that("mrashWeights' slim fit drops the copy of the data", {
    set.seed(42)
    X <- matrix(rnorm(500), nrow = 50)
    colnames(X) <- sprintf("chr1:%d:A:G", seq_len(10))
    y <- X[, 1] * 0.5 + rnorm(50)
    set.seed(1)
    slim <- attr(mrashWeights(X, y, fitRetention = "slim"), "fit")
    expect_s3_class(slim, "mr.ash")
    expect_false(any(is_in(c("data", "update.order"), names(slim))))
    # estimateSparsity() reads the mixture weights it keeps.
    expect_true(is.numeric(slim$pi))
    expect_equal(slim$variantIds, colnames(X))
    set.seed(1)
    full <- attr(mrashWeights(X, y, fitRetention = "full"), "fit")
    expect_true(is.list(full$data))
    expect_null(full$variantIds)
    expect_equal(pecotmr:::.twasFitVariantIds(full), colnames(X))
})
