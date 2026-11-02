context("regularized_regression - bglr")

# ---- bglrWeights, and the alphabet methods it backs ----
test_that("bglrWeights computes weights with BayesB model", {
    skip_if_not_installed("BGLR")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    y <- X[, 1] * 0.5 + rnorm(n)
    result <- bglrWeights(
        X,
        y,
        model = "BayesB",
        nIter = 100,
        burnIn = 20,
        thin = 2,
        etaArgs = list(probIn = 0.05)
    )
    expect_equal(length(result), p)
    expect_true(is.numeric(result))
    expect_true(all(is.finite(result)))
})

test_that("BGLR's BayesB and Bayesian LASSO compute finite length-p weights", {
    skip_if_not_installed("BGLR")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    y <- X[, 1] * 0.5 + rnorm(n)
    fits <- list(
        bayesB = bayesBWeights(X, y, nIter = 100, nBurnin = 20, thin = 2),
        bayesL = bayesLWeights(
            X,
            y,
            backend = "bglr",
            nIter = 100,
            nBurnin = 20,
            thin = 2
        )
    )
    for (nm in names(fits)) {
        expect_equal(length(fits[[nm]]), p, info = nm)
        expect_true(all(is.finite(fits[[nm]])), info = nm)
    }
})

test_that("bayesBWeights passes pi to BGLR's ETA as probIn", {
    skip_if_not_installed("BGLR")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    y <- X[, 1] * 0.5 + rnorm(n)

    captured <- new.env(parent = emptyenv())
    local_mocked_bindings(
        BGLR = function(y, ETA, ...) {
            captured$eta <- ETA
            list(ETA = list(list(b = rep(0, ncol(ETA[[1]]$X)))))
        },
        .package = "BGLR"
    )

    result <- bayesBWeights(
        X,
        y,
        nIter = 100,
        nBurnin = 20,
        thin = 2,
        pi = 0.42
    )
    expect_equal(length(result), p)
    expect_equal(captured$eta[[1]]$model, "BayesB")
    expect_equal(captured$eta[[1]]$probIn, 0.42)
    # probIn is the only ETA setting bayesBWeights adds.
    expect_named(captured$eta[[1]], c("X", "model", "probIn"))
})

test_that("bayesBWeights warns and zero-pads zero-variance columns of X", {
    skip_if_not_installed("BGLR")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    X[, 7] <- 9
    y <- X[, 1] * 0.5 + rnorm(n)
    # The dispatcher drops them before any engine sees X.
    expect_warning(
        result <- bayesBWeights(X, y, nIter = 100, nBurnin = 20, thin = 2),
        "bayesAlphabetWeights: dropping 1 zero-variance column"
    )
    expect_equal(length(result), p)
    expect_equal(result[7], 0)
})

test_that("bglrWeights cleans up its tempdir on exit", {
    skip_if_not_installed("BGLR")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    y <- X[, 1] * 0.5 + rnorm(n)
    before <- list.files(tempdir(), pattern = "^bglr_")
    bayesBWeights(X, y, nIter = 100, nBurnin = 20, thin = 2)
    after <- list.files(tempdir(), pattern = "^bglr_")
    expect_setequal(before, after)
})

test_that("bglrWeights keeps BGLR's fitted object, naming its variants", {
    skip_if_not_installed("BGLR")
    set.seed(42)
    X <- matrix(rnorm(500), nrow = 50)
    colnames(X) <- sprintf("chr1:%d:A:G", seq_len(10))
    y <- X[, 1] * 0.5 + rnorm(50)
    w <- bglrWeights(
        X,
        y,
        model = "BayesB",
        nIter = 100,
        burnIn = 20,
        thin = 5,
        fitRetention = "slim"
    )
    fit <- attr(w, "fit")
    expect_s3_class(fit, "BGLR")
    expect_equal(as.numeric(w), as.numeric(fit$ETA[[1L]]$b))
    expect_null(fit$variantIds)
    expect_equal(pecotmr:::.twasFitVariantIds(fit), colnames(X))
})
