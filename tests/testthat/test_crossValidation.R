context("crossValidation (shared CV engine)")

# The engine is exercised directly through .crossValidateWeights() with a mock
# per-fold fit, so these tests cover the harness mechanics (partitioning, fold
# aggregation, prediction, metrics, key format, subsampling, threading) without
# any real weight-learning. twasWeightsCv() and .fmWeightsCv() are thin callers
# that supply a domain-specific fitFold; their integration is tested elsewhere.

cv <- function(...) pecotmr:::.crossValidateWeights(...)

# One "mock" method whose weights are all 1s over the training columns, so a
# held-out prediction is the row sum of that sample's (training-column) dosages.
mock_fit_fold <- function(Xtr, Ytr, j, ...) {
    list(
        weights = list(
            mock = matrix(
                1,
                ncol(Xtr),
                ncol(Ytr),
                dimnames = list(colnames(Xtr), NULL)
            )
        ),
        fits = list()
    )
}

mk_xy <- function(n = 30, p = 6, k = 1, seed = 1) {
    set.seed(seed)
    X <- matrix(
        rnorm(n * p),
        n,
        p,
        dimnames = list(
            paste0("s", seq_len(n)),
            sprintf("chr1:%d:A:G", 100L * (seq_len(p)))
        )
    )
    Y <- matrix(
        rnorm(n * k),
        n,
        k,
        dimnames = list(rownames(X), paste0("c", seq_len(k)))
    )
    list(X = X, Y = Y)
}

test_that("input is validated", {
    d <- mk_xy()
    expect_error(
        cv(d$X, d$Y, fold = 0, fitFold = mock_fit_fold),
        "positive integer"
    )
    expect_error(
        cv(d$X, d$Y, fold = "a", fitFold = mock_fit_fold),
        "positive integer"
    )
    expect_error(
        cv(as.data.frame(d$X), d$Y, fold = 2, fitFold = mock_fit_fold),
        "X must be a matrix"
    )
    expect_error(
        cv(d$X, d$Y[1:5, , drop = FALSE], fold = 2, fitFold = mock_fit_fold),
        "same"
    )
    expect_error(
        cv(d$X, d$Y, fitFold = mock_fit_fold),
        "Either 'fold' or 'samplePartitions'"
    )
})

test_that("Y as a vector is converted to a matrix with a message", {
    d <- mk_xy(k = 1)
    expect_message(
        cv(
            d$X,
            as.numeric(d$Y),
            fold = 2,
            fitFold = mock_fit_fold,
            verbose = 1
        ),
        "Y converted to matrix"
    )
})

test_that("output keys use the canonical <method>_predicted / _performance", {
    d <- mk_xy()
    r <- suppressMessages(cv(d$X, d$Y, fold = 3, fitFold = mock_fit_fold))
    expect_equal(names(r$prediction), "mock_predicted")
    expect_equal(names(r$performance), "mock_performance")
})

test_that("performance carries the six metric colnames and outcome rownames", {
    d <- mk_xy(k = 2)
    r <- suppressMessages(cv(d$X, d$Y, fold = 3, fitFold = mock_fit_fold))
    perf <- r$performance[["mock_performance"]]
    expect_equal(
        colnames(perf),
        c("corr", "rsq", "adj_rsq", "pval", "RMSE", "MAE")
    )
    expect_equal(rownames(perf), colnames(d$Y))
})

test_that("every sample is predicted exactly once across folds", {
    d <- mk_xy(n = 30)
    r <- suppressMessages(cv(d$X, d$Y, fold = 5, fitFold = mock_fit_fold))
    pred <- r$prediction[["mock_predicted"]]
    expect_equal(nrow(pred), nrow(d$X))
    expect_false(any(is.na(pred)))
})

test_that("a provided samplePartition is reused; a fold mismatch warns", {
    d <- mk_xy(n = 12)
    sp <- data.frame(
        Sample = rownames(d$X),
        Fold = rep(1:3, each = 4),
        stringsAsFactors = FALSE
    )
    r <- suppressMessages(cv(
        d$X,
        d$Y,
        samplePartitions = sp,
        fitFold = mock_fit_fold
    ))
    expect_equal(r$samplePartition, sp)
    expect_message(
        cv(
            d$X,
            d$Y,
            fold = 2,
            samplePartitions = sp,
            fitFold = mock_fit_fold,
            verbose = 1
        ),
        "does not match"
    )
})

test_that("a samplePartition with unknown samples errors", {
    d <- mk_xy()
    sp <- data.frame(
        Sample = paste0("zzz", 1:5),
        Fold = rep(1:2, length.out = 5)
    )
    expect_error(
        cv(d$X, d$Y, samplePartitions = sp, fitFold = mock_fit_fold),
        "do not match"
    )
})

test_that("maxNumVariants subsamples variants with a message", {
    d <- mk_xy(p = 20)
    expect_message(
        cv(
            d$X,
            d$Y,
            fold = 2,
            fitFold = mock_fit_fold,
            maxNumVariants = 8,
            verbose = 1
        ),
        "Randomly selecting 8 out of 20"
    )
})

test_that("maxNumVariants with variantsToKeep retains the specified variants", {
    d <- mk_xy(p = 20)
    expect_message(
        cv(
            d$X,
            d$Y,
            fold = 2,
            fitFold = mock_fit_fold,
            maxNumVariants = 8,
            variantsToKeep = c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G"),
            verbose = 1
        ),
        "Including 3 specified variants"
    )
})

test_that("a degenerate fold (empty train/test) is skipped, not errored", {
    d <- mk_xy(n = 10)
    sp <- data.frame(Sample = rownames(d$X), Fold = rep(1L, nrow(d$X)))
    r <- suppressMessages(cv(
        d$X,
        d$Y,
        samplePartitions = sp,
        fitFold = mock_fit_fold
    ))
    expect_true(all(is.na(r$prediction[["mock_predicted"]])))
})

test_that("zero-variance predictions yield NA metrics with a message", {
    d <- mk_xy()
    zero_fit <- function(Xtr, Ytr, j, ...) {
        list(
            weights = list(
                mock = matrix(
                    0,
                    ncol(Xtr),
                    ncol(Ytr),
                    dimnames = list(colnames(Xtr), NULL)
                )
            ),
            fits = list()
        )
    }
    expect_message(
        r <- cv(d$X, d$Y, fold = 3, fitFold = zero_fit, verbose = 1),
        "zero variance"
    )
    expect_true(all(is.na(r$performance[["mock_performance"]])))
})

test_that("the parallel fold path matches the serial one", {
    d <- mk_xy(n = 40, seed = 7)
    set.seed(1)
    r1 <- suppressMessages(cv(
        d$X,
        d$Y,
        fold = 4,
        fitFold = mock_fit_fold,
        numThreads = 1
    ))
    set.seed(1)
    r2 <- suppressMessages(cv(
        d$X,
        d$Y,
        fold = 4,
        fitFold = mock_fit_fold,
        numThreads = 2
    ))
    expect_equal(r1$prediction, r2$prediction)
})

test_that("retainFits collects per-fold fits only when requested", {
    d <- mk_xy()
    fit_with_model <- function(Xtr, Ytr, j, ...) {
        list(
            weights = list(
                mock = matrix(
                    1,
                    ncol(Xtr),
                    ncol(Ytr),
                    dimnames = list(colnames(Xtr), NULL)
                )
            ),
            fits = list(mock = list(fold = j))
        )
    }
    r_off <- suppressMessages(cv(
        d$X,
        d$Y,
        fold = 3,
        fitFold = fit_with_model,
        retainFits = FALSE
    ))
    expect_true(all(vapply(r_off$foldFits, length, integer(1)) == 0L))
    r_on <- suppressMessages(cv(
        d$X,
        d$Y,
        fold = 3,
        fitFold = fit_with_model,
        retainFits = TRUE
    ))
    expect_equal(r_on$foldFits[["fold_1"]][["mock"]]$fold, 1L)
})

test_that("X without rownames inherits the sample names (from Y)", {
    d <- mk_xy()
    X2 <- d$X
    rownames(X2) <- NULL
    r <- suppressMessages(cv(X2, d$Y, fold = 3, fitFold = mock_fit_fold))
    expect_equal(rownames(r$prediction$mock_predicted), rownames(d$Y))
})

test_that("maxNumVariants subsamples from variantsToKeep when it already exceeds the cap", {
    d <- mk_xy(p = 6)
    expect_message(
        cv(
            d$X,
            d$Y,
            fold = 2,
            fitFold = mock_fit_fold,
            maxNumVariants = 2,
            variantsToKeep = c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G"),
            verbose = 1
        ),
        "Randomly selecting 2 out of 3"
    )
})

test_that("a NULL per-method weight matrix yields an all-NA prediction, not an error", {
    d <- mk_xy()
    fit_with_null <- function(Xtr, Ytr, j, ...) {
        list(
            weights = list(
                mock = matrix(
                    1,
                    ncol(Xtr),
                    ncol(Ytr),
                    dimnames = list(colnames(Xtr), NULL)
                ),
                empty = NULL
            ), # NULL W -> skipped per fold
            fits = list()
        )
    }
    r <- suppressMessages(cv(d$X, d$Y, fold = 3, fitFold = fit_with_null))
    expect_true(all(is.na(r$prediction$empty_predicted)))
    expect_false(all(is.na(r$prediction$mock_predicted)))
})
