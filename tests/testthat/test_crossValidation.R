context("crossValidation (shared CV engine)")

# The engine is exercised directly through .crossValidateWeights() with a mock
# per-fold fit, so these tests cover the harness mechanics (partitioning, fold
# aggregation, prediction, metrics, key format, subsampling, threading) without
# any real weight-learning. twasWeightsCv() and .fmWeightsCv() are thin callers
# that supply a domain-specific fitFold; their integration is tested elsewhere.

cv <- function(...) pecotmr:::.crossValidateWeights(...)

# One "mock" method whose weights are all 1s over the training columns, so a
# held-out prediction is the row sum of that sample's (training-column) dosages.
mockFitFold <- function(Xtr, Ytr, j, ...) {
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

mkXY <- function(n = 30, p = 6, k = 1, seed = 1) {
    set.seed(seed)
    X <- matrix(
        rnorm(n * p),
        n,
        p,
        dimnames = list(
            paste0("s", seq_len(n)),
            sprintf("chr1:%d:A:G", 100L * seq_len(p))
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
    d <- mkXY()
    expect_error(
        cv(d$X, d$Y, fold = 0, fitFold = mockFitFold),
        "Must be >= 1"
    )
    expect_error(
        cv(d$X, d$Y, fold = "a", fitFold = mockFitFold),
        "Must be of type 'count'"
    )
    expect_error(
        cv(as.data.frame(d$X), d$Y, fold = 2, fitFold = mockFitFold),
        "Must be of type 'matrix'"
    )
    expect_error(
        cv(d$X, d$Y[1:5, , drop = FALSE], fold = 2, fitFold = mockFitFold),
        "Must have exactly 30 rows"
    )
    expect_error(
        cv(d$X, d$Y, fitFold = mockFitFold),
        "Either 'fold' or 'samplePartitions'"
    )
})

test_that("Y as a vector is converted to a matrix with a message", {
    d <- mkXY(k = 1)
    expect_message(
        cv(
            d$X,
            as.numeric(d$Y),
            fold = 2,
            fitFold = mockFitFold,
            verbose = 1
        ),
        "Y converted to matrix"
    )
})

test_that("output keys use the canonical <method>_predicted / _performance", {
    d <- mkXY()
    r <- suppressMessages(cv(d$X, d$Y, fold = 3, fitFold = mockFitFold))
    expect_equal(names(r$prediction), "mock_predicted")
    expect_equal(names(r$performance), "mock_performance")
})

test_that("performance carries the six metric colnames and outcome rownames", {
    d <- mkXY(k = 2)
    r <- suppressMessages(cv(d$X, d$Y, fold = 3, fitFold = mockFitFold))
    perf <- r$performance[["mock_performance"]]
    expect_equal(
        colnames(perf),
        c("corr", "rsq", "adj_rsq", "pval", "RMSE", "MAE")
    )
    expect_equal(rownames(perf), colnames(d$Y))
})

test_that("every sample is predicted exactly once across folds", {
    d <- mkXY(n = 30)
    r <- suppressMessages(cv(d$X, d$Y, fold = 5, fitFold = mockFitFold))
    pred <- r$prediction[["mock_predicted"]]
    expect_equal(nrow(pred), nrow(d$X))
    expect_false(any(is.na(pred)))
})

test_that("a provided samplePartition is reused; a fold mismatch warns", {
    d <- mkXY(n = 12)
    sp <- data.frame(
        Sample = rownames(d$X),
        Fold = rep(1:3, each = 4),
        stringsAsFactors = FALSE
    )
    r <- suppressMessages(cv(
        d$X,
        d$Y,
        samplePartitions = sp,
        fitFold = mockFitFold
    ))
    expect_equal(r$samplePartition, sp)
    expect_message(
        cv(
            d$X,
            d$Y,
            fold = 2,
            samplePartitions = sp,
            fitFold = mockFitFold,
            verbose = 1
        ),
        "does not match"
    )
})

test_that("a samplePartition with unknown samples errors", {
    d <- mkXY()
    sp <- data.frame(
        Sample = paste0("zzz", 1:5),
        Fold = rep(1:2, length.out = 5)
    )
    expect_error(
        cv(d$X, d$Y, samplePartitions = sp, fitFold = mockFitFold),
        "do not match"
    )
})

test_that("maxNumVariants subsamples variants with a message", {
    d <- mkXY(p = 20)
    expect_message(
        cv(
            d$X,
            d$Y,
            fold = 2,
            fitFold = mockFitFold,
            maxNumVariants = 8,
            verbose = 1
        ),
        "Randomly selecting 8 out of 20"
    )
})

test_that("maxNumVariants with variantsToKeep retains the specified variants", {
    d <- mkXY(p = 20)
    expect_message(
        cv(
            d$X,
            d$Y,
            fold = 2,
            fitFold = mockFitFold,
            maxNumVariants = 8,
            variantsToKeep = c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G"),
            verbose = 1
        ),
        "Including 3 specified variants"
    )
})

test_that("a degenerate fold (empty train/test) is skipped, not errored", {
    d <- mkXY(n = 10)
    sp <- data.frame(Sample = rownames(d$X), Fold = rep(1L, nrow(d$X)))
    r <- suppressMessages(cv(
        d$X,
        d$Y,
        samplePartitions = sp,
        fitFold = mockFitFold
    ))
    expect_true(all(is.na(r$prediction[["mock_predicted"]])))
})

test_that("zero-variance predictions yield NA metrics with a message", {
    d <- mkXY()
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
    d <- mkXY(n = 40, seed = 7)
    set.seed(1)
    r1 <- suppressMessages(cv(
        d$X,
        d$Y,
        fold = 4,
        fitFold = mockFitFold,
        numThreads = 1
    ))
    set.seed(1)
    r2 <- suppressMessages(cv(
        d$X,
        d$Y,
        fold = 4,
        fitFold = mockFitFold,
        numThreads = 2
    ))
    expect_equal(r1$prediction, r2$prediction)
})

test_that("fold fits are collected only when fitRetention asks", {
    d <- mkXY()
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
        fitRetention = "none"
    ))
    expect_true(all(map_int(r_off$foldFits, length) == 0L))
    r_on <- suppressMessages(cv(
        d$X,
        d$Y,
        fold = 3,
        fitFold = fit_with_model,
        fitRetention = "slim"
    ))
    expect_equal(r_on$foldFits[["fold_1"]][["mock"]]$fold, 1L)
})

test_that("X without rownames inherits the sample names (from Y)", {
    d <- mkXY()
    X2 <- d$X
    rownames(X2) <- NULL
    r <- suppressMessages(cv(X2, d$Y, fold = 3, fitFold = mockFitFold))
    expect_equal(rownames(r$prediction$mock_predicted), rownames(d$Y))
})

test_that("maxNumVariants subsamples from variantsToKeep when it already exceeds the cap", {
    d <- mkXY(p = 6)
    expect_message(
        cv(
            d$X,
            d$Y,
            fold = 2,
            fitFold = mockFitFold,
            maxNumVariants = 2,
            variantsToKeep = c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G"),
            verbose = 1
        ),
        "Randomly selecting 2 out of 3"
    )
})

test_that("a NULL per-method weight matrix yields an all-NA prediction, not an error", {
    d <- mkXY()
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


test_that("sample names are synthesized when neither matrix carries rownames", {
    # The CV split is keyed by sample name, so unnamed inputs still need a
    # stable per-row identity rather than falling back to positions.
    X <- matrix(1:6, 3L, 2L)
    Y <- matrix(1:3, 3L, 1L)
    out <- pecotmr:::.cvSetDimnames(X, Y)
    expect_equal(rownames(out$X), c("sample_1", "sample_2", "sample_3"))
    expect_equal(rownames(out$Y), rownames(out$X))
})

test_that("numThreads = -1 asks BiocParallel for the worker count", {
    # -1 means "all available"; anything else is capped at what is available.
    expect_equal(
        pecotmr:::.cvNumCores(-1),
        BiocParallel::multicoreWorkers()
    )
    expect_equal(pecotmr:::.cvNumCores(1), 1)
})

test_that("a seed makes a CV call reproducible without leaking RNG state", {
    # withr::local_seed, not set.seed: the caller's stream must be untouched
    # after the call returns (Bioconductor asks packages not to set.seed).
    f <- function() {
        pecotmr:::.applySeed(42)
        runif(1)
    }
    expect_equal(f(), f())
    set.seed(1)
    before <- runif(1)
    set.seed(1)
    invisible(f())
    expect_equal(runif(1), before)
})

test_that(".cvPartitionKey answers NULL when there is no partition to key", {
    expect_null(pecotmr:::.cvPartitionKey(NULL))
})

test_that("sumstats inputs ignore residualization but refuse crossValidation", {
    # The asymmetry follows from the defaults, not from taste. CV is off
    # unless asked for (folds = 0), so a non-default value on a
    # summary-statistics run is an explicit request for something the input
    # cannot do -> error. Residualization is ON by default, so refusing a
    # non-default would reject the DEFAULT bundle and force every sumstats
    # caller to unset it -> ignore.
    expect_equal(CrossValidationParam()$folds, 0)
    # Residualization's flags are NULL = unset, which the accessors read as
    # "on"; the default bundle therefore pins nothing and is safe to carry
    # onto an input that cannot residualize.
    expect_length(ResidualizationParam(), 0L)
    # CV is refused on a summary-statistics input ...
    expect_error(
        pecotmr:::.cvRefuseOnSumstats(
            CrossValidationParam(folds = 5),
            "fineMappingPipeline",
            "QtlSumStats"
        ),
        "cross-validation"
    )
    # ... while the DEFAULT CV bundle passes, so one bundle still travels to
    # either input kind.
    expect_silent(pecotmr:::.cvRefuseOnSumstats(
        CrossValidationParam(),
        "fineMappingPipeline",
        "QtlSumStats"
    ))
    # Residualization has no equivalent refusal: the sumstats workers do not
    # take it at all, so it cannot reach an engine and cannot error.
    for (fn in c(
        ".fmPipelineQtlSumStats",
        ".fmPipelineGwas",
        ".twasPipelineQtlSumStats"
    )) {
        f <- get(fn, envir = asNamespace("pecotmr"))
        expect_false("residualization" %in% names(formals(f)), label = fn)
    }
})
