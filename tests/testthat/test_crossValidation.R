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

test_that("the CV codec leaves shapes it does not split as they are", {
    ids <- c("s1", "s2", "s3")
    expect_identical(
        pecotmr:::.cvUnnameRows(c(0.1, 0.2, 0.3), ids),
        c(0.1, 0.2, 0.3)
    )
    named <- matrix(1:3, 3L, 1L, dimnames = list(ids, NULL))
    expect_identical(pecotmr:::.cvNameRows(named, ids), named)
    # A TWAS matrix of predictions, one column per condition, keeps its row
    # names, whatever its width.
    for (k in 1:2) {
        cv <- list(
            predictions = matrix(seq_len(3L * k), 3L, k, dimnames = list(ids))
        )
        expect_identical(pecotmr:::.cvNamePredictions(cv, ids), cv)
    }
    # A collection saved before partitions were shared has no dictionary.
    data(qtlFineMappingExample)
    old <- qtlFineMappingExample
    attr(old, "cvPartition") <- NULL
    expect_identical(pecotmr:::.cvPartitionDict(old), list())
})

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
    # One fold holding every sample leaves nothing to train on: the fold is
    # skipped before any fit, so nothing is predicted or scored.
    d <- mkXY(n = 10)
    sp <- data.frame(Sample = rownames(d$X), Fold = rep(1L, nrow(d$X)))
    fits <- 0L
    countingFit <- function(...) {
        fits <<- fits + 1L
        mockFitFold(...)
    }
    r <- suppressMessages(cv(
        d$X,
        d$Y,
        samplePartitions = sp,
        fitFold = countingFit
    ))
    expect_equal(fits, 0L)
    expect_length(r$prediction, 0L)
    expect_length(r$performance, 0L)
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
    skip_on_os("windows")
    d <- mkXY(n = 40, seed = 7)
    set.seed(1)
    r1 <- suppressMessages(cv(
        d$X,
        d$Y,
        fold = 4,
        fitFold = mockFitFold,
        BPPARAM = BiocParallel::SerialParam()
    ))
    set.seed(1)
    r2 <- suppressMessages(cv(
        d$X,
        d$Y,
        fold = 4,
        fitFold = mockFitFold,
        BPPARAM = BiocParallel::MulticoreParam(2)
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

test_that("maxNumVariants subsamples variantsToKeep when it exceeds the cap", {
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

test_that("a NULL method weight matrix predicts all NA, not an error", {
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

test_that("BPPARAM must be a BiocParallelParam", {
    expect_error(pecotmr:::.assertBpparam(2), "must be a BiocParallelParam")
    expect_error(
        CrossValidationParam(BPPARAM = 4),
        "must be a BiocParallelParam"
    )
})

test_that("the workers draw from a seeded copy of the caller's BPPARAM", {
    bp <- BiocParallel::SerialParam()
    # `seed` first, then the param's own RNGseed, and the caller's param is
    # left as it was.
    expect_identical(
        BiocParallel::bpRNGseed(pecotmr:::.seededBpparam(bp, 7L)),
        7L
    )
    expect_null(BiocParallel::bpRNGseed(bp))
    own <- BiocParallel::SerialParam(RNGseed = 3L)
    expect_identical(
        BiocParallel::bpRNGseed(pecotmr:::.seededBpparam(own)),
        3L
    )
    # Neither: one draw from the session RNG, so set.seed() still governs.
    set.seed(11)
    a <- BiocParallel::bpRNGseed(pecotmr:::.seededBpparam(bp))
    set.seed(11)
    b <- BiocParallel::bpRNGseed(pecotmr:::.seededBpparam(bp))
    expect_identical(a, b)
})

test_that("a worker's error is raised with its own message", {
    boom <- function(i) {
        if (i == 2L) {
            rlang::abort("element two failed")
        }
        i
    }
    expect_error(
        pecotmr:::.bpMap(1:3, boom, BPPARAM = BiocParallel::SerialParam()),
        "^element two failed$"
    )
    expect_identical(
        pecotmr:::.bpMap(1:2, identity, BPPARAM = BiocParallel::SerialParam()),
        list(1L, 2L)
    )
})

test_that("a threaded engine counts only workers on this machine", {
    expect_identical(
        pecotmr:::.bpLocalThreads(BiocParallel::SerialParam(), "f"),
        1L
    )
    expect_identical(
        pecotmr:::.bpLocalThreads(BiocParallel::MulticoreParam(3), "f"),
        3L
    )
    rlang::local_options(rlib_message_verbosity = "verbose")
    expect_message(
        n <- pecotmr:::.bpLocalThreads(BiocParallel::SnowParam(4), "f"),
        "f: the workers of a SnowParam need not run on this machine"
    )
    expect_identical(n, 1L)
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

test_that(".cvPartitionKey ignores row order but not fold assignment", {
    sp <- data.frame(Sample = c("s1", "s2", "s3"), Fold = c(1L, 2L, 1L))
    shuffled <- sp[c(3, 1, 2), ]
    expect_identical(
        pecotmr:::.cvPartitionKey(sp),
        pecotmr:::.cvPartitionKey(shuffled)
    )
    moved <- sp
    moved$Fold <- c(1L, 1L, 2L)
    expect_false(identical(
        pecotmr:::.cvPartitionKey(sp),
        pecotmr:::.cvPartitionKey(moved)
    ))
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

test_that("CrossValidationParam(fitFullData) is validated at construction", {
    expect_true(CrossValidationParam()$fitFullData)
    expect_false(
        CrossValidationParam(folds = 2, fitFullData = FALSE)$fitFullData
    )
    # FALSE keeps only the cross-validation results, so there must be
    # cross-validation to keep: 0 and 1 folds both mean "none".
    expect_error(
        CrossValidationParam(fitFullData = FALSE),
        "fitFullData = FALSE requires cross-validation"
    )
    expect_error(
        CrossValidationParam(folds = 1, fitFullData = FALSE),
        "fitFullData = FALSE requires cross-validation"
    )
    expect_error(
        CrossValidationParam(fitFullData = NA),
        "single TRUE or FALSE"
    )
    expect_error(
        CrossValidationParam(fitFullData = c(TRUE, FALSE)),
        "single TRUE or FALSE"
    )
    expect_true(pecotmr:::.cvResolve(CrossValidationParam())$fitFullData)
    # The default is not a request, so the default bundle still passes the
    # summary-statistics refusal; an explicit FALSE is one, and is refused.
    expect_silent(pecotmr:::.cvRefuseOnSumstats(
        CrossValidationParam(),
        "twasWeightsPipeline",
        "QtlSumStats"
    ))
    expect_error(
        pecotmr:::.cvRefuseOnSumstats(
            CrossValidationParam(folds = 2, fitFullData = FALSE),
            "twasWeightsPipeline",
            "QtlSumStats"
        ),
        "fitFullData"
    )
})

test_that("a BPPARAM other than SerialParam() is a cross-validation request", {
    expect_silent(pecotmr:::.cvRefuseOnSumstats(
        CrossValidationParam(BPPARAM = BiocParallel::SerialParam()),
        "fineMappingPipeline",
        "QtlSumStats"
    ))
    expect_error(
        pecotmr:::.cvRefuseOnSumstats(
            CrossValidationParam(BPPARAM = BiocParallel::MulticoreParam(2)),
            "fineMappingPipeline",
            "QtlSumStats"
        ),
        "Remove CrossValidationParam(BPPARAM)",
        fixed = TRUE
    )
})

test_that("a BiocParallel report without a worker error is raised as is", {
    cnd <- structure(
        class = c("bplist_error", "bperror", "error", "condition"),
        list(message = "no element failed", call = NULL)
    )
    attr(cnd, "result") <- list(1L)
    expect_error(pecotmr:::.bpRaiseFirstError(cnd), "no element failed")
})

# ===========================================================================
# .cvSamplePartition and .cvMetricRow
# ===========================================================================

test_that(".cvSamplePartition puts every sample in the requested folds", {
    part <- pecotmr:::.cvSamplePartition(paste0("s", 1:20), fold = 4L)
    expect_setequal(part$Sample, paste0("s", 1:20))
    expect_setequal(sort(unique(part$Fold)), 1:4)
    expect_equal(nrow(part), 20L)
})

test_that(".cvMetricRow: < 3 usable predictions -> all-NA row", {
    expect_true(all(is.na(pecotmr:::.cvMetricRow(c(1, 2), c(1, 2)))))
    ok <- pecotmr:::.cvMetricRow(c(1, 2, 3, 4, 5), c(1.1, 2, 2.9, 4, 5))
    expect_false(is.na(ok[["rsq"]]))
})

# ---- CV results in a collection ----------------------------------------------
# A CV run's rows share its fold partition and the order its predictions are
# held in; a collection keeps each such pair once (.cvStoredColumn()).

.cvp_partition <- function() {
    tibble::tibble(
        Sample = c("s3", "s1", "s4", "s2"),
        Fold = c(1L, 1L, 2L, 2L)
    )
}

# One method's fine-mapping CV payload (its `prediction` matrices name the
# samples by row), as fineMappingPipeline() attaches it to a row.
.cvp_fm <- function(key = "susie", sp = .cvp_partition()) {
    pred <- matrix(
        c(0.1, 0.2, 0.3, 0.4),
        4L,
        1L,
        dimnames = list(c("s1", "s2", "s3", "s4"), "y1")
    )
    list(
        samplePartition = sp,
        prediction = set_names(list(pred), str_c(key, "_predicted")),
        performance = set_names(
            list(matrix(0.5, 1L, 6L)),
            str_c(key, "_performance")
        )
    )
}

# One condition's TWAS CV payload (its `predictions` vector names them).
.cvp_twas <- function(sp = .cvp_partition()) {
    list(
        samplePartition = sp,
        predictions = c(s1 = 0.1, s2 = 0.2, s3 = 0.3, s4 = 0.4),
        metrics = c(corr = 0.5, rsq = 0.25),
        foldFits = NULL
    )
}

test_that("a CV payload splits from its partition and joins back as it was", {
    for (payload in list(.cvp_fm(), .cvp_twas())) {
        s <- pecotmr:::.cvSplitPartition(payload)
        expect_null(s$stored$samplePartition)
        expect_identical(s$partition$samplePartition, payload$samplePartition)
        expect_identical(s$partition$sampleIds, c("s1", "s2", "s3", "s4"))
        expect_identical(
            pecotmr:::.cvJoinPartition(s$stored, s$partition),
            payload
        )
    }
    # Stored without the sample names the partition entry holds.
    fm <- pecotmr:::.cvSplitPartition(.cvp_fm())$stored
    expect_null(rownames(fm$prediction$susie_predicted))
    expect_equal(colnames(fm$prediction$susie_predicted), "y1")
    twas <- pecotmr:::.cvSplitPartition(.cvp_twas())$stored
    expect_null(names(twas$predictions))
})

test_that("a split keeps the samples outside every fold, in place", {
    # A supplied partition may leave samples out; their predictions stay.
    payload <- .cvp_fm(sp = .cvp_partition()[1:3, ])
    s <- pecotmr:::.cvSplitPartition(payload)
    expect_identical(s$partition$sampleIds, c("s1", "s2", "s3", "s4"))
    expect_identical(pecotmr:::.cvJoinPartition(s$stored, s$partition), payload)
})

test_that("a merged row's regions split and join one by one", {
    payload <- list(region1 = .cvp_fm(), region2 = NULL)
    s <- pecotmr:::.cvSplitPartition(payload)
    expect_named(s$partition, c("region1", "region2"))
    expect_null(s$partition$region2)
    expect_identical(pecotmr:::.cvJoinPartition(s$stored, s$partition), payload)
})

test_that("a payload with no partition is stored as it is", {
    for (payload in list(NULL, list(rsq = 0.5), list(samplePartition = NULL))) {
        s <- pecotmr:::.cvSplitPartition(payload)
        expect_identical(s$stored, payload)
        expect_null(s$partition)
        expect_identical(pecotmr:::.cvJoinPartition(payload, NULL), payload)
    }
})

# Two fine-mapping rows of one CV run: they share its partition.
.cvp_fmr <- function(sp = .cvp_partition()) {
    tl <- data.frame(variant_id = "chr1:100:A:G", pip = 0.5)
    QtlFineMappingResult(
        studyName = c("S", "S"),
        context = c("C", "C"),
        trait = c("T", "T"),
        method = c("susie", "susieInf"),
        entry = list(
            FineMappingRow("chr1:100:A:G", list(), tl, .cvp_fm("susie", sp)),
            FineMappingRow(
                "chr1:100:A:G",
                list(),
                tl,
                .cvp_fm("susie_inf", sp)
            )
        )
    )
}

.cvp_cvOf <- function(x, method) {
    cvResult(x, studyName = "S", context = "C", trait = "T", method = method)
}

test_that("a collection keeps each CV partition once; its rows get it back", {
    fmr <- .cvp_fmr()
    expect_length(pecotmr:::.cvPartitionDict(fmr), 1L)
    expect_equal(fmr$cvPartition, c("cv1", "cv1"))
    stored <- fmr$cvResult[[1L]]
    expect_null(stored$samplePartition)
    expect_null(rownames(stored$prediction$susie_predicted))
    expect_identical(.cvp_cvOf(fmr, "susie"), .cvp_fm("susie"))
    # A subset keeps it, and a combine keeps one copy of what both share.
    expect_identical(
        cvResult(fmr[2L], method = "susieInf"),
        .cvp_fm("susie_inf")
    )
    both <- c(fmr[1L], fmr[2L])
    expect_length(pecotmr:::.cvPartitionDict(both), 1L)
    expect_identical(.cvp_cvOf(both, "susieInf"), .cvp_fm("susie_inf"))
    # A different partition is a second entry.
    other <- .cvp_fmr(sp = .cvp_partition()[c(2L, 1L, 3L, 4L), ])
    expect_length(pecotmr:::.cvPartitionDict(c(fmr[1L], other[2L])), 2L)
    # Every row's key must name an entry.
    expect_error(
        {
            fmr$cvPartition <- c("cv1", "cv9")
            validObject(fmr)
        },
        "cv9"
    )
})

test_that("TwasWeights keeps each CV partition once; its rows get it back", {
    rows <- list(
        TwasWeightsRow("chr1:100:A:G", 0.1, cvResult = .cvp_twas()),
        TwasWeightsRow("chr1:100:A:G", 0.2, cvResult = .cvp_twas())
    )
    tw <- TwasWeights(
        studyName = c("S", "S"),
        context = c("C", "C"),
        trait = c("T", "T"),
        method = c("lasso", "enet"),
        entry = rows
    )
    expect_length(pecotmr:::.cvPartitionDict(tw), 1L)
    expect_null(names(tw$cvResult[[1L]]$predictions))
    expect_identical(.cvp_cvOf(tw, "enet"), .cvp_twas())
})

# Row `i`'s CV as the collection hands it out.
.cvp_rowCv <- function(i, x) {
    cvResult(pecotmr:::.collectionEntry(x, i))
}

test_that("updateObject() keeps a saved collection's CV partitions once", {
    fmr <- .cvp_fmr()
    tw <- TwasWeights(
        studyName = "S",
        context = "C",
        trait = "T",
        method = "lasso",
        entry = list(
            TwasWeightsRow("chr1:100:A:G", 0.1, cvResult = .cvp_twas())
        )
    )
    # As saved before rows shared their partitions: each row's own, in full,
    # and no key or dictionary. A TwasWeights was then saved one row per
    # method, which already holds each row's full cvResult.
    full <- map(seq_along(fmr), .cvp_rowCv, x = fmr)
    savedFmr <- fmr
    savedFmr$cvPartition <- NULL
    savedFmr$cvResult <- S4Vectors::SimpleList(full)
    attr(savedFmr, "cvPartition") <- NULL
    expect_identical(map(seq_along(savedFmr), .cvp_rowCv, x = savedFmr), full)
    savedTw <- .savedTwasPerMethod(tw)
    attr(savedTw, "cvPartition") <- NULL
    for (case in list(list(fmr, savedFmr), list(tw, savedTw))) {
        up <- updateObject(case[[2L]])
        expect_identical(up, case[[1L]])
        expect_identical(updateObject(up), up)
    }
})
