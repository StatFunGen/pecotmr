context("twasWeights")

# ---------------------------------------------------------------------------
# Shared synthetic data generator
# ---------------------------------------------------------------------------
makeData <- function(n = 50, p = 10, seed = 42, add_zero_var_col = FALSE) {
    set.seed(seed)
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * seq_len(p))
    rownames(X) <- paste0("sample_", seq_len(n))

    beta <- rep(0, p)
    beta[1:3] <- c(1.5, -0.8, 0.5)
    noise <- rnorm(n, sd = 0.5)
    Y <- X %*% beta + noise
    Y <- matrix(Y, ncol = 1)
    colnames(Y) <- "outcome_1"
    rownames(Y) <- rownames(X)

    if (add_zero_var_col) {
        # Append a constant column (zero variance)
        X <- cbind(X, zero_var = rep(7, n))
        colnames(X)[p + 1] <- "chr1:9900:A:G"
    }

    list(X = X, Y = Y, beta = beta)
}

makeFakeSusieFit <- function(p = 10, L = 3, inf = FALSE) {
    fit <- list(
        alpha = matrix(1 / p, nrow = L, ncol = p),
        mu = matrix(0, nrow = L, ncol = p),
        lbf_variable = matrix(0, nrow = L, ncol = p),
        X_column_scale_factors = rep(1, p),
        pip = rep(0.1, p),
        V = rep(0.5, L),
        sets = list(cs = NULL, purity = NULL)
    )
    if (inf) {
        fit$theta <- rep(0, p)
    }
    fit
}

mockSusie <- function(...) {
    args <- list(...)
    L <- if (is.null(args$L)) 3 else args$L
    makeFakeSusieFit(
        ncol(args$X),
        L = L,
        inf = identical(args$unmappable_effects, "inf")
    )
}

# Test helper: fetch the weights matrix for a given method token from a
# TwasWeights collection. Accepts either the short token ("lasso") or the
# legacy suffixed name ("lassoWeights" / "lasso_weights"). Single-outcome
# entries store a bare numeric vector internally (drop()'d from the
# learnTwasWeights matrix); promote those back to a 1-column matrix here
# so test assertions on nrow/ncol/rownames keep working.
.weightsByMethod <- function(tw, method) {
    shortName <- sub("_?[Ww]eights$", "", method)
    idx <- which(as.character(tw$method) == shortName)
    if (length(idx) == 0L) {
        return(NULL)
    }
    # The entry is a derived view now: ask the collection for it rather than
    # reaching into a stored `entry` column that no longer exists.
    entry <- getTwasWeights(
        tw,
        study = as.character(tw$study)[[idx[[1L]]]],
        context = as.character(tw$context)[[idx[[1L]]]],
        trait = as.character(tw$trait)[[idx[[1L]]]],
        method = shortName
    )
    w <- getWeights(entry)
    vids <- getVariantIds(entry)
    if (is.numeric(w) && is.null(dim(w))) {
        nm <- names(w)
        if (is.null(nm) && length(vids) == length(w)) {
            nm <- vids
        }
        w <- matrix(w, ncol = 1L, dimnames = list(nm, NULL))
    }
    w
}

# ===========================================================================
#
#  .twas_method_lookup
#
# ===========================================================================

# Simulate a small (X, Y) pair. Local to this file: it was a helper-*.R
# shared with the causal-inference tests, which never adopted it.
generateXY <- function(
    seed = 1,
    numSamples = 10,
    numFeatures = 10,
    xRownames = TRUE,
    yRownames = TRUE
) {
    set.seed(seed)
    X <- scale(
        matrix(rnorm(numSamples * numFeatures), nrow = numSamples),
        center = TRUE,
        scale = TRUE
    )

    if (xRownames) {
        rownames(X) <- paste0("sample", seq_len(numSamples))
    } else {
        rownames(X) <- NULL
    }

    beta <- rep(0, numFeatures)
    beta[1:4] <- 1
    y <- X %*% beta + rnorm(numSamples)
    y <- matrix(y, nrow = numSamples, ncol = 1)
    if (yRownames) {
        rownames(y) <- paste0("sample", seq_len(numSamples))
    } else {
        rownames(y) <- NULL
    }
    colnames(y) <- c("Outcome")

    return(list(X = X, Y = y))
}

test_that(".twas_method_lookup: 'default' preset returns 10 methods", {
    result <- pecotmr:::.twasMethodLookup("default")
    expected_names <- c(
        "susie_weights",
        "susie_inf_weights",
        "mrash_weights",
        "enet_weights",
        "lasso_weights",
        "mcp_weights",
        "scad_weights",
        "l0learn_weights",
        "bayes_r_weights",
        "bayes_c_weights"
    )
    expect_equal(sort(names(result)), sort(expected_names))
})

test_that(".twas_method_lookup: 'fast_default' preset returns 8 methods", {
    result <- pecotmr:::.twasMethodLookup("fastDefault")
    expected_names <- c(
        "susie_weights",
        "susie_inf_weights",
        "mrash_weights",
        "enet_weights",
        "lasso_weights",
        "mcp_weights",
        "scad_weights",
        "l0learn_weights"
    )
    expect_equal(sort(names(result)), sort(expected_names))
})

test_that(".twas_method_lookup: custom vector of short names", {
    result <- pecotmr:::.twasMethodLookup(c("susie", "enet", "dprVb"))
    expect_equal(
        sort(names(result)),
        sort(c("susie_weights", "enet_weights", "dpr_vb_weights"))
    )
})

test_that(".twas_method_lookup: unknown method produces error", {
    expect_error(
        pecotmr:::.twasMethodLookup(c("susie", "nonexistent_method")),
        "unknown method token"
    )
})

test_that(".twas_method_lookup: default args are set for mrash, not susie", {
    result <- pecotmr:::.twasMethodLookup("fastDefault")
    # susie carries NO fitting defaults: susieWeights extracts from a supplied
    # fit and never runs susie, so `refine` / `L` would have nothing to configure.
    expect_length(result$susie_weights, 0L)
    expect_equal(result$mrash_weights$initPriorSd, TRUE)
    expect_equal(result$mrash_weights$max.iter, 100)
})

test_that(".twas_method_lookup: methods with no special args get empty list", {
    result <- pecotmr:::.twasMethodLookup(c("enet", "lasso"))
    expect_equal(length(result$enet_weights), 0L)
    expect_equal(length(result$lasso_weights), 0L)
})

test_that(".twas_method_lookup: all DPR variants can coexist", {
    result <- pecotmr:::.twasMethodLookup(c(
        "dprVb",
        "dprGibbs",
        "dprAdaptiveGibbs"
    ))
    expect_equal(
        sort(names(result)),
        sort(c(
            "dpr_vb_weights",
            "dpr_gibbs_weights",
            "dpr_adaptive_gibbs_weights"
        ))
    )
})

# ===========================================================================
#
#  twasPredict
#
# ===========================================================================

test_that("twasPredict: basic matrix multiplication is correct", {
    d <- makeData(n = 20, p = 5)
    set.seed(99)
    w <- matrix(runif(5), ncol = 1)
    rownames(w) <- colnames(d$X)
    wl <- list(test_weights = w)
    res <- twasPredict(d$X, wl)

    expected <- d$X %*% w
    expect_equal(res[["test_predicted"]], expected)
})

test_that("twasPredict: multiple weight methods in list", {
    d <- makeData(n = 20, p = 5)
    w1 <- matrix(c(1, 0, 0, 0, 0), ncol = 1)
    w2 <- matrix(c(0, 0, 0, 0, 1), ncol = 1)
    wl <- list(method_a_weights = w1, method_b_weights = w2)

    res <- twasPredict(d$X, wl)

    expect_length(res, 2)
    expect_equal(res[["method_a_predicted"]], d$X %*% w1)
    expect_equal(res[["method_b_predicted"]], d$X %*% w2)
})

test_that("twasPredict: name transformation weights -> predicted", {
    set.seed(42)
    wl <- list(
        lassoWeights = matrix(1, nrow = 3, ncol = 1),
        enetWeights = matrix(1, nrow = 3, ncol = 1),
        susieWeights = matrix(1, nrow = 3, ncol = 1)
    )
    X <- matrix(rnorm(9), nrow = 3, ncol = 3)
    res <- twasPredict(X, wl)

    expect_equal(
        names(res),
        c("lassoPredicted", "enetPredicted", "susiePredicted")
    )
})

test_that("twasPredict: names without _weights suffix are kept unchanged", {
    wl <- list(custom_method = matrix(1, nrow = 2, ncol = 1))
    X <- matrix(1:4, nrow = 2, ncol = 2)
    res <- twasPredict(X, wl)

    # gsub("_weights", "_predicted", "custom_method") == "custom_method"
    expect_equal(names(res), "custom_method")
})

test_that("twasPredict: single column Y dimension preserved", {
    d <- makeData(n = 10, p = 4)
    w <- matrix(rep(0.25, 4), ncol = 1)
    wl <- list(avg_weights = w)
    res <- twasPredict(d$X, wl)

    expect_true(is.matrix(res[["avg_predicted"]]))
    expect_equal(nrow(res[["avg_predicted"]]), 10)
    expect_equal(ncol(res[["avg_predicted"]]), 1)
})

test_that("twasPredict: multi-column weights produce multi-column predictions", {
    set.seed(42)
    n <- 15
    p <- 5
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    W <- matrix(rnorm(p * 3), nrow = p, ncol = 3)
    wl <- list(multi_weights = W)
    res <- twasPredict(X, wl)

    expect_equal(ncol(res[["multi_predicted"]]), 3)
    expect_equal(res[["multi_predicted"]], X %*% W)
})

test_that("twasPredict: zero weights give zero predictions", {
    X <- matrix(1:6, nrow = 2, ncol = 3)
    wl <- list(null_weights = matrix(0, nrow = 3, ncol = 1))
    res <- twasPredict(X, wl)
    expect_true(all(res[["null_predicted"]] == 0))
})

# ===========================================================================
#
#  twasWeights  (input validation and basic behavior)
#
# ===========================================================================

test_that("twasWeights: X must be a matrix", {
    d <- makeData()
    expect_error(
        learnTwasWeights(as.data.frame(d$X), d$Y, weightMethods = list()),
        "X.*Must be of type 'matrix'"
    )
})

test_that("twasWeights: Y must be a matrix or vector", {
    d <- makeData()
    # In R, is.vector(list(...)) returns TRUE, so a list passes the initial
    # type check and gets converted via matrix(). The resulting matrix has
    # 1 row which mismatches X's 50 rows, triggering the row count error.
    expect_error(
        learnTwasWeights(d$X, list(d$Y), weightMethods = list()),
        "One of the following must apply"
    )
})

test_that("twasWeights: Y as vector gets converted to matrix internally", {
    d <- makeData()
    y_vec <- as.numeric(d$Y)

    # Mock lassoWeights (an existing package function) to return trivial weights
    local_mocked_bindings(
        lassoWeights = function(X, y, ...) rep(0, ncol(X))
    )
    result <- learnTwasWeights(
        d$X,
        y_vec,
        weightMethods = list(lasso_weights = list())
    )
    expect_true(is(result, "TwasWeights"))
    expect_equal(length(getMethodNames(result)), 1)
    expect_equal(nrow(.weightsByMethod(result, "lassoWeights")), ncol(d$X))
    # Weight vector length must equal number of predictors and be numeric/finite
    w <- .weightsByMethod(result, "lassoWeights")[, 1]
    expect_equal(length(w), ncol(d$X))
    expect_true(is.numeric(w))
    expect_true(all(is.finite(w)))
})

test_that("twasWeights: mismatched row counts error", {
    d <- makeData(n = 50, p = 10)
    Y_short <- d$Y[1:30, , drop = FALSE]
    expect_error(
        learnTwasWeights(d$X, Y_short, weightMethods = list()),
        "Y.*Must have exactly 50 rows"
    )
})

test_that("twasWeights: character weight_methods input is accepted", {
    d <- makeData()
    local_mocked_bindings(
        lassoWeights = function(X, y, ...) rep(0, ncol(X))
    )
    # Short name should be resolved via .twas_method_lookup
    result <- learnTwasWeights(d$X, d$Y, weightMethods = c("lasso"))
    expect_true(is(result, "TwasWeights"))
    expect_equal(getMethodNames(result), "lasso")
})

test_that("twasWeights: zero variance columns are filtered and padded back with zeros", {
    d <- makeData(n = 50, p = 10, add_zero_var_col = TRUE)
    p_with_extra <- ncol(d$X) # 11 columns, last is zero-var

    local_mocked_bindings(
        lassoWeights = function(X, y, ...) {
            # After filtering, the zero-var column should be removed
            # So ncol(X) should be p (10), not p+1 (11)
            rep(1, ncol(X))
        }
    )
    result <- learnTwasWeights(
        d$X,
        d$Y,
        weightMethods = list(lasso_weights = list())
    )

    # The returned weight matrix should have rows equal to total columns (including zero-var)
    expect_equal(nrow(.weightsByMethod(result, "lassoWeights")), p_with_extra)
    # The zero-var column weight should be 0 (padded back)
    expect_equal(
        unname(.weightsByMethod(result, "lassoWeights")["chr1:9900:A:G", 1]),
        0
    )
})

test_that("twasWeights: rownames of result match colnames of X", {
    d <- makeData()
    local_mocked_bindings(
        enetWeights = function(X, y, ...) rep(0.1, ncol(X))
    )
    result <- learnTwasWeights(
        d$X,
        d$Y,
        weightMethods = list(enetWeights = list())
    )
    expect_equal(
        rownames(.weightsByMethod(result, "enetWeights")),
        colnames(d$X)
    )
})

test_that("twasWeights: result dimensions match ncol(X) x ncol(Y)", {
    d <- makeData()
    local_mocked_bindings(
        enetWeights = function(X, y, ...) rep(0, ncol(X))
    )
    result <- learnTwasWeights(
        d$X,
        d$Y,
        weightMethods = list(enetWeights = list())
    )
    expect_equal(
        dim(.weightsByMethod(result, "enetWeights")),
        c(ncol(d$X), ncol(d$Y))
    )
})

test_that("twasWeights: multiple methods return named list with one entry per method", {
    d <- makeData()
    local_mocked_bindings(
        lassoWeights = function(X, y, ...) rep(0.1, ncol(X)),
        enetWeights = function(X, y, ...) rep(0.2, ncol(X))
    )
    result <- learnTwasWeights(
        d$X,
        d$Y,
        weightMethods = list(lasso_weights = list(), enetWeights = list())
    )
    expect_equal(length(getMethodNames(result)), 2)
    expect_true("lasso" %in% getMethodNames(result))
    expect_true("enet" %in% getMethodNames(result))
})

# ===========================================================================
#
#  twasWeights with actual glmnet (lasso/enet) -- skip if not available
#
# ===========================================================================

test_that("twasWeights: lassoWeights produces correct structure with real glmnet", {
    skip_if_not_installed("glmnet")
    d <- makeData(n = 50, p = 10)
    result <- learnTwasWeights(
        d$X,
        d$Y,
        weightMethods = list(lasso_weights = list())
    )

    expect_true(is(result, "TwasWeights"))
    expect_equal(getMethodNames(result), "lasso")
    expect_equal(nrow(.weightsByMethod(result, "lassoWeights")), ncol(d$X))
    expect_equal(ncol(.weightsByMethod(result, "lassoWeights")), 1)
    # At least some weights should be non-zero for this strong signal
    expect_true(any(.weightsByMethod(result, "lassoWeights") != 0))
})

test_that("twasWeights: enetWeights produces correct structure with real glmnet", {
    skip_if_not_installed("glmnet")
    d <- makeData(n = 50, p = 10)
    result <- learnTwasWeights(
        d$X,
        d$Y,
        weightMethods = list(enetWeights = list())
    )

    expect_true(is(result, "TwasWeights"))
    expect_equal(getMethodNames(result), "enet")
    expect_equal(nrow(.weightsByMethod(result, "enetWeights")), ncol(d$X))
})

# ===========================================================================
#
#  twasWeightsCv  (input validation)
#
# ===========================================================================

test_that("twasWeightsCv: NULL weight_methods returns only samplePartition", {
    d <- makeData()
    result <- twasWeightsCv(d$X, d$Y, fold = 3, weightMethods = NULL)
    expect_equal(names(result), "samplePartition")
    expect_true(is.data.frame(result$samplePartition))
})

test_that("twasWeightsCv: samplePartition structure is correct", {
    d <- makeData()
    result <- twasWeightsCv(d$X, d$Y, fold = 5, weightMethods = NULL)
    sp <- result$samplePartition

    expect_true("Sample" %in% colnames(sp))
    expect_true("Fold" %in% colnames(sp))
    expect_equal(nrow(sp), nrow(d$X))
    expect_equal(length(unique(sp$Fold)), 5)
    # All sample names should appear
    expect_true(all(rownames(d$X) %in% sp$Sample))
})

test_that("twasWeightsCv: character weight_methods are accepted", {
    d <- makeData()
    local_mocked_bindings(
        lassoWeights = function(X, y, ...) rep(0, ncol(X))
    )
    set.seed(42)
    result <- twasWeightsCv(
        d$X,
        d$Y,
        fold = 2,
        weightMethods = c("lasso")
    )
    expect_true(is.list(result))
    expect_true("prediction" %in% names(result))
})

# ---------------------------------------------------------------------------
# CV with real lassoWeights (integration test)
# ---------------------------------------------------------------------------

test_that("twasWeightsCv: basic CV with lassoWeights produces correct metrics structure", {
    skip_if_not_installed("glmnet")
    d <- makeData(n = 50, p = 10)

    set.seed(123)
    result <- twasWeightsCv(
        d$X,
        d$Y,
        fold = 3,
        weightMethods = list(lasso_weights = list())
    )

    # Structure checks
    expect_true(is.list(result))
    expect_true("samplePartition" %in% names(result))
    expect_true("prediction" %in% names(result))
    expect_true("performance" %in% names(result))
    expect_true("timeElapsed" %in% names(result))

    # Prediction name transformation
    expect_equal(names(result$prediction), "lasso_predicted")

    # Prediction dimensions should match Y
    pred <- result$prediction[["lasso_predicted"]]
    expect_equal(dim(pred), dim(d$Y))

    # Performance table structure
    perf <- result$performance[["lasso_performance"]]
    expect_true(is.matrix(perf))
    expect_equal(
        colnames(perf),
        c("corr", "rsq", "adj_rsq", "pval", "RMSE", "MAE")
    )
    expect_equal(nrow(perf), ncol(d$Y))

    # With strong signal, correlation should be positive
    expect_true(perf[1, "corr"] > 0)
})

test_that("twasWeightsCv: multiple real methods produce per-method metrics", {
    skip_if_not_installed("glmnet")
    d <- makeData(n = 50, p = 10)

    set.seed(99)
    result <- twasWeightsCv(
        d$X,
        d$Y,
        fold = 3,
        weightMethods = list(
            lassoWeights = list(),
            enetWeights = list()
        )
    )

    expect_equal(length(result$prediction), 2)
    expect_equal(length(result$performance), 2)
    expect_true("lasso_predicted" %in% names(result$prediction))
    expect_true("enet_predicted" %in% names(result$prediction))
    expect_true("lasso_performance" %in% names(result$performance))
    expect_true("enet_performance" %in% names(result$performance))
})

# ===========================================================================
#
#  twasWeightsCv: multivariate Y
#
# ===========================================================================

test_that("twasWeightsCv: multivariate Y with multiple columns", {
    d <- makeData(n = 50, p = 10)
    # Create multi-column Y
    set.seed(42)
    Y_multi <- cbind(
        d$Y,
        d$X %*% c(0, 0, 0, 0, 0, 1, -1, 0, 0, 0) + rnorm(50, sd = 0.5)
    )
    colnames(Y_multi) <- c("outcome_1", "outcome_2")
    rownames(Y_multi) <- rownames(d$X)

    local_mocked_bindings(
        lassoWeights = function(X, y, ...) rep(0.1, ncol(X))
    )
    set.seed(42)
    result <- twasWeightsCv(
        d$X,
        Y_multi,
        fold = 2,
        weightMethods = list(lasso_weights = list())
    )

    pred <- result$prediction[["lasso_predicted"]]
    expect_equal(ncol(pred), 2)
    expect_equal(nrow(pred), 50)

    perf <- result$performance[["lasso_performance"]]
    expect_equal(nrow(perf), 2)
    expect_equal(rownames(perf), c("outcome_1", "outcome_2"))
})

# ===========================================================================
#
#  twasWeightsPipeline  (structure and input validation)
#
# ===========================================================================

test_that("learnTwasWeights refuses the susie + susieInf pair without fits", {
    # The chained susieInf -> susie fit lives in fineMappingPipeline() now;
    # learnTwasWeights never fine-maps, so the pair is an error here.
    # fitSusieInfThenSusie() itself is covered in test_fineMappingWrappers.R.
    d <- makeData(n = 50, p = 10)
    expect_error(
        learnTwasWeights(
            d$X,
            as.numeric(d$Y),
            weightMethods = list(
                susie_weights = list(),
                susie_inf_weights = list()
            )
        ),
        "susie, susieInf"
    )
})


test_that("learnTwasWeights resolves fits under camelCase method names", {
    # `susieWeights` and `susie_weights` name the same method; a fit supplied
    # for one spelling must land on the other's arguments too.
    d <- makeData(n = 50, p = 10)
    seen <- NULL
    local_mocked_bindings(
        susieWeights = function(X, y, susieFit = NULL, ...) {
            seen <<- susieFit
            rep(0, ncol(X))
        }
    )
    learnTwasWeights(
        d$X,
        as.numeric(d$Y),
        weightMethods = list(susie_weights = list()),
        fittedModels = list(susie = makeFakeSusieFit(p = 10, L = 5))
    )
    expect_true("susie" %in% class(seen))
})


test_that("learnTwasWeights runs susie + susieInf from supplied fits", {
    d <- makeData(n = 50, p = 10)
    local_mocked_bindings(
        susieInfWeights = function(X, y, ...) rep(0, ncol(X)),
        susieWeights = function(X, y, ...) rep(0, ncol(X))
    )
    result <- learnTwasWeights(
        d$X,
        as.numeric(d$Y),
        weightMethods = list(
            susie_weights = list(),
            susie_inf_weights = list()
        ),
        fittedModels = list(
            susie = makeFakeSusieFit(p = 10, L = 5),
            susieInf = makeFakeSusieFit(p = 10, L = 7, inf = TRUE)
        )
    )
    expect_equal(getMethodNames(result), c("susie", "susie_inf"))
})


test_that("twasWeightsCv: NA values in Y trigger NA-removal branch in metrics", {
    set.seed(42)
    n <- 30
    p <- 5
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * seq_len(p))
    rownames(X) <- paste0("s", seq_len(n))
    Y <- matrix(rnorm(n), ncol = 1)
    rownames(Y) <- rownames(X)
    colnames(Y) <- "outcome"
    Y[c(3, 11, 17), 1] <- NA # introduce NAs

    # Mock to return non-zero (so prediction has nonzero variance and lm_fit runs)
    local_mocked_bindings(
        lassoWeights = function(X, y, ...) {
            w <- rep(0, ncol(X))
            w[1] <- 0.5
            w
        }
    )
    set.seed(42)
    result <- twasWeightsCv(
        X,
        Y,
        fold = 2,
        weightMethods = list(lasso_weights = list())
    )
    perf <- result$performance[["lasso_performance"]]
    # NA-removal branch ran; metrics should be finite (not all-NA)
    expect_true(is.finite(perf[1, "rsq"]))
})

test_that("twasWeightsCv: dataDrivenPriorMatricesCv is plumbed through", {
    set.seed(42)
    n <- 20
    p <- 4
    X <- matrix(rnorm(n * p), nrow = n)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * seq_len(p))
    rownames(X) <- paste0("s", seq_len(n))
    Y <- matrix(rnorm(n * 2), nrow = n)
    colnames(Y) <- c("y1", "y2")
    rownames(Y) <- rownames(X)

    captured_args <- list()
    local_mocked_bindings(
        mrmashWeights = function(X, Y, ...) {
            captured_args[[length(captured_args) + 1]] <<- list(...)
            matrix(
                0,
                nrow = ncol(X),
                ncol = ncol(Y),
                dimnames = list(colnames(X), colnames(Y))
            )
        }
    )
    prior_cv <- list(matrix(1, 2, 2), matrix(2, 2, 2))
    set.seed(42)
    result <- twasWeightsCv(
        X,
        Y,
        fold = 2,
        weightMethods = list(mrmashWeights = list()),
        dataDrivenPriorMatricesCv = prior_cv
    )
    # mrmashWeights mock should have been called and received the per-fold prior
    # matrix under the camelCase name that actually binds mrmashWrapper's
    # `dataDrivenPriorMatrices` argument (the snake_case form was a latent no-op).
    expect_true(length(captured_args) >= 1)
    expect_true(any(map_lgl(
        captured_args,
        function(a) {
            "dataDrivenPriorMatrices" %in% names(a)
        }
    )))
})

# ===========================================================================
# twasWeightsPipeline: removed_methods warning + max_cv_variants subsampling
# ===========================================================================

# ===========================================================================
# twasWeights: dim-fix branch when nrow(weights_matrix) != length(valid_columns)
# ===========================================================================

test_that("twasWeights: multivariate weights_matrix is reduced to valid_columns when row counts mismatch", {
    set.seed(42)
    n <- 20
    p <- 5
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    # all columns valid (no zero variance)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * seq_len(p))
    Y <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
    colnames(Y) <- c("y1", "y2")

    local_mocked_bindings(
        mrmashWeights = function(X, Y, ...) {
            # Return more rows than valid_columns (length p) so the dim-fix branch
            # subsets the matrix back to names(valid_columns).
            extra_rows <- p + 2
            m <- matrix(
                seq_len(extra_rows * ncol(Y)),
                nrow = extra_rows,
                ncol = ncol(Y)
            )
            rownames(m) <- c(
                sprintf("chr1:%d:A:G", 100L * seq_len(p)),
                "extra1",
                "extra2"
            )
            colnames(m) <- colnames(Y)
            m
        }
    )
    result <- learnTwasWeights(
        X,
        Y,
        weightMethods = list(mrmashWeights = list())
    )
    # After the dim-fix, the weights matrix is restricted to v1..v5 -> shape p x ncol(Y)
    expect_equal(nrow(.weightsByMethod(result, "mrmashWeights")), p)
    expect_equal(ncol(.weightsByMethod(result, "mrmashWeights")), 2)
    expect_equal(
        rownames(.weightsByMethod(result, "mrmashWeights")),
        sprintf("chr1:%d:A:G", 100L * seq_len(p))
    )
})

# ===========================================================================
# Tests migrated from test_twas.R (twasWeightsCv, learnTwasWeights, twasPredict)
# ===========================================================================

test_that("twasWeightsCv is reproducible with seed", {
    sim <- generateXY(seed = 1)
    X <- sim$X
    y <- sim$Y
    local_mocked_bindings(
        enetWeights = function(X, y, ...) rnorm(ncol(X)),
        glmnetWeights = function(X, y, ...) runif(ncol(X))
    )
    # Non-SuSiE methods: this test is about the seeded fold partition, and a
    # SuSiE-family token now requires per-fold fits it has no reason to carry.
    weight_methods_test <- list(enetWeights = list(), glmnetWeights = list())
    set.seed(1)
    result_seed1 <- twasWeightsCv(
        X,
        y,
        fold = 2,
        weightMethods = weight_methods_test
    )
    set.seed(1)
    result_seed2 <- twasWeightsCv(
        X,
        y,
        fold = 2,
        weightMethods = weight_methods_test
    )
    expect_equal(result_seed1$samplePartition, result_seed2$samplePartition)
})


test_that("twasWeightsCv handles errors appropriately", {
    sim <- generateXY(seed = 1)
    X <- sim$X
    y <- sim$Y
    local_mocked_bindings(
        susieWeights = function(X, y, ...) rnorm(ncol(X)),
        glmnetWeights = function(X, y, ...) runif(ncol(X))
    )
    weight_methods_test <- list(susie = list(), glmnetWeights = list())
    expect_error(twasWeightsCv(X, y, fold = NULL), "fold.*samplePartitions")
    expect_error(
        twasWeightsCv(X, y, fold = "invalid"),
        "Must be of type 'count'"
    )
    expect_error(twasWeightsCv(X, y, fold = -1), "Must be >= 1")
    expect_error(twasWeightsCv(2, y, fold = 2), "Must be of type 'matrix'")
    expect_error(twasWeightsCv(X, 2, fold = 2), "Y.*Must have exactly 10 rows")
    expect_error(
        twasWeightsCv(
            matrix(rnorm(4, nrow = 2)),
            matrix(rnorm(2, nrow = 1)),
            fold = 2
        ),
        "unused argument"
    )
    expect_error(twasWeightsCv(X, y), "fold.*samplePartitions")
})


test_that("learnTwasWeights handles errors appropriately", {
    sim <- generateXY(seed = 1)
    X <- sim$X
    y <- sim$Y
    local_mocked_bindings(
        susieWeights = function(X, y, ...) rnorm(ncol(X)),
        glmnetWeights = function(X, y, ...) runif(ncol(X))
    )
    weight_methods_test <- list(susie = list(), glmnetWeights = list())
    expect_error(
        learnTwasWeights(
            matrix(rnorm(4, nrow = 2)),
            matrix(rnorm(2, nrow = 1))
        ),
        "weightMethods.*is missing"
    )
    expect_error(learnTwasWeights(X, y), "weightMethods")
})

# ===========================================================================
# twasZ: mathematical correctness (single-method / vector path)
# ===========================================================================

test_that("twasPredict multiplies X by weights", {
    X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
    weights_list <- list(method1_weights = c(0.5, -0.5))
    result <- twasPredict(X, weights_list)
    expect_length(result, 1)
    expect_equal(names(result), "method1_predicted")
    expected <- X %*% c(0.5, -0.5)
    expect_equal(result[[1]], expected)
})


test_that("twasPredict handles multiple weight methods", {
    set.seed(42)
    X <- matrix(rnorm(30), nrow = 10, ncol = 3)
    weights_list <- list(
        lassoWeights = c(1, 0, -1),
        enetWeights = c(0.5, 0.3, 0.2),
        susieWeights = c(0, 0, 1)
    )
    result <- twasPredict(X, weights_list)
    expect_length(result, 3)
    expect_equal(
        names(result),
        c("lassoPredicted", "enetPredicted", "susiePredicted")
    )
    # Verify computation for one method
    expect_equal(result$lassoPredicted, X %*% c(1, 0, -1))
})


test_that("twasPredict with zero weights gives zero predictions", {
    X <- matrix(1:6, nrow = 2, ncol = 3)
    weights_list <- list(null_weights = c(0, 0, 0))
    result <- twasPredict(X, weights_list)
    expect_true(all(result$null_predicted == 0))
})


test_that("twasPredict with single variant", {
    X <- matrix(c(1, 2, 3), nrow = 3, ncol = 1)
    weights_list <- list(single_weights = 2.0)
    result <- twasPredict(X, weights_list)
    expect_equal(as.numeric(result$single_predicted), c(2, 4, 6))
})


# === Tests migrated from test_s4Constructors.R (TwasWeights) ===

test_that("TwasWeights: builds a collection keyed by 4-tuple", {
    e1 <- .sc_makeTwasWeightsRow()
    e2 <- .sc_makeTwasWeightsRow()
    tw <- TwasWeights(
        study = c("s1", "s1"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        method = c("lasso", "enet"),
        entry = list(e1, e2)
    )
    expect_s4_class(tw, "TwasWeights")
    expect_equal(nrow(tw), 2L)
    expect_setequal(getMethodNames(tw), c("lasso", "enet"))
})


test_that("TwasWeights: getStudy / getContexts / getTraits / getMethodNames", {
    e <- .sc_makeTwasWeightsRow()
    tw <- TwasWeights(
        study = c("s1", "s2"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        method = c("lasso", "lasso"),
        entry = list(e, e)
    )
    expect_setequal(getContexts(tw), c("c1", "c2"))
    expect_equal(getTraits(tw), "t1")
    expect_equal(getMethodNames(tw), "lasso")
})


test_that("TwasWeights: rejects duplicate 4-tuples", {
    e <- .sc_makeTwasWeightsRow()
    expect_error(
        TwasWeights(
            study = c("s1", "s1"),
            context = c("c1", "c1"),
            trait = c("t1", "t1"),
            method = c("lasso", "lasso"),
            entry = list(e, e)
        ),
        "uniqueness violated"
    )
})


test_that("TwasWeights: joint columns work the same as on the FMR class", {
    e <- .sc_makeTwasWeightsRow()
    # Univariate lasso at c1 + the c1 slice of an mr.mash joint over (c1, c2).
    tw <- TwasWeights(
        study = c("s1", "s1"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        method = c("lasso", "mrmash"),
        entry = list(e, e),
        jointContexts = c(NA_character_, "c1;c2")
    )
    expect_true("jointContexts" %in% pecotmr:::.tupleColumnNames(tw))
    expect_identical(tw$jointContexts, c(NA_character_, "c1;c2"))
    # uniqueness: same (s1, c1, t1, mrmash) tuple from two joint fits over
    # (c1, c2) and (c1, c3) -> distinct rows via jointContexts.
    tw2 <- TwasWeights(
        study = c("s1", "s1"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        method = c("mrmash", "mrmash"),
        entry = list(e, e),
        jointContexts = c("c1;c2", "c1;c3")
    )
    expect_equal(nrow(tw2), 2L)
})


test_that("TwasWeights: getTwasWeights extracts the entry for a tuple", {
    e1 <- .sc_makeTwasWeightsRow()
    e2 <- .sc_makeTwasWeightsRow()
    tw <- TwasWeights(
        study = c("s1", "s1"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        method = c("lasso", "enet"),
        entry = list(e1, e2)
    )
    # getTwasWeights() returns the single-row COLLECTION for that tuple, not a
    # detached entry object, so compare what the row carries.
    picked <- getTwasWeights(
        tw,
        study = "s1",
        context = "c1",
        trait = "t1",
        method = "enet"
    )
    expect_s4_class(picked, "TwasWeights")
    expect_equal(nrow(picked), 1L)
    expect_equal(as.character(picked$method), "enet")
    expect_identical(getVariantIds(picked), .twrPartsVariantIds(e2))
    expect_equal(unname(getWeights(picked)), unname(getWeights(e2)))
})

# ===========================================================================
# LdData
# ===========================================================================

# === Tests migrated from test_showMethods.R (TwasWeights) ===

test_that("show.TwasWeights prints entry/study/context/trait/method counts", {
    e <- .sh_makeTwEntry()
    tw <- TwasWeights(
        study = c("s1", "s1"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        method = c("lasso", "enet"),
        entry = list(e, e)
    )
    out <- capture.output(show(tw))
    expect_true(any(grepl("TwasWeights: 2 entries", out)))
    expect_true(any(grepl("1 studies.*2 contexts.*1 traits.*2 methods", out)))
})


test_that("show.TwasWeights reports ldSketch when present", {
    e <- .sh_makeTwEntry()
    tw <- TwasWeights(
        study = "s1",
        context = "c1",
        trait = "t1",
        method = "lasso",
        entry = list(e),
        ldSketch = .sh_makeGenotypeHandle()
    )
    out <- capture.output(show(tw))
    expect_true(any(grepl("LD sketch: gds @ /tmp/test.gds", out)))
})

# ===========================================================================
#
#  .normalizeCvFolds: fold-spec normalization (direct internal calls)
#
# ===========================================================================

test_that(".normalizeCvFolds: list cvFolds + samplePartition are mutually exclusive", {
    expect_error(
        pecotmr:::.normalizeCvFolds(
            cvFolds = list(c("s1", "s2"), c("s3", "s4")),
            samplePartition = data.frame(Sample = "s1", Fold = 1)
        ),
        "not both"
    )
})

test_that(".normalizeCvFolds: samplePartition must have Sample and Fold columns", {
    expect_error(
        pecotmr:::.normalizeCvFolds(samplePartition = data.frame(a = 1, b = 2)),
        "must have columns"
    )
})

test_that(".normalizeCvFolds: a sample assigned to >1 fold errors", {
    sp <- data.frame(
        Sample = c("s1", "s1", "s2"),
        Fold = c(1, 2, 2),
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::.normalizeCvFolds(samplePartition = sp),
        "more than one fold"
    )
})

test_that(".normalizeCvFolds: unknown sample (vs sampleNames) errors", {
    sp <- data.frame(
        Sample = c("s1", "s99"),
        Fold = c(1, 2),
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::.normalizeCvFolds(
            samplePartition = sp,
            sampleNames = c("s1", "s2")
        ),
        "unknown sample"
    )
})

test_that(".normalizeCvFolds: uncovered samples error", {
    sp <- data.frame(Sample = "s1", Fold = 1, stringsAsFactors = FALSE)
    expect_error(
        pecotmr:::.normalizeCvFolds(
            samplePartition = sp,
            sampleNames = c("s1", "s2", "s3")
        ),
        "does not cover"
    )
})

test_that(".normalizeCvFolds: valid samplePartition returns df + nFolds", {
    sp <- data.frame(
        Sample = c("s1", "s2", "s3", "s4"),
        Fold = c(1, 1, 2, 2),
        stringsAsFactors = FALSE
    )
    res <- pecotmr:::.normalizeCvFolds(
        samplePartition = sp,
        sampleNames = c("s1", "s2", "s3", "s4")
    )
    expect_equal(res$nFolds, 2L)
    expect_equal(nrow(res$samplePartition), 4L)
    expect_setequal(unique(res$samplePartition$Fold), c(1, 2))
})

test_that(".normalizeCvFolds: list-form requires at least 2 folds", {
    expect_error(
        pecotmr:::.normalizeCvFolds(cvFolds = list(c("s1", "s2"))),
        "at least 2 folds"
    )
})

test_that(".normalizeCvFolds: numeric fold ids require sampleNames", {
    expect_error(
        pecotmr:::.normalizeCvFolds(cvFolds = list(c(1, 2), c(3, 4))),
        "Numeric fold vectors require"
    )
})

test_that(".normalizeCvFolds: numeric fold ids out of range error", {
    expect_error(
        pecotmr:::.normalizeCvFolds(
            cvFolds = list(c(1, 2), c(3, 99)),
            sampleNames = c("s1", "s2", "s3", "s4")
        ),
        "out-of-range"
    )
})

test_that(".normalizeCvFolds: numeric fold ids resolve via sampleNames", {
    res <- pecotmr:::.normalizeCvFolds(
        cvFolds = list(c(1, 2), c(3, 4)),
        sampleNames = c("s1", "s2", "s3", "s4")
    )
    expect_equal(res$nFolds, 2L)
    expect_setequal(res$samplePartition$Sample, c("s1", "s2", "s3", "s4"))
    expect_equal(
        res$samplePartition$Fold[res$samplePartition$Sample == "s1"],
        1
    )
})

test_that(".normalizeCvFolds: character fold-id vectors are used as-is", {
    res <- pecotmr:::.normalizeCvFolds(
        cvFolds = list(c("s1", "s2"), c("s3", "s4")),
        sampleNames = c("s1", "s2", "s3", "s4")
    )
    expect_equal(res$nFolds, 2L)
    expect_equal(nrow(res$samplePartition), 4L)
})

test_that(".normalizeCvFolds: integer cvFolds returns NULL partition + nFolds = k", {
    res <- pecotmr:::.normalizeCvFolds(cvFolds = 5)
    expect_null(res$samplePartition)
    expect_equal(res$nFolds, 5L)
})

test_that(".normalizeCvFolds: a non-integer scalar cvFolds errors", {
    expect_error(
        pecotmr:::.normalizeCvFolds(cvFolds = "not_a_fold_spec"),
        "single integer, a list of fold vectors"
    )
})

# ===========================================================================
#
#  TwasWeights constructor: length-mismatch guards + accessors
#
# ===========================================================================

test_that("TwasWeights: mismatched core-vector lengths error", {
    e <- .sc_makeTwasWeightsRow()
    expect_error(
        TwasWeights(
            study = c("s1", "s2"),
            context = "c1",
            trait = "t1",
            method = "lasso",
            entry = list(e)
        ),
        "same length"
    )
})

test_that("TwasWeights: joint* column length must match study", {
    e <- .sc_makeTwasWeightsRow()
    expect_error(
        TwasWeights(
            study = "s1",
            context = "c1",
            trait = "t1",
            method = "lasso",
            entry = list(e),
            jointStudies = c("a", "b")
        ),
        "same length as"
    )
})

test_that("TwasWeights: getStandardized/getDataType/getVariantIds delegate to the entry", {
    e1 <- twasWeightsRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
        weights = rnorm(4),
        standardized = TRUE,
        dataType = "expression"
    )
    e2 <- twasWeightsRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
        weights = rnorm(4),
        standardized = FALSE,
        dataType = "splicing"
    )
    tw <- TwasWeights(
        study = c("s1", "s1"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        method = c("lasso", "enet"),
        entry = list(e1, e2)
    )

    expect_true(getStandardized(
        tw,
        study = "s1",
        context = "c1",
        trait = "t1",
        method = "lasso"
    ))
    expect_false(getStandardized(
        tw,
        study = "s1",
        context = "c1",
        trait = "t1",
        method = "enet"
    ))
    expect_equal(
        getDataType(
            tw,
            study = "s1",
            context = "c1",
            trait = "t1",
            method = "enet"
        ),
        "splicing"
    )
    expect_equal(
        getVariantIds(
            tw,
            study = "s1",
            context = "c1",
            trait = "t1",
            method = "lasso"
        ),
        sprintf("chr1:%d:A:G", 100L * (1:4))
    )
})

test_that("TwasWeights: getStudy returns unique study labels", {
    e <- .sc_makeTwasWeightsRow()
    tw <- TwasWeights(
        study = c("s1", "s2"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        method = c("lasso", "lasso"),
        entry = list(e, e)
    )
    expect_setequal(getStudy(tw), c("s1", "s2"))
})

test_that("TwasWeights: getWeights/getCvResult/getFits/getLdSketch delegate per tuple", {
    w1 <- rnorm(4)
    e1 <- twasWeightsRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
        weights = w1,
        fits = list(tag = "fitA"),
        cvResult = list(rsq = 0.42)
    )
    e2 <- twasWeightsRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
        weights = rnorm(4)
    )
    tw <- TwasWeights(
        study = c("s1", "s1"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        method = c("lasso", "enet"),
        entry = list(e1, e2),
        ldSketch = .sh_makeGenotypeHandle()
    )

    expect_equal(
        getWeights(
            tw,
            study = "s1",
            context = "c1",
            trait = "t1",
            method = "lasso"
        ),
        w1
    )
    expect_equal(
        getCvResult(
            tw,
            study = "s1",
            context = "c1",
            trait = "t1",
            method = "lasso"
        )$rsq,
        0.42
    )
    expect_equal(
        getFits(
            tw,
            study = "s1",
            context = "c1",
            trait = "t1",
            method = "lasso"
        )$tag,
        "fitA"
    )
    expect_s4_class(getLdSketch(tw), "RangedSummarizedExperiment")
})

# ===========================================================================
#
#  .resolveMethodFunction: unresolvable-key fallback
#
# ===========================================================================

test_that(".resolveMethodFunction: unresolvable key falls back to the key itself", {
    expect_equal(
        pecotmr:::.resolveMethodFunction("no_such_fn_xyz"),
        "no_such_fn_xyz"
    )
})

# ===========================================================================
#
#  .prepareSusieWeightMethods: seed susie from a supplied susieInf fit
#
# ===========================================================================

test_that(".prepareSusieWeightMethods writes supplied fits onto the method args", {
    infFit <- makeFakeSusieFit(p = 8, L = 3, inf = TRUE)
    susieFit <- makeFakeSusieFit(p = 8, L = 5)

    wm <- pecotmr:::.prepareSusieWeightMethods(
        weightMethods = list(
            susie_weights = list(),
            susie_inf_weights = list()
        ),
        fittedModels = list(susie = susieFit, susieInf = infFit)
    )

    # Each supplied fit is class-tagged and lands on its own method's args.
    # susie's fitting arguments are NOT derived from the inf fit: that
    # prepares a susie fit, which belongs to fineMappingPipeline().
    expect_true("susieInf" %in% class(wm$susie_inf_weights$susieInfFit))
    expect_true("susie" %in% class(wm$susie_weights$susieFit))
    expect_null(wm$susie_weights$model_init)
})

# ===========================================================================
#
#  twasWeightsCv: no-seed warning + multivariate/univariate fitter branches
#
# ===========================================================================

test_that("twasWeightsCv: warns when no random seed has been set", {
    d <- makeData(n = 20, p = 5)
    # makeData() set a seed; drop it just before the call so the unset-seed
    # branch (verbose>=1) is exercised. The fold-sampling at line ~713 restores
    # .Random.seed afterwards, so later tests are unaffected.
    if (exists(".Random.seed", envir = .GlobalEnv)) {
        rm(".Random.seed", envir = .GlobalEnv)
    }
    expect_message(
        twasWeightsCv(d$X, d$Y, fold = 2, weightMethods = NULL),
        "No seed set"
    )
})

test_that("twasWeightsCv: mvsusie per-fold reweighted prior is plumbed (verbose=2)", {
    set.seed(42)
    n <- 24
    p <- 4
    X <- matrix(rnorm(n * p), nrow = n)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * seq_len(p))
    rownames(X) <- paste0("s", seq_len(n))
    Y <- matrix(rnorm(n * 2), nrow = n)
    colnames(Y) <- c("y1", "y2")
    rownames(Y) <- rownames(X)

    captured <- list()
    local_mocked_bindings(
        mvsusieWeights = function(X, Y, ...) {
            captured[[length(captured) + 1]] <<- list(...)
            matrix(
                0,
                nrow = ncol(X),
                ncol = ncol(Y),
                dimnames = list(colnames(X), colnames(Y))
            )
        }
    )
    prior_cv <- list(matrix(1, 2, 2), matrix(2, 2, 2))
    set.seed(1)
    # A SuSiE-family token needs that fold's own fit; the fitter is mocked
    # here, so a stub per fold is enough to reach the per-fold prior path.
    sp <- suppressMessages(twasWeightsCv(X, Y, fold = 2))$samplePartition
    foldFits <- list(fold_1 = "FIT1", fold_2 = "FIT2")
    attr(foldFits, "partitionKey") <- pecotmr:::.cvPartitionKey(sp)
    result <- suppressMessages(twasWeightsCv(
        X,
        Y,
        samplePartitions = sp,
        weightMethods = list(mvsusieWeights = list()),
        reweightedMixturePriorCv = prior_cv,
        fittedModelsCv = list(mvsusie = foldFits),
        verbose = 2
    ))
    expect_true("prediction" %in% names(result))
    # the per-fold prior_variance was forwarded to the multivariate fitter
    expect_true(any(map_lgl(
        captured,
        function(a) {
            "prior_variance" %in% names(a)
        }
    )))
})

test_that("twasWeightsCv forwards fitRetention to a fitter that takes it", {
    set.seed(42)
    n <- 24
    p <- 4
    X <- matrix(rnorm(n * p), nrow = n)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * seq_len(p))
    rownames(X) <- paste0("s", seq_len(n))
    Y <- matrix(rnorm(n * 2), nrow = n)
    colnames(Y) <- c("y1", "y2")
    rownames(Y) <- rownames(X)

    captured <- list()
    local_mocked_bindings(
        mrmashWeights = function(X, Y, fitRetention = "none", ...) {
            captured[[length(captured) + 1]] <<-
                list(fitRetention = fitRetention)
            matrix(
                0,
                nrow = ncol(X),
                ncol = ncol(Y),
                dimnames = list(colnames(X), colnames(Y))
            )
        }
    )
    set.seed(1)
    result <- suppressMessages(twasWeightsCv(
        X,
        Y,
        fold = 2,
        weightMethods = list(mrmashWeights = list()),
        fitRetention = "slim"
    ))
    expect_true("foldFits" %in% names(result))
    expect_true(all(map_lgl(
        captured,
        function(a) identical(a$fitRetention, "slim")
    )))
})

test_that("twasWeightsCv: univariate fitter runs under verbose=2 (no quiet wrapper)", {
    d <- makeData(n = 30, p = 6)
    local_mocked_bindings(
        lassoWeights = function(X, y, ...) {
            w <- rep(0, ncol(X))
            w[1] <- 0.3
            w
        }
    )
    set.seed(1)
    result <- suppressMessages(twasWeightsCv(
        d$X,
        d$Y,
        fold = 2,
        weightMethods = list(lasso_weights = list()),
        verbose = 2
    ))
    expect_true("prediction" %in% names(result))
    expect_equal(nrow(result$prediction[["lasso_predicted"]]), nrow(d$X))
})

# ===========================================================================
#
#  learnTwasWeights: fitRetention plumbing + verbose=2 + parallel
#
# ===========================================================================

test_that("learnTwasWeights: multivariate fitter, fitRetention + verbose=2", {
    set.seed(42)
    n <- 24
    p <- 5
    X <- matrix(rnorm(n * p), nrow = n)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * seq_len(p))
    rownames(X) <- paste0("s", seq_len(n))
    Y <- matrix(rnorm(n * 2), nrow = n)
    colnames(Y) <- c("y1", "y2")

    captured <- list()
    local_mocked_bindings(
        mrmashWeights = function(X, Y, fitRetention = "none", ...) {
            captured[[length(captured) + 1]] <<-
                list(fitRetention = fitRetention)
            matrix(
                0,
                nrow = ncol(X),
                ncol = ncol(Y),
                dimnames = list(colnames(X), colnames(Y))
            )
        }
    )
    result <- suppressMessages(learnTwasWeights(
        X,
        Y,
        weightMethods = list(mrmashWeights = list()),
        fitRetention = "slim",
        verbose = 2
    ))
    expect_true(is(result, "TwasWeights"))
    expect_equal(captured[[1]]$fitRetention, "slim")
})

test_that("a method that takes no fitRetention is left alone", {
    # Retention is passed only to a weight function that declares it; one
    # that does not keeps nothing, and must not be handed the argument.
    d <- makeData(n = 30, p = 6)
    seen <- NULL
    local_mocked_bindings(
        bayesRWeights = function(X, y, ...) {
            seen <<- names(list(...))
            rep(0, ncol(X))
        }
    )
    result <- suppressMessages(learnTwasWeights(
        d$X,
        d$Y,
        weightMethods = list(bayesRWeights = list()),
        fitRetention = "slim"
    ))
    expect_true(is(result, "TwasWeights"))
    expect_false("fitRetention" %in% seen)
})

test_that("learnTwasWeights: univariate fitter runs under verbose=2", {
    d <- makeData(n = 30, p = 6)
    local_mocked_bindings(
        lassoWeights = function(X, y, ...) rep(0.2, ncol(X))
    )
    result <- suppressMessages(learnTwasWeights(
        d$X,
        d$Y,
        weightMethods = list(lasso_weights = list()),
        verbose = 2
    ))
    expect_true(is(result, "TwasWeights"))
    expect_equal(nrow(.weightsByMethod(result, "lassoWeights")), ncol(d$X))
})

test_that("learnTwasWeights: parallel weights path (numThreads = 2)", {
    d <- makeData(n = 30, p = 6)
    local_mocked_bindings(
        lassoWeights = function(X, y, ...) rep(0.1, ncol(X)),
        enetWeights = function(X, y, ...) rep(0.2, ncol(X))
    )
    result <- suppressMessages(learnTwasWeights(
        d$X,
        d$Y,
        weightMethods = list(lasso_weights = list(), enetWeights = list()),
        numThreads = 2
    ))
    expect_true(is(result, "TwasWeights"))
    expect_setequal(getMethodNames(result), c("lasso", "enet"))
})

# ===========================================================================
#
#  twasPredict: TwasWeights S4-collection path
#
# ===========================================================================

test_that("twasPredict: accepts a TwasWeights S4 collection", {
    set.seed(42)
    p <- 5
    w1 <- rnorm(p)
    w2 <- rnorm(p)
    e1 <- twasWeightsRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * seq_len(p)),
        weights = w1
    )
    e2 <- twasWeightsRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * seq_len(p)),
        weights = w2
    )
    tw <- TwasWeights(
        study = c("s1", "s1"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        method = c("lasso", "enet"),
        entry = list(e1, e2)
    )

    X <- matrix(rnorm(8 * p), nrow = 8, ncol = p)
    res <- twasPredict(X, tw)
    expect_equal(names(res), c("lasso_predicted", "enet_predicted"))
    expect_equal(res[["lasso_predicted"]], X %*% matrix(w1, ncol = 1))
    expect_equal(res[["enet_predicted"]], X %*% matrix(w2, ncol = 1))
})


# ===========================================================================
# TwasWeights validity and variant identity
# ===========================================================================

test_that("validity names the missing key and payload columns", {
    data(twasWeightsExample)
    bad <- twasWeightsExample
    S4Vectors::mcols(bad)$method <- NULL
    expect_error(
        methods::validObject(bad),
        "missing elements \\{'method'\\}"
    )
    expect_true(methods::validObject(twasWeightsExample))
    bad2 <- twasWeightsExample
    S4Vectors::mcols(bad2)$cvResult <- NULL
    expect_equal(
        pecotmr:::.twasValidateEntries(bad2),
        "missing entry payload columns: cvResult"
    )
})

test_that("an unnamed genotype matrix is refused rather than given fake ids", {
    # A synthetic "variant_<i>" label would not create identity, only defer
    # the failure to wherever the genomic range is needed.
    X <- matrix(1:4, 2L, 2L)
    expect_error(
        pecotmr:::.twasVariantIds(X),
        "the genotype matrix has no colnames"
    )
    colnames(X) <- c("chr1:1:A:G", "chr1:2:C:T")
    expect_equal(pecotmr:::.twasVariantIds(X), c("chr1:1:A:G", "chr1:2:C:T"))
})

test_that(".twasApplyRownames leaves weights alone when X has no colnames", {
    weightsList <- list(a = matrix(1:4, nrow = 2L))
    noNames <- matrix(0, nrow = 3L, ncol = 2L)
    # Without variant names on X there is nothing to label the rows with.
    expect_identical(
        pecotmr:::.twasApplyRownames(weightsList, noNames),
        weightsList
    )
    named <- noNames
    colnames(named) <- c("v1", "v2")
    out <- pecotmr:::.twasApplyRownames(weightsList, named)
    expect_equal(rownames(out$a), c("v1", "v2"))
})

test_that(".twasBadColMsg names the offending column and its class", {
    expect_equal(
        as.character(pecotmr:::.twasBadColMsg(
            "study",
            data.frame(study = 1:2)
        )),
        "'study' column must be character (got integer)"
    )
})

test_that(".twasMethodRows keeps a per-outcome context vector", {
    vids <- c("chr1:100:A:G", "chr1:200:C:T")
    Y <- matrix(0, nrow = 4L, ncol = 2L, dimnames = list(NULL, c("y1", "y2")))
    wMat <- matrix(
        c(0.1, 0.2, 0.3, 0.4),
        nrow = 2L,
        dimnames = list(vids, c("y1", "y2"))
    )
    mkCtx <- function(contexts) {
        list(
            Y = Y,
            trait = c("t1", "t2"),
            context = contexts,
            study = "s1",
            fitRetention = "none",
            standardized = TRUE,
            dataType = "rnaseq"
        )
    }
    # One context per outcome column: used as-is, not recycled.
    perOutcome <- pecotmr:::.twasMethodRows(
        "lasso_weights",
        wMat,
        vids,
        mkCtx(c("cA", "cB"))
    )
    expect_length(perOutcome, 2L)
    expect_equal(
        map_chr(perOutcome, function(z) z$context),
        c("cA", "cB")
    )
    # A single context is recycled across the outcomes instead.
    recycled <- pecotmr:::.twasMethodRows(
        "lasso_weights",
        wMat,
        vids,
        mkCtx("cOnly")
    )
    expect_equal(
        map_chr(recycled, function(z) z$context),
        c("cOnly", "cOnly")
    )
})

test_that("twasWeightsCv: argument guards fire", {
    d <- generateXY(seed = 1)
    base <- list(X = d$X, Y = d$Y, fold = 2, weightMethods = list())
    expect_error(
        exec(
            twasWeightsCv,
            !!!list_modify(base, !!!list(samplePartitions = 1L))
        ),
        "samplePartitions.*Must be of type 'data.frame'"
    )
    expect_error(
        exec(twasWeightsCv, !!!list_modify(base, !!!list(maxNumVariants = 0))),
        "maxNumVariants.*is not >= 1"
    )
    expect_error(
        exec(twasWeightsCv, !!!list_modify(base, !!!list(numThreads = 1.5))),
        "numThreads.*Must be of type 'single integerish value'"
    )
    expect_error(
        exec(
            twasWeightsCv,
            !!!list_modify(base, !!!list(fitRetention = "sometimes"))
        ),
        "must be one of"
    )
    expect_error(
        exec(twasWeightsCv, !!!list_modify(base, !!!list(seed = "x"))),
        "seed.*Must be of type 'single integerish value'"
    )
    # Inf is the documented "no cap" sentinel and must still be accepted.
    expect_no_error(
        exec(twasWeightsCv, !!!list_modify(base, !!!list(maxNumVariants = Inf)))
    )
})

test_that("learnTwasWeights: argument guards fire", {
    d <- generateXY(seed = 1)
    base <- list(X = d$X, Y = d$Y, weightMethods = list())
    expect_error(
        exec(learnTwasWeights, !!!list_modify(base, !!!list(study = 1L))),
        "study.*Must be of type 'string'"
    )
    # Called directly, not via modifyList(): modifyList() DROPS an element
    # whose value is NULL, so the argument would fall back to its default.
    expect_error(
        learnTwasWeights(
            d$X,
            d$Y,
            weightMethods = list(),
            standardized = NULL
        ),
        "standardized.*Must be of type 'logical flag'"
    )
    expect_error(
        exec(learnTwasWeights, !!!list_modify(base, !!!list(dataType = 1L))),
        "dataType.*Must be of type 'string'"
    )
    expect_error(
        exec(
            learnTwasWeights,
            !!!list_modify(base, !!!list(weightMethods = 1L))
        ),
        "weightMethods.*One of the following must apply"
    )
})

test_that("twasPredict: weightsList must be a list or TwasWeights", {
    expect_error(
        twasPredict(matrix(0, 2, 2), "nope"),
        "weightsList.*One of the following must apply"
    )
})

# ---------------------------------------------------------------------------
# Cross-validating a SuSiE-family method. Those wrappers extract from a
# supplied fit and never fine-map, so CV is only possible when
# fineMappingPipeline's own CV retained each fold's fit.
# ---------------------------------------------------------------------------

.twcv_foldFits <- function(seed = 11, fold = 3) {
    set.seed(seed)
    data(eqtlRegionExample)
    X <- eqtlRegionExample$X[, 1:40]
    y <- eqtlRegionExample$yRes
    Y <- matrix(y, ncol = 1, dimnames = list(rownames(X), "t1"))
    cv <- pecotmr:::.fmWeightsCv(
        X,
        Y,
        tokens = "susie",
        methodArgs = list(),
        fold = fold,
        verbose = 0,
        seed = 1
    )
    list(X = X, Y = Y, cv = cv, slice = pecotmr:::.fmSliceCv(cv, "susie"))
}

test_that("fineMappingPipeline CV retains a lean per-fold fit", {
    skip_if_not_installed("susieR")
    f <- suppressMessages(.twcv_foldFits())
    expect_false(is.null(f$cv$foldFits))
    expect_equal(names(f$cv$foldFits), c("fold_1", "fold_2", "fold_3"))
    # lean: only the fields the weight extractors read
    fit1 <- f$cv$foldFits[["fold_1"]][["susie"]]
    expect_true(all(c("pip", "alpha", "mu") %in% names(fit1)))
    expect_false("lbf_variable" %in% names(fit1))
    # and it slices per method onto the row payload
    expect_equal(names(f$slice$foldFits), c("fold_1", "fold_2", "fold_3"))
})

test_that("twasWeightsCv cannot cross-validate susie without the fold fits", {
    skip_if_not_installed("susieR")
    f <- suppressMessages(.twcv_foldFits())
    expect_error(
        suppressMessages(twasWeightsCv(
            f$X,
            f$Y,
            samplePartitions = f$slice$samplePartition,
            weightMethods = list(susie_weights = list()),
            verbose = 0
        )),
        "never run fine-mapping themselves"
    )
})

test_that("twasWeightsCv cross-validates susie from the retained fold fits", {
    skip_if_not_installed("susieR")
    f <- suppressMessages(.twcv_foldFits())
    out <- suppressMessages(twasWeightsCv(
        f$X,
        f$Y,
        samplePartitions = f$slice$samplePartition,
        weightMethods = list(susie_weights = list()),
        fittedModelsCv = list(susie = f$slice$foldFits),
        verbose = 0
    ))
    expect_true(all(c("prediction", "performance") %in% names(out)))
    expect_false(is.null(out$prediction))
})

test_that(".twasFoldFit injects the fold's fit under the adapter's fit arg", {
    ff <- list(susie = list(fold_1 = "FIT1", fold_2 = "FIT2"))
    a <- pecotmr:::.twasFoldFit(list(), "susie_weights", 2L, ff)
    expect_identical(a$susieFit, "FIT2")
    # a method with no fine-mapping adapter is untouched
    b <- pecotmr:::.twasFoldFit(list(), "lasso_weights", 1L, ff)
    expect_length(b, 0L)
    # and so is the NULL case
    expect_length(pecotmr:::.twasFoldFit(list(), "susie_weights", 1L, NULL), 0L)
})

test_that("twasWeightsCv refuses fold fits from a different partition", {
    skip_if_not_installed("susieR")
    f <- suppressMessages(.twcv_foldFits())
    wm <- list(susie = list())

    # (a) no partition at all: a freshly drawn one would score each fold with
    # a fit that saw its held-out samples.
    expect_error(
        suppressMessages(twasWeightsCv(
            f$X,
            f$Y,
            weightMethods = wm,
            verbose = 0,
            fittedModelsCv = list(susie = f$slice$foldFits)
        )),
        "needs the fold partition"
    )

    # (b) fits with no fingerprint cannot be shown to match.
    unstamped <- f$slice$foldFits
    attr(unstamped, "partitionKey") <- NULL
    expect_error(
        suppressMessages(twasWeightsCv(
            f$X,
            f$Y,
            weightMethods = wm,
            verbose = 0,
            samplePartitions = f$slice$samplePartition,
            fittedModelsCv = list(susie = unstamped)
        )),
        "no partition fingerprint"
    )

    # (c) a genuinely different split is caught by the fingerprint.
    other <- suppressMessages(pecotmr:::.fmWeightsCv(
        f$X,
        f$Y,
        tokens = "susie",
        methodArgs = list(),
        fold = 3,
        verbose = 0,
        seed = 77
    ))
    expect_error(
        suppressMessages(twasWeightsCv(
            f$X,
            f$Y,
            weightMethods = wm,
            verbose = 0,
            samplePartitions = other$samplePartition,
            fittedModelsCv = list(susie = f$slice$foldFits)
        )),
        "trained on a different fold partition"
    )
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

# ---------------------------------------------------------------------------
# Up-front gate: a SuSiE-family method without its fit is refused before any
# fitting work, rather than surfacing from inside the per-fold map.
# ---------------------------------------------------------------------------

test_that("learnTwasWeights refuses a susie token with no fit", {
    skip_if_not_installed("susieR")
    set.seed(11)
    data(eqtlRegionExample)
    X <- eqtlRegionExample$X[, 1:40]
    Y <- matrix(
        eqtlRegionExample$yRes,
        ncol = 1,
        dimnames = list(rownames(X), "t1")
    )
    expect_error(
        suppressMessages(learnTwasWeights(
            X,
            Y,
            weightMethods = list(susie_weights = list()),
            verbose = 0
        )),
        "never run fine-mapping themselves"
    )
    # supplying it through fittedModels satisfies the gate
    fit <- suppressMessages(susieR::susie(X, Y[, 1], L = 5))
    expect_no_error(suppressMessages(learnTwasWeights(
        X,
        Y,
        weightMethods = list(susie_weights = list()),
        fittedModels = list(susie = fit),
        verbose = 0
    )))
    # a method with no fine-mapping adapter is unaffected
    expect_no_error(suppressMessages(learnTwasWeights(
        X,
        Y,
        weightMethods = list(lasso_weights = list()),
        verbose = 0
    )))
})

test_that("twasWeightsCv refuses a susie token with no per-fold fits", {
    skip_if_not_installed("susieR")
    f <- suppressMessages(.twcv_foldFits())
    expect_error(
        suppressMessages(twasWeightsCv(
            f$X,
            f$Y,
            samplePartitions = f$slice$samplePartition,
            weightMethods = list(susie_weights = list()),
            verbose = 0
        )),
        "needs that fold's own fit"
    )
})

test_that(".twasSusieTokensRequested matches both method spellings", {
    expect_equal(
        pecotmr:::.twasSusieTokensRequested(list(susie = list())),
        "susie"
    )
    expect_equal(
        pecotmr:::.twasSusieTokensRequested(c("mvsusieWeights")),
        "mvsusie"
    )
    expect_length(
        pecotmr:::.twasSusieTokensRequested(list(lasso_weights = list())),
        0L
    )
})

test_that("the twas lookup helpers answer NULL when the token is absent", {
    expect_null(pecotmr:::.twasMethodArgsFor(list(), "susie"))
    expect_null(pecotmr:::.twasFoldFitsFor(NULL, "susie"))
})

# --- per-token argument chains ----------------------------------------------
#
# The recorded chains in .twasMethodChains() are the basis for checking a
# method's arguments, so they must not drift silently. These tests re-derive
# them by EXECUTION -- tracing which functions a real fit enters -- which is
# how they were established in the first place; reading the sources hides the
# forwarding behind list_modify()/exec() and imported-without-:: calls.

# Traces the pecotmr functions entered while `expr` runs, restricted to the
# ones a chain could name. Returns their names.
.tw_traceChain <- function(expr) {
    ns <- asNamespace("pecotmr")
    watched <- c(
        "glmnetWeights",
        "ncvregWeights",
        "bglrWeights",
        "bayesAlphabetWeights",
        "dprWeights",
        "mrmashWrapper",
        ".penalizedRssWeights",
        ".rssShrinkGridWeights",
        "buildMrmashPriorMatrices",
        "penalizedRss",
        "lassosumRss",
        "sdpr",
        "prsCs"
    )
    hit <- new.env(parent = emptyenv())
    for (f in watched) {
        ob <- tryCatch(get(f, envir = ns), error = function(e) NULL)
        if (!is.function(ob)) {
            next
        }
        suppressMessages(trace(
            f,
            where = ns,
            print = FALSE,
            tracer = substitute(
                assign(NM, TRUE, envir = HIT),
                list(NM = f, HIT = hit)
            )
        ))
    }
    on.exit(
        for (f in watched) {
            try(suppressMessages(untrace(f, where = ns)), silent = TRUE)
        },
        add = TRUE
    )
    invisible(tryCatch(
        suppressWarnings(suppressMessages(force(expr))),
        error = function(e) NULL
    ))
    # all.names: the internal hops are dot-prefixed and ls() hides those.
    sort(ls(hit, all.names = TRUE))
}

.tw_rssFixture <- function(p = 6L) {
    set.seed(1)
    LD <- diag(p)
    LD[abs(row(LD) - col(LD)) == 1] <- 0.3
    list(
        LD = LD,
        stat = list(bhat = rnorm(p) * 0.1, shat = rep(0.05, p), n = 500L)
    )
}

test_that("the recorded individual-path chains match what a fit enters", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("ncvreg")
    set.seed(1)
    n <- 60L
    p <- 5L
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * seq_len(p))
    y <- as.numeric(X %*% rnorm(p) + rnorm(n))
    recorded <- function(tk) {
        setdiff(
            pecotmr:::.twasMethodChainFor(tk, "QtlDataset"),
            c(paste0(tk, "Weights"), "lassoWeights", "enetWeights")
        )
    }
    expect_true(is_in(
        "glmnetWeights",
        .tw_traceChain(lassoWeights(X, y))
    ))
    expect_true(is_in("glmnetWeights", recorded("lasso")))
    expect_true(is_in(
        "ncvregWeights",
        .tw_traceChain(scadWeights(X, y))
    ))
    expect_true(is_in("ncvregWeights", recorded("scad")))
})

test_that("the recorded RSS chains match what a fit enters", {
    fx <- .tw_rssFixture()
    seen <- .tw_traceChain(scadRssWeights(fx$stat, fx$LD))
    recorded <- pecotmr:::.twasMethodChainFor("scad", "QtlSumStats")
    # Every pecotmr hop the run entered must be in the recorded chain.
    expect_true(all(is_in(seen, recorded)))
    expect_true(is_in("penalizedRss", seen))
    expect_true(is_in(".rssShrinkGridWeights", seen))
})

test_that("every RSS chain yields a checkable argument set", {
    # The point of the chains: each summary-statistics path ends in one of
    # pecotmr's own solvers, which enumerate their formals, so the RSS side
    # is checkable even where the individual side is not.
    for (tk in c("scad", "mcp", "l0learn", "lasso", "prsCs", "dprGibbs")) {
        expect_false(
            is.null(pecotmr:::.twasChainAccepted(tk, "QtlSumStats")),
            label = paste("accepted set for", tk)
        )
    }
})

test_that("a middle-hop argument is accepted and a typo is not", {
    skip_if_not_installed("mr.mashr")
    # canonicalPriorMatrices belongs to buildMrmashPriorMatrices, in the
    # middle of the mrmash chain -- neither the entry wrapper nor the engine.
    acc <- pecotmr:::.twasChainAccepted("mrmash", "QtlDataset")
    expect_true(is_in("canonicalPriorMatrices", acc))
    expect_true(is_in("max_iter", acc))
    expect_false(is_in("zzz", acc))
})

test_that("per-class checking rejects a name only the other path accepts", {
    skip_if_not_installed("ncvreg")
    # `s` is a formal of the RSS shrinkage grid, not of the individual path.
    sumAcc <- pecotmr:::.twasChainAccepted("scad", "QtlSumStats")
    expect_true(is_in("s", sumAcc))
    expect_error(
        pecotmr:::.twasCheckMethodArgsForInput(
            list(scad = list(zzz = 1)),
            "QtlSumStats"
        ),
        "method 'scad': unknown argument\\(s\\) zzz"
    )
    # The individual path reaches ncvreg, which takes `...`, so nothing can
    # be rejected there.
    expect_silent(
        pecotmr:::.twasCheckMethodArgsForInput(
            list(scad = list(zzz = 1)),
            "QtlDataset"
        )
    )
})

test_that("the TWAS Options registry is derived from each implementation", {
    # It used to be 18 hand-written token -> constructor pairs, which had to
    # restate a pairing the implementation already declares through its own
    # `methodArgs` default. The two paths of one method often reach different
    # packages, so a second table is a second thing to keep in step.
    reg <- pecotmr:::.twasMethodCtors()
    expect_setequal(names(reg), names(pecotmr:::.twasMethodCapabilities))
    expect_setequal(names(reg), pecotmr:::.twasConfigurableMethods())
    expect_true(all(map_lgl(reg, is.function)))

    nameOf <- function(f) {
        if (is.null(f)) {
            return(NA_character_)
        }
        hit <- keep(
            grep("Options$", getNamespaceExports("pecotmr"), value = TRUE),
            function(n) identical(get(n, envir = asNamespace("pecotmr")), f)
        )
        if (length(hit) > 0L) hit[[1L]] else "?"
    }
    # Spot-check the pairs the derivation must reproduce.
    expect_equal(nameOf(reg$lasso), "GlmnetOptions")
    expect_equal(nameOf(reg$mcp), "NcvregOptions")
    expect_equal(nameOf(reg$bayesB), "BglrOptions")
    # prsCs has no individual implementation, so the path-blind view falls
    # back to its summary-statistics constructor rather than dropping it.
    expect_equal(nameOf(reg$prsCs), "PrsCsOptions")
})

test_that("each method's Options constructor is resolved per input path", {
    nameOf <- function(f) {
        if (is.null(f)) {
            return(NA_character_)
        }
        hit <- keep(
            grep("Options$", getNamespaceExports("pecotmr"), value = TRUE),
            function(n) identical(get(n, envir = asNamespace("pecotmr")), f)
        )
        if (length(hit) > 0L) hit[[1L]] else "?"
    }
    # The two paths of one method can reach different PACKAGES, which is why
    # one constructor per method was never enough.
    expect_equal(
        nameOf(pecotmr:::.twasMethodCtorFor("lasso", "QtlDataset")),
        "GlmnetOptions"
    )
    expect_equal(
        nameOf(pecotmr:::.twasMethodCtorFor("lasso", "QtlSumStats")),
        "LassosumOptions"
    )
    expect_equal(
        nameOf(pecotmr:::.twasMethodCtorFor("dprGibbs", "QtlSumStats")),
        "SdprOptions"
    )
    expect_equal(
        nameOf(pecotmr:::.twasMethodCtorFor("mrmash", "QtlSumStats")),
        "MrmashRssOptions"
    )
    # A method that does not run on a path has no constructor there.
    expect_null(pecotmr:::.twasMethodCtorFor("enet", "QtlSumStats"))
    expect_null(pecotmr:::.twasMethodCtorFor("prsCs", "QtlDataset"))
    expect_null(pecotmr:::.twasMethodCtorFor("nosuch", "QtlDataset"))
})

# ===========================================================================
# Lookup and per-fold helpers: the absent-input branches
# ===========================================================================

test_that(".twasMethodListArgs answers an empty list for nothing usable", {
    f <- pecotmr:::.twasMethodListArgs
    expect_equal(f(list()), list())
    # Unnamed entries carry no method to key on.
    expect_equal(f(list(1, 2)), list())
    # Named entries are re-keyed with the _weights suffix dropped.
    expect_named(f(list(lasso_weights = list(a = 1))), "lasso")
})

test_that("name lookups answer NA / NULL when the name resolves to nothing", {
    expect_equal(
        pecotmr:::.twasChainHopFormals("noSuchFunctionAnywhere"),
        NA_character_
    )
    expect_null(pecotmr:::.twasImplCtorName("noSuchImplAnywhere"))
    expect_null(pecotmr:::.twasImplCtorName(NULL))
})

test_that(".twasCanonicalShortName maps a function name back, else passes on", {
    f <- pecotmr:::.twasCanonicalShortName
    map <- c(lassoWeights = "lasso")
    expect_equal(f("lassoWeights", map), "lasso")
    expect_equal(f("lasso", map), "lasso")
    expect_equal(f("unknownThing", map), "unknownThing")
})

test_that(".twasFoldFitsFor answers NULL when no spelling is present", {
    f <- pecotmr:::.twasFoldFitsFor
    expect_null(f(NULL, "susie"))
    # A CV payload keyed by something else is not this token's fits.
    expect_null(f(list(somethingElse = list()), "susie"))
})

test_that(".twasFoldFit leaves args alone when there is no fold fit", {
    f <- pecotmr:::.twasFoldFit
    args <- list(X = 1)
    # No CV payload at all.
    expect_identical(f(args, "susie", 1L, NULL), args)
    # A method with no fine-mapping token to look up.
    expect_identical(f(args, "lasso", 1L, list(susie = list())), args)
    # The token is present but carries no folds, so there is nothing to add.
    expect_identical(f(args, "susie", 1L, list(susie = list())), args)
})
