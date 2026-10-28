context("twasWeights")

.sc_makeTwasWeightsRow <- function(
    p = 5L,
    weightStandardized = FALSE,
    weightsDataType = "expression"
) {
    TwasWeightsRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * seq_len(p)),
        weights = rnorm(p),
        weightStandardized = weightStandardized,
        weightsDataType = weightsDataType
    )
}

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

makeFakeSusieFit <- function(p = 10, maxNumSingleEffects = 3, inf = FALSE) {
    fit <- list(
        alpha = matrix(1 / p, nrow = maxNumSingleEffects, ncol = p),
        mu = matrix(0, nrow = maxNumSingleEffects, ncol = p),
        lbf_variable = matrix(0, nrow = maxNumSingleEffects, ncol = p),
        X_column_scale_factors = rep(1, p),
        pip = rep(0.1, p),
        V = rep(0.5, maxNumSingleEffects),
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
        maxNumSingleEffects = L,
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
    entry <- twasWeights(
        tw,
        studyName = as.character(tw$study)[[idx[[1L]]]],
        context = as.character(tw$context)[[idx[[1L]]]],
        trait = as.character(tw$trait)[[idx[[1L]]]],
        method = shortName
    )
    w <- weights(entry)
    vids <- variantIds(entry)
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
#  .twasMethodLookup
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

test_that(".twasMethodLookup: 'default' (10) and 'fastDefault' (8) presets", {
    fastDefault <- c(
        "susie_weights",
        "susie_inf_weights",
        "mrash_weights",
        "enet_weights",
        "lasso_weights",
        "mcp_weights",
        "scad_weights",
        "l0learn_weights"
    )
    presets <- list(
        default = c(fastDefault, "bayes_r_weights", "bayes_c_weights"),
        fastDefault = fastDefault
    )
    for (preset in names(presets)) {
        expect_equal(
            sort(names(pecotmr:::.twasMethodLookup(preset))),
            sort(presets[[preset]]),
            info = preset
        )
    }
})

test_that(".twasMethodLookup: custom vector of short names", {
    result <- pecotmr:::.twasMethodLookup(c("susie", "enet", "dprVb"))
    expect_equal(
        sort(names(result)),
        sort(c("susie_weights", "enet_weights", "dpr_vb_weights"))
    )
})

test_that(".twasMethodLookup: unknown method produces error", {
    expect_error(
        pecotmr:::.twasMethodLookup(c("susie", "nonexistent_method")),
        "unknown method token"
    )
})

test_that(".twasMethodLookup: default args are set for mrash, not susie", {
    result <- pecotmr:::.twasMethodLookup("fastDefault")
    # susie carries NO fitting defaults: susieWeights extracts from a supplied
    # fit and never runs susie, so `refine` / `L` would have nothing to
    # configure.
    expect_length(result$susie_weights, 0L)
    # mrash keeps mr.ash's own prior grid, so only the iteration cap is set.
    expect_named(result$mrash_weights, "max.iter")
    expect_equal(result$mrash_weights$max.iter, 100)
})

test_that(".twasMethodLookup: methods with no special args get empty list", {
    result <- pecotmr:::.twasMethodLookup(c("enet", "lasso"))
    expect_equal(length(result$enet_weights), 0L)
    expect_equal(length(result$lasso_weights), 0L)
})

test_that(".twasMethodLookup: all DPR variants can coexist", {
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

test_that("twasPredict: names without _weights suffix are kept unchanged", {
    wl <- list(custom_method = matrix(1, nrow = 2, ncol = 1))
    X <- matrix(1:4, nrow = 2, ncol = 2)
    res <- twasPredict(X, wl)

    # gsub("_weights", "_predicted", "custom_method") == "custom_method"
    expect_equal(names(res), "custom_method")
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

# ===========================================================================
#
#  learnTwasWeights  (input validation and basic behavior)
#
# ===========================================================================

test_that("learnTwasWeights: X must be a matrix", {
    d <- makeData()
    expect_error(
        learnTwasWeights(as.data.frame(d$X), d$Y, weightMethods = list()),
        "X.*Must be of type 'matrix'"
    )
})

test_that("learnTwasWeights: Y must be a matrix or vector", {
    d <- makeData()
    # In R, is.vector(list(...)) returns TRUE, so a list passes the initial
    # type check and gets converted via matrix(). The resulting matrix has
    # 1 row which mismatches X's 50 rows, triggering the row count error.
    expect_error(
        learnTwasWeights(d$X, list(d$Y), weightMethods = list()),
        "One of the following must apply"
    )
})

test_that("learnTwasWeights: Y as vector gets converted to matrix internally", {
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
    expect_equal(length(methodNames(result)), 1)
    expect_equal(nrow(.weightsByMethod(result, "lassoWeights")), ncol(d$X))
    # Weight vector length must equal number of predictors and be numeric/finite
    w <- .weightsByMethod(result, "lassoWeights")[, 1]
    expect_equal(length(w), ncol(d$X))
    expect_true(is.numeric(w))
    expect_true(all(is.finite(w)))
})

test_that("learnTwasWeights: mismatched row counts error", {
    d <- makeData(n = 50, p = 10)
    Y_short <- d$Y[1:30, , drop = FALSE]
    expect_error(
        learnTwasWeights(d$X, Y_short, weightMethods = list()),
        "Y.*Must have exactly 50 rows"
    )
})

test_that("learnTwasWeights: character weightMethods input is accepted", {
    d <- makeData()
    local_mocked_bindings(
        lassoWeights = function(X, y, ...) rep(0, ncol(X))
    )
    # Short name should be resolved via .twasMethodLookup
    result <- learnTwasWeights(d$X, d$Y, weightMethods = c("lasso"))
    expect_true(is(result, "TwasWeights"))
    expect_equal(methodNames(result), "lasso")
})

test_that("learnTwasWeights: zero-variance columns are filtered, weighted 0", {
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

    # The returned weight matrix should have rows equal to total columns
    # (including zero-var)
    expect_equal(nrow(.weightsByMethod(result, "lassoWeights")), p_with_extra)
    # The zero-var column weight should be 0 (padded back)
    expect_equal(
        unname(.weightsByMethod(result, "lassoWeights")["chr1:9900:A:G", 1]),
        0
    )
})

test_that("learnTwasWeights: weights are ncol(X) x ncol(Y), rows named by X", {
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
        dim(.weightsByMethod(result, "enetWeights")),
        c(ncol(d$X), ncol(d$Y))
    )
    expect_equal(
        rownames(.weightsByMethod(result, "enetWeights")),
        colnames(d$X)
    )
})

test_that("learnTwasWeights: multiple methods return one entry per method", {
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
    expect_equal(length(methodNames(result)), 2)
    expect_true("lasso" %in% methodNames(result))
    expect_true("enet" %in% methodNames(result))
})

# ===========================================================================
#
#  learnTwasWeights with real glmnet (lasso/enet)
#
# ===========================================================================

test_that("learnTwasWeights: lasso and enet fit p x 1 weights with glmnet", {
    skip_if_not_installed("glmnet")
    d <- makeData(n = 50, p = 10)
    # One snake_case and one camelCase key: both spellings resolve.
    keys <- c(lasso = "lasso_weights", enet = "enetWeights")
    for (method in names(keys)) {
        result <- learnTwasWeights(
            d$X,
            d$Y,
            weightMethods = set_names(list(list()), keys[[method]])
        )
        w <- .weightsByMethod(result, method)
        expect_true(is(result, "TwasWeights"), info = method)
        expect_equal(methodNames(result), method, info = method)
        expect_equal(dim(w), c(ncol(d$X), 1L), info = method)
        # At least some weights should be non-zero for this strong signal
        expect_true(any(w != 0), info = method)
    }
})

# ===========================================================================
#
#  twasWeightsCv  (input validation)
#
# ===========================================================================

test_that("twasWeightsCv: NULL methods return only a well-formed partition", {
    d <- makeData()
    result <- twasWeightsCv(d$X, d$Y, fold = 5, weightMethods = NULL)
    expect_equal(names(result), "samplePartition")
    sp <- result$samplePartition

    expect_true(is.data.frame(sp))
    expect_true("Sample" %in% colnames(sp))
    expect_true("Fold" %in% colnames(sp))
    expect_equal(nrow(sp), nrow(d$X))
    expect_equal(length(unique(sp$Fold)), 5)
    # All sample names should appear
    expect_true(all(rownames(d$X) %in% sp$Sample))
})

test_that("twasWeightsCv: character weightMethods are accepted", {
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
#  learnTwasWeights with supplied SuSiE fits; twasWeightsCv plumbing
#
# ===========================================================================

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
        fittedModels = list(
            susie = makeFakeSusieFit(p = 10, maxNumSingleEffects = 5)
        )
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
            susie = makeFakeSusieFit(p = 10, maxNumSingleEffects = 5),
            susieInf = makeFakeSusieFit(
                p = 10,
                maxNumSingleEffects = 7,
                inf = TRUE
            )
        )
    )
    expect_equal(methodNames(result), c("susie", "susie_inf"))
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
# learnTwasWeights: a fit with extra rows is cut back to the fitted variants
# ===========================================================================

test_that("learnTwasWeights: mr.mash rows beyond X's variants are dropped", {
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
            # Return more rows than the p fitted variants, so the result is
            # subset back to those variants' names.
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
    # Restricted to the p fitted variants -> shape p x ncol(Y)
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
    # fold defaults to NULL, so this is also the call that omits it.
    expect_error(twasWeightsCv(X, y, fold = NULL), "fold.*samplePartitions")
    expect_error(
        twasWeightsCv(X, y, fold = "invalid"),
        "Must be of type 'count'"
    )
    expect_error(twasWeightsCv(X, y, fold = -1), "Must be >= 1")
    expect_error(twasWeightsCv(2, y, fold = 2), "Must be of type 'matrix'")
    expect_error(twasWeightsCv(X, 2, fold = 2), "Y.*Must have exactly 10 rows")
})


test_that("learnTwasWeights handles errors appropriately", {
    sim <- generateXY(seed = 1)
    expect_error(
        learnTwasWeights(sim$X, sim$Y),
        "weightMethods.*is missing"
    )
})

# ===========================================================================
# twasPredict: vector weights
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
        studyName = c("s1", "s1"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        method = c("lasso", "enet"),
        entry = list(e1, e2)
    )
    expect_s4_class(tw, "TwasWeights")
    expect_equal(nrow(tw), 2L)
    expect_setequal(methodNames(tw), c("lasso", "enet"))
})


test_that("TwasWeights: studyName / contexts / traitNames / methodNames", {
    e <- .sc_makeTwasWeightsRow()
    tw <- TwasWeights(
        studyName = c("s1", "s2"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        method = c("lasso", "lasso"),
        entry = list(e, e)
    )
    expect_setequal(studyName(tw), c("s1", "s2"))
    expect_setequal(contexts(tw), c("c1", "c2"))
    expect_equal(traitNames(tw), "t1")
    expect_equal(methodNames(tw), "lasso")
})


test_that("TwasWeights: rejects duplicate 4-tuples", {
    e <- .sc_makeTwasWeightsRow()
    expect_error(
        TwasWeights(
            studyName = c("s1", "s1"),
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
        studyName = c("s1", "s1"),
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
        studyName = c("s1", "s1"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        method = c("mrmash", "mrmash"),
        entry = list(e, e),
        jointContexts = c("c1;c2", "c1;c3")
    )
    expect_equal(nrow(tw2), 2L)
})


test_that("TwasWeights: twasWeights extracts the entry for a tuple", {
    e1 <- .sc_makeTwasWeightsRow()
    e2 <- .sc_makeTwasWeightsRow()
    tw <- TwasWeights(
        studyName = c("s1", "s1"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        method = c("lasso", "enet"),
        entry = list(e1, e2)
    )
    # twasWeights() returns the single-row COLLECTION for that tuple, not a
    # detached entry object, so compare what the row carries.
    picked <- twasWeights(
        tw,
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "enet"
    )
    expect_s4_class(picked, "TwasWeights")
    expect_equal(nrow(picked), 1L)
    expect_equal(as.character(picked$method), "enet")
    expect_identical(variantIds(picked), .twrPartsVariantIds(e2))
    expect_equal(unname(weights(picked)), unname(weights(e2)))
})

test_that("TwasWeights: rows are picked by window through subsetting", {
    near <- TwasWeightsRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
        weights = c(0.1, 0.2, 0.3, 0.4)
    )
    far <- TwasWeightsRow(
        variantIds = sprintf("chr2:%d:A:G", 100L * (1:4)),
        weights = c(0.5, 0.6, 0.7, 0.8)
    )
    tw <- TwasWeights(
        studyName = c("s1", "s1"),
        context = c("c1", "c1"),
        trait = c("t1", "t2"),
        method = c("lasso", "lasso"),
        entry = list(near, far)
    )
    # A window selects rows by subsetting the collection on its rows'
    # analysis windows; the one row left needs no other selector.
    inChr2 <- tw[IRanges::overlapsAny(
        genomicRegion(tw),
        GenomicRanges::GRanges("chr2:150-250")
    )]
    expect_identical(as.character(twasWeights(inChr2)$trait), "t2")
    expect_equal(resolveWeights(inChr2)$weights, c(0.5, 0.6, 0.7, 0.8))
    inChr1 <- tw[IRanges::overlapsAny(
        genomicRegion(tw),
        GenomicRanges::GRanges("chr1:1-1000")
    )]
    expect_equal(unname(weights(inChr1)), c(0.1, 0.2, 0.3, 0.4))
    expect_identical(
        variantIds(inChr1),
        sprintf("chr1:%d:A:G", 100L * (1:4))
    )
    # A window over both rows still needs the tuple to choose between them.
    both <- tw[IRanges::overlapsAny(
        genomicRegion(tw),
        GenomicRanges::GRanges(c("chr1:1-1000", "chr2:1-1000"))
    )]
    expect_error(twasWeights(both), "has 2 entries")
})

# === Tests migrated from test_showMethods.R (TwasWeights) ===

test_that("show.TwasWeights prints entry/study/context/trait/method counts", {
    e <- .sh_makeTwEntry()
    tw <- TwasWeights(
        studyName = c("s1", "s1"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        method = c("lasso", "enet"),
        entry = list(e, e)
    )
    out <- capture.output(show(tw))
    expect_true(any(grepl("TwasWeights: 2 entries", out)))
    expect_true(any(grepl("1 studies.*2 contexts.*1 traits.*2 methods", out)))
})


test_that("show.TwasWeights reports ldReference when present", {
    e <- .sh_makeTwEntry()
    tw <- TwasWeights(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "lasso",
        entry = list(e),
        ldReference = .sh_makeGenotypeHandle()
    )
    out <- capture.output(show(tw))
    expect_true(any(grepl("LD reference: plink2 @ /tmp/test", out)))
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
            studyName = c("s1", "s2"),
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
            studyName = "s1",
            context = "c1",
            trait = "t1",
            method = "lasso",
            entry = list(e),
            jointStudies = c("a", "b")
        ),
        "same length as"
    )
})

test_that("TwasWeights: weightStandardized/weightsDataType/variantIds delegate to the entry", {
    e1 <- TwasWeightsRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
        weights = rnorm(4),
        weightStandardized = TRUE,
        weightsDataType = "expression"
    )
    e2 <- TwasWeightsRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
        weights = rnorm(4),
        weightStandardized = FALSE,
        weightsDataType = "splicing"
    )
    tw <- TwasWeights(
        studyName = c("s1", "s1"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        method = c("lasso", "enet"),
        entry = list(e1, e2)
    )

    expect_true(weightStandardized(
        tw,
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "lasso"
    ))
    expect_false(weightStandardized(
        tw,
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "enet"
    ))
    expect_equal(
        weightsDataType(
            tw,
            studyName = "s1",
            context = "c1",
            trait = "t1",
            method = "enet"
        ),
        "splicing"
    )
    expect_equal(
        variantIds(
            tw,
            studyName = "s1",
            context = "c1",
            trait = "t1",
            method = "lasso"
        ),
        sprintf("chr1:%d:A:G", 100L * (1:4))
    )
})

test_that("TwasWeights: per-tuple weights/cvResult/methodFits/ldReference", {
    w1 <- rnorm(4)
    e1 <- TwasWeightsRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
        weights = w1,
        methodFits = list(tag = "fitA"),
        cvResult = list(rsq = 0.42)
    )
    e2 <- TwasWeightsRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
        weights = rnorm(4)
    )
    tw <- TwasWeights(
        studyName = c("s1", "s1"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        method = c("lasso", "enet"),
        entry = list(e1, e2),
        ldReference = .sh_makeGenotypeHandle()
    )

    expect_equal(
        weights(
            tw,
            studyName = "s1",
            context = "c1",
            trait = "t1",
            method = "lasso"
        ),
        w1
    )
    expect_equal(
        cvResult(
            tw,
            studyName = "s1",
            context = "c1",
            trait = "t1",
            method = "lasso"
        )$rsq,
        0.42
    )
    expect_equal(
        methodFits(
            tw,
            studyName = "s1",
            context = "c1",
            trait = "t1",
            method = "lasso"
        )$tag,
        "fitA"
    )
    expect_length(ldReference(tw), 1L)
    expect_s4_class(ldReference(tw)[[1L]], "RangedSummarizedExperiment")
})

test_that("TwasWeights: rejects non-TwasWeightsRow rows", {
    expect_error(
        TwasWeights(
            studyName = "s1",
            context = "c1",
            trait = "t1",
            method = "lasso",
            entry = list("not_an_entry")
        ),
        "must be a TWAS-weight row"
    )
})

test_that("TwasWeights: validity does not recurse on key subset (#546)", {
    # Building the tuple-uniqueness key via `object[, keyCols]` preserves the
    # TwasWeights class while dropping the required `entry` column; older
    # S4Vectors revalidates that intermediate and fails with "missing columns".
    e <- .sc_makeTwasWeightsRow(p = 5L)
    res <- TwasWeights(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "lasso",
        entry = list(e)
    )
    expect_s4_class(res, "TwasWeights")
    expect_true(validObject(res))
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
    infFit <- makeFakeSusieFit(p = 8, maxNumSingleEffects = 3, inf = TRUE)
    susieFit <- makeFakeSusieFit(p = 8, maxNumSingleEffects = 5)

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
    # branch (verbose>=1) is exercised. The fold sampling creates a fresh
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

test_that("learnTwasWeights: parallel weights path (two workers)", {
    skip_on_os("windows")
    d <- makeData(n = 30, p = 6)
    local_mocked_bindings(
        lassoWeights = function(X, y, ...) rep(0.1, ncol(X)),
        enetWeights = function(X, y, ...) rep(0.2, ncol(X))
    )
    result <- suppressMessages(learnTwasWeights(
        d$X,
        d$Y,
        weightMethods = list(lasso_weights = list(), enetWeights = list()),
        BPPARAM = BiocParallel::MulticoreParam(2)
    ))
    expect_true(is(result, "TwasWeights"))
    expect_setequal(methodNames(result), c("lasso", "enet"))
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
    e1 <- TwasWeightsRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * seq_len(p)),
        weights = w1
    )
    e2 <- TwasWeightsRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * seq_len(p)),
        weights = w2
    )
    tw <- TwasWeights(
        studyName = c("s1", "s1"),
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
        exec(twasWeightsCv, !!!list_modify(base, !!!list(BPPARAM = 2))),
        "must be a BiocParallelParam"
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
        exec(
            learnTwasWeights,
            !!!list_modify(base, !!!list(studyName = 1L))
        ),
        "studyName.*Must be of type 'string'"
    )
    # Called directly, not via modifyList(): modifyList() DROPS an element
    # whose value is NULL, so the argument would fall back to its default.
    expect_error(
        learnTwasWeights(
            d$X,
            d$Y,
            weightMethods = list(),
            weightStandardized = NULL
        ),
        "weightStandardized.*Must be of type 'logical flag'"
    )
    expect_error(
        exec(
            learnTwasWeights,
            !!!list_modify(base, !!!list(weightsDataType = 1L))
        ),
        "weightsDataType.*Must be of type 'string'"
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

test_that("twasWeightsCv cross-validates susie from the retained fold fits", {
    f <- .twcvFixture
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

test_that(".twasFoldFit injects the fold's fit, else leaves the args alone", {
    f <- pecotmr:::.twasFoldFit
    ff <- list(susie = list(fold_1 = "FIT1", fold_2 = "FIT2"))
    a <- f(list(), "susie_weights", 2L, ff)
    expect_identical(a$susieFit, "FIT2")
    args <- list(X = 1)
    # No CV payload at all.
    expect_identical(f(args, "susie_weights", 1L, NULL), args)
    # A method with no fine-mapping adapter is untouched.
    expect_identical(f(args, "lasso_weights", 1L, ff), args)
    # The token is present but carries no folds, so there is nothing to add.
    expect_identical(f(args, "susie", 1L, list(susie = list())), args)
})

test_that("twasWeightsCv refuses fold fits from a different partition", {
    f <- .twcvFixture
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
    unkeyed <- f$slice$foldFits
    attr(unkeyed, "partitionKey") <- NULL
    expect_error(
        suppressMessages(twasWeightsCv(
            f$X,
            f$Y,
            weightMethods = wm,
            verbose = 0,
            samplePartitions = f$slice$samplePartition,
            fittedModelsCv = list(susie = unkeyed)
        )),
        "no partition fingerprint"
    )

    # (c) a genuinely different split is caught by the fingerprint. Drawing
    # one needs no fits: with no methods, twasWeightsCv returns just the
    # partition.
    other <- suppressMessages(twasWeightsCv(f$X, f$Y, fold = 3, seed = 77))
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
    # The susie + susieInf pair is refused too, naming both: the chained
    # susieInf -> susie fit lives in fineMappingPipeline(), and
    # fitSusieInfThenSusie() itself is covered in test_fineMappingWrappers.R.
    expect_error(
        learnTwasWeights(
            X,
            Y,
            weightMethods = list(
                susie_weights = list(),
                susie_inf_weights = list()
            ),
            verbose = 0
        ),
        "susie, susieInf"
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
    f <- .twcvFixture
    expect_error(
        suppressMessages(twasWeightsCv(
            f$X,
            f$Y,
            samplePartitions = f$slice$samplePartition,
            weightMethods = list(susie_weights = list()),
            verbose = 0
        )),
        "never run fine-mapping themselves.*needs that fold's own fit"
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
        ".penalizedRssPath",
        "lassosumRss",
        "sdpr",
        ".sdprRun",
        "prsCs",
        ".prsCsRun"
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
    ids <- sprintf("chr1:%d:A:G", 100L * seq_len(p))
    list(
        ldReference = LdMatrix(`dimnames<-`(LD, list(ids, ids))),
        stat = list(
            b = rnorm(p) * 0.1,
            seb = rep(0.05, p),
            n = rep(500L, p),
            variantNames = ids
        )
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
    recordedFor <- function(tk) {
        pecotmr:::.twasMethodChainFor(tk, "QtlSumStats")
    }
    # Every pecotmr hop a run entered must be in the recorded chain, and the
    # chain's last pecotmr hop must be one the run entered.
    seen <- .tw_traceChain(scadRssWeights(fx$stat, fx$ldReference))
    expect_true(all(is_in(seen, recordedFor("scad"))))
    expect_true(is_in(".penalizedRssPath", seen))
    expect_true(is_in(".rssShrinkGridWeights", seen))
    seen <- .tw_traceChain(lassosumRssWeights(fx$stat, fx$ldReference))
    expect_true(all(is_in(seen, recordedFor("lasso"))))
    expect_true(is_in(".penalizedRssPath", seen))
    seen <- .tw_traceChain(enetRssWeights(fx$stat, fx$ldReference))
    expect_true(all(is_in(seen, recordedFor("enet"))))
    expect_true(is_in(".penalizedRssPath", seen))
    seen <- .tw_traceChain(ridgeRssWeights(fx$stat, fx$ldReference))
    expect_true(all(is_in(seen, recordedFor("ridge"))))
    expect_true(is_in(".penalizedRssPath", seen))
    prsCsOpts <- PrsCsOptions(nIter = 20, nBurnin = 5, thin = 1)
    seen <- .tw_traceChain(
        prsCsWeights(fx$stat, fx$ldReference, methodArgs = prsCsOpts)
    )
    expect_equal(seen, ".prsCsRun")
    expect_true(is_in(".prsCsRun", recordedFor("prsCs")))
    sdprOpts <- SdprOptions(iter = 20, burn = 5, thin = 1, verbose = FALSE)
    seen <- .tw_traceChain(
        sdprWeights(fx$stat, fx$ldReference, methodArgs = sdprOpts)
    )
    expect_equal(seen, ".sdprRun")
    expect_true(is_in(".sdprRun", recordedFor("dprGibbs")))
})

test_that("an RSS method accepts its own penalty settings and no other's", {
    # The penalty reaches .penalizedRssPath() as the wrapper's record, so
    # every chain shares that hop without taking each other's settings.
    accepted <- function(tk) pecotmr:::.twasChainAccepted(tk, "QtlSumStats")
    path <- c("lambda", "thr", "maxiter")
    l0Only <- c("lambda0", "lambda2", "maxSwaps")
    expect_true(all(is_in(path, accepted("lasso"))))
    expect_false(any(is_in(
        c("alpha", "gamma", "penalty", l0Only, "spec", "label", "method"),
        accepted("lasso")
    )))
    expect_true(all(is_in(c("alpha", path), accepted("enet"))))
    expect_false(any(is_in(c("gamma", l0Only), accepted("enet"))))
    expect_false(any(is_in(c("alpha", "gamma", l0Only), accepted("ridge"))))
    expect_true(all(is_in(c("alpha", "gamma"), accepted("scad"))))
    expect_false(any(is_in(l0Only, accepted("scad"))))
    expect_true(all(is_in(
        c("penalty", "lambda", l0Only),
        accepted("l0learn")
    )))
    expect_false(any(is_in(c("alpha", "gamma"), accepted("l0learn"))))
})

test_that("every RSS chain yields a checkable argument set", {
    # The point of the chains: each summary-statistics path ends in one of
    # pecotmr's own solvers, which enumerate their formals, so the RSS side
    # is checkable even where the individual side is not.
    tokens <- c(
        "scad",
        "mcp",
        "l0learn",
        "lasso",
        "enet",
        "ridge",
        "prsCs",
        "dprGibbs"
    )
    for (tk in tokens) {
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

# The exported *Options constructor that method `token` resolves to on input
# path `path`, by name: NA when the method does not run on that path, "?"
# when what it resolves to is not an exported *Options constructor.
optionsCtorName <- function(token, path) {
    ctor <- pecotmr:::.twasMethodCtorFor(token, path)
    if (is.null(ctor)) {
        return(NA_character_)
    }
    exports <- getNamespaceExports("pecotmr")
    candidates <- exports[str_ends(exports, "Options")]
    isCtor <- map_lgl(
        mget(candidates, envir = asNamespace("pecotmr")),
        identical,
        ctor
    )
    if (any(isCtor)) candidates[isCtor][[1L]] else "?"
}

test_that("the TWAS Options registry is derived from each implementation", {
    # It used to be 18 hand-written token -> constructor pairs, which had to
    # restate a pairing the implementation already declares through its own
    # `methodArgs` default. The two paths of one method often reach different
    # packages, so a second table is a second thing to keep in step.
    expect_setequal(
        pecotmr:::.twasConfigurableMethods(),
        names(pecotmr:::.twasMethodCapabilities)
    )
    # Every configurable method resolves, on each path it runs on, to an
    # exported *Options constructor a caller can build its options with. The
    # pairs themselves are pinned in the next test.
    for (path in c("QtlDataset", "QtlSumStats")) {
        toks <- intersect(
            pecotmr:::.twasConfigurableMethods(),
            pecotmr:::.twasMethodsFor(path)
        )
        ctorNames <- map_chr(toks, optionsCtorName, path = path)
        expect_equal(
            toks[is.na(ctorNames) | ctorNames == "?"],
            character(0),
            info = path
        )
    }
})

test_that("each method's Options constructor is resolved per input path", {
    # The two paths of one method can reach different PACKAGES, which is why
    # one constructor per method was never enough.
    expect_equal(optionsCtorName("lasso", "QtlDataset"), "GlmnetOptions")
    expect_equal(optionsCtorName("lasso", "QtlSumStats"), "LassosumOptions")
    expect_equal(optionsCtorName("mcp", "QtlDataset"), "NcvregOptions")
    expect_equal(optionsCtorName("bayesB", "QtlDataset"), "BglrOptions")
    expect_equal(optionsCtorName("dprGibbs", "QtlSumStats"), "SdprOptions")
    expect_equal(optionsCtorName("mrmash", "QtlSumStats"), "MrmashRssOptions")
    # enet and ridge run on summary statistics through the penalized solver.
    expect_equal(optionsCtorName("enet", "QtlSumStats"), "PenalizedRssOptions")
    expect_equal(optionsCtorName("ridge", "QtlDataset"), "GlmnetOptions")
    # prsCs has no individual implementation; its one constructor is the
    # summary-statistics one.
    expect_equal(optionsCtorName("prsCs", "QtlSumStats"), "PrsCsOptions")
    # A method that does not run on a path has no constructor there.
    expect_null(pecotmr:::.twasMethodCtorFor("dprVb", "QtlSumStats"))
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
    # A token no weightMethods name spells has no arguments to look up.
    expect_null(pecotmr:::.twasMethodArgsFor(list(), "susie"))
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

# ===========================================================================
# metadata, and the variants each row's fits were made on
# ===========================================================================

.tw_ids <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")

# A TwasWeights of one row per payload in `rows`, named r1, r2, ...
.tw_collection <- function(rows, metadata = list()) {
    n <- length(rows)
    TwasWeights(
        studyName = rep("s1", n),
        context = rep("c1", n),
        trait = sprintf("g%d", seq_len(n)),
        method = rep("lasso", n),
        entry = rows,
        metadata = metadata
    )
}

test_that("TwasWeights() keeps the caller's metadata", {
    tw <- .tw_collection(
        list(TwasWeightsRow(.tw_ids, c(0.1, 0.2, 0.3))),
        metadata = list(run = "r1")
    )
    expect_equal(S4Vectors::metadata(tw), list(run = "r1"))
})

test_that("a weights row records its fit variants only with a fit payload", {
    tw <- .tw_collection(list(
        TwasWeightsRow(.tw_ids, c(0.1, 0.2, 0.3), cvResult = list(r2 = 0.4)),
        TwasWeightsRow(.tw_ids, c(0.1, 0.2, 0.3))
    ))
    expect_s4_class(tw$fitVariants, "CharacterList")
    expect_equal(as.character(tw$fitVariants[[1L]]), .tw_ids)
    expect_equal(as.character(tw$fitVariants[[2L]]), character(0))
})

test_that("cvResult() and methodFits() of a narrowed row warn", {
    fits <- list(lasso = list(lambda = 0.1))
    cv <- list(r2 = 0.4)
    tw <- .tw_collection(list(
        TwasWeightsRow(.tw_ids, c(0.1, 0.2, 0.3), fits, cv)
    ))
    narrowed <- endoapply(tw, utils::head, 2L)
    expect_warning(
        out <- cvResult(narrowed),
        "cvResult\\(\\): the stored fit was made with 3 variants, 1 of"
    )
    expect_identical(out, cv)
    expect_warning(
        outFits <- methodFits(narrowed),
        "methodFits\\(\\): .* Refit on the variants kept"
    )
    expect_identical(outFits, fits)
    expect_no_warning(cvResult(tw))
    expect_no_warning(methodFits(tw))
})

test_that("updateObject() backfills a fitted row's held variants", {
    tw <- .tw_collection(list(
        TwasWeightsRow(.tw_ids, c(0.1, 0.2, 0.3), cvResult = list(r2 = 0.4)),
        TwasWeightsRow(.tw_ids, c(0.1, 0.2, 0.3))
    ))
    saved <- tw
    saved$fitVariants <- NULL
    up <- updateObject(saved)
    expect_equal(as.character(up$fitVariants[[1L]]), .tw_ids)
    expect_equal(as.character(up$fitVariants[[2L]]), character(0))
    # A current collection is left as it is.
    expect_identical(updateObject(tw)$fitVariants, tw$fitVariants)
    expect_message(updateObject(tw, verbose = TRUE), "TwasWeights")
})

test_that("updateObject() backfills held variants and keys LD rows at once", {
    tw <- .tw_collection(list(
        TwasWeightsRow(.tw_ids, c(0.1, 0.2, 0.3), cvResult = list(r2 = 0.4))
    ))
    saved <- tw
    saved$fitVariants <- NULL
    up <- updateObject(.savedWithLdSketch(saved, NULL))
    expect_true(validObject(up))
    expect_equal(as.character(up$fitVariants[[1L]]), .tw_ids)
    expect_identical(ldReference(up), list())
    expect_true(is.na(up$ldReference))
    expect_identical(updateObject(up), up)
})

# ===========================================================================
# estimateSparsity
# ===========================================================================

test_that("estimateSparsity: legacy list input reads attr(.,'fit')$pi", {
    # Build a weight result that mimics learnTwasWeights keeping its fits:
    # element name carries the `_weights` suffix; attr 'fit' carries an mr.ash
    # object whose pi[1] is the spike weight.
    fake_w <- structure(c(0.1, 0, 0.3), fit = list(pi = c(0.6, 0.2, 0.2)))
    weightResults <- list(mrash_weights = fake_w)
    expect_equal(estimateSparsity(weightResults), 1 - 0.6, tolerance = 1e-12)
})

test_that("estimateSparsity: TwasWeights collection input reads from the mrash entry", {
    entry <- TwasWeightsRow(
        variantIds = c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G"),
        weights = c(0.1, 0, 0.3),
        methodFits = list(pi = c(0.7, 0.2, 0.1)),
        weightStandardized = FALSE
    )
    tw <- TwasWeights(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "mrash",
        entry = list(entry)
    )
    expect_equal(estimateSparsity(tw), 1 - 0.7, tolerance = 1e-12)
})

test_that("estimateSparsity: TwasWeights without mrash entry errors", {
    entry <- TwasWeightsRow(variantIds = c("chr1:100:A:G"), weights = c(0.1))
    tw <- TwasWeights(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "lasso",
        entry = list(entry)
    )
    expect_error(
        estimateSparsity(tw),
        "mr.ash entry not found in TwasWeights"
    )
})

test_that("estimateSparsity: TwasWeights with mrash entry but no fit$pi errors", {
    entry <- TwasWeightsRow(
        variantIds = c("chr1:100:A:G"),
        weights = c(0.1),
        methodFits = list(other = 1)
    )
    tw <- TwasWeights(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "mrash",
        entry = list(entry)
    )
    expect_error(
        estimateSparsity(tw),
        "mr.ash fit object not found"
    )
})

test_that("estimateSparsity: legacy list input without mrash_weights errors", {
    expect_error(
        estimateSparsity(list(lasso = c(0.1, 0.2))),
        "'mrash_weights'.*not found"
    )
})

test_that("estimateSparsity: legacy list input without fit attr errors", {
    expect_error(
        estimateSparsity(list(mrash_weights = c(0.1, 0.2))),
        "mr.ash fit object not found"
    )
})

# ===========================================================================
# .twasMethodMap
# ===========================================================================

test_that("ser is reachable as a TWAS method", {
    # It is fit-derived like susie/mvsusie, so it carries no learner
    # arguments and its entry is an empty stub.
    expect_true(is_in("ser", names(pecotmr:::.twasMethodMap)))
    expect_equal(pecotmr:::.twasMethodMap$ser$impl, "serWeights")
    expect_length(pecotmr:::.twasMethodMap$ser$args, 0L)
    expect_setequal(
        pecotmr:::.twasNormalizeMethods(c("ser", "lasso"))$tokens,
        c("ser", "lasso")
    )
    # Nothing to configure, same as the other fit-derived methods.
    expect_error(
        TwasWeightsMethodsParam(
            qtlSumStatsMethods = list(ser = list(nope = 1))
        ),
        "nothing to configure here"
    )
})
