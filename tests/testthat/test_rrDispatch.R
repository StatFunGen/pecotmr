context("regularized_regression - dispatch verification")

# ---- dispatch verification ----
# These tests use local_mocked_bindings to replace the inner function with a
# stub that captures its arguments, then assert the wrapper forwarded the
# correct values. They catch silent dispatch bugs (wrong method, wrong penalty,
# dropped argument) that shape-only tests would not. The glmnet, ncvreg and
# BGLR wrappers are checked at the external solver itself, in
# test_rrE2eSolvers.R and test_rrBglr.R.

test_that("prsCsWeights hands the PRS-CS sampler the right arguments", {
    set.seed(42)
    p <- 10
    bhat <- rnorm(p, sd = 0.1)
    R <- diag(p)
    # Heterogeneous n with median != mean: this lets the test catch a regression
    # that swapped median(stat$n) for mean(stat$n) or sum(stat$n).
    # sorted: 10, 20, 30, 40, 50, 60, 70, 80, 90, 1000 -> median = 55, mean = 145.
    stat <- list(
        b = bhat,
        n = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 1000),
        variantNames = .rrwIds(p)
    )
    captured <- new.env(parent = emptyenv())
    local_mocked_bindings(
        .prsCsRun = function(bhat, ld, n, ...) {
            captured$bhat <- bhat
            captured$R <- ld$R
            captured$n <- n
            captured$dots <- list(...)
            list(betaEst = seq_len(length(bhat)) * 0.01)
        }
    )
    result <- prsCsWeights(
        stat = stat,
        ldReference = .rrwLdMatrix(R),
        methodArgs = PrsCsOptions(maf = rep(0.3, p), nIter = 17)
    )
    expect_equal(captured$bhat, bhat)
    expect_equal(unname(captured$R), R)
    expect_equal(captured$n, 55) # median of stat$n, NOT mean (145)
    expect_equal(captured$dots$maf, rep(0.3, p))
    expect_equal(captured$dots$nIter, 17)
    expect_equal(result, seq_len(p) * 0.01)
})

test_that("sdprWeights dispatches to the SDPR sampler with correct arguments", {
    set.seed(42)
    p <- 10
    bhat <- rnorm(p, sd = 0.1)
    R <- diag(p)
    stat <- list(b = bhat, n = rep(456, p), variantNames = .rrwIds(p))
    captured <- new.env(parent = emptyenv())
    local_mocked_bindings(
        .sdprRun = function(bhat, ld, n, ...) {
            captured$bhat <- bhat
            captured$R <- ld$R
            captured$n <- n
            captured$dots <- list(...)
            list(betaEst = seq_len(length(bhat)) * 0.02)
        }
    )
    result <- sdprWeights(
        stat = stat,
        ldReference = .rrwLdMatrix(R),
        methodArgs = SdprOptions(iter = 19, burn = 3)
    )
    expect_equal(captured$bhat, bhat)
    expect_equal(unname(captured$R), R)
    expect_equal(captured$n, 456)
    expect_equal(captured$dots$iter, 19)
    expect_equal(captured$dots$burn, 3)
    expect_equal(result, seq_len(p) * 0.02)
})

test_that("lassosumRssWeights runs the lassosum path once per s value", {
    set.seed(42)
    p <- 10
    bhat <- rnorm(p, sd = 0.1)
    R <- diag(p)
    for (i in seq_len(p - 1)) {
        R[i, i + 1] <- 0.4
        R[i + 1, i] <- 0.4
    }
    stat <- list(b = bhat, n = rep(100, p), variantNames = .rrwIds(p))
    call_log <- new.env(parent = emptyenv())
    call_log$calls <- list()
    local_mocked_bindings(
        .penalizedRssPath = function(bhat, ld, n, ...) {
            call_log$calls <- c(
                call_log$calls,
                list(list(bhat = bhat, ld = ld, n = n, spec = list(...)$spec))
            )
            list(
                beta = cbind(rep(0, length(bhat)), rep(0.05, length(bhat))),
                lambda = c(0.05, 0.01),
                fbeta = c(1.0, 0.5)
            )
        }
    )
    lassosumRssWeights(
        stat = stat,
        ldReference = .rrwLdMatrix(R),
        s = c(0.2, 0.9),
        selection = "minFbeta"
    )
    expect_equal(length(call_log$calls), 2L)
    expect_equal(call_log$calls[[1]]$n, 100)
    # lassosum is the penalized path at the lasso, alpha = 1.
    expect_equal(call_log$calls[[1]]$spec$penalty, "lasso")
    expect_equal(call_log$calls[[1]]$spec$alpha, 1)
    expect_equal(length(call_log$calls[[1]]$bhat), p)
    # Each call carries its own s for the solver to apply to the one LD
    # matrix, rather than a copy of the matrix shrunk for that s.
    expect_equal(map_dbl(call_log$calls, c("ld", "shrink")), c(0.2, 0.9))
    expect_identical(call_log$calls[[1]]$ld$R, call_log$calls[[2]]$ld$R)
    expect_equal(unname(call_log$calls[[1]]$ld$R), R)
})

test_that("lassosumRssWeights defaults to LD-quadratic selection", {
    p <- 4
    stat <- list(
        cor = c(0.4, 0.1, 0, 0),
        n = rep(100, p),
        variantNames = .rrwIds(p)
    )
    R <- diag(p)
    R[1, 2] <- 0.5
    R[2, 1] <- 0.5
    expected <- c(1, 0, 0, 0)

    local_mocked_bindings(
        .penalizedRssPath = function(bhat, ld, n, ...) {
            list(
                beta = cbind(expected, c(1, 1, 0, 0)),
                lambda = c(0.05, 0.01),
                fbeta = c(5, 1)
            )
        }
    )

    result <- lassosumRssWeights(
        stat = stat,
        ldReference = .rrwLdMatrix(R),
        s = 0.5
    )
    expect_equal(c(result), expected)
    expect_equal(
        unname(attr(result, "lassosum_selection")["mode"]),
        "ldQuadratic"
    )
})

test_that("lassosumRssWeights uses first-max tie behavior for LD-quadratic selection", {
    p <- 3
    stat <- list(
        cor = c(0.3, 0.3, 0),
        n = rep(100, p),
        variantNames = .rrwIds(p)
    )
    R <- diag(p)

    local_mocked_bindings(
        .penalizedRssPath = function(bhat, ld, n, ...) {
            list(
                beta = cbind(c(1, 0, 0), c(0, 1, 0)),
                lambda = c(0.05, 0.01),
                fbeta = c(2, 1)
            )
        }
    )

    result <- lassosumRssWeights(
        stat = stat,
        ldReference = .rrwLdMatrix(R),
        s = 0.5
    )
    expect_equal(c(result), c(1, 0, 0))
    expect_equal(
        unname(attr(result, "lassosum_selection")["mode"]),
        "ldQuadratic"
    )
})

test_that("lassosumRssWeights still supports explicit min(fbeta)", {
    p <- 4
    stat <- list(
        cor = c(0.4, 0.1, 0, 0),
        n = rep(100, p),
        variantNames = .rrwIds(p)
    )
    R <- diag(p)
    expected <- c(1, 1, 0, 0)

    local_mocked_bindings(
        .penalizedRssPath = function(bhat, ld, n, ...) {
            list(
                beta = cbind(c(1, 0, 0, 0), expected),
                lambda = c(0.05, 0.01),
                fbeta = c(5, 1)
            )
        }
    )

    result <- lassosumRssWeights(
        stat = stat,
        ldReference = .rrwLdMatrix(R),
        s = 0.5,
        selection = "minFbeta"
    )
    expect_equal(c(result), expected)
    expect_equal(
        unname(attr(result, "lassosum_selection")["mode"]),
        "minFbeta"
    )
})

test_that("bayes_{n,l,a,c,r}_weights each dispatch to bayesAlphabetWeights with correct method", {
    set.seed(42)
    X <- matrix(rnorm(50), nrow = 10)
    y <- rnorm(10)
    dispatchers <- list(
        list(fn = bayesNWeights, expected = "bayesN"),
        list(fn = bayesLWeights, expected = "bayesL"),
        list(fn = bayesAWeights, expected = "bayesA"),
        list(fn = bayesCWeights, expected = "bayesC"),
        list(fn = bayesRWeights, expected = "bayesR")
    )
    for (d in dispatchers) {
        captured <- new.env(parent = emptyenv())
        local_mocked_bindings(
            bayesAlphabetWeights = function(X, y, method, ...) {
                captured$method <- method
                rep(0, ncol(X))
            }
        )
        d$fn(X, y)
        expect_equal(
            captured$method,
            d$expected,
            label = paste("dispatch for", d$expected)
        )
    }
})


test_that("mrashWeights actually calls lassoWeights for default beta.init", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    y <- X[, 1] * 0.5 + rnorm(n)
    called <- new.env(parent = emptyenv())
    called$lasso <- FALSE
    local_mocked_bindings(
        lassoWeights = function(X, y) {
            called$lasso <- TRUE
            rep(0.01, ncol(X))
        }
    )
    result <- suppressWarnings(mrashWeights(X, y))
    expect_true(called$lasso)
    expect_true(is.numeric(result))
    expect_equal(length(result), p)
})

test_that("mrashWeights leaves the prior grid to mr.ash unless one is given", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    X <- matrix(rnorm(50 * 10), nrow = 50)
    y <- X[, 1] * 0.5 + rnorm(50)
    seen <- NULL
    real <- susieR::mr.ash
    local_mocked_bindings(
        mr.ash = function(...) {
            seen <<- list(...)
            real(...)
        }
    )
    # mr.ash's own default grid, built from the data's column scales.
    suppressWarnings(mrashWeights(X, y))
    expect_false("sa2" %in% names(seen))
    grid <- c(0, 0.1, 0.5)
    suppressWarnings(
        mrashWeights(X, y, methodArgs = MrashOptions(sa2 = grid))
    )
    expect_equal(seen$sa2, grid)
})

gc()
