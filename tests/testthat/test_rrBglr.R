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

# ---- several conditions at once, through BGLR's Multitrait() ----
# Two conditions sharing one effect, the first unmeasured on 15 samples.
.bglrConditions <- function() {
    set.seed(7)
    X <- matrix(rnorm(600), nrow = 60)
    colnames(X) <- sprintf("v%d", seq_len(10))
    Y <- cbind(a = X[, 1] * 0.6 + rnorm(60), b = X[, 1] * 0.5 + rnorm(60))
    Y[1:15, "a"] <- NA
    list(X = X, Y = Y)
}

test_that("a matrix y runs Multitrait(), the conditions' weights its columns", {
    skip_if_not_installed("BGLR")
    d <- .bglrConditions()
    w <- bayesCWeights(
        d$X,
        d$Y,
        backend = "bglr",
        nIter = 300,
        nBurnin = 100,
        thin = 2,
        fitRetention = "slim"
    )
    expect_equal(dim(w), c(10L, 2L))
    expect_equal(dimnames(w), list(colnames(d$X), c("a", "b")))
    expect_true(all(is.finite(w)))
    # Both conditions find the shared effect, the first on its 45 samples.
    expect_equal(unname(apply(abs(w), 2L, which.max)), c(1L, 1L))
    fit <- attr(w, "fit")
    expect_equal(dim(fit$coef), c(10L, 2L))
    expect_equal(dim(fit$pip), c(10L, 2L))
    expect_equal(fit$variantIds, colnames(d$X))
    gaussian <- bayesNWeights(
        d$X,
        d$Y,
        backend = "bglr",
        nIter = 300,
        nBurnin = 100,
        thin = 2,
        fitRetention = "full"
    )
    expect_equal(dim(gaussian), c(10L, 2L))
    expect_equal(attr(gaussian, "fit")$ETA[[1L]]$model, "BRR")
    # A one-column matrix is that column, fitted by BGLR() as before.
    one <- bayesCWeights(
        d$X,
        d$Y[, "b", drop = FALSE],
        backend = "bglr",
        nIter = 100,
        nBurnin = 20
    )
    expect_null(dim(one))
    expect_length(one, 10L)
})

test_that("Multitrait()'s linear term gets Cov and every condition's probIn", {
    skip_if_not_installed("BGLR")
    d <- .bglrConditions()
    # BGLR()'s own residual prior (df0) would be an unused argument to
    # Multitrait(), which therefore never sees it.
    w <- bayesCWeights(
        d$X,
        d$Y,
        backend = "bglr",
        pi = 0.3,
        nIter = 100,
        nBurnin = 20,
        fitRetention = "full",
        methodArgs = BglrOptions(
            Cov = list(type = "DIAG"),
            resCov = list(type = "DIAG"),
            df0 = 4
        )
    )
    fit <- attr(w, "fit")
    term <- fit$ETA[[1L]]
    expect_equal(term$model, "SpikeSlab")
    expect_equal(term$Cov$type, "DIAG")
    expect_equal(fit$resCov$type, "DIAG")
    # Every condition starts at pi; Multitrait() fails on a probIn given
    # without its prior counts.
    expect_equal(term$inclusionProb$counts, c(2, 2))
    expect_equal(
        .bglrMultitraitTerm("SpikeSlab", list(probIn = 0.3), BglrOptions(), 2L),
        list(inclusionProb = list(probIn = c(0.3, 0.3), counts = c(2, 2)))
    )
    expect_equal(
        .bglrMultitraitTerm("BRR", list(probIn = 0.3), BglrOptions(), 2L),
        list()
    )
})

test_that("BglrOptions serves BGLR() and Multitrait(), each its own options", {
    skip_if_not_installed("BGLR")
    opts <- BglrOptions(
        resCov = list(type = "UN"),
        intercept = TRUE,
        Cov = list(type = "DIAG"),
        R2 = 0.4
    )
    expect_s4_class(opts, "MethodOptions")
    expect_error(BglrOptions(resCovv = 1), "unknown argument")
    expect_error(BglrOptions(nIter = 10), "nIter")
    # One phenotype ignores Multitrait()'s settings, which would be unused
    # arguments to BGLR(), and keeps its own.
    d <- .bglrConditions()
    w <- bayesCWeights(
        d$X,
        d$Y[, "b"],
        backend = "bglr",
        nIter = 100,
        nBurnin = 20,
        fitRetention = "full",
        methodArgs = opts
    )
    expect_length(w, 10L)
    expect_equal(attr(w, "fit")$ETA[[1L]]$R2, 0.4)
})

test_that("only BGLR's bayesC and bayesN fit a matrix y", {
    skip_if_not_installed("BGLR")
    d <- .bglrConditions()
    expect_error(
        bayesCWeights(d$X, d$Y, backend = "qgg"),
        "only backend = \"bglr\" fits a matrix"
    )
    expect_error(
        bayesBWeights(d$X, d$Y, nIter = 100, nBurnin = 20),
        "bayesB on bglr takes a vector"
    )
    expect_error(
        bglrWeights(d$X, d$Y, "BL", nIter = 100, burnIn = 20, thin = 1),
        "Multitrait\\(\\) has no form of BL"
    )
    expect_error(bayesCWeights(d$X, d$Y[-1, ], backend = "bglr"), "rows")
})

test_that("a matrix y's heritability estimate is its conditions' mean", {
    skip_if_not_installed("BGLR")
    d <- .bglrConditions()
    perCondition <- vapply(
        c("a", "b"),
        function(k) {
            seen <- !is.na(d$Y[, k])
            pecotmr:::.alphabetConditionH2(
                "fastHer",
                d$X[seen, ],
                d$Y[seen, k]
            )
        },
        numeric(1)
    )
    expect_equal(
        pecotmr:::.alphabetIndividualH2("fastHer", d$X, d$Y),
        mean(perCondition)
    )
    expect_equal(pecotmr:::.alphabetIndividualH2(0.3, d$X, d$Y), 0.3)
    expect_null(pecotmr:::.alphabetIndividualH2(NULL, d$X, d$Y))
    # No condition giving an estimate leaves the engine its own default.
    local_mocked_bindings(.alphabetColumnH2 = function(...) NULL)
    expect_null(pecotmr:::.alphabetIndividualH2("fastHer", d$X, d$Y))
})

test_that("BglrOptions checks nothing where BGLR's arguments are unreadable", {
    local_mocked_bindings(.engineAcceptedNames = function(...) NULL)
    expect_null(pecotmr:::.bglrAcceptedNames())
})

test_that("a joint fit across studies takes diagonal residuals", {
    opts <- BglrOptions(resCov = list(type = "UN", df0 = 5))
    across <- pecotmr:::.bglrStudyResiduals(opts, TRUE)
    expect_equal(as.list(across)$resCov, list(type = "DIAG"))
    expect_identical(pecotmr:::.bglrStudyResiduals(opts, FALSE), opts)
    # A diagonal one keeps its prior settings.
    diagonal <- BglrOptions(resCov = list(type = "DIAG", df0 = c(4, 4)))
    expect_identical(pecotmr:::.bglrStudyResiduals(diagonal, TRUE), diagonal)
    expect_equal(
        as.list(pecotmr:::.bglrStudyResiduals(BglrOptions(), TRUE))$resCov,
        list(type = "DIAG")
    )
})
