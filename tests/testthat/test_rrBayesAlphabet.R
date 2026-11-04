context("regularized_regression - bayes_alphabet")

.baXy <- function(n = 50, p = 10) {
    set.seed(42)
    X <- matrix(rnorm(n * p), nrow = n)
    list(X = X, y = X[, 1] * 0.5 + rnorm(n), p = p)
}

# ---- bayesAlphabetWeights ----
test_that("bayesAlphabetWeights errors on dimension mismatch", {
    skip_if_not_installed("qgg")
    X <- matrix(rnorm(100), nrow = 10)
    y <- rnorm(5)
    expect_error(
        bayesAlphabetWeights(X, y, method = "bayesN"),
        "y.*Must have length 10"
    )
})

test_that("bayesAlphabetWeights runs each method on its default engine", {
    skip_if_not_installed("qgg")
    d <- .baXy()
    result <- bayesAlphabetWeights(
        d$X,
        d$y,
        method = "bayesN",
        nIter = 50,
        nBurnin = 10
    )
    expect_equal(length(result), d$p)
    expect_true(is.numeric(result))
})

test_that("each alphabet wrapper offers its method's backends, default first", {
    # The wrappers spell their choices in their signatures; the pipeline reads
    # them from .alphabetBackends. The two must not drift apart.
    wrappers <- c(
        bayesA = "bayesAWeights",
        bayesB = "bayesBWeights",
        bayesC = "bayesCWeights",
        bayesL = "bayesLWeights",
        bayesN = "bayesNWeights",
        bayesR = "bayesRWeights"
    )
    for (m in names(wrappers)) {
        expect_equal(
            eval(formals(get(wrappers[[m]]))$backend),
            names(pecotmr:::.alphabetBackends$individual[[m]]),
            info = m
        )
    }
})

test_that("bayesAlphabetWeights refuses settings the method or engine lacks", {
    d <- .baXy()
    expect_error(
        bayesAlphabetWeights(d$X, d$y, "bayesA", backend = "bglr"),
        "`backend` for bayesA must be one of qgg, hibayes"
    )
    expect_error(
        bayesAlphabetWeights(d$X, d$y, "bayesN", pi = 0.1),
        "`pi` applies to bayesB and bayesC, not bayesN"
    )
    # The options must come from the chosen engine's own constructor.
    expect_error(
        bayesBWeights(d$X, d$y, methodArgs = QggOptions()),
        "built for engine 'qgg', but backend 'bglr' takes BglrOptions"
    )
    expect_error(
        bayesAWeights(d$X, d$y, methodArgs = list(nub = 4)),
        "must be built with QggOptions"
    )
    expect_error(
        bayesAWeights(d$X, d$y, nIter = 100, nBurnin = 100),
        "must exceed `nBurnin`"
    )
})

test_that("pi reaches qgg once, from the wrapper", {
    skip_if_not_installed("qgg")
    d <- .baXy()
    seen <- NULL
    local_mocked_bindings(
        gbayes = function(...) {
            seen <<- list(...)
            list(bm = rep(0, 10))
        },
        .package = "qgg"
    )
    invisible(bayesCWeights(d$X, d$y, pi = 0.05, nIter = 20, nBurnin = 5))
    expect_equal(seen$pi, 0.05)
    expect_equal(seen$method, "bayesC")
    # Setting it in the engine's options too is the same option twice.
    expect_error(
        bayesCWeights(d$X, d$y, methodArgs = QggOptions(pi = 0.05)),
        "given twice: pi"
    )
})

test_that("the hibayes backend runs ibrm with the shared settings", {
    skip_if_not_installed("hibayes")
    d <- .baXy()
    seen <- NULL
    local_mocked_bindings(
        ibrm = function(...) {
            seen <<- list(...)
            list(alpha = rep(0.5, ncol(seen$M)))
        },
        .package = "hibayes"
    )
    w <- bayesCWeights(
        d$X,
        d$y,
        backend = "hibayes",
        nIter = 30,
        nBurnin = 7,
        thin = 3
    )
    expect_equal(w, rep(0.5, d$p))
    expect_equal(seen$method, "BayesCpi")
    expect_equal(seen$Pi, c(0.9, 0.1))
    expect_equal(c(seen$niter, seen$nburn, seen$thin), c(30, 7, 3))
    expect_equal(seen$threads, 1L)
    # ibrm() matches genotypes to phenotypes by the first column of `data`.
    expect_equal(seen$M.id, seen$data[[1L]])
    expect_error(
        bayesCWeights(
            d$X,
            d$y,
            "hibayes",
            methodArgs = HibayesOptions(seed = 1)
        ),
        "`seed` cannot be set here"
    )
})

test_that("a hibayes fit leaves the caller's random stream where it was", {
    skip_if_not_installed("hibayes")
    d <- .baXy(n = 40, p = 5)
    set.seed(9)
    w <- bayesNWeights(
        d$X,
        d$y,
        backend = "hibayes",
        nIter = 200,
        nBurnin = 50
    )
    after <- runif(1)
    expect_length(w, 5)
    expect_true(all(is.finite(w)))
    # The seed hibayes gets is one draw from R's generator; ibrm()'s own
    # set.seed() is undone afterwards.
    set.seed(9)
    invisible(sample.int(.Machine$integer.max, 1L))
    expect_equal(after, runif(1))
})

# ---- bayesAlphabetWeights zero-variance handling ----
test_that("bayesAlphabetWeights warns, zero-pads X's zero-variance columns", {
    skip_if_not_installed("qgg")
    d <- .baXy()
    X <- d$X
    X[, 4] <- 2
    expect_warning(
        result <- bayesAlphabetWeights(
            X,
            d$y,
            method = "bayesN",
            nIter = 50,
            nBurnin = 10
        ),
        "bayesAlphabetWeights: dropping 1 zero-variance column"
    )
    expect_equal(length(result), d$p)
    expect_equal(result[4], 0)
    expect_true(all(is.finite(result)))
})

# ---- bayesAlphabetRssWeights ----

test_that("bayesAlphabetWeights keeps the posterior or the engine's fit", {
    skip_if_not_installed("qgg")
    skip_if_not_installed("BGLR")
    skip_if_not_installed("hibayes")
    d <- .baXy()
    X <- d$X
    colnames(X) <- sprintf("chr1:%d:A:G", seq_len(d$p))
    runs <- list(
        list(method = "bayesC", backend = "qgg"),
        list(method = "bayesB", backend = "bglr"),
        list(method = "bayesR", backend = "hibayes")
    )
    fits <- list()
    for (r in runs) {
        set.seed(1)
        slim <- bayesAlphabetWeights(
            X,
            d$y,
            r$method,
            backend = r$backend,
            nIter = 200,
            nBurnin = 50,
            fitRetention = "slim"
        )
        fit <- attr(slim, "fit")
        expect_equal(fit$coef, as.numeric(slim))
        expect_equal(fit$variantIds, colnames(X))
        full <- bayesAlphabetWeights(
            X,
            d$y,
            r$method,
            backend = r$backend,
            nIter = 200,
            nBurnin = 50,
            fitRetention = "full"
        )
        fits[[r$backend]] <- attr(full, "fit")
        expect_equal(
            pecotmr:::.twasFitVariantIds(fits[[r$backend]]),
            colnames(X)
        )
    }
    # qgg and BGLR report inclusion probabilities, and name the effects.
    expect_length(attr(slim, "fit")$coef, d$p)
    expect_true(is.numeric(fits$qgg$dm))
    expect_s3_class(fits$bglr, "BGLR")
    expect_null(fits$bglr$variantIds)
    # hibayes names none, so its fit takes them once.
    expect_s3_class(fits$hibayes, "blrMod")
    expect_equal(fits$hibayes$variantIds, colnames(X))
})

test_that("a summary-statistics fit picks its engine by the LD reference", {
    d <- .rrwPanelStat()
    f <- pecotmr:::.rssBackend
    # A genotype reference runs pecotmr's own sampler on its factor; an LD
    # matrix runs qgg, or hibayes for bayesB and bayesL.
    expect_equal(f("bayesA", NULL, d$panel), "pecotmr")
    expect_equal(f("bayesA", NULL, d$ldMatrix), "qgg")
    expect_equal(f("bayesB", NULL, d$ldMatrix), "hibayes")
    expect_equal(f("bayesL", NULL, d$ldMatrix), "hibayes")
    # The matrix engines also run on a genotype reference, through the
    # matrix formed from it.
    expect_equal(f("bayesA", "hibayes", d$panel), "hibayes")
    expect_error(
        f("bayesA", "pecotmr", d$ldMatrix),
        "runs on a genotype reference's factor"
    )
    expect_error(
        f("bayesB", "qgg", d$ldMatrix),
        "must be one of hibayes, pecotmr"
    )
    # LDpred2 defaults as the alphabet does; the lasso's pecotmr engine reads
    # either LD, so it is the default and is never refused.
    expect_equal(f("ldpred2Auto", NULL, d$panel), "pecotmr")
    expect_equal(f("ldpred2Grid", NULL, d$ldMatrix), "bigsnpr")
    expect_error(
        f("ldpred2Grid", "pecotmr", d$ldMatrix),
        "runs on a genotype reference's factor"
    )
    expect_equal(f("lasso", NULL, d$ldMatrix), "pecotmr")
    expect_equal(f("lasso", "pecotmr", d$ldMatrix), "pecotmr")
    expect_equal(f("lasso", "bigsnpr", d$panel), "bigsnpr")
})

test_that("pecotmr's sampler lands where the matrix engines do", {
    skip_if_not_installed("qgg")
    skip_if_not_installed("hibayes")
    d <- .rrwPanelStatLd()
    run <- function(method, backend) {
        native <- identical(backend, "pecotmr")
        set.seed(1)
        suppressMessages(bayesAlphabetRssWeights(
            d$stat,
            if (native) d$panel else d$ldMatrix,
            method,
            backend = backend,
            nIter = 3000,
            nBurnin = 500,
            # Without SDPR's projection, which the matrix engines lack.
            methodArgs = if (native) {
                BayesAlphabetRssOptions(ldProjection = NULL)
            }
        ))
    }
    # Each method's priors follow one engine's, on the same LD. This LD has
    # rank 22 of 40, so the effects of variants in perfect LD are not
    # identified and two runs of one engine split them differently; their
    # fitted marginal effects R b are what the engines must agree on.
    ldFactor <- pecotmr:::.ldFactorFromRef(d$panel, d$stat$variantNames, "test")
    fitted <- function(b) as.vector(crossprod(ldFactor, ldFactor %*% b))
    agree <- function(method, other) {
        cor(fitted(run(method, "pecotmr")), fitted(run(method, other)))
    }
    expect_gt(agree("bayesA", "qgg"), 0.99)
    expect_gt(agree("bayesC", "qgg"), 0.99)
    expect_gt(agree("bayesN", "qgg"), 0.99)
    expect_gt(agree("bayesR", "qgg"), 0.99)
    expect_gt(agree("bayesB", "hibayes"), 0.99)
    expect_gt(agree("bayesL", "hibayes"), 0.99)
})

test_that("pecotmr's Legarra bayesL matches qgg's individual-level fit", {
    # qgg fits the Bayesian LASSO in Legarra et al.'s form on individual
    # data; its summary-statistics sampler uses a different update, so the
    # individual fit is the one the form is checked against.
    skip_if_not_installed("qgg")
    set.seed(31)
    n <- 600
    p <- 30
    X <- matrix(rbinom(n * p, 2, 0.3), n, p)
    for (j in 2:p) {
        X[, j] <- ifelse(runif(n) < 0.5, X[, j - 1], X[, j])
    }
    effects <- replace(numeric(p), c(5, 18), c(0.25, -0.2))
    Xs <- scale(X)
    y <- as.vector(Xs %*% effects + rnorm(n, sd = 0.95))
    ldFactor <- pecotmr:::.ldFactorFromDosage(X)
    native <- pecotmr:::bayesAlphabetRss(
        as.vector(cor(Xs, y)),
        ldFactor,
        n,
        "bayesL",
        nIter = 3000,
        nBurnin = 500,
        bayesLForm = "legarra",
        ldProjection = NULL
    )$betaEst
    individual <- suppressMessages(
        bayesLWeights(Xs, y, nIter = 3000, nBurnin = 500)
    )
    expect_gt(cor(native, individual), 0.95)
})

test_that("pecotmr's sampler is seeded by R's generator", {
    d <- .rrwPanelStatLd()
    run <- function() {
        bayesARssWeights(d$stat, d$panel, nIter = 400, nBurnin = 100)
    }
    set.seed(5)
    a <- run()
    set.seed(5)
    expect_identical(run(), a)
    expect_length(a, length(d$stat$z))
    expect_true(all(is.finite(a)))
})

test_that("pecotmr's sampler checks the options only it takes", {
    d <- .rrwPanelStatLd()
    expect_error(
        bayesARssWeights(
            d$stat,
            d$panel,
            methodArgs = BayesAlphabetRssOptions(bayesLForm = "legarra")
        ),
        "`bayesLForm` applies to bayesL, not bayesA"
    )
    expect_error(BayesAlphabetRssOptions(pi = 0.1), "the wrapper's `pi`")
    expect_error(
        BayesAlphabetRssOptions(shrinkCorr = 0.9),
        "the wrapper's `shrinkCorr`"
    )
    # qgg's per-variant weights run in place of n for every variant.
    w <- bayesNRssWeights(
        d$stat,
        d$panel,
        nIter = 400,
        nBurnin = 100,
        fitRetention = "slim",
        methodArgs = BayesAlphabetRssOptions(perVariantWeights = TRUE)
    )
    expect_true(all(is.finite(w)))
    expect_equal(attr(w, "fit")$backend, "pecotmr")
    # bayesN has no inclusions to report.
    expect_null(attr(w, "fit")$pip)
    pip <- attr(
        bayesCRssWeights(
            d$stat,
            d$panel,
            nIter = 400,
            nBurnin = 100,
            fitRetention = "slim"
        ),
        "fit"
    )$pip
    expect_true(all(pip >= 0 & pip <= 1))
})

test_that("an sbrm fit leaves the caller's random stream where it was", {
    skip_if_not_installed("hibayes")
    d <- .rrwPanelStatLd()
    set.seed(9)
    w <- bayesNRssWeights(
        d$stat,
        d$ldMatrix,
        backend = "hibayes",
        nIter = 300,
        nBurnin = 100
    )
    after <- runif(1)
    expect_true(all(is.finite(w)))
    set.seed(9)
    invisible(sample.int(.Machine$integer.max, 1L))
    expect_equal(after, runif(1))
})

test_that("pecotmr's sampler stays finite on statistics unlike the LD", {
    # .rrwPanelStat()'s noise ignores the LD, so y'y - 2b'X'y + b'X'Xb, the
    # residual sum of squares, turns negative; qgg's own sampler goes NaN
    # there, while pecotmr's keeps the last residual variance. On the LD
    # unshrunk, the chain then diverges, and says so.
    d <- .rrwPanelStat()
    run <- function(fn, shrinkCorr, ldProjection = NULL) {
        fn(
            d$stat,
            d$panel,
            nIter = 400,
            nBurnin = 100,
            shrinkCorr = shrinkCorr,
            methodArgs = BayesAlphabetRssOptions(ldProjection = ldProjection)
        )
    }
    diverged <- "bayesAlphabetRss: the posterior mean effects' squares"
    w <- expect_warning(run(bayesNRssWeights, 1), diverged)
    expect_true(all(is.finite(w)))
    # hibayes' forms fall back on half the genetic variance, as hibayes does.
    w <- expect_warning(run(bayesBRssWeights, 1), diverged)
    expect_true(all(is.finite(w)))
    # Shrinkage keeps the chain from diverging; so does SDPR's projection,
    # the default, which drops the statistics' part the LD cannot explain.
    expect_silent(run(bayesNRssWeights, 0.8))
    expect_silent(run(bayesNRssWeights, 1, ldProjection = 0.1))
    expect_silent(run(bayesBRssWeights, 1, ldProjection = 0.1))
})

test_that("SDPR's projection reads R (R + aI)^-1 R and R (R + aI)^-1 bhat", {
    d <- .rrwPanelStatLd()
    ldFactor <- pecotmr:::.ldFactorFromRef(d$panel, d$stat$variantNames, "test")
    R <- unname(crossprod(ldFactor) / (nrow(ldFactor) - 1))
    bhat <- pecotmr:::.zToCor(d$stat$z, 500)
    toSpan <- R %*% solve(R + 0.1 * diag(nrow(R)))
    projected <- pecotmr:::.alphabetProjected(bhat, ldFactor, 0.1)
    ldOf <- function(x) crossprod(x$ldFactor) / (nrow(x$ldFactor) - 1)
    expect_equal(ldOf(projected), toSpan %*% R, tolerance = 1e-8)
    expect_equal(projected$bhat, as.vector(toSpan %*% bhat), tolerance = 1e-8)
    # At a = 0 it keeps R, and the statistics' part in R's span.
    hard <- pecotmr:::.alphabetProjected(bhat, ldFactor, 0)
    expect_equal(ldOf(hard), R, tolerance = 1e-8)
    expect_equal(
        hard$bhat,
        as.vector(R %*% MASS::ginv(R) %*% bhat),
        tolerance = 1e-8
    )
    expect_identical(
        pecotmr:::.alphabetProjected(bhat, ldFactor, NULL),
        list(bhat = bhat, ldFactor = ldFactor)
    )
})

test_that("BayesAlphabetRssOptions sets SDPR's projection or turns it off", {
    d <- .rrwPanelStatLd()
    run <- function(methodArgs = NULL) {
        set.seed(7)
        bayesNRssWeights(
            d$stat,
            d$panel,
            nIter = 300,
            nBurnin = 50,
            methodArgs = methodArgs
        )
    }
    ldFactor <- pecotmr:::.ldFactorFromRef(d$panel, d$stat$variantNames, "test")
    direct <- function(ldProjection) {
        set.seed(7)
        pecotmr:::bayesAlphabetRss(
            pecotmr:::.zToCor(d$stat$z, 500),
            ldFactor,
            500,
            "bayesN",
            nIter = 300,
            nBurnin = 50,
            shrinkCorr = 0.95,
            ldProjection = ldProjection
        )$betaEst
    }
    # On by default at SDPR's a = 0.1.
    expect_equal(run(), direct(0.1))
    expect_equal(
        run(BayesAlphabetRssOptions(ldProjection = NULL)),
        direct(NULL)
    )
    # A record whose last setting is NULL, which rlang's splice would drop,
    # reaches the sampler too.
    expect_equal(
        run(BayesAlphabetRssOptions(h2 = 0.5, ldProjection = NULL)),
        direct(NULL)
    )
    expect_equal(run(BayesAlphabetRssOptions(ldProjection = 1)), direct(1))
    expect_false(isTRUE(all.equal(direct(0.1), direct(NULL))))
    expect_error(
        run(BayesAlphabetRssOptions(ldProjection = -1)),
        "ldProjection.*>= 0"
    )
})

test_that("pecotmr's sampler reports progress and refuses bad input", {
    d <- .rrwPanelStatLd()
    ldFactor <- pecotmr:::.ldFactorFromRef(d$panel, d$stat$variantNames, "test")
    bhat <- pecotmr:::.zToCor(d$stat$z, 500)
    out <- capture.output(
        fit <- pecotmr:::bayesAlphabetRss(
            bhat,
            ldFactor,
            500,
            "bayesN",
            nIter = 100,
            nBurnin = 10,
            verbose = TRUE
        )
    )
    expect_match(out, "Iteration +100 of 100", all = FALSE)
    expect_error(
        pecotmr:::bayesAlphabetRss(
            bhat,
            ldFactor,
            500,
            "bayesN",
            perVariantWeights = TRUE,
            nPerVariant = 1:3
        ),
        "one sample size per variant"
    )
    # The entry point checks what the R side already does, for a direct call.
    w <- rep(500, ncol(ldFactor))
    expect_error(
        pecotmr:::bayesAlphabetFactorCpp(
            bhat,
            ldFactor,
            w,
            500L,
            "bayesX",
            FALSE,
            0,
            0.5,
            100L,
            10L,
            1L,
            FALSE,
            1L,
            1,
            TRUE,
            FALSE,
            FALSE
        ),
        "Unknown Bayesian alphabet method 'bayesX'"
    )
    expect_error(
        pecotmr:::bayesAlphabetFactorCpp(
            bhat[-1],
            ldFactor,
            w,
            500L,
            "bayesN",
            FALSE,
            0,
            0.5,
            100L,
            10L,
            1L,
            FALSE,
            1L,
            1,
            TRUE,
            FALSE,
            FALSE
        ),
        "The LD factor has"
    )
})

# ---- pecotmr's sampler: the shared options ----

# The factor of c R + (1 - c) I for R = F'F / (k - 1): F's rows scaled by
# sqrt(c) stacked on the identity's scaled by sqrt(1 - c), all rescaled for
# the k + p - 1 that the sampler divides by.
.baShrunkFactor <- function(ldFactor, c) {
    k <- nrow(ldFactor)
    p <- ncol(ldFactor)
    sqrt(k + p - 1) * rbind(sqrt(c / (k - 1)) * ldFactor, sqrt(1 - c) * diag(p))
}

# One run of pecotmr's sampler on the factor `ldFactor`.
.baFactorRun <- function(
    bhat,
    ldFactor,
    method = "bayesA",
    pi = 0,
    shrinkCorr = 1,
    allowJumpSign = TRUE,
    raoBlackwell = FALSE,
    seed = 3L
) {
    pecotmr:::bayesAlphabetFactorCpp(
        bhat,
        ldFactor,
        rep(500, ncol(ldFactor)),
        500L,
        method,
        FALSE,
        pi,
        0.5,
        300L,
        50L,
        1L,
        FALSE,
        seed,
        shrinkCorr,
        allowJumpSign,
        raoBlackwell,
        FALSE
    )
}

test_that("pecotmr's sampler reads the shrunk LD", {
    # The sampler reads the LD only through its products, so shrinkCorr = c
    # on F runs as the unshrunk sampler on a factor of c R + (1 - c) I, draw
    # for draw up to rounding.
    d <- .rrwPanelStatLd()
    ldFactor <- pecotmr:::.ldFactorFromRef(d$panel, d$stat$variantNames, "test")
    bhat <- pecotmr:::.zToCor(d$stat$z, 500)
    expect_equal(
        .baFactorRun(bhat, ldFactor, shrinkCorr = 0.9),
        .baFactorRun(bhat, .baShrunkFactor(ldFactor, 0.9)),
        tolerance = 1e-8
    )
})

test_that("pecotmr's sampler refuses a sign jump only when asked", {
    # A lone variant with a strong statistic never draws the other sign, so
    # refusing jumps changes nothing there; a null one's draws do jump.
    set.seed(16)
    ldFactor <- matrix(rnorm(50), 50, 1)
    expect_identical(
        .baFactorRun(0.5, ldFactor, allowJumpSign = FALSE),
        .baFactorRun(0.5, ldFactor)
    )
    expect_false(identical(
        .baFactorRun(0, ldFactor, allowJumpSign = FALSE)$betaEst,
        .baFactorRun(0, ldFactor)$betaEst
    ))
})

test_that("pecotmr's Rao-Blackwellized sampler runs the same chain", {
    d <- .rrwPanelStatLd()
    ldFactor <- pecotmr:::.ldFactorFromRef(d$panel, d$stat$variantNames, "test")
    bhat <- pecotmr:::.zToCor(d$stat$z, 500)
    draws <- .baFactorRun(bhat, ldFactor, "bayesC", pi = 0.1)
    means <- .baFactorRun(
        bhat,
        ldFactor,
        "bayesC",
        pi = 0.1,
        raoBlackwell = TRUE
    )
    chain <- c("pip", "pi", "residualVariance")
    expect_identical(means[chain], draws[chain])
    expect_false(identical(means$betaEst, draws$betaEst))
    expect_gt(cor(means$betaEst, draws$betaEst), 0.95)
})

test_that("pecotmr's sampler averages several chains over the kept ones", {
    d <- .rrwPanelStatLd()
    ldFactor <- pecotmr:::.ldFactorFromRef(d$panel, d$stat$variantNames, "test")
    bhat <- pecotmr:::.zToCor(d$stat$z, 500)
    set.seed(8)
    fit <- pecotmr:::bayesAlphabetRss(
        bhat,
        ldFactor,
        500,
        "bayesA",
        nIter = 300,
        nBurnin = 50,
        thin = 1,
        shrinkCorr = 1,
        # Against the kernel run on the statistics as given.
        ldProjection = NULL,
        nChains = 2
    )
    # Each chain's seed is drawn from R's generator.
    set.seed(8)
    seeds <- sample.int(.Machine$integer.max, 2L)
    expect_equal(fit$chains$seed, as.numeric(seeds))
    expect_false(any(fit$chains$diverged))
    kept <- fit$chains$kept
    singles <- map(seeds[kept], function(seed) {
        .baFactorRun(bhat, ldFactor, seed = seed)
    })
    expect_equal(
        fit$betaEst,
        reduce(map(singles, "betaEst"), `+`) / sum(kept)
    )
})

test_that("every engine reads the wrapper's shrunk LD", {
    d <- .rrwPanelStatLd()
    R <- pecotmr:::.ldFromRef(d$panel, d$stat$variantNames, "test")
    shrunk <- LdMatrix(0.9 * R + 0.1 * diag(nrow(R)), nSamples = ncol(d$panel))
    # pecotmr's sampler takes the shrinkage itself, on the factor.
    ld <- pecotmr:::.rssLdSource(d$stat, d$panel, "test")
    bhat <- pecotmr:::.rssMarginalCor(d$stat, n = 500, p = ncol(ld$ldFactor))
    set.seed(4)
    want <- pecotmr:::bayesAlphabetRss(
        bhat,
        ld$ldFactor,
        500,
        "bayesA",
        nIter = 300,
        nBurnin = 50,
        nPerVariant = d$stat$n,
        shrinkCorr = 0.9
    )$betaEst
    set.seed(4)
    expect_equal(
        bayesARssWeights(
            d$stat,
            d$panel,
            nIter = 300,
            nBurnin = 50,
            shrinkCorr = 0.9
        ),
        want
    )
    # qgg and hibayes are given the shrunk matrix.
    run <- function(fn, backend, ldReference, shrinkCorr) {
        set.seed(5)
        suppressMessages(fn(
            d$stat,
            ldReference,
            backend = backend,
            nIter = 300,
            nBurnin = 50,
            shrinkCorr = shrinkCorr
        ))
    }
    skip_if_not_installed("hibayes")
    expect_equal(
        run(bayesBRssWeights, "hibayes", d$ldMatrix, 0.9),
        run(bayesBRssWeights, "hibayes", shrunk, 1)
    )
    skip_if_not_installed("qgg")
    expect_equal(
        run(bayesCRssWeights, "qgg", d$ldMatrix, 0.9),
        run(bayesCRssWeights, "qgg", shrunk, 1)
    )
})

test_that("the alphabet wrapper records and checks its LD shrinkage", {
    d <- .rrwPanelStatLd()
    set.seed(6)
    w <- bayesARssWeights(
        d$stat,
        d$panel,
        nIter = 300,
        nBurnin = 50,
        shrinkCorr = 0.7,
        fitRetention = "slim"
    )
    expect_equal(attr(w, "fit")$shrinkCorr, 0.7)
    expect_error(
        bayesARssWeights(d$stat, d$panel, shrinkCorr = "heldOut"),
        "shrinkCorr"
    )
    expect_error(
        bayesARssWeights(d$stat, d$panel, shrinkCorr = 0),
        "`shrinkCorr` must be above 0"
    )
})

test_that("pecotmr's sampler errors when every one of its chains diverges", {
    # Statistics of 1e-6 leave the draws near the prior's, far past twice
    # the statistics' sum of squares.
    set.seed(17)
    ldFactor <- matrix(rnorm(50 * 3), 50, 3)
    run <- function(n) {
        pecotmr:::bayesAlphabetRss(
            rep(1e-6, 3),
            ldFactor,
            n,
            "bayesA",
            nIter = 60,
            nBurnin = 10,
            nChains = 2
        )
    }
    expect_error(run(100), "bayesAlphabetRss: every chain diverged")
    # The samplers divide by n - 1 and the standard errors by n - 2.
    expect_error(run(2), "bayesAlphabetRss: the sample size must exceed 2")
    stat <- list(z = rep(1, 3), n = rep(2, 3), variantNames = .rrwIds(3))
    expect_error(
        bayesARssWeights(stat, .rrwLdMatrix(diag(3), nSamples = 50L)),
        "bayesAlphabetRssWeights: the sample size must exceed 2"
    )
})

test_that("hibayes' form sends a refused jump to the null class", {
    # Effects of variants in LD trade signs between sweeps; in bayesB a
    # refused jump zeroes the effect and drops it from the slab, which
    # changes the chain.
    d <- .rrwPanelStatLd()
    ldFactor <- pecotmr:::.ldFactorFromRef(d$panel, d$stat$variantNames, "test")
    bhat <- pecotmr:::.zToCor(d$stat$z, 500)
    expect_false(identical(
        .baFactorRun(bhat, ldFactor, "bayesB", pi = 0.2, allowJumpSign = FALSE),
        .baFactorRun(bhat, ldFactor, "bayesB", pi = 0.2)
    ))
})
