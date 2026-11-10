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

test_that("a summary-statistics fit picks its engine by method", {
    f <- pecotmr:::.methodBackend
    # qgg, or hibayes for bayesB and bayesL, on either kind of LD reference.
    expect_equal(f("bayesA", NULL, "sumstat"), "qgg")
    expect_equal(f("bayesB", NULL, "sumstat"), "hibayes")
    expect_equal(f("bayesL", NULL, "sumstat"), "hibayes")
    expect_equal(f("bayesA", "hibayes", "sumstat"), "hibayes")
    expect_error(
        f("bayesA", "pecotmr", "sumstat"),
        "must be one of qgg, hibayes"
    )
    expect_error(f("bayesB", "qgg", "sumstat"), "must be one of hibayes")
    # bayesB has hibayes alone there, which leaves nothing to choose; LDpred2
    # has bigsnpr alone; the lasso has two.
    has <- pecotmr:::.methodHasBackends
    expect_false(has("bayesB", "sumstat"))
    expect_true(has("bayesB", "individual"))
    expect_false(has("ldpred2Auto", "sumstat"))
    expect_true(has("lasso", "sumstat"))
    expect_equal(f("lasso", NULL, "sumstat"), "pecotmr")
    expect_equal(f("lasso", "bigsnpr", "sumstat"), "bigsnpr")
})

test_that("qgg and hibayes read a genotype reference as its LD matrix", {
    skip_if_not_installed("qgg")
    skip_if_not_installed("hibayes")
    d <- .rrwPanelStatLd()
    run <- function(fn, ldReference) {
        set.seed(2)
        suppressMessages(fn(
            d$stat,
            ldReference,
            nIter = 300,
            nBurnin = 50,
            ldProjection = NULL
        ))
    }
    # The matrix formed from the panel is the LdMatrix's, so as given one
    # seed gives one fit from either.
    expect_equal(
        run(bayesCRssWeights, d$panel),
        run(bayesCRssWeights, d$ldMatrix)
    )
    expect_equal(
        run(bayesBRssWeights, d$panel),
        run(bayesBRssWeights, d$ldMatrix)
    )
})

# The input the wrapper gives `engine` for .rrwPanelStatLd()'s statistics,
# at the default shrinkage and, unless given, the default projection.
.baDense <- function(d, engine, ldReference = d$ldMatrix, ldProjection = 0.1) {
    pecotmr:::.alphabetRssDense(
        d$stat,
        ldReference,
        list(shrinkCorr = 0.95, ldProjection = ldProjection),
        engine
    )
}

# SDPR's likelihood at a = 0.1, shrunk at c = 0.95, from the eigenvectors of
# the panel's LD: the LD R (R + aI)^-1 R and the effects R (R + aI)^-1 bhat.
.baProjected <- function(d) {
    R <- unname(pecotmr:::.ldFromRef(d$panel, d$stat$variantNames, "test"))
    bhat <- pecotmr:::.rssMarginalCor(d$stat, n = 500, p = nrow(R))
    ev <- eigen(R, symmetric = TRUE)
    kept <- ev$values > sqrt(.Machine$double.eps) * max(ev$values)
    lambda <- ev$values[kept]
    V <- ev$vectors[, kept]
    projected <- V %*% (lambda^2 / (lambda + 0.1) * t(V))
    list(
        R = 0.95 * projected + 0.05 * diag(nrow(R)),
        b = as.vector(V %*% (lambda / (lambda + 0.1) * crossprod(V, bhat))),
        bhat = bhat
    )
}

test_that("the projection gives each engine SDPR's likelihood", {
    d <- .rrwPanelStatLd()
    want <- .baProjected(d)
    n <- 500
    # From the panel's factor or from its LD matrix alike.
    for (ref in list(d$panel, d$ldMatrix)) {
        # qgg reads X'X = sqrt(ww ww') LD, X'y = wy and
        # y'y = median((b^2 + (n - 2) seb^2) ww).
        qgg <- .baDense(d, "qgg", ref)
        expect_equal(
            unname(sqrt(outer(qgg$ww, qgg$ww)) * qgg$R),
            n * want$R,
            tolerance = 1e-8
        )
        expect_equal(qgg$wy, n * want$b, tolerance = 1e-8)
        expect_equal(median((qgg$b^2 + (n - 2) * qgg$seb^2) * qgg$ww), n)
        expect_equal(
            pecotmr:::.alphabetFittedCorr(qgg, want$b),
            as.vector(want$R %*% want$b),
            tolerance = 1e-8
        )
        # hibayes reads X'X = n ldm, X'y = n ldm_jj b and
        # y'y = mean(n ldm_jj (b^2 + (n - 2) se^2)).
        hib <- .baDense(d, "hibayes", ref)
        expect_equal(unname(hib$R), want$R, tolerance = 1e-8)
        expect_equal(diag(hib$R) * hib$b, want$b, tolerance = 1e-8)
        expect_equal(
            mean(n * diag(hib$R) * (hib$b^2 + (n - 2) * hib$seb^2)),
            n
        )
        # A chain's divergence is tested against the effects as given.
        expect_equal(qgg$bhat, want$bhat)
        expect_equal(hib$bhat, want$bhat)
    }
    # Off, both engines read the statistics as given.
    given <- .baDense(d, "qgg", ldProjection = NULL)
    R <- unname(pecotmr:::.ldFromRef(d$panel, d$stat$variantNames, "test"))
    expect_equal(unname(given$R), 0.95 * R + 0.05 * diag(nrow(R)))
    expect_equal(given$b, want$bhat)
    expect_null(given$ww)
    expect_null(given$scale)
})

test_that("hibayes reads the projected likelihood exactly", {
    skip_if_not_installed("hibayes")
    d <- .rrwPanelStatLd()
    hib <- .baDense(d, "hibayes")
    # sbrm()'s conjugate gradient solves ldm g = X'y / n, here R g = b'.
    fit <- withr::with_preserve_seed(hibayes::sbrm(
        sumstat = data.frame(
            SNP = hib$ids,
            A1 = "A",
            A2 = "G",
            freq = 0.5,
            b = hib$b,
            se = hib$seb,
            p = 2 * pnorm(-abs(hib$b / hib$seb)),
            N = hib$n
        ),
        ldm = unname(hib$R),
        method = "CG",
        verbose = FALSE
    ))
    expect_equal(
        as.numeric(fit$g),
        as.vector(solve(hib$R, diag(hib$R) * hib$b)),
        tolerance = 1e-4
    )
})

test_that("qgg reads the projected likelihood exactly", {
    skip_if_not_installed("qgg")
    d <- .rrwPanelStatLd()
    qgg <- .baDense(d, "qgg")
    qin <- pecotmr:::.alphabetQggInputs(qgg)
    n <- qgg$n
    h2 <- 0.01
    # blup with its variances held solves (X'X + (ve / vb) I) b = X'y. qgg
    # 1.1.6 hands its kernel the genetic variance y'y h2 / (n - 1) in vb's
    # place; giving vb that same value makes ve / vb = (1 - h2) / h2 = 99
    # whichever way it is passed.
    fit <- suppressMessages(qgg::gbayes(
        stat = qin$stat,
        Glist = qin$Glist,
        LD = qin$LD,
        method = "blup",
        h2 = h2,
        vb = n / (n - 1) * h2,
        updateB = FALSE,
        updateE = FALSE,
        updateG = FALSE,
        nit = 100,
        nburn = 1000,
        nthin = 1
    ))
    XtX <- sqrt(outer(qgg$ww, qgg$ww)) * qgg$R
    expect_equal(
        unname(fit$stat[qgg$ids, "bm"]),
        as.vector(solve(XtX + 99 * diag(nrow(XtX)), qgg$wy)),
        tolerance = 1e-6
    )
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

# The first draw of the stream `seed` starts.
.baSeedDraw <- function(seed) {
    withr::with_seed(seed, runif(1))
}

# Effects that are the marginal ones scaled by the first draw of the stream
# a chain's seed starts: each chain's fitted correlations are that draw
# times the same vector, so LDpred2's filter on their spans keeps the chains
# whose draws are near the largest. `runaway` scales the draws below it by
# 10, past the divergence bound.
.baChainFitter <- function(runaway = 0) {
    function(solverInput, ...) {
        u <- runif(1)
        list(betaEst = solverInput$bhat * if (u < runaway) 10 else u)
    }
}

test_that("each chain is the engine run from its seed, tested on its mean", {
    d <- .rrwPanelStatLd()
    inp <- .baDense(d, "hibayes")
    local_mocked_bindings(.alphabetHibayesRssFit = .baChainFitter())
    chain <- function(track) {
        pecotmr:::.alphabetRssChain(
            ".alphabetHibayesRssFit",
            inp,
            method = "bayesN",
            seed = 11L,
            track = track
        )
    }
    set.seed(1)
    plain <- chain(FALSE)
    after <- runif(1)
    set.seed(1)
    expect_equal(after, runif(1))
    expect_named(plain, "betaEst")
    expect_equal(plain$betaEst, inp$bhat * .baSeedDraw(11L))
    # Tracked, it carries the fitted correlations and the divergence test of
    # its posterior mean.
    tracked <- chain(TRUE)
    expect_equal(tracked$betaEst, plain$betaEst)
    expect_equal(tracked$corrEst, as.vector(inp$R %*% plain$betaEst))
    expect_false(tracked$diverged)
    local_mocked_bindings(.alphabetHibayesRssFit = .baChainFitter(1))
    expect_true(chain(TRUE)$diverged)
})

test_that("several chains average those LDpred2's filter keeps", {
    d <- .rrwPanelStatLd()
    inp <- .baDense(d, "hibayes")
    local_mocked_bindings(.alphabetHibayesRssFit = .baChainFitter())
    run <- function() {
        set.seed(4)
        bayesBRssWeights(
            d$stat,
            d$ldMatrix,
            nIter = 20,
            nBurnin = 5,
            nChains = 4,
            fitRetention = "slim"
        )
    }
    w <- run()
    fit <- attr(w, "fit")
    expect_equal(nrow(fit$chains), 4L)
    expect_false(any(fit$chains$diverged))
    u <- vapply(fit$chains$seed, .baSeedDraw, numeric(1))
    kept <- u > 0.95 * quantile(u, 0.95)
    expect_equal(fit$chains$kept, kept)
    expect_equal(as.numeric(w), mean(u[kept]) * inp$bhat)
    expect_equal(fit$backend, "hibayes")
    expect_equal(w, run())
})

test_that("a diverged chain is dropped, warned of alone, fatal when all are", {
    d <- .rrwPanelStatLd()
    local_mocked_bindings(.alphabetHibayesRssFit = .baChainFitter(0.5))
    set.seed(21)
    w <- bayesBRssWeights(
        d$stat,
        d$ldMatrix,
        nIter = 20,
        nBurnin = 5,
        nChains = 6,
        fitRetention = "slim"
    )
    chains <- attr(w, "fit")$chains
    u <- vapply(chains$seed, .baSeedDraw, numeric(1))
    expect_equal(chains$diverged, u < 0.5)
    expect_true(any(chains$diverged) && !all(chains$diverged))
    expect_false(any(chains$kept & chains$diverged))
    # Alone, a diverged chain is returned with a warning.
    local_mocked_bindings(.alphabetHibayesRssFit = .baChainFitter(1))
    expect_warning(
        bayesBRssWeights(d$stat, d$ldMatrix, nIter = 20, nBurnin = 5),
        "bayesAlphabetRssWeights: the posterior mean effects' squares"
    )
    expect_error(
        bayesBRssWeights(
            d$stat,
            d$ldMatrix,
            nIter = 20,
            nBurnin = 5,
            nChains = 2
        ),
        "bayesAlphabetRssWeights: every chain diverged"
    )
})

test_that("qgg and hibayes read the wrapper's shrunk LD", {
    d <- .rrwPanelStatLd()
    R <- pecotmr:::.ldFromRef(d$panel, d$stat$variantNames, "test")
    shrunk <- LdMatrix(0.9 * R + 0.1 * diag(nrow(R)), nSamples = ncol(d$panel))
    run <- function(fn, ldReference, shrinkCorr) {
        set.seed(5)
        suppressMessages(fn(
            d$stat,
            ldReference,
            nIter = 300,
            nBurnin = 50,
            shrinkCorr = shrinkCorr,
            ldProjection = NULL
        ))
    }
    skip_if_not_installed("hibayes")
    expect_equal(
        run(bayesBRssWeights, d$ldMatrix, 0.9),
        run(bayesBRssWeights, shrunk, 1)
    )
    skip_if_not_installed("qgg")
    expect_equal(
        run(bayesCRssWeights, d$ldMatrix, 0.9),
        run(bayesCRssWeights, shrunk, 1)
    )
})

test_that("the alphabet wrapper records and checks its settings", {
    skip_if_not_installed("hibayes")
    d <- .rrwPanelStatLd()
    set.seed(6)
    w <- bayesARssWeights(
        d$stat,
        d$panel,
        backend = "hibayes",
        nIter = 300,
        nBurnin = 50,
        shrinkCorr = 0.7,
        fitRetention = "slim"
    )
    expect_equal(attr(w, "fit")$shrinkCorr, 0.7)
    expect_equal(attr(w, "fit")$ldProjection, 0.1)
    expect_equal(attr(w, "fit")$backend, "hibayes")
    expect_null(attr(w, "fit")$chains)
    expect_error(
        bayesARssWeights(d$stat, d$panel, shrinkCorr = "heldOut"),
        "shrinkCorr"
    )
    expect_error(
        bayesARssWeights(d$stat, d$panel, shrinkCorr = 0),
        "`shrinkCorr` must be above 0"
    )
    expect_error(bayesARssWeights(d$stat, d$panel, nChains = 0), "nChains")
    expect_error(
        bayesARssWeights(d$stat, d$panel, ldProjection = -1),
        "ldProjection"
    )
    expect_error(
        bayesARssWeights(d$stat, d$panel, ldProjection = "none"),
        "ldProjection"
    )
    # The engines' standard errors divide by n - 2.
    stat <- list(z = rep(1, 3), n = rep(2, 3), variantNames = .rrwIds(3))
    expect_error(
        bayesARssWeights(stat, .rrwLdMatrix(diag(3), nSamples = 50L)),
        "bayesAlphabetRssWeights: the sample size must exceed 2"
    )
})

# ---- h2 ----
test_that("h2 reaches qgg once, from the wrapper", {
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
    invisible(bayesNWeights(d$X, d$y, h2 = 0.2, nIter = 20, nBurnin = 5))
    expect_equal(seen$h2, 0.2)
    # NULL leaves the engine its own default, so none is sent.
    invisible(bayesNWeights(d$X, d$y, nIter = 20, nBurnin = 5))
    expect_null(seen$h2)
    # Alone, the engine's own spelling still reaches it...
    invisible(bayesNWeights(
        d$X,
        d$y,
        nIter = 20,
        nBurnin = 5,
        methodArgs = QggOptions(h2 = 0.3)
    ))
    expect_equal(seen$h2, 0.3)
    # ...but both together are the same setting twice.
    expect_error(
        bayesNWeights(d$X, d$y, h2 = 0.2, methodArgs = QggOptions(h2 = 0.3)),
        "given twice: h2"
    )
})

test_that("each backend reads h2 as its engine spells it", {
    skip_if_not_installed("BGLR")
    skip_if_not_installed("hibayes")
    withH2 <- pecotmr:::.alphabetWithH2
    expect_equal(
        as.list(withH2(BglrOptions(), "BglrOptions", "bglr", 0.2)),
        list(R2 = 0.2)
    )
    # hibayes fixes its own at 0.5, so it gets the variances its defaults
    # would give at 0.2: (dfvg - 2) / dfvg vary h2 and vary (1 - h2).
    hib <- as.list(
        withH2(HibayesOptions(), "HibayesOptions", "hibayes", 0.2, 2)
    )
    expect_equal(hib[c("vg", "ve")], list(vg = 0.5 * 2 * 0.2, ve = 2 * 0.8))
    six <- as.list(withH2(
        HibayesOptions(dfvg = 6),
        "HibayesOptions",
        "hibayes",
        0.2,
        2
    ))
    expect_equal(six[["vg"]], 4 / 6 * 2 * 0.2)
    expect_identical(
        withH2(QggOptions(), "QggOptions", "qgg", NULL),
        QggOptions()
    )
})

test_that("h2 = 0.5, every engine's own default, changes no fit", {
    skip_if_not_installed("qgg")
    skip_if_not_installed("hibayes")
    skip_if_not_installed("BGLR")
    d <- .baXy()
    run <- function(method, backend, h2 = NULL) {
        set.seed(1)
        suppressMessages(bayesAlphabetWeights(
            d$X,
            d$y,
            method,
            backend = backend,
            h2 = h2,
            nIter = 60,
            nBurnin = 10
        ))
    }
    engines <- list(
        c("bayesN", "qgg"),
        c("bayesB", "bglr"),
        c("bayesN", "hibayes")
    )
    for (b in engines) {
        expect_equal(
            run(b[[1L]], b[[2L]], 0.5),
            run(b[[1L]], b[[2L]]),
            tolerance = 1e-8,
            info = b[[2L]]
        )
    }
    s <- .rrwPanelStatLd()
    runRss <- function(backend, ld, h2 = NULL) {
        set.seed(1)
        suppressMessages(bayesAlphabetRssWeights(
            s$stat,
            ld,
            "bayesN",
            backend = backend,
            h2 = h2,
            nIter = 60,
            nBurnin = 10
        ))
    }
    for (b in list(
        list("qgg", s$ldMatrix),
        list("hibayes", s$panel)
    )) {
        expect_equal(
            runRss(b[[1L]], b[[2L]], 0.5),
            runRss(b[[1L]], b[[2L]]),
            tolerance = 1e-8,
            info = b[[1L]]
        )
    }
})

test_that("h2 can be estimated first, on the data the fit reads", {
    skip_if_not_installed("qgg")
    d <- .baXy(n = 120, p = 30)
    seen <- NULL
    local_mocked_bindings(
        gbayes = function(...) {
            seen <<- list(...)
            list(bm = rep(0, 30))
        },
        .package = "qgg"
    )
    invisible(bayesNWeights(d$X, d$y, h2 = "lava", nIter = 20, nBurnin = 5))
    want <- pecotmr:::.h2IndividualFit("lava", d$X, d$y)$h2
    expect_equal(seen$h2, min(max(want, 0.001), 0.99))
    expect_error(bayesNWeights(d$X, d$y, h2 = "lavaa"), "must be one of")
    s <- .rrwPanelStatLd()
    got <- suppressWarnings(
        pecotmr:::.alphabetRssH2("lava", s$stat, s$panel)
    )
    est <- suppressWarnings(
        pecotmr:::.h2StatFit("lava", s$stat, s$panel, "test")$h2
    )
    expect_equal(got, min(max(est, 0.001), 0.99))
    # GREML needs genotypes.
    expect_error(
        bayesNRssWeights(s$stat, s$panel, h2 = "greml"),
        "must be one of"
    )
})
