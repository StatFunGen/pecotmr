context("regularized_regression - LDpred2")

# ---- pecotmr's sampler against bigsnpr's -----------------------------------

# The panel's factor for pecotmr, the matrix it gives for bigsnpr. Effects in
# tight LD are not identified from the panel's 100 samples, so the fits are
# compared by their fitted correlations R b, averaged over seeds.
.ld2FittedMean <- function(fit, R, seeds) {
    fitted <- map(seeds, function(sd) {
        set.seed(sd)
        as.vector(R %*% fit())
    })
    reduce(fitted, `+`) / length(seeds)
}

test_that("pecotmr's LDpred2-grid sampler agrees with bigsnpr's", {
    skip_if_not_installed("bigsnpr")
    d <- .rrwPanelStatLd()
    R <- pecotmr:::.ldFromRef(d$panel, d$stat$variantNames, "test")
    point <- function(ld, backend) {
        function() {
            ldpred2GridWeights(
                d$stat,
                ld,
                backend = backend,
                p = 0.1,
                h2 = 0.15,
                sparse = FALSE
            )
        }
    }
    ours <- .ld2FittedMean(point(d$panel, "pecotmr"), R, 1:6)
    theirs <- .ld2FittedMean(point(d$ldMatrix, "bigsnpr"), R, 1:6)
    expect_gt(cor(ours, theirs), 0.99)
})

test_that("pecotmr's LDpred2-auto sampler agrees with bigsnpr's", {
    skip_if_not_installed("bigsnpr")
    d <- .rrwPanelStatLd()
    R <- pecotmr:::.ldFromRef(d$panel, d$stat$variantNames, "test")
    run <- function(ld, backend) {
        ldpred2AutoWeights(
            d$stat,
            ld,
            backend = backend,
            h2Init = 0.1,
            pInit = c(0.05, 0.1),
            fitRetention = "slim"
        )
    }
    set.seed(1)
    ours <- run(d$panel, "pecotmr")
    set.seed(1)
    theirs <- run(d$ldMatrix, "bigsnpr")
    expect_gt(cor(as.vector(R %*% ours), as.vector(R %*% theirs)), 0.99)
    expect_equal(attr(ours, "fit")$h2, attr(theirs, "fit")$h2, tolerance = 0.2)
})

# ---- the sampler ------------------------------------------------------------

test_that("an LDpred2 chain stops once its effects outgrow the statistics", {
    # Statistics of 1e-6 on one sample each leave the prior's draws, of
    # variance near 1, far above twice the statistics' sum of squares.
    set.seed(6)
    ldFactor <- scale(matrix(rnorm(40 * 5), 40, 5))
    run <- pecotmr:::ldpred2FactorCpp(
        rep(1e-6, 5),
        ldFactor,
        rep(1, 5),
        FALSE,
        1,
        1,
        FALSE,
        1,
        TRUE,
        1e-5,
        1,
        1,
        100L,
        50L,
        FALSE,
        1L
    )
    expect_true(run$diverged)
    stat <- list(z = rep(1e-6, 5), n = rep(1, 5), variantNames = .rrwIds(5))
    ld <- .rrwLdMatrix(crossprod(ldFactor) / 39, nSamples = 40L)
    skip_if_not_installed("bigsnpr")
    expect_error(
        ldpred2GridWeights(stat, ld, p = 1, h2 = 1, sparse = FALSE),
        "every grid point diverged"
    )
    expect_error(
        ldpred2AutoWeights(stat, ld, h2Init = 1, pInit = 1),
        "every chain diverged"
    )
})

test_that("pecotmr's LDpred2 sampler checks its input and reports progress", {
    set.seed(7)
    ldFactor <- scale(matrix(rnorm(40 * 5), 40, 5))
    run <- function(n, verbose) {
        pecotmr:::ldpred2FactorCpp(
            rep(0.3, 5),
            ldFactor,
            n,
            FALSE,
            0.5,
            0.2,
            FALSE,
            1,
            TRUE,
            1e-5,
            1,
            1,
            100L,
            50L,
            verbose,
            1L
        )
    }
    expect_error(
        run(rep(100, 4), FALSE),
        "There are 4 sample sizes but 5 effects"
    )
    expect_match(
        capture.output(run(rep(100, 5), TRUE)),
        "Iteration +100 of 100",
        all = FALSE
    )
})

test_that("LDpred2-auto keeps the chains the LDpred2 vignette keeps", {
    chain <- function(span, diverged = FALSE) {
        list(corrEst = c(0, span), diverged = diverged)
    }
    chains <- list(chain(1), chain(0.99), chain(0.5), chain(NA, TRUE))
    # The threshold is 0.95 times the spans' 95th percentile, about 0.95.
    expect_equal(
        pecotmr:::.mcmcKeptChains(chains, "test"),
        c(TRUE, TRUE, FALSE, FALSE)
    )
    expect_error(
        pecotmr:::.mcmcKeptChains(list(chain(NA, TRUE)), "test"),
        "test: every chain diverged"
    )
})

# ---- the weight functions --------------------------------------------------

test_that("ldpred2AutoWeights reports its chains at each retention", {
    d <- .rrwPanelStatLd()
    set.seed(2)
    w <- ldpred2AutoWeights(
        d$stat,
        d$panel,
        h2Init = 0.1,
        pInit = c(0.01, 0.1, 0.2),
        fitRetention = "full"
    )
    fit <- attr(w, "fit")
    expect_length(w, 40L)
    expect_named(fit$chains, c("pInit", "p", "h2", "diverged", "kept"))
    expect_equal(fit$chains$pInit, c(0.01, 0.1, 0.2))
    expect_length(fit$postp, 40L)
    expect_equal(dim(fit$beta), c(40L, 3L))
    expect_equal(
        as.numeric(w),
        rowMeans(fit$beta[, fit$chains$kept, drop = FALSE])
    )
    set.seed(2)
    expect_null(attr(
        ldpred2AutoWeights(d$stat, d$panel, h2Init = 0.1, pInit = 0.1),
        "fit"
    ))
})

test_that("ldpred2GridWeights records its pick and its grid", {
    d <- .rrwPanelStatLd()
    set.seed(3)
    w <- ldpred2GridWeights(
        d$stat,
        d$panel,
        p = c(0.01, 0.1, 1),
        h2 = c(0.1, 0.2),
        fitRetention = "full"
    )
    sel <- attr(w, "ldpred2_selection")
    expect_equal(unname(sel[["mode"]]), "heldOut")
    path <- attr(w, "fit")$path
    expect_named(path, c("p", "h2", "sparse", "diverged", "selected"))
    expect_equal(nrow(path), 12L)
    # bigsnpr's order: large p first.
    expect_equal(path$p[[1L]], 1)
    expect_equal(sum(path$selected), 1L)
    expect_equal(as.numeric(sel[["p"]]), path$p[path$selected])
    expect_equal(as.numeric(w), attr(w, "fit")$beta[, path$selected])
    skip_if_not_installed("fdrtool")
    q <- suppressWarnings(ldpred2GridWeights(
        d$stat,
        d$panel,
        p = c(0.1, 1),
        h2 = 0.2,
        selection = "ldQuadratic"
    ))
    expect_equal(unname(attr(q, "ldpred2_selection")[["mode"]]), "ldQuadratic")
})

test_that("LDpred2's heritability defaults to the window's LAVA estimate", {
    d <- .rrwPanelStatLd()
    gridH2 <- function() {
        set.seed(4)
        w <- ldpred2GridWeights(d$stat, d$panel, p = 0.1, fitRetention = "slim")
        sort(unique(attr(w, "fit")$path$h2))
    }
    local_mocked_bindings(.cisH2Lava = function(...) list(h2 = 0.2))
    expect_equal(gridH2(), 0.2 * c(0.3, 0.7, 1, 1.4))
    # Where LAVA gives no estimate, or a negative one, the grid scales
    # bigsnpr's floor.
    for (h2 in c(NA_real_, -0.05)) {
        local_mocked_bindings(.cisH2Lava = function(...) list(h2 = h2))
        expect_equal(gridH2(), 1e-3 * c(0.3, 0.7, 1, 1.4))
    }
    # LAVA reads the reference's sample size.
    unsized <- LdMatrix(
        pecotmr:::.ldFromRef(d$panel, d$stat$variantNames, "test")
    )
    expect_error(
        ldpred2AutoWeights(d$stat, unsized),
        "estimating the heritability by LAVA needs the LD reference's"
    )
})

test_that("the LDpred2 weights check their engine and settings", {
    d <- .rrwPanelStatLd()
    expect_error(
        ldpred2AutoWeights(d$stat, d$ldMatrix, backend = "pecotmr"),
        "runs on a genotype reference's factor"
    )
    expect_error(
        ldpred2GridWeights(
            d$stat,
            d$panel,
            backend = "bigsnpr",
            allowJumpSign = FALSE
        ),
        "always allows a drawn effect to change sign"
    )
    expect_error(
        ldpred2AutoWeights(d$stat, d$panel, h2Init = 0.1, pInit = c(0, 0.1)),
        "`pInit` must be above 0"
    )
    expect_error(
        ldpred2AutoWeights(d$stat, d$panel, h2Init = 0.1, shrinkCorr = 0),
        "`shrinkCorr` must be above 0"
    )
    expect_error(
        ldpred2GridWeights(d$stat, d$panel, h2 = 0.1, nIter = 50),
        "must exceed `nBurnin`"
    )
})

test_that("bigsnpr's grid sampler reads the LD shrunk to shrinkCorr", {
    # bigsnpr's grid has no shrinkage of its own, so the matrix written for it
    # is shrunk; LDpred2-auto's is left to bigsnpr's own shrink_corr.
    skip_if_not_installed("bigsnpr")
    d <- .rrwPanelStatLd()
    inp <- pecotmr:::.rssPerVariantInput(d$stat, d$ldMatrix, "test")
    written <- NULL
    local_mocked_bindings(
        as_SFBM = function(spmat, ...) {
            written <<- as.matrix(spmat)
            NULL
        },
        .package = "bigsnpr"
    )
    R <- unname(inp$ld$R)
    pecotmr:::.ldpred2Engine(inp, "bigsnpr", 0.9)
    expect_equal(unname(written), 0.9 * R + 0.1 * diag(nrow(R)))
    pecotmr:::.ldpred2Engine(inp, "bigsnpr")
    expect_equal(unname(written), R)
    expect_equal(formals(ldpred2GridWeights)$shrinkCorr, 0.95)
})

test_that("bigsnpr's LD file lives only as long as the call", {
    skip_if_not_installed("bigsnpr")
    d <- .rrwPanelStatLd()
    before <- list.files(tempdir(), pattern = "[.]sbk$")
    set.seed(5)
    ldpred2GridWeights(d$stat, d$ldMatrix, p = c(0.1, 1), h2 = 0.2)
    expect_setequal(list.files(tempdir(), pattern = "[.]sbk$"), before)
})
