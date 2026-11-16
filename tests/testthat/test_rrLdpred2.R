context("regularized_regression - LDpred2")

# ---- divergence and the chain filter ----------------------------------------

test_that("an LDpred2 run stops once its effects outgrow the statistics", {
    # Statistics of 1e-6 on one sample each leave the prior's draws, of
    # variance near 1, far above twice the statistics' sum of squares.
    set.seed(6)
    ldFactor <- scale(matrix(rnorm(40 * 5), 40, 5))
    stat <- list(cor = rep(1e-6, 5), n = rep(1, 5), variantNames = .rrwIds(5))
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
    skip_if_not_installed("bigsnpr")
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
    skip_if_not_installed("bigsnpr")
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
    skip_if_not_installed("bigsnpr")
    d <- .rrwPanelStatLd()
    gridH2 <- function() {
        set.seed(4)
        w <- ldpred2GridWeights(d$stat, d$panel, p = 0.1, fitRetention = "slim")
        sort(unique(attr(w, "fit")$path$h2))
    }
    local_mocked_bindings(.lavaFit = function(...) list(h2 = 0.2))
    expect_equal(gridH2(), 0.2 * c(0.3, 0.7, 1, 1.4))
    # Where LAVA gives no estimate, or a negative one, the grid scales
    # bigsnpr's floor.
    for (h2 in c(NA_real_, -0.05)) {
        local_mocked_bindings(.lavaFit = function(...) list(h2 = h2))
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

test_that("the LDpred2 weights check their settings", {
    d <- .rrwPanelStatLd()
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
    pecotmr:::.bigsnprEngine(inp, 0.9)
    expect_equal(unname(written), 0.9 * R + 0.1 * diag(nrow(R)))
    pecotmr:::.bigsnprEngine(inp)
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

test_that("LDpred2's heritability can come from any per-trait estimator", {
    d <- .rrwPanelStatLd()
    inp <- pecotmr:::.rssPerVariantInput(d$stat, d$panel, "test")
    hdl <- pecotmr:::.cisH2HdlL(
        pecotmr:::.corToZ(inp$bhat, inp$n),
        inp$n,
        inp$ld
    )$h2
    expect_equal(
        pecotmr:::.ldpred2H2("hdlL", inp, "test"),
        if (is.na(hdl)) 1e-3 else max(hdl, 1e-3)
    )
    expect_equal(
        pecotmr:::.ldpred2GridH2("lava", inp, "test"),
        pecotmr:::.ldpred2GridH2(NULL, inp, "test")
    )
    expect_error(pecotmr:::.ldpred2H2("greml", inp, "test"), "must be one of")
})
