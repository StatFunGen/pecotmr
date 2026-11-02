context("regularized_regression - lassosumRss")

# ---- lassosumRss ----
test_that("lassosumRss errors on non-matrix R input", {
    expect_error(
        lassosumRss(bhat = rnorm(5), R = "not_a_matrix", n = 100),
        "as a matrix"
    )
})

test_that("lassosumRss errors on non-positive sample size", {
    expect_error(
        lassosumRss(bhat = rnorm(5), R = diag(5), n = -1),
        "valid sample size"
    )
})

test_that("lassosumRss errors on mismatched bhat and R dimensions", {
    expect_error(
        lassosumRss(bhat = rnorm(10), R = diag(5), n = 100),
        "bhat.*Must have length"
    )
})

test_that("lassosumRss runs successfully with valid input", {
    set.seed(42)
    p <- 10
    n <- 100
    bhat <- rnorm(p, sd = 0.1)
    R <- diag(p)
    for (i in seq_len(p - 1)) {
        R[i, i + 1] <- 0.3
        R[i + 1, i] <- 0.3
    }
    result <- lassosumRss(bhat = bhat, R = R, n = n)
    expect_type(result, "list")
    expect_true("betaEst" %in% names(result))
    expect_equal(length(result$betaEst), p)
    expect_true(all(is.finite(result$betaEst)))
    expect_equal(nrow(result$beta), p)
    expect_equal(ncol(result$beta), 20)
    expect_length(result$conv, length(result$lambda))
})

test_that("lassosumRss with large lambda gives all-zero betas", {
    set.seed(42)
    p <- 10
    n <- 100
    bhat <- rnorm(p, sd = 0.1)
    R <- diag(p)
    result <- lassosumRss(bhat = bhat, R = R, n = n, lambda = c(100))
    expect_true(all(result$betaEst == 0))
})

# ---- lassosumRssWeights (wrapper) ----
test_that("lassosumRssWeights calls lassosumRss and returns betaEst", {
    set.seed(42)
    p <- 10
    n <- 100
    bhat <- rnorm(p, sd = 0.1)
    R <- diag(p)
    for (i in seq_len(p - 1)) {
        R[i, i + 1] <- 0.3
        R[i + 1, i] <- 0.3
    }
    stat <- list(b = bhat, n = rep(n, p), variantNames = .rrwIds(p))
    expected <- seq_len(p) * 0.02
    local_mocked_bindings(
        .penalizedRssPath = function(bhat, ld, n, ...) {
            list(
                beta = cbind(rep(0, length(bhat)), expected),
                lambda = c(0.05, 0.01),
                fbeta = c(1, 0.1)
            )
        }
    )
    # The all-zero candidate scores -Inf on every held-out split.
    result <- lassosumRssWeights(
        stat = stat,
        ldReference = .rrwLdMatrix(R, nSamples = n),
        s = 0.5
    )
    expect_equal(length(result), p)
    expect_true(is.numeric(result))
    expect_equal(c(result), expected)
    expect_equal(
        unname(attr(result, "lassosum_selection")["mode"]),
        "heldOut"
    )
})

test_that("lassosumRssWeights clamps correlation input before scaling", {
    set.seed(42)
    p <- 5
    n <- 100
    # Construct bhat with max abs >= 1 to trigger the clamp branch
    bhat <- c(1.5, -0.2, 0.1, 0.05, -0.3)
    R <- diag(p)
    stat <- list(b = bhat, n = rep(n, p), variantNames = .rrwIds(p))
    captured <- NULL
    local_mocked_bindings(
        .penalizedRssPath = function(bhat, ld, n, ...) {
            captured <<- bhat
            list(
                beta = matrix(0, nrow = length(bhat), ncol = 2),
                lambda = c(0.05, 0.01),
                fbeta = c(1, 0.5)
            )
        }
    )
    result <- lassosumRssWeights(
        stat = stat,
        ldReference = .rrwLdMatrix(R),
        s = 0.5,
        selection = "minFbeta"
    )
    scaled_cor <- captured / sqrt(n)
    expect_true(max(abs(scaled_cor)) < 1)
    expect_equal(max(abs(scaled_cor)), 0.9999, tolerance = 1e-9)
    # Sign-preserving rescale
    expect_true(all(sign(scaled_cor) == sign(bhat)))
})

# ---- lassosum on bigsnpr's lassosum2 ----------------------------------------
# snp_lassosum2() times each grid point with system.time(), which runs a full
# garbage collection of the session first, so these tests keep their grids
# to a point or two: held-out selection runs lassosum2 five times (the fit
# and 4 refits).

test_that("lassosum2 fits lassosum's objective at s = delta / (1 + delta)", {
    skip_if_not_installed("bigsnpr")
    d <- .rrwPanelStatLd()
    R <- pecotmr:::.ldFromRef(d$panel, d$stat$variantNames, "test")
    bhat <- d$stat$z / sqrt(500)
    inp <- list(ld = list(R = R), bhat = bhat, n = rep(500, length(bhat)))
    engine <- pecotmr:::.bigsnprEngine(inp)
    on.exit(pecotmr:::.bigsnprRelease(engine))
    fits <- pecotmr:::.lassosum2Run(
        inp,
        engine,
        list(
            delta = c(0.01, 1),
            options = list(nlambda = 2, tol = 1e-10, maxiter = 1e5)
        )
    )
    # At the same lambda, lassosum's fit is lassosum2's times 1 + delta.
    for (j in seq_len(ncol(fits$beta))) {
        ours <- pecotmr:::.penalizedRssPath(
            bhat * sqrt(500),
            list(R = R, shrink = fits$path$s[[j]]),
            500,
            pecotmr:::.rssPenalty("lasso"),
            lambda = fits$path$lambda[[j]],
            thr = 1e-12,
            maxiter = 1e6
        )$beta[, 1]
        expect_equal(ours, fits$beta[, j], tolerance = 1e-6)
    }
})

test_that("lassosumRssWeights runs lassosum2 in lassosum's s", {
    skip_if_not_installed("bigsnpr")
    d <- .rrwPanelStatLd()
    # Without `s`, lassosum2's own delta grid.
    expect_equal(pecotmr:::.lassosum2Delta(NULL), c(0.001, 0.01, 0.1, 1))
    set.seed(1)
    w <- lassosumRssWeights(
        d$stat,
        d$ldMatrix,
        backend = "bigsnpr",
        s = c(0.1, 0.5),
        methodArgs = BigsnprOptions(nlambda = 1),
        fitRetention = "full"
    )
    path <- attr(w, "fit")$path
    expect_named(
        path,
        c("s", "lambda", "iterations", "nparams", "selected")
    )
    # delta = s / (1 - s), so the path reports the s it was given.
    expect_equal(path$s, c(0.1, 0.5))
    sel <- attr(w, "lassosum_selection")
    expect_equal(unname(sel[["mode"]]), "heldOut")
    expect_equal(as.numeric(w), attr(w, "fit")$beta[, path$selected])
    # A genotype reference runs on the matrix formed from it.
    set.seed(1)
    expect_length(
        lassosumRssWeights(
            d$stat,
            d$panel,
            backend = "bigsnpr",
            s = 0.5,
            methodArgs = BigsnprOptions(nlambda = 1)
        ),
        40L
    )
})

test_that("lassosum2 takes lassosum's own pseudovalidation", {
    skip_if_not_installed("bigsnpr")
    skip_if_not_installed("fdrtool")
    d <- .rrwPanelStatLd()
    # fdrtool warns that 40 statistics are too few for reliable FDRs.
    w <- suppressWarnings(lassosumRssWeights(
        d$stat,
        d$ldMatrix,
        backend = "bigsnpr",
        s = c(0.1, 0.5),
        selection = "ldQuadratic",
        methodArgs = BigsnprOptions(nlambda = 1)
    ))
    expect_equal(
        unname(attr(w, "lassosum_selection")[["mode"]]),
        "ldQuadratic"
    )
})

test_that("lassosumRssWeights on bigsnpr refuses what lassosum2 cannot do", {
    skip_if_not_installed("bigsnpr")
    d <- .rrwPanelStatLd()
    expect_error(
        lassosumRssWeights(
            d$stat,
            d$ldMatrix,
            backend = "bigsnpr",
            selection = "minFbeta"
        ),
        "lassosum2 reports no objective"
    )
    expect_error(
        lassosumRssWeights(d$stat, d$ldMatrix, backend = "bigsnpr", s = 1),
        "strictly between 0 and 1"
    )
    expect_error(
        lassosumRssWeights(
            d$stat,
            d$ldMatrix,
            backend = "bigsnpr",
            methodArgs = LassosumOptions()
        ),
        "backend 'bigsnpr' takes BigsnprOptions"
    )
    expect_error(BigsnprOptions(delta = 1), "lassosumRssWeights\\(\\)'s `s`")
    expect_error(BigsnprOptions(ncores = 2), "BiocParallel")
})
