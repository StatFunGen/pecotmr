context("regularized_regression - grpreg")

# ---- grpregWeights on one phenotype ----
test_that("on one phenotype grpreg fits ncvreg's MCP and SCAD", {
    skip_if_not_installed("grpreg")
    skip_if_not_installed("ncvreg")
    d <- .rrwXy(n = 80, p = 10)
    grid <- c(0.4, 0.2, 0.1)
    for (penalty in c("MCP", "SCAD")) {
        gamma <- if (penalty == "SCAD") 3.7 else 3
        viaGrpreg <- grpregWeights(
            d$X,
            d$y,
            paste0("gr", penalty),
            fitRetention = "full",
            methodArgs = GrpregOptions(
                lambda = grid,
                gamma = gamma,
                eps = 1e-10
            )
        )
        viaNcvreg <- ncvregWeights(
            d$X,
            d$y,
            penalty,
            fitRetention = "full",
            methodArgs = NcvregOptions(
                lambda = grid,
                gamma = gamma,
                eps = 1e-10
            )
        )
        expect_equal(
            unname(attr(viaGrpreg, "fit")$fit$beta),
            unname(attr(viaNcvreg, "fit")$fit$beta),
            tolerance = 1e-6
        )
    }
})

test_that("mcp, scad, lasso and enet run on grpreg as their backend", {
    skip_if_not_installed("grpreg")
    d <- .rrwXy(n = 80, p = 8)
    set.seed(1)
    w <- mcpWeights(d$X, d$y, backend = "grpreg", fitRetention = "slim")
    expect_equal(dim(w), c(8L, 1L))
    expect_named(
        attr(w, "fit"),
        c("coef", "intercept", "lambda", "cve", "cvse", "variantIds")
    )
    set.seed(1)
    expect_equal(dim(scadWeights(d$X, d$y, backend = "grpreg")), c(8L, 1L))
    set.seed(1)
    expect_equal(dim(lassoWeights(d$X, d$y, backend = "grpreg")), c(8L, 1L))
    set.seed(1)
    expect_equal(
        dim(enetWeights(d$X, d$y, alpha = 0.3, backend = "grpreg")),
        c(8L, 1L)
    )
    # The lasso and elastic net take their mixing from the wrapper, and
    # grpreg has no ridge.
    expect_error(
        lassoWeights(
            d$X,
            d$y,
            backend = "grpreg",
            methodArgs = GrpregOptions(alpha = 0.5)
        ),
        "enetWeights takes its own `alpha`"
    )
    expect_error(
        enetWeights(d$X, d$y, alpha = 0, backend = "grpreg"),
        "no ridge fit"
    )
})

# ---- grpregWeights across conditions ----
test_that("grpreg fits conditions on the same samples as its multitask model", {
    skip_if_not_installed("grpreg")
    d <- .rrwMulti(n = 80, p = 8, K = 3)
    grid <- c(0.3, 0.15)
    set.seed(1)
    w <- grpregWeights(
        d$X,
        d$Y,
        "grLasso",
        fitRetention = "full",
        methodArgs = GrpregOptions(lambda = grid, eps = 1e-12, max.iter = 1e6)
    )
    expect_equal(dimnames(w), list(colnames(d$X), colnames(d$Y)))
    expect_true(all(rowSums(w != 0) %in% c(0, 3)))
    # grpreg's own multitask fit, on X centred and each phenotype scaled to
    # unit variance, at the lambda the cross-validation picked.
    scaled <- scale(d$Y)
    direct <- grpreg::grpreg(
        scale(d$X, scale = FALSE),
        scaled,
        penalty = "grLasso",
        lambda = grid,
        eps = 1e-12,
        max.iter = 1e6
    )
    at <- match(attr(w, "fit")$lambda.min, grid)
    expect_equal(
        unname(w),
        sweep(t(direct$beta[, -1L, at]), 2L, attr(scaled, "scaled:scale"), "*"),
        tolerance = 1e-5,
        check.attributes = FALSE
    )
})

test_that("grpreg fits conditions on samples of their own, each on its loss", {
    skip_if_not_installed("grpreg")
    set.seed(4)
    n <- 150
    p <- 8
    X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, .rrwIds(p)))
    first <- seq_len(n) <= 70
    Y <- cbind(
        a = ifelse(first, X[, 1] - X[, 2] + rnorm(n), NA),
        b = ifelse(first, NA, 0.7 * X[, 1] + rnorm(n, sd = 2))
    )
    set.seed(1)
    w <- grpregWeights(
        X,
        Y,
        "grLasso",
        fitRetention = "slim",
        methodArgs = GrpregOptions(
            lambda = c(0.3, 0.1, 0.05),
            eps = 1e-12,
            max.iter = 1e6
        )
    )
    expect_true(all(rowSums(w != 0) %in% c(0, 2)))
    # The optimality conditions of each condition's own mean squared error
    # on its phenotype scaled to unit variance, plus the group lasso on its
    # per-SD effects: each fitted on its own samples.
    lam <- 2 * attr(w, "fit")$lambda
    parts <- lapply(1:2, function(k) {
        rows <- !is.na(Y[, k])
        x <- scale(X[rows, ], scale = FALSE)
        sdx <- sqrt(colMeans(x^2))
        spread <- sd(Y[rows, k])
        b <- w[, k] / spread
        r <- (Y[rows, k] - mean(Y[rows, k])) / spread - x %*% b
        list(z = crossprod(x, r) / (sum(rows) * sdx), u = b * sdx)
    })
    z <- sapply(parts, `[[`, "z")
    u <- sapply(parts, `[[`, "u")
    on <- rowSums(u != 0) > 0
    expect_true(any(on))
    expect_equal(
        unname(z[on, , drop = FALSE]),
        unname(
            lam * u[on, , drop = FALSE] / sqrt(rowSums(u[on, , drop = FALSE]^2))
        ),
        tolerance = 1e-5
    )
    expect_true(all(sqrt(rowSums(z[!on, , drop = FALSE]^2)) <= lam + 1e-8))
})

test_that("grpregWeights keeps the selected fit or the path without data", {
    skip_if_not_installed("grpreg")
    d <- .rrwMulti(n = 60, p = 6, K = 2)
    set.seed(1)
    slim <- attr(
        mcpWeights(d$X, d$Y, backend = "grpreg", fitRetention = "slim"),
        "fit"
    )
    expect_named(
        slim,
        c("coef", "intercept", "lambda", "cve", "cvse", "variantIds")
    )
    expect_equal(dim(slim$coef), c(6L, 2L))
    set.seed(1)
    full <- attr(
        mcpWeights(d$X, d$Y, backend = "grpreg", fitRetention = "full"),
        "fit"
    )
    expect_s3_class(full, "cv.grpreg")
    expect_null(full$fit$y)
    expect_null(full$fit$linear.predictors)
    expect_null(full$fit$XG)
    set.seed(1)
    withX <- attr(
        mcpWeights(
            d$X,
            d$Y,
            backend = "grpreg",
            fitRetention = "full",
            retainX = TRUE
        ),
        "fit"
    )
    expect_false(is.null(withX$fit$XG))
})

test_that("ncvreg refuses several conditions' phenotypes", {
    skip_if_not_installed("ncvreg")
    d <- .rrwMulti(n = 40, p = 5, K = 2)
    expect_error(mcpWeights(d$X, d$Y), "backend = \"grpreg\"")
})

# ---- GrpregOptions ----
test_that("GrpregOptions refuses what the wrappers set and checks the rest", {
    skip_if_not_installed("grpreg")
    expect_error(GrpregOptions(group = 1:3), "a variant's effects across")
    expect_error(GrpregOptions(penalty = "grMCP"), "the calling wrapper")
    expect_error(GrpregOptions(nfolds = 3), "the wrapper's `nfolds`")
    expect_error(GrpregOptions(returnX = TRUE), "`retainX`")
    expect_error(GrpregOptions(gama = 2), "unknown argument")
    expect_equal(GrpregOptions(gamma = 2.5)$gamma, 2.5)
    expect_equal(pecotmr:::.optionsInputType(GrpregOptions()), "QtlDataset")
})
