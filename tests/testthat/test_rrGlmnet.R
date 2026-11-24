context("regularized_regression - glmnet")

# ---- glmnetWeights ----
test_that("glmnetWeights computes LASSO weights", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    y <- X[, 1] * 0.5 + rnorm(n)
    result <- glmnetWeights(X, y, alpha = 1)
    expect_equal(nrow(result), p)
    expect_equal(ncol(result), 1)
})

test_that("enetWeights computes elastic net weights", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    y <- X[, 1] * 0.5 + rnorm(n)
    result <- enetWeights(X, y)
    expect_equal(nrow(result), p)
})

test_that("glmnetWeights handles zero-variance columns", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    X[, 3] <- 5
    y <- X[, 1] * 0.5 + rnorm(n)
    expect_warning(
        result <- glmnetWeights(X, y, alpha = 1),
        "glmnetWeights: dropping 1 zero-variance column"
    )
    expect_equal(nrow(result), p)
    expect_equal(result[3, 1], 0)
})

test_that("glmnetWeights errors when all columns are constant", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 50
    p <- 5
    X <- matrix(rep(seq_len(p), each = n), nrow = n, ncol = p)
    y <- rnorm(n)

    expect_warning(
        expect_error(
            glmnetWeights(X, y, alpha = 1),
            "matrix with 2 or more columns"
        ),
        "dropping 5 zero-variance"
    )
})

test_that("glmnetWeights handles NA-producing columns gracefully", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    y <- X[, 1] * 0.5 + rnorm(n)

    X[1, 4] <- NA

    expect_warning(
        result <- glmnetWeights(X, y, alpha = 1),
        "glmnetWeights: dropping 1 zero-variance column"
    )
    expect_equal(nrow(result), p)
    expect_equal(ncol(result), 1)
    expect_equal(result[4, 1], 0)
})

test_that("enetWeights takes its mixing, and ridgeWeights fixes it at 0", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 60
    p <- 8
    X <- matrix(rnorm(n * p), nrow = n)
    y <- X[, 1] * 0.5 + rnorm(n)
    # cv.glmnet draws its folds, so each pair is fitted from the same seed.
    set.seed(11)
    mixed <- enetWeights(X, y, alpha = 0.3)
    set.seed(11)
    expect_equal(mixed, glmnetWeights(X, y, 0.3))
    set.seed(11)
    ridge <- ridgeWeights(X, y)
    set.seed(11)
    expect_equal(ridge, glmnetWeights(X, y, 0))
    expect_equal(nrow(ridge), p)
    # A ridge keeps every variant.
    expect_true(all(ridge != 0))
    expect_error(enetWeights(X, y, alpha = 2), "alpha.*<= 1")
})

test_that("glmnetWeights produces non-zero weights for correlated data", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 100
    p <- 10
    X <- matrix(rnorm(n * p), n, p)
    beta_true <- c(3, -2, rep(0, p - 2))
    y <- X %*% beta_true + rnorm(n)

    w <- glmnetWeights(X, y, alpha = 0.5)
    expect_length(w, p)
    expect_true(any(w != 0))
})

test_that("glmnetWeights keeps the selected fit or the whole path", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    X <- matrix(rnorm(500), nrow = 50)
    colnames(X) <- sprintf("chr1:%d:A:G", seq_len(10))
    X[, 4] <- 1
    y <- X[, 1] * 0.5 + rnorm(50)
    fitted <- colnames(X)[-4]
    set.seed(1)
    none <- suppressWarnings(glmnetWeights(X, y, alpha = 1))
    expect_null(attr(none, "fit"))
    set.seed(1)
    slim <- suppressWarnings(lassoWeights(X, y, fitRetention = "slim"))
    expect_equal(as.numeric(slim), as.numeric(none))
    fit <- attr(slim, "fit")
    expect_named(
        fit,
        c("coef", "intercept", "lambda", "cvm", "cvsd", "variantIds")
    )
    expect_equal(fit$coef, as.numeric(slim)[-4])
    expect_equal(fit$variantIds, fitted)
    full <- attr(
        suppressWarnings(ridgeWeights(X, y, fitRetention = "full")),
        "fit"
    )
    expect_s3_class(full, "cv.glmnet")
    # The path names its variants, so no second copy is added.
    expect_null(full$variantIds)
    expect_equal(pecotmr:::.twasFitVariantIds(full), fitted)
    # Its stored calls name the data rather than holding it...
    expect_identical(full$call$x, quote(x))
    expect_identical(full$glmnet.fit$call$x, quote(x))
    expect_identical(full$call$y, quote(y))
    # ...unless asked to keep it, and the coefficients do not depend on that.
    withX <- attr(
        suppressWarnings(
            ridgeWeights(X, y, fitRetention = "full", retainX = TRUE)
        ),
        "fit"
    )
    expect_true(is.matrix(withX$call$x))
    expect_gt(
        length(serialize(withX$call, NULL)) -
            length(serialize(full$call, NULL)),
        length(serialize(X[, fitted], NULL))
    )
    expect_equal(dim(withX$glmnet.fit$beta), dim(full$glmnet.fit$beta))
})

# ---- glmnetWeights across conditions (mgaussian) ----
test_that("glmnetWeights fits a matrix y with glmnet's mgaussian family", {
    skip_if_not_installed("glmnet")
    d <- .rrwMulti(n = 80, p = 8, K = 3)
    set.seed(1)
    w <- glmnetWeights(d$X, d$Y, alpha = 1, fitRetention = "slim")
    expect_equal(dimnames(w), list(colnames(d$X), colnames(d$Y)))
    # A variant is in every condition's model or in none.
    expect_true(all(rowSums(w != 0) %in% c(0, 3)))
    fit <- attr(w, "fit")
    expect_named(
        fit,
        c("coef", "intercept", "lambda", "cvm", "cvsd", "variantIds")
    )
    # glmnet's fit at the lambda it selected, on each phenotype scaled to
    # unit variance, put back on each one's own scale.
    scaled <- scale(d$Y)
    spread <- attr(scaled, "scaled:scale")
    direct <- glmnet::glmnet(
        d$X,
        scaled,
        family = "mgaussian",
        standardize = FALSE,
        lambda = fit$lambda
    )
    B <- do.call(cbind, lapply(coef(direct), as.numeric))
    expect_equal(
        unname(w),
        sweep(B[-1L, ], 2L, spread, "*"),
        tolerance = 1e-4,
        check.attributes = FALSE
    )
    expect_equal(
        fit$intercept,
        colMeans(d$Y) + spread * B[1L, ],
        tolerance = 1e-4,
        check.attributes = FALSE
    )
})

test_that("glmnet fits conditions on the samples observing all of them", {
    skip_if_not_installed("glmnet")
    d <- .rrwMulti(n = 40, p = 5, K = 2)
    Y <- replace(d$Y, cbind(1:10, 2L), NA)
    withr::local_seed(1L)
    partial <- lassoWeights(d$X, Y)
    withr::local_seed(1L)
    complete <- lassoWeights(d$X[11:40, ], d$Y[11:40, ])
    expect_equal(partial, complete)
    # Conditions on samples of their own share none.
    disjoint <- replace(d$Y, cbind(c(1:20, 21:40), rep(1:2, each = 20L)), NA)
    expect_error(lassoWeights(d$X, disjoint), "backend = \"grpreg\"")
})

test_that("the lasso and elastic net take their backend's options", {
    skip_if_not_installed("glmnet")
    d <- .rrwXy(n = 50, p = 6)
    expect_error(
        lassoWeights(d$X, d$y, methodArgs = NcvregOptions()),
        "backend 'glmnet' takes GlmnetOptions"
    )
    expect_error(
        enetWeights(d$X, d$y, backend = "ncvreg"),
        "must be one of glmnet, grpreg"
    )
})
