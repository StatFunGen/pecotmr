context("regularized_regression - mrash")

# ---- mrashWeights ----
test_that("mrashWeights subsets a full-length beta.init by keep", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    X[, 5] <- 7 # zero-variance column to force subset by keep
    y <- X[, 1] * 0.5 + rnorm(n)
    user_beta_init <- seq(0.01, by = 0.01, length.out = p)
    expect_warning(
        result <- mrashWeights(
            X,
            y,
            methodArgs = MrashOptions(beta.init = user_beta_init)
        ),
        "mrashWeights: dropping 1 zero-variance column"
    )
    expect_equal(length(result), p)
    expect_equal(result[5], 0)
})

test_that("mrashWeights warns and zero-pads when X has zero-variance columns", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 50
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    X[, 3] <- 7
    y <- X[, 1] * 0.5 + rnorm(n)
    local_mocked_bindings(
        lassoWeights = function(X, y) rep(0.01, ncol(X))
    )
    expect_warning(
        result <- mrashWeights(X, y),
        "mrashWeights: dropping 1 zero-variance column"
    )
    expect_equal(length(result), p)
    expect_equal(result[3], 0)
})

gc()

test_that("mrashWeights' slim fit drops the copy of the data", {
    set.seed(42)
    X <- matrix(rnorm(500), nrow = 50)
    colnames(X) <- sprintf("chr1:%d:A:G", seq_len(10))
    y <- X[, 1] * 0.5 + rnorm(50)
    set.seed(1)
    slim <- attr(mrashWeights(X, y, fitRetention = "slim"), "fit")
    expect_s3_class(slim, "mr.ash")
    expect_false(any(is_in(c("data", "update.order"), names(slim))))
    # estimateSparsity() reads the mixture weights it keeps.
    expect_true(is.numeric(slim$pi))
    expect_equal(slim$variantIds, colnames(X))
    set.seed(1)
    full <- attr(mrashWeights(X, y, fitRetention = "full"), "fit")
    # All of the engine's object but its copy of X, so it names its variants.
    expect_true(is.list(full$data))
    expect_null(full$data$X)
    expect_false(is.null(full$update.order))
    expect_equal(full$variantIds, colnames(X))
    expect_equal(pecotmr:::.twasFitVariantIds(full), colnames(X))
    withX <- attr(
        mrashWeights(X, y, fitRetention = "full", retainX = TRUE),
        "fit"
    )
    expect_identical(colnames(withX$data$X), colnames(X))
    expect_null(withX$variantIds)
    expect_equal(pecotmr:::.twasFitVariantIds(withX), colnames(X))
})

# ---- mrashWeights' update order ----

# mr.ash's prior for X and y, built independently of the package: centred
# data, its default grid and equal weights, and the variance of y.
.mrw_logBf <- function(X, y) {
    Xc <- scale(X, center = TRUE, scale = FALSE)
    yc <- y - mean(y)
    w <- colSums(Xc^2)
    sa2 <- (2^((0:24) / 25) - 1)^2 / median(w) * nrow(X)
    sigma2 <- var(y)
    bhat <- as.numeric(crossprod(Xc, yc)) / w
    s2 <- sigma2 / w
    perComponent <- vapply(
        sa2,
        function(s) {
            dnorm(bhat, 0, sqrt(sigma2 * s + s2), log = TRUE) -
                dnorm(bhat, 0, sqrt(s2), log = TRUE)
        },
        numeric(ncol(X))
    )
    apply(perComponent, 1, function(r) log(mean(exp(r))))
}

.mrw_data <- function() {
    set.seed(9)
    X <- matrix(rbinom(120 * 12, 2, 0.3), 120, 12)
    y <- as.numeric(X %*% replace(numeric(12), c(3, 8), c(0.5, -0.4))) +
        rnorm(120)
    list(X = X, y = y)
}

test_that("mrashWeights passes mr.ash its update order", {
    d <- .mrw_data()
    captured <- new.env(parent = emptyenv())
    local_mocked_bindings(
        lassoWeights = function(X, y) rep(0, ncol(X)),
        mr.ash = function(X, y, ...) {
            captured$args <- list(...)
            structure(list(beta = rep(0, ncol(X))), class = "mr.ash")
        },
        predict = function(object, ...) rep(0, length(object$beta) + 1L)
    )
    mrashWeights(d$X, d$y)
    expect_false(is_in("update.order", names(captured$args)))
    mrashWeights(d$X, d$y, updateOrder = "decreasingLogBf")
    expect_identical(
        captured$args$update.order,
        order(-.mrw_logBf(d$X, d$y))
    )
})

test_that("decreasingLogBf is mr.ash run in that order", {
    d <- .mrw_data()
    local_mocked_bindings(lassoWeights = function(X, y) rep(0.01, ncol(X)))
    w <- mrashWeights(d$X, d$y, updateOrder = "decreasingLogBf")
    direct <- susieR::mr.ash(
        d$X,
        d$y,
        beta.init = rep(0.01, ncol(d$X)),
        update.order = order(-.mrw_logBf(d$X, d$y)),
        verbose = FALSE
    )
    expect_equal(w, predict(direct, type = "coefficients")[-1])
    expect_error(
        MrashOptions(update.order = 1:3),
        "update.order.*mrashWeights\\(\\)' own `updateOrder`"
    )
})

test_that(".mixtureLogBf is the marginal likelihood ratio it claims", {
    # Against numerical integration over the effect, for a grid without the
    # point mass integrate() cannot see.
    xtx <- c(50, 200)
    xty <- c(4, -1)
    sigma2 <- 0.8
    scales <- c(0.01, 0.1)
    weights <- c(0.3, 0.7)
    numeric <- vapply(
        seq_along(xtx),
        function(j) {
            bhat <- xty[[j]] / xtx[[j]]
            s2 <- sigma2 / xtx[[j]]
            marginal <- sum(vapply(
                seq_along(scales),
                function(k) {
                    integrate(
                        function(b) {
                            dnorm(bhat, b, sqrt(s2)) *
                                dnorm(b, 0, sqrt(sigma2 * scales[[k]]))
                        },
                        -Inf,
                        Inf
                    )$value *
                        weights[[k]]
                },
                numeric(1)
            ))
            log(marginal / dnorm(bhat, 0, sqrt(s2)))
        },
        numeric(1)
    )
    expect_equal(
        pecotmr:::.mixtureLogBf(xtx, xty, sigma2, scales, weights),
        numeric,
        tolerance = 1e-6
    )
})
