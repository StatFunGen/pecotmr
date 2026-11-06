context("regularized_regression - mrashRss")

# ============================================================================
# mrashRssWeights - dispatch + smoke test
# ============================================================================

test_that("mrashRssWeights forwards arguments to susieR::mr.ash.rss", {
    set.seed(42)
    p <- 10
    bhat <- rnorm(p, sd = 0.1)
    shat <- rep(0.05, p)
    R <- diag(p)
    # Heterogeneous n with median != mean: median = 55, mean = 145.
    stat <- list(
        b = bhat,
        seb = shat,
        n = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 1000),
        variantNames = .rrwIds(p)
    )
    z_vec <- bhat / shat
    captured <- new.env(parent = emptyenv())
    # mr.ash.rss is `@importFrom`'d into pecotmr's namespace, so mock the
    # pecotmr binding (no .package argument), not the susieR binding.
    local_mocked_bindings(
        mr.ash.rss = function(
            bhat,
            shat,
            z,
            R,
            var_y,
            n,
            sigma2_e,
            s0,
            w0,
            ...
        ) {
            captured$bhat <- bhat
            captured$shat <- shat
            captured$z <- z
            captured$R <- R
            captured$var_y <- var_y
            captured$n <- n
            captured$sigma2_e <- sigma2_e
            captured$s0 <- s0
            captured$w0 <- w0
            captured$dots <- list(...)
            list(mu1 = seq_len(length(bhat)) * 0.01)
        }
    )
    result <- mrashRssWeights(
        stat = stat,
        ldReference = .rrwLdMatrix(R),
        varY = 1.5,
        sigma2E = 0.4,
        s0 = c(0, 0.1, 0.2),
        w0 = c(0.5, 0.3, 0.2),
        updateOrder = "consecutive",
        methodArgs = MrashRssOptions(tol = 1e-6)
    )
    expect_equal(captured$bhat, bhat)
    expect_equal(captured$shat, shat)
    expect_equal(captured$z, z_vec)
    expect_equal(unname(captured$R), R)
    expect_equal(captured$var_y, 1.5)
    expect_equal(captured$n, 55) # median of stat$n, NOT mean (145)
    expect_equal(captured$sigma2_e, 0.4)
    expect_equal(captured$s0, c(0, 0.1, 0.2))
    expect_equal(captured$w0, c(0.5, 0.3, 0.2))
    expect_equal(captured$dots$tol, 1e-6)
    expect_equal(result, seq_len(p) * 0.01)
})

test_that("mrashRssWeights returns a numeric vector of expected length", {
    skip_if_not_installed("susieR")
    set.seed(2024)
    n <- 500
    p <- 10
    X <- matrix(rbinom(n * p, 2, 0.3), nrow = n)
    beta_true <- rep(0, p)
    beta_true[c(3, 7)] <- c(0.4, -0.3)
    y <- as.numeric(X %*% beta_true) + rnorm(n)
    bhat <- as.vector(cor(y, X))
    shat <- rep(1 / sqrt(n - 1), p)
    R <- cor(X)
    stat <- list(b = bhat, seb = shat, n = rep(n, p), variantNames = .rrwIds(p))
    w <- mrashRssWeights(
        stat = stat,
        ldReference = .rrwLdMatrix(R),
        varY = var(y),
        sigma2E = NULL,
        s0 = c(0, 0.01, 0.1, 0.5),
        w0 = rep(1 / 4, 4)
    )
    expect_equal(length(w), p)
    expect_true(all(is.finite(w)))
})

# ============================================================================
# mrashRssWeights - update order
# ============================================================================

# Each variant's marginal log Bayes factor under mr.ash.rss's mixture prior,
# from its sufficient statistics on the unstandardized scale, written
# independently as the ratio of the two marginal densities of bhat.
.mrr_logBf <- function(z, shat, n, varY, sigma2, s0, w0) {
    adj <- (n - 1) / (z^2 + n - 2)
    xtx <- varY * adj / shat^2
    bhat <- z * adj * varY / shat / xtx
    s2 <- sigma2 / xtx
    perComponent <- vapply(
        s0,
        function(s) {
            dnorm(bhat, 0, sqrt(sigma2 * s + s2), log = TRUE) -
                dnorm(bhat, 0, sqrt(s2), log = TRUE)
        },
        numeric(length(z))
    )
    apply(perComponent, 1, function(r) log(sum(w0 * exp(r))))
}

# Correlated variants of unequal precision, so a permutation shows in R and
# the Bayes factors' order differs from that of |z|.
.mrr_orderInput <- function(p = 8L) {
    set.seed(5)
    R <- 0.4^abs(outer(seq_len(p), seq_len(p), "-"))
    shat <- seq(0.03, 0.1, length.out = p)
    bhat <- c(0.2, -0.05, 0.3, 0.01, -0.25, 0.08, 0.15, -0.02)[seq_len(p)]
    list(
        R = R,
        stat = list(
            b = bhat,
            seb = shat,
            n = rep(400, p),
            variantNames = .rrwIds(p)
        ),
        z = bhat / shat
    )
}

test_that("decreasingLogBf forwards the variants permuted, then unpermutes", {
    d <- .mrr_orderInput()
    p <- length(d$z)
    s0 <- c(0, 0.01, 0.1)
    w0 <- c(0.5, 0.3, 0.2)
    captured <- new.env(parent = emptyenv())
    local_mocked_bindings(
        mr.ash.rss = function(bhat, shat, z, R, mu1_init, ...) {
            captured$bhat <- bhat
            captured$shat <- shat
            captured$z <- z
            captured$R <- R
            captured$mu1_init <- mu1_init
            list(
                mu1 = matrix(z * 0.01, ncol = 1L),
                beta = z * 0.01,
                sigma2_1 = matrix(z^2, ncol = 1L),
                w1 = cbind(z, 2 * z, 3 * z),
                w0 = w0
            )
        }
    )
    w <- mrashRssWeights(
        d$stat,
        .rrwLdMatrix(d$R),
        varY = 1.2,
        sigma2E = 0.9,
        s0 = s0,
        w0 = w0,
        updateOrder = "decreasingLogBf",
        fitRetention = "full",
        methodArgs = MrashRssOptions(mu1_init = seq_len(p) / 10)
    )
    expected <- order(-.mrr_logBf(d$z, d$stat$seb, 400, 1.2, 0.9, s0, w0))
    expect_false(identical(expected, order(-abs(d$z))))
    perm <- match(captured$z, d$z)
    expect_identical(perm, expected)
    expect_equal(captured$bhat, d$stat$b[perm])
    expect_equal(captured$shat, d$stat$seb[perm])
    expect_equal(unname(captured$R), d$R[perm, perm])
    expect_equal(captured$mu1_init, (seq_len(p) / 10)[perm])
    # The fit comes back in the statistics' order, its shapes kept.
    expect_equal(as.numeric(w), d$z * 0.01)
    fit <- attr(w, "fit")
    expect_equal(dim(fit$mu1), c(p, 1L))
    expect_equal(fit$beta, d$z * 0.01)
    expect_equal(as.numeric(fit$sigma2_1), d$z^2)
    expect_equal(unname(fit$w1), cbind(d$z, 2 * d$z, 3 * d$z))
    expect_equal(fit$w0, w0)
})

test_that("z alone forwards an empty bhat and seb, whatever the order", {
    d <- .mrr_orderInput()
    captured <- new.env(parent = emptyenv())
    local_mocked_bindings(
        mr.ash.rss = function(bhat, shat, z, ...) {
            captured$bhat <- bhat
            captured$shat <- shat
            list(mu1 = z * 0)
        }
    )
    stat <- list(z = d$z, n = rep(400, length(d$z)), variantNames = .rrwIds(8L))
    mrashRssWeights(stat, .rrwLdMatrix(d$R), updateOrder = "decreasingLogBf")
    expect_identical(captured$bhat, numeric(0))
    expect_identical(captured$shat, numeric(0))
})

test_that("the order is mr.ash.rss itself run on the permuted variants", {
    skip_if_not_installed("susieR")
    d <- .mrr_orderInput()
    s0 <- c(0, 0.01, 0.1)
    w <- mrashRssWeights(
        d$stat,
        .rrwLdMatrix(d$R),
        varY = 1.2,
        s0 = s0,
        updateOrder = "decreasingLogBf",
        fitRetention = "full"
    )
    perm <- order(
        -.mrr_logBf(d$z, d$stat$seb, 400, 1.2, 1.2, s0, rep(1 / 3, 3))
    )
    manual <- susieR::mr.ash.rss(
        bhat = d$stat$b[perm],
        shat = d$stat$seb[perm],
        z = d$z[perm],
        R = d$R[perm, perm],
        var_y = 1.2,
        n = 400,
        s0 = s0,
        w0 = rep(1 / 3, 3)
    )
    inverse <- order(perm)
    expect_equal(as.numeric(w), as.numeric(manual$mu1)[inverse])
    fit <- attr(w, "fit")
    expect_equal(fit$w1, manual$w1[inverse, , drop = FALSE])
    expect_equal(fit$sigma2_1, manual$sigma2_1[inverse, , drop = FALSE])
    # Column order is the engine's own call unchanged.
    plain <- mrashRssWeights(
        d$stat,
        .rrwLdMatrix(d$R),
        varY = 1.2,
        s0 = s0,
        updateOrder = "consecutive"
    )
    direct <- susieR::mr.ash.rss(
        bhat = d$stat$b,
        shat = d$stat$seb,
        z = d$z,
        R = d$R,
        var_y = 1.2,
        n = 400,
        s0 = s0,
        w0 = rep(1 / 3, 3)
    )
    expect_equal(as.numeric(plain), as.numeric(direct$mu1))
})

test_that("on a common x'x the Bayes-factor order is the |z| order", {
    d <- .mrr_orderInput()
    inp <- list(
        b = d$stat$b,
        shat = d$stat$seb,
        z = d$z,
        n = 400,
        varY = 1.2,
        sigma2E = NULL
    )
    prior <- list(s0 = c(0, 0.01, 0.1), w0 = rep(1 / 3, 3))
    orderOf <- pecotmr:::.mrashRssOrder
    # The standardized scale: no varY with the seb.
    zOnly <- list_assign(inp, shat = numeric(0), varY = NULL)
    expect_identical(
        orderOf("decreasingLogBf", zOnly, prior, MrashRssOptions()),
        order(-abs(d$z))
    )
    # mr.ash.rss's own standardization makes every x'x n - 1 too.
    expect_identical(
        orderOf(
            "decreasingLogBf",
            inp,
            prior,
            MrashRssOptions(standardize = TRUE)
        ),
        order(-abs(d$z))
    )
    expect_null(orderOf("consecutive", inp, prior, MrashRssOptions()))
})
