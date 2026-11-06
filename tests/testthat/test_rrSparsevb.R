context("regularized_regression - sparsevb")

.sv_data <- function(n = 60L, p = 8L) {
    set.seed(7)
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- sprintf("chr1:%d:A:G", 100L * seq_len(p))
    y <- as.numeric(X[, 1] * 0.6 - X[, 3] * 0.4 + rnorm(n))
    list(X = X, y = y)
}

# ---- sparsevbWeights ----
test_that("sparsevbWeights is svb.fit's posterior mean, gamma * mu", {
    skip_if_not_installed("sparsevb")
    d <- .sv_data()
    set.seed(1)
    slim <- sparsevbWeights(d$X, d$y, fitRetention = "slim")
    # With no identical columns, the noise SD estimated on the distinct
    # columns is the engine's own, drawn at the same point of the stream.
    set.seed(1)
    fit <- sparsevb::svb.fit(d$X, d$y)
    expect_equal(as.numeric(slim), fit$mu * fit$gamma)
    kept <- attr(slim, "fit")
    expect_equal(kept[c("mu", "sigma", "gamma", "intercept")], fit)
    expect_equal(kept$variantIds, colnames(d$X))
    # The posterior holds no copy of the data, so "full" keeps the same.
    set.seed(1)
    full <- sparsevbWeights(d$X, d$y, fitRetention = "full")
    expect_identical(attr(full, "fit"), kept)
    expect_null(attr(sparsevbWeights(d$X, d$y), "fit"))
})

test_that("sparsevbWeights estimates the noise SD on the distinct columns", {
    skip_if_not_installed("sparsevb")
    d <- .sv_data()
    # Variants in perfect LD give identical columns, which sparsevb's own
    # noise estimate refuses.
    X <- cbind(d$X, `chr1:900:A:G` = d$X[, 2])
    expect_error(sparsevb::svb.fit(X, d$y), "duplicate columns")
    set.seed(1)
    w <- sparsevbWeights(X, d$y)
    set.seed(1)
    noise <- selectiveInference::estimateSigma(d$X, d$y)$sigmahat
    fit <- sparsevb::svb.fit(X, d$y, noise_sd = noise)
    expect_equal(as.numeric(w), fit$mu * fit$gamma)
    # A noise SD the caller gives is used as it is.
    set.seed(1)
    given <- sparsevbWeights(
        X,
        d$y,
        methodArgs = SparsevbOptions(noise_sd = 0.8)
    )
    set.seed(1)
    fit <- sparsevb::svb.fit(X, d$y, noise_sd = 0.8)
    expect_equal(as.numeric(given), fit$mu * fit$gamma)
})

test_that("sparsevbWeights warns, zero-pads X's zero-variance columns", {
    skip_if_not_installed("sparsevb")
    d <- .sv_data()
    X <- d$X
    X[, 6] <- 4
    expect_warning(
        w <- sparsevbWeights(X, d$y, fitRetention = "slim"),
        "sparsevbWeights: dropping 1 zero-variance column"
    )
    expect_length(w, ncol(X))
    expect_equal(w[[6]], 0)
    expect_equal(attr(w, "fit")$variantIds, colnames(X)[-6])
})

# ---- SparsevbOptions ----
test_that("SparsevbOptions checks its names against svb.fit", {
    skip_if_not_installed("sparsevb")
    expect_error(SparsevbOptions(X = 1), "supplied by pecotmr")
    expect_error(SparsevbOptions(Y = 1), "supplied by pecotmr")
    expect_error(SparsevbOptions(noiseSd = 1), "unknown argument")
    expect_equal(
        as.list(SparsevbOptions(slab = "gaussian")),
        list(slab = "gaussian")
    )
    d <- .sv_data()
    expect_error(
        sparsevbWeights(d$X, d$y, methodArgs = list(slab = "gaussian")),
        "SparsevbOptions"
    )
})
