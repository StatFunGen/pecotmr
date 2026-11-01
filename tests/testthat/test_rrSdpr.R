context("regularized_regression - sdpr")

# ---- sdpr ----
test_that("sdpr errors on mismatched bhat and R dimensions", {
    expect_error(
        sdpr(bhat = rnorm(10), R = diag(5), n = 100),
        "bhat.*Must have length"
    )
})

test_that("sdpr errors on non-matrix R input", {
    expect_error(
        sdpr(bhat = rnorm(5), R = "not_a_matrix", n = 100),
        "as a matrix"
    )
})

test_that("sdpr errors on non-positive sample size", {
    expect_error(
        sdpr(bhat = rnorm(5), R = diag(5), n = -1),
        "valid sample size"
    )
})

test_that("sdpr errors when M is less than 4", {
    expect_error(
        sdpr(bhat = rnorm(5), R = diag(5), n = 100, M = 3),
        "'M' must be at least 4"
    )
})

test_that("sdpr errors on invalid per_variant_sample_size", {
    expect_error(
        sdpr(
            bhat = rnorm(5),
            R = diag(5),
            n = 100,
            perVariantSampleSize = c(100, -1, 100, 100, 100)
        ),
        "positive values"
    )
})

test_that("sdpr errors on invalid array values", {
    expect_error(
        sdpr(bhat = rnorm(5), R = diag(5), n = 100, array = c(0, 1, 3, 1, 0)),
        "0, 1, or 2"
    )
})

test_that("sdpr runs successfully", {
    set.seed(42)
    p <- 10
    bhat <- rnorm(p, sd = 0.1)
    R <- diag(p)
    result <- sdpr(
        bhat = bhat,
        R = R,
        n = 100,
        iter = 50,
        burn = 10,
        thin = 2,
        verbose = FALSE
    )
    expect_type(result, "list")
    expect_true("betaEst" %in% names(result))
    expect_equal(length(result$betaEst), p)
    expect_true(all(is.finite(result$betaEst)))
})

test_that("sdpr runs with a per-variant sample size or an array vector", {
    set.seed(42)
    p <- 10
    bhat <- rnorm(p, sd = 0.1)
    R <- diag(p)
    extras <- list(
        perVariantSampleSize = list(perVariantSampleSize = rep(100, p)),
        array = list(array = rep(1, p))
    )
    for (nm in names(extras)) {
        result <- exec(
            sdpr,
            bhat = bhat,
            R = R,
            n = 100,
            !!!extras[[nm]],
            iter = 50,
            burn = 10,
            thin = 2,
            verbose = FALSE
        )
        expect_equal(length(result$betaEst), p, info = nm)
        expect_true(all(is.finite(result$betaEst)), info = nm)
    }
})

# ---- sdpr signal recovery ----
test_that("sdpr recovers signal direction on simulated genotype data", {
    set.seed(2024)
    n <- 500
    p <- 20
    # Realistic genotype matrix with non-trivial LD (binomial, MAF=0.3)
    X <- matrix(rbinom(n * p, 2, 0.3), nrow = n)
    beta_true <- rep(0, p)
    beta_true[c(3, 10, 15)] <- c(0.4, -0.3, 0.2)
    y <- X %*% beta_true + rnorm(n)
    bhat <- as.vector(cor(y, X))
    R <- cor(X)
    result <- sdpr(
        bhat = bhat,
        R = R,
        n = n,
        iter = 200,
        burn = 50,
        thin = 5,
        verbose = FALSE,
        seed = 42L
    )
    expect_true("betaEst" %in% names(result))
    expect_equal(length(result$betaEst), p)
    expect_true(all(is.finite(result$betaEst)))
    # Correlation with truth should be positive (signal recovery)
    expect_gt(cor(result$betaEst, beta_true), 0.3)
})

# SDPR on 20 variants (n = 500), two of them causal.
.sdprTwoCausalFit <- function(seed, raoBlackwell = FALSE, nChains = 1) {
    set.seed(2025)
    n <- 500
    p <- 20
    X <- matrix(rbinom(n * p, 2, 0.3), nrow = n)
    y <- X[, c(3, 10)] %*% c(0.4, -0.3) + rnorm(n)
    sdpr(
        bhat = as.vector(cor(y, X)),
        R = cor(X),
        n = n,
        iter = 200,
        burn = 50,
        thin = 5,
        verbose = FALSE,
        seed = seed,
        raoBlackwell = raoBlackwell,
        nChains = nChains
    )
}

test_that("Rao-Blackwellized SDPR runs the same chain", {
    draws <- .sdprTwoCausalFit(42L)
    means <- .sdprTwoCausalFit(42L, raoBlackwell = TRUE)
    expect_identical(means$h2, draws$h2)
    expect_false(identical(means$betaEst, draws$betaEst))
    expect_gt(cor(means$betaEst, draws$betaEst), 0.95)
})

test_that("sdpr errors when every one of several chains diverges", {
    # Statistics of 1e-6 on one sample leave the draws near the prior's, far
    # past twice the statistics' sum of squares.
    expect_error(
        sdpr(
            rep(1e-6, 5),
            diag(5),
            1,
            iter = 30,
            burn = 10,
            thin = 5,
            verbose = FALSE,
            seed = 3L,
            nChains = 2
        ),
        "sdpr: every chain diverged"
    )
})

test_that("sdpr averages several chains over the kept ones", {
    # SDPR's first sweeps draw far past the divergence test's bound before
    # settling, which the test, past the burn-in, leaves alone.
    many <- .sdprTwoCausalFit(4L, nChains = 3)
    expect_equal(many$chains$seed, c(4, 5, 6))
    expect_false(any(many$chains$diverged))
    kept <- many$chains$kept
    singles <- map(4:6, .sdprTwoCausalFit)[kept]
    expect_equal(
        many$betaEst,
        reduce(map(singles, "betaEst"), `+`) / sum(kept)
    )
    expect_equal(many$h2, reduce(map(singles, "h2"), `+`) / sum(kept))
})

# A window of 30 variants in LD (n = 600) whose 15th variant alone has a
# standardized effect explaining `h2` of the trait, fitted by SDPR.
.sdprLoneCausalFit <- function(seed, h2) {
    set.seed(seed)
    n <- 600
    p <- 30
    X <- matrix(rbinom(n * p, 2, 0.3), n, p)
    for (j in 2:p) {
        X[, j] <- ifelse(runif(n) < 0.6, X[, j - 1], X[, j])
    }
    X <- scale(X)
    y <- X[, 15] * sqrt(h2) + rnorm(n, sd = sqrt(1 - h2))
    sdpr(
        bhat = as.vector(cor(X, y)),
        R = cor(X),
        n = n,
        M = 20,
        iter = 600,
        burn = 100,
        thin = 5,
        verbose = FALSE,
        seed = 1L
    )
}

test_that("sdpr recovers a lone causal variant's effect and h2", {
    # A lone causal variant is drawn as several are, so eta's update reads
    # N beta'A'bhat. The original's closed-form draw took eta times that,
    # which centred eta near eta^2: on such a window the effect and h2 ran
    # away (h2 in the tens of thousands here). The truth is an effect of
    # sqrt(0.1) and h2 = 0.1.
    result <- .sdprLoneCausalFit(1001, 0.1)
    expect_gt(result$h2, 0.03)
    expect_lt(result$h2, 0.3)
    expect_gt(result$betaEst[15], 0.15)
    expect_lt(result$betaEst[15], 0.6)
})

test_that("sdpr finds a weak lone causal variant rather than staying null", {
    # Every chain starts with no causal variant, and while a sweep includes
    # none, eta is drawn from its prior alone. The original's flat prior, in
    # effect N(0, 1e6), put |eta| in the hundreds, which made every
    # inclusion unlikely, so the chain stayed null. On this window (causal
    # z = 4.8, h2 = 0.03) all of 500 seeds found the variant under the
    # N(1, 1) prior and 4 of 100 under the original's.
    result <- .sdprLoneCausalFit(2003, 0.03)
    expect_gt(result$h2, 0.01)
    expect_lt(result$h2, 0.2)
})

# ---- sdpr verbose output ----
test_that("sdpr with verbose = TRUE produces output", {
    set.seed(42)
    p <- 10
    bhat <- rnorm(p, sd = 0.1)
    R <- diag(p)
    # iter >= 100 triggers the verbose print inside the MCMC loop
    out <- capture.output(
        result <- sdpr(
            bhat = bhat,
            R = R,
            n = 100,
            iter = 110,
            burn = 10,
            thin = 2,
            verbose = TRUE,
            seed = 42L
        )
    )
    expect_match(out, "^100 iter\\. h2:", all = FALSE)
    expect_type(result, "list")
    expect_equal(length(result$betaEst), p)
    expect_true(all(is.finite(result$betaEst)))
})

# ---- sdpr optLlk = 2 (multi-array) ----
test_that("sdpr runs with optLlk = 2 on mixed arrays and genotype LD", {
    set.seed(2024)
    n <- 500
    p <- 20
    X <- matrix(rbinom(n * p, 2, 0.3), nrow = n)
    beta_true <- rep(0, p)
    beta_true[c(3, 15)] <- c(0.4, 0.2)
    y <- X %*% beta_true + rnorm(n)
    bhat <- as.vector(cor(y, X))
    R <- cor(X)
    # Mixed arrays with varying sample sizes
    arr <- rep(c(1, 2), length.out = p)
    per_n <- rep(c(400, 500), length.out = p)
    result <- sdpr(
        bhat = bhat,
        R = R,
        n = n,
        perVariantSampleSize = per_n,
        array = arr,
        optLlk = 2,
        iter = 100,
        burn = 30,
        thin = 5,
        verbose = FALSE,
        seed = 42L
    )
    expect_type(result, "list")
    expect_equal(length(result$betaEst), p)
    expect_true(all(is.finite(result$betaEst)))
})

# ---- sdprWeights (wrapper) ----
test_that("sdprWeights calls sdpr and returns beta_est", {
    set.seed(42)
    p <- 10
    bhat <- rnorm(p, sd = 0.1)
    R <- diag(p)
    stat <- list(b = bhat, n = rep(100, p), variantNames = .rrwIds(p))
    result <- sdprWeights(
        stat = stat,
        ldReference = .rrwLdMatrix(R),
        methodArgs = SdprOptions(
            iter = 50,
            burn = 10,
            thin = 2,
            verbose = FALSE
        )
    )
    expect_equal(length(result), p)
    expect_true(is.numeric(result))
})

test_that("sdpr runs on a factor that carries no LD at all", {
    # Every variant monomorphic: the factor has rank 0, so B and A vanish and
    # nothing is heritable.
    fit <- pecotmr:::sdprFactorCpp(
        bhatR = rep(0.01, 5),
        factorR = matrix(0, 10, 5),
        n = 100L,
        perVariantSampleSize = NULL,
        array = NULL,
        a = 0.1,
        c = 1,
        M = 20L,
        a0k = 0.5,
        b0k = 0.5,
        iter = 50L,
        burn = 10L,
        thin = 1L,
        optLlk = 1L,
        verbose = FALSE,
        seed = 1L,
        raoBlackwell = FALSE,
        track = FALSE
    )
    expect_equal(fit$h2, 0)
    expect_true(all(is.finite(fit$betaEst)))
})
