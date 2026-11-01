context("regularized_regression - prsCs")

# ---- prsCs ----
test_that("prsCs errors on non-matrix R input", {
    expect_error(
        prsCs(bhat = rnorm(5), R = "not_a_matrix", n = 100),
        "as a matrix"
    )
})

test_that("prsCs errors on non-positive sample size", {
    expect_error(
        prsCs(bhat = rnorm(5), R = diag(5), n = -1),
        "valid sample size"
    )
})

test_that("prsCs errors on mismatched maf length", {
    expect_error(
        prsCs(bhat = rnorm(5), R = diag(5), n = 100, maf = rep(0.3, 3)),
        "same as 'maf'"
    )
})

test_that("prsCs errors on mismatched bhat and R dimensions", {
    expect_error(
        prsCs(bhat = rnorm(10), R = diag(5), n = 100),
        "bhat.*Must have length"
    )
})

test_that("prsCs averages exactly the draws its schedule keeps", {
    # Iterations 5 and 10 are kept, two draws, though (12 - 3) / 5 rounds to
    # one. A fixed phi is the same in every draw, so its mean is phi.
    result <- prsCs(
        bhat = rep(0.05, 5),
        R = diag(5),
        n = 100,
        phi = 0.01,
        nIter = 12,
        nBurnin = 3,
        thin = 5,
        seed = 1L
    )
    expect_equal(result$phiEst, 0.01)
    expect_error(
        prsCs(bhat = rep(0.05, 5), R = diag(5), n = 100, thin = 0),
        "thinning interval positive"
    )
    expect_error(
        prsCs(
            bhat = rep(0.05, 5),
            R = diag(5),
            n = 100,
            nIter = 10,
            nBurnin = 10
        ),
        "keeps no draws"
    )
})

test_that("prsCs runs with valid input and estimates phi when none is given", {
    set.seed(42)
    p <- 10
    n <- 100
    bhat <- rnorm(p, sd = 0.1)
    R <- diag(p)
    for (i in seq_len(p - 1)) {
        R[i, i + 1] <- 0.3
        R[i + 1, i] <- 0.3
    }
    result <- prsCs(
        bhat = bhat,
        R = R,
        n = n,
        maf = rep(0.3, p),
        nIter = 50,
        nBurnin = 10,
        thin = 2
    )
    expect_type(result, "list")
    expect_true("betaEst" %in% names(result))
    expect_equal(length(result$betaEst), p)
    expect_true(all(is.finite(result$betaEst)))
    expect_true("phiEst" %in% names(result))
})

test_that("prsCs with explicit phi value", {
    set.seed(42)
    p <- 10
    n <- 100
    bhat <- rnorm(p, sd = 0.1)
    R <- diag(p)
    result <- prsCs(
        bhat = bhat,
        R = R,
        n = n,
        phi = 0.01,
        maf = rep(0.3, p),
        nIter = 50,
        nBurnin = 10,
        thin = 2
    )
    expect_true("phiEst" %in% names(result))
    expect_true("sigmaEst" %in% names(result))
    expect_true("psiEst" %in% names(result))
})

# ---- prsCs verbose output ----
test_that("prsCs verbose = TRUE reports its iterations and estimated phi", {
    set.seed(42)
    p <- 10
    n <- 100
    bhat <- rnorm(p, sd = 0.1)
    R <- diag(p)
    # n_iter >= 100 triggers the verbose print inside the MCMC loop; phi is
    # left NULL, so the sampler estimates it and prints the estimate.
    out <- capture.output(
        result <- prsCs(
            bhat = bhat,
            R = R,
            n = n,
            maf = rep(0.3, p),
            nIter = 110,
            nBurnin = 10,
            thin = 2,
            verbose = TRUE,
            seed = 42L
        )
    )
    expect_match(out, "Iteration +100 of 110", all = FALSE)
    expect_match(out, "Estimated global shrinkage parameter", all = FALSE)
    expect_type(result, "list")
    expect_equal(length(result$betaEst), p)
    expect_true(all(is.finite(result$betaEst)))
    expect_true(result$phiEst > 0)
})

# ---- prsCs signal recovery ----
test_that("prsCs recovers signal direction on simulated genotype data", {
    set.seed(42)
    n <- 500
    p <- 20
    X <- matrix(rbinom(n * p, 2, 0.3), nrow = n)
    beta_true <- rep(0, p)
    beta_true[c(3, 10, 15)] <- c(0.4, -0.3, 0.2)
    y <- X %*% beta_true + rnorm(n)
    bhat <- as.vector(cor(y, X))
    R <- cor(X)
    result <- prsCs(
        bhat = bhat,
        R = R,
        n = n,
        nIter = 1000,
        nBurnin = 500,
        thin = 5,
        seed = 42
    )
    expect_true("betaEst" %in% names(result))
    expect_equal(length(result$betaEst), p)
    expect_true(all(is.finite(result$betaEst)))
    # Sigma should be reasonable (near 1 for standardized data)
    expect_true(result$sigmaEst > 0.1 && result$sigmaEst < 10)
    # Correlation with truth should be positive (signal recovery)
    expect_gt(cor(result$betaEst, beta_true), 0.5)
})

# ---- prsCsWeights (wrapper) ----
test_that("prsCsWeights calls prsCs and returns betaEst", {
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
    result <- prsCsWeights(
        stat = stat,
        ldReference = .rrwLdMatrix(R),
        methodArgs = PrsCsOptions(
            maf = rep(0.3, p),
            nIter = 50,
            nBurnin = 10,
            thin = 2
        )
    )
    expect_equal(length(result), p)
    expect_true(is.numeric(result))
})

test_that("PRS-CS draws beta from the same posterior off a factor of the LD", {
    # The factor sampler (Woodbury mean, Bhattacharya-Chakraborty-Mallick
    # noise) draws differently from the dense Cholesky one, so the check is
    # the distribution: the first draw, from psi = 1 and sigma = 1, must be
    # N(A^-1 bhat, A^-1 / n) with A = R + I.
    set.seed(12)
    X <- matrix(rbinom(30 * 40, 2, 0.3), 30, 40)
    F <- pecotmr:::.ldFactorFromDosage(X)
    R <- crossprod(F) / (nrow(F) - 1)
    effects <- replace(numeric(40), c(4, 20), c(0.2, -0.15))
    bhat <- as.vector(R %*% effects) + rnorm(40, sd = 0.03)
    n <- 200L
    A <- R + diag(40)
    wantMean <- solve(A, bhat)
    wantCov <- solve(A) / n
    nDraws <- 3000L
    draws <- vapply(
        seq_len(nDraws),
        function(s) {
            pecotmr:::prsCsFactorCpp(
                1,
                0.5,
                1e-3,
                bhat,
                NULL,
                n,
                F,
                1L,
                0L,
                1L,
                FALSE,
                s,
                1,
                FALSE,
                FALSE
            )$betaEst
        },
        numeric(40)
    )
    se <- sqrt(diag(wantCov) / nDraws)
    expect_lt(max(abs(rowMeans(draws) - wantMean) / se), 4.5)
    expect_equal(cov(t(draws)), wantCov, tolerance = 0.1)
    # Each block also returns beta'A beta for sigma's update, so the chains'
    # sigma, steady across seeds, must agree too; their beta means agree to
    # Monte Carlo error.
    chain <- function(fn, ld, seed) {
        fn(
            1,
            0.5,
            1e-3,
            bhat,
            NULL,
            n,
            ld,
            4000L,
            500L,
            5L,
            FALSE,
            seed,
            1,
            FALSE,
            FALSE
        )
    }
    dense <- chain(pecotmr:::prsCsCpp, list(R), 1L)
    factor <- chain(pecotmr:::prsCsFactorCpp, F, 2L)
    expect_equal(factor$sigmaEst, dense$sigmaEst, tolerance = 0.01)
    expect_gt(cor(dense$betaEst, factor$betaEst), 0.9)
})

test_that("prsCs skips an empty LD block", {
    R <- diag(5)
    R[abs(row(R) - col(R)) == 1] <- 0.3
    run <- function(blocks) {
        pecotmr:::prsCsCpp(
            1,
            0.5,
            NULL,
            rep(0.05, 5),
            NULL,
            100L,
            blocks,
            60L,
            10L,
            5L,
            FALSE,
            7L,
            1,
            FALSE,
            FALSE
        )
    }
    expect_identical(run(list(matrix(numeric(0), 0, 0), R)), run(list(R)))
})

test_that("prsCs runs with the local shape a below and at 0.5", {
    # a - 1/2 is the order of the generalized inverse Gaussian psi is drawn
    # from: negative below a = 0.5, which the generator handles by drawing
    # the reciprocal, and zero at it, which takes its own branch.
    set.seed(6)
    R <- diag(10)
    R[abs(row(R) - col(R)) == 1] <- 0.3
    bhat <- rnorm(10, sd = 0.1)
    for (a in c(0.3, 0.5)) {
        fit <- prsCs(
            bhat,
            R,
            200,
            a = a,
            nIter = 60,
            nBurnin = 10,
            thin = 5,
            seed = 1L
        )
        expect_true(all(is.finite(fit$betaEst)), info = a)
    }
})

# ---- the shared sampler options ----

# A factor F of 30 samples' LD R over 20 variants, and statistics with one
# causal variant.
.prsCsShrinkFixture <- function() {
    set.seed(13)
    X <- matrix(rbinom(30 * 20, 2, 0.3), 30, 20)
    F <- pecotmr:::.ldFactorFromDosage(X)
    R <- crossprod(F) / (nrow(F) - 1)
    bhat <- as.vector(R %*% replace(numeric(20), 5, 0.2)) +
        rnorm(20, sd = 0.03)
    list(F = F, R = R, bhat = bhat)
}

# One PRS-CS iteration on `ld`, from psi = 1 and sigma = 1.
.prsCsOneSweep <- function(
    fn,
    ld,
    bhat,
    shrinkCorr,
    raoBlackwell,
    track = FALSE
) {
    fn(
        1,
        0.5,
        1e-3,
        bhat,
        NULL,
        200L,
        ld,
        1L,
        0L,
        1L,
        FALSE,
        9L,
        shrinkCorr,
        raoBlackwell,
        track
    )
}

test_that("PRS-CS reads the shrunk LD and reports conditional means", {
    # Rao-Blackwellized, one iteration reports its draw's conditional mean,
    # (c R + (1 - c) I + I)^-1 bhat, on the dense route and the factor's.
    f <- .prsCsShrinkFixture()
    want <- solve(0.8 * f$R + 0.2 * diag(20) + diag(20), f$bhat)
    dense <- .prsCsOneSweep(pecotmr:::prsCsCpp, list(f$R), f$bhat, 0.8, TRUE)
    factor <- .prsCsOneSweep(pecotmr:::prsCsFactorCpp, f$F, f$bhat, 0.8, TRUE)
    expect_equal(dense$betaEst, want, tolerance = 1e-10)
    expect_equal(factor$betaEst, want, tolerance = 1e-10)
    # Without, the iteration reports its draw: the mean plus noise.
    draw <- .prsCsOneSweep(pecotmr:::prsCsCpp, list(f$R), f$bhat, 0.8, FALSE)
    expect_false(isTRUE(all.equal(draw$betaEst, want)))
})

test_that("a tracked PRS-CS chain fits correlations through the shrunk LD", {
    # With one kept draw, the mean fitted correlations are c R + (1 - c) I
    # times that draw, on both routes.
    f <- .prsCsShrinkFixture()
    shrunk <- 0.8 * f$R + 0.2 * diag(20)
    for (route in list(
        list(pecotmr:::prsCsCpp, list(f$R)),
        list(pecotmr:::prsCsFactorCpp, f$F)
    )) {
        fit <- .prsCsOneSweep(route[[1]], route[[2]], f$bhat, 0.8, FALSE, TRUE)
        expect_false(fit$diverged)
        expect_equal(
            fit$corrEst,
            as.vector(shrunk %*% fit$betaEst),
            tolerance = 1e-10
        )
    }
})

# PRS-CS on a ten-variant window with one causal variant.
.prsCsOptionFit <- function(seed, raoBlackwell = FALSE, nChains = 1) {
    set.seed(14)
    R <- diag(10)
    R[abs(row(R) - col(R)) == 1] <- 0.4
    bhat <- as.vector(R %*% replace(numeric(10), 3, 0.3)) +
        rnorm(10, sd = 0.05)
    prsCs(
        bhat,
        R,
        300,
        nIter = 200,
        nBurnin = 50,
        thin = 5,
        seed = seed,
        raoBlackwell = raoBlackwell,
        nChains = nChains
    )
}

test_that("Rao-Blackwellized PRS-CS runs the same chain", {
    draws <- .prsCsOptionFit(2L)
    means <- .prsCsOptionFit(2L, raoBlackwell = TRUE)
    chain <- c("psiEst", "sigmaEst", "phiEst")
    expect_identical(means[chain], draws[chain])
    expect_false(identical(means$betaEst, draws$betaEst))
    expect_gt(cor(means$betaEst, draws$betaEst), 0.95)
})

test_that("prsCs averages several chains over the kept ones", {
    fit <- .prsCsOptionFit(4L, nChains = 3)
    expect_named(fit$chains, c("seed", "diverged", "kept"))
    expect_equal(fit$chains$seed, c(4, 5, 6))
    expect_false(any(fit$chains$diverged))
    # A tracked chain draws as the untracked one of its seed does.
    kept <- fit$chains$kept
    singles <- map(4:6, .prsCsOptionFit)[kept]
    expect_equal(fit$betaEst, reduce(map(singles, "betaEst"), `+`) / sum(kept))
    expect_null(.prsCsOptionFit(4L)$chains)
})

test_that("a diverged PRS-CS chain warns, and chains that all diverge error", {
    # Statistics of 1e-6 on one sample leave the draws near the prior's, far
    # past twice the statistics' sum of squares.
    run <- function(nChains) {
        prsCs(
            rep(1e-6, 5),
            diag(5),
            1,
            nIter = 30,
            nBurnin = 10,
            thin = 5,
            seed = 3L,
            nChains = nChains
        )
    }
    expect_warning(run(1), "prsCs: the posterior mean effects' squares")
    expect_error(run(2), "prsCs: every chain diverged")
})
