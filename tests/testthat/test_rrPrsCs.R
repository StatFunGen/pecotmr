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
                s
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
        fn(1, 0.5, 1e-3, bhat, NULL, n, ld, 4000L, 500L, 5L, FALSE, seed)
    }
    dense <- chain(pecotmr:::prsCsCpp, list(R), 1L)
    factor <- chain(pecotmr:::prsCsFactorCpp, F, 2L)
    expect_equal(factor$sigmaEst, dense$sigmaEst, tolerance = 0.01)
    expect_gt(cor(dense$betaEst, factor$betaEst), 0.9)
})
