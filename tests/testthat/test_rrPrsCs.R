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
        "number of rows of 'R'"
    )
})

test_that("prsCs runs successfully with valid input", {
    set.seed(42)
    p <- 10
    n <- 100
    bhat <- rnorm(p, sd = 0.1)
    R <- diag(p)
    for (i in 1:(p - 1)) {
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
})

test_that("prsCs with phi = NULL estimates phi automatically", {
    set.seed(42)
    p <- 10
    n <- 100
    bhat <- rnorm(p, sd = 0.1)
    R <- diag(p)
    result <- prsCs(
        bhat = bhat,
        R = R,
        n = n,
        phi = NULL,
        maf = rep(0.3, p),
        nIter = 50,
        nBurnin = 10,
        thin = 2
    )
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

test_that("prsCs works without maf (maf = NULL)", {
    set.seed(42)
    p <- 10
    n <- 100
    bhat <- rnorm(p, sd = 0.1)
    R <- diag(p)
    result <- prsCs(
        bhat = bhat,
        R = R,
        n = n,
        maf = NULL,
        nIter = 50,
        nBurnin = 10,
        thin = 2
    )
    expect_equal(length(result$betaEst), p)
})

# ---- prsCs verbose output ----
test_that("prsCs with verbose = TRUE produces output", {
    set.seed(42)
    p <- 10
    n <- 100
    bhat <- rnorm(p, sd = 0.1)
    R <- diag(p)
    # n_iter >= 100 triggers the verbose print inside the MCMC loop
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
    expect_type(result, "list")
    expect_equal(length(result$betaEst), p)
    expect_true(all(is.finite(result$betaEst)))
})

test_that("prsCs verbose with phi = NULL shows estimated phi", {
    set.seed(42)
    p <- 10
    n <- 100
    bhat <- rnorm(p, sd = 0.1)
    R <- diag(p)
    result <- prsCs(
        bhat = bhat,
        R = R,
        n = n,
        phi = NULL,
        maf = rep(0.3, p),
        nIter = 110,
        nBurnin = 10,
        thin = 2,
        verbose = TRUE,
        seed = 42L
    )
    expect_true("phiEst" %in% names(result))
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
    for (i in 1:(p - 1)) {
        R[i, i + 1] <- 0.3
        R[i + 1, i] <- 0.3
    }
    stat <- list(b = bhat, n = rep(n, p))
    result <- prsCsWeights(
        stat = stat,
        LD = R,
        maf = rep(0.3, p),
        nIter = 50,
        nBurnin = 10,
        thin = 2
    )
    expect_equal(length(result), p)
    expect_true(is.numeric(result))
})
