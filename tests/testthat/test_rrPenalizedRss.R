context("regularized_regression — penalizedRss")

# ---- penalizedRss (low-level solver) ----

test_that("penalizedRss errors on non-matrix R input", {
    expect_error(
        penalizedRss(
            bhat = rnorm(5),
            R = "not_a_matrix",
            n = 100,
            penalty = "MCP"
        ),
        "as a matrix"
    )
})

test_that("penalizedRss errors on non-positive sample size", {
    expect_error(
        penalizedRss(bhat = rnorm(5), R = diag(5), n = -1, penalty = "SCAD"),
        "valid sample size"
    )
})

test_that("penalizedRss errors on mismatched bhat and R dimensions", {
    expect_error(
        penalizedRss(bhat = rnorm(10), R = diag(5), n = 100, penalty = "MCP"),
        "bhat.*Must have length"
    )
})

test_that("penalizedRss with large lambda gives all-zero betas (MCP, SCAD)", {
    set.seed(42)
    p <- 10
    n <- 100
    bhat <- rnorm(p, sd = 0.1)
    for (pen in c("MCP", "SCAD")) {
        result <- penalizedRss(
            bhat = bhat,
            R = diag(p),
            n = n,
            penalty = pen,
            lambda = c(100)
        )
        expect_true(all(result$betaEst == 0), info = pen)
    }
})

test_that("penalizedRss with large lambda0 gives all-zero betas (L0)", {
    set.seed(42)
    p <- 10
    n <- 100
    bhat <- rnorm(p, sd = 0.1)
    result <- penalizedRss(
        bhat = bhat,
        R = diag(p),
        n = n,
        penalty = "L0",
        lambda = c(0),
        lambda0 = 1e6
    )
    expect_true(all(result$betaEst == 0))
})

test_that("penalizedRss runs with MCP and returns correct structure", {
    set.seed(42)
    p <- 10
    n <- 100
    bhat <- rnorm(p, sd = 0.1)
    R <- diag(p)
    for (i in seq_len(p - 1)) {
        R[i, i + 1] <- 0.3
        R[i + 1, i] <- 0.3
    }
    result <- penalizedRss(bhat = bhat, R = R, n = n, penalty = "MCP")
    expect_type(result, "list")
    expect_true("betaEst" %in% names(result))
    expect_equal(length(result$betaEst), p)
    expect_true(all(is.finite(result$betaEst)))
    expect_equal(nrow(result$beta), p)
    expect_equal(ncol(result$beta), 20)
    expect_true(all(result$conv %in% c(0L, 1L)))
})

test_that("penalizedRss runs with SCAD and returns correct structure", {
    set.seed(42)
    p <- 10
    n <- 100
    bhat <- rnorm(p, sd = 0.1)
    R <- diag(p)
    result <- penalizedRss(bhat = bhat, R = R, n = n, penalty = "SCAD")
    expect_type(result, "list")
    expect_equal(length(result$betaEst), p)
    expect_true(all(is.finite(result$betaEst)))
    expect_equal(nrow(result$beta), p)
})

test_that("penalizedRss runs with L0 and returns correct structure", {
    set.seed(42)
    p <- 10
    n <- 100
    bhat <- rnorm(p, sd = 0.1)
    R <- diag(p)
    result <- penalizedRss(
        bhat = bhat,
        R = R,
        n = n,
        penalty = "L0",
        lambda = c(0),
        lambda0 = 0.01
    )
    expect_type(result, "list")
    expect_equal(length(result$betaEst), p)
    expect_true(all(is.finite(result$betaEst)))
    expect_equal(nrow(result$beta), p)
})

# ---- Signal recovery on simulated summary statistics ----

test_that("mcpRssWeights / scadRssWeights recover the signal direction", {
    set.seed(2024)
    n <- 500
    p <- 20
    X <- matrix(rbinom(n * p, 2, 0.3), nrow = n)
    beta_true <- rep(0, p)
    beta_true[c(3, 10, 15)] <- c(0.4, -0.3, 0.2)
    y <- as.numeric(X %*% beta_true) + rnorm(n)
    bhat <- as.vector(cor(y, X))
    R <- cor(X)
    stat <- list(b = bhat, n = rep(n, p), variantNames = .rrwIds(p))
    fns <- list(MCP = mcpRssWeights, SCAD = scadRssWeights)
    for (pen in names(fns)) {
        w <- fns[[pen]](
            stat = stat,
            ldReference = .rrwLdMatrix(R, nSamples = n),
            s = c(0.5, 0.9)
        )
        expect_equal(length(w), p, info = pen)
        expect_true(all(is.finite(w)), info = pen)
        expect_gt(cor(w, beta_true), 0.4, label = str_c(pen, " correlation"))
    }
})

test_that("penalizedRss L0 recovers signal direction", {
    set.seed(2024)
    n <- 500
    p <- 20
    X <- matrix(rbinom(n * p, 2, 0.3), nrow = n)
    beta_true <- rep(0, p)
    beta_true[c(3, 10, 15)] <- c(0.4, -0.3, 0.2)
    y <- as.numeric(X %*% beta_true) + rnorm(n)
    bhat <- as.vector(cor(y, X))
    R <- cor(X)
    stat <- list(b = bhat, n = rep(n, p), variantNames = .rrwIds(p))
    w <- l0learnRssWeights(
        stat = stat,
        ldReference = .rrwLdMatrix(R, nSamples = n),
        penalty = "L0",
        s = c(0.5, 0.9),
        lambda0 = exp(seq(log(0.001), log(0.5), length.out = 5))
    )
    expect_equal(length(w), p)
    expect_true(all(is.finite(w)))
    expect_gt(cor(w, beta_true), 0.3)
})

# ---- Weight wrapper structure ----

test_that("scadRssWeights / mcpRssWeights return p weights with attributes", {
    set.seed(42)
    p <- 10
    n <- 100
    bhat <- rnorm(p, sd = 0.1)
    R <- diag(p)
    stat <- list(b = bhat, n = rep(n, p), variantNames = .rrwIds(p))
    fns <- list(SCAD = scadRssWeights, MCP = mcpRssWeights)
    for (pen in names(fns)) {
        w <- fns[[pen]](
            stat = stat,
            ldReference = .rrwLdMatrix(R, nSamples = n),
            s = 0.5
        )
        expect_equal(length(w), p, info = pen)
        expect_true(is.numeric(w), info = pen)
        sel <- attr(w, "penalized_rss_selection")
        expect_true(!is.null(sel), info = pen)
        expect_equal(unname(sel["penalty"]), pen, info = pen)
    }
})

test_that("l0learnRssWeights returns correct-length vector with attributes", {
    set.seed(42)
    p <- 10
    n <- 100
    bhat <- rnorm(p, sd = 0.1)
    R <- diag(p)
    stat <- list(b = bhat, n = rep(n, p), variantNames = .rrwIds(p))
    w <- l0learnRssWeights(
        stat = stat,
        ldReference = .rrwLdMatrix(R, nSamples = n),
        penalty = "L0",
        s = 0.5,
        lambda0 = c(0.01, 0.1)
    )
    expect_equal(length(w), p)
    expect_true(is.numeric(w))
    sel <- attr(w, "penalized_rss_selection")
    expect_true(!is.null(sel))
    expect_equal(unname(sel["penalty"]), "L0")
})

# ---- Elastic-net mixing (alpha) ----------------------------------------------

# `v` centered and scaled to unit variance with the 1/n denominator.
.prs_scaleN <- function(v) {
    v <- v - mean(v)
    v / sqrt(mean(v^2))
}

# Genotypes and a phenotype standardized with the 1/n denominator, so that
# R = X'X / n is their LD and z = X'y / n the marginal effects: the RSS
# objective is then exactly twice glmnet's gaussian one.
.prs_standardized <- function(n = 200L, p = 15L, seed = 7L) {
    set.seed(seed)
    X <- matrix(rbinom(n * p, 2, 0.3), n, p)
    beta <- replace(numeric(p), c(2, 4, 5), c(0.5, -0.4, 0.3))
    y <- as.numeric(X %*% beta) + rnorm(n)
    Xs <- apply(X, 2, .prs_scaleN)
    ys <- .prs_scaleN(y)
    list(
        X = Xs,
        y = ys,
        R = crossprod(Xs) / n,
        z = as.numeric(crossprod(Xs, ys)) / n,
        n = n
    )
}

# MCP and SCAD at threshold `l` and concavity `g`, for |beta| = `ab`.
.prs_mcp <- function(ab, l, g) {
    inner <- ab <= g * l
    inner * (l * ab - ab^2 / (2 * g)) + (!inner) * (g * l^2 / 2)
}

.prs_scad <- function(ab, l, g) {
    lin <- ab <= l
    mid <- !lin & ab <= g * l
    flat <- !lin & !mid
    lin *
        l *
        ab +
        mid * (2 * g * l * ab - ab^2 - l^2) / (2 * (g - 1)) +
        flat * l^2 * (g + 1) / 2
}

# The penalty penalizedRss() minimizes, written out:
# 2 P(beta; lambda * alpha) + lambda * (1 - alpha) * beta'beta.
.prs_penalty <- function(beta, lambda, alpha, penalty, gamma) {
    l1 <- lambda * alpha
    ab <- abs(beta)
    p1 <- switch(
        penalty,
        lasso = l1 * ab,
        MCP = .prs_mcp(ab, l1, gamma),
        SCAD = .prs_scad(ab, l1, gamma)
    )
    2 * sum(p1) + lambda * (1 - alpha) * sum(beta^2)
}

# The derivative of the L1-type part of the penalty at |beta| = `ab` > 0.
.prs_dpen <- function(ab, l1, penalty, gamma) {
    switch(
        penalty,
        lasso = rep(l1, length(ab)),
        MCP = pmax(l1 - ab / gamma, 0),
        SCAD = (ab <= l1) *
            l1 +
            (ab > l1) * pmax(gamma * l1 - ab, 0) / (gamma - 1)
    )
}

test_that("penalizedRss reproduces glmnet's elastic-net path on exact LD", {
    skip_if_not_installed("glmnet")
    d <- .prs_standardized()
    lam <- c(0.1, 0.05, 0.02, 0.01)
    for (a in c(1, 0.5, 0.2)) {
        rss <- penalizedRss(
            bhat = d$z * sqrt(d$n),
            R = d$R,
            n = d$n,
            penalty = "lasso",
            alpha = a,
            lambda = lam,
            thr = 1e-12,
            maxiter = 1e5
        )
        gl <- glmnet::glmnet(
            d$X,
            d$y,
            alpha = a,
            lambda = lam,
            standardize = FALSE,
            intercept = FALSE,
            thresh = 1e-14
        )
        expect_equal(
            rss$beta,
            unname(as.matrix(gl$beta)),
            tolerance = 1e-6,
            label = str_c("alpha = ", a)
        )
    }
})

test_that("penalizedRss at alpha = 0 is ridge regression", {
    d <- .prs_standardized()
    p <- ncol(d$R)
    rss <- penalizedRss(
        bhat = d$z * sqrt(d$n),
        R = d$R,
        n = d$n,
        penalty = "lasso",
        alpha = 0,
        lambda = c(1, 0.3),
        thr = 1e-13,
        maxiter = 1e5
    )
    expect_equal(
        rss$beta[, 1],
        as.numeric(solve(d$R + diag(p), d$z)),
        tolerance = 1e-8
    )
    expect_equal(
        rss$beta[, 2],
        as.numeric(solve(d$R + 0.3 * diag(p), d$z)),
        tolerance = 1e-8
    )
})

test_that("fbeta is the objective the operators minimize, at any alpha", {
    d <- .prs_standardized()
    lam <- c(0.08, 0.03, 0.01)
    for (pen in c("lasso", "MCP", "SCAD")) {
        gamma <- if (pen == "SCAD") 3.7 else 3
        res <- penalizedRss(
            bhat = d$z * sqrt(d$n),
            R = d$R,
            n = d$n,
            penalty = pen,
            alpha = 0.5,
            lambda = lam,
            gamma = gamma
        )
        for (k in seq_along(lam)) {
            b <- res$beta[, k]
            objective <- sum(b * (d$R %*% b)) -
                2 * sum(b * d$z) +
                .prs_penalty(b, lam[[k]], 0.5, pen, gamma)
            expect_equal(
                res$fbeta[[k]],
                objective,
                tolerance = 1e-10,
                label = str_c(pen, " at lambda ", lam[[k]])
            )
        }
    }
})

test_that("the updates are exact when the LD diagonal is not 1", {
    # The curvature of one coordinate is v + l2. ncvreg writes it v (1 + l2),
    # the same only at v = 1, so with a ridge term and a diagonal of 2 its
    # fixed point is not the optimum.
    d <- .prs_standardized()
    R2 <- 2 * d$R
    lam <- 0.03
    for (pen in c("lasso", "MCP", "SCAD")) {
        gamma <- if (pen == "SCAD") 3.7 else 3
        b <- penalizedRss(
            bhat = d$z * sqrt(d$n),
            R = R2,
            n = d$n,
            penalty = pen,
            alpha = 0.5,
            lambda = lam,
            gamma = gamma,
            thr = 1e-13,
            maxiter = 1e5
        )$beta[, 1]
        l1 <- lam * 0.5
        # Half the gradient of the smooth part and the ridge term.
        grad <- as.numeric(R2 %*% b) - d$z + lam * 0.5 * b
        nz <- b != 0
        expect_gt(sum(nz), 0L)
        stationary <- grad[nz] +
            sign(b[nz]) * .prs_dpen(abs(b[nz]), l1, pen, gamma)
        expect_lt(max(abs(stationary)), 1e-8, label = pen)
        expect_true(all(abs(grad[!nz]) <= l1 + 1e-9), label = pen)
    }
})

test_that("the default lambda path is scaled by the L1 share", {
    d <- .prs_standardized(p = 6L)
    bhat <- d$z * sqrt(d$n)
    base <- exp(seq(log(1e-4), log(0.1), length.out = 20))
    expect_equal(penalizedRss(bhat, d$R, d$n)$lambda, base)
    expect_equal(penalizedRss(bhat, d$R, d$n, alpha = 0.5)$lambda, base / 0.5)
    expect_equal(penalizedRss(bhat, d$R, d$n, alpha = 0)$lambda, base / 1e-3)
    # A path the caller gives is used as given.
    expect_equal(
        penalizedRss(bhat, d$R, d$n, alpha = 0.5, lambda = c(0.2, 0.1))$lambda,
        c(0.2, 0.1)
    )
})

test_that("alpha is a share in [0, 1], above 0 for MCP and SCAD", {
    d <- .prs_standardized(p = 5L)
    bhat <- d$z * sqrt(d$n)
    expect_error(penalizedRss(bhat, d$R, d$n, alpha = 1.5), "alpha.*<= 1")
    expect_error(penalizedRss(bhat, d$R, d$n, alpha = -0.1), "alpha.*>= 0")
    expect_error(
        penalizedRss(bhat, d$R, d$n, penalty = "MCP", alpha = 0),
        "above 0 for MCP"
    )
    expect_error(
        penalizedRss(bhat, d$R, d$n, penalty = "SCAD", alpha = 0),
        "above 0 for SCAD"
    )
    expect_error(penalizedRss(bhat, d$R, d$n, maxSwaps = -1), "maxSwaps")
})

test_that("lassosumRss is penalizedRss's lasso at alpha = 1, bit for bit", {
    d <- .prs_standardized()
    bhat <- d$z * sqrt(d$n)
    lassosum <- lassosumRss(bhat, d$R, d$n)
    lasso <- penalizedRss(bhat, d$R, d$n, penalty = "lasso")
    expect_identical(lassosum$beta, lasso$beta)
    expect_identical(lassosum$conv, lasso$conv)
    expect_identical(lassosum$fbeta, lasso$fbeta)
})

# The ridge's minimizer at shrinkage s and penalty l, from its normal
# equations ((1 - s) R + (s + l) I) beta = z.
.prs_ridgeExact <- function(l, R, z, s) {
    solve((1 - s) * R + (s + l) * diag(nrow(R)), z)
}

# A 30-sample reference over 50 variants: an LD of rank below its size, as
# a small reference gives it, so a factor's thin SVD spans only part of the
# variants.
.prs_lowRankLd <- function() {
    set.seed(5)
    ldFactor <- pecotmr:::.ldFactorFromDosage(matrix(
        rbinom(30 * 50, 2, 0.4),
        30,
        50
    ))
    list(ldFactor = ldFactor, R = crossprod(ldFactor) / (nrow(ldFactor) - 1))
}

test_that("a ridge is solved in closed form, on the matrix or a factor", {
    d <- .prs_lowRankLd()
    z <- as.numeric(d$R %*% rnorm(50, sd = 0.05)) + rnorm(50, sd = 0.02)
    lambda <- c(2, 0.5, 0.1)
    for (s in c(0, 0.3)) {
        exact <- exec(
            cbind,
            !!!map(lambda, .prs_ridgeExact, R = d$R, z = z, s = s)
        )
        for (ld in list(
            list(R = d$R, shrink = s),
            list(ldFactor = d$ldFactor, shrink = s)
        )) {
            fit <- pecotmr:::.rssRidgeSolve(z, lambda, ld)
            expect_equal(fit$beta, exact, tolerance = 1e-10)
            expect_identical(fit$conv, rep(1L, 3L))
            bb <- colSums(fit$beta^2)
            loss <- (1 - s) *
                colSums(fit$beta * (d$R %*% fit$beta)) +
                s * bb -
                2 * colSums(z * fit$beta)
            expect_equal(fit$loss, loss, tolerance = 1e-10)
            expect_equal(fit$fbeta, loss + lambda * bb, tolerance = 1e-10)
        }
    }
    # penalizedRss() at alpha = 0 takes it: the minimum coordinate descent
    # converges to.
    descent <- pecotmr:::penalizedRssCpp(
        z,
        list(d$R),
        0,
        lambda,
        "lasso",
        0,
        0,
        0,
        0,
        1e-13,
        100000L,
        0L
    )
    expect_equal(
        penalizedRss(z * 10, d$R, 100, alpha = 0, lambda = lambda)$beta,
        descent$beta,
        tolerance = 1e-8
    )
})

test_that("a ridge that nothing makes invertible is the minimum-norm fit", {
    d <- .prs_lowRankLd()
    z <- as.numeric(d$R %*% rnorm(50))
    ev <- eigen(d$R, symmetric = TRUE)
    kept <- ev$values > 1e-8
    minNorm <- ev$vectors[, kept] %*%
        (crossprod(ev$vectors[, kept], z) / ev$values[kept])
    # At s = lambda = 0 the LD's null space is unpenalized.
    for (ld in list(list(R = d$R), list(ldFactor = d$ldFactor))) {
        expect_equal(
            pecotmr:::.rssRidgeSolve(z, 0, ld)$beta,
            minNorm,
            tolerance = 1e-8
        )
    }
})

# ---- enetRssWeights / ridgeRssWeights ----------------------------------------

# Twelve variants, two of them causal, as summary statistics with their LD.
.prs_weightsInput <- function() {
    set.seed(42)
    p <- 12
    n <- 300
    X <- matrix(rbinom(n * p, 2, 0.3), n, p)
    y <- as.numeric(X[, c(2, 9)] %*% c(0.4, -0.3)) + rnorm(n)
    list(
        stat = list(
            b = as.vector(cor(y, X)),
            n = rep(n, p),
            variantNames = .rrwIds(p)
        ),
        ld = .rrwLdMatrix(cor(X), nSamples = n)
    )
}

test_that("enetRssWeights mixes at alpha and is lassosum at alpha = 1", {
    f <- .prs_weightsInput()
    w <- enetRssWeights(f$stat, f$ld)
    expect_length(w, 12L)
    sel <- attr(w, "penalized_rss_selection")
    expect_equal(unname(sel[c("penalty", "alpha")]), c("lasso", "0.5"))
    other <- enetRssWeights(f$stat, f$ld, alpha = 0.2)
    expect_equal(
        unname(attr(other, "penalized_rss_selection")["alpha"]),
        "0.2"
    )
    # One seed per call, so both draw the same held-out splits.
    set.seed(3)
    atOne <- enetRssWeights(f$stat, f$ld, alpha = 1)
    set.seed(3)
    expect_equal(
        as.numeric(atOne),
        as.numeric(lassosumRssWeights(f$stat, f$ld))
    )
    expect_error(enetRssWeights(f$stat, f$ld, alpha = 2), "alpha.*<= 1")
    # Its errors name it.
    expect_error(
        enetRssWeights(list(b = 1:3, n = rep(10, 3)), f$ld),
        "enetRssWeights: `stat` must carry `variantNames`"
    )
})

test_that("ridgeRssWeights fits the alpha = 0 mixing over the ridge path", {
    f <- .prs_weightsInput()
    w <- ridgeRssWeights(f$stat, f$ld, fitRetention = "full")
    expect_equal(unname(attr(w, "penalized_rss_selection")["alpha"]), "0")
    # A ridge keeps every variant, along a path from 0.1 up.
    expect_true(all(w != 0))
    fit <- attr(w, "fit")
    expect_gte(min(fit$path$lambda), 0.1)
    expect_equal(ncol(fit$beta), nrow(fit$path))
    expect_error(
        ridgeRssWeights(list(b = 1:3, n = rep(10, 3)), f$ld),
        "ridgeRssWeights: `stat` must carry `variantNames`"
    )
})

test_that("a ridge's fits at s = 1, alike but for scale, are told apart", {
    # At s = 1 each fit is c / (1 + lambda), so all tie on the selection's
    # correlation. Held-out selection lets the variance each explains set
    # the scale, so its pick does not depend on the grid's order; lassosum's
    # pseudovalidation takes the grid's first.
    f <- .prs_weightsInput()
    lambda <- c(0.1, 0.3, 1, 3, 10)
    picked <- function(lambda, selection = "heldOut") {
        set.seed(4)
        w <- ridgeRssWeights(
            f$stat,
            f$ld,
            s = 1,
            selection = selection,
            methodArgs = PenalizedRssOptions(lambda = lambda)
        )
        as.numeric(attr(w, "penalized_rss_selection")[["lambda"]])
    }
    expect_equal(picked(lambda), picked(rev(lambda)))
    skip_if_not_installed("fdrtool")
    # fdrtool warns that 12 statistics are few to fit its null on.
    expect_equal(suppressWarnings(picked(lambda, "ldQuadratic")), 0.1)
    expect_equal(suppressWarnings(picked(rev(lambda), "ldQuadratic")), 10)
})

test_that("a ridge's whole search shares one decomposition of the LD", {
    f <- .prs_weightsInput()
    source <- pecotmr:::.rssLdSource(f$stat, f$ld, "test")
    carried <- pecotmr:::.rssSolverSource(
        source,
        pecotmr:::.rssPenalty("lasso", alpha = 0)
    )
    expect_identical(carried$eigen, eigen(source$R, symmetric = TRUE))
    # Any other penalty reads the source as it is.
    expect_identical(
        pecotmr:::.rssSolverSource(
            source,
            pecotmr:::.rssPenalty("lasso", alpha = 0.5)
        ),
        source
    )
    # Four shrinkage levels, each refit on four held-out splits, and the
    # splits' draws: one decomposition.
    computed <- 0L
    decompose <- pecotmr:::.ldSourceEigen
    local_mocked_bindings(.ldSourceEigen = function(ld) {
        computed <<- computed + is.null(ld$eigen)
        decompose(ld)
    })
    ridgeRssWeights(f$stat, f$ld)
    expect_identical(computed, 1L)
})

test_that("scadRssWeights and mcpRssWeights refuse alpha = 0 before fitting", {
    f <- .prs_weightsInput()
    cnd <- tryCatch(scadRssWeights(f$stat, f$ld, alpha = 0), error = identity)
    expect_match(conditionMessage(cnd), "above 0 for SCAD")
    # Raised up front, not from inside the shrinkage grid.
    expect_null(cnd$parent)
    expect_error(mcpRssWeights(f$stat, f$ld, alpha = 0), "above 0 for MCP")
    expect_error(mcpRssWeights(f$stat, f$ld, alpha = 1.2), "alpha.*<= 1")
})

test_that("l0learnRssWeights takes its lambda path only as its own argument", {
    f <- .prs_weightsInput()
    expect_error(
        l0learnRssWeights(
            f$stat,
            f$ld,
            methodArgs = PenalizedRssOptions(lambda = 0.01)
        ),
        "lambda is supplied by pecotmr.*l0learnRssWeights' own `lambda`"
    )
})

test_that("l0learnRssWeights' L1 grid is the caller's or the penalty's", {
    lambdaFor <- pecotmr:::.l0learnRssLambda
    # L0L1 sweeps ten values on a log scale; the others have no L1 term.
    expect_equal(
        lambdaFor(NULL, "L0L1"),
        exp(seq(log(0.0001), log(0.1), length.out = 10))
    )
    expect_equal(lambdaFor(NULL, "L0"), 0)
    expect_equal(lambdaFor(NULL, "L0L2"), 0)
    expect_equal(lambdaFor(c(0.01, 0.02), "L0L1"), c(0.01, 0.02))
})

test_that(".rssPenalty checks a penalty once and fills MCP / SCAD's gamma", {
    pen <- pecotmr:::.rssPenalty
    expect_equal(pen("SCAD")$gamma, 3.7)
    expect_equal(pen("MCP")$gamma, 3)
    expect_equal(pen("lasso")$gamma, 0)
    expect_equal(pen("MCP", gamma = 5)$gamma, 5)
    expect_equal(pen("lasso", alpha = 0)$alpha, 0)
    expect_error(pen("SCAD", gamma = Inf), "gamma")
    expect_error(pen("L0", lambda0 = -1), "lambda0")
    expect_error(pen("L0L2", lambda2 = -1), "lambda2")
    expect_error(pen("L0", maxSwaps = -1), "maxSwaps")
    expect_error(pen("MCP", alpha = 0), "above 0 for MCP")
    # l0learnRssWeights checks its lambda0 grid before any fit.
    f <- .prs_weightsInput()
    expect_error(l0learnRssWeights(f$stat, f$ld, lambda0 = -1), "lambda0")
})
