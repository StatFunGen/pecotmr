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
        extraArgs <- extras[[nm]]
        result <- exec(
            sdpr,
            bhat = bhat,
            R = R,
            n = 100,
            !!!extraArgs,
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

test_that("sdprWeights on one study ignores the settings across studies", {
    set.seed(42)
    p <- 10
    stat <- list(
        b = rnorm(p, sd = 0.1),
        n = rep(100, p),
        variantNames = .rrwIds(p)
    )
    fit <- function(...) {
        sdprWeights(
            stat,
            .rrwLdMatrix(diag(p)),
            methodArgs = SdprOptions(
                iter = 50,
                burn = 10,
                thin = 2,
                verbose = FALSE,
                seed = 4L,
                ...
            )
        )
    }
    expect_identical(fit(sharing = "all", rho = 0.5), fit())
})

# ---- SDPRX across studies ----

# SDPR's chain across `studies` (.rrwTwoStudies()) on their LD matrices, with
# the kinds of component and rho given.
.sdprx_joint <- function(
    studies,
    kinds,
    rho,
    union = list(1:30, 1:30),
    seed = 11L
) {
    sdprxCpp(
        map(studies, "z"),
        map(studies, "R"),
        map_int(studies, "n"),
        union,
        30L,
        kinds,
        rho,
        0.1,
        1,
        50L,
        0.5,
        0.5,
        200L,
        50L,
        5L,
        FALSE,
        seed,
        FALSE,
        FALSE
    )
}

test_that("SDPR across one study is SDPR, draw for draw", {
    one <- .rrwTwoStudies()[1L]
    joint <- .sdprx_joint(one, list(1L), NULL, list(1:30))
    alone <- sdprCpp(
        one[[1L]]$z,
        list(one[[1L]]$R),
        250L,
        NULL,
        NULL,
        0.1,
        1,
        50L,
        0.5,
        0.5,
        200L,
        50L,
        5L,
        1L,
        FALSE,
        11L,
        FALSE,
        FALSE
    )
    expect_identical(as.numeric(joint$betaEst), alone$betaEst)
    expect_identical(joint$h2, alone$h2)
    # And off a factor of the LD, as sdprFactorCpp() draws, Rao-Blackwellized
    # and tracked.
    factor <- one[[1L]]$X * sqrt((250 - 1) / 250)
    jointF <- sdprxFactorCpp(
        list(one[[1L]]$z),
        list(factor),
        250L,
        list(1:30),
        30L,
        list(1L),
        NULL,
        0.1,
        1,
        50L,
        0.5,
        0.5,
        200L,
        50L,
        5L,
        FALSE,
        11L,
        TRUE,
        TRUE
    )
    aloneF <- sdprFactorCpp(
        one[[1L]]$z,
        factor,
        250L,
        NULL,
        NULL,
        0.1,
        1,
        50L,
        0.5,
        0.5,
        200L,
        50L,
        5L,
        1L,
        FALSE,
        11L,
        TRUE,
        TRUE
    )
    expect_identical(as.numeric(jointF$betaEst), aloneF$betaEst)
    expect_identical(as.numeric(jointF$corrEst), aloneF$corrEst)
})

test_that("SDPR across studies sits over the union, rho given or learned", {
    studies <- .rrwTwoStudies()
    # The second study lacks the last variant.
    studies[[2L]] <- list_modify(
        studies[[2L]],
        z = studies[[2L]]$z[-30L],
        R = studies[[2L]]$R[-30L, -30L]
    )
    union <- list(1:30, 1:29)
    kinds <- list(1L, 2L, 1:2)
    rho <- matrix(c(1, 0.8, 0.8, 1), 2L)
    given <- .sdprx_joint(studies, kinds, rho, union)
    expect_equal(dim(given$betaEst), c(30L, 2L))
    expect_equal(given$betaEst[30L, 2L], 0)
    expect_length(given$h2, 2L)
    expect_equal(given$rhoEst, rho)
    expect_false(given$diverged)
    learned <- .sdprx_joint(studies, kinds, NULL, union)
    expect_true(learned$rhoEst[1L, 2L] >= 0 && learned$rhoEst[1L, 2L] <= 0.99)
    expect_equal(diag(learned$rhoEst), c(1, 1))
})

test_that("components every study shares tie the studies' effects", {
    shared <- .sdprx_joint(
        .rrwTwoStudies(),
        list(1:2),
        matrix(c(1, 0.9, 0.9, 1), 2L)
    )
    expect_gt(cor(shared$betaEst[, 1L], shared$betaEst[, 2L]), 0.8)
})

test_that("the kinds of component follow `sharing`", {
    kinds <- pecotmr:::.sdprKinds
    expect_equal(kinds("eachOrAll", 1L), list(1L))
    expect_equal(kinds("eachOrAll", 3L), list(1L, 2L, 3L, 1:3))
    expect_equal(kinds("everySubset", 2L), kinds("eachOrAll", 2L))
    expect_equal(
        kinds("everySubset", 3L),
        list(1L, 2L, 3L, 1:2, c(1L, 3L), 2:3, 1:3)
    )
    expect_equal(kinds("all", 3L), list(1:3))
})

test_that("SDPR across studies checks its kinds, rho and inputs", {
    studies <- .rrwTwoStudies()
    expect_error(
        .sdprx_joint(studies, list(3L), NULL),
        "Kind 1's studies must lie in 1..2"
    )
    expect_error(
        .sdprx_joint(studies, list(1L, 1L), NULL),
        "Kind 2 must name distinct studies"
    )
    expect_error(
        .sdprx_joint(studies, list(1:2), matrix(c(1, 2, 2, 1), 2L)),
        "positive-definite correlation matrix"
    )
    expect_error(
        .sdprx_joint(studies, list(1:2), diag(3)),
        "positive-definite correlation matrix"
    )
    expect_error(
        .sdprx_joint(studies, list(1:2), NULL, list(1:30, 1:3)),
        "Study 2 has 3 union indices but 30 effects"
    )
    expect_error(
        .sdprx_joint(studies[1L], list(1L), NULL, list(1:30, 1:30)),
        "There are 2 union indices for 1 studies"
    )
})

test_that("SDPR weights across studies sit over the union of variants", {
    inputs <- .rrwStudyInputs(.rrwTwoStudies())
    run <- function(user) {
        pecotmr:::.sdprCoupledWeights(
            inputs$stats,
            inputs$refs,
            inputs$ids,
            user,
            "slim",
            "test"
        )
    }
    W <- run(list(
        iter = 60L,
        burn = 20L,
        thin = 2L,
        M = 20L,
        verbose = FALSE,
        seed = 3L,
        nChains = 2L,
        rho = 0.5
    ))
    expect_equal(dim(W), c(30L, 2L))
    expect_equal(rownames(W), inputs$ids)
    expect_equal(unname(W[30L, 2L]), 0)
    fit <- attr(W, "fit")
    expect_equal(fit$rhoEst, matrix(c(1, 0.5, 0.5, 1), 2L))
    expect_length(fit$h2, 2L)
    expect_equal(nrow(fit$chains), 2L)
    expect_error(
        run(list(perVariantSampleSize = rep(100, 30))),
        "describe one study's variants"
    )
    expect_error(run(list(optLlk = 2)), "describe one study's variants")
    expect_error(run(list(rho = 1)), "a single `rho` must lie in")
    expect_error(run(list(rho = diag(3))), "a 2 x 2 positive-definite")
    expect_error(run(list(sharing = "pairs")), "`sharing` must be one of")
})

test_that("SDPR skips an LD block holding no variant", {
    one <- .rrwTwoStudies()[[1L]]
    run <- function(LD) {
        sdprCpp(
            one$z,
            LD,
            250L,
            NULL,
            NULL,
            0.1,
            1,
            20L,
            0.5,
            0.5,
            60L,
            20L,
            2L,
            1L,
            FALSE,
            6L,
            FALSE,
            TRUE
        )
    }
    withEmpty <- run(list(one$R, matrix(0, 0L, 0L)))
    expect_length(withEmpty$betaEst, 30L)
    expect_true(all(is.finite(withEmpty$betaEst)))
    expect_identical(withEmpty$betaEst, run(list(one$R))$betaEst)
})

test_that("SDPR across three studies learns each pair a kind spans", {
    # Three copies of one study: their effects agree, so the learned
    # correlations run high and some grid values leave rho indefinite.
    one <- .rrwTwoStudies()[[1L]]
    studies <- rep(list(one), 3L)
    run <- function(kinds) {
        sdprxCpp(
            map(studies, "z"),
            map(studies, "R"),
            rep(250L, 3L),
            rep(list(1:30), 3L),
            30L,
            kinds,
            NULL,
            0.1,
            1,
            20L,
            0.5,
            0.5,
            200L,
            50L,
            5L,
            FALSE,
            7L,
            FALSE,
            FALSE
        )
    }
    every <- run(list(1L, 2L, 3L, 1:2, 1:3))
    rho <- every$rhoEst[upper.tri(every$rhoEst)]
    expect_true(all(rho >= 0 & rho <= 0.99))
    # A pair no kind spans keeps its correlation at 0.
    partial <- run(list(1L, 2L, 3L, 1:2))
    expect_equal(partial$rhoEst[1L, 3L], 0)
    expect_equal(partial$rhoEst[2L, 3L], 0)
})
