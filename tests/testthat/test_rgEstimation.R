context("rg_estimation")

# ---- fixtures ---------------------------------------------------------------

# Two traits on the same samples, from real genotypes (eqtlRegionExample):
# 300 variants, the LD from 150 of its 415 people and the statistics from the
# other 265, 30 causal variants with correlated effects (genetic correlation
# `rg`) and residual correlation `re`.
.rgBivariateFixture <- function(seed, h2, rg, re) {
    data(eqtlRegionExample, envir = environment())
    G <- eqtlRegionExample$X
    af <- colMeans(G) / 2
    G <- G[, pmin(af, 1 - af) > 0.05]
    set.seed(seed)
    win <- G[, (seed - 1) * 300 + 1:300]
    ref <- sample(nrow(win), 150)
    ldFactor <- pecotmr:::.ldFactorFromDosage(win[ref, ])
    study <- setdiff(seq_len(nrow(win)), ref)
    Xs <- scale(win[study, ])
    Xs[!is.finite(Xs)] <- 0
    idx <- sample(300, 30)
    b1 <- rnorm(30)
    b2 <- rg * b1 + sqrt(1 - rg^2) * rnorm(30)
    g1 <- as.vector(Xs[, idx] %*% b1)
    g2 <- as.vector(Xs[, idx] %*% b2)
    nS <- length(study)
    e <- matrix(rnorm(2 * nS), nS) %*% chol(matrix(c(1, re, re, 1), 2))
    y1 <- g1 * sqrt(h2[1] / var(g1)) + e[, 1] * sqrt(1 - h2[1])
    y2 <- g2 * sqrt(h2[2] / var(g2)) + e[, 2] * sqrt(1 - h2[2])
    r1 <- as.vector(cor(Xs, y1))
    r2 <- as.vector(cor(Xs, y2))
    list(
        R = crossprod(ldFactor) / (nrow(ldFactor) - 1),
        r1 = replace(r1, !is.finite(r1), 0),
        r2 = replace(r2, !is.finite(r2), 0),
        n = nS
    )
}

# A pair's inputs to the fits from a bivariate fixture.
.rgUnitData <- function(fx, noiseCorrelation = NULL, overlap = c(0, 0)) {
    M <- length(fx$r1)
    list(
        z = list(
            pecotmr:::.corToZ(fx$r1, fx$n),
            pecotmr:::.corToZ(fx$r2, fx$n)
        ),
        nVariant = list(rep(fx$n, M), rep(fx$n, M)),
        ld = list(R = fx$R, nSamples = 150L),
        overlap = overlap,
        heldIntercept = c(NA_real_, NA_real_),
        noiseCorrelation = noiseCorrelation
    )
}

# A noise-correlation matrix of `rho` between two conditions.
.rgRho <- function(rho, labels = c("a", "b")) {
    matrix(c(1, rho, rho, 1), 2, dimnames = list(labels, labels))
}

# The settings each fit is tested under.
.rgSettings <- function() {
    list(
        likelihood = methods::initialize(
            .h2LocalSettings("hdlL"),
            crossNoise = "window"
        ),
        regression = H2EstimateParam(
            criterion = "regression",
            noise = "free",
            crossNoise = "window",
            referenceNoise = "none",
            truncation = H2TruncationParam(
                of = "eigenvalues",
                maxPerSample = 0.75,
                minimum = 2
            ),
            sampleSize = "mean",
            inference = "information"
        ),
        moments = methods::initialize(
            .h2LocalSettings("lava"),
            crossNoise = "zero"
        ),
        variantRegression = H2EstimateParam(
            statistic = "variant",
            criterion = "regression",
            noise = "free",
            crossNoise = "window",
            referenceNoise = "unbiasedR2",
            sampleSize = "median",
            inference = "jackknife",
            procedures = H2RegressionParam(
                initialWeights = "moments",
                ldRedundancy = TRUE,
                jackknifeBlocks = 10
            )
        ),
        gls = H2EstimateParam(
            statistic = "variant",
            criterion = "gls",
            noise = "free",
            crossNoise = "window",
            referenceNoise = "none",
            sampleSize = "mean",
            inference = "information"
        )
    )
}

# P conditions of one gene on the same 415 samples, from eqtlRegionExample's
# genotypes with in-sample LD (an LdMatrix), as a QC'd QtlSumStats: each
# condition's genetic effects correlated `rg` with the first's.
.rgQtlFixture <- function(
    contexts = c("a", "b", "c"),
    traits = rep("g1", length(contexts)),
    studies = rep("s1", length(contexts)),
    rg = 0.7,
    seed = 11L,
    M = 150L
) {
    geno <- .rgFixtureGenotypes(M)
    set.seed(seed)
    base <- rnorm(M)
    zs <- map(
        seq_along(contexts),
        .rgQtlFixtureZ,
        X = geno$X,
        base = base,
        rg = rg
    )
    parts <- stringr::str_split_fixed(geno$ids, ":", 4)
    QtlSumStats(
        studyName = studies,
        context = contexts,
        trait = traits,
        entry = map(zs, .rgEntry, parts = parts, ids = geno$ids, n = geno$n),
        genome = "hg19",
        ldReference = .rgFixtureLd(geno),
        qcOptions = list(qc = "ok")
    )
}

# The first M MAF >= 0.05 variants of eqtlRegionExample, standardized.
.rgFixtureGenotypes <- function(M) {
    data(eqtlRegionExample, envir = environment())
    G <- eqtlRegionExample$X
    af <- colMeans(G) / 2
    G <- G[, pmin(af, 1 - af) > 0.05][, seq_len(M)]
    list(X = scale(G), ids = colnames(G), n = nrow(G), af = colMeans(G) / 2)
}

# The genotypes' in-sample LD as an LdMatrix recording `nSamples`.
.rgFixtureLd <- function(geno, nSamples = geno$n) {
    R <- crossprod(geno$X) / (geno$n - 1)
    LdMatrix(`dimnames<-`(R, list(geno$ids, geno$ids)), nSamples = nSamples)
}

# One condition's z: polygenic effects correlated `rg` with `base`, h2 0.4.
.rgQtlFixtureZ <- function(k, X, base, rg) {
    b <- rg * base + sqrt(1 - rg^2) * rnorm(length(base))
    g <- as.vector(X %*% b)
    y <- g * sqrt(0.4 / var(g)) + rnorm(nrow(X), sd = sqrt(0.6))
    r <- as.vector(cor(X, y))
    pecotmr:::.corToZ(r, nrow(X))
}

# A sumstats entry of `z` at the variants `ids`.
.rgEntry <- function(z, parts, ids, n) {
    gr <- GenomicRanges::GRanges(
        parts[, 1],
        IRanges::IRanges(as.integer(parts[, 2]), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = ids,
        A1 = parts[, 4],
        A2 = parts[, 3],
        Z = z,
        N = as.integer(n)
    )
    gr
}

# A GWAS on the QTL fixture's variants and LD reference, `...` its other
# per-study columns.
.rgGwasFixture <- function(qs, studies = "gw", seed = 3L, ...) {
    set.seed(seed)
    entries <- map(
        seq_along(studies),
        .rgGwasEntry,
        gr = qs[[1L]]
    )
    GwasSumStats(
        studyName = studies,
        entry = entries,
        genome = "hg19",
        ldReference = ldReference(qs)[[1L]],
        ...
    )
}

# @noRd
.rgGwasEntry <- function(k, gr) {
    out <- gr
    S4Vectors::mcols(out)$Z <- S4Vectors::mcols(gr)$Z * 0.5 + rnorm(length(gr))
    S4Vectors::mcols(out)$N <- 20000L
    out
}

# ---- the likelihood: HDL-L's step 2 -----------------------------------------

test_that("the likelihood's covariance step reproduces HDL.L()'s", {
    # Pinned from HDL.L()'s own statements (the HDL repository's
    # HDL/R/HDL.L.R: llfun.gcov.part.2, its covariance start, its loops over
    # starts and step sizes, the likelihood ratio and its 10,000-point grid)
    # run on these loci's rotated statistics with pecotmr's step 1, all
    # samples shared (N0 = N1 = N2), where the two scales of the cross
    # intercept coincide: covariance, cross intercept, p, rg and its interval.
    # HDL.L() reads rg and its interval off that grid, the cross intercept
    # held at the fit; pecotmr reads rg off the fit, within half the grid's
    # step of HDL.L()'s, and the interval off the profile, the cross
    # intercept refitted at each covariance, which takes in the held one's.
    want <- list(
        c(
            -0.014425885636129,
            1.48960457845523,
            0.908578366418015,
            -0.0563056305630563,
            -0.749174917491749,
            0.593759375937594
        ),
        c(
            0.0114035820030154,
            -0.126783528750789,
            0.908333738244179,
            0.0573057305730574,
            -0.635563556355635,
            0.812581258125813
        ),
        c(
            -0.184869604926025,
            -0.216426875893655,
            0.250906478847237,
            -0.635763576357636,
            -0.98019801980198,
            0.306930693069307
        )
    )
    profile <- list(
        c(-0.781554579039782, 0.716436270567631),
        c(-0.668558851829815, 0.874748323431337),
        c(-1, 0.46128643951685)
    )
    cases <- list(
        list(2L, c(0.3, 0.3), 0.6, 0.3),
        list(3L, c(0.2, 0.4), 0, 0.5),
        list(4L, c(0.4, 0.25), -0.5, 0.2)
    )
    param <- .rgSettings()$likelihood
    for (k in seq_along(cases)) {
        case <- cases[[k]]
        fx <- exec(.rgBivariateFixture, !!!case)
        ud <- .rgUnitData(fx)
        got <- pecotmr:::.rgEigenLikelihood(param, ud)[[1L]]$values
        d <- pecotmr:::.rgEigenData(param, ud)
        members <- map(
            1:2,
            pecotmr:::.rgLikelihoodMember,
            d = d,
            floor = exp(-18)
        )
        steps <- map2(
            members,
            c(NA, NA),
            pecotmr:::.rgLikelihoodStep1,
            param = param
        )
        root <- sqrt(steps[[1L]]$intercept * steps[[2L]]$intercept)
        expect_equal(
            c(got$covariance, got$noise_correlation * root, got$p_value),
            want[[k]][1:3],
            tolerance = 1e-6,
            info = str_c("case ", k)
        )
        expect_equal(got$rg, got$covariance / sqrt(got$h2_1 * got$h2_2))
        expect_lte(abs(got$rg - want[[k]][[4L]]), 1 / 9999)
        expect_lte(got$rg_lower, want[[k]][[5L]])
        expect_gte(got$rg_upper, want[[k]][[6L]])
        expect_equal(
            c(got$rg_lower, got$rg_upper),
            profile[[k]],
            tolerance = 1e-6,
            info = str_c("case ", k)
        )
        expect_identical(got$noise_source, "window")
    }
})

test_that("a window fit that fails falls back to the supplied correlation", {
    fx <- .rgBivariateFixture(2L, c(0.3, 0.3), 0.6, 0.3)
    param <- .rgSettings()$likelihood
    real <- pecotmr:::.rgCrossBest
    local_mocked_bindings(
        .rgCrossBest = function(param, ud, d, cross, mode, held) {
            if (identical(mode, "window")) {
                NULL
            } else {
                real(param, ud, d, cross, mode, held)
            }
        }
    )
    fallback <- pecotmr:::.rgEigenLikelihood(
        param,
        .rgUnitData(fx, .rgRho(0.25))
    )[[1L]]$values
    expect_identical(fallback$noise_source, "fallback")
    expect_equal(fallback$noise_correlation, 0.25)
    expect_false(is.na(fallback$covariance))
    # Without one, the pair has no estimate.
    none <- pecotmr:::.rgEigenLikelihood(param, .rgUnitData(fx))[[1L]]$values
    expect_true(is.na(none$covariance))
    expect_false(is.na(none$h2_1))
})

test_that("a pair with a heritability that is not positive has no rg", {
    expect_true(is.na(pecotmr:::.rgCorrelation(0.1, c(0.2, 0))))
    expect_true(is.na(pecotmr:::.rgCorrelation(0.1, c(0.2, NA))))
    expect_true(is.na(pecotmr:::.rgCorrelation(NA_real_, c(0.2, 0.3))))
    expect_equal(pecotmr:::.rgCorrelation(0.1, c(0.25, 0.16)), 0.5)
})

# ---- the likelihood of a group at once ---------------------------------------

# Three traits on eqtlRegionExample's 415 samples (300 variants, in-sample
# LD): genetic correlations 0.6, 0.3 and 0.5, residual 0.3, 0.2 and 0.25, h2
# 0.3 each unless `h2` gives others; a unit's data for them.
.rgTrivariateUnitData <- function(seed = 7L, h2 = rep(0.3, 3L)) {
    data(eqtlRegionExample, envir = environment())
    G <- eqtlRegionExample$X
    af <- colMeans(G) / 2
    X <- scale(G[, pmin(af, 1 - af) > 0.05][, 1:300])
    n <- nrow(X)
    set.seed(seed)
    gcor <- matrix(c(1, 0.6, 0.3, 0.6, 1, 0.5, 0.3, 0.5, 1), 3L)
    ecor <- matrix(c(1, 0.3, 0.2, 0.3, 1, 0.25, 0.2, 0.25, 1), 3L)
    b <- matrix(rnorm(3L * ncol(X)), ncol = 3L) %*% chol(gcor)
    e <- matrix(rnorm(3L * n), ncol = 3L) %*% chol(ecor)
    y <- sweep(scale(X %*% b), 2L, sqrt(h2), "*") +
        sweep(e, 2L, sqrt(1 - h2), "*")
    r <- stats::cor(X, y)
    list(
        z = map(1:3, function(k) pecotmr:::.corToZ(r[, k], n)),
        nVariant = rep(list(rep(n, ncol(X))), 3L),
        ld = list(R = crossprod(X) / (n - 1), nSamples = n),
        overlap = c(1, 1, 1),
        heldIntercept = rep(NA_real_, 3L),
        noiseCorrelation = NULL,
        phenotypicCorrelation = stats::cor(y)
    )
}

# The likelihood's settings in-sample, pair by pair or with jointFit.
.rgJointSettings <- function(jointFit = TRUE) {
    base <- methods::initialize(
        .h2LocalSettings("hdlL"),
        crossNoise = "window",
        referenceNoise = "none"
    )
    methods::initialize(
        base,
        procedures = methods::initialize(base@procedures, jointFit = jointFit)
    )
}

test_that("canonical partial correlations give a correlation matrix", {
    theta <- c(0.5, -0.3, 0.4, 0.1, 0.2, -0.6)
    W <- pecotmr:::.rgCpcFactor(theta, 4L)
    R <- tcrossprod(W)
    expect_equal(diag(R), rep(1, 4L))
    expect_equal(R[[2L, 1L]], 0.5)
    expect_gt(min(eigen(R, only.values = TRUE)$values), 0)
    expect_equal(pecotmr:::.rgCorrelationCpc(R), theta, tolerance = 1e-8)
    # The pairs' values in .rgGroupPairs()'s order land on their entries.
    pairs <- pecotmr:::.rgGroupPairs(4L)
    values <- seq_along(pairs$i) / 10
    M <- pecotmr:::.rgPairsMatrix(values, pairs, 4L)
    expect_identical(M[cbind(pairs$i, pairs$j)], values)
    expect_identical(M[cbind(pairs$j, pairs$i)], values)
    expect_identical(diag(M), rep(0, 4L))
    expect_identical(
        pecotmr:::.rgPairValues(M),
        M[cbind(c(2L, 3L, 3L, 4L, 4L, 4L), c(1L, 1L, 2L, 1L, 2L, 3L))]
    )
})

test_that("the joint likelihood is HDL.L()'s at two conditions", {
    ud <- .rgTrivariateUnitData()
    param <- .rgJointSettings()
    d <- pecotmr:::.rgEigenData(param, ud)
    floor <- param@procedures@varianceFloor
    members <- map(1:3, pecotmr:::.rgLikelihoodMember, d = d, floor = floor)
    steps <- map2(
        members,
        ud$heldIntercept,
        pecotmr:::.rgLikelihoodStep1,
        param = param
    )
    h2 <- map_dbl(steps, "h2")
    intercepts <- map_dbl(steps, "intercept")
    two <- list_modify(d, u = d$u[1:2], n = d$n[1:2], r = d$r[1:2])
    joint <- pecotmr:::.rgJointData(
        two,
        members[1:2],
        h2[1:2],
        intercepts[1:2],
        param,
        matrix(0, 2L, 2L)
    )
    cross <- pecotmr:::.rgCrossData(
        two,
        members[1:2],
        h2[1:2],
        intercepts[1:2],
        floor
    )
    # The first condition's own density holds no free parameter.
    v1 <- pmax(joint$v[, 1L], floor)
    first <- -sum(log(v1) + joint$u[, 1L]^2 / v1) / 2
    for (theta in list(c(-0.4, -0.5), c(0.2, 0.3), c(0.7, 0.3))) {
        expect_equal(
            pecotmr:::.rgJointLl(theta, joint, NULL) - first,
            pecotmr:::.rgCrossObjective(
                c(theta[[1L]] * sqrt(prod(h2[1:2])), theta[[2L]]),
                pecotmr:::.rgCrossScoringData(cross),
                FALSE
            )$value
        )
    }
    # At three, it is the multivariate normal density of each direction.
    three <- pecotmr:::.rgJointData(
        d,
        members,
        h2,
        intercepts,
        param,
        matrix(0, 3L, 3L)
    )
    theta <- c(0.5, -0.2, 0.3, 0.1, 0.25, -0.2)
    R <- tcrossprod(pecotmr:::.rgCpcFactor(theta[1:3], 3L))
    shared <- pecotmr:::.rgPairMatrix(theta[4:6], 3L) * three$crossScale
    direct <- sum(map_dbl(seq_along(d$lam), function(i) {
        sigma <- R *
            sqrt(tcrossprod(h2)) *
            three$lam2M[[i]] +
            shared * three$lam[[i]]
        diag(sigma) <- three$v[i, ]
        u <- three$u[i, ]
        -(as.numeric(determinant(sigma)$modulus) + sum(u * solve(sigma, u))) /
            2
    }))
    expect_equal(pecotmr:::.rgJointLl(theta, three, NULL), direct)
})

test_that("jointFit fits a group of three at once, and a pair as before", {
    ud <- .rgTrivariateUnitData()
    joint <- pecotmr:::.rgFit(.rgJointSettings(), ud)
    expect_length(joint, 3L)
    values <- map(joint, "values")
    rg <- map_dbl(values, "rg")
    expect_true(all(is.finite(rg)))
    # The three correlations make one correlation matrix.
    R <- pecotmr:::.rgPairsMatrix(rg, pecotmr:::.rgGroupPairs(3L), 3L) +
        diag(3L)
    expect_gt(min(eigen(R, only.values = TRUE)$values), -1e-8)
    expect_true(all(map_dbl(values, "rg_lower") <= rg))
    expect_true(all(map_dbl(values, "rg_upper") >= rg))
    expect_true(all(map_dbl(values, "p_value") > 0))
    expect_identical(map_chr(values, "noise_source"), rep("window", 3L))
    expect_equal(
        map_dbl(values, "covariance"),
        rg *
            sqrt(c(
                values[[1L]]$h2_1 * values[[1L]]$h2_2,
                values[[2L]]$h2_1 * values[[2L]]$h2_2,
                values[[3L]]$h2_1 * values[[3L]]$h2_2
            ))
    )
    # Two conditions are fitted as a pair, as HDL.L() fits them.
    two <- pecotmr:::.rgPairUnitData(ud, 1:2)
    expect_identical(
        pecotmr:::.rgFit(.rgJointSettings(), two),
        pecotmr:::.rgFit(.rgJointSettings(FALSE), two)
    )
    # Tied noise holds each pair's shared noise at Ph less the covariance.
    tied <- pecotmr:::.rgFit(
        methods::initialize(.rgJointSettings(), crossNoise = "tied"),
        ud
    )
    expect_identical(
        map_chr(map(tied, "values"), "noise_source"),
        rep("tied", 3L)
    )
})

test_that("estimateLocalRg with jointFit fits each group at once", {
    qs <- .rgQtlFixture()
    res <- estimateLocalRg(qs, param = .rgJointSettings())
    expect_length(res, 3L)
    expect_identical(res$context_1, c("a", "a", "b"))
    expect_identical(res$context_2, c("b", "c", "c"))
    expect_true(all(is.finite(res$rg)))
})

# ---- a pair given its group's other conditions -------------------------------

# The genetic covariance matrix of a group's marginal pairs: their
# covariances off the diagonal and the heritabilities on it.
.rgGroupCovariance <- function(values, P = 3L) {
    pairs <- pecotmr:::.rgGroupPairs(P)
    G <- pecotmr:::.rgPairsMatrix(map_dbl(values, "covariance"), pairs, P)
    h2 <- c(values[[1L]]$h2_1, values[[1L]]$h2_2, values[[2L]]$h2_2)
    G + diag(h2)
}

# Pair (i, j) of `G` given the third condition, and each one's r2 given both
# others, by hand.
.rgByHand <- function(G, i, j) {
    k <- setdiff(1:3, c(i, j))
    covariance <- G[i, j] - G[i, k] * G[j, k] / G[k, k]
    h2 <- c(G[i, i] - G[i, k]^2 / G[k, k], G[j, j] - G[j, k]^2 / G[k, k])
    r2 <- map_dbl(c(i, j), function(m) {
        o <- setdiff(1:3, m)
        sum(G[m, o] * solve(G[o, o], G[o, m])) / G[m, m]
    })
    list(covariance = covariance, h2 = h2, r2 = r2)
}

test_that("the moments given a group's others are LAVA's partial reading", {
    ud <- .rgTrivariateUnitData()
    moments <- methods::initialize(
        .h2LocalSettings("lava"),
        crossNoise = "tied"
    )
    marginal <- map(pecotmr:::.rgFit(moments, ud), "values")
    given <- pecotmr:::.rgFit(moments, list_assign(ud, conditional = TRUE))
    G <- .rgGroupCovariance(marginal)
    for (pair in given) {
        want <- .rgByHand(G, pair$i, pair$j)
        got <- pair$values
        info <- paste(pair$i, pair$j)
        expect_equal(got$covariance, want$covariance, info = info)
        expect_equal(c(got$h2_1, got$h2_2), want$h2, info = info)
        expect_equal(got$rg, want$covariance / sqrt(prod(want$h2)), info = info)
        expect_equal(c(got$r2_1, got$r2_2), want$r2, info = info)
        expect_true(got$p_value > 0 && got$p_value <= 1, info = info)
        expect_identical(got$noise_source, "tied", info = info)
        # The shared noise is the pair's own, not given the others.
        expect_equal(
            got$noise_correlation,
            marginal[[match(
                info,
                map_chr(given, function(x) {
                    paste(x$i, x$j)
                })
            )]]$noise_correlation,
            info = info
        )
    }
    # Given no other condition, a pair is its own, each r2 its rg squared.
    two <- pecotmr:::.rgPairUnitData(ud, 1:2)
    own <- pecotmr:::.rgFit(moments, two)[[1L]]$values
    set.seed(1L)
    alone <- pecotmr:::.rgFit(moments, list_assign(two, conditional = TRUE))
    alone <- alone[[1L]]$values
    expect_equal(alone$rg, own$rg)
    expect_equal(c(alone$r2_1, alone$r2_2), rep(own$rg^2, 2L))
})

test_that("the partial covariance's p-value is LAVA's simulation", {
    Omega <- matrix(
        c(0.010, 0.004, 0.003, 0.004, 0.012, 0.005, 0.003, 0.005, 0.009),
        3L
    )
    Sigma <- matrix(c(1, 0.3, 0.2, 0.3, 1, 0.25, 0.2, 0.25, 1), 3L) / 300
    partial <- pecotmr:::.rgPartial(Omega, 1L, 2L, 3L)
    at <- outer(1:3 == 1L, 1:3 == 2L) | outer(1:3 == 2L, 1:3 == 1L)
    null <- Omega - partial$cov * at
    sim <- list(
        top = sqrt(40) * t(pecotmr:::.rgPsdFactor(null)),
        Sigma = Sigma,
        K = 40L,
        observed = partial$cov
    )
    # Pinned to these draws; it matched a brute-force simulation of whole
    # Wishart draws under the null to its Monte Carlo error.
    set.seed(3L)
    expect_equal(
        pecotmr:::.rgPartialSimulatedP(sim, 1L, 2L, 3L, 2e4),
        0.0843853227921542,
        tolerance = 1e-6
    )
    # The null's Omega_12 is the one the third condition implies.
    expect_equal(pecotmr:::.rgPartial(null, 1L, 2L, 3L)$cov, 0)
    # Undefined where the conditioning covariance is not positive definite.
    singular <- replace(Omega, cbind(3L, 3L), -0.001)
    expect_true(is.na(pecotmr:::.rgPartial(singular, 1L, 2L, 3L)$cov))
    expect_true(is.na(pecotmr:::.rgPartialP(singular, Sigma, 40L, 1L, 2L, 3L)))
    expect_true(is.na(pecotmr:::.rgExplained(singular, 1L)))
})

test_that("the joint likelihood's pairs given their group are partial", {
    ud <- .rgTrivariateUnitData()
    param <- .rgJointSettings()
    marginal <- map(pecotmr:::.rgFit(param, ud), "values")
    given <- pecotmr:::.rgFit(param, list_assign(ud, conditional = TRUE))
    G <- .rgGroupCovariance(marginal)
    for (pair in given) {
        want <- .rgByHand(G, pair$i, pair$j)
        got <- pair$values
        info <- paste(pair$i, pair$j)
        expect_equal(got$covariance, want$covariance, info = info)
        expect_equal(c(got$h2_1, got$h2_2), want$h2, info = info)
        expect_equal(c(got$r2_1, got$r2_2), want$r2, info = info)
        expect_true(got$rg_lower <= got$rg && got$rg <= got$rg_upper)
        expect_true(got$p_value > 0 && got$p_value <= 1, info = info)
        # The likelihood gives the r2s no interval.
        expect_true(is.na(got$r2_1_lower) && is.na(got$r2_2_upper))
    }
    # A pair fitted on its own has each r2 its rg squared.
    two <- pecotmr:::.rgPairUnitData(ud, 1:2)
    alone <- pecotmr:::.rgFit(param, list_assign(two, conditional = TRUE))
    alone <- alone[[1L]]$values
    expect_equal(c(alone$r2_1, alone$r2_2), rep(alone$rg^2, 2L))
})

test_that("estimateLocalRg given the group pairs each condition with a GWAS", {
    qs <- .rgQtlFixture()
    gw <- .rgGwasFixture(qs)
    moments <- methods::initialize(
        .h2LocalSettings("lava"),
        crossNoise = "zero"
    )
    res <- suppressWarnings(
        estimateLocalRg(qs, gwas = gw, param = moments, conditional = TRUE)
    )
    # One unit, the three contexts with the GWAS: each context's pair with it.
    expect_length(res, 3L)
    expect_identical(res$context_1, c("a", "b", "c"))
    expect_identical(res$study_2, rep("gw", 3L))
    # The GWAS's r2 is its share explained by all three contexts, one value.
    expect_length(unique(res$r2_2), 1L)
    # The fixture's GWAS is context a's statistics scaled, so the contexts
    # explain all of it (above 1 under the unbounded moments): given a, b and
    # c have nothing left of it, and their partial rg is NA.
    expect_gt(res$r2_2[[1L]], 1)
    expect_true(is.finite(res$rg[[1L]]))
    expect_true(all(is.na(res$rg[2:3])))
    # Without conditional, the r2s are NA; each object records its reading.
    marginal <- suppressWarnings(
        estimateLocalRg(qs, gwas = gw, param = moments)
    )
    expect_true(all(is.na(marginal$r2_1)) && all(is.na(marginal$r2_2)))
    expect_true(isConditional(res))
    expect_false(isConditional(marginal))
})

test_that("conditional needs each group fitted at once", {
    qs <- .rgQtlFixture(contexts = c("a", "b"))
    expect_error(
        estimateLocalRg(
            qs,
            param = .rgJointSettings(FALSE),
            conditional = TRUE
        ),
        "conditional = TRUE reads each group's conditions fitted at once"
    )
    expect_error(
        estimateLocalRg(
            qs,
            param = .rgSettings()$regression,
            conditional = TRUE
        ),
        "fitted at once"
    )
    expect_error(estimateLocalRg(qs, conditional = NA), "conditional")
})

# ---- every setting -----------------------------------------------------------

test_that("each setting's step 1 is estimateLocalH2()'s fit of the condition", {
    fx <- .rgBivariateFixture(3L, c(0.3, 0.4), 0.5, 0.2)
    ud <- .rgUnitData(fx, overlap = c(0, 0))
    for (setting in names(.rgSettings())) {
        param <- .rgSettings()[[setting]]
        got <- pecotmr:::.rgFit(param, ud)[[1L]]$values
        own <- map_dbl(
            1:2,
            function(k) {
                pecotmr:::.h2LocalFit(
                    param,
                    ud$z[[k]],
                    ud$nVariant[[k]],
                    list_assign(ud$ld, referenceOverlap = 0)
                )$h2
            }
        )
        expect_equal(
            c(got$h2_1, got$h2_2),
            own,
            tolerance = 1e-8,
            info = setting
        )
        expect_identical(got$n_snps, 300L, info = setting)
    }
})

test_that("crossNoise holds the shared noise or drops it in every setting", {
    fx <- .rgBivariateFixture(4L, c(0.4, 0.3), 0.5, 0.2)
    for (setting in setdiff(names(.rgSettings()), "moments")) {
        base <- .rgSettings()[[setting]]
        zero <- pecotmr:::.rgFit(
            methods::initialize(base, crossNoise = "zero"),
            .rgUnitData(fx)
        )[[1L]]$values
        expect_identical(zero$noise_source, "zero", info = setting)
        expect_equal(zero$noise_correlation, 0, info = setting)
        supplied <- pecotmr:::.rgFit(
            methods::initialize(base, crossNoise = "supplied"),
            .rgUnitData(fx, .rgRho(0.25))
        )[[1L]]$values
        expect_identical(supplied$noise_source, "supplied", info = setting)
        expect_equal(supplied$noise_correlation, 0.25, info = setting)
    }
    moments <- methods::initialize(
        .h2LocalSettings("lava"),
        crossNoise = "supplied"
    )
    held <- pecotmr:::.rgFit(moments, .rgUnitData(fx, .rgRho(0.25)))[[1L]]
    expect_equal(held$values$noise_correlation, 0.25)
})

test_that("tied shared noise is the phenotypes' correlation less G12", {
    # A window where every setting's own intercepts are positive, so that each
    # reports a noise correlation.
    fx <- .rgBivariateFixture(6L, c(0.4, 0.3), 0.5, 0.2)
    ud <- list_assign(.rgUnitData(fx), phenotypicCorrelation = .rgRho(0.25))
    ownIntercepts <- function(param, ud) {
        map_dbl(1:2, function(k) {
            pecotmr:::.h2LocalFit(
                param,
                ud$z[[k]],
                ud$nVariant[[k]],
                list_assign(ud$ld, referenceOverlap = 0)
            )$intercept
        })
    }
    for (setting in setdiff(names(.rgSettings()), "moments")) {
        param <- methods::initialize(
            .rgSettings()[[setting]],
            crossNoise = "tied"
        )
        got <- pecotmr:::.rgFit(param, ud)[[1L]]$values
        expect_identical(got$noise_source, "tied", info = setting)
        expect_equal(
            got$noise_correlation * sqrt(prod(ownIntercepts(param, ud))),
            0.25 - got$covariance,
            tolerance = 1e-6,
            info = setting
        )
    }
    # The likelihood's tied fit is its window fit's, held on the line through
    # the window fit's covariance and shared noise.
    likelihood <- .rgSettings()$likelihood
    window <- pecotmr:::.rgFit(likelihood, .rgUnitData(fx))[[1L]]$values
    shared <- window$noise_correlation *
        sqrt(prod(ownIntercepts(likelihood, ud)))
    onLine <- pecotmr:::.rgFit(
        methods::initialize(likelihood, crossNoise = "tied"),
        list_assign(
            .rgUnitData(fx),
            phenotypicCorrelation = .rgRho(window$covariance + shared)
        )
    )[[1L]]$values
    expect_equal(onLine$covariance, window$covariance, tolerance = 1e-3)
    # The moments' shared noise per component is (Ph12 - (1 - phi) G12) /
    # (n - 1), as each condition's own is (1 - (1 - phi) h2) / (n - 1).
    moments <- methods::initialize(
        .h2LocalSettings("lava"),
        crossNoise = "tied"
    )
    got <- pecotmr:::.rgFit(moments, ud)[[1L]]$values
    phi <- pecotmr:::.lavaPhi(fx$n, 150L, 0)
    own <- (1 + c(got$h2_1, got$h2_2) * (phi - 1)) / (fx$n - 1)
    expect_equal(
        got$noise_correlation * sqrt(prod(own)),
        (0.25 - (1 - phi) * got$covariance) / (fx$n - 1)
    )
    expect_identical(got$noise_source, "tied")
})

test_that("the regression's joint sandwich holds each condition's own", {
    fx <- .rgBivariateFixture(4L, c(0.4, 0.3), 0.5, 0.2)
    param <- .rgSettings()$regression
    ud <- .rgUnitData(fx)
    d <- pecotmr:::.rgEigenData(param, ud)
    taper <- rep(1, length(d$lam))
    rows <- map(1:2, pecotmr:::.rgEigenRows, d = d, taper = taper)
    cross <- list(
        y = d$u[[1L]] * d$u[[2L]] / d$lam,
        X = cbind(d$lam / d$M, 1 / sqrt(d$n[[1L]] * d$n[[2L]])),
        taper = taper,
        redundancy = taper
    )
    fit <- pecotmr:::.rgRegressionPair(param, ud, rows, cross)
    sandwich <- pecotmr:::.rgSandwichCovariance(fit)
    own <- map2_dbl(
        rows,
        fit$steps,
        pecotmr:::.h2RegressionSe,
        param = param
    )
    expect_equal(
        sqrt(diag(sandwich$cov)[sandwich$at[1:2]]),
        own,
        tolerance = 1e-10
    )
    result <- pecotmr:::.rgEigenRegression(param, ud)[[1L]]$values
    expect_gt(result$covariance_se, 0)
    expect_gt(result$rg_se, 0)
    expect_equal(
        result$p_value,
        2 * pnorm(-abs(result$covariance / result$covariance_se))
    )
})

test_that("two references' cross covariance is R1 W R2 and the mean noise", {
    set.seed(4)
    A <- crossprod(matrix(rnorm(400), 40, 10)) / 39
    B <- crossprod(matrix(rnorm(400), 40, 10)) / 39
    w <- runif(10, 0.5, 1.5)
    model <- list(
        lds = list(list(R = A), list(R = B)),
        shared = w,
        n = c(400, 900),
        M = 10L
    )
    expect_equal(
        pecotmr:::.rgCrossModelS(model, c(0.2, 0.1)),
        sqrt(400 * 900) * 0.2 / 10 * A %*% diag(w) %*% B + 0.1 * (A + B) / 2
    )
    # On one reference, through g-LDSC's S.
    one <- list_assign(model, lds = list(list(R = A), list(R = A)))
    expect_equal(
        pecotmr:::.rgCrossModelS(one, c(0.2, 0.1)),
        sqrt(400 * 900) * 0.2 / 10 * A %*% diag(w) %*% A + 0.1 * A
    )
})

test_that("the moments fit a group's conditions at once, as LAVA does", {
    fx <- .rgBivariateFixture(4L, c(0.4, 0.3), 0.5, 0.2)
    third <- .rgBivariateFixture(4L, c(0.4, 0.3), 0.5, 0.2)
    ud <- .rgUnitData(fx)
    group <- list_modify(
        ud,
        z = c(ud$z, list(pecotmr:::.corToZ(rev(third$r2), fx$n))),
        nVariant = c(ud$nVariant, ud$nVariant[1L]),
        overlap = c(0, 0, 0)
    )
    param <- .rgSettings()$moments
    pairs <- pecotmr:::.rgEigenMoments(param, group)
    expect_length(pairs, 3L)
    expect_identical(map_int(pairs, "i"), c(1L, 1L, 2L))
    expect_identical(map_int(pairs, "j"), c(2L, 3L, 3L))
    # With one sample size and variant set, the group's (1, 2) pair is the
    # pair's own fit; each heritability is estimateLocalH2()'s.
    alone <- pecotmr:::.rgEigenMoments(param, ud)[[1L]]$values
    expect_equal(
        pairs[[1L]]$values[c("covariance", "rg", "h2_1", "h2_2")],
        alone[c("covariance", "rg", "h2_1", "h2_2")],
        tolerance = 1e-12
    )
    expect_equal(
        alone$h2_1,
        pecotmr:::.h2EigenMoments(
            param,
            ud$z[[1L]],
            ud$nVariant[[1L]],
            list_assign(ud$ld, referenceOverlap = 0)
        )$h2,
        tolerance = 1e-12
    )
})

test_that("the moments' p-value is LAVA's simulation", {
    # Pinned from a brute-force simulation of the components under no
    # covariance (2,000,000 draws, Monte Carlo standard errors 2e-4 and
    # 4e-5): the partial integration matches it.
    K <- 50
    sigma <- matrix(c(0.002, 0.0006, 0.0006, 0.0025), 2)
    withr::with_seed(5, {
        p1 <- pecotmr:::.rgMomentsP(
            matrix(c(0.004, 0.0012, 0.0012, 0.003), 2),
            sigma,
            K
        )
        p2 <- pecotmr:::.rgMomentsP(
            matrix(c(0.004, 0.002, 0.002, 0.003), 2),
            sigma,
            K
        )
    })
    expect_equal(p1, 0.066332, tolerance = 0.006 / 0.066332)
    expect_equal(p2, 0.002609, tolerance = 0.0006 / 0.002609)
    # No test without positive variances, or with perfectly shared noise.
    expect_true(is.na(pecotmr:::.rgMomentsP(diag(c(0.004, -1e-4)), sigma, K)))
    expect_true(is.na(pecotmr:::.rgMomentsP(
        diag(c(0.004, 0.003)),
        matrix(c(0.002, 0.002, 0.002, 0.002), 2),
        K
    )))
})

test_that("the moments' interval comes from noncentral Wishart draws", {
    skip_if_not_installed("matrixsampling")
    omega <- matrix(c(0.004, 0.002, 0.002, 0.003), 2)
    sigma <- matrix(c(0.002, 0.0006, 0.0006, 0.0025), 2)
    interval <- withr::with_seed(
        3,
        pecotmr:::.rgMomentsInterval(omega, sigma, 50)
    )
    rho <- 0.002 / sqrt(0.004 * 0.003)
    expect_lt(interval[[1L]], rho)
    expect_gt(interval[[2L]], rho)
    # No interval where the estimate is not a covariance matrix.
    expect_true(all(is.na(pecotmr:::.rgMomentsInterval(
        matrix(c(0.004, 0.01, 0.01, 0.003), 2),
        sigma,
        50
    ))))
})

test_that("a pair given its group has its intervals from the same draws", {
    skip_if_not_installed("matrixsampling")
    Omega <- matrix(
        c(0.010, 0.004, 0.003, 0.004, 0.012, 0.005, 0.003, 0.005, 0.009),
        3L
    )
    Sigma <- matrix(c(1, 0.3, 0.2, 0.3, 1, 0.25, 0.2, 0.25, 1), 3L) / 300
    draws <- withr::with_seed(
        3,
        pecotmr:::.rgMomentsDraws(Omega, Sigma, 40L)
    )
    drawn <- pecotmr:::.rgConditionalDrawQuantiles(
        draws,
        Sigma,
        40L,
        1L,
        2L,
        3L
    )
    partial <- pecotmr:::.rgPartial(Omega, 1L, 2L, 3L)
    rg <- partial$cov / sqrt(prod(partial$var))
    expect_lt(drawn$rg[[1L]], rg)
    expect_gt(drawn$rg[[2L]], rg)
    r2 <- c(
        pecotmr:::.rgExplained(Omega, 1L),
        pecotmr:::.rgExplained(Omega, 2L)
    )
    expect_true(all(drawn$r2Lower < r2 & r2 < drawn$r2Upper))
})

test_that("the variant GLS keeps its weight a covariance in one window", {
    fx <- .rgBivariateFixture(3L, c(0.3, 0.4), 0.5, 0.2)
    param <- .rgSettings()$gls
    ud <- .rgUnitData(fx)
    g <- pecotmr:::.rgGlsData(param, ud)
    steps <- map(1:2, pecotmr:::.rgGlsStep1, g = g, held = c(NA, NA))
    cross <- pecotmr:::.rgGlsCross(param, g, steps, NA_real_)
    weight <- steps[[1L]]$S * steps[[2L]]$S + cross$S12 * cross$S12
    values <- eigen(weight, symmetric = TRUE, only.values = TRUE)$values
    expect_gt(min(values) / max(values), -1e-10)
    result <- pecotmr:::.rgVariantGls(param, ud)[[1L]]$values
    expect_gt(result$covariance_se, 0)
})

# ---- estimateLocalRg --------------------------------------------------------

test_that("estimateLocalRg pairs each gene's contexts by default", {
    qs <- .rgQtlFixture()
    res <- suppressWarnings(estimateLocalRg(qs))
    expect_s4_class(res, "RgEstimate")
    expect_length(res, 3L)
    expect_identical(res$context_1, c("a", "a", "b"))
    expect_identical(res$context_2, c("b", "c", "c"))
    expect_identical(res$trait_1, rep("g1", 3L))
    expect_identical(h2EstimateParam(res), H2EstimateParam())
    expect_identical(res$noise_source, rep("window", 3L))
    # The context axis named, and an empty list of settings, are the default.
    expect_identical(
        suppressWarnings(estimateLocalRg(qs, jointSpecification = "context")),
        res
    )
    expect_identical(suppressWarnings(estimateLocalRg(qs, param = list())), res)
})

test_that("estimateLocalRg pairs along the axes and within the scope given", {
    qs <- .rgQtlFixture()
    scoped <- suppressWarnings(estimateLocalRg(
        qs,
        jointSpecification = list(list(
            axes = "context",
            scope = list(context = c("a", "c"))
        ))
    ))
    expect_length(scoped, 1L)
    expect_identical(c(scoped$context_1, scoped$context_2), c("a", "c"))
    traits <- .rgQtlFixture(contexts = c("a", "a"), traits = c("g1", "g2"))
    byTrait <- suppressWarnings(estimateLocalRg(
        traits,
        jointSpecification = "trait"
    ))
    expect_length(byTrait, 1L)
    expect_identical(c(byTrait$trait_1, byTrait$trait_2), c("g1", "g2"))
    # The default context axis finds no pair across two genes.
    expect_error(
        estimateLocalRg(traits),
        "no two selected rows share a group"
    )
})

test_that("estimateLocalRg's moments fit each group at once", {
    qs <- .rgQtlFixture()
    moments <- methods::initialize(
        .h2LocalSettings("lava"),
        crossNoise = "zero"
    )
    res <- suppressWarnings(
        estimateLocalRg(qs, param = moments, referenceOverlap = 1)
    )
    expect_length(res, 3L)
    expect_identical(res$noise_source, rep("zero", 3L))
    expect_true(all(is.finite(res$p_value)))
    # Each pair's heritabilities are the conditions' own under the moments.
    own <- estimateLocalH2(qs, moments, referenceOverlap = 1)
    expect_equal(res$h2_1[1:2], rep(own$h2[[1L]], 2L), tolerance = 1e-10)
})

test_that("estimateLocalRg reads the supplied noise correlation by label", {
    qs <- .rgQtlFixture()
    rho <- matrix(
        c(1, 0.2, 0.1, 0.2, 1, 0.3, 0.1, 0.3, 1),
        3,
        dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
    )
    supplied <- methods::initialize(
        H2EstimateParam(),
        crossNoise = "supplied"
    )
    res <- suppressWarnings(
        estimateLocalRg(qs, param = supplied, noiseCorrelation = rho)
    )
    expect_equal(res$noise_correlation, c(0.2, 0.1, 0.3))
    expect_identical(res$noise_source, rep("supplied", 3L))
    expect_error(
        estimateLocalRg(qs, param = supplied, noiseCorrelation = rho[1:2, 1:2]),
        "names no condition c"
    )
})

test_that("estimateLocalRg pairs each QTL row with each GWAS study", {
    qs <- .rgQtlFixture(contexts = c("a", "b"))
    gw <- .rgGwasFixture(qs, studies = c("gw1", "gw2"))
    res <- suppressWarnings(estimateLocalRg(qs, gwas = gw))
    expect_length(res, 4L)
    expect_identical(res$study_2, c("gw1", "gw1", "gw2", "gw2"))
    expect_true(all(is.na(res$context_2)) && all(is.na(res$trait_2)))
    expect_identical(res$context_1, c("a", "b", "a", "b"))
    # A GWAS side's intercept held at its genome-wide value.
    h2 <- pecotmr:::.h2EstimateNew(
        list(GenomicRanges::GRanges(), GenomicRanges::GRanges()),
        list(
            study = c("gw1", "gw2"),
            context = NA_character_,
            trait = NA_character_,
            h2 = c(0.2, 0.2),
            h2_se = c(0.01, 0.01),
            h2_lower = NA_real_,
            h2_upper = NA_real_,
            h2_liability = NA_real_,
            h2_liability_se = NA_real_,
            h2_liability_lower = NA_real_,
            h2_liability_upper = NA_real_,
            intercept = c(1.02, 0.98),
            intercept_se = c(0.01, 0.01),
            p_value = NA_real_,
            n_snps = c(100L, 100L),
            n_components = NA_integer_,
            scale = NA_real_,
            reference_overlap = NA_real_
        ),
        pecotmr:::.h2NoObjects(2L),
        GenomicRanges::GRanges(),
        .h2GenomeWideSettings("lder", FALSE, FALSE)
    )
    held <- suppressWarnings(estimateLocalRg(
        qs,
        gwas = gw,
        gwasH2Estimate = h2
    ))
    expect_length(held, 4L)
    expect_error(
        estimateLocalRg(qs, gwas = gw, gwasH2Estimate = h2[1L]),
        "no genome-wide intercept for GWAS study 'gw2'"
    )
    expect_error(
        estimateLocalRg(qs, gwas = gw, jointSpecification = "context"),
        "give one or the other"
    )
})

test_that("a case/control GWAS side is on the liability scale as well", {
    qs <- .rgQtlFixture(contexts = c("a", "b"))
    studies <- c("gw1", "gw2")
    # gw1 a case/control study, 4,000 cases of the 20,000 its N counts; gw2
    # a quantitative trait.
    gw <- .rgGwasFixture(
        qs,
        studies = studies,
        nCase = c(4000, NA),
        nControl = c(16000, NA),
        prevalence = c(0.05, NA)
    )
    res <- suppressWarnings(estimateLocalRg(qs, gwas = gw))
    # Lee et al. (2011) equation 23 at K = 0.05 and P = 0.2.
    factor <- 0.05^2 * 0.95^2 / (0.2 * 0.8 * dnorm(qnorm(0.95))^2)
    cc <- res$study_2 == "gw1"
    expect_true(all(is.finite(res$covariance)))
    expect_equal(res$h2_2_liability[cc], res$h2_2[cc] * factor)
    expect_equal(
        res$covariance_liability[cc],
        res$covariance[cc] * sqrt(factor)
    )
    expect_equal(
        res$covariance_liability_se[cc],
        res$covariance_se[cc] * sqrt(factor)
    )
    # The QTL side and the quantitative GWAS have nothing to convert.
    expect_true(all(is.na(res$h2_1_liability)))
    expect_true(all(is.na(res$h2_2_liability[!cc])))
    expect_true(all(is.na(res$covariance_liability[!cc])))
    # rg is the same on either scale.
    plain <- suppressWarnings(
        estimateLocalRg(qs, gwas = .rgGwasFixture(qs, studies = studies))
    )
    expect_identical(res$rg, plain$rg)
    expect_true(all(is.na(plain$covariance_liability)))
})

test_that("estimateLocalRg refuses settings no pair can be fitted under", {
    qs <- .rgQtlFixture(contexts = c("a", "b"))
    base <- H2EstimateParam()
    expect_error(
        estimateLocalRg(
            qs,
            param = methods::initialize(base, crossNoise = "tied")
        ),
        "tied shared noise reads each pair's phenotypic correlation"
    )
    expect_error(
        estimateLocalRg(
            qs,
            param = methods::initialize(base, crossNoise = "supplied")
        ),
        "needs `noiseCorrelation`"
    )
    expect_error(
        estimateLocalRg(
            qs,
            param = methods::initialize(
                .h2LocalSettings("fastHer"),
                crossNoise = "window"
            )
        ),
        "so it fits free noise"
    )
    expect_error(
        estimateLocalRg(qs, gwasH2Estimate = makeTestH2Estimate()),
        "so it needs `gwas`"
    )
    expect_error(
        estimateLocalRg(
            qs,
            gwas = .rgGwasFixture(qs),
            param = methods::initialize(
                .h2LocalSettings("lava"),
                crossNoise = "zero"
            ),
            gwasH2Estimate = makeTestH2Estimate()
        ),
        "the moments tie each condition's noise"
    )
    expect_error(
        estimateLocalRg(qs, noiseCorrelation = matrix(c(1, 0.2, 0.3, 1), 2)),
        "must be a symmetric correlation matrix"
    )
    expect_error(
        estimateLocalRg(qs, param = list(noise = "free")),
        "must be built with H2EstimateParam\\(\\)"
    )
    expect_error(
        estimateLocalRg(qs, param = .h2LocalSettings("greml")),
        "the phenotype statistic needs individual-level data"
    )
})

test_that("rows on different LD references are fitted each on its own", {
    qs <- .rgQtlFixture(contexts = c("a", "a"), studies = c("s1", "s2"))
    geno <- .rgFixtureGenotypes(150L)
    two <- QtlSumStats(
        studyName = c("s1", "s2"),
        context = c("a", "a"),
        trait = c("g1", "g1"),
        entry = list(qs[[1L]], qs[[2L]]),
        genome = "hg19",
        ldReference = list(.rgFixtureLd(geno), .rgFixtureLd(geno, 400L)),
        qcOptions = list(qc = "ok")
    )
    # The default fits their shared noise in the window, as on one reference.
    est <- suppressWarnings(estimateLocalRg(two, jointSpecification = "study"))
    expect_identical(h2EstimateParam(est)@crossNoise, "window")
    expect_identical(est$noise_source, "window")
    expect_true(is.finite(est$rg) && is.finite(est$noise_correlation))
    held <- suppressWarnings(estimateLocalRg(
        two,
        param = methods::initialize(H2EstimateParam(), crossNoise = "zero"),
        jointSpecification = "study"
    ))
    expect_identical(held$noise_source, "zero")
    shared <- suppressWarnings(estimateLocalRg(
        qs,
        jointSpecification = "study"
    ))
    expect_identical(c(shared$study_1, shared$study_2), c("s1", "s2"))
    expect_identical(h2EstimateParam(shared)@crossNoise, "window")
})

test_that("estimateLocalRg warns about the pairs without an rg", {
    qs <- .rgQtlFixture(contexts = c("a", "b"), rg = 0.5)
    null <- QtlSumStats(
        studyName = c("s1", "s1"),
        context = c("a", "b"),
        trait = c("g1", "g1"),
        entry = map(
            list(qs[[1L]], qs[[2L]]),
            function(gr) {
                S4Vectors::mcols(gr)$Z <- 0
                gr
            }
        ),
        genome = "hg19",
        ldReference = ldReference(qs)[[1L]],
        qcOptions = list(qc = "ok")
    )
    expect_warning(
        res <- estimateLocalRg(null),
        "1 of 1 pair\\(s\\) have no rg"
    )
    # No statistics: no step-1 fit, so nothing to bound.
    expect_true(is.na(res$rg) && is.na(res$covariance) && is.na(res$h2_1))
})

test_that("a heritability at its bound of 0 holds the covariance at 0", {
    fx <- .rgBivariateFixture(4L, c(0.3, 0.001), 0.5, 0.2)
    got <- pecotmr:::.rgFit(.rgSettings()$likelihood, .rgUnitData(fx))
    got <- got[[1L]]$values
    expect_identical(got$h2_2, 0)
    expect_identical(got$covariance, 0)
    expect_true(is.na(got$rg))
    expect_gt(got$h2_1, 0)
})

test_that("a group fits its other conditions without one at its bound of 0", {
    # The third trait's heritability sits at its bound of 0 at step 1: its
    # pairs' covariances are 0, as a pair's own fit gives, and the others are
    # fitted at once without it, here the first two as a pair, on one
    # reference and on three copies of it alike.
    ud <- .rgTrivariateUnitData(h2 = c(0.3, 0.3, 0.001))
    param <- .rgJointSettings()
    three <- list_assign(ud, ld = NULL, lds = rep(list(ud$ld), 3L), ldOf = 1:3)
    own <- pecotmr:::.rgFit(param, pecotmr:::.rgPairUnitData(ud, 1:2))
    cols <- c("rg", "rg_lower", "rg_upper", "covariance", "p_value")
    for (group in list(ud, three)) {
        got <- map(pecotmr:::.rgFit(param, group), "values")
        expect_identical(map_dbl(got[2:3], "h2_2"), c(0, 0))
        expect_identical(map_dbl(got[2:3], "covariance"), c(0, 0))
        expect_true(all(is.na(map_dbl(got[2:3], "rg"))))
        expect_equal(
            unlist(got[[1L]][cols]),
            unlist(own[[1L]]$values[cols]),
            tolerance = 1e-6
        )
    }
})

# ---- estimateLocalRg on a QtlDataset -----------------------------------------

# Three contexts of qtlDatasetExample's trait on its real genotypes (165
# donors, 200 variants): polygenic effects correlated across the contexts,
# correlated residuals, h2 0.4, and context c without its first 20 donors.
.rgDatasetFixture <- function() {
    data(qtlDatasetExample, envir = environment())
    panel <- pecotmr:::.qtlGenotypeSe(qtlDatasetExample)
    brain <- pecotmr:::.qtlContextSe(qtlDatasetExample, "brain")
    G <- scale(t(as.matrix(SummarizedExperiment::assay(panel))))
    G[!is.finite(G)] <- 0
    set.seed(4L)
    genetic <- matrix(c(1, 0.6, -0.3, 0.6, 1, 0, -0.3, 0, 1), 3L)
    residual <- matrix(c(1, 0.3, 0.2, 0.3, 1, 0.3, 0.2, 0.3, 1), 3L)
    b <- matrix(rnorm(3L * ncol(G)), ncol = 3L) %*% chol(genetic)
    e <- matrix(rnorm(3L * nrow(G)), ncol = 3L) %*% chol(residual)
    y <- scale(G %*% b) * sqrt(0.4) + e * sqrt(0.6)
    drops <- list(a = integer(0), b = integer(0), c = 1:20)
    contexts <- stats::setNames(nm = c("a", "b", "c"))
    QtlDataset(
        studyName = "s1",
        genotypes = panel,
        phenotypes = map(
            contexts,
            .rgDatasetContext,
            y = y,
            drops = drops,
            brain = brain
        )
    )
}

# One context's phenotypes: its column of `y`, missing for its dropped donors.
.rgDatasetContext <- function(cx, y, drops, brain) {
    k <- match(cx, c("a", "b", "c"))
    se <- brain[, rownames(y)]
    SummarizedExperiment::assay(se) <- matrix(
        replace(y[, k], drops[[cx]], NA),
        1L,
        dimnames = dimnames(se)
    )
    se
}

# A GWAS of 20,000 on the fixture's genotyped variants whose local signal is
# context a's, coded the other way round: each variant's alleles swapped and
# its z negated, so that matching flips them back, and drops the
# strand-ambiguous ones, whose strand a swap leaves undecided.
.rgDatasetGwas <- function(qd) {
    panel <- pecotmr:::.qtlGenotypeSe(qd)
    rr <- SummarizedExperiment::rowRanges(panel)
    G <- scale(t(as.matrix(SummarizedExperiment::assay(panel))))
    a <- pecotmr:::.qtlContextSe(qd, "a")
    y <- SummarizedExperiment::assay(a)[1L, rownames(G)]
    r <- suppressWarnings(as.vector(stats::cor(G, y)))
    set.seed(5L)
    z <- sqrt(20000) * 0.8 * replace(r, !is.finite(r), 0) + rnorm(length(r))
    gr <- GenomicRanges::GRanges(
        GenomicRanges::seqnames(rr),
        IRanges::IRanges(GenomicRanges::start(rr), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = paste(
            GenomicRanges::seqnames(rr),
            GenomicRanges::start(rr),
            rr$A1,
            rr$A2,
            sep = ":"
        ),
        A1 = rr$A2,
        A2 = rr$A1,
        Z = -z,
        N = 20000L
    )
    GwasSumStats(studyName = "gw", entry = list(gr), genome = "hg19")
}

test_that("a QtlDataset pair is fitted on the donors both conditions have", {
    skip_if_not_installed("PlinkMatrix")
    qd <- .rgDatasetFixture()
    param <- methods::initialize(
        .rgSettings()$likelihood,
        referenceNoise = "none"
    )
    res <- estimateLocalRg(qd, param = param, cisWindow = 1e6)
    expect_s4_class(res, "RgEstimate")
    expect_identical(res$context_1, c("a", "a", "b"))
    expect_identical(res$context_2, c("b", "c", "c"))
    expect_identical(unique(c(res$trait_1, res$trait_2)), "ENSG_example")
    expect_true(all(is.finite(res$covariance)))
    # The a-b pair is on all 165 donors, so it is the same fitted alone.
    alone <- estimateLocalRg(
        qd,
        param = param,
        contexts = c("a", "b"),
        cisWindow = 1e6
    )
    expect_identical(alone$covariance, res$covariance[[1L]])
    # Its statistics are each condition's marginal correlations with the
    # window's genotypes on those donors, and those genotypes its LD.
    pf <- pecotmr:::.qtlResolvePhenotypeFilter(
        PhenotypeFilterParam(),
        naDefault = "drop",
        naAllowed = c("drop", "impute"),
        caller = "test"
    )
    y <- map(c("a", "b"), function(cx) {
        pecotmr:::.fmResidPheno(
            qd,
            ResidualizationParam(),
            contexts = cx,
            traitId = "ENSG_example",
            phenotypeFilter = pf
        )
    })
    donors <- intersect(rownames(y[[1L]]), rownames(y[[2L]]))
    X <- pecotmr:::.fmResidGeno(
        qd,
        ResidualizationParam(),
        pecotmr:::.qtlResolveFilter(GenotypeFilterParam()),
        contexts = c("a", "b"),
        traitId = "ENSG_example",
        cisWindow = 1e6,
        samples = donors
    )[donors, ]
    Xs <- scale(X[, matrixStats::colSds(X) > 0])
    Y <- cbind(y[[1L]][donors, 1L], y[[2L]][donors, 1L])
    r <- stats::cor(Xs, Y)
    n <- length(donors)
    M <- ncol(Xs)
    want <- pecotmr:::.rgFit(
        param,
        list(
            z = list(
                pecotmr:::.corToZ(r[, 1L], n),
                pecotmr:::.corToZ(r[, 2L], n)
            ),
            nVariant = list(rep(n, M), rep(n, M)),
            ld = list(ldFactor = Xs, nSamples = n),
            overlap = c(1, 1),
            heldIntercept = c(NA_real_, NA_real_),
            noiseCorrelation = NULL,
            phenotypicCorrelation = stats::cor(Y)
        )
    )[[1L]]$values
    expect_identical(n, 165L)
    expect_identical(res$n_snps[[1L]], M)
    expect_equal(res$covariance[[1L]], want$covariance)
    expect_equal(res$rg[[1L]], want$rg)
    expect_equal(res$noise_correlation[[1L]], want$noise_correlation)
})

test_that("a QtlDataset group fit is on the donors all its conditions have", {
    skip_if_not_installed("PlinkMatrix")
    qd <- .rgDatasetFixture()
    moments <- methods::initialize(
        .h2LocalSettings("lava"),
        crossNoise = "tied"
    )
    group <- suppressWarnings(
        estimateLocalRg(qd, param = moments, cisWindow = 1e6)
    )
    expect_identical(group$noise_source, rep("tied", 3L))
    expect_true(all(is.finite(group$covariance)))
    # c lacks 20 donors, so the group's a-b pair is not the pair alone.
    alone <- suppressWarnings(estimateLocalRg(
        qd,
        param = moments,
        contexts = c("a", "b"),
        cisWindow = 1e6
    ))
    expect_false(isTRUE(all.equal(alone$covariance, group$covariance[[1L]])))
})

test_that("estimateLocalRg pairs each QtlDataset condition with a GWAS", {
    skip_if_not_installed("PlinkMatrix")
    qd <- .rgDatasetFixture()
    gw <- .rgDatasetGwas(qd)
    param <- methods::initialize(
        .rgSettings()$likelihood,
        crossNoise = "zero"
    )
    res <- estimateLocalRg(qd, gwas = gw, param = param, cisWindow = 1e6)
    expect_length(res, 3L)
    expect_identical(res$context_1, c("a", "b", "c"))
    expect_identical(res$study_2, rep("gw", 3L))
    expect_true(all(is.na(res$context_2)) && all(is.na(res$trait_2)))
    expect_identical(res$noise_source, rep("zero", 3L))
    # The GWAS's swapped alleles are flipped back, so it correlates with a.
    expect_gt(res$rg[[1L]], 0.5)
    # Its strand-ambiguous variants, whose strand a swap leaves undecided, are
    # dropped.
    rr <- SummarizedExperiment::rowRanges(pecotmr:::.qtlGenotypeSe(qd))
    ambiguous <- paste0(rr$A1, rr$A2) %in% c("AT", "TA", "CG", "GC")
    expect_true(all(res$n_snps <= sum(!ambiguous)))
    expect_error(
        estimateLocalRg(
            qd,
            gwas = gw,
            param = methods::initialize(param, crossNoise = "tied")
        ),
        "a GWAS side does not have"
    )
})

test_that("a QtlDataset's GWAS side reads its own reference", {
    skip_if_not_installed("PlinkMatrix")
    qd <- .rgDatasetFixture()
    gw <- .rgDatasetGwas(qd)
    # The GWAS on a panel of 100 of the donors, named as the GWAS names its
    # variants (every variant's alleles swapped, which leaves LD as it is): a
    # reference of its own. A variant constant on those donors is
    # independent.
    panel <- pecotmr:::.qtlGenotypeSe(qd)
    G <- scale(t(as.matrix(SummarizedExperiment::assay(panel)))[1:100, ])
    G[!is.finite(G)] <- 0
    R <- stats::cov2cor(
        crossprod(G) + diag(as.numeric(colSums(G^2) == 0), ncol(G))
    )
    ids <- S4Vectors::mcols(gw[[1L]])$SNP
    held <- GwasSumStats(
        studyName = "gw",
        entry = list(gw[[1L]]),
        genome = "hg19",
        ldReference = LdMatrix(`dimnames<-`(R, list(ids, ids)), nSamples = 100L)
    )
    est <- suppressWarnings(estimateLocalRg(
        qd,
        gwas = held,
        contexts = "a",
        cisWindow = 1e6
    ))
    expect_identical(h2EstimateParam(est)@crossNoise, "window")
    expect_identical(est$noise_source, "window")
    # Without one, the GWAS side reads the genotypes as before.
    plain <- suppressWarnings(estimateLocalRg(
        qd,
        gwas = gw,
        contexts = "a",
        cisWindow = 1e6
    ))
    expect_identical(h2EstimateParam(plain)@crossNoise, "window")
})

test_that("a QtlDataset's conditions pair with a GWAS given their group", {
    skip_if_not_installed("PlinkMatrix")
    qd <- .rgDatasetFixture()
    gw <- .rgDatasetGwas(qd)
    moments <- methods::initialize(
        .h2LocalSettings("lava"),
        crossNoise = "zero"
    )
    res <- suppressWarnings(estimateLocalRg(
        qd,
        gwas = gw,
        param = moments,
        cisWindow = 1e6,
        conditional = TRUE
    ))
    expect_length(res, 3L)
    expect_identical(res$context_1, c("a", "b", "c"))
    expect_identical(res$study_2, rep("gw", 3L))
    # One unit, the GWAS's r2 one value (NA here: the fixture's GWAS is
    # context a's statistics scaled, so the four conditions' Omega is not
    # positive definite).
    expect_length(unique(res$r2_2), 1L)
    # The group is on the donors all three contexts have (c lacks 20), so
    # every pair has one K; alone, each context is on its own donors.
    marginal <- suppressWarnings(estimateLocalRg(
        qd,
        gwas = gw,
        param = moments,
        cisWindow = 1e6
    ))
    expect_length(unique(res$n_components), 1L)
    expect_identical(res$n_components, rep(marginal$n_components[[3L]], 3L))
    expect_gt(marginal$n_components[[1L]], marginal$n_components[[3L]])
    expect_error(
        estimateLocalRg(
            qd,
            param = .rgJointSettings(FALSE),
            conditional = TRUE
        ),
        "fitted at once"
    )
})

test_that("estimateLocalRg on a QtlDataset refuses what it cannot fit", {
    data(qtlDatasetExample, envir = environment())
    expect_error(
        estimateLocalRg(
            qtlDatasetExample,
            param = methods::initialize(
                .h2LocalSettings("greml"),
                crossNoise = "zero"
            )
        ),
        "REML has no bivariate fit here"
    )
    expect_error(
        estimateLocalRg(qtlDatasetExample, gwas = H2EstimateParam()),
        "gwas"
    )
    expect_error(
        estimateLocalRg(qtlDatasetExample, jointSpecification = "study"),
        "the study axis pairs nothing"
    )
    expect_error(
        estimateLocalRg(qtlDatasetExample),
        "no two selected conditions share a group"
    )
    two <- GenomicRanges::GRanges(c("chr22:1-1000", "chr22:5000-6000"))
    expect_error(
        estimateLocalRg(qtlDatasetExample, region = two),
        "estimate one region at a time"
    )
})

# ---- summarizeLocalRg --------------------------------------------------------

# An RgEstimate of one pair of contexts (a, b) over `n` genes, and of (a, c)
# over the first three of them, from per-gene covariances, heritabilities and
# log-likelihood curves.
.rgSummaryFixture <- function(
    covariance,
    h1,
    h2,
    curves = rep(list(NULL), length(covariance)),
    param = methods::initialize(.h2LocalSettings("lava"), crossNoise = "zero"),
    conditional = FALSE
) {
    n <- length(covariance)
    genes <- paste0("g", seq_len(n))
    at <- c(seq_len(n), 1:3)
    rg <- covariance / sqrt(pmax(h1 * h2, 0))
    window <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 1000))
    na <- rep(NA_real_, length(at))
    pecotmr:::.rgEstimateNew(
        rep(list(window), length(at)),
        list(
            study_1 = rep("s1", length(at)),
            context_1 = rep("a", length(at)),
            trait_1 = genes[at],
            study_2 = rep("s1", length(at)),
            context_2 = rep(c("b", "c"), c(n, 3L)),
            trait_2 = genes[at],
            rg = ifelse(h1[at] > 0 & h2[at] > 0, rg[at], NA_real_),
            rg_se = na,
            rg_lower = na,
            rg_upper = na,
            covariance = covariance[at],
            covariance_se = na,
            h2_1 = h1[at],
            h2_2 = h2[at],
            covariance_liability = na,
            covariance_liability_se = na,
            h2_1_liability = na,
            h2_2_liability = na,
            noise_correlation = na,
            noise_source = rep("zero", length(at)),
            p_value = na,
            r2_1 = na,
            r2_1_lower = na,
            r2_1_upper = na,
            r2_2 = na,
            r2_2_lower = na,
            r2_2_upper = na,
            n_snps = rep(100L, length(at)),
            n_components = rep(40L, length(at)),
            rg_loglik = curves[at]
        ),
        window,
        param,
        conditional
    )
}

# A gene's log-likelihood curve -c (r - mu)^2 / 2 on the RgEstimate's grid.
.rgQuadraticCurve <- function(mu, c) {
    grid <- seq(-1, 1, length.out = 101L)
    -c * (grid - mu)^2 / 2
}

test_that("the pooled summary is the ratio of sums, with a gene jackknife", {
    covariance <- c(0.05, 0.02, 0.08, -0.01, 0.03)
    h1 <- c(0.1, 0.05, 0.2, 0.02, 0.1)
    h2 <- c(0.08, 0.04, 0.15, 0.05, 0.06)
    out <- summarizeLocalRg(.rgSummaryFixture(covariance, h1, h2), "pooled")
    expect_s4_class(out, "DataFrame")
    expect_identical(nrow(out), 2L)
    expect_identical(out$context_2, c("b", "c"))
    expect_identical(out$summary, rep("pooled", 2L))
    ratio <- function(k) sum(covariance[k]) / sqrt(sum(h1[k]) * sum(h2[k]))
    expect_equal(out$rg[[1L]], ratio(1:5))
    expect_equal(out$rg[[2L]], ratio(1:3))
    loo <- map_dbl(1:5, function(g) ratio(setdiff(1:5, g)))
    expect_equal(out$rg_se[[1L]], sqrt(4 / 5 * sum((loo - mean(loo))^2)))
    expect_identical(out$n_genes, c(5L, 3L))
    expect_equal(
        out$rg_lower[[1L]],
        out$rg[[1L]] - qnorm(0.975) * out$rg_se[[1L]]
    )
    expect_true(all(is.na(out$mash_rg)))
    expect_false(S4Vectors::metadata(out)$conditional)
})

test_that("the mean summary averages the genes whose heritabilities pass", {
    covariance <- c(0.05, 0.02, 0.08, -0.01, 0.03)
    h1 <- c(0.1, 0.05, 0.2, 0.02, 0.1)
    h2 <- c(0.08, 0.04, 0.15, 0.05, 0.06)
    x <- .rgSummaryFixture(covariance, h1, h2)
    rg <- covariance / sqrt(pmax(h1 * h2, 0))
    every <- summarizeLocalRg(x, "mean")
    expect_equal(every$rg[[1L]], mean(rg))
    expect_equal(every$rg_se[[1L]], sd(rg) / sqrt(5))
    strong <- summarizeLocalRg(x, "mean", minH2 = 0.045)
    kept <- h1 > 0.045 & h2 > 0.045
    expect_equal(strong$rg[[1L]], mean(rg[kept]))
    expect_identical(strong$n_genes[[1L]], sum(kept))
    # A gene whose heritability is not positive has no rg to average.
    negative <- .rgSummaryFixture(covariance, replace(h1, 2L, -0.01), h2)
    expect_identical(summarizeLocalRg(negative, "mean")$n_genes[[1L]], 4L)
    expect_error(summarizeLocalRg(x, "average"), "must be one of")
})

test_that("the shared summary maximizes the genes' summed log-likelihoods", {
    mu <- c(0.3, 0.5, 0.7, 0.4, 0.6)
    curvature <- c(10, 20, 15, 5, 8)
    curves <- map2(mu, curvature, .rgQuadraticCurve)
    x <- .rgSummaryFixture(
        rep(0.01, 5),
        rep(0.1, 5),
        rep(0.1, 5),
        curves
    )
    out <- summarizeLocalRg(x, "shared")
    total <- sum(curvature)
    shared <- sum(curvature * mu) / total
    expect_equal(out$rg[[1L]], shared)
    # The sandwich: each gene's slope at the maximum, over the curvature.
    expect_equal(
        out$rg_se[[1L]],
        sqrt(sum((curvature * (shared - mu))^2)) / total
    )
    half <- sqrt(qchisq(0.95, 1) / total)
    expect_equal(out$rg_lower[[1L]], shared - half, tolerance = 1e-3)
    expect_equal(out$rg_upper[[1L]], shared + half, tolerance = 1e-3)
    expect_equal(
        out$p_value[[1L]],
        pchisq(total * shared^2, 1, lower.tail = FALSE)
    )
    # No curves, as under the moments: nothing to add up.
    none <- summarizeLocalRg(.rgSummaryFixture(mu, mu, mu), "shared")
    expect_true(is.na(none$rg[[1L]]))
    expect_identical(none$n_genes[[1L]], 0L)
})

test_that("mash's prior correlation cross-checks each pair of contexts", {
    x <- .rgSummaryFixture(c(0.05, 0.02, 0.08), rep(0.1, 3), rep(0.1, 3))
    conds <- c("a", "b", "c")
    U1 <- matrix(
        c(1, 0.5, 0.2, 0.5, 1, 0.1, 0.2, 0.1, 1),
        3,
        dimnames = list(conds, conds)
    )
    U2 <- diag(3)
    dimnames(U2) <- list(conds, conds)
    model <- list(
        fitted_g = list(
            pi = c(0.5, 0.1, 0.2, 0.15, 0.05),
            Ulist = list(U1 = U1, U2 = U2),
            grid = c(0.5, 2),
            usepointmass = TRUE
        )
    )
    out <- summarizeLocalRg(x, "pooled", mashModel = model)
    weights <- c(0.1 * 0.25 + 0.15 * 4, 0.2 * 0.25 + 0.05 * 4)
    implied <- stats::cov2cor(weights[[1L]] * U1 + weights[[2L]] * U2)
    expect_equal(out$mash_rg, c(implied[["a", "b"]], implied[["a", "c"]]))
    expect_error(
        summarizeLocalRg(
            x,
            mashModel = list(fitted_g = list(Ulist = list(diag(3))))
        ),
        "must be a fitted mash model"
    )
    given <- .rgSummaryFixture(
        c(0.05, 0.02, 0.08),
        rep(0.1, 3),
        rep(0.1, 3),
        conditional = TRUE
    )
    expect_true(
        S4Vectors::metadata(summarizeLocalRg(given, "mean"))$conditional
    )
    expect_error(
        summarizeLocalRg(given, mashModel = model),
        "does not cross-check a conditional RgEstimate"
    )
})

test_that("estimateLocalRg keeps each pair's curve for the shared summary", {
    qs <- .rgQtlFixture(contexts = c("a", "b"))
    res <- estimateLocalRg(qs)
    curve <- res$rg_loglik[[1L]]
    expect_length(curve, 101L)
    expect_equal(max(curve), 0)
    # One gene's shared rg is its own, read off its curve.
    shared <- summarizeLocalRg(res, "shared")
    expect_equal(shared$rg[[1L]], res$rg[[1L]], tolerance = 0.01)
    # The likelihood bounds each pair, which a pooled ratio reads low.
    expect_warning(summarizeLocalRg(res, "pooled"), "pooled ratio reads low")
    # The moments keep none.
    moments <- methods::initialize(
        .h2LocalSettings("lava"),
        crossNoise = "zero"
    )
    none <- suppressWarnings(estimateLocalRg(qs, param = moments))
    expect_null(none$rg_loglik[[1L]])
})


# ---- conditions on different LD references -----------------------------------

# Each setting with no shared noise; the variant regression without the
# unbiased r^2, which a cross-population LD score from two independent panels
# does not take.
.rgZeroSettings <- function() {
    s <- .rgSettings()
    zero <- list(
        likelihood = s$likelihood,
        regression = s$regression,
        moments = s$moments,
        variantRegression = methods::initialize(
            s$variantRegression,
            referenceNoise = "none",
            inference = "information"
        ),
        gls = s$gls
    )
    map(zero, methods::initialize, crossNoise = "zero")
}

# A unit's data on the first `M` variants, with its one LD given as one
# reference per condition.
.rgAsTwoLd <- function(ud, M = 120L) {
    at <- seq_len(M)
    one <- list(
        z = map(ud$z, `[`, at),
        nVariant = map(ud$nVariant, `[`, at),
        ld = list(R = ud$ld$R[at, at], nSamples = ud$ld$nSamples),
        overlap = ud$overlap,
        heldIntercept = ud$heldIntercept,
        noiseCorrelation = NULL
    )
    two <- c(
        one[setdiff(names(one), "ld")],
        list(
            ld = NULL,
            lds = rep(list(one$ld), length(one$z)),
            ldOf = seq_along(one$z)
        )
    )
    list(one = one, two = two)
}

test_that("each setting on two references is its own fit when they agree", {
    # The two-reference forms reduce to the one-reference fits where the two
    # references are one LD. With no shared noise: the likelihood's reduced
    # pairs to HDL.L()'s step 2, the eigen regression's to its normal
    # equations and sandwich, the variant forms' cross scores to the LD
    # scores. With shared noise, whose shape (R1 + R2) / 2 is then that LD,
    # each to its fit of it. The moments' estimate is LAVA's, its inference
    # simulated.
    pair <- .rgAsTwoLd(pecotmr:::.rgPairUnitData(.rgTrivariateUnitData(), 1:2))
    cols <- c("rg", "rg_lower", "rg_upper", "covariance", "h2_1", "h2_2")
    for (mode in c("zero", "window", "supplied")) {
        rho <- if (mode == "supplied") .rgRho(0.2)
        # The moments have no window fit of the shared noise.
        skipped <- if (mode == "window") "moments"
        for (name in setdiff(names(.rgZeroSettings()), skipped)) {
            param <- methods::initialize(
                .rgZeroSettings()[[name]],
                crossNoise = mode
            )
            one <- pecotmr:::.rgFit(
                param,
                list_assign(pair$one, noiseCorrelation = rho)
            )[[1L]]$values
            two <- pecotmr:::.rgFit(
                param,
                list_assign(pair$two, noiseCorrelation = rho)
            )[[1L]]$values
            compared <- if (name == "moments") {
                c("rg", "covariance", "h2_1", "h2_2", "noise_correlation")
            } else {
                c(cols, "rg_se", "p_value", "noise_correlation")
            }
            expect_equal(
                unlist(two[compared]),
                unlist(one[compared]),
                tolerance = 1e-5,
                info = paste(name, mode)
            )
            expect_identical(two$noise_source, mode, info = paste(name, mode))
        }
    }
})

test_that("the variant regression's information reads LD between variants", {
    # Its three regressions' sandwich takes the responses' covariances under
    # the model over every pair of variants (.rgLdMeat()). With flat weights
    # the covariance's estimate is linear in the cross-products, so the
    # sandwich at the model's own covariances is its variance, which draws
    # from that model reproduce.
    R <- .rgBivariateFixture(4L, c(0.4, 0.3), 0.5, 0.2)$R
    M <- nrow(R)
    n <- 500
    RR <- n / M * R %*% R
    S <- list(s11 = 0.3 * RR + R, s22 = 0.2 * RR + R, s12 = 0.1 * RR + 0.3 * R)
    X <- cbind(n * rowSums(R^2) / M, 1)
    bread <- pecotmr:::.rgBlockDiagonal(rep(list(crossprod(X)), 3L))
    meat <- pecotmr:::.rgLdMeat(list(X, X, X), S)
    variance <- solve(bread, t(solve(bread, meat)))[[5L, 5L]]
    a <- solve(crossprod(X), t(X))[1L, ]
    cross <- S$s11 * S$s22 + S$s12 * t(S$s12)
    expect_equal(variance, sum(a * (cross %*% a)))
    sigma <- rbind(cbind(S$s11, S$s12), cbind(t(S$s12), S$s22))
    e <- eigen(sigma, symmetric = TRUE)
    root <- e$vectors %*% (sqrt(pmax(e$values, 0)) * t(e$vectors))
    flat <- H2EstimateParam(
        statistic = "variant",
        criterion = "regression",
        noise = "free",
        crossNoise = "window",
        referenceNoise = "none",
        inference = "information",
        procedures = H2RegressionParam(weights = "fixed")
    )
    set.seed(21)
    estimates <- map_dbl(seq_len(400), function(i) {
        z <- pecotmr:::.corToZ(as.vector(root %*% rnorm(2L * M)) / sqrt(n), n)
        ud <- list(
            z = list(z[seq_len(M)], z[M + seq_len(M)]),
            nVariant = list(rep(n, M), rep(n, M)),
            ld = list(R = R, nSamples = 150L),
            overlap = c(0, 0),
            heldIntercept = c(NA_real_, NA_real_),
            noiseCorrelation = NULL
        )
        pecotmr:::.rgFit(flat, ud)[[1L]]$values$covariance
    })
    expect_equal(sd(estimates), sqrt(variance), tolerance = 0.1)
    # Two copies of one reference read as the one.
    pair <- .rgAsTwoLd(.rgUnitData(.rgBivariateFixture(
        4L,
        c(0.4, 0.3),
        0.5,
        0.2
    )))
    information <- methods::initialize(flat, procedures = H2RegressionParam())
    expect_equal(
        pecotmr:::.rgFit(information, pair$two)[[1L]]$values$covariance_se,
        pecotmr:::.rgFit(information, pair$one)[[1L]]$values$covariance_se,
        tolerance = 1e-8
    )
})

test_that("the moments read the effects as fixed or random", {
    # On one reference both readings give LAVA's estimate; random effects
    # draw the components' genetic variance too, as on two references, where
    # the reading is the same when they agree. Fixed effects on two
    # references are LAVA's reading, their p-value and interval within the
    # draws' error of its simulation's (at these seeds, the p-values within
    # 35% and the bounds within 0.005).
    pair <- .rgAsTwoLd(pecotmr:::.rgPairUnitData(.rgTrivariateUnitData(), 1:2))
    estimates <- c("rg", "covariance", "h2_1", "h2_2", "noise_correlation")
    inference <- c("rg_lower", "rg_upper", "p_value")
    for (mode in c("zero", "supplied")) {
        rho <- if (mode == "supplied") .rgRho(0.2)
        one <- list_assign(pair$one, noiseCorrelation = rho)
        two <- list_assign(pair$two, noiseCorrelation = rho)
        fixed <- methods::initialize(
            .rgZeroSettings()$moments,
            crossNoise = mode
        )
        random <- methods::initialize(fixed, effects = "random")
        fits <- map(
            list(
                lava = list(fixed, one),
                oneRandom = list(random, one),
                twoRandom = list(random, two),
                twoFixed = list(fixed, two)
            ),
            function(at) {
                set.seed(1L)
                pecotmr:::.rgFit(at[[1L]], at[[2L]])[[1L]]$values
            }
        )
        for (name in c("oneRandom", "twoRandom", "twoFixed")) {
            expect_equal(
                unlist(fits[[name]][estimates]),
                unlist(fits$lava[estimates]),
                info = paste(mode, name)
            )
        }
        expect_identical(
            fits$twoRandom[inference],
            fits$oneRandom[inference],
            info = mode
        )
        expect_equal(
            fits$twoFixed$p_value,
            fits$lava$p_value,
            tolerance = 0.5,
            info = mode
        )
        expect_equal(
            c(fits$twoFixed$rg_lower, fits$twoFixed$rg_upper),
            c(fits$lava$rg_lower, fits$lava$rg_upper),
            tolerance = 0.05,
            info = mode
        )
        # The genetic variance's draws widen the random reading's interval.
        expect_gt(
            fits$oneRandom$rg_upper - fits$oneRandom$rg_lower,
            fits$lava$rg_upper - fits$lava$rg_lower
        )
    }
})

test_that("a group's moments given its others read the effects alike", {
    three <- .rgAsTwoLd(.rgTrivariateUnitData())
    for (effects in c("fixed", "random")) {
        param <- methods::initialize(
            .rgZeroSettings()$moments,
            effects = effects
        )
        set.seed(3L)
        one <- pecotmr:::.rgFit(param, c(three$one, list(conditional = TRUE)))
        set.seed(3L)
        mixed <- pecotmr:::.rgFit(
            param,
            c(three$two, list(conditional = TRUE))
        )
        for (col in c("rg", "covariance", "r2_1", "r2_2")) {
            expect_equal(
                map_dbl(mixed, c("values", col)),
                map_dbl(one, c("values", col)),
                info = paste(effects, col)
            )
        }
        # Random effects draw alike on one reference and on two; fixed
        # effects on two are within the draws' error of LAVA's Wishart.
        tolerance <- if (effects == "random") 1e-10 else 0.25
        for (col in c("rg_lower", "rg_upper", "p_value")) {
            expect_equal(
                map_dbl(mixed, c("values", col)),
                map_dbl(one, c("values", col)),
                tolerance = tolerance,
                info = paste(effects, col)
            )
        }
    }
})

test_that("fixed-effects means hold a group's cross-products", {
    # Three conditions on two bases of a 20-variant window, six components
    # each: the means' cross-products through C are H's, and the first
    # condition built keeps its observed direction.
    set.seed(9L)
    basis <- function(a) {
        list(
            lam = sort(stats::rexp(6L), decreasing = TRUE),
            vectors = qr.Q(qr(matrix(stats::rnorm(120L), 20L, 6L)))
        )
    }
    bases <- map(1:2, basis)
    shapes <- pecotmr:::.rgMixedMomentShapes(bases, "fixed")
    own <- map(1:3, function(k) list(delta = stats::rnorm(6L)))
    basisOf <- c(1L, 1L, 2L)
    H <- matrix(c(0.3, 0.1, 0.05, 0.1, 0.2, 0.02, 0.05, 0.02, 0.25), 3L)
    crossOf <- function(mus, k, l) {
        sum(mus[[k]] * (shapes[[basisOf[[k]], basisOf[[l]]]]$C %*% mus[[l]]))
    }
    for (order in list(1:3, 3:1)) {
        mus <- pecotmr:::.rgFixedMeans(H, own, shapes, basisOf, order)
        cells <- expand.grid(k = 1:3, l = 1:3)
        expect_equal(
            matrix(map2_dbl(cells$k, cells$l, crossOf, mus = mus), 3L),
            H
        )
        first <- order[[1L]]
        expect_equal(
            mus[[first]] / sqrt(H[[first, first]]),
            own[[first]]$delta / sqrt(sum(own[[first]]$delta^2))
        )
    }
})

# The first two traits of .rgTrivariateUnitData(), the second read against
# the LD of the first 200 samples alone: two references that differ, whose
# eigenvectors overlap densely.
.rgDenseTwoLd <- function() {
    ud <- pecotmr:::.rgPairUnitData(.rgTrivariateUnitData(), 1:2)
    data(eqtlRegionExample, envir = environment())
    G <- eqtlRegionExample$X
    af <- colMeans(G) / 2
    X <- scale(G[, pmin(af, 1 - af) > 0.05][, 1:300])
    other <- list(R = stats::cor(X[1:200, ]), nSamples = 200L)
    list_assign(ud, ld = NULL, lds = list(ud$ld, other), ldOf = 1:2)
}

# The log-likelihood's central difference in parameter j.
.rgJointSlope <- function(j, theta, data, h = 1e-6) {
    f <- pecotmr:::.rgJointObjective
    up <- f(replace(theta, j, theta[[j]] + h), data, FALSE)$value
    down <- f(replace(theta, j, theta[[j]] - h), data, FALSE)$value
    (up - down) / (2 * h)
}

test_that("the joint likelihood's score on two references is its slope", {
    # On references that differ the decomposition's blocks are dense, and a
    # block's quadratic term is r' V^-1 dV V^-1 r, which a' (V^-1 dV) r is
    # not where V^-1 dV is not symmetric.
    ud <- .rgDenseTwoLd()
    param <- .rgJointSettings()
    bases <- map(
        seq_along(ud$lds),
        pecotmr:::.rgBasisEigen,
        param = param,
        ud = ud
    )
    step1 <- pecotmr:::.rgMixedStep1(param, ud, bases)
    joint <- pecotmr:::.rgMixedData(
        bases,
        step1,
        ud$ldOf,
        param@procedures@varianceFloor,
        pecotmr:::.rgTieOf(param, ud, 2L)
    )
    bounds <- pecotmr:::.rgJointBounds(joint, param@procedures, NULL)
    data <- pecotmr:::.rgJointScoringData(joint, NULL, bounds)
    for (theta in list(c(0.2, 0.1), c(-0.3, 0.4))) {
        score <- pecotmr:::.rgJointObjective(theta, data, TRUE)$score
        slope <- map_dbl(seq_along(theta), .rgJointSlope, theta, data)
        expect_equal(score, slope, tolerance = 1e-5)
    }
})

test_that("a group on several references is the joint fit when they agree", {
    three <- .rgAsTwoLd(.rgTrivariateUnitData())
    cases <- expand.grid(
        mode = c("zero", "window"),
        given = c(FALSE, TRUE),
        stringsAsFactors = FALSE
    )
    for (k in seq_len(nrow(cases))) {
        joint <- methods::initialize(
            .rgJointSettings(),
            crossNoise = cases$mode[[k]]
        )
        given <- list(conditional = cases$given[[k]])
        one <- pecotmr:::.rgFit(joint, c(three$one, given))
        mixed <- pecotmr:::.rgFit(joint, c(three$two, given))
        for (col in c(
            "rg",
            "rg_lower",
            "rg_upper",
            "p_value",
            "noise_correlation"
        )) {
            expect_equal(
                map_dbl(mixed, c("values", col)),
                map_dbl(one, c("values", col)),
                tolerance = 1e-5,
                info = paste(col, cases$mode[[k]], cases$given[[k]])
            )
        }
    }
})

# Two cohorts of eqtlRegionExample's people (207 and 208, disjoint), each a
# study of context a on its own in-sample LD reference, sharing polygenic
# effects; `swap` codes cohort B's variants with their alleles the other way
# round (z negated, LD's sign flipped), as another population's panel may.
.rgTwoCohortFixture <- function(swap = FALSE, M = 150L) {
    data(eqtlRegionExample, envir = environment())
    G <- eqtlRegionExample$X
    set.seed(21L)
    perm <- sample.int(nrow(G))
    cohorts <- list(perm[1:207], perm[208:415])
    af <- map(cohorts, .rgCohortAlleleFrequency, G = G)
    usable <- reduce(map(af, .rgCommonEnough), `&`) &
        !.rgStrandAmbiguous(colnames(G))
    keep <- which(usable)[seq_len(M)]
    effects <- rnorm(M)
    sides <- map2(
        cohorts,
        c(FALSE, swap),
        .rgCohortSide,
        G = G[, keep],
        effects = effects
    )
    QtlSumStats(
        studyName = c("A", "B"),
        context = c("a", "a"),
        trait = c("g1", "g1"),
        entry = map(sides, "entry"),
        genome = "hg19",
        ldReference = map(sides, "ld"),
        qcOptions = list(qc = "ok")
    )
}

.rgCohortAlleleFrequency <- function(rows, G) {
    colMeans(G[rows, ]) / 2
}

.rgCommonEnough <- function(af) {
    pmin(af, 1 - af) > 0.05
}

# Variants whose two alleles pair on the other strand (A/T, C/G): across two
# references a swap of one cannot be told from a strand flip, and it is dropped.
.rgStrandAmbiguous <- function(ids) {
    parts <- stringr::str_split_fixed(ids, ":", 4)
    pairs <- paste0(parts[, 3], parts[, 4])
    pairs %in% c("AT", "TA", "CG", "GC")
}

# One cohort's entry and in-sample LdMatrix.
.rgCohortSide <- function(rows, flip, G, effects) {
    X <- scale(G[rows, ])
    g <- as.vector(X %*% effects)
    y <- g * sqrt(0.4 / var(g)) + rnorm(nrow(X), sd = sqrt(0.6))
    z <- pecotmr:::.corToZ(as.vector(cor(X, y)), nrow(X))
    parts <- stringr::str_split_fixed(colnames(G), ":", 4)
    if (flip) {
        parts <- parts[, c(1L, 2L, 4L, 3L)]
    }
    ids <- paste(parts[, 1], parts[, 2], parts[, 3], parts[, 4], sep = ":")
    sign <- if (flip) -1 else 1
    R <- crossprod(X) / (nrow(X) - 1)
    entry <- .rgEntry(sign * z, parts, ids, nrow(X))
    list(
        entry = entry,
        ld = LdMatrix(`dimnames<-`(R, list(ids, ids)), nSamples = nrow(X))
    )
}

test_that("two studies on their own references are matched by allele", {
    est <- suppressWarnings(estimateLocalRg(
        .rgTwoCohortFixture(),
        jointSpecification = "study",
        referenceOverlap = 1
    ))
    expect_identical(est$noise_source, "window")
    expect_true(is.finite(est$rg) && est$rg > 0)
    # Cohort B coded the other way round reads the same once matched.
    swapped <- suppressWarnings(estimateLocalRg(
        .rgTwoCohortFixture(swap = TRUE),
        jointSpecification = "study",
        referenceOverlap = 1
    ))
    expect_equal(swapped$rg, est$rg, tolerance = 1e-8)
    expect_equal(swapped$h2_2, est$h2_2, tolerance = 1e-8)
    # The moments on two references need no matrixsampling for their interval.
    moments <- suppressWarnings(estimateLocalRg(
        .rgTwoCohortFixture(),
        param = .rgZeroSettings()$moments,
        jointSpecification = "study",
        referenceOverlap = 1
    ))
    expect_true(is.finite(moments$rg_lower) && is.finite(moments$rg_upper))
})

test_that("a GWAS side reads its own reference where it names one", {
    qs <- .rgQtlFixture(contexts = "a")
    geno <- .rgFixtureGenotypes(150L)
    own <- GwasSumStats(
        studyName = "gw",
        entry = list(.rgGwasEntry(1L, qs[[1L]])),
        genome = "hg19",
        ldReference = .rgFixtureLd(geno, 400L)
    )
    est <- suppressWarnings(estimateLocalRg(qs, gwas = own))
    expect_identical(h2EstimateParam(est)@crossNoise, "window")
    expect_identical(est$noise_source, "window")
    expect_true(is.finite(est$rg))
    # A GWAS on the row's own reference stays on that one LD.
    shared <- suppressWarnings(estimateLocalRg(qs, gwas = .rgGwasFixture(qs)))
    expect_identical(h2EstimateParam(shared)@crossNoise, "window")
})

# Two studies of qtlDatasetExample's 165 donors (82 and 83), contexts a and b
# with polygenic effects shared across the studies (a with b at 0.6).
.rgTwoStudyFixture <- function() {
    data(qtlDatasetExample, envir = environment())
    panel <- pecotmr:::.qtlGenotypeSe(qtlDatasetExample)
    brain <- pecotmr:::.qtlContextSe(qtlDatasetExample, "brain")
    G <- scale(t(as.matrix(SummarizedExperiment::assay(panel))))
    G[!is.finite(G)] <- 0
    set.seed(4L)
    b <- matrix(rnorm(2L * ncol(G)), ncol = 2L) %*%
        chol(matrix(c(1, 0.6, 0.6, 1), 2L))
    e <- matrix(rnorm(2L * nrow(G)), ncol = 2L)
    y <- scale(G %*% b) * sqrt(0.5) + e * sqrt(0.5)
    contexts <- stats::setNames(nm = c("a", "b"))
    full <- map(
        c(s1 = "s1", s2 = "s2"),
        .rgStudyOf,
        panel = panel,
        phenotypes = map(contexts, .rgContextOf, y = y, brain = brain)
    )
    donors <- rownames(G)
    MultiStudyQtlDataset(
        qtlDatasets = list(
            s1 = full$s1[, donors[1:82], ],
            s2 = full$s2[, donors[83:165], ]
        )
    )
}

.rgStudyOf <- function(name, panel, phenotypes) {
    QtlDataset(studyName = name, genotypes = panel, phenotypes = phenotypes)
}

.rgContextOf <- function(cx, y, brain) {
    se <- brain[, rownames(y)]
    SummarizedExperiment::assay(se) <- matrix(
        y[, match(cx, c("a", "b"))],
        1L,
        dimnames = dimnames(se)
    )
    se
}

test_that("a MultiStudyQtlDataset pairs studies on their own genotypes", {
    skip_if_not_installed("PlinkMatrix")
    ms <- .rgTwoStudyFixture()
    # Within each study: its contexts on its genotypes, the noise fitted.
    within <- suppressWarnings(estimateLocalRg(ms, cisWindow = 1e6))
    expect_identical(within$study_1, c("s1", "s2"))
    expect_identical(within$study_1, within$study_2)
    expect_identical(h2EstimateParam(within)@crossNoise, "window")
    # Across studies: each context on the two studies' own genotypes.
    across <- suppressWarnings(estimateLocalRg(
        ms,
        cisWindow = 1e6,
        jointSpecification = "study"
    ))
    expect_identical(across$context_1, c("a", "b"))
    expect_identical(
        c(across$study_1[[1L]], across$study_2[[1L]]),
        c("s1", "s2")
    )
    expect_identical(across$noise_source, c("window", "window"))
    expect_true(all(across$rg > 0.5))
})

test_that("tied noise ties a two-study group's pairs that share donors", {
    skip_if_not_installed("PlinkMatrix")
    ms <- .rgTwoStudyFixture()
    tied <- function(param) methods::initialize(param, crossNoise = "tied")
    fits <- list(
        pairs = list(tied(H2EstimateParam()), FALSE),
        joint = list(tied(.rgJointSettings()), FALSE),
        given = list(tied(.rgSettings()$moments), TRUE)
    )
    for (name in names(fits)) {
        est <- suppressWarnings(estimateLocalRg(
            ms,
            cisWindow = 1e6,
            param = fits[[name]][[1L]],
            jointSpecification = list(c("study", "context")),
            conditional = fits[[name]][[2L]]
        ))
        # A study's two contexts share its donors and tie their noise; two
        # studies' conditions share none.
        same <- est$study_1 == est$study_2
        expect_identical(sum(same), 2L, info = name)
        expect_identical(est$noise_source[same], c("tied", "tied"), info = name)
        expect_true(all(est$noise_source[!same] == "zero"), info = name)
        expect_true(all(est$noise_correlation[!same] == 0), info = name)
        expect_true(all(is.finite(est$noise_correlation[same])), info = name)
    }
    # A summary row carries no phenotypes to tie.
    expect_error(
        pecotmr:::.rgAssertTiedDonors(tied(H2EstimateParam()), list(), "u"),
        "summary row"
    )
})

# ---- the frequency exponent -------------------------------------------------

test_that("frequency weights carry to two references", {
    # With alpha = 0 and the same LD and frequencies given as two references,
    # each two-reference form is its one-reference fit, as at alpha = -1.
    pair <- .rgAsTwoLd(pecotmr:::.rgPairUnitData(.rgTrivariateUnitData(), 1:2))
    M <- length(pair$one$z[[1L]])
    p <- seq(0.05, 0.5, length.out = M)
    het <- 2 * p * (1 - p)
    one <- list_assign(pair$one, ld = list_assign(pair$one$ld, het = het))
    two <- list_assign(
        pair$two,
        lds = map(pair$two$lds, list_assign, het = het)
    )
    cols <- c("rg", "rg_lower", "rg_upper", "covariance", "h2_1", "h2_2")
    for (name in names(.rgZeroSettings())) {
        param <- methods::initialize(.rgZeroSettings()[[name]], alpha = 0)
        a <- pecotmr:::.rgFit(param, one)[[1L]]$values
        b <- pecotmr:::.rgFit(param, two)[[1L]]$values
        compared <- if (name == "moments") {
            c("rg", "covariance", "h2_1", "h2_2")
        } else {
            c(cols, "p_value")
        }
        expect_equal(
            unlist(b[compared]),
            unlist(a[compared]),
            tolerance = 1e-4,
            info = name
        )
    }
})

# An LdMatrix of the fixture's genotypes, recording each variant's A1
# frequency.
.rgFixtureLdAf <- function(geno) {
    parts <- stringr::str_split_fixed(geno$ids, ":", 4)
    R <- crossprod(geno$X) / (geno$n - 1)
    LdMatrix(
        `dimnames<-`(R, list(geno$ids, geno$ids)),
        variants = data.frame(
            chrom = parts[, 1],
            pos = as.integer(parts[, 2]),
            A1 = parts[, 4],
            A2 = parts[, 3],
            AF = geno$af
        ),
        nSamples = geno$n
    )
}

test_that("estimateLocalRg weighs each reference's variants by alpha", {
    qs <- .rgQtlFixture(contexts = c("a", "b"))
    alpha0 <- methods::initialize(H2EstimateParam(), alpha = 0)
    expect_error(
        estimateLocalRg(qs, param = alpha0),
        "records none"
    )
    geno <- .rgFixtureGenotypes(150L)
    withAf <- QtlSumStats(
        studyName = c("s1", "s1"),
        context = c("a", "b"),
        trait = c("g1", "g1"),
        entry = list(qs[[1L]], qs[[2L]]),
        genome = "hg19",
        ldReference = .rgFixtureLdAf(geno),
        qcOptions = list(qc = "ok")
    )
    est <- suppressWarnings(estimateLocalRg(withAf, param = alpha0))
    expect_identical(h2EstimateParam(est)@alpha, 0)
    expect_true(is.finite(est$rg))
})

test_that("a MultiStudyQtlDataset weighs each study's variants by alpha", {
    skip_if_not_installed("PlinkMatrix")
    est <- suppressWarnings(estimateLocalRg(
        .rgTwoStudyFixture(),
        cisWindow = 1e6,
        param = methods::initialize(H2EstimateParam(), alpha = 0),
        jointSpecification = "study"
    ))
    expect_true(all(is.finite(est$rg)))
})

# ---- edge cases and fallbacks ------------------------------------------------

# Whether no pair of `fits` has a covariance.
.rgNoCovariance <- function(fits) {
    all(is.na(map_dbl(fits, c("values", "covariance"))))
}

# `param` asking for more directions than any window here has.
.rgTooFew <- function(param) {
    methods::initialize(
        param,
        truncation = methods::initialize(param@truncation, minimum = 1000)
    )
}

test_that("a window too small for a setting leaves its pairs no estimate", {
    ud <- .rgUnitData(.rgBivariateFixture(2L, c(0.3, 0.3), 0.6, 0.3))
    s <- .rgSettings()
    two <- list_assign(
        ud,
        z = map(ud$z, `[`, 1:2),
        nVariant = map(ud$nVariant, `[`, 1:2),
        ld = list(R = ud$ld$R[1:2, 1:2], nSamples = ud$ld$nSamples)
    )
    for (name in c("variantRegression", "gls")) {
        expect_true(
            .rgNoCovariance(pecotmr:::.rgFit(s[[name]], two)),
            info = name
        )
    }
    eigen <- list(
        likelihood = s$likelihood,
        regression = s$regression,
        fixed = s$moments,
        random = methods::initialize(s$moments, effects = "random")
    )
    for (name in names(eigen)) {
        fits <- pecotmr:::.rgFit(.rgTooFew(eigen[[name]]), ud)
        expect_true(.rgNoCovariance(fits), info = name)
    }
    ud3 <- .rgTrivariateUnitData()
    expect_true(.rgNoCovariance(
        pecotmr:::.rgFit(.rgTooFew(.rgJointSettings()), ud3)
    ))
    # On two references, each side's own basis.
    pair <- .rgAsTwoLd(pecotmr:::.rgPairUnitData(ud3, 1:2))
    small <- .rgAsTwoLd(pecotmr:::.rgPairUnitData(ud3, 1:2), M = 2L)
    zero <- .rgZeroSettings()
    for (name in c("likelihood", "regression", "moments")) {
        fits <- pecotmr:::.rgFit(.rgTooFew(zero[[name]]), pair$two)
        expect_true(.rgNoCovariance(fits), info = name)
    }
    for (name in c("variantRegression", "gls")) {
        fits <- pecotmr:::.rgFit(zero[[name]], small$two)
        expect_true(.rgNoCovariance(fits), info = name)
    }
    three <- .rgAsTwoLd(ud3)
    mixed <- methods::initialize(.rgJointSettings(), crossNoise = "window")
    expect_true(.rgNoCovariance(
        pecotmr:::.rgFit(.rgTooFew(mixed), three$two)
    ))
    given <- list_assign(three$two, conditional = TRUE)
    expect_true(.rgNoCovariance(
        pecotmr:::.rgFit(.rgTooFew(zero$moments), given)
    ))
    # Fewer samples than components (here 42 hold the trace) leave the
    # moments no noise to read; LAVA caps the components at 0.75 of the
    # samples, so this keeps them all.
    uncapped <- methods::initialize(
        zero$moments,
        truncation = H2TruncationParam(
            fraction = 0.9999,
            of = "eigenvalues",
            minimum = 2
        )
    )
    scarce <- list_assign(
        given,
        nVariant = map(given$nVariant, `*`, 0.05)
    )
    expect_true(.rgNoCovariance(pecotmr:::.rgFit(uncapped, scarce)))
})

test_that("two references fit a group's pairs one at a time where asked", {
    three <- .rgAsTwoLd(.rgTrivariateUnitData())
    zero <- .rgZeroSettings()
    pairs <- pecotmr:::.rgGroupPairs(3L)
    for (name in c("likelihood", "moments")) {
        two <- map_dbl(
            pecotmr:::.rgFit(zero[[name]], three$two),
            c("values", "rg")
        )
        # Each pair as its own unit on the one reference.
        one <- map2_dbl(
            pairs$i,
            pairs$j,
            function(i, j) {
                unit <- pecotmr:::.rgPairUnitData(three$one, c(i, j))
                pecotmr:::.rgFit(zero[[name]], unit)[[1L]]$values$rg
            }
        )
        expect_length(two, 3L)
        expect_equal(two, one, tolerance = 1e-5, info = name)
    }
    # Each condition's noise read without the reference's overlap.
    s <- .rgSettings()
    plainOne <- s$moments
    plainOne@referenceNoise <- "none"
    ud <- .rgUnitData(.rgBivariateFixture(2L, c(0.3, 0.3), 0.6, 0.3))
    expect_true(is.finite(pecotmr:::.rgFit(plainOne, ud)[[1L]]$values$rg))
    plainTwo <- zero$moments
    plainTwo@referenceNoise <- "none"
    pair <- .rgAsTwoLd(pecotmr:::.rgPairUnitData(.rgTrivariateUnitData(), 1:2))
    expect_true(is.finite(
        pecotmr:::.rgFit(plainTwo, pair$two)[[1L]]$values$rg
    ))
})

test_that("a step-2 fit with no null, or none beating it, has no estimate", {
    fx <- .rgBivariateFixture(2L, c(0.3, 0.3), 0.6, 0.3)
    param <- .rgSettings()$likelihood
    param@procedures@mustBeatNull <- TRUE
    # The supplied noise correlation starts the window's fit too.
    started <- pecotmr:::.rgEigenLikelihood(
        param,
        .rgUnitData(fx, .rgRho(0.25))
    )[[1L]]$values
    expect_identical(started$noise_source, "window")
    expect_true(is.finite(started$covariance))
    for (null in list(NULL, list(value = 1e300, par = 0))) {
        local_mocked_bindings(.rgCrossNullFit = function(...) null)
        got <- pecotmr:::.rgEigenLikelihood(param, .rgUnitData(fx))
        expect_true(.rgNoCovariance(got))
    }
})

test_that("a step-1 fit at a held intercept that fails has no estimate", {
    fx <- .rgBivariateFixture(2L, c(0.3, 0.3), 0.6, 0.3)
    held <- list_assign(.rgUnitData(fx), heldIntercept = c(1, 1))
    local_mocked_bindings(.h2ScoringBest = function(...) NULL)
    got <- pecotmr:::.rgEigenLikelihood(.rgSettings()$likelihood, held)
    expect_true(is.na(got[[1L]]$values$h2_1))
    expect_true(.rgNoCovariance(got))
})

test_that("a joint fit that fails falls back to the supplied noise", {
    ud <- .rgTrivariateUnitData()
    param <- .rgJointSettings()
    real <- pecotmr:::.rgJointBest
    local_mocked_bindings(
        .rgJointBest = function(param, ud, d, members, joint, mode, held) {
            if (identical(mode, "window")) {
                return(NULL)
            }
            real(param, ud, d, members, joint, mode, held)
        }
    )
    rho <- matrix(0.2, 3L, 3L) + diag(0.8, 3L)
    fallback <- map(
        pecotmr:::.rgFit(param, list_assign(ud, noiseCorrelation = rho)),
        "values"
    )
    expect_identical(map_chr(fallback, "noise_source"), rep("fallback", 3L))
    expect_true(all(is.finite(map_dbl(fallback, "rg"))))
    expect_true(.rgNoCovariance(pecotmr:::.rgFit(param, ud)))
})

test_that("a joint fit with no null, or none beating it, is no fit", {
    ud <- .rgTrivariateUnitData()
    param <- .rgJointSettings()
    param@procedures@mustBeatNull <- TRUE
    for (none in c(TRUE, FALSE)) {
        local_mocked_bindings(
            .rgJointCoordinateNull = function(
                theta,
                at,
                joint,
                proc,
                held,
                values = 0
            ) {
                if (none) {
                    return(NULL)
                }
                list(value = 1e300, par = replace(theta, at, values))
            }
        )
        expect_true(.rgNoCovariance(pecotmr:::.rgFit(param, ud)), info = none)
    }
})

test_that("a joint fit starts each pair at 0 where its own fit fails", {
    local_mocked_bindings(.rgLikelihoodCrossFit = function(...) NULL)
    fits <- pecotmr:::.rgFit(.rgJointSettings(), .rgTrivariateUnitData())
    expect_true(all(is.finite(map_dbl(fits, c("values", "rg")))))
})

test_that("a regression or GLS cross fit that fails falls back or has none", {
    fx <- .rgBivariateFixture(2L, c(0.3, 0.3), 0.6, 0.3)
    s <- .rgSettings()
    irwls <- pecotmr:::.rgCrossIrwls
    gls <- pecotmr:::.rgGlsCross
    local_mocked_bindings(
        .rgCrossIrwls = function(param, rows, cross, steps, held) {
            if (is.na(held)) {
                stop("singular")
            }
            irwls(param, rows, cross, steps, held)
        },
        .rgGlsCross = function(param, g, steps, held) {
            if (is.na(held)) {
                stop("singular")
            }
            gls(param, g, steps, held)
        }
    )
    for (name in c("regression", "variantRegression", "gls")) {
        fallback <- pecotmr:::.rgFit(
            s[[name]],
            .rgUnitData(fx, .rgRho(0.25))
        )[[1L]]$values
        expect_identical(fallback$noise_source, "fallback", info = name)
        expect_true(is.finite(fallback$covariance), info = name)
        none <- pecotmr:::.rgFit(s[[name]], .rgUnitData(fx))
        expect_true(.rgNoCovariance(none), info = name)
    }
    pair <- .rgAsTwoLd(pecotmr:::.rgPairUnitData(.rgTrivariateUnitData(), 1:2))
    expect_true(.rgNoCovariance(pecotmr:::.rgFit(s$gls, pair$two)))
})

test_that("the regressions hold an intercept and give SEs as asked", {
    ud <- .rgUnitData(.rgBivariateFixture(2L, c(0.3, 0.3), 0.6, 0.3))
    s <- .rgSettings()
    held <- pecotmr:::.rgFit(
        s$regression,
        list_assign(ud, heldIntercept = c(NA, 1))
    )[[1L]]$values
    expect_true(is.finite(held$covariance))
    for (name in c("regression", "gls")) {
        param <- s[[name]]
        param@inference <- "lrt"
        got <- pecotmr:::.rgFit(param, ud)[[1L]]$values
        expect_true(is.finite(got$rg) && is.na(got$rg_se), info = name)
    }
    pair <- .rgAsTwoLd(pecotmr:::.rgPairUnitData(.rgTrivariateUnitData(), 1:2))
    zero <- .rgZeroSettings()$regression
    zero@inference <- "lrt"
    expect_true(is.na(pecotmr:::.rgFit(zero, pair$two)[[1L]]$values$rg_se))
    noBlocks <- s$variantRegression
    noBlocks@procedures@jackknifeBlocks <- NULL
    expect_error(pecotmr:::.rgFit(noBlocks, ud), "jackknifeBlocks")
    # With no rg, the GLS gives the covariance's error alone.
    local_mocked_bindings(.rgCorrelation = function(...) NA_real_)
    alone <- pecotmr:::.rgFit(s$gls, ud)[[1L]]$values
    expect_true(is.finite(alone$covariance_se) && is.na(alone$rg_se))
})

test_that("the moments' draws and tails take their edge cases", {
    set.seed(1L)
    few <- pecotmr:::.rgCentralGram(4L, 1L, diag(2L))
    expect_identical(dim(few), c(2L, 2L, 4L))
    expect_true(all(apply(few, 3L, det) < 1e-12))
    expect_identical(
        pecotmr:::.rgCentralGram(3L, 0L, diag(2L)),
        array(0, c(2L, 2L, 3L))
    )
    H <- matrix(c(1, 0.4, 0.4, 1), 2L)
    expect_equal(pecotmr:::.rgMixedNullH(H, 1L, 2L, integer(0)), diag(2L))
    expect_identical(pecotmr:::.rgUnitVector(c(0, 0)), c(0, 0))
    expect_identical(
        pecotmr:::.rgJoinedColumns(NULL, 1:3, character(0), "|"),
        rep("all", 3L)
    )
    expect_identical(
        pecotmr:::.rgStatisticCorrelation(list(rep(9, 5), rep(9, 5))),
        0
    )
    # Non-positive heritabilities give the fixed effects no p or interval.
    flat <- list(effects = "fixed", H = diag(c(0, 1)))
    expect_true(is.na(pecotmr:::.rgFixedPairP(flat)))
    inference <- pecotmr:::.rgGroupPairInference(flat)
    expect_true(is.na(inference$p) && all(is.na(inference$interval)))
})

test_that("draws with no covariance matrix give no p-value or interval", {
    local_mocked_bindings(
        .rgMixedStacked = function(...) -diag(2L),
        .rgTwoLdJoint = function(...) -diag(2L),
        .rgSideComponentVariance = function(...) 1
    )
    H <- matrix(c(1, 0.4, 0.4, 1), 2L)
    expect_null(pecotmr:::.rgMixedDrawSources(list(), H))
    expect_true(is.na(
        pecotmr:::.rgMixedDrawP(list(), H, 1L, 2L, integer(0), 0.1, 10L)
    ))
    drawn <- pecotmr:::.rgMixedDrawn(list(H = H), 1L, 2L, integer(0), 0.1)
    expect_true(is.na(drawn$p) && all(is.na(drawn$rg)))
    cross <- list(
        sigma12 = 0.1,
        noiseShape = 1,
        covariance = 0,
        genetic = 0
    )
    own <- list(list(lam = 1, K = 1L), list(lam = 1, K = 1L))
    expect_true(is.na(pecotmr:::.rgTwoLdDrawP(cross, own, list(1, 1))))
    expect_true(all(is.na(pecotmr:::.rgTwoLdMomentsInterval(cross, own))))
})

test_that("small p-values are read off more draws", {
    own <- list(list(lam = 1, K = 1L), list(lam = 1, K = 1L))
    zeroNoise <- list(sigma12 = 0, C = matrix(0.5), statistic = 3)
    local_mocked_bindings(.rgSideComponentVariance = function(...) 1)
    for (deep in c(5e-7, 5e-6)) {
        local_mocked_bindings(
            .rgBilinearTail = function(rho, observed, nDraws, ...) {
                if (nDraws < 1e5) {
                    1e-5
                } else if (nDraws < 1e6) {
                    deep
                } else {
                    2e-7
                }
            }
        )
        want <- if (deep < 1e-6) 2e-7 else deep
        expect_identical(pecotmr:::.rgTwoLdMomentsP(zeroNoise, own), want)
    }
    local_mocked_bindings(
        .rgTwoLdJoint = function(...) diag(2L),
        .rgTwoLdDrawShare = function(nDraws, ...) {
            if (nDraws < 1e5) 1e-4 else 3e-5
        }
    )
    shared <- list(sigma12 = 0.1, noiseShape = 1)
    expect_identical(pecotmr:::.rgTwoLdDrawP(shared, own, list(1, 1)), 3e-5)
    # A fixed-effects pair with shared noise draws its components together;
    # without, it reads its tails.
    pair <- list(
        H = matrix(c(1, 0.3, 0.3, 1), 2L),
        sigma = matrix(c(1, 0.2, 0.2, 1), 2L)
    )
    local_mocked_bindings(
        .rgMixedDrawP = function(group, H, i, j, rest, observed, nDraws) {
            if (nDraws < 1e5) 1e-4 else 3e-5
        },
        .rgFixedTails = function(...) NULL,
        .rgFixedTail = function(tails, observed, nDraws) {
            if (nDraws < 1e5) 1e-5 else 5e-6
        }
    )
    expect_identical(pecotmr:::.rgFixedPairP(pair), 3e-5)
    pair$sigma <- diag(2L)
    expect_identical(pecotmr:::.rgFixedPairP(pair), 5e-6)
})

test_that("the summary's crossings and mash's correlation read their edges", {
    excess <- rep(-1, 101L)
    expect_identical(pecotmr:::.rgCurveCrossing(51:101, excess, 1), 1)
    crossed <- replace(excess, 51:52, 1)
    expect_identical(
        pecotmr:::.rgCurveCrossing(52:101, crossed, 1),
        pecotmr:::.RG_CURVE_GRID[[52L]]
    )
    implied <- matrix(1, 2L, 2L, dimnames = list(c("x", "y"), c("x", "y")))
    identity <- list(
        study_1 = "s",
        context_1 = "a",
        study_2 = "s",
        context_2 = "b"
    )
    expect_true(is.na(pecotmr:::.rgMashCorrelation(implied, identity)))
})

test_that("a reduced fit that does not converge has no estimate", {
    local_mocked_bindings(.h2ScoringBest = function(...) NULL)
    got <- pecotmr:::.rgReducedResult(
        .rgSettings()$likelihood,
        list(bound = 0.5),
        list(nSnps = 10L)
    )
    expect_true(is.na(got$values$covariance))
    expect_identical(got$values$n_snps, 10L)
})

test_that("a unit whose conditions share too few variants has no pairs", {
    qs <- .rgQtlFixture(contexts = c("a", "b"))
    local_mocked_bindings(
        .rgCommonHeldBy = function(common, refs, label) {
            list(ids = character(0))
        }
    )
    res <- suppressWarnings(estimateLocalRg(qs))
    expect_length(res, 1L)
    expect_true(is.na(res$rg))
    expect_identical(res$n_snps, 0L)
})

test_that("a GWAS side outside a window is empty, on two references refused", {
    qs <- .rgQtlFixture(contexts = c("a", "b"))
    gw <- .rgGwasFixture(qs)
    expect_identical(
        nrow(pecotmr:::.rgGwasMember(gw, "gw", NULL)$df),
        0L
    )
    expect_null(pecotmr:::.rgGwasRef(gw, "gw", NULL))
    entry <- gw[[1L]]
    half <- seq_len(length(entry) %/% 2L)
    geno <- .rgFixtureGenotypes(150L)
    split <- GwasSumStats(
        studyName = c("gw", "gw"),
        entry = list(entry[half], entry[-half]),
        genome = "hg19",
        ldReference = list(
            .rgFixtureLd(geno),
            .rgFixtureLd(geno, nSamples = 100L)
        ),
        blockId = c("b1", "b2")
    )
    expect_error(
        pecotmr:::.rgGwasRef(split, "gw", range(entry)),
        "different LD references"
    )
})

test_that("a QtlDataset pair on fewer than three shared donors has none", {
    skip_if_not_installed("PlinkMatrix")
    data(qtlDatasetExample, envir = environment())
    panel <- pecotmr:::.qtlGenotypeSe(qtlDatasetExample)
    brain <- pecotmr:::.qtlContextSe(qtlDatasetExample, "brain")
    donors <- colnames(panel)
    set.seed(2L)
    y <- matrix(
        rnorm(3L * length(donors)),
        ncol = 3L,
        dimnames = list(donors, NULL)
    )
    drops <- list(a = 3:length(donors), b = 1:2)
    qd <- QtlDataset(
        studyName = "s1",
        genotypes = panel,
        phenotypes = map(
            c(a = "a", b = "b"),
            .rgDatasetContext,
            y = y,
            drops = drops,
            brain = brain
        )
    )
    res <- suppressWarnings(estimateLocalRg(qd, cisWindow = 1e6))
    expect_true(all(is.na(res$rg)))
})

test_that("a MultiStudyQtlDataset pairs its conditions with a GWAS", {
    skip_if_not_installed("PlinkMatrix")
    ms <- .rgTwoStudyFixture()
    gw <- .rgDatasetGwas(as.list(ms)$s1)
    each <- suppressWarnings(estimateLocalRg(ms, gwas = gw, cisWindow = 1e6))
    expect_length(each, 4L)
    expect_true(all(each$study_2 == "gw"))
    # The GWAS carries s1's signal; on s2's donors its heritability sits at
    # its bound of 0, which holds the covariance at 0.
    own <- each$study_1 == "s1"
    expect_true(all(is.finite(each$rg[own])))
    expect_true(all(each$covariance[!own] == 0))
    given <- suppressWarnings(estimateLocalRg(
        ms,
        gwas = gw,
        cisWindow = 1e6,
        param = .rgSettings()$moments,
        conditional = TRUE
    ))
    expect_true(all(given$study_2 == "gw"))
    expect_error(
        estimateLocalRg(ms, cisWindow = 1e6, jointSpecification = "trait"),
        "no two selected conditions share a group"
    )
    expect_error(estimateLocalRg(ms, cisWindow = 1e6, traitId = "none"))
})

# A summary-statistics study "s3" of context a's trait, on study s1's donors
# with their LD as its reference, beside a MultiStudyQtlDataset's studies;
# its trait anchored `shift` bases from the studies'.
.rgSummaryStudy <- function(ms, shift = 0L) {
    qd <- as.list(ms)$s1
    panel <- pecotmr:::.qtlGenotypeSe(qd)
    rr <- SummarizedExperiment::rowRanges(panel)
    a <- pecotmr:::.qtlContextSe(qd, "a")
    y <- SummarizedExperiment::assay(a)[1L, ]
    G <- t(as.matrix(SummarizedExperiment::assay(panel)))
    donors <- intersect(names(y)[!is.na(y)], rownames(G))
    X <- scale(G[donors, ])
    X[!is.finite(X)] <- 0
    r <- suppressWarnings(as.vector(stats::cor(X, y[donors])))
    r <- replace(r, !is.finite(r), 0)
    ids <- paste(
        GenomicRanges::seqnames(rr),
        GenomicRanges::start(rr),
        rr$A2,
        rr$A1,
        sep = ":"
    )
    gr <- GenomicRanges::GRanges(
        GenomicRanges::seqnames(rr),
        IRanges::IRanges(GenomicRanges::start(rr), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = ids,
        A1 = rr$A1,
        A2 = rr$A2,
        Z = pecotmr:::.corToZ(r, length(donors)),
        N = length(donors)
    )
    R <- stats::cov2cor(
        crossprod(X) + diag(as.numeric(colSums(X^2) == 0), ncol(X))
    )
    QtlSumStats(
        studyName = "s3",
        context = "a",
        trait = rownames(a),
        entry = list(gr),
        genome = "hg19",
        ldReference = LdMatrix(
            `dimnames<-`(R, list(ids, ids)),
            nSamples = length(donors)
        ),
        traitPos = IRanges::shift(SummarizedExperiment::rowRanges(a), shift),
        qcOptions = list(qc = "ok")
    )
}

test_that("a MultiStudyQtlDataset pairs a summary study with its own", {
    skip_if_not_installed("PlinkMatrix")
    ms <- .rgTwoStudyFixture()
    withSummary <- MultiStudyQtlDataset(
        qtlDatasets = as.list(ms),
        sumStats = .rgSummaryStudy(ms)
    )
    across <- suppressWarnings(estimateLocalRg(
        withSummary,
        cisWindow = 1e6,
        jointSpecification = "study"
    ))
    summary <- across$study_1 == "s3" | across$study_2 == "s3"
    expect_identical(sum(summary), 2L)
    expect_true(all(across$context_1[summary] == "a"))
    expect_true(all(is.finite(across$rg[summary])))
    # A selection the summary study holds none of leaves it out.
    onlyB <- suppressWarnings(estimateLocalRg(
        withSummary,
        cisWindow = 1e6,
        contexts = "b",
        jointSpecification = "study"
    ))
    expect_false(any(c(onlyB$study_1, onlyB$study_2) == "s3"))
    expect_identical(onlyB$context_1, "b")
})

test_that("a summary study in a window of its own pairs with a GWAS alone", {
    skip_if_not_installed("PlinkMatrix")
    ms <- .rgTwoStudyFixture()
    # Its trait anchored 50 kb off the studies': a window of its own, which
    # neither genotyped study has a condition in.
    shifted <- MultiStudyQtlDataset(
        qtlDatasets = as.list(ms),
        sumStats = .rgSummaryStudy(ms, shift = 50000L)
    )
    gw <- .rgDatasetGwas(as.list(ms)$s1)
    each <- suppressWarnings(estimateLocalRg(
        shifted,
        gwas = gw,
        cisWindow = 1e6
    ))
    expect_length(each, 5L)
    summary <- each$study_1 == "s3"
    expect_identical(sum(summary), 1L)
    expect_identical(each$study_2[summary], "gw")
    expect_gt(each$n_snps[summary], 2L)
    # A unit with no window keys as none.
    expect_identical(pecotmr:::.rgMsWindowKey(list(window = NULL)), "none")
})

# Two studies of contexts a and b on the example's genotypes: s1's two
# contexts measured on disjoint donors (two and the rest), s2's on all.
.rgSparseStudies <- function() {
    data(qtlDatasetExample, envir = environment())
    panel <- pecotmr:::.qtlGenotypeSe(qtlDatasetExample)
    brain <- pecotmr:::.qtlContextSe(qtlDatasetExample, "brain")
    donors <- colnames(panel)
    set.seed(6L)
    y <- matrix(
        rnorm(3L * length(donors)),
        ncol = 3L,
        dimnames = list(donors, NULL)
    )
    drops <- list(
        s1 = list(a = 3:length(donors), b = 1:2),
        s2 = list(a = integer(0), b = integer(0))
    )
    MultiStudyQtlDataset(qtlDatasets = imap(drops, .rgSparseStudy, y = y))
}

# @noRd
.rgSparseStudy <- function(drops, name, y) {
    data(qtlDatasetExample, envir = environment())
    QtlDataset(
        studyName = name,
        genotypes = pecotmr:::.qtlGenotypeSe(qtlDatasetExample),
        phenotypes = map(
            c(a = "a", b = "b"),
            .rgDatasetContext,
            y = y,
            drops = drops,
            brain = pecotmr:::.qtlContextSe(qtlDatasetExample, "brain")
        )
    )
}

test_that("a study's conditions sharing too few donors give no pair", {
    skip_if_not_installed("PlinkMatrix")
    est <- suppressWarnings(
        estimateLocalRg(.rgSparseStudies(), cisWindow = 1e6)
    )
    expect_identical(est$study_1, c("s1", "s2"))
    expect_identical(est$n_snps[[1L]], 0L)
    expect_true(is.na(est$rg[[1L]]))
    expect_gt(est$n_snps[[2L]], 2L)
})

test_that("a QtlDataset group within a scope pairs only its conditions", {
    skip_if_not_installed("PlinkMatrix")
    est <- suppressWarnings(estimateLocalRg(
        .rgDatasetFixture(),
        cisWindow = 1e6,
        jointSpecification = list(list(
            axes = "context",
            scope = list(context = c("a", "c"))
        ))
    ))
    expect_length(est, 1L)
    expect_identical(c(est$context_1, est$context_2), c("a", "c"))
})

test_that("random-effects moments need more samples than components", {
    ud <- .rgUnitData(.rgBivariateFixture(2L, c(0.3, 0.3), 0.6, 0.3))
    # Every component kept (LAVA would cap them at 0.75 of the samples),
    # and fewer samples than they number.
    random <- methods::initialize(
        .rgSettings()$moments,
        effects = "random",
        truncation = H2TruncationParam(
            fraction = 0.9999,
            of = "eigenvalues",
            minimum = 2
        )
    )
    scarce <- list_assign(ud, nVariant = map(ud$nVariant, `*`, 0.05))
    expect_true(.rgNoCovariance(pecotmr:::.rgFit(random, scarce)))
})

test_that("two references' step 1 without a positive heritability gives none", {
    pair <- .rgAsTwoLd(pecotmr:::.rgPairUnitData(.rgTrivariateUnitData(), 1:2))
    zero <- .rgZeroSettings()$likelihood
    for (h2 in c(NA_real_, 0)) {
        local_mocked_bindings(
            .rgLikelihoodStep1 = function(member, held, param) {
                list(h2 = h2, intercept = 1)
            }
        )
        got <- pecotmr:::.rgFit(zero, pair$two)[[1L]]$values
        expect_true(is.na(got$rg), info = h2)
        expect_identical(
            got$covariance,
            if (is.na(h2)) NA_real_ else 0,
            info = h2
        )
    }
})

test_that("a joint fit on several references that fails leaves no pair", {
    three <- .rgAsTwoLd(.rgTrivariateUnitData())
    mixed <- methods::initialize(.rgJointSettings(), crossNoise = "window")
    local_mocked_bindings(.rgJointFit = function(...) NULL)
    expect_true(.rgNoCovariance(pecotmr:::.rgFit(mixed, three$two)))
})
