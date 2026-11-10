# Shared fixtures for the TWAS weight-learner wrappers. Used by both
# test_regularizedRegressionWrappers.R (lasso/enet/scad/... + RSS solvers) and
# test_fineMappingWrappers.R (the SuSiE/mvSuSiE/fSuSiE weight extractors), so
# they live in a helper that testthat auto-loads for every test file.

# Individual-level (X, y) fixture: n samples x p variants with a sparse signal.
.rrwXy <- function(n = 50, p = 6, seed = 1) {
    set.seed(seed)
    X <- matrix(
        rnorm(n * p),
        n,
        p,
        dimnames = list(NULL, sprintf("chr1:%d:A:G", 100L * seq_len(p)))
    )
    b <- rnorm(p)
    b[-(1:2)] <- 0
    y <- as.numeric(X %*% b + rnorm(n))
    list(X = X, y = y, n = n, p = p)
}

# Single-context summary statistics + LD for the RSS weights contract:
# `stat` names its variants, and `ldReference` holds the same LD as an
# LdMatrix for the weight functions to match them to.
.rrwStatLd <- function(n = 50, p = 6, seed = 1) {
    d <- .rrwXy(n, p, seed)
    bhat <- map_dbl(
        seq_len(p),
        function(j) summary(lm(d$y ~ d$X[, j]))$coefficients[2, 1]
    )
    sehat <- map_dbl(
        seq_len(p),
        function(j) summary(lm(d$y ~ d$X[, j]))$coefficients[2, 2]
    )
    zhat <- bhat / sehat
    LD <- cor(d$X)
    stat <- list(
        b = bhat,
        seb = sehat,
        z = zhat,
        cor = bhat,
        n = rep(n, p),
        var_y = 1,
        variantNames = colnames(LD)
    )
    list(
        stat = stat,
        LD = LD,
        ldReference = LdMatrix(LD, nSamples = n),
        p = p,
        n = n
    )
}

# Multi-context (variants x conditions) fixture for mr.mash / mvSuSiE wrappers.
.rrwMulti <- function(n = 60, p = 6, K = 3, seed = 2) {
    set.seed(seed)
    X <- matrix(
        rnorm(n * p),
        n,
        p,
        dimnames = list(NULL, sprintf("chr1:%d:A:G", 100L * seq_len(p)))
    )
    B <- matrix(0, p, K)
    B[1, ] <- rnorm(K, sd = 2)
    B[2, ] <- rnorm(K, sd = 2)
    Y <- X %*% B + matrix(rnorm(n * K), n, K)
    colnames(Y) <- paste0("ctx", seq_len(K))
    marginalZ <- map(seq_len(K), function(k) {
        map_dbl(seq_len(p), function(j) {
            summary(lm(Y[, k] ~ X[, j]))$coefficients[2, 3]
        })
    })
    Z <- exec(cbind, !!!marginalZ)
    colnames(Z) <- colnames(Y)
    list(
        X = X,
        Y = Y,
        LD = cor(X),
        ldReference = LdMatrix(cor(X)),
        stat = list(z = Z, n = n, variantNames = colnames(X)),
        p = p,
        K = K
    )
}

# Variant ids chr1:100:A:G, chr1:200:A:G, ... for p variants, and an LdMatrix
# over `R` named by them, for tests that build their own statistics and LD.
# `nSamples` is the reference size held-out selection reads.
.rrwIds <- function(p) {
    sprintf("chr1:%d:A:G", 100L * seq_len(p))
}

.rrwLdMatrix <- function(R, nSamples = NA_integer_) {
    ids <- .rrwIds(nrow(R))
    LdMatrix(`dimnames<-`(R, list(ids, ids)), nSamples = nSamples)
}

# A genotype panel's first `p` variants, the LD matrix it gives, and
# statistics with two strong effects, for the RSS weights that read either.
.rrwPanelStat <- function(p = 40) {
    panel <- readGenotypes(test_path("test_data", "test_variants.bed"))
    ids <- names(SummarizedExperiment::rowRanges(panel))[seq_len(p)]
    R <- pecotmr:::.ldFromRef(panel, ids, "test")
    set.seed(3)
    effects <- replace(numeric(p), c(3, 20), c(4, -3))
    z <- as.vector(R %*% effects) + rnorm(p)
    list(
        panel = panel,
        ldMatrix = LdMatrix(R, nSamples = ncol(panel)),
        stat = list(z = z, n = rep(500, p), variantNames = ids)
    )
}

# .rrwPanelStat()'s panel and LD, with statistics whose noise has the LD's
# covariance, as real marginal statistics do: the samplers that estimate a
# residual variance or diverge from them need statistics consistent with the
# LD.
.rrwPanelStatLd <- function() {
    d <- .rrwPanelStat()
    # Noise drawn through the factor, F'e / sqrt(k - 1), has covariance R
    # even where R is singular, as it is here.
    ldFactor <- pecotmr:::.ldFactorFromRef(d$panel, d$stat$variantNames, "test")
    k <- nrow(ldFactor)
    set.seed(3)
    effects <- replace(numeric(ncol(ldFactor)), c(3, 20), c(4, -3))
    signal <- crossprod(ldFactor, ldFactor %*% effects) / (k - 1)
    noise <- crossprod(ldFactor, rnorm(k)) / sqrt(k - 1)
    d$stat$z <- as.vector(signal + noise)
    d
}

# Two studies on disjoint halves of eqtlRegionExample's people: each one's
# marginal correlations and in-sample LD over 30 common variants, scaled to
# population SD, with an effect on four of them.
.rrwTwoStudies <- function() {
    data(eqtlRegionExample, envir = environment())
    X0 <- eqtlRegionExample$X
    p <- colMeans(X0) / 2
    cols <- which(pmin(p, 1 - p) >= 0.05)[301:330]
    withr::local_seed(5L)
    effects <- stats::rnorm(4L)
    map(list(1:250, 251:415), function(rows) {
        X <- scale(X0[rows, cols]) * sqrt(length(rows) / (length(rows) - 1))
        g <- drop(X[, c(3, 10, 17, 25)] %*% effects)
        y <- scale(g + stats::rnorm(length(rows)))[, 1L]
        list(
            X = X,
            z = as.numeric(crossprod(X, y)) / length(rows),
            R = crossprod(X) / length(rows),
            n = length(rows)
        )
    })
}

# .rrwTwoStudies()'s studies as the samplers across studies take them: each
# one's statistics (z-scores) and LdMatrix over the variants it keeps, the
# second lacking the last, with the union's ids.
.rrwStudyInputs <- function(studies, keep = list(1:30, 1:29)) {
    ids <- .rrwIds(30L)
    list(
        ids = ids,
        stats = map2(studies, keep, function(s, k) {
            list(
                z = (s$z * sqrt(s$n))[k],
                n = rep(s$n, length(k)),
                variantNames = ids[k]
            )
        }),
        refs = map2(studies, keep, function(s, k) {
            LdMatrix(
                `dimnames<-`(s$R[k, k], list(ids[k], ids[k])),
                nSamples = s$n
            )
        })
    )
}
