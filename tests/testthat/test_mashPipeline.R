# =============================================================================
# Tests for utility helpers exported from R/mashPipeline.R that don't
# require mashr / flashier installs:
#   sanitizeMashData, makePairwiseContrastCol, sliceMashData,
#   metaAnalysisPerCondition
# =============================================================================

# ---------------------------------------------------------------------------
# sanitizeMashData
# ---------------------------------------------------------------------------

test_that("sanitizeMashData replaces NaN in bhat with 0", {
    d <- list(
        bhat = matrix(c(1, NaN, 3, 4), nrow = 2),
        sbhat = matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2)
    )
    out <- sanitizeMashData(d)
    expect_equal(out$bhat[1, 2], 3) # untouched
    expect_equal(out$bhat[2, 1], 0) # NaN -> 0
    expect_equal(out$sbhat, d$sbhat) # sbhat untouched
})

test_that("sanitizeMashData replaces NaN/Inf in sbhat with 1e3", {
    d <- list(
        bhat = matrix(c(1, 2, 3, 4), nrow = 2),
        sbhat = matrix(c(0.1, NaN, Inf, 0.4), nrow = 2)
    )
    out <- sanitizeMashData(d)
    expect_equal(out$sbhat[2, 1], 1e3)
    expect_equal(out$sbhat[1, 2], 1e3)
    expect_equal(out$sbhat[1, 1], 0.1)
    expect_equal(out$sbhat[2, 2], 0.4)
    expect_equal(out$bhat, d$bhat) # bhat untouched (no NaN there)
})

test_that("sanitizeMashData is idempotent on already-clean data", {
    d <- list(
        bhat = matrix(c(1, 2, 3, 4), nrow = 2),
        sbhat = matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2)
    )
    expect_equal(sanitizeMashData(d), d)
    expect_equal(sanitizeMashData(sanitizeMashData(d)), d)
})

test_that("sanitizeMashData leaves -Inf in bhat alone (only NaN is replaced)", {
    d <- list(
        bhat = matrix(c(-Inf, Inf, 3, 4), nrow = 2),
        sbhat = matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2)
    )
    out <- sanitizeMashData(d)
    expect_true(is.infinite(out$bhat[1, 1]))
    expect_true(is.infinite(out$bhat[2, 1]))
})

# ---------------------------------------------------------------------------
# makePairwiseContrastCol
# ---------------------------------------------------------------------------

test_that("makePairwiseContrastCol sets +1/-1 at named pair positions", {
    tmpl <- setNames(rep(0, 4), c("a", "b", "c", "d"))
    out <- makePairwiseContrastCol(c("b", "d"), tmpl)
    expect_equal(out[["a"]], 0)
    expect_equal(out[["b"]], 1)
    expect_equal(out[["c"]], 0)
    expect_equal(out[["d"]], -1)
})

test_that("makePairwiseContrastCol preserves template names", {
    tmpl <- setNames(rep(0, 3), c("x", "y", "z"))
    out <- makePairwiseContrastCol(c("x", "z"), tmpl)
    expect_equal(names(out), c("x", "y", "z"))
})

test_that("makePairwiseContrastCol overwrites pre-existing values in template", {
    tmpl <- setNames(c(5, -3, 7), c("a", "b", "c"))
    out <- makePairwiseContrastCol(c("a", "c"), tmpl)
    expect_equal(out[["a"]], 1) # was 5, now 1
    expect_equal(out[["b"]], -3) # untouched
    expect_equal(out[["c"]], -1) # was 7, now -1
})

# ---------------------------------------------------------------------------
# sliceMashData
# ---------------------------------------------------------------------------

test_that("sliceMashData subsets bhat / sbhat / Z by SNP and sample", {
    snps <- c("s1", "s2", "s3")
    samples <- c("ctxA", "ctxB", "ctxC")
    data <- list(
        bhat = matrix(seq_len(9), 3, 3, dimnames = list(snps, samples)),
        sbhat = matrix(seq_len(9) / 10, 3, 3, dimnames = list(snps, samples)),
        Z = matrix(seq_len(9) * 2, 3, 3, dimnames = list(snps, samples)),
        snp = snps
    )
    vhat <- diag(1, 3, 3)
    dimnames(vhat) <- list(samples, samples)

    out <- sliceMashData(
        data,
        vhat,
        snps = c("s1", "s3"),
        samples = c("ctxA", "ctxC")
    )
    expect_equal(dim(out$data$bhat), c(2, 2))
    expect_equal(dim(out$vhat), c(2, 2))
    expect_equal(colnames(out$data$bhat), c("ctxA", "ctxC"))
    expect_equal(colnames(out$data$sbhat), c("ctxA", "ctxC"))
    expect_equal(colnames(out$data$Z), c("ctxA", "ctxC"))
    expect_equal(colnames(out$vhat), c("ctxA", "ctxC"))
    expect_equal(out$data$snp, c("s1", "s3"))
})

test_that("sliceMashData restricts data$snp to intersection of snps argument", {
    snps <- c("s1", "s2", "s3", "s4")
    samples <- c("ctxA", "ctxB")
    data <- list(
        bhat = matrix(1, 4, 2, dimnames = list(snps, samples)),
        sbhat = matrix(1, 4, 2, dimnames = list(snps, samples)),
        Z = matrix(1, 4, 2, dimnames = list(snps, samples)),
        snp = snps
    )
    vhat <- diag(1, 2, 2)
    dimnames(vhat) <- list(samples, samples)
    out <- sliceMashData(data, vhat, snps = c("s2", "s4"), samples = samples)
    expect_equal(out$data$snp, c("s2", "s4"))
})

# ---------------------------------------------------------------------------
# metaAnalysisPerCondition
# ---------------------------------------------------------------------------

test_that("metaAnalysisPerCondition returns single-effect p-value when only one feature passes filter", {
    feat <- "var1"
    cols <- c("mean_contrast_brain_vs_blood")
    es <- matrix(0.5, nrow = 1, ncol = 1, dimnames = list(feat, cols))
    se <- matrix(0.1, nrow = 1, ncol = 1, dimnames = list(feat, cols))
    out <- metaAnalysisPerCondition(es, se)
    # 2 conditions: brain, blood -> 2 rows
    expect_equal(nrow(out), 2L)
    expect_true(all(
        c(
            "condition",
            "contrast",
            "meta_pvalue",
            "meta_effect",
            "meta_se",
            "tau2",
            "I2"
        ) %in%
            names(out)
    ))
    # With a single effect both rows return single-effect p-values (not NA)
    expect_false(any(is.na(out$meta_pvalue)))
    # tau2 / I2 only meaningful with >= 2 effects
    expect_true(all(is.na(out$tau2)))
    expect_true(all(is.na(out$I2)))
})

test_that("metaAnalysisPerCondition returns NA pvalue when SE cutoff drops everything", {
    feat <- c("v1", "v2")
    cols <- "mean_contrast_brain_vs_blood"
    es <- matrix(c(0.1, 0.2), nrow = 2, ncol = 1, dimnames = list(feat, cols))
    se <- matrix(c(0.01, 0.02), nrow = 2, ncol = 1, dimnames = list(feat, cols))
    # seCutoff = 0.5 drops both rows
    out <- metaAnalysisPerCondition(es, se, seCutoff = 0.5)
    expect_true(all(is.na(out$meta_pvalue)))
    expect_true(all(is.na(out$meta_effect)))
})

test_that("metaAnalysisPerCondition runs DerSimonian-Laird when >=2 effects survive", {
    feat <- c("v1", "v2", "v3")
    cols <- "mean_contrast_brain_vs_blood"
    es <- matrix(
        c(0.3, 0.5, 0.4),
        nrow = 3,
        ncol = 1,
        dimnames = list(feat, cols)
    )
    se <- matrix(
        c(0.1, 0.1, 0.1),
        nrow = 3,
        ncol = 1,
        dimnames = list(feat, cols)
    )
    out <- metaAnalysisPerCondition(es, se)
    expect_equal(nrow(out), 2L)
    # All three effects survive (SE = 0.1 > 0): meta_effect and meta_se populated
    expect_false(any(is.na(out$meta_effect)))
    expect_false(any(is.na(out$meta_se)))
    expect_false(any(is.na(out$tau2)))
    expect_false(any(is.na(out$I2)))
    # I2 in [0, 1]
    expect_true(all(out$I2 >= 0 & out$I2 <= 1))
})

test_that("metaAnalysisPerCondition unique-condition extraction handles >2 conditions", {
    feat <- "v1"
    cols <- c(
        "mean_contrast_brain_vs_blood",
        "mean_contrast_brain_vs_muscle",
        "mean_contrast_blood_vs_muscle"
    )
    es <- matrix(c(0.3, 0.4, 0.1), nrow = 1, dimnames = list(feat, cols))
    se <- matrix(c(0.1, 0.1, 0.1), nrow = 1, dimnames = list(feat, cols))
    out <- metaAnalysisPerCondition(es, se)
    # 3 conditions (brain, blood, muscle), each with 2 vs-comparisons -> 6 rows
    expect_setequal(unique(out$condition), c("brain", "blood", "muscle"))
    expect_equal(nrow(out), 6L)
})

# ---------------------------------------------------------------------------
# updateMashModelCov — works on a hand-built mock fitted_g (no mashr dep)
# ---------------------------------------------------------------------------

test_that("updateMashModelCov drops dropped conditions + resizes remaining cov matrices", {
    R <- 3L
    samples <- c("brain", "blood", "muscle")
    U <- list(
        brain = diag(c(1, 0, 0)),
        blood = diag(c(0, 1, 0)),
        muscle = diag(c(0, 0, 1)),
        identity = diag(1, R),
        PCA_1 = matrix(seq_len(R * R), R, R)
    ) # no dimnames -> last branch
    pi <- setNames(
        rep(0.2, 5L),
        c(
            "brain.scale1",
            "blood.scale1",
            "muscle.scale1",
            "identity.scale1",
            "PCA_1.scale1"
        )
    )
    m <- list(fitted_g = list(Ulist = U, pi = pi))
    m2 <- updateMashModelCov(
        m,
        allSamples = samples,
        samples = c("brain", "blood")
    )
    expect_false("muscle" %in% names(m2$fitted_g$Ulist))
    expect_true(all(vapply(
        m2$fitted_g$Ulist,
        function(x) all(dim(x) == c(2L, 2L)),
        logical(1)
    )))
    expect_false(any(grepl("muscle", names(m2$fitted_g$pi))))
    # Brain matrix has a single 1 at the brain position (the first of the
    # retained `samples` ordering).
    expect_equal(m2$fitted_g$Ulist$brain[1, 1], 1)
    expect_equal(sum(m2$fitted_g$Ulist$brain), 1)
})

# ---------------------------------------------------------------------------
# fitMashContrast — fabricated posterior inputs (no mashr dep)
# ---------------------------------------------------------------------------

test_that("fitMashContrast returns NULL when fewer than 2 tested conditions", {
    origMean <- matrix(
        0,
        nrow = 1,
        ncol = 3,
        dimnames = list("v1", c("a", "b", "c"))
    )
    origMean[1, "b"] <- 1
    pm <- matrix(0, nrow = 1, ncol = 3, dimnames = list("v1", c("a", "b", "c")))
    pv <- array(diag(3), dim = c(3, 3, 1))
    dimnames(pv) <- list(c("a", "b", "c"), c("a", "b", "c"), NULL)
    expect_null(fitMashContrast(1L, origMean, pm, pv))
})

test_that("fitMashContrast: 2-tested-conditions fast path yields one pairwise contrast", {
    origMean <- matrix(
        c(0.5, 0.3, 0),
        nrow = 1,
        dimnames = list("v1", c("a", "b", "c"))
    )
    pm <- matrix(
        c(0.5, 0.3, 0),
        nrow = 1,
        dimnames = list("v1", c("a", "b", "c"))
    )
    pv <- array(diag(3) * 0.1, dim = c(3, 3, 1))
    dimnames(pv) <- list(c("a", "b", "c"), c("a", "b", "c"), NULL)
    out <- fitMashContrast(1L, origMean, pm, pv)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 1L)
    # feature_id + 2 tested -> 1 pairwise contrast (mean, se, p) = 4 columns
    expect_equal(ncol(out), 4L)
    expect_setequal(
        names(out),
        c(
            "feature_id",
            "mean_contrast_a_vs_b",
            "se_contrast_a_vs_b",
            "p_contrast_a_vs_b"
        )
    )
    expect_equal(out[["mean_contrast_a_vs_b"]], 0.5 - 0.3)
})

test_that("fitMashContrast: 3-tested-conditions yields deviation + pairwise contrasts", {
    origMean <- matrix(
        c(0.5, 0.3, -0.2),
        nrow = 1,
        dimnames = list("v1", c("a", "b", "c"))
    )
    pm <- matrix(
        c(0.5, 0.3, -0.2),
        nrow = 1,
        dimnames = list("v1", c("a", "b", "c"))
    )
    pv <- array(diag(3) * 0.1, dim = c(3, 3, 1))
    dimnames(pv) <- list(c("a", "b", "c"), c("a", "b", "c"), NULL)
    out <- fitMashContrast(1L, origMean, pm, pv)
    expect_s3_class(out, "data.frame")
    # feature_id + 3 deviation + choose(3,2)=3 pairwise = 6 contrasts -> 19 cols
    expect_equal(ncol(out), 19L)
    contrastSuffix <- sub("^(mean|se|p)_contrast_", "", names(out))
    expect_true(any(grepl("_deviation$", contrastSuffix)))
    expect_true(any(grepl("_vs_", contrastSuffix)))
})

# ---------------------------------------------------------------------------
# mashPipeline — end-to-end on the bundled multi-context example
# ---------------------------------------------------------------------------

test_that("mashPipeline runs end-to-end on qtlSumStatsMulticontextExample", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    # Use the same fixture for strong/random; nPcs <= ncol - 1 (3 contexts).
    res <- suppressMessages(suppressWarnings(
        mashPipeline(
            sumStatsList = list(strong = ss, random = ss),
            alpha = 0,
            nPcs = 2L,
            setSeed = 1L
        )
    ))
    expect_named(res, c("U", "w"))
    expect_type(res$U, "list")
    expect_gt(length(res$U), 0L)
    # Every covariance matrix is 3x3 (one row/col per context)
    expect_true(all(vapply(
        res$U,
        function(m) all(dim(m) == c(3L, 3L)),
        logical(1)
    )))
    expect_type(res$w, "double")
    expect_equal(sum(res$w), 1, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# fitMashContrast — condition grouping (>2 conditions, grouped replicates)
# ---------------------------------------------------------------------------

test_that("fitMashContrast applies deviation + pairwise group adjustments", {
    conds <- c("a", "b", "c", "d")
    origMean <- matrix(
        c(0.5, 0.3, -0.2, 0.4),
        nrow = 1,
        dimnames = list("v1", conds)
    )
    pm <- matrix(c(0.5, 0.3, -0.2, 0.4), nrow = 1, dimnames = list("v1", conds))
    pv <- array(0, dim = c(4, 4, 1), dimnames = list(conds, conds, NULL))
    pv[,, 1] <- diag(4) * 0.1
    # a,b share group 1 (replicates); c is its own group 2; d ungrouped (0).
    # Non-NULL grouping triggers `grouping <- grouping[tested]`, the >2-condition
    # deviation re-weighting loop, and the pairwise group-adjustment loop.
    grouping <- setNames(c(1L, 1L, 2L, 0L), conds)
    out <- fitMashContrast(1L, origMean, pm, pv, grouping = grouping)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 1L)
    # feature_id + 4 deviation + choose(4,2)=6 pairwise = 10 contrasts -> 31 cols.
    expect_equal(ncol(out), 31L)
    contrastSuffix <- sub("^(mean|se|p)_contrast_", "", names(out))
    expect_true(any(grepl("_deviation$", contrastSuffix)))
    expect_true(any(grepl("_vs_", contrastSuffix)))
    contrastVals <- unlist(out[grepl("_contrast_", names(out))])
    expect_true(all(is.finite(contrastVals)))
})

# ---------------------------------------------------------------------------
# updateMashModelCov — named data-driven cov matrices (the `[samples, samples]`
# else branch, distinct from the no-dimnames positional-slice branch)
# ---------------------------------------------------------------------------

test_that("updateMashModelCov slices named data-driven cov matrices by sample", {
    allSamples <- c("brain", "blood", "muscle")
    ddMat <- matrix(seq_len(9), 3, 3, dimnames = list(allSamples, allSamples))
    U <- list(identity = diag(1, 3), dataDriven = ddMat)
    pi <- setNames(c(0.5, 0.5), c("identity.scale1", "dataDriven.scale1"))
    m <- list(fitted_g = list(Ulist = U, pi = pi))
    m2 <- updateMashModelCov(
        m,
        allSamples = allSamples,
        samples = c("brain", "muscle")
    )
    # The named data-driven matrix is sliced by name: cov[[d]][samples, samples].
    expect_equal(dim(m2$fitted_g$Ulist$dataDriven), c(2L, 2L))
    expect_equal(
        m2$fitted_g$Ulist$dataDriven,
        ddMat[c("brain", "muscle"), c("brain", "muscle")]
    )
    # identity collapses to a single 1 in the top-left corner.
    expect_equal(m2$fitted_g$Ulist$identity[1, 1], 1)
    expect_equal(sum(m2$fitted_g$Ulist$identity), 1)
})

# ---------------------------------------------------------------------------
# metaAnalysisPerCondition — conditions matching no contrast column are skipped
# ---------------------------------------------------------------------------

test_that("metaAnalysisPerCondition skips conditions whose name matches no column", {
    # The condition "x$y" is a derived condition name; used as a grep() pattern
    # the embedded `$` anchor matches nothing, exercising the
    # `if (length(idx) == 0) next` skip branch.
    cols <- "mean_contrast_x$y_vs_z"
    es <- matrix(c(0.3, 0.5), nrow = 2, dimnames = list(c("v1", "v2"), cols))
    se <- matrix(c(0.1, 0.1), nrow = 2, dimnames = list(c("v1", "v2"), cols))
    out <- metaAnalysisPerCondition(es, se)
    # "x$y" is skipped; only the "z" condition survives.
    expect_false("x$y" %in% out$condition)
    expect_true("z" %in% out$condition)
    expect_equal(nrow(out), 1L)
})

# ---------------------------------------------------------------------------
# mashPipeline — input validation (errors fire before any mashr call).
# Guarded by skips because the requireNamespace() checks run first; without
# mashr/flashier the function stops with an install message instead.
# ---------------------------------------------------------------------------

test_that("mashPipeline rejects a sumStatsList that is not a named list", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    expect_error(mashPipeline(list(1, 2), alpha = 0), "must be a named list")
})

test_that("mashPipeline errors when a required entry is missing", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    expect_error(
        mashPipeline(list(strong = 1), alpha = 0),
        "missing required entr"
    )
})

test_that("mashPipeline errors on unrecognised entries", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    expect_error(
        mashPipeline(list(strong = 1, random = 1, bogus = 1), alpha = 0),
        "unrecognised entries"
    )
})

test_that("mashPipeline coerces a SimpleList before validating its names", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    skip_if_not_installed("S4Vectors")
    # A SimpleList is converted to a base list first (the as.list branch), then
    # validation runs; here it is missing both required entries.
    expect_error(
        mashPipeline(S4Vectors::SimpleList(bogus = 1), alpha = 0),
        "missing required entr"
    )
})

# ---------------------------------------------------------------------------
# mashPipeline — priorCovariances validation + bypass path.
# Supplying residualCorrelation makes random/null optional and short-circuits
# null-correlation estimation, so the supplied Vhat branch is also exercised.
# ---------------------------------------------------------------------------

test_that("mashPipeline rejects priorCovariances not a non-empty named list", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    vhat <- diag(3)
    # Empty list.
    expect_error(
        suppressMessages(suppressWarnings(
            mashPipeline(
                list(strong = ss),
                alpha = 0,
                residualCorrelation = vhat,
                priorCovariances = list()
            )
        )),
        "non-empty named"
    )
    # Unnamed list.
    expect_error(
        suppressMessages(suppressWarnings(
            mashPipeline(
                list(strong = ss),
                alpha = 0,
                residualCorrelation = vhat,
                priorCovariances = list(diag(3))
            )
        )),
        "non-empty named"
    )
})

test_that("mashPipeline rejects priorCovariances with wrong dimensions", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    vhat <- diag(3)
    expect_error(
        suppressMessages(suppressWarnings(
            mashPipeline(
                list(strong = ss),
                alpha = 0,
                residualCorrelation = vhat,
                priorCovariances = list(myU = diag(2))
            )
        )),
        "3 x 3 matrix"
    )
})

test_that("mashPipeline passes supplied residualCorrelation + priorCovariances through", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    vhat <- diag(3)
    U0 <- list(identity = diag(3), effectA = diag(c(1, 0, 0)))
    res <- suppressMessages(suppressWarnings(
        mashPipeline(
            list(strong = ss),
            alpha = 0,
            residualCorrelation = vhat,
            priorCovariances = U0
        )
    ))
    expect_named(res, c("U", "w"))
    # priorCovariances passed straight through as the covariance list (bypass of
    # the cov_canonical / cov_pca / cov_flash / cov_ed chain).
    expect_identical(res$U, U0)
    expect_type(res$w, "double")
    expect_equal(sum(res$w), 1, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# mashPipeline — null-based Vhat estimation + default nPcs
# ---------------------------------------------------------------------------

test_that("mashPipeline estimates Vhat from a null set and defaults nPcs", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    # Supplying `null` triggers estimate_null_correlation_simple; leaving nPcs
    # NULL exercises the `nPcs <- ncol(Bhat) - 1` default in the cov_* chain.
    res <- suppressMessages(suppressWarnings(
        mashPipeline(
            list(strong = ss, random = ss, null = ss),
            alpha = 0,
            setSeed = 1L
        )
    ))
    expect_named(res, c("U", "w"))
    expect_gt(length(res$U), 0L)
    expect_true(all(vapply(
        res$U,
        function(m) all(dim(m) == c(3L, 3L)),
        logical(1)
    )))
    expect_equal(sum(res$w), 1, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# mashResidualCorrelation — the Vhat estimator extracted from mashPipeline.
# ---------------------------------------------------------------------------

test_that("mashResidualCorrelation(identity) is an identity of the right size", {
    skip_if_not_installed("mashr")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    v <- mashResidualCorrelation(
        list(strong = ss),
        alpha = 0,
        method = "identity"
    )
    expect_equal(dim(v), c(3L, 3L))
    expect_equal(v, diag(3))
})

test_that("mashResidualCorrelation(simple) returns a null correlation matrix", {
    skip_if_not_installed("mashr")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    v <- suppressMessages(suppressWarnings(
        mashResidualCorrelation(
            list(strong = ss, null = ss),
            alpha = 0,
            method = "simple"
        )
    ))
    expect_equal(dim(v), c(3L, 3L))
    expect_equal(unname(diag(v)), rep(1, 3), tolerance = 1e-8) # a correlation matrix
    expect_true(isSymmetric(unname(v)))
})

test_that("mashResidualCorrelation(simple) errors without a null entry", {
    skip_if_not_installed("mashr")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    expect_error(
        mashResidualCorrelation(
            list(strong = ss),
            alpha = 0,
            method = "simple"
        ),
        "requires a 'null' entry"
    )
})

test_that("mashResidualCorrelation(simple_specific) returns a null correlation", {
    skip_if_not_installed("mashr")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    v <- suppressMessages(suppressWarnings(
        mashResidualCorrelation(
            list(strong = ss, null = ss),
            alpha = 0,
            method = "simple_specific"
        )
    ))
    expect_equal(dim(v), c(3L, 3L))
    expect_equal(unname(diag(v)), rep(1, 3), tolerance = 1e-6)
})

test_that("mashResidualCorrelation(corshrink) returns a 3x3 correlation matrix", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("CorShrink")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    v <- suppressMessages(suppressWarnings(
        mashResidualCorrelation(
            list(strong = ss, null = ss),
            alpha = 0,
            method = "corshrink"
        )
    ))
    expect_equal(dim(v), c(3L, 3L))
    expect_equal(unname(diag(v)), rep(1, 3), tolerance = 1e-6)
})

test_that("mashResidualCorrelation(mle) refines V against a supplied prior", {
    skip_if_not_installed("mashr")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    v <- suppressMessages(suppressWarnings(
        mashResidualCorrelation(
            list(strong = ss, random = ss),
            alpha = 0,
            method = "mle",
            priorCovariances = list(identity = diag(3)),
            nSubset = 100L,
            maxIter = 3L,
            setSeed = 1L
        )
    ))
    expect_equal(dim(v), c(3L, 3L))
})

test_that("mashResidualCorrelation errors when a method's inputs are missing", {
    skip_if_not_installed("mashr")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    expect_error(
        mashResidualCorrelation(
            list(strong = ss),
            alpha = 0,
            method = "corshrink"
        ),
        "requires a 'null'"
    )
    expect_error(
        mashResidualCorrelation(
            list(strong = ss, random = ss),
            alpha = 0,
            method = "mle"
        ),
        "priorCovariances"
    )
})

# ---------------------------------------------------------------------------
# mashPriorCovariances — the covariance + weight estimator extracted from
# mashPipeline. Default builds every non-udr component (canonical + pca +
# flash + flash_nonneg) refined by cov_ed.
# ---------------------------------------------------------------------------

test_that("mashPriorCovariances computes the default (all-but-udr) prior", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    pc <- suppressMessages(suppressWarnings(
        mashPriorCovariances(
            list(strong = ss),
            alpha = 0,
            vhat = diag(3),
            nPcs = 2L,
            setSeed = 1L
        )
    ))
    expect_named(pc, c("U", "w", "loglik"))
    expect_gt(length(pc$U), 0L)
    expect_true(all(vapply(
        pc$U,
        function(m) all(dim(m) == c(3L, 3L)),
        logical(1)
    )))
    expect_equal(sum(pc$w), 1, tolerance = 1e-6)
    expect_null(pc$loglik)
})

test_that("mashPriorCovariances passes a supplied prior through unchanged", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    U0 <- list(identity = diag(3), effectA = diag(c(1, 0, 0)))
    pc <- suppressMessages(suppressWarnings(
        mashPriorCovariances(
            list(strong = ss),
            alpha = 0,
            vhat = diag(3),
            priorCovariances = U0
        )
    ))
    expect_identical(pc$U, U0)
    expect_equal(sum(pc$w), 1, tolerance = 1e-6)
})

test_that("mashPriorCovariances validates a supplied prior", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    expect_error(
        suppressMessages(suppressWarnings(
            mashPriorCovariances(
                list(strong = ss),
                alpha = 0,
                vhat = diag(3),
                priorCovariances = list()
            )
        )),
        "non-empty named"
    )
    expect_error(
        suppressMessages(suppressWarnings(
            mashPriorCovariances(
                list(strong = ss),
                alpha = 0,
                vhat = diag(3),
                priorCovariances = list(myU = diag(2))
            )
        )),
        "3 x 3 matrix"
    )
})

test_that("mashPriorCovariances(flash_nonneg) adds components vs flash-only", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    base <- suppressMessages(suppressWarnings(
        mashPriorCovariances(
            list(strong = ss),
            alpha = 0,
            vhat = diag(3),
            components = c("canonical", "pca", "flash"),
            nPcs = 2L,
            setSeed = 1L
        )
    ))
    wide <- suppressMessages(suppressWarnings(
        mashPriorCovariances(
            list(strong = ss),
            alpha = 0,
            vhat = diag(3),
            components = c("canonical", "pca", "flash", "flash_nonneg"),
            nPcs = 2L,
            setSeed = 1L
        )
    ))
    expect_gt(length(wide$U), length(base$U))
})

test_that("mashPriorCovariances engine 'ud' (udr) produces U + weights", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    skip_if_not_installed("udr")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    # A toy-sized udr config: `n_unconstrained` dominates the cost, and the default
    # (50, sized for many-condition data) is pathological on a 3-condition fixture
    # (it drove a 5+ minute fit). 2 unconstrained matrices suffice to exercise the
    # udr path.
    pc <- suppressMessages(suppressWarnings(
        mashPriorCovariances(
            list(strong = ss),
            alpha = 0,
            vhat = diag(3),
            engine = "ud",
            udControl = list(n_unconstrained = 2L, maxiter = 20L),
            setSeed = 1L
        )
    ))
    expect_gt(length(pc$U), 0L)
    expect_true(all(vapply(
        pc$U,
        function(m) all(dim(m) == c(3L, 3L)),
        logical(1)
    )))
    expect_equal(sum(pc$w), 1, tolerance = 1e-6)
})

test_that("mashPriorCovariances engine 'ud_ted' errors clearly on non-i.i.d. data", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    skip_if_not_installed("udr")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    expect_error(
        suppressMessages(suppressWarnings(
            mashPriorCovariances(
                list(strong = ss),
                alpha = 0,
                vhat = diag(3),
                engine = "ud_ted",
                setSeed = 1L
            )
        )),
        "i.i.d"
    )
})

test_that("mashPriorCovariances rejects an unknown component", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    expect_error(
        mashPriorCovariances(
            list(strong = ss),
            alpha = 0,
            components = "bogus"
        ),
        "unknown component"
    )
})

# ---------------------------------------------------------------------------
# mashCovarianceComponents — the raw per-method component builder that
# mashPriorCovariances refines (and the mixture-prior notebook demonstrates).
# ---------------------------------------------------------------------------

test_that("mashCovarianceComponents builds a single requested component", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    fl <- suppressMessages(suppressWarnings(
        mashCovarianceComponents(
            list(strong = ss),
            alpha = 0,
            vhat = diag(3),
            components = "flash",
            setSeed = 1L
        )
    ))
    expect_gt(length(fl), 0L)
    expect_true(all(vapply(
        fl,
        function(m) all(dim(m) == c(3L, 3L)),
        logical(1)
    )))
})

test_that("mashCovarianceComponents default builds all non-udr components", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    one <- suppressMessages(suppressWarnings(
        mashCovarianceComponents(
            list(strong = ss),
            alpha = 0,
            vhat = diag(3),
            components = "canonical",
            setSeed = 1L
        )
    ))
    all4 <- suppressMessages(suppressWarnings(
        mashCovarianceComponents(
            list(strong = ss),
            alpha = 0,
            vhat = diag(3),
            nPcs = 2L,
            setSeed = 1L
        )
    ))
    expect_gt(length(all4), length(one))
})

test_that("mashCovarianceComponents feeds mashPriorCovariances (same components)", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    comps <- suppressMessages(suppressWarnings(
        mashCovarianceComponents(
            list(strong = ss),
            alpha = 0,
            vhat = diag(3),
            nPcs = 2L,
            setSeed = 1L
        )
    ))
    prior <- suppressMessages(suppressWarnings(
        mashPriorCovariances(
            list(strong = ss),
            alpha = 0,
            vhat = diag(3),
            nPcs = 2L,
            setSeed = 1L
        )
    ))
    # mashPriorCovariances refines these components, so its U contains them all.
    expect_true(all(names(comps) %in% names(prior$U)))
})

test_that("mashCovarianceComponents rejects unknown components", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    expect_error(
        mashCovarianceComponents(
            list(strong = ss),
            alpha = 0,
            components = "bogus"
        ),
        "unknown component"
    )
})

test_that("mashPriorCovariances refines supplied priorComponents (pipeline mode)", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    comps <- suppressMessages(suppressWarnings(
        mashCovarianceComponents(
            list(strong = ss),
            alpha = 0,
            vhat = diag(3),
            nPcs = 2L,
            setSeed = 1L
        )
    ))
    pr <- suppressMessages(suppressWarnings(
        mashPriorCovariances(
            list(strong = ss),
            alpha = 0,
            vhat = diag(3),
            priorComponents = comps,
            engine = "cov_ed",
            setSeed = 1L
        )
    ))
    # the supplied components are refined into the returned U (not rebuilt)
    expect_true(all(names(comps) %in% names(pr$U)))
    expect_equal(sum(pr$w), 1, tolerance = 1e-6)
})

test_that("mashPriorCovariances validates priorComponents", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    expect_error(
        mashPriorCovariances(
            list(strong = ss),
            alpha = 0,
            priorComponents = list()
        ),
        "non-empty named"
    )
})

test_that("mashPipeline result == composing the two extracted building blocks", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    full <- suppressMessages(suppressWarnings(
        mashPipeline(
            list(strong = ss, random = ss, null = ss),
            alpha = 0,
            setSeed = 1L
        )
    ))
    # Same seed discipline mashPipeline uses: seed once, then delegate with
    # setSeed = NULL so the RNG stream stays continuous.
    set.seed(1L)
    vhat <- mashResidualCorrelation(
        list(strong = ss, null = ss),
        alpha = 0,
        method = "simple",
        setSeed = NULL
    )
    prior <- suppressMessages(suppressWarnings(
        mashPriorCovariances(
            list(strong = ss),
            alpha = 0,
            vhat = vhat,
            setSeed = NULL
        )
    ))
    expect_equal(full$U, prior$U)
    expect_equal(full$w, prior$w)
})

# ---------------------------------------------------------------------------
# mashModelFit + mashPosterior — the fit -> posterior chain (mash_fit /
# mash_posterior). A tiny 2-component prior keeps these fast.
# ---------------------------------------------------------------------------

.mashTestModel <- function(ss) {
    suppressMessages(suppressWarnings(
        mashModelFit(
            list(random = ss),
            alpha = 0,
            priorCovariances = list(
                identity = diag(3),
                effectA = diag(c(1, 0, 0))
            ),
            vhat = diag(3),
            setSeed = 1L
        )
    ))
}

test_that("mashModelFit returns a fitted mash model", {
    skip_if_not_installed("mashr")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    m <- .mashTestModel(ss)
    expect_s3_class(m, "mash")
    expect_false(is.null(m$fitted_g))
})

test_that("mashModelFit validates the prior and the fitOn entry", {
    skip_if_not_installed("mashr")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    expect_error(
        mashModelFit(list(random = ss), alpha = 0, priorCovariances = list()),
        "non-empty named list"
    )
    expect_error(
        mashModelFit(
            list(strong = ss),
            alpha = 0,
            priorCovariances = list(identity = diag(3))
        ),
        "no 'random' entry"
    )
})

test_that("mashPosterior returns posterior matrices with covariance", {
    skip_if_not_installed("mashr")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    post <- suppressMessages(suppressWarnings(
        mashPosterior(.mashTestModel(ss), ss, alpha = 0, vhat = diag(3))
    ))
    expect_true(all(
        c("PosteriorMean", "PosteriorSD", "lfsr", "PosteriorCov") %in%
            names(post)
    ))
    expect_equal(ncol(post$PosteriorMean), 3L)
    expect_equal(nrow(post$PosteriorMean), nrow(post$PosteriorSD))
})

test_that("mashPosterior outputPosteriorCov = FALSE omits PosteriorCov", {
    skip_if_not_installed("mashr")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    post <- suppressMessages(suppressWarnings(
        mashPosterior(
            .mashTestModel(ss),
            ss,
            alpha = 0,
            vhat = diag(3),
            outputPosteriorCov = FALSE
        )
    ))
    expect_false("PosteriorCov" %in% names(post))
    expect_equal(ncol(post$PosteriorMean), 3L)
})

test_that("mashPosterior(excludeCondition) drops the condition from model + output", {
    skip_if_not_installed("mashr")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    conds <- colnames(.mashSumStatsToMatrices(ss, "strong")$b)
    post <- suppressMessages(suppressWarnings(
        mashPosterior(
            .mashTestModel(ss),
            ss,
            alpha = 0,
            vhat = diag(3),
            excludeCondition = conds[3]
        )
    ))
    expect_equal(ncol(post$PosteriorMean), 2L)
    expect_false(conds[3] %in% colnames(post$PosteriorMean))
})

test_that("mashPosterior errors on an unknown excludeCondition", {
    skip_if_not_installed("mashr")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    expect_error(
        mashPosterior(
            .mashTestModel(ss),
            ss,
            alpha = 0,
            vhat = diag(3),
            excludeCondition = "not_a_condition"
        ),
        "not found in the target"
    )
})

test_that("fitMashContrast consumes a mashPosterior result", {
    skip_if_not_installed("mashr")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    post <- suppressMessages(suppressWarnings(
        mashPosterior(.mashTestModel(ss), ss, alpha = 0, vhat = diag(3))
    ))
    # Force all 3 conditions "tested" (non-zero orig) so a contrast is returned.
    om <- post$PosteriorMean
    om[om == 0] <- 0.01
    fc <- fitMashContrast(1L, om, post$PosteriorMean, post$PosteriorCov)
    expect_s3_class(fc, "data.frame")
    expect_equal(nrow(fc), 1L)
})

# ---------------------------------------------------------------------------
# .rmaMeta: thin metafor::rma adapter (DL default; REML/ML/... pass through)
# ---------------------------------------------------------------------------

test_that(".rmaMeta: DL is the default and reshapes the metafor fit", {
    set.seed(1)
    m <- rnorm(6)
    s <- abs(rnorm(6)) + 0.2
    expect_equal(
        pecotmr:::.rmaMeta(m, s),
        pecotmr:::.rmaMeta(m, s, method = "DL")
    )
    out <- pecotmr:::.rmaMeta(m, s)
    expect_true(all(c("mean", "se", "tau2", "I2", "Q") %in% names(out)))
    expect_true(out$I2 >= 0 && out$I2 <= 1) # metafor's % rescaled to [0,1]
})

test_that(".rmaMeta: DL matches the closed-form DerSimonian-Laird", {
    means <- c(0.5, 0.8, 0.3)
    ses <- c(0.2, 0.3, 0.15)
    res <- pecotmr:::.rmaMeta(means, ses)
    wFe <- 1 / ses^2
    muFe <- sum(wFe * means) / sum(wFe)
    Q <- sum(wFe * (means - muFe)^2)
    tau2 <- max(0, (Q - 2) / (sum(wFe) - sum(wFe^2) / sum(wFe)))
    wRe <- 1 / (ses^2 + tau2)
    expect_equal(res$Q, Q, tolerance = 1e-6)
    expect_equal(res$tau2, tau2, tolerance = 1e-6)
    expect_equal(res$mean, sum(wRe * means) / sum(wRe), tolerance = 1e-6)
    expect_equal(res$se, sqrt(1 / sum(wRe)), tolerance = 1e-6)
})

test_that(".rmaMeta: non-DL estimator (REML) runs via metafor and differs from DL", {
    set.seed(3)
    m <- rnorm(8, sd = 1.5)
    s <- abs(rnorm(8)) + 0.3
    reml <- pecotmr:::.rmaMeta(m, s, method = "REML")
    expect_true(all(c("mean", "se", "tau2", "I2", "Q") %in% names(reml)))
    expect_true(reml$I2 >= 0 && reml$I2 <= 1)
    # REML tau2 generally differs from DL tau2 on heterogeneous data
    expect_false(isTRUE(all.equal(reml$tau2, pecotmr:::.rmaMeta(m, s)$tau2)))
})

test_that(".rmaMeta: non-DL estimator falls back to DL on non-convergence (never errors)", {
    # metafor's iterative estimators can fail to converge on small / degenerate
    # inputs; .rmaMeta must still return a finite fit (via the DL fallback).
    set.seed(1)
    for (i in seq_len(25)) {
        m <- rnorm(sample(3:8, 1L))
        s <- abs(rnorm(length(m))) + 0.1
        r <- suppressWarnings(pecotmr:::.rmaMeta(m, s, method = "REML"))
        expect_true(is.finite(r$mean) && is.finite(r$se))
    }
})

test_that("metaAnalysisPerCondition: threads metaMethod (DL default unchanged)", {
    es <- matrix(
        rnorm(12),
        6,
        2,
        dimnames = list(NULL, c("mean_contrast_A_vs_B", "mean_contrast_A_vs_C"))
    )
    sv <- matrix(abs(rnorm(12)) + 0.1, 6, 2, dimnames = dimnames(es))
    expect_equal(
        metaAnalysisPerCondition(es, sv),
        metaAnalysisPerCondition(es, sv, metaMethod = "DL")
    )
})

# ---------------------------------------------------------------------------
# mashPosteriorContrast (orchestrates fitMashContrast over all features)
# ---------------------------------------------------------------------------

.mpc_fixture <- function(nf = 5, conds = c("Ast", "Mic", "Oli"), seed = 1) {
    set.seed(seed)
    pm <- matrix(
        rnorm(nf * length(conds)),
        nf,
        length(conds),
        dimnames = list(paste0("chr1:", seq_len(nf), ":A:G"), conds)
    )
    pv <- array(
        0,
        c(length(conds), length(conds), nf),
        dimnames = list(conds, conds, NULL)
    )
    for (i in seq_len(nf)) {
        A <- matrix(rnorm(length(conds)^2), length(conds))
        pv[,, i] <- crossprod(A)
    }
    orig <- matrix(
        rnorm(nf * length(conds)),
        nf,
        length(conds),
        dimnames = list(rownames(pm), conds)
    )
    list(pm = pm, pv = pv, orig = orig)
}

test_that("mashPosteriorContrast: deviation + pairwise columns, rownames preserved", {
    f <- .mpc_fixture()
    res <- mashPosteriorContrast(f$pm, f$pv, f$orig)
    expect_equal(nrow(res), 5L)
    expect_identical(res$feature_id, rownames(f$pm))
    # 3 conditions -> 3 deviation + 3 pairwise = 6 contrasts x (mean/se/p)
    expect_equal(sum(grepl("^mean_contrast", names(res))), 6L)
    expect_equal(sum(grepl("^se_contrast", names(res))), 6L)
    expect_equal(sum(grepl("^p_contrast", names(res))), 6L)
    expect_true(
        any(grepl("deviation", names(res))) && any(grepl("_vs_", names(res)))
    )
    # column order: all mean_* precede se_*, which precede p_*
    contrastCols <- names(res)[grepl("_contrast_", names(res))]
    kinds <- sub("_contrast_.*", "", contrastCols)
    expect_false(is.unsorted(match(kinds, c("mean", "se", "p"))))
})

test_that("mashPosteriorContrast: features with < 2 tested conditions are dropped", {
    f <- .mpc_fixture()
    f$orig[2, ] <- 0 # feature 2 has no tested condition
    f$orig[4, c(2, 3)] <- 0 # feature 4 has only 1 tested condition
    res <- mashPosteriorContrast(f$pm, f$pv, f$orig)
    expect_equal(nrow(res), 3L)
    expect_false(any(res$feature_id %in% rownames(f$pm)[c(2, 4)]))
})

test_that("mashPosteriorContrast: grouping is forwarded to fitMashContrast", {
    f <- .mpc_fixture()
    res <- mashPosteriorContrast(
        f$pm,
        f$pv,
        f$orig,
        grouping = c(Ast = 1L, Mic = 1L, Oli = 0L)
    )
    expect_equal(nrow(res), 5L)
    expect_true(any(grepl("deviation", names(res))))
})

# ---------------------------------------------------------------------------
# Feature scores (calculateFeatureScores / nSignificantScore / scoreFromCs)
# ---------------------------------------------------------------------------

.fs_contrast <- function(nf = 8, conds = c("Ast", "Mic", "Oli"), seed = 1) {
    f <- .mpc_fixture(nf = nf, conds = conds, seed = seed)
    mashPosteriorContrast(f$pm, f$pv, f$orig)
}

test_that("calculateFeatureScores: one Z per condition from deviation contrasts", {
    cr <- .fs_contrast()
    fs <- calculateFeatureScores(cr, metaMethod = "REML")
    expect_setequal(fs$condition, c("Ast", "Mic", "Oli"))
    expect_equal(names(fs), c("condition", "zScore"))
    expect_true(all(is.finite(fs$zScore)))
    # empty when no deviation columns
    expect_equal(nrow(calculateFeatureScores(data.frame(x = 1))), 0L)
})

test_that("nSignificantScore: fraction of significant deviation contrasts in [0,1]", {
    cr <- .fs_contrast()
    ns <- nSignificantScore(cr, pCutoff = 0.5)
    expect_setequal(ns$condition, c("Ast", "Mic", "Oli"))
    expect_true(all(ns$ratio >= 0 & ns$ratio <= 1, na.rm = TRUE))
    # a hand-built p column: 2 of 4 below cutoff -> 0.5
    df <- data.frame(p_contrast_A_deviation = c(1e-8, 1e-9, 0.2, 0.3))
    expect_equal(nSignificantScore(df, pCutoff = 1e-5)$ratio, 0.5)
})

test_that("scoreFromCs: CS lead intersection score; NA on empty / no-overlap", {
    cr <- .fs_contrast()
    fm <- data.frame(
        variants = cr$feature_id,
        cs_order = c(1, 1, 1, 2, 2, 0, 0, 0),
        pip = c(0.6, 0.3, 0.1, 0.7, 0.3, 0.02, 0.02, 0.02),
        stringsAsFactors = FALSE
    )
    sc <- scoreFromCs(fm, cr, "Ast")
    expect_true(is.finite(sc) && sc >= 0)
    # no credible sets -> NA
    expect_true(is.na(scoreFromCs(
        data.frame(
            variants = character(0),
            cs_order = integer(0),
            pip = numeric(0)
        ),
        cr,
        "Ast"
    )))
    # CS leads that don't overlap the contrast variants -> NA
    fmNoOverlap <- data.frame(
        variants = c("chrX:1:A:G", "chrX:2:A:G"),
        cs_order = c(1, 1),
        pip = c(0.9, 0.1),
        stringsAsFactors = FALSE
    )
    expect_true(is.na(scoreFromCs(fmNoOverlap, cr, "Ast")))
})

# ---------------------------------------------------------------------------
# Coverage: feature-score edge branches + posterior-contrast empty result
# ---------------------------------------------------------------------------

test_that("calculateFeatureScores: NA when a condition's se column is absent", {
    cr <- data.frame(mean_contrast_Ast_deviation = c(1, 2))
    out <- calculateFeatureScores(cr)
    expect_equal(out$condition, "Ast")
    expect_true(is.na(out$zScore))
})

test_that("calculateFeatureScores: NA when no finite (effect, se) pair remains", {
    cr <- data.frame(
        mean_contrast_Ast_deviation = c(1, 2),
        se_contrast_Ast_deviation = c(0, -1)
    ) # se <= 0 -> all dropped
    expect_true(is.na(calculateFeatureScores(cr)$zScore))
})

test_that("nSignificantScore: empty result when there are no deviation p-value columns", {
    out <- nSignificantScore(data.frame(foo = 1:3, bar = 4:6))
    expect_equal(nrow(out), 0L)
    expect_named(out, c("condition", "ratio"))
})

test_that("scoreFromCs: falls back to the single pairwise contrast when no deviation column", {
    fm <- data.frame(
        cs_order = c(1, 1, 0),
        pip = c(0.9, 0.5, 0.1),
        variants = c("v1", "v2", "v3"),
        stringsAsFactors = FALSE
    )
    cr <- data.frame(
        p_contrast_Ast_vs_Mic = c(0.5, 0.2),
        mean_contrast_Ast_vs_Mic = c(1.0, 0.5),
        se_contrast_Ast_vs_Mic = c(0.2, 0.1),
        feature_id = c("v1", "v2"),
        stringsAsFactors = FALSE
    )
    # condition "Xyz" has no deviation column; exactly one pairwise contrast exists.
    expect_true(is.finite(scoreFromCs(fm, cr, condition = "Xyz")))
})

test_that("scoreFromCs: NA when neither a deviation nor a single pairwise column exists", {
    fm <- data.frame(
        cs_order = c(1, 0),
        pip = c(0.9, 0.1),
        variants = c("v1", "v2"),
        stringsAsFactors = FALSE
    )
    # Lead variant overlaps, but cr carries no p_contrast_*_deviation for the
    # condition and no pairwise p_contrast_*_vs_* column -> NA.
    cr <- data.frame(some_other_col = 1, feature_id = "v1")
    expect_true(is.na(scoreFromCs(fm, cr, condition = "Xyz")))
})

test_that("mashPosteriorContrast: empty frame when every feature is dropped", {
    f <- .mpc_fixture()
    f$orig[] <- 0 # no feature has any tested condition
    expect_equal(nrow(mashPosteriorContrast(f$pm, f$pv, f$orig)), 0L)
})

# ---------------------------------------------------------------------------
# Coverage: mash* helpers — NULL vhat default, SimpleList input, error guards
# ---------------------------------------------------------------------------

test_that("mashResidualCorrelation(mle): errors without a 'random' entry", {
    skip_if_not_installed("mashr")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    expect_error(
        mashResidualCorrelation(list(strong = ss), alpha = 0, method = "mle"),
        "requires a 'random' entry"
    )
})

test_that("mashResidualCorrelation accepts a SimpleList sumStatsList", {
    skip_if_not_installed("mashr")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    V <- suppressMessages(suppressWarnings(
        mashResidualCorrelation(
            S4Vectors::SimpleList(null = ss),
            alpha = 0,
            method = "simple"
        )
    ))
    expect_equal(dim(V), c(3L, 3L))
})

test_that("mashCovarianceComponents: SimpleList input + default (NULL) vhat", {
    skip_if_not_installed("mashr")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    cc <- suppressMessages(suppressWarnings(
        mashCovarianceComponents(
            S4Vectors::SimpleList(strong = ss),
            alpha = 0,
            components = "canonical",
            setSeed = 1L
        )
    ))
    expect_gt(length(cc), 0L)
})

test_that("mashPriorCovariances: SimpleList input (canonical only)", {
    skip_if_not_installed("mashr")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    pc <- suppressMessages(suppressWarnings(
        mashPriorCovariances(
            S4Vectors::SimpleList(strong = ss),
            alpha = 0,
            vhat = diag(3),
            components = "canonical"
        )
    ))
    expect_named(pc, c("U", "w", "loglik"))
})

test_that("mashModelFit: SimpleList input + default (NULL) vhat", {
    skip_if_not_installed("mashr")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    m <- suppressMessages(suppressWarnings(
        mashModelFit(
            S4Vectors::SimpleList(random = ss),
            alpha = 0,
            priorCovariances = list(identity = diag(3)),
            setSeed = 1L
        )
    ))
    expect_s3_class(m, "mash")
})

test_that("mashPosterior: default (NULL) vhat and excludeCondition dropping every condition", {
    skip_if_not_installed("mashr")
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    model <- .mashTestModel(ss)
    conds <- colnames(.mashSumStatsToMatrices(ss, "strong")$b)
    # NULL vhat path: omit vhat so mashPosterior fills the identity default.
    post <- suppressMessages(suppressWarnings(mashPosterior(
        model,
        ss,
        alpha = 0
    )))
    expect_equal(ncol(post$PosteriorMean), 3L)
    # excludeCondition removing every condition errors.
    expect_error(
        suppressMessages(suppressWarnings(
            mashPosterior(model, ss, alpha = 0, excludeCondition = conds)
        )),
        "drops every condition"
    )
})
