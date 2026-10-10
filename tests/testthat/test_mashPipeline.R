# =============================================================================
# Tests for utility helpers exported from R/mashPipeline.R that don't
# require mashr / flashier installs:
#   sanitizeMashData, makePairwiseContrastCol, sliceMashData,
#   metaAnalysisPerCondition

# Fixtures for the input-preparation tests, moved here with them.

# Helper: tiny QtlSumStats with 2 contexts on one trait. The mcols
# layout is controlled so each test can vary which columns are present.
.mssm_makeQtlSumStats <- function(
    mcolsBuilder,
    contexts = c("brain", "liver"),
    nSnp = 5L
) {
    set.seed(13L)
    gh <- new(
        "GenotypeHandle",
        path = "/tmp/sketch.gds",
        format = "gds",
        snpInfo = data.frame(
            SNP = sprintf("chr1:%d:A:G", 100L * seq_len(nSnp)),
            CHR = "1",
            BP = seq(100L, by = 100L, length.out = nSnp),
            A1 = "A",
            A2 = "G",
            stringsAsFactors = FALSE
        ),
        nSamples = 50L,
        sampleIds = paste0("s", seq_len(50L)),
        pgenPtr = NULL
    )
    ranges <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(
            start = seq(100L, by = 100L, length.out = nSnp),
            width = 1L
        )
    )
    entries <- map(seq_along(contexts), function(i) {
        gr <- ranges
        S4Vectors::mcols(gr) <- S4Vectors::DataFrame(mcolsBuilder(i, nSnp))
        gr
    })
    QtlSumStats(
        studyName = rep("s1", length(contexts)),
        context = contexts,
        trait = rep("g1", length(contexts)),
        entry = entries,
        genome = "hg19",
        ldSketch = gh,
        qcInfo = list(prebuilt = "synthetic")
    )
}

.qszm_gh <- function() {
    new(
        "GenotypeHandle",
        path = "/tmp/sketch.gds",
        format = "gds",
        snpInfo = data.frame(
            SNP = "chr1:100:A:G",
            CHR = "1",
            BP = 100L,
            A1 = "A",
            A2 = "G",
            stringsAsFactors = FALSE
        ),
        nSamples = 10L,
        sampleIds = paste0("s", seq_len(10L)),
        pgenPtr = NULL
    )
}

.qszmBeta <- function() {
    rn <- c("chr1:100:A:G", "chr1:200:C:T", "chr2:50:A:T")
    cn <- c("brain", "liver")
    list(
        bhat = matrix(
            c(0.5, -0.3, 0.2, 0.1, 0.4, -0.6),
            nrow = 3,
            dimnames = list(rn, cn)
        ),
        shat = matrix(
            c(0.1, 0.2, 0.15, 0.12, 0.09, 0.2),
            nrow = 3,
            dimnames = list(rn, cn)
        )
    )
}

# A multi-context QtlFineMappingResult fixture: two contexts sharing the same
# 6 variants, each with one credible set whose lead (max PIP) is a distinct
# variant, plus low-signal background for random/null sampling.
.mi_makeFmr <- function() {
    mkTL <- function(zvec, csIdx) {
        n <- length(zvec)
        data.frame(
            variant_id = paste0("chr1:", 100 * seq_len(n), ":A:G"),
            chrom = "1",
            pos = as.integer(100 * seq_len(n)),
            A1 = "G",
            A2 = "A",
            N = 1000,
            MAF = 0.2,
            marginal_beta = zvec * 0.05,
            marginal_se = 0.05,
            marginal_z = zvec,
            marginal_p = 2 * pnorm(-abs(zvec)),
            pip = {
                p <- rep(0.05, n)
                p[csIdx] <- seq(0.9, by = -0.2, length.out = length(csIdx))
                p
            },
            posterior_mean = zvec * 0.05,
            posterior_sd = 0.02,
            cs_95 = {
                cc <- rep("susie_0", n)
                cc[csIdx] <- "susie_1"
                cc
            },
            stringsAsFactors = FALSE
        )
    }
    vids <- paste0("chr1:", 100 * seq_len(6), ":A:G")
    e1 <- fineMappingRow(
        variantIds = vids,
        susieFit = list(x = 1),
        topLoci = mkTL(c(0.5, -1, 6.0, 0.2, 1.1, -0.3), c(3, 2))
    )
    e2 <- fineMappingRow(
        variantIds = vids,
        susieFit = list(x = 1),
        topLoci = mkTL(c(-0.4, 0.7, 0.1, 5.0, -1.2, 0.6), c(4, 5))
    )
    QtlFineMappingResult(
        studyName = c("s1", "s1"),
        context = c("brain", "blood"),
        trait = c("t1", "t1"),
        method = c("susie", "susie"),
        entry = list(e1, e2)
    )
}

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
    feat <- c("chr1:100:A:G", "chr1:200:A:G")
    cols <- "mean_contrast_brain_vs_blood"
    es <- matrix(c(0.1, 0.2), nrow = 2, ncol = 1, dimnames = list(feat, cols))
    se <- matrix(c(0.01, 0.02), nrow = 2, ncol = 1, dimnames = list(feat, cols))
    # seCutoff = 0.5 drops both rows
    out <- metaAnalysisPerCondition(es, se, seCutoff = 0.5)
    expect_true(all(is.na(out$meta_pvalue)))
    expect_true(all(is.na(out$meta_effect)))
})

test_that("metaAnalysisPerCondition runs DerSimonian-Laird when >=2 effects survive", {
    feat <- c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G")
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
    feat <- "chr1:100:A:G"
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
    expect_true(all(map_lgl(
        m2$fitted_g$Ulist,
        function(x) all(dim(x) == c(2L, 2L))
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
        dimnames = list("chr1:100:A:G", c("a", "b", "c"))
    )
    origMean[1, "b"] <- 1
    pm <- matrix(
        0,
        nrow = 1,
        ncol = 3,
        dimnames = list("chr1:100:A:G", c("a", "b", "c"))
    )
    pv <- array(diag(3), dim = c(3, 3, 1))
    dimnames(pv) <- list(c("a", "b", "c"), c("a", "b", "c"), NULL)
    expect_null(fitMashContrast(1L, origMean, pm, pv))
})

test_that("fitMashContrast: 2-tested-conditions fast path yields one pairwise contrast", {
    origMean <- matrix(
        c(0.5, 0.3, 0),
        nrow = 1,
        dimnames = list("chr1:100:A:G", c("a", "b", "c"))
    )
    pm <- matrix(
        c(0.5, 0.3, 0),
        nrow = 1,
        dimnames = list("chr1:100:A:G", c("a", "b", "c"))
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
        dimnames = list("chr1:100:A:G", c("a", "b", "c"))
    )
    pm <- matrix(
        c(0.5, 0.3, -0.2),
        nrow = 1,
        dimnames = list("chr1:100:A:G", c("a", "b", "c"))
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
    ss <- mashFixture()
    # Use the same fixture for strong/random; nPcs <= ncol - 1 (3 contexts).
    res <- suppressMessages(suppressWarnings(
        mashPipeline(
            sumStatsList = list(strong = ss, random = ss),
            alpha = 0,
            prior = MashPriorParam(nPcs = 2L),
            setSeed = 1L
        )
    ))
    expect_named(res, c("U", "w"))
    expect_type(res$U, "list")
    expect_gt(length(res$U), 0L)
    # Every covariance matrix is 3x3 (one row/col per context)
    expect_true(all(map_lgl(
        res$U,
        function(m) all(dim(m) == c(3L, 3L))
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
        dimnames = list("chr1:100:A:G", conds)
    )
    pm <- matrix(
        c(0.5, 0.3, -0.2, 0.4),
        nrow = 1,
        dimnames = list("chr1:100:A:G", conds)
    )
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
    es <- matrix(
        c(0.3, 0.5),
        nrow = 2,
        dimnames = list(c("chr1:100:A:G", "chr1:200:A:G"), cols)
    )
    se <- matrix(
        c(0.1, 0.1),
        nrow = 2,
        dimnames = list(c("chr1:100:A:G", "chr1:200:A:G"), cols)
    )
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
    ss <- mashFixture()
    vhat <- diag(3)
    # Empty list.
    expect_error(
        suppressMessages(suppressWarnings(
            mashPipeline(
                list(strong = ss),
                alpha = 0,
                residualCorrelation = vhat,
                prior = MashPriorParam(priorCovariances = list())
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
                prior = MashPriorParam(priorCovariances = list(diag(3)))
            )
        )),
        "non-empty named"
    )
})

test_that("mashPipeline rejects priorCovariances with wrong dimensions", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    ss <- mashFixture()
    vhat <- diag(3)
    expect_error(
        suppressMessages(suppressWarnings(
            mashPipeline(
                list(strong = ss),
                alpha = 0,
                residualCorrelation = vhat,
                prior = MashPriorParam(priorCovariances = list(myU = diag(2)))
            )
        )),
        "3 x 3 matrix"
    )
})

test_that("mashPipeline passes supplied residualCorrelation + priorCovariances through", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    ss <- mashFixture()
    vhat <- diag(3)
    U0 <- list(identity = diag(3), effectA = diag(c(1, 0, 0)))
    res <- suppressMessages(suppressWarnings(
        mashPipeline(
            list(strong = ss),
            alpha = 0,
            residualCorrelation = vhat,
            prior = MashPriorParam(priorCovariances = U0)
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
    ss <- mashFixture()
    # `residualCorrelationMethod = "simple"` runs
    # estimate_null_correlation_simple on the null set; leaving nPcs NULL
    # exercises the `nPcs <- ncol(Bhat) - 1` default in the cov_* chain.
    res <- suppressMessages(suppressWarnings(
        mashPipeline(
            list(strong = ss, random = ss, null = ss),
            alpha = 0,
            residualCorrelationMethod = "simple",
            setSeed = 1L
        )
    ))
    expect_named(res, c("U", "w"))
    expect_gt(length(res$U), 0L)
    expect_true(all(map_lgl(
        res$U,
        function(m) all(dim(m) == c(3L, 3L))
    )))
    expect_equal(sum(res$w), 1, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# residualCorrelationMethod
#
# Vhat used to be chosen by inspecting which partitions were present, so a
# caller without a null set silently assumed zero residual correlation. The
# choice is now always explicit.
# ---------------------------------------------------------------------------

test_that("the default is identity regardless of which partitions are given", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    ss <- mashFixture()
    prior <- mashTinyPrior()
    # The default must not depend on data shape -- that was the old silent
    # behaviour this argument replaced. `w` is what carries the evidence: it
    # is fitted against V, so an identical `w` across the three calls is the
    # same V reaching mash() each time. (`U` is the supplied prior passed
    # through, so comparing it would prove nothing here.)
    withNull <- suppressMessages(suppressWarnings(mashPipeline(
        list(strong = ss, random = ss, null = ss),
        alpha = 0,
        prior = MashPriorParam(priorCovariances = prior),
        setSeed = 1L
    )))
    withoutNull <- suppressMessages(suppressWarnings(mashPipeline(
        list(strong = ss, random = ss),
        alpha = 0,
        prior = MashPriorParam(priorCovariances = prior),
        setSeed = 1L
    )))
    explicit <- suppressMessages(suppressWarnings(mashPipeline(
        list(strong = ss, random = ss, null = ss),
        alpha = 0,
        residualCorrelationMethod = "identity",
        prior = MashPriorParam(priorCovariances = prior),
        setSeed = 1L
    )))
    expect_equal(withNull$w, withoutNull$w)
    expect_equal(withNull$w, explicit$w)
})

test_that("an unused null partition is reported", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    ss <- mashFixture()
    # Assembling a null set and still getting identity is more often an
    # oversight than an intent, so it must not pass silently.
    expect_message(
        suppressWarnings(mashPipeline(
            list(strong = ss, random = ss, null = ss),
            alpha = 0,
            prior = MashPriorParam(priorCovariances = mashTinyPrior()),
            setSeed = 1L
        )),
        "does not use it"
    )
    # Naming identity explicitly does not change that the set is unused.
    expect_message(
        suppressWarnings(mashPipeline(
            list(strong = ss, random = ss, null = ss),
            alpha = 0,
            residualCorrelationMethod = "identity",
            prior = MashPriorParam(priorCovariances = mashTinyPrior()),
            setSeed = 1L
        )),
        "'null' partition"
    )
})

test_that("no unused-null notice when there is nothing to ignore", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    ss <- mashFixture()
    noNull <- function() {
        suppressWarnings(mashPipeline(
            list(strong = ss, random = ss),
            alpha = 0,
            prior = MashPriorParam(priorCovariances = mashTinyPrior()),
            setSeed = 1L
        ))
    }
    expect_no_message(noNull(), message = "'null' partition")
    # Nor when the null set is actually being consumed.
    usesNull <- function() {
        suppressWarnings(mashPipeline(
            list(strong = ss, random = ss, null = ss),
            alpha = 0,
            residualCorrelationMethod = "simple",
            prior = MashPriorParam(priorCovariances = mashTinyPrior()),
            setSeed = 1L
        ))
    }
    expect_no_message(usesNull(), message = "does not use it")
})

test_that("mashPipeline honours a data-driven residualCorrelationMethod", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    ss <- mashFixture()
    sl <- list(strong = ss, random = ss, null = ss)
    prior <- mashTinyPrior()
    identityFit <- suppressMessages(suppressWarnings(mashPipeline(
        sl,
        alpha = 0,
        residualCorrelationMethod = "identity",
        prior = MashPriorParam(priorCovariances = prior),
        setSeed = 1L
    )))
    simpleFit <- suppressMessages(suppressWarnings(mashPipeline(
        sl,
        alpha = 0,
        residualCorrelationMethod = "simple",
        prior = MashPriorParam(priorCovariances = prior),
        setSeed = 1L
    )))
    # A different V has to move the fit, or the argument is not reaching it.
    expect_false(isTRUE(all.equal(identityFit$w, simpleFit$w)))
})

test_that("mashPipeline rejects an unknown residualCorrelationMethod", {
    ss <- mashFixture()
    expect_error(
        mashPipeline(
            list(strong = ss, random = ss),
            alpha = 0,
            residualCorrelationMethod = "bogus"
        ),
        "bogus"
    )
})

test_that("a method needing a partition it lacks errors through mashPipeline", {
    skip_if_not_installed("mashr")
    ss <- mashFixture()
    expect_error(
        mashPipeline(
            list(strong = ss, random = ss),
            alpha = 0,
            residualCorrelationMethod = "simple"
        ),
        "requires a 'null' entry"
    )
})

test_that("a supplied residualCorrelation wins over the method", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    ss <- mashFixture()
    # 'simple' would error without a null set; the supplied matrix means the
    # estimator is never consulted.
    res <- suppressMessages(suppressWarnings(mashPipeline(
        list(strong = ss, random = ss),
        alpha = 0,
        residualCorrelation = diag(3),
        residualCorrelationMethod = "simple",
        prior = MashPriorParam(priorCovariances = mashTinyPrior()),
        setSeed = 1L
    )))
    expect_named(res, c("U", "w"))
})

# ---------------------------------------------------------------------------
# mashResidualCorrelation — the Vhat estimator extracted from mashPipeline.
# ---------------------------------------------------------------------------

test_that("mashResidualCorrelation(identity) is an identity of the right size", {
    skip_if_not_installed("mashr")
    ss <- mashFixture()
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
    ss <- mashFixture()
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
    ss <- mashFixture()
    expect_error(
        mashResidualCorrelation(
            list(strong = ss),
            alpha = 0,
            method = "simple"
        ),
        "requires a 'null' entry"
    )
})

test_that("mashResidualCorrelation(simpleSpecific) returns a null correlation", {
    skip_if_not_installed("mashr")
    ss <- mashFixture()
    v <- suppressMessages(suppressWarnings(
        mashResidualCorrelation(
            list(strong = ss, null = ss),
            alpha = 0,
            method = "simpleSpecific"
        )
    ))
    expect_equal(dim(v), c(3L, 3L))
    expect_equal(unname(diag(v)), rep(1, 3), tolerance = 1e-6)
})

test_that("mashResidualCorrelation(corshrink) returns a 3x3 correlation matrix", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("CorShrink")
    ss <- mashFixture()
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
    ss <- mashFixture()
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
    ss <- mashFixture()
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
# flash + flashNonneg) refined by cov_ed.
# ---------------------------------------------------------------------------

test_that("mashPriorCovariances computes the default (all-but-udr) prior", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    ss <- mashFixture()
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
    expect_true(all(map_lgl(
        pc$U,
        function(m) all(dim(m) == c(3L, 3L))
    )))
    expect_equal(sum(pc$w), 1, tolerance = 1e-6)
    expect_null(pc$loglik)
})

test_that("mashPriorCovariances passes a supplied prior through unchanged", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    ss <- mashFixture()
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
    ss <- mashFixture()
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

test_that("mashPriorCovariances(flashNonneg) adds components vs flash-only", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    ss <- mashFixture()
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
            components = c("canonical", "pca", "flash", "flashNonneg"),
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
    ss <- mashFixture()
    # A toy-sized udr config: `n_unconstrained` dominates the cost, and the default
    # (50, sized for many-condition data) is pathological on a 3-condition fixture
    # (it drove a 5+ minute fit). 2 unconstrained matrices suffice to exercise the
    # udr path.
    pc <- suppressMessages(suppressWarnings(
        mashPriorCovariances(
            list(strong = ss),
            alpha = 0,
            vhat = diag(3),
            engine = CovUdrOptions(
                init = UdInitOptions(n_unconstrained = 2L),
                fit = UdFitOptions(maxiter = 20L)
            ),
            setSeed = 1L
        )
    ))
    expect_gt(length(pc$U), 0L)
    expect_true(all(map_lgl(
        pc$U,
        function(m) all(dim(m) == c(3L, 3L))
    )))
    expect_equal(sum(pc$w), 1, tolerance = 1e-6)
})

test_that("mashPriorCovariances TED update errors clearly on non-i.i.d. data", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    skip_if_not_installed("udr")
    ss <- mashFixture()
    expect_error(
        suppressMessages(suppressWarnings(
            mashPriorCovariances(
                list(strong = ss),
                alpha = 0,
                vhat = diag(3),
                engine = CovUdrOptions(
                    fit = UdFitOptions(unconstrained.update = "ted")
                ),
                setSeed = 1L
            )
        )),
        "i.i.d"
    )
})

test_that("mashPriorCovariances rejects an unknown component", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    ss <- mashFixture()
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
    ss <- mashFixture()
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
    expect_true(all(map_lgl(
        fl,
        function(m) all(dim(m) == c(3L, 3L))
    )))
})

test_that("mashCovarianceComponents default builds all non-udr components", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    ss <- mashFixture()
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
    ss <- mashFixture()
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
    # Canonical components are structural hypotheses and reach the prior
    # unrefined, so their names survive verbatim. The data-driven ones are
    # replaced by their Extreme Deconvolution refinements (mashr names those
    # "ED_<source>"), rather than appearing twice as they used to.
    # Derived rather than pattern-matched: mashr names singleton components
    # after the conditions themselves (brain, blood, ...), so only a
    # canonical-only build says reliably which names are canonical.
    canonicalNames <- names(suppressMessages(suppressWarnings(
        mashCovarianceComponents(
            list(strong = ss),
            alpha = 0,
            vhat = diag(3),
            components = "canonical",
            setSeed = 1L
        )
    )))
    expect_true(length(canonicalNames) > 0L)
    expect_true(all(canonicalNames %in% names(prior$U)))
    dataDrivenNames <- setdiff(names(comps), canonicalNames)
    expect_true(length(dataDrivenNames) > 0L)
    expect_false(any(dataDrivenNames %in% names(prior$U)))
    expect_true(any(grepl("^ED", names(prior$U))))
})

test_that("mashCovarianceComponents rejects unknown components", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    ss <- mashFixture()
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
    ss <- mashFixture()
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
            engine = "covEd",
            setSeed = 1L
        )
    ))
    # Supplied components are treated as the engine's input, so what comes
    # back are their refinements, not the originals.
    expect_false(any(names(comps) %in% names(pr$U)))
    expect_true(any(grepl("^ED", names(pr$U))))
    expect_equal(length(pr$U), length(comps))
    expect_equal(sum(pr$w), 1, tolerance = 1e-6)
})

test_that("mashPriorCovariances validates priorComponents", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    ss <- mashFixture()
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
    ss <- mashFixture()
    full <- suppressMessages(suppressWarnings(
        mashPipeline(
            list(strong = ss, random = ss, null = ss),
            alpha = 0,
            residualCorrelationMethod = "simple",
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

# ---------------------------------------------------------------------------
# Prior-shape normalisation (.mashAsUlist)
#
# The two producers disagree on shape: mashCovarianceComponents() returns a
# bare Ulist, mashPriorCovariances() wraps it as list(U, w, loglik). Every
# consumer must take either, so a chain can be written without reaching for
# `$U`.
# ---------------------------------------------------------------------------

test_that(".mashAsUlist unwraps the wrapper and passes a bare Ulist through", {
    ulist <- list(identity = diag(3), effectA = diag(c(1, 0, 0)))
    wrapped <- list(U = ulist, w = c(0.5, 0.5), loglik = -10)

    expect_equal(pecotmr:::.mashAsUlist(wrapped), ulist)
    expect_equal(pecotmr:::.mashAsUlist(ulist), ulist)
    expect_null(pecotmr:::.mashAsUlist(NULL))
})

test_that(".mashAsUlist does not mistake a covariance named U for the wrapper", {
    # A genuine Ulist entry named "U" is a matrix, not a list, which is what
    # keeps the two shapes distinguishable.
    ulist <- list(U = diag(3), identity = diag(3))
    expect_equal(pecotmr:::.mashAsUlist(ulist), ulist)
})

test_that("mashPriorCovariances accepts either shape for its two prior args", {
    skip_if_not_installed("mashr")
    ss <- mashFixture()
    sl <- list(strong = ss, random = ss)
    ulist <- list(identity = diag(3), effectA = diag(c(1, 0, 0)))
    wrapped <- list(U = ulist, w = NULL, loglik = NULL)

    bare <- suppressMessages(suppressWarnings(mashPriorCovariances(
        sl,
        alpha = 0,
        vhat = diag(3),
        priorCovariances = ulist
    )))
    wrap <- suppressMessages(suppressWarnings(mashPriorCovariances(
        sl,
        alpha = 0,
        vhat = diag(3),
        priorCovariances = wrapped
    )))
    expect_equal(names(wrap$U), names(bare$U))
    expect_equal(wrap$w, bare$w)
})

test_that("fitMashContrast consumes a mashPosterior result", {
    skip_if_not_installed("mashr")
    ss <- mashFixture()
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
        variants = c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G"),
        stringsAsFactors = FALSE
    )
    cr <- data.frame(
        p_contrast_Ast_vs_Mic = c(0.5, 0.2),
        mean_contrast_Ast_vs_Mic = c(1.0, 0.5),
        se_contrast_Ast_vs_Mic = c(0.2, 0.1),
        feature_id = c("chr1:100:A:G", "chr1:200:A:G"),
        stringsAsFactors = FALSE
    )
    # condition "Xyz" has no deviation column; exactly one pairwise contrast exists.
    expect_true(is.finite(scoreFromCs(fm, cr, condition = "Xyz")))
})

test_that("scoreFromCs: NA when neither a deviation nor a single pairwise column exists", {
    fm <- data.frame(
        cs_order = c(1, 0),
        pip = c(0.9, 0.1),
        variants = c("chr1:100:A:G", "chr1:200:A:G"),
        stringsAsFactors = FALSE
    )
    # Lead variant overlaps, but cr carries no p_contrast_*_deviation for the
    # condition and no pairwise p_contrast_*_vs_* column -> NA.
    cr <- data.frame(some_other_col = 1, feature_id = "chr1:100:A:G")
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
    ss <- mashFixture()
    expect_error(
        mashResidualCorrelation(list(strong = ss), alpha = 0, method = "mle"),
        "requires a 'random' entry"
    )
})

test_that("mashResidualCorrelation accepts a SimpleList sumStatsList", {
    skip_if_not_installed("mashr")
    ss <- mashFixture()
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
    ss <- mashFixture()
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
    ss <- mashFixture()
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

test_that("contrast rows fall back to the positional index when unnamed", {
    # posteriorMean usually carries feature rownames, but an unnamed matrix
    # still needs a feature_id -- the position, not NA, so the contrast table
    # stays joinable.
    f <- pecotmr:::.mashContrastDf
    pm <- matrix(1:4, 2L, 2L)
    d <- f(
        1L,
        pm,
        "c1",
        matrix(0.1, 2L, 1L),
        matrix(0.2, 2L, 1L),
        matrix(0.3, 2L, 1L)
    )
    expect_equal(d$feature_id, "1")
    rownames(pm) <- c("f1", "f2")
    d2 <- f(
        2L,
        pm,
        "c1",
        matrix(0.1, 2L, 1L),
        matrix(0.2, 2L, 1L),
        matrix(0.3, 2L, 1L)
    )
    expect_equal(d2$feature_id, "f2")
})

test_that("mashPipeline helpers: argument guards fire", {
    expect_error(
        fitMashContrast(
            0L,
            matrix(0, 2, 2),
            matrix(0, 2, 2),
            array(0, c(2, 2, 2))
        ),
        "index.*Must be >= 1"
    )
    expect_error(
        updateMashModelCov(list(), allSamples = 1L, samples = "a"),
        "allSamples.*Must be of type 'character'"
    )
    expect_error(
        sliceMashData("not-a-list", vhat = diag(2), snps = 1L, samples = NULL),
        "data.*Must be of type 'list'"
    )
    expect_error(
        calculateFeatureScores(data.frame(), metaMethod = 1L),
        "metaMethod.*Must be of type 'string'"
    )
    expect_error(
        nSignificantScore(data.frame(), pCutoff = 2),
        "pCutoff.*is not <= 1"
    )
    expect_error(
        makePairwiseContrastCol(c("a", "b", "c"), template = c(a = 0)),
        "pair.*Must have length 2"
    )
    expect_error(sanitizeMashData("nope"), "data.*Must be of type 'list'")
})

test_that("MashComponentParam takes plain lists or constructors per entry", {
    skip_if_not_installed("mashr")
    ca <- MashComponentParam(
        pca = list(subset = 1:5),
        canonical = CovCanonicalOptions()
    )
    expect_s4_class(ca, "MashComponentParam")
    expect_s4_class(ca, "MethodParam")
    # Naming a component selects it. The entries live in one named-list
    # slot rather than one nullable slot per component: a record selecting
    # pca used to show three NULL slots, which told a reader nothing.
    expect_setequal(names(ca$components), c("canonical", "pca"))
    expect_null(ca$components$flash)
    expect_setequal(names(ca), "components")
    # Each entry is still an engine argument bag, whichever form it arrived
    # in: a plain list is spliced into that component's constructor.
    expect_s4_class(ca$components$pca, "MethodOptions")
    expect_equal(ca$components$pca$subset, 1:5)
    expect_error(MashComponentParam(pca = list(nope = 1)), "unknown argument")
    expect_error(
        MashComponentParam(pca = CovCanonicalOptions()),
        "was built with the constructor for 'canonical'"
    )
    # The four components are formals now, so a misspelled one is R's own
    # error rather than a hand-rolled "unknown method" check.
    expect_error(MashComponentParam(nosuch = list()), "unused argument")
})

test_that("MashPriorParam takes only a MashComponentParam for components", {
    skip_if_not_installed("mashr")
    # The slot was a character|MethodOptions union while the aggregator was
    # still a MethodOptions, so any engine bag was accepted there. Now it is
    # exact: a component record, or the component names.
    expect_s4_class(
        MashPriorParam(components = MashComponentParam(pca = list())),
        "MashPriorParam"
    )
    expect_s4_class(
        MashPriorParam(components = c("canonical", "pca")),
        "MashPriorParam"
    )
    expect_error(
        MashPriorParam(components = CovPcaOptions()),
        "must be component names or a MashComponentParam"
    )
})

test_that("mashPriorCovariances takes components as names or a record", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    ss <- mashFixture()
    byRecord <- suppressMessages(suppressWarnings(
        mashPriorCovariances(
            list(strong = ss),
            alpha = 0,
            vhat = diag(3),
            components = MashComponentParam(canonical = CovCanonicalOptions()),
            engine = "none",
            setSeed = 1L
        )
    ))
    # Naming a component in the record selects it, so canonical-only here.
    expect_true(length(byRecord$U) > 0L)
    expect_false(any(grepl("^ED", names(byRecord$U))))
})

test_that("engine 'none' leaves the data-driven components unrefined", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("flashier")
    ss <- mashFixture()
    raw <- suppressMessages(suppressWarnings(
        mashPriorCovariances(
            list(strong = ss),
            alpha = 0,
            vhat = diag(3),
            components = c("pca"),
            engine = "none",
            nPcs = 2L,
            setSeed = 1L
        )
    ))
    # mashr supports using generator output directly as the prior; pecotmr
    # used to refine unconditionally.
    expect_true(any(grepl("^PCA", names(raw$U))))
    expect_false(any(grepl("^ED", names(raw$U))))
})

test_that("mashPriorCovariances rejects an unknown engine", {
    skip_if_not_installed("mashr")
    ss <- mashFixture()
    expect_error(
        mashPriorCovariances(
            list(strong = ss),
            alpha = 0,
            vhat = diag(3),
            engine = "ud_ted"
        ),
        "`engine` must be one of"
    )
})

test_that("mashPipeline forwards mashDataArgs and mashArgs downstream", {
    skip_if_not_installed("mashr")
    dataSeen <- list()
    mashSeen <- list()
    realData <- mashr::mash_set_data
    realMash <- mashr::mash
    ss <- mashFixture(20L)
    suppressMessages(suppressWarnings(with_mocked_bindings(
        mashPipeline(
            list(strong = ss, random = ss, null = ss),
            alpha = 0,
            residualCorrelationMethod = "simple",
            prior = MashPriorParam(priorCovariances = mashTinyPrior()),
            mashDataArgs = MashDataOptions(zero_Shat_reset = 0.25),
            mashArgs = MashOptions(nullweight = 3)
        ),
        mash_set_data = function(...) {
            dataSeen[[length(dataSeen) + 1L]] <<- list(...)
            realData(...)
        },
        mash = function(...) {
            mashSeen[[length(mashSeen) + 1L]] <<- list(...)
            realMash(...)
        },
        .package = "mashr"
    )))
    # Both the Vhat estimator and the prior/weight fit build mash data, and
    # every one of those calls must carry the caller's settings.
    expect_gte(length(dataSeen), 2L)
    expect_true(all(map_dbl(dataSeen, "zero_Shat_reset") == 0.25))
    expect_gte(length(mashSeen), 1L)
    expect_true(all(map_dbl(mashSeen, "nullweight") == 3))
})

test_that("mashResidualCorrelation takes method as a name or a constructor", {
    skip_if_not_installed("mashr")
    seen <- NULL
    real <- mashr::estimate_null_correlation_simple
    ss <- mashFixture(20L)
    byName <- suppressMessages(suppressWarnings(mashResidualCorrelation(
        list(strong = ss, null = ss),
        alpha = 0,
        method = "simple"
    )))
    byCtor <- suppressMessages(suppressWarnings(with_mocked_bindings(
        mashResidualCorrelation(
            list(strong = ss, null = ss),
            alpha = 0,
            method = MashCorSimpleOptions(z_thresh = 1)
        ),
        estimate_null_correlation_simple = function(...) {
            seen <<- list(...)
            real(...)
        },
        .package = "mashr"
    )))
    expect_equal(seen$z_thresh, 1)
    expect_equal(dim(byName), dim(byCtor))
    expect_error(
        mashResidualCorrelation(
            list(strong = ss),
            alpha = 0,
            method = MashOptions()
        ),
        "unknown engine 'mash'"
    )
})

test_that("mashResidualCorrelation routes CorShrinkOptions to CorShrink", {
    skip_if_not_installed("mashr")
    skip_if_not_installed("CorShrink")
    seen <- NULL
    real <- CorShrink::CorShrinkData
    ss <- mashFixture(20L)
    suppressMessages(suppressWarnings(with_mocked_bindings(
        mashResidualCorrelation(
            list(strong = ss, null = ss),
            alpha = 0,
            method = CorShrinkOptions(cor_method = "pearson")
        ),
        CorShrinkData = function(...) {
            seen <<- list(...)
            real(...)
        },
        .package = "CorShrink"
    )))
    expect_equal(seen$cor_method, "pearson")
    # The pecotmr defaults ride along.
    expect_equal(seen$image, "null")
    expect_equal(seen$ash.control, list(mixcompdist = "halfuniform"))
})

test_that(".mashUdControl defers to udr only where udr chooses well", {
    skip_if_not_installed("udr")
    # udr resolves its NA default as ifelse(is.matrix(fit$V), "ted", "none").
    # On the z scale that is TED, which is what we want, so the field is
    # omitted and udr decides.
    iid <- pecotmr:::.mashUdControl(UdFitOptions(), 3L, iid = TRUE)
    expect_false(is_in("unconstrained.update", names(iid)))
    expect_equal(
        udr::ud_fit_control_default()$unconstrained.update,
        NA
    )
    # On the beta scale ud_init() leaves a per-variant V, so udr would pick
    # "none" and never refine the unconstrained components. Pin ED, which
    # does not need i.i.d. data.
    beta <- pecotmr:::.mashUdControl(UdFitOptions(), 3L, iid = FALSE)
    expect_equal(beta$unconstrained.update, "ed")
    # Either way the caller's own value wins.
    expect_equal(
        pecotmr:::.mashUdControl(
            UdFitOptions(unconstrained.update = "ed"),
            3L,
            iid = TRUE
        )$unconstrained.update,
        "ed"
    )
    expect_equal(
        pecotmr:::.mashUdControl(
            UdFitOptions(unconstrained.update = "ted"),
            3L,
            iid = FALSE
        )$unconstrained.update,
        "ted"
    )
    # The settings pecotmr owns outright are unaffected by the scale.
    expect_equal(iid$scaled.update, "fa")
    expect_equal(iid$lambda, 3L)
})

test_that("the mash region partitions forward every cfg field", {
    # mashInput() builds the record inline before handing it on.
    fields <- .cfgFieldNames("mashInput")
    accepted <- setdiff(
        names(formals(get(
            ".mashObjectPartitions",
            envir = asNamespace("pecotmr")
        ))),
        c("...", "obj")
    )
    passed <- .argNamesPassed(".mashRegionPartitions")
    expect_equal(setdiff(intersect(fields, accepted), passed), character(0))
})

test_that("filterInvalidSummaryStat replaces NaN/Inf in bhat", {
    dat <- list(
        bhat = data.frame(a = c(1, NaN, 3), b = c(Inf, 2, -Inf)),
        sbhat = data.frame(a = c(0.1, 0.2, 0.3), b = c(0.1, NA, 0.3))
    )
    result <- filterInvalidSummaryStat(dat, bhat = "bhat", sbhat = "sbhat")
    expect_true(all(!is.nan(result$bhat)))
    expect_true(all(!is.infinite(result$bhat)))
    # NaN/Inf in bhat replaced with 0
    expect_equal(unname(result$bhat[1, 2]), 0) # Inf -> 0
})

test_that("filterInvalidSummaryStat replaces NaN/Inf in sbhat", {
    dat <- list(
        bhat = data.frame(a = c(1, 2, 3)),
        sbhat = data.frame(a = c(0.1, NaN, Inf))
    )
    result <- filterInvalidSummaryStat(dat, bhat = "bhat", sbhat = "sbhat")
    # NaN/Inf in sbhat replaced with 1000
    expect_equal(unname(result$sbhat[1, "a"]), 0.1)
    expect_equal(unname(result$sbhat[2, "a"]), 1000)
    expect_equal(unname(result$sbhat[3, "a"]), 1000)
})

test_that("filterInvalidSummaryStat filters by missing_rate when null.b present", {
    dat <- list(
        bhat = data.frame(a = c(0, 0, 1, 2), b = c(0, 0, 0, 3)),
        sbhat = data.frame(a = c(1, 1, 1, 1), b = c(1, 1, 1, 1)),
        null.b = TRUE
    )
    result <- filterInvalidSummaryStat(
        dat,
        bhat = "bhat",
        sbhat = "sbhat",
        filterByMissingRate = 0.5
    )
    expect_equal(nrow(result$bhat), 2) # rows 3 and 4 survive
})

test_that("filterInvalidSummaryStat filters by missing_rate when random.b present", {
    dat <- list(
        bhat = data.frame(a = c(0, 1, 2), b = c(0, 1, 3)),
        sbhat = data.frame(a = c(1, 1, 1), b = c(1, 1, 1)),
        random.b = TRUE
    )
    result <- filterInvalidSummaryStat(
        dat,
        bhat = "bhat",
        sbhat = "sbhat",
        filterByMissingRate = 0.5
    )
    expect_true(nrow(result$bhat) < 3)
})

test_that("filterInvalidSummaryStat btoz with .b and .s pattern creates condition.z", {
    dat <- list(
        strong.b = data.frame(a = c(1, 2), b = c(3, 4)),
        strong.s = data.frame(a = c(0.1, 0.2), b = c(0.3, 0.4))
    )
    result <- filterInvalidSummaryStat(
        dat,
        bhat = "strong.b",
        sbhat = "strong.s",
        btoz = TRUE,
        sigPCutoff = NULL
    )
    expect_true("strong.z" %in% names(result))
    expect_true(is.matrix(result$strong.z))
    expect_equal(nrow(result$strong.z), 2)
    expect_equal(ncol(result$strong.z), 2)
    expect_equal(
        as.numeric(result$strong.z),
        c(10, 10, 10, 10),
        tolerance = 1e-10
    )
})

test_that("filterInvalidSummaryStat btoz when bhat/sbhat data is NULL creates NULL z", {
    dat <- list(
        strong.b = NULL,
        strong.s = data.frame(a = c(0.1))
    )
    result <- filterInvalidSummaryStat(
        dat,
        bhat = "strong.b",
        sbhat = "strong.s",
        btoz = TRUE,
        sigPCutoff = NULL
    )
    expect_true("strong.z" %in% names(result))
    expect_null(result$strong.z)
})

test_that("filterInvalidSummaryStat btoz without .b/.s pattern creates generic z", {
    dat <- list(
        bhat = data.frame(a = c(1, 2, 3)),
        sbhat = data.frame(a = c(0.5, 1, 0.5))
    )
    result <- filterInvalidSummaryStat(
        dat,
        bhat = "bhat",
        sbhat = "sbhat",
        btoz = TRUE
    )
    expect_true("z" %in% names(result))
    expect_equal(as.numeric(result$z[, 1]), c(2, 2, 6))
})

test_that("filterInvalidSummaryStat btoz creates NULL z when bhat is NULL (no .b suffix)", {
    dat <- list(
        bhat = NULL,
        sbhat = data.frame(a = c(0.1))
    )
    result <- filterInvalidSummaryStat(
        dat,
        bhat = "bhat",
        sbhat = "sbhat",
        btoz = TRUE
    )
    expect_true("z" %in% names(result))
    expect_null(result$z)
})

test_that("filterInvalidSummaryStat btoz filters strong.z by significance cutoff", {
    dat <- list(
        strong.b = data.frame(
            a = c(1, 0.01, 0.02, 2),
            b = c(0.01, 0.01, 0.01, 0.01)
        ),
        strong.s = data.frame(
            a = c(0.1, 0.1, 0.1, 0.1),
            b = c(0.1, 0.1, 0.1, 0.1)
        ),
        strong.z = NULL
    )
    result <- filterInvalidSummaryStat(
        dat,
        bhat = "strong.b",
        sbhat = "strong.s",
        btoz = TRUE,
        sigPCutoff = 1E-6
    )
    expect_true("strong.z" %in% names(result))
    expect_equal(nrow(result$strong.z), 2)
    expect_equal(nrow(result$strong.b), 2)
    expect_equal(nrow(result$strong.s), 2)
})

test_that("filterInvalidSummaryStat processes z directly with null component", {
    dat <- list(
        strong = list(z = data.frame(a = c(5, NaN, 0.1), b = c(1, 2, Inf))),
        random = list(z = data.frame(a = c(0.5, 0.2), b = c(0.3, 0.4))),
        null = list(z = data.frame(a = c(0.01, NaN), b = c(Inf, 0.02)))
    )
    result <- filterInvalidSummaryStat(dat, z = "z")
    expect_true(all(!is.nan(result$strong$z)))
    expect_true(all(!is.infinite(result$strong$z)))
    expect_true(all(!is.nan(result$null$z)))
    expect_true(all(!is.infinite(result$null$z)))
})

test_that("filterInvalidSummaryStat z path applies significance cutoff to strong.z", {
    dat <- list(
        strong = list(
            z = data.frame(a = c(10, 0.1, 0.2), b = c(0.1, 0.1, 0.1))
        ),
        random = list(
            z = data.frame(a = c(0.5, 0.2, 0.3), b = c(0.3, 0.4, 0.2))
        )
    )
    result <- filterInvalidSummaryStat(dat, z = "z", sigPCutoff = 1E-6)
    expect_equal(nrow(result$strong$z), 1)
})

test_that("filterInvalidSummaryStat z path with filterByMissingRate", {
    dat <- list(
        random = list(z = data.frame(a = c(0, 0, 5), b = c(0, 3, 4)))
    )
    result <- filterInvalidSummaryStat(dat, z = "z", filterByMissingRate = 0.5)
    expect_equal(nrow(result$random$z), 2)
})

test_that("filterInvalidSummaryStat processes bhat/sbhat without filterByMissingRate when no null.b/random.b", {
    dat <- list(
        bhat = data.frame(a = c(1, NaN, 3), b = c(Inf, 2, -Inf)),
        sbhat = data.frame(a = c(0.1, NA, 0.3), b = c(0.1, 0.2, NaN))
    )
    result <- filterInvalidSummaryStat(
        dat,
        bhat = "bhat",
        sbhat = "sbhat",
        filterByMissingRate = 0.5
    )
    expect_equal(nrow(result$bhat), 3)
    expect_equal(unname(result$bhat[2, "a"]), 0)
    expect_equal(unname(result$sbhat[3, "b"]), 1000)
})

test_that("filterInvalidSummaryStat with filterByMissingRate=NULL keeps all rows even with null.b", {
    dat <- list(
        bhat = data.frame(a = c(0, 0, 1), b = c(0, 0, 1)),
        sbhat = data.frame(a = c(1, 1, 1), b = c(1, 1, 1)),
        null.b = TRUE
    )
    result <- filterInvalidSummaryStat(
        dat,
        bhat = "bhat",
        sbhat = "sbhat",
        filterByMissingRate = NULL
    )
    expect_equal(nrow(result$bhat), 3)
})

test_that("filterInvalidSummaryStat z path handles NULL strong component", {
    dat <- list(
        strong = NULL,
        random = list(z = data.frame(a = c(0.5, 0.2), b = c(0.3, 0.4)))
    )
    result <- filterInvalidSummaryStat(dat, z = "z")
    expect_null(result$strong)
    expect_true(!is.null(result$random$z))
})

test_that("filterMixtureComponents removes zero matrices", {
    U <- list(
        comp1 = matrix(
            c(1, 0, 0, 2),
            2,
            2,
            dimnames = list(c("A", "B"), c("A", "B"))
        ),
        comp2 = matrix(
            c(0, 0, 0, 0),
            2,
            2,
            dimnames = list(c("A", "B"), c("A", "B"))
        ),
        A = matrix(
            c(3, 0, 0, 4),
            2,
            2,
            dimnames = list(c("A", "B"), c("A", "B"))
        )
    )
    w <- c(comp1 = 0.5, comp2 = 0.3, A = 0.2)
    result <- filterMixtureComponents(c("A", "B"), U, w)
    expect_false("comp2" %in% names(result$U))
})

test_that("filterMixtureComponents removes matrices below weight cutoff", {
    U <- list(
        comp1 = matrix(
            c(1, 0, 0, 2),
            2,
            2,
            dimnames = list(c("A", "B"), c("A", "B"))
        ),
        comp2 = matrix(
            c(0.1, 0, 0, 0.1),
            2,
            2,
            dimnames = list(c("A", "B"), c("A", "B"))
        )
    )
    w <- c(comp1 = 0.9, comp2 = 0.00001)
    result <- filterMixtureComponents(c("A", "B"), U, w, wCutoff = 1e-4)
    expect_false("comp2" %in% names(result$U))
})

test_that("filterMixtureComponents errors on missing conditions", {
    U <- list(
        comp1 = matrix(
            c(1, 0, 0, 2),
            2,
            2,
            dimnames = list(c("A", "B"), c("A", "B"))
        )
    )
    expect_error(
        filterMixtureComponents(c("A", "C"), U),
        "not found in matrix"
    )
})

test_that("filterMixtureComponents removes components named as filtered conditions", {
    U <- list(
        comp1 = matrix(
            c(1, 0, 0, 0, 2, 0, 0, 0, 3),
            3,
            3,
            dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
        ),
        C = matrix(
            c(4, 0, 0, 0, 5, 0, 0, 0, 6),
            3,
            3,
            dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
        )
    )
    w <- c(comp1 = 0.6, C = 0.4)
    result <- filterMixtureComponents(c("A", "B"), U, w)
    expect_false("C" %in% names(result$U))
    expect_equal(nrow(result$U$comp1), 2)
    expect_equal(ncol(result$U$comp1), 2)
})

test_that("filterMixtureComponents renormalizes weights to preserve original sum", {
    U <- list(
        comp1 = matrix(
            c(1, 0, 0, 2),
            2,
            2,
            dimnames = list(c("A", "B"), c("A", "B"))
        ),
        comp2 = matrix(
            c(3, 0, 0, 4),
            2,
            2,
            dimnames = list(c("A", "B"), c("A", "B"))
        ),
        comp3 = matrix(
            c(0, 0, 0, 0),
            2,
            2,
            dimnames = list(c("A", "B"), c("A", "B"))
        )
    )
    w <- c(comp1 = 0.5, comp2 = 0.3, comp3 = 0.2)
    original_sum <- sum(w)

    result <- filterMixtureComponents(c("A", "B"), U, w)
    expect_false("comp3" %in% names(result$U))
    expect_equal(sum(result$w), original_sum, tolerance = 1e-10)
    expect_true(result$w["comp1"] > 0.5)
    expect_true(result$w["comp2"] > 0.3)
})

test_that("filterMixtureComponents subsets 3x3 matrices to 2x2 and removes filtered condition names", {
    U <- list(
        comp1 = matrix(
            c(1, 0.1, 0, 0.1, 2, 0, 0, 0, 3),
            3,
            3,
            dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
        ),
        B = matrix(
            c(4, 0, 0, 0, 5, 0, 0, 0, 6),
            3,
            3,
            dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
        )
    )
    w <- c(comp1 = 0.7, B = 0.3)

    result <- filterMixtureComponents(c("A", "C"), U, w)
    expect_false("B" %in% names(result$U))
    expect_equal(nrow(result$U$comp1), 2)
    expect_equal(ncol(result$U$comp1), 2)
    expect_equal(rownames(result$U$comp1), c("A", "C"))
})

test_that("filterMixtureComponents handles NULL weights gracefully", {
    U <- list(
        comp1 = matrix(
            c(1, 0, 0, 2),
            2,
            2,
            dimnames = list(c("A", "B"), c("A", "B"))
        ),
        comp2 = matrix(
            c(0, 0, 0, 0),
            2,
            2,
            dimnames = list(c("A", "B"), c("A", "B"))
        )
    )
    result <- filterMixtureComponents(c("A", "B"), U, w = NULL)
    expect_false("comp2" %in% names(result$U))
    expect_true("comp1" %in% names(result$U))
})

test_that("mergeMashData combines two datasets with identical columns", {
    d1 <- list(
        random = data.frame(
            a = 1:3,
            b = 4:6,
            row.names = c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G")
        )
    )
    d2 <- list(
        random = data.frame(
            a = 7:8,
            b = 9:10,
            row.names = c("chr1:400:A:G", "chr1:500:A:G")
        )
    )
    result <- mergeMashData(d1, d2)
    expect_equal(nrow(result$random), 5)
    expect_equal(ncol(result$random), 2)
    expect_equal(colnames(result$random), c("a", "b"))
    expect_equal(result$random$a, c(1, 2, 3, 7, 8))
})

test_that("mergeMashData handles NULL input", {
    d1 <- NULL
    d2 <- list(random = data.frame(a = 1:3))
    result <- mergeMashData(d1, d2)
    expect_equal(nrow(result$random), 3)
})

test_that("mergeMashData aligns different column names correctly", {
    d1 <- list(
        random = data.frame(
            a = 1:2,
            b = 3:4,
            row.names = c("chr1:100:A:G", "chr1:200:A:G")
        )
    )
    d2 <- list(
        random = data.frame(
            a = 5:6,
            c = 7:8,
            row.names = c("chr1:300:A:G", "chr1:400:A:G")
        )
    )
    result <- mergeMashData(d1, d2)
    expect_equal(nrow(result$random), 4)
    expect_true(all(c("a", "b", "c") %in% colnames(result$random)))
    expect_true(is.nan(result$random[3, "b"]))
    expect_true(is.nan(result$random[1, "c"]))
})

test_that("mergeMashData preserves data when one side is empty data.frame", {
    d1 <- list(random = data.frame(a = 1:3))
    d2 <- list(random = data.frame())
    result <- mergeMashData(d1, d2)
    expect_equal(nrow(result$random), 3)
})

test_that("mergeMashData preserves data when one side is NULL element", {
    d1 <- list(random = NULL)
    d2 <- list(random = data.frame(a = 1:3))
    result <- mergeMashData(d1, d2)
    expect_equal(nrow(result$random), 3)
})

test_that("mergeMashData handles multiple named elements", {
    d1 <- list(
        random = data.frame(
            a = 1:2,
            b = 3:4,
            row.names = c("chr1:100:A:G", "chr1:200:A:G")
        ),
        null = data.frame(x = 10:11, row.names = c("n1", "n2"))
    )
    d2 <- list(
        random = data.frame(
            a = 5:6,
            b = 7:8,
            row.names = c("chr1:300:A:G", "chr1:400:A:G")
        ),
        null = data.frame(x = 12:13, row.names = c("n3", "n4"))
    )
    result <- mergeMashData(d1, d2)
    expect_equal(nrow(result$random), 4)
    expect_equal(nrow(result$null), 4)
})

test_that("mergeMashData errors on duplicate variant ids across sides", {
    d1 <- list(
        random = data.frame(
            a = 1:2,
            row.names = c("chr1:100:A:G", "chr1:200:A:G")
        )
    )
    d2 <- list(
        random = data.frame(
            a = 3:4,
            row.names = c("chr1:200:A:G", "chr1:900:A:G")
        )
    )
    expect_error(mergeMashData(d1, d2), "duplicate variant ids")
})

test_that("mergeMashData uses one_data when res_data element has zero rows", {
    d1 <- list(random = data.frame(a = numeric(0), b = numeric(0)))
    d2 <- list(random = data.frame(a = 1:3, b = 4:6))
    result <- mergeMashData(d1, d2)
    expect_true(nrow(result$random) >= 3)
})

test_that("mashRandNullSample with z scores returns random and null", {
    set.seed(42)
    dat <- list(
        z = data.frame(
            cond1 = c(5, 0.1, 0.2, 0.3, 0.5, 6, 0.1, 0.2, 0.4, 0.3),
            cond2 = c(0.2, 0.3, 0.1, 0.5, 0.4, 0.1, 0.3, 0.2, 0.1, 0.5)
        )
    )
    result <- mashRandNullSample(
        dat,
        nRandom = 5,
        nNull = 3,
        excludeCondition = c(),
        seed = 123
    )
    expect_type(result, "list")
    expect_true("random" %in% names(result))
    expect_true("null" %in% names(result))
    expect_true("z" %in% names(result$random))
    expect_equal(nrow(result$random$z), 5)
})

test_that("mashRandNullSample with seed is reproducible", {
    dat <- list(
        z = data.frame(
            cond1 = c(0.1, 0.2, 0.3, 0.4, 0.5),
            cond2 = c(0.5, 0.4, 0.3, 0.2, 0.1)
        )
    )
    result1 <- mashRandNullSample(
        dat,
        nRandom = 3,
        nNull = 2,
        excludeCondition = c(),
        seed = 42
    )
    result2 <- mashRandNullSample(
        dat,
        nRandom = 3,
        nNull = 2,
        excludeCondition = c(),
        seed = 42
    )
    expect_equal(result1$random$z, result2$random$z)
})

test_that("mashRandNullSample NULL input returns NULL", {
    result <- mashRandNullSample(
        NULL,
        nRandom = 5,
        nNull = 3,
        excludeCondition = c()
    )
    expect_null(result)
})

test_that("mashRandNullSample warns when no null variants found (all abs_z > 2)", {
    dat <- list(
        z = data.frame(
            cond1 = c(5, 6, 7, 8, 9),
            cond2 = c(5, 6, 7, 8, 9)
        )
    )
    expect_warning(
        result <- mashRandNullSample(
            dat,
            nRandom = 3,
            nNull = 2,
            excludeCondition = c(),
            seed = 42
        ),
        "no variants are included in the null"
    )
    expect_equal(length(result$null), 0)
})

test_that("mashRandNullSample warns when not enough null data", {
    dat <- list(
        z = data.frame(
            cond1 = c(5, 6, 0.1),
            cond2 = c(5, 6, 0.1),
            cond3 = c(5, 6, 0.1)
        )
    )
    expect_warning(
        result <- mashRandNullSample(
            dat,
            nRandom = 2,
            nNull = 1,
            excludeCondition = c(),
            seed = 42
        ),
        "not enough null data"
    )
    expect_equal(length(result$null), 0)
})

test_that("mashRandNullSample with bhat/sbhat processes random and null samples", {
    dat <- list(
        bhat = data.frame(
            cond1 = c(0.1, 0.05, 0.02, 0.01, 0.03),
            cond2 = c(0.02, 0.01, 0.03, 0.05, 0.04)
        ),
        sbhat = data.frame(
            cond1 = c(0.1, 0.1, 0.1, 0.1, 0.1),
            cond2 = c(0.1, 0.1, 0.1, 0.1, 0.1)
        )
    )
    result <- mashRandNullSample(
        dat,
        nRandom = 3,
        nNull = 3,
        excludeCondition = c(),
        seed = 42
    )
    expect_equal(ncol(result$random$bhat), 2)
    expect_equal(nrow(result$random$bhat), 3)
    expect_true(length(result$null) > 0)
})

test_that("mashRandNullSample errors when excludeCondition not found (z path)", {
    dat <- list(
        z = data.frame(cond1 = 1:5, cond2 = 1:5)
    )
    expect_error(
        mashRandNullSample(
            dat,
            nRandom = 3,
            nNull = 2,
            excludeCondition = "nonexistent",
            seed = 42
        ),
        "excludeCondition are not present"
    )
})

test_that("mashRandNullSample errors when excludeCondition not found (bhat path)", {
    dat <- list(
        bhat = data.frame(cond1 = 1:5, cond2 = 1:5),
        sbhat = data.frame(cond1 = rep(1, 5), cond2 = rep(1, 5))
    )
    expect_error(
        mashRandNullSample(
            dat,
            nRandom = 3,
            nNull = 2,
            excludeCondition = "nonexistent",
            seed = 42
        ),
        "excludeCondition are not present"
    )
})

test_that("mashRandNullSample drops excluded condition by column name", {
    dat <- list(
        z = data.frame(
            cond1 = c(0.1, 0.2, 0.3, 0.4, 0.5),
            cond2 = c(0.5, 0.4, 0.3, 0.2, 0.1),
            cond3 = c(0.3, 0.3, 0.3, 0.3, 0.3)
        )
    )
    result <- mashRandNullSample(
        dat,
        nRandom = 3,
        nNull = 2,
        excludeCondition = "cond3",
        seed = 42
    )
    expect_equal(colnames(result$random$z), c("cond1", "cond2"))
    expect_equal(colnames(result$null$z), c("cond1", "cond2"))
    expect_false("cond3" %in% colnames(result$random$z))
})

test_that("mashRandNullSample extracts null data with z scores when enough null variants exist", {
    dat <- list(
        z = data.frame(
            cond1 = c(0.1, 0.3, 0.2, 0.5, 0.4, 0.1, 0.3, 0.2, 0.5, 0.4),
            cond2 = c(0.2, 0.1, 0.4, 0.3, 0.5, 0.2, 0.1, 0.4, 0.3, 0.5)
        )
    )
    result <- mashRandNullSample(
        dat,
        nRandom = 5,
        nNull = 4,
        excludeCondition = c(),
        seed = 42
    )
    expect_true("null" %in% names(result))
    expect_true("z" %in% names(result$null))
    expect_equal(nrow(result$null$z), 4)
    expect_equal(ncol(result$null$z), 2)
    expect_equal(nrow(result$random$z), 5)
})

test_that("mashRandNullSample null data capped at available null variants", {
    dat <- list(
        z = data.frame(
            cond1 = c(0.1, 0.2, 0.3),
            cond2 = c(0.1, 0.2, 0.3)
        )
    )
    result <- mashRandNullSample(
        dat,
        nRandom = 2,
        nNull = 100,
        excludeCondition = c(),
        seed = 42
    )
    expect_true(length(result$null) > 0)
    expect_equal(nrow(result$null$z), 3)
})

test_that("mashRandNullSample extracts null data with bhat/sbhat when enough null variants", {
    dat <- list(
        bhat = data.frame(
            cond1 = c(0.01, 0.02, 0.01, 0.03, 0.02, 0.01),
            cond2 = c(0.02, 0.01, 0.03, 0.01, 0.02, 0.01)
        ),
        sbhat = data.frame(
            cond1 = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
            cond2 = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
        )
    )
    result <- mashRandNullSample(
        dat,
        nRandom = 3,
        nNull = 4,
        excludeCondition = c(),
        seed = 42
    )
    expect_true("null" %in% names(result))
    expect_true("bhat" %in% names(result$null))
    expect_true("sbhat" %in% names(result$null))
    expect_equal(nrow(result$null$bhat), 4)
    expect_equal(nrow(result$null$sbhat), 4)
})

test_that("mashRandNullSample excludeCondition with numeric index errors on z path", {
    dat <- list(
        z = data.frame(
            cond1 = c(0.1, 0.2, 0.3, 0.4, 0.5),
            cond2 = c(0.5, 0.4, 0.3, 0.2, 0.1),
            cond3 = c(0.3, 0.3, 0.3, 0.3, 0.3)
        )
    )
    expect_error(
        mashRandNullSample(
            dat,
            nRandom = 3,
            nNull = 2,
            excludeCondition = 3,
            seed = 42
        ),
        "excludeCondition are not present"
    )
})

test_that("mashRandNullSample excludeCondition with numeric index errors on bhat path", {
    dat <- list(
        bhat = data.frame(
            cond1 = c(0.01, 0.02, 0.01, 0.03, 0.02),
            cond2 = c(0.02, 0.01, 0.03, 0.01, 0.02),
            cond3 = c(0.01, 0.01, 0.01, 0.01, 0.01)
        ),
        sbhat = data.frame(
            cond1 = c(0.1, 0.1, 0.1, 0.1, 0.1),
            cond2 = c(0.1, 0.1, 0.1, 0.1, 0.1),
            cond3 = c(0.1, 0.1, 0.1, 0.1, 0.1)
        )
    )
    expect_error(
        mashRandNullSample(
            dat,
            nRandom = 3,
            nNull = 2,
            excludeCondition = 3,
            seed = 42
        ),
        "excludeCondition are not present"
    )
})

test_that("mashRandNullSample caps random sample at available rows", {
    dat <- list(
        z = data.frame(
            cond1 = c(0.1, 0.2, 0.3),
            cond2 = c(0.2, 0.1, 0.3)
        )
    )
    result <- mashRandNullSample(
        dat,
        nRandom = 100,
        nNull = 2,
        excludeCondition = c(),
        seed = 42
    )
    expect_equal(nrow(result$random$z), 3)
})

test_that("filterMixtureComponents filters zero matrices", {
    U <- list(
        mat1 = matrix(
            c(1, 0.5, 0.5, 1),
            2,
            2,
            dimnames = list(c("A", "B"), c("A", "B"))
        ),
        mat2 = matrix(0, 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
        mat3 = matrix(
            c(0.8, 0.3, 0.3, 0.9),
            2,
            2,
            dimnames = list(c("A", "B"), c("A", "B"))
        )
    )
    w <- c(mat1 = 0.5, mat2 = 0.3, mat3 = 0.2)
    conditions_to_keep <- c("A", "B")

    result <- filterMixtureComponents(conditions_to_keep, U, w)

    # mat2 should be removed (all zeros)
    expect_true(!"mat2" %in% names(result$U))
    # weights should be rescaled to maintain sum
    expect_equal(sum(result$w), sum(w), tolerance = 1e-10)
})

test_that("filterMixtureComponents removes low weight components", {
    U <- list(
        mat1 = matrix(
            c(1, 0.5, 0.5, 1),
            2,
            2,
            dimnames = list(c("A", "B"), c("A", "B"))
        ),
        mat2 = matrix(
            c(0.8, 0.3, 0.3, 0.9),
            2,
            2,
            dimnames = list(c("A", "B"), c("A", "B"))
        )
    )
    w <- c(mat1 = 0.999, mat2 = 0.00001) # mat2 below default cutoff
    conditions_to_keep <- c("A", "B")

    result <- filterMixtureComponents(conditions_to_keep, U, w, wCutoff = 1e-04)
    expect_true(!"mat2" %in% names(result$U))
})

test_that("filterMixtureComponents errors on missing condition", {
    U <- list(
        mat1 = matrix(
            c(1, 0.5, 0.5, 1),
            2,
            2,
            dimnames = list(c("A", "B"), c("A", "B"))
        )
    )
    w <- c(mat1 = 1.0)
    expect_error(
        filterMixtureComponents(c("A", "C"), U, w),
        "not found in matrix"
    )
})

test_that("filterMixtureComponents subsets conditions", {
    U <- list(
        mat1 = matrix(
            c(1, 0.5, 0.2, 0.5, 1, 0.3, 0.2, 0.3, 1),
            3,
            3,
            dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
        )
    )
    w <- c(mat1 = 1.0)

    result <- filterMixtureComponents(c("A", "B"), U, w)
    expect_equal(nrow(result$U[[1]]), 2)
    expect_equal(ncol(result$U[[1]]), 2)
})

test_that(".mashSumStatsToMatrices: auto picks BETA+SE when present", {
    ss <- .mssm_makeQtlSumStats(function(i, n) {
        list(
            SNP = sprintf("chr1:%d:A:G", 100L * seq_len(n)),
            A1 = "A",
            A2 = "G",
            Z = rnorm(n),
            BETA = rnorm(n, sd = 0.1),
            SE = abs(rnorm(n, sd = 0.05)) + 0.01
        )
    })
    out <- pecotmr:::.mashSumStatsToMatrices(ss, "strong", inputScale = "auto")
    expect_equal(ncol(out$b), 2L) # 2 contexts
    # On BETA scale, Shat values should be the small SEs we generated.
    expect_true(all(out$s[out$s < 1000] < 1))
})

test_that(".mashSumStatsToMatrices: auto falls back to Z when no BETA/SE", {
    ss <- .mssm_makeQtlSumStats(function(i, n) {
        list(
            SNP = sprintf("chr1:%d:A:G", 100L * seq_len(n)),
            A1 = "A",
            A2 = "G",
            Z = rnorm(n)
        )
    })
    out <- pecotmr:::.mashSumStatsToMatrices(ss, "strong", inputScale = "auto")
    # Shat should be 1 on the Z scale.
    expect_true(all(out$s == 1 | out$s == 1000))
})

test_that(".mashSumStatsToMatrices: inputScale='beta' errors when BETA missing", {
    ss <- .mssm_makeQtlSumStats(function(i, n) {
        list(
            SNP = sprintf("chr1:%d:A:G", 100L * seq_len(n)),
            A1 = "A",
            A2 = "G",
            Z = rnorm(n)
        )
    })
    expect_error(
        pecotmr:::.mashSumStatsToMatrices(ss, "strong", inputScale = "beta"),
        "BETA and SE"
    )
})

test_that(".mashSumStatsToMatrices: inputScale='z' forces Z+1 even when BETA present", {
    ss <- .mssm_makeQtlSumStats(function(i, n) {
        list(
            SNP = sprintf("chr1:%d:A:G", 100L * seq_len(n)),
            A1 = "A",
            A2 = "G",
            Z = rnorm(n),
            BETA = rnorm(n, sd = 0.1),
            SE = abs(rnorm(n, sd = 0.05)) + 0.01
        )
    })
    out <- pecotmr:::.mashSumStatsToMatrices(ss, "strong", inputScale = "z")
    # Forced Z scale: Shat must be 1 everywhere except the NA fill (1000).
    expect_true(all(out$s == 1 | out$s == 1000))
})

test_that(".mashSumStatsToMatrices: errors when no usable scale", {
    ss <- .mssm_makeQtlSumStats(function(i, n) {
        list(
            SNP = sprintf("chr1:%d:A:G", 100L * seq_len(n)),
            A1 = "A",
            A2 = "G",
            N = rep(1000L, n)
        )
    }) # only N — no Z, no BETA/SE
    expect_error(
        pecotmr:::.mashSumStatsToMatrices(ss, "strong", inputScale = "auto"),
        "no usable scale"
    )
})

test_that(".mashSumStatsToMatrices: inputScale='z' errors when Z missing", {
    ss <- .mssm_makeQtlSumStats(function(i, n) {
        list(
            SNP = sprintf("chr1:%d:A:G", 100L * seq_len(n)),
            A1 = "A",
            A2 = "G",
            BETA = rnorm(n, sd = 0.1),
            SE = abs(rnorm(n, sd = 0.05)) + 0.01
        )
    }) # BETA/SE but no Z
    expect_error(
        pecotmr:::.mashSumStatsToMatrices(ss, "strong", inputScale = "z"),
        "carry a Z mcol"
    )
})

test_that(".mashObjectMatrices errors when marginal effects lack the required columns", {
    res <- QtlFineMappingResult(
        studyName = "s",
        context = "brain",
        trait = "t",
        method = "susie",
        entry = list(.sc_makeFineMappingRow(3))
    )
    # Force a marginal-effects table with no `context` column (as an mv/f-SuSiE
    # result trimmed of marginal sumstats would yield): mash cannot pivot it.
    testthat::local_mocked_bindings(
        marginalEffects = function(x, ...) {
            data.frame(variant_id = "v", beta = 1, se = 1)
        },
        .package = "pecotmr"
    )
    expect_error(
        pecotmr:::.mashObjectMatrices(
            res,
            inputScale = "auto",
            coverage = 0.95
        ),
        ">= 2 contexts"
    )
})

test_that(".mashObjectMatrices warns and pins the first method on a multi-method result", {
    res <- QtlFineMappingResult(
        studyName = c("s", "s"),
        context = c("brain", "liver"),
        trait = c("t", "t"),
        method = c("susie", "mvsusie"),
        entry = list(.sc_makeFineMappingRow(3), .sc_makeFineMappingRow(3))
    )
    expect_warning(
        pecotmr:::.mashObjectMatrices(
            res,
            inputScale = "auto",
            coverage = 0.95
        ),
        "multiple methods"
    )
})

test_that(".mashObjectPartitions errors when < 2 conditions remain after excludeCondition", {
    ss <- .mssm_makeQtlSumStats(function(i, n) {
        list(
            SNP = sprintf("chr1:%d:A:G", 100L * seq_len(n)),
            A1 = "A",
            A2 = "G",
            BETA = rnorm(n, sd = 0.1),
            SE = abs(rnorm(n, sd = 0.05)) + 0.01
        )
    })
    expect_error(
        pecotmr:::.mashObjectPartitions(
            ss,
            nRandom = 3,
            nNull = 3,
            excludeCondition = "liver",
            coverage = 0.95,
            inputScale = "auto",
            seed = 1
        ),
        "fewer than 2 conditions"
    )
})

test_that(".mashObjectPartitions warns when no variants match the independent-variant list", {
    ss <- .mssm_makeQtlSumStats(function(i, n) {
        list(
            SNP = sprintf("chr1:%d:A:G", 100L * seq_len(n)),
            A1 = "A",
            A2 = "G",
            BETA = rnorm(n, sd = 0.1),
            SE = abs(rnorm(n, sd = 0.05)) + 0.01
        )
    })
    # No independent variant matches the panel, so the random/null background is
    # empty for this object; the warning fires before the (empty) sampling.
    expect_warning(
        tryCatch(
            pecotmr:::.mashObjectPartitions(
                ss,
                nRandom = 2,
                nNull = 2,
                excludeCondition = character(0),
                coverage = 0.95,
                inputScale = "auto",
                seed = 1,
                independentVariants = c("chrX:999999:N:N")
            ),
            error = function(e) invisible(NULL)
        ),
        "no variants matched"
    )
})

test_that(".mashSumStatsToMatrices on bundled multicontext fixture: shape and rowname format", {
    data(qtlSumStatsMulticontextExample)
    out <- pecotmr:::.mashSumStatsToMatrices(
        qtlSumStatsMulticontextExample,
        "strong",
        inputScale = "auto"
    )
    expect_equal(ncol(out$b), 3L)
    expect_equal(colnames(out$b), c("brain", "blood", "muscle"))
    # One (study, trait) block, 200 variants -> 200 rows
    expect_equal(nrow(out$b), 200L)
    # Rownames are disambiguated by (study::trait::variant)
    expect_true(all(grepl("^study1::ENSG_example::", rownames(out$b))))
    # On the BETA scale, sbhat values are small; no NA fill needed since every
    # context has every variant
    expect_true(all(out$s < 1))
})

test_that(".mashSumStatsToMatrices fills missing variants with bhat=0 / sbhat=1000", {
    set.seed(42L)
    gh <- new(
        "GenotypeHandle",
        path = "/tmp/sketch.gds",
        format = "gds",
        snpInfo = data.frame(
            SNP = sprintf("chr1:%d:A:G", 100L * (1:5)),
            CHR = "1",
            BP = seq(100L, by = 100L, length.out = 5L),
            A1 = "A",
            A2 = "G",
            stringsAsFactors = FALSE
        ),
        nSamples = 50L,
        sampleIds = paste0("s", seq_len(50L)),
        pgenPtr = NULL
    )
    mkGr <- function(snpIds) {
        gr <- GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(
                start = seq(100L, by = 100L, length.out = length(snpIds)),
                width = 1L
            )
        )
        S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
            SNP = snpIds,
            A1 = "A",
            A2 = "G",
            Z = rnorm(length(snpIds)),
            BETA = rnorm(length(snpIds), sd = 0.1),
            SE = rep(0.05, length(snpIds))
        )
        gr
    }
    # ctx1 has all 5 variants; ctx2 has only the first 3
    ss <- QtlSumStats(
        studyName = c("s1", "s1"),
        context = c("ctx1", "ctx2"),
        trait = c("g1", "g1"),
        entry = list(
            mkGr(sprintf("chr1:%d:A:G", 100L * (1:5))),
            mkGr(sprintf("chr1:%d:A:G", 100L * (1:3)))
        ),
        genome = "hg19",
        ldSketch = gh,
        qcInfo = list(prebuilt = "synthetic")
    )
    out <- pecotmr:::.mashSumStatsToMatrices(ss, "strong", inputScale = "auto")
    expect_equal(dim(out$b), c(5L, 2L))
    expect_setequal(colnames(out$b), c("ctx1", "ctx2"))
    # In ctx2, the last 2 variants are missing -> bhat NA -> 0, shat NA -> 1000
    expect_equal(unname(out$b[4:5, "ctx2"]), c(0, 0))
    expect_equal(unname(out$s[4:5, "ctx2"]), c(1000, 1000))
    # ctx1 has them present
    expect_true(all(abs(out$b[, "ctx1"]) < 1))
    expect_true(all(out$s[, "ctx1"] < 1))
})

test_that(".mashSumStatsToMatrices disambiguates rownames across (study, trait) blocks", {
    set.seed(7L)
    gh <- new(
        "GenotypeHandle",
        path = "/tmp/sketch.gds",
        format = "gds",
        snpInfo = data.frame(
            SNP = sprintf("chr1:%d:A:G", 100L * (1:3)),
            CHR = "1",
            BP = c(100L, 200L, 300L),
            A1 = "A",
            A2 = "G",
            stringsAsFactors = FALSE
        ),
        nSamples = 50L,
        sampleIds = paste0("s", seq_len(50L)),
        pgenPtr = NULL
    )
    mkGr <- function(snpIds) {
        gr <- GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(
                start = seq(100L, by = 100L, length.out = length(snpIds)),
                width = 1L
            )
        )
        S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
            SNP = snpIds,
            A1 = "A",
            A2 = "G",
            Z = rnorm(length(snpIds)),
            BETA = rnorm(length(snpIds), sd = 0.1),
            SE = rep(0.05, length(snpIds))
        )
        gr
    }
    # Two (study, trait) blocks but they share SNP IDs v1, v2, v3 — without
    # the prefix the rbind would silently merge them.
    ss <- QtlSumStats(
        studyName = c("s1", "s1"),
        context = c("ctx1", "ctx1"),
        trait = c("g1", "g2"),
        entry = list(
            mkGr(sprintf("chr1:%d:A:G", 100L * (1:3))),
            mkGr(sprintf("chr1:%d:A:G", 100L * (1:3)))
        ),
        genome = "hg19",
        ldSketch = gh,
        qcInfo = list(prebuilt = "synthetic")
    )
    out <- pecotmr:::.mashSumStatsToMatrices(ss, "strong", inputScale = "auto")
    # 3 variants per block × 2 blocks = 6 rows
    expect_equal(nrow(out$b), 6L)
    expect_setequal(
        rownames(out$b),
        c(
            paste0("s1::g1::", sprintf("chr1:%d:A:G", 100L * (1:3))),
            paste0("s1::g2::", sprintf("chr1:%d:A:G", 100L * (1:3)))
        )
    )
})

test_that(".mashSumStatsToMatrices errors when entry lacks SNP mcol", {
    set.seed(8L)
    gh <- new(
        "GenotypeHandle",
        path = "/tmp/sketch.gds",
        format = "gds",
        snpInfo = data.frame(
            SNP = sprintf("chr1:%d:A:G", 100L * (1:3)),
            CHR = "1",
            BP = c(100L, 200L, 300L),
            A1 = "A",
            A2 = "G",
            stringsAsFactors = FALSE
        ),
        nSamples = 50L,
        sampleIds = paste0("s", seq_len(50L)),
        pgenPtr = NULL
    )
    gr <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = c(100L, 200L, 300L), width = 1L)
    )
    # NO SNP mcol — should trigger the variant-alignment error
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        A1 = "A",
        A2 = "G",
        Z = rnorm(3),
        BETA = rnorm(3, sd = 0.1),
        SE = rep(0.05, 3)
    )
    ss <- QtlSumStats(
        studyName = "s1",
        context = "ctx1",
        trait = "g1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = gh,
        qcInfo = list(prebuilt = "synthetic")
    )
    expect_error(
        pecotmr:::.mashSumStatsToMatrices(ss, "strong", inputScale = "auto"),
        "SNP"
    )
})

test_that(".mashSumStatsToMatrices on GwasSumStats: studies become columns", {
    set.seed(11L)
    gh <- new(
        "GenotypeHandle",
        path = "/tmp/sketch.gds",
        format = "gds",
        snpInfo = data.frame(
            SNP = sprintf("chr1:%d:A:G", 100L * (1:3)),
            CHR = "1",
            BP = c(100L, 200L, 300L),
            A1 = "A",
            A2 = "G",
            stringsAsFactors = FALSE
        ),
        nSamples = 50L,
        sampleIds = paste0("s", seq_len(50L)),
        pgenPtr = NULL
    )
    mkGr <- function(snpIds) {
        gr <- GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(
                start = seq(100L, by = 100L, length.out = length(snpIds)),
                width = 1L
            )
        )
        S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
            SNP = snpIds,
            A1 = "A",
            A2 = "G",
            Z = rnorm(length(snpIds)),
            BETA = rnorm(length(snpIds), sd = 0.1),
            SE = rep(0.05, length(snpIds))
        )
        gr
    }
    # Each study is its own (study) block; columns of the mash matrix are the
    # studies, so the result is block-diagonal with NA-fill off the diagonal.
    ss <- GwasSumStats(
        studyName = c("studyA", "studyB"),
        entry = list(
            mkGr(sprintf("chr1:%d:A:G", 100L * (1:3))),
            mkGr(sprintf("chr1:%d:A:G", 100L * (1:3)))
        ),
        genome = "hg19",
        ldSketch = gh,
        qcInfo = list(prebuilt = "synthetic")
    )
    out <- pecotmr:::.mashSumStatsToMatrices(ss, "strong", inputScale = "auto")
    expect_equal(ncol(out$b), 2L)
    expect_equal(colnames(out$b), c("studyA", "studyB"))
    # 2 studies x 3 variants = 6 rows; rownames prefixed by the study block key.
    expect_equal(nrow(out$b), 6L)
    expect_setequal(
        rownames(out$b),
        c(
            paste0("studyA::", sprintf("chr1:%d:A:G", 100L * (1:3))),
            paste0("studyB::", sprintf("chr1:%d:A:G", 100L * (1:3)))
        )
    )
    # studyA's rows are absent from studyB's column -> bhat 0 / shat 1000 fill.
    studyArows <- grep("^studyA::", rownames(out$b))
    expect_equal(unname(out$b[studyArows, "studyB"]), rep(0, 3))
    expect_equal(unname(out$s[studyArows, "studyB"]), rep(1000, 3))
    # On the BETA scale, the present cells carry the small generated SEs.
    expect_true(all(out$s[studyArows, "studyA"] < 1))
})

test_that(".mashSumStatsToMatrices errors on a non-SumStats input", {
    expect_error(
        pecotmr:::.mashSumStatsToMatrices(list(a = 1), "strong"),
        "must be a QtlSumStats or GwasSumStats"
    )
})

test_that(".mashSumStatsToMatrices errors when SumStats has empty QC info", {
    gh <- new(
        "GenotypeHandle",
        path = "/tmp/sketch.gds",
        format = "gds",
        snpInfo = data.frame(
            SNP = sprintf("chr1:%d:A:G", 100L * (1:3)),
            CHR = "1",
            BP = c(100L, 200L, 300L),
            A1 = "A",
            A2 = "G",
            stringsAsFactors = FALSE
        ),
        nSamples = 50L,
        sampleIds = paste0("s", seq_len(50L)),
        pgenPtr = NULL
    )
    gr <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = c(100L, 200L, 300L), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = sprintf("chr1:%d:A:G", 100L * (1:3)),
        A1 = "A",
        A2 = "G",
        Z = rnorm(3),
        BETA = rnorm(3, sd = 0.1),
        SE = rep(0.05, 3)
    )
    ss <- QtlSumStats(
        studyName = "s1",
        context = "c1",
        trait = "g1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = gh,
        qcInfo = list()
    ) # empty QC info
    expect_error(
        pecotmr:::.mashSumStatsToMatrices(ss, "strong"),
        "no QC info"
    )
})

test_that(".mashSumStatsToMatrices errors when SumStats has zero entries", {
    gh <- new(
        "GenotypeHandle",
        path = "/tmp/sketch.gds",
        format = "gds",
        snpInfo = data.frame(
            SNP = sprintf("chr1:%d:A:G", 100L * (1:3)),
            CHR = "1",
            BP = c(100L, 200L, 300L),
            A1 = "A",
            A2 = "G",
            stringsAsFactors = FALSE
        ),
        nSamples = 50L,
        sampleIds = paste0("s", seq_len(50L)),
        pgenPtr = NULL
    )
    ss <- QtlSumStats(
        studyName = character(0),
        context = character(0),
        trait = character(0),
        entry = list(),
        genome = "hg19",
        ldSketch = gh,
        varY = numeric(0),
        qcInfo = list(prebuilt = "synthetic")
    )
    expect_error(
        pecotmr:::.mashSumStatsToMatrices(ss, "strong"),
        "no entries"
    )
})

test_that("mergeMashData returns resData when oneData is NULL", {
    d1 <- list(random = data.frame(a = 1:3, b = 4:6))
    expect_equal(mergeMashData(d1, NULL), d1)
})

test_that("mergeMashData returns resData when oneData is an empty list", {
    d1 <- list(random = data.frame(a = 1:3))
    expect_equal(mergeMashData(d1, list()), d1)
})

test_that("qtlSumStatsFromZMatrix: one row per context, Z preserved verbatim", {
    z <- matrix(
        c(1.1, -2.2, 0.3, 0.4, -0.5, 0.6),
        nrow = 3,
        dimnames = list(
            c("chr1:100:A:G", "chr1:200:C:T", "chr2:50:A:T"),
            c("brain", "liver")
        )
    )
    qss <- qtlSumStatsFromZMatrix(z, studyName = "s1", ldSketch = .qszm_gh())
    expect_s4_class(qss, "QtlSumStats")
    # The variants span chr1 and chr2, so each context becomes one ELEMENT per
    # chromosome. The tuple is what stays 1:1 with a matrix column.
    expect_equal(nrow(qss), 4L)
    expect_setequal(as.character(qss$context), c("brain", "liver"))
    expect_equal(unique(as.character(qss$study)), "s1")
    expect_equal(unique(as.character(qss$trait)), "mash")
    # sumStats stitches a tuple's elements back together, so the Z column
    # still matches the input matrix column verbatim (values only; the mcols
    # column is unnamed whereas z[, j] carries the row ids as names).
    brain <- sumStats(qss, studyName = "s1", context = "brain", trait = "mash")
    liver <- sumStats(qss, studyName = "s1", context = "liver", trait = "mash")
    expect_equal(S4Vectors::mcols(brain)$Z, unname(z[, 1]))
    expect_equal(S4Vectors::mcols(liver)$Z, unname(z[, 2]))
})

test_that("qtlSumStatsFromZMatrix: decodes chrom/pos from ids, synthesises when they don't parse", {
    z <- matrix(
        rnorm(2),
        ncol = 1,
        dimnames = list(c("chr1:250:A:G", "not_a_variant"), "ctx")
    )
    qss <- qtlSumStatsFromZMatrix(z, studyName = "s1", ldSketch = .qszm_gh())
    # Stitched across whatever seqnames the ids decoded to, so the variant
    # order still matches the matrix rows.
    e <- sumStats(qss, studyName = "s1", context = "ctx", trait = "mash")
    expect_equal(GenomicRanges::start(e)[1], 250L) # decoded
    expect_true(GenomicRanges::start(e)[2] >= 1L) # synthetic fallback
    # un-parseable chrom falls back to chr1 (never NA)
    expect_false(any(is.na(as.character(GenomicRanges::seqnames(e)))))
    # SNP ids are carried through unchanged
    expect_equal(S4Vectors::mcols(e)$SNP, rownames(z))
})

test_that("qtlSumStatsFromZMatrix: NULL rownames get synthetic variant ids", {
    z <- matrix(rnorm(4), nrow = 2, dimnames = list(NULL, c("a", "b")))
    qss <- qtlSumStatsFromZMatrix(z, studyName = "s1", ldSketch = .qszm_gh())
    expect_equal(S4Vectors::mcols(qss[[1]])$SNP, c("var1", "var2"))
})

test_that("qtlSumStatsFromZMatrix: placeholders and pass-through qcInfo are set", {
    z <- matrix(rnorm(6), nrow = 3, dimnames = list(NULL, c("x", "y")))
    qss <- qtlSumStatsFromZMatrix(
        z,
        studyName = "s1",
        ldSketch = .qszm_gh(),
        n = 500L,
        a1 = "T",
        a2 = "C",
        role = "strong"
    )
    mc <- S4Vectors::mcols(qss[[1]])
    expect_equal(unique(mc$A1), "T")
    expect_equal(unique(mc$A2), "C")
    expect_equal(unique(mc$N), 500L)
    expect_equal(qcInfo(qss)$role, "strong")
    expect_equal(length(qcInfo(qss)$entryAudit), 2L) # one slot per context
})

test_that("qtlSumStatsFromZMatrix: columns can map to traits or context x trait pairs", {
    z <- matrix(rnorm(6), nrow = 3, dimnames = list(NULL, c("geneA", "geneB")))
    # columns as traits: constant context, one trait per column
    qss <- qtlSumStatsFromZMatrix(
        z,
        studyName = "s1",
        ldSketch = .qszm_gh(),
        context = "brain",
        trait = colnames(z)
    )
    expect_equal(as.character(qss$context), c("brain", "brain"))
    expect_equal(as.character(qss$trait), c("geneA", "geneB"))
    # columns as (context, trait) pairs
    qss2 <- qtlSumStatsFromZMatrix(
        z,
        studyName = "s1",
        ldSketch = .qszm_gh(),
        context = c("brain", "liver"),
        trait = c("geneA", "geneA")
    )
    expect_equal(as.character(qss2$context), c("brain", "liver"))
    expect_equal(as.character(qss2$trait), c("geneA", "geneA"))
})

test_that("qtlSumStatsFromZMatrix: a condition label of the wrong length errors", {
    z <- matrix(rnorm(6), nrow = 3, dimnames = list(NULL, c("a", "b")))
    expect_error(
        qtlSumStatsFromZMatrix(
            z,
            studyName = "s1",
            ldSketch = .qszm_gh(),
            trait = c("t1", "t2", "t3")
        ),
        "must be length 1 or ncol"
    )
})

test_that("qtlSumStatsFromZMatrix: rejects non-matrix input and unlabelled conditions", {
    expect_error(
        qtlSumStatsFromZMatrix(1:5, studyName = "s1", ldSketch = .qszm_gh()),
        "variants x conditions matrix"
    )
    # no colnames -> the default context = colnames(z) is NULL
    z <- matrix(rnorm(4), nrow = 2)
    expect_error(
        qtlSumStatsFromZMatrix(z, studyName = "s1", ldSketch = .qszm_gh()),
        "column names"
    )
})

test_that("qtlSumStatsFromBetaMatrix: one entry per context, BETA/SE/Z mcols set", {
    d <- .qszmBeta()
    qss <- qtlSumStatsFromBetaMatrix(
        d$bhat,
        d$shat,
        studyName = "s1",
        ldSketch = .qszm_gh()
    )
    expect_s4_class(qss, "QtlSumStats")
    expect_setequal(as.character(qss$context), c("brain", "liver"))
    # One element per (context, chromosome); sumStats stitches a context
    # back into the single GRanges matching the matrix column.
    mc <- S4Vectors::mcols(
        sumStats(qss, studyName = "s1", context = "brain", trait = "mash")
    )
    expect_true(all(c("BETA", "SE", "Z") %in% colnames(mc)))
    expect_equal(mc$BETA, unname(d$bhat[, 1]))
    expect_equal(mc$SE, unname(d$shat[, 1]))
    expect_equal(mc$Z, unname(d$bhat[, 1] / d$shat[, 1]))
})

test_that("qtlSumStatsFromBetaMatrix: feeds .mashSumStatsToMatrices on both scales", {
    d <- .qszmBeta()
    qss <- qtlSumStatsFromBetaMatrix(
        d$bhat,
        d$shat,
        studyName = "s1",
        ldSketch = .qszm_gh()
    )
    mb <- .mashSumStatsToMatrices(qss, "strong", inputScale = "beta")
    expect_equal(unname(mb$b), unname(d$bhat))
    expect_equal(unname(mb$s), unname(d$shat))
    mz <- .mashSumStatsToMatrices(qss, "strong", inputScale = "z")
    expect_equal(unname(mz$b), unname(d$bhat / d$shat))
})

test_that("qtlSumStatsFromBetaMatrix: validates matrices and matching dimensions", {
    d <- .qszmBeta()
    expect_error(
        qtlSumStatsFromBetaMatrix(1:5, d$shat, "s1", .qszm_gh()),
        "`bhat` must be a numeric"
    )
    expect_error(
        qtlSumStatsFromBetaMatrix(d$bhat, 1:5, "s1", .qszm_gh()),
        "`shat` must be a numeric"
    )
    expect_error(
        qtlSumStatsFromBetaMatrix(d$bhat, d$shat[1:2, ], "s1", .qszm_gh()),
        "identical dimensions"
    )
})

test_that("qtlSumStatsFromBetaMatrix: NULL rownames -> synthetic ids; placeholders set", {
    bhat <- matrix(rnorm(4), nrow = 2, dimnames = list(NULL, c("x", "y")))
    shat <- matrix(
        abs(rnorm(4)) + 0.1,
        nrow = 2,
        dimnames = list(NULL, c("x", "y"))
    )
    qss <- qtlSumStatsFromBetaMatrix(
        bhat,
        shat,
        studyName = "s1",
        ldSketch = .qszm_gh(),
        n = 500L,
        a1 = "T",
        a2 = "C",
        role = "strong"
    )
    mc <- S4Vectors::mcols(qss[[1]])
    expect_equal(mc$SNP, c("var1", "var2"))
    expect_equal(unique(mc$A1), "T")
    expect_equal(unique(mc$N), 500L)
    expect_equal(qcInfo(qss)$role, "strong")
})

test_that("mashInput: QtlSumStats path returns the flat b/s/z + XtX contract", {
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    out <- mashInput(list(geneA = ss), nRandom = 5, nNull = 5, seed = 1)
    for (k in c(
        "strong.b",
        "strong.s",
        "strong.z",
        "random.b",
        "random.s",
        "random.z",
        "null.b",
        "null.s",
        "null.z",
        "XtX"
    )) {
        expect_true(k %in% names(out), info = k)
    }
    nCond <- ncol(out$strong.z)
    expect_equal(nrow(out$random.z), 5L)
    expect_equal(nrow(out$null.z), 5L)
    # XtX is conditions x conditions and symmetric
    expect_equal(dim(out$XtX), c(nCond, nCond))
    expect_equal(unname(out$XtX), unname(t(out$XtX)))
})

test_that("mashInput: FineMappingResult strong = CS lead (max PIP) per condition", {
    fmr <- .mi_makeFmr()
    out <- mashInput(
        list(geneA = fmr),
        nRandom = 4,
        nNull = 4,
        coverage = 0.95,
        sigPCutoff = 0.5,
        seed = 7
    )
    expect_equal(colnames(out$strong.z), c("brain", "blood"))
    # brain CS lead = variant 3 (chr1:300), blood CS lead = variant 4 (chr1:400)
    expect_setequal(
        sub("_geneA$", "", rownames(out$strong.z)),
        c("chr1:300:A:G", "chr1:400:A:G")
    )
    expect_equal(nrow(out$random.z), 4L)
})

test_that("mashInput: a FineMappingResult with no credible set yields no strong", {
    vids <- paste0("chr1:", 100 * seq_len(6), ":A:G")
    noCsTL <- function(zvec) {
        data.frame(
            variant_id = vids,
            chrom = "1",
            pos = as.integer(100 * seq_len(6)),
            A1 = "G",
            A2 = "A",
            N = 1000,
            MAF = 0.2,
            marginal_beta = zvec * 0.05,
            marginal_se = 0.05,
            marginal_z = zvec,
            marginal_p = 2 * pnorm(-abs(zvec)),
            pip = rep(0.05, 6),
            posterior_mean = zvec * 0.05,
            posterior_sd = 0.02,
            cs_95 = rep("susie_0", 6),
            stringsAsFactors = FALSE
        )
    }
    e1 <- fineMappingRow(
        variantIds = vids,
        susieFit = list(x = 1),
        topLoci = noCsTL(rnorm(6))
    )
    e2 <- fineMappingRow(
        variantIds = vids,
        susieFit = list(x = 1),
        topLoci = noCsTL(rnorm(6))
    )
    fmr <- QtlFineMappingResult(
        studyName = c("s1", "s1"),
        context = c("brain", "blood"),
        trait = c("t1", "t1"),
        method = c("susie", "susie"),
        entry = list(e1, e2)
    )
    out <- mashInput(list(g = fmr), nRandom = 3, nNull = 3, seed = 2)
    expect_null(out$strong.z)
    expect_equal(nrow(out$random.z), 3L)
})

test_that("mashInput: multiple regions accumulate rows and disambiguate names", {
    fmr <- .mi_makeFmr()
    one <- mashInput(
        list(a = fmr),
        nRandom = 4,
        nNull = 4,
        sigPCutoff = 0.5,
        seed = 7
    )
    two <- mashInput(
        list(a = fmr, b = fmr),
        nRandom = 4,
        nNull = 4,
        sigPCutoff = 0.5,
        seed = 7
    )
    expect_equal(nrow(two$strong.z), 2L * nrow(one$strong.z))
    expect_equal(nrow(two$random.z), 2L * nrow(one$random.z))
    expect_false(any(duplicated(rownames(two$random.z))))
})

test_that("mashInput: zOnly = TRUE drops the .b/.s matrices", {
    fmr <- .mi_makeFmr()
    out <- mashInput(
        list(g = fmr),
        nRandom = 3,
        nNull = 3,
        zOnly = TRUE,
        sigPCutoff = 0.5,
        seed = 7
    )
    expect_false(any(grepl("\\.(b|s)$", names(out))))
    expect_true(all(c("strong.z", "random.z", "null.z", "XtX") %in% names(out)))
})

test_that("mashInput: excludeCondition drops the condition column everywhere", {
    # Exclude one of three conditions (mash needs >= 2), leaving brain + blood.
    data(qtlSumStatsMulticontextExample)
    ss <- qtlSumStatsMulticontextExample
    out <- mashInput(
        list(g = ss),
        nRandom = 4,
        nNull = 4,
        excludeCondition = "muscle",
        seed = 7
    )
    expect_false("muscle" %in% colnames(out$random.z))
    expect_setequal(colnames(out$random.z), c("brain", "blood"))
})

test_that("mashInput: rejects a non-SumStats/non-FineMapping element", {
    expect_error(
        mashInput(list(1L)),
        "must be a QtlSumStats, GwasSumStats, or FineMappingResult"
    )
    expect_error(mashInput(list()), "non-empty list")
})

test_that("mashInput: independentVariants restricts random/null pool (strong untouched)", {
    fmr <- .mi_makeFmr() # 6 variants chr1:100..600, 2 contexts
    vids <- paste0("chr1:", 100 * seq_len(6), ":A:G")
    indep <- vids[1:3] # only the first three are "independent"
    out <- mashInput(
        list(g = fmr),
        nRandom = 3,
        nNull = 3,
        independentVariants = indep,
        sigPCutoff = 0.5,
        seed = 7
    )
    # random / null variant ids (strip the "_g" region suffix) must lie in indep.
    rn_ids <- sub("_g$", "", rownames(out$random.z))
    expect_true(length(rn_ids) > 0L && all(rn_ids %in% indep))
    if (!is.null(out$null.z)) {
        expect_true(all(sub("_g$", "", rownames(out$null.z)) %in% indep))
    }
    # strong is NOT filtered -- its CS lead may lie outside the independent set.
    expect_true(nrow(out$strong.z) >= 1L)
})

test_that("mashInput: independentVariants matches across allele flip + chr prefix", {
    fmr <- .mi_makeFmr()
    # same positions as the first three variants, but ref/alt swapped and the
    # chr prefix dropped -- matchVariants must still match them (not a string cmp).
    indep_flipped <- c("1:100:G:A", "1:200:G:A", "1:300:G:A")
    out <- mashInput(
        list(g = fmr),
        nRandom = 3,
        nNull = 3,
        independentVariants = indep_flipped,
        sigPCutoff = 0.5,
        seed = 7
    )
    rn_ids <- sub("_g$", "", rownames(out$random.z))
    expect_true(length(rn_ids) > 0L)
    expect_true(all(rn_ids %in% paste0("chr1:", c(100, 200, 300), ":A:G")))
})

test_that("a NULL partition passes through as NULL", {
    expect_null(pecotmr:::.mashAsDataFrameOrNull(NULL))
    # A matrix becomes a base data.frame, keeping the variant-id rownames the
    # downstream combine step relies on.
    m <- matrix(
        1:4,
        2L,
        2L,
        dimnames = list(c("v1", "v2"), c("a", "b"))
    )
    out <- pecotmr:::.mashAsDataFrameOrNull(m)
    expect_s3_class(out, "data.frame")
    expect_equal(rownames(out), c("v1", "v2"))
})

test_that(".qtlSumStatsFromMatrix synthesises coords for unparseable ids", {
    local_mocked_bindings(
        parseVariantId = function(...) stop("bad"),
        .package = "pecotmr"
    )
    out <- pecotmr:::.qtlSumStatsFromMatrix(
        vids = c("weird1", "weird2"),
        nCond = 1L,
        studyName = "s1",
        ldSketch = NULL,
        context = "cA",
        trait = "g1",
        genome = "hg19",
        role = "input",
        mcolFn = function(i, k) list(Z = c(1, 2))
    )
    expect_s4_class(out, "QtlSumStats")
    gr <- out[[1L]]
    # Ids that carry no coordinates land on chr1 at their own row index, so
    # the object stays well-formed instead of failing to build.
    expect_equal(as.character(GenomicRanges::seqnames(gr)), c("chr1", "chr1"))
    expect_equal(GenomicRanges::start(gr), c(1L, 2L))
})

test_that("mashWrapper: argument guards fire", {
    expect_error(
        filterInvalidSummaryStat("not-a-list"),
        "datList.*Must be of type 'list'"
    )
    expect_error(
        filterInvalidSummaryStat(list(), sigPCutoff = 2),
        "sigPCutoff.*is not <= 1"
    )
    expect_error(
        filterInvalidSummaryStat(list(), filterByMissingRate = -1),
        "filterByMissingRate.*is not >= 0"
    )
    expect_error(
        filterMixtureComponents(conditionsToKeep = 1L, U = list()),
        "conditionsToKeep.*Must be of type 'character'"
    )
    expect_error(
        filterMixtureComponents("a", U = list(), wCutoff = -1),
        "wCutoff.*is not >= 0"
    )
    expect_error(mergeMashData("nope", list()), "Must be of type 'list'")
})

test_that(".mashConcatChr answers character(0) for no pieces", {
    expect_identical(pecotmr:::.mashConcatChr(list()), character(0))
})

test_that(".mashPrefixRownames leaves an absent or empty matrix alone", {
    expect_null(pecotmr:::.mashPrefixRownames(NULL, "ctx"))
    empty <- matrix(numeric(0), nrow = 0L, ncol = 1L)
    expect_identical(pecotmr:::.mashPrefixRownames(empty, "ctx"), empty)
})

test_that(".mashProcessZ skips the missing-rate filter when none is asked for", {
    # A NULL filterByMissingRate means "clean the values, keep every row".
    z <- matrix(c(1, NA, 3, 4), nrow = 2, dimnames = list(c("v1", "v2"), NULL))
    out <- pecotmr:::.mashProcessZ(z, NULL)
    expect_equal(nrow(out), 2L)
    expect_equal(unname(out[2, 1]), 0) # NA replaced by 0, row retained
})

test_that(".mashFilterZ leaves a list with no z-carrying partition alone", {
    # list_assign() would CREATE an absent component as NULL, so a list with
    # nothing to rewrite has to come back identical.
    datList <- list(other = list(bhat = matrix(1)))
    expect_identical(pecotmr:::.mashFilterZ(datList, 0.5, NULL), datList)
})
