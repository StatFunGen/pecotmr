context("h2_estimate_wrappers")

# ===========================================================================
# Helper functions to build test S4 objects
# ===========================================================================

# Two baseline annotations (binary and continuous) and one candidate over the
# variants of `ref`.
.h2PanelAnnotations <- function(ref) {
    set.seed(11)
    p <- length(ref)
    mat <- cbind(
        coding = rbinom(p, 1, 0.3),
        cons = runif(p),
        cand = rbinom(p, 1, 0.2)
    )
    AnnotationMatrix(
        mat,
        GenomicRanges::granges(ref, use.mcols = FALSE),
        data.frame(
            name = colnames(mat),
            tier = c("baseline", "baseline", "candidate"),
            type = c("binary", "continuous", "binary")
        ),
        genome = "hg19"
    )
}

# GWAS z-scores over the fixture's blocks with heritability `h2`:
# z = sqrt(n) R beta + e, e ~ N(0, R), block by block.
.h2PanelZ <- function(ld, n, h2, seed = 3) {
    set.seed(seed)
    M <- sum(map_int(ld, length))
    list_c(map(ld, function(x) {
        R <- ldMatrix(x)
        beta <- rnorm(nrow(R), sd = sqrt(h2 / M))
        e <- eigen(R, symmetric = TRUE)
        noise <- e$vectors %*% (sqrt(pmax(e$values, 0)) * rnorm(nrow(R)))
        as.vector(sqrt(n) * R %*% beta + noise)
    }))
}

# The fixture's LD blocks recording the A1 frequencies `af` (one per
# variant, in block order) as their AF.
.h2PanelLdWithAf <- function(fx, af) {
    sizes <- map_int(fx$ld, function(x) nrow(ldMatrix(x)))
    map2(fx$ld, split(af, rep(seq_along(sizes), sizes)), function(x, f) {
        variants <- as(x, "GRanges")
        LdMatrix(
            unname(ldMatrix(x)),
            variants = S4Vectors::`mcols<-`(
                variants,
                value = S4Vectors::DataFrame(
                    S4Vectors::mcols(variants),
                    AF = f
                )
            ),
            nSamples = ncol(fx$panel)
        )
    })
}

makeTestEigenRef <- function(nSnps = 20, nBlocks = 2) {
    snps_per_block <- nSnps / nBlocks
    blocks_gr <- GenomicRanges::GRanges(
        seqnames = rep("chr1", nBlocks),
        ranges = IRanges::IRanges(
            start = seq(1, by = snps_per_block * 100, length.out = nBlocks),
            width = snps_per_block * 100 - 1
        )
    )
    ld_blocks <- blocks_gr

    snp_info <- data.frame(
        SNP = paste0("rs", seq_len(nSnps)),
        CHR = rep("chr1", nSnps),
        BP = as.integer(seq(50, by = 100, length.out = nSnps)),
        A1 = rep("A", nSnps),
        A2 = rep("G", nSnps),
        stringsAsFactors = FALSE
    )

    eigen_list <- map(seq_len(nBlocks), function(b) {
        idx <- seq((b - 1) * snps_per_block + 1, b * snps_per_block)
        p <- length(idx)
        set.seed(42 + b)
        R <- matrix(0.3, p, p)
        diag(R) <- 1
        e <- eigen(R)
        list(values = e$values, vectors = e$vectors, snpIdx = as.integer(idx))
    })

    LdEigen(
        ldBlocks = ld_blocks,
        snpInfo = snp_info,
        nSamples = 500000L,
        inSample = FALSE,
        genome = "hg19",
        eigenList = eigen_list,
        eigenvalueTruncation = 1.0
    )
}

makeTestScoreRef <- function(
    nSnps = 20,
    nBlocks = 2,
    with_ld_matrices = FALSE
) {
    snps_per_block <- nSnps / nBlocks
    blocks_gr <- GenomicRanges::GRanges(
        seqnames = rep("chr1", nBlocks),
        ranges = IRanges::IRanges(
            start = seq(1, by = snps_per_block * 100, length.out = nBlocks),
            width = snps_per_block * 100 - 1
        )
    )
    ld_blocks <- blocks_gr

    snp_info <- data.frame(
        SNP = paste0("rs", seq_len(nSnps)),
        CHR = rep("chr1", nSnps),
        BP = as.integer(seq(50, by = 100, length.out = nSnps)),
        A1 = rep("A", nSnps),
        A2 = rep("G", nSnps),
        stringsAsFactors = FALSE
    )

    set.seed(99)
    ld_scores <- matrix(
        runif(nSnps, 1, 10),
        ncol = 1,
        dimnames = list(NULL, "base_l2")
    )
    ld_score_weights <- rep(1 / nSnps, nSnps)

    ld_matrix_list <- if (with_ld_matrices) {
        map(seq_len(nBlocks), function(b) {
            idx <- seq((b - 1) * snps_per_block + 1, b * snps_per_block)
            p <- length(idx)
            # AR(1) LD so per-SNP LD scores vary (equicorrelated R gives a
            # constant LD score, which makes the [1 | ldsc] design singular).
            R <- 0.5^abs(outer(seq_len(p), seq_len(p), "-"))
            list(R = R, snpIdx = as.integer(idx))
        })
    } else {
        list()
    }

    LdScore(
        ldBlocks = ld_blocks,
        snpInfo = snp_info,
        nSamples = 500000L,
        inSample = FALSE,
        genome = "hg19",
        ldScores = ld_scores,
        ldScoreWeights = ld_score_weights,
        ldMatrixList = ld_matrix_list
    )
}

makeTestAnnotations <- function(nSnps = 20) {
    set.seed(77)
    snp_ranges <- GenomicRanges::GRanges(
        seqnames = rep("chr1", nSnps),
        ranges = IRanges::IRanges(
            start = as.integer(seq(50, by = 100, length.out = nSnps)),
            width = 1L
        )
    )
    annot_mat <- matrix(
        c(rbinom(nSnps, 1, 0.5), rbinom(nSnps, 1, 0.3)),
        nrow = nSnps,
        ncol = 2
    )
    colnames(annot_mat) <- c("annot_A", "annot_B")
    annotation_meta <- data.frame(
        name = c("annot_A", "annot_B"),
        tier = c("baseline", "baseline"),
        type = c("binary", "binary"),
        stringsAsFactors = FALSE
    )
    AnnotationMatrix(annot_mat, snp_ranges, annotation_meta, genome = "hg19")
}

# ===========================================================================
# The LD statistic each setting reads (internal)
# ===========================================================================

test_that("the eigen statistic's settings refuse an LdScore", {
    score_ref <- makeTestScoreRef()
    for (method in c("lder", "hdl")) {
        expect_error(
            pecotmr:::.h2AssertSettingsRef(
                .h2GenomeWideSettings(method, FALSE, FALSE),
                score_ref,
                "estimateH2"
            ),
            "estimateH2: the eigen statistic requires an LdEigen",
            info = method
        )
    }
})

test_that("the variants' settings refuse an LdEigen", {
    eigen_ref <- makeTestEigenRef()
    for (method in c("gldsc", "sldsc")) {
        expect_error(
            pecotmr:::.h2AssertSettingsRef(
                .h2GenomeWideSettings(method, FALSE, FALSE),
                eigen_ref,
                "estimateH2"
            ),
            "the variant statistic requires an LdScore object, got LdEigen",
            info = method
        )
    }
})

test_that("each setting accepts the LD statistic it reads", {
    eigen_ref <- makeTestEigenRef()
    score_ref <- makeTestScoreRef()
    ok <- function(method, ref) {
        pecotmr:::.h2AssertSettingsRef(
            .h2GenomeWideSettings(method, FALSE, FALSE),
            ref,
            "estimateH2"
        )
    }
    expect_true(ok("lder", eigen_ref))
    expect_true(ok("hdl", eigen_ref))
    expect_true(ok("gldsc", score_ref))
    expect_true(ok("sldsc", score_ref))
})

# ===========================================================================
# computeLdScores for LdEigen
# ===========================================================================

test_that("computeLdScores LdEigen returns a non-negative base_l2 column", {
    eigen_ref <- makeTestEigenRef()
    result <- computeLdScores(eigen_ref)
    expect_true(is.matrix(result))
    expect_equal(ncol(result), 1)
    expect_equal(colnames(result), "base_l2")
    expect_equal(nrow(result), length(eigen_ref))
    expect_true(all(result[, "base_l2"] >= 0))
})

test_that("computeLdScores(LdEigen) with annotations: base + annotations", {
    eigen_ref <- makeTestEigenRef()
    annot <- makeTestAnnotations()
    result <- computeLdScores(eigen_ref, annotations = annot)
    expect_true(is.matrix(result))
    # base_l2 + 2 annotation columns
    expect_equal(ncol(result), 3)
    expect_equal(nrow(result), length(eigen_ref))
    expect_equal(colnames(result), c("base_l2", "annot_A", "annot_B"))
})

test_that("computeLdScores gives the same scores from LdEigen and LdScore", {
    # Both are ldsc's unbiased r^2 summed against each annotation. Weighted
    # by an annotation, the cross terms between eigenvectors do not cancel,
    # so an eigen shortcut that drops them fails here -- the LdEigen route
    # rebuilds R for those columns instead.
    ld <- makeTestLdMatrixBlocks(sizes = c(5L, 4L), rho = 0.6)
    scoreRef <- buildLdScore(ld, genome = "hg19")
    eigenRef <- buildLdEigen(ld, genome = "hg19")
    set.seed(5)
    mat <- cbind(binary = c(1, 0, 1, 1, 0, 0, 1, 0, 1), cont = runif(9))
    annot <- AnnotationMatrix(
        mat,
        GenomicRanges::granges(scoreRef, use.mcols = FALSE),
        data.frame(
            name = c("binary", "cont"),
            tier = "baseline",
            type = c("binary", "continuous")
        ),
        genome = "hg19"
    )
    viaScore <- computeLdScores(scoreRef, annotations = annot)
    expect_equal(
        computeLdScores(eigenRef, annotations = annot),
        viaScore,
        tolerance = 1e-12
    )
    # Both are the definition, block by block, at the fixture's n = 500.
    blockScores <- function(R, idx) {
        u <- R^2 - (1 - R^2) / (500 - 2)
        cbind(rowSums(u), u %*% mat[idx, ])
    }
    want <- rbind(
        blockScores(unname(ldMatrix(ld[[1]])), 1:5),
        blockScores(unname(ldMatrix(ld[[2]])), 6:9)
    )
    expect_equal(unname(viaScore), unname(want), tolerance = 1e-12)
})

test_that("LD scores read from a panel match its blocks' stored LD", {
    # The stored route squares each block's matrix; the panel route forms R's
    # rows from the panel's factor a stripe at a time, and the LdEigen route
    # from each block's SVD.
    fx <- .h2PanelFixture()
    lazyScore <- buildLdScore(fx$panel, ldBlocks = fx$blocks, genome = "hg19")
    storedScore <- buildLdScore(fx$ld, genome = "hg19")
    lazyEigen <- buildLdEigen(fx$panel, ldBlocks = fx$blocks, genome = "hg19")
    annot <- .h2PanelAnnotations(storedScore)
    want <- computeLdScores(storedScore, annotations = annot)
    expect_equal(
        computeLdScores(lazyScore, annotations = annot),
        want,
        tolerance = 1e-10
    )
    expect_equal(
        computeLdScores(lazyEigen, annotations = annot),
        want,
        tolerance = 1e-10
    )
    expect_equal(
        computeLdScores(lazyEigen)[, 1],
        mcols(lazyScore)$ldScores[, 1],
        tolerance = 1e-10
    )
})

test_that("estimators read annotations stored sparse as they read dense", {
    fx <- .h2PanelFixture()
    storedScore <- buildLdScore(fx$ld, genome = "hg19")
    storedEigen <- buildLdEigen(fx$ld, genome = "hg19")
    annot <- .h2PanelAnnotations(storedScore)
    sparse <- annot
    SummarizedExperiment::assay(sparse, "annotations") <- methods::as(
        SummarizedExperiment::assay(annot, "annotations"),
        "CsparseMatrix"
    )
    n <- 2000
    z <- .h2PanelZ(fx$ld, n, h2 = 0.5)
    expect_equal(
        computeLdScores(storedScore, annotations = sparse),
        computeLdScores(storedScore, annotations = annot)
    )
    runs <- list(
        lder = list(pecotmr:::eigenRegressionUnivariate, storedEigen),
        gldsc = list(pecotmr:::variantGlsUnivariate, storedScore),
        sldsc = list(pecotmr:::variantRegressionUnivariate, storedScore),
        hdl = list(pecotmr:::eigenLikelihoodUnivariate, storedEigen)
    )
    # HDL warns on this small panel that its intercept sits at its bound;
    # the estimates agree either way, which is what is checked here.
    for (method in names(runs)) {
        fn <- runs[[method]][[1L]]
        ref <- runs[[method]][[2L]]
        expect_equal(
            suppressWarnings(fn(z, n, ref, sparse)),
            suppressWarnings(fn(z, n, ref, annot)),
            tolerance = 1e-10,
            label = method
        )
    }
})

test_that("estimators read a panel's LD on demand as they read stored LD", {
    # The same blocks two ways: read from the genotype panel one block at a
    # time, and stored as LdMatrix blocks. The middle block has more variants
    # than the panel has samples, so its LD has 40 zero eigenvalues that the
    # SVD of the panel's factor does not list; LDER and HDL count those
    # directions without them. g-LDSC forms each block's Omega from the
    # factor rather than from R.
    fx <- .h2PanelFixture()
    lazyScore <- buildLdScore(fx$panel, ldBlocks = fx$blocks, genome = "hg19")
    storedScore <- buildLdScore(fx$ld, genome = "hg19")
    lazyEigen <- buildLdEigen(fx$panel, ldBlocks = fx$blocks, genome = "hg19")
    storedEigen <- buildLdEigen(fx$ld, genome = "hg19")
    annot <- .h2PanelAnnotations(storedScore)
    n <- 2000
    z <- .h2PanelZ(fx$ld, n, h2 = 0.5)
    runs <- list(
        lder = list(
            pecotmr:::eigenRegressionUnivariate,
            lazyEigen,
            storedEigen
        ),
        gldsc = list(pecotmr:::variantGlsUnivariate, lazyScore, storedScore),
        sldsc = list(
            pecotmr:::variantRegressionUnivariate,
            lazyScore,
            storedScore
        )
    )
    for (method in names(runs)) {
        fn <- runs[[method]][[1L]]
        for (an in list(NULL, annot)) {
            lazy <- fn(z, n, runs[[method]][[2L]], an)
            stored <- fn(z, n, runs[[method]][[3L]], an)
            expect_equal(lazy, stored, tolerance = 1e-8, label = method)
        }
    }
    # The eigen likelihood is the same function either way, and so its fit
    # agrees to rounding: here h2 to 1e-13 and its jackknife SE to 1e-9,
    # well inside the tolerance.
    eigenDesign <- function(ref) {
        data <- map(
            pecotmr:::.ldEigenBlocks(ref),
            pecotmr:::.h2EigenBlockData,
            eigenRef = ref,
            z = z,
            n = n,
            baselineMat = NULL,
            candMat = NULL
        )
        pecotmr:::.h2EigenDesign(
            data,
            length(ref),
            NULL,
            .h2GenomeWideSettings("hdl", FALSE, FALSE)@procedures
        )
    }
    lazyDesign <- eigenDesign(lazyEigen)
    storedDesign <- eigenDesign(storedEigen)
    expect_equal(length(lazyDesign$lam), length(storedDesign$lam))
    proc <- .h2GenomeWideSettings("hdl", FALSE, FALSE)@procedures
    ll <- function(design, par) {
        pecotmr:::.h2Likelihood(
            "variance",
            pecotmr:::.h2VarianceData(design, n, 100, proc),
            par
        )$value
    }
    for (par in list(c(0.3, 1), c(0.5, 1.2))) {
        expect_equal(
            ll(lazyDesign, par),
            ll(storedDesign, par),
            tolerance = 1e-10
        )
    }
    eigenFit <- pecotmr:::eigenLikelihoodUnivariate
    expect_equal(
        suppressWarnings(eigenFit(z, n, lazyEigen, NULL)),
        suppressWarnings(eigenFit(z, n, storedEigen, NULL)),
        tolerance = 1e-6
    )
    # Stratified, this fixture (a 100-sample panel, so the intercept sits at
    # its bound) leaves the likelihood flat along the base / annotation
    # ridge: the two optima have the same likelihood to 1e-7 but sit apart
    # along it, so only the total heritability is pinned down.
    stratifiedH2 <- function(ref) {
        suppressWarnings(eigenFit(z, n, ref, annot))$h2
    }
    expect_equal(
        stratifiedH2(lazyEigen),
        stratifiedH2(storedEigen),
        tolerance = 1e-3
    )
    # The truncated decompositions HDL conventionally uses agree too.
    lazy9 <- buildLdEigen(
        fx$panel,
        ldBlocks = fx$blocks,
        genome = "hg19",
        eigenvalueTruncation = 0.9
    )
    stored9 <- buildLdEigen(fx$ld, genome = "hg19", eigenvalueTruncation = 0.9)
    expect_equal(stratifiedH2(lazy9), stratifiedH2(stored9), tolerance = 1e-3)
})

# ===========================================================================
# computeLdScores for LdScore
# ===========================================================================

test_that("computeLdScores(LdScore) without annotations: stored ld_scores", {
    score_ref <- makeTestScoreRef()
    result <- computeLdScores(score_ref)
    expect_identical(result, mcols(score_ref)$ldScores)
})

test_that("computeLdScores(LdScore) errors on annotations without matrices", {
    score_ref <- makeTestScoreRef(with_ld_matrices = FALSE)
    annot <- makeTestAnnotations()
    expect_error(
        computeLdScores(score_ref, annotations = annot),
        "need each LD block's LD"
    )
})

# ===========================================================================
# h2EstimateToSldscTrait
# ===========================================================================

test_that("h2EstimateToSldscTrait returns correct list structure", {
    h2_obj <- makeTestH2Estimate(with_enrichment = TRUE)
    result <- h2EstimateToSldscTrait(h2_obj)

    expected_names <- c(
        "categories",
        "tau",
        "tauSe",
        "enrichment",
        "enrichmentSe",
        "enrichmentP",
        "propH2",
        "propSnps",
        "h2g",
        "annotationJackknifeCoefs",
        "nBlocks"
    )
    expect_true(all(expected_names %in% names(result)))
    expect_equal(length(result), length(expected_names))

    expect_equal(result$categories, c("annot1", "annot2"))
    expect_equal(result$h2g, 0.3)
    expect_true(is.matrix(result$annotationJackknifeCoefs))
    expect_equal(result$nBlocks, 10L)
    expect_equal(names(result$tau), c("annot1", "annot2"))
    expect_equal(names(result$enrichment), c("annot1", "annot2"))
})

test_that("h2EstimateToSldscTrait errors on non-H2Estimate input", {
    expect_error(
        h2EstimateToSldscTrait(list(h2 = 0.3)),
        "Must inherit from class 'H2Estimate'"
    )
})

test_that("h2EstimateToSldscTrait errors when enrichment is NULL", {
    h2_obj <- makeTestH2Estimate(with_enrichment = FALSE)
    expect_error(
        h2EstimateToSldscTrait(h2_obj),
        "no enrichment"
    )
})

test_that("h2EstimateToSldscTrait converts one estimate", {
    two <- c(makeTestH2Estimate(), `$<-`(makeTestH2Estimate(), "study", "b"))
    expect_error(h2EstimateToSldscTrait(two), "holds 2 estimates")
    expect_identical(
        h2EstimateToSldscTrait(two[two$study == "b"]),
        h2EstimateToSldscTrait(makeTestH2Estimate())
    )
})

test_that("h2EstimateToSldscTrait makes a 1-row dummy for NULL jackknife", {
    h2_obj <- makeTestH2Estimate(with_enrichment = TRUE)
    h2_obj$jackknife <- list(NULL)
    result <- h2EstimateToSldscTrait(h2_obj)
    expect_equal(nrow(result$annotationJackknifeCoefs), 1)
    expect_equal(result$nBlocks, 1L)
    expect_equal(
        colnames(result$annotationJackknifeCoefs),
        c("annot1", "annot2")
    )
})

# ===========================================================================
# Helper: GwasSumStats matched to a reference panel
# ===========================================================================

makeTestSumstatsForRef <- function(
    ref,
    traitName = "test",
    varY = NA_real_,
    ...
) {
    n_snps <- length(ref)
    set.seed(123)
    # An LdStatistic is a GRanges now: coordinates come from the ranges and
    # the alleles from mcols, rather than from a parallel snpInfo frame.
    md <- S4Vectors::mcols(ref, use.names = FALSE)
    df <- data.frame(
        SNP = names(ref),
        CHR = sub("^chr", "", as.character(GenomicRanges::seqnames(ref))),
        BP = GenomicRanges::start(ref),
        A1 = as.character(md$A1),
        A2 = as.character(md$A2),
        Z = rnorm(n_snps),
        N = rep(50000, n_snps),
        stringsAsFactors = FALSE
    )
    GwasSumStats(
        studyName = traitName,
        entry = list(.dfToGwasGr(df)),
        genome = "hg19",
        ldReference = makeTestGwasGenotypeHandle(),
        varY = varY,
        ...
    )
}

# ===========================================================================
# estimateH2 end-to-end tests
# ===========================================================================

test_that("estimateH2 returns an H2Estimate for lder, gldsc and hdl", {
    eigen_ref <- makeTestEigenRef()
    refs <- list(
        lder = eigen_ref,
        gldsc = makeTestScoreRef(with_ld_matrices = TRUE),
        hdl = eigen_ref
    )
    for (method in names(refs)) {
        ss <- makeTestSumstatsForRef(refs[[method]])
        result <- estimateH2(
            ss,
            refs[[method]],
            param = .h2MethodSettings(method, NULL, refs[[method]])
        )

        expect_true(is(result, "H2Estimate"), info = method)
        expect_true(is.numeric(heritabilityEstimate(result)), info = method)
        expect_true(is.numeric(heritabilityEstimateSe(result)), info = method)
        expect_true(is.numeric(heritabilityIntercept(result)), info = method)
        expect_true(is.numeric(heritabilityInterceptSe(result)), info = method)
        expect_identical(
            h2EstimateParam(result),
            .h2GenomeWideSettings(
                method,
                FALSE,
                identical(method, "lder") && inSample(refs[[method]])
            ),
            info = method
        )
        expect_equal(nSnps(result), nSnps(ss), info = method)
        expect_equal(studyName(result), "test", info = method)
        expect_identical(traitNames(result), character(0), info = method)
    }
})

test_that("a genome-wide estimate is one row for the study", {
    ref <- makeTestEigenRef()
    ss <- makeTestSumstatsForRef(ref)
    res <- estimateH2(ss, ref)
    expect_identical(length(res), 1L)
    # It covers the whole reference, so its element holds no ranges; the row
    # takes the reference's genome.
    expect_length(res[[1]], 0L)
    expect_identical(
        unname(GenomeInfoDb::genome(res)),
        unname(GenomeInfoDb::genome(ref))
    )
    expect_equal(
        heritabilityEstimate(res),
        eigenLikelihoodUnivariate(z(ss), 50000, ref)$h2
    )
})

# ===========================================================================
# estimateH2 input alignment
#
# The estimators index z and the annotation rows by each LD block's snpIdx,
# matching nothing on variant id, so a misaligned input used to return a
# silently wrong h2 instead of failing.
# ===========================================================================

# Rebuild `ss`'s single entry with the mcols kept and the ranges replaced.
.reposition <- function(ss, gr) {
    entry <- sumStats(ss)
    S4Vectors::mcols(gr) <- S4Vectors::mcols(entry)
    GwasSumStats(
        studyName = "test",
        entry = list(gr),
        genome = "hg19",
        ldReference = makeTestGwasGenotypeHandle()
    )
}

test_that("estimateH2 rejects a sumstats with a different variant count", {
    ref <- makeTestEigenRef()
    ss <- makeTestSumstatsForRef(ref)
    short <- GwasSumStats(
        studyName = "test",
        entry = list(sumStats(ss)[1:10]),
        genome = "hg19",
        ldReference = makeTestGwasGenotypeHandle()
    )
    expect_error(
        estimateH2(short, ref),
        "carries 10 variant\\(s\\) but the LD reference carries 20"
    )
})

test_that("estimateH2 rejects a sumstats at different coordinates", {
    ref <- makeTestEigenRef()
    ss <- makeTestSumstatsForRef(ref)
    shifted <- .reposition(
        ss,
        GenomicRanges::GRanges(
            "chr1",
            IRanges::IRanges(GenomicRanges::start(ref) + 1L, width = 1L)
        )
    )
    expect_error(
        estimateH2(shifted, ref),
        "disagree at 20 of 20 position"
    )
})

test_that("estimateH2 refuses sumstats with the right variants out of order", {
    ref <- makeTestEigenRef()
    ss <- makeTestSumstatsForRef(ref)
    reordered <- .reposition(
        ss,
        GenomicRanges::GRanges(
            "chr1",
            IRanges::IRanges(rev(GenomicRanges::start(ref)), width = 1L)
        )
    )
    expect_error(
        estimateH2(reordered, ref),
        "same variants in the same order"
    )
})

test_that("estimateStratifiedH2 rejects annotations short of the reference", {
    ref <- makeTestEigenRef()
    ss <- makeTestSumstatsForRef(ref)
    annot <- makeTestAnnotations(nSnps = 10)
    expect_error(
        estimateStratifiedH2(ss, ref, annot),
        "estimateStratifiedH2: `annotations` carries 10 variant\\(s\\)"
    )
})

test_that("estimateStratifiedH2 refuses annotations elsewhere or reordered", {
    ref <- makeTestEigenRef()
    ss <- makeTestSumstatsForRef(ref)
    annot <- makeTestAnnotations()
    # The right count of rows is not enough: the estimators read annotation
    # rows by the reference's variant order.
    relocate <- function(starts) {
        AnnotationMatrix(
            SummarizedExperiment::assay(annot, "annotations"),
            GenomicRanges::GRanges(
                "chr1",
                IRanges::IRanges(starts, width = 1L)
            ),
            as.data.frame(SummarizedExperiment::colData(annot)),
            genome = "hg19"
        )
    }
    expect_error(
        estimateStratifiedH2(
            ss,
            ref,
            relocate(GenomicRanges::start(ref) + 1L)
        ),
        "`annotations` and the LD reference disagree at 20 of 20 position"
    )
    expect_error(
        estimateStratifiedH2(
            ss,
            ref,
            relocate(rev(GenomicRanges::start(ref)))
        ),
        "`annotations` and the LD reference disagree"
    )
})

test_that("estimateH2 rejects a genome build mismatch", {
    ref <- makeTestEigenRef()
    ss <- makeTestSumstatsForRef(ref)
    hg38 <- GwasSumStats(
        studyName = "test",
        entry = list(sumStats(ss)),
        genome = "hg38",
        ldReference = makeTestGwasGenotypeHandle()
    )
    expect_error(estimateH2(hg38, ref), "build mismatch")
})

test_that("estimateH2 accepts an LD reference that records no genome build", {
    ref <- makeTestEigenRef()
    ss <- makeTestSumstatsForRef(ref)
    # An LD reference often records no build, so one built from it
    # names none. An unknown build must not read as a mismatch.
    GenomeInfoDb::genome(ref) <- NA_character_
    expect_s4_class(estimateH2(ss, ref), "H2Estimate")
})

# ===========================================================================
# estimateStratifiedH2: each genome-wide fit from its settings
# ===========================================================================

test_that("estimateStratifiedH2 defaults to stratified HDL's settings", {
    ref <- makeTestEigenRef()
    ss <- makeTestSumstatsForRef(ref)
    annot <- makeTestAnnotations()
    res <- suppressWarnings(estimateStratifiedH2(ss, ref, annot))
    expect_identical(
        h2EstimateParam(res),
        .h2GenomeWideSettings("hdl", TRUE, inSample(ref))
    )
    expect_s3_class(stratifiedHeritabilityEnrichment(res), "data.frame")
    # An empty list carries no settings, as everywhere a Param is taken.
    expect_identical(
        suppressWarnings(estimateStratifiedH2(ss, ref, annot, param = list())),
        res
    )
})

test_that("estimateStratifiedH2 runs the fit its settings choose", {
    eigenRef <- makeTestEigenRef()
    scoreRef <- makeTestScoreRef(with_ld_matrices = TRUE)
    annot <- makeTestAnnotations()
    fits <- list(
        lder = eigenRegressionUnivariate,
        hdl = eigenLikelihoodUnivariate,
        gldsc = variantGlsUnivariate,
        sldsc = variantRegressionUnivariate
    )
    for (method in names(fits)) {
        ref <- if (method %in% c("lder", "hdl")) eigenRef else scoreRef
        ss <- makeTestSumstatsForRef(ref)
        settings <- .h2MethodSettings(method, annot, ref)
        res <- suppressWarnings(
            estimateStratifiedH2(ss, ref, annot, param = settings)
        )
        direct <- suppressWarnings(fits[[method]](
            z(ss),
            median(sampleSize(ss)),
            ref,
            annot,
            param = settings
        ))
        expect_identical(h2EstimateParam(res), settings, info = method)
        expect_identical(heritabilityEstimate(res), direct$h2, info = method)
        expect_identical(
            heritabilityEstimateSe(res),
            direct$h2Se,
            info = method
        )
    }
})

test_that("a setting that is not the published one reaches the fit", {
    ref <- makeTestScoreRef(with_ld_matrices = TRUE)
    ss <- makeTestSumstatsForRef(ref)
    annot <- makeTestAnnotations()
    sldsc <- .h2MethodSettings("sldsc", annot, ref)
    plain <- initialize(
        sldsc,
        procedures = initialize(sldsc@procedures, ldRedundancy = FALSE)
    )
    h2 <- function(param) {
        heritabilityEstimate(suppressWarnings(
            estimateStratifiedH2(ss, ref, annot, param = param)
        ))
    }
    expect_false(isTRUE(all.equal(h2(plain), h2(sldsc))))
})

test_that("estimateStratifiedH2 refuses settings no genome-wide fit takes", {
    ref <- makeTestEigenRef()
    ss <- makeTestSumstatsForRef(ref)
    annot <- makeTestAnnotations()
    lder <- .h2GenomeWideSettings("lder", TRUE, FALSE)
    hdl <- .h2GenomeWideSettings("hdl", TRUE, FALSE)
    sldsc <- .h2GenomeWideSettings("sldsc", TRUE, FALSE)
    refused <- function(param, pattern) {
        expect_error(
            estimateStratifiedH2(ss, ref, annot, param = param),
            pattern
        )
    }
    # LAVA's moments are a local setting.
    refused(
        .h2LocalSettings("lava"),
        "estimateStratifiedH2: a genome-wide estimate fits the variant"
    )
    refused(
        initialize(hdl, noise = "tied"),
        "fits a free or floored noise term"
    )
    refused(
        initialize(lder, inference = "information"),
        "reports the jackknife over the reference's LD blocks"
    )
    refused(
        initialize(sldsc, referenceNoise = "none"),
        "variant regression reads referenceNoise = \"unbiasedR2\""
    )
    refused(
        initialize(lder, truncation = H2TruncationParam()),
        "keeps every direction its LdEigen holds"
    )
    # A least number of directions, which a block's fit would not read.
    refused(
        initialize(
            lder,
            truncation = H2TruncationParam(fraction = 1, minimum = 3)
        ),
        "keeps every direction its LdEigen holds"
    )
    refused(
        initialize(
            lder,
            truncation = H2TruncationParam(fraction = 1, maxPerSample = 0.75)
        ),
        "maxPerSample = NULL"
    )
    refused(
        initialize(
            lder,
            procedures = initialize(lder@procedures, jackknifeBlocks = 10)
        ),
        "leave `jackknifeBlocks` unset"
    )
    refused(
        initialize(
            hdl,
            procedures = initialize(
                hdl@procedures,
                mustBeatNull = TRUE
            )
        ),
        "does not take mustBeatNull: set it to FALSE"
    )
})

test_that("estimateStratifiedH2 checks its reference and arguments", {
    eigenRef <- makeTestEigenRef()
    scoreRef <- makeTestScoreRef(with_ld_matrices = TRUE)
    annot <- makeTestAnnotations()
    # LDER's settings, the default, read an LdEigen.
    expect_error(
        estimateStratifiedH2(makeTestSumstatsForRef(scoreRef), scoreRef, annot),
        "the eigen statistic requires an LdEigen object, got LdScore"
    )
    ss <- makeTestSumstatsForRef(eigenRef)
    expect_error(
        estimateStratifiedH2(ss, eigenRef, annot, param = list(noise = "free")),
        "must be built with H2EstimateParam\\(\\)"
    )
    expect_error(estimateStratifiedH2(ss, eigenRef, NULL), "annotations")
})

test_that("estimateH2 no longer takes annotations, a method or its options", {
    ref <- makeTestEigenRef()
    ss <- makeTestSumstatsForRef(ref)
    expect_error(
        estimateH2(ss, ref, annotations = makeTestAnnotations()),
        "unused argument"
    )
    expect_error(estimateH2(ss, ref, method = "lder"), "unused argument")
    expect_error(
        estimateH2(ss, ref, estimatorArgs = list(nIter = 2)),
        "unused argument"
    )
    expect_error(
        estimateH2(ss, ref, param = list(criterion = "moments")),
        "must be built with H2EstimateParam\\(\\)"
    )
})

test_that("estimateH2 fits the settings it is given and records them", {
    ref <- makeTestScoreRef(with_ld_matrices = TRUE)
    ss <- makeTestSumstatsForRef(ref)
    sldsc <- .h2GenomeWideSettings("sldsc", FALSE, FALSE)
    settings <- methods::initialize(
        sldsc,
        procedures = methods::initialize(sldsc@procedures, iterations = 5)
    )
    res <- estimateH2(ss, ref, param = settings)
    expect_identical(h2EstimateParam(res), settings)
    expect_identical(
        heritabilityEstimate(res),
        variantRegressionUnivariate(
            z(ss),
            median(sampleSize(ss)),
            ref,
            param = settings
        )$h2
    )
    # Without settings, HDL's.
    eigenRef <- makeTestEigenRef()
    byDefault <- estimateH2(makeTestSumstatsForRef(eigenRef), eigenRef)
    expect_identical(
        h2EstimateParam(byDefault),
        .h2GenomeWideSettings("hdl", FALSE, inSample(eigenRef))
    )
})

# ===========================================================================
# S-LDSC (Finucane et al. 2015)
#
# The original stratified LD score regression, added alongside the three
# newer estimators. Unlike g-LDSC it regresses chi2 on per-variant LD scores
# rather than working from the per-block LD matrices, so an unstratified run
# needs only the scores.
# ===========================================================================

test_that("unstratified S-LDSC runs without the per-block LD matrices", {
    # This is the property that separates it from g-LDSC: a scores-only
    # LdScore is enough, which is why it is the cheap option.
    noMat <- makeTestScoreRef(with_ld_matrices = FALSE)
    ss <- makeTestSumstatsForRef(noMat)
    expect_error(
        estimateH2(ss, noMat, param = .h2MethodSettings("gldsc", NULL, noMat)),
        "The variant GLS needs each LD block's LD"
    )
    res <- estimateH2(
        ss,
        noMat,
        param = .h2MethodSettings("sldsc", NULL, noMat)
    )
    expect_s4_class(res, "H2Estimate")
    expect_identical(
        h2EstimateParam(res),
        .h2GenomeWideSettings("sldsc", FALSE, FALSE)
    )
    expect_true(is.finite(heritabilityEstimate(res)))
})

test_that("stratified S-LDSC needs the per-block LD matrices", {
    noMat <- makeTestScoreRef(with_ld_matrices = FALSE)
    ss <- makeTestSumstatsForRef(noMat)
    annot <- makeTestAnnotations(nSnps = length(noMat))
    expect_error(
        estimateStratifiedH2(
            ss,
            noMat,
            annot,
            param = .h2MethodSettings("sldsc", annot, noMat)
        ),
        "stratified variant regression needs each LD block's LD"
    )
})

test_that("S-LDSC reports an enrichment table when stratified", {
    ref <- makeTestScoreRef(with_ld_matrices = TRUE)
    ss <- makeTestSumstatsForRef(ref)
    annot <- makeTestAnnotations(nSnps = length(ref))
    res <- estimateStratifiedH2(
        ss,
        ref,
        annot,
        param = .h2MethodSettings("sldsc", annot, ref)
    )

    enr <- stratifiedHeritabilityEnrichment(res)
    expect_s3_class(enr, "data.frame")
    expect_equal(nrow(enr), 2L)
    expect_true(all(
        c(
            "annotation",
            "tau",
            "tauSe",
            "enrichment",
            "enrichmentSe",
            "enrichmentP",
            "propH2",
            "propSnps"
        ) %in%
            names(enr)
    ))
    expect_true(is.matrix(annotationJackknifeCoefs(res)[[1]]))
})

test_that("S-LDSC requires an LdScore, not an LdEigen", {
    eigenRef <- makeTestEigenRef()
    ss <- makeTestSumstatsForRef(eigenRef)
    expect_error(
        estimateH2(
            ss,
            eigenRef,
            param = .h2MethodSettings("sldsc", NULL, eigenRef)
        ),
        "requires an LdScore"
    )
})

test_that("the S-LDSC weights reproduce the original Hsq.weights", {
    # het_w = 1 / (2 (intercept + hsq*N/M * ld)^2), oc_w = 1 / w_ld, both
    # LD scores floored at 1 and hsq clamped to [0, 1].
    w <- pecotmr:::.h2LdScoreWeights(
        ld = c(4, 1, 0.5),
        wLd = c(10, 1, 0.2),
        n = 1000,
        M = 100,
        hsq = 0.5,
        intercept = 1
    )
    cc <- 0.5 * 1000 / 100
    expected <- (1 / (2 * (1 + cc * c(4, 1, 1))^2)) / c(10, 1, 1)
    expect_equal(w, expected)
    expect_true(all(w > 0))

    # hsq is clamped, so an out-of-range estimate cannot blow the weights up.
    expect_equal(
        pecotmr:::.h2LdScoreWeights(2, 1, 1000, 100, hsq = 5, intercept = 1),
        pecotmr:::.h2LdScoreWeights(2, 1, 1000, 100, hsq = 1, intercept = 1)
    )
    expect_equal(
        pecotmr:::.h2LdScoreWeights(2, 1, 1000, 100, hsq = -3, intercept = 1),
        pecotmr:::.h2LdScoreWeights(2, 1, 1000, 100, hsq = 0, intercept = 1)
    )
})

# ===========================================================================
# Rank-deficient LD
#
# A real LD block is routinely rank-deficient, and its decomposition then
# carries a tail of tiny negative eigenvalues (order -1e-15) that are pure
# floating-point noise. Those used to reach the eigen regression's weights as
# negative regression weights, so IRWLS took sqrt() of a negative and lm.fit
# aborted with "NA/NaN/Inf in 'x'".
# ===========================================================================

# An LD block with `dup` variants repeated, so R is singular by construction
# and eigen() returns a negative-noise tail.
.rankDeficientEigenRef <- function(nBlocks = 6L, p = 12L, dup = 4L) {
    blocks <- map(seq_len(nBlocks), function(b) {
        base <- 0.5^abs(outer(seq_len(p - dup), seq_len(p - dup), "-"))
        keep <- c(seq_len(p - dup), seq_len(dup))
        R <- base[keep, keep]
        e <- eigen(R, symmetric = TRUE)
        list(
            values = e$values,
            vectors = e$vectors,
            snpIdx = as.integer(((b - 1L) * p + 1L):(b * p))
        )
    })
    M <- nBlocks * p
    snpInfo <- data.frame(
        SNP = paste0("rs", seq_len(M)),
        CHR = "chr1",
        BP = as.integer(seq_len(M) * 100L),
        A1 = "A",
        A2 = "G",
        stringsAsFactors = FALSE
    )
    starts <- as.integer(map_dbl(blocks, function(b) min(b$snpIdx))) * 100L
    ends <- as.integer(map_dbl(blocks, function(b) max(b$snpIdx))) * 100L
    LdEigen(
        ldBlocks = GenomicRanges::GRanges(
            "chr1",
            IRanges::IRanges(start = starts, end = ends)
        ),
        snpInfo = snpInfo,
        nSamples = 500000L,
        inSample = FALSE,
        genome = "hg19",
        eigenvalueTruncation = 1.0,
        eigenList = blocks
    )
}

test_that("a rank-deficient LD block really does yield negative eigenvalues", {
    ref <- .rankDeficientEigenRef()
    values <- unlist(map(eigenList(ref), function(b) b$values))
    expect_true(any(values < 0))
    expect_true(min(values) > -1e-6) # noise, not structure
})

test_that("estimateH2 survives negative eigenvalues from rank-deficient LD", {
    ref <- .rankDeficientEigenRef()
    ss <- makeTestSumstatsForRef(ref)
    res <- estimateH2(ss, ref)
    expect_s4_class(res, "H2Estimate")
    expect_false(is.na(heritabilityEstimate(res)))
})

test_that("the eigen regression's weights are never negative", {
    lam <- c(2, 1, 0.5, 0, -1e-15, -3e-15)
    # The in-sample (linear) taper, min(1, lambda), is the one a negative
    # eigenvalue would turn negative.
    param <- .h2GenomeWideSettings("lder", FALSE, TRUE)
    rows <- pecotmr:::.h2EigenRegressionRows(
        list(list(x = rep(1, 6), lam = lam, ldAnnot = matrix(1, 6, 1))),
        NULL,
        1e5,
        1000,
        param
    )
    w <- pecotmr:::.h2RegressionWeights(param@procedures, rows, c(0.1, 1))
    expect_false(any(is.na(w)))
    expect_true(all(w >= 0))
    # sqrt() of these is what IRWLS forms the design from.
    expect_false(any(is.na(sqrt(w))))
})

# ===========================================================================
# buildLdEigen / buildLdScore -> estimateH2 round trip
#
# The builders exist so that a loaded LD reference can reach estimateH2
# without hand-assembling S4 slots; these pin that the objects they produce
# are actually accepted, alignment guard included.
# ===========================================================================

test_that("estimateH2 runs on references built by buildLdEigen/buildLdScore", {
    ld <- makeTestLdMatrixBlocks(sizes = rep(10L, 5), rho = 0.6)
    eigenRef <- buildLdEigen(ld, genome = "hg19")
    scoreRef <- buildLdScore(ld, genome = "hg19")
    ss <- makeTestSumstatsForRef(eigenRef)

    expect_s4_class(estimateH2(ss, eigenRef), "H2Estimate")
    suppressWarnings(
        expect_s4_class(
            estimateH2(
                ss,
                eigenRef,
                param = .h2MethodSettings("hdl", NULL, eigenRef)
            ),
            "H2Estimate"
        )
    )
    expect_s4_class(
        estimateH2(
            ss,
            scoreRef,
            param = .h2MethodSettings("gldsc", NULL, scoreRef)
        ),
        "H2Estimate"
    )
})

test_that("estimateH2 rejects a single-block LD reference", {
    # One block used to fail inside the jackknife fold fit with "'x' must be
    # a matrix" (LDER, g-LDSC) or return an SE of exactly 0 (HDL).
    ld <- makeTestLdMatrix(n = 30L, rho = 0.6)
    eigenRef <- buildLdEigen(ld, genome = "hg19")
    ss <- makeTestSumstatsForRef(eigenRef)
    expect_error(
        estimateH2(ss, eigenRef),
        "delete-one-block jackknife"
    )
    expect_error(
        estimateH2(
            ss,
            buildLdScore(ld, genome = "hg19"),
            param = .h2GenomeWideSettings("gldsc", FALSE, FALSE)
        ),
        "at least two are required"
    )
})

test_that("buildLdScore(keepLdMatrices = FALSE) is rejected by g-LDSC", {
    ld <- makeTestLdMatrixBlocks(sizes = rep(10L, 5))
    scoreRef <- buildLdScore(ld, genome = "hg19", keepLdMatrices = FALSE)
    ss <- makeTestSumstatsForRef(scoreRef)
    expect_error(
        estimateH2(
            ss,
            scoreRef,
            param = .h2MethodSettings("gldsc", NULL, scoreRef)
        ),
        "The variant GLS needs each LD block's LD"
    )
})

test_that("estimateH2 accepts a chr-prefix difference between the two sides", {
    ref <- makeTestEigenRef()
    ss <- makeTestSumstatsForRef(ref)
    # makeTestSumstatsForRef builds the entry on bare names ("1" vs the
    # reference's "chr1"); the constructor reads them in UCSC style, so the
    # two sides meet exactly.
    expect_equal(
        as.character(GenomicRanges::seqnames(sumStats(ss)))[[1]],
        "chr1"
    )
    expect_s4_class(estimateH2(ss, ref), "H2Estimate")
})

test_that("a study's recorded varY leaves its genome-wide h2 as it is", {
    eigenRef <- makeTestEigenRef()
    plain <- estimateH2(makeTestSumstatsForRef(eigenRef), eigenRef)
    withVarY <- makeTestSumstatsForRef(eigenRef, varY = 4)
    expect_equal(varY(withVarY), 4)
    expect_identical(
        heritabilityEstimate(estimateH2(withVarY, eigenRef)),
        heritabilityEstimate(plain)
    )
    # A study that records no prevalence has no liability-scale estimate.
    expect_true(is.na(S4Vectors::mcols(plain)$h2_liability))
})

test_that("a case/control study's prevalence gives h2 on the liability scale", {
    eigenRef <- makeTestEigenRef()
    # Lee et al. (2011) equation 23, as ldsc's h2_obs_to_liab() and LAVA's
    # h2.latent read it: K^2 (1 - K)^2 / (P (1 - P) phi(t)^2).
    lee <- function(K, P) {
        K^2 * (1 - K)^2 / (P * (1 - P) * dnorm(qnorm(1 - K))^2)
    }
    # 10,000 cases of the 50,000 samples N counts: P = 0.2.
    ss <- makeTestSumstatsForRef(
        eigenRef,
        nCase = 10000,
        nControl = 40000,
        prevalence = 0.01
    )
    expect_equal(pecotmr:::.h2LiabilityFactor(ss, "test"), lee(0.01, 0.2))
    # Where the QC set N to the effective sample size, that size counts a
    # balanced study's samples: summaryStatsQc()'s record, one audit per row,
    # which qcInfo() reads as nSource.
    qcd <- `$<-`(
        ss,
        "qcAudit",
        S4Vectors::SimpleList(list(list(nSource = "effective")))
    )
    expect_identical(qcInfo(qcd)$nSource, "effective")
    expect_equal(pecotmr:::.h2LiabilityFactor(qcd, "test"), lee(0.01, 0.5))
    # Per-variant counts give the share where the study records them so.
    gr <- sumStats(ss, studyName = "test")
    S4Vectors::mcols(gr)$N_CASE <- rep(5000, length(gr))
    S4Vectors::mcols(gr)$N_CONTROL <- rep(45000, length(gr))
    perVariant <- GwasSumStats(
        studyName = "test",
        entry = list(gr),
        genome = "hg19",
        ldReference = makeTestGwasGenotypeHandle(),
        prevalence = 0.01
    )
    expect_equal(
        pecotmr:::.h2LiabilityFactor(perVariant, "test"),
        lee(0.01, 0.1)
    )
    expect_true(is.na(pecotmr:::.h2LiabilityFactor(
        makeTestSumstatsForRef(eigenRef),
        "test"
    )))
    # The estimate holds its heritability and standard error times the
    # factor; the fit here gives a known one.
    local_mocked_bindings(
        .h2GenomeWideFit = function(...) {
            list(h2 = 0.3, h2Se = 0.05, intercept = 1.02, interceptSe = 0.01)
        }
    )
    est <- S4Vectors::mcols(estimateH2(ss, eigenRef))
    expect_equal(est$h2, 0.3)
    expect_equal(est$h2_liability, 0.3 * lee(0.01, 0.2))
    expect_equal(est$h2_liability_se, 0.05 * lee(0.01, 0.2))
    expect_true(is.na(est$h2_liability_lower) && is.na(est$h2_liability_upper))
    eff <- S4Vectors::mcols(estimateH2(qcd, eigenRef))
    expect_equal(eff$h2_liability, 0.3 * lee(0.01, 0.5))
})

test_that("a study's rows that disagree on its N's source are refused", {
    eigenRef <- makeTestEigenRef()
    # The study in two rows, one per LD block.
    ss <- makeTestSumstatsForRef(
        eigenRef,
        nCase = 10000,
        nControl = 40000,
        prevalence = 0.01,
        region = GenomicRanges::GRanges(
            "chr1",
            IRanges::IRanges(c(1, 1001), width = 999)
        )
    )
    expect_identical(nrow(ss), 2L)
    mixed <- `$<-`(
        ss,
        "qcAudit",
        S4Vectors::SimpleList(
            list(nSource = "effective"),
            list(nSource = "total")
        )
    )
    expect_error(
        pecotmr:::.h2StudyNSource(mixed, "test"),
        "disagree on how the QC set N"
    )
    expect_identical(pecotmr:::.h2StudyNSource(ss, "test"), NA_character_)
    # Rows disagreeing on a per-study value leave it no one value either.
    twoK <- ss
    S4Vectors::mcols(twoK)$prevalence <- c(0.01, 0.02)
    expect_error(
        pecotmr:::.h2StudyValue(twoK, "prevalence", "test"),
        "disagree on prevalence"
    )
})

# ===========================================================================
# computeLdScores for LdScore WITH annotations
# ===========================================================================

test_that("computeLdScores(LdScore) with annotations and matrices is right", {
    score_ref <- makeTestScoreRef(with_ld_matrices = TRUE)
    annot <- makeTestAnnotations()
    result <- computeLdScores(score_ref, annotations = annot)

    expect_true(is.matrix(result))
    # base_l2 + 2 annotation columns
    expect_equal(ncol(result), 3)
    expect_equal(nrow(result), length(score_ref))
    expect_equal(colnames(result), c("base_l2", "annot_A", "annot_B"))
    # First column should be the stored base LD scores
    expect_equal(result[, 1], mcols(score_ref)$ldScores[, 1])
})


context("h2_utils")

# =============================================================================
# .jackknifeSeFromLoo
# =============================================================================

test_that(".jackknifeSeFromLoo is zero when every LOO estimate is equal", {
    expect_equal(pecotmr:::.jackknifeSeFromLoo(rep(2, 5)), 0)
})

test_that(".jackknifeSeFromLoo equals the pseudo-value jackknife", {
    # The delete-one form sqrt((n - 1) / n * sum((loo - mean)^2)) and the
    # pseudo-value form sqrt(var(n * full - (n - 1) * loo) / n) are the same
    # quantity -- the full estimate cancels -- which is why one helper serves
    # every estimator. Leave-one-out means of a known vector:
    vals <- c(2, 4, 6, 8, 10)
    n <- length(vals)
    looMeans <- map_dbl(seq_len(n), function(i) mean(vals[-i]))
    pseudo <- n * mean(vals) - (n - 1) * looMeans
    expect_equal(
        pecotmr:::.jackknifeSeFromLoo(looMeans),
        sqrt(var(pseudo) / n)
    )
})

# =============================================================================
# standardizeTauStar
# =============================================================================

test_that("standardizeTauStar computes tauStar = tau * sdAnnot * MRef / h2g", {
    tau <- c(0.01, 0.02)
    sd_annot <- c(0.5, 0.3)
    M_ref <- 1000L
    h2g <- 0.4
    nBlocks <- 5
    annotationJackknifeCoefs <- matrix(
        rep(tau, each = nBlocks),
        nrow = nBlocks,
        ncol = 2
    )
    res <- pecotmr:::standardizeTauStar(
        tau,
        annotationJackknifeCoefs,
        sd_annot,
        M_ref,
        h2g
    )
    expected <- tau * sd_annot * M_ref / h2g
    expect_true(is.list(res))
    expect_named(res, c("tauStar", "tauStarSe"))
    expect_equal(res$tauStar, expected)
    expect_length(res$tauStarSe, 2)
})

test_that("standardizeTauStar errors when h2g == 0", {
    expect_error(
        pecotmr:::standardizeTauStar(
            c(0.01),
            matrix(0.01, nrow = 3, ncol = 1),
            c(0.5),
            1000L,
            0
        ),
        "h2g must be non-zero"
    )
})

test_that("standardizeTauStar errors when tau and sdAnnot differ in length", {
    expect_error(
        pecotmr:::standardizeTauStar(
            c(0.01, 0.02),
            matrix(0, nrow = 3, ncol = 2),
            c(0.5),
            1000L,
            0.5
        ),
        "tau and sdAnnot must have the same length"
    )
})

# =============================================================================
# .rmaMeta
# =============================================================================

test_that(".rmaMeta returns all NA with k=0", {
    res <- pecotmr:::.rmaMeta(numeric(0), numeric(0))
    expect_true(is.na(res$mean))
    expect_true(is.na(res$se))
    expect_true(is.na(res$tau2))
    expect_true(is.na(res$I2))
    expect_true(is.na(res$Q))
})

test_that(".rmaMeta with k=1 returns input values", {
    res <- pecotmr:::.rmaMeta(5.0, 1.0)
    expect_equal(res$mean, 5.0)
    expect_equal(res$se, 1.0)
    expect_equal(res$tau2, 0)
})

test_that(".rmaMeta with identical means gives tau2=0", {
    means <- rep(3.0, 5)
    ses <- rep(1.0, 5)
    res <- pecotmr:::.rmaMeta(means, ses)
    expect_equal(res$tau2, 0)
    expect_equal(res$mean, 3.0)
})

test_that(".rmaMeta errors with non-positive ses", {
    expect_error(
        pecotmr:::.rmaMeta(c(1, 2), c(1, 0)),
        "all ses must be positive and finite"
    )
    expect_error(
        pecotmr:::.rmaMeta(c(1, 2), c(1, -1)),
        "all ses must be positive and finite"
    )
})

test_that(".rmaMeta defaults to DL and matches a known example", {
    # Three studies with known values
    means <- c(0.5, 0.8, 0.3)
    ses <- c(0.2, 0.3, 0.15)
    res <- pecotmr:::.rmaMeta(means, ses)
    expect_equal(res, pecotmr:::.rmaMeta(means, ses, method = "DL"))
    expect_true(is.list(res))
    expect_named(res, c("mean", "se", "tau2", "I2", "Q"))
    # metafor reports I2 as a percentage; .rmaMeta rescales it to [0, 1].
    expect_true(res$I2 >= 0 && res$I2 <= 1)

    # Fixed-effect weights
    w_fe <- 1 / ses^2
    mu_fe <- sum(w_fe * means) / sum(w_fe)
    Q <- sum(w_fe * (means - mu_fe)^2)
    c_dl <- sum(w_fe) - sum(w_fe^2) / sum(w_fe)
    tau2 <- max(0, (Q - 2) / c_dl)
    w_re <- 1 / (ses^2 + tau2)
    mu_re <- sum(w_re * means) / sum(w_re)
    se_re <- sqrt(1 / sum(w_re))

    expect_equal(res$Q, Q)
    expect_equal(res$tau2, tau2)
    expect_equal(res$mean, mu_re)
    expect_equal(res$se, se_re)
})

test_that(".rmaMeta runs REML through metafor, differing from DL", {
    set.seed(3)
    m <- rnorm(8, sd = 1.5)
    s <- abs(rnorm(8)) + 0.3
    reml <- pecotmr:::.rmaMeta(m, s, method = "REML")
    expect_true(all(c("mean", "se", "tau2", "I2", "Q") %in% names(reml)))
    expect_true(reml$I2 >= 0 && reml$I2 <= 1)
    # REML tau2 generally differs from DL tau2 on heterogeneous data
    expect_false(isTRUE(all.equal(reml$tau2, pecotmr:::.rmaMeta(m, s)$tau2)))
})

test_that(".rmaMeta falls back to DL when REML does not converge", {
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

# =============================================================================
# snpsPerBlock
# =============================================================================

test_that("snpsPerBlock assigns SNPs to correct blocks", {
    # 10 SNPs across 2 blocks
    snp_info <- data.frame(
        CHR = rep("chr1", 10),
        BP = c(100, 500, 1000, 3000, 4500, 5500, 6000, 7500, 9000, 9999)
    )
    blocks_gr <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1"),
        ranges = IRanges::IRanges(start = c(1, 5001), end = c(5000, 10000))
    )
    ld_blocks <- blocks_gr
    res <- pecotmr:::snpsPerBlock(snp_info, ld_blocks)
    # Block 1 covers 1-5000: SNPs at 100, 500, 1000, 3000, 4500
    expect_equal(sort(res[["1"]]), c(1L, 2L, 3L, 4L, 5L))
    # Block 2 covers 5001-10000: SNPs at 5500, 6000, 7500, 9000, 9999
    expect_equal(sort(res[["2"]]), c(6L, 7L, 8L, 9L, 10L))
})

test_that("snpsPerBlock returns empty for blocks with no SNPs", {
    snp_info <- data.frame(
        CHR = rep("chr1", 3),
        BP = c(100, 200, 300)
    )
    blocks_gr <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr2"),
        ranges = IRanges::IRanges(start = c(1, 1), end = c(5000, 5000))
    )
    ld_blocks <- blocks_gr
    res <- pecotmr:::snpsPerBlock(snp_info, ld_blocks)
    # All SNPs on chr1, so block 2 (chr2) should have no entries
    expect_true("1" %in% names(res))
    expect_false("2" %in% names(res))
    expect_equal(sort(res[["1"]]), c(1L, 2L, 3L))
})

# =============================================================================
# checkGenomeBuild
# =============================================================================

test_that("checkGenomeBuild returns TRUE when all objects match", {
    blocks_gr_a <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 1, end = 5000)
    )
    ld_blocks_a <- blocks_gr_a
    GenomeInfoDb::genome(ld_blocks_a) <- "hg19"
    ld_blocks_b <- ld_blocks_a

    expect_true(pecotmr:::checkGenomeBuild(ld_blocks_a, ld_blocks_b))
})

test_that("checkGenomeBuild errors when genome builds mismatch", {
    blocks_gr_19 <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 1, end = 5000)
    )
    ld_blocks_19 <- blocks_gr_19
    GenomeInfoDb::genome(ld_blocks_19) <- "hg19"

    blocks_gr_38 <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 1, end = 5000)
    )
    ld_blocks_38 <- blocks_gr_38
    GenomeInfoDb::genome(ld_blocks_38) <- "hg38"

    expect_error(
        pecotmr:::checkGenomeBuild(ld_blocks_19, ld_blocks_38),
        "Genome build mismatch"
    )
})

test_that("checkGenomeBuild works with AnnotationMatrix", {
    snp_gr <- GenomicRanges::GRanges(
        seqnames = rep("chr1", 3),
        ranges = IRanges::IRanges(start = 1:3, width = 1)
    )
    annot_mat <- matrix(c(1, 0, 1), nrow = 3, ncol = 1)
    annot_meta <- data.frame(
        name = "annot1",
        tier = "baseline",
        type = "binary",
        stringsAsFactors = FALSE
    )
    am <- AnnotationMatrix(
        snpRanges = snp_gr,
        annotations = annot_mat,
        annotationMeta = annot_meta,
        genome = "hg19"
    )

    blocks_gr <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 1, end = 5000)
    )
    ld_blocks <- blocks_gr
    GenomeInfoDb::genome(ld_blocks) <- "hg19"

    expect_true(pecotmr:::checkGenomeBuild(am, ld_blocks))
})

# Tests for the heritability estimators eigenRegressionUnivariate,
# variantGlsUnivariate and eigenLikelihoodUnivariate.

set.seed(42)

# =============================================================================
# Helpers: simulate test data
# =============================================================================

simulateH2Data <- function(
    n_snps = 100,
    nBlocks = 2,
    n_gwas = 50000,
    h2_true = 0.3
) {
    set.seed(42)
    snps_per_block <- n_snps / nBlocks

    # Block structure
    blocks_gr <- GenomicRanges::GRanges(
        seqnames = rep("chr1", nBlocks),
        ranges = IRanges::IRanges(
            start = seq(1, by = snps_per_block * 100, length.out = nBlocks),
            width = snps_per_block * 100 - 1
        )
    )
    ld_blocks <- blocks_gr

    # SNP info
    snp_info <- data.frame(
        SNP = paste0("rs", seq_len(n_snps)),
        CHR = rep("chr1", n_snps),
        BP = as.integer(seq(50, by = 100, length.out = n_snps)),
        A1 = rep("A", n_snps),
        A2 = rep("G", n_snps),
        stringsAsFactors = FALSE
    )

    # Build block-diagonal LD (AR(1) with rho = 0.5)
    eigen_list <- list()
    R_blocks <- list()
    for (b in seq_len(nBlocks)) {
        idx <- seq((b - 1) * snps_per_block + 1, b * snps_per_block)
        p <- length(idx)
        R <- 0.5^abs(outer(seq_len(p), seq_len(p), "-"))
        e <- eigen(R, symmetric = TRUE)
        eigen_list[[b]] <- list(
            values = e$values,
            vectors = e$vectors,
            snpIdx = as.integer(idx)
        )
        R_blocks[[b]] <- list(R = R, snpIdx = as.integer(idx))
    }

    # Simulate z-scores under the LD-aware infinitesimal model: per block
    # Cov(z) = (N h2 / M) R R + R (marginal-effect LD covariance + sampling
    # covariance). White-noise z (ignoring LD) is inconsistent with the LD
    # reference and makes every eigenvalue-based estimator peg at a boundary.
    z <- numeric(n_snps)
    for (b in seq_len(nBlocks)) {
        idx <- R_blocks[[b]]$snpIdx
        R <- R_blocks[[b]]$R
        sigma <- (n_gwas * h2_true / n_snps) * (R %*% R) + R
        z[idx] <- as.vector(crossprod(chol(sigma), rnorm(length(idx))))
    }

    # LdEigen
    eigen_ref <- LdEigen(
        ldBlocks = ld_blocks,
        snpInfo = snp_info,
        nSamples = 500000L,
        inSample = FALSE,
        genome = "hg19",
        eigenList = eigen_list,
        eigenvalueTruncation = 1.0
    )

    # LD scores: l2_j = sum_k r^2_{jk}
    ld_scores_vec <- numeric(n_snps)
    for (b in seq_len(nBlocks)) {
        idx <- R_blocks[[b]]$snpIdx
        R <- R_blocks[[b]]$R
        ld_scores_vec[idx] <- rowSums(R^2)
    }

    ld_score_ref <- LdScore(
        ldBlocks = ld_blocks,
        snpInfo = snp_info,
        nSamples = 500000L,
        inSample = FALSE,
        genome = "hg19",
        ldScores = matrix(
            ld_scores_vec,
            ncol = 1,
            dimnames = list(NULL, "base_l2")
        ),
        ldScoreWeights = rep(1, n_snps),
        ldMatrixList = R_blocks
    )

    list(
        z = z,
        n = n_gwas,
        eigen_ref = eigen_ref,
        ld_score_ref = ld_score_ref,
        h2_true = h2_true,
        n_snps = n_snps
    )
}

# Flexible annotation builder: one baseline column (all-ones or spatial
# half-on/half-off) plus `nCand` random binary candidate columns.
makeCoverageAnnot <- function(
    n_snps,
    baseline = "ones",
    nCand = 1L,
    seed = 321
) {
    set.seed(seed)
    snp_gr <- GenomicRanges::GRanges(
        seqnames = rep("chr1", n_snps),
        ranges = IRanges::IRanges(
            start = seq(50, by = 100, length.out = n_snps),
            width = 1L
        )
    )
    base_col <- if (identical(baseline, "spatial")) {
        as.numeric(seq_len(n_snps) <= n_snps / 2)
    } else {
        rep(1, n_snps)
    }
    mat <- matrix(base_col, ncol = 1)
    meta <- data.frame(
        name = "baseline1",
        tier = "baseline",
        type = "binary",
        stringsAsFactors = FALSE
    )
    if (nCand > 0L) {
        candidateCols <- map(
            seq_len(nCand),
            function(k) as.numeric(rbinom(n_snps, 1, 0.3 + 0.1 * k))
        )
        cand <- exec(cbind, !!!candidateCols)
        mat <- cbind(mat, cand)
        meta <- rbind(
            meta,
            data.frame(
                name = paste0("cand", seq_len(nCand)),
                tier = "candidate",
                type = "binary",
                stringsAsFactors = FALSE
            )
        )
    }
    colnames(mat) <- meta$name
    AnnotationMatrix(mat, snp_gr, meta, genome = "hg19")
}

test_that(".h2Stratified reads the width of the baseline design", {
    expect_false(.h2Stratified(NULL))
    # An all-ones baseline alone is the unstratified design.
    expect_false(.h2Stratified(makeCoverageAnnot(100, "ones", nCand = 0L)))
    expect_false(.h2Stratified(makeCoverageAnnot(100, "ones", nCand = 2L)))
    # One annotation that does not span the genome gets a base column too.
    expect_true(.h2Stratified(makeCoverageAnnot(100, "spatial", nCand = 0L)))
})

# The four estimators, each with the reference it reads and the baseline
# its annotations take. An all-ones baseline makes the LDER and HDL
# eigenvalue-score column the eigenvalue itself; for g-LDSC and S-LDSC it
# would repeat the base LD score column and leave the design singular, so
# they take a half-on/half-off baseline instead.
.h2Runs <- list(
    lder = list(
        fn = pecotmr:::eigenRegressionUnivariate,
        ref = "eigen_ref",
        baseline = "ones"
    ),
    gldsc = list(
        fn = pecotmr:::variantGlsUnivariate,
        ref = "ld_score_ref",
        baseline = "spatial"
    ),
    hdl = list(
        fn = pecotmr:::eigenLikelihoodUnivariate,
        ref = "eigen_ref",
        baseline = "ones"
    ),
    sldsc = list(
        fn = pecotmr:::variantRegressionUnivariate,
        ref = "ld_score_ref",
        baseline = "spatial"
    )
)

# Run `method` on the simulated data `d`, passing `...` to the estimator.
.h2Run <- function(method, d, ...) {
    run <- .h2Runs[[method]]
    run$fn(d$z, d$n, d[[run$ref]], ...)
}

# Annotations over `d`'s variants for `method`, with `nCand` candidates.
.h2RunAnnot <- function(method, d, nCand) {
    makeCoverageAnnot(
        d$n_snps,
        baseline = .h2Runs[[method]]$baseline,
        nCand = nCand
    )
}

# Pre-compute shared test data
dat <- simulateH2Data()
# Five blocks and 500 SNPs for the annotation tests, so that the stratified
# designs are not singular.
dat_annot <- simulateH2Data(n_snps = 500, nBlocks = 5)

# =============================================================================
# Unstratified estimates
# =============================================================================

test_that("each univariate estimator returns a finite h2 and a positive SE", {
    fields <- c(
        "h2",
        "h2Se",
        "intercept",
        "interceptSe",
        "enrichment",
        "annotationJackknifeCoefs",
        "scoreStats"
    )
    for (method in names(.h2Runs)) {
        res <- .h2Run(method, dat)
        expect_true(is.list(res), info = method)
        expect_true(all(fields %in% names(res)), info = method)
        expect_true(is.finite(res$h2), info = method)
        expect_true(res$h2 > -0.5 && res$h2 < 1.0, info = method)
        expect_true(is.finite(res$h2Se), info = method)
        expect_true(res$h2Se > 0, info = method)
        # Without annotations there is no enrichment to report.
        expect_null(res$enrichment, info = method)
    }
})

test_that("the eigen regression reports its intercept on the others' scale", {
    # The fit's inflation a enters as E[chi2Rot] = n * h2/M * d + 1 + n * a. The
    # estimators report that intercept, 1 + n * a, about 1 without confounding.
    res <- pecotmr:::eigenRegressionUnivariate(dat$z, dat$n, dat$eigen_ref)
    expect_true(is.finite(res$intercept))
    expect_gte(res$intercept, 1)
    expect_lt(res$intercept, 1 + dat$n * 0.01)
    expect_gt(res$interceptSe, 0)
})

# =============================================================================
# Stratified estimates
# =============================================================================

test_that("each univariate estimator reports enrichment and score statistics", {
    expected_cols <- c(
        "annotation",
        "tau",
        "tauSe",
        "enrichment",
        "enrichmentSe",
        "enrichmentP",
        "propH2",
        "propSnps"
    )
    for (method in names(.h2Runs)) {
        annot <- .h2RunAnnot(method, dat_annot, 1L)
        res <- .h2Run(method, dat_annot, annotations = annot)
        # One enrichment row per baseline annotation
        expect_true(is.data.frame(res$enrichment), info = method)
        expect_true(
            all(expected_cols %in% colnames(res$enrichment)),
            info = method
        )
        expect_equal(nrow(res$enrichment), 1L, info = method)
        expect_true(all(is.finite(res$enrichment$tau)), info = method)
        expect_true(all(is.finite(res$enrichment$enrichment)), info = method)
        # One jackknife row per LD block (5), one column per baseline
        # annotation
        jk <- res$annotationJackknifeCoefs
        expect_true(is.matrix(jk), info = method)
        expect_equal(dim(jk), c(5L, 1L), info = method)
        expect_true(all(is.finite(jk)), info = method)
        # One score z per candidate annotation
        ss <- res$scoreStats
        expect_true(is.list(ss), info = method)
        expect_true(all(c("z", "R") %in% names(ss)), info = method)
        expect_equal(length(ss$z), 1L, info = method)
        expect_true(is.finite(ss$z[1]), info = method)
        expect_true(is.matrix(ss$R), info = method)
        expect_equal(dim(ss$R), c(1L, 1L), info = method)
        expect_equal(ss$annotationNames, "cand1", info = method)
    }
})

# ---------------------------------------------------------------------------
# Stratified score helpers — multiple-candidate (cor) branch
# ---------------------------------------------------------------------------

test_that("each estimator with 2 candidates returns a 2x2 score correlation", {
    for (method in names(.h2Runs)) {
        annot <- .h2RunAnnot(method, dat_annot, 2L)
        res <- .h2Run(method, dat_annot, annotations = annot)
        expect_equal(length(res$scoreStats$z), 2L, info = method)
        expect_equal(dim(res$scoreStats$R), c(2L, 2L), info = method)
    }
})

# ---------------------------------------------------------------------------
# Stratified score helpers — no-candidate (NULL scoreStats) branch
# ---------------------------------------------------------------------------

test_that("each estimator with baseline-only annotations has NULL scoreStats", {
    for (method in names(.h2Runs)) {
        annot <- .h2RunAnnot(method, dat_annot, 0L)
        res <- .h2Run(method, dat_annot, annotations = annot)
        expect_null(res$scoreStats, info = method)
    }
})

# ===========================================================================
# Error branches and degenerate inputs: an unknown object type for the
# genome check, mismatched .rmaMeta lengths, and unnamed jackknife
# coefficients for h2EstimateToSldscTrait.
# ===========================================================================

# ---------------------------------------------------------------------------
# checkGenomeBuild — unknown object type branch
# ---------------------------------------------------------------------------

test_that("checkGenomeBuild errors on an unrecognized object type", {
    expect_error(
        pecotmr:::checkGenomeBuild(list(genome = "hg19")),
        "Unknown object type"
    )
})

# ---------------------------------------------------------------------------
# .rmaMeta — length-mismatch error
# ---------------------------------------------------------------------------

test_that(".rmaMeta errors when means and ses differ in length", {
    expect_error(
        pecotmr:::.rmaMeta(c(1, 2), c(1, 2, 3)),
        "must have the same length"
    )
})

# ---------------------------------------------------------------------------
# h2EstimateToSldscTrait — annotationJackknifeCoefs present but without column
# names
# ---------------------------------------------------------------------------

test_that("h2EstimateToSldscTrait names unnamed jackknife categories", {
    h2_obj <- makeTestH2Estimate(with_enrichment = TRUE)
    tb <- annotationJackknifeCoefs(h2_obj)[[1]]
    colnames(tb) <- NULL
    h2_obj$jackknife <- list(tb)
    result <- h2EstimateToSldscTrait(h2_obj)
    expect_equal(
        colnames(result$annotationJackknifeCoefs),
        c("annot1", "annot2")
    )
    expect_equal(result$nBlocks, 10L)
})

# ===========================================================================
# Deterministic h2 recovery + regression tests (correctness validation)
#
# These validate the LDER / gLDSC / HDL estimators against data simulated with
# a KNOWN heritability, using fixed seeds so every run is deterministic and
# self-contained (no dependence on the original LDER/gLDSC/HDL packages).
# Recovery tests check the mean estimate recovers the true h2; the stored means
# pin the known-input -> known-output values. A behavioral test locks in the
# LDER inSample wiring.
# ===========================================================================

# AR(1) LD blocks + a chol-based per-block MVN simulator for the standard
# model z ~ N(0, N R diag(perSnpVar) R + R). Returns the eigen and score
# references plus a simulator closure. No external dependencies.
.h2RecoveryRefs <- function(nBlocks = 30L, p = 40L, N = 80000) {
    M <- nBlocks * p
    blocks <- map(seq_len(nBlocks), function(b) {
        R <- 0.6^abs(outer(seq_len(p), seq_len(p), "-"))
        e <- eigen(R, symmetric = TRUE)
        list(
            R = R,
            values = e$values,
            vectors = e$vectors,
            idx = ((b - 1L) * p + 1L):(b * p)
        )
    })
    snpInfo <- data.frame(
        SNP = paste0("rs", seq_len(M)),
        CHR = "chr1",
        BP = as.integer(seq_len(M) * 100L),
        A1 = "A",
        A2 = "G",
        stringsAsFactors = FALSE
    )
    starts <- as.integer(map_dbl(blocks, function(b) min(b$idx))) * 100L
    ends <- as.integer(map_dbl(blocks, function(b) max(b$idx))) * 100L
    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(start = starts, end = ends)
    )
    lb <- gr
    eigenRef <- LdEigen(
        ldBlocks = lb,
        snpInfo = snpInfo,
        nSamples = 500000L,
        inSample = FALSE,
        genome = "hg19",
        eigenvalueTruncation = 1.0,
        eigenList = map(blocks, function(b) {
            list(
                values = b$values,
                vectors = b$vectors,
                snpIdx = as.integer(b$idx)
            )
        })
    )
    lsv <- numeric(M)
    for (b in blocks) {
        lsv[b$idx] <- rowSums(b$R^2)
    }
    scoreRef <- LdScore(
        ldBlocks = lb,
        snpInfo = snpInfo,
        nSamples = 500000L,
        inSample = FALSE,
        genome = "hg19",
        ldScores = matrix(lsv, ncol = 1, dimnames = list(NULL, "base_l2")),
        ldScoreWeights = rep(1, M),
        ldMatrixList = map(blocks, function(b) {
            list(R = b$R, snpIdx = as.integer(b$idx))
        })
    )
    simZ <- function(perSnpVar) {
        unlist(map(blocks, function(b) {
            S <- N * (b$R %*% diag(perSnpVar[b$idx], p) %*% b$R) + b$R
            as.vector(crossprod(chol(S), rnorm(p)))
        }))
    }
    list(eigenRef = eigenRef, scoreRef = scoreRef, simZ = simZ, M = M, N = N)
}

.h2RecoveryAnnot <- function(M, split) {
    annotA <- as.integer(seq_len(M) <= split)
    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(start = as.integer(seq_len(M) * 100L), width = 1L)
    )
    mat <- cbind(annotA = annotA, annotB = 1L - annotA)
    meta <- data.frame(
        name = c("annotA", "annotB"),
        tier = c("baseline", "baseline"),
        type = "binary",
        stringsAsFactors = FALSE
    )
    AnnotationMatrix(mat, gr, meta, genome = "hg19")
}

# The HDL recovery tests use 10 blocks rather than 30: HDL refits once per
# block for its jackknife, so its cost grows with the square of the block
# count, and on 10 blocks the mean estimate still lands well inside the
# tolerances.
test_that("eigenLikelihoodUnivariate recovers a known h2 (seed = 1)", {
    d <- .h2RecoveryRefs(nBlocks = 10L)
    uni <- rep(0.4 / d$M, d$M)
    set.seed(1)
    ests <- replicate(
        30,
        eigenLikelihoodUnivariate(d$simZ(uni), d$N, d$eigenRef)$h2
    )
    expect_lt(abs(mean(ests) - 0.4), 0.05) # recovers true h2 = 0.4
    expect_equal(mean(ests), 0.395, tolerance = 0.04) # stored known output
})

test_that("eigenRegressionUnivariate recovers a known h2 (seed = 1)", {
    d <- .h2RecoveryRefs()
    uni <- rep(0.4 / d$M, d$M)
    set.seed(1)
    ests <- replicate(
        30,
        eigenRegressionUnivariate(d$simZ(uni), d$N, d$eigenRef)$h2
    )
    expect_lt(abs(mean(ests) - 0.4), 0.05)
    expect_equal(mean(ests), 0.389, tolerance = 0.04)
})

# ===========================================================================
# Candidate-annotation score tests: calibration and power
#
# A candidate's score is calibrated by the spread of its LD blocks' scores, so
# under the null its z is close to N(0, 1) whatever the regression weights
# get wrong. Both LD-score estimators are checked through
# estimateStratifiedH2().
# ===========================================================================

# A half-on/half-off baseline annotation plus the candidate `cand`.
.h2ScoreAnnot <- function(M, cand) {
    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(start = as.integer(seq_len(M) * 100L), width = 1L)
    )
    mat <- cbind(half = as.integer(seq_len(M) <= M / 2), cand = cand)
    meta <- data.frame(
        name = c("half", "cand"),
        tier = c("baseline", "candidate"),
        type = "binary"
    )
    AnnotationMatrix(mat, gr, meta, genome = "hg19")
}

# Simulated z-scores over the recovery reference as a GwasSumStats.
.h2ScoreSumstats <- function(d, z) {
    df <- data.frame(
        SNP = paste0("rs", seq_len(d$M)),
        CHR = "chr1",
        BP = as.integer(seq_len(d$M) * 100L),
        A1 = "A",
        A2 = "G",
        Z = z,
        N = d$N
    )
    GwasSumStats(
        studyName = "sim",
        entry = list(.dfToGwasGr(df)),
        genome = "hg19",
        ldReference = makeTestGwasGenotypeHandle()
    )
}

# The reference each method reads.
.h2MethodRef <- function(d, method) {
    if (method %in% c("lder", "hdl")) d$eigenRef else d$scoreRef
}

# Score z of a random 30% candidate whose variants carry `enrich` times the
# per-variant heritability of the rest, for each method, over `nRep` traits.
.h2CandidateZ <- function(
    d,
    enrich,
    nRep,
    methods = c("sldsc", "gldsc", "lder", "hdl")
) {
    t(replicate(nRep, {
        cand <- stats::rbinom(d$M, 1, 0.3)
        perSnp <- 0.3 / d$M * ifelse(cand == 1L, enrich, 1)
        ss <- .h2ScoreSumstats(d, d$simZ(perSnp))
        annot <- .h2ScoreAnnot(d$M, cand)
        vapply(
            methods,
            function(m) {
                ref <- .h2MethodRef(d, m)
                scoreStats(suppressWarnings(estimateStratifiedH2(
                    ss,
                    ref,
                    annot,
                    param = .h2MethodSettings(m, annot, ref)
                )))[[1L]]$z
            },
            numeric(1)
        )
    }))
}

test_that("candidate score z is calibrated under the null", {
    d <- .h2RecoveryRefs(nBlocks = 20L)
    set.seed(11)
    z <- .h2CandidateZ(d, enrich = 1, nRep = 60L)
    for (method in colnames(z)) {
        expect_lt(abs(mean(z[, method])), 0.25)
        expect_gt(stats::sd(z[, method]), 0.7)
        expect_lt(stats::sd(z[, method]), 1.35)
    }
})

test_that("candidate score z detects an enriched candidate", {
    d <- .h2RecoveryRefs(nBlocks = 20L)
    set.seed(12)
    z <- .h2CandidateZ(d, enrich = 4, nRep = 15L)
    expect_gt(mean(z[, "gldsc"]), 3)
    # LDER and HDL read the same eigen-decomposition, and on this reference
    # detect it less sharply than g-LDSC's GLS does.
    expect_gt(mean(z[, "lder"]), 1.5)
    expect_gt(mean(z[, "hdl"]), 1.5)
    # S-LDSC spends its weight on the base score and so tests less sharply
    # on a reference this small, but still points the right way.
    expect_gt(mean(z[, "sldsc"]), 0.5)
})

# ===========================================================================
# The baseline design: an all-ones base column unless the baseline spans it
# ===========================================================================

# A baseline of the given columns over the recovery reference.
.h2BaselineOf <- function(d, cols) {
    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(start = as.integer(seq_len(d$M) * 100L), width = 1L)
    )
    AnnotationMatrix(
        cols,
        gr,
        data.frame(name = colnames(cols), tier = "baseline", type = "binary"),
        genome = "hg19"
    )
}

# Each method's fit, stratified by `annotations`, or unstratified by
# estimateH2() where they are NULL.
.h2FitAll <- function(d, ss, annotations) {
    lapply(
        set_names(c("lder", "hdl", "gldsc", "sldsc")),
        function(m) {
            ref <- .h2MethodRef(d, m)
            suppressWarnings(
                if (is.null(annotations)) {
                    estimateH2(
                        ss,
                        ref,
                        param = .h2MethodSettings(m, NULL, ref)
                    )
                } else {
                    estimateStratifiedH2(
                        ss,
                        ref,
                        annotations,
                        param = .h2MethodSettings(m, annotations, ref)
                    )
                }
            )
        }
    )
}

test_that("every estimator fits a lone annotation beside a base column", {
    # One annotation and the partition it makes with its complement span the
    # same designs once the base is added; without it, LDER and HDL used to
    # model the uncovered half as carrying no heritability.
    d <- .h2RecoveryRefs(nBlocks = 10L)
    half <- as.integer(seq_len(d$M) <= d$M / 2)
    set.seed(9)
    ss <- .h2ScoreSumstats(d, d$simZ(ifelse(half == 1L, 0.5, 0.2) / d$M))
    lone <- .h2FitAll(d, ss, .h2BaselineOf(d, cbind(half = half)))
    pair <- .h2FitAll(
        d,
        ss,
        .h2BaselineOf(d, cbind(half = half, rest = 1L - half))
    )
    for (m in c("lder", "hdl")) {
        expect_equal(
            heritabilityEstimate(lone[[m]]),
            heritabilityEstimate(pair[[m]]),
            tolerance = 1e-3,
            info = m
        )
        lonePart <- stratifiedHeritabilityEnrichment(lone[[m]])
        pairPart <- stratifiedHeritabilityEnrichment(pair[[m]])
        expect_equal(
            lonePart$enrichment,
            pairPart$enrichment[pairPart$annotation == "half"],
            tolerance = 1e-3,
            info = m
        )
    }
    # Only the user's annotations are reported, never the added base.
    for (m in names(lone)) {
        expect_identical(
            stratifiedHeritabilityEnrichment(lone[[m]])$annotation,
            "half",
            info = m
        )
    }
})

test_that("an all-ones baseline fits as the unstratified model does", {
    d <- .h2RecoveryRefs(nBlocks = 10L)
    set.seed(9)
    ss <- .h2ScoreSumstats(d, d$simZ(rep(0.35 / d$M, d$M)))
    ones <- .h2FitAll(d, ss, .h2BaselineOf(d, cbind(ones = rep(1L, d$M))))
    none <- .h2FitAll(d, ss, NULL)
    for (m in names(ones)) {
        # S-LDSC's one design column is then the all-ones annotation's
        # stratified score, formed from the LD, rather than the stored base
        # score: the two agree to rounding of the LD.
        expect_equal(
            heritabilityEstimate(ones[[m]]),
            heritabilityEstimate(none[[m]]),
            tolerance = 1e-4,
            info = m
        )
        expect_equal(
            stratifiedHeritabilityEnrichment(ones[[m]])$enrichment,
            1,
            tolerance = 1e-6,
            info = m
        )
    }
})

test_that("HDL can estimate a depleted annotation", {
    # Beside the base an annotation's parameter is its excess over the
    # genome-wide rate; a lower bound of 0 used to forbid depletion.
    d <- .h2RecoveryRefs(nBlocks = 10L)
    half <- as.integer(seq_len(d$M) <= d$M / 2)
    set.seed(4)
    ss <- .h2ScoreSumstats(d, d$simZ(ifelse(half == 1L, 0.1, 0.6) / d$M))
    annot <- .h2BaselineOf(d, cbind(half = half))
    fit <- suppressWarnings(estimateStratifiedH2(
        ss,
        d$eigenRef,
        annot,
        param = .h2MethodSettings("hdl", annot, d$eigenRef)
    ))
    part <- stratifiedHeritabilityEnrichment(fit)
    expect_lt(part$tau, 0)
    expect_lt(part$enrichment, 1)
})

test_that("enrichment is reported only for columns within [0, 1]", {
    set.seed(3)
    M <- 40L
    A <- `attr<-`(
        cbind(
            base = 1,
            inside = runif(M),
            above = runif(M, 0, 2),
            below = rnorm(M)
        ),
        "userCols",
        2:4
    )
    tau <- c(1e-3, 2e-3, -1e-3, 5e-4)
    folds <- matrix(tau, 6L, 4L, byrow = TRUE) + rnorm(24L, sd = 1e-4)
    out <- pecotmr:::.h2StratifiedReport(A, tau, folds)$enrichment
    # The fit's coefficients stay meaningful outside the domain.
    expect_equal(out$tau, tau[2:4])
    expect_true(all(is.finite(out$tauSe)))
    for (col in c(
        "enrichment",
        "enrichmentSe",
        "enrichmentP",
        "propH2",
        "propSnps"
    )) {
        expect_true(is.finite(out[[col]][[1L]]), info = col)
        expect_true(all(is.na(out[[col]][2:3])), info = col)
    }
})

test_that("the enrichment p-value is ldsc's t test of inside against outside", {
    set.seed(5)
    M <- 50L
    inside <- rbinom(M, 1L, 0.3)
    A <- `attr<-`(cbind(base = 1, annot = inside), "userCols", 2L)
    tau <- c(1e-3, 2e-3)
    folds <- matrix(tau, 10L, 2L, byrow = TRUE) + rnorm(20L, sd = 2e-4)
    # regressions.py, _overlap_output: per-variant h2 inside minus outside.
    difference <- function(tt) {
        v <- as.vector(A %*% tt)
        sum(v[inside == 1L]) /
            sum(inside) -
            sum(v[inside == 0L]) / sum(1L - inside)
    }
    se <- pecotmr:::.jackknifeSeFromLoo(apply(folds, 1L, difference))
    out <- pecotmr:::.h2StratifiedReport(A, tau, folds)$enrichment
    expect_equal(
        out$enrichmentP,
        2 * pt(abs(difference(tau)) / se, 10L, lower.tail = FALSE)
    )
    # A column covering every variant has no outside to compare with.
    ones <- `attr<-`(cbind(ones = 1, annot = inside), "userCols", 1L)
    whole <- pecotmr:::.h2StratifiedReport(ones, tau, folds)$enrichment
    expect_equal(whole$enrichment, 1)
    expect_true(is.na(whole$enrichmentP))
})

test_that(".h2AnnotationBasis centers only columns whose zero is unreached", {
    M <- 20L
    bin <- rep(0:1, each = 10L)
    x <- seq(-1, 1, length.out = M)
    # A column reaching 0 keeps its heritability parametrization; the base
    # is a real variant's rate and an annotation beside it an excess.
    b <- pecotmr:::.h2AnnotationBasis(cbind(base = 1, bin = bin))
    expect_equal(b$transform, diag(2L))
    expect_equal(b$size, c(M, 10))
    expect_equal(b$share, c(1, 1))
    expect_equal(b$lower, c(0, -1))
    expect_equal(b$upper, c(1, 1))
    # Each level of a partition is its own variants' rate.
    expect_equal(
        pecotmr:::.h2AnnotationBasis(cbind(a = bin, b = 1L - bin))$lower,
        c(0, 0)
    )
    # Annotations covering every variant leave the base no variant of its own.
    expect_equal(
        pecotmr:::.h2AnnotationBasis(cbind(
            base = 1,
            a = bin,
            b = 1L - bin
        ))$lower,
        c(-1, -1, -1)
    )
    # A column above 0 throughout moves its location into the base, which
    # then reads as the total heritability.
    shifted <- pecotmr:::.h2AnnotationBasis(cbind(base = 1, x = x + 3))
    expect_equal(shifted$size, c(M, sum(abs(x))))
    expect_equal(shifted$share, c(1, 0))
    expect_equal(shifted$lower, c(0, -2))
    expect_equal(shifted$upper, c(1, 2))
    # The procedures' h2 bounds scale them: a heritability within them, an
    # excess within their span, a centered column's within twice it.
    capped <- pecotmr:::.h2AnnotationBasis(
        cbind(base = 1, x = x + 3),
        NULL,
        c(0, 0.5)
    )
    expect_equal(capped$lower, c(0, -1))
    expect_equal(capped$upper, c(0.5, 1))
    excess <- pecotmr:::.h2AnnotationBasis(
        cbind(base = 1, bin = bin),
        NULL,
        c(0, 0.5)
    )
    expect_equal(excess$lower, c(0, -0.5))
    expect_equal(excess$upper, c(0.5, 0.5))
    expect_equal(
        as.vector(shifted$transform %*% c(0.5, 0.2)),
        c(0.5 - 3 * 0.2, 0.2)
    )
    # Columns that themselves span the all-ones vector cannot all be centered.
    u <- seq(0.1, 0.9, length.out = M)
    expect_equal(
        pecotmr:::.h2AnnotationBasis(cbind(u = u, v = 1 - u))$transform,
        diag(2L)
    )
})

# A block-level annotation, which the eigen-methods can identify, with values
# that keep every per-variant rate positive, and z-scores under it.
.h2ContinuousFixture <- function() {
    d <- .h2RecoveryRefs(nBlocks = 12L)
    set.seed(21)
    raw <- rep(seq(-1.5, 1.5, length.out = 12L), each = 40L)
    x <- raw - mean(raw)
    rate <- 0.4 / d$M * (1 + 0.5 * x)
    list(d = d, x = x, ss = .h2ScoreSumstats(d, d$simZ(rate)))
}

test_that("HDL reads a continuous annotation whatever its zero", {
    f <- .h2ContinuousFixture()
    fit <- function(col) {
        annot <- .h2BaselineOf(f$d, cbind(x = col))
        suppressWarnings(estimateStratifiedH2(
            f$ss,
            f$d$eigenRef,
            annot,
            param = .h2MethodSettings("hdl", annot, f$d$eigenRef)
        ))
    }
    codings <- list(
        f$x + 3,
        f$x + 10,
        f$x,
        f$x + 0.5,
        f$x - min(f$x),
        (f$x + 3) / 100
    )
    fits <- lapply(codings, fit)
    taus <- map_dbl(fits, function(e) {
        stratifiedHeritabilityEnrichment(e)$tau
    }) *
        c(1, 1, 1, 1, 1, 1 / 100)
    # Moving the zero leaves the coefficient of x and the total alone, where
    # it used to pin the base at 0 and take tau from the shift alone.
    expect_gt(taus[[1L]], 0)
    expect_equal(taus, rep(taus[[1L]], 6L), tolerance = 0.02)
    expect_equal(
        map_dbl(fits, heritabilityEstimate),
        rep(heritabilityEstimate(fits[[1L]]), 6L),
        tolerance = 0.02
    )
    # Enrichment only for the coding inside [0, 1], and with a real p-value.
    part <- map(fits, stratifiedHeritabilityEnrichment)
    expect_true(all(is.na(map_dbl(part[1:5], "enrichment"))))
    expect_gt(part[[6L]]$enrichment, 1)
    expect_gt(part[[6L]]$enrichmentP, 0)
})

test_that("a column of mean 0 fits by every estimator, without an enrichment", {
    f <- .h2ContinuousFixture()
    fits <- .h2FitAll(f$d, f$ss, .h2BaselineOf(f$d, cbind(x = f$x)))
    for (m in names(fits)) {
        part <- stratifiedHeritabilityEnrichment(fits[[m]])
        expect_true(is.finite(heritabilityEstimate(fits[[m]])), info = m)
        expect_true(is.finite(part$tau), info = m)
        expect_true(is.na(part$enrichment), info = m)
    }
})

test_that("gLDSC's Omega floors a negative plug-in variance", {
    R <- 0.5^abs(outer(1:4, 1:4, "-"))
    negative <- pecotmr:::.h2GlsOmegaEigen(
        list(R = R),
        cbind(c(1, 1, 1, 1), c(-3, 0.5, 0, 0)),
        2
    )
    expect_true(all(is.finite(negative$values)))
    # The floor reads the variants at 0, and only those.
    floored <- pecotmr:::.h2GlsOmegaEigen(
        list(R = R),
        cbind(c(0, 1.5, 1, 1)),
        2
    )
    expect_equal(negative$values, floored$values)
})

test_that("sHDL's score test reads a candidate whatever its zero", {
    f <- .h2ContinuousFixture()
    gr <- GenomicRanges::granges(f$d$eigenRef, use.mcols = FALSE)
    candidate <- function(col) {
        cols <- cbind(base = 1, cand = col)
        AnnotationMatrix(
            cols,
            gr,
            data.frame(
                name = colnames(cols),
                tier = c("baseline", "candidate"),
                type = "continuous"
            ),
            genome = "hg19"
        )
    }
    z <- function(col) {
        annot <- candidate(col)
        fit <- suppressWarnings(estimateStratifiedH2(
            f$ss,
            f$d$eigenRef,
            annot,
            param = .h2MethodSettings("hdl", annot, f$d$eigenRef)
        ))
        scoreStats(fit)[[1L]]$z
    }
    # A mean of 0 used to scale the score by next to nothing, of either sign.
    zero <- z(f$x)
    expect_gt(zero, 0)
    expect_equal(zero, z(f$x + 3), tolerance = 0.01)
})

test_that("S-LDSC's candidate score test needs the reference's LD", {
    # The recovery reference's scores, without the LD they came from.
    d <- .h2RecoveryRefs(nBlocks = 4L)
    ref <- LdScore(
        ldBlocks = ldBlocks(d$scoreRef),
        snpInfo = data.frame(
            SNP = paste0("rs", seq_len(d$M)),
            CHR = "chr1",
            BP = as.integer(seq_len(d$M) * 100L),
            A1 = "A",
            A2 = "G"
        ),
        nSamples = 500000L,
        inSample = FALSE,
        genome = "hg19",
        ldScores = S4Vectors::mcols(d$scoreRef)$ldScores,
        ldScoreWeights = rep(1, d$M)
    )
    # Candidates only: the unstratified fit needs no LD, the score does.
    cand <- matrix(stats::rbinom(d$M, 1, 0.3), dimnames = list(NULL, "cand"))
    annot <- AnnotationMatrix(
        cand,
        GenomicRanges::GRanges(
            "chr1",
            IRanges::IRanges(as.integer(seq_len(d$M) * 100L), width = 1L)
        ),
        data.frame(name = "cand", tier = "candidate", type = "binary"),
        genome = "hg19"
    )
    expect_error(
        estimateStratifiedH2(
            .h2ScoreSumstats(d, stats::rnorm(d$M)),
            ref,
            annot,
            param = .h2MethodSettings("sldsc", annot, ref)
        ),
        "candidate-annotation score test needs each LD block's LD"
    )
})

test_that("variantGlsUnivariate recovers a known h2 (seed = 1)", {
    d <- .h2RecoveryRefs()
    uni <- rep(0.4 / d$M, d$M)
    set.seed(1)
    ests <- replicate(30, variantGlsUnivariate(d$simZ(uni), d$N, d$scoreRef)$h2)
    expect_lt(abs(mean(ests) - 0.4), 0.1) # g-LDSC has higher finite-block var
    expect_equal(mean(ests), 0.365, tolerance = 0.08)
})

test_that("eigenLikelihoodUnivariate recovers block-separable partitioned h2", {
    d <- .h2RecoveryRefs(nBlocks = 10L)
    annot <- .h2RecoveryAnnot(d$M, d$M / 2)
    perSnp <- ifelse(seq_len(d$M) <= d$M / 2, 0.6 / d$M, 0.2 / d$M)
    M_a <- c(d$M / 2, d$M / 2)
    set.seed(1)
    h2 <- rowMeans(replicate(15, {
        r <- eigenLikelihoodUnivariate(
            d$simZ(perSnp),
            d$N,
            d$eigenRef,
            annotations = annot
        )
        r$enrichment$tau * M_a
    }))
    expect_lt(abs(h2[1] - 0.3), 0.06) # annotA true h2 = 0.3
    expect_lt(abs(h2[2] - 0.1), 0.05) # annotB true h2 = 0.1
})

test_that("variantGlsUnivariate recovers enrichment's direction (A > B)", {
    d <- .h2RecoveryRefs()
    annot <- .h2RecoveryAnnot(d$M, d$M / 2)
    perSnp <- ifelse(seq_len(d$M) <= d$M / 2, 0.6 / d$M, 0.2 / d$M)
    set.seed(1)
    e <- rowMeans(replicate(20, {
        variantGlsUnivariate(
            d$simZ(perSnp),
            d$N,
            d$scoreRef,
            annotations = annot
        )$enrichment$enrichment
    }))
    expect_gt(e[1], 1.2) # annotA enriched
    expect_lt(e[2], 0.8) # annotB depleted
})

test_that("LDER uses the in-sample vs out-of-sample weight switch", {
    d <- .h2RecoveryRefs()
    set.seed(1)
    z <- d$simZ(rep(0.4 / d$M, d$M))
    refOut <- d$eigenRef
    refIn <- d$eigenRef
    refIn@inSample <- TRUE
    expect_false(isTRUE(all.equal(
        eigenRegressionUnivariate(z, d$N, refOut)$h2,
        eigenRegressionUnivariate(z, d$N, refIn)$h2
    )))
})

# The S-LDSC checks below live here rather than with the rest of the S-LDSC
# block because they need .h2RecoveryRefs(), defined just above.
# S-LDSC is a genome-scale method: real analyses regress over ~1M variants.
# Every reference in this suite is a few thousand, and at that size the
# intercept is not identifiable -- an unweighted OLS on the same data also
# returns an intercept above 3 against a true 1, and attenuates h2 to ~0.30.
# So these tests check that the estimator is in the right neighbourhood and
# that its machinery behaves, not that it hits the truth; the attenuation is
# a property of the data, not of the fit.
test_that("variantRegressionUnivariate lands near h2 on a small reference", {
    d <- .h2RecoveryRefs(N = 5000)
    uni <- rep(0.4 / d$M, d$M)
    set.seed(1)
    ests <- replicate(
        20,
        pecotmr:::variantRegressionUnivariate(d$simZ(uni), d$N, d$scoreRef)$h2
    )
    expect_gt(mean(ests), 0.15)
    expect_lt(mean(ests), 0.65)
})

test_that("S-LDSC tracks h2 upward when the true h2 rises", {
    # The check the point estimate cannot make on this scale: whatever the
    # attenuation, a bigger true h2 must give a bigger estimate.
    d <- .h2RecoveryRefs(N = 5000)
    set.seed(1)
    low <- replicate(
        10,
        pecotmr:::variantRegressionUnivariate(
            d$simZ(rep(0.1 / d$M, d$M)),
            d$N,
            d$scoreRef
        )$h2
    )
    set.seed(1)
    high <- replicate(
        10,
        pecotmr:::variantRegressionUnivariate(
            d$simZ(rep(0.6 / d$M, d$M)),
            d$N,
            d$scoreRef
        )$h2
    )
    expect_gt(mean(high), mean(low))
})

test_that("the S-LDSC intercept is estimated on the low-chi2 variants", {
    # Two-step: step 1 fits the intercept on chi2 < 30, step 2 holds it fixed
    # while fitting the slopes. The value is reported as estimated -- the
    # original ldsc does not floor it at 1, and one below 1 is informative
    # rather than an error to be clamped away.
    d <- .h2RecoveryRefs(N = 5000)
    set.seed(1)
    res <- pecotmr:::variantRegressionUnivariate(
        d$simZ(rep(0.4 / d$M, d$M)),
        d$N,
        d$scoreRef
    )
    expect_true(is.finite(res$intercept))
    expect_gt(res$intercept, 0)
})

test_that("S-LDSC recovers enrichment direction for one annotation", {
    # A single annotation, not a full partition: a set of baseline columns
    # that exactly partitions the variants is rank-deficient against the
    # model's implicit all-ones base column, and the split between them is
    # then arbitrary. Real baseline models overlap and do not partition.
    d <- .h2RecoveryRefs(N = 5000)
    M <- d$M
    inA <- as.integer(seq_len(M) <= M / 2)
    annot <- AnnotationMatrix(
        matrix(inA, ncol = 1, dimnames = list(NULL, "annotA")),
        GenomicRanges::GRanges(
            "chr1",
            IRanges::IRanges(as.integer(seq_len(M) * 100L), width = 1L)
        ),
        data.frame(name = "annotA", tier = "baseline", type = "binary"),
        genome = "hg19"
    )
    perSnp <- ifelse(seq_len(M) <= M / 2, 0.6 / M, 0.2 / M)
    set.seed(1)
    e <- mean(replicate(10, {
        pecotmr:::variantRegressionUnivariate(
            d$simZ(perSnp),
            d$N,
            d$scoreRef,
            annotations = annot
        )$enrichment$enrichment
    }))
    expect_gt(e, 1) # annotA carries 3x the per-variant heritability
})


# =============================================================================
# Numerical agreement with the original ldsc
#
# Values below were produced by running ldsc's own `ldscore.regressions.Hsq`
# (github.com/CBIIT/ldsc) on exactly the data these tests rebuild, in the
# two configurations `ldscore/sumstats.py` selects. The LD scores, base and
# stratified alike, sum ldsc's unbiased r^2 at n = 1000, as ldsc's own scores
# do, so .h2ParityRef() builds them at that n:
#
#   univariate  -> Hsq(..., twostep = 30)
#   partitioned -> chi2 < max(0.001 N, 80) filter, then Hsq(..., old_weights
#                  = True)
#
# Agreement is to about 1e-6 relative; the residual is the small ridge in
# .h2WlsCoef(), which the original has no equivalent of. Case 2's
# intercept sits below 1 and case 3's partitioned h2 is negative -- the
# original returns those too, and reproducing them is the point: the
# implementation should track the reference, pathologies included.
# =============================================================================

# The z-scores the original implementations were run on:
# z ~ N(0, N R diag(v) R + R), drawn block by block over the parity LD below
# and PINNED as data.
#
# Read from a fixture rather than re-simulated because the draw used
# `eigen()`: z = V sqrt(D) w. eigen() fixes eigenvectors only up to sign, and
# the sign LAPACK returns differs between builds, so the same seed produced
# different z -- and therefore different estimates -- on linux-64 than on
# osx-arm64 / linux-aarch64. The comparisons above are only meaningful
# against the exact vectors the original implementations were run on, so
# those are stored.
#
# Rebuild with inst/scripts/build_h2_parity_z.R, and ONLY together with
# re-deriving the expected values: the data and the expected numbers are a
# matched pair.
.h2ParityZ <- function(h2, n, seed) {
    key <- paste(h2, n, seed, sep = "_")
    zAll <- readRDS(test_path("test_data", "h2ParityZ.rds"))
    if (!(key %in% names(zAll))) {
        stop("no pinned z-scores for case '", key, "'")
    }
    zAll[[key]]
}

# The LD every pinned value in this file was derived from: the 14 chr22 LD
# blocks (1,300 variants) the bundled ldScoreExample / ldEigenExample held
# until 2026-10-20, frozen in test_data/h2ParityLd.rds as each block's strict
# upper triangle. Keeping them apart from the bundled examples lets those be
# rebuilt without moving a single pin. The S-LDSC scores are ldsc's unbiased
# r^2 at `nSamples`; eigen() fixes each eigenvector only up to sign, which
# none of the estimators depends on.
.h2ParityRef <- function(nSamples = 1000L) {
    fx <- readRDS(test_path("test_data", "h2ParityLd.rds"))
    blocks <- map(fx$blocks, function(b) {
        R <- diag(length(b$snpIdx))
        R[upper.tri(R)] <- b$upper
        R[lower.tri(R)] <- t(R)[lower.tri(R)]
        ids <- fx$snpInfo$SNP[b$snpIdx]
        list(R = `dimnames<-`(R, list(ids, ids)), snpIdx = b$snpIdx)
    })
    base <- numeric(nrow(fx$snpInfo))
    for (b in blocks) {
        base[b$snpIdx] <- rowSums(
            pecotmr:::.ldScoreUnbiasedR2(b$R^2, nSamples)
        )
    }
    eigens <- map(blocks, function(b) {
        e <- eigen(b$R, symmetric = TRUE)
        list(values = e$values, vectors = e$vectors, snpIdx = b$snpIdx)
    })
    # The frozen table holds the blocks 1-based, as they were always read.
    ldBlocks <- GenomicRanges::GRanges(
        fx$ldBlocks$chrom,
        IRanges::IRanges(fx$ldBlocks$start, fx$ldBlocks$end)
    )
    list(
        score = LdScore(
            snpInfo = fx$snpInfo,
            ldScores = matrix(base, ncol = 1, dimnames = list(NULL, "base_l2")),
            ldScoreWeights = base,
            ldBlocks = ldBlocks,
            nSamples = nSamples,
            genome = fx$genome,
            ldMatrixList = blocks
        ),
        eigen = LdEigen(
            snpInfo = fx$snpInfo,
            eigenList = eigens,
            ldBlocks = ldBlocks,
            nSamples = nSamples,
            genome = fx$genome
        )
    )
}

# The variant regression's fit through .h2Regression(), from its LD scores
# `scores` (the stored base score `baseScore` the weights divide by), the
# design `A` and the settings `param`.
.variantRegressionRowsFit <- function(chi2, scores, baseScore, n, A, param) {
    rows <- pecotmr:::.h2VariantRegressionRows(
        chi2,
        scores,
        A,
        baseScore,
        n,
        param,
        list(seq_along(chi2))
    )
    pecotmr:::.h2Regression(param, rows)
}

test_that("univariate S-LDSC reproduces the original ldsc to 1e-5", {
    ref <- .h2ParityRef()$score
    M <- length(ref)
    scores <- matrix(as.vector(mcols(ref)$ldScores[, 1]), ncol = 1)
    baseScore <- as.vector(mcols(ref)$ldScores[, 1])

    reference <- list(
        list(
            h2 = 0.4,
            n = 10000,
            seed = 42,
            h2Hat = 0.307612483,
            int = 5.724347641
        ),
        list(
            h2 = 0.15,
            n = 4000,
            seed = 7,
            h2Hat = 0.152458832,
            int = 0.826813787
        ),
        list(
            h2 = 0.6,
            n = 25000,
            seed = 99,
            h2Hat = 0.515391513,
            int = 8.460347119
        )
    )
    for (g in reference) {
        z <- .h2ParityZ(g$h2, g$n, g$seed)
        fit <- .variantRegressionRowsFit(
            z^2,
            scores,
            baseScore,
            g$n,
            matrix(1, M, 1),
            .h2GenomeWideSettings("sldsc", FALSE, FALSE)
        )
        expect_equal(fit$h2, g$h2Hat, tolerance = 1e-5)
        expect_equal(fit$intercept, g$int, tolerance = 1e-5)
    }
})

test_that("partitioned S-LDSC reproduces the original ldsc to 1e-5", {
    ref <- .h2ParityRef()$score
    M <- length(ref)
    inA <- as.integer(seq_len(M) <= M / 2)
    A <- cbind(base = 1, annotA = inA)
    scores <- pecotmr:::.h2VariantScoreMatrix(ref, A)

    reference <- list(
        list(
            h2 = 0.4,
            n = 10000,
            seed = 42,
            h2Hat = 0.366278912,
            int = 2.143763992
        ),
        list(
            h2 = 0.15,
            n = 4000,
            seed = 7,
            h2Hat = 0.146933464,
            int = 0.894710259
        ),
        list(
            h2 = 0.6,
            n = 25000,
            seed = 99,
            h2Hat = -0.042445193,
            int = 22.018356554
        )
    )
    for (g in reference) {
        z <- .h2ParityZ(g$h2, g$n, g$seed)
        fit <- .variantRegressionRowsFit(
            z^2,
            scores,
            scores[, 1],
            g$n,
            A,
            .h2GenomeWideSettings("sldsc", TRUE, FALSE)
        )
        expect_equal(fit$h2, g$h2Hat, tolerance = 1e-5)
        expect_equal(fit$intercept, g$int, tolerance = 1e-5)
    }
})

test_that("S-LDSC's two sets of settings differ, as the original ldsc's do", {
    # Unstratified, S-LDSC's settings take the two-step estimator;
    # stratified, the chi2-filtered single-pass `old_weights` route. The
    # original ldsc refuses to run two-step on a partitioned design at all,
    # so these must not converge to one shared code path.
    ref <- .h2ParityRef()$score
    M <- length(ref)
    z <- .h2ParityZ(0.4, 10000, 42)
    scores <- matrix(as.vector(mcols(ref)$ldScores[, 1]), ncol = 1)
    uni <- .variantRegressionRowsFit(
        z^2,
        scores,
        scores[, 1],
        10000,
        matrix(1, M, 1),
        .h2GenomeWideSettings("sldsc", FALSE, FALSE)
    )
    inA <- as.integer(seq_len(M) <= M / 2)
    part <- .variantRegressionRowsFit(
        z^2,
        pecotmr:::.h2VariantScoreMatrix(ref, cbind(base = 1, annotA = inA)),
        scores[, 1],
        10000,
        cbind(base = 1, annotA = inA),
        .h2GenomeWideSettings("sldsc", TRUE, FALSE)
    )
    expect_false(isTRUE(all.equal(uni$intercept, part$intercept)))
})

test_that("S-LDSC's regression settings each reach its fit", {
    dat <- simulateH2Data(n_gwas = 2000, h2_true = 0.1)
    base <- .h2GenomeWideSettings("sldsc", FALSE, FALSE)
    with <- function(...) {
        initialize(base, procedures = initialize(base@procedures, ...))
    }
    h2 <- function(param) {
        variantRegressionUnivariate(
            dat$z,
            dat$n,
            dat$ld_score_ref,
            param = param
        )$h2
    }
    published <- h2(base)
    expect_false(isTRUE(all.equal(h2(with(ldRedundancy = FALSE)), published)))
    expect_false(isTRUE(all.equal(
        h2(with(fitWeightIntercept = FALSE)),
        published
    )))
    # Updated to convergence the first weights barely matter; with no
    # updates they are the weights.
    expect_equal(
        h2(with(iterations = NULL, initialWeights = "flat")),
        h2(with(iterations = NULL)),
        tolerance = 1e-6
    )
    expect_false(isTRUE(all.equal(
        h2(with(iterations = 0, initialWeights = "flat")),
        h2(with(iterations = 0))
    )))
})

test_that("floored noise holds S-LDSC's intercept at the floor", {
    dat <- simulateH2Data(n_gwas = 2000, h2_true = 0.1)
    # Deflated statistics: the free intercept falls below 1.
    z <- dat$z * 0.8
    for (stratified in c(FALSE, TRUE)) {
        free <- .h2GenomeWideSettings("sldsc", stratified, FALSE)
        floored <- initialize(
            free,
            noise = "floored",
            procedures = initialize(free@procedures, interceptFloor = 1)
        )
        annot <- if (stratified) {
            makeCoverageAnnot(dat$n_snps, "spatial", nCand = 0L)
        }
        fit <- function(param) {
            variantRegressionUnivariate(
                z,
                dat$n,
                dat$ld_score_ref,
                annot,
                param = param
            )
        }
        expect_lt(fit(free)$intercept, 1)
        expect_equal(fit(floored)$intercept, 1)
    }
})

test_that("S-LDSC's two stages take LDER's threshold when none is set", {
    dat <- simulateH2Data()
    scores <- matrix(as.vector(mcols(dat$ld_score_ref)$ldScores[, 1]))
    param <- .h2GenomeWideSettings("sldsc", FALSE, FALSE)
    m <- dat$n_snps
    n <- 2000
    threshold <- function(chi2) {
        rows <- pecotmr:::.h2VariantRegressionRows(
            chi2,
            scores,
            matrix(1, m, 1),
            scores[, 1],
            n,
            param,
            list(seq_len(m))
        )
        pecotmr:::.h2Chi2Threshold(param, rows)
    }
    # Chi-squares the model fits exactly, so the free fit's intercept is the
    # one they are built with. No inflation: LDER's rule on the variants.
    expect_equal(
        threshold(1 + 0.05 * scores[, 1]),
        3.5 * sqrt(m / n) + 2.5
    )
    # An intercept of 2, past 1.05: the top 5 / m of the chi-squares.
    inflated <- 2 + 0.05 * scores[, 1]
    expect_equal(threshold(inflated), stats::quantile(inflated, 1 - 5 / m))
})

test_that("LDER counts a truncated LdEigen's variants, not its directions", {
    # Holding the leading directions that carry 99% of each block's LD keeps
    # about half of the parity reference's directions. The heritability is
    # per variant, so LDER's M is the variant count, as HDL's is; counting
    # the directions read it at about half. It still reads somewhat lower
    # (about 10% here), its intercept losing the near-null directions it
    # rests on.
    ld <- readRDS(test_path("test_data", "h2ParityLd.rds"))
    info <- ld$snpInfo
    gr <- GenomicRanges::GRanges(
        paste0("chr", info$CHR),
        IRanges::IRanges(info$BP, width = 1L)
    )
    names(gr) <- info$SNP
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(A1 = info$A1, A2 = info$A2)
    blocks <- map(ld$blocks, function(b) {
        R <- diag(length(b$snpIdx))
        R[upper.tri(R)] <- b$upper
        R[lower.tri(R)] <- t(R)[lower.tri(R)]
        LdMatrix(R, variants = gr[b$snpIdx], nSamples = 20000L)
    })
    z <- .h2ParityZ(0.4, 1e5, 2)
    full <- buildLdEigen(blocks)
    truncated <- buildLdEigen(blocks, eigenvalueTruncation = 0.99)
    kept <- sum(map_int(eigenList(truncated), function(e) length(e$values)))
    expect_lt(kept, 0.6 * length(truncated))
    expect_gt(
        pecotmr:::eigenRegressionUnivariate(z, 1e5, truncated)$h2 /
            pecotmr:::eigenRegressionUnivariate(z, 1e5, full)$h2,
        0.85
    )
})

test_that("LDER's threshold leaves the caller's settings unchanged", {
    # Its rough fit takes the linear taper on a copy: `slot<-` called as a
    # function would change the caller's record in place, and the LDER fit
    # after it would weigh with the wrong taper.
    dat <- simulateH2Data()
    blockStats <- .h2EigenRegressionAllBlockStats(
        dat$eigen_ref,
        dat$z,
        NULL,
        NULL
    )
    param <- .h2GenomeWideSettings("lder", FALSE, FALSE)
    rows <- pecotmr:::.h2EigenRegressionRows(
        blockStats,
        NULL,
        dat$n,
        length(dat$eigen_ref),
        param
    )
    pecotmr:::.h2Chi2Threshold(param, rows)
    expect_identical(param@procedures@eigenvalueTaper, "logistic")
    expect_equal(rows$taper, pecotmr:::.h2Taper("logistic", pmax(rows$lam, 0)))
})


# =============================================================================
# HDL: the finite-reference bound
# =============================================================================

test_that("HDL reports when its fit lands on an optimiser bound", {
    # HDL subtracts lam * h2 / nRef from the modelled variance. Once the
    # reference is small enough relative to the GWAS, that term outgrows the
    # signal, the variance is driven negative across most directions, and the
    # optimiser walks the intercept out to its upper bound -- returning a
    # heritability that reflects the bound rather than the data. It used to do
    # that silently.
    z <- .h2ParityZ(0.4, 1e5, 2)

    small <- .h2ParityRef(nSamples = 1000L)$eigen
    expect_warning(
        est <- pecotmr:::eigenLikelihoodUnivariate(z, 1e5, small)$h2,
        "sits at the top of its range"
    )

    big <- .h2ParityRef(nSamples = 20000L)$eigen
    expect_no_warning(
        ok <- pecotmr:::eigenLikelihoodUnivariate(z, 1e5, big)$h2,
        message = "sits at the top of its range"
    )
    expect_gt(ok, 0.3) # recovers the simulated 0.4
    # The bound answer is simply not the data's answer -- here it overshoots,
    # elsewhere it undershoots, which is the point of flagging it rather than
    # asserting a direction.
    expect_gt(abs(est - ok), 0.2)
})

test_that("HDL recovers a known h2 given an adequate reference", {
    ref <- .h2ParityRef(nSamples = 20000L)$eigen
    ests <- map_dbl(
        1:5,
        function(seed) {
            suppressWarnings(
                pecotmr:::eigenLikelihoodUnivariate(
                    .h2ParityZ(0.4, 1e5, seed),
                    1e5,
                    ref
                )$h2
            )
        }
    )
    expect_lt(abs(mean(ests) - 0.4), 0.1)
})


# =============================================================================
# Numerical agreement with the other three original implementations
#
# Values reproduce the reference implementations run on the data these tests
# rebuild:
#   HDL   github.com/zhenin/HDL          (HDL/R/llfun.R, optimised the same way)
#   LDER  github.com/shuangsong0110/LDER (R/get.res.R two-stage + ldscore.R)
#   gLDSC github.com/xzw20046/gldsc      (R/middle_fun.R gls.left.right +
#                                         gls.estimator)
#
# Measured agreement when this was pinned: HDL exact (0e+00), LDER machine
# precision (~1e-16), gLDSC ~3e-8 -- the last from the small ridge on the GLS
# solve, which the original has no equivalent of.
#
# Two deliberate departures, both no-ops wherever the originals are defined:
#   * LDER clamps negative eigenvalues before forming weights. The original
#     does not, and dies with "NA/NaN/Inf in 'x'" on a rank-deficient
#     reference -- verified by running it on this same data.
#   * gLDSC and S-LDSC carry a 1e-8 ridge so an exactly singular design
#     returns a value rather than NA.
# =============================================================================

test_that("HDL reproduces the original HDL's model at its maximum", {
    # The original HDL's L-BFGS-B (optim()'s defaults, one start) stops short
    # of its likelihood's maximum, which every likelihood here reaches by
    # Fisher scoring: the estimates lie within 2e-4 of the original's, at a
    # likelihood no lower than that at the original's estimate (its intercept
    # refitted). Each maximum is confirmed by Nelder-Mead (relative
    # tolerance 1e-16), at the same log-likelihood to 10 digits; the floor
    # puts kinks in the likelihood, on which a fit can stall short of it.
    # The second's lies on one, along which the log-likelihood is within its
    # rounding of the maximum over about 1e-5 of h2, so that estimate is
    # pinned to 2e-5, the others to 1e-6.
    reference <- list(
        list(
            h2 = 0.4,
            n = 1e5,
            seed = 2,
            nRef = 20000L,
            want = 0.405108188,
            within = 1e-6,
            original = 0.405111844
        ),
        list(
            h2 = 0.3,
            n = 1e4,
            seed = 5,
            nRef = 20000L,
            want = 0.286726683,
            within = 2e-5,
            original = 0.286758612
        ),
        list(
            h2 = 0.6,
            n = 5e4,
            seed = 9,
            nRef = 50000L,
            want = 0.646853841,
            within = 1e-6,
            original = 0.646809004
        )
    )
    proc <- .h2GenomeWideSettings("hdl", FALSE, FALSE)@procedures
    for (g in reference) {
        ref <- .h2ParityRef(nSamples = g$nRef)$eigen
        z <- .h2ParityZ(g$h2, g$n, g$seed)
        fit <- suppressWarnings(pecotmr:::eigenLikelihoodUnivariate(
            z,
            g$n,
            ref
        ))
        expect_equal(fit$h2, g$want, tolerance = g$within)
        expect_equal(fit$h2, g$original, tolerance = 2e-4)
        blocks <- map(
            pecotmr:::.ldEigenBlocks(ref),
            pecotmr:::.h2EigenBlockData,
            eigenRef = ref,
            z = z,
            n = g$n,
            baselineMat = NULL,
            candMat = NULL
        )
        data <- pecotmr:::.h2VarianceData(
            pecotmr:::.h2EigenDesign(blocks, length(ref), NULL, proc),
            g$n,
            g$nRef,
            proc
        )
        at <- function(theta) {
            pecotmr:::.h2Likelihood("variance", data, theta)$value
        }
        original <- pecotmr:::.h2ScoringBest(pecotmr:::.h2Scoring(
            "variance",
            data,
            list(1),
            0,
            10,
            100,
            full = c(g$original, 1),
            free = 2L
        ))
        expect_gte(at(c(fit$h2, fit$intercept)), original$value - 1e-8)
    }
})

test_that("the HDL variance floor matches the original's exp(-10)", {
    # The original llfun floors the modelled per-direction variance at
    # exp(-10). On a real reference the floor is active for a large share of
    # the eigen-directions, and the bstar^2/lamh2 term it feeds dominates the
    # likelihood -- a different floor silently moves the optimum.
    expect_equal(
        .h2GenomeWideSettings("hdl", FALSE, FALSE)@procedures@varianceFloor,
        exp(-10)
    )
})

test_that("LDER reproduces the original LDER", {
    ref <- .h2ParityRef()$eigen
    reference <- list(
        list(h2 = 0.4, n = 1e5, seed = 2, want = 0.416906429),
        list(h2 = 0.3, n = 1e4, seed = 5, want = 0.289677426),
        list(h2 = 0.6, n = 5e4, seed = 9, want = 0.626712459)
    )
    for (g in reference) {
        z <- .h2ParityZ(g$h2, g$n, g$seed)
        est <- pecotmr:::eigenRegressionUnivariate(z, g$n, ref)$h2
        expect_equal(est, g$want, tolerance = 1e-6)
    }
})

test_that("gLDSC reproduces the original gldsc", {
    ref <- .h2ParityRef()$score
    reference <- list(
        list(h2 = 0.4, n = 5e4, seed = 3, want = 0.407716183),
        list(h2 = 0.3, n = 1e4, seed = 5, want = 0.290373743)
    )
    for (g in reference) {
        z <- .h2ParityZ(g$h2, g$n, g$seed)
        est <- pecotmr:::variantGlsUnivariate(z, g$n, ref)$h2
        expect_equal(est, g$want, tolerance = 1e-6)
    }
})

test_that("stratified gLDSC reproduces the original gldsc", {
    # Pinned from gldsc's own statements (R/middle_fun.R: gls.left.right,
    # gls.estimator, oneout, with its crossprod helper as crossprod()) on a
    # binary and a continuous annotation, and its jackknife's variances.
    ref <- .h2ParityRef()$score
    M <- length(ref)
    cols <- cbind(
        halfA = as.integer(seq_len(M) %% 2L == 0L),
        probB = (seq_len(M) %% 5L) / 4
    )
    annot <- AnnotationMatrix(
        cols,
        GenomicRanges::granges(ref, use.mcols = FALSE),
        data.frame(
            name = colnames(cols),
            tier = "baseline",
            type = c("binary", "continuous")
        ),
        genome = unique(GenomeInfoDb::genome(ref))
    )
    fit <- pecotmr:::variantGlsUnivariate(
        .h2ParityZ(0.4, 5e4, 3),
        5e4,
        ref,
        annot
    )
    expect_equal(
        fit$enrichment$tau,
        c(-6.50492608301308e-06, 0.000122544026681053),
        tolerance = 1e-5
    )
    expect_equal(
        fit$enrichment$tauSe,
        c(0.000105543613297095, 5.97338520034948e-05),
        tolerance = 1e-5
    )
    expect_equal(
        fit$enrichment$enrichment,
        c(0.989676161122372, 1.09724359437187),
        tolerance = 1e-6
    )
    expect_equal(
        fit$enrichment$enrichmentSe,
        c(0.17181162331594, 0.0548583368764506),
        tolerance = 1e-6
    )
    expect_equal(fit$h2, 0.409557142849348, tolerance = 1e-6)
    expect_equal(fit$intercept, -8.37517497402338, tolerance = 1e-5)
})


# ===========================================================================
# Study resolution, genome checking and the chi2 filter fallback
# ===========================================================================

test_that("an explicit study is used without inspecting the collection", {
    expect_equal(pecotmr:::.estimateH2ResolveStudy(NULL, "myStudy"), "myStudy")
})

test_that("the genome check skips an unlabelled build and refuses a mixture", {
    gr <- GenomicRanges::GRanges(
        c("chr1", "chr2"),
        IRanges::IRanges(c(1L, 2L), width = 1L)
    )
    # No build recorded: skipped rather than guessed at or rejected.
    expect_true(pecotmr:::checkGenomeBuild(gr))
    # Two builds inside ONE object are as much a mismatch as two objects
    # disagreeing, and report through the same path.
    GenomeInfoDb::genome(gr) <- c("hg19", "hg38")
    expect_error(
        pecotmr:::checkGenomeBuild(gr),
        "Genome build mismatch"
    )
    GenomeInfoDb::genome(gr) <- "hg38"
    expect_true(pecotmr:::checkGenomeBuild(gr))
    # An unlabelled object alongside a labelled one still passes.
    bare <- GenomicRanges::GRanges("chr3", IRanges::IRanges(1L, width = 1L))
    expect_true(pecotmr:::checkGenomeBuild(gr, bare))
})

test_that("the chi2 filter is abandoned when too few variants survive it", {
    # The original ldsc drops chi2 >= 30, but with fewer than 4 survivors the
    # regression has nothing to fit, so the filter is dropped rather than
    # returning an estimate from 3 points.
    out <- .variantRegressionRowsFit(
        rep(100, 6L),
        cbind(seq_len(6L)),
        rep(1, 6L),
        1000,
        matrix(1, 6L, 1L),
        .h2GenomeWideSettings("sldsc", FALSE, FALSE)
    )
    expect_length(out$kept$y, 6L)
    expect_false(is.na(out$h2))
})


# ===========================================================================
# HDL bound diagnostics and the partitioned chi2 filter
# ===========================================================================

test_that("HDL warns only when a parameter sits at its UPPER bound", {
    # Only upper bounds are diagnostic: an intercept resting at 0 just means
    # no confounding and is routine, so a lower bound must not warn.
    f <- pecotmr:::.h2WarnIfAtBound
    expect_silent(f(c(0.1, 0.5), c(1, 2), 1000, 500))
    # At a lower bound: still silent.
    expect_silent(f(c(0, 0), c(1, 2), 1000, 500))
})

test_that("the HDL bound warning names which parameter pinned", {
    # The last parameter is the intercept; anything else is heritability.
    f <- pecotmr:::.h2WarnIfAtBound
    expect_warning(
        f(c(0.1, 2), c(1, 2), 1000, 500),
        "the fitted intercept sits at the top of its range"
    )
    expect_warning(
        f(c(1, 0.5), c(1, 2), 1000, 500),
        "the fitted heritability sits at the top of its range"
    )
})

test_that("the partitioned chi2 filter is abandoned when too few survive", {
    # Stratified S-LDSC needs more surviving variants than annotations + 2;
    # below that the filter is dropped rather than solving an underdetermined
    # system.
    out <- .variantRegressionRowsFit(
        rep(1e6, 6L),
        cbind(seq_len(6L), rev(seq_len(6L))),
        rep(1, 6L),
        1000,
        matrix(1, 6L, 2L),
        .h2GenomeWideSettings("sldsc", TRUE, FALSE)
    )
    expect_length(out$kept$y, 6L)
})

# ---------------------------------------------------------------------------
# Fallbacks and degenerate inputs: each returns a defined answer rather than
# propagating an error out of an estimator.
# ---------------------------------------------------------------------------

test_that(".rmaMetaFallbackToDL falls back to the closed-form estimator", {
    skip_if_not_installed("metafor")
    # DL never iterates, so it is a safe landing point when REML/ML fail on
    # small or near-homogeneous inputs.
    expect_warning(
        fit <- pecotmr:::.rmaMetaFallbackToDL(
            simpleError("did not converge"),
            means = c(0.1, 0.2, 0.15),
            ses = c(0.05, 0.06, 0.04),
            method = "REML"
        ),
        "method = 'REML'\\) failed .* falling back to DL"
    )
    expect_equal(fit$method, "DL")
})

test_that(".h2BaselineMat is NULL when no annotation is baseline tier", {
    n <- 5L
    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(
            start = as.integer(seq(50, by = 100, length.out = n)),
            width = 1L
        )
    )
    mat <- matrix(c(rep(1L, n), rep(0L, n)), nrow = n, ncol = 2L)
    colnames(mat) <- c("cand_A", "cand_B")
    meta <- data.frame(
        name = c("cand_A", "cand_B"),
        tier = c("candidate", "candidate"),
        type = c("binary", "binary"),
        stringsAsFactors = FALSE
    )
    am <- AnnotationMatrix(mat, gr, meta, genome = "hg19")
    # The object has annotations, but none of them are baseline.
    expect_equal(ncol(assay(baselineAnnotations(am), "annotations")), 0L)
    expect_null(pecotmr:::.h2BaselineMat(am))
    expect_null(pecotmr:::.h2BaselineMat(NULL))
})

test_that("LDER in one stage with free noise returns the free fit", {
    set.seed(4)
    m <- 60L
    lder <- .h2GenomeWideSettings("lder", FALSE, TRUE)
    param <- initialize(
        lder,
        noise = "free",
        procedures = initialize(lder@procedures, twoStage = FALSE)
    )
    rows <- pecotmr:::.h2EigenRegressionRows(
        list(list(
            x = rnorm(m),
            lam = abs(rnorm(m, 3, 0.5)),
            ldAnnot = matrix(1, nrow = m, ncol = 1L)
        )),
        NULL,
        1000,
        5000,
        param
    )
    res <- pecotmr:::.h2Regression(param, rows)
    # One pass, no chi2 outlier filtering, intercept left free.
    expect_length(res$kept$y, m)
    expect_true(is.finite(res$intercept))
})

test_that(".h2VariantBlockIndex errors when no block covers a variant", {
    n <- 6L
    snpInfo <- data.frame(
        CHR = rep("chr22", n),
        BP = as.integer(seq(1000, by = 100, length.out = n)),
        SNP = str_c("chr22:", seq(1000, by = 100, length.out = n), ":A:G"),
        A1 = rep("A", n),
        A2 = rep("G", n),
        stringsAsFactors = FALSE
    )
    mk <- function(blocks) {
        LdScore(
            snpInfo = snpInfo,
            ldScores = matrix(1, nrow = n, ncol = 1L),
            ldScoreWeights = rep(1, n),
            ldBlocks = blocks,
            nSamples = 100L,
            genome = "hg19",
            ldMatrixList = list()
        )
    }
    # With no precomputed matrices the blocks are intersected with the
    # variants; blocks on another chromosome cover nothing.
    disjoint <- suppressWarnings(
        mk(GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10000)))
    )
    expect_error(
        suppressWarnings(pecotmr:::.h2VariantBlockIndex(disjoint)),
        "no LD block covers any reference variant"
    )
    covering <- mk(GenomicRanges::GRanges("chr22", IRanges::IRanges(1, 10000)))
    idx <- pecotmr:::.h2VariantBlockIndex(covering)
    expect_equal(unname(lengths(idx)), n)
})

test_that("RmaOptions refuses the estimates and the estimator choice", {
    expect_error(RmaOptions(yi = 1), "the per-study estimates")
    expect_error(RmaOptions(sei = 1), "the per-study standard errors")
    expect_error(RmaOptions(method = "REML"), "the caller's `metaMethod`")
    # metafor::rma() takes dots, so nothing else can be rejected.
    expect_s4_class(RmaOptions(nosuch = 1), "MethodOptions")
    expect_equal(RmaOptions(test = "knha")$test, "knha")
})

test_that(".rmaMeta forwards metaArgs to metafor::rma", {
    seen <- NULL
    real <- metafor::rma
    means <- c(0.2, 0.3, 0.25, 0.4)
    ses <- c(0.1, 0.12, 0.09, 0.15)
    suppressWarnings(with_mocked_bindings(
        pecotmr:::.rmaMeta(means, ses, metaArgs = RmaOptions(test = "knha")),
        rma = function(...) {
            seen <<- list(...)
            real(...)
        },
        .package = "metafor"
    ))
    expect_equal(seen$test, "knha")
    expect_equal(seen$method, "DL")
})

test_that("the DL fallback carries metaArgs too", {
    seen <- list()
    real <- metafor::rma
    means <- c(0.2, 0.3, 0.25, 0.4)
    ses <- c(0.1, 0.12, 0.09, 0.15)
    suppressWarnings(with_mocked_bindings(
        pecotmr:::.rmaMeta(
            means,
            ses,
            method = "REML",
            metaArgs = RmaOptions(level = 90)
        ),
        rma = function(...) {
            args <- list(...)
            seen[[length(seen) + 1L]] <<- args
            # Fail the iterative estimator so the DL fallback runs.
            if (!identical(args$method, "DL")) {
                stop("forced iterative failure")
            }
            real(...)
        },
        .package = "metafor"
    ))
    expect_length(seen, 2L)
    expect_equal(seen[[2]]$method, "DL")
    expect_true(all(map_dbl(seen, "level") == 90))
})

test_that("estimateH2 reads a study split across LD blocks as one study", {
    ref <- makeTestEigenRef()
    ss <- makeTestSumstatsForRef(ref)
    gr <- sumStats(ss)
    mid <- GenomicRanges::start(gr)[[10L]]
    split <- GwasSumStats(
        studyName = "test",
        entry = list(gr),
        genome = "hg19",
        ldReference = makeTestGwasGenotypeHandle(),
        region = GenomicRanges::GRanges(
            as.character(GenomicRanges::seqnames(gr))[[1L]],
            IRanges::IRanges(c(1L, mid + 1L), c(mid, .Machine$integer.max))
        )
    )
    expect_equal(nrow(split), 2L)
    est <- estimateH2(split, ref)
    # One estimate, one row, for the study.
    expect_identical(length(est), 1L)
    expect_equal(
        heritabilityEstimate(est),
        heritabilityEstimate(estimateH2(ss, ref))
    )
    # Two studies still need one named.
    two <- GwasSumStats(
        studyName = c("a", "b"),
        entry = list(gr, gr),
        genome = "hg19",
        ldReference = makeTestGwasGenotypeHandle()
    )
    expect_error(
        estimateH2(two, ref),
        "`studyName` is required when the GwasSumStats has 2 studies"
    )
})

# ---- HDL-L's settings: one locus's heritability with its own intercept -----

# Real genotypes from eqtlRegionExample: the LD from 150 of its 415 people and
# statistics simulated on the other 265, over 300 variants.
.windowLikelihoodFixture <- function(seed, h2) {
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
    gamma <- replace(numeric(300), sample(300, 5), rnorm(5))
    g <- as.numeric(Xs %*% gamma)
    y <- g * sqrt(h2 / var(g)) + rnorm(length(g), sd = sqrt(1 - h2))
    corInput <- as.numeric(cor(Xs, y))
    list(
        ldFactor = ldFactor,
        R = crossprod(ldFactor) / (nrow(ldFactor) - 1),
        corInput = replace(corInput, !is.finite(corInput), 0),
        n = length(study)
    )
}

# The engine's fit under `settings` of the correlations `r` from a study of
# `n`, against `ld`.
.localFit <- function(settings, r, n, ld) {
    pecotmr:::.h2LocalFit(
        settings,
        pecotmr:::.corToZ(r, n),
        rep(n, length(r)),
        ld
    )
}

test_that("HDL-L's settings reproduce HDL.L()'s per-locus heritability", {
    # Pinned reference values from HDL.L()'s own statements (the HDL
    # repository's HDL/R/HDL.L.R: its log-likelihoods, LD score regression
    # start, eigen cut and optimization loops) run on these three loci.
    want <- c(0.0563749163361593, 0.1002563602000716, 0.1277057381167920)
    h2 <- c(0.05, 0.1, 0.2)
    hdl <- .h2LocalSettings("hdlL")
    for (i in 1:3) {
        fx <- .windowLikelihoodFixture(i + 1L, h2[[i]])
        expect_equal(
            .localFit(
                hdl,
                fx$corInput,
                fx$n,
                list(R = fx$R, nSamples = 150L)
            )$h2,
            want[[i]],
            tolerance = 1e-6
        )
        # The factor's eigenpairs come from its SVD, equal to rounding.
        expect_equal(
            .localFit(
                hdl,
                fx$corInput,
                fx$n,
                list(ldFactor = fx$ldFactor, nSamples = 150L)
            )$h2,
            want[[i]],
            tolerance = 1e-6
        )
    }
})

test_that("HDL-L's settings report HDL.L()'s p-value, halved, and intercept", {
    # Pinned reference values from the same HDL.L() statements on the same
    # loci: its likelihood-ratio p-value against the intercept-only null,
    # and the fitted intercept. The settings refer the ratio to the boundary
    # mixture, halving HDL.L()'s p.
    wantP <- c(5.50193943944049e-03, 1.17185324412135e-04, 7.62771871683717e-05)
    wantIntercept <- c(0.781751722851809, 1.413136925604808, 1.165862452737627)
    h2 <- c(0.05, 0.1, 0.2)
    hdl <- .h2LocalSettings("hdlL")
    full <- hdl
    full@procedures@halfBoundaryP <- FALSE
    for (i in 1:3) {
        fx <- .windowLikelihoodFixture(i + 1L, h2[[i]])
        ld <- list(R = fx$R, nSamples = 150L)
        fit <- .localFit(hdl, fx$corInput, fx$n, ld)
        expect_equal(fit$p, wantP[[i]] / 2, tolerance = 1e-6)
        expect_equal(fit$intercept, wantIntercept[[i]], tolerance = 1e-6)
        expect_equal(.localFit(full, fx$corInput, fx$n, ld)$p, wantP[[i]])
    }
})

test_that("HDL-L's settings give NA without signal, start a flat locus at 0", {
    fx <- .windowLikelihoodFixture(2L, 0.05)
    hdl <- .h2LocalSettings("hdlL")
    # No fit beats the intercept-only null when every effect is zero.
    expect_true(is.na(
        .localFit(hdl, numeric(300), fx$n, list(R = fx$R, nSamples = 150L))$h2
    ))
    # One direction: the eigen regression's start cannot tell h2 from the
    # intercept, and holds what it cannot fit at 0.
    start <- pecotmr:::.h2VarianceRegressionStart(
        list(G = cbind(1, 1), b2 = 0.01, offset = 0, floor = 1e-18),
        c(0, 0),
        c(1, 20)
    )
    expect_true(all(is.finite(start)) && any(start == 0))
    # A monomorphic variant leaves the trace short of the variant count, so
    # no eigenvalue reaches 99% of it and every one is kept.
    kept <- pecotmr:::.h2Truncate(
        pecotmr:::.h2Eigen(list(R = diag(c(1, 1, 0)))),
        hdl$truncation,
        100,
        3L
    )
    expect_length(kept$values, 3L)
})

test_that("a free intercept absorbs the reference term into the noise", {
    # The reference's -h2 lam / nRef is proportional to lam, as the noise
    # term is, so with the intercept free only the intercept moves: by h2 n /
    # nRef, the variance the reference term takes away.
    fx <- .windowLikelihoodFixture(3L, 0.1)
    hdl <- .h2LocalSettings("hdlL")
    noRef <- hdl
    noRef@referenceNoise <- "none"
    ld <- list(R = fx$R, nSamples = 150L)
    a <- .localFit(hdl, fx$corInput, fx$n, ld)
    b <- .localFit(noRef, fx$corInput, fx$n, ld)
    expect_equal(a$h2, b$h2, tolerance = 1e-6)
    expect_equal(
        a$intercept - b$intercept,
        a$h2 * fx$n / 150,
        tolerance = 1e-4
    )
})

# ---- the moments (LAVA's settings) and the null space (FastHer's) ----

test_that("the moments follow the reference's distance from the study's LD", {
    fx <- .windowLikelihoodFixture(3L, 0.1)
    z <- fx$corInput * sqrt(fx$n)
    e <- eigen(fx$R, symmetric = TRUE)
    lam <- e$values[e$values > 0]
    K <- which(cumsum(lam) / sum(lam) >= 0.99)[[1L]]
    r <- z / sqrt(z^2 + fx$n - 2)
    delta <- crossprod(e$vectors[, seq_len(K)], r) / sqrt(lam[seq_len(K)])
    dtd <- sum(delta^2)
    m <- K / (fx$n - 1)
    lava <- .h2LocalSettings("lava")
    # A reference of N holding a share rho of the study's samples, phi being
    # 1 + (1 - 2 rho) n / N. Half the study's samples give 1 at any N.
    cases <- list(
        list(N = fx$n, rho = 1, phi = 0),
        list(N = 2 * fx$n, rho = 1, phi = 0.5),
        list(N = 400L, rho = 0.5, phi = 1),
        list(N = fx$n, rho = 0, phi = 2),
        list(N = fx$n, rho = 0.25, phi = 1.5)
    )
    for (case in cases) {
        fit <- pecotmr:::.h2LocalFit(
            lava,
            z,
            rep(fx$n, length(z)),
            list(
                ldFactor = fx$ldFactor,
                nSamples = case$N,
                referenceOverlap = case$rho
            )
        )
        v <- 1 + m * (case$phi - 1)
        s <- (1 - m) * (fx$n + 1) / (2 * v) - 1
        expect_equal(fit$nComponents, K)
        expect_equal(fit$h2, (dtd - m) / v, tolerance = 1e-8)
        expect_equal(
            fit$p,
            pbeta(dtd, m * s, (1 - m) * s, lower.tail = FALSE),
            tolerance = 1e-8
        )
    }
    # No overlap is the default, and what it assumes is flagged only while the
    # reference is small enough for that to matter.
    bare <- function(N) {
        pecotmr:::.h2LocalFit(
            lava,
            z,
            rep(fx$n, length(z)),
            list(ldFactor = fx$ldFactor, nSamples = N)
        )
    }
    expect_equal(bare(fx$n)$h2, (dtd - m) / (1 + m), tolerance = 1e-8)
    expect_true(bare(2L * fx$n)$smallReference)
    expect_false(bare(2L * fx$n + 1L)$smallReference)
    # Without the overlap factor the reference has no noise of its own: phi
    # is 1, and the estimate delta'delta - m.
    noFactor <- lava
    noFactor@referenceNoise <- "none"
    expect_equal(
        pecotmr:::.h2LocalFit(
            noFactor,
            z,
            rep(fx$n, length(z)),
            list(ldFactor = fx$ldFactor, nSamples = fx$n)
        )$h2,
        dtd - m,
        tolerance = 1e-8
    )
})

test_that("LAVA's test has no Beta to fit at a very small n", {
    expect_true(is.na(pecotmr:::.h2MomentsP(0.5, 0.8, 1, 6)))
    expect_true(is.finite(pecotmr:::.h2MomentsP(0.5, 0.8, 1, 100)))
})

test_that("the moments hold on real genotypes where delta'delta passes 1", {
    # A panel of 3000 from one population, on real genotypes, with K / n near
    # 0.6: the noise then fills most of delta'delta and the signal the rest.
    data(eqtlRegionExample, envir = environment())
    af <- colMeans(eqtlRegionExample$X) / 2
    common <- eqtlRegionExample$X[, pmin(af, 1 - af) > 0.05]
    set.seed(6L)
    win <- common[, 301:1500]
    n <- 150L
    Xpop <- scale(win)
    Rpop <- cor(win)
    b <- replace(numeric(ncol(win)), sample.int(ncol(win), 3L), rnorm(3L))
    b <- b * sqrt(0.5 / sum(b * (Rpop %*% b)))
    qtl <- sample.int(nrow(win), n, replace = TRUE)
    y <- as.vector(Xpop[qtl, ] %*% b) + rnorm(n, sd = sqrt(0.5))
    r <- as.vector(cor(win[qtl, ], y))
    ld <- list(
        R = cor(win[sample.int(nrow(win), 3000L, replace = TRUE), ]),
        nSamples = 3000L
    )
    fit <- .localFit(.h2LocalSettings("lava"), r, n, ld)
    e <- eigen(ld$R, symmetric = TRUE)
    top <- seq_len(fit$nComponents)
    dtd <- sum((crossprod(e$vectors[, top], r) / sqrt(e$values[top]))^2)
    expect_gt(dtd, 1)
    expect_lt(abs(fit$h2 - 0.5), 0.2)
    expect_lt(fit$p, 0.05)
})

test_that("the moments are NA with too few LD components or samples", {
    lava <- .h2LocalSettings("lava")
    # Two copies of one variant hold all the variance in one component.
    x <- c(1, -1, 2, 0, -2, 1, -1, 0)
    fit <- pecotmr:::.h2LocalFit(
        lava,
        c(2, 2),
        c(100, 100),
        list(ldFactor = cbind(x, x), nSamples = 8L)
    )
    expect_true(is.na(fit$h2) && is.na(fit$p))
    expect_equal(fit$nComponents, 1L)
    # Two components leave no degrees of freedom in a study of 3.
    ldFactor <- cbind(c(1, 1, -1, -1), c(1, -1, 1, -1))
    fit <- pecotmr:::.h2LocalFit(
        lava,
        c(1, 1),
        c(3, 3),
        list(ldFactor = ldFactor, nSamples = 4L)
    )
    expect_true(is.na(fit$h2) && is.na(fit$p))
    expect_equal(fit$nComponents, 2L)
})

# Standardized so that X'X = n R, X'y = n r and y'y = n, the
# summary-statistics likelihood is that of y ~ N(0, d (h2 X X' / m + (1 - h2)
# I)) itself, r being the correlation behind each t statistic.
.nullSpaceFixture <- function() {
    set.seed(9)
    n <- 300
    m <- 8
    X <- matrix(rnorm(n * m), n, m) %*% chol(0.5^abs(outer(1:m, 1:m, "-")))
    X <- scale(X) * sqrt(n / (n - 1))
    y <- as.vector(X %*% rnorm(m, sd = 0.15)) + rnorm(n)
    y <- (y - mean(y)) / sqrt(mean((y - mean(y))^2))
    list(X = X, y = y, n = n, m = m, r = as.vector(crossprod(X, y)) / n)
}

test_that("FastHer's settings are the individual-level maximum likelihood", {
    fx <- .nullSpaceFixture()
    fit <- .localFit(
        .h2LocalSettings("fastHer"),
        fx$r,
        fx$n,
        list(R = crossprod(fx$X) / fx$n)
    )
    e <- eigen(tcrossprod(fx$X) / fx$m, symmetric = TRUE)
    yt <- as.vector(crossprod(e$vectors, fx$y))
    nll <- function(par) {
        v <- par[[2L]] * (par[[1L]] * pmax(e$values, 0) + 1 - par[[1L]])
        sum(log(v) + yt^2 / v) / 2
    }
    opt <- optim(
        c(0.3, 1),
        nll,
        method = "L-BFGS-B",
        lower = c(1e-9, 1e-6),
        upper = c(1 - 1e-9, 10),
        control = list(factr = 1e2, pgtol = 0)
    )
    expect_equal(fit$h2, opt$par[[1L]], tolerance = 1e-5)
    expect_equal(fit$scale, opt$par[[2L]], tolerance = 1e-5)
    expect_equal(fit$nComponents, fx$m)
    expect_gt(fit$h2Se, 0)
    # Held at d = 1 the scale is not fitted, and the test is the ratio.
    held <- .h2LocalSettings("fastHer")
    held@procedures@profileScale <- FALSE
    held@inference <- "lrt"
    fixed <- .localFit(held, fx$r, fx$n, list(R = crossprod(fx$X) / fx$n))
    nllFixed <- function(h2) nll(c(h2, 1))
    want <- optimize(nllFixed, c(0, 1 - 1e-9), tol = 1e-10)$minimum
    expect_equal(fixed$h2, want, tolerance = 1e-6)
    expect_true(is.na(fixed$scale) && is.finite(fixed$p))
})

# A strong QTL on real genotypes, with its LD estimated in the same samples:
# the correlations, their t statistics and the LD.
.cisH2StrongQtl <- function(seed, h2 = 0.6) {
    data(eqtlRegionExample, envir = environment())
    G <- eqtlRegionExample$X
    af <- colMeans(G) / 2
    set.seed(seed)
    win <- G[, pmin(af, 1 - af) > 0.05][, 301:600]
    Xs <- scale(win) * sqrt(nrow(win) / (nrow(win) - 1))
    g <- Xs[, sample(ncol(Xs), 1L)]
    y <- sqrt(h2) * g + rnorm(length(g), sd = sqrt(1 - h2))
    r <- as.vector(cor(Xs, y))
    n <- nrow(Xs)
    list(
        r = r,
        z = r * sqrt((n - 2) / (1 - r^2)),
        n = n,
        Xs = Xs,
        y = y,
        R = cor(win)
    )
}

test_that("every setting reads a t statistic as its exact correlation", {
    fx <- .cisH2StrongQtl(11L)
    ld <- list(R = fx$R, nSamples = fx$n)
    hdl <- .h2LocalSettings("hdlL")
    fit <- pecotmr:::.h2LocalFit(hdl, fx$z, rep(fx$n, length(fx$z)), ld)
    expect_equal(fit$h2, .localFit(hdl, fx$r, fx$n, ld)$h2)
    # z / sqrt(n), as HDL.L() and FastHer read it, overstates r^2 by 1 / (1 -
    # r^2): past 1 at this QTL's strongest variant.
    expect_equal(pecotmr:::.zToCor(fx$z, fx$n), fx$r)
    expect_gt(max(abs(fx$z)) / sqrt(fx$n), 1)
})

test_that("the null space's likelihood is defined for a strong QTL in-sample", {
    # Reading t / sqrt(n) as the correlation overstates r^2 by 1 / (1 - r^2),
    # which took the statistics' norm past the sample size.
    fx <- .cisH2StrongQtl(12L)
    nVariant <- rep(fx$n, length(fx$z))
    fastHer <- .h2LocalSettings("fastHer")
    expect_gt(fx$n * max(fx$r^2) / (1 - max(fx$r^2)), 0.5 * fx$n)
    fit <- pecotmr:::.h2LocalFit(fastHer, fx$z, nVariant, list(R = fx$R))
    expect_true(is.finite(fit$h2))
    undefined <- pecotmr:::.h2LocalFit(
        fastHer,
        fx$r * fx$n,
        nVariant,
        list(R = diag(length(fx$z)))
    )
    expect_true(is.na(undefined$h2) && undefined$undefined)
})

test_that("the in-sample moments are the adjusted R^2 of OLS on components", {
    fx <- .cisH2StrongQtl(13L, h2 = 0.2)
    nVariant <- rep(fx$n, length(fx$z))
    fit <- pecotmr:::.h2LocalFit(
        .h2LocalSettings("lava"),
        fx$z,
        nVariant,
        list(R = fx$R, nSamples = fx$n, referenceOverlap = 1)
    )
    e <- eigen(fx$R, symmetric = TRUE)
    K <- fit$nComponents
    pcs <- fx$Xs %*% e$vectors[, seq_len(K)]
    ols <- summary(lm(fx$y ~ pcs))
    expect_equal(fit$h2, ols$adj.r.squared, tolerance = 1e-6)
    expect_equal(
        fit$p,
        pf(
            ols$fstatistic[[1L]],
            ols$fstatistic[[2L]],
            ols$fstatistic[[3L]],
            lower.tail = FALSE
        ),
        tolerance = 1e-6
    )
})

test_that("the null space's likelihood is NA where it is undefined", {
    set.seed(10)
    ldFactor <- matrix(rnorm(60 * 30), 60, 30)
    fastHer <- .h2LocalSettings("fastHer")
    # A reference of rank 30 for a study of 20.
    expect_true(is.na(
        pecotmr:::.h2LocalFit(
            fastHer,
            rnorm(30),
            rep(20, 30),
            list(ldFactor = ldFactor)
        )$h2
    ))
    # Statistics larger than the sample size allows.
    expect_true(is.na(
        pecotmr:::.h2LocalFit(
            fastHer,
            rep(20, 30),
            rep(100, 30),
            list(ldFactor = ldFactor)
        )$h2
    ))
})

# ---- the local settings that follow the genome-wide estimators ----

# The local eigen regression, the variant GLS and the variant regression, as
# the documentation describes them.
.localSettings <- function() {
    list(
        eigenRegression = H2EstimateParam(
            criterion = "regression",
            noise = "free",
            crossNoise = "window",
            referenceNoise = "none",
            truncation = H2TruncationParam(
                of = "variants",
                maxPerSample = NULL,
                minimum = 1
            ),
            sampleSize = "median",
            inference = "information"
        ),
        variantGls = H2EstimateParam(
            statistic = "variant",
            criterion = "gls",
            noise = "free",
            crossNoise = "window",
            referenceNoise = "none",
            sampleSize = "median",
            inference = "information"
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
                iterations = 2,
                jackknifeBlocks = 10
            )
        )
    )
}

test_that("the eigen regression to convergence is the free-noise likelihood", {
    # At convergence its weights are the inverse variances, so its normal
    # equations are the likelihood's score: with the same directions and
    # sample size, HDL-L's estimate, the reference term going to the noise.
    for (seed in 2:4) {
        fx <- .windowLikelihoodFixture(seed, 0.1)
        ld <- list(R = fx$R, nSamples = 150L)
        reg <- .localFit(
            .localSettings()$eigenRegression,
            fx$corInput,
            fx$n,
            ld
        )
        hdl <- .localFit(.h2LocalSettings("hdlL"), fx$corInput, fx$n, ld)
        expect_equal(reg$h2, hdl$h2, tolerance = 1e-5, info = seed)
        expect_true(is.finite(reg$h2Se) && reg$h2Se > 0)
        expect_true(is.na(reg$p))
    }
})

test_that("an eigen regression to convergence halves steps that cycle", {
    # One causal variant on real LD (p >> n) whose signal sits on the top
    # eigenvalues: plain weight updates cycle between 0.336, from flat weights,
    # and -0.043, whose weights (h2 held at 0) are flat again. A fit to
    # convergence halves a step that raises the deviance, and reaches the
    # free-noise likelihood's estimate; a fixed count of updates, as the
    # published procedures run, still takes each update whole.
    data(eqtlRegionExample)
    G <- eqtlRegionExample$X
    set.seed(300001L)
    idx <- sample.int(nrow(G), 300L, replace = TRUE)
    af <- colMeans(G[idx, ]) / 2
    Xs <- scale(G[idx, pmin(af, 1 - af) >= 0.05])
    causal <- sample.int(ncol(Xs), sample(1:2, 1L))
    g <- as.vector(Xs[, causal, drop = FALSE] %*% rnorm(length(causal)))
    y <- scale(g * sqrt(0.1 / var(g)) + rnorm(300L) * sqrt(0.9))
    r <- as.vector(cor(Xs, y))
    ld <- list(ldFactor = Xs, nSamples = 300L, referenceOverlap = 1)
    reg <- H2EstimateParam(
        criterion = "regression",
        noise = "free",
        referenceNoise = "eigenTerm",
        truncation = H2TruncationParam(
            of = "eigenvalues",
            maxPerSample = 0.75,
            minimum = 2
        ),
        inference = "information"
    )
    capped <- function(k) {
        updates <- H2RegressionParam(iterations = k)
        .localFit(
            methods::initialize(reg, procedures = updates),
            r,
            300L,
            ld
        )$h2
    }
    expect_gt(abs(capped(100) - capped(101)), 0.3)
    expect_equal(
        .localFit(reg, r, 300L, ld)$h2,
        .localFit(.h2LocalSettings("hdlL"), r, 300L, ld)$h2,
        tolerance = 1e-4
    )
})

test_that("the local GLS and variant regression estimate a window", {
    fx <- .windowLikelihoodFixture(4L, 0.2)
    ld <- list(R = fx$R, nSamples = 150L)
    settings <- .localSettings()
    gls <- .localFit(settings$variantGls, fx$corInput, fx$n, ld)
    expect_true(is.finite(gls$h2) && is.finite(gls$intercept))
    expect_gt(gls$h2Se, 0)
    reg <- .localFit(settings$variantRegression, fx$corInput, fx$n, ld)
    expect_true(is.finite(reg$h2) && is.finite(reg$intercept))
    expect_gt(reg$h2Se, 0)
    # A jackknife needs blocks to leave out.
    noBlocks <- settings$variantRegression
    noBlocks@procedures@jackknifeBlocks <- NULL
    expect_error(
        .localFit(noBlocks, fx$corInput, fx$n, ld),
        "needs `jackknifeBlocks`"
    )
})

test_that("the regression's noise can be tied, floored or fitted in stages", {
    fx <- .windowLikelihoodFixture(3L, 0.1)
    ld <- list(R = fx$R, nSamples = 150L)
    base <- .localSettings()$eigenRegression
    free <- .localFit(base, fx$corInput, fx$n, ld)
    tied <- base
    tied@noise <- "tied"
    tiedFit <- .localFit(tied, fx$corInput, fx$n, ld)
    expect_true(is.finite(tiedFit$h2) && is.na(tiedFit$intercept))
    # A floor above the free intercept holds it there.
    floored <- base
    floored@noise <- "floored"
    floored@procedures@interceptFloor <- free$intercept + 0.5
    flooredFit <- .localFit(floored, fx$corInput, fx$n, ld)
    expect_equal(flooredFit$intercept, free$intercept + 0.5)
    # Below it, the free fit stands.
    floored@procedures@interceptFloor <- free$intercept - 0.5
    expect_equal(.localFit(floored, fx$corInput, fx$n, ld)$h2, free$h2)
    # Two stages hold the noise of a first fit below the threshold.
    staged <- base
    staged@procedures@twoStage <- TRUE
    staged@procedures@chi2Threshold <- 1e6
    expect_equal(
        .localFit(staged, fx$corInput, fx$n, ld)$h2,
        free$h2,
        tolerance = 1e-6
    )
    # LDER's own threshold, and fixed weights, still give an estimate.
    staged@procedures@chi2Threshold <- NULL
    expect_true(is.finite(.localFit(staged, fx$corInput, fx$n, ld)$h2))
    fixed <- base
    fixed@procedures@weights <- "fixed"
    expect_true(is.finite(.localFit(fixed, fx$corInput, fx$n, ld)$h2))
    # A filter dropping every chi-square above 1 changes the fit.
    filtered <- base
    filtered@procedures@chi2Filter <- c(0, 1)
    expect_false(isTRUE(all.equal(
        .localFit(filtered, fx$corInput, fx$n, ld)$h2,
        free$h2
    )))
})

test_that("the weights read the intercept and start as the settings say", {
    fx <- .windowLikelihoodFixture(3L, 0.1)
    ld <- list(R = fx$R, nSamples = 150L)
    eigen <- .localSettings()$eigenRegression
    free <- .localFit(eigen, fx$corInput, fx$n, ld)
    # Weights holding the intercept at 1, as LDER's do, move the estimate.
    held <- eigen
    held@procedures@fitWeightIntercept <- FALSE
    expect_gt(abs(.localFit(held, fx$corInput, fx$n, ld)$h2 - free$h2), 1e-3)
    # To convergence the first weights do not matter.
    moments <- eigen
    moments@procedures@initialWeights <- "moments"
    expect_equal(
        .localFit(moments, fx$corInput, fx$n, ld)$h2,
        free$h2,
        tolerance = 1e-6
    )
    # Without an update they are the fit's own.
    variant <- .localSettings()$variantRegression
    variant@inference <- "information"
    variant@procedures@iterations <- 0
    flat <- variant
    flat@procedures@initialWeights <- "flat"
    expect_gt(
        abs(
            .localFit(variant, fx$corInput, fx$n, ld)$h2 -
                .localFit(flat, fx$corInput, fx$n, ld)$h2
        ),
        1e-3
    )
})

test_that("the eigen statistic keeps the directions its truncation says", {
    dec <- list(values = c(4, 3, 2, 1, 0), vectors = diag(5))
    trunc <- function(..., M = 5L) {
        pecotmr:::.h2Truncate(dec, H2TruncationParam(...), 100, M)$values
    }
    # A share of the positive eigenvalues' sum, of the variant count, every
    # direction.
    expect_identical(trunc(fraction = 0.6, of = "eigenvalues"), c(4, 3))
    expect_identical(trunc(fraction = 0.6, of = "variants"), 4)
    expect_identical(trunc(fraction = 1), c(4, 3, 2, 1, 0))
    # A window of more variants than the directions reach keeps them all.
    expect_identical(
        trunc(fraction = 0.6, of = "variants", M = 20L),
        c(4, 3, 2, 1, 0)
    )
    # A cap at a share of the sample size, never under the minimum.
    expect_identical(
        pecotmr:::.h2Truncate(
            dec,
            H2TruncationParam(
                of = "eigenvalues",
                maxPerSample = 0.01,
                minimum = 2
            ),
            100,
            5L
        )$values,
        c(4, 3)
    )
    # Ties with the last direction kept.
    tied <- list(values = c(4, 2, 2, 1), vectors = diag(4))
    expect_identical(
        pecotmr:::.h2Truncate(
            tied,
            H2TruncationParam(
                fraction = 0.6,
                of = "eigenvalues",
                keepTies = TRUE
            ),
            100,
            4L
        )$values,
        c(4, 2, 2)
    )
})

test_that("a local estimate refuses the settings it cannot be made under", {
    refuse <- pecotmr:::.h2RefuseLocalSetting
    expect_error(
        refuse(
            H2EstimateParam(
                criterion = "regression",
                noise = "supplied",
                referenceNoise = "none",
                inference = "information"
            ),
            FALSE,
            "test"
        ),
        "supplied noise needs a value"
    )
    expect_error(
        refuse(
            H2EstimateParam(
                criterion = "likelihood",
                noise = "free",
                referenceNoise = "eigenTerm",
                inference = "jackknife"
            ),
            FALSE,
            "test"
        ),
        "only the variant regression has blocks"
    )
    reml <- .h2LocalSettings("greml")
    expect_error(refuse(reml, FALSE, "test"), "needs individual-level data")
    expect_null(refuse(reml, TRUE, "test"))
    expect_error(
        pecotmr:::.h2LocalFit(reml, 1, 10, list()),
        "needs individual-level data"
    )
})

# ---- estimateLocalH2 on QTL summary statistics: one cis estimate per gene ----

# qtlSumStatsExample rebuilt from its parts, with the entry's sample sizes,
# the LD reference, the QC record or the trait position replaced.
.cisH2Example <- function(
    N = NULL,
    ref = NULL,
    qcOptions = list(qc = "ok"),
    traitPos = NULL
) {
    data(qtlSumStatsExample, envir = environment())
    gr <- qtlSumStatsExample[[1]]
    if (!is.null(N)) {
        S4Vectors::mcols(gr)$N <- N
    }
    QtlSumStats(
        studyName = "study1",
        context = "brain",
        trait = "ENSG_example",
        entry = list(gr),
        genome = "hg19",
        ldReference = ref %||% ldReference(qtlSumStatsExample)[[1]],
        qcOptions = qcOptions,
        traitPos = traitPos
    )
}

test_that("estimateLocalH2 gives each QtlSumStats row its cis estimate", {
    skip_if_not_installed("PlinkMatrix")
    x <- .cisH2Example()
    methods <- c("lava", "hdlL", "fastHer")
    fits <- map(
        methods,
        function(m) {
            estimateLocalH2(x, .h2LocalSettings(m), referenceOverlap = 1)
        }
    )
    for (i in seq_along(fits)) {
        fit <- fits[[i]]
        expect_s4_class(fit, "H2Estimate")
        expect_identical(length(fit), 1L)
        expect_identical(as.character(fit$trait), "ENSG_example")
        expect_identical(nSnps(fit), 200L)
        expect_true(is.finite(heritabilityEstimate(fit)))
        expect_identical(h2EstimateParam(fit), .h2LocalSettings(methods[[i]]))
        # The element is the window the estimate covers: the row's region.
        expect_identical(
            IRanges::ranges(fit[[1]]),
            IRanges::ranges(x$region)
        )
        expect_identical(unname(GenomeInfoDb::genome(fit)), "hg19")
    }
    lava <- fits[[1L]]
    # The moments report their estimate past K = 0.1 n (68 of 165 here) with
    # their test, and no standard error.
    expect_gt(lava$n_components / 165, 0.1)
    expect_true(is.finite(lava$p_value) && is.na(lava$h2_se))
    hdl <- fits[[2L]]
    expect_true(is.finite(hdl$p_value) && is.finite(hdl$intercept))
    expect_true(is.na(hdl$h2_se) && is.na(hdl$scale))
    fastHer <- fits[[3L]]
    expect_true(is.finite(fastHer$h2_se) && is.finite(fastHer$scale))
    expect_true(is.na(fastHer$p_value))
    # The default is HDL-L's settings.
    expect_identical(estimateLocalH2(x, referenceOverlap = 1), hdl)
    # One estimate per selected row.
    two <- c(x, `$<-`(x, "trait", "ENSG_other"))
    both <- estimateLocalH2(two, referenceOverlap = 1)
    expect_identical(as.character(both$trait), c("ENSG_example", "ENSG_other"))
    expect_equal(heritabilityEstimate(both), rep(heritabilityEstimate(hdl), 2))
    expect_identical(
        as.character(
            estimateLocalH2(
                two,
                traitId = "ENSG_other",
                referenceOverlap = 1
            )$trait
        ),
        "ENSG_other"
    )
    # Any setting that fits a window: the local ones too.
    for (settings in .localSettings()) {
        est <- estimateLocalH2(x, settings, referenceOverlap = 1)
        expect_true(is.finite(heritabilityEstimate(est)))
        expect_true(is.finite(heritabilityEstimateSe(est)))
        expect_true(is.na(est$reference_overlap))
    }
})

test_that("estimateLocalH2 reads the overlap, per row, for the factor", {
    skip_if_not_installed("PlinkMatrix")
    x <- .cisH2Example()
    two <- c(x, `$<-`(x, "trait", "ENSG_other"))
    # LAVA's settings read it; the default's do not.
    factor <- .h2LocalSettings("lava")
    lava <- estimateLocalH2(x, factor, referenceOverlap = 1)
    own <- heritabilityEstimate(lava)
    half <- heritabilityEstimate(
        estimateLocalH2(x, factor, referenceOverlap = 0.5)
    )
    # This reference is as large as the study: it is the study's own samples at
    # 1 (phi 0), and half of them, with others, at 0.5 (phi 1).
    expect_equal(own, half / (1 - lava$n_components / (165 - 1)))
    expect_identical(lava$reference_overlap, 1)
    # One per selected row.
    expect_equal(
        heritabilityEstimate(estimateLocalH2(
            two,
            factor,
            referenceOverlap = c(1, 0.5)
        )),
        c(own, half)
    )
    expect_equal(
        heritabilityEstimate(
            estimateLocalH2(
                two,
                factor,
                traitId = "ENSG_other",
                referenceOverlap = 0.5
            )
        ),
        half
    )
    expect_error(
        estimateLocalH2(two, referenceOverlap = c(1, 0.5, 0)),
        "has 3 values for 2 selected"
    )
    expect_error(estimateLocalH2(x, referenceOverlap = 2), "referenceOverlap")
    expect_error(
        estimateLocalH2(x, referenceOverlap = NA_real_),
        "referenceOverlap"
    )
    # A reference cannot hold more of the study than it has samples.
    big <- .cisH2Example(N = rep(1000L, 200L))
    expect_error(
        estimateLocalH2(big, factor, referenceOverlap = 1),
        "puts 1000 of the 1000 study samples in an LD reference of 165"
    )
    # HDL-L's and FastHer's settings do not read it.
    for (m in c("hdlL", "fastHer")) {
        settings <- .h2LocalSettings(m)
        expect_equal(
            heritabilityEstimate(estimateLocalH2(x, settings)),
            heritabilityEstimate(
                estimateLocalH2(x, settings, referenceOverlap = 1)
            )
        )
    }
})

test_that("estimateLocalH2 warns when the factor assumes no overlap", {
    skip_if_not_installed("PlinkMatrix")
    x <- .cisH2Example()
    two <- c(x, `$<-`(x, "trait", "ENSG_other"))
    factor <- .h2LocalSettings("lava")
    expect_warning(
        estimateLocalH2(two, factor),
        "2 of 2 gene\\(s\\) has at most twice"
    )
    expect_warning(estimateLocalH2(x, factor), "ENSG_example")
    # Saying so, or a reference big enough, ends it; the other settings, the
    # default's included, never had it.
    expect_silent(estimateLocalH2(x, factor, referenceOverlap = 1))
    expect_silent(estimateLocalH2(x, factor, referenceOverlap = 0.5))
    expect_silent(estimateLocalH2(x))
    smallStudy <- .cisH2Example(N = rep(80L, 200L))
    expect_silent(estimateLocalH2(smallStudy, factor))
})

test_that("estimateLocalH2 analyses a row over its cis window", {
    skip_if_not_installed("PlinkMatrix")
    tp <- GenomicRanges::GRanges(
        "chr22",
        IRanges::IRanges(15200000L, width = 1L)
    )
    x <- .cisH2Example(traitPos = tp)
    fit <- estimateLocalH2(x, cisWindow = 1e5, referenceOverlap = 1)
    expect_identical(
        IRanges::ranges(fit[[1]]),
        IRanges::IRanges(15100000L, 15300000L)
    )
    expect_lt(nSnps(fit), 200L)
    expect_error(
        estimateLocalH2(.cisH2Example(), cisWindow = 1e5),
        "trait position"
    )
    # A window between two variants holds none: NA, with a warning.
    pos <- GenomicRanges::start(x[[1]])
    gap <- which.max(diff(pos))
    between <- GenomicRanges::GRanges(
        "chr22",
        IRanges::IRanges(pos[[gap]] + diff(pos)[[gap]] %/% 2L, width = 1L)
    )
    empty <- expect_warning(
        estimateLocalH2(.cisH2Example(traitPos = between), cisWindow = 10),
        "1 of 1 gene\\(s\\) hold no variants in their window"
    )
    expect_identical(nSnps(empty), 0L)
    expect_true(is.na(heritabilityEstimate(empty)))
    expect_identical(IRanges::width(empty[[1]]), 21L)
})

test_that("estimateLocalH2 marks the genes the null space cannot fit", {
    # Twenty study samples against a 165-sample reference.
    skip_if_not_installed("PlinkMatrix")
    x <- .cisH2Example(N = 20L)
    out <- expect_warning(
        estimateLocalH2(x, .h2LocalSettings("fastHer")),
        "the null space is undefined for 1 of 1 gene"
    )
    expect_true(is.na(heritabilityEstimate(out)))
})

test_that("estimateLocalH2 refuses what a cis estimate cannot be made from", {
    skip_if_not_installed("PlinkMatrix")
    x <- .cisH2Example()
    expect_error(
        estimateLocalH2(.cisH2Example(qcOptions = list())),
        "has no QC record"
    )
    expect_error(
        estimateLocalH2(x, traitId = "ENSG_other"),
        "no rows matched"
    )
    expect_error(estimateLocalH2(x, CredibleSetParam()), "H2EstimateParam")
    # An empty list carries no settings: the defaults.
    expect_identical(
        h2EstimateParam(estimateLocalH2(x, list(), referenceOverlap = 1)),
        H2EstimateParam()
    )
    expect_error(
        estimateLocalH2(x, .h2LocalSettings("greml")),
        "needs individual-level data"
    )
    data(qtlSumStatsExample, envir = environment())
    gr <- qtlSumStatsExample[[1]]
    noRef <- QtlSumStats(
        studyName = "study1",
        context = "brain",
        trait = "ENSG_example",
        entry = list(gr),
        genome = "hg19",
        qcOptions = list(qc = "ok")
    )
    expect_error(estimateLocalH2(noRef), "the row has no LD reference")
    # The overlap factor and the eigenvalue term read the reference's sample
    # size; the null space does not.
    R <- pecotmr:::.ldFromRef(
        ldReference(qtlSumStatsExample)[[1]],
        gr$SNP,
        "test"
    )
    unsized <- .cisH2Example(ref = LdMatrix(R))
    for (m in c("lava", "hdlL")) {
        expect_error(
            estimateLocalH2(unsized, .h2LocalSettings(m)),
            "reads the LD reference's sample size"
        )
    }
    expect_true(is.finite(
        heritabilityEstimate(
            estimateLocalH2(unsized, .h2LocalSettings("fastHer"))
        )
    ))
})

# ---- estimateLocalH2 on a QtlDataset; h2 for the methods starting from it
test_that("estimateLocalH2 gives each QtlDataset trait an in-sample estimate", {
    skip_if_not_installed("PlinkMatrix")
    data(qtlDatasetExample, envir = environment())
    fit <- estimateLocalH2(qtlDatasetExample, cisWindow = 1e6)
    expect_s4_class(fit, "H2Estimate")
    expect_identical(length(fit), 1L)
    expect_identical(h2EstimateParam(fit), H2EstimateParam())
    expect_identical(as.character(fit$study), "study1")
    expect_identical(as.character(fit$context), "brain")
    expect_identical(as.character(fit$trait), "ENSG_example")
    expect_true(is.finite(heritabilityEstimate(fit)))
    # The genotypes are the study's own LD, which the overlap factor records.
    expect_identical(
        estimateLocalH2(
            qtlDatasetExample,
            .h2LocalSettings("lava"),
            cisWindow = 1e6
        )$reference_overlap,
        1
    )
    # The element is the trait's analysis window.
    expect_gt(sum(IRanges::width(fit[[1]])), 0L)
    for (m in c("hdlL", "fastHer")) {
        expect_true(
            is.finite(heritabilityEstimate(
                estimateLocalH2(
                    qtlDatasetExample,
                    .h2LocalSettings(m),
                    cisWindow = 1e6
                )
            )),
            info = m
        )
    }
})

test_that("an in-sample estimate reads a trait's genotypes as their own LD", {
    fx <- .cisH2StrongQtl(13L, h2 = 0.2)
    lava <- .h2LocalSettings("lava")
    fit <- pecotmr:::.h2IndividualFitParam(lava, fx$Xs, fx$y)
    want <- pecotmr:::.h2LocalFit(
        .h2LocalSettings("lava"),
        fx$z,
        rep(fx$n, length(fx$z)),
        list(R = fx$R, nSamples = fx$n, referenceOverlap = 1)
    )
    expect_equal(fit$h2, want$h2, tolerance = 1e-6)
    expect_identical(fit$nSnps, ncol(fx$Xs))
    # A constant column carries nothing and is left out.
    constant <- pecotmr:::.h2IndividualFitParam(lava, cbind(fx$Xs, 1), fx$y)
    expect_equal(constant$h2, fit$h2)
    expect_identical(constant$nSnps, ncol(fx$Xs))
    # A window where nothing varies has no estimate.
    none <- pecotmr:::.h2IndividualFitParam(
        lava,
        matrix(1, 20, 3),
        rnorm(20)
    )
    expect_true(is.na(none$h2))
    expect_identical(none$nSnps, 0L)
})

test_that("REML is qgg's on the window's relationship matrix", {
    skip_if_not_installed("qgg")
    fx <- .cisH2StrongQtl(13L, h2 = 0.4)
    fit <- pecotmr:::.h2IndividualFitParam(
        .h2LocalSettings("greml"),
        fx$Xs,
        fx$y
    )
    Xs <- scale(fx$Xs)
    theta <- as.numeric(
        qgg::greml(y = fx$y, GRM = list(tcrossprod(Xs) / ncol(Xs)))$theta
    )
    expect_equal(fit$h2, theta[[1L]] / sum(theta))
    expect_true(is.finite(fit$h2Se) && fit$h2Se > 0)
    expect_true(is.na(fit$p))
    # GremlOptions() reach qgg::greml(), which cannot run a single iteration.
    expect_error(pecotmr:::.h2IndividualFitParam(
        H2EstimateParam(
            statistic = "phenotype",
            criterion = "reml",
            referenceNoise = "none",
            inference = "information",
            procedures = GremlOptions(maxit = 1)
        ),
        fx$Xs,
        fx$y
    ))
})

test_that("estimateLocalH2 on a QtlDataset refuses what it cannot estimate", {
    data(qtlDatasetExample, envir = environment())
    two <- GenomicRanges::GRanges(c("chr22:1-1000", "chr22:5000-6000"))
    expect_error(
        estimateLocalH2(qtlDatasetExample, region = two),
        "estimate one region at a time"
    )
    expect_error(
        estimateLocalH2(qtlDatasetExample, CredibleSetParam()),
        "H2EstimateParam"
    )
    # No trait with samples in common with its genotypes leaves nothing.
    expect_error(
        pecotmr:::.h2IndividualEstimate(
            list(),
            H2EstimateParam(),
            qtlDatasetExample
        ),
        "no selected trait shares two or more samples"
    )
})

test_that("a supplied or estimated heritability is settled for the methods", {
    fakeFit <- function(settings, value) list(h2 = value)
    resolve <- function(
        h2Estimate = NULL,
        h2EstimateParam = NULL,
        value = 0.4
    ) {
        pecotmr:::.h2Resolve(
            h2Estimate,
            h2EstimateParam,
            "test",
            fakeFit,
            list(value = value),
            individual = TRUE
        )
    }
    lava <- .h2LocalSettings("lava")
    expect_null(resolve())
    expect_identical(resolve(0.3), 0.3)
    # Floored and capped, keeping every engine's variances positive.
    expect_identical(resolve(0), 0.001)
    expect_identical(resolve(0.995), 0.99)
    expect_identical(
        resolve(h2EstimateParam = .h2LocalSettings("greml"), value = 2),
        0.99
    )
    expect_identical(resolve(h2EstimateParam = lava), 0.4)
    expect_error(resolve(1), "must be below 1")
    expect_error(resolve(-0.1))
    expect_error(resolve("lava"), "h2Estimate")
    expect_error(
        resolve(0.3, lava),
        "give the heritability \\(`h2Estimate`\\) or the settings"
    )
    expect_error(
        resolve(h2EstimateParam = list(criterion = "moments")),
        "must be built with H2EstimateParam\\(\\)"
    )
    # A one-row H2Estimate gives its estimate.
    one <- pecotmr:::.h2EstimateNew(
        list(GenomicRanges::GRanges("chr1:1-10")),
        list(
            study = "s",
            context = "c",
            trait = "t",
            h2 = 0.2,
            h2_se = NA_real_,
            h2_lower = NA_real_,
            h2_upper = NA_real_,
            h2_liability = NA_real_,
            h2_liability_se = NA_real_,
            h2_liability_lower = NA_real_,
            h2_liability_upper = NA_real_,
            intercept = NA_real_,
            intercept_se = NA_real_,
            p_value = NA_real_,
            n_snps = 10L,
            n_components = NA_integer_,
            scale = NA_real_,
            reference_overlap = NA_real_
        ),
        pecotmr:::.h2NoObjects(1L),
        GenomicRanges::GRanges(),
        lava
    )
    expect_identical(resolve(one), 0.2)
    two <- pecotmr:::.h2EstimateNew(
        rep(list(GenomicRanges::GRanges("chr1:1-10")), 2L),
        list(
            study = c("s", "s"),
            context = c("c", "c"),
            trait = c("t1", "t2"),
            h2 = c(0.2, 0.3),
            h2_se = NA_real_,
            h2_lower = NA_real_,
            h2_upper = NA_real_,
            h2_liability = NA_real_,
            h2_liability_se = NA_real_,
            h2_liability_lower = NA_real_,
            h2_liability_upper = NA_real_,
            intercept = NA_real_,
            intercept_se = NA_real_,
            p_value = NA_real_,
            n_snps = 10L,
            n_components = NA_integer_,
            scale = NA_real_,
            reference_overlap = NA_real_
        ),
        pecotmr:::.h2NoObjects(2L),
        GenomicRanges::GRanges(),
        lava
    )
    expect_error(resolve(two), "holds one row, not 2")
    # No estimate leaves the method its own default.
    expect_warning(
        expect_null(resolve(h2EstimateParam = lava, value = NA_real_)),
        "no heritability estimate could be made"
    )
    # Summary statistics refuse the phenotype statistic.
    expect_error(
        pecotmr:::.h2Resolve(
            NULL,
            .h2LocalSettings("greml"),
            "test",
            fakeFit,
            list(value = 0.4)
        ),
        "the phenotype statistic needs individual-level data"
    )
})

test_that("an H2Estimate gives a run each selected trait's value, once each", {
    skip_if_not_installed("PlinkMatrix")
    x <- .cisH2Example()
    est <- estimateLocalH2(x, referenceOverlap = 1)
    tuples <- tibble::tibble(
        study = "study1",
        context = "brain",
        trait = "ENSG_example"
    )
    spec <- pecotmr:::.h2PipelineSpec(est, NULL, tuples, "test")
    key <- pecotmr:::.h2TupleKey("study1", "brain", "ENSG_example")
    expect_equal(unname(spec$values), heritabilityEstimate(est))
    expect_equal(
        pecotmr:::.h2ForTrait(spec, key, "test", NULL, list()),
        min(max(heritabilityEstimate(est), 0.001), 0.99)
    )
    # Each trait reads its own row.
    byTrait <- list(values = c(a = 0.2, b = 0.4))
    expect_identical(pecotmr:::.h2ForTrait(byTrait, "a", "t", NULL), 0.2)
    expect_identical(pecotmr:::.h2ForTrait(byTrait, "b", "t", NULL), 0.4)
    # A condition the run derives, such as a principal component, has no row
    # and keeps its default.
    expect_null(pecotmr:::.h2ForTrait(
        spec,
        pecotmr:::.h2TupleKey("study1", "brain", "PC1"),
        "test",
        NULL,
        list()
    ))
    other <- tibble::tibble(
        study = "study1",
        context = "brain",
        trait = "ENSG_other"
    )
    expect_error(
        pecotmr:::.h2PipelineSpec(est, NULL, rbind(tuples, other), "test"),
        "the `h2Estimate` H2Estimate holds no row for 1 selected"
    )
    # One trait's estimates under two settings cannot share an object.
    expect_error(
        c(est, estimateLocalH2(x, .h2LocalSettings("fastHer"))),
        "fitted under different settings"
    )
    fastHer <- .h2LocalSettings("fastHer")
    expect_identical(
        pecotmr:::.h2PipelineSpec(NULL, fastHer, tuples, "t"),
        list(param = fastHer)
    )
    # One number serves every trait.
    one <- pecotmr:::.h2PipelineSpec(0.25, NULL, tuples, "t")
    expect_identical(one, list(value = 0.25))
    expect_identical(pecotmr:::.h2ForTrait(one, key, "t", NULL, list()), 0.25)
    expect_null(pecotmr:::.h2PipelineSpec(NULL, NULL, tuples, "test"))
    expect_error(
        pecotmr:::.h2PipelineSpec(0.25, fastHer, tuples, "t"),
        "not both"
    )
})

# ===========================================================================
# The frequency exponent (alpha)
# ===========================================================================

test_that("the frequency exponent's pseudo-eigenpairs keep the dense model", {
    # Var(r) = h2 R W R / M + c R / n: rotated onto the pseudo-eigenpairs,
    # every direction is independent, and their likelihood is the dense one's
    # up to a constant.
    fx <- .windowLikelihoodFixture(1L, 0.1)
    M <- ncol(fx$R)
    p <- seq(0.05, 0.5, length.out = M)
    w <- pecotmr:::.h2FrequencyWeights(2 * p * (1 - p), 0)
    expect_equal(mean(w), 1)
    dec <- pecotmr:::.h2Truncate(
        pecotmr:::.h2Eigen(list(R = fx$R)),
        H2TruncationParam(),
        fx$n,
        M
    )
    pseudo <- pecotmr:::.h2Weighted(dec, w)
    u <- as.vector(crossprod(dec$vectors, fx$corInput))
    ut <- as.vector(crossprod(pseudo$vectors, fx$corInput))
    lam <- dec$values
    dense <- function(h, c) {
        S <- h *
            crossprod(dec$vectors * w, dec$vectors) *
            tcrossprod(lam) /
            M +
            diag(c * lam / fx$n, length(lam))
        -(as.numeric(determinant(S)$modulus) + sum(u * solve(S, u))) / 2
    }
    diagonal <- function(h, c) {
        v <- h * pseudo$values^2 / M + c * pseudo$values / fx$n
        -sum(log(v) + ut^2 / v) / 2
    }
    grid <- expand.grid(h = c(0.05, 0.2, 0.5), c = c(0.8, 1.2))
    gap <- purrr::pmap_dbl(grid, dense) - purrr::pmap_dbl(grid, diagonal)
    expect_equal(gap - gap[[1L]], rep(0, nrow(grid)), tolerance = 1e-8)
    # Equal frequencies weigh every variant alike.
    expect_equal(pecotmr:::.h2Weighted(dec, rep(1, M))$values, dec$values)
})

test_that("alpha weighs nothing where every frequency is alike", {
    fx <- .windowLikelihoodFixture(1L, 0.1)
    ld <- list(ldFactor = fx$ldFactor, nSamples = 150L, het = rep(0.42, 300))
    settings <- c(
        list(
            lava = .h2LocalSettings("lava"),
            hdlL = .h2LocalSettings("hdlL")
        ),
        .localSettings()
    )
    for (name in names(settings)) {
        perSd <- .localFit(settings[[name]], fx$corInput, fx$n, ld)
        alike <- .localFit(
            methods::initialize(settings[[name]], alpha = 0),
            fx$corInput,
            fx$n,
            ld
        )
        expect_equal(alike$h2, perSd$h2, tolerance = 1e-6, info = name)
    }
})

test_that("alpha reads frequencies where they are recorded, or refuses", {
    lava0 <- methods::initialize(.h2LocalSettings("lava"), alpha = 0)
    # A bare genotype matrix, which may be residualized, records none.
    fx <- .windowLikelihoodFixture(1L, 0.1)
    expect_error(
        pecotmr:::.h2IndividualFitParam(
            lava0,
            fx$ldFactor[1:20, 1:50],
            rnorm(20)
        ),
        "genotype matrix alone"
    )
    # A genome-wide reference that records no frequencies, and carries no
    # panel to read them from, is refused.
    ref <- makeTestEigenRef()
    expect_error(
        estimateH2(
            makeTestSumstatsForRef(ref),
            ref,
            param = methods::initialize(
                .h2MethodSettings("lder", NULL, ref),
                alpha = 0
            )
        ),
        "records none"
    )
})

test_that("a genome-wide block's pseudo-eigenpairs keep its weighted model", {
    # Under the weights W a block's statistics have covariance
    # h2 R W R / M + c R: on its pseudo-eigenvectors the two parts are
    # diagonal, b^2 and b, and a direction's annotation scores on the
    # loadings are its genetic variance per unit rate over b^2.
    fx <- .h2PanelFixture()
    R <- ldMatrix(fx$ld[[1L]])
    p <- nrow(R)
    e <- eigen(R, symmetric = TRUE)
    set.seed(5)
    af <- runif(p, 0.05, 0.5)
    w <- pecotmr:::.h2FrequencyWeights(2 * af * (1 - af), 0)
    block <- pecotmr:::.h2WeightedBlock(
        list(
            values = e$values,
            vectors = e$vectors,
            snpIdx = seq_len(p),
            nNull = 0L
        ),
        w
    )
    V <- block$vectors
    k <- length(block$values)
    expect_equal(k + block$nNull, p)
    expect_equal(crossprod(V, R %*% V), diag(block$values), tolerance = 1e-8)
    expect_equal(
        crossprod(V, R %*% (w * R) %*% V),
        diag(block$values^2),
        tolerance = 1e-8
    )
    a <- rbinom(p, 1, 0.3)
    scores <- pecotmr:::.h2EigenAnnot(block, cbind(w * a, w))
    expect_equal(
        diag(crossprod(V, R %*% (w * a * R) %*% V)),
        block$values^2 * scores[seq_len(k), 1L],
        tolerance = 1e-8
    )
    # The weights' own column scores 1 along every direction.
    expect_equal(scores[seq_len(k), 2L], rep(1, k), tolerance = 1e-8)
})

test_that("a genome-wide reference's frequencies come from its MAF or panel", {
    fx <- .h2PanelFixture()
    lazyEigen <- buildLdEigen(fx$panel, ldBlocks = fx$blocks, genome = "hg19")
    p <- colMeans(
        pecotmr:::.ldRefDosage(fx$panel, seq_len(300), meanImpute = TRUE)
    ) /
        2
    het <- 2 * p * (1 - p)
    # The panel records none, so its dosages are read, a block at a time.
    expect_equal(pecotmr:::.h2RefHeterozygosity(lazyEigen), unname(het))
    # Blocks with an AF column record its MAF, which gives the same.
    stored <- buildLdScore(.h2PanelLdWithAf(fx, p), genome = "hg19")
    expect_equal(S4Vectors::mcols(stored)[["MAF"]], unname(pmin(p, 1 - p)))
    expect_equal(pecotmr:::.h2RefHeterozygosity(stored), unname(het))
    expect_error(
        pecotmr:::.h2RefHeterozygosity(buildLdEigen(fx$ld, genome = "hg19")),
        "records none"
    )
})

test_that("genome-wide alpha weighs nothing where every frequency is alike", {
    fx <- .h2PanelFixture()
    alike <- .h2PanelLdWithAf(fx, rep(0.3, 300))
    storedEigen <- buildLdEigen(alike, genome = "hg19")
    storedScore <- buildLdScore(alike, genome = "hg19")
    annot <- .h2PanelAnnotations(storedScore)
    n <- 2000
    z <- .h2PanelZ(fx$ld, n, h2 = 0.5)
    refs <- list(
        lder = storedEigen,
        hdl = storedEigen,
        gldsc = storedScore,
        sldsc = storedScore
    )
    fields <- c("h2", "h2Se", "intercept", "tau", "tauSe", "enrichment")
    runs <- c(
        map(names(refs), function(m) list(method = m, an = NULL)),
        map(c("lder", "gldsc", "sldsc"), function(m) {
            list(method = m, an = annot)
        })
    )
    for (run in runs) {
        perSd <- .h2GenomeWideSettings(
            run$method,
            .h2Stratified(run$an),
            FALSE
        )
        fit <- function(param) {
            suppressWarnings(pecotmr:::.h2GenomeWideFit(
                param,
                z,
                n,
                refs[[run$method]],
                run$an
            ))[fields]
        }
        expect_equal(
            fit(methods::initialize(perSd, alpha = 0)),
            fit(perSd),
            tolerance = 1e-5,
            label = run$method
        )
    }
    # sHDL's single start stops on a failed line search in one jackknife fold
    # of this small panel, wherever its likelihood's rounding takes it, so
    # its likelihood is compared rather than its folds.
    proc <- .h2GenomeWideSettings("hdl", TRUE, FALSE)@procedures
    A <- pecotmr:::.h2DesignMat(annot)
    w <- pecotmr:::.h2RefWeights(
        methods::initialize(
            .h2GenomeWideSettings("hdl", TRUE, FALSE),
            alpha = 0
        ),
        storedEigen
    )
    eigenDesign <- function(weights) {
        data <- map(
            pecotmr:::.ldEigenBlocks(storedEigen),
            pecotmr:::.h2EigenBlockData,
            eigenRef = storedEigen,
            z = z,
            n = n,
            baselineMat = pecotmr:::.h2WeightedDesign(A, weights),
            candMat = NULL,
            w = weights
        )
        pecotmr:::.h2EigenDesign(
            data,
            length(storedEigen),
            pecotmr:::.h2AnnotationBasis(A, weights),
            proc
        )
    }
    perSdDesign <- eigenDesign(NULL)
    alikeDesign <- eigenDesign(w)
    nll <- function(design, par) {
        pecotmr:::.h2Likelihood(
            "variance",
            pecotmr:::.h2VarianceData(design, n, ncol(fx$panel), proc),
            par
        )$value
    }
    for (par in list(c(0.2, 0.1, 0.05, 1), c(0.4, 0, -0.1, 1.2))) {
        expect_equal(nll(alikeDesign, par), nll(perSdDesign, par))
    }
})

test_that("the stratified report reads heritability under the weights", {
    # Each variant's variance is w_j (A tau)_j: an annotation's heritability
    # is its variants' sum, A'WA tau, and its share of variants is unweighted.
    set.seed(8)
    M <- 200
    A <- structure(
        cbind(base = 1, a = rbinom(M, 1, 0.3), b = runif(M)),
        userCols = 2:3
    )
    af <- runif(M, 0.05, 0.5)
    w <- pecotmr:::.h2FrequencyWeights(2 * af * (1 - af), 0)
    tau <- c(1e-3, 2e-3, 5e-4)
    folds <- matrix(tau, 10, 3, byrow = TRUE) *
        (1 + matrix(rnorm(30, sd = 0.05), 10, 3))
    out <- pecotmr:::.h2StratifiedReport(A, tau, folds, w)$enrichment
    d <- w * as.vector(A %*% tau)
    expect_equal(
        out$propH2,
        c(sum(d * A[, "a"]), sum(d * A[, "b"])) / sum(d)
    )
    expect_equal(out$propSnps, unname(colSums(A[, 2:3]) / M))
    expect_equal(out$enrichment, out$propH2 / out$propSnps)
    # Weights of 1 are none.
    expect_equal(
        pecotmr:::.h2StratifiedReport(A, tau, folds, rep(1, M)),
        pecotmr:::.h2StratifiedReport(A, tau, folds)
    )
})

test_that("estimateH2 weighs a genome-wide reference's variants", {
    fx <- .h2PanelFixture()
    lazyEigen <- buildLdEigen(fx$panel, ldBlocks = fx$blocks, genome = "hg19")
    hdl0 <- methods::initialize(
        .h2GenomeWideSettings("hdl", FALSE, FALSE),
        alpha = 0
    )
    est <- suppressWarnings(
        estimateH2(makeTestSumstatsForRef(lazyEigen), lazyEigen, param = hdl0)
    )
    expect_identical(h2EstimateParam(est)@alpha, 0)
    expect_true(is.finite(heritabilityEstimate(est)))
    # S-LDSC's weighted LD scores need each block's LD.
    scoresOnly <- buildLdScore(
        .h2PanelLdWithAf(fx, rep(0.3, 300)),
        genome = "hg19",
        keepLdMatrices = FALSE
    )
    expect_error(
        pecotmr:::variantRegressionUnivariate(
            rnorm(300),
            2000,
            scoresOnly,
            param = methods::initialize(
                .h2GenomeWideSettings("sldsc", FALSE, FALSE),
                alpha = 0
            )
        ),
        "weighted LD scores"
    )
})

test_that("estimateLocalH2 weighs a panel's and a QtlDataset's variants", {
    skip_if_not_installed("PlinkMatrix")
    hdl0 <- methods::initialize(.h2LocalSettings("hdlL"), alpha = 0)
    x <- .cisH2Example()
    est <- estimateLocalH2(x, hdl0, referenceOverlap = 1)
    expect_identical(h2EstimateParam(est)@alpha, 0)
    expect_true(is.finite(heritabilityEstimate(est)))
    expect_false(isTRUE(all.equal(
        heritabilityEstimate(est),
        heritabilityEstimate(estimateLocalH2(
            x,
            .h2LocalSettings("hdlL"),
            referenceOverlap = 1
        ))
    )))
    data(qtlDatasetExample, envir = environment())
    individual <- estimateLocalH2(
        qtlDatasetExample,
        hdl0,
        contexts = "brain",
        cisWindow = 1e6
    )
    expect_true(all(is.finite(heritabilityEstimate(individual))))
})

test_that("the variant regression's information reads LD between variants", {
    # Its sandwich's middle is the responses' covariance under the model over
    # every pair of variants, 2 S o S with S = (n h2 / M) R R + c R. With flat
    # weights h2's estimate is linear in the responses, so at the model's own
    # S the sandwich is its variance, which draws from that model reproduce.
    # A signal small enough that every draw's sqrt(n) r stays a correlation's.
    fx <- .windowLikelihoodFixture(1L, 0.1)
    M <- ncol(fx$R)
    n <- 300
    S <- n * 0.1 / M * fx$R %*% fx$R + fx$R
    X <- cbind(n * rowSums(fx$R^2) / M, 1)
    rows <- list(
        X = X,
        index = seq_len(M),
        ldModel = list(ld = list(R = fx$R), w = NULL, n = n, M = M)
    )
    meat <- pecotmr:::.h2RegressionMeat(
        rows,
        list(design = list(X = X), w = rep(1, M), coef = c(0.1, 1))
    )
    a <- solve(crossprod(X), t(X))[1L, ]
    variance <- sum(a * ((2 * S^2) %*% a))
    bread <- solve(crossprod(X))
    expect_equal((bread %*% meat %*% bread)[[1L, 1L]], variance)
    e <- eigen(S, symmetric = TRUE)
    root <- e$vectors %*% (sqrt(pmax(e$values, 0)) * t(e$vectors))
    param <- H2EstimateParam(
        statistic = "variant",
        criterion = "regression",
        noise = "free",
        referenceNoise = "none",
        inference = "information",
        procedures = H2RegressionParam(weights = "fixed")
    )
    ld <- list(R = fx$R, nSamples = 150L)
    set.seed(12)
    h2 <- map_dbl(seq_len(400), function(i) {
        z <- pecotmr:::.corToZ(as.vector(root %*% rnorm(M)) / sqrt(n), n)
        pecotmr:::.h2LocalFit(param, z, rep(n, M), ld)$h2
    })
    expect_true(all(is.finite(h2)))
    expect_equal(sd(h2), sqrt(variance), tolerance = 0.1)
})

test_that("HDL-L's interval inverts its likelihood ratio", {
    # Each end of the interval lies where the profile log-likelihood, the
    # intercept refitted, has fallen half the chi-square's 95% point below
    # its value at the fit; settings without the likelihood-ratio test give
    # none.
    fx <- .windowLikelihoodFixture(2L, 0.1)
    ld <- list(ldFactor = fx$ldFactor, nSamples = 150L)
    hdl <- .h2LocalSettings("hdlL")
    fit <- .localFit(hdl, fx$corInput, fx$n, ld)
    expect_true(fit$h2Lower <= fit$h2 && fit$h2 <= fit$h2Upper)
    M <- length(fx$corInput)
    dec <- pecotmr:::.h2Basis(hdl, ld, fx$n, M)
    r <- fx$corInput
    data <- list(
        lam = dec$values,
        bstar = as.vector(crossprod(dec$vectors, r)),
        n = fx$n,
        nRef = 150L,
        M = M,
        floor = hdl@procedures@varianceFloor
    )
    profile <- function(h) {
        pecotmr:::.h2WindowProfile(h, data, hdl@procedures)
    }
    top <- profile(fit$h2)
    for (end in c(fit$h2Lower, fit$h2Upper)) {
        if (end > 0 && end < 1) {
            expect_equal(
                2 * (top - profile(end)),
                qchisq(0.95, 1),
                tolerance = 1e-4
            )
        }
    }
    regression <- .localSettings()$eigenRegression
    none <- .localFit(regression, fx$corInput, fx$n, ld)
    expect_true(is.na(none$h2Lower) && is.na(none$h2Upper))
})

test_that("estimateLocalH2's default reports HDL-L's interval", {
    skip_if_not_installed("PlinkMatrix")
    est <- estimateLocalH2(.cisH2Example())
    md <- S4Vectors::mcols(est)
    expect_true(is.finite(md$h2_lower) && is.finite(md$h2_upper))
    expect_true(md$h2_lower <= md$h2 && md$h2 <= md$h2_upper)
    lava <- estimateLocalH2(
        .cisH2Example(),
        .h2LocalSettings("lava"),
        referenceOverlap = 1
    )
    expect_true(is.na(S4Vectors::mcols(lava)$h2_lower))
})

# ---- edge cases of the local and genome-wide fits ---------------------------

test_that("tied noise without the null space fits the rotated correlations", {
    fx <- .nullSpaceFixture()
    R <- crossprod(fx$X) / fx$n
    tied <- .h2LocalSettings("fastHer")
    tied@procedures@nullSpace <- FALSE
    fit <- .localFit(tied, fx$r, fx$n, list(R = R))
    # The likelihood of b = Q'r, each direction's variance the scale times
    # lambda / n + h2 (lambda^2 / m - lambda / n).
    e <- eigen(R, symmetric = TRUE)
    b <- as.vector(crossprod(e$vectors, fx$r))
    lam <- e$values
    nll <- function(par) {
        v <- par[[2L]] *
            (lam / fx$n + par[[1L]] * (lam^2 / fx$m - lam / fx$n))
        sum(log(v) + b^2 / v) / 2
    }
    opt <- optim(
        c(0.3, 1),
        nll,
        method = "L-BFGS-B",
        lower = c(1e-9, 1e-6),
        upper = c(1 - 1e-9, 10),
        control = list(factr = 1e2, pgtol = 0)
    )
    # Over eight directions the likelihood is flat in h2 and the scale: the
    # fit is its maximum, no lower than L-BFGS-B's, and Nelder-Mead started
    # from it finds nothing higher.
    at <- nll(c(fit$h2, fit$scale))
    expect_lte(at, opt$value + 1e-10)
    polish <- suppressWarnings(optim(
        c(fit$h2, fit$scale),
        nll,
        control = list(reltol = 1e-14, maxit = 5000)
    ))
    expect_gte(polish$value, at - 1e-8)
    expect_equal(fit$nComponents, fx$m)
    # Held at a scale of 1, h2 alone.
    held <- tied
    held@procedures@profileScale <- FALSE
    one <- .localFit(held, fx$r, fx$n, list(R = R))
    nllOne <- function(h2) nll(c(h2, 1))
    want <- optimize(nllOne, c(0, 1 - 1e-9), tol = 1e-10)$minimum
    expect_equal(one$h2, want, tolerance = 1e-6)
    expect_true(is.na(one$scale))
})

test_that("a tied fit that never converges has no estimate", {
    fx <- .nullSpaceFixture()
    R <- crossprod(fx$X) / fx$n
    tied <- .h2LocalSettings("fastHer")
    tied@procedures@nullSpace <- FALSE
    local_mocked_bindings(.h2TiedScoring = function(...) list(best = NULL))
    nullSpace <- .localFit(.h2LocalSettings("fastHer"), fx$r, fx$n, list(R = R))
    plain <- .localFit(tied, fx$r, fx$n, list(R = R))
    expect_true(is.na(nullSpace$h2) && is.na(plain$h2))
    expect_equal(plain$nComponents, fx$m)
})

test_that("a window with too few directions or variants has no estimate", {
    fx <- .nullSpaceFixture()
    R <- crossprod(fx$X) / fx$n
    few <- function(param) {
        methods::initialize(
            param,
            truncation = methods::initialize(param@truncation, minimum = 100)
        )
    }
    regression <- H2EstimateParam(
        criterion = "regression",
        noise = "free",
        referenceNoise = "none",
        truncation = H2TruncationParam(of = "eigenvalues", minimum = 2),
        sampleSize = "mean",
        inference = "information"
    )
    for (param in list(.h2LocalSettings("hdlL"), regression)) {
        fit <- .localFit(few(param), fx$r, fx$n, list(R = R, nSamples = fx$n))
        expect_true(is.na(fit$h2))
        expect_equal(fit$nComponents, fx$m)
    }
    variants <- H2EstimateParam(
        statistic = "variant",
        criterion = "regression",
        noise = "free",
        referenceNoise = "none",
        sampleSize = "median",
        inference = "information"
    )
    two <- .localFit(variants, fx$r[1:2], fx$n, list(R = R[1:2, 1:2]))
    expect_true(is.na(two$h2))
    # HDL-L's interval is the likelihood ratio's, under its test alone.
    hdl <- .h2LocalSettings("hdlL")
    hdl@inference <- "information"
    noTest <- .localFit(hdl, fx$r, fx$n, list(R = R, nSamples = fx$n))
    expect_true(is.finite(noTest$h2))
    expect_true(is.na(noTest$h2Lower) && is.na(noTest$h2Upper))
})

# ---- the eigen likelihoods in compiled code ---------------------------------

# Each one-window likelihood's data and a point inside its bounds, small
# enough to check by hand; the variance model's last two directions sit at
# their floor.
.likelihoodModels <- function() {
    set.seed(3)
    K <- 12L
    lam <- sort(rexp(K), decreasing = TRUE)
    b1 <- rnorm(K, sd = 0.05)
    b2 <- rnorm(K, sd = 0.05)
    genetic <- lam^2 / K
    noise <- lam / 1000
    tied <- list(
        offset = noise,
        g = genetic - noise,
        b2 = b1^2,
        floor = 1e-6
    )
    nullSpace <- list(a = rexp(K), ck = 500 * lam / K - 1, n = 500)
    list(
        variance = list(
            data = list(
                G = cbind(genetic, noise),
                b2 = b1^2,
                offset = c(rep(0, K - 2L), -1, -1),
                floor = 1e-6
            ),
            theta = c(0.3, 1.2)
        ),
        tied = list(data = c(tied, scaled = FALSE), theta = 0.3),
        tiedScaled = list(data = c(tied, scaled = TRUE), theta = c(0.3, 1.1)),
        nullSpace = list(data = c(nullSpace, scaled = FALSE), theta = 0.3),
        nullSpaceScaled = list(
            data = c(nullSpace, scaled = TRUE),
            theta = c(0.3, 0.9)
        ),
        cross = list(
            data = list(
                G = cbind(genetic, noise),
                b1 = b1,
                b2 = b2,
                v11 = 0.3 * genetic + noise,
                v22 = 0.4 * genetic + noise,
                offset = 0,
                floor = 1e-8
            ),
            theta = c(0.1, 0.2)
        ),
        reduced = list(
            data = list(
                s = runif(K),
                a = rnorm(K),
                b = rnorm(K),
                floor = 1e-6
            ),
            theta = 0.4
        )
    )
}

# A likelihood's central difference in each parameter.
.likelihoodSlope <- function(model, data, theta, h = 1e-6) {
    map_dbl(seq_along(theta), function(j) {
        at <- function(x) pecotmr:::.h2Likelihood(model, data, x)$value
        (at(replace(theta, j, theta[[j]] + h)) -
            at(replace(theta, j, theta[[j]] - h))) /
            (2 * h)
    })
}

test_that("every one-window likelihood's score is its slope", {
    models <- .likelihoodModels()
    kind <- c(
        variance = "variance",
        tied = "tied",
        tiedScaled = "tied",
        nullSpace = "nullSpace",
        nullSpaceScaled = "nullSpace",
        cross = "cross",
        reduced = "reduced"
    )
    for (m in names(models)) {
        x <- models[[m]]
        at <- pecotmr:::.h2Likelihood(kind[[m]], x$data, x$theta, TRUE)
        expect_equal(
            at$score,
            .likelihoodSlope(kind[[m]], x$data, x$theta),
            tolerance = 1e-6,
            info = m
        )
        expect_equal(at$information, t(at$information), info = m)
        values <- eigen(at$information, TRUE, only.values = TRUE)$values
        expect_gt(min(values), -1e-10)
    }
    # The variance model's terms, by hand: a direction at its floor carries
    # no parameter.
    x <- models$variance
    raw <- x$data$offset + as.vector(x$data$G %*% x$theta)
    v <- pmax(raw, x$data$floor)
    G <- x$data$G * (raw > x$data$floor)
    at <- pecotmr:::.h2Likelihood("variance", x$data, x$theta, TRUE)
    expect_equal(at$value, -sum(log(v) + x$data$b2 / v) / 2)
    expect_equal(at$information, unname(crossprod(G / v)) / 2)
})

test_that("a fit is the maximum within its bounds, or its free part's", {
    x <- .likelihoodModels()$variance
    lower <- c(0, 0.5)
    upper <- c(1, 10)
    fits <- pecotmr:::.h2Scoring(
        "variance",
        x$data,
        list(c(0.5, 1), c(0.05, 5)),
        lower,
        upper,
        1000
    )
    expect_true(all(map_lgl(fits, "converged")))
    best <- pecotmr:::.h2ScoringBest(fits)
    score <- pecotmr:::.h2Likelihood("variance", x$data, best$par, TRUE)$score
    inside <- best$par > lower & best$par < upper
    expect_lt(max(abs(score[inside])), 1e-4)
    expect_true(all(score[best$par == lower] <= 0))
    expect_true(all(score[best$par == upper] >= 0))
    # With h2 held, the intercept alone, and its value the full model's.
    held <- pecotmr:::.h2Scoring(
        "variance",
        x$data,
        list(1),
        0.5,
        10,
        1000,
        full = c(0.3, 1),
        free = 2L
    )[[1L]]
    expect_length(held$par, 1L)
    expect_equal(dim(held$information), c(1L, 1L))
    expect_equal(
        held$value,
        pecotmr:::.h2Likelihood("variance", x$data, c(0.3, held$par))$value
    )
    expect_gte(best$value, held$value)
})

test_that("a fit that runs out of steps or reads no finite value fails", {
    models <- .likelihoodModels()
    x <- models$variance
    short <- pecotmr:::.h2Scoring(
        "variance",
        x$data,
        list(c(0.9, 9)),
        c(0, 0.5),
        c(1, 10),
        1
    )
    expect_false(short[[1L]]$converged)
    # At h2 = 1 the null space's variance is 0: no finite log-likelihood.
    y <- models$nullSpace
    edge <- pecotmr:::.h2Scoring("nullSpace", y$data, list(1), 0, 1, 100)
    expect_false(edge[[1L]]$converged)
    expect_equal(edge[[1L]]$value, -Inf)
})

test_that("a parameter no direction informs takes no step", {
    # Every direction at its floor: a flat likelihood, read from a design
    # given as a vector.
    flat <- list(G = rep(1, 3), b2 = c(1, 2, 3), offset = -10, floor = 0.1)
    fit <- pecotmr:::.h2Scoring("variance", flat, list(0.5), 0, 1, 100)
    expect_true(fit[[1L]]$converged)
    expect_equal(fit[[1L]]$par, 0.5)
})

test_that("the compiled likelihoods refuse data they cannot read", {
    x <- .likelihoodModels()$variance
    expect_error(
        pecotmr:::.h2Likelihood("other", x$data, x$theta),
        "Unknown likelihood"
    )
    expect_error(
        pecotmr:::.h2Likelihood("variance", list(1, 2), 0),
        "has no 'G'"
    )
    expect_error(
        pecotmr:::.h2Likelihood("variance", x$data[c("G", "b2", "offset")], 0),
        "has no 'floor'"
    )
    expect_error(
        pecotmr:::.h2Scoring(
            "variance",
            x$data,
            list(c(0.1, 1, 2)),
            c(0, 0, 0),
            c(1, 1, 1),
            10
        ),
        "incompatible"
    )
})

test_that("the fits' small helpers take each kind of input", {
    # Settings without the eigenvalue term read no reference size.
    hdl <- .h2GenomeWideSettings("hdl", FALSE, FALSE)
    hdl@referenceNoise <- "none"
    expect_identical(pecotmr:::.h2PanelSize(hdl, NULL), Inf)
    # g-LDSC's S from eigenpairs is the one from their matrix.
    fx <- .nullSpaceFixture()
    R <- crossprod(fx$X) / fx$n
    e <- eigen(R, symmetric = TRUE)
    d <- seq(0.5, 1.5, length.out = nrow(R))
    expect_equal(
        pecotmr:::.h2GlsS(list(vectors = e$vectors, values = e$values), d, 2),
        pecotmr:::.h2GlsS(list(R = R), d, 2)
    )
    # REML's genotypes carry the frequency weights as column scales.
    het <- c(0.1, 0.2, 0.3, 0.4, 0.42, 0.45, 0.48, 0.5)
    alpha <- methods::initialize(.h2LocalSettings("greml"), alpha = 0)
    w <- pecotmr:::.h2FrequencyWeights(het, 0)
    expect_equal(
        pecotmr:::.h2WeightedColumns(fx$X, alpha, list(het = het)),
        sweep(fx$X, 2L, sqrt(w), `*`)
    )
})

test_that("a genome-wide HDL fit with no converged start has no estimate", {
    eigenRef <- makeTestEigenRef()
    local_mocked_bindings(.h2ScoringBest = function(...) NULL)
    est <- suppressWarnings(
        estimateH2(makeTestSumstatsForRef(eigenRef), eigenRef)
    )
    expect_true(is.na(heritabilityEstimate(est)))
})

test_that("a panel that records its frequencies gives them to the weights", {
    fx <- .h2PanelFixture()
    lazyEigen <- buildLdEigen(fx$panel, ldBlocks = fx$blocks, genome = "hg19")
    # The reference records no MAF, so its panel's frequencies are read where
    # the panel records them (here 0.3 for every variant, which no dosages of
    # it give), its dosages only where it does not.
    local_mocked_bindings(.ldRefAf = function(x) rep(0.3, nrow(x)))
    expect_equal(
        pecotmr:::.h2RefHeterozygosity(lazyEigen),
        rep(2 * 0.3 * 0.7, length(lazyEigen))
    )
})

test_that("estimateLocalH2 on a MultiStudyQtlDataset is each study's own", {
    skip_if_not_installed("PlinkMatrix")
    ms <- .twoStudyFixture()
    got <- estimateLocalH2(ms, cisWindow = 1e6)
    own <- c(
        estimateLocalH2(as.list(ms)$s1, cisWindow = 1e6),
        estimateLocalH2(as.list(ms)$s2, cisWindow = 1e6)
    )
    expect_identical(got$study, own$study)
    expect_identical(got$context, own$context)
    expect_identical(heritabilityEstimate(got), heritabilityEstimate(own))
    # Its summary statistics' rows are estimated as a QtlSumStats's.
    withSummary <- MultiStudyQtlDataset(
        qtlDatasets = as.list(ms),
        sumStats = .summaryStudy(ms)
    )
    all <- estimateLocalH2(withSummary, cisWindow = 1e6, referenceOverlap = 1)
    rows <- estimateLocalH2(
        sumStats(withSummary),
        cisWindow = 1e6,
        referenceOverlap = 1
    )
    expect_identical(all$study, c(own$study, "s3"))
    expect_identical(
        heritabilityEstimate(all),
        c(heritabilityEstimate(own), heritabilityEstimate(rows))
    )
    # A context named selects within the studies that hold it.
    onlyB <- estimateLocalH2(withSummary, cisWindow = 1e6, contexts = "b")
    expect_identical(onlyB$study, c("s1", "s2"))
    expect_identical(onlyB$context, c("b", "b"))
    expect_error(
        estimateLocalH2(ms, cisWindow = 1e6, traitId = "none"),
        "no study holds a selected context and trait"
    )
})

test_that("a MultiStudyQtlDataset's REML needs no summary rows selected", {
    skip_if_not_installed("PlinkMatrix")
    skip_if_not_installed("qgg")
    ms <- .twoStudyFixture()
    withSummary <- MultiStudyQtlDataset(
        qtlDatasets = as.list(ms),
        sumStats = .summaryStudy(ms)
    )
    greml <- .h2LocalSettings("greml")
    expect_error(
        estimateLocalH2(withSummary, greml, cisWindow = 1e6),
        "needs individual-level data"
    )
    # The summary study holds only context a, so b selects none of it.
    onlyB <- estimateLocalH2(
        withSummary,
        greml,
        cisWindow = 1e6,
        contexts = "b"
    )
    expect_identical(onlyB$study, c("s1", "s2"))
    expect_true(all(is.finite(heritabilityEstimate(onlyB))))
})
