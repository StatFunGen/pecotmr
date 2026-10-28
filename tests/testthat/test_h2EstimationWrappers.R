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
# .validateMethodRef (internal)
# ===========================================================================

test_that(".validateMethodRef errors when lder or hdl is given an LdScore", {
    score_ref <- makeTestScoreRef()
    for (method in c("lder", "hdl")) {
        expect_error(
            pecotmr:::.validateMethodRef(method, score_ref),
            "requires an LdEigen",
            info = method
        )
    }
})

test_that(".validateMethodRef errors when method='gldsc' but ref is LdEigen", {
    eigen_ref <- makeTestEigenRef()
    expect_error(
        pecotmr:::.validateMethodRef("gldsc", eigen_ref),
        "requires an LdScore"
    )
})

test_that(".validateMethodRef returns TRUE for valid combinations", {
    eigen_ref <- makeTestEigenRef()
    score_ref <- makeTestScoreRef()

    expect_true(pecotmr:::.validateMethodRef("lder", eigen_ref))
    expect_true(pecotmr:::.validateMethodRef("hdl", eigen_ref))
    expect_true(pecotmr:::.validateMethodRef("gldsc", score_ref))
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

test_that("computeLdScores LdEigen with annotations returns base + annotation columns", {
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
        lder = list(pecotmr:::lderUnivariate, lazyEigen, storedEigen),
        gldsc = list(pecotmr:::gldscUnivariate, lazyScore, storedScore),
        sldsc = list(pecotmr:::sldscUnivariate, lazyScore, storedScore)
    )
    for (method in names(runs)) {
        fn <- runs[[method]][[1L]]
        for (an in list(NULL, annot)) {
            lazy <- fn(z, n, runs[[method]][[2L]], an, TRUE)
            stored <- fn(z, n, runs[[method]][[3L]], an, TRUE)
            expect_equal(lazy, stored, tolerance = 1e-8, label = method)
        }
    }
    # HDL's likelihood is the same function either way. Its fit then agrees
    # to the optimizer's tolerance: L-BFGS-B turns the likelihood's rounding
    # (1e-13) into differences near 1e-7 in the optimum.
    hdlDesign <- function(ref) {
        data <- map(
            pecotmr:::.ldEigenBlocks(ref),
            pecotmr:::.hdlBlockData,
            eigenRef = ref,
            z = z,
            n = n,
            baselineMat = NULL,
            candMat = NULL
        )
        pecotmr:::.hdlDesign(data, length(ref), NULL)
    }
    lazyDesign <- hdlDesign(lazyEigen)
    storedDesign <- hdlDesign(storedEigen)
    expect_equal(length(lazyDesign$lam), length(storedDesign$lam))
    for (par in list(c(0.3, 1), c(0.5, 1.2))) {
        expect_equal(
            pecotmr:::.hdlNll(
                par,
                n,
                100,
                lazyDesign$lam,
                lazyDesign$bstar,
                lazyDesign$ldAnnotScaled
            ),
            pecotmr:::.hdlNll(
                par,
                n,
                100,
                storedDesign$lam,
                storedDesign$bstar,
                storedDesign$ldAnnotScaled
            ),
            tolerance = 1e-10
        )
    }
    hdl <- pecotmr:::hdlUnivariate
    for (an in list(NULL, annot)) {
        expect_equal(
            suppressWarnings(hdl(z, n, lazyEigen, an, TRUE)),
            suppressWarnings(hdl(z, n, storedEigen, an, TRUE)),
            tolerance = 1e-6
        )
    }
    # The truncated decompositions HDL conventionally uses agree too.
    lazy9 <- buildLdEigen(
        fx$panel,
        ldBlocks = fx$blocks,
        genome = "hg19",
        eigenvalueTruncation = 0.9
    )
    stored9 <- buildLdEigen(fx$ld, genome = "hg19", eigenvalueTruncation = 0.9)
    expect_equal(
        suppressWarnings(hdl(z, n, lazy9, annot, TRUE)),
        suppressWarnings(hdl(z, n, stored9, annot, TRUE)),
        tolerance = 1e-6
    )
})

# ===========================================================================
# computeLdScores for LdScore
# ===========================================================================

test_that("computeLdScores LdScore without annotations returns stored ld_scores", {
    score_ref <- makeTestScoreRef()
    result <- computeLdScores(score_ref)
    expect_identical(result, mcols(score_ref)$ldScores)
})

test_that("computeLdScores LdScore with annotations but no ld_matrix_list errors", {
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

test_that("h2EstimateToSldscTrait creates 1-row dummy when annotationJackknifeCoefs is NULL", {
    h2_obj <- makeTestH2Estimate(with_enrichment = TRUE)
    h2_obj@annotationJackknifeCoefs <- NULL
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

makeTestGwasGenotypeHandle <- function() {
    new(
        "GenotypeHandle",
        path = "/tmp/test",
        format = "plink2",
        snpInfo = data.frame(),
        nSamples = 0L,
        sampleIds = character(),
        pgenPtr = NULL
    )
}

.dfToGwasGr <- function(df) {
    gr <- GenomicRanges::GRanges(
        seqnames = df$CHR,
        ranges = IRanges::IRanges(start = df$BP, width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = df$SNP,
        A1 = df$A1,
        A2 = df$A2,
        Z = df$Z,
        N = df$N
    )
    gr
}

makeTestSumstatsForRef <- function(
    ref,
    traitName = "test",
    varY = NA_real_
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
        varY = varY
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
        result <- estimateH2(ss, refs[[method]], method = method)

        expect_true(is(result, "H2Estimate"), info = method)
        expect_true(is.numeric(heritabilityEstimate(result)), info = method)
        expect_true(is.numeric(heritabilityEstimateSe(result)), info = method)
        expect_true(is.numeric(heritabilityIntercept(result)), info = method)
        expect_true(is.numeric(heritabilityInterceptSe(result)), info = method)
        expect_equal(methodNames(result), method, info = method)
        expect_equal(nSnps(result), nSnps(ss), info = method)
        expect_equal(traitNames(result), "test", info = method)
    }
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
        estimateH2(short, ref, method = "lder"),
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
        estimateH2(shifted, ref, method = "lder"),
        "disagree at 20 of 20 position"
    )
})

test_that("estimateH2 rejects a sumstats holding the right variants in the wrong order", {
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
        estimateH2(reordered, ref, method = "lder"),
        "same variants in the same order"
    )
})

test_that("estimateH2 rejects annotations that do not cover the reference", {
    ref <- makeTestEigenRef()
    ss <- makeTestSumstatsForRef(ref)
    annot <- makeTestAnnotations(nSnps = 10)
    expect_error(
        estimateH2(ss, ref, method = "lder", annotations = annot),
        "`annotations` carries 10 variant\\(s\\)"
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
    expect_error(estimateH2(hg38, ref, method = "lder"), "build mismatch")
})

test_that("estimateH2 accepts an LD reference that records no genome build", {
    ref <- makeTestEigenRef()
    ss <- makeTestSumstatsForRef(ref)
    # An LD reference often records no build, so one built from it
    # names none. An unknown build must not read as a mismatch.
    GenomeInfoDb::genome(ref) <- NA_character_
    expect_s4_class(estimateH2(ss, ref, method = "lder"), "H2Estimate")
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
        estimateH2(ss, noMat, method = "gldsc"),
        "g-LDSC needs each LD block's LD"
    )
    res <- estimateH2(ss, noMat, method = "sldsc")
    expect_s4_class(res, "H2Estimate")
    expect_equal(methodNames(res), "sldsc")
    expect_true(is.finite(heritabilityEstimate(res)))
})

test_that("stratified S-LDSC needs the per-block LD matrices", {
    noMat <- makeTestScoreRef(with_ld_matrices = FALSE)
    ss <- makeTestSumstatsForRef(noMat)
    annot <- makeTestAnnotations(nSnps = length(noMat))
    expect_error(
        estimateH2(ss, noMat, method = "sldsc", annotations = annot),
        "stratified S-LDSC needs each LD block's LD"
    )
})

test_that("S-LDSC reports an enrichment table when stratified", {
    ref <- makeTestScoreRef(with_ld_matrices = TRUE)
    ss <- makeTestSumstatsForRef(ref)
    annot <- makeTestAnnotations(nSnps = length(ref))
    res <- estimateH2(ss, ref, method = "sldsc", annotations = annot)

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
    expect_false(is.null(annotationJackknifeCoefs(res)))
})

test_that("S-LDSC requires an LdScore, not an LdEigen", {
    eigenRef <- makeTestEigenRef()
    ss <- makeTestSumstatsForRef(eigenRef)
    expect_error(
        estimateH2(ss, eigenRef, method = "sldsc"),
        "requires an LdScore"
    )
})

test_that("the S-LDSC weights reproduce the original Hsq.weights", {
    # het_w = 1 / (2 (intercept + hsq*N/M * ld)^2), oc_w = 1 / w_ld, both
    # LD scores floored at 1 and hsq clamped to [0, 1].
    w <- pecotmr:::.sldscWeights(
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
        pecotmr:::.sldscWeights(2, 1, 1000, 100, hsq = 5, intercept = 1),
        pecotmr:::.sldscWeights(2, 1, 1000, 100, hsq = 1, intercept = 1)
    )
    expect_equal(
        pecotmr:::.sldscWeights(2, 1, 1000, 100, hsq = -3, intercept = 1),
        pecotmr:::.sldscWeights(2, 1, 1000, 100, hsq = 0, intercept = 1)
    )
})

# ===========================================================================
# Rank-deficient LD
#
# A real LD block is routinely rank-deficient, and its decomposition then
# carries a tail of tiny negative eigenvalues (order -1e-15) that are pure
# floating-point noise. Those used to reach .lderWeights() as negative
# regression weights, so IRWLS took sqrt() of a negative and lm.fit aborted
# with "NA/NaN/Inf in 'x'".
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
    res <- estimateH2(ss, ref, method = "lder")
    expect_s4_class(res, "H2Estimate")
    expect_false(is.na(heritabilityEstimate(res)))
})

test_that(".lderWeights never returns a negative weight", {
    lam <- c(2, 1, 0.5, 0, -1e-15, -3e-15)
    w <- pecotmr:::.lderWeights(
        lam,
        matrix(1, length(lam), 1),
        slopes = 0.1,
        N = 1e5,
        M = 1000,
        rough = FALSE
    )
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

    expect_s4_class(estimateH2(ss, eigenRef, method = "lder"), "H2Estimate")
    suppressWarnings(
        expect_s4_class(estimateH2(ss, eigenRef, method = "hdl"), "H2Estimate")
    )
    expect_s4_class(estimateH2(ss, scoreRef, method = "gldsc"), "H2Estimate")
})

test_that("estimateH2 rejects a single-block LD reference", {
    # One block used to fail inside the jackknife fold fit with "'x' must be
    # a matrix" (LDER, g-LDSC) or return an SE of exactly 0 (HDL).
    ld <- makeTestLdMatrix(n = 30L, rho = 0.6)
    eigenRef <- buildLdEigen(ld, genome = "hg19")
    ss <- makeTestSumstatsForRef(eigenRef)
    expect_error(
        estimateH2(ss, eigenRef, method = "lder"),
        "delete-one-block jackknife"
    )
    expect_error(
        estimateH2(ss, buildLdScore(ld, genome = "hg19"), method = "gldsc"),
        "at least two are required"
    )
})

test_that("buildLdScore(keepLdMatrices = FALSE) is rejected by g-LDSC", {
    ld <- makeTestLdMatrixBlocks(sizes = rep(10L, 5))
    scoreRef <- buildLdScore(ld, genome = "hg19", keepLdMatrices = FALSE)
    ss <- makeTestSumstatsForRef(scoreRef)
    expect_error(
        estimateH2(ss, scoreRef, method = "gldsc"),
        "g-LDSC needs each LD block's LD"
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
    expect_s4_class(estimateH2(ss, ref, method = "lder"), "H2Estimate")
})

test_that("estimateH2 with a varY correction runs without error", {
    eigen_ref <- makeTestEigenRef()
    ss <- makeTestSumstatsForRef(
        eigen_ref,
        traitName = "cc_trait",
        varY = 4.0
    )
    expect_equal(varY(ss), 4.0)

    result <- estimateH2(ss, eigen_ref, method = "lder")
    expect_s4_class(result, "H2Estimate")
    expect_true(is.numeric(heritabilityEstimate(result)))
})

# ===========================================================================
# computeLdScores for LdScore WITH annotations
# ===========================================================================

test_that("computeLdScores LdScore with annotations and ld_matrix_list returns correct matrix", {
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
# computeBaselineEnrichment
# =============================================================================

test_that("computeBaselineEnrichment computes enrichment and propSnps", {
    M <- 100
    baseline_mat <- matrix(0, nrow = M, ncol = 2)
    baseline_mat[1:50, 1] <- 1
    baseline_mat[, 2] <- 1
    tau <- c(0.01, 0.005)
    tauSe <- c(0.002, 0.001)
    h2 <- 0.5
    annot_names <- c("annot1", "base")
    res <- pecotmr:::computeBaselineEnrichment(
        tau,
        tauSe,
        NULL,
        baseline_mat,
        annot_names,
        h2
    )
    expect_s3_class(res, "data.frame")
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
    expect_named(res, expected_cols)
    expect_equal(nrow(res), 2)
    expect_equal(res$enrichment, tau * M / h2)
    expect_equal(res$propSnps, c(0.5, 1.0))
})

test_that("computeBaselineEnrichment uses jackknife blocks when provided", {
    M <- 100
    nBlocks <- 5
    baseline_mat <- matrix(0, nrow = M, ncol = 2)
    baseline_mat[1:50, 1] <- 1
    baseline_mat[, 2] <- 1
    tau <- c(0.01, 0.005)
    tauSe <- c(0.002, 0.001)
    h2 <- 0.5
    annot_names <- c("annot1", "base")
    annotationJackknifeCoefs <- matrix(
        rep(tau, each = nBlocks),
        nrow = nBlocks,
        ncol = 2
    )
    res <- pecotmr:::computeBaselineEnrichment(
        tau,
        tauSe,
        annotationJackknifeCoefs,
        baseline_mat,
        annot_names,
        h2
    )
    # With constant annotationJackknifeCoefs, enrichmentSe should be 0
    expect_equal(res$enrichmentSe, c(0, 0))
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

# Tests for the heritability estimators lderUnivariate, gldscUnivariate and
# hdlUnivariate.

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
        cand <- exec(
            cbind,
            !!!map(
                seq_len(nCand),
                function(k) as.numeric(rbinom(n_snps, 1, 0.3 + 0.1 * k))
            )
        )
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

# The three estimators, each with the reference it reads and the baseline
# its annotations take. An all-ones baseline makes the LDER and HDL
# eigenvalue-score column the eigenvalue itself; for g-LDSC it would repeat
# the base LD score column and leave the design singular, so g-LDSC takes a
# half-on/half-off baseline instead.
.h2Runs <- list(
    lder = list(
        fn = pecotmr:::lderUnivariate,
        ref = "eigen_ref",
        baseline = "ones"
    ),
    gldsc = list(
        fn = pecotmr:::gldscUnivariate,
        ref = "ld_score_ref",
        baseline = "spatial"
    ),
    hdl = list(
        fn = pecotmr:::hdlUnivariate,
        ref = "eigen_ref",
        baseline = "ones"
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
        "local",
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

test_that("lderUnivariate intercept is near zero for well-calibrated data", {
    # In LDER the intercept parameter a represents confounding deviation:
    # the model is E[chi2Rot - 1] = n * h2/M * d + n * a
    # so a ~ 0 when there is no confounding.
    res <- pecotmr:::lderUnivariate(dat$z, dat$n, dat$eigen_ref)
    expect_true(is.finite(res$intercept))
    expect_true(abs(res$intercept) < 0.01)
})

test_that("each univariate estimator with local = TRUE returns per-block h2", {
    for (method in names(.h2Runs)) {
        res <- .h2Run(method, dat, local = TRUE)
        expect_true(is.data.frame(res$local), info = method)
        expect_true(
            all(c("blockId", "h2Local", "h2LocalSe") %in% colnames(res$local)),
            info = method
        )
    }
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

test_that("each estimator with annotations and local = TRUE gives local h2", {
    for (method in names(.h2Runs)) {
        annot <- .h2RunAnnot(method, dat_annot, 1L)
        res <- .h2Run(method, dat_annot, annotations = annot, local = TRUE)
        expect_true(is.data.frame(res$local), info = method)
        expect_true(
            all(c("blockId", "h2Local", "h2LocalSe") %in% colnames(res$local)),
            info = method
        )
        # One row per block
        expect_equal(nrow(res$local), 5L, info = method)
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
# h2EstimateToSldscTrait — annotationJackknifeCoefs present but without column names
# ---------------------------------------------------------------------------

test_that("h2EstimateToSldscTrait assigns category names to unnamed annotationJackknifeCoefs", {
    h2_obj <- makeTestH2Estimate(with_enrichment = TRUE)
    tb <- annotationJackknifeCoefs(h2_obj)
    colnames(tb) <- NULL
    h2_obj@annotationJackknifeCoefs <- tb
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
test_that("hdlUnivariate recovers a known h2 (deterministic, seed = 1)", {
    d <- .h2RecoveryRefs(nBlocks = 10L)
    uni <- rep(0.4 / d$M, d$M)
    set.seed(1)
    ests <- replicate(30, hdlUnivariate(d$simZ(uni), d$N, d$eigenRef)$h2)
    expect_lt(abs(mean(ests) - 0.4), 0.05) # recovers true h2 = 0.4
    expect_equal(mean(ests), 0.395, tolerance = 0.04) # stored known output
})

test_that("lderUnivariate recovers a known h2 (deterministic, seed = 1)", {
    d <- .h2RecoveryRefs()
    uni <- rep(0.4 / d$M, d$M)
    set.seed(1)
    ests <- replicate(30, lderUnivariate(d$simZ(uni), d$N, d$eigenRef)$h2)
    expect_lt(abs(mean(ests) - 0.4), 0.05)
    expect_equal(mean(ests), 0.389, tolerance = 0.04)
})

test_that("gldscUnivariate recovers a known h2 (deterministic, seed = 1)", {
    d <- .h2RecoveryRefs()
    uni <- rep(0.4 / d$M, d$M)
    set.seed(1)
    ests <- replicate(30, gldscUnivariate(d$simZ(uni), d$N, d$scoreRef)$h2)
    expect_lt(abs(mean(ests) - 0.4), 0.1) # g-LDSC has higher finite-block var
    expect_equal(mean(ests), 0.365, tolerance = 0.08)
})

test_that("hdlUnivariate recovers partitioned h2 (block-separable annots)", {
    d <- .h2RecoveryRefs(nBlocks = 10L)
    annot <- .h2RecoveryAnnot(d$M, d$M / 2)
    perSnp <- ifelse(seq_len(d$M) <= d$M / 2, 0.6 / d$M, 0.2 / d$M)
    M_a <- c(d$M / 2, d$M / 2)
    set.seed(1)
    h2 <- rowMeans(replicate(15, {
        r <- hdlUnivariate(d$simZ(perSnp), d$N, d$eigenRef, annotations = annot)
        r$enrichment$tau * M_a
    }))
    expect_lt(abs(h2[1] - 0.3), 0.06) # annotA true h2 = 0.3
    expect_lt(abs(h2[2] - 0.1), 0.05) # annotB true h2 = 0.1
})

test_that("gldscUnivariate recovers enrichment direction (annotA > annotB)", {
    d <- .h2RecoveryRefs()
    annot <- .h2RecoveryAnnot(d$M, d$M / 2)
    perSnp <- ifelse(seq_len(d$M) <= d$M / 2, 0.6 / d$M, 0.2 / d$M)
    set.seed(1)
    e <- rowMeans(replicate(20, {
        gldscUnivariate(
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
        lderUnivariate(z, d$N, refOut)$h2,
        lderUnivariate(z, d$N, refIn)$h2
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
test_that("sldscUnivariate lands near a known h2 on a small reference", {
    d <- .h2RecoveryRefs(N = 5000)
    uni <- rep(0.4 / d$M, d$M)
    set.seed(1)
    ests <- replicate(
        20,
        pecotmr:::sldscUnivariate(d$simZ(uni), d$N, d$scoreRef)$h2
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
        pecotmr:::sldscUnivariate(
            d$simZ(rep(0.1 / d$M, d$M)),
            d$N,
            d$scoreRef
        )$h2
    )
    set.seed(1)
    high <- replicate(
        10,
        pecotmr:::sldscUnivariate(
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
    res <- pecotmr:::sldscUnivariate(
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
        pecotmr:::sldscUnivariate(
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
# .sldscWlsCoef(), which the original has no equivalent of. Case 2's
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
            pecotmr:::.ldscUnbiasedR2(b$R^2, nSamples)
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
        fit <- pecotmr:::.sldscFit(
            z^2,
            scores,
            baseScore,
            g$n,
            matrix(1, M, 1),
            2L
        )
        expect_equal(fit$h2, g$h2Hat, tolerance = 1e-5)
        expect_equal(fit$intercept, g$int, tolerance = 1e-5)
    }
})

test_that("partitioned S-LDSC reproduces the original ldsc to 1e-5", {
    ref <- .h2ParityRef()$score
    M <- length(ref)
    inA <- as.integer(seq_len(M) <= M / 2)
    baselineMat <- matrix(inA, ncol = 1, dimnames = list(NULL, "annotA"))
    scores <- pecotmr:::.sldscScoreMatrix(ref, baselineMat, M)
    A <- cbind(base = 1, annotA = inA)

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
        fit <- pecotmr:::.sldscFit(z^2, scores, scores[, 1], g$n, A, 2L)
        expect_equal(fit$h2, g$h2Hat, tolerance = 1e-5)
        expect_equal(fit$intercept, g$int, tolerance = 1e-5)
    }
})

test_that("the two S-LDSC routes differ, as the original ldsc's do", {
    # Univariate takes the two-step estimator; partitioned takes the
    # chi2-filtered single-pass `old_weights` route. The original ldsc refuses
    # to run two-step on a partitioned design at all, so these must not
    # converge to one shared code path.
    ref <- .h2ParityRef()$score
    M <- length(ref)
    z <- .h2ParityZ(0.4, 10000, 42)
    scores <- matrix(as.vector(mcols(ref)$ldScores[, 1]), ncol = 1)
    uni <- pecotmr:::.sldscFit(
        z^2,
        scores,
        scores[, 1],
        10000,
        matrix(1, M, 1),
        2L
    )
    inA <- as.integer(seq_len(M) <= M / 2)
    part <- pecotmr:::.sldscFit(
        z^2,
        pecotmr:::.sldscScoreMatrix(
            ref,
            matrix(inA, ncol = 1, dimnames = list(NULL, "annotA")),
            M
        ),
        scores[, 1],
        10000,
        cbind(base = 1, annotA = inA),
        2L
    )
    expect_false(isTRUE(all.equal(uni$intercept, part$intercept)))
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
        est <- pecotmr:::hdlUnivariate(z, 1e5, small)$h2,
        "sits at the top of its range"
    )

    big <- .h2ParityRef(nSamples = 20000L)$eigen
    expect_no_warning(
        ok <- pecotmr:::hdlUnivariate(z, 1e5, big)$h2,
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
                pecotmr:::hdlUnivariate(.h2ParityZ(0.4, 1e5, seed), 1e5, ref)$h2
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

test_that("HDL reproduces the original HDL", {
    reference <- list(
        list(h2 = 0.4, n = 1e5, seed = 2, nRef = 20000L, want = 0.405111844),
        list(h2 = 0.3, n = 1e4, seed = 5, nRef = 20000L, want = 0.286758612),
        list(h2 = 0.6, n = 5e4, seed = 9, nRef = 50000L, want = 0.646809004)
    )
    for (g in reference) {
        ref <- .h2ParityRef(nSamples = g$nRef)$eigen
        z <- .h2ParityZ(g$h2, g$n, g$seed)
        est <- suppressWarnings(pecotmr:::hdlUnivariate(z, g$n, ref)$h2)
        expect_equal(est, g$want, tolerance = 1e-6)
    }
})

test_that("the HDL variance floor matches the original's exp(-10)", {
    # The original llfun floors the modelled per-direction variance at
    # exp(-10). On a real reference the floor is active for a large share of
    # the eigen-directions, and the bstar^2/lamh2 term it feeds dominates the
    # likelihood -- a different floor silently moves the optimum.
    expect_equal(formals(pecotmr:::.hdlNll)$lim, quote(exp(-10)))
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
        est <- pecotmr:::lderUnivariate(z, g$n, ref)$h2
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
        est <- pecotmr:::gldscUnivariate(z, g$n, ref)$h2
        expect_equal(est, g$want, tolerance = 1e-6)
    }
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
    A <- matrix(1, 6L, 1L)
    sc <- matrix(1, 6L, 1L)
    out <- pecotmr:::.sldscFitUnivariate(
        rep(100, 6L),
        sc,
        rep(1, 6L),
        1000,
        A,
        2L
    )
    expect_setequal(
        names(out),
        c("tau", "h2", "estH", "intercept", "coef")
    )
    expect_false(is.na(out$h2))
})


# ===========================================================================
# HDL bound diagnostics and the partitioned chi2 filter
# ===========================================================================

test_that("HDL warns only when a parameter sits at its UPPER bound", {
    # Only upper bounds are diagnostic: an intercept resting at 0 just means
    # no confounding and is routine, so a lower bound must not warn.
    f <- pecotmr:::.hdlWarnIfAtBound
    expect_silent(f(c(0.1, 0.5), c(1, 2), 1000, 500))
    # At a lower bound: still silent.
    expect_silent(f(c(0, 0), c(1, 2), 1000, 500))
})

test_that("the HDL bound warning names which parameter pinned", {
    # The last parameter is the intercept; anything else is heritability.
    f <- pecotmr:::.hdlWarnIfAtBound
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
    A <- matrix(1, 6L, 1L)
    sc <- matrix(1, 6L, 2L)
    out <- pecotmr:::.sldscFitPartitioned(
        rep(1e6, 6L),
        sc,
        rep(1, 6L),
        1000,
        A
    )
    expect_setequal(
        names(out),
        c("tau", "h2", "estH", "intercept", "coef")
    )
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

test_that(".lderGetRes(twostage = FALSE) returns the free-intercept fit", {
    set.seed(4)
    m <- 60L
    res <- pecotmr:::.lderGetRes(
        x = rnorm(m),
        lam = abs(rnorm(m, 3, 0.5)),
        ldAnnot = matrix(1, nrow = m, ncol = 1L),
        N = 1000,
        M = 5000,
        M_a = 5000,
        rough = FALSE,
        twostage = FALSE
    )
    # One pass, no chi2 outlier filtering, intercept left free.
    expect_setequal(names(res), c("h2", "h2a", "tau", "a"))
    expect_true(is.finite(res$a))
})

test_that(".lderLocalBlock reports NA for a block too small to fit", {
    out <- pecotmr:::.lderLocalBlock(
        list(lam = c(1, 2), x = c(0.5, 0.6)),
        b = 7L,
        N = 1000,
        a = 1
    )
    # Fewer than 3 eigenvalues: a slope cannot be identified.
    expect_equal(out$blockId, 7L)
    expect_true(is.na(out$h2Local))
    expect_true(is.na(out$h2LocalSe))
})

test_that(".h2Enrichment names unnamed annotation columns positionally", {
    baselineMat <- matrix(c(1, 0, 1, 1), nrow = 2L)
    expect_null(colnames(baselineMat))
    out <- pecotmr:::.h2Enrichment(
        baselineMat,
        tau = c(0.1, 0.2),
        tauSe = c(0.01, 0.02),
        annotationJackknifeCoefs = matrix(c(0.1, 0.2, 0.1, 0.2), nrow = 2L),
        h2 = 0.3
    )
    expect_equal(out$annotation, c("annot_1", "annot_2"))
})

test_that(".sldscBlockIndex errors when no block covers a variant", {
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
        suppressWarnings(pecotmr:::.sldscBlockIndex(disjoint)),
        "no LD block covers any reference variant"
    )
    covering <- mk(GenomicRanges::GRanges("chr22", IRanges::IRanges(1, 10000)))
    idx <- pecotmr:::.sldscBlockIndex(covering)
    expect_equal(unname(lengths(idx)), n)
})

test_that(".sldscIrwls derives its own starting weights when none are given", {
    set.seed(3)
    nv <- 40L
    scores <- matrix(abs(rnorm(nv * 2L, 5, 1)), nrow = nv)
    y <- 1 + rowSums(scores) * 1e-4 * 500 + rnorm(nv, sd = 0.05)
    fit <- pecotmr:::.sldscIrwls(
        y = y,
        scores = scores,
        baseScore = rowSums(scores),
        n = 500,
        M = 1000,
        nIter = 2L,
        freeIntercept = FALSE,
        intercept = 1,
        w0 = NULL
    )
    expect_setequal(names(fit), c("coef", "intercept"))
    expect_true(all(is.finite(fit$coef)))
})

test_that(".sldscLocal clamps a negative LD score to zero", {
    out <- pecotmr:::.sldscLocal(
        chi2 = c(1.2, 3.4),
        baseScore = c(10, -2),
        fit = list(tau = c(1e-5))
    )
    expect_equal(out$variantIdx, 1:2)
    expect_equal(out$ldScore, c(10, -2))
    # A negative LD score cannot contribute negative heritability.
    expect_equal(out$h2Local, c(1e-4, 0))
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
    expect_equal(
        heritabilityEstimate(estimateH2(split, ref, method = "lder")),
        heritabilityEstimate(estimateH2(ss, ref, method = "lder"))
    )
    # Two studies still need one named.
    two <- GwasSumStats(
        studyName = c("a", "b"),
        entry = list(gr, gr),
        genome = "hg19",
        ldReference = makeTestGwasGenotypeHandle()
    )
    expect_error(
        estimateH2(two, ref, method = "lder"),
        "`studyName` is required when the GwasSumStats has 2 studies"
    )
})
