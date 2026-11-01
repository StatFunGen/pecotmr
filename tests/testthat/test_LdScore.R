# Tests migrated from test_h2ClassesSumstats.R

# === Tests migrated from test_h2ClassesSumstats.R (LdScore) ===

test_that("LdScore constructs and validates correctly", {
    ldblocks <- makeTestLdBlocks()
    n <- 10
    snp_info <- makeTestSnpInfo(n)

    obj <- LdScore(
        ldBlocks = ldblocks,
        snpInfo = snp_info,
        nSamples = 500L,
        inSample = FALSE,
        genome = "hg19",
        ldScores = matrix(runif(n), nrow = n, ncol = 1),
        ldScoreWeights = runif(n),
        ldMatrixList = list()
    )
    expect_s4_class(obj, "LdScore")
    expect_true(methods::validObject(obj))
})


test_that("LdScore rejects ld_scores row mismatch with snp_info", {
    ldblocks <- makeTestLdBlocks()
    snp_info <- makeTestSnpInfo(10)

    expect_error(
        methods::validObject(
            LdScore(
                ldBlocks = ldblocks,
                snpInfo = snp_info,
                nSamples = 500L,
                inSample = FALSE,
                genome = "hg19",
                ldScores = matrix(0, nrow = 5, ncol = 1), # wrong rows
                ldScoreWeights = runif(10),
                ldMatrixList = list()
            )
        ),
        "ldScores.*Must have exactly 10 rows, but has 5 rows"
    )
})

# show() smoke test, moved here from test_showMethods.R so the test
# tree mirrors R/.
test_that("show(LdScore) does not error", {
    n <- 10
    lsr <- LdScore(
        ldBlocks = makeTestLdBlocks(),
        snpInfo = makeTestSnpInfo(n),
        nSamples = 500L,
        inSample = FALSE,
        genome = "hg19",
        ldScores = matrix(1, nrow = n, ncol = 1),
        ldScoreWeights = rep(1, n),
        ldMatrixList = list()
    )
    expect_output(show(lsr), "LdScore")
})

test_that("LdScore rejects weights that are not parallel to the variants", {
    n <- 10
    expect_error(
        LdScore(
            ldBlocks = makeTestLdBlocks(),
            snpInfo = makeTestSnpInfo(n),
            nSamples = 500L,
            inSample = FALSE,
            genome = "hg19",
            ldScores = matrix(runif(n), nrow = n, ncol = 1),
            ldScoreWeights = runif(n - 1L),
            ldMatrixList = list()
        ),
        "ldScoreWeights.*Must have length 10, but has length 9"
    )
})

test_that("mcols exposes the per-variant LD-score weights", {
    n <- 10
    w <- runif(n)
    obj <- LdScore(
        ldBlocks = makeTestLdBlocks(),
        snpInfo = makeTestSnpInfo(n),
        nSamples = 500L,
        inSample = FALSE,
        genome = "hg19",
        ldScores = matrix(runif(n), nrow = n, ncol = 1),
        ldScoreWeights = w,
        ldMatrixList = list()
    )
    expect_equal(mcols(obj)$ldScoreWeights, w)
})

test_that("validity requires the score columns to be present in mcols", {
    # The columns live in mcols now, so dropping one is an ordinary mcols
    # edit rather than a slot edit -- which is exactly why validity checks it.
    n <- 10
    obj <- LdScore(
        ldBlocks = makeTestLdBlocks(),
        snpInfo = makeTestSnpInfo(n),
        nSamples = 500L,
        inSample = FALSE,
        genome = "hg19",
        ldScores = matrix(runif(n), nrow = n, ncol = 1),
        ldScoreWeights = runif(n),
        ldMatrixList = list()
    )
    bad <- obj
    S4Vectors::mcols(bad)$ldScoreWeights <- NULL
    expect_error(methods::validObject(bad), "ldScoreWeights")
})


# =============================================================================
# buildLdScore
# =============================================================================

test_that("buildLdScore sums ldsc's unbiased r^2 within each block", {
    ld <- makeTestLdMatrix(n = 6L)
    ref <- buildLdScore(ld)
    r2 <- unname(ldMatrix(ld))^2

    expect_s4_class(ref, "LdScore")
    expect_equal(length(ref), 6L)
    expect_equal(colnames(mcols(ref)$ldScores), "base_l2")
    # makeTestLdMatrix records a 500-sample panel.
    expect_equal(
        as.vector(mcols(ref)$ldScores[, 1]),
        rowSums(r2 - (1 - r2) / (500 - 2))
    )
})

test_that("buildLdScore scores each block against only its own variants", {
    ld <- makeTestLdMatrixBlocks(sizes = c(4L, 3L))
    scores <- as.vector(mcols(buildLdScore(ld))$ldScores[, 1])
    perBlock <- unlist(map(
        map(ld, ldMatrix),
        function(R) rowSums(R^2 - (1 - R^2) / (500 - 2))
    ))

    expect_equal(length(scores), 7L)
    expect_equal(scores, unname(perBlock))
})

# The unbiased r^2 is ldsc's (ldscore.py, __l2_unbiased__). These inputs are
# exact in binary: r = 0.5 gives 0.25 - 0.75 / 8 = 0.15625, and r = 0 gives
# -1 / 8, which ldsc keeps rather than flooring. A variant's own term stays 1.
test_that("the unbiased r^2 is r^2 - (1 - r^2) / (n - 2), as in ldsc", {
    R <- matrix(c(1, 0.5, 0, 0.5, 1, -0.5, 0, -0.5, 1), 3L, 3L)
    ld <- LdMatrix(
        R,
        variants = .testLdVariants(3L, "chr1", 1000L),
        nSamples = 10L
    )
    # The raw sums would be c(1.25, 1.5, 1.25).
    expect_identical(
        as.vector(mcols(buildLdScore(ld))$ldScores[, 1]),
        c(1.03125, 1.3125, 1.03125)
    )
})

# Summed over a block with weights a_k, the unbiased r^2 is affine in the raw
# weighted sum: ((n - 1) * sum_k a_k r2_jk - sum_k a_k) / (n - 2). The LdEigen
# route relies on that to correct scores it never forms r^2 for, so it must
# hold exactly, for binary and continuous weights alike.
test_that("the summed form of the unbiased r^2 is exact under any weights", {
    set.seed(11)
    R <- cov2cor(crossprod(matrix(rnorm(40 * 8), 40, 8)))
    a <- cbind(base = 1, binary = c(1, 0, 0, 1, 1, 0, 1, 0), cont = runif(8))
    expect_equal(
        pecotmr:::.ldscUnbiasedSum((R^2) %*% a, colSums(a), 40L),
        pecotmr:::.ldscUnbiasedR2(R^2, 40L) %*% a,
        tolerance = 1e-12
    )
})

test_that("the unbiased r^2 needs a reference panel of more than 2", {
    # An LdMatrix that records no panel size leaves nothing to divide by.
    expect_error(
        buildLdScore(makeTestLdMatrix(nSamples = NA_integer_)),
        "`nSamples` is required"
    )
    expect_error(
        buildLdScore(makeTestLdMatrix(), nSamples = 2L),
        "needs more than 2 samples, and it records 2"
    )
    # Naming the panel size rescues an LdMatrix that records none.
    expect_s4_class(
        buildLdScore(makeTestLdMatrix(nSamples = NA_integer_), nSamples = 10L),
        "LdScore"
    )
})

# buildLdScore and computeLdScores(LdEigen) are two routes to the same
# quantity -- the first sums the unbiased r^2 over R, the second corrects
# sum_i V[j,i]^2 d[i]^2 (= sum_k r^2_jk) through the summed form -- so a
# disagreement means one of them has drifted.
test_that("buildLdScore agrees with computeLdScores on the same reference", {
    for (ld in list(makeTestLdMatrix(n = 6L), makeTestLdMatrixBlocks())) {
        expect_equal(
            as.vector(mcols(buildLdScore(ld))$ldScores[, 1]),
            as.vector(computeLdScores(buildLdEigen(ld))[, 1]),
            tolerance = 1e-12
        )
    }
})

test_that("buildLdScore keeps per-block LD matrices for g-LDSC by default", {
    ld <- makeTestLdMatrixBlocks(sizes = c(4L, 3L))
    mats <- ldMatrixList(buildLdScore(ld))

    expect_equal(length(mats), 2L)
    expect_equal(dim(mats[[1]]$R), c(4L, 4L))
    expect_equal(mats[[1]]$snpIdx, 1:4)
    expect_equal(mats[[2]]$snpIdx, 5:7)
    expect_equal(
        length(ldMatrixList(buildLdScore(ld, keepLdMatrices = FALSE))),
        0L
    )
})

test_that("an LdScore stores block matrices packed and hands them out whole", {
    ld <- makeTestLdMatrixBlocks(sizes = c(4L, 3L))
    scores <- buildLdScore(ld)
    stored <- methods::slot(scores, "ldMatrixList")
    expect_s4_class(stored[[1L]]$R, "dspMatrix")
    mats <- ldMatrixList(scores)
    expect_true(is.matrix(mats[[1L]]$R))
    expect_identical(mats[[1L]]$R, unname(ldMatrix(ld[[1L]])))
})

test_that("updateObject packs the block matrices of a saved LdScore", {
    scores <- buildLdScore(makeTestLdMatrixBlocks(sizes = c(4L, 3L)))
    saved <- scores
    methods::slot(saved, "ldMatrixList") <- ldMatrixList(scores)
    expect_true(is.matrix(methods::slot(saved, "ldMatrixList")[[1L]]$R))
    expect_identical(updateObject(saved), scores)
    expect_identical(updateObject(scores), scores)
})

test_that("buildLdScore defaults weights to 1/max(l2, 1)", {
    ld <- makeTestLdMatrix(n = 6L)
    ref <- buildLdScore(ld)
    l2 <- as.vector(mcols(ref)$ldScores[, 1])
    expect_equal(mcols(ref)$ldScoreWeights, 1 / pmax(l2, 1))

    custom <- buildLdScore(ld, ldScoreWeights = rep(2, 6))
    expect_equal(mcols(custom)$ldScoreWeights, rep(2, 6))
    expect_error(
        buildLdScore(ld, ldScoreWeights = rep(2, 3)),
        "ldScoreWeights.*Must have length 6, but has length 3"
    )
})

test_that("buildLdScore: argument guards fire", {
    ld <- makeTestLdMatrix()
    expect_error(buildLdScore(ld, nSamples = 0L), "nSamples.*Must be >= 1")
    expect_error(buildLdScore(ld, inSample = NA), "inSample.*May not be NA")
    expect_error(
        buildLdScore(ld, keepLdMatrices = NA),
        "keepLdMatrices.*May not be NA"
    )
})

test_that("buildLdScore reads a panel's blocks as LdMatrix blocks give them", {
    fx <- .h2PanelFixture()
    lazy <- buildLdScore(fx$panel, ldBlocks = fx$blocks, genome = "hg19")
    stored <- buildLdScore(fx$ld, genome = "hg19")
    expect_equal(
        S4Vectors::mcols(lazy)$ldScores,
        S4Vectors::mcols(stored)$ldScores,
        tolerance = 1e-10
    )
    expect_equal(
        S4Vectors::mcols(lazy)$ldScoreWeights,
        S4Vectors::mcols(stored)$ldScoreWeights,
        tolerance = 1e-10
    )
    # The panel is never copied into matrices, whatever keepLdMatrices says.
    expect_length(methods::slot(lazy, "ldMatrixList"), 0L)
    expect_length(methods::slot(stored, "ldMatrixList"), 3L)
    expect_null(ldReference(stored))
    out <- capture.output(show(lazy))
    expect_true(any(grepl("read from the genotype panel", out)))
    out <- capture.output(show(stored))
    expect_true(any(grepl("stored matrices", out)))
    bare <- buildLdScore(fx$ld, genome = "hg19", keepLdMatrices = FALSE)
    out <- capture.output(show(bare))
    expect_true(any(grepl("scores only", out)))
})

test_that("an LdScore storing block matrices refuses every subset", {
    # Each stored block indexes its variants by position, so a subset or a
    # reorder would leave the matrices describing other variants. head() and
    # rev() reach the rows without `[`, so they are refused as well.
    fx <- .h2PanelFixture()
    stored <- buildLdScore(fx$ld, genome = "hg19")
    msg <- "an LdScore storing per-block LD matrices cannot be subset"
    expect_error(stored[1:3], msg)
    expect_error(head(stored, 3), msg)
    expect_error(rev(stored), msg)
})

test_that("an LdScore without stored matrices subsets freely", {
    # Each block's variants are found afresh from the ranges, so the scores
    # and the block LD of the variants kept stay in step.
    fx <- .h2PanelFixture()
    lazy <- buildLdScore(fx$panel, ldBlocks = fx$blocks, genome = "hg19")
    bare <- buildLdScore(fx$ld, genome = "hg19", keepLdMatrices = FALSE)
    keep <- c(5:40, 100:150)
    for (obj in list(lazy, bare)) {
        sub <- obj[keep]
        expect_length(sub, length(keep))
        expect_equal(
            S4Vectors::mcols(sub)$ldScores,
            S4Vectors::mcols(obj)$ldScores[keep, , drop = FALSE]
        )
    }
    # The panel's LD over the kept variants is the full blocks' restricted to
    # them; the third block keeps no variant and so has no entry.
    full <- ldMatrixList(lazy)
    sub <- ldMatrixList(lazy[keep])
    expect_length(sub, 2L)
    expect_equal(sub[[1L]]$R, full[[1L]]$R[5:40, 5:40], tolerance = 1e-12)
    expect_equal(sub[[2L]]$R, full[[2L]]$R[40:90, 40:90], tolerance = 1e-12)
    expect_equal(map(sub, "snpIdx"), list(1:36, 37:87))
})
