# Tests migrated from test_h2ClassesSumstats.R

# === Tests migrated from test_h2ClassesSumstats.R (LdScore) ===

test_that("LdScore constructs and validates correctly", {
    ldblocks <- make_test_ldblocks()
    n <- 10
    snp_info <- make_test_snp_info(n)

    obj <- LdScore(
        ldBlocks = ldblocks,
        snpInfo = snp_info,
        nRef = 500L,
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
    ldblocks <- make_test_ldblocks()
    snp_info <- make_test_snp_info(10)

    expect_error(
        methods::validObject(
            LdScore(
                ldBlocks = ldblocks,
                snpInfo = snp_info,
                nRef = 500L,
                inSample = FALSE,
                genome = "hg19",
                ldScores = matrix(0, nrow = 5, ncol = 1), # wrong rows
                ldScoreWeights = runif(10),
                ldMatrixList = list()
            )
        ),
        "ldScores.*must be parallel"
    )
})

# show() smoke test, moved here from test_showMethods.R so the test
# tree mirrors R/.
test_that("show(LdScore) does not error", {
    n <- 10
    lsr <- LdScore(
        ldBlocks = make_test_ldblocks(),
        snpInfo = make_test_snp_info(n),
        nRef = 500L,
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
            ldBlocks = make_test_ldblocks(),
            snpInfo = make_test_snp_info(n),
            nRef = 500L,
            inSample = FALSE,
            genome = "hg19",
            ldScores = matrix(runif(n), nrow = n, ncol = 1),
            ldScoreWeights = runif(n - 1L),
            ldMatrixList = list()
        ),
        "they must be parallel"
    )
})

test_that("getLdScoreWeights returns the per-variant weights", {
    n <- 10
    w <- runif(n)
    obj <- LdScore(
        ldBlocks = make_test_ldblocks(),
        snpInfo = make_test_snp_info(n),
        nRef = 500L,
        inSample = FALSE,
        genome = "hg19",
        ldScores = matrix(runif(n), nrow = n, ncol = 1),
        ldScoreWeights = w,
        ldMatrixList = list()
    )
    expect_equal(getLdScoreWeights(obj), w)
})

test_that("validity requires the score columns to be present in mcols", {
    # The columns live in mcols now, so dropping one is an ordinary mcols
    # edit rather than a slot edit -- which is exactly why validity checks it.
    n <- 10
    obj <- LdScore(
        ldBlocks = make_test_ldblocks(),
        snpInfo = make_test_snp_info(n),
        nRef = 500L,
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

test_that("buildLdScore computes per-block sums of r^2", {
    ld <- makeTestLdData(n = 6L)
    ref <- buildLdScore(ld)

    expect_s4_class(ref, "LdScore")
    expect_equal(length(ref), 6L)
    expect_equal(colnames(getLdScores(ref)), "base_l2")
    expect_equal(
        as.vector(getLdScores(ref)[, 1]),
        rowSums(unname(getCorrelation(ld))^2)
    )
})

test_that("buildLdScore scores each block against only its own variants", {
    ld <- makeTestLdDataMultiBlock(sizes = c(4L, 3L))
    scores <- as.vector(getLdScores(buildLdScore(ld))[, 1])
    perBlock <- unlist(lapply(getCorrelation(ld), function(R) rowSums(R^2)))

    expect_equal(length(scores), 7L)
    expect_equal(scores, perBlock)
})

# buildLdScore and computeLdScores(LdEigen) are two routes to the same
# quantity -- sum_k r^2_jk versus sum_i V[j,i]^2 d[i]^2 -- so a disagreement
# means one of them has drifted.
test_that("buildLdScore agrees with computeLdScores on the same reference", {
    for (ld in list(makeTestLdData(n = 6L), makeTestLdDataMultiBlock())) {
        expect_equal(
            as.vector(getLdScores(buildLdScore(ld))[, 1]),
            as.vector(computeLdScores(buildLdEigen(ld))[, 1])
        )
    }
})

test_that("buildLdScore keeps per-block LD matrices for g-LDSC by default", {
    ld <- makeTestLdDataMultiBlock(sizes = c(4L, 3L))
    mats <- getLdMatrixList(buildLdScore(ld))

    expect_equal(length(mats), 2L)
    expect_equal(dim(mats[[1]]$R), c(4L, 4L))
    expect_equal(mats[[1]]$snpIdx, 1:4)
    expect_equal(mats[[2]]$snpIdx, 5:7)
    expect_equal(
        length(getLdMatrixList(buildLdScore(ld, keepLdMatrices = FALSE))),
        0L
    )
})

test_that("buildLdScore defaults weights to 1/max(l2, 1)", {
    ld <- makeTestLdData(n = 6L)
    ref <- buildLdScore(ld)
    l2 <- as.vector(getLdScores(ref)[, 1])
    expect_equal(getLdScoreWeights(ref), 1 / pmax(l2, 1))

    custom <- buildLdScore(ld, ldScoreWeights = rep(2, 6))
    expect_equal(getLdScoreWeights(custom), rep(2, 6))
    expect_error(
        buildLdScore(ld, ldScoreWeights = rep(2, 3)),
        "3 value\\(s\\) for 6 variant\\(s\\)"
    )
})
