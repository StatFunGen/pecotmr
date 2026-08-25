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
