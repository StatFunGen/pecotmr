# Tests migrated from test_h2ClassesSumstats.R

# === Tests migrated from test_h2ClassesSumstats.R (LdEigen) ===

test_that("LdEigen constructs and validates correctly", {
    ldblocks <- make_test_ldblocks()
    snp_info <- make_test_snp_info()
    eigen_list <- list(
        list(
            values = c(1, 0.5),
            vectors = matrix(rnorm(20), 10, 2),
            snpIdx = 1:10
        ),
        list(values = c(0.8), vectors = matrix(rnorm(10), 10, 1), snpIdx = 1:10)
    )

    obj <- LdEigen(
        ldBlocks = ldblocks,
        snpInfo = snp_info,
        nRef = 500L,
        inSample = FALSE,
        genome = "hg19",
        eigenList = eigen_list,
        eigenvalueTruncation = 0.9
    )
    expect_s4_class(obj, "LdEigen")
    expect_true(methods::validObject(obj))
})


test_that("LdEigen rejects eigen_list length mismatch", {
    ldblocks <- make_test_ldblocks() # 2 blocks
    # Only 1 element in eigen_list
    expect_error(
        methods::validObject(
            LdEigen(
                ldBlocks = ldblocks,
                snpInfo = make_test_snp_info(),
                nRef = 500L,
                inSample = FALSE,
                genome = "hg19",
                eigenList = list(list(values = 1)),
                eigenvalueTruncation = 0.9
            )
        ),
        "eigenList.*must match"
    )
})


test_that("LdEigen rejects invalid eigenvalue_truncation", {
    ldblocks <- make_test_ldblocks()
    expect_error(
        methods::validObject(
            LdEigen(
                ldBlocks = ldblocks,
                snpInfo = make_test_snp_info(),
                nRef = 500L,
                inSample = FALSE,
                genome = "hg19",
                eigenList = list(list(), list()),
                eigenvalueTruncation = 0
            )
        ),
        "eigenvalueTruncation"
    )
})

# show() smoke test, moved here from test_showMethods.R so the test
# tree mirrors R/.
test_that("show(LdEigen) does not error", {
    eig <- LdEigen(
        ldBlocks = make_test_ldblocks(),
        snpInfo = make_test_snp_info(),
        nRef = 500L,
        inSample = FALSE,
        genome = "hg19",
        eigenList = list(list(), list()),
        eigenvalueTruncation = 0.9
    )
    expect_output(show(eig), "LdEigen")
})

test_that("subsetting an LdEigen is refused, not silently allowed", {
    # The guard exists because eigenList is per LD BLOCK while the ranges are
    # per variant: narrowing the ranges would leave decompositions describing
    # variants the object no longer carries.
    obj <- LdEigen(
        ldBlocks = make_test_ldblocks(),
        snpInfo = make_test_snp_info(),
        nRef = 500L,
        inSample = FALSE,
        genome = "hg19",
        eigenList = list(
            list(
                values = c(1, 0.5),
                vectors = matrix(rnorm(20), 10, 2),
                snpIdx = 1:10
            ),
            list(
                values = c(0.8),
                vectors = matrix(rnorm(10), 10, 1),
                snpIdx = 1:10
            )
        ),
        eigenvalueTruncation = 0.9
    )
    expect_error(obj[1:5], "cannot be subset")
    expect_error(obj[1], "Recompute over the subset")
})


# =============================================================================
# buildLdEigen
# =============================================================================

test_that("buildLdEigen decomposes a single-block LdData", {
    ld <- makeTestLdData(n = 6L)
    ref <- buildLdEigen(ld)

    expect_s4_class(ref, "LdEigen")
    expect_equal(length(ref), 6L)
    expect_equal(length(getEigenList(ref)), 1L)
    expect_equal(length(getLdBlocks(ref)), 1L)
    expect_equal(getNRef(ref), 500L)
    expect_false(getInSample(ref))
    expect_equal(getEigenList(ref)[[1]]$snpIdx, 1:6)
})

test_that("buildLdEigen gives one block per matrix of a multi-block LdData", {
    ref <- buildLdEigen(makeTestLdDataMultiBlock(sizes = c(4L, 3L)))

    expect_equal(length(ref), 7L)
    blocks <- getEigenList(ref)
    expect_equal(length(blocks), 2L)
    # snpIdx must index the concatenated variant order, not each block's own.
    expect_equal(blocks[[1]]$snpIdx, 1:4)
    expect_equal(blocks[[2]]$snpIdx, 5:7)
    expect_equal(length(getLdBlocks(ref)), 2L)
})

test_that("buildLdEigen reconstructs the correlation it was given", {
    ld <- makeTestLdData(n = 6L)
    block <- getEigenList(buildLdEigen(ld))[[1]]
    rebuilt <- block$vectors %*% diag(block$values) %*% t(block$vectors)
    expect_equal(rebuilt, unname(getCorrelation(ld)), tolerance = 1e-10)
})

test_that("buildLdEigen carries variant identity across from the LdData", {
    ref <- buildLdEigen(makeTestLdData(n = 6L))
    md <- S4Vectors::mcols(ref, use.names = FALSE)

    # An LdData names variants `variant_id` and reports `allele_freq`; an
    # LdStatistic wants SNP and MAF.
    expect_equal(names(ref)[[1]], "chr1:1000:C:T")
    expect_equal(as.character(md$A1)[[1]], "T")
    expect_equal(as.character(md$A2)[[1]], "C")
    expect_equal(
        md$MAF,
        pmin(seq(0.2, 0.8, length.out = 6), seq(0.8, 0.2, length.out = 6))
    )
})

test_that("buildLdEigen truncates to the requested eigenvalue mass", {
    ld <- makeTestLdData(n = 6L)
    full <- getEigenList(buildLdEigen(ld))[[1]]
    cut <- getEigenList(buildLdEigen(ld, eigenvalueTruncation = 0.9))[[1]]

    expect_equal(length(full$values), 6L)
    expect_lt(length(cut$values), 6L)
    expect_equal(ncol(cut$vectors), length(cut$values))
    # The retained components are the leading ones, unchanged.
    expect_equal(cut$values, full$values[seq_along(cut$values)])
})

test_that("buildLdEigen prefers an explicit nRef, inSample and genome", {
    ref <- buildLdEigen(
        makeTestLdData(),
        nRef = 12345L,
        inSample = TRUE,
        genome = "hg38"
    )
    expect_equal(getNRef(ref), 12345L)
    expect_true(getInSample(ref))
    expect_equal(getGenome(ref), "hg38")
})


test_that("eigenvalue truncation keeps everything when there is no mass", {
    # All-nonpositive eigenvalues carry no cumulative mass to threshold
    # against, so every index is kept rather than none.
    expect_equal(pecotmr:::.ldEigenKeep(c(0, -1, -2), 0.9), 1:3)
    # ...and a truncation of 1 keeps everything by definition.
    expect_equal(pecotmr:::.ldEigenKeep(c(3, 2, 1), 1), 1:3)
})
