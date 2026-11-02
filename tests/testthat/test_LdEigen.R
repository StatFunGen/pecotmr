# Tests migrated from test_h2ClassesSumstats.R

# === Tests migrated from test_h2ClassesSumstats.R (LdEigen) ===

test_that("LdEigen constructs and validates correctly", {
    ldblocks <- makeTestLdBlocks()
    snp_info <- makeTestSnpInfo()
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
        nSamples = 500L,
        inSample = FALSE,
        genome = "hg19",
        eigenList = eigen_list,
        eigenvalueTruncation = 0.9
    )
    expect_s4_class(obj, "LdEigen")
    expect_true(methods::validObject(obj))
})


test_that("LdEigen rejects eigen_list length mismatch", {
    ldblocks <- makeTestLdBlocks() # 2 blocks
    # Only 1 element in eigen_list
    expect_error(
        methods::validObject(
            LdEigen(
                ldBlocks = ldblocks,
                snpInfo = makeTestSnpInfo(),
                nSamples = 500L,
                inSample = FALSE,
                genome = "hg19",
                eigenList = list(list(values = 1)),
                eigenvalueTruncation = 0.9
            )
        ),
        "Variable 'eigenList'.*Must have length"
    )
})


test_that("LdEigen rejects invalid eigenvalue_truncation", {
    ldblocks <- makeTestLdBlocks()
    expect_error(
        methods::validObject(
            LdEigen(
                ldBlocks = ldblocks,
                snpInfo = makeTestSnpInfo(),
                nSamples = 500L,
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
        ldBlocks = makeTestLdBlocks(),
        snpInfo = makeTestSnpInfo(),
        nSamples = 500L,
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
        ldBlocks = makeTestLdBlocks(),
        snpInfo = makeTestSnpInfo(),
        nSamples = 500L,
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
    # head() and tail() reach the rows without `[`; they used to return the
    # first variants with every block's decomposition still attached.
    expect_error(head(obj, 3), "cannot be subset")
    expect_error(tail(obj, 3), "cannot be subset")
})


# =============================================================================
# buildLdEigen
# =============================================================================

test_that("buildLdEigen decomposes a single LdMatrix as one block", {
    ld <- makeTestLdMatrix(n = 6L)
    ref <- buildLdEigen(ld)

    expect_s4_class(ref, "LdEigen")
    expect_equal(length(ref), 6L)
    expect_equal(length(eigenList(ref)), 1L)
    expect_equal(length(ldBlocks(ref)), 1L)
    expect_equal(ldPanelNSamples(ref), 500L)
    expect_false(inSample(ref))
    expect_equal(eigenList(ref)[[1]]$snpIdx, 1:6)
})

test_that("buildLdEigen gives one block per LdMatrix in a list", {
    ref <- buildLdEigen(makeTestLdMatrixBlocks(sizes = c(4L, 3L)))

    expect_equal(length(ref), 7L)
    blocks <- eigenList(ref)
    expect_equal(length(blocks), 2L)
    # snpIdx must index the concatenated variant order, not each block's own.
    expect_equal(blocks[[1]]$snpIdx, 1:4)
    expect_equal(blocks[[2]]$snpIdx, 5:7)
    expect_equal(length(ldBlocks(ref)), 2L)
})

test_that("buildLdEigen reconstructs the correlation it was given", {
    ld <- makeTestLdMatrix(n = 6L)
    block <- eigenList(buildLdEigen(ld))[[1]]
    rebuilt <- block$vectors %*% diag(block$values) %*% t(block$vectors)
    expect_equal(rebuilt, unname(ldMatrix(ld)), tolerance = 1e-10)
})

test_that("buildLdEigen carries variant identity across from an LdMatrix", {
    ref <- buildLdEigen(makeTestLdMatrix(n = 6L))
    md <- S4Vectors::mcols(ref, use.names = FALSE)

    # An LdMatrix records the frequency of A1 as `AF`; an LdStatistic wants
    # it folded, as MAF.
    expect_equal(names(ref)[[1]], "chr1:1000:C:T")
    expect_equal(as.character(md$A1)[[1]], "T")
    expect_equal(as.character(md$A2)[[1]], "C")
    expect_equal(
        md$MAF,
        pmin(seq(0.2, 0.8, length.out = 6), seq(0.8, 0.2, length.out = 6))
    )
})

test_that("buildLdEigen truncates to the requested eigenvalue mass", {
    ld <- makeTestLdMatrix(n = 6L)
    full <- eigenList(buildLdEigen(ld))[[1]]
    cut <- eigenList(buildLdEigen(ld, eigenvalueTruncation = 0.9))[[1]]

    expect_equal(length(full$values), 6L)
    expect_lt(length(cut$values), 6L)
    expect_equal(ncol(cut$vectors), length(cut$values))
    # The retained components are the leading ones, unchanged.
    expect_equal(cut$values, full$values[seq_along(cut$values)])
})

test_that("buildLdEigen prefers an explicit nRef, inSample and genome", {
    ref <- buildLdEigen(
        makeTestLdMatrix(),
        nSamples = 12345L,
        inSample = TRUE,
        genome = "hg38"
    )
    expect_equal(ldPanelNSamples(ref), 12345L)
    expect_true(inSample(ref))
    expect_equal(unname(GenomeInfoDb::genome(ref)), "hg38")
})


test_that("eigenvalue truncation keeps everything when there is no mass", {
    # All-nonpositive eigenvalues carry no cumulative mass to threshold
    # against, so every index is kept rather than none.
    expect_equal(pecotmr:::.ldEigenKeep(c(0, -1, -2), 0.9), 1:3)
    # ...and a truncation of 1 keeps everything by definition.
    expect_equal(pecotmr:::.ldEigenKeep(c(3, 2, 1), 1), 1:3)
})

test_that("LdEigen: argument guards fire", {
    ldblocks <- makeTestLdBlocks()
    snp_info <- makeTestSnpInfo(10)
    eigen_list <- list(
        list(values = c(1, 0.5), vectors = matrix(0, 10, 2), snpIdx = 1:10),
        list(values = 0.8, vectors = matrix(0, 10, 1), snpIdx = 1:10)
    )
    ok <- list(
        ldBlocks = ldblocks,
        snpInfo = snp_info,
        nSamples = 500L,
        eigenList = eigen_list
    )
    badSnpInfo <- list_modify(ok, snpInfo = "not-a-df")
    expect_error(
        exec(LdEigen, !!!badSnpInfo),
        "snpInfo.*Must be of type 'data.frame'"
    )
    badEigenList <- list_modify(ok, eigenList = "not-a-list")
    expect_error(
        exec(LdEigen, !!!badEigenList),
        "eigenList.*Must be of type 'list'"
    )
    zeroSamples <- list_modify(ok, nSamples = 0L)
    expect_error(
        exec(LdEigen, !!!zeroSamples),
        "nSamples.*Must be >= 1"
    )
    textSamples <- list_modify(ok, nSamples = "500")
    expect_error(
        exec(LdEigen, !!!textSamples),
        "nSamples.*Must be of type 'count'"
    )
})

test_that("buildLdEigen: argument guards fire", {
    ld <- makeTestLdMatrix()
    expect_error(buildLdEigen(ld, nSamples = 0L), "nSamples.*Must be >= 1")
    expect_error(
        buildLdEigen(ld, inSample = NA),
        "inSample.*May not be NA"
    )
})

test_that("buildLdEigen from a panel decomposes each block when asked", {
    fx <- .h2PanelFixture()
    lazy <- buildLdEigen(fx$panel, ldBlocks = fx$blocks, genome = "hg19")
    stored <- buildLdEigen(fx$ld, genome = "hg19")
    expect_length(methods::slot(lazy, "eigenList"), 0L)
    expect_identical(ldReference(lazy), fx$panel)
    blocks <- eigenList(lazy)
    expect_length(blocks, 3L)
    # Blocks of 60, 140 and 100 variants from 100 samples: the middle one
    # lists 100 directions and counts its other 40 zero-eigenvalue ones.
    expect_equal(map_int(blocks, "nNull"), c(0L, 40L, 0L))
    expect_equal(lengths(map(blocks, "values")), c(60L, 100L, 100L))
    expect_equal(
        blocks[[2L]]$values[1:90],
        eigenList(stored)[[2L]]$values[1:90]
    )
    expect_equal(map(blocks, "snpIdx"), map(eigenList(stored), "snpIdx"))
    out <- capture.output(show(lazy))
    expect_true(any(grepl("across 3 blocks", out)))
    expect_true(any(grepl("computed from the genotype panel", out)))
    expect_error(lazy[1:2], "cannot be subset")
})

test_that("an LdEigen stores one decomposition per block or reads them all", {
    fx <- .h2PanelFixture()
    lazy <- buildLdEigen(fx$panel, ldBlocks = fx$blocks, genome = "hg19")
    stored <- buildLdEigen(fx$ld, genome = "hg19")
    mixed <- lazy
    mixed@eigenList <- eigenList(stored)
    expect_error(validObject(mixed), "Variable 'eigenList'.*Must have length")
})
