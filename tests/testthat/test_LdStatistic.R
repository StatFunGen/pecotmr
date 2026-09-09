# Tests for R/LdStatistic.R (virtual base class)
# getGenome() is defined on the virtual LdStatistic and inherited by its
# concrete subclasses (LdEigen / LdScore); exercise it through a concrete
# LdScore instance. Fixtures (make_test_ldblocks / make_test_snp_info) come
# from helper-h2Classes.R.

test_that("getGenome returns the genome build string (via an LdScore subclass)", {
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
    expect_equal(getGenome(obj), "hg19")
})

# ===========================================================================
# Validity and the shared ranges builder
#
# Every test above builds a VALID object, so the validity function's error
# branches and the ranges builder's missing-column guard were never executed.
# ===========================================================================

# Overrides REPLACE the defaults rather than being appended: forwarding `...`
# alongside a fixed `nRef = 500L` makes .lds_score(nRef = 0L) fail with
# "formal argument matched by multiple actual arguments", which an
# expect_error(., "nRef") happily matches without ever reaching validity.
.lds_score <- function(n = 10, ...) {
    args <- list(
        ldBlocks = make_test_ldblocks(),
        snpInfo = make_test_snp_info(n),
        nRef = 500L,
        inSample = FALSE,
        genome = "hg19",
        ldScores = matrix(runif(n), nrow = n, ncol = 1),
        ldScoreWeights = runif(n),
        ldMatrixList = list()
    )
    over <- list(...)
    args[names(over)] <- over
    exec(LdScore, !!!args)
}

test_that("validity rejects an nRef that is not a single positive integer", {
    # Matched on the validity text, not just "nRef": the argument name alone
    # also appears in R's own argument-matching errors.
    msg <- "'nRef' must be a single positive integer"
    expect_error(.lds_score(nRef = 0L), msg)
    expect_error(.lds_score(nRef = -1L), msg)
    expect_error(.lds_score(nRef = c(10L, 20L)), msg)
})

test_that("validity rejects an inSample that is not a single flag", {
    # Through new(), not LdScore(): the constructor writes
    # `inSample = isTRUE(inSample)`, so a multi-element or empty flag is
    # collapsed before validity can see it. (`nRef = as.integer(nRef)` does
    # not collapse length, which is why the nRef branch IS reachable above.)
    # This one guards direct new() use -- migration scripts, deserialization.
    gr <- as(.lds_score(), "GRanges")
    msg <- "'inSample' must be a single logical value"
    mk <- function(flag) {
        methods::new(
            "LdScore",
            gr,
            ldBlocks = make_test_ldblocks(),
            nRef = 500L,
            inSample = flag,
            ldMatrixList = list()
        )
    }
    expect_error(mk(c(TRUE, FALSE)), msg)
    expect_error(mk(logical(0)), msg)
})

test_that("validity rejects a statistic carrying no variants", {
    # Reached through new() rather than LdScore(): the constructor cannot
    # produce an empty statistic at all, because .ldStatRanges() builds
    # IRanges(start = integer(0), width = 1L) and those lengths do not
    # recycle. So this branch guards direct new() use, which is what the
    # migration scripts and deserialization do.
    gr <- GenomicRanges::GRanges()
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = character(0),
        A1 = character(0),
        A2 = character(0),
        ldScores = I(matrix(numeric(0), nrow = 0, ncol = 1)),
        ldScoreWeights = numeric(0)
    )
    expect_error(
        methods::new(
            "LdScore",
            gr,
            ldBlocks = make_test_ldblocks(),
            nRef = 500L,
            inSample = FALSE,
            ldMatrixList = list()
        ),
        "at least one variant"
    )
})

test_that(".ldStatRanges names the snpInfo columns it is missing", {
    si <- make_test_snp_info(4)
    expect_error(
        LdScore(
            ldBlocks = make_test_ldblocks(),
            snpInfo = si[, setdiff(colnames(si), c("A1", "A2"))],
            nRef = 500L,
            inSample = FALSE,
            genome = "hg19",
            ldScores = matrix(runif(4), nrow = 4, ncol = 1),
            ldScoreWeights = runif(4),
            ldMatrixList = list()
        ),
        "missing column\\(s\\): A1, A2"
    )
})

test_that("getLdBlocks returns the blocks the statistic was built against", {
    obj <- .lds_score()
    expect_s4_class(getLdBlocks(obj), "GRanges")
    expect_equal(length(getLdBlocks(obj)), 2L)
})


# =============================================================================
# .ldRefPrepare -- the input handling shared by buildLdEigen/buildLdScore
# =============================================================================

test_that("the builders accept a bare LdData as well as a list of them", {
    one <- makeTestLdData(n = 4L, startBp = 1000L)
    two <- makeTestLdData(n = 3L, startBp = 5000L)

    expect_equal(length(buildLdEigen(one)), 4L)
    expect_equal(length(buildLdEigen(list(one))), 4L)

    # A list of LdData concatenates: one LD block per element here, and the
    # variant order follows the list order.
    joined <- buildLdEigen(list(one, two))
    expect_equal(length(joined), 7L)
    expect_equal(length(getEigenList(joined)), 2L)
    expect_equal(getEigenList(joined)[[2]]$snpIdx, 5:7)
})

test_that("the builders reject input that is not LdData", {
    expect_error(buildLdEigen("not-ld"), "must be an LdData")
    expect_error(buildLdEigen(list()), "must be an LdData")
    expect_error(
        buildLdScore(list(makeTestLdData(), "not-ld")),
        "element\\(s\\) 2 are not LdData"
    )
})

test_that("nRef is taken from the LdData, and disagreement is an error", {
    a <- makeTestLdData(n = 4L, nRef = 500L)
    b <- makeTestLdData(n = 3L, startBp = 5000L, nRef = 900L)

    expect_equal(getNRef(buildLdEigen(list(a, a))), 500L)
    expect_error(buildLdEigen(list(a, b)), "differing reference panel sizes")
    # An explicit nRef settles it.
    expect_equal(getNRef(buildLdEigen(list(a, b), nRef = 700L)), 700L)
})

test_that("the genome build comes from the LdData when the caller names none", {
    ld <- makeTestLdData()
    # loadLdMatrix() leaves the build unset, which must stay NA rather than
    # becoming a made-up default.
    expect_true(is.na(getGenome(buildLdEigen(ld))))

    tagged <- ld
    GenomeInfoDb::genome(tagged) <- "hg38"
    expect_equal(getGenome(buildLdEigen(tagged)), "hg38")
    # An explicit argument still wins.
    expect_equal(getGenome(buildLdEigen(tagged, genome = "hg19")), "hg19")
})

test_that("LD block ranges span the variants each block covers", {
    ref <- buildLdEigen(makeTestLdDataMultiBlock(sizes = c(4L, 3L)))
    blocks <- getLdBlocks(ref)
    variants <- GenomicRanges::start(ref)

    expect_equal(
        as.character(GenomicRanges::seqnames(blocks)),
        c("chr1", "chr1")
    )
    expect_equal(GenomicRanges::start(blocks), c(variants[[1]], variants[[5]]))
    expect_equal(GenomicRanges::end(blocks), c(variants[[4]], variants[[7]]))
})

test_that("an LD block spanning chromosomes is rejected", {
    ld <- makeTestLdData(n = 4L)
    gr <- as(ld, "GRanges")
    GenomeInfoDb::seqlevels(gr) <- c("chr1", "chr2")
    GenomicRanges::seqnames(gr) <- factor(
        c("chr1", "chr1", "chr2", "chr2"),
        levels = c("chr1", "chr2")
    )
    crossChrom <- LdData(
        correlation = getCorrelation(ld),
        variants = gr,
        blockMetadata = getBlockMetadata(ld),
        nRef = getNRef(ld)
    )
    expect_error(buildLdEigen(crossChrom), "spans 2 chromosomes")
})
