# Tests for R/LdStatistic.R (virtual base class)
# The build is read from seqinfo() via GenomeInfoDb::genome(), which the
# virtual LdStatistic and its concrete subclasses (LdEigen / LdScore) all
# inherit from GRanges; exercise it through a concrete
# LdScore instance. Fixtures (makeTestLdBlocks / makeTestSnpInfo) come
# from helper-h2Classes.R.

test_that("genome() reports the build from seqinfo (via an LdScore)", {
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
    expect_equal(unname(GenomeInfoDb::genome(obj)), "hg19")
})

# ===========================================================================
# Validity and the shared ranges builder
#
# Every test above builds a VALID object, so the validity function's error
# branches and the ranges builder's missing-column guard were never executed.
# ===========================================================================

# Overrides REPLACE the defaults rather than being appended: forwarding `...`
# alongside a fixed `nSamples = 500L` makes .lds_score(nSamples = 0L) fail
# with "formal argument matched by multiple actual arguments", which an
# expect_error(., "nSamples") happily matches without ever reaching validity.
.lds_score <- function(n = 10, ...) {
    args <- list(
        ldBlocks = makeTestLdBlocks(),
        snpInfo = makeTestSnpInfo(n),
        nSamples = 500L,
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

test_that("validity rejects an nSamples that is not a positive integer", {
    # Matched on the validity text, not just the name: the argument name
    # alone
    # also appears in R's own argument-matching errors.
    msg <- "Variable 'nSamples'"
    expect_error(.lds_score(nSamples = 0L), msg)
    expect_error(.lds_score(nSamples = -1L), msg)
    expect_error(.lds_score(nSamples = c(10L, 20L)), msg)
})

test_that("validity rejects an inSample that is not a single flag", {
    # Through new(), not LdScore(): the constructor writes
    # `inSample = isTRUE(inSample)`, so a multi-element or empty flag is
    # collapsed before validity can see it. (`nSamples = as.integer(nSamples)`
    # does not collapse length, which is why the nSamples branch IS reachable
    # above.)
    # This one guards direct new() use -- migration scripts, deserialization.
    gr <- as(.lds_score(), "GRanges")
    msg <- "Variable 'inSample': Must have length 1"
    mk <- function(flag) {
        methods::new(
            "LdScore",
            gr,
            ldBlocks = makeTestLdBlocks(),
            nSamples = 500L,
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
            ldBlocks = makeTestLdBlocks(),
            nSamples = 500L,
            inSample = FALSE,
            ldMatrixList = list()
        ),
        "at least one variant"
    )
})

test_that(".ldStatRanges names the snpInfo columns it is missing", {
    si <- makeTestSnpInfo(4)
    expect_error(
        LdScore(
            ldBlocks = makeTestLdBlocks(),
            snpInfo = si[, setdiff(colnames(si), c("A1", "A2"))],
            nSamples = 500L,
            inSample = FALSE,
            genome = "hg19",
            ldScores = matrix(runif(4), nrow = 4, ncol = 1),
            ldScoreWeights = runif(4),
            ldMatrixList = list()
        ),
        "missing column\\(s\\): A1, A2"
    )
})

test_that("ldBlocks returns the blocks the statistic was built against", {
    obj <- .lds_score()
    expect_s4_class(ldBlocks(obj), "GRanges")
    expect_equal(length(ldBlocks(obj)), 2L)
})


# =============================================================================
# .h2RefPrepare -- the input handling shared by buildLdEigen/buildLdScore
# =============================================================================

test_that("the builders accept one LdMatrix as well as a list of them", {
    one <- makeTestLdMatrix(n = 4L, startBp = 1000L)
    two <- makeTestLdMatrix(n = 3L, startBp = 5000L)

    expect_equal(length(buildLdEigen(one)), 4L)
    expect_equal(length(buildLdEigen(list(one))), 4L)

    # A list concatenates, one LD block per element, and the variant order
    # follows the list order.
    joined <- buildLdEigen(list(one, two))
    expect_equal(length(joined), 7L)
    expect_equal(length(eigenList(joined)), 2L)
    expect_equal(eigenList(joined)[[2]]$snpIdx, 5:7)
})

test_that("the builders reject anything but LdMatrix blocks or a panel", {
    msg <- "must be a genotype panel, an LdMatrix, or a non-empty list"
    expect_error(buildLdEigen("not-ld"), msg)
    expect_error(buildLdEigen(list()), msg)
    expect_error(buildLdScore(list(makeTestLdMatrix(), "not-ld")), msg)
    expect_error(
        buildLdScore(makeTestLdMatrix(), ldBlocks = makeTestLdBlocks()),
        "`ldBlocks` cuts a genotype panel"
    )
})

test_that("nRef is taken from the LdMatrix blocks; disagreement is an error", {
    a <- makeTestLdMatrix(n = 4L, nSamples = 500L)
    b <- makeTestLdMatrix(n = 3L, startBp = 5000L, nSamples = 900L)

    expect_equal(ldPanelNSamples(buildLdEigen(list(a, a))), 500L)
    expect_error(buildLdEigen(list(a, b)), "differing reference panel sizes")
    # An explicit nSamples settles it.
    expect_equal(
        ldPanelNSamples(buildLdEigen(list(a, b), nSamples = 700L)),
        700L
    )
})

test_that("nSamples is required when no LdMatrix records one", {
    unknown <- makeTestLdMatrix(n = 3L, nSamples = NA_integer_)
    expect_error(buildLdScore(unknown), "`nSamples` is required")
    expect_equal(
        ldPanelNSamples(buildLdScore(unknown, nSamples = 250L)),
        250L
    )
    # Alongside a block that does record one, that size is used.
    known <- makeTestLdMatrix(n = 3L, startBp = 5000L, nSamples = 900L)
    expect_equal(ldPanelNSamples(buildLdEigen(list(unknown, known))), 900L)
})

test_that("the genome build comes from the reference when none is named", {
    ld <- makeTestLdMatrix()
    # An unlabelled reference stays NA rather than becoming a made-up default.
    expect_true(is.na(unname(GenomeInfoDb::genome(buildLdEigen(ld)))))

    tagged <- ld
    GenomeInfoDb::genome(tagged) <- "hg38"
    expect_equal(unname(GenomeInfoDb::genome(buildLdEigen(tagged))), "hg38")
    # An explicit argument still wins.
    expect_equal(
        unname(GenomeInfoDb::genome(buildLdEigen(tagged, genome = "hg19"))),
        "hg19"
    )
})

test_that("LD block ranges span the variants each block covers", {
    ref <- buildLdEigen(makeTestLdMatrixBlocks(sizes = c(4L, 3L)))
    blocks <- ldBlocks(ref)
    variants <- GenomicRanges::start(ref)

    expect_equal(
        as.character(GenomicRanges::seqnames(blocks)),
        c("chr1", "chr1")
    )
    expect_equal(GenomicRanges::start(blocks), c(variants[[1]], variants[[5]]))
    expect_equal(GenomicRanges::end(blocks), c(variants[[4]], variants[[7]]))
})

test_that("an LD block spanning chromosomes is rejected", {
    tbl <- .testLdVariants(4L, "chr1", 1000L)
    tbl$chrom <- c("chr1", "chr1", "chr2", "chr2")
    crossChrom <- LdMatrix(.testArBlock(4L, 0.5), variants = tbl, nSamples = 50)
    expect_error(buildLdEigen(crossChrom), "spans 2 chromosomes")
})

test_that("an AF column becomes a folded MAF column; none leaves it off", {
    tbl <- .testLdVariants(3L, "chr1", 1000L)
    tbl$AF <- c(0.1, 0.9, 0.3)
    got <- S4Vectors::mcols(
        buildLdScore(LdMatrix(diag(3), variants = tbl, nSamples = 50)),
        use.names = FALSE
    )
    expect_equal(as.numeric(got$MAF), c(0.1, 0.1, 0.3))

    tbl$AF <- NULL
    bare <- S4Vectors::mcols(
        buildLdScore(LdMatrix(diag(3), variants = tbl, nSamples = 50)),
        use.names = FALSE
    )
    expect_false(is_in("MAF", colnames(bare)))
    expect_identical(as.character(bare$SNP), tbl$SNP)
})

test_that("a genotype panel is cut into LD blocks at ldBlocks", {
    panel <- readGenotypes(test_path("test_data", "test_variants.bed"))
    pos <- GenomicRanges::start(SummarizedExperiment::rowRanges(panel))
    cut <- pos[[200L]]
    # The blocks may name chromosomes either way; a panel's are
    # chr-prefixed. Variants outside every block are left out.
    blocks <- GenomicRanges::GRanges(
        "21",
        IRanges::IRanges(c(pos[[1L]], cut), c(cut - 1L, pos[[300L]]))
    )
    ref <- buildLdScore(panel, ldBlocks = blocks)
    expect_equal(length(ref), 300L)
    expect_equal(ldPanelNSamples(ref), ncol(panel))
    # It keeps the panel and each variant's row in it, not the blocks' LD...
    expect_identical(ldReference(ref), panel)
    expect_identical(S4Vectors::mcols(ref)$panelRow, 1:300)
    expect_length(methods::slot(ref, "ldMatrixList"), 0L)
    # ...and computes a block's LD from the panel when asked for it.
    expect_equal(length(ldMatrixList(ref)), 2L)
    expect_equal(
        unname(ldMatrixList(ref)[[2]]$R),
        unname(computeLd(panel, snpIdx = 200:300))
    )
    expect_error(buildLdScore(panel), "`ldBlocks` is required")
    overlapping <- GenomicRanges::GRanges(
        "chr21",
        IRanges::IRanges(
            c(pos[[1L]], pos[[100L]]),
            c(pos[[150L]], pos[[300L]])
        )
    )
    expect_error(
        buildLdEigen(panel, ldBlocks = overlapping),
        "`ldBlocks` overlap"
    )
})

test_that("a stochastic sketch's blocks carry its effective size", {
    skip_if_not_installed("pgenlibr")
    sketch <- readGenotypes(
        test_path("test_data", "test_harmonize_regions.pgen")
    )
    gr <- SummarizedExperiment::rowRanges(sketch)
    blocks <- GenomicRanges::GRanges(
        "chr21",
        IRanges::IRanges(
            min(GenomicRanges::start(gr)),
            max(GenomicRanges::start(gr))
        )
    )
    ref <- buildLdEigen(sketch, ldBlocks = blocks)
    # 50 real samples projected onto 100 pseudo-samples: 1 / (1/50 + 1/100).
    expect_equal(ldPanelNSamples(ref), 33L)
    # The sketch's real allele frequencies are carried as MAF.
    af <- SummarizedExperiment::rowRanges(sketch)$AF
    expect_equal(as.numeric(S4Vectors::mcols(ref)$MAF), pmin(af, 1 - af))
})

test_that("a statistic read from a panel must name each variant's panel row", {
    fx <- .h2PanelFixture()
    ref <- buildLdScore(fx$panel, ldBlocks = fx$blocks)
    expect_true(validObject(ref))
    offPanel <- ref
    S4Vectors::mcols(offPanel)$panelRow <- ref$panelRow + 1000L
    expect_error(validObject(offPanel), "needs an integer `panelRow`")
    unread <- ref
    S4Vectors::mcols(unread)$panelRow <- NULL
    expect_error(validObject(unread), "needs an integer `panelRow`")
    # An LdMatrix reference's LD is stored per block, never read on demand.
    asMatrix <- ref
    asMatrix@ldReference <- fx$ld[[1L]]
    expect_error(
        validObject(asMatrix),
        "must be a genotype panel or a GenotypeMixture"
    )
})

test_that("a block's LD and decomposition are read as stored", {
    fx <- .h2PanelFixture()
    idx <- 61:200
    want <- unname(computeLd(fx$panel, snpIdx = idx))
    stored <- list(idx = seq_along(idx), R = want)
    lazy <- list(idx = seq_along(idx), rows = idx)
    expect_equal(pecotmr:::.ldStatBlockLd(stored, NULL), list(R = want))
    F <- pecotmr:::.ldStatBlockLd(lazy, fx$panel)$F
    expect_equal(unname(crossprod(F)) / (nrow(F) - 1), want)
    # 140 variants and 100 samples: the SVD lists 100 directions and counts
    # the other 40, whose eigenvalue is exactly zero, unless they would be
    # truncated away anyway.
    dec <- pecotmr:::.ldStatBlockEigen(lazy, fx$panel, 1)
    expect_equal(length(dec$values), 100L)
    expect_equal(dec$nNull, 40L)
    expect_equal(dec$snpIdx, seq_along(idx))
    full <- eigen(want, symmetric = TRUE)$values
    expect_equal(dec$values[1:90], full[1:90])
    expect_equal(
        unname(dec$vectors %*% (dec$values * t(dec$vectors))),
        want
    )
    cut <- pecotmr:::.ldStatBlockEigen(lazy, fx$panel, 0.9)
    expect_equal(cut$nNull, 0L)
    expect_equal(
        length(cut$values),
        length(pecotmr:::.ldEigenKeep(full, 0.9))
    )
    storedEigen <- list(idx = 1:3, eigen = list(values = 1, vectors = 1))
    expect_equal(
        pecotmr:::.ldStatBlockEigen(storedEigen, NULL, 1)$nNull,
        0L
    )
})
