# Tests for R/MultiStudyQtlDataset.R

test_that("MultiStudyQtlDataset: combines two QtlDatasets", {
    qd1 <- QtlDataset(
        studyName = "s1",
        genotypes = .sc_makeGenotypeHandle(),
        phenotypes = list(brain = .sc_makeSe())
    )
    qd2 <- QtlDataset(
        studyName = "s2",
        genotypes = .sc_makeGenotypeHandle(),
        phenotypes = list(brain = .sc_makeSe())
    )
    mt <- MultiStudyQtlDataset(qtlDatasets = list(s1 = qd1, s2 = qd2))
    expect_s4_class(mt, "MultiStudyQtlDataset")
    expect_setequal(studyName(mt), c("s1", "s2"))
})


test_that("MultiStudyQtlDataset: rejects single dataset with no sumStats", {
    qd <- QtlDataset(
        studyName = "s1",
        genotypes = .sc_makeGenotypeHandle(),
        phenotypes = list(brain = .sc_makeSe())
    )
    expect_error(
        MultiStudyQtlDataset(qtlDatasets = list(s1 = qd)),
        "at least 2 studies"
    )
})


test_that("MultiStudyQtlDataset: rejects unnamed qtlDatasets list", {
    qd <- QtlDataset(
        studyName = "s1",
        genotypes = .sc_makeGenotypeHandle(),
        phenotypes = list(brain = .sc_makeSe())
    )
    expect_error(
        MultiStudyQtlDataset(qtlDatasets = list(qd, qd)),
        "Must have names"
    )
})


test_that("MultiStudyQtlDataset: rejects non-QtlDataset entries", {
    qd <- QtlDataset(
        studyName = "s1",
        genotypes = .sc_makeGenotypeHandle(),
        phenotypes = list(brain = .sc_makeSe())
    )
    expect_error(
        MultiStudyQtlDataset(qtlDatasets = list(s1 = qd, s2 = "not a dataset")),
        "must be a QtlDataset"
    )
})

test_that("MultiStudyQtlDataset: rejects a trait at two positions", {
    se1 <- .sc_makeSe(traits = "ENSG1")
    # Build se2 from scratch (see note in the QtlDataset trait-conflict test).
    rng2 <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 9999L, width = 500L)
    )
    names(rng2) <- "ENSG1"
    expr2 <- matrix(
        rnorm(10),
        nrow = 1,
        ncol = 10,
        dimnames = list("ENSG1", paste0("s", 1:10))
    )
    cd2 <- S4Vectors::DataFrame(
        sex = rep(c("M", "F"), 5),
        row.names = paste0("s", 1:10)
    )
    se2 <- SummarizedExperiment::SummarizedExperiment(
        assays = list(expression = expr2),
        rowRanges = rng2,
        colData = cd2
    )
    qd1 <- QtlDataset(
        studyName = "s1",
        genotypes = .sc_makeGenotypeHandle(),
        phenotypes = list(brain = se1)
    )
    qd2 <- QtlDataset(
        studyName = "s2",
        genotypes = .sc_makeGenotypeHandle(),
        phenotypes = list(brain = se2)
    )
    expect_error(
        MultiStudyQtlDataset(qtlDatasets = list(s1 = qd1, s2 = qd2)),
        "inconsistent rowRanges"
    )
})


test_that("MultiStudyQtlDataset: allows studies differing only in chr prefix", {
    # se1 uses "chr1"; se2 is the SAME locus (chr1:1000-1499) but seqnames "1".
    se1 <- .sc_makeSe(traits = "ENSG1")
    rng2 <- GenomicRanges::GRanges(
        seqnames = "1",
        ranges = IRanges::IRanges(start = 1000L, width = 500L)
    )
    names(rng2) <- "ENSG1"
    expr2 <- matrix(
        rnorm(10),
        nrow = 1,
        ncol = 10,
        dimnames = list("ENSG1", paste0("s", 1:10))
    )
    cd2 <- S4Vectors::DataFrame(
        sex = rep(c("M", "F"), 5),
        row.names = paste0("s", 1:10)
    )
    se2 <- SummarizedExperiment::SummarizedExperiment(
        assays = list(expression = expr2),
        rowRanges = rng2,
        colData = cd2
    )
    qd1 <- QtlDataset(
        studyName = "s1",
        genotypes = .sc_makeGenotypeHandle(),
        phenotypes = list(brain = se1)
    )
    qd2 <- QtlDataset(
        studyName = "s2",
        genotypes = .sc_makeGenotypeHandle(),
        phenotypes = list(brain = se2)
    )
    expect_s4_class(
        MultiStudyQtlDataset(qtlDatasets = list(s1 = qd1, s2 = qd2)),
        "MultiStudyQtlDataset"
    )
})


test_that("sumStats(MultiStudyQtlDataset) rejects selection arguments", {
    # Compose one individual-level QtlDataset with a QtlSumStats of
    # summary-statistic-only studies (1 + 1 = 2 studies total).
    gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100L, width = 1L))
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = "rs1",
        A1 = "A",
        A2 = "G",
        Z = 1.0,
        N = 1000L
    )
    ss <- QtlSumStats(
        studyName = "s3",
        context = "c1",
        trait = "t1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .sc_makeGenotypeHandle()
    )
    qd1 <- QtlDataset(
        studyName = "s1",
        genotypes = .sc_makeGenotypeHandle(),
        phenotypes = list(brain = .sc_makeSe())
    )
    mt <- MultiStudyQtlDataset(qtlDatasets = list(s1 = qd1), sumStats = ss)

    # Bare call returns the embedded QtlSumStats collection ...
    expect_s4_class(sumStats(mt), "QtlSumStats")
    # ... but any selection argument is rejected.
    # Taking no extra formals is the rejection: R names the argument.
    expect_error(sumStats(mt, studyName = "s1"), "unused argument")
})


# ===========================================================================
# Validity branches
# ===========================================================================

test_that("MultiStudyQtlDataset: rejects an empty qtlDatasets list", {
    expect_error(
        MultiStudyQtlDataset(qtlDatasets = list()),
        "qtlDatasets.*Must have length >= 1"
    )
})

test_that("MultiStudyQtlDataset: rejects duplicated study names", {
    qd <- QtlDataset(
        studyName = "s1",
        genotypes = .sc_makeGenotypeHandle(),
        phenotypes = list(brain = .sc_makeSe())
    )
    expect_error(
        MultiStudyQtlDataset(qtlDatasets = list(s1 = qd, s1 = qd)),
        "qtlDatasets.*Must have unique names"
    )
})

test_that("MultiStudyQtlDataset: rejects a non-QtlSumStats sumStats", {
    # The slot is typed QtlSumStats_OR_NULL, so S4 refuses the value when the
    # object is built.
    qd <- QtlDataset(
        studyName = "s1",
        genotypes = .sc_makeGenotypeHandle(),
        phenotypes = list(brain = .sc_makeSe())
    )
    expect_error(
        MultiStudyQtlDataset(
            qtlDatasets = list(s1 = qd, s2 = qd),
            sumStats = "not a QtlSumStats"
        ),
        'invalid object for slot "sumStats"'
    )
})


test_that("studyName unions the individual-level and summary-only studies", {
    data(qtlSumStatsExample)
    qd <- QtlDataset(
        studyName = "s1",
        genotypes = .sc_makeGenotypeHandle(),
        phenotypes = list(brain = .sc_makeSe())
    )
    mt <- MultiStudyQtlDataset(
        qtlDatasets = list(s1 = qd),
        sumStats = qtlSumStatsExample
    )
    expect_setequal(studyName(mt), c("s1", "study1"))
})

test_that("MultiStudyQtlDataset: qtlDatasets guard fires on a non-list", {
    expect_error(
        MultiStudyQtlDataset(qtlDatasets = "not-a-list"),
        "qtlDatasets.*Must be of type 'list'"
    )
})

test_that(".msqdConcat answers an empty list for no pieces", {
    expect_identical(pecotmr:::.msqdConcat(list()), list())
})

test_that("updateObject rebuilds a MultiStudyQtlDataset from its studies", {
    data(multiStudyQtlDatasetExample)
    ms <- multiStudyQtlDatasetExample
    expect_equal(updateObject(ms), ms)
    # A study set with no summary statistics keeps none: the NULL slot is
    # read through slot(), not as the placeholder attr() would return.
    expect_null(pecotmr:::.msqdSavedSumStats(ms))
    expect_null(sumStats(updateObject(ms)))
    noSlot <- ms
    attr(noSlot, "sumStats") <- NULL
    expect_null(pecotmr:::.msqdSavedSumStats(noSlot))
    expect_message(
        updateObject(ms, verbose = TRUE),
        "updateObject\\(object = 'MultiStudyQtlDataset'\\)"
    )
    # Saved before it became a SimpleList, the studies sat in a
    # `qtlDatasets` slot, which the object still carries as an attribute.
    old <- ms
    attr(old, "qtlDatasets") <- rev(as.list(ms))
    expect_identical(names(updateObject(old)), rev(names(ms)))
    # Its summary-statistic studies are updated along with it.
    qd <- QtlDataset(
        studyName = "s1",
        genotypes = .sc_makeGenotypeHandle(),
        phenotypes = list(brain = .sc_makeSe())
    )
    withSs <- MultiStudyQtlDataset(
        qtlDatasets = list(s1 = qd),
        sumStats = qtlSumStatsExample
    )
    expect_s4_class(sumStats(updateObject(withSs)), "QtlSumStats")
})

test_that("subsetting a MultiStudyQtlDataset re-checks what is left", {
    data(multiStudyQtlDatasetExample)
    ms <- multiStudyQtlDatasetExample
    swapped <- ms[rev(names(ms))]
    expect_s4_class(swapped, "MultiStudyQtlDataset")
    expect_identical(names(swapped), rev(names(ms)))
    expect_error(ms[1L], "at least 2 studies")
})

# ===========================================================================
# show()
# ===========================================================================

test_that("show.MultiStudyQtlDataset reports per-source study counts", {
    qd1 <- .sh_makeQtlDataset(studyName = "s1")
    qd2 <- .sh_makeQtlDataset(studyName = "s2")
    mt <- MultiStudyQtlDataset(qtlDatasets = list(s1 = qd1, s2 = qd2))
    out <- capture.output(show(mt))
    expect_true(any(grepl(
        "MultiStudyQtlDataset: 2 individual-level \\+ 0 sumstats",
        out
    )))
    expect_true(any(grepl("Individual-level studies: s1, s2", out)))
})

test_that("show.MultiStudyQtlDataset reports sumstats studies when present", {
    qd <- .sh_makeQtlDataset(studyName = "s1")
    ss <- QtlSumStats(
        studyName = "s2",
        context = "c1",
        trait = "t1",
        entry = list(.sh_makeQtlSumstatsGr()),
        genome = "hg19",
        ldReference = .sh_makeGenotypeHandle()
    )
    mt <- MultiStudyQtlDataset(qtlDatasets = list(s1 = qd), sumStats = ss)
    out <- capture.output(show(mt))
    expect_true(any(grepl("Sumstats studies: s2", out)))
})
