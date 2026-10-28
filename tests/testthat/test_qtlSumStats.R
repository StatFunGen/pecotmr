context("qtlSumStats")

# ===========================================================================
# Constructor / validity
# ===========================================================================

test_that("QtlSumStats: minimal single-tuple object builds and validates", {
    obj <- .qtlMakeOne()
    expect_s4_class(obj, "QtlSumStats")
    expect_true(methods::validObject(obj))
    expect_equal(nrow(obj), 1L)
})

test_that("QtlSumStats: multi-tuple object is keyed by (study, context, trait)", {
    obj <- QtlSumStats(
        studyName = c("s1", "s1", "s2"),
        context = c("c1", "c2", "c1"),
        trait = c("t1", "t1", "t1"),
        entry = list(
            .qtlMakeEntryGr(3),
            .qtlMakeEntryGr(4),
            .qtlMakeEntryGr(5)
        ),
        genome = "hg38",
        ldReference = .qtlMakeGenotypeHandle()
    )
    expect_equal(nrow(obj), 3L)
    expect_setequal(contexts(obj), c("c1", "c2"))
    expect_equal(traitNames(obj), "t1")
})

test_that("QtlSumStats: traitPosition is empty when no position was supplied", {
    obj <- .qtlMakeOne(trait = "ENSG001")
    # traitPos cannot be inferred from summary statistics, so an unsupplied
    # position is an empty GRanges rather than a fabricated range.
    expect_s4_class(traitPosition(obj), "GRanges")
    expect_length(traitPosition(obj), 0L)
    expect_length(traitPosition(obj, "ENSG001"), 0L)
})

test_that("QtlSumStats: traitPosition returns the supplied trait position", {
    tpos <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1000L, 1500L))
    obj <- QtlSumStats(
        studyName = "s1",
        context = "c1",
        trait = "ENSGX",
        entry = list(.qtlMakeEntryGr(5)),
        genome = "hg19",
        ldReference = .qtlMakeGenotypeHandle(),
        traitPos = tpos
    )
    tp <- traitPosition(obj)
    expect_s4_class(tp, "GRanges")
    expect_equal(GenomicRanges::start(tp), 1000L)
    expect_equal(GenomicRanges::end(tp), 1500L)
    expect_equal(GenomicRanges::start(traitPosition(obj, "ENSGX")), 1000L)
    # a trait not in the collection -> an empty GRanges
    expect_s4_class(traitPosition(obj, "nope"), "GRanges")
    expect_length(traitPosition(obj, "nope"), 0L)
})

test_that("QtlSumStats: errors when required args are missing", {
    expect_error(QtlSumStats(studyName = "s1", context = "c1"), "all required")
})

test_that("QtlSumStats: errors on length mismatch among study/context/trait/entry", {
    expect_error(
        QtlSumStats(
            studyName = c("s1", "s2"),
            context = c("c1"),
            trait = c("t1"),
            entry = list(.qtlMakeEntryGr(1)),
            genome = "hg19",
            ldReference = .qtlMakeGenotypeHandle()
        ),
        "same length"
    )
})

test_that("QtlSumStats: errors when entry is not a list", {
    expect_error(
        QtlSumStats(
            studyName = "s1",
            context = "c1",
            trait = "t1",
            entry = "not_a_list",
            genome = "hg19",
            ldReference = .qtlMakeGenotypeHandle()
        ),
        "must be a list"
    )
})

test_that("QtlSumStats: errors on non-singleton genome", {
    expect_error(
        QtlSumStats(
            studyName = "s1",
            context = "c1",
            trait = "t1",
            entry = list(.qtlMakeEntryGr(1)),
            genome = c("hg19", "hg38"),
            ldReference = .qtlMakeGenotypeHandle()
        ),
        "single character string"
    )
})

test_that("QtlSumStats: errors on duplicate (study, context, trait) tuple", {
    expect_error(
        QtlSumStats(
            studyName = c("s1", "s1"),
            context = c("c1", "c1"),
            trait = c("t1", "t1"),
            entry = list(.qtlMakeEntryGr(1), .qtlMakeEntryGr(1)),
            genome = "hg19",
            ldReference = .qtlMakeGenotypeHandle()
        ),
        "uniqueness violated"
    )
})

test_that("QtlSumStats: errors when entry list contains non-GRanges", {
    expect_error(
        QtlSumStats(
            studyName = "s1",
            context = "c1",
            trait = "t1",
            entry = list("not_a_granges"),
            genome = "hg19",
            ldReference = .qtlMakeGenotypeHandle()
        ),
        "must be a GRanges"
    )
})

test_that("QtlSumStats: scalar varY is recycled across tuples", {
    obj <- QtlSumStats(
        studyName = c("s1", "s2"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        entry = list(.qtlMakeEntryGr(2), .qtlMakeEntryGr(2)),
        genome = "hg19",
        ldReference = .qtlMakeGenotypeHandle(),
        varY = 1.5
    )
    expect_equal(obj$varY, c(1.5, 1.5))
})

test_that("QtlSumStats: errors when varY or nSample length is not 1 or n", {
    for (arg in c("varY", "nSample")) {
        expect_error(
            exec(
                QtlSumStats,
                studyName = c("s1", "s2"),
                context = c("c1", "c1"),
                trait = c("t1", "t1"),
                entry = list(.qtlMakeEntryGr(2), .qtlMakeEntryGr(2)),
                genome = "hg19",
                ldReference = .qtlMakeGenotypeHandle(),
                !!arg := c(1, 2, 3)
            ),
            "length 1 or length\\(study\\)",
            info = arg
        )
    }
})

test_that("QtlSumStats: optional nSample column is absent by default, attached and recycled when supplied", {
    bare <- .qtlMakeOne()
    expect_false("nSample" %in% names(bare)) # schema unchanged by default
    obj <- QtlSumStats(
        studyName = c("s1", "s2"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        entry = list(.qtlMakeEntryGr(2), .qtlMakeEntryGr(2)),
        genome = "hg19",
        ldReference = .qtlMakeGenotypeHandle(),
        nSample = 838
    ) # scalar recycled to both tuples
    expect_equal(as.numeric(obj$nSample), c(838, 838))
    expect_true(methods::validObject(obj))
})

test_that("QtlSumStats: keepSeqlevels preserves the nSample column", {
    gr <- .qtlMakeEntryGr(4, chr = "chr1")
    obj <- QtlSumStats(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .qtlMakeGenotypeHandle(),
        nSample = 838
    )
    sub <- GenomeInfoDb::keepSeqlevels(obj, "chr1", pruning.mode = "coarse")
    expect_equal(as.numeric(sub$nSample), 838)
})

test_that("QtlSumStats: accepts and stores extra per-tuple columns via ...", {
    obj <- QtlSumStats(
        studyName = c("s1", "s2"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        entry = list(.qtlMakeEntryGr(2), .qtlMakeEntryGr(2)),
        genome = "hg19",
        ldReference = .qtlMakeGenotypeHandle(),
        extraCols = list(cohort = c("UKB", "FinnGen"))
    )
    expect_equal(as.character(obj$cohort), c("UKB", "FinnGen"))
    # Not variadic: the base class's GRangesList(...) means ELEMENTS,
    # so a bare `cohort =` here would read as adding an entry.
    expect_error(
        QtlSumStats(
            studyName = "s1",
            context = "c1",
            trait = "t1",
            entry = list(),
            genome = "hg19",
            cohort = "UKB"
        ),
        "unused argument"
    )
})

# ===========================================================================
# Accessors
# ===========================================================================

test_that("sumStats / z / sampleSize read the single-tuple entry", {
    obj <- .qtlMakeOne(n = 4)
    expect_s4_class(sumStats(obj), "GRanges")
    expect_equal(length(sumStats(obj)), 4L)
    expect_equal(z(obj), seq(1.0, by = 0.5, length.out = 4))
    expect_equal(sampleSize(obj), rep(1000L, 4))
})

test_that("nSnps reports the number of variants in the selected entry", {
    obj <- .qtlMakeOne(n = 7)
    expect_equal(nSnps(obj), 7L)
})

test_that("varY returns numeric value or NULL when NA", {
    obj_na <- .qtlMakeOne(varY = NA_real_)
    expect_null(varY(obj_na))

    obj_v <- .qtlMakeOne(varY = 2.5)
    expect_equal(varY(obj_v), 2.5)
})

test_that("contexts / traitNames return unique values across tuples", {
    obj <- QtlSumStats(
        studyName = c("s1", "s1", "s2"),
        context = c("c1", "c2", "c1"),
        trait = c("t1", "t1", "t2"),
        entry = list(
            .qtlMakeEntryGr(1),
            .qtlMakeEntryGr(1),
            .qtlMakeEntryGr(1)
        ),
        genome = "hg19",
        ldReference = .qtlMakeGenotypeHandle()
    )
    expect_setequal(contexts(obj), c("c1", "c2"))
    expect_setequal(traitNames(obj), c("t1", "t2"))
})

test_that("accessors require (study, context, trait) when collection has >1 entry", {
    obj <- QtlSumStats(
        studyName = c("s1", "s2"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        entry = list(.qtlMakeEntryGr(1), .qtlMakeEntryGr(1)),
        genome = "hg19",
        ldReference = .qtlMakeGenotypeHandle()
    )
    expect_error(
        sumStats(obj),
        "Pass `studyName`, `context`, and `trait` to select one"
    )

    expect_s4_class(
        sumStats(obj, studyName = "s1", context = "c1", trait = "t1"),
        "GRanges"
    )
})

test_that("accessors error on unknown tuple", {
    obj <- .qtlMakeOne(studyName = "s1", context = "c1", trait = "t1")
    expect_error(
        sumStats(obj, studyName = "ghost", context = "c1", trait = "t1"),
        "No entry"
    )
})

test_that("accessors require length-1 selection args", {
    obj <- .qtlMakeOne()
    expect_error(
        sumStats(obj, studyName = c("s1", "s2"), context = "c1", trait = "t1"),
        "Must have length 1"
    )
})

# ===========================================================================
# Narrowing with the standard verbs
# ===========================================================================

test_that("keepSeqlevels keeps only the requested chromosome's rows", {
    # Multi-chromosome entry: 3 on chr1, 2 on chr2. Built in one shot so
    # the resulting GRanges already shares seqlevels across rows.
    gr <- GenomicRanges::GRanges(
        seqnames = c(rep("chr1", 3), rep("chr2", 2)),
        ranges = IRanges::IRanges(
            start = c(100L, 200L, 300L, 1000L, 1100L),
            width = 1L
        )
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = paste0("rs", seq_along(gr)),
        A1 = rep("A", length(gr)),
        A2 = rep("G", length(gr)),
        Z = seq(1.0, by = 0.5, length.out = length(gr)),
        N = rep(1000L, length(gr))
    )
    obj <- QtlSumStats(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .qtlMakeGenotypeHandle()
    )

    chr1Only <- GenomeInfoDb::keepSeqlevels(
        obj,
        "chr1",
        pruning.mode = "coarse"
    )
    expect_equal(nSnps(chr1Only), 3L)
    chr2Only <- GenomeInfoDb::keepSeqlevels(
        obj,
        "chr2",
        pruning.mode = "coarse"
    )
    expect_equal(nSnps(chr2Only), 2L)
})

test_that("a range subset that keeps no row leaves an empty collection", {
    obj <- .qtlMakeOne(n = 4)
    # A window on the collection's own chromosome, past its variants.
    empty <- IRanges::subsetByOverlaps(
        obj,
        GenomicRanges::GRanges("chr1:90000000-90000100")
    )
    expect_equal(nrow(empty), 0L)
    expect_equal(length(empty), 0L)
    # The collection-level state survives even when nothing is left; no row
    # is left to use the reference, so none is reported.
    expect_equal(
        unname(GenomeInfoDb::genome(empty)),
        unname(GenomeInfoDb::genome(obj))
    )
    expect_identical(empty@ldReference, obj@ldReference)
    expect_identical(ldReference(empty), list())
})

test_that("keepSeqlevels keeps the collection-level slots", {
    obj <- .qtlMakeOne(qcOptions = list(step1 = "ok"))
    res <- GenomeInfoDb::keepSeqlevels(obj, "chr1", pruning.mode = "coarse")
    expect_equal(
        unname(GenomeInfoDb::genome(res)),
        unname(GenomeInfoDb::genome(obj))
    )
    expect_equal(qcOptions(res), qcOptions(obj))
    expect_identical(ldReference(res), ldReference(obj))
})

# ===========================================================================
# Show method
# ===========================================================================

test_that("show prints entry count and genome build", {
    obj <- .qtlMakeOne()
    expect_output(show(obj), "QtlSumStats: 1 entries, genome build hg19")
})

# ===========================================================================
# Empty-collection + duplicate-tuple selector error branches
# ===========================================================================

test_that("sumStats errors on an empty QtlSumStats collection", {
    empty <- QtlSumStats(
        studyName = character(0),
        context = character(0),
        trait = character(0),
        entry = list(),
        genome = "hg19",
        ldReference = .qtlMakeGenotypeHandle(),
        varY = numeric(0)
    )
    expect_equal(nrow(empty), 0L)
    expect_error(sumStats(empty), "has no rows")
})

test_that("QtlSumStats: ldReference is optional (NULL for LD-free workflows)", {
    obj <- QtlSumStats(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        entry = list(.qtlMakeEntryGr(1)),
        genome = "hg19"
    ) # ldReference omitted
    expect_identical(ldReference(obj), list())
    expect_true(is.na(obj$ldReference))
    expect_output(show(obj), "none \\(LD-free\\)")
})

test_that("QtlSumStats: a non-GenotypeHandle ldReference is rejected", {
    expect_error(
        QtlSumStats(
            studyName = "s1",
            context = "c1",
            trait = "t1",
            entry = list(.qtlMakeEntryGr(1)),
            genome = "hg19",
            ldReference = "not_a_handle"
        ),
        "must be a genotype panel"
    )
})

# ===========================================================================
# First-class summary-statistic accessors: pval / coef / se
# ===========================================================================

test_that("pval / coef / se read the optional P/BETA/SE mcols", {
    gr <- .qtlMakeEntryGr(4)
    S4Vectors::mcols(gr)$P <- c(0.1, 0.01, 1e-4, 0.5)
    S4Vectors::mcols(gr)$BETA <- c(0.2, -0.3, 0.4, 0.05)
    S4Vectors::mcols(gr)$SE <- rep(0.1, 4)
    x <- QtlSumStats(
        studyName = "s",
        context = "c",
        trait = "g",
        entry = list(gr),
        genome = "hg19"
    )
    expect_equal(pval(x), c(0.1, 0.01, 1e-4, 0.5))
    expect_equal(coef(x), c(0.2, -0.3, 0.4, 0.05))
    expect_equal(se(x), rep(0.1, 4))
})

test_that("pval / coef / se return NULL when the entry omits them", {
    y <- .qtlMakeOne(n = 4) # Z/N-only entry
    expect_null(pval(y))
    expect_null(coef(y))
    expect_null(se(y))
    expect_equal(z(y), seq(1.0, by = 0.5, length.out = 4)) # Z still works
})

# ===========================================================================
# Constructor auto-fills tss_distance / tes_distance from traitPos
# ===========================================================================

test_that("constructor derives tss/tes distance from a point trait position", {
    gr <- .qtlMakeEntryGr(5, start_at = 100L, step = 100L) # 100..500
    tp <- GenomicRanges::GRanges("chr1", IRanges::IRanges(250L, width = 1L))
    x <- QtlSumStats(
        studyName = "s",
        context = "c",
        trait = "g",
        entry = list(gr),
        genome = "hg19",
        traitPos = tp
    )
    mc <- S4Vectors::mcols(sumStats(x))
    expect_equal(mc$tss_distance, c(100, 200, 300, 400, 500) - 250)
    # 1bp trait -> tss and tes distance are identical
    expect_equal(mc$tss_distance, mc$tes_distance)
})

test_that("constructor derives distinct tss/tes distance for a multi-base trait", {
    gr <- .qtlMakeEntryGr(3, start_at = 100L, step = 100L) # 100,200,300
    tp <- GenomicRanges::GRanges("chr1", IRanges::IRanges(200L, 400L)) # TSS 200, TES 400
    x <- QtlSumStats(
        studyName = "s",
        context = "c",
        trait = "g",
        entry = list(gr),
        genome = "hg19",
        traitPos = tp
    )
    mc <- S4Vectors::mcols(sumStats(x))
    expect_equal(mc$tss_distance, c(100, 200, 300) - 200)
    expect_equal(mc$tes_distance, c(100, 200, 300) - 400)
})

test_that("no traitPos -> no distance columns; existing distance is preserved", {
    z <- .qtlMakeOne(n = 3)
    expect_null(S4Vectors::mcols(sumStats(z))$tss_distance)

    gr <- .qtlMakeEntryGr(2, start_at = 100L, step = 100L)
    S4Vectors::mcols(gr)$tss_distance <- c(-999L, -888L) # pre-existing
    w <- QtlSumStats(
        studyName = "s",
        context = "c",
        trait = "g",
        entry = list(gr),
        genome = "hg19",
        traitPos = GenomicRanges::GRanges(
            "chr1",
            IRanges::IRanges(250L, width = 1L)
        )
    )
    mc <- S4Vectors::mcols(sumStats(w))
    expect_equal(mc$tss_distance, c(-999L, -888L)) # not clobbered
    expect_equal(mc$tes_distance, c(100, 200) - 250) # absent -> filled
})


# === combineQtlSumStats() ===

test_that("combineQtlSumStats() row-binds collections and LD-free rows", {
    data(qtlSumStatsExample)
    a <- qtlSumStatsExample
    b <- a
    mcols(b)$study <- paste0(mcols(a)$study, "_b")

    out <- combineQtlSumStats(a, b)
    expect_s4_class(out, "QtlSumStats")
    expect_equal(nrow(out), nrow(a) + nrow(b))
    expect_identical(unname(as.list(out)), unname(c(as.list(a), as.list(b))))
    # The combined rows are named from their own tuples; mcols<- left `b`
    # with the names of the study it was copied from.
    expect_identical(
        names(out),
        c("study1|brain|ENSG_example", "study1_b|brain|ENSG_example")
    )
    expect_equal(
        unname(GenomeInfoDb::genome(out)),
        unname(GenomeInfoDb::genome(a))
    )

    # Both parts read the one panel, so it is stored once for every row.
    expect_identical(ldReference(out), ldReference(a))
    expect_equal(unique(out$ldReference), names(ldReference(out)))

    ldFree <- function(x) {
        QtlSumStats(
            studyName = x$study,
            context = x$context,
            trait = x$trait,
            entry = as.list(x),
            genome = unname(GenomeInfoDb::genome(x))[[1L]]
        )
    }
    expect_identical(
        ldReference(combineQtlSumStats(ldFree(a), ldFree(b))),
        list()
    )
})


test_that("combineQtlSumStats() passes a lone input through and validates", {
    data(qtlSumStatsExample)
    expect_identical(combineQtlSumStats(qtlSumStatsExample), qtlSumStatsExample)
    expect_error(combineQtlSumStats(list()), "nothing to combine")
    expect_error(
        combineQtlSumStats(qtlSumStatsExample, 1L),
        "every input must be a QtlSumStats"
    )
})


# ===========================================================================
# Trait-distance annotation and validity messages
# ===========================================================================

test_that("trait distances are added to an entry that carries no mcols", {
    # Regression: mcols() on a bare GRanges is a zero-COLUMN DataFrame that
    # still has one row per variant. Swapping it for DataFrame(row.names =
    # NULL) discarded the row count and the assignment died with
    # "2 elements in value to replace 0 elements".
    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(10L, 20L), width = 1L)
    )
    tp <- GenomicRanges::GRanges("chr1", IRanges::IRanges(5L, width = 1L))
    out <- pecotmr:::.qssAppendTraitDist(1L, list(gr), tp)
    expect_equal(S4Vectors::mcols(out)$tss_distance, c(5L, 15L))
    expect_equal(S4Vectors::mcols(out)$tes_distance, c(5L, 15L))
})

test_that("existing mcols survive the trait-distance annotation", {
    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(10L, 20L), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(z = c(1.5, 2.5))
    tp <- GenomicRanges::GRanges("chr1", IRanges::IRanges(5L, width = 1L))
    out <- pecotmr:::.qssAppendTraitDist(1L, list(gr), tp)
    expect_equal(
        colnames(S4Vectors::mcols(out)),
        c("z", "tss_distance", "tes_distance")
    )
    expect_equal(S4Vectors::mcols(out)$z, c(1.5, 2.5))
})

test_that("an entry with no variants is returned untouched", {
    # A variant set with no ranges has no positions to take a distance from.
    gr0 <- GenomicRanges::GRanges()
    tp <- GenomicRanges::GRanges("chr1", IRanges::IRanges(5L, width = 1L))
    out <- pecotmr:::.qssAppendTraitDist(1L, list(gr0), tp)
    expect_s4_class(out, "GRanges")
    expect_length(out, 0L)
})

test_that("validity names the identity columns that are missing", {
    data(qtlSumStatsExample)
    bad <- qtlSumStatsExample
    S4Vectors::mcols(bad)$trait <- NULL
    expect_error(
        methods::validObject(bad),
        "missing elements \\{'trait'\\}"
    )
    expect_true(methods::validObject(qtlSumStatsExample))
})

test_that(".appendTraitDistances is a no-op without a trait position", {
    entry <- list(
        GenomicRanges::GRanges("chr1", IRanges::IRanges(100, width = 1))
    )
    # No anchor to measure TSS/TES distance from, so entries pass through.
    expect_identical(pecotmr:::.appendTraitDistances(entry, NULL), entry)
})

test_that("QtlSumStats records each tuple's window and its origin", {
    data(qtlSumStatsExample)
    q <- qtlSumStatsExample
    g <- q[[1L]]
    chrom <- as.character(GenomicRanges::seqnames(g))[[1L]]
    lo <- GenomicRanges::start(g)[[10L]]
    hi <- GenomicRanges::start(g)[[150L]]
    build <- function(...) {
        QtlSumStats(
            studyName = q$study,
            context = q$context,
            trait = q$trait,
            entry = list(g),
            genome = unique(GenomeInfoDb::genome(q)),
            ...
        )
    }
    # Without one, a tuple's window is the span it was built with.
    plain <- build()
    expect_false(plain$regionGiven)
    expect_identical(as.character(plain$region), as.character(range(g)))
    # A window drops the variants outside it, with a warning, and is given.
    window <- str_c(chrom, ":", lo, "-", hi)
    expect_warning(
        inside <- build(region = GenomicRanges::GRanges(window)),
        "variant(s) lie outside their tuple's `region`",
        fixed = TRUE
    )
    expect_true(inside$regionGiven)
    expect_identical(as.character(inside$region), window)
    expect_identical(length(inside[[1L]]), 141L)
    # Read from a BED file, its start is 0-based.
    bedPath <- withr::local_tempfile(fileext = ".bed")
    writeLines(str_c(chrom, "\t", lo - 1L, "\t", hi), bedPath)
    bed <- suppressWarnings(build(region = bedPath))
    expect_identical(as.character(bed$region), window)
    # One window for every tuple, or one per tuple.
    expect_error(
        build(region = GenomicRanges::GRanges(c(window, window))),
        "one window, or one per tuple"
    )
})

test_that("a cisWindow narrows each selected tuple to its trait's window", {
    data(qtlSumStatsExample)
    q <- qtlSumStatsExample
    g <- q[[1L]]
    chrom <- as.character(GenomicRanges::seqnames(g))[[1L]]
    tp <- GenomicRanges::GRanges(
        chrom,
        IRanges::IRanges(GenomicRanges::start(g)[[80L]], width = 1L)
    )
    build <- function(...) {
        QtlSumStats(
            studyName = q$study,
            context = q$context,
            trait = q$trait,
            entry = list(g),
            genome = unique(GenomeInfoDb::genome(q)),
            ...
        )
    }
    x <- build(traitPos = tp)
    expect_identical(pecotmr:::.qssApplyCisWindow(x, NULL, "f", 1L), x)
    cis <- pecotmr:::.qssApplyCisWindow(x, 5000, "f", 1L)
    expect_true(cis$regionGiven)
    expect_identical(
        GenomicRanges::start(cis$region),
        GenomicRanges::start(tp) - 5000L
    )
    expect_true(all(
        abs(GenomicRanges::start(cis[[1L]]) - GenomicRanges::start(tp)) <= 5000L
    ))
    # The cis window keeps inside the window the tuple already had.
    edge <- GenomicRanges::GRanges(
        chrom,
        IRanges::IRanges(
            GenomicRanges::start(tp),
            GenomicRanges::end(g)[[200L]]
        )
    )
    clipped <- pecotmr:::.qssApplyCisWindow(
        suppressWarnings(build(traitPos = tp, region = edge)),
        5000,
        "f",
        1L
    )
    expect_identical(
        GenomicRanges::start(clipped$region),
        GenomicRanges::start(tp)
    )
    # A selected tuple with no trait position has no cis window; one the
    # selection leaves out is left alone.
    bare <- build()
    expect_error(
        pecotmr:::.qssApplyCisWindow(
            bare,
            5000,
            "fineMappingPipeline",
            1L
        ),
        "fineMappingPipeline: `cisWindow` is laid around each tuple's trait",
        fixed = TRUE
    )
    expect_identical(
        pecotmr:::.qssApplyCisWindow(bare, 5000, "f", integer(0)),
        bare
    )
})

test_that(".clipWindows keeps nothing of a window that misses its row's", {
    windows <- GenomicRanges::GRanges(
        c("chr1", "chr1", "chr2"),
        IRanges::IRanges(c(100L, 100L, 100L), c(200L, 200L, 200L))
    )
    regions <- GenomicRanges::GRanges(
        c("chr1", "chr1", "chr1"),
        IRanges::IRanges(c(150L, 300L, 100L), c(400L, 400L, 200L))
    )
    clipped <- pecotmr:::.clipWindows(windows, regions)
    expect_identical(as.character(clipped)[[1L]], "chr1:150-200")
    expect_identical(GenomicRanges::width(clipped)[2:3], c(0L, 0L))
})

# ===========================================================================
# metadata, qcOptions and a tabular extraCols
# ===========================================================================

test_that("QtlSumStats() takes extra columns as a DataFrame", {
    obj <- QtlSumStats(
        studyName = c("s1", "s2"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        entry = list(.qtlMakeEntryGr(2), .qtlMakeEntryGr(2)),
        genome = "hg19",
        ldReference = .qtlMakeGenotypeHandle(),
        extraCols = S4Vectors::DataFrame(cohort = c("UKB", "FinnGen"))
    )
    expect_equal(obj$cohort, c("UKB", "FinnGen"))
})

test_that("QtlSumStats() keeps the caller's metadata and QC options", {
    obj <- .qtlMakeOne(qcOptions = list(step1 = "ok"))
    expect_equal(qcOptions(obj), list(step1 = "ok"))
    withMeta <- QtlSumStats(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        entry = list(.qtlMakeEntryGr(2)),
        genome = "hg19",
        ldReference = .qtlMakeGenotypeHandle(),
        metadata = list(source = "unit test")
    )
    expect_equal(S4Vectors::metadata(withMeta), list(source = "unit test"))
})

# ===========================================================================
# sumStats(annotateSignificance =)
# ===========================================================================

test_that("sumStats(annotateSignificance=) adds a derived logical mcol", {
    x <- .qapDefaultFixture
    r <- qtlAssociationPostprocess(x, methods = "permutation")
    gr <- sumStats(
        r,
        studyName = "s",
        context = "brain",
        trait = "g1",
        annotateSignificance = "permutation"
    )
    expect_true("significant" %in% names(S4Vectors::mcols(gr)))
    expect_type(S4Vectors::mcols(gr)$significant, "logical")
    # matches the permutation rule P < p_nominal_threshold[gene 1]
    thr1 <- as.numeric(r$p_nominal_threshold)[1]
    exp1 <- if (is.na(thr1)) {
        rep(FALSE, length(gr))
    } else {
        S4Vectors::mcols(gr)$P < thr1
    }
    expect_equal(S4Vectors::mcols(gr)$significant, exp1)
})
