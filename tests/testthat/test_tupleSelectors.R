context("tupleSelectors (internal row-selector helpers)")

# These helpers read identity columns off `mcols(x)`, so the tests below use
# a real collection rather than a plain data.frame. `.ts_coll()` builds the
# smallest one that exercises them: one empty GRanges per row, identity
# columns in mcols, no payload class.

setClass("TsTestCollection", contains = "RangedTupleList")

.ts_coll <- function(..., stringsAsFactors = FALSE) {
    cols <- list(...)
    n <- if (length(cols) == 0L) 0L else length(cols[[1L]])
    grl <- GenomicRanges::GRangesList(
        rep(list(GenomicRanges::GRanges()), n)
    )
    if (length(cols) > 0L) {
        S4Vectors::mcols(grl) <- exec(
            S4Vectors::DataFrame,
            !!!c(cols, list(check.names = FALSE))
        )
    }
    methods::new("TsTestCollection", grl)
}

# ===========================================================================
# .matchTupleRows
# ===========================================================================

test_that(".matchTupleRows: empty keys returns every row index", {
    df <- .ts_coll(study = c("s1", "s2"), method = c("susie", "lasso"))
    expect_equal(pecotmr:::.matchTupleRows(df, list()), c(1L, 2L))
})

test_that(".matchTupleRows: AND-matches across multiple (column, value) pairs", {
    df <- .ts_coll(
        study = c("s1", "s1", "s2"),
        context = c("c1", "c2", "c1"),
        stringsAsFactors = FALSE
    )
    expect_equal(pecotmr:::.matchTupleRows(df, list(study = "s1")), c(1L, 2L))
    expect_equal(
        pecotmr:::.matchTupleRows(df, list(study = "s1", context = "c2")),
        2L
    )
    expect_equal(
        pecotmr:::.matchTupleRows(df, list(study = "ghost", context = "c1")),
        integer(0)
    )
})

# ===========================================================================
# .tupleSelectRow (QtlFineMappingResult / TwasWeights shape)
# ===========================================================================

test_that(".tupleSelectRow: zero-row input errors with the class label", {
    empty <- .ts_coll(
        study = character(0),
        context = character(0),
        trait = character(0),
        method = character(0),
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::.tupleSelectRow(
            empty,
            studyName = "s1",
            context = "c1",
            trait = "t1",
            method = "susie",
            cls = "TwasWeights"
        ),
        "TwasWeights has no rows"
    )
})

test_that(".tupleSelectRow: single-row collection returns 1L without selectors", {
    one <- .ts_coll(
        study = "s1",
        context = "c1",
        trait = "t1",
        method = "susie",
        stringsAsFactors = FALSE
    )
    expect_equal(pecotmr:::.tupleSelectRow(one), 1L)
})

test_that(".tupleSelectRow: multi-row + missing selectors errors with row count", {
    multi <- .ts_coll(
        study = c("s1", "s1"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        method = c("susie", "susie"),
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::.tupleSelectRow(multi, cls = "QtlFineMappingResult"),
        "QtlFineMappingResult has 2 entries"
    )
})

test_that(".tupleSelectRow: non-scalar selectors error", {
    multi <- .ts_coll(
        study = c("s1", "s2"),
        context = c("c1", "c2"),
        trait = c("t1", "t2"),
        method = c("susie", "susie"),
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::.tupleSelectRow(
            multi,
            studyName = c("s1", "s2"),
            context = "c1",
            trait = "t1",
            method = "susie"
        ),
        "Must have length 1"
    )
})

test_that(".tupleSelectRow: matching tuple returns first row index", {
    multi <- .ts_coll(
        study = c("s1", "s1"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        method = c("susie", "susie"),
        stringsAsFactors = FALSE
    )
    expect_equal(
        pecotmr:::.tupleSelectRow(
            multi,
            studyName = "s1",
            context = "c2",
            trait = "t1",
            method = "susie"
        ),
        2L
    )
})

test_that(".tupleSelectRow: missing tuple errors with the 4-tuple in the message", {
    multi <- .ts_coll(
        study = c("s1", "s1"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        method = c("susie", "susie"),
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::.tupleSelectRow(
            multi,
            studyName = "ghost",
            context = "c1",
            trait = "t1",
            method = "susie"
        ),
        "No entry for"
    )
})

# ===========================================================================
# .tupleSelectRowGwasFmr
# ===========================================================================

test_that(".tupleSelectRowGwasFmr: zero-row input errors", {
    empty <- .ts_coll(
        study = character(0),
        method = character(0),
        blockId = character(0),
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::.tupleSelectRowGwasFmr(
            empty,
            studyName = "g1",
            method = "susie"
        ),
        "GwasFineMappingResult has no rows"
    )
})

test_that(".tupleSelectRowGwasFmr: single-row collection returns 1L", {
    one <- .ts_coll(
        study = "g1",
        method = "susie",
        blockId = "region_1",
        stringsAsFactors = FALSE
    )
    expect_equal(pecotmr:::.tupleSelectRowGwasFmr(one), 1L)
})

test_that(".tupleSelectRowGwasFmr: missing selectors on multi-row errors", {
    multi <- .ts_coll(
        study = c("g1", "g2"),
        method = c("susie", "susie"),
        blockId = c("region_1", "region_1"),
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::.tupleSelectRowGwasFmr(multi),
        "Pass `studyName` and `method`"
    )
})

test_that(".tupleSelectRowGwasFmr: non-scalar blockId errors", {
    multi <- .ts_coll(
        study = c("g1", "g1"),
        method = c("susie", "susie"),
        blockId = c("r1", "r2"),
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::.tupleSelectRowGwasFmr(
            multi,
            studyName = "g1",
            method = "susie",
            blockId = c("r1", "r2")
        ),
        "blockId.*Must have length 1"
    )
})

test_that(".tupleSelectRowGwasFmr: blockId disambiguates per-block rows", {
    # Same (study, method) across two blocks; blockId picks the right row.
    multi <- .ts_coll(
        study = c("g1", "g1"),
        method = c("susie", "susie"),
        blockId = c("chr22_1_100", "chr22_500_600"),
        stringsAsFactors = FALSE
    )
    expect_equal(
        pecotmr:::.tupleSelectRowGwasFmr(
            multi,
            studyName = "g1",
            method = "susie",
            blockId = "chr22_500_600"
        ),
        2L
    )
})

test_that(".tupleSelectRowGwasFmr: a missing tuple's error names its blockId", {
    one <- .ts_coll(
        study = "g1",
        method = "susie",
        blockId = "r1",
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::.tupleSelectRowGwasFmr(
            one,
            studyName = "g1",
            method = "susie",
            blockId = "ghost"
        ),
        "blockId='ghost'"
    )
})

test_that(".tupleSelectRowGwasFmr: ambiguous multi-match lists the blocks", {
    # Two rows share (study, method); .tupleSelectRowGwasFmr should error
    # listing the available blockIds since the caller didn't disambiguate.
    multi <- .ts_coll(
        study = c("g1", "g1"),
        method = c("susie", "susie"),
        blockId = c("region_A", "region_B"),
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::.tupleSelectRowGwasFmr(
            multi,
            studyName = "g1",
            method = "susie"
        ),
        "pass `blockId` to disambiguate"
    )
})

# ===========================================================================
# .fmrRowsMatching (aggregate row selector -- never errors on ambiguity)
# ===========================================================================

test_that(".fmrRowsMatching: no selectors returns every row", {
    df <- .ts_coll(
        study = c("s1", "s1"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        method = c("susie", "susie"),
        stringsAsFactors = FALSE
    )
    expect_equal(pecotmr:::.fmrRowsMatching(df), c(1L, 2L))
})

test_that(".fmrRowsMatching: matches a subset without erroring on ambiguity", {
    df <- .ts_coll(
        study = c("s1", "s1", "s2"),
        context = c("c1", "c2", "c1"),
        trait = c("t1", "t1", "t1"),
        method = c("susie", "susie", "susie"),
        stringsAsFactors = FALSE
    )
    expect_equal(pecotmr:::.fmrRowsMatching(df, studyName = "s1"), c(1L, 2L))
    expect_equal(pecotmr:::.fmrRowsMatching(df, context = "c1"), c(1L, 3L))
})

test_that(".fmrRowsMatching: selectors on absent columns are ignored", {
    # GWAS-shaped frame has no context/trait column; passing context must not
    # error and must not constrain the result.
    gwas <- .ts_coll(
        study = c("g1", "g1"),
        method = c("susie", "susie"),
        blockId = c("r1", "r2"),
        stringsAsFactors = FALSE
    )
    expect_equal(pecotmr:::.fmrRowsMatching(gwas, context = "c1"), c(1L, 2L))
})

test_that(".fmrRowsMatching: `blockId` matches the blockId column", {
    gwas <- .ts_coll(
        study = c("g1", "g1"),
        method = c("susie", "susie"),
        blockId = c("r1", "r2"),
        stringsAsFactors = FALSE
    )
    expect_equal(pecotmr:::.fmrRowsMatching(gwas, blockId = "r2"), 2L)
})

test_that(".fmrRowsMatching: a vector selector matches any listed value", {
    df <- .ts_coll(
        study = c("s1", "s2", "s3"),
        context = c("c1", "c2", "c3"),
        trait = c("t1", "t1", "t1"),
        method = c("susie", "susie", "susie"),
        stringsAsFactors = FALSE
    )
    expect_equal(
        pecotmr:::.fmrRowsMatching(df, studyName = c("s1", "s3")),
        c(1L, 3L)
    )
})

# ===========================================================================
# .fmrRowMetadata (stable 8-column identity frame)
# ===========================================================================

test_that(".fmrRowMetadata: emits all eight identity columns, NA-filling absent ones", {
    qtl <- .ts_coll(
        study = c("s1", "s1"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        method = c("susie", "susie"),
        stringsAsFactors = FALSE
    )
    m <- pecotmr:::.fmrRowMetadata(qtl)
    expect_equal(
        names(m),
        c(
            "study",
            "context",
            "trait",
            "blockId",
            "method",
            "dataType",
            "dataSource",
            "standardized"
        )
    )
    expect_equal(m$context, c("c1", "c2"))
    expect_true(all(is.na(m$blockId))) # QTL frame has no blockId
})

test_that(".fmrRowMetadata: GWAS frame NA-fills context/trait, keeps blockId", {
    gwas <- .ts_coll(
        study = c("g1", "g1"),
        method = c("susie", "susie"),
        blockId = c("r1", "r2"),
        stringsAsFactors = FALSE
    )
    m <- pecotmr:::.fmrRowMetadata(gwas)
    expect_equal(m$blockId, c("r1", "r2"))
    expect_true(all(is.na(m$context)))
    expect_true(all(is.na(m$trait)))
})

test_that(".fmrRowMetadata: zero-row input yields a zero-row 8-column frame", {
    empty <- .ts_coll(
        study = character(0),
        method = character(0),
        blockId = character(0),
        stringsAsFactors = FALSE
    )
    m <- pecotmr:::.fmrRowMetadata(empty)
    expect_equal(nrow(m), 0L)
    expect_equal(
        names(m),
        c(
            "study",
            "context",
            "trait",
            "blockId",
            "method",
            "dataType",
            "dataSource",
            "standardized"
        )
    )
})

# ===========================================================================
# .rbindAligned (union-of-columns rbind)
# ===========================================================================

test_that(".rbindAligned: a single part is returned unchanged", {
    df <- data.frame(a = 1:2, b = c("x", "y"), stringsAsFactors = FALSE)
    expect_identical(pecotmr:::.rbindAligned(list(df)), df)
})

test_that(".rbindAligned: identical columns rbind and preserve type", {
    p1 <- data.frame(a = 1L, b = "x", stringsAsFactors = FALSE)
    p2 <- data.frame(a = 2L, b = "y", stringsAsFactors = FALSE)
    out <- pecotmr:::.rbindAligned(list(p1, p2))
    expect_equal(nrow(out), 2L)
    expect_type(out$a, "integer")
    expect_equal(out$b, c("x", "y"))
})

test_that(".rbindAligned: differing columns align on the union, NA-filling gaps", {
    p1 <- data.frame(a = 1L, b = "x", stringsAsFactors = FALSE)
    p2 <- data.frame(a = 2L, c = "z", stringsAsFactors = FALSE)
    out <- pecotmr:::.rbindAligned(list(p1, p2))
    expect_equal(sort(names(out)), c("a", "b", "c"))
    expect_equal(out$a, c(1L, 2L))
    expect_equal(out$b, c("x", NA))
    expect_equal(out$c, c(NA, "z"))
})

# ===========================================================================
# Coverage: region-column + rbind helpers
# ===========================================================================
test_that(".naLikeColumn: a non-character/non-GRanges exemplar pads with logical NA", {
    expect_identical(pecotmr:::.naLikeColumn(TRUE, 3L), rep(NA, 3L))
})

test_that(".rbindCollections: all-NULL input returns NULL", {
    expect_null(pecotmr:::.rbindCollections(list(NULL, NULL)))
})

test_that(".getRegionColumn: a zero-row collection yields an empty GRanges", {
    gr <- pecotmr:::.getRegionColumn(.ts_coll(study = character(0)))
    expect_s4_class(gr, "GRanges")
    expect_length(gr, 0L)
})

test_that(".appendBlockIdCol: blockId must match the row count", {
    # blockId is a label keying the external block manifest, so subsetting
    # cannot invalidate it.
    expect_equal(pecotmr:::.appendBlockIdCol(list(), NULL, 2L), list())
    expect_equal(
        pecotmr:::.appendBlockIdCol(list(), c("b1", "b2"), 2L)$blockId,
        c("b1", "b2")
    )
    expect_error(
        pecotmr:::.appendBlockIdCol(list(), "only_one", 2L),
        "same length as"
    )
})

test_that(".appendTraitPosCol: traitPos must be a GRanges of matching length", {
    expect_error(
        pecotmr:::.appendTraitPosCol(list(), "notgranges", 1L),
        "must be a GRanges"
    )
    expect_error(
        pecotmr:::.appendTraitPosCol(
            list(),
            GenomicRanges::GRanges(
                c("chr1", "chr1"),
                IRanges::IRanges(1:2, 1:2)
            ),
            1L
        ),
        "same length as"
    )
})

test_that(".validateTraitPosColumn: reports non-GRanges and wrong-length traitPos", {
    expect_equal(
        pecotmr:::.validateTraitPosColumn(.ts_coll(traitPos = c("x", "y"))),
        str_c(
            "'traitPos' column Must inherit from class 'GRanges', ",
            "but has class 'character'"
        )
    )
    # A one-range traitPos beside two rows: assigned through the mcols
    # listData because the parallel-length check would reject it otherwise,
    # which is exactly the state the validator has to catch.
    bad <- .ts_coll(a = 1:2)
    S4Vectors::mcols(bad)@listData$traitPos <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(1, 1)
    )
    expect_equal(
        pecotmr:::.validateTraitPosColumn(bad),
        "'traitPos' column must have one range per row"
    )
    expect_length(
        pecotmr:::.validateTraitPosColumn(.ts_coll(a = 1L)),
        0L
    )
})


test_that("an exemplar column absent from every part is NULL", {
    # Callers draw the name from allCols so in practice some part has it; the
    # fallthrough is the honest answer when none does, not an error.
    expect_null(pecotmr:::.exemplarColumn("nosuchcol", list()))
})


# ===========================================================================
# Combining collections: QC records and chromosome-path unions
# ===========================================================================

test_that("combining inputs that carry no QC record yields no QC options", {
    # Empty, not a fabricated record: claiming QC that never ran would be
    # worse than reporting none.
    expect_equal(pecotmr:::.ssCombineQcOptions(list(), "combineX"), list())
})

test_that("chromosome paths union when they agree and abort when they clash", {
    # A sharded sketch routes a row to its file by CHR, so one chromosome
    # mapping to two files means the parts are not the same panel at all.
    mk <- function(paths) {
        new(
            "GenotypeHandle",
            path = "/tmp/x",
            format = "plink2",
            snpInfo = data.frame(
                SNP = "a",
                CHR = "chr1",
                BP = 1L,
                A1 = "A",
                A2 = "G",
                fileIdx = 1L,
                stringsAsFactors = FALSE
            ),
            nSamples = 1L,
            sampleIds = "s1",
            pgenPtr = NULL,
            chromPaths = paths
        )
    }
    h1 <- mk(c("1" = "/data/chr1.pgen"))
    h2 <- mk(c("1" = "/data/chr1.pgen", "2" = "/data/chr2.pgen"))
    out <- pecotmr:::.ssUnionChromPaths(list(h1, h2), "combineX")
    expect_equal(sort(names(out)), c("1", "2"))
    expect_equal(unname(out[["1"]]), "/data/chr1.pgen")
    # Same chromosome, different file -> refused.
    h3 <- mk(c("1" = "/other/chr1.pgen"))
    expect_error(
        pecotmr:::.ssUnionChromPaths(list(h1, h3), "combineX"),
        "chromosome '1' maps to two different genotype"
    )
})


test_that("panels read through no handle union their rows", {
    # A PLINK1, VCF or GDS panel has no handle whose snpInfo could be
    # widened, so the parts' assays are bound: each variant once, in genomic
    # order. A panel from another file is another source, and kept apart.
    skip_if_not_installed("PlinkMatrix")
    skip_if_not_installed("pgenlibr")
    stem <- test_path("test_data", "test_variants")
    panel <- readGenotypes(plink1Prefix = stem)
    a <- pecotmr:::.panelSubsetRows(panel, 4:8)
    b <- pecotmr:::.panelSubsetRows(panel, 1:5)
    merged <- pecotmr:::.ldRefMergeViews(
        pecotmr:::.ldRefDict(list(a, b)),
        "combineX"
    )
    expect_identical(merged$keys, c("ld1", "ld1"))
    out <- merged$refs$ld1
    expect_identical(rownames(out), rownames(panel)[1:8])
    expect_equal(
        as.matrix(SummarizedExperiment::assay(out, "dosage")),
        as.matrix(SummarizedExperiment::assay(panel, "dosage"))[1:8, ]
    )
    expect_identical(
        pecotmr:::.ldRefSource(out)$path,
        pecotmr:::.ldRefSource(panel)$path
    )
    # The same panel twice is that panel.
    expect_identical(pecotmr:::.ssUnionPanels(list(a, a), character(0)), a)
    apart <- pecotmr:::.ldRefMergeViews(
        pecotmr:::.ldRefDict(list(a, readGenotypes(plink2Prefix = stem))),
        "combineX"
    )
    expect_length(apart$refs, 2L)
})


# ===========================================================================
# Row-payload coercion and unmerged collection slots
# ===========================================================================

test_that("an empty collection coerces to itself, not to a row", {
    # fineMappingResult() hands over a single-row COLLECTION, so the
    # coercion has to accept that shape -- but an empty one has no row to
    # extract and must come back untouched.
    data(qtlFineMappingExample, twasWeightsExample)
    expect_s4_class(
        pecotmr:::.asFmRowPayload(qtlFineMappingExample[0]),
        "QtlFineMappingResult"
    )
    expect_s4_class(
        pecotmr:::.asTwRowPayload(twasWeightsExample[0]),
        "TwasWeights"
    )
    # A populated one becomes the row payload.
    expect_s4_class(
        pecotmr:::.asFmRowPayload(qtlFineMappingExample),
        "FineMappingRow"
    )
    expect_s4_class(
        pecotmr:::.asTwRowPayload(twasWeightsExample),
        "TwasWeightsRow"
    )
    # Anything that is not a collection passes straight through.
    expect_equal(pecotmr:::.asFmRowPayload(42), 42)
    expect_equal(pecotmr:::.asTwRowPayload("x"), "x")
})

test_that("a collection-level slot with no merge rule is refused", {
    # Silently dropping a slot when combining would lose provenance, so an
    # unrecognised one is an error naming the slot(s) and the class.
    setClass(
        "TsExtraKid",
        contains = "RangedTupleList",
        representation(
            genome = "character",
            ldReference = "ANY",
            surprise = "character"
        )
    )
    k <- new(
        "TsExtraKid",
        GenomicRanges::GRangesList(
            a = GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, width = 1L))
        ),
        genome = "hg38",
        ldReference = NULL,
        surprise = "boo"
    )
    expect_error(
        pecotmr:::.rtlExtraSlots(list(k, k), "combineX"),
        "no merge rule for the collection-level slot"
    )
    # The message names the offending slots and the class.
    expect_error(
        pecotmr:::.rtlExtraSlots(list(k, k), "combineX"),
        "surprise"
    )
})

test_that(".twrRowResolveWeights drops rows whose weights are all NA", {
    allNa <- TwasWeightsRow(
        variantIds = c("chr1:1:A:G", "chr1:2:C:T"),
        weights = c(NA_real_, NA_real_)
    )
    out <- pecotmr:::.twrRowResolveWeights(allNa)
    # Nothing usable survives, so the row contributes no variants at all.
    expect_length(out$variantIds, 0L)
    expect_length(out$weights, 0L)
    # Control: a single usable weight is kept, and only that one.
    partial <- TwasWeightsRow(
        variantIds = c("chr1:1:A:G", "chr1:2:C:T"),
        weights = c(NA_real_, 0.5)
    )
    kept <- pecotmr:::.twrRowResolveWeights(partial)
    expect_equal(kept$variantIds, "chr1:2:C:T")
    expect_equal(kept$weights, 0.5)
})

test_that(".ssConcatChr answers character(0) for no pieces", {
    expect_identical(pecotmr:::.ssConcatChr(list()), character(0))
})

# ---- the stored analysis window ---------------------------------------------

test_that("genomicRegion() is the analysis window, not the variant span", {
    data(qtlFineMappingExample)
    x <- qtlFineMappingExample
    # Built without a window: fall back to the derived span, so the accessor
    # never returns nothing where it used to return a range.
    expect_s4_class(genomicRegion(x), "GRanges")

    want <- GenomicRanges::GRanges("chr22:1-50000000")
    withWindow <- QtlFineMappingResult(
        studyName = as.character(.tupleColumn(x, "study")),
        context = as.character(.tupleColumn(x, "context")),
        trait = as.character(.tupleColumn(x, "trait")),
        method = as.character(.tupleColumn(x, "method")),
        entry = list(.fmrRowParts(x, 1L)),
        region = rep(want, nrow(x))
    )
    expect_identical(
        as.character(genomicRegion(withWindow)),
        rep(as.character(want), nrow(x))
    )
    # The window is WIDER than the variants that survived -- that gap is the
    # whole reason it has to be stored rather than derived.
    span <- unlist(range(withWindow), use.names = FALSE)
    expect_gt(
        GenomicRanges::start(span)[[1L]],
        GenomicRanges::start(want)[[1L]]
    )
    expect_lt(GenomicRanges::end(span)[[1L]], GenomicRanges::end(want)[[1L]])
})

test_that("the stored window survives subsetting, the derived span does not", {
    data(qtlFineMappingExample)
    x <- qtlFineMappingExample
    want <- GenomicRanges::GRanges("chr22:1-50000000")
    obj <- QtlFineMappingResult(
        studyName = as.character(.tupleColumn(x, "study")),
        context = as.character(.tupleColumn(x, "context")),
        trait = as.character(.tupleColumn(x, "trait")),
        method = as.character(.tupleColumn(x, "method")),
        entry = list(.fmrRowParts(x, 1L)),
        region = rep(want, nrow(x))
    )
    # Provenance: a subset of a fit's variants was still produced by a fit
    # over the original window, so this must NOT move.
    expect_identical(
        as.character(genomicRegion(obj[1L])),
        as.character(want)
    )
})

test_that("a region column must be a GRanges with one range per row", {
    data(qtlFineMappingExample)
    x <- qtlFineMappingExample
    build <- function(rg) {
        QtlFineMappingResult(
            studyName = as.character(.tupleColumn(x, "study")),
            context = as.character(.tupleColumn(x, "context")),
            trait = as.character(.tupleColumn(x, "trait")),
            method = as.character(.tupleColumn(x, "method")),
            entry = list(.fmrRowParts(x, 1L)),
            region = rg
        )
    }
    expect_error(build("chr22:1-2"), "must be a GRanges")
    expect_error(
        build(GenomicRanges::GRanges(c("chr22:1-2", "chr22:3-4"))),
        "same length as"
    )
})

test_that(".asAnalysisWindow keeps a window and spans a multi-range one", {
    expect_null(.asAnalysisWindow(NULL))
    expect_null(.asAnalysisWindow(GenomicRanges::GRanges()))
    one <- GenomicRanges::GRanges("chr1:100-200")
    expect_identical(.asAnalysisWindow(one), one)
    # A window is a GRanges; a string is refused rather than read.
    expect_error(.asAnalysisWindow("chr1:100-200"), "must be a GRanges")
    # A multi-range window collapses to the span it covered.
    span <- .asAnalysisWindow(GenomicRanges::GRanges(c(
        "chr1:100-200",
        "chr1:500-900"
    )))
    expect_length(span, 1L)
    expect_identical(GenomicRanges::start(span)[[1L]], 100L)
    expect_identical(GenomicRanges::end(span)[[1L]], 900L)
})


test_that("LdMatrix views are stored as the one covering them, in any order", {
    ids <- c("chr22:100:A:G", "chr22:200:C:T", "chr22:300:G:A")
    vars <- data.frame(
        chrom = "chr22",
        pos = c(100L, 200L, 300L),
        SNP = ids,
        A1 = c("G", "T", "A"),
        A2 = c("A", "C", "G")
    )
    R <- matrix(c(1, 0.5, 0.2, 0.5, 1, 0.3, 0.2, 0.3, 1), 3L, 3L)
    ld <- LdMatrix(R, variants = vars, nSamples = 100)
    # Two views neither of which covers the other are both kept: the LD
    # between their variants is held by neither.
    apart <- pecotmr:::.ldRefMergeViews(
        pecotmr:::.ldRefDict(list(ld[1:2], ld[2:3])),
        "combine"
    )
    expect_length(apart$refs, 2L)
    # With the matrix covering both among them, both are stored as it,
    # wherever it comes.
    orders <- list(
        list(ld[1:2], ld, ld[2:3]),
        list(ld[1:2], ld[2:3], ld),
        list(ld, ld[1:2], ld[2:3])
    )
    for (parts in orders) {
        merged <- pecotmr:::.ldRefMergeViews(
            pecotmr:::.ldRefDict(parts),
            "combine"
        )
        expect_identical(merged$refs, list(ld1 = ld))
        expect_identical(merged$keys, rep("ld1", 3L))
    }
})

test_that("binding stores an LdMatrix another part's covers as that one", {
    ids <- sprintf("chr1:%d:A:G", 100L * (1:6))
    r <- 0.3^abs(outer(1:6, 1:6, "-"))
    wide <- LdMatrix(`dimnames<-`(r, list(ids, ids)))
    narrow <- LdMatrix(`dimnames<-`(r[4:6, 4:6], list(ids[4:6], ids[4:6])))
    merged <- pecotmr:::.ldRefMergeViews(
        pecotmr:::.ldRefDict(list(wide, narrow, wide)),
        "probe"
    )
    expect_identical(merged$refs, list(ld1 = wide))
    expect_identical(merged$keys, c("ld1", "ld1", "ld1"))
    # Matrices that do not nest are different references, and stay apart.
    other <- LdMatrix(
        `dimnames<-`(0.9^abs(outer(1:6, 1:6, "-")), list(ids, ids))
    )
    apart <- pecotmr:::.ldRefMergeViews(
        pecotmr:::.ldRefDict(list(wide, other)),
        "probe"
    )
    expect_length(apart$refs, 2L)
    # A panel and a matrix are never views of one source.
    panel <- readGenotypes(test_path("test_data", "test_variants.bed"))
    kinds <- pecotmr:::.ldRefMergeViews(
        pecotmr:::.ldRefDict(list(wide, panel)),
        "probe"
    )
    expect_length(kinds$refs, 2L)
})

test_that("binding stores a mixture another part's covers as that one", {
    panel <- readGenotypes(test_path("test_data", "test_variants.bed"))
    halves <- list(a = panel[, 1:40], b = panel[, 41:100])
    mix <- GenotypeMixture(halves, mixtureWeights = c(0.5, 0.5))
    # Two views that overlap only in part are both held by the whole mixture.
    merged <- pecotmr:::.ldRefMergeViews(
        pecotmr:::.ldRefDict(list(mix[1:200], mix, mix[100:349])),
        "probe"
    )
    expect_identical(merged$refs, list(ld1 = mix))
    expect_identical(merged$keys, c("ld1", "ld1", "ld1"))
    # The same panels weighted differently are another mixture.
    reweighted <- GenotypeMixture(halves, mixtureWeights = c(0.4, 0.6))
    apart <- pecotmr:::.ldRefMergeViews(
        pecotmr:::.ldRefDict(list(mix, reweighted[1:200])),
        "probe"
    )
    expect_length(apart$refs, 2L)
})

# ===========================================================================
# GWAS block ranges and the QTL variants past them
# ===========================================================================

# Ten GWAS variants at chr1:100-1000, 100 apart.
.ts_gwasVariants <- function(pos = 100L * (1:10)) {
    gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(pos, width = 1L))
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = sprintf("chr1:%d:A:G", pos),
        A1 = "A",
        A2 = "G",
        Z = 0,
        N = 100L
    )
    gr
}

.ts_positions <- function(pos) {
    GenomicRanges::GRanges("chr1", IRanges::IRanges(pos, width = 1L))
}

test_that(".rowSpanRanges spans each element, an empty one a placeholder", {
    v <- .ts_gwasVariants()
    spans <- pecotmr:::.rowSpanRanges(
        GenomicRanges::GRangesList(v[1:3], v[0], v[5:10])
    )
    expect_identical(
        as.character(spans),
        c("chr1:100-300", "chrUn:1-0", "chr1:500-1000")
    )
    expect_identical(GenomicRanges::width(spans)[[2L]], 0L)
})

test_that(".asRegionRanges takes a GRanges, or reads a BED file as BED", {
    gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(10L, 20L))
    expect_identical(pecotmr:::.asRegionRanges(gr), gr)
    # A BED file's start is 0-based, and its fourth column names a window.
    path <- withr::local_tempfile(fileext = ".bed")
    writeLines("chr1\t9\t20\tw1", path)
    w <- pecotmr:::.asRegionRanges(path)
    expect_identical(as.character(w), "chr1:10-20")
    expect_identical(S4Vectors::mcols(w)$blockId, "w1")
    # Given in R, a window is a GRanges: a table or a range string is
    # refused rather than read, and a lone string is taken for a path.
    table <- data.frame(chrom = "chr1", start = 9L, end = 20L)
    expect_error(
        pecotmr:::.asRegionRanges(table),
        "must be a GRanges or the path of a BED file"
    )
    expect_error(
        pecotmr:::.asRegionRanges(c("chr1:10-20", "chr2:5-8")),
        "must be a GRanges or the path of a BED file"
    )
    expect_error(
        pecotmr:::.asRegionRanges("chr1:10-20"),
        "no file exists at 'chr1:10-20'"
    )
    expect_error(pecotmr:::.asRegionRanges(42), "must be a GRanges")
})

test_that(".assertRegionGRanges takes a GRanges and nothing else", {
    gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(10L, 20L))
    expect_identical(pecotmr:::.assertRegionGRanges(gr), gr)
    expect_error(
        pecotmr:::.assertRegionGRanges("chr1:10-20"),
        "`region` must be a GRanges"
    )
})

test_that(".rowRegionAt reads one row's window, or its span", {
    v <- .ts_gwasVariants()
    x <- GwasSumStats(
        studyName = c("s", "s"),
        entry = list(v[1:5], v[0L]),
        genome = "hg19",
        blockId = c("a", "b")
    )
    expect_identical(
        as.character(pecotmr:::.rowRegionAt(x, 1L)),
        "chr1:100-500"
    )
    # An empty row built empty has only a placeholder, so no window.
    expect_null(pecotmr:::.rowRegionAt(x, 2L))
    old <- x
    old$region <- NULL
    expect_identical(
        as.character(pecotmr:::.rowRegionAt(old, 1L)),
        "chr1:100-500"
    )
})

test_that(".regionColumns keeps a given range and spans the rest", {
    given <- GenomicRanges::GRanges(
        c("chr1", "chrUn"),
        IRanges::IRanges(1L, c(5000L, 0L))
    )
    spans <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100L, 300L))
    cols <- pecotmr:::.regionColumns(given, c(TRUE, FALSE), rep(spans, 2L))
    expect_identical(cols$regionGiven, c(TRUE, FALSE))
    expect_identical(
        as.character(cols$region),
        c("chr1:1-5000", "chr1:100-300")
    )
    # A row with no variants keeps the range it was given.
    empty <- pecotmr:::.regionColumns(
        given[1L],
        TRUE,
        pecotmr:::.placeholderRanges(1L)
    )
    expect_true(empty$regionGiven)
    expect_identical(as.character(empty$region), "chr1:1-5000")
})

test_that(".refuseUncontainedRegions wants a given range to hold its row", {
    pieces <- GenomicRanges::GRangesList(
        GenomicRanges::GRanges(
            "chr1",
            IRanges::IRanges(c(100L, 300L), width = 1L)
        ),
        GenomicRanges::GRanges(),
        GenomicRanges::GRanges(
            "chr1",
            IRanges::IRanges(c(100L, 300L), width = 1L)
        )
    )
    given <- function(chrom, start, end) {
        GenomicRanges::GRanges(chrom, IRanges::IRanges(start, end))
    }
    refuse <- function(g, has = c(TRUE, TRUE, FALSE)) {
        pecotmr:::.refuseUncontainedRegions(g, has, pieces, "f")
    }
    # Held; an empty row holds nothing, and a row given no range is not
    # checked.
    expect_null(refuse(given(c("chr1", "chr9", "chr9"), 50L, 400L)))
    expect_error(
        refuse(given(c("chr1", "chr1", "chr1"), 150L, 400L)),
        "row 1: range chr1:150-400, variants at chr1:100-300",
        fixed = TRUE
    )
    # Another chromosome holds none of them.
    expect_error(
        refuse(given(c("chr2", "chr1", "chr1"), 1L, 5000L)),
        "variants at chr1:100-300"
    )
    # A variant sits at its start, so an indel's width cannot carry it past
    # the end of the block it starts in.
    indel <- GenomicRanges::GRangesList(
        GenomicRanges::GRanges("chr1", IRanges::IRanges(300L, 310L))
    )
    expect_null(pecotmr:::.refuseUncontainedRegions(
        given("chr1", 1L, 305L),
        TRUE,
        indel,
        "f"
    ))
})

test_that(".rowRegions falls back to the span where no range is stored", {
    x <- GwasSumStats(
        studyName = "s",
        entry = list(.ts_gwasVariants()),
        genome = "hg19",
        region = GenomicRanges::GRanges("chr1:1-5000")
    )
    r <- pecotmr:::.rowRegions(x)
    expect_identical(as.character(r$region), "chr1:1-5000")
    expect_true(r$regionGiven)
    # Saved before the columns existed: the span, computed.
    old <- x
    old$region <- NULL
    old$regionGiven <- NULL
    r <- pecotmr:::.rowRegions(old)
    expect_identical(as.character(r$region), "chr1:100-1000")
    expect_false(r$regionGiven)
    # A range without its origin is not taken as given.
    noOrigin <- x
    noOrigin$regionGiven <- NULL
    expect_false(pecotmr:::.rowRegions(noOrigin)$regionGiven)
    # A placeholder a combine padded in falls back to the span as well.
    padded <- x
    padded$region <- pecotmr:::.placeholderRanges(1L)
    padded$regionGiven <- NA
    r <- pecotmr:::.rowRegions(padded)
    expect_identical(as.character(r$region), "chr1:100-1000")
    expect_false(r$regionGiven)
})

test_that(".checkGwasBlockReach words what it finds", {
    v <- .ts_gwasVariants()
    computed <- GwasSumStats(
        studyName = c("s", "s"),
        entry = list(v[1:5], v[6:10]),
        genome = "hg19",
        blockId = c("bA", "bB")
    )
    # Two variants before the first block's span count against it.
    expect_warning(
        pecotmr:::.checkGwasBlockReach(
            list(u = .ts_positions(c(10L, 20L, 100L))),
            computed,
            "f",
            "units"
        ),
        "u: 2 variant(s) past the start of block bA (chr1:100-500) of GWAS",
        fixed = TRUE
    )
    # Past five lines the rest are counted.
    many <- set_names(
        rep(list(.ts_positions(c(1000L, 1100L, 1200L))), 7L),
        str_c("u", 1:7)
    )
    msg <- tryCatch(
        pecotmr:::.checkGwasBlockReach(many, computed, "f", "units"),
        warning = conditionMessage
    )
    expect_match(msg, "u5: 2 variant(s)", fixed = TRUE)
    expect_match(msg, "... and 2 more.", fixed = TRUE)
    expect_no_match(msg, "u6: ", fixed = TRUE)
    # A fine-mapping result's fit is a study and a method, and its given
    # window is checked strictly.
    fmr <- GwasFineMappingResult(
        studyName = "s",
        method = "susie",
        entry = list(.sc_makeFineMappingRow(3)),
        region = GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, 350L))
    )
    expect_error(
        pecotmr:::.checkGwasBlockReach(
            list(u = .ts_positions(c(200L, 400L))),
            fmr,
            "f",
            "units"
        ),
        "u: 1 variant(s) past the end of block chr1_100_300 (chr1:1-350) of",
        fixed = TRUE
    )
    # Nothing to check against, or with.
    expect_null(
        pecotmr:::.checkGwasBlockReach(list(), computed, "f", "units")
    )
    expect_null(
        pecotmr:::.checkGwasBlockReach(
            list(u = .ts_positions(100L)),
            computed[0L],
            "f",
            "units"
        )
    )
})

# ===========================================================================
# Fit coverage: what a narrowed row's stored fit still describes
# ===========================================================================

.ts_fitted <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")

# A fine-mapping row whose fit was made on `.ts_fitted`, holding `keep` of
# them.
.ts_fitRow <- function(keep = seq_along(.ts_fitted)) {
    pip <- c(0.5, 0.3, 0.2)
    full <- FineMappingRow(
        .ts_fitted,
        list(pip = set_names(pip, .ts_fitted)),
        data.frame(variant_id = .ts_fitted, pip = pip)
    )
    methods::new(
        "FineMappingRow",
        variants = variants(full)[keep],
        susieFit = susieFit(full)
    )
}

# A collection of `n` rows, each fitted on `.ts_fitted` and narrowed to its
# first two.
.ts_narrowedCollection <- function(n) {
    pip <- c(0.5, 0.3, 0.2)
    row <- FineMappingRow(
        .ts_fitted,
        list(pip = pip),
        data.frame(variant_id = .ts_fitted, pip = pip)
    )
    x <- QtlFineMappingResult(
        studyName = sprintf("s%d", seq_len(n)),
        context = rep("c1", n),
        trait = rep("g1", n),
        method = rep("susie", n),
        entry = rep(list(row), n)
    )
    endoapply(x, utils::head, 2L)
}

test_that(".lostFitVariantCount counts the fit's variants a row lost", {
    f <- pecotmr:::.lostFitVariantCount
    expect_equal(f(.ts_fitted, .ts_fitted), 0L)
    expect_equal(f(.ts_fitted, .ts_fitted[1:2]), 1L)
    expect_equal(f(character(0), .ts_fitted), 0L)
    # A variant held under its swapped alleles is still held.
    expect_equal(f(.ts_fitted, c("chr1:100:G:A", .ts_fitted[2:3])), 0L)
})

test_that(".fitColumnOf finds each held variant's column in the fit", {
    f <- pecotmr:::.fitColumnOf
    expect_equal(f(.ts_fitted[c(3L, 1L)], .ts_fitted), c(3L, 1L))
    # Allele-aware for a swapped id; NA for a variant the fit never had.
    expect_equal(f(c("chr1:200:T:C", "chr1:900:A:G"), .ts_fitted), c(2L, NA))
})

test_that(".mdFitVariants reads a group's first row, nothing for old rows", {
    md <- S4Vectors::DataFrame(
        fitVariants = IRanges::CharacterList(list(c("a", "b"), c("a", "b")))
    )
    # Rows split from one entry repeat its fit and its variants alike.
    expect_equal(pecotmr:::.mdFitVariants(md, c(1L, 2L)), c("a", "b"))
    expect_equal(
        pecotmr:::.mdFitVariants(S4Vectors::DataFrame(x = 1:2), 1L),
        character(0)
    )
})

test_that(".partsFitIds: the row's own ids, else the fit's names", {
    # A row holding all of its fit is its fit's variables, read by position.
    full <- .ts_fitRow()
    expect_equal(pecotmr:::.partsFitIds(full), .ts_fitted)
    # A narrowed row's fit names them, wherever it keeps the names.
    expect_equal(pecotmr:::.partsFitIds(.ts_fitRow(1:2)), .ts_fitted)
    named <- matrix(0, 1L, 3L, dimnames = list(NULL, rev(.ts_fitted)))
    fitless <- methods::new(
        "FineMappingRow",
        variants = variants(full)[1:2],
        susieFit = NULL
    )
    expect_equal(pecotmr:::.partsFitIds(fitless, named), rev(.ts_fitted))
    # A narrowed row whose fit names none cannot be placed in it.
    unnamed <- methods::new(
        "FineMappingRow",
        variants = variants(full)[1:2],
        susieFit = list(pip = c(0.5, 0.3, 0.2))
    )
    expect_error(
        pecotmr:::.partsFitIds(unnamed),
        "holds 2 variants but its fit was made with 3"
    )
})

test_that(".fitColumnsFor reads a row holding its whole fit by position", {
    f <- pecotmr:::.fitColumnsFor
    gr <- variants(.ts_fitRow())
    fit <- list(pip = c(0.5, 0.3, 0.2))
    expect_equal(f(gr, fit), 1:3)
    # Narrowed: by the fit's names, or not at all without them.
    expect_equal(f(gr[2:3], list(pip = set_names(fit$pip, .ts_fitted))), 2:3)
    expect_null(f(gr[2:3], fit))
    # A fit of unknown width describes the row as it is.
    expect_equal(f(gr[2:3], NULL), 1:2)
})

test_that(".fitWidth counts the fit's variables, never the null column", {
    f <- pecotmr:::.fitWidth
    alpha <- matrix(0.25, 1L, 4L)
    expect_equal(f(list(pip = c(0.1, 0.2, 0.3), alpha = alpha)), 3L)
    expect_equal(f(list(alpha = alpha)), 4L)
    expect_equal(f(list(), alpha), 4L)
    expect_true(is.na(f(NULL)))
    expect_true(is.na(f(list(list(pip = 1)))))
})

test_that("a TWAS row records fit variants only when it carries a fit", {
    w <- c(0.1, 0.2, 0.3)
    withCv <- TwasWeightsRow(.ts_fitted, w, cvResult = list(r2 = 0.5))
    expect_equal(pecotmr:::.twRowFitVariants(withCv), .ts_fitted)
    # Weights alone are no fit payload.
    expect_equal(
        pecotmr:::.twRowFitVariants(TwasWeightsRow(.ts_fitted, w)),
        character(0)
    )
})

test_that(".fitCoverageChecked warns about a narrowed row and returns it", {
    narrowed <- .ts_fitRow(1:2)
    expect_warning(
        out <- pecotmr:::.fitCoverageChecked(narrowed, "f()"),
        "f\\(\\): the stored fit was made with 3 variants, 1 of which"
    )
    expect_identical(out, narrowed)
    cnd <- tryCatch(
        pecotmr:::.fitCoverageChecked(narrowed, "f()"),
        warning = identity
    )
    expect_s3_class(cnd, "pecotmrFitCoverage")
    expect_no_warning(pecotmr:::.fitCoverageChecked(.ts_fitRow(), "f()"))
})

test_that(".fitCoverageCheckedRows names at most five narrowed rows", {
    x <- .ts_narrowedCollection(6L)
    cnd <- tryCatch(
        pecotmr:::.fitCoverageCheckedRows(x, "f()"),
        warning = identity
    )
    expect_s3_class(cnd, "pecotmrFitCoverage")
    msg <- conditionMessage(cnd)
    expect_match(msg, "f\\(\\): 6 row\\(s\\) hold fewer variants")
    expect_equal(stringr::str_count(msg, "\\|c1\\|g1\\|susie"), 5L)
    expect_invisible(suppressWarnings(
        pecotmr:::.fitCoverageCheckedRows(x, "f()")
    ))
})

test_that(".fitCoverageCheckedRows is silent without a fit or rows", {
    x <- .ts_narrowedCollection(1L)
    expect_no_warning(pecotmr:::.fitCoverageCheckedRows(x[0], "f()"))
    fitless <- x
    fitless$susieFit <- S4Vectors::SimpleList(list(NULL))
    expect_no_warning(pecotmr:::.fitCoverageCheckedRows(fitless, "f()"))
    # An optional input left out is nothing to check.
    expect_no_warning(pecotmr:::.fitCoverageCheckedRows(NULL, "f()"))
})

# TWAS weights fitted on `.ts_fitted`, holding `keep` of them.
.ts_weights <- function(keep = seq_along(.ts_fitted)) {
    tw <- TwasWeights(
        studyName = "w1",
        context = "c1",
        trait = "g1",
        method = "lasso",
        entry = list(
            TwasWeightsRow(.ts_fitted, c(0.1, 0.2, 0.3), cvResult = list())
        )
    )
    endoapply(tw, `[`, keep)
}

test_that(".narrowedRowLabels labels the narrowed rows and no others", {
    labels <- pecotmr:::.narrowedRowLabels(.ts_narrowedCollection(2L))
    expect_length(labels, 2L)
    expect_match(labels, "^s[12]\\|c1\\|g1\\|susie")
    expect_identical(pecotmr:::.narrowedRowLabels(.ts_weights()), character(0))
    expect_identical(
        pecotmr:::.narrowedRowLabels(.ts_weights(1:2)),
        pecotmr:::.pubRowLabels(.ts_weights(1:2), 1L)
    )
})

test_that(".fitCoverageCheckedAll names every input's narrowed rows at once", {
    fm <- .ts_narrowedCollection(2L)
    cnd <- tryCatch(
        pecotmr:::.fitCoverageCheckedAll(
            list(NULL, fm, .ts_weights(1:2)),
            "f()"
        ),
        warning = identity
    )
    expect_s3_class(cnd, "pecotmrFitCoverage")
    expect_match(conditionMessage(cnd), "f\\(\\): 3 row\\(s\\) hold fewer")
    expect_match(conditionMessage(cnd), "w1\\|c1\\|g1\\|lasso")
    expect_null(suppressWarnings(
        pecotmr:::.fitCoverageCheckedAll(list(fm), "f()")
    ))
    expect_no_warning(
        pecotmr:::.fitCoverageCheckedAll(list(NULL, .ts_weights()), "f()")
    )
    # Anything but a collection is left to its caller's own validation.
    expect_identical(
        pecotmr:::.narrowedRowLabels(data.frame(a = 1)),
        character(0)
    )
})

test_that("the joint warning's hint is that of the first narrowed input", {
    fm <- .ts_narrowedCollection(1L)
    weightsFirst <- tryCatch(
        pecotmr:::.fitCoverageCheckedAll(list(.ts_weights(1:2), fm), "f()"),
        warning = identity
    )
    expect_match(conditionMessage(weightsFirst), "Refit on the variants kept")
    expect_false(str_detect(
        conditionMessage(weightsFirst),
        "intersectVariants"
    ))
    # An input with nothing narrowed does not choose the hint.
    fmFirst <- tryCatch(
        pecotmr:::.fitCoverageCheckedAll(list(.ts_weights(), fm), "f()"),
        warning = identity
    )
    expect_match(conditionMessage(fmFirst), "intersectVariants")
})

test_that(".muffleFitCoverage muffles fit-coverage warnings and no others", {
    expect_no_warning(
        out <- pecotmr:::.muffleFitCoverage({
            warn("narrowed", class = "pecotmrFitCoverage")
            "value"
        })
    )
    expect_equal(out, "value")
    expect_warning(pecotmr:::.muffleFitCoverage(warn("other")), "other")
})

test_that("the coverage hint names intersectVariants() only where it applies", {
    expect_match(
        pecotmr:::.fitCoverageHint(.ts_fitRow()),
        "intersectVariants\\(\\) restricts a fit"
    )
    expect_match(
        pecotmr:::.fitCoverageHint(.ts_narrowedCollection(1L)),
        "intersectVariants"
    )
    weightsHint <- pecotmr:::.fitCoverageHint(
        TwasWeightsRow(.ts_fitted, c(0.1, 0.2, 0.3))
    )
    expect_false(str_detect(weightsHint, "intersectVariants"))
    expect_match(weightsHint, "Refit on the variants kept")
})

test_that(".naLikeColumn pads a list column with empty entries of its kind", {
    cl <- pecotmr:::.naLikeColumn(IRanges::CharacterList("a"), 2L)
    expect_s4_class(cl, "CharacterList")
    expect_equal(unname(lengths(cl)), c(0L, 0L))
    sl <- pecotmr:::.naLikeColumn(S4Vectors::SimpleList(list(1)), 2L)
    expect_s4_class(sl, "SimpleList")
    expect_length(sl, 2L)
    expect_null(sl[[1L]])
})

test_that("binding a collection saved before derivedColumns pads the column", {
    x <- .fmRebuildExample(.fmLoadExample("qtlFineMappingExample"))
    recorded <- which(lengths(x$derivedColumns) > 0L)[[1L]]
    old <- x[-recorded]
    old$derivedColumns <- NULL
    out <- pecotmr:::.rbindCollections(list(old, x[recorded]))
    expect_s4_class(out$derivedColumns, "CharacterList")
    expect_equal(
        unname(lengths(out$derivedColumns)),
        c(rep(0L, nrow(old)), length(x$derivedColumns[[recorded]]))
    )
})

# ===========================================================================
# .fmrRowTopLoci
# ===========================================================================

test_that(".fmrRowTopLoci() tables convert to ranges", {
    tl <- data.frame(
        variant_id = c("1:100:A:G", "1:200:C:T"),
        pip = c(0.9, 0.1),
        betahat = c(0.5, -0.2),
        sebetahat = c(0.1, 0.2),
        cs = c(1L, 0L),
        method = "susie",
        stringsAsFactors = FALSE
    )
    ent <- FineMappingRow(
        variantIds = tl$variant_id,
        susieFit = list(),
        topLoci = tl
    )
    # .fmrRowTopLoci() is the per-row table; topLoci() converts once, after
    # the aggregate view has stacked the rows, so the conversion is tested
    # through the converter it now lives in.
    out <- .fmrRowTopLoci(ent)
    expect_s3_class(out, "data.frame")
    gr <- pecotmr:::.variantTableGRanges(out)
    expect_s4_class(gr, "GRanges")
    expect_equal(length(gr), 2)
    expect_equal(S4Vectors::mcols(gr)$pip, c(0.9, 0.1))
})

test_that(".fmrRowTopLoci() empty input converts to empty ranges", {
    ent <- FineMappingRow(
        variantIds = character(0),
        susieFit = list(),
        topLoci = data.frame()
    )
    gr <- pecotmr:::.variantTableGRanges(.fmrRowTopLoci(ent))
    expect_s4_class(gr, "GRanges")
    expect_equal(length(gr), 0)
})

test_that(".fmrRowTopLoci(raw = TRUE) returns the stored table verbatim", {
    vids <- paste0("chr1:", 1:4, ":A:G")
    entry <- .makeAdjustEntry(vids)
    raw <- .fmrRowTopLoci(entry, raw = TRUE)
    expect_equal(nrow(raw), length(vids))
    expect_true(all(c("betahat", "sebetahat") %in% names(raw)))
    # The default view filters by pip and projects; raw does neither.
    expect_false(identical(names(raw), names(.fmrRowTopLoci(entry))))
})

test_that("topLoci surfaces the mvsusie lfsr through the projector", {
    tl <- data.frame(
        variant_id = c("chr1:100:A:G", "chr1:200:C:T"),
        chrom = "1",
        pos = c(100L, 200L),
        A1 = c("A", "C"),
        A2 = c("G", "T"),
        pip = c(0.6, 0.4),
        posterior_mean = c(0.1, 0.2),
        posterior_sd = c(0.05, 0.05),
        cs_95 = c("mvsusie_1", "mvsusie_1"),
        lfsr = c(0.001, 0.02),
        stringsAsFactors = FALSE
    )
    e <- FineMappingRow(
        variantIds = tl$variant_id,
        susieFit = list(pip = tl$pip),
        topLoci = tl
    )
    out <- .fmrRowTopLoci(e, signalCutoff = 0)
    expect_true(is_in("lfsr", names(out)))
    expect_equal(out$lfsr, c(0.001, 0.02))
})
