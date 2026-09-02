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
        S4Vectors::mcols(grl) <- do.call(
            S4Vectors::DataFrame,
            c(cols, list(check.names = FALSE))
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
            study = "s1",
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
            study = c("s1", "s2"),
            context = "c1",
            trait = "t1",
            method = "susie"
        ),
        "must each be length 1"
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
            study = "s1",
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
            study = "ghost",
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
        pecotmr:::.tupleSelectRowGwasFmr(empty, study = "g1", method = "susie"),
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
        "Pass `study` and `method`"
    )
})

test_that(".tupleSelectRowGwasFmr: non-scalar region errors", {
    multi <- .ts_coll(
        study = c("g1", "g1"),
        method = c("susie", "susie"),
        blockId = c("r1", "r2"),
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::.tupleSelectRowGwasFmr(
            multi,
            study = "g1",
            method = "susie",
            region = c("r1", "r2")
        ),
        "`region` must be length 1"
    )
})

test_that(".tupleSelectRowGwasFmr: region disambiguates per-block rows", {
    # Same (study, method) across two regions; region picks the right row.
    multi <- .ts_coll(
        study = c("g1", "g1"),
        method = c("susie", "susie"),
        blockId = c("chr22_1_100", "chr22_500_600"),
        stringsAsFactors = FALSE
    )
    expect_equal(
        pecotmr:::.tupleSelectRowGwasFmr(
            multi,
            study = "g1",
            method = "susie",
            region = "chr22_500_600"
        ),
        2L
    )
})

test_that(".tupleSelectRowGwasFmr: missing tuple errors and includes region in message", {
    one <- .ts_coll(
        study = "g1",
        method = "susie",
        blockId = "r1",
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::.tupleSelectRowGwasFmr(
            one,
            study = "g1",
            method = "susie",
            region = "ghost"
        ),
        "region='ghost'"
    )
})

test_that(".tupleSelectRowGwasFmr: ambiguous multi-match (no region) lists candidates", {
    # Two rows share (study, method); .tupleSelectRowGwasFmr should error
    # listing the available blockIds since the caller didn't disambiguate.
    multi <- .ts_coll(
        study = c("g1", "g1"),
        method = c("susie", "susie"),
        blockId = c("region_A", "region_B"),
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::.tupleSelectRowGwasFmr(multi, study = "g1", method = "susie"),
        "pass `region` to disambiguate"
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
    expect_equal(pecotmr:::.fmrRowsMatching(df, study = "s1"), c(1L, 2L))
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

test_that(".fmrRowsMatching: `region` matches the blockId column", {
    gwas <- .ts_coll(
        study = c("g1", "g1"),
        method = c("susie", "susie"),
        blockId = c("r1", "r2"),
        stringsAsFactors = FALSE
    )
    expect_equal(pecotmr:::.fmrRowsMatching(gwas, region = "r2"), 2L)
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
        pecotmr:::.fmrRowsMatching(df, study = c("s1", "s3")),
        c(1L, 3L)
    )
})

# ===========================================================================
# .fmrRowMetadata (stable 5-column identity frame)
# ===========================================================================

test_that(".fmrRowMetadata: emits all five identity columns, NA-filling absent ones", {
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
        c("study", "context", "trait", "blockId", "method")
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

test_that(".fmrRowMetadata: zero-row input yields a zero-row 5-column frame", {
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
        c("study", "context", "trait", "blockId", "method")
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
    # Replaces the old .appendRegionCol / .validateRegionColumn pair. The
    # stored `region` GRanges column retired because it had no correct update
    # rule under subsetRegion() -- it would simply go stale. blockId is a
    # label keying the external block manifest, so subsetting cannot
    # invalidate it.
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

test_that(".getRegionColumn derives the span rather than reading a column", {
    # The one range concept is the intrinsic element range; there is no stored
    # window to fall back to any more.
    data(gwasFineMappingExample, envir = environment())
    gr <- pecotmr:::.getRegionColumn(gwasFineMappingExample)
    expect_s4_class(gr, "GRanges")
    expect_length(gr, nrow(gwasFineMappingExample))
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
        "'traitPos' column must be a GRanges"
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
