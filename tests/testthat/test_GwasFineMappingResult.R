# Tests for R/GwasFineMappingResult.R

# === Tests migrated from test_s4Constructors.R (GwasFineMappingResult) ===

test_that("GwasFineMappingResult: builds a collection keyed by 2-tuple", {
    e1 <- .sc_makeFineMappingRow(3)
    e2 <- .sc_makeFineMappingRow(3)
    res <- GwasFineMappingResult(
        studyName = c("g1", "g2"),
        method = c("susie", "susie"),
        entry = list(e1, e2)
    )
    expect_s4_class(res, "GwasFineMappingResult")
    expect_equal(nrow(res), 2L)
})

test_that("GwasFineMappingResult: region is derived from the variants", {
    # genomicRegion() reports the element's REALIZED variant span, not a nominal
    # window parsed out of blockId. A stored window had no correct update
    # rule under subsetRegion() and would quietly go stale; the span is in
    # sync by construction.
    e1 <- .sc_makeFineMappingRow(3)
    e2 <- .sc_makeFineMappingRow(2)
    res <- GwasFineMappingResult(
        studyName = c("g1", "g1"),
        method = c("susie", "susie"),
        blockId = c("chr1_1_500", "chr2_600_900"),
        entry = list(e1, e2)
    )
    reg <- genomicRegion(res)
    expect_equal(as.character(GenomicRanges::seqnames(reg)), c("chr1", "chr1"))
    expect_equal(GenomicRanges::start(reg), c(100L, 100L))
    expect_equal(GenomicRanges::end(reg), c(300L, 200L))
    # A synthetic blockId no longer produces a chrUn sentinel: with the span
    # derived from the variants there is nothing to fabricate.
    res2 <- GwasFineMappingResult(
        studyName = "g1",
        method = "susie",
        entry = list(e1)
    )
    expect_equal(
        as.character(GenomicRanges::seqnames(genomicRegion(res2))),
        "chr1"
    )
})


test_that("GwasFineMappingResult: validity does not recurse on key subset (#546)", {
    # The validity method builds a key-column data.frame to check tuple
    # uniqueness. Doing so via `object[, keyCols]` preserves the
    # GwasFineMappingResult class while dropping the required `entry` column;
    # older S4Vectors revalidates that intermediate and fails with
    # "missing columns: entry".
    e <- .sc_makeFineMappingRow(3)
    res <- GwasFineMappingResult(
        studyName = "g1",
        method = "susie",
        entry = list(e)
    )
    expect_s4_class(res, "GwasFineMappingResult")
    expect_true(validObject(res))

    # blockId is optional provenance and none was supplied, so the identity
    # columns are just (study, method); the range supplies the rest.
    sub <- S4Vectors::mcols(res)[, c("study", "method")]
    expect_false("entry" %in% names(sub))
})


test_that("GwasFineMappingResult: errors on length mismatch", {
    e <- .sc_makeFineMappingRow(3)
    expect_error(
        GwasFineMappingResult(
            studyName = c("g1", "g2"),
            method = c("susie"),
            entry = list(e)
        ),
        "same length"
    )
})


test_that("GwasFineMappingResult: same (study, method), other ranges", {
    # Row identity is (study, method, range), so (study, method) may repeat as
    # long as the two rows cover different genomic spans. This is the
    # genome-wide-across-blocks shape that qtlEnrichmentPipeline +
    # colocPipeline expect, and it now needs no label to disambiguate it.
    e1 <- .sc_makeFineMappingRow(3)
    e2 <- .sc_makeFineMappingRow(3, offset = 10000L)
    res <- GwasFineMappingResult(
        studyName = c("g1", "g1"),
        method = c("susie", "susie"),
        entry = list(e1, e2)
    )
    expect_s4_class(res, "GwasFineMappingResult")
    expect_equal(nrow(res), 2L)
    # No blockId was supplied, so the column is simply absent -- nothing
    # synthesises a placeholder id any more.
    expect_false(is_in("blockId", colnames(res)))
    expect_equal(length(unique(.rtlRangeKeys(res))), 2L)
})


test_that("GwasFineMappingResult: rejects two rows covering the same range", {
    # Same (study, method) AND the same span: these are the same block, and no
    # blockId label can make them distinct -- which is the point of keying on
    # the range rather than on a label.
    e1 <- .sc_makeFineMappingRow(3)
    e2 <- .sc_makeFineMappingRow(3)
    expect_error(
        GwasFineMappingResult(
            studyName = c("g1", "g1"),
            method = c("susie", "susie"),
            blockId = c("chr22_1_100", "chr22_200_300"),
            entry = list(e1, e2)
        ),
        "uniqueness violated"
    )
})


test_that("GwasFineMappingResult: errors when blockId length mismatches", {
    e1 <- .sc_makeFineMappingRow(3)
    e2 <- .sc_makeFineMappingRow(3)
    expect_error(
        GwasFineMappingResult(
            studyName = c("g1", "g1"),
            method = c("susie", "susie"),
            blockId = "only_one",
            entry = list(e1, e2)
        ),
        "same length"
    )
})


test_that("GwasFineMappingResult: show prints summary", {
    e <- .sc_makeFineMappingRow(3)
    res <- GwasFineMappingResult(
        studyName = "g1",
        method = "susie",
        entry = list(e)
    )
    expect_output(show(res), "GwasFineMappingResult")
})

# ===========================================================================
# TwasWeights collection
# ===========================================================================

# === Tests migrated from test_showMethods.R (GwasFineMappingResult) ===

test_that("show.GwasFineMappingResult prints (study, method) summary", {
    res <- GwasFineMappingResult(
        studyName = c("g1", "g1"),
        method = c("susie", "susieRss"),
        entry = list(.sh_makeFmEntry(), .sh_makeFmEntry())
    )
    out <- capture.output(show(res))
    expect_true(any(grepl("GwasFineMappingResult: 2 entries", out)))
    expect_true(any(grepl("1 studies.*2 methods", out)))
    expect_true(any(grepl("LD sketch: NULL", out)))
})


test_that("show.GwasFineMappingResult reports the ldSketch source when present", {
    res <- GwasFineMappingResult(
        studyName = "g1",
        method = "susie",
        entry = list(.sh_makeFmEntry()),
        ldSketch = .sh_makeGenotypeHandle()
    )
    out <- capture.output(show(res))
    expect_true(any(grepl("LD sketch: gds @ /tmp/test.gds", out)))
})


# === Tests migrated from test_collectionAccessors.R (GwasFineMappingResult) ===

test_that("GwasFineMappingResult: pip with study/method selectors", {
    e1 <- .ca_makeFmEntry(3)
    e2 <- .ca_makeFmEntry(4)
    res <- GwasFineMappingResult(
        studyName = c("g1", "g2"),
        method = c("susie", "susie"),
        entry = list(e1, e2)
    )
    pip <- pip(res, studyName = "g2", method = "susie")
    expect_equal(length(pip), 4L)
})


test_that("GwasFineMappingResult: contexts/traitNames return NULL", {
    e <- .ca_makeFmEntry(3)
    res <- GwasFineMappingResult(
        studyName = "g1",
        method = "susie",
        entry = list(e)
    )
    expect_null(contexts(res))
    expect_null(traitNames(res))
})


test_that("GwasFineMappingResult: credibleSets/topLoci/susieFit/variantIds dispatch", {
    e <- .ca_makeFmEntry(3)
    res <- GwasFineMappingResult(
        studyName = "g1",
        method = "susie",
        entry = list(e)
    )
    expect_equal(nrow(credibleSets(res)), 2L)
    # topLoci returns the projected posterior view (filtered by default
    # signalCutoff = 0.025; .ca_makeTopLoci sets all pip > 0.025 so all rows
    # survive). Compare on the projected shape, not the slot's raw shape.
    tl <- topLoci(res, signalCutoff = 0)
    expect_equal(length(tl), 3L)
    expect_equal(tl$variant_id, .ca_makeTopLoci(3)$variant_id)
    expect_equal(susieFit(res), list(payload = "fit_n=3"))
    expect_equal(length(variantIds(res)), 3L)
})

test_that("GwasFineMappingResult: topLoci aggregates per-block rows genome-wide", {
    # A genome-wide collection: same (study, method) across two region blocks.
    # With no selectors topLoci now stacks both blocks, tagging each variant
    # with its blockId; context/trait are NA-filled (GWAS keys on region).
    e1 <- .sc_makeFineMappingRow(3)
    e2 <- .sc_makeFineMappingRow(2)
    res <- GwasFineMappingResult(
        studyName = c("g1", "g1"),
        method = c("susie", "susie"),
        blockId = c("chr1:1-100", "chr1:200-300"),
        entry = list(e1, e2)
    )
    agg <- topLoci(res, signalCutoff = 0)
    expect_equal(length(agg), 5L)
    expect_equal(
        agg$blockId,
        c(
            "chr1:1-100",
            "chr1:1-100",
            "chr1:1-100",
            "chr1:200-300",
            "chr1:200-300"
        )
    )
    expect_true(all(is.na(agg$context)))
    expect_true(all(is.na(agg$trait)))
    expect_true("variant_id" %in% names(S4Vectors::mcols(agg)))
})

test_that("GwasFineMappingResult: topLoci region= selects a single block", {
    e1 <- .sc_makeFineMappingRow(3)
    e2 <- .sc_makeFineMappingRow(2)
    res <- GwasFineMappingResult(
        studyName = c("g1", "g1"),
        method = c("susie", "susie"),
        blockId = c("chr1:1-100", "chr1:200-300"),
        entry = list(e1, e2)
    )
    # region= pins one row, so this hits the single-entry fast path (bare table).
    tl <- topLoci(
        res,
        studyName = "g1",
        method = "susie",
        region = "chr1:200-300",
        signalCutoff = 0
    )
    expect_equal(length(tl), 2L)
    expect_false("blockId" %in% names(S4Vectors::mcols(tl)))
})

test_that("GwasFineMappingResult: credibleSets aggregates CS across blocks", {
    e1 <- .sc_makeFineMappingRow(3)
    e2 <- .sc_makeFineMappingRow(2)
    res <- GwasFineMappingResult(
        studyName = c("g1", "g1"),
        method = c("susie", "susie"),
        blockId = c("chr1:1-100", "chr1:200-300"),
        entry = list(e1, e2)
    )
    cs <- credibleSets(res)
    expect_equal(nrow(cs), 4L) # 2 CS members per block
    expect_equal(
        cs$blockId,
        c("chr1:1-100", "chr1:1-100", "chr1:200-300", "chr1:200-300")
    )
    expect_true(all(is.na(cs$context)))
    # region= pins one block -> bare table
    bare <- credibleSets(res, studyName = "g1", method = "susie", region = "chr1:1-100")
    expect_false("blockId" %in% names(bare))
})


test_that("GwasFineMappingResult: .tupleSelectRowGwasFmr requires both selectors for multi-row", {
    e <- .ca_makeFmEntry(3)
    res <- GwasFineMappingResult(
        studyName = c("g1", "g2"),
        method = c("susie", "susie"),
        entry = list(e, e)
    )
    expect_error(pip(res), "Pass `study` and `method`")
    expect_error(
        pip(res, studyName = c("g1", "g2"), method = "susie"),
        "Must have length 1"
    )
    expect_error(pip(res, studyName = "ghost", method = "susie"), "No entry for")
})


test_that("GwasFineMappingResult: studyName/methodNames inherit from base", {
    e <- .ca_makeFmEntry(3)
    res <- GwasFineMappingResult(
        studyName = c("g1", "g2"),
        method = c("susie", "susieRss"),
        entry = list(e, e)
    )
    expect_setequal(studyName(res), c("g1", "g2"))
    expect_setequal(methodNames(res), c("susie", "susieRss"))
})


test_that("GwasFineMappingResult: marginalEffects with study/method selectors", {
    e1 <- .ca_makeFmEntry(3)
    e2 <- .ca_makeFmEntry(4)
    res <- GwasFineMappingResult(
        studyName = c("g1", "g2"),
        method = c("susie", "susie"),
        entry = list(e1, e2)
    )
    # Collection-level selection picks the g2 entry, then delegates to the
    # entry-level getMarginalEffects.
    me <- marginalEffects(res, studyName = "g2", method = "susie")
    expect_s3_class(me, "data.frame")
    expect_equal(nrow(me), 4L)
    expect_true(all(c("variant_id", "beta", "se", "z", "p") %in% names(me)))
})

# ===========================================================================
# Validity: the messages that name what is missing
#
# Every test above builds a valid collection, so the validity function's
# early-return branches were never executed. They are reachable by dropping a
# column from a built object, which is what a careless mcols edit would do.
# ===========================================================================

test_that("validity names a missing identity column", {
    res <- GwasFineMappingResult(
        studyName = c("g1", "g2"),
        method = c("susie", "susie"),
        entry = list(.sc_makeFineMappingRow(3), .sc_makeFineMappingRow(3))
    )
    bad <- res
    mcols(bad)$method <- NULL
    expect_error(
        methods::validObject(bad),
        "missing elements \\{'method'\\}"
    )
})

test_that("validity names a missing entry payload column", {
    res <- GwasFineMappingResult(
        studyName = "g1",
        method = "susie",
        entry = list(.sc_makeFineMappingRow(3))
    )
    bad <- res
    mcols(bad)$cvResult <- NULL
    expect_error(
        methods::validObject(bad),
        "missing entry payload columns: .*missing elements \\{'cvResult'\\}"
    )
})

test_that("GwasFineMappingResult: pip(returnList) keys by study|method", {
    # Documented as "a per-entry list keyed by identity tuple". The QTL
    # method always honoured it; this one used to ignore the flag and hand
    # back the flat vector, so the list branch had no test.
    res <- GwasFineMappingResult(
        studyName = c("g1", "g2"),
        method = c("susie", "susie"),
        entry = list(.ca_makeFmEntry(3), .ca_makeFmEntry(4))
    )
    got <- pip(res, studyName = "g2", method = "susie", returnList = TRUE)
    expect_type(got, "list")
    expect_named(got, "g2|susie")
    expect_equal(length(got[["g2|susie"]]), 4L)
    # The flat vector is still what you get without the flag.
    expect_false(is.list(pip(res, studyName = "g2", method = "susie")))
})

test_that("GwasFineMappingResult: context / trait selectors are refused", {
    # They belong to the QTL axis. The shared accessor signature keeps them,
    # but a GWAS collection has no such axis, so asking is a mistake rather
    # than a silent no-op.
    res <- GwasFineMappingResult(
        studyName = "g1",
        method = "susie",
        entry = list(.ca_makeFmEntry(3))
    )
    expect_error(
        pip(res, studyName = "g1", method = "susie", context = "brain"),
        "has no context or trait axis"
    )
    expect_error(
        pip(res, studyName = "g1", method = "susie", context = "brain"),
        "`context` does not select anything"
    )
    expect_error(
        pip(res, studyName = "g1", method = "susie", trait = "ENSG_A"),
        "`trait` does not select anything"
    )
    # Both at once: the message pluralises rather than naming one.
    expect_error(
        pip(
            res,
            studyName = "g1",
            method = "susie",
            context = "brain",
            trait = "ENSG_A"
        ),
        "`context` and `trait` do not select anything"
    )
    # NA is an absence, not a request: callers threading a whole
    # (study, context, trait, method) record fill absent axes with NA.
    expect_equal(
        length(pip(
            res,
            studyName = "g1",
            method = "susie",
            context = NA,
            trait = NA
        )),
        3L
    )
})
