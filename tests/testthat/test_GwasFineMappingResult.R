# Tests for R/GwasFineMappingResult.R

# === Tests migrated from test_s4Constructors.R (GwasFineMappingResult) ===

test_that("GwasFineMappingResult: builds a collection keyed by 2-tuple", {
    e1 <- .sc_makeFineMappingEntry(3)
    e2 <- .sc_makeFineMappingEntry(3)
    res <- GwasFineMappingResult(
        study = c("g1", "g2"),
        method = c("susie", "susie"),
        entry = list(e1, e2)
    )
    expect_s4_class(res, "GwasFineMappingResult")
    expect_equal(nrow(res), 2L)
})

test_that("GwasFineMappingResult: region is auto-derived from region_id", {
    e1 <- .sc_makeFineMappingEntry(3)
    e2 <- .sc_makeFineMappingEntry(2)
    res <- GwasFineMappingResult(
        study = c("g1", "g1"),
        method = c("susie", "susie"),
        region_id = c("chr1_1_500", "chr2_600_900"),
        entry = list(e1, e2)
    )
    reg <- getRegion(res)
    expect_equal(as.character(GenomicRanges::seqnames(reg)), c("chr1", "chr2"))
    expect_equal(GenomicRanges::start(reg), c(1L, 600L))
    expect_equal(GenomicRanges::end(reg), c(500L, 900L))
    # synthetic region_id (no coordinates) -> a chrUn 0-width sentinel
    res2 <- GwasFineMappingResult(
        study = "g1",
        method = "susie",
        entry = list(e1)
    )
    expect_equal(
        as.character(GenomicRanges::seqnames(getRegion(res2))),
        "chrUn"
    )
})


test_that("GwasFineMappingResult: validity does not recurse on key subset (#546)", {
    # The validity method builds a key-column data.frame to check tuple
    # uniqueness. Doing so via `object[, keyCols]` preserves the
    # GwasFineMappingResult class while dropping the required `entry` column;
    # older S4Vectors revalidates that intermediate and fails with
    # "missing columns: entry".
    e <- .sc_makeFineMappingEntry(3)
    res <- GwasFineMappingResult(
        study = "g1",
        method = "susie",
        entry = list(e)
    )
    expect_s4_class(res, "GwasFineMappingResult")
    expect_true(validObject(res))

    sub <- res[, c("study", "method", "region_id")]
    expect_false("entry" %in% names(sub))
})


test_that("GwasFineMappingResult: errors on length mismatch", {
    e <- .sc_makeFineMappingEntry(3)
    expect_error(
        GwasFineMappingResult(
            study = c("g1", "g2"),
            method = c("susie"),
            entry = list(e)
        ),
        "same length"
    )
})


test_that("GwasFineMappingResult: same (study, method) across distinct region_ids OK", {
    # The default constructor synthesises region_1, region_2, ... so the
    # (study, method, region_id) triple is unique even when (study, method)
    # repeats. This is the genome-wide-across-blocks shape that
    # qtlEnrichmentPipeline + colocPipeline expect.
    e1 <- .sc_makeFineMappingEntry(3)
    e2 <- .sc_makeFineMappingEntry(3)
    res <- GwasFineMappingResult(
        study = c("g1", "g1"),
        method = c("susie", "susie"),
        entry = list(e1, e2)
    )
    expect_s4_class(res, "GwasFineMappingResult")
    expect_equal(nrow(res), 2L)
    expect_setequal(as.character(res$region_id), c("region_1", "region_2"))
})


test_that("GwasFineMappingResult: rejects duplicate (study, method, region_id) tuples", {
    e1 <- .sc_makeFineMappingEntry(3)
    e2 <- .sc_makeFineMappingEntry(3)
    expect_error(
        GwasFineMappingResult(
            study = c("g1", "g1"),
            method = c("susie", "susie"),
            region_id = c("chr22_1_100", "chr22_1_100"),
            entry = list(e1, e2)
        ),
        "uniqueness violated"
    )
})


test_that("GwasFineMappingResult: errors when region_id length mismatches", {
    e1 <- .sc_makeFineMappingEntry(3)
    e2 <- .sc_makeFineMappingEntry(3)
    expect_error(
        GwasFineMappingResult(
            study = c("g1", "g1"),
            method = c("susie", "susie"),
            region_id = "only_one",
            entry = list(e1, e2)
        ),
        "same length"
    )
})


test_that("GwasFineMappingResult: show prints summary", {
    e <- .sc_makeFineMappingEntry(3)
    res <- GwasFineMappingResult(
        study = "g1",
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
        study = c("g1", "g1"),
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
        study = "g1",
        method = "susie",
        entry = list(.sh_makeFmEntry()),
        ldSketch = .sh_makeGenotypeHandle()
    )
    out <- capture.output(show(res))
    expect_true(any(grepl("LD sketch: gds @ /tmp/test.gds", out)))
})


# === Tests migrated from test_collectionAccessors.R (GwasFineMappingResult) ===

test_that("GwasFineMappingResult: getPip with study/method selectors", {
    e1 <- .ca_makeFmEntry(3)
    e2 <- .ca_makeFmEntry(4)
    res <- GwasFineMappingResult(
        study = c("g1", "g2"),
        method = c("susie", "susie"),
        entry = list(e1, e2)
    )
    pip <- getPip(res, study = "g2", method = "susie")
    expect_equal(length(pip), 4L)
})


test_that("GwasFineMappingResult: getContexts/getTraits return NULL", {
    e <- .ca_makeFmEntry(3)
    res <- GwasFineMappingResult(
        study = "g1",
        method = "susie",
        entry = list(e)
    )
    expect_null(getContexts(res))
    expect_null(getTraits(res))
})


test_that("GwasFineMappingResult: getCs/getTopLoci/getSusieFit/getVariantIds dispatch", {
    e <- .ca_makeFmEntry(3)
    res <- GwasFineMappingResult(
        study = "g1",
        method = "susie",
        entry = list(e)
    )
    expect_equal(nrow(getCs(res)), 2L)
    # getTopLoci returns the projected posterior view (filtered by default
    # signalCutoff = 0.025; .ca_makeTopLoci sets all pip > 0.025 so all rows
    # survive). Compare on the projected shape, not the slot's raw shape.
    tl <- getTopLoci(res, signalCutoff = 0)
    expect_equal(nrow(tl), 3L)
    expect_equal(tl$variant_id, .ca_makeTopLoci(3)$variant_id)
    expect_equal(getSusieFit(res), list(payload = "fit_n=3"))
    expect_equal(length(getVariantIds(res)), 3L)
})

test_that("GwasFineMappingResult: getTopLoci aggregates per-block rows genome-wide", {
    # A genome-wide collection: same (study, method) across two region blocks.
    # With no selectors getTopLoci now stacks both blocks, tagging each variant
    # with its region_id; context/trait are NA-filled (GWAS keys on region).
    e1 <- .sc_makeFineMappingEntry(3)
    e2 <- .sc_makeFineMappingEntry(2)
    res <- GwasFineMappingResult(
        study = c("g1", "g1"),
        method = c("susie", "susie"),
        region_id = c("chr1:1-100", "chr1:200-300"),
        entry = list(e1, e2)
    )
    agg <- getTopLoci(res, signalCutoff = 0)
    expect_equal(nrow(agg), 5L)
    expect_equal(
        agg$region_id,
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
    expect_true("variant_id" %in% names(agg))
})

test_that("GwasFineMappingResult: getTopLoci region= selects a single block", {
    e1 <- .sc_makeFineMappingEntry(3)
    e2 <- .sc_makeFineMappingEntry(2)
    res <- GwasFineMappingResult(
        study = c("g1", "g1"),
        method = c("susie", "susie"),
        region_id = c("chr1:1-100", "chr1:200-300"),
        entry = list(e1, e2)
    )
    # region= pins one row, so this hits the single-entry fast path (bare table).
    tl <- getTopLoci(
        res,
        study = "g1",
        method = "susie",
        region = "chr1:200-300",
        signalCutoff = 0
    )
    expect_equal(nrow(tl), 2L)
    expect_false("region_id" %in% names(tl))
})

test_that("GwasFineMappingResult: getCs aggregates per-block credible sets tagged by region_id", {
    e1 <- .sc_makeFineMappingEntry(3)
    e2 <- .sc_makeFineMappingEntry(2)
    res <- GwasFineMappingResult(
        study = c("g1", "g1"),
        method = c("susie", "susie"),
        region_id = c("chr1:1-100", "chr1:200-300"),
        entry = list(e1, e2)
    )
    cs <- getCs(res)
    expect_equal(nrow(cs), 4L) # 2 CS members per block
    expect_equal(
        cs$region_id,
        c("chr1:1-100", "chr1:1-100", "chr1:200-300", "chr1:200-300")
    )
    expect_true(all(is.na(cs$context)))
    # region= pins one block -> bare table
    bare <- getCs(res, study = "g1", method = "susie", region = "chr1:1-100")
    expect_false("region_id" %in% names(bare))
})


test_that("GwasFineMappingResult: .tupleSelectRowGwasFmr requires both selectors for multi-row", {
    e <- .ca_makeFmEntry(3)
    res <- GwasFineMappingResult(
        study = c("g1", "g2"),
        method = c("susie", "susie"),
        entry = list(e, e)
    )
    expect_error(getPip(res), "Pass `study` and `method`")
    expect_error(
        getPip(res, study = c("g1", "g2"), method = "susie"),
        "must each be length 1"
    )
    expect_error(getPip(res, study = "ghost", method = "susie"), "No entry for")
})


test_that("GwasFineMappingResult: getStudy/getMethodNames inherit from base", {
    e <- .ca_makeFmEntry(3)
    res <- GwasFineMappingResult(
        study = c("g1", "g2"),
        method = c("susie", "susieRss"),
        entry = list(e, e)
    )
    expect_setequal(getStudy(res), c("g1", "g2"))
    expect_setequal(getMethodNames(res), c("susie", "susieRss"))
})


test_that("GwasFineMappingResult: getMarginalEffects with study/method selectors", {
    e1 <- .ca_makeFmEntry(3)
    e2 <- .ca_makeFmEntry(4)
    res <- GwasFineMappingResult(
        study = c("g1", "g2"),
        method = c("susie", "susie"),
        entry = list(e1, e2)
    )
    # Collection-level selection picks the g2 entry, then delegates to the
    # entry-level getMarginalEffects.
    me <- getMarginalEffects(res, study = "g2", method = "susie")
    expect_s3_class(me, "data.frame")
    expect_equal(nrow(me), 4L)
    expect_true(all(c("variant_id", "beta", "se", "z", "p") %in% names(me)))
})
