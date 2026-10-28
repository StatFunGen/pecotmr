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

test_that("GwasFineMappingResult: a window comes from region, not blockId", {
    # A blockId is a name, however it is spelled: it gives no window, so the
    # rows keep their spans and nothing was given.
    e1 <- .sc_makeFineMappingRow(3)
    e2 <- .sc_makeFineMappingRow(2, offset = 1000L)
    res <- GwasFineMappingResult(
        studyName = c("g1", "g1"),
        method = c("susie", "susie"),
        blockId = c("chr1_1_500", "b2"),
        entry = list(e1, e2)
    )
    expect_identical(res$blockId, c("chr1_1_500", "b2"))
    expect_identical(
        as.character(genomicRegion(res)),
        c("chr1:100-300", "chr1:1100-1200")
    )
    expect_identical(res$regionGiven, c(FALSE, FALSE))
    # A window that leaves out its row's variants is not that row's window.
    expect_error(
        GwasFineMappingResult(
            studyName = c("g1", "g1"),
            method = c("susie", "susie"),
            entry = list(e1, e2),
            region = GenomicRanges::GRanges(
                c("chr1", "chr2"),
                IRanges::IRanges(c(1L, 600L), c(500L, 900L))
            )
        ),
        "row 2: range chr2:600-900, variants at chr1:1100-1200",
        fixed = TRUE
    )
    expect_error(
        GwasFineMappingResult(
            studyName = "g1",
            method = "susie",
            entry = list(e1),
            region = GenomicRanges::GRanges(
                "chr1",
                IRanges::IRanges(150L, 400L)
            )
        ),
        "must contain the row's variants"
    )
    # With no range given there is no region column: genomicRegion() reads
    # the span of the variants, which cannot drift out of step with them.
    res2 <- GwasFineMappingResult(
        studyName = "g1",
        method = "susie",
        entry = list(e1)
    )
    expect_false("region" %in% colnames(S4Vectors::mcols(res2)))
    expect_false(res2$regionGiven)
    expect_identical(as.character(genomicRegion(res2)), "chr1:100-300")
    # A `region` is given, whatever the blockId.
    res3 <- GwasFineMappingResult(
        studyName = "g1",
        method = "susie",
        entry = list(e1),
        blockId = "b1",
        region = GenomicRanges::GRanges("chr1", IRanges::IRanges(50L, 400L))
    )
    expect_true(res3$regionGiven)
    expect_identical(as.character(genomicRegion(res3)), "chr1:50-400")
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
            blockId = c("blockA", "blockB"),
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
    expect_true(any(grepl("LD reference: NULL", out)))
})


test_that("show.GwasFineMappingResult reports the ldReference source", {
    res <- GwasFineMappingResult(
        studyName = "g1",
        method = "susie",
        entry = list(.sh_makeFmEntry()),
        ldReference = .sh_makeGenotypeHandle()
    )
    out <- capture.output(show(res))
    expect_true(any(grepl("LD reference: plink2 @ /tmp/test", out)))
})


# === Tests migrated from test_collectionAccessors.R (GwasFineMappingResult) ===

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
    expect_equal(length(credibleSets(res)), 2L)
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
    e2 <- .sc_makeFineMappingRow(2, offset = 1000L)
    res <- GwasFineMappingResult(
        studyName = c("g1", "g1"),
        method = c("susie", "susie"),
        blockId = c("chr1:1-500", "chr1:501-2000"),
        entry = list(e1, e2)
    )
    agg <- topLoci(res, signalCutoff = 0)
    expect_equal(length(agg), 5L)
    expect_equal(
        agg$blockId,
        c(
            "chr1:1-500",
            "chr1:1-500",
            "chr1:1-500",
            "chr1:501-2000",
            "chr1:501-2000"
        )
    )
    expect_true(all(is.na(agg$context)))
    expect_true(all(is.na(agg$trait)))
    expect_true("variant_id" %in% names(S4Vectors::mcols(agg)))
})

test_that("GwasFineMappingResult: topLoci picks one block by name or window", {
    e1 <- .sc_makeFineMappingRow(3)
    e2 <- .sc_makeFineMappingRow(2, offset = 1000L)
    res <- GwasFineMappingResult(
        studyName = c("g1", "g1"),
        method = c("susie", "susie"),
        blockId = c("chr1:1-500", "chr1:501-2000"),
        entry = list(e1, e2)
    )
    # blockId= pins one row, so this hits the single-entry fast path (bare
    # table).
    tl <- topLoci(
        res,
        studyName = "g1",
        method = "susie",
        blockId = "chr1:501-2000",
        signalCutoff = 0
    )
    expect_equal(length(tl), 2L)
    expect_false("blockId" %in% names(S4Vectors::mcols(tl)))
    # A window selects rows by subsetting the collection on its rows'
    # analysis windows -- here, alone, the second block, whose variants span
    # chr1:1100-1200.
    byWindow <- topLoci(
        res[IRanges::overlapsAny(
            genomicRegion(res),
            GenomicRanges::GRanges("chr1:1150-1500")
        )],
        signalCutoff = 0
    )
    expect_identical(as.character(byWindow), as.character(tl))
    # A window over both blocks keeps both, which then aggregate.
    both <- topLoci(
        res[IRanges::overlapsAny(
            genomicRegion(res),
            GenomicRanges::GRanges("chr1:1-5000")
        )],
        signalCutoff = 0
    )
    expect_equal(length(both), 5L)
    expect_true("blockId" %in% names(S4Vectors::mcols(both)))
})

test_that("GwasFineMappingResult: credibleSets aggregates CS across blocks", {
    e1 <- .sc_makeFineMappingRow(3)
    e2 <- .sc_makeFineMappingRow(2, offset = 1000L)
    res <- GwasFineMappingResult(
        studyName = c("g1", "g1"),
        method = c("susie", "susie"),
        blockId = c("chr1:1-500", "chr1:501-2000"),
        entry = list(e1, e2)
    )
    cs <- credibleSets(res)
    expect_equal(length(cs), 4L) # 2 CS members per block
    expect_equal(
        cs$blockId,
        c("chr1:1-500", "chr1:1-500", "chr1:501-2000", "chr1:501-2000")
    )
    expect_true(all(is.na(cs$context)))
    # blockId= pins one block -> bare table
    bare <- credibleSets(
        res,
        studyName = "g1",
        method = "susie",
        blockId = "chr1:1-500"
    )
    expect_false("blockId" %in% names(bare))
})


test_that("GwasFineMappingResult: .tupleSelectRowGwasFmr requires both selectors for multi-row", {
    e <- .ca_makeFmEntry(3)
    res <- GwasFineMappingResult(
        studyName = c("g1", "g2"),
        method = c("susie", "susie"),
        entry = list(e, e)
    )
    expect_error(pip(res), "Pass `studyName` and `method`")
    expect_error(
        pip(res, studyName = c("g1", "g2"), method = "susie"),
        "Must have length 1"
    )
    expect_error(
        pip(res, studyName = "ghost", method = "susie"),
        "No entry for"
    )
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
    expect_s4_class(me, "GRanges")
    expect_equal(length(me), 4L)
    expect_true(all(
        c("variant_id", "beta", "se", "z", "p") %in% names(S4Vectors::mcols(me))
    ))
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

test_that("GwasFineMappingResult: pip() reads the row the selectors pick", {
    res <- GwasFineMappingResult(
        studyName = c("g1", "g2"),
        method = c("susie", "susie"),
        entry = list(.ca_makeFmEntry(3), .ca_makeFmEntry(4))
    )
    got <- pip(res, studyName = "g2", method = "susie")
    expect_type(got, "double")
    expect_length(got, 4L)
    # One return type: the list form is gone with its switch.
    expect_error(pip(res, studyName = "g2", returnList = TRUE), "unused")
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

# ===========================================================================
# A study fine-mapped block by block
# ===========================================================================

# One study fine-mapped in two LD blocks with two methods, plus a second
# study in one block.
.gfr_makeBlocked <- function() {
    GwasFineMappingResult(
        studyName = c("g1", "g1", "g1", "g1", "g2"),
        method = c("susie", "susie", "ser", "ser", "susie"),
        entry = list(
            .sc_makeFineMappingRow(3),
            .sc_makeFineMappingRow(3, offset = 10000L),
            .sc_makeFineMappingRow(3),
            .sc_makeFineMappingRow(3, offset = 10000L),
            .sc_makeFineMappingRow(3)
        ),
        blockId = c(
            "chr1_1_5000",
            "chr1_5001_20000",
            "chr1_1_5000",
            "chr1_5001_20000",
            "chr1_1_5000"
        )
    )
}

test_that("pip(GwasFineMappingResult) picks a block by name or window", {
    x <- .gfr_makeBlocked()
    named <- pip(
        x,
        studyName = "g1",
        method = "susie",
        blockId = "chr1_5001_20000"
    )
    near <- pip(
        x[IRanges::overlapsAny(
            genomicRegion(x),
            GenomicRanges::GRanges("chr1:10100-10150")
        )],
        studyName = "g1",
        method = "susie"
    )
    expect_identical(near, named)
    expect_length(named, 3L)
    expect_error(
        pip(x, studyName = "g1", method = "susie"),
        "pass `blockId` to disambiguate"
    )
})

test_that("fineMappingResult returns every block of a fit, or the one named", {
    x <- .gfr_makeBlocked()
    fit <- fineMappingResult(x, studyName = "g1", method = "susie")
    expect_s4_class(fit, "GwasFineMappingResult")
    expect_equal(nrow(fit), 2L)
    expect_identical(fit$blockId, c("chr1_1_5000", "chr1_5001_20000"))
    one <- fineMappingResult(
        x,
        studyName = "g1",
        method = "susie",
        blockId = "chr1_5001_20000"
    )
    expect_equal(nrow(one), 1L)
    expect_identical(one$blockId, "chr1_5001_20000")
    # Subsetting on the rows' windows keeps the blocks a window overlaps.
    near <- fineMappingResult(
        x[IRanges::overlapsAny(
            genomicRegion(x),
            GenomicRanges::GRanges("chr1:10100-10150")
        )],
        studyName = "g1",
        method = "susie"
    )
    expect_identical(near$blockId, "chr1_5001_20000")
    # Selectors left unset are fine while they still name one fit.
    expect_equal(nrow(fineMappingResult(x, studyName = "g2")), 1L)
    expect_error(
        fineMappingResult(x, studyName = "g1"),
        "holds 2 \\(study, method\\) fits matching the selection"
    )
    expect_error(
        fineMappingResult(x, studyName = "g1", method = "susie", blockId = "z"),
        "No entry for \\(study='g1', method='susie', blockId='z'\\)"
    )
    expect_error(
        fineMappingResult(x, method = "lasso"),
        "No entry for \\(study='<any>', method='lasso'\\)"
    )
    expect_error(
        fineMappingResult(x, studyName = "g1", context = "brain"),
        "no context or trait axis"
    )
    expect_error(
        fineMappingResult(x[0L], studyName = "g1"),
        "has no rows"
    )
})

test_that("a fit carried by a multi-chromosome entry is refused", {
    # Splitting the entry by chromosome would copy its one fit onto every
    # piece, and each piece would then be read as the whole fit.
    ids <- c("chr1:100:A:G", "chr2:100:A:G")
    tl <- data.frame(variant_id = ids, pip = c(0.6, 0.4))
    fitted <- FineMappingRow(
        variantIds = ids,
        susieFit = list(pip = setNames(c(0.6, 0.4), ids)),
        topLoci = tl
    )
    expect_error(
        GwasFineMappingResult(
            studyName = "g1",
            method = "susie",
            entry = list(fitted)
        ),
        "`entry` element\\(s\\) 1 span more than one chromosome but carry"
    )
    # Without a fit there is nothing to copy, and the entry splits as before.
    bare <- FineMappingRow(variantIds = ids, susieFit = list(), topLoci = tl)
    split <- GwasFineMappingResult(
        studyName = "g1",
        method = "susie",
        entry = list(bare)
    )
    expect_equal(nrow(split), 2L)
})

test_that("filtering every variant out of a fit's blocks keeps them apart", {
    skip_if_not_installed("plyranges")
    x <- .gfr_makeBlocked()
    # Each block becomes an empty row; its blockId keeps it distinct.
    emptied <- dplyr::filter(x, pip > 2)
    expect_equal(unname(lengths(emptied)), rep(0L, 5L))
    expect_identical(emptied$blockId, x$blockId)
})

test_that("a row without a blockId is named by its own span", {
    # Row identity is (study, method, range), so one (study, method) may hold
    # two blocks over other spans -- the genome-wide-across-blocks shape that
    # qtlEnrichmentPipeline + colocPipeline expect. The caller named neither
    # row; each is named by the stretch its variants cover, the range its
    # identity already keys on.
    x <- GwasFineMappingResult(
        studyName = c("g1", "g1"),
        method = c("susie", "susie"),
        entry = list(
            .sc_makeFineMappingRow(3),
            .sc_makeFineMappingRow(3, offset = 10000L)
        )
    )
    expect_identical(x$blockId, c("chr1_100_300", "chr1_10100_10300"))
    expect_identical(x$blockId, .rtlRangeKeys(x))
    # A span filled in as the id was computed, not given.
    expect_identical(x$regionGiven, c(FALSE, FALSE))
    # A supplied id is kept, and an NA one filled.
    y <- GwasFineMappingResult(
        studyName = c("g1", "g1"),
        method = c("susie", "susie"),
        entry = list(
            .sc_makeFineMappingRow(3),
            .sc_makeFineMappingRow(3, offset = 10000L)
        ),
        blockId = c("blockA", NA)
    )
    expect_identical(y$blockId, c("blockA", "chr1_10100_10300"))
    # Each chromosome piece of a fit-less entry gets its own span.
    ids <- c("chr1:100:A:G", "chr2:500:A:G")
    split <- GwasFineMappingResult(
        studyName = "g1",
        method = "susie",
        entry = list(FineMappingRow(
            ids,
            susieFit = list(),
            topLoci = data.frame(variant_id = ids, pip = c(0.6, 0.4))
        ))
    )
    expect_identical(split$blockId, c("chr1_100_100", "chr2_500_500"))
})

test_that("GwasFineMappingResult() keeps the caller's metadata", {
    res <- GwasFineMappingResult(
        studyName = "g1",
        method = "susie",
        entry = list(.sc_makeFineMappingRow(3)),
        metadata = list(run = "r1")
    )
    expect_equal(S4Vectors::metadata(res), list(run = "r1"))
})

test_that("a fit-less entry split across chromosomes records no fit variants", {
    # Its rows hold one chromosome each; with no fit there is nothing for a
    # row to have narrowed, so none of them warns.
    ids <- c("chr1:100:A:G", "chr2:200:C:T")
    tl <- data.frame(variant_id = ids, pip = c(0.6, 0.4))
    res <- GwasFineMappingResult(
        studyName = "g1",
        method = "susie",
        entry = list(FineMappingRow(ids, NULL, tl))
    )
    expect_equal(nrow(res), 2L)
    expect_equal(unname(lengths(res$fitVariants)), c(0L, 0L))
    expect_no_warning(lbf(res))
})
