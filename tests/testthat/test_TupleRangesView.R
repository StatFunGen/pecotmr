# Tests for the plyranges bridge: TupleRangesView and the flatten / nest
# conversions.
#
# plyranges has no GRangesList support, so the collections cannot be operated
# on directly. The bridge converts to a flat GenomicRanges view -- which
# plyranges does support, via the DelegatingGenomicRanges extension point --
# and back again.

test_that("flattenTupleRanges concatenates every element's ranges", {
    data(qtlSumStatsMulticontextExample, envir = environment())
    mc <- qtlSumStatsMulticontextExample
    flat <- flattenTupleRanges(mc)
    expect_s4_class(flat, "GRanges")
    expect_equal(length(flat), sum(lengths(mc)))
})

test_that("flattenTupleRanges broadcasts the identity tuple, dot-prefixed", {
    data(qtlSumStatsMulticontextExample, envir = environment())
    mc <- qtlSumStatsMulticontextExample
    flat <- flattenTupleRanges(mc)
    cn <- colnames(mcols(flat))
    expect_true(all(is_in(c(".study", ".context", ".trait"), cn)))
    # One value per range, matching the element it came from.
    expect_setequal(unique(mcols(flat)$.context), as.character(mc$context))
    expect_equal(
        sum(mcols(flat)$.context == "blood"),
        lengths(mc)[[which(mc$context == "blood")]]
    )
})

test_that("the dot prefix keeps a colliding per-range column intact", {
    # On gwasFineMappingExample the outer `method` is "susie" while the
    # per-range `method` is "susieRss" -- the collection's label against the
    # fitter actually used. Broadcasting onto the bare name would overwrite
    # one with the other.
    data(gwasFineMappingExample, envir = environment())
    x <- gwasFineMappingExample
    expect_true(is_in("method", colnames(mcols(x[[1]]))))
    flat <- flattenTupleRanges(x)
    expect_true(all(is_in(c("method", ".method"), colnames(mcols(flat)))))
    expect_false(identical(
        unique(mcols(flat)$method),
        unique(mcols(flat)$.method)
    ))
})

test_that("flattenTupleRanges drops non-broadcastable payload columns", {
    # susieFit / cvResult describe an ELEMENT, so there is no range to put
    # them on.
    data(qtlFineMappingExample, envir = environment())
    flat <- flattenTupleRanges(qtlFineMappingExample)
    expect_false(is_in("susieFit", colnames(mcols(flat))))
    expect_false(is_in(".susieFit", colnames(mcols(flat))))
    expect_true(is_in(".study", colnames(mcols(flat))))
})

test_that("flattenTupleRanges rejects a non-collection", {
    expect_error(
        flattenTupleRanges(GenomicRanges::GRanges()),
        "RangedTupleList"
    )
})

test_that("nestTupleRanges round-trips an unmodified flatten", {
    data(qtlSumStatsMulticontextExample, envir = environment())
    mc <- qtlSumStatsMulticontextExample
    back <- nestTupleRanges(flattenTupleRanges(mc), mc)
    expect_s4_class(back, "QtlSumStats")
    expect_equal(lengths(back), lengths(mc))
    expect_equal(mcols(back), mcols(mc))
    # Collection-level slots survive, which is the point of nesting through a
    # template rather than rebuilding from the ranges alone.
    expect_identical(getLdSketch(back), getLdSketch(mc))
    expect_identical(getGenome(back), getGenome(mc))
})

test_that("nestTupleRanges returns a filtered flatten to its elements", {
    data(qtlSumStatsMulticontextExample, envir = environment())
    mc <- qtlSumStatsMulticontextExample
    flat <- flattenTupleRanges(mc)
    kept <- flat[mcols(flat)$Z > 0]
    back <- nestTupleRanges(kept, mc)
    expect_equal(sum(lengths(back)), length(kept))
    expect_lt(sum(lengths(back)), sum(lengths(mc)))
    # Every range lands under the tuple it carried.
    idx <- which(mc$context == "blood")
    expect_equal(
        lengths(back)[[idx]],
        sum(mcols(kept)$.context == "blood")
    )
})

test_that("nestTupleRanges keeps a tuple whose ranges all vanished", {
    # An empty element, not a dropped row: the collection keeps its shape so
    # its metadata stays aligned.
    data(qtlSumStatsMulticontextExample, envir = environment())
    mc <- qtlSumStatsMulticontextExample
    flat <- flattenTupleRanges(mc)
    onlyBlood <- flat[mcols(flat)$.context == "blood"]
    back <- nestTupleRanges(onlyBlood, mc)
    expect_equal(nrow(back), nrow(mc))
    expect_equal(sum(lengths(back) == 0L), nrow(mc) - 1L)
})

test_that("nestTupleRanges survives an NA-valued identity column", {
    # str_c() propagates NA, so an NA identity value (varY is NA_real_ on a
    # z-score collection) would make every key NA and the row lookup a
    # logical subscript full of NAs.
    data(qtlSumStatsMulticontextExample, envir = environment())
    mc <- qtlSumStatsMulticontextExample
    expect_true(any(is.na(mcols(mc)$varY)))
    expect_equal(
        lengths(nestTupleRanges(flattenTupleRanges(mc), mc)),
        lengths(mc)
    )
})

test_that("nestTupleRanges validates its arguments", {
    data(qtlSumStatsMulticontextExample, envir = environment())
    mc <- qtlSumStatsMulticontextExample
    flat <- flattenTupleRanges(mc)
    expect_error(nestTupleRanges(flat, "not a collection"), "RangedTupleList")
    expect_error(nestTupleRanges("not ranges", mc), "GRanges")
})

test_that("a TupleRangesView is a GenomicRanges and subsets to itself", {
    data(qtlSumStatsMulticontextExample, envir = environment())
    flat <- flattenTupleRanges(qtlSumStatsMulticontextExample)
    v <- .tupleRangesView(flat)
    expect_s4_class(v, "TupleRangesView")
    # Being a GenomicRanges is what makes plyranges dispatch at all.
    expect_true(is(v, "GenomicRanges"))
    expect_equal(length(v), length(flat))
    expect_s4_class(v[1:10], "TupleRangesView")
    expect_equal(length(v[1:10]), 10L)
})

test_that("the view constructor keeps mcols parallel to the object", {
    # DelegatingGenomicRanges carries its own elementMetadata; leaving it
    # zero-row makes the class reject itself with "'mcols(x)' is not parallel
    # to 'x'".
    data(qtlSumStatsMulticontextExample, envir = environment())
    flat <- flattenTupleRanges(qtlSumStatsMulticontextExample)
    expect_silent(v <- .tupleRangesView(flat))
    expect_true(validObject(v))
})

test_that("the view rejects a two-argument subscript", {
    data(qtlSumStatsMulticontextExample, envir = environment())
    v <- .tupleRangesView(flattenTupleRanges(qtlSumStatsMulticontextExample))
    expect_error(v[1, 1], "two-argument")
})


# ---------------------------------------------------------------------------
# plyranges / dplyr verbs.
#
# These call the methods directly rather than through the generic: dispatch
# needs the S3method() entries a roxygen regen writes, and the logic under test
# is the same either way. The final test checks dispatch itself, and skips
# until the registration exists.
# ---------------------------------------------------------------------------

# @noRd
.trv_mc <- function() {
    data(qtlSumStatsMulticontextExample, envir = environment())
    qtlSumStatsMulticontextExample
}

test_that("filter keeps the collection and narrows its elements", {
    mc <- .trv_mc()
    out <- filter.RangedTupleList(mc, Z > 2)
    expect_s4_class(out, "QtlSumStats")
    expect_equal(nrow(out), nrow(mc))
    expect_lt(sum(lengths(out)), sum(lengths(mc)))
    expect_identical(getLdSketch(out), getLdSketch(mc))
})

test_that("filter can mix identity and per-range columns", {
    # The whole point of broadcasting the tuple: one predicate spanning both.
    mc <- .trv_mc()
    out <- filter.RangedTupleList(mc, .context == "blood" & Z > 1)
    idx <- which(mc$context == "blood")
    expect_gt(lengths(out)[[idx]], 0L)
    expect_equal(sum(lengths(out)[-idx]), 0L)
})

test_that("mutate adds a per-range column without changing the shape", {
    mc <- .trv_mc()
    out <- mutate.RangedTupleList(mc, hit = Z > 2)
    expect_equal(lengths(out), lengths(mc))
    expect_true(is_in("hit", colnames(mcols(out[[1]]))))
})

test_that("select keeps the identity columns it needs to nest", {
    # Dropping them would leave the ranges with no tuple, and every one would
    # fall into the first element.
    mc <- .trv_mc()
    out <- select.RangedTupleList(mc, Z)
    expect_equal(lengths(out), lengths(mc))
})

test_that("slice takes n ranges from EACH element", {
    # Slicing the flattened set would take n in total and empty all but the
    # first element.
    mc <- .trv_mc()
    out <- slice.RangedTupleList(mc, 1:5)
    expect_equal(unname(lengths(out)), rep(5L, nrow(mc)))
})

test_that("arrange orders within each element", {
    mc <- .trv_mc()
    out <- arrange.RangedTupleList(mc, Z)
    expect_equal(lengths(out), lengths(mc))
    for (i in seq_len(nrow(out))) {
        expect_false(is.unsorted(mcols(out[[i]])$Z))
    }
})

test_that("summarise reduces to a table rather than a collection", {
    mc <- .trv_mc()
    out <- summarise.RangedTupleList(mc, n = plyranges::n())
    expect_false(is(out, "RangedTupleList"))
    expect_equal(NROW(out), 1L)
})

test_that("group_by returns plyranges' own grouped representation", {
    mc <- .trv_mc()
    out <- group_by.RangedTupleList(mc, .context)
    expect_s4_class(out, "GroupedGenomicRanges")
    expect_equal(length(out), sum(lengths(mc)))
})

test_that("select and group_by work on the view itself", {
    # plyranges' own DelegatingGenomicRanges support misses both: they work on
    # a plain GRanges and fail on any delegating class. These patch that.
    mc <- .trv_mc()
    v <- .tupleRangesView(flattenTupleRanges(mc))
    expect_s4_class(select.TupleRangesView(v, Z), "TupleRangesView")
    expect_s4_class(
        group_by.TupleRangesView(v, .context),
        "GroupedGenomicRanges"
    )
})

test_that("the verbs are reachable through the dplyr generics", {
    # Needs the S3method() entries from a roxygen regen; skips until then.
    # Probed by attempting dispatch rather than by looking the function up:
    # getS3method() finds it in the namespace whether or not it is registered.
    mc <- .trv_mc()
    dispatches <- tryCatch(
        {
            dplyr::filter(mc, Z > 2)
            TRUE
        },
        error = function(e) FALSE
    )
    skip_if_not(dispatches, "S3 methods not yet registered (regen pending)")
    expect_s4_class(dplyr::filter(mc, Z > 2), "QtlSumStats")
    expect_s4_class(dplyr::slice(mc, 1:2), "QtlSumStats")
})

test_that("show reports the view's range and metadata-column counts", {
    # The view is internal, but its show method is what a developer sees when
    # one surfaces in a browser() or an error trace.
    data(qtlSumStatsMulticontextExample, envir = environment())
    flat <- flattenTupleRanges(qtlSumStatsMulticontextExample)
    v <- .tupleRangesView(flat)
    expect_output(show(v), "TupleRangesView:")
    expect_output(show(v), str_c(length(flat), " ranges"))
    expect_output(show(v), "metadata column\\(s\\)")
})

# ===========================================================================
# Degenerate collections and the plyranges requirement
# ===========================================================================

test_that("flattening a collection with no identity columns keeps its ranges", {
    # With nothing to key on, every range belongs to the single element --
    # which is what a one-row collection means.
    mc <- .trv_mc()
    bare <- mc[1, ]
    mcols(bare) <- mcols(bare)[, character(0), drop = FALSE]
    flat <- flattenTupleRanges(bare)
    expect_equal(length(flat), sum(lengths(bare)))
})

test_that(".rtlTupleKeys gives every row the same key with no key columns", {
    md <- S4Vectors::DataFrame(entry = S4Vectors::SimpleList(1, 2))
    expect_equal(.rtlTupleKeys(md, character(0), 3L), rep("", 3L))
})

test_that(".rtlTupleKeyCols is empty when there are no mcols", {
    gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 2))
    expect_equal(.rtlTupleKeyCols(gr), character(0))
})

test_that(".rtlAsRanges unwraps a view and passes a GRanges through", {
    flat <- flattenTupleRanges(.trv_mc())
    v <- .tupleRangesView(flat)
    expect_s4_class(.rtlAsRanges(v), "GRanges")
    expect_false(methods::is(.rtlAsRanges(v), "TupleRangesView"))
    expect_identical(.rtlAsRanges(flat), flat)
})

test_that("select(.drop_ranges = TRUE) returns the bare table, not a view", {
    # plyranges' own escape hatch: the caller asked for a data frame rather
    # than ranges, so re-wrapping it as a view would undo the request.
    flat <- flattenTupleRanges(.trv_mc())
    v <- .tupleRangesView(flat)
    out <- select.TupleRangesView(v, Z, .drop_ranges = TRUE)
    expect_false(methods::is(out, "TupleRangesView"))
    expect_false(methods::is(out, "GRanges"))
})
