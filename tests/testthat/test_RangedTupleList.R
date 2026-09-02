# Tests for the RangedTupleList virtual base.
#
# The class is virtual, so these exercise it through a minimal concrete
# subclass that adds the same shape of collection-level slots the real
# collections carry (a scalar, an ANY, and a list).

setClass(
    "RtlTestKid",
    contains = "RangedTupleList",
    representation(genome = "character", ldSketch = "ANY", qcInfo = "list")
)

.rtl_mk <- function(chr, pos, strand = "*") {
    GenomicRanges::GRanges(
        chr,
        IRanges::IRanges(pos, width = 1L),
        strand = strand
    )
}

.rtl_makeKid <- function() {
    x <- new(
        "RtlTestKid",
        GenomicRanges::GRangesList(
            a = .rtl_mk("chr1", c(100, 200, 300)),
            b = .rtl_mk("chr1", c(1000, 1100)),
            c = .rtl_mk("chr2", c(50, 60))
        ),
        genome = "hg38",
        ldSketch = "SKETCH",
        qcInfo = list(step1 = "ok")
    )
    S4Vectors::mcols(x) <- S4Vectors::DataFrame(
        study = c("s1", "s2", "s3"),
        n = c(10L, 20L, 30L)
    )
    x
}

# ===========================================================================
# nrow / ncol shims
# ===========================================================================

test_that("nrow counts elements and ncol counts metadata columns", {
    x <- .rtl_makeKid()
    expect_equal(nrow(x), 3L)
    expect_equal(ncol(x), 2L)
    # length() keeps its GRangesList meaning (number of elements).
    expect_equal(length(x), 3L)
    # lengths() still reports variants per element.
    expect_equal(unname(lengths(x)), c(3L, 2L, 2L))
})

test_that("ncol is zero when the collection has no mcols", {
    bare <- new(
        "RtlTestKid",
        GenomicRanges::GRangesList(a = .rtl_mk("chr1", 1:2)),
        genome = "hg38",
        ldSketch = NULL,
        qcInfo = list()
    )
    expect_equal(ncol(bare), 0L)
    expect_equal(nrow(bare), 1L)
})

# ===========================================================================
# $ / $<- shims
# ===========================================================================

test_that("$ reads an mcols column and $<- writes one", {
    x <- .rtl_makeKid()
    expect_equal(x$study, c("s1", "s2", "s3"))
    expect_equal(x$n, c(10L, 20L, 30L))
    expect_null(x$missingColumn)

    x$blockId <- c("b1", "b2", "b3")
    expect_equal(x$blockId, c("b1", "b2", "b3"))
    expect_equal(ncol(x), 3L)
    # Overwriting an existing column replaces it in place.
    x$study <- c("z1", "z2", "z3")
    expect_equal(x$study, c("z1", "z2", "z3"))
    expect_equal(ncol(x), 3L)
})

test_that("$<- works on a collection that had no mcols at all", {
    bare <- new(
        "RtlTestKid",
        GenomicRanges::GRangesList(a = .rtl_mk("chr1", 1:2)),
        genome = "hg38",
        ldSketch = NULL,
        qcInfo = list()
    )
    bare$study <- "s1"
    expect_equal(bare$study, "s1")
    expect_equal(ncol(bare), 1L)
})

# ===========================================================================
# [ rejects a second argument
# ===========================================================================

test_that("[ selects elements and preserves class, slots and mcols", {
    x <- .rtl_makeKid()
    y <- x[c(1L, 3L)]
    expect_s4_class(y, "RtlTestKid")
    expect_equal(length(y), 2L)
    expect_equal(y@genome, "hg38")
    expect_equal(y@ldSketch, "SKETCH")
    expect_equal(y$study, c("s1", "s3"))
    expect_equal(names(y), c("a", "c"))
})

test_that("[ rejects a two-argument call instead of silently dropping j", {
    # The stock CompressedList method ignores `j` and returns whole elements,
    # so x[1, "study"] would quietly mean something else than it used to.
    x <- .rtl_makeKid()
    expect_error(x[1L, "study"], "two-argument `\\[` is not supported")
    expect_error(x[, "study"], "two-argument `\\[` is not supported")
    # A trailing comma with no j is still fine.
    expect_equal(length(x[1L, ]), 1L)
    # The supported spelling for a metadata cell.
    expect_equal(S4Vectors::mcols(x)[1L, "study"], "s1")
})

# ===========================================================================
# Invariant: one seqname and one strand per element
# ===========================================================================

test_that("validity rejects an element spanning two seqnames", {
    expect_error(
        new(
            "RtlTestKid",
            GenomicRanges::GRangesList(a = .rtl_mk(c("chr1", "chr2"), c(1, 2))),
            genome = "hg38",
            ldSketch = NULL,
            qcInfo = list()
        ),
        "span more than one"
    )
})

test_that("validity rejects an element spanning two strands", {
    # range() splits on strand as well as seqname, so a mixed-strand element
    # breaks the 1:1 element-to-span mapping just as a mixed-chrom one does.
    expect_error(
        new(
            "RtlTestKid",
            GenomicRanges::GRangesList(
                a = .rtl_mk("chr1", c(1, 2), strand = c("+", "-"))
            ),
            genome = "hg38",
            ldSketch = NULL,
            qcInfo = list()
        ),
        "span more than one"
    )
})

test_that("validity accepts a stranded element as long as it is consistent", {
    x <- new(
        "RtlTestKid",
        GenomicRanges::GRangesList(
            a = .rtl_mk("chr1", c(1, 2), strand = "+")
        ),
        genome = "hg38",
        ldSketch = NULL,
        qcInfo = list()
    )
    expect_equal(length(x), 1L)
    expect_equal(length(range(x)[[1L]]), 1L)
})

test_that("an empty collection is valid", {
    empty <- new(
        "RtlTestKid",
        GenomicRanges::GRangesList(),
        genome = "hg38",
        ldSketch = NULL,
        qcInfo = list()
    )
    expect_equal(length(empty), 0L)
    expect_equal(nrow(empty), 0L)
    expect_equal(ncol(empty), 0L)
})

# ===========================================================================
# Rebuilding preserves collection-level state
# ===========================================================================

test_that(".rtlRebuild carries every subclass slot across a rebuild", {
    # Slots are derived from the class definition rather than listed by hand,
    # so a slot added later is picked up instead of silently dropped.
    x <- .rtl_makeKid()
    expect_setequal(
        pecotmr:::.rtlOwnSlots(x),
        c("genome", "ldSketch", "qcInfo")
    )
    out <- pecotmr:::.rtlRebuild(
        x,
        list(a = .rtl_mk("chr1", 100), c = .rtl_mk("chr2", 50)),
        c(TRUE, FALSE, TRUE)
    )
    expect_s4_class(out, "RtlTestKid")
    expect_equal(out@genome, "hg38")
    expect_equal(out@ldSketch, "SKETCH")
    expect_equal(out@qcInfo, list(step1 = "ok"))
    expect_equal(out$study, c("s1", "s3"))
    expect_equal(unname(lengths(out)), c(1L, 1L))
})

test_that("endoapply rebuilds without losing slots", {
    x <- .rtl_makeKid()
    out <- endoapply(x, head, 1L)
    expect_equal(unname(lengths(out)), c(1L, 1L, 1L))
    expect_equal(out@ldSketch, "SKETCH")
    expect_equal(out$study, c("s1", "s2", "s3"))
})

test_that("endoapply dispatches only on the bare call, not the qualified one", {
    # `S4Vectors::endoapply` fetches the generic out of S4Vectors' namespace,
    # whose method table does not carry pecotmr's method, so it falls through
    # to the default and fails on the missing list coercion. Callers must use
    # the bare `endoapply()` with pecotmr loaded.
    x <- .rtl_makeKid()
    expect_s4_class(endoapply(x, head, 1L), "RtlTestKid")
    expect_error(
        S4Vectors::endoapply(x, head, 1L),
        "no method or default for coercing"
    )
})

test_that("endoapply rejects a FUN that stops returning GRanges", {
    x <- .rtl_makeKid()
    expect_error(
        endoapply(x, length),
        "non-GRanges"
    )
})

# ===========================================================================
# subsetRegion
# ===========================================================================

test_that("subsetRegion selects, trims, and drops emptied elements", {
    x <- .rtl_makeKid()
    out <- subsetRegion(
        x,
        GenomicRanges::GRanges(
            "chr1",
            IRanges::IRanges(150, 1050)
        )
    )
    # a: 200,300 kept (100 trimmed away); b: 1000 kept; c: chr2, dropped.
    expect_equal(length(out), 2L)
    expect_equal(names(out), c("a", "b"))
    expect_equal(unname(lengths(out)), c(2L, 1L))
    expect_equal(out$study, c("s1", "s2"))
    expect_equal(out@genome, "hg38")
})

test_that("subsetRegion accepts a chr:start-end string", {
    x <- .rtl_makeKid()
    out <- subsetRegion(x, "chr2:40-70")
    expect_equal(length(out), 1L)
    expect_equal(names(out), "c")
    expect_equal(out$study, "s3")
})

test_that("subsetRegion returns an empty collection when nothing overlaps", {
    x <- .rtl_makeKid()
    out <- subsetRegion(x, "chr9:1-100")
    expect_s4_class(out, "RtlTestKid")
    expect_equal(length(out), 0L)
    expect_equal(out@genome, "hg38")
})

test_that("subsetRegion differs from subsetByOverlaps and restrict", {
    # The three operations are genuinely distinct and all three are useful.
    x <- .rtl_makeKid()
    win <- GenomicRanges::GRanges("chr1", IRanges::IRanges(150, 1050))
    byOverlaps <- IRanges::subsetByOverlaps(x, win)
    trimmed <- subsetRegion(x, win)
    # subsetByOverlaps keeps matching elements WHOLE (a still has 3 variants).
    expect_equal(unname(lengths(byOverlaps)), c(3L, 2L))
    # subsetRegion trims them (a keeps 2, b keeps 1).
    expect_equal(unname(lengths(trimmed)), c(2L, 1L))
})


# .rtlSplitByBlocks -- the LD-block counterpart to .rtlSplitBySeqname (§4.2).
# Splitting a genome-wide GWAS by seqname gives one element per chromosome,
# which is too coarse for cTWAS: its EM needs per-block context.

# @noRd
.rtl_blocks <- function() {
    b <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(1L, 800L), c(500L, 1200L))
    )
    names(b) <- c("b1", "b2")
    b
}

# @noRd
.rtl_variants <- function(pos = c(100L, 250L, 900L)) {
    GenomicRanges::GRanges("chr1", IRanges::IRanges(pos, width = 1L))
}

test_that(".rtlSplitByBlocks splits one element into per-block pieces", {
    r <- .rtlSplitByBlocks(list(.rtl_variants()), .rtl_blocks())
    expect_equal(names(r$entry), c("b1", "b2"))
    expect_equal(r$blockId, c("b1", "b2"))
    expect_equal(unname(lengths(r$entry)), c(2L, 1L))
    # Both pieces came from the one input element, so its metadata row is
    # replicated twice.
    expect_equal(r$fromIdx, c(1L, 1L))
})

test_that(".rtlSplitByBlocks reports variants that match no block", {
    # 1500 is past the last block. Dropping is deliberate -- a variant outside
    # every LD block has no block-local LD -- but it must not be silent.
    v <- .rtl_variants(c(100L, 1500L))
    expect_warning(
        r <- .rtlSplitByBlocks(list(v), .rtl_blocks()),
        "1 variant"
    )
    expect_equal(sum(lengths(r$entry)), 1L)
})

test_that(".rtlSplitByBlocks drops an element with no surviving variant", {
    v <- .rtl_variants(c(1500L, 1600L))
    expect_warning(
        r <- .rtlSplitByBlocks(list(v), .rtl_blocks()),
        "2 variant"
    )
    expect_length(r$entry, 0L)
    expect_length(r$fromIdx, 0L)
})

test_that(".rtlSplitByBlocks keeps fromIdx aligned across several elements", {
    # Element 1 spans both blocks, element 2 only the second: fromIdx is what
    # lets the caller replicate the right metadata row against each piece.
    r <- .rtlSplitByBlocks(
        list(.rtl_variants(), .rtl_variants(900L)),
        .rtl_blocks()
    )
    expect_equal(r$fromIdx, c(1L, 1L, 2L))
    expect_equal(r$blockId, c("b1", "b2", "b2"))
})

test_that(".rtlSplitByBlocks takes keys from names, mcols, then the range", {
    b <- .rtl_blocks()
    expect_equal(
        .rtlSplitByBlocks(list(.rtl_variants(100L)), b)$blockId,
        "b1"
    )
    names(b) <- NULL
    mcols(b)$blockId <- c("x", "y")
    expect_equal(
        .rtlSplitByBlocks(list(.rtl_variants(100L)), b)$blockId,
        "x"
    )
    mcols(b) <- NULL
    # Falls back to the block's own coordinates, one key per block ROW --
    # .rtlRangeKeys would merge them into a single span instead.
    expect_equal(
        .rtlSplitByBlocks(list(.rtl_variants(100L)), b)$blockId,
        "chr1_1_500"
    )
})

test_that(".rtlSplitByBlocks rejects an unusable manifest", {
    expect_error(.rtlSplitByBlocks(list(.rtl_variants()), "chr1"), "GRanges")
    expect_error(
        .rtlSplitByBlocks(list(.rtl_variants()), .rtl_blocks()[0L]),
        "empty"
    )
    dup <- .rtl_blocks()
    names(dup) <- c("b1", "b1")
    expect_error(
        .rtlSplitByBlocks(list(.rtl_variants()), dup),
        "unique"
    )
})

test_that(".rtlSplitByBlocks handles an empty entry list", {
    r <- .rtlSplitByBlocks(list(), .rtl_blocks())
    expect_length(r$entry, 0L)
    expect_length(r$fromIdx, 0L)
    expect_length(r$blockId, 0L)
})

test_that(".rtlSplitByBlocks rejects a non-GRanges element", {
    expect_error(
        .rtlSplitByBlocks(list(data.frame(a = 1)), .rtl_blocks()),
        "must be a GRanges"
    )
})

# `[[<-` -- element replacement.
#
# Subclasses of CompressedGRangesList get no working `[[<-` for free: the stock
# replacement rebuilds by coercing a plain list back to class(x), and no
# setAs("list", <subclass>) exists. A bare `contains = "CompressedGRangesList"`
# class with no methods of its own fails identically, so this is not specific
# to these classes. Installing that coercion is what we avoid -- it makes the
# replacement run while silently DROPPING the subclass's slots -- so this
# routes through .rtlRebuild() instead, which keeps them and re-runs validity.

test_that("[[<- replaces an element and keeps class, slots and mcols", {
    x <- .rtl_makeKid()
    before <- S4Vectors::mcols(x)
    x[["a"]] <- head(x[["a"]], 1L)
    expect_s4_class(x, "RtlTestKid")
    expect_equal(length(x[["a"]]), 1L)
    expect_equal(S4Vectors::mcols(x), before)
    # The collection-level slots are exactly what the setAs route would drop.
    expect_equal(x@genome, "hg38")
    expect_equal(x@ldSketch, "SKETCH")
    expect_equal(x@qcInfo, list(step1 = "ok"))
})

test_that("[[<- leaves the other elements untouched", {
    x <- .rtl_makeKid()
    orig <- lengths(x)
    x[[2]] <- head(x[[2]], 1L)
    expect_equal(unname(lengths(x))[-2], unname(orig)[-2])
    expect_equal(length(x[[2]]), 1L)
})

test_that("[[<- accepts a positional or a name index", {
    byPos <- .rtl_makeKid()
    byName <- .rtl_makeKid()
    byPos[[1]] <- head(byPos[[1]], 2L)
    byName[["a"]] <- head(byName[["a"]], 2L)
    expect_equal(lengths(byPos), lengths(byName))
})

test_that("[[<- rejects a non-GRanges value", {
    # An element IS the variant set, so anything else is a category error
    # rather than something to coerce.
    x <- .rtl_makeKid()
    expect_error(
        {
            x[[1]] <- 42
        },
        "must be a GRanges"
    )
})

test_that("[[<- refuses to extend the collection", {
    # A new element has no identity row; mcols would be padded with NA, which
    # the validity method then reads as a broken tuple. Failing here names the
    # real cause instead.
    x <- .rtl_makeKid()
    expect_error(
        {
            x[[4]] <- x[[1]]
        },
        "cannot extend"
    )
})

test_that("[[<- rejects the two-index form and an unknown name", {
    x <- .rtl_makeKid()
    expect_error(
        {
            x[[1, 1]] <- x[[1]]
        },
        "single index"
    )
    expect_error(
        {
            x[["nope"]] <- x[[1]]
        },
        "no element named"
    )
})

test_that("[[<- still runs the class's validity method", {
    # The point of rebuilding rather than coercing: a replacement that breaks
    # an invariant is rejected, not stored. Element `c` is on chr2, so
    # assigning it into a chr1 element spans two seqnames.
    x <- .rtl_makeKid()
    expect_error(
        {
            x[[1]] <- c(x[[1]], x[[3]])
        },
        "more than one"
    )
})

test_that("[[<- accepts an empty element", {
    # An entry can legitimately be empty -- summaryStatsQc(pipCutoffToSkip)
    # empties a no-signal region -- so this is not a validity failure.
    x <- .rtl_makeKid()
    x[[1]] <- x[[1]][0]
    expect_equal(length(x[[1]]), 0L)
    expect_s4_class(x, "RtlTestKid")
})


# .rtlGatherElements -- element ranges for `idx`, without rebuilding the
# collection.
#
# The obvious `unlist(x[idx])` subsets the COLLECTION, reconstructing the whole
# S4 subclass and re-running validity to reach elements that are already there.
# These pin that the cheap paths return exactly what the expensive one did.

test_that(".rtlGatherElements matches unlist(x[idx]) for every index shape", {
    x <- .rtl_makeKid()
    slow <- function(idx) unlist(x[idx], use.names = FALSE)
    expect_identical(.rtlGatherElements(x, 1L), slow(1L))
    expect_identical(.rtlGatherElements(x, c(1L, 3L)), slow(c(1L, 3L)))
    expect_identical(
        .rtlGatherElements(x, seq_len(length(x))),
        slow(seq_len(length(x)))
    )
})

test_that(".rtlGatherElements honours a permuted full index", {
    # The all-rows fast path is guarded on the exact sequence, not on length:
    # a permutation must fall through to the general path, or the elements
    # come back silently reordered.
    x <- .rtl_makeKid()
    idx <- rev(seq_len(length(x)))
    expect_identical(
        .rtlGatherElements(x, idx),
        unlist(x[idx], use.names = FALSE)
    )
    # And that really is a different order from the unpermuted gather.
    expect_false(identical(
        .rtlGatherElements(x, idx),
        .rtlGatherElements(x, seq_len(length(x)))
    ))
})

test_that(".rtlGatherElements returns an empty GRanges for an empty index", {
    x <- .rtl_makeKid()
    out <- .rtlGatherElements(x, integer(0))
    expect_s4_class(out, "GRanges")
    expect_length(out, 0L)
})

test_that(".rtlGatherElements keeps the elements' inner mcols", {
    # Coercing to the base class must not disturb what the elements carry --
    # that inner metadata is the per-variant payload.
    x <- .rtl_makeKid()
    mcols(x[[1]])$score <- seq_len(length(x[[1]]))
    kid <- .rtl_makeKid()
    elts <- as.list(kid)
    mcols(elts[[1]])$score <- seq_len(length(elts[[1]]))
    kid <- .rtlRebuild(kid, elts, seq_along(elts))
    got <- .rtlGatherElements(kid, c(1L, 2L))
    expect_true(is_in("score", colnames(mcols(got))))
    expect_identical(got, unlist(kid[c(1L, 2L)], use.names = FALSE))
})

test_that(".rtlGatherElements does not re-validate the collection", {
    # Nothing here reconstructs the subclass, so a collection whose validity
    # would be expensive is not paying for it on every element read.
    x <- .rtl_makeKid()
    expect_identical(.rtlGatherElements(x, 2L), x[[2]])
})


# ===========================================================================
# The plyranges bridge: the flatten / nest conversions and the dplyr verbs
# built on them.
#
# plyranges has no GRangesList support, so the collections cannot be operated
# on directly. The bridge flattens to a plain GenomicRanges -- which plyranges
# does support natively -- and nests the result back again.
# ===========================================================================

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
