# =============================================================================
# RangedTupleList
# -----------------------------------------------------------------------------
# Virtual base for the pecotmr collection classes. Each element is a GRanges of
# variants; the identity tuple (study / context / trait / method / ...) and any
# other per-element metadata live in `mcols(x)`, which is itself a DataFrame.
# That gives arbitrary multi-key tables AND the Bioconductor range machinery in
# one object, so `subsetByOverlaps()` / `range()` come for free.
#
# Two API layers, deliberately distinct:
#
#   x$study     the mcols COLUMN named "study"   (the old DFrame feel)
#   x[[i]]      the i-th ELEMENT, a GRanges      (GRangesList semantics)
#
# `[[` is left alone on purpose: on a GRangesList it addresses elements, so
# overloading it for columns would collide. Use `mcols(x)[i, j]` for a cell.
#
# Slots hold COLLECTION-LEVEL scalars only (ldSketch, genome, qcInfo). Anything
# with one value per element belongs in mcols -- a per-element slot silently
# desyncs on subsetting, because `x[i]` narrows the elements and the mcols but
# has no idea the slot is parallel to them.
# =============================================================================

#' @include AllGenerics.R
NULL

#' @title RangedTupleList methods
#' @description Container and metadata methods shared by every
#'   \code{RangedTupleList} subclass.
#'
#'   \code{nrow()} counts elements, \code{ncol()} counts metadata columns and
#'   \code{colnames()} names them, so the collections read exactly as the
#'   \code{DFrame}-shaped ones they replace. \code{$} and \code{$<-}
#'   read and write \code{mcols} columns. \code{[} selects elements and
#'   \strong{rejects} a second argument, because the stock method would
#'   silently discard it. \code{endoapply()} maps over the elements and
#'   rebuilds the object with its collection-level slots intact.
#' @param x,X A \code{RangedTupleList} collection.
#' @param name An \code{mcols} column name.
#' @param do.NULL,prefix Unused; present for generic compatibility.
#' @param value Replacement value for that column.
#' @param i Element selector.
#' @param j Not supported; supplying it is an error.
#' @param drop Passed through to the inherited method.
#' @param ... Additional arguments passed on.
#' @param FUN A function applied to each element, returning a \code{GRanges}.
#' @return \code{nrow()} and \code{ncol()} integers; \code{colnames()} a
#'   character vector; \code{$} a metadata column; \code{$<-}, \code{[} and
#'   \code{endoapply()} a collection of the same class.
#' @examples
#' data(qtlSumStatsMulticontextExample)
#' mc <- qtlSumStatsMulticontextExample
#' # A collection reads as a table of tuples: one row per element, with the
#' # identity in mcols and the variants in the element itself.
#' nrow(mc)
#' ncol(mc)
#' colnames(mc)
#' mc$context
#' length(mc[[1]])
#' # Row subsetting and endoapply() both return the same class, carrying the
#' # collection-level slots with them.
#' class(mc[1:2])
#' lengths(endoapply(mc, head, n = 5L))
#' @name RangedTupleList-methods
#' @rdname RangedTupleList-methods
NULL

#' @title Ranged Tuple List
#' @description Virtual base class for pecotmr's per-element collections. Each
#'   element is a \code{GRanges} of variants; the identity tuple and any other
#'   per-element metadata live in \code{mcols(x)}.
#'
#'   Concrete subclasses add their own collection-level slots. Per-element
#'   values must go in \code{mcols}, never in a slot: subsetting narrows the
#'   elements and \code{mcols} together, but leaves slots untouched, so a
#'   per-element slot silently ends up the wrong length and positionally wrong.
#'
#' @section Invariant: One seqname and one strand per element. \code{range()}
#'   splits on both, so a mixed element would return two ranges and break the
#'   1:1 element-to-span mapping the rest of the package relies on. Variants
#'   carry strand \code{*} unless the source data says otherwise.
#' @export
setClass("RangedTupleList", contains = c("VIRTUAL", "CompressedGRangesList"))

methods::setValidity("RangedTupleList", function(object) {
    .rtlCheckSingleSpan(object)
})

# Each element must occupy exactly one (seqname, strand) span. range() already
# splits on both, so an element that yields more than one range has violated
# the invariant -- no need to inspect seqnames and strand separately.
# @noRd
.rtlCheckSingleSpan <- function(object) {
    if (length(object) == 0L) {
        return(TRUE)
    }
    spans <- lengths(range(object))
    bad <- which(spans > 1L)
    if (length(bad) == 0L) {
        return(TRUE)
    }
    labels <- names(object)[bad] %||% as.character(bad)
    glue(
        "element(s) {str_flatten(labels, ', ')} span more than one ",
        "(seqname, strand) combination; each element must cover exactly one. ",
        "Split by seqname at construction, and set strand to '*' unless the ",
        "source data is stranded."
    )
}

# -----------------------------------------------------------------------------
# Rebuilding
# -----------------------------------------------------------------------------

# The slots a concrete subclass adds on top of the container itself. Used to
# carry collection-level state across a rebuild; deriving it from the class
# definition means a slot added later is picked up automatically instead of
# being silently dropped.
# @noRd
.rtlOwnSlots <- function(x) {
    setdiff(
        names(methods::getSlots(class(x))),
        names(methods::getSlots("CompressedGRangesList"))
    )
}

# Rebuild `x` from a new list of GRanges elements, carrying every
# collection-level slot over and narrowing mcols by `keep`.
#
# This is the ONLY supported way to produce a modified collection. `endoapply`
# via a `setAs("list", <class>)` coercion looks like it works but silently
# drops the subclass's slots -- verified -- which is exactly the desync this
# class exists to prevent.
# @noRd
# Put `x`'s seqinfo back onto a rebuilt GRangesList. Only the seqlevels the
# rebuild actually spans can be kept, so this restores the genome (and any
# lengths) for those, and supplies the whole seqinfo when the rebuild spans
# nothing.
# @noRd
.rtlRestoreSeqinfo <- function(grl, x) {
    si <- GenomeInfoDb::seqinfo(x)
    if (length(GenomeInfoDb::seqlevels(si)) == 0L) {
        return(grl)
    }
    if (length(GenomeInfoDb::seqlevels(grl)) == 0L) {
        GenomeInfoDb::seqinfo(grl) <- si
        return(grl)
    }
    keepLevels <- intersect(
        GenomeInfoDb::seqlevels(si),
        GenomeInfoDb::seqlevels(grl)
    )
    GenomeInfoDb::seqinfo(
        grl,
        new2old = match(keepLevels, GenomeInfoDb::seqlevels(grl))
    ) <- si[keepLevels]
    grl
}

.rtlRebuild <- function(x, elements, keep) {
    # mcols and slots are attached BEFORE new(), not after: new() validates
    # during initialize(), and a subclass's validity method reads its identity
    # columns and slots. Building the object bare and filling it in afterwards
    # trips that check on the way past.
    grl <- GenomicRanges::GRangesList(elements)
    # seqinfo is collection-level state, exactly like the slots below: it
    # carries the genome build. A rebuild from bare elements starts with the
    # seqlevels those elements happen to span -- none at all when everything
    # was dropped -- so the original seqinfo is merged back in, or subsetting
    # to nothing would silently discard the build.
    grl <- .rtlRestoreSeqinfo(grl, x)
    md <- mcols(x, use.names = FALSE)
    if (!is.null(md)) {
        mcols(grl) <- md[keep, , drop = FALSE]
    }
    ownSlots <- .rtlOwnSlots(x)
    slotArgs <- set_names(map(ownSlots, .rtlGetSlot, x = x), ownSlots)
    exec(methods::new, class(x), grl, !!!slotArgs)
}

# @noRd
.rtlGetSlot <- function(nm, x) {
    methods::slot(x, nm)
}

# The element ranges for `idx`, concatenated, without rebuilding the
# collection.
#
# The obvious `unlist(x[idx])` subsets the COLLECTION, which reconstructs the
# whole S4 subclass and re-runs its validity method just to reach elements that
# are already there. Measured on an 11-row TwasWeights: 14.9 ms for a single
# index against 0.9 ms for `x[[i]]`, and 16.2 ms for every row against 0.2 ms
# for a bare `unlist(x)`. Validity is not the cost (0.1 ms) -- the rebuild is.
#
# Nothing here re-validates: the elements came out of a collection that was
# validated when it was built.
# @noRd
.rtlGatherElements <- function(x, idx) {
    if (length(idx) == 0L) {
        return(unlist(x, use.names = FALSE)[0L])
    }
    if (length(idx) == 1L) {
        return(x[[idx[[1L]]]])
    }
    # Guarded on the exact sequence, not on length: a PERMUTED full index must
    # fall through to the general path, which honours the caller's order, or
    # the rows would come back silently reordered.
    if (identical(as.integer(idx), seq_len(length(x)))) {
        return(unlist(x, use.names = FALSE))
    }
    # Coercing to the base class drops the subclass rebuild and its validity
    # pass while leaving the elements (and their inner mcols) untouched.
    unlist(
        as(x, "CompressedGRangesList")[idx],
        use.names = FALSE
    )
}

# -----------------------------------------------------------------------------
# API shims: make the mcols table feel like the DFrame it replaced
# -----------------------------------------------------------------------------

#' @rdname RangedTupleList-methods
#' @importFrom BiocGenerics nrow ncol
#' @export
setMethod("nrow", "RangedTupleList", function(x) length(x))

#' @rdname RangedTupleList-methods
#' @export
setMethod("ncol", "RangedTupleList", function(x) {
    md <- mcols(x, use.names = FALSE)
    if (is.null(md)) 0L else ncol(md)
})

#' @rdname RangedTupleList-methods
#' @importFrom BiocGenerics colnames
#' @export
setMethod("colnames", "RangedTupleList", function(x, do.NULL = TRUE, prefix) {
    # Completes the nrow / ncol / `$` set. Without it `colnames(x)` returns
    # NULL rather than erroring, so code that asks a collection for its column
    # names gets a silent wrong answer instead of a signal.
    md <- mcols(x, use.names = FALSE)
    if (is.null(md)) NULL else colnames(md)
})

#' @rdname RangedTupleList-methods
#' @export
setMethod("$", "RangedTupleList", function(x, name) {
    md <- mcols(x, use.names = FALSE)
    if (is.null(md)) {
        return(NULL)
    }
    md[[name]]
})

#' @rdname RangedTupleList-methods
#' @export
setMethod("$<-", "RangedTupleList", function(x, name, value) {
    md <- mcols(x, use.names = FALSE)
    if (is.null(md)) {
        md <- S4Vectors::DataFrame(set_names(list(value), name))
    } else {
        md[[name]] <- value
    }
    mcols(x) <- md
    x
})

#' @rdname RangedTupleList-methods
#' @export
setMethod("[", "RangedTupleList", function(x, i, j, ..., drop = TRUE) {
    # A two-argument `[` is rejected rather than honoured or ignored. The
    # stock CompressedList method silently DISCARDS `j` and returns whole
    # elements, so `x[1, "study"]` would quietly change meaning from "that
    # cell" to "that element" -- a silent flip is worse than an error.
    if (!missing(j)) {
        msg <- glue(
            "two-argument `[` is not supported on a {class(x)}; use ",
            "mcols(x)[i, j] for a metadata cell, or x[i] for elements."
        )
        abort(msg)
    }
    methods::callNextMethod(x, i, ..., drop = drop)
})

#' @rdname RangedTupleList-methods
#' @export
setMethod("[[<-", "RangedTupleList", function(x, i, j, ..., value) {
    # Subclasses of CompressedGRangesList get no working `[[<-`: the stock
    # replacement rebuilds by coercing a plain list back to class(x), and no
    # `setAs("list", <subclass>)` exists, so it fails with an opaque "failed to
    # coerce 'list(value)'" (a bare `contains = "CompressedGRangesList"` class
    # with no methods fails identically -- this is not specific to pecotmr).
    #
    # The coercion is deliberately NOT installed: adding it makes the
    # replacement run while silently dropping the subclass's slots. Rebuilding
    # through .rtlRebuild() instead keeps every slot and re-runs validity, so a
    # replacement that would break the identity invariants is rejected rather
    # than stored.
    if (!missing(j)) {
        abort(glue("`[[<-` on a {class(x)} takes a single index."))
    }
    if (!methods::is(value, "GRanges")) {
        msg <- glue(
            "`[[<-` value must be a GRanges; a {class(x)} element is the ",
            "variant set itself (got {class(value)[[1L]]})."
        )
        abort(msg)
    }
    idx <- .rtlAssignIndex(x, i)
    elements <- as.list(x)
    elements[[idx]] <- value
    # Growing the collection has no defined identity row for the new element,
    # and mcols would be padded with NA -- which the validity method reads as a
    # broken tuple. Rejected here so the failure names the real cause.
    if (idx > length(x)) {
        msg <- glue(
            "`[[<-` cannot extend a {class(x)}: a new element has no identity ",
            "row. Build the collection with its full tuple set, or use c()."
        )
        abort(msg)
    }
    .rtlRebuild(x, elements, seq_along(x))
})

# Resolve `[[` index forms (positive integer or element name) to a position.
# @noRd
.rtlAssignIndex <- function(x, i) {
    if (is.character(i)) {
        pos <- match(i, names(x))
        if (is.na(pos)) {
            abort(glue("`[[<-`: no element named '{i}'."))
        }
        return(pos)
    }
    if (length(i) != 1L || is.na(i)) {
        abort("`[[<-` takes a single non-NA index.")
    }
    as.integer(i)
}

#' @rdname RangedTupleList-methods
#' @importFrom S4Vectors endoapply
#' @export
setMethod("endoapply", "RangedTupleList", function(X, FUN, ...) {
    # Defined explicitly rather than left to the inherited machinery. Stock
    # `endoapply` on a subclass errors unless a `setAs("list", <class>)`
    # coercion exists, and adding that coercion makes it run while silently
    # DROPPING the subclass's slots (verified: genome / ldSketch come back
    # unset). Rebuilding here keeps every collection-level slot.
    #
    # Call this as `endoapply(x, ...)`, NOT `S4Vectors::endoapply(x, ...)`:
    # the qualified form fetches the generic from S4Vectors' namespace, whose
    # method table does not carry this method, and silently falls through to
    # the default.
    elements <- map(as.list(X), FUN, ...)
    notRanges <- !map_lgl(elements, methods::is, "GRanges")
    if (any(notRanges)) {
        msg <- glue(
            "endoapply(): FUN returned a non-GRanges for ",
            "{sum(notRanges)} of {length(elements)} element(s); a ",
            "{class(X)} can only hold GRanges."
        )
        abort(msg)
    }
    .rtlRebuild(X, elements, seq_along(X))
})

# A stable per-element key for the element's own genomic span, used where an
# identity tuple needs a positional component. Derived from range() rather than
# from a stored label so it cannot drift out of step with the variants, and so
# it stays correct after subsetRegion() narrows an element.
# @noRd
.rtlRangeKeys <- function(x) {
    if (length(x) == 0L) {
        return(character(0))
    }
    spans <- range(x)
    map_chr(as.list(spans), .rtlOneRangeKey)
}

# @noRd
.rtlOneRangeKey <- function(g) {
    if (length(g) == 0L) {
        return(NA_character_)
    }
    str_c(
        as.character(seqnames(g))[[1L]],
        start(g)[[1L]],
        end(g)[[1L]],
        sep = "_"
    )
}

# -----------------------------------------------------------------------------
# subsetRegion
# -----------------------------------------------------------------------------

#' @rdname subsetRegion
#' @export
setMethod("subsetRegion", "RangedTupleList", function(x, region, ...) {
    win <- .asGRegion(region)
    elements <- map(as.list(x), .rtlRestrictOne, win = win)
    keep <- lengths(elements) > 0L
    .rtlRebuild(x, elements[keep], keep)
})

# -----------------------------------------------------------------------------
# Construction helper: enforce one seqname per element
# -----------------------------------------------------------------------------

# Split any multi-seqname element into one element per seqname, and report the
# index each piece came from so the caller can replicate that element's
# metadata row alongside it.
#
# A genome-wide GWAS arrives as a single GRanges spanning every chromosome,
# which the one-seqname invariant forbids; splitting is how it becomes a
# per-chromosome collection. Elements already on one seqname pass through
# untouched, so this is a no-op for the QTL side (a cis window is
# single-chromosome by construction).
# @noRd
.rtlSplitBySeqname <- function(entry) {
    entry <- as.list(entry)
    if (length(entry) == 0L) {
        return(list(entry = list(), fromIdx = integer(0)))
    }
    # Checked here rather than in class validity: once the elements are the
    # container, a non-GRanges cannot even be stored, so the error has to come
    # before construction or it surfaces as an unrelated seqnames() failure.
    notRanges <- !map_lgl(entry, methods::is, "GRanges")
    if (any(notRanges)) {
        msg <- glue(
            "every `entry` element must be a GRanges; element(s) ",
            "{str_flatten(which(notRanges), ', ')} are not."
        )
        abort(msg)
    }
    pieces <- map(entry, .rtlSplitOne)
    list(
        entry = unlist(pieces, recursive = FALSE, use.names = TRUE),
        fromIdx = rep(seq_along(entry), lengths(pieces))
    )
}

# One element -> a named list of single-seqname GRanges, in the order the
# seqnames first appear so a caller's ordering is not silently permuted.
# @noRd
.rtlSplitOne <- function(g) {
    chroms <- as.character(seqnames(g))
    levels <- unique(chroms)
    if (length(levels) <= 1L) {
        return(list(g))
    }
    set_names(map(levels, .rtlPickSeqname, g = g, chroms = chroms), levels)
}

# @noRd
.rtlPickSeqname <- function(level, g, chroms) {
    g[chroms == level]
}

# Split each element by an LD-block manifest instead of by seqname, returning
# the same (entry, fromIdx) shape plus the block key each piece belongs to.
#
# A genome-wide GWAS split by seqname gives one element per chromosome, which
# is too coarse for cTWAS: its EM needs per-block context (§4.7). `blocks` is a
# GRanges of LD blocks whose names supply the keys.
#
# A variant matching no block is dropped rather than pooled into a catch-all:
# a variant outside every LD block has no block-local LD to be fine-mapped
# against, so carrying it would put variants in a block whose LD does not
# describe them. The count is reported so the loss is never silent.
# @noRd
.rtlSplitByBlocks <- function(entry, blocks) {
    entry <- as.list(entry)
    if (length(entry) == 0L) {
        return(list(
            entry = list(),
            fromIdx = integer(0),
            blockId = character(0)
        ))
    }
    notRanges <- !map_lgl(entry, methods::is, "GRanges")
    if (any(notRanges)) {
        msg <- glue(
            "every `entry` element must be a GRanges; element(s) ",
            "{str_flatten(which(notRanges), ', ')} are not."
        )
        abort(msg)
    }
    blocks <- .rtlValidateBlocks(blocks)
    pieces <- map(entry, .rtlSplitOneByBlocks, blocks = blocks)
    dropped <- sum(lengths(entry)) - sum(map_int(pieces, .rtlPieceSize))
    if (dropped > 0L) {
        msg <- glue(
            "{dropped} variant(s) overlap no LD block and were dropped; ",
            "a variant outside every block has no block-local LD to be ",
            "fine-mapped against."
        )
        warn(msg)
    }
    list(
        entry = unlist(pieces, recursive = FALSE, use.names = TRUE),
        fromIdx = rep(seq_along(entry), lengths(pieces)),
        blockId = unname(list_c(map(pieces, names)))
    )
}

# The manifest must name its blocks: the names become `blockId`, which is what
# downstream code keys regions by.
# @noRd
.rtlValidateBlocks <- function(blocks) {
    if (!methods::is(blocks, "GRanges")) {
        msg <- glue(
            "`ldBlocks` must be a GRanges of LD blocks ",
            "(got {class(blocks)[[1L]]})."
        )
        abort(msg)
    }
    if (length(blocks) == 0L) {
        abort("`ldBlocks` is empty; it must describe at least one LD block.")
    }
    keys <- .rtlBlockKeys(blocks)
    if (anyDuplicated(keys) > 0L) {
        dup <- unique(keys[duplicated(keys)])
        msg <- glue(
            "`ldBlocks` keys must be unique; repeated: ",
            "{str_flatten(head(dup, 5L), ', ')}."
        )
        abort(msg)
    }
    names(blocks) <- keys
    blocks
}

# Block keys come from names(), a `blockId` mcol, or the rendered range, in
# that order -- so a manifest read straight from a BED-like file still keys
# stably without the caller having to name it.
# @noRd
.rtlBlockKeys <- function(blocks) {
    if (!is.null(names(blocks)) && !any(is.na(names(blocks)))) {
        return(as.character(names(blocks)))
    }
    md <- mcols(blocks)
    if (!is.null(md) && is_in("blockId", colnames(md))) {
        return(as.character(md$blockId))
    }
    # Not .rtlRangeKeys(): that merges via range() to give one key per
    # collection ELEMENT, whereas a manifest needs one key per block row.
    str_c(
        as.character(seqnames(blocks)),
        start(blocks),
        end(blocks),
        sep = "_"
    )
}

# One element -> a named list of per-block GRanges, in manifest order.
# @noRd
.rtlSplitOneByBlocks <- function(g, blocks) {
    hits <- IRanges::findOverlaps(g, blocks, select = "first")
    if (!any(!is.na(hits))) {
        return(list())
    }
    idx <- sort(unique(hits[!is.na(hits)]))
    set_names(
        map(idx, .rtlPickBlock, g = g, hits = hits),
        names(blocks)[idx]
    )
}

# Variants kept across one element's per-block pieces.
# @noRd
.rtlPieceSize <- function(pieces) {
    as.integer(sum(lengths(pieces)))
}

# @noRd
.rtlPickBlock <- function(i, g, hits) {
    g[!is.na(hits) & hits == i]
}

# Keep only the variants of one element that overlap the window.
#
# The seqname is checked first because an element on another chromosome cannot
# overlap by definition, and going straight to overlapsAny() would emit
# GenomicRanges' "no sequence levels in common" warning for every such element
# on every call -- pure noise for a collection spanning several chromosomes.
# @noRd
.rtlRestrictOne <- function(g, win) {
    onWindowChrom <- as.character(seqnames(g)) %in% as.character(seqnames(win))
    if (!any(onWindowChrom)) {
        return(g[0L])
    }
    g[onWindowChrom & IRanges::overlapsAny(g, win)]
}


# =============================================================================
# Combining collections
# =============================================================================

#' @rdname RangedTupleList-methods
#' @importFrom S4Vectors bindROWS
#' @export
setMethod(
    "bindROWS",
    "RangedTupleList",
    function(
        x,
        objects = list(),
        use.names = TRUE,
        ignore.mcols = FALSE,
        check = TRUE
    ) {
        # bindROWS() is the one hook `c()` and `append()` share, so defining it
        # here fixes both. The inherited CompressedList version rbinds the mcols,
        # which requires identical columns and so fails whenever the parts came
        # from runs with different optional columns (traitPos / jointContexts /
        # blockId); it also carries the FIRST part's collection-level slots
        # silently, making the result depend on argument order.
        .combineTupleCollections(compact(c(list(x), objects)), NULL, "c")
    }
)


# =============================================================================
# The plyranges / dplyr bridge for RangedTupleList collections
#
# plyranges has no GRangesList support at all: its verbs dispatch on
# GenomicRanges / Ranges, and a plain CompressedGRangesList fails exactly as
# these collections do. That is not a deficiency in the collections, and it is
# not fixable by changing them.
#
# The bridge is a round trip through a flat GRanges instead: flatten the
# collection (broadcasting its identity tuple onto every range), let plyranges
# work on the shape it already supports, and nest the result back into the
# collection's elements. Nothing about the collections has to change, and no
# wrapper class is needed -- plyranges acts on a real GRanges natively.
# =============================================================================

# -----------------------------------------------------------------------------
# Conversions
# -----------------------------------------------------------------------------

#' Flatten a tuple collection to a plain GRanges
#'
#' Concatenates every element's ranges and broadcasts the collection's identity
#' tuple onto each range, so a plyranges predicate can mix tuple and per-range
#' columns: \code{filter(x, .context == "blood" & pip > 0.9)}.
#'
#' Broadcast columns are prefixed with a dot -- \code{.study},
#' \code{.context}, \code{.trait}, \code{.method} -- following the
#' tidySummarizedExperiment convention for framework-injected columns
#' (\code{.sample}, \code{.feature}). The prefix is not cosmetic: a
#' collection's identity column can collide with a per-range column of the same
#' name carrying different information. On
#' \code{\link{gwasFineMappingExample}} the outer \code{method} is
#' \code{"susie"} while the per-range \code{method} is \code{"susieRss"} --
#' the collection's label against the fitter actually used. Prefixing keeps
#' both, and keeps the name stable so a predicate does not change meaning
#' between objects.
#'
#' Non-atomic \code{mcols} columns -- the per-element \code{susieFit} /
#' \code{cvResult} payloads -- are dropped. They describe an element, not a
#' range, so there is no row to broadcast them onto. Use the accessors
#' (\code{\link{getSusieFit}}, \code{\link{getCvResult}}) for those.
#'
#' @param x A \code{RangedTupleList}.
#' @return A \code{GRanges}.
#' @examples
#' data(qtlFineMappingExample)
#' flat <- flattenTupleRanges(qtlFineMappingExample)
#' head(names(S4Vectors::mcols(flat)))
#' @export
flattenTupleRanges <- function(x) {
    if (!methods::is(x, "RangedTupleList")) {
        msg <- glue(
            "`x` must be a RangedTupleList (got {class(x)[[1L]]})."
        )
        abort(msg)
    }
    gr <- .rtlGatherElements(x, seq_len(length(x)))
    md <- mcols(x, use.names = FALSE)
    if (is.null(md) || ncol(md) == 0L || length(gr) == 0L) {
        return(gr)
    }
    tupleCols <- names(md)[map_lgl(as.list(md), is.atomic)]
    reps <- rep(seq_len(length(x)), lengths(x))
    for (nm in tupleCols) {
        mcols(gr)[[.rtlDotName(nm)]] <- md[[nm]][reps]
    }
    gr
}

# The broadcast name for an identity column. Dotted so it cannot collide with
# a per-range column, and so the name is the same whatever the object holds.
# @noRd
.rtlDotName <- function(nm) {
    str_c(".", nm)
}

#' Re-nest a flattened GRanges into a tuple collection
#'
#' The inverse of \code{\link{flattenTupleRanges}}: ranges are grouped by the
#' identity tuple they carry and returned to \code{template}'s elements, so the
#' collection's class, its per-element payloads and its collection-level slots
#' all survive a plyranges round trip.
#'
#' A tuple with no surviving ranges becomes an empty element rather than being
#' dropped, so the collection keeps its shape and its metadata stays aligned.
#'
#' @param flat A \code{GRanges} produced by \code{flattenTupleRanges} (possibly
#'   filtered or mutated).
#' @param template The collection it came from, supplying the tuple grid, the
#'   payload columns and the slots.
#' @return An object of \code{template}'s class.
#' @examples
#' data(qtlFineMappingExample)
#' flat <- flattenTupleRanges(qtlFineMappingExample)
#' nestTupleRanges(flat, qtlFineMappingExample)
#' @export
nestTupleRanges <- function(flat, template) {
    if (!methods::is(template, "RangedTupleList")) {
        msg <- glue(
            "`template` must be a RangedTupleList (got ",
            "{class(template)[[1L]]})."
        )
        abort(msg)
    }
    if (!methods::is(flat, "GRanges")) {
        msg <- glue("`flat` must be a GRanges (got {class(flat)[[1L]]}).")
        abort(msg)
    }
    keyCols <- .rtlTupleKeyCols(template)
    dotted <- map_chr(keyCols, .rtlDotName)
    wanted <- .rtlTupleKeys(
        mcols(template, use.names = FALSE),
        keyCols,
        n = length(template)
    )
    have <- .rtlTupleKeys(
        mcols(flat, use.names = FALSE),
        dotted,
        n = length(flat)
    )
    elements <- map(wanted, .rtlPickByKey, flat = flat, have = have)
    .rtlRebuild(template, elements, seq_len(length(template)))
}

# The identity columns to group by: the atomic mcols the flattener broadcasts.
# @noRd
.rtlTupleKeyCols <- function(x) {
    md <- mcols(x, use.names = FALSE)
    if (is.null(md)) {
        return(character(0))
    }
    names(md)[map_lgl(as.list(md), is.atomic)]
}

# One key string per row. With no identity columns every range belongs to the
# single element, which is what a one-row collection means.
# @noRd
.rtlTupleKeys <- function(md, keyCols, n) {
    present <- intersect(keyCols, colnames(md))
    if (length(present) == 0L) {
        return(rep("", n))
    }
    exec(str_c, !!!map(present, .rtlKeyPart, md = md), sep = "\r")
}

# NA is mapped to a sentinel rather than left alone: str_c() propagates NA, so
# a single NA-valued identity column (varY is NA_real_ on a z-score collection)
# would turn every key into NA and the subsequent `have == key` into a logical
# subscript full of NAs.
# @noRd
.rtlKeyPart <- function(nm, md) {
    v <- as.character(md[[nm]])
    if_else(is.na(v), "\u0001NA", v)
}

# @noRd
.rtlPickByKey <- function(key, flat, have) {
    flat[have == key]
}


# -----------------------------------------------------------------------------
# Verbs on the collections themselves
# -----------------------------------------------------------------------------

# S3 dispatch reaches a method registered on an S4 VIRTUAL base, so one set of
# methods here serves every RangedTupleList subclass -- the fine-mapping,
# sumstats, TWAS-weight and coloc collections alike.
#
# Shape-preserving verbs flatten, delegate and re-nest, so the caller gets the
# collection back. Reducing verbs return what the verb produces, because a
# summary has no per-element shape to nest into.

# These verbs delegate to plyranges' GRanges methods, which only exist once
# plyranges' namespace is loaded. Without this the failure is an opaque "no
# applicable method for 'filter' applied to an object of class GRanges" --
# pointing at the flattened form rather than at the missing package.
# requireNamespace() both checks and loads, so the check is also the fix.
# @noRd
.rtlRequirePlyranges <- function(verb) {
    if (!requireNamespace("plyranges", quietly = TRUE)) {
        msg <- glue(
            "`{verb}()` on a tuple collection needs the plyranges package; ",
            "install it, or work on flattenTupleRanges(x) directly."
        )
        abort(msg)
    }
    invisible(NULL)
}

# `...` is forwarded directly rather than captured with enquos() and spliced:
# splicing hands the verb quosure OBJECTS instead of expressions to evaluate,
# and plyranges then fails with "Argument to filter condition must evaluate to
# a logical vector".

#' @exportS3Method dplyr::filter
filter.RangedTupleList <- function(.data, ...) {
    .rtlRequirePlyranges("filter")
    out <- dplyr::filter(flattenTupleRanges(.data), ...)
    nestTupleRanges(out, .data)
}

#' @exportS3Method dplyr::mutate
mutate.RangedTupleList <- function(.data, ...) {
    .rtlRequirePlyranges("mutate")
    out <- dplyr::mutate(flattenTupleRanges(.data), ...)
    nestTupleRanges(out, .data)
}

#' @exportS3Method dplyr::select
select.RangedTupleList <- function(.data, ...) {
    .rtlRequirePlyranges("select")
    # The identity columns are kept whatever the selection, or the result
    # could not be nested back into its elements.
    out <- dplyr::select(flattenTupleRanges(.data), ...)
    nestTupleRanges(
        .rtlRestoreKeys(out, flattenTupleRanges(.data), .data),
        .data
    )
}

# Takes n ranges FROM EACH element, matching arrange()'s within-element rule.
#
# Unlike the other verbs this does not flatten: slicing the flat set would take
# n ranges in TOTAL, which for a many-element collection empties all but the
# first. Applying the slice per element is both the right semantics and simpler
# than reconstructing the grouping after a flatten.
#' @exportS3Method dplyr::slice
slice.RangedTupleList <- function(.data, ...) {
    .rtlRequirePlyranges("slice")
    elements <- map(as.list(.data), .rtlSliceOne, ...)
    .rtlRebuild(.data, elements, seq_len(length(.data)))
}

# @noRd
.rtlSliceOne <- function(g, ...) {
    dplyr::slice(g, ...)
}

# Orders WITHIN each element. A global sort of the flattened set followed by
# partitioning on the identity tuple leaves each element internally sorted --
# the two are the same thing for the partitioned result -- so no grouping is
# needed here. Ordering ACROSS elements would permute the elements themselves,
# which the identity tuple cannot express.
#' @exportS3Method dplyr::arrange
arrange.RangedTupleList <- function(.data, ...) {
    .rtlRequirePlyranges("arrange")
    out <- dplyr::arrange(flattenTupleRanges(.data), ...)
    nestTupleRanges(out, .data)
}

#' @exportS3Method dplyr::summarise
summarise.RangedTupleList <- function(.data, ...) {
    .rtlRequirePlyranges("summarise")
    dplyr::summarise(flattenTupleRanges(.data), ...)
}

# No count() method: plyranges defines none for GRanges either, so there is
# nothing to delegate to. Use summarise(group_by(x, ...), n = n()).

#' @exportS3Method dplyr::group_by
group_by.RangedTupleList <- function(.data, ...) {
    .rtlRequirePlyranges("group_by")
    dplyr::group_by(flattenTupleRanges(.data), ...)
}

# Put back any identity column a selection dropped. Without them the ranges
# carry no tuple and every one would fall into the first element.
# @noRd
.rtlRestoreKeys <- function(out, flat, template) {
    dotted <- map_chr(.rtlTupleKeyCols(template), .rtlDotName)
    missing <- setdiff(
        intersect(dotted, colnames(mcols(flat))),
        colnames(mcols(out))
    )
    for (nm in missing) {
        mcols(out)[[nm]] <- mcols(flat)[[nm]]
    }
    out
}
