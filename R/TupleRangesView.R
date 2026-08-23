# =============================================================================
# TupleRangesView -- the plyranges bridge for RangedTupleList collections
#
# plyranges has no GRangesList support at all: its verbs dispatch on
# GenomicRanges / Ranges, and a plain CompressedGRangesList fails exactly as
# these collections do. That is not a deficiency in the collections, and it is
# not fixable by changing them.
#
# It IS fixable through the extension point GenomicRanges provides:
# `DelegatingGenomicRanges` is a virtual class that wraps a flat `delegate` and
# *is* a GenomicRanges, and plyranges ships methods for it. A subclass supplies
# a `[` hook and inherits the verb surface instead of reimplementing it.
#
# The conversion is the right shape rather than a workaround: plyranges' own
# `GroupedGenomicRanges` (slots group_keys / group_indices / n / delegate)
# represents grouping as flat ranges plus group metadata, which is exactly a
# RangedTupleList transposed. Same information, other representation.
# =============================================================================

#' Flat view of a tuple collection for plyranges verbs
#'
#' A \code{\link[GenomicRanges]{GenomicRanges}} view over a
#' \code{RangedTupleList}: every element's ranges concatenated, with the
#' collection's identity tuple broadcast onto each range. Because it is a
#' \code{GenomicRanges}, plyranges verbs act on it directly.
#'
#' Users do not normally build one. The verb methods on \code{RangedTupleList}
#' convert to a view, delegate, and convert back.
#'
#' @slot delegate The flattened \code{GRanges}.
#' @name TupleRangesView-class
#' @keywords internal
#' @export
setClass("TupleRangesView", contains = "DelegatingGenomicRanges")

# Build a view over `flat`.
#
# DelegatingGenomicRanges carries its OWN elementMetadata alongside the
# delegate's, and validity requires it to be parallel to the object. Supplying
# only `delegate=` leaves it zero-row against a longer object and the class
# rejects itself with "'mcols(x)' is not parallel to 'x'", so it is set here.
# @noRd
.tupleRangesView <- function(flat) {
    methods::new(
        "TupleRangesView",
        delegate = flat,
        elementMetadata = S4Vectors::make_zero_col_DFrame(length(flat))
    )
}

# The one hook a DelegatingGenomicRanges subclass has to supply. Without it,
# plyranges' subsetting verbs (slice, filter_by_overlaps, join_overlap_*) fail
# with "subscript is a NSBS object that is incompatible with the current
# subsetting operation".
#' @rdname TupleRangesView-class
#' @export
setMethod("[", "TupleRangesView", function(x, i, j, ..., drop = TRUE) {
    if (!missing(j)) {
        abort("two-argument `[` is not supported on a TupleRangesView.")
    }
    .tupleRangesView(x@delegate[i])
})

#' @rdname show-methods
#' @export
setMethod("show", "TupleRangesView", function(object) {
    cat(glue(
        "TupleRangesView: {length(object@delegate)} ranges, ",
        "{ncol(mcols(object@delegate))} metadata column(s)\n",
        .trim = FALSE
    ))
    invisible(NULL)
})

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
# Patches for the two verbs plyranges' DelegatingGenomicRanges support misses
# -----------------------------------------------------------------------------

# Both work on a plain GRanges and fail on any DelegatingGenomicRanges:
# `select` reports "Cannot select/rename the following columns: seqnames,
# start, end, width, strand", and `group_by` -- for which plyranges defines no
# DelegatingGenomicRanges method at all -- reports "Can't select columns with
# `dots`". Both are upstream gaps, not something about this class.
#
# The fix in each case is to run the verb on the delegate, where plyranges
# works, and re-wrap.

#' @exportS3Method dplyr::select
select.TupleRangesView <- function(.data, ..., .drop_ranges = FALSE) {
    out <- dplyr::select(.data@delegate, ..., .drop_ranges = .drop_ranges)
    if (isTRUE(.drop_ranges)) {
        return(out)
    }
    .tupleRangesView(out)
}

# Not re-wrapped: grouping is plyranges' own representation, and the result is
# a GroupedGenomicRanges rather than a view of a collection.
#' @exportS3Method dplyr::group_by
group_by.TupleRangesView <- function(.data, ..., .add = FALSE) {
    dplyr::group_by(.data@delegate, ...)
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

# Unwrap a view back to its delegate; pass a plain GRanges through.
# @noRd
.rtlAsRanges <- function(out) {
    if (methods::is(out, "TupleRangesView")) out@delegate else out
}

# `...` is forwarded directly rather than captured with enquos() and spliced:
# splicing hands the verb quosure OBJECTS instead of expressions to evaluate,
# and plyranges then fails with "Argument to filter condition must evaluate to
# a logical vector".

#' @exportS3Method dplyr::filter
filter.RangedTupleList <- function(.data, ...) {
    .rtlRequirePlyranges("filter")
    out <- dplyr::filter(flattenTupleRanges(.data), ...)
    nestTupleRanges(.rtlAsRanges(out), .data)
}

#' @exportS3Method dplyr::mutate
mutate.RangedTupleList <- function(.data, ...) {
    .rtlRequirePlyranges("mutate")
    out <- dplyr::mutate(flattenTupleRanges(.data), ...)
    nestTupleRanges(.rtlAsRanges(out), .data)
}

#' @exportS3Method dplyr::select
select.RangedTupleList <- function(.data, ...) {
    .rtlRequirePlyranges("select")
    # The identity columns are kept whatever the selection, or the result
    # could not be nested back into its elements.
    out <- dplyr::select(flattenTupleRanges(.data), ...)
    nestTupleRanges(
        .rtlRestoreKeys(.rtlAsRanges(out), flattenTupleRanges(.data), .data),
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
    nestTupleRanges(.rtlAsRanges(out), .data)
}

#' @exportS3Method dplyr::summarise
summarise.RangedTupleList <- function(.data, ...) {
    .rtlRequirePlyranges("summarise")
    .rtlAsRanges(dplyr::summarise(flattenTupleRanges(.data), ...))
}

# No count() method: plyranges defines none for GRanges either, so there is
# nothing to delegate to. Use summarise(group_by(x, ...), n = n()).

#' @exportS3Method dplyr::group_by
group_by.RangedTupleList <- function(.data, ...) {
    .rtlRequirePlyranges("group_by")
    .rtlAsRanges(dplyr::group_by(flattenTupleRanges(.data), ...))
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
