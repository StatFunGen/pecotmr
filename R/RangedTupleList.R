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
# Slots hold COLLECTION-LEVEL scalars only (ldReference, genome, qcOptions).
# Anything with one value per element belongs in mcols -- a per-element slot
# silently desyncs on subsetting, because `x[i]` narrows the elements and the
# mcols but has no idea the slot is parallel to them.
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
#'   read and write \code{mcols} columns. Storing metadata, whether by
#'   \code{$<-}, \code{mcols<-} or its older spellings \code{values<-} and
#'   \code{elementMetadata<-}, refuses an identity column (\code{study},
#'   \code{context}, \code{trait}, \code{method}, their \code{gwas}
#'   counterparts and \code{blockId}) that is not character; \code{NA} values
#'   are allowed. \code{[} selects elements and
#'   \strong{rejects} a second argument, because the stock method would
#'   silently discard it. \code{endoapply()} maps over the elements and
#'   rebuilds the object with its collection-level slots intact;
#'   \code{sort()}, \code{unique()} and \code{trim()} work on each element
#'   the same way, so they too return the collection's own class rather than
#'   a bare \code{GRangesList}.
#' @param x,X A \code{RangedTupleList} collection.
#' @param name An \code{mcols} column name.
#' @param do.NULL,prefix Unused; present for generic compatibility.
#' @param value For \code{$<-}, the replacement column; for
#'   \code{elementMetadata<-}, the replacement metadata table.
#' @param i Element selector.
#' @param j Not supported; supplying it is an error.
#' @param drop Passed through to the inherited method.
#' @param ... Additional arguments passed on.
#' @param FUN A function applied to each element, returning a \code{GRanges}.
#' @param objects Further collections to append, as \code{c()} and
#'   \code{append()} pass them.
#' @param use.names For \code{trim()}, passed to each element's
#'   \code{trim()}; otherwise unused, present for generic compatibility.
#' @param ignore.mcols,check Unused; present for generic compatibility.
#' @param decreasing Passed to each element's \code{sort()}.
#' @param incomparables Passed to each element's \code{unique()}.
#' @return \code{nrow()} and \code{ncol()} integers; \code{colnames()} a
#'   character vector; \code{$} a metadata column; \code{$<-},
#'   \code{elementMetadata<-}, \code{[}, \code{endoapply()}, \code{sort()},
#'   \code{unique()} and \code{trim()} a collection of the same class.
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
#' @section Sequence names and genome build: Every range pecotmr holds is on
#'   UCSC sequence names (\code{chr1}, \code{chrX}, \code{chrM}), the style
#'   \code{seqlevelsStyle(x) <- "UCSC"} gives. Ranges read from a file,
#'   placed by variant id, passed to a constructor or given as a range
#'   argument are renamed on the way in (\code{"1"} and \code{"chr1"} both
#'   become \code{chr1}), so ranges from different inputs compare as they
#'   are. A query handed straight to an inherited Bioconductor method, such
#'   as \code{subsetByOverlaps()}, is not renamed: give it in UCSC style.
#'   A pipeline's result carries the genome build its inputs name in
#'   \code{seqinfo()}, as do the ranged views of it (\code{topLoci()},
#'   \code{colocVariants()}, ...); inputs naming different builds are
#'   refused. When no input names one, \code{genome()} is \code{NA} until it
#'   is set with \code{genome(x) <- "hg38"}.
#'
#' @section Invariant: One seqname and one strand per element. \code{range()}
#'   splits on both, so a mixed element would return two ranges and break the
#'   1:1 element-to-span mapping the rest of the package relies on. Variants
#'   carry strand \code{*} unless the source data says otherwise.
#' @export
setClass("RangedTupleList", contains = c("VIRTUAL", "CompressedGRangesList"))

methods::setValidity("RangedTupleList", function(object) {
    span <- .rtlCheckSingleSpan(object)
    problems <- c(
        if (!isTRUE(span)) span,
        .rtlIdentityProblems(mcols(object, use.names = FALSE))
    )
    if (length(problems) == 0L) TRUE else problems
})

# The identity columns of the metadata `md` name a row, so each is character.
# NA stays allowed: a GWAS side has no context or trait, and a row off the
# LD-block axis has no blockId.
# @noRd
.rtlIdentityProblems <- function(md) {
    cols <- intersect(
        c(
            "study",
            "context",
            "trait",
            "method",
            "gwasStudy",
            "gwasContext",
            "gwasTrait",
            "gwasMethod",
            "blockId"
        ),
        colnames(md)
    )
    list_c(map(cols, .rtlIdentityProblem, md = md), ptype = character())
}

# @noRd
.rtlIdentityProblem <- function(cn, md) {
    v <- md[[cn]]
    if (is.character(v)) {
        return(character(0))
    }
    as.character(glue(
        "the `{cn}` column must be character (got {class(v)[[1L]]})."
    ))
}

#' @rdname genomicRegion
#' @examples
#' data(qtlFineMappingExample)
#' genomicRegion(qtlFineMappingExample)
#' data(gwasSumStatsS4Example)
#' genomicRegion(gwasSumStatsS4Example)
#' @export
setMethod("genomicRegion", "RangedTupleList", function(x) {
    .getRegionColumn(x)
})

# `obj` carrying the caller's `metadata`: a free-form list, as metadata()
# is on any Bioconductor vector.
# @noRd
#' @importFrom checkmate assertList
.withMetadata <- function(obj, metadata) {
    assertList(metadata, .var.name = "metadata")
    if (length(metadata) == 0L) {
        return(obj)
    }
    S4Vectors::`metadata<-`(obj, value = metadata)
}

# Each element must occupy exactly one (seqname, strand) span. range() already
# splits on both, so an element that yields more than one range has violated
# the invariant -- no need to inspect seqnames and strand separately.
# @noRd
.rtlCheckSingleSpan <- function(object) {
    if (length(object) == 0L) {
        return(TRUE)
    }
    # Compared as integer (seqname, strand) codes against each element's first
    # range: range() would paste a string key for every range of a
    # genome-wide collection to find the same thing.
    flat <- unlist(object, use.names = FALSE)
    n <- lengths(object)
    key <- as.integer(seqnames(flat)) *
        3L +
        as.integer(BiocGenerics::strand(flat))
    firstKey <- rep(key[(cumsum(n) - n + 1L)[n > 0L]], n[n > 0L])
    bad <- unique(rep(seq_along(object), n)[key != firstKey])
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
        return(GenomeInfoDb::`seqinfo<-`(grl, value = si))
    }
    keepLevels <- intersect(
        GenomeInfoDb::seqlevels(si),
        GenomeInfoDb::seqlevels(grl)
    )
    GenomeInfoDb::`seqinfo<-`(
        grl,
        new2old = match(keepLevels, GenomeInfoDb::seqlevels(grl)),
        value = si[keepLevels]
    )
}

# -----------------------------------------------------------------------------
# Element names
# -----------------------------------------------------------------------------

# The identity columns an element's name is read from, in reading order. A
# name is the present columns' values joined by "|", missing parts skipped, so
# it reads as the row's key: "Q1|brain|ENSG1|susie".
# @noRd
.RTL_NAME_COLUMNS <- c(
    "study",
    "context",
    "trait",
    "method",
    "gwasStudy",
    "gwasContext",
    "gwasTrait",
    "gwasMethod",
    "analysis",
    "cosId",
    "blockId",
    "qtlCs",
    "gwasCs",
    "jointStudies",
    "jointContexts",
    "jointTraits"
)

# Element names from the identity columns, as Bioconductor keys a List's
# elements. Rows sharing an identity -- one tuple split across chromosomes --
# are told apart by their span, the key the uniqueness checks already use.
# NULL when the collection carries no identity columns.
# @noRd
.rtlElementNames <- function(grl) {
    md <- mcols(grl, use.names = FALSE)
    cols <- intersect(.RTL_NAME_COLUMNS, colnames(md))
    if (length(cols) == 0L || length(grl) == 0L) {
        return(NULL)
    }
    parts <- exec(cbind, !!!map(cols, .rtlNameColumn, md = md))
    keys <- apply(parts, 1L, .rtlJoinNameParts)
    shared <- is_in(keys, keys[duplicated(keys)])
    if (!any(shared)) {
        return(keys)
    }
    spans <- .rtlRangeKeys(grl)
    spans <- if_else(is.na(spans), str_c("row", seq_along(keys)), spans)
    if_else(shared, str_c(keys, spans, sep = "|"), keys)
}

# @noRd
.rtlNameColumn <- function(nm, md) {
    as.character(md[[nm]])
}

# @noRd
.rtlJoinNameParts <- function(parts) {
    str_flatten(parts[!is.na(parts) & nzchar(parts)], "|")
}

# `grl` with its element names set from its identity columns, left as it is
# when it carries none.
# @noRd
.rtlNamed <- function(grl) {
    nms <- .rtlElementNames(grl)
    if (is.null(nms)) {
        return(grl)
    }
    `names<-`(grl, nms)
}

# Rebuild `x` from a new list of GRanges elements, carrying every
# collection-level slot over and narrowing mcols by `keep`.
#
# This is the ONLY supported way to produce a modified collection. `endoapply`
# via a `setAs("list", <class>)` coercion looks like it works but silently
# drops the subclass's slots -- verified -- which is exactly the desync this
# class exists to prevent.
# @noRd
.rtlRebuild <- function(x, elements, keep) {
    # mcols and slots are attached BEFORE new(), not after: new() validates
    # during initialize(), and a subclass's validity method reads its identity
    # columns and slots. Building the object bare and filling it in afterwards
    # trips that check on the way past.
    bare <- .rtlBindElements(elements)
    # seqinfo is collection-level state, exactly like the slots below: it
    # carries the genome build. A rebuild from bare elements starts with the
    # seqlevels those elements happen to span -- none at all when everything
    # was dropped -- so the original seqinfo is merged back in, or subsetting
    # to nothing would silently discard the build.
    withSeqinfo <- .rtlRestoreSeqinfo(bare, x)
    md <- mcols(x, use.names = FALSE)
    grl <- if (is.null(md)) {
        withSeqinfo
    } else {
        S4Vectors::`mcols<-`(withSeqinfo, value = md[keep, , drop = FALSE])
    }
    # The rebuilt elements keep the names they had -- the caller may have set
    # their own -- and an unnamed collection is named off its identity.
    grl <- if (is.null(names(x))) {
        .rtlNamed(grl)
    } else {
        `names<-`(grl, names(x)[keep])
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

# One method on the shared base rather than one per flavour: every collection
# that keys rows by LD block stores the id in the same mcols column, and the
# ones that have no block axis want the same NA-padded answer. Serialized
# collections predating the column take that path too.
#' @rdname blockId
#' @export
setMethod("blockId", "RangedTupleList", function(x, ...) {
    .fmrIdentityColumn(x, "blockId")
})

#' @rdname RangedTupleList-methods
#' @export
setMethod("ncol", "RangedTupleList", function(x) {
    # No is.null() arm: `elementMetadata` is a typed DataFrame slot, so a
    # collection with no metadata carries a zero-COLUMN DataFrame rather than
    # NULL, and ncol() already answers 0 for it.
    ncol(mcols(x, use.names = FALSE))
})

#' @rdname RangedTupleList-methods
#' @importFrom BiocGenerics colnames
#' @export
setMethod("colnames", "RangedTupleList", function(x, do.NULL = TRUE, prefix) {
    # Completes the nrow / ncol / `$` set. Without it `colnames(x)` returns
    # NULL rather than erroring, so code that asks a collection for its column
    # names gets a silent wrong answer instead of a signal.
    colnames(mcols(x, use.names = FALSE))
})

#' @rdname RangedTupleList-methods
#' @export
setMethod("$", "RangedTupleList", function(x, name) {
    # A zero-column DataFrame answers NULL for any name, which is the same
    # thing the removed is.null() arm returned.
    mcols(x, use.names = FALSE)[[name]]
})

#' @rdname RangedTupleList-methods
#' @export
setMethod("$<-", "RangedTupleList", function(x, name, value) {
    derived <- is.null(names(x)) || identical(names(x), .rtlElementNames(x))
    md <- `[[<-`(mcols(x, use.names = FALSE), name, value = value)
    out <- S4Vectors::`mcols<-`(x, value = md)
    # Names read off the identity follow an edit to it; names the caller set
    # are the caller's to keep.
    if (derived) .rtlNamed(out) else out
})

# mcols<- and values<- both store their table through this generic, so one
# method checks all three spellings. The inherited method stores the table
# without running validity, which would keep an identity column of the wrong
# type.
#' @rdname RangedTupleList-methods
#' @importFrom S4Vectors elementMetadata<-
#' @export
setMethod("elementMetadata<-", "RangedTupleList", function(x, ..., value) {
    out <- methods::callNextMethod()
    problems <- .rtlIdentityProblems(mcols(out, use.names = FALSE))
    if (length(problems) > 0L) {
        abort(problems)
    }
    out
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
    elements <- `[[<-`(as.list(x), idx, value = value)
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
#' @importFrom checkmate assertScalar
.rtlAssignIndex <- function(x, i) {
    if (is.character(i)) {
        pos <- match(i, names(x))
        if (is.na(pos)) {
            abort(glue("`[[<-`: no element named '{i}'."))
        }
        return(pos)
    }
    assertScalar(i, na.ok = FALSE, .var.name = "`[[<-` index")
    as.integer(i)
}

#' @rdname RangedTupleList-methods
#' @importFrom S4Vectors endoapply
#' @export
setMethod("endoapply", "RangedTupleList", function(X, FUN, ...) {
    # Defined explicitly rather than left to the inherited machinery. Stock
    # `endoapply` on a subclass errors unless a `setAs("list", <class>)`
    # coercion exists, and adding that coercion makes it run while silently
    # DROPPING the subclass's slots (verified: genome / ldReference come back
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

# The inherited GRangesList sort() / unique() / trim() work per element too,
# but rebuild the result as a bare CompressedGRangesList, dropping the class
# and every collection-level slot. Routed through endoapply() they keep both.
#' @rdname RangedTupleList-methods
#' @importFrom BiocGenerics sort
#' @export
setMethod("sort", "RangedTupleList", function(x, decreasing = FALSE, ...) {
    endoapply(x, sort, decreasing = decreasing, ...)
})

#' @rdname RangedTupleList-methods
#' @importFrom BiocGenerics unique
#' @export
setMethod(
    "unique",
    "RangedTupleList",
    function(x, incomparables = FALSE, ...) {
        endoapply(x, unique, incomparables = incomparables, ...)
    }
)

#' @rdname RangedTupleList-methods
#' @importFrom IRanges trim
#' @export
setMethod("trim", "RangedTupleList", function(x, use.names = TRUE, ...) {
    endoapply(x, trim, use.names = use.names, ...)
})

# A stable per-element key for the element's own genomic span, used where an
# identity tuple needs a positional component. Derived from range() rather than
# from a stored label so it cannot drift out of step with the variants, and so
# it stays correct after an element is narrowed. The keys are values, so they
# do not carry the collection's element names. An empty element has no span
# and keys NA.
# @noRd
.rtlRangeKeys <- function(x) {
    if (length(x) == 0L) {
        return(character(0))
    }
    spans <- .rtlFirstSpans(x)
    keys <- rep(NA_character_, length(x))
    keys[spans$has] <- .rtlSpanKey(spans$first)
    keys
}

# Each element's span -- its first range under range(), which every element
# of a collection has exactly one of (.rtlCheckSingleSpan()) -- in element
# order, and which elements have one. One range() over the whole list: one
# per element made building the keys cost a millisecond a row.
# @noRd
.rtlFirstSpans <- function(x) {
    # Each element sits on one (seqname, strand) (.rtlCheckSingleSpan()), so
    # its span is its smallest start to its largest end there: taken from
    # integer vectors, where range() would paste a key for every range.
    n <- lengths(x)
    has <- n > 0L
    flat <- unlist(x, use.names = FALSE)
    group <- rep(seq_along(x), n)
    first <- (cumsum(n) - n + 1L)[has]
    starts <- vapply(split(start(flat), group), min, integer(1))
    ends <- vapply(split(end(flat), group), max, integer(1))
    spans <- GenomicRanges::GRanges(
        seqnames(flat)[first],
        IRanges::IRanges(starts, ends),
        strand = BiocGenerics::strand(flat)[first],
        seqinfo = GenomeInfoDb::seqinfo(flat)
    )
    list(first = spans, has = has)
}

# The key of the first range of `g`, a span: callers pass range() of a
# GRanges that has variants (an empty GWAS entry is skipped before it is
# keyed).
# @noRd
.rtlOneRangeKey <- function(g) {
    .rtlSpanKey(g[1L])
}

# "chrom_start_end" for each range of `gr`: the one formatter every span key
# goes through, so a key minted for one range and the keys of a collection
# agree by construction.
# @noRd
.rtlSpanKey <- function(gr) {
    str_c(as.character(seqnames(gr)), start(gr), end(gr), sep = "_")
}

# The seqname each element of a GRangesList sits on, NA for an empty one, as
# .rtlElementSeqname() reads one element -- read off the joined elements
# rather than element by element.
# @noRd
.rtlElementSeqnames <- function(x) {
    n <- lengths(x)
    has <- n > 0L
    flat <- unlist(x, use.names = FALSE)
    out <- rep(NA_character_, length(x))
    out[has] <- as.character(seqnames(flat))[(cumsum(n) - n + 1L)[has]]
    out
}

# Whether the elements' variant metadata has column `cn`. A collection's
# elements share one variant metadata table, so all of them have it or none
# does.
# @noRd
.rtlHasVariantColumn <- function(x, cn) {
    md <- mcols(unlist(x, use.names = FALSE), use.names = FALSE)
    is_in(cn, colnames(md))
}

# -----------------------------------------------------------------------------
# updateObject
# -----------------------------------------------------------------------------

#' @title Update an Object Saved by an Earlier pecotmr
#' @description Bring an object saved by an earlier version of pecotmr up to
#'   date. An up-to-date object is returned unchanged; an older one gets every
#'   update below that applies to it, in the one call.
#'
#'   \itemize{
#'     \item Collections keep their LD references per row: each distinct
#'       reference once in the \code{ldReference} slot, and an
#'       \code{ldReference} column naming the one every row uses (in
#'       \code{outcomeInfo} for a \code{ColocBoostResult}). An object saved
#'       when it held a single reference for all its rows is converted by
#'       keying every row to that reference (for a \code{ColocBoostResult},
#'       every summary-statistic outcome), and one that held none by marking
#'       every row LD-free. An object saved when the slot was named
#'       \code{ldSketch} has it renamed first.
#'     \item A genotype panel saved before panels recorded where their
#'       genotypes come from (an LD reference, or a \code{QtlDataset}'s
#'       genotypes) is read again from its files, as
#'       \code{\link{readGenotypes}} reads them, and cut to the variants and
#'       samples it held: a VCF, GDS or PLINK1 panel then reads through its
#'       format's own backend, and a stochastic sketch records its real
#'       sample count and allele frequencies. Only the files' variant and
#'       sample lists, and a sketch's \code{.afreq}, are read; no genotypes
#'       are, and nothing is written. An LD reference whose files cannot be
#'       read again is kept as saved, with a warning, and updating again once
#'       they are back finishes the job; a \code{QtlDataset} whose genotypes
#'       cannot be is an error.
#'     \item Collections name their elements by their identity tuple; one
#'       saved without names gets them.
#'     \item Summary statistics keep their QC settings in \code{qcOptions}
#'       and each row's QC audit in a \code{qcAudit} column; ones saved with
#'       both in one \code{qcInfo} record have it split between the two.
#'       Their genome build is in \code{seqinfo()}; ones saved with it in a
#'       \code{genome} slot have it moved there.
#'     \item TWAS weights record the variants each row was fitted on in a
#'       \code{fitVariants} column; ones saved without it have it filled in
#'       from their fits. Fine-mapping results no longer record them (a row's
#'       variants are its fit's variables); ones saved with the column have
#'       it dropped.
#'     \item A \code{CtwasResult} holds one element per region; one saved as a
#'       table of whole runs is rebuilt.
#'     \item A \code{QtlDataset} no longer stores run settings (genotype
#'       QC cutoffs, \code{scaleResiduals}); one saved with them has them
#'       dropped, as the pipelines' \code{\link{GenotypeFilterParam}} and
#'       \code{\link{ResidualizationParam}} now carry them. One saved with
#'       its study in a \code{study} slot has it moved to \code{studyName}.
#'     \item A \code{MultiStudyQtlDataset} is now a list of its studies; one
#'       saved when it held them in a slot is rebuilt as one, its studies
#'       updated too.
#'     \item An \code{LdScore} or \code{LdEigen} saved before it could keep
#'       an LD reference stores its blocks' LD itself, and gets none; one
#'       saved with the reference's sample size as \code{nRef} has it moved
#'       to \code{nSamples}.
#'     \item An \code{H2Estimate} saved with its estimates under the slots'
#'       earlier names (\code{h2}, \code{h2Se}, \code{local},
#'       \code{tauBlocks}) has them moved to the current ones.
#'   }
#' @param object A collection (summary statistics, fine-mapping results, TWAS
#'   weights, a colocalization or cTWAS result), a \code{QtlDataset}, a
#'   \code{MultiStudyQtlDataset}, an \code{LdScore} or \code{LdEigen}, or an
#'   \code{H2Estimate}.
#' @param ... Passed on to the inherited method.
#' @param verbose Report the update.
#' @return \code{object} in the current representation.
#' @examples
#' data(qtlSumStatsExample)
#' # An object built by this version needs nothing.
#' identical(updateObject(qtlSumStatsExample), qtlSumStatsExample)
#' @name updateObject-methods
#' @rdname updateObject-methods
NULL

#' @rdname updateObject-methods
#' @importFrom BiocGenerics updateObject
#' @export
setMethod(
    "updateObject",
    "RangedTupleList",
    function(object, ..., verbose = FALSE) {
        if (isTRUE(verbose)) {
            message("updateObject(object = '", class(object)[[1L]], "')")
        }
        methods::callNextMethod(
            .rtlNamed(.updateRowLdRefs(object)),
            ...,
            verbose = verbose
        )
    }
)

# `object` with each slot an earlier version saved under another name moved
# to its current one: `renames` maps each old name to the new. The inherited
# update keeps only the slots a class still has, so a value left under its
# old name would be dropped. A name the object was not saved under is
# skipped.
# @noRd
#' @importFrom purrr reduce2
.renameSavedSlots <- function(object, renames) {
    reduce2(names(renames), unname(renames), .renameSavedSlot, .init = object)
}

# slot() rather than attr(): R stores a NULL slot as a placeholder symbol,
# which slot() reads back as NULL and attr() hands back unconverted.
# @noRd
.renameSavedSlot <- function(object, from, to) {
    if (is.null(attr(object, from, exact = TRUE))) {
        return(object)
    }
    value <- methods::slot(object, from)
    attr(object, from) <- NULL
    `slot<-`(object, to, check = FALSE, value = value)
}

# -----------------------------------------------------------------------------
# Construction helper: enforce one seqname per element
# -----------------------------------------------------------------------------

# Each entry's variants become one ELEMENT of the GRangesList, with that
# entry's metadata row as the element's outer mcols. A multi-seqname entry
# splits into one element per chromosome, its metadata row replicated
# alongside -- which is why the metadata is indexed by `fromIdx` rather
# than used as-is. Shared by the TwasWeights and *FineMappingResult
# constructors, whose assembly differs only in the columns handed in.
# @noRd
.rtlGrlWithMetadata <- function(payloads, cols) {
    split <- .rtlSplitBySeqname(map(payloads, variants))
    md <- exec(S4Vectors::DataFrame, !!!c(cols, list(check.names = FALSE)))
    S4Vectors::`mcols<-`(
        .rtlBindElements(split$entry),
        value = md[split$fromIdx, , drop = FALSE]
    )
}

# Build a RangedTupleList subclass from its payloads and metadata columns --
# the shared tail of the TwasWeights and QtlFineMappingResult constructors.
# new() validates it, once. `ldReference` is one reference for every payload
# or a list with one per payload; each distinct one is stored once, and every
# row names its own in the `ldReference` column.
# @noRd
.rtlNewValidated <- function(Class, payloads, cols, ldReference) {
    dict <- .ldRefDict(.ldRefPerTuple(ldReference, length(payloads)))
    new(
        Class,
        .rtlNamed(.rtlGrlWithMetadata(
            payloads,
            c(cols, list(ldReference = dict$keys))
        )),
        ldReference = dict$refs
    )
}

# -----------------------------------------------------------------------------
# Joining elements into a list
# -----------------------------------------------------------------------------

# GenomicRanges::GRangesList(elements), built faster. GRangesList() joins
# the elements' variant metadata through S4Vectors, which reads every column
# of every element by S4 dispatch, about 0.25 ms a cell: 3 s for a 1,000-row
# collection of fine-mapping fits. When every element is a plain GRanges
# whose metadata holds the same plain vector columns, in the same order and
# of the same types -- what the collection constructors pass -- they are
# joined from their vectors instead, giving the identical object. Anything
# else takes GRangesList(). Like GRangesList(), which builds its result inside
# suppressWarnings(), the join reports no warnings: merging the elements'
# seqinfo warns about every two on different chromosomes.
# @noRd
.rtlBindElements <- function(elements) {
    flat <- suppressWarnings(.rtlJoinPlainElements(elements))
    if (is.null(flat)) {
        return(GenomicRanges::GRangesList(elements))
    }
    IRanges::relist(
        flat,
        IRanges::PartitioningByEnd(
            cumsum(lengths(elements)),
            names = names(elements)
        )
    )
}

# The elements of collections `parts`, in order, as one GRangesList: what
# GRangesList() makes of their elements, keeping their names. Each part is
# joined whole -- its elements are already one GRanges underneath -- rather
# than taken apart with as.list() and put back together, which cost about
# 0.6 ms an element. c() joins the parts' variant metadata as GRangesList()
# would, filling a column one part lacks with NA, and its warnings are
# silenced as GRangesList() silences them (.rtlBindElements()).
# @noRd
.rtlJoinParts <- function(parts) {
    flats <- unname(map(parts, unlist, use.names = FALSE))
    flat <- suppressWarnings(exec(c, !!!flats))
    ends <- cumsum(list_c(map(parts, .rtlElementLengths), ptype = integer()))
    IRanges::relist(
        flat,
        IRanges::PartitioningByEnd(ends, names = .rtlJoinedNames(parts))
    )
}

# @noRd
.rtlElementLengths <- function(x) {
    unname(lengths(x))
}

# The elements joined into one GRanges as c() joins them, or NULL when one of
# them falls outside the plain case .rtlBindElements() builds.
# @noRd
.rtlJoinPlainElements <- function(elements) {
    elements <- unname(as.list(elements))
    if (
        length(elements) == 0L ||
            !all(map_lgl(elements, .rtlIsPlainGRanges))
    ) {
        return(NULL)
    }
    columns <- map(elements, .rtlMcolList)
    if (!.rtlSameVectorColumns(columns)) {
        return(NULL)
    }
    flat <- GenomicRanges::GRanges(
        seqnames = list_c(
            map(elements, .rtlSeqnameValues),
            ptype = character()
        ),
        ranges = IRanges::IRanges(
            start = list_c(map(elements, start), ptype = integer()),
            end = list_c(map(elements, end), ptype = integer()),
            names = .rtlJoinedNames(elements)
        ),
        strand = list_c(map(elements, .rtlStrandValues), ptype = character()),
        seqinfo = .rtlMergedSeqinfo(elements)
    )
    .rtlWithJoinedMcols(flat, columns)
}

# A GRanges itself, not a subclass, with no metadata() of its own, which c()
# would carry over from the first element.
# @noRd
.rtlIsPlainGRanges <- function(g) {
    methods::is(g, "GRanges") &&
        identical(class(g)[[1L]], "GRanges") &&
        length(S4Vectors::metadata(g)) == 0L
}

# One element's variant metadata as a plain list of its columns.
# @noRd
.rtlMcolList <- function(g) {
    as.list(mcols(g, use.names = FALSE))
}

# TRUE when every element's metadata holds the first one's plain vector
# columns: the same names, in the same order, of the same types.
# @noRd
.rtlSameVectorColumns <- function(columns) {
    first <- columns[[1L]]
    all(map_lgl(
        columns,
        .rtlColumnsMatch,
        names = names(first) %||% character(0),
        types = map_chr(first, typeof)
    ))
}

# @noRd
.rtlColumnsMatch <- function(cols, names, types) {
    identical(names(cols) %||% character(0), names) &&
        identical(map_chr(cols, typeof), types) &&
        all(map_lgl(cols, .rtlIsPlainVector))
}

# An atomic vector with no attributes beyond its names -- no dimensions or
# class -- so joining it is plain concatenation, names included.
# @noRd
.rtlIsPlainVector <- function(v) {
    is.atomic(v) && all(is_in(names(attributes(v)), "names"))
}

# @noRd
.rtlSeqnameValues <- function(g) {
    as.character(seqnames(g))
}

# @noRd
.rtlStrandValues <- function(g) {
    as.character(GenomicRanges::strand(g))
}

# The elements' variant names joined, "" for an element without any; NULL
# when none has any, as c() leaves them.
# @noRd
.rtlJoinedNames <- function(elements) {
    nms <- map(elements, names)
    if (all(map_lgl(nms, is.null))) {
        return(NULL)
    }
    list_c(map2(nms, elements, .rtlNamesOrBlank), ptype = character())
}

# @noRd
.rtlNamesOrBlank <- function(nm, g) {
    nm %||% rep("", length(g))
}

# The elements' seqinfo merged as c() merges it, each distinct one once, in
# the order they appear.
# @noRd
.rtlMergedSeqinfo <- function(elements) {
    infos <- unique(map(elements, GenomeInfoDb::seqinfo))
    reduce(infos, S4Vectors::merge)
}

# `flat` carrying the elements' metadata columns, each joined once.
# @noRd
.rtlWithJoinedMcols <- function(flat, columns) {
    cn <- names(columns[[1L]]) %||% character(0)
    if (length(cn) == 0L) {
        return(flat)
    }
    joined <- set_names(map(cn, .rtlJoinColumn, columns = columns), cn)
    S4Vectors::`mcols<-`(
        flat,
        value = S4Vectors::DataFrame(joined, check.names = FALSE)
    )
}

# @noRd
.rtlJoinColumn <- function(cn, columns) {
    list_c(map(columns, cn))
}

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
        # unname(): the flattened names were an artifact of base unlist()'s
        # `outer.inner` mangling and were never coherent -- a multi-seqname
        # entry got `a.chr1`, a single-seqname one got bare `a`, an unnamed
        # input got `chr1` or "". Nothing reads them, and the seqname path
        # already had to unname() derived values to keep them out of mcols.
        entry = unname(list_flatten(pieces)),
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

# Split each element across a set of windows (a GWAS `region`, typically LD
# blocks) instead of by seqname, returning the same (entry, fromIdx) shape
# plus the key and the interval of the window each piece belongs to.
#
# A genome-wide GWAS split by seqname gives one element per chromosome, which
# is too coarse for cTWAS: its EM needs per-block context (§4.7). `blocks` is a
# GRanges of windows whose names supply the keys.
#
# A variant matching no window is dropped rather than pooled into a
# catch-all: the windows are what the caller asked to analyse, and a variant
# outside every LD block has no block-local LD to be fine-mapped against. The
# count is reported so the loss is never silent.
# @noRd
.rtlSplitByBlocks <- function(entry, blocks) {
    entry <- as.list(entry)
    if (length(entry) == 0L) {
        return(list(
            entry = list(),
            fromIdx = integer(0),
            blockId = character(0),
            blockRange = GenomicRanges::GRanges()
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
            "{dropped} variant(s) lie outside every `region` window and ",
            "were dropped."
        )
        warn(msg)
    }
    keys <- unname(list_c(map(pieces, names), ptype = character()))
    list(
        # unname(): the flattened names were an artifact of base unlist()'s
        # `outer.inner` mangling and were never coherent -- a multi-seqname
        # entry got `a.chr1`, a single-seqname one got bare `a`, an unnamed
        # input got `chr1` or "". Nothing reads them, and the seqname path
        # already had to unname() derived values to keep them out of mcols.
        entry = unname(list_flatten(pieces)),
        fromIdx = rep(seq_along(entry), lengths(pieces)),
        blockId = keys,
        blockRange = unname(blocks[keys])
    )
}

# Every window needs a name: the names become `blockId`, which is what
# downstream code keys regions by.
# @noRd
.rtlValidateBlocks <- function(blocks) {
    if (!methods::is(blocks, "GRanges")) {
        msg <- glue(
            "`region` must be a GRanges of windows ",
            "(got {class(blocks)[[1L]]})."
        )
        abort(msg)
    }
    if (length(blocks) == 0L) {
        abort("`region` is empty; it must describe at least one window.")
    }
    keys <- .rtlBlockKeys(blocks)
    if (anyDuplicated(keys) > 0L) {
        dup <- unique(keys[duplicated(keys)])
        msg <- glue(
            "`region` window names must be unique; repeated: ",
            "{str_flatten(head(dup, 5L), ', ')}."
        )
        abort(msg)
    }
    # `names<-` not set_names(): `blocks` is a GRanges, and set_names() is
    # vector-only -- it rejects an S4 object.
    `names<-`(blocks, keys)
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
    hits <- IRanges::findOverlaps(
        .locusRanges(g),
        .locusRanges(blocks),
        select = "first"
    )
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
        # bindROWS() is the one hook `c()` and `append()` share, so defining
        # it here fixes both. The inherited CompressedList version rbinds the
        # mcols, requiring identical columns, so it fails whenever the parts
        # came from runs with different optional columns (traitPos /
        # jointContexts / blockId); it also carries the FIRST part's
        # collection-level slots silently, making the result depend on
        # argument order.
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
#' (\code{\link{susieFit}}, \code{\link{cvResult}}) for those.
#'
#' \code{as(x, "GRanges")} is the same conversion.
#'
#' @param x A \code{RangedTupleList}.
#' @return A \code{GRanges}.
#' @examples
#' data(qtlFineMappingExample)
#' flat <- flattenTupleRanges(qtlFineMappingExample)
#' head(names(S4Vectors::mcols(flat)))
#' identical(as(qtlFineMappingExample, "GRanges"), flat)
#' @export
flattenTupleRanges <- function(x) {
    if (!methods::is(x, "RangedTupleList")) {
        msg <- glue(
            "`x` must be a RangedTupleList (got {class(x)[[1L]]})."
        )
        abort(msg)
    }
    gr <- if (methods::is(x, "FineMappingResultBase")) {
        .fmrGatherVariants(x, seq_len(length(x)))
    } else {
        .rtlGatherElements(x, seq_len(length(x)))
    }
    md <- mcols(x, use.names = FALSE)
    if (is.null(md) || ncol(md) == 0L || length(gr) == 0L) {
        return(gr)
    }
    tupleCols <- names(md)[map_lgl(as.list(md), is.atomic)]
    if (length(tupleCols) == 0L) {
        return(gr)
    }
    reps <- rep(seq_len(length(x)), lengths(x))
    broadcast <- set_names(
        map(md[tupleCols], .rtlBroadcastColumn, reps = reps),
        map_chr(tupleCols, .rtlDotName)
    )
    `mcols<-`(
        gr,
        value = cbind(
            mcols(gr, use.names = FALSE),
            DataFrame(broadcast, check.names = FALSE)
        )
    )
}

#' @rdname flattenTupleRanges
#' @name coerce,RangedTupleList,GRanges-method
#' @aliases coerce,RangedTupleList,GRanges-method
#' @usage NULL
#' @importFrom methods coerce
#' @exportMethod coerce
setAs("RangedTupleList", "GRanges", function(from) {
    flattenTupleRanges(from)
})

# @noRd
.rtlBroadcastColumn <- function(col, reps) {
    col[reps]
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
#' identity tuple they carry -- and, for a tuple split across chromosomes, by
#' chromosome -- and returned to \code{template}'s elements, so the
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
    keys <- .rtlNestKeys(flat, template, .rtlTupleKeyCols(template))
    elements <- map(keys$wanted, .rtlPickByKey, flat = flat, have = keys$have)
    if (methods::is(template, "FineMappingResultBase")) {
        elements <- .fmrDropRestoredColumns(elements, template)
    }
    .rtlRebuild(template, elements, seq_len(length(template)))
}

# Re-nested fine-mapping elements without the columns flattening derived back
# (see flattenTupleRanges()): the pip, unless the template stored it, and the
# credible-set columns each row records as derived -- each only while the
# row's fit still gives exactly those values, so a column a verb changed is
# kept.
# @noRd
.fmrDropRestoredColumns <- function(elements, template) {
    storedPip <- .rtlHasVariantColumn(template, "pip")
    map2(
        elements,
        seq_along(elements),
        .fmrDropRestoredAt,
        template = template,
        storedPip = storedPip
    )
}

# @noRd
.fmrDropRestoredAt <- function(g, i, template, storedPip) {
    parts <- .fmrRowParts(template, i)
    recorded <- .fmRowDerivedColumns(parts)
    if (length(g) == 0L) {
        md <- mcols(g, use.names = FALSE)
        labelCols <- .fmRowCsRecord(recorded)$labels
        derived <- c(
            if (storedPip) character(0) else "pip",
            labelCols,
            str_c(labelCols, "_purity")
        )
        return(`mcols<-`(
            g,
            value = md[, setdiff(colnames(md), derived), drop = FALSE]
        ))
    }
    fit <- susieFit(parts)
    cols <- .fitColumnsFor(g, fit)
    if (!storedPip) {
        g <- .fmRowDropFitPip(g, fit, cols)
    }
    if (length(recorded) == 0L) {
        return(g)
    }
    cs <- .fmRowDropFitCs(g, fit, cols)
    if (identical(cs$derived, recorded)) cs$variants else g
}

# Element and range keys for re-nesting. The identity tuple alone keys an
# element unless several elements share it -- one tuple split across
# chromosomes -- and then the one seqname each of those spans is added on both
# sides. An empty element spans none, so its key is NA and it stays empty.
# @noRd
.rtlNestKeys <- function(flat, template, keyCols) {
    wanted <- .rtlTupleKeys(
        mcols(template, use.names = FALSE),
        keyCols,
        n = length(template)
    )
    have <- .rtlTupleKeys(
        mcols(flat, use.names = FALSE),
        map_chr(keyCols, .rtlDotName),
        n = length(flat)
    )
    shared <- unique(wanted[duplicated(wanted)])
    if (length(shared) == 0L) {
        return(list(wanted = wanted, have = have))
    }
    list(
        wanted = .rtlKeyWithSeqname(
            wanted,
            is_in(wanted, shared),
            map_chr(as.list(template), .rtlElementSeqname)
        ),
        have = .rtlKeyWithSeqname(
            have,
            is_in(have, shared),
            as.character(seqnames(flat))
        )
    )
}

# @noRd
.rtlKeyWithSeqname <- function(keys, add, seqs) {
    if_else(add, str_c(keys, seqs, sep = "\r"), keys)
}

# The identity columns to group by: the atomic mcols the flattener broadcasts.
# @noRd
.rtlTupleKeyCols <- function(x) {
    # No is.null() arm: see ncol() above -- mcols() answers a zero-column
    # DataFrame, and names() of one is already character(0).
    md <- mcols(x, use.names = FALSE)
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
    flat[which(have == key)]
}

# The one seqname an element sits on. NA for an empty element, which carries
# no coordinate to name.
# @noRd
.rtlElementSeqname <- function(g) {
    if (length(g) == 0L) {
        return(NA_character_)
    }
    as.character(seqnames(g))[[1L]]
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
    # map() forwards `...` to the function, so the user's slice expressions
    # reach dplyr::slice directly -- no wrapper needed.
    elements <- map(as.list(.data), dplyr::slice, ...)
    .rtlRebuild(.data, elements, seq_len(length(.data)))
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
    if (length(missing) == 0L) {
        return(out)
    }
    `mcols<-`(
        out,
        value = cbind(
            mcols(out, use.names = FALSE),
            mcols(flat, use.names = FALSE)[, missing, drop = FALSE]
        )
    )
}
