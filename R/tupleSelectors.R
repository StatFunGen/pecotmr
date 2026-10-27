# =============================================================================
# Tuple row matchers
# -----------------------------------------------------------------------------
# Internal helpers shared by the fine-mapping result / TwasWeights / SumStats
# collections to resolve a tuple-keyed selection to a single row index. Every
# collection is a RangedTupleList, so these read identity columns off mcols
# and payloads off the elements -- there is no second shape to branch on.
# Pure R helpers -- no S4 dispatch, no exports.
# =============================================================================

# Internal: return integer row indices of `x` where every (column, value)
# pair in `keys` matches as.character(x[[column]]) == value. Shared
# building block for tuple-keyed row selectors and cache lookups
# (.tupleSelectRow, .qtlSumStatsSelectRow, .gwasSelectStudy,
# .fmCacheLookup, .cipFmrHasTuple, etc.). Pure vectorised AND-match with
# character coercion -- no validation, no error reporting.
.matchTupleRows <- function(x, keys) {
    if (length(keys) == 0L) {
        return(seq_len(nrow(x)))
    }
    which(reduce(
        map(names(keys), .tupleKeyEquals, x = x, keys = keys),
        `&`,
        .init = rep(TRUE, nrow(x))
    ))
}

# @noRd
.tupleKeyEquals <- function(k, x, keys) {
    as.character(.tupleColumn(x, k)) == keys[[k]]
}

# Read an identity column by name. Every collection keeps its identity
# columns in mcols; `x[[k]]` addresses an ELEMENT, not a column, so mcols is
# the only correct read.
# @noRd
.tupleColumn <- function(x, k) {
    mcols(x)[[k]]
}

# ---- flavour-agnostic identity reads ----------------------------------------
# `context` and `trait` are QTL-only axes: a GwasFineMappingResult is keyed by
# (study, method) alone. Pipelines that pair the two flavours -- colocPipeline,
# qtlEnrichmentPipeline -- need "this side has no context" as a VALUE, so these
# report an absent axis as NA instead of erroring on a missing column.

# One identity column of a collection, as character, NA-filled when the
# collection has no such axis (length always matches the collection).
# @noRd
.fmrIdentityColumn <- function(x, column) {
    values <- .tupleColumn(x, column)
    if (is.null(values)) {
        return(rep(NA_character_, nrow(x)))
    }
    as.character(values)
}

# The same, for a single row. Indexed before coercion so reading one row does
# not convert the whole column.
# @noRd
.fmrIdentityAt <- function(x, column, i) {
    values <- .tupleColumn(x, column)
    if (is.null(values)) {
        return(NA_character_)
    }
    as.character(values[[i]])
}

# Which flavour of fine-mapping collection a side is, for messages.
# @noRd
.fmrSideName <- function(x) {
    if (methods::is(x, "QtlFineMappingResult")) "QTL" else "GWAS"
}

# A human-readable identity for warnings: the side's flavour plus the identity
# fields it actually carries (the axes a GWAS side does not have are dropped
# rather than reported as NA).
# @noRd
.fmrTupleLabel <- function(side, ident, block = NULL) {
    supplied <- compact(c(ident, list(block = block)))
    fields <- discard(supplied, is.na)
    shown <- str_c(names(fields), "='", map_chr(fields, as.character), "'")
    glue("{side} ({str_flatten(shown, ', ')})")
}

# One row's payload, whichever collection `x` is -- what its constructor
# takes as `entry`: for TWAS weights and fine-mapping results, the row parts
# rebuilt from its element and the collection's mcols, where the fits sit;
# for the summary-statistics collections, the element itself. A caller
# taking either weight source reads its rows through here.
# @noRd
.collectionEntry <- function(x, i) {
    if (methods::is(x, "TwasWeights")) {
        return(.twrRowParts(x, i))
    }
    if (methods::is(x, "FineMappingResultBase")) {
        return(.fmrRowParts(x, i))
    }
    x[[i]]
}

# Every row's payload (.collectionEntry()) as a plain list.
# @noRd
.collectionEntries <- function(x) {
    # Weight and fine-mapping rows are rebuilt one by one: their elements
    # alone are bare GRanges, without the fits.
    if (
        methods::is(x, "TwasWeights") ||
            methods::is(x, "FineMappingResultBase")
    ) {
        return(map(seq_len(nrow(x)), .collectionEntry, x = x))
    }
    as.list(x)
}

# The identity column NAMES. `names(x)` would be the ELEMENT names, not the
# columns, so this reads mcols.
# @noRd
.tupleColumnNames <- function(x) {
    colnames(mcols(x))
}

# Internal: resolve a tuple-keyed selection (study, context, trait,
# method) to a single row index. Used by the QtlFineMappingResult and
# TwasWeights accessors. Returns an error when no row matches; returns
# the single row index when the collection has exactly one row and any
# selector argument was omitted.
#' @importFrom checkmate assertVector
.tupleSelectRow <- function(
    x,
    studyName,
    context,
    trait,
    method,
    cls = "QtlFineMappingResult"
) {
    if (nrow(x) == 0L) {
        msg <- glue("{cls} has no rows.")
        abort(msg)
    }
    anyUnset <- missing(studyName) ||
        is.null(studyName) ||
        missing(context) ||
        is.null(context) ||
        missing(trait) ||
        is.null(trait) ||
        missing(method) ||
        is.null(method)
    if (anyUnset) {
        if (nrow(x) == 1L) {
            return(1L)
        }
        msg <- glue(
            "{cls} has {nrow(x)} entries. Pass `studyName`, `context`, ",
            "`trait`, and `method` to select one."
        )
        abort(msg)
    }
    assertVector(studyName, len = 1L)
    assertVector(context, len = 1L)
    assertVector(trait, len = 1L)
    assertVector(method, len = 1L)
    .tupleMatchQtl(x, studyName, context, trait, method)
}

# Resolve a (study, context, trait, method) tuple to a single row index.
# @noRd
.tupleMatchQtl <- function(x, studyName, context, trait, method) {
    idx <- .matchTupleRows(
        x,
        list(
            study = studyName,
            context = context,
            trait = trait,
            method = method
        )
    )
    if (length(idx) == 0L) {
        msg <- glue(
            "No entry for (study='{studyName}', context='{context}', ",
            "trait='{trait}', method='{method}')."
        )
        abort(msg)
    }
    idx[[1L]]
}

# Internal: resolve a (study, method, blockId) tuple to a single row
# index of a GwasFineMappingResult collection. `blockId` may be NULL when the
# (study, method) pair maps to a single row; otherwise it names one of a
# genome-wide collection's per-block rows.
# @noRd
#' @importFrom checkmate assertVector
.tupleSelectRowGwasFmr <- function(
    x,
    studyName = NULL,
    method = NULL,
    blockId = NULL
) {
    if (nrow(x) == 0L) {
        abort("GwasFineMappingResult has no rows.")
    }
    if (is.null(studyName) || is.null(method)) {
        if (nrow(x) == 1L) {
            return(1L)
        }
        msg <- glue(
            "GwasFineMappingResult has {nrow(x)} entries. Pass `studyName` ",
            "and `method` to select one."
        )
        abort(msg)
    }
    assertVector(studyName, len = 1L)
    assertVector(method, len = 1L)
    assertVector(blockId, len = 1L, null.ok = TRUE)
    .tupleMatchGwas(x, studyName, method, blockId)
}

# Internal: every row of one (study, method) fit of a GwasFineMappingResult
# -- one per LD block -- or the block `blockId` names. `study` and `method`
# may be left NULL only while they still match a single (study, method)
# pair; a selection spanning several fits is refused rather than returned as
# if it were one.
# @noRd
.tupleSelectRowsGwasFmr <- function(x, studyName, method, blockId = NULL) {
    if (nrow(x) == 0L) {
        abort("GwasFineMappingResult has no rows.")
    }
    assertString(studyName, null.ok = TRUE)
    assertString(method, null.ok = TRUE)
    assertString(blockId, null.ok = TRUE)
    idx <- .matchTupleRows(
        x,
        compact(list(study = studyName, method = method, blockId = blockId))
    )
    if (length(idx) == 0L) {
        abort(glue(
            "No entry for ",
            "({.tupleGwasSelectionLabel(studyName, method, blockId)})."
        ))
    }
    fits <- unique(str_c(
        .fmrIdentityColumn(x, "study")[idx],
        .fmrIdentityColumn(x, "method")[idx],
        sep = "/"
    ))
    if (length(fits) > 1L) {
        abort(glue(
            "GwasFineMappingResult holds {length(fits)} (study, method) fits ",
            "matching the selection ({str_flatten(fits, ', ')}); pass ",
            "`studyName` and `method` to select one."
        ))
    }
    idx
}

# The selectors as they appear in a no-match message; an unset study or
# method matched anything.
# @noRd
.tupleGwasSelectionLabel <- function(studyName, method, blockId) {
    parts <- c(
        study = studyName %||% "<any>",
        method = method %||% "<any>",
        blockId = blockId
    )
    str_flatten(str_c(names(parts), "='", parts, "'"), ", ")
}

# Resolve a (study, method[, blockId]) tuple to a single row index.
# @noRd
.tupleMatchGwas <- function(x, studyName, method, blockId) {
    keys <- c(
        list(study = studyName, method = method),
        compact(list(blockId = blockId))
    )
    idx <- .matchTupleRows(x, keys)
    if (length(idx) == 0L) {
        blockPart <- if (is.null(blockId)) {
            ""
        } else {
            glue(", blockId='{blockId}'")
        }
        msg <- glue(
            "No entry for (study='{studyName}', method='{method}'{blockPart})."
        )
        abort(msg)
    }
    if (length(idx) > 1L) {
        .tupleGwasAmbiguous(x, studyName, method, idx)
    }
    idx[[1L]]
}

# Multiple (study, method) rows matched: report the blocks to choose from.
# @noRd
.tupleGwasAmbiguous <- function(x, studyName, method, idx) {
    blocks <- str_flatten(
        shQuote(as.character(x$blockId[idx])),
        ", "
    )
    msg <- glue(
        "GwasFineMappingResult has {length(idx)} rows matching ",
        "(study='{studyName}', method='{method}'); pass `blockId` to ",
        "disambiguate (available: {blocks})."
    )
    abort(msg)
}

# Internal: row indices of a FineMappingResultBase collection matching the
# given selectors. NULL selectors, and selectors naming a column the
# collection lacks, are ignored -- QTL rows carry study/context/trait/method,
# GWAS rows carry study/method/blockId -- so the same call works on either.
# Returns all rows when nothing constrains. Unlike the single-row selectors
# above this never errors on ambiguity; it is the aggregate counterpart used
# by topLoci().
.fmrRowsMatching <- function(
    x,
    studyName = NULL,
    context = NULL,
    trait = NULL,
    method = NULL,
    blockId = NULL
) {
    supplied <- compact(list(
        study = studyName,
        context = context,
        trait = trait,
        method = method,
        blockId = blockId
    ))
    keys <- supplied[is_in(names(supplied), .tupleColumnNames(x))]
    if (length(keys) == 0L) {
        return(seq_len(nrow(x)))
    }
    which(reduce(
        map(names(keys), .tupleKeyIsIn, x = x, keys = keys),
        `&`,
        .init = rep(TRUE, nrow(x))
    ))
}

# @noRd
.tupleKeyIsIn <- function(k, x, keys) {
    is_in(as.character(.tupleColumn(x, k)), as.character(keys[[k]]))
}

# Internal: per-row identity metadata for a FineMappingResultBase collection,
# as a data.frame with a stable column set (study, context, trait, blockId,
# method). Columns the collection lacks are NA-filled so QTL and GWAS results
# yield the same metadata shape. Safe for zero-row collections.
.fmrRowMetadata <- function(x) {
    cols <- c(
        "study",
        "context",
        "trait",
        "blockId",
        "method",
        "dataType",
        "dataSource",
        "standardized"
    )
    n <- nrow(x)
    as_tibble(set_names(map(cols, .fmrMetadataCol, x = x, n = n), cols))
}

# Internal: row-bind view tibbles whose column sets may differ; bind_rows
# unions the columns (missing cells become NA) and preserves column types
# when every input already shares the same columns (the common case).
.rbindAligned <- function(parts) {
    if (length(parts) == 1L) {
        return(parts[[1L]])
    }
    bind_rows(parts)
}

# Internal: the per-row analysis WINDOW -- the genomic interval each row's fit
# was run over, as stored on the `region` column. Shared by genomicRegion() on
# TwasWeights and FineMappingResultBase.
#
# Stored rather than derived, because it is the half that cannot be recovered:
# QC and screening drop variants before the fit, so the realized span of the
# variants present is narrower than the window that was analysed (a fit over
# chr22:10,000,000-18,000,000 whose variants span 10,414,272-16,916,492). The
# span IS recoverable -- it is just range(x) -- so nothing is lost by deriving
# it on demand, and that is what callers wanting current content should use.
#
# Spec 4.4 removed a stored column on the grounds that it "would quietly go
# stale" under subsetting. That argument holds for a SPAN, which is
# current state and must track subsetting; it does not hold for a WINDOW,
# which is provenance -- subsetting a fit's variants does not retroactively
# change which interval was analysed.
#
# Falls back to the derived span for collections built before the column
# existed (older saved objects, hand-built collections), so the accessor never
# returns nothing where it used to return a range. So does a row a combine
# padded with a placeholder (.naLikeColumn()), having no window of its own.
# @noRd
.getRegionColumn <- function(x) {
    if (nrow(x) == 0L) {
        return(GenomicRanges::GRanges())
    }
    if (!is_in("region", .tupleColumnNames(x))) {
        return(unlist(range(x), use.names = FALSE))
    }
    stored <- .tupleColumn(x, "region")
    padded <- width(stored) == 0L
    if (!any(padded)) {
        return(stored)
    }
    .pickRegions(!padded, stored, .rowSpanRanges(x))
}

# Internal: append the optional `blockId` provenance column (no-op when NULL).
#
# blockId keys the EXTERNAL LD block manifest. It exists because block
# BOUNDARIES are the one thing not recoverable from the variants -- adjacent
# blocks' spans leave gaps, so the realized span of the variants present is
# narrower than the block. Everything else positional comes from range().
.appendBlockIdCol <- function(cols, blockId, n) {
    if (is.null(blockId)) {
        return(cols)
    }
    if (length(blockId) != n) {
        msg <- glue(
            "`blockId` must have the same length as `studyName` ",
            "(got {length(blockId)} vs {n})."
        )
        abort(msg)
    }
    c(cols, list(blockId = as.character(blockId)))
}

# Internal: normalize a caller-supplied window -- a GRanges, as the
# pipeline's `region` is -- into the single range the `region` column
# stores, or NULL when there is none. A multi-range window (jointRegions) is
# recorded as its overall span, which is the interval the fit covered.
.asAnalysisWindow <- function(region) {
    if (is.null(region) || length(region) == 0L) {
        return(NULL)
    }
    gr <- .assertRegionGRanges(region)
    if (length(gr) == 1L) gr else range(gr)[1L]
}

# Internal: `.asAnalysisWindow()` replicated to one range per row, for the
# constructors that take a whole collection's worth of rows at once.
.analysisWindowForRows <- function(region, n) {
    gr <- .asAnalysisWindow(region)
    if (is.null(gr) || n == 0L) {
        return(NULL)
    }
    rep(gr, n)
}

# Internal: an optional per-row provenance column. One helper for all three
# (`dataType`, `dataSource`, `standardized`) since they differ only in name and
# coercion.
#
# `dataSource` is per ROW, not per collection: it used to be derived from the
# collection's single `ldReference` slot, so after combineFineMappingResults()
# folded individual-level and RSS-derived fits together every row got the
# combined object's answer.
#
# `standardized` gets the same treatment fine-mapping needed anyway: the
# posterior effect IS the prediction weight, so it carries the same genotype-
# scale ambiguity as a TWAS weight -- and the fit cannot answer it (susie
# records whether IT scaled, not whether the supplied X was already
# standardized), so it has to be recorded rather than derived.
# @noRd
.appendProvenanceCol <- function(cols, value, nm, n, coerce = as.character) {
    if (is.null(value)) {
        return(cols)
    }
    v <- coerce(value)
    if (!is_in(length(v), c(1L, n))) {
        abort(glue("`{nm}` must be length 1 or the same length as `study`."))
    }
    c(cols, set_names(list(rep(v, length.out = n)), nm))
}

# @noRd
.appendFmProvenanceCols <- function(
    cols,
    dataType,
    dataSource,
    standardized,
    n
) {
    cols <- .appendProvenanceCol(cols, dataType, "dataType", n)
    cols <- .appendProvenanceCol(cols, dataSource, "dataSource", n)
    .appendProvenanceCol(cols, standardized, "standardized", n, as.logical)
}

# Internal: append a validated `region` GRanges to a constructor's column list
# (no-op when NULL). Shared by the TwasWeights / FineMappingResultBase
# constructors so the provenance column is added identically everywhere.
# @noRd
.appendRegionCol <- function(cols, region, n) {
    if (is.null(region)) {
        return(cols)
    }
    if (!methods::is(region, "GRanges")) {
        abort("`region` must be a GRanges (one range per row) or NULL.")
    }
    if (length(region) != n) {
        abort("`region` must have the same length as `study`.")
    }
    c(cols, list(region = .ucscRanges(region)))
}

# Internal: validate the optional `region` GRanges column (one range per row),
# shared by the TwasWeights / FineMappingResultBase validity methods. Returns a
# character vector of errors (empty when valid or the column is absent).
.validateRegionColumn <- function(object) {
    if (!is_in("region", .tupleColumnNames(object))) {
        return(character(0))
    }
    region <- .tupleColumn(object, "region")
    res <- checkClass(region, "GRanges")
    if (!isTRUE(res)) {
        return(str_c("'region' column ", res))
    }
    if (length(region) != nrow(object)) {
        return("'region' column must have one range per row")
    }
    character(0)
}

# Internal: `n` placeholder ranges, zero-width on "chrUn", for the rows of a
# range column that have no range of their own: an empty element, or the rows
# of a combined part that lacked the column. No real range is zero-width, so
# that is what marks one.
# @noRd
.placeholderRanges <- function(n) {
    GenomicRanges::GRanges(
        rep("chrUn", n),
        IRanges::IRanges(start = rep(1L, n), width = rep(0L, n))
    )
}

# Internal: each element's span as one bare range (no seqinfo, strand or
# metadata), or a placeholder for an empty element. An element sits on one
# seqname, so strand aside its span is one range.
# @noRd
.rowSpanRanges <- function(x) {
    spans <- range(x, ignore.strand = TRUE)
    has <- lengths(spans) > 0L
    flat <- unlist(spans, use.names = FALSE)
    n <- length(spans)
    GenomicRanges::GRanges(
        replace(rep("chrUn", n), has, as.character(seqnames(flat))),
        IRanges::IRanges(
            start = replace(rep(1L, n), has, start(flat)),
            end = replace(rep(0L, n), has, end(flat))
        )
    )
}

# Internal: per row, the range of `a` where `use` is TRUE, else that of `b`,
# as bare ranges. `a` and `b` have one range per row.
# @noRd
.pickRegions <- function(use, a, b) {
    GenomicRanges::GRanges(
        if_else(use, as.character(seqnames(a)), as.character(seqnames(b))),
        IRanges::IRanges(
            start = if_else(use, start(a), start(b)),
            end = if_else(use, end(a), end(b))
        )
    )
}

# Internal: a window argument -- a `region`, or h2's `ldBlocks` -- as the
# GRanges of the windows it names. Given in an R call, a window is a GRanges,
# 1-based like everything in Bioconductor; a table or a "chr:start-end"
# string is refused rather than read, so no window is misread by one base. A
# file follows its own standard: the path of a BED file is read as BED (see
# .readBedWindows()).
# @noRd
.asRegionRanges <- function(x, what = "`region`") {
    if (methods::is(x, "GRanges")) {
        return(.ucscRanges(x))
    }
    if (!is.character(x) || length(x) != 1L) {
        abort(c(
            glue(
                "{what} must be a GRanges or the path of a BED file ",
                "(got {class(x)[[1L]]})."
            ),
            i = .regionArgHint()
        ))
    }
    if (!file.exists(x)) {
        abort(c(
            glue("{what}: no file exists at '{x}'."),
            i = .regionArgHint()
        ))
    }
    .ucscRanges(.readBedWindows(x, what))
}

# Internal: a window argument that takes no file -- a pipeline's `region`,
# or a QtlDataset's -- which must be a GRanges. It comes back in UCSC style,
# so everything it is compared with inside is compared exactly.
# @noRd
.assertRegionGRanges <- function(x, what = "`region`") {
    if (!methods::is(x, "GRanges")) {
        abort(c(
            glue("{what} must be a GRanges (got {class(x)[[1L]]})."),
            i = .regionArgHint()
        ))
    }
    .ucscRanges(x)
}

# Internal: how to give a window in an R call, for the errors refusing a
# window argument of the wrong kind.
# @noRd
.regionArgHint <- function() {
    str_c(
        "A window given in R is a GRanges, which is 1-based: for example ",
        "GenomicRanges::GRanges('chr22:1-100')."
    )
}

# Internal: the variants of `g` inside the one window `win`. A variant sits at
# its start.
# @noRd
.withinWindow <- function(g, win) {
    same <- as.character(seqnames(g)) == as.character(seqnames(win))[[1L]]
    g[same & start(g) >= start(win)[[1L]] & start(g) <= end(win)[[1L]]]
}

# Internal: `gr` as bare ranges, its seqnames and intervals alone.
# @noRd
.bareRanges <- function(gr) {
    GenomicRanges::GRanges(
        as.character(seqnames(gr)),
        IRanges::IRanges(start = start(gr), end = end(gr))
    )
}

# Internal: the `region` and `regionGiven` columns of a GWAS collection. A row
# whose block range was given keeps it, and every other row takes the span
# of its own variants. `given` holds one range per row and `has` says which
# of them were given.
# @noRd
.regionColumns <- function(given, has, spans) {
    list(region = .pickRegions(has, given, spans), regionGiven = has)
}

# Internal: refuse a window a fine-mapping result was given that leaves out
# its own row's variants. The check on QTL variants (.checkGwasBlockReach())
# takes a given window as the block's true boundary, so one that does not hold
# its own variants would send it after a block that is not missing. A variant
# sits at its start, so an indel's width cannot carry it past a window end.
# `pieces` holds each row's variants, and `has` says which rows were given
# `given`.
# @noRd
.refuseUncontainedRegions <- function(given, has, pieces, caller) {
    check <- has & lengths(pieces) > 0L
    if (!any(check)) {
        return(invisible(NULL))
    }
    starts <- start(pieces[check])
    chrom <- .rtlElementSeqnames(pieces[check])
    lo <- min(starts)
    hi <- max(starts)
    inside <- as.character(seqnames(given[check])) == chrom &
        start(given[check]) <= lo &
        hi <= end(given[check])
    if (all(inside)) {
        return(invisible(NULL))
    }
    out <- !inside
    rows <- glue(
        "row {which(check)[out]}: range {as.character(given[check][out])}, ",
        "variants at {chrom[out]}:{lo[out]}-{hi[out]}"
    )
    abort(c(
        glue(
            "{caller}: each row's `region` must contain the row's variants, ",
            "since it is the window the row was fit over."
        ),
        set_names(utils::head(rows, 5L), rep("*", min(5L, length(rows))))
    ))
}

# Internal: each row's analysis window and whether the user gave it, from the
# `region` / `regionGiven` columns of a summary-statistic collection or a
# fine-mapping result. A collection saved before the columns existed has no
# stored window, and neither has a row a combine padded; such a row falls back
# to its span and counts as computed.
# @noRd
.rowRegions <- function(x) {
    spans <- .rowSpanRanges(x)
    cols <- .tupleColumnNames(x)
    if (!is_in("region", cols)) {
        return(list(region = spans, regionGiven = rep(FALSE, length(spans))))
    }
    stored <- .tupleColumn(x, "region")
    has <- width(stored) > 0L
    given <- if (is_in("regionGiven", cols)) {
        as.logical(.tupleColumn(x, "regionGiven"))
    } else {
        rep(FALSE, length(spans))
    }
    list(
        region = .pickRegions(has, stored, spans),
        regionGiven = has & !is.na(given) & given
    )
}

# Internal: row `i`'s analysis window alone (see .rowRegions()), read off its
# own column entry rather than every row's, for per-row lookups. NULL for an
# empty row with no window stored.
# @noRd
.rowRegionAt <- function(x, i) {
    if (is_in("region", .tupleColumnNames(x))) {
        stored <- .tupleColumn(x, "region")[i]
        if (width(stored) > 0L) {
            return(.bareRanges(stored))
        }
    }
    g <- x[[i]]
    if (length(g) == 0L) {
        return(NULL)
    }
    .bareRanges(range(g, ignore.strand = TRUE)[1L])
}

# Internal: set a GWAS collection's `region` / `regionGiven` columns, for the
# rebuilds that know them better than the constructor can: a QC pass, a union
# of blocks, a fine-mapping run, a file read back.
# @noRd
.withRegionColumns <- function(x, cols) {
    md <- mcols(x, use.names = FALSE)
    md$region <- cols$region
    md$regionGiven <- as.logical(cols$regionGiven)
    S4Vectors::`mcols<-`(x, value = md)
}

# ---- QTL variants past the GWAS blocks --------------------------------------
# A QTL region is scored against every block of a GWAS study it overlaps,
# merged into one. Its variants outside all of the study's loaded blocks are
# left out of that merge, which is a sign a block it reaches was not loaded.
# Past a block range the caller gave, they do lie in a block that was not
# loaded, so that is an error. Past a range computed from the GWAS variants
# they may only be variants the GWAS lacks, at a chromosome end or untyped,
# so that is a warning, and only for more than one variant in a row.

# Internal: check `units` -- a named list of variant positions (GRanges), one
# per region the caller scores against `gws`, named as a message should call
# it -- against the block ranges of each GWAS fit in `gws`. Each variant past
# a fit's blocks counts against its nearest block, on the side it falls.
# `caller` and `what` (the units, as a plural noun) word the messages.
# @noRd
.checkGwasBlockReach <- function(units, gws, caller, what) {
    if (length(units) == 0L || nrow(gws) == 0L) {
        return(invisible(NULL))
    }
    blocks <- .gwasBlockTable(gws)
    issues <- list_rbind(imap(units, .unitBlockReach, blocks = blocks))
    .reportBlockReach(issues, caller, what)
}

# Internal: the block ranges of `gws`, one per row with a range, each with
# the GWAS fit it belongs to (its study, and method for a fine-mapping
# result) and a label naming the block.
# @noRd
.gwasBlockTable <- function(gws) {
    regions <- .rowRegions(gws)
    keep <- width(regions$region) > 0L
    region <- regions$region[keep]
    chrom <- as.character(seqnames(region))
    ids <- as.character(.tupleColumn(gws, "blockId"))[keep]
    key <- str_c(chrom, start(region), end(region), sep = "_")
    span <- str_c(chrom, ":", start(region), "-", end(region))
    list(
        ranges = .locusRanges(region),
        given = regions$regionGiven[keep],
        fit = .gwasFitLabels(gws)[keep],
        # An id that only spells the range is not repeated beside it.
        label = if_else(
            is.na(ids) | ids == key,
            span,
            str_c(ids, " (", span, ")")
        )
    )
}

# @noRd
.gwasFitLabels <- function(gws) {
    study <- str_c("'", as.character(gws$study), "'")
    if (!is_in("method", .tupleColumnNames(gws))) {
        return(study)
    }
    str_c(study, " (", as.character(gws$method), ")")
}

# Internal: one unit's variants past each GWAS fit whose blocks they meet.
# @noRd
.unitBlockReach <- function(variants, unit, blocks) {
    v <- .locusRanges(variants)
    list_rbind(map(
        unique(blocks$fit),
        .unitFitReach,
        v = v,
        unit = unit,
        blocks = blocks
    ))
}

# @noRd
#' @importFrom dplyr count
.unitFitReach <- function(fit, v, unit, blocks) {
    own <- which(blocks$fit == fit)
    r <- blocks$ranges[own]
    inside <- .overlapsOnShared(v, r)
    if (!any(inside)) {
        return(NULL)
    }
    chroms <- unique(as.character(seqnames(v))[inside])
    out <- v[!inside & is_in(as.character(seqnames(v)), chroms)]
    if (length(out) == 0L) {
        return(NULL)
    }
    onChroms <- which(is_in(as.character(seqnames(r)), chroms))
    hit <- own[onChroms][GenomicRanges::nearest(
        .onSeqlevels(out, chroms),
        .onSeqlevels(r[onChroms], chroms),
        ignore.strand = TRUE
    )]
    count(
        tibble(
            unit = unit,
            fit = fit,
            block = blocks$label[hit],
            given = blocks$given[hit],
            side = if_else(
                start(out) < start(blocks$ranges[hit]),
                "start",
                "end"
            )
        ),
        .data$unit,
        .data$fit,
        .data$block,
        .data$given,
        .data$side
    )
}

# Internal: abort over variants past a given range, else warn over runs of
# more than one variant past a computed one.
# @noRd
.reportBlockReach <- function(issues, caller, what) {
    if (is.null(issues) || nrow(issues) == 0L) {
        return(invisible(NULL))
    }
    given <- filter(issues, .data$given)
    if (nrow(given) > 0L) {
        abort(c(
            glue(
                "{caller}: {what} reach past the LD blocks loaded for a GWAS ",
                "study they overlap. Those blocks' ranges were given, so the ",
                "variants past them lie in blocks that were not loaded."
            ),
            .blockReachLines(given),
            i = "Load the neighbouring LD blocks of those studies too."
        ))
    }
    computed <- filter(issues, !.data$given, .data$n > 1L)
    if (nrow(computed) > 0L) {
        warn(c(
            glue(
                "{caller}: {what} reach past the variants of a GWAS study ",
                "they overlap. Those blocks' ranges are the span of the GWAS ",
                "variants, so the variants past them may only be ones the ",
                "GWAS lacks, or a block may not have been loaded."
            ),
            .blockReachLines(computed),
            i = str_c(
                "Give the GWAS its windows as a `region` to have this ",
                "checked strictly."
            )
        ))
    }
    invisible(NULL)
}

# @noRd
.blockReachLines <- function(issues) {
    shown <- utils::head(issues, 5L)
    lines <- glue(
        "{shown$unit}: {shown$n} variant(s) past the {shown$side} of block ",
        "{shown$block} of GWAS study {shown$fit}."
    )
    more <- nrow(issues) - nrow(shown)
    set_names(
        c(lines, if (more > 0L) glue("... and {more} more.")),
        rep("*", length(lines) + (more > 0L))
    )
}

# Internal: trait-position provenance. One GRanges per trait carrying the
# molecular feature's OWN genomic coordinates (gene/peak range; TSS = start()),
# distinct from the fine-mapping `region`. The true trait position cannot be
# inferred from summary statistics, so it is threaded from QtlDataset rowRanges
# / an explicit QtlSumStats trait-position and carried onto QtlFineMappingResult
# and TwasWeights. Mirrors the `region` helpers above; absent for GWAS.
# traitPos is optional provenance: always known for a QtlDataset, but only known
# for a QtlSumStats when the caller supplied it (it cannot be inferred from
# summary statistics). When the column is absent the answer is an empty
# GRanges -- no position, in the type every other answer has.
.getTraitPosColumn <- function(x) {
    if (is_in("traitPos", .tupleColumnNames(x))) {
        .tupleColumn(x, "traitPos")
    } else {
        GenomicRanges::GRanges()
    }
}

.appendTraitPosCol <- function(cols, traitPos, n) {
    if (is.null(traitPos)) {
        return(cols)
    }
    if (!methods::is(traitPos, "GRanges")) {
        abort("`traitPos` must be a GRanges (one range per row) or NULL.")
    }
    if (length(traitPos) != n) {
        abort("`traitPos` must have the same length as `study`.")
    }
    c(cols, list(traitPos = .ucscRanges(traitPos)))
}

#' @importFrom checkmate checkClass
.validateTraitPosColumn <- function(object) {
    if (!is_in("traitPos", .tupleColumnNames(object))) {
        return(character(0))
    }
    traitPos <- .tupleColumn(object, "traitPos")
    res <- checkClass(traitPos, "GRanges")
    if (!isTRUE(res)) {
        return(str_c("'traitPos' column ", res))
    }
    if (length(traitPos) != nrow(object)) {
        return("'traitPos' column must have one range per row")
    }
    character(0)
}

# Internal: a length-n "unknown" column matching the type of `exemplar`, used
# to pad a collection that lacks an optional column before a union row-bind.
# Character -> NA; GRanges -> placeholder ranges; else NA.
.naLikeColumn <- function(exemplar, n) {
    if (methods::is(exemplar, "GRanges")) {
        return(.placeholderRanges(n))
    }
    if (methods::is(exemplar, "CharacterList")) {
        return(IRanges::CharacterList(rep(list(character(0)), n)))
    }
    if (methods::is(exemplar, "List")) {
        return(S4Vectors::SimpleList(rep(list(NULL), n)))
    }
    if (is.character(exemplar)) {
        return(rep(NA_character_, n))
    }
    rep(NA, n)
}

# First existing value for column `cn` across `parts` (a type exemplar used for
# NA-filling); NULL if no part has it.
# @noRd
.exemplarColumn <- function(cn, parts) {
    for (p in parts) {
        if (is_in(cn, .tupleColumnNames(p))) {
            return(.tupleColumn(p, cn))
        }
    }
    # Callers always draw `cn` from allCols, so in practice some part has it;
    # the NULL is the honest answer when none does.
    NULL
}

# Internal: row-bind two or more per-tuple collection objects (TwasWeights /
# FineMappingResultBase subclasses), carrying forward EVERY column. Columns are
# unioned, and a collection lacking an optional column (e.g. jointContexts or
# region) is padded with `.naLikeColumn` -- so adding a column to a class flows
# through combines automatically rather than being hand-listed at each site.
# Preserves the concrete class and rebuilds the LD references: every row keeps
# its own (.ldRefCombine()), unless `ldReference` names one for all of them.
# Returns NULL when given no inputs. Every collection is a RangedTupleList now,
# so this only drops the NULL parts and hands off; the old DFrame path went
# with the rebase onto GRangesList.
#
# `slots` carries the collection-level state a concrete subclass adds beyond
# `ldReference` (SumStatsBase's `genome` / `qcOptions`). It is passed in rather
# than read off `parts[[1L]]` because those slots describe the WHOLE collection,
# so merging them is the caller's decision -- first-wins would silently drop the
# other parts' QC audit.
.rbindCollections <- function(
    parts,
    ldReference = NULL,
    slots = list(),
    genome = NULL,
    fn = ".rbindCollections"
) {
    parts <- compact(parts)
    if (length(parts) == 0L) {
        return(NULL)
    }
    dict <- .ldRefCombine(parts, ldReference, fn)
    .rbindRangedCollections(parts, dict, slots, genome)
}

# Internal: aggregate a per-entry accessor across every row of a
# FineMappingResultBase collection that matches the given selectors, prefixing
# each entry's rows with the row identity (study/context/trait/blockId/method)
# so rows stay attributable to their source entry. Shared by topLoci /
# credibleSets
# / marginalEffects on FineMappingResultBase.
#
# When at least one selector is supplied AND resolves to exactly one entry, this
# short-circuits to `onSingle(entry)` -- the bare per-entry view with no
# identity columns -- preserving the single-entry contract. A no-selector call
# always
# aggregates (even a one-row collection) so the output shape is uniform.
#
#   perEntry(entry) -> data.frame : per-entry view used when aggregating. Any
#     column whose name collides with an identity column is dropped before the
#     metadata prefix is attached (so an added `method` never duplicates).
#   onSingle(entry) -> any        : returned verbatim when selectors pin one
#     entry (defaults to perEntry).
.fmrAggregateView <- function(
    x,
    studyName = NULL,
    context = NULL,
    trait = NULL,
    method = NULL,
    blockId = NULL,
    perEntry,
    onSingle = perEntry,
    viewArgs = list()
) {
    sel <- list(
        studyName = studyName,
        context = context,
        trait = trait,
        method = method,
        blockId = blockId
    )
    single <- .fmrTrySingle(x, sel, onSingle, viewArgs)
    if (single$done) {
        return(single$result)
    }
    idx <- exec(.fmrRowsMatching, x, !!!sel)
    meta <- .fmrRowMetadata(x)
    if (length(idx) == 0L) {
        return(.fmrEmptyOrError(single, meta))
    }
    parts <- compact(map(
        idx,
        .fmrEntryPart,
        x = x,
        perEntry = perEntry,
        meta = meta,
        viewArgs = viewArgs
    ))
    if (length(parts) == 0L) {
        return(slice(meta, integer(0)))
    }
    .rbindAligned(parts)
}

# Empty-match handling for an aggregate view: an explicit selector that matched
# nothing re-raises its informative selection error; otherwise return an empty
# identity-only frame.
# @noRd
.fmrEmptyOrError <- function(single, meta) {
    if (inherits(single$sel, "error")) {
        cnd_signal(single$sel)
    }
    slice(meta, integer(0))
}

# Attempt the single-entry fast path when any of `selectors` -- the named row
# selectors -- is supplied. Returns list(done, result) on a successful
# single-select, else list(done = FALSE, sel) carrying the selection error
# (or NULL when no selector was given).
# @noRd
#' @importFrom rlang try_fetch
.fmrTrySingle <- function(x, selectors, onSingle, viewArgs) {
    if (length(compact(selectors)) == 0L) {
        return(list(done = FALSE, sel = NULL))
    }
    sel <- try_fetch(
        exec(.fmrSelectEntry, x, !!!selectors),
        error = function(cnd) cnd
    )
    if (!inherits(sel, "error")) {
        return(list(done = TRUE, result = exec(onSingle, sel, !!!viewArgs)))
    }
    list(done = FALSE, sel = sel)
}

# Project one entry to its per-entry view with the row's identity metadata
# prepended; NULL when the entry yields no rows.
# @noRd
.fmrEntryPart <- function(i, x, perEntry, meta, viewArgs) {
    # `perEntry` is an internal per-row function taking the stored payload, not
    # an S4 generic dispatching on a rebuilt entry. That indirection is what
    # the retired FineMappingEntry class existed for.
    raw <- exec(perEntry, .fmrRowParts(x, i), !!!viewArgs)
    if (is.null(raw) || nrow(raw) == 0L) {
        return(NULL)
    }
    v <- select(raw, all_of(setdiff(names(raw), names(meta))))
    bind_cols(slice(meta, rep(i, nrow(v))), v)
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# One identity-metadata column `cc` of a collection (NA-filled when absent).
# @noRd
.fmrMetadataCol <- function(cc, x, n) {
    if (is_in(cc, .tupleColumnNames(x))) {
        as.character(.tupleColumn(x, cc))
    } else {
        rep(NA_character_, n)
    }
}

# Column `cn` of part `p`, or a NA-like column of an exemplar type when absent.
# @noRd
.rbindColPiece <- function(p, cn, parts) {
    if (is_in(cn, .tupleColumnNames(p))) {
        .tupleColumn(p, cn)
    } else {
        .naLikeColumn(.exemplarColumn(cn, parts), nrow(p))
    }
}

# Row-binding a collection: the elements append and the mcols rbind.
#
# S4Vectors::combineRows() does the union-with-NA-fill this needs, but it
# fails on a GRanges-valued mcols column (`traitPos`) with "GRanges objects
# don't support [[, as.list(), lapply(), or unlist()". .rbindColumn() binds
# each column with c() instead, which GRanges does support, so the union is
# done per column here rather than delegated.
# @noRd
.rbindRangedCollections <- function(
    parts,
    dict,
    slots = list(),
    genome = NULL
) {
    cls <- class(parts[[1L]])[[1L]]
    allCols <- reduce(map(parts, .tupleColumnNames), union)
    # Each part keyed its rows into its own references, so the bound rows take
    # the keys the combined references were built with.
    combined <- list_assign(
        set_names(map(allCols, .rbindColumn, parts = parts), allCols),
        ldReference = dict$keys
    )
    md <- exec(
        S4Vectors::DataFrame,
        !!!c(combined, list(check.names = FALSE))
    )
    # Joining the parts merges each part's seqinfo, so a part whose elements
    # never had a build set would leave the result carrying both the real
    # build and NA. The agreed build is written back explicitly to keep
    # seqinfo single-valued, which validity requires.
    built <- .withGenomeBuild(.rtlJoinParts(parts), !is.null(genome), genome)
    grl <- .rtlNamed(S4Vectors::`mcols<-`(built, value = md))
    exec(methods::new, cls, grl, !!!c(list(ldReference = dict$refs), slots))
}

# Concatenate column `cn` across all parts (NA-filling parts that lack it).
# @noRd
.rbindColumn <- function(cn, parts) {
    pieces <- map(parts, .rbindColPiece, cn = cn, parts = parts)
    # GRanges concat may warn on seqinfo.
    suppressWarnings(exec(c, !!!pieces))
}

# =============================================================================
# Per-row primitives
# -----------------------------------------------------------------------------
# The per-row views a collection projects. These replace the FineMappingEntry
# S4 class, which existed only to be the DISPATCH TARGET for exactly this work:
# a collection method would rebuild an entry for each row and call the exported
# generic on it. That indirection cost a whole class, a validity method and a
# hand-written alignment check (.fmeCheckTopLociOrder) whose only reason to
# exist was that the rebuild re-split the element into `variantIds` and
# `topLoci` slots and could get them out of step.
#
# Reading the stored form directly removes all of it. `parts` is a plain list,
# not a class: it adds no public surface, no dispatch and nothing to validate.
# =============================================================================

# The stored payload of one row: the element (variants, with topLoci as its
# inner mcols) plus the outer-mcols fit payload.
# @noRd
.fmrRowParts <- function(x, idx) {
    md <- mcols(x, use.names = FALSE)
    # `idx` may name a group of rows that form one logical entry (a
    # multi-block sweep). Their elements concatenate; the fit payload comes
    # from the first, matching what the entry rebuild did.
    first <- idx[[1L]]
    new(
        "FineMappingRow",
        variants = .rtlGatherElements(x, idx),
        susieFit = md$susieFit[[first]],
        cvResult = md$cvResult[[first]],
        fitVariants = .mdFitVariants(md, idx)
    )
}

# The variants the fit of rows `idx` was made with, in its order: the first
# row's, as the fit payload is, since rows split from one entry repeat both.
# None for a collection saved before rows recorded them.
# @noRd
.mdFitVariants <- function(md, idx) {
    fv <- md$fitVariants
    if (is.null(fv)) {
        return(character(0))
    }
    as.character(fv[[idx[[1L]]]])
}

# @noRd
.fmrPartsVariantIds <- function(parts) {
    variantIds(.asFmRowPayload(parts))
}

# @noRd
.fmrPartsSusieFit <- function(parts) {
    susieFit(.asFmRowPayload(parts))
}

# @noRd
.fmrPartsCvResult <- function(parts) {
    cvResult(.asFmRowPayload(parts))
}

# @noRd
.fmrPartsTopLoci <- function(parts) {
    .fmeTopLociFromElement(variants(.asFmRowPayload(parts)))
}

# ---- per-row views ----------------------------------------------------------
# Each of these is the body of what used to be a FineMappingEntry method,
# reading the stored form through `parts` instead of through S4 slots. They are
# what .fmrAggregateView() maps over.

# @noRd
.fmrRowTopLoci <- function(
    parts,
    signalCutoff = 0.025,
    minPurity = NULL,
    raw = FALSE
) {
    tl <- .fmrPartsTopLoci(parts)
    # raw = TRUE returns the stored canonical table verbatim: every variant,
    # every column, no posterior-view projection.
    if (isTRUE(raw)) {
        return(tl)
    }
    # Always a table: .fmrAggregateView() stacks these, and topLoci()
    # converts the stacked result to ranges once.
    .fmeFilterTopLoci(tl, signalCutoff, minPurity)
}

# @noRd
.fmrRowPip <- function(parts) {
    tl <- .fmrPartsTopLoci(parts)
    if (nrow(tl) == 0L || !is_in("pip", names(tl))) {
        return(numeric(0))
    }
    set_names(tl$pip, tl$variant_id)
}

# @noRd
.fmrRowMarginalEffects <- function(parts) {
    .projectMarginalView(.fmrPartsTopLoci(parts))
}

# @noRd
.fmrRowCs <- function(parts, coverage = 0.95, minPurity = NULL) {
    tl <- .fmrPartsTopLoci(parts)
    if (nrow(tl) == 0L) {
        return(.projectPosteriorView(tl))
    }
    csCol <- names(tl)[str_detect(
        names(tl),
        str_c("^cs_", coverage * 100, "$")
    )]
    if (length(csCol) == 0L) {
        return(.projectPosteriorView(slice(tl, 0)))
    }
    inCs <- !is.na(tl[[csCol[1L]]]) &
        str_length(tl[[csCol[1L]]]) > 0L &
        !str_detect(tl[[csCol[1L]]], "_0$")
    # Independent purity filter (min.abs.corr), orthogonal to `coverage`: drop
    # CS members whose credible set at THIS coverage is below `minPurity`.
    purCol <- str_c(csCol[1L], "_purity")
    hasPurity <- is_in(purCol, names(tl))
    if (!is.null(minPurity) && !hasPurity) {
        msg <- glue(
            "credibleSets: no purity column '{purCol}' for coverage ",
            "{coverage}; minPurity filter skipped."
        )
        warn(msg)
    }
    keep <- if (is.null(minPurity) || !hasPurity) {
        inCs
    } else {
        pur <- as.numeric(tl[[purCol]])
        inCs & !is.na(pur) & pur >= minPurity
    }
    .projectPosteriorView(tl[keep, , drop = FALSE])
}

# @noRd
# The fit's own log Bayes factors of the variants the row still holds, each
# found in the fit by id; a row narrowed since the fit loses the columns of
# the variants it dropped, and nothing is rescaled.
.fmrRowLbf <- function(parts) {
    lbf <- .asLbfMatrix(susieFit(parts))
    if (is.null(lbf)) {
        return(tibble(variant_id = character(0)))
    }
    fitted <- .partsFitIds(parts, lbf)
    if (ncol(lbf) != length(fitted)) {
        return(tibble(variant_id = character(0)))
    }
    vids <- .fmrPartsVariantIds(parts)
    col <- .fitColumnOf(vids, fitted)
    held <- !is.na(col)
    wide <- as_tibble(
        t(as.matrix(lbf)[, col[held], drop = FALSE]),
        .name_repair = "minimal"
    )
    w <- `names<-`(wide, str_c("lbf_L", seq_len(ncol(wide))))
    bind_cols(tibble(variant_id = as.character(vids[held])), w)
}

# @noRd
.fmrRowCredibleSetSummary <- function(parts, coverage = 0.95) {
    .csSummaryFit(.fmrPartsTopLoci(parts), susieFit(parts), coverage)
}

# @noRd
.fmrRowFsusieCredibleBand <- function(parts) {
    .fsusieCredibleBandFit(susieFit(parts))
}

# @noRd
.fmrRowFsusieAffectedRegions <- function(parts) {
    .fsusieAffectedRegionsFit(
        susieFit(parts),
        topLoci = .fmrPartsTopLoci(parts)
    )
}

# @noRd
.fmrRowResolveWeights <- function(parts) {
    empty <- list(variantIds = character(0), weights = numeric(0))
    # The topLoci posterior view projects the effect to `posterior_effect`,
    # and that effect IS the per-variant weight -- a susie posterior mean is
    # the TWAS coefficient. Aligned with variant_id.
    tl <- .fmrRowTopLoci(parts)
    if (
        is.null(tl) ||
            nrow(tl) == 0L ||
            !all(is_in(c("variant_id", "posterior_effect"), names(tl)))
    ) {
        return(empty)
    }
    vids <- as.character(tl$variant_id)
    w <- as.numeric(tl$posterior_effect)
    ok <- !is.na(vids) & !is.na(w)
    if (!any(ok)) {
        return(empty)
    }
    list(variantIds = vids[ok], weights = w[ok])
}

# The stored payload of one TwasWeights row. Mirrors .fmrRowParts: weights are
# the element's inner mcols, the rest is outer mcols.
# @noRd
.twrRowParts <- function(x, idx) {
    md <- mcols(x, use.names = FALSE)
    first <- idx[[1L]]
    gr <- .rtlGatherElements(x, idx)
    new(
        "TwasWeightsRow",
        variants = gr,
        weights = mcols(gr, use.names = FALSE)$weight,
        methodFits = md$fits[[first]],
        cvResult = md$cvResult[[first]],
        weightStandardized = isTRUE(md$standardized[[first]]),
        weightsDataType = md$dataType[[first]],
        fitVariants = .mdFitVariants(md, idx)
    )
}

# @noRd
.twrPartsVariantIds <- function(parts) {
    variantIds(parts)
}

# The weight vector aligned to the row's variant ids -- what resolveWeights
# produced from a TwasWeightsRow.
# @noRd
.twrRowResolveWeights <- function(parts) {
    empty <- list(variantIds = character(0), weights = numeric(0))
    vids <- .twrPartsVariantIds(parts)
    raw <- weights(parts)
    if (length(vids) == 0L || is.null(raw)) {
        return(empty)
    }
    w <- as.numeric(raw)
    if (length(w) != length(vids)) {
        return(empty)
    }
    ok <- !is.na(vids) & !is.na(w)
    if (!any(ok)) {
        return(empty)
    }
    list(variantIds = vids[ok], weights = w[ok])
}

# ---- fit coverage ------------------------------------------------------------
# A row keeps the fit it was built with when its variants are narrowed (to
# those a GWAS also carries, say). The fit, and everything read off it, still
# describes the run on all of its variants, so reading it from a row that
# holds fewer warns. Nothing is rescaled: intersectVariants() is the explicit
# way to restrict a fit.

# The variants a row's stored fit was made with.
# @noRd
.rowFitVariants <- function(parts) {
    if (methods::is(parts, "TwasWeightsRow")) {
        return(.twRowFitVariants(parts))
    }
    .fmRowFitVariants(parts)
}

# The fit's variables in its own order: the ids the row recorded when it was
# built, else the names the fit's matrix `mat` carries, else -- an unnamed fit
# of a row saved before rows recorded them -- the row's variants as they are.
# @noRd
.partsFitIds <- function(parts, mat = NULL) {
    recorded <- .rowFitVariants(parts)
    if (length(recorded) > 0L) {
        return(recorded)
    }
    colnames(mat) %||% variantIds(parts)
}

# How many of `fitted` are gone from `held`. A variant whose alleles
# harmonization swapped is held under another id, so ids that do not match
# exactly are matched allele-aware before they count as gone.
# @noRd
.lostFitVariantCount <- function(fitted, held) {
    gone <- fitted[!is_in(fitted, held)]
    if (length(gone) == 0L) {
        return(0L)
    }
    matched <- matchVariants(gone, held, removeStrandAmbiguous = FALSE)
    length(gone) - length(unique(matched$idxA))
}

# For each of `held`, its column in a fit made on `fitted` (NA when the fit
# has none), matched the same allele-aware way.
# @noRd
.fitColumnOf <- function(held, fitted) {
    col <- match(held, fitted)
    miss <- which(is.na(col))
    if (length(miss) == 0L) {
        return(col)
    }
    m <- matchVariants(held[miss], fitted, removeStrandAmbiguous = FALSE)
    replace(col, miss[m$idxA], m$idxB)
}

# `parts` unchanged, after warning when its row holds fewer variants than its
# stored fit was made with.
# @noRd
.fitCoverageChecked <- function(parts, caller) {
    fitted <- .rowFitVariants(parts)
    lost <- .lostFitVariantCount(fitted, variantIds(parts))
    if (lost > 0L) {
        warn(
            c(
                glue(
                    "{caller}: the stored fit was made with ",
                    "{length(fitted)} variants, {lost} of which this row no ",
                    "longer holds."
                ),
                i = .fitCoverageHint(parts)
            ),
            class = "pecotmrFitCoverage"
        )
    }
    parts
}

# `x` unchanged, after warning once about every row that holds fewer variants
# than its stored fit was made with.
# @noRd
.fitCoverageCheckedRows <- function(x, caller) {
    .fitCoverageCheckedAll(list(x), caller)
    invisible(x)
}

# Warn once about the narrowed rows of every collection in `xs`, so a
# pipeline reading the fits of two inputs names all of them in one warning.
# @noRd
.fitCoverageCheckedAll <- function(xs, caller) {
    labels <- map(xs, .narrowedRowLabels)
    rows <- list_c(labels, ptype = character())
    if (length(rows) == 0L) {
        return(invisible(NULL))
    }
    warn(
        c(
            glue(
                "{caller}: {length(rows)} row(s) hold fewer variants than ",
                "their stored fit was made with."
            ),
            set_names(utils::head(rows, 5L), rep("*", min(5L, length(rows)))),
            i = .fitCoverageHint(xs[[which(lengths(labels) > 0L)[[1L]]]])
        ),
        class = "pecotmrFitCoverage"
    )
    invisible(NULL)
}

# The labels of the rows of `x` that hold fewer variants than their stored
# fit was made with. None when `x` records no fit variants, or is not a
# collection at all: an optional input left NULL, or one the caller has yet
# to validate, which this advisory check must not pre-empt.
# @noRd
.narrowedRowLabels <- function(x) {
    if (!methods::is(x, "RangedTupleList")) {
        return(character(0))
    }
    fv <- mcols(x, use.names = FALSE)$fitVariants
    if (is.null(fv) || length(x) == 0L) {
        return(character(0))
    }
    lost <- map_int(seq_along(x), .rowLostFitVariants, x = x, fv = fv)
    short <- which(lost > 0L)
    if (length(short) == 0L) {
        return(character(0))
    }
    .pubRowLabels(x, short)
}

# @noRd
.rowLostFitVariants <- function(i, x, fv) {
    .lostFitVariantCount(as.character(fv[[i]]), .grVariantIds(x[[i]]))
}

# `expr` with its fit-coverage warnings muffled, for a caller that has
# already warned about the rows `expr` reads -- a pipeline that runs itself
# once per study, say.
# @noRd
.muffleFitCoverage <- function(expr) {
    withCallingHandlers(
        expr,
        pecotmrFitCoverage = function(cnd) invokeRestart("muffleWarning")
    )
}

# What to do about a narrowed row of `x`. Only fine-mapping results have a
# way to restrict their fits, intersectVariants(); weights are refitted.
# @noRd
.fitCoverageHint <- function(x) {
    restrict <- if (.isFineMappingPayload(x)) {
        str_c(
            "intersectVariants() restricts a fit to the variants two results ",
            "share; otherwise refit on the variants kept."
        )
    } else {
        "Refit on the variants kept to have fits that describe them."
    }
    str_c(
        "The fit, and the values read from it, describe the original run on ",
        "all of its variants, unchanged. ",
        restrict
    )
}

# @noRd
.isFineMappingPayload <- function(x) {
    methods::is(x, "FineMappingRow") ||
        methods::is(x, "FineMappingResultBase")
}

# The per-variant weight vector of one row, whichever weight source it came
# from.
# @noRd
.rowResolveWeights <- function(parts) {
    if (methods::is(parts, "TwasWeightsRow")) {
        return(.twrRowResolveWeights(parts))
    }
    .fmrRowResolveWeights(parts)
}

# The cross-validated / per-method fits of one row, or NULL for a
# fine-mapping row.
#
# A fine-mapping row deliberately reports NULL: cTWAS would otherwise
# renormalize alpha into an UNstandardized weight, which is inconsistent with
# the standardized posterior effect the topLoci view already carries. Preserved
# verbatim from the retired FineMappingEntry methodFits method.
# @noRd
.rowFits <- function(parts) {
    if (methods::is(parts, "TwasWeightsRow")) methodFits(parts) else NULL
}

# Whether the row's weights are already on the standardized scale.
#
# TRUE for a fine-mapping row: its posterior effect is colSums(alpha * mu),
# which does not divide by the column scale factors, so cTWAS's per-variant
# variance scaling would double-standardize it.
# @noRd
.rowStandardized <- function(parts) {
    if (methods::is(parts, "TwasWeightsRow")) {
        return(weightStandardized(parts))
    }
    TRUE
}

# ---- row-payload normalizers ------------------------------------------------
# The per-row views read a FineMappingRow / TwasWeightsRow. Callers may also
# hold a single-row collection -- what fineMappingResult() returns -- so these
# accept either and hand back the row.

# A FineMappingRow from a row or from a fine-mapping collection's rows.
# @noRd
.asFmRowPayload <- function(e) {
    # A single-row collection is the other shape callers hand over -- it is
    # what fineMappingResult() returns, so anything accepting "one row's
    # worth of fine-mapping" has to take it too.
    if (methods::is(e, "FineMappingResultBase")) {
        if (length(e) == 0L) {
            return(e)
        }
        return(.fmrRowParts(e, seq_len(length(e))))
    }
    e
}

# @noRd
.asTwRowPayload <- function(e) {
    if (methods::is(e, "TwasWeights")) {
        if (length(e) == 0L) {
            return(e)
        }
        return(.twrRowParts(e, seq_len(length(e))))
    }
    e
}

# Every element must be a row payload once normalization has run. Checked
# before the payload is unpacked, so a wrong type reports itself rather than
# surfacing later as `$ operator is invalid for atomic vectors`.
# @noRd
.checkRowPayloads <- function(entry, cls, what) {
    bad <- which(!map_lgl(entry, methods::is, cls))
    if (length(bad) == 0L) {
        return(invisible(NULL))
    }
    msg <- glue(
        "every element of `entry` must be a {what} row ({cls}); element(s) ",
        "{str_flatten(utils::head(bad, 5L), ', ')} are not."
    )
    abort(msg)
}

# The `entry` elements whose variants lie on more than one chromosome. A
# collection element spans exactly one, so such an entry would be split into
# one row per chromosome, each repeating the entry's identity and payload.
# @noRd
.multiChromosomeEntries <- function(entry) {
    which(map_int(entry, .fmRowSeqnameCount) > 1L)
}

# @noRd
.fmRowSeqnameCount <- function(e) {
    n_distinct(as.character(seqnames(variants(e))))
}

# Shape-tolerant readers for values that may still be a legacy entry object
# while the flip is in progress. Delete the entry branch (and these) once
# nothing produces entries.
# @noRd
.rowCvResult <- function(e) {
    cvResult(e)
}

# @noRd
.rowWeights <- function(e) {
    weights(e)
}

# @noRd
.rowVariantIds <- function(e) {
    variantIds(e)
}

# The non-variant half of each TWAS row payload, as outer mcols columns.
# Objects of arbitrary shape (fits, cvResult, dataType) ride in SimpleLists.
# @noRd
.twRowPayloadCols <- function(entry) {
    list(
        fits = S4Vectors::SimpleList(map(entry, methodFits)),
        cvResult = S4Vectors::SimpleList(map(entry, cvResult)),
        standardized = map_lgl(entry, .twPayloadStandardized),
        dataType = S4Vectors::SimpleList(map(entry, weightsDataType)),
        fitVariants = IRanges::CharacterList(map(entry, .twRowFitVariants))
    )
}

# =============================================================================
# Combining summary-statistics collections
# -----------------------------------------------------------------------------
# A SumStatsBase subclass carries three COLLECTION-level slots -- ldReference,
# genome and qcOptions -- on top of the elements and their mcols, so row-binding
# per-block pieces back into one collection has to say what happens to each.
# The rules live here once, shared by combineGwasSumStats() and
# combineQtlSumStats(), so the two cannot drift:
#
#   genome    must agree; a mismatch is an error
#   qcOptions carried from the first part after checking the QC options
#             agree; each row's audit is a column, so it binds with the rows
#   ldReference  the panels union, read through the parts' shared GenotypeHandle
#
# The ldReference rule is the one that matters for correctness. A block-parallel
# pipeline narrows each piece's panel to that block, so first-wins would hand
# back a collection spanning every block whose LD reference covers only the
# first -- and cTWAS computes its full-panel LD from exactly that slot.
# =============================================================================

# The concrete class of `x`, for map_chr over a list of collections.
# @noRd
.firstClass <- function(x) class(x)[[1L]]

# Every part must be the same concrete class: the combined object is built as
# `parts[[1L]]`'s class, so a mixed set would silently coerce the rest.
# @noRd
.rtlRequireSameClass <- function(parts, fn) {
    classes <- unique(map_chr(parts, .firstClass))
    if (length(classes) > 1L) {
        msg <- glue(
            "{fn}: inputs must be the same concrete class (got ",
            "{str_flatten(classes, ', ')})."
        )
        abort(msg)
    }
    invisible(NULL)
}

# The collection-level slots BEYOND ldReference that a combine has to merge, for
# whatever concrete class `parts` are. A slot describes the whole collection,
# so first-wins would silently drop the other parts' state; each gets its own
# rule and an unrecognised slot is an error rather than a quiet default, so
# adding one to a class cannot slip through a combine unnoticed.
# @noRd
.rtlExtraSlots <- function(parts, fn) {
    own <- setdiff(.rtlOwnSlots(parts[[1L]]), "ldReference")
    # discard(is.null), not compact(): parts no QC has run on merge to an
    # empty qcOptions, which is a merged value -- compact() would drop it and
    # report the slot as having no merge rule.
    out <- discard(
        list(
            qcOptions = if (is_in("qcOptions", own)) {
                .ssCombineQcOptions(parts, fn)
            }
        ),
        is.null
    )
    unknown <- setdiff(own, names(out))
    if (length(unknown) > 0L) {
        cls <- class(parts[[1L]])[[1L]]
        msg <- glue(
            "{fn}: no merge rule for the collection-level slot(s) ",
            "{str_flatten(unknown, ', ')} on class '{cls}'."
        )
        abort(msg)
    }
    out
}

# Merge `parts` into one collection, reconciling every collection-level slot.
# The single path behind `c()` / `append()` and the four combine*() wrappers.
# Every row keeps its own LD reference; a non-NULL `ldReference` is every
# row's instead.
# @noRd
.combineTupleCollections <- function(parts, ldReference, fn) {
    .rtlRequireSameClass(parts, fn)
    # Forced before the rebuild: GRangesList() merges seqinfo and rejects
    # mismatched builds itself, with a message about sequence-level genomes
    # that says nothing about which inputs disagreed.
    genome <- .rtlCombineGenome(parts, fn)
    .rbindCollections(
        parts,
        ldReference = ldReference,
        slots = .rtlExtraSlots(parts, fn),
        genome = genome,
        fn = fn
    )
}

# The single build every part must agree on, or NULL for a family that does
# not carry one. Read seqinfo directly via GenomeInfoDb::genome().
# @noRd
.rtlCombineGenome <- function(parts, fn) {
    if (!methods::is(parts[[1L]], "SumStatsBase")) {
        return(NULL)
    }
    .ssCombineGenome(parts, fn)
}

# Row-bind summary-statistics parts, merging the three collection-level slots.
# `ldReference` (when non-NULL) overrides the unioned panel.
# @noRd
.rbindSumStats <- function(parts, ldReference, fn) {
    .combineTupleCollections(parts, ldReference, fn)
}

# One genome build per collection (every entry shares the LD sketch), so the
# parts must agree.
# @noRd
.ssCombineGenome <- function(parts, fn) {
    # genome() is per-seqlevel, so flatten across parts before uniquing
    # rather than map_chr (which would reject a multi-seqlevel part).
    builds <- discard(list_c(map(parts, GenomeInfoDb::genome)), is.na)
    genomes <- unique(unname(builds))
    if (length(genomes) == 0L) {
        return(NA_character_)
    }
    if (length(genomes) > 1L) {
        msg <- glue(
            "{fn}: every input must share one genome build (got ",
            "{str_flatten(genomes, ', ')})."
        )
        abort(msg)
    }
    genomes
}

# qcOptions describes the whole collection, so it is carried from the first
# part once every part is shown to have run the same QC. Each row's audit is
# a column of mcols(), so it binds with the rows it describes.
# @noRd
.ssCombineQcOptions <- function(parts, fn) {
    populated <- keep(parts, .ssHasQcOptions)
    if (length(populated) == 0L) {
        return(list())
    }
    if (length(populated) < length(parts)) {
        msg <- glue(
            "{fn}: {length(parts) - length(populated)} of {length(parts)} ",
            "inputs carry no QC record. Run summaryStatsQc() on every input ",
            "before combining, or the combined qcOptions would claim QC ",
            "that only some elements went through."
        )
        abort(msg)
    }
    .ssCheckQcOptions(populated, fn)
    qcOptions(populated[[1L]])
}

# @noRd
.ssHasQcOptions <- function(x) length(qcOptions(x)) > 0L

# The QC options must match: a collection stitched from blocks filtered at
# different MAF / INFO cutoffs has no single honest answer for "what QC ran",
# and every downstream consumer reads that answer off one record.
# @noRd
.ssCheckQcOptions <- function(parts, fn) {
    first <- qcOptions(parts[[1L]])$options
    same <- map_lgl(parts, .ssSameQcOptions, first = first)
    if (all(same)) {
        return(invisible(NULL))
    }
    msg <- glue(
        "{fn}: inputs were QC'd with different summaryStatsQc() options ",
        "(input {which(!same)[[1L]]} differs from the first). Re-run QC with ",
        "one set of options before combining."
    )
    abort(msg)
}

# @noRd
.ssSameQcOptions <- function(x, first) {
    identical(qcOptions(x)$options, first)
}

# The LD references of `parts`, row-bound: each row keeps the reference it
# had, and a reference that is a view of another is stored as that one
# (.ldRefMergeViews()). `override`, when given, is every row's reference
# instead.
# @noRd
.ldRefCombine <- function(parts, override, fn) {
    if (!is.null(override)) {
        n <- sum(map_int(parts, length))
        return(.ldRefDict(.ldRefPerTuple(override, n)))
    }
    perRow <- list_flatten(map(parts, .rowLdRefs))
    .ldRefMergeViews(.ldRefDict(perRow), fn)
}

# `dict` with every entry that is a view of another stored as that one, and
# the rows re-keyed to it. Block-parallel pipelines narrow each piece's panel
# to its own block, so binding the pieces back together would otherwise hold
# one panel per block where there is one panel.
#
# Entries are taken widest first, and each joins the first group whose head
# (its widest entry) it is a view of, or heads a new one. An entry so meets
# the widest reference holding it before any narrower one, whatever order the
# rows came in. Groups are numbered by their first row, so `ld1` stays the
# first row's reference.
# @noRd
.ldRefMergeViews <- function(dict, fn) {
    refs <- dict$refs
    if (length(refs) <= 1L) {
        return(dict)
    }
    grouped <- reduce(
        order(map_int(refs, .ldRefNVariants), decreasing = TRUE),
        .ldRefAssignGroup,
        refs = refs,
        .init = list(heads = integer(0), group = integer(length(refs)))
    )
    seen <- unique(grouped$group)
    group <- match(grouped$group, seen)
    merged <- map(
        seq_along(seen),
        .ldRefGroupRef,
        group = group,
        heads = grouped$heads[seen],
        refs = refs,
        fn = fn
    )
    keys <- str_c("ld", seq_along(merged))
    list(
        refs = set_names(merged, keys),
        keys = keys[group[match(dict$keys, names(refs))]]
    )
}

# `grouped` extended by entry `i`: it joins the first group whose head it is
# a view of, or heads a new group.
# @noRd
#' @importFrom purrr detect_index
.ldRefAssignGroup <- function(grouped, i, refs) {
    hit <- detect_index(grouped$heads, .ldRefHeadHolds, i = i, refs = refs)
    if (hit == 0L) {
        grouped$heads <- c(grouped$heads, i)
        hit <- length(grouped$heads)
    }
    grouped$group[[i]] <- hit
    grouped
}

# @noRd
.ldRefHeadHolds <- function(h, i, refs) {
    .ldRefIsViewOf(refs[[i]], refs[[h]])
}

# TRUE when reference `x` can be stored as `of`: panels reading the same files
# over the same samples, or an LdMatrix or mixture all of whose variants `of`
# holds, with the same LD.
# @noRd
.ldRefIsViewOf <- function(x, of) {
    if (.ldIsPanel(x) && .ldIsPanel(of)) {
        return(.ssSameGenotypeSource(x, of))
    }
    if (.ldIsMatrixRef(x) && .ldIsMatrixRef(of)) {
        return(.ldmCoveredBy(x, of))
    }
    if (.ldIsMixtureRef(x) && .ldIsMixtureRef(of)) {
        return(.ldRefVariantsCoveredBy(x, of) && .gmSameMixture(x, of))
    }
    FALSE
}

# The one reference group `g` is stored as. Panels from one source union. An
# LdMatrix or mixture cannot be widened -- the LD between variants that
# different members hold is in none of them -- so the group's head, which
# holds every other member, stands for it.
# @noRd
.ldRefGroupRef <- function(g, group, heads, refs, fn) {
    members <- unname(refs[group == g])
    if (length(members) == 1L || !.ldIsPanel(members[[1L]])) {
        return(refs[[heads[[g]]]])
    }
    .unionPanelRefs(members, fn)
}

# One panel standing for `panels`, which read the same files over the same
# samples: a handle union keeps the first handle and widens its snpInfo, each
# row's `fileIdx` addressing a position in THAT file, and a panel union binds
# the parts over the first one's samples.
# @noRd
.unionPanelRefs <- function(panels, fn) {
    chromPaths <- .ssUnionChromPaths(panels, fn)
    handles <- map(panels, .ldRefHandle)
    if (any(map_lgl(handles, is.null))) {
        return(.ssUnionPanels(panels, chromPaths))
    }
    handle <- methods::initialize(
        handles[[1L]],
        snpInfo = .ssUnionSnpInfo(handles),
        chromPaths = chromPaths
    )
    .asLdRef(handle)
}

# The union of panels from one source that read through no handle -- a VCF,
# GDS or PLINK1 panel, or a per-chromosome set of them: their rows, each
# once, in genomic order, over the first panel's samples. The parts' assays
# are bound, not read. Parts that are all the same panel union to it.
# @noRd
.ssUnionPanels <- function(panels, chromPaths) {
    first <- panels[[1L]]
    if (all(map_lgl(map(panels, rownames), identical, rownames(first)))) {
        return(first)
    }
    rr <- suppressWarnings(
        exec(c, !!!unname(map(panels, SummarizedExperiment::rowRanges)))
    )
    stacked <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            dosage = exec(
                BiocGenerics::rbind,
                !!!unname(map(panels, .seDosage))
            )
        ),
        rowRanges = rr,
        colData = SummarizedExperiment::colData(first),
        metadata = S4Vectors::metadata(first)
    )
    kept <- which(!duplicated(names(rr)))
    ord <- kept[order(
        as.character(seqnames(rr)[kept]),
        start(rr)[kept]
    )]
    md <- S4Vectors::metadata(stacked)
    if (!is.null(md$genotypeSource)) {
        md$genotypeSource$chromPaths <- chromPaths
    }
    S4Vectors::`metadata<-`(stacked[ord, ], value = md)
}

# TRUE when every variant of `x` is in `by`, coded the same way.
# @noRd
.ldRefVariantsCoveredBy <- function(x, by) {
    m <- matchVariants(
        .ldRefMatchIds(x),
        .ldRefMatchIds(by),
        removeStrandAmbiguous = FALSE
    )
    length(m$idxA) == .ldRefNVariants(x) && all(m$sign > 0)
}

# TRUE when every variant of `x` is in `by`, coded the same way, with the
# same LD.
# @noRd
.ldmCoveredBy <- function(x, by) {
    m <- matchVariants(
        .ldRefMatchIds(x),
        .ldRefMatchIds(by),
        removeStrandAmbiguous = FALSE
    )
    length(m$idxA) == length(x) &&
        all(m$sign > 0) &&
        .ldmAgreesOn(x, by, m$idxA, m$idxB)
}

# Union the per-chromosome shard paths across LD references (panels, or
# bare handles). Two sketches trimmed from one genoMeta panel to different
# chromosomes share path/format/samples and differ ONLY here, so unioning is
# what lets their rows keep routing: a row reaches its file by its own CHR.
# One chromosome mapping to two different files means the parts are not the
# same panel after all.
# @noRd
.ssUnionChromPaths <- function(refs, fn) {
    combined <- .ssConcatChr(map(refs, .ldRefChromPaths))
    byChrom <- split(unname(combined), names(combined))
    walk2(names(byChrom), byChrom, .ssCheckOneFilePerChrom, fn = fn)
    # Every chromosome agrees, so the first mapping stands for all of them.
    combined[!duplicated(names(combined))]
}

# @noRd
.ssConcatChr <- function(pieces) {
    if (length(pieces) == 0L) {
        return(character(0))
    }
    list_c(pieces)
}

# A chromosome resolving to two different files means the inputs are not the
# same panel, so their sketches cannot be unioned.
# @noRd
.ssCheckOneFilePerChrom <- function(ch, paths, fn) {
    if (n_distinct(paths) <= 1L) {
        return(invisible(NULL))
    }
    msg <- glue(
        "{fn}: chromosome '{ch}' maps to two different genotype ",
        "files across the inputs, so their LD references cannot ",
        "be unioned."
    )
    abort(msg)
}

# @noRd
.ssSameGenotypeSource <- function(x, first) {
    # chromPaths is deliberately NOT compared: two sketches trimmed from the
    # same genoMeta panel to different chromosomes differ there and nowhere
    # else, and .ssUnionChromPaths() merges them (erroring on a genuine
    # conflict). Everything below must match for the union to mean anything --
    # rows address the same files, over the same sample axis.
    a <- .ldRefSource(x)
    b <- .ldRefSource(first)
    identical(a$path, b$path) &&
        identical(a$format, b$format) &&
        identical(.ssSampleIds(x), .ssSampleIds(first))
}

# An LD reference's samples: a panel's columns, a bare handle's sample ids.
# @noRd
.ssSampleIds <- function(x) {
    if (methods::is(x, "GenotypeHandle")) sampleIds(x) else colnames(x)
}

# The parts' snpInfo rows, de-duplicated by variant id and returned in genomic
# order so the unioned panel reads like one built over the whole span. The
# first handle's data-frame class is preserved: a handle built with a base
# data.frame keeps one, because some downstream readers index snpInfo
# positionally and a tibble column is a one-column list there.
# @noRd
.ssUnionSnpInfo <- function(handles) {
    first <- snpInfo(handles[[1L]])
    combined <- list_rbind(map(handles, .ssHandleSnpInfo))
    keyed <- all(is_in(c("SNP", "CHR", "BP"), names(combined)))
    si <- if (nrow(combined) > 0L && keyed) {
        deduped <- combined[!duplicated(combined$SNP), , drop = FALSE]
        deduped[
            order(deduped$CHR, as.integer(deduped$BP)),
            ,
            drop = FALSE
        ]
    } else {
        combined
    }
    if (inherits(first, "tbl_df")) si else as.data.frame(si)
}

# @noRd
.ssHandleSnpInfo <- function(h) as_tibble(snpInfo(h))
