# =============================================================================
# Tuple row matchers
# -----------------------------------------------------------------------------
# Internal helpers shared by the FineMappingResult / TwasWeights /
# SumStats DFrame-subclass collections to resolve a tuple-keyed selection
# to a single row index. Pure R helpers -- no S4 dispatch, no exports.
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
    ok <- rep(TRUE, nrow(x))
    for (k in names(keys)) {
        ok <- ok & as.character(.tupleColumn(x, k)) == keys[[k]]
    }
    which(ok)
}

# Read an identity column by name, whichever collection shape `x` has.
# `x[[k]]` is a column on the DFrame-backed collections but an ELEMENT on a
# RangedTupleList, where the identity columns live in mcols. The two shapes
# coexist until every collection has migrated.
# @noRd
.tupleColumn <- function(x, k) {
    if (methods::is(x, "RangedTupleList")) {
        return(mcols(x)[[k]])
    }
    x[[k]]
}

# One collection ELEMENT, whichever shape `x` has. On a RangedTupleList the
# elements are the container itself; the DFrame-backed collections still keep
# them in an `entry` column. Both shapes coexist until every collection has
# migrated, so polymorphic call sites go through here.
# @noRd
.collectionEntry <- function(x, i) {
    if (methods::is(x, "TwasWeights")) {
        return(.twrRowParts(x, i))
    }
    if (methods::is(x, "FineMappingResultBase")) {
        return(.fmrRowParts(x, i))
    }
    if (methods::is(x, "RangedTupleList")) {
        return(x[[i]])
    }
    x$entry[[i]]
}

# ALL elements as a plain list, whichever shape `x` has.
# @noRd
.collectionEntries <- function(x) {
    # TwasWeights entries are derived views, so they have to be rebuilt one by
    # one rather than read off as elements (which are bare GRanges).
    if (
        methods::is(x, "TwasWeights") ||
            methods::is(x, "FineMappingResultBase")
    ) {
        return(map(seq_len(nrow(x)), .collectionEntry, x = x))
    }
    if (methods::is(x, "RangedTupleList")) {
        return(as.list(x))
    }
    as.list(x$entry)
}

# The identity column NAMES, whichever shape `x` has. `names(x)` is element
# names on a RangedTupleList, not columns.
# @noRd
.tupleColumnNames <- function(x) {
    if (methods::is(x, "RangedTupleList")) {
        return(colnames(mcols(x)))
    }
    names(x)
}

# Internal: resolve a tuple-keyed selection (study, context, trait,
# method) to a single row index. Used by the QtlFineMappingResult and
# TwasWeights accessors. Returns an error when no row matches; returns
# the single row index when the collection has exactly one row and any
# selector argument was omitted.
.tupleSelectRow <- function(
    x,
    study,
    context,
    trait,
    method,
    cls = "QtlFineMappingResult"
) {
    if (nrow(x) == 0L) {
        msg <- glue("{cls} has no rows.")
        abort(msg)
    }
    anyUnset <- missing(study) ||
        is.null(study) ||
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
            "{cls} has {nrow(x)} entries. Pass `study`, `context`, ",
            "`trait`, and `method` to select one."
        )
        abort(msg)
    }
    if (
        length(study) != 1L ||
            length(context) != 1L ||
            length(trait) != 1L ||
            length(method) != 1L
    ) {
        abort(
            "`study`, `context`, `trait`, and `method` must each be length 1."
        )
    }
    .tupleMatchQtl(x, study, context, trait, method)
}

# Resolve a (study, context, trait, method) tuple to a single row index.
# @noRd
.tupleMatchQtl <- function(x, study, context, trait, method) {
    idx <- .matchTupleRows(
        x,
        list(study = study, context = context, trait = trait, method = method)
    )
    if (length(idx) == 0L) {
        msg <- glue(
            "No entry for (study='{study}', context='{context}', ",
            "trait='{trait}', method='{method}')."
        )
        abort(msg)
    }
    idx[[1L]]
}

# Internal: resolve a (study, method, blockId) tuple to a single row
# index of a GwasFineMappingResult collection. `region` may be NULL when
# the (study, method) pair maps to a single row; otherwise it disambiguates
# among per-block rows of a genome-wide collection.
.tupleSelectRowGwasFmr <- function(x, study, method, region = NULL) {
    if (nrow(x) == 0L) {
        abort("GwasFineMappingResult has no rows.")
    }
    anyUnset <- missing(study) ||
        is.null(study) ||
        missing(method) ||
        is.null(method)
    if (anyUnset) {
        if (nrow(x) == 1L) {
            return(1L)
        }
        msg <- glue(
            "GwasFineMappingResult has {nrow(x)} entries. Pass `study` ",
            "and `method` to select one."
        )
        abort(msg)
    }
    if (length(study) != 1L || length(method) != 1L) {
        abort("`study` and `method` must each be length 1.")
    }
    if (!is.null(region) && length(region) != 1L) {
        abort("`region` must be length 1 when supplied.")
    }
    .tupleMatchGwas(x, study, method, region)
}

# Resolve a (study, method[, blockId]) tuple to a single row index.
# @noRd
.tupleMatchGwas <- function(x, study, method, region) {
    keys <- list(study = study, method = method)
    if (!is.null(region)) {
        keys$blockId <- region
    }
    idx <- .matchTupleRows(x, keys)
    if (length(idx) == 0L) {
        regionPart <- if (is.null(region)) {
            ""
        } else {
            glue(", region='{region}'")
        }
        msg <- glue(
            "No entry for (study='{study}', method='{method}'{regionPart})."
        )
        abort(msg)
    }
    if (length(idx) > 1L) {
        .tupleGwasAmbiguous(x, study, method, idx)
    }
    idx[[1L]]
}

# Multiple (study, method) rows matched: report the disambiguating regions.
# @noRd
.tupleGwasAmbiguous <- function(x, study, method, idx) {
    regions <- str_flatten(
        shQuote(as.character(x$blockId[idx])),
        ", "
    )
    msg <- glue(
        "GwasFineMappingResult has {length(idx)} rows matching ",
        "(study='{study}', method='{method}'); pass `region` to ",
        "disambiguate (available: {regions})."
    )
    abort(msg)
}

# Internal: row indices of a FineMappingResultBase collection matching the
# given selectors. NULL selectors, and selectors naming a column the
# collection lacks, are ignored -- QTL rows carry study/context/trait/method,
# GWAS rows carry study/method/blockId -- so the same call works on either.
# `region` matches the `blockId` column. Returns all rows when nothing
# constrains. Unlike the single-row selectors above this never errors on
# ambiguity; it is the aggregate counterpart used by getTopLoci.
.fmrRowsMatching <- function(
    x,
    study = NULL,
    context = NULL,
    trait = NULL,
    method = NULL,
    region = NULL
) {
    keys <- list(
        study = study,
        context = context,
        trait = trait,
        method = method,
        blockId = region
    )
    keys <- keys[!map_lgl(keys, is.null)]
    keys <- keys[is_in(names(keys), .tupleColumnNames(x))]
    if (length(keys) == 0L) {
        return(seq_len(nrow(x)))
    }
    ok <- rep(TRUE, nrow(x))
    for (k in names(keys)) {
        ok <- ok &
            is_in(as.character(.tupleColumn(x, k)), as.character(keys[[k]]))
    }
    which(ok)
}

# Internal: per-row identity metadata for a FineMappingResultBase collection,
# as a data.frame with a stable column set (study, context, trait, blockId,
# method). Columns the collection lacks are NA-filled so QTL and GWAS results
# yield the same metadata shape. Safe for zero-row collections.
.fmrRowMetadata <- function(x) {
    cols <- c("study", "context", "trait", "blockId", "method")
    n <- nrow(x)
    vals <- map(cols, .fmrMetadataCol, x = x, n = n)
    names(vals) <- cols
    as_tibble(vals)
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

# Internal: read the per-row `region` GRanges column of a collection, or an
# empty GRanges when the collection carries none. Shared by getRegion on
# TwasWeights and FineMappingResultBase (region is a uniform column across the
# family; no derivation from topLoci / blockId).
# The per-element span, DERIVED from the ranges rather than read from a stored
# `region` column (spec 4.4). range() is always in sync by construction, where
# a stored window had no correct update rule under subsetRegion() and would
# quietly go stale. `blockId` carries the nominal block identity instead.
# @noRd
.getRegionColumn <- function(x) {
    if (nrow(x) == 0L) {
        return(GenomicRanges::GRanges())
    }
    if (methods::is(x, "RangedTupleList")) {
        return(unlist(range(x), use.names = FALSE))
    }
    # No fallback to a stored column: nothing carries one any more. A DFrame
    # collection (CtwasResult) has no region column either, so an empty
    # GRanges is the honest answer rather than a fabricated span.
    GenomicRanges::GRanges()
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
            "`blockId` must have the same length as `study` ",
            "(got {length(blockId)} vs {n})."
        )
        abort(msg)
    }
    cols[["blockId"]] <- as.character(blockId)
    cols
}

# Internal: append a validated `region` GRanges to a constructor's column list
# (no-op when region is NULL). Shared by the TwasWeights / FineMappingResult
# constructors so the provenance column is added identically everywhere.

# Internal: validate the optional `region` GRanges column (one range per row),
# shared by the TwasWeights / FineMappingResult validity methods. Returns a
# character vector of errors (empty when valid or the column is absent).

# Internal: trait-position provenance. One GRanges per trait carrying the
# molecular feature's OWN genomic coordinates (gene/peak range; TSS = start()),
# distinct from the fine-mapping `region`. The true trait position cannot be
# inferred from summary statistics, so it is threaded from QtlDataset rowRanges
# / an explicit QtlSumStats trait-position and carried onto QtlFineMappingResult
# and TwasWeights. Mirrors the `region` helpers above; absent for GWAS.
# traitPos is optional provenance: always known for a QtlDataset, but only known
# for a QtlSumStats when the caller supplied it (it cannot be inferred from
# summary statistics). When the column is absent we return a scalar NA rather
# than an empty GRanges, so getTraitPosition() reports "no trait position"
# honestly instead of a zero-length range.
.getTraitPosColumn <- function(x) {
    if (is_in("traitPos", .tupleColumnNames(x))) {
        .tupleColumn(x, "traitPos")
    } else {
        NA
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
    cols[["traitPos"]] <- traitPos
    cols
}

.validateTraitPosColumn <- function(object) {
    if (!is_in("traitPos", .tupleColumnNames(object))) {
        return(character(0))
    }
    traitPos <- .tupleColumn(object, "traitPos")
    if (!methods::is(traitPos, "GRanges")) {
        return("'traitPos' column must be a GRanges")
    }
    if (length(traitPos) != nrow(object)) {
        return("'traitPos' column must have one range per row")
    }
    character(0)
}

# Internal: a length-n "unknown" column matching the type of `exemplar`, used
# to pad a collection that lacks an optional column before a union row-bind.
# Character -> NA; GRanges -> a 0-width chrUn sentinel range; else NA.
.naLikeColumn <- function(exemplar, n) {
    if (methods::is(exemplar, "GRanges")) {
        return(GenomicRanges::GRanges(
            rep("chrUn", n),
            IRanges::IRanges(start = 1L, width = 0L)
        ))
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
    # unreachable: cn is always drawn from allCols, so some part has it
    NULL # nocov
}

# Internal: row-bind two or more per-tuple collection objects (TwasWeights /
# FineMappingResult subclasses), carrying forward EVERY column. Columns are
# unioned, and a collection lacking an optional column (e.g. jointContexts or
# region) is padded with `.naLikeColumn` -- so adding a column to a class flows
# through combines automatically rather than being hand-listed at each site.
# Preserves the concrete class and sets the `ldSketch` slot explicitly (rbind
# does not reliably carry it). Returns NULL when given no inputs.
.rbindCollections <- function(parts, ldSketch = NULL) {
    parts <- compact(parts)
    if (length(parts) == 0L) {
        return(NULL)
    }
    if (methods::is(parts[[1L]], "RangedTupleList")) {
        return(.rbindRangedCollections(parts, ldSketch))
    }
    cls <- class(parts[[1L]])[[1L]]
    allCols <- reduce(map(parts, names), union)
    combined <- map(allCols, .rbindColumn, parts = parts)
    names(combined) <- allCols
    dfArgs <- c(combined, list(check.names = FALSE))
    df <- exec(S4Vectors::DataFrame, !!!dfArgs)
    out <- new(cls, df, ldSketch = ldSketch)
    validObject(out)
    out
}

# Internal: aggregate a per-entry accessor across every row of a
# FineMappingResultBase collection that matches the given selectors, prefixing
# each entry's rows with the row identity (study/context/trait/blockId/method)
# so rows stay attributable to their source entry. Shared by getTopLoci / getCs
# / getMarginalEffects on FineMappingResultBase.
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
    study = NULL,
    context = NULL,
    trait = NULL,
    method = NULL,
    region = NULL,
    perEntry,
    onSingle = perEntry,
    ...
) {
    single <- .fmrTrySingle(
        x,
        study,
        context,
        trait,
        method,
        region,
        onSingle,
        ...
    )
    if (single$done) {
        return(single$result)
    }
    idx <- .fmrRowsMatching(
        x,
        study = study,
        context = context,
        trait = trait,
        method = method,
        region = region
    )
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
        ...
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

# Attempt the single-entry fast path when any selector is supplied. Returns
# list(done, result) on a successful single-select, else list(done = FALSE, sel)
# carrying the selection error (or NULL when no selector was given).
# @noRd
.fmrTrySingle <- function(
    x,
    study,
    context,
    trait,
    method,
    region,
    onSingle,
    ...
) {
    anySelector <- !is.null(study) ||
        !is.null(context) ||
        !is.null(trait) ||
        !is.null(method) ||
        !is.null(region)
    if (!anySelector) {
        return(list(done = FALSE, sel = NULL))
    }
    sel <- tryCatch(
        .fmrSelectEntry(
            x,
            study = study,
            context = context,
            trait = trait,
            method = method,
            region = region
        ),
        error = function(e) e
    )
    if (!inherits(sel, "error")) {
        return(list(done = TRUE, result = onSingle(sel, ...)))
    }
    list(done = FALSE, sel = sel)
}

# Project one entry to its per-entry view with the row's identity metadata
# prepended; NULL when the entry yields no rows.
# @noRd
.fmrEntryPart <- function(i, x, perEntry, meta, ...) {
    # `perEntry` is an internal per-row function taking the stored payload, not
    # an S4 generic dispatching on a rebuilt entry. That indirection is what
    # the FineMappingRow class existed for.
    v <- perEntry(.fmrRowParts(x, i), ...)
    if (is.null(v) || nrow(v) == 0L) {
        return(NULL)
    }
    v <- select(v, all_of(setdiff(names(v), names(meta))))
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

# Concatenate column `cn` across all parts (NA-filling parts that lack it).
# @noRd
# Row-binding a RANGED collection: the elements append and the mcols rbind.
# The DFrame path cannot be reused because there is no column holding the
# payload any more -- the payload is the container.
# @noRd
.rbindRangedCollections <- function(parts, ldSketch) {
    cls <- class(parts[[1L]])[[1L]]
    allCols <- reduce(map(parts, .tupleColumnNames), union)
    combined <- set_names(map(allCols, .rbindColumn, parts = parts), allCols)
    md <- exec(
        S4Vectors::DataFrame,
        !!!c(combined, list(check.names = FALSE))
    )
    elements <- list_flatten(map(parts, as.list))
    grl <- GenomicRanges::GRangesList(elements)
    mcols(grl) <- md
    out <- new(cls, grl, ldSketch = ldSketch)
    validObject(out)
    out
}

.rbindColumn <- function(cn, parts) {
    pieces <- map(parts, .rbindColPiece, cn = cn, parts = parts)
    # GRanges concat may warn on seqinfo.
    suppressWarnings(exec(c, !!!pieces))
}

# =============================================================================
# Per-row primitives
# -----------------------------------------------------------------------------
# The per-row views a collection projects. These replace the FineMappingRow
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
        cvResult = md$cvResult[[first]]
    )
}

# The same, for a set of rows collapsed into one logical entry (a multi-block
# row group). The fit payload comes from the first row, matching the previous
# .fmeEntryFromRows behaviour.
# @noRd

# @noRd
.fmrPartsVariantIds <- function(parts) {
    getVariantIds(.asFmRowPayload(parts))
}

# @noRd
.fmrPartsSusieFit <- function(parts) {
    getSusieFit(.asFmRowPayload(parts))
}

# @noRd
.fmrPartsCvResult <- function(parts) {
    getCvResult(.asFmRowPayload(parts))
}

# @noRd
.fmrPartsTopLoci <- function(parts) {
    .fmeTopLociFromElement(rowVariants(.asFmRowPayload(parts)))
}

# ---- per-row views ----------------------------------------------------------
# Each of these is the body of what used to be a FineMappingRow method,
# reading the stored form through `parts` instead of through S4 slots. They are
# what .fmrAggregateView() maps over.

# @noRd
.fmrRowTopLoci <- function(
    parts,
    type = c("data.frame", "GRanges"),
    signalCutoff = 0.025,
    minPurity = NULL,
    raw = FALSE,
    ...
) {
    tl <- .fmrPartsTopLoci(parts)
    # raw = TRUE returns the stored canonical table verbatim: every variant,
    # every column, no posterior-view projection.
    if (isTRUE(raw)) {
        return(tl)
    }
    type <- arg_match(type)
    out <- .fmeFilterTopLoci(tl, signalCutoff, minPurity)
    if (type == "data.frame") {
        return(out)
    }
    .fmeTopLociGRanges(out)
}

# @noRd
.fmrRowPip <- function(parts, ...) {
    tl <- .fmrPartsTopLoci(parts)
    if (nrow(tl) == 0L || !is_in("pip", names(tl))) {
        return(numeric(0))
    }
    set_names(tl$pip, tl$variant_id)
}

# @noRd
.fmrRowMarginalEffects <- function(parts, maxPval = NULL, ...) {
    tl <- .fmrPartsTopLoci(parts)
    if (nrow(tl) == 0L) {
        return(.projectMarginalView(tl))
    }
    out <- .projectMarginalView(tl)
    if (!is.null(maxPval) && nrow(out) > 0L) {
        keep <- !is.na(out$p) & out$p <= maxPval
        out <- filter(out, keep)
    }
    out
}

# @noRd
.fmrRowCs <- function(parts, coverage = 0.95, minPurity = NULL, ...) {
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
    keep <- !is.na(tl[[csCol[1L]]]) &
        str_length(tl[[csCol[1L]]]) > 0L &
        !str_detect(tl[[csCol[1L]]], "_0$")
    # Independent purity filter (min.abs.corr), orthogonal to `coverage`: drop
    # CS members whose credible set at THIS coverage is below `minPurity`.
    if (!is.null(minPurity)) {
        purCol <- str_c(csCol[1L], "_purity")
        if (is_in(purCol, names(tl))) {
            pur <- as.numeric(tl[[purCol]])
            keep <- keep & !is.na(pur) & pur >= minPurity
        } else {
            msg <- glue(
                "getCs: no purity column '{purCol}' for coverage ",
                "{coverage}; minPurity filter skipped."
            )
            warn(msg)
        }
    }
    .projectPosteriorView(tl[keep, , drop = FALSE])
}

# @noRd
.fmrRowLbf <- function(parts, ...) {
    lbf <- .asLbfMatrix(getSusieFit(parts))
    vids <- .fmrPartsVariantIds(parts)
    if (is.null(lbf) || ncol(lbf) != length(vids)) {
        return(tibble(variant_id = character(0)))
    }
    w <- as_tibble(t(as.matrix(lbf)), .name_repair = "minimal")
    names(w) <- str_c("lbf_L", seq_len(ncol(w)))
    bind_cols(tibble(variant_id = as.character(vids)), w)
}

# @noRd
.fmrRowCredibleSetSummary <- function(parts, coverage = 0.95, ...) {
    .csSummaryFit(.fmrPartsTopLoci(parts), getSusieFit(parts), coverage)
}

# @noRd
.fmrRowFsusieCredibleBand <- function(parts, ...) {
    .fsusieCredibleBandFit(getSusieFit(parts))
}

# @noRd
.fmrRowFsusieAffectedRegions <- function(parts, ...) {
    .fsusieAffectedRegionsFit(
        getSusieFit(parts),
        topLoci = .fmrPartsTopLoci(parts)
    )
}

# @noRd
.fmrRowResolveWeights <- function(parts, ...) {
    empty <- list(variantIds = character(0), weights = numeric(0))
    # The topLoci posterior view projects the effect to `beta`; use it as the
    # per-variant weight, aligned with variant_id.
    tl <- .fmrRowTopLoci(parts)
    if (
        is.null(tl) ||
            nrow(tl) == 0L ||
            !all(is_in(c("variant_id", "beta"), names(tl)))
    ) {
        return(empty)
    }
    vids <- as.character(tl$variant_id)
    w <- as.numeric(tl$beta)
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
        fits = md$fits[[first]],
        cvResult = md$cvResult[[first]],
        standardized = isTRUE(md$standardized[[first]]),
        dataType = md$dataType[[first]]
    )
}

# @noRd
.twrPartsVariantIds <- function(parts) {
    getVariantIds(parts)
}

# The weight vector aligned to the row's variant ids -- what resolveWeights
# produced from a TwasWeightsRow.
# @noRd
.twrRowResolveWeights <- function(parts, ...) {
    empty <- list(variantIds = character(0), weights = numeric(0))
    vids <- .twrPartsVariantIds(parts)
    w <- getWeights(parts)
    if (length(vids) == 0L || is.null(w)) {
        return(empty)
    }
    w <- as.numeric(w)
    if (length(w) != length(vids)) {
        return(empty)
    }
    ok <- !is.na(vids) & !is.na(w)
    if (!any(ok)) {
        return(empty)
    }
    list(variantIds = vids[ok], weights = w[ok])
}

# Class-aware row payload: the TwasWeights and FineMappingResult shapes differ
# (weights vs a susie fit), so callers that accept either weight source route
# through here. Replaces the old .collectionEntry bridge.
# @noRd
.rowParts <- function(x, i) {
    if (methods::is(x, "TwasWeights")) {
        return(.twrRowParts(x, i))
    }
    .fmrRowParts(x, i)
}

# The per-variant weight vector of one row, whichever weight source it came
# from.
# @noRd
.rowResolveWeights <- function(parts, ...) {
    if (methods::is(parts, "TwasWeightsRow")) {
        return(.twrRowResolveWeights(parts, ...))
    }
    .fmrRowResolveWeights(parts, ...)
}

# The cross-validated / per-method fits of one row, or NULL for a
# fine-mapping row.
#
# A fine-mapping row deliberately reports NULL: cTWAS would otherwise
# renormalize alpha into an UNstandardized weight, which is inconsistent with
# the standardized posterior effect the topLoci view already carries. Preserved
# verbatim from the FineMappingRow getFits method.
# @noRd
.rowFits <- function(parts) {
    if (methods::is(parts, "TwasWeightsRow")) getFits(parts) else NULL
}

# Whether the row's weights are already on the standardized scale.
#
# TRUE for a fine-mapping row: its posterior effect is colSums(alpha * mu),
# which does not divide by the column scale factors, so cTWAS's per-variant
# variance scaling would double-standardize it.
# @noRd
.rowStandardized <- function(parts) {
    if (methods::is(parts, "TwasWeightsRow")) {
        return(getStandardized(parts))
    }
    TRUE
}

# ---- row-payload builders ---------------------------------------------------
# What callers use in place of the retired FineMappingRow / TwasWeightsRow
# constructors. Same arguments, but they return the STORED form directly -- a
# plain list holding one GRanges plus the fit payload -- so nothing has to be
# converted on the way into a collection.

# Accept either a row payload (the stored form) or a legacy entry object, so
# the constructors can be flipped independently of their ~300 call sites.
# Delete the entry branch once nothing constructs entries.
# @noRd
.asFmRowPayload <- function(e) {
    # A single-row collection is the other shape callers hand over -- it is
    # what getFineMappingResult() returns, so anything accepting "one row's
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

# @noRd

# Shape-tolerant readers for values that may still be a legacy entry object
# while the flip is in progress. Delete the entry branch (and these) once
# nothing produces entries.
# @noRd
.rowCvResult <- function(e) {
    getCvResult(e)
}

# @noRd
.rowWeights <- function(e) {
    getWeights(e)
}

# @noRd
.rowVariantIds <- function(e) {
    getVariantIds(e)
}

# The non-variant half of each TWAS row payload, as outer mcols columns.
# Objects of arbitrary shape (fits, cvResult, dataType) ride in SimpleLists.
# @noRd
.twRowPayloadCols <- function(entry) {
    list(
        fits = S4Vectors::SimpleList(map(entry, getFits)),
        cvResult = S4Vectors::SimpleList(map(entry, getCvResult)),
        standardized = map_lgl(entry, .twPayloadStandardized),
        dataType = S4Vectors::SimpleList(map(entry, getDataType))
    )
}
