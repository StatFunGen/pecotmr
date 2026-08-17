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
        ok <- ok & as.character(x[[k]]) == keys[[k]]
    }
    which(ok)
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

# Internal: resolve a (study, method, region_id) tuple to a single row
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

# Resolve a (study, method[, region_id]) tuple to a single row index.
# @noRd
.tupleMatchGwas <- function(x, study, method, region) {
    keys <- list(study = study, method = method)
    if (!is.null(region)) {
        keys$region_id <- region
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
        shQuote(as.character(x$region_id[idx])),
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
# GWAS rows carry study/method/region_id -- so the same call works on either.
# `region` matches the `region_id` column. Returns all rows when nothing
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
        region_id = region
    )
    keys <- keys[!map_lgl(keys, is.null)]
    keys <- keys[is_in(names(keys), names(x))]
    if (length(keys) == 0L) {
        return(seq_len(nrow(x)))
    }
    ok <- rep(TRUE, nrow(x))
    for (k in names(keys)) {
        ok <- ok & is_in(as.character(x[[k]]), as.character(keys[[k]]))
    }
    which(ok)
}

# Internal: per-row identity metadata for a FineMappingResultBase collection,
# as a data.frame with a stable column set (study, context, trait, region_id,
# method). Columns the collection lacks are NA-filled so QTL and GWAS results
# yield the same metadata shape. Safe for zero-row collections.
.fmrRowMetadata <- function(x) {
    cols <- c("study", "context", "trait", "region_id", "method")
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
# family; no derivation from topLoci / region_id).
.getRegionColumn <- function(x) {
    if (is_in("region", names(x))) x[["region"]] else GenomicRanges::GRanges()
}

# Internal: append a validated `region` GRanges to a constructor's column list
# (no-op when region is NULL). Shared by the TwasWeights / FineMappingResult
# constructors so the provenance column is added identically everywhere.
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
    cols[["region"]] <- region
    cols
}

# Internal: validate the optional `region` GRanges column (one range per row),
# shared by the TwasWeights / FineMappingResult validity methods. Returns a
# character vector of errors (empty when valid or the column is absent).
.validateRegionColumn <- function(object) {
    if (!is_in("region", names(object))) {
        return(character(0))
    }
    if (!methods::is(object[["region"]], "GRanges")) {
        return("'region' column must be a GRanges")
    }
    if (length(object[["region"]]) != nrow(object)) {
        return("'region' column must have one range per row")
    }
    character(0)
}

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
    if (is_in("traitPos", names(x))) x[["traitPos"]] else NA
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
    if (!is_in("traitPos", names(object))) {
        return(character(0))
    }
    if (!methods::is(object[["traitPos"]], "GRanges")) {
        return("'traitPos' column must be a GRanges")
    }
    if (length(object[["traitPos"]]) != nrow(object)) {
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
        if (is_in(cn, names(p))) return(p[[cn]])
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
# each entry's rows with the row identity (study/context/trait/region_id/method)
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
#     metadata prefix is attached (so a stamped `method` never duplicates).
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
    v <- perEntry(x$entry[[i]], ...)
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
    if (is_in(cc, names(x))) as.character(x[[cc]]) else rep(NA_character_, n)
}

# Column `cn` of part `p`, or a NA-like column of an exemplar type when absent.
# @noRd
.rbindColPiece <- function(p, cn, parts) {
    if (is_in(cn, names(p))) {
        p[[cn]]
    } else {
        .naLikeColumn(.exemplarColumn(cn, parts), nrow(p))
    }
}

# Concatenate column `cn` across all parts (NA-filling parts that lack it).
# @noRd
.rbindColumn <- function(cn, parts) {
    pieces <- map(parts, .rbindColPiece, cn = cn, parts = parts)
    # GRanges concat may warn on seqinfo.
    suppressWarnings(exec(c, !!!pieces))
}
