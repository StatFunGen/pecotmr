# =============================================================================
# Tuple row matchers
# -----------------------------------------------------------------------------
# Internal helpers shared by the FineMappingResult / TwasWeights / SumStats
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
    ok <- rep(TRUE, nrow(x))
    for (k in names(keys)) {
        ok <- ok & as.character(.tupleColumn(x, k)) == keys[[k]]
    }
    which(ok)
}

# Read an identity column by name. Every collection keeps its identity
# columns in mcols; `x[[k]]` addresses an ELEMENT, not a column, so mcols is
# the only correct read.
# @noRd
.tupleColumn <- function(x, k) {
    mcols(x)[[k]]
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
    x[[i]]
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
    # No stored region column: nothing carries one any more, so the span is
    # always derived from the elements' own ranges.
    unlist(range(x), use.names = FALSE)
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
# does not reliably carry it). Returns NULL when given no inputs. Every
# collection is a RangedTupleList now, so this only drops the NULL parts and
# hands off; the old DFrame path went with the rebase onto GRangesList.
#
# `slots` carries the collection-level state a concrete subclass adds beyond
# `ldSketch` (SumStatsBase's `genome` / `qcInfo`). It is passed in rather than
# read off `parts[[1L]]` because those slots describe the WHOLE collection, so
# merging them is the caller's decision -- first-wins would silently drop the
# other parts' QC audit.
.rbindCollections <- function(
    parts,
    ldSketch = NULL,
    slots = list(),
    genome = NULL
) {
    parts <- compact(parts)
    if (length(parts) == 0L) {
        return(NULL)
    }
    .rbindRangedCollections(parts, ldSketch, slots, genome)
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
    ldSketch,
    slots = list(),
    genome = NULL
) {
    cls <- class(parts[[1L]])[[1L]]
    allCols <- reduce(map(parts, .tupleColumnNames), union)
    combined <- set_names(map(allCols, .rbindColumn, parts = parts), allCols)
    md <- exec(
        S4Vectors::DataFrame,
        !!!c(combined, list(check.names = FALSE))
    )
    elements <- list_flatten(map(parts, as.list))
    grl <- GenomicRanges::GRangesList(elements)
    # Rebuilding from as.list() merges each part's seqinfo, so a part whose
    # elements never had a build set would leave the result carrying both the
    # real build and NA. The agreed build is written back explicitly to keep
    # seqinfo single-valued, which validity requires.
    if (!is.null(genome)) {
        GenomeInfoDb::genome(grl) <- genome
    }
    mcols(grl) <- md
    out <- exec(
        methods::new,
        cls,
        grl,
        !!!c(list(ldSketch = ldSketch), slots)
    )
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

# =============================================================================
# Combining summary-statistics collections
# -----------------------------------------------------------------------------
# A SumStatsBase subclass carries three COLLECTION-level slots -- ldSketch,
# genome and qcInfo -- on top of the elements and their mcols, so row-binding
# per-block pieces back into one collection has to say what happens to each.
# The rules live here once, shared by combineGwasSumStats() and
# combineQtlSumStats(), so the two cannot drift:
#
#   genome    must agree; a mismatch is an error
#   qcInfo    `entryAudit` concatenates in element order (it is indexed BY
#             element); the rest is carried from the first part after checking
#             the QC options agree
#   ldSketch  the panels union, read through the parts' shared GenotypeHandle
#
# The ldSketch rule is the one that matters for correctness. A block-parallel
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

# The collection-level slots BEYOND ldSketch that a combine has to merge, for
# whatever concrete class `parts` are. A slot describes the whole collection,
# so first-wins would silently drop the other parts' state; each gets its own
# rule and an unrecognised slot is an error rather than a quiet default, so
# adding one to a class cannot slip through a combine unnoticed.
# @noRd
.rtlExtraSlots <- function(parts, fn) {
    own <- setdiff(.rtlOwnSlots(parts[[1L]]), "ldSketch")
    out <- list()
    if (is_in("qcInfo", own)) {
        out$qcInfo <- .ssCombineQcInfo(parts, fn)
    }
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
# A non-NULL `ldSketch` overrides the unioned panel.
# @noRd
.combineTupleCollections <- function(parts, ldSketch, fn) {
    .rtlRequireSameClass(parts, fn)
    # Forced before the rebuild: GRangesList() merges seqinfo and rejects
    # mismatched builds itself, with a message about sequence-level genomes
    # that says nothing about which inputs disagreed.
    genome <- .rtlCombineGenome(parts, fn)
    .rbindCollections(
        parts,
        ldSketch = ldSketch %||% .combineLdSketch(parts, fn),
        slots = .rtlExtraSlots(parts, fn),
        genome = genome
    )
}

# The single build every part must agree on, or NULL for a family that does
# not carry one. Read through getGenome() so it works off seqinfo.
# @noRd
.rtlCombineGenome <- function(parts, fn) {
    if (!methods::is(parts[[1L]], "SumStatsBase")) {
        return(NULL)
    }
    .ssCombineGenome(parts, fn)
}

# Row-bind summary-statistics parts, merging the three collection-level slots.
# `ldSketch` (when non-NULL) overrides the unioned panel.
# @noRd
.rbindSumStats <- function(parts, ldSketch, fn) {
    .combineTupleCollections(parts, ldSketch, fn)
}

# One genome build per collection (every entry shares the LD sketch), so the
# parts must agree.
# @noRd
.ssCombineGenome <- function(parts, fn) {
    genomes <- unique(map_chr(parts, getGenome))
    if (length(genomes) > 1L) {
        msg <- glue(
            "{fn}: every input must share one genome build (got ",
            "{str_flatten(genomes, ', ')})."
        )
        abort(msg)
    }
    genomes
}

# qcInfo is collection-level EXCEPT `entryAudit`, which getQcDiagnostics()
# addresses BY element index -- so the audit concatenates in element order
# while the rest is carried from the first part. Concatenation pads a part
# whose audit is absent or short, because an audit that slipped out of step
# with the elements would report another block's diagnostics.
# @noRd
.ssCombineQcInfo <- function(parts, fn) {
    populated <- keep(parts, .ssHasQcInfo)
    if (length(populated) == 0L) {
        return(list())
    }
    if (length(populated) < length(parts)) {
        msg <- glue(
            "{fn}: {length(parts) - length(populated)} of {length(parts)} ",
            "inputs carry no QC record. Run summaryStatsQc() on every input ",
            "before combining, or the combined qcInfo would claim QC that ",
            "only some elements went through."
        )
        abort(msg)
    }
    .ssCheckQcOptions(populated, fn)
    out <- getQcInfo(populated[[1L]])
    out$entryAudit <- list_flatten(map(parts, .ssEntryAudit))
    out
}

# @noRd
.ssHasQcInfo <- function(x) length(getQcInfo(x)) > 0L

# One part's per-element audit, padded to its element count so positions stay
# aligned with the elements they describe.
# @noRd
.ssEntryAudit <- function(x) {
    audit <- getQcInfo(x)$entryAudit %||% list()
    n <- nrow(x)
    if (length(audit) >= n) {
        return(as.list(audit)[seq_len(n)])
    }
    c(as.list(audit), vector("list", n - length(audit)))
}

# The QC options must match: a collection stitched from blocks filtered at
# different MAF / INFO cutoffs has no single honest answer for "what QC ran",
# and every downstream consumer reads that answer off one record.
# @noRd
.ssCheckQcOptions <- function(parts, fn) {
    first <- getQcInfo(parts[[1L]])$options
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
    identical(getQcInfo(x)$options, first)
}

# The union of the parts' LD panels. All-NULL stays NULL (an individual-level
# collection has no panel); a mix is an error, because silently keeping the
# panels that exist would leave elements harmonized against nothing.
# @noRd
.combineLdSketch <- function(parts, fn) {
    sketches <- map(parts, getLdSketch)
    present <- !map_lgl(sketches, is.null)
    if (!any(present)) {
        return(NULL)
    }
    if (!all(present)) {
        msg <- glue(
            "{fn}: {sum(!present)} of {length(parts)} inputs carry no ",
            "ldSketch. Combine collections that all have an LD reference, or ",
            "all have none."
        )
        abort(msg)
    }
    handles <- map(sketches, .ldSketchHandle)
    .ssCheckSameSource(handles, fn)
    handle <- handles[[1L]]
    handle@snpInfo <- .ssUnionSnpInfo(handles)
    handle@chromPaths <- .ssUnionChromPaths(handles, fn)
    .asLdSketch(handle)
}

# Union the per-chromosome shard paths across handles. Two sketches trimmed
# from one genoMeta panel to different chromosomes share path/format/samples
# and differ ONLY here, so unioning is what lets their snpInfo rows keep
# routing: a row reaches its file by its own CHR. One chromosome mapping to
# two different files means the parts are not the same panel after all.
# @noRd
.ssUnionChromPaths <- function(handles, fn) {
    out <- character(0)
    for (h in handles) {
        cp <- h@chromPaths
        for (ch in names(cp)) {
            if (is_in(ch, names(out)) && !identical(out[[ch]], cp[[ch]])) {
                msg <- glue(
                    "{fn}: chromosome '{ch}' maps to two different genotype ",
                    "files across the inputs, so their LD sketches cannot ",
                    "be unioned."
                )
                abort(msg)
            }
            out[[ch]] <- cp[[ch]]
        }
    }
    out
}

# Every panel must read from the same file(s): the union keeps the first
# handle and widens its snpInfo, and each row's `fileIdx` addresses a position
# in THAT file, so rows from another panel would read the wrong variants.
# @noRd
.ssCheckSameSource <- function(handles, fn) {
    same <- map_lgl(handles, .ssSameGenotypeSource, first = handles[[1L]])
    if (all(same)) {
        return(invisible(NULL))
    }
    msg <- glue(
        "{fn}: inputs reference different genotype panels (input ",
        "{which(!same)[[1L]]} differs from the first). Their LD sketches ",
        "cannot be unioned."
    )
    abort(msg)
}

# @noRd
.ssSameGenotypeSource <- function(h, first) {
    # chromPaths is deliberately NOT compared: two sketches trimmed from the
    # same genoMeta panel to different chromosomes differ there and nowhere
    # else, and .ssUnionChromPaths() merges them (erroring on a genuine
    # conflict). Everything below must match for the union to mean anything --
    # snpInfo rows carry file POSITIONS, which only make sense against the
    # same files and the same sample axis.
    identical(h@path, first@path) &&
        identical(h@format, first@format) &&
        identical(h@nSamples, first@nSamples) &&
        identical(h@sampleIds, first@sampleIds)
}

# The parts' snpInfo rows, de-duplicated by variant id and returned in genomic
# order so the unioned panel reads like one built over the whole span. The
# first handle's data-frame class is preserved: a handle built with a base
# data.frame keeps one, because some downstream readers index snpInfo
# positionally and a tibble column is a one-column list there.
# @noRd
.ssUnionSnpInfo <- function(handles) {
    first <- getSnpInfo(handles[[1L]])
    si <- list_rbind(map(handles, .ssHandleSnpInfo))
    if (nrow(si) > 0L && all(is_in(c("SNP", "CHR", "BP"), names(si)))) {
        si <- si[!duplicated(si$SNP), , drop = FALSE]
        si <- si[order(canonChrom(si$CHR), as.integer(si$BP)), , drop = FALSE]
    }
    if (inherits(first, "tbl_df")) si else as.data.frame(si)
}

# @noRd
.ssHandleSnpInfo <- function(h) as_tibble(getSnpInfo(h))
