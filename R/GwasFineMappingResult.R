# =============================================================================
# GwasFineMappingResult S4 class
# -----------------------------------------------------------------------------
# Collection keyed by the identity tuple (study, method) plus the element's
# own genomic RANGE. Each row holds a FineMappingRow payload
# for one GWAS study at one fine-mapping method over one LD block.
# Multiple rows per (study, method) are allowed when they differ on
# range -- the genome-wide-across-blocks shape that
# qtlEnrichmentPipeline / colocPipeline expect when sweeping a study
# across LD blocks. Class-level slots:
#   * ldReference   the LD references the fits used, each stored once; the
#                ldReference mcols column names each row's. Required for the
#                LD-block-indexed susieRSS workflow.
# Methods that take per-row selectors accept (study, method) and ignore
# context/trait (GWAS has no per-tuple context or trait axis).
# =============================================================================

#' @include AllClasses.R tupleSelectors.R
NULL

#' @title GWAS Fine-Mapping Result Collection
#' @description S4 collection of fine-mapping fits for one or more GWAS studies
#'   on a single LD block. Keyed by the identity tuple \code{(study, method)};
#'   each entry is a \code{FineMappingRow}.
#'
#' Required columns: \code{study}, \code{method}, \code{entry}. The 2-tuple is
#' unique. The caller is expected to construct one \code{GwasFineMappingResult}
#' per LD block (no in-class block indexing).
#' @export
setClass(
    "GwasFineMappingResult",
    contains = "FineMappingResultBase",
    validity = function(object) .validateGwasFineMappingResult(object)
)

# ---- GwasFineMappingResult validity helpers --------------------------------

# @noRd
#' @importFrom checkmate makeAssertCollection assertNames checkNames
.validateGwasFineMappingResult <- function(object) {
    coll <- makeAssertCollection()
    assertNames(
        .tupleColumnNames(object),
        must.include = c("study", "method"),
        what = "colnames",
        .var.name = "mcols",
        add = coll
    )
    # The checks below read those columns; running them on an object missing
    # them reports the consequence rather than the cause.
    if (!coll$isEmpty()) {
        return(coll$getMessages())
    }
    coll$push(.gfmrCheckEntries(object))
    coll$getMessages()
}

# @noRd
.gfmrCheckEntries <- function(object) {
    c(
        .gfmrCheckEntryLength(object),
        .validateTraitPosColumn(object),
        .validateRegionColumn(object),
        .gfmrCheckTupleUniqueness(object)
    )
}

# @noRd
.gfmrCheckEntryLength <- function(object) {
    # The variants and their topLoci ARE the elements now, so what is left to
    # check is that the fit payload columns are present and parallel.
    res <- checkNames(
        .tupleColumnNames(object),
        must.include = c("susieFit", "cvResult"),
        what = "colnames"
    )
    if (isTRUE(res)) {
        return(NULL)
    }
    str_c("missing entry payload columns: ", res)
}

# @noRd
.gfmrCheckTupleUniqueness <- function(object) {
    # Keyed on the element's RANGE rather than on a region_id label. The range
    # is intrinsic -- it cannot drift out of step with the variants the way a
    # stored label can, and it stays correct after an element is narrowed.
    # Extract key columns directly rather than via `object[, keyCols]`:
    # column-subsetting preserves the class, and older S4Vectors revalidates
    # that intermediate, spuriously failing.
    if (length(object) == 0L) {
        return(NULL)
    }
    keyCols <- c("study", "method")
    keyTbl <- as_tibble(set_names(
        map(keyCols, .gfmrColOf, object = object),
        keyCols
    )) |>
        mutate(range = .keyEmptyElements(object, .rtlRangeKeys(object)))
    if (nrow(distinct(keyTbl)) < nrow(keyTbl)) {
        return("(study, method, range) tuple uniqueness violated")
    }
    NULL
}

# @noRd
.gfmrColOf <- function(cn, object) {
    .tupleColumn(object, cn)
}


#' @title Create a GwasFineMappingResult Collection
#' @description Construct a \code{GwasFineMappingResult} collection from
#'   per-(study, method) tuples and a list of \code{FineMappingRow}
#'   payloads. The collection can represent either a single LD block (one row
#'   per (study, method)) or a genome-wide sweep across blocks (multiple rows
#'   per (study, method), each covering its own genomic range).
#'
#'   Rows are identified by (\code{study}, \code{method}) plus the element's
#'   own range, which is derived from the variants rather than stored, so it
#'   cannot drift out of step with them and stays correct when an element's
#'   variants are narrowed.
#' @param studyName Character vector of study identifiers (per tuple).
#' @param method Character vector of fine-mapping method names (per tuple).
#' @param entry List / \code{SimpleList} of \code{FineMappingRow} objects.
#' @param blockId Optional character vector (per tuple) naming each row's
#'   block, e.g. \code{"chr22_10516173_17414263"}. Provenance only, not part
#'   of the identity key, and a name only: a block's boundaries, which the
#'   variants cannot give since adjacent blocks' spans leave gaps, go in
#'   \code{region}. Where it is \code{NULL}, or \code{NA} for a row, the
#'   row's \code{blockId} is its own span, as \code{"chrom_start_end"}, so
#'   every row is named by its block.
#' @param dataType Optional assay / modality label ("bulk_rnaseq",
#'   "methylation"), recorded as a per-row provenance column. Length 1 or
#'   one per row; \code{NULL} (default) omits the column.
#' @param dataSource Optional per-row label for the kind of input the fit saw,
#'   \code{"individual"} or \code{"sumstat"}. Where it is absent, the result
#'   writers read a row's kind from its LD reference: a row with none saw
#'   individual-level data.
#' @param weightStandardized Optional per-row flag: whether the posterior
#'   effect is on the STANDARDIZED genotype scale. The posterior effect IS the
#'   prediction weight, so it carries the same scale question as a TWAS
#'   weight. \code{\link{fineMappingPipeline}} always passes \code{TRUE}:
#'   it standardizes genotypes before fitting, and the posterior mean is
#'   reported on the scale of the \code{X} the fit was given. The argument
#'   exists for a collection assembled by hand, where only the caller knows.
#' @param ldReference The LD reference the fits were computed against -- a
#'   genotype panel (see \code{\link{readGenotypes}}), a
#'   \code{GenotypeMixture} or an \code{LdMatrix} -- used for every row, or
#'   a list with one per row; \code{NULL} for none. Each distinct reference is
#'   stored once.
#' @param traitPos Optional per-row trait genomic anchor (a \code{GRanges} or
#'   \code{NULL}), carried forward as provenance; not part of the identity key.
#'   \code{NULL} (default) omits the column.
#' @param region Optional per-row analysis WINDOW (a \code{GRanges} or
#'   \code{NULL}) -- the interval each row's fit was run over, carried
#'   forward as provenance; not part of the identity key. It is recorded as
#'   given, so each range must contain its row's variants.
#'   \code{\link{fineMappingPipeline}} always sets it, from the window of the
#'   input row. \code{NULL} (default) omits the column, and
#'   \code{\link{genomicRegion}} then derives the window from the span of the
#'   variants present.
#' @param metadata A list of free-form metadata about the collection,
#'   as \code{metadata()} holds on any Bioconductor vector. Default
#'   \code{list()}.
#' @return A \code{GwasFineMappingResult} object. Its logical
#'   \code{regionGiven} column says, per row, whether the window was given
#'   or is the span of the variants; \code{\link{fineMappingPipeline}}
#'   carries it over from the input. See \code{\linkS4class{GwasSumStats}}
#'   for how the pipelines use it.
#' @examples
#' tl <- data.frame(variant_id = paste0("chr1:", 100 * 1:3, ":A:G"),
#'   pip = c(0.9, 0.5, 0.1), cs = c(1L, 1L, NA))
#' fe <- FineMappingRow(
#'   variantIds = tl$variant_id, susieFit = list(), topLoci = tl)
#' GwasFineMappingResult(studyName = "t1", method = "susie", entry = list(fe))
#' @export
GwasFineMappingResult <- function(
    studyName,
    method,
    entry,
    blockId = NULL,
    traitPos = NULL,
    region = NULL,
    dataType = NULL,
    dataSource = NULL,
    weightStandardized = NULL,
    ldReference = NULL,
    metadata = list()
) {
    n <- length(studyName)
    if (length(method) != n || length(entry) != n) {
        abort(
            "`studyName`, `method`, and `entry` must all have the same length."
        )
    }
    entry <- map(entry, .asFmRowPayload)
    .checkRowPayloads(entry, "FineMappingRow", "fine-mapping")
    .gfmrRefuseSplitFits(entry)
    cols <- list(
        study = as.character(studyName),
        method = as.character(method),
        susieFit = S4Vectors::SimpleList(map(entry, susieFit)),
        cvResult = S4Vectors::SimpleList(map(entry, cvResult)),
        derivedColumns = IRanges::CharacterList(
            map(entry, .fmRowDerivedColumns)
        )
    ) |>
        .appendBlockIdCol(blockId %||% rep(NA_character_, n), n) |>
        .appendTraitPosCol(traitPos, n) |>
        .appendRegionCol(region, n) |>
        .appendFmProvenanceCols(dataType, dataSource, weightStandardized, n)
    .withMetadata(.gfmrBuild(entry, cols, region, ldReference), metadata)
}

# The validated collection from its entries and metadata columns. Each
# entry's variants become one ELEMENT, its topLoci that element's inner
# mcols, and its fit/cv payload outer mcols. A multi-seqname entry -- only a
# fit-less one gets this far -- splits by chromosome with its metadata row
# replicated.
# @noRd
.gfmrBuild <- function(entry, cols, region, ldReference) {
    n <- length(entry)
    dict <- .ldRefDict(.ldRefPerTuple(ldReference, n))
    dfArgs <- c(cols, list(ldReference = dict$keys, check.names = FALSE))
    split <- .rtlSplitBySeqname(map(entry, variants))
    pieces <- .rtlBindElements(split$entry)
    grl <- S4Vectors::`mcols<-`(
        pieces,
        value = .gfmrPieceMetadata(
            exec(S4Vectors::DataFrame, !!!dfArgs),
            split,
            .gfmrGivenRegions(region, n),
            pieces
        )
    )
    new("GwasFineMappingResult", .rtlNamed(grl), ldReference = dict$refs)
}

# The window each entry was given: its `region`, when the caller passed one.
# A `blockId` is a name, and gives no window.
# @noRd
.gfmrGivenRegions <- function(region, n) {
    if (is.null(region)) {
        return(list(ranges = .placeholderRanges(n), given = rep(FALSE, n)))
    }
    list(ranges = .bareRanges(region), given = rep(TRUE, n))
}

# Each piece's metadata row: its entry's, replicated across the pieces a
# multi-chromosome entry splits into, and a blockId for every piece -- the
# caller's, else the piece's own span as "chrom_start_end" -- so each row is
# named by its block whether or not the caller named it. An empty piece spans
# nothing and keeps NA. `regionGiven` records which pieces were given their
# range. A row given none stores no `region`, where genomicRegion() reads its
# span, unless another row was given one and the column has to exist.
# @noRd
.gfmrPieceMetadata <- function(md, split, given, pieces) {
    out <- md[split$fromIdx, , drop = FALSE]
    out$blockId <- coalesce(as.character(out$blockId), .rtlRangeKeys(pieces))
    ranges <- given$ranges[split$fromIdx]
    has <- given$given[split$fromIdx]
    .refuseUncontainedRegions(ranges, has, pieces, "GwasFineMappingResult")
    cols <- .regionColumns(ranges, has, .rowSpanRanges(pieces))
    if (any(cols$regionGiven)) {
        out$region <- cols$region
    }
    out$regionGiven <- cols$regionGiven
    out
}

# A multi-chromosome entry is split into one row per chromosome, its metadata
# copied onto each. Variants and their top-loci split cleanly, but a fit does
# not: its effects span every variant it was given, so each piece would carry
# the whole fit, and coloc and the enrichment would read it once per piece. An
# entry carrying a fit must therefore lie on one chromosome -- one LD block,
# as fineMappingPipeline's entries do.
# @noRd
.gfmrRefuseSplitFits <- function(entry) {
    multi <- .multiChromosomeEntries(entry)
    bad <- multi[map_lgl(entry[multi], .fmRowCarriesFit)]
    if (length(bad) == 0L) {
        return(invisible(NULL))
    }
    abort(glue(
        "GwasFineMappingResult: `entry` element(s) ",
        "{str_flatten(utils::head(bad, 5L), ', ')} span more than one ",
        "chromosome but carry a single fit, which cannot be split across ",
        "rows. Fine-map each LD block separately and pass one entry per block."
    ))
}

# @noRd
.fmRowCarriesFit <- function(e) {
    length(susieFit(e)) > 0L || !is.null(cvResult(e))
}

# Internal: return integer row indices of `x` where every (column, value)
# pair in `keys` matches as.character(x[[column]]) == value. Shared
# building block for all of pecotmr's tuple-keyed row selectors and
# cache lookups (.tupleSelectRow, .qtlSumStatsSelectRow,

# GwasFineMappingResult has no context / trait columns; the generic
# returns NULL so callers can write generic code that handles either
# class without conditionals.
#' @rdname contexts
#' @export
setMethod("contexts", "GwasFineMappingResult", function(x) NULL)

#' @rdname traitNames
#' @export
setMethod("traitNames", "GwasFineMappingResult", function(x) NULL)

# One fit, keyed by (study, method): a study fine-mapped block by block is one
# row per LD block, so the fit is all of them; `blockId` picks one block by
# name. The generic accepts
# the full set of selectors; `context` / `trait` are QTL-only and a
# GwasFineMappingResult refuses them rather than ignoring them, the same way
# sumStats(GwasSumStats) does.
#' @rdname fineMappingResult
#' @export
setMethod(
    "fineMappingResult",
    "GwasFineMappingResult",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        blockId = NULL
    ) {
        .gwasFmrRefuseQtlSelectors(context, trait)
        x[.tupleSelectRowsGwasFmr(x, studyName, method, blockId)]
    }
)

#' @rdname pip
#' @export
setMethod(
    "pip",
    "GwasFineMappingResult",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        blockId = NULL
    ) {
        .gwasFmrRefuseQtlSelectors(context, trait)
        idx <- .tupleSelectRowGwasFmr(x, studyName, method, blockId)
        .fmrRowPip(.fmrRowParts(x, idx))
    }
)

# `context` / `trait` belong to the QTL axis. Naming them keeps the shared
# accessor signature, but a GWAS collection has no such axis, so asking for
# one is a mistake rather than a no-op.
# @noRd
.gwasFmrRefuseQtlSelectors <- function(context, trait) {
    # NA counts as unset, not as a selection: callers that thread a whole
    # (study, context, trait, method) record fill the axes a collection does
    # not have with NA, and that is an absence rather than a request.
    given <- c(
        context = .gwasFmrSelectorGiven(context),
        trait = .gwasFmrSelectorGiven(trait)
    )
    if (!any(given)) {
        return(invisible(NULL))
    }
    bad <- names(given)[given]
    abort(glue(
        "GwasFineMappingResult has no context or trait axis, so ",
        "{str_flatten(sprintf('`%s`', bad), ' and ')} ",
        "{if (length(bad) == 1L) 'does' else 'do'} not select anything. ",
        "Select with `studyName`, `method` and `blockId`."
    ))
}

# @noRd
.gwasFmrSelectorGiven <- function(v) {
    !is.null(v) && !all(is.na(v))
}

# Row selector for the base delegating accessors (credibleSets / topLoci /
# marginalEffects / susieFit / variantIds live on
# FineMappingResultBase, AllClasses.R). The Gwas selector threads `blockId`.
setMethod(
    ".fmrSelectEntry",
    "GwasFineMappingResult",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        blockId = NULL
    ) {
        # No refusal here, unlike the public accessors: this is the blind
        # delegation point for the shared base accessors, and callers such
        # as the VCF writer thread a whole selector record through it
        # without knowing the class. The user-facing entry points are where
        # a QTL-only selector is a mistake worth reporting.
        .fmrRowParts(
            x,
            .tupleSelectRowGwasFmr(x, studyName, method, blockId)
        )
    }
)


#' @rdname show-methods
#' @export
setMethod("show", "GwasFineMappingResult", function(object) {
    cat(glue("GwasFineMappingResult: {nrow(object)} entries\n", .trim = FALSE))
    if (nrow(object) > 0L) {
        cat(glue(
            "  {n_distinct(object$study)} studies, ",
            "{n_distinct(object$method)} methods\n",
            .trim = FALSE
        ))
    }
    ldSrc <- .ldRefsLabel(ldReference(object), "NULL")
    cat(glue("  LD reference: {ldSrc}\n", .trim = FALSE))
})
