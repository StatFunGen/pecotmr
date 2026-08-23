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
#   * ldSketch   GenotypeHandle for the LD reference; required for the
#                LD-block-indexed susieRSS workflow.
# Methods that take per-row selectors accept (study, method) and ignore
# context/trait (GWAS has no per-tuple context or trait axis).
# =============================================================================

#' @include AllClasses.R tupleSelectors.R
NULL

setClass(
    "GwasFineMappingResult",
    contains = "FineMappingResultBase",
    validity = function(object) .validateGwasFineMappingResult(object)
)

# ---- GwasFineMappingResult validity helpers --------------------------------

# @noRd
.validateGwasFineMappingResult <- function(object) {
    errors <- .gfmrCheckRequiredCols(object)
    if (length(errors) == 0L) {
        errors <- .gfmrCheckEntries(object)
    }
    errors <- c(errors, .gfmrCheckLdSketch(object))
    if (length(errors) == 0L) TRUE else errors
}

# @noRd
.gfmrCheckRequiredCols <- function(object) {
    required <- c("study", "method")
    missingCols <- setdiff(required, .tupleColumnNames(object))
    if (length(missingCols) > 0L) {
        return(str_c("missing columns: ", str_flatten(missingCols, ", ")))
    }
    NULL
}

# @noRd
.gfmrCheckEntries <- function(object) {
    c(
        .gfmrCheckEntryLength(object),
        .validateTraitPosColumn(object),
        .gfmrCheckTupleUniqueness(object)
    )
}

# @noRd
.gfmrCheckEntryLength <- function(object) {
    # The variants and their topLoci ARE the elements now, so what is left to
    # check is that the fit payload columns are present and parallel.
    missingCols <- setdiff(
        c("susieFit", "cvResult"),
        .tupleColumnNames(object)
    )
    if (length(missingCols) > 0L) {
        return(str_c(
            "missing entry payload columns: ",
            str_flatten(missingCols, ", ")
        ))
    }
    NULL
}

# @noRd

# @noRd

# @noRd
.gfmrCheckTupleUniqueness <- function(object) {
    # Keyed on the element's RANGE rather than on a region_id label. The range
    # is intrinsic -- it cannot drift out of step with the variants the way a
    # stored label can, and it stays correct after subsetRegion() narrows an
    # element. Extract key columns directly rather than via `object[, keyCols]`:
    # column-subsetting preserves the class, and older S4Vectors revalidates
    # that intermediate, spuriously failing.
    if (length(object) == 0L) {
        return(NULL)
    }
    keyCols <- c("study", "method")
    keyTbl <- as_tibble(set_names(
        map(keyCols, .gfmrColOf, object = object),
        keyCols
    ))
    keyTbl$range <- .rtlRangeKeys(object)
    if (nrow(distinct(keyTbl)) < nrow(keyTbl)) {
        return("(study, method, range) tuple uniqueness violated")
    }
    NULL
}

# @noRd
.gfmrColOf <- function(cn, object) {
    .tupleColumn(object, cn)
}

# @noRd
.gfmrCheckLdSketch <- function(object) {
    if (
        !is.null(object@ldSketch) &&
            !methods::is(object@ldSketch, "GenotypeHandle")
    ) {
        return("'ldSketch' must be a GenotypeHandle or NULL")
    }
    NULL
}

# =============================================================================
# TWAS Weights
# =============================================================================

#' @title TWAS Weights Collection
#' @description S4 collection of TWAS weights keyed by the identity tuple
#'   \code{(study, context, trait, method)}. Each entry is a
#'   \code{TwasWeightsRow} carrying one method's weights for one
#'   trait/context/study. Implements the \code{DFrame}-subclass collection
#'   pattern.
#'
#' Required columns: \code{study}, \code{context}, \code{trait}, \code{method},
#' \code{entry}. Each \code{entry} is a \code{TwasWeightsRow}.
#'
#' Optional columns \code{jointStudies}, \code{jointContexts},
#' \code{jointTraits} appear when the collection contains rows produced by a
#' \code{jointSpecification}-driven joint fit. For such a row, the corresponding
#' identity-tuple column carries the sentinel \code{"joint"} and the joint
#' column lists the semicolon-joined members of the joined axis. For non-joint
#' rows the joint columns are \code{NA_character_}. Tuple uniqueness is enforced
#' jointly across the identity-tuple columns and any present joint columns.
#' @slot ldSketch The LD reference \code{GenotypeHandle} the weights were
#'   derived against, or \code{NULL} when the weights were learned from
#'   individual-level data. Used downstream for cross-pipeline LD-sketch
#'   identity validation.
#' @export

#' @title Create a GwasFineMappingResult Collection
#' @description Construct a \code{GwasFineMappingResult} collection from
#'   per-(study, method) tuples and a list of \code{FineMappingRow}
#'   payloads. The collection can represent either a single LD block (one row
#'   per (study, method)) or a genome-wide sweep across blocks (multiple rows
#'   per (study, method), each covering its own genomic range).
#'
#'   Rows are identified by (\code{study}, \code{method}) plus the element's
#'   own range, which is derived from the variants rather than stored, so it
#'   cannot drift out of step with them and stays correct after
#'   \code{\link{subsetRegion}}.
#' @param study Character vector of study identifiers (per tuple).
#' @param method Character vector of fine-mapping method names (per tuple).
#' @param entry List / \code{SimpleList} of \code{FineMappingRow} objects.
#' @param blockId Optional character vector (per tuple) keying the external LD
#'   block manifest, e.g. \code{"chr22_10516173_17414263"}. Provenance only,
#'   not part of the identity key. Supply it when the block BOUNDARIES matter:
#'   those are the one thing not recoverable from the variants, since adjacent
#'   blocks' spans leave gaps.
#' @param ldSketch An optional \code{GenotypeHandle}.
#' @param traitPos Optional per-row trait genomic anchor (a \code{GRanges} or
#'   \code{NULL}), carried forward as provenance; not part of the identity key.
#'   \code{NULL} (default) omits the column.
#' @return A \code{GwasFineMappingResult} object.
#' @examples
#' tl <- data.frame(variant_id = paste0("chr1:", 100 * 1:3, ":A:G"),
#'   pip = c(0.9, 0.5, 0.1), cs = c(1L, 1L, NA))
#' fe <- fineMappingRow(
#'   variantIds = tl$variant_id, susieFit = list(), topLoci = tl)
#' GwasFineMappingResult(study = "t1", method = "susie", entry = list(fe))
#' @export
GwasFineMappingResult <- function(
    study,
    method,
    entry,
    blockId = NULL,
    traitPos = NULL,
    ldSketch = NULL
) {
    n <- length(study)
    if (length(method) != n || length(entry) != n) {
        abort("`study`, `method`, and `entry` must all have the same length.")
    }
    entry <- map(entry, .asFmRowPayload)
    .checkRowPayloads(entry, "FineMappingRow", "fine-mapping")
    cols <- list(
        study = as.character(study),
        method = as.character(method),
        susieFit = S4Vectors::SimpleList(map(entry, getSusieFit)),
        cvResult = S4Vectors::SimpleList(map(entry, getCvResult))
    )
    cols <- .appendBlockIdCol(cols, blockId, n)
    cols <- .appendTraitPosCol(cols, traitPos, n)
    dfArgs <- c(cols, list(check.names = FALSE))
    # Each entry's variants become one ELEMENT, its topLoci that element's
    # inner mcols, and its fit/cv payload outer mcols. A multi-seqname entry
    # splits by chromosome with its metadata row replicated.
    split <- .rtlSplitBySeqname(map(entry, rowVariants))
    grl <- GenomicRanges::GRangesList(split$entry)
    md <- exec(S4Vectors::DataFrame, !!!dfArgs)
    mcols(grl) <- md[split$fromIdx, , drop = FALSE]
    obj <- new("GwasFineMappingResult", grl, ldSketch = ldSketch)
    validObject(obj)
    obj
}

# Internal: return integer row indices of `x` where every (column, value)
# pair in `keys` matches as.character(x[[column]]) == value. Shared
# building block for all of pecotmr's tuple-keyed row selectors and
# cache lookups (.tupleSelectRow, .qtlSumStatsSelectRow,

# GwasFineMappingResult has no context / trait columns; the generic
# returns NULL so callers can write generic code that handles either
# class without conditionals.
#' @rdname getContexts
#' @export
setMethod("getContexts", "GwasFineMappingResult", function(x) NULL)

#' @rdname getTraits
#' @export
setMethod("getTraits", "GwasFineMappingResult", function(x) NULL)

# Per-tuple lookup keyed by (study, method, blockId). The generic
# accepts the full set of selectors; context/trait args are ignored for
# GwasFineMappingResult. `region` (passed via `...`) is the per-block
# disambiguator for multi-block genome-wide collections.
#' @rdname getFineMappingResult
#' @export
setMethod(
    "getFineMappingResult",
    "GwasFineMappingResult",
    function(x, study = NULL, context = NULL, trait = NULL, method = NULL) {
        x[.tupleSelectRowGwasFmr(x, study, method)]
    }
)

#' @rdname getPip
#' @export
setMethod(
    "getPip",
    "GwasFineMappingResult",
    function(
        x,
        study = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        region = NULL,
        returnList = FALSE,
        ...
    ) {
        idx <- .tupleSelectRowGwasFmr(x, study, method, region)
        .fmrRowPip(.fmrRowParts(x, idx))
    }
)

# Row selector for the base delegating accessors (getCs / getTopLoci /
# getMarginalEffects / getSusieFit / getVariantIds live on
# FineMappingResultBase, AllClasses.R). The Gwas selector threads `region`.
setMethod(
    ".fmrSelectEntry",
    "GwasFineMappingResult",
    function(
        x,
        study = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        region = NULL,
        ...
    ) {
        .fmrRowParts(x, .tupleSelectRowGwasFmr(x, study, method, region))
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
    ldSrc <- if (is.null(object@ldSketch)) {
        "NULL"
    } else {
        glue("{object@ldSketch@format} @ {object@ldSketch@path}")
    }
    cat(glue("  LD sketch: {ldSrc}\n", .trim = FALSE))
})
