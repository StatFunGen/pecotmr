# =============================================================================
# QtlFineMappingResult S4 class
# -----------------------------------------------------------------------------
# DFrame-subclass collection keyed by the identity tuple (study, context,
# trait, method). Each row holds a FineMappingRow payload. Class-level
# slots:
#   * ldSketch   GenotypeHandle (NULL for individual-level fits, the
#                LD-sketch handle for RSS-derived fits).
# Optional columns jointStudies / jointContexts / jointTraits carry
# semicolon-joined member identities for cross-axis joint fits emitted by
# the joint dispatchers.
# =============================================================================

#' @include AllClasses.R tupleSelectors.R
NULL

setClass(
    "QtlFineMappingResult",
    contains = "FineMappingResultBase",
    validity = function(object) .validateQtlFineMappingResult(object)
)

# ---- QtlFineMappingResult validity helpers ---------------------------------

# Collect all contract violations (empty vector = valid). The per-entry checks
# run only once the required columns are present; the ldSketch check always
# runs.
# @noRd
.validateQtlFineMappingResult <- function(object) {
    errors <- .qfmrCheckRequiredCols(object)
    if (length(errors) == 0L) {
        errors <- .qfmrCheckEntries(object)
    }
    errors <- c(errors, .qfmrCheckLdSketch(object))
    if (length(errors) == 0L) TRUE else errors
}

# The study/context/trait/method/entry columns must be present.
# @noRd
.qfmrCheckRequiredCols <- function(object) {
    required <- c("study", "context", "trait", "method")
    missingCols <- setdiff(required, .tupleColumnNames(object))
    if (length(missingCols) > 0L) {
        return(str_c("missing columns: ", str_flatten(missingCols, ", ")))
    }
    NULL
}

# Entry-column + region/traitPos + joint-column + tuple-uniqueness contract.
# @noRd
.qfmrCheckEntries <- function(object) {
    jointCols <- intersect(
        c("jointStudies", "jointContexts", "jointTraits"),
        .tupleColumnNames(object)
    )
    c(
        .qfmrCheckEntryLength(object),
        .validateTraitPosColumn(object),
        .qfmrCheckJointCols(object, jointCols),
        .qfmrCheckTupleUniqueness(object, jointCols)
    )
}

# The variants and their topLoci ARE the elements now, so what is left
# to check is that the fit payload columns are present and parallel.
# @noRd
.qfmrCheckEntryLength <- function(object) {
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

# Each present joint* column must be character.
# @noRd
.qfmrCheckJointCols <- function(object, jointCols) {
    unlist(compact(map(jointCols, .qfmrJointColError, object = object)))
}

# @noRd
.qfmrJointColError <- function(jc, object) {
    vals <- .tupleColumn(object, jc)
    if (!is.character(vals)) {
        return(glue(
            "'{jc}' column must be character (got {class(vals)[[1L]]})"
        ))
    }
    NULL
}

# (study, context, trait, method[, joint*]) tuple uniqueness.
# @noRd
.qfmrCheckTupleUniqueness <- function(object, jointCols) {
    keyCols <- c("study", "context", "trait", "method", jointCols)
    # Extract key columns directly rather than via `object[, keyCols]`: column-
    # subsetting preserves the class while dropping the required `entry` column,
    # and older S4Vectors revalidates that intermediate, spuriously failing.
    keyTbl <- as_tibble(set_names(
        map(keyCols, .qfmrColOf, object = object),
        keyCols
    ))
    if (nrow(distinct(keyTbl)) < nrow(keyTbl)) {
        return(str_c(
            "(study, context, trait, method[, joint*]) tuple uniqueness ",
            "violated"
        ))
    }
    NULL
}

# @noRd
.qfmrColOf <- function(cn, object) {
    .tupleColumn(object, cn)
}

# ldSketch must be a GenotypeHandle or NULL.
# @noRd
.qfmrCheckLdSketch <- function(object) {
    if (
        !is.null(object@ldSketch) &&
            !methods::is(object@ldSketch, "GenotypeHandle")
    ) {
        return("'ldSketch' must be a GenotypeHandle or NULL")
    }
    NULL
}

#' @title GWAS Fine-Mapping Result Collection
#' @description S4 collection of fine-mapping fits for one or more GWAS studies
#'   on a single LD block. Keyed by the identity tuple \code{(study, method)};
#'   each entry is a \code{FineMappingRow}.
#'
#' Required columns: \code{study}, \code{method}, \code{entry}. The 2-tuple is
#' unique. The caller is expected to construct one \code{GwasFineMappingResult}
#' per LD block (no in-class block indexing).
#' @export

#' @title Create a QtlFineMappingResult Collection
#' @description Construct a \code{QtlFineMappingResult} DFrame-subclass
#'   collection from per-tuple vectors and a list of \code{FineMappingRow}
#'   payloads (one per tuple). The optional \code{ldSketch} slot records the LD
#'   reference used for RSS-derived fits; pass \code{NULL} (the default) for
#'   individual-level fits.
#' @param study Character vector of study identifiers (per tuple). Use the
#'   sentinel \code{"joint"} for rows produced by a cross-study joint fit.
#' @param context Character vector of context labels (per tuple). Use
#'   \code{"joint"} for rows produced by a cross-context joint fit.
#' @param trait Character vector of trait identifiers (per tuple). Use
#'   \code{"joint"} for rows produced by a cross-trait joint fit.
#' @param method Character vector of fine-mapping method names (per tuple).
#' @param entry List / \code{SimpleList} of \code{FineMappingRow} objects.
#' @param jointStudies Optional character vector (length \code{length(study)})
#'   listing the semicolon-joined studies participating in each row's
#'   cross-study joint fit, or \code{NA_character_} for non-joint rows. When
#'   \code{NULL} (default) the column is omitted.
#' @param jointContexts Optional character vector for cross-context joints. Same
#'   shape as \code{jointStudies}.
#' @param jointTraits Optional character vector for cross-trait joints. Same
#'   shape as \code{jointStudies}.
#' @param ldSketch An optional \code{GenotypeHandle} (the LD reference for
#'   RSS-derived fits), or \code{NULL} for individual-level fits.
#' @param traitPos Optional per-row trait genomic anchor (a \code{GRanges} or
#'   \code{NULL}), carried forward as provenance; not part of the identity key.
#'   \code{NULL} (default) omits the column.
#' @return A \code{QtlFineMappingResult} object.
#' @examples
#' tl <- data.frame(variant_id = paste0("chr1:", 100 * 1:3, ":A:G"),
#'   pip = c(0.9, 0.5, 0.1), cs = c(1L, 1L, NA))
#' fe <- fineMappingRow(
#'   variantIds = tl$variant_id, susieFit = list(), topLoci = tl)
#' QtlFineMappingResult(study = "s1", context = "brain", trait = "g1",
#'   method = "susie", entry = list(fe))
#' @export
QtlFineMappingResult <- function(
    study,
    context,
    trait,
    method,
    entry,
    jointStudies = NULL,
    jointContexts = NULL,
    jointTraits = NULL,
    traitPos = NULL,
    ldSketch = NULL
) {
    n <- length(study)
    if (
        length(context) != n ||
            length(trait) != n ||
            length(method) != n ||
            length(entry) != n
    ) {
        msg <- glue(
            "`study`, `context`, `trait`, `method`, and `entry` must all ",
            "have the same length."
        )
        abort(msg)
    }
    entry <- map(entry, .asFmRowPayload)
    .checkRowPayloads(entry, "FineMappingRow", "fine-mapping")
    cols <- list(
        study = as.character(study),
        context = as.character(context),
        trait = as.character(trait),
        method = as.character(method),
        susieFit = S4Vectors::SimpleList(map(entry, getSusieFit)),
        cvResult = S4Vectors::SimpleList(map(entry, getCvResult))
    )
    cols <- .qfmrAppendJointCols(
        cols,
        jointStudies,
        jointContexts,
        jointTraits,
        n
    )
    cols <- .appendTraitPosCol(cols, traitPos, n)
    dfArgs <- c(cols, list(check.names = FALSE))
    # Each entry's variants become one ELEMENT, its topLoci that element's
    # inner mcols, and its fit/cv payload outer mcols. A multi-seqname entry
    # splits by chromosome with its metadata row replicated.
    split <- .rtlSplitBySeqname(map(entry, rowVariants))
    grl <- GenomicRanges::GRangesList(split$entry)
    md <- exec(S4Vectors::DataFrame, !!!dfArgs)
    mcols(grl) <- md[split$fromIdx, , drop = FALSE]
    obj <- new("QtlFineMappingResult", grl, ldSketch = ldSketch)
    validObject(obj)
    obj
}

# Append any supplied joint* provenance columns (each must match length(study)).
# @noRd
.qfmrAppendJointCols <- function(
    cols,
    jointStudies,
    jointContexts,
    jointTraits,
    n
) {
    joints <- list(
        jointStudies = jointStudies,
        jointContexts = jointContexts,
        jointTraits = jointTraits
    )
    for (nm in names(joints)) {
        val <- joints[[nm]]
        if (is.null(val)) {
            next
        }
        if (length(val) != n) {
            msg <- glue("`{nm}` must have the same length as `study`.")
            abort(msg)
        }
        cols[[nm]] <- as.character(val)
    }
    cols
}

# The single row a (study, context, trait, method) selector pins.
# @noRd
.qfmrSelectRowIndex <- function(x, study, context, trait, method) {
    .tupleSelectRow(
        x,
        study,
        context,
        trait,
        method,
        cls = "QtlFineMappingResult"
    )
}

#' @rdname getFineMappingResult
setMethod(
    "getFineMappingResult",
    "QtlFineMappingResult",
    function(x, study = NULL, context = NULL, trait = NULL, method = NULL) {
        x[.qfmrSelectRowIndex(x, study, context, trait, method)]
    }
)

# Derived collection-level accessors (delegate to entry-level methods).

#' @rdname getPip
#' @export
setMethod(
    "getPip",
    "QtlFineMappingResult",
    function(
        x,
        study = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        returnList = FALSE,
        ...
    ) {
        pip <- .fmrRowPip(
            .fmrRowParts(
                x,
                .qfmrSelectRowIndex(x, study, context, trait, method)
            )
        )
        if (isTRUE(returnList)) {
            nm <- glue(
                "{as.character(x$study)[1L]}|{as.character(x$context)[1L]}|",
                "{as.character(x$trait)[1L]}|{as.character(x$method)[1L]}"
            )
            out <- list()
            out[[nm]] <- pip
            return(out)
        }
        pip
    }
)

# Row selector for the base delegating accessors (getCs / getTopLoci /
# getMarginalEffects / getSusieFit / getVariantIds now live on
# FineMappingResultBase, AllClasses.R).
setMethod(
    ".fmrSelectEntry",
    "QtlFineMappingResult",
    function(
        x,
        study = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        ...
    ) {
        .fmrRowParts(
            x,
            .qfmrSelectRowIndex(x, study, context, trait, method)
        )
    }
)

#' @rdname getCvResult
#' @export
setMethod(
    "getCvResult",
    "QtlFineMappingResult",
    function(
        x,
        study = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        ...
    ) {
        getCvResult(.fmrRowParts(
            x,
            .qfmrSelectRowIndex(x, study, context, trait, method)
        ))
    }
)

# getTopLoci / getMarginalEffects / getSusieFit / getVariantIds are defined once
# on FineMappingResultBase (AllClasses.R), routing through .fmrSelectEntry.

#' @rdname getContexts
#' @export
setMethod("getContexts", "QtlFineMappingResult", function(x) {
    unique(as.character(x$context))
})

#' @rdname getTraits
#' @export
setMethod("getTraits", "QtlFineMappingResult", function(x) {
    unique(as.character(x$trait))
})
#' @rdname show-methods
#' @export
setMethod("show", "QtlFineMappingResult", function(object) {
    cat(glue("QtlFineMappingResult: {nrow(object)} entries\n", .trim = FALSE))
    if (nrow(object) > 0L) {
        cat(glue(
            "  {n_distinct(object$study)} studies, ",
            "{n_distinct(object$context)} contexts, ",
            "{n_distinct(object$trait)} traits, ",
            "{n_distinct(object$method)} methods\n",
            .trim = FALSE
        ))
    }
    ldSrc <- if (is.null(object@ldSketch)) {
        "NULL (individual-level fit)"
    } else {
        glue("{object@ldSketch@format} @ {object@ldSketch@path}")
    }
    cat(glue("  LD sketch: {ldSrc}\n", .trim = FALSE))
})
