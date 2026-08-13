# =============================================================================
# QtlFineMappingResult S4 class
# -----------------------------------------------------------------------------
# DFrame-subclass collection keyed by the identity tuple (study, context,
# trait, method). Each row holds a FineMappingEntry payload. Class-level
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
    required <- c("study", "context", "trait", "method", "entry")
    missingCols <- setdiff(required, names(object))
    if (length(missingCols) > 0L) {
        return(paste("missing columns:", paste(missingCols, collapse = ", ")))
    }
    NULL
}

# Entry-column + region/traitPos + joint-column + tuple-uniqueness contract.
# @noRd
.qfmrCheckEntries <- function(object) {
    jointCols <- intersect(
        c("jointStudies", "jointContexts", "jointTraits"),
        names(object)
    )
    c(
        .qfmrCheckEntryLength(object),
        .qfmrCheckEntryTypes(object),
        .validateRegionColumn(object),
        .validateTraitPosColumn(object),
        .qfmrCheckJointCols(object, jointCols),
        .qfmrCheckTupleUniqueness(object, jointCols)
    )
}

# length(entry) must equal nrow.
# @noRd
.qfmrCheckEntryLength <- function(object) {
    if (length(object$entry) != nrow(object)) {
        return("length(entry) must equal nrow(.) for QtlFineMappingResult")
    }
    NULL
}

# Every `entry` element must be a FineMappingEntry.
# @noRd
.qfmrCheckEntryTypes <- function(object) {
    if (!all(map_lgl(object$entry, .qfmrIsEntry))) {
        return("every element of the `entry` column must be a FineMappingEntry")
    }
    NULL
}

# @noRd
.qfmrIsEntry <- function(e) {
    methods::is(e, "FineMappingEntry")
}

# Each present joint* column must be character.
# @noRd
.qfmrCheckJointCols <- function(object, jointCols) {
    unlist(compact(map(jointCols, .qfmrJointColError, object = object)))
}

# @noRd
.qfmrJointColError <- function(jc, object) {
    vals <- object[[jc]]
    if (!is.character(vals)) {
        return(sprintf(
            "'%s' column must be character (got %s)",
            jc,
            class(vals)[[1L]]
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
    keyDf <- as.data.frame(
        map(keyCols, .qfmrColOf, object = object),
        col.names = keyCols,
        stringsAsFactors = FALSE
    )
    if (anyDuplicated(keyDf)) {
        return(paste0(
            "(study, context, trait, method[, joint*]) tuple uniqueness ",
            "violated"
        ))
    }
    NULL
}

# @noRd
.qfmrColOf <- function(cn, object) {
    object[[cn]]
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
#'   each entry is a \code{FineMappingEntry}.
#'
#' Required columns: \code{study}, \code{method}, \code{entry}. The 2-tuple is
#' unique. The caller is expected to construct one \code{GwasFineMappingResult}
#' per LD block (no in-class block indexing).
#' @export

#' @title Create a QtlFineMappingResult Collection
#' @description Construct a \code{QtlFineMappingResult} DFrame-subclass
#'   collection from per-tuple vectors and a list of \code{FineMappingEntry}
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
#' @param entry List / \code{SimpleList} of \code{FineMappingEntry} objects.
#' @param jointStudies Optional character vector (length \code{length(study)})
#'   listing the semicolon-joined studies participating in each row's
#'   cross-study joint fit, or \code{NA_character_} for non-joint rows. When
#'   \code{NULL} (default) the column is omitted.
#' @param jointContexts Optional character vector for cross-context joints. Same
#'   shape as \code{jointStudies}.
#' @param jointTraits Optional character vector for cross-trait joints. Same
#'   shape as \code{jointStudies}.
#' @param region Optional \code{GRanges} (length \code{length(study)}) giving
#'   the genomic anchor of each row's trait (its own coordinates). Carried
#'   forward as provenance (e.g. for cTWAS LD-block placement); not part of the
#'   identity key. \code{NULL} (default) omits the column.
#' @param ldSketch An optional \code{GenotypeHandle} (the LD reference for
#'   RSS-derived fits), or \code{NULL} for individual-level fits.
#' @param traitPos Optional per-row trait genomic anchor (a \code{GRanges} or
#'   \code{NULL}), carried forward as provenance; not part of the identity key.
#'   \code{NULL} (default) omits the column.
#' @return A \code{QtlFineMappingResult} object.
#' @examples
#' tl <- data.frame(variant_id = paste0("chr1:", 100 * 1:3, ":A:G"),
#'   pip = c(0.9, 0.5, 0.1), cs = c(1L, 1L, NA))
#' fe <- FineMappingEntry(
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
    region = NULL,
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
        stop(
            "`study`, `context`, `trait`, `method`, and `entry` must all ",
            "have the same length."
        )
    }
    cols <- list(
        study = as.character(study),
        context = as.character(context),
        trait = as.character(trait),
        method = as.character(method),
        entry = S4Vectors::SimpleList(entry)
    )
    cols <- .qfmrAppendJointCols(
        cols,
        jointStudies,
        jointContexts,
        jointTraits,
        n
    )
    cols <- .appendRegionCol(cols, region, n)
    cols <- .appendTraitPosCol(cols, traitPos, n)
    df <- do.call(S4Vectors::DataFrame, c(cols, list(check.names = FALSE)))
    obj <- new("QtlFineMappingResult", df, ldSketch = ldSketch)
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
            stop("`", nm, "` must have the same length as `study`.")
        }
        cols[[nm]] <- as.character(val)
    }
    cols
}

#' @rdname getFineMappingResult
setMethod(
    "getFineMappingResult",
    "QtlFineMappingResult",
    function(x, study = NULL, context = NULL, trait = NULL, method = NULL) {
        idx <- .tupleSelectRow(
            x,
            study,
            context,
            trait,
            method,
            cls = "QtlFineMappingResult"
        )
        x$entry[[idx]]
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
        entry <- getFineMappingResult(x, study, context, trait, method)
        pip <- getPip(entry)
        if (isTRUE(returnList)) {
            nm <- sprintf(
                "%s|%s|%s|%s",
                as.character(x$study)[1L],
                as.character(x$context)[1L],
                as.character(x$trait)[1L],
                as.character(x$method)[1L]
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
        getFineMappingResult(x, study, context, trait, method)
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
        entry <- getFineMappingResult(x, study, context, trait, method)
        getCvResult(entry)
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
    cat(sprintf("QtlFineMappingResult: %d entries\n", nrow(object)))
    if (nrow(object) > 0L) {
        cat(sprintf(
            "  %d studies, %d contexts, %d traits, %d methods\n",
            length(unique(object$study)),
            length(unique(object$context)),
            length(unique(object$trait)),
            length(unique(object$method))
        ))
    }
    ldSrc <- if (is.null(object@ldSketch)) {
        "NULL (individual-level fit)"
    } else {
        sprintf("%s @ %s", object@ldSketch@format, object@ldSketch@path)
    }
    cat(sprintf("  LD sketch: %s\n", ldSrc))
})
