# =============================================================================
# CtwasResult S4 class
# -----------------------------------------------------------------------------
# DFrame-subclass collection of cTWAS runs, keyed by the identity tuple
# (gwasStudy, study, context, method). Each row holds a CtwasResultEntry
# payload (fine-mapping posteriors + group priors + region metadata) for one
# run. Unlike the QTL family, `trait` is NOT a key -- cTWAS is multi-gene per
# run, so genes live inside the payload. Optional joint columns (jointStudies,
# jointContexts) tag rows born from a multi-study / multi-context cTWAS run and
# participate in the uniqueness key, exactly as in the TwasWeights /
# FineMappingResult families.
# =============================================================================

#' @include AllGenerics.R tupleSelectors.R
NULL

setClass("CtwasResult", contains = "DFrame", validity = function(object) {
    .validateCtwasResult(object)
})

# ---- CtwasResult validity helpers ------------------------------------------

# Collect all contract violations (empty vector = valid). Entry checks run only
# once the required columns are present.
# @noRd
.validateCtwasResult <- function(object) {
    errors <- .ctwasResCheckRequiredCols(object)
    if (length(errors) == 0L) {
        errors <- .ctwasResCheckEntries(object)
    }
    if (length(errors) == 0L) TRUE else errors
}

# @noRd
.ctwasResCheckRequiredCols <- function(object) {
    required <- c("gwasStudy", "study", "context", "method", "entry")
    missingCols <- setdiff(required, names(object))
    if (length(missingCols) > 0L) {
        return(paste("missing columns:", paste(missingCols, collapse = ", ")))
    }
    NULL
}

# @noRd
.ctwasResCheckEntries <- function(object) {
    jointCols <- intersect(c("jointStudies", "jointContexts"), names(object))
    c(
        .ctwasResCheckEntryLength(object),
        .ctwasResCheckEntryTypes(object),
        .ctwasResCheckJointCols(object, jointCols),
        .ctwasResCheckTupleUniqueness(object, jointCols)
    )
}

# @noRd
.ctwasResCheckEntryLength <- function(object) {
    if (length(object$entry) != nrow(object)) {
        return("length(entry) must equal nrow(.) for CtwasResult")
    }
    NULL
}

# @noRd
.ctwasResCheckEntryTypes <- function(object) {
    if (!all(map_lgl(object$entry, .ctwasResIsEntry))) {
        return("every element of the `entry` column must be a CtwasResultEntry")
    }
    NULL
}

# @noRd
.ctwasResIsEntry <- function(e) {
    methods::is(e, "CtwasResultEntry")
}

# @noRd
.ctwasResCheckJointCols <- function(object, jointCols) {
    unlist(compact(map(jointCols, .ctwasResJointColError, object = object)))
}

# @noRd
.ctwasResJointColError <- function(jc, object) {
    if (!is.character(object[[jc]])) {
        return(sprintf(
            "'%s' column must be character (got %s)",
            jc,
            class(object[[jc]])[[1L]]
        ))
    }
    NULL
}

# @noRd
.ctwasResCheckTupleUniqueness <- function(object, jointCols) {
    keyCols <- c("gwasStudy", "study", "context", "method", jointCols)
    keyDf <- as.data.frame(
        map(keyCols, .ctwasResColOf, object = object),
        col.names = keyCols,
        stringsAsFactors = FALSE
    )
    if (anyDuplicated(keyDf)) {
        return(paste0(
            "(gwasStudy, study, context, method[, joint*]) ",
            "tuple uniqueness violated"
        ))
    }
    NULL
}

# @noRd
.ctwasResColOf <- function(cn, object) {
    object[[cn]]
}

#' @title Create a CtwasResult Collection
#' @description Construct a \code{CtwasResult} DFrame-subclass collection of
#'   cTWAS runs from per-tuple vectors and a list of
#'   \code{\link{CtwasResultEntry}} payloads (one per run). A single-context run
#'   is one row; a multi-context run emits per-context rows sharing the same
#'   \code{jointContexts} set and jointly-estimated group priors.
#' @param gwasStudy Character vector of GWAS (disease) study identifiers.
#' @param study Character vector of QTL study identifiers.
#' @param context Character vector of context labels (scalar per row).
#' @param method Character vector of weight-method names.
#' @param entry List / \code{SimpleList} of \code{CtwasResultEntry} objects.
#' @param jointStudies Optional character vector (length
#'   \code{length(gwasStudy)}) listing the joined studies for a multi-study run,
#'   \code{NA_character_} otherwise. \code{NULL} (default) omits the column.
#' @param jointContexts Optional character vector for multi-context runs. Same
#'   shape as \code{jointStudies}.
#' @return A \code{CtwasResult} object.
#' @examples
#' cre <- CtwasResultEntry(
#'   finemap = data.frame(id = c("g1", "g2"), susie_pip = c(0.9, 0.1)),
#'   susieAlpha = data.frame(id = c("g1", "g2"), alpha = c(0.9, 0.1)))
#' cr <- CtwasResult(gwasStudy = "gwas1", study = "s1", context = "brain",
#'   method = "susie", entry = list(cre))
#' cr
#' @export
CtwasResult <- function(
    gwasStudy,
    study,
    context,
    method,
    entry,
    jointStudies = NULL,
    jointContexts = NULL
) {
    n <- length(gwasStudy)
    if (
        length(study) != n ||
            length(context) != n ||
            length(method) != n ||
            length(entry) != n
    ) {
        stop(
            "`gwasStudy`, `study`, `context`, `method`, and `entry` must all ",
            "have the same length."
        )
    }
    cols <- list(
        gwasStudy = as.character(gwasStudy),
        study = as.character(study),
        context = as.character(context),
        method = as.character(method),
        entry = S4Vectors::SimpleList(entry)
    )
    for (nm in c("jointStudies", "jointContexts")) {
        val <- get(nm)
        if (is.null(val)) {
            next
        }
        if (length(val) != n) {
            stop("`", nm, "` must have the same length as `gwasStudy`.")
        }
        cols[[nm]] <- as.character(val)
    }
    df <- do.call(S4Vectors::DataFrame, c(cols, list(check.names = FALSE)))
    obj <- new("CtwasResult", df)
    validObject(obj)
    obj
}

#' @rdname getMethodNames
#' @export
setMethod("getMethodNames", "CtwasResult", function(x) {
    unique(as.character(x$method))
})

#' @rdname getStudy
#' @export
setMethod("getStudy", "CtwasResult", function(x) unique(as.character(x$study)))

#' @rdname getContexts
#' @export
setMethod("getContexts", "CtwasResult", function(x) {
    unique(as.character(x$context))
})

# Aggregate a per-entry table (finemap or susieAlpha) across all rows, prefixing
# each with the row's (gwasStudy, study, context, method) run identity.
# @noRd
.ctwasAggregateRows <- function(x, getter) {
    n <- nrow(x)
    if (n == 0L) {
        return(NULL)
    }
    parts <- lapply(seq_len(n), function(i) {
        tb <- getter(x$entry[[i]])
        if (is.null(tb) || nrow(tb) == 0L) {
            return(NULL)
        }
        tb <- as.data.frame(tb)
        id <- data.frame(
            gwasStudy = as.character(x$gwasStudy)[[i]],
            study = as.character(x$study)[[i]],
            context = as.character(x$context)[[i]],
            method = as.character(x$method)[[i]],
            stringsAsFactors = FALSE
        )
        cbind(id[rep(1L, nrow(tb)), , drop = FALSE], tb, row.names = NULL)
    })
    parts <- parts[!vapply(parts, is.null, logical(1L))]
    if (length(parts) == 0L) {
        return(NULL)
    }
    .rbindAligned(parts)
}

#' @rdname getFinemap
#' @export
setMethod("getFinemap", "CtwasResult", function(x, ...) {
    .ctwasAggregateRows(x, getFinemap)
})

#' @rdname getSusieAlpha
#' @export
setMethod("getSusieAlpha", "CtwasResult", function(x, ...) {
    .ctwasAggregateRows(x, getSusieAlpha)
})

#' @rdname show-methods
#' @export
setMethod("show", "CtwasResult", function(object) {
    cat(sprintf("CtwasResult: %d run(s)\n", nrow(object)))
    if (nrow(object) > 0L) {
        cat(sprintf(
            "  %d GWAS stud(y/ies), %d method(s): %s\n",
            length(unique(object$gwasStudy)),
            length(unique(object$method)),
            paste(unique(as.character(object$method)), collapse = ", ")
        ))
    }
})
