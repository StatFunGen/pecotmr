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

setClass("QtlFineMappingResult",
  contains = "FineMappingResultBase",
  validity = function(object) {
    errors <- character()
    required <- c("study", "context", "trait", "method", "entry")
    missingCols <- setdiff(required, names(object))
    if (length(missingCols) > 0L)
      errors <- c(errors, paste("missing columns:",
                                paste(missingCols, collapse = ", ")))
    if (length(errors) == 0L) {
      if (length(object$entry) != nrow(object))
        errors <- c(errors,
          "length(entry) must equal nrow(.) for QtlFineMappingResult")
      entryTypes <- vapply(object$entry,
                          function(e) methods::is(e, "FineMappingEntry"),
                          logical(1))
      if (!all(entryTypes))
        errors <- c(errors,
          "every element of the `entry` column must be a FineMappingEntry")
      jointCols <- intersect(
        c("jointStudies", "jointContexts", "jointTraits"), names(object))
      for (jc in jointCols) {
        vals <- object[[jc]]
        if (!is.character(vals))
          errors <- c(errors, sprintf(
            "'%s' column must be character (got %s)", jc, class(vals)[[1L]]))
      }
      keyCols <- c("study", "context", "trait", "method", jointCols)
      keyDf <- as.data.frame(object[, keyCols, drop = FALSE])
      if (anyDuplicated(keyDf))
        errors <- c(errors,
          "(study, context, trait, method[, joint*]) tuple uniqueness violated")
    }
    if (!is.null(object@ldSketch) &&
        !methods::is(object@ldSketch, "GenotypeHandle")) {
      errors <- c(errors,
        "'ldSketch' must be a GenotypeHandle or NULL")
    }
    if (length(errors) == 0L) TRUE else errors
  }
)

#' @title GWAS Fine-Mapping Result Collection
#' @description S4 collection of fine-mapping fits for one or more GWAS
#'   studies on a single LD block. Keyed by the identity tuple
#'   \code{(study, method)}; each entry is a \code{FineMappingEntry}.
#'
#'   Required columns: \code{study}, \code{method}, \code{entry}. The
#'   2-tuple is unique. The caller is expected to construct one
#'   \code{GwasFineMappingResult} per LD block (no in-class block
#'   indexing).
#' @export

#' @title Create a QtlFineMappingResult Collection
#' @description Construct a \code{QtlFineMappingResult} DFrame-subclass
#'   collection from per-tuple vectors and a list of
#'   \code{FineMappingEntry} payloads (one per tuple). The optional
#'   \code{ldSketch} slot records the LD reference used for RSS-derived
#'   fits; pass \code{NULL} (the default) for individual-level fits.
#' @param study Character vector of study identifiers (per tuple). Use
#'   the sentinel \code{"joint"} for rows produced by a cross-study
#'   joint fit.
#' @param context Character vector of context labels (per tuple). Use
#'   \code{"joint"} for rows produced by a cross-context joint fit.
#' @param trait Character vector of trait identifiers (per tuple). Use
#'   \code{"joint"} for rows produced by a cross-trait joint fit.
#' @param method Character vector of fine-mapping method names (per tuple).
#' @param entry List / \code{SimpleList} of \code{FineMappingEntry} objects.
#' @param jointStudies Optional character vector (length \code{length(study)})
#'   listing the semicolon-joined studies participating in each row's
#'   cross-study joint fit, or \code{NA_character_} for non-joint rows.
#'   When \code{NULL} (default) the column is omitted.
#' @param jointContexts Optional character vector for cross-context joints.
#'   Same shape as \code{jointStudies}.
#' @param jointTraits Optional character vector for cross-trait joints.
#'   Same shape as \code{jointStudies}.
#' @param ldSketch An optional \code{GenotypeHandle} (the LD reference for
#'   RSS-derived fits), or \code{NULL} for individual-level fits.
#' @return A \code{QtlFineMappingResult} object.
#' @export
QtlFineMappingResult <- function(study, context, trait, method, entry,
                                 jointStudies = NULL,
                                 jointContexts = NULL,
                                 jointTraits = NULL,
                                 ldSketch = NULL) {
  n <- length(study)
  if (length(context) != n || length(trait) != n || length(method) != n ||
      length(entry) != n) {
    stop("`study`, `context`, `trait`, `method`, and `entry` must all ",
         "have the same length.")
  }
  cols <- list(
    study   = as.character(study),
    context = as.character(context),
    trait   = as.character(trait),
    method  = as.character(method),
    entry   = S4Vectors::SimpleList(entry)
  )
  for (pair in list(c("jointStudies", "jointStudies"),
                    c("jointContexts", "jointContexts"),
                    c("jointTraits", "jointTraits"))) {
    val <- get(pair[[1L]])
    if (is.null(val)) next
    if (length(val) != n)
      stop("`", pair[[1L]], "` must have the same length as `study`.")
    cols[[pair[[2L]]]] <- as.character(val)
  }
  df <- do.call(S4Vectors::DataFrame,
                c(cols, list(check.names = FALSE)))
  obj <- new("QtlFineMappingResult", df, ldSketch = ldSketch)
  validObject(obj)
  obj
}

setMethod("getFineMappingResult", "QtlFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL) {
    idx <- .tupleSelectRow(x, study, context, trait, method,
                           cls = "QtlFineMappingResult")
    x$entry[[idx]]
  })

# Derived collection-level accessors (delegate to entry-level methods).

#' @rdname getPip
#' @export
setMethod("getPip", "QtlFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL,
           returnList = FALSE, ...) {
    entry <- getFineMappingResult(x, study, context, trait, method)
    pip <- getPip(entry)
    if (isTRUE(returnList)) {
      nm <- sprintf("%s|%s|%s|%s",
                    as.character(x$study)[1L],
                    as.character(x$context)[1L],
                    as.character(x$trait)[1L],
                    as.character(x$method)[1L])
      out <- list(); out[[nm]] <- pip
      return(out)
    }
    pip
  })

#' @rdname getCs
#' @export
setMethod("getCs", "QtlFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL,
           coverage = 0.95, ...) {
    entry <- getFineMappingResult(x, study, context, trait, method)
    getCs(entry, coverage = coverage)
  })

#' @rdname getCvResult
#' @export
setMethod("getCvResult", "QtlFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL, ...) {
    entry <- getFineMappingResult(x, study, context, trait, method)
    getCvResult(entry)
  })

#' @rdname getTopLoci
#' @export
setMethod("getTopLoci", "QtlFineMappingResult",
  function(x, type = c("data.frame", "GRanges"),
           signalCutoff = 0.025,
           study = NULL, context = NULL, trait = NULL, method = NULL,
           ...) {
    entry <- getFineMappingResult(x, study, context, trait, method)
    getTopLoci(entry, type = match.arg(type), signalCutoff = signalCutoff)
  })

#' @rdname getMarginalEffects
#' @export
setMethod("getMarginalEffects", "QtlFineMappingResult",
  function(x, maxPval = NULL,
           study = NULL, context = NULL, trait = NULL, method = NULL, ...) {
    entry <- getFineMappingResult(x, study, context, trait, method)
    getMarginalEffects(entry, maxPval = maxPval)
  })

#' @rdname getSusieFit
#' @export
setMethod("getSusieFit", "QtlFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL,
           ...) {
    entry <- getFineMappingResult(x, study, context, trait, method)
    getSusieFit(entry)
  })

#' @rdname getVariantIds
#' @export
setMethod("getVariantIds", "QtlFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL,
           ...) {
    entry <- getFineMappingResult(x, study, context, trait, method)
    getVariantIds(entry)
  })

#' @rdname getContexts
#' @export
setMethod("getContexts", "QtlFineMappingResult",
          function(x) unique(as.character(x$context)))

#' @rdname getTraits
#' @export
setMethod("getTraits", "QtlFineMappingResult",
          function(x) unique(as.character(x$trait)))
#' @export
setMethod("show", "QtlFineMappingResult", function(object) {
  cat(sprintf("QtlFineMappingResult: %d entries\n", nrow(object)))
  if (nrow(object) > 0L) {
    cat(sprintf("  %d studies, %d contexts, %d traits, %d methods\n",
                length(unique(object$study)),
                length(unique(object$context)),
                length(unique(object$trait)),
                length(unique(object$method))))
  }
  ldSrc <- if (is.null(object@ldSketch)) "NULL (individual-level fit)"
           else sprintf("%s @ %s",
                         object@ldSketch@format,
                         object@ldSketch@path)
  cat(sprintf("  LD sketch: %s\n", ldSrc))
})
