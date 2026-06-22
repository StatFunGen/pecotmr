# =============================================================================
# GwasFineMappingResult S4 class
# -----------------------------------------------------------------------------
# DFrame-subclass collection keyed by the identity tuple
# (study, method, region_id). Each row holds a FineMappingEntry payload
# for one GWAS study at one fine-mapping method over one LD block.
# Multiple rows per (study, method) are allowed when they differ on
# region_id — the genome-wide-across-blocks shape that
# qtlEnrichmentPipeline / colocPipeline expect when sweeping a study
# across LD blocks. Class-level slots:
#   * ldSketch   GenotypeHandle for the LD reference; required for the
#                LD-block-indexed susieRSS workflow.
# Methods that take per-row selectors accept (study, method) and ignore
# context/trait (GWAS has no per-tuple context or trait axis).
# =============================================================================

#' @include AllClasses.R tupleSelectors.R
NULL

setClass("GwasFineMappingResult",
  contains = "FineMappingResultBase",
  validity = function(object) {
    errors <- character()
    required <- c("study", "method", "region_id", "entry")
    missingCols <- setdiff(required, names(object))
    if (length(missingCols) > 0L)
      errors <- c(errors, paste("missing columns:",
                                paste(missingCols, collapse = ", ")))
    if (length(errors) == 0L) {
      if (length(object$entry) != nrow(object))
        errors <- c(errors,
          "length(entry) must equal nrow(.) for GwasFineMappingResult")
      entryTypes <- vapply(object$entry,
                          function(e) methods::is(e, "FineMappingEntry"),
                          logical(1))
      if (!all(entryTypes))
        errors <- c(errors,
          "every element of the `entry` column must be a FineMappingEntry")
      keyDf <- as.data.frame(object[, c("study", "method", "region_id")])
      if (anyDuplicated(keyDf))
        errors <- c(errors,
          "(study, method, region_id) tuple uniqueness violated")
    }
    if (!is.null(object@ldSketch) &&
        !methods::is(object@ldSketch, "GenotypeHandle")) {
      errors <- c(errors,
        "'ldSketch' must be a GenotypeHandle or NULL")
    }
    if (length(errors) == 0L) TRUE else errors
  }
)

# =============================================================================
# TWAS Weights
# =============================================================================

#' @title TWAS Weights Collection
#' @description S4 collection of TWAS weights keyed by the identity tuple
#'   \code{(study, context, trait, method)}. Each entry is a
#'   \code{TwasWeightsEntry} carrying one method's weights for one
#'   trait/context/study. Implements the \code{DFrame}-subclass
#'   collection pattern.
#'
#'   Required columns: \code{study}, \code{context}, \code{trait},
#'   \code{method}, \code{entry}. Each \code{entry} is a
#'   \code{TwasWeightsEntry}.
#'
#'   Optional columns \code{jointStudies}, \code{jointContexts},
#'   \code{jointTraits} appear when the collection contains rows
#'   produced by a \code{jointSpecification}-driven joint fit. For such
#'   a row, the corresponding identity-tuple column carries the sentinel
#'   \code{"joint"} and the joint column lists the semicolon-joined
#'   members of the joined axis. For non-joint rows the joint columns
#'   are \code{NA_character_}. Tuple uniqueness is enforced jointly
#'   across the identity-tuple columns and any present joint columns.
#' @slot ldSketch The LD reference \code{GenotypeHandle} the weights were
#'   derived against, or \code{NULL} when the weights were learned from
#'   individual-level data. Used downstream for cross-pipeline
#'   LD-sketch identity validation.
#' @export

#' @title Create a GwasFineMappingResult Collection
#' @description Construct a \code{GwasFineMappingResult} DFrame-subclass
#'   collection from per-(study, method, region_id) tuples and a list of
#'   \code{FineMappingEntry} payloads. The collection can represent
#'   either a single LD block (one row per (study, method)) or a
#'   genome-wide sweep across blocks (multiple rows per (study, method),
#'   each tagged with its own region_id).
#' @param study Character vector of study identifiers (per tuple).
#' @param method Character vector of fine-mapping method names (per tuple).
#' @param entry List / \code{SimpleList} of \code{FineMappingEntry} objects.
#' @param region_id Character vector of LD-block identifiers (per
#'   tuple). When omitted (\code{NULL}), defaults to a per-row synthetic
#'   id (\code{"region_1"}, \code{"region_2"}, ...) so the
#'   (\code{study}, \code{method}, \code{region_id}) triple is unique.
#'   Supplying meaningful labels (e.g. \code{"chr22_10516173_17414263"})
#'   is preferred for downstream consumers that join on region.
#' @param ldSketch An optional \code{GenotypeHandle}.
#' @return A \code{GwasFineMappingResult} object.
#' @export
GwasFineMappingResult <- function(study, method, entry,
                                  region_id = NULL,
                                  ldSketch = NULL) {
  n <- length(study)
  if (length(method) != n || length(entry) != n) {
    stop("`study`, `method`, and `entry` must all have the same length.")
  }
  if (is.null(region_id)) {
    region_id <- paste0("region_", seq_len(n))
  } else if (length(region_id) != n) {
    stop("`region_id` must have the same length as `study` (got ",
         length(region_id), " vs ", n, ").")
  }
  cols <- list(
    study     = as.character(study),
    method    = as.character(method),
    region_id = as.character(region_id),
    entry     = S4Vectors::SimpleList(entry)
  )
  df <- do.call(S4Vectors::DataFrame,
                c(cols, list(check.names = FALSE)))
  obj <- new("GwasFineMappingResult", df, ldSketch = ldSketch)
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

# Per-tuple lookup keyed by (study, method, region_id). The generic
# accepts the full set of selectors; context/trait args are ignored for
# GwasFineMappingResult. `region` (passed via `...`) is the per-block
# disambiguator for multi-block genome-wide collections.
#' @rdname getFineMappingResult
#' @export
setMethod("getFineMappingResult", "GwasFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL) {
    idx <- .tupleSelectRowGwasFmr(x, study, method)
    x$entry[[idx]]
  })

#' @rdname getPip
#' @export
setMethod("getPip", "GwasFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL,
           region = NULL, returnList = FALSE, ...) {
    idx <- .tupleSelectRowGwasFmr(x, study, method, region)
    getPip(x$entry[[idx]])
  })

#' @rdname getCs
#' @export
setMethod("getCs", "GwasFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL,
           region = NULL, ...) {
    idx <- .tupleSelectRowGwasFmr(x, study, method, region)
    getCs(x$entry[[idx]])
  })

#' @rdname getTopLoci
#' @export
setMethod("getTopLoci", "GwasFineMappingResult",
  function(x, type = c("data.frame", "GRanges"),
           signalCutoff = 0.025,
           study = NULL, context = NULL, trait = NULL, method = NULL,
           region = NULL, ...) {
    idx <- .tupleSelectRowGwasFmr(x, study, method, region)
    getTopLoci(x$entry[[idx]], type = match.arg(type),
               signalCutoff = signalCutoff)
  })

#' @rdname getMarginalEffects
#' @export
setMethod("getMarginalEffects", "GwasFineMappingResult",
  function(x, maxPval = NULL,
           study = NULL, context = NULL, trait = NULL, method = NULL,
           region = NULL, ...) {
    idx <- .tupleSelectRowGwasFmr(x, study, method, region)
    getMarginalEffects(x$entry[[idx]], maxPval = maxPval)
  })

#' @rdname getSusieFit
#' @export
setMethod("getSusieFit", "GwasFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL,
           region = NULL, ...) {
    idx <- .tupleSelectRowGwasFmr(x, study, method, region)
    getSusieFit(x$entry[[idx]])
  })

#' @rdname getVariantIds
#' @export
setMethod("getVariantIds", "GwasFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL,
           region = NULL, ...) {
    idx <- .tupleSelectRowGwasFmr(x, study, method, region)
    getVariantIds(x$entry[[idx]])
  })


#' @export
setMethod("show", "GwasFineMappingResult", function(object) {
  cat(sprintf("GwasFineMappingResult: %d entries\n", nrow(object)))
  if (nrow(object) > 0L) {
    cat(sprintf("  %d studies, %d methods\n",
                length(unique(object$study)),
                length(unique(object$method))))
  }
  ldSrc <- if (is.null(object@ldSketch)) "NULL"
           else sprintf("%s @ %s",
                         object@ldSketch@format,
                         object@ldSketch@path)
  cat(sprintf("  LD sketch: %s\n", ldSrc))
})

