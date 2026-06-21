# =============================================================================
# TwasWeightsEntry S4 class
# -----------------------------------------------------------------------------
# Per-tuple TWAS weight payload: variant ids + per-variant weight
# vector/matrix, optional retained fit object, optional CV-performance
# metrics, and the standardized + dataType flags. One entry sits in
# every row of a TwasWeights collection.
# =============================================================================

#' @include allGenerics.R
NULL

setClass("TwasWeightsEntry",
  representation(
    variantIds    = "character",
    weights       = "ANY",
    fits          = "ANY",
    cvPerformance = "ANY",
    standardized  = "logical",
    dataType      = "ANY"
  ),
  validity = function(object) {
    errors <- character()
    if (length(object@standardized) != 1L)
      errors <- c(errors, "'standardized' must be a single logical")
    if (is.matrix(object@weights) &&
        nrow(object@weights) != length(object@variantIds))
      errors <- c(errors,
        "nrow(weights) must equal length(variantIds)")
    if (length(errors) == 0L) TRUE else errors
  }
)

# =============================================================================

# =============================================================================
# Summary Statistics collection classes (post-refactor; DFrame subclasses)
# =============================================================================

#' @title QTL Summary Statistics Collection
#' @description S4 collection of QTL summary statistics keyed by the
#'   identity tuple \code{(study, context, trait)}. Each entry holds a
#'   \code{GRanges} of summary statistics for that tuple. Class-level
#'   slots \code{ldSketch} (the LD reference \code{GenotypeHandle}) and
#'   \code{genome} (the genome build, a single character string) apply
#'   to every entry; the genome build must be uniform because all
#'   entries necessarily share the LD reference.
#'
#'   Required columns: \code{study}, \code{context}, \code{trait},
#'   \code{entry}. Optional columns include \code{varY} (numeric,
#'   per-tuple phenotype variance; \code{NA_real_} when unused). The
#'   3-tuple \code{(study, context, trait)} is unique. Each \code{entry}
#'   is a \code{GRanges} whose mcols carry the per-variant statistics
#'   (\code{SNP}, \code{A1}, \code{A2}, \code{Z}, \code{N}; plus
#'   optional \code{MAF}, \code{INFO}, \code{BETA}, \code{SE}, \code{P}).
#' @slot ldSketch A \code{GenotypeHandle} carrying the LD reference for
#'   downstream QC and RSS analysis.
#' @slot genome A single character string giving the genome build that
#'   the LD sketch and every entry are aligned to.
#' @title Summary Statistics Base Class
#' @description Virtual base class for summary-statistic collections.
#'   Concrete subclasses (\code{QtlSumStats}, \code{GwasSumStats}) carry
#'   a \code{DFrame} of per-entry rows plus shared slots \code{ldSketch},
#'   \code{genome}, and \code{qcInfo}. Downstream pipelines should
#'   dispatch on \code{SumStatsBase} for behaviors that apply to either
#'   flavour and on the concrete subclass when the tuple shape matters.
#' @slot ldSketch A \code{GenotypeHandle} carrying the LD reference for
#'   downstream QC and RSS analysis.
#' @slot genome A single character string giving the genome build that
#'   the LD sketch and every entry are aligned to.
#' @slot qcInfo A \code{list} recording which QC steps ran. Empty
#'   \code{list()} on construction; populated by \code{summaryStatsQc()}.
#'   Fine-mapping and TWAS-weights pipelines reject inputs where
#'   \code{length(getQcInfo(x)) == 0L} — the slot serves as both the
#'   gating flag and the audit trail.
#' @export

#' @title Create a TwasWeightsEntry Object
#' @description Construct a \code{TwasWeightsEntry} payload for one
#'   \code{(study, context, trait, method)} row of a \code{TwasWeights}
#'   collection.
#' @param variantIds Character vector of variant IDs.
#' @param weights Numeric vector or matrix.
#' @param fits Optional method-specific fit object.
#' @param cvPerformance Optional list of CV metrics.
#' @param standardized Logical (length 1).
#' @param dataType Optional data-type tag.
#' @return A \code{TwasWeightsEntry} object.
#' @export
TwasWeightsEntry <- function(variantIds, weights, fits = NULL,
                             cvPerformance = NULL, standardized = FALSE,
                             dataType = NULL) {
  obj <- new("TwasWeightsEntry",
             variantIds    = as.character(variantIds),
             weights       = weights,
             fits          = fits,
             cvPerformance = cvPerformance,
             standardized  = isTRUE(standardized),
             dataType      = dataType)
  validObject(obj)
  obj
}

#' @rdname getWeights
#' @export
setMethod("getWeights", "TwasWeightsEntry",
          function(x, ...) x@weights)

#' @rdname getVariantIds
#' @export
setMethod("getVariantIds", "TwasWeightsEntry",
          function(x, ...) x@variantIds)

#' @rdname getFits
#' @export
setMethod("getFits", "TwasWeightsEntry",
          function(x, ...) x@fits)

#' @rdname getCvPerformance
#' @export
setMethod("getCvPerformance", "TwasWeightsEntry",
          function(x, ...) x@cvPerformance)

#' @rdname getStandardized
#' @export
setMethod("getStandardized", "TwasWeightsEntry",
          function(x, ...) x@standardized)

#' @rdname getDataType
#' @export
setMethod("getDataType", "TwasWeightsEntry",
          function(x, ...) x@dataType)

#' @export
setMethod("show", "TwasWeightsEntry", function(object) {
  cat(sprintf("TwasWeightsEntry: %d variants, standardized=%s\n",
              length(object@variantIds), object@standardized))
  hasCv <- !is.null(object@cvPerformance)
  cat(sprintf("  CV performance: %s\n", hasCv))
})
