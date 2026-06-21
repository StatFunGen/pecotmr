# =============================================================================
# FineMappingResultBase S4 virtual class
# -----------------------------------------------------------------------------
# Shared parent of the QTL and GWAS fine-mapping result collections.
# Concrete subclasses (QtlFineMappingResult, GwasFineMappingResult) carry
# a DFrame of per-fit rows plus a shared ldSketch slot. Downstream
# pipelines dispatch on FineMappingResultBase for behaviors that apply to
# either flavour, and on the concrete subclass when the tuple shape
# matters.
# =============================================================================

#' @include allGenerics.R
#' @importFrom methods setClass setMethod new is validObject
NULL

#' @title Fine-Mapping Result Base Class
#' @description Virtual base class for fine-mapping result collections.
#'   Concrete subclasses (\code{QtlFineMappingResult},
#'   \code{GwasFineMappingResult}) carry a \code{DFrame} of per-fit rows
#'   and a shared \code{ldSketch} slot. Downstream pipelines should
#'   dispatch on \code{FineMappingResultBase} for behaviors that apply to
#'   either flavour, and on the concrete subclass when the tuple shape
#'   matters.
#' @slot ldSketch The LD reference \code{GenotypeHandle} the fits were
#'   computed against, or \code{NULL} when the fits were derived from
#'   individual-level data (no LD reference). Used downstream for
#'   cross-pipeline LD-sketch identity validation.
#' @export
setClass("FineMappingResultBase",
  contains = c("VIRTUAL", "DFrame"),
  representation(ldSketch = "ANY"))

#' @rdname getStudy
#' @export
setMethod("getStudy", "FineMappingResultBase",
          function(x) unique(as.character(x$study)))

#' @rdname getLdSketch
#' @export
setMethod("getLdSketch", "FineMappingResultBase",
          function(x, ...) x@ldSketch)

#' @rdname getMethodNames
#' @export
setMethod("getMethodNames", "FineMappingResultBase",
          function(x) unique(as.character(x$method)))

#' @rdname adjustPips
#' @export
setMethod("adjustPips", "FineMappingResultBase",
  function(x, keepVariants, ...) {
    if (nrow(x) == 0L) return(x)
    entries <- x@listData$entry
    for (i in seq_along(entries)) {
      adj <- tryCatch(adjustPips(entries[[i]], keepVariants, ...),
                      error = function(err) NULL)
      if (!is.null(adj)) entries[[i]] <- adj
    }
    x@listData$entry <- entries
    x
  })
