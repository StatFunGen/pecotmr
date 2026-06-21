# =============================================================================
# SumStatsBase S4 virtual class
# -----------------------------------------------------------------------------
# Shared parent of the QTL and GWAS summary statistics collections.
# Concrete subclasses (QtlSumStats, GwasSumStats) inherit from DFrame and
# share the ldSketch / genome / qcInfo slots. Class-specific accessors
# (getZ / getN / getMaf / nSnps / subsetChr / getVarY / getSumStats)
# stay on the concrete subclass because they rely on the tuple shape
# (3-tuple for QtlSumStats, 1-tuple for GwasSumStats).
# =============================================================================

#' @include allGenerics.R GenotypeHandle.R
#' @importFrom methods setClass setMethod new is validObject
NULL

#' @title Summary Statistics Base Class
#' @description Virtual base class for QTL and GWAS summary statistics
#'   collections. Concrete subclasses (\code{QtlSumStats},
#'   \code{GwasSumStats}) inherit from \code{DFrame} and share the
#'   \code{ldSketch} / \code{genome} / \code{qcInfo} slots.
#' @slot ldSketch The \code{GenotypeHandle} the QC pipeline harmonized
#'   against. Required: \code{summaryStatsQc()} sets it.
#' @slot genome Character, genome build label.
#' @slot qcInfo A \code{list} recording which QC steps ran. Empty
#'   \code{list()} on construction; populated by \code{summaryStatsQc()}
#'   with a per-step audit record (filter names, drop counts, liftover
#'   target, RAISS settings, etc.). Fine-mapping and TWAS-weights
#'   pipelines reject inputs where \code{length(getQcInfo(x)) == 0L} — the
#'   slot serves as both the gating flag and the audit trail.
#' @export
setClass("SumStatsBase",
  contains = c("VIRTUAL", "DFrame"),
  representation(
    ldSketch = "GenotypeHandle",
    genome   = "character",
    qcInfo   = "list"
  ))

#' @rdname getGenome
#' @export
setMethod("getGenome", "SumStatsBase", function(x, ...) x@genome)

#' @rdname getQcInfo
#' @export
setMethod("getQcInfo", "SumStatsBase", function(x, ...) x@qcInfo)

#' @rdname getLdSketch
#' @export
setMethod("getLdSketch", "SumStatsBase", function(x, ...) x@ldSketch)

#' @rdname getStudy
#' @export
setMethod("getStudy", "SumStatsBase",
          function(x) unique(as.character(x$study)))
