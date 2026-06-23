# =============================================================================
# AllClasses.R
# -----------------------------------------------------------------------------
# Virtual base classes shared across the package. Concrete subclasses live
# in their own per-class files (QtlSumStats.R, GwasSumStats.R, QtlDataset.R,
# QtlFineMappingResult.R, GwasFineMappingResult.R, etc.).
#
# Per Bioconductor convention this file is loaded first in the Collate
# ordering (the "AllClasses.R" filename sorts to the top of the alphabet),
# and every method-bearing file uses `@include AllClasses.R` so roxygen
# topologically orders the Collate field for us.
# =============================================================================

#' @include AllGenerics.R GenotypeHandle.R
#' @importFrom methods setClass setMethod new is validObject
NULL

# =============================================================================
# SumStatsBase
# -----------------------------------------------------------------------------
# Shared parent of the QTL and GWAS summary statistics collections.
# Concrete subclasses (QtlSumStats, GwasSumStats) inherit from DFrame and
# share the ldSketch / genome / qcInfo slots. Class-specific accessors
# (getZ / getN / getMaf / nSnps / subsetChr / getVarY / getSumStats)
# stay on the concrete subclass because they rely on the tuple shape
# (3-tuple for QtlSumStats, 1-tuple for GwasSumStats).
# =============================================================================

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

#' @rdname getQcDiagnostics
#' @export
setMethod("getQcDiagnostics", "SumStatsBase",
  function(x, entry = 1L, ...) {
    qc <- x@qcInfo
    if (length(qc) == 0L) return(NULL)
    audits <- qc$entryAudit
    if (is.null(audits)) return(NULL)
    if (is.null(entry)) {
      out <- lapply(audits, function(a) a$ldMismatchDiagnostics)
      keep <- !vapply(out, is.null, logical(1L))
      if (!any(keep)) return(NULL)
      setNames(out[keep], seq_along(audits)[keep])
    } else {
      if (!is.numeric(entry) || length(entry) != 1L ||
          entry < 1L || entry > length(audits)) {
        stop("`entry` must be a single integer in 1:", length(audits), ".")
      }
      audits[[as.integer(entry)]]$ldMismatchDiagnostics
    }
  })

#' @rdname getLdSketch
#' @export
setMethod("getLdSketch", "SumStatsBase", function(x, ...) x@ldSketch)

#' @rdname getStudy
#' @export
setMethod("getStudy", "SumStatsBase",
          function(x) unique(as.character(x$study)))

# =============================================================================
# FineMappingResultBase
# -----------------------------------------------------------------------------
# Shared parent of the QTL and GWAS fine-mapping result collections.
# Concrete subclasses (QtlFineMappingResult, GwasFineMappingResult) carry
# a DFrame of per-fit rows plus a shared ldSketch slot. Downstream
# pipelines dispatch on FineMappingResultBase for behaviors that apply to
# either flavour, and on the concrete subclass when the tuple shape
# matters.
# =============================================================================

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
