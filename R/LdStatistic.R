# =============================================================================
# LdStatistic S4 virtual class
# -----------------------------------------------------------------------------
# Abstract container for pre-computed LD statistics. Subclasses
# (LdEigen, LdScore) provide method-specific representations: LdEigen
# for eigendecomposition-based methods (LDER/HDL/sHDL), LdScore for
# LD-score-based methods (S-LDSC/g-LDSC).
# =============================================================================

#' @include AllGenerics.R LdBlocks.R
NULL

#' @title LD Statistic (Virtual Base Class)
#' @description Abstract container for pre-computed LD statistics. Subclasses
#'   provide method-specific representations: eigendecompositions (for
#'   LDER/HDL/sHDL) and LD score matrices (for S-LDSC/g-LDSC).
#' @slot ldBlocks An \code{LdBlocks} object defining the block structure.
#' @slot snpInfo A \code{data.frame} with columns \code{SNP}, \code{CHR},
#'   \code{BP}, \code{A1}, \code{A2}, and optionally \code{MAF}.
#' @slot nRef Integer, sample size of the LD reference panel.
#' @slot inSample Logical, whether the LD reference is from the same
#'   cohort as the GWAS (affects bias correction).
#' @slot genome Character string for genome build.
#' @export
setClass("LdStatistic",
  contains = "VIRTUAL",
  representation(
    ldBlocks = "LdBlocks",
    snpInfo = "data.frame",
    nRef = "integer",
    inSample = "logical",
    genome = "character"
  ),
  validity = function(object) {
    errors <- character()
    if (length(object@nRef) != 1L || object@nRef <= 0L)
      errors <- c(errors, "'nRef' must be a single positive integer")
    if (length(object@inSample) != 1L)
      errors <- c(errors, "'inSample' must be a single logical value")
    if (length(object@genome) != 1L || !nzchar(object@genome))
      errors <- c(errors, "'genome' must be a single non-empty character string")
    if (nrow(object@snpInfo) == 0L)
      errors <- c(errors, "'snpInfo' must have at least one row")
    if (length(errors) == 0) TRUE else errors
  }
)

#' @rdname getSnpInfo
#' @export
setMethod("getSnpInfo", "LdStatistic", function(x) x@snpInfo)

#' @rdname getNRef
#' @export
setMethod("getNRef", "LdStatistic", function(x) x@nRef)

#' @rdname getInSample
#' @export
setMethod("getInSample", "LdStatistic", function(x) x@inSample)

#' @rdname getLdBlocks
#' @export
setMethod("getLdBlocks", "LdStatistic", function(x) x@ldBlocks)

#' @rdname getGenome
#' @export
setMethod("getGenome", "LdStatistic", function(x, ...) x@genome)
