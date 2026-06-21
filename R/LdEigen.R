# =============================================================================
# LdEigen S4 class
# -----------------------------------------------------------------------------
# Pre-computed per-block eigendecompositions of the LD correlation matrix.
# Consumed by LDER / HDL / sHDL h2 estimators.
# =============================================================================

#' @include LdStatistic.R
NULL

#' @title Eigendecomposition-Based LD Statistic
#' @description Pre-computed per-block eigendecompositions of the LD
#'   correlation matrix. Used by LDER, HDL, and sHDL.
#' @slot eigenList A list of length \code{nBlocks}, each element a list
#'   with components:
#'   \describe{
#'     \item{values}{Numeric vector of eigenvalues}
#'     \item{vectors}{Numeric matrix of eigenvectors (SNPs x retained components)}
#'     \item{snpIdx}{Integer vector of SNP indices in \code{snpInfo}}
#'   }
#' @slot eigenvalueTruncation Numeric, proportion of variance retained
#'   (e.g., 0.9 for HDL's default). If 1.0, no truncation.
#' @export
setClass("LdEigen",
  contains = "LdStatistic",
  representation(
    eigenList = "list",
    eigenvalueTruncation = "numeric"
  ),
  validity = function(object) {
    parent_check <- getValidity(getClass("LdStatistic"))(object)
    errors <- if (isTRUE(parent_check)) character() else parent_check
    nBlocks <- length(object@ldBlocks@blocks)
    if (length(object@eigenList) != nBlocks)
      errors <- c(errors,
        "Length of 'eigenList' must match number of LD blocks")
    if (length(object@eigenvalueTruncation) != 1L ||
        object@eigenvalueTruncation <= 0 ||
        object@eigenvalueTruncation > 1)
      errors <- c(errors,
        "'eigenvalueTruncation' must be a single value in (0, 1]")
    if (length(errors) == 0) TRUE else errors
  }
)

#' @export
setMethod("show", "LdEigen", function(object) {
  cat(sprintf("LdEigen: %d SNPs across %d blocks\n",
              nrow(object@snpInfo), length(object@eigenList)))
  cat(sprintf("  Eigenvalue truncation: %.2f\n",
              object@eigenvalueTruncation))
  cat(sprintf("  Reference N: %d, In-sample: %s\n",
              object@nRef, object@inSample))
})

#' @rdname getEigenList
#' @export
setMethod("getEigenList", "LdEigen", function(x) x@eigenList)
