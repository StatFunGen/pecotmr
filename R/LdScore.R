# =============================================================================
# LdScore S4 class
# -----------------------------------------------------------------------------
# Pre-computed LD scores (sum of r^2) per SNP. Consumed by S-LDSC and
# g-LDSC. Holds the optional per-block LD matrices needed for g-LDSC's
# FGLS residual covariance.
# =============================================================================

#' @include LdStatistic.R
NULL

#' @title LD Score-Based LD Statistic
#' @description Pre-computed LD scores for each SNP. Used by S-LDSC and
#'   g-LDSC. Supports both standard LD scores and annotation-stratified
#'   LD scores.
#' @slot ldScores A numeric matrix (SNPs x annotations+1). The first
#'   column is the base LD score (sum of r^2). Additional columns are
#'   annotation-stratified LD scores if annotations are provided.
#' @slot ldScoreWeights A numeric vector of regression weights for each SNP.
#' @slot ldMatrixList For g-LDSC: a list of per-block LD (R^2) matrices
#'   used to compute the FGLS residual covariance. NULL for S-LDSC.
#' @export
setClass("LdScore",
  contains = "LdStatistic",
  representation(
    ldScores = "matrix",
    ldScoreWeights = "numeric",
    ldMatrixList = "list"
  ),
  validity = function(object) {
    parent_check <- getValidity(getClass("LdStatistic"))(object)
    errors <- if (isTRUE(parent_check)) character() else parent_check
    if (nrow(object@ldScores) != nrow(object@snpInfo))
      errors <- c(errors,
        "Number of rows in 'ldScores' must match 'snpInfo'")
    if (length(object@ldScoreWeights) != nrow(object@snpInfo))
      errors <- c(errors,
        "Length of 'ldScoreWeights' must match 'snpInfo'")
    if (length(errors) == 0) TRUE else errors
  }
)

#' @export
setMethod("show", "LdScore", function(object) {
  n_scores <- ncol(object@ldScores)
  has_matrix <- length(object@ldMatrixList) > 0
  cat(sprintf("LdScore: %d SNPs, %d LD score columns\n",
              nrow(object@snpInfo), n_scores))
  cat(sprintf("  Full LD matrices: %s (needed for g-LDSC)\n", has_matrix))
  cat(sprintf("  Reference N: %d, In-sample: %s\n",
              object@nRef, object@inSample))
})

#' @rdname getLdScores
#' @export
setMethod("getLdScores", "LdScore", function(x) x@ldScores)

#' @rdname getLdScoreWeights
#' @export
setMethod("getLdScoreWeights", "LdScore", function(x) x@ldScoreWeights)

#' @rdname getLdMatrixList
#' @export
setMethod("getLdMatrixList", "LdScore", function(x) x@ldMatrixList)
