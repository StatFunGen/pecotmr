# =============================================================================
# AnnotationMatrix S4 class
# -----------------------------------------------------------------------------
# Container for SNP-level annotations used in stratified heritability
# analysis. Supports binary (0/1) and continuous annotations classified as
# baseline (always jointly fitted) or candidate (score-tested).
# =============================================================================

#' @include allGenerics.R
NULL

#' @title Genomic Annotation Matrix
#' @description Container for SNP-level annotations used in stratified
#'   heritability analysis. Supports binary (0/1) and continuous annotations.
#'   Annotations are classified as baseline (always jointly fitted) or
#'   candidate (evaluated via score statistics).
#' @slot snpRanges A \code{GRanges} object with one range per SNP,
#'   defining genomic positions.
#' @slot annotations A numeric matrix (SNPs x annotations). Dense for
#'   small annotation counts, can be sparse (\code{dgCMatrix}) for large
#'   binary annotation sets.
#' @slot annotationMeta A \code{data.frame} with columns:
#'   \describe{
#'     \item{name}{Character, annotation name}
#'     \item{tier}{Character, one of "baseline" or "candidate"}
#'     \item{type}{Character, one of "binary" or "continuous"}
#'   }
#' @slot genome Character string for genome build.
#' @export
setClass("AnnotationMatrix",
  representation(
    snpRanges = "GRanges",
    annotations = "ANY",
    annotationMeta = "data.frame",
    genome = "character"
  ),
  validity = function(object) {
    errors <- character()
    n_snp <- length(object@snpRanges)
    n_annot <- ncol(object@annotations)
    if (nrow(object@annotations) != n_snp)
      errors <- c(errors,
        "Number of rows in 'annotations' must match length of 'snpRanges'")
    required_meta_cols <- c("name", "tier", "type")
    if (!all(required_meta_cols %in% colnames(object@annotationMeta)))
      errors <- c(errors,
        "annotationMeta must have columns: name, tier, type")
    if (nrow(object@annotationMeta) != n_annot)
      errors <- c(errors,
        "Number of rows in 'annotationMeta' must match annotation count")
    valid_tiers <- c("baseline", "candidate")
    if (!all(object@annotationMeta$tier %in% valid_tiers))
      errors <- c(errors,
        "annotationMeta$tier must be 'baseline' or 'candidate'")
    valid_types <- c("binary", "continuous")
    if (!all(object@annotationMeta$type %in% valid_types))
      errors <- c(errors,
        "annotationMeta$type must be 'binary' or 'continuous'")
    if (length(errors) == 0) TRUE else errors
  }
)

#' @export
setMethod("show", "AnnotationMatrix", function(object) {
  n_base <- sum(object@annotationMeta$tier == "baseline")
  n_cand <- sum(object@annotationMeta$tier == "candidate")
  n_bin <- sum(object@annotationMeta$type == "binary")
  n_cont <- sum(object@annotationMeta$type == "continuous")
  cat(sprintf("AnnotationMatrix: %d SNPs x %d annotations\n",
              nrow(object@annotations), ncol(object@annotations)))
  cat(sprintf("  Baseline: %d, Candidate: %d\n", n_base, n_cand))
  cat(sprintf("  Binary: %d, Continuous: %d\n", n_bin, n_cont))
  cat(sprintf("  Genome build: %s\n", object@genome))
})

#' @rdname getAnnotations
#' @export
setMethod("getAnnotations", "AnnotationMatrix", function(x) x@annotations)

#' @rdname getAnnotationMeta
#' @export
setMethod("getAnnotationMeta", "AnnotationMatrix",
          function(x) x@annotationMeta)

#' @rdname getSnpRanges
#' @export
setMethod("getSnpRanges", "AnnotationMatrix", function(x) x@snpRanges)

#' @rdname getGenome
#' @export
setMethod("getGenome", "AnnotationMatrix", function(x, ...) x@genome)
