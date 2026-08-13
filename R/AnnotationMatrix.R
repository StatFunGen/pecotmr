# =============================================================================
# AnnotationMatrix S4 class
# -----------------------------------------------------------------------------
# Container for SNP-level annotations used in stratified heritability
# analysis. Supports binary (0/1) and continuous annotations classified as
# baseline (always jointly fitted) or candidate (score-tested).
# =============================================================================

#' @include AllGenerics.R
NULL

#' @title Genomic Annotation Matrix
#' @description Container for SNP-level annotations used in stratified
#'   heritability analysis. Supports binary (0/1) and continuous annotations.
#'   Annotations are classified as baseline (always jointly fitted) or candidate
#'   (evaluated via score statistics).
#' @slot snpRanges A \code{GRanges} object with one range per SNP, defining
#'   genomic positions.
#' @slot annotations A numeric matrix (SNPs x annotations). Dense for small
#'   annotation counts, can be sparse (\code{dgCMatrix}) for large binary
#'   annotation sets.
#' @slot annotationMeta A \code{data.frame} with columns:
#'   \describe{
#'     \item{name}{Character, annotation name}
#'     \item{tier}{Character, one of "baseline" or "candidate"}
#'     \item{type}{Character, one of "binary" or "continuous"}
#'   }
#' @slot genome Character string for genome build.
#' @export
setClass(
    "AnnotationMatrix",
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
        if (nrow(object@annotations) != n_snp) {
            errors <- c(
                errors,
                paste0(
                    "Number of rows in 'annotations' must match ",
                    "length of 'snpRanges'"
                )
            )
        }
        required_meta_cols <- c("name", "tier", "type")
        if (!all(required_meta_cols %in% colnames(object@annotationMeta))) {
            errors <- c(
                errors,
                "annotationMeta must have columns: name, tier, type"
            )
        }
        if (nrow(object@annotationMeta) != n_annot) {
            errors <- c(
                errors,
                "Number of rows in 'annotationMeta' must match annotation count"
            )
        }
        valid_tiers <- c("baseline", "candidate")
        if (!all(object@annotationMeta$tier %in% valid_tiers)) {
            errors <- c(
                errors,
                "annotationMeta$tier must be 'baseline' or 'candidate'"
            )
        }
        valid_types <- c("binary", "continuous")
        if (!all(object@annotationMeta$type %in% valid_types)) {
            errors <- c(
                errors,
                "annotationMeta$type must be 'binary' or 'continuous'"
            )
        }
        if (length(errors) == 0) TRUE else errors
    }
)

#' @rdname show-methods
#' @importFrom methods show
#' @export
setMethod("show", "AnnotationMatrix", function(object) {
    n_base <- sum(object@annotationMeta$tier == "baseline")
    n_cand <- sum(object@annotationMeta$tier == "candidate")
    n_bin <- sum(object@annotationMeta$type == "binary")
    n_cont <- sum(object@annotationMeta$type == "continuous")
    cat(sprintf(
        "AnnotationMatrix: %d SNPs x %d annotations\n",
        nrow(object@annotations),
        ncol(object@annotations)
    ))
    cat(sprintf("  Baseline: %d, Candidate: %d\n", n_base, n_cand))
    cat(sprintf("  Binary: %d, Continuous: %d\n", n_bin, n_cont))
    cat(sprintf("  Genome build: %s\n", object@genome))
})

#' @rdname getAnnotations
#' @export
setMethod("getAnnotations", "AnnotationMatrix", function(x) x@annotations)

#' @rdname getAnnotationMeta
#' @export
setMethod("getAnnotationMeta", "AnnotationMatrix", function(x) x@annotationMeta)

#' @rdname getSnpRanges
#' @export
setMethod("getSnpRanges", "AnnotationMatrix", function(x) x@snpRanges)

#' @rdname getGenome
#' @export
setMethod("getGenome", "AnnotationMatrix", function(x, ...) x@genome)

# =============================================================================
# Constructor
# =============================================================================

#' @title Create an AnnotationMatrix Object
#' @description Construct an \code{AnnotationMatrix} from a matrix and metadata.
#' @param annotations A numeric matrix or sparse matrix (SNPs x annotations).
#' @param snpRanges A \code{GRanges} object with SNP positions.
#' @param annotationMeta A data.frame with columns: name, tier, type.
#' @param genome Character, genome build.
#' @return An \code{AnnotationMatrix} object.
#' @examples
#' snpRanges <- GenomicRanges::GRanges(
#'   "22", IRanges::IRanges((1:10) * 100, width = 1))
#' annotations <- matrix(rbinom(50, 1, 0.3), 10, 5,
#'   dimnames = list(NULL, paste0("annot", 1:5)))
#' meta <- data.frame(name = paste0("annot", 1:5), tier = "baseline",
#'   type = "binary")
#' am <- AnnotationMatrix(annotations, snpRanges, annotationMeta = meta)
#' am
#' @export
AnnotationMatrix <- function(
    annotations,
    snpRanges,
    annotationMeta,
    genome = "hg19"
) {
    # Validate annotationMeta
    if (!is.data.frame(annotationMeta)) {
        stop("annotationMeta must be a data.frame")
    }

    requiredCols <- c("name", "tier", "type")
    if (!all(requiredCols %in% colnames(annotationMeta))) {
        stop("annotationMeta must have columns: name, tier, type")
    }

    # Set column names on matrix
    if (is.null(colnames(annotations))) {
        colnames(annotations) <- annotationMeta$name
    }

    new(
        "AnnotationMatrix",
        snpRanges = snpRanges,
        annotations = annotations,
        annotationMeta = annotationMeta,
        genome = genome
    )
}

# =============================================================================
# Tier accessors
# =============================================================================

#' @title Get Baseline Annotations
#' @description Extract only baseline-tier annotations from an
#'   \code{AnnotationMatrix}.
#' @param annot An \code{AnnotationMatrix} object.
#' @return An \code{AnnotationMatrix} with only baseline annotations.
#' @examples
#' snpRanges <- GenomicRanges::GRanges(
#'   "22", IRanges::IRanges((1:10) * 100, width = 1))
#' annotations <- matrix(rbinom(50, 1, 0.3), 10, 5,
#'   dimnames = list(NULL, paste0("annot", 1:5)))
#' meta <- data.frame(name = paste0("annot", 1:5), tier = "baseline",
#'   type = "binary")
#' am <- AnnotationMatrix(annotations, snpRanges, annotationMeta = meta)
#' getBaseline(am)
#' @export
getBaseline <- function(annot) {
    meta <- getAnnotationMeta(annot)
    idx <- meta$tier == "baseline"
    AnnotationMatrix(
        annotations = getAnnotations(annot)[, idx, drop = FALSE],
        snpRanges = getSnpRanges(annot),
        annotationMeta = meta[idx, , drop = FALSE],
        genome = getGenome(annot)
    )
}

#' @title Get Candidate Annotations
#' @description Extract only candidate-tier annotations from an
#'   \code{AnnotationMatrix}.
#' @param annot An \code{AnnotationMatrix} object.
#' @return An \code{AnnotationMatrix} with only candidate annotations.
#' @examples
#' snpRanges <- GenomicRanges::GRanges(
#'   "22", IRanges::IRanges((1:10) * 100, width = 1))
#' annotations <- matrix(rbinom(50, 1, 0.3), 10, 5,
#'   dimnames = list(NULL, paste0("annot", 1:5)))
#' meta <- data.frame(name = paste0("annot", 1:5), tier = "baseline",
#'   type = "binary")
#' am <- AnnotationMatrix(annotations, snpRanges, annotationMeta = meta)
#' getCandidates(am)
#' @export
getCandidates <- function(annot) {
    meta <- getAnnotationMeta(annot)
    idx <- meta$tier == "candidate"
    AnnotationMatrix(
        annotations = getAnnotations(annot)[, idx, drop = FALSE],
        snpRanges = getSnpRanges(annot),
        annotationMeta = meta[idx, , drop = FALSE],
        genome = getGenome(annot)
    )
}
