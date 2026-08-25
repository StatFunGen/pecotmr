# =============================================================================
# AnnotationMatrix S4 class
# -----------------------------------------------------------------------------
# Container for SNP-level annotations used in stratified heritability
# analysis. Supports binary (0/1) and continuous annotations classified as
# baseline (always jointly fitted) or candidate (score-tested).
#
# A RangedSummarizedExperiment with nothing added: the SNP positions are
# rowRanges, the annotation matrix is an assay, and the per-annotation table
# (name / tier / type) is colData -- which is what it always was, written out
# by hand. Tier selection is therefore ordinary column subsetting, and the
# rows stay aligned to it without anything rebuilding the object.
# =============================================================================

#' @include AllGenerics.R
NULL

#' @title Genomic Annotation Matrix
#' @description Container for SNP-level annotations used in stratified
#'   heritability analysis. Supports binary (0/1) and continuous annotations,
#'   classified as baseline (always jointly fitted) or candidate (evaluated
#'   via score statistics).
#'
#'   It \strong{is} a \code{RangedSummarizedExperiment}: SNP positions in
#'   \code{rowRanges}, the annotation matrix in the \code{annotations}
#'   assay, and the per-annotation table in \code{colData}. Subsetting works
#'   in both directions -- \code{x[i, j]} narrows SNPs and annotations
#'   together -- and \code{\link{getBaseline}} / \code{\link{getCandidates}}
#'   are column subsets rather than reconstructions.
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @export
setClass(
    "AnnotationMatrix",
    contains = "RangedSummarizedExperiment",
    validity = function(object) .validateAnnotationMatrix(object)
)

# The tier/type vocabulary is the only thing left to check: the SNP-by-
# annotation shape is now enforced by SummarizedExperiment itself.
# @noRd
.validateAnnotationMatrix <- function(object) {
    errors <- character()
    cd <- SummarizedExperiment::colData(object)
    required <- c("name", "tier", "type")
    if (!all(is_in(required, colnames(cd)))) {
        return("annotationMeta must have columns: name, tier, type")
    }
    if (!all(is_in(cd$tier, c("baseline", "candidate")))) {
        errors <- c(
            errors,
            "annotationMeta$tier must be 'baseline' or 'candidate'"
        )
    }
    if (!all(is_in(cd$type, c("binary", "continuous")))) {
        errors <- c(
            errors,
            "annotationMeta$type must be 'binary' or 'continuous'"
        )
    }
    if (length(errors) == 0) TRUE else errors
}

#' @rdname show-methods
#' @importFrom methods show
#' @export
setMethod("show", "AnnotationMatrix", function(object) {
    meta <- SummarizedExperiment::colData(object)
    cat(glue(
        "AnnotationMatrix: {nrow(object)} SNPs x ",
        "{ncol(object)} annotations\n",
        .trim = FALSE
    ))
    cat(glue(
        "  Baseline: {sum(meta$tier == 'baseline')}, ",
        "Candidate: {sum(meta$tier == 'candidate')}\n",
        .trim = FALSE
    ))
    cat(glue(
        "  Binary: {sum(meta$type == 'binary')}, ",
        "Continuous: {sum(meta$type == 'continuous')}\n",
        .trim = FALSE
    ))
    cat(glue("  Genome build: {getGenome(object)}\n", .trim = FALSE))
})

#' @rdname getGenome
#' @export
setMethod("getGenome", "AnnotationMatrix", function(x, ...) {
    build <- unique(GenomeInfoDb::genome(SummarizedExperiment::rowRanges(x)))
    build <- build[!is.na(build)]
    if (length(build) == 0L) NA_character_ else build[[1L]]
})

# =============================================================================
# Constructor
# =============================================================================

# Shape checks that do not involve snpRanges.
#
# Checked here rather than left to SummarizedExperiment, whose message for a
# mismatch is "'x@assays' is not parallel to 'x'" -- accurate but no help in
# finding which of the three inputs is the odd one.
#
# The column count is checked BEFORE the caller defaults colnames from
# annotationMeta$name: naming an n-column matrix from an m-row table fails
# inside `colnames<-` with "length of 'dimnames' [2] not equal to array
# extent", which is the same unhelpful shape this guard exists to replace.
# @noRd
.amCheckInputs <- function(annotations, annotationMeta) {
    if (!is.data.frame(annotationMeta)) {
        abort("annotationMeta must be a data.frame")
    }
    requiredCols <- c("name", "tier", "type")
    if (!all(is_in(requiredCols, colnames(annotationMeta)))) {
        abort("annotationMeta must have columns: name, tier, type")
    }
    if (ncol(annotations) != nrow(annotationMeta)) {
        abort(glue(
            "`annotations` has {ncol(annotations)} column(s) for ",
            "{nrow(annotationMeta)} metadata row(s); they must match."
        ))
    }
    invisible(NULL)
}


#' @title Create an AnnotationMatrix Object
#' @description Construct an \code{AnnotationMatrix} from a matrix and
#'   metadata.
#' @param annotations A numeric matrix or sparse matrix (SNPs x annotations).
#' @param snpRanges A \code{GRanges} object with SNP positions.
#' @param annotationMeta A data.frame with columns: name, tier, type.
#' @param genome Character, genome build; recorded in \code{seqinfo()}.
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
#' dim(am)
#' @export
AnnotationMatrix <- function(
    annotations,
    snpRanges,
    annotationMeta,
    genome = "hg19"
) {
    .amCheckInputs(annotations, annotationMeta)
    if (is.null(colnames(annotations))) {
        colnames(annotations) <- annotationMeta$name
    }
    if (nrow(annotations) != length(snpRanges)) {
        abort(glue(
            "`annotations` has {nrow(annotations)} row(s) for ",
            "{length(snpRanges)} SNP range(s); the rows must match."
        ))
    }
    if (!is.null(genome) && length(genome) == 1L && !is.na(genome)) {
        GenomeInfoDb::genome(snpRanges) <- genome
    }
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(annotations = annotations),
        rowRanges = snpRanges,
        # Keyed off the assay's own colnames, not annotationMeta$name:
        # SummarizedExperiment requires the two to agree, and the previous
        # class let them diverge, so taking the name column would reject
        # objects that were legal before.
        colData = S4Vectors::DataFrame(
            annotationMeta,
            row.names = colnames(annotations)
        )
    )
    obj <- methods::new("AnnotationMatrix", se)
    validObject(obj)
    obj
}

# =============================================================================
# Tier accessors
# =============================================================================

# Annotations of one tier. A column subset: SummarizedExperiment keeps the
# rows, the assay and the per-annotation table aligned, where the previous
# implementation rebuilt the object from three separately-subset pieces.
# @noRd
.annotTier <- function(annot, tier) {
    annot[, SummarizedExperiment::colData(annot)$tier == tier]
}

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
    .annotTier(annot, "baseline")
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
#' meta <- data.frame(name = paste0("annot", 1:5),
#'   tier = c(rep("baseline", 3), rep("candidate", 2)), type = "binary")
#' am <- AnnotationMatrix(annotations, snpRanges, annotationMeta = meta)
#' getCandidates(am)
#' @export
getCandidates <- function(annot) {
    .annotTier(annot, "candidate")
}
