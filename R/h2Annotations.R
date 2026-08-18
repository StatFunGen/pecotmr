#' @title Annotation Handling for Stratified Heritability
#' @description Read and manage genomic annotations for stratified heritability
#'   analysis. Supports BED, BigWig, and LDSC .annot formats.
#' @name pecotmr-h2-annotations
#' @keywords internal
#' @importFrom tools file_ext
#' @importFrom stringr str_detect str_to_lower
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @include AllGenerics.R
NULL

# =============================================================================
# Reader method
# =============================================================================

#' @rdname readAnnotations
#' @export
setMethod(
    "readAnnotations",
    signature(paths = "character"),
    function(paths, snpRanges, annotationMeta = NULL, genome = "hg19", ...) {
        if (is.null(names(paths))) {
            msg <- glue(
                "'paths' must be a named character vector (names = ",
                "annotation names)"
            )
            abort(msg)
        }
        annotNames <- names(paths)
        nAnnots <- length(paths)
        annotMat <- .readAnnotMatrix(
            paths,
            snpRanges,
            annotNames,
            length(snpRanges),
            nAnnots
        )
        if (is.null(annotationMeta)) {
            annotationMeta <- tibble(
                name = annotNames,
                tier = rep("candidate", nAnnots),
                type = .readAnnotTypes(paths)
            )
        }
        AnnotationMatrix(annotMat, snpRanges, annotationMeta, genome)
    }
)

# Auto-detected annotation types (continuous for BigWig, else binary). Unnamed
# so the tibble column carries no element names (tibble, unlike data.frame,
# would otherwise keep the paths' names on the column).
# @noRd
.readAnnotTypes <- function(paths) {
    unname(map_chr(paths, .annotType))
}

# @noRd
.annotType <- function(p) {
    if (.annotDetectFormat(p) == "bigwig") "continuous" else "binary"
}

# Build the (SNP x annotation) matrix by reading each annotation column.
# @noRd
.readAnnotMatrix <- function(paths, snpRanges, annotNames, nSnps, nAnnots) {
    annotMat <- matrix(0, nrow = nSnps, ncol = nAnnots)
    colnames(annotMat) <- annotNames
    for (i in seq_along(paths)) {
        annotMat[, i] <- .readAnnotColumn(
            paths[i],
            snpRanges,
            annotNames[i]
        )
    }
    annotMat
}

# One annotation column, dispatched by detected format (BigWig / .annot / BED).
# @noRd
.readAnnotColumn <- function(path, snpRanges, annotName) {
    fmt <- .annotDetectFormat(path)
    if (fmt == "bigwig") {
        .readBigwigAtSnps(path, snpRanges)
    } else if (fmt == "ldsc_annot") {
        .readLdscAnnot(path, snpRanges, annotName)
    } else {
        .readBedAnnotation(path, snpRanges)
    }
}

# =============================================================================
# Internal helpers
# =============================================================================

#' @title Detect Annotation File Format
#' @description Detect annotation file format from extension. This is separate
#'   from \code{.h2DetectFormat} because BED annotation files (genomic intervals
#'   for rtracklayer) must be distinguished from plink BED files.
#' @param path Character, file path.
#' @return Character, one of "bigwig", "ldsc_annot", or "bed".
#' @keywords internal
.annotDetectFormat <- function(path) {
    lpath <- str_to_lower(path)
    if (str_detect(lpath, "\\.annot\\.gz$")) {
        return("ldsc_annot")
    }

    ext <- str_to_lower(file_ext(path))
    switch(
        ext,
        "bw" = ,
        "bigwig" = "bigwig",
        "annot" = "ldsc_annot",
        # Default: treat as BED (genomic interval file for rtracklayer)
        "bed"
    )
}

#' @title Read BigWig Scores at SNP Positions
#' @description Import scores from a BigWig file at specified SNP positions.
#' @param bwPath Character, path to a BigWig file.
#' @param snpRanges A \code{GRanges} object with SNP positions.
#' @return Numeric vector of scores (length = number of SNPs).
#' @keywords internal
.readBigwigAtSnps <- function(bwPath, snpRanges) {
    bw <- rtracklayer::BigWigFile(bwPath)
    scores <- rtracklayer::import(bw, which = snpRanges, as = "NumericList")
    # Take mean score at each SNP position
    map_dbl(scores, .bwMeanScore)
}

#' @title Read BED Annotation
#' @description Read a BED file and compute binary overlap with SNP positions.
#' @param bedPath Character, path to a BED file.
#' @param snpRanges A \code{GRanges} object with SNP positions.
#' @return Numeric vector of 0/1 values (length = number of SNPs).
#' @keywords internal
.readBedAnnotation <- function(bedPath, snpRanges) {
    regions <- rtracklayer::import(bedPath)
    hits <- findOverlaps(snpRanges, regions)
    result <- rep(0L, length(snpRanges))
    result[queryHits(hits)] <- 1L
    as.numeric(result)
}

#' @title Read LDSC Annotation File
#' @description Read an S-LDSC .annot[.gz] file and extract a named annotation
#'   column, matched to SNP positions.
#' @param annotPath Character, path to an .annot or .annot.gz file.
#' @param snpRanges A \code{GRanges} object with SNP positions.
#' @param annotName Character, name of the annotation column to extract.
#' @return Numeric vector of annotation values (length = number of SNPs).
#' @keywords internal
.readLdscAnnot <- function(annotPath, snpRanges, annotName) {
    # S-LDSC .annot files are tab-separated with columns: CHR, BP, SNP, CM, ...
    dt <- vroom(annotPath, show_col_types = FALSE)

    if (!is_in(annotName, colnames(dt))) {
        msg <- glue("Annotation column '{annotName}' not found in {annotPath}")
        abort(msg)
    }

    if (!all(is_in(c("CHR", "BP"), colnames(dt)))) {
        abort("LDSC annot file must contain CHR and BP columns")
    }

    # Build GRanges from the annot file positions
    annotGr <- GRanges(
        seqnames = withChrPrefix(dt$CHR),
        ranges = IRanges(start = dt$BP, width = 1L)
    )

    # Match SNPs by genomic position
    hits <- findOverlaps(snpRanges, annotGr)

    # Initialize result with default 0
    result <- rep(0, length(snpRanges))
    result[queryHits(hits)] <-
        as.numeric(dt[[annotName]][subjectHits(hits)])

    result
}

# (The AnnotationMatrix() constructor and the getBaseline / getCandidates tier
# accessors now live in R/AnnotationMatrix.R alongside the class definition.)

# The mean BigWig score at one SNP (0 when the SNP has no overlapping
# intervals).
# @noRd
.bwMeanScore <- function(x) {
    if (length(x) > 0) mean(x) else 0
}
