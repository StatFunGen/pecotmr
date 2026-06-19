#' @title Annotation Handling for Stratified Heritability
#' @description Read and manage genomic annotations for stratified
#'   heritability analysis. Supports BED, BigWig, and LDSC .annot formats.
#' @name pecotmr-h2-annotations
#' @keywords internal
#' @importFrom tools file_ext
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @include allGenerics.R
NULL

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
#' @export
AnnotationMatrix <- function(annotations, snpRanges, annotationMeta,
                             genome = "hg19") {
  # Validate annotationMeta
  if (!is.data.frame(annotationMeta))
    stop("annotationMeta must be a data.frame")

  requiredCols <- c("name", "tier", "type")
  if (!all(requiredCols %in% colnames(annotationMeta)))
    stop("annotationMeta must have columns: name, tier, type")

  # Set column names on matrix
  if (is.null(colnames(annotations)))
    colnames(annotations) <- annotationMeta$name

  new("AnnotationMatrix",
    snpRanges = snpRanges,
    annotations = annotations,
    annotationMeta = annotationMeta,
    genome = genome
  )
}

# =============================================================================
# Reader method
# =============================================================================

#' @rdname readAnnotations
#' @export
setMethod("readAnnotations",
  signature(paths = "character"),
  function(paths, snpRanges, annotationMeta = NULL, genome = "hg19", ...) {

    if (is.null(names(paths)))
      stop("'paths' must be a named character vector (names = annotation names)")

    annotNames <- names(paths)
    nSnps <- length(snpRanges)
    nAnnots <- length(paths)

    # Auto-detect types from file extensions
    types <- vapply(paths, function(p) {
      fmt <- .annotDetectFormat(p)
      if (fmt == "bigwig") "continuous"
      else "binary"
    }, character(1))

    # Initialize annotation matrix
    annotMat <- matrix(0, nrow = nSnps, ncol = nAnnots)
    colnames(annotMat) <- annotNames

    for (i in seq_along(paths)) {
      fmt <- .annotDetectFormat(paths[i])

      if (fmt == "bigwig") {
        # Continuous annotation from BigWig
        annotMat[, i] <- .readBigwigAtSnps(paths[i], snpRanges)
      } else if (fmt == "ldsc_annot") {
        # S-LDSC .annot format
        annotMat[, i] <- .readLdscAnnot(paths[i], snpRanges,
                                        annotNames[i])
      } else {
        # Binary annotation from BED or similar
        annotMat[, i] <- .readBedAnnotation(paths[i], snpRanges)
      }
    }

    # Build annotationMeta if not provided
    if (is.null(annotationMeta)) {
      annotationMeta <- data.frame(
        name = annotNames,
        tier = rep("candidate", nAnnots),
        type = types,
        stringsAsFactors = FALSE
      )
    }

    AnnotationMatrix(annotMat, snpRanges, annotationMeta, genome)
  }
)

# =============================================================================
# Internal helpers
# =============================================================================

#' @title Detect Annotation File Format
#' @description Detect annotation file format from extension. This is separate
#'   from \code{.h2DetectFormat} because BED annotation files (genomic
#'   intervals for rtracklayer) must be distinguished from plink BED files.
#' @param path Character, file path.
#' @return Character, one of "bigwig", "ldsc_annot", or "bed".
#' @keywords internal
.annotDetectFormat <- function(path) {
  lpath <- tolower(path)
  if (grepl("\\.annot\\.gz$", lpath))
    return("ldsc_annot")

  ext <- tolower(file_ext(path))
  switch(ext,
    "bw" = , "bigwig" = "bigwig",
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
  vapply(scores, function(x) if (length(x) > 0) mean(x) else 0,
         numeric(1))
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
  dt <- as.data.frame(vroom(annotPath, show_col_types = FALSE))

  if (!annotName %in% colnames(dt))
    stop("Annotation column '", annotName, "' not found in ", annotPath)

  if (!all(c("CHR", "BP") %in% colnames(dt)))
    stop("LDSC annot file must contain CHR and BP columns")

  # Build GRanges from the annot file positions
  annotGr <- GRanges(
    seqnames = paste0("chr", sub("^chr", "", dt$CHR)),
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

# =============================================================================
# Annotation subsetting
# =============================================================================

#' @title Get Baseline Annotations
#' @description Extract only baseline-tier annotations from an
#'   \code{AnnotationMatrix}.
#' @param annot An \code{AnnotationMatrix} object.
#' @return An \code{AnnotationMatrix} with only baseline annotations.
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
