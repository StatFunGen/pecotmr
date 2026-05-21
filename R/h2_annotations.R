#' @title Annotation Handling for Stratified Heritability
#' @description Read and manage genomic annotations for stratified
#'   heritability analysis. Supports BED, BigWig, and LDSC .annot formats.
#' @importFrom tools file_ext
#' @importFrom GenomicRanges GRanges
#' @include AllGenerics.R
NULL

# =============================================================================
# Constructor
# =============================================================================

#' @title Create an AnnotationMatrix Object
#' @description Construct an \code{AnnotationMatrix} from a matrix and metadata.
#' @param annotations A numeric matrix or sparse matrix (SNPs x annotations).
#' @param snp_ranges A \code{GRanges} object with SNP positions.
#' @param annotation_meta A data.frame with columns: name, tier, type.
#' @param genome Character, genome build.
#' @return An \code{AnnotationMatrix} object.
#' @export
AnnotationMatrix <- function(annotations, snp_ranges, annotation_meta,
                             genome = "hg19") {
  # Validate annotation_meta
  if (!is.data.frame(annotation_meta))
    stop("annotation_meta must be a data.frame")

  required_cols <- c("name", "tier", "type")
  if (!all(required_cols %in% colnames(annotation_meta)))
    stop("annotation_meta must have columns: name, tier, type")

  # Set column names on matrix
  if (is.null(colnames(annotations)))
    colnames(annotations) <- annotation_meta$name

  new("AnnotationMatrix",
    snp_ranges = snp_ranges,
    annotations = annotations,
    annotation_meta = annotation_meta,
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
  function(paths, snp_ranges, annotation_meta = NULL, genome = "hg19", ...) {

    if (is.null(names(paths)))
      stop("'paths' must be a named character vector (names = annotation names)")

    annot_names <- names(paths)
    n_snps <- length(snp_ranges)
    n_annots <- length(paths)

    # Auto-detect types from file extensions
    types <- vapply(paths, function(p) {
      fmt <- .annot_detect_format(p)
      if (fmt == "bigwig") "continuous"
      else "binary"
    }, character(1))

    # Initialize annotation matrix
    annot_mat <- matrix(0, nrow = n_snps, ncol = n_annots)
    colnames(annot_mat) <- annot_names

    for (i in seq_along(paths)) {
      fmt <- .annot_detect_format(paths[i])

      if (fmt == "bigwig") {
        # Continuous annotation from BigWig
        annot_mat[, i] <- .read_bigwig_at_snps(paths[i], snp_ranges)
      } else if (fmt == "ldsc_annot") {
        # S-LDSC .annot format
        annot_mat[, i] <- .read_ldsc_annot(paths[i], snp_ranges,
                                            annot_names[i])
      } else {
        # Binary annotation from BED or similar
        annot_mat[, i] <- .read_bed_annotation(paths[i], snp_ranges)
      }
    }

    # Build annotation_meta if not provided
    if (is.null(annotation_meta)) {
      annotation_meta <- data.frame(
        name = annot_names,
        tier = rep("candidate", n_annots),
        type = types,
        stringsAsFactors = FALSE
      )
    }

    AnnotationMatrix(annot_mat, snp_ranges, annotation_meta, genome)
  }
)

# =============================================================================
# Internal helpers
# =============================================================================

#' @title Detect Annotation File Format
#' @description Detect annotation file format from extension. This is separate
#'   from \code{.h2_detect_format} because BED annotation files (genomic
#'   intervals for rtracklayer) must be distinguished from plink BED files.
#' @param path Character, file path.
#' @return Character, one of "bigwig", "ldsc_annot", or "bed".
#' @keywords internal
.annot_detect_format <- function(path) {
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
#' @param bw_path Character, path to a BigWig file.
#' @param snp_ranges A \code{GRanges} object with SNP positions.
#' @return Numeric vector of scores (length = number of SNPs).
#' @keywords internal
.read_bigwig_at_snps <- function(bw_path, snp_ranges) {
  bw <- rtracklayer::BigWigFile(bw_path)
  scores <- rtracklayer::import(bw, which = snp_ranges, as = "NumericList")
  # Take mean score at each SNP position
  vapply(scores, function(x) if (length(x) > 0) mean(x) else 0,
         numeric(1))
}

#' @title Read BED Annotation
#' @description Read a BED file and compute binary overlap with SNP positions.
#' @param bed_path Character, path to a BED file.
#' @param snp_ranges A \code{GRanges} object with SNP positions.
#' @return Numeric vector of 0/1 values (length = number of SNPs).
#' @keywords internal
.read_bed_annotation <- function(bed_path, snp_ranges) {
  regions <- rtracklayer::import(bed_path)
  hits <- findOverlaps(snp_ranges, regions)
  result <- rep(0L, length(snp_ranges))
  result[queryHits(hits)] <- 1L
  as.numeric(result)
}

#' @title Read LDSC Annotation File
#' @description Read an S-LDSC .annot[.gz] file and extract a named annotation
#'   column, matched to SNP positions.
#' @param annot_path Character, path to an .annot or .annot.gz file.
#' @param snp_ranges A \code{GRanges} object with SNP positions.
#' @param annot_name Character, name of the annotation column to extract.
#' @return Numeric vector of annotation values (length = number of SNPs).
#' @keywords internal
.read_ldsc_annot <- function(annot_path, snp_ranges, annot_name) {
  # S-LDSC .annot files are tab-separated with columns: CHR, BP, SNP, CM, ...
  dt <- as.data.frame(vroom(annot_path, show_col_types = FALSE))

  if (!annot_name %in% colnames(dt))
    stop("Annotation column '", annot_name, "' not found in ", annot_path)

  if (!all(c("CHR", "BP") %in% colnames(dt)))
    stop("LDSC annot file must contain CHR and BP columns")

  # Build GRanges from the annot file positions
  annot_gr <- GRanges(
    seqnames = paste0("chr", sub("^chr", "", dt$CHR)),
    ranges = IRanges(start = dt$BP, width = 1L)
  )

  # Match SNPs by genomic position
  hits <- findOverlaps(snp_ranges, annot_gr)

  # Initialize result with default 0
  result <- rep(0, length(snp_ranges))
  result[queryHits(hits)] <-
    as.numeric(dt[[annot_name]][subjectHits(hits)])

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
  idx <- annot@annotation_meta$tier == "baseline"
  AnnotationMatrix(
    annotations = annot@annotations[, idx, drop = FALSE],
    snp_ranges = annot@snp_ranges,
    annotation_meta = annot@annotation_meta[idx, , drop = FALSE],
    genome = annot@genome
  )
}

#' @title Get Candidate Annotations
#' @description Extract only candidate-tier annotations from an
#'   \code{AnnotationMatrix}.
#' @param annot An \code{AnnotationMatrix} object.
#' @return An \code{AnnotationMatrix} with only candidate annotations.
#' @export
getCandidates <- function(annot) {
  idx <- annot@annotation_meta$tier == "candidate"
  AnnotationMatrix(
    annotations = annot@annotations[, idx, drop = FALSE],
    snp_ranges = annot@snp_ranges,
    annotation_meta = annot@annotation_meta[idx, , drop = FALSE],
    genome = annot@genome
  )
}
