#' @title Genotype I/O via GenotypeHandle
#' @description Read genotype data from various formats (VCF, plink1,
#'   plink2, GDS) and provide block-level genotype extraction without
#'   requiring format conversion.
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom S4Vectors DataFrame mcols mcols<-
#' @importFrom tools file_ext
#' @importFrom methods as
#' @include AllGenerics.R
NULL

# =============================================================================
# Main reader method — returns a GenotypeHandle
# =============================================================================

#' @rdname readGenotypes
#' @export
setMethod("readGenotypes",
  signature(path = "character"),
  function(path, format = NULL, ...) {
    if (is.null(format)) {
      format <- .h2_detect_format(path)
    }
    switch(format,
      "gds" = .make_gds_handle(path),
      "vcf" = .make_vcf_handle(path, ...),
      "plink1" = .make_plink1_handle(path, ...),
      "plink2" = .make_plink2_handle(path, ...),
      stop("Unsupported genotype format: ", format)
    )
  }
)

# =============================================================================
# Handle constructors — read metadata, defer genotype loading
# =============================================================================

#' @keywords internal
.make_gds_handle <- function(path) {
  if (!requireNamespace("SNPRelate", quietly = TRUE))
    stop("Package 'SNPRelate' is required for reading GDS files.")
  if (!requireNamespace("gdsfmt", quietly = TRUE))
    stop("Package 'gdsfmt' is required for reading GDS files.")
  if (!file.exists(path))
    stop("GDS file not found: ", path)

  snp_info <- .gds_snp_info(path)

  gds <- SNPRelate::snpgdsOpen(path, readonly = TRUE)
  on.exit(SNPRelate::snpgdsClose(gds))
  sample_ids <- as.character(gdsfmt::read.gdsn(
    gdsfmt::index.gdsn(gds, "sample.id")))
  n_samples <- length(sample_ids)

  new("GenotypeHandle",
    path = path,
    format = "gds",
    snp_info = snp_info,
    n_samples = as.integer(n_samples),
    sample_ids = sample_ids,
    pgen_ptr = NULL
  )
}

#' @keywords internal
.make_vcf_handle <- function(path, ...) {
  if (!requireNamespace("VariantAnnotation", quietly = TRUE))
    stop("Package 'VariantAnnotation' is required for reading VCF files.")
  if (!file.exists(path))
    stop("VCF file not found: ", path)

  hdr <- VariantAnnotation::scanVcfHeader(path)
  sample_ids <- as.character(VariantAnnotation::samples(hdr))
  n_samples <- length(sample_ids)

  param <- VariantAnnotation::ScanVcfParam(fixed = c("ALT"), info = NA,
                                            geno = NA)
  vcf <- VariantAnnotation::readVcf(path, param = param, ...)
  rd <- rowRanges(vcf)

  # pecotmr convention: A1 = ALT (effect), A2 = REF
  snp_info <- data.frame(
    SNP = names(rd),
    CHR = as.character(seqnames(rd)),
    BP = as.integer(start(rd)),
    A1 = vapply(rd$ALT, function(x) as.character(x)[1], character(1)),
    A2 = as.character(rd$REF),
    stringsAsFactors = FALSE
  )

  new("GenotypeHandle",
    path = normalizePath(path),
    format = "vcf",
    snp_info = snp_info,
    n_samples = as.integer(n_samples),
    sample_ids = sample_ids,
    pgen_ptr = NULL
  )
}

#' @keywords internal
.make_plink1_handle <- function(path, ...) {
  if (!requireNamespace("snpStats", quietly = TRUE))
    stop("Package 'snpStats' is required for reading plink1 files.")

  stem <- .plink_stem(path)
  bed_file <- paste0(stem, ".bed")
  bim_file <- paste0(stem, ".bim")
  fam_file <- paste0(stem, ".fam")

  for (f in c(bed_file, bim_file, fam_file)) {
    if (!file.exists(f))
      stop("Plink file not found: ", f)
  }

  bim <- read.table(bim_file, header = FALSE,
                            stringsAsFactors = FALSE,
                            col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))
  fam <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
  sample_ids <- as.character(fam[, 2])  # IID column
  n_samples <- nrow(fam)

  # plink1 bim: col5 = A1 (minor/effect), col6 = A2 (major/ref)
  # Matches pecotmr convention directly
  snp_info <- data.frame(
    SNP = bim$SNP,
    CHR = as.character(bim$CHR),
    BP = as.integer(bim$BP),
    A1 = bim$A1,
    A2 = bim$A2,
    stringsAsFactors = FALSE
  )

  new("GenotypeHandle",
    path = stem,
    format = "plink1",
    snp_info = snp_info,
    n_samples = as.integer(n_samples),
    sample_ids = sample_ids,
    pgen_ptr = NULL
  )
}

#' @keywords internal
.make_plink2_handle <- function(path, ...) {
  if (!requireNamespace("pgenlibr", quietly = TRUE))
    stop("Package 'pgenlibr' is required for reading plink2 files.")

  stem <- .plink_stem(path)

  # Use pecotmr's resolve_plink2_paths for robust path detection (.pvar.zst)
  paths <- resolve_plink2_paths(stem)

  # Use pecotmr's read_pvar for robust .pvar/.pvar.zst handling via pgenlibr
  vi <- read_pvar(paths$pvar)
  # read_pvar returns: chrom, id, pos, A2 (REF), A1 (ALT) — pecotmr convention
  snp_info <- data.frame(
    SNP = vi$id,
    CHR = as.character(vi$chrom),
    BP = as.integer(vi$pos),
    A1 = vi$A1,
    A2 = vi$A2,
    stringsAsFactors = FALSE
  )

  # Read sample IDs from .psam
  psam <- as.data.frame(vroom(paths$psam, delim = "\t",
                                      show_col_types = FALSE))
  names(psam) <- sub("^#", "", names(psam))
  sample_ids <- as.character(psam$IID)

  pgen <- pgenlibr::NewPgen(paths$pgen)
  n_samples <- pgenlibr::GetRawSampleCt(pgen)

  new("GenotypeHandle",
    path = stem,
    format = "plink2",
    snp_info = snp_info,
    n_samples = as.integer(n_samples),
    sample_ids = sample_ids,
    pgen_ptr = pgen
  )
}

# =============================================================================
# Block genotype extraction — dispatches by format
# =============================================================================

#' @title Extract Block Genotypes
#' @description Extract a genotype matrix for a subset of SNPs from a
#'   \code{GenotypeHandle}. Returns a \code{RangedSummarizedExperiment}
#'   with dosage assay in Bioconductor convention (variants x samples),
#'   variant metadata as \code{rowRanges} (GRanges), and sample IDs as
#'   \code{colData}.
#' @param handle A \code{GenotypeHandle} object.
#' @param snp_idx Integer vector of 1-based SNP indices into
#'   \code{handle@@snp_info}.
#' @param mean_impute Logical, whether to mean-impute missing values.
#'   Default TRUE.
#' @return A \code{RangedSummarizedExperiment} with:
#'   \describe{
#'     \item{assay("dosage")}{Numeric matrix (variants x samples)}
#'     \item{rowRanges}{GRanges with A1, A2 metadata}
#'     \item{colData}{DataFrame with sample_id column}
#'   }
#' @export
extractBlockGenotypes <- function(handle, snp_idx, mean_impute = TRUE) {
  geno <- switch(handle@format,
    "gds" = .extract_block_gds(handle, snp_idx),
    "vcf" = .extract_block_vcf(handle, snp_idx),
    "plink1" = .extract_block_plink1(handle, snp_idx),
    "plink2" = .extract_block_plink2(handle, snp_idx),
    stop("Unsupported format in extractBlockGenotypes: ", handle@format)
  )
  if (is.null(geno)) return(NULL)
  if (mean_impute) geno <- .mean_impute_geno(geno)

  # geno is samples x variants from the format extractors
  si <- handle@snp_info[snp_idx, , drop = FALSE]
  chr <- as.character(si$CHR)
  chr <- sub("^chr", "", chr, ignore.case = TRUE)
  chr <- paste0("chr", chr)

  row_ranges <- GRanges(
    seqnames = chr,
    ranges = IRanges(start = as.integer(si$BP), width = 1L)
  )
  mcols(row_ranges) <- DataFrame(
    SNP = si$SNP, A1 = si$A1, A2 = si$A2
  )

  col_data <- DataFrame(
    sample_id = handle@sample_ids,
    row.names = handle@sample_ids
  )

  # Transpose to Bioc convention: variants x samples
  dosage <- t(geno)
  rownames(dosage) <- si$SNP
  colnames(dosage) <- handle@sample_ids

  SummarizedExperiment(
    assays = list(dosage = dosage),
    rowRanges = row_ranges,
    colData = col_data
  )
}

#' @keywords internal
.extract_block_gds <- function(handle, snp_idx) {
  gds <- SNPRelate::snpgdsOpen(handle@path, readonly = TRUE,
                                allow.fork = TRUE)
  on.exit(SNPRelate::snpgdsClose(gds))

  snp_ids <- handle@snp_info$SNP[snp_idx]
  # Use snpgdsGetGeno for proper non-contiguous SNP selection
  geno <- SNPRelate::snpgdsGetGeno(gds, snp.id = snp_ids,
                                    with.id = FALSE, verbose = FALSE)
  if (is.null(geno) || length(geno) == 0) return(NULL)

  # snpgdsGetGeno returns count of the first allele in snp.allele,
  # which we label A1 in .gds_snp_info. No flip needed.
  storage.mode(geno) <- "double"
  geno
}

#' @keywords internal
.extract_block_vcf <- function(handle, snp_idx) {
  si <- handle@snp_info[snp_idx, ]
  gr <- GRanges(
    seqnames = si$CHR,
    ranges = IRanges(start = si$BP, end = si$BP)
  )
  param <- VariantAnnotation::ScanVcfParam(which = gr, fixed = NA,
                                            info = NA, geno = "GT")
  vcf <- VariantAnnotation::readVcf(handle@path, genome = "", param = param)
  gt <- VariantAnnotation::geno(vcf)$GT

  # Convert GT strings to ALT dosage (A1 dosage)
  geno <- matrix(NA_real_, nrow = ncol(gt), ncol = nrow(gt))
  for (j in seq_len(nrow(gt))) {
    g <- gt[j, ]
    geno[, j] <- vapply(g, function(x) {
      if (is.na(x) || x == "./.") return(NA_real_)
      alleles <- strsplit(x, "[/|]")[[1]]
      sum(alleles != "0")
    }, numeric(1))
  }

  geno
}

#' @keywords internal
.extract_block_plink1 <- function(handle, snp_idx) {
  snp_ids <- handle@snp_info$SNP[snp_idx]
  plink_data <- snpStats::read.plink(
    bed = paste0(handle@path, ".bed"),
    bim = paste0(handle@path, ".bim"),
    fam = paste0(handle@path, ".fam"),
    select.snps = snp_ids
  )
  # snpStats as(x, "numeric") gives count of B allele (A2/bim col 6).
  # Flip to count A1 (bim col 5 / effect allele).
  geno <- 2 - as(plink_data$genotypes, "numeric")
  storage.mode(geno) <- "double"
  geno
}

#' @keywords internal
.extract_block_plink2 <- function(handle, snp_idx) {
  # pgenlibr::ReadList returns ALT dosage = A1 dosage in pecotmr convention
  geno <- pgenlibr::ReadList(handle@pgen_ptr, variant_subset = snp_idx,
                              meanimpute = FALSE)
  storage.mode(geno) <- "double"
  geno
}

# =============================================================================
# LD correlation computation
# =============================================================================

#' @title Compute Block LD Correlation
#' @description Compute the LD correlation matrix for a block of SNPs.
#'   Delegates to \code{\link{compute_LD}} for the actual computation,
#'   with automatic backend selection based on file format unless
#'   overridden.
#' @param handle A \code{GenotypeHandle} object.
#' @param snp_idx Integer vector of 1-based SNP indices.
#' @param backend Character, one of \code{"internal"} (default),
#'   \code{"snprelate"}, or \code{"snpstats"}. When \code{"internal"},
#'   GDS-format handles automatically use \code{SNPRelate::snpgdsLDMat}
#'   via the native GDS path; other formats use the internal correlator.
#' @param method Character, LD computation method passed to
#'   \code{\link{compute_LD}}. Default \code{"sample"}.
#' @param ... Additional arguments passed to \code{\link{compute_LD}}
#'   (e.g., \code{shrinkage}, \code{trim_samples}).
#' @return Numeric correlation matrix (p x p).
#' @export
computeBlockLdCor <- function(handle, snp_idx, backend = "internal",
                              method = "sample", ...) {
  # For GDS format with internal backend, use the native SNPRelate path
  # which avoids extracting genotypes into memory
  if (handle@format == "gds" && backend == "internal") {
    return(.compute_block_ld_gds(handle, snp_idx))
  }

  # Extract genotypes via the unified GenotypeHandle pipeline
  rse <- extractBlockGenotypes(handle, snp_idx)
  if (is.null(rse)) return(diag(length(snp_idx)))
  geno <- t(SummarizedExperiment::assay(rse, "dosage"))
  if (ncol(geno) < 2) return(diag(length(snp_idx)))

  # Delegate to compute_LD for all computation
  compute_LD(geno, method = method, backend = backend, ...)
}

#' @keywords internal
.compute_block_ld_gds <- function(handle, snp_idx) {
  gds <- SNPRelate::snpgdsOpen(handle@path, readonly = TRUE,
                                allow.fork = TRUE)
  on.exit(SNPRelate::snpgdsClose(gds))

  snp_ids <- handle@snp_info$SNP[snp_idx]
  ld_mat <- SNPRelate::snpgdsLDMat(
    gds, snp.id = snp_ids, method = "corr",
    slide = -1, verbose = FALSE
  )
  R <- ld_mat$LD
  R[is.na(R)] <- 0
  R
}

# =============================================================================
# Region filtering helper
# =============================================================================

#' @title Filter SNP Info by Region
#' @description Return 1-based indices into snp_info for SNPs within a
#'   genomic region string.
#' @param snp_info data.frame with CHR and BP columns.
#' @param region Character region string "chr:start-end".
#' @return Integer vector of matching SNP indices.
#' @keywords internal
.region_to_snp_idx <- function(snp_info, region) {
  parsed <- parse_region(region)
  chr_match <- strip_chr_prefix(as.character(snp_info$CHR)) == parsed$chrom
  pos_match <- snp_info$BP >= parsed$start & snp_info$BP <= parsed$end
  which(chr_match & pos_match)
}

#' @title Convert SNP Info to Variant Info
#' @description Convert GenotypeHandle snp_info (uppercase columns) to
#'   pecotmr variant_info format (lowercase columns).
#' @param snp_info data.frame with SNP, CHR, BP, A1, A2 columns.
#' @return data.frame with chrom, id, pos, A2, A1 columns.
#' @keywords internal
.snp_info_to_variant_info <- function(snp_info) {
  data.frame(
    chrom = snp_info$CHR,
    id = snp_info$SNP,
    pos = snp_info$BP,
    A2 = snp_info$A2,
    A1 = snp_info$A1,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# Helpers
# =============================================================================

#' @keywords internal
.mean_impute_geno <- function(geno) {
  na_cols <- which(colSums(is.na(geno)) > 0L)
  for (j in na_cols) {
    col_mean <- mean(geno[, j], na.rm = TRUE)
    geno[is.na(geno[, j]), j] <- col_mean
  }
  geno
}

#' @keywords internal
.gds_snp_info <- function(gds_path) {
  gds <- SNPRelate::snpgdsOpen(gds_path, readonly = TRUE)
  on.exit(SNPRelate::snpgdsClose(gds))

  snp_id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.id"))
  chr <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.chromosome"))
  pos <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.position"))
  allele <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.allele"))

  alleles_split <- strsplit(allele, "/")
  # snpgdsGetGeno counts copies of the first allele in snp.allele.
  # Label the first allele as A1 so dosage = count of A1.
  a1 <- vapply(alleles_split, `[`, character(1), 1L)
  a2 <- vapply(alleles_split, `[`, character(1), 2L)

  data.frame(
    SNP = snp_id,
    CHR = as.character(chr),
    BP = as.integer(pos),
    A1 = a1,
    A2 = a2,
    stringsAsFactors = FALSE
  )
}

#' @title Detect File Format from Extension
#' @description Infer file format from the file extension.
#' @param path Character, file path.
#' @return Character, detected format.
#' @keywords internal
.h2_detect_format <- function(path) {
  lpath <- tolower(path)
  if (grepl("\\.vcf\\.gz$", lpath) || grepl("\\.vcf\\.bgz$", lpath))
    return("vcf")
  if (grepl("\\.annot\\.gz$", lpath))
    return("ldsc_annot")

  ext <- tolower(file_ext(path))
  if (nzchar(ext)) {
    return(switch(ext,
      "vcf" = "vcf",
      "bcf" = "vcf",
      "bed" = "plink1",
      "bim" = "plink1",
      "fam" = "plink1",
      "pgen" = "plink2",
      "pvar" = "plink2",
      "psam" = "plink2",
      "gds" = "gds",
      "rds" = "rds",
      "rdata" = "rds",
      "annot" = "ldsc_annot",
      "bw" = "bigwig",
      "bigwig" = "bigwig",
      stop("Cannot detect format from extension: ", ext)
    ))
  }
  # No extension — check for plink stem files
  if (file.exists(paste0(path, ".pgen")) || file.exists(paste0(path, ".pvar")))
    return("plink2")
  if (file.exists(paste0(path, ".bed")) || file.exists(paste0(path, ".bim")))
    return("plink1")
  if (file.exists(paste0(path, ".gds")))
    return("gds")
  stop("Cannot detect genotype format for path: ", path)
}

#' @title Detect Plink File Stem
#' @description Given any plink file path, return the stem.
#' @param path Character, path to any plink file.
#' @return Character, file stem without extension.
#' @keywords internal
.plink_stem <- function(path) {
  # Only strip known plink extensions; leave other paths as-is (they may
  # already be the stem, e.g. "prefix.genotype" -> "prefix.genotype.bed")
  ext <- file_ext(path)
  plink_exts <- c("bed", "bim", "fam", "pgen", "pvar", "psam")
  if (tolower(ext) %in% plink_exts) {
    file_path_sans_ext(path)
  } else {
    path
  }
}
