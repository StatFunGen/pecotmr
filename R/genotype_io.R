#' @title Genotype I/O via GenotypeHandle
#' @description Read genotype data from various formats (VCF, plink1,
#'   plink2, GDS) and provide block-level genotype extraction without
#'   requiring format conversion.
#' @name pecotmr-genotype-io
#' @keywords internal
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
      format <- .h2DetectFormat(path)
    }
    switch(format,
      "gds" = .makeGdsHandle(path),
      "vcf" = .makeVcfHandle(path, ...),
      "plink1" = .makePlink1Handle(path, ...),
      "plink2" = .makePlink2Handle(path, ...),
      stop("Unsupported genotype format: ", format)
    )
  }
)

# =============================================================================
# Handle constructors — read metadata, defer genotype loading
# =============================================================================

#' @keywords internal
.makeGdsHandle <- function(path) {
  if (!requireNamespace("SNPRelate", quietly = TRUE))
    stop("Package 'SNPRelate' is required for reading GDS files.")
  if (!requireNamespace("gdsfmt", quietly = TRUE))
    stop("Package 'gdsfmt' is required for reading GDS files.")
  if (!file.exists(path))
    stop("GDS file not found: ", path)

  snpInfo <- .gdsSnpInfo(path)

  gds <- SNPRelate::snpgdsOpen(path, readonly = TRUE)
  on.exit(SNPRelate::snpgdsClose(gds))
  sampleIds <- as.character(gdsfmt::read.gdsn(
    gdsfmt::index.gdsn(gds, "sample.id")))
  nSamples <- length(sampleIds)

  new("GenotypeHandle",
    path = path,
    format = "gds",
    snpInfo = snpInfo,
    nSamples = as.integer(nSamples),
    sampleIds = sampleIds,
    pgenPtr = NULL
  )
}

#' @keywords internal
.makeVcfHandle <- function(path, ...) {
  if (!requireNamespace("VariantAnnotation", quietly = TRUE))
    stop("Package 'VariantAnnotation' is required for reading VCF files.")
  if (!file.exists(path))
    stop("VCF file not found: ", path)

  hdr <- VariantAnnotation::scanVcfHeader(path)
  sampleIds <- as.character(VariantAnnotation::samples(hdr))
  nSamples <- length(sampleIds)

  param <- VariantAnnotation::ScanVcfParam(fixed = c("ALT"), info = NA,
                                            geno = NA)
  vcf <- VariantAnnotation::readVcf(path, param = param, ...)
  rd <- rowRanges(vcf)

  # pecotmr convention: A1 = ALT (effect), A2 = REF
  snpInfo <- data.frame(
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
    snpInfo = snpInfo,
    nSamples = as.integer(nSamples),
    sampleIds = sampleIds,
    pgenPtr = NULL
  )
}

#' @keywords internal
.makePlink1Handle <- function(path, ...) {
  if (!requireNamespace("snpStats", quietly = TRUE))
    stop("Package 'snpStats' is required for reading plink1 files.")

  stem <- .plinkStem(path)
  bedFile <- paste0(stem, ".bed")
  bimFile <- paste0(stem, ".bim")
  famFile <- paste0(stem, ".fam")

  for (f in c(bedFile, bimFile, famFile)) {
    if (!file.exists(f))
      stop("Plink file not found: ", f)
  }

  bim <- read.table(bimFile, header = FALSE,
                            stringsAsFactors = FALSE,
                            col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))
  fam <- read.table(famFile, header = FALSE, stringsAsFactors = FALSE)
  sampleIds <- as.character(fam[, 2])  # IID column
  nSamples <- nrow(fam)

  # plink1 bim: col5 = A1 (minor/effect), col6 = A2 (major/ref)
  # Matches pecotmr convention directly
  snpInfo <- data.frame(
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
    snpInfo = snpInfo,
    nSamples = as.integer(nSamples),
    sampleIds = sampleIds,
    pgenPtr = NULL
  )
}

#' @keywords internal
.makePlink2Handle <- function(path, ...) {
  if (!requireNamespace("pgenlibr", quietly = TRUE))
    stop("Package 'pgenlibr' is required for reading plink2 files.")

  stem <- .plinkStem(path)

  # Use pecotmr's resolvePlink2Paths for robust path detection (.pvar.zst)
  paths <- resolvePlink2Paths(stem)

  # Use pecotmr's readPvar for robust .pvar/.pvar.zst handling via pgenlibr
  vi <- readPvar(paths$pvar)
  # readPvar returns: chrom, id, pos, A2 (REF), A1 (ALT) — pecotmr convention
  snpInfo <- data.frame(
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
  sampleIds <- as.character(psam$IID)

  pgen <- pgenlibr::NewPgen(paths$pgen)
  nSamples <- pgenlibr::GetRawSampleCt(pgen)

  new("GenotypeHandle",
    path = stem,
    format = "plink2",
    snpInfo = snpInfo,
    nSamples = as.integer(nSamples),
    sampleIds = sampleIds,
    pgenPtr = pgen
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
#' @param snpIdx Integer vector of 1-based SNP indices into
#'   \code{handle@@snpInfo}.
#' @param meanImpute Logical, whether to mean-impute missing values.
#'   Default TRUE.
#' @return A \code{RangedSummarizedExperiment} with:
#'   \describe{
#'     \item{assay("dosage")}{Numeric matrix (variants x samples)}
#'     \item{rowRanges}{GRanges with A1, A2 metadata}
#'     \item{colData}{DataFrame with sampleId column}
#'   }
#' @export
extractBlockGenotypes <- function(handle, snpIdx, meanImpute = TRUE) {
  geno <- switch(handle@format,
    "gds" = .extractBlockGds(handle, snpIdx),
    "vcf" = .extractBlockVcf(handle, snpIdx),
    "plink1" = .extractBlockPlink1(handle, snpIdx),
    "plink2" = .extractBlockPlink2(handle, snpIdx),
    stop("Unsupported format in extractBlockGenotypes: ", handle@format)
  )
  if (is.null(geno)) return(NULL)
  if (meanImpute) geno <- .meanImputeGeno(geno)

  # geno is samples x variants from the format extractors
  si <- handle@snpInfo[snpIdx, , drop = FALSE]
  chr <- as.character(si$CHR)
  chr <- sub("^chr", "", chr, ignore.case = TRUE)
  chr <- paste0("chr", chr)

  rowRanges <- GRanges(
    seqnames = chr,
    ranges = IRanges(start = as.integer(si$BP), width = 1L)
  )
  mcols(rowRanges) <- DataFrame(
    SNP = si$SNP, A1 = si$A1, A2 = si$A2
  )

  colData <- DataFrame(
    sampleId = handle@sampleIds,
    row.names = handle@sampleIds
  )

  # Transpose to Bioc convention: variants x samples
  dosage <- t(geno)
  rownames(dosage) <- si$SNP
  colnames(dosage) <- handle@sampleIds

  SummarizedExperiment(
    assays = list(dosage = dosage),
    rowRanges = rowRanges,
    colData = colData
  )
}

#' @keywords internal
.extractBlockGds <- function(handle, snpIdx) {
  gds <- SNPRelate::snpgdsOpen(handle@path, readonly = TRUE,
                                allow.fork = TRUE)
  on.exit(SNPRelate::snpgdsClose(gds))

  snpIds <- handle@snpInfo$SNP[snpIdx]
  # Use snpgdsGetGeno for proper non-contiguous SNP selection
  geno <- SNPRelate::snpgdsGetGeno(gds, snp.id = snpIds,
                                    with.id = FALSE, verbose = FALSE)
  if (is.null(geno) || length(geno) == 0) return(NULL)

  # snpgdsGetGeno returns count of the first allele in snp.allele,
  # which we label A1 in .gdsSnpInfo. No flip needed.
  storage.mode(geno) <- "double"
  geno
}

#' @keywords internal
.extractBlockVcf <- function(handle, snpIdx) {
  si <- handle@snpInfo[snpIdx, ]
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
.extractBlockPlink1 <- function(handle, snpIdx) {
  snpIds <- handle@snpInfo$SNP[snpIdx]
  plinkData <- snpStats::read.plink(
    bed = paste0(handle@path, ".bed"),
    bim = paste0(handle@path, ".bim"),
    fam = paste0(handle@path, ".fam"),
    select.snps = snpIds
  )
  # snpStats as(x, "numeric") gives count of B allele (A2/bim col 6).
  # Flip to count A1 (bim col 5 / effect allele).
  geno <- 2 - as(plinkData$genotypes, "numeric")
  storage.mode(geno) <- "double"
  geno
}

#' @keywords internal
.extractBlockPlink2 <- function(handle, snpIdx) {
  # pgenlibr::ReadList returns ALT dosage = A1 dosage in pecotmr convention
  geno <- pgenlibr::ReadList(handle@pgenPtr, variant_subset = snpIdx,
                              meanimpute = FALSE)
  storage.mode(geno) <- "double"
  geno
}

# =============================================================================
# LD correlation computation
# =============================================================================

#' @title Compute Block LD Correlation
#' @description Compute the LD correlation matrix for a block of SNPs.
#'   Delegates to \code{\link{computeLd}} for the actual computation,
#'   with automatic backend selection based on file format unless
#'   overridden.
#' @param handle A \code{GenotypeHandle} object.
#' @param snpIdx Integer vector of 1-based SNP indices.
#' @param backend Character, one of \code{"internal"} (default),
#'   \code{"snprelate"}, or \code{"snpstats"}. When \code{"internal"},
#'   GDS-format handles automatically use \code{SNPRelate::snpgdsLDMat}
#'   via the native GDS path; other formats use the internal correlator.
#' @param method Character, LD computation method passed to
#'   \code{\link{computeLd}}. Default \code{"sample"}.
#' @param ... Additional arguments passed to \code{\link{computeLd}}
#'   (e.g., \code{shrinkage}, \code{trimSamples}).
#' @return Numeric correlation matrix (p x p).
#' @export
computeBlockLdCor <- function(handle, snpIdx, backend = "internal",
                              method = "sample", ...) {
  # For GDS format with internal backend, use the native SNPRelate path
  # which avoids extracting genotypes into memory
  if (handle@format == "gds" && backend == "internal") {
    return(.computeBlockLdGds(handle, snpIdx))
  }

  # Extract genotypes via the unified GenotypeHandle pipeline
  rse <- extractBlockGenotypes(handle, snpIdx)
  if (is.null(rse)) return(diag(length(snpIdx)))
  geno <- t(SummarizedExperiment::assay(rse, "dosage"))
  if (ncol(geno) < 2) return(diag(length(snpIdx)))

  # Delegate to computeLd for all computation
  computeLd(geno, method = method, backend = backend, ...)
}

#' @keywords internal
.computeBlockLdGds <- function(handle, snpIdx) {
  gds <- SNPRelate::snpgdsOpen(handle@path, readonly = TRUE,
                                allow.fork = TRUE)
  on.exit(SNPRelate::snpgdsClose(gds))

  snpIds <- handle@snpInfo$SNP[snpIdx]
  ldMat <- SNPRelate::snpgdsLDMat(
    gds, snp.id = snpIds, method = "corr",
    slide = -1, verbose = FALSE
  )
  R <- ldMat$LD
  R[is.na(R)] <- 0
  R
}

# =============================================================================
# Region filtering helper
# =============================================================================

#' @title Filter SNP Info by Region
#' @description Return 1-based indices into snpInfo for SNPs within a
#'   genomic region string.
#' @param snpInfo data.frame with CHR and BP columns.
#' @param region Character region string "chr:start-end".
#' @return Integer vector of matching SNP indices.
#' @keywords internal
.regionToSnpIdx <- function(snpInfo, region) {
  parsed <- parseRegion(region)
  chrMatch <- stripChrPrefix(as.character(snpInfo$CHR)) == parsed$chrom
  posMatch <- snpInfo$BP >= parsed$start & snpInfo$BP <= parsed$end
  which(chrMatch & posMatch)
}

#' @title Convert SNP Info to Variant Info
#' @description Convert GenotypeHandle snpInfo (uppercase columns) to
#'   pecotmr variant_info format (lowercase columns).
#' @param snpInfo data.frame with SNP, CHR, BP, A1, A2 columns.
#' @return data.frame with chrom, id, pos, A2, A1 columns.
#' @keywords internal
.snpInfoToVariantInfo <- function(snpInfo) {
  data.frame(
    chrom = snpInfo$CHR,
    id = snpInfo$SNP,
    pos = snpInfo$BP,
    A2 = snpInfo$A2,
    A1 = snpInfo$A1,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# Helpers
# =============================================================================

#' @keywords internal
.meanImputeGeno <- function(geno) {
  naCols <- which(colSums(is.na(geno)) > 0L)
  for (j in naCols) {
    colMean <- mean(geno[, j], na.rm = TRUE)
    geno[is.na(geno[, j]), j] <- colMean
  }
  geno
}

#' @keywords internal
.gdsSnpInfo <- function(gdsPath) {
  gds <- SNPRelate::snpgdsOpen(gdsPath, readonly = TRUE)
  on.exit(SNPRelate::snpgdsClose(gds))

  snpId <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.id"))
  chr <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.chromosome"))
  pos <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.position"))
  allele <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.allele"))

  allelesSplit <- strsplit(allele, "/")
  # snpgdsGetGeno counts copies of the first allele in snp.allele.
  # Label the first allele as A1 so dosage = count of A1.
  a1 <- vapply(allelesSplit, `[`, character(1), 1L)
  a2 <- vapply(allelesSplit, `[`, character(1), 2L)

  data.frame(
    SNP = snpId,
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
.h2DetectFormat <- function(path) {
  lpath <- tolower(path)
  if (grepl("\\.vcf\\.gz$", lpath) || grepl("\\.vcf\\.bgz$", lpath))
    return("vcf")
  if (grepl("\\.annot\\.gz$", lpath))
    return("ldsc_annot")

  ext <- tolower(file_ext(path))
  if (nzchar(ext)) {
    detected <- switch(ext,
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
      NULL
    )
    if (!is.null(detected)) return(detected)
  }
  # Check for file stems, including dotted prefixes such as sample.EUR.chr21.
  if (file.exists(paste0(path, ".pgen")) || file.exists(paste0(path, ".pvar")))
    return("plink2")
  if (file.exists(paste0(path, ".bed")) || file.exists(paste0(path, ".bim")))
    return("plink1")
  if (file.exists(paste0(path, ".gds")))
    return("gds")
  if (nzchar(ext))
    stop("Cannot detect format from extension: ", ext)
  stop("Cannot detect genotype format for path: ", path)
}

#' @title Detect Plink File Stem
#' @description Given any plink file path, return the stem.
#' @param path Character, path to any plink file.
#' @return Character, file stem without extension.
#' @keywords internal
.plinkStem <- function(path) {
  # Only strip known plink extensions; leave other paths as-is (they may
  # already be the stem, e.g. "prefix.genotype" -> "prefix.genotype.bed")
  ext <- file_ext(path)
  plinkExts <- c("bed", "bim", "fam", "pgen", "pvar", "psam")
  if (tolower(ext) %in% plinkExts) {
    file_path_sans_ext(path)
  } else {
    path
  }
}
