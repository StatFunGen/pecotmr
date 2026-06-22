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
  fmt <- getFormat(handle)
  geno <- switch(fmt,
    "gds" = .extractBlockGds(handle, snpIdx),
    "vcf" = .extractBlockVcf(handle, snpIdx),
    "plink1" = .extractBlockPlink1(handle, snpIdx),
    "plink2" = .extractBlockPlink2(handle, snpIdx),
    stop("Unsupported format in extractBlockGenotypes: ", fmt)
  )
  if (is.null(geno)) return(NULL)
  if (meanImpute) geno <- .meanImputeGeno(geno)

  # geno is samples x variants from the format extractors
  si <- getSnpInfo(handle)[snpIdx, , drop = FALSE]
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

  sampleIds <- getSampleIds(handle)
  colData <- DataFrame(
    sampleId = sampleIds,
    row.names = sampleIds
  )

  # Transpose to Bioc convention: variants x samples
  dosage <- t(geno)
  rownames(dosage) <- si$SNP
  colnames(dosage) <- sampleIds

  SummarizedExperiment(
    assays = list(dosage = dosage),
    rowRanges = rowRanges,
    colData = colData
  )
}

#' @keywords internal
.extractBlockGds <- function(handle, snpIdx) {
  gds <- SNPRelate::snpgdsOpen(getPath(handle), readonly = TRUE,
                                allow.fork = TRUE)
  on.exit(SNPRelate::snpgdsClose(gds))

  snpIds <- getSnpInfo(handle)$SNP[snpIdx]
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
  si <- getSnpInfo(handle)[snpIdx, ]
  gr <- GRanges(
    seqnames = si$CHR,
    ranges = IRanges(start = si$BP, end = si$BP)
  )
  param <- VariantAnnotation::ScanVcfParam(which = gr, fixed = NA,
                                            info = NA, geno = "GT")
  vcf <- VariantAnnotation::readVcf(getPath(handle), genome = "", param = param)
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
  snpIds <- getSnpInfo(handle)$SNP[snpIdx]
  pathStem <- getPath(handle)
  plinkData <- snpStats::read.plink(
    bed = paste0(pathStem, ".bed"),
    bim = paste0(pathStem, ".bim"),
    fam = paste0(pathStem, ".fam"),
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
  geno <- pgenlibr::ReadList(getPgenPtr(handle), variant_subset = snpIdx,
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
  if (getFormat(handle) == "gds" && backend == "internal") {
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
  gds <- SNPRelate::snpgdsOpen(getPath(handle), readonly = TRUE,
                                allow.fork = TRUE)
  on.exit(SNPRelate::snpgdsClose(gds))

  snpIds <- getSnpInfo(handle)$SNP[snpIdx]
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


# =============================================================================
# Format-specific file readers
# -----------------------------------------------------------------------------
# Low-level PLINK / VCF / GDS variant-metadata readers, the stochastic
# genotype sidecar helpers (.afreq / .stochastic_meta.tsv), and the
# top-level dispatchers loadGenotypeRegion + getRefVariantInfo that
# auto-detect the underlying format and route to the correct reader.
# =============================================================================

# read PLINK files

#' @importFrom vroom vroom
#' @importFrom tools file_path_sans_ext
#' @importFrom Rsamtools TabixFile seqnamesTabix scanTabix headerTabix
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom SummarizedExperiment assay
#' @importFrom MungeSumstats standardise_header
readBim <- function(bed) {
  bimf <- paste0(file_path_sans_ext(bed), ".bim")
  bim <- vroom(bimf, col_names = FALSE)
  colnames(bim) <- c("chrom", "id", "gpos", "pos", "a1", "a0")
  return(bim)
}

#' @importFrom vroom vroom
#' @importFrom tools file_path_sans_ext
readFam <- function(bed) {
  famf <- paste0(file_path_sans_ext(bed), ".fam")
  return(vroom(famf, col_names = FALSE))
}

# open bed/bim/fam: A PLINK 1 .bed is a valid .pgen
openBed <- function(bed) {
  if (!requireNamespace("pgenlibr", quietly = TRUE)) {
    stop("To use this function, please install pgenlibr: https://cran.r-project.org/web/packages/pgenlibr/index.html")
  }
  rawSCt <- nrow(readFam(bed))
  return(pgenlibr::NewPgen(bed, raw_sample_ct = rawSCt))
}

#' Read a PLINK2 allele frequency file (.afreq or .afreq.zst)
#'
#' @param prefix File prefix (without .afreq extension).
#' @return A data.frame with columns: chrom, id, A2 (REF), A1 (ALT), alt_freq, obs_ct.
#'   alt_freq is the frequency of the A1 (ALT/effect) allele.
#' @importFrom vroom vroom
#' @importFrom dplyr rename select
#' @export
readAfreq <- function(prefix) {
  afreqZst <- paste0(prefix, ".afreq.zst")
  afreqPlain <- paste0(prefix, ".afreq")
  if (file.exists(afreqZst)) {
    if (Sys.which("zstd") == "") stop("zstd CLI is required to read .afreq.zst files")
    af <- as.data.frame(vroom(pipe(paste0("zstd -dcq ", shQuote(afreqZst))),
                              delim = "\t", show_col_types = FALSE))
  } else if (file.exists(afreqPlain)) {
    af <- as.data.frame(vroom(afreqPlain, delim = "\t", show_col_types = FALSE))
  } else {
    return(NULL)
  }
  # PLINK2 .afreq: REF = A2, ALT = A1, ALT_FREQS = A1 (effect allele) frequency
  af <- rename(af,
    "chrom" = "#CHROM", "id" = "ID",
    "A2" = "REF", "A1" = "ALT",
    "alt_freq" = "ALT_FREQS", "obs_ct" = "OBS_CT"
  )
  cols <- c("chrom", "id", "A2", "A1", "alt_freq", "obs_ct")
  # Stochastic genotype .afreq includes U_MIN/U_MAX for exact min-max inversion
  if ("U_MIN" %in% colnames(af)) {
    af <- rename(af, "u_min" = "U_MIN", "u_max" = "U_MAX")
    cols <- c(cols, "u_min", "u_max")
  }
  af <- select(af, all_of(cols))
  return(af)
}

#' Read stochastic genotype sidecar metadata (U_MIN/U_MAX).
#'
#' Reads per-variant min/max values used to invert min-max [0,2] scaling
#' of stochastic genotype data. Supports two formats:
#' \itemize{
#'   \item \strong{afreq}: PLINK2 .afreq/.afreq.zst with U_MIN/U_MAX columns
#'     (read via \code{readAfreq}, which also returns allele frequencies).
#'   \item \strong{generic}: Tab-delimited file with columns id, u_min, u_max.
#' }
#'
#' @param path Path to the sidecar metadata file.
#' @param format One of \code{NULL} (auto-detect from extension), \code{"afreq"},
#'   or \code{"generic"}. When \code{NULL}, files ending in \code{.afreq} or
#'   \code{.afreq.zst} are parsed as afreq; all others as generic.
#' @return A data.frame with columns \code{id}, \code{u_min}, \code{u_max},
#'   or \code{NULL} if the file lacks U_MIN/U_MAX columns (afreq format) or
#'   doesn't exist.
#' @importFrom vroom vroom
#' @noRd
readStochasticMeta <- function(path, format = NULL) {
  if (!file.exists(path)) return(NULL)

  if (is.null(format)) {
    format <- if (grepl("\\.afreq(\\.zst)?$", path)) "afreq" else "generic"
  }
  format <- match.arg(format, c("afreq", "generic"))

  if (format == "afreq") {
    # readAfreq expects a prefix, not a full path - strip the .afreq[.zst] suffix
    prefix <- sub("\\.afreq(\\.zst)?$", "", path)
    af <- readAfreq(prefix)
    if (is.null(af) || !all(c("u_min", "u_max") %in% colnames(af))) return(NULL)
    return(af[, c("id", "u_min", "u_max"), drop = FALSE])
  }

  # Generic: expect tab-delimited with columns id, u_min, u_max
  meta <- as.data.frame(vroom(path, delim = "\t", show_col_types = FALSE))
  required <- c("id", "u_min", "u_max")
  if (!all(required %in% colnames(meta))) {
    stop("Stochastic metadata file '", path, "' must contain columns: ",
         paste(required, collapse = ", "))
  }
  meta[, required, drop = FALSE]
}

#' Search for a stochastic genotype sidecar file alongside a genotype path.
#'
#' Looks for \code{.afreq}, \code{.afreq.zst}, and
#' \code{.stochastic_meta.tsv} files next to the given genotype path.
#' For extension-based paths (VCF, GDS), the extension is stripped first.
#' For prefix-based paths (PLINK1/2), the prefix is used directly.
#'
#' @param genotypePath Path to the genotype data (prefix or file path).
#' @return Path to the first sidecar file found, or \code{NULL}.
#' @noRd
findStochasticMeta <- function(genotypePath) {
  # Strip known genotype extensions to get the stem
  stem <- sub("\\.(vcf|vcf\\.gz|bcf|gds|bed|bim|fam|pgen|pvar|psam)$", "",
              genotypePath)
  candidates <- c(
    paste0(stem, ".afreq"),
    paste0(stem, ".afreq.zst"),
    paste0(stem, ".stochastic_meta.tsv")
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) > 0) found[1] else NULL
}


#' Invert min-max [0,2] scaling to recover the original U matrix.
#'
#' Stochastic genotype data is stored after min-max scaling:
#' U_scaled = 2 * (U - u_min) / (u_max - u_min).
#' This function exactly inverts that transform using the stored per-variant
#' u_min and u_max values from a companion sidecar file (.afreq or
#' .stochastic_meta.tsv).
#'
#' The recovered U satisfies U'U/B ~ Wishart(B, R)/B, the correct distributional
#' property for LD-based fine-mapping with dynamic variance tracking.
#'
#' @param X Numeric matrix (B x p) of min-max scaled values in [0, 2].
#' @param uMin Numeric vector of per-variant minimum values before scaling.
#' @param uMax Numeric vector of per-variant maximum values before scaling.
#' @return Matrix of original U values with same dimensions.
#' @export
invertMinmaxScaling <- function(X, uMin, uMax) {
  if (length(uMin) != ncol(X) || length(uMax) != ncol(X)) {
    stop("Length of u_min/u_max (", length(uMin), ") must equal ncol(X) (", ncol(X), ")")
  }
  denom <- uMax - uMin
  denom[denom == 0] <- 1  # monomorphic: scaling was identity
  # Invert: U_original = U_scaled * (u_max - u_min) / 2 + u_min
  sweep(sweep(X, 2, denom / 2, "*"), 2, uMin, "+")
}

# ---------- Internal helpers for PLINK2 format ----------

#' Resolve and validate PLINK2 file paths for a given prefix.
#' @return Named list with pgen, pvar, psam paths.
#' @noRd
resolvePlink2Paths <- function(prefix) {
  pgen <- paste0(prefix, ".pgen")
  if (!file.exists(pgen)) {
    stop("PLINK2 .pgen file not found at: ", pgen,
         "\n  Note: .pgen must be uncompressed (plink2 does not compress .pgen).")
  }
  # Prefer plain .pvar (fast, no extra deps); fall back to .pvar.zst
  pvar <- if (file.exists(paste0(prefix, ".pvar"))) {
    paste0(prefix, ".pvar")
  } else if (file.exists(paste0(prefix, ".pvar.zst"))) {
    paste0(prefix, ".pvar.zst")
  } else {
    stop("PLINK2 .pvar[.zst] file not found at prefix: ", prefix)
  }
  psam <- paste0(prefix, ".psam")
  if (!file.exists(psam)) {
    stop("PLINK2 .psam file not found at: ", psam,
         "\n  Note: .psam must be uncompressed (plink2 does not compress .psam).")
  }
  list(pgen = pgen, pvar = pvar, psam = psam)
}

#' Read .pvar or .pvar.zst into a data.frame via pgenlibr.
#'
#' Uses pgenlibr::NewPvar() to parse the file (handles both plain .pvar and
#' zstd-compressed .pvar.zst natively, no external CLI required).
#'
#' @param pvarPath Path to .pvar or .pvar.zst file.
#' @return data.frame with columns: chrom, id, pos, A2 (REF), A1 (ALT).
#' @noRd
readPvar <- function(pvarPath) {
  if (!requireNamespace("pgenlibr", quietly = TRUE)) {
    stop("pgenlibr is required. Install from https://cran.r-project.org/web/packages/pgenlibr/index.html")
  }
  pvar <- pgenlibr::NewPvar(pvarPath)
  on.exit(pgenlibr::ClosePvar(pvar), add = TRUE)
  n <- pgenlibr::GetVariantCt(pvar)
  idx <- seq_len(n)
  data.frame(
    chrom = vapply(idx, function(i) pgenlibr::GetVariantChrom(pvar, i), character(1)),
    id    = vapply(idx, function(i) pgenlibr::GetVariantId(pvar, i), character(1)),
    pos   = vapply(idx, function(i) pgenlibr::GetVariantPos(pvar, i), integer(1)),
    A2    = vapply(idx, function(i) pgenlibr::GetAlleleCode(pvar, i, 1L), character(1)),
    A1    = vapply(idx, function(i) pgenlibr::GetAlleleCode(pvar, i, 2L), character(1)),
    stringsAsFactors = FALSE
  )
}

#' Read variant metadata from either .bim or .pvar/.pvar.zst file.
#'
#' Auto-detects the format by extension and header, then returns a
#' standardized data.frame. For PLINK1 .bim files, assigns column names
#' based on the number of columns (6 or 9). For PLINK2 .pvar files,
#' delegates to \code{readPvar()}.
#'
#' @param snpFilePath Path to .bim, .pvar, or .pvar.zst file.
#' @return data.frame with at minimum columns: chrom, id, pos, A2, A1.
#'   Extended .bim files (9 columns) also include: variance, allele_freq, n_nomiss.
#' @importFrom utils read.table
#' @noRd
readVariantMetadata <- function(snpFilePath) {
  isPvar <- grepl("\\.(pvar|pvar\\.zst)$", snpFilePath)
  if (!isPvar) {
    firstLine <- readLines(snpFilePath, n = 1)
    isPvar <- grepl("^#CHROM", firstLine)
  }

  if (isPvar) {
    readPvar(snpFilePath)
  } else {
    df <- read.table(snpFilePath, stringsAsFactors = FALSE)
    n <- ncol(df)
    if (n == 6) {
      names(df) <- c("chrom", "id", "gpos", "pos", "A1", "A2")
    } else if (n == 9) {
      names(df) <- c("chrom", "id", "gpos", "pos", "A1", "A2", "variance", "allele_freq", "n_nomiss")
    } else {
      stop("Unexpected number of columns (", n, ") in variant file: ", snpFilePath)
    }
    df
  }
}

#' Get variant information from any LD reference source.
#'
#' Auto-detects the source type (PLINK2, PLINK1, VCF, GDS, or pre-computed
#' LD metadata) and returns variant metadata. For PLINK2, opens only the
#' .pvar file. For PLINK1, reads only the .bim file. For VCF and GDS,
#' loads the full file and extracts variant info.
#'
#' @param source Genotype file path/prefix or LD metadata file path.
#' @param region Region of interest: "chr:start-end" string or data.frame with
#'   chrom/start/end. If NULL, returns all variants.
#' @return A data.frame with columns: chrom, id, pos, A2, A1.
#'   May also include allele_freq, variance, n_nomiss depending on source.
#'
#' @importFrom vroom vroom
#' @export
getRefVariantInfo <- function(source, region = NULL) {
  resolved <- resolveLdSource(source)

  # For genotype sources via metadata, resolve per-chromosome path
  if (resolved$type %in% c("plink2", "plink1", "vcf", "gds") && !is.null(resolved$metaPath) && !is.null(region)) {
    dataPath <- resolveGenotypePathForRegion(resolved$metaPath, region)
  } else {
    dataPath <- resolved$dataPath
  }

  if (resolved$type == "plink2") {
    paths <- resolvePlink2Paths(dataPath)
    info <- readPvar(paths$pvar)
    afreq <- readAfreq(dataPath)
    if (!is.null(afreq)) {
      info$allele_freq <- afreq$alt_freq[match(info$id, afreq$id)]
    }
  } else if (resolved$type == "plink1") {
    bim <- readBim(paste0(dataPath, ".bed"))
    info <- data.frame(
      chrom = bim$chrom, id = bim$id, pos = bim$pos,
      A2 = bim$a0, A1 = bim$a1,
      stringsAsFactors = FALSE
    )
  } else if (resolved$type %in% c("vcf", "gds")) {
    # VCF/GDS: load via the genotype loader and extract variant_info
    result <- loadGenotypeRegion(dataPath, region = region,
                                 returnVariantInfo = TRUE)
    info <- result$variant_info
    # Compute allele frequency from the genotype matrix
    info$allele_freq <- colMeans(result$X, na.rm = TRUE) / 2
    return(info)  # Already region-filtered by the loader
  } else {
    # Pre-computed LD: read bim/pvar files via metadata
    bimPaths <- getRegionalLdMeta(resolved$metaPath, region)$intersections$bimFilePaths
    info <- do.call(rbind, lapply(bimPaths, function(path) {
      df <- readVariantMetadata(path)
      out <- data.frame(
        chrom = df$chrom, id = df$id, pos = df$pos,
        A2 = df$A2, A1 = df$A1,
        stringsAsFactors = FALSE
      )
      if ("variance" %in% names(df)) out$variance <- df$variance
      if ("allele_freq" %in% names(df)) out$allele_freq <- df$allele_freq
      if ("n_nomiss" %in% names(df)) out$n_nomiss <- df$n_nomiss
      out
    }))
    info$id <- normalizeVariantId(info$id)
    return(info)  # Already region-filtered by getRegionalLdMeta
  }

  # Region filter for plink2/plink1
  if (!is.null(region)) {
    parsed <- parseRegion(region)
    infoChrom <- stripChrPrefix(info$chrom)
    # Handle multi-row region data.frame (one row per chrom)
    if (is.data.frame(parsed) && nrow(parsed) > 1) {
      inRegion <- rep(FALSE, nrow(info))
      for (r in seq_len(nrow(parsed))) {
        inRegion <- inRegion | (infoChrom == as.character(parsed$chrom[r]) &
                                info$pos >= parsed$start[r] & info$pos <= parsed$end[r])
      }
    } else {
      inRegion <- infoChrom == as.character(parsed$chrom) &
                  info$pos >= parsed$start & info$pos <= parsed$end
    }
    info <- info[inRegion, , drop = FALSE]
  }
  info
}

#' Match variant_info against a whitelist file, returning logical index.
#' Uses parse_variant_id() from misc.R to handle all variant ID formats.
#' @importFrom vroom vroom
#' @importFrom readr read_lines
#' @noRd
matchVariantsToKeep <- function(variantInfo, keepVariantsPath) {
  keepRaw <- tryCatch(
    as.data.frame(vroom(keepVariantsPath, show_col_types = FALSE)),
    error = function(e) NULL
  )
  if (!is.null(keepRaw) && "chrom" %in% names(keepRaw) && "pos" %in% names(keepRaw)) {
    keepVariants <- parseVariantId(keepRaw)
  } else {
    # Fall back to reading as single-column variant IDs
    ids <- read_lines(keepVariantsPath)
    keepVariants <- parseVariantId(ids)
  }
  viChrom <- as.integer(stripChrPrefix(variantInfo$chrom))
  hasAlleles <- "A1" %in% names(keepVariants) && "A2" %in% names(keepVariants) &&
    !any(is.na(keepVariants$A1)) && !any(is.na(keepVariants$A2))
  if (hasAlleles) {
    paste0(viChrom, ":", variantInfo$pos, ":", variantInfo$A2, ":", variantInfo$A1) %in%
      paste0(keepVariants$chrom, ":", keepVariants$pos, ":", keepVariants$A2, ":", keepVariants$A1)
  } else {
    paste0(viChrom, ":", variantInfo$pos) %in%
      paste0(keepVariants$chrom, ":", keepVariants$pos)
  }
}

NoSNPsError <- function(message) {
  structure(list(message = message), class = c("NoSNPsError", "error", "condition"))
}




#' Load genotype data for a specific region
#'
#' Auto-detects PLINK2 (.pgen/.pvar[.zst]/.psam), PLINK1 (.bed/.bim/.fam),
#' VCF (.vcf/.vcf.gz/.bcf), or GDS (.gds) format and loads genotype data
#' via \code{\link{readGenotypes}} and \code{\link{extractBlockGenotypes}}.
#' If a stochastic genotype sidecar file (.afreq or
#' .stochastic_meta.tsv) is found alongside the genotype file, non-integer
#' dosages are automatically rescaled using the stored U_MIN/U_MAX values.
#'
#' @param genotype Path to the genotype data file (without extension).
#' @param region The target region in the format "chr:start-end".
#' @param keep_indel Whether to keep indel SNPs.
#' @param keep_variants_path Path to a file listing variants to keep.
#' @param return_variant_info If TRUE, return a list with X (dosage matrix) and
#'   variant_info (data.frame). If FALSE (default), return only the dosage matrix.
#' @param stochastic_meta_path Optional explicit path to a stochastic genotype
#'   sidecar file. If NULL (default), auto-detected via \code{findStochasticMeta}.
#' @param stochastic_meta_format Optional format override for the sidecar file:
#'   \code{"afreq"} or \code{"generic"}. If NULL (default), auto-detected from
#'   file extension.
#' @return If return_variant_info is FALSE, a numeric dosage matrix (rows=samples,
#'   cols=variants). If TRUE, a list with elements X and variant_info.
#'
#' @export
loadGenotypeRegion <- function(genotype, region = NULL, keepIndel = TRUE,
                               keepVariantsPath = NULL,
                               returnVariantInfo = FALSE,
                               stochasticMetaPath = NULL,
                               stochasticMetaFormat = NULL) {
  # --- Detect format and create GenotypeHandle ---
  if (grepl("\\.(vcf|vcf\\.gz|bcf)$", genotype)) {
    handle <- readGenotypes(genotype, format = "vcf")
  } else if (grepl("\\.gds$", genotype)) {
    handle <- readGenotypes(genotype, format = "gds")
  } else if (hasPlink2Files(genotype)) {
    handle <- readGenotypes(genotype, format = "plink2")
  } else if (hasPlink1Files(genotype)) {
    handle <- readGenotypes(genotype, format = "plink1")
  } else {
    stop("Genotype files not found at: ", genotype,
         "\n  Expected: .vcf/.vcf.gz/.bcf, .gds, or PLINK prefix (.pgen/.pvar[.zst]/.psam or .bed/.bim/.fam)")
  }

  # --- Region filter ---
  handleSnpInfo <- getSnpInfo(handle)
  if (!is.null(region)) {
    snpIdx <- .regionToSnpIdx(handleSnpInfo, region)
    if (length(snpIdx) == 0) {
      stop(NoSNPsError(paste("No SNPs found in the specified region", region)))
    }
  } else {
    snpIdx <- seq_len(nrow(handleSnpInfo))
  }

  # --- Extract genotypes (no mean imputation — callers handle missing) ---
  rse <- extractBlockGenotypes(handle, snpIdx, meanImpute = FALSE)
  # Convert RSE to samples x variants matrix for pecotmr convention
  X <- t(assay(rse, "dosage"))
  variantInfo <- .snpInfoToVariantInfo(
    handleSnpInfo[snpIdx, , drop = FALSE])

  # --- Attach allele frequency from .afreq sidecar (plink2 only) ---
  if (getFormat(handle) == "plink2") {
    afreq <- readAfreq(getPath(handle))
    if (!is.null(afreq)) {
      afreqCols <- intersect(c("id", "alt_freq", "obs_ct"), colnames(afreq))
      variantInfo <- merge(variantInfo, afreq[, afreqCols, drop = FALSE],
                           by = "id", all.x = TRUE, sort = FALSE)
    }
  }

  result <- list(X = X, variant_info = variantInfo)

  # --- Post-filters: indels and variant whitelist ---
  if (!keepIndel) {
    snpMask <- isSnpAlleles(result$variant_info$A1, result$variant_info$A2)
    result$X <- result$X[, snpMask, drop = FALSE]
    result$variant_info <- result$variant_info[snpMask, , drop = FALSE]
  }
  if (!is.null(keepVariantsPath)) {
    keepIdx <- matchVariantsToKeep(result$variant_info, keepVariantsPath)
    result$X <- result$X[, keepIdx, drop = FALSE]
    result$variant_info <- result$variant_info[keepIdx, , drop = FALSE]
  }

  # --- Detect and invert stochastic genotype scaling ---
  metaPath <- stochasticMetaPath %||% findStochasticMeta(genotype)
  if (!is.null(metaPath)) {
    smeta <- readStochasticMeta(metaPath, format = stochasticMetaFormat)
    if (!is.null(smeta)) {
      idx <- match(colnames(result$X), smeta$id)
      matched <- !is.na(idx)
      if (any(matched)) {
        result$X[, matched] <- invertMinmaxScaling(
          result$X[, matched, drop = FALSE],
          smeta$u_min[idx[matched]],
          smeta$u_max[idx[matched]]
        )
        result$variant_info$u_min <- smeta$u_min[idx]
        result$variant_info$u_max <- smeta$u_max[idx]
        message("Stochastic genotype detected: restored original scale via ", basename(metaPath))
      }
    }
  } else {
    isStochastic <- !all(result$X == round(result$X), na.rm = TRUE)
    if (isStochastic) {
      warning("Non-integer genotype values detected but no stochastic metadata sidecar found. ",
              "Place a .afreq or .stochastic_meta.tsv file with u_min/u_max columns ",
              "alongside the genotype files to restore the original scale.")
    }
  }

  if (returnVariantInfo) result else result$X
}
