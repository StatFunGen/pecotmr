# =============================================================================
# GenotypeHandle S4 class
# -----------------------------------------------------------------------------
# Lazy file handle for genotype data (PLINK1, PLINK2, VCF/BCF, GDS). Opens
# the file for metadata only (sample IDs and SNP info); dosage extraction
# is deferred until extractBlockGenotypes() is called.
# =============================================================================

#' @include AllGenerics.R
NULL

#' @title Genotype File Handle
#' @description S4 container holding a path + format + metadata for lazy
#'   genotype access. Supports PLINK1 (.bed/.bim/.fam), PLINK2
#'   (.pgen/.pvar/.psam), VCF/BCF, and GDS.
#' @slot path Character, file path.
#' @slot format Character, one of \code{"plink1"}, \code{"plink2"},
#'   \code{"vcf"}, \code{"gds"}.
#' @slot snpInfo data.frame, SNP metadata read from the index/sidecar.
#' @slot nSamples Integer, number of samples.
#' @slot sampleIds Character vector of sample identifiers.
#' @slot pgenPtr Opaque pointer for PLINK2 reader state (NULL otherwise).
#' @export
setClass("GenotypeHandle",
  representation(
    path = "character",
    format = "character",
    snpInfo = "data.frame",
    nSamples = "integer",
    sampleIds = "character",
    pgenPtr = "ANY"
  ),
  validity = function(object) {
    errors <- character()
    if (length(object@path) != 1L)
      errors <- c(errors, "'path' must be a single character string")
    valid_formats <- c("gds", "vcf", "plink1", "plink2")
    if (!object@format %in% valid_formats)
      errors <- c(errors, paste("'format' must be one of:",
                                paste(valid_formats, collapse = ", ")))
    if (length(errors) == 0) TRUE else errors
  }
)

#' @export
setMethod("show", "GenotypeHandle", function(object) {
  cat(sprintf("GenotypeHandle [%s]\n", object@format))
  cat(sprintf("  Path: %s\n", object@path))
  cat(sprintf("  %d samples, %d SNPs\n",
              object@nSamples, nrow(object@snpInfo)))
})

#' @title Create a GenotypeHandle Object
#' @description Construct a \code{GenotypeHandle} from one of several input
#'   forms. Exactly one of the following must be specified:
#'   \describe{
#'     \item{\code{path}}{A single file path with a recognized extension:
#'       \code{.vcf}, \code{.vcf.gz}, \code{.vcf.bgz}, \code{.bcf}, or
#'       \code{.gds}. Format is auto-detected from the extension.}
#'     \item{\code{plink1Prefix}}{A path prefix; the constructor appends
#'       \code{.bed}, \code{.bim}, and \code{.fam} to locate the triplet.}
#'     \item{\code{plink2Prefix}}{A path prefix; the constructor appends
#'       \code{.pgen}, \code{.pvar} (or \code{.pvar.zst}), and \code{.psam}
#'       to locate the triplet.}
#'     \item{\code{bed} + \code{bim} + \code{fam}}{Explicit PLINK1 triplet.
#'       The three files must share a stem; if they don't, use
#'       \code{plink1Prefix} or arrange symlinks at a common stem.}
#'     \item{\code{pgen} + \code{pvar} + \code{psam}}{Explicit PLINK2
#'       triplet. Same shared-stem requirement.}
#'   }
#'   The constructor opens the file for metadata only (sample IDs and SNP
#'   info); dosage extraction is deferred until \code{extractBlockGenotypes()}
#'   is called.
#'
#' @param path Single file path (.vcf/.vcf.gz/.vcf.bgz/.bcf/.gds), or
#'   \code{NULL}.
#' @param plink1Prefix Path prefix for a PLINK1 triplet, or \code{NULL}.
#' @param plink2Prefix Path prefix for a PLINK2 triplet, or \code{NULL}.
#' @param bed,bim,fam Explicit paths to the PLINK1 triplet, or all \code{NULL}.
#' @param pgen,pvar,psam Explicit paths to the PLINK2 triplet, or all
#'   \code{NULL}.
#' @param ldMeta Path to an LD-meta TSV file (columns \code{chrom},
#'   \code{start}, \code{end}, \code{path}; the \code{path} column may be
#'   comma-separated as \code{ld_file,bim_file}). Requires \code{region}.
#'   The constructor resolves the row covering \code{region}, then
#'   delegates to the appropriate file-based handler. When the resolved
#'   row points at PLINK1 / PLINK2 / VCF / GDS files, the corresponding
#'   reader is used; \code{.cor.xz} (pre-computed LD-matrix) rows are not
#'   supported here — use \code{\link{loadLdMatrix}} for that case.
#' @param region Region specification for \code{ldMeta} lookup:
#'   \code{"chr:start-end"} string or a one-row data.frame with
#'   \code{chrom}, \code{start}, \code{end}.
#' @param ... Additional arguments forwarded to the format-specific reader.
#' @return A \code{GenotypeHandle} object.
#' @export
GenotypeHandle <- function(path = NULL,
                           plink1Prefix = NULL, plink2Prefix = NULL,
                           bed = NULL, bim = NULL, fam = NULL,
                           pgen = NULL, pvar = NULL, psam = NULL,
                           ldMeta = NULL, region = NULL,
                           ...) {
  bedTrioGiven <- !is.null(bed) || !is.null(bim) || !is.null(fam)
  bedTrioComplete <- !is.null(bed) && !is.null(bim) && !is.null(fam)
  if (bedTrioGiven && !bedTrioComplete) {
    stop("If specifying the bed/bim/fam triplet, all three must be provided.")
  }
  pgenTrioGiven <- !is.null(pgen) || !is.null(pvar) || !is.null(psam)
  pgenTrioComplete <- !is.null(pgen) && !is.null(pvar) && !is.null(psam)
  if (pgenTrioGiven && !pgenTrioComplete) {
    stop("If specifying the pgen/pvar/psam triplet, all three must be provided.")
  }
  if (!is.null(ldMeta) && is.null(region)) {
    stop("`ldMeta` requires a `region` (a 'chr:start-end' string or a ",
         "one-row data.frame with chrom/start/end).")
  }
  if (is.null(ldMeta) && !is.null(region)) {
    stop("`region` is only meaningful when `ldMeta` is supplied.")
  }

  sources <- c(
    path           = !is.null(path),
    plink1Prefix   = !is.null(plink1Prefix),
    plink2Prefix   = !is.null(plink2Prefix),
    plink1Triplet  = bedTrioComplete,
    plink2Triplet  = pgenTrioComplete,
    ldMeta         = !is.null(ldMeta)
  )
  nSources <- sum(sources)
  if (nSources != 1L) {
    stop("Exactly one of `path`, `plink1Prefix`, `plink2Prefix`, the ",
         "bed/bim/fam triplet, the pgen/pvar/psam triplet, or `ldMeta` ",
         "must be specified (got ", nSources, ").")
  }

  if (sources[["path"]]) {
    return(readGenotypes(path, ...))
  }
  if (sources[["plink1Prefix"]]) {
    return(.makePlink1Handle(plink1Prefix, ...))
  }
  if (sources[["plink2Prefix"]]) {
    return(.makePlink2Handle(plink2Prefix, ...))
  }
  if (sources[["plink1Triplet"]]) {
    return(.genotypeHandleFromPlink1Triplet(bed, bim, fam, ...))
  }
  if (sources[["plink2Triplet"]]) {
    return(.genotypeHandleFromPlink2Triplet(pgen, pvar, psam, ...))
  }
  if (sources[["ldMeta"]]) {
    return(.genotypeHandleFromLdMeta(ldMeta, region, ...))
  }
}

.genotypeHandleFromLdMeta <- function(ldMeta, region, ...) {
  resolved <- getRegionalLdMeta(ldMeta, region)
  ldPaths  <- resolved$intersections$LD_file_paths
  bimPaths <- resolved$intersections$bimFilePaths
  if (length(ldPaths) == 0L) {
    stop("GenotypeHandle: no LD-meta row covers region ", deparse(region),
         " in ", ldMeta, ".")
  }
  if (length(ldPaths) > 1L) {
    stop("GenotypeHandle: region ", deparse(region), " spans multiple LD-meta ",
         "rows; the GenotypeHandle constructor only resolves single-row ",
         "regions. Use loadLdMatrix() for multi-row regions, or restrict the ",
         "region to a single LD block.")
  }
  ldPath  <- ldPaths[[1L]]
  bimPath <- if (length(bimPaths) > 0L) bimPaths[[1L]] else NULL

  if (grepl("\\.cor(\\.xz)?$", ldPath, ignore.case = TRUE)) {
    stop("GenotypeHandle: the LD-meta row for region ", deparse(region),
         " points at a pre-computed correlation matrix (", ldPath,
         "). Use loadLdMatrix() / loadLdSketch() for .cor.xz inputs; ",
         "GenotypeHandle accepts only genotype payloads (VCF/GDS/PLINK).")
  }

  lower <- tolower(ldPath)
  if (grepl("\\.vcf(\\.b?gz)?$", lower) || endsWith(lower, ".bcf")) {
    return(readGenotypes(ldPath, format = "vcf", ...))
  }
  if (endsWith(lower, ".gds")) {
    return(readGenotypes(ldPath, format = "gds", ...))
  }
  if (endsWith(lower, ".bed")) {
    return(.makePlink1Handle(sub("\\.bed$", "", ldPath, ignore.case = TRUE), ...))
  }
  if (endsWith(lower, ".pgen")) {
    return(.makePlink2Handle(sub("\\.pgen$", "", ldPath, ignore.case = TRUE), ...))
  }
  stop("GenotypeHandle: unsupported LD-meta file extension on '", ldPath,
       "'. Expected one of .vcf/.vcf.gz/.vcf.bgz/.bcf/.gds/.bed/.pgen.")
}

.genotypeHandleFromPlink1Triplet <- function(bed, bim, fam, ...) {
  for (f in list(bed = bed, bim = bim, fam = fam)) {
    if (!is.character(f) || length(f) != 1L) {
      stop("Each of `bed`, `bim`, `fam` must be a single file path.")
    }
  }
  stems <- c(
    bed = file_path_sans_ext(bed),
    bim = file_path_sans_ext(bim),
    fam = file_path_sans_ext(fam)
  )
  if (length(unique(stems)) != 1L) {
    stop("`bed`, `bim`, and `fam` must share a common path stem. Got:\n",
         paste0("  ", names(stems), ": ", stems, collapse = "\n"), "\n",
         "If your files are at different paths, either rename them to share ",
         "a stem or arrange symlinks at a common prefix and pass ",
         "`plink1Prefix` instead.")
  }
  .makePlink1Handle(unname(stems[1L]), ...)
}

.genotypeHandleFromPlink2Triplet <- function(pgen, pvar, psam, ...) {
  for (f in list(pgen = pgen, pvar = pvar, psam = psam)) {
    if (!is.character(f) || length(f) != 1L) {
      stop("Each of `pgen`, `pvar`, `psam` must be a single file path.")
    }
  }
  pvarStem <- sub("\\.zst$", "", pvar, ignore.case = TRUE)
  stems <- c(
    pgen = file_path_sans_ext(pgen),
    pvar = file_path_sans_ext(pvarStem),
    psam = file_path_sans_ext(psam)
  )
  if (length(unique(stems)) != 1L) {
    stop("`pgen`, `pvar`, and `psam` must share a common path stem. Got:\n",
         paste0("  ", names(stems), ": ", stems, collapse = "\n"), "\n",
         "If your files are at different paths, either rename them to share ",
         "a stem or arrange symlinks at a common prefix and pass ",
         "`plink2Prefix` instead.")
  }
  .makePlink2Handle(unname(stems[1L]), ...)
}

#' @rdname getSnpInfo
#' @export
setMethod("getSnpInfo", "GenotypeHandle", function(x) x@snpInfo)

#' @rdname getFormat
#' @export
setMethod("getFormat", "GenotypeHandle", function(x) x@format)

#' @rdname getPath
#' @export
setMethod("getPath", "GenotypeHandle", function(x) x@path)

#' @rdname getSampleIds
#' @export
setMethod("getSampleIds", "GenotypeHandle", function(x) x@sampleIds)

#' @rdname getPgenPtr
#' @export
setMethod("getPgenPtr", "GenotypeHandle", function(x) x@pgenPtr)

#' @rdname getNSamples
#' @export
setMethod("getNSamples", "GenotypeHandle", function(x) x@nSamples)
