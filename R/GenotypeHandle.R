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
#' @slot path Character, file path. For a one-file-per-chromosome handle
#'   this is the chrom-meta file path (a display/provenance value); the
#'   per-chromosome payload files live in \code{chromPaths}.
#' @slot format Character, one of \code{"plink1"}, \code{"plink2"},
#'   \code{"vcf"}, \code{"gds"}.
#' @slot snpInfo data.frame, SNP metadata read from the index/sidecar. For a
#'   sharded handle this is the union across chromosomes (row-bound in the
#'   order the shards were supplied), so \code{snpIdx} stays a single global
#'   index space.
#' @slot nSamples Integer, number of samples.
#' @slot sampleIds Character vector of sample identifiers.
#' @slot pgenPtr Opaque pointer for PLINK2 reader state (NULL otherwise).
#' @slot chromPaths Named character vector mapping canonical chromosome
#'   (e.g. \code{"21"}, \code{"X"}) to the per-chromosome payload path/prefix.
#'   Empty (\code{character(0)}) for a single-file handle; non-empty marks a
#'   one-file-per-chromosome (sharded) handle whose extraction is routed by
#'   chromosome.
#' @export
setClass("GenotypeHandle",
  representation(
    path = "character",
    format = "character",
    snpInfo = "data.frame",
    nSamples = "integer",
    sampleIds = "character",
    pgenPtr = "ANY",
    chromPaths = "character"
  ),
  prototype = prototype(
    chromPaths = character(0)
  ),
  validity = function(object) {
    errors <- character()
    if (length(object@path) != 1L)
      errors <- c(errors, "'path' must be a single character string")
    valid_formats <- c("gds", "vcf", "plink1", "plink2")
    if (!object@format %in% valid_formats)
      errors <- c(errors, paste("'format' must be one of:",
                                paste(valid_formats, collapse = ", ")))
    if (length(object@chromPaths) > 0L) {
      nm <- names(object@chromPaths)
      if (is.null(nm) || any(!nzchar(nm)) || anyDuplicated(nm))
        errors <- c(errors, paste("'chromPaths' must be a uniquely-named",
                                  "character vector (names = chromosomes)"))
    }
    if (length(errors) == 0) TRUE else errors
  }
)

#' @export
setMethod("show", "GenotypeHandle", function(object) {
  cat(sprintf("GenotypeHandle [%s]\n", object@format))
  chromPaths <- .genotypeChromPaths(object)
  if (length(chromPaths) > 0L) {
    cat(sprintf("  Chrom-meta: %s\n", object@path))
    cat(sprintf("  %d per-chromosome files: %s\n",
                length(chromPaths), paste(names(chromPaths), collapse = ", ")))
  } else {
    cat(sprintf("  Path: %s\n", object@path))
  }
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
#'     \item{\code{genoMeta}}{One genotype file per chromosome. Either a path
#'       to a whitespace/TSV meta file whose first column is the chromosome
#'       (\code{#chr}) and second column the per-chromosome payload
#'       (a \code{.bed}/\code{.pgen}/\code{.vcf[.gz]}/\code{.bcf}/\code{.gds}
#'       file or a PLINK prefix; relative paths resolve against the meta
#'       file's directory), or a named character vector mapping chromosome to
#'       payload. All shards must share one format and identical sample IDs in
#'       the same order. Extraction is routed to the correct file by the
#'       requested region's chromosome.}
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
#' @param genoMeta One-file-per-chromosome specification: a path to a
#'   \code{#chr,path} meta file or a named character vector
#'   (names = chromosomes, values = payload paths/prefixes). Optionally pass
#'   \code{format} via \code{...} to force a single backend for every shard.
#' @param ... Additional arguments forwarded to the format-specific reader.
#' @return A \code{GenotypeHandle} object.
#' @export
GenotypeHandle <- function(path = NULL,
                           plink1Prefix = NULL, plink2Prefix = NULL,
                           bed = NULL, bim = NULL, fam = NULL,
                           pgen = NULL, pvar = NULL, psam = NULL,
                           ldMeta = NULL, region = NULL,
                           genoMeta = NULL,
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
    ldMeta         = !is.null(ldMeta),
    genoMeta       = !is.null(genoMeta)
  )
  nSources <- sum(sources)
  if (nSources != 1L) {
    stop("Exactly one of `path`, `plink1Prefix`, `plink2Prefix`, the ",
         "bed/bim/fam triplet, the pgen/pvar/psam triplet, `ldMeta`, or ",
         "`genoMeta` must be specified (got ", nSources, ").")
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
  if (sources[["genoMeta"]]) {
    return(.genotypeHandleFromChromMeta(genoMeta, ...))
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

# ---------------------------------------------------------------------------
# One-file-per-chromosome (sharded) handle support
# ---------------------------------------------------------------------------

# Canonicalize a chromosome label to a bare token (strip a leading "chr").
# Used as the routing key in @chromPaths and when matching @snpInfo$CHR.
#' @keywords internal
.canonChr <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)

# Per-chromosome shard map, tolerant of GenotypeHandle objects deserialized
# from before the `chromPaths` slot existed (e.g. an RDS saved by an older
# pecotmr). Such objects have no `chromPaths` slot, so a direct `@` access
# errors; treat them as single-file handles.
#' @keywords internal
.genotypeChromPaths <- function(handle) {
  tryCatch(handle@chromPaths, error = function(e) character(0))
}

# Parse the genoMeta input into a data.frame(chrom, path). Accepts either a
# path to a #chr/path meta file (whitespace- or tab-delimited, with header)
# or a named character vector (names = chromosomes). Relative payload paths
# in a meta file are resolved against the meta file's own directory.
#' @keywords internal
.parseChromMeta <- function(genoMeta) {
  isMetaFile <- is.character(genoMeta) && length(genoMeta) == 1L &&
    is.null(names(genoMeta)) && file.exists(genoMeta)
  if (isMetaFile) {
    meta <- utils::read.table(genoMeta, header = TRUE, sep = "",
                              comment.char = "", stringsAsFactors = FALSE,
                              check.names = FALSE)
    if (ncol(meta) < 2L)
      stop("GenotypeHandle(genoMeta): meta file '", genoMeta,
           "' must have at least 2 columns (chromosome, path).")
    chrom <- as.character(meta[[1L]])
    pth   <- as.character(meta[[2L]])
    base  <- dirname(normalizePath(genoMeta))
    pth <- vapply(pth, function(p) {
      if (grepl("^(/|[A-Za-z]:)", p) || file.exists(p)) p else file.path(base, p)
    }, character(1), USE.NAMES = FALSE)
    return(data.frame(chrom = chrom, path = pth, stringsAsFactors = FALSE))
  }
  if (is.character(genoMeta) && length(genoMeta) >= 1L &&
      !is.null(names(genoMeta))) {
    return(data.frame(chrom = names(genoMeta), path = unname(genoMeta),
                      stringsAsFactors = FALSE))
  }
  stop("GenotypeHandle(genoMeta): expected a path to a `#chr,path` meta file ",
       "or a named character vector (names = chromosomes, values = paths).")
}

# Build a single-file GenotypeHandle for one shard payload, dispatching to the
# right reader. `format` (optional) forces a backend; otherwise it is detected
# from the file extension, falling back to PLINK prefix probing.
#' @keywords internal
.resolveGenotypeShard <- function(p, format = NULL) {
  lower <- tolower(p)
  if (!is.null(format)) {
    if (format == "plink1") return(.makePlink1Handle(p))
    if (format == "plink2") return(.makePlink2Handle(p))
    return(readGenotypes(p, format = format))
  }
  if (grepl("\\.vcf(\\.b?gz)?$", lower) || endsWith(lower, ".bcf"))
    return(readGenotypes(p, format = "vcf"))
  if (endsWith(lower, ".gds")) return(readGenotypes(p, format = "gds"))
  if (endsWith(lower, ".bed"))
    return(.makePlink1Handle(sub("\\.bed$", "", p, ignore.case = TRUE)))
  if (endsWith(lower, ".pgen"))
    return(.makePlink2Handle(sub("\\.pgen$", "", p, ignore.case = TRUE)))
  # No recognized extension: treat as a PLINK prefix, probe for the sidecar.
  if (file.exists(paste0(p, ".bed"))) return(.makePlink1Handle(p))
  if (file.exists(paste0(p, ".pgen"))) return(.makePlink2Handle(p))
  stop("GenotypeHandle(genoMeta): cannot determine genotype format for '", p,
       "'. Use a recognized extension (.bed/.pgen/.vcf[.gz]/.bcf/.gds), a ",
       "PLINK prefix, or pass `format=`.")
}

# Assemble a sharded handle from a per-chromosome meta. Reads each shard's
# metadata via the existing single-file readers, validates a single shared
# format and identical sample IDs (same order, required for cross-shard
# cbind), and row-binds the per-shard snpInfo into one global index space.
#' @keywords internal
.genotypeHandleFromChromMeta <- function(genoMeta, ...) {
  dots   <- list(...)
  format <- dots$format
  parsed <- .parseChromMeta(genoMeta)
  if (nrow(parsed) == 0L)
    stop("GenotypeHandle(genoMeta): no chromosomes found in the meta input.")

  shards <- lapply(parsed$path, .resolveGenotypeShard, format = format)

  formats <- vapply(shards, function(h) h@format, character(1))
  if (length(unique(formats)) != 1L)
    stop("GenotypeHandle(genoMeta): all per-chromosome files must share one ",
         "format; got: ", paste(unique(formats), collapse = ", "), ".")

  sample0 <- shards[[1L]]@sampleIds
  for (i in seq_along(shards)[-1L]) {
    if (!identical(shards[[i]]@sampleIds, sample0))
      stop("GenotypeHandle(genoMeta): all per-chromosome files must have ",
           "identical sample IDs in the same order (mismatch at '",
           parsed$path[[i]], "').")
  }

  unifiedSnpInfo <- do.call(rbind, lapply(shards, function(h) h@snpInfo))
  rownames(unifiedSnpInfo) <- NULL

  chromPaths <- character(0)
  for (i in seq_along(shards)) {
    chs <- unique(.canonChr(shards[[i]]@snpInfo$CHR))
    for (ch in chs) {
      if (ch %in% names(chromPaths))
        stop("GenotypeHandle(genoMeta): chromosome '", ch, "' appears in ",
             "more than one per-chromosome file.")
      chromPaths[[ch]] <- shards[[i]]@path
    }
  }

  metaPath <- if (is.character(genoMeta) && length(genoMeta) == 1L &&
                  is.null(names(genoMeta)) && file.exists(genoMeta))
    normalizePath(genoMeta) else "<chrom-meta>"

  new("GenotypeHandle",
    path       = metaPath,
    format     = formats[[1L]],
    snpInfo    = unifiedSnpInfo,
    nSamples   = shards[[1L]]@nSamples,
    sampleIds  = sample0,
    pgenPtr    = NULL,
    chromPaths = chromPaths
  )
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
