#' @title S4 Method Implementations
#' @description Constructors and accessor method implementations for all
#'   S4 classes: LdData, RegionalData, FineMappingResult, TwasWeights.
#' @name pecotmr-methods
#' @keywords internal
#' @include allGenerics.R
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors DataFrame mcols mcols<-
#' @importFrom GenomicRanges seqnames GRanges
#' @importFrom tools file_path_sans_ext
NULL

# =============================================================================
# GenotypeHandle constructor
# =============================================================================

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
  # Detect partial triplets and reject them early with a helpful message.
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

# Resolve an LD-meta TSV row covering `region`, then construct a
# GenotypeHandle from whichever file it points to (auto-detecting format
# by extension). Errors when the region spans multiple rows (mixture
# panels are handled at the LdData layer via list-valued genotypeHandle;
# see task #32).
# @noRd
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

  # If the row points at a precomputed correlation file (.cor.xz), that's
  # an LD-matrix payload, not a genotype payload — out of scope here.
  if (grepl("\\.cor(\\.xz)?$", ldPath, ignore.case = TRUE)) {
    stop("GenotypeHandle: the LD-meta row for region ", deparse(region),
         " points at a pre-computed correlation matrix (", ldPath,
         "). Use loadLdMatrix() / loadLdSketch() for .cor.xz inputs; ",
         "GenotypeHandle accepts only genotype payloads (VCF/GDS/PLINK).")
  }

  # Auto-detect format from the resolved file (mirrors the dispatch in
  # readGenotypesFromFile / fileUtils.R::readGenotypes auto path).
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
  # If ldPath is a stem and a bim sibling was specified, treat as PLINK1.
  if (!is.null(bimPath) && endsWith(tolower(bimPath), ".bim")) {
    return(.makePlink1Handle(sub("\\.[^.]+$", "", ldPath), ...))
  }
  stop("GenotypeHandle: could not detect genotype format from LD-meta row ",
       "path '", ldPath, "'. Expected .vcf(.gz|.bgz)/.bcf/.gds/.bed/.pgen.")
}

# Validate that three explicit PLINK1 paths share a common stem and delegate
# to the standard prefix-based handle constructor.
#' @keywords internal
#' @noRd
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

# Validate that three explicit PLINK2 paths share a common stem and delegate
# to the standard prefix-based handle constructor. Accepts .pvar or .pvar.zst.
#' @keywords internal
#' @noRd
.genotypeHandleFromPlink2Triplet <- function(pgen, pvar, psam, ...) {
  for (f in list(pgen = pgen, pvar = pvar, psam = psam)) {
    if (!is.character(f) || length(f) != 1L) {
      stop("Each of `pgen`, `pvar`, `psam` must be a single file path.")
    }
  }
  # .pvar may be .pvar.zst, so strip both extensions if present
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

# =============================================================================
# LdData constructor and accessors
# =============================================================================

#' @title Create an LdData Object
#' @description Construct an \code{LdData} from a correlation matrix and/or
#'   genotype handle, plus variant metadata as a GRanges.
#' @param correlation A correlation matrix, list of matrices, or NULL.
#' @param genotypeHandle A GenotypeHandle, list of GenotypeHandles, or NULL.
#' @param snpIdx Integer vector of SNP indices, or NULL.
#' @param variants A GRanges with variant metadata (must have variant_id in
#'   mcols, plus A1, A2).
#' @param blockMetadata LdBlocks or data.frame with block info.
#' @param nRef Integer, reference panel sample size.
#' @param mixtureWeights Optional numeric vector of mixing proportions,
#'   one per panel in \code{genotypeHandle} when it is a list. Must be
#'   non-negative and sum to 1. Required whenever
#'   \code{genotypeHandle} is a list and downstream code will call
#'   \code{getCorrelation()}.
#' @return An \code{LdData} object.
#' @export
LdData <- function(correlation = NULL, genotypeHandle = NULL,
                   snpIdx = NULL, variants, blockMetadata,
                   nRef = 0L, mixtureWeights = NULL) {
  obj <- new("LdData",
    correlation = correlation,
    genotypeHandle = genotypeHandle,
    snpIdx = snpIdx,
    variants = variants,
    blockMetadata = blockMetadata,
    nRef = as.integer(nRef),
    mixtureWeights = mixtureWeights
  )
  validObject(obj)
  obj
}

#' @rdname getCorrelation
#' @export
setMethod("getCorrelation", "LdData", function(x) {
  if (!is.null(x@correlation)) return(x@correlation)
  if (is.null(x@genotypeHandle)) {
    stop("No correlation matrix or genotype handle available")
  }
  if (is.list(x@genotypeHandle)) {
    if (is.null(x@mixtureWeights))
      stop("Cannot compute mixture LD: `mixtureWeights` is NULL. ",
           "Construct LdData with mixtureWeights = <numeric vector> ",
           "when supplying a list of GenotypeHandles.")
    perPanel <- lapply(x@genotypeHandle, function(h) {
      geno <- extractBlockGenotypes(h, x@snpIdx)
      X <- t(assay(geno, "dosage"))
      computeLd(X, method = "sample")
    })
    # All panels must agree on dimension; weighted sum element-wise.
    dims <- vapply(perPanel, function(R) nrow(R), integer(1))
    if (length(unique(dims)) != 1L)
      stop("Mixture panels yielded LD matrices of differing dimensions: ",
           paste(dims, collapse = ", "),
           ". All panels must be aligned on the same variant subset.")
    w <- x@mixtureWeights
    R <- matrix(0, nrow = dims[[1L]], ncol = dims[[1L]])
    for (k in seq_along(perPanel)) R <- R + w[[k]] * perPanel[[k]]
    dimnames(R) <- dimnames(perPanel[[1L]])
    return(R)
  }
  geno <- extractBlockGenotypes(x@genotypeHandle, x@snpIdx)
  X <- t(assay(geno, "dosage"))
  computeLd(X, method = "sample")
})

#' @rdname getGenotypes
#' @export
setMethod("getGenotypes", "LdData", function(x, ...) {
  if (is.null(x@genotypeHandle)) return(NULL)
  # Plain matrix stored directly (e.g. from loadLdSketch after filtering)
  if (is.matrix(x@genotypeHandle)) return(x@genotypeHandle)
  if (is.list(x@genotypeHandle)) {
    lapply(x@genotypeHandle, function(h) {
      geno <- extractBlockGenotypes(h, x@snpIdx)
      t(assay(geno, "dosage"))
    })
  } else {
    geno <- extractBlockGenotypes(x@genotypeHandle, x@snpIdx)
    t(assay(geno, "dosage"))
  }
})

#' @rdname hasGenotypes
#' @export
setMethod("hasGenotypes", "LdData", function(x) {
  !is.null(x@genotypeHandle)
})

#' @rdname getVariantIds
#' @export
setMethod("getVariantIds", "LdData", function(x, ...) {
  mcols(x@variants)$variant_id
})

#' @rdname getVariantInfo
#' @export
setMethod("getVariantInfo", "LdData", function(x) {
  x@variants
})

#' @rdname getBlockMetadata
#' @export
setMethod("getBlockMetadata", "LdData", function(x) {
  x@blockMetadata
})

#' @rdname getRefPanel
#' @export
setMethod("getRefPanel", "LdData", function(x) {
  mc <- as.data.frame(mcols(x@variants))
  mc$chrom <- as.character(seqnames(x@variants))
  mc$pos <- start(x@variants)
  mc
})

# =============================================================================
# Helper: build variant GRanges from refPanel data.frame
# =============================================================================

#' @title Build Variant GRanges
#' @description Convert a refPanel data.frame to a GRanges object suitable
#'   for the LdData variants slot.
#' @param refPanel data.frame with columns: chrom, pos, A1, A2, variant_id.
#'   May also have allele_freq, variance, n_nomiss.
#' @return A GRanges object.
#' @keywords internal
#' @noRd
.refPanelToGranges <- function(refPanel) {
  chr <- as.character(refPanel$chrom)
  chr <- sub("^chr", "", chr, ignore.case = TRUE)
  chr <- paste0("chr", chr)
  pos <- as.integer(refPanel$pos)

  gr <- GRanges(
    seqnames = chr,
    ranges = IRanges(start = pos, width = 1L)
  )

  mcolsData <- DataFrame(
    variant_id = refPanel$variant_id,
    A1 = refPanel$A1,
    A2 = refPanel$A2
  )

  optional <- c("allele_freq", "variance", "n_nomiss")
  for (col in optional) {
    if (col %in% names(refPanel)) {
      mcolsData[[col]] <- refPanel[[col]]
    }
  }
  mcols(gr) <- mcolsData
  gr
}

# =============================================================================
# FineMappingResult collection constructor and accessors
# =============================================================================

#' @title Create a QtlFineMappingResult Collection
#' @description Construct a \code{QtlFineMappingResult} DFrame-subclass
#'   collection from per-tuple vectors and a list of
#'   \code{FineMappingEntry} payloads (one per tuple). The optional
#'   \code{ldSketch} slot records the LD reference used for RSS-derived
#'   fits; pass \code{NULL} (the default) for individual-level fits.
#' @param study Character vector of study identifiers (per tuple).
#' @param context Character vector of context labels (per tuple).
#' @param trait Character vector of trait identifiers (per tuple).
#' @param method Character vector of fine-mapping method names (per tuple).
#' @param entry List / \code{SimpleList} of \code{FineMappingEntry} objects.
#' @param ldSketch An optional \code{GenotypeHandle} (the LD reference for
#'   RSS-derived fits), or \code{NULL} for individual-level fits.
#' @return A \code{QtlFineMappingResult} object.
#' @export
QtlFineMappingResult <- function(study, context, trait, method, entry,
                                 ldSketch = NULL) {
  n <- length(study)
  if (length(context) != n || length(trait) != n || length(method) != n ||
      length(entry) != n) {
    stop("`study`, `context`, `trait`, `method`, and `entry` must all ",
         "have the same length.")
  }
  cols <- list(
    study   = as.character(study),
    context = as.character(context),
    trait   = as.character(trait),
    method  = as.character(method),
    entry   = S4Vectors::SimpleList(entry)
  )
  df <- do.call(S4Vectors::DataFrame,
                c(cols, list(check.names = FALSE)))
  obj <- new("QtlFineMappingResult", df, ldSketch = ldSketch)
  validObject(obj)
  obj
}

#' @title Create a GwasFineMappingResult Collection
#' @description Construct a \code{GwasFineMappingResult} DFrame-subclass
#'   collection from per-(study, method) tuples and a list of
#'   \code{FineMappingEntry} payloads. The collection represents one LD
#'   block of GWAS fine-mapping fits; build a separate collection per
#'   block when sweeping the genome.
#' @param study Character vector of study identifiers (per tuple).
#' @param method Character vector of fine-mapping method names (per tuple).
#' @param entry List / \code{SimpleList} of \code{FineMappingEntry} objects.
#' @param ldSketch An optional \code{GenotypeHandle}.
#' @return A \code{GwasFineMappingResult} object.
#' @export
GwasFineMappingResult <- function(study, method, entry,
                                  ldSketch = NULL) {
  n <- length(study)
  if (length(method) != n || length(entry) != n) {
    stop("`study`, `method`, and `entry` must all have the same length.")
  }
  cols <- list(
    study  = as.character(study),
    method = as.character(method),
    entry  = S4Vectors::SimpleList(entry)
  )
  df <- do.call(S4Vectors::DataFrame,
                c(cols, list(check.names = FALSE)))
  obj <- new("GwasFineMappingResult", df, ldSketch = ldSketch)
  validObject(obj)
  obj
}

# Internal: return integer row indices of `x` where every (column, value)
# pair in `keys` matches as.character(x[[column]]) == value. Shared
# building block for all of pecotmr's tuple-keyed row selectors and
# cache lookups (.tupleSelectRow, .qtlSumStatsSelectRow,
# .gwasSelectStudy, .fmCacheLookup, .cipFmrHasTuple, etc.). Pure
# vectorised AND-match with character coercion -- no validation, no
# error reporting.
.matchTupleRows <- function(x, keys) {
  if (length(keys) == 0L) return(seq_len(nrow(x)))
  ok <- rep(TRUE, nrow(x))
  for (k in names(keys)) {
    ok <- ok & as.character(x[[k]]) == keys[[k]]
  }
  which(ok)
}

# Internal: resolve a tuple-keyed selection (study, context, trait,
# method) to a single row index. Used by the QtlFineMappingResult and
# TwasWeights accessors. Returns an error when no row matches; returns
# the single row index when the collection has exactly one row and any
# selector argument was omitted.
.tupleSelectRow <- function(x, study, context, trait, method,
                            cls = "QtlFineMappingResult") {
  if (nrow(x) == 0L) stop(cls, " has no rows.")
  if (missing(study) || is.null(study) ||
      missing(context) || is.null(context) ||
      missing(trait) || is.null(trait) ||
      missing(method) || is.null(method)) {
    if (nrow(x) == 1L) return(1L)
    stop(cls, " has ", nrow(x), " entries. Pass `study`, `context`, ",
         "`trait`, and `method` to select one.")
  }
  if (length(study) != 1L || length(context) != 1L ||
      length(trait) != 1L || length(method) != 1L) {
    stop("`study`, `context`, `trait`, and `method` must each be length 1.")
  }
  idx <- .matchTupleRows(x, list(study = study, context = context,
                                  trait = trait, method = method))
  if (length(idx) == 0L) {
    stop(sprintf(
      "No entry for (study='%s', context='%s', trait='%s', method='%s').",
      study, context, trait, method))
  }
  idx[[1L]]
}

# Internal: resolve a (study, method) tuple to a single row index of a
# GwasFineMappingResult collection.
.tupleSelectRowGwasFmr <- function(x, study, method) {
  if (nrow(x) == 0L) stop("GwasFineMappingResult has no rows.")
  if (missing(study) || is.null(study) ||
      missing(method) || is.null(method)) {
    if (nrow(x) == 1L) return(1L)
    stop("GwasFineMappingResult has ", nrow(x), " entries. Pass `study` ",
         "and `method` to select one.")
  }
  if (length(study) != 1L || length(method) != 1L)
    stop("`study` and `method` must each be length 1.")
  idx <- .matchTupleRows(x, list(study = study, method = method))
  if (length(idx) == 0L)
    stop(sprintf("No entry for (study='%s', method='%s').", study, method))
  idx[[1L]]
}

#' @title Get a Single Fine-Mapping Entry
#' @description Return the \code{FineMappingEntry} for one
#'   \code{(study, context, trait, method)} row of a
#'   \code{FineMappingResult} collection.
#' @param x A \code{FineMappingResult} object.
#' @param study,context,trait,method Single character identifiers. All
#'   required when the collection has more than one row; optional when
#'   the collection has a single row.
#' @return A \code{FineMappingEntry} object.
#' @export
setGeneric("getFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL)
    standardGeneric("getFineMappingResult"))

#' @rdname getFineMappingResult
#' @export
setMethod("getFineMappingResult", "QtlFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL) {
    idx <- .tupleSelectRow(x, study, context, trait, method,
                           cls = "QtlFineMappingResult")
    x$entry[[idx]]
  })

# Derived collection-level accessors (delegate to entry-level methods).

#' @rdname getPip
#' @export
setMethod("getPip", "QtlFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL,
           returnList = FALSE, ...) {
    entry <- getFineMappingResult(x, study, context, trait, method)
    pip <- getPip(entry)
    if (isTRUE(returnList)) {
      nm <- sprintf("%s|%s|%s|%s",
                    as.character(x$study)[1L],
                    as.character(x$context)[1L],
                    as.character(x$trait)[1L],
                    as.character(x$method)[1L])
      out <- list(); out[[nm]] <- pip
      return(out)
    }
    pip
  })

#' @rdname getCs
#' @export
setMethod("getCs", "QtlFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL,
           coverage = 0.95, ...) {
    entry <- getFineMappingResult(x, study, context, trait, method)
    getCs(entry, coverage = coverage)
  })

#' @rdname getTopLoci
#' @export
setMethod("getTopLoci", "QtlFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL,
           ...) {
    entry <- getFineMappingResult(x, study, context, trait, method)
    getTopLoci(entry)
  })

#' @rdname getTrimmedFit
#' @export
setMethod("getTrimmedFit", "QtlFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL,
           ...) {
    entry <- getFineMappingResult(x, study, context, trait, method)
    getTrimmedFit(entry)
  })

#' @rdname getVariantIds
#' @export
setMethod("getVariantIds", "QtlFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL,
           ...) {
    entry <- getFineMappingResult(x, study, context, trait, method)
    getVariantIds(entry)
  })

#' @rdname getStudy
#' @export
setMethod("getStudy", "FineMappingResultBase",
          function(x) unique(as.character(x$study)))

#' @rdname getLdSketch
#' @export
setMethod("getLdSketch", "FineMappingResultBase",
          function(x, ...) x@ldSketch)

#' @rdname getMethodNames
#' @export
setMethod("getMethodNames", "FineMappingResultBase",
          function(x) unique(as.character(x$method)))

#' @rdname getContexts
#' @export
setMethod("getContexts", "QtlFineMappingResult",
          function(x) unique(as.character(x$context)))

#' @rdname getTraits
#' @export
setMethod("getTraits", "QtlFineMappingResult",
          function(x) unique(as.character(x$trait)))

# GwasFineMappingResult has no context / trait columns; the generic
# returns NULL so callers can write generic code that handles either
# class without conditionals.
#' @rdname getContexts
#' @export
setMethod("getContexts", "GwasFineMappingResult", function(x) NULL)

#' @rdname getTraits
#' @export
setMethod("getTraits", "GwasFineMappingResult", function(x) NULL)

# Per-tuple lookup for the GWAS variant (2-tuple instead of 4-tuple).
# The generic accepts the full set of selectors; context/trait args are
# ignored for GwasFineMappingResult.
#' @rdname getFineMappingResult
#' @export
setMethod("getFineMappingResult", "GwasFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL) {
    idx <- .tupleSelectRowGwasFmr(x, study, method)
    x$entry[[idx]]
  })

#' @rdname getPip
#' @export
setMethod("getPip", "GwasFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL,
           returnList = FALSE, ...) {
    entry <- getFineMappingResult(x, study = study, method = method)
    getPip(entry)
  })

#' @rdname getCs
#' @export
setMethod("getCs", "GwasFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL, ...) {
    entry <- getFineMappingResult(x, study = study, method = method)
    getCs(entry)
  })

#' @rdname getTopLoci
#' @export
setMethod("getTopLoci", "GwasFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL, ...) {
    entry <- getFineMappingResult(x, study = study, method = method)
    getTopLoci(entry)
  })

#' @rdname getTrimmedFit
#' @export
setMethod("getTrimmedFit", "GwasFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL, ...) {
    entry <- getFineMappingResult(x, study = study, method = method)
    getTrimmedFit(entry)
  })

#' @rdname getVariantIds
#' @export
setMethod("getVariantIds", "GwasFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL, ...) {
    entry <- getFineMappingResult(x, study = study, method = method)
    getVariantIds(entry)
  })

# =============================================================================
# SumStatsBase shared accessors
# =============================================================================
#
# Methods defined once on the virtual SumStatsBase parent so QtlSumStats
# and GwasSumStats share the implementation. Class-specific accessors
# (getZ / getN / getMaf / nSnps / subsetChr / getVarY / getSumStats) stay
# on the concrete subclass because they rely on the tuple shape (3-tuple
# for QtlSumStats, 1-tuple for GwasSumStats).

#' @rdname getGenome
#' @export
setMethod("getGenome", "SumStatsBase", function(x, ...) x@genome)

#' @rdname getQcInfo
#' @export
setMethod("getQcInfo", "SumStatsBase", function(x, ...) x@qcInfo)

#' @rdname getLdSketch
#' @export
setMethod("getLdSketch", "SumStatsBase", function(x, ...) x@ldSketch)

#' @rdname getStudy
#' @export
setMethod("getStudy", "SumStatsBase",
          function(x) unique(as.character(x$study)))

# =============================================================================
# TwasWeights collection constructor and accessors
# =============================================================================

#' @title Create a TwasWeights Collection Object
#' @description Construct a \code{TwasWeights} DFrame-subclass collection
#'   from per-tuple vectors and a list of \code{TwasWeightsEntry}
#'   payloads (one per tuple).
#' @param study Character vector of study identifiers.
#' @param context Character vector of context labels.
#' @param trait Character vector of trait identifiers.
#' @param method Character vector of TWAS weight method names.
#' @param entry List / \code{SimpleList} of \code{TwasWeightsEntry} objects.
#' @param ldSketch An optional \code{GenotypeHandle}, or \code{NULL} for
#'   individual-level fits.
#' @return A \code{TwasWeights} object.
#' @export
TwasWeights <- function(study, context, trait, method, entry,
                        ldSketch = NULL) {
  n <- length(study)
  if (length(context) != n || length(trait) != n || length(method) != n ||
      length(entry) != n) {
    stop("`study`, `context`, `trait`, `method`, and `entry` must all ",
         "have the same length.")
  }
  cols <- list(
    study   = as.character(study),
    context = as.character(context),
    trait   = as.character(trait),
    method  = as.character(method),
    entry   = S4Vectors::SimpleList(entry)
  )
  df <- do.call(S4Vectors::DataFrame,
                c(cols, list(check.names = FALSE)))
  obj <- new("TwasWeights", df, ldSketch = ldSketch)
  validObject(obj)
  obj
}

#' @title Get a Single TWAS Weights Entry
#' @description Return the \code{TwasWeightsEntry} for one
#'   \code{(study, context, trait, method)} row of a \code{TwasWeights}
#'   collection.
#' @param x A \code{TwasWeights} object.
#' @param study,context,trait,method Single character identifiers. All
#'   required when the collection has more than one row; optional when
#'   the collection has a single row.
#' @return A \code{TwasWeightsEntry} object.
#' @export
setGeneric("getTwasWeights",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL)
    standardGeneric("getTwasWeights"))

#' @rdname getTwasWeights
#' @export
setMethod("getTwasWeights", "TwasWeights",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL) {
    idx <- .tupleSelectRow(x, study, context, trait, method,
                           cls = "TwasWeights")
    x$entry[[idx]]
  })

#' @rdname getWeights
#' @export
setMethod("getWeights", "TwasWeights",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL,
           ...) {
    entry <- getTwasWeights(x, study, context, trait, method)
    getWeights(entry)
  })

#' @rdname getCvPerformance
#' @export
setMethod("getCvPerformance", "TwasWeights",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL,
           ...) {
    entry <- getTwasWeights(x, study, context, trait, method)
    getCvPerformance(entry)
  })

#' @rdname getFits
#' @export
setMethod("getFits", "TwasWeights",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL,
           ...) {
    entry <- getTwasWeights(x, study, context, trait, method)
    getFits(entry)
  })

#' @rdname getStandardized
#' @export
setMethod("getStandardized", "TwasWeights",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL) {
    entry <- getTwasWeights(x, study, context, trait, method)
    getStandardized(entry)
  })

#' @rdname getDataType
#' @export
setMethod("getDataType", "TwasWeights",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL) {
    entry <- getTwasWeights(x, study, context, trait, method)
    getDataType(entry)
  })

#' @rdname getVariantIds
#' @export
setMethod("getVariantIds", "TwasWeights",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL,
           ...) {
    entry <- getTwasWeights(x, study, context, trait, method)
    getVariantIds(entry)
  })

#' @rdname getStudy
#' @export
setMethod("getStudy", "TwasWeights",
          function(x) unique(as.character(x$study)))

#' @rdname getLdSketch
#' @export
setMethod("getLdSketch", "TwasWeights",
          function(x, ...) x@ldSketch)

#' @rdname getContexts
#' @export
setMethod("getContexts", "TwasWeights",
          function(x) unique(as.character(x$context)))

#' @rdname getTraits
#' @export
setMethod("getTraits", "TwasWeights",
          function(x) unique(as.character(x$trait)))

#' @rdname getMethodNames
#' @export
setMethod("getMethodNames", "TwasWeights",
          function(x) unique(as.character(x$method)))


# =============================================================================
# QtlDataset constructor and accessors
# =============================================================================

#' @title Create a QtlDataset Object
#' @description Construct a \code{QtlDataset} S4 object containing one
#'   study's individual-level QTL data: a genotype handle and a named list
#'   of \code{SummarizedExperiment} objects (one per QTL context), plus
#'   genotype-derived covariates and a residual-scaling flag.
#' @param study Character (length 1). Study identifier.
#' @param genotypes A \code{GenotypeHandle}.
#' @param phenotypes Named list of \code{SummarizedExperiment} objects,
#'   keyed by context. Each SE must have \code{rowRanges} carrying trait
#'   positions and \code{colData} carrying per-context phenotype covariates.
#' @param genotypeCovariates Numeric matrix of genotype-derived covariates
#'   (e.g., ancestry PCs); rows are samples.
#' @param scaleResiduals Logical (length 1). Default \code{TRUE}.
#' @return A \code{QtlDataset} object.
#' @export
QtlDataset <- function(study, genotypes, phenotypes,
                       genotypeCovariates = matrix(numeric(0), nrow = 0, ncol = 0),
                       scaleResiduals = TRUE,
                       mafCutoff = 0,
                       macCutoff = 0,
                       xvarCutoff = 0,
                       imissCutoff = 0,
                       keepSamples = character(0),
                       keepVariants = character(0)) {
  obj <- new("QtlDataset",
             study              = as.character(study),
             genotypes          = genotypes,
             phenotypes         = phenotypes,
             genotypeCovariates = as.matrix(genotypeCovariates),
             scaleResiduals     = isTRUE(scaleResiduals),
             mafCutoff          = as.numeric(mafCutoff),
             macCutoff          = as.numeric(macCutoff),
             xvarCutoff         = as.numeric(xvarCutoff),
             imissCutoff        = as.numeric(imissCutoff),
             keepSamples        = as.character(keepSamples),
             keepVariants       = as.character(keepVariants))
  validObject(obj)
  obj
}

#' @rdname getStudy
#' @export
setMethod("getStudy", "QtlDataset", function(x) x@study)

#' @rdname getContexts
#' @export
setMethod("getContexts", "QtlDataset", function(x) names(x@phenotypes))

#' @rdname getGenotypeCovariates
#' @export
setMethod("getGenotypeCovariates", "QtlDataset",
          function(x) x@genotypeCovariates)

#' @rdname getScaleResiduals
#' @export
setMethod("getScaleResiduals", "QtlDataset", function(x) x@scaleResiduals)

# --- Internal: resolve the variant-selection region for the genotype handle.
# Returns a single GRanges. When `traitId` is supplied, expand each trait's
# rowRange by `cisWindow` bp and take the union span (per the multi-trait rule:
# `[min(start) - cisWindow, max(end) + cisWindow]`). When `region` is supplied,
# extend by `cisWindow` if given. Exactly one of (traitId, region) may be
# supplied; if neither is, return NULL meaning "all variants in handle".
.qtlResolveVariantRegion <- function(x, traitId = NULL, region = NULL,
                                     cisWindow = NULL) {
  if (!is.null(traitId) && !is.null(region)) {
    stop("Specify either `traitId` or `region`, not both.")
  }
  if (is.null(traitId) && is.null(region)) {
    return(NULL)
  }
  if (!is.null(traitId)) {
    if (is.null(cisWindow) || length(cisWindow) != 1L || cisWindow < 0) {
      stop("`cisWindow` is required (and must be non-negative) when ",
           "`traitId` is specified.")
    }
    # Build a GRanges from each trait's rowRanges across all contexts;
    # take the union span +/- cisWindow.
    perTraitRanges <- list()
    for (ctxIdx in seq_along(x@phenotypes)) {
      se <- x@phenotypes[[ctxIdx]]
      rr <- SummarizedExperiment::rowRanges(se)
      hits <- match(traitId, rownames(se))
      hits <- hits[!is.na(hits)]
      if (length(hits) > 0) {
        perTraitRanges[[length(perTraitRanges) + 1L]] <- rr[hits]
      }
    }
    if (length(perTraitRanges) == 0L) {
      stop("None of the requested traitId values were found in any context.")
    }
    allRanges <- do.call(c, perTraitRanges)
    chrs <- unique(as.character(GenomicRanges::seqnames(allRanges)))
    if (length(chrs) != 1L) {
      stop("Multi-trait variant extraction requires all selected traits to ",
           "share a chromosome (got: ",
           paste(chrs, collapse = ", "), ").")
    }
    spanStart <- max(1L, min(GenomicRanges::start(allRanges)) - cisWindow)
    spanEnd   <- max(GenomicRanges::end(allRanges)) + cisWindow
    return(GenomicRanges::GRanges(
      seqnames = chrs,
      ranges   = IRanges::IRanges(start = spanStart, end = spanEnd)
    ))
  }
  # region path
  if (!methods::is(region, "GRanges")) {
    stop("`region` must be a GRanges object.")
  }
  if (length(region) != 1L) {
    stop("`region` must be a single range.")
  }
  if (!is.null(cisWindow)) {
    if (length(cisWindow) != 1L || cisWindow < 0) {
      stop("`cisWindow` must be a single non-negative value.")
    }
    region <- GenomicRanges::GRanges(
      seqnames = GenomicRanges::seqnames(region),
      ranges   = IRanges::IRanges(
        start = max(1L, GenomicRanges::start(region) - cisWindow),
        end   = GenomicRanges::end(region) + cisWindow
      )
    )
  }
  region
}

# Internal: map a GRanges region into 1-based snpIdx into handle@snpInfo.
.qtlVariantIndices <- function(x, region = NULL) {
  handle <- x@genotypes
  if (is.null(region)) {
    return(seq_len(nrow(handle@snpInfo)))
  }
  chr <- as.character(GenomicRanges::seqnames(region))[[1L]]
  chrCanon <- sub("^chr", "", chr, ignore.case = TRUE)
  snpInfo <- handle@snpInfo
  siChr <- sub("^chr", "", as.character(snpInfo$CHR), ignore.case = TRUE)
  bp <- as.integer(snpInfo$BP)
  start <- GenomicRanges::start(region)
  end   <- GenomicRanges::end(region)
  which(siChr == chrCanon & bp >= start & bp <= end)
}

# Internal: extract the panel dosage block (samples x variants) for the
# requested region, narrow to the requested sample set, and apply lazy QC
# (per-sample imiss filter, then per-variant max(mafCutoff,
# macCutoff / (2 * n)) and xvarCutoff filters). Used by getGenotypes,
# getResidualizedGenotypes (via getGenotypes), and getMaf so all three
# share a single variant/sample selection result.
#
# Returns a list:
#   geno       : numeric matrix (kept samples x kept variants)
#   variantIds : character vector of kept variant IDs (= colnames(geno))
#   sampleIds  : character vector of kept sample IDs (= rownames(geno))
#   maf        : numeric vector of per-variant MAF for kept variants
.qtlExtractBlock <- function(x, traitId = NULL, region = NULL,
                             cisWindow = NULL, samples = NULL) {
  gr <- .qtlResolveVariantRegion(x, traitId = traitId, region = region,
                                 cisWindow = cisWindow)
  snpIdx <- .qtlVariantIndices(x, gr)
  if (length(snpIdx) == 0L) {
    return(list(
      geno       = matrix(numeric(0), nrow = x@genotypes@nSamples, ncol = 0L,
                          dimnames = list(x@genotypes@sampleIds, character(0))),
      variantIds = character(0),
      sampleIds  = x@genotypes@sampleIds,
      maf        = numeric(0)
    ))
  }

  # Apply keepVariants restriction before materialization so we do not
  # extract dosage we will immediately drop.
  if (length(x@keepVariants) > 0L) {
    snpAll <- as.character(x@genotypes@snpInfo$SNP[snpIdx])
    keepMask <- snpAll %in% x@keepVariants
    snpIdx <- snpIdx[keepMask]
    if (length(snpIdx) == 0L) {
      return(list(
        geno       = matrix(numeric(0), nrow = 0L, ncol = 0L,
                            dimnames = list(character(0), character(0))),
        variantIds = character(0),
        sampleIds  = character(0),
        maf        = numeric(0)
      ))
    }
  }

  block <- extractBlockGenotypes(x@genotypes, snpIdx, meanImpute = FALSE)
  # `block` is variants x samples (Bioc convention); transpose to
  # samples x variants for analysis-style operations.
  dosage <- t(SummarizedExperiment::assay(block, "dosage"))

  # Resolve the requested sample set: keepSamples (panel-level) intersected
  # with the per-call samples arg, then intersected with the panel sample IDs.
  panelSamples <- rownames(dosage)
  keep <- panelSamples
  if (length(x@keepSamples) > 0L) {
    keep <- intersect(keep, x@keepSamples)
  }
  if (!is.null(samples)) {
    keep <- intersect(keep, as.character(samples))
  }
  if (length(keep) == 0L) {
    return(list(
      geno       = dosage[integer(0), , drop = FALSE],
      variantIds = colnames(dosage),
      sampleIds  = character(0),
      maf        = rep(NA_real_, ncol(dosage))
    ))
  }
  dosage <- dosage[keep, , drop = FALSE]

  # Per-sample missingness filter.
  if (x@imissCutoff > 0 && nrow(dosage) > 0L && ncol(dosage) > 0L) {
    imiss <- rowMeans(is.na(dosage))
    keepSampleMask <- imiss <= x@imissCutoff
    dosage <- dosage[keepSampleMask, , drop = FALSE]
  }

  # Per-variant MAF / MAC / X-variance filters against the post-narrowing
  # sample count. We mean-impute internally for the variance / dosage
  # returned but compute MAF from the un-imputed values so missingness is
  # handled correctly.
  nSamp <- nrow(dosage)
  if (ncol(dosage) > 0L) {
    nObs <- colSums(!is.na(dosage))
    sumD <- colSums(dosage, na.rm = TRUE)
    p <- ifelse(nObs > 0L, sumD / (2 * nObs), NA_real_)
    mafVec <- pmin(p, 1 - p)
    effectiveMaf <- max(x@mafCutoff, if (nSamp > 0L)
      x@macCutoff / (2 * nSamp) else 0)
    keepVarMask <- !is.na(mafVec) & mafVec >= effectiveMaf
    if (x@xvarCutoff > 0 && nSamp > 1L) {
      # Compute variance with mean imputation per column so variance is
      # defined when missingness is present.
      mu <- ifelse(nObs > 0L, sumD / nObs, 0)
      centered <- sweep(dosage, 2L, mu, FUN = "-")
      centered[is.na(centered)] <- 0
      varVec <- colSums(centered * centered) / (nSamp - 1L)
      keepVarMask <- keepVarMask & varVec >= x@xvarCutoff
    }
    dosage <- dosage[, keepVarMask, drop = FALSE]
    mafVec <- mafVec[keepVarMask]
  } else {
    mafVec <- numeric(0)
  }

  # Mean-impute remaining missing dosage cells so downstream linear
  # algebra is well-defined; MAF was computed before imputation.
  if (anyNA(dosage)) {
    for (j in seq_len(ncol(dosage))) {
      col <- dosage[, j]
      na <- is.na(col)
      if (any(na)) {
        col[na] <- mean(col[!na])
        dosage[, j] <- col
      }
    }
  }

  list(
    geno       = dosage,
    variantIds = colnames(dosage),
    sampleIds  = rownames(dosage),
    maf        = mafVec
  )
}

#' @rdname getGenotypes
#' @export
setMethod("getGenotypes", "QtlDataset",
  function(x, traitId = NULL, region = NULL, cisWindow = NULL,
           samples = NULL, ...) {
    .qtlExtractBlock(x, traitId = traitId, region = region,
                     cisWindow = cisWindow, samples = samples)$geno
  })

#' @rdname getMaf
#' @export
setMethod("getMaf", "QtlDataset",
  function(x, region = NULL, cisWindow = NULL, samples = NULL, ...) {
    block <- .qtlExtractBlock(x, traitId = NULL, region = region,
                              cisWindow = cisWindow, samples = samples)
    out <- block$maf
    names(out) <- block$variantIds
    out
  })

#' @rdname getPhenotypes
#' @export
setMethod("getPhenotypes", "QtlDataset",
  function(x, contexts, traitId = NULL, region = NULL, ...) {
    if (missing(contexts) || is.null(contexts) || length(contexts) == 0L) {
      stop("`contexts` is required for getPhenotypes(QtlDataset). ",
           "Pass a character vector of one or more context names; ",
           "use getContexts(x) to list the available contexts.")
    }
    available <- names(x@phenotypes)
    bad <- setdiff(contexts, available)
    if (length(bad) > 0L) {
      stop("Unknown context(s): ", paste(bad, collapse = ", "),
           ". Available: ", paste(available, collapse = ", "))
    }
    out <- x@phenotypes[contexts]
    if (!is.null(traitId)) {
      out <- lapply(seq_along(out), function(i) {
        se <- out[[i]]
        ctx <- names(out)[[i]]
        present <- intersect(traitId, rownames(se))
        missing <- setdiff(traitId, rownames(se))
        if (length(missing) > 0L) {
          warning(sprintf("context '%s' is missing trait(s): %s",
                          ctx, paste(missing, collapse = ", ")))
        }
        se[present, , drop = FALSE]
      })
      names(out) <- contexts
    }
    if (!is.null(region)) {
      out <- lapply(out, function(se) {
        rr <- SummarizedExperiment::rowRanges(se)
        keep <- IRanges::overlapsAny(rr, region)
        se[keep, , drop = FALSE]
      })
      names(out) <- contexts
    }
    out
  })

#' @rdname getPhenotypeCovariates
#' @export
setMethod("getPhenotypeCovariates", "QtlDataset",
  function(x, contexts) {
    if (missing(contexts) || is.null(contexts) || length(contexts) == 0L) {
      stop("`contexts` is required.")
    }
    available <- names(x@phenotypes)
    bad <- setdiff(contexts, available)
    if (length(bad) > 0L) {
      stop("Unknown context(s): ", paste(bad, collapse = ", "))
    }
    out <- lapply(contexts, function(ctx) {
      se <- x@phenotypes[[ctx]]
      cd <- SummarizedExperiment::colData(se)
      as.matrix(as.data.frame(cd))
    })
    names(out) <- contexts
    out
  })

# Internal: residualize a numeric matrix Y (n x k) against a covariate
# matrix C (n x p) via pivoted QR decomposition. Adds an intercept column
# to C. When C is rank-deficient (e.g., union of all contexts' phenotype
# covariates includes collinear / duplicate columns), the pivoted QR drops
# the redundant columns automatically. Optionally rescales each residual
# column to unit standard deviation; constant-valued columns are left
# unchanged.
.qtlResidualizeQR <- function(Y, C, scaleResiduals = TRUE) {
  X <- if (is.null(C) || ncol(C) == 0L) {
    matrix(1, nrow = nrow(Y), ncol = 1L,
           dimnames = list(rownames(Y), "intercept"))
  } else {
    cbind(intercept = 1, C)
  }
  # `qr.resid` does not support LAPACK pivoted QR, so use `lm.fit`. It
  # handles rank-deficient designs gracefully via base-R's pivoted QR
  # internally — same effect the LAPACK path was meant to deliver.
  res <- stats::lm.fit(x = X, y = Y)$residuals
  res <- as.matrix(res)
  rownames(res) <- rownames(Y)
  colnames(res) <- colnames(Y)
  if (isTRUE(scaleResiduals)) {
    sds <- apply(res, 2L, function(v) stats::sd(v, na.rm = TRUE))
    # `sds == 0` exact-zero test is unreliable for residuals coming out of
    # lm.fit on a constant Y: roundoff gives sd ~ 1e-16 instead of 0, and
    # dividing the (also-tiny) residuals by it amplifies floating-point
    # noise to unit-scale. Treat anything below sqrt(.Machine$double.eps)
    # as effectively zero (column is constant) and skip rescaling.
    nearZero <- !is.finite(sds) | sds < sqrt(.Machine$double.eps)
    sds[nearZero] <- 1
    res[, nearZero] <- 0
    res <- sweep(res, 2L, sds, FUN = "/")
  }
  res
}

# Internal: validate and resolve the `*ToResidualize` argument against a
# set of contexts and the covariates actually present in those contexts'
# colData. Accepts either NULL (use all), a character vector (apply to all
# listed contexts), or a named list keyed by context. Returns a named list
# keyed by context giving the actual character vector of covariate names
# to use for that context (or character(0) if none). Errors when:
#   - a named-list key is not in `contexts`
#   - `contexts` contains entries missing from a supplied named-list
#     (per the rule: named-list keys must equal `contexts`)
#   - an explicitly requested name matches no actual covariate
.qtlResolvePhenoSelection <- function(x, contexts, toResidualize) {
  resolveOne <- function(ctx, requested) {
    se <- x@phenotypes[[ctx]]
    avail <- colnames(SummarizedExperiment::colData(se))
    if (is.null(requested)) return(avail)
    keep <- intersect(requested, avail)
    if (length(keep) != length(requested)) {
      missingNames <- setdiff(requested, avail)
      stop(sprintf(
        "phenotypeCovariatesToResidualize: context '%s' has no covariate(s) named: %s",
        ctx, paste(missingNames, collapse = ", ")))
    }
    keep
  }
  if (is.null(toResidualize)) {
    out <- lapply(contexts, resolveOne, requested = NULL)
    names(out) <- contexts
    return(out)
  }
  if (is.list(toResidualize)) {
    if (is.null(names(toResidualize)) ||
        any(!nzchar(names(toResidualize)))) {
      stop("phenotypeCovariatesToResidualize: when supplied as a list, ",
           "it must be named with context names.")
    }
    badKeys <- setdiff(names(toResidualize), contexts)
    if (length(badKeys) > 0L) {
      stop("phenotypeCovariatesToResidualize: list key(s) not in `contexts`: ",
           paste(badKeys, collapse = ", "))
    }
    missingKeys <- setdiff(contexts, names(toResidualize))
    if (length(missingKeys) > 0L) {
      stop("phenotypeCovariatesToResidualize: list does not cover all ",
           "`contexts`. Per-context lists must have exactly the same ",
           "context set as `contexts`. Missing keys: ",
           paste(missingKeys, collapse = ", "))
    }
    out <- lapply(contexts, function(ctx) resolveOne(ctx, toResidualize[[ctx]]))
    names(out) <- contexts
    return(out)
  }
  if (is.character(toResidualize)) {
    out <- lapply(contexts, resolveOne, requested = toResidualize)
    names(out) <- contexts
    return(out)
  }
  stop("phenotypeCovariatesToResidualize must be NULL, a character vector, ",
       "or a named list keyed by context.")
}

# Internal: validate the genotype-covariate selection vector. Returns
# character(0) when nothing selected, the resolved set otherwise.
.qtlResolveGenoSelection <- function(x, toResidualize) {
  avail <- colnames(x@genotypeCovariates)
  if (is.null(avail)) avail <- character(0)
  if (is.null(toResidualize)) return(avail)
  keep <- intersect(toResidualize, avail)
  if (length(keep) != length(toResidualize)) {
    missingNames <- setdiff(toResidualize, avail)
    stop("genotypeCovariatesToResidualize: no covariate(s) named: ",
         paste(missingNames, collapse = ", "))
  }
  keep
}

# Internal: build the covariate matrix used for residualization, given a
# set of contexts, the resolved per-context phenotype selections, and the
# resolved genotype-covariate selection. Honors the inclusion flags. For
# `length(contexts) == 1` (per-context mode), the per-context phenotype
# covariates and the genotype covariates are taken with no cross-context
# alignment. For `length(contexts) >= 2` (joint mode), per-context
# phenotype covariates from all listed contexts are concatenated
# (prefixed with "{context}." to keep same-named columns distinct) and
# the sample set is intersected across all contributing matrices.
# Returns a single matrix (with rownames = sample IDs) or NULL.
.qtlBuildResidualizationDesign <- function(x, contexts,
                                           phenoSelection,
                                           genoSelection,
                                           includePheno, includeGeno) {
  perContext <- list()
  if (includePheno) {
    for (ctx in contexts) {
      keep <- phenoSelection[[ctx]]
      if (length(keep) == 0L) next
      se <- x@phenotypes[[ctx]]
      cd <- as.matrix(as.data.frame(SummarizedExperiment::colData(se)))
      cdMat <- cd[, keep, drop = FALSE]
      colnames(cdMat) <- paste0(ctx, ".", colnames(cdMat))
      if (is.null(rownames(cdMat))) {
        rownames(cdMat) <- as.character(
          rownames(SummarizedExperiment::colData(se)))
      }
      perContext[[ctx]] <- cdMat
    }
  }
  gCov <- if (includeGeno && length(genoSelection) > 0L) {
    x@genotypeCovariates[, genoSelection, drop = FALSE]
  } else {
    matrix(numeric(0), nrow = 0, ncol = 0)
  }
  haveAny <- length(perContext) > 0L ||
    (!is.null(gCov) && ncol(gCov) > 0L)
  if (!haveAny) return(NULL)

  sampleSets <- list()
  for (mat in perContext) {
    if (!is.null(rownames(mat))) {
      sampleSets[[length(sampleSets) + 1L]] <- rownames(mat)
    }
  }
  if (!is.null(gCov) && ncol(gCov) > 0L && !is.null(rownames(gCov))) {
    sampleSets[[length(sampleSets) + 1L]] <- rownames(gCov)
  }
  common <- if (length(sampleSets) == 0L) character(0)
            else Reduce(intersect, sampleSets)
  if (length(common) == 0L) return(NULL)

  blocks <- list()
  for (mat in perContext) {
    blocks[[length(blocks) + 1L]] <- mat[common, , drop = FALSE]
  }
  if (!is.null(gCov) && ncol(gCov) > 0L) {
    blocks[[length(blocks) + 1L]] <- gCov[common, , drop = FALSE]
  }
  do.call(cbind, blocks)
}

# Internal: resolve a (convenience, precise) flag pair to a single boolean.
# `missing*` arguments are passed as the result of `missing()` evaluated in
# the calling method to detect whether the user explicitly set the value.
# Rules:
#   - both missing: returns TRUE (the documented default)
#   - only convenience set: returns convenience
#   - only precise set: returns precise
#   - both set: must agree, else error
.qtlResolveResidualizationFlag <- function(conveniencePassed, convenienceMissing,
                                           precisePassed, preciseMissing,
                                           convenienceName, preciseName) {
  if (preciseMissing && convenienceMissing) return(TRUE)
  if (preciseMissing) return(isTRUE(conveniencePassed))
  if (convenienceMissing) return(isTRUE(precisePassed))
  if (isTRUE(conveniencePassed) != isTRUE(precisePassed)) {
    stop(sprintf(
      "Conflicting values: `%s` = %s and `%s` = %s. Set only one, or ",
      convenienceName, conveniencePassed,
      preciseName, precisePassed),
      "pass consistent values.")
  }
  isTRUE(precisePassed)
}

#' @rdname getResidualizedGenotypes
#' @export
setMethod("getResidualizedGenotypes", "QtlDataset",
  function(x, contexts, traitId = NULL, region = NULL, cisWindow = NULL,
           samples = NULL,
           phenotypeCovariatesToResidualize = NULL,
           genotypeCovariatesToResidualize = NULL,
           residualizePhenotypeCovariates = TRUE,
           residualizeGenotypeCovariates  = TRUE,
           residualizePhenotypeCovariatesFromGenotypes = NULL,
           residualizeGenotypeCovariatesFromGenotypes  = NULL,
           ...) {
    if (missing(contexts) || is.null(contexts) || length(contexts) == 0L) {
      stop("`contexts` is required for getResidualizedGenotypes(QtlDataset). ",
           "Use getContexts(x) to list the available contexts. ",
           "Pass a single context for per-context mode or multiple ",
           "contexts for joint mode (sample intersection).")
    }
    bad <- setdiff(contexts, names(x@phenotypes))
    if (length(bad) > 0L) {
      stop("Unknown context(s): ", paste(bad, collapse = ", "))
    }

    # Resolve inclusion flags (convenience vs precise).
    convPhenoMissing <- missing(residualizePhenotypeCovariates)
    convGenoMissing  <- missing(residualizeGenotypeCovariates)
    precPhenoMissing <- missing(residualizePhenotypeCovariatesFromGenotypes) ||
                       is.null(residualizePhenotypeCovariatesFromGenotypes)
    precGenoMissing  <- missing(residualizeGenotypeCovariatesFromGenotypes) ||
                       is.null(residualizeGenotypeCovariatesFromGenotypes)
    includePheno <- .qtlResolveResidualizationFlag(
      residualizePhenotypeCovariates, convPhenoMissing,
      residualizePhenotypeCovariatesFromGenotypes, precPhenoMissing,
      "residualizePhenotypeCovariates",
      "residualizePhenotypeCovariatesFromGenotypes")
    includeGeno <- .qtlResolveResidualizationFlag(
      residualizeGenotypeCovariates, convGenoMissing,
      residualizeGenotypeCovariatesFromGenotypes, precGenoMissing,
      "residualizeGenotypeCovariates",
      "residualizeGenotypeCovariatesFromGenotypes")

    # Resolve the covariate selections.
    phenoSel <- .qtlResolvePhenoSelection(x, contexts,
                                          phenotypeCovariatesToResidualize)
    genoSel  <- .qtlResolveGenoSelection(x, genotypeCovariatesToResidualize)

    G <- getGenotypes(x, traitId = traitId, region = region,
                      cisWindow = cisWindow, samples = samples)
    if (ncol(G) == 0L) return(G)

    C <- .qtlBuildResidualizationDesign(
      x, contexts = contexts,
      phenoSelection = phenoSel,
      genoSelection  = genoSel,
      includePheno   = includePheno,
      includeGeno    = includeGeno)
    if (!is.null(C)) {
      common <- intersect(rownames(G), rownames(C))
      if (length(common) == 0L) {
        stop("No samples in common between the genotype matrix and the ",
             "covariate matrix for contexts: ",
             paste(contexts, collapse = ", "))
      }
      G <- G[common, , drop = FALSE]
      C <- C[common, , drop = FALSE]
    }
    .qtlResidualizeQR(G, C, scaleResiduals = x@scaleResiduals)
  })

#' @rdname getResidualizedPhenotypes
#' @export
setMethod("getResidualizedPhenotypes", "QtlDataset",
  function(x, contexts, traitId = NULL, region = NULL,
           phenotypeCovariatesToResidualize = NULL,
           genotypeCovariatesToResidualize = NULL,
           residualizePhenotypeCovariates = TRUE,
           residualizeGenotypeCovariates  = TRUE,
           residualizePhenotypeCovariatesFromPhenotypes = NULL,
           residualizeGenotypeCovariatesFromPhenotypes  = NULL,
           ...) {
    if (missing(contexts) || is.null(contexts) || length(contexts) == 0L) {
      stop("`contexts` is required for getResidualizedPhenotypes().")
    }
    bad <- setdiff(contexts, names(x@phenotypes))
    if (length(bad) > 0L) {
      stop("Unknown context(s): ", paste(bad, collapse = ", "))
    }

    convPhenoMissing <- missing(residualizePhenotypeCovariates)
    convGenoMissing  <- missing(residualizeGenotypeCovariates)
    precPhenoMissing <- missing(residualizePhenotypeCovariatesFromPhenotypes) ||
                       is.null(residualizePhenotypeCovariatesFromPhenotypes)
    precGenoMissing  <- missing(residualizeGenotypeCovariatesFromPhenotypes) ||
                       is.null(residualizeGenotypeCovariatesFromPhenotypes)
    includePheno <- .qtlResolveResidualizationFlag(
      residualizePhenotypeCovariates, convPhenoMissing,
      residualizePhenotypeCovariatesFromPhenotypes, precPhenoMissing,
      "residualizePhenotypeCovariates",
      "residualizePhenotypeCovariatesFromPhenotypes")
    includeGeno <- .qtlResolveResidualizationFlag(
      residualizeGenotypeCovariates, convGenoMissing,
      residualizeGenotypeCovariatesFromPhenotypes, precGenoMissing,
      "residualizeGenotypeCovariates",
      "residualizeGenotypeCovariatesFromPhenotypes")

    phenoSel <- .qtlResolvePhenoSelection(x, contexts,
                                          phenotypeCovariatesToResidualize)
    genoSel  <- .qtlResolveGenoSelection(x, genotypeCovariatesToResidualize)

    Yraw <- getPhenotypes(x, contexts = contexts, traitId = traitId,
                          region = region)
    C <- .qtlBuildResidualizationDesign(
      x, contexts = contexts,
      phenoSelection = phenoSel,
      genoSelection  = genoSel,
      includePheno   = includePheno,
      includeGeno    = includeGeno)

    out <- lapply(contexts, function(ctx) {
      se <- Yraw[[ctx]]
      Y <- t(SummarizedExperiment::assay(se))  # samples x traits
      if (!is.null(C)) {
        common <- intersect(rownames(Y), rownames(C))
        if (length(common) == 0L) {
          stop(sprintf(
            "context '%s': no samples shared between phenotype data and ",
            "the resolved covariate matrix.", ctx))
        }
        Y <- Y[common, , drop = FALSE]
        Cctx <- C[common, , drop = FALSE]
      } else {
        Cctx <- NULL
      }
      .qtlResidualizeQR(Y, Cctx, scaleResiduals = x@scaleResiduals)
    })
    names(out) <- contexts
    out
  })

# =============================================================================
# Per-entry payload constructors and accessors
# =============================================================================

#' @title Create a FineMappingEntry Object
#' @description Construct a \code{FineMappingEntry} payload for one
#'   \code{(study, context, trait, method)} row of a
#'   \code{FineMappingResult} collection.
#' @param variantIds Character vector of variant IDs.
#' @param trimmedFit Method-specific fit object.
#' @param topLoci Long-format \code{data.frame}.
#' @param sumstats Optional list of summary statistics, or \code{NULL}.
#' @return A \code{FineMappingEntry} object.
#' @export
FineMappingEntry <- function(variantIds, trimmedFit, topLoci,
                             sumstats = NULL) {
  obj <- new("FineMappingEntry",
             variantIds = as.character(variantIds),
             trimmedFit = trimmedFit,
             topLoci    = as.data.frame(topLoci),
             sumstats   = sumstats)
  validObject(obj)
  obj
}

#' @title Create a TwasWeightsEntry Object
#' @description Construct a \code{TwasWeightsEntry} payload for one
#'   \code{(study, context, trait, method)} row of a \code{TwasWeights}
#'   collection.
#' @param variantIds Character vector of variant IDs.
#' @param weights Numeric vector or matrix.
#' @param fits Optional method-specific fit object.
#' @param cvPerformance Optional list of CV metrics.
#' @param standardized Logical (length 1).
#' @param dataType Optional data-type tag.
#' @return A \code{TwasWeightsEntry} object.
#' @export
TwasWeightsEntry <- function(variantIds, weights, fits = NULL,
                             cvPerformance = NULL, standardized = FALSE,
                             dataType = NULL) {
  obj <- new("TwasWeightsEntry",
             variantIds    = as.character(variantIds),
             weights       = weights,
             fits          = fits,
             cvPerformance = cvPerformance,
             standardized  = isTRUE(standardized),
             dataType      = dataType)
  validObject(obj)
  obj
}

# Per-entry accessors (reuse the existing generics; these methods read
# slots from the payload classes directly).

#' @rdname getVariantIds
#' @export
setMethod("getVariantIds", "FineMappingEntry",
          function(x, ...) x@variantIds)

#' @rdname getTrimmedFit
#' @export
setMethod("getTrimmedFit", "FineMappingEntry",
          function(x, ...) x@trimmedFit)

#' @rdname getTopLoci
#' @export
setMethod("getTopLoci", "FineMappingEntry",
          function(x, ...) x@topLoci)

#' @rdname getPip
#' @export
setMethod("getPip", "FineMappingEntry", function(x, ...) {
  tl <- x@topLoci
  if (nrow(tl) == 0L || !"pip" %in% names(tl)) return(numeric(0))
  setNames(tl$pip, tl$variant_id)
})

#' @rdname getCs
#' @export
setMethod("getCs", "FineMappingEntry",
  function(x, coverage = 0.95, ...) {
    tl <- x@topLoci
    if (nrow(tl) == 0L) return(data.frame())
    csCol <- grep(paste0("^cs.*", coverage * 100), names(tl), value = TRUE)
    if (length(csCol) == 0L && "cs" %in% names(tl)) csCol <- "cs"
    if (length(csCol) == 0L) return(data.frame())
    tl[tl[[csCol[1L]]] > 0, , drop = FALSE]
  })

#' @rdname getWeights
#' @export
setMethod("getWeights", "TwasWeightsEntry",
          function(x, ...) x@weights)

#' @rdname getVariantIds
#' @export
setMethod("getVariantIds", "TwasWeightsEntry",
          function(x, ...) x@variantIds)

#' @rdname getFits
#' @export
setMethod("getFits", "TwasWeightsEntry",
          function(x, ...) x@fits)

#' @rdname getCvPerformance
#' @export
setMethod("getCvPerformance", "TwasWeightsEntry",
          function(x, ...) x@cvPerformance)

#' @rdname getStandardized
#' @export
setMethod("getStandardized", "TwasWeightsEntry",
          function(x, ...) x@standardized)

#' @rdname getDataType
#' @export
setMethod("getDataType", "TwasWeightsEntry",
          function(x, ...) x@dataType)

# =============================================================================
# MultiTaskQtlDataset constructor and accessors
# =============================================================================

#' @title Create a MultiTaskQtlDataset Object
#' @description Construct a \code{MultiTaskQtlDataset} S4 object from a
#'   named list of \code{QtlDataset} objects (individual-level studies)
#'   and an optional \code{QtlSumStats} of summary-statistic-only
#'   studies. The total study count must be at least two, satisfied by
#'   either (a) at least two \code{qtlDatasets} entries, or (b) at least
#'   one \code{qtlDatasets} entry plus a non-empty \code{sumStats}.
#' @param qtlDatasets A named list of \code{QtlDataset} objects, keyed
#'   by study identifier.
#' @param sumStats An optional \code{QtlSumStats} collection. Default
#'   \code{NULL}.
#' @return A \code{MultiTaskQtlDataset} object.
#' @export
MultiTaskQtlDataset <- function(qtlDatasets, sumStats = NULL) {
  obj <- new("MultiTaskQtlDataset",
             qtlDatasets = qtlDatasets,
             sumStats    = sumStats)
  validObject(obj)
  obj
}

#' @rdname getQtlDatasets
#' @export
setMethod("getQtlDatasets", "MultiTaskQtlDataset",
          function(x) x@qtlDatasets)

#' @rdname getSumStats
#' @export
setMethod("getSumStats", "MultiTaskQtlDataset",
  function(x, ...) {
    if (length(list(...)) > 0L) {
      stop("getSumStats(MultiTaskQtlDataset) does not accept selection ",
           "arguments; it returns the embedded QtlSumStats collection ",
           "(use getSumStats() on that result to fetch one entry).")
    }
    x@sumStats
  })

#' @rdname getStudy
#' @export
setMethod("getStudy", "MultiTaskQtlDataset", function(x) {
  fromQtl <- names(x@qtlDatasets)
  fromSs  <- if (is.null(x@sumStats)) character(0)
             else unique(as.character(x@sumStats$study))
  unique(c(fromQtl, fromSs))
})

# =============================================================================
# topLoci GRanges conversion
# =============================================================================

#' @title Convert topLoci to GRanges
#' @description Convert a long-format topLoci data.frame (with variant_id
#'   in chr:pos:A2:A1 format) to a GRanges object with all metadata columns.
#' @param topLoci A data.frame with a variant_id column encoding
#'   chr:pos:A2:A1.
#' @return A GRanges object with all original columns as metadata.
#' @export
topLociToGranges <- function(topLoci) {
  if (is.null(topLoci) || nrow(topLoci) == 0) {
    return(GRanges())
  }
  parsed <- parseVariantId(topLoci$variant_id)
  chr <- paste0("chr", parsed$chrom)
  gr <- GRanges(
    seqnames = chr,
    ranges = IRanges(start = parsed$pos, width = 1L)
  )
  mcols(gr) <- DataFrame(topLoci)
  gr
}

# =============================================================================
# GenotypeHandle accessors
# =============================================================================

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

# =============================================================================
# LdStatistic / LdEigen / LdScore accessors
# =============================================================================

#' @rdname getSnpInfo
#' @export
setMethod("getSnpInfo", "LdStatistic", function(x) x@snpInfo)

#' @rdname getNRef
#' @export
setMethod("getNRef", "LdStatistic", function(x) x@nRef)

#' @rdname getInSample
#' @export
setMethod("getInSample", "LdStatistic", function(x) x@inSample)

#' @rdname getLdBlocks
#' @export
setMethod("getLdBlocks", "LdStatistic", function(x) x@ldBlocks)

#' @rdname getEigenList
#' @export
setMethod("getEigenList", "LdEigen", function(x) x@eigenList)

#' @rdname getLdScores
#' @export
setMethod("getLdScores", "LdScore", function(x) x@ldScores)

#' @rdname getLdScoreWeights
#' @export
setMethod("getLdScoreWeights", "LdScore", function(x) x@ldScoreWeights)

#' @rdname getLdMatrixList
#' @export
setMethod("getLdMatrixList", "LdScore", function(x) x@ldMatrixList)

# =============================================================================
# AnnotationMatrix accessors
# =============================================================================

#' @rdname getAnnotations
#' @export
setMethod("getAnnotations", "AnnotationMatrix", function(x) x@annotations)

#' @rdname getAnnotationMeta
#' @export
setMethod("getAnnotationMeta", "AnnotationMatrix",
          function(x) x@annotationMeta)

#' @rdname getSnpRanges
#' @export
setMethod("getSnpRanges", "AnnotationMatrix", function(x) x@snpRanges)

# =============================================================================
# LdBlocks accessor
# =============================================================================

#' @rdname getBlocks
#' @export
setMethod("getBlocks", "LdBlocks", function(x) x@blocks)

# =============================================================================
# LdData accessors (genotypeHandle, snpIdx, nRef)
# =============================================================================

#' @rdname getGenotypeHandle
#' @export
setMethod("getGenotypeHandle", "LdData", function(x) x@genotypeHandle)

#' @rdname getMixtureWeights
#' @export
setMethod("getMixtureWeights", "LdData", function(x) x@mixtureWeights)

#' @rdname getSnpIdx
#' @export
setMethod("getSnpIdx", "LdData", function(x) x@snpIdx)

#' @rdname getNRef
#' @export
setMethod("getNRef", "LdData", function(x) x@nRef)

# =============================================================================
# H2Estimate accessors (tauBlocks, h2)
# =============================================================================

#' @rdname getTauBlocks
#' @export
setMethod("getTauBlocks", "H2Estimate", function(x) x@tauBlocks)

#' @rdname getH2
#' @export
setMethod("getH2", "H2Estimate", function(x) x@h2)

# =============================================================================
# Cross-class getGenome methods (LdStatistic / AnnotationMatrix / LdBlocks)
# =============================================================================

#' @rdname getGenome
#' @export
setMethod("getGenome", "LdStatistic", function(x, ...) x@genome)

#' @rdname getGenome
#' @export
setMethod("getGenome", "AnnotationMatrix", function(x, ...) x@genome)

#' @rdname getGenome
#' @export
setMethod("getGenome", "LdBlocks", function(x, ...) x@genome)
