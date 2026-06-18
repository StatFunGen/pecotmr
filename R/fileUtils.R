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

#' @importFrom vroom vroom
#' @importFrom dplyr as_tibble mutate filter
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom stringr str_detect

# Internal helper: read a region from a tabix-indexed file via Rsamtools
readTabixRegion <- function(file, region, useColNames) {
  tbx <- TabixFile(file)
  parsed <- parseRegion(region)
  # Match chromosome naming convention in the tabix index
  chrom <- as.character(parsed$chrom)
  tbxSeqnames <- seqnamesTabix(tbx)
  if (any(grepl("^chr", tbxSeqnames))) {
    chrom <- paste0("chr", chrom)
  }
  gr <- GRanges(
    seqnames = chrom,
    ranges = IRanges(start = parsed$start, end = parsed$end)
  )
  lines <- scanTabix(tbx, param = gr)[[1]]
  if (length(lines) == 0) return(NULL)

  # Get header for column names
  colNamesVec <- NULL
  if (useColNames) {
    hdr <- headerTabix(tbx)$header
    if (length(hdr) > 0) {
      lastHdr <- hdr[length(hdr)]
      colNamesVec <- strsplit(sub("^#", "", lastHdr), "\t")[[1]]
    }
  }

  # Parse tab-delimited lines
  txt <- paste(lines, collapse = "\n")
  if (!is.null(colNamesVec)) {
    as.data.frame(vroom(I(txt), delim = "\t", col_names = colNamesVec,
                               show_col_types = FALSE))
  } else {
    as.data.frame(vroom(I(txt), delim = "\t", col_names = useColNames,
                               show_col_types = FALSE))
  }
}

tabixRegion <- function(file, region, tabixHeader = "auto", target = "", targetColumnIndex = "") {
  if (!file.exists(file)) {
    stop("Input file does not exist: ", file)
  }

  useColNames <- if (identical(tabixHeader, FALSE)) FALSE else TRUE

  cmdOutput <- tryCatch(
    readTabixRegion(file, region, useColNames),
    error = function(e) NULL
  )

  if (!is.null(cmdOutput) && target != "" && targetColumnIndex != "") {
    cmdOutput <- cmdOutput %>%
      filter(str_detect(.[[targetColumnIndex]], target))
  } else if (!is.null(cmdOutput) && target != "") {
    cmdOutput <- cmdOutput %>%
      mutate(text = apply(., 1, function(row) paste(row, collapse = "_"))) %>%
      filter(str_detect(text, target)) %>%
      select(-text)
  }

  if (is.null(cmdOutput) || nrow(cmdOutput) == 0) {
    return(tibble())
  }

  cmdOutput %>%
    as_tibble() %>%
    mutate(
      !!names(.)[1] := as.character(.[[1]]),
      !!names(.)[2] := as.numeric(.[[2]])
    )
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

#' @importFrom purrr map
#' @importFrom readr read_delim cols
#' @importFrom dplyr select mutate across everything
#' @importFrom magrittr %>%
#' @noRd
readSingleCovariate <- function(path) {
  rawDf <- read_delim(path, "\t", col_types = cols(.default = "c")) %>% select(-1)
  df <- rawDf
  nonNumeric <- character()
  for (nm in names(df)) {
    values <- trimws(as.character(df[[nm]]))
    converted <- suppressWarnings(as.numeric(values))
    bad <- !is.na(values) & values != "" & is.na(converted)
    if (any(bad)) {
      nonNumeric <- c(nonNumeric, nm)
    } else {
      df[[nm]] <- converted
    }
  }
  if (length(nonNumeric) > 0) {
    stop("Non-numeric columns found in covariate file ", path, ": ",
         paste(nonNumeric, collapse = ", "),
         ". All columns except the first (sample ID) must be numeric.")
  }
  df %>% mutate(across(everything(), as.numeric)) %>% t()
}

#' @noRd
loadCovariateData <- function(covariatePath) {
  # Validate all covariate files exist
  missing <- covariatePath[!file.exists(covariatePath)]
  if (length(missing) > 0) {
    stop("Covariate file(s) not found: ", paste(missing, collapse = ", "))
  }
  return(map(covariatePath, readSingleCovariate))
}

NoPhenotypeError <- function(message) {
  structure(list(message = message), class = c("NoPhenotypeError", "error", "condition"))
}

#' @importFrom purrr map2 compact
#' @importFrom readr read_delim cols
#' @importFrom dplyr filter select mutate across everything
#' @importFrom magrittr %>%
#' @noRd
loadPhenotypeData <- function(phenotypePath, region, extractRegionName = NULL, regionNameCol = NULL, tabixHeader = TRUE) {
  if (is.null(extractRegionName)) {
    extractRegionName <- rep(list(NULL), length(phenotypePath))
  } else if (is.list(extractRegionName) && length(extractRegionName) != length(phenotypePath)) {
    stop("extract_region_name must be NULL or a list with the same length as phenotype_path.")
  } else if (!is.null(extractRegionName) && !is.list(extractRegionName)) {
    stop("extract_region_name must be NULL or a list.")
  }

  # Use `map2` to iterate over `phenotype_path` and `extract_region_name` simultaneously
  phenotypeDataRaw <- map2(phenotypePath, extractRegionName, ~ {
    tabixData <- if (!is.null(region)) tabixRegion(.x, region, tabixHeader = tabixHeader) else read_delim(.x, "\t", col_types = cols())
    if (nrow(tabixData) == 0) {
      message(paste("Phenotype file ", .x, " is empty for the specified region", if (is.null(region)) "" else region))
      return(NULL)
    }
    if (!is.null(.y) && is.vector(.y) && !is.null(regionNameCol) && (regionNameCol %% 1 == 0)) {
      if (regionNameCol <= ncol(tabixData)) {
        regionColName <- colnames(tabixData)[regionNameCol]
        tabixData <- tabixData %>%
          filter(.data[[regionColName]] %in% .y) %>%
          t()
        colnames(tabixData) <- tabixData[regionNameCol, ]
        return(tabixData)
      } else {
        stop("region_name_col is out of bounds for the number of columns in tabix_data.")
      }
    } else {
      result <- tabixData %>% t()
      # Assign region names from region_name_col if available
      if (!is.null(regionNameCol) && (regionNameCol %% 1 == 0) && regionNameCol <= ncol(tabixData)) {
        colnames(result) <- tabixData[[regionNameCol]]
      }
      return(result)
    }
  })

  # Track which indices had non-NULL data, then remove NULLs
  keptIndices <- which(vapply(phenotypeDataRaw, Negate(is.null), logical(1)))
  phenotypeData <- phenotypeDataRaw[keptIndices]

  # Check if all phenotype files are empty
  if (length(phenotypeData) == 0) {
    stop(NoPhenotypeError(paste("All phenotype files are empty for the specified region", if (!is.null(region)) "" else region)))
  }
  # Store kept indices as attribute so callers can align covariates/conditions
  attr(phenotypeData, "kept_indices") <- keptIndices
  return(phenotypeData)
}

#' @importFrom purrr map
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @noRd
extractPhenotypeCoordinates <- function(phenotypeList) {
  return(map(phenotypeList, ~ t(.x[1:3, ]) %>%
    as_tibble() %>%
    mutate(start = as.numeric(start), end = as.numeric(end))))
}

#' @importFrom magrittr %>%
#' @noRd
filterByCommonSamples <- function(dat, commonSamples) {
  dat[commonSamples, , drop = FALSE] %>% .[order(rownames(.)), , drop = FALSE]
}

#' @importFrom tibble tibble
#' @importFrom dplyr mutate select
#' @importFrom purrr map map2
#' @importFrom magrittr %>%
#' @noRd
prepareDataList <- function(genoBed, phenotype, covariate, imissCutoff, mafCutoff, macCutoff, xvarCutoff, phenotypeHeader = 4, keepSamples = NULL) {
  dataList <- tibble(
    covar = covariate,
    Y = lapply(phenotype, function(x) apply(x[-c(1:phenotypeHeader), , drop = FALSE], c(1, 2), as.numeric))
  ) %>%
    mutate(
      # Determine common complete samples across Y, covar, and geno_bed, considering missing values
      common_complete_samples = map2(covar, Y, ~ {
        covar_non_na <- rownames(.x)[!apply(.x, 1, function(row) all(is.na(row)))]
        y_non_na <- rownames(.y)[!apply(.y, 1, function(row) all(is.na(row)))]
        if (length(intersect(intersect(covar_non_na, y_non_na), rownames(genoBed))) == 0) {
          stop("No common complete samples between genotype and phenotype/covariate data")
        }
        intersect(intersect(covar_non_na, y_non_na), rownames(genoBed))
      }),
      # Further intersect with keep_samples if provided
      common_complete_samples = if (!is.null(keepSamples) && length(keepSamples) > 0) {
        map(common_complete_samples, ~ intersect(.x, keepSamples))
      } else {
        common_complete_samples
      },
      # Determine dropped samples before filtering
      dropped_samples_covar = map2(covar, common_complete_samples, ~ setdiff(rownames(.x), .y)),
      dropped_samples_Y = map2(Y, common_complete_samples, ~ setdiff(rownames(.x), .y)),
      dropped_samples_X = map(common_complete_samples, ~ setdiff(rownames(genoBed), .x)),
      # Filter data based on common complete samples
      Y = map2(Y, common_complete_samples, ~ filterByCommonSamples(.x, .y)),
      covar = map2(covar, common_complete_samples, ~ filterByCommonSamples(.x, .y)),
      # Apply filter_X on the geno_bed data filtered by common complete samples and then format column names
      X = map(common_complete_samples, ~ {
        filteredGenoBed <- filterByCommonSamples(genoBed, .x)
        macVal <- if (nrow(filteredGenoBed) == 0) 0 else (macCutoff / (2 * nrow(filteredGenoBed)))
        mafVal <- max(mafCutoff, macVal)
        filteredData <- filterX(filteredGenoBed, imissCutoff, mafVal, varThresh = xvarCutoff)
        colnames(filteredData) <- normalizeVariantId(colnames(filteredData)) # Normalize to canonical format
        filteredData
      })
    ) %>%
    select(covar, Y, X, dropped_samples_Y, dropped_samples_X, dropped_samples_covar)
  return(dataList)
}

#' @importFrom purrr map
#' @importFrom dplyr intersect
#' @importFrom stringr str_split_fixed
#' @importFrom magrittr %>%
#' @noRd
prepareXMatrix <- function(genoBed, dataList, imissCutoff, mafCutoff, macCutoff, xvarCutoff) {
  # Calculate the union of all samples from data_list: any of X, covar and Y would do
  allSamplesUnion <- map(dataList$covar, ~ rownames(.x)) %>%
    unlist() %>%
    unique()
  # Find the intersection of these samples with the samples in geno_bed
  commonSamples <- intersect(allSamplesUnion, rownames(genoBed))
  # Filter geno_bed using common_samples
  XFiltered <- filterByCommonSamples(genoBed, commonSamples)
  # Calculate MAF cutoff considering the number of common samples
  mafVal <- max(mafCutoff, macCutoff / (2 * length(commonSamples)))
  # Apply further filtering on X
  XFiltered <- filterX(XFiltered, imissCutoff, mafVal, xvarCutoff)
  colnames(XFiltered) <- normalizeVariantId(colnames(XFiltered))

  # To keep a log message
  variants <- str_split_fixed(colnames(XFiltered), ":", 3)
  message(paste0("Dimension of input genotype data is ", nrow(XFiltered), " rows and ", ncol(XFiltered), " columns for genomic region of ", variants[1, 1], ":", min(as.integer(variants[, 2])), "-", max(as.integer(variants[, 2]))))
  return(XFiltered)
}

#' @importFrom purrr map map2
#' @importFrom dplyr mutate
#' @importFrom stats lm.fit sd
#' @importFrom magrittr %>%
#' @noRd
addXResiduals <- function(dataList, scaleResiduals = FALSE) {
  # Compute residuals for X and add them to data_list
  dataList <- dataList %>%
    mutate(
      lm_res_X = map2(X, covar, ~ .lm.fit(x = cbind(1, .y), y = .x)$residuals %>% as.matrix()),
      X_resid_mean = map(lm_res_X, ~ apply(.x, 2, mean)),
      X_resid_sd = map(lm_res_X, ~ apply(.x, 2, sd)),
      X_resid = map(lm_res_X, ~ {
        if (scaleResiduals) {
          scale(.x)
        } else {
          .x
        }
      })
    )

  return(dataList)
}

#' @importFrom purrr map map2
#' @importFrom dplyr mutate
#' @importFrom stats lm.fit sd
#' @importFrom magrittr %>%
#' @noRd
addYResiduals <- function(dataList, conditions, scaleResiduals = FALSE) {
  # Compute residuals, their mean, and standard deviation, and add them to data_list
  dataList <- dataList %>%
    mutate(
      lm_res = map2(Y, covar, ~ {
        res <- .lm.fit(x = cbind(1, .y), y = .x)$residuals %>% as.matrix()
        colnames(res) <- colnames(.x)
        res
      }),
      Y_resid_mean = map(lm_res, ~ apply(.x, 2, mean)),
      Y_resid_sd = map(lm_res, ~ apply(.x, 2, sd)),
      Y_resid = map(lm_res, ~ {
        if (scaleResiduals) {
          scale(.x)
        } else {
          .x
        }
      })
    )

  names(dataList$Y_resid) <- conditions

  return(dataList)
}

#' (Deprecated) Load regional association data
#'
#' \strong{Deprecated.} The on-disk regional loaders have been removed.
#' Construct a \code{\link{QtlDataset}} directly from in-memory inputs
#' (a \code{\link{GenotypeHandle}} plus a named list of per-context
#' \code{SummarizedExperiment} phenotypes). Variant- and sample-level QC
#' (\code{mafCutoff}, \code{macCutoff}, \code{xvarCutoff},
#' \code{imissCutoff}, \code{keepSamples}, \code{keepVariants}) are
#' constructor arguments on \code{QtlDataset()} and are applied lazily at
#' extraction time inside \code{getGenotypes()} and
#' \code{getResidualizedGenotypes()}.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
loadRegionalAssociationData <- function(...) {
  .Deprecated(new = "QtlDataset", package = "pecotmr",
    msg = paste(
      "loadRegionalAssociationData() has been removed.",
      "Build a QtlDataset() directly from a GenotypeHandle and a named",
      "list of per-context SummarizedExperiment phenotypes; pass",
      "mafCutoff / macCutoff / xvarCutoff / imissCutoff /",
      "keepSamples / keepVariants to the constructor."))
  invisible(NULL)
}

#' (Deprecated) Load regional univariate data
#'
#' \strong{Deprecated.} Use \code{\link{QtlDataset}} with a single context.
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
loadRegionalUnivariateData <- function(...) {
  .Deprecated(new = "QtlDataset", package = "pecotmr",
    msg = paste(
      "loadRegionalUnivariateData() has been removed. Build a QtlDataset()",
      "with a single context entry in the phenotypes list."))
  invisible(NULL)
}

#' (Deprecated) Load regional data for regression modeling
#'
#' \strong{Deprecated.} Use \code{\link{QtlDataset}}.
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
loadRegionalRegressionData <- function(...) {
  .Deprecated(new = "QtlDataset", package = "pecotmr",
    msg = paste(
      "loadRegionalRegressionData() has been removed. Build a QtlDataset()",
      "directly; per-condition residualized genotype/phenotype views are",
      "available via getResidualizedGenotypes() and",
      "getResidualizedPhenotypes()."))
  invisible(NULL)
}

#' (Deprecated) Load and preprocess regional multivariate data
#'
#' \strong{Deprecated.} Use \code{\link{MultiTaskQtlDataset}} (passing a
#' list of per-study \code{QtlDataset} objects).
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
loadRegionalMultivariateData <- function(...) {
  .Deprecated(new = "MultiTaskQtlDataset", package = "pecotmr",
    msg = paste(
      "loadRegionalMultivariateData() has been removed. Build per-study",
      "QtlDataset() objects and combine them with MultiTaskQtlDataset();",
      "the multivariate-Y join is now a pipeline-side concern (mvSuSiE",
      "and mr.mash wrappers form the joint Y matrix from the QtlDataset",
      "list at use-time)."))
  invisible(NULL)
}

#' (Deprecated) Load regional functional association data
#'
#' \strong{Deprecated.} Use \code{\link{QtlDataset}}.
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
loadRegionalFunctionalData <- function(...) {
  .Deprecated(new = "QtlDataset", package = "pecotmr",
    msg = paste(
      "loadRegionalFunctionalData() has been removed. Build a QtlDataset()",
      "directly; the previous `minMarkers` filter can be applied to the",
      "phenotype SummarizedExperiment list before constructor entry."))
  invisible(NULL)
}

# Function to remove gene name at the end of context name
#' @export
cleanContextNames <- function(context, gene) {
  # Remove gene name if it matches the last part of the context
  gene <- gene[order(-nchar(unique(gene)))]
  for (geneId in gene) {
    context <- gsub(paste0("_", geneId), "", context)
  }
  return(context)
}

#' (Deprecated) Load TWAS Weights from RDS Files
#'
#' \strong{Deprecated.} File-path loaders have been removed from pecotmr.
#' Construct a \code{TwasWeights} collection directly with the
#' \code{TwasWeights()} constructor (passing per-tuple vectors and a
#' \code{SimpleList} of \code{TwasWeightsEntry} payloads). If you need to
#' assemble a collection from on-disk inputs, read those files in your
#' own code and pass the resulting in-memory objects to the constructor.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
loadTwasWeights <- function(...) {
  .Deprecated(
    new = "TwasWeights",
    package = "pecotmr",
    msg = paste(
      "loadTwasWeights() has been removed. Construct TwasWeights",
      "collections directly with TwasWeights() and TwasWeightsEntry();",
      "file-path loaders are no longer part of pecotmr."))
  invisible(NULL)
}


#' Standardize GWAS summary statistics column names
#'
#' Uses MungeSumstats' comprehensive column name mapping to standardize
#' column names from various GWAS formats, then renames to pecotmr conventions.
#' Optionally applies an additional custom column mapping file.
#'
#' @param sumstats A data frame of summary statistics.
#' @param columnFilePath Optional file path to a custom column mapping file
#'   (format: standard_name:original_name, one per line). Applied after
#'   MungeSumstats standardization.
#' @param commentString Comment character in columnFilePath. Default is "#".
#' @return A data frame with standardized column names.
#' @export
standardiseSumstatsColumns <- function(sumstats, columnFilePath = NULL, commentString = "#") {
  # MungeSumstats standard names -> pecotmr conventions
  msToPecotmr <- c(
    CHR = "chrom", BP = "pos", SNP = "variant_id",
    BETA = "beta", SE = "se", Z = "z", P = "p",
    N = "n_sample", N_CAS = "n_case", N_CON = "n_control",
    FRQ = "maf"
  )
  # Make a copy to avoid in-place modification by MungeSumstats
  sumstatsCopy <- data.frame(sumstats, check.names = FALSE)

  # Read the explicit user column mapping first. User declarations are
  # AUTHORITATIVE: a column the user mapped (e.g. `af:effect_allele_frequency`)
  # must not be silently overridden by MungeSumstats (which would otherwise
  # absorb `effect_allele_frequency` into `FRQ` -> `maf` before the custom map
  # could run). We therefore shield each declared source column behind a unique
  # placeholder, let MungeSumstats standardize everything else, then restore the
  # declared columns to their requested standard names last.
  placeholders <- character(0)
  if (!is.null(columnFilePath)) {
    if (!file.exists(columnFilePath)) {
      stop("Column mapping file not found: ", columnFilePath)
    }
    columnData <- read.table(columnFilePath,
      header = FALSE, sep = ":",
      comment.char = if (is.null(commentString)) "" else commentString,
      stringsAsFactors = FALSE
    )
    colnames(columnData) <- c("standard", "original")
    for (i in seq_len(nrow(columnData))) {
      idx <- which(colnames(sumstatsCopy) == columnData$original[i])
      if (length(idx) > 0) {
        ph <- paste0(".pecotmr_decl_", i)
        colnames(sumstatsCopy)[idx] <- ph
        placeholders[[ph]] <- columnData$standard[i]
      }
    }
  }

  # Use MungeSumstats for comprehensive column standardization (shielded
  # declared columns pass through untouched as unmapped placeholders).
  sumstatsCopy <- standardise_header(
    sumstatsCopy, return_list = FALSE, uppercase_unmapped = FALSE
  )
  # Rename MungeSumstats standard names to pecotmr conventions
  for (msName in names(msToPecotmr)) {
    idx <- which(colnames(sumstatsCopy) == msName)
    if (length(idx) > 0) {
      colnames(sumstatsCopy)[idx] <- msToPecotmr[msName]
    }
  }
  # Restore user-declared columns to their requested standard names (last word).
  for (ph in names(placeholders)) {
    idx <- which(colnames(sumstatsCopy) == ph)
    if (length(idx) > 0) {
      colnames(sumstatsCopy)[idx] <- placeholders[[ph]]
    }
  }
  as.data.frame(sumstatsCopy)
}

#' (Deprecated) Load Summary Statistic Data
#'
#' \strong{Deprecated.} File-path summary-statistic loading has been
#' removed from pecotmr. Read your summary statistics in your own code
#' (see \code{vignette('rss-qc')} for examples using MungeSumstats for
#' text files, Rsamtools for tabix-indexed files, and VariantAnnotation
#' for VCFs), build a \code{GRanges} with the required mcols
#' (\code{SNP}, \code{A1}, \code{A2}, \code{Z}, \code{N}; plus optional
#' \code{MAF}, \code{INFO}, \code{BETA}, \code{SE}, \code{P}), and pass
#' it to the \code{\link{GwasSumStats}} constructor along with the
#' \code{ldSketch} and \code{genome}. Per-trait phenotype variance
#' \code{varY} (including the binary-trait OLS formula
#' \code{n / (n - 1) * phi * (1 - phi)}) is a constructor argument.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
loadRssData <- function(...) {
  .Deprecated(new = "GwasSumStats", package = "pecotmr",
    msg = paste(
      "loadRssData() has been removed. Build a GRanges of summary",
      "statistics in your own code (vignette('rss-qc') shows examples",
      "with MungeSumstats / Rsamtools / VariantAnnotation), then pass it",
      "to GwasSumStats()."))
  invisible(NULL)
}

#' (Deprecated) Load mixture regional data across multiple cohorts
#'
#' \strong{Deprecated.} Build per-study \code{\link{QtlDataset}} objects
#' for the individual-level cohorts and a \code{\link{QtlSumStats}} for
#' the summary-statistic cohorts, then combine them with
#' \code{\link{MultiTaskQtlDataset}}. The multivariate-Y join across
#' cohorts is now a pipeline-side concern (mvSuSiE / mr.mash wrappers
#' form the joint Y matrix at use-time).
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
loadMultitaskRegionalData <- function(...) {
  .Deprecated(new = "MultiTaskQtlDataset", package = "pecotmr",
    msg = paste(
      "loadMultitaskRegionalData() has been removed. Build per-study",
      "QtlDataset() objects for individual-level cohorts and a",
      "QtlSumStats() for summary-statistic cohorts, then combine with",
      "MultiTaskQtlDataset()."))
  invisible(NULL)
}

#' (Deprecated) Convert loaded regional data to individual-level inputs
#'
#' \strong{Deprecated.} The \code{RegionalData} class has been removed;
#' build individual-level inputs directly via the \code{\link{QtlDataset}}
#' constructor.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
regionDataToIndInput <- function(...) {
  .Deprecated(new = "QtlDataset", package = "pecotmr",
    msg = paste(
      "regionDataToIndInput() has been removed alongside RegionalData.",
      "Build individual-level inputs directly via the QtlDataset()",
      "constructor."))
  invisible(NULL)
}

#' (Deprecated) Convert loaded regional data to RSS inputs
#'
#' \strong{Deprecated.} The \code{RegionalData} / \code{QcResult} classes
#' have been removed; build RSS inputs directly via the
#' \code{\link{QtlSumStats}} or \code{\link{GwasSumStats}} constructors,
#' then run \code{\link{summaryStatsQc}}.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
regionDataToRssInput <- function(...) {
  .Deprecated(new = "QtlSumStats", package = "pecotmr",
    msg = paste(
      "regionDataToRssInput() has been removed alongside RegionalData.",
      "Build RSS inputs directly via QtlSumStats() / GwasSumStats() and",
      "run summaryStatsQc()."))
  invisible(NULL)
}

#' Load and filter tabular data with optional region subsetting
#'
#' This function loads summary statistics data from tabular files (TSV, TXT).
#' For compressed (.gz) and tabix-indexed files, it can subset data by genomic region.
#' Additionally, it can filter results by a specified target value in a designated column.
#'
#' @param filePath Path to the summary statistics file.
#' @param region Genomic region for subsetting tabix-indexed files. Format: chr:start-end (e.g., "9:10000-50000").
#' @param extractRegionName Value to filter for in the specified filter column.
#' @param regionNameCol Index of the column to apply the extract_region_name against.
#'
#' @return A dataframe containing the filtered summary statistics.
#'
#' @importFrom vroom vroom
#' @export
loadTsvRegion <- function(filePath, region = NULL, extractRegionName = NULL, regionNameCol = NULL) {
  sumstats <- NULL

  if (grepl("\\.gz$", filePath)) {
    if (!is.null(region)) {
      # Use Rsamtools to query the tabix-indexed file by region
      sumstats <- tryCatch({
        tbx <- TabixFile(filePath)
        parsed <- parseRegion(region)
        # Match chromosome naming convention in the tabix index
        chrom <- as.character(parsed$chrom)
        tbxSeqnames <- seqnamesTabix(tbx)
        if (any(grepl("^chr", tbxSeqnames))) {
          chrom <- paste0("chr", chrom)
        }
        gr <- GRanges(
          seqnames = chrom,
          ranges = IRanges(start = parsed$start, end = parsed$end)
        )
        lines <- scanTabix(tbx, param = gr)[[1]]
        if (length(lines) == 0) return(NULL)

        # Get header for column names
        hdr <- headerTabix(tbx)$header
        colNamesVec <- NULL
        if (length(hdr) > 0) {
          lastHdr <- hdr[length(hdr)]
          colNamesVec <- strsplit(sub("^#", "", lastHdr), "\t")[[1]]
        } else {
          headerCon <- gzfile(filePath, "rt")
          firstLine <- readLines(headerCon, n = 1)
          close(headerCon)
          firstFields <- strsplit(sub("^#", "", firstLine), "\t")[[1]]
          headerTokens <- c("chrom", "chr", "#chrom", "pos", "bp", "snp",
                            "variant_id", "a1", "a2", "beta", "se", "z",
                            "p", "pvalue")
          if (any(tolower(firstFields) %in% headerTokens)) {
            colNamesVec <- firstFields
          }
        }

        txt <- paste(lines, collapse = "\n")
        if (!is.null(colNamesVec)) {
          as.data.frame(vroom(I(txt), delim = "\t", col_names = colNamesVec,
                                     show_col_types = FALSE))
        } else {
          as.data.frame(vroom(I(txt), delim = "\t", col_names = TRUE,
                                     show_col_types = FALSE))
        }
      }, error = function(e) {
        stop("Data read error. Please make sure this gz file is tabix-indexed and the specified filter column exists.")
      })
    } else {
      # No region specified - read the whole gz file
      sumstats <- as.data.frame(vroom(filePath, show_col_types = FALSE))
    }
  } else {
    warning("Not a tabix-indexed gz file, loading the entire dataset.")
    sumstats <- as.data.frame(vroom(filePath, show_col_types = FALSE))
  }

  # Apply name-based filter if specified
  if (!is.null(sumstats) && !is.null(extractRegionName) && !is.null(regionNameCol)) {
    keepIndex <- which(str_detect(sumstats[[regionNameCol]], extractRegionName))
    sumstats <- sumstats[keepIndex, ]
  }

  return(sumstats)
}

#' Split loaded twas_weights_results into batches based on maximum memory usage
#'
#' @param twasWeightsResults List of loaded gene data by loadTwasWeights()
#' @param metaDataDf Dataframe containing gene metadata with region_id and TSS columns
#' @param maxMemoryPerBatch Maximum memory per batch in MB (default: 750)
#' @return List of batches, where each batch contains a subset of twas_weights_results
#' @export
batchLoadTwasWeights <- function(twasWeightsResults, metaDataDf, maxMemoryPerBatch = 750) {
  geneNames <- names(twasWeightsResults)
  if (length(geneNames) == 0) {
    message("No genes in twas_weights_results.")
    return(list())
  }

  geneMemoryDf <- data.frame(
    geneName = geneNames, memoryMb = sapply(geneNames, function(gene) {
      as.numeric(object.size(twasWeightsResults[[gene]])) / (1024^2) # Get object size in bytes and convert to MB
    })
  )

  # Merge with meta_data_df to get TSS information
  metaDataDf <- metaDataDf[!duplicated(metaDataDf[, c("region_id", "TSS")]), ]
  geneMemoryDf <- merge(geneMemoryDf, metaDataDf[, c("region_id", "TSS")],
    by.x = "geneName",
    by.y = "region_id", all.x = TRUE
  )
  geneMemoryDf <- geneMemoryDf[order(geneMemoryDf$TSS), ]

  # Check if we need to split into batches
  totalMemoryMb <- sum(geneMemoryDf$memoryMb)
  message("Total memory usage: ", round(totalMemoryMb, 2), " MB")
  if (totalMemoryMb <= maxMemoryPerBatch) {
    message("All genes fit within the memory limit. No need to split into batches.")
    return(list(allGenes = twasWeightsResults))
  }

  # Create batches by adding genes until we reach the memory limit
  batches <- list()
  currentBatchGenes <- character(0)
  currentBatchMemory <- 0
  batchIndex <- 1

  for (i in 1:nrow(geneMemoryDf)) {
    gene <- geneMemoryDf$geneName[i]
    geneMemory <- geneMemoryDf$memoryMb[i]
    # If a single gene exceeds the memory limit, include it in its own batch
    if (geneMemory > maxMemoryPerBatch) {
      batches[[paste0("batch_", batchIndex)]] <- twasWeightsResults[gene]
      batchIndex <- batchIndex + 1
      next
    }
    # If adding this gene would exceed the memory limit, start a new batch
    if (currentBatchMemory + geneMemory > maxMemoryPerBatch && length(currentBatchGenes) > 0) {
      batches[[paste0("batch_", batchIndex)]] <- twasWeightsResults[currentBatchGenes]
      currentBatchGenes <- character(0)
      currentBatchMemory <- 0
      batchIndex <- batchIndex + 1
    }
    currentBatchGenes <- c(currentBatchGenes, gene)
    currentBatchMemory <- currentBatchMemory + geneMemory
  }
  # Add the last batch if not empty
  if (length(currentBatchGenes) > 0) {
    batches[[batchIndex]] <- twasWeightsResults[currentBatchGenes]
  }
  message("Split into ", length(batches), " batches")
  names(batches) <- NULL
  return(batches)
}

# Function to filter a single credible set based on coverage and purity
#' @importFrom susieR susie_get_cs
#' @importFrom purrr map_lgl
#' @export
getFilterLbfIndex <- function(susieObj, coverage = 0.5, sizeFactor = 0.5) {
  susieObj$V <- NULL  # ensure no filtering by estimated prior

  # Get CS list with coverage
  csList <- susie_get_cs(susieObj, coverage = coverage, dedup = FALSE)

  # Total number of variants
  totalVariants <- ncol(susieObj$alpha)

  # Maximum allowed CS size to be considered 'concentrated'
  maxSize <- totalVariants * coverage * sizeFactor

  # Identify which CSs are 'concentrated enough'
  keepIdx <- map_lgl(csList$cs, ~ length(.x) < maxSize)

  # Extract the CS indices that pass the filter
  csIndex <- which(keepIdx) %>% names %>% gsub("L","", .) %>% as.numeric

  # Return filtered lbf_variable rows (one per CS)
  return(csIndex)
}
