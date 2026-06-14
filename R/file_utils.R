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
  if (!is.null(region)) {
    snpIdx <- .regionToSnpIdx(handle@snpInfo, region)
    if (length(snpIdx) == 0) {
      stop(NoSNPsError(paste("No SNPs found in the specified region", region)))
    }
  } else {
    snpIdx <- seq_len(nrow(handle@snpInfo))
  }

  # --- Extract genotypes (no mean imputation — callers handle missing) ---
  rse <- extractBlockGenotypes(handle, snpIdx, meanImpute = FALSE)
  # Convert RSE to samples x variants matrix for pecotmr convention
  X <- t(assay(rse, "dosage"))
  variantInfo <- .snpInfoToVariantInfo(
    handle@snpInfo[snpIdx, , drop = FALSE])

  # --- Attach allele frequency from .afreq sidecar (plink2 only) ---
  if (handle@format == "plink2") {
    afreq <- readAfreq(handle@path)
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

#' Load regional association data
#'
#' This function loads genotype, phenotype, and covariate data for a specific region and performs data preprocessing.
#'
#' @param genotype PLINK bed file containing genotype data.
#' @param phenotype A vector of phenotype file names.
#' @param covariate A vector of covariate file names corresponding to the phenotype file vector.
#' @param region A string of chr:start-end for the phenotype region.
#' @param conditions A vector of strings representing different conditions or groups.
#' @param maf_cutoff Minimum minor allele frequency (MAF) cutoff. Default is 0.
#' @param mac_cutoff Minimum minor allele count (MAC) cutoff. Default is 0.
#' @param xvar_cutoff Minimum variance cutoff. Default is 0.
#' @param imiss_cutoff Maximum individual missingness cutoff. Default is 0.
#' @param association_window A string of chr:start-end for the association analysis window (cis or trans). If not provided, all genotype data will be loaded.
#' @param extract_region_name A list of vectors of strings (e.g., gene ID ENSG00000269699) to subset the information when there are multiple regions available. Default is NULL.
#' @param region_name_col Column name containing the region name. Default is NULL.
#' @param keep_indel Logical indicating whether to keep insertions/deletions (INDELs). Default is TRUE.
#' @param keep_samples A vector of sample names to keep. Default is NULL.
#' @param phenotype_header Number of rows to skip at the beginning of the transposed phenotype file (default is 4 for chr, start, end, and ID).
#' @param scale_residuals Logical indicating whether to scale residuals. Default is FALSE.
#' @param tabix_header Logical indicating whether the tabix file has a header. Default is TRUE.
#'
#' @return A \code{RegionalData} S4 object. Per-condition residualized
#'   phenotypes, residualized genotypes, and their scaling factors are
#'   computed on demand via accessors (\code{getResidualX()},
#'   \code{getResidualY()}, \code{getResidualXScalar()},
#'   \code{getResidualYScalar()}, \code{getXVariance()}). Region metadata is
#'   available via \code{getChrom()} and \code{getGrange()}.
#'
#' @export
loadRegionalAssociationData <- function(genotype, # PLINK file
                                        phenotype, # a vector of phenotype file names
                                        covariate, # a vector of covariate file names corresponding to the phenotype file vector
                                        region, # a string of chr:start-end for phenotype region
                                        conditions, # a vector of strings
                                        mafCutoff = 0,
                                        macCutoff = 0,
                                        xvarCutoff = 0,
                                        imissCutoff = 0,
                                        associationWindow = NULL,
                                        extractRegionName = NULL,
                                        regionNameCol = NULL,
                                        keepIndel = TRUE,
                                        keepSamples = NULL,
                                        keepVariants = NULL,
                                        phenotypeHeader = 4, # skip first 4 rows of transposed phenotype for chr, start, end and ID
                                        scaleResiduals = FALSE,
                                        tabixHeader = TRUE) {
  ## Load genotype
  geno <- loadGenotypeRegion(genotype, associationWindow, keepIndel, keepVariantsPath = keepVariants)
  ## Load phenotype and covariates and perform some pre-processing
  covar <- loadCovariateData(covariate)
  pheno <- loadPhenotypeData(phenotype, region, extractRegionName = extractRegionName, regionNameCol = regionNameCol, tabixHeader = tabixHeader)
  # Align covariates and conditions with phenotypes after filtering
  # loadPhenotypeData removes empty phenotypes and stores which indices survived
  keptIdx <- attr(pheno, "kept_indices")
  if (!is.null(keptIdx) && length(keptIdx) < length(covar)) {
    covar <- covar[keptIdx]
    if (!is.null(conditions)) conditions <- conditions[keptIdx]
  }
  ### including Y ( cov ) and specific X and covar match, filter X variants based on the overlapped samples.
  dataList <- prepareDataList(geno, pheno, covar, imissCutoff,
    mafCutoff, macCutoff, xvarCutoff,
    phenotypeHeader = phenotypeHeader, keepSamples = keepSamples
  )
  mafList <- setNames(lapply(dataList$X, function(x) apply(x, 2, computeMaf)), conditions)
  ## Get residue Y for each of condition and its mean and sd
  dataList <- addYResiduals(dataList, conditions, scaleResiduals)
  ## Get residue X for each of condition and its mean and sd
  dataList <- addXResiduals(dataList, scaleResiduals)
  # Get X matrix for union of samples.
  # Short-circuit when there's only one condition: the per-condition X computed in
  # prepareDataList already operates on the same sample set (the single condition's
  # common_complete_samples, which is itself a subset of rownames(geno)) and applies
  # the same MAF/imiss/var thresholds with the same MAC cutoff scaling, so the union
  # X is bit-equivalent to data_list$X[[1]]. Skipping the redundant filter_X call saves
  # work and avoids a duplicate "N out of M total variants dropped" log line.
  if (length(dataList$X) == 1) {
    X <- dataList$X[[1]]
    variants <- str_split_fixed(colnames(X), ":", 3)
    message(paste0(
      "Dimension of input genotype data is ", nrow(X), " rows and ",
      ncol(X), " columns for genomic region of ",
      variants[1, 1], ":", min(as.integer(variants[, 2])), "-",
      max(as.integer(variants[, 2]))
    ))
  } else {
    X <- prepareXMatrix(geno, dataList, imissCutoff, mafCutoff, macCutoff, xvarCutoff)
  }
  parsedRegion <- if (!is.null(region)) parseRegion(region) else NULL
  regionGr <- if (!is.null(parsedRegion)) {
    GRanges(
      seqnames = paste0("chr", parsedRegion$chrom),
      ranges = IRanges(
        start = as.integer(parsedRegion$start),
        end = as.integer(parsedRegion$end)
      )
    )
  } else NULL

  phenoList <- dataList$Y
  covarList <- dataList$covar
  if (!is.null(conditions)) {
    names(phenoList) <- conditions
    names(covarList) <- conditions
  }
  RegionalData(
    genotypeMatrix = X,
    phenotypes = phenoList,
    covariates = covarList,
    scaleResiduals = scaleResiduals,
    maf = mafList,
    region = regionGr,
    droppedSamples = list(
      X = dataList$dropped_samples_X,
      Y = dataList$dropped_samples_Y,
      covar = dataList$dropped_samples_covar
    ),
    coordinates = if (!is.null(region)) extractPhenotypeCoordinates(pheno) else NULL
  )
}

#' Load Regional Univariate Association Data
#'
#' Loads regional association data for univariate analysis. Returns a
#' \code{RegionalData} S4 object; derived quantities (residuals, scalars,
#' per-variant variance) are computed lazily via accessors
#' (\code{getResidualX}, \code{getResidualY}, \code{getResidualXScalar},
#' \code{getResidualYScalar}, \code{getXVariance}, \code{getChrom},
#' \code{getGrange}).
#'
#' @return A \code{RegionalData} object.
#' @export
loadRegionalUnivariateData <- function(...) {
  loadRegionalAssociationData(...)
}

#' Load Regional Data for Regression Modeling
#'
#' Loads regional association data formatted for regression modeling.
#' Returns a \code{RegionalData} S4 object; the per-condition \code{X_data}
#' previously returned in a list is available as
#' \code{getResidualX(rd, i)} (residualized) or by subsetting
#' \code{rd@@genotypeMatrix} by condition rownames.
#'
#' @return A \code{RegionalData} object.
#' @export
loadRegionalRegressionData <- function(...) {
  loadRegionalAssociationData(...)
}

# return matrix of R conditions, with column names being the names of the conditions (phenotypes) and row names being sample names. Even for one condition it has to be a matrix with just one column.
#' @noRd
phenoListToMat <- function(dataList) {
  allRowNames <- unique(unlist(lapply(dataList$residual_Y, rownames)))
  # Step 2: Align matrices and fill with NA where necessary
  alignedMats <- lapply(dataList$residual_Y, function(mat) {
    ### change the ncol of each matrix
    expandedMat <- matrix(NA, nrow = length(allRowNames), ncol = ncol(mat), dimnames = list(allRowNames, colnames(mat)))
    commonRows <- intersect(rownames(mat), allRowNames)
    expandedMat[commonRows, ] <- mat[commonRows, ]
    return(expandedMat)
  })
  YResidMatrix <- do.call(cbind, alignedMats)
  if (!is.null(names(dataList$residual_Y))) {
    colnames(YResidMatrix) <- names(dataList$residual_Y)
  }
  dataList$residual_Y <- YResidMatrix
  return(dataList)
}

#' Load and Preprocess Regional Multivariate Data
#'
#' Loads regional association data and packages it for multivariate modeling.
#' Phenotypes across conditions are joined into a single multivariate matrix
#' (samples x conditions). When \code{matrix_y_min_complete} is supplied,
#' samples with fewer than that many non-missing condition values are dropped.
#' Per-variant MAF and variance are computed on the (post-filter) genotype
#' matrix and exposed via \code{getMAF()} / \code{getXVariance()} on the
#' returned object.
#'
#' @return A \code{MultivariateRegionalData} object.
#' @export
loadRegionalMultivariateData <- function(matrixYMinComplete = NULL,
                                         ...) {
  rd <- loadRegionalAssociationData(...)
  nCond <- length(rd@phenotypes)
  residualYList <- lapply(seq_len(nCond), function(i) getResidualY(rd, i))
  names(residualYList) <- names(rd@phenotypes)
  residualYScalarList <- lapply(seq_len(nCond), function(i) getResidualYScalar(rd, i))
  dat <- list(residual_Y = residualYList)
  dat <- phenoListToMat(dat)

  X <- rd@genotypeMatrix
  YScalar <- unlist(residualYScalarList)
  droppedSample <- rd@droppedSamples
  regionGr <- rd@region
  Y <- dat$residual_Y

  if (!is.null(matrixYMinComplete)) {
    filt <- filterY(Y, matrixYMinComplete)
    if (length(filt$rm_rows) > 0) {
      X <- X[-filt$rm_rows, , drop = FALSE]
      Y <- filt$Y
      droppedSample <- rownames(dat$residual_Y)[filt$rm_rows]
    }
  }

  MultivariateRegionalData(
    genotypeMatrix = X,
    Y = as.matrix(Y),
    scaling = YScalar,
    droppedSamples = droppedSample,
    region = regionGr,
    coordinates = rd@coordinates
  )
}

#' Load Regional Functional Association Data
#'
#' Loads precomputed regional functional association data. Returns a
#' \code{RegionalData} S4 object; derived quantities are computed lazily
#' via accessors. When \code{min_markers} is supplied, conditions whose
#' \code{Y_coordinates} have fewer than \code{min_markers} rows are
#' dropped from the returned \code{RegionalData}.
#'
#' @param minMarkers Minimum number of phenotype markers required for a study.
#'   If \code{NULL}, no marker-count filtering is applied.
#' @return A \code{RegionalData} object.
#' @export
loadRegionalFunctionalData <- function(..., minMarkers = NULL) {
  rd <- loadRegionalAssociationData(...)
  if (!is.null(minMarkers)) rd <- .filterRegionalDataByMarkerCount(rd, minMarkers)
  rd
}

# Subset per-condition slots of a RegionalData by the marker counts in
# Y_coordinates. The genotype_matrix, region, and dropped_samples are
# preserved (those are not per-condition or are panel-wide).
.filterRegionalDataByMarkerCount <- function(rd, minMarkers) {
  if (is.null(rd@coordinates)) return(rd)
  keep <- vapply(rd@coordinates, function(x) nrow(x) >= minMarkers, logical(1))
  if (all(keep)) return(rd)
  RegionalData(
    genotypeMatrix = rd@genotypeMatrix,
    phenotypes = rd@phenotypes[keep],
    covariates = rd@covariates[keep],
    scaleResiduals = rd@scaleResiduals,
    maf = if (length(rd@maf) == length(keep)) rd@maf[keep] else rd@maf,
    region = rd@region,
    droppedSamples = rd@droppedSamples,
    coordinates = rd@coordinates[keep]
  )
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

#' Load, Validate, and Consolidate TWAS Weights from Multiple RDS Files
#'
#' This function loads TWAS weight data from multiple RDS files, checks for the presence
#' of specified region and condition. If variable_name_obj is provided, it aligns and
#' consolidates weight matrices based on the object's variant names, filling missing data
#' with zeros. If variable_name_obj is NULL, it checks that all files have the same row
#' numbers for the condition and consolidates weights accordingly.
#'
#' @param weight_db_file weight_db_files Vector of file paths for RDS files containing TWAS weights..
#' Each element organized as region/condition/weights
#' @param condition The specific condition to be checked and consolidated across all files.
#' @param variable_name_obj The name of the variable/object to fetch from each file, if not NULL.
#' @return A consolidated list of weights for the specified condition and a list of SuSiE results.
#' @examples
#' # Example usage (replace with actual file paths, condition, region, and variable_name_obj):
#' weight_db_files <- c("path/to/file1.rds", "path/to/file2.rds")
#' condition <- "example_condition"
#' region <- "example_region"
#' variable_name_obj <- "example_variable" # or NULL for standard processing
#' consolidated_weights <- loadTwasWeights(weight_db_files, condition, region, variable_name_obj)
#' print(consolidated_weights)
#' @export
loadTwasWeights <- function(weightDbFiles, conditions = NULL,
                            variableNameObj = c("preset_variants_result", "variantNames"),
                            susieObj = c("preset_variants_result", "susie_result_trimmed"),
                            twasWeightsTable = "twas_weights") {
  ## Internal function to load and validate data from RDS files
  loadAndValidateData <- function(weightDbFiles, conditions, variableNameObj) {
    allData <- do.call(c, lapply(unname(weightDbFiles), function(rdsFile) {
      # Validate file before loading
      if (!file.exists(rdsFile)) {
        warning(paste0("Skipping weight file '", rdsFile, "': file does not exist."))
        return(NULL)
      }
      if (file.size(rdsFile) <= 200) {
        warning(paste0("Skipping weight file '", rdsFile, "': file too small (", file.size(rdsFile), " bytes), likely empty or corrupt."))
        return(NULL)
      }
      db <- tryCatch(readRDS(rdsFile), error = function(e) {
        warning(paste0("Skipping weight file '", rdsFile, "': failed to read RDS — ", conditionMessage(e)))
        return(NULL)
      })
      if (!is.list(db) || length(db) == 0) {
        warning(paste0("Skipping weight file '", rdsFile, "': unexpected structure (not a non-empty list)."))
        return(NULL)
      }
      gene <- names(db)
      # Filter by conditions if specified
      if (!is.null(conditions)) {
        # Split contexts if specified and trim whitespace, cen handle single or multiple conditions
        conditions <- trimws(strsplit(conditions, ",")[[1]])

        # Filter the gene's data to only include specified context layers
        if (length(gene) == 1 && gene != "mnm_rs") { # Need check
          availableContexts <- names(db[[gene]])
          matchingContexts <- availableContexts[availableContexts %in% paste0(conditions, "_", gene)]
          if (length(matchingContexts) == 0) {
            warning(paste0("No matching context layers found in ", rdsFile, ". Skipping this file."))
            return(NULL)
          }

          db[[gene]] <- db[[gene]][matchingContexts]
        }
      } else {
        # Set default for 'conditions' if they are not specified
        conditions <- names(db[[gene]])
      }
      if (any(unique(names(findData(db, c(3, "twas_weights")))) %in% c("mrmash_weights", "mvsusie_weights"))) {
        names(db[[1]]) <- cleanContextNames(names(db[[1]]), gene = gene)
        db <- list(mnm_rs = db[[1]])
      } else {
        # Check if region from all RDS files are the same
        if (length(gene) != 1) {
          stop("More than one region provided in the RDS file. ")
        } else {
          names(db[[gene]]) <- cleanContextNames(names(db[[gene]]), gene = gene)
        }
      }
      return(db)
    }))
    # Remove NULL entries (from files that had no matching context layers)
    allData <- allData[!sapply(allData, is.null)]

    if (length(allData) == 0) {
      stop("No data loaded. Check that conditions match available context layers in the RDS files.")
    }
    # Combine the lists with the same region name
    gene <- unique(names(allData)[!names(allData) %in% "mnm_rs"])
    if (length(gene) > 1) stop("More than one region of twas weights data provided. ")
    combinedAllData <- lapply(split(allData, names(allData)), function(lst) {
      if (length(lst) > 1) {
        lst <- do.call(c, unname(lst))
      }
      if (isTRUE(names(lst) == "mnm_rs")) lst <- lst[[1]]
      if (gene %in% names(lst)) lst <- do.call(c, lapply(unname(lst), function(x) x))
      return(lst)
    })

    # merge univariate and multivariate results for same gene-context pair
    if ("mnm_rs" %in% names(combinedAllData)) {
      # gene <- names(combinedAllData)[!names(combinedAllData) %in% "mnm_rs"]
      overlContexts <- names(combinedAllData[["mnm_rs"]])[names(combinedAllData[["mnm_rs"]]) %in% names(combinedAllData[[gene]])]
      multiVariants <- unique(findData(combinedAllData$mnm_rs, c(2, variableNameObj)))
      for (context in overlContexts) {
        uniVariants <- getNestedElement(combinedAllData[[gene]][[context]], variableNameObj)
        # Harmonize chr prefix convention between multivariate and univariate variant IDs
        chrMatched <- ensureChrMatch(multiVariants, uniVariants)
        multiVariantsH <- chrMatched$idsA
        uniVariantsH <- chrMatched$idsB
        multiWeights <- setNames(rep(0, length(uniVariantsH)), uniVariantsH)
        multiWeights <- lapply(combinedAllData[["mnm_rs"]][[context]]$twasWeights, function(weightList) {
          alignedWeights <- setNames(rep(0, length(uniVariantsH)), uniVariantsH)
          weightVals <- unlist(weightList)
          names(weightVals) <- ensureChrMatch(names(weightVals), uniVariantsH)$idsA
          methodWeightVariants <- names(weightVals)
          overlapVariants <- methodWeightVariants[methodWeightVariants %in% multiVariantsH[multiVariantsH %in% uniVariantsH]]
          alignedWeights[overlapVariants] <- weightVals[overlapVariants]
          alignedWeights <- as.matrix(alignedWeights)
        })
        combinedAllData[[gene]][[context]]$twasWeights <- c(combinedAllData[[gene]][[context]]$twasWeights, multiWeights)
        combinedAllData[[gene]][[context]]$twasCvResult$performance <- c(
          combinedAllData[[gene]][[context]]$twasCvResult$performance,
          combinedAllData[["mnm_rs"]][[context]]$twasCvResult$performance
        )
      }
      combinedAllData[["mnm_rs"]] <- NULL
    }
    if (gene %in% names(combinedAllData)) combinedAllData <- do.call(c, unname(combinedAllData))
    if (gene %in% names(combinedAllData)) combinedAllData <- combinedAllData[[1]]

    # ## Check if the specified condition and variable_name_obj are available in all files
    # if (!all(conditions %in% names(combinedAllData))) {
    #   stop("The specified condition is not available in all RDS files.")
    # }
    return(combinedAllData)
  }

  # Internal function to align and merge weight matrices
  alignAndMerge <- function(weightsList, variableObjs) {
    if (length(weightsList) != length(variableObjs)) {
      stop("The length of the weights_list and variable_objs must be the same.")
    }
    # Loop through each weight matrix and assign variant names as rownames
    for (i in seq_along(weightsList)) {
      # Ensure dimensions match
      if (nrow(weightsList[[i]]) != length(variableObjs[[i]])) {
        stop(paste("Number of rows in weights_list[[", i, "]] does not match the length of variable_objs[[", i, "]]", sep = ""))
      }
      # Apply variant names to the row names of the weight matrix
      rownames(weightsList[[i]]) <- variableObjs[[i]]
    }
    return(weightsList)
  }

  # Internal function to consolidate weights for given condition
  consolidateWeightsList <- function(combinedAllData, conditions, variableNameObj, twasWeightsTable) {
    combinedWeightsByCondition <- lapply(conditions, function(condition) {
      tempList <- getNestedElement(combinedAllData, c(condition, twasWeightsTable))
      sapply(tempList, cbind)
    })
    names(combinedWeightsByCondition) <- conditions
    if (is.null(variableNameObj)) {
      # Standard processing: Check for identical row numbers and consolidate
      rowNumbers <- sapply(combinedWeightsByCondition, function(data) nrow(data))
      if (length(unique(rowNumbers)) > 1) {
        stop("Not all files have the same number of rows for the specified condition.")
      }
      weights <- combinedWeightsByCondition
    } else {
      # Processing with variable_name_obj: Align and merge data, fill missing with zeros
      variableObjs <- lapply(conditions, function(condition) {
        getNestedElement(combinedAllData, c(condition, variableNameObj))
      })
      weights <- alignAndMerge(combinedWeightsByCondition, variableObjs)
    }
    names(weights) <- conditions
    return(weights)
  }

  ## Load, validate, and consolidate data
  try(
    {
      combinedAllData <- loadAndValidateData(weightDbFiles, conditions, variableNameObj)
      if (is.null(combinedAllData)) {
        return(NULL)
      }
      # update condition in case of merging rds files
      conditions <- names(combinedAllData)
      weights <- consolidateWeightsList(combinedAllData, conditions, variableNameObj, twasWeightsTable)
      combinedSusieResult <- lapply(combinedAllData, function(context) getNestedElement(context, susieObj))
      performanceTables <- lapply(conditions, function(condition) {
        getNestedElement(combinedAllData, c(condition, "twasCvResult", "performance"))
      })
      names(performanceTables) <- conditions
      # Extract variant_ids from weight matrices (union across all contexts)
      allVariantIds <- Reduce(union, lapply(weights, function(w) {
        if (is.matrix(w)) rownames(w) else names(w)
      }))
      if (is.null(allVariantIds)) allVariantIds <- character(0)
      # Pad weight matrices to common variant set
      if (length(allVariantIds) > 0) {
        weights <- lapply(weights, function(w) {
          if (!is.matrix(w)) return(w)
          missing <- setdiff(allVariantIds, rownames(w))
          if (length(missing) > 0) {
            pad <- matrix(0, nrow = length(missing), ncol = ncol(w),
                          dimnames = list(missing, colnames(w)))
            w <- rbind(w, pad)[allVariantIds, , drop = FALSE]
          }
          w
        })
      }
      return(TwasWeights(
        weights = weights,
        variantIds = allVariantIds,
        fits = combinedSusieResult,
        cvPerformance = performanceTables
      ))
    },
    silent = FALSE
  )
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

#' Load summary statistic data
#'
#' This function formats the input summary statistics dataframe with uniform column names
#' to fit into the SuSiE pipeline. Column standardization is performed via
#' MungeSumstats::standardise_header(), with an optional custom column mapping file
#' for additional non-standard names.
#' Additionally, it extracts sample size, case number, control number, and
#' optionally variance of Y for observed-scale binary OLS summary statistics.
#' Missing values in n_sample, n_case, and n_control are backfilled with median values.
#'
#' @param sumstatPath File path to the summary statistics.
#' @param columnFilePath Optional file path to a custom column mapping file for
#'   non-standard column names not recognized by MungeSumstats.
#' @param n_sample User-specified sample size. If unknown, set as 0 to retrieve from the sumstat file.
#' @param n_case User-specified number of cases.
#' @param n_control User-specified number of controls.
#' @param binary_trait_model How to treat case-control sample counts. The
#'   default \code{"rss"} uses counts only to infer \code{n} and leaves
#'   \code{var_y = NULL}, so \code{susieRss()} uses its z-score RSS interface
#'   on the standardized phenotype scale. Use \code{"ols"} only when
#'   \code{beta} and \code{se} come from ordinary least squares on a 0/1
#'   phenotype and the full \code{bhat/shat/var_y} sufficient-statistic
#'   interface is desired. In that case, if \code{phi = n_case / n}, centering gives
#'   \code{sum((y - phi)^2) = n * phi * (1 - phi)}, so the \code{susieR}
#'   \code{var_y = sum(y^2) / (n - 1)} input is
#'   \code{n / (n - 1) * phi * (1 - phi)}.
#' @param region The region where tabix use to subset the input dataset.
#' @param extract_region_name User-specified gene/phenotype name used to further subset the phenotype data.
#' @param region_name_col Filter this specific column for the extract_region_name.
#' @param comment_string Comment sign in the column_mapping file, default is #
#' @return A list of rss_input, including the column-name-formatted summary statistics,
#' sample size (n), and var_y.
#'
#' @importFrom dplyr mutate group_by summarise
#' @importFrom magrittr %>%
#' @export
loadRssData <- function(sumstatPath, columnFilePath = NULL, nSample = 0, nCase = 0, nControl = 0, region = NULL,
                        extractRegionName = NULL, regionNameCol = NULL, commentString = "#",
                        binaryTraitModel = c("rss", "ols")) {
  binaryTraitModel <- match.arg(binaryTraitModel)
  nSample <- if (length(nSample) == 1L && is.na(nSample)) 0 else nSample
  nCase <- if (length(nCase) == 1L && is.na(nCase)) 0 else nCase
  nControl <- if (length(nControl) == 1L && is.na(nControl)) 0 else nControl
  # Validate input files exist
  if (!file.exists(sumstatPath)) {
    stop("Summary statistics file not found: ", sumstatPath)
  }
  if (!is.null(columnFilePath) && !file.exists(columnFilePath)) {
    stop("Column mapping file not found: ", columnFilePath)
  }
  varY <- NULL
  sumstats <- loadTsvRegion(filePath = sumstatPath, region = region, extractRegionName = extractRegionName, regionNameCol = regionNameCol)

  # To keep a log message
  nVariants <- if (is.null(sumstats)) 0L else nrow(sumstats)
  if (length(nVariants) == 0 || is.na(nVariants)) {
    nVariants <- 0L
  }
  if (nVariants == 0) {
    message(paste0("No variants in region ", region, "."))
    if (is.null(sumstats)) {
      sumstats <- data.frame()
    }
    return(list(sumstats = sumstats, n = NULL, varY = NULL))
  } else {
    message(paste0("Region ", region, " include ", nVariants, " in input sumstats."))
  }

  # Standardize column names via MungeSumstats + optional custom mapping
  sumstats <- standardiseSumstatsColumns(sumstats, columnFilePath, commentString)

  # ---- Effect-allele frequency (af) propagation -------------------------------
  # `af` is the frequency of the effect allele / A1 and is exported (after
  # harmonization) as top_loci$af. It becomes available ONLY through an explicit
  # column-file mapping to the standard name `af` (MungeSumstats never emits
  # `af`; it maps FRQ -> `maf`, which stays an internal QC quantity). Ambiguous
  # frequency headers therefore never silently become `af`. The effect allele
  # must also be resolvable (an A1 column or an allele-bearing variant id), or
  # the declared frequency cannot be tied to a direction and is not exported.
  afDeclared <- "af" %in% colnames(sumstats)
  hasEffectAllele <- "A1" %in% colnames(sumstats) ||
    any(c("variant_id", "variant") %in% colnames(sumstats))
  if (afDeclared && hasEffectAllele) {
    sumstats$af <- suppressWarnings(as.numeric(sumstats$af))
    if (all(is.na(sumstats$af))) {
      warning("Effect-allele frequency column 'af' was declared but its values ",
              "are missing/unusable; top_loci$af will be NA and MAF filtering ",
              "will be skipped.")
    }
  } else {
    if (afDeclared && !hasEffectAllele) {
      warning("Effect-allele frequency 'af' was declared but no effect allele ",
              "(A1 / allele-bearing variant id) is available to tie it to a ",
              "direction; it will not be exported. top_loci$af will be NA and ",
              "MAF filtering will be skipped.")
    } else {
      warning("Effect-allele frequency (af) was not declared in the column ",
              "file; top_loci$af will be NA and MAF filtering will be skipped. ",
              "Generic frequency headers (FRQ/AF/allele_frequency) are kept as ",
              "internal MAF only and are never exported as af.")
    }
    sumstats$af <- NA_real_
  }
  # ----------------------------------------------------------------------------

  hasObservedBetaSe <- all(c("beta", "se") %in% colnames(sumstats))
  if (binaryTraitModel == "ols" && !hasObservedBetaSe) {
    stop("binaryTraitModel = 'ols' requires observed beta and se columns ",
         "from ordinary least squares on a 0/1 phenotype; z-only summary ",
         "statistics should use binaryTraitModel = 'rss'.")
  }
  if (!"z" %in% colnames(sumstats) && all(c("beta", "se") %in%
    colnames(sumstats))) {
    sumstats$z <- sumstats$beta / sumstats$se
  }
  if (!"beta" %in% colnames(sumstats) && "z" %in% colnames(sumstats)) {
    sumstats$beta <- sumstats$z
    sumstats$se <- 1
    attr(sumstats, "pecotmr_beta_se_from_z") <- TRUE
  } else {
    attr(sumstats, "pecotmr_beta_se_from_z") <- FALSE
  }
  for (col in c("n_sample", "n_case", "n_control")) {
    if (col %in% colnames(sumstats)) {
      sumstats[[col]][is.na(sumstats[[col]])] <- median(sumstats[[col]],
        na.rm = TRUE
      )
    }
  }
  binaryVarY <- function(phi, n) {
    if (length(phi) != 1 || is.na(phi) || !is.finite(phi) ||
        phi <= 0 || phi >= 1) {
      stop("Invalid case fraction for binaryTraitModel = 'ols': ", phi,
           ". Expected 0 < nCase / n < 1.")
    }
    if (is.null(n) || length(n) != 1 || is.na(n) || !is.finite(n) || n <= 1) {
      stop("Invalid sample size for binaryTraitModel = 'ols': ", n,
           ". Expected n > 1.")
    }
    n / (n - 1) * phi * (1 - phi)
  }
  if (nSample != 0 && (nCase + nControl) != 0) {
    stop("Please provide sample size, or case number with control number, but not both")
  } else if (nSample != 0) {
    n <- nSample
  } else if ((nCase + nControl) != 0) {
    n <- nCase + nControl
    if (binaryTraitModel == "ols") {
      phi <- nCase / n
      varY <- binaryVarY(phi, n)
    }
  } else {
    if ("n_sample" %in% colnames(sumstats)) {
      n <- median(sumstats$n_sample)
    } else if (all(c("n_case", "n_control") %in% colnames(sumstats))) {
      nByVariant <- sumstats$n_case + sumstats$n_control
      n <- median(nByVariant)
      if (binaryTraitModel == "ols") {
        phi <- median(sumstats$n_case / nByVariant)
        varY <- binaryVarY(phi, n)
      }
    } else {
      warning("Sample size could not be determined from the summary statistics.")
      n <- NULL
    }
  }
  # Validate determined sample size
  if (!is.null(n)) {
    if (length(n) != 1) {
      stop("Sample size must be a single value, got length ", length(n), ".")
    }
    if (is.na(n) || !is.finite(n) || n <= 0) {
      stop("Invalid sample size determined: ", n,
           ". Sample size must be a positive finite number.",
           "\n  Hint: check n_sample, n_case, n_control parameters or the ",
           "n_sample/n_case/n_control columns in your summary statistics.")
    }
  }
  return(list(sumstats = sumstats, n = n, varY = varY))
}


#' This function loads a mixture data sets for a specific region, including individual-level data (genotype, phenotype, covariate data)
#' or summary statistics (sumstats, LD). Run \code{loadRegionalUnivariateData} and \code{loadRssData} multiple times for different datasets
#'
#' @section Loading individual level data from multiple corhorts
#' @param region A string of chr:start-end for the phenotype region.
#' @param genotype_list a vector of PLINK bed file containing genotype data.
#' @param phenotype_list A vector of phenotype file names.
#' @param covariate_list A vector of covariate file names corresponding to the phenotype file vector.
#' @param conditions_list_individual A vector of strings representing different conditions or groups.
#' @param match_geno_pheno A vector of index of phentoypes matched to genotype if mulitple genotype PLINK files
#' @param maf_cutoff Minimum minor allele frequency (MAF) cutoff. Default is 0.
#' @param mac_cutoff Minimum minor allele count (MAC) cutoff. Default is 0.
#' @param xvar_cutoff Minimum variance cutoff. Default is 0.
#' @param imiss_cutoff Maximum individual missingness cutoff. Default is 0.
#' @param association_window A string of chr:start-end for the association analysis window (cis or trans). If not provided, all genotype data will be loaded.
#' @param extract_region_name A list of vectors of strings (e.g., gene ID ENSG00000269699) to subset the information when there are multiple regions available. Default is NULL.
#' @param region_name_col Column name containing the region name. Default is NULL.
#' @param keep_indel Logical indicating whether to keep insertions/deletions (INDELs). Default is TRUE.
#' @param keep_samples A vector of sample names to keep. Default is NULL.
#' @param phenotype_header Number of rows to skip at the beginning of the transposed phenotype file (default is 4 for chr, start, end, and ID).
#' @param scale_residuals Logical indicating whether to scale residuals. Default is FALSE.
#' @param tabix_header Logical indicating whether the tabix file has a header. Default is TRUE.
#'
#' @section Loading summary statistics from multiple corhorts or data set
#' @param sumstat_path_list A vector of file path to the summary statistics.
#' @param column_file_path_list A vector of file path to the column file for mapping.
#' @param LD_meta_file_path_list A vector of path of LD_metadata, LD_metadata is a data frame specifying LD blocks with columns "chrom", "start", "end", and "path". "start" and "end" denote the positions of LD blocks. "path" is the path of each LD block, optionally including bim file paths.
#' @param match_LD_sumstat A vector of index of sumstat matched to LD if mulitple sumstat files
#' @param conditions_list_sumstat A vector of strings representing different sumstats.
#' @param n_samples User-specified sample size. If unknown, set as 0 to retrieve from the sumstat file.
#' @param n_cases User-specified number of cases.
#' @param n_controls User-specified number of controls.
#' @param binary_trait_model How to treat case-control sample counts for summary
#'   statistics. Passed to \code{load_rss_data()}; the default \code{"rss"}
#'   keeps \code{var_y = NULL}, while \code{"ols"} computes the observed-scale
#'   OLS \code{var_y} from \code{n_case / n} when true OLS \code{beta/se}
#'   columns are available.
#' @param region The region where tabix use to subset the input dataset.
#' @param extract_sumstats_region_name User-specified gene/phenotype name used to further subset the phenotype data.
#' @param sumstats_region_name_col Filter this specific column for the extract_sumstats_region_name.
#' @param comment_string comment sign in the column_mapping file, default is #
#' @param extract_coordinates Optional data frame with columns "chrom" and "pos" for specific coordinates extraction.
#'
#' @return A list containing the individualData and sumstatData:
#' individualData contains the following components if exist
#' \itemize{
#'   \item residual_Y: A list of residualized phenotype values (either a vector or a matrix).
#'   \item residual_X: A list of residualized genotype matrices for each condition.
#'   \item residual_Y_scalar: Scaling factor for residualized phenotype values.
#'   \item residual_X_scalar: Scaling factor for residualized genotype values.
#'   \item dropped_sample: A list of dropped samples for X, Y, and covariates.
#'   \item covar: Covariate data.
#'   \item Y: Original phenotype data.
#'   \item X_data: Original genotype data.
#'   \item X: Filtered genotype matrix.
#'   \item maf: Minor allele frequency (MAF) for each variant.
#'   \item chrom: Chromosome of the region.
#'   \item grange: Genomic range of the region (start and end positions).
#'   \item Y_coordinates: Phenotype coordinates if a region is specified.
#' }
#' sumstatData contains the following components if exist
#' \itemize{
#'   \item sumstats: A list of summary statistics for the matched ldInfo, each sublist contains sumstats, n, var_y from \code{loadRssData}.
#'   \item ldInfo: A list of \code{LdData} S4 objects (one per LD reference), as returned by \code{load_LD_matrix}.
#' }
#'
#' @export
loadMultitaskRegionalData <- function(region, # a string of chr:start-end for phenotype region
                                      genotypeList = NULL, # PLINK file
                                      phenotypeList = NULL, # a vector of phenotype file names
                                      covariateList = NULL, # a vector of covariate file names corresponding to the phenotype file vector
                                      conditionsListIndividual = NULL, # a vector of strings
                                      matchGenoPheno = NULL, # a vector of index of phentoypes matched to genotype if mulitple genotype files
                                      mafCutoff = 0,
                                      macCutoff = 0,
                                      xvarCutoff = 0,
                                      imissCutoff = 0,
                                      associationWindow = NULL,
                                      extractRegionName = NULL,
                                      regionNameCol = NULL,
                                      keepIndel = TRUE,
                                      keepSamples = NULL,
                                      keepVariants = NULL,
                                      phenotypeHeader = 4, # skip first 4 rows of transposed phenotype for chr, start, end and ID
                                      scaleResiduals = FALSE,
                                      tabixHeader = TRUE,
                                      # sumstat if need
                                      sumstatPathList = NULL,
                                      columnFilePathList = NULL,
                                      ldMetaFilePathList = NULL,
                                      matchLdSumstat = NULL, # a vector of index of sumstat matched to LD if mulitple sumstat files
                                      conditionsListSumstat = NULL,
                                      nSamples = 0,
                                      nCases = 0,
                                      nControls = 0,
                                      binaryTraitModel = c("rss", "ols"),
                                      extractSumstatsRegionName = NULL,
                                      sumstatsRegionNameCol = NULL,
                                      commentString = "#",
                                      extractCoordinates = NULL) {
  binaryTraitModel <- match.arg(binaryTraitModel)
  if (is.null(genotypeList) & is.null(sumstatPathList)) {
    stop("Data load error. Please make sure at least one data set (sumstat_path_list or genotype_list) exists.")
  }

  # - if exist individual level data
  individualData <- NULL
  if (!is.null(genotypeList)) {
    if (length(phenotypeList) != length(covariateList)) {
      stop("Data load error. 'phenotype_list' and 'covariate_list' must have the same length.")
    }
    if (is.null(conditionsListIndividual)) {
      conditionsListIndividual <- paste0("condition", seq_along(phenotypeList))
      warning("Data load warning. 'conditions_list_individual' is not supplied; using default condition names. ",
              "Provide 'conditions_list_individual' to preserve cohort or cell-type labels.")
    }
    if (length(conditionsListIndividual) != length(phenotypeList)) {
      stop("Data load error. 'conditions_list_individual' must have the same length as 'phenotype_list'.")
    }
    #### FIXME: later if we have mulitple genotype list
    if (length(genotypeList) != 1 & is.null(matchGenoPheno)) {
      stop("Data load error. Please make sure 'match_geno_pheno' exists if you load data from multiple individual-level data.")
    } else if (length(genotypeList) == 1 & is.null(matchGenoPheno)) {
      matchGenoPheno <- rep(1, length(phenotypeList))
    }
    if (length(matchGenoPheno) != length(phenotypeList)) {
      stop("Data load error. 'match_geno_pheno' must have the same length as 'phenotype_list'.")
    }
    if (any(is.na(matchGenoPheno)) ||
        any(matchGenoPheno < 1 | matchGenoPheno > length(genotypeList))) {
      stop("Data load error. 'match_geno_pheno' must contain valid indices into 'genotype_list'.")
    }

    # - load individual data from multiple datasets
    nDataset <- unique(matchGenoPheno)
    for (iData in nDataset) {
      # extract genotype file name
      genotype <- genotypeList[iData]
      # extract phenotype and covariate file names
      pos <- which(matchGenoPheno == iData)
      phenotype <- phenotypeList[pos]
      covariate <- covariateList[pos]
      conditions <- conditionsListIndividual[pos]
      extractRegionNameI <- extractRegionName
      if (is.list(extractRegionName) && length(extractRegionName) == length(phenotypeList)) {
        extractRegionNameI <- extractRegionName[pos]
      }
      dat <- loadRegionalUnivariateData(
        genotype = genotype, phenotype = phenotype,
        covariate = covariate,
        region = region,
        associationWindow = associationWindow,
        conditions = conditions, xvarCutoff = xvarCutoff,
        mafCutoff = mafCutoff, macCutoff = macCutoff,
        imissCutoff = imissCutoff, keepIndel = keepIndel,
        keepSamples = keepSamples, keepVariants = keepVariants,
        extractRegionName = extractRegionNameI,
        phenotypeHeader = phenotypeHeader,
        regionNameCol = regionNameCol,
        scaleResiduals = scaleResiduals
      )
      if (is.null(individualData)) {
        individualData <- dat
      } else {
        individualData <- c(individualData, dat)
      }
    }
  }

  # - if exist summstat data
  sumstatData <- NULL
  if (!is.null(sumstatPathList)) {
    if (length(matchLdSumstat) == 0) {
      matchLdSumstat[[1]] <- conditionsListSumstat
    }
    if (length(matchLdSumstat) != length(ldMetaFilePathList)) {
      stop("Please make sure 'match_LD_sumstat' matched 'LD_meta_file_path_list' if you load data from multiple sumstats.")
    }
    # - load sumstat data from multiple datasets
    nLd <- length(matchLdSumstat)
    for (iLd in 1:nLd) {
      # extract LD meta file path name
      ldMetaFilePath <- ldMetaFilePathList[iLd]
      ldInfo <- loadLdMatrix(ldMetaFilePath,
        region = associationWindow,
        extractCoordinates = extractCoordinates,
        returnGenotype = "auto"
      )
      # extract sumstat information
      conditions <- matchLdSumstat[[iLd]]
      pos <- match(conditions, conditionsListSumstat)
      sumstats <- lapply(pos, function(ii) {
        sumstatPath <- sumstatPathList[ii]
        columnFilePath <- columnFilePathList[ii]
        # Load sumstat for this study (multiple LD references handled by outer loop)
        tmp <- loadRssData(
          sumstatPath = sumstatPath, columnFilePath = columnFilePath,
          nSample = nSamples[ii], nCase = nCases[ii], nControl = nControls[ii],
          region = associationWindow, extractRegionName = extractSumstatsRegionName,
          regionNameCol = sumstatsRegionNameCol, commentString = commentString,
          binaryTraitModel = binaryTraitModel
        )
        if (nrow(tmp$sumstats) == 0){ return(NULL) }
        if (!("variant_id" %in% colnames(tmp$sumstats))) {
          tmp$sumstats <- tmp$sumstats %>%
            mutate(variant_id = formatVariantId(chrom, pos, A2, A1))
        }
        return(tmp)
      })
      names(sumstats) <- conditions
      ifNoVariants <- sapply(sumstats, is.null)
      if (sum(ifNoVariants)!=0){
        posNoVariants <- which(ifNoVariants)
        sumstats <- sumstats[-posNoVariants]
      }
      sumstatData$sumstats <- c(sumstatData$sumstats, list(sumstats))
      sumstatData$ldInfo <- c(sumstatData$ldInfo, list(ldInfo))
    }
    names(sumstatData$sumstats) <- names(sumstatData$ldInfo) <- names(matchLdSumstat)
  }

  return(list(
    individualData = individualData,
    sumstatData = sumstatData
  ))
}

#' Convert loaded regional data to individual-level inputs
#'
#' @param regionData A list returned by \code{loadMultitaskRegionalData()}.
#' @return A list containing \code{X}, \code{Y}, \code{maf},
#'   \code{X_variance}, and source information.
#' @export
regionDataToIndInput <- function(regionData) {
  firstNonNull <- function(...) {
    values <- list(...)
    for (value in values) {
      if (!is.null(value)) return(value)
    }
    NULL
  }

  alignIndividualContexts <- function(X, Y) {
    cbindY <- function(Y, fallbackNames) {
      keep <- !vapply(Y, is.null, logical(1))
      if (!any(keep)) return(NULL)
      Y <- Y[keep]
      fallbackNames <- fallbackNames[keep]
      mats <- Map(function(y, nm) {
        if (is.null(dim(y))) y <- matrix(y, ncol = 1)
        if (is.null(colnames(y))) colnames(y) <- nm
        y
      }, Y, fallbackNames)
      do.call(cbind, mats)
    }

    if (!is.list(X) || is.matrix(X) || is.data.frame(X) ||
        !is.list(Y) || is.matrix(Y) || is.data.frame(Y) ||
        is.null(names(X)) || is.null(names(Y)) ||
        length(intersect(names(X), names(Y))) > 0) {
      return(list(X = X, Y = Y))
    }
    xNames <- names(X)
    yNames <- names(Y)
    grouped <- list()
    for (context in xNames) {
      matched <- yNames[yNames == context | startsWith(yNames, paste0(context, "_"))]
      if (length(matched) > 0) {
        yGroup <- cbindY(Y[matched], matched)
        if (!is.null(yGroup)) grouped[[context]] <- yGroup
      }
    }
    if (length(grouped) == 0 && length(X) == 1 && length(Y) > 0) {
      yGroup <- cbindY(Y, yNames)
      if (!is.null(yGroup)) grouped[[xNames[[1]]]] <- yGroup
    }
    if (length(grouped) == 0) {
      return(list(X = X, Y = Y))
    }
    list(X = X[names(grouped)], Y = grouped)
  }

  individualData <- regionData$individualData
  if (is.null(individualData)) {
    return(list(X = NULL, Y = NULL, maf = NULL, xVariance = NULL,
                sourceInfo = list(hasIndividual = FALSE, contexts = character())))
  }

  if (is(individualData, "RegionalData")) {
    contexts <- names(individualData@phenotypes)
    nCond <- length(contexts)
    X <- stats::setNames(
      lapply(seq_len(nCond), function(i) getResidualX(individualData, i)),
      contexts
    )
    Y <- stats::setNames(
      lapply(seq_len(nCond), function(i) getResidualY(individualData, i)),
      contexts
    )
    aligned <- alignIndividualContexts(X, Y)
    X <- aligned$X
    Y <- aligned$Y
    maf <- individualData@maf
    XVariance <- stats::setNames(
      lapply(seq_len(nCond), function(i) getXVariance(individualData, i)),
      contexts
    )
    return(list(
      X = X,
      Y = Y,
      maf = maf,
      xVariance = XVariance,
      sourceInfo = list(hasIndividual = !is.null(X) && !is.null(Y),
                        contexts = contexts)
    ))
  }

  # Post-QC shape: list(X = list_of_matrices, Y = list_of_matrices, ...)
  if (is.list(individualData) &&
      (!is.null(individualData$X) || !is.null(individualData$Y))) {
    X <- individualData$X
    Y <- individualData$Y
    aligned <- alignIndividualContexts(X, Y)
    X <- aligned$X
    Y <- aligned$Y
    maf <- individualData$maf
    XVariance <- individualData$X_variance
    contexts <- if (!is.null(X) && is.list(X) && !is.matrix(X)) names(X) else character()
    return(list(
      X = X,
      Y = Y,
      maf = maf,
      xVariance = XVariance,
      sourceInfo = list(hasIndividual = !is.null(X) && !is.null(Y),
                        contexts = contexts)
    ))
  }

  stop("region_data$individualData must be a RegionalData object or a post-QC list with X/Y entries")
}

#' Convert loaded regional data to RSS inputs
#'
#' @param regionData A list returned by \code{loadMultitaskRegionalData()}.
#' @return A list containing named RSS inputs, matched LD data, and source
#'   information.
#' @export
regionDataToRssInput <- function(regionData) {
  rssInputFromQcedSumstat <- function(sumstatData) {
    rssInput <- sumstatData$sumstats
    ldDataIn <- sumstatData$ldData
    ldMatch <- sumstatData$ldMatch
    studies <- names(rssInput)
    ldData <- list()
    ldGroup <- character()
    for (i in seq_along(studies)) {
      study <- studies[[i]]
      ldName <- if (!is.null(ldMatch) && length(ldMatch) >= i) ldMatch[[i]] else study
      if (is.null(ldName) || is.na(ldName) || !ldName %in% names(ldDataIn)) {
        ldName <- names(ldDataIn)[min(i, length(ldDataIn))]
      }
      ld <- ldDataIn[[ldName]]
      if (!is.null(ld) && !is(ld, "LdData")) {
        stop("region_data$sumstatData$ldData entries must be LdData objects.")
      }
      ldData[[study]] <- ld
      ldGroup[[study]] <- ldName
    }
    list(
      rssInput = rssInput,
      ldData = ldData,
      sourceInfo = list(hasSumstat = length(rssInput) > 0,
                         studies = names(rssInput),
                         ldGroup = ldGroup)
    )
  }

  sumstatData <- regionData$sumstatData
  if (is.null(sumstatData) || is.null(sumstatData$sumstats)) {
    return(list(rssInput = list(), ldData = list(),
                sourceInfo = list(hasSumstat = FALSE, studies = character(),
                                   ldGroup = character())))
  }
  if (!is.null(sumstatData$ldData)) {
    return(rssInputFromQcedSumstat(sumstatData))
  }

  rssInput <- list()
  ldData <- list()
  ldGroup <- character()

  for (i in seq_along(sumstatData$sumstats)) {
    studies <- sumstatData$sumstats[[i]]
    ldIndex <- min(i, length(sumstatData$ldInfo))
    groupName <- names(sumstatData$ldInfo)[ldIndex]
    if (is.null(groupName) || is.na(groupName) || groupName == "") {
      groupName <- paste0("LD", ldIndex)
    }
    ldEntry <- sumstatData$ldInfo[[ldIndex]]
    if (!is.null(ldEntry) && !is(ldEntry, "LdData")) {
      stop("region_data$sumstatData$ldInfo entries must be LdData objects.")
    }
    for (study in names(studies)) {
      outputName <- study
      if (outputName %in% names(rssInput)) {
        outputName <- make.unique(c(names(rssInput), outputName))[length(rssInput) + 1]
      }
      rssInput[[outputName]] <- studies[[study]]
      ldData[[outputName]] <- ldEntry
      ldGroup[[outputName]] <- groupName
    }
  }

  list(
    rssInput = rssInput,
    ldData = ldData,
    sourceInfo = list(hasSumstat = length(rssInput) > 0,
                       studies = names(rssInput),
                       ldGroup = ldGroup)
  )
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
