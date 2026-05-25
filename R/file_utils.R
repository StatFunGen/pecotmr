# read PLINK files

#' @importFrom vroom vroom
#' @importFrom tools file_path_sans_ext
read_bim <- function(bed) {
  bimf <- paste0(file_path_sans_ext(bed), ".bim")
  bim <- vroom(bimf, col_names = FALSE)
  colnames(bim) <- c("chrom", "id", "gpos", "pos", "a1", "a0")
  return(bim)
}

#' @importFrom vroom vroom
#' @importFrom tools file_path_sans_ext
read_fam <- function(bed) {
  famf <- paste0(file_path_sans_ext(bed), ".fam")
  return(vroom(famf, col_names = FALSE))
}

# open bed/bim/fam: A PLINK 1 .bed is a valid .pgen
open_bed <- function(bed) {
  if (!requireNamespace("pgenlibr", quietly = TRUE)) {
    stop("To use this function, please install pgenlibr: https://cran.r-project.org/web/packages/pgenlibr/index.html")
  }
  raw_s_ct <- nrow(read_fam(bed))
  return(pgenlibr::NewPgen(bed, raw_sample_ct = raw_s_ct))
}

#' Read a PLINK2 allele frequency file (.afreq or .afreq.zst)
#'
#' @param prefix File prefix (without .afreq extension).
#' @return A data.frame with columns: chrom, id, A2 (REF), A1 (ALT), alt_freq, obs_ct.
#'   alt_freq is the frequency of the A1 (ALT/effect) allele.
#' @importFrom vroom vroom
#' @importFrom dplyr rename select
#' @export
read_afreq <- function(prefix) {
  afreq_zst <- paste0(prefix, ".afreq.zst")
  afreq_plain <- paste0(prefix, ".afreq")
  if (file.exists(afreq_zst)) {
    if (Sys.which("zstd") == "") stop("zstd CLI is required to read .afreq.zst files")
    af <- as.data.frame(vroom(pipe(paste0("zstd -dcq ", shQuote(afreq_zst))),
                              delim = "\t", show_col_types = FALSE))
  } else if (file.exists(afreq_plain)) {
    af <- as.data.frame(vroom(afreq_plain, delim = "\t", show_col_types = FALSE))
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
#'     (read via \code{read_afreq}, which also returns allele frequencies).
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
read_stochastic_meta <- function(path, format = NULL) {
  if (!file.exists(path)) return(NULL)

  if (is.null(format)) {
    format <- if (grepl("\\.afreq(\\.zst)?$", path)) "afreq" else "generic"
  }
  format <- match.arg(format, c("afreq", "generic"))

  if (format == "afreq") {
    # read_afreq expects a prefix, not a full path - strip the .afreq[.zst] suffix
    prefix <- sub("\\.afreq(\\.zst)?$", "", path)
    af <- read_afreq(prefix)
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
#' @param genotype_path Path to the genotype data (prefix or file path).
#' @return Path to the first sidecar file found, or \code{NULL}.
#' @noRd
find_stochastic_meta <- function(genotype_path) {
  # Strip known genotype extensions to get the stem
  stem <- sub("\\.(vcf|vcf\\.gz|bcf|gds|bed|bim|fam|pgen|pvar|psam)$", "",
              genotype_path)
  candidates <- c(
    paste0(stem, ".afreq"),
    paste0(stem, ".afreq.zst"),
    paste0(stem, ".stochastic_meta.tsv")
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) > 0) found[1] else NULL
}

#' Load PLINK2 genotype data via pgenlibr
#'
#' Loads genotype data from PLINK2 format files using pgenlibr directly
#' (no plink2 CLI required). Supports uncompressed .pgen with .pvar or
#' .pvar.zst (the standard plink2 layout produced by \code{--make-pgen vzs}).
#' The .psam file must be uncompressed (plink2 never compresses it).
#'
#' Dosage convention: X contains ALT/A1 (effect allele) dosage counts (0, 1, 2),
#' consistent with the PLINK1 path in \code{load_genotype_region()} which returns
#' A1 dosage via \code{2 - as(geno_bed, "numeric")}.
#'
#' @param prefix File prefix (without extension). Expected layout:
#'   prefix.pgen (uncompressed), prefix.pvar or prefix.pvar.zst,
#'   prefix.psam (uncompressed).
#' @param region Target region in format "chr:start-end" (e.g., "chr1:1000-2000").
#'   If NULL, loads all variants.
#' @param keep_indel Whether to keep indel variants. Default TRUE.
#' @param keep_variants_path Path to a file listing variants to keep. Default NULL.
#' @return A list with:
#'   \item{X}{Numeric ALT/A1 dosage matrix. Rows are samples named by IID
#'     (individual ID from .psam). Columns are variants named by the ID column
#'     from .pvar. Values 0/1/2 count copies of the A1 (ALT/effect) allele.}
#'   \item{variant_info}{Data.frame with columns: chrom, id, pos, A2 (REF allele),
#'     A1 (ALT/effect allele). If a .afreq[.zst] file exists at the same prefix,
#'     also includes alt_freq (A1 frequency) and obs_ct.}
#'
#' @importFrom vroom vroom
#' @noRd
load_plink2_data <- function(prefix, region = NULL, keep_indel = TRUE, keep_variants_path = NULL) {
  if (!requireNamespace("pgenlibr", quietly = TRUE)) {
    stop("pgenlibr is required. Install from https://cran.r-project.org/web/packages/pgenlibr/index.html")
  }

  paths <- resolve_plink2_paths(prefix)

  # --- Read variant info from .pvar as text (pgenlibr::NewPvar is unreliable) ---
  all_variant_info <- read_pvar(paths$pvar)

  variant_idx <- seq_len(nrow(all_variant_info))
  if (!is.null(region)) {
    parsed <- parse_region(region)
    in_region <- strip_chr_prefix(all_variant_info$chrom) == parsed$chrom &
                 all_variant_info$pos >= parsed$start &
                 all_variant_info$pos <= parsed$end
    variant_idx <- which(in_region)
    if (length(variant_idx) == 0) {
      stop(NoSNPsError(paste("No variants found in region", region)))
    }
  }
  variant_info <- all_variant_info[variant_idx, , drop = FALSE]

  # --- Read samples from .psam ---
  psam <- as.data.frame(vroom(paths$psam, delim = "\t", show_col_types = FALSE))
  names(psam) <- sub("^#", "", names(psam))

  # --- Read genotype dosage via pgenlibr ---
  pgen <- pgenlibr::NewPgen(paths$pgen)
  on.exit(pgenlibr::ClosePgen(pgen), add = TRUE)
  X <- pgenlibr::ReadList(pgen, variant_subset = variant_idx, meanimpute = FALSE)
  rownames(X) <- psam$IID
  colnames(X) <- variant_info$id

  # --- Attach allele frequency from .afreq sidecar ---
  afreq <- read_afreq(prefix)
  if (!is.null(afreq)) {
    afreq_cols <- intersect(c("id", "alt_freq", "obs_ct"), colnames(afreq))
    variant_info <- merge(variant_info, afreq[, afreq_cols, drop = FALSE],
                          by = "id", all.x = TRUE, sort = FALSE)
  }

  # --- Post-filters: indels and variant whitelist ---
  if (!keep_indel) {
    snp_mask <- is_snp_alleles(variant_info$A1, variant_info$A2)
    X <- X[, snp_mask, drop = FALSE]
    variant_info <- variant_info[snp_mask, , drop = FALSE]
  }
  if (!is.null(keep_variants_path)) {
    keep_idx <- match_variants_to_keep(variant_info, keep_variants_path)
    X <- X[, keep_idx, drop = FALSE]
    variant_info <- variant_info[keep_idx, , drop = FALSE]
  }

  list(X = X, variant_info = variant_info)
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
#' @param u_min Numeric vector of per-variant minimum values before scaling.
#' @param u_max Numeric vector of per-variant maximum values before scaling.
#' @return Matrix of original U values with same dimensions.
#' @export
invert_minmax_scaling <- function(X, u_min, u_max) {
  if (length(u_min) != ncol(X) || length(u_max) != ncol(X)) {
    stop("Length of u_min/u_max (", length(u_min), ") must equal ncol(X) (", ncol(X), ")")
  }
  denom <- u_max - u_min
  denom[denom == 0] <- 1  # monomorphic: scaling was identity
  # Invert: U_original = U_scaled * (u_max - u_min) / 2 + u_min
  sweep(sweep(X, 2, denom / 2, "*"), 2, u_min, "+")
}

# ---------- Internal helpers for load_plink2_data ----------

#' Resolve and validate PLINK2 file paths for a given prefix.
#' @return Named list with pgen, pvar, psam paths.
#' @noRd
resolve_plink2_paths <- function(prefix) {
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
#' @param pvar_path Path to .pvar or .pvar.zst file.
#' @return data.frame with columns: chrom, id, pos, A2 (REF), A1 (ALT).
#' @noRd
read_pvar <- function(pvar_path) {
  if (!requireNamespace("pgenlibr", quietly = TRUE)) {
    stop("pgenlibr is required. Install from https://cran.r-project.org/web/packages/pgenlibr/index.html")
  }
  pvar <- pgenlibr::NewPvar(pvar_path)
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
#' delegates to \code{read_pvar()}.
#'
#' @param snp_file_path Path to .bim, .pvar, or .pvar.zst file.
#' @return data.frame with at minimum columns: chrom, id, pos, A2, A1.
#'   Extended .bim files (9 columns) also include: variance, allele_freq, n_nomiss.
#' @importFrom utils read.table
#' @noRd
read_variant_metadata <- function(snp_file_path) {
  is_pvar <- grepl("\\.(pvar|pvar\\.zst)$", snp_file_path)
  if (!is_pvar) {
    first_line <- readLines(snp_file_path, n = 1)
    is_pvar <- grepl("^#CHROM", first_line)
  }

  if (is_pvar) {
    read_pvar(snp_file_path)
  } else {
    df <- read.table(snp_file_path, stringsAsFactors = FALSE)
    n <- ncol(df)
    if (n == 6) {
      names(df) <- c("chrom", "id", "gpos", "pos", "A1", "A2")
    } else if (n == 9) {
      names(df) <- c("chrom", "id", "gpos", "pos", "A1", "A2", "variance", "allele_freq", "n_nomiss")
    } else {
      stop("Unexpected number of columns (", n, ") in variant file: ", snp_file_path)
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
get_ref_variant_info <- function(source, region = NULL) {
  resolved <- resolve_ld_source(source)

  # For genotype sources via metadata, resolve per-chromosome path
  if (resolved$type %in% c("plink2", "plink1", "vcf", "gds") && !is.null(resolved$meta_path) && !is.null(region)) {
    data_path <- resolve_genotype_path_for_region(resolved$meta_path, region)
  } else {
    data_path <- resolved$data_path
  }

  if (resolved$type == "plink2") {
    paths <- resolve_plink2_paths(data_path)
    info <- read_pvar(paths$pvar)
    afreq <- read_afreq(data_path)
    if (!is.null(afreq)) {
      info$allele_freq <- afreq$alt_freq[match(info$id, afreq$id)]
    }
  } else if (resolved$type == "plink1") {
    bim <- read_bim(paste0(data_path, ".bed"))
    info <- data.frame(
      chrom = bim$chrom, id = bim$id, pos = bim$pos,
      A2 = bim$a0, A1 = bim$a1,
      stringsAsFactors = FALSE
    )
  } else if (resolved$type %in% c("vcf", "gds")) {
    # VCF/GDS: load via the genotype loader and extract variant_info
    result <- load_genotype_region(data_path, region = region,
                                   return_variant_info = TRUE)
    info <- result$variant_info
    # Compute allele frequency from the genotype matrix
    info$allele_freq <- colMeans(result$X, na.rm = TRUE) / 2
    return(info)  # Already region-filtered by the loader
  } else {
    # Pre-computed LD: read bim/pvar files via metadata
    bim_paths <- get_regional_ld_meta(resolved$meta_path, region)$intersections$bim_file_paths
    info <- do.call(rbind, lapply(bim_paths, function(path) {
      df <- read_variant_metadata(path)
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
    info$id <- normalize_variant_id(info$id)
    return(info)  # Already region-filtered by get_regional_ld_meta
  }

  # Region filter for plink2/plink1
  if (!is.null(region)) {
    parsed <- parse_region(region)
    info_chrom <- strip_chr_prefix(info$chrom)
    # Handle multi-row region data.frame (one row per chrom)
    if (is.data.frame(parsed) && nrow(parsed) > 1) {
      in_region <- rep(FALSE, nrow(info))
      for (r in seq_len(nrow(parsed))) {
        in_region <- in_region | (info_chrom == as.character(parsed$chrom[r]) &
                                  info$pos >= parsed$start[r] & info$pos <= parsed$end[r])
      }
    } else {
      in_region <- info_chrom == as.character(parsed$chrom) &
                   info$pos >= parsed$start & info$pos <= parsed$end
    }
    info <- info[in_region, , drop = FALSE]
  }
  info
}

#' Match variant_info against a whitelist file, returning logical index.
#' Uses parse_variant_id() from misc.R to handle all variant ID formats.
#' @importFrom vroom vroom
#' @importFrom readr read_lines
#' @noRd
match_variants_to_keep <- function(variant_info, keep_variants_path) {
  keep_raw <- tryCatch(
    as.data.frame(vroom(keep_variants_path, show_col_types = FALSE)),
    error = function(e) NULL
  )
  if (!is.null(keep_raw) && "chrom" %in% names(keep_raw) && "pos" %in% names(keep_raw)) {
    keep_variants <- parse_variant_id(keep_raw)
  } else {
    # Fall back to reading as single-column variant IDs
    ids <- read_lines(keep_variants_path)
    keep_variants <- parse_variant_id(ids)
  }
  vi_chrom <- as.integer(strip_chr_prefix(variant_info$chrom))
  has_alleles <- "A1" %in% names(keep_variants) && "A2" %in% names(keep_variants) &&
    !any(is.na(keep_variants$A1)) && !any(is.na(keep_variants$A2))
  if (has_alleles) {
    paste0(vi_chrom, ":", variant_info$pos, ":", variant_info$A2, ":", variant_info$A1) %in%
      paste0(keep_variants$chrom, ":", keep_variants$pos, ":", keep_variants$A2, ":", keep_variants$A1)
  } else {
    paste0(vi_chrom, ":", variant_info$pos) %in%
      paste0(keep_variants$chrom, ":", keep_variants$pos)
  }
}

#' @importFrom vroom vroom
#' @importFrom dplyr as_tibble mutate filter
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom stringr str_detect

# Internal helper: read a region from a tabix-indexed file via Rsamtools
read_tabix_region <- function(file, region, use_col_names) {
  tbx <- Rsamtools::TabixFile(file)
  parsed <- parse_region(region)
  # Match chromosome naming convention in the tabix index
  chrom <- as.character(parsed$chrom)
  tbx_seqnames <- Rsamtools::seqnamesTabix(tbx)
  if (any(grepl("^chr", tbx_seqnames))) {
    chrom <- paste0("chr", chrom)
  }
  gr <- GenomicRanges::GRanges(
    seqnames = chrom,
    ranges = IRanges::IRanges(start = parsed$start, end = parsed$end)
  )
  lines <- Rsamtools::scanTabix(tbx, param = gr)[[1]]
  if (length(lines) == 0) return(NULL)

  # Get header for column names
  col_names_vec <- NULL
  if (use_col_names) {
    hdr <- Rsamtools::headerTabix(tbx)$header
    if (length(hdr) > 0) {
      last_hdr <- hdr[length(hdr)]
      col_names_vec <- strsplit(sub("^#", "", last_hdr), "\t")[[1]]
    }
  }

  # Parse tab-delimited lines
  txt <- paste(lines, collapse = "\n")
  if (!is.null(col_names_vec)) {
    as.data.frame(vroom::vroom(I(txt), delim = "\t", col_names = col_names_vec,
                               show_col_types = FALSE))
  } else {
    as.data.frame(vroom::vroom(I(txt), delim = "\t", col_names = use_col_names,
                               show_col_types = FALSE))
  }
}

tabix_region <- function(file, region, tabix_header = "auto", target = "", target_column_index = "") {
  if (!file.exists(file)) {
    stop("Input file does not exist: ", file)
  }

  use_col_names <- if (identical(tabix_header, FALSE)) FALSE else TRUE

  cmd_output <- tryCatch(
    read_tabix_region(file, region, use_col_names),
    error = function(e) NULL
  )

  if (!is.null(cmd_output) && target != "" && target_column_index != "") {
    cmd_output <- cmd_output %>%
      filter(str_detect(.[[target_column_index]], target))
  } else if (!is.null(cmd_output) && target != "") {
    cmd_output <- cmd_output %>%
      mutate(text = apply(., 1, function(row) paste(row, collapse = "_"))) %>%
      filter(str_detect(text, target)) %>%
      select(-text)
  }

  if (is.null(cmd_output) || nrow(cmd_output) == 0) {
    return(tibble())
  }

  cmd_output %>%
    as_tibble() %>%
    mutate(
      !!names(.)[1] := as.character(.[[1]]),
      !!names(.)[2] := as.numeric(.[[2]])
    )
}


NoSNPsError <- function(message) {
  structure(list(message = message), class = c("NoSNPsError", "error", "condition"))
}

#' Load PLINK1 genotype data via snpStats
#'
#' Loads genotype data from PLINK1 format files (.bed/.bim/.fam) using snpStats.
#' Returns the same structure as \code{load_plink2_data()} for consistency.
#'
#' Dosage convention: X contains A1 (effect allele) dosage counts (0, 1, 2),
#' computed as \code{2 - as(geno_bed, "numeric")} from snpStats encoding.
#'
#' @param prefix File prefix (without extension). Expected: prefix.bed, prefix.bim, prefix.fam.
#' @param region Target region in format "chr:start-end" (e.g., "chr1:1000-2000").
#'   If NULL, loads all variants.
#' @param keep_indel Whether to keep indel variants. Default TRUE.
#' @param keep_variants_path Path to a file listing variants to keep. Default NULL.
#' @return A list with:
#'   \item{X}{Numeric A1 dosage matrix. Rows are samples, columns are variants.}
#'   \item{variant_info}{Data.frame with columns: chrom, id, pos, A2 (allele.2),
#'     A1 (allele.1/effect allele).}
#'
#' @importFrom vroom vroom
#' @importFrom readr col_character col_guess col_integer
#' @noRd
load_plink1_data <- function(prefix, region = NULL, keep_indel = TRUE, keep_variants_path = NULL) {
  bed_file <- paste0(prefix, ".bed")
  bim_file <- paste0(prefix, ".bim")
  fam_file <- paste0(prefix, ".fam")
  if (!all(file.exists(bed_file, bim_file, fam_file))) {
    stop("PLINK1 fileset (.bed/.bim/.fam) not found at prefix: ", prefix)
  }
  if (!requireNamespace("snpStats", quietly = TRUE)) {
    stop("snpStats is required. Install from https://bioconductor.org/packages/release/bioc/html/snpStats.html")
  }

  # --- Region filter via bim ---
  if (!is.null(region)) {
    parsed <- parse_region(region)
    bim_data <- read_bim(bed_file)
    bim_data$chrom <- strip_chr_prefix(bim_data$chrom)
    snp_ids <- bim_data$id[bim_data$chrom == parsed$chrom &
                            bim_data$pos >= parsed$start &
                            bim_data$pos <= parsed$end]
    if (length(snp_ids) == 0) {
      stop(NoSNPsError(paste("No SNPs found in the specified region", region)))
    }
  } else {
    snp_ids <- NULL
  }

  geno <- snpStats::read.plink(prefix, select.snps = snp_ids)
  geno_map <- geno$map

  # --- Build variant_info from snpStats $map ---
  # snpStats: allele.1 = effect allele (A1), allele.2 = reference allele (A2)
  variant_info <- data.frame(
    chrom = as.character(geno_map$chromosome),
    id    = rownames(geno_map),
    pos   = geno_map$position,
    A2    = as.character(geno_map$allele.2),
    A1    = as.character(geno_map$allele.1),
    stringsAsFactors = FALSE
  )

  # --- Dosage matrix: 2 - snpStats encoding gives A1 dosage ---
  X <- 2 - as(geno$genotypes, "numeric")
  rownames(X) <- rownames(geno$genotypes)
  colnames(X) <- variant_info$id

  # --- Post-filters: indels and variant whitelist ---
  if (!keep_indel) {
    snp_mask <- is_snp_alleles(variant_info$A1, variant_info$A2)
    X <- X[, snp_mask, drop = FALSE]
    variant_info <- variant_info[snp_mask, , drop = FALSE]
  }
  if (!is.null(keep_variants_path)) {
    keep_idx <- match_variants_to_keep(variant_info, keep_variants_path)
    X <- X[, keep_idx, drop = FALSE]
    variant_info <- variant_info[keep_idx, , drop = FALSE]
  }

  return(list(X = X, variant_info = variant_info))
}

#' Load genotype data from a VCF file via VariantAnnotation
#'
#' Reads biallelic SNP genotypes from a VCF (or VCF.gz/BCF) file using
#' VariantAnnotation::readVcf(). Extracts GT field and converts to 0/1/2
#' ALT dosage. Returns the same structure as load_plink2_data().
#'
#' @param path Path to VCF file (.vcf, .vcf.gz, or .bcf).
#' @param region Target region in format "chr:start-end". If NULL, loads all.
#' @param keep_indel Whether to keep indel variants. Default TRUE.
#' @param keep_variants_path Path to a file listing variants to keep.
#' @return A list with X (dosage matrix) and variant_info (data.frame).
#' @noRd
load_vcf_data <- function(path, region = NULL, keep_indel = TRUE,
                          keep_variants_path = NULL) {
  if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
    stop("VariantAnnotation is required for VCF loading. ",
         "Install from Bioconductor: BiocManager::install('VariantAnnotation')")
  }

  # Build scan parameters
  param <- if (!is.null(region)) {
    parsed <- parse_region(region)
    chrom_name <- as.character(parsed$chrom)
    # Match chromosome naming convention in the VCF header
    vcf_seqnames <- Rsamtools::seqnamesTabix(Rsamtools::TabixFile(path))
    if (any(grepl("^chr", vcf_seqnames))) {
      chrom_name <- paste0("chr", chrom_name)
    }
    gr <- GenomicRanges::GRanges(
      seqnames = chrom_name,
      ranges = IRanges::IRanges(start = parsed$start, end = parsed$end)
    )
    VariantAnnotation::ScanVcfParam(which = gr, geno = "GT")
  } else {
    VariantAnnotation::ScanVcfParam(geno = "GT")
  }

  # Read VCF
  vcf <- VariantAnnotation::readVcf(path, param = param)
  if (length(vcf) == 0) {
    stop(NoSNPsError(paste("No variants found in VCF", path,
                           if (!is.null(region)) paste("for region", region))))
  }

  # Extract GT matrix and convert to ALT dosage (0/1/2)
  gt <- VariantAnnotation::geno(vcf)$GT
  # GT is a character matrix: "0/0", "0/1", "1/1", "0|0", "0|1", "1|0", "1|1"
  dosage <- matrix(NA_real_, nrow = ncol(gt), ncol = nrow(gt))
  rownames(dosage) <- colnames(gt)
  for (j in seq_len(nrow(gt))) {
    g <- gt[j, ]
    alleles <- strsplit(g, "[/|]")
    dosage[, j] <- vapply(alleles, function(a) {
      a <- as.integer(a)
      if (any(is.na(a))) NA_real_ else sum(a)
    }, numeric(1))
  }

  # Build variant_info from rowRanges
  rr <- SummarizedExperiment::rowRanges(vcf)
  variant_info <- data.frame(
    chrom = as.character(GenomicRanges::seqnames(rr)),
    id    = names(rr),
    pos   = as.integer(GenomicRanges::start(rr)),
    A2    = as.character(VariantAnnotation::ref(vcf)),
    A1    = vapply(VariantAnnotation::alt(vcf),
                   function(x) as.character(x)[1], character(1)),
    stringsAsFactors = FALSE
  )
  colnames(dosage) <- variant_info$id

  # Post-filters
  if (!keep_indel) {
    snp_mask <- is_snp_alleles(variant_info$A1, variant_info$A2)
    dosage <- dosage[, snp_mask, drop = FALSE]
    variant_info <- variant_info[snp_mask, , drop = FALSE]
  }
  if (!is.null(keep_variants_path)) {
    keep_idx <- match_variants_to_keep(variant_info, keep_variants_path)
    dosage <- dosage[, keep_idx, drop = FALSE]
    variant_info <- variant_info[keep_idx, , drop = FALSE]
  }

  list(X = dosage, variant_info = variant_info)
}

#' Load genotype data from a GDS file via SNPRelate
#'
#' Reads genotype data from a CoreArray GDS file using SNPRelate.
#' Returns the same structure as load_plink2_data().
#'
#' @param path Path to GDS file (.gds).
#' @param region Target region in format "chr:start-end". If NULL, loads all.
#' @param keep_indel Whether to keep indel variants. Default TRUE.
#' @param keep_variants_path Path to a file listing variants to keep.
#' @return A list with X (dosage matrix) and variant_info (data.frame).
#' @noRd
load_gds_data <- function(path, region = NULL, keep_indel = TRUE,
                          keep_variants_path = NULL) {
  if (!requireNamespace("SNPRelate", quietly = TRUE)) {
    stop("SNPRelate is required for GDS loading. ",
         "Install from Bioconductor: BiocManager::install('SNPRelate')")
  }

  gds <- SNPRelate::snpgdsOpen(path, readonly = TRUE, allow.duplicate = TRUE)
  on.exit(SNPRelate::snpgdsClose(gds), add = TRUE)

  # Read variant metadata
  snp_chrom <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.chromosome"))
  snp_pos <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.position"))
  snp_id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.id"))
  snp_allele <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.allele"))

  # Parse alleles: "REF/ALT" format
  allele_split <- strsplit(snp_allele, "/")
  a2 <- vapply(allele_split, `[`, character(1), 1L)  # REF
  a1 <- vapply(allele_split, `[`, character(1), 2L)  # ALT

  # Region filter
  snp_subset <- NULL
  if (!is.null(region)) {
    parsed <- parse_region(region)
    in_region <- strip_chr_prefix(as.character(snp_chrom)) == parsed$chrom &
                 snp_pos >= parsed$start & snp_pos <= parsed$end
    if (!any(in_region)) {
      stop(NoSNPsError(paste("No variants found in region", region)))
    }
    snp_subset <- snp_id[in_region]
  }

  # Read genotype matrix (samples x variants, ALT dosage 0/1/2)
  geno <- SNPRelate::snpgdsGetGeno(gds, snp.id = snp_subset,
                                    with.id = TRUE, verbose = FALSE)
  X <- geno$genotype
  rownames(X) <- geno$sample.id

  # Build variant_info for the subset
  if (!is.null(snp_subset)) {
    idx <- match(geno$snp.id, snp_id)
  } else {
    idx <- seq_along(snp_id)
  }
  variant_info <- data.frame(
    chrom = as.character(snp_chrom[idx]),
    id    = as.character(snp_id[idx]),
    pos   = as.integer(snp_pos[idx]),
    A2    = a2[idx],
    A1    = a1[idx],
    stringsAsFactors = FALSE
  )
  colnames(X) <- variant_info$id

  # Post-filters
  if (!keep_indel) {
    snp_mask <- is_snp_alleles(variant_info$A1, variant_info$A2)
    X <- X[, snp_mask, drop = FALSE]
    variant_info <- variant_info[snp_mask, , drop = FALSE]
  }
  if (!is.null(keep_variants_path)) {
    keep_idx <- match_variants_to_keep(variant_info, keep_variants_path)
    X <- X[, keep_idx, drop = FALSE]
    variant_info <- variant_info[keep_idx, , drop = FALSE]
  }

  list(X = X, variant_info = variant_info)
}

#' Load genotype data for a specific region
#'
#' Auto-detects PLINK2 (.pgen/.pvar[.zst]/.psam), PLINK1 (.bed/.bim/.fam),
#' VCF (.vcf/.vcf.gz/.bcf), or GDS (.gds) format and loads genotype data
#' accordingly. If a stochastic genotype sidecar file (.afreq or
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
#'   sidecar file. If NULL (default), auto-detected via \code{find_stochastic_meta}.
#' @param stochastic_meta_format Optional format override for the sidecar file:
#'   \code{"afreq"} or \code{"generic"}. If NULL (default), auto-detected from
#'   file extension.
#' @return If return_variant_info is FALSE, a numeric dosage matrix (rows=samples,
#'   cols=variants). If TRUE, a list with elements X and variant_info.
#'
#' @export
load_genotype_region <- function(genotype, region = NULL, keep_indel = TRUE,
                                 keep_variants_path = NULL,
                                 return_variant_info = FALSE,
                                 stochastic_meta_path = NULL,
                                 stochastic_meta_format = NULL) {
  result <- NULL
  # VCF and GDS: detect by file extension on the path itself
  if (grepl("\\.(vcf|vcf\\.gz|bcf)$", genotype)) {
    result <- load_vcf_data(genotype, region = region, keep_indel = keep_indel,
                            keep_variants_path = keep_variants_path)
  } else if (grepl("\\.gds$", genotype)) {
    result <- load_gds_data(genotype, region = region, keep_indel = keep_indel,
                            keep_variants_path = keep_variants_path)
  } else if (has_plink2_files(genotype)) {
    # PLINK prefix detection (no extension on the path)
    result <- load_plink2_data(genotype, region = region, keep_indel = keep_indel,
                               keep_variants_path = keep_variants_path)
  } else if (has_plink1_files(genotype)) {
    result <- load_plink1_data(genotype, region = region, keep_indel = keep_indel,
                               keep_variants_path = keep_variants_path)
  } else {
    stop("Genotype files not found at: ", genotype,
         "\n  Expected: .vcf/.vcf.gz/.bcf, .gds, or PLINK prefix (.pgen/.pvar[.zst]/.psam or .bed/.bim/.fam)")
  }

  # --- Detect and invert stochastic genotype scaling ---
  meta_path <- stochastic_meta_path %||% find_stochastic_meta(genotype)
  if (!is.null(meta_path)) {
    smeta <- read_stochastic_meta(meta_path, format = stochastic_meta_format)
    if (!is.null(smeta)) {
      idx <- match(colnames(result$X), smeta$id)
      matched <- !is.na(idx)
      if (any(matched)) {
        result$X[, matched] <- invert_minmax_scaling(
          result$X[, matched, drop = FALSE],
          smeta$u_min[idx[matched]],
          smeta$u_max[idx[matched]]
        )
        result$variant_info$u_min <- smeta$u_min[idx]
        result$variant_info$u_max <- smeta$u_max[idx]
        message("Stochastic genotype detected: restored original scale via ", basename(meta_path))
      }
    }
  } else {
    is_stochastic <- !all(result$X == round(result$X), na.rm = TRUE)
    if (is_stochastic) {
      warning("Non-integer genotype values detected but no stochastic metadata sidecar found. ",
              "Place a .afreq or .stochastic_meta.tsv file with u_min/u_max columns ",
              "alongside the genotype files to restore the original scale.")
    }
  }

  if (return_variant_info) result else result$X
}

#' @importFrom purrr map
#' @importFrom readr read_delim cols
#' @importFrom dplyr select mutate across everything
#' @importFrom magrittr %>%
#' @noRd
read_single_covariate <- function(path) {
  raw_df <- read_delim(path, "\t", col_types = cols(.default = "c")) %>% select(-1)
  df <- raw_df
  non_numeric <- character()
  for (nm in names(df)) {
    values <- trimws(as.character(df[[nm]]))
    converted <- suppressWarnings(as.numeric(values))
    bad <- !is.na(values) & values != "" & is.na(converted)
    if (any(bad)) {
      non_numeric <- c(non_numeric, nm)
    } else {
      df[[nm]] <- converted
    }
  }
  if (length(non_numeric) > 0) {
    stop("Non-numeric columns found in covariate file ", path, ": ",
         paste(non_numeric, collapse = ", "),
         ". All columns except the first (sample ID) must be numeric.")
  }
  df %>% mutate(across(everything(), as.numeric)) %>% t()
}

#' @noRd
load_covariate_data <- function(covariate_path) {
  # Validate all covariate files exist
  missing <- covariate_path[!file.exists(covariate_path)]
  if (length(missing) > 0) {
    stop("Covariate file(s) not found: ", paste(missing, collapse = ", "))
  }
  return(map(covariate_path, read_single_covariate))
}

NoPhenotypeError <- function(message) {
  structure(list(message = message), class = c("NoPhenotypeError", "error", "condition"))
}

#' @importFrom purrr map2 compact
#' @importFrom readr read_delim cols
#' @importFrom dplyr filter select mutate across everything
#' @importFrom magrittr %>%
#' @noRd
load_phenotype_data <- function(phenotype_path, region, extract_region_name = NULL, region_name_col = NULL, tabix_header = TRUE) {
  if (is.null(extract_region_name)) {
    extract_region_name <- rep(list(NULL), length(phenotype_path))
  } else if (is.list(extract_region_name) && length(extract_region_name) != length(phenotype_path)) {
    stop("extract_region_name must be NULL or a list with the same length as phenotype_path.")
  } else if (!is.null(extract_region_name) && !is.list(extract_region_name)) {
    stop("extract_region_name must be NULL or a list.")
  }

  # Use `map2` to iterate over `phenotype_path` and `extract_region_name` simultaneously
  phenotype_data_raw <- map2(phenotype_path, extract_region_name, ~ {
    tabix_data <- if (!is.null(region)) tabix_region(.x, region, tabix_header = tabix_header) else read_delim(.x, "\t", col_types = cols())
    if (nrow(tabix_data) == 0) {
      message(paste("Phenotype file ", .x, " is empty for the specified region", if (is.null(region)) "" else region))
      return(NULL)
    }
    if (!is.null(.y) && is.vector(.y) && !is.null(region_name_col) && (region_name_col %% 1 == 0)) {
      if (region_name_col <= ncol(tabix_data)) {
        region_col_name <- colnames(tabix_data)[region_name_col]
        tabix_data <- tabix_data %>%
          filter(.data[[region_col_name]] %in% .y) %>%
          t()
        colnames(tabix_data) <- tabix_data[region_name_col, ]
        return(tabix_data)
      } else {
        stop("region_name_col is out of bounds for the number of columns in tabix_data.")
      }
    } else {
      result <- tabix_data %>% t()
      # Assign region names from region_name_col if available
      if (!is.null(region_name_col) && (region_name_col %% 1 == 0) && region_name_col <= ncol(tabix_data)) {
        colnames(result) <- tabix_data[[region_name_col]]
      }
      return(result)
    }
  })

  # Track which indices had non-NULL data, then remove NULLs
  kept_indices <- which(vapply(phenotype_data_raw, Negate(is.null), logical(1)))
  phenotype_data <- phenotype_data_raw[kept_indices]

  # Check if all phenotype files are empty
  if (length(phenotype_data) == 0) {
    stop(NoPhenotypeError(paste("All phenotype files are empty for the specified region", if (!is.null(region)) "" else region)))
  }
  # Store kept indices as attribute so callers can align covariates/conditions
  attr(phenotype_data, "kept_indices") <- kept_indices
  return(phenotype_data)
}

#' @importFrom purrr map
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @noRd
extract_phenotype_coordinates <- function(phenotype_list) {
  return(map(phenotype_list, ~ t(.x[1:3, ]) %>%
    as_tibble() %>%
    mutate(start = as.numeric(start), end = as.numeric(end))))
}

#' @importFrom magrittr %>%
#' @noRd
filter_by_common_samples <- function(dat, common_samples) {
  dat[common_samples, , drop = FALSE] %>% .[order(rownames(.)), , drop = FALSE]
}

#' @importFrom tibble tibble
#' @importFrom dplyr mutate select
#' @importFrom purrr map map2
#' @importFrom magrittr %>%
#' @noRd
prepare_data_list <- function(geno_bed, phenotype, covariate, imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff, phenotype_header = 4, keep_samples = NULL) {
  data_list <- tibble(
    covar = covariate,
    Y = lapply(phenotype, function(x) apply(x[-c(1:phenotype_header), , drop = FALSE], c(1, 2), as.numeric))
  ) %>%
    mutate(
      # Determine common complete samples across Y, covar, and geno_bed, considering missing values
      common_complete_samples = map2(covar, Y, ~ {
        covar_non_na <- rownames(.x)[!apply(.x, 1, function(row) all(is.na(row)))]
        y_non_na <- rownames(.y)[!apply(.y, 1, function(row) all(is.na(row)))]
        if (length(intersect(intersect(covar_non_na, y_non_na), rownames(geno_bed))) == 0) {
          stop("No common complete samples between genotype and phenotype/covariate data")
        }
        intersect(intersect(covar_non_na, y_non_na), rownames(geno_bed))
      }),
      # Further intersect with keep_samples if provided
      common_complete_samples = if (!is.null(keep_samples) && length(keep_samples) > 0) {
        map(common_complete_samples, ~ intersect(.x, keep_samples))
      } else {
        common_complete_samples
      },
      # Determine dropped samples before filtering
      dropped_samples_covar = map2(covar, common_complete_samples, ~ setdiff(rownames(.x), .y)),
      dropped_samples_Y = map2(Y, common_complete_samples, ~ setdiff(rownames(.x), .y)),
      dropped_samples_X = map(common_complete_samples, ~ setdiff(rownames(geno_bed), .x)),
      # Filter data based on common complete samples
      Y = map2(Y, common_complete_samples, ~ filter_by_common_samples(.x, .y)),
      covar = map2(covar, common_complete_samples, ~ filter_by_common_samples(.x, .y)),
      # Apply filter_X on the geno_bed data filtered by common complete samples and then format column names
      X = map(common_complete_samples, ~ {
        filtered_geno_bed <- filter_by_common_samples(geno_bed, .x)
        mac_val <- if (nrow(filtered_geno_bed) == 0) 0 else (mac_cutoff / (2 * nrow(filtered_geno_bed)))
        maf_val <- max(maf_cutoff, mac_val)
        filtered_data <- filter_X(filtered_geno_bed, imiss_cutoff, maf_val, var_thresh = xvar_cutoff)
        colnames(filtered_data) <- normalize_variant_id(colnames(filtered_data)) # Normalize to canonical format
        filtered_data
      })
    ) %>%
    select(covar, Y, X, dropped_samples_Y, dropped_samples_X, dropped_samples_covar)
  return(data_list)
}

#' @importFrom purrr map
#' @importFrom dplyr intersect
#' @importFrom stringr str_split_fixed
#' @importFrom magrittr %>%
#' @noRd
prepare_X_matrix <- function(geno_bed, data_list, imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff) {
  # Calculate the union of all samples from data_list: any of X, covar and Y would do
  all_samples_union <- map(data_list$covar, ~ rownames(.x)) %>%
    unlist() %>%
    unique()
  # Find the intersection of these samples with the samples in geno_bed
  common_samples <- intersect(all_samples_union, rownames(geno_bed))
  # Filter geno_bed using common_samples
  X_filtered <- filter_by_common_samples(geno_bed, common_samples)
  # Calculate MAF cutoff considering the number of common samples
  maf_val <- max(maf_cutoff, mac_cutoff / (2 * length(common_samples)))
  # Apply further filtering on X
  X_filtered <- filter_X(X_filtered, imiss_cutoff, maf_val, xvar_cutoff)
  colnames(X_filtered) <- normalize_variant_id(colnames(X_filtered))

  # To keep a log message
  variants <- str_split_fixed(colnames(X_filtered), ":", 3)
  message(paste0("Dimension of input genotype data is ", nrow(X_filtered), " rows and ", ncol(X_filtered), " columns for genomic region of ", variants[1, 1], ":", min(as.integer(variants[, 2])), "-", max(as.integer(variants[, 2]))))
  return(X_filtered)
}

#' @importFrom purrr map map2
#' @importFrom dplyr mutate
#' @importFrom stats lm.fit sd
#' @importFrom magrittr %>%
#' @noRd
add_X_residuals <- function(data_list, scale_residuals = FALSE) {
  # Compute residuals for X and add them to data_list
  data_list <- data_list %>%
    mutate(
      lm_res_X = map2(X, covar, ~ .lm.fit(x = cbind(1, .y), y = .x)$residuals %>% as.matrix()),
      X_resid_mean = map(lm_res_X, ~ apply(.x, 2, mean)),
      X_resid_sd = map(lm_res_X, ~ apply(.x, 2, sd)),
      X_resid = map(lm_res_X, ~ {
        if (scale_residuals) {
          scale(.x)
        } else {
          .x
        }
      })
    )

  return(data_list)
}

#' @importFrom purrr map map2
#' @importFrom dplyr mutate
#' @importFrom stats lm.fit sd
#' @importFrom magrittr %>%
#' @noRd
add_Y_residuals <- function(data_list, conditions, scale_residuals = FALSE) {
  # Compute residuals, their mean, and standard deviation, and add them to data_list
  data_list <- data_list %>%
    mutate(
      lm_res = map2(Y, covar, ~ {
        res <- .lm.fit(x = cbind(1, .y), y = .x)$residuals %>% as.matrix()
        colnames(res) <- colnames(.x)
        res
      }),
      Y_resid_mean = map(lm_res, ~ apply(.x, 2, mean)),
      Y_resid_sd = map(lm_res, ~ apply(.x, 2, sd)),
      Y_resid = map(lm_res, ~ {
        if (scale_residuals) {
          scale(.x)
        } else {
          .x
        }
      })
    )

  names(data_list$Y_resid) <- conditions

  return(data_list)
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
#' @return A list containing the following components:
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
#'
#' @export
load_regional_association_data <- function(genotype, # PLINK file
                                           phenotype, # a vector of phenotype file names
                                           covariate, # a vector of covariate file names corresponding to the phenotype file vector
                                           region, # a string of chr:start-end for phenotype region
                                           conditions, # a vector of strings
                                           maf_cutoff = 0,
                                           mac_cutoff = 0,
                                           xvar_cutoff = 0,
                                           imiss_cutoff = 0,
                                           association_window = NULL,
                                           extract_region_name = NULL,
                                           region_name_col = NULL,
                                           keep_indel = TRUE,
                                           keep_samples = NULL,
                                           keep_variants = NULL,
                                           phenotype_header = 4, # skip first 4 rows of transposed phenotype for chr, start, end and ID
                                           scale_residuals = FALSE,
                                           tabix_header = TRUE) {
  ## Load genotype
  geno <- load_genotype_region(genotype, association_window, keep_indel, keep_variants_path = keep_variants)
  ## Load phenotype and covariates and perform some pre-processing
  covar <- load_covariate_data(covariate)
  pheno <- load_phenotype_data(phenotype, region, extract_region_name = extract_region_name, region_name_col = region_name_col, tabix_header = tabix_header)
  # Align covariates and conditions with phenotypes after filtering
  # load_phenotype_data removes empty phenotypes and stores which indices survived
  kept_idx <- attr(pheno, "kept_indices")
  if (!is.null(kept_idx) && length(kept_idx) < length(covar)) {
    covar <- covar[kept_idx]
    if (!is.null(conditions)) conditions <- conditions[kept_idx]
  }
  ### including Y ( cov ) and specific X and covar match, filter X variants based on the overlapped samples.
  data_list <- prepare_data_list(geno, pheno, covar, imiss_cutoff,
    maf_cutoff, mac_cutoff, xvar_cutoff,
    phenotype_header = phenotype_header, keep_samples = keep_samples
  )
  maf_list <- setNames(lapply(data_list$X, function(x) apply(x, 2, compute_maf)), colnames(data_list$X))
  ## Get residue Y for each of condition and its mean and sd
  data_list <- add_Y_residuals(data_list, conditions, scale_residuals)
  ## Get residue X for each of condition and its mean and sd
  data_list <- add_X_residuals(data_list, scale_residuals)
  # Get X matrix for union of samples.
  # Short-circuit when there's only one condition: the per-condition X computed in
  # prepare_data_list already operates on the same sample set (the single condition's
  # common_complete_samples, which is itself a subset of rownames(geno)) and applies
  # the same MAF/imiss/var thresholds with the same MAC cutoff scaling, so the union
  # X is bit-equivalent to data_list$X[[1]]. Skipping the redundant filter_X call saves
  # work and avoids a duplicate "N out of M total variants dropped" log line.
  if (length(data_list$X) == 1) {
    X <- data_list$X[[1]]
    variants <- str_split_fixed(colnames(X), ":", 3)
    message(paste0(
      "Dimension of input genotype data is ", nrow(X), " rows and ",
      ncol(X), " columns for genomic region of ",
      variants[1, 1], ":", min(as.integer(variants[, 2])), "-",
      max(as.integer(variants[, 2]))
    ))
  } else {
    X <- prepare_X_matrix(geno, data_list, imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff)
  }
  parsed_region <- if (!is.null(region)) parse_region(region) else NULL
  ## residual_Y: a list of y either vector or matrix (CpG for example), and they need to match with residual_X in terms of which samples are missing.
  ## residual_X: is a list of R conditions each is a matrix, with list names being the names of conditions, column names being SNP names and row names being sample names.
  ## X: is the somewhat original genotype matrix output from `filter_X`, with column names being SNP names and row names being sample names. Sample names of X should match example sample names of residual_Y matrix form (not list); but the matrices inside residual_X would be subsets of sample name of residual_Y matrix form (not list).
  return(list(
    residual_Y = data_list$Y_resid,
    residual_X = data_list$X_resid,
    residual_Y_scalar = if (scale_residuals) data_list$Y_resid_sd else rep(1, length(data_list$Y_resid)),
    residual_X_scalar = if (scale_residuals) data_list$X_resid_sd else rep(1, length(data_list$X_resid)),
    dropped_sample = list(X = data_list$dropped_samples_X, Y = data_list$dropped_samples_Y, covar = data_list$dropped_samples_covar),
    covar = data_list$covar,
    Y = data_list$Y,
    X_data = data_list$X,
    X = X,
    maf = maf_list,
    chrom = if (!is.null(parsed_region)) paste0("chr", parsed_region$chrom) else NULL,
    grange = if (!is.null(parsed_region)) as.character(c(parsed_region$start, parsed_region$end)) else NULL,
    Y_coordinates = if (!is.null(region)) extract_phenotype_coordinates(pheno) else NULL
  ))
}

#' Load Regional Univariate Association Data
#'
#' This function loads regional association data for univariate analysis.
#' It includes residual matrices, original genotype data, and additional metadata.
#'
#' @importFrom matrixStats colVars
#' @return A list
#' @export
load_regional_univariate_data <- function(...) {
  dat <- load_regional_association_data(...)
  return(list(
    residual_Y = dat$residual_Y,
    residual_X = dat$residual_X,
    residual_Y_scalar = dat$residual_Y_scalar,
    residual_X_scalar = dat$residual_X_scalar,
    dropped_sample = dat$dropped_sample,
    maf = dat$maf,
    X = dat$X, # X unadjusted by covariate
    chrom = dat$chrom,
    grange = dat$grange,
    X_variance = lapply(dat$residual_X, function(x) colVars(x))
  ))
}

#' Load Regional Data for Regression Modeling
#'
#' This function loads regional association data formatted for regression modeling.
#' It includes phenotype, genotype, and covariate matrices along with metadata.
#'
#' @return A list
#' @export
load_regional_regression_data <- function(...) {
  dat <- load_regional_association_data(...)
  return(list(
    Y = dat$Y,
    X_data = dat$X_data,
    covar = dat$covar,
    dropped_sample = dat$dropped_sample,
    maf = dat$maf,
    chrom = dat$chrom,
    grange = dat$grange
  ))
}

# return matrix of R conditions, with column names being the names of the conditions (phenotypes) and row names being sample names. Even for one condition it has to be a matrix with just one column.
#' @noRd
pheno_list_to_mat <- function(data_list) {
  all_row_names <- unique(unlist(lapply(data_list$residual_Y, rownames)))
  # Step 2: Align matrices and fill with NA where necessary
  aligned_mats <- lapply(data_list$residual_Y, function(mat) {
    ### change the ncol of each matrix
    expanded_mat <- matrix(NA, nrow = length(all_row_names), ncol = ncol(mat), dimnames = list(all_row_names, colnames(mat)))
    common_rows <- intersect(rownames(mat), all_row_names)
    expanded_mat[common_rows, ] <- mat[common_rows, ]
    return(expanded_mat)
  })
  Y_resid_matrix <- do.call(cbind, aligned_mats)
  if (!is.null(names(data_list$residual_Y))) {
    colnames(Y_resid_matrix) <- names(data_list$residual_Y)
  }
  data_list$residual_Y <- Y_resid_matrix
  return(data_list)
}

#' Load and Preprocess Regional Multivariate Data
#'
#' This function loads regional association data and processes it into a multivariate format.
#' It optionally filters out samples based on missingness thresholds in the response matrix.
#'
#' @importFrom matrixStats colVars
#' @return A list
#' @export
load_regional_multivariate_data <- function(matrix_y_min_complete = NULL, # when Y is saved as matrix, remove those with non-missing counts less than this cutoff
                                            ...) {
  dat <- pheno_list_to_mat(load_regional_association_data(...))
  if (!is.null(matrix_y_min_complete)) {
    Y <- filter_Y(dat$residual_Y, matrix_y_min_complete)
    if (length(Y$rm_rows) > 0) {
      X <- dat$X[-Y$rm_rows, ]
      Y_scalar <- dat$residual_Y_scalar[-Y$rm_rows]
      dropped_sample <- rownames(dat$residual_Y)[Y$rm_rows]
    } else {
      X <- dat$X
      Y_scalar <- dat$residual_Y_scalar
      dropped_sample <- dat$dropped_sample
    }
  } else {
    Y <- dat$residual_Y
    X <- dat$X
    Y_scalar <- dat$residual_Y_scalar
    dropped_sample <- dat$dropped_sample
  }
  return(list(
    residual_Y = Y,
    residual_Y_scalar = Y_scalar,
    dropped_sample = dropped_sample,
    X = X,
    maf = apply(X, 2, compute_maf),
    chrom = dat$chrom,
    grange = dat$grange,
    X_variance = colVars(X)
  ))
}

#' Load Regional Functional Association Data
#'
#' This function loads precomputed regional functional association data.
#'
#' @param min_markers Minimum number of phenotype markers required for a study.
#'   If \code{NULL}, no marker-count filtering is applied.
#' @return A list
#' @export
load_regional_functional_data <- function(..., min_markers = NULL) {
  dat <- load_regional_association_data(...)
  if (!is.null(min_markers)) {
    dat <- .filter_functional_data_by_marker_count(dat, min_markers)
  }
  dat
}

.filter_functional_data_by_marker_count <- function(fdat, min_markers,
                                                    always_keep = c("dropped_sample", "dropped_samples", "X", "chrom", "grange")) {
  if (is.null(fdat$Y_coordinates)) return(fdat)
  keep <- vapply(fdat$Y_coordinates, function(x) nrow(x) >= min_markers, logical(1))
  filter_names <- setdiff(names(fdat), always_keep)
  fdat[filter_names] <- lapply(fdat[filter_names], function(x) {
    if (length(x) == length(keep)) x[keep] else x
  })
  fdat
}



# Function to remove gene name at the end of context name
#' @export
clean_context_names <- function(context, gene) {
  # Remove gene name if it matches the last part of the context
  gene <- gene[order(-nchar(unique(gene)))]
  for (gene_id in gene) {
    context <- gsub(paste0("_", gene_id), "", context)
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
#' consolidated_weights <- load_twas_weights(weight_db_files, condition, region, variable_name_obj)
#' print(consolidated_weights)
#' @export
load_twas_weights <- function(weight_db_files, conditions = NULL,
                              variable_name_obj = c("preset_variants_result", "variant_names"),
                              susie_obj = c("preset_variants_result", "susie_result_trimmed"),
                              twas_weights_table = "twas_weights") {
  ## Internal function to load and validate data from RDS files
  load_and_validate_data <- function(weight_db_files, conditions, variable_name_obj) {
    all_data <- do.call(c, lapply(unname(weight_db_files), function(rds_file) {
      db <- readRDS(rds_file)
      gene <- names(db)
      # Filter by conditions if specified
      if (!is.null(conditions)) {
        # Split contexts if specified and trim whitespace, cen handle single or multiple conditions
        conditions <- trimws(strsplit(conditions, ",")[[1]])

        # Filter the gene's data to only include specified context layers
        if (length(gene) == 1 && gene != "mnm_rs") { # Need check
          available_contexts <- names(db[[gene]])
          matching_contexts <- available_contexts[available_contexts %in% paste0(conditions, "_", gene)]
          if (length(matching_contexts) == 0) {
            warning(paste0("No matching context layers found in ", rds_file, ". Skipping this file."))
            return(NULL)
          }

          db[[gene]] <- db[[gene]][matching_contexts]
        }
      } else {
        # Set default for 'conditions' if they are not specified
        conditions <- names(db[[gene]])
      }
      if (any(unique(names(find_data(db, c(3, "twas_weights")))) %in% c("mrmash_weights", "mvsusie_weights"))) {
        names(db[[1]]) <- clean_context_names(names(db[[1]]), gene = gene)
        db <- list(mnm_rs = db[[1]])
      } else {
        # Check if region from all RDS files are the same
        if (length(gene) != 1) {
          stop("More than one region provided in the RDS file. ")
        } else {
          names(db[[gene]]) <- clean_context_names(names(db[[gene]]), gene = gene)
        }
      }
      return(db)
    }))
    # Remove NULL entries (from files that had no matching context layers)
    all_data <- all_data[!sapply(all_data, is.null)]

    if (length(all_data) == 0) {
      stop("No data loaded. Check that conditions match available context layers in the RDS files.")
    }
    # Combine the lists with the same region name
    gene <- unique(names(all_data)[!names(all_data) %in% "mnm_rs"])
    if (length(gene) > 1) stop("More than one region of twas weights data provided. ")
    combined_all_data <- lapply(split(all_data, names(all_data)), function(lst) {
      if (length(lst) > 1) {
        lst <- do.call(c, unname(lst))
      }
      if (isTRUE(names(lst) == "mnm_rs")) lst <- lst[[1]]
      if (gene %in% names(lst)) lst <- do.call(c, lapply(unname(lst), function(x) x))
      return(lst)
    })

    # merge univariate and multivariate results for same gene-context pair
    if ("mnm_rs" %in% names(combined_all_data)) {
      # gene <- names(combined_all_data)[!names(combined_all_data) %in% "mnm_rs"]
      overl_contexts <- names(combined_all_data[["mnm_rs"]])[names(combined_all_data[["mnm_rs"]]) %in% names(combined_all_data[[gene]])]
      multi_variants <- unique(find_data(combined_all_data$mnm_rs, c(2, variable_name_obj)))
      for (context in overl_contexts) {
        uni_variants <- get_nested_element(combined_all_data[[gene]][[context]], variable_name_obj)
        multi_weights <- setNames(rep(0, length(uni_variants)), uni_variants)
        multi_weights <- lapply(combined_all_data[["mnm_rs"]][[context]]$twas_weights, function(weight_list) {
          aligned_weights <- setNames(rep(0, length(uni_variants)), uni_variants)
          method_weight_variants <- names(unlist(weight_list))
          overlap_variants <- method_weight_variants[method_weight_variants %in% multi_variants[multi_variants %in% uni_variants]] # overlapping variants from method, multivariate, univariate
          aligned_weights[overlap_variants] <- unlist(weight_list)[overlap_variants]
          aligned_weights <- as.matrix(aligned_weights)
        })
        combined_all_data[[gene]][[context]]$twas_weights <- c(combined_all_data[[gene]][[context]]$twas_weights, multi_weights)
        combined_all_data[[gene]][[context]]$twas_cv_result$performance <- c(
          combined_all_data[[gene]][[context]]$twas_cv_result$performance,
          combined_all_data[["mnm_rs"]][[context]]$twas_cv_result$performance
        )
      }
      combined_all_data[["mnm_rs"]] <- NULL
    }
    if (gene %in% names(combined_all_data)) combined_all_data <- do.call(c, unname(combined_all_data))
    if (gene %in% names(combined_all_data)) combined_all_data <- combined_all_data[[1]]

    # ## Check if the specified condition and variable_name_obj are available in all files
    # if (!all(conditions %in% names(combined_all_data))) {
    #   stop("The specified condition is not available in all RDS files.")
    # }
    return(combined_all_data)
  }

  # Internal function to align and merge weight matrices
  align_and_merge <- function(weights_list, variable_objs) {
    if (length(weights_list) != length(variable_objs)) {
      stop("The length of the weights_list and variable_objs must be the same.")
    }
    # Loop through each weight matrix and assign variant names as rownames
    for (i in seq_along(weights_list)) {
      # Ensure dimensions match
      if (nrow(weights_list[[i]]) != length(variable_objs[[i]])) {
        stop(paste("Number of rows in weights_list[[", i, "]] does not match the length of variable_objs[[", i, "]]", sep = ""))
      }
      # Apply variant names to the row names of the weight matrix
      rownames(weights_list[[i]]) <- variable_objs[[i]]
    }
    return(weights_list)
  }

  # Internal function to consolidate weights for given condition
  consolidate_weights_list <- function(combined_all_data, conditions, variable_name_obj, twas_weights_table) {
    combined_weights_by_condition <- lapply(conditions, function(condition) {
      temp_list <- get_nested_element(combined_all_data, c(condition, twas_weights_table))
      sapply(temp_list, cbind)
    })
    names(combined_weights_by_condition) <- conditions
    if (is.null(variable_name_obj)) {
      # Standard processing: Check for identical row numbers and consolidate
      row_numbers <- sapply(combined_weights_by_condition, function(data) nrow(data))
      if (length(unique(row_numbers)) > 1) {
        stop("Not all files have the same number of rows for the specified condition.")
      }
      weights <- combined_weights_by_condition
    } else {
      # Processing with variable_name_obj: Align and merge data, fill missing with zeros
      variable_objs <- lapply(conditions, function(condition) {
        get_nested_element(combined_all_data, c(condition, variable_name_obj))
      })
      weights <- align_and_merge(combined_weights_by_condition, variable_objs)
    }
    names(weights) <- conditions
    return(weights)
  }

  ## Load, validate, and consolidate data
  try(
    {
      combined_all_data <- load_and_validate_data(weight_db_files, conditions, variable_name_obj)
      if (is.null(combined_all_data)) {
        return(NULL)
      }
      # update condition in case of merging rds files
      conditions <- names(combined_all_data)
      weights <- consolidate_weights_list(combined_all_data, conditions, variable_name_obj, twas_weights_table)
      combined_susie_result <- lapply(combined_all_data, function(context) get_nested_element(context, susie_obj))
      performance_tables <- lapply(conditions, function(condition) {
        get_nested_element(combined_all_data, c(condition, "twas_cv_result", "performance"))
      })
      names(performance_tables) <- conditions
      return(list(susie_results = combined_susie_result, weights = weights, twas_cv_performance = performance_tables))
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
#' @param column_file_path Optional file path to a custom column mapping file
#'   (format: standard_name:original_name, one per line). Applied after
#'   MungeSumstats standardization.
#' @param comment_string Comment character in column_file_path. Default is "#".
#' @return A data frame with standardized column names.
#' @export
standardise_sumstats_columns <- function(sumstats, column_file_path = NULL, comment_string = "#") {
  # MungeSumstats standard names -> pecotmr conventions
  ms_to_pecotmr <- c(
    CHR = "chrom", BP = "pos", SNP = "variant_id",
    BETA = "beta", SE = "se", Z = "z", P = "p",
    N = "n_sample", N_CAS = "n_case", N_CON = "n_control",
    FRQ = "maf"
  )
  # Make a copy to avoid in-place modification by MungeSumstats
  sumstats_copy <- data.frame(sumstats, check.names = FALSE)
  # Use MungeSumstats for comprehensive column standardization
  sumstats_copy <- MungeSumstats::standardise_header(
    sumstats_copy, return_list = FALSE, uppercase_unmapped = FALSE
  )
  # Rename MungeSumstats standard names to pecotmr conventions
  for (ms_name in names(ms_to_pecotmr)) {
    idx <- which(colnames(sumstats_copy) == ms_name)
    if (length(idx) > 0) {
      colnames(sumstats_copy)[idx] <- ms_to_pecotmr[ms_name]
    }
  }
  # Apply additional custom column mapping if provided
  if (!is.null(column_file_path)) {
    if (!file.exists(column_file_path)) {
      stop("Column mapping file not found: ", column_file_path)
    }
    column_data <- read.table(column_file_path,
      header = FALSE, sep = ":",
      comment.char = if (is.null(comment_string)) "" else comment_string,
      stringsAsFactors = FALSE
    )
    colnames(column_data) <- c("standard", "original")
    for (i in seq_len(nrow(column_data))) {
      idx <- which(colnames(sumstats_copy) == column_data$original[i])
      if (length(idx) > 0) {
        colnames(sumstats_copy)[idx] <- column_data$standard[i]
      }
    }
  }
  as.data.frame(sumstats_copy)
}

#' Load summary statistic data
#'
#' This function formats the input summary statistics dataframe with uniform column names
#' to fit into the SuSiE pipeline. Column standardization is performed via
#' MungeSumstats::standardise_header(), with an optional custom column mapping file
#' for additional non-standard names.
#' Additionally, it extracts sample size, case number, control number, and variance of Y.
#' Missing values in n_sample, n_case, and n_control are backfilled with median values.
#'
#' @param sumstat_path File path to the summary statistics.
#' @param column_file_path Optional file path to a custom column mapping file for
#'   non-standard column names not recognized by MungeSumstats.
#' @param n_sample User-specified sample size. If unknown, set as 0 to retrieve from the sumstat file.
#' @param n_case User-specified number of cases.
#' @param n_control User-specified number of controls.
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
load_rss_data <- function(sumstat_path, column_file_path = NULL, n_sample = 0, n_case = 0, n_control = 0, region = NULL,
                          extract_region_name = NULL, region_name_col = NULL, comment_string = "#") {
  # Validate input files exist
  if (!file.exists(sumstat_path)) {
    stop("Summary statistics file not found: ", sumstat_path)
  }
  if (!is.null(column_file_path) && !file.exists(column_file_path)) {
    stop("Column mapping file not found: ", column_file_path)
  }
  var_y <- NULL
  sumstats <- load_tsv_region(file_path = sumstat_path, region = region, extract_region_name = extract_region_name, region_name_col = region_name_col)

  # To keep a log message
  n_variants <- nrow(sumstats)
  if (n_variants == 0) {
    message(paste0("No variants in region ", region, "."))
    return(list(sumstats = sumstats, n = NULL, var_y = NULL))
  } else {
    message(paste0("Region ", region, " include ", n_variants, " in input sumstats."))
  }

  # Standardize column names via MungeSumstats + optional custom mapping
  sumstats <- standardise_sumstats_columns(sumstats, column_file_path, comment_string)
  if (!"z" %in% colnames(sumstats) && all(c("beta", "se") %in%
    colnames(sumstats))) {
    sumstats$z <- sumstats$beta / sumstats$se
  }
  if (!"beta" %in% colnames(sumstats) && "z" %in% colnames(sumstats)) {
    sumstats$beta <- sumstats$z
    sumstats$se <- 1
  }
  for (col in c("n_sample", "n_case", "n_control")) {
    if (col %in% colnames(sumstats)) {
      sumstats[[col]][is.na(sumstats[[col]])] <- median(sumstats[[col]],
        na.rm = TRUE
      )
    }
  }
  if (n_sample != 0 && (n_case + n_control) != 0) {
    stop("Please provide sample size, or case number with control number, but not both")
  } else if (n_sample != 0) {
    n <- n_sample
  } else if ((n_case + n_control) != 0) {
    n <- n_case + n_control
    phi <- n_case / n
    var_y <- 1 / (phi * (1 - phi))
  } else {
    if ("n_sample" %in% colnames(sumstats)) {
      n <- median(sumstats$n_sample)
    } else if (all(c("n_case", "n_control") %in% colnames(sumstats))) {
      n <- median(sumstats$n_case + sumstats$n_control)
      phi <- median(sumstats$n_case / n)
      var_y <- 1 / (phi * (1 - phi))
    } else {
      warning("Sample size and variance of Y could not be determined from the summary statistics.")
      n <- NULL
    }
  }
  return(list(sumstats = sumstats, n = n, var_y = var_y))
}


#' This function loads a mixture data sets for a specific region, including individual-level data (genotype, phenotype, covariate data)
#' or summary statistics (sumstats, LD). Run \code{load_regional_univariate_data} and \code{load_rss_data} multiple times for different datasets
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
#' @param region The region where tabix use to subset the input dataset.
#' @param extract_sumstats_region_name User-specified gene/phenotype name used to further subset the phenotype data.
#' @param sumstats_region_name_col Filter this specific column for the extract_sumstats_region_name.
#' @param comment_string comment sign in the column_mapping file, default is #
#' @param extract_coordinates Optional data frame with columns "chrom" and "pos" for specific coordinates extraction.
#'
#' @return A list containing the individual_data and sumstat_data:
#' individual_data contains the following components if exist
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
#' sumstat_data contains the following components if exist
#' \itemize{
#'   \item sumstats: A list of summary statistics for the matched LD_info, each sublist contains sumstats, n, var_y from \code{load_rss_data}.
#'   \item LD_info: A list of LD information, each sublist contains LD_variants, LD_matrix, ref_panel  \code{load_LD_matrix}.
#' }
#'
#' @export
load_multitask_regional_data <- function(region, # a string of chr:start-end for phenotype region
                                         genotype_list = NULL, # PLINK file
                                         phenotype_list = NULL, # a vector of phenotype file names
                                         covariate_list = NULL, # a vector of covariate file names corresponding to the phenotype file vector
                                         conditions_list_individual = NULL, # a vector of strings
                                         match_geno_pheno = NULL, # a vector of index of phentoypes matched to genotype if mulitple genotype files
                                         maf_cutoff = 0,
                                         mac_cutoff = 0,
                                         xvar_cutoff = 0,
                                         imiss_cutoff = 0,
                                         association_window = NULL,
                                         extract_region_name = NULL,
                                         region_name_col = NULL,
                                         keep_indel = TRUE,
                                         keep_samples = NULL,
                                         keep_variants = NULL,
                                         phenotype_header = 4, # skip first 4 rows of transposed phenotype for chr, start, end and ID
                                         scale_residuals = FALSE,
                                         tabix_header = TRUE,
                                         # sumstat if need
                                         sumstat_path_list = NULL,
                                         column_file_path_list = NULL,
                                         LD_meta_file_path_list = NULL,
                                         match_LD_sumstat = NULL, # a vector of index of sumstat matched to LD if mulitple sumstat files
                                         conditions_list_sumstat = NULL,
                                         n_samples = 0,
                                         n_cases = 0,
                                         n_controls = 0,
                                         extract_sumstats_region_name = NULL,
                                         sumstats_region_name_col = NULL,
                                         comment_string = "#",
                                         extract_coordinates = NULL) {
  if (is.null(genotype_list) & is.null(sumstat_path_list)) {
    stop("Data load error. Please make sure at least one data set (sumstat_path_list or genotype_list) exists.")
  }

  # - if exist individual level data
  individual_data <- NULL
  if (!is.null(genotype_list)) {
    if (length(phenotype_list) != length(covariate_list)) {
      stop("Data load error. 'phenotype_list' and 'covariate_list' must have the same length.")
    }
    if (is.null(conditions_list_individual)) {
      conditions_list_individual <- paste0("condition", seq_along(phenotype_list))
      warning("Data load warning. 'conditions_list_individual' is not supplied; using default condition names. ",
              "Provide 'conditions_list_individual' to preserve cohort or cell-type labels.")
    }
    if (length(conditions_list_individual) != length(phenotype_list)) {
      stop("Data load error. 'conditions_list_individual' must have the same length as 'phenotype_list'.")
    }
    #### FIXME: later if we have mulitple genotype list
    if (length(genotype_list) != 1 & is.null(match_geno_pheno)) {
      stop("Data load error. Please make sure 'match_geno_pheno' exists if you load data from multiple individual-level data.")
    } else if (length(genotype_list) == 1 & is.null(match_geno_pheno)) {
      match_geno_pheno <- rep(1, length(phenotype_list))
    }
    if (length(match_geno_pheno) != length(phenotype_list)) {
      stop("Data load error. 'match_geno_pheno' must have the same length as 'phenotype_list'.")
    }
    if (any(is.na(match_geno_pheno)) ||
        any(match_geno_pheno < 1 | match_geno_pheno > length(genotype_list))) {
      stop("Data load error. 'match_geno_pheno' must contain valid indices into 'genotype_list'.")
    }

    # - load individual data from multiple datasets
    n_dataset <- unique(match_geno_pheno)
    for (i_data in n_dataset) {
      # extract genotype file name
      genotype <- genotype_list[i_data]
      # extract phenotype and covariate file names
      pos <- which(match_geno_pheno == i_data)
      phenotype <- phenotype_list[pos]
      covariate <- covariate_list[pos]
      conditions <- conditions_list_individual[pos]
      extract_region_name_i <- extract_region_name
      if (is.list(extract_region_name) && length(extract_region_name) == length(phenotype_list)) {
        extract_region_name_i <- extract_region_name[pos]
      }
      dat <- load_regional_univariate_data(
        genotype = genotype, phenotype = phenotype,
        covariate = covariate,
        region = region,
        association_window = association_window,
        conditions = conditions, xvar_cutoff = xvar_cutoff,
        maf_cutoff = maf_cutoff, mac_cutoff = mac_cutoff,
        imiss_cutoff = imiss_cutoff, keep_indel = keep_indel,
        keep_samples = keep_samples, keep_variants = keep_variants,
        extract_region_name = extract_region_name_i,
        phenotype_header = phenotype_header,
        region_name_col = region_name_col,
        scale_residuals = scale_residuals
      )
      if (is.null(individual_data)) {
        individual_data <- dat
      } else {
        individual_data <- stats::setNames(lapply(names(dat), function(k) {
          c(individual_data[[k]], dat[[k]])
        }), names(dat))
        individual_data$chrom <- dat$chrom
        individual_data$grange <- dat$grange
      }
    }
  }

  # - if exist summstat data
  sumstat_data <- NULL
  if (!is.null(sumstat_path_list)) {
    if (length(match_LD_sumstat) == 0) {
      match_LD_sumstat[[1]] <- conditions_list_sumstat
    }
    if (length(match_LD_sumstat) != length(LD_meta_file_path_list)) {
      stop("Please make sure 'match_LD_sumstat' matched 'LD_meta_file_path_list' if you load data from multiple sumstats.")
    }
    # - load sumstat data from multiple datasets
    n_LD <- length(match_LD_sumstat)
    for (i_ld in 1:n_LD) {
      # extract LD meta file path name
      LD_meta_file_path <- LD_meta_file_path_list[i_ld]
      LD_info <- load_LD_matrix(LD_meta_file_path,
        region = association_window,
        extract_coordinates = extract_coordinates,
        return_genotype = "auto"
      )
      # extract sumstat information
      conditions <- match_LD_sumstat[[i_ld]]
      pos <- match(conditions, conditions_list_sumstat)
      sumstats <- lapply(pos, function(ii) {
        sumstat_path <- sumstat_path_list[ii]
        column_file_path <- column_file_path_list[ii]
        # Load sumstat for this study (multiple LD references handled by outer loop)
        tmp <- load_rss_data(
          sumstat_path = sumstat_path, column_file_path = column_file_path,
          n_sample = n_samples[ii], n_case = n_cases[ii], n_control = n_controls[ii],
          region = association_window, extract_region_name = extract_sumstats_region_name,
          region_name_col = sumstats_region_name_col, comment_string = comment_string
        )
        if (nrow(tmp$sumstats) == 0){ return(NULL) }
        if (!("variant_id" %in% colnames(tmp$sumstats))) {
          tmp$sumstats <- tmp$sumstats %>%
            mutate(variant_id = format_variant_id(chrom, pos, A2, A1))
        }
        return(tmp)
      })
      names(sumstats) <- conditions
      if_no_variants <- sapply(sumstats, is.null)
      if (sum(if_no_variants)!=0){
        pos_no_variants <- which(if_no_variants)
        sumstats <- sumstats[-pos_no_variants]
      }
      sumstat_data$sumstats <- c(sumstat_data$sumstats, list(sumstats))
      sumstat_data$LD_info <- c(sumstat_data$LD_info, list(LD_info))
    }
    names(sumstat_data$sumstats) <- names(sumstat_data$LD_info) <- names(match_LD_sumstat)
  }

  return(list(
    individual_data = individual_data,
    sumstat_data = sumstat_data
  ))
}

#' Convert loaded regional data to individual-level inputs
#'
#' @param region_data A list returned by \code{load_multitask_regional_data()}.
#' @return A list containing \code{X}, \code{Y}, \code{maf},
#'   \code{X_variance}, and source information.
#' @export
region_data_to_ind_input <- function(region_data) {
  first_non_null <- function(...) {
    values <- list(...)
    for (value in values) {
      if (!is.null(value)) return(value)
    }
    NULL
  }

  align_individual_contexts <- function(X, Y) {
    cbind_y <- function(Y, fallback_names) {
      keep <- !vapply(Y, is.null, logical(1))
      if (!any(keep)) return(NULL)
      Y <- Y[keep]
      fallback_names <- fallback_names[keep]
      mats <- Map(function(y, nm) {
        if (is.null(dim(y))) y <- matrix(y, ncol = 1)
        if (is.null(colnames(y))) colnames(y) <- nm
        y
      }, Y, fallback_names)
      do.call(cbind, mats)
    }

    if (!is.list(X) || is.matrix(X) || is.data.frame(X) ||
        !is.list(Y) || is.matrix(Y) || is.data.frame(Y) ||
        is.null(names(X)) || is.null(names(Y)) ||
        length(intersect(names(X), names(Y))) > 0) {
      return(list(X = X, Y = Y))
    }
    x_names <- names(X)
    y_names <- names(Y)
    grouped <- list()
    for (context in x_names) {
      matched <- y_names[y_names == context | startsWith(y_names, paste0(context, "_"))]
      if (length(matched) > 0) {
        y_group <- cbind_y(Y[matched], matched)
        if (!is.null(y_group)) grouped[[context]] <- y_group
      }
    }
    if (length(grouped) == 0 && length(X) == 1 && length(Y) > 0) {
      y_group <- cbind_y(Y, y_names)
      if (!is.null(y_group)) grouped[[x_names[[1]]]] <- y_group
    }
    if (length(grouped) == 0) {
      return(list(X = X, Y = Y))
    }
    list(X = X[names(grouped)], Y = grouped)
  }

  individual_data <- region_data$individual_data
  if (is.null(individual_data)) {
    return(list(X = NULL, Y = NULL, maf = NULL, X_variance = NULL,
                source_info = list(has_individual = FALSE, contexts = character())))
  }

  X <- first_non_null(individual_data$residual_X, individual_data$X)
  Y <- first_non_null(individual_data$residual_Y, individual_data$Y)
  if (is.list(X) && !is.matrix(X) && !is.data.frame(X) &&
      is.null(names(X)) && !is.null(names(Y)) && length(X) == length(Y)) {
    names(X) <- names(Y)
  }
  if (is.list(Y) && !is.matrix(Y) && !is.data.frame(Y) &&
      is.null(names(Y)) && !is.null(names(X)) && length(Y) == length(X)) {
    names(Y) <- names(X)
  }
  if (is.matrix(X) && is.list(Y) && !is.null(names(Y))) {
    X <- stats::setNames(rep(list(X), length(Y)), names(Y))
  }
  aligned <- align_individual_contexts(X, Y)
  X <- aligned$X
  Y <- aligned$Y
  maf <- individual_data$maf
  X_variance <- individual_data$X_variance
  if (is.list(maf) && is.null(names(maf)) && !is.null(names(X)) && length(maf) == length(X)) {
    names(maf) <- names(X)
  }
  if (is.list(X_variance) && is.null(names(X_variance)) && !is.null(names(X)) &&
      length(X_variance) == length(X)) {
    names(X_variance) <- names(X)
  }
  contexts <- unique(c(names(X), names(Y)))
  list(
    X = X,
    Y = Y,
    maf = maf,
    X_variance = X_variance,
    source_info = list(has_individual = !is.null(X) && !is.null(Y),
                       contexts = contexts)
  )
}

#' Convert loaded regional data to RSS inputs
#'
#' @param region_data A list returned by \code{load_multitask_regional_data()}.
#' @return A list containing named RSS inputs, matched LD data, and source
#'   information.
#' @export
region_data_to_rss_input <- function(region_data) {
  make_ld_data_from_matrix <- function(ld, variant_ids = NULL) {
    is_genotype <- is.matrix(ld) && nrow(ld) != ncol(ld)
    if (!is.null(variant_ids) && is.matrix(ld)) {
      if (is.null(colnames(ld)) && length(variant_ids) == ncol(ld)) {
        colnames(ld) <- variant_ids
      }
      if (!is_genotype && is.null(rownames(ld)) && length(variant_ids) == nrow(ld)) {
        rownames(ld) <- variant_ids
      }
    }
    ids <- if (is.matrix(ld) && !is_genotype) rownames(ld) else colnames(ld)
    parsed <- NULL
    if (!is.null(ids) && length(ids) > 0) {
      parsed <- tryCatch(parse_variant_id(ids), error = function(e) NULL)
      if (!is.null(parsed)) {
        ids <- format_variant_id(parsed$chrom, parsed$pos, parsed$A2, parsed$A1)
        if (!is_genotype && is.matrix(ld)) rownames(ld) <- colnames(ld) <- ids
        if (is_genotype && is.matrix(ld)) colnames(ld) <- ids
        parsed$variant_id <- ids
      }
    }
    list(
      LD_matrix = ld,
      LD_variants = ids,
      ref_panel = parsed,
      block_metadata = if (!is_genotype && !is.null(parsed)) .infer_single_ld_block_metadata(parsed) else NULL,
      is_genotype = isTRUE(is_genotype)
    )
  }

  rss_input_from_qced_sumstat <- function(sumstat_data) {
    variant_ids_from_rss <- function(rss) {
      ss <- rss$sumstats
      if (is.null(ss)) return(character())
      if ("variant_id" %in% colnames(ss)) return(normalize_variant_id(as.character(ss$variant_id)))
      if (all(c("chrom", "pos", "A2", "A1") %in% colnames(ss))) {
        return(format_variant_id(ss$chrom, ss$pos, ss$A2, ss$A1))
      }
      character()
    }

    rss_input <- sumstat_data$sumstats
    LD_mat <- sumstat_data$LD_mat
    LD_match <- sumstat_data$LD_match
    studies <- names(rss_input)
    LD_data <- list()
    ld_group <- character()
    for (i in seq_along(studies)) {
      study <- studies[[i]]
      ld_name <- if (!is.null(LD_match) && length(LD_match) >= i) LD_match[[i]] else study
      if (is.null(ld_name) || is.na(ld_name) || !ld_name %in% names(LD_mat)) {
        ld_name <- names(LD_mat)[min(i, length(LD_mat))]
      }
      ld <- LD_mat[[ld_name]]
      rss <- rss_input[[study]]
      variant_ids <- variant_ids_from_rss(rss)
      LD_data[[study]] <- make_ld_data_from_matrix(ld, variant_ids)
      ld_group[[study]] <- ld_name
    }
    list(
      rss_input = rss_input,
      LD_data = LD_data,
      source_info = list(has_sumstat = length(rss_input) > 0,
                         studies = names(rss_input),
                         ld_group = ld_group)
    )
  }

  sumstat_data <- region_data$sumstat_data
  if (is.null(sumstat_data) || is.null(sumstat_data$sumstats)) {
    return(list(rss_input = list(), LD_data = list(),
                source_info = list(has_sumstat = FALSE, studies = character(),
                                   ld_group = character())))
  }
  if (!is.null(sumstat_data$LD_mat)) {
    return(rss_input_from_qced_sumstat(sumstat_data))
  }

  rss_input <- list()
  LD_data <- list()
  ld_group <- character()

  for (i in seq_along(sumstat_data$sumstats)) {
    studies <- sumstat_data$sumstats[[i]]
    ld_index <- min(i, length(sumstat_data$LD_info))
    group_name <- names(sumstat_data$LD_info)[ld_index]
    if (is.null(group_name) || is.na(group_name) || group_name == "") {
      group_name <- paste0("LD", ld_index)
    }
    for (study in names(studies)) {
      output_name <- study
      if (output_name %in% names(rss_input)) {
        output_name <- make.unique(c(names(rss_input), output_name))[length(rss_input) + 1]
      }
      rss_input[[output_name]] <- studies[[study]]
      LD_data[[output_name]] <- sumstat_data$LD_info[[ld_index]]
      ld_group[[output_name]] <- group_name
    }
  }

  list(
    rss_input = rss_input,
    LD_data = LD_data,
    source_info = list(has_sumstat = length(rss_input) > 0,
                       studies = names(rss_input),
                       ld_group = ld_group)
  )
}

#' Load and filter tabular data with optional region subsetting
#'
#' This function loads summary statistics data from tabular files (TSV, TXT).
#' For compressed (.gz) and tabix-indexed files, it can subset data by genomic region.
#' Additionally, it can filter results by a specified target value in a designated column.
#'
#' @param file_path Path to the summary statistics file.
#' @param region Genomic region for subsetting tabix-indexed files. Format: chr:start-end (e.g., "9:10000-50000").
#' @param extract_region_name Value to filter for in the specified filter column.
#' @param region_name_col Index of the column to apply the extract_region_name against.
#'
#' @return A dataframe containing the filtered summary statistics.
#'
#' @importFrom vroom vroom
#' @export
load_tsv_region <- function(file_path, region = NULL, extract_region_name = NULL, region_name_col = NULL) {
  sumstats <- NULL

  if (grepl("\\.gz$", file_path)) {
    if (!is.null(region)) {
      # Use Rsamtools to query the tabix-indexed file by region
      sumstats <- tryCatch({
        tbx <- Rsamtools::TabixFile(file_path)
        parsed <- parse_region(region)
        # Match chromosome naming convention in the tabix index
        chrom <- as.character(parsed$chrom)
        tbx_seqnames <- Rsamtools::seqnamesTabix(tbx)
        if (any(grepl("^chr", tbx_seqnames))) {
          chrom <- paste0("chr", chrom)
        }
        gr <- GenomicRanges::GRanges(
          seqnames = chrom,
          ranges = IRanges::IRanges(start = parsed$start, end = parsed$end)
        )
        lines <- Rsamtools::scanTabix(tbx, param = gr)[[1]]
        if (length(lines) == 0) return(NULL)

        # Get header for column names
        hdr <- Rsamtools::headerTabix(tbx)$header
        col_names_vec <- NULL
        if (length(hdr) > 0) {
          last_hdr <- hdr[length(hdr)]
          col_names_vec <- strsplit(sub("^#", "", last_hdr), "\t")[[1]]
        } else {
          header_con <- gzfile(file_path, "rt")
          first_line <- readLines(header_con, n = 1)
          close(header_con)
          first_fields <- strsplit(sub("^#", "", first_line), "\t")[[1]]
          header_tokens <- c("chrom", "chr", "#chrom", "pos", "bp", "snp",
                             "variant_id", "a1", "a2", "beta", "se", "z",
                             "p", "pvalue")
          if (any(tolower(first_fields) %in% header_tokens)) {
            col_names_vec <- first_fields
          }
        }

        txt <- paste(lines, collapse = "\n")
        if (!is.null(col_names_vec)) {
          as.data.frame(vroom::vroom(I(txt), delim = "\t", col_names = col_names_vec,
                                     show_col_types = FALSE))
        } else {
          as.data.frame(vroom::vroom(I(txt), delim = "\t", col_names = TRUE,
                                     show_col_types = FALSE))
        }
      }, error = function(e) {
        stop("Data read error. Please make sure this gz file is tabix-indexed and the specified filter column exists.")
      })
    } else {
      # No region specified - read the whole gz file
      sumstats <- as.data.frame(vroom::vroom(file_path, show_col_types = FALSE))
    }
  } else {
    warning("Not a tabix-indexed gz file, loading the entire dataset.")
    sumstats <- as.data.frame(vroom::vroom(file_path, show_col_types = FALSE))
  }

  # Apply name-based filter if specified
  if (!is.null(sumstats) && !is.null(extract_region_name) && !is.null(region_name_col)) {
    keep_index <- which(str_detect(sumstats[[region_name_col]], extract_region_name))
    sumstats <- sumstats[keep_index, ]
  }

  return(sumstats)
}

#' Split loaded twas_weights_results into batches based on maximum memory usage
#'
#' @param twas_weights_results List of loaded gene data by load_twas_weights()
#' @param meta_data_df Dataframe containing gene metadata with region_id and TSS columns
#' @param max_memory_per_batch Maximum memory per batch in MB (default: 750)
#' @return List of batches, where each batch contains a subset of twas_weights_results
#' @export
batch_load_twas_weights <- function(twas_weights_results, meta_data_df, max_memory_per_batch = 750) {
  gene_names <- names(twas_weights_results)
  if (length(gene_names) == 0) {
    message("No genes in twas_weights_results.")
    return(list())
  }

  gene_memory_df <- data.frame(
    gene_name = gene_names, memory_mb = sapply(gene_names, function(gene) {
      as.numeric(object.size(twas_weights_results[[gene]])) / (1024^2) # Get object size in bytes and convert to MB
    })
  )

  # Merge with meta_data_df to get TSS information
  meta_data_df <- meta_data_df[!duplicated(meta_data_df[, c("region_id", "TSS")]), ]
  gene_memory_df <- merge(gene_memory_df, meta_data_df[, c("region_id", "TSS")],
    by.x = "gene_name",
    by.y = "region_id", all.x = TRUE
  )
  gene_memory_df <- gene_memory_df[order(gene_memory_df$TSS), ]

  # Check if we need to split into batches
  total_memory_mb <- sum(gene_memory_df$memory_mb)
  message("Total memory usage: ", round(total_memory_mb, 2), " MB")
  if (total_memory_mb <= max_memory_per_batch) {
    message("All genes fit within the memory limit. No need to split into batches.")
    return(list(all_genes = twas_weights_results))
  }

  # Create batches by adding genes until we reach the memory limit
  batches <- list()
  current_batch_genes <- character(0)
  current_batch_memory <- 0
  batch_index <- 1

  for (i in 1:nrow(gene_memory_df)) {
    gene <- gene_memory_df$gene_name[i]
    gene_memory <- gene_memory_df$memory_mb[i]
    # If a single gene exceeds the memory limit, include it in its own batch
    if (gene_memory > max_memory_per_batch) {
      batches[[paste0("batch_", batch_index)]] <- twas_weights_results[gene]
      batch_index <- batch_index + 1
      next
    }
    # If adding this gene would exceed the memory limit, start a new batch
    if (current_batch_memory + gene_memory > max_memory_per_batch && length(current_batch_genes) > 0) {
      batches[[paste0("batch_", batch_index)]] <- twas_weights_results[current_batch_genes]
      current_batch_genes <- character(0)
      current_batch_memory <- 0
      batch_index <- batch_index + 1
    }
    current_batch_genes <- c(current_batch_genes, gene)
    current_batch_memory <- current_batch_memory + gene_memory
  }
  # Add the last batch if not empty
  if (length(current_batch_genes) > 0) {
    batches[[batch_index]] <- twas_weights_results[current_batch_genes]
  }
  message("Split into ", length(batches), " batches")
  names(batches) <- NULL
  return(batches)
}

# Function to filter a single credible set based on coverage and purity
#' @importFrom susieR susie_get_cs
#' @importFrom purrr map_lgl
#' @export
get_filter_lbf_index <- function(susie_obj, coverage = 0.5, size_factor = 0.5) {
  susie_obj$V <- NULL  # ensure no filtering by estimated prior

  # Get CS list with coverage
  cs_list <- susie_get_cs(susie_obj, coverage = coverage, dedup = FALSE)

  # Total number of variants
  total_variants <- ncol(susie_obj$alpha)

  # Maximum allowed CS size to be considered 'concentrated'
  max_size <- total_variants * coverage * size_factor

  # Identify which CSs are 'concentrated enough'
  keep_idx <- map_lgl(cs_list$cs, ~ length(.x) < max_size)

  # Extract the CS indices that pass the filter
  cs_index <- which(keep_idx) %>% names %>% gsub("L","", .) %>% as.numeric

  # Return filtered lbf_variable rows (one per CS)
  return(cs_index)
}
