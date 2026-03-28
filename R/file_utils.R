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
  af <- select(af, chrom, id, A2, A1, alt_freq, obs_ct)
  return(af)
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

  # --- Read variant info and determine region indices ---
  pvar <- pgenlibr::NewPvar(paths$pvar)
  on.exit(pgenlibr::ClosePvar(pvar), add = TRUE)
  all_variant_info <- read_pvar_info(pvar)

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
  colnames(psam)[1:2] <- c("FID", "IID")

  # --- Read genotype dosage via pgenlibr ---
  pgen <- pgenlibr::NewPgen(paths$pgen, pvar = pvar)
  on.exit(pgenlibr::ClosePgen(pgen), add = TRUE)
  X <- pgenlibr::ReadList(pgen, variant_subset = variant_idx, meanimpute = FALSE)
  rownames(X) <- psam$IID
  colnames(X) <- variant_info$id

  # --- Attach allele frequency from .afreq if available ---
  afreq <- read_afreq(prefix)
  if (!is.null(afreq)) {
    variant_info <- merge(variant_info, afreq[, c("id", "alt_freq", "obs_ct")],
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
  # pgenlibr::NewPvar handles both .pvar and .pvar.zst natively
  pvar <- if (file.exists(paste0(prefix, ".pvar.zst"))) {
    paste0(prefix, ".pvar.zst")
  } else if (file.exists(paste0(prefix, ".pvar"))) {
    paste0(prefix, ".pvar")
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

#' Build variant_info data.frame from a pgenlibr pvar object.
#' @return data.frame with columns: chrom, id, pos, A2 (REF), A1 (ALT).
#' @noRd
read_pvar_info <- function(pvar) {
  n <- pgenlibr::GetVariantCt(pvar)
  idx <- seq_len(n)
  data.frame(
    chrom = vapply(idx, function(i) pgenlibr::GetVariantChrom(pvar, i), character(1)),
    id    = vapply(idx, function(i) pgenlibr::GetVariantId(pvar, i), character(1)),
    pos   = vapply(idx, function(i) pgenlibr::GetVariantPos(pvar, i), integer(1)),
    # PLINK2: allele 1 = REF (A2), allele 2 = first ALT (A1)
    A2    = vapply(idx, function(i) pgenlibr::GetAlleleCode(pvar, i, 1L), character(1)),
    A1    = vapply(idx, function(i) pgenlibr::GetAlleleCode(pvar, i, 2L), character(1)),
    stringsAsFactors = FALSE
  )
}

#' Get variant information from any LD reference source without loading genotypes.
#'
#' Auto-detects the source type (PLINK2, PLINK1, or pre-computed LD metadata)
#' and returns variant metadata. For PLINK2, opens only the .pvar file.
#' For PLINK1, reads only the .bim file. No genotype data is loaded.
#'
#' @param source PLINK prefix or LD metadata file path.
#' @param region Region of interest: "chr:start-end" string or data.frame with
#'   chrom/start/end. If NULL, returns all variants.
#' @return A data.frame with columns: chrom, id, pos, A2, A1.
#'   May also include allele_freq, variance, n_nomiss depending on source.
#'
#' @importFrom vroom vroom
#' @export
get_ref_variant_info <- function(source, region = NULL) {
  resolved <- resolve_ld_source(source)

  if (resolved$type == "plink2") {
    paths <- resolve_plink2_paths(resolved$data_path)
    pvar <- pgenlibr::NewPvar(paths$pvar)
    on.exit(pgenlibr::ClosePvar(pvar), add = TRUE)
    info <- read_pvar_info(pvar)
    # Attach allele frequency from .afreq if available
    afreq <- read_afreq(resolved$data_path)
    if (!is.null(afreq)) {
      info$allele_freq <- afreq$alt_freq[match(info$id, afreq$id)]
    }
  } else if (resolved$type == "plink1") {
    bim <- read_bim(paste0(resolved$data_path, ".bed"))
    info <- data.frame(
      chrom = bim$chrom, id = bim$id, pos = bim$pos,
      A2 = bim$a0, A1 = bim$a1,
      stringsAsFactors = FALSE
    )
  } else {
    # Pre-computed LD: read bim files via metadata
    bim_paths <- get_regional_ld_meta(resolved$meta_path, region)$intersections$bim_file_paths
    info <- do.call(rbind, lapply(bim_paths, function(path) {
      df <- as.data.frame(vroom(path, col_names = FALSE, show_col_types = FALSE))
      out <- data.frame(
        chrom = df$X1, id = df$X2, pos = df$X4,
        A2 = df$X5, A1 = df$X6,
        stringsAsFactors = FALSE
      )
      if (ncol(df) >= 8) { out$variance <- df$X7; out$allele_freq <- df$X8 }
      if (ncol(df) >= 9) { out$n_nomiss <- df$X9 }
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
#' @noRd
match_variants_to_keep <- function(variant_info, keep_variants_path) {
  keep_raw <- as.data.frame(vroom(keep_variants_path, show_col_types = FALSE))
  # parse_variant_id handles chr prefix, underscore/colon formats, build suffixes,
  # and returns data.frame with integer chrom/pos and character A2/A1
  keep_variants <- parse_variant_id(
    if ("chrom" %in% names(keep_raw) & "pos" %in% names(keep_raw)) keep_raw else keep_raw[[1]]
  )
  vi_chrom <- as.integer(strip_chr_prefix(variant_info$chrom))
  paste0(vi_chrom, ":", variant_info$pos) %in% paste0(keep_variants$chrom, ":", keep_variants$pos)
}

#' @importFrom vroom vroom
#' @importFrom dplyr as_tibble mutate filter
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom stringr str_detect

tabix_region <- function(file, region, tabix_header = "auto", target = "", target_column_index = "") {
  if (!file.exists(file)) {
    stop("Input file does not exist: ", file)
  }
  cmd_output <- tryCatch(
    {
      use_col_names <- if (identical(tabix_header, FALSE)) FALSE else TRUE
      as.data.frame(vroom(pipe(paste0("tabix -h ", file, " ", region)),
                          delim = "\t", col_names = use_col_names, show_col_types = FALSE))
    },
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
    col_types <- list(col_character(), col_character(), col_guess(), col_integer(), col_guess(), col_guess())
    bim_data <- vroom(bim_file, col_names = FALSE, col_types = col_types)
    bim_data$X1 <- strip_chr_prefix(bim_data$X1)
    snp_ids <- bim_data$X2[bim_data$X1 == parsed$chrom &
                            bim_data$X4 >= parsed$start &
                            bim_data$X4 <= parsed$end]
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

#' Load genotype data for a specific region
#'
#' Auto-detects PLINK2 (.pgen/.pvar[.zst]/.psam) or PLINK1 (.bed/.bim/.fam) format
#' and loads genotype data accordingly. Returns a numeric dosage matrix.
#'
#' @param genotype Path to the genotype data file (without extension).
#' @param region The target region in the format "chr:start-end".
#' @param keep_indel Whether to keep indel SNPs.
#' @param keep_variants_path Path to a file listing variants to keep.
#' @return A numeric dosage matrix (rows=samples, cols=variants).
#'
#' @export
load_genotype_region <- function(genotype, region = NULL, keep_indel = TRUE, keep_variants_path = NULL) {
  resolved <- resolve_ld_source(genotype)
  if (resolved$type == "plink2") {
    return(load_plink2_data(resolved$data_path, region = region, keep_indel = keep_indel,
                            keep_variants_path = keep_variants_path)$X)
  }
  if (resolved$type == "plink1") {
    return(load_plink1_data(resolved$data_path, region = region, keep_indel = keep_indel,
                            keep_variants_path = keep_variants_path)$X)
  }
  stop("Genotype files not found. Expected PLINK2 (.pgen/.pvar[.zst]/.psam), ",
       "PLINK1 (.bed/.bim/.fam), or metadata TSV pointing to PLINK files at: ", genotype)
}

#' @importFrom purrr map
#' @importFrom readr read_delim cols
#' @importFrom dplyr select mutate across everything
#' @importFrom magrittr %>%
#' @noRd
load_covariate_data <- function(covariate_path) {
  # Validate all covariate files exist
  missing <- covariate_path[!file.exists(covariate_path)]
  if (length(missing) > 0) {
    stop("Covariate file(s) not found: ", paste(missing, collapse = ", "))
  }
  return(map(covariate_path, ~ read_delim(.x, "\t", col_types = cols()) %>%
    select(-1) %>%
    mutate(across(everything(), as.numeric)) %>%
    t()))
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
  variants <- as.data.frame(do.call(rbind, lapply(colnames(X_filtered), function(x) strsplit(x, ":")[[1]][1:2])), stringsAsFactors = FALSE)
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
  # Get X matrix for union of samples
  X <- prepare_X_matrix(geno, data_list, imiss_cutoff, maf_cutoff, mac_cutoff, xvar_cutoff)
  region <- if (!is.null(region)) unlist(strsplit(region, ":", fixed = TRUE))
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
    chrom = region[1],
    grange = if (!is.null(region)) unlist(strsplit(region[2], "-", fixed = TRUE)) else NULL,
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
#' @return A list
#' @export
load_regional_functional_data <- function(...) {
  dat <- load_regional_association_data(...)
  return(dat)
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

#' Load summary statistic data
#'
#' This function formats the input summary statistics dataframe with uniform column names
#' to fit into the SuSiE pipeline. The mapping is performed through the specified column file.
#' Additionally, it extracts sample size, case number, control number, and variance of Y.
#' Missing values in n_sample, n_case, and n_control are backfilled with median values.
#'
#' @param sumstat_path File path to the summary statistics.
#' @param column_file_path File path to the column file for mapping.
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
load_rss_data <- function(sumstat_path, column_file_path, n_sample = 0, n_case = 0, n_control = 0, region = NULL,
                          extract_region_name = NULL, region_name_col = NULL, comment_string = "#") {
  # Validate input files exist
  if (!file.exists(sumstat_path)) {
    stop("Summary statistics file not found: ", sumstat_path)
  }
  if (!file.exists(column_file_path)) {
    stop("Column mapping file not found: ", column_file_path)
  }
  # Read and preprocess column mapping
  if (is.null(comment_string)) {
    column_data <- read.table(column_file_path,
      header = FALSE, sep = ":",
      comment.char = "", # This tells R not to treat any character as comment
      stringsAsFactors = FALSE
    ) %>%
      rename(standard = V1, original = V2)
  } else {
    column_data <- read.table(column_file_path,
      header = FALSE, sep = ":",
      comment.char = comment_string,
      stringsAsFactors = FALSE
    ) %>%
      rename(standard = V1, original = V2)
  }
  # Initialize sumstats variable
  sumstats <- NULL
  var_y <- NULL
  sumstats <- load_tsv_region(file_path = sumstat_path, region = region, extract_region_name = extract_region_name, region_name_col = region_name_col)
  
  # To keep a log message
  n_variants <- nrow(sumstats)
  if (n_variants == 0){
      message(paste0("No variants in region ", region, "."))
      return(list(sumstats = sumstats, n = NULL, var_y = NULL))
  } else {
      message(paste0("Region ", region, " include ", n_variants, " in input sumstats."))
  }
  
  # Standardize column names based on mapping
  for (name in colnames(sumstats)) {
    if (name %in% column_data$original) {
      index <- which(column_data$original == name)
      colnames(sumstats)[colnames(sumstats) == name] <- column_data$standard[index]
    }
  }
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
    #### FIXME: later if we have mulitple genotype list
    if (length(genotype_list) != 1 & is.null(match_geno_pheno)) {
      stop("Data load error. Please make sure 'match_geno_pheno' exists if you load data from multiple individual-level data.")
    } else if (length(genotype_list) == 1 & is.null(match_geno_pheno)) {
      match_geno_pheno <- rep(1, length(phenotype_list))
    }

    # - load individual data from multiple datasets
    n_dataset <- unique(match_geno_pheno) ### FIXME
    for (i_data in 1:n_dataset) {
      # extract genotype file name
      genotype <- genotype_list[i_data]
      # extract phenotype and covariate file names
      pos <- which(match_geno_pheno == i_data)
      phenotype <- phenotype_list[pos]
      covariate <- covariate_list[pos]
      conditions <- conditions_list_individual[pos]
      dat <- load_regional_univariate_data(
        genotype = genotype, phenotype = phenotype,
        covariate = covariate, region = region,
        association_window = association_window,
        conditions = conditions, xvar_cutoff = xvar_cutoff,
        maf_cutoff = maf_cutoff, mac_cutoff = mac_cutoff,
        imiss_cutoff = imiss_cutoff, keep_indel = keep_indel,
        keep_samples = keep_samples, keep_variants = keep_variants,
        extract_region_name = extract_region_name,
        phenotype_header = phenotype_header,
        region_name_col = region_name_col,
        scale_residuals = scale_residuals
      )
      if (is.null(individual_data)) {
        individual_data <- dat
      } else {
        individual_data <- lapply(names(dat), function(k) {
          c(individual_data[[k]], dat[[k]])
        })
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
        extract_coordinates = extract_coordinates
      )
      # extract sumstat information
      conditions <- match_LD_sumstat[[i_ld]]
      pos <- match(conditions, conditions_list_sumstat)
      sumstats <- lapply(pos, function(ii) {
        sumstat_path <- sumstat_path_list[ii]
        column_file_path <- column_file_path_list[ii]
        # FIXME later: when consider multiple LD reference
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
  cmd <- NULL

  if (!is.null(region)) {
    if (grepl("^chr", region)) {
      region <- strip_chr_prefix(region)
    }
  }

  if (grepl("\\.gz$", file_path)) {
    if (is.null(sumstats) || nrow(sumstats) == 0) {
      # Determine the appropriate command based on provided parameters
      if (!is.null(extract_region_name) && !is.null(region) && !is.null(region_name_col)) {
        # Both region and filter specified
        cmd <- paste0(
          "zcat ", file_path, " | head -1 && tabix ", file_path, " ", region,
          " | awk '$", region_name_col, " ~ /", extract_region_name, "/'"
        )
      } else if (!is.null(extract_region_name) && is.null(region) && !is.null(region_name_col)) {
        # Only filter specified, no region
        cmd <- paste0("zcat ", file_path, " | awk '$", region_name_col, " ~ /", extract_region_name, "/'")
      } else if (!is.null(region) && (is.null(region_name_col) || is.null(extract_region_name))) {
        # Only region specified, no filter
        cmd <- paste0("zcat ", file_path, " | head -1 && tabix ", file_path, " ", region)
      } else {
        # Neither region nor filter specified - read gz directly
        cmd <- NULL
      }

      sumstats <- tryCatch(
        {
          if (is.null(cmd)) {
            as.data.frame(vroom(file_path, show_col_types = FALSE))
          } else {
            as.data.frame(vroom(pipe(cmd), delim = "\t", show_col_types = FALSE))
          }
        },
        error = function(e) {
          stop("Data read error. Please make sure this gz file is tabix-indexed and the specified filter column exists.")
        }
      )
    }
  } else {
    warning("Not a tabix-indexed gz file, loading the entire dataset.")
    sumstats <- as.data.frame(vroom(file_path, show_col_types = FALSE))
    # Apply filter if specified
    if (!is.null(extract_region_name) && !is.null(region_name_col)) {
      keep_index <- which(str_detect(sumstats[[region_name_col]], extract_region_name))
      sumstats <- sumstats[keep_index, ]
    }
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

