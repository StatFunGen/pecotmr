#' Deduplicate and sort genomic regions by chromosome and start position.
#' @importFrom dplyr distinct arrange
#' @importFrom magrittr %>%
#' @noRd
order_dedup_regions <- function(df) {
  df$chrom <- as.integer(strip_chr_prefix(df$chrom))
  df <- distinct(df, chrom, start, .keep_all = TRUE) %>%
    arrange(chrom, start)
  df
}

#' Find the first and last rows of genomic_data that overlap a query region.
#' Clamps the query to the available data range before searching.
#' @importFrom dplyr filter arrange slice
#' @noRd
find_intersection_rows <- function(genomic_data, region_chrom, region_start, region_end) {
  chrom_data <- genomic_data %>% filter(chrom == region_chrom)
  if (nrow(chrom_data) == 0) stop("No data for chromosome ", region_chrom)

  # Clamp query to available range
  region_start <- max(region_start, min(chrom_data$start))
  region_end   <- min(region_end,   max(chrom_data$end))

  start_row <- genomic_data %>%
    filter(chrom == region_chrom, start <= region_start, end > region_start) %>%
    slice(1)
  end_row <- genomic_data %>%
    filter(chrom == region_chrom, start < region_end, end >= region_end) %>%
    arrange(desc(end)) %>%
    slice(1)

  if (nrow(start_row) == 0 || nrow(end_row) == 0) {
    stop("Region ", region_chrom, ":", region_start, "-", region_end,
         " is not covered by any rows in the LD metadata.")
  }
  list(start_row = start_row, end_row = end_row)
}

#' Validate that start_row..end_row fully covers [region_start, region_end].
#' @noRd
validate_selected_region <- function(start_row, end_row, region_start, region_end) {
  if (start_row$start > region_start || end_row$end < region_end) {
    stop("Region ", region_start, "-", region_end, " is not fully covered by the LD metadata ",
         "(available: ", start_row$start, "-", end_row$end, ").")
  }
}

#' Extract values of a column for rows spanning the intersection range.
#' @noRd
extract_file_paths <- function(genomic_data, intersection_rows, column_to_extract) {
  if (!column_to_extract %in% names(genomic_data)) {
    stop("Column '", column_to_extract, "' not found in genomic data.")
  }
  idx <- which(genomic_data$chrom == intersection_rows$start_row$chrom &
               genomic_data$start >= intersection_rows$start_row$start &
               genomic_data$start <= intersection_rows$end_row$start)
  genomic_data[[column_to_extract]][idx]
}

#' Find LD blocks overlapping a query region from a metadata TSV file.
#'
#' @param ld_reference_meta_file TSV with columns chrom, start, end, path.
#'   The path column may be comma-separated: "ld_file,bim_file".
#' @param region "chr:start-end" string or data.frame with chrom/start/end.
#' @param complete_coverage_required If TRUE, error when the region extends
#'   beyond available LD blocks.
#' @return A list with: intersections (LD_file_paths, bim_file_paths),
#'   ld_meta_data, and parsed region.
#' @importFrom stringr str_split
#' @importFrom dplyr select
#' @importFrom vroom vroom
#' @noRd
get_regional_ld_meta <- function(ld_reference_meta_file, region, complete_coverage_required = FALSE) {
  genomic_data <- vroom(ld_reference_meta_file)
  region <- parse_region(region)
  # Set column names
  names(genomic_data) <- c("chrom", "start", "end", "path")
  names(region) <- c("chrom", "start", "end")

  # Treat start=0, end=0 as "covers all regions" (used for whole-chromosome PLINK files)
  whole_chrom <- genomic_data$start == 0 & genomic_data$end == 0
  if (any(whole_chrom)) genomic_data$end[whole_chrom] <- Inf

  # Order and deduplicate regions
  genomic_data <- order_dedup_regions(genomic_data)
  region <- order_dedup_regions(region)

  # Process file paths
  file_path <- genomic_data$path %>%
    str_split(",", simplify = TRUE) %>%
    data.frame() %>%
    `colnames<-`(if (ncol(.) == 2) c("LD_file_path", "bim_file_path") else c("LD_file_path"))

  genomic_data <- cbind(genomic_data, file_path) %>% select(-path)

  # Find intersection rows
  intersection_rows <- find_intersection_rows(genomic_data, region$chrom, region$start, region$end)

  # Validate region
  if (complete_coverage_required) {
    validate_selected_region(intersection_rows$start_row, intersection_rows$end_row, region$start, region$end)
  }

  # Extract file paths
  LD_paths <- find_valid_file_paths(ld_reference_meta_file, extract_file_paths(genomic_data, intersection_rows, "LD_file_path"))
  bim_paths <- if ("bim_file_path" %in% names(genomic_data)) {
    find_valid_file_paths(ld_reference_meta_file, extract_file_paths(genomic_data, intersection_rows, "bim_file_path"))
  } else {
    NULL
  }

  return(list(
    intersections = list(
      start_index = intersection_rows$start_row,
      end_index = intersection_rows$end_row,
      LD_file_paths = LD_paths,
      bim_file_paths = bim_paths
    ),
    ld_meta_data = genomic_data,
    region = region
  ))
}

#' Read a pre-computed LD matrix (.cor.xz) and its bim file, returning a
#' symmetric matrix with variants ordered by position.
#' @importFrom dplyr mutate
#' @importFrom utils read.table
#' @importFrom stats setNames
#' @noRd
process_LD_matrix <- function(LD_file_path, snp_file_path = NULL) {
  # Read .cor.xz matrix
  LD_file_con <- xzfile(LD_file_path)
  LD_matrix <- scan(LD_file_con, quiet = TRUE)
  close(LD_file_con)
  LD_matrix <- matrix(LD_matrix, ncol = sqrt(length(LD_matrix)), byrow = TRUE)

  # Auto-detect variant metadata file: .bim (PLINK1) or .pvar/.pvar.zst (PLINK2)
  if (is.null(snp_file_path)) {
    candidates <- paste0(LD_file_path, c(".bim", ".pvar", ".pvar.zst"))
    found <- candidates[file.exists(candidates)]
    if (length(found) == 0) stop("No variant file found for: ", LD_file_path,
                                  " (tried .bim, .pvar, .pvar.zst)")
    snp_file_path <- found[1]
  }

  LD_variants <- read_variant_metadata(snp_file_path)
  is_pvar <- !("gpos" %in% names(LD_variants))
  LD_variants <- LD_variants %>%
    mutate(chrom = as.character(as.integer(strip_chr_prefix(chrom))),
           variants = normalize_variant_id(id))
  if (is_pvar) {
    LD_variants <- rename(LD_variants, GD = pos)
    LD_variants$GD <- LD_variants$pos <- as.integer(
      sapply(LD_variants$variants, function(v) strsplit(v, ":")[[1]][2]))
  } else {
    LD_variants <- rename(LD_variants, GD = gpos)
  }

  # Label and symmetrize the matrix
  colnames(LD_matrix) <- rownames(LD_matrix) <- LD_variants$variants
  if (all(LD_matrix[lower.tri(LD_matrix)] == 0)) {
    LD_matrix[lower.tri(LD_matrix)] <- t(LD_matrix)[lower.tri(LD_matrix)]
  } else {
    LD_matrix[upper.tri(LD_matrix)] <- t(LD_matrix)[upper.tri(LD_matrix)]
  }

  # Order variants by genomic position
  pos_order <- order(sapply(LD_variants$variants, function(v) as.integer(strsplit(v, ":")[[1]][2])))
  LD_variants <- LD_variants[pos_order, ]
  LD_matrix <- LD_matrix[LD_variants$variants, LD_variants$variants]

  list(LD_matrix = LD_matrix, LD_variants = LD_variants)
}

#' Subset an LD matrix and variant info to a genomic region, optionally
#' further restricted to specific coordinates.
#' @importFrom dplyr mutate select
#' @importFrom magrittr %>%
#' @noRd
extract_LD_for_region <- function(LD_matrix, variants, region, extract_coordinates) {
  extracted <- subset(variants, chrom == region$chrom & pos >= region$start & pos <= region$end)

  if (!is.null(extract_coordinates)) {
    extract_coordinates <- extract_coordinates %>%
      mutate(chrom = as.integer(strip_chr_prefix(chrom))) %>%
      select(chrom, pos)
    extracted <- extracted %>%
      mutate(chrom = as.integer(strip_chr_prefix(chrom))) %>%
      merge(extract_coordinates, by = c("chrom", "pos"))
    keep_cols <- intersect(c("chrom", "variants", "pos", "GD", "A1", "A2",
                             "variance", "allele_freq", "n_nomiss"), names(extracted))
    extracted <- select(extracted, all_of(keep_cols))
  }

  mat <- LD_matrix[extracted$variants, extracted$variants, drop = FALSE]
  list(extracted_LD_matrix = mat, extracted_LD_variants = extracted)
}

#' Combine multiple block-level LD matrices into one, handling boundary overlaps.
#' @importFrom utils tail
#' @noRd
create_LD_matrix <- function(LD_matrices, variants) {
  # Merge variant lists, deduplicating boundary overlaps
  merge_variants <- function(variant_list) {
    merged <- character(0)
    for (v in variant_list) {
      ids <- if (is.list(v) && !is.null(v$variants)) v$variants else v
      if (length(ids) == 0) next
      if (length(merged) > 0 && tail(merged, 1) == ids[1]) ids <- ids[-1]
      merged <- c(merged, ids)
    }
    merged
  }

  all_variants <- merge_variants(variants)
  combined <- matrix(0, nrow = length(all_variants), ncol = length(all_variants),
                     dimnames = list(all_variants, all_variants))

  # Place each block into the combined matrix
  for (i in seq_along(LD_matrices)) {
    v <- rownames(LD_matrices[[i]])
    idx <- match(v, all_variants)
    combined[idx, idx] <- LD_matrices[[i]]
  }
  combined
}

#' Load and Process Linkage Disequilibrium (LD) Matrix
#'
#' Unified entry point for loading LD data from a metadata TSV file.
#'
#' The metadata TSV must have columns: chrom, start, end, path. Two formats:
#' \itemize{
#'   \item Pre-computed LD blocks: many rows per chromosome with block boundaries
#'     in start/end and path pointing to .cor.xz files (optionally comma-separated
#'     with a .bim path).
#'   \item PLINK genotype files: one row per chromosome with start=0, end=0, and
#'     path pointing to a per-chromosome PLINK prefix (.pgen/.pvar[.zst]/.psam or
#'     .bed/.bim/.fam). LD is computed on the fly via \code{compute_LD()}.
#' }
#'
#' @param LD_meta_file_path Path to the LD metadata TSV file.
#' @param region Region of interest: "chr:start-end" string or data.frame with chrom/start/end.
#' @param extract_coordinates Optional data.frame with columns "chrom" and "pos" for
#'   specific coordinates extraction (only for pre-computed LD blocks).
#' @param return_genotype Controls what LD_matrix contains in the return value.
#'   FALSE (default): always return correlation matrix R.
#'   TRUE: return genotype matrix X (only valid for PLINK sources).
#'   "auto": return X for PLINK sources, R for pre-computed sources.
#' @param n_sample Optional sample size for computing variance (= 2*p*(1-p)*n/(n-1)).
#'   If NULL, ref_panel will not include variance or n_nomiss columns.
#'   Only used for PLINK genotype sources.
#'
#' @return A list with:
#' \describe{
#'   \item{LD_variants}{Character vector of variant IDs (canonical format).}
#'   \item{LD_matrix}{LD correlation matrix R (or genotype matrix X when return_genotype is TRUE or "auto" with PLINK source).}
#'   \item{ref_panel}{Data.frame with variant metadata (chrom, pos, A2, A1, variant_id,
#'     and optionally allele_freq, variance, n_nomiss).}
#'   \item{is_genotype}{Logical: TRUE if LD_matrix contains genotype X, FALSE if correlation R.}
#'   \item{block_metadata}{Data.frame with region/block info. For pre-computed LD: one row per block.
#'     For PLINK: a single row spanning the loaded region.}
#' }
#' @export
load_LD_matrix <- function(LD_meta_file_path, region, extract_coordinates = NULL,
                           return_genotype = FALSE, n_sample = NULL) {
  source <- resolve_ld_source(LD_meta_file_path)
  is_geno <- source$type %in% c("plink2", "plink1", "vcf", "gds")

  # "auto": return X for genotype sources, R for pre-computed
  if (identical(return_genotype, "auto")) return_genotype <- is_geno

  if (is_geno) {
    geno_path <- resolve_genotype_path_for_region(source$meta_path, region)
    result <- load_LD_from_genotype(geno_path, region,
                                    return_genotype = return_genotype,
                                    n_sample = n_sample)
  } else {
    # Pre-computed LD blocks (.cor.xz)
    if (return_genotype) {
      stop("return_genotype=TRUE requires genotype files, not pre-computed LD matrices.")
    }
    result <- load_LD_from_blocks(source$meta_path, region, extract_coordinates, n_sample = n_sample)
  }

  # Remove any duplicate variant IDs (safety net for boundary overlaps)
  variant_ids <- getVariantIds(result)
  if (!is.null(variant_ids)) {
    dup_idx <- which(duplicated(variant_ids))
    if (length(dup_idx) > 0) {
      variant_ids_clean <- variant_ids[-dup_idx]
      corr <- getCorrelation(result)
      if (!is.null(corr)) {
        corr <- corr[-dup_idx, -dup_idx, drop = FALSE]
      }
      variants_gr <- result@variants[-dup_idx]
      result <- LDData(
        correlation = corr,
        genotype_handle = result@genotype_handle,
        snp_idx = result@snp_idx,
        variants = variants_gr,
        block_metadata = result@block_metadata,
        n_ref = result@n_ref
      )
    }
  }

  result
}

# ---------- Internal: resolve LD source type ----------

#' @noRd
has_plink2_files <- function(prefix) {
  file.exists(paste0(prefix, ".pgen")) &&
    (file.exists(paste0(prefix, ".pvar")) || file.exists(paste0(prefix, ".pvar.zst"))) &&
    file.exists(paste0(prefix, ".psam"))
}

#' @noRd
has_plink1_files <- function(prefix) {
  file.exists(paste0(prefix, ".bed")) &&
    file.exists(paste0(prefix, ".bim")) &&
    file.exists(paste0(prefix, ".fam"))
}

#' @noRd
is_vcf_path <- function(path) {
  grepl("\\.(vcf|vcf\\.gz|bcf)$", path) && file.exists(path)
}

#' @noRd
is_gds_path <- function(path) {
  grepl("\\.gds$", path) && file.exists(path)
}

#' Check whether a path points to a genotype source (PLINK, VCF, or GDS).
#' @noRd
is_genotype_source <- function(path) {
  has_plink2_files(path) || has_plink1_files(path) || is_vcf_path(path) || is_gds_path(path)
}

#' Resolve an LD source metadata TSV to its actual data type.
#'
#' The metadata TSV has columns: chrom, start, end, path. Three categories are
#' supported:
#' \itemize{
#'   \item Pre-computed LD blocks (.cor.xz): many rows per chromosome, each with
#'     specific start/end block boundaries and path pointing to .cor.xz files.
#'   \item Genotype files (PLINK2, PLINK1, VCF, or GDS): one row per chromosome
#'     with start=0, end=0, and path pointing to a per-chromosome genotype file
#'     or prefix. The actual region filter is applied by the genotype loader.
#' }
#'
#' This function peeks at the first row to determine the data type.
#' The actual per-chromosome path is resolved later by
#' \code{resolve_genotype_path_for_region()} at load time.
#'
#' @param path Path to a metadata TSV file with columns chrom, start, end, path.
#' @return A list with:
#'   \item{type}{"plink2", "plink1", "vcf", "gds", or "precomputed"}
#'   \item{data_path}{Genotype path from first row (for type detection only; actual
#'     per-chromosome path is resolved at load time)}
#'   \item{meta_path}{The metadata TSV path (always set)}
#' @importFrom vroom vroom
#' @noRd
resolve_ld_source <- function(path) {
  if (!file.exists(path)) {
    stop("LD metadata file not found: ", path,
         "\n  Expected: a TSV file with columns chrom, start, end, path.")
  }

  # Peek at first row to determine underlying data type
  meta <- as.data.frame(vroom(path, show_col_types = FALSE, n_max = 1))
  if (ncol(meta) < 4) stop("LD metadata file must have at least 4 columns (chrom, start, end, path): ", path)
  colnames(meta)[1:4] <- c("chrom", "start", "end", "path")
  raw_path <- gsub(",.*$", "", meta$path[1])  # strip comma-separated bim path
  resolved <- file.path(dirname(path), raw_path)

  if (has_plink2_files(resolved)) return(list(type = "plink2", data_path = resolved, meta_path = path))
  if (has_plink1_files(resolved)) return(list(type = "plink1", data_path = resolved, meta_path = path))
  if (is_vcf_path(resolved)) return(list(type = "vcf", data_path = resolved, meta_path = path))
  if (is_gds_path(resolved)) return(list(type = "gds", data_path = resolved, meta_path = path))

  # Pre-computed .cor.xz blocks - verify not using 0:0 sentinel
  if (!is.na(meta$start) && !is.na(meta$end) && meta$start == 0 && meta$end == 0) {
    stop("Metadata has start=0, end=0 but path does not resolve to genotype files: ", resolved,
         "\n  The 0:0 sentinel is only valid for whole-chromosome genotype files.")
  }

  list(type = "precomputed", meta_path = path)
}

#' Resolve the correct genotype path for a given region from a metadata TSV.
#' Reads the TSV, finds the row matching the query region's chromosome,
#' and returns the resolved genotype file path or prefix.
#' @importFrom vroom vroom
#' @noRd
resolve_genotype_path_for_region <- function(meta_path, region) {
  parsed <- parse_region(region)
  meta <- as.data.frame(vroom(meta_path, show_col_types = FALSE))
  colnames(meta) <- c("chrom", "start", "end", "path")
  meta$chrom <- as.integer(strip_chr_prefix(meta$chrom))
  query_chrom <- as.integer(strip_chr_prefix(parsed$chrom))

  matching <- meta[meta$chrom == query_chrom, , drop = FALSE]
  if (nrow(matching) == 0) {
    stop("No entry for chromosome ", query_chrom, " in metadata file: ", meta_path)
  }
  raw_path <- gsub(",.*$", "", matching$path[1])
  file.path(dirname(meta_path), raw_path)
}

# ---------- Internal: load LD from genotype files ----------

#' Load genotype data and compute LD or return genotype matrix.
#' @noRd
load_LD_from_genotype <- function(genotype_path, region,
                                  return_genotype = FALSE, n_sample = NULL) {
  # Load genotype matrix and variant info via the unified loader
  result <- load_genotype_region(genotype_path, region = region,
                                 return_variant_info = TRUE)
  X <- result$X
  variant_info <- result$variant_info

  # Normalize variant IDs to canonical format (chr:pos:A2:A1)
  variant_ids <- normalize_variant_id(
    format_variant_id(variant_info$chrom, variant_info$pos, variant_info$A2, variant_info$A1)
  )
  colnames(X) <- variant_ids

  # Build ref_panel
  ref_panel <- parse_variant_id(variant_ids)
  ref_panel$variant_id <- variant_ids

  # Load allele frequency from .afreq file if available, otherwise compute from genotypes
  afreq <- read_afreq(genotype_path)
  if (!is.null(afreq)) {
    freq_match <- match(variant_info$id, afreq$id)
    n_unmatched <- sum(is.na(freq_match))
    if (n_unmatched > 0) {
      warning(n_unmatched, " out of ", length(freq_match),
              " variants have no allele frequency in .afreq file.")
    }
    ref_panel$allele_freq <- afreq$alt_freq[freq_match]
  } else {
    # Compute ALT allele frequency directly from the dosage matrix
    ref_panel$allele_freq <- colMeans(X, na.rm = TRUE) / 2
  }

  # Compute variance if sample size provided
  if (!is.null(n_sample)) {
    p <- ref_panel$allele_freq
    ref_panel$variance <- 2 * p * (1 - p) * n_sample / (n_sample - 1)
    ref_panel$n_nomiss <- n_sample
  }

  # Block metadata (single block spanning the loaded region)
  positions <- variant_info$pos
  block_metadata <- data.frame(
    block_id = 1L,
    chrom = as.character(variant_info$chrom[1]),
    block_start = min(positions),
    block_end = max(positions),
    size = length(variant_ids),
    start_idx = 1L,
    end_idx = length(variant_ids),
    stringsAsFactors = FALSE
  )

  # Build variant GRanges for LDData
  variants_gr <- .ref_panel_to_granges(ref_panel)

  if (return_genotype) {
    # Store genotype handle + snp_idx for lazy access
    handle <- readGenotypes(genotype_path)
    snp_idx <- .region_to_snp_idx(handle@snp_info, region)
    return(LDData(
      correlation = NULL,
      genotype_handle = handle,
      snp_idx = snp_idx,
      variants = variants_gr,
      block_metadata = block_metadata,
      n_ref = as.integer(nrow(X))
    ))
  }

  R <- compute_LD(X, method = "sample")

  LDData(
    correlation = R,
    genotype_handle = NULL,
    snp_idx = NULL,
    variants = variants_gr,
    block_metadata = block_metadata,
    n_ref = as.integer(nrow(X))
  )
}

# ---------- LD sketch: genotype loading ----------

#' HWE-based standardization of a genotype matrix
#'
#' Centers by 2*allele_freq, scales by sqrt(2*allele_freq*(1-allele_freq)).
#' Assumes monomorphic variants have already been removed.
#'
#' @param X Numeric genotype matrix (n x p).
#' @param allele_freq Numeric vector of allele frequencies (length p).
#' @return Standardized matrix (n x p).
#' @noRd
standardize_genotype_hwe <- function(X, allele_freq) {
  X_std <- sweep(X, 2, 2 * allele_freq)
  sweep(X_std, 2, sqrt(2 * allele_freq * (1 - allele_freq)), "/")
}

#' Load LD sketch genotypes for a region
#'
#' Loads genotype data for a region via \code{load_LD_matrix(return_genotype=TRUE)}
#' and removes monomorphic variants. Returns the raw genotype matrix and metadata,
#' which callers can use to derive either a correlation matrix R (for summary-based
#' weight training or fine-mapping) or an SVD (for TWAS z-score computation).
#'
#' @param ld_meta_file_path Path to the LD metadata TSV file.
#' @param region Region of interest: "chr:start-end" string or data.frame with chrom/start/end.
#' @param n_sample Optional original panel sample size for computing variance
#'   (= 2*p*(1-p)*n/(n-1)). Passed through to \code{load_LD_matrix()}.
#'
#' @return An \code{LDData} S4 object with monomorphic variants removed.
#'   Consumers should use S4 accessors: \code{getGenotypes()}, \code{getRefPanel()},
#'   \code{getVariantIds()}. The number of sketch samples is
#'   \code{nrow(getGenotypes(result))}.
#' @export
load_ld_sketch <- function(ld_meta_file_path, region, n_sample = NULL) {
  result <- load_LD_matrix(ld_meta_file_path, region, return_genotype = TRUE, n_sample = n_sample)
  if (!is(result, "LDData")) {
    stop("load_LD_matrix must return an LDData object")
  }
  X <- getGenotypes(result)
  ref_panel <- getRefPanel(result)

  # Remove monomorphic variants (zero variance under HWE)
  p <- ref_panel$allele_freq
  polymorphic <- p > 0 & p < 1
  if (!all(polymorphic)) {
    X <- X[, polymorphic, drop = FALSE]
    ref_panel <- ref_panel[polymorphic, , drop = FALSE]
  }

  # Rebuild LDData with the extracted (and filtered) genotype matrix stored
  # directly in genotype_handle so getGenotypes() returns it without needing
  # the original file handle.
  variants_gr <- .ref_panel_to_granges(ref_panel)
  LDData(
    correlation = NULL,
    genotype_handle = X,
    snp_idx = NULL,
    variants = variants_gr,
    block_metadata = getBlockMetadata(result),
    n_ref = result@n_ref
  )
}

# ---------- Internal: load LD from pre-computed blocks ----------

#' Load pre-computed LD from block-based metadata files.
#' @noRd
load_LD_from_blocks <- function(LD_meta_file_path, region, extract_coordinates = NULL, n_sample = NULL) {
  # Intersect LD metadata with specified regions
  intersected_LD_files <- get_regional_ld_meta(LD_meta_file_path, region)

  LD_file_paths <- intersected_LD_files$intersections$LD_file_paths
  bim_file_paths <- intersected_LD_files$intersections$bim_file_paths

  extracted_LD_matrices_list <- list()
  extracted_LD_variants_list <- list()
  block_chroms <- character(length(LD_file_paths))

  for (j in seq_along(LD_file_paths)) {
    LD_matrix_processed <- process_LD_matrix(LD_file_paths[j], bim_file_paths[j])
    extracted_LD_list <- extract_LD_for_region(
      LD_matrix = LD_matrix_processed$LD_matrix,
      variants = LD_matrix_processed$LD_variants,
      region = intersected_LD_files$region,
      extract_coordinates = extract_coordinates
    )
    extracted_LD_matrices_list[[j]] <- extracted_LD_list$extracted_LD_matrix
    extracted_LD_variants_list[[j]] <- extracted_LD_list$extracted_LD_variants
    if (nrow(extracted_LD_variants_list[[j]]) > 0) {
      block_chroms[j] <- as.character(extracted_LD_variants_list[[j]]$chrom[1])
    } else {
      block_chroms[j] <- as.character(intersected_LD_files$region$chrom)
    }
  }

  # Filter out empty blocks before combining
  non_empty <- sapply(extracted_LD_variants_list, function(v) nrow(v) > 0)
  if (!any(non_empty)) {
    stop("No variants found in any LD block for the specified region.")
  }
  if (any(!non_empty)) {
    message(paste(
      "Removing", sum(!non_empty), "empty LD block(s) with no variants in the region."
    ))
    extracted_LD_matrices_list <- extracted_LD_matrices_list[non_empty]
    extracted_LD_variants_list <- extracted_LD_variants_list[non_empty]
    block_chroms <- block_chroms[non_empty]
    LD_file_paths <- LD_file_paths[non_empty]
  }

  LD_matrix <- create_LD_matrix(
    LD_matrices = extracted_LD_matrices_list,
    variants = extracted_LD_variants_list
  )
  LD_variants <- rownames(LD_matrix)

  block_variants <- lapply(extracted_LD_variants_list, function(v) v$variants)
  block_positions <- lapply(extracted_LD_variants_list, function(v) v$pos)
  block_metadata <- data.frame(
    block_id = seq_along(LD_file_paths),
    chrom = block_chroms,
    block_start = sapply(block_positions, min),
    block_end = sapply(block_positions, max),
    size = sapply(block_variants, length),
    start_idx = sapply(block_variants, function(v) min(match(v, LD_variants))),
    end_idx = sapply(block_variants, function(v) max(match(v, LD_variants))),
    stringsAsFactors = FALSE
  )

  rm(extracted_LD_matrices_list)

  ref_panel <- parse_variant_id(rownames(LD_matrix))
  merged_variant_list <- do.call(rbind, extracted_LD_variants_list)
  ref_panel$variant_id <- rownames(LD_matrix)

  if ("allele_freq" %in% colnames(merged_variant_list)) {
    ref_panel$allele_freq <- merged_variant_list$allele_freq[match(rownames(LD_matrix), merged_variant_list$variants)]
  }
  if ("variance" %in% colnames(merged_variant_list)) {
    ref_panel$variance <- merged_variant_list$variance[match(rownames(LD_matrix), merged_variant_list$variants)]
  }
  if ("n_nomiss" %in% colnames(merged_variant_list)) {
    ref_panel$n_nomiss <- merged_variant_list$n_nomiss[match(rownames(LD_matrix), merged_variant_list$variants)]
  }

  # Compute variance from n_sample + allele_freq if not already present
  if (!is.null(n_sample) && (!"variance" %in% colnames(ref_panel) || all(is.na(ref_panel$variance)))) {
    if ("allele_freq" %in% colnames(ref_panel)) {
      p <- ref_panel$allele_freq
      ref_panel$variance <- 2 * p * (1 - p) * n_sample / (n_sample - 1)
      ref_panel$n_nomiss <- n_sample
    }
  }

  variants_gr <- .ref_panel_to_granges(ref_panel)

  LDData(
    correlation = LD_matrix,
    genotype_handle = NULL,
    snp_idx = NULL,
    variants = variants_gr,
    block_metadata = block_metadata,
    n_ref = if (is.null(n_sample)) 0L else as.integer(n_sample)
  )
}

#' Filter variants by LD Reference
#'
#' Filters a vector of variant IDs to those present in the LD reference panel.
#' Auto-detects the reference type (PLINK2, PLINK1, or pre-computed LD metadata).
#'
#' @param variant_ids variant names in the format chr:pos:ref:alt.
#' @param ld_reference_meta_file Path to LD metadata file or PLINK prefix.
#' @param keep_indel Whether to keep indel variants. Default TRUE.
#' @return A list with:
#'   \item{data}{Character vector of filtered variant IDs.}
#'   \item{idx}{Integer vector of indices into the original variant_ids.}
#' @importFrom dplyr group_by summarise
#' @importFrom vroom vroom
#' @importFrom magrittr %>%
#' @export
filter_variants_by_ld_reference <- function(variant_ids, ld_reference_meta_file, keep_indel = TRUE) {
  variants_df <- parse_variant_id(variant_ids)

  # Derive region to scope the reference lookup
  region_df <- variants_df %>%
    group_by(chrom) %>%
    summarise(start = min(pos), end = max(pos))

  # Use shared helper -- no genotype loading
  ref_info <- get_ref_variant_info(ld_reference_meta_file, region_df)
  ref_chrom <- as.integer(strip_chr_prefix(ref_info$chrom))
  ref_key <- paste0(ref_chrom, ":", ref_info$pos)

  variant_key <- paste0(variants_df$chrom, ":", variants_df$pos)
  keep_indices <- which(variant_key %in% ref_key)

  if (!keep_indel) {
    snp_idx <- which(is_snp_alleles(variants_df$A1, variants_df$A2))
    keep_indices <- intersect(keep_indices, snp_idx)
  }

  n_dropped <- length(variant_ids) - length(keep_indices)
  if (n_dropped > 0) {
    message(n_dropped, " out of ", length(variant_ids),
            " total variants dropped due to absence on the reference LD panel.")
  }

  list(data = variant_ids[keep_indices], idx = keep_indices)
}

#' Partition LD Matrix into Block-Specific Matrices
#'
#' This function takes the output from load_LD_matrix and partitions the combined LD matrix
#' into a list of smaller matrices based on the block_indices, making it easier to work with
#' large LD matrices that span multiple blocks.
#'
#' @param ld_data An \code{LDData} S4 object as returned by \code{load_LD_matrix()}.
#' @param merge_small_blocks Logical, whether to merge blocks smaller than min_merged_block_size (default: TRUE).
#' @param min_merged_block_size Integer, minimum number of variants for a block after merging (default: 500).
#' @param max_merged_block_size Integer, maximum number of variants in a block after merging (default: 10000).
#'
#' @return returns a list containing:
#' \describe{
#' \item{ld_matrices}{A list of matrices, each representing LD for a specific block.}
#' \item{variant_indices}{A data frame that maps variant IDs to their corresponding block.}
#' \item{block_metadata}{Information about each block including size, chromosome, start and end positions.}
#' }
#' @noRd
partition_LD_matrix <- function(ld_data, merge_small_blocks = TRUE,
                                min_merged_block_size = 500, max_merged_block_size = 10000) {
  if (!is(ld_data, "LDData")) {
    stop("ld_data must be an LDData object")
  }
  combined_matrix <- getCorrelation(ld_data)
  block_metadata <- ld_data@block_metadata
  if (is(block_metadata, "LDBlocks")) {
    block_metadata <- as.data.frame(block_metadata@blocks)
  }
  variant_ids <- getVariantIds(ld_data)

  # Error if matrix is empty
  if (is.null(combined_matrix) || nrow(combined_matrix) == 0 || ncol(combined_matrix) == 0) {
    stop("Empty or NULL LD matrix provided.")
  }

  # Ensure the row and column names of the matrix match the variant_ids
  if (is.null(rownames(combined_matrix)) || is.null(colnames(combined_matrix)) ||
    !identical(rownames(combined_matrix), variant_ids) || !identical(colnames(combined_matrix), variant_ids)) {
    rownames(combined_matrix) <- variant_ids
    colnames(combined_matrix) <- variant_ids
  }

  # Filter out blocks with invalid indices (empty blocks, out-of-range, NA, Inf)
  n_variants <- length(variant_ids)
  valid_blocks <- sapply(seq_len(nrow(block_metadata)), function(i) {
    s <- block_metadata$start_idx[i]
    e <- block_metadata$end_idx[i]
    sz <- block_metadata$size[i]
    # Block is valid if: size > 0, indices are finite integers, and within range
    !is.na(s) && !is.na(e) && is.finite(s) && is.finite(e) &&
      sz > 0 && s >= 1 && e >= s && e <= n_variants
  })

  if (!any(valid_blocks)) {
    stop("No valid LD blocks found. All block indices are out of range or empty.")
  }

  if (any(!valid_blocks)) {
    message(paste(
      "Removing", sum(!valid_blocks),
      "LD block(s) with invalid or out-of-range indices."
    ))
    block_metadata <- block_metadata[valid_blocks, , drop = FALSE]
    block_metadata$block_id <- seq_len(nrow(block_metadata))
  }

  # Validate the block structure of the matrix (skip if only one block)
  if (nrow(block_metadata) > 1) {
    validate_block_structure(combined_matrix, block_metadata, variant_ids)
  }

  # Optionally merge small blocks
  if (merge_small_blocks && any(block_metadata$size < min_merged_block_size) && nrow(block_metadata) > 1) {
    block_metadata <- merge_blocks(block_metadata, min_merged_block_size, max_merged_block_size)
  }

  # Partition the matrix based on block metadata
  result <- extract_block_matrices(combined_matrix, block_metadata, variant_ids)
  return(result)
}

#' Validate that cross-block entries are zero (excluding boundary variants).
#' @noRd
validate_block_structure <- function(matrix, block_metadata, variant_ids) {
  msgs <- character(0)
  n <- length(variant_ids)

  for (i in 1:(nrow(block_metadata) - 1)) {
    for (j in (i + 1):nrow(block_metadata)) {
      si <- block_metadata$start_idx[i]; ei <- block_metadata$end_idx[i]
      sj <- block_metadata$start_idx[j]; ej <- block_metadata$end_idx[j]
      if (si > n || ei > n || sj > n || ej > n) {
        msgs <- c(msgs, paste("Block indices out of range for blocks", i, "and", j))
        next
      }
      # Exclude boundary variants (potential overlaps)
      vi <- variant_ids[si:(ei - 1)]
      vj <- variant_ids[(sj + 1):ej]
      if (length(vi) > 0 && length(vj) > 0) {
        max_val <- max(abs(matrix[vi, vj, drop = FALSE]))
        if (max_val > 1e-10) {
          msgs <- c(msgs, paste("Non-zero correlation between blocks", i, "and", j,
                                "- max:", max_val))
        }
      }
    }
  }
  if (length(msgs) > 0) stop("Matrix lacks expected block structure:\n", paste(msgs, collapse = "\n"))
}

#' @noRd
can_merge <- function(block1, block2, max_size) {
  block1$chrom == block2$chrom && (block1$size + block2$size) <= max_size
}

#' @noRd
merge_two_blocks <- function(block_metadata, idx1, idx2) {
  if (idx1 > idx2) { tmp <- idx1; idx1 <- idx2; idx2 <- tmp }
  result <- block_metadata
  result$end_idx[idx1] <- block_metadata$end_idx[idx2]
  result$size[idx1] <- block_metadata$size[idx1] + block_metadata$size[idx2]
  result <- result[-idx2, ]
  result$block_id <- seq_len(nrow(result))
  result
}

#' Find blocks below min_size and identify the best neighbor to merge with.
#' @noRd
find_merge_candidates <- function(block_metadata, min_size, max_size) {
  candidates <- data.frame(block_idx = integer(), merge_with = integer(), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(block_metadata))) {
    if (block_metadata$size[i] >= min_size) next
    prev_ok <- i > 1 && can_merge(block_metadata[i, ], block_metadata[i - 1, ], max_size)
    next_ok <- i < nrow(block_metadata) && can_merge(block_metadata[i, ], block_metadata[i + 1, ], max_size)
    merge_with <- if (prev_ok && next_ok) {
      if (block_metadata$size[i - 1] <= block_metadata$size[i + 1]) i - 1 else i + 1
    } else if (prev_ok) i - 1
      else if (next_ok) i + 1
      else next
    candidates <- rbind(candidates, data.frame(block_idx = i, merge_with = merge_with))
  }
  candidates
}

#' Iteratively merge blocks below min_size with their smallest neighbor.
#' @noRd
merge_blocks <- function(block_metadata, min_size, max_size) {
  if (nrow(block_metadata) <= 1) return(block_metadata)
  repeat {
    candidates <- find_merge_candidates(block_metadata, min_size, max_size)
    if (nrow(candidates) == 0) break
    block_metadata <- merge_two_blocks(block_metadata, candidates$block_idx[1], candidates$merge_with[1])
  }
  block_metadata
}

# Helper function to extract block matrices
extract_block_matrices <- function(matrix, block_metadata, variant_ids) {
  ld_matrices <- list()
  variant_mapping <- data.frame(
    variant_id = character(),
    block_id = integer(),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(block_metadata))) {
    start_idx <- block_metadata$start_idx[i]
    end_idx <- block_metadata$end_idx[i]

    # Skip empty blocks
    if (end_idx < start_idx) next

    # Ensure indices are within bounds
    if (start_idx > length(variant_ids) || end_idx > length(variant_ids)) {
      warning(paste("Block", i, "has indices outside the range of variant_ids. Skipping."))
      next
    }

    # Extract variant IDs for this block
    block_variants <- variant_ids[start_idx:end_idx]

    # Extract submatrix for this block - use named indexing
    block_matrix <- matrix[block_variants, block_variants, drop = FALSE]

    # Store in list
    ld_matrices[[i]] <- block_matrix

    # Update variant mapping
    block_mapping <- data.frame(
      variant_id = block_variants,
      block_id = i,
      stringsAsFactors = FALSE
    )
    variant_mapping <- rbind(variant_mapping, block_mapping)

  }

  return(list(
    ld_matrices = ld_matrices,
    variant_indices = variant_mapping,
    block_metadata = block_metadata
  ))
}


#' Check and optionally repair LD matrix quality
#'
#' Diagnoses positive-definiteness of an LD correlation matrix and optionally
#' repairs it. Downstream methods like PRS-CS require positive-definite LD
#' (Cholesky decomposition), while others (lassosum, SDPR) handle non-PD
#' matrices internally via their own regularization.
#'
#' Three modes are available:
#' \describe{
#'   \item{\code{"check"}}{Diagnostic only - returns eigenvalue statistics
#'     without modifying the matrix.}
#'   \item{\code{"shrink"}}{Apply shrinkage toward identity:
#'     \code{R_s = (1 - shrinkage) * R + shrinkage * I}. Simple and fast;
#'     always produces a positive-definite matrix when \code{shrinkage > 0}.}
#'   \item{\code{"eigenfix"}}{Set negative eigenvalues to zero and
#'     reconstruct the matrix. Matches the approach used in susieR's
#'     \code{rss_lambda_constructor} and is the closest positive
#'     semidefinite matrix in the Frobenius norm. Does not inflate the
#'     diagonal like shrinkage does.}
#' }
#'
#' @param R Symmetric correlation matrix.
#' @param method One of \code{"check"}, \code{"shrink"}, or \code{"eigenfix"}.
#' @param r_tol Eigenvalue tolerance. Eigenvalues with absolute value below
#'   \code{r_tol} are treated as zero. Default: \code{1e-8}.
#' @param shrinkage Shrinkage parameter for \code{method = "shrink"}.
#'   Default: \code{0.01}.
#'
#' @return A list with components:
#' \describe{
#'   \item{R}{The (possibly repaired) LD matrix.}
#'   \item{is_pd}{Logical: is the matrix positive definite?}
#'   \item{is_psd}{Logical: is the matrix positive semidefinite (within r_tol)?}
#'   \item{min_eigenvalue}{Smallest eigenvalue of the original matrix.}
#'   \item{n_negative}{Number of negative eigenvalues (below -r_tol).}
#'   \item{condition_number}{Ratio of largest to smallest positive eigenvalue
#'     (\code{Inf} if any eigenvalue is zero).}
#'   \item{method_applied}{Character: \code{"none"}, \code{"shrink"}, or
#'     \code{"eigenfix"}.}
#' }
#'
#' @examples
#' # A well-conditioned matrix
#' R_good <- diag(5)
#' check_ld(R_good)$is_pd  # TRUE
#'
#' # A matrix with negative eigenvalues
#' R_bad <- matrix(0.9, 3, 3); diag(R_bad) <- 1; R_bad[1,3] <- R_bad[3,1] <- -0.5
#' check_ld(R_bad)$is_psd  # FALSE
#' R_fixed <- check_ld(R_bad, method = "eigenfix")$R
#' check_ld(R_fixed)$is_psd  # TRUE
#'
#' @export
check_ld <- function(R,
                     method = c("check", "shrink", "eigenfix"),
                     r_tol = 1e-8,
                     shrinkage = 0.01) {
  method <- match.arg(method)
  p <- nrow(R)

  # Eigen decomposition (symmetric)
  eig <- eigen(R, symmetric = TRUE)
  vals <- eig$values

  # Diagnostics
  min_eval <- min(vals)
  n_neg <- sum(vals < -r_tol)
  pos_vals <- vals[vals > r_tol]
  cond_num <- if (length(pos_vals) > 0) max(pos_vals) / min(pos_vals) else Inf
  is_psd <- !any(vals < -r_tol)
  is_pd <- all(vals > r_tol)

  method_applied <- "none"
  R_out <- R

  if (method == "shrink" && !is_pd) {
    R_out <- (1 - shrinkage) * R + shrinkage * diag(p)
    method_applied <- "shrink"
  } else if (method == "eigenfix" && !is_pd) {
    # Set negative eigenvalues to a small positive value and reconstruct.
    # Using r_tol (not zero) ensures the result is strictly positive
    # definite, which is required by methods that use Cholesky decomposition
    # (PRS-CS, SDPR). Setting to exactly zero would produce PSD but not PD.
    vals_fixed <- pmax(vals, r_tol)
    R_out <- eig$vectors %*% diag(vals_fixed) %*% t(eig$vectors)
    # Restore exact symmetry and unit diagonal
    R_out <- (R_out + t(R_out)) / 2
    diag(R_out) <- 1
    method_applied <- "eigenfix"
  }

  list(
    R = R_out,
    is_pd = is_pd,
    is_psd = is_psd,
    min_eigenvalue = min_eval,
    n_negative = n_neg,
    condition_number = cond_num,
    method_applied = method_applied
  )
}

#' Prune columns by pairwise correlation (LD-style prune)
#'
#' Performs LD pruning using one of two backends. The default \code{"hclust"}
#' backend computes the full correlation matrix, builds a single-linkage
#' hierarchical clustering on the distance (1 - |cor|), and keeps one
#' representative column per cluster. The \code{"snprelate"} backend delegates
#' to \code{SNPRelate::snpgdsLDpruning}, which performs a sliding-window
#' greedy prune directly on a temporary GDS file.
#'
#' @param X Numeric matrix. Columns are the variables to prune (typically SNP
#'   genotype dosages); rows are observations.
#' @param cor_thres Numeric in (0, 1). Absolute correlation threshold.
#'   Columns whose pairwise |cor| exceeds this are grouped; one survivor is
#'   kept per group. Default 0.8.
#' @param backend Character, one of \code{"hclust"} (default) or
#'   \code{"snprelate"}. Controls the pruning algorithm:
#'   \describe{
#'     \item{\code{"hclust"}}{Uses the internal hierarchical-clustering approach
#'       with \code{Rfast::cora} (if available) or base \code{cor()}.}
#'     \item{\code{"snprelate"}}{Requires \pkg{SNPRelate} and \pkg{gdsfmt}.
#'       Creates a temporary GDS file and runs
#'       \code{SNPRelate::snpgdsLDpruning(method = "corr")}.}
#'   }
#' @param verbose Logical. If TRUE, print progress messages. Default FALSE.
#'
#' @return A list with:
#'   \describe{
#'     \item{X.new}{Matrix containing the retained columns of \code{X}.}
#'     \item{filter.id}{Integer vector of the column indices of \code{X} that
#'       were retained (in original order).}
#'   }
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(100 * 5), 100, 5)
#' X[, 2] <- X[, 1] + rnorm(100, sd = 0.01)   # near-duplicate of col 1
#' res <- ld_prune_by_correlation(X, cor_thres = 0.9)
#' ncol(res$X.new)
#'
#' @importFrom stats as.dist hclust cutree cor
#' @export
ld_prune_by_correlation <- function(X, cor_thres = 0.8,
                                    backend = c("hclust", "snprelate"),
                                    verbose = FALSE) {
  backend <- match.arg(backend)
  p <- ncol(X)

  if (backend == "snprelate") {
    return(.ld_prune_snprelate(X, cor_thres = cor_thres, verbose = verbose))
  }

  # ---- hclust backend (default) ----
  if (requireNamespace("Rfast", quietly = TRUE)) {
    cor.X <- Rfast::cora(X, large = TRUE)
  } else {
    cor.X <- cor(X)
  }
  Sigma.distance <- as.dist(1 - abs(cor.X))
  fit <- hclust(Sigma.distance, method = "single")
  clusters <- cutree(fit, h = 1 - cor_thres)
  groups <- unique(clusters)
  ind.delete <- NULL
  X.new <- X
  filter.id <- seq_len(p)
  for (ig in seq_along(groups)) {
    temp.group <- which(clusters == groups[ig])
    if (length(temp.group) > 1) {
      ind.delete <- c(ind.delete, temp.group[-1])
    }
  }
  ind.delete <- unique(ind.delete)
  if (length(ind.delete) > 0) {
    X.new <- as.matrix(X[, -ind.delete])
    filter.id <- filter.id[-ind.delete]
    if (verbose) {
      message("ld_prune_by_correlation: pruned ", length(ind.delete),
              " of ", p, " columns at |cor| > ", cor_thres)
    }
  } else if (verbose) {
    message("ld_prune_by_correlation: no columns pruned at |cor| > ", cor_thres)
  }

  if (ncol(X.new) == 1) {
    colnames(X.new) <- colnames(X)[-ind.delete]
  }

  list(X.new = X.new, filter.id = filter.id)
}

#' SNPRelate-based LD pruning helper
#' @noRd
.ld_prune_snprelate <- function(X, cor_thres, verbose) {
  if (!requireNamespace("SNPRelate", quietly = TRUE) ||
      !requireNamespace("gdsfmt", quietly = TRUE)) {
    stop("Packages 'SNPRelate' and 'gdsfmt' are required for backend='snprelate'.")
  }
  p <- ncol(X)
  snp_names <- colnames(X) %||% paste0("snp", seq_len(p))

  # Round dosages to integer genotype codes for GDS
  geno_int <- round(X)
  storage.mode(geno_int) <- "integer"

  tmp_gds <- tempfile(fileext = ".gds")
  on.exit(unlink(tmp_gds), add = TRUE)

  SNPRelate::snpgdsCreateGeno(
    gds.fn = tmp_gds,
    genmat = t(geno_int),
    sample.id = seq_len(nrow(X)),
    snp.id = seq_len(p),
    snp.rs.id = snp_names,
    snp.chromosome = rep(1L, p),
    snp.position = seq_len(p),
    snpfirstdim = TRUE
  )

  gds <- SNPRelate::snpgdsOpen(tmp_gds, allow.duplicate = TRUE)
  on.exit(SNPRelate::snpgdsClose(gds), add = TRUE)

  keep_list <- SNPRelate::snpgdsLDpruning(
    gds,
    method = "corr",
    ld.threshold = cor_thres,
    verbose = verbose
  )

  keep_ids <- sort(unlist(keep_list, use.names = FALSE))
  filter.id <- keep_ids
  X.new <- X[, keep_ids, drop = FALSE]

  if (verbose) {
    message("ld_prune_by_correlation (snprelate): kept ", length(keep_ids),
            " of ", p, " columns at |cor| > ", cor_thres)
  }

  list(X.new = X.new, filter.id = filter.id)
}

#' Drop collinear columns from a design matrix by a chosen strategy
#'
#' Given a numeric matrix \code{X} and a set of column names known to be
#' involved in linear dependencies, remove one column using one of three
#' strategies. Designed to be called iteratively by
#' \code{\link{enforce_design_full_rank}}, but can be used standalone.
#'
#' @param X Numeric matrix. Must have column names covering
#'   \code{problematic_cols}.
#' @param problematic_cols Character vector of column names in \code{X} that
#'   are candidates for removal. If empty, \code{X} is returned unchanged.
#' @param strategy One of \code{"correlation"} (remove the column with the
#'   largest sum of absolute pairwise correlations among the candidates;
#'   when only two candidates, one is picked at random), \code{"variance"}
#'   (remove the lowest-variance candidate), or \code{"response_correlation"}
#'   (remove the candidate whose correlation with \code{response} has the
#'   smallest magnitude).
#' @param response Numeric vector required when \code{strategy =
#'   "response_correlation"}; the outcome to correlate against.
#' @param verbose Logical. If TRUE, print which column was removed. Default
#'   FALSE.
#'
#' @return \code{X} with exactly one column removed (or unchanged if
#'   \code{problematic_cols} is empty).
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(100 * 3), 100, 3)
#' X[, 3] <- X[, 1] + X[, 2]
#' colnames(X) <- c("a", "b", "c")
#' drop_collinear_columns(X, problematic_cols = c("a", "b", "c"),
#'                        strategy = "variance")
#'
#' @importFrom stats var cor
#' @keywords internal
#' @noRd
drop_collinear_columns <- function(X, problematic_cols,
                                   strategy = c("correlation", "variance", "response_correlation"),
                                   response = NULL, verbose = FALSE) {
  strategy <- match.arg(strategy)

  if (length(problematic_cols) == 0) {
    return(X)
  }

  if (length(problematic_cols) == 1) {
    col_to_remove <- problematic_cols[1]
    if (verbose) message("drop_collinear_columns: removing single column ", col_to_remove)
    X <- X[, !(colnames(X) %in% col_to_remove), drop = FALSE]
    return(X)
  }

  if (strategy == "variance") {
    variances <- apply(X[, problematic_cols, drop = FALSE], 2, var)
    col_to_remove <- problematic_cols[which.min(variances)]
    if (verbose) message("drop_collinear_columns: smallest variance -> removing ", col_to_remove)
  } else if (strategy == "correlation") {
    cor_matrix <- abs(cor(X[, problematic_cols, drop = FALSE]))
    diag(cor_matrix) <- 0

    if (length(problematic_cols) == 2) {
      col_to_remove <- sample(problematic_cols, 1)
      if (verbose) message("drop_collinear_columns: two candidates, randomly removing ", col_to_remove)
    } else {
      cor_sums <- colSums(cor_matrix)
      col_to_remove <- problematic_cols[which.max(cor_sums)]
      if (verbose) message("drop_collinear_columns: highest sum |cor| -> removing ", col_to_remove)
    }
  } else if (strategy == "response_correlation") {
    if (is.null(response)) {
      stop("response must be supplied for strategy = 'response_correlation'")
    }
    cor_with_response <- apply(X[, problematic_cols, drop = FALSE], 2,
                               function(col) cor(col, response))
    col_to_remove <- problematic_cols[which.min(abs(cor_with_response))]
    if (verbose) message("drop_collinear_columns: smallest |cor| with response -> removing ", col_to_remove)
  }

  X[, !(colnames(X) %in% col_to_remove), drop = FALSE]
}

#' Iteratively enforce full column rank on a design matrix
#'
#' Given a candidate predictor matrix \code{X} and an optional unnamed
#' covariate matrix \code{C}, builds the design \code{[1, X, C]} and removes
#' rank-deficient columns from \code{X} until the design has full column rank.
#' Rank-deficient columns are identified via the pivot of
#' \code{qr([1, X, C])}. On each iteration, one problematic column is dropped
#' using \code{\link{drop_collinear_columns}}. If iterative pruning does not
#' achieve full rank, falls back to \code{\link{ld_prune_by_correlation}} at a
#' descending sequence of correlation thresholds.
#'
#' @param X Numeric matrix with column names (the predictors subject to
#'   pruning).
#' @param C Numeric matrix of covariates (can be unnamed) that will be kept.
#'   Pass \code{NULL} or a zero-column matrix when there are no covariates.
#' @param strategy Passed through to \code{\link{drop_collinear_columns}}.
#' @param response Passed through to \code{\link{drop_collinear_columns}}
#'   when \code{strategy = "response_correlation"}.
#' @param max_iterations Integer. Hard cap on the iterative-prune loop.
#'   Default 300.
#' @param corr_thresholds Numeric vector of |cor| thresholds used for the
#'   \code{\link{ld_prune_by_correlation}} fallback, tried in order.
#'   Default \code{seq(0.75, 0.5, by = -0.05)}.
#' @param verbose Logical. If TRUE, print per-iteration progress. Default
#'   FALSE.
#'
#' @return The pruned predictor matrix \code{X} (covariates \code{C} are not
#'   modified).
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(100 * 4), 100, 4)
#' X[, 4] <- X[, 1] + X[, 2]          # rank-deficient
#' colnames(X) <- c("a", "b", "c", "d")
#' C <- matrix(rnorm(100), 100, 1)
#' X2 <- enforce_design_full_rank(X, C, strategy = "variance")
#' qr(cbind(1, X2, C))$rank == ncol(cbind(1, X2, C))
#'
#' @export
enforce_design_full_rank <- function(X, C,
                                     strategy = c("correlation", "variance", "response_correlation"),
                                     response = NULL,
                                     max_iterations = 300L,
                                     corr_thresholds = seq(0.75, 0.5, by = -0.05),
                                     verbose = FALSE) {
  strategy <- match.arg(strategy)
  original_colnames <- colnames(X)
  initial_ncol <- ncol(X)
  iteration <- 0L

  build_design <- function(X) {
    XD <- cbind(1, X, C)
    colnames(XD)[seq_len(ncol(X) + 1L)] <- c("Intercept", colnames(X))
    XD
  }

  X_design <- build_design(X)
  matrix_rank <- qr(X_design)$rank
  if (verbose) {
    message("enforce_design_full_rank: initial rank ", matrix_rank,
            " / ", ncol(X_design))
  }

  skip_iterative <- FALSE

  # Fast path: try removing all QR-pivot-flagged columns at once.
  if (matrix_rank < ncol(X_design)) {
    qrd <- qr(X_design)
    problematic_cols <- qrd$pivot[(qrd$rank + 1L):ncol(X_design)]
    problematic_colnames <- colnames(X_design)[problematic_cols]
    problematic_colnames <- problematic_colnames[problematic_colnames %in% colnames(X)]

    if (length(problematic_colnames) > 0) {
      X_temp <- X[, !(colnames(X) %in% problematic_colnames), drop = FALSE]
      if (qr(build_design(X_temp))$rank == ncol(build_design(X_temp))) {
        if (verbose) {
          message("enforce_design_full_rank: full rank after batch-removing ",
                  length(problematic_colnames), " column(s)")
        }
      } else {
        skip_iterative <- TRUE
        if (verbose) {
          message("enforce_design_full_rank: batch removal insufficient, ",
                  "skipping to correlation-pruning fallback")
        }
      }
    }
  }

  # Iterative path.
  if (!skip_iterative) {
    while (matrix_rank < ncol(X_design) && iteration < max_iterations) {
      qrd <- qr(X_design)
      problematic_cols <- qrd$pivot[(qrd$rank + 1L):ncol(X_design)]
      problematic_colnames <- colnames(X_design)[problematic_cols]
      problematic_colnames <- problematic_colnames[problematic_colnames %in% colnames(X)]

      if (length(problematic_colnames) == 0) break

      X <- drop_collinear_columns(X, problematic_colnames, strategy = strategy,
                                  response = response, verbose = verbose)

      X_design <- build_design(X)
      matrix_rank <- qr(X_design)$rank
      iteration <- iteration + 1L
      if (verbose) {
        message("enforce_design_full_rank: iter ", iteration,
                " rank ", matrix_rank, " / ", ncol(X_design))
      }
    }

    if (iteration == max_iterations) {
      warning("enforce_design_full_rank: max_iterations reached; design may still be rank-deficient")
    }
  }

  # Correlation-threshold fallback.
  X_design <- build_design(X)
  matrix_rank <- qr(X_design)$rank
  if (matrix_rank < ncol(X_design)) {
    if (verbose) {
      message("enforce_design_full_rank: applying ld_prune_by_correlation fallback")
    }
    for (threshold in corr_thresholds) {
      filter_result <- ld_prune_by_correlation(X, cor_thres = threshold,
                                               verbose = verbose)
      X <- filter_result$X.new
      X_design <- build_design(X)
      matrix_rank <- qr(X_design)$rank
      if (verbose) {
        message("enforce_design_full_rank: threshold ", threshold,
                " -> rank ", matrix_rank, " / ", ncol(X_design))
      }
      if (matrix_rank == ncol(X_design)) break
    }
  }

  if (ncol(X) == 1L && initial_ncol == 1L) {
    colnames(X) <- original_colnames
  }
  X
}

#' LD clumping by a per-variant score using bigsnpr
#'
#' Wraps \code{bigsnpr::snp_clumping} with the boilerplate of wrapping a
#' numeric dosage matrix into a \code{bigstatsr::FBM.code256} object and of
#' handling the common pitfall of a single-variant input.
#'
#' @param X Numeric matrix of 0/1/2 allele dosages, n rows by p variants.
#'   Column names are expected to be variant IDs but are not required.
#' @param score Numeric vector of length \code{ncol(X)}. Higher values favour
#'   retention during clumping (e.g. -log10 p, |Z|, MAF). May be \code{NULL},
#'   in which case bigsnpr falls back to minor allele frequency computed from
#'   \code{X}.
#' @param chr Integer or character vector of length \code{ncol(X)} giving the
#'   chromosome for each variant.
#' @param pos Integer vector of length \code{ncol(X)} giving the base-pair
#'   position for each variant.
#' @param r2 Numeric in (0, 1]. r-squared threshold for clumping (variants
#'   within \code{window_kb} whose r2 exceeds \code{r2} and have lower
#'   \code{score} are removed). Default 0.2.
#' @param window_kb Numeric. Window size in kilobases. Default is
#'   \code{100 / r2}, matching the common "ld-clump size = 100/r2" heuristic
#'   used in many GWAS pipelines.
#' @param verbose Logical. If TRUE, print the number of retained variants.
#'   Default FALSE.
#'
#' @return An integer vector of indices (into \code{X} columns) kept after
#'   clumping. For a single-column \code{X}, returns \code{1L}.
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   n <- 500; p <- 20
#'   X <- matrix(rbinom(n * p, 2, 0.3), n, p)
#'   colnames(X) <- paste0("chr1:", seq_len(p) * 1000, ":A:G")
#'   s <- runif(p)
#'   chr <- rep(1L, p); pos <- seq_len(p) * 1000L
#'   keep <- ld_clump_by_score(X, score = s, chr = chr, pos = pos, r2 = 0.2)
#' }
#'
#' @export
ld_clump_by_score <- function(X, score, chr, pos, r2 = 0.2,
                              window_kb = 100 / r2, verbose = FALSE) {
  if (!requireNamespace("bigsnpr", quietly = TRUE)) {
    stop("Package 'bigsnpr' is required. Install from CRAN: install.packages('bigsnpr')")
  }
  if (!requireNamespace("bigstatsr", quietly = TRUE)) {
    stop("Package 'bigstatsr' is required. Install from CRAN: install.packages('bigstatsr')")
  }

  if (ncol(X) < 1L) stop("ld_clump_by_score: X must have at least one column")
  if (!is.null(score) && length(score) != ncol(X)) {
    stop("ld_clump_by_score: length(score) must equal ncol(X)")
  }
  if (length(chr) != ncol(X) || length(pos) != ncol(X)) {
    stop("ld_clump_by_score: chr and pos must have length equal to ncol(X)")
  }

  if (ncol(X) == 1L) {
    if (verbose) message("ld_clump_by_score: single variant, skipping clumping")
    return(1L)
  }

  if (inherits(X, "FBM")) {
    G <- X
  } else {
    code_vec <- c(0, 1, 2, rep(NA, 256L - 3L))
    G <- bigstatsr::FBM.code256(
      nrow = nrow(X), ncol = ncol(X),
      init = X, code = code_vec
    )
  }

  keep <- bigsnpr::snp_clumping(
    G = G,
    infos.chr = as.integer(chr),
    infos.pos = as.integer(pos),
    S = score,
    thr.r2 = r2,
    size = window_kb
  )

  if (verbose) {
    message("ld_clump_by_score: ", length(keep), " / ", ncol(X),
            " variants retained at r2 <= ", r2)
  }
  keep
}
