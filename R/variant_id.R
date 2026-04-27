#' @title Variant ID Parsing and Formatting Utilities
#' @description Functions for parsing, formatting, normalizing, and detecting
#'   the naming conventions of variant IDs (e.g., "chr1:100:A:G") and genomic
#'   region strings (e.g., "chr1:100-200").
#' @importFrom stringr str_split
NULL

#' Strip "chr" prefix from chromosome identifiers.
#' @param x Character vector of chromosome identifiers (e.g., "chr1", "chrX").
#' @return Character vector with "chr" prefix removed (e.g., "1", "X").
#' @noRd
strip_chr_prefix <- function(x) sub("^chr", "", x)

#' Strip build suffix from variant IDs (e.g., ":b38" or "_b38").
#' @param x Character vector of variant IDs.
#' @return Character vector with build suffix removed.
#' @noRd
strip_build_suffix <- function(x) sub("(:|_)b[0-9]+$", "", x)

#' Test whether allele pairs are single-nucleotide (SNP, not indel).
#'
#' Returns TRUE for each pair where both alleles are exactly one of A, T, C, G.
#' @param a1 Character vector of first alleles.
#' @param a2 Character vector of second alleles.
#' @return Logical vector, TRUE if the variant is a SNP.
#' @noRd
is_snp_alleles <- function(a1, a2) {
  nchar(a1) == 1L & nchar(a2) == 1L &
    grepl("^[ATCG]$", a1) & grepl("^[ATCG]$", a2)
}

#' Detect the naming convention of variant IDs
#'
#' Examines variant ID strings to detect their format: whether they have a "chr"
#' prefix, what separator is used between allele fields, and whether they include
#' a genome build suffix (e.g., ":b38" or "_b38").
#'
#' Supported formats:
#' \itemize{
#'   \item All colons: \code{"chr1:100:A:G"} or \code{"1:100:A:G"}
#'   \item Mixed colon/underscore: \code{"chr1:100_A_G"} or \code{"1:100_A_G"}
#'   \item All underscores (PLINK BIM): \code{"chr1_100_A_G"} or \code{"1_100_A_G"}
#' }
#'
#' @param ids A character vector of variant IDs.
#' @return A list with components:
#'   \describe{
#'     \item{has_chr}{Logical, whether the IDs have a "chr" prefix.}
#'     \item{allele_sep}{Character, the separator between allele fields (":" or "_").
#'       For mixed format \code{"chr1:100_A_G"}, this is \code{"_"}.}
#'     \item{has_build}{Logical, whether a build suffix is present.}
#'     \item{example}{Character, the first non-NA ID for reference.}
#'   }
#' @noRd
detect_variant_convention <- function(ids) {
  # Find first non-NA element
  first_id <- ids[!is.na(ids)][1]
  if (is.na(first_id) || length(first_id) == 0) {
    return(list(has_chr = FALSE, allele_sep = ":", has_build = FALSE, example = NA_character_))
  }
  has_chr <- grepl("^chr", first_id)
  # Detect build suffix like :b38 or _b38 at end
  has_build <- grepl("(:|_)b[0-9]+$", first_id)
  id_clean <- strip_build_suffix(first_id)
  # Detect allele separator: check if variant uses underscores between allele fields
  # This catches both full underscore ("1_100_A_G") and mixed ("chr1:100_A_G") formats
  allele_sep <- if (grepl("_[ATCGID*]+_[ATCGID*]+$", id_clean)) "_" else ":"
  list(has_chr = has_chr, allele_sep = allele_sep, has_build = has_build, example = first_id)
}

#' Parse variant IDs into a data frame
#'
#' Converts variant IDs from any supported string format or data.frame into a
#' standardized data.frame with integer chrom, integer pos, and character allele
#' columns (A2, A1). Supports colon-separated ("chr1:100:A:G"), underscore-separated
#' ("1_100_A_G"), with or without "chr" prefix, and with optional build suffix
#' (":b38" or "_b38"). The detected input convention is stored as an attribute.
#'
#' @param ids A character vector of variant IDs, or a data.frame with columns
#'   "chrom", "pos", and allele columns (A2/A1 or ref/alt or any 4-column layout).
#' @return A data.frame with columns "chrom" (integer), "pos" (integer), "A2"
#'   (character), "A1" (character). The detected convention is stored as
#'   \code{attr(result, "convention")}.
#' @export
parse_variant_id <- function(ids) {
  # Handle data.frame input
  if (is.data.frame(ids)) {
    if (all(c("chrom", "pos", "A2", "A1") %in% names(ids))) {
      # Already has correct column names
    } else if (all(c("chrom", "pos", "A1", "A2") %in% names(ids))) {
      # Has A1/A2 but need to check they're in the right semantic order
      # (A2 = ref, A1 = alt/effect) -- keep as-is since column names are explicit
    } else if (ncol(ids) >= 4) {
      # Assume positional: chrom, pos, A2, A1
      names(ids)[1:4] <- c("chrom", "pos", "A2", "A1")
    }
    # Detect convention from chrom column before converting
    conv <- list(
      has_chr = any(grepl("^chr", as.character(ids$chrom))),
      allele_sep = ":", has_build = FALSE, example = NA_character_
    )
    ids$chrom <- as.integer(strip_chr_prefix(as.character(ids$chrom)))
    ids$pos <- as.integer(ids$pos)
    attr(ids, "convention") <- conv
    return(ids)
  }

  # Detect convention before parsing
  convention <- detect_variant_convention(ids)

  # Normalize: convert underscores to colons, strip build suffix
  normalized <- gsub("_", ":", ids)
  normalized <- strip_build_suffix(normalized)

  # Split into exactly 4 fields using strcapture (vectorized, no list overhead)
  data <- strcapture(
    "^([^:]+):([^:]+):([^:]+):([^:]+)",
    normalized,
    proto = data.frame(chrom = character(), pos = character(),
                       A2 = character(), A1 = character(),
                       stringsAsFactors = FALSE)
  )

  data$chrom <- as.integer(strip_chr_prefix(data$chrom))
  data$pos <- as.integer(data$pos)

  attr(data, "convention") <- convention
  return(data)
}

#' Format variant ID strings from component columns
#'
#' Constructs variant ID strings from chrom, pos, A2, A1 columns.
#' The chrom:pos separator is always a colon. The allele separator can be
#' either colon (canonical: \code{"chr1:100:A:G"}) or underscore
#' (mixed: \code{"chr1:100_A_G"}).
#'
#' When a \code{convention} object (from \code{detect_variant_convention}) is
#' provided, the output format is driven automatically by the detected
#' convention, so callers do not need to specify \code{chr_prefix} or
#' \code{allele_sep} manually.
#'
#' @param chrom Integer or character chromosome (e.g., 1 or "chr1").
#' @param pos Integer position.
#' @param A2 Character reference allele.
#' @param A1 Character alternate/effect allele.
#' @param chr_prefix Logical, whether to add "chr" prefix. Default TRUE.
#'   Ignored if \code{convention} is provided.
#' @param allele_sep Character, separator between pos/A2 and A2/A1 fields.
#'   Default \code{":"} produces canonical \code{"chr1:100:A:G"};
#'   \code{"_"} produces mixed \code{"chr1:100_A_G"}.
#'   Ignored if \code{convention} is provided.
#' @param convention Optional list from \code{detect_variant_convention}.
#'   When provided, \code{has_chr} and \code{allele_sep} are read from the
#'   convention automatically. This is the preferred way to preserve the
#'   user's input format.
#' @return A character vector of formatted variant IDs.
#' @noRd
format_variant_id <- function(chrom, pos, A2, A1, chr_prefix = TRUE, allele_sep = ":", convention = NULL) {
  # If convention is provided, use it to determine format automatically
  if (!is.null(convention)) {
    chr_prefix <- convention$has_chr
    allele_sep <- if (!is.null(convention$allele_sep)) convention$allele_sep else ":"
  }
  # Strip any existing chr prefix to normalize, then re-add if requested
  chrom_clean <- strip_chr_prefix(as.character(chrom))
  if (chr_prefix) {
    paste0("chr", chrom_clean, ":", pos, allele_sep, A2, allele_sep, A1)
  } else {
    paste0(chrom_clean, ":", pos, allele_sep, A2, allele_sep, A1)
  }
}

#' Normalize variant IDs to canonical format
#'
#' One-step convenience function: parses variant IDs in any supported format
#' and re-formats them. By default, outputs the canonical format
#' (\code{"chr{N}:{pos}:{A2}:{A1}"}). When a \code{convention} object is
#' provided, the output preserves the user's original format automatically.
#'
#' @param ids A character vector of variant IDs in any supported format.
#' @param chr_prefix Logical, whether to include "chr" prefix. Default TRUE.
#'   Ignored if \code{convention} is provided.
#' @param convention Optional list from \code{detect_variant_convention} or
#'   \code{attr(parse_variant_id(ids), "convention")}. When provided, the
#'   output format is driven automatically by the detected convention.
#' @return A character vector of normalized variant IDs.
#' @export
normalize_variant_id <- function(ids, chr_prefix = TRUE, convention = NULL) {
  parsed <- parse_variant_id(ids)
  if (!is.null(convention)) {
    format_variant_id(parsed$chrom, parsed$pos, parsed$A2, parsed$A1, convention = convention)
  } else {
    format_variant_id(parsed$chrom, parsed$pos, parsed$A2, parsed$A1, chr_prefix = chr_prefix)
  }
}

# Internal convenience wrapper around parse_variant_id.
variant_id_to_df <- function(variant_id) {
  parse_variant_id(variant_id)
}

#' @importFrom stringr str_split
#' @export
parse_region <- function(region) {
  if (!is.character(region) || length(region) != 1) {
    return(region)
  }

  if (!grepl("^chr[0-9XY]+:[0-9]+-[0-9]+$", region)) {
    stop("Input string format must be 'chr:start-end'.")
  }
  parts <- str_split(region, "[:-]")[[1]]
  df <- data.frame(
    chrom = strip_chr_prefix(parts[1]),
    start = as.integer(parts[2]),
    end = as.integer(parts[3])
  )

  return(df)
}

#' Utility function to convert LD region_ids to `region of interest` dataframe
#' @param ld_region_id A string of region in the format of chrom_start_end.
#' @export
region_to_df <- function(ld_region_id, colnames = c("chrom", "start", "end")) {
  region_of_interest <- as.data.frame(do.call(rbind, lapply(strsplit(ld_region_id, "[_:-]"), function(x) as.integer(strip_chr_prefix(x)))))
  colnames(region_of_interest) <- colnames
  return(region_of_interest)
}

#' Ensure two sets of variant IDs use matching chr prefix convention
#'
#' Detects whether \code{ids_a} and \code{ids_b} have mismatched chr prefixes.
#' If mismatched, normalizes both to canonical format (with "chr" prefix) using
#' \code{\link{normalize_variant_id}}. If already matching, returns inputs
#' unchanged.
#'
#' @param ids_a Character vector of variant IDs.
#' @param ids_b Character vector of variant IDs.
#' @return A list with components \code{ids_a} and \code{ids_b}, both normalized
#'   to canonical chr-prefix format if they were mismatched.
#' @noRd
ensure_chr_match <- function(ids_a, ids_b) {
  has_chr_a <- any(grepl("^chr", ids_a[!is.na(ids_a)][1:min(5, sum(!is.na(ids_a)))]))
  has_chr_b <- any(grepl("^chr", ids_b[!is.na(ids_b)][1:min(5, sum(!is.na(ids_b)))]))
  if (has_chr_a == has_chr_b) {
    return(list(ids_a = ids_a, ids_b = ids_b))
  }
  list(
    ids_a = normalize_variant_id(ids_a, chr_prefix = TRUE),
    ids_b = normalize_variant_id(ids_b, chr_prefix = TRUE)
  )
}

#' Convert region specifications to a GRanges object
#'
#' Accepts region strings ("chr1:100-200", "1_100_200"), character vectors of
#' such strings, or data.frames with chrom/start/end columns. Returns a
#' \code{\link[GenomicRanges]{GRanges}} object.
#'
#' @param regions A region string, character vector, or data.frame with
#'   chrom/start/end columns.
#' @return A \code{GRanges} object.
#' @noRd
as_granges <- function(regions) {
  if (is.character(regions)) {
    df <- region_to_df(regions)
  } else if (is.data.frame(regions)) {
    if (!all(c("chrom", "start", "end") %in% names(regions))) {
      stop("data.frame must have columns: chrom, start, end")
    }
    df <- regions
  } else {
    stop("regions must be a character vector or data.frame with chrom/start/end columns")
  }
  # GRanges expects character seqnames; prefix with "chr" if numeric
  seqnames <- as.character(df$chrom)
  if (!any(grepl("^chr", seqnames))) {
    seqnames <- paste0("chr", seqnames)
  }
  GenomicRanges::GRanges(
    seqnames = seqnames,
    ranges = IRanges::IRanges(start = as.integer(df$start), end = as.integer(df$end))
  )
}

#' Test whether two genomic regions overlap
#'
#' @param region_a A region string ("chr1:100-200" or "1_100_200") or a
#'   single-row data.frame with chrom/start/end columns.
#' @param region_b A region string or single-row data.frame.
#' @return Logical scalar: TRUE if the regions share at least one base pair.
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges findOverlaps
#' @export
regions_overlap <- function(region_a, region_b) {
  gr_a <- as_granges(region_a)
  gr_b <- as_granges(region_b)
  length(IRanges::findOverlaps(gr_a, gr_b)) > 0
}

#' Find which target regions overlap a query region
#'
#' @param query A single region string or single-row data.frame with
#'   chrom/start/end columns.
#' @param targets A character vector of region strings, or a multi-row
#'   data.frame with chrom/start/end columns.
#' @return Integer vector of 1-based indices into \code{targets} that overlap
#'   the query. Empty integer vector if no overlaps.
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges findOverlaps
#' @importFrom S4Vectors subjectHits
#' @export
find_overlapping_regions <- function(query, targets) {
  gr_query <- as_granges(query)
  gr_targets <- as_granges(targets)
  hits <- IRanges::findOverlaps(gr_query, gr_targets)
  unique(S4Vectors::subjectHits(hits))
}

#' Classify variant type from allele strings
#'
#' Determines whether each variant is a SNP, insertion, deletion, or
#' multi-nucleotide polymorphism (MNP) based on the allele lengths.
#'
#' @param ids A character vector of variant IDs in "chr:pos:ref:alt" format,
#'   or a data.frame with A2 (ref) and A1 (alt) columns (e.g., from
#'   \code{\link{parse_variant_id}}).
#' @return A character vector with one of "SNP", "insertion", "deletion", or
#'   "MNP" for each variant.
#' @export
classify_variant_type <- function(ids) {
  if (is.character(ids)) {
    ids <- parse_variant_id(ids)
  }
  if (!is.data.frame(ids) || !all(c("A2", "A1") %in% names(ids))) {
    stop("Input must be a character vector of variant IDs or a data.frame with A2 and A1 columns.")
  }
  len_ref <- nchar(ids$A2)
  len_alt <- nchar(ids$A1)
  type <- character(nrow(ids))
  type[len_ref == 1L & len_alt == 1L & grepl("^[ATCG]$", ids$A2) & grepl("^[ATCG]$", ids$A1)] <- "SNP"
  type[len_ref == len_alt & (len_ref > 1L | !grepl("^[ATCG]$", ids$A2) | !grepl("^[ATCG]$", ids$A1)) & type == ""] <- "MNP"
  type[len_ref > len_alt] <- "deletion"
  type[len_alt > len_ref] <- "insertion"
  type
}
