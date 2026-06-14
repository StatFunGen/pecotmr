#' @title Variant ID Parsing and Formatting Utilities
#' @description Functions for parsing, formatting, normalizing, and detecting
#'   the naming conventions of variant IDs (e.g., "chr1:100:A:G") and genomic
#'   region strings (e.g., "chr1:100-200").
#' @name pecotmr-variant-id
#' @keywords internal
#' @importFrom stringr str_split
NULL

#' Strip "chr" prefix from chromosome identifiers.
#' @param x Character vector of chromosome identifiers (e.g., "chr1", "chrX").
#' @return Character vector with "chr" prefix removed (e.g., "1", "X").
#' @noRd
stripChrPrefix <- function(x) sub("^chr", "", x)

# Backwards-compat alias

#' Strip build suffix from variant IDs (e.g., ":b38" or "_b38").
#' @param x Character vector of variant IDs.
#' @return Character vector with build suffix removed.
#' @noRd
stripBuildSuffix <- function(x) sub("(:|_)b[0-9]+$", "", x)

# Backwards-compat alias

#' Test whether allele pairs are single-nucleotide (SNP, not indel).
#'
#' Returns TRUE for each pair where both alleles are exactly one of A, T, C, G.
#' @param a1 Character vector of first alleles.
#' @param a2 Character vector of second alleles.
#' @return Logical vector, TRUE if the variant is a SNP.
#' @noRd
isSnpAlleles <- function(a1, a2) {
  nchar(a1) == 1L & nchar(a2) == 1L &
    grepl("^[ATCG]$", a1) & grepl("^[ATCG]$", a2)
}

# Backwards-compat alias

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
#'     \item{hasChr}{Logical, whether the IDs have a "chr" prefix.}
#'     \item{alleleSep}{Character, the separator between allele fields (":" or "_").
#'       For mixed format \code{"chr1:100_A_G"}, this is \code{"_"}.}
#'     \item{hasBuild}{Logical, whether a build suffix is present.}
#'     \item{example}{Character, the first non-NA ID for reference.}
#'   }
#' @noRd
detectVariantConvention <- function(ids) {
  # Find first non-NA element
  firstId <- ids[!is.na(ids)][1]
  if (is.na(firstId) || length(firstId) == 0) {
    return(list(hasChr = FALSE, alleleSep = ":", hasBuild = FALSE, example = NA_character_))
  }
  hasChr <- grepl("^chr", firstId)
  # Detect build suffix like :b38 or _b38 at end
  hasBuild <- grepl("(:|_)b[0-9]+$", firstId)
  idClean <- stripBuildSuffix(firstId)
  # Detect allele separator: check if variant uses underscores between allele fields
  # This catches both full underscore ("1_100_A_G") and mixed ("chr1:100_A_G") formats
  alleleSep <- if (grepl("_[ATCGID*]+_[ATCGID*]+$", idClean)) "_" else ":"
  list(hasChr = hasChr, alleleSep = alleleSep, hasBuild = hasBuild, example = firstId)
}

# Backwards-compat alias for external callers

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
parseVariantId <- function(ids) {
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
      hasChr = any(grepl("^chr", as.character(ids$chrom))),
      alleleSep = ":", hasBuild = FALSE, example = NA_character_
    )
    ids$chrom <- as.integer(stripChrPrefix(as.character(ids$chrom)))
    ids$pos <- as.integer(ids$pos)
    attr(ids, "convention") <- conv
    return(ids)
  }

  # Detect convention before parsing
  convention <- detectVariantConvention(ids)

  # Normalize: convert underscores to colons, strip build suffix
  normalized <- gsub("_", ":", ids)
  normalized <- stripBuildSuffix(normalized)

  # Split into exactly 4 fields using strcapture (vectorized, no list overhead)
  data <- strcapture(
    "^([^:]+):([^:]+):([^:]+):([^:]+)",
    normalized,
    proto = data.frame(chrom = character(), pos = character(),
                       A2 = character(), A1 = character(),
                       stringsAsFactors = FALSE)
  )

  data$chrom <- as.integer(stripChrPrefix(data$chrom))
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
#' \code{alleleSep} manually.
#'
#' @param chrom Integer or character chromosome (e.g., 1 or "chr1").
#' @param pos Integer position.
#' @param A2 Character reference allele.
#' @param A1 Character alternate/effect allele.
#' @param chrPrefix Logical, whether to add "chr" prefix. Default TRUE.
#'   Ignored if \code{convention} is provided.
#' @param alleleSep Character, separator between pos/A2 and A2/A1 fields.
#'   Default \code{":"} produces canonical \code{"chr1:100:A:G"};
#'   \code{"_"} produces mixed \code{"chr1:100_A_G"}.
#'   Ignored if \code{convention} is provided.
#' @param convention Optional list from \code{detectVariantConvention}.
#'   When provided, \code{hasChr} and \code{alleleSep} are read from the
#'   convention automatically. This is the preferred way to preserve the
#'   user's input format.
#' @return A character vector of formatted variant IDs.
#' @noRd
formatVariantId <- function(chrom, pos, A2, A1, chrPrefix = TRUE, alleleSep = ":", convention = NULL) {
  # If convention is provided, use it to determine format automatically
  if (!is.null(convention)) {
    chrPrefix <- convention$hasChr
    alleleSep <- if (!is.null(convention$alleleSep)) convention$alleleSep else ":"
  }
  # Strip any existing chr prefix to normalize, then re-add if requested
  chromClean <- stripChrPrefix(as.character(chrom))
  if (chrPrefix) {
    paste0("chr", chromClean, ":", pos, alleleSep, A2, alleleSep, A1)
  } else {
    paste0(chromClean, ":", pos, alleleSep, A2, alleleSep, A1)
  }
}

# Backwards-compat alias for external callers

#' Normalize variant IDs to canonical format
#'
#' One-step convenience function: parses variant IDs in any supported format
#' and re-formats them. By default, outputs the canonical format
#' (\code{"chr{N}:{pos}:{A2}:{A1}"}). When a \code{convention} object is
#' provided, the output preserves the user's original format automatically.
#'
#' @param ids A character vector of variant IDs in any supported format.
#' @param chrPrefix Logical, whether to include "chr" prefix. Default TRUE.
#'   Ignored if \code{convention} is provided.
#' @param convention Optional list from \code{detectVariantConvention} or
#'   \code{attr(parseVariantId(ids), "convention")}. When provided, the
#'   output format is driven automatically by the detected convention.
#' @return A character vector of normalized variant IDs.
#' @export
normalizeVariantId <- function(ids, chrPrefix = TRUE, convention = NULL) {
  parsed <- parseVariantId(ids)
  if (!is.null(convention)) {
    formatVariantId(parsed$chrom, parsed$pos, parsed$A2, parsed$A1, convention = convention)
  } else {
    formatVariantId(parsed$chrom, parsed$pos, parsed$A2, parsed$A1, chrPrefix = chrPrefix)
  }
}

#' @export

# Internal convenience wrapper around parseVariantId.
variantIdToDf <- function(variantId) {
  parseVariantId(variantId)
}

# Backwards-compat alias for external callers

#' @importFrom stringr str_split
#' @export
parseRegion <- function(region) {
  if (!is.character(region) || length(region) != 1) {
    return(region)
  }

  if (!grepl("^chr[0-9XY]+:[0-9]+-[0-9]+$", region)) {
    stop("Input string format must be 'chr:start-end'.")
  }
  parts <- str_split(region, "[:-]")[[1]]
  df <- data.frame(
    chrom = stripChrPrefix(parts[1]),
    start = as.integer(parts[2]),
    end = as.integer(parts[3])
  )

  return(df)
}

#' Utility function to convert LD region_ids to `region of interest` dataframe
#' @param ldRegionId A string of region in the format of chrom_start_end.
#' @export
regionToDf <- function(ldRegionId, colnames = c("chrom", "start", "end")) {
  regionOfInterest <- as.data.frame(do.call(rbind, lapply(strsplit(ldRegionId, "[_:-]"), function(x) as.integer(stripChrPrefix(x)))))
  colnames(regionOfInterest) <- colnames
  return(regionOfInterest)
}

#' Ensure two sets of variant IDs use matching chr prefix convention
#'
#' Detects whether \code{idsA} and \code{idsB} have mismatched chr prefixes.
#' If mismatched, normalizes both to canonical format (with "chr" prefix) using
#' \code{\link{normalizeVariantId}}. If already matching, returns inputs
#' unchanged.
#'
#' @param idsA Character vector of variant IDs.
#' @param idsB Character vector of variant IDs.
#' @return A list with components \code{idsA} and \code{idsB}, both normalized
#'   to canonical chr-prefix format if they were mismatched.
#' @noRd
ensureChrMatch <- function(idsA, idsB) {
  hasChrA <- any(grepl("^chr", idsA[!is.na(idsA)][1:min(5, sum(!is.na(idsA)))]))
  hasChrB <- any(grepl("^chr", idsB[!is.na(idsB)][1:min(5, sum(!is.na(idsB)))]))
  if (hasChrA == hasChrB) {
    return(list(idsA = idsA, idsB = idsB))
  }
  list(
    idsA = normalizeVariantId(idsA, chrPrefix = TRUE),
    idsB = normalizeVariantId(idsB, chrPrefix = TRUE)
  )
}

# Backwards-compat alias for external callers

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
asGranges <- function(regions) {
  if (is.character(regions)) {
    df <- regionToDf(regions)
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

# Backwards-compat alias for external callers

#' Test whether two genomic regions overlap
#'
#' @param regionA A region string ("chr1:100-200" or "1_100_200") or a
#'   single-row data.frame with chrom/start/end columns.
#' @param regionB A region string or single-row data.frame.
#' @return Logical scalar: TRUE if the regions share at least one base pair.
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges findOverlaps
#' @export
regionsOverlap <- function(regionA, regionB) {
  grA <- asGranges(regionA)
  grB <- asGranges(regionB)
  length(IRanges::findOverlaps(grA, grB)) > 0
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
findOverlappingRegions <- function(query, targets) {
  grQuery <- asGranges(query)
  grTargets <- asGranges(targets)
  hits <- IRanges::findOverlaps(grQuery, grTargets)
  unique(S4Vectors::subjectHits(hits))
}

#' Classify variant type from allele strings
#'
#' Determines whether each variant is a SNP, insertion, deletion, or
#' multi-nucleotide polymorphism (MNP) based on the allele lengths.
#'
#' @param ids A character vector of variant IDs in "chr:pos:ref:alt" format,
#'   or a data.frame with A2 (ref) and A1 (alt) columns (e.g., from
#'   \code{\link{parseVariantId}}).
#' @return A character vector with one of "SNP", "insertion", "deletion", or
#'   "MNP" for each variant.
#' @export
classifyVariantType <- function(ids) {
  if (is.character(ids)) {
    ids <- parseVariantId(ids)
  }
  if (!is.data.frame(ids) || !all(c("A2", "A1") %in% names(ids))) {
    stop("Input must be a character vector of variant IDs or a data.frame with A2 and A1 columns.")
  }
  lenRef <- nchar(ids$A2)
  lenAlt <- nchar(ids$A1)
  type <- character(nrow(ids))
  type[lenRef == 1L & lenAlt == 1L & grepl("^[ATCG]$", ids$A2) & grepl("^[ATCG]$", ids$A1)] <- "SNP"
  type[lenRef == lenAlt & (lenRef > 1L | !grepl("^[ATCG]$", ids$A2) | !grepl("^[ATCG]$", ids$A1)) & type == ""] <- "MNP"
  type[lenRef > lenAlt] <- "deletion"
  type[lenAlt > lenRef] <- "insertion"
  type
}

