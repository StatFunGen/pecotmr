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

#' Canonicalize chromosome identifiers to a normalized string.
#'
#' Strips a leading "chr"/"ch" prefix (case-insensitive), uppercases, and maps
#' common synonyms to a single form (23 -> X, 24 -> Y, M -> MT). Returns a
#' character vector and never coerces to integer, so non-autosomal chromosomes
#' (X, Y, MT) survive instead of collapsing to NA. This is the single chromosome
#' normalizer used for variant matching and identity throughout the package;
#' it mirrors the chromosome handling previously inlined in summary-statistics
#' QC.
#'
#' @param x Character or numeric vector of chromosome identifiers.
#' @return Character vector of normalized chromosome names.
#' @noRd
canonChrom <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x <- sub("^ch", "", x, ignore.case = TRUE)
  x <- toupper(x)
  ok <- !is.na(x)
  x[ok & x == "23"] <- "X"
  x[ok & x == "24"] <- "Y"
  x[ok & x == "M"] <- "MT"
  x
}

#' Ensure a leading chr prefix on a chromosome identifier.
#'
#' Case-insensitive additive companion to \code{canonChrom} (which strips the
#' prefix): \code{"1"}, \code{"chr1"}, \code{"CHR1"} all become \code{"chr1"}.
#' Does not remap 23/24/M.
#' @param x Character or numeric vector of chromosome identifiers.
#' @return Character vector with a lowercase \code{"chr"} prefix.
#' @noRd
withChrPrefix <- function(x) {
  paste0("chr", sub("^chr", "", as.character(x), ignore.case = TRUE))
}

#' Order key for chromosome identifiers.
#'
#' Returns an ordered factor placing autosomes 1..22 first, then X, Y, XY, MT,
#' and finally any non-standard contigs (in order of appearance). Use this in
#' place of \code{as.integer(chrom)} wherever chromosomes must be sorted, so
#' that X / Y / MT sort sensibly instead of collapsing to NA.
#'
#' @param x Vector of chromosome identifiers.
#' @return An ordered factor suitable for \code{order()} / \code{arrange()}.
#' @noRd
chromOrder <- function(x) {
  x <- canonChrom(x)
  standard <- c(as.character(1:22), "X", "Y", "XY", "MT")
  extra <- setdiff(unique(x[!is.na(x)]), standard)
  factor(x, levels = c(standard, extra), ordered = TRUE)
}

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
#' standardized data.frame with a normalized character chrom, integer pos, and
#' character allele columns (A2, A1). The chrom is canonicalized (chr stripped,
#' uppercased, 23 -> X / 24 -> Y / M -> MT) but kept as a string so X/Y/MT
#' survive. Supports colon-separated ("chr1:100:A:G"), underscore-separated
#' ("1_100_A_G"), with or without "chr" prefix, and with optional build suffix
#' (":b38" or "_b38"). The detected input convention is stored as an attribute.
#'
#' @param ids A character vector of variant IDs, or a data.frame with columns
#'   "chrom", "pos", and allele columns (A2/A1 or ref/alt or any 4-column layout).
#' @return A data.frame with columns "chrom" (character, normalized), "pos"
#'   (integer), "A2" (character), "A1" (character). The detected convention is
#'   stored as
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
    ids$chrom <- canonChrom(ids$chrom)
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

  data$chrom <- canonChrom(data$chrom)
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
  # Normalize the chromosome (strip prefix, uppercase, map synonyms) then re-add
  # the "chr" prefix if requested. canonChrom keeps X/Y/MT as strings.
  chromClean <- canonChrom(chrom)
  if (chrPrefix) {
    paste0("chr", chromClean, ":", pos, alleleSep, A2, alleleSep, A1)
  } else {
    paste0(chromClean, ":", pos, alleleSep, A2, alleleSep, A1)
  }
}

# Backwards-compat alias for external callers

#' Re-format variant IDs to a chosen output convention
#'
#' Output-only convenience: parses variant IDs and re-emits them in a single
#' format. By default produces the canonical format
#' (\code{"chr{N}:{pos}:{A2}:{A1}"}); pass a \code{convention} to preserve the
#' input format instead.
#'
#' This normalizes only the textual \emph{format} (chr prefix and field
#' separators); it does NOT reorder alleles, so it is not a matching/identity
#' operation -- two records that differ by a ref/alt swap remain distinct
#' strings. For identity, match on the (chrom, pos, ref, alt) tuple via the
#' variant matcher; use this only for display, file output, or relabeling for
#' name-based downstream consumers.
#'
#' @param ids A character vector of variant IDs in any supported format.
#' @param chrPrefix Logical, whether to include "chr" prefix. Default TRUE.
#'   Ignored if \code{convention} is provided.
#' @param convention Optional list from \code{detectVariantConvention} or
#'   \code{attr(parseVariantId(ids), "convention")}. When provided, the
#'   output format is driven automatically by the detected convention.
#' @return A character vector of re-formatted variant IDs. Unparseable ids
#'   (e.g. rsIDs) are returned unchanged.
#' @export
normalizeVariantId <- function(ids, chrPrefix = TRUE, convention = NULL) {
  parsed <- parseVariantId(ids)
  out <- as.character(ids)
  # Only re-format ids that parsed into a chrom + pos; leave unparseable ids
  # (e.g. rsIDs) unchanged rather than emitting "chrNA:..." garbage.
  ok <- !is.na(parsed$chrom) & !is.na(parsed$pos)
  if (any(ok)) {
    out[ok] <- if (!is.null(convention))
      formatVariantId(parsed$chrom[ok], parsed$pos[ok], parsed$A2[ok],
                      parsed$A1[ok], convention = convention)
    else
      formatVariantId(parsed$chrom[ok], parsed$pos[ok], parsed$A2[ok],
                      parsed$A1[ok], chrPrefix = chrPrefix)
  }
  out
}

#' @export

# Internal convenience wrapper around parseVariantId.
variantIdToDf <- function(variantId) {
  parseVariantId(variantId)
}

#' Harmonize variant alleles against a reference
#'
#' The allele-harmonization engine for the package (used by summary-statistics
#' QC, ctwasPipeline, and -- via \code{matchVariants} -- the pipeline join
#' sites). Matches a target against a reference by ("chrom", "pos") and the
#' allele pair ("A2", "A1"), accounting for strand flips and major/minor allele
#' (ref/alt) swaps, and sign-flips the specified columns when alleles are
#' swapped relative to the reference.
#'
#' @param targetData A data frame with columns "chrom", "pos", "A2", "A1" (and optionally other columns like "beta" or "z"),
#'   or a vector of strings in the format of "chr:pos:A2:A1"/"chr:pos_A2_A1". Can be automatically converted to a data frame if a vector.
#' @param refVariants A data frame with columns "chrom", "pos", "A2", "A1" or strings in the format of "chr:pos:A2:A1"/"chr:pos_A2_A1".
#' @param colToFlip The name of the column in targetData where flips are to be applied.
#'   On an allele swap these columns are sign-flipped (multiplied by -1), the
#'   correct operation for signed quantities like \code{beta} and \code{z}.
#' @param colToComplement Names of columns in targetData to complement
#'   (\code{1 - x}) on an allele swap, the correct operation for an
#'   effect-allele frequency like \code{af}. Default \code{character()} does no
#'   complementing, so non-RSS callers are unchanged. Distinct from
#'   \code{colToFlip}: frequencies are complemented, signed effects are
#'   sign-flipped.
#' @param matchMinProp Minimum proportion of variants in the smallest data
#'   to be matched, otherwise stops with an error. Default is 20%.
#' @param removeDups Whether to remove duplicates, default is TRUE.
#' @param removeIndels Whether to remove INDELs, default is FALSE.
#' @param flip Whether the alleles must be flipped: A <--> T & C <--> G, in which case
#'   corresponding `colToFlip` are multiplied by -1. Default is `TRUE`.
#' @param removeStrandAmbiguous Whether to remove strand SNPs (if any). Default is `TRUE`.
#' @param flipStrand Whether to output the variants after strand flip. Default is `FALSE`.
#' @param removeUnmatched Whether to remove unmatched variants. Default is `TRUE`.
#' @return An \code{AlleleQcResult} S4 object. Use
#'   \code{$harmonizedData} to recover the post-QC variant
#'   data.frame and \code{$qcSummary} to inspect the per-variant
#'   merge/flip/strand diagnostics.
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate inner_join filter pull select everything row_number if_else any_of all_of rename
#' @importFrom vctrs vec_duplicate_detect
#' @importFrom tidyr separate
#' @keywords internal
#' @noRd
#' @details
#' Pure panel-vs-sumstats allele harmonization: match by (chrom, pos),
#' detect A1/A2 swap, sign-flip \code{colToFlip} columns and complement
#' \code{colToComplement} columns on swap. Variant-allele filters
#' (indels, strand-ambiguous, duplicates) are applied here directly when
#' the corresponding \code{removeIndels} / \code{removeStrandAmbiguous} /
#' \code{removeDups} flags are set; MAF / INFO / N column-numeric filters
#' run in \code{.applyContentFilters()} before this function.
harmonizeAlleles <- function(targetData, refVariants, colToFlip = NULL,
                           matchMinProp = 0.2, flipStrand = FALSE,
                           removeUnmatched = TRUE,
                           removeIndels = FALSE,
                           removeStrandAmbiguous = TRUE,
                           removeDups = FALSE,
                           colToComplement = character(), ...) {
  strandFlip <- function(ref) chartr("ATCG", "TAGC", ref)

  sanitizeNames <- function(df) {
    nm <- colnames(df)
    if (is.null(nm)) nm <- rep("unnamed", ncol(df))
    emptyIdx <- is.na(nm) | nm == ""
    if (any(emptyIdx))
      nm[emptyIdx] <- paste0("unnamed_", seq_len(sum(emptyIdx)))
    colnames(df) <- make.unique(nm, sep = "_")
    df
  }

  if (is.data.frame(targetData)) {
    if (ncol(targetData) > 4 &&
        all(c("chrom", "pos", "A2", "A1") %in% names(targetData))) {
      variantCols <- c("chrom", "pos", "A2", "A1")
      variantDf <- targetData %>% select(all_of(variantCols))
      otherCols <- targetData %>% select(-all_of(variantCols))
      targetData <- cbind(variantIdToDf(variantDf), otherCols)
    } else {
      targetData <- variantIdToDf(targetData)
    }
  } else {
    targetData <- variantIdToDf(targetData)
  }
  refVariants <- variantIdToDf(refVariants)

  # Strip merge-conflicting columns; keep target A1/A2. `variant_id` is also
  # stripped from `refVariants` because the post-harmonization variant_id is
  # rebuilt from the QC'd alleles further down (`variants_id_qced`), and
  # leaving the input variant_id on either side causes the final rename to
  # collide on duplicate names.
  columnsToRemove <- c("chromosome", "position", "ref", "alt", "variant_id")
  if (any(columnsToRemove %in% colnames(targetData)))
    targetData <- select(targetData, -any_of(columnsToRemove))
  if ("variant_id" %in% colnames(refVariants))
    refVariants <- select(refVariants, -any_of("variant_id"))

  matchResult <- inner_join(targetData, refVariants,
                            by = c("chrom", "pos"),
                            suffix = c(".target", ".ref")) %>%
                 as.data.frame() %>%
                 sanitizeNames()

  if (nrow(matchResult) == 0) {
    warning("No matching variants found between target data and reference variants.")
    emptyOut <- list(harmonizedData = matchResult, qcSummary = matchResult)
    attr(emptyOut, "qcCounts") <- list(
      considered = 0L, signFlip = 0L, strandFlip = 0L, kept = 0L,
      dropped = 0L, droppedIndel = 0L, droppedAmbiguous = 0L,
      droppedOther = 0L)
    return(emptyOut)
  }

  matchResult <- matchResult %>%
    mutate(variants_id_original = formatVariantId(chrom, pos, A2.target, A1.target),
           variants_id_qced     = formatVariantId(chrom, pos, A2.ref, A1.ref)) %>%
    mutate(across(c(A1.target, A2.target, A1.ref, A2.ref), toupper)) %>%
    mutate(flip1.ref = strandFlip(A1.ref),
           flip2.ref = strandFlip(A2.ref)) %>%
    # AT / CG pairs cannot be distinguished from strand-flip without external
    # context; the keep rule below relies on this flag as a safety guard for
    # callers that may not have removed strand-ambiguous variants upstream.
    mutate(strand_unambiguous = if_else(
      (A1.target == "A" & A2.target == "T") |
      (A1.target == "T" & A2.target == "A") |
      (A1.target == "C" & A2.target == "G") |
      (A1.target == "G" & A2.target == "C"),
      FALSE, TRUE)) %>%
    mutate(exact_match = A1.target == A1.ref & A2.target == A2.ref) %>%
    mutate(sign_flip   = ((A1.target == A2.ref & A2.target == A1.ref) |
                         (A1.target == flip2.ref & A2.target == flip1.ref)) &
                        (A1.target != A1.ref & A2.target != A2.ref)) %>%
    mutate(strand_flip = ((A1.target == flip1.ref & A2.target == flip2.ref) |
                         (A1.target == flip2.ref & A2.target == flip1.ref)) &
                        (A1.target != A1.ref & A2.target != A2.ref)) %>%
    # INDEL detection: explicit "I"/"D" notation, or any allele wider than 1bp.
    mutate(INDEL = (A2.target == "I" | A2.target == "D" |
                   nchar(A2.target) > 1L | nchar(A1.target) > 1L)) %>%
    # ID_match: an indel encoded as I/D on the target side matches an indel
    # on the reference side (where the reference uses multi-base alleles).
    mutate(ID_match = ((A2.target == "D" | A2.target == "I") &
                      (nchar(A1.ref) > 1L | nchar(A2.ref) > 1L)))

  # When removeStrandAmbiguous = FALSE, the A/T - C/G safety guard is
  # disabled: ambiguous variants are treated as exact/sign-flip cases.
  if (!removeStrandAmbiguous)
    matchResult$strand_unambiguous <- TRUE

  # If no strand_flip survives the unambiguous test, the remaining ambiguous
  # variants can be treated as exact/sign-flip cases rather than dropped.
  if (!any(matchResult$strand_flip & matchResult$strand_unambiguous))
    matchResult$strand_unambiguous <- TRUE

  matchResult <- matchResult %>%
    mutate(keep = if_else(strand_flip,
                          true  = strand_unambiguous | exact_match | ID_match,
                          false = exact_match | sign_flip | ID_match))

  if (removeIndels)
    matchResult <- matchResult %>%
      mutate(keep = if_else(INDEL, FALSE, keep))

  if (!is.null(colToFlip)) {
    missing <- setdiff(colToFlip, colnames(matchResult))
    if (length(missing) > 0L)
      stop("Column(s) '", paste(missing, collapse = "', '"),
           "' not found in targetData.")
    matchResult[matchResult$sign_flip, colToFlip] <-
      -1 * matchResult[matchResult$sign_flip, colToFlip]
  }
  # A frequency tracks the effect allele, so an allele swap takes af -> 1 - af
  # (not a sign flip). Kept independent of colToFlip so signed columns are
  # untouched here.
  if (length(colToComplement) > 0L) {
    missing <- setdiff(colToComplement, colnames(matchResult))
    if (length(missing) > 0L)
      stop("Column(s) '", paste(missing, collapse = "', '"),
           "' not found in targetData.")
    matchResult[matchResult$sign_flip, colToComplement] <-
      1 - matchResult[matchResult$sign_flip, colToComplement]
  }
  if (flipStrand) {
    sIdx <- which(matchResult$strand_flip)
    matchResult[sIdx, "A1.target"] <- strandFlip(matchResult[sIdx, "A1.target"])
    matchResult[sIdx, "A2.target"] <- strandFlip(matchResult[sIdx, "A2.target"])
  }

  # Per-step QC counts (used by .runEntrySummaryStatsQc for "kept N of M
  # (corrected: sign-flipped A, strand-flipped B; dropped C)" logging).
  # Computed from the per-variant flags before they are stripped from the
  # returned data frame so callers reading the data frame are unaffected.
  qcCounts <- list(
    considered = nrow(matchResult),
    signFlip   = sum(matchResult$sign_flip   & matchResult$keep, na.rm = TRUE),
    strandFlip = sum(matchResult$strand_flip & matchResult$keep, na.rm = TRUE),
    kept       = sum(matchResult$keep,  na.rm = TRUE),
    dropped    = sum(!matchResult$keep, na.rm = TRUE))
  if ("INDEL" %in% colnames(matchResult)) {
    qcCounts$droppedIndel <- sum(!matchResult$keep & matchResult$INDEL,
                                  na.rm = TRUE)
  } else {
    qcCounts$droppedIndel <- 0L
  }
  qcCounts$droppedAmbiguous <- sum(
    !matchResult$keep & matchResult$strand_flip &
    !matchResult$strand_unambiguous &
    if ("INDEL" %in% colnames(matchResult)) !matchResult$INDEL else TRUE,
    na.rm = TRUE)
  qcCounts$droppedOther <- qcCounts$dropped - qcCounts$droppedIndel -
                           qcCounts$droppedAmbiguous

  result <- matchResult[matchResult$keep, , drop = FALSE]

  qcCols <- c("flip1.ref", "flip2.ref", "strand_unambiguous",
              "exact_match", "sign_flip", "strand_flip", "INDEL",
              "ID_match", "keep")
  result <- result %>%
    select(-any_of(qcCols), -A1.target, -A2.target) %>%
    rename(A1 = A1.ref, A2 = A2.ref, variant_id = variants_id_qced)

  # removeDups: drop duplicate variant rows (same chrom/pos/qced ID).
  # Default FALSE keeps the existing strict behavior (error on dups).
  if (removeDups) {
    dups <- duplicated(result[, c("chrom", "pos", "variant_id")])
    if (any(dups)) {
      warning(sprintf("Removed %d duplicate variant(s), keeping first occurrence.",
                      sum(dups)))
      result <- result[!dups, , drop = FALSE]
    }
  }

  if (!removeUnmatched) {
    matchVariant <- result %>% pull(variants_id_original)
    matchResult <- matchResult %>%
      select(-any_of(qcCols), -variants_id_original, -A1.target, -A2.target) %>%
      rename(A1 = A1.ref, A2 = A2.ref, variant_id = variants_id_qced)
    targetData <- targetData %>%
      mutate(variant_id = formatVariantId(chrom, pos, A2, A1))
    if (length(setdiff(targetData %>% pull(variant_id), matchVariant)) > 0L) {
      unmatchData <- targetData %>% filter(!variant_id %in% matchVariant)
      result <- rbind(result,
                      unmatchData %>% mutate(variants_id_original = variant_id))
      result <- result[match(targetData$variant_id,
                             result$variants_id_original), ] %>%
                select(-variants_id_original)
    }
  }

  if (nrow(result) < matchMinProp * nrow(refVariants))
    stop("Not enough variants have been matched.")
  if (any(duplicated(result$variant_id)))
    stop("Duplicated variant IDs remain after harmonization; pass ",
         "removeDups = TRUE or deduplicate upstream before calling ",
         "harmonizeAlleles.")

  out <- list(harmonizedData = result, qcSummary = matchResult)
  attr(out, "qcCounts") <- qcCounts
  out
}

#' Allele-aware variant matcher (match by chrom/pos/ref/alt, not id string)
#'
#' The single matching primitive for the package: given two sets of variant
#' identifiers, match them on the (chrom, pos) position and the allele pair
#' (exact, ref/alt swap, strand flip) rather than by raw id-string equality.
#' This makes matching robust to chr-prefix, field-separator, and allele-order
#' differences -- the id string is a serialization, not an identity.
#'
#' Allele semantics mirror summary-statistics QC (it delegates the allele
#' classification to \code{harmonizeAlleles}): an allele swap is matched and
#' reported with \code{sign = -1} so callers can sign-flip a paired effect / z
#' / weight; strand flips are resolved; and strand-ambiguous (A/T, C/G) matches
#' are dropped when \code{removeStrandAmbiguous = TRUE}.
#'
#' When ids cannot be parsed into chrom/pos/ref/alt (e.g. rsIDs), the matcher
#' falls back to exact string identity with \code{sign = 1} (no allele
#' awareness is possible without alleles).
#'
#' @param idsA,idsB Character vectors of variant IDs, or data.frames with
#'   chrom/pos/A2/A1 columns (A2 = ref, A1 = alt).
#' @param allowFlip When TRUE (default) match by (chrom, pos) with ref/alt swap
#'   and strand handling, reporting \code{sign = -1} for a swap. When FALSE,
#'   match on exact alleles only (canonicalized format; \code{sign} always +1)
#'   so a ref/alt swap does not match.
#' @param removeStrandAmbiguous Drop strand-ambiguous (A/T, C/G) matches.
#'   Default TRUE (mirrors summaryStatsQc).
#' @param removeIndels Drop indels. Default FALSE.
#' @return A list of three equal-length vectors describing the matched pairs:
#'   \code{idxA} / \code{idxB} (indices into \code{idsA} / \code{idsB}) and
#'   \code{sign} (+1 exact, -1 allele-swap). At most one pair per \code{idsA}
#'   entry; unmatched entries are omitted.
#' @noRd
matchVariants <- function(idsA, idsB, allowFlip = TRUE,
                          removeStrandAmbiguous = TRUE,
                          removeIndels = FALSE) {
  empty <- list(idxA = integer(0), idxB = integer(0), sign = numeric(0))
  if (!allowFlip) {
    # Exact-allele matching: canonicalize the id FORMAT (chr-prefix + separator,
    # rsID-safe) but keep allele order, then match by exact string identity, so
    # a ref/alt swap does NOT match (sign is always +1).
    keyOf <- function(x) if (is.data.frame(x)) {
        p <- parseVariantId(x)
        formatVariantId(p$chrom, p$pos, p$A2, p$A1)
      } else normalizeVariantId(x)
    idxB <- match(keyOf(idsA), keyOf(idsB))
    matched <- which(!is.na(idxB))
    return(list(idxA = matched, idxB = idxB[matched],
                sign = rep(1, length(matched))))
  }
  dfA <- parseVariantId(idsA)
  dfB <- parseVariantId(idsB)
  unparseable <- nrow(dfA) == 0L || nrow(dfB) == 0L ||
    anyNA(dfA$chrom) || anyNA(dfA$pos) ||
    anyNA(dfB$chrom) || anyNA(dfB$pos)
  if (unparseable) {
    # rsID / unparseable ids: exact string identity only (no allele awareness).
    if (is.data.frame(idsA) || is.data.frame(idsB)) return(empty)
    idxB <- match(as.character(idsA), as.character(idsB))
    matched <- which(!is.na(idxB))
    return(list(idxA = matched, idxB = idxB[matched],
                sign = rep(1, length(matched))))
  }
  # Inject sentinel index/sign columns so the matched pairs and the swap sign
  # can be read straight back out of harmonizeAlleles without re-deriving them.
  dfA$.mvTidx <- seq_len(nrow(dfA))
  dfA$.mvSign <- 1
  dfB$.mvRidx <- seq_len(nrow(dfB))
  res <- suppressWarnings(harmonizeAlleles(
    targetData = dfA, refVariants = dfB,
    colToFlip = ".mvSign", matchMinProp = 0, removeDups = TRUE,
    flipStrand = FALSE, removeIndels = removeIndels,
    removeStrandAmbiguous = removeStrandAmbiguous, removeUnmatched = TRUE))
  h <- res$harmonizedData
  if (is.null(h) || nrow(h) == 0L) return(empty)
  keep <- !duplicated(h$.mvTidx)          # at most one ref per A id
  list(idxA = as.integer(h$.mvTidx[keep]),
       idxB = as.integer(h$.mvRidx[keep]),
       sign = as.numeric(h$.mvSign[keep]))
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
    chrom = canonChrom(parts[1]),
    start = as.integer(parts[2]),
    end = as.integer(parts[3])
  )

  return(df)
}

#' Utility function to convert LD region_ids to `region of interest` dataframe
#'
#' The first field is treated as the chromosome and kept as a normalized string
#' (via \code{canonChrom}) so X/Y/MT survive; the remaining fields are genomic
#' positions and are returned as integers.
#' @param ldRegionId A string of region in the format of chrom_start_end.
#' @export
regionToDf <- function(ldRegionId, colnames = c("chrom", "start", "end")) {
  parts <- do.call(rbind, strsplit(ldRegionId, "[_:-]"))
  regionOfInterest <- as.data.frame(parts, stringsAsFactors = FALSE)
  colnames(regionOfInterest) <- colnames
  regionOfInterest[[1]] <- canonChrom(regionOfInterest[[1]])
  for (j in seq_along(colnames)[-1]) {
    regionOfInterest[[j]] <- as.integer(regionOfInterest[[j]])
  }
  regionOfInterest
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

