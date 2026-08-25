#' @title Variant ID Parsing and Formatting Utilities
#' @description Functions for parsing, formatting, normalizing, and detecting
#'   the naming conventions of variant IDs (e.g., "chr1:100:A:G") and genomic
#'   region strings (e.g., "chr1:100-200").
#' @name pecotmr-variant-id
#' @keywords internal
#' @importFrom stringr str_split str_split_1 str_remove str_detect str_to_upper
#'   str_length str_c str_replace_all str_match str_flatten regex
#' @importFrom tidyr replace_na
#' @importFrom tibble as_tibble
#' @importFrom dplyr bind_cols
#' @importFrom magrittr is_in
NULL

#' Strip "chr" prefix from chromosome identifiers.
#' @param x Character vector of chromosome identifiers (e.g., "chr1", "chrX").
#' @return Character vector with "chr" prefix removed (e.g., "1", "X").
#' @noRd
stripChrPrefix <- function(x) str_remove(x, "^chr")

#' Canonicalize chromosome identifiers to a normalized string.
#'
#' Strips a leading "chr"/"ch" prefix (case-insensitive), uppercases, and maps
#' common synonyms to a single form (23 -> X, 24 -> Y, M -> MT). Returns a
#' character vector and never coerces to integer, so non-autosomal chromosomes
#' (X, Y, MT) survive instead of collapsing to NA. This is the single chromosome
#' normalizer used for variant matching and identity throughout the package; it
#' mirrors the chromosome handling previously inlined in summary-statistics QC.
#'
#' @param x Character or numeric vector of chromosome identifiers.
#' @return Character vector of normalized chromosome names.
#' @noRd
canonChrom <- function(x) {
    x <- as.character(x)
    x <- str_remove(x, regex("^chr", ignore_case = TRUE))
    x <- str_remove(x, regex("^ch", ignore_case = TRUE))
    x <- str_to_upper(x)
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
    str_c("chr", str_remove(as.character(x), regex("^chr", ignore_case = TRUE)))
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
    standard <- c(as.character(seq_len(22)), "X", "Y", "XY", "MT")
    extra <- setdiff(unique(x[!is.na(x)]), standard)
    factor(x, levels = c(standard, extra), ordered = TRUE)
}

# Backwards-compat alias

#' Strip build suffix from variant IDs (e.g., ":b38" or "_b38").
#' @param x Character vector of variant IDs.
#' @return Character vector with build suffix removed.
#' @noRd
stripBuildSuffix <- function(x) str_remove(x, "(:|_)b[0-9]+$")

# Backwards-compat alias

#' Test whether allele pairs are single-nucleotide (SNP, not indel).
#'
#' Returns TRUE for each pair where both alleles are exactly one of A, T, C, G.
#' @param a1 Character vector of first alleles.
#' @param a2 Character vector of second alleles.
#' @return Logical vector, TRUE if the variant is a SNP.
#' @noRd
isSnpAlleles <- function(a1, a2) {
    isSnp <- str_length(a1) == 1L &
        str_length(a2) == 1L &
        str_detect(a1, "^[ATCG]$") &
        str_detect(a2, "^[ATCG]$")
    # str_length(NA)/str_detect(NA) yield NA where base nchar/grepl gave FALSE;
    # preserve the base "NA allele is not a SNP" behavior.
    replace_na(isSnp, FALSE)
}

# Backwards-compat alias

#' Detect the naming convention of variant IDs
#'
#' Examines variant ID strings to detect their format: whether they have a "chr"
#' prefix, what separator is used between allele fields, and whether they
#' include a genome build suffix (e.g., ":b38" or "_b38").
#'
#' Supported formats:
#' \itemize{
#'   \item All colons: \code{"chr1:100:A:G"} or \code{"1:100:A:G"}
#'   \item Mixed colon/underscore: \code{"chr1:100_A_G"} or \code{"1:100_A_G"}
#'   \item All underscores (PLINK BIM): \code{"chr1_100_A_G"} or
#'   \code{"1_100_A_G"}
#' }
#'
#' @param ids A character vector of variant IDs.
#' @return A list with components:
#'   \describe{
#'     \item{hasChr}{Logical, whether the IDs have a "chr" prefix.}
#'     \item{alleleSep}{Character, the separator between allele fields (":" or
#'     "_").
#'       For mixed format \code{"chr1:100_A_G"}, this is \code{"_"}.}
#'     \item{hasBuild}{Logical, whether a build suffix is present.}
#'     \item{example}{Character, the first non-NA ID for reference.}
#'   }
#' @noRd
detectVariantConvention <- function(ids) {
    # Find first non-NA element
    firstId <- ids[!is.na(ids)][1]
    if (is.na(firstId) || length(firstId) == 0) {
        return(list(
            hasChr = FALSE,
            alleleSep = ":",
            hasBuild = FALSE,
            example = NA_character_
        ))
    }
    hasChr <- str_detect(firstId, "^chr")
    # Detect build suffix like :b38 or _b38 at end
    hasBuild <- str_detect(firstId, "(:|_)b[0-9]+$")
    idClean <- stripBuildSuffix(firstId)
    # Detect allele separator: check if variant uses underscores between allele
    # fields This catches both full underscore ("1_100_A_G") and mixed
    # ("chr1:100_A_G") formats
    alleleSep <- if (str_detect(idClean, "_[ATCGID*]+_[ATCGID*]+$")) {
        "_"
    } else {
        ":"
    }
    list(
        hasChr = hasChr,
        alleleSep = alleleSep,
        hasBuild = hasBuild,
        example = firstId
    )
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
#'   "chrom", "pos", and allele columns (A2/A1 or ref/alt or any 4-column
#'   layout).
#' @return A tibble with columns "chrom" (character, normalized), "pos"
#'   (integer), "A2" (character), "A1" (character). The detected convention is
#'   stored as \code{attr(result, "convention")}.
#' @examples
#' parseVariantId(c("chr1:100:A:G", "chr2:200:T:C"))
#' @export
parseVariantId <- function(ids) {
    if (is.data.frame(ids)) {
        return(.parseVariantIdDf(ids))
    }
    convention <- detectVariantConvention(ids)
    # Normalize: convert underscores to colons, strip build suffix.
    normalized <- stripBuildSuffix(str_replace_all(ids, "_", ":"))
    # Split into exactly 4 fields via a single vectorized regex match; a
    # non-matching id yields an all-NA capture row (mirrors strcapture).
    m <- str_match(normalized, "^([^:]+):([^:]+):([^:]+):([^:]+)")
    data <- tibble(
        chrom = canonChrom(m[, 2L]),
        pos = as.integer(m[, 3L]),
        A2 = m[, 4L],
        A1 = m[, 5L]
    )
    attr(data, "convention") <- convention
    data
}

# Parse a data.frame of already-split ids: resolve the 4 identity columns
# (explicit chrom/pos/A2/A1 or chrom/pos/A1/A2 kept as-is; otherwise positional
# chrom/pos/A2/A1), then canonicalize chrom/pos.
# @noRd
.parseVariantIdDf <- function(ids) {
    # minimal repair: preserve any empty/duplicate extra-column names (e.g. an
    # unnamed passthrough column) for .sanitizeNames() to canonicalize later.
    ids <- as_tibble(ids, .name_repair = "minimal")
    hasA2A1 <- all(is_in(c("chrom", "pos", "A2", "A1"), names(ids)))
    hasA1A2 <- all(is_in(c("chrom", "pos", "A1", "A2"), names(ids)))
    if (!hasA2A1 && !hasA1A2 && ncol(ids) >= 4) {
        names(ids)[seq_len(4)] <- c("chrom", "pos", "A2", "A1")
    }
    conv <- list(
        hasChr = any(str_detect(as.character(ids$chrom), "^chr")),
        alleleSep = ":",
        hasBuild = FALSE,
        example = NA_character_
    )
    ids$chrom <- canonChrom(ids$chrom)
    ids$pos <- as.integer(ids$pos)
    attr(ids, "convention") <- conv
    ids
}

#' Format variant ID strings from component columns
#'
#' Constructs variant ID strings from chrom, pos, A2, A1 columns. The chrom:pos
#' separator is always a colon. The allele separator can be either colon
#' (canonical: \code{"chr1:100:A:G"}) or underscore (mixed:
#' \code{"chr1:100_A_G"}).
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
#' @param chrPrefix Logical, whether to add "chr" prefix. Default TRUE. Ignored
#'   if \code{convention} is provided.
#' @param alleleSep Character, separator between pos/A2 and A2/A1 fields.
#'   Default \code{":"} produces canonical \code{"chr1:100:A:G"}; \code{"_"}
#'   produces mixed \code{"chr1:100_A_G"}. Ignored if \code{convention} is
#'   provided.
#' @param convention Optional list from \code{detectVariantConvention}. When
#'   provided, \code{hasChr} and \code{alleleSep} are read from the convention
#'   automatically. This is the preferred way to preserve the user's input
#'   format.
#' @return A character vector of formatted variant IDs.
#' @noRd
formatVariantId <- function(
    chrom,
    pos,
    A2,
    A1,
    chrPrefix = TRUE,
    alleleSep = ":",
    convention = NULL
) {
    # If convention is provided, use it to determine format automatically
    if (!is.null(convention)) {
        chrPrefix <- convention$hasChr
        alleleSep <- if (!is.null(convention$alleleSep)) {
            convention$alleleSep
        } else {
            ":"
        }
    }
    # Normalize the chromosome (strip prefix, uppercase, map synonyms) then
    # re-add the "chr" prefix if requested. canonChrom keeps X/Y/MT as strings.
    chromClean <- canonChrom(chrom)
    if (chrPrefix) {
        str_c("chr", chromClean, ":", pos, alleleSep, A2, alleleSep, A1)
    } else {
        str_c(chromClean, ":", pos, alleleSep, A2, alleleSep, A1)
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
#'   \code{attr(parseVariantId(ids), "convention")}. When provided, the output
#'   format is driven automatically by the detected convention.
#' @return A character vector of re-formatted variant IDs. Unparseable ids (e.g.
#'   rsIDs) are returned unchanged.
#' @examples
#' normalizeVariantId(c("1:100:A:G", "2:200:T:C"))
#' @export
normalizeVariantId <- function(ids, chrPrefix = TRUE, convention = NULL) {
    parsed <- parseVariantId(ids)
    out <- as.character(ids)
    # Only re-format ids that parsed into a chrom + pos; leave unparseable ids
    # (e.g. rsIDs) unchanged rather than emitting "chrNA:..." garbage.
    ok <- !is.na(parsed$chrom) & !is.na(parsed$pos)
    if (any(ok)) {
        out[ok] <- if (!is.null(convention)) {
            formatVariantId(
                parsed$chrom[ok],
                parsed$pos[ok],
                parsed$A2[ok],
                parsed$A1[ok],
                convention = convention
            )
        } else {
            formatVariantId(
                parsed$chrom[ok],
                parsed$pos[ok],
                parsed$A2[ok],
                parsed$A1[ok],
                chrPrefix = chrPrefix
            )
        }
    }
    out
}

#' Parse variant IDs into a data frame
#'
#' Convenience wrapper around \code{\link{parseVariantId}} returning the parsed
#' \code{chrom}/\code{pos}/\code{A2}/\code{A1} data frame for a vector of
#' variant IDs.
#'
#' @param variantId A character vector of variant IDs.
#' @return A data.frame with columns \code{chrom}, \code{pos}, \code{A2},
#'   \code{A1}.
#' @seealso \code{\link{parseVariantId}}
#' @examples
#' variantIdToDf("chr1:100:A:G")
#' @export
variantIdToDf <- function(variantId) {
    parseVariantId(variantId)
}

# Complement a DNA allele string (A<->T, C<->G) for strand flipping.
# @noRd
.strandFlip <- function(ref) chartr("ATCG", "TAGC", ref)

# Ensure a data.frame has unique, non-empty column names (blanks become
# `unnamed_<i>`, duplicates de-duplicated with make.unique).
# @noRd
.sanitizeNames <- function(df) {
    nm <- colnames(df)
    if (is.null(nm)) {
        nm <- rep("unnamed", ncol(df))
    }
    emptyIdx <- is.na(nm) | nm == ""
    if (any(emptyIdx)) {
        nm[emptyIdx] <- str_c("unnamed_", seq_len(sum(emptyIdx)))
    }
    colnames(df) <- make.unique(nm, sep = "_")
    df
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
#' @param targetData A data frame with columns "chrom", "pos", "A2", "A1" (and
#'   optionally other columns like "beta" or "z"), or a vector of strings in the
#'   format of "chr:pos:A2:A1"/"chr:pos_A2_A1". Can be automatically converted
#'   to a data frame if a vector.
#' @param refVariants A data frame with columns "chrom", "pos", "A2", "A1" or
#'   strings in the format of "chr:pos:A2:A1"/"chr:pos_A2_A1".
#' @param colToFlip The name of the column in targetData where flips are to be
#'   applied. On an allele swap these columns are sign-flipped (multiplied by
#'   -1), the correct operation for signed quantities like \code{beta} and
#'   \code{z}.
#' @param colToComplement Names of columns in targetData to complement (\code{1
#'   - x}) on an allele swap, the correct operation for an effect-allele
#'   frequency like \code{af}. Default \code{character()} does no complementing,
#'   so non-RSS callers are unchanged. Distinct from \code{colToFlip}:
#'   frequencies are complemented, signed effects are sign-flipped.
#' @param matchMinProp Minimum proportion of variants in the smallest data to be
#'   matched, otherwise stops with an error. Default is 20%.
#' @param removeDups Whether to remove duplicates, default is TRUE.
#' @param removeIndels Whether to remove INDELs, default is FALSE.
#' @param flip Whether the alleles must be flipped: A <--> T & C <--> G, in
#'   which case corresponding `colToFlip` are multiplied by -1. Default is
#'   `TRUE`.
#' @param removeStrandAmbiguous Whether to remove strand SNPs (if any). Default
#'   is `TRUE`.
#' @param flipStrand Whether to output the variants after strand flip. Default
#'   is `FALSE`.
#' @param removeUnmatched Whether to remove unmatched variants. Default is
#'   `TRUE`.
#' @return An \code{AlleleQcResult} S4 object. Use \code{$harmonizedData} to
#'   recover the post-QC variant data.frame and \code{$qcSummary} to inspect the
#'   per-variant merge/flip/strand diagnostics.
#' @importFrom dplyr mutate inner_join filter pull select everything row_number
#' @importFrom dplyr if_else any_of all_of rename across
#' @importFrom vctrs vec_duplicate_detect
#' @importFrom tidyr separate
#' @keywords internal
#' @noRd
#' @details Pure panel-vs-sumstats allele harmonization: match by (chrom, pos),
#'   detect A1/A2 swap, sign-flip \code{colToFlip} columns and complement
#'   \code{colToComplement} columns on swap. Variant-allele filters (indels,
#'   strand-ambiguous, duplicates) are applied here directly when the
#'   corresponding \code{removeIndels} / \code{removeStrandAmbiguous} /
#'   \code{removeDups} flags are set; MAF / INFO / N column-numeric filters run
#'   in \code{.applyContentFilters()} before this function.
harmonizeAlleles <- function(
    targetData,
    refVariants,
    colToFlip = NULL,
    matchMinProp = 0.2,
    flipStrand = FALSE,
    removeUnmatched = TRUE,
    removeIndels = FALSE,
    removeStrandAmbiguous = TRUE,
    removeDups = FALSE,
    colToComplement = character(),
    ...
) {
    coerced <- .harmonizeCoerceInputs(targetData, refVariants)
    targetData <- coerced$targetData
    refVariants <- coerced$refVariants
    matchResult <- .harmonizeJoin(targetData, refVariants)
    if (nrow(matchResult) == 0) {
        return(.harmonizeEmptyResult(matchResult))
    }
    matchResult <- .harmonizeFlags(matchResult)
    matchResult <- .harmonizeResolveAmbiguity(
        matchResult,
        removeStrandAmbiguous
    )
    matchResult <- .harmonizeKeepRule(matchResult, removeIndels)
    matchResult <- .harmonizeApplyFlips(
        matchResult,
        colToFlip,
        colToComplement,
        flipStrand
    )
    qcCounts <- .harmonizeQcCounts(matchResult)
    qcSummary <- matchResult
    result <- .harmonizeCleanResult(matchResult)
    if (removeDups) {
        result <- .harmonizeRemoveDups(result)
    }
    if (!removeUnmatched) {
        restored <- .harmonizeRestoreUnmatched(result, matchResult, targetData)
        result <- restored$result
        qcSummary <- restored$qcSummary
    }
    .harmonizeFinalChecks(result, refVariants, matchMinProp)
    out <- list(harmonizedData = result, qcSummary = qcSummary)
    attr(out, "qcCounts") <- qcCounts
    out
}

# QC / flag columns stripped from the harmonized result before it is returned.
.harmonizeQcCols <- c(
    "flip1.ref",
    "flip2.ref",
    "strand_unambiguous",
    "exact_match",
    "sign_flip",
    "strand_flip",
    "INDEL",
    "ID_match",
    "keep"
)

# Coerce both sides to canonical variant data.frames and strip merge-conflicting
# columns (variant_id is dropped on both sides because it is rebuilt from the
# QC'd alleles below; leaving it in collides on the final rename).
# @noRd
.harmonizeCoerceInputs <- function(targetData, refVariants) {
    if (
        is.data.frame(targetData) &&
            ncol(targetData) > 4 &&
            all(is_in(c("chrom", "pos", "A2", "A1"), names(targetData)))
    ) {
        variantCols <- c("chrom", "pos", "A2", "A1")
        variantDf <- targetData |> select(all_of(variantCols))
        otherCols <- targetData |> select(-all_of(variantCols))
        targetData <- bind_cols(
            variantIdToDf(variantDf),
            otherCols,
            .name_repair = "minimal"
        )
    } else {
        targetData <- variantIdToDf(targetData)
    }
    refVariants <- variantIdToDf(refVariants)
    dropCols <- c("chromosome", "position", "ref", "alt", "variant_id")
    if (any(is_in(dropCols, colnames(targetData)))) {
        targetData <- select(targetData, -any_of(dropCols))
    }
    if (is_in("variant_id", colnames(refVariants))) {
        refVariants <- select(refVariants, -any_of("variant_id"))
    }
    list(targetData = targetData, refVariants = refVariants)
}

# Inner-join target + reference on (chrom, pos).
# @noRd
.harmonizeJoin <- function(targetData, refVariants) {
    inner_join(
        targetData,
        refVariants,
        by = c("chrom", "pos"),
        suffix = c(".target", ".ref")
    ) |>
        as_tibble(.name_repair = "minimal") |>
        .sanitizeNames()
}

# Empty-match early return (warning + zeroed qcCounts).
# @noRd
.harmonizeEmptyResult <- function(matchResult) {
    msg <- glue(
        "No matching variants found between target data and ",
        "reference variants."
    )
    warn(msg)
    emptyOut <- list(harmonizedData = matchResult, qcSummary = matchResult)
    attr(emptyOut, "qcCounts") <- list(
        considered = 0L,
        signFlip = 0L,
        strandFlip = 0L,
        kept = 0L,
        dropped = 0L,
        droppedIndel = 0L,
        droppedAmbiguous = 0L,
        droppedOther = 0L
    )
    emptyOut
}

# Per-variant harmonization flags: original/QC'd ids, uppercased alleles, strand
# complements, and the exact/sign-flip/strand-flip/INDEL/ID-match indicators.
# @noRd
.harmonizeFlags <- function(matchResult) {
    matchResult |>
        mutate(
            variants_id_original = formatVariantId(
                .data$chrom,
                .data$pos,
                .data$A2.target,
                .data$A1.target
            ),
            variants_id_qced = formatVariantId(
                .data$chrom,
                .data$pos,
                .data$A2.ref,
                .data$A1.ref
            )
        ) |>
        mutate(across(
            c("A1.target", "A2.target", "A1.ref", "A2.ref"),
            str_to_upper
        )) |>
        mutate(
            flip1.ref = .strandFlip(.data$A1.ref),
            flip2.ref = .strandFlip(.data$A2.ref)
        ) |>
        .harmonizeMatchFlags()
}

# Exact / sign-flip / strand-flip / INDEL / ID-match indicators, consuming the
# uppercased alleles and strand complements from the earlier pipe stages.
# @noRd
.harmonizeMatchFlags <- function(matchResult) {
    matchResult |>
        mutate(
            strand_unambiguous = if_else(
                (.data$A1.target == "A" & .data$A2.target == "T") |
                    (.data$A1.target == "T" & .data$A2.target == "A") |
                    (.data$A1.target == "C" & .data$A2.target == "G") |
                    (.data$A1.target == "G" & .data$A2.target == "C"),
                FALSE,
                TRUE
            )
        ) |>
        mutate(
            exact_match = .data$A1.target == .data$A1.ref &
                .data$A2.target == .data$A2.ref
        ) |>
        mutate(
            sign_flip = ((.data$A1.target == .data$A2.ref &
                .data$A2.target == .data$A1.ref) |
                (.data$A1.target == .data$flip2.ref &
                    .data$A2.target == .data$flip1.ref)) &
                (.data$A1.target != .data$A1.ref &
                    .data$A2.target != .data$A2.ref)
        ) |>
        mutate(
            strand_flip = ((.data$A1.target == .data$flip1.ref &
                .data$A2.target == .data$flip2.ref) |
                (.data$A1.target == .data$flip2.ref &
                    .data$A2.target == .data$flip1.ref)) &
                (.data$A1.target != .data$A1.ref &
                    .data$A2.target != .data$A2.ref)
        ) |>
        mutate(
            INDEL = (.data$A2.target == "I" |
                .data$A2.target == "D" |
                str_length(.data$A2.target) > 1L |
                str_length(.data$A1.target) > 1L)
        ) |>
        mutate(
            ID_match = ((.data$A2.target == "D" | .data$A2.target == "I") &
                (str_length(.data$A1.ref) > 1L |
                    str_length(.data$A2.ref) > 1L))
        )
}

# Strand-ambiguity resolution: disable the A/T-C/G guard when the caller opted
# out, and when no unambiguous strand flip survives (ambiguous variants then
# fall through as exact / sign-flip cases rather than being dropped).
# @noRd
.harmonizeResolveAmbiguity <- function(matchResult, removeStrandAmbiguous) {
    if (!removeStrandAmbiguous) {
        matchResult$strand_unambiguous <- TRUE
    }
    if (!any(matchResult$strand_flip & matchResult$strand_unambiguous)) {
        matchResult$strand_unambiguous <- TRUE
    }
    matchResult
}

# Compute the keep flag (strand-flip vs non-strand-flip rules); drop indels when
# requested.
# @noRd
.harmonizeKeepRule <- function(matchResult, removeIndels) {
    matchResult <- matchResult |>
        mutate(
            keep = if_else(
                .data$strand_flip,
                true = .data$strand_unambiguous |
                    .data$exact_match |
                    .data$ID_match,
                false = .data$exact_match |
                    .data$sign_flip |
                    .data$ID_match
            )
        )
    if (removeIndels) {
        matchResult <- matchResult |>
            mutate(keep = if_else(.data$INDEL, FALSE, .data$keep))
    }
    matchResult
}

# Named per-row conditional column transforms for the harmonize across() calls
# (the `flip` condition vector is passed through across's `...`, so no anonymous
# functions are needed).
# @noRd
.negateWhere <- function(x, flip) if_else(flip, -x, x)
# @noRd
.complementWhere <- function(x, flip) if_else(flip, 1 - x, x)
# @noRd
.strandFlipWhere <- function(x, flip) if_else(flip, .strandFlip(x), x)

# Apply signed-column flips (colToFlip), effect-allele-frequency complements
# (colToComplement, af -> 1 - af on a swap), and optional target strand flips.
# @noRd
.harmonizeApplyFlips <- function(
    matchResult,
    colToFlip,
    colToComplement,
    flipStrand
) {
    if (!is.null(colToFlip)) {
        .harmonizeCheckCols(colToFlip, matchResult)
        matchResult <- matchResult |>
            mutate(across(
                all_of(colToFlip),
                .negateWhere,
                matchResult$sign_flip
            ))
    }
    if (length(colToComplement) > 0L) {
        .harmonizeCheckCols(colToComplement, matchResult)
        matchResult <- matchResult |>
            mutate(across(
                all_of(colToComplement),
                .complementWhere,
                matchResult$sign_flip
            ))
    }
    if (flipStrand) {
        matchResult <- .harmonizeFlipStrandCols(matchResult)
    }
    matchResult
}

# Assert the named columns exist in matchResult.
# @noRd
.harmonizeCheckCols <- function(cols, matchResult) {
    missing <- setdiff(cols, colnames(matchResult))
    if (length(missing) > 0L) {
        joined <- str_flatten(missing, "', '")
        msg <- glue("Column(s) '{joined}' not found in targetData.")
        abort(msg)
    }
    invisible(NULL)
}

# Strand-flip the target alleles of the strand-flipped rows.
# @noRd
.harmonizeFlipStrandCols <- function(matchResult) {
    matchResult |>
        mutate(across(
            c("A1.target", "A2.target"),
            .strandFlipWhere,
            matchResult$strand_flip
        ))
}

# Per-step QC counts (for the "kept N of M (corrected: ...; dropped ...)" logs),
# computed before the flag columns are stripped from the returned frame.
# @noRd
.harmonizeQcCounts <- function(matchResult) {
    hasIndel <- is_in("INDEL", colnames(matchResult))
    qcCounts <- list(
        considered = nrow(matchResult),
        signFlip = sum(matchResult$sign_flip & matchResult$keep, na.rm = TRUE),
        strandFlip = sum(
            matchResult$strand_flip & matchResult$keep,
            na.rm = TRUE
        ),
        kept = sum(matchResult$keep, na.rm = TRUE),
        dropped = sum(!matchResult$keep, na.rm = TRUE),
        droppedIndel = if (hasIndel) {
            sum(!matchResult$keep & matchResult$INDEL, na.rm = TRUE)
        } else {
            0L
        }
    )
    qcCounts$droppedAmbiguous <- sum(
        !matchResult$keep &
            matchResult$strand_flip &
            !matchResult$strand_unambiguous &
            if (hasIndel) !matchResult$INDEL else TRUE,
        na.rm = TRUE
    )
    qcCounts$droppedOther <- qcCounts$dropped -
        qcCounts$droppedIndel -
        qcCounts$droppedAmbiguous
    qcCounts
}

# Kept rows with QC/flag + target-allele columns stripped and ref alleles /
# QC'd id renamed to the canonical A1 / A2 / variant_id.
# @noRd
.harmonizeCleanResult <- function(matchResult) {
    matchResult |>
        filter(.data$keep) |>
        select(
            -any_of(.harmonizeQcCols),
            -any_of(c("A1.target", "A2.target"))
        ) |>
        rename(
            A1 = "A1.ref",
            A2 = "A2.ref",
            variant_id = "variants_id_qced"
        )
}

# Drop duplicate (chrom, pos, variant_id) rows, keeping the first occurrence.
# @noRd
.harmonizeRemoveDups <- function(result) {
    deduped <- distinct(
        result,
        .data$chrom,
        .data$pos,
        .data$variant_id,
        .keep_all = TRUE
    )
    nDropped <- nrow(result) - nrow(deduped)
    if (nDropped > 0) {
        msg <- glue(
            "Removed {nDropped} duplicate variant(s), keeping first ",
            "occurrence."
        )
        warn(msg)
    }
    deduped
}

# removeUnmatched = FALSE path: re-append the unmatched target variants in the
# original target order. Returns list(result, qcSummary) -- qcSummary is the
# cleaned/renamed matchResult (matching the original's returned qcSummary here).
# @noRd
.harmonizeRestoreUnmatched <- function(result, matchResult, targetData) {
    matchVariant <- result |> pull("variants_id_original")
    qcSummary <- matchResult |>
        select(
            -any_of(.harmonizeQcCols),
            -any_of(c("variants_id_original", "A1.target", "A2.target"))
        ) |>
        rename(
            A1 = "A1.ref",
            A2 = "A2.ref",
            variant_id = "variants_id_qced"
        )
    targetData <- targetData |>
        mutate(
            variant_id = formatVariantId(
                .data$chrom,
                .data$pos,
                .data$A2,
                .data$A1
            )
        )
    if (length(setdiff(targetData |> pull("variant_id"), matchVariant)) == 0L) {
        return(list(result = result, qcSummary = qcSummary))
    }
    unmatchData <- targetData |> filter(!is_in(.data$variant_id, matchVariant))
    result <- bind_rows(
        result,
        unmatchData |> mutate(variants_id_original = .data$variant_id)
    )
    result <- result |>
        slice(match(targetData$variant_id, .data$variants_id_original)) |>
        select(-any_of("variants_id_original"))
    list(result = result, qcSummary = qcSummary)
}

# Final guards: enough variants matched, and no duplicate ids remain.
# @noRd
.harmonizeFinalChecks <- function(result, refVariants, matchMinProp) {
    if (nrow(result) < matchMinProp * nrow(refVariants)) {
        abort("Not enough variants have been matched.")
    }
    if (any(duplicated(result$variant_id))) {
        msg <- glue(
            "Duplicated variant IDs remain after harmonization; pass ",
            "removeDups = TRUE or deduplicate upstream before calling ",
            "harmonizeAlleles."
        )
        abort(msg)
    }
    invisible(NULL)
}

# Canonical per-id key for allele-aware matching: data.frame -> chrom/pos/A2/A1
# formatted id (keeps allele order, so a ref/alt swap does NOT match); character
# vector -> normalized id.
# @noRd
.matchVariantKeyOf <- function(x) {
    if (is.data.frame(x)) {
        p <- parseVariantId(x)
        formatVariantId(p$chrom, p$pos, p$A2, p$A1)
    } else {
        normalizeVariantId(x)
    }
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
#' reported with \code{sign = -1} so callers can sign-flip a paired effect / z /
#' weight; strand flips are resolved; and strand-ambiguous (A/T, C/G) matches
#' are dropped when \code{removeStrandAmbiguous = TRUE}.
#'
#' When ids cannot be parsed into chrom/pos/ref/alt (e.g. rsIDs), the matcher
#' falls back to exact string identity with \code{sign = 1} (no allele awareness
#' is possible without alleles).
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
matchVariants <- function(
    idsA,
    idsB,
    allowFlip = TRUE,
    removeStrandAmbiguous = TRUE,
    removeIndels = FALSE
) {
    if (!allowFlip) {
        # Exact-allele matching: canonicalize the id FORMAT (chr-prefix +
        # separator, rsID-safe) but keep allele order, then match by exact
        # string identity, so a ref/alt swap does NOT match (sign is always +1).
        return(.matchVariantsExact(
            .matchVariantKeyOf(idsA),
            .matchVariantKeyOf(idsB)
        ))
    }
    dfA <- parseVariantId(idsA)
    dfB <- parseVariantId(idsB)
    if (.matchVariantsUnparseable(dfA, dfB)) {
        # rsID / unparseable ids: exact string identity only (no allele
        # awareness). data.frame inputs have no string form to fall back to.
        if (is.data.frame(idsA) || is.data.frame(idsB)) {
            return(.matchVariantsEmpty())
        }
        return(.matchVariantsExact(as.character(idsA), as.character(idsB)))
    }
    .matchVariantsHarmonized(dfA, dfB, removeIndels, removeStrandAmbiguous)
}

# The empty match result.
# @noRd
.matchVariantsEmpty <- function() {
    list(idxA = integer(0), idxB = integer(0), sign = numeric(0))
}

# Exact string match (sign always +1) between two pre-computed key vectors.
# @noRd
.matchVariantsExact <- function(keyA, keyB) {
    idxB <- match(keyA, keyB)
    matched <- which(!is.na(idxB))
    list(
        idxA = matched,
        idxB = idxB[matched],
        sign = rep(1, length(matched))
    )
}

# TRUE when either parsed side lacks usable chrom/pos (rsID / unparseable ids).
# @noRd
.matchVariantsUnparseable <- function(dfA, dfB) {
    nrow(dfA) == 0L ||
        nrow(dfB) == 0L ||
        anyNA(dfA$chrom) ||
        anyNA(dfA$pos) ||
        anyNA(dfB$chrom) ||
        anyNA(dfB$pos)
}

# Allele-aware match via harmonizeAlleles, reading the matched pairs + swap sign
# back out of the injected sentinel columns.
# @noRd
.matchVariantsHarmonized <- function(
    dfA,
    dfB,
    removeIndels,
    removeStrandAmbiguous
) {
    # Inject sentinel index/sign columns so the matched pairs and the swap sign
    # can be read straight back out of harmonizeAlleles without re-deriving
    # them.
    dfA$.mvTidx <- seq_len(nrow(dfA))
    dfA$.mvSign <- 1
    dfB$.mvRidx <- seq_len(nrow(dfB))
    res <- suppressWarnings(harmonizeAlleles(
        targetData = dfA,
        refVariants = dfB,
        colToFlip = ".mvSign",
        matchMinProp = 0,
        removeDups = TRUE,
        flipStrand = FALSE,
        removeIndels = removeIndels,
        removeStrandAmbiguous = removeStrandAmbiguous,
        removeUnmatched = TRUE
    ))
    h <- res$harmonizedData
    if (is.null(h) || nrow(h) == 0L) {
        return(.matchVariantsEmpty())
    }
    h <- .matchDropSignConflicts(h)
    if (nrow(h) == 0L) {
        return(.matchVariantsEmpty())
    }
    keep <- !duplicated(h$.mvTidx) # at most one ref per A id
    list(
        idxA = as.integer(h$.mvTidx[keep]),
        idxB = as.integer(h$.mvRidx[keep]),
        sign = as.numeric(h$.mvSign[keep])
    )
}

# Drop target variants whose reference matches disagree about the sign.
#
# A reference panel can legitimately carry a variant and its own allele flip as
# two separate entries -- two distinct indels at one position that happen to be
# each other's flip is rare but biologically real. A target variant then
# matches BOTH: once exactly (sign +1) and once as an allele swap (sign -1),
# and there is no way to tell which entry it is. The caller below keeps the
# first match, so without this the answer would depend on the panel's row
# order, silently flipping the effect direction for that variant.
#
# Only a sign DISAGREEMENT is ambiguous. Duplicate reference entries that agree
# (the same variant listed twice) still resolve to the first match, so this
# drops nothing it does not have to.
# @noRd
.matchDropSignConflicts <- function(h) {
    ambiguous <- h |>
        group_by(.data$.mvTidx) |>
        summarise(nSigns = n_distinct(.data$.mvSign), .groups = "drop") |>
        filter(.data$nSigns > 1L) |>
        pull(".mvTidx")
    if (length(ambiguous) == 0L) {
        return(h)
    }
    filter(h, !is_in(.data$.mvTidx, ambiguous))
}

# Backwards-compat alias for external callers

#' Parse a region string into its components
#'
#' Split a region identifier (e.g. \code{"chr1:100-200"} or
#' \code{"chr1_100_200"}) into chromosome, start and end. Non-character or
#' non-scalar input is returned unchanged.
#'
#' @param region A single region string.
#' @return A vector/list with the chromosome, start and end, or the input
#'   unchanged when it is not a scalar string.
#' @importFrom stringr str_split
#' @examples
#' parseRegion("chr1:1000000-2000000")
#' @export
parseRegion <- function(region) {
    if (!is.character(region) || length(region) != 1) {
        return(region)
    }

    if (!str_detect(region, "^chr[0-9XY]+:[0-9]+-[0-9]+$")) {
        abort("Input string format must be 'chr:start-end'.")
    }
    parts <- str_split_1(region, "[:-]")
    tibble(
        chrom = canonChrom(parts[1]),
        start = as.integer(parts[2]),
        end = as.integer(parts[3])
    )
}

#' Utility function to convert LD region_ids to `region of interest` dataframe
#'
#' The first field is treated as the chromosome and kept as a normalized string
#' (via \code{canonChrom}) so X/Y/MT survive; the remaining fields are genomic
#' positions and are returned as integers.
#' @param ldRegionId A string of region in the format of chrom_start_end.
#' @param colnames Character vector of length 3 giving output column names for
#'   chromosome, start and end. Default \code{c("chrom", "start", "end")}.
#' @return A tibble with one row per input region and columns named by
#'   \code{colnames}: a normalized character chromosome plus integer start/end.
#' @examples
#' regionToDf(c("1_100_200", "2_300_400"))
#' @export
regionToDf <- function(ldRegionId, colnames = c("chrom", "start", "end")) {
    parts <- str_split(ldRegionId, "[_:-]", simplify = TRUE)
    regionOfInterest <- as_tibble(parts, .name_repair = "minimal")
    colnames(regionOfInterest) <- colnames
    regionOfInterest |>
        mutate(
            across(all_of(colnames[1]), canonChrom),
            across(all_of(colnames[-1]), as.integer)
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
#' @examples
#' regions <- c("chr1:100-200", "chr1:300-400")
#' asGranges(regions = regions)
#' @export
asGranges <- function(regions) {
    if (is.character(regions)) {
        df <- regionToDf(regions)
    } else if (is.data.frame(regions)) {
        if (!all(is_in(c("chrom", "start", "end"), names(regions)))) {
            abort("data.frame must have columns: chrom, start, end")
        }
        df <- regions
    } else {
        msg <- glue(
            "regions must be a character vector or data.frame with ",
            "chrom/start/end columns"
        )
        abort(msg)
    }
    # GRanges expects character seqnames; prefix with "chr" if numeric
    seqnames <- as.character(df$chrom)
    if (!any(str_detect(seqnames, "^chr"))) {
        seqnames <- str_c("chr", seqnames)
    }
    GenomicRanges::GRanges(
        seqnames = seqnames,
        ranges = IRanges::IRanges(
            start = as.integer(df$start),
            end = as.integer(df$end)
        )
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
#' @examples
#' regionsOverlap("chr1:100-200", "chr1:150-250")
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
#' @examples
#' query <- "chr1:100-200"
#' targets <- c("chr1:150-250", "chr1:300-400")
#' findOverlappingRegions(query = query, targets = targets)
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
#' @param ids A character vector of variant IDs in "chr:pos:ref:alt" format, or
#'   a data.frame with A2 (ref) and A1 (alt) columns (e.g., from
#'   \code{\link{parseVariantId}}).
#' @return A character vector with one of "SNP", "insertion", "deletion", or
#'   "MNP" for each variant.
#' @examples
#' classifyVariantType(c("chr1:100:A:G", "chr1:200:AT:A"))
#' @export
classifyVariantType <- function(ids) {
    if (is.character(ids)) {
        ids <- parseVariantId(ids)
    }
    if (!is.data.frame(ids) || !all(is_in(c("A2", "A1"), names(ids)))) {
        msg <- glue(
            "Input must be a character vector of variant IDs or a ",
            "data.frame with A2 and A1 columns."
        )
        abort(msg)
    }
    lenRef <- str_length(ids$A2)
    lenAlt <- str_length(ids$A1)
    isSnpRef <- str_detect(ids$A2, "^[ATCG]$")
    isSnpAlt <- str_detect(ids$A1, "^[ATCG]$")
    case_when(
        lenRef > lenAlt ~ "deletion",
        lenAlt > lenRef ~ "insertion",
        lenRef == 1L & lenAlt == 1L & isSnpRef & isSnpAlt ~ "SNP",
        lenRef == lenAlt ~ "MNP",
        .default = ""
    )
}

# =============================================================================
# Variant identity <-> ranges
# -----------------------------------------------------------------------------
# A variant id is a RENDERING of (seqname, position, REF, ALT), not stored
# state: the identity lives in a one-width GRanges row plus its A2 (REF) / A1
# (ALT) mcols. These two helpers are the boundary conversions -- render on the
# way out, parse on the way in. See spec 4.1a.
# =============================================================================

# Render the variant ids of a GRanges from its coordinates and alleles.
# @noRd
.grVariantIds <- function(gr) {
    if (length(gr) == 0L) {
        return(character(0))
    }
    mc <- mcols(gr)
    formatVariantId(
        as.character(seqnames(gr)),
        start(gr),
        as.character(mc$A2),
        as.character(mc$A1)
    )
}

# Build the element GRanges for a set of variant ids.
#
# Ids that do not encode coordinates are a CONSTRUCTION ERROR. Substituting a
# synthetic range would make range() lie, and range() is the block identity that
# subsetRegion() and the (tuple, range) uniqueness key are both built on -- a
# fabricated coordinate corrupts both silently.
# @noRd
.variantIdsToGRanges <- function(ids, what = "variant ids") {
    ids <- as.character(ids)
    if (length(ids) == 0L) {
        gr <- GenomicRanges::GRanges()
        mcols(gr) <- S4Vectors::DataFrame(
            A1 = character(0),
            A2 = character(0)
        )
        return(gr)
    }
    parsed <- parseVariantId(ids)
    bad <- is.na(parsed$chrom) | is.na(parsed$pos)
    if (any(bad)) {
        msg <- glue(
            "{what}: {sum(bad)} of {length(ids)} do not encode ",
            "coordinates ",
            "(expected chrom:pos:ref:alt), e.g. ",
            "{str_flatten(ids[bad][seq_len(min(3L, sum(bad)))], ', ')}. ",
            "A variant id renders its range and alleles, so an id without ",
            "coordinates has no variant identity to store."
        )
        abort(msg)
    }
    gr <- GenomicRanges::GRanges(
        withChrPrefix(parsed$chrom),
        IRanges::IRanges(start = parsed$pos, width = 1L)
    )
    mcols(gr) <- S4Vectors::DataFrame(
        A1 = as.character(parsed$A1),
        A2 = as.character(parsed$A2)
    )
    gr
}
