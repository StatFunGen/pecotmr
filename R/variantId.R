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

#' Canonicalize chromosome identifiers to a normalized string.
#'
#' Strips a leading "chr"/"ch" prefix (case-insensitive), uppercases, drops
#' the leading zeros of a numeric label (01 -> 1, as UK Biobank's BGEN files
#' write chromosomes) and maps common synonyms to a single form (23 -> X,
#' 24 -> Y, M -> MT). Returns a
#' character vector and never coerces to integer, so non-autosomal chromosomes
#' (X, Y, MT) survive instead of collapsing to NA. This is the single chromosome
#' normalizer used for variant matching and identity throughout the package; it
#' mirrors the chromosome handling previously inlined in summary-statistics QC.
#'
#' @param x Character or numeric vector of chromosome identifiers.
#' @return Character vector of normalized chromosome names.
#' @noRd
canonChrom <- function(x) {
    bare <- as.character(x) |>
        str_remove(regex("^chr", ignore_case = TRUE)) |>
        str_remove(regex("^ch", ignore_case = TRUE)) |>
        str_to_upper() |>
        str_remove("^0+(?=[0-9]+$)")
    canonChromLabel(bare)
}

# Fold the numeric / abbreviated synonyms of the non-autosomes onto their
# canonical labels. NA passes through as NA.
# @noRd
canonChromLabel <- function(x) {
    ok <- !is.na(x)
    x |>
        replace(ok & x == "23", "X") |>
        replace(ok & x == "24", "Y") |>
        replace(ok & x == "M", "MT")
}

#' Chromosome identifiers in the UCSC style.
#'
#' The style \code{seqlevelsStyle(x) <- "UCSC"} gives, and the only one pecotmr
#' holds: every chromosome label that enters -- read from a file, parsed from a
#' variant id, or passed to a constructor or as a range argument -- goes
#' through here once, so comparisons inside the package are exact.
#' \code{"1"}, \code{"chr1"},
#' \code{"CHR1"} and \code{"01"} all become \code{"chr1"}; the numeric codes
#' fold as in \code{canonChrom} (\code{23 -> chrX}, \code{24 -> chrY}); and
#' the mitochondrion, \code{"M"} or \code{"MT"}, becomes \code{"chrM"}. A
#' label outside the standard set keeps its spelling behind the prefix.
#' @param x Character or numeric vector of chromosome identifiers.
#' @return Character vector of UCSC-style chromosome names.
#' @noRd
withChrPrefix <- function(x) {
    # Worked out once per distinct label: a variant table has millions of
    # rows but a few dozen chromosomes.
    labels <- unique(as.character(x))
    bare <- str_remove(labels, regex("^chr", ignore_case = TRUE))
    canon <- canonChrom(bare)
    standard <- is_in(canon, .STANDARD_CHROMS)
    out <- str_c(
        "chr",
        if_else(standard, str_replace(canon, "^MT$", "M"), bare)
    )
    out[match(as.character(x), labels)]
}

# The chromosomes withChrPrefix() spells as UCSC does, in canonChrom() form.
# @noRd
.STANDARD_CHROMS <- c(as.character(seq_len(22L)), "X", "Y", "MT")

# `gr` with its chromosomes in UCSC style, for ranges entering from outside:
# a constructor's ranges, a range argument. Strand, metadata columns and the
# genome build are kept, and ranges already in UCSC style come back as they
# are. Anything but ranges (NULL, a stray value a constructor is about to
# refuse) passes through untouched.
# @noRd
.ucscRanges <- function(gr) {
    if (!methods::is(gr, "GenomicRanges") && !methods::is(gr, "GRangesList")) {
        return(gr)
    }
    old <- GenomeInfoDb::seqlevels(gr)
    new <- withChrPrefix(old)
    if (identical(old, new)) {
        return(gr)
    }
    if (!anyDuplicated(new)) {
        return(GenomeInfoDb::renameSeqlevels(gr, set_names(new, old)))
    }
    # Two spellings of one chromosome ("1" and "chr1") cannot both be
    # renamed onto one name, so the ranges are rebuilt on the merged names.
    out <- S4Vectors::`mcols<-`(
        GenomicRanges::GRanges(
            withChrPrefix(as.character(GenomicRanges::seqnames(gr))),
            IRanges::ranges(gr),
            strand = GenomicRanges::strand(gr)
        ),
        value = S4Vectors::mcols(gr)
    )
    .withGenome(out, .oneBuild(gr))
}

# A collection constructor's `entry` list with every element's chromosomes in
# UCSC style; a NULL element stays NULL.
# @noRd
.ucscEntries <- function(entry) {
    map(as.list(entry), .ucscRanges)
}

# A context or genotype experiment with its rowRanges in UCSC style, so a
# dataset compares them exactly with every other range pecotmr holds.
# @noRd
.ucscExperiment <- function(se) {
    if (!methods::is(se, "RangedSummarizedExperiment")) {
        return(se)
    }
    old <- GenomeInfoDb::seqlevels(se)
    if (identical(withChrPrefix(old), old)) {
        return(se)
    }
    SummarizedExperiment::`rowRanges<-`(
        se,
        value = .ucscRanges(SummarizedExperiment::rowRanges(se))
    )
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

#' Strip build suffix from variant IDs (e.g., ":b38" or "_b38").
#' @param x Character vector of variant IDs.
#' @return Character vector with build suffix removed.
#' @noRd
stripBuildSuffix <- function(x) {
    # The regex only where a suffix could be: nearly no id carries one, and
    # a fixed-string scan costs a fraction of a regex over millions of ids.
    candidate <- which(
        str_detect(x, fixed(":b")) | str_detect(x, fixed("_b"))
    )
    if (length(candidate) == 0L) {
        return(x)
    }
    replace(x, candidate, str_remove(x[candidate], "(:|_)b[0-9]+$"))
}

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

#' Re-render variant ids whose allele slot holds a tag instead of an allele.
#'
#' Some LD panels (e.g. ADSP R5) name a variant with a tag where an allele
#' belongs -- \code{chr21:13988152:INS:T} whose real alleles are
#' \code{A}/\code{AT} -- so the id string stops encoding its alleles even
#' though the panel's own A1/A2 columns are correct. The id-string matchers
#' (\code{.ldFromRefMatch}, \code{.subsetLdRefToIds}) key on the id and
#' cannot interpret such a record. The allele COLUMNS are authoritative, so the
#' id is rebuilt from them; the variant is never dropped for how its id happens
#' to spell an allele.
#'
#' Only an id that parses to \code{chr:pos:allele:allele} but whose allele
#' slots are not pure DNA is rewritten, and only where both supplied alleles are
#' present. An id that does not parse to the allele form -- an rsID, say -- is a
#' different, valid convention that \code{matchVariants} handles by exact
#' string, so it is left alone; so is a clean allele-encoded id, which keeps the
#' panel's own spelling. The DNA test is case-insensitive, so a lower-case id
#' (\code{chr1:100:a:g}) is clean, not a tag. Rewritten ids follow the
#' convention (chr prefix, field separator) detected on \code{ids}, so a
#' repaired id looks like its neighbours.
#'
#' @param ids Character vector of variant ids (\code{chr:pos:A2:A1}).
#' @param A2,A1 Character vectors of authoritative reference / alternate
#'   alleles, parallel to \code{ids}.
#' @return \code{ids} with the tag-allele entries re-rendered from
#'   \code{A2}/\code{A1}.
#' @noRd
.repairVariantIds <- function(ids, A2, A1) {
    ids <- as.character(ids)
    if (length(A2) != length(ids) || length(A1) != length(ids)) {
        return(ids)
    }
    p <- parseVariantId(ids)
    dna <- regex("^[ACGT]+$", ignore_case = TRUE)
    tagged <- !is.na(p$A1) &
        !is.na(p$A2) &
        !replace_na(str_detect(p$A1, dna) & str_detect(p$A2, dna), FALSE)
    repair <- tagged & !is.na(A2) & !is.na(A1) & A2 != "" & A1 != ""
    if (!any(repair)) {
        return(ids)
    }
    replace(
        ids,
        repair,
        formatVariantId(
            p$chrom[repair],
            p$pos[repair],
            A2[repair],
            A1[repair],
            convention = attr(p, "convention")
        )
    )
}

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

#' Parse variant IDs into a data frame
#'
#' Converts variant IDs from any supported string format or data.frame into a
#' standardized data.frame with a character chrom, integer pos, and character
#' allele columns (A2, A1). The chrom is spelled in UCSC style, as every
#' chromosome pecotmr holds is ("chr1", "chrX", "chrM"; 23 -> chrX,
#' 24 -> chrY). Supports colon-separated ("chr1:100:A:G"),
#' underscore-separated ("1_100_A_G"), with or without "chr" prefix, and with
#' optional build suffix (":b38" or "_b38"). The detected input convention is
#' stored as an attribute.
#'
#' @param ids A character vector of variant IDs, or a data.frame with columns
#'   "chrom", "pos", and allele columns (A2/A1 or ref/alt or any 4-column
#'   layout).
#' @return A tibble with columns "chrom" (character, UCSC style), "pos"
#'   (integer), "A2" (character), "A1" (character). The detected convention is
#'   stored as \code{attr(result, "convention")}.
#' @examples
#' parseVariantId(c("chr1:100:A:G", "chr2:200:T:C"))
#' @importFrom checkmate assert checkCharacter checkDataFrame
#' @export
parseVariantId <- function(ids) {
    # `ids` is documented as a character vector OR a data.frame.
    assert(
        checkCharacter(ids),
        checkDataFrame(ids),
        .var.name = "ids"
    )
    if (is.data.frame(ids)) {
        return(.parseVariantIdDf(ids))
    }
    convention <- detectVariantConvention(ids)
    m <- .variantIdFields(ids)
    data <- tibble(
        chrom = withChrPrefix(m[, 2L]),
        pos = as.integer(m[, 3L]),
        A2 = m[, 4L],
        A1 = m[, 5L]
    )
    `attr<-`(data, "convention", convention)
}

# The fields of variant ids as every parser reads them: underscores read as
# colons and a build suffix dropped, then exactly 4 fields split by a single
# vectorized regex match. A non-matching id yields an all-NA capture row
# (mirrors strcapture).
# @noRd
.variantIdFields <- function(ids) {
    normalized <- stripBuildSuffix(str_replace_all(ids, fixed("_"), ":"))
    str_match(normalized, "^([^:]+):([^:]+):([^:]+):([^:]+)")
}

# `[[<-` rather than mutate(): callers hand this frame straight through from
# a user-supplied table, which may carry NA or "" column names -- and
# dplyr refuses to transform such a frame at all.
# @noRd
.withCanonicalCoords <- function(ids) {
    # `value =` named: tibble's `[[<-` method carries an extra `j` formal, so a
    # third POSITIONAL argument binds to `j` and the write silently misfires.
    coerced <- `[[<-`(ids, "chrom", value = withChrPrefix(ids$chrom))
    `[[<-`(coerced, "pos", value = as.integer(ids$pos))
}

# Parse a data.frame of already-split ids: resolve the 4 identity columns
# (explicit chrom/pos/A2/A1 or chrom/pos/A1/A2 kept as-is; otherwise positional
# chrom/pos/A2/A1), then canonicalize chrom/pos.
# @noRd
.parseVariantIdDf <- function(ids) {
    # minimal repair: preserve any empty/duplicate extra-column names (e.g. an
    # unnamed passthrough column) for .sanitizeNames() to canonicalize later.
    tbl <- as_tibble(ids, .name_repair = "minimal")
    hasA2A1 <- all(is_in(c("chrom", "pos", "A2", "A1"), names(tbl)))
    hasA1A2 <- all(is_in(c("chrom", "pos", "A1", "A2"), names(tbl)))
    # With neither allele naming present the frame is positional: its first
    # four columns are the variant key.
    named <- if (!hasA2A1 && !hasA1A2 && ncol(tbl) >= 4) {
        `names<-`(
            tbl,
            replace(names(tbl), seq_len(4), c("chrom", "pos", "A2", "A1"))
        )
    } else {
        tbl
    }
    conv <- list(
        hasChr = any(str_detect(as.character(named$chrom), "^chr")),
        alleleSep = ":",
        hasBuild = FALSE,
        example = NA_character_
    )
    `attr<-`(.withCanonicalCoords(named), "convention", conv)
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
#' @param chrPrefix Logical. \code{TRUE} (default) spells the chromosome in UCSC
#'   style (\code{chr1}, \code{chrX}, \code{chrM}); \code{FALSE} writes it
#'   bare (\code{1}, \code{X}, \code{MT}). Ignored if \code{convention} is
#'   provided.
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
    # UCSC style is pecotmr's own; the bare form is for a consumer that wants
    # it (canonChrom keeps X/Y/MT as strings).
    # Over the distinct labels only: millions of variants, a few dozen
    # chromosomes (withChrPrefix() does the same).
    chromOut <- if (chrPrefix) {
        withChrPrefix(chrom)
    } else {
        labels <- unique(chrom)
        canonChrom(labels)[match(chrom, labels)]
    }
    str_c(chromOut, ":", pos, alleleSep, A2, alleleSep, A1)
}

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
#' @param chrPrefix Logical. \code{TRUE} (default) writes the chromosome in
#'   UCSC style (\code{chr1}, \code{chrM}), the style pecotmr uses
#'   throughout; \code{FALSE} writes it bare. Ignored if \code{convention} is
#'   provided.
#' @param convention Optional list from \code{detectVariantConvention} or
#'   \code{attr(parseVariantId(ids), "convention")}. When provided, the output
#'   format is driven automatically by the detected convention.
#' @return A character vector of re-formatted variant IDs. Unparseable ids (e.g.
#'   rsIDs) are returned unchanged.
#' @examples
#' normalizeVariantId(c("1:100:A:G", "2:200:T:C"))
#' @importFrom checkmate assertCharacter assertFlag assertList
#' @export
normalizeVariantId <- function(ids, chrPrefix = TRUE, convention = NULL) {
    assertCharacter(ids)
    assertFlag(chrPrefix)
    assertList(convention, null.ok = TRUE)
    out <- as.character(ids)
    # An id already in the canonical form re-formats to itself, and a
    # genotype panel's millions of ids usually all are: one regex finds them,
    # where parsing and re-rendering every one costs several times as much.
    todo <- if (chrPrefix && is.null(convention)) {
        which(!str_detect(out, .CANONICAL_VARIANT_ID))
    } else {
        seq_along(out)
    }
    if (length(todo) == 0L) {
        return(out)
    }
    parsed <- parseVariantId(out[todo])
    # Only re-format ids that parsed into a chrom + pos; leave unparseable ids
    # (e.g. rsIDs) unchanged rather than emitting "chrNA:..." garbage.
    ok <- !is.na(parsed$chrom) & !is.na(parsed$pos)
    if (!any(ok)) {
        return(out)
    }
    reformatted <- if (!is.null(convention)) {
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
    replace(out, todo[ok], reformatted)
}

# The ids normalizeVariantId() returns unchanged in its default convention: a
# UCSC chromosome, a position without leading zeros, and two alleles, with
# no build suffix and no underscore.
# @noRd
.CANONICAL_VARIANT_ID <- "^chr(1[0-9]|2[0-2]|[1-9]|X|Y|M):[1-9][0-9]*:[^:_]+:[^:_]+$"

# Complement a DNA allele string (A<->T, C<->G) for strand flipping.
# @noRd
.strandFlip <- function(ref) chartr("ATCG", "TAGC", ref)

# Ensure a data.frame has unique, non-empty column names (blanks become
# `unnamed_<i>`, duplicates de-duplicated with make.unique).
# @noRd
.sanitizeNames <- function(df) {
    nm <- colnames(df) %||% rep("unnamed", ncol(df))
    emptyIdx <- is.na(nm) | nm == ""
    named <- replace(nm, emptyIdx, str_c("unnamed_", seq_len(sum(emptyIdx))))
    `colnames<-`(df, make.unique(named, sep = "_"))
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
#' @param removeStrandAmbiguous Whether to drop a strand-ambiguous match (A/T
#'   or C/G alleles, where a strand flip and an allele swap look the same).
#'   Each target-reference pair is judged on its own, whatever else is being
#'   harmonized; an exact A/T match is kept. Default is `TRUE`.
#' @param flipStrand Whether to output the variants after strand flip. Default
#'   is `FALSE`.
#' @param removeUnmatched Whether to remove unmatched variants. Default is
#'   `TRUE`.
#' @return A list with \code{harmonizedData}, the post-QC variant
#'   data.frame, and \code{qcSummary}, the per-variant merge/flip/strand
#'   diagnostics.
#' @importFrom dplyr mutate inner_join filter pull select everything row_number
#' @importFrom dplyr if_else any_of all_of rename across
#' @importFrom purrr partial
#' @importFrom vctrs vec_duplicate_detect
#' @importFrom tidyr separate
#' @keywords internal
#' @details Pure panel-vs-sumstats allele harmonization: match by (chrom, pos),
#'   detect A1/A2 swap, sign-flip \code{colToFlip} columns and complement
#'   \code{colToComplement} columns on swap. Variant-allele filters (indels,
#'   strand-ambiguous, duplicates) are applied here directly when the
#'   corresponding \code{removeIndels} / \code{removeStrandAmbiguous} /
#'   \code{removeDups} flags are set; MAF / INFO / N column-numeric filters run
#'   in \code{.applyContentFilters()} before this function.
#' @noRd
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
    colToComplement = character()
) {
    coerced <- .harmonizeCoerceInputs(targetData, refVariants)
    target <- coerced$targetData
    panel <- coerced$refVariants
    # The index rides along on the join copy only: the restore path below
    # re-reads `target`, and a sentinel column there would bind_rows its
    # way into the returned frame.
    joined <- .harmonizeJoin(.haIndexTargets(target), panel)
    if (nrow(joined) == 0) {
        return(.harmonizeEmptyResult(joined))
    }
    matchResult <- .harmonizeDecideMatches(
        joined,
        colToFlip,
        colToComplement,
        flipStrand,
        removeIndels,
        removeStrandAmbiguous
    )
    qcCounts <- .harmonizeQcCounts(matchResult)
    cleaned <- .harmonizeCleanResult(matchResult)
    deduped <- if (removeDups) .harmonizeRemoveDups(cleaned) else cleaned
    # Putting the unmatched rows back also re-derives the QC summary over
    # them, so the two travel together.
    final <- if (removeUnmatched) {
        list(result = deduped, qcSummary = matchResult)
    } else {
        .harmonizeRestoreUnmatched(deduped, matchResult, target)
    }
    .harmonizeFinalChecks(final$result, panel, matchMinProp)
    out <- list(harmonizedData = final$result, qcSummary = final$qcSummary)
    `attr<-`(out, "qcCounts", qcCounts)
}

# Turn the raw (chrom, pos) join into a decided match table: classify each
# pair's allele relationship, drop the strand-ambiguous and indel rows the
# caller refused, settle each target on a single reference row, then apply the
# sign / frequency flips the chosen orientation implies. Everything before this
# builds the join; everything after it shapes the output.
# @noRd
.harmonizeDecideMatches <- function(
    matchResult,
    colToFlip,
    colToComplement,
    flipStrand,
    removeIndels,
    removeStrandAmbiguous
) {
    decided <- .harmonizeFlags(matchResult) |>
        .harmonizeResolveAmbiguity(removeStrandAmbiguous) |>
        .harmonizeKeepRule(removeIndels) |>
        .harmonizeResolveTargets()
    .harmonizeApplyFlips(
        decided,
        colToFlip,
        colToComplement,
        flipStrand
    )
}

# Target-row index carried through the (chrom, pos) join so the resolution
# below can group the reference rows one target variant matched.
.haTargetIdx <- ".haTidx"

# Assigned directly rather than through mutate(): the target frame is not
# name-repaired until the join (.sanitizeNames), and dplyr's data mask refuses
# a frame carrying an NA or empty column name.
# @noRd
.haIndexTargets <- function(targetData) {
    `[[<-`(targetData, .haTargetIdx, value = seq_len(nrow(targetData)))
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
    "keep",
    .haTargetIdx
)

# Coerce both sides to canonical variant data.frames and strip merge-conflicting
# columns (variant_id is dropped on both sides because it is rebuilt from the
# QC'd alleles below; leaving it in collides on the final rename).
# @noRd
.harmonizeCoerceInputs <- function(targetData, refVariants) {
    parsedTarget <- if (
        is.data.frame(targetData) &&
            ncol(targetData) > 4 &&
            all(is_in(c("chrom", "pos", "A2", "A1"), names(targetData)))
    ) {
        variantCols <- c("chrom", "pos", "A2", "A1")
        variantDf <- targetData |> select(all_of(variantCols))
        otherCols <- targetData |> select(-all_of(variantCols))
        bind_cols(
            parseVariantId(variantDf),
            otherCols,
            .name_repair = "minimal"
        )
    } else {
        parseVariantId(targetData)
    }
    parsedRef <- parseVariantId(refVariants)
    dropCols <- c("chromosome", "position", "ref", "alt", "variant_id")
    # select() rejects a frame with an empty column name, which a passthrough
    # column legitimately has, so only touch a frame that has something to
    # drop.
    trimmedTarget <- if (any(is_in(dropCols, colnames(parsedTarget)))) {
        select(parsedTarget, -any_of(dropCols))
    } else {
        parsedTarget
    }
    trimmedRef <- if (is_in("variant_id", colnames(parsedRef))) {
        select(parsedRef, -any_of("variant_id"))
    } else {
        parsedRef
    }
    list(targetData = trimmedTarget, refVariants = trimmedRef)
}

# Target + reference rows paired at the same (chrom, pos) by the matching
# engine (.variantPairs()), laid out as an inner join on (chrom, pos) lays
# them out: the target's columns, then the reference's other columns, the ones
# both carry suffixed ".target" / ".ref".
# @noRd
#' @importFrom dplyr bind_cols
.harmonizeJoin <- function(targetData, refVariants) {
    pairs <- .variantPairs(
        .frameVariantTable(targetData),
        .frameVariantTable(refVariants)
    )
    keys <- c("chrom", "pos")
    left <- targetData[pairs$query, , drop = FALSE]
    right <- refVariants[pairs$subject, , drop = FALSE]
    right <- right[, !is_in(colnames(right), keys), drop = FALSE]
    both <- setdiff(intersect(colnames(left), colnames(right)), keys)
    colnames(left) <- .suffixed(colnames(left), both, ".target")
    colnames(right) <- .suffixed(colnames(right), both, ".ref")
    bind_cols(left, right, .name_repair = "minimal") |>
        as_tibble(.name_repair = "minimal") |>
        .sanitizeNames()
}

# `nm` with `suffix` appended to the names in `which`.
# @noRd
.suffixed <- function(nm, which, suffix) {
    hit <- !is.na(nm) & is_in(nm, which)
    replace(nm, hit, str_c(nm[hit], suffix))
}

# The variant table (.variantTable()) of a parsed chrom/pos/A2/A1 frame.
# @noRd
.frameVariantTable <- function(df) {
    list(
        chrom = as.character(df$chrom),
        pos = as.integer(df$pos),
        A2 = as.character(df$A2),
        A1 = as.character(df$A1)
    )
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
    `attr<-`(
        emptyOut,
        "qcCounts",
        list(
            considered = 0L,
            signFlip = 0L,
            strandFlip = 0L,
            kept = 0L,
            dropped = 0L,
            droppedIndel = 0L,
            droppedAmbiguous = 0L,
            droppedOther = 0L
        )
    )
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

# Exact / sign-flip / strand-flip / INDEL / ID-match indicators, as the
# matching engine decides them (.variantPairFlags()), on the uppercased
# alleles from the earlier pipe stages.
# @noRd
.harmonizeMatchFlags <- function(matchResult) {
    f <- .variantPairFlags(
        matchResult$A1.target,
        matchResult$A2.target,
        matchResult$A1.ref,
        matchResult$A2.ref
    )
    mutate(
        matchResult,
        strand_unambiguous = f$strandUnambiguous,
        exact_match = f$exact,
        sign_flip = f$signFlip,
        strand_flip = f$strandFlip,
        INDEL = f$indel,
        ID_match = f$idMatch
    )
}

# Strand-ambiguity resolution: the A/T-C/G guard is off when the caller opted
# out. Otherwise each pair is judged on its own: an ambiguous pair is not
# rescued by the rest of the batch having no strand flip.
# @noRd
.harmonizeResolveAmbiguity <- function(matchResult, removeStrandAmbiguous) {
    if (removeStrandAmbiguous) {
        return(matchResult)
    }
    mutate(matchResult, strand_unambiguous = TRUE)
}

# The keep flag, as the matching engine decides it (.variantPairKeep()); the
# strand guard is already in `strand_unambiguous`.
# @noRd
.harmonizeKeepRule <- function(matchResult, removeIndels) {
    keep <- .variantPairKeep(
        list(
            strandUnambiguous = matchResult$strand_unambiguous,
            exact = matchResult$exact_match,
            signFlip = matchResult$sign_flip,
            strandFlip = matchResult$strand_flip,
            indel = matchResult$INDEL,
            idMatch = matchResult$ID_match
        ),
        removeStrandAmbiguous = TRUE,
        removeIndels = removeIndels
    )
    mutate(matchResult, keep = keep)
}

# Reduce the join to at most one reference row per TARGET variant.
#
# The join is on (chrom, pos), so a target matches every reference row at its
# position whose alleles reconcile. Usually that is one row. It is more when
# the reference carries the variant AND its own flip as separate entries --
# two distinct indels at one position that happen to be each other's flip is
# rare but biologically real -- or the same variant once per strand.
#
# Matches that DISAGREE about the allele-swap sign are undecidable: the target
# is either the exact entry or the swapped one and nothing says which. Keeping
# the first would make the effect direction depend on the reference's row
# order; keeping both would return the target twice with opposite signs, which
# is not harmonization but duplication. Those targets are dropped. Matches that
# agree (the same variant listed twice, or once per strand) resolve to the
# first row, so this costs nothing where there is no ambiguity.
#
# Deciding it here, before the keep flags are counted, is what lets every
# caller -- matchVariants, the summaryStatsQc panel harmonization, the cTWAS
# weight harmonization -- agree on the surviving variant set without repeating
# the rule.
# @noRd
.harmonizeResolveTargets <- function(matchResult) {
    # Named apart from the column it replaces: `keep = keep` inside mutate()
    # reads the column, not this vector, and silently changes nothing.
    resolvedKeep <- .variantPairResolve(
        matchResult[[.haTargetIdx]],
        matchResult$keep,
        matchResult$sign_flip
    )
    mutate(matchResult, keep = resolvedKeep)
}

# Named per-row conditional column transforms for the harmonize across() calls.
# The `flip` condition vector is bound with purrr::partial() rather than passed
# through across()'s `...`, which dplyr deprecated in 1.1.0 -- partial() fixes
# the argument without an anonymous function, so the no-lambda rule the `...`
# form was chosen for still holds. partial() is LAZY, so callers bind the
# condition to a local first: what it captures must not be reassigned before
# across() calls it.
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
    signFlip <- matchResult$sign_flip
    negated <- if (is.null(colToFlip)) {
        matchResult
    } else {
        .harmonizeCheckCols(colToFlip, matchResult)
        matchResult |>
            mutate(across(
                all_of(colToFlip),
                partial(.negateWhere, flip = signFlip)
            ))
    }
    complemented <- if (length(colToComplement) == 0L) {
        negated
    } else {
        .harmonizeCheckCols(colToComplement, negated)
        negated |>
            mutate(across(
                all_of(colToComplement),
                partial(.complementWhere, flip = signFlip)
            ))
    }
    if (!flipStrand) {
        return(complemented)
    }
    .harmonizeFlipStrandCols(complemented)
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
    strandFlip <- matchResult$strand_flip
    matchResult |>
        mutate(across(
            c("A1.target", "A2.target"),
            partial(.strandFlipWhere, flip = strandFlip)
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
    droppedAmbiguous <- sum(
        !matchResult$keep &
            matchResult$strand_flip &
            !matchResult$strand_unambiguous &
            if (hasIndel) !matchResult$INDEL else TRUE,
        na.rm = TRUE
    )
    list_assign(
        qcCounts,
        droppedAmbiguous = droppedAmbiguous,
        droppedOther = qcCounts$dropped -
            qcCounts$droppedIndel -
            droppedAmbiguous
    )
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
    withUnmatched <- bind_rows(
        result,
        unmatchData |> mutate(variants_id_original = .data$variant_id)
    )
    # Restore the caller's input order, which bind_rows put the unmatched
    # rows behind.
    restored <- withUnmatched |>
        slice(match(targetData$variant_id, .data$variants_id_original)) |>
        select(-any_of("variants_id_original"))
    list(result = restored, qcSummary = qcSummary)
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

# =============================================================================
# The variant matching engine
# -----------------------------------------------------------------------------
# Variants are paired by position first -- overlaps of width-1 ranges, so a
# window (`maxDistance`) costs nothing extra -- and their alleles decided
# second, pair by pair: exact, ref/alt swap (sign -1), strand flip, the I/D
# codes. matchVariants() and harmonizeAlleles() both decide through it, so the
# two cannot disagree about which variants match.
# =============================================================================

# A variant table: canonical chromosome labels, positions and both alleles,
# from a GRanges (A1/A2 metadata columns), variant ids, or a chrom/pos/A2/A1
# data.frame (parseVariantId()). Chromosomes are relabelled over their distinct
# values only (withChrPrefix()), so a genome-wide GRanges costs one pass.
# @noRd
.variantTable <- function(x) {
    if (methods::is(x, "GRanges")) {
        mc <- mcols(x, use.names = FALSE)
        if (!all(is_in(c("A1", "A2"), colnames(mc)))) {
            abort("A GRanges of variants needs A1 and A2 metadata columns.")
        }
        labels <- withChrPrefix(as.character(GenomeInfoDb::seqlevels(x)))
        return(list(
            chrom = labels[as.integer(seqnames(x))],
            pos = start(x),
            A2 = as.character(mc$A2),
            A1 = as.character(mc$A1)
        ))
    }
    p <- parseVariantId(x)
    list(
        chrom = as.character(p$chrom),
        pos = as.integer(p$pos),
        A2 = as.character(p$A2),
        A1 = as.character(p$A1)
    )
}

# The pairs of variants of tables `a` and `b` on one chromosome at most
# `maxDistance` bp apart, as list(query, subject) indices, in query then
# subject order -- the order a join on (chrom, pos) gives. A variant with no
# chromosome or position pairs with none.
# @noRd
#' @importFrom IRanges findOverlaps IRanges
.variantPairs <- function(a, b, maxDistance = 0L) {
    levels <- unique(c(a$chrom, b$chrom))
    levels <- levels[!is.na(levels)]
    ra <- .variantPairRanges(a, levels)
    rb <- .variantPairRanges(b, levels)
    hits <- findOverlaps(
        ra$gr,
        rb$gr,
        maxgap = as.integer(maxDistance) - 1L,
        ignore.strand = TRUE
    )
    q <- ra$rows[S4Vectors::queryHits(hits)]
    sbj <- rb$rows[S4Vectors::subjectHits(hits)]
    o <- order(q, sbj)
    list(query = q[o], subject = sbj[o])
}

# @noRd
.variantPairRanges <- function(x, levels) {
    rows <- which(!is.na(x$chrom) & !is.na(x$pos))
    gr <- GenomicRanges::GRanges(
        factor(x$chrom[rows], levels = levels),
        IRanges(x$pos[rows], width = 1L)
    )
    list(gr = gr, rows = rows)
}

# `x` upper-cased, over its distinct values (alleles are few).
# @noRd
.upperAlleles <- function(x) {
    u <- unique(x)
    str_to_upper(u)[match(x, u)]
}

# The allele relationship of each pair of a target variant (alleles `a1t`,
# `a2t`) with a reference variant (`a1r`, `a2r`): exact; a ref/alt swap, or a
# swap of the complements (`signFlip`); a strand flip; an indel; an I/D-coded
# target meeting a reference indel (`idMatch`); and whether the target's
# alleles leave the strand undecidable (A/T, C/G). Case-insensitive.
# @noRd
.variantPairFlags <- function(a1t, a2t, a1r, a2r) {
    a1t <- .upperAlleles(a1t)
    a2t <- .upperAlleles(a2t)
    a1r <- .upperAlleles(a1r)
    a2r <- .upperAlleles(a2r)
    flip1 <- .strandFlip(a1r)
    flip2 <- .strandFlip(a2r)
    differ <- a1t != a1r & a2t != a2r
    list(
        strandUnambiguous = !((a1t == "A" & a2t == "T") |
            (a1t == "T" & a2t == "A") |
            (a1t == "C" & a2t == "G") |
            (a1t == "G" & a2t == "C")),
        exact = a1t == a1r & a2t == a2r,
        signFlip = ((a1t == a2r & a2t == a1r) |
            (a1t == flip2 & a2t == flip1)) &
            differ,
        strandFlip = ((a1t == flip1 & a2t == flip2) |
            (a1t == flip2 & a2t == flip1)) &
            differ,
        indel = a2t == "I" |
            a2t == "D" |
            str_length(a2t) > 1L |
            str_length(a1t) > 1L,
        idMatch = (a2t == "D" | a2t == "I") &
            (str_length(a1r) > 1L | str_length(a2r) > 1L)
    )
}

# Whether each pair matches: a strand flip only when the target's alleles
# decide the strand (each pair judged on its own), an exact pair, or an I/D
# match; any other pair when exact, swapped or I/D-matched. Indels are refused
# when `removeIndels`. NA where an allele is missing.
# @noRd
.variantPairKeep <- function(flags, removeStrandAmbiguous, removeIndels) {
    decided <- if (removeStrandAmbiguous) flags$strandUnambiguous else TRUE
    keep <- ifelse(
        flags$strandFlip,
        decided | flags$exact | flags$idMatch,
        flags$exact | flags$signFlip | flags$idMatch
    )
    if (removeIndels) {
        keep <- ifelse(flags$indel, FALSE, keep)
    }
    keep
}

# `keep` with each target settled on at most one reference variant. A target
# whose kept pairs disagree about the swap sign is undecidable and dropped --
# its effect direction would depend on the reference's row order --; one whose
# pairs agree keeps the first. Unchanged when no target has two.
# @noRd
#' @importFrom vctrs vec_unique_loc
.variantPairResolve <- function(query, keep, signFlip) {
    hit <- which(keep)
    target <- query[hit]
    if (anyDuplicated(target) == 0L) {
        return(keep)
    }
    distinctSigns <- vec_unique_loc(data.frame(
        target = target,
        sign = signFlip[hit]
    ))
    nSigns <- tabulate(target[distinctSigns], nbins = max(target))
    resolved <- hit[nSigns[target] == 1L & !duplicated(target)]
    replace(rep(FALSE, length(keep)), resolved, TRUE)
}

# The matched (a, b) pairs of variant tables `a` and `b`: paired within
# `maxDistance`, decided by their alleles, each target settled on one
# reference variant, and each reference variant (its alleles as given)
# claimed by the first target only. list(idxA, idxB, sign), by idxA.
# @noRd
.matchVariantTables <- function(
    a,
    b,
    maxDistance = 0L,
    removeStrandAmbiguous = TRUE,
    removeIndels = FALSE
) {
    pairs <- .variantPairs(a, b, maxDistance)
    q <- pairs$query
    sbj <- pairs$subject
    if (length(q) == 0L) {
        return(.matchVariantsEmpty())
    }
    flags <- .variantPairFlags(a$A1[q], a$A2[q], b$A1[sbj], b$A2[sbj])
    keep <- .variantPairKeep(flags, removeStrandAmbiguous, removeIndels)
    kept <- which(.variantPairResolve(q, keep, flags$signFlip) %in% TRUE)
    claimed <- kept[vec_unique_loc(data.frame(
        chrom = b$chrom[sbj[kept]],
        pos = b$pos[sbj[kept]],
        A2 = b$A2[sbj[kept]],
        A1 = b$A1[sbj[kept]]
    ))]
    list(
        idxA = q[claimed],
        idxB = sbj[claimed],
        sign = as.numeric(ifelse(flags$signFlip[claimed], -1, 1))
    )
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
#' are dropped when \code{removeStrandAmbiguous = TRUE}, each pair judged on
#' its own.
#'
#' When ids cannot be parsed into chrom/pos/ref/alt (e.g. rsIDs), the matcher
#' falls back to exact string identity with \code{sign = 1} (no allele awareness
#' is possible without alleles).
#'
#' @param idsA,idsB Character vectors of variant IDs, data.frames with
#'   chrom/pos/A2/A1 columns (A2 = ref, A1 = alt), or GRanges of variants with
#'   A1/A2 metadata columns.
#' @param allowFlip When TRUE (default) match by (chrom, pos) with ref/alt swap
#'   and strand handling, reporting \code{sign = -1} for a swap. When FALSE,
#'   match on exact alleles only (canonicalized format; \code{sign} always +1)
#'   so a ref/alt swap does not match.
#' @param removeStrandAmbiguous Drop strand-ambiguous (A/T, C/G) matches.
#'   Default TRUE (mirrors summaryStatsQc).
#' @param removeIndels Drop indels. Default FALSE.
#' @param maxDistance Largest distance in bp between the positions of two
#'   variants that may match; their alleles must still agree. Default 0: the
#'   same position.
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
    removeIndels = FALSE,
    maxDistance = 0L
) {
    ranged <- methods::is(idsA, "GRanges") || methods::is(idsB, "GRanges")
    if (!allowFlip && !ranged && maxDistance == 0L) {
        # Exact-allele matching: canonicalize the id FORMAT (chr-prefix +
        # separator, rsID-safe) but keep allele order, then match by exact
        # string identity, so a ref/alt swap does NOT match (sign is always +1).
        return(.matchVariantsExact(
            .matchVariantKeyOf(idsA),
            .matchVariantKeyOf(idsB)
        ))
    }
    tabA <- .variantTable(idsA)
    tabB <- .variantTable(idsB)
    if (!ranged && .matchVariantsUnparseable(tabA, tabB)) {
        # rsID / unparseable ids: exact string identity only (no allele
        # awareness). data.frame inputs have no string form to fall back to.
        if (is.data.frame(idsA) || is.data.frame(idsB)) {
            return(.matchVariantsEmpty())
        }
        return(.matchVariantsExact(as.character(idsA), as.character(idsB)))
    }
    if (!allowFlip) {
        return(.matchVariantTablesExact(tabA, tabB, maxDistance))
    }
    .matchVariantTables(
        tabA,
        tabB,
        maxDistance,
        removeStrandAmbiguous,
        removeIndels
    )
}

# Pairs with exactly the same alleles in the same order (sign +1): the first
# reference variant for each target, within `maxDistance`.
# @noRd
.matchVariantTablesExact <- function(a, b, maxDistance = 0L) {
    pairs <- .variantPairs(a, b, maxDistance)
    q <- pairs$query
    sbj <- pairs$subject
    same <- which(a$A2[q] == b$A2[sbj] & a$A1[q] == b$A1[sbj])
    first <- same[!duplicated(q[same])]
    list(
        idxA = q[first],
        idxB = sbj[first],
        sign = rep(1, length(first))
    )
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
    length(dfA$pos) == 0L ||
        length(dfB$pos) == 0L ||
        anyNA(dfA$chrom) ||
        anyNA(dfA$pos) ||
        anyNA(dfB$chrom) ||
        anyNA(dfB$pos)
}

# The variant table (.variantTable()) of variant ids, read as parseVariantId()
# reads them (NA fields for an id that does not parse). A canonical id
# ("chr1:100:A:G") is split on its separator, at a fraction of the cost of the
# regex every other spelling takes: a genome-wide panel holds millions.
# @noRd
#' @importFrom stringi stri_split_fixed
.idsVariantTable <- function(ids) {
    ids <- as.character(ids)
    canon <- str_detect(ids, .CANONICAL_VARIANT_ID)
    canon[is.na(canon)] <- FALSE
    if (all(canon)) {
        return(.canonicalVariantTable(ids))
    }
    out <- list(
        chrom = rep(NA_character_, length(ids)),
        pos = rep(NA_integer_, length(ids)),
        A2 = rep(NA_character_, length(ids)),
        A1 = rep(NA_character_, length(ids))
    )
    if (any(canon)) {
        fields <- .canonicalVariantTable(ids[canon])
        for (nm in names(out)) {
            out[[nm]][canon] <- fields[[nm]]
        }
    }
    p <- parseVariantId(ids[!canon])
    out$chrom[!canon] <- as.character(p$chrom)
    out$pos[!canon] <- as.integer(p$pos)
    out$A2[!canon] <- as.character(p$A2)
    out$A1[!canon] <- as.character(p$A1)
    .unparsedAsSpelled(out, ids)
}

# A variant table's rows whose id did not parse (an rsID), each made a variant
# of its own: on a chromosome named after the id, at position 1, with the id
# as its A2 and no A1. Such a variant pairs only with the same spelling, and
# matches it exactly, so an id without coordinates is matched by its spelling
# while the rest are matched by their variants.
# @noRd
.unparsedAsSpelled <- function(tab, ids) {
    bad <- which(is.na(tab$chrom) | is.na(tab$pos))
    if (length(bad) == 0L) {
        return(tab)
    }
    spelled <- !is.na(ids[bad])
    bad <- bad[spelled]
    tab$chrom[bad] <- str_c(.UNPARSED_CHROM, ids[bad])
    tab$pos[bad] <- 1L
    tab$A2[bad] <- ids[bad]
    tab$A1[bad] <- ""
    tab
}

# The chromosome-name prefix of the variants .unparsedAsSpelled() makes.
.UNPARSED_CHROM <- "<id>"

# The ids of a variant table's variants: rendered from their fields, or the
# spelling of an id that did not parse (.unparsedAsSpelled()).
# @noRd
.variantTableIds <- function(tab) {
    ids <- formatVariantId(tab$chrom, tab$pos, tab$A2, tab$A1)
    spelled <- which(str_starts(tab$chrom, fixed(.UNPARSED_CHROM)))
    replace(ids, spelled, tab$A2[spelled])
}

# The number of variants in a variant table.
# @noRd
.variantTableSize <- function(tab) {
    length(tab$pos)
}

# A variant table as a data.frame keyed by every field, for hashing variants
# (vctrs) without rendering ids.
# @noRd
.variantKeyFrame <- function(tab) {
    data.frame(chrom = tab$chrom, pos = tab$pos, A2 = tab$A2, A1 = tab$A1)
}

# Variant tables joined, in order.
# @noRd
.bindVariantTables <- function(tables) {
    fields <- c("chrom", "pos", "A2", "A1")
    set_names(
        map(fields, function(f) unlist(map(tables, f), use.names = FALSE)),
        fields
    )
}

# The fields of canonical ids, split on ":".
# @noRd
.canonicalVariantTable <- function(ids) {
    m <- stri_split_fixed(ids, ":", n = 4L, simplify = TRUE)
    list(
        chrom = m[, 1L],
        pos = as.integer(m[, 2L]),
        A2 = m[, 3L],
        A1 = m[, 4L]
    )
}

# A variant table subset to rows `i`.
# @noRd
.variantTableRows <- function(tab, i) {
    map(tab, function(col) col[i])
}

# A panel's variant table, computed on first use: parsing a genome-wide
# panel's ids costs seconds, and when every id matched against it is spelled
# as a panel id, none is needed.
# @noRd
.lazyVariantTable <- function(panel) {
    tab <- NULL
    function() {
        if (is.null(tab)) {
            tab <<- .idsVariantTable(panel)
        }
        tab
    }
}

# Match variant tables `query` to a `panel`: each query variant's index in the
# panel, or NA. A variant with exactly a panel variant's alleles takes that
# one; the rest are matched allele-aware by the matching engine
# (.matchVariantTables(): ref/alt swap, strand), once per distinct variant.
# Strand-ambiguous pairs are matched, since the callers carry orientation-free
# values such as PIPs.
# @noRd
#' @importFrom vctrs vec_match vec_group_id
.alignTablesToPanel <- function(query, panel) {
    qf <- .variantKeyFrame(query)
    idx <- vec_match(qf, .variantKeyFrame(panel))
    unmatched <- which(is.na(idx) & !is.na(query$pos))
    if (length(unmatched) == 0L) {
        return(idx)
    }
    group <- vec_group_id(qf[unmatched, , drop = FALSE])
    first <- unmatched[match(seq_len(attr(group, "n")), group)]
    m <- .matchVariantTables(
        .variantTableRows(query, first),
        panel,
        removeStrandAmbiguous = FALSE
    )
    hit <- rep(NA_integer_, length(first))
    hit[m$idxA] <- m$idxB
    idx[unmatched] <- hit[as.vector(group)]
    idx
}

# Match variant ids to a panel of variant ids: each id's index in `panel`, or
# NA. An id spelled exactly as a panel id takes that one; the rest that parse
# are matched by their variants (.alignTablesToPanel()), and an id that does
# not parse (an rsID) only by its spelling. `panelTable` supplies the panel's
# variant table (see .lazyVariantTable()).
# @noRd
.alignIdsToPanel <- function(
    ids,
    panel,
    panelTable = .lazyVariantTable(panel)
) {
    idx <- match(ids, panel)
    unmatched <- which(is.na(idx))
    if (length(unmatched) == 0L) {
        return(idx)
    }
    rest <- .idsVariantTable(ids[unmatched])
    parsed <- which(!is.na(rest$chrom) & !is.na(rest$pos))
    if (length(parsed) == 0L) {
        return(idx)
    }
    at <- .alignTablesToPanel(.variantTableRows(rest, parsed), panelTable())
    idx[unmatched[parsed]] <- at
    idx
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
# synthetic range would make range() lie, and range() is the block identity the
# (tuple, range) uniqueness key is built on -- a
# fabricated coordinate corrupts both silently.
# @noRd
.variantIdsToGRanges <- function(ids, what = "variant ids") {
    ids <- as.character(ids)
    if (length(ids) == 0L) {
        gr <- S4Vectors::`mcols<-`(
            GenomicRanges::GRanges(),
            value = S4Vectors::DataFrame(
                A1 = character(0),
                A2 = character(0)
            )
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
    gr <- S4Vectors::`mcols<-`(
        GenomicRanges::GRanges(
            parsed$chrom,
            IRanges::IRanges(start = parsed$pos, width = 1L)
        ),
        value = S4Vectors::DataFrame(
            A1 = as.character(parsed$A1),
            A2 = as.character(parsed$A2)
        )
    )
    gr
}
