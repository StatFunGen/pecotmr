#' @importFrom GenomicRanges seqnames
#' @importFrom S4Vectors mcols
NULL

#' Match target data alleles against a reference panel
#'
#' Match by ("chrom", "A1", "A2" and "pos"), accounting for possible
#' strand flips and major/minor allele flips (opposite effects and zscores).
#' Flips specified columns when alleles are swapped relative to the reference.
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
#'   \code{getHarmonizedData()} to recover the post-QC variant
#'   data.frame and \code{getQcSummary()} to inspect the per-variant
#'   merge/flip/strand diagnostics.
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate inner_join filter pull select everything row_number if_else any_of all_of rename
#' @importFrom vctrs vec_duplicate_detect
#' @importFrom tidyr separate
#' @export
matchRefPanel <- function(targetData, refVariants, colToFlip = NULL,
                          matchMinProp = 0.2, removeDups = TRUE,
                          removeIndels = FALSE, removeStrandAmbiguous = TRUE,
                          flipStrand = FALSE, removeUnmatched = TRUE,
                          colToComplement = character(), ...) {
	strandFlip <- function(ref) {
	  chartr("ATCG", "TAGC", ref)
	}

  # helper to sanitize column names to avoid NA/empty names that break dplyr verbs
  sanitizeNames <- function(df) {
	 nm <- colnames(df)
	 if (is.null(nm)) {
	   nm <- rep("unnamed", ncol(df))
	 }
	 emptyIdx <- is.na(nm) | nm == ""
	 if (any(emptyIdx)) {
	   # assign stable placeholder names for empties
	   nm[emptyIdx] <- paste0("unnamed_", seq_len(sum(emptyIdx)))
	 }
	 # ensure names are unique and syntactic
	 nm <- make.unique(nm, sep = "_")
	 colnames(df) <- nm
	 df
  }

  # check if the pattern is ATCG/DI
  checkAtcg <- function(vec) {
	grepl("^[ATCGDI]+$", vec)
  }

  # transform all inputs to dataframe
  if (is.data.frame(targetData)) {
	 if (ncol(targetData) > 4 && all(c("chrom", "pos", "A2", "A1") %in% names(targetData))) {
		# Extract variant columns and standardize
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

  # Remove redundant columns that would conflict with the merge.
  # Keep A2/A1 in targetData - they become A2.target/A1.target after merge.
  columnsToRemove <- c("chromosome", "position", "ref", "alt", "variant_id")
  if (any(columnsToRemove %in% colnames(targetData))) {
	 targetData <- select(targetData, -any_of(columnsToRemove))
  }
  matchResult <- inner_join(targetData, refVariants, by = c("chrom", "pos"), suffix = c(".target", ".ref")) %>% as.data.frame()

  # sanitize names after merge as well (merge can introduce empty names in edge cases)
  matchResult <- sanitizeNames(matchResult)

  if (nrow(matchResult) == 0) {
	warning("No matching variants found between target data and reference variants.")
	return(AlleleQcResult(harmonizedData = matchResult, qcSummary = matchResult))
  }
    # match target & ref by chrom and position
  matchResult = matchResult %>%
	mutate(variants_id_original = formatVariantId(chrom, pos, A2.target, A1.target)) %>%
	mutate(variants_id_qced = formatVariantId(chrom, pos, A2.ref, A1.ref)) %>%
	# upper case target/reference A1 A2
	mutate(across(c(A1.target, A2.target, A1.ref, A2.ref), toupper)) %>%
	mutate(flip1.ref = strandFlip(A1.ref), flip2.ref = strandFlip(A2.ref)) %>%
	# these pairings are ambiguous: because we cannot tell it's an sign flip / strand flip
	mutate(strand_unambiguous = if_else((A1.target == "A" & A2.target == "T") | (A1.target == "T" & A2.target == "A") |
	  (A1.target == "C" & A2.target == "G") | (A1.target == "G" & A2.target == "C"), FALSE, TRUE)) %>%
	# filter out non-ATCG coded alleles
	mutate(non_ATCG = !(checkAtcg(A1.target) & checkAtcg(A2.target))) %>%
	# exact match should be kept all the time
	mutate(exact_match = A1.target == A1.ref & A2.target == A2.ref) %>%
	mutate(sign_flip = ((A1.target == A2.ref & A2.target == A1.ref) | (A1.target == flip2.ref & A2.target == flip1.ref)) & (A1.target != A1.ref & A2.target != A2.ref)) %>%
	mutate(strand_flip = ((A1.target == flip1.ref & A2.target == flip2.ref) | (A1.target == flip2.ref & A2.target == flip1.ref)) & (A1.target != A1.ref & A2.target != A2.ref)) %>%
	mutate(INDEL = (A2.target == "I" | A2.target == "D" | nchar(A2.target) > 1 | nchar(A1.target) > 1)) %>%
	mutate(ID_match = ((A2.target == "D" | A2.target == "I") & (nchar(A1.ref) > 1 | nchar(A2.ref) > 1)))

  # if not remove, then this should'nt be a condition to filter out any variants
  if (!removeStrandAmbiguous) {
	matchResult <- matchResult %>% mutate(strand_unambiguous = TRUE)
  }
  # if all strand flip is un-ambigous, then we know ambigous cases are indeed a strand flip
  # not a sign flip, then we infer there is no ambigous in the whole dataset, and keep those ambigous ones
  if (!any(matchResult$strand_flip & matchResult$strand_unambiguous)) {
	matchResult <- matchResult %>% mutate(strand_unambiguous = TRUE)
  }

  # To keep variants: if it's a strand flip, we will keep those unambiguous (because if ambiguous, cannot know if it's strand / sign flip, so discard all)
  # or exact match or indel match (ID_match)
  # If not a strand flip, then we will keep those that are exact match / those are sign flip / INDEL matched
  matchResult <- matchResult %>% mutate(keep = if_else(strand_flip, true = (strand_unambiguous | exact_match | ID_match), false =
	(exact_match | sign_flip | ID_match)
  ))

  if (removeIndels) {
	matchResult <- matchResult %>% mutate(keep = if_else(INDEL, FALSE, keep))
  }

  # Flip the signs of colToFlip for sign-flipped variants
  if (!is.null(colToFlip)) {
	missing <- setdiff(colToFlip, colnames(matchResult))
	if (length(missing) > 0) {
	  stop("Column(s) '", paste(missing, collapse = "', '"), "' not found in targetData.")
	}
	matchResult[matchResult$sign_flip, colToFlip] <- -1 * matchResult[matchResult$sign_flip, colToFlip]
  }
  # Complement (1 - x) colToComplement for the same swapped variants. A
  # frequency tracks the effect allele, so an allele swap takes af -> 1 - af
  # (not a sign flip). Kept independent of colToFlip so signed columns are
  # untouched here.
  if (length(colToComplement) > 0) {
	missing <- setdiff(colToComplement, colnames(matchResult))
	if (length(missing) > 0) {
	  stop("Column(s) '", paste(missing, collapse = "', '"), "' not found in targetData.")
	}
	matchResult[matchResult$sign_flip, colToComplement] <-
	  1 - matchResult[matchResult$sign_flip, colToComplement]
  }
  # flip the strands if there is a strand flip
  if (flipStrand) {
	strandFlippedIndices <- which(matchResult$strand_flip)
	matchResult[strandFlippedIndices, "A1.target"] <- strandFlip(matchResult[strandFlippedIndices, "A1.target"])
	matchResult[strandFlippedIndices, "A2.target"] <- strandFlip(matchResult[strandFlippedIndices, "A2.target"])
  }

  # Remove all unnecessary columns used to determine qc status
  # Finally keep those variants with FLAG keep = TRUE
  result <- matchResult[matchResult$keep, , drop = FALSE]

  if (removeDups) {
	dups <- duplicated(result[, c("chrom", "pos", "variants_id_qced")])
	if (any(dups)) {
	  nRemoved <- sum(dups)
	  warning(sprintf("Removed %d duplicate variant(s), keeping first occurrence.", nRemoved))
	  result <- result[!dups, , drop = FALSE]
	}
  }

  qcCols <- c("flip1.ref", "flip2.ref", "strand_unambiguous", "non_ATCG",
              "exact_match", "sign_flip", "strand_flip", "INDEL", "ID_match", "keep")
  result <- result %>%
	select(-any_of(qcCols), -A1.target, -A2.target) %>%
	rename(A1 = A1.ref, A2 = A2.ref, variant_id = variants_id_qced)

  if (!removeUnmatched) {
	matchVariant <- result %>% pull(variants_id_original)
	matchResult <- matchResult %>%
	  select(-any_of(qcCols), -variants_id_original, -A1.target, -A2.target) %>%
	  rename(A1 = A1.ref, A2 = A2.ref, variant_id = variants_id_qced)
	targetData <- targetData %>% mutate(variant_id = formatVariantId(chrom, pos, A2, A1))
	if (length(setdiff(targetData %>% pull(variant_id), matchVariant)) > 0) {
	  unmatchData <- targetData %>% filter(!variant_id %in% matchVariant)
	  result <- rbind(result, unmatchData %>% mutate(variants_id_original = variant_id))
	  result <- result[match(targetData$variant_id, result$variants_id_original), ] %>% select(-variants_id_original)
	}
  }

  minMatch <- matchMinProp * nrow(refVariants)
  if (nrow(result) < minMatch) {
	stop("Not enough variants have been matched.")
  }

  # Error if duplicated variant IDs remain (same variant with different values)
  if (any(duplicated(result$variant_id))) {
	stop("Duplicated variants with different values found. Please check the input data and determine which to keep.")
  }

  return(AlleleQcResult(harmonizedData = result, qcSummary = matchResult))
}

#' @rdname matchRefPanel
#' @export
alleleQc <- matchRefPanel

#' Align Variant Names
#'
#' This function aligns variant names from two strings containing variant names in the format of
#' "chr:pos:A1:A2" or "chr:pos_A1_A2". The first string should be the "source" and the second
#' should be the "reference".
#'
#' @param source A character vector of variant names in the format "chr:pos:A2:A1" or "chr:pos_A2_A1".
#' @param reference A character vector of variant names in the format "chr:pos:A2:A1" or "chr:pos_A2_A1".
#' @param removeBuildSuffix Whether to strip trailing genome build suffixes like ":b38" or "_b38" before alignment. Default TRUE.
#'
#' @return A list with two elements:
#' - aligned_variants: A character vector of aligned variant names.
#' - unmatched_indices: A vector of indices for the variants in the source that could not be matched.
#'
#' @examples
#' source <- c("1:123:A:C", "2:456:G:T", "3:789:C:A")
#' reference <- c("1:123:A:C", "2:456:T:G", "4:101:G:C")
#' alignVariantNames(source, reference)
#'
#' @export
alignVariantNames <- function(source, reference, removeIndels = FALSE, removeBuildSuffix = TRUE) {
  # Optionally strip build suffix like :b38 or _b38 from both sides for robust alignment
  if (removeBuildSuffix) {
    source <- gsub("(:|_)b[0-9]+$", "", source)
    reference <- gsub("(:|_)b[0-9]+$", "", reference)
  }
  # Check if source and reference follow the expected pattern
  sourcePattern <- grepl("^(chr)?[0-9]+[_:][0-9]+[_:][ATCG*]+[_:][ATCG*]+$", source)
  referencePattern <- grepl("^(chr)?[0-9]+[_:][0-9]+[_:][ATCG*]+[_:][ATCG*]+$", reference)

  if (!all(sourcePattern) && !all(referencePattern)) {
    warning("Cannot unify variant names because they do not follow the expected variant naming convention chr:pos:A2:A1 or chr:pos_A2_A1.")
    return(list(aligned_variants = source, unmatched_indices = integer(0)))
  }

  if ((!all(sourcePattern) && all(referencePattern)) || (all(sourcePattern) && !all(referencePattern))) {
    stop("Source and reference have different variant naming conventions. They cannot be aligned.")
  }

  # Detect reference convention to preserve in output
  refConvention <- detectVariantConvention(reference)

  sourceDf <- parseVariantId(source)
  referenceDf <- parseVariantId(reference)

  qcResult <- matchRefPanel(
    targetData = sourceDf,
    refVariants = referenceDf,
    colToFlip = NULL,
    matchMinProp = 0,
    removeDups = FALSE,
    flipStrand = TRUE,
    removeIndels = removeIndels,
    removeStrandAmbiguous = FALSE,
    removeUnmatched = FALSE
  )

  alignedDf <- getHarmonizedData(qcResult)

  # Format output using reference convention (preserving user's format automatically)
  alignedVariants <- formatVariantId(
    alignedDf$chrom, alignedDf$pos, alignedDf$A2, alignedDf$A1,
    convention = refConvention
  )
  names(alignedVariants) <- NULL

  # Normalize reference to the same output format for accurate matching
  refNormalized <- normalizeVariantId(reference, convention = refConvention)
  unmatchedIndices <- which(match(alignedVariants, refNormalized, nomatch = 0) == 0)

  list(
    aligned_variants = alignedVariants,
    unmatched_indices = unmatchedIndices
  )
}

#' Merge variant info from two sources with allele-flip-aware matching
#'
#' Merges variant metadata (chromosome, position, ref, alt) from two sources,
#' detecting and correcting allele flips (where alt/ref are swapped). Creates
#' a canonical key from sorted alleles to match across datasets.
#'
#' @param variants1 A data.frame with columns \code{chrom}, \code{pos},
#'   \code{alt}, \code{ref}, or a \code{GRanges} with corresponding metadata
#'   columns.
#' @param variants2 A data.frame or \code{GRanges} with the same columns.
#' @param all Logical. If TRUE (default), returns the union of both sets.
#'   If FALSE, returns only variants from \code{variants2} (flipped to match
#'   \code{variants1}'s allele orientation).
#' @return A data.frame with columns \code{chrom}, \code{pos}, \code{alt},
#'   \code{ref}, deduplicated by position and alleles.
#' @export
mergeVariantInfo <- function(variants1, variants2, all = TRUE) {
  # Convert GRanges to data.frame if needed
  toDf <- function(x) {
    if (is(x, "GRanges")) {
      mc <- as.data.frame(mcols(x))
      mc$chrom <- as.character(seqnames(x))
      mc$pos <- start(x)
      mc[, c("chrom", "pos", "alt", "ref")]
    } else {
      as.data.frame(x)[, c("chrom", "pos", "alt", "ref")]
    }
  }

  df1 <- toDf(variants1)
  df2 <- toDf(variants2)

  # Create a canonical key from sorted alleles so flipped pairs match
  makeKey <- function(df) {
    aMin <- pmin(df$alt, df$ref)
    aMax <- pmax(df$alt, df$ref)
    paste(df$chrom, df$pos, aMin, aMax)
  }

  key1 <- makeKey(df1)
  key2 <- makeKey(df2)

  # Detect flips: where df2's alt matches df1's ref at the same key
  matchIdx <- match(key2, key1)
  hasMatch <- !is.na(matchIdx)

  flip <- hasMatch &
    df2$alt[hasMatch] == df1$ref[matchIdx[hasMatch]] &
    df2$ref[hasMatch] == df1$alt[matchIdx[hasMatch]]

  # Apply flips to df2
  flipRows <- which(hasMatch)[flip[hasMatch]]
  if (length(flipRows) > 0) {
    tmp <- df2$alt[flipRows]
    df2$alt[flipRows] <- df2$ref[flipRows]
    df2$ref[flipRows] <- tmp
  }

  if (all) {
    combined <- rbind(
      df1[, c("chrom", "pos", "alt", "ref")],
      df2[, c("chrom", "pos", "alt", "ref")])
    combined[!duplicated(paste(combined$chrom, combined$pos,
                               combined$alt, combined$ref)), ]
  } else {
    df2[, c("chrom", "pos", "alt", "ref")]
  }
}

