#' Match alleles between target data and reference variants
#'
#' Match by ("chrom", "A1", "A2" and "pos"), accounting for possible
#' strand flips and major/minor allele flips (opposite effects and zscores).
#'
#' @param target_data A data frame with columns "chrom", "pos", "A2", "A1" (and optionally other columns like "beta" or "z"),
#'   or a vector of strings in the format of "chr:pos:A2:A1"/"chr:pos_A2_A1". Can be automatically converted to a data frame if a vector.
#' @param ref_variants A data frame with columns "chrom", "pos", "A2", "A1" or strings in the format of "chr:pos:A2:A1"/"chr:pos_A2_A1".
#' @param col_to_flip The name of the column in target_data where flips are to be applied.
#' @param match_min_prop Minimum proportion of variants in the smallest data
#'   to be matched, otherwise stops with an error. Default is 20%.
#' @param remove_dups Whether to remove duplicates, default is TRUE.
#' @param remove_indels Whether to remove INDELs, default is FALSE.
#' @param flip Whether the alleles must be flipped: A <--> T & C <--> G, in which case
#'   corresponding `col_to_flip` are multiplied by -1. Default is `TRUE`.
#' @param remove_strand_ambiguous Whether to remove strand SNPs (if any). Default is `TRUE`.
#' @param flip_strand Whether to output the variants after strand flip. Default is `FALSE`.
#' @param remove_unmatched Whether to remove unmatched variants. Default is `TRUE`.
#' @return A single data frame with matched variants.
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate inner_join filter pull select everything row_number if_else any_of all_of rename
#' @importFrom vctrs vec_duplicate_detect
#' @importFrom tidyr separate
#' @export
allele_qc <- function(target_data, ref_variants, col_to_flip = NULL,
                     match_min_prop = 0.2, remove_dups = TRUE,
                     remove_indels = FALSE, remove_strand_ambiguous = TRUE,
                     flip_strand = FALSE, remove_unmatched = TRUE, ...) {
	strand_flip <- function(ref) {
	# Define a mapping for complementary bases
	base_mapping <- c("A" = "T", "T" = "A", "G" = "C", "C" = "G")

	# Function to complement a single base
	complement_base <- function(base) {
	  complement <- base_mapping[base]
	  return(complement)
	}

	# Function to reverse and complement a DNA sequence
	reverse_complement <- function(sequence) {
	  reversed <- rev(strsplit(sequence, NULL)[[1]])
	  complemented <- sapply(reversed, complement_base, USE.NAMES = FALSE)
	  return(paste(complemented, collapse = ""))
	}

	complemented_sequence <- sapply(ref, reverse_complement, USE.NAMES = FALSE)

	return(complemented_sequence)
  }

  # helper to sanitize column names to avoid NA/empty names that break dplyr verbs
  sanitize_names <- function(df) {
	 nm <- colnames(df)
	 if (is.null(nm)) {
	   nm <- rep("unnamed", ncol(df))
	 }
	 empty_idx <- is.na(nm) | nm == ""
	 if (any(empty_idx)) {
	   # assign stable placeholder names for empties
	   nm[empty_idx] <- paste0("unnamed_", seq_len(sum(empty_idx)))
	 }
	 # ensure names are unique and syntactic
	 nm <- base::make.unique(nm, sep = "_")
	 colnames(df) <- nm
	 df
  }

  # check if the pattern is ATCG/DI
  check_ATCG <- function(vec) {
	pattern <- "^[ATCGDI]+$"

	# Function to check if a single element matches the pattern
	check_element <- function(element) {
	  grepl(pattern, element)
	}

	result <- sapply(vec, check_element, USE.NAMES = FALSE)
	return(result)
  }

  # transform all inputs to dataframe
  if (is.data.frame(target_data)) {
	 if (ncol(target_data) > 4 && all(c("chrom", "pos", "A2", "A1") %in% names(target_data))) {
		# Extract variant columns and standardize
		variant_cols <- c("chrom", "pos", "A2", "A1")
		variant_df <- target_data %>% select(all_of(variant_cols))
		other_cols <- target_data %>% select(-all_of(variant_cols))
		target_data <- cbind(variant_id_to_df(variant_df), other_cols)
	 } else {
		target_data <- variant_id_to_df(target_data)
	 }
  } else {
		target_data <- variant_id_to_df(target_data)
  }
  ref_variants <- variant_id_to_df(ref_variants)

  columns_to_remove <- c("chromosome", "position", "A2", "A1", "ref", "alt", "variant_id")

  # Check if any of the specified columns are present
  if (any(columns_to_remove %in% colnames(target_data))) {
	 target_data <- select(target_data, -any_of(columns_to_remove))
  }
  match_result <- merge(target_data, ref_variants, by = c("chrom", "pos"), all = FALSE, suffixes = c(".target", ".ref")) %>% as.data.frame()

  # sanitize names after merge as well (merge can introduce empty names in edge cases)
  match_result <- sanitize_names(match_result)

  if (nrow(match_result) == 0) {
	warning("No matching variants found between target data and reference variants.") 
	return(list(target_data_qced = match_result, qc_summary = match_result))
  }
    # match target & ref by chrom and position
  match_result = match_result %>%
	mutate(variants_id_original = format_variant_id(chrom, pos, A2.target, A1.target)) %>%
	mutate(variants_id_qced = format_variant_id(chrom, pos, A2.ref, A1.ref)) %>%
	# upper case target/reference A1 A2
	mutate(across(c(A1.target, A2.target, A1.ref, A2.ref), toupper)) %>%
	mutate(flip1.ref = strand_flip(A1.ref), flip2.ref = strand_flip(A2.ref)) %>%
	# these pairings are ambiguous: because we cannot tell it's an sign flip / strand flip
	mutate(strand_unambiguous = if_else((A1.target == "A" & A2.target == "T") | (A1.target == "T" & A2.target == "A") |
	  (A1.target == "C" & A2.target == "G") | (A1.target == "G" & A2.target == "C"), FALSE, TRUE)) %>%
	# filter out non-ATCG coded alleles
	mutate(non_ATCG = !(check_ATCG(A1.target) & check_ATCG(A2.target))) %>%
	# exact match should be kept all the time
	mutate(exact_match = A1.target == A1.ref & A2.target == A2.ref) %>%
	mutate(sign_flip = ((A1.target == A2.ref & A2.target == A1.ref) | (A1.target == flip2.ref & A2.target == flip1.ref)) & (A1.target != A1.ref & A2.target != A2.ref)) %>%
	mutate(strand_flip = ((A1.target == flip1.ref & A2.target == flip2.ref) | (A1.target == flip2.ref & A2.target == flip1.ref)) & (A1.target != A1.ref & A2.target != A2.ref)) %>%
	mutate(INDEL = (A2.target == "I" | A2.target == "D" | nchar(A2.target) > 1 | nchar(A1.target) > 1)) %>%
	mutate(ID_match = ((A2.target == "D" | A2.target == "I") & (nchar(A1.ref) > 1 | nchar(A2.ref) > 1)))

  # if not remove, then this should'nt be a condition to filter out any variants
  if (!remove_strand_ambiguous) {
	match_result <- match_result %>% mutate(strand_unambiguous = TRUE)
  }
  # if all strand flip is un-ambigous, then we know ambigous cases are indeed a strand flip
  # not a sign flip, then we infer there is no ambigous in the whole dataset, and keep those ambigous ones
  if (nrow(match_result %>% filter(strand_flip == TRUE) %>% filter(strand_unambiguous == TRUE)) == 0) {
	match_result <- match_result %>% mutate(strand_unambiguous = TRUE)
  }

  # To keep variants: if it's a strand flip, we will keep those unambigous (because if ambigous, cannot know it's trand / sign flip, so discard all)
  # or exact match or indel match (ID_match)
  # If not a strand flip, then we will keep those that are exact match / those are sign flip / INDEL matched
  match_result <- match_result %>% mutate(keep = if_else(strand_flip, true = (strand_unambiguous | exact_match | ID_match), false =
	(exact_match | sign_flip | ID_match)
  ))

  if (remove_indels) {
	match_result <- match_result %>% mutate(keep = if_else(INDEL, FALSE, keep))
  }

  # Flip the signs of col_to_flip for sign-flipped variants
  if (!is.null(col_to_flip)) {
	missing <- setdiff(col_to_flip, colnames(match_result))
	if (length(missing) > 0) {
	  stop("Column(s) '", paste(missing, collapse = "', '"), "' not found in target_data.")
	}
	match_result[match_result$sign_flip, col_to_flip] <- -1 * match_result[match_result$sign_flip, col_to_flip]
  }
  # flip the strands if there is a strand flip
  if (flip_strand) {
	strand_flipped_indices <- which(match_result$strand_flip)
	match_result[strand_flipped_indices, "A1.target"] <- strand_flip(match_result[strand_flipped_indices, "A1.target"])
	match_result[strand_flipped_indices, "A2.target"] <- strand_flip(match_result[strand_flipped_indices, "A2.target"])
  }

  # Remove all unnecessary columns used to determine qc status
  # Finally keep those variants with FLAG keep = TRUE
  result <- match_result[match_result$keep, , drop = FALSE]
	
  if (remove_dups) {
	dups <- vec_duplicate_detect(result[, c("chrom", "pos", "variants_id_qced")])
	if (any(dups)) {
	  result <- result[!dups, , drop = FALSE]
	  warning("Duplicate variants were removed.")
	}
  }

  qc_cols <- c("flip1.ref", "flip2.ref", "strand_unambiguous", "non_ATCG",
               "exact_match", "sign_flip", "strand_flip", "INDEL", "ID_match", "keep")
  result <- result %>%
	select(-any_of(qc_cols), -A1.target, -A2.target) %>%
	rename(A1 = A1.ref, A2 = A2.ref, variant_id = variants_id_qced)

  if (!remove_unmatched) {
	match_variant <- result %>% pull(variants_id_original)
	match_result <- match_result %>%
	  select(-any_of(qc_cols), -variants_id_original, -A1.target, -A2.target) %>%
	  rename(A1 = A1.ref, A2 = A2.ref, variant_id = variants_id_qced)
	target_data <- target_data %>% mutate(variant_id = format_variant_id(chrom, pos, A2, A1))
	if (length(setdiff(target_data %>% pull(variant_id), match_variant)) > 0) {
	  unmatch_data <- target_data %>% filter(!variant_id %in% match_variant)
	  result <- rbind(result, unmatch_data %>% mutate(variants_id_original = variant_id))
	  result <- result[match(target_data$variant_id, result$variants_id_original), ] %>% select(-variants_id_original)
	}
  }

  min_match <- match_min_prop * nrow(ref_variants)
  if (nrow(result) < min_match) {
	stop("Not enough variants have been matched.")
  }

  # Error if duplicated variant IDs remain (same variant with different values)
  if (any(duplicated(result$variant_id))) {
	stop("Duplicated variants with different values found. Please check the input data and determine which to keep.")
  }

  return(list(target_data_qced = result, qc_summary = match_result))
}

#' Align Variant Names
#'
#' This function aligns variant names from two strings containing variant names in the format of
#' "chr:pos:A1:A2" or "chr:pos_A1_A2". The first string should be the "source" and the second
#' should be the "reference".
#'
#' @param source A character vector of variant names in the format "chr:pos:A2:A1" or "chr:pos_A2_A1".
#' @param reference A character vector of variant names in the format "chr:pos:A2:A1" or "chr:pos_A2_A1".
#' @param remove_build_suffix Whether to strip trailing genome build suffixes like ":b38" or "_b38" before alignment. Default TRUE.
#'
#' @return A list with two elements:
#' - aligned_variants: A character vector of aligned variant names.
#' - unmatched_indices: A vector of indices for the variants in the source that could not be matched.
#'
#' @examples
#' source <- c("1:123:A:C", "2:456:G:T", "3:789:C:A")
#' reference <- c("1:123:A:C", "2:456:T:G", "4:101:G:C")
#' align_variant_names(source, reference)
#'
#' @export
align_variant_names <- function(source, reference, remove_indels = FALSE, remove_build_suffix = TRUE) {
  # Optionally strip build suffix like :b38 or _b38 from both sides for robust alignment
  if (remove_build_suffix) {
    source <- gsub("(:|_)b[0-9]+$", "", source)
    reference <- gsub("(:|_)b[0-9]+$", "", reference)
  }
  # Check if source and reference follow the expected pattern
  source_pattern <- grepl("^(chr)?[0-9]+[_:][0-9]+[_:][ATCG*]+[_:][ATCG*]+$", source)
  reference_pattern <- grepl("^(chr)?[0-9]+[_:][0-9]+[_:][ATCG*]+[_:][ATCG*]+$", reference)

  if (!all(source_pattern) && !all(reference_pattern)) {
    warning("Cannot unify variant names because they do not follow the expected variant naming convention chr:pos:A2:A1 or chr:pos_A2_A1.")
    return(list(aligned_variants = source, unmatched_indices = integer(0)))
  }

  if ((!all(source_pattern) && all(reference_pattern)) || (all(source_pattern) && !all(reference_pattern))) {
    stop("Source and reference have different variant naming conventions. They cannot be aligned.")
  }

  # Detect reference convention to preserve in output
  ref_convention <- detect_variant_convention(reference)

  source_df <- parse_variant_id(source)
  reference_df <- parse_variant_id(reference)

  qc_result <- allele_qc(
    target_data = source_df,
    ref_variants = reference_df,
    col_to_flip = NULL,
    match_min_prop = 0,
    remove_dups = FALSE,
    flip_strand = TRUE,
    remove_indels = remove_indels,
    remove_strand_ambiguous = FALSE,
    remove_unmatched = FALSE
  )

  aligned_df <- qc_result$target_data_qced

  # Format output using reference convention (preserving user's format automatically)
  aligned_variants <- format_variant_id(
    aligned_df$chrom, aligned_df$pos, aligned_df$A2, aligned_df$A1,
    convention = ref_convention
  )
  names(aligned_variants) <- NULL

  # Normalize reference to the same output format for accurate matching
  ref_normalized <- normalize_variant_id(reference, convention = ref_convention)
  unmatched_indices <- which(match(aligned_variants, ref_normalized, nomatch = 0) == 0)

  list(
    aligned_variants = aligned_variants,
    unmatched_indices = unmatched_indices
  )
}