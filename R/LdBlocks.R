# =============================================================================
# LdBlocks S4 class
# -----------------------------------------------------------------------------
# Container for genome-partitioned LD block boundaries. Holds a GRanges
# of block intervals plus a `genome` build label. Consumed by h2
# estimation (per-block jackknife), LD score computation, and LD-block-
# indexed GWAS fine-mapping.
# =============================================================================

#' @include AllGenerics.R
NULL

setClass("LdBlocks",
  representation(
    blocks = "GRanges",
    genome = "character"
  ),
  validity = function(object) {
    errors <- character()
    if (length(object@genome) != 1L)
      errors <- c(errors, "'genome' must be a single character string")
    if (length(errors) == 0) TRUE else errors
  }
)

# =============================================================================
# Genotype Handle
# =============================================================================

#' @title Genotype Handle
#' @description Lightweight handle to genotype data in any supported format.
#'   Stores the file path, detected format, and cached SNP metadata. Used to
#'   defer reading genotypes until block-level extraction is needed.
#' @slot path Character, path to the genotype file (or stem for plink).
#' @slot format Character, one of "gds", "vcf", "plink1", "plink2".
#' @slot snpInfo A \code{data.frame} with columns \code{SNP}, \code{CHR},
#'   \code{BP}, \code{A1}, \code{A2}. Cached on first access.
#' @slot nSamples Integer, number of samples.
#' @slot sampleIds Character vector of sample identifiers.
#' @slot pgenPtr An external pointer for plink2 pgen handle, or NULL.
#' @export


#' @rdname getBlocks
#' @export
setMethod("getBlocks", "LdBlocks", function(x) x@blocks)

#' @rdname getGenome
#' @export
setMethod("getGenome", "LdBlocks", function(x, ...) x@genome)


#' @export
setMethod("show", "LdBlocks", function(object) {
  cat(sprintf("LdBlocks: %d blocks, genome build: %s\n",
              length(object@blocks), object@genome))
})
