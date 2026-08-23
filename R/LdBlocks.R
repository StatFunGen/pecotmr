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

setClass(
    "LdBlocks",
    representation(
        blocks = "GRanges",
        genome = "character"
    ),
    validity = function(object) {
        errors <- character()
        if (length(object@genome) != 1L) {
            errors <- c(errors, "'genome' must be a single character string")
        }
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


#' @rdname show-methods
#' @export
setMethod("show", "LdBlocks", function(object) {
    cat(glue(
        "LdBlocks: {length(object@blocks)} blocks, ",
        "genome build: {object@genome}\n",
        .trim = FALSE
    ))
})


# Coerce an LD-block specification to the GRanges the splitters consume.
#
# Accepts the canonical `LdBlocks`, a bare `GRanges`, a data.frame with
# chrom/start/end (plus an optional `blockId`), or a path to a delimited file
# holding those columns. `LdBlocks` is listed first because it is the type the
# rest of the package passes around; the looser forms exist so a caller with a
# plain block table does not have to construct one.
# @noRd
.asLdBlockRanges <- function(x) {
    if (methods::is(x, "LdBlocks")) {
        return(getBlocks(x))
    }
    if (methods::is(x, "GRanges")) {
        return(x)
    }
    if (is.character(x) && length(x) == 1L) {
        x <- .readManifest(x)
    }
    if (!is.data.frame(x)) {
        msg <- glue(
            "`ldBlocks` must be an LdBlocks, a GRanges, a data.frame with ",
            "chrom/start/end, or a path to one (got {class(x)[[1L]]})."
        )
        abort(msg)
    }
    .ldBlockTableToRanges(x)
}

# One block per row. `blockId` is carried into mcols when present so the key
# survives; without it the splitter falls back to the block's coordinates.
# @noRd
.ldBlockTableToRanges <- function(df) {
    missingCols <- setdiff(c("chrom", "start", "end"), colnames(df))
    if (length(missingCols) > 0L) {
        msg <- glue(
            "`ldBlocks` table is missing required column(s): ",
            "{str_flatten(missingCols, ', ')}."
        )
        abort(msg)
    }
    gr <- GenomicRanges::GRanges(
        seqnames = as.character(df$chrom),
        ranges = IRanges::IRanges(
            start = as.integer(df$start),
            end = as.integer(df$end)
        )
    )
    if (is_in("blockId", colnames(df))) {
        mcols(gr)$blockId <- as.character(df$blockId)
    }
    gr
}
