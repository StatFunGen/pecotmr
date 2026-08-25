# =============================================================================
# LdStatistic S4 virtual class
# -----------------------------------------------------------------------------
# Abstract container for pre-computed LD statistics. Subclasses
# (LdEigen, LdScore) provide method-specific representations: LdEigen
# for eigendecomposition-based methods (LDER/HDL/sHDL), LdScore for
# LD-score-based methods (S-LDSC/g-LDSC).
# =============================================================================

#' @include AllGenerics.R LdBlocks.R
NULL

#' @title LD Statistic (Virtual Base Class)
#' @description Abstract container for pre-computed LD statistics. Subclasses
#'   provide method-specific representations: eigendecompositions (for
#'   LDER/HDL/sHDL) and LD score matrices (for S-LDSC/g-LDSC).
#'
#'   An \code{LdStatistic} \strong{is} a \code{GRanges} of the reference
#'   variants: one range per SNP, with \code{SNP} / \code{A1} / \code{A2}
#'   (and optionally \code{MAF}) in \code{mcols}. Per-variant statistics
#'   live in \code{mcols} too, so subsetting narrows the ranges and the
#'   statistics together. The genome build comes from \code{seqinfo()} rather
#'   than a slot of its own.
#' @slot ldBlocks A \code{GRanges} of LD block intervals.
#' @slot nRef Integer, sample size of the LD reference panel.
#' @slot inSample Logical, whether the LD reference is from the same cohort as
#'   the GWAS (affects bias correction).
#' @export
setClass(
    "LdStatistic",
    contains = c("VIRTUAL", "GRanges"),
    representation(
        ldBlocks = "GRanges",
        nRef = "integer",
        inSample = "logical"
    ),
    validity = function(object) .validateLdStatistic(object)
)

# @noRd
.validateLdStatistic <- function(object) {
    errors <- character()
    if (length(object@nRef) != 1L || object@nRef <= 0L) {
        errors <- c(errors, "'nRef' must be a single positive integer")
    }
    if (length(object@inSample) != 1L) {
        errors <- c(errors, "'inSample' must be a single logical value")
    }
    if (length(object) == 0L) {
        errors <- c(errors, "an LdStatistic must carry at least one variant")
    }
    if (length(errors) == 0) TRUE else errors
}

# The variants of an LD reference, as the GRanges every subclass is built on.
# Shared by the LdScore() and LdEigen() constructors so the two agree on what
# a reference panel's variant table looks like.
# @noRd
.ldStatRanges <- function(snpInfo, genome) {
    required <- c("SNP", "CHR", "BP", "A1", "A2")
    missingCols <- setdiff(required, colnames(snpInfo))
    if (length(missingCols) > 0L) {
        abort(glue(
            "`snpInfo` is missing column(s): ",
            "{str_flatten(missingCols, ', ')}."
        ))
    }
    gr <- GenomicRanges::GRanges(
        seqnames = withChrPrefix(as.character(snpInfo$CHR)),
        ranges = IRanges::IRanges(as.integer(snpInfo$BP), width = 1L)
    )
    names(gr) <- as.character(snpInfo$SNP)
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        select(snpInfo, -any_of(c("CHR", "BP"))),
        row.names = NULL
    )
    if (!is.null(genome) && length(genome) == 1L && nzchar(genome)) {
        GenomeInfoDb::genome(gr) <- genome
    }
    gr
}

#' @rdname getNRef
#' @export
setMethod("getNRef", "LdStatistic", function(x) x@nRef)

#' @rdname getInSample
#' @export
setMethod("getInSample", "LdStatistic", function(x) x@inSample)

#' @rdname getLdBlocks
#' @export
setMethod("getLdBlocks", "LdStatistic", function(x) x@ldBlocks)

#' @rdname getGenome
#' @export
setMethod("getGenome", "LdStatistic", function(x, ...) {
    # The build lives in seqinfo, not a slot: a GRanges already has somewhere
    # to keep it, and storing it twice is what the retired LdBlocks class did.
    build <- unique(GenomeInfoDb::genome(x))
    build <- build[!is.na(build)]
    if (length(build) == 0L) NA_character_ else build[[1L]]
})
