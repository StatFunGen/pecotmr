# =============================================================================
# LdData S4 class
# -----------------------------------------------------------------------------
# Container for LD information used by fine-mapping / colocalization. Holds
# a pre-computed correlation matrix (or a list of per-block matrices) and/or
# a GenotypeHandle (or list of handles for mixture panels) for on-demand
# correlation computation.
# =============================================================================

#' @include GenotypeHandle.R
NULL

#' @title LD Data Container
#' @description S4 container for LD information. Stores either a pre-computed
#'   correlation matrix or a \code{GenotypeHandle} (or list of handles for
#'   mixture panels) for lazy genotype/correlation access.
#'
#' @slot correlation A correlation matrix, a list of per-block matrices
#'   (block-diagonal LD), or NULL if genotypes are available and R should be
#'   computed on demand.
#' @slot genotypeHandle A \code{GenotypeHandle}, a list of
#'   \code{GenotypeHandle}s (for mixture panels), or NULL when only pre-computed
#'   R is available.
#' @slot snpIdx Integer vector of 1-based SNP indices into the handle's
#'   \code{snpInfo}. NULL when correlation is pre-computed.
#' @slot variants A \code{GRanges} object with variant metadata (A1, A2,
#'   variant_id, and optionally allele_freq, variance, n_nomiss).
#' @slot blockMetadata An \code{LdBlocks} object or a \code{data.frame} with
#'   block boundary information.
#' @slot nRef Integer, reference panel sample size.
#' @slot mixtureWeights NULL when \code{genotypeHandle} is a single
#'   \code{GenotypeHandle}; a numeric vector of mixing proportions (one per
#'   panel, summing to 1) when \code{genotypeHandle} is a list of
#'   \code{GenotypeHandle}s. Used by \code{getCorrelation()} to compute a
#'   weighted-average mixture LD matrix; required whenever \code{genotypeHandle}
#'   is a list and \code{getCorrelation()} will be called.
#' @export
setClass(
    "LdData",
    representation(
        correlation = "ANY",
        genotypeHandle = "ANY",
        snpIdx = "ANY",
        variants = "GRanges",
        blockMetadata = "ANY",
        nRef = "integer",
        mixtureWeights = "ANY"
    ),
    validity = function(object) {
        errors <- character()
        if (is.null(object@correlation) && is.null(object@genotypeHandle)) {
            errors <- c(
                errors,
                str_c(
                    "At least one of 'correlation' or ",
                    "'genotypeHandle' must be non-NULL"
                )
            )
        }
        if (length(object@variants) == 0) {
            errors <- c(errors, "'variants' must not be empty")
        }
        if (!is.null(object@mixtureWeights)) {
            if (!is.list(object@genotypeHandle)) {
                errors <- c(
                    errors,
                    str_c(
                        "'mixtureWeights' may only be set when ",
                        "'genotypeHandle' is a list of GenotypeHandles"
                    )
                )
            } else {
                w <- object@mixtureWeights
                if (
                    !is.numeric(w) || length(w) != length(object@genotypeHandle)
                ) {
                    errors <- c(
                        errors,
                        str_c(
                            "'mixtureWeights' must be numeric of length ",
                            "equal to the genotypeHandle list"
                        )
                    )
                } else if (any(w < 0) || abs(sum(w) - 1) > 1e-6) {
                    errors <- c(
                        errors,
                        "'mixtureWeights' must be non-negative and sum to 1"
                    )
                }
            }
        }
        if (length(errors) == 0) TRUE else errors
    }
)

#' @rdname show-methods
#' @export
setMethod("show", "LdData", function(object) {
    n_var <- length(object@variants)
    has_R <- !is.null(object@correlation)
    has_geno <- !is.null(object@genotypeHandle)
    r_type <- if (has_R && is.list(object@correlation)) {
        "block-diagonal"
    } else {
        "single"
    }
    cat(glue("LdData: {n_var} variants\n", .trim = FALSE))
    cat(glue(
        "  Correlation: {if (has_R) r_type else 'NULL'}, ",
        "Genotype handle: {if (has_geno) 'available' else 'NULL'}\n",
        .trim = FALSE
    ))
    cat(glue("  Reference N: {object@nRef}\n", .trim = FALSE))
})

#' @title Create an LdData Object
#' @description Construct an \code{LdData} from a correlation matrix and/or
#'   genotype handle, plus variant metadata as a GRanges.
#' @param correlation A correlation matrix, list of matrices, or NULL.
#' @param genotypeHandle A GenotypeHandle, list of GenotypeHandles, or NULL.
#' @param snpIdx Integer vector of SNP indices, or NULL.
#' @param variants A GRanges with variant metadata (must have variant_id in
#'   mcols, plus A1, A2).
#' @param blockMetadata LdBlocks or data.frame with block info.
#' @param nRef Integer, reference panel sample size.
#' @param mixtureWeights Optional numeric vector of mixing proportions, one per
#'   panel in \code{genotypeHandle} when it is a list. Must be non-negative and
#'   sum to 1. Required whenever \code{genotypeHandle} is a list and downstream
#'   code will call \code{getCorrelation()}.
#' @return An \code{LdData} object.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:8]
#' gr <- GenomicRanges::GRanges("22",
#'   IRanges::IRanges(seq(1L, by = 100L, length.out = 8), width = 1L))
#' ld <- LdData(correlation = cor(X), variants = gr,
#'   blockMetadata = S4Vectors::DataFrame(
#'     chrom = "22", start = 1L, end = 1000L))
#' ld
#' @export
LdData <- function(
    correlation = NULL,
    genotypeHandle = NULL,
    snpIdx = NULL,
    variants,
    blockMetadata,
    nRef = 0L,
    mixtureWeights = NULL
) {
    obj <- new(
        "LdData",
        correlation = correlation,
        genotypeHandle = genotypeHandle,
        snpIdx = snpIdx,
        variants = variants,
        blockMetadata = blockMetadata,
        nRef = as.integer(nRef),
        mixtureWeights = mixtureWeights
    )
    validObject(obj)
    obj
}

# Internal: convert a refPanel data.frame (chrom/pos/A1/A2/variant_id, with
# optional allele_freq/variance/n_nomiss) into the GRanges form used by the
# LdData `variants` slot.
.refPanelToGranges <- function(refPanel) {
    chr <- withChrPrefix(refPanel$chrom)
    pos <- as.integer(refPanel$pos)

    gr <- GRanges(
        seqnames = chr,
        ranges = IRanges(start = pos, width = 1L)
    )

    mcolsData <- DataFrame(
        variant_id = refPanel$variant_id,
        A1 = refPanel$A1,
        A2 = refPanel$A2
    )

    optional <- c("allele_freq", "variance", "n_nomiss")
    for (col in optional) {
        if (is_in(col, names(refPanel))) {
            mcolsData[[col]] <- refPanel[[col]]
        }
    }
    mcols(gr) <- mcolsData
    gr
}

#' @rdname getCorrelation
#' @export
setMethod("getCorrelation", "LdData", function(x) {
    if (!is.null(x@correlation)) {
        return(x@correlation)
    }
    if (is.null(x@genotypeHandle)) {
        abort("No correlation matrix or genotype handle available")
    }
    if (is.list(x@genotypeHandle)) {
        if (is.null(x@mixtureWeights)) {
            msg <- glue(
                "Cannot compute mixture LD: `mixtureWeights` is NULL. ",
                "Construct LdData with mixtureWeights = <numeric vector> ",
                "when supplying a list of GenotypeHandles."
            )
            abort(msg)
        }
        perPanel <- map(x@genotypeHandle, .ldPanelLd, snpIdx = x@snpIdx)
        dims <- map_int(perPanel, nrow)
        if (length(unique(dims)) != 1L) {
            msg <- glue(
                "Mixture panels yielded LD matrices of differing ",
                "dimensions: {str_flatten(dims, ', ')}. All panels ",
                "must be aligned on the same variant subset."
            )
            abort(msg)
        }
        w <- x@mixtureWeights
        R <- matrix(0, nrow = dims[[1L]], ncol = dims[[1L]])
        for (k in seq_along(perPanel)) {
            R <- R + w[[k]] * perPanel[[k]]
        }
        dimnames(R) <- dimnames(perPanel[[1L]])
        return(R)
    }
    computeLd(.dosageMatrix(x@genotypeHandle, x@snpIdx), method = "sample")
})

#' @rdname getGenotypes
#' @export
setMethod("getGenotypes", "LdData", function(x, ...) {
    if (is.null(x@genotypeHandle)) {
        return(NULL)
    }
    if (is.matrix(x@genotypeHandle)) {
        return(x@genotypeHandle)
    }
    if (is.list(x@genotypeHandle)) {
        map(x@genotypeHandle, .dosageMatrix, x@snpIdx)
    } else {
        .dosageMatrix(x@genotypeHandle, x@snpIdx)
    }
})

#' @rdname hasGenotypes
#' @export
setMethod("hasGenotypes", "LdData", function(x) {
    !is.null(getGenotypeHandle(x))
})

#' @rdname getVariantIds
#' @export
setMethod("getVariantIds", "LdData", function(x, ...) {
    mcols(x@variants)$variant_id
})

#' @rdname getVariantInfo
#' @export
setMethod("getVariantInfo", "LdData", function(x) {
    x@variants
})

#' @rdname getBlockMetadata
#' @export
setMethod("getBlockMetadata", "LdData", function(x) {
    x@blockMetadata
})

#' @rdname getRefPanel
#' @export
setMethod("getRefPanel", "LdData", function(x) {
    mc <- as_tibble(as.data.frame(mcols(x@variants)))
    mc$chrom <- as.character(seqnames(x@variants))
    mc$pos <- start(x@variants)
    mc
})

#' @rdname getGenotypeHandle
#' @export
setMethod("getGenotypeHandle", "LdData", function(x) x@genotypeHandle)

#' @rdname getMixtureWeights
#' @export
setMethod("getMixtureWeights", "LdData", function(x) x@mixtureWeights)

#' @rdname getSnpIdx
#' @export
setMethod("getSnpIdx", "LdData", function(x) x@snpIdx)

#' @rdname getNRef
#' @export
setMethod("getNRef", "LdData", function(x) x@nRef)

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# Sample-LD matrix for one mixture panel's genotype handle (over `snpIdx`).
# @noRd
.ldPanelLd <- function(h, snpIdx) {
    computeLd(.dosageMatrix(h, snpIdx), method = "sample")
}
