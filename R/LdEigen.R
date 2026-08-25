# =============================================================================
# LdEigen S4 class
# -----------------------------------------------------------------------------
# Pre-computed per-block eigendecompositions of the LD correlation matrix.
# Consumed by LDER / HDL / sHDL h2 estimators.
# =============================================================================

#' @include LdStatistic.R
NULL

#' @title Eigendecomposition-Based LD Statistic
#' @description Pre-computed per-block eigendecompositions of the LD
#'   correlation matrix. Used by LDER, HDL and sHDL.
#'
#'   Unlike \code{\linkS4class{LdScore}}, whose statistics are per variant,
#'   an eigendecomposition describes a whole LD \emph{block}. It therefore
#'   stays in a slot rather than \code{mcols} -- and because a slot does not
#'   narrow when the ranges do, subsetting an \code{LdEigen} is refused
#'   rather than allowed to produce an object whose decompositions no longer
#'   describe its variants.
#' @slot eigenList A list of length \code{nBlocks}, each element a list with
#'   components:
#'   \describe{
#'     \item{values}{Numeric vector of eigenvalues}
#'     \item{vectors}{Numeric matrix of eigenvectors (SNPs x retained
#'     components)}
#'     \item{snpIdx}{Integer vector of variant indices}
#'   }
#' @slot eigenvalueTruncation Numeric, proportion of variance retained (e.g.
#'   0.9 for HDL's default). If 1.0, no truncation.
#' @export
setClass(
    "LdEigen",
    contains = "LdStatistic",
    representation(
        eigenList = "list",
        eigenvalueTruncation = "numeric"
    ),
    validity = function(object) .validateLdEigen(object)
)

# @noRd
.validateLdEigen <- function(object) {
    parentCheck <- .validateLdStatistic(object)
    errors <- if (isTRUE(parentCheck)) character() else parentCheck
    if (length(object@eigenList) != length(object@ldBlocks)) {
        errors <- c(
            errors,
            "Length of 'eigenList' must match number of LD blocks"
        )
    }
    if (
        length(object@eigenvalueTruncation) != 1L ||
            object@eigenvalueTruncation <= 0 ||
            object@eigenvalueTruncation > 1
    ) {
        errors <- c(
            errors,
            "'eigenvalueTruncation' must be a single value in (0, 1]"
        )
    }
    if (length(errors) == 0) TRUE else errors
}

#' @title Create an LdEigen
#' @description Bundle per-block eigendecompositions with the variants they
#'   were computed over.
#' @param snpInfo A \code{data.frame} with columns \code{SNP}, \code{CHR},
#'   \code{BP}, \code{A1}, \code{A2} (and optionally \code{MAF}).
#' @param eigenList A list with one entry per LD block.
#' @param ldBlocks A \code{GRanges} of LD block intervals.
#' @param nRef Integer, sample size of the LD reference panel.
#' @param inSample Logical, whether the reference is the GWAS cohort.
#' @param genome Character, genome build; recorded in \code{seqinfo()}.
#' @param eigenvalueTruncation Numeric in (0, 1]; proportion of variance
#'   retained.
#' @return An \code{LdEigen}.
#' @examples
#' snpInfo <- data.frame(SNP = paste0("rs", 1:4), CHR = "chr1",
#'   BP = c(50L, 150L, 250L, 350L), A1 = "A", A2 = "G")
#' blocks <- GenomicRanges::GRanges("chr1",
#'   IRanges::IRanges(c(1L, 200L), c(199L, 400L)))
#' mkBlock <- function(idx) list(values = rep(1, length(idx)),
#'   vectors = diag(length(idx)), snpIdx = idx)
#' le <- LdEigen(snpInfo = snpInfo,
#'   eigenList = list(mkBlock(1:2), mkBlock(3:4)),
#'   ldBlocks = blocks, nRef = 100L, genome = "hg19")
#' length(le)
#' length(getEigenList(le))
#' @export
LdEigen <- function(
    snpInfo,
    eigenList,
    ldBlocks,
    nRef,
    inSample = FALSE,
    genome = NA_character_,
    eigenvalueTruncation = 1
) {
    obj <- methods::new(
        "LdEigen",
        .ldStatRanges(snpInfo, genome),
        ldBlocks = .asLdBlockRanges(ldBlocks),
        nRef = as.integer(nRef),
        inSample = isTRUE(inSample),
        eigenList = eigenList,
        eigenvalueTruncation = as.numeric(eigenvalueTruncation)
    )
    validObject(obj)
    obj
}

#' @describeIn LdEigen-class Refused. Subsetting would narrow the variants
#'   while \code{eigenList} -- a slot, because it is per block rather than
#'   per variant -- stayed as it was, leaving decompositions that describe
#'   variants the object no longer has. Recompute over the subset instead.
#' @param x An \code{LdEigen}.
#' @param i,j,... Subscripts; any use is an error.
#' @param drop Ignored.
#' @return Nothing: this method always signals an error.
#' @export
setMethod("[", "LdEigen", function(x, i, j, ..., drop = TRUE) {
    abort(glue(
        "an LdEigen cannot be subset: its eigendecompositions are per LD ",
        "block, so narrowing the variants would leave them describing ",
        "variants that are no longer present. Recompute over the subset."
    ))
})

#' @rdname show-methods
#' @export
setMethod("show", "LdEigen", function(object) {
    cat(glue(
        "LdEigen: {length(object)} SNPs across ",
        "{length(object@eigenList)} blocks\n",
        .trim = FALSE
    ))
    cat(sprintf("  Eigenvalue truncation: %.2f\n", object@eigenvalueTruncation))
    cat(glue(
        "  Reference N: {object@nRef}, In-sample: {object@inSample}\n",
        .trim = FALSE
    ))
})

#' @rdname getEigenList
#' @export
setMethod("getEigenList", "LdEigen", function(x) x@eigenList)
