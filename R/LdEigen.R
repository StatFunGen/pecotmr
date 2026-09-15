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

#' @title Build an LdEigen from loaded LD
#' @description Eigendecompose already-loaded LD, block by block, into the
#'   \code{LdEigen} that \code{\link{estimateH2}} consumes for
#'   \code{method = "lder"} and \code{method = "hdl"}. This is the supported
#'   route from an LD reference to an h2 input: \code{\link{LdEigen}} itself
#'   is the low-level constructor and expects the decompositions to have been
#'   computed already.
#'
#'   No file I/O happens here. Read the LD first with
#'   \code{\link{loadLdMatrix}} -- passing a vector of regions returns one
#'   \code{LdData} per block -- then pass the result in.
#' @param ldBlockData An \code{\link{LdData}}, or a list of them, covering
#'   the reference variants. An \code{LdData} whose correlation is a list of
#'   per-block matrices contributes one LD block per matrix; one whose
#'   correlation is a single matrix contributes one block. Blocks are kept in
#'   the order given, and that order defines the variant order of the result.
#'
#'   Block structure matters downstream: \code{\link{estimateH2}} takes its
#'   standard error from a delete-one-block jackknife, so it requires at
#'   least two blocks and wants many more. Loading a whole region as one
#'   dense matrix yields a single block.
#' @param nRef Integer, the LD reference panel sample size. Defaults to the
#'   size recorded by the supplied \code{LdData}; required when they record
#'   none, and an error when they disagree.
#' @param inSample Logical, whether the reference is the GWAS cohort itself.
#'   Selects LDER's weighting and HDL's finite-reference correction.
#' @param genome Character, genome build. Defaults to the build recorded by
#'   the supplied \code{LdData}, which \code{\link{loadLdMatrix}} leaves
#'   unset.
#' @param eigenvalueTruncation Numeric in (0, 1]; the proportion of each
#'   block's eigenvalue mass to retain. \code{1} (the default) keeps every
#'   component; HDL conventionally uses \code{0.9}.
#' @return An \code{LdEigen} over every variant in \code{ldBlockData}.
#' @seealso \code{\link{buildLdScore}} for the g-LDSC input,
#'   \code{\link{estimateH2}}, \code{\link{loadLdMatrix}}
#' @examples
#' meta <- system.file("extdata", "ld_reference", "ld_meta_file.tsv",
#'   package = "pecotmr")
#' ld <- loadLdMatrix(meta, region = "chr22:10000000-19000000")
#' ldEigen <- buildLdEigen(ld, genome = "hg38")
#' ldEigen
#' length(getEigenList(ldEigen))
#' @export
buildLdEigen <- function(
    ldBlockData,
    nRef = NULL,
    inSample = FALSE,
    genome = NA_character_,
    eigenvalueTruncation = 1
) {
    prep <- .ldRefPrepare(ldBlockData, nRef, genome)
    eigenList <- map2(
        prep$blocks,
        prep$snpIdx,
        .ldEigenOneBlock,
        truncation = eigenvalueTruncation
    )
    LdEigen(
        snpInfo = prep$snpInfo,
        eigenList = eigenList,
        ldBlocks = prep$ldBlocks,
        nRef = prep$nRef,
        inSample = inSample,
        genome = prep$genome,
        eigenvalueTruncation = eigenvalueTruncation
    )
}

# One block's decomposition, truncated to the leading components carrying
# `truncation` of its eigenvalue mass.
# @noRd
.ldEigenOneBlock <- function(block, snpIdx, truncation) {
    e <- eigen(block$R, symmetric = TRUE)
    keep <- .ldEigenKeep(e$values, truncation)
    list(
        values = e$values[keep],
        vectors = e$vectors[, keep, drop = FALSE],
        snpIdx = as.integer(snpIdx)
    )
}

# Rounding can leave a near-singular block with slightly negative trailing
# eigenvalues, so the cumulative mass is taken over the clamped values --
# otherwise the running total can dip and pick the wrong cut point.
# @noRd
.ldEigenKeep <- function(values, truncation) {
    if (truncation >= 1) {
        return(seq_along(values))
    }
    positive <- pmax(values, 0)
    total <- sum(positive)
    if (total <= 0) {
        return(seq_along(values))
    }
    seq_len(which(cumsum(positive) / total >= truncation)[[1L]])
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
