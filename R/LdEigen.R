# =============================================================================
# LdEigen S4 class
# -----------------------------------------------------------------------------
# Per-block eigendecompositions of the LD correlation matrix: stored for
# LdMatrix input, computed from the genotype panel one block at a time for
# panel input. Consumed by LDER / HDL / sHDL h2 estimators.
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
#'
#'   An \code{LdEigen} built from a genotype panel stores no
#'   decompositions: each block's is computed from the panel's factor (its
#'   standardized genotypes) when an estimator reaches the block. That gives
#'   the block's nonzero eigenvalues and their eigenvectors; when the panel
#'   has fewer samples than the block has variants, the block's remaining
#'   eigenvalues are exactly zero, and the estimators account for those
#'   directions without listing them.
#' @slot eigenList A list of length \code{nBlocks} (empty for an
#'   \code{LdEigen} built from a genotype panel), each element a list with
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
#' @importFrom checkmate makeAssertCollection assertList
.validateLdEigen <- function(object) {
    coll <- makeAssertCollection()
    parentCheck <- .validateLdStatistic(object)
    if (!isTRUE(parentCheck)) {
        coll$push(parentCheck)
    }
    # One decomposition per block, or none for an LdEigen that computes them
    # from its genotype panel.
    assertList(
        object@eigenList,
        len = if (is.null(object@ldReference)) length(object@ldBlocks) else 0L,
        .var.name = "eigenList",
        add = coll
    )
    # checkmate has no exclusive lower bound, so (0, 1] stays a plain check.
    if (
        length(object@eigenvalueTruncation) != 1L ||
            object@eigenvalueTruncation <= 0 ||
            object@eigenvalueTruncation > 1
    ) {
        coll$push("'eigenvalueTruncation' must be a single value in (0, 1]")
    }
    coll$getMessages()
}

#' @title Create an LdEigen
#' @description Bundle per-block eigendecompositions with the variants they
#'   were computed over.
#' @param snpInfo A \code{data.frame} with columns \code{SNP}, \code{CHR},
#'   \code{BP}, \code{A1}, \code{A2} (and optionally \code{MAF}).
#' @param eigenList A list with one entry per LD block.
#' @param ldBlocks The LD block intervals: a \code{GRanges}, or the path
#'   of a BED file, read as BED with a 0-based start.
#' @param nSamples Integer, sample size of the LD reference panel.
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
#'   ldBlocks = blocks, nSamples = 100L, genome = "hg19")
#' length(le)
#' length(eigenList(le))
#' @importFrom checkmate assertDataFrame assertList assertCount
#' @export
LdEigen <- function(
    snpInfo,
    eigenList,
    ldBlocks,
    nSamples,
    inSample = FALSE,
    genome = NA_character_,
    eigenvalueTruncation = 1
) {
    assertDataFrame(snpInfo)
    assertList(eigenList)
    assertCount(nSamples, positive = TRUE)
    methods::new(
        "LdEigen",
        .ldStatRanges(snpInfo, genome),
        ldBlocks = .asLdBlockRanges(ldBlocks),
        nSamples = as.integer(nSamples),
        inSample = isTRUE(inSample),
        eigenList = eigenList,
        eigenvalueTruncation = as.numeric(eigenvalueTruncation)
    )
}

#' @importFrom checkmate assertCount assertFlag
#' @title Build an LdEigen from an LD reference
#' @description Eigendecompose an LD reference, block by block, into the
#'   \code{LdEigen} that \code{\link{estimateH2}} consumes for
#'   \code{method = "lder"} and \code{method = "hdl"}. This is the supported
#'   route from an LD reference to an h2 input: \code{\link{LdEigen}} itself
#'   is the low-level constructor and expects the decompositions to have been
#'   computed already.
#'
#' @param ldReference The LD reference, block by block: a list of
#'   \code{\link{LdMatrix}} objects, one per LD block (a single
#'   \code{LdMatrix} is one block), or a genotype panel (see
#'   \code{\link{readGenotypes}}) together with \code{ldBlocks}. Blocks are
#'   kept in the order given, and that order defines the variant order of the
#'   result.
#'
#'   Block structure matters downstream: \code{\link{estimateH2}} takes its
#'   standard error from a delete-one-block jackknife, so it requires at
#'   least two blocks and wants many more.
#' @param ldBlocks With a genotype panel, the LD blocks to cut it into: a
#'   \code{GRanges}, or the path of a BED file, read as BED with a 0-based
#'   start. Each block's LD is computed from the panel variants
#'   inside it; variants outside every block are left out. Not used with
#'   \code{LdMatrix} blocks.
#' @param nSamples Integer, the number of samples the reference's LD was
#'   estimated from. Defaults to a panel's sample count (a stochastic
#'   sketch's effective size), or the size the \code{LdMatrix} blocks record;
#'   required when they record none, and an error when they disagree.
#' @param inSample Logical, whether the reference is the GWAS cohort itself.
#'   Selects LDER's weighting and HDL's finite-reference correction.
#' @param genome Character, genome build. Defaults to the build recorded by
#'   the reference's variants.
#' @param eigenvalueTruncation Numeric in (0, 1]; the proportion of each
#'   block's eigenvalue mass to retain. \code{1} (the default) keeps every
#'   component; HDL conventionally uses \code{0.9}.
#' @return An \code{LdEigen} over every variant in the reference's blocks.
#'   One built from a genotype panel keeps the panel and decomposes each
#'   block when it is needed, rather than storing the decompositions.
#' @seealso \code{\link{buildLdScore}} for the g-LDSC input,
#'   \code{\link{estimateH2}}, \code{\link{LdMatrix}}
#' @examples
#' panel <- readGenotypes(plink2Prefix = file.path(system.file("extdata",
#'   "ld_reference", "chr22", package = "pecotmr"),
#'   "protocol_example.LD.chr22"))
#' blocks <- GenomicRanges::GRanges("chr22",
#'   IRanges::IRanges(c(10000000, 15000001), c(15000000, 19000000)))
#' ldEigen <- buildLdEigen(panel, ldBlocks = blocks, genome = "hg38")
#' ldEigen
#' length(eigenList(ldEigen))
#' @export
buildLdEigen <- function(
    ldReference,
    ldBlocks = NULL,
    nSamples = NULL,
    inSample = FALSE,
    genome = NA_character_,
    eigenvalueTruncation = 1
) {
    assertCount(nSamples, positive = TRUE, null.ok = TRUE)
    assertFlag(inSample)
    prep <- .h2RefPrepare(ldReference, ldBlocks, nSamples, genome)
    if (!is.null(prep$ldReference)) {
        return(.ldEigenFromPanel(prep, inSample, eigenvalueTruncation))
    }
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
        nSamples = prep$nRef,
        inSample = inSample,
        genome = prep$genome,
        eigenvalueTruncation = eigenvalueTruncation
    )
}

# An LdEigen over a genotype panel: the panel and each variant's row in it,
# and no decompositions, which are computed a block at a time when needed.
# Built in one step because an empty eigenList is valid only alongside the
# panel it is read from.
# @noRd
.ldEigenFromPanel <- function(prep, inSample, truncation) {
    gr <- .ldStatRanges(prep$snpInfo, prep$genome)
    S4Vectors::mcols(gr)$panelRow <- prep$panelRow
    methods::new(
        "LdEigen",
        gr,
        ldBlocks = .asLdBlockRanges(prep$ldBlocks),
        nSamples = as.integer(prep$nRef),
        inSample = isTRUE(inSample),
        eigenList = list(),
        eigenvalueTruncation = as.numeric(truncation),
        ldReference = prep$ldReference
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
    .ldEigenRefuseSubset()
})

# head(), tail(), window() and subset() reach the rows through extractROWS()
# without passing `[`, so the refusal is made there as well.
setMethod("extractROWS", "LdEigen", function(x, i) {
    .ldEigenRefuseSubset()
})

# @noRd
.ldEigenRefuseSubset <- function() {
    abort(glue(
        "an LdEigen cannot be subset: its eigendecompositions are per LD ",
        "block, so narrowing the variants would leave them describing ",
        "variants that are no longer present. Recompute over the subset."
    ))
}

#' @rdname show-methods
#' @export
setMethod("show", "LdEigen", function(object) {
    cat(glue(
        "LdEigen: {length(object)} SNPs across ",
        "{length(object@ldBlocks)} blocks\n",
        .trim = FALSE
    ))
    cat(sprintf("  Eigenvalue truncation: %.2f\n", object@eigenvalueTruncation))
    if (!is.null(object@ldReference)) {
        cat("  Decompositions: computed from the genotype panel when needed\n")
    }
    cat(glue(
        "  Reference N: {object@nSamples}, In-sample: {object@inSample}\n",
        .trim = FALSE
    ))
})

# The block records of an LdEigen (see .ldStatLazyBlocks()).
# @noRd
.ldEigenBlocks <- function(x) {
    if (length(x@eigenList) > 0L) {
        return(map(x@eigenList, .ldEigenStoredBlock))
    }
    .ldStatLazyBlocks(x)
}

# @noRd
.ldEigenStoredBlock <- function(block) {
    list(idx = block$snpIdx, eigen = block)
}

# One block's decomposition, at the truncation the LdEigen records.
# @noRd
.ldEigenBlockAt <- function(block, x) {
    .ldStatBlockEigen(block, x@ldReference, x@eigenvalueTruncation)
}

#' @rdname eigenList
#' @export
setMethod("eigenList", "LdEigen", function(x) {
    if (is.null(x@ldReference)) {
        return(x@eigenList)
    }
    map(.ldEigenBlocks(x), .ldEigenBlockAt, x = x)
})
