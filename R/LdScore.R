# =============================================================================
# LdScore S4 class
# -----------------------------------------------------------------------------
# Pre-computed LD scores (sum of ldsc's unbiased r^2) per SNP. Consumed by
# S-LDSC and g-LDSC. g-LDSC's FGLS residual covariance also needs each block's
# LD: stored per block for LdMatrix input, read from the genotype panel on
# demand for panel input.
# =============================================================================

#' @include LdStatistic.R
NULL

#' @title LD Score-Based LD Statistic
#' @description Pre-computed LD scores for each SNP. Used by S-LDSC and
#'   g-LDSC. Supports both standard LD scores and annotation-stratified ones.
#'
#'   The scores and their regression weights are \code{mcols} rather than
#'   slots, because they are parallel to the variants: as slots they would
#'   survive a subset unchanged while the ranges narrowed, leaving scores
#'   describing variants that are no longer there.
#'
#'   An \code{LdScore} storing per-block LD matrices (\code{ldMatrixList})
#'   cannot be subset or reordered: each block's matrix indexes its variants
#'   by position, so narrowing the variants would leave it describing other
#'   variants. A scores-only object, or one reading its LD from a genotype
#'   panel, subsets freely.
#' @slot ldMatrixList For g-LDSC: a list of per-block LD correlation (R)
#'   matrices used to compute the FGLS residual covariance, kept when the
#'   scores were built from \code{LdMatrix} blocks. Empty for a scores-only
#'   object, and for one built from a genotype panel, which reads each
#'   block's LD from the panel instead (see \code{\linkS4class{LdStatistic}}).
#' @export
setClass(
    "LdScore",
    contains = "LdStatistic",
    representation(ldMatrixList = "list"),
    validity = function(object) .validateLdScore(object)
)

# @noRd
#' @importFrom checkmate makeAssertCollection assertNames
.validateLdScore <- function(object) {
    coll <- makeAssertCollection()
    parentCheck <- .validateLdStatistic(object)
    if (!isTRUE(parentCheck)) {
        coll$push(parentCheck)
    }
    assertNames(
        colnames(S4Vectors::mcols(object, use.names = FALSE)),
        must.include = c("ldScores", "ldScoreWeights"),
        what = "colnames",
        .var.name = "mcols",
        add = coll
    )
    coll$getMessages()
}

#' @title Create an LdScore
#' @description Bundle pre-computed LD scores with the variants they describe.
#' @param snpInfo A \code{data.frame} with columns \code{SNP}, \code{CHR},
#'   \code{BP}, \code{A1}, \code{A2} (and optionally \code{MAF}).
#' @param ldScores A numeric matrix, one row per variant. The first column is
#'   the base LD score (sum of r^2, which \code{\link{buildLdScore}} takes as
#'   ldsc's unbiased r^2, as ldsc's own score files do); further columns are
#'   annotation-stratified scores.
#' @param ldScoreWeights Numeric regression weights, one per variant.
#' @param ldBlocks The LD block intervals: a \code{GRanges}, or the path
#'   of a BED file, read as BED with a 0-based start.
#' @param nSamples Integer, sample size of the LD reference panel.
#' @param inSample Logical, whether the reference is the GWAS cohort.
#' @param genome Character, genome build; recorded in \code{seqinfo()}.
#' @param ldMatrixList Optional list of per-block LD matrices (g-LDSC only).
#' @return An \code{LdScore}.
#' @examples
#' snpInfo <- data.frame(SNP = paste0("rs", 1:4), CHR = "chr1",
#'   BP = c(50L, 150L, 250L, 350L), A1 = "A", A2 = "G")
#' blocks <- GenomicRanges::GRanges("chr1",
#'   IRanges::IRanges(c(1L, 200L), c(199L, 400L)))
#' ls <- LdScore(snpInfo = snpInfo,
#'   ldScores = matrix(runif(4), ncol = 1, dimnames = list(NULL, "base_l2")),
#'   ldScoreWeights = rep(1, 4), ldBlocks = blocks, nSamples = 100L,
#'   inSample = FALSE, genome = "hg19")
#' length(ls)
#' head(S4Vectors::mcols(ls)$ldScores)
#' @importFrom checkmate assertMatrix assertNumeric
#' @export
LdScore <- function(
    snpInfo,
    ldScores,
    ldScoreWeights,
    ldBlocks,
    nSamples,
    inSample = FALSE,
    genome = NA_character_,
    ldMatrixList = list()
) {
    gr <- .ldStatRanges(snpInfo, genome)
    ldScores <- as.matrix(ldScores)
    assertMatrix(ldScores, nrows = length(gr))
    assertNumeric(ldScoreWeights, len = length(gr))
    scored <- S4Vectors::`mcols<-`(
        gr,
        value = cbind(
            S4Vectors::mcols(gr, use.names = FALSE),
            S4Vectors::DataFrame(
                ldScores = I(ldScores),
                ldScoreWeights = as.numeric(ldScoreWeights)
            )
        )
    )
    methods::new(
        "LdScore",
        scored,
        ldBlocks = .asLdBlockRanges(ldBlocks),
        nSamples = as.integer(nSamples),
        inSample = isTRUE(inSample),
        ldMatrixList = ldMatrixList
    )
}

#' @importFrom checkmate assertCount assertFlag
#' @title Build an LdScore from an LD reference
#' @description Compute per-variant LD scores from an LD reference, block by
#'   block, into the \code{LdScore} that \code{\link{estimateH2}} consumes
#'   for \code{method = "gldsc"}. This is the supported route from an LD
#'   reference to an h2 input: \code{\link{LdScore}} itself is the low-level
#'   constructor and expects the scores to have been computed already.
#'
#'   The score is \eqn{\ell_j = \sum_k r^2_{jk}} within each block, with
#'   each term ldsc's unbiased estimate of \eqn{r^2} (Bulik-Sullivan et al.
#'   2015, \code{ldscore.py}): \code{r^2 - (1 - r^2) / (n - 2)}, where
#'   \code{n} is the reference panel sample size. A sample \eqn{r^2}
#'   overstates the population value by about \code{(1 - r^2) / (n - 2)},
#'   and a block sums many of them. A variant's own term stays 1, and
#'   negative estimates are kept, as in ldsc. \code{\link{computeLdScores}}
#'   gives the same scores from an \code{\link{LdEigen}} built on the same
#'   reference.
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
#'   required when they record none, and an error when they disagree. It is
#'   also the \code{n} of the unbiased \eqn{r^2}, and it must exceed 2.
#' @param inSample Logical, whether the reference is the GWAS cohort itself.
#' @param genome Character, genome build. Defaults to the build recorded by
#'   the reference's variants.
#' @param ldScoreWeights Optional numeric regression weights, one per variant.
#'   Defaults to the conventional \eqn{1/\max(\ell_j, 1)}.
#' @param keepLdMatrices Logical, for \code{LdMatrix} blocks: keep each
#'   block's LD matrix on the result. \code{TRUE} (the default) is required
#'   for \code{method = "gldsc"} and for annotation-stratified scores, which
#'   need each block's LD. \code{FALSE} drops them, leaving a scores-only
#'   object. A genotype panel is never copied into matrices: the result keeps
#'   the panel and reads each block's LD from it when it is needed, so it
#'   serves every method either way.
#' @return An \code{LdScore} over every variant in the reference's blocks.
#' @seealso \code{\link{buildLdEigen}} for the LDER / HDL input,
#'   \code{\link{estimateH2}}, \code{\link{computeLdScores}}
#' @examples
#' panel <- readGenotypes(plink2Prefix = file.path(system.file("extdata",
#'   "ld_reference", "chr22", package = "pecotmr"),
#'   "protocol_example.LD.chr22"))
#' blocks <- GenomicRanges::GRanges("chr22",
#'   IRanges::IRanges(c(10000000, 15000001), c(15000000, 19000000)))
#' ldScore <- buildLdScore(panel, ldBlocks = blocks, genome = "hg38")
#' ldScore
#' head(S4Vectors::mcols(ldScore)$ldScores)
#' @export
buildLdScore <- function(
    ldReference,
    ldBlocks = NULL,
    nSamples = NULL,
    inSample = FALSE,
    genome = NA_character_,
    ldScoreWeights = NULL,
    keepLdMatrices = TRUE
) {
    assertCount(nSamples, positive = TRUE, null.ok = TRUE)
    assertFlag(inSample)
    assertFlag(keepLdMatrices)
    prep <- .h2RefPrepare(ldReference, ldBlocks, nSamples, genome)
    l2 <- .ldScoreVector(prep, .ldscAssertN(prep$nRef))
    ldMatrixList <- if (isTRUE(keepLdMatrices) && is.null(prep$ldReference)) {
        map2(prep$blocks, prep$snpIdx, .ldScoreKeepMatrix)
    } else {
        list()
    }
    scores <- LdScore(
        snpInfo = prep$snpInfo,
        ldScores = matrix(l2, ncol = 1, dimnames = list(NULL, "base_l2")),
        ldScoreWeights = .ldScoreResolveWeights(ldScoreWeights, l2),
        ldBlocks = prep$ldBlocks,
        nSamples = prep$nRef,
        inSample = inSample,
        genome = prep$genome,
        ldMatrixList = ldMatrixList
    )
    .ldStatAttachReference(scores, prep)
}

# Per-variant sum of ldsc's unbiased r^2 within the variant's own block, in
# reference order. The blocks partition the variants, which are stacked
# block by block, so the blocks' sums concatenate into place; a genotype
# panel's blocks are read one at a time.
# @noRd
.ldScoreVector <- function(prep, nRef) {
    blocks <- map2(prep$blocks, prep$snpIdx, .ldScoreBuildBlock)
    list_c(map(
        blocks,
        .ldScoreBlockSums,
        ldReference = prep$ldReference,
        nRef = nRef
    ))
}

# A block of the reference being built, as the record the block readers take.
# @noRd
.ldScoreBuildBlock <- function(block, snpIdx) {
    list(idx = snpIdx, R = block$R, rows = block$rows)
}

# One block's base LD scores: the raw sums of r^2 over the block, corrected to
# ldsc's unbiased r^2 by its summed form.
# @noRd
.ldScoreBlockSums <- function(block, ldReference, nRef) {
    ld <- .ldStatBlockLd(block, ldReference)
    p <- .ldSourceSize(ld)
    raw <- .ldSquaredProduct(ld, matrix(1, p, 1L))
    as.vector(.ldscUnbiasedSum(raw, p, nRef))
}

# ldsc's unbiased estimator of r^2 (ldscore.py, __l2_unbiased__): a sample
# r^2 overstates the population value by about (1 - r^2) / (n - 2). The self
# term (r^2 = 1) stays 1, and negative estimates are kept, as in ldsc. Every
# LD score pecotmr computes goes through this or its summed form below.
# @noRd
.ldscUnbiasedR2 <- function(r2, nRef) {
    r2 - (1 - r2) / (nRef - 2)
}

# .ldscUnbiasedR2() summed over a block with weights a_k:
#   sum_k a_k (r2_jk - (1 - r2_jk) / (n - 2))
#     = ((n - 1) * sum_k a_k r2_jk - sum_k a_k) / (n - 2).
# Exact, so a route that has the raw weighted sums but not r^2 itself -- an
# eigendecomposition -- corrects them without rebuilding the matrix.
# `rawSums` is variants x columns; `weightSums` holds each column's sum_k a_k
# over the block, which for the unweighted base score is the block size.
# @noRd
.ldscUnbiasedSum <- function(rawSums, weightSums, nRef) {
    sweep((nRef - 1) * rawSums, 2L, weightSums) / (nRef - 2)
}

# ldsc's unbiased r^2 divides by n - 2, so the panel must exceed 2 samples.
# @noRd
.ldscAssertN <- function(nRef) {
    if (nRef > 2L) {
        return(nRef)
    }
    abort(glue(
        "ldsc's unbiased r^2 divides by n - 2, so the LD reference panel ",
        "needs more than 2 samples, and it records {nRef}. Give the panel's ",
        "size as `nSamples`."
    ))
}

# The n of ldsc's unbiased r^2 for a built LD reference: its panel size.
# @noRd
.ldscPanelN <- function(ldRef) {
    .ldscAssertN(ldPanelNSamples(ldRef))
}

# @noRd
.ldScoreKeepMatrix <- function(block, snpIdx) {
    list(R = block$R, snpIdx = as.integer(snpIdx))
}

# The conventional LDSC heteroskedasticity weight, with the score floored at
# 1 so a variant in near-perfect linkage equilibrium cannot dominate.
# @noRd
.ldScoreResolveWeights <- function(ldScoreWeights, l2) {
    if (is.null(ldScoreWeights)) {
        return(1 / pmax(l2, 1))
    }
    assertNumeric(ldScoreWeights, len = length(l2))
    as.numeric(ldScoreWeights)
}

#' @rdname show-methods
#' @export
setMethod("show", "LdScore", function(object) {
    scores <- mcols(object)$ldScores
    cat(glue(
        "LdScore: {length(object)} SNPs, ",
        "{ncol(scores)} LD score columns\n",
        .trim = FALSE
    ))
    cat(glue("  Block LD: {.ldScoreLdSource(object)}\n", .trim = FALSE))
    cat(glue(
        "  Reference N: {object@nSamples}, In-sample: {object@inSample}\n",
        .trim = FALSE
    ))
})

# Where an LdScore's per-block LD comes from, as show() describes it.
# @noRd
.ldScoreLdSource <- function(x) {
    if (!is.null(x@ldReference)) {
        return("read from the genotype panel when needed")
    }
    if (length(x@ldMatrixList) > 0L) {
        return("stored matrices")
    }
    "none (scores only; g-LDSC and stratified scores need it)"
}

# TRUE when each block's LD can be had: stored, or read from the panel.
# @noRd
.ldScoreHasLd <- function(x) {
    !is.null(x@ldReference) || length(x@ldMatrixList) > 0L
}

# The block records of an LdScore (see .ldStatLazyBlocks()).
# @noRd
.ldScoreBlocks <- function(x) {
    if (length(x@ldMatrixList) > 0L) {
        return(map(x@ldMatrixList, .ldScoreStoredBlock))
    }
    .ldStatLazyBlocks(x)
}

# @noRd
.ldScoreStoredBlock <- function(block) {
    list(idx = block$snpIdx, R = block$R)
}

# ---- subsetting -------------------------------------------------------------

# Stored block matrices index their variants by position (`snpIdx`), so a
# subset or reorder would leave them describing other variants. Without them
# each block's variants are found afresh from the ranges, so the subset is
# safe. Every subset -- `[`, head(), tail(), window(), subset() -- reaches the
# rows through extractROWS().
setMethod("extractROWS", "LdScore", function(x, i) {
    if (length(x@ldMatrixList) > 0L) {
        abort(glue(
            "an LdScore storing per-block LD matrices cannot be subset: each ",
            "block's matrix indexes its variants by position, so narrowing ",
            "or reordering the variants would leave it describing other ",
            "variants. Rebuild the scores over the subset instead."
        ))
    }
    methods::callNextMethod()
})

#' @rdname ldMatrixList
#' @export
setMethod("ldMatrixList", "LdScore", function(x) {
    if (is.null(x@ldReference)) {
        return(x@ldMatrixList)
    }
    map(.ldScoreBlocks(x), .ldScorePanelMatrix, ldReference = x@ldReference)
})

# One block's LD matrix computed from the panel an LdScore reads it from.
# @noRd
.ldScorePanelMatrix <- function(block, ldReference) {
    ldFactor <- .ldRefFactor(ldReference, block$rows)
    list(
        R = crossprod(ldFactor) / (nrow(ldFactor) - 1),
        snpIdx = block$idx
    )
}
