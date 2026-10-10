# =============================================================================
# LdScore S4 class
# -----------------------------------------------------------------------------
# Pre-computed LD scores (sum of r^2) per SNP. Consumed by S-LDSC and
# g-LDSC. Holds the optional per-block LD matrices needed for g-LDSC's
# FGLS residual covariance.
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
#' @slot ldMatrixList For g-LDSC: a list of per-block LD (R^2) matrices used to
#'   compute the FGLS residual covariance. Empty for S-LDSC.
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
#'   the base LD score (sum of r^2); further columns are
#'   annotation-stratified scores.
#' @param ldScoreWeights Numeric regression weights, one per variant.
#' @param ldBlocks A \code{GRanges} of LD block intervals.
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
    obj <- methods::new(
        "LdScore",
        scored,
        ldBlocks = .asLdBlockRanges(ldBlocks),
        nSamples = as.integer(nSamples),
        inSample = isTRUE(inSample),
        ldMatrixList = ldMatrixList
    )
    validObject(obj)
    obj
}

#' @importFrom checkmate assertCount assertFlag
#' @title Build an LdScore from loaded LD
#' @description Compute per-variant LD scores from already-loaded LD, block by
#'   block, into the \code{LdScore} that \code{\link{estimateH2}} consumes
#'   for \code{method = "gldsc"}. This is the supported route from an LD
#'   reference to an h2 input: \code{\link{LdScore}} itself is the low-level
#'   constructor and expects the scores to have been computed already.
#'
#'   The score is \eqn{\ell_j = \sum_k r^2_{jk}} within each block, the same
#'   quantity \code{\link{computeLdScores}} reconstructs from an
#'   \code{\link{LdEigen}}, so the two agree on a shared reference.
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
#' @param nSamples Integer, the LD reference panel sample size. Defaults to the
#'   size recorded by the supplied \code{LdData}; required when they record
#'   none, and an error when they disagree.
#' @param inSample Logical, whether the reference is the GWAS cohort itself.
#' @param genome Character, genome build. Defaults to the build recorded by
#'   the supplied \code{LdData}, which \code{\link{loadLdMatrix}} leaves
#'   unset.
#' @param ldScoreWeights Optional numeric regression weights, one per variant.
#'   Defaults to the conventional \eqn{1/\max(\ell_j, 1)}.
#' @param keepLdMatrices Logical. Keep each block's LD matrix on the result.
#'   \code{TRUE} (the default) is required for \code{method = "gldsc"},
#'   which forms its FGLS residual covariance from them and errors without
#'   them. \code{FALSE} drops them, leaving a scores-only object.
#' @return An \code{LdScore} over every variant in \code{ldBlockData}.
#' @seealso \code{\link{buildLdEigen}} for the LDER / HDL input,
#'   \code{\link{estimateH2}}, \code{\link{computeLdScores}}
#' @examples
#' meta <- system.file("extdata", "ld_reference", "ld_meta_file.tsv",
#'   package = "pecotmr")
#' ld <- loadLdMatrix(meta, region = "chr22:10000000-19000000")
#' ldScore <- buildLdScore(ld, genome = "hg38")
#' ldScore
#' head(S4Vectors::mcols(ldScore)$ldScores)
#' @export
buildLdScore <- function(
    ldBlockData,
    nSamples = NULL,
    inSample = FALSE,
    genome = NA_character_,
    ldScoreWeights = NULL,
    keepLdMatrices = TRUE
) {
    assertCount(nSamples, positive = TRUE, null.ok = TRUE)
    assertFlag(inSample)
    assertFlag(keepLdMatrices)
    prep <- .ldRefPrepare(ldBlockData, nSamples, genome)
    l2 <- .ldScoreVector(prep$blocks, prep$snpIdx, nrow(prep$snpInfo))
    ldMatrixList <- if (isTRUE(keepLdMatrices)) {
        map2(prep$blocks, prep$snpIdx, .ldScoreKeepMatrix)
    } else {
        list()
    }
    LdScore(
        snpInfo = prep$snpInfo,
        ldScores = matrix(l2, ncol = 1, dimnames = list(NULL, "base_l2")),
        ldScoreWeights = .ldScoreResolveWeights(ldScoreWeights, l2),
        ldBlocks = prep$ldBlocks,
        nSamples = prep$nRef,
        inSample = inSample,
        genome = prep$genome,
        ldMatrixList = ldMatrixList
    )
}

# Per-variant sum of r^2 within the variant's own block, scattered back into
# reference order.
# @noRd
.ldScoreVector <- function(blocks, snpIdx, nVariants) {
    if (length(blocks) == 0L) {
        return(numeric(nVariants))
    }
    # The blocks partition the variants, so each variant's score lands in
    # exactly one place; variants in no block stay zero.
    replace(
        numeric(nVariants),
        list_c(snpIdx),
        list_c(map(blocks, .ldScoreBlockSums))
    )
}

# @noRd
.ldScoreBlockSums <- function(block) {
    rowSums(block$R^2)
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
    cat(glue(
        "  Full LD matrices: {length(object@ldMatrixList) > 0} ",
        "(needed for g-LDSC)\n",
        .trim = FALSE
    ))
    cat(glue(
        "  Reference N: {object@nSamples}, In-sample: {object@inSample}\n",
        .trim = FALSE
    ))
})

#' @rdname ldMatrixList
#' @export
setMethod("ldMatrixList", "LdScore", function(x) x@ldMatrixList)
