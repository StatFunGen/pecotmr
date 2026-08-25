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
.validateLdScore <- function(object) {
    parentCheck <- .validateLdStatistic(object)
    errors <- if (isTRUE(parentCheck)) character() else parentCheck
    md <- S4Vectors::mcols(object, use.names = FALSE)
    for (col in c("ldScores", "ldScoreWeights")) {
        if (!is_in(col, colnames(md))) {
            errors <- c(errors, glue("mcols must carry an '{col}' column"))
        }
    }
    if (length(errors) == 0) TRUE else errors
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
#' @param nRef Integer, sample size of the LD reference panel.
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
#'   ldScoreWeights = rep(1, 4), ldBlocks = blocks, nRef = 100L,
#'   inSample = FALSE, genome = "hg19")
#' length(ls)
#' head(getLdScores(ls))
#' @export
LdScore <- function(
    snpInfo,
    ldScores,
    ldScoreWeights,
    ldBlocks,
    nRef,
    inSample = FALSE,
    genome = NA_character_,
    ldMatrixList = list()
) {
    gr <- .ldStatRanges(snpInfo, genome)
    ldScores <- as.matrix(ldScores)
    if (nrow(ldScores) != length(gr)) {
        abort(glue(
            "`ldScores` has {nrow(ldScores)} row(s) for {length(gr)} ",
            "variant(s); they must be parallel."
        ))
    }
    if (length(ldScoreWeights) != length(gr)) {
        abort(glue(
            "`ldScoreWeights` has {length(ldScoreWeights)} value(s) for ",
            "{length(gr)} variant(s); they must be parallel."
        ))
    }
    md <- S4Vectors::mcols(gr, use.names = FALSE)
    md$ldScores <- ldScores
    md$ldScoreWeights <- as.numeric(ldScoreWeights)
    S4Vectors::mcols(gr) <- md
    obj <- methods::new(
        "LdScore",
        gr,
        ldBlocks = .asLdBlockRanges(ldBlocks),
        nRef = as.integer(nRef),
        inSample = isTRUE(inSample),
        ldMatrixList = ldMatrixList
    )
    validObject(obj)
    obj
}

#' @rdname show-methods
#' @export
setMethod("show", "LdScore", function(object) {
    scores <- getLdScores(object)
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
        "  Reference N: {object@nRef}, In-sample: {object@inSample}\n",
        .trim = FALSE
    ))
})

#' @rdname getLdScores
#' @export
setMethod("getLdScores", "LdScore", function(x) {
    S4Vectors::mcols(x, use.names = FALSE)$ldScores
})

#' @rdname getLdScoreWeights
#' @export
setMethod("getLdScoreWeights", "LdScore", function(x) {
    S4Vectors::mcols(x, use.names = FALSE)$ldScoreWeights
})

#' @rdname getLdMatrixList
#' @export
setMethod("getLdMatrixList", "LdScore", function(x) x@ldMatrixList)
