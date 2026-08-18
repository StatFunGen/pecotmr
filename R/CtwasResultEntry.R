# =============================================================================
# CtwasResultEntry S4 class
# -----------------------------------------------------------------------------
# Per-run cTWAS payload: the per-gene (+ per-SNP) posterior table from
# ctwas::ctwas_sumstats, the jointly-estimated group prior + prior variance,
# and per-region metadata. One entry sits in every row of a CtwasResult
# collection.
# =============================================================================

#' @include AllGenerics.R
NULL

setClass(
    "CtwasResultEntry",
    representation(
        finemap = "ANY", # per-gene/SNP posterior summary (ctwas finemap_res)
        # per-effect susie alpha table (ctwas susie_alpha_res)
        susieAlpha = "ANY",
        param = "ANY", # group_prior / group_prior_var for this run
        regionInfo = "ANY" # per-region metadata (optional)
    ),
    prototype = prototype(
        finemap = NULL,
        susieAlpha = NULL,
        param = NULL,
        regionInfo = NULL
    )
)

#' @title Create a CtwasResultEntry
#' @description Per-run cTWAS payload wrapping the fine-mapping posterior table,
#'   the full per-effect susie alpha table, the estimated group prior(s), and
#'   region metadata. Held in every row of a \code{\link{CtwasResult}}
#'   collection.
#' @param finemap The per-gene (and, when SNPs are retained, per-SNP) posterior
#'   summary table (\code{ctwas::finemap_regions} \code{finemap_res} shape), or
#'   \code{NULL}.
#' @param susieAlpha The per-effect susie alpha table
#'   (\code{ctwas::finemap_regions} \code{susie_alpha_res} shape) -- the fuller
#'   cTWAS output retained so the raw run is reconstructable, or \code{NULL}.
#' @param param The estimated \code{group_prior} / \code{group_prior_var} for
#'   this run, or \code{NULL}.
#' @param regionInfo Per-region metadata, or \code{NULL}.
#' @return A \code{CtwasResultEntry} object.
#' @examples
#' cre <- CtwasResultEntry(
#'   finemap = data.frame(id = c("g1", "g2"), susie_pip = c(0.9, 0.1)),
#'   susieAlpha = data.frame(id = c("g1", "g2"), alpha = c(0.9, 0.1)))
#' cre
#' @export
CtwasResultEntry <- function(
    finemap = NULL,
    susieAlpha = NULL,
    param = NULL,
    regionInfo = NULL
) {
    new(
        "CtwasResultEntry",
        finemap = finemap,
        susieAlpha = susieAlpha,
        param = param,
        regionInfo = regionInfo
    )
}

#' @rdname getFinemap
#' @export
setMethod("getFinemap", "CtwasResultEntry", function(x, ...) x@finemap)

#' @rdname getSusieAlpha
#' @export
setMethod("getSusieAlpha", "CtwasResultEntry", function(x, ...) x@susieAlpha)

#' @rdname getCtwasParam
#' @export
setMethod("getCtwasParam", "CtwasResultEntry", function(x, ...) x@param)
