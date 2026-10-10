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

#' @title cTWAS Per-Run Payload
#' @description Per-run cTWAS payload: the fine-mapping posterior table, the
#'   full per-effect susie alpha table, the jointly-estimated group prior(s),
#'   and per-region metadata. One entry sits in every row of a
#'   \code{\linkS4class{CtwasResult}} collection.
#' @slot posteriors The per-gene (and, when SNPs are retained,
#'   per-SNP) posterior
#'   summary table (\code{ctwas::finemap_regions} \code{finemap_res} shape),
#'   or \code{NULL}.
#' @slot susieAlpha The per-effect susie alpha table
#'   (\code{ctwas::finemap_regions} \code{susie_alpha_res} shape) -- the
#'   fuller cTWAS output retained so the raw run is reconstructable, or
#'   \code{NULL}.
#' @slot groupPriors The estimated \code{group_prior} /
#'   \code{group_prior_var} for
#'   this run, or \code{NULL}.
#' @slot regionInfo Per-region metadata, or \code{NULL}.
#' @seealso \code{\link{CtwasResultEntry}} for the constructor.
#' @export
setClass(
    "CtwasResultEntry",
    representation(
        posteriors = "ANY", # per-gene/SNP posterior summary (ctwas finemap_res)
        # per-effect susie alpha table (ctwas susie_alpha_res)
        susieAlpha = "ANY",
        groupPriors = "ANY", # group_prior / group_prior_var for this run
        regionInfo = "ANY" # per-region metadata (optional)
    ),
    prototype = prototype(
        posteriors = NULL,
        susieAlpha = NULL,
        groupPriors = NULL,
        regionInfo = NULL
    )
)

#' @title Create a CtwasResultEntry
#' @description Per-run cTWAS payload wrapping the fine-mapping posterior table,
#'   the full per-effect susie alpha table, the estimated group prior(s), and
#'   region metadata. Held in every row of a \code{\link{CtwasResult}}
#'   collection.
#' @param posteriors The per-gene (and, when SNPs are retained,
#'   per-SNP) posterior
#'   summary table (\code{ctwas::finemap_regions} \code{finemap_res} shape), or
#'   \code{NULL}.
#' @param susieAlpha The per-effect susie alpha table
#'   (\code{ctwas::finemap_regions} \code{susie_alpha_res} shape) -- the fuller
#'   cTWAS output retained so the raw run is reconstructable, or \code{NULL}.
#' @param groupPriors The estimated \code{group_prior} /
#'   \code{group_prior_var} for
#'   this run, or \code{NULL}.
#' @param regionInfo Per-region metadata, or \code{NULL}.
#' @return A \code{CtwasResultEntry} object.
#' @examples
#' cre <- CtwasResultEntry(
#'   posteriors = data.frame(id = c("g1", "g2"), susie_pip = c(0.9, 0.1)),
#'   susieAlpha = data.frame(id = c("g1", "g2"), alpha = c(0.9, 0.1)))
#' cre
#' @export
CtwasResultEntry <- function(
    posteriors = NULL,
    susieAlpha = NULL,
    groupPriors = NULL,
    regionInfo = NULL
) {
    new(
        "CtwasResultEntry",
        posteriors = posteriors,
        susieAlpha = susieAlpha,
        groupPriors = groupPriors,
        regionInfo = regionInfo
    )
}

#' @rdname ctwasPosteriors
#' @export
setMethod("ctwasPosteriors", "CtwasResultEntry", function(x) x@posteriors)

#' @rdname susieAlpha
#' @export
setMethod("susieAlpha", "CtwasResultEntry", function(x) x@susieAlpha)

#' @rdname ctwasGroupPriors
#' @export
setMethod(
    "ctwasGroupPriors",
    "CtwasResultEntry",
    function(x) x@groupPriors
)
