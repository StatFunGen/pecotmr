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
#' @section Coercion: \code{as(x, "GRanges")} gives the per-gene (and, with
#'   SNPs kept, per-SNP) posteriors as ranges -- an error when they could not
#'   be placed on the genome and are still a table.
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

#' @rdname CtwasResultEntry-class
#' @name coerce,CtwasResultEntry,GRanges-method
#' @aliases coerce,CtwasResultEntry,GRanges-method
#' @usage NULL
#' @exportMethod coerce
setAs("CtwasResultEntry", "GRanges", function(from) {
    .ctwasEntryRanges(from@posteriors)
})

# Posteriors as ranges: an absent table is no ranges; one the pipeline could
# not place on the genome has no ranges to give.
# @noRd
.ctwasEntryRanges <- function(posteriors) {
    if (is.null(posteriors)) {
        return(GenomicRanges::GRanges())
    }
    if (methods::is(posteriors, "GRanges")) {
        return(posteriors)
    }
    abort(str_c(
        "These cTWAS posteriors were not placed on the genome, so they have ",
        "no ranges; ctwasPosteriors(x) returns them as a table."
    ))
}

#' @rdname ctwasGroupPriors
#' @export
setMethod(
    "ctwasGroupPriors",
    "CtwasResultEntry",
    function(x) x@groupPriors
)

#' @rdname show-methods
#' @export
setMethod("show", "CtwasResultEntry", function(object) {
    post <- .ctwasPosteriorTable(object@posteriors)
    priors <- if (is.null(object@groupPriors)) "not set" else "set"
    cat(glue(
        "CtwasResultEntry: {.ctwasGeneCount(post)} gene(s) in ",
        "{.ctwasRegionCount(post)} region(s); group priors {priors}\n",
        .trim = FALSE
    ))
    invisible(NULL)
})

# The posterior table whichever form it is held in: a ranged payload's
# columns, or the table itself.
# @noRd
.ctwasPosteriorTable <- function(posteriors) {
    if (methods::is(posteriors, "GRanges")) {
        return(as.data.frame(mcols(posteriors)))
    }
    posteriors %||% data.frame()
}

# The genes a posterior table scores: its rows that are not SNPs.
# @noRd
.ctwasGeneCount <- function(post) {
    if (!is_in("id", colnames(post))) {
        return(0L)
    }
    genes <- if (is_in("type", colnames(post))) post$type != "SNP" else TRUE
    length(unique(as.character(post$id[genes])))
}

# @noRd
.ctwasRegionCount <- function(post) {
    if (!is_in("region_id", colnames(post))) {
        return(0L)
    }
    length(unique(as.character(post$region_id)))
}
