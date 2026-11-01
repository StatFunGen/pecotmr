# =============================================================================
# H2Estimate S4 class
# -----------------------------------------------------------------------------
# Container for univariate heritability estimation results: global h2,
# optional per-block local estimates, optional annotation-stratified
# enrichment with jackknife blocks, and method-specific score statistics.
# Produced by `estimateH2()` and consumed by `h2EstimateToSldscTrait()` for
# integration with the sLDSC postprocessing pipeline.
# =============================================================================

#' @include MethodParam.R
NULL

#' @title Heritability Estimate
#' @description Container for univariate heritability estimation results. Holds
#'   global, local, and annotation-stratified estimates.
#' @slot estimate Numeric, global SNP heritability estimate.
#' @slot estimateSe Numeric, standard error of global h2.
#' @slot intercept Numeric, confounding intercept estimate (NA if method does
#'   not estimate one).
#' @slot interceptSe Numeric, SE of intercept.
#' @slot localBlocks A \code{data.frame} with per-block local
#'   heritability estimates
#'   (columns: \code{blockId}, \code{h2Local}, \code{h2LocalSe}). NULL if
#'   \code{local = FALSE}.
#' @slot enrichment A \code{data.frame} of per-annotation enrichment results,
#'   NULL if unstratified. Two shapes are possible, depending on which route
#'   produced the estimate:
#'   \describe{
#'     \item{partitioned}{enrichment of the user's baseline annotations
#'       from their partitioned heritabilities \eqn{A'A\tau}, with
#'       columns \code{annotation}, \code{tau}, \code{tauSe},
#'       \code{enrichment}, \code{enrichmentSe}, \code{enrichmentP},
#'       \code{propH2}, \code{propSnps}.}
#'     \item{score test}{the score test over the tested (candidate)
#'       annotations, with columns \code{annotation}, \code{scoreZ},
#'       \code{scoreP}.}
#'   }
#' @slot annotationJackknifeCoefs A numeric matrix
#'   (nBlocks x n_annotations) of per-block
#'   jackknife tau values. Required for Gazal tauStar standardization
#'   downstream. NULL if not available (e.g., unstratified analysis).
#' @slot scoreStats A list with score statistics for candidate annotations,
#'   suitable for input to \code{susieRss}. Contains:
#'   \describe{
#'     \item{z}{Numeric vector of z-scores for each candidate annotation}
#'     \item{R}{Correlation matrix of the score statistics}
#'     \item{annotationNames}{Character vector of candidate annotation names}
#'   }
#'   NULL if no candidate annotations provided.
#' @slot method Character string identifying the estimation method.
#' @slot nSnps Integer, number of SNPs used in estimation.
#' @slot traitName Character string for trait identifier.
#' @export
setClass(
    "H2Estimate",
    representation(
        estimate = "numeric",
        estimateSe = "numeric",
        intercept = "numeric",
        interceptSe = "numeric",
        localBlocks = "data.frame_OR_NULL",
        enrichment = "data.frame_OR_NULL",
        annotationJackknifeCoefs = "matrix_OR_NULL",
        scoreStats = "list_OR_NULL",
        method = "character",
        nSnps = "integer",
        traitName = "character"
    )
)


# =============================================================================
# Accessors
# =============================================================================

#' @rdname heritabilityEstimate
#' @export
setMethod("heritabilityEstimate", "H2Estimate", function(x) x@estimate)

#' @rdname heritabilityEstimateSe
#' @export
setMethod("heritabilityEstimateSe", "H2Estimate", function(x) x@estimateSe)

#' @rdname heritabilityIntercept
#' @export
setMethod("heritabilityIntercept", "H2Estimate", function(x) x@intercept)

#' @rdname heritabilityInterceptSe
#' @export
setMethod("heritabilityInterceptSe", "H2Estimate", function(x) x@interceptSe)

#' @rdname nSnps
#' @export
setMethod("nSnps", "H2Estimate", function(x) x@nSnps)

#' @rdname traitNames
#' @export
setMethod("traitNames", "H2Estimate", function(x) x@traitName)

#' @rdname methodNames
#' @export
setMethod("methodNames", "H2Estimate", function(x) x@method)

#' @rdname annotationJackknifeCoefs
#' @export
setMethod("annotationJackknifeCoefs", "H2Estimate", function(x) {
    x@annotationJackknifeCoefs
})

#' @rdname localH2Blocks
#' @export
setMethod("localH2Blocks", "H2Estimate", function(object) {
    object@localBlocks
})

#' @rdname stratifiedHeritabilityEnrichment
#' @export
setMethod("stratifiedHeritabilityEnrichment", "H2Estimate", function(object) {
    object@enrichment
})

#' @rdname scoreStats
#' @export
setMethod("scoreStats", "H2Estimate", function(object) {
    object@scoreStats
})

# Each slot an H2Estimate was saved under before its slots were renamed, and
# its name now.
# @noRd
.H2_RENAMED_SLOTS <- c(
    h2 = "estimate",
    h2Se = "estimateSe",
    local = "localBlocks",
    tauBlocks = "annotationJackknifeCoefs"
)

# An H2Estimate saved before the rename holds its estimates under the old
# names, which the inherited update would drop as slots the class no longer
# has, leaving the estimate empty; they are moved to the new names first.
#' @rdname updateObject-methods
#' @export
setMethod(
    "updateObject",
    "H2Estimate",
    function(object, ..., verbose = FALSE) {
        if (isTRUE(verbose)) {
            message("updateObject(object = 'H2Estimate')")
        }
        methods::callNextMethod(
            .renameSavedSlots(object, .H2_RENAMED_SLOTS),
            ...,
            verbose = verbose
        )
    }
)


# =============================================================================
# Show
# =============================================================================

#' @rdname show-methods
#' @export
setMethod("show", "H2Estimate", function(object) {
    cat(glue(
        "H2Estimate for '{object@traitName}' ",
        "(method: {object@method})\n",
        .trim = FALSE
    ))
    cat(sprintf(
        "  h2 = %.4f (SE = %.4f)\n",
        object@estimate,
        object@estimateSe
    ))
    if (!is.na(object@intercept)) {
        cat(sprintf(
            "  intercept = %.4f (SE = %.4f)\n",
            object@intercept,
            object@interceptSe
        ))
    }
    has_local <- !is.null(object@localBlocks)
    has_enrich <- !is.null(object@enrichment)
    has_tau_blocks <- !is.null(object@annotationJackknifeCoefs)
    cat(glue(
        "  Local: {has_local}, Enrichment: {has_enrich}, ",
        "annotationJackknifeCoefs: {has_tau_blocks}\n",
        .trim = FALSE
    ))
    cat(glue("  N SNPs: {object@nSnps}\n", .trim = FALSE))
})
