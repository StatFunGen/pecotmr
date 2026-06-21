# =============================================================================
# H2Estimate S4 class
# -----------------------------------------------------------------------------
# Container for univariate heritability estimation results: global h2,
# optional per-block local estimates, optional annotation-stratified
# enrichment with jackknife blocks, and method-specific score statistics.
# Produced by `estimateH2()` and consumed by `h2EstimateToSldscTrait()` for
# integration with the sLDSC postprocessing pipeline.
# =============================================================================

#' @title Heritability Estimate
#' @description Container for univariate heritability estimation results.
#'   Holds global, local, and annotation-stratified estimates.
#' @slot h2 Numeric, global SNP heritability estimate.
#' @slot h2Se Numeric, standard error of global h2.
#' @slot intercept Numeric, confounding intercept estimate (NA if method
#'   does not estimate one).
#' @slot interceptSe Numeric, SE of intercept.
#' @slot local A \code{data.frame} with per-block local heritability
#'   estimates (columns: \code{blockId}, \code{h2Local}, \code{h2LocalSe}).
#'   NULL if \code{local = FALSE}.
#' @slot enrichment A \code{data.frame} with baseline annotation enrichment
#'   estimates (columns: \code{annotation}, \code{tau}, \code{tauSe},
#'   \code{enrichment}, \code{enrichmentSe}, \code{enrichmentP},
#'   \code{propH2}, \code{propSnps}). NULL if unstratified.
#' @slot tauBlocks A numeric matrix (nBlocks x n_annotations) of per-block
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
setClass("H2Estimate",
  representation(
    h2 = "numeric",
    h2Se = "numeric",
    intercept = "numeric",
    interceptSe = "numeric",
    local = "ANY",        # data.frame or NULL
    enrichment = "ANY",   # data.frame or NULL
    tauBlocks = "ANY",    # matrix or NULL
    scoreStats = "ANY",   # list or NULL
    method = "character",
    nSnps = "integer",
    traitName = "character"
  )
)


# =============================================================================
# Accessors
# =============================================================================

#' @rdname getH2
#' @export
setMethod("getH2", "H2Estimate", function(x) x@h2)

#' @rdname getTauBlocks
#' @export
setMethod("getTauBlocks", "H2Estimate", function(x) x@tauBlocks)

#' @rdname getLocal
#' @export
setMethod("getLocal", "H2Estimate", function(object) {
  object@local
})

#' @rdname getEnrichment
#' @export
setMethod("getEnrichment", "H2Estimate", function(object) {
  object@enrichment
})

#' @rdname getScoreStats
#' @export
setMethod("getScoreStats", "H2Estimate", function(object) {
  object@scoreStats
})


# =============================================================================
# Show
# =============================================================================

#' @export
setMethod("show", "H2Estimate", function(object) {
  cat(sprintf("H2Estimate for '%s' (method: %s)\n",
              object@traitName, object@method))
  cat(sprintf("  h2 = %.4f (SE = %.4f)\n", object@h2, object@h2Se))
  if (!is.na(object@intercept))
    cat(sprintf("  intercept = %.4f (SE = %.4f)\n",
                object@intercept, object@interceptSe))
  has_local <- !is.null(object@local)
  has_enrich <- !is.null(object@enrichment)
  has_tau_blocks <- !is.null(object@tauBlocks)
  cat(sprintf("  Local: %s, Enrichment: %s, tauBlocks: %s\n",
              has_local, has_enrich, has_tau_blocks))
  cat(sprintf("  N SNPs: %d\n", object@nSnps))
})
