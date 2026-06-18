#' @title QTL Enrichment Pipeline (Genome-Wide)
#' @description Genome-wide pipeline that computes per-pair (GWAS study,
#'   QTL context) enrichment estimates by passing the GWAS PIP vector
#'   and the QTL credible-set posteriors to
#'   \code{\link{computeQtlEnrichment}}. The returned table feeds
#'   \code{\link{enlocPipeline}}.
#'
#'   \strong{Not gene-parallelisable}: the enrichment estimator runs
#'   over the full genome of GWAS PIPs and the full collection of QTL
#'   fits at once.
#'
#' @section Inputs:
#' \itemize{
#'   \item \code{gwasFineMappingResult}: a genome-wide
#'     \code{\link{GwasFineMappingResult}} (one row per (study, LD
#'     block) tuple). Each entry's \code{FineMappingEntry$trimmedFit}
#'     must carry a \code{pip} vector.
#'   \item \code{qtlFineMappingResult}: the genome-wide
#'     \code{\link{QtlFineMappingResult}}. Each entry's
#'     \code{trimmedFit} must carry \code{alpha}, \code{pip}, and
#'     prior-variance fields (\code{V}).
#' }
#'
#' @section LD-sketch identity check:
#' The GWAS \code{\link{FineMappingResultBase}} must have a non-NULL
#' \code{ldSketch} (RSS-derived). If the QTL FMR also has a non-NULL
#' \code{ldSketch}, the two must match exactly. When the QTL FMR's
#' \code{ldSketch} is NULL (individual-level QTL fit), validation is
#' skipped on the QTL side.
#'
#' @param gwasFineMappingResult See above.
#' @param qtlFineMappingResult See above.
#' @param numGwas Number of GWAS variants used to estimate \code{piGwas}.
#'   When \code{NULL} (default) it is estimated from the data — bias
#'   warning applies if the input PIP vector is not genome-wide.
#' @param piQtl Per-variant prior of being a QTL causal variant.
#'   \code{NULL} (default) estimates from the data.
#' @param lambda Shrinkage parameter for the enrichment estimator.
#'   Default \code{1.0}.
#' @param impN Number of imputed samples used by the estimator. Default
#'   \code{25}.
#' @param numThreads Number of threads used by
#'   \code{computeQtlEnrichment}. Default \code{1}.
#' @param ... Additional arguments forwarded to
#'   \code{\link{computeQtlEnrichment}}.
#' @return A data frame with one row per (gwasStudy, qtlContext) pair
#'   and columns \code{gwasStudy}, \code{qtlContext},
#'   \code{enrichment}, \code{enrichmentSe}, \code{enrichmentLogOdds},
#'   plus any extras the underlying estimator emits. Suitable as the
#'   \code{enrichment} argument to \code{\link{enlocPipeline}}.
#' @export
qtlEnrichmentPipeline <- function(gwasFineMappingResult,
                                  qtlFineMappingResult,
                                  numGwas = NULL,
                                  piQtl = NULL,
                                  lambda = 1.0,
                                  impN = 25,
                                  numThreads = 1L,
                                  ...) {
  if (!methods::is(gwasFineMappingResult, "GwasFineMappingResult"))
    stop("`gwasFineMappingResult` must be a GwasFineMappingResult.")
  if (!methods::is(qtlFineMappingResult, "QtlFineMappingResult"))
    stop("`qtlFineMappingResult` must be a QtlFineMappingResult.")

  gwasLd <- getLdSketch(gwasFineMappingResult)
  if (is.null(gwasLd))
    stop("qtlEnrichmentPipeline: the GWAS FineMappingResult must have a ",
         "non-NULL ldSketch (it should be RSS-derived).")
  qtlLd <- getLdSketch(qtlFineMappingResult)
  .colocRequireMatchingLdSketches(qtlLd, gwasLd)

  # Per-study genome-wide GWAS PIP vector (named by variant id).
  gwasStudies <- unique(as.character(gwasFineMappingResult$study))
  qtlContexts <- unique(as.character(qtlFineMappingResult$context))

  if (length(gwasStudies) == 0L || length(qtlContexts) == 0L)
    stop("qtlEnrichmentPipeline: no (gwasStudy, qtlContext) pairs ",
         "to compute (one of the inputs has zero rows).")

  results <- list()
  for (gStudy in gwasStudies) {
    gwasPip <- .enrBuildGwasPipVector(gwasFineMappingResult, gStudy)
    if (length(gwasPip) == 0L) {
      warning(sprintf(
        "qtlEnrichmentPipeline: no usable PIPs for gwasStudy='%s'; skipping.",
        gStudy))
      next
    }
    for (qContext in qtlContexts) {
      qtlRegions <- .enrBuildQtlRegionsList(qtlFineMappingResult, qContext)
      if (length(qtlRegions) == 0L) {
        warning(sprintf(
          "qtlEnrichmentPipeline: no usable QTL regions for qtlContext='%s'; skipping.",
          qContext))
        next
      }
      enr <- tryCatch(
        computeQtlEnrichment(
          gwasPip          = gwasPip,
          susieQtlRegions  = qtlRegions,
          numGwas          = numGwas,
          piQtl            = piQtl,
          lambda           = lambda,
          impN             = impN,
          numThreads       = numThreads,
          ...),
        error = function(e) {
          warning(sprintf(
            "qtlEnrichmentPipeline: computeQtlEnrichment failed for (gwasStudy='%s', qtlContext='%s'): %s",
            gStudy, qContext, conditionMessage(e)))
          NULL
        })
      if (is.null(enr)) next
      row <- .enrFlattenEnrichment(enr)
      row$gwasStudy  <- gStudy
      row$qtlContext <- qContext
      results[[length(results) + 1L]] <- row
    }
  }

  if (length(results) == 0L) {
    return(data.frame(
      gwasStudy  = character(0),
      qtlContext = character(0),
      enrichment = numeric(0),
      enrichmentSe = numeric(0),
      enrichmentLogOdds = numeric(0),
      stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, lapply(results, as.data.frame,
                              stringsAsFactors = FALSE))
  rownames(out) <- NULL
  idCols <- c("gwasStudy", "qtlContext")
  other  <- setdiff(colnames(out), idCols)
  out[, c(idCols, other), drop = FALSE]
}

# =============================================================================
# Internal helpers
# =============================================================================

# Build a named GWAS PIP vector for one study. Walks every row of the
# GwasFineMappingResult tagged with that study, extracts the per-row
# pip from each FineMappingEntry, and concatenates with variant-id
# names. Errors if any single variant appears with conflicting PIP
# values across rows.
# @noRd
.enrBuildGwasPipVector <- function(gwasFmr, gStudy) {
  idx <- which(as.character(gwasFmr$study) == gStudy)
  if (length(idx) == 0L) return(numeric(0))
  pieces <- list()
  for (i in idx) {
    entry <- gwasFmr$entry[[i]]
    fit <- getTrimmedFit(entry)
    if (is.null(fit) || is.null(fit$pip)) next
    pip <- as.numeric(fit$pip)
    ids <- if (!is.null(names(fit$pip))) names(fit$pip)
           else getVariantIds(entry)
    if (length(ids) != length(pip)) next
    pieces[[length(pieces) + 1L]] <-
      stats::setNames(pip, as.character(ids))
  }
  if (length(pieces) == 0L) return(numeric(0))
  all <- unlist(pieces)
  uniqIds <- unique(names(all))
  if (length(uniqIds) != length(all)) {
    # Same variant id in multiple blocks. Verify the values agree.
    dedup <- tapply(all, names(all), function(v) {
      if (length(unique(round(v, 12))) > 1L) {
        stop("qtlEnrichmentPipeline: variant '", names(v)[[1L]],
             "' appears with conflicting PIPs across GWAS blocks for ",
             "study '", "?", "'; the GWAS fine-mapping must produce a ",
             "consistent PIP per variant.")
      }
      v[[1L]]
    })
    all <- as.numeric(dedup)
    names(all) <- names(dedup)
  }
  all
}

# Build the per-(qtl context) list of region fits in the shape that
# computeQtlEnrichment expects: list(d) where each d carries
# alpha, pip, prior_variance (V).
# @noRd
.enrBuildQtlRegionsList <- function(qtlFmr, qContext) {
  idx <- which(as.character(qtlFmr$context) == qContext)
  if (length(idx) == 0L) return(list())
  out <- list()
  for (i in idx) {
    entry <- qtlFmr$entry[[i]]
    fit <- getTrimmedFit(entry)
    if (is.null(fit) || is.null(fit$alpha) || is.null(fit$pip)) next
    pV <- if (!is.null(fit$V)) fit$V
          else if (!is.null(fit$prior_variance)) fit$prior_variance
          else NULL
    if (is.null(pV)) next
    if (is.null(names(fit$pip)))
      names(fit$pip) <- getVariantIds(entry)
    out[[length(out) + 1L]] <- list(
      alpha          = fit$alpha,
      pip            = fit$pip,
      prior_variance = pV)
  }
  out
}

# Coerce computeQtlEnrichment's variable-shape output into a single-row
# named list with the canonical columns the caller documents. The
# underlying estimator returns either a list with named numeric scalars
# (enrichment, enrichmentSe, enrichmentLogOdds, ...) or a matrix/df —
# this helper handles both.
# @noRd
.enrFlattenEnrichment <- function(enr) {
  if (is.list(enr) && is.null(dim(enr))) {
    pickScalar <- function(field) {
      v <- enr[[field]]
      if (is.null(v)) NA_real_
      else as.numeric(v[[1L]])
    }
    list(
      enrichment        = pickScalar("enrichment"),
      enrichmentSe      = pickScalar("enrichmentSe"),
      enrichmentLogOdds = pickScalar("enrichmentLogOdds"))
  } else if (is.matrix(enr) || is.data.frame(enr)) {
    df <- as.data.frame(enr, stringsAsFactors = FALSE)
    if (nrow(df) == 0L) {
      list(enrichment = NA_real_, enrichmentSe = NA_real_,
           enrichmentLogOdds = NA_real_)
    } else {
      list(
        enrichment        = .enrPickColumn(df, c("enrichment", "Enrichment")),
        enrichmentSe      = .enrPickColumn(df, c("enrichmentSe", "se", "stderr")),
        enrichmentLogOdds = .enrPickColumn(df, c("enrichmentLogOdds", "logOdds", "log_odds")))
    }
  } else {
    list(enrichment = NA_real_, enrichmentSe = NA_real_,
         enrichmentLogOdds = NA_real_)
  }
}

# @noRd
.enrPickColumn <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0L) return(NA_real_)
  as.numeric(df[[hit[[1L]]]][[1L]])
}
