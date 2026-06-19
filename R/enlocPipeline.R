#' @title Enrichment-Informed Colocalization Pipeline (enloc)
#' @description Variant of \code{\link{colocPipeline}} that uses per-pair
#'   QTL-GWAS enrichment estimates (from \code{\link{qtlEnrichmentPipeline}})
#'   to inform the colocalization priors. For each (QTL tuple, GWAS
#'   tuple) pair the per-variant shared-signal prior \code{p12} is
#'   scaled by the enrichment factor for that (gwasStudy, qtlContext)
#'   pair before being passed to \code{coloc::coloc.bf_bf}; the single-
#'   trait priors \code{p1} and \code{p2} remain at the supplied
#'   baseline.
#'
#' @section Prior adjustment:
#' For a baseline \code{p12} and enrichment factor \code{r}, the adjusted
#' prior is \code{p12 * (1 + r)} (capped at \code{p12Max}). When the
#' enrichment table has no row for a given (gwasStudy, qtlContext) pair
#' the baseline \code{p12} is used unchanged and a warning is emitted.
#'
#' @section Why coloc.bf_bf and not coloc.susie:
#' The earlier coloc.susie-driven implementation depended on rebuilding a
#' full SuSiE fit object from a \code{FineMappingEntry}; that path proved
#' brittle because the SuSiE post-processing trimming \code{pecotmr}
#' applies (drop-by-V, fSuSiE-shape handling, missing CS index) does not
#' round-trip cleanly to what \code{coloc.susie} expects. The
#' \code{coloc.bf_bf} path used by \code{\link{colocPipeline}} consumes
#' the already-trimmed LBF matrix directly and is the supported entry
#' point for both pipelines.
#'
#' @section LD-sketch identity check:
#' Same as \code{\link{colocPipeline}}.
#'
#' @param qtlFineMappingResult A \code{\link{QtlFineMappingResult}}
#'   (required).
#' @param gwasInput Either a \code{\link{GwasSumStats}} or a
#'   \code{\link{GwasFineMappingResult}}.
#' @param enrichment Data frame from
#'   \code{\link{qtlEnrichmentPipeline}} with at least \code{gwasStudy},
#'   \code{qtlContext}, and \code{enrichment} columns.
#' @param filterLbfCs,filterLbfCsSecondary,priorTol LBF filtering knobs
#'   shared with \code{\link{colocPipeline}}.
#' @param p1,p2,p12 Baseline priors for \code{coloc::coloc.bf_bf}.
#'   Default \code{1e-4}, \code{1e-4}, \code{5e-6}.
#' @param p12Max Maximum allowed value for the enrichment-adjusted
#'   \code{p12} prior. Default \code{1e-3}.
#' @param finemappingMethods See \code{\link{colocPipeline}}.
#' @param returnGwasFineMapping See \code{\link{colocPipeline}}.
#' @param ... Forwarded to \code{coloc::coloc.bf_bf}.
#' @return Data frame as in \code{\link{colocPipeline}}, with additional
#'   \code{enrichment} and \code{p12Used} columns recording the per-pair
#'   factor and the prior actually passed to \code{coloc::coloc.bf_bf}.
#' @export
enlocPipeline <- function(qtlFineMappingResult,
                          gwasInput,
                          enrichment,
                          filterLbfCs           = FALSE,
                          filterLbfCsSecondary  = NULL,
                          priorTol              = 1e-9,
                          p1 = 1e-4, p2 = 1e-4, p12 = 5e-6,
                          p12Max = 1e-3,
                          finemappingMethods = "susie",
                          returnGwasFineMapping = FALSE,
                          ...) {
  if (!requireNamespace("coloc", quietly = TRUE)) {
    stop("Package 'coloc' is required for enlocPipeline. ",
         "Install with: install.packages('coloc').")
  }
  if (!methods::is(qtlFineMappingResult, "QtlFineMappingResult")) {
    stop("`qtlFineMappingResult` must be a QtlFineMappingResult ",
         "(got class '", class(qtlFineMappingResult)[[1L]], "').")
  }
  if (!methods::is(gwasInput, "GwasSumStats") &&
      !methods::is(gwasInput, "GwasFineMappingResult")) {
    stop("`gwasInput` must be a GwasSumStats or a ",
         "GwasFineMappingResult (got class '",
         class(gwasInput)[[1L]], "').")
  }
  if (missing(enrichment) || is.null(enrichment) ||
      !is.data.frame(enrichment))
    stop("`enrichment` must be a data.frame with at least gwasStudy, ",
         "qtlContext, enrichment columns (output of ",
         "qtlEnrichmentPipeline).")
  required <- c("gwasStudy", "qtlContext", "enrichment")
  missingCols <- setdiff(required, colnames(enrichment))
  if (length(missingCols) > 0L)
    stop("`enrichment` is missing column(s): ",
         paste(missingCols, collapse = ", "))

  # Resolve gwas fine-mapping (shared with colocPipeline contract).
  gwasFmr <- if (methods::is(gwasInput, "GwasFineMappingResult")) {
    gwasInput
  } else {
    if (length(getQcInfo(gwasInput)) == 0L)
      stop("enlocPipeline: gwasInput (GwasSumStats) has no QC record. ",
           "Call summaryStatsQc() first.")
    fineMappingPipeline(gwasInput, methods = finemappingMethods)
  }
  .colocRequireMatchingLdSketches(getLdSketch(qtlFineMappingResult),
                                  getLdSketch(gwasFmr))

  # Pre-extract per-GWAS-tuple LBF matrices (same path colocPipeline uses).
  gwasLbfByPair <- .colocPreextractGwasLbf(
    gwasFmr, filterLbfCs, filterLbfCsSecondary, priorTol)
  if (length(gwasLbfByPair) == 0L) {
    out <- .enlocEmptyResult()
    if (returnGwasFineMapping && methods::is(gwasInput, "GwasSumStats")) {
      attr(out, "gwasFineMapping") <- gwasFmr
    }
    return(out)
  }

  results <- list()
  for (qi in seq_len(nrow(qtlFineMappingResult))) {
    qStudy   <- as.character(qtlFineMappingResult$study)[[qi]]
    qContext <- as.character(qtlFineMappingResult$context)[[qi]]
    qTrait   <- as.character(qtlFineMappingResult$trait)[[qi]]
    qMethod  <- as.character(qtlFineMappingResult$method)[[qi]]
    qEntry   <- qtlFineMappingResult$entry[[qi]]
    qLbfInfo <- .colocExtractLbfFromEntry(
      qEntry, filterLbfCs, filterLbfCsSecondary, priorTol,
      label = sprintf("QTL (study='%s', context='%s', trait='%s', method='%s')",
                      qStudy, qContext, qTrait, qMethod))
    if (is.null(qLbfInfo)) next
    qLbf <- qLbfInfo$lbf

    for (gKey in names(gwasLbfByPair)) {
      gInfo  <- gwasLbfByPair[[gKey]]
      gLbf   <- gInfo$lbf
      gStudy <- gInfo$study
      gMethod <- gInfo$method

      aligned <- .colocAlignLbf(qLbf, gLbf)
      if (is.null(aligned)) next

      enRow <- .enlocLookupEnrichment(enrichment, gStudy, qContext)
      if (is.na(enRow)) {
        warning(sprintf(
          "enlocPipeline: no enrichment entry for (gwasStudy='%s', qtlContext='%s'); using baseline p12.",
          gStudy, qContext))
        enRow <- 0
      }
      p12Used <- min(p12 * (1 + enRow), p12Max)

      pairRes <- tryCatch(
        coloc::coloc.bf_bf(aligned$qtl, aligned$gwas,
                           p1 = p1, p2 = p2, p12 = p12Used, ...),
        error = function(e) {
          warning(sprintf(
            "enlocPipeline: coloc.bf_bf failed for QTL (study='%s', context='%s', trait='%s', method='%s') x GWAS (study='%s', method='%s'): %s",
            qStudy, qContext, qTrait, qMethod,
            gStudy, gMethod, conditionMessage(e)))
          NULL
        })
      if (is.null(pairRes) || is.null(pairRes$summary)) next

      sm <- as.data.frame(pairRes$summary, stringsAsFactors = FALSE)
      sm$qtlStudy   <- qStudy
      sm$qtlContext <- qContext
      sm$qtlTrait   <- qTrait
      sm$qtlMethod  <- qMethod
      sm$gwasStudy  <- gStudy
      sm$gwasMethod <- gMethod
      sm$enrichment <- enRow
      sm$p12Used    <- p12Used
      results[[length(results) + 1L]] <- sm
    }
  }

  if (length(results) == 0L) {
    out <- .enlocEmptyResult()
  } else {
    out <- do.call(rbind, lapply(results, .colocStandardiseRow))
    rownames(out) <- NULL
    idCols <- c("qtlStudy", "qtlContext", "qtlTrait", "qtlMethod",
                "gwasStudy", "gwasMethod", "enrichment", "p12Used")
    other  <- setdiff(colnames(out), idCols)
    out    <- out[, c(idCols, other), drop = FALSE]
  }

  if (returnGwasFineMapping && methods::is(gwasInput, "GwasSumStats")) {
    attr(out, "gwasFineMapping") <- gwasFmr
  }
  out
}

# Internal: look up the enrichment factor for a (gwasStudy, qtlContext)
# pair. Returns NA when the pair is not in the table; the caller emits
# the warning and falls back to baseline.
# @noRd
.enlocLookupEnrichment <- function(enrichment, gwasStudy, qtlContext) {
  idx <- which(as.character(enrichment$gwasStudy)  == gwasStudy &
               as.character(enrichment$qtlContext) == qtlContext)
  if (length(idx) == 0L) return(NA_real_)
  as.numeric(enrichment$enrichment[[idx[[1L]]]])
}

# Empty-result frame for the no-pair case. Mirrors .colocEmptyResult
# but with the enloc-specific extra columns.
# @noRd
.enlocEmptyResult <- function() {
  base <- .colocEmptyResult()
  base$enrichment <- numeric(0)
  base$p12Used    <- numeric(0)
  base
}
