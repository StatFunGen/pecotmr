#' (Deprecated) Univariate Analysis Pipeline
#'
#' \strong{Deprecated.} Use \code{\link{fineMappingPipeline}} with
#' \code{methods = "susie"} (and \code{methods = c("susie", "susieInf")}
#' for the SuSiE-inf + SuSiE initialization). The new pipeline accepts
#' a \code{\link{QtlDataset}} (individual-level) and dispatches on
#' \code{contexts} / \code{traitId} / \code{region}.
#'
#' For the TWAS-weights portion of the old pipeline, use
#' \code{\link{twasWeightsPipeline}} on the same \code{QtlDataset}.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
univariateAnalysisPipeline <- function(...) {
  .Deprecated(new = "fineMappingPipeline", package = "pecotmr",
    msg = paste(
      "univariateAnalysisPipeline() has been removed. Use",
      "fineMappingPipeline(qtlDataset, methods = 'susie') for",
      "fine-mapping and twasWeightsPipeline(qtlDataset, ...) for",
      "TWAS weights."))
  invisible(NULL)
}

#' (Deprecated) Load LD for a study, supporting single or mixture panels
#'
#' \strong{Deprecated.} Mixture-LD support is deferred from the new
#' \code{\link{GenotypeHandle}}-based LD-sketch design. For single-panel
#' studies, build a \code{GenotypeHandle} (or load via
#' \code{\link{loadLdMatrix}}) and pass it as the \code{ldSketch} slot of
#' your \code{QtlSumStats} / \code{GwasSumStats}. Mixture-LD support
#' will be reintroduced in a follow-up.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
loadStudyLd <- function(...) {
  .Deprecated(new = "GenotypeHandle", package = "pecotmr",
    msg = paste(
      "loadStudyLd() has been removed. For single-panel LD, build a",
      "GenotypeHandle (or use loadLdMatrix()) and pass it as the",
      "ldSketch slot of your SumStats. Mixture-LD support is deferred",
      "and will be reintroduced in a follow-up."))
  invisible(NULL)
}

#' (Deprecated) RSS Analysis Pipeline
#'
#' \strong{Deprecated.} Use \code{\link{fineMappingPipeline}} with
#' \code{methods = "susieRSS"}. The new pipeline accepts a
#' \code{\link{QtlSumStats}} (per-trait RSS) or a
#' \code{\link{GwasSumStats}} (per-LD-block GWAS RSS) and runs
#' \code{\link{summaryStatsQc}}-checked harmonization before
#' fine-mapping; reject inputs where \code{getQcInfo(x)} is empty.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
rssAnalysisPipeline <- function(...) {
  .Deprecated(new = "fineMappingPipeline", package = "pecotmr",
    msg = paste(
      "rssAnalysisPipeline() has been removed. Use",
      "fineMappingPipeline(sumStats, methods = 'susieRSS') after",
      "running summaryStatsQc() on the SumStats input."))
  invisible(NULL)
}
