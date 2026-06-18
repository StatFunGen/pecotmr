#' (Deprecated) Multivariate Analysis Pipeline
#'
#' \strong{Deprecated.} Use \code{\link{fineMappingPipeline}} with
#' \code{methods = "mvsusie"} (and \code{methods = "mvsusieRSS"} for
#' the RSS variant). The new pipeline accepts a
#' \code{\link{MultiTaskQtlDataset}} (individual-level studies) or a
#' \code{\link{QtlSumStats}} (summary-statistic studies) and handles
#' the joint multi-trait / multi-context dispatch internally.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
multivariateAnalysisPipeline <- function(...) {
  .Deprecated(new = "fineMappingPipeline", package = "pecotmr",
    msg = paste(
      "multivariateAnalysisPipeline() has been removed. Use",
      "fineMappingPipeline(data, methods = 'mvsusie') for",
      "individual-level joint analyses or",
      "fineMappingPipeline(data, methods = 'mvsusieRSS') for",
      "RSS-based joint analyses."))
  invisible(NULL)
}
