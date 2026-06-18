#' (Deprecated) xQTL GWAS Enrichment Analysis
#'
#' \strong{Deprecated.} The file-path enrichment wrapper has been
#' removed. Use \code{\link{qtlEnrichmentPipeline}}: build
#' genome-wide \code{\link{FineMappingResult}} collections for both
#' the GWAS (LD-block-indexed via \code{susieRSS}) and the QTLs in
#' your own code, then pass them to \code{qtlEnrichmentPipeline()}
#' along with their shared \code{ldSketch}.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
xqtlEnrichmentWrapper <- function(...) {
  .Deprecated(new = "qtlEnrichmentPipeline", package = "pecotmr",
    msg = paste(
      "xqtlEnrichmentWrapper() has been removed. Use",
      "qtlEnrichmentPipeline() with FineMappingResult collections for",
      "the GWAS and the QTLs."))
  invisible(NULL)
}

#' (Deprecated) Colocalization Wrapper
#'
#' \strong{Deprecated.} The file-path colocalization wrapper has been
#' removed. Use \code{\link{colocPipeline}}: build a
#' \code{FineMappingResult} for the QTLs in your own code, then pass
#' it along with a \code{GwasSumStats} (window-level GWAS
#' fine-mapping is run internally) or an existing genome-wide GWAS
#' \code{FineMappingResult} to \code{colocPipeline()}.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
colocWrapper <- function(...) {
  .Deprecated(new = "colocPipeline", package = "pecotmr",
    msg = paste(
      "colocWrapper() has been removed. Use colocPipeline() with a",
      "QTL FineMappingResult plus a GwasSumStats or GWAS",
      "FineMappingResult."))
  invisible(NULL)
}

#' (Deprecated) Colocalization Post-Processor
#'
#' \strong{Deprecated.} Post-processing now happens inside
#' \code{\link{colocPipeline}} and \code{\link{enlocPipeline}}, which
#' return the cleaned, ranked colocalization output directly.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
colocPostProcessor <- function(...) {
  .Deprecated(new = "colocPipeline", package = "pecotmr",
    msg = paste(
      "colocPostProcessor() has been removed. Post-processing is now",
      "internal to colocPipeline() / enlocPipeline()."))
  invisible(NULL)
}
