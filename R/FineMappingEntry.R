# =============================================================================
# FineMappingEntry S4 class
# -----------------------------------------------------------------------------
# Per-tuple fine-mapping payload: variant ids + a method-specific trimmed
# fit object + a long-format topLoci data.frame (per-variant PIP, CS
# membership, beta/se). One entry sits in every row of a
# FineMappingResult collection. Accessors read directly from the payload
# slots (no further lookups required).
# =============================================================================

#' @include allGenerics.R
NULL

setClass("FineMappingEntry",
  representation(
    variantIds = "character",
    trimmedFit = "ANY",
    topLoci    = "data.frame",
    sumstats   = "ANY"
  ),
  validity = function(object) {
    errors <- character()
    if (nrow(object@topLoci) > 0L) {
      required <- c("variant_id", "pip")
      missingCols <- setdiff(required, colnames(object@topLoci))
      if (length(missingCols) > 0L)
        errors <- c(errors, paste("topLoci missing columns:",
                                  paste(missingCols, collapse = ", ")))
    }
    if (length(errors) == 0L) TRUE else errors
  }
)

#' @title TWAS Weights Entry (per-tuple payload)
#' @description S4 container for one method's TWAS weights, attached to
#'   a \code{TwasWeights} row. One entry corresponds to one
#'   \code{(study, context, trait, method)} tuple.
#' @slot variantIds Character vector of variant IDs that have weights.
#' @slot weights Numeric vector (single-method, single-outcome) or
#'   matrix (multi-outcome).
#' @slot fits Optional method-specific fit object.
#' @slot cvPerformance Optional named list of CV metrics (\code{rsq},
#'   \code{pval}, etc.).
#' @slot standardized Logical (length 1). Whether the weights are on the
#'   standardized scale.
#' @slot dataType Data-type tag for downstream usage (e.g.,
#'   \code{"expression"}, \code{"splicing"}); may be \code{NULL}.
#' @export


#' @title Create a FineMappingEntry Object
#' @description Construct a \code{FineMappingEntry} payload for one
#'   \code{(study, context, trait, method)} row of a
#'   \code{FineMappingResult} collection.
#' @param variantIds Character vector of variant IDs.
#' @param trimmedFit Method-specific fit object.
#' @param topLoci Long-format \code{data.frame}.
#' @param sumstats Optional list of summary statistics, or \code{NULL}.
#' @return A \code{FineMappingEntry} object.
#' @export
FineMappingEntry <- function(variantIds, trimmedFit, topLoci,
                             sumstats = NULL) {
  obj <- new("FineMappingEntry",
             variantIds = as.character(variantIds),
             trimmedFit = trimmedFit,
             topLoci    = as.data.frame(topLoci),
             sumstats   = sumstats)
  validObject(obj)
  obj
}


# Per-entry accessors (reuse the existing generics; these methods read
# slots from the payload classes directly).

#' @rdname getVariantIds
#' @export
setMethod("getVariantIds", "FineMappingEntry",
          function(x, ...) x@variantIds)

#' @rdname getTrimmedFit
#' @export
setMethod("getTrimmedFit", "FineMappingEntry",
          function(x, ...) x@trimmedFit)

#' @rdname getTopLoci
#' @export
setMethod("getTopLoci", "FineMappingEntry",
  function(x, type = c("data.frame", "GRanges"), ...) {
    type <- match.arg(type)
    tl <- x@topLoci
    if (type == "data.frame") return(tl)
    if (is.null(tl) || nrow(tl) == 0L) return(GRanges())
    parsed <- parseVariantId(tl$variant_id)
    gr <- GRanges(
      seqnames = paste0("chr", parsed$chrom),
      ranges = IRanges(start = parsed$pos, width = 1L)
    )
    mcols(gr) <- DataFrame(tl)
    gr
  })

#' @rdname getPip
#' @export
setMethod("getPip", "FineMappingEntry", function(x, ...) {
  tl <- x@topLoci
  if (nrow(tl) == 0L || !"pip" %in% names(tl)) return(numeric(0))
  setNames(tl$pip, tl$variant_id)
})

#' @rdname getCs
#' @export
setMethod("getCs", "FineMappingEntry",
  function(x, coverage = 0.95, ...) {
    tl <- x@topLoci
    if (nrow(tl) == 0L) return(data.frame())
    csCol <- grep(paste0("^cs.*", coverage * 100), names(tl), value = TRUE)
    if (length(csCol) == 0L && "cs" %in% names(tl)) csCol <- "cs"
    if (length(csCol) == 0L) return(data.frame())
    tl[tl[[csCol[1L]]] > 0, , drop = FALSE]
  })

#' @rdname adjustPips
#' @export
setMethod("adjustPips", "FineMappingEntry",
  function(x, keepVariants, ...) {
    keepVariants <- as.character(keepVariants)
    common <- intersect(x@variantIds, keepVariants)
    if (!length(common))
      stop("adjustPips: intersection of entry variants with `keepVariants` ",
           "is empty.")
    keepIdx <- match(common, x@variantIds)
    fit <- x@trimmedFit
    if (is.null(fit$lbf_variable))
      stop("adjustPips: entry's trimmedFit has no `lbf_variable` matrix; ",
           "PIP renormalization requires lbf_variable.")
    lbfSub <- fit$lbf_variable[, keepIdx, drop = FALSE]
    fit$lbf_variable <- lbfSub
    fit$alpha <- lbfToAlpha(lbfSub)
    fit$pip <- as.numeric(1 - apply(1 - fit$alpha, 2, prod))
    if (!is.null(fit$mu))
      fit$mu <- if (length(dim(fit$mu)) == 3)
                  fit$mu[, keepIdx, , drop = FALSE]
                else fit$mu[, keepIdx, drop = FALSE]
    if (!is.null(fit$mu2))
      fit$mu2 <- if (length(dim(fit$mu2)) == 3)
                   fit$mu2[, keepIdx, , drop = FALSE]
                 else fit$mu2[, keepIdx, drop = FALSE]
    if (!is.null(fit$X_column_scale_factors))
      fit$X_column_scale_factors <- fit$X_column_scale_factors[keepIdx]
    newTopLoci <- x@topLoci
    if (nrow(newTopLoci) > 0L)
      newTopLoci <- newTopLoci[newTopLoci$variant_id %in% common, ,
                                drop = FALSE]
    if (nrow(newTopLoci) > 0L && "pip" %in% names(newTopLoci)) {
      pipByVid <- setNames(fit$pip, common)
      newTopLoci$pip <- unname(pipByVid[newTopLoci$variant_id])
    }
    new("FineMappingEntry",
        variantIds = common,
        trimmedFit = fit,
        topLoci = newTopLoci,
        sumstats = x@sumstats)
  })

#' @export
setMethod("show", "FineMappingEntry", function(object) {
  nCs <- if (nrow(object@topLoci) > 0L && "cs" %in% names(object@topLoci))
           length(unique(object@topLoci$cs[object@topLoci$cs > 0]))
         else 0L
  cat(sprintf("FineMappingEntry: %d variants, %d credible sets\n",
              length(object@variantIds), nCs))
})
