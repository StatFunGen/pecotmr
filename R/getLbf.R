#' @include AllGenerics.R AllClasses.R FineMappingEntry.R tupleSelectors.R
NULL

#' Per-variant per-effect log Bayes factors (wide)
#'
#' The variant x effect matrix of single-effect log Bayes factors
#' (\code{lbf_variable}, or fSuSiE's \code{lBF}) from the stored fit: one row
#' per variant, one \code{lbf_L<k>} column per effect. The per-variant scalar
#' summary (max across effects) is the \code{logBF} column of
#' \code{\link{getTopLoci}}; this accessor keeps the full per-effect breakdown.
#' Across entries with different effect counts the collection method NA-fills
#' the ragged columns.
#'
#' @param x A \code{FineMappingEntry} or \code{FineMappingResultBase}.
#' @param ... Ignored.
#' @return A \code{data.frame}: \code{variant_id} + \code{lbf_L1..lbf_LL} (the
#'   collection method also carries the entry identity columns).
#' @seealso \code{\link{getTopLoci}} (the scalar \code{logBF} column)
#' @examples
#' data(qtlFineMappingExample)
#' getLbf(qtlFineMappingExample)
#' @export
setGeneric("getLbf", function(x, ...) standardGeneric("getLbf"))

#' @rdname getLbf
#' @export
setMethod("getLbf", "FineMappingEntry", function(x, ...) {
    lbf <- .asLbfMatrix(getSusieFit(x))
    vids <- x@variantIds
    if (is.null(lbf) || ncol(lbf) != length(vids)) {
        return(data.frame(variant_id = character(0), stringsAsFactors = FALSE))
    }
    w <- as.data.frame(t(as.matrix(lbf)), stringsAsFactors = FALSE)
    names(w) <- paste0("lbf_L", seq_len(ncol(w)))
    cbind(
        data.frame(variant_id = as.character(vids), stringsAsFactors = FALSE),
        w
    )
})

#' @rdname getLbf
#' @export
setMethod("getLbf", "FineMappingResultBase", function(x, ...) {
    .fmrAggregateView(x, perEntry = function(e) getLbf(e))
})
