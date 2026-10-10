# =============================================================================
# TwasWeightsRow S4 class
# -----------------------------------------------------------------------------
# One row's TWAS payload: the variants (with their per-variant weight as a
# metadata column) plus the fit payload. Sibling of FineMappingRow; see that
# class for why neither carries view methods.
# =============================================================================

#' @include AllGenerics.R
NULL

#' @title TWAS Weights Row
#' @description One row's worth of TWAS weights. Build one with
#'   \code{\link{twasWeightsRow}} and pass a list of them as
#'   \code{\link{TwasWeights}}' \code{entry} argument.
#' @slot variants A \code{GRanges} of the weighted variants.
#' @slot weights Per-variant weights, aligned to \code{variants}.
#' @slot methodFits Per-method fit payload, or \code{NULL}.
#' @slot cvResult Cross-validation payload, or \code{NULL}.
#' @slot weightStandardized Whether the weights are on the standardized scale.
#' @slot weightsDataType Optional data-type label.
#' @seealso \code{\link{twasWeightsRow}},
#'   \code{\linkS4class{FineMappingRow}}
#' @export
setClass(
    "TwasWeightsRow",
    representation(
        variants = "GRanges",
        weights = "ANY",
        methodFits = "ANY",
        cvResult = "ANY",
        weightStandardized = "logical",
        weightsDataType = "ANY"
    ),
    prototype(
        weights = NULL,
        methodFits = NULL,
        cvResult = NULL,
        weightStandardized = FALSE,
        weightsDataType = NULL
    )
)

#' @importFrom checkmate makeAssertCollection assert assertFlag
#' @importFrom checkmate checkNumeric checkMatrix
methods::setValidity("TwasWeightsRow", function(object) {
    coll <- makeAssertCollection()
    w <- object@weights
    n <- length(object@variants)
    if (!is.null(w)) {
        # A matrix carries one ROW per variant (columns are conditions), so
        # the two shapes need separate checks -- validating only the vector
        # case would let a mis-sized matrix through.
        assert(
            checkNumeric(w, len = n),
            checkMatrix(w, nrows = n),
            .var.name = "weights",
            add = coll
        )
    }
    assertFlag(
        object@weightStandardized,
        .var.name = "weightStandardized",
        add = coll
    )
    coll$getMessages()
})

#' @title Build One TWAS-Weight Row
#' @description Assemble a single row's payload for
#'   \code{\link{TwasWeights}}: the variants, their weights and the fit
#'   payload. Pass a list of these as the collection's \code{entry} argument.
#' @param variantIds Character vector of variant ids, each encoding
#'   coordinates (\code{chrom:pos:ref:alt}).
#' @param weights Per-variant weights, aligned to \code{variantIds} (a vector,
#'   or a matrix with one column per condition).
#' @param methodFits Optional per-method fit payload.
#' @param cvResult Optional cross-validation payload.
#' @param weightStandardized Whether the weights are already on the standardized
#'   scale.
#' @param weightsDataType Optional data-type label.
#' @return A \code{\linkS4class{TwasWeightsRow}}.
#' @seealso \code{\link{fineMappingRow}}
#' @examples
#' row <- twasWeightsRow(
#'     variantIds = c("chr1:100:A:G", "chr1:200:C:T"),
#'     weights = c(0.4, -0.2)
#' )
#' TwasWeights(
#'     studyName = "s1", context = "brain", trait = "g1", method = "lasso",
#'     entry = list(row)
#' )
#' @export
twasWeightsRow <- function(
    variantIds,
    weights,
    methodFits = NULL,
    cvResult = NULL,
    weightStandardized = FALSE,
    weightsDataType = NULL
) {
    vids <- as.character(variantIds)
    gr <- .variantIdsToGRanges(vids, "variantIds")
    w <- weights
    if (!is.null(w) && is.null(dim(w)) && length(w) != length(gr)) {
        msg <- glue(
            "length(weights) is {length(w)} but {length(gr)} variants were ",
            "supplied."
        )
        abort(msg)
    }
    if (!is.null(w) && !is.null(dim(w)) && nrow(w) != length(gr)) {
        msg <- glue(
            "nrow(weights) must equal length(variantIds) ",
            "(got {nrow(w)} vs {length(gr)})."
        )
        abort(msg)
    }
    withWeight <- S4Vectors::`mcols<-`(
        gr,
        value = `[[<-`(mcols(gr, use.names = FALSE), "weight", value = w)
    )
    obj <- new(
        "TwasWeightsRow",
        variants = withWeight,
        weights = w,
        methodFits = methodFits,
        cvResult = cvResult,
        weightStandardized = isTRUE(weightStandardized),
        weightsDataType = weightsDataType
    )
    validObject(obj)
    obj
}


# ---- field accessors --------------------------------------------------------

#' @rdname variantIds
#' @export
setMethod("variantIds", "TwasWeightsRow", function(x) {
    .grVariantIds(x@variants)
})

#' @rdname weights-methods
#' @export
setMethod("weights", "TwasWeightsRow", function(object, ...) object@weights)

#' @rdname methodFits
#' @export
setMethod("methodFits", "TwasWeightsRow", function(x) x@methodFits)

#' @rdname cvResult
#' @export
setMethod("cvResult", "TwasWeightsRow", function(x) x@cvResult)

#' @rdname weightStandardized
#' @export
setMethod("weightStandardized", "TwasWeightsRow", function(x) {
    isTRUE(x@weightStandardized)
})

#' @rdname weightsDataType
#' @export
setMethod("weightsDataType", "TwasWeightsRow", function(x) x@weightsDataType)

# @noRd
setMethod("variants", "TwasWeightsRow", function(x) x@variants)

#' @rdname show-methods
#' @export
setMethod("show", "TwasWeightsRow", function(object) {
    cat(glue(
        "TwasWeightsRow: {length(object@variants)} variants, ",
        "standardized={object@weightStandardized}\n",
        .trim = FALSE
    ))
    hasCv <- !is.null(object@cvResult)
    cat(glue("  CV performance: {hasCv}\n", .trim = FALSE))
    invisible(NULL)
})
