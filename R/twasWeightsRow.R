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
#' @slot fits Per-method fit payload, or \code{NULL}.
#' @slot cvResult Cross-validation payload, or \code{NULL}.
#' @slot standardized Whether the weights are on the standardized scale.
#' @slot dataType Optional data-type label.
#' @seealso \code{\link{twasWeightsRow}},
#'   \code{\linkS4class{FineMappingRow}}
#' @export
setClass(
    "TwasWeightsRow",
    representation(
        variants = "GRanges",
        weights = "ANY",
        fits = "ANY",
        cvResult = "ANY",
        standardized = "logical",
        dataType = "ANY"
    ),
    prototype(
        weights = NULL,
        fits = NULL,
        cvResult = NULL,
        standardized = FALSE,
        dataType = NULL
    )
)

methods::setValidity("TwasWeightsRow", function(object) {
    errors <- character(0)
    w <- object@weights
    n <- length(object@variants)
    if (!is.null(w)) {
        # A matrix carries one ROW per variant (columns are conditions), so
        # the two shapes need separate checks -- validating only the vector
        # case would let a mis-sized matrix through.
        if (is.null(dim(w)) && length(w) != n) {
            errors <- c(
                errors,
                "length(weights) must equal length(variantIds)"
            )
        } else if (!is.null(dim(w)) && nrow(w) != n) {
            errors <- c(
                errors,
                "nrow(weights) must equal length(variantIds)"
            )
        }
    }
    if (length(object@standardized) != 1L || is.na(object@standardized)) {
        errors <- c(errors, "'standardized' must be a single logical value")
    }
    if (length(errors) == 0L) TRUE else errors
})

#' @title Build One TWAS-Weight Row
#' @description Assemble a single row's payload for
#'   \code{\link{TwasWeights}}: the variants, their weights and the fit
#'   payload. Pass a list of these as the collection's \code{entry} argument.
#' @param variantIds Character vector of variant ids, each encoding
#'   coordinates (\code{chrom:pos:ref:alt}).
#' @param weights Per-variant weights, aligned to \code{variantIds} (a vector,
#'   or a matrix with one column per condition).
#' @param fits Optional per-method fit payload.
#' @param cvResult Optional cross-validation payload.
#' @param standardized Whether the weights are already on the standardized
#'   scale.
#' @param dataType Optional data-type label.
#' @return A \code{\linkS4class{TwasWeightsRow}}.
#' @seealso \code{\link{fineMappingRow}}
#' @examples
#' row <- twasWeightsRow(
#'     variantIds = c("chr1:100:A:G", "chr1:200:C:T"),
#'     weights = c(0.4, -0.2)
#' )
#' TwasWeights(
#'     study = "s1", context = "brain", trait = "g1", method = "lasso",
#'     entry = list(row)
#' )
#' @export
twasWeightsRow <- function(
    variantIds,
    weights,
    fits = NULL,
    cvResult = NULL,
    standardized = FALSE,
    dataType = NULL
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
    mcols(gr)$weight <- w
    obj <- new(
        "TwasWeightsRow",
        variants = gr,
        weights = w,
        fits = fits,
        cvResult = cvResult,
        standardized = isTRUE(standardized),
        dataType = dataType
    )
    validObject(obj)
    obj
}


# Variant ids must encode coordinates: the id is a rendering of the variant's
# range and alleles, so an id that carries neither (e.g. "v1") names nothing.
# Enforced here so the requirement bites before the entry becomes a ranged
# element, rather than surfacing later as an unparseable-coordinate failure.
# @noRd
.tweCheckVariantIdentity <- function(object) {
    ids <- object@variantIds
    if (length(ids) == 0L) {
        return(NULL)
    }
    parsed <- parseVariantId(ids)
    bad <- is.na(parsed$chrom) | is.na(parsed$pos)
    if (!any(bad)) {
        return(NULL)
    }
    glue(
        "variantIds: {sum(bad)} of {length(ids)} do not encode coordinates ",
        "(expected chrom:pos:ref:alt), e.g. ",
        "{str_flatten(ids[bad][seq_len(min(3L, sum(bad)))], ', ')}"
    )
}

# ---- field accessors --------------------------------------------------------

#' @rdname getVariantIds
#' @export
setMethod("getVariantIds", "TwasWeightsRow", function(x, ...) {
    .grVariantIds(x@variants)
})

#' @rdname getWeights
#' @export
setMethod("getWeights", "TwasWeightsRow", function(x, ...) x@weights)

#' @rdname getFits
#' @export
setMethod("getFits", "TwasWeightsRow", function(x, ...) x@fits)

#' @rdname getCvResult
#' @export
setMethod("getCvResult", "TwasWeightsRow", function(x, ...) x@cvResult)

#' @rdname getStandardized
#' @export
setMethod("getStandardized", "TwasWeightsRow", function(x, ...) {
    isTRUE(x@standardized)
})

#' @rdname getDataType
#' @export
setMethod("getDataType", "TwasWeightsRow", function(x, ...) x@dataType)

# @noRd
setMethod("rowVariants", "TwasWeightsRow", function(x, ...) x@variants)

#' @rdname show-methods
#' @export
setMethod("show", "TwasWeightsRow", function(object) {
    cat(glue(
        "TwasWeightsRow: {length(object@variants)} variants, ",
        "standardized={object@standardized}\n",
        .trim = FALSE
    ))
    hasCv <- !is.null(object@cvResult)
    cat(glue("  CV performance: {hasCv}\n", .trim = FALSE))
    invisible(NULL)
})
