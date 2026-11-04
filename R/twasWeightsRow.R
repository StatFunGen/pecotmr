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
#'   \code{\link{TwasWeightsRow}} and pass a list of them as
#'   \code{\link{TwasWeights}}' \code{entry} argument.
#' @slot variants A \code{GRanges} of the weighted variants.
#' @slot weights Per-variant weights, aligned to \code{variants}.
#' @slot methodFits Per-method fit payload, or \code{NULL}.
#' @slot cvResult Cross-validation payload, or \code{NULL}.
#' @slot weightStandardized Whether the weights are on the standardized scale.
#' @slot weightsDataType Optional molecular type.
#' @section Coercion: \code{as(x, "GRanges")} gives the row's variants with
#'   their \code{weight} column: the element a collection stores for it.
#' @seealso \code{\link{TwasWeightsRow}},
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
#' @param weightsDataType Optional molecular type (\code{"expression"},
#'   \code{"splicing"}).
#' @return A \code{\linkS4class{TwasWeightsRow}}.
#' @seealso \code{\link{FineMappingRow}}
#' @examples
#' row <- TwasWeightsRow(
#'     variantIds = c("chr1:100:A:G", "chr1:200:C:T"),
#'     weights = c(0.4, -0.2)
#' )
#' TwasWeights(
#'     studyName = "s1", context = "brain", trait = "g1", method = "lasso",
#'     entry = list(row)
#' )
#' @export
TwasWeightsRow <- function(
    variantIds,
    weights,
    methodFits = NULL,
    cvResult = NULL,
    weightStandardized = FALSE,
    weightsDataType = NULL
) {
    vids <- as.character(variantIds)
    gr <- .variantIdsToGRanges(vids, "variantIds")
    w <- .unnamedWeights(weights)
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
    new(
        "TwasWeightsRow",
        variants = withWeight,
        weights = w,
        methodFits = methodFits,
        cvResult = cvResult,
        weightStandardized = isTRUE(weightStandardized),
        weightsDataType = weightsDataType
    )
}

# A row saved while rows recorded the variants their weights were fitted on
# carries a `fitVariants` slot; the fits name their variants now, so the slot
# is dropped.
#' @rdname updateObject-methods
#' @export
setMethod(
    "updateObject",
    "TwasWeightsRow",
    function(object, ..., verbose = FALSE) {
        if (isTRUE(verbose)) {
            message("updateObject(object = 'TwasWeightsRow')")
        }
        if (!is.null(attr(object, "fitVariants", exact = TRUE))) {
            object <- `attr<-`(object, "fitVariants", NULL)
        }
        methods::callNextMethod(object, ..., verbose = verbose)
    }
)


# ---- field accessors --------------------------------------------------------

#' @rdname variantIds
#' @export
setMethod("variantIds", "TwasWeightsRow", function(x) {
    .grVariantIds(x@variants)
})

#' @rdname weights-methods
#' @export
setMethod("weights", "TwasWeightsRow", function(object, ...) {
    .namedWeights(object@weights, .grVariantIds(object@variants))
})

# Weights as a row stores them: without the variant ids as names, which the
# row's variants already hold -- every id would otherwise be stored a third
# time. A matrix keeps its column labels (conditions, or fSuSiE's features).
# @noRd
.unnamedWeights <- function(w) {
    if (is.null(w)) {
        return(w)
    }
    if (is.null(dim(w))) {
        return(unname(w))
    }
    `rownames<-`(w, NULL)
}

# Stored weights named by their variants again, as the accessors return them.
# @noRd
.namedWeights <- function(w, ids) {
    if (is.null(w)) {
        return(w)
    }
    if (is.null(dim(w))) {
        return(if (length(w) == length(ids)) set_names(w, ids) else w)
    }
    if (nrow(w) == length(ids)) `rownames<-`(w, ids) else w
}

#' @rdname methodFits
#' @export
setMethod("methodFits", "TwasWeightsRow", function(x) x@methodFits)

#' @rdname cvResult
#' @export
setMethod("cvResult", "TwasWeightsRow", function(x) x@cvResult)

# The variants a row's stored fit was made on, read off the fit itself; none
# for a row kept without a fit.
# @noRd
.twRowFitVariants <- function(x) .twasFitVariantIds(x@methodFits)

# ---- the variants a stored fit was made on ----------------------------------

# The variants a TWAS method's stored fit was made on. A fit records them once:
# where its engine already names them, or in the `variantIds` element
# .withFitVariantIds() gives a fit that names none. A row joined from several
# regions stores a list of the regions' fits, whose variants run on in turn.
# None for no fit, or one that names nothing.
# @noRd
.twasFitVariantIds <- function(fit) {
    if (!is.list(fit)) {
        return(character(0))
    }
    ids <- .twasFitOwnIds(fit)
    if (length(ids) > 0L || !all(map_lgl(fit, is.list))) {
        return(ids)
    }
    list_c(map(unname(fit), .twasFitOwnIds), ptype = character())
}

# One fit's variants: its `variantIds`, else where its engine names them --
# glmnet's and ncvreg's coefficient paths, qgg's and BGLR's effect vectors,
# mr.ash's and DPR's copies of the genotypes, a SuSiE-family fit's own names.
# @noRd
.twasFitOwnIds <- function(fit) {
    ids <- fit[["variantIds"]] %||%
        rownames(pluck(fit, "glmnet.fit", "beta")) %||%
        .twasNcvregFitIds(fit) %||%
        names(fit[["bm"]]) %||%
        names(pluck(fit, "ETA", 1L, "b")) %||%
        colnames(pluck(fit, "data", "X")) %||%
        colnames(fit[["x"]])
    if (is.null(ids)) {
        return(.fitNamedVariants(fit))
    }
    as.character(ids)
}

# ncvreg names its coefficient path's rows, the intercept's among them.
# @noRd
.twasNcvregFitIds <- function(fit) {
    ids <- rownames(pluck(fit, "fit", "beta"))
    ids[ids != "(Intercept)"]
}

# `fit` naming the variants it was made on: as it is when its engine already
# names them (or it is no list to name them in, or there are no ids to give
# it), else with `ids` as one `variantIds` element -- a single copy, however
# many per-variant values the fit holds.
# @noRd
.withFitVariantIds <- function(fit, ids) {
    if (!is.list(fit) || is.null(ids) || length(.twasFitOwnIds(fit)) > 0L) {
        return(fit)
    }
    `[[<-`(fit, "variantIds", value = as.character(ids))
}

# Weights `w` with the fit they carry, if any, naming its variants.
# @noRd
.withWeightsFitIds <- function(w, ids) {
    fit <- attr(w, "fit", exact = TRUE)
    if (is.null(fit)) {
        return(w)
    }
    `attr<-`(w, "fit", .withFitVariantIds(fit, ids))
}

#' @rdname weightStandardized
#' @export
setMethod("weightStandardized", "TwasWeightsRow", function(x) {
    isTRUE(x@weightStandardized)
})

#' @rdname dataType
#' @export
setMethod("dataType", "TwasWeightsRow", function(x) x@weightsDataType)

# @noRd
setMethod("variants", "TwasWeightsRow", function(x) x@variants)

#' @rdname TwasWeightsRow-class
#' @name coerce,TwasWeightsRow,GRanges-method
#' @aliases coerce,TwasWeightsRow,GRanges-method
#' @usage NULL
#' @exportMethod coerce
setAs("TwasWeightsRow", "GRanges", function(from) variants(from))

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
