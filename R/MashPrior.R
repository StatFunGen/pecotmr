# =============================================================================
# MashPrior S4 class
# -----------------------------------------------------------------------------
# Data-driven (mash) prior bundle consumed by twasWeightsPipeline (mr.mash) and,
# transitively, by fineMappingPipeline (mvSuSiE). It is an INPUT-only container:
# mashPipeline() produces the prior payload(s); this class packages a full-data
# prior together with optional per-fold (cross-validated) priors plus the fold
# partition they were computed on, so honest CV reuses the SAME folds.
#
#   fullFit  the full-data data-driven prior payload: the mashPipeline() output
#            list(U = <covariance list>, w = <mixture weights>), handed to
#            mr.mash as `dataDrivenPriorMatrices` for the full-data fit.
#   cvFits   NULL, or a list with:
#              samplePartition  data.frame(Sample, Fold) -- the CV folds the
#                               per-fold priors were computed on.
#              perFoldFits      list of per-fold prior payloads (each the same
#                               shape as `fullFit`); perFoldFits[[j]] is the
#                               prior for fold j, ordered to match
#                               sort(unique(Fold)).
# =============================================================================

#' @include AllGenerics.R
NULL

#' @title Data-Driven (mash) Prior Bundle
#' @description Input container packaging a full-data data-driven prior with
#'   optional per-fold (cross-validated) priors and the fold partition they were
#'   computed on. Produced (eventually) by \code{mashPipeline()} and consumed by
#'   \code{twasWeightsPipeline()} (mr.mash); the per-fold fits then flow to
#'   \code{fineMappingPipeline()} (mvSuSiE) via the resulting
#'   \code{\link{TwasWeights}}.
#' @slot fullFit The full-data data-driven prior payload -- the
#'   \code{mashPipeline()} output \code{list(U, w)}; fed to mr.mash as
#'   \code{dataDrivenPriorMatrices} for the full-data fit. \code{NULL} when only
#'   per-fold priors are supplied (a CV-only run).
#' @slot cvFits \code{NULL}, or a list with \code{samplePartition}
#'   (\code{data.frame(Sample, Fold)}) and \code{perFoldFits} (a list of
#'   per-fold prior payloads, \code{perFoldFits[[j]]} for fold \code{j}).
#' @export
setClass(
    "MashPrior",
    representation(
        fullFit = "ANY",
        cvFits = "ANY"
    ),
    validity = function(object) .validateMashPrior(object)
)

# ---- MashPrior validity helpers --------------------------------------------

# @noRd
.validateMashPrior <- function(object) {
    errors <- c(
        .mashPriorCheckPresence(object),
        .mashPriorCheckCvFits(object@cvFits)
    )
    if (length(errors) == 0L) TRUE else errors
}

# @noRd
.mashPriorCheckPresence <- function(object) {
    if (is.null(object@fullFit) && is.null(object@cvFits)) {
        return("a MashPrior must carry at least one of `fullFit` or `cvFits`")
    }
    NULL
}

# @noRd
.mashPriorCheckCvFits <- function(cv) {
    if (is.null(cv)) {
        return(NULL)
    }
    if (!is.list(cv) || is.null(cv$perFoldFits)) {
        return("`cvFits` must be a list with a `perFoldFits` element")
    }
    c(
        .mashPriorCheckPerFoldFits(cv$perFoldFits),
        .mashPriorCheckPartition(cv)
    )
}

# @noRd
.mashPriorCheckPerFoldFits <- function(perFoldFits) {
    if (!is.list(perFoldFits) || length(perFoldFits) == 0L) {
        return("`cvFits$perFoldFits` must be a non-empty list")
    }
    NULL
}

# Fold-partition consistency (only meaningful when perFoldFits is a list).
# @noRd
.mashPriorCheckPartition <- function(cv) {
    sp <- cv$samplePartition
    if (is.null(sp)) {
        return(NULL)
    }
    if (!is.data.frame(sp) || !all(is_in(c("Sample", "Fold"), names(sp)))) {
        return(str_c(
            "`cvFits$samplePartition` must be a data.frame ",
            "with `Sample` and `Fold` columns"
        ))
    }
    if (!is.list(cv$perFoldFits)) {
        return(NULL)
    }
    nF <- n_distinct(sp$Fold)
    if (length(cv$perFoldFits) != nF) {
        return(glue(
            "`cvFits$perFoldFits` has {length(cv$perFoldFits)} element(s) ",
            "but the partition defines {nF} fold(s)"
        ))
    }
    NULL
}

#' @title Create a MashPrior Object
#' @description Construct a \code{\link{MashPrior}} bundling a full-data
#'   data-driven prior with optional per-fold (cross-validated) priors.
#' @param fullFit Full-data data-driven prior payload (the \code{mashPipeline()}
#'   \code{list(U, w)} output), or \code{NULL} for a CV-only bundle.
#' @param cvFits \code{NULL}, or a list with \code{perFoldFits} (a non-empty
#'   list of per-fold prior payloads) and optionally \code{samplePartition}
#'   (\code{data.frame(Sample, Fold)}).
#' @return A \code{MashPrior} object.
#' @examples
#' U <- list(shared = diag(3), singleton = matrix(0.3, 3, 3) + diag(0.7, 3))
#' mp <- MashPrior(fullFit = list(U = U, w = c(0.5, 0.5)))
#' mp
#' @export
MashPrior <- function(fullFit = NULL, cvFits = NULL) {
    obj <- new("MashPrior", fullFit = fullFit, cvFits = cvFits)
    validObject(obj)
    obj
}

#' @rdname getFullFit
#' @export
setMethod("getFullFit", "MashPrior", function(x, ...) x@fullFit)

#' @rdname getCvFits
#' @export
setMethod("getCvFits", "MashPrior", function(x, ...) x@cvFits)

#' @rdname show-methods
#' @export
setMethod("show", "MashPrior", function(object) {
    cat("MashPrior\n")
    cat(glue(
        "  fullFit: {if (is.null(object@fullFit)) 'none' else 'present'}\n",
        .trim = FALSE
    ))
    cv <- object@cvFits
    if (is.null(cv)) {
        cat("  cvFits: none\n")
    } else {
        nF <- if (!is.null(cv$perFoldFits)) length(cv$perFoldFits) else 0L
        cat(glue(
            "  cvFits: {nF} per-fold prior(s)",
            "{if (!is.null(cv$samplePartition)) ' + samplePartition' else ''}",
            "\n",
            .trim = FALSE
        ))
    }
})
