# =============================================================================
# MashPrior S4 class
# -----------------------------------------------------------------------------
# Data-driven (mash) prior bundle consumed by twasWeightsPipeline (mr.mash) and,
# transitively, by fineMappingPipeline (mvSuSiE). It is an INPUT-only container:
# mashPipeline() produces the priors; this class packages the full-data priors
# together with optional per-fold (cross-validated) priors plus the fold
# partition they were computed on, so honest CV reuses the SAME folds.
#
#   fullFit  the full-data priors: the mashPipeline() output, a named list with
#            one prior per joint spec, each list(U, w, residualCorrelation,
#            alpha, inputScale, axes, scope, conditions). Each mr.mash group
#            takes the one over its spec's axes whose conditions include the
#            group's (.mashPriorForGroup()), as `dataDrivenPriorMatrices`.
#   cvFits   NULL, or a list with:
#              samplePartition  data.frame(Sample, Fold) -- the CV folds the
#                               per-fold priors were computed on.
#              perFoldFits      list of per-fold priors (each the same shape
#                               as `fullFit`); perFoldFits[[j]] is the prior
#                               for fold j, ordered to match
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
#' @slot fullFit The full-data priors -- the \code{\link{mashPipeline}}
#'   output: a named list with one prior per joint spec, each a list of
#'   \code{U}, \code{w}, \code{residualCorrelation}, \code{alpha},
#'   \code{inputScale}, \code{axes}, \code{scope} and \code{conditions}.
#'   Each mr.mash group takes the prior over its joint spec's axes whose
#'   conditions include the group's (and whose scope does); its \code{U} is
#'   taken on the group's conditions, in their order, and put on the per-SD
#'   scale mr.mash fits on by dividing entry \eqn{(k, l)} by
#'   \eqn{(n_k n_l)^{e/2}}, where \eqn{n} is the conditions' sample size and
#'   \eqn{e} is 1 on the z scale and \code{alpha} otherwise. Its
#'   \code{residualCorrelation}, a correlation, is taken on the group's
#'   conditions unscaled. A prior without one (from an earlier
#'   \code{mashPipeline()}) is accepted, and the mr.mash fits then take
#'   their residual covariance from the data. No such prior, or several, is
#'   an error. \code{NULL} when only per-fold priors are supplied (a CV-only
#'   run).
#' @slot cvFits \code{NULL}, or a list with \code{samplePartition}
#'   (\code{data.frame(Sample, Fold)}) and \code{perFoldFits} (a list of
#'   per-fold priors, each the shape of \code{fullFit},
#'   \code{perFoldFits[[j]]} for fold \code{j}).
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
        .mashPriorCheckPriors(object@fullFit, "fullFit"),
        .mashPriorCheckCvFits(object@cvFits)
    )
    if (length(errors) == 0L) TRUE else errors
}

# A set of priors as mashPipeline() returns them: a named list with one prior
# per joint spec. NULL is no set.
# @noRd
.mashPriorCheckPriors <- function(priors, what) {
    if (is.null(priors)) {
        return(NULL)
    }
    if (!is.list(priors) || length(priors) == 0L || !is_named(priors)) {
        return(glue(
            "`{what}` must be a named list with one prior per joint spec, ",
            "as mashPipeline() returns"
        ))
    }
    list_c(imap(priors, .mashPriorProblem, what = what))
}

# What is wrong with one prior, or NULL: it carries its covariances, the
# conditions their rows are, its axes and the scale it was fitted on.
# @noRd
.mashPriorProblem <- function(prior, name, what) {
    label <- glue("`{what}${name}`")
    fields <- c("U", "alpha", "inputScale", "axes", "conditions")
    if (!is.list(prior) || !all(is_in(fields, names(prior)))) {
        return(glue("{label} must carry {str_flatten(fields, ', ')}"))
    }
    conditionCols <- c("study", "context", "trait", "n")
    if (
        !is.data.frame(prior$conditions) ||
            !all(is_in(conditionCols, names(prior$conditions)))
    ) {
        return(glue(
            "{label}$conditions must be a data.frame with columns ",
            "{str_flatten(conditionCols, ', ')}"
        ))
    }
    k <- nrow(prior$conditions)
    if (
        !is.list(prior$U) ||
            length(prior$U) == 0L ||
            !all(map_lgl(prior$U, .mashCovIsSquareN, nCond = k))
    ) {
        return(glue(
            "{label}$U must be a non-empty list of {k} x {k} matrices, one ",
            "row per condition"
        ))
    }
    residualOk <- is.null(prior$residualCorrelation) ||
        .mashCovIsSquareN(prior$residualCorrelation, nCond = k)
    if (!residualOk) {
        return(glue(
            "{label}$residualCorrelation must be a {k} x {k} matrix, one row ",
            "per condition"
        ))
    }
    .mashPriorScaleProblem(prior, label)
}

# A prior's axes and scale: the joint-spec axes, a numeric alpha and the
# scale mashPipeline() resolved.
# @noRd
.mashPriorScaleProblem <- function(prior, label) {
    axesOk <- is.character(prior$axes) &&
        length(prior$axes) > 0L &&
        all(is_in(prior$axes, c("study", "context", "trait")))
    if (!axesOk) {
        return(glue("{label}$axes must name study, context or trait"))
    }
    scaleOk <- is.numeric(prior$alpha) &&
        length(prior$alpha) == 1L &&
        is_in(prior$inputScale, c("beta", "z", "standardized"))
    if (!scaleOk) {
        return(glue(
            "{label} must have a numeric `alpha` and an `inputScale` of ",
            "'beta', 'z' or 'standardized'"
        ))
    }
    NULL
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
    list_c(map(
        seq_along(perFoldFits),
        .mashPriorCheckFold,
        perFoldFits = perFoldFits
    ))
}

# @noRd
.mashPriorCheckFold <- function(j, perFoldFits) {
    .mashPriorCheckPriors(perFoldFits[[j]], glue("cvFits$perFoldFits[[{j}]]"))
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
#' @description Construct a \code{\link{MashPrior}} bundling full-data
#'   data-driven priors with optional per-fold (cross-validated) priors.
#' @param fullFit Full-data priors (the \code{\link{mashPipeline}} output: a
#'   named list with one prior per joint spec), or \code{NULL} for a CV-only
#'   bundle.
#' @param cvFits \code{NULL}, or a list with \code{perFoldFits} (a non-empty
#'   list of per-fold priors, each the shape of \code{fullFit}) and
#'   optionally \code{samplePartition} (\code{data.frame(Sample, Fold)}).
#' @return A \code{MashPrior} object.
#' @examples
#' U <- list(shared = matrix(1, 2, 2), brain = diag(c(1, 0)))
#' conditions <- data.frame(
#'   study = NA, context = c("brain", "liver"), trait = NA, n = 100
#' )
#' prior <- list(
#'   U = U, w = c(0.5, 0.5), alpha = 0, inputScale = "standardized",
#'   axes = "context", scope = NULL, conditions = conditions
#' )
#' mp <- MashPrior(fullFit = list(context = prior))
#' mp
#' @importFrom checkmate assertList
#' @export
MashPrior <- function(fullFit = NULL, cvFits = NULL) {
    assertList(cvFits, null.ok = TRUE)
    new("MashPrior", fullFit = fullFit, cvFits = cvFits)
}

#' @rdname fullFit
#' @export
setMethod("fullFit", "MashPrior", function(x) x@fullFit)

#' @rdname cvFits
#' @export
setMethod("cvFits", "MashPrior", function(x) x@cvFits)

#' @rdname show-methods
#' @export
setMethod("show", "MashPrior", function(object) {
    cat("MashPrior\n")
    full <- object@fullFit
    cat(glue(
        "  fullFit: ",
        "{if (is.null(full)) 'none' else str_flatten(names(full), ', ')}\n",
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
