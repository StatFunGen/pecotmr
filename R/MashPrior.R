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
#              samplePartition  data.frame(Sample, Fold) — the CV folds the
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
#' @slot fullFit The full-data data-driven prior payload — the
#'   \code{mashPipeline()} output \code{list(U, w)}; fed to mr.mash as
#'   \code{dataDrivenPriorMatrices} for the full-data fit. \code{NULL} when only
#'   per-fold priors are supplied (a CV-only run).
#' @slot cvFits \code{NULL}, or a list with \code{samplePartition}
#'   (\code{data.frame(Sample, Fold)}) and \code{perFoldFits} (a list of per-fold
#'   prior payloads, \code{perFoldFits[[j]]} for fold \code{j}).
#' @export
setClass("MashPrior",
  representation(
    fullFit = "ANY",
    cvFits  = "ANY"),
  validity = function(object) {
    errors <- character()
    if (is.null(object@fullFit) && is.null(object@cvFits)) {
      errors <- c(errors,
        "a MashPrior must carry at least one of `fullFit` or `cvFits`")
    }
    cv <- object@cvFits
    if (!is.null(cv)) {
      if (!is.list(cv) || is.null(cv$perFoldFits)) {
        errors <- c(errors,
          "`cvFits` must be a list with a `perFoldFits` element")
      } else {
        if (!is.list(cv$perFoldFits) || length(cv$perFoldFits) == 0L) {
          errors <- c(errors, "`cvFits$perFoldFits` must be a non-empty list")
        }
        sp <- cv$samplePartition
        if (!is.null(sp)) {
          if (!is.data.frame(sp) || !all(c("Sample", "Fold") %in% names(sp))) {
            errors <- c(errors,
              "`cvFits$samplePartition` must be a data.frame with `Sample` and `Fold` columns")
          } else if (is.list(cv$perFoldFits)) {
            nF <- length(unique(sp$Fold))
            if (length(cv$perFoldFits) != nF) {
              errors <- c(errors, sprintf(
                "`cvFits$perFoldFits` has %d element(s) but the partition defines %d fold(s)",
                length(cv$perFoldFits), nF))
            }
          }
        }
      }
    }
    if (length(errors) == 0L) TRUE else errors
  }
)

#' @title Create a MashPrior Object
#' @description Construct a \code{\link{MashPrior}} bundling a full-data
#'   data-driven prior with optional per-fold (cross-validated) priors.
#' @param fullFit Full-data data-driven prior payload (the
#'   \code{mashPipeline()} \code{list(U, w)} output), or \code{NULL} for a
#'   CV-only bundle.
#' @param cvFits \code{NULL}, or a list with \code{perFoldFits} (a non-empty
#'   list of per-fold prior payloads) and optionally \code{samplePartition}
#'   (\code{data.frame(Sample, Fold)}).
#' @return A \code{MashPrior} object.
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

#' @export
setMethod("show", "MashPrior", function(object) {
  cat("MashPrior\n")
  cat(sprintf("  fullFit: %s\n",
              if (is.null(object@fullFit)) "none" else "present"))
  cv <- object@cvFits
  if (is.null(cv)) {
    cat("  cvFits: none\n")
  } else {
    nF <- if (!is.null(cv$perFoldFits)) length(cv$perFoldFits) else 0L
    cat(sprintf("  cvFits: %d per-fold prior(s)%s\n", nF,
                if (!is.null(cv$samplePartition)) " + samplePartition" else ""))
  }
})
