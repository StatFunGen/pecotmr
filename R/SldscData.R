#' @title S-LDSC input data container
#' @description An in-memory bundle of the loaded S-LDSC inputs, produced from
#'   the reader functions (\code{\link{readSldscAnnot}}, \code{\link{readSldscFrq}},
#'   \code{\link{readSldscTrait}}) and consumed by
#'   \code{\link{sldscPostprocessingPipeline}}. The class itself performs no
#'   file I/O: the user runs the readers, then constructs an \code{SldscData}
#'   from those in-memory objects, and the pipeline does all computation on it.
#' @slot annot A \code{data.frame} of target annotations with at least
#'   \code{CHR} and \code{SNP} columns plus one or more annotation columns
#'   (\code{BP}/\code{CM} optional).
#' @slot frq A \code{data.frame} of reference-panel allele frequencies with
#'   \code{SNP} and \code{MAF} columns (a 0-row frame when no \code{.frq} data
#'   was supplied).
#' @slot traits A named list, one entry per trait, each a list with a
#'   \code{single} element (list of per-target \code{\link{readSldscTrait}}
#'   runs) and an optional \code{joint} element (a single run, or \code{NULL}).
#' @name SldscData-class
#' @include AllGenerics.R
#' @importFrom methods new validObject is
#' @exportClass SldscData
NULL

setClass("SldscData",
  slots = c(
    annot  = "data.frame",
    frq    = "data.frame",
    traits = "list"
  ),
  prototype = list(
    annot  = data.frame(),
    frq    = data.frame(),
    traits = list()
  ))

setValidity("SldscData", function(object) {
  errs <- character(0)

  annot <- object@annot
  if (!all(c("CHR", "SNP") %in% names(annot)))
    errs <- c(errs, "`annot` must have columns CHR and SNP.")
  annotCols <- setdiff(names(annot), c("CHR", "SNP", "BP", "CM"))
  if (length(annotCols) == 0L)
    errs <- c(errs, "`annot` must have at least one annotation column beyond CHR/SNP/BP/CM.")

  frq <- object@frq
  if (nrow(frq) > 0L && !all(c("SNP", "MAF") %in% names(frq)))
    errs <- c(errs, "non-empty `frq` must have columns SNP and MAF.")

  tr <- object@traits
  if (length(tr) > 0L) {
    if (is.null(names(tr)) || any(!nzchar(names(tr))))
      errs <- c(errs, "`traits` must be a named list (one entry per trait).")
    for (nm in names(tr)) {
      t <- tr[[nm]]
      if (!is.list(t) || !("single" %in% names(t)))
        errs <- c(errs, sprintf(
          "traits[['%s']] must be a list with a `single` element.", nm))
      else if (!is.list(t$single))
        errs <- c(errs, sprintf(
          "traits[['%s']]$single must be a list of runs.", nm))
    }
  }

  if (length(errs)) errs else TRUE
})

#' Construct an SldscData object
#'
#' Bundles the in-memory outputs of the S-LDSC readers into a single object for
#' \code{\link{sldscPostprocessingPipeline}}. Performs no file I/O.
#'
#' @param annot A target-annotation \code{data.frame} (e.g. from
#'   \code{\link{readSldscAnnot}}): \code{CHR}, \code{SNP}, and one or more
#'   annotation columns.
#' @param frq Optional reference-panel allele-frequency \code{data.frame} (e.g.
#'   from \code{\link{readSldscFrq}}): \code{SNP}, \code{MAF}. \code{NULL} (the
#'   default) stores an empty frame, which disables MAF-based filtering.
#' @param traits A named list of per-trait runs; each entry a list with a
#'   \code{single} list (per-target \code{\link{readSldscTrait}} outputs) and an
#'   optional \code{joint} run.
#' @return An \code{SldscData} object.
#' @seealso \code{\link{readSldscAnnot}}, \code{\link{readSldscFrq}},
#'   \code{\link{readSldscTrait}}, \code{\link{sldscPostprocessingPipeline}}
#' @rdname SldscData
#' @export
SldscData <- function(annot, frq = NULL, traits = list()) {
  if (missing(annot)) stop("SldscData: `annot` is required.")
  if (is.null(frq)) frq <- data.frame()
  obj <- new("SldscData",
             annot  = as.data.frame(annot),
             frq    = as.data.frame(frq),
             traits = traits)
  validObject(obj)
  obj
}

# ---- accessors ----

#' @rdname getAnnotData
#' @export
setMethod("getAnnotData", "SldscData", function(x) x@annot)

#' @rdname getFrqData
#' @export
setMethod("getFrqData", "SldscData", function(x) x@frq)

#' @rdname getTraitRuns
#' @export
setMethod("getTraitRuns", "SldscData", function(x) x@traits)

#' @rdname getTraitNames
#' @export
setMethod("getTraitNames", "SldscData", function(x) names(x@traits))

#' @rdname getAnnotCols
#' @export
setMethod("getAnnotCols", "SldscData",
  function(x) setdiff(names(x@annot), c("CHR", "SNP", "BP", "CM")))

#' @rdname getTraitRun
#' @export
setMethod("getTraitRun", "SldscData",
  function(x, trait, mode = c("single", "joint"), idx = NULL) {
    mode <- match.arg(mode)
    t <- x@traits[[trait]]
    if (is.null(t)) return(NULL)
    if (mode == "joint") return(t$joint)
    if (is.null(idx)) return(t$single)
    if (idx > length(t$single)) return(NULL)
    t$single[[idx]]
  })

#' @rdname SldscData
setMethod("show", "SldscData", function(object) {
  cat("SldscData\n")
  cat("  annotations (", length(getAnnotCols(object)), "): ",
      paste(getAnnotCols(object), collapse = ", "), "\n", sep = "")
  cat("  annot SNPs: ", nrow(object@annot),
      " | frq SNPs: ", nrow(object@frq), "\n", sep = "")
  cat("  traits (", length(object@traits), "): ",
      paste(names(object@traits), collapse = ", "), "\n", sep = "")
  invisible(object)
})
