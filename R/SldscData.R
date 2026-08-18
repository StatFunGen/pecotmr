#' @title S-LDSC input data container
#' @description An in-memory bundle of the loaded S-LDSC inputs, produced from
#'   the reader functions (\code{\link{readSldscAnnot}},
#'   \code{\link{readSldscFrq}}, \code{\link{readSldscTrait}}) and consumed by
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
#' @include AllGenerics.R
#' @importFrom methods new validObject is
#' @exportClass SldscData
setClass(
    "SldscData",
    slots = c(
        annot = "data.frame",
        frq = "data.frame",
        traits = "list"
    ),
    prototype = list(
        annot = tibble(),
        frq = tibble(),
        traits = list()
    )
)

setValidity("SldscData", function(object) .validateSldscData(object))

# ---- SldscData validity helpers --------------------------------------------

# @noRd
.validateSldscData <- function(object) {
    errs <- c(
        .sldscDataCheckAnnot(object@annot),
        .sldscDataCheckFrq(object@frq),
        .sldscDataCheckTraits(object@traits)
    )
    if (length(errs)) errs else TRUE
}

# @noRd
.sldscDataCheckAnnot <- function(annot) {
    errs <- character(0)
    if (!all(is_in(c("CHR", "SNP"), names(annot)))) {
        errs <- c(errs, "`annot` must have columns CHR and SNP.")
    }
    annotCols <- setdiff(names(annot), c("CHR", "SNP", "BP", "CM"))
    if (length(annotCols) == 0L) {
        errs <- c(
            errs,
            str_c(
                "`annot` must have at least one annotation ",
                "column beyond CHR/SNP/BP/CM."
            )
        )
    }
    errs
}

# @noRd
.sldscDataCheckFrq <- function(frq) {
    if (nrow(frq) > 0L && !all(is_in(c("SNP", "MAF"), names(frq)))) {
        return("non-empty `frq` must have columns SNP and MAF.")
    }
    NULL
}

# @noRd
.sldscDataCheckTraits <- function(tr) {
    if (length(tr) == 0L) {
        return(NULL)
    }
    errs <- character(0)
    if (is.null(names(tr)) || any(str_length(names(tr)) == 0L, na.rm = TRUE)) {
        errs <- c(errs, "`traits` must be a named list (one entry per trait).")
    }
    c(errs, unlist(compact(map(names(tr), .sldscDataCheckOneTrait, tr = tr))))
}

# @noRd
.sldscDataCheckOneTrait <- function(nm, tr) {
    t <- tr[[nm]]
    if (!is.list(t) || !is_in("single", names(t))) {
        return(glue(
            "traits[['{nm}']] must be a list with a `single` element."
        ))
    }
    if (!is.list(t$single)) {
        return(glue("traits[['{nm}']]$single must be a list of runs."))
    }
    NULL
}

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
#' @param object An \code{SldscData} object (used by the \code{show} method).
#' @return An \code{SldscData} object.
#' @seealso \code{\link{readSldscAnnot}}, \code{\link{readSldscFrq}},
#'   \code{\link{readSldscTrait}}, \code{\link{sldscPostprocessingPipeline}}
#' @rdname SldscData
#' @examples
#' mkRun <- function(cats) {
#'   n <- length(cats)
#'   list(categories = cats, tau = setNames(rep(1e-7, n), cats),
#'     tauSe = setNames(rep(3e-8, n), cats),
#'     enrichment = setNames(rep(2, n), cats),
#'     enrichmentSe = setNames(rep(0.4, n), cats),
#'     enrichmentP = setNames(rep(0.01, n), cats),
#'     propH2 = setNames(rep(0.2, n), cats),
#'     propSnps = setNames(rep(0.1, n), cats), h2g = 0.3,
#'     tauBlocks = matrix(1e-7, 10, n, dimnames = list(NULL, cats)),
#'     nBlocks = 10L)
#' }
#' annot <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   annot_A = c(1, 0, 1, 0, 1, 0), annot_B = c(2.1, 1.8, 2.5, 1.9, 2.3, 2))
#' frq <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   MAF = rep(0.2, 6))
#' mkTrait <- function() {
#'   list(single = list(mkRun(c("annot_A_0", "baselineLD_0")),
#'     mkRun(c("annot_B_0", "baselineLD_0"))),
#'     joint = mkRun(c("annot_A_0", "annot_B_0", "baselineLD_0")))
#' }
#' traits <- setNames(list(mkTrait(), mkTrait()), c("traitX", "traitY"))
#' sd <- SldscData(annot = annot, frq = frq, traits = traits)
#' sd
#' @export
SldscData <- function(annot, frq = NULL, traits = list()) {
    if (missing(annot)) {
        abort("SldscData: `annot` is required.")
    }
    if (is.null(frq)) {
        frq <- tibble()
    }
    obj <- new(
        "SldscData",
        annot = as_tibble(annot),
        frq = as_tibble(frq),
        traits = traits
    )
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
setMethod("getAnnotCols", "SldscData", function(x) {
    setdiff(names(x@annot), c("CHR", "SNP", "BP", "CM"))
})

#' @rdname getTraitRun
#' @export
setMethod(
    "getTraitRun",
    "SldscData",
    function(x, trait, mode = c("single", "joint"), idx = NULL) {
        mode <- arg_match(mode)
        t <- x@traits[[trait]]
        if (is.null(t)) {
            return(NULL)
        }
        if (mode == "joint") {
            return(t$joint)
        }
        if (is.null(idx)) {
            return(t$single)
        }
        if (idx > length(t$single)) {
            return(NULL)
        }
        t$single[[idx]]
    }
)

#' @rdname SldscData
setMethod("show", "SldscData", function(object) {
    cat("SldscData\n")
    cat(glue(
        "  annotations ({length(getAnnotCols(object))}): ",
        "{str_flatten(getAnnotCols(object), ', ')}\n",
        .trim = FALSE
    ))
    cat(glue(
        "  annot SNPs: {nrow(object@annot)} | ",
        "frq SNPs: {nrow(object@frq)}\n",
        .trim = FALSE
    ))
    cat(glue(
        "  traits ({length(object@traits)}): ",
        "{str_flatten(names(object@traits), ', ')}\n",
        .trim = FALSE
    ))
    invisible(object)
})
