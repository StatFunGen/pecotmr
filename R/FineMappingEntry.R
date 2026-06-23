# =============================================================================
# FineMappingEntry S4 class
# -----------------------------------------------------------------------------
# Per-tuple fine-mapping payload backing one row of a FineMappingResult
# collection. Three slots:
#
#   variantIds : character vector, variant IDs in fit order
#   susieFit   : the SuSiE fit (full or trimmed; controlled by the
#                pipeline's `trim` parameter)
#   topLoci    : unfiltered per-variant data.frame carrying BOTH marginal
#                univariate effects and posterior fine-mapping output in
#                a single wide table. Stored in canonical schema (column
#                names with `marginal_*` / `posterior_*` prefixes);
#                accessors project + rename to user-facing column names.
#
# Accessors:
#   getVariantIds(x)
#   getSusieFit(x)
#   getTopLoci(x, signalCutoff = 0.025) ........... posterior view (PIP filter)
#   getMarginalEffects(x, maxPval = NULL) ......... marginal view (p-value filter)
#   getPip(x), getCs(x, coverage)
#   adjustPips(x, keepVariants)
# =============================================================================

#' @include AllGenerics.R
NULL

setClass("FineMappingEntry",
  representation(
    variantIds = "character",
    susieFit   = "ANY",
    topLoci    = "data.frame"
  ),
  validity = function(object) {
    errors <- character()
    n <- length(object@variantIds)
    if (nrow(object@topLoci) > 0L) {
      # Minimal contract: variant_id + pip. Canonical projector columns
      # (marginal_*, posterior_*, etc.) are pipeline-populated; tests
      # and downstream consumers building skeletal entries can omit
      # them, in which case accessor projections return NA-filled cols.
      required <- c("variant_id", "pip")
      missingCols <- setdiff(required, colnames(object@topLoci))
      if (length(missingCols) > 0L)
        errors <- c(errors,
          paste("topLoci missing required columns:",
                paste(missingCols, collapse = ", ")))
      if (n > 0L && nrow(object@topLoci) != n)
        errors <- c(errors,
          sprintf("topLoci has %d rows but variantIds has %d entries; ",
                  nrow(object@topLoci), n))
      if (n > 0L && nrow(object@topLoci) == n &&
          !identical(as.character(object@topLoci$variant_id),
                     as.character(object@variantIds)))
        errors <- c(errors,
          "topLoci$variant_id must equal variantIds in order")
      # Drift check: if susieFit carries its own pip vector, it must match
      # the topLoci pip column. Catches the case where adjustPips() (or
      # any future mutator) updates one and forgets the other.
      sf <- object@susieFit
      if (!is.null(sf) && is.list(sf) && !is.null(sf$pip) &&
          length(sf$pip) == n && "pip" %in% colnames(object@topLoci)) {
        if (!isTRUE(all.equal(as.numeric(sf$pip),
                              as.numeric(object@topLoci$pip),
                              tolerance = 1e-10)))
          errors <- c(errors,
            "susieFit$pip and topLoci$pip have drifted out of sync")
      }
    }
    if (length(errors) == 0L) TRUE else errors
  }
)

#' @title Create a FineMappingEntry Object
#' @description Construct a \code{FineMappingEntry} payload for one
#'   \code{(study, context, trait, method)} row of a
#'   \code{FineMappingResult} collection.
#' @param variantIds Character vector of variant IDs in fit order.
#' @param susieFit The SuSiE fit object (full or trimmed; controlled by
#'   the pipeline's \code{trim} parameter).
#' @param topLoci Per-variant \code{data.frame} in canonical schema:
#'   identity columns (\code{variant_id, chrom, pos, A1, A2}), context
#'   (\code{N, af}; effect-allele frequency, never MAF), marginal columns (\code{marginal_beta,
#'   marginal_se, marginal_z, marginal_p}), posterior columns
#'   (\code{pip, posterior_mean, posterior_sd, cs_*, cs_*_purity}),
#'   pipeline stamps (\code{method, gene, event, grange_start,
#'   grange_end}). Unfiltered: one row per variant in the fit.
#' @return A \code{FineMappingEntry} object.
#' @export
FineMappingEntry <- function(variantIds, susieFit, topLoci) {
  obj <- new("FineMappingEntry",
             variantIds = as.character(variantIds),
             susieFit   = susieFit,
             topLoci    = as.data.frame(topLoci))
  validObject(obj)
  obj
}

# =============================================================================
# Accessors
# =============================================================================

#' @rdname getVariantIds
#' @export
setMethod("getVariantIds", "FineMappingEntry",
          function(x, ...) x@variantIds)

#' @rdname getSusieFit
#' @export
setMethod("getSusieFit", "FineMappingEntry",
          function(x, ...) x@susieFit)

#' @rdname getTopLoci
#' @export
setMethod("getTopLoci", "FineMappingEntry",
  function(x, type = c("data.frame", "GRanges"),
           signalCutoff = 0.025, ...) {
    type <- match.arg(type)
    tl <- x@topLoci
    if (nrow(tl) == 0L) {
      out <- tl
    } else {
      keep <- if (is.null(signalCutoff) || signalCutoff <= 0) {
        rep(TRUE, nrow(tl))
      } else {
        !is.na(tl$pip) & tl$pip > signalCutoff
      }
      out <- .projectPosteriorView(tl[keep, , drop = FALSE])
    }
    if (type == "data.frame") return(out)
    if (is.null(out) || nrow(out) == 0L) return(GenomicRanges::GRanges())
    parsed <- parseVariantId(out$variant_id)
    gr <- GenomicRanges::GRanges(
      seqnames = paste0("chr", parsed$chrom),
      ranges = IRanges::IRanges(start = parsed$pos, width = 1L))
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(out)
    gr
  })

#' @rdname getMarginalEffects
#' @export
setMethod("getMarginalEffects", "FineMappingEntry",
  function(x, maxPval = NULL, ...) {
    tl <- x@topLoci
    if (nrow(tl) == 0L) return(.projectMarginalView(tl))
    out <- .projectMarginalView(tl)
    if (!is.null(maxPval) && nrow(out) > 0L) {
      keep <- !is.na(out$p) & out$p <= maxPval
      out <- out[keep, , drop = FALSE]
    }
    out
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
    if (nrow(tl) == 0L) return(.projectPosteriorView(tl))
    csCol <- grep(paste0("^cs_", coverage * 100, "$"), names(tl), value = TRUE)
    if (length(csCol) == 0L) return(.projectPosteriorView(tl[FALSE, , drop = FALSE]))
    keep <- !is.na(tl[[csCol[1L]]]) & nzchar(tl[[csCol[1L]]]) &
            !grepl("_0$", tl[[csCol[1L]]])
    .projectPosteriorView(tl[keep, , drop = FALSE])
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
    fit <- x@susieFit
    if (is.null(fit$lbf_variable))
      stop("adjustPips: entry's susieFit has no `lbf_variable` matrix; ",
           "PIP renormalization requires lbf_variable. Re-run the ",
           "pipeline with trim = FALSE to retain it.")
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
    # Rebuild topLoci consistently from the new fit + the existing
    # marginal columns (which are per-variant and just need subsetting).
    newTopLoci <- x@topLoci
    if (nrow(newTopLoci) > 0L) {
      newTopLoci <- newTopLoci[newTopLoci$variant_id %in% common, ,
                               drop = FALSE]
      newTopLoci$pip <- as.numeric(fit$pip)
      # Posterior mean / SD computed from the fit when alpha + mu/mu2
      # are matrix-shaped. When either is missing or shapes don't
      # match, leave the existing column values in place.
      alphaMat <- if (!is.null(fit$alpha)) as.matrix(fit$alpha) else NULL
      muMat    <- if (!is.null(fit$mu))    as.matrix(fit$mu)    else NULL
      mu2Mat   <- if (!is.null(fit$mu2))   as.matrix(fit$mu2)   else NULL
      if (!is.null(alphaMat) && !is.null(muMat) &&
          all(dim(alphaMat) == dim(muMat))) {
        newTopLoci$posterior_mean <- as.numeric(colSums(alphaMat * muMat))
        if (!is.null(mu2Mat) && all(dim(alphaMat) == dim(mu2Mat))) {
          newTopLoci$posterior_sd <- as.numeric(sqrt(pmax(
            colSums(alphaMat * mu2Mat) - newTopLoci$posterior_mean^2, 0)))
        }
      }
    }
    new("FineMappingEntry",
        variantIds = common,
        susieFit   = fit,
        topLoci    = newTopLoci)
  })

#' @export
setMethod("show", "FineMappingEntry", function(object) {
  tl <- object@topLoci
  nCs <- if (nrow(tl) > 0L) {
    csCols <- grep("^cs_[0-9]+$", names(tl), value = TRUE)
    if (length(csCols) > 0L) {
      vals <- unique(unlist(lapply(csCols, function(cc) {
        v <- tl[[cc]]; v <- v[!grepl("_0$", v)]; v
      })))
      length(vals)
    } else 0L
  } else 0L
  cat(sprintf("FineMappingEntry: %d variants, %d credible sets\n",
              length(object@variantIds), nCs))
})

# =============================================================================
# Internal column projectors used by accessors
# =============================================================================

# Read a column from the canonical wide topLoci, returning NAs of the
# given type when the column is absent. Lets accessor projectors tolerate
# skeletal entries that lack optional columns.
# @noRd
.tlCol <- function(tl, name, type = c("character", "integer", "numeric")) {
  type <- match.arg(type)
  if (name %in% colnames(tl)) {
    return(switch(type,
      character = as.character(tl[[name]]),
      integer   = as.integer(tl[[name]]),
      numeric   = as.numeric(tl[[name]])))
  }
  switch(type,
    character = rep(NA_character_, nrow(tl)),
    integer   = rep(NA_integer_,   nrow(tl)),
    numeric   = rep(NA_real_,      nrow(tl)))
}

# Project the canonical wide topLoci to the posterior view: identity +
# N/af + (beta=posterior_mean, se=posterior_sd) + pip + cs_* + signal_cluster
# + pipeline stamps. Renames `posterior_mean`/`posterior_sd` to `beta`/`se`.
# Exports effect-allele frequency as `af` (never MAF). Missing optional
# columns are NA-filled.
# @noRd
.projectPosteriorView <- function(tl) {
  if (nrow(tl) == 0L) {
    return(data.frame(
      variant_id = character(0), chrom = character(0), pos = integer(0),
      A1 = character(0), A2 = character(0),
      N = numeric(0), af = numeric(0),
      beta = numeric(0), se = numeric(0), pip = numeric(0),
      stringsAsFactors = FALSE))
  }
  out <- data.frame(
    variant_id      = .tlCol(tl, "variant_id", "character"),
    chrom           = .tlCol(tl, "chrom",      "character"),
    pos             = .tlCol(tl, "pos",        "integer"),
    A1              = .tlCol(tl, "A1",         "character"),
    A2              = .tlCol(tl, "A2",         "character"),
    N               = .tlCol(tl, "N",          "numeric"),
    af              = .tlCol(tl, "af",         "numeric"),
    beta            = .tlCol(tl, "posterior_mean", "numeric"),
    se              = .tlCol(tl, "posterior_sd",   "numeric"),
    pip             = .tlCol(tl, "pip",        "numeric"),
    stringsAsFactors = FALSE)
  # Pass through CS columns (cs_95, cs_70, cs_50, cs_95_purity) and
  # pipeline stamps (method, gene, event, grange_*) when present.
  extraCols <- intersect(
    c("cs_95", "cs_70", "cs_50", "cs_95_purity",
      "method", "gene", "event", "grange_start", "grange_end"),
    colnames(tl))
  for (cc in extraCols) out[[cc]] <- tl[[cc]]
  rownames(out) <- NULL
  out
}

# Project to the marginal view: identity + N/af + (beta, se, z, p) where
# beta/se/z/p are the marginal univariate columns renamed from their
# `marginal_*` storage names. Exports effect-allele frequency as `af`
# (never MAF). Missing optional columns are NA-filled.
# @noRd
.projectMarginalView <- function(tl) {
  if (nrow(tl) == 0L) {
    return(data.frame(
      variant_id = character(0), chrom = character(0), pos = integer(0),
      A1 = character(0), A2 = character(0),
      N = numeric(0), af = numeric(0),
      beta = numeric(0), se = numeric(0), z = numeric(0), p = numeric(0),
      stringsAsFactors = FALSE))
  }
  data.frame(
    variant_id = .tlCol(tl, "variant_id",    "character"),
    chrom      = .tlCol(tl, "chrom",         "character"),
    pos        = .tlCol(tl, "pos",           "integer"),
    A1         = .tlCol(tl, "A1",            "character"),
    A2         = .tlCol(tl, "A2",            "character"),
    N          = .tlCol(tl, "N",             "numeric"),
    af         = .tlCol(tl, "af",            "numeric"),
    beta       = .tlCol(tl, "marginal_beta", "numeric"),
    se         = .tlCol(tl, "marginal_se",   "numeric"),
    z          = .tlCol(tl, "marginal_z",    "numeric"),
    p          = .tlCol(tl, "marginal_p",    "numeric"),
    stringsAsFactors = FALSE)
}
