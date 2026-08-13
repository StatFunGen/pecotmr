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
#   getMarginalEffects(x, maxPval = NULL) ......... marginal view (p-value
#   filter)
#   getPip(x), getCs(x, coverage)
#   adjustPips(x, keepVariants)
# =============================================================================

#' @include AllGenerics.R
NULL

#' @title Fine-Mapping Entry (per-tuple payload)
#' @description S4 container for a single fine-mapping fit attached to a
#'   \code{FineMappingResult} row. One entry corresponds to one \code{(study,
#'   context, trait, method)} tuple.
#'
#' For joint fits (e.g., multi-trait mvSuSiE or fSuSiE), multiple
#' \code{FineMappingEntry} objects in the same \code{FineMappingResult}
#' collection may carry references to the same underlying R fit object (R's
#' copy-on-modify semantics keep this memory-efficient).
#' @slot variantIds Character vector of variant IDs in the fit.
#' @slot susieFit The method-specific fit object (SuSiE list, mvSuSiE object,
#'   fSuSiE object, etc.). May be shared by reference across joint-fit entries.
#' @slot topLoci A long-format \code{data.frame} with at minimum
#'   \code{variant_id} and \code{pip} columns; optional \code{cs},
#'   \code{coverage}, \code{betahat}, \code{sd}, \code{csLog10bf}, \code{z}.
#' @slot cvResult Optional cross-validation payload (list with
#'   \code{samplePartition}, \code{predictions}, \code{performance}) recorded
#'   when fine-mapping is run with \code{cvFolds > 1}. \code{NULL} otherwise.
#' @export
setClass(
    "FineMappingEntry",
    representation(
        variantIds = "character",
        susieFit = "ANY",
        topLoci = "data.frame",
        cvResult = "ANY"
    ),
    validity = function(object) .validateFineMappingEntry(object)
)

# ---- FineMappingEntry validity helpers -------------------------------------

# Full validity check: collect all contract violations (empty vector = valid).
# @noRd
.validateFineMappingEntry <- function(object) {
    errors <- c(
        .fmeCheckCvResult(object),
        .fmeCheckTopLoci(object)
    )
    if (length(errors) == 0L) TRUE else errors
}

# cvResult must be NULL or a list.
# @noRd
.fmeCheckCvResult <- function(object) {
    if (!is.null(object@cvResult) && !is.list(object@cvResult)) {
        return(paste0(
            "cvResult must be NULL or a list ",
            "(samplePartition/predictions/performance)"
        ))
    }
    NULL
}

# topLoci contract checks (only meaningful when topLoci has rows).
# @noRd
.fmeCheckTopLoci <- function(object) {
    if (nrow(object@topLoci) == 0L) {
        return(NULL)
    }
    c(
        .fmeCheckTopLociCols(object),
        .fmeCheckTopLociRows(object),
        .fmeCheckTopLociOrder(object),
        .fmeCheckPipDrift(object)
    )
}

# Minimal column contract: variant_id + pip. Canonical projector columns
# (marginal_*, posterior_*, ...) are pipeline-populated; skeletal entries may
# omit them (accessor projections then return NA-filled cols).
# @noRd
.fmeCheckTopLociCols <- function(object) {
    missingCols <- setdiff(c("variant_id", "pip"), colnames(object@topLoci))
    if (length(missingCols) > 0L) {
        return(paste(
            "topLoci missing required columns:",
            paste(missingCols, collapse = ", ")
        ))
    }
    NULL
}

# topLoci row count must match variantIds length.
# @noRd
.fmeCheckTopLociRows <- function(object) {
    n <- length(object@variantIds)
    if (n > 0L && nrow(object@topLoci) != n) {
        return(sprintf(
            "topLoci has %d rows but variantIds has %d entries; ",
            nrow(object@topLoci),
            n
        ))
    }
    NULL
}

# topLoci$variant_id must equal variantIds in order.
# @noRd
.fmeCheckTopLociOrder <- function(object) {
    n <- length(object@variantIds)
    if (n == 0L || nrow(object@topLoci) != n) {
        return(NULL)
    }
    sameOrder <- identical(
        as.character(object@topLoci$variant_id),
        as.character(object@variantIds)
    )
    if (!sameOrder) {
        return("topLoci$variant_id must equal variantIds in order")
    }
    NULL
}

# Drift check: a susieFit$pip vector must match the topLoci pip column. Catches
# adjustPips() (or any mutator) updating one and forgetting the other.
# @noRd
.fmeCheckPipDrift <- function(object) {
    n <- length(object@variantIds)
    sf <- object@susieFit
    carriesPip <- !is.null(sf) &&
        is.list(sf) &&
        !is.null(sf$pip) &&
        length(sf$pip) == n &&
        "pip" %in% colnames(object@topLoci)
    if (!carriesPip) {
        return(NULL)
    }
    inSync <- isTRUE(all.equal(
        as.numeric(sf$pip),
        as.numeric(object@topLoci$pip),
        tolerance = 1e-10
    ))
    if (!inSync) {
        return("susieFit$pip and topLoci$pip have drifted out of sync")
    }
    NULL
}

#' @title Create a FineMappingEntry Object
#' @description Construct a \code{FineMappingEntry} payload for one
#'   \code{(study, context, trait, method)} row of a \code{FineMappingResult}
#'   collection.
#' @param variantIds Character vector of variant IDs in fit order.
#' @param susieFit The SuSiE fit object (full or trimmed; controlled by the
#'   pipeline's \code{trim} parameter).
#' @param topLoci Per-variant \code{data.frame} in canonical schema: identity
#'   columns (\code{variant_id, chrom, pos, A1, A2}), context (\code{N, af};
#'   effect-allele frequency, never MAF), marginal columns (\code{marginal_beta,
#'   marginal_se, marginal_z, marginal_p}), posterior columns (\code{pip,
#'   posterior_mean, posterior_sd, cs_*, cs_*_purity}), pipeline stamps
#'   (\code{method, gene, event, grange_start, grange_end}). Unfiltered: one row
#'   per variant in the fit.
#' @param cvResult Optional cross-validation payload (list with
#'   \code{samplePartition}, \code{predictions}, \code{performance}) recorded
#'   when fine-mapping is run with \code{cvFolds > 1}. \code{NULL} otherwise.
#' @return A \code{FineMappingEntry} object.
#' @examples
#' tl <- data.frame(variant_id = paste0("chr1:", 100 * 1:3, ":A:G"),
#'   pip = c(0.9, 0.5, 0.1), cs = c(1L, 1L, NA))
#' FineMappingEntry(variantIds = tl$variant_id, susieFit = list(), topLoci = tl)
#' @export
FineMappingEntry <- function(variantIds, susieFit, topLoci, cvResult = NULL) {
    obj <- new(
        "FineMappingEntry",
        variantIds = as.character(variantIds),
        susieFit = susieFit,
        topLoci = as.data.frame(topLoci),
        cvResult = cvResult
    )
    validObject(obj)
    obj
}

# =============================================================================
# Accessors
# =============================================================================

#' @rdname getVariantIds
#' @export
setMethod("getVariantIds", "FineMappingEntry", function(x, ...) x@variantIds)

# cTWAS weight-source accessors for a FineMappingEntry. resolveWeights returns
# the topLoci posterior effect colSums(alpha * mu) (buildTopLoci), which is on
# the STANDARDIZED scale -- unlike susieWeights / the cTWAS renorm, it does NOT
# divide by X_column_scale_factors. So:
#   getStandardized -> TRUE   skip cTWAS's per-variant variance scaling, which
#                             would double-standardize the effect (RSS and
#                             individual fits alike -- posterior_mean is always
#                             colSums(alpha*mu)).
# getFits -> NULL skip cTWAS's alpha renormalisation, which recomputes an
# UNstandardized weight (/ scale factors) that is inconsistent with the
# standardized posterior effect. This lands on the same scale as the
# susie-TwasWeights path: an unstandardized susie weight (mu / scale) x
# sqrt(var) = mu = the standardized posterior effect. (getSusieFit still exposes
# the fit for non-weight uses.)
#' @rdname getFits
#' @export
setMethod("getFits", "FineMappingEntry", function(x, ...) NULL)

#' @rdname getStandardized
#' @export
setMethod("getStandardized", "FineMappingEntry", function(x, ...) TRUE)

#' @rdname resolveWeights
#' @export
setMethod("resolveWeights", "FineMappingEntry", function(x, ...) {
    empty <- list(variantIds = character(0), weights = numeric(0))
    # getTopLoci projects the posterior effect to `beta` (= posterior_mean);
    # use it as the per-variant weight, aligned with variant_id.
    tl <- getTopLoci(x)
    if (
        is.null(tl) ||
            nrow(tl) == 0L ||
            !all(c("variant_id", "beta") %in% names(tl))
    ) {
        return(empty)
    }
    vids <- as.character(tl$variant_id)
    w <- as.numeric(tl$beta)
    ok <- !is.na(vids) & !is.na(w)
    if (!any(ok)) {
        return(empty)
    }
    list(variantIds = vids[ok], weights = w[ok])
})

#' @rdname getSusieFit
#' @export
setMethod("getSusieFit", "FineMappingEntry", function(x, ...) x@susieFit)

#' @rdname getCvResult
#' @export
setMethod("getCvResult", "FineMappingEntry", function(x, ...) x@cvResult)

#' @rdname getTopLoci
#' @export
setMethod(
    "getTopLoci",
    "FineMappingEntry",
    function(
        x,
        type = c("data.frame", "GRanges"),
        signalCutoff = 0.025,
        minPurity = NULL,
        ...
    ) {
        type <- match.arg(type)
        out <- .fmeFilterTopLoci(x@topLoci, signalCutoff, minPurity)
        if (type == "data.frame") {
            return(out)
        }
        .fmeTopLociGRanges(out)
    }
)

# Apply the pip `signalCutoff` + optional purity filter, then project to the
# posterior view. An empty topLoci passes through unchanged.
# @noRd
.fmeFilterTopLoci <- function(tl, signalCutoff, minPurity) {
    if (nrow(tl) == 0L) {
        return(tl)
    }
    keep <- if (is.null(signalCutoff) || signalCutoff <= 0) {
        rep(TRUE, nrow(tl))
    } else {
        !is.na(tl$pip) & tl$pip > signalCutoff
    }
    if (!is.null(minPurity)) {
        keep <- keep & .fmePurityKeep(tl, minPurity)
    }
    .projectPosteriorView(tl[keep, , drop = FALSE])
}

# Purity keep-mask (orthogonal to the pip cutoff): drop variants belonging to a
# below-minPurity primary (highest-coverage) credible set. Variants not in a CS
# are unaffected. All-TRUE when the entry carries no purity columns.
# @noRd
.fmePurityKeep <- function(tl, minPurity) {
    purCols <- grep("^cs_[0-9.]+_purity$", names(tl), value = TRUE)
    if (length(purCols) == 0L) {
        return(rep(TRUE, nrow(tl)))
    }
    covs <- suppressWarnings(as.numeric(sub(
        "_purity$",
        "",
        sub("^cs_", "", purCols)
    )))
    primaryPur <- purCols[which.max(covs)]
    csCol <- sub("_purity$", "", primaryPur)
    inCs <- if (csCol %in% names(tl)) {
        !grepl("_0$", tl[[csCol]])
    } else {
        rep(FALSE, nrow(tl))
    }
    pur <- as.numeric(tl[[primaryPur]])
    !(inCs & !is.na(pur) & pur < minPurity)
}

# Project a posterior-view data.frame to GRanges (empty GRanges when empty).
# @noRd
.fmeTopLociGRanges <- function(out) {
    if (is.null(out) || nrow(out) == 0L) {
        return(GenomicRanges::GRanges())
    }
    parsed <- parseVariantId(out$variant_id)
    gr <- GenomicRanges::GRanges(
        seqnames = paste0("chr", parsed$chrom),
        ranges = IRanges::IRanges(start = parsed$pos, width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(out)
    gr
}

#' @rdname getMarginalEffects
#' @export
setMethod(
    "getMarginalEffects",
    "FineMappingEntry",
    function(x, maxPval = NULL, ...) {
        tl <- x@topLoci
        if (nrow(tl) == 0L) {
            return(.projectMarginalView(tl))
        }
        out <- .projectMarginalView(tl)
        if (!is.null(maxPval) && nrow(out) > 0L) {
            keep <- !is.na(out$p) & out$p <= maxPval
            out <- out[keep, , drop = FALSE]
        }
        out
    }
)

#' @rdname getPip
#' @export
setMethod("getPip", "FineMappingEntry", function(x, ...) {
    tl <- x@topLoci
    if (nrow(tl) == 0L || !"pip" %in% names(tl)) {
        return(numeric(0))
    }
    setNames(tl$pip, tl$variant_id)
})

#' @rdname getCs
#' @export
setMethod(
    "getCs",
    "FineMappingEntry",
    function(x, coverage = 0.95, minPurity = NULL, ...) {
        tl <- x@topLoci
        if (nrow(tl) == 0L) {
            return(.projectPosteriorView(tl))
        }
        csCol <- grep(
            paste0("^cs_", coverage * 100, "$"),
            names(tl),
            value = TRUE
        )
        if (length(csCol) == 0L) {
            return(.projectPosteriorView(tl[FALSE, , drop = FALSE]))
        }
        keep <- !is.na(tl[[csCol[1L]]]) &
            nzchar(tl[[csCol[1L]]]) &
            !grepl("_0$", tl[[csCol[1L]]])
        # Independent purity filter (min.abs.corr), orthogonal to `coverage`:
        # drop CS members whose credible set at THIS coverage is below
        # `minPurity`.
        if (!is.null(minPurity)) {
            purCol <- paste0(csCol[1L], "_purity")
            if (purCol %in% names(tl)) {
                pur <- as.numeric(tl[[purCol]])
                keep <- keep & !is.na(pur) & pur >= minPurity
            } else {
                warning(
                    "getCs: no purity column '",
                    purCol,
                    "' for coverage ",
                    coverage,
                    "; minPurity filter skipped."
                )
            }
        }
        .projectPosteriorView(tl[keep, , drop = FALSE])
    }
)

#' @rdname adjustPips
#' @export
setMethod("adjustPips", "FineMappingEntry", function(x, keepVariants, ...) {
    keepIdx <- .adjustPipsKeepIdx(x@variantIds, keepVariants)
    common <- x@variantIds[keepIdx]
    fit <- .adjustPipsSubsetFit(x@susieFit, keepIdx)
    newTopLoci <- .adjustPipsRebuildTopLoci(x@topLoci, fit, common)
    # cvResult is sample-indexed (partition + held-out predictions), not
    # variant-indexed, so it carries through a variant-subsetting unchanged.
    new(
        "FineMappingEntry",
        variantIds = common,
        susieFit = fit,
        topLoci = newTopLoci,
        cvResult = x@cvResult
    )
})

# Match the entry's variants to `keepVariants` by (chrom, pos, allele) so a
# chr-prefix / separator difference does not read as no-overlap. PIP / LBF are
# coding-invariant (no sign). The index is sorted to keep the entry's variant
# order so fit-matrix subsetting is unchanged.
# @noRd
.adjustPipsKeepIdx <- function(variantIds, keepVariants) {
    m <- matchVariants(variantIds, as.character(keepVariants))
    if (!length(m$idxA)) {
        stop(
            "adjustPips: intersection of entry variants with `keepVariants` ",
            "is empty."
        )
    }
    sort(m$idxA)
}

# Subset the SuSiE fit to `keepIdx` and renormalize PIPs from the retained
# lbf_variable columns.
# @noRd
.adjustPipsSubsetFit <- function(fit, keepIdx) {
    if (is.null(fit$lbf_variable)) {
        stop(
            "adjustPips: entry's susieFit has no `lbf_variable` matrix; ",
            "PIP renormalization requires lbf_variable. Re-run the ",
            "pipeline with trim = FALSE to retain it."
        )
    }
    lbfSub <- fit$lbf_variable[, keepIdx, drop = FALSE]
    fit$lbf_variable <- lbfSub
    fit$alpha <- lbfToAlpha(lbfSub)
    fit$pip <- as.numeric(1 - apply(1 - fit$alpha, 2, prod))
    fit$mu <- .adjustPipsSubsetCols(fit$mu, keepIdx)
    fit$mu2 <- .adjustPipsSubsetCols(fit$mu2, keepIdx)
    if (!is.null(fit$X_column_scale_factors)) {
        fit$X_column_scale_factors <- fit$X_column_scale_factors[keepIdx]
    }
    fit
}

# Column-subset a (possibly 3-D) per-variant fit matrix; NULL passes through.
# @noRd
.adjustPipsSubsetCols <- function(m, keepIdx) {
    if (is.null(m)) {
        return(NULL)
    }
    if (length(dim(m)) == 3) {
        m[, keepIdx, , drop = FALSE]
    } else {
        m[, keepIdx, drop = FALSE]
    }
}

# Rebuild topLoci from the subset fit + existing (per-variant) marginal columns.
# @noRd
.adjustPipsRebuildTopLoci <- function(topLoci, fit, common) {
    if (nrow(topLoci) == 0L) {
        return(topLoci)
    }
    newTopLoci <- topLoci[topLoci$variant_id %in% common, , drop = FALSE]
    newTopLoci$pip <- as.numeric(fit$pip)
    .adjustPipsPosterior(newTopLoci, fit)
}

# Recompute posterior_mean / posterior_sd from the fit when alpha + mu/mu2 are
# matrix-shaped and conformable; otherwise leave the existing values in place.
# @noRd
.adjustPipsPosterior <- function(newTopLoci, fit) {
    alphaMat <- if (!is.null(fit$alpha)) as.matrix(fit$alpha) else NULL
    muMat <- if (!is.null(fit$mu)) as.matrix(fit$mu) else NULL
    if (
        is.null(alphaMat) ||
            is.null(muMat) ||
            !all(dim(alphaMat) == dim(muMat))
    ) {
        return(newTopLoci)
    }
    newTopLoci$posterior_mean <- as.numeric(colSums(alphaMat * muMat))
    mu2Mat <- if (!is.null(fit$mu2)) as.matrix(fit$mu2) else NULL
    if (!is.null(mu2Mat) && all(dim(alphaMat) == dim(mu2Mat))) {
        newTopLoci$posterior_sd <- as.numeric(sqrt(pmax(
            colSums(alphaMat * mu2Mat) - newTopLoci$posterior_mean^2,
            0
        )))
    }
    newTopLoci
}

#' @rdname show-methods
#' @export
setMethod("show", "FineMappingEntry", function(object) {
    tl <- object@topLoci
    nCs <- if (nrow(tl) > 0L) {
        csCols <- grep("^cs_[0-9]+$", names(tl), value = TRUE)
        if (length(csCols) > 0L) {
            vals <- unique(unlist(lapply(csCols, function(cc) {
                v <- tl[[cc]]
                v <- v[!grepl("_0$", v)]
                v
            })))
            length(vals)
        } else {
            0L
        }
    } else {
        0L
    }
    cat(sprintf(
        "FineMappingEntry: %d variants, %d credible sets\n",
        length(object@variantIds),
        nCs
    ))
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
        return(switch(
            type,
            character = as.character(tl[[name]]),
            integer = as.integer(tl[[name]]),
            numeric = as.numeric(tl[[name]])
        ))
    }
    switch(
        type,
        character = rep(NA_character_, nrow(tl)),
        integer = rep(NA_integer_, nrow(tl)),
        numeric = rep(NA_real_, nrow(tl))
    )
}

# Project the canonical wide topLoci to the posterior view: identity +
# N/af + (beta=posterior_mean, se=posterior_sd) + pip + cs_* + signal_cluster
# + pipeline stamps. Renames `posterior_mean`/`posterior_sd` to `beta`/`se`.
# Exports effect-allele frequency as `af` (never MAF). Missing optional
# columns are NA-filled.
# @noRd
.projectPosteriorView <- function(tl) {
    if (nrow(tl) == 0L) {
        return(.projectPosteriorEmpty())
    }
    out <- data.frame(
        variant_id = .tlCol(tl, "variant_id", "character"),
        chrom = .tlCol(tl, "chrom", "character"),
        pos = .tlCol(tl, "pos", "integer"),
        A1 = .tlCol(tl, "A1", "character"),
        A2 = .tlCol(tl, "A2", "character"),
        N = .tlCol(tl, "N", "numeric"),
        af = .tlCol(tl, "af", "numeric"),
        beta = .tlCol(tl, "posterior_mean", "numeric"),
        se = .tlCol(tl, "posterior_sd", "numeric"),
        pip = .tlCol(tl, "pip", "numeric"),
        logBF = .tlCol(tl, "logBF", "numeric"),
        stringsAsFactors = FALSE
    )
    for (cc in .projectPosteriorExtraCols(tl)) {
        out[[cc]] <- tl[[cc]]
    }
    rownames(out) <- NULL
    out
}

# Empty posterior-view frame (canonical columns, zero rows).
# @noRd
.projectPosteriorEmpty <- function() {
    data.frame(
        variant_id = character(0),
        chrom = character(0),
        pos = integer(0),
        A1 = character(0),
        A2 = character(0),
        N = numeric(0),
        af = numeric(0),
        beta = numeric(0),
        se = numeric(0),
        pip = numeric(0),
        logBF = numeric(0),
        stringsAsFactors = FALSE
    )
}

# Optional pass-through columns for the posterior view, in canonical order: the
# dynamic CS columns (cs_<C> membership / cs_<C>_purity), the per-CS variant-
# level fullFit columns (within_cs_pip[_<lab>] / cs_logbf_<lab> /
# cs_effect[_var]_<lab>), the per-condition posterior quantities
# (conditional_effect, lfsr), and pipeline stamps -- whichever are present.
# @noRd
.projectPosteriorExtraCols <- function(tl) {
    csCols <- grep("^cs_[0-9.]+(_purity)?$", colnames(tl), value = TRUE)
    fullFitCols <- grep(
        "^(within_cs_pip|cs_logbf_|cs_effect_)",
        colnames(tl),
        value = TRUE
    )
    intersect(
        c(
            csCols,
            fullFitCols,
            "conditional_effect",
            "lfsr",
            "method",
            "gene",
            "event",
            "grange_start",
            "grange_end"
        ),
        colnames(tl)
    )
}

# Project to the marginal view: identity + N/af + (beta, se, z, p) where
# beta/se/z/p are the marginal univariate columns renamed from their
# `marginal_*` storage names. Exports effect-allele frequency as `af`
# (never MAF). Missing optional columns are NA-filled.
# @noRd
.projectMarginalView <- function(tl) {
    if (nrow(tl) == 0L) {
        return(data.frame(
            variant_id = character(0),
            chrom = character(0),
            pos = integer(0),
            A1 = character(0),
            A2 = character(0),
            N = numeric(0),
            af = numeric(0),
            beta = numeric(0),
            se = numeric(0),
            z = numeric(0),
            p = numeric(0),
            stringsAsFactors = FALSE
        ))
    }
    data.frame(
        variant_id = .tlCol(tl, "variant_id", "character"),
        chrom = .tlCol(tl, "chrom", "character"),
        pos = .tlCol(tl, "pos", "integer"),
        A1 = .tlCol(tl, "A1", "character"),
        A2 = .tlCol(tl, "A2", "character"),
        N = .tlCol(tl, "N", "numeric"),
        af = .tlCol(tl, "af", "numeric"),
        beta = .tlCol(tl, "marginal_beta", "numeric"),
        se = .tlCol(tl, "marginal_se", "numeric"),
        z = .tlCol(tl, "marginal_z", "numeric"),
        p = .tlCol(tl, "marginal_p", "numeric"),
        stringsAsFactors = FALSE
    )
}
