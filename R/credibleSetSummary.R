#' @include AllGenerics.R AllClasses.R FineMappingEntry.R tupleSelectors.R
NULL

.emptyCsSummary <- function() {
    data.frame(
        cs = character(),
        effect_id = character(),
        coverage = numeric(),
        n_variants = integer(),
        purity_min = numeric(),
        purity_mean = numeric(),
        V = numeric(),
        cs_log10bf = numeric(),
        cs_log_bf = numeric(),
        cs_pip = numeric(),
        cs_mean_effect = numeric(),
        lead_variant = character(),
        lead_pip = numeric(),
        stringsAsFactors = FALSE
    )
}

# The per-effect single-effect log Bayes factor: log of the average (over the
# flat prior) of the effect's per-variant Bayes factors, i.e.
# logSumExp(lbf_variable[L, ]) - log(p). This is the principled per-CS evidence,
# unlike the max-member per-variant logBF. NA when the fit carries no lbf
# matrix.
.csEffectLogBf <- function(fit, Lidx) {
    lbf <- fit$lbf_variable
    if (is.null(lbf) || !is.matrix(lbf) || is.na(Lidx) || Lidx > nrow(lbf)) {
        return(NA_real_)
    }
    lv <- as.numeric(lbf[Lidx, ])
    lv <- lv[is.finite(lv)]
    if (length(lv) == 0L) {
        return(NA_real_)
    }
    mx <- max(lv)
    mx + log(sum(exp(lv - mx))) - log(length(lv))
}

# One row per credible set at `coverage`, summarising size, purity (min from the
# canonical cs_<cov>_purity column, mean from the fit's sets$purity when
# present),
# the per-effect prior variance V, the CS log Bayes factor (max member logBF),
# and the lead (highest-PIP) variant. Sourced primarily from the entry's topLoci
# so it stays consistent with getCs / getTopLoci; V / mean-purity come from the
# stored fit. The credible-set index in the `cs` label doubles as the effect (L)
# index for susie-family fits (unfiltered CS are ordered by effect).
# @noRd
.csSummaryFit <- function(tl, fit, coverage) {
    csCol <- paste0("cs_", coverage * 100)
    purCol <- paste0(csCol, "_purity")
    if (is.null(tl) || nrow(tl) == 0L || !(csCol %in% names(tl))) {
        return(.emptyCsSummary())
    }
    labels <- unique(tl[[csCol]][
        !is.na(tl[[csCol]]) & nzchar(tl[[csCol]]) & !grepl("_0$", tl[[csCol]])
    ])
    if (length(labels) == 0L) {
        return(.emptyCsSummary())
    }
    ctx <- .csSummaryContext(tl, fit, coverage, csCol, purCol)
    out <- do.call(rbind, map(labels, function(lab) .csSummaryRow(lab, ctx)))
    out$cs_log10bf[!is.finite(out$cs_log10bf)] <- NA_real_
    out
}

# Per-fit context for the CS summary rows: prior variances (V) and per-effect
# mean-purity, plus the resolved column names.
# @noRd
.csSummaryContext <- function(tl, fit, coverage, csCol, purCol) {
    Vvec <- if (!is.null(fit$V)) as.numeric(fit$V) else NULL
    meanPur <- NULL
    sp <- fit$sets$purity
    if (!is.null(sp)) {
        mc <- intersect(c("meanAbsCorr", "mean.abs.corr"), names(sp))
        if (length(mc) > 0L) {
            meanPur <- stats::setNames(as.numeric(sp[[mc[1L]]]), rownames(sp))
        }
    }
    list(
        tl = tl,
        fit = fit,
        coverage = coverage,
        csCol = csCol,
        purCol = purCol,
        Vvec = Vvec,
        meanPur = meanPur
    )
}

# One credible-set summary row for label `lab`.
# @noRd
.csSummaryRow <- function(lab, ctx) {
    m <- ctx$tl[ctx$tl[[ctx$csCol]] == lab, , drop = FALSE]
    Lidx <- suppressWarnings(as.integer(sub("^.*_", "", lab)))
    Lname <- if (!is.na(Lidx)) paste0("L", Lidx) else NA_character_
    lead <- which.max(m$pip)
    data.frame(
        cs = lab,
        effect_id = Lname,
        coverage = ctx$coverage,
        n_variants = nrow(m),
        purity_min = .csPurityMin(m, ctx$purCol),
        purity_mean = .csPurityMean(ctx$meanPur, Lname),
        V = .csEffectV(ctx$Vvec, Lidx),
        # cs_log10bf: the strongest MEMBER variant's logBF (max over the CS).
        cs_log10bf = .csLog10Bf(m),
        # cs_log_bf: the true per-EFFECT single-effect log Bayes factor.
        cs_log_bf = .csEffectLogBf(ctx$fit, Lidx),
        # cs_pip: total posterior inclusion mass captured by the CS.
        cs_pip = suppressWarnings(sum(as.numeric(m$pip), na.rm = TRUE)),
        # cs_mean_effect: mean posterior conditional effect over the CS (NA when
        # the entry has no conditional_effect column, e.g. univariate susie).
        cs_mean_effect = .csMeanEffect(m),
        lead_variant = as.character(m$variant_id[lead]),
        lead_pip = as.numeric(m$pip[lead]),
        stringsAsFactors = FALSE
    )
}

# purity_min: first per-CS purity value (NA when the column is absent).
# @noRd
.csPurityMin <- function(m, purCol) {
    if (purCol %in% names(m)) as.numeric(m[[purCol]][1]) else NA_real_
}

# purity_mean: per-effect mean absolute correlation (NA when unavailable).
# @noRd
.csPurityMean <- function(meanPur, Lname) {
    if (!is.null(meanPur) && !is.na(Lname) && Lname %in% names(meanPur)) {
        as.numeric(meanPur[[Lname]])
    } else {
        NA_real_
    }
}

# V: the effect's prior variance (NA when out of range / unavailable).
# @noRd
.csEffectV <- function(Vvec, Lidx) {
    if (!is.null(Vvec) && !is.na(Lidx) && Lidx <= length(Vvec)) {
        Vvec[Lidx]
    } else {
        NA_real_
    }
}

# cs_log10bf: strongest member logBF (NA when the column is absent).
# @noRd
.csLog10Bf <- function(m) {
    if ("logBF" %in% names(m)) {
        suppressWarnings(max(m$logBF, na.rm = TRUE))
    } else {
        NA_real_
    }
}

# cs_mean_effect: mean conditional effect over the CS (NA when absent).
# @noRd
.csMeanEffect <- function(m) {
    if ("conditional_effect" %in% names(m)) {
        suppressWarnings(mean(as.numeric(m$conditional_effect), na.rm = TRUE))
    } else {
        NA_real_
    }
}

#' Per-credible-set summary of a fine-mapping result
#'
#' One row per credible set (at a given coverage) with its size, purity, prior
#' variance, log Bayes factor, and lead variant -- the per-CS complement to the
#' per-variant \code{\link{getCs}}. Replaces the legacy per-effect `effect.tsv`.
#'
#' @param x A \code{FineMappingEntry} or \code{FineMappingResultBase}.
#' @param coverage Credible-set coverage to summarise. Default 0.95.
#' @param ... Ignored.
#' @return A \code{data.frame}: \code{cs, effect_id, coverage, n_variants,
#'   purity_min, purity_mean, V, cs_log10bf} (strongest member logBF),
#'   \code{cs_log_bf} (true per-effect single-effect log Bayes factor),
#'   \code{cs_pip} (summed member PIP = inclusion mass captured),
#'   \code{cs_mean_effect} (mean posterior conditional effect),
#'   \code{lead_variant, lead_pip} (the collection method also carries the entry
#'   identity columns).
#' @seealso \code{\link{getCs}}, \code{\link{getTopLoci}}
#' @examples
#' data(qtlFineMappingExample)
#' getCredibleSetSummary(qtlFineMappingExample)
#' @export
setGeneric("getCredibleSetSummary", function(x, ...) {
    standardGeneric("getCredibleSetSummary")
})

#' @rdname getCredibleSetSummary
#' @export
setMethod(
    "getCredibleSetSummary",
    "FineMappingEntry",
    function(x, coverage = 0.95, ...) {
        .csSummaryFit(x@topLoci, getSusieFit(x), coverage)
    }
)

#' @rdname getCredibleSetSummary
#' @export
setMethod(
    "getCredibleSetSummary",
    "FineMappingResultBase",
    function(x, coverage = 0.95, ...) {
        .fmrAggregateView(x, perEntry = function(e) {
            getCredibleSetSummary(e, coverage = coverage)
        })
    }
)
