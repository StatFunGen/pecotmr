#' @include AllGenerics.R AllClasses.R FineMappingEntry.R tupleSelectors.R
NULL

.emptyCsSummary <- function() {
    tibble(
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
        lead_pip = numeric()
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

# One row per credible set at `coverage`. Prefers the fit's own per-coverage
# credible sets (`susieFit$sets$cs`): sourcing from the fit makes the table
# COMPLETE -- every set the fit found gets a row, including one whose variants
# were all relabelled to a smaller overlapping set in the per-variant column --
# with per-set stats over the true full membership. When the stored fit carries
# no per-coverage sets (a minimal / hand-built fit), it falls back to
# enumerating the per-variant cs_<cov> column (the labelled members only). Rows
# carry the fit's true effect index (cs = <method>_<k>, effect_id = L<k>).
# @noRd
.csSummaryFit <- function(tl, fit, coverage) {
    specs <- .csSummarySpecs(tl, fit, coverage)
    if (length(specs) == 0L) {
        return(.emptyCsSummary())
    }
    ctx <- list(tl = tl, fit = fit, coverage = coverage)
    map(specs, .csSummaryRow, ctx = ctx) |>
        bind_rows() |>
        mutate(
            cs_log10bf = if_else(
                is.finite(.data$cs_log10bf),
                .data$cs_log10bf,
                NA_real_
            )
        )
}

# Per-CS specs (member variant ids, true effect index, purity, method tag) for
# `coverage`: from the fit's sets$cs when present (complete membership), else
# from the per-variant cs_<cov> column (the labelled members only).
# @noRd
.csSummarySpecs <- function(tl, fit, coverage) {
    setsObj <- .csSetsForCoverage(fit, coverage)
    if (!is.null(setsObj) && !is.null(setsObj$cs) && length(setsObj$cs) > 0L) {
        return(.csSpecsFromFit(setsObj, tl, fit, coverage))
    }
    .csSpecsFromColumn(tl, fit, coverage)
}

# Complete per-CS specs sourced from the fit's credible sets. Member positions
# (setsObj$cs) map to variant ids via the alpha column names; the effect index
# comes from the sets$cs "L<k>" names; purity is read by set position.
# @noRd
.csSpecsFromFit <- function(setsObj, tl, fit, coverage) {
    vn <- colnames(fit$alpha)
    eff <- .fmEffectIndices(setsObj$cs)
    methodTag <- .csMethodTag(tl, fit, coverage)
    map(seq_along(setsObj$cs), function(k) {
        members <- as.integer(setsObj$cs[[k]])
        list(
            ids = if (!is.null(vn)) vn[members] else NA_character_,
            n = length(members),
            eff = eff[k],
            purity_min = .csPurityAtK(setsObj$purity, k, "min.abs.corr"),
            purity_mean = .csPurityAtK(
                setsObj$purity, k, c("mean.abs.corr", "meanAbsCorr")
            ),
            methodTag = methodTag
        )
    })
}

# Fallback per-CS specs from the per-variant cs_<cov> column (labelled members
# only); used when the stored fit carries no per-coverage credible sets.
# @noRd
.csSpecsFromColumn <- function(tl, fit, coverage) {
    csCol <- str_c("cs_", round(coverage * 100))
    purCol <- str_c(csCol, "_purity")
    if (is.null(tl) || nrow(tl) == 0L || !is_in(csCol, names(tl))) {
        return(list())
    }
    vals <- tl[[csCol]]
    labels <- unique(vals[
        !is.na(vals) & str_length(vals) > 0L & !str_detect(vals, "_0$")
    ])
    if (length(labels) == 0L) {
        return(list())
    }
    meanPur <- .csMeanPurLookup(fit)
    map(labels, function(lab) {
        m <- tl[!is.na(vals) & vals == lab, , drop = FALSE]
        eff <- suppressWarnings(as.integer(str_remove(lab, "^.*_")))
        Lname <- if (!is.na(eff)) str_c("L", eff) else NA_character_
        list(
            ids = m$variant_id,
            n = nrow(m),
            eff = eff,
            purity_min = if (is_in(purCol, names(m))) {
                as.numeric(m[[purCol]][1])
            } else {
                NA_real_
            },
            purity_mean = if (!is.null(meanPur) && !is.na(Lname) &&
                is_in(Lname, names(meanPur))) {
                as.numeric(meanPur[[Lname]])
            } else {
                NA_real_
            },
            methodTag = str_remove(lab, "_[0-9]+$")
        )
    })
}

# The credible-set object (sets$cs + sets$purity) for `coverage`: the fit's
# primary `sets` when it matches the fit's requested coverage, else the stored
# secondary set for that coverage (`fit$sets_secondary`, named
# "CS_<cov*100>_<method>"). NULL when that coverage was not computed / stored.
# @noRd
.csSetsForCoverage <- function(fit, coverage) {
    pc <- fit$sets$requested_coverage
    if (!is.null(pc) &&
        isTRUE(all.equal(as.numeric(pc[1L]), as.numeric(coverage)))) {
        return(fit$sets)
    }
    sec <- fit$sets_secondary
    if (!is.null(sec) && length(sec) > 0L) {
        nums <- suppressWarnings(
            as.integer(str_match(names(sec), "(?i)cs_(\\d+)_")[, 2L])
        )
        h <- which(nums == round(coverage * 100))
        if (length(h) > 0L) {
            return(sec[[h[1L]]]$sets)
        }
    }
    # Untrimmed fit with no secondary store and no recorded coverage: the
    # primary sets are the only ones present -- use them.
    if ((is.null(sec) || length(sec) == 0L) &&
        is.null(pc) && !is.null(fit$sets$cs)) {
        return(fit$sets)
    }
    NULL
}

# The method tag for the cs label: prefer the tag already in the entry's
# cs_<cov> column (keeps the summary consistent with the per-variant labels);
# fall back to the fit's class when the column carries no CS label.
# @noRd
.csMethodTag <- function(tl, fit, coverage) {
    csCol <- str_c("cs_", round(coverage * 100))
    if (!is.null(tl) && is_in(csCol, names(tl))) {
        v <- tl[[csCol]]
        v <- v[!is.na(v) & str_length(v) > 0L & !str_detect(v, "_0$")]
        if (length(v) > 0L) {
            return(str_remove(v[1L], "_[0-9]+$"))
        }
    }
    .camelToSnakeMethod(class(fit)[1L])
}

# Per-effect mean-abs-corr lookup (named by L<k>) from the fit's sets$purity.
# @noRd
.csMeanPurLookup <- function(fit) {
    sp <- fit$sets$purity
    if (is.null(sp)) {
        return(NULL)
    }
    mc <- intersect(c("meanAbsCorr", "mean.abs.corr"), names(sp))
    if (length(mc) == 0L) {
        return(NULL)
    }
    set_names(as.numeric(sp[[mc[1L]]]), rownames(sp))
}

# One summary row from a per-CS spec. Size / effect index / purity / method tag
# come from the spec; member statistics (pip, member logBF, conditional effect,
# lead) from the entry topLoci; prior variance / per-effect logBF from the fit.
# @noRd
.csSummaryRow <- function(spec, ctx) {
    m <- .csMemberRows(ctx$tl, spec$ids)
    pipVec <- if (!is.null(m) && is_in("pip", names(m))) {
        as.numeric(m$pip)
    } else {
        NULL
    }
    lead <- if (!is.null(pipVec) && any(is.finite(pipVec))) {
        which.max(pipVec)
    } else {
        integer(0)
    }
    eff <- spec$eff
    tibble(
        cs = str_c(spec$methodTag, "_", eff),
        effect_id = if (!is.na(eff)) str_c("L", eff) else NA_character_,
        coverage = ctx$coverage,
        n_variants = spec$n,
        purity_min = spec$purity_min,
        purity_mean = spec$purity_mean,
        # V: the effect's prior variance, keyed by true effect index.
        V = .csEffectV(if (!is.null(ctx$fit$V)) as.numeric(ctx$fit$V) else NULL, eff),
        # cs_log10bf: the strongest MEMBER variant's logBF (max over the CS).
        cs_log10bf = .csLog10Bf(m),
        # cs_log_bf: the true per-EFFECT single-effect log Bayes factor.
        cs_log_bf = .csEffectLogBf(ctx$fit, eff),
        # cs_pip: total posterior inclusion mass captured by the CS.
        cs_pip = if (!is.null(pipVec)) {
            suppressWarnings(sum(pipVec, na.rm = TRUE))
        } else {
            NA_real_
        },
        # cs_mean_effect: mean posterior conditional effect over the CS (NA when
        # the entry has no conditional_effect column, e.g. univariate susie).
        cs_mean_effect = .csMeanEffect(m),
        lead_variant = if (length(lead) > 0L) {
            as.character(m$variant_id[lead])
        } else {
            NA_character_
        },
        lead_pip = if (length(lead) > 0L) as.numeric(pipVec[lead]) else NA_real_
    )
}

# The entry-topLoci rows for the given member variant ids, matched by
# variant_id (robust to any PIP filtering / reordering of the table). NULL when
# the table lacks a variant_id column.
# @noRd
.csMemberRows <- function(tl, ids) {
    if (is.null(tl) || !is_in("variant_id", names(tl)) || all(is.na(ids))) {
        return(NULL)
    }
    tl[match(ids, tl$variant_id), , drop = FALSE]
}

# Per-set purity at set POSITION k, first available of `cols` (NA when absent).
# @noRd
.csPurityAtK <- function(purity, k, cols) {
    if (is.null(purity)) {
        return(NA_real_)
    }
    nm <- intersect(cols, colnames(purity))
    if (length(nm) == 0L || k > nrow(purity)) {
        return(NA_real_)
    }
    as.numeric(purity[[nm[1L]]][k])
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
    if (is_in("logBF", names(m))) {
        suppressWarnings(max(m$logBF, na.rm = TRUE))
    } else {
        NA_real_
    }
}

# cs_mean_effect: mean conditional effect over the CS (NA when absent).
# @noRd
.csMeanEffect <- function(m) {
    if (is_in("conditional_effect", names(m))) {
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
#' @return A \code{tibble}: \code{cs, effect_id, coverage, n_variants,
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
        .fmrAggregateView(
            x,
            perEntry = getCredibleSetSummary,
            coverage = coverage
        )
    }
)
