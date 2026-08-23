# =============================================================================
# FineMappingRow S4 class
# -----------------------------------------------------------------------------
# One row's fine-mapping payload: the variants (with topLoci as their metadata
# columns) plus the stored fit. This is what a FineMappingResult collection is
# built FROM and what its per-row views are computed OVER.
#
# It replaces the retired FineMappingRow, but is a different object in the
# two ways that made that class expensive:
#
#   * It is NOT a dispatch target. The per-row views (.fmrRowTopLoci, .fmrRowCs,
#     ...) take it as data; it carries only field accessors, so there is no
#     parallel method surface duplicating the collection's.
#   * It holds ONE GRanges, not `variantIds` and `topLoci` in separate slots,
#     so the two cannot drift apart and no per-rebuild alignment check is
#     needed. The alignment is verified once, in the constructor, where the
#     caller still supplies them separately.
# =============================================================================

#' @include AllGenerics.R
NULL

#' @title Fine-Mapping Row
#' @description One row's worth of fine-mapping: the variants and the stored
#'   fit. Build one with \code{\link{fineMappingRow}} and pass a list of them
#'   as a collection's \code{entry} argument.
#' @slot variants A \code{GRanges} of the fit's variants whose metadata
#'   columns are the per-variant \code{topLoci} values.
#' @slot susieFit The stored fine-mapping fit, or \code{NULL}.
#' @slot cvResult Cross-validation payload, or \code{NULL}.
#' @seealso \code{\link{fineMappingRow}},
#'   \code{\linkS4class{TwasWeightsRow}}
#' @export
setClass(
    "FineMappingRow",
    representation(
        variants = "GRanges",
        susieFit = "ANY",
        cvResult = "ANY"
    ),
    prototype(susieFit = NULL, cvResult = NULL)
)

methods::setValidity("FineMappingRow", function(object) {
    errors <- character(0)
    if (!is.null(object@cvResult) && !is.list(object@cvResult)) {
        errors <- c(errors, "cvResult must be NULL or a list")
    }
    md <- mcols(object@variants, use.names = FALSE)
    if (!is.null(md) && nrow(md) != length(object@variants)) {
        errors <- c(
            errors,
            "variants' metadata columns must have one row per variant"
        )
    }
    if (length(errors) == 0L) TRUE else errors
})

#' @title Build One Fine-Mapping Row
#' @description Assemble a single row's payload for
#'   \code{\link{QtlFineMappingResult}} / \code{\link{GwasFineMappingResult}}:
#'   the variants, the stored fit and the per-variant \code{topLoci} table.
#'   Pass a list of these as the collections' \code{entry} argument.
#'
#'   The alignment between \code{variantIds} and \code{topLoci} is checked
#'   here, once, at the boundary where the two are still separate --
#'   \code{topLoci}'s columns are taken positionally, so a table listing the
#'   same variants in a different order would silently mis-assign every one.
#'   Past this point they are a single \code{GRanges} and cannot disagree.
#' @param variantIds Character vector of variant ids. Each must encode
#'   coordinates (\code{chrom:pos:ref:alt}): an id renders a range and alleles,
#'   so one without coordinates has no variant identity to store.
#' @param susieFit The stored fine-mapping fit, or \code{NULL}.
#' @param topLoci A per-variant table aligned row-for-row with
#'   \code{variantIds}.
#' @param cvResult Optional cross-validation payload.
#' @return A \code{\linkS4class{FineMappingRow}}.
#' @seealso \code{\link{twasWeightsRow}}
#' @examples
#' tl <- data.frame(
#'     variant_id = c("chr1:100:A:G", "chr1:200:C:T"),
#'     pip = c(0.9, 0.1)
#' )
#' row <- fineMappingRow(tl$variant_id, susieFit = list(), topLoci = tl)
#' QtlFineMappingResult(
#'     study = "s1", context = "c1", trait = "g1", method = "susie",
#'     entry = list(row)
#' )
#' @export
fineMappingRow <- function(variantIds, susieFit, topLoci, cvResult = NULL) {
    vids <- as.character(variantIds)
    gr <- .variantIdsToGRanges(vids, "variantIds")
    tl <- as_tibble(topLoci)
    missingCols <- setdiff(c("variant_id", "pip"), colnames(tl))
    if (nrow(tl) > 0L && length(missingCols) > 0L) {
        msg <- glue(
            "topLoci missing required columns: ",
            "{str_flatten(missingCols, ', ')}"
        )
        abort(msg)
    }
    if (nrow(tl) != length(gr)) {
        msg <- glue(
            "topLoci has {nrow(tl)} rows but {length(gr)} variants were ",
            "supplied; they must be aligned row-for-row."
        )
        abort(msg)
    }
    if (nrow(tl) > 0L && is_in("variant_id", names(tl))) {
        .fmRowCheckIdOrder(vids, as.character(tl$variant_id))
    }
    for (nm in setdiff(names(tl), .fmeIdentityCols)) {
        mcols(gr)[[nm]] <- tl[[nm]]
    }
    obj <- new(
        "FineMappingRow",
        variants = gr,
        susieFit = susieFit,
        cvResult = cvResult
    )
    validObject(obj)
    obj
}


# topLoci's mcols are taken POSITIONALLY, so a table listing the same variants
# in a different order would silently mis-assign every column.
# @noRd
.fmRowCheckIdOrder <- function(vids, tlIds) {
    if (identical(vids, tlIds)) {
        return(invisible(NULL))
    }
    bad <- which(vids != tlIds)
    msg <- glue(
        "topLoci$variant_id must be in the same order as `variantIds`; they ",
        "first differ at position {bad[[1L]]} ",
        "('{vids[bad[[1L]]]}' vs '{tlIds[bad[[1L]]]}')."
    )
    abort(msg)
}


# ---- FineMappingRow validity helpers -------------------------------------

# Full validity check: collect all contract violations (empty vector = valid).
# @noRd

# cvResult must be NULL or a list.
# @noRd

# topLoci contract checks (only meaningful when topLoci has rows).
# @noRd

# Minimal column contract: variant_id + pip. Canonical projector columns
# (marginal_*, posterior_*, ...) are pipeline-populated; skeletal entries may
# omit them (accessor projections then return NA-filled cols).
# @noRd

# topLoci row count must match variantIds length.
# @noRd

# topLoci$variant_id must equal variantIds in order.
#
# Still earns its place after the migration. The collection stores topLoci AS
# the element's inner mcols, where alignment is definitional -- but the entry
# VIEW keeps `variantIds` and `topLoci` in separate slots, and
# .fmeEntryToGRanges() builds the range from the former while taking mcols
# positionally from the latter. Without this check a misaligned entry would
# convert into silently mis-assigned mcols. It becomes deletable only if the
# entry stops carrying the two separately.
# @noRd

# Drift check: a susieFit$pip vector must match the topLoci pip column. Catches
# adjustPips() (or any mutator) updating one and forgetting the other.
# @noRd

# ---- fSuSiE credible-band / affected-region helpers --------------------

.isFsusieFit <- function(fit) {
    !is.null(fit) &&
        is.list(fit) &&
        !is.null(fit$fitted_wc2) &&
        !is.null(fit$fitted_func) &&
        !is.null(fit$outing_grid) &&
        !is.null(fit$alpha)
}

# Populate `obj$cred_band` via fsusieR's wavethresh/GenW band computation. That
# function is registered as an S3 method but NOT exported, and the exported
# affected_reg() depends on cred_band already being populated, so this internal
# call is the only path that works across all post_processing modes. Guarded so
# an upstream fsusieR change surfaces as a clear error, not a silent NULL.
# @noRd
.fsusiePopulateCredibleBand <- function(fit) {
    fn <- tryCatch(
        get("update_cal_credible_band.susiF", envir = asNamespace("fsusieR")),
        error = function(e) NULL
    )
    if (is.null(fn)) {
        # nocov start  (defensive guard against an upstream fsusieR rename; only
        # reachable if fsusieR drops this unexported S3 method)
        msg <- glue(
            "fsusieR's internal update_cal_credible_band.susiF not found; ",
            "cannot compute the fSuSiE credible band (upstream fsusieR API ",
            "changed)."
        )
        abort(msg)
        # nocov end
    }
    indxLst <- fsusieR::gen_wavelet_indx(log2(length(fit$outing_grid)))
    fn(fit, indxLst)
}

# Chromosome of the fit's variants (matches getTopLoci's `chrom`, no "chr").
# @noRd
.fsusieChrom <- function(fit) {
    vid <- names(fit$csd_X)
    if (is.null(vid) || length(vid) == 0L) {
        return(NA_character_)
    }
    tryCatch(parseVariantId(vid[1])$chrom, error = function(e) NA_character_)
}

.emptyCredibleBand <- function() {
    tibble(
        cs = character(),
        chrom = character(),
        pos = numeric(),
        effect = numeric(),
        lower = numeric(),
        upper = numeric()
    )
}

# @noRd
.fsusieCredibleBandFit <- function(fit) {
    if (!.isFsusieFit(fit)) {
        return(.emptyCredibleBand())
    }
    fit <- .fsusiePopulateCredibleBand(fit)
    grid <- as.numeric(fit$outing_grid)
    chrom <- .fsusieChrom(fit)
    parts <- map(
        seq_along(fit$cred_band),
        .fsusieBandRow,
        fit = fit,
        chrom = chrom,
        grid = grid
    )
    parts <- compact(parts)
    if (length(parts) == 0L) {
        return(.emptyCredibleBand())
    }
    bind_rows(parts)
}

# Map each fit-effect index (1..L, as affected_reg's `CS`) to pecotmr's
# canonical CS label + purity. The raw fSuSiE fit stores no purity (fsusieR
# computes it from the genotype matrix, which the fit does not retain), so the
# canonical value is `cs_95_purity` in the entry's topLoci. Matching is by the
# CS's variant
# MEMBERSHIP (not by index), which is robust to pecotmr's CS
# filtering/renumbering.
# @noRd
.fsusieCsMapFromTopLoci <- function(fit, topLoci) {
    L <- length(fit$cs)
    label <- set_names(str_c("fsusie_", seq_len(L)), seq_len(L))
    purity <- set_names(rep(NA_real_, L), seq_len(L))
    vids <- names(fit$csd_X)
    if (
        is.null(topLoci) ||
            is.null(vids) ||
            !all(is_in(
                c("cs_95", "cs_95_purity", "variant_id"),
                names(topLoci)
            ))
    ) {
        return(list(label = label, purity = purity))
    }
    retained <- unique(topLoci$cs_95[!str_detect(topLoci$cs_95, "_0$")])
    for (lab in retained) {
        inLab <- topLoci$cs_95 == lab
        kVars <- as.character(topLoci$variant_id[inLab])
        kPurit <- as.numeric(topLoci$cs_95_purity[inLab][1])
        for (l in seq_len(L)) {
            fv <- fit$cs[[l]]
            if (is.numeric(fv)) {
                fv <- vids[as.integer(fv)]
            }
            fv <- as.character(fv)
            if (length(fv) > 0L && setequal(fv, kVars)) {
                label[[as.character(l)]] <- lab
                purity[[as.character(l)]] <- kPurit
                break
            }
        }
    }
    list(label = label, purity = purity)
}

# @noRd
.fsusieAffectedRegionsFit <- function(fit, topLoci = NULL) {
    if (!.isFsusieFit(fit)) {
        return(GenomicRanges::GRanges())
    }
    fit <- .fsusiePopulateCredibleBand(fit)
    reg <- tryCatch(fsusieR::affected_reg(fit), error = function(e) NULL)
    if (is.null(reg) || nrow(reg) == 0L) {
        return(GenomicRanges::GRanges())
    }
    reg <- as_tibble(reg)
    chrom <- .fsusieChrom(fit)
    grid <- as.numeric(fit$outing_grid)
    csMap <- .fsusieCsMapFromTopLoci(fit, topLoci)
    csKey <- as.character(reg$CS)
    # Effect direction over each region (sign of the fitted effect curve), which
    # upstream affected_reg() collapses away.
    direction <- map_chr(
        seq_len(nrow(reg)),
        .fsusieRegionDirection,
        grid = grid,
        reg = reg,
        fit = fit
    )
    GenomicRanges::GRanges(
        seqnames = str_c("chr", str_remove(chrom, "^chr")),
        ranges = IRanges::IRanges(
            start = as.integer(reg$Start),
            end = as.integer(reg$End)
        ),
        cs = unname(csMap$label[csKey]),
        purity = unname(csMap$purity[csKey]),
        direction = direction
    )
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# One credible band's variant-level data.frame (effect + lower/upper), or NULL
# when band `l` has no fitted band/effect.
# @noRd
.fsusieBandRow <- function(l, fit, chrom, grid) {
    band <- fit$cred_band[[l]]
    eff <- fit$fitted_func[[l]]
    if (is.null(band) || is.null(eff)) {
        return(NULL)
    }
    tibble(
        cs = str_c("fsusie_", l),
        chrom = chrom,
        pos = grid,
        effect = as.numeric(eff),
        lower = as.numeric(band[2, ]),
        upper = as.numeric(band[1, ])
    )
}

# The effect direction ("pos"/"neg"/NA) of affected region `i` (sign of the
# fitted effect curve over the region's grid range).
# @noRd
.fsusieRegionDirection <- function(i, grid, reg, fit) {
    inRange <- which(grid >= reg$Start[i] & grid <= reg$End[i])
    eff <- fit$fitted_func[[reg$CS[i]]]
    if (length(inRange) == 0L || is.null(eff)) {
        return(NA_character_)
    }
    s <- sign(mean(as.numeric(eff)[inRange], na.rm = TRUE))
    if (is.na(s) || s == 0) {
        NA_character_
    } else if (s > 0) {
        "pos"
    } else {
        "neg"
    }
}

# One entry's affected regions, labelled with the row's study/context/trait
# mcols.
# @noRd
.fsusieEntryAffectedRegions <- function(i, x) {
    gr <- .fmrRowFsusieAffectedRegions(.fmrRowParts(x, i))
    if (length(gr) > 0L) {
        S4Vectors::mcols(gr)$study <- as.character(x$study[[i]])
        if (is_in("context", .tupleColumnNames(x))) {
            S4Vectors::mcols(gr)$context <- as.character(x$context[[i]])
        }
        if (is_in("trait", .tupleColumnNames(x))) {
            S4Vectors::mcols(gr)$trait <- as.character(x$trait[[i]])
        }
    }
    gr
}


# ---- getCredibleSetSummary helpers -------------------------------------

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
    csCol <- str_c("cs_", coverage * 100)
    purCol <- str_c(csCol, "_purity")
    if (is.null(tl) || nrow(tl) == 0L || !is_in(csCol, names(tl))) {
        return(.emptyCsSummary())
    }
    csValues <- tl[[csCol]]
    labels <- unique(csValues[
        !is.na(csValues) &
            str_length(csValues) > 0L &
            !str_detect(csValues, "_0$")
    ])
    if (length(labels) == 0L) {
        return(.emptyCsSummary())
    }
    ctx <- .csSummaryContext(tl, fit, coverage, csCol, purCol)
    map(labels, .csSummaryRow, ctx = ctx) |>
        bind_rows() |>
        mutate(
            cs_log10bf = if_else(
                is.finite(.data$cs_log10bf),
                .data$cs_log10bf,
                NA_real_
            )
        )
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
            meanPur <- set_names(as.numeric(sp[[mc[1L]]]), rownames(sp))
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
    m <- filter(ctx$tl, .data[[ctx$csCol]] == lab)
    Lidx <- suppressWarnings(as.integer(str_remove(lab, "^.*_")))
    Lname <- if (!is.na(Lidx)) str_c("L", Lidx) else NA_character_
    lead <- which.max(m$pip)
    tibble(
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
        lead_pip = as.numeric(m$pip[lead])
    )
}

# purity_min: first per-CS purity value (NA when the column is absent).
# @noRd
.csPurityMin <- function(m, purCol) {
    if (is_in(purCol, names(m))) as.numeric(m[[purCol]][1]) else NA_real_
}

# purity_mean: per-effect mean absolute correlation (NA when unavailable).
# @noRd
.csPurityMean <- function(meanPur, Lname) {
    if (!is.null(meanPur) && !is.na(Lname) && is_in(Lname, names(meanPur))) {
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
    purCols <- names(tl)[str_detect(names(tl), "^cs_[0-9.]+_purity$")]
    if (length(purCols) == 0L) {
        return(rep(TRUE, nrow(tl)))
    }
    covs <- suppressWarnings(as.numeric(str_remove(
        str_remove(purCols, "^cs_"),
        "_purity$"
    )))
    primaryPur <- purCols[which.max(covs)]
    csCol <- str_remove(primaryPur, "_purity$")
    inCs <- if (is_in(csCol, names(tl))) {
        !str_detect(tl[[csCol]], "_0$")
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
        seqnames = str_c("chr", parsed$chrom),
        ranges = IRanges::IRanges(start = parsed$pos, width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(out)
    gr
}


# Match the entry's variants to `keepVariants` by (chrom, pos, allele) so a
# chr-prefix / separator difference does not read as no-overlap. PIP / LBF are
# coding-invariant (no sign). The index is sorted to keep the entry's variant
# order so fit-matrix subsetting is unchanged.
# @noRd
.adjustPipsKeepIdx <- function(variantIds, keepVariants) {
    m <- matchVariants(variantIds, as.character(keepVariants))
    if (!length(m$idxA)) {
        msg <- glue(
            "adjustPips: intersection of entry variants with `keepVariants` ",
            "is empty."
        )
        abort(msg)
    }
    sort(m$idxA)
}

# Restrict the SuSiE fit to `keepIdx` and renormalize every single-effect
# posterior over the retained variants.
# @noRd
.adjustPipsSubsetFit <- function(fit, keepIdx, nVariants) {
    .adjustPipsRejectWeighted(fit)
    spec <- .adjustPipsAlphaSpec(fit, nVariants)
    cols <- c(keepIdx, spec$nullIdx)
    alpha <- .adjustPipsRenormAlpha(spec$alpha, cols)
    nKeep <- length(keepIdx)
    varAlpha <- alpha[, seq_len(nKeep), drop = FALSE]
    fit$alpha <- alpha
    fit$pip <- as.numeric(1 - apply(1 - varAlpha, 2L, prod))
    if (length(spec$nullIdx) > 0L) {
        fit$null_index <- nKeep + 1L
    }
    # Recorded BEFORE renormalization: after it every row sums to 1, so the
    # share that survived is only knowable at this point. This is the §3.6
    # diagnostic -- always reported, never a filter.
    fit$retained_mass <- .adjustPipsRetainedMass(spec$alpha, cols)
    fit <- .adjustPipsSubsetVariantSlots(fit, keepIdx, cols, nVariants)
    .adjustPipsRebuildAllSets(fit, varAlpha)
}

# susieInf / susieAsh carry Omega-weighted per-variant terms whose posterior
# moments are not separable by variant, so no honest subset exists without
# refitting on the retained set. Refuse rather than return a fit whose `theta`
# still spans the pre-subset variants -- coef.susie would then recycle it
# against a shorter alpha and silently emit the wrong number of weights.
# @noRd
.adjustPipsRejectWeighted <- function(fit) {
    weighted <- c("theta", "omega_weights")
    present <- weighted[map_lgl(weighted, .adjustPipsHasField, fit = fit)]
    if (length(present) == 0L) {
        return(invisible(NULL))
    }
    msg <- glue(
        "adjustPips: fit carries Omega-weighted per-variant term(s) ",
        "({str_flatten(present, ', ')}); a susieInf / susieAsh posterior ",
        "cannot be honestly restricted to a variant subset. Refit on the ",
        "retained variants instead."
    )
    abort(msg)
}

# TRUE when `fit` carries a non-NULL field named `nm`.
# @noRd
.adjustPipsHasField <- function(nm, fit) {
    !is.null(fit[[nm]])
}

# Resolve the alpha matrix and the position of a null_weight column. susieR
# appends the null as alpha column p + 1 while `pip` stays length p, so the
# column count -- not the pip length -- is what reveals it.
# @noRd
.adjustPipsAlphaSpec <- function(fit, nVariants) {
    if (is.null(fit[["alpha"]])) {
        msg <- glue(
            "adjustPips: entry's susieFit has no `alpha` matrix; ",
            "renormalization restricts the stored per-effect alpha."
        )
        abort(msg)
    }
    alpha <- as.matrix(fit[["alpha"]])
    nCol <- ncol(alpha)
    if (nCol == nVariants) {
        return(list(alpha = alpha, nullIdx = integer(0)))
    }
    if (nCol == nVariants + 1L) {
        return(list(alpha = alpha, nullIdx = .adjustPipsNullIdx(fit, nCol)))
    }
    msg <- glue(
        "adjustPips: susieFit$alpha has {nCol} columns but the entry has ",
        "{nVariants} variants (a null_weight fit would have ",
        "{nVariants + 1L})."
    )
    abort(msg)
}

# susieR appends the null column last, so a non-trailing null_index means a
# layout this code cannot slice safely.
# @noRd
.adjustPipsNullIdx <- function(fit, nCol) {
    nullIdx <- fit[["null_index"]]
    if (is.null(nullIdx) || length(nullIdx) != 1L || !is.finite(nullIdx)) {
        return(nCol)
    }
    if (as.integer(nullIdx) != nCol) {
        msg <- glue(
            "adjustPips: null_index is {nullIdx} but the null column must be ",
            "the last of {nCol} alpha columns."
        )
        abort(msg)
    }
    nCol
}

# Renormalize each single-effect posterior over `cols`. alpha already carries
# the fit's prior (alpha is proportional to pi * BF), so restricting it and
# renormalizing is exact for any prior weights -- the fit's `pi` is not needed,
# and rebuilding alpha from lbf_variable would silently substitute a uniform
# prior. Working in log space keeps an effect whose retained mass has
# underflowed from collapsing to NaN. A null_weight column stays in `cols` so
# the "no signal here" mass survives.
# @noRd
.adjustPipsRenormAlpha <- function(alpha, cols) {
    logSub <- log(alpha[, cols, drop = FALSE])
    rowMax <- apply(logSub, 1L, max)
    dead <- !is.finite(rowMax)
    if (any(dead)) {
        msg <- glue(
            "adjustPips: effect(s) {str_flatten(which(dead), ', ')} carry no ",
            "posterior mass on the retained variants, so their single-effect ",
            "posterior cannot be renormalized."
        )
        abort(msg)
    }
    w <- exp(sweep(logSub, 1L, rowMax, `-`))
    w / rowSums(w)
}

# The share of each single effect's posterior that the retained variants carry.
#
# Read off the PRE-renormalization alpha, whose rows sum to 1 over the full
# variant set, so the row-sum over the retained columns is exactly the surviving
# share. A null_weight column counts as retained: its mass is the "no signal
# here" hypothesis, which the subset does not remove.
# @noRd
.adjustPipsRetainedMass <- function(alpha, cols) {
    as.numeric(rowSums(alpha[, cols, drop = FALSE]))
}

# Restrict every per-variant fit slot. Slots spanning alpha's columns take
# `cols` (which carries any null_weight column); variant-only slots take
# `keepIdx`. A slot whose variant axis matches neither width is a contract
# violation, not something to pass through silently.
#
# Fields are READ with `[[`, never `$`: `$` on a list falls back to prefix
# matching when the exact name is absent, so `fit$pi` silently returns `pip`,
# `fit$mu` returns `mu2`, and `fit$sets` returns `sets_secondary`.
# @noRd
.adjustPipsSubsetVariantSlots <- function(fit, keepIdx, cols, nVariants) {
    p <- list(keepIdx = keepIdx, cols = cols, nVariants = nVariants)
    fit$lbf_variable <- .adjustPipsCols(
        fit[["lbf_variable"]],
        p,
        "lbf_variable"
    )
    fit$mu <- .adjustPipsCols(fit[["mu"]], p, "mu")
    fit$mu2 <- .adjustPipsCols(fit[["mu2"]], p, "mu2")
    fit$mu2_diag <- .adjustPipsCols(fit[["mu2_diag"]], p, "mu2_diag")
    fit$clfsr <- .adjustPipsCols(fit[["clfsr"]], p, "clfsr")
    fit$coef <- .adjustPipsRows(fit[["coef"]], p, "coef")
    fit$X_column_scale_factors <- .adjustPipsVec(
        fit[["X_column_scale_factors"]],
        p,
        "X_column_scale_factors"
    )
    fit$XtXr <- .adjustPipsVec(fit[["XtXr"]], p, "XtXr")
    fit$pi <- .adjustPipsVec(fit[["pi"]], p, "pi")
    fit
}

# Index set for a slot measured along its variant axis: `nVariants` selects the
# retained variants, `nVariants + 1` those plus the null_weight column.
# @noRd
.adjustPipsIdxFor <- function(n, p, what) {
    if (n == p$nVariants) {
        return(p$keepIdx)
    }
    if (n == p$nVariants + 1L && length(p$cols) > length(p$keepIdx)) {
        return(p$cols)
    }
    msg <- glue(
        "adjustPips: fit slot `{what}` spans {n} variants but the entry has ",
        "{p$nVariants}."
    )
    abort(msg)
}

# Column-subset a (possibly 3-D) per-variant fit matrix; NULL passes through
# and a dimensionless slot is treated as a per-variant vector.
# @noRd
.adjustPipsCols <- function(m, p, what) {
    if (is.null(m)) {
        return(NULL)
    }
    d <- dim(m)
    if (is.null(d)) {
        return(.adjustPipsVec(m, p, what))
    }
    idx <- .adjustPipsIdxFor(d[[2L]], p, what)
    if (length(d) == 3L) {
        m[, idx, , drop = FALSE]
    } else {
        m[, idx, drop = FALSE]
    }
}

# Row-subset a variants-by-something fit matrix (mvsusie / fsusie `coef`).
# @noRd
.adjustPipsRows <- function(m, p, what) {
    if (is.null(m)) {
        return(NULL)
    }
    m <- as.matrix(m)
    m[.adjustPipsIdxFor(nrow(m), p, what), , drop = FALSE]
}

# Subset a per-variant vector slot; NULL passes through.
# @noRd
.adjustPipsVec <- function(v, p, what) {
    if (is.null(v)) {
        return(NULL)
    }
    v[.adjustPipsIdxFor(length(v), p, what)]
}

# Rebuild the credible sets from the renormalized alpha. Renormalization
# changes CS composition, so remapping the stored variant indices would point
# at the wrong variants. Purity is LD-dependent and cannot be recomputed
# without a reference panel, so it is dropped rather than carried forward
# stale.
# @noRd
.adjustPipsRebuildAllSets <- function(fit, varAlpha) {
    if (!is.null(fit[["sets"]])) {
        fit$sets <- .adjustPipsCsAtCoverage(
            varAlpha,
            fit[["V"]],
            .adjustPipsCoverage(fit[["sets"]])
        )
    }
    if (!is.null(fit[["sets_secondary"]])) {
        fit$sets_secondary <- map(
            fit[["sets_secondary"]],
            .adjustPipsRebuildSecondary,
            varAlpha = varAlpha,
            V = fit[["V"]]
        )
    }
    fit
}

# Rebuild one secondary-coverage credible-set table at its own coverage.
# @noRd
.adjustPipsRebuildSecondary <- function(entry, varAlpha, V) {
    if (is.null(entry$sets)) {
        return(entry)
    }
    entry$sets <- .adjustPipsCsAtCoverage(
        varAlpha,
        V,
        .adjustPipsCoverage(entry$sets)
    )
    entry
}

# Requested coverage recorded on a susie_get_cs() result (0.95 when absent).
# @noRd
.adjustPipsCoverage <- function(sets) {
    coverage <- sets$requested_coverage
    if (
        is.null(coverage) ||
            length(coverage) != 1L ||
            !is.finite(coverage)
    ) {
        return(0.95)
    }
    as.numeric(coverage)
}

# Smallest-prefix credible sets from a renormalized alpha, in susie_get_cs's
# shape. Effects at or below `priorTol` prior variance are skipped, matching
# susie_get_cs's own V filter; purity is left NULL (no LD available here).
# @noRd
.adjustPipsCsAtCoverage <- function(varAlpha, V, coverage, priorTol = 1e-9) {
    idx <- if (is.null(V)) seq_len(nrow(varAlpha)) else which(V > priorTol)
    cs <- set_names(
        map(idx, .adjustPipsOneCs, varAlpha = varAlpha, coverage = coverage),
        str_c("L", idx)
    )
    list(
        cs = cs,
        cs_index = as.integer(idx),
        coverage = map_dbl(
            seq_along(idx),
            .adjustPipsCsMass,
            cs = cs,
            varAlpha = varAlpha,
            idx = idx
        ),
        requested_coverage = coverage,
        purity = NULL
    )
}

# Variant indices of one effect's smallest prefix reaching `coverage`.
# @noRd
.adjustPipsOneCs <- function(l, varAlpha, coverage) {
    a <- varAlpha[l, ]
    ord <- order(a, decreasing = TRUE)
    hit <- which(cumsum(a[ord]) >= coverage)
    n <- if (length(hit) == 0L) length(ord) else hit[[1L]]
    sort(ord[seq_len(n)])
}

# Realized posterior mass of the i-th recomputed credible set.
# @noRd
.adjustPipsCsMass <- function(i, cs, varAlpha, idx) {
    sum(varAlpha[idx[[i]], cs[[i]]])
}

# Alpha restricted to the variant columns (any null_weight column dropped).
# @noRd
.adjustPipsVariantAlpha <- function(fit) {
    alpha <- as.matrix(fit[["alpha"]])
    nullIdx <- fit[["null_index"]]
    hasNull <- !is.null(nullIdx) &&
        length(nullIdx) == 1L &&
        is.finite(nullIdx) &&
        nullIdx >= 1L &&
        nullIdx <= ncol(alpha)
    if (hasNull) {
        alpha[, -as.integer(nullIdx), drop = FALSE]
    } else {
        alpha
    }
}

# Rebuild topLoci from the subset fit + existing (per-variant) marginal columns.
# @noRd
.adjustPipsRebuildTopLoci <- function(topLoci, fit, common) {
    if (nrow(topLoci) == 0L) {
        return(topLoci)
    }
    newTopLoci <- filter(topLoci, is_in(.data$variant_id, common))
    newTopLoci$pip <- as.numeric(fit[["pip"]])
    .adjustPipsRelabelCs(.adjustPipsPosterior(newTopLoci, fit), fit)
}

# Rewrite the cs_<C> membership labels from the recomputed credible sets, and
# blank cs_<C>_purity: purity is LD-dependent, so the stored value describes a
# credible set that no longer exists.
# @noRd
.adjustPipsRelabelCs <- function(tl, fit) {
    csCols <- names(tl)[str_detect(names(tl), "^cs_[0-9]+$")]
    if (length(csCols) == 0L || is.null(fit[["alpha"]])) {
        return(tl)
    }
    varAlpha <- .adjustPipsVariantAlpha(fit)
    label <- .adjustPipsMethodLabel(tl)
    for (csCol in csCols) {
        tl[[csCol]] <- .adjustPipsCsLabels(
            varAlpha,
            fit[["V"]],
            csCol,
            label,
            nrow(tl)
        )
        purCol <- str_c(csCol, "_purity")
        if (is_in(purCol, names(tl))) {
            tl[[purCol]] <- NA_real_
        }
    }
    tl
}

# Membership labels at the coverage encoded in the column name: "<method>_<k>"
# for members of effect k, "<method>_0" for variants in no credible set.
# @noRd
.adjustPipsCsLabels <- function(varAlpha, V, csCol, label, n) {
    coverage <- as.numeric(str_remove(csCol, "^cs_")) / 100
    sets <- .adjustPipsCsAtCoverage(varAlpha, V, coverage)
    out <- rep(str_c(label, "_0"), n)
    for (k in seq_along(sets$cs)) {
        out[sets$cs[[k]]] <- str_c(label, "_", sets$cs_index[[k]])
    }
    out
}

# The method label topLoci uses to prefix credible-set names.
# @noRd
.adjustPipsMethodLabel <- function(tl) {
    if (!is_in("method", names(tl))) {
        return("cs")
    }
    labels <- unique(as.character(tl$method))
    labels <- labels[!is.na(labels)]
    if (length(labels) == 0L) {
        return("cs")
    }
    labels[[1L]]
}

# Recompute posterior_mean / posterior_sd from the fit when alpha + mu/mu2 are
# matrix-shaped and conformable; otherwise leave the existing values in place.
# @noRd
.adjustPipsPosterior <- function(newTopLoci, fit) {
    alphaMat <- if (!is.null(fit[["alpha"]])) {
        .adjustPipsVariantAlpha(fit)
    } else {
        NULL
    }
    muMat <- if (!is.null(fit[["mu"]])) as.matrix(fit[["mu"]]) else NULL
    if (
        is.null(alphaMat) ||
            is.null(muMat) ||
            !all(dim(alphaMat) == dim(muMat))
    ) {
        return(newTopLoci)
    }
    newTopLoci$posterior_mean <- as.numeric(colSums(alphaMat * muMat))
    mu2Mat <- if (!is.null(fit[["mu2"]])) as.matrix(fit[["mu2"]]) else NULL
    if (!is.null(mu2Mat) && all(dim(alphaMat) == dim(mu2Mat))) {
        newTopLoci$posterior_sd <- as.numeric(sqrt(pmax(
            colSums(alphaMat * mu2Mat) - newTopLoci$posterior_mean^2,
            0
        )))
    }
    newTopLoci
}


# =============================================================================
# Internal column projectors used by accessors
# =============================================================================

# Read a column from the canonical wide topLoci, returning NAs of the
# given type when the column is absent. Lets accessor projectors tolerate
# skeletal entries that lack optional columns.
# @noRd
.tlCol <- function(tl, name, type = c("character", "integer", "numeric")) {
    type <- arg_match(type)
    if (is_in(name, colnames(tl))) {
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
# + pipeline-added annotations. Renames `posterior_mean`/`posterior_sd` to
# `beta`/`se`.
# Exports effect-allele frequency as `af` (never MAF). Missing optional
# columns are NA-filled.
# @noRd
.projectPosteriorView <- function(tl) {
    if (nrow(tl) == 0L) {
        return(.projectPosteriorEmpty())
    }
    out <- tibble(
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
        logBF = .tlCol(tl, "logBF", "numeric")
    )
    for (cc in .projectPosteriorExtraCols(tl)) {
        out[[cc]] <- tl[[cc]]
    }
    out
}

# Empty posterior-view frame (canonical columns, zero rows).
# @noRd
.projectPosteriorEmpty <- function() {
    tibble(
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
        logBF = numeric(0)
    )
}

# Optional pass-through columns for the posterior view, in canonical order: the
# dynamic CS columns (cs_<C> membership / cs_<C>_purity), the per-CS variant-
# level fullFit columns (within_cs_pip[_<lab>] / cs_logbf_<lab> /
# cs_effect[_var]_<lab>), the per-condition posterior quantities
# (conditional_effect, lfsr), and pipeline-added annotations -- whichever are
# present.
# @noRd
.projectPosteriorExtraCols <- function(tl) {
    csCols <- colnames(tl)[str_detect(colnames(tl), "^cs_[0-9.]+(_purity)?$")]
    fullFitCols <- colnames(tl)[str_detect(
        colnames(tl),
        "^(within_cs_pip|cs_logbf_|cs_effect_)"
    )]
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
        return(tibble(
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
            p = numeric(0)
        ))
    }
    tibble(
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
        p = .tlCol(tl, "marginal_p", "numeric")
    )
}

# The non-null credible-set labels in one cs_<coverage> column (dropping the
# "_0" not-captured sentinel).
# @noRd
.fmeNonNullCsLabels <- function(cc, tl) {
    v <- tl[[cc]]
    v[!str_detect(v, "_0$")]
}

# =============================================================================
# FineMappingRow <-> collection element
# -----------------------------------------------------------------------------
# The entry is a DERIVED VIEW of one collection element, not stored state:
#
#   element GRanges   the fit's variants (coordinates + A1 / A2)
#   inner mcols       topLoci, minus the identity columns the range now carries
#   outer mcols       susieFit, cvResult
#
# `variant_id` / `chrom` / `pos` are dropped on the way in and re-rendered on
# the way out, so topLoci can no longer drift out of order relative to the
# variants -- the alignment is definitional rather than checked.
# =============================================================================

# Identity columns that the range and its alleles already carry.
# @noRd
.fmeIdentityCols <- c("variant_id", "chrom", "pos", "A1", "A2")

# Every `entry` must be a FineMappingRow. Checked before the payload is
# unpacked, so a wrong type names itself rather than failing as a missing
# method on whatever was passed instead.
# @noRd

# One entry -> its element.
# @noRd

# The entry view for a collection row (or rows). A tuple split across
# chromosomes owns several elements; they are stitched back into the single
# entry callers expect, the same contract getSumStats() keeps.
# @noRd

# One element (plus its outer-mcols payload) -> the entry view.
# @noRd

# Rebuild the canonical topLoci table from an element: the identity columns are
# re-rendered from the range and alleles, the rest come straight from mcols.
# @noRd
.fmeTopLociFromElement <- function(gr) {
    if (length(gr) == 0L) {
        return(as_tibble(list(variant_id = character(0), pip = numeric(0))))
    }
    mc <- as.list(mcols(gr))
    identity <- list(
        variant_id = .grVariantIds(gr),
        chrom = as.character(seqnames(gr)),
        pos = as.integer(start(gr))
    )
    as_tibble(c(identity, mc), .name_repair = "minimal")
}

#' @rdname adjustPips
#' @noRd
setMethod("adjustPips", "FineMappingRow", function(x, keepVariants, ...) {
    vids <- getVariantIds(x)
    keepIdx <- .adjustPipsKeepIdx(vids, keepVariants)
    common <- vids[keepIdx]
    fit <- .adjustPipsSubsetFit(getSusieFit(x), keepIdx, length(vids))
    newTopLoci <- .adjustPipsRebuildTopLoci(
        .fmeTopLociFromElement(rowVariants(x)),
        fit,
        common
    )
    # cvResult is sample-indexed (partition + held-out predictions), not
    # variant-indexed, so it carries through a variant-subsetting unchanged.
    fineMappingRow(
        variantIds = common,
        susieFit = fit,
        topLoci = newTopLoci,
        cvResult = getCvResult(x)
    )
})

# ---- field accessors --------------------------------------------------------
# The only methods this class carries. `@` is confined to these bodies, which
# is what keeps the slot-access rule satisfied everywhere else.

#' @rdname getVariantIds
#' @export
setMethod("getVariantIds", "FineMappingRow", function(x, ...) {
    .grVariantIds(x@variants)
})

#' @rdname getSusieFit
#' @export
setMethod("getSusieFit", "FineMappingRow", function(x, ...) x@susieFit)

#' @rdname getCvResult
#' @export
setMethod("getCvResult", "FineMappingRow", function(x, ...) x@cvResult)

#' @rdname show-methods
#' @export
setMethod("show", "FineMappingRow", function(object) {
    tl <- .fmeTopLociFromElement(object@variants)
    nCs <- if (nrow(tl) > 0L) {
        csCols <- names(tl)[str_detect(names(tl), "^cs_[0-9]+$")]
        if (length(csCols) > 0L) {
            length(unique(unlist(map(csCols, .fmeNonNullCsLabels, tl = tl))))
        } else {
            0L
        }
    } else {
        0L
    }
    cat(glue(
        "FineMappingRow: {length(object@variants)} variants, ",
        "{nCs} credible sets\n",
        .trim = FALSE
    ))
    invisible(NULL)
})

# The element itself. Used by the collection constructors, which store it
# directly as the row's element.
# @noRd
setGeneric("rowVariants", function(x, ...) standardGeneric("rowVariants"))

# @noRd
setMethod("rowVariants", "FineMappingRow", function(x, ...) x@variants)
