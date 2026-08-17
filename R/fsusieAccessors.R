#' @include AllGenerics.R AllClasses.R FineMappingEntry.R tupleSelectors.R
NULL

# Detect an untrimmed fSuSiE (susiF) fit carrying the wavelet slots that the
# credible band needs. A trimmed fit drops these, so we degrade to empty.
# @noRd
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

#' fSuSiE credible band (fitted effect + uncertainty band)
#'
#' For a functional-SuSiE (\code{fsusieR::susiF}) fine-mapping entry, returns
#' the fitted effect curve and its credible band over the functional grid, one
#' row per (credible set, grid position). Wraps the upstream fSuSiE band
#' computation (fsusieR's internal \code{update_cal_credible_band.susiF} +
#' \code{get_fitted_effect}). Requires an UNtrimmed fit (the wavelet slots are
#' dropped by trimming); degrades to zero rows for non-fSuSiE or trimmed fits.
#'
#' @param x A \code{FineMappingEntry} or \code{FineMappingResultBase}.
#' @param ... Ignored.
#' @return A long \code{tibble}: \code{cs, chrom, pos, effect, lower, upper}
#'   (the collection method additionally carries the entry identity columns).
#' @seealso \code{\link{fsusieAffectedRegions}}
#' @examples
#' data(fsusieFineMappingExample)
#' fsusieCredibleBand(fsusieFineMappingExample)
#' @export
setGeneric("fsusieCredibleBand", function(x, ...) {
    standardGeneric("fsusieCredibleBand")
})

#' @rdname fsusieCredibleBand
#' @export
setMethod("fsusieCredibleBand", "FineMappingEntry", function(x, ...) {
    .fsusieCredibleBandFit(getSusieFit(x))
})

#' @rdname fsusieCredibleBand
#' @export
setMethod("fsusieCredibleBand", "FineMappingResultBase", function(x, ...) {
    .fmrAggregateView(x, perEntry = fsusieCredibleBand)
})

#' fSuSiE affected genomic regions
#'
#' The sub-intervals of the functional grid where a credible set's band excludes
#' zero (the genomic footprint of each fSuSiE effect), from the upstream
#' \code{fsusieR::affected_reg}. Returns a \code{GRanges} with the credible-set
#' label, its purity, and the effect \code{direction}
#' (\code{"pos"}/\code{"neg"}; upstream discards it). Requires an UNtrimmed fit;
#' empty otherwise.
#'
#' @param x A \code{FineMappingEntry} or \code{FineMappingResultBase}.
#' @param ... Ignored.
#' @return A \code{GRanges} of affected intervals with mcols \code{cs},
#'   \code{purity}, \code{direction} (plus entry identity for the collection).
#' @seealso \code{\link{fsusieCredibleBand}}
#' @examples
#' data(fsusieFineMappingExample)
#' fsusieAffectedRegions(fsusieFineMappingExample)
#' @export
setGeneric("fsusieAffectedRegions", function(x, ...) {
    standardGeneric("fsusieAffectedRegions")
})

#' @rdname fsusieAffectedRegions
#' @export
setMethod("fsusieAffectedRegions", "FineMappingEntry", function(x, ...) {
    .fsusieAffectedRegionsFit(getSusieFit(x), topLoci = x@topLoci)
})

#' @rdname fsusieAffectedRegions
#' @export
setMethod("fsusieAffectedRegions", "FineMappingResultBase", function(x, ...) {
    grs <- map(seq_len(nrow(x)), .fsusieEntryAffectedRegions, x = x)
    grs <- grs[lengths(grs) > 0L]
    if (length(grs) == 0L) {
        return(GenomicRanges::GRanges())
    }
    exec(c, !!!grs)
})

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

# One entry's affected regions, stamped with the row's study/context/trait
# mcols.
# @noRd
.fsusieEntryAffectedRegions <- function(i, x) {
    gr <- fsusieAffectedRegions(x$entry[[i]])
    if (length(gr) > 0L) {
        S4Vectors::mcols(gr)$study <- as.character(x$study[[i]])
        if (is_in("context", names(x))) {
            S4Vectors::mcols(gr)$context <- as.character(x$context[[i]])
        }
        if (is_in("trait", names(x))) {
            S4Vectors::mcols(gr)$trait <- as.character(x$trait[[i]])
        }
    }
    gr
}
