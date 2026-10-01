# =============================================================================
# Hierarchical multiple-testing correction for cis-QTL association statistics.
#
# QtlSumStats-native port of the (retired) inst/code/tensorqtl_postprocessor.R
# engine. The input is a per-gene QtlSumStats: one ROW per gene/trait carrying
# the regional summary as row columns (n_variants [+ n_variants_filtered]; for
# permutation p_beta / q_beta / beta_shape1 / beta_shape2), and each ROW's
# ENTRY a per-variant GRanges with `P` (p-value) + optional `af` /
# tss_distance / tes_distance / qvalue mcols. It returns the SAME object
# enriched with the corrected-statistic columns. No file I/O, no globbing (that
# stays in the wrapper); no method key and no new class (the result is enriched
# summary statistics). Significance is never stored -- see getSignificantQtls().
#
# ALL multiple-testing math uses established package functions:
#   Bonferroni -> stats::p.adjust(method="bonferroni", n = <per-gene test
#   count>)
#   BH-FDR     -> stats::p.adjust(method="fdr")
#   Storey q   -> qvalue::qvalue  (NO hand-rolled BH-as-qvalue fallback)
#   perm thresh-> stats::qbeta    (empirical bracketing feeds the beta quantile)
# =============================================================================

#' @title Options for Storey q-value Estimation
#' @description Build a record of extra arguments for
#'   \code{qvalue::qvalue()}, the Storey q-value engine behind
#'   \code{\link{qtlAssociationPostprocess}}.
#' @param ... Arguments for \code{qvalue::qvalue()}: \code{pi0}, the
#'   \code{pi0.method} / \code{lambda} pair that sets how the null
#'   proportion is estimated, \code{pfdr}, \code{fdr.level} and
#'   \code{lfdr.out}. \code{qvalue()}'s signature ends in \code{...}, so
#'   names cannot be checked here and a misspelling is passed through.
#'   \code{p} is the p-value vector pecotmr assembles and is refused.
#'
#'   pecotmr retries \code{qvalue()} on its two documented degenerate cases
#'   ("missing or infinite" with \code{lambda = 0}, "pi0 <= 0" with a
#'   bootstrap \code{pi0}); setting \code{lambda} or \code{pi0.method}
#'   here replaces the first attempt, and the retries still override them
#'   when that attempt fails.
#' @return A \code{MethodConfig} record for
#'   \code{qtlAssociationPostprocess(qvalueArgs =)}.
#' @seealso \code{\link{qtlAssociationPostprocess}}
#' @examples
#' qvalueConfig(pi0.method = "bootstrap")
#' @export
qvalueConfig <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(p = "the p-value vector pecotmr assembles"),
        "qvalueConfig"
    )
    .newMethodConfig(
        "qvalue::qvalue",
        defaults = list(),
        extra = extra,
        label = "qvalueConfig",
        engine = "qvalue"
    )
}

# Storey q-values via the Bioconductor qvalue package, with the qvalue-NATIVE
# edge-case retries only (never a hand-rolled substitute): lambda=0 for
# missing/infinite handling, bootstrap pi0 for the "pi0 <= 0" degenerate case.
# Returns a numeric vector aligned to `p`.
#' @importFrom rlang try_fetch
.qapSafeQvalue <- function(p, qvalueArgs = list()) {
    if (!requireNamespace("qvalue", quietly = TRUE)) {
        # Optional-package guard; qvalue is Suggests-only.
        msg <- glue(
            "qtlAssociationPostprocess: the 'qvalue' package is required for ",
            "Storey q-values."
        )
        abort(msg)
    }
    extra <- as.list(qvalueArgs)
    try_fetch(
        exec(qvalue::qvalue, p, !!!extra)$qvalues,
        error = function(cnd) {
            if (str_detect(conditionMessage(cnd), "missing or infinite")) {
                exec(
                    qvalue::qvalue,
                    p,
                    !!!list_assign(extra, lambda = 0)
                )$qvalues
            } else if (str_detect(conditionMessage(cnd), "pi0 <= 0")) {
                maxP <- max(p, na.rm = TRUE)
                lambdaSeq <- seq(0, min(0.9, maxP * 0.95), length.out = 10)
                exec(
                    qvalue::qvalue,
                    p,
                    !!!list_assign(
                        extra,
                        lambda = lambdaSeq,
                        pi0.method = "bootstrap"
                    )
                )$qvalues
            } else {
                # The cause is chained via `parent`, so it is no longer
                # interpolated into the message.
                msg <- glue(
                    "qtlAssociationPostprocess: qvalue::qvalue failed. ",
                    "Not substituting a hand-rolled q-value."
                )
                abort(msg, parent = cnd)
            }
        }
    )
}

# Event-level global adjustment of a per-gene p-value vector: Benjamini-Hochberg
# FDR (stats::p.adjust) and Storey q (qvalue). Returns list(fdr=, q=).
.qapGlobalAdjust <- function(eventP, qvalueArgs = list()) {
    list(
        fdr = stats::p.adjust(eventP, method = "fdr"),
        q = .qapSafeQvalue(eventP, qvalueArgs)
    )
}

# Per-gene permutation nominal p-value threshold (FastQTL/TensorQTL): find the
# empirical p_beta cutoff that separates q_beta-significant from non-significant
# genes, then map it through each gene's beta distribution with stats::qbeta.
# The lb/ub bracketing is the FastQTL algorithm (no package equivalent); the
# distribution math is stats::qbeta. Returns a per-gene numeric vector.
.qapPermutationNominalThreshold <- function(
    pBeta,
    qBeta,
    shape1,
    shape2,
    fdrThreshold
) {
    out <- rep(NA_real_, length(pBeta))
    ok <- !is.na(pBeta) & !is.na(qBeta)
    lb <- sort(pBeta[ok & qBeta <= fdrThreshold]) # passing p_beta values
    ub <- sort(pBeta[ok & qBeta > fdrThreshold]) # failing p_beta values
    if (length(lb) == 0L) {
        return(out)
    } # no significant events
    lbVal <- lb[length(lb)] # max passing p_beta
    thr <- if (length(ub) > 0L) (lbVal + ub[1L]) / 2 else lbVal
    stats::qbeta(thr, shape1, shape2) # vectorised over per-gene shapes
}

# Per-variant MAF + cis-window keep mask (the FILTERED-flavour restriction) from
# an entry's af / tss_distance / tes_distance mcols. Shared by the correction
# and the significance derivation so the filtered set is defined once.
.qapFilterKeep <- function(mc, mafCutoff, cisWindow, afCol, nVar) {
    byMaf <- if (mafCutoff > 0 && !is.null(mc[[afCol]])) {
        pmin(mc[[afCol]], 1 - mc[[afCol]]) > mafCutoff
    } else {
        rep(TRUE, nVar)
    }
    inCis <- cisWindow > 0 &&
        !is.null(mc$tss_distance) &&
        !is.null(mc$tes_distance)
    if (!inCis) {
        return(byMaf)
    }
    byMaf & (mc$tss_distance >= -cisWindow & mc$tes_distance <= cisWindow)
}

# Per-gene logical masks of significant variants under a correction method (the
# derived significance, never stored). `method` is one of permutation,
# bonferroni_original, bonferroni_filtered, qvalue; `threshold` defaults to the
# fdrThreshold stashed by qtlAssociationPostprocess. Uses the SAME package math
# as the correction (p.adjust for the Bonferroni-adjusted per-variant p).
.qapSignificanceMask <- function(x, method, threshold = NULL) {
    recipe <- getQcInfo(x)$associationPostprocess
    if (is.null(recipe)) {
        msg <- glue(
            "This QtlSumStats was not produced by ",
            "qtlAssociationPostprocess(); ",
            "no significance recipe to derive from."
        )
        abort(msg)
    }
    if (is.null(threshold)) {
        threshold <- recipe$fdrThreshold
    }
    pcol <- recipe$pvalueCol
    masks <- map(seq_len(nrow(x)), .qapEmptyMask, x = x, pcol = pcol)
    if (method == "permutation") {
        return(.qapMaskPermutation(x, masks, pcol))
    }
    if (is_in(method, c("bonferroni_original", "bonferroni_filtered"))) {
        return(.qapMaskBonferroni(x, masks, method, threshold, recipe))
    }
    if (method == "qvalue") {
        return(.qapMaskQvalue(x, masks, threshold))
    }
    masks
}

# Permutation significance: per-entry mask of nominal p below the gene's
# permutation threshold.
# @noRd
.qapMaskPermutation <- function(x, masks, pcol) {
    if (is.null(x$p_nominal_threshold)) {
        msg <- glue(
            "getSignificantQtls: no p_nominal_threshold (run the ",
            "permutation method)."
        )
        abort(msg)
    }
    map(
        seq_len(nrow(x)),
        .qapPermutationMaskAt,
        x = x,
        masks = masks,
        thr = as.numeric(x$p_nominal_threshold),
        pcol = pcol
    )
}

# Entry `i`'s permutation mask, or the incoming mask when the gene has no
# threshold to apply.
# @noRd
.qapPermutationMaskAt <- function(i, x, masks, thr, pcol) {
    if (is.na(thr[i])) {
        return(masks[[i]])
    }
    pv <- S4Vectors::mcols(x[[i]])[[pcol]]
    !is.na(pv) & pv < thr[i]
}

# Bonferroni significance (original / filtered flavour): a global variant-level
# p threshold derived from the FDR-significant genes, applied per entry.
# @noRd
.qapMaskBonferroni <- function(x, masks, method, threshold, recipe) {
    flav <- str_remove(method, "bonferroni_")
    fdrCol <- str_c("fdr_bonferroni_min_", flav)
    pMinCol <- str_c("p_bonferroni_min_", flav)
    nCol <- if (flav == "filtered") "n_variants_filtered" else "n_variants"
    if (is.null(.tupleColumn(x, fdrCol))) {
        msg <- glue(
            "getSignificantQtls: '{method}' columns absent; recompute with ",
            "the matching flavour."
        )
        abort(msg)
    }
    sig <- replace_na(as.numeric(.tupleColumn(x, fdrCol)) < threshold, FALSE)
    if (!any(sig)) {
        return(masks)
    }
    # global scalar
    varThr <- max(as.numeric(.tupleColumn(x, pMinCol))[sig], na.rm = TRUE)
    map(
        seq_len(nrow(x)),
        .qapBonferroniMaskAt,
        x = x,
        recipe = recipe,
        flav = flav,
        nVar = as.numeric(.tupleColumn(x, nCol)),
        varThr = varThr
    )
}

# @noRd
.qapBonferroniMaskAt <- function(i, x, recipe, flav, nVar, varThr) {
    .qapBonferroniRowMask(x[[i]], recipe, flav, nVar[i], varThr)
}

# One entry's Bonferroni keep-mask: Bonferroni-adjusted p <= the global
# threshold, intersected with the filtered-flavour variant filter.
# @noRd
.qapBonferroniRowMask <- function(entry, recipe, flav, nVarI, varThr) {
    mc <- S4Vectors::mcols(entry)
    pv <- mc[[recipe$pvalueCol]]
    bySignificance <- pmin(1, pv * nVarI) <= varThr
    keep <- if (flav != "filtered") {
        bySignificance
    } else {
        bySignificance &
            .qapFilterKeep(
                mc,
                recipe$mafCutoff,
                recipe$cisWindow,
                recipe$afCol,
                length(pv)
            )
    }
    keep & !is.na(keep)
}

# Q-value significance: per-entry mask of variant q-values below the threshold
# for the FDR-significant genes.
# @noRd
# Entry `i`'s q-value mask, or the incoming mask when the gene is not
# FDR-significant or carries no per-variant q-values.
# @noRd
.qapQvalueMaskAt <- function(i, x, masks, sig, threshold) {
    mc <- S4Vectors::mcols(x[[i]])
    if (!isTRUE(sig[i]) || is.null(mc$qvalue)) {
        return(masks[[i]])
    }
    !is.na(mc$qvalue) & mc$qvalue < threshold
}

.qapMaskQvalue <- function(x, masks, threshold) {
    qCol <- if (!is.null(x$q_beta)) "q_beta" else "q_bonferroni_min_original"
    if (is.null(.tupleColumn(x, qCol))) {
        msg <- glue(
            "getSignificantQtls: no event q-value column (q_beta / ",
            "q_bonferroni_min_original)."
        )
        abort(msg)
    }
    map(
        seq_len(nrow(x)),
        .qapQvalueMaskAt,
        x = x,
        masks = masks,
        sig = replace_na(as.numeric(.tupleColumn(x, qCol)) < threshold, FALSE),
        threshold = threshold
    )
}

#' @rdname getSignificantQtls
#' @param method Correction method whose significant variants to extract:
#'   \code{"permutation"}, \code{"bonferroni_original"},
#'   \code{"bonferroni_filtered"}, or \code{"qvalue"}.
#' @param threshold FDR threshold (defaults to the value stashed by
#'   \code{qtlAssociationPostprocess}).
#' @export
setMethod(
    "getSignificantQtls",
    "QtlSumStats",
    function(
        x,
        method = c(
            "permutation",
            "bonferroni_original",
            "bonferroni_filtered",
            "qvalue"
        ),
        threshold = NULL
    ) {
        method <- arg_match(method)
        masks <- .qapSignificanceMask(x, method, threshold)
        pieces <- compact(
            map(seq_len(nrow(x)), .qapMaskedEntry, x = x, masks = masks)
        )
        if (length(pieces) == 0L) {
            return(x[[1L]][0L])
        }
        exec(c, !!!pieces)
    }
)

# Rebuild a QtlSumStats with extra row columns added + a new qcInfo. Column
# assignment on the S4 subclass itself coerces per-element, so we coerce to a
# plain DFrame, add the columns, and re-wrap (the same rebuild idiom as
# subsetChr).
.qapRebuild <- function(x, newCols, qcInfo) {
    # The collection is a GRangesList now, not a DFrame: the per-variant
    # GRanges are the elements and `newCols` are per-tuple metadata, so the
    # rebuild keeps the elements as-is and only rewrites mcols.
    existing <- mcols(x, use.names = FALSE)
    # A NULL entry in `newCols` removed that column, which is what dropping
    # every named column and re-adding only the non-NULL ones reproduces.
    kept <- existing[,
        setdiff(colnames(existing), names(newCols)),
        drop = FALSE
    ]
    added <- compact(newCols)
    md <- if (length(added) == 0L) {
        kept
    } else {
        cbind(kept, S4Vectors::DataFrame(added, check.names = FALSE))
    }
    # Rebuilding from as.list() starts from the elements' own seqinfo, so the
    # build is written back explicitly -- it is collection-level state, and
    # there is no genome slot to carry it any more.
    built <- .withGenomeBuild(
        GenomicRanges::GRangesList(as.list(x)),
        TRUE,
        getGenome(x)
    )
    grl <- S4Vectors::`mcols<-`(built, value = md)
    methods::new(
        "QtlSumStats",
        grl,
        ldSketch = getLdSketch(x),
        qcInfo = qcInfo
    )
}

#' @rdname qtlAssociationPostprocess
#' @param fdrThreshold Event- and variant-level FDR threshold (default 0.05).
#' @param mafCutoff Minor-allele-frequency cutoff for the FILTERED Bonferroni
#'   flavour (fold of the per-variant \code{af} mcol). \code{0} disables it.
#' @param cisWindow cis-window (bp) for the FILTERED Bonferroni flavour, applied
#'   to the per-variant \code{tss_distance}/\code{tes_distance} mcols. \code{0}
#'   disables it. When either \code{mafCutoff} or \code{cisWindow} is > 0 the
#'   \code{*_filtered} columns are produced (requires the
#'   \code{n_variants_filtered}
#'   row column).
#' @param methods Correction families to compute: any of \code{"permutation"}
#'   (needs \code{p_beta}/\code{beta_shape1}/\code{beta_shape2}) and
#'   \code{"bonferroni"} (needs \code{n_variants}). The q-value SNP method adds
#'   no stored column -- it is a significance query (see
#'   \code{getSignificantQtls}).
#' @param pvalueCol,afCol Entry mcol names for the per-variant p-value / allele
#'   frequency (defaults \code{"P"} / \code{"af"}).
#' @param qvalueArgs Extra arguments for \code{qvalue::qvalue()}, built with
#'   \code{\link{qvalueConfig}} -- the \code{pi0.method} / \code{lambda}
#'   pair in particular, which sets how the null proportion is estimated.
#' @export
setMethod(
    "qtlAssociationPostprocess",
    "QtlSumStats",
    function(
        x,
        fdrThreshold = 0.05,
        mafCutoff = 0,
        cisWindow = 0,
        methods = c("permutation", "bonferroni"),
        pvalueCol = "P",
        afCol = "af",
        qvalueArgs = qvalueConfig()
    ) {
        .assertMethodConfig(qvalueArgs, "qvalueConfig", "qvalueArgs")
        methods <- arg_match(
            methods,
            c("permutation", "bonferroni"),
            multiple = TRUE
        )
        filtering <- (mafCutoff > 0 || cisWindow > 0)
        newCols <- c(
            if (is_in("bonferroni", methods)) {
                .qapBonferroniCols(
                    x,
                    mafCutoff,
                    cisWindow,
                    pvalueCol,
                    afCol,
                    filtering,
                    qvalueArgs
                )
            },
            if (is_in("permutation", methods) && !is.null(x$p_beta)) {
                .qapPermutationCols(x, fdrThreshold, qvalueArgs)
            }
        ) %||%
            list()
        # Stash the correction recipe so getSignificantQtls /
        # annotateSignificance can reproduce significance cheaply (thresholds,
        # not flags).
        qc <- list_assign(
            getQcInfo(x),
            associationPostprocess = .qapRecipe(
                fdrThreshold,
                mafCutoff,
                cisWindow,
                methods,
                pvalueCol,
                afCol
            )
        )
        .qapRebuild(x, newCols, qc)
    }
)

# Local Bonferroni correction columns: per-gene min adjusted p (original and,
# when filtering, filtered) plus their global FDR / q-value adjustments.
# @noRd
.qapBonferroniCols <- function(
    x,
    mafCutoff,
    cisWindow,
    pvalueCol,
    afCol,
    filtering,
    qvalueArgs = list()
) {
    nVar <- .qapBonferroniNVar(x)
    nVarFilt <- if (!is.null(x$n_variants_filtered)) {
        as.numeric(x$n_variants_filtered)
    } else {
        NULL
    }
    if (filtering && is.null(nVarFilt)) {
        msg <- glue(
            "qtlAssociationPostprocess: `n_variants_filtered` row column is ",
            "required when mafCutoff/cisWindow > 0."
        )
        abort(msg)
    }
    perGene <- .qapBonferroniPerGene(
        x,
        nVar,
        nVarFilt,
        pvalueCol,
        afCol,
        mafCutoff,
        cisWindow,
        filtering
    )
    gaO <- .qapGlobalAdjust(perGene$orig, qvalueArgs)
    gaF <- if (filtering) {
        .qapGlobalAdjust(perGene$filt, qvalueArgs)
    } else {
        NULL
    }
    c(
        list(
            p_bonferroni_min_original = perGene$orig,
            fdr_bonferroni_min_original = gaO$fdr,
            q_bonferroni_min_original = gaO$q
        ),
        compact(list(
            p_bonferroni_min_filtered = if (filtering) perGene$filt,
            fdr_bonferroni_min_filtered = gaF$fdr,
            q_bonferroni_min_filtered = gaF$q
        ))
    )
}

# The n_variants row column is mandatory for the Bonferroni correction.
# @noRd
.qapBonferroniNVar <- function(x) {
    if (is.null(x$n_variants)) {
        msg <- glue(
            "qtlAssociationPostprocess: `n_variants` row column is required ",
            "for the Bonferroni correction."
        )
        abort(msg)
    }
    as.numeric(x$n_variants)
}

# Per-gene min Bonferroni-adjusted p (original + filtered). The upstream p-value
# pre-filter always retains the global-min variant, so the min over the entry is
# the exact per-gene min. Returns list(orig, filt).
# @noRd
# One gene's Bonferroni-adjusted minimum p-value, before and after the
# variant filter. NA on either side means the gene had nothing to adjust.
# @noRd
.qapBonferroniForGene <- function(
    i,
    x,
    pvalueCol,
    nVar,
    nVarFilt,
    afCol,
    mafCutoff,
    cisWindow,
    filtering
) {
    mc <- S4Vectors::mcols(x[[i]])
    pv <- mc[[pvalueCol]]
    if (is.null(pv) || length(pv) == 0L) {
        return(list(orig = NA_real_, filt = NA_real_))
    }
    orig <- min(stats::p.adjust(pv, method = "bonferroni", n = nVar[i]))
    if (!filtering) {
        return(list(orig = orig, filt = NA_real_))
    }
    keep <- .qapFilterKeep(mc, mafCutoff, cisWindow, afCol, length(pv))
    if (!any(keep)) {
        return(list(orig = orig, filt = NA_real_))
    }
    list(
        orig = orig,
        filt = min(stats::p.adjust(
            pv[keep],
            method = "bonferroni",
            n = nVarFilt[i]
        ))
    )
}

.qapBonferroniPerGene <- function(
    x,
    nVar,
    nVarFilt,
    pvalueCol,
    afCol,
    mafCutoff,
    cisWindow,
    filtering
) {
    perGene <- map(
        seq_len(nrow(x)),
        .qapBonferroniForGene,
        x = x,
        pvalueCol = pvalueCol,
        nVar = nVar,
        nVarFilt = nVarFilt,
        afCol = afCol,
        mafCutoff = mafCutoff,
        cisWindow = cisWindow,
        filtering = filtering
    )
    list(
        orig = map_dbl(perGene, "orig"),
        filt = map_dbl(perGene, "filt")
    )
}

# Permutation columns: BH-FDR of p_beta, the Storey q-value (q_beta, when
# absent), and the per-gene nominal p threshold from the beta shape params.
# @noRd
.qapPermutationCols <- function(x, fdrThreshold, qvalueArgs = list()) {
    pBeta <- as.numeric(x$p_beta)
    qBeta <- if (!is.null(x$q_beta)) {
        as.numeric(x$q_beta)
    } else {
        .qapSafeQvalue(pBeta, qvalueArgs)
    }
    c(
        compact(list(q_beta = if (is.null(x$q_beta)) qBeta)),
        list(fdr_beta = stats::p.adjust(pBeta, method = "fdr")),
        compact(list(
            p_nominal_threshold = if (
                !is.null(x$beta_shape1) && !is.null(x$beta_shape2)
            ) {
                .qapPermutationNominalThreshold(
                    pBeta,
                    qBeta,
                    as.numeric(x$beta_shape1),
                    as.numeric(x$beta_shape2),
                    fdrThreshold
                )
            }
        ))
    )
}

# The significance recipe stashed for cheap downstream re-derivation.
# @noRd
.qapRecipe <- function(
    fdrThreshold,
    mafCutoff,
    cisWindow,
    methods,
    pvalueCol,
    afCol
) {
    list(
        fdrThreshold = fdrThreshold,
        mafCutoff = mafCutoff,
        cisWindow = cisWindow,
        methods = methods,
        pvalueCol = pvalueCol,
        afCol = afCol
    )
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# An all-FALSE significance mask sized to entry `i`'s p-value column.
# @noRd
.qapEmptyMask <- function(i, x, pcol) {
    rep(FALSE, length(S4Vectors::mcols(x[[i]])[[pcol]]))
}

# Entry `i`'s significant variants (masks[[i]]), labelled with identity mcols;
# NULL when nothing is significant.
# @noRd
.qapMaskedEntry <- function(i, x, masks) {
    gr <- x[[i]][masks[[i]]]
    if (length(gr) == 0L) {
        return(NULL)
    }
    S4Vectors::`mcols<-`(
        gr,
        value = cbind(
            mcols(gr, use.names = FALSE),
            S4Vectors::DataFrame(
                study = as.character(x$study)[i],
                context = as.character(x$context)[i],
                trait = as.character(x$trait)[i]
            )
        )
    )
}
