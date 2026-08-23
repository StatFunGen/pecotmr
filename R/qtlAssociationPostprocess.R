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

# Storey q-values via the Bioconductor qvalue package, with the qvalue-NATIVE
# edge-case retries only (never a hand-rolled substitute): lambda=0 for
# missing/infinite handling, bootstrap pi0 for the "pi0 <= 0" degenerate case.
# Returns a numeric vector aligned to `p`.
.qapSafeQvalue <- function(p) {
    if (!requireNamespace("qvalue", quietly = TRUE)) {
        # nocov start  (optional-package guard; qvalue is Suggests-only)
        msg <- glue(
            "qtlAssociationPostprocess: the 'qvalue' package is required for ",
            "Storey q-values. Install Bioconductor 'qvalue'."
        )
        abort(msg)
        # nocov end
    }
    tryCatch(
        qvalue::qvalue(p)$qvalues,
        error = function(e) {
            if (str_detect(conditionMessage(e), "missing or infinite")) {
                qvalue::qvalue(p, lambda = 0)$qvalues
            } else if (str_detect(conditionMessage(e), "pi0 <= 0")) {
                maxP <- max(p, na.rm = TRUE)
                lambdaSeq <- seq(0, min(0.9, maxP * 0.95), length.out = 10)
                qvalue::qvalue(
                    p,
                    lambda = lambdaSeq,
                    pi0.method = "bootstrap"
                )$qvalues
            } else {
                msg <- glue(
                    "qtlAssociationPostprocess: qvalue::qvalue failed ",
                    "({conditionMessage(e)}). Not substituting a ",
                    "hand-rolled q-value."
                )
                abort(msg)
            }
        }
    )
}

# Event-level global adjustment of a per-gene p-value vector: Benjamini-Hochberg
# FDR (stats::p.adjust) and Storey q (qvalue). Returns list(fdr=, q=).
.qapGlobalAdjust <- function(eventP) {
    list(
        fdr = stats::p.adjust(eventP, method = "fdr"),
        q = .qapSafeQvalue(eventP)
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
    keep <- rep(TRUE, nVar)
    if (mafCutoff > 0 && !is.null(mc[[afCol]])) {
        af <- mc[[afCol]]
        keep <- keep & (pmin(af, 1 - af) > mafCutoff)
    }
    if (
        cisWindow > 0 && !is.null(mc$tss_distance) && !is.null(mc$tes_distance)
    ) {
        keep <- keep &
            (mc$tss_distance >= -cisWindow & mc$tes_distance <= cisWindow)
    }
    keep
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
    thr <- as.numeric(x$p_nominal_threshold)
    for (i in seq_len(nrow(x))) {
        if (is.na(thr[i])) {
            next
        }
        pv <- S4Vectors::mcols(x[[i]])[[pcol]]
        masks[[i]] <- !is.na(pv) & pv < thr[i]
    }
    masks
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
    nVar <- as.numeric(.tupleColumn(x, nCol))
    for (i in seq_len(nrow(x))) {
        masks[[i]] <- .qapBonferroniRowMask(
            x[[i]],
            recipe,
            flav,
            nVar[i],
            varThr
        )
    }
    masks
}

# One entry's Bonferroni keep-mask: Bonferroni-adjusted p <= the global
# threshold, intersected with the filtered-flavour variant filter.
# @noRd
.qapBonferroniRowMask <- function(entry, recipe, flav, nVarI, varThr) {
    mc <- S4Vectors::mcols(entry)
    pv <- mc[[recipe$pvalueCol]]
    keep <- pmin(1, pv * nVarI) <= varThr
    if (flav == "filtered") {
        keep <- keep &
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
.qapMaskQvalue <- function(x, masks, threshold) {
    qCol <- if (!is.null(x$q_beta)) "q_beta" else "q_bonferroni_min_original"
    if (is.null(.tupleColumn(x, qCol))) {
        msg <- glue(
            "getSignificantQtls: no event q-value column (q_beta / ",
            "q_bonferroni_min_original)."
        )
        abort(msg)
    }
    sig <- replace_na(as.numeric(.tupleColumn(x, qCol)) < threshold, FALSE)
    for (i in seq_len(nrow(x))) {
        mc <- S4Vectors::mcols(x[[i]])
        if (isTRUE(sig[i]) && !is.null(mc$qvalue)) {
            masks[[i]] <- !is.na(mc$qvalue) & mc$qvalue < threshold
        }
    }
    masks
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
        pieces <- map(seq_len(nrow(x)), .qapMaskedEntry, x = x, masks = masks)
        pieces <- pieces[!map_lgl(pieces, is.null)]
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
    md <- mcols(x, use.names = FALSE)
    for (nm in names(newCols)) {
        md[[nm]] <- newCols[[nm]]
    }
    grl <- GenomicRanges::GRangesList(as.list(x))
    mcols(grl) <- md
    methods::new(
        "QtlSumStats",
        grl,
        ldSketch = getLdSketch(x),
        genome = getGenome(x),
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
        afCol = "af"
    ) {
        methods <- arg_match(
            methods,
            c("permutation", "bonferroni"),
            multiple = TRUE
        )
        filtering <- (mafCutoff > 0 || cisWindow > 0)
        newCols <- list()
        if (is_in("bonferroni", methods)) {
            newCols <- c(
                newCols,
                .qapBonferroniCols(
                    x,
                    mafCutoff,
                    cisWindow,
                    pvalueCol,
                    afCol,
                    filtering
                )
            )
        }
        if (is_in("permutation", methods) && !is.null(x$p_beta)) {
            newCols <- c(newCols, .qapPermutationCols(x, fdrThreshold))
        }
        # Stash the correction recipe so getSignificantQtls /
        # annotateSignificance can reproduce significance cheaply (thresholds,
        # not flags).
        qc <- getQcInfo(x)
        qc$associationPostprocess <- .qapRecipe(
            fdrThreshold,
            mafCutoff,
            cisWindow,
            methods,
            pvalueCol,
            afCol
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
    filtering
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
    cols <- list(p_bonferroni_min_original = perGene$orig)
    gaO <- .qapGlobalAdjust(perGene$orig)
    cols$fdr_bonferroni_min_original <- gaO$fdr
    cols$q_bonferroni_min_original <- gaO$q
    if (filtering) {
        cols$p_bonferroni_min_filtered <- perGene$filt
        gaF <- .qapGlobalAdjust(perGene$filt)
        cols$fdr_bonferroni_min_filtered <- gaF$fdr
        cols$q_bonferroni_min_filtered <- gaF$q
    }
    cols
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
    n <- nrow(x)
    pBonfOrig <- rep(NA_real_, n)
    pBonfFilt <- rep(NA_real_, n)
    for (i in seq_len(n)) {
        mc <- S4Vectors::mcols(x[[i]])
        pv <- mc[[pvalueCol]]
        if (is.null(pv) || length(pv) == 0L) {
            next
        }
        pBonfOrig[i] <- min(stats::p.adjust(
            pv,
            method = "bonferroni",
            n = nVar[i]
        ))
        if (filtering) {
            keep <- .qapFilterKeep(mc, mafCutoff, cisWindow, afCol, length(pv))
            if (any(keep)) {
                pBonfFilt[i] <- min(stats::p.adjust(
                    pv[keep],
                    method = "bonferroni",
                    n = nVarFilt[i]
                ))
            }
        }
    }
    list(orig = pBonfOrig, filt = pBonfFilt)
}

# Permutation columns: BH-FDR of p_beta, the Storey q-value (q_beta, when
# absent), and the per-gene nominal p threshold from the beta shape params.
# @noRd
.qapPermutationCols <- function(x, fdrThreshold) {
    pBeta <- as.numeric(x$p_beta)
    qBeta <- if (!is.null(x$q_beta)) {
        as.numeric(x$q_beta)
    } else {
        .qapSafeQvalue(pBeta)
    }
    cols <- list()
    if (is.null(x$q_beta)) {
        cols$q_beta <- qBeta
    }
    cols$fdr_beta <- stats::p.adjust(pBeta, method = "fdr")
    if (!is.null(x$beta_shape1) && !is.null(x$beta_shape2)) {
        cols$p_nominal_threshold <- .qapPermutationNominalThreshold(
            pBeta,
            qBeta,
            as.numeric(x$beta_shape1),
            as.numeric(x$beta_shape2),
            fdrThreshold
        )
    }
    cols
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
    gr$study <- as.character(x$study)[i]
    gr$context <- as.character(x$context)[i]
    gr$trait <- as.character(x$trait)[i]
    gr
}
