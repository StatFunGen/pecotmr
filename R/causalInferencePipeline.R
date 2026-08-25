#' @title Causal Inference Pipeline (TWAS-Z + Mendelian Randomization)
#' @description Per-region pipeline that pairs QTL-derived weight vectors
#'   (\code{\link{TwasWeights}} and/or a QTL \code{\link{QtlFineMappingResult}})
#'   with one or more GWAS studies (\code{\link{GwasSumStats}}) to produce
#'   per-tuple TWAS Z-scores and, when fine-mapping is supplied, Wald-ratio
#'   Mendelian Randomization estimates over the QTL credible sets.
#'
#'   Input combinations:
#'   \itemize{
#'     \item \code{twasWeights} alone (no \code{fineMappingResult}):
#'           TWAS Z only, no MR.
#'     \item \code{fineMappingResult} alone (no \code{twasWeights}):
#'           TWAS Z derived from SuSiE-style coefficients carried on the
#'           \code{topLoci} slot of each FineMappingRow; plus MR.
#'     \item both: TWAS Z computed from \code{twasWeights};
#'           MR computed from \code{fineMappingResult}.
#'   }
#'
#' @section LD-sketch identity check: If a QTL input (TwasWeights or
#'   QtlFineMappingResult) carries a non-\code{NULL} \code{ldSketch}, it must
#'   match the \code{ldSketch} on \code{gwasSumStats}. Mismatch is a hard error.
#'   A QTL input with \code{ldSketch = NULL} (the fit was learned from
#'   individual-level data) skips the validation for that input.
#'
#' @section Output shape:
#' A long-format \code{GRanges} with one row per
#' \code{(qtlStudy, context, trait, method, gwasStudy)} tuple,
#' positioned at the variant span of the QTL weight set, with mcols:
#' \describe{
#'   \item{\code{qtlStudy}, \code{context}, \code{trait},
#'         \code{method}, \code{gwasStudy}}{Identity columns.}
#'   \item{\code{twasZ}, \code{twasPval}}{Per-tuple TWAS Z and p-value.}
#'   \item{\code{waldRatio}, \code{waldRatioSe}, \code{mrPval}}{Per-tuple
#'         IVW-aggregated Wald-ratio MR estimate, standard error, and
#'         p-value. Present only when MR was computed; \code{NA}
#'         otherwise.}
#'   \item{\code{nIV}}{Number of instrumental variables used in the
#'         MR aggregation.}
#' }
#'
#' @param gwasSumStats A \code{\link{GwasSumStats}} object. Must be QC'd
#'   (\code{length(getQcInfo(x)) > 0L}).
#' @param twasWeights Optional \code{\link{TwasWeights}} carrying per-(study,
#'   context, trait, method) weights. When supplied, drives the TWAS-Z
#'   computation.
#' @param fineMappingResult Optional \code{\link{QtlFineMappingResult}}. When
#'   supplied, drives the MR computation and (when \code{twasWeights = NULL})
#'   the TWAS-Z weights via the SuSiE-style coefficients on each entry's
#'   \code{topLoci}.
#' @param rsqCutoff Numeric (length 1). When \code{> 0}, performs CV weight
#'   selection (ports the legacy \code{twas_pipeline} \code{pick_best_model} +
#'   \code{update_twas_method}): per \code{(study, context, trait, gwasStudy)}
#'   keep only the method whose \code{cvResult} \code{rsqOption} metric is
#'   highest among methods that clear both \code{rsqCutoff} and the
#'   \code{rsqPvalCutoff} gate AND that produced a finite TWAS Z (the NA/Inf
#'   re-selection); groups where no method clears the cutoffs are dropped. A
#'   group whose methods carry no usable \code{cvResult} (the SS-TWAS path)
#'   keeps all methods. Needs the \code{twasWeights} \code{cvResult}, so
#'   selection is a no-op on the fineMappingResult-only path. Default \code{0}
#'   (no selection; score every method).
#' @param rsqPvalCutoff Numeric (length 1). CV-p-value gate for weight selection
#'   (ports legacy \code{rsq_pval_cutoff}): a method is eligible only when its
#'   \code{cvResult} \code{rsqPvalOption} metric is \code{< rsqPvalCutoff}.
#'   Default \code{Inf} (no p-value gate). A finite value activates selection
#'   even when \code{rsqCutoff = 0}.
#' @param rsqOption Character. Which \code{cvResult} metric is the "r-squared"
#'   used for the cutoff and ranking (ports legacy \code{rsq_option}); typically
#'   \code{"rsq"} or \code{"adj_rsq"}. Default \code{"rsq"}.
#' @param rsqPvalOption Character vector of candidate \code{cvResult} metric
#'   names for the p-value gate (ports legacy \code{rsq_pval_option}); the first
#'   one present in a tuple's metrics is used. Default \code{c("adj_rsq_pval",
#'   "pval")}.
#' @param mrPipCutoff Numeric (length 1). PIP threshold for a \code{topLoci}
#'   variant to be used as an instrumental variable. Used only when
#'   \code{mrMethod = "ivwPerVariant"}. Default \code{0.5}.
#' @param mrMethod One of \code{"ivwPerVariant"} (default) or \code{"csAware"}.
#'   The IVW-per-variant method filters topLoci variants by \code{pip >
#'   mrPipCutoff} and IVW-pools Wald ratios across variants. The CS-aware method
#'   groups variants by credible set (column \code{cs} in topLoci), computes a
#'   PIP-weighted composite Wald ratio per CS using \code{mrCpipCutoff} on the
#'   per-CS cumulative PIP, then IVW-pools across CSs and reports Cochran's Q +
#'   I-squared in the output columns \code{Q}, \code{I2}.
#' @param mrCpipCutoff Numeric (length 1). Cumulative-PIP cutoff for retaining a
#'   credible set. Used only when \code{mrMethod = "csAware"}. Default
#'   \code{0.5}.
#' @param mrPvalCutoff Numeric (length 1). TWAS-p-value gate for running MR
#'   (ports the legacy \code{twas_pipeline} \code{mr_pval_cutoff}): MR is
#'   computed for a \code{(qtl tuple, gwas)} only when its \code{twasPval <
#'   mrPvalCutoff}; otherwise the MR output columns are \code{NA}. Default
#'   \code{1} (no gate; MR runs wherever a \code{fineMappingResult} entry
#'   exists).
#' @param combineMethods Optional character vector forwarded to
#'   \code{\link{combinePValues}} for cross-method combination per
#'   \code{(qtlStudy, context, trait, gwasStudy)} group. \code{NULL} (default)
#'   skips combination.
#' @param alleleFlip Logical, default \code{TRUE}. When TRUE, match QTL variants
#'   to the GWAS by (chrom, pos) with ref/alt swaps recognized and the exposure
#'   effect / weight sign-flipped accordingly; when FALSE, match on exact
#'   alleles only, so a ref/alt swap is treated as a distinct variant.
#' @param ... Reserved.
#' @return A \code{GRanges} as described above.
#' @examples
#' data(qtlDatasetExample)
#' data(gwasSumStatsS4Example)
#' tw <- twasWeightsPipeline(qtlDatasetExample,
#'   methods = "lasso", cisWindow = 1e6)
#' causalInferencePipeline(
#'   gwasSumStats = gwasSumStatsS4Example, twasWeights = tw)
#' @export
causalInferencePipeline <- function(
    gwasSumStats,
    twasWeights = NULL,
    fineMappingResult = NULL,
    rsqCutoff = 0,
    rsqPvalCutoff = Inf,
    rsqOption = "rsq",
    rsqPvalOption = c("adj_rsq_pval", "pval"),
    mrPipCutoff = 0.5,
    mrMethod = c("ivwPerVariant", "csAware"),
    mrCpipCutoff = 0.5,
    mrPvalCutoff = 1,
    combineMethods = NULL,
    alleleFlip = TRUE,
    ...
) {
    mrMethod <- arg_match(mrMethod)
    p <- as.list(environment())
    p$dots <- list(...)
    .cipRun(p)
}

# Orchestrate the causal-inference pipeline over a parameter bundle `p` (from
# as.list(environment())): validate, resolve the QTL work list, apply optional
# CV selection, score every (qtl, gwas) pair, and finalize.
# @noRd
.cipRun <- function(p) {
    .cipValidateInputs(p$gwasSumStats, p$twasWeights, p$fineMappingResult)
    p$gwasLd <- .cipCheckLdSketches(
        p$gwasSumStats,
        p$twasWeights,
        p$fineMappingResult
    )
    p$qtlRows <- .cipResolveWorkList(p$twasWeights, p$fineMappingResult)
    sel <- .cipCvSelection(p)
    p$qtlRows <- sel$qtlRows
    outRows <- list_flatten(map(
        seq_len(nrow(p$qtlRows)),
        .cipScoreQtlTuple,
        p = p
    ))
    if (length(outRows) == 0L) {
        abort(
            "causalInferencePipeline: no (qtl, gwas) tuples produced a result."
        )
    }
    .cipFinalize(outRows, sel, p$combineMethods)
}

# Validate the object classes + QC state of the pipeline inputs.
# @noRd
.cipValidateInputs <- function(gwasSumStats, twasWeights, fineMappingResult) {
    if (!methods::is(gwasSumStats, "GwasSumStats")) {
        abort("`gwasSumStats` must be a GwasSumStats object.")
    }
    if (length(getQcInfo(gwasSumStats)) == 0L) {
        msg <- glue(
            "causalInferencePipeline: gwasSumStats has no QC record ",
            "(getQcInfo() is empty). Call summaryStatsQc() first."
        )
        abort(msg)
    }
    if (is.null(twasWeights) && is.null(fineMappingResult)) {
        msg <- glue(
            "causalInferencePipeline: at least one of `twasWeights` or ",
            "`fineMappingResult` must be supplied."
        )
        abort(msg)
    }
    if (!is.null(twasWeights) && !methods::is(twasWeights, "TwasWeights")) {
        abort("`twasWeights` must be a TwasWeights object or NULL.")
    }
    if (
        !is.null(fineMappingResult) &&
            !methods::is(fineMappingResult, "QtlFineMappingResult")
    ) {
        msg <- glue(
            "`fineMappingResult` must be a QtlFineMappingResult or NULL ",
            "(causalInferencePipeline does not accept GWAS-side fine ",
            "mapping for the QTL slot)."
        )
        abort(msg)
    }
    invisible(NULL)
}

# Require the twasWeights / fineMappingResult LD sketches to match the GWAS one;
# returns the GWAS LD sketch.
# @noRd
.cipCheckLdSketches <- function(gwasSumStats, twasWeights, fineMappingResult) {
    gwasLd <- getLdSketch(gwasSumStats)
    if (!is.null(twasWeights)) {
        .cipRequireMatchingLdSketches(
            getLdSketch(twasWeights),
            gwasLd,
            label = "twasWeights"
        )
    }
    if (!is.null(fineMappingResult)) {
        .cipRequireMatchingLdSketches(
            getLdSketch(fineMappingResult),
            gwasLd,
            label = "fineMappingResult"
        )
    }
    gwasLd
}

# Build the (qtlStudy, context, trait, method) work list + the per-row
# useFmrForWeights flag (error when empty).
# @noRd
.cipResolveWorkList <- function(twasWeights, fineMappingResult) {
    qtlRows <- .cipBuildQtlWorkList(twasWeights, fineMappingResult)
    if (nrow(qtlRows) == 0L) {
        msg <- glue(
            "causalInferencePipeline: no QTL tuples to score (the supplied ",
            "twasWeights / fineMappingResult collections are empty)."
        )
        abort(msg)
    }
    qtlRows$useFmrForWeights <- is.null(twasWeights)
    qtlRows
}

# Optional CV weight selection (legacy pick_best_model + update_twas_method):
# filter to eligible methods now, deferring the final best-method pick to after
# the TWAS Z. Returns list(qtlRows, rsqLookup, selectionActive).
# @noRd
.cipCvSelection <- function(p) {
    selectionActive <- !is.null(p$twasWeights) &&
        (p$rsqCutoff > 0 || is.finite(p$rsqPvalCutoff))
    if (!selectionActive) {
        return(list(
            qtlRows = p$qtlRows,
            rsqLookup = NULL,
            selectionActive = FALSE
        ))
    }
    metricTab <- .cipMethodMetrics(
        p$qtlRows,
        p$twasWeights,
        p$rsqOption,
        p$rsqPvalOption
    )
    rsqLookup <- set_names(
        metricTab$rsq,
        str_c(
            metricTab$qtlStudy,
            metricTab$context,
            metricTab$trait,
            metricTab$method,
            sep = "\r"
        )
    )
    qtlRows <- .cipFilterEligibleMethods(
        p$qtlRows,
        metricTab,
        p$rsqCutoff,
        p$rsqPvalCutoff
    )
    if (nrow(qtlRows) == 0L) {
        msg <- glue(
            "causalInferencePipeline: every QTL tuple was filtered out by ",
            "rsqCutoff = {p$rsqCutoff} / rsqPvalCutoff = {p$rsqPvalCutoff} ",
            "(no method cleared the CV cutoffs)."
        )
        abort(msg)
    }
    list(qtlRows = qtlRows, rsqLookup = rsqLookup, selectionActive = TRUE)
}

# Score one QTL tuple against every GWAS study -> a list of result records
# (empty when the tuple has no usable weights).
# @noRd
.cipScoreQtlTuple <- function(qi, p) {
    qStudy <- p$qtlRows$qtlStudy[[qi]]
    qContext <- p$qtlRows$context[[qi]]
    qTrait <- p$qtlRows$trait[[qi]]
    qMethod <- p$qtlRows$method[[qi]]
    weightsInfo <- .cipExtractWeights(
        twasWeights = p$twasWeights,
        fineMappingResult = p$fineMappingResult,
        study = qStudy,
        context = qContext,
        trait = qTrait,
        method = qMethod,
        useFmr = p$qtlRows$useFmrForWeights[[qi]]
    )
    if (is.null(weightsInfo)) {
        return(list())
    }
    fmrEntry <- .cipResolveFmrEntry(
        p$fineMappingResult,
        qStudy,
        qContext,
        qTrait,
        qMethod
    )
    tuple <- list(
        qStudy = qStudy,
        qContext = qContext,
        qTrait = qTrait,
        qMethod = qMethod
    )
    compact(map(
        seq_len(nrow(p$gwasSumStats)),
        .cipScoreGwasPair,
        tuple = tuple,
        weightsInfo = weightsInfo,
        fmrEntry = fmrEntry,
        p = p
    ))
}

# Resolve the FineMappingRow for a tuple (NULL when the fmr lacks it).
# @noRd
.cipResolveFmrEntry <- function(
    fineMappingResult,
    qStudy,
    qContext,
    qTrait,
    qMethod
) {
    hasTuple <- !is.null(fineMappingResult) &&
        .cipFmrHasTuple(fineMappingResult, qStudy, qContext, qTrait, qMethod)
    if (!hasTuple) {
        return(NULL)
    }
    getFineMappingResult(
        fineMappingResult,
        study = qStudy,
        context = qContext,
        trait = qTrait,
        method = qMethod
    )
}

# Score one (qtl tuple, gwas study) pair -> a result record, or NULL when the
# TWAS Z cannot be computed (too little overlap).
# @noRd
.cipScoreGwasPair <- function(gi, tuple, weightsInfo, fmrEntry, p) {
    gStudy <- as.character(p$gwasSumStats$study)[[gi]]
    gdf <- getSumstatDf(p$gwasSumStats, study = gStudy, require = c("SNP", "Z"))
    twasOut <- .cipComputeTwasZ(
        weights = weightsInfo$weights,
        variantIds = weightsInfo$variantIds,
        gwasDf = gdf,
        gwasLd = p$gwasLd,
        alleleFlip = p$alleleFlip,
        label = .cipPairLabel(tuple, gStudy)
    )
    if (is.null(twasOut)) {
        return(NULL)
    }
    mrOut <- .cipRunMr(fmrEntry, gdf, twasOut, p)
    .cipResultRow(tuple, gStudy, twasOut, mrOut)
}

# Identify a (QTL tuple, GWAS study) pair in diagnostics.
# @noRd
.cipPairLabel <- function(tuple, gStudy) {
    glue(
        "QTL (study='{tuple$qStudy}', context='{tuple$qContext}', ",
        "trait='{tuple$qTrait}', method='{tuple$qMethod}') x GWAS ",
        "(study='{gStudy}')"
    )
}

# Run MR for a pair, gated on the TWAS p-value (mrPvalCutoff >= 1 disables the
# gate) and the presence of a fine-mapping entry.
# @noRd
.cipRunMr <- function(fmrEntry, gdf, twasOut, p) {
    mrGateOpen <- p$mrPvalCutoff >= 1 ||
        (!is.na(twasOut$pval) && twasOut$pval < p$mrPvalCutoff)
    if (is.null(fmrEntry) || !mrGateOpen) {
        return(.cipEmptyMr())
    }
    if (p$mrMethod == "csAware") {
        .cipComputeMrCsAware(
            fmrEntry = fmrEntry,
            gwasDf = gdf,
            cpipCutoff = p$mrCpipCutoff,
            alleleFlip = p$alleleFlip
        )
    } else {
        .cipComputeMr(
            fmrEntry = fmrEntry,
            gwasDf = gdf,
            pipCutoff = p$mrPipCutoff,
            alleleFlip = p$alleleFlip
        )
    }
}

# The all-NA MR result used when MR is skipped for a pair.
# @noRd
.cipEmptyMr <- function() {
    list(
        waldRatio = NA_real_,
        waldRatioSe = NA_real_,
        mrPval = NA_real_,
        nIV = NA_integer_,
        Q = NA_real_,
        I2 = NA_real_,
        nCs = NA_integer_
    )
}

# Assemble one result record from a tuple + its TWAS / MR outputs.
# @noRd
.cipResultRow <- function(tuple, gStudy, twasOut, mrOut) {
    list(
        qtlStudy = tuple$qStudy,
        context = tuple$qContext,
        trait = tuple$qTrait,
        method = tuple$qMethod,
        gwasStudy = gStudy,
        twasZ = twasOut$Z,
        twasPval = twasOut$pval,
        waldRatio = mrOut$waldRatio,
        waldRatioSe = mrOut$waldRatioSe,
        mrPval = mrOut$mrPval,
        nIV = mrOut$nIV,
        Q = mrOut$Q %||% NA_real_,
        I2 = mrOut$I2 %||% NA_real_,
        nCs = mrOut$nCs %||% NA_integer_,
        chrom = twasOut$chrom,
        startPos = twasOut$startPos,
        endPos = twasOut$endPos
    )
}

# Assemble the result data.frame, apply the final best-method pick, convert to
# GRanges, and (optionally) combine across methods.
# @noRd
.cipFinalize <- function(outRows, sel, combineMethods) {
    resultDf <- .cipRowsToDf(outRows)
    # Final best-method pick + NA/Inf re-selection (legacy update_twas_method):
    # per (qtlStudy, context, trait, gwasStudy) keep the highest-rsqOption
    # eligible method whose TWAS Z is finite, falling back to the top-rsq method
    # when none is finite. SS-TWAS groups (no usable rsq) keep all methods.
    if (sel$selectionActive) {
        resultDf <- .cipSelectBestMethod(resultDf, sel$rsqLookup)
    }
    out <- .cipDfToGranges(resultDf)
    if (!is.null(combineMethods)) {
        out <- .cipCombineAcrossMethods(out, methods = combineMethods)
    }
    out
}

# =============================================================================
# Internal helpers
# =============================================================================

# Compare two GenotypeHandles for LD-sketch identity. Thin wrapper over
# the shared `.requireMatchingLdSketches` helper (R/ld.R).
.cipRequireMatchingLdSketches <- function(qtlLd, gwasLd, label) {
    .requireMatchingLdSketches(
        qtlLd,
        gwasLd,
        pipelineName = "causalInferencePipeline",
        label = label
    )
}

# Build the (qtlStudy, context, trait, method) work list from
# whichever input was supplied. When both are supplied, prefer the
# TwasWeights tuples and only retain those that also appear in the FMR
# (so MR has something to attach).
.cipBuildQtlWorkList <- function(twasWeights, fineMappingResult) {
    if (!is.null(twasWeights)) {
        df <- tibble(
            qtlStudy = as.character(twasWeights$study),
            context = as.character(twasWeights$context),
            trait = as.character(twasWeights$trait),
            method = as.character(twasWeights$method)
        )
    } else {
        df <- tibble(
            qtlStudy = as.character(fineMappingResult$study),
            context = as.character(fineMappingResult$context),
            trait = as.character(fineMappingResult$trait),
            method = as.character(fineMappingResult$method)
        )
    }
    df
}

# Resolve one CV metric (rsqOption / rsqPvalOption) for a single tuple from the
# TwasWeights cvResult, which the individual-level CV path stores as a list
# with a named $metrics vector (corr, rsq, adj_rsq, pval, RMSE, MAE); a bare
# metrics vector / data frame is tolerated too. `which` is a vector of candidate
# metric names; the first present is used. Returns NA when no usable metric.
# @noRd
.cipCvMetric <- function(twasWeights, study, context, trait, method, which) {
    perf <- tryCatch(
        getCvResult(
            twasWeights,
            study = study,
            context = context,
            trait = trait,
            method = method
        ),
        error = function(e) NULL
    )
    if (is.null(perf)) {
        return(NA_real_)
    }
    metrics <- if (is.list(perf) && !is.null(perf[["metrics"]])) {
        perf[["metrics"]]
    } else {
        perf
    }
    nm <- intersect(which, names(metrics))
    if (length(nm) == 0L) {
        return(NA_real_)
    }
    val <- suppressWarnings(as.numeric(metrics[[nm[[1L]]]]))
    if (length(val) == 0L) NA_real_ else val[[1L]]
}

# Tabulate the rsqOption (rsq) and rsqPvalOption (pval) CV metrics for every
# tuple in the work-list. Returns the identity columns plus `rsq`, `pval`.
# @noRd
.cipMethodMetrics <- function(qtlRows, twasWeights, rsqOption, rsqPvalOption) {
    n <- nrow(qtlRows)
    rsq <- map_dbl(
        seq_len(n),
        .cipRowCvMetric,
        twasWeights = twasWeights,
        qtlRows = qtlRows,
        which = rsqOption
    )
    pval <- map_dbl(
        seq_len(n),
        .cipRowCvMetric,
        twasWeights = twasWeights,
        qtlRows = qtlRows,
        which = rsqPvalOption
    )
    tibble(
        qtlStudy = qtlRows$qtlStudy,
        context = qtlRows$context,
        trait = qtlRows$trait,
        method = qtlRows$method,
        rsq = rsq,
        pval = pval
    )
}

# Eligibility filter (legacy pick_best_model gate): per (study, context, trait)
# keep methods whose rsq >= rsqCutoff and (when the gate is finite) whose CV
# p-value < rsqPvalCutoff. A group whose methods carry no usable rsq (all NA) is
# the SS-TWAS path and keeps all its methods. Groups where no method clears the
# cutoffs contribute nothing. Returns the filtered work-list (same columns).
# @noRd
.cipFilterEligibleMethods <- function(
    qtlRows,
    metricTab,
    rsqCutoff,
    rsqPvalCutoff
) {
    grp <- str_c(
        metricTab$qtlStudy,
        metricTab$context,
        metricTab$trait,
        sep = "\r"
    )
    keep <- logical(nrow(qtlRows))
    pvalGate <- is.finite(rsqPvalCutoff)
    for (g in unique(grp)) {
        idx <- which(grp == g)
        rsq <- metricTab$rsq[idx]
        if (all(is.na(rsq))) {
            keep[idx] <- TRUE
            next
        } # SS-TWAS: keep all
        elig <- !is.na(rsq) & rsq >= rsqCutoff
        if (pvalGate) {
            elig <- elig &
                !is.na(metricTab$pval[idx]) &
                metricTab$pval[idx] < rsqPvalCutoff
        }
        keep[idx[elig]] <- TRUE
    }
    filter(qtlRows, keep)
}

# Final best-method pick + NA/Inf re-selection (legacy update_twas_method): per
# (qtlStudy, context, trait, gwasStudy) rank the (already-eligible) methods by
# rsqLookup descending and keep the first whose twasZ is finite; if none is
# finite, keep the top-rsq method. A group with no usable rsq (SS-TWAS) keeps
# all its rows. `rsqLookup` is keyed by study\rcontext\rtrait\rmethod.
# @noRd
.cipSelectBestMethod <- function(df, rsqLookup) {
    if (nrow(df) == 0L) {
        return(df)
    }
    key <- str_c(df$qtlStudy, df$context, df$trait, df$method, sep = "\r")
    rsq <- unname(rsqLookup[key])
    grp <- str_c(df$qtlStudy, df$context, df$trait, df$gwasStudy, sep = "\r")
    keepRow <- logical(nrow(df))
    for (g in unique(grp)) {
        idx <- which(grp == g)
        r <- rsq[idx]
        if (all(is.na(r))) {
            keepRow[idx] <- TRUE
            next
        } # SS-TWAS: keep all
        ord <- idx[order(r, decreasing = TRUE)] # NA sorts last
        z <- suppressWarnings(as.numeric(df$twasZ[ord]))
        fin <- which(is.finite(z))
        sel <- if (length(fin) > 0L) ord[[fin[[1L]]]] else ord[[1L]]
        keepRow[sel] <- TRUE
    }
    filter(df, keepRow)
}

.cipFmrHasTuple <- function(fmr, study, context, trait, method) {
    length(.matchTupleRows(
        fmr,
        list(study = study, context = context, trait = trait, method = method)
    )) >
        0L
}

# Fetch the per-tuple weight-source entry (a TwasWeightsRow from
# `twasWeights`, or a FineMappingRow from `fineMappingResult`) and resolve it
# to an aligned (variantIds, weights) pair via the shared `resolveWeights` --
# the SuSiE-style posterior effect for the FMR path. Returns NULL when the tuple
# is absent or has no usable weights.
.cipExtractWeights <- function(
    twasWeights,
    fineMappingResult,
    study,
    context,
    trait,
    method,
    useFmr
) {
    ent <- if (!useFmr) {
        .cipWeightsFromTwas(twasWeights, study, context, trait, method)
    } else {
        .cipWeightsFromFmr(fineMappingResult, study, context, trait, method)
    }
    if (is.null(ent)) {
        return(NULL)
    }
    wr <- resolveWeights(ent)
    if (length(wr$variantIds) == 0L) {
        return(NULL)
    }
    wr
}

# TwasWeights entry for a tuple (NULL when absent).
# @noRd
.cipWeightsFromTwas <- function(twasWeights, study, context, trait, method) {
    tuple <- list(
        study = study,
        context = context,
        trait = trait,
        method = method
    )
    if (length(.matchTupleRows(twasWeights, tuple)) == 0L) {
        return(NULL)
    }
    getTwasWeights(
        twasWeights,
        study = study,
        context = context,
        trait = trait,
        method = method
    )
}

# FineMappingResult entry for a tuple (NULL when absent).
# @noRd
.cipWeightsFromFmr <- function(
    fineMappingResult,
    study,
    context,
    trait,
    method
) {
    if (!.cipFmrHasTuple(fineMappingResult, study, context, trait, method)) {
        return(NULL)
    }
    getFineMappingResult(
        fineMappingResult,
        study = study,
        context = context,
        trait = trait,
        method = method
    )
}

# Compute the per-tuple TWAS Z from a single GwasSumStats tuple's
# unpacked data.frame (produced by getSumstatDf upstream). Returns
# NULL when the overlap is too small.
.cipComputeTwasZ <- function(
    weights,
    variantIds,
    gwasDf,
    gwasLd,
    alleleFlip = TRUE,
    label = "this (QTL, GWAS) pair"
) {
    # Match QTL weights to the GWAS sumstats by (chrom, pos, allele) rather than
    # by raw id string, so chr-prefix / separator / allele-swap differences are
    # reconciled instead of silently dropped. sign = -1 marks an allele swap;
    # weights[idxA] * sign brings the weight into the GWAS allele coding, which
    # the GWAS z and the GWAS-panel LD already share.
    m <- matchVariants(variantIds, gwasDf$variant_id, allowFlip = alleleFlip)
    # Reconciliation here is one-sided and unbiased -- subsetting both w and R
    # keeps Z = w'z / sqrt(w'Rw) honest, it only costs power -- but the caller
    # cannot see how much power was lost unless the drop is reported.
    .cipWarnVariantDrop(length(variantIds), length(m$idxA), label)
    if (length(m$idxA) < 2L) {
        return(NULL)
    }

    wSub <- weights[m$idxA] * m$sign
    zSub <- gwasDf$z[m$idxB]
    gwasIds <- gwasDf$variant_id[m$idxB]
    ldMat <- .cipLdFromSketch(gwasLd, gwasIds)

    res <- twasZ(weights = wSub, z = zSub, R = ldMat)
    zMat <- res$Z
    zVal <- as.numeric(zMat[1L, "Z"])
    pVal <- as.numeric(zMat[1L, "pval"])
    # Position the row at the variant span.
    chrom <- gwasDf$chrom[[m$idxB[[1L]]]]
    startPos <- min(gwasDf$pos[m$idxB])
    endPos <- max(gwasDf$pos[m$idxB])
    list(
        Z = zVal,
        pval = pVal,
        chrom = chrom,
        startPos = startPos,
        endPos = endPos
    )
}

# Report how many QTL weight variants failed to match the GWAS sumstats.
# Dropping below two matches yields no TWAS Z at all, which used to be a silent
# NULL, so that case is called out separately.
# @noRd
.cipWarnVariantDrop <- function(nWeights, nMatched, label) {
    if (nMatched >= nWeights) {
        return(invisible(NULL))
    }
    dropped <- nWeights - nMatched
    if (nMatched < 2L) {
        msg <- glue(
            "causalInferencePipeline: {label} matched only {nMatched} of ",
            "{nWeights} weight variant(s) to the GWAS sumstats; ",
            "at least 2 are needed for a TWAS Z, so no result is reported."
        )
    } else {
        msg <- glue(
            "causalInferencePipeline: {label} dropped {dropped} of ",
            "{nWeights} weight variant(s) absent from the GWAS sumstats; ",
            "the TWAS Z uses the remaining {nMatched}."
        )
    }
    warn(msg)
}

# Build an LD correlation matrix for a given variant subset from an
# LD sketch. Thin wrapper over the shared `.ldFromSketch` helper
# (R/ld.R).
.cipLdFromSketch <- function(ldSketch, variantIds) {
    .ldFromSketch(ldSketch, variantIds, label = "causalInferencePipeline")
}

# Compute the Wald-ratio IVW MR estimate for a single tuple. Uses the
# FineMappingRow's topLoci as the instrumental variable source: each
# variant with PIP > pipCutoff contributes one ratio = beta_y / beta_x.
# Returns list(waldRatio, waldRatioSe, mrPval, nIV) with NA fields when
# no IVs survive.
.cipComputeMr <- function(fmrEntry, gwasDf, pipCutoff, alleleFlip = TRUE) {
    iv <- .cipMrInstruments(fmrEntry, pipCutoff)
    if (is.null(iv)) {
        return(.cipEmptyMrRatio())
    }
    m <- matchVariants(iv$ivVars, gwasDf$variant_id, allowFlip = alleleFlip)
    if (length(m$idxA) == 0L) {
        return(.cipEmptyMrRatio())
    }
    # Align the QTL exposure effect to the GWAS allele coding (sign = -1 on an
    # allele swap) so the Wald ratio betaY / betaX has the correct sign.
    betaX <- iv$betaX[m$idxA] * m$sign
    seX <- iv$seX[m$idxA]
    gIdx <- m$idxB
    gZ <- gwasDf$z[gIdx]
    gN <- .cipGwasCol(gwasDf$N, gIdx)
    gMaf <- .cipGwasCol(gwasDf$maf, gIdx)
    betaY <- .cipZToBeta(gZ, gMaf, gN)
    seY <- .cipZToSe(gZ, gMaf, gN)
    ratio <- betaY / betaX
    # Standard Wald-ratio SE via delta method.
    rSe <- sqrt((seY / betaX)^2 + (betaY * seX / betaX^2)^2)
    # IVW pooling of the per-instrument Wald ratios.
    pooled <- .ivwPool(ratio, rSe)
    if (pooled$n == 0L) {
        return(.cipEmptyMrRatio())
    }
    list(
        waldRatio = pooled$effect,
        waldRatioSe = pooled$se,
        mrPval = pooled$pval,
        nIV = pooled$n
    )
}

# The all-NA per-variant Wald-ratio MR result (4 fields).
# @noRd
.cipEmptyMrRatio <- function() {
    list(
        waldRatio = NA_real_,
        waldRatioSe = NA_real_,
        mrPval = NA_real_,
        nIV = 0L
    )
}

# Instrumental variables from a FineMappingRow's topLoci: variants with PIP >
# pipCutoff. Returns list(ivVars, betaX, seX) or NULL when none qualify.
# @noRd
.cipMrInstruments <- function(fmrEntry, pipCutoff) {
    tl <- getTopLoci(fmrEntry)
    if (is.null(tl) || nrow(tl) == 0L) {
        return(NULL)
    }
    cols <- .cipTlCols(tl)
    if (
        length(cols$pip) == 0L ||
            length(cols$beta) == 0L ||
            length(cols$se) == 0L
    ) {
        return(NULL)
    }
    keep <- !is.na(tl[[cols$pip[[1L]]]]) & tl[[cols$pip[[1L]]]] > pipCutoff
    ivVars <- as.character(tl$variant_id)[keep]
    if (length(ivVars) == 0L) {
        return(NULL)
    }
    list(
        ivVars = ivVars,
        betaX = as.numeric(tl[[cols$beta[[1L]]]])[keep],
        seX = as.numeric(tl[[cols$se[[1L]]]])[keep]
    )
}

# Index a possibly-absent GWAS column, filling NA when the column is missing.
# @noRd
.cipGwasCol <- function(col, gIdx) {
    if (!is.null(col)) col[gIdx] else rep(NA_real_, length(gIdx))
}

# Cochran's Q-based I-squared heterogeneity statistic. Ported from
# the legacy mr.R::calcI2. Q = 0 or near-zero -> I2 = 0; clipped to
# [0, 1] to stay in the usual heterogeneity convention.
# @noRd
.cipCalcI2 <- function(Q, nGroups) {
    if (!is.finite(Q) || Q <= 1e-3 || nGroups <= 1L) {
        return(0)
    }
    i2 <- (Q - (nGroups - 1)) / Q
    max(0, min(1, i2))
}

# CS-aware Wald-ratio MR estimate for a single tuple. Adapted from the
# legacy mr.R::mrAnalysis but operating on a `FineMappingRow`'s topLoci
# data.frame + a GWAS `GRanges` instead of the old data.frame-only API.
#
# Logic:
#   1. Group topLoci variants by credible set (column `cs` in topLoci, or
#      first column matching ^cs).
#   2. Per CS, compute cumulative PIP; drop CSs with cpip < cpipCutoff.
#   3. Per surviving CS, compute a PIP-weighted composite Wald ratio
#      `composite_bhat = sum((bhatY/bhatX * pip) / cpip)` with the
#      delta-method composite SE.
#   4. IVW-pool composite Walds across CSs with weights wv = 1/se^2.
#   5. Report (metaEff, metaSe, metaPval, Q, I2, nCs, nIv) where the
#      result list keys are renamed (waldRatio/waldRatioSe/mrPval/nIV) so
#      the calling pipeline emits consistent column names.
#
# Returns list(waldRatio, waldRatioSe, mrPval, nIV, Q, I2, nCs) with NA
# fields when no usable CS survives.
# @noRd
.cipComputeMrCsAware <- function(
    fmrEntry,
    gwasDf,
    cpipCutoff,
    alleleFlip = TRUE
) {
    inst <- .cipCsAwareInstruments(fmrEntry, gwasDf, alleleFlip)
    if (is.null(inst)) {
        return(.cipEmptyMrCs())
    }
    comp <- .cipCsComposites(inst, cpipCutoff)
    if (length(comp$bhat) == 0L) {
        return(.cipEmptyMrCs())
    }
    .cipCsAwareResult(comp, inst)
}

# The all-NA CS-aware MR result (7 fields).
# @noRd
.cipEmptyMrCs <- function() {
    list(
        waldRatio = NA_real_,
        waldRatioSe = NA_real_,
        mrPval = NA_real_,
        nIV = 0L,
        Q = NA_real_,
        I2 = NA_real_,
        nCs = 0L
    )
}

# Extract + GWAS-align the CS-aware exposure instruments. Returns the aligned,
# unit-SE-standardized (cs, pip, bhatX, sbhatX, bhatY, sbhatY) or NULL.
# @noRd
.cipCsAwareInstruments <- function(fmrEntry, gwasDf, alleleFlip) {
    tl <- getTopLoci(fmrEntry)
    if (is.null(tl) || nrow(tl) == 0L) {
        return(NULL)
    }
    cols <- .cipTlCols(tl)
    if (
        length(cols$pip) == 0L ||
            length(cols$beta) == 0L ||
            length(cols$se) == 0L ||
            length(cols$cs) == 0L
    ) {
        return(NULL)
    }
    raw <- .cipCsAwareRaw(tl, cols)
    if (is.null(raw)) {
        return(NULL)
    }
    .cipCsAwareAlign(raw, gwasDf, alleleFlip)
}

# Pull the CS / PIP / beta / se / variant columns from topLoci, dropping rows
# with NA cs / non-positive cs / NA pip / NA beta. Returns NULL when none valid.
# @noRd
.cipCsAwareRaw <- function(tl, cols) {
    cs <- as.integer(tl[[cols$cs[[1L]]]])
    pip <- as.numeric(tl[[cols$pip[[1L]]]])
    bhatX <- as.numeric(tl[[cols$beta[[1L]]]])
    sbhatX <- as.numeric(tl[[cols$se[[1L]]]])
    vids <- as.character(tl$variant_id)
    ok <- !is.na(cs) & cs > 0L & !is.na(pip) & !is.na(bhatX) & !is.na(sbhatX)
    if (!any(ok)) {
        return(NULL)
    }
    list(
        cs = cs[ok],
        pip = pip[ok],
        bhatX = bhatX[ok],
        sbhatX = sbhatX[ok],
        vids = vids[ok]
    )
}

# Match the exposure variants to the GWAS sumstats (allele-aware), derive the
# outcome beta/se, and standardize the exposure to unit SE. Returns NULL when
# nothing matches.
# @noRd
.cipCsAwareAlign <- function(raw, gwasDf, alleleFlip) {
    m <- matchVariants(raw$vids, gwasDf$variant_id, allowFlip = alleleFlip)
    if (length(m$idxA) == 0L) {
        return(NULL)
    }
    gIdx <- m$idxB
    gZ <- gwasDf$z[gIdx]
    gN <- .cipGwasCol(gwasDf$N, gIdx)
    gMaf <- .cipGwasCol(gwasDf$maf, gIdx)
    # Standardize bhatX -> z; sbhatX -> 1 (legacy mrAnalysis rescaling).
    bhatX <- (raw$bhatX[m$idxA] * m$sign) / raw$sbhatX[m$idxA]
    list(
        cs = raw$cs[m$idxA],
        pip = raw$pip[m$idxA],
        bhatX = bhatX,
        sbhatX = rep(1, length(bhatX)),
        bhatY = .cipZToBeta(gZ, gMaf, gN),
        sbhatY = .cipZToSe(gZ, gMaf, gN)
    )
}

# Per-CS PIP-weighted composite Wald ratios (dropping CSs below cpipCutoff or
# with a degenerate SE). Returns list(bhat, sbhat) vectors.
# @noRd
.cipCsComposites <- function(inst, cpipCutoff) {
    csIds <- sort(unique(inst$cs))
    comps <- compact(map(
        csIds,
        .cipCsComposite,
        inst = inst,
        cpipCutoff = cpipCutoff
    ))
    list(bhat = map_dbl(comps, "bhat"), sbhat = map_dbl(comps, "sbhat"))
}

# One CS's composite Wald ratio + delta-method SE, or NULL when the CS is
# dropped (cpip below cutoff, or non-finite / non-positive SE).
# @noRd
.cipCsComposite <- function(cid, inst, cpipCutoff) {
    inCs <- which(inst$cs == cid)
    cpip <- sum(inst$pip[inCs])
    if (!is.finite(cpip) || cpip < cpipCutoff) {
        return(NULL)
    }
    pNorm <- inst$pip[inCs] / cpip
    bx <- inst$bhatX[inCs]
    sx <- inst$sbhatX[inCs]
    by <- inst$bhatY[inCs]
    sy <- inst$sbhatY[inCs]
    cBhat <- sum((by / bx) * pNorm)
    second <- sum(
        ((by / bx)^2 + (sy^2 / bx^2) + ((by^2 * sx^2) / bx^4)) * pNorm
    )
    cSe <- sqrt(max(0, second - cBhat^2))
    if (!is.finite(cBhat) || !is.finite(cSe) || cSe <= 0) {
        return(NULL)
    }
    list(bhat = cBhat, sbhat = cSe)
}

# IVW-meta the per-CS composite Walds into the final CS-aware MR result.
# @noRd
.cipCsAwareResult <- function(comp, inst) {
    pooled <- .ivwPool(comp$bhat, comp$sbhat)
    list(
        waldRatio = pooled$effect,
        waldRatioSe = pooled$se,
        mrPval = pooled$pval,
        nIV = length(inst$cs),
        Q = pooled$Q,
        I2 = .cipCalcI2(pooled$Q, pooled$n),
        nCs = length(comp$bhat)
    )
}


# Derive beta / se from z using maf + n via the shared zToBetaSe() (model-exact
# se = 1/sqrt(2*p*q*(N + z^2)), beta = z*se). Fall back to z as a beta surrogate
# / se = 1 when maf or n is unavailable.
.cipZToBeta <- function(z, maf, n) {
    if (any(is.na(maf)) || any(is.na(n))) {
        return(z)
    } # fall back to z as a beta surrogate when no maf/n
    .zToBetaSe(z, maf, n)$beta
}
.cipZToSe <- function(z, maf, n) {
    if (any(is.na(maf)) || any(is.na(n))) {
        return(rep(1, length(z)))
    }
    .zToBetaSe(z, maf, n)$se
}

# Fixed-effect inverse-variance-weighted (IVW) pooling of per-instrument
# effects. Weights 1/se^2, drops non-finite / non-positive weights, and
# returns the pooled effect, its SE (1/sqrt(sum w)), two-tailed p-value,
# Cochran's Q, and the number pooled (n = 0 when nothing is poolable).
.ivwPool <- function(effect, se) {
    w <- 1 / se^2
    ok <- is.finite(w) & w > 0
    if (!any(ok)) {
        return(list(
            effect = NA_real_,
            se = NA_real_,
            pval = NA_real_,
            Q = NA_real_,
            n = 0L
        ))
    }
    effect <- effect[ok]
    w <- w[ok]
    metaEff <- sum(w * effect) / sum(w)
    metaSe <- 1 / sqrt(sum(w))
    list(
        effect = metaEff,
        se = metaSe,
        pval = .zToPvalue(metaEff / metaSe),
        Q = sum(w * (effect - metaEff)^2),
        n = length(effect)
    )
}

# Resolve the topLoci column aliases used across the CIP MR helpers. Returns the
# first matching column name for each role (character(0) if absent); `cs` also
# falls back to any column whose name starts with "cs" (e.g. cs_0.95).
.cipTlCols <- function(tl) {
    cn <- colnames(tl)
    cs <- intersect("cs", cn)
    if (length(cs) == 0L) {
        csCand <- cn[str_detect(cn, "^cs")]
        if (length(csCand) > 0L) cs <- csCand[[1L]]
    }
    list(
        pip = intersect(c("pip", "PIP"), cn),
        beta = intersect(c("betahat", "beta", "bhat_x"), cn),
        se = intersect(c("sebetahat", "se", "sbhat_x"), cn),
        cs = cs
    )
}

# Convert the accumulated list of row records to a flat data.frame.
.cipRowsToDf <- function(rows) {
    bind_rows(rows)
}

# Convert the assembled result data.frame to a GRanges with mcols.
.cipDfToGranges <- function(df) {
    chr <- withChrPrefix(df$chrom)
    gr <- GenomicRanges::GRanges(
        seqnames = chr,
        ranges = IRanges::IRanges(
            start = as.integer(df$startPos),
            end = as.integer(df$endPos)
        )
    )
    mcols <- select(
        df,
        all_of(c(
            "qtlStudy",
            "context",
            "trait",
            "method",
            "gwasStudy",
            "twasZ",
            "twasPval",
            "waldRatio",
            "waldRatioSe",
            "mrPval",
            "nIV",
            "Q",
            "I2",
            "nCs"
        ))
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(mcols)
    gr
}

# Combine TWAS p-values across method for each (qtlStudy, context,
# trait, gwasStudy) group. Appends one row per group with a
# combined p-value and methodName = "combined.<methodToken>". Uses
# combinePValues() with the cross-method correlation set to the identity
# (we have no cross-method covariance available downstream).
.cipCombineAcrossMethods <- function(gr, methods) {
    mc <- as_tibble(as.data.frame(S4Vectors::mcols(gr)))
    key <- str_c(mc$qtlStudy, mc$context, mc$trait, mc$gwasStudy, sep = "||")
    groups <- split(seq_len(nrow(mc)), key)
    extras <- list()
    for (gkey in names(groups)) {
        rows <- groups[[gkey]]
        if (length(rows) < 2L) {
            next
        }
        pvals <- as.numeric(mc$twasPval[rows])
        zvec <- as.numeric(mc$twasZ[rows])
        cp <- combinePValues(pvals = pvals, zScores = zvec, methods = methods)
        for (m in methods) {
            newRow <- slice(mc, rows[[1L]])
            newRow$method <- str_c("combined.", m)
            newRow$twasZ <- NA_real_
            newRow$twasPval <- as.numeric(cp$results[[m]]$pval)
            newRow$waldRatio <- NA_real_
            newRow$waldRatioSe <- NA_real_
            newRow$mrPval <- NA_real_
            newRow$nIV <- NA_integer_
            extras[[length(extras) + 1L]] <- newRow
        }
    }
    if (length(extras) == 0L) {
        return(gr)
    }
    newMcs <- bind_rows(extras)
    newGr <- gr[rep(1L, nrow(newMcs))]
    S4Vectors::mcols(newGr) <- S4Vectors::DataFrame(newMcs)
    c(gr, newGr)
}


# =============================================================================
# Unified TWAS Z-statistic
# =============================================================================

# Internal: build the K x K covariance W^T R W. Uses the SVD path when the
# triplet (V, D, nSketch) is supplied, otherwise the R / X path. Aligns LD
# rows/cols to the rownames of W when both are named; falls back to
# positional alignment otherwise.
.twasZCovY <- function(
    weights,
    R = NULL,
    X = NULL,
    V = NULL,
    D = NULL,
    nSketch = NULL
) {
    rn <- rownames(weights)
    if (!is.null(V) && !is.null(D) && !is.null(nSketch)) {
        return(.twasZCovYSvd(weights, V, D, nSketch, rn))
    }
    .twasZCovYR(weights, R, X, rn)
}

# SVD path: covY = (V^TW * sqrt(Lambda))^T (V^TW * sqrt(Lambda)) with
# Lambda_i = D_i^2 / (nSketch - 1).
# @noRd
.twasZCovYSvd <- function(weights, V, D, nSketch, rn) {
    vSub <- .twasZAlignV(V, rn, nrow(weights))
    Lambda <- D^2 / (nSketch - 1)
    VtW <- crossprod(vSub, weights)
    list(covY = crossprod(VtW * sqrt(Lambda)))
}

# Row-align V to the weights (by name when available, else positionally).
# @noRd
.twasZAlignV <- function(V, rn, nW) {
    if (!is.null(rownames(V)) && !is.null(rn)) {
        idx <- match(rn, rownames(V))
        if (anyNA(idx)) {
            msg <- glue(
                "twasZ: V is missing rows for {sum(is.na(idx))} ",
                "variant(s) named in weights."
            )
            abort(msg)
        }
        return(V[idx, , drop = FALSE])
    }
    if (nrow(V) != nW) {
        abort(
            "twasZ: positional alignment requires nrow(V) == nrow(weights)."
        )
    }
    V
}

# R / X path: covY = W^T R W (R computed from X when R is absent).
# @noRd
.twasZCovYR <- function(weights, R, X, rn) {
    if (is.null(R)) {
        if (is.null(X)) {
            abort("twasZ: provide R, X, or the (V, D, nSketch) SVD triplet.")
        }
        R <- computeLd(X)
    }
    rSub <- .twasZAlignR(R, rn, nrow(weights))
    list(covY = crossprod(weights, rSub) %*% weights)
}

# Symmetric-align R to the weights (by name when available, else positionally).
# @noRd
.twasZAlignR <- function(R, rn, nW) {
    if (!is.null(rownames(R)) && !is.null(rn)) {
        idx <- match(rn, rownames(R))
        if (anyNA(idx)) {
            msg <- glue(
                "twasZ: R is missing rows for {sum(is.na(idx))} ",
                "variant(s) named in weights."
            )
            abort(msg)
        }
        return(R[idx, idx, drop = FALSE])
    }
    if (nrow(R) != nW) {
        abort(
            "twasZ: positional alignment requires nrow(R) == nrow(weights)."
        )
    }
    R
}

#' Calculate TWAS Z-Statistics for One or More Methods / Contexts
#'
#' Unified TWAS Z-statistic: accepts a weight vector (single method/context) or
#' a (variants x K) weight matrix, computes the per-tuple TWAS Z-score and
#' two-sided p-value, and optionally delegates cross-tuple p-value combination
#' to \code{\link{combinePValues}}.
#'
#' For each column k of \code{weights}:
#' \itemize{
#'   \item \code{stat_k = w_k^T z}
#'   \item \code{denom_k = w_k^T R w_k}
#'   \item \code{Z_k = stat_k / sqrt(denom_k)}, \code{p_k = 2 * (1 -
#'   Phi(|Z_k|))}
#' }
#' When \code{combineMethods} is non-NULL and K >= 2, the cross-tuple
#' correlation matrix \code{rho_{i,j} = covY_{i,j} / sqrt(covY_{i,i} *
#' covY_{j,j})} is constructed once and forwarded to
#' \code{combinePValues} as the \code{R} argument. When K == 1, the
#' combined p-value trivially equals the per-tuple p-value.
#'
#' The SVD path (\code{V}, \code{D}, \code{nSketch}) lets the caller avoid
#' materializing the full LD matrix: \code{covY = (V^TW * sqrt(Lambda))^T (V^TW
#' * sqrt(Lambda))} with \code{Lambda_i = D_i^2 / (nSketch - 1)}. Use the
#' \code{R} path when an LD correlation matrix is already available; use the
#' \code{X} path to compute \code{R} from a genotype matrix.
#'
#' @param weights Numeric vector of weights (single tuple) or a numeric matrix
#'   with one column per tuple (method / context). When a vector, the column
#'   name defaults to \code{"method1"}.
#' @param z Numeric vector of GWAS Z-scores aligned to the rows of
#'   \code{weights}.
#' @param R Optional LD correlation matrix.
#' @param X Optional genotype matrix used to compute \code{R} when \code{R} is
#'   missing.
#' @param V,D,nSketch SVD components of the LD sketch (right-singular vectors,
#'   singular values, panel sample size). Supplying all three selects the SVD
#'   path.
#' @param combineMethods Optional character vector of method names to forward to
#'   \code{\link{combinePValues}} for cross-tuple combination. \code{NULL}
#'   (default) skips combination.
#' @return A list with:
#' \describe{
#'   \item{Z}{A \code{K x 2} numeric matrix with columns
#'     \code{c("Z", "pval")}; rownames are the column names of
#'     \code{weights}.}
#'   \item{combined}{Output of \code{combinePValues} (or its trivial
#'     K=1 equivalent) when \code{combineMethods} is non-NULL,
#'     otherwise \code{NULL}.}
#' }
#' @seealso \code{\link{combinePValues}} for the combination method menu.
#' @importFrom stats pnorm
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:20]
#' w <- setNames(rnorm(20) * 0.1, colnames(X))
#' twasZ(weights = w, z = rnorm(20), R = cor(X))
#' @export
twasZ <- function(
    weights,
    z,
    R = NULL,
    X = NULL,
    V = NULL,
    D = NULL,
    nSketch = NULL,
    combineMethods = NULL
) {
    weights <- .twasZPrepWeights(weights, z)
    covY <- .twasZCovY(
        weights = weights,
        R = R,
        X = X,
        V = V,
        D = D,
        nSketch = nSketch
    )$covY
    ySd <- sqrt(diag(covY))
    stats <- as.numeric(crossprod(weights, as.numeric(z)))
    zVec <- stats / ySd
    pVec <- .zToPvalue(zVec)
    zMatrix <- cbind(Z = zVec, pval = pVec)
    rownames(zMatrix) <- colnames(weights)
    combined <- .twasZCombine(
        combineMethods,
        ncol(weights),
        pVec,
        zVec,
        covY,
        ySd,
        zMatrix
    )
    list(Z = zMatrix, combined = combined)
}

# Coerce a weight vector to a one-column matrix, validate the class / dims, and
# default the column names.
# @noRd
.twasZPrepWeights <- function(weights, z) {
    if (is.numeric(weights) && is.null(dim(weights))) {
        nm <- if (!is.null(names(weights))) names(weights) else NULL
        weights <- matrix(weights, ncol = 1L, dimnames = list(nm, "method1"))
    }
    if (!is.matrix(weights)) {
        abort("`weights` must be a numeric vector or a matrix.")
    }
    if (is.null(colnames(weights))) {
        colnames(weights) <- str_c("method", seq_len(ncol(weights)))
    }
    if (nrow(weights) != length(z)) {
        abort("nrow(weights) must equal length(z).")
    }
    weights
}

# Optional cross-tuple p-value combination via combinePValues (K == 1 uses the
# trivial identity form). Returns NULL when combineMethods is NULL.
# @noRd
.twasZCombine <- function(combineMethods, K, pVec, zVec, covY, ySd, zMatrix) {
    if (is.null(combineMethods)) {
        return(NULL)
    }
    combineMethods <- as.character(combineMethods)
    if (K == 1L) {
        return(.twasZCombineSingle(combineMethods, pVec, zMatrix))
    }
    sig <- covY / tcrossprod(ySd, ySd)
    rownames(sig) <- colnames(sig) <- rownames(zMatrix)
    names(pVec) <- rownames(zMatrix)
    names(zVec) <- rownames(zMatrix)
    combinePValues(
        pvals = pVec,
        zScores = zVec,
        methods = combineMethods,
        R = sig
    )
}

# K == 1 combination: every requested method trivially takes the single tuple's
# p-value, with a 1x1 aligned correlation matrix.
# @noRd
.twasZCombineSingle <- function(combineMethods, pVec, zMatrix) {
    perMethod <- set_names(
        map(combineMethods, .cipSingleMethodPval, pVec = pVec),
        combineMethods
    )
    list(
        input = list(
            nPvalsIn = 1L,
            nZScoresIn = 1L,
            nValid = 1L,
            Raligned = matrix(
                1.0,
                1L,
                1L,
                dimnames = list(rownames(zMatrix), rownames(zMatrix))
            )
        ),
        results = perMethod
    )
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# The CV `which` metric for work-list row `i` of a (study, context, trait,
# method) table. Shared by the rsq and pval columns (differ only by `which`).
# @noRd
.cipRowCvMetric <- function(i, twasWeights, qtlRows, which) {
    .cipCvMetric(
        twasWeights,
        qtlRows$qtlStudy[[i]],
        qtlRows$context[[i]],
        qtlRows$trait[[i]],
        qtlRows$method[[i]],
        which = which
    )
}

# One method's single-tuple record (K == 1: it just takes the lone p-value).
# @noRd
.cipSingleMethodPval <- function(m, pVec) {
    list(method = m, pval = as.numeric(pVec[[1L]]))
}
