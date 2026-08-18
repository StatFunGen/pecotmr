#' @title Colocalization Pipeline (coloc.bf_bf over QTL + GWAS LBF matrices)
#' @description Per-region pipeline that pairs a QTL
#'   \code{\link{QtlFineMappingResult}} with a GWAS fine-mapping result (either
#'   supplied directly as a \code{\link{GwasFineMappingResult}} or computed
#'   inline from a \code{\link{GwasSumStats}}) and runs
#'   \code{coloc::coloc.bf_bf} per (QTL tuple, GWAS tuple) pair to produce
#'   per-pair colocalization posterior probabilities PP.H0-PP.H4.
#'
#' @section Why \code{coloc.bf_bf} and not \code{coloc.susie}:
#' The prior \code{colocWrapper} (now stubbed) used
#' \code{coloc::coloc.bf_bf} on the SuSiE \code{lbf_variable} matrices
#' directly. That choice carries three behaviours that
#' \code{coloc::coloc.susie} does not expose:
#' \itemize{
#'   \item \strong{fSuSiE support}: the LBF matrix lives at a different
#'     slot for fSuSiE fits (\code{fsusie_result$lBF}) and gets stacked
#'     into a single matrix.
#'   \item \strong{Effect filtering}: \code{filterLbfCs} keeps only
#'     effects that produced a credible set; \code{filterLbfCsSecondary}
#'     keeps effects at a secondary coverage; otherwise the default
#'     filter drops effects whose prior variance is below
#'     \code{priorTol}.
#'   \item \strong{Multiple-GWAS batching}: when several GWAS
#'     fine-mapping rows fall in the same region they are merged into
#'     one combined LBF matrix per QTL pair (one \code{coloc.bf_bf}
#'     call covers them all).
#' }
#' This pipeline preserves all three.
#'
#'   GWAS input dispatch:
#'   \itemize{
#'     \item \code{gwasInput} is a \code{\link{GwasSumStats}}: GWAS
#'           fine-mapping is performed inline by
#'           \code{\link{fineMappingPipeline}} with the supplied
#'           \code{finemappingMethods} (default \code{"susie"}).
#'     \item \code{gwasInput} is a \code{\link{GwasFineMappingResult}}:
#'           used directly; no inline fine-mapping.
#'   }
#'
#' @section LD-sketch identity check: If
#'   \code{getLdSketch(qtlFineMappingResult)} is non-\code{NULL}, it must match
#'   the LD sketch on \code{gwasInput}. Mismatch is a hard error. When the QTL
#'   FMR's \code{ldSketch} is \code{NULL} (individual-level fit), the validation
#'   is skipped on the QTL side.
#'
#' @param qtlFineMappingResult A \code{\link{QtlFineMappingResult}} (required).
#' @param gwasInput Either a \code{\link{GwasSumStats}} or a
#'   \code{\link{GwasFineMappingResult}}.
#' @param filterLbfCs Logical. When \code{TRUE} (and \code{filterLbfCsSecondary}
#'   is \code{NULL}), keep only effects that produced a credible set
#'   (\code{trimmedFit$sets$cs_index}). Default \code{FALSE}.
#' @param filterLbfCsSecondary Optional secondary coverage (numeric in \eqn{(0,
#'   1)}). When supplied, run a credible-set concentration filter at this
#'   coverage level instead of \code{filterLbfCs}: each L-effect's credible set
#'   must span fewer than \code{nVariants * filterLbfCsSecondary *
#'   filterLbfCsConcentration} variants to be kept. Effects with diffuse
#'   credible sets are dropped before the LBF matrix is passed to
#'   \code{coloc::coloc.bf_bf}. Overrides \code{filterLbfCs} when set.
#' @param filterLbfCsConcentration Numeric in \eqn{(0, 1)}; the concentration
#'   factor in the cutoff above. With the default \code{0.5} a 50\% credible set
#'   is kept only if it spans fewer than 25\% of the locus's variants. Only
#'   consulted when \code{filterLbfCsSecondary} is non-NULL. Default \code{0.5}.
#' @param priorTol Prior-variance cutoff for the default filter: effects with
#'   \code{V <= priorTol} are dropped. Ignored when either \code{filterLbfCs} or
#'   \code{filterLbfCsSecondary} is in use. Default \code{1e-9}.
#' @param p1 Prior probability of QTL signal per variant. Default \code{1e-4}.
#' @param p2 Prior probability of GWAS signal per variant. Default \code{1e-4}.
#' @param p12 Prior probability of shared signal per variant. Default
#'   \code{5e-6}.
#' @param finemappingMethods Character vector forwarded to
#'   \code{\link{fineMappingPipeline}} when \code{gwasInput} is a
#'   \code{GwasSumStats}. Default \code{"susie"}.
#' @param returnGwasFineMapping Logical. When \code{TRUE}, attach the computed
#'   \code{GwasFineMappingResult} on the returned data frame as attribute
#'   \code{"gwasFineMapping"}. Default \code{FALSE}.
#' @param enrichment Optional data.frame of per-(gwasStudy, qtlStudy,
#'   qtlContext) enrichment factors with columns \code{gwasStudy},
#'   \code{qtlStudy}, \code{qtlContext}, \code{enrichment}. Output of
#'   \code{\link{qtlEnrichmentPipeline}}. When non-\code{NULL}, each pair's
#'   \code{p12} prior is scaled to \code{min(p12 * (1 + enrichment), p12Max)}
#'   (the enrichment-informed colocalization variant, "enloc"). Pairs without a
#'   matching enrichment row fall back to the baseline \code{p12} with a
#'   warning. Default \code{NULL} (baseline coloc).
#' @param p12Max Numeric scalar. Maximum value for the enrichment-adjusted
#'   \code{p12} prior. Default \code{1e-3}. Ignored when \code{enrichment =
#'   NULL}.
#' @param adjustPips Logical, default \code{TRUE}. When TRUE, before any
#'   per-pair inference the QTL and GWAS fine-mapping result collections are
#'   passed through \code{\link{adjustPips}} so each entry's PIPs are
#'   renormalized to the intersection of its variants with the union of the
#'   other side's variant IDs. This matters in two scenarios: (1) the user
#'   declined to impute missing variants in the GWAS \code{SumStats} and the QTL
#'   fine-mapping input has additional variants; (2) the GWAS fine-mapping
#'   result contains variants not present in the QTL fine-mapping result. Pass
#'   \code{FALSE} to use the FMRs as supplied.
#' @param alleleFlip Logical, default \code{TRUE}. When TRUE, align LBF columns
#'   between the QTL and GWAS by (chrom, pos) with ref/alt swaps recognized (LBF
#'   is coding-invariant, so no sign change is needed); when FALSE, match on
#'   exact alleles only, so a ref/alt swap is treated as a distinct variant.
#' @param ... Additional arguments forwarded to \code{coloc::coloc.bf_bf}.
#' @return A data frame with one row per (QTL tuple, GWAS tuple, credible-set
#'   pair) combination. Identity columns: \code{study}, \code{context},
#'   \code{trait}, \code{method}, \code{gwasStudy}, \code{gwasMethod}, plus the
#'   standard coloc fields (\code{idx1}, \code{idx2}, \code{nSnps},
#'   \code{PP.H0.abf} \ldots \code{PP.H4.abf}). When \code{enrichment} is
#'   supplied, two additional columns \code{enrichment} and \code{p12Used}
#'   report the per-pair factor and the prior actually passed to
#'   \code{coloc::coloc.bf_bf}.
#' @examples
#' data(qtlFineMappingLbfExample)
#' data(gwasFineMappingLbfExample)
#' colocPipeline(qtlFineMappingLbfExample,
#'   gwasInput = gwasFineMappingLbfExample)
#' @export
colocPipeline <- function(
    qtlFineMappingResult,
    gwasInput,
    filterLbfCs = FALSE,
    filterLbfCsSecondary = NULL,
    filterLbfCsConcentration = 0.5,
    priorTol = 1e-9,
    p1 = 1e-4,
    p2 = 1e-4,
    p12 = 5e-6,
    finemappingMethods = "susie",
    returnGwasFineMapping = FALSE,
    enrichment = NULL,
    p12Max = 1e-3,
    adjustPips = TRUE,
    alleleFlip = TRUE,
    ...
) {
    p <- as.list(environment())
    p$dots <- list(...)
    p$useEnrichment <- !is.null(enrichment)
    .colocValidateInputs(p)
    p$gwasFmr <- .colocResolveGwasFmr(gwasInput, finemappingMethods)
    .colocRequireMatchingLdSketches(
        getLdSketch(qtlFineMappingResult),
        getLdSketch(p$gwasFmr)
    )
    p <- .colocMaybeAdjustPips(p)
    # Pre-extract per-GWAS-tuple LBF matrices: group the GWAS FMR by study,
    # stack each study's LBF rows, and store per-(study, method) batched
    # matrices (reproduces the legacy row-wise combine per xQTL).
    p$gwasLbfByPair <- .colocPreextractGwasLbf(
        p$gwasFmr,
        filterLbfCs,
        filterLbfCsSecondary,
        filterLbfCsConcentration,
        priorTol
    )
    if (length(p$gwasLbfByPair) == 0L) {
        return(.colocEarlyReturn(p))
    }
    results <- list_flatten(map(
        seq_len(nrow(p$qtlFineMappingResult)),
        .colocScoreQtlTuple,
        p = p
    ))
    .colocFinalize(results, p)
}

# Validate the enrichment table (when supplied), the coloc package, and the
# input object classes.
# @noRd
.colocValidateInputs <- function(p) {
    if (p$useEnrichment) {
        .colocValidateEnrichment(p$enrichment)
    }
    if (!requireNamespace("coloc", quietly = TRUE)) {
        # nocov start
        msg <- glue(
            "Package 'coloc' is required for colocPipeline. ",
            "Install with: install.packages('coloc')."
        )
        abort(msg)
        # nocov end
    }
    if (!methods::is(p$qtlFineMappingResult, "QtlFineMappingResult")) {
        msg <- glue(
            "`qtlFineMappingResult` must be a QtlFineMappingResult ",
            "(got class '{class(p$qtlFineMappingResult)[[1L]]}')."
        )
        abort(msg)
    }
    if (
        !methods::is(p$gwasInput, "GwasSumStats") &&
            !methods::is(p$gwasInput, "GwasFineMappingResult")
    ) {
        msg <- glue(
            "`gwasInput` must be a GwasSumStats or a GwasFineMappingResult ",
            "(got class '{class(p$gwasInput)[[1L]]}')."
        )
        abort(msg)
    }
    invisible(NULL)
}

# The enrichment table must be a data.frame carrying the required id + value
# columns.
# @noRd
.colocValidateEnrichment <- function(enrichment) {
    if (!is.data.frame(enrichment)) {
        msg <- glue(
            "`enrichment` must be a data.frame with at least gwasStudy, ",
            "qtlStudy, qtlContext, enrichment columns (output of ",
            "qtlEnrichmentPipeline)."
        )
        abort(msg)
    }
    required <- c("gwasStudy", "qtlStudy", "qtlContext", "enrichment")
    missingCols <- setdiff(required, colnames(enrichment))
    if (length(missingCols) > 0L) {
        msg <- glue(
            "`enrichment` is missing column(s): ",
            "{str_flatten(missingCols, ', ')}"
        )
        abort(msg)
    }
    invisible(NULL)
}

# Resolve the GWAS side to a GwasFineMappingResult (fine-map QC'd sumstats when
# a GwasSumStats is passed).
# @noRd
.colocResolveGwasFmr <- function(gwasInput, finemappingMethods) {
    if (methods::is(gwasInput, "GwasFineMappingResult")) {
        return(gwasInput)
    }
    if (length(getQcInfo(gwasInput)) == 0L) {
        msg <- glue(
            "colocPipeline: gwasInput (GwasSumStats) has no QC record. ",
            "Call summaryStatsQc() first."
        )
        abort(msg)
    }
    fineMappingPipeline(gwasInput, methods = finemappingMethods)
}

# Optional PIP renormalization: adjust each side's entry PIPs to the
# intersection of its own variants with the union of the other side's.
# @noRd
.colocMaybeAdjustPips <- function(p) {
    if (!isTRUE(p$adjustPips)) {
        return(p)
    }
    qtlVids <- unique(unlist(map(
        p$qtlFineMappingResult$entry,
        .colocEntryVariantIds
    )))
    gwasVids <- unique(unlist(map(p$gwasFmr$entry, .colocEntryVariantIds)))
    if (length(qtlVids) > 0L && length(gwasVids) > 0L) {
        p$qtlFineMappingResult <- adjustPips(p$qtlFineMappingResult, gwasVids)
        p$gwasFmr <- adjustPips(p$gwasFmr, qtlVids)
    }
    p
}

# Empty-result early return (attaching the GWAS fine-mapping when requested).
# @noRd
.colocEarlyReturn <- function(p) {
    out <- .colocEmptyResult(enriched = p$useEnrichment)
    if (p$returnGwasFineMapping && methods::is(p$gwasInput, "GwasSumStats")) {
        attr(out, "gwasFineMapping") <- p$gwasFmr
    }
    out
}

# Score one QTL tuple against every pre-extracted GWAS pair -> summary rows.
# @noRd
.colocScoreQtlTuple <- function(qi, p) {
    q <- .colocQtlTupleInfo(qi, p)
    qLbfInfo <- .colocExtractLbfFromEntry(
        q$entry,
        p$filterLbfCs,
        p$filterLbfCsSecondary,
        p$filterLbfCsConcentration,
        p$priorTol,
        label = q$label
    )
    if (is.null(qLbfInfo)) {
        return(list())
    }
    compact(map(
        names(p$gwasLbfByPair),
        .colocScorePairAt,
        qLbfInfo = qLbfInfo,
        p = p,
        q = q
    ))
}

# Identity + entry + log label for a QTL tuple.
# @noRd
.colocQtlTupleInfo <- function(qi, p) {
    fmr <- p$qtlFineMappingResult
    study <- as.character(fmr$study)[[qi]]
    context <- as.character(fmr$context)[[qi]]
    trait <- as.character(fmr$trait)[[qi]]
    method <- as.character(fmr$method)[[qi]]
    list(
        study = study,
        context = context,
        trait = trait,
        method = method,
        entry = fmr$entry[[qi]],
        label = glue(
            "QTL (study='{study}', context='{context}', ",
            "trait='{trait}', method='{method}')"
        )
    )
}

# Score one (QTL, GWAS) pair via coloc.bf_bf -> a summary row, or NULL when the
# variants don't align or coloc fails / returns no summary.
# @noRd
.colocScorePair <- function(qLbf, gInfo, q, p) {
    # Align variants between the QTL and GWAS LBF matrices by (chrom, pos,
    # allele) tuple via matchVariants (see .colocAlignLbf).
    aligned <- .colocAlignLbf(qLbf, gInfo$lbf, alleleFlip = p$alleleFlip)
    if (is.null(aligned)) {
        return(NULL)
    }
    p12Info <- .colocResolveP12(p, gInfo$study, q$study, q$context)
    pairRes <- .colocRunPair(aligned, p, p12Info$p12Used, q, gInfo)
    if (is.null(pairRes) || is.null(pairRes$summary)) {
        return(NULL)
    }
    .colocSummaryRow(pairRes, q, gInfo, p, p12Info)
}

# Enrichment-informed p12 (per-(gwasStudy, qtlStudy, qtlContext) scaling capped
# at p12Max; baseline p12 with no enrichment table / no matching row).
# @noRd
.colocResolveP12 <- function(p, gwasStudy, qStudy, qContext) {
    if (!p$useEnrichment) {
        return(list(enRow = NA_real_, p12Used = p$p12))
    }
    enRow <- .colocLookupEnrichment(p$enrichment, gwasStudy, qStudy, qContext)
    if (is.na(enRow)) {
        msg <- glue(
            "colocPipeline: no enrichment entry for ",
            "(gwasStudy='{gwasStudy}', qtlStudy='{qStudy}', ",
            "qtlContext='{qContext}'); using baseline p12."
        )
        warn(msg)
        enRow <- 0
    }
    list(enRow = enRow, p12Used = min(p$p12 * (1 + enRow), p$p12Max))
}

# Run coloc.bf_bf for an aligned pair, warning + NULL on failure.
# @noRd
.colocRunPair <- function(aligned, p, p12Used, q, gInfo) {
    colocArgs <- c(
        list(
            aligned$qtl,
            aligned$gwas,
            p1 = p$p1,
            p2 = p$p2,
            p12 = p12Used
        ),
        p$dots
    )
    tryCatch(
        exec(coloc::coloc.bf_bf, !!!colocArgs),
        error = function(e) {
            msg <- glue(
                "colocPipeline: coloc.bf_bf failed for QTL ",
                "(study='{q$study}', context='{q$context}', ",
                "trait='{q$trait}', method='{q$method}') x GWAS ",
                "(study='{gInfo$study}', method='{gInfo$method}'): ",
                "{conditionMessage(e)}"
            )
            warn(msg)
            NULL
        }
    )
}

# Build a coloc summary row stamped with the QTL / GWAS identity + enrichment.
# @noRd
.colocSummaryRow <- function(pairRes, q, gInfo, p, p12Info) {
    sm <- as.data.frame(pairRes$summary, stringsAsFactors = FALSE)
    sm$study <- q$study
    sm$context <- q$context
    sm$trait <- q$trait
    sm$method <- q$method
    sm$gwasStudy <- gInfo$study
    sm$gwasMethod <- gInfo$method
    if (p$useEnrichment) {
        sm$enrichment <- p12Info$enRow
        sm$p12Used <- p12Info$p12Used
    }
    sm
}

# Assemble the result table + attach the GWAS fine-mapping when requested.
# @noRd
.colocFinalize <- function(results, p) {
    out <- .colocAssemble(results, p$useEnrichment)
    if (p$returnGwasFineMapping && methods::is(p$gwasInput, "GwasSumStats")) {
        attr(out, "gwasFineMapping") <- p$gwasFmr
    }
    out
}

# Row-bind + column-order the per-pair summary rows (empty result when none).
# @noRd
.colocAssemble <- function(results, useEnrichment) {
    if (length(results) == 0L) {
        return(.colocEmptyResult(enriched = useEnrichment))
    }
    out <- bind_rows(map(results, .colocStandardiseRow))
    idCols <- c(
        "study",
        "context",
        "trait",
        "method",
        "gwasStudy",
        "gwasMethod"
    )
    if (useEnrichment) {
        idCols <- c(idCols, "enrichment", "p12Used")
    }
    select(out, all_of(idCols), everything())
}

# =============================================================================
# Internal helpers
# =============================================================================

# LD-sketch identity check. Thin wrapper over the shared
# `.requireMatchingLdSketches` helper (R/ld.R). Shared with
# qtlEnrichmentPipeline.
# @noRd
.colocRequireMatchingLdSketches <- function(qtlLd, gwasLd) {
    .requireMatchingLdSketches(qtlLd, gwasLd, pipelineName = "colocPipeline")
}

# SuSiE credible-set concentration filter. Given a trimmed SuSiE fit
# and a coverage level, return the L-effect indices whose credible set
# is "narrow enough" to be informative: |CS| < nVariants * coverage *
# concentration. With concentration = 0.5 a 50% CS is kept only if it
# spans fewer than 25% of the locus variants -- this prunes diffuse
# signals before they reach coloc.bf_bf.
#
# Returns an integer vector of kept effect indices (empty when nothing
# survives), or errors when susieR is unavailable.
# @noRd
#' @importFrom susieR susie_get_cs
#' @importFrom purrr map_lgl
.colocFilterCsByConcentration <- function(
    fit,
    coverage = 0.5,
    concentration = 0.5
) {
    fit$V <- NULL # disable V-based filtering inside susie_get_cs
    csList <- susie_get_cs(fit, coverage = coverage, dedup = FALSE)
    totalVariants <- ncol(fit$alpha)
    maxSize <- totalVariants * coverage * concentration
    keep <- map_lgl(csList$cs, .colocCsUnderMax, maxSize = maxSize)
    as.numeric(str_remove_all(names(which(keep)), "L"))
}

# Extract an LBF matrix (effects x variants) from a FineMappingEntry,
# applying the same filtering knobs as the legacy .extractLbfMatrix:
#   - filterLbfCs (CS-only)
#   - filterLbfCsSecondary (secondary coverage CS, with concentration cutoff)
#   - priorTol drop on V (default)
# Handles the fSuSiE shape (where the LBF lives at a different slot).
# Returns list(lbf = <matrix>, variantIds = <character>) or NULL when
# the entry has no usable LBF matrix.
# @noRd
.colocExtractLbfFromEntry <- function(
    entry,
    filterLbfCs,
    filterLbfCsSecondary,
    filterLbfCsConcentration,
    priorTol,
    label = "entry"
) {
    fit <- getSusieFit(entry)
    if (is.null(fit)) {
        msg <- glue("colocPipeline: {label} has no trimmedFit; skipping.")
        warn(msg)
        return(NULL)
    }
    lbfMatrix <- .colocLbfMatrix(fit, label)
    if (is.null(lbfMatrix)) {
        return(NULL)
    }
    lbfMatrix <- .colocFilterLbfRows(
        lbfMatrix,
        fit,
        filterLbfCs,
        filterLbfCsSecondary,
        filterLbfCsConcentration,
        priorTol
    )
    if (nrow(lbfMatrix) == 0L) {
        return(NULL)
    }
    lbfMatrix <- .colocAssignLbfColnames(lbfMatrix, entry)
    if (ncol(lbfMatrix) == 0L) {
        return(NULL)
    }
    list(lbf = lbfMatrix, variantIds = colnames(lbfMatrix))
}

# Extract the (effects x variants) LBF matrix from a fit: susie lbf_variable, or
# a stacked fSuSiE lBF (direct or nested), or NULL (with a warning) when absent
# / empty.
# @noRd
.colocLbfMatrix <- function(fit, label) {
    lbfMatrix <- if (!is.null(fit$lbf_variable)) {
        as.matrix(fit$lbf_variable)
    } else if (
        !is.null(fit$fsusie_result) &&
            is.list(fit$fsusie_result$lBF)
    ) {
        # fSuSiE path: stack per-trait lBF lists into a single matrix.
        lbfList <- fit$fsusie_result$lBF
        exec(rbind, !!!lbfList)
    } else if (
        is.list(fit) &&
            length(fit) >= 1L &&
            !is.null(fit[[1L]]$fsusie_result$lBF)
    ) {
        lbfList <- fit[[1L]]$fsusie_result$lBF
        exec(rbind, !!!lbfList)
    } else {
        msg <- glue(
            "colocPipeline: {label} trimmedFit has no lbf_variable / fsusie ",
            "lBF; skipping."
        )
        warn(msg)
        return(NULL)
    }
    if (is.null(lbfMatrix) || nrow(lbfMatrix) == 0L) {
        msg <- glue("colocPipeline: {label} LBF matrix is empty.")
        warn(msg)
        return(NULL)
    }
    lbfMatrix
}

# Row (effect) filtering, in the original priority order: primary CS index,
# else secondary-CS-by-concentration, else prior-variance threshold.
# @noRd
.colocFilterLbfRows <- function(
    lbfMatrix,
    fit,
    filterLbfCs,
    filterLbfCsSecondary,
    filterLbfCsConcentration,
    priorTol
) {
    if (isTRUE(filterLbfCs) && is.null(filterLbfCsSecondary)) {
        csIdx <- fit$sets$cs_index
        if (!is.null(csIdx) && length(csIdx) > 0L) {
            lbfMatrix <- lbfMatrix[csIdx, , drop = FALSE]
        }
    } else if (!is.null(filterLbfCsSecondary)) {
        secIdx <- tryCatch(
            .colocFilterCsByConcentration(
                fit,
                coverage = filterLbfCsSecondary,
                concentration = filterLbfCsConcentration
            ),
            error = function(e) NULL
        )
        if (!is.null(secIdx) && length(secIdx) > 0L) {
            lbfMatrix <- lbfMatrix[secIdx, , drop = FALSE]
        }
    } else if (!is.null(fit$V)) {
        lbfMatrix <- lbfMatrix[fit$V > priorTol, , drop = FALSE]
    }
    lbfMatrix
}

# Assign variant-id column names (fit-provided, else the entry's variantIds
# slot) and drop columns with an NA id.
# @noRd
.colocAssignLbfColnames <- function(lbfMatrix, entry) {
    if (is.null(colnames(lbfMatrix)) || any(is.na(colnames(lbfMatrix)))) {
        vids <- getVariantIds(entry)
        if (length(vids) == ncol(lbfMatrix)) {
            colnames(lbfMatrix) <- vids
        }
    }
    lbfMatrix[, !is.na(colnames(lbfMatrix)), drop = FALSE]
}

# Build a per-GWAS-tuple LBF matrix list, keyed by "study|method".
# Within each key we stack multiple FMR rows row-wise (the legacy
# "combined GWAS LBF" pattern), drop NA columns, and replace NAs with
# 0 so a fresh QTL pairing always lands on the same coordinate frame.
# @noRd
.colocPreextractGwasLbf <- function(
    gwasFmr,
    filterLbfCs,
    filterLbfCsSecondary,
    filterLbfCsConcentration,
    priorTol
) {
    groupKey <- str_c(
        as.character(gwasFmr$study),
        as.character(gwasFmr$method),
        sep = "||"
    )
    groups <- split(seq_len(nrow(gwasFmr)), groupKey)
    out <- list()
    for (gkey in names(groups)) {
        rows <- groups[[gkey]]
        pieces <- list()
        for (ri in rows) {
            info <- .colocExtractLbfFromEntry(
                gwasFmr$entry[[ri]],
                filterLbfCs,
                filterLbfCsSecondary,
                filterLbfCsConcentration,
                priorTol,
                label = glue(
                    "GWAS (study='{as.character(gwasFmr$study)[[ri]]}', ",
                    "method='{as.character(gwasFmr$method)[[ri]]}', row={ri})"
                )
            )
            if (!is.null(info)) pieces[[length(pieces) + 1L]] <- info$lbf
        }
        if (length(pieces) == 0L) {
            next
        }
        combined <- .colocRbindLbf(pieces)
        parts <- str_split(gkey, "\\|\\|")[[1L]]
        out[[gkey]] <- list(
            lbf = combined,
            study = parts[[1L]],
            method = if (length(parts) >= 2L) {
                parts[[2L]]
            } else {
                NA_character_
            }
        )
    }
    out
}

# Row-bind a list of LBF matrices with potentially-different column
# sets. Uses the union of columns; missing cells fill with 0 (matching
# the legacy `replace_na(., 0)` step inside colocWrapper).
# @noRd
.colocRbindLbf <- function(mats) {
    allCols <- unique(unlist(map(mats, colnames)))
    padded <- map(mats, .colocPadCols, allCols = allCols)
    exec(rbind, !!!padded)
}

# Align column names between a QTL and a (combined) GWAS LBF matrix by
# (chrom, pos, allele) tuple (via matchVariants), returning both restricted to
# the common variants under one shared id so coloc.bf_bf can align them.
# @noRd
.colocAlignLbf <- function(qtlLbf, gwasLbf, alleleFlip = TRUE) {
    qids <- colnames(qtlLbf)
    gids <- colnames(gwasLbf)
    # Match LBF columns by (chrom, pos, allele) tuple rather than by raw id
    # string, so chr-prefix / separator / allele-order differences resolve. LBF
    # is allele-coding-invariant, so this is pure identity alignment (no sign);
    # with alleleFlip = FALSE, ref/alt-swapped columns are not treated as
    # shared.
    m <- matchVariants(qids, gids, allowFlip = alleleFlip)
    if (length(m$idxA) == 0L) {
        return(NULL)
    }
    # Relabel both matrices to one shared id so coloc.bf_bf sees identical
    # names.
    sharedIds <- gids[m$idxB]
    qSub <- qtlLbf[, m$idxA, drop = FALSE]
    gSub <- gwasLbf[, m$idxB, drop = FALSE]
    colnames(qSub) <- sharedIds
    colnames(gSub) <- sharedIds
    list(qtl = qSub, gwas = gSub)
}

# A blank result data frame for the no-pair case so callers downstream
# do not have to special-case a NULL return. When `enriched = TRUE` the
# enrichment + p12Used columns are appended (the enloc-mode schema).
# @noRd
.colocEmptyResult <- function(enriched = FALSE) {
    base <- tibble(
        study = character(0),
        context = character(0),
        trait = character(0),
        method = character(0),
        gwasStudy = character(0),
        gwasMethod = character(0),
        idx1 = integer(0),
        idx2 = integer(0),
        nSnps = integer(0),
        PP.H0.abf = numeric(0),
        PP.H1.abf = numeric(0),
        PP.H2.abf = numeric(0),
        PP.H3.abf = numeric(0),
        PP.H4.abf = numeric(0)
    )
    if (enriched) {
        base$enrichment <- numeric(0)
        base$p12Used <- numeric(0)
    }
    base
}

# Look up the enrichment factor for a (gwasStudy, qtlStudy, qtlContext)
# triple in the user-supplied enrichment table. Returns NA when the
# triple is not present; the caller falls back to the baseline p12 and
# emits a warning.
# @noRd
.colocLookupEnrichment <- function(
    enrichment,
    gwasStudy,
    qtlStudy,
    qtlContext
) {
    idx <- which(
        as.character(enrichment$gwasStudy) == gwasStudy &
            as.character(enrichment$qtlStudy) == qtlStudy &
            as.character(enrichment$qtlContext) == qtlContext
    )
    if (length(idx) == 0L) {
        return(NA_real_)
    }
    as.numeric(enrichment$enrichment[[idx[[1L]]]])
}

# Ensure each row data.frame from coloc.bf_bf carries the standard PP
# columns even when the underlying call produced a slightly different
# shape.
# @noRd
.colocStandardiseRow <- function(sm) {
    for (col in c(
        "idx1",
        "idx2",
        "nSnps",
        "PP.H0.abf",
        "PP.H1.abf",
        "PP.H2.abf",
        "PP.H3.abf",
        "PP.H4.abf"
    )) {
        if (!is_in(col, colnames(sm))) sm[[col]] <- NA
    }
    sm
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# The variant ids of one fine-mapping entry (S4 slot; not pluckable by name).
# @noRd
.colocEntryVariantIds <- function(e) {
    e@variantIds
}

# Score the QTL LBF against GWAS pair `gKey` -> a summary row (or NULL).
# @noRd
.colocScorePairAt <- function(gKey, qLbfInfo, p, q) {
    .colocScorePair(qLbfInfo$lbf, p$gwasLbfByPair[[gKey]], q, p)
}

# TRUE when a credible set has fewer than `maxSize` variants.
# @noRd
.colocCsUnderMax <- function(x, maxSize) {
    length(x) < maxSize
}

# Pad one LBF matrix to the union column set `allCols` (missing cells -> 0).
# @noRd
.colocPadCols <- function(m, allCols) {
    miss <- setdiff(allCols, colnames(m))
    if (length(miss) > 0L) {
        pad <- matrix(
            0,
            nrow = nrow(m),
            ncol = length(miss),
            dimnames = list(NULL, miss)
        )
        m <- cbind(m, pad)
    }
    m[, allCols, drop = FALSE]
}
