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
#'   reconciled with \code{\link{intersectVariants}}, so both sides are scored
#'   on the variants they share and their PIPs are renormalized to that shared
#'   set. This matters in two scenarios: (1) the user declined to impute
#'   missing variants in the GWAS \code{SumStats} and the QTL fine-mapping
#'   input has additional variants; (2) the GWAS fine-mapping result contains
#'   variants not present in the QTL fine-mapping result. Pass \code{FALSE} to
#'   use the FMRs as supplied.
#'
#'   How much of each effect's posterior survived that restriction is reported
#'   per result row as \code{qtlRetainedMass} / \code{gwasRetainedMass}; a low
#'   value means the effect was largely built on variants the other side does
#'   not carry, so its colocalization evidence rests on little retained signal.
#'   Both are \code{NA} when \code{adjustPips = FALSE}, since nothing was
#'   measured.
#' @param alleleFlip Logical, default \code{TRUE}. When TRUE, align LBF columns
#'   between the QTL and GWAS by (chrom, pos) with ref/alt swaps recognized (LBF
#'   is coding-invariant, so no sign change is needed); when FALSE, match on
#'   exact alleles only, so a ref/alt swap is treated as a distinct variant.
#' @param ... Additional arguments forwarded to \code{coloc::coloc.bf_bf}.
#' @return A \code{\linkS4class{ColocResult}}: one element per tested (QTL
#'   credible set, GWAS credible set, block) pair, holding that pair's aligned
#'   variants with their \code{SNP.PP.H4}. Pair-level metadata carries the
#'   identity columns (\code{study}, \code{context}, \code{trait},
#'   \code{method}, \code{gwasStudy}, \code{gwasMethod}), the block and stable
#'   credible-set ids (\code{blockId}, \code{qtlCs}, \code{gwasCs}), the
#'   standard coloc fields (\code{idx1}, \code{idx2}, \code{nSnps},
#'   \code{hit1}, \code{hit2}, \code{PP.H0.abf} \ldots \code{PP.H4.abf}) and
#'   the reconciliation diagnostics \code{qtlRetainedMass} /
#'   \code{gwasRetainedMass}. When \code{enrichment} is supplied, two
#'   additional columns \code{enrichment} and \code{p12Used} report the
#'   per-pair factor and the prior actually passed to
#'   \code{coloc::coloc.bf_bf}.
#'
#'   Project it with \code{\link{getColocPairs}} (the flat table this
#'   pipeline used to return, also available as \code{as.data.frame}),
#'   \code{\link{getColocVariants}}, \code{\link{getColocCredibleSets}} or
#'   \code{\link{getColocGenes}}.
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

# Optional PIP renormalization, via the symmetric reconciliation verb.
#
# Coloc is the symmetric case: both sides have to end up scored on the SAME
# variant set, or a pair's two halves are not comparable. That is what
# intersectVariants() gives, and it replaces the older hand-rolled pass which
# adjusted each side to the UNION of the other's variants -- a union is not an
# intersection, so the two sides could still end up on different sets.
# @noRd
.colocMaybeAdjustPips <- function(p) {
    if (!isTRUE(p$adjustPips)) {
        return(p)
    }
    if (nrow(p$qtlFineMappingResult) == 0L || nrow(p$gwasFmr) == 0L) {
        return(p)
    }
    both <- intersectVariants(p$qtlFineMappingResult, p$gwasFmr)
    p$qtlFineMappingResult <- both$x
    p$gwasFmr <- both$y
    p
}

# The LD reference the result carries forward, so getColocCredibleSets() can
# recompute purity (section 3.7) without being handed a sketch separately. The
# QTL and GWAS sketches are already required to match by
# .colocRequireMatchingLdSketches, so either one identifies the panel.
# @noRd
.colocLdSketch <- function(p) {
    getLdSketch(p$qtlFineMappingResult)
}

# Empty-result early return (attaching the GWAS fine-mapping when requested).
# @noRd
.colocEarlyReturn <- function(p) {
    out <- .colocEmptyResult(
        enriched = p$useEnrichment,
        ldSketch = .colocLdSketch(p)
    )
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
        q$parts,
        p$filterLbfCs,
        p$filterLbfCsSecondary,
        p$filterLbfCsConcentration,
        p$priorTol,
        label = q$label
    )
    if (is.null(qLbfInfo)) {
        return(list())
    }
    q$retainedMass <- qLbfInfo$retainedMass
    q$effect <- qLbfInfo$effect
    compact(map(
        names(p$gwasLbfByPair),
        .colocScorePairAt,
        qLbfInfo = qLbfInfo,
        p = p,
        q = q
    ))
}

# Identity + row payload + log label for a QTL tuple.
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
        parts = .fmrRowParts(fmr, qi),
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
    rows <- .colocSummaryRow(pairRes, q, gInfo, p, p12Info)
    # $results is the per-variant layer that process_coloc_results() used to
    # consume and that this pipeline silently dropped. It is pivoted here, the
    # only place that knows which results column belongs to which summary row.
    list(
        rows = rows,
        variants = .crPivotColocResults(pairRes$results, nrow(rows))
    )
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

# Build a coloc summary row carrying the QTL / GWAS identity + enrichment.
# @noRd
.colocSummaryRow <- function(pairRes, q, gInfo, p, p12Info) {
    sm <- as.data.frame(pairRes$summary, stringsAsFactors = FALSE)
    sm <- .colocRenameNsnps(sm)
    sm$study <- q$study
    sm$context <- q$context
    sm$trait <- q$trait
    sm$method <- q$method
    sm$gwasStudy <- gInfo$study
    sm$gwasMethod <- gInfo$method
    # idx1 / idx2 index the LBF rows handed to coloc.bf_bf, which is exactly
    # what retainedMass runs parallel to -- so the mass reported here is the
    # mass of the two effects this row actually scores, not a per-entry
    # average.
    sm$qtlRetainedMass <- .colocPickAt(
        q$retainedMass,
        sm[["idx1"]],
        nrow(sm)
    )
    sm$gwasRetainedMass <- .colocPickAt(
        gInfo$retainedMass,
        sm[["idx2"]],
        nrow(sm)
    )
    # coloc's idx1 / idx2 number the rows of THIS call, so they are not
    # comparable across blocks. The fit's own effect indices are, and they are
    # what the credible-set and gene views group on.
    sm$qtlCs <- .colocPickAt(
        q$effect,
        sm[["idx1"]],
        nrow(sm),
        fill = NA_integer_
    )
    sm$gwasCs <- .colocPickAt(
        gInfo$effect,
        sm[["idx2"]],
        nrow(sm),
        fill = NA_integer_
    )
    sm$blockId <- gInfo$blockId %||% NA_character_
    if (p$useEnrichment) {
        sm$enrichment <- p12Info$enRow
        sm$p12Used <- p12Info$p12Used
    }
    sm
}

# Assemble the result table + attach the GWAS fine-mapping when requested.
# @noRd
.colocFinalize <- function(results, p) {
    out <- .colocAssemble(results, p$useEnrichment, .colocLdSketch(p))
    if (p$returnGwasFineMapping && methods::is(p$gwasInput, "GwasSumStats")) {
        attr(out, "gwasFineMapping") <- p$gwasFmr
    }
    out
}

# Row-bind + column-order the per-pair summary rows (empty result when none).
# @noRd
.colocAssemble <- function(results, useEnrichment, ldSketch = NULL) {
    if (length(results) == 0L) {
        return(.colocEmptyResult(enriched = useEnrichment, ldSketch = ldSketch))
    }
    rows <- bind_rows(map(map(results, "rows"), .colocStandardiseRow))
    variants <- unlist(map(results, "variants"), recursive = FALSE)
    ColocResult(.colocOrderColumns(rows, useEnrichment), variants, ldSketch)
}

# Identity columns first, then everything coloc produced.
# @noRd
.colocOrderColumns <- function(out, useEnrichment) {
    idCols <- c(
        "study",
        "context",
        "trait",
        "method",
        "gwasStudy",
        "gwasMethod",
        "blockId",
        "qtlCs",
        "gwasCs"
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

# Extract an LBF matrix (effects x variants) from a FineMappingRow,
# applying the same filtering knobs as the legacy .extractLbfMatrix:
#   - filterLbfCs (CS-only)
#   - filterLbfCsSecondary (secondary coverage CS, with concentration cutoff)
#   - priorTol drop on V (default)
# Handles the fSuSiE shape (where the LBF lives at a different slot).
# Returns list(lbf = <matrix>, variantIds = <character>) or NULL when
# the entry has no usable LBF matrix.
# @noRd
.colocExtractLbfFromEntry <- function(
    parts,
    filterLbfCs,
    filterLbfCsSecondary,
    filterLbfCsConcentration,
    priorTol,
    label = "entry"
) {
    fit <- getSusieFit(parts)
    if (is.null(fit)) {
        msg <- glue("colocPipeline: {label} has no trimmedFit; skipping.")
        warn(msg)
        return(NULL)
    }
    lbfMatrix <- .colocLbfMatrix(fit, label)
    if (is.null(lbfMatrix)) {
        return(NULL)
    }
    mass <- .colocEffectRetainedMass(fit, nrow(lbfMatrix))
    keep <- .colocSelectLbfRows(
        lbfMatrix,
        fit,
        filterLbfCs,
        filterLbfCsSecondary,
        filterLbfCsConcentration,
        priorTol
    )
    lbfMatrix <- lbfMatrix[keep, , drop = FALSE]
    mass <- mass[keep]
    if (nrow(lbfMatrix) == 0L) {
        return(NULL)
    }
    lbfMatrix <- .colocAssignLbfColnames(lbfMatrix, parts)
    if (ncol(lbfMatrix) == 0L) {
        return(NULL)
    }
    list(
        lbf = lbfMatrix,
        variantIds = colnames(lbfMatrix),
        retainedMass = mass,
        # The fit's own effect indices, which -- unlike coloc's idx1 / idx2 --
        # are stable across calls and blocks, so they are what identifies a
        # credible set in the result.
        effect = keep
    )
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

# Row (effect) selection, in the original priority order: primary CS index,
# else secondary-CS-by-concentration, else prior-variance threshold.
#
# Returns the surviving row INDICES rather than the subsetted matrix. The
# caller has a second vector running parallel to those rows (the per-effect
# retained mass), and handing back a selector keeps the two narrowed by one
# decision instead of two hand-synchronised subsets.
#
# `fit` is read with `[[`, never `$`: `$` on a list falls back to prefix
# matching when the exact name is absent, so a fit carrying `sets_secondary`
# but no `sets` would silently filter on the wrong element.
# @noRd
.colocSelectLbfRows <- function(
    lbfMatrix,
    fit,
    filterLbfCs,
    filterLbfCsSecondary,
    filterLbfCsConcentration,
    priorTol
) {
    allRows <- seq_len(nrow(lbfMatrix))
    if (isTRUE(filterLbfCs) && is.null(filterLbfCsSecondary)) {
        csIdx <- fit[["sets"]][["cs_index"]]
        if (!is.null(csIdx) && length(csIdx) > 0L) {
            return(csIdx)
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
            return(secIdx)
        }
    } else if (!is.null(fit[["V"]])) {
        return(allRows[fit[["V"]] > priorTol])
    }
    allRows
}

# Per-effect retained posterior mass, parallel to the LBF rows.
#
# `retained_mass` is written by adjustPips() at reconciliation time and records
# how much of each effect's alpha survived the shared-variant restriction. It
# is absent when no reconciliation ran (adjustPips = FALSE), which is reported
# as NA -- "not measured", distinct from a measured 0.
# @noRd
.colocEffectRetainedMass <- function(fit, nEffects) {
    mass <- fit[["retained_mass"]]
    if (is.null(mass) || length(mass) != nEffects) {
        return(rep(NA_real_, nEffects))
    }
    as.numeric(mass)
}

# Assign variant-id column names (fit-provided, else the row's rendered
# variant ids) and drop columns with an NA id.
# @noRd
.colocAssignLbfColnames <- function(lbfMatrix, parts) {
    if (is.null(colnames(lbfMatrix)) || any(is.na(colnames(lbfMatrix)))) {
        vids <- .fmrPartsVariantIds(parts)
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
    keys <- str_c(
        as.character(gwasFmr$study),
        as.character(gwasFmr$method),
        as.character(.colocGwasBlockIds(gwasFmr)),
        sep = "||"
    )
    out <- list()
    for (ri in seq_len(nrow(gwasFmr))) {
        parts <- str_split(keys[[ri]], "\\|\\|")[[1L]]
        info <- .colocExtractLbfFromEntry(
            .fmrRowParts(gwasFmr, ri),
            filterLbfCs,
            filterLbfCsSecondary,
            filterLbfCsConcentration,
            priorTol,
            label = glue(
                "GWAS (study='{parts[[1L]]}', method='{parts[[2L]]}', ",
                "block='{parts[[3L]]}')"
            )
        )
        if (is.null(info)) {
            next
        }
        out[[keys[[ri]]]] <- list(
            lbf = info$lbf,
            retainedMass = info$retainedMass,
            effect = info$effect,
            study = parts[[1L]],
            method = parts[[2L]],
            blockId = parts[[3L]]
        )
    }
    out
}

# The LD block each GWAS fine-mapping row was computed on. The element's own
# range is the identity; an explicit `blockId` (which keys the external block
# manifest, and so carries the true block BOUNDARIES rather than the span of
# the variants that survived) is preferred when present.
# @noRd
.colocGwasBlockIds <- function(gwasFmr) {
    md <- mcols(gwasFmr, use.names = FALSE)
    if (is_in("blockId", colnames(md))) {
        return(as.character(md$blockId))
    }
    .rtlRangeKeys(gwasFmr)
}

# Align column names between a QTL and a GWAS LBF matrix by
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
.colocEmptyResult <- function(enriched = FALSE, ldSketch = NULL) {
    base <- tibble(
        study = character(0),
        context = character(0),
        trait = character(0),
        method = character(0),
        gwasStudy = character(0),
        gwasMethod = character(0),
        blockId = character(0),
        qtlCs = integer(0),
        gwasCs = integer(0),
        idx1 = integer(0),
        idx2 = integer(0),
        nSnps = integer(0),
        PP.H0.abf = numeric(0),
        PP.H1.abf = numeric(0),
        PP.H2.abf = numeric(0),
        PP.H3.abf = numeric(0),
        PP.H4.abf = numeric(0),
        qtlRetainedMass = numeric(0),
        gwasRetainedMass = numeric(0)
    )
    if (enriched) {
        base$enrichment <- numeric(0)
        base$p12Used <- numeric(0)
    }
    ColocResult(base, list(), ldSketch = ldSketch)
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

# coloc.bf_bf names the aligned-variant count `nsnps`; this package publishes
# it as `nSnps`. Without the rename .colocStandardiseRow() invents an all-NA
# `nSnps` and the real count is left sitting in a second column beside it.
# @noRd
.colocRenameNsnps <- function(sm) {
    if (is_in("nsnps", colnames(sm)) && !is_in("nSnps", colnames(sm))) {
        sm <- rename(sm, nSnps = "nsnps")
    }
    sm
}

# Index a per-effect vector (retained mass, effect id) by coloc's effect index,
# tolerating an index coloc did not report (NA) or one past the end.
# @noRd
.colocPickAt <- function(values, idx, n, fill = NA_real_) {
    if (is.null(idx) || is.null(values)) {
        return(rep(fill, n))
    }
    idx <- as.integer(idx)
    ok <- !is.na(idx) & idx >= 1L & idx <= length(values)
    out <- rep(fill, length(idx))
    out[ok] <- values[idx[ok]]
    out
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# The variant ids of one fine-mapping entry (S4 slot; not pluckable by name).
# @noRd

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
