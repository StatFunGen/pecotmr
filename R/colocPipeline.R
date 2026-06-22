#' @title Colocalization Pipeline (coloc.bf_bf over QTL + GWAS LBF matrices)
#' @description Per-region pipeline that pairs a QTL
#'   \code{\link{QtlFineMappingResult}} with a GWAS fine-mapping result
#'   (either supplied directly as a \code{\link{GwasFineMappingResult}}
#'   or computed inline from a \code{\link{GwasSumStats}}) and runs
#'   \code{coloc::coloc.bf_bf} per (QTL tuple, GWAS tuple) pair to
#'   produce per-pair colocalization posterior probabilities
#'   PP.H0-PP.H4.
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
#' @section LD-sketch identity check:
#' If \code{getLdSketch(qtlFineMappingResult)} is non-\code{NULL}, it
#' must match the LD sketch on \code{gwasInput}. Mismatch is a hard
#' error. When the QTL FMR's \code{ldSketch} is \code{NULL}
#' (individual-level fit), the validation is skipped on the QTL side.
#'
#' @param qtlFineMappingResult A \code{\link{QtlFineMappingResult}}
#'   (required).
#' @param gwasInput Either a \code{\link{GwasSumStats}} or a
#'   \code{\link{GwasFineMappingResult}}.
#' @param filterLbfCs Logical. When \code{TRUE} (and
#'   \code{filterLbfCsSecondary} is \code{NULL}), keep only effects
#'   that produced a credible set (\code{trimmedFit$sets$cs_index}).
#'   Default \code{FALSE}.
#' @param filterLbfCsSecondary Optional secondary coverage (numeric in
#'   \eqn{(0, 1)}). When supplied, run a credible-set concentration
#'   filter at this coverage level instead of \code{filterLbfCs}: each
#'   L-effect's credible set must span fewer than \code{nVariants *
#'   filterLbfCsSecondary * filterLbfCsConcentration} variants to be
#'   kept. Effects with diffuse credible sets are dropped before the
#'   LBF matrix is passed to \code{coloc::coloc.bf_bf}. Overrides
#'   \code{filterLbfCs} when set.
#' @param filterLbfCsConcentration Numeric in \eqn{(0, 1)}; the
#'   concentration factor in the cutoff above. With the default
#'   \code{0.5} a 50\% credible set is kept only if it spans fewer
#'   than 25\% of the locus's variants. Only consulted when
#'   \code{filterLbfCsSecondary} is non-NULL. Default \code{0.5}.
#' @param priorTol Prior-variance cutoff for the default filter:
#'   effects with \code{V <= priorTol} are dropped. Ignored when
#'   either \code{filterLbfCs} or \code{filterLbfCsSecondary} is in
#'   use. Default \code{1e-9}.
#' @param p1 Prior probability of QTL signal per variant. Default
#'   \code{1e-4}.
#' @param p2 Prior probability of GWAS signal per variant. Default
#'   \code{1e-4}.
#' @param p12 Prior probability of shared signal per variant. Default
#'   \code{5e-6}.
#' @param finemappingMethods Character vector forwarded to
#'   \code{\link{fineMappingPipeline}} when \code{gwasInput} is a
#'   \code{GwasSumStats}. Default \code{"susie"}.
#' @param returnGwasFineMapping Logical. When \code{TRUE}, attach the
#'   computed \code{GwasFineMappingResult} on the returned data frame
#'   as attribute \code{"gwasFineMapping"}. Default \code{FALSE}.
#' @param enrichment Optional data.frame of per-(gwasStudy, qtlContext)
#'   enrichment factors with columns \code{gwasStudy}, \code{qtlContext},
#'   \code{enrichment}. Output of \code{\link{qtlEnrichmentPipeline}}.
#'   When non-\code{NULL}, each pair's \code{p12} prior is scaled to
#'   \code{min(p12 * (1 + enrichment), p12Max)} (the enrichment-informed
#'   colocalization variant, "enloc"). Pairs without a matching
#'   enrichment row fall back to the baseline \code{p12} with a warning.
#'   Default \code{NULL} (baseline coloc).
#' @param p12Max Numeric scalar. Maximum value for the enrichment-adjusted
#'   \code{p12} prior. Default \code{1e-3}. Ignored when
#'   \code{enrichment = NULL}.
#' @param adjustPips Logical, default \code{TRUE}. When TRUE, before any
#'   per-pair inference the QTL and GWAS fine-mapping result collections
#'   are passed through \code{\link{adjustPips}} so each entry's PIPs are
#'   renormalized to the intersection of its variants with the union of
#'   the other side's variant IDs. This matters in two scenarios: (1) the
#'   user declined to impute missing variants in the GWAS \code{SumStats}
#'   and the QTL fine-mapping input has additional variants; (2) the GWAS
#'   fine-mapping result contains variants not present in the QTL
#'   fine-mapping result. Pass \code{FALSE} to use the FMRs as supplied.
#' @param ... Additional arguments forwarded to
#'   \code{coloc::coloc.bf_bf}.
#' @return A data frame with one row per (QTL tuple, GWAS tuple,
#'   credible-set pair) combination. Identity columns: \code{study},
#'   \code{context}, \code{trait}, \code{method}, \code{gwasStudy},
#'   \code{gwasMethod}, plus the standard coloc fields
#'   (\code{idx1}, \code{idx2}, \code{nSnps},
#'   \code{PP.H0.abf} \ldots \code{PP.H4.abf}). When \code{enrichment}
#'   is supplied, two additional columns \code{enrichment} and
#'   \code{p12Used} report the per-pair factor and the prior actually
#'   passed to \code{coloc::coloc.bf_bf}.
#' @export
colocPipeline <- function(qtlFineMappingResult,
                          gwasInput,
                          filterLbfCs              = FALSE,
                          filterLbfCsSecondary     = NULL,
                          filterLbfCsConcentration = 0.5,
                          priorTol                 = 1e-9,
                          p1                       = 1e-4,
                          p2                       = 1e-4,
                          p12                      = 5e-6,
                          finemappingMethods       = "susie",
                          returnGwasFineMapping    = FALSE,
                          enrichment               = NULL,
                          p12Max                   = 1e-3,
                          adjustPips               = TRUE,
                          ...) {
  useEnrichment <- !is.null(enrichment)
  if (useEnrichment) {
    if (!is.data.frame(enrichment))
      stop("`enrichment` must be a data.frame with at least gwasStudy, ",
           "qtlContext, enrichment columns (output of ",
           "qtlEnrichmentPipeline).")
    required <- c("gwasStudy", "qtlContext", "enrichment")
    missingCols <- setdiff(required, colnames(enrichment))
    if (length(missingCols) > 0L)
      stop("`enrichment` is missing column(s): ",
           paste(missingCols, collapse = ", "))
  }
  if (!requireNamespace("coloc", quietly = TRUE)) {
    stop("Package 'coloc' is required for colocPipeline. ",
         "Install with: install.packages('coloc').")
  }
  if (!methods::is(qtlFineMappingResult, "QtlFineMappingResult")) {
    stop("`qtlFineMappingResult` must be a QtlFineMappingResult ",
         "(got class '", class(qtlFineMappingResult)[[1L]], "').")
  }
  if (!methods::is(gwasInput, "GwasSumStats") &&
      !methods::is(gwasInput, "GwasFineMappingResult")) {
    stop("`gwasInput` must be a GwasSumStats or a ",
         "GwasFineMappingResult (got class '",
         class(gwasInput)[[1L]], "').")
  }

  # --- Resolve gwas fine-mapping --------------------------------------
  gwasFmr <- if (methods::is(gwasInput, "GwasFineMappingResult")) {
    gwasInput
  } else {
    if (length(getQcInfo(gwasInput)) == 0L) {
      stop("colocPipeline: gwasInput (GwasSumStats) has no QC record. ",
           "Call summaryStatsQc() first.")
    }
    fineMappingPipeline(gwasInput, methods = finemappingMethods)
  }

  # --- LD-sketch identity check ---------------------------------------
  qtlLd <- getLdSketch(qtlFineMappingResult)
  gwasLd <- getLdSketch(gwasFmr)
  .colocRequireMatchingLdSketches(qtlLd, gwasLd)

  # --- Optional PIP renormalization to the cross-FMR variant union -----
  # Renormalize each side's entry PIPs to the intersection of its own
  # variants with the union of the other side's variants. Handles two
  # cases: (a) GWAS sumstats missing variants present in the QTL FMR
  # (no imputation requested); (b) GWAS FMR carrying variants not in
  # the QTL FMR.
  if (isTRUE(adjustPips)) {
    qtlVids  <- unique(unlist(lapply(qtlFineMappingResult$entry,
                                     function(e) e@variantIds)))
    gwasVids <- unique(unlist(lapply(gwasFmr$entry,
                                     function(e) e@variantIds)))
    if (length(qtlVids) > 0L && length(gwasVids) > 0L) {
      qtlFineMappingResult <- adjustPips(qtlFineMappingResult, gwasVids)
      gwasFmr              <- adjustPips(gwasFmr, qtlVids)
    }
  }

  # --- Pre-extract per-GWAS-tuple LBF matrices ------------------------
  # The legacy colocWrapper combined multiple GWAS files' LBF matrices
  # row-wise before a single coloc.bf_bf call per xQTL. Reproduce that
  # by grouping the GWAS FMR by study, stacking each study's LBF rows,
  # and storing per-(study, method) batched matrices.
  gwasLbfByPair <- .colocPreextractGwasLbf(
    gwasFmr, filterLbfCs, filterLbfCsSecondary,
    filterLbfCsConcentration, priorTol)
  if (length(gwasLbfByPair) == 0L) {
    out <- .colocEmptyResult(enriched = useEnrichment)
    if (returnGwasFineMapping && methods::is(gwasInput, "GwasSumStats")) {
      attr(out, "gwasFineMapping") <- gwasFmr
    }
    return(out)
  }

  # --- Per (QTL tuple, GWAS tuple) pair iteration ---------------------
  results <- list()
  for (qi in seq_len(nrow(qtlFineMappingResult))) {
    qStudy   <- as.character(qtlFineMappingResult$study)[[qi]]
    qContext <- as.character(qtlFineMappingResult$context)[[qi]]
    qTrait   <- as.character(qtlFineMappingResult$trait)[[qi]]
    qMethod  <- as.character(qtlFineMappingResult$method)[[qi]]
    qEntry   <- qtlFineMappingResult$entry[[qi]]
    qLbfInfo <- .colocExtractLbfFromEntry(
      qEntry, filterLbfCs, filterLbfCsSecondary,
      filterLbfCsConcentration, priorTol,
      label = sprintf("QTL (study='%s', context='%s', trait='%s', method='%s')",
                      qStudy, qContext, qTrait, qMethod))
    if (is.null(qLbfInfo)) next
    qLbf <- qLbfInfo$lbf

    for (gKey in names(gwasLbfByPair)) {
      gInfo  <- gwasLbfByPair[[gKey]]
      gLbf   <- gInfo$lbf

      # Align variants between QTL and GWAS LBF matrices via the legacy
      # alignVariantNames + intersect-and-drop pattern.
      aligned <- .colocAlignLbf(qLbf, gLbf)
      if (is.null(aligned)) next
      qAligned <- aligned$qtl
      gAligned <- aligned$gwas

      # Enrichment-informed p12: per-pair scaling capped at p12Max.
      # Baseline p12 used when no enrichment table or no matching row.
      if (useEnrichment) {
        enRow <- .colocLookupEnrichment(enrichment, gInfo$study, qContext)
        if (is.na(enRow)) {
          warning(sprintf(
            "colocPipeline: no enrichment entry for (gwasStudy='%s', qtlContext='%s'); using baseline p12.",
            gInfo$study, qContext))
          enRow <- 0
        }
        p12Used <- min(p12 * (1 + enRow), p12Max)
      } else {
        enRow   <- NA_real_
        p12Used <- p12
      }

      pairRes <- tryCatch(
        coloc::coloc.bf_bf(qAligned, gAligned,
                           p1 = p1, p2 = p2, p12 = p12Used, ...),
        error = function(e) {
          warning(sprintf(
            "colocPipeline: coloc.bf_bf failed for QTL (study='%s', context='%s', trait='%s', method='%s') x GWAS (study='%s', method='%s'): %s",
            qStudy, qContext, qTrait, qMethod,
            gInfo$study, gInfo$method, conditionMessage(e)))
          NULL
        })
      if (is.null(pairRes) || is.null(pairRes$summary)) next

      sm <- as.data.frame(pairRes$summary, stringsAsFactors = FALSE)
      sm$study      <- qStudy
      sm$context    <- qContext
      sm$trait      <- qTrait
      sm$method     <- qMethod
      sm$gwasStudy  <- gInfo$study
      sm$gwasMethod <- gInfo$method
      if (useEnrichment) {
        sm$enrichment <- enRow
        sm$p12Used    <- p12Used
      }
      results[[length(results) + 1L]] <- sm
    }
  }

  if (length(results) == 0L) {
    out <- .colocEmptyResult(enriched = useEnrichment)
  } else {
    out <- do.call(rbind, lapply(results, .colocStandardiseRow))
    rownames(out) <- NULL
    idCols <- c("study", "context", "trait", "method",
                "gwasStudy", "gwasMethod")
    if (useEnrichment) idCols <- c(idCols, "enrichment", "p12Used")
    other  <- setdiff(colnames(out), idCols)
    out    <- out[, c(idCols, other), drop = FALSE]
  }

  if (returnGwasFineMapping && methods::is(gwasInput, "GwasSumStats")) {
    attr(out, "gwasFineMapping") <- gwasFmr
  }
  out
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
.colocFilterCsByConcentration <- function(fit, coverage = 0.5,
                                          concentration = 0.5) {
  fit$V <- NULL  # disable V-based filtering inside susie_get_cs
  csList <- susie_get_cs(fit, coverage = coverage, dedup = FALSE)
  totalVariants <- ncol(fit$alpha)
  maxSize <- totalVariants * coverage * concentration
  keep <- map_lgl(csList$cs, ~ length(.x) < maxSize)
  as.numeric(gsub("L", "", names(which(keep))))
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
.colocExtractLbfFromEntry <- function(entry, filterLbfCs,
                                      filterLbfCsSecondary,
                                      filterLbfCsConcentration,
                                      priorTol,
                                      label = "entry") {
  fit <- getSusieFit(entry)
  if (is.null(fit)) {
    warning(sprintf("colocPipeline: %s has no trimmedFit; skipping.",
                    label))
    return(NULL)
  }

  lbfMatrix <- if (!is.null(fit$lbf_variable)) {
    as.matrix(fit$lbf_variable)
  } else if (!is.null(fit$fsusie_result) &&
             is.list(fit$fsusie_result$lBF)) {
    # fSuSiE path: stack per-trait lBF lists into a single matrix.
    do.call(rbind, fit$fsusie_result$lBF)
  } else if (is.list(fit) && length(fit) >= 1L &&
             !is.null(fit[[1L]]$fsusie_result$lBF)) {
    do.call(rbind, fit[[1L]]$fsusie_result$lBF)
  } else {
    warning(sprintf(
      "colocPipeline: %s trimmedFit has no lbf_variable / fsusie lBF; skipping.",
      label))
    return(NULL)
  }
  if (is.null(lbfMatrix) || nrow(lbfMatrix) == 0L) {
    warning(sprintf("colocPipeline: %s LBF matrix is empty.", label))
    return(NULL)
  }

  # Row (effect) filtering, in priority order of the original code.
  if (isTRUE(filterLbfCs) && is.null(filterLbfCsSecondary)) {
    csIdx <- fit$sets$cs_index
    if (!is.null(csIdx) && length(csIdx) > 0L)
      lbfMatrix <- lbfMatrix[csIdx, , drop = FALSE]
  } else if (!is.null(filterLbfCsSecondary)) {
    secIdx <- tryCatch(
      .colocFilterCsByConcentration(fit,
                                    coverage = filterLbfCsSecondary,
                                    concentration = filterLbfCsConcentration),
      error = function(e) NULL)
    if (!is.null(secIdx) && length(secIdx) > 0L)
      lbfMatrix <- lbfMatrix[secIdx, , drop = FALSE]
  } else {
    if (!is.null(fit$V))
      lbfMatrix <- lbfMatrix[fit$V > priorTol, , drop = FALSE]
  }
  if (nrow(lbfMatrix) == 0L) return(NULL)

  # Column-name (variant id) assignment. Prefer fit-provided names;
  # fall back to the entry's variantIds slot.
  if (is.null(colnames(lbfMatrix)) || any(is.na(colnames(lbfMatrix)))) {
    vids <- getVariantIds(entry)
    if (length(vids) == ncol(lbfMatrix)) {
      colnames(lbfMatrix) <- vids
    }
  }
  lbfMatrix <- lbfMatrix[, !is.na(colnames(lbfMatrix)), drop = FALSE]
  if (ncol(lbfMatrix) == 0L) return(NULL)

  list(lbf = lbfMatrix, variantIds = colnames(lbfMatrix))
}

# Build a per-GWAS-tuple LBF matrix list, keyed by "study|method".
# Within each key we stack multiple FMR rows row-wise (the legacy
# "combined GWAS LBF" pattern), drop NA columns, and replace NAs with
# 0 so a fresh QTL pairing always lands on the same coordinate frame.
# @noRd
.colocPreextractGwasLbf <- function(gwasFmr, filterLbfCs,
                                    filterLbfCsSecondary,
                                    filterLbfCsConcentration,
                                    priorTol) {
  groupKey <- paste(as.character(gwasFmr$study),
                    as.character(gwasFmr$method),
                    sep = "||")
  groups <- split(seq_len(nrow(gwasFmr)), groupKey)
  out <- list()
  for (gkey in names(groups)) {
    rows <- groups[[gkey]]
    pieces <- list()
    for (ri in rows) {
      info <- .colocExtractLbfFromEntry(
        gwasFmr$entry[[ri]],
        filterLbfCs, filterLbfCsSecondary,
        filterLbfCsConcentration, priorTol,
        label = sprintf("GWAS (study='%s', method='%s', row=%d)",
                        as.character(gwasFmr$study)[[ri]],
                        as.character(gwasFmr$method)[[ri]], ri))
      if (!is.null(info)) pieces[[length(pieces) + 1L]] <- info$lbf
    }
    if (length(pieces) == 0L) next
    combined <- .colocRbindLbf(pieces)
    parts <- strsplit(gkey, "\\|\\|")[[1L]]
    out[[gkey]] <- list(lbf = combined,
                        study  = parts[[1L]],
                        method = if (length(parts) >= 2L) parts[[2L]]
                                 else NA_character_)
  }
  out
}

# Row-bind a list of LBF matrices with potentially-different column
# sets. Uses the union of columns; missing cells fill with 0 (matching
# the legacy `replace_na(., 0)` step inside colocWrapper).
# @noRd
.colocRbindLbf <- function(mats) {
  allCols <- unique(unlist(lapply(mats, colnames)))
  padded <- lapply(mats, function(m) {
    miss <- setdiff(allCols, colnames(m))
    if (length(miss) > 0L) {
      pad <- matrix(0, nrow = nrow(m), ncol = length(miss),
                    dimnames = list(NULL, miss))
      m <- cbind(m, pad)
    }
    m[, allCols, drop = FALSE]
  })
  do.call(rbind, padded)
}

# Align column names between a QTL and a (combined) GWAS LBF matrix,
# intersect to the common variants, and return both restricted to that
# common set in the same order. Mirrors alignVariantNames +
# intersect-and-drop in the legacy code.
# @noRd
.colocAlignLbf <- function(qtlLbf, gwasLbf) {
  qids <- colnames(qtlLbf)
  gids <- colnames(gwasLbf)
  aligned <- tryCatch(alignVariantNames(qids, gids),
                      error = function(e) NULL)
  if (!is.null(aligned)) {
    qids <- aligned$alignedVariants
    colnames(qtlLbf) <- qids
  }
  common <- intersect(qids, gids)
  if (length(common) == 0L) return(NULL)
  list(qtl  = qtlLbf[, common, drop = FALSE],
       gwas = gwasLbf[, common, drop = FALSE])
}

# A blank result data frame for the no-pair case so callers downstream
# do not have to special-case a NULL return. When `enriched = TRUE` the
# enrichment + p12Used columns are appended (the enloc-mode schema).
# @noRd
.colocEmptyResult <- function(enriched = FALSE) {
  base <- data.frame(
    study      = character(0),
    context    = character(0),
    trait      = character(0),
    method     = character(0),
    gwasStudy  = character(0),
    gwasMethod = character(0),
    idx1       = integer(0),
    idx2       = integer(0),
    nSnps      = integer(0),
    PP.H0.abf  = numeric(0),
    PP.H1.abf  = numeric(0),
    PP.H2.abf  = numeric(0),
    PP.H3.abf  = numeric(0),
    PP.H4.abf  = numeric(0),
    stringsAsFactors = FALSE)
  if (enriched) {
    base$enrichment <- numeric(0)
    base$p12Used    <- numeric(0)
  }
  base
}

# Look up the enrichment factor for a (gwasStudy, qtlContext) pair in
# the user-supplied enrichment table. Returns NA when the pair is not
# present; the caller falls back to the baseline p12 and emits a
# warning.
# @noRd
.colocLookupEnrichment <- function(enrichment, gwasStudy, qtlContext) {
  idx <- which(as.character(enrichment$gwasStudy)  == gwasStudy &
               as.character(enrichment$qtlContext) == qtlContext)
  if (length(idx) == 0L) return(NA_real_)
  as.numeric(enrichment$enrichment[[idx[[1L]]]])
}

# Ensure each row data.frame from coloc.bf_bf carries the standard PP
# columns even when the underlying call produced a slightly different
# shape.
# @noRd
.colocStandardiseRow <- function(sm) {
  for (col in c("idx1", "idx2", "nSnps",
                "PP.H0.abf", "PP.H1.abf", "PP.H2.abf",
                "PP.H3.abf", "PP.H4.abf")) {
    if (!col %in% colnames(sm)) sm[[col]] <- NA
  }
  sm
}
