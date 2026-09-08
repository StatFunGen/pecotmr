#' @title QTL Enrichment Pipeline (Genome-Wide)
#' @description Genome-wide pipeline that computes per-pair enrichment
#'   estimates by passing an outcome PIP vector and a set of annotation
#'   credible-set posteriors to \code{\link{qtlEnrichment}}. The returned table
#'   feeds \code{\link{colocPipeline}} via its \code{enrichment} argument.
#'
#'   Either side may be a \code{\link{QtlFineMappingResult}} or a
#'   \code{\link{GwasFineMappingResult}}, so QTL-in-GWAS, QTL-in-QTL and
#'   GWAS-in-GWAS enrichment all run through one code path. The argument names
#'   keep the QTL / GWAS wording of the common case; what they mean generally
#'   is that \code{gwasFineMappingResult} is the \strong{outcome} whose PIPs
#'   are scanned and \code{qtlFineMappingResult} is the \strong{annotation}
#'   whose region fits are tested for enrichment within them.
#'
#' \strong{Not gene-parallelisable}: the enrichment estimator runs over the full
#' genome of outcome PIPs and the full collection of annotation fits at once.
#'
#' @section Inputs:
#' \itemize{
#'   \item \code{gwasFineMappingResult}: the genome-wide outcome collection.
#'     Each entry's \code{FineMappingRow$trimmedFit} must carry a \code{pip}
#'     vector. One PIP vector is built per outcome trait -- keyed by
#'     \code{study} for a GWAS collection (one trait per study) and by
#'     (\code{study}, \code{context}, \code{trait}) for a QTL one, since a
#'     variant's PIP differs between molecular traits and pooling them would
#'     collide.
#'   \item \code{qtlFineMappingResult}: the annotation collection. Each entry's
#'     \code{trimmedFit} must carry \code{alpha}, \code{pip}, and
#'     prior-variance fields (\code{V}). Its region fits are pooled per
#'     (\code{study}, \code{context}) -- per \code{study} alone for a GWAS
#'     collection, which has no context axis.
#' }
#'
#' @section LD-sketch identity check: A GWAS outcome collection must have a
#'   non-NULL \code{ldSketch} (it should be RSS-derived). Where both sides
#'   carry one, the two must match exactly; a \code{NULL} on either side (an
#'   individual-level fit) skips the check.
#'
#' @param gwasFineMappingResult The outcome side; see above.
#' @param qtlFineMappingResult The annotation side; see above.
#' @param numGwas Number of GWAS variants used to estimate \code{piGwas}. When
#'   \code{NULL} (default) it is estimated from the data -- bias warning applies
#'   if the input PIP vector is not genome-wide.
#' @param piQtl Per-variant prior of being a QTL causal variant. \code{NULL}
#'   (default) estimates from the data.
#' @param lambda Shrinkage parameter for the enrichment estimator. Default
#'   \code{1.0}.
#' @param impN Number of imputed samples used by the estimator. Default
#'   \code{25}.
#' @param numThreads Number of threads used by \code{qtlEnrichment}. Default
#'   \code{1}.
#' @param seed Integer or \code{NULL}. Base random seed forwarded to
#'   \code{\link{qtlEnrichment}} for reproducible multiple imputation.
#'   \code{NULL} (default) draws a nondeterministic seed.
#' @param ... Additional arguments forwarded to \code{\link{qtlEnrichment}}.
#' @return A tibble with one row per (outcome trait, annotation unit) pair.
#'   The identity columns are \code{gwasStudy}, \code{gwasContext},
#'   \code{gwasTrait}, \code{qtlStudy}, \code{qtlContext}; the axes a side does
#'   not have are \code{NA} (\code{gwasContext} / \code{gwasTrait} for a GWAS
#'   outcome, \code{qtlContext} for a GWAS annotation). Suitable as the
#'   \code{enrichment} argument to \code{\link{colocPipeline}}, which joins on
#'   those columns.
#'
#'   The estimates are \code{enrichmentLogOdds}, the enrichment parameter
#'   \eqn{a_1} on the log-odds scale, with its standard error
#'   \code{enrichmentSe}; \code{enrichment} is the same quantity as a
#'   multiplicative factor, \eqn{e^{a_1} - 1}, which is what
#'   \code{colocPipeline} scales \code{p12} by (so \eqn{a_1 = 0} leaves the
#'   prior untouched). \code{enrichmentLogOddsNoShrinkage} and
#'   \code{enrichmentSeNoShrinkage} are the same estimate before shrinkage,
#'   \code{intercept} / \code{interceptSe} are \eqn{a_0}, and \code{colocP1},
#'   \code{colocP2}, \code{colocP12} are the enrichment-informed coloc priors
#'   the estimator derives from \eqn{(a_0, a_1)} -- an alternative to scaling
#'   a baseline \code{p12}. \code{effectiveMiRounds} is how many
#'   multiple-imputation rounds survived outlier filtering.
#' @examples
#' data(gwasFineMappingExample)
#' data(qtlFineMappingExample)
#' qtlEnrichmentPipeline(
#'   gwasFineMappingResult = gwasFineMappingExample,
#'   qtlFineMappingResult = qtlFineMappingExample
#' )
#' @export
qtlEnrichmentPipeline <- function(
    gwasFineMappingResult,
    qtlFineMappingResult,
    numGwas = NULL,
    piQtl = NULL,
    lambda = 1.0,
    impN = 25,
    numThreads = 1L,
    seed = NULL,
    ...
) {
    .enrValidateInputs(gwasFineMappingResult, qtlFineMappingResult)
    p <- as.list(environment())
    p$dots <- list(...)
    p <- .enrPrepare(p)
    results <- list_flatten(map(
        seq_len(nrow(p$gwasTuples)),
        .enrScoreOutcomeTuple,
        p = p
    ))
    .enrAssemble(results)
}

# Validate the input classes + LD-sketch presence / identity.
# @noRd
.enrValidateInputs <- function(gwasFineMappingResult, qtlFineMappingResult) {
    if (!methods::is(gwasFineMappingResult, "FineMappingResultBase")) {
        msg <- glue(
            "`gwasFineMappingResult` must be a GwasFineMappingResult or a ",
            "QtlFineMappingResult ",
            "(got class '{class(gwasFineMappingResult)[[1L]]}')."
        )
        abort(msg)
    }
    if (!methods::is(qtlFineMappingResult, "FineMappingResultBase")) {
        msg <- glue(
            "`qtlFineMappingResult` must be a QtlFineMappingResult or a ",
            "GwasFineMappingResult ",
            "(got class '{class(qtlFineMappingResult)[[1L]]}')."
        )
        abort(msg)
    }
    outcomeLd <- getLdSketch(gwasFineMappingResult)
    if (
        is.null(outcomeLd) &&
            methods::is(gwasFineMappingResult, "GwasFineMappingResult")
    ) {
        msg <- glue(
            "qtlEnrichmentPipeline: the GWAS FineMappingResult must have a ",
            "non-NULL ldSketch (it should be RSS-derived)."
        )
        abort(msg)
    }
    # Lenient rather than qtl-required: a QTL outcome side may be an
    # individual-level fit carrying no panel, and the requirement above already
    # covers the RSS-derived GWAS case.
    .requireMatchingLdSketches(
        getLdSketch(qtlFineMappingResult),
        outcomeLd,
        pipelineName = "qtlEnrichmentPipeline",
        nullPolicy = "lenient"
    )
    invisible(NULL)
}

# Hoist the outcome-independent work out of the double loop: per-trait outcome
# PIP vectors, the union variant-name panel, per-tuple annotation regions, and
# each tuple's one-time alignment to the union panel (errors captured as
# values).
# @noRd
.enrPrepare <- function(p) {
    p$gwasTuples <- .enrOutcomeTuples(p$gwasFineMappingResult)
    p$qtlTuples <- .enrAnnotationTuples(p$qtlFineMappingResult)
    if (nrow(p$gwasTuples) == 0L || nrow(p$qtlTuples) == 0L) {
        msg <- glue(
            "qtlEnrichmentPipeline: no (outcome, annotation) pairs to ",
            "compute (one of the inputs has zero rows)."
        )
        abort(msg)
    }
    p$gwasPipByTuple <- map(
        seq_len(nrow(p$gwasTuples)),
        .enrGwasPipForRow,
        gwasTuples = p$gwasTuples,
        fmr = p$gwasFineMappingResult
    )
    unionGwasNames <- unique(unlist(
        map(p$gwasPipByTuple, names),
        use.names = FALSE
    ))
    p$qtlRegionsByTuple <- map(
        seq_len(nrow(p$qtlTuples)),
        .enrQtlRegionsForRow,
        qtlTuples = p$qtlTuples,
        fmr = p$qtlFineMappingResult
    )
    p$alignedByTuple <- map(
        p$qtlRegionsByTuple,
        .enrAlignRegionsSafe,
        unionGwasNames = unionGwasNames
    )
    p
}

# The outcome side's per-trait keys: one PIP vector is built per key. A GWAS
# collection carries one trait per study, so its QTL-only axes read NA and the
# key collapses to (study).
# @noRd
.enrOutcomeTuples <- function(fmr) {
    distinct(tibble(
        gwasStudy = .fmrIdentityColumn(fmr, "study"),
        gwasContext = .fmrIdentityColumn(fmr, "context"),
        gwasTrait = .fmrIdentityColumn(fmr, "trait")
    ))
}

# The annotation side's keys: (study, context), the unit whose region fits are
# pooled into one enrichment estimate. The joint key matters -- context alone
# would silently merge two studies sharing a context label -- and a GWAS
# annotation pools its blocks under one NA-context key.
# @noRd
.enrAnnotationTuples <- function(fmr) {
    distinct(tibble(
        qtlStudy = .fmrIdentityColumn(fmr, "study"),
        qtlContext = .fmrIdentityColumn(fmr, "context")
    ))
}

# One tuple-table row as an identity list addressed by the collection's OWN
# column names; the table prefixes them so the result can name both sides.
# @noRd
.enrOutcomeIdent <- function(gwasTuples, k) {
    list(
        study = gwasTuples$gwasStudy[[k]],
        context = gwasTuples$gwasContext[[k]],
        trait = gwasTuples$gwasTrait[[k]]
    )
}

# @noRd
.enrAnnotationIdent <- function(qtlTuples, k) {
    list(
        study = qtlTuples$qtlStudy[[k]],
        context = qtlTuples$qtlContext[[k]]
    )
}

# Row indices of `fmr` matching an identity tuple. An axis the collection does
# not have is NA on both sides and matches, rather than excluding every row.
# @noRd
.enrMatchRows <- function(fmr, ident) {
    hits <- map(names(ident), .enrColumnMatches, fmr = fmr, ident = ident)
    which(reduce(hits, `&`, .init = rep(TRUE, nrow(fmr))))
}

# @noRd
.enrColumnMatches <- function(column, fmr, ident) {
    values <- .fmrIdentityColumn(fmr, column)
    wanted <- ident[[column]]
    if (is.na(wanted)) {
        return(is.na(values))
    }
    !is.na(values) & values == wanted
}

# The outcome PIP vector for the k-th outcome tuple.
# @noRd
.enrGwasPipForRow <- function(k, gwasTuples, fmr) {
    .enrBuildGwasPipVector(fmr, .enrOutcomeIdent(gwasTuples, k))
}

# Annotation SuSiE region list for the k-th (study, context) tuple.
# @noRd
.enrQtlRegionsForRow <- function(k, qtlTuples, fmr) {
    .enrBuildQtlRegionsList(fmr, .enrAnnotationIdent(qtlTuples, k))
}

# Align one tuple's regions to the union GWAS panel, capturing any error as a
# value (re-raised + skipped per (gwas, tuple) below, never aborting).
# @noRd
.enrAlignRegionsSafe <- function(regions, unionGwasNames) {
    tryCatch(
        .enrAlignRegions(regions, unionGwasNames),
        error = function(e) e
    )
}

# Score one outcome trait against every annotation tuple -> enrichment records
# (empty when the outcome has no usable PIPs).
# @noRd
.enrScoreOutcomeTuple <- function(gi, p) {
    gwasPip <- p$gwasPipByTuple[[gi]]
    if (length(gwasPip) == 0L) {
        msg <- glue(
            "qtlEnrichmentPipeline: no usable PIPs for ",
            "{.enrOutcomeLabel(p, gi)}; skipping."
        )
        warn(msg)
        return(list())
    }
    compact(map(
        seq_len(nrow(p$qtlTuples)),
        .enrScoreTuple,
        gi = gi,
        gwasPip = gwasPip,
        p = p
    ))
}

# Score one (outcome trait, annotation tuple) pair -> an enrichment record, or
# NULL when the tuple has no regions or qtlEnrichment fails.
# @noRd
.enrScoreTuple <- function(k, gi, gwasPip, p) {
    if (length(p$qtlRegionsByTuple[[k]]) == 0L) {
        msg <- glue(
            "qtlEnrichmentPipeline: no usable regions for ",
            "{.enrAnnotationLabel(p, k)}; skipping."
        )
        warn(msg)
        return(NULL)
    }
    enr <- .enrRunEnrichment(gi, gwasPip, k, p)
    if (is.null(enr)) {
        return(NULL)
    }
    c(
        .enrFlattenEnrichment(enr),
        as.list(p$gwasTuples[gi, , drop = FALSE]),
        as.list(p$qtlTuples[k, , drop = FALSE])
    )
}

# Human-readable identities for the warnings above, naming each side by its own
# flavour and only the axes it has.
# @noRd
.enrOutcomeLabel <- function(p, gi) {
    .fmrTupleLabel(
        .fmrSideName(p$gwasFineMappingResult),
        .enrOutcomeIdent(p$gwasTuples, gi)
    )
}

# @noRd
.enrAnnotationLabel <- function(p, k) {
    .fmrTupleLabel(
        .fmrSideName(p$qtlFineMappingResult),
        .enrAnnotationIdent(p$qtlTuples, k)
    )
}

# Run qtlEnrichment for a pair (with the pre-aligned regions), warning + NULL on
# failure. alignNames = FALSE reuses the shared per-tuple alignment.
# @noRd
.enrRunEnrichment <- function(gi, gwasPip, k, p) {
    aligned <- p$alignedByTuple[[k]]
    tryCatch(
        {
            if (inherits(aligned, "condition")) {
                cnd_signal(aligned)
            }
            enrichArgs <- c(
                list(
                    gwasPip = gwasPip,
                    susieQtlRegions = aligned,
                    numGwas = p$numGwas,
                    piQtl = p$piQtl,
                    lambda = p$lambda,
                    impN = p$impN,
                    numThreads = p$numThreads,
                    seed = p$seed,
                    alignNames = FALSE
                ),
                p$dots
            )
            exec(qtlEnrichment, !!!enrichArgs)
        },
        error = function(e) {
            eMsg <- conditionMessage(e)
            msg <- glue(
                "qtlEnrichmentPipeline: qtlEnrichment failed for ",
                "{.enrOutcomeLabel(p, gi)} x ",
                "{.enrAnnotationLabel(p, k)}: {eMsg}"
            )
            warn(msg)
            NULL
        }
    )
}

# Row-bind the enrichment records into the id-first result table.
# @noRd
.enrAssemble <- function(results) {
    if (length(results) == 0L) {
        return(.enrEmptyResult())
    }
    out <- bind_rows(results)
    select(out, all_of(.enrIdCols()), everything())
}

# The identity columns of a result row: the outcome trait's tuple, then the
# annotation unit's.
# @noRd
.enrIdCols <- function() {
    c("gwasStudy", "gwasContext", "gwasTrait", "qtlStudy", "qtlContext")
}

# The empty enrichment result table, with the same columns a populated one has.
# @noRd
.enrEmptyResult <- function() {
    idCols <- set_names(
        rep(list(character(0)), length(.enrIdCols())),
        .enrIdCols()
    )
    valueNames <- names(.enrNaEnrichment())
    valueCols <- set_names(
        rep(list(numeric(0)), length(valueNames)),
        valueNames
    )
    as_tibble(c(idCols, valueCols))
}

# =============================================================================
# Internal helpers
# =============================================================================

# Align one tuple's QTL regions to the GWAS naming convention: relabel matched
# pip names to the union GWAS panel via the shared matcher (unmatched names
# kept as-is). Pure -- the caller precomputes one result per tuple and shares
# it across the outer GWAS loop.
# @noRd
.enrAlignRegions <- function(regions, unionGwasNames) {
    map(regions, .enrAlignRegion, unionGwasNames = unionGwasNames)
}

# Build a named outcome PIP vector for one trait. Walks every row of the
# collection carrying that identity, extracts the per-row pip from each
# FineMappingRow, and concatenates with variant-id names. Errors if any single
# variant appears with conflicting PIP values across rows -- which is why the
# identity is the full trait tuple: pooling two molecular traits of one study
# would collide on every variant they share.
#' @importFrom dplyr add_count
#' @noRd
.enrBuildGwasPipVector <- function(gwasFmr, ident) {
    idx <- .enrMatchRows(gwasFmr, ident)
    if (length(idx) == 0L) {
        return(numeric(0))
    }
    pieces <- list()
    for (i in idx) {
        parts <- .fmrRowParts(gwasFmr, i)
        fit <- getSusieFit(parts)
        if (is.null(fit) || is.null(fit$pip)) {
            next
        }
        pip <- as.numeric(fit$pip)
        ids <- if (!is.null(names(fit$pip))) {
            names(fit$pip)
        } else {
            .fmrPartsVariantIds(parts)
        }
        if (length(ids) != length(pip)) {
            next
        }
        pieces[[length(pieces) + 1L]] <-
            set_names(pip, as.character(ids))
    }
    if (length(pieces) == 0L) {
        return(numeric(0))
    }
    all <- unlist(pieces)
    if (n_distinct(names(all)) < length(all)) {
        all <- .enrCollapseDuplicatePips(all)
    }
    all
}

# Collapse duplicate variant ids across GWAS blocks: agreeing PIPs (rounded to
# 12 digits) merge; a variant with conflicting PIPs aborts. distinct() keeps one
# row per (id, rounded-pip), so a conflicting variant survives with >1 row for
# add_count() to flag.
# @noRd
.enrCollapseDuplicatePips <- function(all) {
    byId <- tibble(id = names(all), pip = as.numeric(all)) |>
        mutate(pipR = round(.data$pip, 12)) |>
        distinct(.data$id, .data$pipR, .keep_all = TRUE) |>
        add_count(.data$id)
    conflict <- filter(byId, .data$n > 1L)
    if (nrow(conflict) > 0L) {
        vName <- conflict$id[[1L]]
        msg <- glue(
            "qtlEnrichmentPipeline: variant '{vName}' appears with ",
            "conflicting PIPs across GWAS blocks; the GWAS fine-mapping ",
            "must produce a consistent PIP per variant."
        )
        abort(msg)
    }
    set_names(byId$pip, byId$id)
}

# Build the per-(study, context) list of region fits in the shape that
# qtlEnrichment expects: list(d) where each d carries alpha, pip,
# prior_variance (V). Filters on BOTH study and context so entries from
# different studies that happen to share a context label are not pooled into
# one enrichment estimate.
# @noRd
.enrBuildQtlRegionsList <- function(qtlFmr, ident) {
    idx <- .enrMatchRows(qtlFmr, ident)
    if (length(idx) == 0L) {
        return(list())
    }
    out <- list()
    for (i in idx) {
        parts <- .fmrRowParts(qtlFmr, i)
        fit <- getSusieFit(parts)
        if (is.null(fit) || is.null(fit$alpha) || is.null(fit$pip)) {
            next
        }
        pV <- if (!is.null(fit$V)) {
            fit$V
        } else if (!is.null(fit$prior_variance)) {
            fit$prior_variance
        } else {
            NULL
        }
        if (is.null(pV)) {
            next
        }
        if (is.null(names(fit$pip))) {
            names(fit$pip) <- .fmrPartsVariantIds(parts)
        }
        out[[length(out) + 1L]] <- list(
            alpha = fit$alpha,
            pip = fit$pip,
            prior_variance = pV
        )
    }
    out
}

# Pull one enrichment field from qtlEnrichment's list output as a scalar numeric
# (NA when absent).
# @noRd
.enrPickScalar <- function(field, enr) {
    v <- enr[[field]]
    if (is.null(v)) {
        NA_real_
    } else {
        as.numeric(v[[1L]])
    }
}

# Project qtlEnrichment()'s output onto the columns this pipeline publishes.
#
# The field names are the estimator's own, written in src/qtl_enrichment.h and
# shared verbatim with upstream fastenloc's enloc.enrich.out, so they are
# matched literally rather than guessed at. The shrinkage estimates are the
# ones reported, matching upstream, which likewise switches to the shrunk a1
# before deriving the coloc priors.
#
# `enrichment` is expm1 of the log-odds a1 rather than a1 itself, because
# colocPipeline consumes it as `p12 * (1 + enrichment)`, which is then exactly
# the enloc-adjusted prior `p12 * exp(a1)` -- and leaves p12 untouched for an
# unenriched annotation (a1 = 0). expm1 is bounded below by -1, so a depleted
# annotation shrinks p12 towards 0 rather than past it.
# @noRd
.enrFlattenEnrichment <- function(enr) {
    if (!is.list(enr) || !is_in("Enrichment (w/ shrinkage)", names(enr))) {
        msg <- glue(
            "qtlEnrichmentPipeline: the enrichment estimator returned no ",
            "'Enrichment (w/ shrinkage)' field, so every estimate is ",
            "reported as NA. Expected the field names written by ",
            "qtl_enrichment.h."
        )
        warn(msg)
        return(.enrNaEnrichment())
    }
    logOdds <- .enrPickScalar("Enrichment (w/ shrinkage)", enr)
    list(
        enrichment = expm1(logOdds),
        enrichmentSe = .enrPickScalar("sd (w/ shrinkage)", enr),
        enrichmentLogOdds = logOdds,
        enrichmentLogOddsNoShrinkage = .enrPickScalar(
            "Enrichment (no shrinkage)",
            enr
        ),
        enrichmentSeNoShrinkage = .enrPickScalar("sd (no shrinkage)", enr),
        intercept = .enrPickScalar("Intercept", enr),
        interceptSe = .enrPickScalar("sd (intercept)", enr),
        colocP1 = .enrPickScalar("Alternative (coloc) p1", enr),
        colocP2 = .enrPickScalar("Alternative (coloc) p2", enr),
        colocP12 = .enrPickScalar("Alternative (coloc) p12", enr),
        effectiveMiRounds = .enrPickScalar("Effective MI rounds", enr)
    )
}

# The value columns, all unmeasured. Also the single source of the value-column
# schema, so the empty result cannot drift from the populated one.
# @noRd
.enrNaEnrichment <- function() {
    list(
        enrichment = NA_real_,
        enrichmentSe = NA_real_,
        enrichmentLogOdds = NA_real_,
        enrichmentLogOddsNoShrinkage = NA_real_,
        enrichmentSeNoShrinkage = NA_real_,
        intercept = NA_real_,
        interceptSe = NA_real_,
        colocP1 = NA_real_,
        colocP2 = NA_real_,
        colocP12 = NA_real_,
        effectiveMiRounds = NA_real_
    )
}


# =============================================================================
# qtlEnrichment: low-level enrichment estimation
# -----------------------------------------------------------------------------
# Per-(GWAS, QTL-region-list) enrichment estimator. Called per-(gwasStudy,
# qtlContext) pair by qtlEnrichmentPipeline above. Uses the fastenloc-style
# C++ kernel (qtlEnrichmentRcpp) under the hood.
# =============================================================================
#' @title Implementation of enrichment analysis described in
#'   https://doi.org/10.1371/journal.pgen.1006646
#'
#' @description Largely follows from fastenloc
#'   https://github.com/xqwen/fastenloc but uses `susieR` fitted objects as
#'   input to estimate prior for use with `coloc` package (coloc v5, aka
#'   SuSiE-coloc). The main differences are 1) now enrichment is based on all
#'   QTL variants whether or not they are inside signal clusters; 2) Causal QTL
#'   are sampled from SuSiE single effects, not signal clusters; 3) Allow a
#'   variant to be QTL for not only multiple conditions (eg cell types) but also
#'   multiple regions (eg genes). Other minor improvements include 1) Make GSL
#'   RNG thread-safe; 2) Release memory from QTL binary annotation samples
#'   immediately after they are used.
#' @details Uses output of \code{\link[susieR]{susie}} from the \code{susieR}
#'   package.
#'
#' @param gwasPip This is a vector of GWAS PIP, genome-wide.
#' @param susieQtlRegions This is a list of SuSiE fitted objects per QTL unit
#'   analyzed
#' @param numGwas This parameter is highly important if GWAS input does not
#'   contain all SNPs interrogated (e.g., in some cases, only fine-mapped geomic
#'   regions are included). Then users must pick a value of total_variants and
#'   estimate piGwas beforehand by: sum(gwasPip$pip)/numGwas. If numGwas is
#'   null, piGwas would be sum(gwasPip$pip)/total_variants.
#' @param piQtl This parameter can be safely left to default if your input QTL
#'   data has enough regions to estimate it.
#' @param lambda Similar to the shrinkage parameter used in ridge regression. It
#'   takes any non-negative value and shrinks the enrichment estimate towards 0.
#'   When it is set to 0, no shrinkage will be applied. A large value indicates
#'   strong shrinkage. The default value is set to 1.0.
#' @param impN Rounds of multiple imputation to draw QTL from, default is 25.
#' @param numThreads Number of Simultaneous running CPU threads for multiple
#'   imputation, default is 1.
#' @param alignNames Logical; when TRUE (default) QTL pip names are aligned to
#'   the GWAS variant-naming convention via \code{matchVariants}. Set FALSE when
#'   the caller has already aligned them (e.g. \code{qtlEnrichmentPipeline}
#'   aligns each QTL tuple once against the union GWAS panel rather than
#'   re-aligning per GWAS study); only the cheap per-study unmatched set is then
#'   recomputed, skipping the costly \code{harmonizeAlleles} pass.
#' @param doubleShrinkage Logical. Apply the double-shrinkage correction to the
#'   enrichment estimate. Default \code{FALSE}.
#' @param besselCorrection Logical. Apply Bessel's correction when estimating
#'   the sampling variance. Default \code{TRUE}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#' @param seed Integer or \code{NULL}. Base random seed for the multiple-
#'   imputation sampler; each imputation round derives its own seed from it, so
#'   a fixed \code{seed} gives reproducible results. \code{NULL} (default) draws
#'   a nondeterministic seed.
#' @return A named list of enrichment parameter estimates, carrying the
#'   fields the C++ estimator writes -- \code{Intercept},
#'   \code{sd (intercept)}, \code{Enrichment (no shrinkage)},
#'   \code{Enrichment (w/ shrinkage)}, \code{sd (no shrinkage)},
#'   \code{sd (w/ shrinkage)}, \code{Alternative (coloc) p1} / \code{p2} /
#'   \code{p12} and \code{Effective MI rounds} -- plus
#'   \code{unused_xqtl_variants}, the QTL variants of each region that no GWAS
#'   variant matched. The names are upstream fastenloc's;
#'   \code{\link{qtlEnrichmentPipeline}} is what renames them to a tidy table.
#'
#' @examples
#'
#' # Simulate fake data for gwasPip
#' nGwasPip <- 1000
#' gwasPip <- runif(nGwasPip)
#' names(gwasPip) <- paste0("snp", 1:nGwasPip)
#' # Simulate fake data for a single SuSiEFit object
#' simulateSusiefit <- function(n, p) {
#'   pip <- runif(n)
#'   names(pip) <- paste0("snp", 1:n)
#'   alpha <- t(matrix(runif(n * p), nrow = n))
#'   alpha <- t(apply(alpha, 1, function(row) row / sum(row)))
#'   list(
#'     pip = pip,
#'     alpha = alpha,
#'     prior_variance = runif(p)
#'   )
#' }
#' # Simulate multiple SuSiEFit objects
#' nSusieFits <- 2
#' susieFits <- replicate(
#'   nSusieFits, simulateSusiefit(nGwasPip, 10), simplify = FALSE)
#' # Add these fits to a list, providing names to each element
#' names(susieFits) <- paste0("fit", seq_along(susieFits))
#' # Set other parameters
#' impN <- 10
#' lambda <- 1
#' numThreads <- 1
#' library(pecotmr)
#' en <- qtlEnrichment(
#'   gwasPip, susieFits, lambda = lambda, impN = impN,
#'   numThreads = numThreads)
#'
#' @seealso \code{\link[susieR]{susie}}
#' @useDynLib pecotmr, .registration = TRUE
#' @export
#'
qtlEnrichment <- function(
    gwasPip,
    susieQtlRegions,
    numGwas = NULL,
    piQtl = NULL,
    lambda = 1.0,
    impN = 25,
    doubleShrinkage = FALSE,
    besselCorrection = TRUE,
    numThreads = 1,
    verbose = TRUE,
    alignNames = TRUE,
    seed = NULL
) {
    piGwas <- .enrEstimatePiGwas(gwasPip, numGwas, verbose)
    piQtl <- .enrEstimatePiQtl(susieQtlRegions, piQtl, verbose)
    .enrValidatePi(piGwas, piQtl)
    .enrValidateNames(gwasPip, susieQtlRegions)
    # Align each region's pip names to the GWAS convention + record unmatched.
    aligned <- .enrAlignPipNames(susieQtlRegions, gwasPip, alignNames)
    unmatchedVariants <- map(aligned, "unmatched_variants")
    susieQtlRegions <- map(aligned, .enrStripUnmatched)
    # cpp11 requires exact integer types for int parameters.
    en <- qtlEnrichmentRcpp(
        rGwasPip = gwasPip,
        rQtlSusieFit = susieQtlRegions,
        piGwas = piGwas,
        piQtl = piQtl,
        ImpN = as.integer(impN),
        shrinkageLambda = lambda,
        doubleShrinkage = doubleShrinkage,
        besselCorrection = besselCorrection,
        numThreads = as.integer(numThreads),
        seed = if (is.null(seed)) NULL else as.integer(seed)
    )
    en$unused_xqtl_variants <- unmatchedVariants
    en
}

# piGwas = sum(gwasPip) / numGwas (estimated from the data, with a warning, when
# numGwas is absent).
# @noRd
.enrEstimatePiGwas <- function(gwasPip, numGwas, verbose) {
    if (!is.null(numGwas)) {
        return(sum(gwasPip) / numGwas)
    }
    msg <- glue(
        "numGwas is not provided. Estimating piGwas from the data. Note ",
        "that this estimate may be biased if the input gwasPip does not ",
        "contain genome-wide variants."
    )
    warn(msg)
    piGwas <- sum(gwasPip) / length(gwasPip)
    if (verbose) {
        piGwasR <- round(piGwas, 5)
        msg <- glue("Estimated piGwas: {piGwasR}\n", .trim = FALSE)
        inform(msg)
    }
    piGwas
}

# piQtl = total signal / total variants across regions (estimated, with a
# warning, when piQtl is absent).
# @noRd
.enrEstimatePiQtl <- function(susieQtlRegions, piQtl, verbose) {
    if (!is.null(piQtl)) {
        return(piQtl)
    }
    msg <- glue(
        "piQtl is not provided. Estimating piQtl from the data. Note that ",
        "this estimate may be biased if either 1) the input susieQtlRegions ",
        "does not have enough data, or 2) the single effects only include ",
        "variables inside of credible sets or signal clusters."
    )
    warn(msg)
    allPips <- unlist(map(susieQtlRegions, "pip"))
    piQtl <- sum(allPips) / length(allPips)
    if (verbose) {
        piQtlR <- round(piQtl, 5)
        msg <- glue("Estimated piQtl: {piQtlR}\n", .trim = FALSE)
        inform(msg)
    }
    piQtl
}

# Neither prior probability may be zero.
# @noRd
.enrValidatePi <- function(piGwas, piQtl) {
    if (piGwas == 0) {
        msg <- glue(
            "Cannot perform enrichment analysis. No association signal found ",
            "in GWAS data."
        )
        abort(msg)
    }
    if (piQtl == 0) {
        msg <- glue(
            "Cannot perform enrichment analysis. No QTL associated with the ",
            "molecular phenotype."
        )
        abort(msg)
    }
    invisible(NULL)
}

# Both gwasPip and every region's pip must carry variant names.
# @noRd
.enrValidateNames <- function(gwasPip, susieQtlRegions) {
    if (is.null(names(gwasPip))) {
        msg <- glue(
            "Variant names are missing in gwasPip. Please provide named ",
            "gwasPip data."
        )
        abort(msg)
    }
    if (!all(map_lgl(susieQtlRegions, .enrHasPipNames))) {
        msg <- glue(
            "Variant names are missing in susieQtlRegions$pip. Please provide ",
            "susieQtlRegions with named pip data."
        )
        abort(msg)
    }
    invisible(NULL)
}

# TRUE when a region's pip vector is named.
# @noRd
.enrHasPipNames <- function(x) {
    !is.null(names(x$pip))
}

# Align each region's pip names to gwasPip (relabel + record unmatched), or --
# when the caller already aligned (alignNames = FALSE) -- only recompute the
# cheap per-study unmatched set.
# @noRd
.enrAlignPipNames <- function(susieQtlRegions, gwasPip, alignNames) {
    if (alignNames) {
        return(map(susieQtlRegions, .enrAlignRegionByMatch, gwasPip = gwasPip))
    }
    map(susieQtlRegions, .enrMarkUnmatched, gwasNameSet = names(gwasPip))
}

# Relabel a region's matched pip names to the GWAS convention via the shared
# matcher, recording the unmatched variant names.
# @noRd
.enrAlignRegionByMatch <- function(x, gwasPip) {
    mm <- matchVariants(names(x$pip), names(gwasPip))
    nm <- names(x$pip)
    nm[mm$idxA] <- names(gwasPip)[mm$idxB]
    names(x$pip) <- nm
    unmatchedIdx <- setdiff(seq_along(x$pip), mm$idxA)
    if (length(unmatchedIdx) > 0) {
        x$unmatched_variants <- names(x$pip)[unmatchedIdx]
    }
    x
}

# Record the region's variants absent from the GWAS name set (cheap membership
# test; names already aligned by the caller).
# @noRd
.enrMarkUnmatched <- function(x, gwasNameSet) {
    unmatchedIdx <- which(!is_in(names(x$pip), gwasNameSet))
    if (length(unmatchedIdx) > 0) {
        x$unmatched_variants <- names(x$pip)[unmatchedIdx]
    }
    x
}

# Drop the transient unmatched_variants field from a region.
# @noRd
.enrStripUnmatched <- function(x) {
    x$unmatched_variants <- NULL
    x
}

# Relabel one region's matched pip names to the union GWAS panel (unmatched
# names kept as-is).
# @noRd
.enrAlignRegion <- function(x, unionGwasNames) {
    if (!is.null(names(x$pip)) && length(unionGwasNames) > 0L) {
        mm <- matchVariants(names(x$pip), unionGwasNames)
        nm <- names(x$pip)
        nm[mm$idxA] <- unionGwasNames[mm$idxB]
        names(x$pip) <- nm
    }
    x
}
