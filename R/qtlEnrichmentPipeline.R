#' @title QTL Enrichment Pipeline (Genome-Wide)
#' @description Genome-wide pipeline that computes per-pair (GWAS study, QTL
#'   context) enrichment estimates by passing the GWAS PIP vector and the QTL
#'   credible-set posteriors to \code{\link{qtlEnrichment}}. The returned table
#'   feeds \code{\link{colocPipeline}} via its \code{enrichment} argument.
#'
#' \strong{Not gene-parallelisable}: the enrichment estimator runs over the full
#' genome of GWAS PIPs and the full collection of QTL fits at once.
#'
#' @section Inputs:
#' \itemize{
#'   \item \code{gwasFineMappingResult}: a genome-wide
#'     \code{\link{GwasFineMappingResult}} (one row per (study, LD
#'     block) tuple). Each entry's \code{FineMappingEntry$trimmedFit}
#'     must carry a \code{pip} vector.
#'   \item \code{qtlFineMappingResult}: the genome-wide
#'     \code{\link{QtlFineMappingResult}}. Each entry's
#'     \code{trimmedFit} must carry \code{alpha}, \code{pip}, and
#'     prior-variance fields (\code{V}).
#' }
#'
#' @section LD-sketch identity check: The GWAS \code{FineMappingResultBase} must
#'   have a non-NULL \code{ldSketch} (RSS-derived). If the QTL FMR also has a
#'   non-NULL \code{ldSketch}, the two must match exactly. When the QTL FMR's
#'   \code{ldSketch} is NULL (individual-level QTL fit), validation is skipped
#'   on the QTL side.
#'
#' @param gwasFineMappingResult See above.
#' @param qtlFineMappingResult See above.
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
#' @return A tibble with one row per (gwasStudy, qtlStudy, qtlContext)
#'   triple and columns \code{gwasStudy}, \code{qtlStudy}, \code{qtlContext},
#'   \code{enrichment}, \code{enrichmentSe}, \code{enrichmentLogOdds}, plus any
#'   extras the underlying estimator emits. Suitable as the \code{enrichment}
#'   argument to \code{\link{colocPipeline}} (which joins on the same triple).
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
    results <- list_flatten(map(p$gwasStudies, .enrScoreGwasStudy, p = p))
    .enrAssemble(results)
}

# Validate the input classes + LD-sketch presence / identity.
# @noRd
.enrValidateInputs <- function(gwasFineMappingResult, qtlFineMappingResult) {
    if (!methods::is(gwasFineMappingResult, "GwasFineMappingResult")) {
        abort("`gwasFineMappingResult` must be a GwasFineMappingResult.")
    }
    if (!methods::is(qtlFineMappingResult, "QtlFineMappingResult")) {
        abort("`qtlFineMappingResult` must be a QtlFineMappingResult.")
    }
    gwasLd <- getLdSketch(gwasFineMappingResult)
    if (is.null(gwasLd)) {
        msg <- glue(
            "qtlEnrichmentPipeline: the GWAS FineMappingResult must have a ",
            "non-NULL ldSketch (it should be RSS-derived)."
        )
        abort(msg)
    }
    .colocRequireMatchingLdSketches(getLdSketch(qtlFineMappingResult), gwasLd)
    invisible(NULL)
}

# Hoist the GWAS-study-independent work out of the double loop: per-study GWAS
# PIP vectors, the union variant-name panel, per-tuple QTL regions, and each
# tuple's one-time alignment to the union panel (errors captured as values).
# @noRd
.enrPrepare <- function(p) {
    p$gwasStudies <- unique(as.character(p$gwasFineMappingResult$study))
    # Iterate the joint (study, context) key: context alone would silently merge
    # two studies sharing a context label, giving wrong enrichment estimates.
    p$qtlTuples <- distinct(tibble(
        qtlStudy = as.character(p$qtlFineMappingResult$study),
        qtlContext = as.character(p$qtlFineMappingResult$context)
    ))
    if (length(p$gwasStudies) == 0L || nrow(p$qtlTuples) == 0L) {
        msg <- glue(
            "qtlEnrichmentPipeline: no (gwasStudy, qtlStudy, qtlContext) ",
            "triples to compute (one of the inputs has zero rows)."
        )
        abort(msg)
    }
    p$gwasPipByStudy <- .enrGwasPipByStudy(
        p$gwasFineMappingResult,
        p$gwasStudies
    )
    unionGwasNames <- unique(unlist(
        map(p$gwasPipByStudy, names),
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

# Per-study genome-wide GWAS PIP vectors (named by variant id). A loop rather
# than map because the varying study is `.enrBuildGwasPipVector`'s second arg.
# @noRd
.enrGwasPipByStudy <- function(gwasFineMappingResult, gwasStudies) {
    out <- vector("list", length(gwasStudies))
    for (i in seq_along(gwasStudies)) {
        out[[i]] <- .enrBuildGwasPipVector(
            gwasFineMappingResult,
            gwasStudies[[i]]
        )
    }
    set_names(out, gwasStudies)
}

# QTL SuSiE region list for the k-th (study, context) tuple.
# @noRd
.enrQtlRegionsForRow <- function(k, qtlTuples, fmr) {
    .enrBuildQtlRegionsList(
        fmr,
        qtlTuples$qtlStudy[[k]],
        qtlTuples$qtlContext[[k]]
    )
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

# Score one GWAS study against every QTL tuple -> enrichment records (empty when
# the study has no usable PIPs).
# @noRd
.enrScoreGwasStudy <- function(gStudy, p) {
    gwasPip <- p$gwasPipByStudy[[gStudy]]
    if (length(gwasPip) == 0L) {
        msg <- glue(
            "qtlEnrichmentPipeline: no usable PIPs for ",
            "gwasStudy='{gStudy}'; skipping."
        )
        warn(msg)
        return(list())
    }
    compact(map(
        seq_len(nrow(p$qtlTuples)),
        .enrScoreTuple,
        gStudy = gStudy,
        gwasPip = gwasPip,
        p = p
    ))
}

# Score one (gwasStudy, qtl tuple) pair -> an enrichment record, or NULL when
# the tuple has no regions or qtlEnrichment fails.
# @noRd
.enrScoreTuple <- function(k, gStudy, gwasPip, p) {
    qStudy <- p$qtlTuples$qtlStudy[[k]]
    qContext <- p$qtlTuples$qtlContext[[k]]
    if (length(p$qtlRegionsByTuple[[k]]) == 0L) {
        msg <- glue(
            "qtlEnrichmentPipeline: no usable QTL regions for ",
            "(qtlStudy='{qStudy}', qtlContext='{qContext}'); skipping."
        )
        warn(msg)
        return(NULL)
    }
    enr <- .enrRunEnrichment(gStudy, gwasPip, k, qStudy, qContext, p)
    if (is.null(enr)) {
        return(NULL)
    }
    row <- .enrFlattenEnrichment(enr)
    row$gwasStudy <- gStudy
    row$qtlStudy <- qStudy
    row$qtlContext <- qContext
    row
}

# Run qtlEnrichment for a pair (with the pre-aligned regions), warning + NULL on
# failure. alignNames = FALSE reuses the shared per-tuple alignment.
# @noRd
.enrRunEnrichment <- function(gStudy, gwasPip, k, qStudy, qContext, p) {
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
                "(gwasStudy='{gStudy}', qtlStudy='{qStudy}', ",
                "qtlContext='{qContext}'): {eMsg}"
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
    idCols <- c("gwasStudy", "qtlStudy", "qtlContext")
    select(out, all_of(idCols), everything())
}

# The empty enrichment result table.
# @noRd
.enrEmptyResult <- function() {
    tibble(
        gwasStudy = character(0),
        qtlStudy = character(0),
        qtlContext = character(0),
        enrichment = numeric(0),
        enrichmentSe = numeric(0),
        enrichmentLogOdds = numeric(0)
    )
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

# Build a named GWAS PIP vector for one study. Walks every row of the
# GwasFineMappingResult tagged with that study, extracts the per-row
# pip from each FineMappingEntry, and concatenates with variant-id
# names. Errors if any single variant appears with conflicting PIP
# values across rows.
#' @importFrom dplyr add_count
#' @noRd
.enrBuildGwasPipVector <- function(gwasFmr, gStudy) {
    idx <- which(as.character(gwasFmr$study) == gStudy)
    if (length(idx) == 0L) {
        return(numeric(0))
    }
    pieces <- list()
    for (i in idx) {
        entry <- gwasFmr$entry[[i]]
        fit <- getSusieFit(entry)
        if (is.null(fit) || is.null(fit$pip)) {
            next
        }
        pip <- as.numeric(fit$pip)
        ids <- if (!is.null(names(fit$pip))) {
            names(fit$pip)
        } else {
            getVariantIds(entry)
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

# Build the per-(qtlStudy, qtlContext) list of region fits in the shape
# that qtlEnrichment expects: list(d) where each d carries
# alpha, pip, prior_variance (V). Filters on BOTH study and context so
# entries from different studies that happen to share a context label
# are not pooled into one enrichment estimate.
# @noRd
.enrBuildQtlRegionsList <- function(qtlFmr, qStudy, qContext) {
    idx <- which(
        as.character(qtlFmr$study) == qStudy &
            as.character(qtlFmr$context) == qContext
    )
    if (length(idx) == 0L) {
        return(list())
    }
    out <- list()
    for (i in idx) {
        entry <- qtlFmr$entry[[i]]
        fit <- getSusieFit(entry)
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
            names(fit$pip) <- getVariantIds(entry)
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

# Coerce qtlEnrichment's variable-shape output into a single-row
# named list with the canonical columns the caller documents. The
# underlying estimator returns either a list with named numeric scalars
# (enrichment, enrichmentSe, enrichmentLogOdds, ...) or a matrix/df --
# this helper handles both.
# @noRd
.enrFlattenEnrichment <- function(enr) {
    if (is.list(enr) && is.null(dim(enr))) {
        list(
            enrichment = .enrPickScalar("enrichment", enr),
            enrichmentSe = .enrPickScalar("enrichmentSe", enr),
            enrichmentLogOdds = .enrPickScalar("enrichmentLogOdds", enr)
        )
    } else if (is.matrix(enr) || is.data.frame(enr)) {
        df <- as_tibble(enr, .name_repair = "minimal")
        if (nrow(df) == 0L) {
            list(
                enrichment = NA_real_,
                enrichmentSe = NA_real_,
                enrichmentLogOdds = NA_real_
            )
        } else {
            list(
                enrichment = .enrPickColumn(df, c("enrichment", "Enrichment")),
                enrichmentSe = .enrPickColumn(
                    df,
                    c("enrichmentSe", "se", "stderr")
                ),
                enrichmentLogOdds = .enrPickColumn(
                    df,
                    c("enrichmentLogOdds", "logOdds", "log_odds")
                )
            )
        }
    } else {
        list(
            enrichment = NA_real_,
            enrichmentSe = NA_real_,
            enrichmentLogOdds = NA_real_
        )
    }
}

# @noRd
.enrPickColumn <- function(df, candidates) {
    hit <- intersect(candidates, colnames(df))
    if (length(hit) == 0L) {
        return(NA_real_)
    }
    as.numeric(df[[hit[[1L]]]][[1L]])
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
#' @return A list of enrichment parameter estimates
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
    en <- list(qtlEnrichmentRcpp(
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
    ))
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
