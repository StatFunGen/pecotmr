#' @include qtlSumStats.R
#' @title ColocBoost multi-trait colocalization pipeline (S4)
#' @description Protocol-level multi-trait colocalization analysis using
#'   \pkg{colocboost}. Dispatches on the QTL input type:
#'   \itemize{
#'     \item \code{QtlDataset} -- single-study, individual-level
#'           multi-context data. Per-context residualized X / Y are
#'           extracted from the dataset (filtering knobs on the
#'           constructor apply lazily inside the accessors).
#'     \item \code{QtlSumStats} -- summary-statistic-only QTL data with a
#'           shared LD reference (\code{ldSketch}). Must already have
#'           been passed through \code{\link{summaryStatsQc}} (the
#'           pipeline rejects inputs whose \code{getQcInfo()} is empty).
#'     \item \code{MultiStudyQtlDataset} -- a mixture of one or more
#'           individual-level \code{QtlDataset} studies and an optional
#'           \code{QtlSumStats} collection.
#'   }
#'   GWAS is optional and always passed separately as a
#'   \code{GwasSumStats} object (must also be QC'd).
#'
#' \code{colocboostPipeline} does \strong{not} accept a \code{FineMappingResult}
#' for either side; colocboost has its own variable-selection algorithm.
#'
#' @section QC contract:
#'   \itemize{
#'     \item Individual-level QC (MAF / MAC / X-variance / per-sample
#'           missingness, sample / variant restrictions) lives on the
#'           \code{QtlDataset} constructor and is applied lazily inside
#'           \code{getGenotypes()} / \code{getResidualizedGenotypes()}.
#'           The pipeline does \emph{not} run a separate
#'           individual-level QC pass.
#'     \item \code{mafCutoff} / \code{macCutoff} / \code{imissCutoff} act
#'           at ANALYSIS time on the sumstat sides, measured against the
#'           \strong{LD reference panel}: a variant whose panel genotypes fall
#'           below the cutoffs is dropped before the LD matrix is built. They
#'           apply to a \code{QtlSumStats} QTL side and to a
#'           \code{GwasSumStats} GWAS side (so they still bite when the QTL
#'           side is individual-level). Defaults (\code{0}, \code{0},
#'           \code{1}) filter nothing. Unlike
#'           \code{summaryStatsQc(imputeOpts = ...)}, which only bounds what
#'           RAISS will impute, these discard \emph{observed} variants.
#'     \item All summary-statistic QC (variant filters, harmonization
#'           against the \code{ldSketch}, LD-mismatch detection, RAISS
#'           imputation, etc.) lives in
#'           \code{\link{summaryStatsQc}}. The pipeline rejects any
#'           \code{QtlSumStats} or \code{GwasSumStats} where
#'           \code{length(getQcInfo(x)) == 0L}.
#'   }
#'
#' @section Analysis variants:
#'   \itemize{
#'     \item \code{xqtlColoc} (default \code{TRUE}): run a colocboost
#'           model over the QTL contexts only (individual-level inputs).
#'     \item \code{jointGwas} (default \code{FALSE}): run a non-focal
#'           colocboost model that combines all QTL contexts/studies
#'           with the supplied \code{gwasSumStats} studies.
#'     \item \code{separateGwas} (default \code{FALSE}): run one focal
#'           colocboost model per GWAS study, where the GWAS is the
#'           focal outcome.
#'   }
#'
#' @param qtlData One of \code{QtlDataset}, \code{QtlSumStats}, or
#'   \code{MultiStudyQtlDataset}.
#' @param gwasSumStats Optional \code{GwasSumStats} with the GWAS studies to
#'   colocalize against. \code{NULL} to skip GWAS colocalization.
#' @param contexts Optional character vector of context names to restrict the
#'   individual-level / QtlSumStats QTL analysis to. When \code{NULL} (default),
#'   every context present is used.
#' @param traitId Optional character vector of trait identifiers to restrict the
#'   analysis to. When supplied with an individual-level \code{QtlDataset}
#'   input, \code{cisWindow} is required (passed to
#'   \code{getResidualizedGenotypes} / \code{getPhenotypes} for the
#'   variant-window selection).
#' @param region Optional single-range \code{GRanges} describing the analysis
#'   window. Mutually exclusive with \code{traitId} (see the \code{QtlDataset}
#'   accessors).
#' @param cisWindow Optional cis window in basepairs; required with
#'   \code{traitId}, optional with \code{region}.
#' @param samples Optional character vector of sample IDs to restrict the
#'   analysis to; \code{NULL} (default) uses all samples.
#' @param focalTrait Optional trait name; when supplied and present in the
#'   assembled outcome list, the colocboost xQTL-only run uses it as the focal
#'   outcome.
#' @param xqtlColoc,jointGwas,separateGwas Logical flags selecting which
#'   colocboost variants to run.
#' @param pipCutoffToSkip Individual-level pre-filter (ports the legacy
#'   \code{pip_cutoff_to_skip_ind}). Scalar (applied to every context) or a
#'   context-named numeric vector. For each context, every outcome is fit with a
#'   single-effect SuSiE (\code{L = 1}) and dropped unless some variant's PIP
#'   exceeds the cutoff; a context with no surviving outcome is skipped.
#'   \code{0} (default) disables it; a negative value uses \code{3 /
#'   n_variants}. (Summary-statistic skipping is handled upstream by
#'   \code{\link{summaryStatsQc}}'s own \code{pipCutoffToSkip}.)
#' @param absZCutoffToSkip,bfCutoffToSkip,logBfCutoffToSkip Alternative
#'   individual-level pre-filter metrics used in place of
#'   \code{pipCutoffToSkip}: drop an outcome unless its maximum marginal
#'   \code{|z|} (\code{absZCutoffToSkip}), or its maximum per-variant
#'   single-effect Bayes factor (\code{bfCutoffToSkip}) / log Bayes factor
#'   (\code{logBfCutoffToSkip}) from the \code{L = 1} fit, exceeds the cutoff.
#'   Scalars, each defaulting to 0 (off). Exactly one screening metric may be
#'   enabled: setting any of these requires \code{pipCutoffToSkip = 0}.
#' @param mafCutoff,macCutoff,imissCutoff Analysis-time filters applied to the
#'   summary-statistic sides against the LD reference panel: a variant whose
#'   panel genotypes fall below the cutoffs is dropped before the LD matrix is
#'   built. Defaults (\code{0}, \code{0}, \code{1}) filter nothing. See the
#'   Details section for how these relate to the individual-level cutoffs
#'   recorded on a \code{QtlDataset}.
#' @param alleleFlip Logical, default \code{TRUE}. When TRUE, harmonize variants
#'   across the individual X, sumstats, and LD by (chrom, pos) with ref/alt
#'   swaps recognized (flipping z / residualized dosage / LD to a shared
#'   coding); when FALSE, match on exact alleles only (names-only), so a ref/alt
#'   swap is treated as a distinct variant.
#' @param ... Additional arguments forwarded to
#'   \code{\link[colocboost]{colocboost}} (e.g., \code{M}, \code{L},
#'   \code{output_level}).
#' @return A \code{\linkS4class{ColocBoostResult}}: one element per
#'   confidence set (CoS) across every analysis that ran, holding that set's
#'   member variants with their \code{vcp}. The \code{analysis} column marks
#'   which run each set came from (\code{xqtl_coloc}, \code{joint_gwas} or
#'   \code{separate_gwas}), and \code{gwasStudy} distinguishes the per-study
#'   \code{separate_gwas} runs. Outcome-specific (uncolocalized) sets are
#'   included with \code{isColocalized = FALSE}.
#'
#'   Project it with \code{\link{getColocPairs}} (also available as
#'   \code{as.data.frame}), \code{\link{getColocVariants}} or
#'   \code{\link{getColocBoostOutcomes}}; timings are on
#'   \code{\link{getComputingTime}} and the region-wide marginal
#'   probabilities on \code{\link{getRegionVcp}}.
#' @name colocboostPipeline
#' @importFrom methods is setGeneric setMethod
#' @importFrom S4Vectors mcols
#' @importFrom GenomicRanges seqnames start end GRanges
#' @importFrom IRanges IRanges
#' @export
NULL

# =============================================================================
# Generic
# =============================================================================

#' @rdname colocboostPipeline
#' @examples
#' data(qtlDatasetExample)
#' colocboostPipeline(qtlDatasetExample, xqtlColoc = TRUE)
#' @export
setGeneric("colocboostPipeline", function(qtlData, gwasSumStats = NULL, ...) {
    standardGeneric("colocboostPipeline")
})

# =============================================================================
# Helpers (private)
# =============================================================================

# Run colocboost() with tryCatch + timing.
.cbRun <- function(label, args) {
    if (!requireNamespace("colocboost", quietly = TRUE)) {
        # nocov start
        abort("The colocboost package is required for colocboostPipeline().")
        # nocov end
    }
    t1 <- Sys.time()
    args <- compact(args)
    res <- tryCatch(
        exec(colocboost::colocboost, !!!args),
        error = function(e) {
            eMsg <- conditionMessage(e)
            msg <- glue("{label} failed: {eMsg}")
            inform(msg)
            NULL
        }
    )
    list(result = res, time = Sys.time() - t1)
}

# Build the LD / X_ref slot of the colocboost call from a list of LD
# matrices. When any matrix is non-square it is treated as a samples x
# variants genotype reference and routed to X_ref; otherwise routed to LD.
.cbBuildLdArgs <- function(ldList) {
    ldList <- compact(ldList)
    if (length(ldList) == 0L) {
        return(list())
    }
    isGeno <- any(map_lgl(ldList, .cbIsNonSquare))
    if (isGeno) list(X_ref = ldList) else list(LD = ldList)
}

# Reject SumStats objects that have not been passed through
# summaryStatsQc(). Both QtlSumStats and GwasSumStats expose getQcInfo();
# an empty list (the constructor default) signals "no QC run".
.cbRequireSumStatsQc <- function(x, what) {
    if (is.null(x)) {
        return(invisible(NULL))
    }
    if (length(getQcInfo(x)) == 0L) {
        msg <- glue(
            "{what} must be passed through summaryStatsQc() before ",
            "reaching colocboostPipeline (getQcInfo() returned an empty ",
            "list). Call summaryStatsQc(x, ...) and pass the result."
        )
        abort(msg)
    }
    invisible(NULL)
}

# Resolve the per-context screen spec. The `pipCutoffToSkip` channel carries
# EITHER a resolved screen object (list(metric, cutoff) for absZ/bf/logBf/pip --
# applies uniformly to every context) OR the legacy PIP cutoff as a scalar (all
# contexts) or a named vector keyed by context. Default 0 (no screen).
.cbResolveCutoff <- function(pipCutoffToSkip, ctx) {
    if (is.null(pipCutoffToSkip) || length(pipCutoffToSkip) == 0L) {
        return(0)
    }
    if (is.list(pipCutoffToSkip)) {
        return(pipCutoffToSkip)
    } # uniform screen object
    if (!is.null(names(pipCutoffToSkip))) {
        if (is_in(ctx, names(pipCutoffToSkip))) {
            return(pipCutoffToSkip[[ctx]])
        }
        return(0)
    }
    pipCutoffToSkip[[1L]]
}

# Combine the four colocboost screen cutoffs into a single spec to thread
# through the pipCutoffToSkip channel: a resolved screen object when a new
# metric (absZ / bf / logBf) is set, otherwise the (possibly context-named)
# legacy pipCutoffToSkip. Enforces one screening metric at a time.
.cbScreenSpec <- function(
    pipCutoffToSkip,
    absZCutoffToSkip,
    bfCutoffToSkip,
    logBfCutoffToSkip
) {
    newScreen <- .resolveScreenMetric(
        0,
        absZCutoffToSkip,
        bfCutoffToSkip,
        logBfCutoffToSkip
    )
    pipOn <- !is.null(pipCutoffToSkip) &&
        length(pipCutoffToSkip) > 0L &&
        any(as.numeric(pipCutoffToSkip) != 0, na.rm = TRUE)
    if (!is.null(newScreen) && pipOn) {
        msg <- glue(
            "colocboostPipeline: only one signal screen may be enabled ",
            "at a time; unset pipCutoffToSkip to use absZCutoffToSkip / ",
            "bfCutoffToSkip / logBfCutoffToSkip."
        )
        abort(msg)
    }
    if (!is.null(newScreen)) newScreen else pipCutoffToSkip
}

# Per-outcome single-trait skip (ports the legacy qc_individual_data
# pip_cutoff_to_skip): for each outcome column of Y, fit a single-effect
# SuSiE (L = 1, max_iter = 100) on (X, Y[, j]) and keep the outcome only if
# any variant's PIP exceeds the cutoff. A cutoff < 0 means 3 / n_variants.
# Returns the retained Y (NULL when no outcome clears the threshold).
.cbPipSkipOutcomes <- function(X, Y, spec) {
    if (is.null(.asScreen(spec))) {
        return(Y)
    }
    # Single-effect screen per outcome, sharing the L = 1 SuSiE pre-screen
    # (.fmSerScreen) with the fine-mapping pipeline. fallback = FALSE: an
    # outcome that cannot be screened (too few samples / fit failure) is
    # dropped.
    keep <- map_lgl(
        seq_len(ncol(Y)),
        .cbScreenOutcome,
        X = X,
        Y = Y,
        spec = spec
    )
    if (!any(keep)) {
        return(NULL)
    }
    Y[, keep, drop = FALSE]
}

# Materialise an individual-level QtlDataset into the colocboost
# (X, Y, dict_YX, outcome_names) bundle. Each context becomes one X /
# Y pair; the YA matrices are split into single-trait columns and
# dict_YX maps each split column back to its X. Returns NULL when no
# context survives selection. pipCutoffToSkip (scalar or context-named
# vector) optionally drops weak-signal outcomes / contexts up front.
.cbIndividualBundle <- function(
    qd,
    contexts = NULL,
    traitId = NULL,
    region = NULL,
    cisWindow = NULL,
    samples = NULL,
    pipCutoffToSkip = 0
) {
    contexts <- .cbResolveBundleContexts(qd, contexts)
    p <- list(
        qd = qd,
        traitId = traitId,
        region = region,
        cisWindow = cisWindow,
        samples = samples,
        pipCutoffToSkip = pipCutoffToSkip
    )
    built <- compact(map(contexts, .cbBuildContextXY, p = p))
    if (length(built) == 0L) {
        return(NULL)
    }
    XperCtx <- set_names(map(built, "X"), map_chr(built, "ctx"))
    YperCtx <- set_names(map(built, "Y"), map_chr(built, "ctx"))
    dedup <- .cbDedupX(XperCtx)
    split <- .cbSplitY(YperCtx, dedup$xMatch)
    outcomeInfo <- split$outcomeInfo
    outcomeInfo$study <- getStudy(qd)
    outcomeInfo$dataForm <- "individual"
    list(
        X = dedup$uniqueX,
        Y = split$YSplit,
        dict_YX = split$dict,
        outcomeNames = names(split$YSplit),
        outcomeInfo = outcomeInfo
    )
}

# Resolve the contexts to bundle: all of them when unspecified, else the
# supplied set validated against the dataset.
# @noRd
.cbResolveBundleContexts <- function(qd, contexts) {
    if (is.null(contexts) || length(contexts) == 0L) {
        return(getContexts(qd))
    }
    available <- getContexts(qd)
    bad <- setdiff(contexts, available)
    if (length(bad) > 0L) {
        qdStudy <- getStudy(qd)
        badStr <- str_flatten(bad, ", ")
        availStr <- str_flatten(available, ", ")
        msg <- glue(
            "Unknown context(s) for QtlDataset '{qdStudy}': {badStr}. ",
            "Available: {availStr}"
        )
        abort(msg)
    }
    contexts
}

# Residualized (X, Y) for one context (sample-aligned + signal-screened), or
# NULL when the context should be skipped.
# @noRd
.cbBuildContextXY <- function(ctx, p) {
    Y <- .cbResidualizedY(p$qd, ctx, p$traitId, p$region)
    if (is.null(Y) || ncol(Y) == 0L) {
        return(NULL)
    }
    X <- .cbResidualizedX(
        p$qd,
        ctx,
        p$traitId,
        p$region,
        p$cisWindow,
        p$samples
    )
    if (is.null(X) || ncol(X) == 0L) {
        return(NULL)
    }
    # Canonicalize variant colnames (chr-prefix + separator, allele order
    # preserved; rsIDs passed through) so colocboost's name-based matching
    # aligns
    # them with the sumstat / LD ids and across studies. A genuine ref/alt swap
    # stays a distinct id -- names are aligned here, allele *coding* is not.
    colnames(X) <- normalizeVariantId(colnames(X))
    common <- intersect(rownames(X), rownames(Y))
    if (length(common) == 0L) {
        msg <- glue(
            "colocboostPipeline: skipping context '{ctx}' ",
            "(no samples shared between residualized X and Y)."
        )
        inform(msg)
        return(NULL)
    }
    X <- X[common, , drop = FALSE]
    Y <- Y[common, , drop = FALSE]
    Y <- .cbApplyScreen(X, Y, ctx, p$pipCutoffToSkip)
    if (is.null(Y)) {
        return(NULL)
    }
    list(ctx = ctx, X = X, Y = Y)
}

# Residualized phenotypes for one context (message + NULL on failure).
# @noRd
.cbResidualizedY <- function(qd, ctx, traitId, region) {
    tryCatch(
        getResidualizedPhenotypes(
            qd,
            contexts = ctx,
            traitId = traitId,
            region = region
        ),
        error = function(e) {
            eMsg <- conditionMessage(e)
            msg <- glue(
                "colocboostPipeline: skipping context '{ctx}' ",
                "(residualized phenotypes unavailable: {eMsg})."
            )
            inform(msg)
            NULL
        }
    )
}

# Residualized genotypes for one context (message + NULL on failure).
# @noRd
.cbResidualizedX <- function(qd, ctx, traitId, region, cisWindow, samples) {
    tryCatch(
        getResidualizedGenotypes(
            qd,
            contexts = ctx,
            traitId = traitId,
            region = region,
            cisWindow = cisWindow,
            samples = samples
        ),
        error = function(e) {
            eMsg <- conditionMessage(e)
            msg <- glue(
                "colocboostPipeline: skipping context '{ctx}' ",
                "(residualized genotypes unavailable: {eMsg})."
            )
            inform(msg)
            NULL
        }
    )
}

# Apply the optional signal screen to a context's outcomes; NULL when nothing
# clears it (context is then skipped).
# @noRd
.cbApplyScreen <- function(X, Y, ctx, pipCutoffToSkip) {
    cutoffCtx <- .cbResolveCutoff(pipCutoffToSkip, ctx)
    if (is.null(.asScreen(cutoffCtx))) {
        return(Y)
    }
    Y <- .cbPipSkipOutcomes(X, Y, cutoffCtx)
    if (is.null(Y) || ncol(Y) == 0L) {
        msg <- glue(
            "colocboostPipeline: skipping context '{ctx}' ",
            "(no outcome cleared the signal screen)."
        )
        inform(msg)
        return(NULL)
    }
    Y
}

# Deduplicate X matrices identical across contexts so dict_YX can fan out to a
# smaller X set. Returns list(uniqueX, xMatch) where xMatch[i] is the uniqueX
# index for the i-th context. (Sequential back-reference dedup -- kept as a
# loop.)
# @noRd
.cbDedupX <- function(XperCtx) {
    uniqueX <- list()
    xMatch <- integer(length(XperCtx))
    for (i in seq_along(XperCtx)) {
        matched <- names(uniqueX)[
            map_lgl(uniqueX, identical, XperCtx[[i]])
        ]
        if (length(matched) > 0L) {
            xMatch[[i]] <- match(matched[[1L]], names(uniqueX))
        } else {
            uniqueX[[names(XperCtx)[i]]] <- XperCtx[[i]]
            xMatch[[i]] <- length(uniqueX)
        }
    }
    list(uniqueX = uniqueX, xMatch = xMatch)
}

# Split each context's Y into single-trait columns (context-qualifying duplicate
# trait names) and build the dict_YX (Y-index, X-index) map. Returns
# list(YSplit, dict). (Sequential make.unique naming -- kept as a loop.)
# @noRd
.cbSplitY <- function(YperCtx, xMatch) {
    allTraitNames <- unlist(map(YperCtx, colnames), use.names = FALSE)
    dupTraits <- unique(allTraitNames[
        duplicated(allTraitNames) | duplicated(allTraitNames, fromLast = TRUE)
    ])
    YSplit <- list()
    dict <- matrix(integer(0), ncol = 2L)
    # The outcome NAME is a display label: context-qualified only when a trait
    # is ambiguous, then made unique. That makes it lossy -- a bare name could
    # be a trait seen in one context or a trait literally called that -- so the
    # (context, trait) it was minted from is recorded alongside it here, the
    # only point where both are still in hand.
    info <- list()
    for (i in seq_along(YperCtx)) {
        Y <- YperCtx[[i]]
        ctx <- names(YperCtx)[i]
        for (j in seq_len(ncol(Y))) {
            tname <- .cbTraitName(
                colnames(Y)[j],
                ctx,
                dupTraits,
                names(YSplit)
            )
            YSplit[[tname]] <- Y[, j, drop = FALSE]
            dict <- rbind(dict, c(length(YSplit), xMatch[[i]]))
            info[[length(info) + 1L]] <- tibble(
                name = tname,
                context = ctx,
                trait = as.character(colnames(Y)[j] %||% NA_character_)
            )
        }
    }
    colnames(dict) <- c("Y", "X")
    list(
        YSplit = YSplit,
        dict = dict,
        outcomeInfo = bind_rows(info)
    )
}

# Resolve a unique outcome name for a single trait column: default unnamed
# columns, context-qualify duplicates, and make.unique against existing names.
# @noRd
.cbTraitName <- function(tname, ctx, dupTraits, existing) {
    if (is.null(tname) || is.na(tname) || tname == "") {
        tname <- str_c("outcome", length(existing) + 1L)
    }
    if (is_in(tname, dupTraits)) {
        tname <- str_c(ctx, "_", tname)
    }
    if (is_in(tname, existing)) {
        tname <- make.unique(c(existing, tname))[length(existing) + 1L]
    }
    tname
}

# Build a single (sumstat data.frame, LD correlation matrix) pair from a
# QtlSumStats / GwasSumStats entry. Returns NULL when the entry has no
# variants overlapping the ldSketch panel.
.cbSumstatPair <- function(
    df,
    ldSketch,
    varY = NULL,
    nCase = NULL,
    nControl = NULL,
    cutoffs = NULL
) {
    if (is.null(df) || nrow(df) == 0L) {
        return(NULL)
    }
    df$variant_id <- .cbSumstatVariantIds(df)
    # Canonicalize ids (chr-prefix + separator, allele order preserved; rsIDs
    # passed through) so the sumstat `variant` column and the LD dimnames align
    # by name with the individual X colnames and across studies/sumstats.
    variantIds <- normalizeVariantId(df$variant_id)
    # Panel MAF / MAC / missingness cutoffs narrow what is asked of the LD
    # panel. `variantIds` itself stays full-length and aligned to `df`: the
    # row filter below indexes `df` by it, so shrinking it here would misalign
    # the two.
    wanted <- variantIds[.panelKeepMask(
        variantIds,
        ldSketch,
        cutoffs,
        ".cbSumstatPair"
    )]
    if (length(wanted) == 0L) {
        return(NULL)
    }
    # Use the shared `.ldFromSketch` helper in "drop" mode so missing-from-panel
    # variants are silently filtered (the colocboost path expects to operate
    # only on the overlap).
    R <- .ldFromSketch(
        ldSketch,
        wanted,
        label = ".cbSumstatPair",
        onMissing = "drop"
    )
    if (is.null(R)) {
        return(NULL)
    }
    keptIds <- attr(R, "keptVariantIds")
    attr(R, "keptVariantIds") <- NULL
    df <- filter(df, is_in(variantIds, keptIds))
    ss <- .cbSumstatFrame(df, keptIds, nCase, nControl, varY)
    list(sumstat = ss, LD = R, variantIds = keptIds)
}

# The colocboost sumstat frame for one entry, over the variants the LD panel
# kept. `var_y` is added only when supplied: colocboost treats its absence as
# "z-scale", so writing NA would change the model rather than say nothing.
# @noRd
.cbSumstatFrame <- function(df, keptIds, nCase, nControl, varY) {
    ss <- data.frame(
        z = df$z,
        n = .cbEffectiveN(df, nCase, nControl),
        variant = keptIds,
        stringsAsFactors = FALSE
    )
    if (!is.null(varY) && !is.na(varY)) {
        ss$var_y <- varY
    }
    ss
}

# Resolve a sumstat entry's variant ids, falling back to the canonical
# chr:pos:A2:A1 form when the entry had no SNP mcol (variant_id left NA).
# @noRd
.cbSumstatVariantIds <- function(df) {
    variantIds <- df$variant_id
    if (!anyNA(variantIds)) {
        return(variantIds)
    }
    formatVariantId(df$chrom, df$pos, df$A2, df$A1)
}

# Effective sample size: case/control -> 4 / (1/nCase + 1/nControl); otherwise
# the per-variant N (NA when neither is available).
# @noRd
.cbEffectiveN <- function(df, nCase, nControl) {
    okCC <- !is.null(nCase) &&
        !is.null(nControl) &&
        !is.na(nCase) &&
        !is.na(nControl) &&
        nCase > 0 &&
        nControl > 0
    if (okCC) {
        return(effectiveN(nCase, nControl))
    }
    if (!is.null(df$N)) {
        return(df$N)
    }
    NA_real_
}

# Build the colocboost sumstat / LD bundle from a QtlSumStats and an
# optional contexts / traitId filter. Returns a list keyed by sumstat
# study label, where each entry has (sumstat, LD).
.cbQtlSumStatsBundle <- function(
    ss,
    contexts = NULL,
    traitId = NULL,
    cutoffs = NULL
) {
    if (is.null(ss) || nrow(ss) == 0L) {
        return(list())
    }
    ldSketch <- getLdSketch(ss)
    keepRow <- rep(TRUE, nrow(ss))
    if (!is.null(contexts) && length(contexts) > 0L) {
        keepRow <- keepRow & is_in(as.character(ss$context), contexts)
    }
    if (!is.null(traitId) && length(traitId) > 0L) {
        keepRow <- keepRow & is_in(as.character(ss$trait), traitId)
    }
    if (!any(keepRow)) {
        return(list())
    }
    rows <- which(keepRow)

    bundle <- list()
    for (i in rows) {
        st <- as.character(ss$study)[[i]]
        ctx <- as.character(ss$context)[[i]]
        tr <- as.character(ss$trait)[[i]]
        label <- str_c(st, ctx, tr, sep = ":")
        pair <- .cbSumstatPair(
            df = getSumstatDf(
                ss,
                study = st,
                context = ctx,
                trait = tr,
                require = "Z"
            ),
            ldSketch = ldSketch,
            varY = if (is_in("varY", .tupleColumnNames(ss))) {
                ss$varY[[i]]
            } else {
                NA_real_
            },
            cutoffs = cutoffs
        )
        if (!is.null(pair)) bundle[[label]] <- pair
    }
    bundle
}

# Same as .cbQtlSumStatsBundle for a GwasSumStats collection, keyed by
# study label.
.cbGwasSumStatsBundle <- function(gws, cutoffs = NULL) {
    if (is.null(gws) || nrow(gws) == 0L) {
        return(list())
    }
    ldSketch <- getLdSketch(gws)
    bundle <- list()
    for (i in seq_len(nrow(gws))) {
        st <- as.character(gws$study)[[i]]
        pair <- .cbSumstatPair(
            df = getSumstatDf(gws, study = st, require = "Z"),
            ldSketch = ldSketch,
            varY = if (is_in("varY", .tupleColumnNames(gws))) {
                gws$varY[[i]]
            } else {
                NA_real_
            },
            nCase = if (is_in("nCase", .tupleColumnNames(gws))) {
                gws$nCase[[i]]
            } else {
                NA_real_
            },
            nControl = if (is_in("nControl", .tupleColumnNames(gws))) {
                gws$nControl[[i]]
            } else {
                NA_real_
            },
            cutoffs = cutoffs
        )
        if (!is.null(pair)) bundle[[st]] <- pair
    }
    bundle
}

# Compare two GenotypeHandles for the LD-sketch equality contract. Thin
# wrapper over the shared `.requireMatchingLdSketches` helper (R/ld.R)
# using the "lenient" null policy: a NULL on either side skips the check
# (only colocboostPipeline allows that, since some bundles only have a
# QTL side or only a GWAS side).
.cbRequireMatchingLdSketches <- function(qtlLd, gwasLd) {
    .requireMatchingLdSketches(
        qtlLd,
        gwasLd,
        pipelineName = "colocboostPipeline",
        nullPolicy = "lenient"
    )
}

# Combine sumstat bundles into the (sumstat-list, LD-list, dict) shape
# colocboost expects. Deduplicates identical LD matrices so dict_sumstatLD
# can point multiple sumstats at one LD.
.cbMergeSumstatBundles <- function(bundles) {
    if (length(bundles) == 0L) {
        return(list(
            sumstat = list(),
            LD = list(),
            dict_sumstatLD = matrix(integer(0), ncol = 2L)
        ))
    }
    ldUnique <- list()
    ldMatch <- integer(length(bundles))
    for (i in seq_along(bundles)) {
        ld <- bundles[[i]]$LD
        matched <- which(map_lgl(ldUnique, identical, ld))
        if (length(matched) > 0L) {
            ldMatch[[i]] <- matched[[1L]]
        } else {
            ldUnique[[length(ldUnique) + 1L]] <- ld
            ldMatch[[i]] <- length(ldUnique)
        }
    }
    names(ldUnique) <- str_c("LD", seq_along(ldUnique))
    sumstat <- map(bundles, "sumstat")
    names(sumstat) <- names(bundles)
    dict <- cbind(seq_along(bundles), ldMatch)
    colnames(dict) <- c("sumstat", "LD")
    list(sumstat = sumstat, LD = ldUnique, dict_sumstatLD = dict)
}

# Build an empty result skeleton consistent with what the per-method
# dispatch fills in.
.cbEmptyResult <- function() {
    list(
        xqtl_coloc = NULL,
        joint_gwas = NULL,
        separate_gwas = NULL,
        computing_time = list(
            Analysis = list(
                xqtl_coloc = NULL,
                joint_gwas = NULL,
                separate_gwas = NULL
            )
        )
    )
}

# Shared dispatch: accepts a fully-prepared individual bundle (possibly
# NULL) plus a sumstat bundle (possibly empty) and runs the three
# colocboost variants the user requested.
.cbRunVariants <- function(
    individualBundle,
    sumstatBundle,
    xqtlColoc,
    jointGwas,
    separateGwas,
    focalTrait,
    dotArgs,
    qtlLdSketch = NULL
) {
    results <- .cbEmptyResult()
    hasInd <- !is.null(individualBundle)
    hasSs <- length(sumstatBundle$sumstat) > 0L
    if (!hasInd && !hasSs) {
        msg <- glue(
            "colocboostPipeline: no QTL inputs remain after selection. ",
            "Nothing to run."
        )
        inform(msg)
        return(.cbEmptyResultObject())
    }
    if (isTRUE(xqtlColoc) && hasInd) {
        run <- .cbRunXqtlOnly(individualBundle, focalTrait, dotArgs)
        results$xqtl_coloc <- run$result
        results$computing_time$Analysis$xqtl_coloc <- run$time
    }
    if (isTRUE(jointGwas) && hasSs) {
        run <- .cbRunJointGwas(individualBundle, sumstatBundle, hasInd, dotArgs)
        results$joint_gwas <- run$result
        results$computing_time$Analysis$joint_gwas <- run$time
    }
    if (isTRUE(separateGwas) && hasSs) {
        run <- .cbRunSeparateGwas(
            individualBundle,
            sumstatBundle,
            hasInd,
            dotArgs
        )
        results$separate_gwas <- run$result
        results$computing_time$Analysis$separate_gwas <- run$time
    }
    .cbToResultObject(
        results,
        .cbOutcomeInfo(individualBundle, sumstatBundle, hasInd),
        qtlLdSketch
    )
}

# An empty ColocBoostResult: no sets, but the full column schema, so a caller
# reading `result$cosNpc` gets an empty column rather than NULL exactly when
# there is nothing to report.
# @noRd
.cbEmptyResultObject <- function() {
    ColocBoostResult(
        results = list(),
        analysis = character(0),
        gwasStudy = character(0),
        outcomeInfo = .cbEmptyOutcomeInfo()
    )
}

# Flatten the three raw colocboost runs into one ColocBoostResult.
#
# The three analyses are NOT parallel structures: xqtl_coloc and joint_gwas are
# each a single colocboost object, while separate_gwas is a NAMED LIST of them,
# one per GWAS study -- one level deeper. That is why `gwasStudy` exists as a
# key: without it the separate-GWAS sets would be indistinguishable from each
# other once flattened.
# @noRd
.cbToResultObject <- function(results, outcomeInfo, ldSketch = NULL) {
    runs <- list()
    analysis <- character(0)
    gwasStudy <- character(0)
    for (nm in c("xqtl_coloc", "joint_gwas")) {
        if (!is.null(results[[nm]])) {
            runs[[length(runs) + 1L]] <- results[[nm]]
            analysis <- c(analysis, nm)
            gwasStudy <- c(gwasStudy, NA_character_)
        }
    }
    sep <- results$separate_gwas
    if (!is.null(sep) && length(sep) > 0L) {
        keys <- names(sep) %||% as.character(seq_along(sep))
        for (i in seq_along(sep)) {
            if (is.null(sep[[i]])) {
                next
            }
            runs[[length(runs) + 1L]] <- sep[[i]]
            analysis <- c(analysis, "separate_gwas")
            gwasStudy <- c(gwasStudy, keys[[i]])
        }
    }
    ColocBoostResult(
        results = runs,
        analysis = analysis,
        gwasStudy = gwasStudy,
        outcomeInfo = outcomeInfo,
        ldSketch = ldSketch,
        computingTime = results$computing_time %||% list()
    )
}

# xQTL-only ColocBoost run -> list(result, time).
# @noRd
.cbRunXqtlOnly <- function(individualBundle, focalTrait, dotArgs) {
    traits <- individualBundle$outcomeNames
    focalIdx <- if (!is.null(focalTrait) && is_in(focalTrait, traits)) {
        which(traits == focalTrait)
    } else {
        NULL
    }
    nCtx <- length(individualBundle$Y)
    msg <- glue(
        "====== Performing xQTL-only ColocBoost on {nCtx} contexts. ====="
    )
    inform(msg)
    args <- c(
        list(
            X = individualBundle$X,
            Y = individualBundle$Y,
            dict_YX = individualBundle$dict_YX,
            outcome_names = traits,
            focal_outcome_idx = focalIdx,
            output_level = 2
        ),
        dotArgs
    )
    run <- .cbRun("xQTL-only ColocBoost", args)
    list(result = run$result, time = run$time)
}

# The (outcome name -> study / context / trait / dataForm) lookup for one run.
#
# Outcome names reach colocboost as bare labels; this is what lets a
# ColocBoostResult answer "which contexts colocalized" instead of handing the
# caller a string to parse. Sumstat outcomes are keyed by study name and carry
# no context or trait of their own, which is recorded as NA rather than
# invented.
# @noRd
.cbOutcomeInfo <- function(individualBundle, sumstatBundle, hasInd) {
    parts <- list()
    if (isTRUE(hasInd) && !is.null(individualBundle$outcomeInfo)) {
        parts[[length(parts) + 1L]] <- individualBundle$outcomeInfo
    }
    ssNames <- names(sumstatBundle$sumstat)
    if (length(ssNames) > 0L) {
        parts[[length(parts) + 1L]] <- tibble(
            name = ssNames,
            context = NA_character_,
            trait = NA_character_,
            study = ssNames,
            dataForm = "sumstats"
        )
    }
    if (length(parts) == 0L) {
        return(.cbEmptyOutcomeInfo())
    }
    bind_rows(parts)
}

# @noRd
.cbEmptyOutcomeInfo <- function() {
    tibble(
        name = character(0),
        context = character(0),
        trait = character(0),
        study = character(0),
        dataForm = character(0)
    )
}

# Joint (non-focal) QTL + GWAS run -> list(result, time).
# @noRd
.cbRunJointGwas <- function(individualBundle, sumstatBundle, hasInd, dotArgs) {
    traits <- c(
        if (hasInd) individualBundle$outcomeNames else character(),
        names(sumstatBundle$sumstat)
    )
    ldArgs <- .cbBuildLdArgs(sumstatBundle$LD)
    nContexts <- if (hasInd) length(individualBundle$Y) else 0L
    nGwas <- length(sumstatBundle$sumstat)
    msg <- glue(
        "====== Performing non-focal GWAS-xQTL ColocBoost on ",
        "{nContexts} contexts and {nGwas} GWAS. ====="
    )
    inform(msg)
    args <- c(
        list(
            X = if (hasInd) individualBundle$X else NULL,
            Y = if (hasInd) individualBundle$Y else NULL,
            sumstat = sumstatBundle$sumstat,
            dict_YX = if (hasInd) individualBundle$dict_YX else NULL,
            dict_sumstatLD = sumstatBundle$dict_sumstatLD,
            outcome_names = traits,
            focal_outcome_idx = NULL,
            output_level = 2
        ),
        ldArgs,
        dotArgs
    )
    run <- .cbRun("Joint GWAS ColocBoost", args)
    list(result = run$result, time = run$time)
}

# Separate (focal) per-GWAS runs -> list(result, time) with per-study timing.
# @noRd
.cbRunSeparateGwas <- function(
    individualBundle,
    sumstatBundle,
    hasInd,
    dotArgs
) {
    ssNames <- names(sumstatBundle$sumstat)
    t1 <- Sys.time()
    separate <- set_names(
        map(
            seq_along(ssNames),
            .cbSeparateGwasAt,
            ssNames = ssNames,
            individualBundle = individualBundle,
            sumstatBundle = sumstatBundle,
            hasInd = hasInd,
            dotArgs = dotArgs
        ),
        ssNames
    )
    t2 <- Sys.time()
    list(
        result = separate,
        time = list(
            total = t2 - t1,
            n_studies = length(ssNames),
            average = if (length(ssNames) > 0L) {
                (t2 - t1) / length(ssNames)
            } else {
                NA
            }
        )
    )
}

# One focal GWAS-xQTL ColocBoost run -> its result.
# @noRd
.cbRunOneSeparateGwas <- function(
    i,
    study,
    individualBundle,
    sumstatBundle,
    hasInd,
    dotArgs
) {
    ldIdx <- sumstatBundle$dict_sumstatLD[i, 2L]
    ldArgs <- .cbBuildLdArgs(sumstatBundle$LD[ldIdx])
    traits <- c(
        if (hasInd) individualBundle$outcomeNames else character(),
        study
    )
    nContexts <- if (hasInd) length(individualBundle$Y) else 0L
    msg <- glue(
        "====== Performing focal GWAS-xQTL ColocBoost on {nContexts} ",
        "contexts and {study} GWAS. ====="
    )
    inform(msg)
    args <- c(
        list(
            X = if (hasInd) individualBundle$X else NULL,
            Y = if (hasInd) individualBundle$Y else NULL,
            sumstat = sumstatBundle$sumstat[i],
            dict_YX = if (hasInd) individualBundle$dict_YX else NULL,
            outcome_names = traits,
            focal_outcome_idx = length(traits),
            output_level = 2
        ),
        ldArgs,
        dotArgs
    )
    .cbRun(str_c("Separate GWAS ColocBoost for ", study), args)$result
}

# =============================================================================
# Allele harmonization (the alleleFlip = TRUE path)
# =============================================================================

# Relabel a residualized genotype matrix's columns to the canonical coding,
# negating columns whose allele order is swapped relative to canonical. For a
# centered / residualized dosage, negation IS the allele flip (flipping the
# counted allele negates the centered genotype). Columns with no canonical
# match (e.g. a multi-allelic secondary allele) are dropped.
# @noRd
.cbFlipMatrixToCanonical <- function(m, canonical) {
    mm <- matchVariants(colnames(m), canonical)
    if (length(mm$idxA) == 0L) {
        return(m[, integer(0), drop = FALSE])
    }
    o <- order(mm$idxA)
    ia <- mm$idxA[o]
    ib <- mm$idxB[o]
    sgn <- mm$sign[o]
    out <- m[, ia, drop = FALSE]
    flip <- which(sgn == -1)
    if (length(flip) > 0L) {
        out[, flip] <- -out[, flip]
    }
    colnames(out) <- canonical[ib]
    out
}

# Relabel a (sumstat, LD) pair to the canonical coding: flip the z-score and the
# LD sign-submatrix for variants swapped relative to canonical
# (LD_ij -> sign_i * sign_j * LD_ij), drop unmatched variants, and relabel to
# canonical. The sumstat and its LD share one sign vector, so they stay
# consistent. Returns NULL when no variant matches.
# @noRd
.cbFlipPairToCanonical <- function(p, canonical) {
    mm <- matchVariants(p$sumstat$variant, canonical)
    if (length(mm$idxA) == 0L) {
        return(NULL)
    }
    o <- order(mm$idxA)
    ia <- mm$idxA[o]
    ib <- mm$idxB[o]
    sgn <- mm$sign[o]
    ss <- p$sumstat[ia, , drop = FALSE]
    ss$z <- ss$z * sgn
    ss$variant <- canonical[ib]
    ld <- p$LD[ia, ia, drop = FALSE] * outer(sgn, sgn)
    dimnames(ld) <- list(canonical[ib], canonical[ib])
    list(sumstat = ss, LD = ld, variantIds = canonical[ib])
}

# Harmonize the allele coding of every colocboost input (individual X columns,
# each sumstat's z / variant, and its LD) to a single per-locus canonical
# coding, so a variant that appears with opposite ref/alt across sources is
# combined with a consistent sign rather than dropped or silently mis-signed.
# The canonical id for a (chrom, pos) locus is the first-seen id across all
# sources (its ref/alt order becomes the shared coding); rsIDs are their own
# locus. Only invoked when alleleFlip = TRUE.
# @noRd
.cbHarmonizeAlleles <- function(individualBundle, pairs) {
    ids <- character(0)
    if (!is.null(individualBundle)) {
        ids <- c(
            ids,
            unlist(map(individualBundle$X, colnames), use.names = FALSE)
        )
    }
    ids <- c(
        ids,
        unlist(map(pairs, list("sumstat", "variant")), use.names = FALSE)
    )
    ids <- unique(ids[!is.na(ids)])
    if (length(ids) == 0L) {
        return(list(individualBundle = individualBundle, pairs = pairs))
    }
    parsed <- parseVariantId(ids)
    ok <- !is.na(parsed$chrom) & !is.na(parsed$pos)
    locus <- if_else(ok, str_c(parsed$chrom, parsed$pos, sep = ":"), ids)
    canonical <- ids[!duplicated(locus)]
    if (!is.null(individualBundle)) {
        individualBundle$X <- map(
            individualBundle$X,
            .cbFlipMatrixToCanonical,
            canonical = canonical
        )
    }
    pairs <- map(pairs, .cbFlipPairToCanonical, canonical = canonical)
    pairs <- pairs[!map_lgl(pairs, is.null)]
    list(individualBundle = individualBundle, pairs = pairs)
}

# Top-level driver shared by all input methods. qtlPairs and gwasPairs
# are per-tuple lists of `list(sumstat, LD)` produced by the per-class
# bundle helpers; they are merged here so dict_sumstatLD can dedupe
# identical LD matrices across QTL and GWAS sides.
.cbDriver <- function(
    individualBundle,
    qtlPairs,
    gwasSumStats,
    xqtlColoc,
    jointGwas,
    separateGwas,
    focalTrait,
    dotArgs,
    qtlLdSketch = NULL,
    alleleFlip = TRUE,
    cutoffs = NULL
) {
    if (!isTRUE(xqtlColoc) && !isTRUE(jointGwas) && !isTRUE(separateGwas)) {
        inform("colocboostPipeline: no analysis flag is TRUE; nothing to do.")
        return(.cbEmptyResultObject())
    }
    combinedPairs <- .cbAppendGwasPairs(
        qtlPairs,
        gwasSumStats,
        qtlLdSketch,
        cutoffs = cutoffs
    )
    # Harmonize allele coding across all sources to a shared per-locus canonical
    # so swapped variants are combined with a consistent sign (alleleFlip =
    # TRUE); alleleFlip = FALSE leaves the names-only canonicalization done at
    # the source builders, which keeps swapped variants distinct.
    if (isTRUE(alleleFlip)) {
        harmonized <- .cbHarmonizeAlleles(individualBundle, combinedPairs)
        individualBundle <- harmonized$individualBundle
        combinedPairs <- harmonized$pairs
    }
    sumstatBundle <- .cbMergeSumstatBundles(combinedPairs)
    .cbRunVariants(
        individualBundle,
        sumstatBundle,
        xqtlColoc,
        jointGwas,
        separateGwas,
        focalTrait,
        dotArgs,
        qtlLdSketch = qtlLdSketch
    )
}

# Append the GWAS sumstat pairs to the QTL pairs, QC-gating gwasSumStats and
# make.unique-ing colliding labels. (Sequential key resolution -- kept as a
# loop.)
# @noRd
.cbAppendGwasPairs <- function(
    qtlPairs,
    gwasSumStats,
    qtlLdSketch,
    cutoffs = NULL
) {
    if (is.null(gwasSumStats)) {
        return(qtlPairs)
    }
    .cbRequireSumStatsQc(gwasSumStats, "gwasSumStats")
    if (!is.null(qtlLdSketch)) {
        .cbRequireMatchingLdSketches(qtlLdSketch, getLdSketch(gwasSumStats))
    }
    gwasPairs <- .cbGwasSumStatsBundle(gwasSumStats, cutoffs = cutoffs)
    combinedPairs <- qtlPairs
    for (label in names(gwasPairs)) {
        key <- label
        if (is_in(key, names(combinedPairs))) {
            key <- make.unique(
                c(names(combinedPairs), key)
            )[length(combinedPairs) + 1L]
        }
        combinedPairs[[key]] <- gwasPairs[[label]]
    }
    combinedPairs
}

# =============================================================================
# Methods
# =============================================================================

# The QtlDataset colocboost run: screen spec, individual-level bundle, shared
# driver. Split out so the method is only its signature -- twenty formals plus
# a body is past the length a reviewer (or BiocCheck) will accept, and the
# signature is the part that cannot be shortened.
# `p` is the method's named arguments; `dots` is its `...`, which
# as.list(environment()) does not capture.
# @noRd
.cbQtlDatasetDrive <- function(p, dots) {
    screenSpec <- .cbScreenSpec(
        p$pipCutoffToSkip,
        p$absZCutoffToSkip,
        p$bfCutoffToSkip,
        p$logBfCutoffToSkip
    )
    indBundle <- .cbIndividualBundle(
        p$qtlData,
        contexts = p$contexts,
        traitId = p$traitId,
        region = p$region,
        cisWindow = p$cisWindow,
        samples = p$samples,
        pipCutoffToSkip = screenSpec
    )
    .cbDriver(
        indBundle,
        qtlPairs = list(),
        p$gwasSumStats,
        p$xqtlColoc,
        p$jointGwas,
        p$separateGwas,
        p$focalTrait,
        dots,
        alleleFlip = p$alleleFlip,
        # The QTL side is individual-level here, but the GWAS side may still
        # be sumstats, so the panel cutoffs still apply to it.
        cutoffs = .panelCutoffs(list(
            mafCutoff = p$mafCutoff,
            macCutoff = p$macCutoff,
            imissCutoff = p$imissCutoff
        ))
    )
}


#' @rdname colocboostPipeline
#' @export
setMethod(
    "colocboostPipeline",
    "QtlDataset",
    function(
        qtlData,
        gwasSumStats = NULL,
        contexts = NULL,
        traitId = NULL,
        region = NULL,
        cisWindow = NULL,
        focalTrait = NULL,
        xqtlColoc = TRUE,
        jointGwas = FALSE,
        separateGwas = FALSE,
        samples = NULL,
        mafCutoff = 0,
        macCutoff = 0,
        imissCutoff = 1,
        pipCutoffToSkip = 0,
        absZCutoffToSkip = 0,
        bfCutoffToSkip = 0,
        logBfCutoffToSkip = 0,
        alleleFlip = TRUE,
        ...
    ) {
        .cbQtlDatasetDrive(as.list(environment()), list(...))
    }
)

#' @rdname colocboostPipeline
#' @export
setMethod(
    "colocboostPipeline",
    "QtlSumStats",
    function(
        qtlData,
        gwasSumStats = NULL,
        contexts = NULL,
        traitId = NULL,
        region = NULL,
        cisWindow = NULL,
        focalTrait = NULL,
        xqtlColoc = TRUE,
        jointGwas = FALSE,
        separateGwas = FALSE,
        alleleFlip = TRUE,
        mafCutoff = 0,
        macCutoff = 0,
        imissCutoff = 1,
        ...
    ) {
        .cbRequireSumStatsQc(qtlData, "qtlData")
        dotArgs <- list(...)
        cutoffs <- .panelCutoffs(list(
            mafCutoff = mafCutoff,
            macCutoff = macCutoff,
            imissCutoff = imissCutoff
        ))
        qtlPairs <- .cbQtlSumStatsBundle(
            qtlData,
            contexts = contexts,
            traitId = traitId,
            cutoffs = cutoffs
        )
        .cbDriver(
            individualBundle = NULL,
            qtlPairs = qtlPairs,
            gwasSumStats = gwasSumStats,
            xqtlColoc = xqtlColoc,
            jointGwas = jointGwas,
            separateGwas = separateGwas,
            focalTrait = focalTrait,
            dotArgs = dotArgs,
            qtlLdSketch = getLdSketch(qtlData),
            alleleFlip = alleleFlip,
            cutoffs = cutoffs
        )
    }
)

#' @rdname colocboostPipeline
#' @export
setMethod(
    "colocboostPipeline",
    "MultiStudyQtlDataset",
    function(
        qtlData,
        gwasSumStats = NULL,
        contexts = NULL,
        traitId = NULL,
        region = NULL,
        cisWindow = NULL,
        focalTrait = NULL,
        xqtlColoc = TRUE,
        jointGwas = FALSE,
        separateGwas = FALSE,
        samples = NULL,
        mafCutoff = 0,
        macCutoff = 0,
        imissCutoff = 1,
        pipCutoffToSkip = 0,
        absZCutoffToSkip = 0,
        bfCutoffToSkip = 0,
        logBfCutoffToSkip = 0,
        alleleFlip = TRUE,
        ...
    ) {
        p <- as.list(environment())
        p$dotArgs <- list(...)
        .cbPipelineMultiStudy(p)
    }
)

# MultiStudyQtlDataset colocboost worker: aggregate the per-study individual
# bundles + embedded sumstats, then dispatch to the shared driver.
# @noRd
.cbPipelineMultiStudy <- function(p) {
    screenSpec <- .cbScreenSpec(
        p$pipCutoffToSkip,
        p$absZCutoffToSkip,
        p$bfCutoffToSkip,
        p$logBfCutoffToSkip
    )
    indBundle <- .cbMultiStudyIndBundle(p, screenSpec)
    cutoffs <- .panelCutoffs(p)
    ss <- .cbMultiStudySumstats(
        p$qtlData,
        p$contexts,
        p$traitId,
        cutoffs = cutoffs
    )
    .cbDriver(
        indBundle,
        ss$qtlPairs,
        p$gwasSumStats,
        p$xqtlColoc,
        p$jointGwas,
        p$separateGwas,
        p$focalTrait,
        p$dotArgs,
        qtlLdSketch = ss$qtlLdSketch,
        alleleFlip = p$alleleFlip,
        cutoffs = cutoffs
    )
}

# Aggregate the individual-level bundles across all QtlDataset members. Per-
# study trait names are prefixed with "{study}:" so colocboost sees distinct
# outcomes when two studies share a trait. Returns the combined bundle or NULL.
# (Sequential offset-shifted merge -- kept as a loop.)
# @noRd
.cbMultiStudyIndBundle <- function(p, screenSpec) {
    qtlDatasets <- getQtlDatasets(p$qtlData)
    combinedX <- list()
    combinedY <- list()
    combinedDict <- matrix(integer(0), ncol = 2L)
    colnames(combinedDict) <- c("Y", "X")
    combinedOutcomes <- character()
    combinedInfo <- list()
    for (study in names(qtlDatasets)) {
        sub <- .cbIndividualBundle(
            qtlDatasets[[study]],
            contexts = p$contexts,
            traitId = p$traitId,
            region = p$region,
            cisWindow = p$cisWindow,
            samples = p$samples,
            pipCutoffToSkip = screenSpec
        )
        if (is.null(sub)) {
            next
        }
        sub <- .cbPrefixStudyNames(sub, study)
        xOffset <- length(combinedX)
        yOffset <- length(combinedY)
        combinedX <- c(combinedX, sub$X)
        combinedY <- c(combinedY, sub$Y)
        shifted <- sub$dict_YX
        shifted[, "Y"] <- shifted[, "Y"] + yOffset
        shifted[, "X"] <- shifted[, "X"] + xOffset
        combinedDict <- rbind(combinedDict, shifted)
        combinedOutcomes <- c(combinedOutcomes, sub$outcomeNames)
        combinedInfo[[length(combinedInfo) + 1L]] <- sub$outcomeInfo
    }
    if (length(combinedX) == 0L) {
        return(NULL)
    }
    list(
        X = combinedX,
        Y = combinedY,
        dict_YX = combinedDict,
        outcomeNames = combinedOutcomes,
        outcomeInfo = bind_rows(combinedInfo)
    )
}

# Prefix a study's X names + outcome names with "{study}:" (Y names track the
# outcome names).
# @noRd
.cbPrefixStudyNames <- function(sub, study) {
    names(sub$X) <- str_c(study, names(sub$X), sep = ":")
    sub$outcomeNames <- str_c(study, sub$outcomeNames, sep = ":")
    names(sub$Y) <- sub$outcomeNames
    # The lookup table is keyed on the outcome NAME, so it has to be renamed in
    # the same pass -- a stale key silently drops every outcome of this study
    # from the identity join.
    if (!is.null(sub$outcomeInfo) && nrow(sub$outcomeInfo) > 0L) {
        sub$outcomeInfo$name <- str_c(study, sub$outcomeInfo$name, sep = ":")
        sub$outcomeInfo$study <- study
    }
    sub
}

# Sumstat side of a MultiStudyQtlDataset: bundle any embedded QtlSumStats.
# Returns list(qtlPairs, qtlLdSketch).
# @noRd
.cbMultiStudySumstats <- function(
    qtlData,
    contexts,
    traitId,
    cutoffs = NULL
) {
    embeddedSs <- getSumStats(qtlData)
    if (is.null(embeddedSs)) {
        return(list(qtlPairs = list(), qtlLdSketch = NULL))
    }
    .cbRequireSumStatsQc(embeddedSs, "MultiStudyQtlDataset@sumStats")
    list(
        qtlPairs = .cbQtlSumStatsBundle(
            embeddedSs,
            contexts = contexts,
            traitId = traitId,
            cutoffs = cutoffs
        ),
        qtlLdSketch = getLdSketch(embeddedSs)
    )
}

#' @rdname colocboostPipeline
#' @export
setMethod(
    "colocboostPipeline",
    "ANY",
    function(qtlData, gwasSumStats = NULL, ...) {
        cls <- class(qtlData)[[1L]]
        msg <- glue(
            "colocboostPipeline does not accept inputs of class ",
            "'{cls}'. Pass a QtlDataset, QtlSumStats, or ",
            "MultiStudyQtlDataset for QTL data."
        )
        abort(msg)
    }
)

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# TRUE when a matrix is non-square (a samples x variants genotype reference).
# @noRd
.cbIsNonSquare <- function(m) {
    nrow(m) != ncol(m)
}

# Single-effect screen decision for outcome column `j` (drop on failure).
# @noRd
.cbScreenOutcome <- function(j, X, Y, spec) {
    .fmSerScreen(X, Y[, j], spec, fallback = FALSE)
}

# Run one separate-GWAS colocboost fit for sumstat index `i`.
# @noRd
.cbSeparateGwasAt <- function(
    i,
    ssNames,
    individualBundle,
    sumstatBundle,
    hasInd,
    dotArgs
) {
    .cbRunOneSeparateGwas(
        i,
        ssNames[[i]],
        individualBundle,
        sumstatBundle,
        hasInd,
        dotArgs
    )
}
