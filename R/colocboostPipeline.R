#' @include qtlSumStats.R
#' @title ColocBoost multi-trait colocalization pipeline (S4)
#' @description Protocol-level multi-trait colocalization analysis using
#'   \pkg{colocboost}. Dispatches on the QTL input type:
#'   \itemize{
#'     \item \code{QtlDataset} -- single-study, individual-level
#'           multi-context data. Per-context residualized X / Y are
#'           extracted from the dataset (filtering knobs on the
#'           constructor apply lazily inside the accessors).
#'     \item \code{QtlSumStats} -- summary-statistic-only QTL data with an
#'           LD reference per row (\code{ldReference}). Must already have
#'           been passed through \code{\link{summaryStatsQc}} (the
#'           pipeline rejects inputs whose \code{qcOptions()} is empty).
#'           Every analysis variant is available on this input, including
#'           \code{xqtlColoc}, which colocalizes the QTL studies against
#'           each other with no GWAS involved.
#'     \item \code{MultiStudyQtlDataset} -- a mixture of one or more
#'           individual-level \code{QtlDataset} studies and an optional
#'           \code{QtlSumStats} collection.
#'   }
#'   GWAS is optional and always passed separately as a
#'   \code{GwasSumStats} object (must also be QC'd).
#'
#' \code{colocboostPipeline} does \strong{not} accept a fine-mapping result
#' for either side; colocboost has its own variable-selection algorithm.
#'
#' @section QC contract:
#'   \itemize{
#'     \item Individual-level QC (MAF / MAC / X-variance / per-sample
#'           missingness, sample / variant restrictions) is
#'           \code{genotypeFilterParam}, applied to each block of genotypes
#'           the pipeline extracts. A \code{QtlDataset} records no filter
#'           of its own.
#'     \item \code{panelFilterParam} acts at ANALYSIS time on the sumstat
#'           sides, measured against the \strong{LD reference panel}: a
#'           variant whose panel genotypes fall below its cutoffs is dropped
#'           before the LD matrix is built. It applies to a
#'           \code{QtlSumStats} QTL side and to a \code{GwasSumStats} GWAS
#'           side (so it still bites when the QTL side is
#'           individual-level). Its defaults filter nothing. Unlike
#'           \code{summaryStatsQc(raissParam = ...)}, which only bounds what
#'           RAISS will impute, these discard \emph{observed} variants.
#'     \item All summary-statistic QC (variant filters, harmonization
#'           against the \code{ldReference}, LD-mismatch detection, RAISS
#'           imputation, etc.) lives in
#'           \code{\link{summaryStatsQc}}. The pipeline rejects any
#'           \code{QtlSumStats} or \code{GwasSumStats} where
#'           \code{length(qcOptions(x)) == 0L}.
#'   }
#'
#' @section Analysis variants:
#'   \itemize{
#'     \item \code{xqtlColoc} (default \code{TRUE}): run a colocboost
#'           model over the QTL outcomes only, excluding
#'           \code{gwasSumStats}. Works for either data form -- the
#'           individual-level contexts, the summary-level QTL studies, or
#'           both together when a \code{MultiStudyQtlDataset} carries a
#'           mixture. Focal only on a QTL outcome named by
#'           \code{focalOutcome}.
#'     \item \code{jointGwas} (default \code{FALSE}): run one colocboost
#'           model that combines all QTL contexts/studies with the
#'           supplied \code{gwasSumStats} studies. By default it is focal on
#'           the GWAS when \code{gwasSumStats} holds a single study, and
#'           non-focal when it holds several; \code{focalOutcome} chooses.
#'     \item \code{separateGwas} (default \code{FALSE}): run one focal
#'           colocboost model per GWAS study, where the GWAS is the
#'           focal outcome.
#'   }
#'
#' @param qtlData One of \code{QtlDataset}, \code{QtlSumStats}, or
#'   \code{MultiStudyQtlDataset}.
#' @param gwasSumStats Optional \code{GwasSumStats} with the GWAS studies to
#'   colocalize against. \code{NULL} to skip GWAS colocalization. Each study
#'   is one outcome: the union of its rows that the QTL data overlap, so a
#'   study split by chromosome or LD block contributes just the one or two
#'   rows the QTL region reaches, and a study none of whose rows it reaches is
#'   left out. QTL variants outside every loaded block of a study the QTL data
#'   overlap are an error where the blocks' ranges were given
#'   (\code{regionGiven}, see \code{\linkS4class{GwasSumStats}}), since they
#'   lie in a block that was not loaded, and a warning for more than one in a
#'   row past a range computed from the GWAS variants. Those rows must share
#'   one LD reference, but it need not be the
#'   QTL side's: every outcome's summary statistics are matched against their
#'   own reference. When some references are genotype panels and others are
#'   \code{LdMatrix} objects, the runs that combine them pass colocboost the
#'   panels' LD as matrices, since colocboost takes LD matrices or reference
#'   genotypes but not both.
#' @param contexts Optional character vector of context names to restrict
#'   the individual-level / QtlSumStats QTL analysis to. When \code{NULL}
#'   (default), every context present is used. A MultiStudyQtlDataset applies
#'   it to every study: a context no study has is an error, and one only some
#'   studies lack is skipped for those, with a warning.
#' @param traitId Optional character vector of trait identifiers to restrict
#'   the analysis to. When supplied with an individual-level
#'   \code{QtlDataset} input, \code{cisWindow} is required: each trait's
#'   position \eqn{\pm} \code{cisWindow} is the variant window. In a
#'   MultiStudyQtlDataset a study carrying none of them sits out; the run is
#'   refused only when no study or summary-statistics row carries one.
#' @param region Optional single-range \code{GRanges} describing the analysis
#'   window. Mutually exclusive with \code{traitId}. Individual-level input
#'   only: a summary-statistics tuple
#'   keeps the window it was built with (its \code{region}), so the
#'   \code{QtlSumStats} method does not take one --- matching
#'   \code{\link{fineMappingPipeline}} and \code{\link{twasWeightsPipeline}}.
#' @param cisWindow Optional cis window in basepairs. Individual-level
#'   input: required with \code{traitId}, optional with \code{region}.
#'   QtlSumStats input: each selected tuple is analysed over its
#'   \code{traitPos} \eqn{\pm} \code{cisWindow}, within its \code{region};
#'   every selected tuple needs a \code{traitPos}. A MultiStudyQtlDataset
#'   applies it to both kinds of study it holds.
#' @param samples Optional character vector of sample IDs to restrict the
#'   analysis to; \code{NULL} (default) uses all samples.
#' @param focalOutcome The outcome a run treats as focal (colocboost's
#'   \code{focal_outcome_idx}). One of:
#'   \itemize{
#'     \item \code{NULL} (default): the GWAS study, in the \code{jointGwas}
#'           run, when \code{gwasSumStats} holds exactly one. With several
#'           there is no one study to choose, so no run is focal.
#'     \item A GWAS study name, as a plain string: that study is focal in
#'           the \code{jointGwas} run. It must be a study of
#'           \code{gwasSumStats}.
#'     \item A named character vector with any of \code{study},
#'           \code{context} and \code{trait} that identifies exactly one QTL
#'           outcome, e.g. \code{c(context = "brain", trait = "ENSG1")}: that
#'           outcome is focal in every run holding it, the \code{xqtlColoc}
#'           run and the \code{jointGwas} run.
#'     \item \code{NA}: no run is focal.
#'   }
#'   The \code{separateGwas} runs are always focal on their own GWAS study.
#' @param xqtlColoc,jointGwas,separateGwas Logical flags selecting which
#'   colocboost variants to run.
#' @param signalScreenParam Individual-level pre-filter (ports the legacy
#'   \code{pip_cutoff_to_skip_ind}), built with
#'   \code{\link{SignalScreenParam}}. For each context every outcome is fit
#'   with a single-effect SuSiE (\code{L = 1}) and dropped unless the chosen
#'   metric exceeds its cutoff; a context with no surviving outcome is
#'   skipped:
#'   \itemize{
#'     \item \code{pip} --- the maximum PIP. A negative cutoff uses
#'       \code{3 / n_variants}. Uniquely here it may also be a
#'       \strong{context-named} numeric vector rather than a scalar,
#'       applying a different cutoff per context.
#'     \item \code{absZ} --- the maximum marginal \code{|z|}.
#'     \item \code{bf} / \code{logBf} --- the maximum per-variant Bayes
#'       factor or log Bayes factor from the \code{L = 1} fit.
#'   }
#'   Only one metric may be enabled, which
#'   \code{\link{SignalScreenParam}} enforces when you build it. Unset (the
#'   default) screens nothing. Summary-statistic skipping is handled upstream
#'   by \code{\link{summaryStatsQc}}.
#' @param panelFilterParam Analysis-time filters applied to the
#'   summary-statistic sides against the LD reference panel, built with
#'   \code{\link{PanelFilterParam}}: a variant whose panel genotypes fall
#'   below its cutoffs is dropped before the LD matrix is built. The defaults
#'   filter nothing. Individual-level genotypes are filtered by
#'   \code{genotypeFilterParam} instead.
#' @param genotypeFilterParam Which individual-level variants and samples to
#'   keep, built with \code{\link{GenotypeFilterParam}}: its
#'   \code{keepSamples} narrows the dataset once, and its MAF, MAC,
#'   variance, missingness, \code{keepVariants} and \code{keepIndel}
#'   settings apply to each extracted block, as in
#'   \code{\link{fineMappingPipeline}}. The default filters nothing.
#'   \code{QtlDataset} and \code{MultiStudyQtlDataset} input only.
#' @param phenotypeFilterParam The run's phenotype QC, built with
#'   \code{\link{PhenotypeFilterParam}}. colocboost takes missing phenotype
#'   values as they are, so unset keeps them; \code{"drop"} or
#'   \code{"impute"} applies that instead, and \code{outlierAction =
#'   "drop"} removes outlier samples. \code{QtlDataset} and
#'   \code{MultiStudyQtlDataset} input only.
#' @param residualizationParam Covariate residualization settings for
#'   individual-level input, built with \code{\link{ResidualizationParam}}.
#'   Reaches both \code{\link{residualizedGenotypes}} and
#'   \code{\link{residualizedPhenotypes}}, so it governs what is regressed
#'   out of each side and --- through \code{scaleResiduals} --- whether the
#'   columns are scaled to unit variance, as in
#'   \code{\link{fineMappingPipeline}} and
#'   \code{\link{twasWeightsPipeline}}. \strong{Ignored} on
#'   \code{QtlSumStats} input, which has nothing to residualize.
#' @param alleleFlip Logical, default \code{TRUE}. When TRUE, harmonize variants
#'   across the individual X, sumstats, and LD by (chrom, pos) with ref/alt
#'   swaps recognized (flipping z / residualized dosage / LD to a shared
#'   coding); when FALSE, match on exact alleles only (names-only), so a ref/alt
#'   swap is treated as a distinct variant.
#' @param methodArgs Additional arguments forwarded to
#'   \code{colocboost::colocboost}, built with
#'   \code{\link{ColocboostOptions}}. Names are checked against that
#'   function's formals, so a misspelled option errors instead of being
#'   dropped inside the run. The data slots the pipeline assembles itself
#'   (\code{X}, \code{Y}, \code{sumstat}, \code{LD}, \code{X_ref},
#'   \code{focal_outcome_idx}) are refused.
#' @param ... Required by the generic; the methods take no further arguments.
#' @return A \code{\linkS4class{ColocBoostResult}}: one element per
#'   confidence set (CoS) across every analysis that ran, holding that set's
#'   member variants with their \code{vcp}. The \code{analysis} column marks
#'   which run each set came from (\code{xqtl_coloc}, \code{joint_gwas} or
#'   \code{separate_gwas}), and \code{gwasStudy} distinguishes the per-study
#'   \code{separate_gwas} runs. Outcome-specific (uncolocalized) sets are
#'   included with \code{isColocalized = FALSE}.
#'
#'   Project it with \code{\link{colocPairs}}, \code{\link{colocVariants}}
#'   or \code{\link{colocBoostOutcomes}}; timings are on
#'   \code{\link{computingTime}} and the region-wide marginal
#'   probabilities on \code{\link{regionVcp}}. \code{\link{ldReference}}
#'   returns the LD references the summary statistics were matched against,
#'   each stored once, and the \code{ldReference} column of
#'   \code{colocBoostOutcomes()} names the one each outcome used. Its
#'   \code{seqinfo()} names the genome build the inputs name; inputs naming
#'   different builds are refused.
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
    .cbAssertInputClass(qtlData)
    standardGeneric("colocboostPipeline")
})

# At least one analysis has to be selected. The three flags are independent
# -- any combination is a valid request -- so they stay separate arguments;
# all-FALSE is the one combination that is not a request at all, and it used
# to run silently and return nothing.
# @noRd
.cbAssertAnyRun <- function(xqtlColoc, jointGwas, separateGwas) {
    if (any(c(isTRUE(xqtlColoc), isTRUE(jointGwas), isTRUE(separateGwas)))) {
        return(invisible(NULL))
    }
    abort(glue(
        "colocboostPipeline: no analysis selected. Set at least one of ",
        "`xqtlColoc` (xQTL-only), `jointGwas` (all GWAS jointly) or ",
        "`separateGwas` (one run per GWAS study) to TRUE."
    ))
}

# @noRd
.cbAssertInputClass <- function(qtlData) {
    ok <- c("QtlDataset", "QtlSumStats", "MultiStudyQtlDataset")
    if (any(map_lgl(ok, is, object = qtlData))) {
        return(invisible(NULL))
    }
    abort(glue(
        "colocboostPipeline does not accept inputs of class ",
        "'{class(qtlData)[[1L]]}'. Pass a QtlDataset, QtlSumStats, or ",
        "MultiStudyQtlDataset for QTL data."
    ))
}

# focalOutcome's shape: NULL, NA, one GWAS study name, or a named character
# vector of `study` / `context` / `trait` naming one QTL outcome. Whether it
# names a real outcome is settled once the outcomes are built
# (.cbResolveFocal()).
# @noRd
.cbAssertFocalOutcome <- function(focalOutcome) {
    if (is.null(focalOutcome) || .cbFocalIsNa(focalOutcome)) {
        return(invisible(NULL))
    }
    if (.cbFocalShapeOk(focalOutcome)) {
        return(invisible(NULL))
    }
    abort(c(
        str_c(
            "colocboostPipeline: `focalOutcome` must be NULL, NA, one GWAS ",
            "study name, or a named character vector of `study`, `context` ",
            "and/or `trait` naming one QTL outcome."
        ),
        i = str_c(
            "e.g. focalOutcome = \"MyGWAS\", or ",
            "focalOutcome = c(context = \"brain\", trait = \"ENSG1\")."
        )
    ))
}

# @noRd
.cbFocalIsNa <- function(x) {
    is.atomic(x) && length(x) == 1L && is.null(names(x)) && is.na(x)
}

# One unnamed string, or a named character vector whose names are distinct
# fields of a QTL outcome.
# @noRd
.cbFocalShapeOk <- function(x) {
    if (!is.character(x) || length(x) == 0L || anyNA(x)) {
        return(FALSE)
    }
    nm <- names(x)
    if (is.null(nm)) {
        return(length(x) == 1L)
    }
    !anyDuplicated(nm) && all(is_in(nm, c("study", "context", "trait")))
}

# The outcome each run is focal on (NULL for none), from `focalOutcome`.
# `qtlInfo` lists the QTL outcomes (.cbQtlOutcomeInfo()) and `gwasKeys` maps
# each GWAS study that reached the runs to its outcome name (.cbGwasKeys()).
# A GWAS study can be focal only where it is an outcome, the joint run; a QTL
# outcome is focal in both runs that hold it.
# @noRd
.cbResolveFocal <- function(focalOutcome, gwasSumStats, qtlInfo, gwasKeys) {
    if (.cbFocalIsNa(focalOutcome)) {
        return(list(xqtl = NULL, joint = NULL))
    }
    studies <- unique(as.character(gwasSumStats$study))
    if (is.null(focalOutcome)) {
        return(list(xqtl = NULL, joint = .cbDefaultFocal(studies, gwasKeys)))
    }
    if (is.null(names(focalOutcome))) {
        return(list(
            xqtl = NULL,
            joint = .cbFocalGwasKey(focalOutcome, studies, gwasKeys)
        ))
    }
    name <- .cbFocalQtlName(focalOutcome, qtlInfo)
    list(xqtl = name, joint = name)
}

# The default focal outcome: the GWAS, when gwasSumStats holds one study and
# it reached the runs. Several studies leave no one study to choose.
# @noRd
.cbDefaultFocal <- function(studies, gwasKeys) {
    if (length(studies) != 1L || !is_in(studies, names(gwasKeys))) {
        return(NULL)
    }
    unname(gwasKeys[[studies]])
}

# The outcome name of GWAS study `study`, which must be a study of
# gwasSumStats and must have reached the runs.
# @noRd
.cbFocalGwasKey <- function(study, studies, gwasKeys) {
    if (!is_in(study, studies)) {
        known <- if (length(studies) == 0L) {
            "no gwasSumStats was given"
        } else {
            str_c("its studies are ", str_flatten(studies, ", "))
        }
        abort(c(
            glue(
                "colocboostPipeline: focalOutcome '{study}' names no GWAS ",
                "study ({known})."
            ),
            i = str_c(
                "Name a QTL outcome by a named vector, e.g. ",
                "c(context = \"brain\", trait = \"ENSG1\")."
            )
        ))
    }
    if (!is_in(study, names(gwasKeys))) {
        abort(glue(
            "colocboostPipeline: focalOutcome names GWAS study '{study}', ",
            "but none of its rows overlap the QTL data, so it is not an ",
            "outcome of the joint run."
        ))
    }
    unname(gwasKeys[[study]])
}

# The one QTL outcome a focalOutcome vector identifies.
# @noRd
.cbFocalQtlName <- function(focal, qtlInfo) {
    hit <- reduce(
        names(focal),
        .cbFocalField,
        focal = focal,
        qtlInfo = qtlInfo,
        .init = rep(TRUE, nrow(qtlInfo))
    )
    if (sum(hit) == 1L) {
        return(qtlInfo$name[hit])
    }
    .cbAbortFocalQtl(focal, qtlInfo, hit)
}

# @noRd
.cbFocalField <- function(hit, field, focal, qtlInfo) {
    hit & is_in(qtlInfo[[field]], focal[[field]])
}

# A focalOutcome vector that names no QTL outcome, or several.
# @noRd
.cbAbortFocalQtl <- function(focal, qtlInfo, hit) {
    asked <- str_flatten(str_c(names(focal), " = '", focal, "'"), ", ")
    if (any(hit)) {
        abort(c(
            glue(
                "colocboostPipeline: focalOutcome ({asked}) matches ",
                "{sum(hit)} QTL outcomes: ",
                "{.cbDescribeOutcomes(filter(qtlInfo, hit))}."
            ),
            i = "Add `study`, `context` or `trait` until it names one."
        ))
    }
    abort(c(
        glue(
            "colocboostPipeline: focalOutcome ({asked}) matches no QTL ",
            "outcome."
        ),
        i = glue("The QTL outcomes are: {.cbDescribeOutcomes(qtlInfo)}.")
    ))
}

# "study / context / trait" for up to ten outcomes, "-" for a part one lacks.
# @noRd
.cbDescribeOutcomes <- function(info) {
    if (nrow(info) == 0L) {
        return("none")
    }
    shown <- slice(info, seq_len(min(10L, nrow(info))))
    listed <- str_c(
        coalesce(shown$study, "-"),
        coalesce(shown$context, "-"),
        coalesce(shown$trait, "-"),
        sep = " / "
    )
    more <- if (nrow(info) > 10L) str_c(", and ", nrow(info) - 10L, " more")
    str_c(str_flatten(listed, "; "), more %||% "")
}

# The data slots the pipeline builds from its inputs. They are real
# colocboost formals, so the union check would accept them -- but the
# pipeline supplies them, so a value given here would be overwritten.
# @noRd
.cbPipelineOwnedArgs <- function() {
    c(
        "X",
        "Y",
        "sumstat",
        "LD",
        "X_ref",
        "focal_outcome_idx",
        "dict_YX",
        "dict_sumstatLD"
    )
}

# @noRd
.cbRefusePipelineOwned <- function(extra) {
    clash <- intersect(names(extra), .cbPipelineOwnedArgs())
    if (length(clash) == 0L) {
        return(invisible(NULL))
    }
    abort(glue(
        "ColocboostOptions: {str_flatten(clash, ', ')} ",
        "{if (length(clash) == 1L) 'is' else 'are'} assembled by ",
        "colocboostPipeline() from its inputs, so a value set here would be ",
        "overwritten. Use the pipeline's own arguments instead."
    ))
}

# =============================================================================
# Helpers (private)
# =============================================================================

# Build the LD / X_ref slot of the colocboost call from a bundle's LD list.
# `ldArg` says which the list holds -- LD matrices ("LD") or factors of them
# ("X_ref"), which colocboost correlates itself -- as .cbMergeSumstatBundles()
# records it; a factor can be square by coincidence, so its shape cannot say.
.cbBuildLdArgs <- function(ldList, ldArg = "LD") {
    ldList <- compact(ldList)
    if (length(ldList) == 0L) {
        return(list())
    }
    set_names(list(ldList), ldArg)
}

# Reject SumStats objects that have not been passed through
# summaryStatsQc(). Both QtlSumStats and GwasSumStats expose qcOptions();
# an empty list (the constructor default) signals "no QC run".
.cbRequireSumStatsQc <- function(x, what) {
    if (is.null(x)) {
        return(invisible(NULL))
    }
    if (length(qcOptions(x)) == 0L) {
        msg <- glue(
            "{what} must be passed through summaryStatsQc() before ",
            "reaching colocboostPipeline (qcOptions() returned an empty ",
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
    pipCutoffToSkip = 0,
    residualizationParam = ResidualizationParam(),
    genotypeFilter,
    phenotypeFilter = .qtlKeepAllPhenotypes()
) {
    contexts <- .cbResolveBundleContexts(qd, contexts)
    p <- list(
        qd = qd,
        genotypeFilter = genotypeFilter,
        phenotypeFilter = phenotypeFilter,
        traitId = traitId,
        region = region,
        cisWindow = cisWindow,
        samples = samples,
        pipCutoffToSkip = pipCutoffToSkip,
        residualizationParam = residualizationParam
    )
    built <- compact(map(contexts, .cbBuildContextXY, p = p))
    if (length(built) == 0L) {
        return(NULL)
    }
    XperCtx <- set_names(map(built, "X"), map_chr(built, "ctx"))
    YperCtx <- set_names(map(built, "Y"), map_chr(built, "ctx"))
    dedup <- .cbDedupX(XperCtx)
    split <- .cbSplitY(YperCtx, dedup$xMatch)
    # An individual-level outcome's LD comes from its own genotypes, so it has
    # no LD-reference entry to name.
    outcomeInfo <- mutate(
        split$outcomeInfo,
        study = studyName(qd),
        dataForm = "individual",
        ldReference = NA_character_
    )
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
        return(contexts(qd))
    }
    available <- contexts(qd)
    bad <- setdiff(contexts, available)
    if (length(bad) > 0L) {
        qdStudy <- studyName(qd)
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

# A `trait` no selected context carries is refused, as fineMappingPipeline and
# twasWeightsPipeline refuse it, rather than every context being skipped into
# an empty result.
# @noRd
.cbAssertTraitsSelected <- function(qd, contexts, traitId) {
    ctxs <- .cbResolveBundleContexts(qd, contexts)
    carried <- map_lgl(ctxs, .cbContextCarries, qd = qd, traitId = traitId)
    if (length(traitId) > 0L && !any(carried)) {
        abort("colocboostPipeline(QtlDataset): no traits selected.")
    }
    invisible(NULL)
}

# The QtlSumStats form: a contexts / trait selection that matches no row.
# @noRd
.cbAssertQssRowsSelected <- function(ss, contexts, traitId) {
    if (nrow(ss) > 0L && !any(.cbQtlRowsSelected(ss, contexts, traitId))) {
        abort(glue(
            "colocboostPipeline(QtlSumStats): no entries matched the ",
            "supplied contexts / traitId filters."
        ))
    }
    invisible(NULL)
}

# The multi-study form: refused only when no embedded dataset's contexts and
# no sumstat row carry a requested trait. A study that lacks it sits out, as
# it does in fineMappingPipeline and twasWeightsPipeline.
# @noRd
.cbMsAssertTraitsSelected <- function(msqd, contexts, traitId) {
    if (length(traitId) == 0L) {
        return(invisible(NULL))
    }
    inData <- map_lgl(
        as.list(msqd),
        .cbDatasetCarries,
        contexts = contexts,
        traitId = traitId
    )
    ss <- sumStats(msqd)
    inSumStats <- !is.null(ss) &&
        nrow(ss) > 0L &&
        any(.cbQtlRowsSelected(ss, contexts, traitId))
    if (!any(inData) && !inSumStats) {
        abort("colocboostPipeline(MultiStudyQtlDataset): no traits selected.")
    }
    invisible(NULL)
}

# Whether any of a dataset's selected contexts carries a requested trait;
# contexts it lacks are left out (a study need not carry every context).
# @noRd
.cbDatasetCarries <- function(qd, contexts, traitId) {
    ctxs <- intersect(contexts %||% contexts(qd), contexts(qd))
    any(map_lgl(ctxs, .cbContextCarries, qd = qd, traitId = traitId))
}

# Whether one context carries any requested trait (every context does when no
# trait was requested).
# @noRd
.cbContextCarries <- function(ctx, qd, traitId) {
    length(traitId) == 0L ||
        any(is_in(traitId, rownames(.qtlContextSe(qd, context = ctx))))
}

# Residualized (X, Y) for one context (sample-aligned + signal-screened), or
# NULL when the context should be skipped.
# @noRd
.cbBuildContextXY <- function(ctx, p) {
    rawY <- .cbResidualizedY(
        p$qd,
        ctx,
        p$traitId,
        p$region,
        p$residualizationParam,
        p$phenotypeFilter
    )
    if (is.null(rawY) || ncol(rawY) == 0L) {
        return(NULL)
    }
    rawX <- .cbResidualizedX(
        p$qd,
        ctx,
        p$traitId,
        p$region,
        p$cisWindow,
        p$samples,
        p$residualizationParam,
        p$genotypeFilter
    )
    if (is.null(rawX) || ncol(rawX) == 0L) {
        return(NULL)
    }
    .cbAlignContextXY(ctx, rawX, rawY, p$pipCutoffToSkip)
}

# One context's X and Y on the samples they share, Y signal-screened; NULL,
# with a message, when they share none.
# @noRd
.cbAlignContextXY <- function(ctx, rawX, rawY, pipCutoffToSkip) {
    # Canonicalize variant colnames (chr-prefix + separator, allele order
    # preserved; rsIDs passed through) so colocboost's name-based matching
    # aligns them with the sumstat / LD ids and across studies. A genuine
    # ref/alt swap stays a distinct id -- names are aligned here, allele
    # *coding* is not.
    named <- `colnames<-`(rawX, normalizeVariantId(colnames(rawX)))
    common <- intersect(rownames(named), rownames(rawY))
    if (length(common) == 0L) {
        msg <- glue(
            "colocboostPipeline: skipping context '{ctx}' ",
            "(no samples shared between residualized X and Y)."
        )
        inform(msg)
        return(NULL)
    }
    X <- named[common, , drop = FALSE]
    Y <- .cbApplyScreen(X, rawY[common, , drop = FALSE], ctx, pipCutoffToSkip)
    if (is.null(Y)) {
        return(NULL)
    }
    list(ctx = ctx, X = X, Y = Y)
}

# colocboost's phenotype QC: colocboost takes missing phenotypes as they are,
# so an unset naAction keeps them, and every choice is available.
# @noRd
.cbResolvePhenotypeFilter <- function(phenotypeFilterParam) {
    .qtlResolvePhenotypeFilter(
        phenotypeFilterParam,
        naDefault = "keep",
        naAllowed = c("keep", "drop", "impute"),
        caller = "colocboostPipeline"
    )
}

# Residualized phenotypes for one context (message + NULL on failure).
# @noRd
.cbResidualizedY <- function(
    qd,
    ctx,
    traitId,
    region,
    residualizationParam = ResidualizationParam(),
    phenotypeFilter = .qtlKeepAllPhenotypes()
) {
    try_fetch(
        .qtlResidualizedPhenotypes(
            qd,
            contexts = ctx,
            traitId = traitId,
            window = region,
            phenotypeFilter = phenotypeFilter,
            residualizationParam = residualizationParam
        )[[ctx]],
        error = function(cnd) {
            msg <- glue(
                "colocboostPipeline: skipping context '{ctx}' ",
                "(residualized phenotypes unavailable)."
            )
            inform(msg, parent = cnd)
            NULL
        }
    )
}

# Residualized genotypes for one context (message + NULL on failure).
# @noRd
.cbResidualizedX <- function(
    qd,
    ctx,
    traitId,
    region,
    cisWindow,
    samples,
    residualizationParam = ResidualizationParam(),
    genotypeFilter
) {
    try_fetch(
        .qtlResidualizedGenotypes(
            qd,
            contexts = ctx,
            window = .qtlResolveVariantRegion(qd, traitId, region, cisWindow),
            samples = samples,
            filter = genotypeFilter,
            residualizationParam = residualizationParam
        ),
        error = function(cnd) {
            msg <- glue(
                "colocboostPipeline: skipping context '{ctx}' ",
                "(residualized genotypes unavailable)."
            )
            inform(msg, parent = cnd)
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

# The index of the first element of `xs` identical to `xs[[i]]`.
# @noRd
.cbFirstIdentical <- function(i, xs) {
    which(map_lgl(xs, identical, xs[[i]]))[[1L]]
}

# Deduplicate X matrices identical across contexts so dict_YX can fan out to a
# smaller X set. Returns list(uniqueX, xMatch) where xMatch[i] is the uniqueX
# index for the i-th context.
#
# Each context resolves to the first context it matches, so nothing has to be
# carried between iterations: the ones that resolve to themselves are the
# unique set, and every context's slot is its representative's position in it.
# @noRd
.cbDedupX <- function(XperCtx) {
    firstIdx <- map_int(seq_along(XperCtx), .cbFirstIdentical, xs = XperCtx)
    keepIdx <- which(firstIdx == seq_along(XperCtx))
    list(
        uniqueX = XperCtx[keepIdx],
        xMatch = match(firstIdx, keepIdx)
    )
}

# Split each context's Y into single-trait columns (context-qualifying duplicate
# trait names) and build the dict_YX (Y-index, X-index) map. Returns
# list(YSplit, dict). (Sequential make.unique naming -- kept as a loop.)
# @noRd
.cbSplitY <- function(YperCtx, xMatch) {
    allTraitNames <- unname(list_c(map(YperCtx, colnames)))
    dupTraits <- unique(allTraitNames[
        duplicated(allTraitNames) | duplicated(allTraitNames, fromLast = TRUE)
    ])
    # Every (context, column) pair, flattened once so the rest is positional.
    pairs <- .cbConcat(map(
        seq_along(YperCtx),
        .cbYPairsForContext,
        YperCtx = YperCtx,
        xMatch = xMatch
    ))
    # The outcome NAME is a display label: context-qualified only when a trait
    # is ambiguous, then made unique. Making it unique depends on the names
    # already minted, so this one step is a fold -- the only sequential part.
    outcomeNames <- reduce(
        pairs,
        .cbAccumOutcomeName,
        dupTraits = dupTraits,
        .init = character(0)
    )
    # The name is lossy -- a bare name could be a trait seen in one context or
    # a trait literally called that -- so the (context, trait) it was minted
    # from is recorded alongside it, the only point where both are in hand.
    list(
        YSplit = set_names(map(pairs, "col"), outcomeNames),
        dict = cbind(Y = seq_along(pairs), X = map_int(pairs, "xIdx")),
        outcomeInfo = bind_rows(map2(outcomeNames, pairs, .cbOutcomeInfoRow))
    )
}

# Concatenate a list of lists into one, empty-safe.
# @noRd
.cbConcat <- function(pieces) {
    if (length(pieces) == 0L) {
        return(list())
    }
    list_c(pieces)
}

# One (context, column) pair: the column itself plus what it was minted from.
# @noRd
.cbYPair <- function(j, Y, ctx, xIdx) {
    list(
        col = Y[, j, drop = FALSE],
        ctx = ctx,
        raw = colnames(Y)[j],
        trait = as.character(colnames(Y)[j] %||% NA_character_),
        xIdx = xIdx
    )
}

# @noRd
.cbYPairsForContext <- function(i, YperCtx, xMatch) {
    Y <- YperCtx[[i]]
    map(
        seq_len(ncol(Y)),
        .cbYPair,
        Y = Y,
        ctx = names(YperCtx)[i],
        xIdx = xMatch[[i]]
    )
}

# @noRd
.cbAccumOutcomeName <- function(existing, pair, dupTraits) {
    c(existing, .cbTraitName(pair$raw, pair$ctx, dupTraits, existing))
}

# @noRd
.cbOutcomeInfoRow <- function(name, pair) {
    tibble(name = name, context = pair$ctx, trait = pair$trait)
}

# Resolve a unique outcome name for a single trait column: default unnamed
# columns, context-qualify duplicates, and make.unique against existing names.
# @noRd
.cbTraitName <- function(tname, ctx, dupTraits, existing) {
    named <- if (is.null(tname) || is.na(tname) || tname == "") {
        str_c("outcome", length(existing) + 1L)
    } else {
        tname
    }
    qualified <- if (is_in(named, dupTraits)) {
        str_c(ctx, "_", named)
    } else {
        named
    }
    if (!is_in(qualified, existing)) {
        return(qualified)
    }
    make.unique(c(existing, qualified))[length(existing) + 1L]
}

# Build a single (sumstat data.frame, LD correlation matrix) pair from a
# QtlSumStats / GwasSumStats entry. Returns NULL when the entry has no
# variants overlapping the ldReference panel.
.cbSumstatPair <- function(
    df,
    ldReference,
    varY = NULL,
    nCase = NULL,
    nControl = NULL,
    cutoffs = NULL
) {
    if (is.null(df) || nrow(df) == 0L) {
        return(NULL)
    }
    df <- mutate(df, variant_id = .cbSumstatVariantIds(df))
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
        ldReference,
        cutoffs,
        ".cbSumstatPair"
    )]
    if (length(wanted) == 0L) {
        return(NULL)
    }
    ld <- .cbPairLd(ldReference, wanted)
    if (is.null(ld)) {
        return(NULL)
    }
    keptIds <- attr(ld, "keptVariantIds")
    kept <- filter(df, is_in(variantIds, keptIds))
    ss <- .cbSumstatFrame(kept, keptIds, nCase, nControl, varY)
    # The reference rides with the pair, so the result can record the one
    # each outcome was matched against.
    list(
        sumstat = ss,
        LD = `attr<-`(ld, "keptVariantIds", NULL),
        variantIds = keptIds,
        isFactor = .ldIsGenotypeRef(ldReference),
        ldReference = ldReference
    )
}

# The LD of the `wanted` variants the panel holds, the rest dropped (the
# colocboost path expects to operate only on the overlap): a genotype
# reference's factor, which colocboost takes as X_ref and correlates itself,
# or an LdMatrix's matrix.
# @noRd
.cbPairLd <- function(ldReference, wanted) {
    ldFrom <- if (.ldIsGenotypeRef(ldReference)) {
        .ldFactorFromRef
    } else {
        .ldFromRef
    }
    ldFrom(ldReference, wanted, label = ".cbSumstatPair", onMissing = "drop")
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
    mutate(
        ss,
        !!!compact(list(
            var_y = if (!is.null(varY) && !is.na(varY)) varY
        ))
    )
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
# optional contexts / traitId filter, with `cisWindow` laid on the tuples it
# keeps. Returns a list keyed by sumstat study label, where each entry has
# (sumstat, LD).
.cbQtlSumStatsBundle <- function(
    ss,
    contexts = NULL,
    traitId = NULL,
    cutoffs = NULL,
    cisWindow = NULL
) {
    if (is.null(ss) || nrow(ss) == 0L) {
        return(list())
    }
    rows <- which(.cbQtlRowsSelected(ss, contexts, traitId))
    if (length(rows) == 0L) {
        return(list())
    }
    ss <- .qssApplyCisWindow(ss, cisWindow, "colocboostPipeline", rows)
    compact(set_names(
        map(rows, .cbQtlEntryPair, ss = ss, cutoffs = cutoffs),
        map_chr(rows, .cbQtlEntryLabel, ss = ss)
    ))
}

# The rows of `ss` colocboost's selection keeps: `contexts` and `traitId` are
# plain vectors here, and an empty one keeps everything.
# @noRd
.cbQtlRowsSelected <- function(ss, contexts, traitId) {
    byContext <- if (!is.null(contexts) && length(contexts) > 0L) {
        is_in(as.character(ss$context), contexts)
    } else {
        rep(TRUE, nrow(ss))
    }
    byTrait <- if (!is.null(traitId) && length(traitId) > 0L) {
        is_in(as.character(ss$trait), traitId)
    } else {
        TRUE
    }
    byContext & byTrait
}

# @noRd
.cbQtlEntryLabel <- function(i, ss) {
    str_c(
        as.character(ss$study)[[i]],
        as.character(ss$context)[[i]],
        as.character(ss$trait)[[i]],
        sep = ":"
    )
}

# One QTL entry's (sumstat, LD) pair, against the entry's own LD reference,
# or NULL when it has no panel overlap. It carries the tuple it came from,
# which a focalOutcome vector is matched against, and the tuple's span, which
# decides the GWAS rows a study contributes (.cbQtlWindow()).
# @noRd
.cbQtlEntryPair <- function(i, ss, cutoffs) {
    tuple <- list(
        study = as.character(ss$study)[[i]],
        context = as.character(ss$context)[[i]],
        trait = as.character(ss$trait)[[i]]
    )
    pair <- .cbSumstatPair(
        df = .sumStatsDf(
            ss,
            studyName = tuple$study,
            context = tuple$context,
            trait = tuple$trait,
            require = "Z"
        ),
        ldReference = .rowLdRef(ss, i),
        varY = .cbTupleValue(ss, "varY", i),
        cutoffs = cutoffs
    )
    if (is.null(pair)) {
        return(NULL)
    }
    list_assign(pair, !!!tuple, span = range(ss[[i]]))
}

# A tuple column's value for row `i`, or NA when the collection lacks it.
# `[[` on one of these collections extracts an ELEMENT, not a column, so this
# goes through mcols() the way the `$` method does.
# @noRd
.cbTupleValue <- function(collection, column, i) {
    if (!is_in(column, .tupleColumnNames(collection))) {
        return(NA_real_)
    }
    S4Vectors::mcols(collection, use.names = FALSE)[[column]][[i]]
}

# Same as .cbQtlSumStatsBundle for a GwasSumStats collection: one pair per
# STUDY, keyed by study name. A study split by chromosome or LD block owns
# several rows, but it is one trait, so it reaches colocboost as one outcome
# -- the union of the rows that overlap `window`, where the QTL data sit
# (.cbQtlWindow()); its other rows are not read. `window = NULL` takes every
# row.
.cbGwasSumStatsBundle <- function(gws, window = NULL, cutoffs = NULL) {
    if (is.null(gws) || nrow(gws) == 0L) {
        return(list())
    }
    index <- .rowSpans(gws)
    hit <- if (is.null(window)) {
        rep(TRUE, length(index$rows))
    } else {
        .overlapsOnShared(index$spans, window)
    }
    rows <- unique(index$rows[hit])
    studies <- as.character(gws$study)
    .cbInformGwasLeftOut(setdiff(unique(studies), studies[rows]))
    compact(imap(
        split(rows, factor(studies[rows], levels = unique(studies[rows]))),
        .cbGwasStudyPair,
        gws = gws,
        cutoffs = cutoffs
    ))
}

# @noRd
.cbInformGwasLeftOut <- function(studies) {
    if (length(studies) == 0L) {
        return(invisible(NULL))
    }
    inform(glue(
        "colocboostPipeline: leaving out GWAS ",
        "{str_flatten(studies, ', ')}: no row overlaps the QTL data."
    ))
}

# One GWAS study's (sumstat, LD) pair over `rows`, its rows the QTL data
# overlap, stitched into one outcome; NULL when it has no panel overlap. The
# rows' statistics are read against one LD matrix, so they must share an LD
# reference. The study-level columns (varY, nCase, nControl) repeat on every
# row, so the first row's stand for the study.
# @noRd
.cbGwasStudyPair <- function(rows, study, gws, cutoffs) {
    pair <- .cbSumstatPair(
        df = .entryToSumstatDf(
            .ssStitchElements(gws, rows),
            require = "Z",
            derive = "none",
            keepChrPrefix = TRUE,
            label = glue("GwasSumStats[{study}]")
        ),
        ldReference = .rowsSharedLdRef(
            gws,
            rows,
            glue(
                "colocboostPipeline: the rows of GWAS study '{study}' the ",
                "QTL data overlap"
            )
        ),
        varY = .cbTupleValue(gws, "varY", rows[[1L]]),
        nCase = .cbTupleValue(gws, "nCase", rows[[1L]]),
        nControl = .cbTupleValue(gws, "nControl", rows[[1L]]),
        cutoffs = cutoffs
    )
    if (is.null(pair)) {
        return(NULL)
    }
    list_assign(pair, study = study)
}

# Combine sumstat bundles into the (sumstat-list, LD-list, dict) shape
# colocboost expects. Deduplicates identical LD matrices so dict_sumstatLD
# can point multiple sumstats at one LD, and records in `ldArg` whether the
# list holds LD matrices ("LD") or factors of them ("X_ref").
#
# The QTL and GWAS sides each read their own LD reference, so a run can mix a
# genotype panel's factors with an LdMatrix's matrices. colocboost takes LD or
# X_ref but not both, so such a run gets every factor as its LD matrix -- the
# matrix colocboost would compute from it anyway once a panel has at least as
# many samples as variants.
.cbMergeSumstatBundles <- function(bundles) {
    if (length(bundles) == 0L) {
        return(list(
            sumstat = list(),
            LD = list(),
            dict_sumstatLD = matrix(integer(0), ncol = 2L),
            ldArg = "LD"
        ))
    }
    isFactor <- map_lgl(bundles, .cbPairIsFactor)
    if (any(isFactor) && !all(isFactor)) {
        bundles <- map(bundles, .cbPairAsLdMatrix)
    }
    # Same shape as .cbDedupX: every bundle resolves to the first bundle
    # carrying an identical LD matrix, so the unique set and the back-mapping
    # both fall out without carrying state.
    lds <- map(bundles, "LD")
    firstIdx <- map_int(seq_along(lds), .cbFirstIdentical, xs = lds)
    keepIdx <- which(firstIdx == seq_along(lds))
    list(
        sumstat = map(bundles, "sumstat"),
        LD = set_names(lds[keepIdx], str_c("LD", seq_along(keepIdx))),
        dict_sumstatLD = cbind(
            sumstat = seq_along(bundles),
            LD = match(firstIdx, keepIdx)
        ),
        ldArg = if (all(map_lgl(bundles, .cbPairIsFactor))) "X_ref" else "LD"
    )
}

# A (sumstat, LD) pair with its LD held as the matrix. A factor F of a
# genotype reference gives crossprod(F) / (nrow(F) - 1), which is that
# reference's LD over the same variants (see .ldFactorFromRef()).
# @noRd
.cbPairAsLdMatrix <- function(p) {
    if (!.cbPairIsFactor(p)) {
        return(p)
    }
    list_assign(
        p,
        LD = crossprod(p$LD) / (nrow(p$LD) - 1L),
        isFactor = FALSE
    )
}

# One focal separate-GWAS run's sumstat bundle: every QTL pair plus the one
# GWAS under test, in that order, so the GWAS is the LAST outcome and its
# focal index is simply the outcome count.
# @noRd
.cbGwasRunBundle <- function(gwasKey, pairs, qtlKeys) {
    .cbMergeSumstatBundles(pairs[c(qtlKeys, gwasKey)])
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

# TRUE when an analysis was requested and has the data to run. When it was
# requested but cannot run, warn and answer FALSE, so the caller's compact()
# drops it instead of skipping it silently.
# @noRd
.cbCanRun <- function(requested, hasData, flag, needed) {
    if (!isTRUE(requested)) {
        return(FALSE)
    }
    if (isTRUE(hasData)) {
        return(TRUE)
    }
    .cbWarnNoData(flag, needed)
    FALSE
}

# A requested analysis with nothing to run on would otherwise be a silent
# no-op: the caller gets an empty ColocBoostResult whose computingTime()
# has no row for it, with no indication of why. Say so.
# @noRd
.cbWarnNoData <- function(flag, needed) {
    msg <- glue(
        "colocboostPipeline: {flag} = TRUE was requested, but there is no ",
        "{needed} to run it on. Skipping it -- the returned ",
        "ColocBoostResult will hold no confidence sets from this analysis ",
        "and its computingTime() will have no row for it."
    )
    warn(msg)
}

# Shared dispatch: accepts a fully-prepared individual bundle (possibly
# NULL) plus a sumstat bundle (possibly empty) and runs the colocboost
# variants `run` (from .cbResolveRuns) switched on. `qtlSumstatBundle` is the
# QTL-side subset of `sumstatBundle` (see .cbDriver); only the xQTL-only run
# uses it, so that a GWAS study never becomes an xQTL-only outcome. `focal`
# names each run's focal outcome (.cbResolveFocal()). The result is keyed by
# the result slot each run fills; compact() drops what did not run.
# @noRd
.cbDispatchRuns <- function(
    individualBundle,
    sumstatBundle,
    qtlSumstatBundle,
    run,
    hasInd,
    focal,
    methodArgs,
    gwasRunBundles = list()
) {
    compact(list(
        xqtl_coloc = if (run$xqtl) {
            .cbRunXqtlOnly(
                individualBundle,
                qtlSumstatBundle,
                hasInd,
                focal$xqtl,
                methodArgs
            )
        },
        joint_gwas = if (run$joint) {
            .cbRunJointGwas(
                individualBundle,
                sumstatBundle,
                hasInd,
                methodArgs,
                gwasNames = names(gwasRunBundles),
                focal = focal$joint
            )
        },
        separate_gwas = if (run$separate) {
            .cbRunSeparateGwas(
                individualBundle,
                gwasRunBundles,
                hasInd,
                methodArgs
            )
        }
    ))
}

# Which analyses were requested AND have data to run on. .cbCanRun warns in
# the "requested but cannot" case so the skip is never silent; the flags are
# resolved in the order the runs are listed, so those warnings keep that
# order. separateGwas needs a GWAS, not merely "some summary statistics": a
# QTL-only sumstat input used to satisfy a `hasSs` guard and then run one
# focal pass per QTL study.
# @noRd
.cbResolveRuns <- function(
    xqtlColoc,
    jointGwas,
    separateGwas,
    hasInd,
    hasSs,
    qtlSumstatBundle,
    gwasRunBundles
) {
    hasQtlSs <- length(qtlSumstatBundle$sumstat) > 0L
    list(
        xqtl = .cbCanRun(
            xqtlColoc,
            hasInd || hasQtlSs,
            "xqtlColoc",
            "QTL data"
        ),
        joint = .cbCanRun(
            jointGwas,
            hasSs,
            "jointGwas",
            "summary-statistic data"
        ),
        separate = .cbCanRun(
            separateGwas,
            length(gwasRunBundles) > 0L,
            "separateGwas",
            "GWAS summary statistics"
        )
    )
}

# Each run's result and timing folded into the empty result skeleton, so a
# slot no analysis filled keeps the skeleton's own (empty) value.
# @noRd
.cbMergeRuns <- function(empty, runs) {
    list_assign(
        list_assign(empty, !!!map(runs, "result")),
        computing_time = list_assign(
            empty$computing_time,
            Analysis = list_assign(
                empty$computing_time$Analysis,
                !!!map(runs, "time")
            )
        )
    )
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

# Run the requested colocboost analyses over the prepared bundles and fold
# them into a ColocBoostResult (an empty one when no QTL input remains).
# `focal` names each run's focal outcome; `resultInfo` carries the outcome
# table and LD references the result records (.cbPrepareInputs()).
# @noRd
.cbRunVariants <- function(
    individualBundle,
    sumstatBundle,
    xqtlColoc,
    jointGwas,
    separateGwas,
    focal,
    methodArgs,
    resultInfo = list(outcomeInfo = .cbEmptyOutcomeInfo(), ldRefs = list()),
    qtlSumstatBundle = NULL,
    gwasRunBundles = list()
) {
    hasInd <- !is.null(individualBundle)
    hasSs <- length(sumstatBundle$sumstat) > 0L
    qtlSumstatBundle <- qtlSumstatBundle %||% .cbMergeSumstatBundles(list())
    if (!hasInd && !hasSs) {
        msg <- glue(
            "colocboostPipeline: no QTL inputs remain after selection. ",
            "Nothing to run."
        )
        inform(msg)
        return(.cbEmptyResultObject())
    }
    run <- .cbResolveRuns(
        xqtlColoc,
        jointGwas,
        separateGwas,
        hasInd,
        hasSs,
        qtlSumstatBundle,
        gwasRunBundles
    )
    runs <- .cbDispatchRuns(
        individualBundle,
        sumstatBundle,
        qtlSumstatBundle,
        run,
        hasInd,
        focal,
        methodArgs,
        gwasRunBundles = gwasRunBundles
    )
    .cbToResultObject(
        .cbMergeRuns(.cbEmptyResult(), runs),
        resultInfo$outcomeInfo,
        resultInfo$ldRefs
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

# @noRd
.cbJointRun <- function(nm, results) {
    if (is.null(results[[nm]])) {
        return(NULL)
    }
    list(run = results[[nm]], analysis = nm, gwasStudy = NA_character_)
}

# @noRd
.cbSeparateRun <- function(i, sep, keys) {
    if (is.null(sep[[i]])) {
        return(NULL)
    }
    list(run = sep[[i]], analysis = "separate_gwas", gwasStudy = keys[[i]])
}

# One character field across run records, empty-safe.
# @noRd
.cbEntryField <- function(entries, field) {
    if (length(entries) == 0L) {
        return(character(0))
    }
    map_chr(entries, field)
}

# Flatten the three raw colocboost runs into one ColocBoostResult.
#
# The three analyses are NOT parallel structures: xqtl_coloc and joint_gwas are
# each a single colocboost object, while separate_gwas is a NAMED LIST of them,
# one per GWAS study -- one level deeper. That is why `gwasStudy` exists as a
# key: without it the separate-GWAS sets would be indistinguishable from each
# other once flattened.
# @noRd
.cbToResultObject <- function(results, outcomeInfo, ldReference = list()) {
    sep <- results$separate_gwas
    sepRuns <- if (is.null(sep) || length(sep) == 0L) {
        list()
    } else {
        compact(map(
            seq_along(sep),
            .cbSeparateRun,
            sep = sep,
            keys = names(sep) %||% as.character(seq_along(sep))
        ))
    }
    # One record per run, so the three parallel vectors are read off a single
    # list instead of being grown in step with each other.
    entries <- c(
        compact(map(
            c("xqtl_coloc", "joint_gwas"),
            .cbJointRun,
            results = results
        )),
        sepRuns
    )
    ColocBoostResult(
        results = map(entries, "run"),
        analysis = .cbEntryField(entries, "analysis"),
        gwasStudy = .cbEntryField(entries, "gwasStudy"),
        outcomeInfo = outcomeInfo,
        ldReference = ldReference,
        computingTime = results$computing_time %||% list()
    )
}

# The xQTL-only run's progress line, which names the summary-statistic
# studies only when there are any.
# @noRd
.cbAnnounceXqtl <- function(nCtx, nSs) {
    msg <- if (nSs > 0L) {
        glue(
            "====== Performing xQTL-only ColocBoost on {nCtx} contexts ",
            "and {nSs} summary-statistic studies. ====="
        )
    } else {
        glue(
            "====== Performing xQTL-only ColocBoost on {nCtx} contexts. ====="
        )
    }
    inform(msg)
}

# xQTL-only ColocBoost run -> list(result, time).
#
# Either side may be absent: `individualBundle` is NULL for a summary-level
# QTL input, and `sumstatBundle` is empty for an individual-level one. It holds
# the QTL-side sumstats ONLY -- a GWAS study must never be pulled into the
# xQTL-only analysis, which is why this does not take the merged bundle the
# joint / separate runs use.
# @noRd
.cbRunXqtlOnly <- function(
    individualBundle,
    sumstatBundle,
    hasInd,
    focal,
    methodArgs
) {
    traits <- c(
        if (hasInd) individualBundle$outcomeNames else character(),
        names(sumstatBundle$sumstat)
    )
    focalIdx <- .cbFocalIdx(focal, traits)
    nCtx <- if (hasInd) length(individualBundle$Y) else 0L
    nSs <- length(sumstatBundle$sumstat)
    .cbAnnounceXqtl(nCtx, nSs)
    ldArgs <- if (nSs > 0L) {
        .cbBuildLdArgs(sumstatBundle$LD, sumstatBundle$ldArg)
    } else {
        list()
    }
    args <- c(
        list(
            X = if (hasInd) individualBundle$X else NULL,
            Y = if (hasInd) individualBundle$Y else NULL,
            dict_YX = if (hasInd) individualBundle$dict_YX else NULL,
            sumstat = if (nSs > 0L) sumstatBundle$sumstat else NULL,
            dict_sumstatLD = if (nSs > 0L) {
                sumstatBundle$dict_sumstatLD
            } else {
                NULL
            },
            outcome_names = traits,
            focal_outcome_idx = focalIdx,
            output_level = 2
        ),
        ldArgs,
        # as.list(): c() would append the MethodOptions record as one opaque
        # element instead of splicing its entries, and exec() would then hand
        # that object to colocboost as a positional argument.
        as.list(methodArgs)
    )
    run <- .cbRun("xQTL-only ColocBoost", args)
    list(result = run$result, time = run$time)
}

# The (outcome name -> study / context / trait / dataForm) lookup for one run.
#
# Outcome names reach colocboost as bare labels; this is what lets a
# ColocBoostResult answer "which contexts colocalized" instead of handing the
# caller a string to parse. `ldKeys` names, per sumstat outcome, the result's
# LD-reference entry it was matched against (.cbResultLdRefs()), and
# `sumstatInfo` (.cbPairsInfo()) the tuple each sumstat outcome came from.
# @noRd
.cbOutcomeInfo <- function(
    individualBundle,
    sumstatBundle,
    hasInd,
    ldKeys = character(0),
    sumstatInfo = NULL
) {
    ssNames <- names(sumstatBundle$sumstat)
    parts <- compact(list(
        if (isTRUE(hasInd)) individualBundle$outcomeInfo,
        if (length(ssNames) > 0L) {
            .cbSumstatOutcomeInfo(ssNames, ldKeys, sumstatInfo)
        }
    ))
    if (length(parts) == 0L) {
        return(.cbEmptyOutcomeInfo())
    }
    bind_rows(parts)
}

# The summary-statistic outcomes' rows. A QTL tuple records its study,
# context and trait; a GWAS study has no context or trait of its own, which
# is recorded as NA rather than invented. Without `sumstatInfo` each outcome
# stands as its own study.
# @noRd
#' @importFrom dplyr coalesce
.cbSumstatOutcomeInfo <- function(ssNames, ldKeys, sumstatInfo) {
    info <- sumstatInfo %||% .cbPairsInfo(set_names(list(), character(0)))
    hit <- match(ssNames, info$name)
    tibble(
        name = ssNames,
        context = info$context[hit],
        trait = info$trait[hit],
        study = coalesce(info$study[hit], ssNames),
        dataForm = "sumstats",
        ldReference = unname(ldKeys[ssNames])
    )
}

# The (name, study, context, trait) of each summary-statistic pair, as the
# pair records them: a QTL tuple's three, a GWAS study's study alone.
# @noRd
.cbPairsInfo <- function(pairs) {
    tibble(
        name = names(pairs) %||% character(0),
        study = map_chr(unname(pairs), "study", .default = NA_character_),
        context = map_chr(unname(pairs), "context", .default = NA_character_),
        trait = map_chr(unname(pairs), "trait", .default = NA_character_)
    )
}

# @noRd
.cbEmptyOutcomeInfo <- function() {
    tibble(
        name = character(0),
        context = character(0),
        trait = character(0),
        study = character(0),
        dataForm = character(0),
        ldReference = character(0)
    )
}

# Joint QTL + GWAS run -> list(result, time). Focal on `focal`, an outcome
# name, when one was resolved (.cbResolveFocal()); non-focal otherwise.
# @noRd
.cbRunJointGwas <- function(
    individualBundle,
    sumstatBundle,
    hasInd,
    methodArgs,
    gwasNames = character(),
    focal = NULL
) {
    ssNames <- names(sumstatBundle$sumstat)
    traits <- c(
        if (hasInd) individualBundle$outcomeNames else character(),
        ssNames
    )
    focalIdx <- .cbFocalIdx(focal, traits)
    ldArgs <- .cbBuildLdArgs(sumstatBundle$LD, sumstatBundle$ldArg)
    # The merged bundle holds QTL sumstats as well as GWAS, so split the
    # labels rather than calling every sumstat a GWAS: a context is any
    # outcome that is not one of the GWAS studies, whether it arrived as
    # individual-level data or as summary statistics.
    nGwas <- length(intersect(ssNames, gwasNames))
    inform(.cbJointMessage(length(traits) - nGwas, nGwas, traits[focalIdx]))
    args <- c(
        list(
            X = if (hasInd) individualBundle$X else NULL,
            Y = if (hasInd) individualBundle$Y else NULL,
            sumstat = sumstatBundle$sumstat,
            dict_YX = if (hasInd) individualBundle$dict_YX else NULL,
            dict_sumstatLD = sumstatBundle$dict_sumstatLD,
            outcome_names = traits,
            focal_outcome_idx = focalIdx,
            output_level = 2
        ),
        ldArgs,
        # as.list(): c() would append the MethodOptions record as one opaque
        # element instead of splicing its entries, and exec() would then hand
        # that object to colocboost as a positional argument.
        as.list(methodArgs)
    )
    run <- .cbRun("Joint GWAS ColocBoost", args)
    list(result = run$result, time = run$time)
}

# The joint run's announcement; `focal` is the focal outcome's name, or empty
# for a non-focal run.
# @noRd
.cbJointMessage <- function(nContexts, nGwas, focal) {
    if (length(focal) == 0L) {
        return(glue(
            "====== Performing non-focal GWAS-xQTL ColocBoost on ",
            "{nContexts} contexts and {nGwas} GWAS. ====="
        ))
    }
    glue(
        "====== Performing focal GWAS-xQTL ColocBoost on {nContexts} ",
        "contexts and {nGwas} GWAS, focal on {focal}. ====="
    )
}

# The position of outcome `focal` among `traits`, the form colocboost's
# focal_outcome_idx takes; NULL when there is no focal outcome or the run
# does not hold it.
# @noRd
.cbFocalIdx <- function(focal, traits) {
    idx <- if (is.null(focal)) NA_integer_ else match(focal, traits)
    if (is.na(idx)) NULL else idx
}

# Separate (focal) per-GWAS runs -> list(result, time) with per-study timing.
# @noRd
.cbRunSeparateGwas <- function(
    individualBundle,
    gwasRunBundles,
    hasInd,
    methodArgs
) {
    # One run per GWAS study -- NOT per merged sumstat label. Each bundle
    # already holds the QTL side plus its own GWAS.
    ssNames <- names(gwasRunBundles)
    t1 <- Sys.time()
    separate <- set_names(
        map(
            ssNames,
            .cbSeparateGwasAt,
            gwasRunBundles = gwasRunBundles,
            individualBundle = individualBundle,
            hasInd = hasInd,
            methodArgs = methodArgs
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
    studyName,
    runBundle,
    individualBundle,
    hasInd,
    methodArgs
) {
    # `runBundle` carries every QTL sumstat followed by this one GWAS, so the
    # whole LD list and its dict come along -- the same shape as the joint
    # run, but focal on the last outcome.
    ldArgs <- .cbBuildLdArgs(runBundle$LD, runBundle$ldArg)
    traits <- c(
        if (hasInd) individualBundle$outcomeNames else character(),
        names(runBundle$sumstat)
    )
    # Every outcome but the focal GWAS is a context, whether it arrived as
    # individual-level data or as summary statistics. Counting only
    # individualBundle$Y reported "0 contexts" for a summary-statistic QTL
    # run -- which was accurate about the call, because the QTL sumstats were
    # not being passed at all.
    nContexts <- length(traits) - 1L
    msg <- glue(
        "====== Performing focal GWAS-xQTL ColocBoost on {nContexts} ",
        "contexts and {studyName} GWAS. ====="
    )
    inform(msg)
    args <- c(
        list(
            X = if (hasInd) individualBundle$X else NULL,
            Y = if (hasInd) individualBundle$Y else NULL,
            sumstat = runBundle$sumstat,
            dict_YX = if (hasInd) individualBundle$dict_YX else NULL,
            dict_sumstatLD = runBundle$dict_sumstatLD,
            outcome_names = traits,
            focal_outcome_idx = length(traits),
            output_level = 2
        ),
        ldArgs,
        # as.list(): c() would append the MethodOptions record as one opaque
        # element instead of splicing its entries, and exec() would then hand
        # that object to colocboost as a positional argument.
        as.list(methodArgs)
    )
    .cbRun(str_c("Separate GWAS ColocBoost for ", studyName), args)$result
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
    # sgn is +1 / -1 per matched column, so one sweep applies every flip.
    out <- sweep(m[, ia, drop = FALSE], 2L, sgn, FUN = "*")
    `colnames<-`(out, canonical[ib])
}

# Relabel a (sumstat, LD) pair to the canonical coding: flip the z-score and the
# LD sign-submatrix for variants swapped relative to canonical
# (LD_ij -> sign_i * sign_j * LD_ij; on a factor, the swapped columns are
# negated), drop unmatched variants, and relabel to canonical. The sumstat
# and its LD share one sign vector, so they stay consistent. Returns NULL
# when no variant matches.
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
    ss <- mutate(
        p$sumstat[ia, , drop = FALSE],
        z = .data$z * sgn,
        variant = canonical[ib]
    )
    isFactor <- .cbPairIsFactor(p)
    ld <- if (isFactor) {
        `colnames<-`(
            sweep(p$LD[, ia, drop = FALSE], 2L, sgn, FUN = "*"),
            canonical[ib]
        )
    } else {
        `dimnames<-`(
            p$LD[ia, ia, drop = FALSE] * outer(sgn, sgn),
            list(canonical[ib], canonical[ib])
        )
    }
    list(
        sumstat = ss,
        LD = ld,
        variantIds = canonical[ib],
        isFactor = isFactor,
        ldReference = p$ldReference
    )
}

# TRUE when a (sumstat, LD) pair carries a factor of its LD rather than the
# matrix.
# @noRd
.cbPairIsFactor <- function(p) {
    isTRUE(p$isFactor)
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
    allIds <- c(
        if (!is.null(individualBundle)) {
            unname(list_c(map(individualBundle$X, colnames)))
        },
        unname(list_c(map(pairs, list("sumstat", "variant"))))
    )
    ids <- unique(allIds[!is.na(allIds)])
    if (length(ids) == 0L) {
        return(list(individualBundle = individualBundle, pairs = pairs))
    }
    parsed <- parseVariantId(ids)
    ok <- !is.na(parsed$chrom) & !is.na(parsed$pos)
    locus <- if_else(ok, str_c(parsed$chrom, parsed$pos, sep = ":"), ids)
    canonical <- ids[!duplicated(locus)]
    flipped <- if (is.null(individualBundle)) {
        NULL
    } else {
        list_assign(
            individualBundle,
            X = map(
                individualBundle$X,
                .cbFlipMatrixToCanonical,
                canonical = canonical
            )
        )
    }
    list(
        individualBundle = flipped,
        pairs = compact(
            map(pairs, .cbFlipPairToCanonical, canonical = canonical)
        )
    )
}

# Harmonize allele coding across all sources to a shared per-locus
# canonical so swapped variants are combined with a consistent sign
# (alleleFlip = TRUE); alleleFlip = FALSE leaves the names-only
# canonicalization done at the source builders, which keeps swapped
# variants distinct.
# @noRd
.cbMaybeHarmonize <- function(individualBundle, combinedPairs, alleleFlip) {
    if (!isTRUE(alleleFlip)) {
        return(list(individualBundle = individualBundle, pairs = combinedPairs))
    }
    .cbHarmonizeAlleles(individualBundle, combinedPairs)
}

# Top-level driver shared by all input methods. qtlPairs and gwasPairs
# are per-tuple lists of `list(sumstat, LD)` produced by the per-class
# bundle helpers; they are merged here so dict_sumstatLD can dedupe
# identical LD matrices across QTL and GWAS sides.
#
# The xQTL-only run gets its own bundle over just the QTL-side pairs, so a
# GWAS study is never treated as an xQTL outcome. Rebuilding it through
# .cbMergeSumstatBundles (rather than subsetting the merged one) keeps the
# deduplicated LD list and dict_sumstatLD consistent for the subset.
# @noRd
.cbDriver <- function(
    individualBundle,
    qtlPairs,
    gwasSumStats,
    xqtlColoc,
    jointGwas,
    separateGwas,
    focalOutcome,
    methodArgs,
    alleleFlip = TRUE,
    cutoffs = NULL
) {
    if (!isTRUE(xqtlColoc) && !isTRUE(jointGwas) && !isTRUE(separateGwas)) {
        inform("colocboostPipeline: no analysis flag is TRUE; nothing to do.")
        return(.cbEmptyResultObject())
    }
    inputs <- .cbPrepareInputs(
        individualBundle,
        qtlPairs,
        gwasSumStats,
        focalOutcome,
        alleleFlip,
        cutoffs
    )
    .cbRunVariants(
        inputs$individualBundle,
        inputs$sumstatBundle,
        xqtlColoc,
        jointGwas,
        separateGwas,
        inputs$focal,
        methodArgs,
        resultInfo = inputs$resultInfo,
        qtlSumstatBundle = inputs$qtlSumstatBundle,
        gwasRunBundles = inputs$gwasRunBundles
    )
}

# Everything the runs read, from the QTL side and gwasSumStats: each
# outcome's (sumstat, LD) pair -- a GWAS study once, over its rows the QTL
# data overlap -- allele-harmonized when asked; the bundles each run takes;
# the outcome each run is focal on; and the outcome table and LD references
# the result records. The tuples and keys are read before harmonization,
# which keeps every pair's name.
# @noRd
.cbPrepareInputs <- function(
    individualBundle,
    qtlPairs,
    gwasSumStats,
    focalOutcome,
    alleleFlip,
    cutoffs
) {
    pairs <- .cbAppendGwasPairs(
        qtlPairs,
        gwasSumStats,
        window = .cbQtlWindow(individualBundle, qtlPairs),
        cutoffs = cutoffs,
        qtlVariants = .cbQtlVariants(individualBundle, qtlPairs)
    )
    harmonized <- .cbMaybeHarmonize(individualBundle, pairs, alleleFlip)
    sumstatBundle <- .cbMergeSumstatBundles(harmonized$pairs)
    qtlKeys <- intersect(names(harmonized$pairs), names(qtlPairs))
    ldRefs <- .cbResultLdRefs(harmonized$pairs)
    list(
        individualBundle = harmonized$individualBundle,
        sumstatBundle = sumstatBundle,
        qtlSumstatBundle = .cbMergeSumstatBundles(harmonized$pairs[qtlKeys]),
        gwasRunBundles = .cbGwasRunBundles(harmonized$pairs, qtlKeys),
        focal = .cbResolveFocal(
            focalOutcome,
            gwasSumStats,
            .cbQtlOutcomeInfo(individualBundle, qtlPairs),
            .cbGwasKeys(pairs, names(qtlPairs))
        ),
        resultInfo = list(
            outcomeInfo = .cbOutcomeInfo(
                harmonized$individualBundle,
                sumstatBundle,
                !is.null(harmonized$individualBundle),
                ldRefs$keys,
                .cbPairsInfo(pairs)
            ),
            ldRefs = ldRefs$refs
        )
    )
}

# Where the QTL data sit, one span per chromosome: the GWAS rows a study
# contributes are those overlapping it. Individual-level
# variants are placed by their ids, which the panel normalizes to chrom:pos
# form (an id that does not parse, such as an rsID, places nothing), and
# summary-level ones by their tuple's range. NULL when nothing is placed,
# which leaves every GWAS row in.
# @noRd
.cbQtlWindow <- function(individualBundle, qtlPairs) {
    ids <- unique(list_c(
        map(individualBundle$X, colnames),
        ptype = character()
    ))
    parsed <- parseVariantId(ids)
    ok <- !is.na(parsed$chrom) & !is.na(parsed$pos)
    placed <- exec(
        c,
        GenomicRanges::GRanges(
            parsed$chrom[ok],
            IRanges::IRanges(parsed$pos[ok], width = 1L)
        ),
        !!!unname(map(compact(map(qtlPairs, "span")), .locusRanges))
    )
    if (length(placed) == 0L) {
        return(NULL)
    }
    range(placed)
}

# The QTL variants a run reads, placed by their ids, each position once: the
# individual-level ones and those each summary-statistic pair kept. NULL when
# none places (an rsID places nothing).
# @noRd
.cbQtlVariants <- function(individualBundle, qtlPairs) {
    ids <- unique(c(
        list_c(map(individualBundle$X, colnames), ptype = character()),
        list_c(map(unname(qtlPairs), "variantIds"), ptype = character())
    ))
    parsed <- parseVariantId(ids)
    ok <- !is.na(parsed$chrom) & !is.na(parsed$pos)
    if (!any(ok)) {
        return(NULL)
    }
    unique(GenomicRanges::GRanges(
        parsed$chrom[ok],
        IRanges::IRanges(parsed$pos[ok], width = 1L)
    ))
}

# colocboost reads every QTL outcome together, so the QTL data are one
# region, checked against the LD blocks of the GWAS studies they overlap
# (.checkGwasBlockReach()).
# @noRd
.cbCheckBlockReach <- function(qtlVariants, gwasSumStats) {
    if (is.null(qtlVariants)) {
        return(invisible(NULL))
    }
    .checkGwasBlockReach(
        list(`the QTL data` = qtlVariants),
        gwasSumStats,
        "colocboostPipeline",
        "the QTL data"
    )
}

# Each GWAS study that reached the runs, mapped to its outcome name.
# @noRd
.cbGwasKeys <- function(pairs, qtlNames) {
    gwas <- pairs[setdiff(names(pairs), qtlNames)]
    set_names(names(gwas), map_chr(unname(gwas), "study"))
}

# The QTL outcomes a focalOutcome vector is matched against: the
# individual-level ones with the (study, context, trait) they were built
# from, then each QTL summary-statistic tuple.
# @noRd
.cbQtlOutcomeInfo <- function(individualBundle, qtlPairs) {
    ind <- individualBundle$outcomeInfo %||% .cbEmptyOutcomeInfo()
    bind_rows(
        select(ind, all_of(c("name", "study", "context", "trait"))),
        .cbPairsInfo(qtlPairs)
    )
}

# The LD references a result records -- each distinct one the summary
# statistics were matched against, once -- and each summary-statistic
# outcome's key into them, read off the pairs that carry them.
# @noRd
.cbResultLdRefs <- function(pairs) {
    dict <- .ldRefDict(map(unname(pairs), "ldReference"))
    list(refs = dict$refs, keys = set_names(dict$keys, names(pairs)))
}

# The separate-GWAS run bundles, keyed by GWAS label. The runs are FOCAL: each
# is one GWAS tested against the whole QTL side, so each needs its own bundle
# of (every QTL pair + that one GWAS). Built through .cbMergeSumstatBundles
# for the same reason as the xQTL bundle in .cbDriver -- the deduplicated LD
# list and dict_sumstatLD have to match the subset, which subsetting the
# merged bundle would not give. Without this the separate loop ran over EVERY
# merged label and treated each QTL study as a focal GWAS in its own right.
# @noRd
.cbGwasRunBundles <- function(pairs, qtlKeys) {
    gwasKeys <- setdiff(names(pairs), qtlKeys)
    set_names(
        map(gwasKeys, .cbGwasRunBundle, pairs = pairs, qtlKeys = qtlKeys),
        gwasKeys
    )
}

# Append the GWAS sumstat pairs -- one per study, over its rows that overlap
# `window` (.cbGwasSumStatsBundle()) -- to the QTL pairs, QC-gating
# gwasSumStats and make.unique-ing colliding labels. The GWAS side's LD comes
# from its own reference, which need not be the QTL side's: each pair carries
# the LD it was matched against.
# @noRd
.cbAppendGwasPairs <- function(
    qtlPairs,
    gwasSumStats,
    window = NULL,
    cutoffs = NULL,
    qtlVariants = NULL
) {
    if (is.null(gwasSumStats)) {
        return(qtlPairs)
    }
    .cbRequireSumStatsQc(gwasSumStats, "gwasSumStats")
    .cbCheckBlockReach(qtlVariants, gwasSumStats)
    gwasPairs <- .cbGwasSumStatsBundle(
        gwasSumStats,
        window = window,
        cutoffs = cutoffs
    )
    if (length(gwasPairs) == 0L) {
        return(qtlPairs)
    }
    # Disambiguating one label at a time against the names so far gives the
    # same answer as one make.unique() over the whole sequence, since it too
    # renames left to right and leaves the first occurrence alone.
    keys <- make.unique(c(names(qtlPairs), names(gwasPairs)))
    c(
        qtlPairs,
        set_names(gwasPairs, keys[length(qtlPairs) + seq_along(gwasPairs)])
    )
}

# =============================================================================
# Methods
# =============================================================================

# The QtlDataset colocboost run: screen spec, individual-level bundle, shared
# driver. Split out so the method body stays short; the arguments are passed
# explicitly rather than bundled, so a rename is a compile-time error and the
# dependency of each helper is visible.
# @noRd
.cbQtlDatasetDrive <- function(
    qtlData,
    gwasSumStats,
    contexts,
    traitId,
    region,
    cisWindow,
    focalOutcome,
    xqtlColoc,
    jointGwas,
    separateGwas,
    samples,
    panelFilterParam,
    signalScreenParam,
    alleleFlip,
    methodArgs,
    genotypeFilterParam,
    phenotypeFilterParam,
    residualizationParam
) {
    screenSpec <- .screenResolve(signalScreenParam)
    genotypeFilter <- .qtlResolveFilter(genotypeFilterParam)
    indBundle <- .cbIndividualBundle(
        .qtlRestrictSamples(qtlData, genotypeFilter$keepSamples),
        contexts = contexts,
        traitId = traitId,
        region = region,
        cisWindow = cisWindow,
        samples = samples,
        pipCutoffToSkip = screenSpec,
        residualizationParam = residualizationParam,
        genotypeFilter = genotypeFilter,
        phenotypeFilter = .cbResolvePhenotypeFilter(phenotypeFilterParam)
    )
    .cbDriver(
        indBundle,
        qtlPairs = list(),
        gwasSumStats,
        xqtlColoc,
        jointGwas,
        separateGwas,
        focalOutcome,
        methodArgs,
        alleleFlip = alleleFlip,
        # The QTL side is individual-level here, but the GWAS side may still
        # be sumstats, so the panel cutoffs still apply to it.
        cutoffs = .panelCutoffs(panelFilterParam)
    )
}


# The genome build colocboostPipeline()'s inputs share; NA when none names
# one, an error when two disagree.
# @noRd
.cbBuild <- function(...) {
    .inputsBuild(list(...), "colocboostPipeline")
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
        focalOutcome = NULL,
        xqtlColoc = TRUE,
        jointGwas = FALSE,
        separateGwas = FALSE,
        samples = NULL,
        panelFilterParam = PanelFilterParam(),
        signalScreenParam = SignalScreenParam(),
        alleleFlip = TRUE,
        genotypeFilterParam = GenotypeFilterParam(),
        phenotypeFilterParam = PhenotypeFilterParam(),
        residualizationParam = ResidualizationParam(),
        methodArgs = ColocboostOptions()
    ) {
        .assertMethodOptions(methodArgs, "ColocboostOptions", "methodArgs")
        .cbAssertAnyRun(xqtlColoc, jointGwas, separateGwas)
        .cbAssertFocalOutcome(focalOutcome)
        .cbAssertTraitsSelected(qtlData, contexts, traitId)
        build <- .cbBuild(qtlData, gwasSumStats, region)
        res <- .cbQtlDatasetDrive(
            qtlData = qtlData,
            gwasSumStats = gwasSumStats,
            contexts = contexts,
            traitId = traitId,
            region = region,
            cisWindow = cisWindow,
            focalOutcome = focalOutcome,
            xqtlColoc = xqtlColoc,
            jointGwas = jointGwas,
            separateGwas = separateGwas,
            samples = samples,
            panelFilterParam = panelFilterParam,
            signalScreenParam = signalScreenParam,
            alleleFlip = alleleFlip,
            methodArgs = methodArgs,
            genotypeFilterParam = genotypeFilterParam,
            phenotypeFilterParam = phenotypeFilterParam,
            residualizationParam = residualizationParam
        )
        .withGenome(res, build)
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
        cisWindow = NULL,
        focalOutcome = NULL,
        xqtlColoc = TRUE,
        jointGwas = FALSE,
        separateGwas = FALSE,
        alleleFlip = TRUE,
        panelFilterParam = PanelFilterParam(),
        methodArgs = ColocboostOptions()
    ) {
        .assertMethodOptions(methodArgs, "ColocboostOptions", "methodArgs")
        .cbAssertAnyRun(xqtlColoc, jointGwas, separateGwas)
        .cbAssertFocalOutcome(focalOutcome)
        .cbRequireSumStatsQc(qtlData, "qtlData")
        .cbAssertQssRowsSelected(qtlData, contexts, traitId)
        build <- .cbBuild(qtlData, gwasSumStats)
        cutoffs <- .panelCutoffs(panelFilterParam)
        qtlPairs <- .cbQtlSumStatsBundle(
            qtlData,
            contexts = contexts,
            traitId = traitId,
            cutoffs = cutoffs,
            cisWindow = cisWindow
        )
        res <- .cbDriver(
            individualBundle = NULL,
            qtlPairs = qtlPairs,
            gwasSumStats = gwasSumStats,
            xqtlColoc = xqtlColoc,
            jointGwas = jointGwas,
            separateGwas = separateGwas,
            focalOutcome = focalOutcome,
            methodArgs = methodArgs,
            alleleFlip = alleleFlip,
            cutoffs = cutoffs
        )
        .withGenome(res, build)
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
        focalOutcome = NULL,
        xqtlColoc = TRUE,
        jointGwas = FALSE,
        separateGwas = FALSE,
        samples = NULL,
        panelFilterParam = PanelFilterParam(),
        signalScreenParam = SignalScreenParam(),
        alleleFlip = TRUE,
        genotypeFilterParam = GenotypeFilterParam(),
        phenotypeFilterParam = PhenotypeFilterParam(),
        residualizationParam = ResidualizationParam(),
        methodArgs = ColocboostOptions()
    ) {
        .assertMethodOptions(methodArgs, "ColocboostOptions", "methodArgs")
        .cbAssertAnyRun(xqtlColoc, jointGwas, separateGwas)
        .cbAssertFocalOutcome(focalOutcome)
        # As in fineMappingPipeline: a context no study has is refused, and
        # one only some lack is skipped for those, with a warning.
        parseContexts(contexts, qtlData)
        .cbMsAssertTraitsSelected(qtlData, contexts, traitId)
        build <- .cbBuild(qtlData, gwasSumStats, region)
        res <- .cbPipelineMultiStudy(
            qtlData = qtlData,
            gwasSumStats = gwasSumStats,
            contexts = contexts,
            traitId = traitId,
            region = region,
            cisWindow = cisWindow,
            focalOutcome = focalOutcome,
            xqtlColoc = xqtlColoc,
            jointGwas = jointGwas,
            separateGwas = separateGwas,
            samples = samples,
            panelFilterParam = panelFilterParam,
            signalScreenParam = signalScreenParam,
            alleleFlip = alleleFlip,
            methodArgs = methodArgs,
            genotypeFilterParam = genotypeFilterParam,
            phenotypeFilterParam = phenotypeFilterParam,
            residualizationParam = residualizationParam
        )
        .withGenome(res, build)
    }
)

# MultiStudyQtlDataset colocboost worker: aggregate the per-study individual
# bundles + embedded sumstats, then dispatch to the shared driver.
# @noRd
.cbPipelineMultiStudy <- function(
    qtlData,
    gwasSumStats,
    contexts,
    traitId,
    region,
    cisWindow,
    focalOutcome,
    xqtlColoc,
    jointGwas,
    separateGwas,
    samples,
    panelFilterParam,
    signalScreenParam,
    alleleFlip,
    methodArgs,
    genotypeFilterParam,
    phenotypeFilterParam,
    residualizationParam
) {
    cutoffs <- .panelCutoffs(panelFilterParam)
    qtlSide <- .cbMultiStudyQtlSide(
        qtlData,
        contexts = contexts,
        traitId = traitId,
        region = region,
        cisWindow = cisWindow,
        samples = samples,
        screenSpec = .screenResolve(signalScreenParam),
        residualizationParam = residualizationParam,
        genotypeFilter = .qtlResolveFilter(genotypeFilterParam),
        phenotypeFilter = .cbResolvePhenotypeFilter(phenotypeFilterParam),
        cutoffs = cutoffs
    )
    .cbDriver(
        qtlSide$indBundle,
        qtlSide$qtlPairs,
        gwasSumStats,
        xqtlColoc,
        jointGwas,
        separateGwas,
        focalOutcome,
        methodArgs,
        alleleFlip = alleleFlip,
        cutoffs = cutoffs
    )
}

# The QTL side of a MultiStudyQtlDataset run: the individual-level bundle
# across its studies, and the summary-statistic pairs alongside them.
# @noRd
.cbMultiStudyQtlSide <- function(
    qtlData,
    contexts,
    traitId,
    region,
    cisWindow,
    samples,
    screenSpec,
    residualizationParam,
    genotypeFilter,
    phenotypeFilter,
    cutoffs
) {
    list(
        indBundle = .cbMultiStudyIndBundle(
            qtlData = qtlData,
            contexts = contexts,
            traitId = traitId,
            region = region,
            cisWindow = cisWindow,
            samples = samples,
            screenSpec = screenSpec,
            residualizationParam = residualizationParam,
            genotypeFilter = genotypeFilter,
            phenotypeFilter = phenotypeFilter
        ),
        qtlPairs = .cbMultiStudySumstats(
            qtlData,
            contexts,
            traitId,
            cisWindow,
            cutoffs = cutoffs
        )$qtlPairs
    )
}

# One study's bundle with its names study-prefixed, or NULL when it has
# nothing to contribute. A `contexts` vector applies to every study: each
# runs on the ones it has, and a study with none of them sits out (the entry
# point has warned about it).
# @noRd
.cbStudyBundle <- function(
    studyName,
    qtlDatasets,
    contexts,
    traitId,
    region,
    cisWindow,
    samples,
    screenSpec,
    residualizationParam = ResidualizationParam(),
    genotypeFilter,
    phenotypeFilter = .qtlKeepAllPhenotypes()
) {
    qd <- qtlDatasets[[studyName]]
    share <- if (is.null(contexts)) NULL else intersect(contexts, contexts(qd))
    if (!is.null(share) && length(share) == 0L) {
        return(NULL)
    }
    sub <- .cbIndividualBundle(
        .qtlRestrictSamples(qd, genotypeFilter$keepSamples),
        contexts = share,
        traitId = traitId,
        region = region,
        cisWindow = cisWindow,
        samples = samples,
        pipCutoffToSkip = screenSpec,
        residualizationParam = residualizationParam,
        genotypeFilter = genotypeFilter,
        phenotypeFilter = phenotypeFilter
    )
    if (is.null(sub)) {
        return(NULL)
    }
    .cbPrefixStudyNames(sub, studyName)
}

# @noRd
.cbFieldLength <- function(sub, field) {
    length(sub[[field]])
}

# How many entries precede each element, given per-element counts.
# @noRd
.cbOffsets <- function(counts) {
    cumsum(c(0L, counts))[seq_along(counts)]
}

# Study `i`'s dict, renumbered into the combined X / Y ordering.
# @noRd
.cbShiftStudyDict <- function(i, subs, xOffsets, yOffsets) {
    d <- subs[[i]]$dict_YX
    cbind(Y = d[, "Y"] + yOffsets[[i]], X = d[, "X"] + xOffsets[[i]])
}

# The individual-level bundle across every study of a MultiStudyQtlDataset.
# Each study's dict indexes its own X and Y, so it shifts by however many
# entries the studies before it contributed -- a cumulative count, known once
# all the bundles are in hand rather than tracked while looping.
# NULL when no study has anything to contribute.
# @noRd
.cbMultiStudyIndBundle <- function(
    qtlData,
    contexts,
    traitId,
    region,
    cisWindow,
    samples,
    screenSpec,
    residualizationParam = ResidualizationParam(),
    genotypeFilter,
    phenotypeFilter = .qtlKeepAllPhenotypes()
) {
    qtlDatasets <- as.list(qtlData)
    subs <- compact(map(
        names(qtlDatasets),
        .cbStudyBundle,
        qtlDatasets = qtlDatasets,
        contexts = contexts,
        traitId = traitId,
        region = region,
        cisWindow = cisWindow,
        samples = samples,
        screenSpec = screenSpec,
        residualizationParam = residualizationParam,
        genotypeFilter = genotypeFilter,
        phenotypeFilter = phenotypeFilter
    ))
    if (length(subs) == 0L) {
        return(NULL)
    }
    .cbCombineStudyBundles(subs)
}

# The per-study bundles as one: X, Y and the outcomes concatenated, and each
# study's dict shifted past the entries the studies before it contributed.
# NULL when no study contributes any X.
# @noRd
.cbCombineStudyBundles <- function(subs) {
    xOffsets <- .cbOffsets(map_int(subs, .cbFieldLength, field = "X"))
    yOffsets <- .cbOffsets(map_int(subs, .cbFieldLength, field = "Y"))
    combinedX <- .cbConcat(map(subs, "X"))
    if (length(combinedX) == 0L) {
        return(NULL)
    }
    list(
        X = combinedX,
        Y = .cbConcat(map(subs, "Y")),
        dict_YX = exec(
            rbind,
            !!!map(
                seq_along(subs),
                .cbShiftStudyDict,
                subs = subs,
                xOffsets = xOffsets,
                yOffsets = yOffsets
            )
        ),
        outcomeNames = .cbConcat(map(subs, "outcomeNames")),
        outcomeInfo = bind_rows(map(subs, "outcomeInfo"))
    )
}

# Prefix a study's X names + outcome names with "{study}:" (Y names track the
# outcome names), so colocboost sees distinct outcomes when two studies share
# a trait.
# @noRd
.cbPrefixStudyNames <- function(sub, studyName) {
    outcomeNames <- str_c(studyName, sub$outcomeNames, sep = ":")
    renamed <- list_assign(
        sub,
        X = set_names(sub$X, str_c(studyName, names(sub$X), sep = ":")),
        Y = set_names(sub$Y, outcomeNames),
        outcomeNames = outcomeNames
    )
    # The lookup table is keyed on the outcome NAME, so it has to be renamed in
    # the same pass -- a stale key silently drops every outcome of this study
    # from the identity join.
    if (is.null(sub$outcomeInfo) || nrow(sub$outcomeInfo) == 0L) {
        return(renamed)
    }
    # `studyLabel` so the data mask cannot shadow it with the `study` column.
    studyLabel <- studyName
    list_assign(
        renamed,
        outcomeInfo = mutate(
            sub$outcomeInfo,
            name = str_c(studyLabel, .data$name, sep = ":"),
            study = studyLabel
        )
    )
}

# Sumstat side of a MultiStudyQtlDataset: bundle any embedded QtlSumStats,
# with `cisWindow` laid on the tuples selected as a QtlSumStats run lays it.
# Returns list(qtlPairs); each pair carries its own row's LD reference.
# @noRd
.cbMultiStudySumstats <- function(
    qtlData,
    contexts,
    traitId,
    cisWindow,
    cutoffs = NULL
) {
    embeddedSs <- sumStats(qtlData)
    if (is.null(embeddedSs)) {
        return(list(qtlPairs = list()))
    }
    .cbRequireSumStatsQc(embeddedSs, "MultiStudyQtlDataset@sumStats")
    list(
        qtlPairs = .cbQtlSumStatsBundle(
            embeddedSs,
            contexts = contexts,
            traitId = traitId,
            cutoffs = cutoffs,
            cisWindow = cisWindow
        )
    )
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# Single-effect screen decision for outcome column `j` (drop on failure).
# @noRd
.cbScreenOutcome <- function(j, X, Y, spec) {
    .fmSerScreen(X, Y[, j], spec, fallback = FALSE)
}

# Run one separate-GWAS colocboost fit for GWAS study `gwasName`.
# @noRd
.cbSeparateGwasAt <- function(
    gwasName,
    gwasRunBundles,
    individualBundle,
    hasInd,
    methodArgs
) {
    .cbRunOneSeparateGwas(
        gwasName,
        gwasRunBundles[[gwasName]],
        individualBundle,
        hasInd,
        methodArgs
    )
}

#' @title colocboost Engine Settings
#' @description Additional arguments forwarded to
#'   \code{colocboost::colocboost}. Names are checked against that
#'   function's own formals, so a misspelled option is an error here rather
#'   than a setting silently dropped inside the run.
#'
#'   The data arguments \code{colocboostPipeline} assembles itself ---
#'   \code{X}, \code{Y}, \code{sumstat}, \code{LD}, \code{X_ref},
#'   \code{focal_outcome_idx} --- are not settable here.
#' @param ... Named \code{colocboost::colocboost} arguments.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' ColocboostOptions(M = 5)
#' @export
ColocboostOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            X = "supplied from the data by the pipeline",
            Y = "supplied from the data by the pipeline",
            dict_YX = "supplied from the data by the pipeline",
            dict_sumstatLD = "supplied from the data by the pipeline",
            sumstat = "supplied from the data by the pipeline",
            LD = "supplied from the data by the pipeline",
            outcome_names = "supplied from the data by the pipeline",
            focal_outcome_idx = "supplied from the data by the pipeline",
            output_level = "supplied from the data by the pipeline"
        ),
        "ColocboostOptions"
    )
    .cbRefusePipelineOwned(extra)
    .newMethodOptions(
        "colocboost::colocboost",
        defaults = list(),
        extra = extra,
        label = "ColocboostOptions",
        engine = "colocboost"
    )
}

# Run colocboost() with tryCatch + timing.
#' @importFrom rlang try_fetch
.cbRun <- function(label, args) {
    # Engine-call exception: colocboost's whole interface is two blocks --
    # this and ColocboostOptions -- so a colocboostWrapper.R would hold ~50
    # lines against this file's 1900. Kept inline deliberately.
    if (!requireNamespace("colocboost", quietly = TRUE)) {
        abort("The colocboost package is required for colocboostPipeline().")
    }
    t1 <- Sys.time()
    # Callers splice the flattened bundle in with the assembled data slots;
    # drop the slots that resolved to NULL before handing them over.
    args <- compact(args)
    res <- try_fetch(
        exec(colocboost::colocboost, !!!args),
        error = function(cnd) {
            msg <- glue("{label} failed")
            inform(msg, parent = cnd)
            NULL
        }
    )
    list(result = res, time = Sys.time() - t1)
}
