# =============================================================================
# Pipeline-stage and data-filtering argument bundles
# -----------------------------------------------------------------------------
# Settings that pecotmr itself reads, as opposed to the MethodConfig records
# that forward options to an outside engine: every constructor here is built
# with a NULL callee. The one exception is susieRssControlConfig, which does
# wrap susieR::susie_rss_control() and lives here because it is the nested
# value of rssConfig(control =).
#
# Three groups, on axes that genuinely cross -- genotypes are filtered
# whether the run fine-maps or fits TWAS weights, and credible-set settings
# apply whichever filters ran:
#
#   STAGE: how a stage behaves, keyed on where in the pipeline it acts
#     signalScreenConfig     -- before any fit: is this block worth fitting?
#     residualizationConfig  -- before the fit: what is regressed out
#     rssConfig              -- during the fit: the summary-statistics solver
#     credibleSetConfig      -- after the fit: how credible sets are built
#                               and reported
#     ensembleConfig         -- after the fits: how they are stacked
#
#   FILTER: which rows and columns of the data survive, keyed on the object
#     genotypeFilterConfig, panelFilterConfig, sumstatsFilterConfig,
#     sumstatsCleaningConfig -- see that section's header below
#
#   SHARED: crossValidationConfig, one bundle both fineMappingPipeline() and
#     twasWeightsPipeline() take, spanning the other two axes -- it is a
#     stage setting and it carries its own variant cap
#
# `L` / `Lgreedy` are credibleSetConfig fields: they bound how many credible
# sets can exist. They reach the engine by a different route from the rest of
# the bundle -- .fmSeededArgNames() is c("L", "L_greedy"), seeded onto each
# SuSiE-family token's own arguments, whereas the other fields are read by
# postprocessFinemappingFits() -- so changing them needs a refit while a
# stored fit can be re-post-processed under a different coverage.
#
# `addSusieInf` is a standalone argument, not a bundle field: it selects a
# chained initialisation between methods rather than describing credible sets.
#
# `fitRetention` is a plain arg_match() enum on each pipeline, not a bundle:
# see colocPipeline(), crossValidation() and fineMappingPipeline().
# =============================================================================

#' @include AllGenerics.R
NULL

#' @title Signal Pre-Screen Settings
#' @description Whether to skip a block before fitting it, and on which
#'   metric. Shared by \code{\link{fineMappingPipeline}} and
#'   \code{\link{colocboostPipeline}}.
#' @section One metric at a time:
#'   The four metrics are mutually exclusive and the constructor enforces it,
#'   so a conflicting pair fails where it is written rather than part-way
#'   into a pipeline run. Leave them all unset (the default) to screen
#'   nothing.
#' @param pip Posterior-inclusion-probability cutoff. A negative value uses
#'   the adaptive \code{3 / nVariants} threshold.
#' @param absZ Cutoff on the maximum absolute z-score.
#' @param bf Cutoff on the maximum per-variant Bayes factor.
#' @param logBf Cutoff on the maximum per-variant log Bayes factor.
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' signalScreenConfig(absZ = 5)
#' @export
signalScreenConfig <- function(
    pip = NULL,
    absZ = NULL,
    bf = NULL,
    logBf = NULL
) {
    res <- .newMethodConfig(
        NULL,
        defaults = list(pip = pip, absZ = absZ, bf = bf, logBf = logBf),
        extra = list(),
        label = "signalScreenConfig",
        engine = "signalScreen"
    )
    .screenAssertOneMetric(res)
    .screenAssertPositive(res)
    res
}

# absZ and bf have no meaningful negative cutoff: absZ screens on max|Z| and
# a Bayes factor is positive. (pip keeps `< 0 => 3 / nVariants`, and logBf is
# a log, so both may legitimately be negative.)
# @noRd
.screenAssertPositive <- function(screen) {
    for (m in c("absZ", "bf")) {
        v <- screen[[m]]
        if (!is.null(v) && length(v) == 1L && !is.na(v) && v < 0) {
            abort(glue(
                "signalScreenConfig: `{m}` must be > 0, got {v}. ",
                "absZ screens on max|Z| and Bayes factors are positive."
            ))
        }
    }
    invisible(NULL)
}

# Refuse more than one enabled metric. Enabled means set AND non-zero: 0 is
# the long-standing "off" spelling for these cutoffs, so signalScreenConfig(pip
# = 0, absZ = 5) is one screen, not two.
# @noRd
.screenAssertOneMetric <- function(screen) {
    on <- names(screen)[map_lgl(as.list(screen), .screenIsOn)]
    if (length(on) > 1L) {
        abort(glue(
            "signalScreenConfig: only one screening metric may be enabled at ",
            "a time; got {str_flatten(on, ', ')}."
        ))
    }
    invisible(NULL)
}

# @noRd
.screenIsOn <- function(x) {
    !is.null(x) && length(x) > 0L && any(as.numeric(x) != 0, na.rm = TRUE)
}

# The one enabled metric as the polymorphic screen spec the pipelines already
# understand (see .asScreen): a bare numeric for `pip`, a list(metric, cutoff)
# for the others, NULL when nothing is enabled. A field set to 0 is off, so
# signalScreenConfig(pip = 0, absZ = 5) resolves to the absZ screen.
# @noRd
.screenResolve <- function(screen) {
    on <- names(screen)[map_lgl(as.list(screen), .screenIsOn)]
    if (length(on) == 0L) {
        return(NULL)
    }
    metric <- on[[1L]]
    if (metric == "pip") {
        return(screen$pip)
    }
    list(metric = metric, cutoff = screen[[metric]])
}

#' @title Credible-Set Construction And Reporting
#' @description How credible sets are built from a fit and which of them are
#'   reported. Every field is read by
#'   \code{\link{postprocessFinemappingFits}}, so these can be changed and a
#'   stored fit re-summarized without refitting.
#' @param coverage Primary credible-set coverage. Default \code{0.95}.
#' @param secondaryCoverage Additional coverages to report credible sets at.
#'   Default \code{c(0.7, 0.5)}.
#' @param signalCutoff PIP cutoff for including a non-credible-set variant in
#'   the top-loci table. Default \code{0.025}.
#' @param minAbsCorr Minimum absolute within-set correlation for purity.
#'   Default \code{0.8}.
#' @param medianAbsCorr Median absolute within-set correlation threshold, or
#'   \code{NULL} (default) to judge purity on \code{minAbsCorr} alone.
#' @param includeAllCs Logical. Report every credible set rather than only the
#'   top one. Default \code{FALSE}.
#' @param perCsColumns Which per-credible-set variant-level columns are added
#'   to the \code{topLoci} table. \code{"none"} (default) adds only the
#'   always-on scalar \code{within_cs_pip}; \code{"alpha"} also widens
#'   \code{alpha} into one \code{within_cs_pip_<lab>} column per set;
#'   \code{"full"} additionally widens \code{cs_logbf_}, \code{cs_effect_}
#'   and \code{cs_effect_var_}. \code{includeAllCs} decides the \code{<lab>}
#'   in those names. This governs the \emph{table}; the stored fit is
#'   governed by \code{fitRetention}.
#' @param L Integer. Maximum number of single effects the fit may carry, and
#'   so the maximum number of credible sets that can exist. Default
#'   \code{10}. Seeded onto every SuSiE-family token in \code{methods} that
#'   did not set it through \code{\link{fineMappingMethodsConfig}}.
#' @param Lgreedy Integer or \code{NULL}. Maximum number of single effects
#'   for the greedy initialization stage, where the engine has one.
#'   \code{NULL} (the default) leaves the engine's own default in place.
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' credibleSetConfig(coverage = 0.9, includeAllCs = TRUE)
#' @export
credibleSetConfig <- function(
    coverage = 0.95,
    secondaryCoverage = c(0.7, 0.5),
    signalCutoff = 0.025,
    minAbsCorr = 0.8,
    medianAbsCorr = NULL,
    includeAllCs = FALSE,
    perCsColumns = c("none", "alpha", "full"),
    L = 10L,
    Lgreedy = NULL
) {
    perCsColumns <- arg_match(perCsColumns)
    .newMethodConfig(
        NULL,
        defaults = list(
            coverage = coverage,
            secondaryCoverage = secondaryCoverage,
            signalCutoff = signalCutoff,
            minAbsCorr = minAbsCorr,
            medianAbsCorr = medianAbsCorr,
            includeAllCs = includeAllCs,
            perCsColumns = perCsColumns,
            L = L,
            Lgreedy = Lgreedy
        ),
        extra = list(),
        label = "credibleSetConfig",
        engine = "credibleSet"
    )
}

#' @title Covariate Residualization Settings
#' @description What is regressed out of the phenotype and genotype before
#'   fitting. Shared by \code{\link{fineMappingPipeline}} and
#'   \code{\link{twasWeightsPipeline}}, which hold the same four settings.
#'
#'   \code{fineMappingPipeline}'s \code{usePCA} / \code{nPCs} are
#'   \strong{not} here. They do not residualize anything: they PCA-reduce a
#'   multi-trait context's phenotype matrix and fine-map each top principal
#'   component \emph{as a trait}. That is an analysis mode, a sibling of the
#'   univariate and multivariate dispatch paths, not a covariate setting.
#' @section Summary-statistics inputs:
#'   A \code{QtlSumStats} / \code{GwasSumStats} input carries no genotypes
#'   or covariates, so nothing here applies and the whole bundle is
#'   \strong{ignored} --- not refused.
#'
#'   That is the opposite of \code{\link{crossValidationConfig}}, which
#'   \emph{is} refused on those inputs
#'   (\code{fineMappingPipeline(QtlSumStats, crossValidation = ...)} errors).
#'   The difference is the default: CV is off unless asked for
#'   (\code{folds = 0}), so a non-default value on a summary-statistics run
#'   is an explicit request for something impossible. Residualization is on
#'   by default (\code{residualizePhenotype} and \code{residualizeGenotype}
#'   are both \code{TRUE}), so refusing a non-default would reject the
#'   \emph{default} bundle and force every sumstats caller to unset it.
#'
#'   One bundle therefore travels to either input kind unchanged, the same
#'   way \code{crossValidationConfig}'s \code{weightMethods} /
#'   \code{maxVariants} are carried but ignored by
#'   \code{fineMappingPipeline}.
#' @param phenotypeCovariates Covariates to residualize the phenotype on, or
#'   \code{NULL} (default) for the dataset's own.
#' @param genotypeCovariates Covariates to residualize the genotype on, or
#'   \code{NULL} (default) for the dataset's own.
#' @param residualizePhenotype Logical. Residualize the phenotype. Default
#'   \code{TRUE}.
#' @param residualizeGenotype Logical. Residualize the genotype. Default
#'   \code{TRUE}.
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' residualizationConfig(residualizeGenotype = FALSE)
#' @export
residualizationConfig <- function(
    phenotypeCovariates = NULL,
    genotypeCovariates = NULL,
    residualizePhenotype = TRUE,
    residualizeGenotype = TRUE
) {
    .newMethodConfig(
        NULL,
        defaults = list(
            phenotypeCovariates = phenotypeCovariates,
            genotypeCovariates = genotypeCovariates,
            residualizePhenotype = residualizePhenotype,
            residualizeGenotype = residualizeGenotype
        ),
        extra = list(),
        label = "residualizationConfig",
        engine = "residualization"
    )
}

#' @title Summary-Statistics Solver Settings
#' @description How the SuSiE-RSS fit is conditioned and what it falls back
#'   to. Applies to \code{QtlSumStats} / \code{GwasSumStats} input.
#' @param serFallback Logical. Fall back to a single-effect fit when the
#'   multi-effect fit is judged unreliable. Default \code{FALSE}.
#' @param keepFullFit Which pre-fallback multi-effect fits to retain:
#'   \code{"fallback"} (default) keeps them only for regions that fell back,
#'   \code{"all"} for every region, \code{"none"} for none. Only meaningful
#'   when \code{serFallback} is \code{TRUE}.
#' @param rFinite Sample size used for the LD matrix's finiteness check, or
#'   \code{NULL} (default) to use the panel's own.
#' @param rMismatch How to treat an LD / z-score mismatch. Default
#'   \code{"none"}.
#' @param control Options for \code{susieR::susie_rss_control()}, forwarded as
#'   \code{susie_rss()}'s \code{control} argument. Built with
#'   \code{\link{susieRssControlConfig}}; \code{NULL} (default) leaves that
#'   function's own defaults in place.
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' rssConfig(serFallback = TRUE, keepFullFit = "all")
#' @export
rssConfig <- function(
    serFallback = FALSE,
    keepFullFit = c("fallback", "all", "none"),
    rFinite = NULL,
    rMismatch = "none",
    control = NULL
) {
    keepFullFit <- arg_match(keepFullFit)
    if (!is.null(control)) {
        .assertMethodConfig(control, "susieRssControlConfig", "control")
    }
    .newMethodConfig(
        NULL,
        defaults = list(
            serFallback = serFallback,
            keepFullFit = keepFullFit,
            rFinite = rFinite,
            rMismatch = rMismatch,
            control = control
        ),
        extra = list(),
        label = "rssConfig",
        engine = "rss"
    )
}

#' @title Arguments For susieR's RSS Control Block
#' @description Options forwarded to \code{susieR::susie_rss_control()} and
#'   from there to \code{susie_rss()}'s \code{control} argument. Names are
#'   checked against that function's live formals, so a misspelling fails at
#'   the call site rather than being dropped into an ignored list --- which is
#'   what a bare named list here used to do.
#' @param ... Arguments for \code{susieR::susie_rss_control()}, under its own
#'   names (\code{check_prior}, \code{mismatch_estimator}, ...).
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' susieRssControlConfig(check_prior = TRUE)
#' @export
susieRssControlConfig <- function(...) {
    .newMethodConfig(
        "susieR::susie_rss_control",
        defaults = list(),
        extra = list(...),
        label = "susieRssControlConfig",
        engine = "susieRssControl"
    )
}

#' @title SR-TWAS Ensemble Settings
#' @description Whether and how \code{\link{twasWeightsPipeline}} stacks its
#'   per-method weights into an SR-TWAS ensemble.
#' @section Requires cross-validation:
#'   Stacking combines each method's \strong{out-of-fold} predictions, so it
#'   cannot run without cross-validation. \code{enabled = TRUE} together with
#'   \code{crossValidationConfig(folds < 2)} is an error. It used to be
#'   neither: the ensemble row was simply absent from the result, with
#'   \code{ensemble = TRUE} still reading as on.
#'
#'   \code{enabled} defaults to \code{FALSE} because
#'   \code{\link{crossValidationConfig}} defaults to no folds. To get an
#'   ensemble, ask for both.
#' @param enabled Logical. Compute SR-TWAS ensemble weights. Default
#'   \code{FALSE}.
#' @param r2Threshold Minimum cross-validated \eqn{R^2} for a method to enter
#'   the stack. Default \code{0.01}. Stacking needs at least two methods to
#'   clear it.
#' @param solver Stacking solver, \code{"quadprog"} (default) or
#'   \code{"glmnet"}.
#' @param alpha Elastic-net mixing parameter, used only when
#'   \code{solver = "glmnet"}. Default \code{1}.
#' @return A \code{\link{MethodConfig}} object.
#' @seealso \code{\link{crossValidationConfig}}
#' @examples
#' ensembleConfig(enabled = TRUE, r2Threshold = 0.05)
#' @export
ensembleConfig <- function(
    enabled = FALSE,
    r2Threshold = 0.01,
    solver = c("quadprog", "glmnet"),
    alpha = 1
) {
    solver <- arg_match(solver)
    .newMethodConfig(
        NULL,
        defaults = list(
            enabled = enabled,
            r2Threshold = r2Threshold,
            solver = solver,
            alpha = alpha
        ),
        extra = list(),
        label = "ensembleConfig",
        engine = "ensemble"
    )
}

# Refuse an ensemble that cannot be built. Stacking reads out-of-fold
# predictions, so without folds there is nothing to stack -- and the old
# behaviour was to return a result silently missing its ensemble row.
# @noRd
.ensembleAssertCv <- function(ensembleArgs, crossValidationArgs) {
    if (!isTRUE(ensembleArgs$enabled) || .cvEnabled(crossValidationArgs)) {
        return(invisible(NULL))
    }
    abort(glue(
        "twasWeightsPipeline: ensembleConfig(enabled = TRUE) needs ",
        "out-of-fold predictions to stack, so it requires ",
        "crossValidationConfig(folds >= 2); got folds = ",
        "{crossValidation$folds %||% 0}."
    ))
}

# susieR's `control` argument is a plain named list, so the constructor
# result is flattened on the way out. NULL is preserved rather than becoming
# list(): to susie_rss() an absent control means "use susie_rss_control()'s
# own defaults", which an empty list does not.
# @noRd
.rssControlList <- function(control) {
    if (is.null(control) || length(control) == 0L) {
        return(NULL)
    }
    as.list(control)
}

# =============================================================================
# Variant / sample filtering bundles
# -----------------------------------------------------------------------------
# Which rows and columns of the DATA survive, as opposed to the stage bundles
# above, which describe how a stage behaves. Three filters run in this
# package, on three different objects, and they are deliberately three
# constructors rather than one:
#
#   genotypeFilterConfig  -- a study's own genotype matrix (QtlDataset and the
#                          pipelines that build one)
#   panelFilterConfig     -- an LD reference panel's variants
#   sumstatsFilterConfig  -- rows of a summary-statistics table
#
# Their fields overlap without meaning the same thing, and their defaults
# genuinely differ: `imissCutoff` resolves to 0 on the genotype path and 1 on
# the panel path. One union bundle would have to pick one of those, and would
# accept `removeIndels` where nothing reads it -- the failure mode these
# constructors exist to prevent.
#
# genotypeFilterConfig defaults every field to NULL ("not set") because its
# fields are a SPECIFICATION to QtlDataset and an OVERRIDE to the pipelines;
# absence is what distinguishes "pin this value" from "leave it alone". The
# other two have one meaning each and keep ordinary defaults.
#
# sumstatsCleaningConfig sits with them because it runs on the same pass over
# the same table, but it is coercion and normalisation (coerceNumeric,
# normalizeChr, clampSmallP) that also drops rows -- not a filter.
# =============================================================================

#' @title Genotype Filtering Options
#' @description Which variants and samples survive when a genotype matrix is
#'   assembled, as one checked bundle. Used by \code{\link{QtlDataset}}, the
#'   manifest loaders, and the pipelines that build a dataset.
#' @section Unset versus set:
#'   Every field defaults to \code{NULL}, meaning \strong{not set}, and an
#'   unset field is absent from the result rather than carried as \code{NULL}.
#'   That is what lets one bundle serve two roles:
#'   \itemize{
#'     \item to \code{\link{QtlDataset}} it is a \emph{specification}, and
#'       an unset field takes that constructor's own default --- \code{0} for
#'       each cutoff, \code{character(0)} for each \code{keep}, \code{TRUE}
#'       for \code{keepIndel}.
#'     \item to a pipeline it is an \emph{override}, and an unset field
#'       leaves the dataset's construct-time value alone.
#'   }
#'   So \code{genotypeFilterConfig(mafCutoff = 0)} and
#'   \code{genotypeFilterConfig()} differ: the first pins the cutoff at zero,
#'   the second defers. \code{\link{panelFilterConfig}} needs no such
#'   distinction and keeps ordinary defaults.
#' @param mafCutoff Minor-allele-frequency floor.
#' @param macCutoff Minor-allele-count floor; the stricter of this and
#'   \code{mafCutoff} applies.
#' @param xvarCutoff Genotype-variance floor.
#' @param imissCutoff Per-variant missingness ceiling.
#' @param keepSamples Sample ids to restrict to; \code{character(0)} keeps
#'   all samples.
#' @param keepVariants Variant ids to restrict to; \code{character(0)} keeps
#'   all variants.
#' @param keepIndel Logical. Retain insertions and deletions.
#' @return A \code{\link{MethodConfig}} object.
#' @seealso \code{\link{panelFilterConfig}}, \code{\link{sumstatsFilterConfig}}
#' @examples
#' genotypeFilterConfig(mafCutoff = 0.01, keepIndel = FALSE)
#' @export
genotypeFilterConfig <- function(
    mafCutoff = NULL,
    macCutoff = NULL,
    xvarCutoff = NULL,
    imissCutoff = NULL,
    keepSamples = NULL,
    keepVariants = NULL,
    keepIndel = NULL
) {
    .newMethodConfig(
        NULL,
        defaults = list(
            mafCutoff = mafCutoff,
            macCutoff = macCutoff,
            xvarCutoff = xvarCutoff,
            imissCutoff = imissCutoff,
            keepSamples = keepSamples,
            keepVariants = keepVariants,
            keepIndel = keepIndel
        ),
        extra = list(),
        label = "genotypeFilterConfig",
        engine = "genotypeFilter"
    )
}

#' @title LD Reference Panel Filtering Options
#' @description Which of an LD panel's variants are kept before the panel is
#'   used, as one checked bundle. Narrower than
#'   \code{\link{genotypeFilterConfig}}: a panel is not a study's genotypes, so
#'   there is no sample restriction and no variance cutoff, and the
#'   missingness default is \code{1} (no filter) rather than \code{0}.
#' @param mafCutoff Minor-allele-frequency floor. Default \code{0}.
#' @param macCutoff Minor-allele-count floor; the stricter of this and
#'   \code{mafCutoff} applies, using the panel's own sample count. Default
#'   \code{0}.
#' @param imissCutoff Per-variant missingness ceiling. Default \code{1}, which
#'   filters nothing and lets the allele-frequency sidecar be read instead of
#'   materializing dosage.
#' @return A \code{\link{MethodConfig}} object.
#' @seealso \code{\link{genotypeFilterConfig}},
#'   \code{\link{sumstatsFilterConfig}}
#' @examples
#' panelFilterConfig(mafCutoff = 0.001)
#' @export
panelFilterConfig <- function(
    mafCutoff = 0,
    macCutoff = 0,
    imissCutoff = 1
) {
    .newMethodConfig(
        NULL,
        defaults = list(
            mafCutoff = mafCutoff,
            macCutoff = macCutoff,
            imissCutoff = imissCutoff
        ),
        extra = list(),
        label = "panelFilterConfig",
        engine = "panelFilter"
    )
}

#' @title Summary-Statistic Row Filtering Options
#' @description Which rows of a summary-statistics table survive QC, as one
#'   checked bundle. These act on the sumstats themselves, not on any genotype
#'   matrix, which is why they are separate from
#'   \code{\link{genotypeFilterConfig}}.
#' @param removeIndels Logical. Drop insertions and deletions. Default
#'   \code{FALSE}.
#' @param removeStrandAmbiguous Logical. Drop strand-ambiguous variants (A/T
#'   and C/G). Default \code{TRUE}.
#' @param infoCutoff Imputation-INFO floor. Default \code{0}.
#' @param nCutoff Per-variant sample-size floor. Default \code{5}.
#' @return A \code{\link{MethodConfig}} object.
#' @seealso \code{\link{genotypeFilterConfig}}, \code{\link{panelFilterConfig}}
#' @examples
#' sumstatsFilterConfig(removeIndels = TRUE, infoCutoff = 0.8)
#' @export
sumstatsFilterConfig <- function(
    removeIndels = FALSE,
    removeStrandAmbiguous = TRUE,
    infoCutoff = 0,
    nCutoff = 5
) {
    .newMethodConfig(
        NULL,
        defaults = list(
            removeIndels = removeIndels,
            removeStrandAmbiguous = removeStrandAmbiguous,
            infoCutoff = infoCutoff,
            nCutoff = nCutoff
        ),
        extra = list(),
        label = "sumstatsFilterConfig",
        engine = "sumstatsFilter"
    )
}

#' @title Summary-Statistics Cleaning Settings
#' @description What makes a summary-statistics row a well-formed record:
#'   type coercion, chromosome-label normalization, dropping malformed or
#'   impossible rows, and flooring underflowed p-values.
#'
#'   Separate from \code{\link{sumstatsFilterConfig}} on purpose. That one
#'   applies quality thresholds to \emph{valid} records (INFO, N, indels,
#'   strand-ambiguous variants); this one decides whether a row is a valid
#'   record at all. All of it is applied by \code{\link{summaryStatsQc}}
#'   before any filter or screen runs.
#' @param coerceNumeric Logical. Coerce the signed columns
#'   (Z/BETA/SE/OR/LOG_ODDS/SIGNED_SUMSTAT/P/MAF/FRQ/INFO/N) to numeric.
#'   Default \code{TRUE}.
#' @param normalizeChr Logical. Strip the \code{"chr"} prefix, uppercase the
#'   chromosome label, and map 23->X, 24->Y, M->MT. Default \code{TRUE}.
#' @param dropNonstandardChr Logical. Drop variants whose CHR (after
#'   normalization) is outside 1..22, X, Y, MT. Default \code{TRUE}.
#' @param dropMissData Logical. Drop rows with NA in any vital column (chrom,
#'   pos, A1, A2, and at least one of Z / BETA). Default \code{TRUE}.
#' @param dropPOutOfRange Logical. Drop rows where \code{P < 0} or
#'   \code{P > 1}. Default \code{TRUE}.
#' @param clampSmallP Logical. Floor non-negative P values below
#'   \code{smallPFloor} to \code{smallPFloor}, so \code{-log10(P)} stays
#'   finite. Applied to both input and Z-derived P values. Default
#'   \code{TRUE}.
#' @param smallPFloor Numeric (length 1). The floor \code{clampSmallP}
#'   applies. Default \code{5e-324} (R's smallest positive double).
#' @param dropZeroEffect Logical. Drop rows where any effect column is exactly
#'   0 (\code{BETA}, \code{LOG_ODDS}, \code{SIGNED_SUMSTAT}) or \code{OR}
#'   is exactly 1. Default \code{TRUE}.
#' @param dropNonpositiveSe Logical. Drop rows where \code{SE <= 0}. Default
#'   \code{TRUE}.
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' sumstatsCleaningConfig(clampSmallP = FALSE, dropZeroEffect = FALSE)
#' @export
sumstatsCleaningConfig <- function(
    coerceNumeric = TRUE,
    normalizeChr = TRUE,
    dropNonstandardChr = TRUE,
    dropMissData = TRUE,
    dropPOutOfRange = TRUE,
    clampSmallP = TRUE,
    smallPFloor = 5e-324,
    dropZeroEffect = TRUE,
    dropNonpositiveSe = TRUE
) {
    .newMethodConfig(
        NULL,
        defaults = list(
            coerceNumeric = coerceNumeric,
            normalizeChr = normalizeChr,
            dropNonstandardChr = dropNonstandardChr,
            dropMissData = dropMissData,
            dropPOutOfRange = dropPOutOfRange,
            clampSmallP = clampSmallP,
            smallPFloor = smallPFloor,
            dropZeroEffect = dropZeroEffect,
            dropNonpositiveSe = dropNonpositiveSe
        ),
        extra = list(),
        label = "sumstatsCleaningConfig",
        engine = "sumstatsCleaning"
    )
}

# The cleaning settings as a plain list with every default applied.
# .assertMethodConfig accepts `list()` as "no options", and .applySanityChecks
# tests these with `if (!flag)`, so a missing field must not arrive as NULL.
# @noRd
.sumstatsCleaningResolve <- function(sumstatsCleaningArgs) {
    .assertMethodConfig(
        sumstatsCleaningArgs,
        "sumstatsCleaningConfig",
        "sumstatsCleaning"
    )
    list(
        coerceNumeric = sumstatsCleaningArgs$coerceNumeric %||% TRUE,
        normalizeChr = sumstatsCleaningArgs$normalizeChr %||% TRUE,
        dropNonstandardChr = sumstatsCleaningArgs$dropNonstandardChr %||% TRUE,
        dropMissData = sumstatsCleaningArgs$dropMissData %||% TRUE,
        dropPOutOfRange = sumstatsCleaningArgs$dropPOutOfRange %||% TRUE,
        clampSmallP = sumstatsCleaningArgs$clampSmallP %||% TRUE,
        smallPFloor = sumstatsCleaningArgs$smallPFloor %||% 5e-324,
        dropZeroEffect = sumstatsCleaningArgs$dropZeroEffect %||% TRUE,
        dropNonpositiveSe = sumstatsCleaningArgs$dropNonpositiveSe %||% TRUE
    )
}

# =============================================================================
# Cross-validation settings
# -----------------------------------------------------------------------------
# One bundle for both fineMappingPipeline() and twasWeightsPipeline(). The two
# had drifted on a setting that means the same thing in each: `cvFolds`
# defaulted to 0 in fine-mapping and 5 in TWAS. Folding them into one
# constructor is the rule the joint specification already follows.
#
# `folds` is 0 everywhere. twasWeightsPipeline used to default to 5, so the
# two pipelines disagreed on what an unspecified `cvFolds` meant; converging
# on "off unless asked" makes one bundle mean one thing. The knock-on is that
# the SR-TWAS ensemble, which needs out-of-fold predictions, no longer runs by
# default -- so `ensemble` defaults to FALSE too, and asking for an ensemble
# without folds is an error rather than a silently missing row.
#
# Two of the five fields are honoured by twasWeightsPipeline() only, and are
# IGNORED rather than rejected by fineMappingPipeline() -- one bundle the
# caller can hand to either pipeline without rewriting it:
#
#   weightMethods -- selects among TWAS *weight* methods, matching tokens like
#     `lasso` / `mrmash` and their `<token>_weights` spellings. Fine mapping's
#     `methods=` are fine-mapping methods and its CV refits all of them, so
#     there is no weight-method set to select from. It is also the knob that
#     decides whether .jointTwasCvBlocked() reaches back for a fold's
#     fine-mapping fit -- a TWAS-to-fine-mapping handoff with no counterpart
#     in the other direction.
#   maxVariants -- caps the CV design matrix for twasWeightsCv(). Fine mapping
#     does not use that engine; .fmWeightsCv() has its own loop and only
#     mirrors twasWeightsCv()'s output shape.
#
# Both are documented as TWAS-only so the silence is stated, not discovered.
#
# `seed` is not here either. It seeds the whole call in both pipelines
# (fine-mapping wraps the call in withr::local_seed; TWAS also seeds the
# BiocParallel RNG for method fitting), so it is a call-level knob that CV
# happens to benefit from, not a CV setting.
# =============================================================================

#' @title Cross-Validation Settings
#' @description How cross-validation runs, shared by
#'   \code{\link{fineMappingPipeline}} and \code{\link{twasWeightsPipeline}}.
#'
#'   \code{folds}, \code{numThreads} and \code{samplePartition} mean the same
#'   thing in both. \code{maxVariants} and \code{weightMethods} are honoured
#'   by \code{\link{twasWeightsPipeline}} only and are \strong{ignored} by
#'   \code{\link{fineMappingPipeline}}, so one bundle can be handed to
#'   either pipeline unchanged.
#' @param folds Integer. Number of folds; \code{0} (the default) or
#'   \code{1} skips cross-validation. The same value in both pipelines ---
#'   \code{twasWeightsPipeline} used to default to \code{5}, which made an
#'   unspecified setting mean two different things.
#' @param numThreads Integer. Parallel workers for the per-fold refits.
#'   \code{1} (default) is serial; \code{-1} uses all cores. Only consulted
#'   when \code{folds > 1}.
#' @param samplePartition Optional pre-defined partition \code{data.frame}
#'   with columns \code{Sample} and \code{Fold}. When supplied, every method
#'   reuses this exact partition instead of generating a fresh one.
#' @param maxVariants Integer. Cap on the number of variants used for CV;
#'   unset means no limit. \strong{\code{twasWeightsPipeline} only} ---
#'   \code{fineMappingPipeline} ignores it, because it does not use the
#'   \code{twasWeightsCv} engine this configures.
#' @param weightMethods Optional override of which methods are
#'   cross-validated, as a character vector of tokens or a named method list.
#'   Unset cross-validates every method that produced non-zero weights.
#'   \strong{\code{twasWeightsPipeline} only} ---
#'   \code{fineMappingPipeline} ignores it, because its \code{methods=} are
#'   fine-mapping methods and its CV refits all of them.
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' crossValidationConfig(folds = 10, numThreads = 4)
#' @export
crossValidationConfig <- function(
    folds = 0,
    numThreads = 1,
    samplePartition = NULL,
    maxVariants = NULL,
    weightMethods = NULL
) {
    .newMethodConfig(
        NULL,
        defaults = list(
            folds = folds,
            numThreads = numThreads,
            samplePartition = samplePartition,
            maxVariants = maxVariants,
            weightMethods = weightMethods
        ),
        extra = list(),
        label = "crossValidationConfig",
        engine = "crossValidation"
    )
}

# The fields a pipeline reads, as a plain list so call sites index it without
# caring that the caller passed a constructor result.
# @noRd
.cvResolve <- function(crossValidationArgs) {
    list(
        folds = crossValidationArgs$folds %||% 0,
        numThreads = crossValidationArgs$numThreads %||% 1,
        samplePartition = crossValidationArgs$samplePartition,
        # -1 is twasWeightsCv()'s "no cap" sentinel. fineMappingPipeline
        # resolves these two as well but never reads them.
        maxVariants = crossValidationArgs$maxVariants %||% -1,
        weightMethods = crossValidationArgs$weightMethods
    )
}

# TRUE when cross-validation will actually run. Both 0 and 1 mean "no CV",
# and anything that needs out-of-fold predictions has to ask this rather than
# testing `folds` itself.
# @noRd
.cvEnabled <- function(crossValidationArgs) {
    folds <- crossValidationArgs$folds %||% 0
    !is.null(folds) && length(folds) == 1L && !is.na(folds) && folds >= 2L
}

# Refuse a cross-validation request on an input that has no samples to hold
# out. Summary statistics carry no individual-level data, so CV there is not
# unimplemented but meaningless -- saying so beats running to completion and
# returning results that were never cross-validated.
# @noRd
.cvRefuseOnSumstats <- function(crossValidationArgs, pipeline, cls) {
    set <- names(crossValidationArgs)[
        map_lgl(as.list(crossValidationArgs), .cvFieldIsSet)
    ]
    if (length(set) == 0L) {
        return(invisible(NULL))
    }
    abort(glue(
        "{pipeline}: cross-validation is not possible on {cls} input -- ",
        "it holds out samples, and summary statistics carry none. ",
        "Remove crossValidationConfig({str_flatten(set, ', ')})."
    ))
}

# A field counts as "set" when it differs from the constructor's own default,
# so passing crossValidationConfig() unchanged is not mistaken for a request.
# @noRd
.cvFieldIsSet <- function(x) {
    if (is.null(x)) {
        return(FALSE)
    }
    # numThreads = 1 and folds = 0 are the no-op defaults.
    !(length(x) == 1L && !is.na(x) && is.numeric(x) && (x == 0 || x == 1))
}

# The ensemble settings as a plain list, so a cfg record can carry the group
# as one named field instead of four loose ones.
# @noRd
.ensembleResolve <- function(ensembleArgs) {
    list(
        enabled = isTRUE(ensembleArgs$enabled),
        r2Threshold = ensembleArgs$r2Threshold %||% 0.01,
        solver = ensembleArgs$solver %||% "quadprog",
        alpha = ensembleArgs$alpha %||% 1
    )
}
