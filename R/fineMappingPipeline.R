#' @include qtlSumStats.R gwasSumStats.R MethodParam.R
#' @title Fine-Mapping Pipeline
#' @description S4-dispatched per-region fine-mapping entry point that
#'   replaces the deprecated \code{univariateAnalysisPipeline},
#'   \code{multivariateAnalysisPipeline}, \code{rssAnalysisPipeline},
#'   and \code{susieRssPipeline} pipelines. Accepts:
#'   \itemize{
#'     \item a \code{\link{QtlDataset}} for individual-level cohort
#'           fits (per-context / per-trait univariate SuSiE; joint
#'           multi-trait or multi-context mvSuSiE; joint multi-trait
#'           fSuSiE per context);
#'     \item a \code{\link{MultiStudyQtlDataset}} which recurses through
#'           each embedded \code{QtlDataset} per study and processes
#'           the optional embedded \code{QtlSumStats} via the
#'           sumstat method;
#'     \item a \code{\link{QtlSumStats}} for per-trait SuSiE-RSS fits
#'           and per-(study, trait) multi-context mvSuSiE-RSS fits;
#'     \item a \code{\link{GwasSumStats}} for per-(study, LD-block)
#'           SuSiE-RSS fine-mapping (used by
#'           \code{\link{qtlEnrichmentPipeline}} downstream).
#'   }
#'
#'   Method tokens are unified across input classes; auto-dispatch
#'   picks the individual-level vs RSS implementation based on the
#'   input class. The supported tokens are:
#'   \describe{
#'     \item{\code{susie}}{\code{susieR::susie} with
#'           \code{unmappable_effects = "none"} on individual-level
#'           input; \code{susieR::susie_rss} (same) on RSS.}
#'     \item{\code{susieInf}}{\code{unmappable_effects = "inf"}
#'           variant of the same.}
#'     \item{\code{susieAsh}}{\code{unmappable_effects = "ash"}
#'           variant of the same.}
#'     \item{\code{ser}}{Single-effect regression via
#'           \code{susieR::susie_ser} on summary statistics
#'           (\code{z}, \code{n}); LD-free (no \code{R}, no \code{L}),
#'           so distinct from \code{susie} with \code{L = 1}. Sumstat
#'           input only (\code{QtlSumStats} / \code{GwasSumStats}, or the
#'           sumstat side of a \code{MultiStudyQtlDataset}); rejected on
#'           individual-level \code{QtlDataset}.}
#'     \item{\code{mvsusie}}{\code{mvsusieR::mvsusie} on individual-
#'           level input (requires multi-trait OR multi-context Y),
#'           \code{mvsusieR::mvsusie_rss} on sumstat input (requires
#'           multi-context within a single (study, trait) group).
#'           Errors on \code{GwasSumStats} input.}
#'     \item{\code{fsusie}}{\code{fsusieR::susiF} joint multi-trait fit
#'           per context. Individual-level only; errors on any
#'           SumStats input.}
#'     \item{\code{mrmash}}{Always rejected here. \code{mr.mash} is
#'           a TWAS weight-oriented method and lives in
#'           \code{\link{twasWeightsPipeline}}.}
#'   }
#'
#' @section Chained initialisation: When \code{susieInf} is requested alongside
#'   \code{susie} and/or \code{susieAsh} and
#'   \code{addSusieInf = TRUE} (the default), the SuSiE-inf fit is computed
#'   first and used as initialisation
#'   for the SuSiE / SuSiE-ash fits, mirroring the legacy
#'   \code{univariateAnalysisPipeline} / \code{susieRssPipeline} chained init
#'   behaviour. SuSiE-inf is dropped from the final result when the caller did
#'   not explicitly request it (only used as init).
#'
#' @section QC contract: All \code{QtlSumStats} and \code{GwasSumStats} inputs
#'   must have been QC'd via \code{\link{summaryStatsQc}}; the pipeline errors
#'   on inputs where \code{length(getQcInfo(x)) == 0L}. \code{summaryStatsQc}
#'   also drops variants absent from the \code{ldSketch}, so by the time
#'   per-entry processing runs every variant is guaranteed to be present in the
#'   LD panel and a local LD matrix can be built with
#'   \code{extractBlockGenotypes} + \code{computeLd("sample")}.
#'
#' @section Optional resume cache: Supplying \code{fineMappingResult} of an
#'   existing \code{FineMappingResult} skips re-fitting any \code{(study,
#'   context, trait, method)} tuple that already has a matching row; cached
#'   entries are merged with the newly-fit entries in the returned collection.
#'
#' @section Intentional behaviours dropped from the pre-stub pipelines:
#' The four pre-stub pipelines (\code{univariateAnalysisPipeline} /
#' \code{multivariateAnalysisPipeline} / \code{rssAnalysisPipeline} /
#' \code{susieRssPipeline}) carried several behaviours that are
#' deliberately not ported here:
#' \itemize{
#'   \item TWAS weights computation (\code{twasWeights = TRUE} path):
#'         lives in \code{\link{twasWeightsPipeline}} now.
#'   \item Filtering knobs (now \code{genotypeFilter} /
#'         \code{panelFilter}, \code{ldReferenceMetaFile}): individual-
#'         level QC lives on the \code{QtlDataset} constructor; sumstat
#'         QC lives in \code{summaryStatsQc()}. No filtering happens
#'         inside this pipeline.
#'   \item Diagnostic re-analysis paths
#'         (\code{singleEffect} / \code{bayesianConditionalRegression}
#'         reanalysis on the RSS path): these are not exposed as
#'         dedicated method tokens. Callers who want a single-effect
#'         fit can request it via per-method kwargs, e.g.
#'         \code{methods = FineMappingMethodsParam(susie = list(L = 1))}
#'         \code{methods} parameter).
#'   \item \code{loadRssData} and explicit
#'         \code{ldReferenceMetaFile} arguments: the new
#'         \code{QtlSumStats} / \code{GwasSumStats} carry the
#'         (already-QC'd) sumstats and \code{ldSketch} directly.
#'   \item Verbose \code{methodName} suffixing (e.g.
#'         \code{"susie_rss_NO_QC"}, \code{"susie_rss_SLALOM_RAISS_imputed"}):
#'         the method column on the returned \code{FineMappingResult}
#'         carries the bare token (\code{"susie"},
#'         \code{"susieInf"}, \code{"mvsusie"}, ...) only. QC
#'         provenance is recorded on the sumstats' \code{qcInfo}.
#'   \item An entry that \code{summaryStatsQc(pipCutoffToSkip = ...)} screened
#'         out (recorded as \code{qcInfo$entryAudit[[i]]$pipScreenSkipped}, and
#'         emptied to 0 variants) is \strong{skipped}, not fit: it produces no
#'         row and a message with the screen reason. An all-screened collection
#'         yields a valid empty result rather than an error.
#' }
#'
#' @param data A \code{QtlDataset}, \code{MultiStudyQtlDataset},
#'   \code{QtlSumStats}, or \code{GwasSumStats}.
#' @param methods Method specification. Accepts either:
#'   \itemize{
#'     \item A character vector of method tokens, e.g.
#'           \code{c("susie", "susieInf", "mvsusie")} (any subset of
#'           \code{c("susie", "susieInf", "susieAsh", "mvsusie", "fsusie")},
#'           subject to per-class compatibility).
#'     \item A \code{\link{FineMappingMethodsParam}} record keyed by method
#'           token, where each value is a list of per-method kwargs to
#'           splice into the underlying fitter, e.g.
#'           \code{FineMappingMethodsParam(susie = list(L = 1, refine =
#'                      FALSE), mvsusie = list(max_iter = 500))}. Each
#'           entry is checked against the engine that receives it. Mirrors the
#'           convention of \code{\link{twasWeightsPipeline}}'s
#'           \code{methods} argument. User-supplied kwargs override the
#'           capability-table defaults and any base / chained args set
#'           by the pipeline (e.g. you can override \code{model_init}
#'           even when fitting from a susieInf chain).
#'   }
#' @param contexts Optional character vector of context names. Default
#'   \code{NULL} (all contexts).
#' @param traitId Optional character vector of trait names to restrict
#'   processing to.
#' @param region Optional variant window for QtlDataset trait selection: a
#'   \code{GRanges}, a \code{"chr:start-end"} string, or a one-row data.frame
#'   with \code{chrom}/\code{start}/\code{end}. Mutually exclusive with
#'   \code{traitId}.
#' @param cisWindow For QtlDataset: cis-window (bp) around each trait's genomic
#'   position when extracting variants. Required when \code{traitId} is
#'   supplied. Mutually exclusive with \code{region}.
#' @param jointRegions For QtlDataset with a multi-range \code{region}:
#'   \code{FALSE} (default) fits each range independently and merges the
#'   per-range results into one entry per (study, context, trait, method) -- the
#'   merged \code{susieFit} is a named list of per-region fits and credible-set
#'   labels are renumbered to stay unique. \code{TRUE} concatenates the ranges'
#'   genotypes into one joint fit. Ignored for a single-range / cis
#'   (\code{traitId} + \code{cisWindow}) request.
#' @param addSusieInf Logical. When \code{susieInf} is in \code{methods}
#'   alongside \code{susie} and/or \code{susieAsh}, whether the SuSiE-inf
#'   fit initialises the chained downstream method(s). Default \code{TRUE}.
#'   A choice between methods rather than a property of the credible sets,
#'   which is why it stands alone.
#' @param credibleSetArgs How credible sets are built and reported, built with
#'   \code{\link{CredibleSetParam}}: \code{coverage} (default \code{0.95}),
#'   \code{secondaryCoverage} (\code{c(0.7, 0.5)}), \code{signalCutoff}
#'   (the PIP cutoff for top-loci selection, \code{0.025}),
#'   \code{minAbsCorr} (\code{0.8}) and \code{medianAbsCorr}
#'   (\code{NULL}) for purity --- a set is kept if it passes either, OR-logic
#'   --- and \code{includeAllCs}. \code{perCsColumns} adds the per-set
#'   variant-level columns to \code{topLoci}. \code{L} bounds how many
#'   credible sets can exist (default \code{10}) and \code{Lgreedy} the
#'   greedy-L loop of the SuSiE-inf refinement.
#'
#'   Most fields are read by \code{\link{postprocessFinemappingFits}}, so a
#'   stored fit can be re-summarized under different settings without
#'   refitting. \code{L} / \code{Lgreedy} are the exception: they are seeded
#'   onto each SuSiE-family token's own arguments, so changing them needs a
#'   refit.
#' @param fineMappingResult Optional existing \code{FineMappingResult} to use as
#'   a resume cache; tuples already present are not refit.
#' @param crossValidationArgs Cross-validation settings, built with
#'   \code{\link{CrossValidationParam}}: \code{folds} (default \code{0}, no
#'   CV), \code{numThreads} for the per-fold refits, and \code{samplePartition}
#'   to reuse a fixed partition across methods. When \code{folds > 1} each
#'   method is refit on every fold's training samples and used to predict the
#'   held-out ones; the partition, the out-of-fold predictions and the metrics
#'   are stored on each \code{FineMappingRow}'s \code{cvResult} slot.
#'
#'   Cross-validation holds out \strong{samples}, so it is not possible on
#'   \code{QtlSumStats} / \code{GwasSumStats} input. Those methods accept
#'   the argument and \strong{error} if it is set, rather than accepting it
#'   and quietly returning results that were never cross-validated.
#'
#'   \code{maxVariants} and \code{weightMethods} are honoured by
#'   \code{\link{twasWeightsPipeline}} only and are ignored here.
#' @param seed Optional integer. When non-NULL,
#'   \code{withr::local_seed(seed)} is called at the start of the call for
#'   reproducible fits; the global RNG state is restored on return. Default
#'   \code{NULL}
#'   (no seeding).
#' @param signalScreenArgs Individual-level single-effect (SER) pre-screen
#'   applied to each residualized \code{(X, y)} block before a full fit,
#'   built with \code{\link{SignalScreenParam}}. A susie model with
#'   \code{L = 1} is fit and the block is skipped unless the chosen metric
#'   exceeds its cutoff:
#'   \itemize{
#'     \item \code{pip} --- the maximum PIP. A negative cutoff uses the
#'       adaptive \code{3 / nVariants} threshold.
#'     \item \code{absZ} --- the maximum marginal \code{|z|}, which needs
#'       no fit.
#'     \item \code{bf} / \code{logBf} --- the maximum per-variant Bayes
#'       factor or log Bayes factor from the \code{L = 1} fit.
#'   }
#'   Only one metric may be enabled, which
#'   \code{\link{SignalScreenParam}} enforces when you build it. Unset (the
#'   default) screens nothing. The summary-statistics analog lives in
#'   \code{\link{summaryStatsQc}}.
#' @param usePCA Logical (length 1). \code{QtlDataset} only. When \code{TRUE}
#'   (default \code{FALSE}), each multi-trait context's PCA-reduced phenotype is
#'   fine-mapped with univariate SuSiE on its top principal components (ports
#'   the legacy \code{fsusie.R} \code{susie_on_top_pc}). Each PC becomes a
#'   pseudo-trait row keyed \code{trait = "topPC\{i\}"}, \code{method =
#'   "susie"}. Single-trait contexts have no PCA and are skipped.
#' @param nPCs Integer (length 1). \code{QtlDataset} only. Caps the number of
#'   top principal components fine-mapped per context when \code{usePCA = TRUE}
#'   (default \code{10}). The effective count is \code{min(nPCs, usable
#'   traits)}.
#' @param jointSpecification Optional joint-fit specification (NULL by default).
#'   When NULL, the pipeline runs the implicit multi-context / multi-trait
#'   mvSuSiE / fSuSiE branches as before. When non-NULL, the argument is parsed
#'   and validated via the joint-spec grammar documented under
#'   \code{parseJointSpecification} (a character vector of axes, or a list of
#'   \code{list(axes, scope)} specs); the per-spec axis dispatcher
#'   implementation is in progress and a non-NULL value currently errors with an
#'   informative message. See the design notes in \code{R/jointSpecification.R}
#'   for the accepted grammar.
#' @section Panel filters on the RSS path: On \code{QtlSumStats} /
#'   \code{GwasSumStats} input there is no genotype matrix to filter, so
#'   \code{panelFilter}'s \code{mafCutoff} / \code{macCutoff} /
#'   \code{imissCutoff} are measured against the \strong{LD reference panel}
#'   instead: a variant whose panel genotypes fall below the cutoffs is
#'   dropped before the z-scores and LD matrix are built. The thresholds mean
#'   the same thing as \code{genotypeFilter}'s on the \code{QtlDataset} path
#'   (MAC is converted to a MAF equivalent and the stricter of the two
#'   applies), so one number carries across input types. The defaults filter
#'   nothing.
#'
#'   This discards \emph{observed} variants, unlike
#'   \code{summaryStatsQc(imputeArgs = ...)}, which only bounds which variants
#'   RAISS will impute. Use it when the panel cannot support the LD estimate a
#'   rare variant would need.
#'
#' @param genotypeFilterArgs Per-call overrides of the \code{QtlDataset}'s own
#'   filters, built with \code{\link{GenotypeFilterParam}}. A field left unset
#'   keeps the dataset's construct-time value, so
#'   \code{GenotypeFilterParam(mafCutoff = 0)} pins the cutoff at zero while
#'   \code{GenotypeFilterParam()} changes nothing. Applies to the
#'   \code{QtlDataset} method only.
#' @param panelFilterArgs LD-reference-panel filters for the summary-statistics
#'   and GWAS methods, built with \code{\link{PanelFilterParam}}. See
#'   \emph{Panel filters on the RSS path} above.
#' @param mrmashPrior Optional \code{\link{TwasWeights}} from a previous
#'   mr.mash \code{\link{twasWeightsPipeline}} run, mined for the mvSuSiE
#'   data-driven prior: the retained \code{mrmash} fits supply the prior
#'   mixture, and their cross-validation payload supplies honest per-fold
#'   priors. \code{NULL} uses the canonical prior instead.
#' @param dataDrivenPriorWeightsCutoff Numeric or \code{NULL}. Cutoff below
#'   which data-driven prior weights are zeroed; \code{NULL} disables the
#'   cutoff.
#' @param naAction Character. How to handle missing values in the extracted
#'   phenotype/genotype data.
#' @param verbose Verbosity (0 silent, 1 default). Default \code{1}.
#' @param residualizationArgs Covariate residualization settings, built with
#'   \code{\link{ResidualizationParam}}: \code{phenotypeCovariates} and
#'   \code{genotypeCovariates} name which covariates to regress out
#'   (\code{NULL}, the default, uses every available one), and
#'   \code{residualizePhenotype} / \code{residualizeGenotype} turn each side
#'   off.
#'
#'   Only meaningful for \code{QtlDataset} / \code{MultiStudyQtlDataset}
#'   input. A \code{QtlSumStats} / \code{GwasSumStats} run has no covariates
#'   to regress out, so the bundle is \strong{ignored} there rather than
#'   refused --- unlike \code{crossValidation}, which errors. Residualization
#'   is on by default, so refusing a non-default value would reject the
#'   default bundle; CV is off by default, so a non-default value there is an
#'   explicit request for something the input cannot do.
#'
#'   The marginal univariate effects stored on each \code{FineMappingRow}
#'   obey the same choice as the SuSiE fit itself --- they are computed
#'   against the same residualized \code{X} / \code{Y}.
#' @param fitRetention How much of each fit the \code{susieFit} slot keeps:
#'   \code{"slim"} (default) a trimmed view, or \code{"full"} the whole
#'   \code{susie()} return, so \code{getSusieFit()} and
#'   non-default-coverage \code{getCs()} queries can read the full posterior
#'   matrices. The per-variant \code{topLoci} table is fully populated
#'   either way --- its per-credible-set columns are governed by
#'   \code{credibleSet}'s \code{perCsColumns}, not by this.
#' @param rssArgs Summary-statistics solver settings for \code{QtlSumStats} /
#'   \code{GwasSumStats} input, built with \code{\link{SusieRssParam}}:
#'   \itemize{
#'     \item \code{serFallback} (default \code{FALSE}) --- after each
#'       multi-effect SuSiE-RSS fit, read susieR's finite-sample R
#'       diagnostics and, when
#'       \code{fit$R_finite_diagnostics$R_reliability_flag} is \code{TRUE},
#'       report the single-effect (\code{ser_model}) result for that region
#'       instead.
#'     \item \code{keepFullFit} (default \code{"fallback"}) --- which
#'       pre-fallback multi-effect fits to retain: only the regions that fell
#'       back, \code{"all"}, or \code{"none"}. Meaningful only with
#'       \code{serFallback}. The retained fit and the decision are read via
#'       \code{getSusieFit(res)$multiEffectFit},
#'       \code{$R_reliability_flag} and \code{$serFallbackUsed}.
#'     \item \code{rFinite} --- finite-sample size for susieR's
#'       \code{R_finite} correction. \code{NULL} (default) uses susieR's
#'       own, except when a finite/EB mode is active (\code{serFallback} or
#'       \code{rMismatch != "none"}), where it falls back to the LD panel's
#'       \code{getNSamples(ldSketch)}.
#'     \item \code{rMismatch} (default \code{"none"}) --- LD-mismatch
#'       correction forwarded as \code{R_mismatch}: also \code{"eb"} or
#'       \code{"eb_mix"}.
#'     \item \code{control} --- options for
#'       \code{susieR::susie_rss_control()}, built with
#'       \code{\link{SusieRssControlOptions}}. This was a bare named list,
#'       which meant a misspelled option was dropped silently; the
#'       constructor checks the names against that function's live formals.
#'   }
#' @param ... Reserved for future per-method arguments.
#'
#' @return A \code{FineMappingResult} collection keyed by \code{(study, context,
#'   trait, method)}. The \code{ldSketch} slot is set automatically: \code{NULL}
#'   for individual-level (QtlDataset / all-individual-level
#'   MultiStudyQtlDataset) fits, the input's \code{ldSketch} for RSS-derived
#'   fits.
#' @examples
#' data(qtlDatasetExample)
#' fineMappingPipeline(qtlDatasetExample, methods = "susie", cisWindow = 1e6)
#' @export
setGeneric("fineMappingPipeline", function(data, ...) {
    # Checked here rather than in an ANY method: the accepted classes are
    # known, so this is a validity check, and running it before dispatch
    # means the message survives whatever else the call passed.
    .fmAssertInputClass(data)
    standardGeneric("fineMappingPipeline")
})

# @noRd
.fmAssertInputClass <- function(data) {
    ok <- c("QtlDataset", "MultiStudyQtlDataset", "QtlSumStats", "GwasSumStats")
    if (any(map_lgl(ok, is, object = data))) {
        return(invisible(NULL))
    }
    abort(glue(
        "fineMappingPipeline does not accept inputs of class ",
        "'{class(data)[[1L]]}'. Pass a QtlDataset, MultiStudyQtlDataset, ",
        "QtlSumStats, or GwasSumStats. Use summaryStatsQc() on SumStats ",
        "inputs first."
    ))
}

# =============================================================================
# Method capability table -- unified naming, individual vs sumstat dispatch
# =============================================================================

# `individualImpl`  : function-call symbol used when input is QtlDataset /
#                     MultiStudyQtlDataset (NULL = not supported).
# `sumstatImpl`     : function-call symbol used when input is QtlSumStats /
#                     GwasSumStats (NULL = not supported).
# `multivariate`    : requires a multi-trait or multi-context joint Y
#                     (mvsusie / mvsusie_rss / fsusie).
# `gwasAllowed`     : whether the method is permitted on a GwasSumStats
#                     input. Only the SuSiE-RSS family supports per-LD-block
#                     GWAS fine-mapping.
# `unmappableEffects`: the value passed to susieR::susie /
#                     susieR::susie_rss to switch between susie / susieInf /
#                     susieAsh variants. NA for non-SuSiE-family methods.
#
# This table lists ONLY fine-mapping methods. TWAS-weight-oriented tokens (e.g.
# mr.mash) are not here -- they live in .fmTwasOnlyTokens and are rejected with
# a clear pointer to twasWeightsPipeline() (see .fmCheckMethodCapabilities).
#
# @noRd
.fineMappingMethodCapabilities <- list(
    susie = list(
        individualImpl = "susieR::susie",
        sumstatImpl = "susieR::susie_rss",
        multivariate = FALSE,
        gwasAllowed = TRUE,
        unmappableEffects = "none",
        args = list()
    ),
    susieInf = list(
        individualImpl = "susieR::susie",
        sumstatImpl = "susieR::susie_rss",
        multivariate = FALSE,
        gwasAllowed = TRUE,
        unmappableEffects = "inf",
        args = list()
    ),
    susieAsh = list(
        individualImpl = "susieR::susie",
        sumstatImpl = "susieR::susie_rss",
        multivariate = FALSE,
        gwasAllowed = TRUE,
        unmappableEffects = "ash",
        args = list()
    ),
    # Single-effect regression (SER) on summary statistics via
    # susieR::susie_ser.
    # LD-free (z + n; no R, no L), so distinct from susie with L = 1.
    # Sumstat-only (individualImpl = NULL): runs on QtlSumStats / GwasSumStats
    # (and the sumstat side of a MultiStudyQtlDataset); rejected on
    # individual-level QtlDataset.
    ser = list(
        individualImpl = NULL,
        sumstatImpl = "susieR::susie_ser",
        multivariate = FALSE,
        gwasAllowed = TRUE,
        unmappableEffects = NA_character_,
        args = list()
    ),
    mvsusie = list(
        individualImpl = "mvsusieR::mvsusie",
        sumstatImpl = "mvsusieR::mvsusie_rss",
        multivariate = TRUE,
        gwasAllowed = FALSE,
        unmappableEffects = NA_character_,
        args = list()
    ),
    fsusie = list(
        individualImpl = "fsusieR::susiF",
        sumstatImpl = NULL,
        multivariate = TRUE,
        gwasAllowed = FALSE,
        unmappableEffects = NA_character_,
        args = list()
    )
)

# TWAS-weight-oriented method tokens. NOT fine-mapping methods (they belong to
# twasWeightsPipeline); enumerated only so .fmCheckMethodCapabilities rejects
# them with a clear pointer rather than an "unknown token" error.
.fmTwasOnlyTokens <- c("mrmash")

# Normalize a user-supplied `methods` argument into a character vector of
# canonical tokens. Mirrors `.twasNormalizeMethods` but the fine-mapping
# pipeline takes only a character vector (no preset strings, no list form).
# @noRd
# Normalize a user-supplied `methods` argument into `(tokens, methodArgs)`.
#
# Accepts: * character vector c("susie", "susieInf") -> empty kwargs per token
# * a FineMappingMethodsParam() record -> per-token kwargs, already checked
#
# Names of the returned `methodArgs` always equal `tokens` (one entry per
# token, empty list when the user supplied none). The fitters then
# `modifyList`-merge each entry into the base arg list before do.call.
#
# Mirrors the convention of .twasNormalizeMethods so the two pipelines
# expose the same shape on the user side.
#
# A plain named list is routable HERE even though the constructor refuses one
# under `methods`: this call knows its input class, so the overrides go into
# that path's slot and are checked against the engine that will actually
# receive them. A retired FineMappingMethodsParam() record is still accepted
# too.
# @noRd
.fmNormalizeMethods <- function(
    methods,
    inputKind,
    L = 10L,
    Lgreedy = NULL
) {
    if (is.null(methods) || length(methods) == 0L) {
        msg <- glue(
            "fineMappingPipeline: `methods` must be a non-empty character ",
            "vector, a named list of per-method options, or a ",
            "FineMappingMethodsParam() record."
        )
        abort(msg)
    }
    parsed <- if (is.character(methods)) {
        .fmMethodsFromTokens(methods)
    } else if (is(methods, "MethodsSelectionParam") || is.list(methods)) {
        .methodsParamResolve(
            .methodsParamFor(
                methods,
                inputKind,
                "FineMappingMethodsParam",
                "fineMappingPipeline"
            ),
            inputKind
        )
    } else if (.isMethodOptions(methods)) {
        list(
            tokens = names(methods),
            methodArgs = map(as.list(methods), as.list)
        )
    } else {
        cls <- class(methods)[[1L]]
        msg <- glue(
            "fineMappingPipeline: `methods` must be a character vector, a ",
            "named list of per-method options, or a ",
            "FineMappingMethodsParam() record. Got class '{cls}'."
        )
        abort(msg)
    }
    list(
        tokens = parsed$tokens,
        methodArgs = .fmSeedSusieDefaults(
            parsed$methodArgs,
            parsed$tokens,
            L,
            Lgreedy
        )
    )
}

# A bare character vector of tokens: every token gets empty kwargs.
# @noRd
.fmMethodsFromTokens <- function(methods) {
    tokens <- unique(methods)
    list(
        tokens = tokens,
        methodArgs = set_names(rep(list(list()), length(tokens)), tokens)
    )
}

# SuSiE-family fit defaults live here (the single source of truth), not in CLI
# wrappers: seed L / L_greedy on every susie-family token whose kwargs did not
# already set them.
# @noRd
.fmSeedSusieDefaults <- function(methodArgs, tokens, L, Lgreedy) {
    seeded <- intersect(tokens, c("susie", "susieInf", "susieAsh"))
    if (length(seeded) == 0L) {
        return(methodArgs)
    }
    list_assign(
        methodArgs,
        !!!set_names(
            map(
                seeded,
                .fmSeedTokenDefaults,
                methodArgs = methodArgs,
                L = L,
                Lgreedy = Lgreedy
            ),
            seeded
        )
    )
}

# One token's kwargs with L / L_greedy filled in where the caller left them
# unset. A token with no kwargs at all still gets the defaults.
# @noRd
.fmSeedTokenDefaults <- function(tk, methodArgs, L, Lgreedy) {
    args <- methodArgs[[tk]] %||% list()
    list_assign(
        args,
        L = args[["L"]] %||% L,
        L_greedy = args[["L_greedy"]] %||% Lgreedy
    )
}

# Which engine a token actually reaches for a given input class. A method may
# have one entry point for individual-level data and another for summary
# statistics, and any single run uses exactly one of them -- so once the input
# class is known the argument names can be checked against that one engine
# rather than against the union of both.
# @noRd
.fmMethodCalleeFor <- function(token, inputKind) {
    info <- .fineMappingMethodCapabilities[[token]]
    if (is.null(info)) {
        return(NULL)
    }
    if (identical(inputKind, "QtlDataset")) {
        return(info$individualImpl)
    }
    info$sumstatImpl
}

# Check each token's user arguments against the single engine this run will
# reach. Stricter than the constructor's check, which cannot know the input
# class: an argument that only `susie_rss` accepts is caught here on an
# individual-level run, and mvsusie becomes checkable on the QtlDataset path
# even though mvsusie_rss takes `...`.
# @noRd
.fmCheckMethodArgsForInput <- function(methodArgs, inputKind) {
    for (token in names(methodArgs)) {
        callee <- .fmMethodCalleeFor(token, inputKind)
        if (is.null(callee)) {
            next
        }
        accepted <- .engineAcceptedNames(callee)
        # Defaults pecotmr seeds itself are not the caller's to answer for.
        given <- setdiff(names(methodArgs[[token]]), .fmSeededArgNames())
        .engineCheckExtra(
            methodArgs[[token]][given],
            accepted,
            glue("fineMappingPipeline: method '{token}'"),
            callee
        )
    }
    invisible(NULL)
}

# Argument names pecotmr seeds on the caller's behalf, which therefore need no
# checking against the engine.
# @noRd
.fmSeededArgNames <- function() {
    c("L", "L_greedy")
}

# Enforce input-class / method compatibility against the fine-mapping
# capability table. Rejects TWAS-weight-oriented tokens (.fmTwasOnlyTokens,
# e.g. mr.mash) with a clear pointer to twasWeightsPipeline(). Routes the input
# class through individual / sumstat / GWAS branches and emits a single error
# listing every offending token.
# @noRd
.fmCheckMethodCapabilities <- function(tokens, inputKind) {
    if (length(tokens) == 0L) {
        return(invisible(NULL))
    }
    caps <- .fineMappingMethodCapabilities
    unknown <- setdiff(tokens, c(names(caps), .fmTwasOnlyTokens))
    if (length(unknown) > 0L) {
        unknownStr <- str_flatten(unknown, ", ")
        knownStr <- str_flatten(names(caps), ", ")
        msg <- glue(
            "fineMappingPipeline: unknown method token(s): {unknownStr}. ",
            "Known tokens: {knownStr}."
        )
        abort(msg)
    }
    issues <- compact(map(
        tokens,
        .fmTokenCapabilityIssue,
        inputKind = inputKind,
        caps = caps
    ))
    if (length(issues) == 0L) {
        return(invisible(NULL))
    }
    bad <- map_chr(issues, "token")
    detail <- map_chr(issues, .fmIssueDetail)
    badStr <- str_flatten(unique(bad), ", ")
    detailStr <- str_flatten(detail, "; ")
    msg <- glue(
        "fineMappingPipeline: the following method(s) are not ",
        "available for input class '{inputKind}': {badStr}. {detailStr}."
    )
    abort(msg)
}

# Capability issue for one method token under `inputKind`: NULL when the token
# is usable, else list(token, reason) describing why it is not available. A
# TWAS-weight token is always rejected (use twasWeightsPipeline instead).
# @noRd
.fmTokenCapabilityIssue <- function(tk, inputKind, caps) {
    if (is_in(tk, .fmTwasOnlyTokens)) {
        return(list(
            token = tk,
            reason = str_c(
                "is a TWAS-weight-oriented method; ",
                "use twasWeightsPipeline()"
            )
        ))
    }
    info <- caps[[tk]]
    reason <- switch(
        inputKind,
        QtlDataset = if (is.null(info$individualImpl)) {
            "is sumstat-only (use a QtlSumStats input)"
        },
        MultiStudyQtlDataset = if (
            is.null(info$individualImpl) && is.null(info$sumstatImpl)
        ) {
            "has no individual or sumstat implementation"
        },
        QtlSumStats = if (is.null(info$sumstatImpl)) {
            "is individual-only (use a QtlDataset input)"
        },
        GwasSumStats = if (
            !isTRUE(info$gwasAllowed) || is.null(info$sumstatImpl)
        ) {
            str_c(
                "is not supported on GwasSumStats (only the SuSiE-RSS ",
                "family is)"
            )
        },
        NULL
    )
    if (is.null(reason)) {
        NULL
    } else {
        list(token = tk, reason = reason)
    }
}

# TRUE if method token `tk` is unknown (kept; validated elsewhere) or its
# capability advertises a non-NULL `capField`.
# @noRd
.fmMethodOk <- function(tk, capField, caps) {
    info <- caps[[tk]]
    is.null(info) || !is.null(info[[capField]])
}

# Keep only the tokens in `methods` whose capability has a non-NULL `capField`
# (individualImpl / sumstatImpl), so a sumstat-only method (e.g. ser) is dropped
# from the individual-level recursion and an individual-only method from the
# sumstat recursion. `methods` is a character vector of tokens or a named list
# of per-token args; unknown tokens pass through (handled elsewhere).
.fmFilterMethodsForKind <- function(methods, capField) {
    caps <- .fineMappingMethodCapabilities
    if (is.character(methods)) {
        methods[map_lgl(methods, .fmMethodOk, capField, caps)]
    } else if (is.list(methods)) {
        methods[map_lgl(names(methods), .fmMethodOk, capField, caps)]
    } else {
        methods
    }
}

# Reject SumStats inputs that have not been QC'd via summaryStatsQc.
# @noRd
.fmAssertQcd <- function(sumstats) {
    if (length(getQcInfo(sumstats)) == 0L) {
        cls <- class(sumstats)[[1L]]
        msg <- glue(
            "fineMappingPipeline: the supplied {cls} has no QC record ",
            "(qcInfo is empty). Call summaryStatsQc() first and pass the ",
            "QC-applied result."
        )
        abort(msg)
    }
}

# Given a `methods` vector, decide whether the SuSiE-inf chained-init
# shortcut applies. Returns a list of (chainSusie, chainAsh, runInf,
# keepInf): runInf is TRUE when susieInf must be fitted (either user
# requested it OR a chained init needs it); keepInf is TRUE when the
# user asked for "susieInf" in `methods` directly.
# @noRd
.fmResolveSusieChain <- function(tokens, addSusieInf) {
    hasInf <- is_in("susieInf", tokens)
    hasSu <- is_in("susie", tokens)
    hasAsh <- is_in("susieAsh", tokens)
    chainSusie <- isTRUE(addSusieInf) && hasInf && hasSu
    chainAsh <- isTRUE(addSusieInf) && hasInf && hasAsh
    runInf <- hasInf || chainSusie || chainAsh
    keepInf <- hasInf
    list(
        chainSusie = chainSusie,
        chainAsh = chainAsh,
        runInf = runInf,
        keepInf = keepInf
    )
}

# Optional resume-cache lookup. Returns the matching FineMappingRow from
# `fineMappingResult` for the tuple (study, context, trait, method), or
# NULL when there is no hit. Returns NULL silently when fineMappingResult
# is NULL or not a QtlFineMappingResult.
# @noRd
.fmCacheLookup <- function(fineMappingResult, study, context, trait, method) {
    if (is.null(fineMappingResult)) {
        return(NULL)
    }
    if (!is(fineMappingResult, "QtlFineMappingResult")) {
        return(NULL)
    }
    idx <- .matchTupleRows(
        fineMappingResult,
        list(study = study, context = context, trait = trait, method = method)
    )
    if (length(idx) == 0L) {
        return(NULL)
    }
    .fmrRowParts(fineMappingResult, idx[[1L]])
}

# GwasFineMappingResult cache lookup using the (study, method, range) identity.
# Multi-block FMRs carry one entry per block, so the key has to include the
# block -- but the block is now the element's own RANGE rather than a stored
# label, which means the cache cannot miss because a label was absent or
# spelled differently.
# @noRd
.fmCacheLookupGwas <- function(fineMappingResult, study, method, blockId) {
    if (is.null(fineMappingResult)) {
        return(NULL)
    }
    if (!is(fineMappingResult, "GwasFineMappingResult")) {
        return(NULL)
    }
    matched <- .matchTupleRows(
        fineMappingResult,
        list(study = study, method = method)
    )
    if (length(matched) == 0L) {
        return(NULL)
    }
    # Several rows share (study, method) when the study was split into blocks;
    # the block key picks the one this lookup wants.
    idx <- if (length(matched) == 1L) {
        matched
    } else {
        matched[.rtlRangeKeys(fineMappingResult)[matched] == blockId]
    }
    if (length(idx) == 0L) {
        return(NULL)
    }
    .fmrRowParts(fineMappingResult, idx[[1L]])
}

# Build a QtlFineMappingResult collection from per-tuple parallel vectors.
# `jointStudies`, `jointContexts`, `jointTraits` are optional character
# vectors (length matches `studies`) describing semicolon-joined joint
# members for cross-study / cross-context / cross-trait joint fits; pass
# `NULL` (default) to omit the column entirely.
# @noRd
.fmBuildQtlResult <- function(
    studies,
    contexts,
    traits,
    methods,
    entries,
    jointStudies = NULL,
    jointContexts = NULL,
    jointTraits = NULL,
    traitPos = NULL,
    ldSketch = NULL,
    allowEmpty = FALSE
) {
    if (length(entries) == 0L && !allowEmpty) {
        msg <- glue(
            "fineMappingPipeline: no (study, context, trait, method) tuples ",
            "produced a fine-mapping result."
        )
        abort(msg)
    }
    QtlFineMappingResult(
        study = studies,
        context = contexts,
        trait = traits,
        method = methods,
        entry = entries,
        jointStudies = jointStudies,
        jointContexts = jointContexts,
        jointTraits = jointTraits,
        traitPos = traitPos,
        ldSketch = ldSketch
    )
}

# Build a GwasFineMappingResult collection from per-row vectors. `blockIds` is
# optional provenance keying the external LD block manifest; row identity comes
# from (study, method) plus the element's own range either way.
# @noRd
.fmBuildGwasResult <- function(
    studies,
    methods,
    entries,
    blockIds = NULL,
    ldSketch = NULL,
    allowEmpty = FALSE
) {
    if (length(entries) == 0L && !allowEmpty) {
        msg <- glue(
            "fineMappingPipeline: no (study, method) tuples produced a ",
            "fine-mapping result."
        )
        abort(msg)
    }
    GwasFineMappingResult(
        study = studies,
        method = methods,
        blockId = blockIds,
        entry = entries,
        ldSketch = ldSketch
    )
}

# One QTL-side result row as an immutable record (study/context/trait/method
# + the FineMappingRow). Dispatch helpers RETURN these; the orchestrator
# flattens them and extracts the parallel vectors -- no mutable accumulator.
# @noRd
.fmQtlRow <- function(study, context, trait, method, entry) {
    list(
        study = study,
        context = context,
        trait = trait,
        method = method,
        entry = entry
    )
}

# One GWAS-side result row as an immutable record (study/method/blockId + the
# FineMappingRow). blockId is provenance for the external block manifest;
# row identity comes from (study, method) plus the element's own range.
# @noRd
.fmGwasRow <- function(study, method, blockId, entry) {
    list(study = study, method = method, blockId = blockId, entry = entry)
}

# Effect-allele frequency vector aligned to `variantIds` from an entry's
# DIRECTIONAL AF mcol (post-QC harmonized/complemented to the final effect
# allele). NULL when the entry carries no declared af -- a directionless MAF
# is NOT used here (it is QC-only), so an undeclared-af entry exports af = NA
# rather than a silently mislabelled minor-allele frequency.
# @noRd
.fmAfByVar <- function(entry, variantIds) {
    mc <- S4Vectors::mcols(entry)
    if (!is_in("AF", colnames(mc))) {
        return(NULL)
    }
    set_names(as.numeric(mc$AF), as.character(mc$SNP))[variantIds]
}

# Block label derived from a GwasSumStats entry's GRanges. Built through the
# same helper the collection uses for its range key, so a label minted here and
# the identity of the row it ends up on agree by construction rather than by
# two formatters happening to match.
# @noRd
.fmGwasBlockId <- function(gr) {
    .rtlOneRangeKey(range(gr))
}

# GWAS resume lookup using the GwasFineMappingResult (study, method, range)
# identity; NULL when no compatible cache was supplied.
# @noRd
.fmCacheLookupGwasResume <- function(fineMappingResult, st, tk, blockId) {
    if (
        !is.null(fineMappingResult) &&
            is(fineMappingResult, "GwasFineMappingResult")
    ) {
        .fmCacheLookupGwas(fineMappingResult, st, tk, blockId)
    } else {
        NULL
    }
}

# Fit the still-to-run RSS tokens for one GWAS region and return one row-record
# per fitted token.
# @noRd
.fmGwasFitRows <- function(gr, zn, st, blockId, toRun, cfg) {
    ents <- .fmFitRssBlockP(
        set_names(zn$z, zn$variantIds),
        zn,
        toRun,
        glue("GWAS (study='{st}', region='{blockId}')"),
        .fmAfByVar(gr, zn$variantIds),
        cfg
    )
    map(names(ents), .fmGwasRowFor, st = st, blockId = blockId, ents = ents)
}

# All result rows for one GwasSumStats entry: cache hits + freshly-fitted
# tokens, or an empty set when the region was screened out. Returns
# list(rows, skipped) so the caller sums the skip flags functionally.
# @noRd
.fmGwasEntryRows <- function(i, studyCol, tokens, cfg) {
    st <- studyCol[[i]]
    gr <- .collectionEntry(cfg$data, i)
    skip <- .fmEntrySkipInfo(cfg$data, i)
    if (isTRUE(skip$skipped)) {
        if (cfg$verbose >= 1) {
            reason <- skip$reason
            inform(glue(
                "fineMappingPipeline(GwasSumStats): study='{st}' region ",
                "skipped: {reason}"
            ))
        }
        return(list(rows = list(), skipped = TRUE))
    }
    zn <- .fmExtractZn(
        gr,
        glue("fineMappingPipeline(GwasSumStats): study='{st}'"),
        ldSketch = cfg$ldSketch,
        cutoffs = .panelCutoffs(cfg$panelFilterArgs)
    )
    blockId <- .fmGwasBlockId(gr)
    lookups <- map(
        tokens,
        .fmGwasLookup,
        fineMappingResult = cfg$fineMappingResult,
        st = st,
        blockId = blockId
    )
    cachedRows <- map(
        keep(lookups, .fmHasCached),
        .fmGwasRowFromLookup,
        st = st,
        blockId = blockId
    )
    toRun <- map_chr(keep(lookups, .fmNotCached), "tk")
    if (length(toRun) == 0L) {
        return(list(rows = cachedRows, skipped = FALSE))
    }
    computed <- .fmGwasFitRows(gr, zn, st, blockId, toRun, cfg)
    list(rows = c(cachedRows, computed), skipped = FALSE)
}

# The fine-mapping settings shared by every summary-statistics path: the
# QtlSumStats joint phase and the RSS / GWAS row builders all act on these.
# Named once here so the records built from it cannot disagree about the same
# setting. Explicitly constructed, not captured from the calling frame.
# @noRd
.fmSsCommonCfg <- function(
    data,
    credibleSetArgs,
    fineMappingResult,
    fitRetention,
    verbose,
    panelFilterArgs
) {
    list(
        data = data,
        credibleSetArgs = credibleSetArgs,
        fineMappingResult = fineMappingResult,
        fitRetention = fitRetention,
        verbose = verbose,
        panelFilterArgs = panelFilterArgs
    )
}

# The per-run settings the RSS / GWAS row builders all need: the common
# summary-statistics record plus the LD and RSS-control fields only the row
# builders use. The chain threads one record instead of twenty arguments, and
# every field is named deliberately -- this is not a capture of the calling
# frame.
# @noRd
.fmRssRunCfg <- function(
    common,
    ldSketch,
    addSusieInf,
    methodArgs,
    rssArgs
) {
    c(
        common,
        list(
            ldSketch = ldSketch,
            addSusieInf = addSusieInf,
            methodArgs = methodArgs,
            rssArgs = rssArgs
        )
    )
}

# .fmFitRssBlock with the per-run knobs supplied from `cfg`; callers pass only
# the block-specific arguments. The RSS analog of .fmFitXBlockP.
# @noRd
.fmFitRssBlockP <- function(z, zn, toRun, label, af, cfg) {
    ldMat <- .ldFromSketch(
        cfg$ldSketch,
        zn$variantIds,
        label = "fineMappingPipeline"
    )
    .fmFitRssBlock(
        z,
        ldMat,
        zn$n,
        toRun,
        cfg$addSusieInf,
        cfg$methodArgs,
        cfg$verbose,
        label = label,
        af = af,
        nVar = zn$nVar,
        credibleSetArgs = cfg$credibleSetArgs,
        fitRetention = cfg$fitRetention,
        rssArgs = cfg$rssArgs
    )
}

# Fit the still-to-run RSS tokens for one QtlSumStats entry and return one
# row-record per fitted token.
# @noRd
.fmRssFitRows <- function(i, st, ctx, tr, toRun, cfg) {
    entry <- .collectionEntry(cfg$data, i)
    zn <- .fmExtractZn(
        entry,
        glue(
            "fineMappingPipeline(QtlSumStats): entry {i} (study='{st}', ",
            "context='{ctx}', trait='{tr}')"
        ),
        ldSketch = cfg$ldSketch,
        cutoffs = .panelCutoffs(cfg$panelFilterArgs)
    )
    ents <- .fmFitRssBlockP(
        set_names(zn$z, zn$variantIds),
        zn,
        toRun,
        glue("(study='{st}', context='{ctx}', trait='{tr}')"),
        .fmAfByVar(entry, zn$variantIds),
        cfg
    )
    map(names(ents), .fmQtlRowFor, st = st, ctx = ctx, tr = tr, ents = ents)
}

# All result rows for one QtlSumStats entry: cache hits first, then (only when
# tokens remain to run and the entry was not screened out) freshly-fitted
# tokens. Returns list(rows, skipped) so the caller sums the skip flags.
# @noRd
.fmRssEntryRows <- function(
    i,
    studyCol,
    contextCol,
    traitCol,
    univTokens,
    cfg
) {
    st <- studyCol[i]
    ctx <- contextCol[i]
    tr <- traitCol[i]
    lookups <- map(
        univTokens,
        .fmQtlLookup,
        fineMappingResult = cfg$fineMappingResult,
        st = st,
        ctx = ctx,
        tr = tr
    )
    cachedRows <- map(
        keep(lookups, .fmHasCached),
        .fmQtlRowFromLookup,
        st = st,
        ctx = ctx,
        tr = tr
    )
    toRun <- map_chr(keep(lookups, .fmNotCached), "tk")
    if (length(toRun) == 0L) {
        return(list(rows = cachedRows, skipped = FALSE))
    }
    skip <- .fmEntrySkipInfo(cfg$data, i)
    if (isTRUE(skip$skipped)) {
        if (cfg$verbose >= 1) {
            reason <- skip$reason
            inform(glue(
                "fineMappingPipeline(QtlSumStats): entry {i} ",
                "(study='{st}', context='{ctx}', trait='{tr}') ",
                "skipped: {reason}"
            ))
        }
        return(list(rows = cachedRows, skipped = TRUE))
    }
    computed <- .fmRssFitRows(i, st, ctx, tr, toRun, cfg)
    list(rows = c(cachedRows, computed), skipped = FALSE)
}

# Concatenate two same-class FineMappingResult collections row-wise, carrying
# forward every column (delegates to the generic `.rbindCollections`).
# @noRd
#' @importFrom checkmate assertClass
.rbindFineMappingResult <- function(a, b, ldSketch = NULL) {
    assertClass(a, "FineMappingResultBase")
    assertClass(b, "FineMappingResultBase")
    # Carry forward every column (blockId / joint* / ...) and reconcile the
    # collection-level slots via the shared combine; the concrete class
    # (QTL vs GWAS) is preserved and checked there.
    .combineTupleCollections(list(a, b), ldSketch, ".rbindFineMappingResult")
}

#' Combine FineMappingResult collections
#'
#' Row-bind two or more fine-mapping result collections of the SAME concrete
#' class (all \code{\link{QtlFineMappingResult}} or all
#' \code{\link{GwasFineMappingResult}}) into one -- e.g. per-block GWAS results
#' into a genome-wide collection for cTWAS. Mixing the two concrete classes is
#' an error.
#'
#' @param ... Two or more \code{FineMappingResultBase} objects, or a single
#'   \code{list} of them.
#' @param ldSketch Optional genotype panel (see
#'   \code{\link{readGenotypes}}) to attach to the
#'   combined collection. Default \code{NULL}. Applied when combining two or
#'   more inputs; a single input is returned unchanged.
#' @return A single combined fine-mapping result of the shared concrete class.
#' @seealso \code{\link{combineTwasWeights}}
#' @examples
#' data(qtlFineMappingExample)
#' combineFineMappingResults(qtlFineMappingExample)
#' @export
combineFineMappingResults <- function(..., ldSketch = NULL) {
    parts <- .asCombineList(
        list(...),
        "FineMappingResultBase",
        "combineFineMappingResults"
    )
    .combineTupleCollections(parts, ldSketch, "combineFineMappingResults")
}


# The selector arguments are named rather than carried in `...`: an omitted
# one must still fall through to the extractor's own default, so NULL entries
# are dropped instead of forwarded. discard(is.null), not compact(): compact()
# would also drop a legitimately zero-length `contexts` / `samples`.
.fmResidPheno <- function(
    x,
    residualizationArgs,
    contexts = NULL,
    traitId = NULL,
    naAction = NULL
) {
    sel <- discard(
        list(contexts = contexts, traitId = traitId, naAction = naAction),
        is.null
    )
    exec(
        getResidualizedPhenotypes,
        !!!c(
            list(
                x = x,
                residualizationArgs = residualizationArgs %||%
                    ResidualizationParam()
            ),
            sel
        )
    )
}

.fmResidGeno <- function(
    x,
    residualizationArgs,
    contexts = NULL,
    traitId = NULL,
    region = NULL,
    cisWindow = NULL,
    samples = NULL
) {
    sel <- discard(
        list(
            contexts = contexts,
            traitId = traitId,
            region = region,
            cisWindow = cisWindow,
            samples = samples
        ),
        is.null
    )
    exec(
        getResidualizedGenotypes,
        !!!c(
            list(
                x = x,
                residualizationArgs = residualizationArgs %||%
                    ResidualizationParam()
            ),
            sel
        )
    )
}

# Directional effect-allele (A1) frequency for the variants in a fitted
# genotype block `X` (samples x variants, post-residualization and post
# sample-intersection). Re-extracts the allele frequency from the dataset
# `data` over the SAME selection used to build `X` and aligns it to
# `colnames(X)`; variants `getAf` does not return (e.g. dropped by a
# borderline MAF re-check on the final sample set) come back as NA. Returns
# NULL when `X` is empty or the dataset exposes no `getAf` (non-QtlDataset
# sources whose entries already carry `af`). The branch mirrors the
# `.fmResidGeno` call that built `X`: `region`-driven when a joint range is
# given, else `traitId` + `cisWindow` for the cis window.
#' @importFrom rlang try_fetch
.fmAfForX <- function(
    data,
    X,
    traitId = NULL,
    region = NULL,
    cisWindow = NULL
) {
    if (is.null(X) || ncol(X) == 0L || nrow(X) == 0L) {
        return(NULL)
    }
    if (!is(data, "QtlDataset")) {
        return(NULL)
    }
    afAll <- try_fetch(
        if (is.null(region)) {
            getAf(
                data,
                traitId = traitId,
                cisWindow = cisWindow,
                samples = rownames(X)
            )
        } else {
            getAf(data, region = region, samples = rownames(X))
        },
        error = function(cnd) NULL
    )
    if (is.null(afAll) || length(afAll) == 0L) {
        return(NULL)
    }
    unname(afAll[colnames(X)])
}

.fmPostprocessOne <- function(
    fit,
    method,
    dataX,
    dataY,
    credibleSetArgs,
    fitRetention,
    csInput = NULL,
    af = NULL,
    n = NULL,
    region = NULL,
    conditionIdx = NULL
) {
    .fmRunPostprocess(
        fit,
        method,
        dataX,
        dataY,
        region,
        af,
        n,
        csInput,
        conditionIdx,
        credibleSetArgs = credibleSetArgs,
        fitRetention = fitRetention
    )
}

# Run postprocessFinemappingFits for a single (method -> fit) mapping with the
# resolved defaults `d`, then format + validate the FineMappingRow payload.
# @noRd
.fmRunPostprocess <- function(
    fit,
    method,
    dataX,
    dataY,
    region,
    af,
    n,
    csInput,
    conditionIdx,
    credibleSetArgs,
    fitRetention
) {
    post <- postprocessFinemappingFits(
        fits = set_names(list(fit), method),
        dataX = dataX,
        dataY = dataY,
        af = af,
        n = n,
        region = region,
        csInput = csInput,
        conditionIdx = conditionIdx,
        # Both were threaded all the way down here and then not passed on,
        # so every caller's credible-set settings (coverage, purity,
        # signalCutoff, perCsColumns, ...) and retention level were replaced
        # by postprocessFinemappingFits()'s own defaults.
        credibleSetArgs = credibleSetArgs,
        fitRetention = fitRetention
    )
    out <- formatFinemappingOutput(post, primaryMethod = method)
    # `formatFinemappingOutput` returns $finemappingEntry as a bare row
    # payload (variants + susieFit + cvResult) per the helper's contract.
    if (!methods::is(out$finemappingEntry, "FineMappingRow")) {
        msg <- glue(
            ".fmPostprocessOne: postprocess output did not carry a ",
            "fine-mapping row - check pecotmr internal contract."
        )
        abort(msg)
    }
    out$finemappingEntry
}

# --- Multi-region (jointRegions) helpers ------------------------------------

# Resolve the per-trait X windows from a (region, jointRegions) pair. The cis
# path (region NULL) is a single trait-derived block; an explicit `region` is
# taken literally as one joint block (jointRegions=TRUE -> concatenated
# genotypes) or one block per range (jointRegions=FALSE -> independent fits
# merged downstream). Shared by the QtlDataset / MultiStudyQtlDataset
# fineMapping & twas methods.
#' @keywords internal
.makeXRegions <- function(region, jointRegions) {
    # Accept a "chr:start-end" string / one-row data.frame as well as a GRanges
    # (a GRanges passes through unchanged), so pipeline callers need not
    # pre-parse.
    if (!is.null(region)) {
        region <- .asGRegion(region)
    }
    if (is.null(region)) {
        list(NULL)
    } else if (isTRUE(jointRegions)) {
        list(region)
    } else {
        map(seq_along(region), .fmNthRegion, region = region)
    }
}

# Rows of a TwasWeights whose method is mr.mash and whose identity matches
# every axis the caller FIXED. A NULL axis means "do not filter on it".
# @noRd
.fmMrmashSelector <- function(mrmashPrior, study, trait, context) {
    axes <- compact(list(study = study, trait = trait, context = context))
    reduce(
        names(axes),
        .fmMrmashNarrow,
        mrmashPrior = mrmashPrior,
        axes = axes,
        .init = as.character(mrmashPrior$method) == "mrmash"
    )
}

# `.tupleColumn()` not `[[`: on a RangedTupleList `[[` extracts an ELEMENT
# (the variant set), while the identity axes live in mcols.
# @noRd
.fmMrmashNarrow <- function(sel, axis, mrmashPrior, axes) {
    sel & as.character(.tupleColumn(mrmashPrior, axis)) == axes[[axis]]
}

# Locate the retained mr.mash fit payload {dataDrivenPriorMatrices, w0, V} for
# one (study, trait[, context]) inside a `TwasWeights` collection from a prior
# mr.mash twasWeightsPipeline run (the producer side of the mvSuSiE data-driven
# prior). The joint fit is attached to a single mrmash row of the group (the
# other rows carry fits = NULL), so scan the matching mrmash rows and return the
# first non-NULL payload. The fit may span more conditions than the mvsusie
# block fits -- `.buildMvsusieReweightedPrior(include_indices=)` subsets it.
#
# Each axis is optional: a NULL axis is NOT filtered (match-any). A joint fit is
# shared across all its per-context rows, so the consumer fixes the constant
# axes and leaves the jointed (varying) axis NULL -- e.g. cross-context mvsusie
# keys on (study, trait) with context = NULL; cross-trait keys on (study,
# context) with trait = NULL (see .jointPriorKey). Returns NULL when no
# TwasWeights is supplied or it carries no matching mr.mash fit (caller falls
# back to the canonical prior).
# @noRd
.fmLookupMrmashFit <- function(
    mrmashPrior,
    study = NULL,
    trait = NULL,
    context = NULL
) {
    if (is.null(mrmashPrior)) {
        return(NULL)
    }
    # Each per-context mr.mash row of a joint group carries the SHARED joint
    # fit, so the consumer matches the FIXED axes and leaves the jointed axis
    # NULL (match-any). study/trait/context = NULL means "skip that axis".
    sel <- .fmMrmashSelector(mrmashPrior, study, trait, context)
    for (i in which(sel)) {
        f <- getFits(.twrRowParts(mrmashPrior, i))
        if (!is.null(f)) return(f)
    }
    NULL
}

# Locate the retained per-fold mr.mash CV payload for one (study, trait[,
# context]) inside a `TwasWeights` collection: the mrmash entry's `cvResult`,
# carrying `foldFits` (per-fold lean payloads) + `samplePartition` (the folds
# the per-fold priors were computed on). These let the mvSuSiE CV use an honest
# per-fold prior instead of reusing the full-data prior on every fold. Returns
# NULL when no TwasWeights / no matching mr.mash CV result with fold fits.
# @noRd
.fmLookupMrmashCv <- function(
    mrmashPrior,
    study = NULL,
    trait = NULL,
    context = NULL
) {
    if (is.null(mrmashPrior)) {
        return(NULL)
    }
    sel <- .fmMrmashSelector(mrmashPrior, study, trait, context)
    for (i in which(sel)) {
        cv <- getCvResult(.twrRowParts(mrmashPrior, i))
        if (!is.null(cv) && !is.null(cv$foldFits)) return(cv)
    }
    NULL
}

# Build the per-fold mvSuSiE reweighted priors for cross-validation from a
# TwasWeights mr.mash CV payload (`mvCv` from .fmLookupMrmashCv). For each fold:
# * full per-fold fit (carries its own w0) -> reweight that fit [mode B] *
# prior-only stub (U but no w0) -> reuse `fullFitParts` w0/V with the fold's U
# via overrideU [mode C]
# Returns a list named by fold id (as character, matching samplePartition$Fold);
# each element a list(priorVariance, residualVariance). NULL if no fold fits.
# @noRd
.fmBuildMvsusiePriorCv <- function(
    mvCv,
    fullFitParts,
    conditionNames,
    weightsTol = 1e-10
) {
    if (is.null(mvCv) || is.null(mvCv$foldFits)) {
        return(NULL)
    }
    foldFits <- mvCv$foldFits
    sp <- mvCv$samplePartition
    foldIds <- if (!is.null(sp)) sort(unique(sp$Fold)) else seq_along(foldFits)
    # A fold with no fit keeps its NULL slot, so the result stays aligned
    # with `foldIds` -- which is what `out[[i]] <- ...` left behind.
    set_names(
        map(
            seq_along(foldIds),
            .fmFoldPriorAt,
            foldIds = foldIds,
            foldFits = foldFits,
            fullFitParts = fullFitParts,
            conditionNames = conditionNames,
            weightsTol = weightsTol
        ),
        as.character(foldIds)
    )
}

# Fold `i`'s fit, matched by name ("fold_<id>") when available, else by
# position.
# @noRd
.fmFoldFitAt <- function(i, foldIds, foldFits) {
    nm <- str_c("fold_", foldIds[[i]])
    if (!is.null(names(foldFits)) && is_in(nm, names(foldFits))) {
        return(foldFits[[nm]])
    }
    if (length(foldFits) >= i) {
        return(foldFits[[i]])
    }
    NULL
}

# Fold `i`'s reweighted mixture prior, or NULL when the fold has no fit. A
# fold that carries its own w0 defines the prior; otherwise the full-data fit
# does, with the fold's data-driven matrices substituted in.
# @noRd
.fmFoldPriorAt <- function(
    i,
    foldIds,
    foldFits,
    fullFitParts,
    conditionNames,
    weightsTol
) {
    ff <- .fmFoldFitAt(i, foldIds, foldFits)
    if (is.null(ff)) {
        return(NULL)
    }
    if (!is.null(ff$w0)) {
        return(.buildMvsusieReweightedPrior(ff, conditionNames, weightsTol))
    }
    .buildMvsusieReweightedPrior(
        fullFitParts,
        conditionNames,
        weightsTol,
        overrideU = ff$dataDrivenPriorMatrices
    )
}

# PCA-reduce a (samples x traits) phenotype matrix to its top `nPCs` principal
# component scores, for the `usePCA` top-PC susie path. Centers + scales
# (matching the legacy fsusie.R susie_on_top_pc), dropping incomplete rows and
# zero-variance traits first (prcomp requires complete, non-degenerate columns).
# Returns a (samples x k) score matrix, k = min(nPCs, usable traits), columns
# named topPC1..topPCk and rows keyed by sample; NULL when < 2 usable traits or
# samples (single-trait -> PCA undefined, so the caller skips).
# @noRd
.fmTopPcScores <- function(Y, nPCs) {
    if (is.null(dim(Y)) || ncol(Y) < 2L) {
        return(NULL)
    }
    complete <- Y[stats::complete.cases(Y), , drop = FALSE]
    if (nrow(complete) < 2L) {
        return(NULL)
    }
    varying <- complete[,
        apply(complete, 2L, stats::var) > 0,
        drop = FALSE
    ]
    if (ncol(varying) < 2L) {
        return(NULL)
    }
    scores <- stats::prcomp(varying, center = TRUE, scale. = TRUE)$x
    k <- min(as.integer(nPCs), ncol(scores))
    if (k < 1L) {
        return(NULL)
    }
    `colnames<-`(scores[, seq_len(k), drop = FALSE], str_c("topPC", seq_len(k)))
}

# Is a signal screen enabled? Any spec that .asScreen resolves to a screen
# object (a non-zero pip cutoff or a resolved metric) activates it; this gates
# the extra screening extraction so the default (no screen) costs nothing.
# @noRd
.fmScreenActive <- function(screen) {
    !is.null(.asScreen(screen))
}

# Per-condition SER pre-screen for a joint (multi-context / multi-trait) fit:
# returns a logical vector over the columns of `Y` (the conditions) marking
# which show single-effect signal. The multivariate analog of `.fmSerScreen`
# and a port of the deleted `skipConditions`: callers drop the FALSE columns
# (null contexts / traits) before the joint mvSuSiE fit.
# @noRd
.fmSerScreenColumns <- function(X, Y, screen) {
    map_lgl(seq_len(ncol(Y)), .fmSerScreenColumn, X = X, Y = Y, screen = screen)
}

# .fmFitXBlock / .fmFitRssBlock (per-block SuSiE dispatch) now live in
# fineMappingWrappers.R with the other method-fitting wrappers.

# Extract integer credible-set indices from a "<method>_<idx>" vector.
.fmCsIdx <- function(csVec) {
    suppressWarnings(as.integer(str_replace(
        as.character(csVec),
        "^.*_([0-9]+)$",
        "\\1"
    )))
}

# Re-number credible-set membership labels by `offset`, preserving the
# "<method>_0" (not-in-any-CS) sentinel.
.fmRelabelCs <- function(csVec, offset) {
    csVec <- as.character(csVec)
    if (offset == 0L) {
        return(csVec)
    }
    parts <- str_match(csVec, "^(.*)_([0-9]+)$")
    map_chr(
        seq_along(csVec),
        .fmRelabelCsOne,
        csVec = csVec,
        parts = parts,
        offset = offset
    )
}

# Merge per-region FineMappingRow payloads (same study/context/trait/method,
# independent fits) into one entry: concatenate variants and topLoci rows,
# renumber credible sets so per-region indices do not collide, and keep the
# per-region SuSiE fits as a named list in `susieFit` (consumers needing a
# single fit must iterate the list).
.fmMergeEntries <- function(entries) {
    entries <- entries[!map_lgl(entries, is.null)]
    if (length(entries) == 0L) {
        return(NULL)
    }
    if (length(entries) == 1L) {
        return(entries[[1L]])
    }
    variantIds <- list_c(map(entries, .fmEntryVariantIds))
    tls <- map(entries, .fmEntryTopLoci)
    allNames <- unique(list_c(map(tls, names)))
    csCols <- allNames[str_detect(allNames, "^cs_[0-9]+$")]
    topLoci <- bind_rows(.fmRenumberCs(tls, csCols))
    susieFit <- set_names(
        map(entries, .fmEntrySusieFit),
        str_c("region", seq_along(entries))
    )
    # Per-region CV partitions/predictions kept under region* names so a
    # multi-region entry retains each block's CV (NULL when no region had CV).
    cvList <- set_names(
        map(entries, .fmEntryCvResult),
        str_c("region", seq_along(entries))
    )
    cvResult <- if (all(map_lgl(cvList, is.null))) NULL else cvList
    fineMappingRow(
        variantIds = variantIds,
        susieFit = susieFit,
        topLoci = topLoci,
        cvResult = cvResult
    )
}

# Renumber per-region credible-set columns so region-local cs indices do not
# collide once concatenated: each region's indices are shifted by the running
# max of the regions before it (a sequential offset fold, per cs_<coverage>
# column).
# @noRd
# The highest credible-set index table `tl` uses in column `cc`, or 0 when it
# has no such column.
# @noRd
.fmCsMaxIn <- function(tl, cc) {
    if (!is_in(cc, names(tl))) {
        return(0L)
    }
    max(c(0L, .fmCsIdx(tl[[cc]])), na.rm = TRUE)
}

# How many credible sets precede each element, given per-element counts.
# @noRd
.fmExclusiveCumsum <- function(counts) {
    cumsum(c(0L, counts))[seq_along(counts)]
}

# @noRd
.fmColumnOffsets <- function(cc, tls) {
    .fmExclusiveCumsum(map_int(tls, .fmCsMaxIn, cc = cc))
}

# @noRd
.fmTableHasCol <- function(cc, tl) {
    is_in(cc, names(tl))
}

# @noRd
.fmRelabelColumn <- function(cc, tl, offsets, i) {
    .fmRelabelCs(tl[[cc]], offsets[[cc]][[i]])
}

# @noRd
.fmRenumberTable <- function(i, tls, csCols, offsets) {
    tl <- tls[[i]]
    present <- keep(csCols, .fmTableHasCol, tl = tl)
    if (length(present) == 0L) {
        return(tl)
    }
    mutate(
        tl,
        !!!set_names(
            map(present, .fmRelabelColumn, tl = tl, offsets = offsets, i = i),
            present
        )
    )
}

.fmRenumberCs <- function(tls, csCols) {
    # A table's offset in a column is how many credible sets the tables
    # before it contributed there -- a cumulative count, so it is known up
    # front instead of being carried through the walk.
    offsets <- set_names(map(csCols, .fmColumnOffsets, tls = tls), csCols)
    map(
        seq_along(tls),
        .fmRenumberTable,
        tls = tls,
        csCols = csCols,
        offsets = offsets
    )
}

# Merge per-method user kwargs onto a base arg list. `userArgs` is the
# per-token kwargs supplied by the caller (e.g. `list(L = 1, refine =
# FALSE)`); the capability table's `args` default fills in any keys the
# user did not set. User-supplied values always win over base, capability
# defaults, and chain-derived args. Returns the merged list.
# @noRd
# Split a merged argument list into the target's own formals plus a
# `methodArgs` remainder. Callers build one flat list (base args, capability
# defaults, user overrides); this routes each name to the right place, so a
# tool option reaches the tool and an unknown name errors at the call instead
# of vanishing. A target with no `methodArgs` formal gets the list unchanged.
# @noRd
.splitMethodArgs <- function(fn, args) {
    target <- match.fun(fn)
    fm <- names(formals(target))
    if (!is_in("methodArgs", fm)) {
        return(args)
    }
    isFormal <- is_in(names(args), setdiff(fm, "methodArgs"))
    c(
        args[isFormal],
        list(methodArgs = .methodConfigLike(target, args[!isFormal]))
    )
}

# Wrap the remainder in the record the target's `methodArgs` expects, so a
# target that demands a constructor gets one rather than the bare list it
# would refuse. The constructor is read off the target's own default --
# `methodArgs = MvsusieOptions()` names it -- so this needs no registry and
# cannot drift from the signature. A target still defaulting to `list()`
# gets the plain list.
# @noRd
.methodConfigLike <- function(target, extra) {
    d <- formals(target)$methodArgs
    if (!is.call(d) || !is.name(d[[1L]])) {
        return(extra)
    }
    exec(match.fun(as.character(d[[1L]])), !!!extra)
}

.fmMergeUserArgs <- function(baseArgs, token, userArgs = NULL) {
    if (is.null(userArgs)) {
        userArgs <- list()
    }
    info <- .fineMappingMethodCapabilities[[token]]
    capDefaults <- if (!is.null(info) && !is.null(info$args)) {
        info$args
    } else {
        list()
    }
    # Order matters: base < capability defaults < user overrides.
    withCaps <- if (length(capDefaults) > 0L) {
        list_modify(baseArgs, !!!compact(capDefaults))
    } else {
        baseArgs
    }
    if (length(userArgs) == 0L) {
        return(withCaps)
    }
    list_modify(withCaps, !!!compact(userArgs))
}

# .fmFitSusieIndiv / .fmFitSusieRss / .fmFitSusieSer (single SuSiE fits) now
# live in fineMappingWrappers.R with the other method-fitting wrappers.

# Extract variant ids + Z + (median) N from a single QtlSumStats /
# GwasSumStats entry GRanges. Errors when Z or N is missing. Wraps the
# shared `.entryToSumstatDf` helper (R/sumstatsQc.R).
# @noRd
.fmExtractZn <- function(gr, label, ldSketch = NULL, cutoffs = NULL) {
    allDf <- .entryToSumstatDf(gr, require = c("SNP", "Z", "N"), label = label)
    # Filtered HERE rather than at the LD build: z, the LD matrix and the
    # allele frequencies are all keyed off `variantIds`, so narrowing the id
    # set at its source keeps them aligned by construction instead of by three
    # subsetting steps staying in step with one another.
    keep <- .panelKeepMask(allDf$variant_id, ldSketch, cutoffs, label)
    df <- allDf[keep, , drop = FALSE]
    list(
        variantIds = df$variant_id,
        z = df$z,
        # Scalar block N the RSS fit consumes (susie_rss takes a single n) ...
        n = stats::median(df$N, na.rm = TRUE),
        # ... and the per-variant effective N, aligned to `variantIds`, used
        # only for the reporting-only top_loci$N column (never the fit).
        nVar = df$N
    )
}

# Whether entry `i` of a QC'd SumStats was deliberately screened out (and why),
# so fineMappingPipeline can skip it gracefully instead of erroring on a
# 0-variant entry. `summaryStatsQc(pipCutoffToSkip = ...)` empties a no-signal
# region and records qcInfo$entryAudit[[i]]$pipScreenSkipped (+
# pipScreenReason);
# an entry may also be empty for other reasons. Returns list(skipped, reason).
.fmEntrySkipInfo <- function(data, i) {
    ea <- try_fetch(getQcInfo(data)$entryAudit[[i]], error = function(cnd) NULL)
    screened <- isTRUE(ea$pipScreenSkipped)
    entry <- .collectionEntry(data, i)
    empty <- is.null(entry) || length(entry) == 0L
    reason <-
        if (
            !is.null(ea$pipScreenReason) &&
                str_length(as.character(ea$pipScreenReason)) > 0L
        ) {
            as.character(ea$pipScreenReason)
        } else if (screened) {
            "no signal above the PIP pre-screen cutoff"
        } else if (empty) {
            "empty entry (no variants)"
        } else {
            NA_character_
        }
    list(skipped = isTRUE(screened || empty), reason = reason)
}

# =============================================================================
# Per-fold cross-validation of fine-mapping methods
# -----------------------------------------------------------------------------
# fineMappingPipeline mirrors twasWeightsPipeline's cross-validation: when
# cvFolds > 1, each fine-mapping method is refit on the training samples of
# every fold, its weights extracted and used to predict the held-out samples,
# yielding out-of-fold predictions + per-outcome metrics. The partition and
# predictions are stored on each FineMappingRow's cvResult slot so
# twasWeightsPipeline can (a) reuse the identical fold partition and (b) feed
# fine-mapping's own cross-validated predictions straight into the SR-TWAS
# ensemble instead of recomputing them. Output shape mirrors twasWeightsCv()
# (samplePartition + per-method <key>_predicted / <key>_performance), keyed by
# the TWAS snake method name (adapter methodKey) for a drop-in merge.
# =============================================================================

# Snake method key (e.g. "susie_inf") for a fine-mapping token, taken from the
# shared adapter registry so fineMapping CV keys match the TwasWeights `method`
# column and twasWeightsCv()'s prediction keys.
# @noRd
.fmTwasMethodKey <- function(token) {
    adapter <- .twasFineMappingMethodAdapters[[token]]
    if (is.null(adapter)) {
        return(token)
    }
    str_remove(adapter$methodKey, "_weights$")
}

# Coerce a weight vector to a single-column matrix (rows named by the vector's
# names); pass matrices through unchanged.
# @noRd
.fmAsMat <- function(w) {
    if (is.matrix(w)) {
        return(w)
    }
    matrix(w, ncol = 1L, dimnames = list(names(w), NULL))
}

# Fit one fine-mapping method on (Xtr, Ytr) for a CV fold and return a
# variants x outcomes weight matrix (rownames = colnames(Xtr)). susie-family
# tokens are fit independently (no chained init) per fold, matching
# twasWeightsCv's per-fold refit. Returns NULL on failure (caller skips it).
# @noRd
.fmFoldWeights <- function(
    token,
    Xtr,
    Ytr,
    coverage,
    userArgs,
    pos,
    mvPrior = NULL
) {
    if (is_in(token, c("susie", "susieInf", "susieAsh"))) {
        return(.fmFoldWeightsSusie(token, Xtr, Ytr, coverage, userArgs))
    }
    if (token == "mvsusie") {
        return(.fmFoldWeightsMv(Xtr, Ytr, coverage, userArgs, mvPrior))
    }
    if (token == "fsusie") {
        return(.fmFoldWeightsFsusie(Xtr, Ytr, pos, userArgs))
    }
    NULL
}

# The minimal fold-fit payload: exactly the fields the weight extractors
# read, so K folds x (study, context, trait) does not retain K full SuSiE fits.
# `trimFinemappingFit()` is the credible-set-driven trimmer for the fit stored
# on a row; a fold fit has no credible sets to key on and is only ever used to
# re-extract weights, so it keeps a smaller set.
# @noRd
.fmLeanFoldFit <- function(fit, token) {
    if (is.null(fit)) {
        return(NULL)
    }
    keep <- c(
        "pip",
        "alpha",
        "mu",
        "X_column_scale_factors",
        "theta",
        "coef",
        "V"
    )
    lean <- fit[intersect(keep, names(fit))]
    `class<-`(lean, unique(c(token, class(fit))))
}

# Per-fold univariate-susie-family weights (susie / susieInf / susieAsh).
# @noRd
.fmFoldWeightsSusie <- function(token, Xtr, Ytr, coverage, userArgs) {
    y <- if (is.matrix(Ytr)) Ytr[, 1L] else Ytr
    fit <- .fmFitSusieIndiv(
        Xtr,
        y,
        token,
        coverage = coverage,
        userArgs = userArgs
    )
    w <- switch(
        token,
        susie = susieWeights(susieFit = fit),
        susieInf = susieInfWeights(susieInfFit = fit),
        susieAsh = susieAshWeights(susieAshFit = fit)
    )
    out <- .fmAsMat(set_names(as.numeric(w), colnames(Xtr)))
    `attr<-`(out, "fit", .fmLeanFoldFit(fit, token))
}

# Per-fold fsusie weights.
# @noRd
.fmFoldWeightsFsusie <- function(Xtr, Ytr, pos, userArgs) {
    fsArgs <- .fmMergeUserArgs(
        list(X = Xtr, Y = Ytr, pos = pos),
        "fsusie",
        userArgs
    )
    fit <- exec(fitFsusie, !!!.splitMethodArgs(fitFsusie, fsArgs))
    W <- fsusieWeights(fsusieFit = fit, variantIds = colnames(Xtr))
    out <- as.matrix(W)
    # fSuSiE cannot be re-extracted from a trimmed fit, so the precomputed
    # weight matrix rides along as `coef` the same way the row-level trimmer
    # does it.
    withCoef <- list_assign(fit, coef = out)
    `attr<-`(out, "fit", .fmLeanFoldFit(withCoef, "fsusie"))
}

# Per-fold fine-mapping fit for the CV engine. `ctx` carries mvPrior, mvPriorCv,
# tokens, coverage, methodArgs, pos, verbose. Weights keyed by canonical method
# key.
# @noRd
.fmFitFold <- function(Xtr, Ytr, j, ctx) {
    mvPrior <- ctx$mvPrior
    mvPriorCv <- ctx$mvPriorCv
    tokens <- ctx$tokens
    coverage <- ctx$coverage
    methodArgs <- ctx$methodArgs
    pos <- ctx$pos
    verbose <- ctx$verbose
    # Honest per-fold mvSuSiE prior when supplied (the fold's own
    # mr.mash-derived prior); otherwise the single full-data prior is reused on
    # every fold.
    mvPriorThisFold <- if (!is.null(mvPriorCv)) {
        p <- mvPriorCv[[as.character(j)]]
        if (is.null(p)) mvPrior else p
    } else {
        mvPrior
    }
    keys <- map_chr(tokens, .fmTwasMethodKey)
    results <- map(
        tokens,
        .fmFoldTokenResult,
        Xtr = Xtr,
        Ytr = Ytr,
        coverage = coverage,
        methodArgs = methodArgs,
        pos = pos,
        mvPriorThisFold = mvPriorThisFold,
        j = j,
        verbose = verbose
    )
    list(
        weights = set_names(map(results, "weights"), keys),
        fits = compact(set_names(map(results, "fit"), keys))
    )
}

# One token's fold result. The CV engine's fitFold contract is
# list(weights = <key -> matrix>, fits = <key -> fitted model>); the fold
# helpers attach their (lean) fit to the weight matrix, so it is split back
# out here rather than stripped in place.
# @noRd
.fmFoldTokenResult <- function(
    tk,
    Xtr,
    Ytr,
    coverage,
    methodArgs,
    pos,
    mvPriorThisFold,
    j,
    verbose
) {
    w <- try_fetch(
        .fmFoldWeights(
            tk,
            Xtr,
            Ytr,
            coverage,
            methodArgs[[tk]],
            pos,
            mvPriorThisFold
        ),
        error = function(cnd) {
            if (verbose >= 1) {
                eMsg <- conditionMessage(cnd)
                msg <- glue(
                    "  CV fold {j}, method {tk} failed: {eMsg}",
                    .trim = FALSE
                )
                inform(msg)
            }
            NULL
        }
    )
    if (is.null(w)) {
        return(list(weights = NULL, fit = NULL))
    }
    list(
        weights = `attr<-`(w, "fit", NULL),
        fit = attr(w, "fit")
    )
}

# Cross-validate a homogeneous set of fine-mapping `tokens` over (X, Y) via the
# shared .crossValidateWeights() engine. For univariate tokens Y is a single
# column; for mvsusie/fsusie Y carries one column per condition/feature (and
# fsusie additionally needs `pos`). Each token's per-fold fit is refit here; the
# engine owns partitioning, the (optionally parallel) fold loop, prediction, and
# the metric block. Returns list(samplePartition, prediction, performance),
# keyed identically to twasWeightsCv().
# @noRd
.fmWeightsCv <- function(
    X,
    Y,
    tokens,
    methodArgs,
    fold,
    samplePartition = NULL,
    coverage = 0.95,
    pos = NULL,
    verbose = 1,
    mvPrior = NULL,
    mvPriorCv = NULL,
    numThreads = 1,
    seed = NULL
) {
    if (length(tokens) == 0L) {
        return(NULL)
    }
    # Per-fold fit context passed to the shared engine's top-level fitter
    # (.fmFitFold). Weights are keyed by the canonical method key so
    # <key>_predicted / <key>_performance line up with the TwasWeights method
    # column.
    cvFitCtx <- list(
        mvPrior = mvPrior,
        mvPriorCv = mvPriorCv,
        tokens = tokens,
        coverage = coverage,
        methodArgs = methodArgs,
        pos = pos,
        verbose = verbose
    )
    res <- .crossValidateWeights(
        X,
        Y,
        fold = fold,
        samplePartitions = samplePartition,
        fitFold = .fmFitFold,
        fitFoldCtx = cvFitCtx,
        numThreads = numThreads,
        verbose = verbose,
        seed = seed,
        # Retained so twasWeightsCv can cross-validate a SuSiE-family method:
        # those wrappers extract from a supplied fit and never refit, so the
        # fold's own fit is the only thing that makes per-fold extraction
        # possible.
        fitRetention = "slim"
    )
    .fmCvResult(res)
}

# The cross-validation result, with the sample partition recorded on the
# fold fits. Recorded at the producer so a consumer can refuse fits trained
# on a different split; without it a mismatched partition would leak
# held-out samples into training and inflate the CV metrics.
# @noRd
.fmCvResult <- function(res) {
    foldFits <- if (is.null(res$foldFits)) {
        NULL
    } else {
        `attr<-`(
            res$foldFits,
            "partitionKey",
            .cvPartitionKey(res$samplePartition)
        )
    }
    list(
        samplePartition = res$samplePartition,
        prediction = res$prediction,
        performance = res$performance,
        foldFits = foldFits
    )
}

# Slice a full .fmWeightsCv() result down to one method's payload, keeping
# the shared samplePartition. Stored on that method's FineMappingRow.
# @noRd
.fmSliceCv <- function(cv, token) {
    if (is.null(cv)) {
        return(NULL)
    }
    key <- .fmTwasMethodKey(token)
    pk <- str_c(key, "_predicted")
    mk <- str_c(key, "_performance")
    if (!is_in(pk, names(cv$prediction))) {
        return(NULL)
    }
    list(
        samplePartition = cv$samplePartition,
        prediction = cv$prediction[pk],
        performance = cv$performance[mk],
        foldFits = .fmSliceFoldFits(cv$foldFits, key)
    )
}

# One method's per-fold fits, keyed by fold, from the full CV payload.
# @noRd
.fmSliceFoldFits <- function(foldFits, key) {
    if (is.null(foldFits) || length(foldFits) == 0L) {
        return(NULL)
    }
    out <- map(foldFits, key)
    if (all(map_lgl(out, is.null))) {
        return(NULL)
    }
    `attr<-`(out, "partitionKey", attr(foldFits, "partitionKey"))
}

# Rebuild a FineMappingRow with a cvResult attached (the class is immutable).
# @noRd
.fmAttachCv <- function(entry, cvResult) {
    if (is.null(entry) || is.null(cvResult)) {
        return(entry)
    }
    fineMappingRow(
        variantIds = .fmrPartsVariantIds(entry),
        susieFit = .fmrPartsSusieFit(entry),
        topLoci = .fmrPartsTopLoci(entry),
        cvResult = cvResult
    )
}

# =============================================================================
# Residualized genotype for one X window: the trait-derived cis block when
# `rg` is NULL, else the explicit region. Shared by the univariate and PCA
# dispatch paths.
# @noRd
.fmResidGenoBlock <- function(
    ctx,
    traitId,
    rg,
    samples,
    data,
    cisWindow,
    residualizationArgs
) {
    if (is.null(rg)) {
        .fmResidGeno(
            data,
            contexts = ctx,
            traitId = traitId,
            cisWindow = cisWindow,
            samples = samples,
            residualizationArgs = residualizationArgs
        )
    } else {
        .fmResidGeno(
            data,
            contexts = ctx,
            region = rg,
            samples = samples,
            residualizationArgs = residualizationArgs
        )
    }
}

# Merge per-window block entries into one row-record per token: a token is
# dropped when any window failed to fit it (NULL), otherwise the windows are
# merged via .fmMergeEntries. Reused for univariate tokens and the single PCA
# "susie" token (trait = the PC name).
# @noRd
.fmMergeTokenRows <- function(study, ctx, trait, tokens, blockEntries) {
    compact(map(
        tokens,
        .fmMergeTokenRow,
        study = study,
        ctx = ctx,
        trait = trait,
        blockEntries = blockEntries
    ))
}

# .fmFitXBlock with the many per-run knobs supplied from the config bundle `p`;
# callers pass only the block-specific arguments (design, response, tokens,
# addSusieInf, context label, trait/PC label, allele frequencies).
# @noRd
.fmFitXBlockP <- function(X, y, tokens, addSusieInf, ctx, label, afVec, cfg) {
    .fmFitXBlock(
        X,
        y,
        tokens,
        addSusieInf,
        cfg$methodArgs,
        cfg$verbose,
        ctx,
        label,
        cvFolds = cfg$crossValidationArgs$folds,
        cvThreads = cfg$crossValidationArgs$numThreads,
        samplePartition = cfg$crossValidationArgs$samplePartition,
        af = afVec,
        credibleSetArgs = cfg$credibleSetArgs,
        fitRetention = cfg$fitRetention,
        seed = cfg$seed
    )
}

# Fit one univariate X window for (ctx, trait): residualize genotype, align to
# Y, SER pre-screen, then .fmFitXBlock. Errors when too few shared samples
# (a hard data problem), returns list() when the window screens out.
# @noRd
.fmUnivBlockFit <- function(rg, ctx, tid, Y, toRun, cfg) {
    allX <- .fmResidGenoBlock(
        ctx,
        tid,
        rg,
        rownames(Y),
        data = cfg$data,
        cisWindow = cfg$cisWindow,
        residualizationArgs = cfg$residualizationArgs
    )
    common <- intersect(rownames(allX), rownames(Y))
    if (length(common) < 2L) {
        abort(glue(
            "fineMappingPipeline: too few shared samples between ",
            "residualized X and Y for (context='{ctx}', trait='{tid}')."
        ))
    }
    X <- allX[common, , drop = FALSE]
    yBlock <- Y[common, , drop = FALSE]
    y <- if (ncol(yBlock) > 1L) {
        yBlock[, 1L, drop = TRUE]
    } else {
        drop(yBlock)
    }
    if (!.fmSerScreen(X, y, cfg$screen)) {
        if (cfg$verbose >= 1) {
            inform(glue(
                "Skipping (context='{ctx}', trait='{tid}'): SER ",
                "pre-screen found no signal above the cutoff."
            ))
        }
        return(list())
    }
    afVec <- .fmAfForX(
        cfg$data,
        X,
        traitId = tid,
        region = rg,
        cisWindow = cfg$cisWindow
    )
    .fmFitXBlockP(X, y, toRun, cfg$addSusieInf, ctx, tid, afVec, cfg)
}

# All univariate row-records for one (context, trait): cache hits, then (when
# tokens remain) per-window fits merged per token.
# @noRd
.fmUnivTraitRows <- function(tid, ctx, cfg) {
    lookups <- map(
        cfg$univTokens,
        .fmUnivLookup,
        fineMappingResult = cfg$fineMappingResult,
        study = cfg$study,
        ctx = ctx,
        tid = tid
    )
    cachedRows <- map(
        keep(lookups, .fmHasCached),
        .fmUnivCachedRow,
        study = cfg$study,
        ctx = ctx,
        tid = tid
    )
    toRun <- map_chr(keep(lookups, .fmNotCached), "tk")
    if (length(toRun) == 0L) {
        return(cachedRows)
    }
    Y <- .fmResidPheno(
        cfg$data,
        contexts = ctx,
        traitId = tid,
        naAction = cfg$naAction,
        residualizationArgs = cfg$residualizationArgs
    )
    blockEntries <- map(
        cfg$xRegions,
        .fmUnivBlockFit,
        ctx = ctx,
        tid = tid,
        Y = Y,
        toRun = toRun,
        cfg = cfg
    )
    computed <- .fmMergeTokenRows(cfg$study, ctx, tid, toRun, blockEntries)
    c(cachedRows, computed)
}

# Fit one PCA X window for (ctx, pcName): residualize genotype, align to the
# PC scores, SER pre-screen, then univariate-susie .fmFitXBlock. Returns
# list() when too few shared samples or the window screens out (soft skip).
# @noRd
.fmPcaBlockFit <- function(rg, ctx, traits, pcName, pcY, samples, cfg) {
    X <- .fmResidGenoBlock(
        ctx,
        traits,
        rg,
        samples,
        data = cfg$data,
        cisWindow = cfg$cisWindow,
        residualizationArgs = cfg$residualizationArgs
    )
    common <- intersect(rownames(X), names(pcY))
    if (length(common) < 2L) {
        return(list())
    }
    Xb <- X[common, , drop = FALSE]
    if (!.fmSerScreen(Xb, pcY[common], cfg$screen)) {
        return(list())
    }
    afVec <- .fmAfForX(
        cfg$data,
        Xb,
        traitId = traits,
        region = rg,
        cisWindow = cfg$cisWindow
    )
    # The PC path never chains susieInf: a PC is not a trait whose fit the
    # susieInf stage is defined for.
    .fmFitXBlockP(Xb, pcY[common], "susie", FALSE, ctx, pcName, afVec, cfg)
}

# Row-records for one PC pseudo-trait: a cache hit, else per-window fits merged
# into the single "susie" token (trait = the PC name).
# @noRd
.fmPcaScoreRows <- function(pcName, ctx, traits, scores, cfg) {
    cached <- .fmCacheLookup(
        cfg$fineMappingResult,
        cfg$study,
        ctx,
        pcName,
        "susie"
    )
    if (!is.null(cached)) {
        return(list(.fmQtlRow(cfg$study, ctx, pcName, "susie", cached)))
    }
    pcY <- scores[, pcName]
    samples <- rownames(scores)
    blockEntries <- map(
        cfg$xRegions,
        .fmPcaBlockFit,
        ctx = ctx,
        traits = traits,
        pcName = pcName,
        pcY = pcY,
        samples = samples,
        cfg = cfg
    )
    .fmMergeTokenRows(cfg$study, ctx, pcName, "susie", blockEntries)
}

# All usePCA row-records for one context: PCA-reduce the multi-trait phenotype
# and fine-map each top PC. A single-trait context (or one with no usable PC
# scores) contributes nothing.
# @noRd
.fmPcaContextRows <- function(ctx, cfg) {
    traits <- cfg$perCtxTraits[[ctx]]
    if (length(traits) < 2L) {
        return(list())
    }
    Yctx <- .fmResidPheno(
        cfg$data,
        contexts = ctx,
        traitId = traits,
        naAction = cfg$naAction,
        residualizationArgs = cfg$residualizationArgs
    )
    scores <- .fmTopPcScores(Yctx, cfg$nPCs)
    if (is.null(scores)) {
        return(list())
    }
    if (cfg$verbose >= 1) {
        nPc <- ncol(scores)
        nTr <- length(traits)
        inform(glue(
            "usePCA: fine-mapping {nPc} top PC(s) of context='{ctx}' ",
            "({nTr} traits) ..."
        ))
    }
    list_flatten(map(
        colnames(scores),
        .fmPcaScoreRows,
        ctx = ctx,
        traits = traits,
        scores = scores,
        cfg = cfg
    ))
}

# QtlDataset method
# =============================================================================

# The settings every QtlDataset joint-phase function forwards to the shared
# joint dispatcher: the common fine-mapping record plus the scope and
# TWAS-weight fields only the joint engine uses.
# @noRd
.fmQdsJointCfg <- function(
    common,
    contexts,
    traitId,
    mrmashPrior,
    dataDrivenPriorWeightsCutoff
) {
    c(
        common,
        list(
            contexts = contexts,
            traitId = traitId,
            mrmashPrior = mrmashPrior,
            dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff
        )
    )
}

# Run the joint engine for a QtlDataset with the given spec + token set. Shared
# by the explicit-jointSpecification path and the auto-detected multivariate
# path; only the spec and token arguments differ between them.
# @noRd
.fmQdsJointDispatch <- function(jointSpec, tokens, methodArgs, cfg) {
    .fmDispatchJointSpecsQtlDataset(
        jointSpec,
        cfg$data,
        tokens,
        cfg$contexts,
        cfg$traitId,
        cfg$cisWindow,
        cfg$verbose,
        methodArgs = methodArgs,
        xRegions = cfg$xRegions,
        mrmashPrior = cfg$mrmashPrior,
        dataDrivenPriorWeightsCutoff = cfg$dataDrivenPriorWeightsCutoff,
        crossValidationArgs = cfg$crossValidationArgs,
        residualizationArgs = cfg$residualizationArgs,
        pipCutoffToSkip = cfg$screen,
        fineMappingResult = cfg$fineMappingResult,
        credibleSetArgs = cfg$credibleSetArgs,
        fitRetention = cfg$fitRetention,
        seed = cfg$seed
    )
}

# Resolve the signal screen, apply per-call filter overrides to a validated
# copy of the dataset, and derive the X windows. Rejects the region + cisWindow
# combination. Returns `p` extended with data / screen / xRegions.
# @noRd
.fmQdsResolveInputs <- function(
    data,
    signalScreenArgs,
    genotypeFilterArgs,
    region,
    cisWindow,
    jointRegions
) {
    .assertMethodParam(signalScreenArgs, "SignalScreenParam", "signalScreen")
    screen <- .screenResolve(signalScreenArgs)
    data <- .qtlApplyFilterOverrides(data, genotypeFilterArgs)
    if (!is.null(region) && !is.null(cisWindow)) {
        msg <- glue(
            "fineMappingPipeline(QtlDataset): specify either `region` or ",
            "`cisWindow`, not both. `cisWindow` expands each trait's own ",
            "coordinates, whereas `region` is the literal variant window."
        )
        abort(msg)
    }
    list(
        data = data,
        screen = screen,
        xRegions = .makeXRegions(region, jointRegions)
    )
}

# Normalize methods, capability-check them, then hand off to the joint phase.
# @noRd
# The joint-specification phase for the QtlDataset path: dispatch whatever the
# spec covers, narrow the tokens the joint fits consumed, and report whether
# anything is left for the per-tuple loop. `done = TRUE` means the joint fits
# are the whole answer.
#
# Deliberately mirrors .twasQdsJointPhase -- same responsibilities, same
# done / result shape -- so the two pipelines' joint handling reads as one
# pattern. The argument sets differ because the outputs do (credible sets
# here, weights there); the structure should not.
# @noRd
.fmQdsJointPhase <- function(parsedJointSpec, norm, cfg) {
    if (length(parsedJointSpec) == 0L) {
        return(list(
            done = FALSE,
            result = NULL,
            tokens = norm$tokens,
            methodArgs = norm$methodArgs
        ))
    }
    jointResult <- .fmQdsJointDispatch(
        parsedJointSpec,
        intersect(norm$tokens, c("mvsusie", "fsusie")),
        norm$methodArgs,
        cfg
    )
    # The joint dispatch consumed the multivariate tokens; the univariate
    # phase runs on whatever is left.
    narrowed <- .fmNarrowAfterJoint(TRUE, norm$tokens, norm$methodArgs)
    if (length(narrowed$tokens) == 0L) {
        if (is.null(jointResult)) {
            abort(glue(
                "fineMappingPipeline(QtlDataset): no joint fits produced. ",
                "Check that the jointSpecification scope intersects the ",
                "available studies / contexts / traits."
            ))
        }
        return(list(done = TRUE, result = jointResult))
    }
    list(
        done = FALSE,
        result = jointResult,
        tokens = narrowed$tokens,
        methodArgs = narrowed$methodArgs
    )
}

.fmQdsResolveTokens <- function(jointSpecification, methods, L, Lgreedy, cfg) {
    parsedJointSpec <- parseJointSpecification(jointSpecification, cfg$data)
    norm <- .fmNormalizeMethods(
        methods,
        inputKind = "QtlDataset",
        L = L,
        Lgreedy = Lgreedy
    )
    .fmCheckMethodCapabilities(norm$tokens, "QtlDataset")
    .fmCheckMethodArgsForInput(norm$methodArgs, "QtlDataset")
    .fmQdsJointPhase(parsedJointSpec, norm, cfg)
}

# The tokens (and their args) still owed a univariate fit after a joint
# dispatch ran. A run with no joint specification keeps everything.
# @noRd
.fmNarrowAfterJoint <- function(hadJointSpec, tokens, methodArgs) {
    if (!hadJointSpec) {
        return(list(tokens = tokens, methodArgs = methodArgs))
    }
    remaining <- setdiff(tokens, c("mvsusie", "fsusie"))
    list(tokens = remaining, methodArgs = methodArgs[remaining])
}

# Trait ids available in one context, intersected with the requested traitId or
# overlapping the requested region (mirrors twasWeightsPipeline).
# @noRd
.fmQdsTraitsForContext <- function(ctx, data, traitId, region) {
    se <- getPhenotypes(data, contexts = ctx)
    allIds <- rownames(se)
    ids <- if (!is.null(traitId)) {
        intersect(allIds, traitId)
    } else if (!is.null(region)) {
        allIds[IRanges::overlapsAny(
            SummarizedExperiment::rowRanges(se),
            region
        )]
    } else {
        allIds
    }
    ids
}

# Resolve the study, the contexts to use (validating any requested ones), and
# the per-context trait lists. Returns `p` extended with study / useCtx /
# perCtxTraits / nCtx / nTraits.
# @noRd
.fmQdsResolveContexts <- function(data, contexts, traitId, region) {
    allCtx <- getContexts(data)
    useCtx <- if (is.null(contexts)) {
        allCtx
    } else {
        bad <- setdiff(contexts, allCtx)
        if (length(bad) > 0L) {
            badStr <- str_flatten(bad, ", ")
            msg <- glue(
                "fineMappingPipeline(QtlDataset): unknown context(s): ",
                "{badStr}"
            )
            abort(msg)
        }
        contexts
    }
    perCtxTraits <- set_names(
        map(
            useCtx,
            .fmQdsTraitsForContext,
            data = data,
            traitId = traitId,
            region = region
        ),
        useCtx
    )
    allTraits <- unique(list_c(perCtxTraits))
    if (length(allTraits) == 0L) {
        abort("fineMappingPipeline(QtlDataset): no traits selected.")
    }
    list(
        study = getStudy(data),
        useCtx = useCtx,
        perCtxTraits = perCtxTraits,
        nCtx = length(useCtx),
        nTraits = length(allTraits)
    )
}

# Partition tokens into univariate / mvsusie / fsusie sets and validate the
# multivariate requirements (mvsusie needs multi-trait OR multi-context; fsusie
# needs multi-trait per context). Returns `p` extended with the three sets.
# @noRd
.fmQdsSplitTokens <- function(tokens, nCtx, nTraits) {
    univTokens <- tokens[!is_in(tokens, c("mvsusie", "fsusie"))]
    mvTokens <- tokens[tokens == "mvsusie"]
    fsTokens <- tokens[tokens == "fsusie"]
    if (length(mvTokens) > 0L && nCtx < 2L && nTraits < 2L) {
        msg <- glue(
            "fineMappingPipeline(QtlDataset): mvsusie requires multi-trait ",
            "or multi-context input (got {nTraits} trait(s) x ",
            "{nCtx} context(s))."
        )
        abort(msg)
    }
    if (length(fsTokens) > 0L && nTraits < 2L) {
        msg <- glue(
            "fineMappingPipeline(QtlDataset): fsusie requires multi-trait ",
            "input within a context (got {nTraits} trait(s))."
        )
        abort(msg)
    }
    list(
        univTokens = univTokens,
        mvTokens = mvTokens,
        fsTokens = fsTokens
    )
}

# Univariate + usePCA dispatch: each (context, trait) -> merged-per-token
# row-records; each multi-trait context's top PCs -> pseudo-trait rows.
# @noRd
# The per-run settings the univariate / PCA row builders all need: the common
# fine-mapping record plus the grid and method fields only the per-tuple
# dispatch uses. The chain threads a single record instead of two dozen
# arguments, and every field is named deliberately -- this is not a capture of
# the calling frame, so what travels is exactly what these lists declare.
# @noRd
.fmQdsRunCfg <- function(
    common,
    addSusieInf,
    methodArgs,
    nPCs,
    naAction,
    perCtxTraits,
    study,
    univTokens,
    useCtx,
    usePCA
) {
    c(
        common,
        list(
            addSusieInf = addSusieInf,
            methodArgs = methodArgs,
            nPCs = nPCs,
            naAction = naAction,
            perCtxTraits = perCtxTraits,
            study = study,
            univTokens = univTokens,
            useCtx = useCtx,
            usePCA = usePCA
        )
    )
}

.fmQdsDispatchRows <- function(cfg) {
    univRows <- if (length(cfg$univTokens) > 0L) {
        list_flatten(map(cfg$useCtx, .fmUnivContextRows, cfg = cfg))
    } else {
        list()
    }
    pcaRows <- if (isTRUE(cfg$usePCA)) {
        list_flatten(map(cfg$useCtx, .fmPcaContextRows, cfg = cfg))
    } else {
        list()
    }
    c(univRows, pcaRows)
}

# Multivariate dispatch via the joint engine (auto-detected shape) for
# mvsusie / fsusie WITHOUT an explicit jointSpecification, merged with any
# explicit-spec result already in `p$jointResult`.
# @noRd
.fmQdsAutoJoint <- function(
    cfg,
    mvTokens,
    fsTokens,
    jointResult,
    methodArgs,
    nCtx,
    nTraits
) {
    if (length(mvTokens) == 0L && length(fsTokens) == 0L) {
        return(jointResult)
    }
    autoJoint <- .fmQdsJointDispatch(
        .fmSynthesizeJointSpec(nCtx, nTraits),
        c(mvTokens, fsTokens),
        methodArgs,
        cfg
    )
    if (is.null(jointResult)) {
        autoJoint
    } else if (is.null(autoJoint)) {
        jointResult
    } else {
        .rbindFineMappingResult(jointResult, autoJoint, ldSketch = NULL)
    }
}

# Assemble the QtlDataset result from the per-tuple row-records (region =
# trait-anchored cis span) combined with the joint result. Errors only when
# neither path produced anything.
# @noRd
.fmQdsAssemble <- function(data, rows, jointResult) {
    rowContext <- map_chr(rows, "context")
    rowTrait <- map_chr(rows, "trait")
    perTupleResult <- if (length(rows) > 0L) {
        .fmBuildQtlResult(
            map_chr(rows, "study"),
            rowContext,
            rowTrait,
            map_chr(rows, "method"),
            map(rows, "entry"),
            traitPos = try_fetch(
                .anchorVector(data, rowContext, rowTrait, "traitPos"),
                error = function(cnd) NULL
            ),
            ldSketch = NULL
        )
    } else {
        NULL
    }
    if (is.null(jointResult)) {
        if (is.null(perTupleResult)) {
            msg <- glue(
                "fineMappingPipeline: no (study, context, trait, method) ",
                "tuples produced a fine-mapping result."
            )
            abort(msg)
        }
        return(perTupleResult)
    }
    if (is.null(perTupleResult)) {
        return(jointResult)
    }
    .rbindFineMappingResult(perTupleResult, jointResult, ldSketch = NULL)
}

# Resolve the context/trait tuples, split the tokens across the univariate,
# multivariate and functional paths, dispatch one row per tuple and assemble
# the result. Carved out of .fmPipelineQtlDataset purely to shorten it --
# every argument is one the entry point already holds, so the two still have
# to be read together.
# @noRd
.fmQdsRunAndAssemble <- function(
    data,
    commonCfg,
    jointCfg,
    contexts,
    traitId,
    region,
    tokens,
    methodArgs,
    jointResult,
    addSusieInf,
    nPCs,
    naAction,
    usePCA
) {
    ctxInfo <- .fmQdsResolveContexts(
        data = data,
        contexts = contexts,
        traitId = traitId,
        region = region
    )
    split <- .fmQdsSplitTokens(tokens, ctxInfo$nCtx, ctxInfo$nTraits)
    rows <- .fmQdsDispatchRows(.fmQdsRunCfg(
        commonCfg,
        addSusieInf = addSusieInf,
        methodArgs = methodArgs,
        nPCs = nPCs,
        naAction = naAction,
        perCtxTraits = ctxInfo$perCtxTraits,
        study = ctxInfo$study,
        univTokens = split$univTokens,
        useCtx = ctxInfo$useCtx,
        usePCA = usePCA
    ))
    .fmQdsAssemble(
        data,
        rows,
        .fmQdsAutoJoint(
            jointCfg,
            mvTokens = split$mvTokens,
            fsTokens = split$fsTokens,
            jointResult = jointResult,
            methodArgs = methodArgs,
            nCtx = ctxInfo$nCtx,
            nTraits = ctxInfo$nTraits
        )
    )
}

# QtlDataset fine-mapping worker. Each stage returns only the values it
# derives, and the next stage takes them as named arguments -- nothing rides
# on a captured environment. The susieInf chaining is applied downstream
# inside .fmFitXBlock / .fmFitRssBlock (which recompute
# .fmResolveSusieChain), so no chain config is threaded here.
# @noRd
.fmPipelineQtlDataset <- function(
    data,
    methods,
    contexts = NULL,
    traitId = NULL,
    region = NULL,
    cisWindow = NULL,
    # Per-call genotype-filter overrides; NULL = use the QtlDataset's
    # construct-time slot value (applied lazily at extraction).
    genotypeFilterArgs = GenotypeFilterParam(),
    jointRegions = FALSE,
    jointSpecification = NULL,
    addSusieInf = TRUE,
    credibleSetArgs = CredibleSetParam(),
    fineMappingResult = NULL,
    crossValidationArgs = CrossValidationParam(),
    residualizationArgs = ResidualizationParam(),
    signalScreenArgs = SignalScreenParam(),
    usePCA = FALSE,
    nPCs = 10L,
    seed = NULL,
    mrmashPrior = NULL,
    dataDrivenPriorWeightsCutoff = 1e-10,
    naAction = c("drop", "impute"),
    verbose = 1,
    fitRetention = c("slim", "full")
) {
    fitRetention <- arg_match(fitRetention)
    naAction <- arg_match(naAction, c("drop", "impute"))
    if (!is.null(seed)) {
        withr::local_seed(as.integer(seed))
    }
    # Each stage returns only the values it derives; nothing is grafted onto a
    # captured environment.
    resolved <- .fmQdsResolveInputs(
        data = data,
        signalScreenArgs = signalScreenArgs,
        genotypeFilterArgs = genotypeFilterArgs,
        region = region,
        cisWindow = cisWindow,
        jointRegions = jointRegions
    )
    data <- resolved$data
    screen <- resolved$screen
    xRegions <- resolved$xRegions
    # The bundle is the user-facing form; below this line the existing cfg
    # plumbing (and the joint engine that reads it) keeps plain scalars.
    cvCfg <- .cvResolve(crossValidationArgs)
    # The settings the joint phase and the per-tuple dispatch share, named
    # once so the two records built from it cannot disagree.
    commonCfg <- list(
        data = data,
        cisWindow = cisWindow,
        credibleSetArgs = credibleSetArgs,
        screen = screen,
        xRegions = xRegions,
        crossValidationArgs = cvCfg,
        residualizationArgs = residualizationArgs,
        fineMappingResult = fineMappingResult,
        fitRetention = fitRetention,
        verbose = verbose,
        seed = seed
    )
    jointCfg <- .fmQdsJointCfg(
        commonCfg,
        contexts = contexts,
        traitId = traitId,
        mrmashPrior = mrmashPrior,
        dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff
    )
    rt <- .fmQdsResolveTokens(
        jointSpecification,
        methods,
        credibleSetArgs$L %||% 10L,
        credibleSetArgs$Lgreedy,
        jointCfg
    )
    if (rt$done) {
        return(rt$result)
    }
    .fmQdsRunAndAssemble(
        data,
        commonCfg,
        jointCfg,
        contexts = contexts,
        traitId = traitId,
        region = region,
        tokens = rt$tokens,
        methodArgs = rt$methodArgs,
        jointResult = rt$result,
        addSusieInf = addSusieInf,
        nPCs = nPCs,
        naAction = naAction,
        usePCA = usePCA
    )
}

#' @rdname fineMappingPipeline
#' @importFrom purrr list_c list_flatten list_modify list_rbind list_assign
#' @importFrom purrr discard walk2 detect zap
#' @export
setMethod(
    "fineMappingPipeline",
    "QtlDataset",
    .fmPipelineQtlDataset
)

# =============================================================================
# MultiStudyQtlDataset method
# =============================================================================

# Per-embedded-study fine-mapping worker for .multiStudyPipelineDriver: recurse
# fineMappingPipeline on one QtlDataset with the individual-capable methods.
# `cfg` bundles the parent call's forwarded arguments.
# @noRd
.fmPerStudy <- function(qd, cfg) {
    m <- .fmFilterMethodsForKind(cfg$methods, "individualImpl")
    if (length(m) == 0L) {
        return(NULL)
    }
    # Every setting is named. `cfg` also carries panelFilterArgs and rssArgs,
    # which belong to the summary-statistics path: a study's own genotypes
    # are filtered by genotypeFilterArgs, so those two are deliberately not
    # forwarded here. A name this entry point does not accept would be an
    # "unused argument" error, so the only way this can go wrong is by
    # OMITTING a setting -- which test_fineMappingPipeline.R pins.
    fineMappingPipeline(
        data = qd,
        methods = m,
        jointSpecification = NULL,
        contexts = cfg$contexts,
        traitId = cfg$traitId,
        region = cfg$region,
        cisWindow = cfg$cisWindow,
        jointRegions = cfg$jointRegions,
        addSusieInf = cfg$addSusieInf,
        credibleSetArgs = cfg$credibleSetArgs,
        fineMappingResult = cfg$fineMappingResult,
        verbose = cfg$verbose,
        crossValidationArgs = cfg$crossValidationArgs,
        residualizationArgs = cfg$residualizationArgs,
        seed = cfg$seed,
        naAction = cfg$naAction,
        signalScreenArgs = cfg$signalScreenArgs,
        genotypeFilterArgs = cfg$genotypeFilterArgs,
        usePCA = cfg$usePCA,
        nPCs = cfg$nPCs,
        fitRetention = cfg$fitRetention,
        mrmashPrior = cfg$mrmashPrior,
        dataDrivenPriorWeightsCutoff = cfg$dataDrivenPriorWeightsCutoff
    )
}

# Embedded-sumstats fine-mapping worker for .multiStudyPipelineDriver: recurse
# fineMappingPipeline on the QtlSumStats with the sumstat-capable methods.
# @noRd
.fmSumStats <- function(ss, cfg) {
    m <- .fmFilterMethodsForKind(cfg$methods, "sumstatImpl")
    if (length(m) == 0L) {
        return(NULL)
    }
    # Summary statistics carry no genotypes and no per-study region
    # selection: the region is already baked into the input, and there is
    # nothing to residualize, screen, cross-validate or PCA. `cfg` carries
    # those settings for the individual-level sibling and they are not
    # forwarded here.
    fineMappingPipeline(
        data = ss,
        methods = m,
        jointSpecification = NULL,
        contexts = cfg$contexts,
        traitId = cfg$traitId,
        addSusieInf = cfg$addSusieInf,
        credibleSetArgs = cfg$credibleSetArgs,
        fineMappingResult = cfg$fineMappingResult,
        verbose = cfg$verbose,
        panelFilterArgs = cfg$panelFilterArgs,
        rssArgs = cfg$rssArgs,
        fitRetention = cfg$fitRetention,
        mrmashPrior = cfg$mrmashPrior,
        dataDrivenPriorWeightsCutoff = cfg$dataDrivenPriorWeightsCutoff
    )
}

# Validate the naAction choice and the region + cisWindow combination, then
# derive the X windows. Mirrors .fmQdsResolveInputs, minus the filter
# overrides the components apply themselves.
# @noRd
.fmMsResolveInputs <- function(naAction, region, cisWindow, jointRegions) {
    if (!is.null(region) && !is.null(cisWindow)) {
        msg <- glue(
            "fineMappingPipeline(MultiStudyQtlDataset): specify either ",
            "`region` or `cisWindow`, not both."
        )
        abort(msg)
    }
    list(
        naAction = arg_match(naAction, c("drop", "impute")),
        xRegions = .makeXRegions(region, jointRegions)
    )
}

# Run the joint engine for a MultiStudyQtlDataset with the given spec + token
# set.
# @noRd
.fmMsJointDispatch <- function(jointSpec, tokens, methodArgs, cfg) {
    .fmDispatchJointSpecsMultiStudy(
        jointSpec,
        cfg$data,
        tokens,
        cfg$contexts,
        cfg$traitId,
        cfg$cisWindow,
        cfg$verbose,
        methodArgs = methodArgs,
        xRegions = cfg$xRegions,
        mrmashPrior = cfg$mrmashPrior,
        dataDrivenPriorWeightsCutoff = cfg$dataDrivenPriorWeightsCutoff
    )
}

# The joint-specification phase for the MultiStudyQtlDataset path: dispatch
# whatever the spec covers, narrow the tokens the joint fits consumed, and
# report whether anything is left for the per-component recursion.
# `done = TRUE` means the joint fits are the whole answer.
#
# Deliberately mirrors .fmQdsJointPhase and .twasQdsJointPhase -- same
# responsibilities, same done / result shape.
# @noRd
.fmMsJointPhase <- function(parsedJointSpec, norm, cfg) {
    if (length(parsedJointSpec) == 0L) {
        return(list(
            done = FALSE,
            result = NULL,
            tokens = norm$tokens,
            methodArgs = norm$methodArgs
        ))
    }
    jointResult <- .fmMsJointDispatch(
        parsedJointSpec,
        intersect(norm$tokens, c("mvsusie", "fsusie")),
        norm$methodArgs,
        cfg
    )
    # The joint dispatch consumed the multivariate tokens; the univariate
    # phase runs on whatever is left.
    narrowed <- .fmNarrowAfterJoint(TRUE, norm$tokens, norm$methodArgs)
    if (length(narrowed$tokens) == 0L) {
        if (is.null(jointResult)) {
            abort(glue(
                "fineMappingPipeline(MultiStudyQtlDataset): no joint fits ",
                "produced. Check that the jointSpecification scope ",
                "intersects the available data."
            ))
        }
        return(list(done = TRUE, result = jointResult))
    }
    list(
        done = FALSE,
        result = jointResult,
        tokens = narrowed$tokens,
        methodArgs = narrowed$methodArgs
    )
}

# Resolve method tokens for a MultiStudyQtlDataset run and run any EXPLICIT
# jointSpecification (per-component axis dispatcher), removing the joint methods
# from the per-tuple recursion. Returns the still-pending tokens, the forwarded
# `methods` (kwargs-preserving), the joint result, and `done`.
# @noRd
.fmMsResolveTokens <- function(jointSpecification, methods, L, Lgreedy, cfg) {
    parsedJointSpec <- parseJointSpecification(jointSpecification, cfg$data)
    # A MultiStudyQtlDataset may hold individual-level studies AND summary
    # statistics, so `methods` is translated here -- where the dataset is in
    # hand -- and the Param is what travels into the per-component
    # recursion. Each component then dispatches on its own class and reads
    # the slot for its own path.
    methods <- .methodsParamForMulti(
        methods,
        "FineMappingMethodsParam",
        "fineMappingPipeline",
        !is.null(getSumStats(cfg$data))
    )
    norm <- .fmNormalizeMethods(
        methods,
        inputKind = "QtlDataset",
        L = L,
        Lgreedy = Lgreedy
    )
    .fmCheckMethodCapabilities(norm$tokens, "MultiStudyQtlDataset")
    .fmCheckMethodArgsForInput(norm$methodArgs, "MultiStudyQtlDataset")
    phase <- .fmMsJointPhase(parsedJointSpec, norm, cfg)
    if (phase$done) {
        return(phase)
    }
    # The per-component recursion re-enters fineMappingPipeline through
    # `methods=`, so a run with no joint spec forwards the caller's own
    # spelling untouched rather than the normalized tokens.
    list(
        done = FALSE,
        result = phase$result,
        tokens = phase$tokens,
        methods = if (length(parsedJointSpec) == 0L) {
            methods
        } else if (length(phase$methodArgs) > 0L) {
            phase$methodArgs
        } else {
            phase$tokens
        }
    )
}

# MultiStudyQtlDataset fine-mapping worker. After resolving the explicit
# joint spec it routes each remaining method to the components it supports,
# via the shared multi-study driver.
# @noRd
.fmPipelineMultiStudy <- function(
    data,
    methods,
    contexts = NULL,
    traitId = NULL,
    region = NULL,
    cisWindow = NULL,
    genotypeFilterArgs = GenotypeFilterParam(),
    panelFilterArgs = PanelFilterParam(),
    jointRegions = FALSE,
    jointSpecification = NULL,
    addSusieInf = TRUE,
    credibleSetArgs = CredibleSetParam(),
    fineMappingResult = NULL,
    mrmashPrior = NULL,
    dataDrivenPriorWeightsCutoff = 1e-10,
    crossValidationArgs = CrossValidationParam(),
    residualizationArgs = ResidualizationParam(),
    signalScreenArgs = SignalScreenParam(),
    rssArgs = SusieRssParam(),
    usePCA = FALSE,
    nPCs = 10L,
    seed = NULL,
    naAction = c("drop", "impute"),
    verbose = 1,
    fitRetention = c("slim", "full")
) {
    fitRetention <- arg_match(fitRetention)
    inputs <- .fmMsResolveInputs(naAction, region, cisWindow, jointRegions)
    naAction <- inputs$naAction
    # One record of the settings the joint phase and the shared dispatcher
    # both need, so each names them once.
    jointCfg <- list(
        data = data,
        contexts = contexts,
        traitId = traitId,
        cisWindow = cisWindow,
        credibleSetArgs = credibleSetArgs,
        verbose = verbose,
        xRegions = inputs$xRegions,
        mrmashPrior = mrmashPrior,
        dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff
    )
    rt <- .fmMsResolveTokens(
        jointSpecification,
        methods,
        credibleSetArgs$L %||% 10L,
        credibleSetArgs$Lgreedy,
        jointCfg
    )
    if (rt$done) {
        return(rt$result)
    }
    methods <- rt$methods
    .multiStudyPipelineDriver(
        data,
        rt$result,
        .fmPerStudy,
        .fmSumStats,
        list(
            methods = methods,
            contexts = contexts,
            traitId = traitId,
            region = region,
            cisWindow = cisWindow,
            jointRegions = jointRegions,
            addSusieInf = addSusieInf,
            credibleSetArgs = credibleSetArgs,
            fineMappingResult = fineMappingResult,
            verbose = verbose,
            crossValidationArgs = crossValidationArgs,
            residualizationArgs = residualizationArgs,
            seed = seed,
            naAction = naAction,
            signalScreenArgs = signalScreenArgs,
            genotypeFilterArgs = genotypeFilterArgs,
            panelFilterArgs = panelFilterArgs,
            rssArgs = rssArgs,
            usePCA = usePCA,
            nPCs = nPCs,
            fitRetention = fitRetention,
            mrmashPrior = mrmashPrior,
            dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff
        ),
        .rbindFineMappingResult,
        QtlFineMappingResult,
        "fineMappingPipeline"
    )
}

#' @rdname fineMappingPipeline
#' @export
setMethod(
    "fineMappingPipeline",
    "MultiStudyQtlDataset",
    .fmPipelineMultiStudy
)

# =============================================================================
# QtlSumStats method
# =============================================================================

# The settings every QtlSumStats joint-phase function forwards to the shared
# joint dispatcher: the common summary-statistics record plus the scope and
# TWAS-weight fields only the joint engine uses.
# @noRd
.fmQssJointCfg <- function(
    common,
    contexts,
    traitId,
    mrmashPrior,
    dataDrivenPriorWeightsCutoff
) {
    c(
        common,
        list(
            contexts = contexts,
            traitId = traitId,
            mrmashPrior = mrmashPrior,
            dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff
        )
    )
}

# Run the joint engine for a QtlSumStats with the given spec + token set.
# Shared by the explicit-jointSpecification path and the auto-detected
# cross-context path; only the spec and token arguments differ between them.
# @noRd
.fmQssJointDispatch <- function(jointSpec, tokens, methodArgs, cfg) {
    .fmDispatchJointSpecsQtlSumStats(
        jointSpec,
        cfg$data,
        tokens,
        cfg$contexts,
        cfg$traitId,
        cfg$verbose,
        methodArgs = methodArgs,
        mrmashPrior = cfg$mrmashPrior,
        dataDrivenPriorWeightsCutoff = cfg$dataDrivenPriorWeightsCutoff,
        fineMappingResult = cfg$fineMappingResult,
        credibleSetArgs = cfg$credibleSetArgs,
        fitRetention = cfg$fitRetention,
        panelFilterArgs = cfg$panelFilterArgs %||% PanelFilterParam()
    )
}

# The joint-specification phase for the QtlSumStats path: dispatch whatever
# the spec covers, narrow the tokens the joint fits consumed, and report
# whether anything is left for the per-tuple loop. `done = TRUE` means the
# joint fits are the whole answer.
#
# Deliberately mirrors .fmQdsJointPhase and .twasQdsJointPhase -- same
# responsibilities, same done / result shape.
# @noRd
.fmQssJointPhase <- function(parsedJointSpec, norm, cfg) {
    if (length(parsedJointSpec) == 0L) {
        return(list(
            done = FALSE,
            result = NULL,
            tokens = norm$tokens,
            methodArgs = norm$methodArgs
        ))
    }
    jointResult <- .fmQssJointDispatch(
        parsedJointSpec,
        intersect(norm$tokens, "mvsusie"),
        norm$methodArgs,
        cfg
    )
    # The joint dispatch consumed the multivariate tokens; the univariate
    # phase runs on whatever is left.
    narrowed <- .fmNarrowAfterJoint(TRUE, norm$tokens, norm$methodArgs)
    if (length(narrowed$tokens) == 0L) {
        if (is.null(jointResult)) {
            abort(glue(
                "fineMappingPipeline(QtlSumStats): no joint fits produced. ",
                "Check that the jointSpecification scope intersects the ",
                "available data."
            ))
        }
        return(list(done = TRUE, result = jointResult))
    }
    list(
        done = FALSE,
        result = jointResult,
        tokens = narrowed$tokens,
        methodArgs = narrowed$methodArgs
    )
}

# --- fine-mapping method-argument constructors ------------------------------
#
# A token can reach two susieR/mvsusieR entry points depending on whether the
# input is individual-level or summary statistics, so each constructor checks
# against BOTH, taken live from the capability table rather than transcribed.

# The engine(s) a fine-mapping token forwards to, from the capability table.
# @noRd
.fmMethodCallees <- function(token) {
    info <- .fineMappingMethodCapabilities[[token]]
    unique(compact(list(info$individualImpl, info$sumstatImpl)))
}

# Shared body for the per-engine constructors: pecotmr keeps no defaults of
# its own here (they are seeded per run from `L` / `Lgreedy`), so everything
# the caller gives is checked against that ONE engine and carried through.
#
# One constructor per engine, not per token. A method whose two input paths
# run different engines gets one for each: checking the union instead would
# accept a susie_rss-only name on an individual-level run and then drop it,
# and mvsusie could not be checked at all, since mvsusie_rss takes `...` and
# a union containing it accepts anything.
# @noRd
.fmMethodOptions <- function(callee, label, engine, extra) {
    .newMethodOptions(
        callee,
        defaults = list(),
        extra = extra,
        label = label,
        engine = engine
    )
}

# The union constructor behind a token, for splicing a PLAIN LIST in the
# aggregator. A plain list declares no engine and the run's input class is
# not known yet, so it must stay acceptable for either path; the per-input
# check in .fmCheckMethodArgsForInput narrows it once the class is known.
# @noRd
#' @rdname FineMappingMethodsParam
#' @aliases FineMappingMethodsParam-class
#' @exportClass FineMappingMethodsParam
setClass("FineMappingMethodsParam", contains = "MethodsSelectionParam")

#' @title Which Fine-Mapping Methods To Run, And How
#' @description Selects the methods \code{\link{fineMappingPipeline}} runs and
#'   carries each one's engine arguments.
#'
#'   A method that runs on both individual-level and summary-statistic data
#'   reaches a different engine on each path, so for a
#'   \code{MultiStudyQtlDataset} carrying both there is no single set of
#'   arguments per method. Name such a method under
#'   \code{qtlDatasetMethods} or \code{qtlSumStatsMethods} to say which
#'   path its options are for; name it under \code{methods} when the path
#'   need not be stated --- a single-type run, or a method with nothing
#'   path-specific to configure.
#'
#'   Naming a method selects it. A method belongs in exactly one of the
#'   three slots.
#' @param methods Named list of per-method options whose input path need not
#'   be stated. Each entry is that method's \code{*Options()} record, or
#'   \code{list()} to run it with its defaults.
#' @param qtlDatasetMethods Named list of per-method options for the
#'   individual-level path.
#' @param qtlSumStatsMethods Named list of per-method options for the
#'   summary-statistics path.
#' @return A \code{FineMappingMethodsParam} object, a \code{\link{MethodParam}}.
#' @examples
#' FineMappingMethodsParam(methods = "susie")
#' FineMappingMethodsParam(methods = list(susie = SusieOptions(L = 20)))
#' FineMappingMethodsParam(
#'     qtlDatasetMethods = list(susie = SusieOptions(L = 20)),
#'     qtlSumStatsMethods = list(susie = SusieRssOptions(L = 20))
#' )
#' @export
FineMappingMethodsParam <- function(
    methods = NULL,
    qtlDatasetMethods = NULL,
    qtlSumStatsMethods = NULL
) {
    new(
        "FineMappingMethodsParam",
        methods = .methodsNormalizeSlot(
            methods,
            "methods",
            "FineMappingMethodsParam"
        ),
        qtlDatasetMethods = .methodsNormalizeSlot(
            qtlDatasetMethods,
            "qtlDatasetMethods",
            "FineMappingMethodsParam"
        ),
        qtlSumStatsMethods = .methodsNormalizeSlot(
            qtlSumStatsMethods,
            "qtlSumStatsMethods",
            "FineMappingMethodsParam"
        )
    )
}

# Every fine-mapping method that runs on one input path. GwasSumStats is a
# summary-statistics class, narrowed to the methods that accept one trait.
# @noRd
.fmMethodsFor <- function(inputKind) {
    caps <- .fineMappingMethodCapabilities
    keep(names(caps), .fmMethodRunsOn, caps = caps, inputKind = inputKind)
}

# @noRd
.fmMethodRunsOn <- function(token, caps, inputKind) {
    info <- caps[[token]]
    if (identical(inputKind, "QtlDataset")) {
        return(!is.null(info$individualImpl))
    }
    if (is.null(info$sumstatImpl)) {
        return(FALSE)
    }
    !identical(inputKind, "GwasSumStats") || isTRUE(info$gwasAllowed)
}

# One method's Options constructor for ONE input path, validating against
# just that path's callee.
#
# Deliberately not the union of both callees: validating against the union
# would accept a name valid on only one of them -- `susie` takes 41
# arguments on individual data and 49 on summary statistics, sharing 37, so
# the union admits 16 names that are wrong for whichever path actually runs.
# Once the caller knows its input class there is no reason to be that
# lenient.
# @noRd
.fmTokenPathConfig <- function(token, inputKind) {
    force(token)
    force(inputKind)
    callee <- .fmMethodCalleeFor(token, inputKind)
    if (is.null(callee)) {
        return(NULL)
    }
    function(...) {
        .newMethodOptions(
            callee,
            defaults = list(),
            extra = list(...),
            label = glue("method '{token}' on {inputKind}"),
            engine = token
        )
    }
}

# Every method's constructor for one input path, omitting the methods that
# do not run on it -- `fsusie` has no summary-statistics implementation,
# `ser` no individual one.
# @noRd
.fmMethodCtorsFor <- function(inputKind) {
    toks <- names(.fineMappingMethodCapabilities)
    compact(set_names(
        map(toks, .fmTokenPathConfig, inputKind = inputKind),
        toks
    ))
}

#' @title Arguments For A Fine-Mapping Engine
#' @description Options forwarded to one fine-mapping engine, checked
#'   against that engine's live formals.
#'
#'   There is one constructor per ENGINE, not per method. A method whose
#'   individual-level and summary-statistic paths run different functions
#'   has one for each --- \code{SusieOptions} for \code{susieR::susie} and
#'   \code{SusieRssOptions} for \code{susieR::susie_rss}, \code{MvsusieOptions}
#'   and \code{MvsusieRssOptions} likewise. Checking the union of the two
#'   instead would accept an \code{susie_rss}-only name on an
#'   individual-level run and then silently drop it, and would leave
#'   mvsusie unchecked entirely, since \code{mvsusie_rss} takes \code{...}
#'   and a union containing it accepts any name at all.
#'
#'   \code{SerOptions} and \code{FsusieOptions} have no sibling: SER is
#'   summary-statistics-only and fSuSiE individual-only.
#'
#'   Pass the result as an element of \code{\link{FineMappingMethodsParam}},
#'   or directly to the matching single-engine function. In the aggregator a
#'   method accepts either of its engines' constructors, since the input
#'   class is not known at that point. A plain list works there too and
#'   stays checked against both; the exact engine is checked once the run's
#'   input class is known.
#' @param ... Arguments for the engine, under its own names (\code{L},
#'   \code{coverage}, \code{max_iter}, ...).
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' SusieOptions(L = 5)
#' SusieRssOptions(maf_thresh = 0.01)
#' MvsusieOptions(prior_tol = 1e-9)
#' @name fineMappingMethodOptions
NULL

# Token -> the engine constructor for each of its input paths. A method with
# both an individual-level and a summary-statistic engine accepts either
# one's constructor, since the aggregator does not know the input class.
# @noRd
.fmTokenCtorsByPath <- function() {
    list(
        susie = list(SusieOptions, SusieRssOptions),
        susieInf = list(SusieInfOptions, SusieInfRssOptions),
        susieAsh = list(SusieAshOptions, SusieAshRssOptions),
        ser = list(SerOptions),
        mvsusie = list(MvsusieOptions, MvsusieRssOptions),
        fsusie = list(FsusieOptions)
    )
}

# Resolve the method tokens for a QtlSumStats run and run any EXPLICIT
# jointSpecification up front (mvsusie via the axis dispatcher), removing the
# joint methods from the per-tuple token set. `done` flags the case
# where an explicit spec consumed every method, so the caller returns the joint
# result directly (or errors when it produced nothing).
# @noRd
.fmQssResolveTokens <- function(jointSpecification, methods, L, Lgreedy, cfg) {
    parsedJointSpec <- parseJointSpecification(jointSpecification, cfg$data)
    norm <- .fmNormalizeMethods(
        methods,
        inputKind = "QtlSumStats",
        L = L,
        Lgreedy = Lgreedy
    )
    .fmCheckMethodCapabilities(norm$tokens, "QtlSumStats")
    .fmCheckMethodArgsForInput(norm$methodArgs, "QtlSumStats")
    .fmQssJointPhase(parsedJointSpec, norm, cfg)
}

# Resolve the study/context/trait columns and the selected row indices for a
# QtlSumStats run, applying the optional contexts / traitId filters.
# @noRd
.fmQssSelectRows <- function(
    data,
    contexts,
    traitId
) {
    studyCol <- as.character(data$study)
    contextCol <- as.character(data$context)
    traitCol <- as.character(data$trait)
    byContext <- if (is.null(contexts)) {
        seq_len(nrow(data))
    } else {
        which(is_in(contextCol, contexts))
    }
    selRows <- if (is.null(traitId)) {
        byContext
    } else {
        byContext[is_in(traitCol[byContext], traitId)]
    }
    if (length(selRows) == 0L) {
        msg <- glue(
            "fineMappingPipeline(QtlSumStats): no entries matched the ",
            "supplied contexts / traitId filters."
        )
        abort(msg)
    }
    list(
        studyCol = studyCol,
        contextCol = contextCol,
        traitCol = traitCol,
        selRows = selRows
    )
}

# Partition tokens into the univariate RSS family (everything that isn't
# multivariate) and mvsusie, validating that mvsusie has >= 2 contexts for at
# least one (study, trait) group.
# @noRd
.fmQssSplitTokens <- function(tokens, sel) {
    univTokens <- tokens[!is_in(tokens, c("mvsusie", "fsusie"))]
    mvTokens <- tokens[tokens == "mvsusie"]
    if (length(mvTokens) > 0L) {
        groupKey <- str_c(
            sel$studyCol[sel$selRows],
            sel$traitCol[sel$selRows],
            sep = "||"
        )
        perGroupNCtx <- map_int(
            split(sel$contextCol[sel$selRows], groupKey),
            length
        )
        if (all(perGroupNCtx < 2L)) {
            msg <- glue(
                "fineMappingPipeline(QtlSumStats): mvsusie requires at ",
                "least two contexts per (study, trait); the supplied ",
                "collection has only one context per trait."
            )
            abort(msg)
        }
    }
    list(univTokens = univTokens, mvTokens = mvTokens)
}

# Multivariate dispatch via the joint engine (auto-detected cross-context RSS
# joint) for mvsusie WITHOUT an explicit jointSpecification, merged with any
# explicit-spec result the joint phase already produced.
# @noRd
.fmQssAutoJoint <- function(cfg, mvTokens, methodArgs, jointResult, ldSketch) {
    if (length(mvTokens) == 0L) {
        return(jointResult)
    }
    autoJoint <- .fmQssJointDispatch(
        list(list(axes = "context", scope = NULL)),
        mvTokens,
        methodArgs,
        cfg
    )
    if (is.null(jointResult)) {
        autoJoint
    } else if (is.null(autoJoint)) {
        jointResult
    } else {
        .rbindFineMappingResult(jointResult, autoJoint, ldSketch = ldSketch)
    }
}

# Assemble the QtlSumStats result: build the per-tuple QtlFineMappingResult
# from the collected row-records (region = entry variant span, no cis-window)
# and combine with the joint result. An all-screened collection yields a valid
# empty result rather than an error.
# @noRd
.fmQssAssemble <- function(
    data,
    ldSketch,
    rows,
    nSkipped,
    jointResult
) {
    rowContext <- map_chr(rows, "context")
    rowTrait <- map_chr(rows, "trait")
    perTupleResult <- if (length(rows) > 0L) {
        .fmBuildQtlResult(
            map_chr(rows, "study"),
            rowContext,
            rowTrait,
            map_chr(rows, "method"),
            map(rows, "entry"),
            traitPos = try_fetch(
                .anchorVector(data, rowContext, rowTrait, "traitPos"),
                error = function(cnd) NULL
            ),
            ldSketch = ldSketch
        )
    } else {
        NULL
    }
    if (is.null(jointResult)) {
        if (!is.null(perTupleResult)) {
            return(perTupleResult)
        }
        if (nSkipped > 0L) {
            return(.fmBuildQtlResult(
                character(0),
                character(0),
                character(0),
                character(0),
                list(),
                ldSketch = ldSketch,
                allowEmpty = TRUE
            ))
        }
        abort(
            "fineMappingPipeline(QtlSumStats): no entries produced a result."
        )
    }
    if (is.null(perTupleResult)) {
        return(jointResult)
    }
    .rbindFineMappingResult(perTupleResult, jointResult, ldSketch = ldSketch)
}

# The settings the joint phase and the RSS row builders share, named once
# so the two records built from them cannot disagree.
# @noRd
.fmQssConfigs <- function(
    data,
    credibleSetArgs,
    fineMappingResult,
    fitRetention,
    verbose,
    panelFilterArgs,
    contexts,
    traitId,
    mrmashPrior,
    dataDrivenPriorWeightsCutoff
) {
    common <- .fmSsCommonCfg(
        data = data,
        credibleSetArgs = credibleSetArgs,
        fineMappingResult = fineMappingResult,
        fitRetention = fitRetention,
        verbose = verbose,
        panelFilterArgs = panelFilterArgs
    )
    list(
        common = common,
        joint = .fmQssJointCfg(
            common,
            contexts = contexts,
            traitId = traitId,
            mrmashPrior = mrmashPrior,
            dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff
        )
    )
}

# QtlSumStats fine-mapping worker. The resolved tokens, row selection and LD
# sketch are ordinary locals, handed to the per-entry dispatch helpers as
# named arguments.
# @noRd
.fmPipelineQtlSumStats <- function(
    data,
    methods,
    addSusieInf,
    contexts,
    traitId,
    jointSpecification,
    credibleSetArgs,
    fineMappingResult,
    mrmashPrior,
    dataDrivenPriorWeightsCutoff,
    verbose,
    fitRetention,
    rssArgs,
    panelFilterArgs
) {
    .fmAssertQcd(data)
    cfgs <- .fmQssConfigs(
        data = data,
        credibleSetArgs = credibleSetArgs,
        fineMappingResult = fineMappingResult,
        fitRetention = fitRetention,
        verbose = verbose,
        panelFilterArgs = panelFilterArgs,
        contexts = contexts,
        traitId = traitId,
        mrmashPrior = mrmashPrior,
        dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff
    )
    commonCfg <- cfgs$common
    jointCfg <- cfgs$joint
    rt <- .fmQssResolveTokens(
        jointSpecification,
        methods,
        credibleSetArgs$L %||% 10L,
        credibleSetArgs$Lgreedy,
        jointCfg
    )
    if (rt$done) {
        return(rt$result)
    }
    .fmQssAfterJoint(
        data,
        rt = rt,
        commonCfg = commonCfg,
        jointCfg = jointCfg,
        contexts = contexts,
        traitId = traitId,
        addSusieInf = addSusieInf,
        rssArgs = rssArgs
    )
}

# Everything after the joint phase has had its say: select the rows the
# univariate tokens run on, resolve the LD-finiteness setting against the
# sketch, and fit. `rt` is .fmQssResolveTokens()'s record -- it carries the
# remaining tokens, their options, and whatever the joint phase produced.
# @noRd
.fmQssAfterJoint <- function(
    data,
    rt,
    commonCfg,
    jointCfg,
    contexts,
    traitId,
    addSusieInf,
    rssArgs
) {
    sel <- .fmQssSelectRows(data = data, contexts = contexts, traitId = traitId)
    ldSketch <- getLdSketch(data)
    cfg <- .fmRssRunCfg(
        commonCfg,
        ldSketch = ldSketch,
        addSusieInf = addSusieInf,
        methodArgs = rt$methodArgs,
        rssArgs = .fmRssArgsResolved(rssArgs, ldSketch)
    )
    .fmQssRunRows(
        data,
        sel = sel,
        split = .fmQssSplitTokens(rt$tokens, sel),
        cfg = cfg,
        jointCfg = jointCfg,
        ldSketch = ldSketch,
        methodArgs = rt$methodArgs,
        jointResult = rt$result
    )
}

# Fine-map the selected rows with the univariate tokens, then assemble that
# with whatever the joint phase produced. `sel` and `split` travel whole:
# the row selection and the token split are each one answer, and splitting
# them back into loose columns here is what this function exists to avoid.
# @noRd
.fmQssRunRows <- function(
    data,
    sel,
    split,
    cfg,
    jointCfg,
    ldSketch,
    methodArgs,
    jointResult
) {
    univOut <- if (length(split$univTokens) > 0L) {
        map(
            sel$selRows,
            .fmRssEntryRows,
            studyCol = sel$studyCol,
            contextCol = sel$contextCol,
            traitCol = sel$traitCol,
            univTokens = split$univTokens,
            cfg = cfg
        )
    } else {
        list()
    }
    .fmQssAssemble(
        data = data,
        ldSketch = ldSketch,
        rows = list_flatten(map(univOut, "rows")),
        nSkipped = sum(map_lgl(univOut, "skipped")),
        jointResult = .fmQssAutoJoint(
            cfg = jointCfg,
            mvTokens = split$mvTokens,
            methodArgs = methodArgs,
            jointResult = jointResult,
            ldSketch = ldSketch
        )
    )
}

#' @rdname fineMappingPipeline
#' @export
setMethod(
    "fineMappingPipeline",
    "QtlSumStats",
    function(
        data,
        methods,
        contexts = NULL,
        traitId = NULL,
        jointSpecification = NULL,
        addSusieInf = TRUE,
        credibleSetArgs = CredibleSetParam(),
        fineMappingResult = NULL,
        mrmashPrior = NULL,
        dataDrivenPriorWeightsCutoff = 1e-10,
        verbose = 1,
        fitRetention = c("slim", "full"),
        rssArgs = SusieRssParam(),
        panelFilterArgs = PanelFilterParam(),
        crossValidationArgs = CrossValidationParam(),
        residualizationArgs = ResidualizationParam()
    ) {
        fitRetention <- arg_match(fitRetention)
        .cvRefuseOnSumstats(
            crossValidationArgs,
            "fineMappingPipeline",
            "QtlSumStats"
        )
        .fmPipelineQtlSumStats(
            data = data,
            methods = methods,
            addSusieInf = addSusieInf,
            contexts = contexts,
            traitId = traitId,
            jointSpecification = jointSpecification,
            credibleSetArgs = credibleSetArgs,
            fineMappingResult = fineMappingResult,
            mrmashPrior = mrmashPrior,
            dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
            verbose = verbose,
            fitRetention = fitRetention,
            rssArgs = rssArgs,
            panelFilterArgs = panelFilterArgs
        )
    }
)

# =============================================================================
# GwasSumStats method
# =============================================================================

# Default the finite-sample LD reference size when an EB / SER-fallback LD-
# mismatch mode is active (serFallback on, or rMismatch other than "none") and
# rFinite is unset: use the LD-panel sample size (the notebook's `B`). Shared by
# the QtlSumStats and GwasSumStats workers.
# @noRd
# The resolved record, still a SusieRssParam. Previously the resolved value was
# carried alongside the record as a second cfg field and spliced in with
# list_assign(as.list(...)) at the point of use, which handed everything
# downstream a bare list exactly where the value became concrete.
# @noRd
.fmRssArgsResolved <- function(rssArgs, ldSketch) {
    resolved <- .fmResolveRFinite(
        rssArgs$rFinite,
        rssArgs$serFallback,
        rssArgs$rMismatch,
        ldSketch
    )
    if (identical(resolved, rssArgs$rFinite)) {
        return(rssArgs)
    }
    .rssParamWithRFinite(rssArgs, resolved)
}

.fmResolveRFinite <- function(rFinite, serFallback, rMismatch, ldSketch) {
    if (
        is.null(rFinite) &&
            (isTRUE(serFallback) || !identical(rMismatch, "none"))
    ) {
        .ldSketchNSamples(ldSketch)
    } else {
        rFinite
    }
}

# GwasSumStats fine-mapping worker. The resolved tokens, LD sketch and
# finite-sample size are ordinary locals, handed to the per-entry dispatch
# helpers as named arguments. One GwasSumStats is one LD block (the caller
# builds one collection per block when sweeping the genome); we fine-map each
# (study, method) tuple across the whole entry, with no in-pipeline block
# partitioning.
# @noRd
.fmPipelineGwas <- function(
    data,
    methods,
    addSusieInf,
    credibleSetArgs,
    fineMappingResult,
    verbose,
    fitRetention,
    rssArgs,
    panelFilterArgs
) {
    .fmAssertQcd(data)
    norm <- .fmNormalizeMethods(
        methods,
        L = credibleSetArgs$L %||% 10L,
        Lgreedy = credibleSetArgs$Lgreedy
    )
    .fmCheckMethodCapabilities(norm$tokens, "GwasSumStats")
    .fmCheckMethodArgsForInput(norm$methodArgs, "GwasSumStats")
    ldSketch <- getLdSketch(data)
    # Derived values are ordinary locals now, not fields grafted onto a bundle.
    tokens <- norm$tokens
    methodArgs <- norm$methodArgs
    rssArgs <- .fmRssArgsResolved(rssArgs, ldSketch)
    studyCol <- as.character(data$study)
    commonCfg <- .fmSsCommonCfg(
        data = data,
        credibleSetArgs = credibleSetArgs,
        fineMappingResult = fineMappingResult,
        fitRetention = fitRetention,
        verbose = verbose,
        panelFilterArgs = panelFilterArgs
    )
    .fmGwasRunEntries(
        data,
        studyCol = studyCol,
        tokens = tokens,
        commonCfg = commonCfg,
        ldSketch = ldSketch,
        addSusieInf = addSusieInf,
        methodArgs = methodArgs,
        rssArgs = rssArgs
    )
}

# Fine-map every row of the GWAS collection and assemble the result.
# @noRd
.fmGwasRunEntries <- function(
    data,
    studyCol,
    tokens,
    commonCfg,
    ldSketch,
    addSusieInf,
    methodArgs,
    rssArgs
) {
    entryOut <- map(
        seq_len(nrow(data)),
        .fmGwasEntryRows,
        studyCol = studyCol,
        tokens = tokens,
        cfg = .fmRssRunCfg(
            commonCfg,
            ldSketch = ldSketch,
            addSusieInf = addSusieInf,
            methodArgs = methodArgs,
            rssArgs = rssArgs
        )
    )
    rows <- list_flatten(map(entryOut, "rows"))
    nSkipped <- sum(map_lgl(entryOut, "skipped"))
    # An all-screened (or empty-input) collection legitimately yields a 0-row
    # result -- allow it instead of erroring "no ... tuples produced a result".
    .fmBuildGwasResult(
        map_chr(rows, "study"),
        map_chr(rows, "method"),
        map(rows, "entry"),
        blockIds = map_chr(rows, "blockId"),
        ldSketch = ldSketch,
        allowEmpty = (nSkipped > 0L || nrow(data) == 0L)
    )
}

#' @rdname fineMappingPipeline
#' @export
setMethod(
    "fineMappingPipeline",
    "GwasSumStats",
    function(
        data,
        methods,
        addSusieInf = TRUE,
        credibleSetArgs = CredibleSetParam(),
        fineMappingResult = NULL,
        verbose = 1,
        fitRetention = c("slim", "full"),
        rssArgs = SusieRssParam(),
        panelFilterArgs = PanelFilterParam(),
        crossValidationArgs = CrossValidationParam(),
        residualizationArgs = ResidualizationParam()
    ) {
        fitRetention <- arg_match(fitRetention)
        .cvRefuseOnSumstats(
            crossValidationArgs,
            "fineMappingPipeline",
            "GwasSumStats"
        )
        .fmPipelineGwas(
            data = data,
            methods = methods,
            addSusieInf = addSusieInf,
            credibleSetArgs = credibleSetArgs,
            fineMappingResult = fineMappingResult,
            verbose = verbose,
            fitRetention = fitRetention,
            rssArgs = rssArgs,
            panelFilterArgs = panelFilterArgs
        )
    }
)

# =============================================================================
# Named helpers for map/apply call sites (no inline lambdas)
# =============================================================================

# @noRd
.fmIssueDetail <- function(x) {
    glue("{x$token} {x$reason}")
}

# @noRd
.fmHasCached <- function(l) {
    !is.null(l$cached)
}

# @noRd
.fmNotCached <- function(l) {
    is.null(l$cached)
}

# @noRd
.fmGwasLookup <- function(tk, fineMappingResult, st, blockId) {
    list(
        tk = tk,
        cached = .fmCacheLookupGwasResume(
            fineMappingResult,
            st,
            tk,
            blockId
        )
    )
}

# @noRd
.fmGwasRowFor <- function(tk, st, blockId, ents) {
    .fmGwasRow(st, tk, blockId, ents[[tk]])
}

# @noRd
.fmGwasRowFromLookup <- function(l, st, blockId) {
    .fmGwasRow(st, l$tk, blockId, l$cached)
}

# @noRd
.fmQtlLookup <- function(tk, fineMappingResult, st, ctx, tr) {
    list(
        tk = tk,
        cached = .fmCacheLookup(fineMappingResult, st, ctx, tr, tk)
    )
}

# @noRd
.fmQtlRowFor <- function(tk, st, ctx, tr, ents) {
    .fmQtlRow(st, ctx, tr, tk, ents[[tk]])
}

# @noRd
.fmQtlRowFromLookup <- function(l, st, ctx, tr) {
    .fmQtlRow(st, ctx, tr, l$tk, l$cached)
}

# The i-th element of a region set (kept as a length-1 region).
# @noRd
.fmNthRegion <- function(i, region) {
    region[i]
}

# SER pre-screen for the j-th response column.
# @noRd
.fmSerScreenColumn <- function(j, X, Y, screen) {
    .fmSerScreen(X, Y[, j], screen)
}

# Re-label the j-th credible-set entry, preserving the "_0" sentinel.
# @noRd
.fmRelabelCsOne <- function(j, csVec, parts, offset) {
    if (is.na(parts[j, 1L])) {
        return(csVec[[j]])
    }
    idx <- as.integer(parts[j, 3L])
    if (idx == 0L) csVec[[j]] else str_c(parts[j, 2L], "_", idx + offset)
}

# FineMappingRow slot accessors (S4 slots can't be plucked by name).
# @noRd
.fmEntryVariantIds <- function(e) {
    .fmrPartsVariantIds(e)
}

# @noRd
.fmEntryTopLoci <- function(e) {
    .fmrPartsTopLoci(e)
}

# @noRd
.fmEntrySusieFit <- function(e) {
    .fmrPartsSusieFit(e)
}

# @noRd
.fmEntryCvResult <- function(e) {
    .fmrPartsCvResult(e)
}

# One merged row-record for token `tk` across window block entries (NULL when
# any window failed to fit the token).
# @noRd
.fmMergeTokenRow <- function(tk, study, ctx, trait, blockEntries) {
    ents <- map(blockEntries, .fmBlockToken, tk = tk)
    if (any(map_lgl(ents, is.null))) {
        return(NULL)
    }
    entry <- if (length(ents) == 1L) ents[[1L]] else .fmMergeEntries(ents)
    .fmQtlRow(study, ctx, trait, tk, entry)
}

# @noRd
.fmBlockToken <- function(be, tk) {
    be[[tk]]
}

# @noRd
.fmUnivLookup <- function(tk, fineMappingResult, study, ctx, tid) {
    list(
        tk = tk,
        cached = .fmCacheLookup(fineMappingResult, study, ctx, tid, tk)
    )
}

# @noRd
.fmUnivCachedRow <- function(l, study, ctx, tid) {
    .fmQtlRow(study, ctx, tid, l$tk, l$cached)
}

# All univariate row-records for one context (over its per-context traits).
# @noRd
.fmUnivContextRows <- function(ctx, cfg) {
    list_flatten(map(
        cfg$perCtxTraits[[ctx]],
        .fmUnivTraitRows,
        ctx = ctx,
        cfg = cfg
    ))
}

# =============================================================================
# Fine-mapping Param bundles
# -----------------------------------------------------------------------------
# Moved here from the former pipelineParams.R: these describe fine-mapping
# stages, and are read here and in fineMappingWrappers.R / jointEngine.R even
# when another pipeline forwards them.
#
# `L` / `Lgreedy` are CredibleSetParam fields: they bound how many credible
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

# --- CredibleSetParam accessors --------------------------------------------
# Setters revalidate: slot assignment checks the declared type on its own but
# does not re-run a validity method, and this package does write those (see
# SignalScreenParam).

#' @rdname CredibleSetParam
setMethod("getCoverage", "CredibleSetParam", function(x) x@coverage)

#' @rdname CredibleSetParam
setMethod("setCoverage", "CredibleSetParam", function(x, value) {
    x@coverage <- value
    validObject(x)
    x
})

#' @rdname CredibleSetParam
setMethod("getSecondaryCoverage", "CredibleSetParam", function(x) {
    x@secondaryCoverage
})

#' @rdname CredibleSetParam
setMethod("setSecondaryCoverage", "CredibleSetParam", function(x, value) {
    x@secondaryCoverage <- value
    validObject(x)
    x
})

#' @rdname CredibleSetParam
setMethod("getSignalCutoff", "CredibleSetParam", function(x) x@signalCutoff)

#' @rdname CredibleSetParam
setMethod("setSignalCutoff", "CredibleSetParam", function(x, value) {
    x@signalCutoff <- value
    validObject(x)
    x
})

#' @rdname CredibleSetParam
setMethod("getMinAbsCorr", "CredibleSetParam", function(x) x@minAbsCorr)

#' @rdname CredibleSetParam
setMethod("setMinAbsCorr", "CredibleSetParam", function(x, value) {
    x@minAbsCorr <- value
    validObject(x)
    x
})

#' @rdname CredibleSetParam
setMethod("getMedianAbsCorr", "CredibleSetParam", function(x) x@medianAbsCorr)

#' @rdname CredibleSetParam
setMethod("setMedianAbsCorr", "CredibleSetParam", function(x, value) {
    x@medianAbsCorr <- value
    validObject(x)
    x
})

#' @rdname CredibleSetParam
setMethod("getIncludeAllCs", "CredibleSetParam", function(x) x@includeAllCs)

#' @rdname CredibleSetParam
setMethod("setIncludeAllCs", "CredibleSetParam", function(x, value) {
    x@includeAllCs <- value
    validObject(x)
    x
})

#' @rdname CredibleSetParam
setMethod("getPerCsColumns", "CredibleSetParam", function(x) x@perCsColumns)

#' @rdname CredibleSetParam
setMethod("setPerCsColumns", "CredibleSetParam", function(x, value) {
    x@perCsColumns <- value
    validObject(x)
    x
})

#' @rdname CredibleSetParam
setMethod("getL", "CredibleSetParam", function(x) x@L)

#' @rdname CredibleSetParam
setMethod("getLgreedy", "CredibleSetParam", function(x) x@Lgreedy)

#' @rdname SignalScreenParam
#' @aliases SignalScreenParam-class
#' @exportClass SignalScreenParam
setClass(
    "SignalScreenParam",
    contains = "MethodParam",
    slots = c(
        pip = "numeric_OR_NULL",
        absZ = "numeric_OR_NULL",
        bf = "numeric_OR_NULL",
        logBf = "numeric_OR_NULL"
    )
)

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
#' @return A \code{SignalScreenParam} object, a \code{\link{MethodParam}}.
#' @examples
#' SignalScreenParam(absZ = 5)
#' @export
SignalScreenParam <- function(
    pip = NULL,
    absZ = NULL,
    bf = NULL,
    logBf = NULL
) {
    new("SignalScreenParam", pip = pip, absZ = absZ, bf = bf, logBf = logBf)
}

# absZ and bf have no meaningful negative cutoff: absZ screens on max|Z| and
# a Bayes factor is positive. (pip keeps `< 0 => 3 / nVariants`, and logBf is
# a log, so both may legitimately be negative.)
# @noRd
.screenPositiveProblem <- function(screen) {
    bad <- keep(c("absZ", "bf"), .screenIsNegative, screen = screen)
    if (length(bad) == 0L) {
        return(character(0))
    }
    m <- bad[[1L]]
    glue(
        "`{m}` must be > 0, got {screen[[m]]}. ",
        "absZ screens on max|Z| and Bayes factors are positive."
    )
}

# @noRd
.screenIsNegative <- function(m, screen) {
    v <- screen[[m]]
    !is.null(v) && length(v) == 1L && !is.na(v) && v < 0
}

# Refuse more than one enabled metric. Enabled means set AND non-zero: 0 is
# the long-standing "off" spelling for these cutoffs, so SignalScreenParam(pip
# = 0, absZ = 5) is one screen, not two.
# @noRd
.screenOneMetricProblem <- function(screen) {
    on <- names(screen)[map_lgl(as.list(screen), .screenIsOn)]
    if (length(on) <= 1L) {
        return(character(0))
    }
    glue(
        "only one screening metric may be enabled at a time; ",
        "got {str_flatten(on, ', ')}."
    )
}

setValidity("SignalScreenParam", function(object) {
    problems <- c(
        .screenOneMetricProblem(object),
        .screenPositiveProblem(object)
    )
    if (length(problems) == 0L) TRUE else as.character(problems)
})

# @noRd
.screenIsOn <- function(x) {
    !is.null(x) && length(x) > 0L && any(as.numeric(x) != 0, na.rm = TRUE)
}

# The one enabled metric as the polymorphic screen spec the pipelines already
# understand (see .asScreen): a bare numeric for `pip`, a list(metric, cutoff)
# for the others, NULL when nothing is enabled. A field set to 0 is off, so
# SignalScreenParam(pip = 0, absZ = 5) resolves to the absZ screen.
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

#' @rdname CredibleSetParam
#' @aliases CredibleSetParam-class
#' @exportClass CredibleSetParam
setClass(
    "CredibleSetParam",
    contains = "MethodParam",
    slots = c(
        coverage = "numeric",
        secondaryCoverage = "numeric",
        signalCutoff = "numeric",
        minAbsCorr = "numeric",
        medianAbsCorr = "numeric_OR_NULL",
        includeAllCs = "logical",
        perCsColumns = "character",
        L = "numeric",
        Lgreedy = "numeric_OR_NULL"
    )
)

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
#'   did not set it through \code{\link{FineMappingMethodsParam}}.
#' @param Lgreedy Integer or \code{NULL}. Maximum number of single effects
#'   for the greedy initialization stage, where the engine has one.
#'   \code{NULL} (the default) leaves the engine's own default in place.
#' @return \code{CredibleSetParam} returns a \code{CredibleSetParam}
#'   object, a \code{\link{MethodParam}}. Each \code{get*} returns that
#'   setting's value; each \code{set*} returns a modified copy. \code{L}
#'   and \code{Lgreedy} have no setter: they are seeded onto the fitting
#'   tokens, so changing one needs a refit rather than a re-summary.
#' @examples
#' CredibleSetParam(coverage = 0.9, includeAllCs = TRUE)
#' @export
CredibleSetParam <- function(
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
    new(
        "CredibleSetParam",
        coverage = coverage,
        secondaryCoverage = secondaryCoverage,
        signalCutoff = signalCutoff,
        minAbsCorr = minAbsCorr,
        medianAbsCorr = medianAbsCorr,
        includeAllCs = includeAllCs,
        perCsColumns = perCsColumns,
        L = L,
        Lgreedy = Lgreedy
    )
}

#' @rdname SusieRssParam
#' @aliases SusieRssParam-class
#' @exportClass SusieRssParam
setClass(
    "SusieRssParam",
    contains = "MethodParam",
    slots = c(
        serFallback = "logical",
        keepFullFit = "character",
        rFinite = "numeric_OR_NULL",
        rMismatch = "character_OR_NULL",
        control = "MethodOptions_OR_NULL"
    )
)

#' @title SuSiE-RSS Solver Settings
#' @description How the SuSiE-RSS fit is conditioned and what it falls back
#'   to, on \code{QtlSumStats} / \code{GwasSumStats} input.
#'
#'   SuSiE-family only, which is what the name says: every field is a
#'   \code{susieR} concept --- \code{serFallback} drops to
#'   \code{susie_ser}, \code{rFinite} / \code{rMismatch} are
#'   \code{susie_rss}'s \code{R_finite} / \code{R_mismatch}, and
#'   \code{control} is \code{susie_rss_control}. The other
#'   summary-statistics engine, \code{mvsusieR::mvsusie_rss}, has none of
#'   these formals and is not configured by this.
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
#'   \code{\link{SusieRssControlOptions}}; \code{NULL} (default) leaves that
#'   function's own defaults in place.
#' @return A \code{SusieRssParam} object, a \code{\link{MethodParam}}.
#' @examples
#' SusieRssParam(serFallback = TRUE, keepFullFit = "all")
#' @export
SusieRssParam <- function(
    serFallback = FALSE,
    keepFullFit = c("fallback", "all", "none"),
    rFinite = NULL,
    rMismatch = "none",
    control = NULL
) {
    keepFullFit <- arg_match(keepFullFit)
    if (!is.null(control)) {
        .assertMethodOptions(control, "SusieRssControlOptions", "control")
    }
    new(
        "SusieRssParam",
        serFallback = serFallback,
        keepFullFit = keepFullFit,
        rFinite = rFinite,
        rMismatch = rMismatch,
        control = control
    )
}

# `rFinite` is the one field pecotmr itself fills in: it defaults to NULL
# meaning "use the LD panel's sample size", which is not knowable when the
# caller builds the record. Resolution happens once the sketch is in hand
# (see .fmRssArgsResolved). Internal, not an exported setter -- no user sets
# this, it is derived.
# @noRd
.rssParamWithRFinite <- function(x, value) {
    x@rFinite <- value
    validObject(x)
    x
}
