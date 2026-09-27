#' @include qtlSumStats.R gwasSumStats.R
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
#'   \code{susie} and/or \code{susieAsh} and \code{addSusieInf = TRUE}
#'   (default), the SuSiE-inf fit is computed first and used as initialisation
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
#'   \item Filtering knobs (\code{mafCutoff}, \code{imissCutoff},
#'         \code{xvarCutoff}, \code{ldReferenceMetaFile}): individual-
#'         level QC lives on the \code{QtlDataset} constructor; sumstat
#'         QC lives in \code{summaryStatsQc()}. No filtering happens
#'         inside this pipeline.
#'   \item Diagnostic re-analysis paths
#'         (\code{singleEffect} / \code{bayesianConditionalRegression}
#'         reanalysis on the RSS path): these are not exposed as
#'         dedicated method tokens. Callers who want a single-effect
#'         fit can request it via per-method kwargs, e.g.
#'         \code{methods = list(susie = list(L = 1))} (see the
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
#'     \item A named list keyed by method token, where each value is a
#'           list of per-method kwargs to splice into the underlying
#'           fitter, e.g.
#'           \code{list(susie = list(L = 1, refine = FALSE),
#'                      mvsusie = list(max_iter = 500))}. Mirrors the
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
#'   alongside \code{susie} and/or \code{susieAsh}, controls whether the
#'   SuSiE-inf fit initialises the chained downstream method(s). Default
#'   \code{TRUE}.
#' @param coverage Primary credible-set coverage (numeric, length 1). Default
#'   \code{0.95}.
#' @param secondaryCoverage Secondary coverages forwarded to
#'   \code{postprocessFinemappingFits}. Default \code{c(0.7, 0.5)}.
#' @param signalCutoff PIP cutoff for top-loci selection. Default \code{0.025}.
#' @param minAbsCorr Minimum absolute correlation for credible-set purity.
#'   Default \code{0.8}.
#' @param medianAbsCorr Optional median absolute correlation for credible-set
#'   purity, routed to \code{susieR::susie_get_cs}. A set is kept if it passes
#'   either \code{minAbsCorr} or \code{medianAbsCorr} (OR-logic). Default
#'   \code{NULL} (off).
#' @param fineMappingResult Optional existing \code{FineMappingResult} to use as
#'   a resume cache; tuples already present are not refit.
#' @param cvFolds Integer. Number of cross-validation folds. Default \code{0}
#'   (no CV). When \code{> 1}, each method is refit on the training samples of
#'   every fold and used to predict the held-out samples; the fold partition
#'   plus per-fold out-of-fold predictions and metrics are stored on each
#'   \code{FineMappingRow}'s \code{cvResult} slot (see
#'   \code{\link{getCvResult}}). \code{twasWeightsPipeline} reuses this
#'   partition and feeds these predictions into the SR-TWAS ensemble.
#'   Individual-level (\code{QtlDataset} / \code{MultiStudyQtlDataset}) input
#'   only; ignored for sumstat inputs.
#' @param cvThreads Integer. Number of parallel workers for the cross-validation
#'   fold refits (passed to the shared CV engine's \code{numThreads}). Default
#'   \code{1} (serial); \code{-1} uses all cores. Only consulted when
#'   \code{cvFolds > 1}.
#' @param samplePartition Optional pre-defined CV partition \code{data.frame}
#'   with columns \code{Sample} and \code{Fold}. When supplied (and
#'   \code{cvFolds > 1}), every method reuses this exact partition; otherwise a
#'   fresh partition is generated per \code{(study, context, trait)}.
#' @param seed Optional integer. When non-NULL,
#'   \code{withr::local_seed(seed)} is called at the start of the call for
#'   reproducible fits; the global RNG state is restored on return. Default
#'   \code{NULL}
#'   (no seeding).
#' @param pipCutoffToSkip Numeric (length 1). Individual-level single-effect
#'   (SER) pre-screen applied to each residualized \code{(X, y)} block before a
#'   full fit: a susie model with \code{L = 1} is fit and the block is skipped
#'   when no PIP exceeds the cutoff (no potentially significant variant). The
#'   summary-statistics analog lives in \code{summaryStatsQc()}. \code{0}
#'   (default) disables the screen; a negative value uses the adaptive \code{3 /
#'   nVariants} threshold.
#' @param absZCutoffToSkip,bfCutoffToSkip,logBfCutoffToSkip Numeric (length 1).
#'   Alternative individual-level pre-screen metrics, in place of
#'   \code{pipCutoffToSkip}: skip a block unless its maximum marginal \code{|z|}
#'   (\code{absZCutoffToSkip}), or its maximum per-variant single-effect Bayes
#'   factor (\code{bfCutoffToSkip}) / log Bayes factor
#'   (\code{logBfCutoffToSkip}) from the \code{L = 1} fit, exceeds the cutoff.
#'   Each defaults to 0 (off). Exactly one of the four \code{*CutoffToSkip}
#'   arguments may be non-zero (one screening metric at a time).
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
#'   \code{mafCutoff} / \code{macCutoff} / \code{imissCutoff} are measured
#'   against the \strong{LD reference panel} instead: a variant whose panel
#'   genotypes fall below the cutoffs is dropped before the z-scores and LD
#'   matrix are built. The thresholds mean the same thing as on the
#'   \code{QtlDataset} path (MAC is converted to a MAF equivalent and the
#'   stricter of the two applies), so one number carries across input types.
#'   Defaults (\code{0}, \code{0}, \code{1}) filter nothing.
#'
#'   This discards \emph{observed} variants, unlike
#'   \code{summaryStatsQc(imputeOpts = ...)}, which only bounds which variants
#'   RAISS will impute. Use it when the panel cannot support the LD estimate a
#'   rare variant would need.
#'
#' @param mafCutoff Numeric or \code{NULL}. Per-call override of the
#'   \code{QtlDataset} minor-allele-frequency filter; \code{NULL} uses the
#'   dataset's stored value.
#' @param macCutoff Numeric or \code{NULL}. Per-call override of the
#'   \code{QtlDataset} minor-allele-count filter; \code{NULL} uses the dataset's
#'   stored value.
#' @param xvarCutoff Numeric or \code{NULL}. Per-call override of the
#'   \code{QtlDataset} per-variant variance filter; \code{NULL} uses the
#'   dataset's stored value.
#' @param imissCutoff Numeric or \code{NULL}. Per-call override of the
#'   \code{QtlDataset} per-sample missingness filter; \code{NULL} uses the
#'   dataset's stored value.
#' @param keepIndel Logical or \code{NULL}. Per-call override of the
#'   \code{QtlDataset} indel-retention flag; \code{NULL} uses the dataset's
#'   stored value.
#' @param keepSamples Character vector or \code{NULL}. Per-call override of the
#'   \code{QtlDataset} sample allow-list; \code{NULL} uses the dataset's stored
#'   value.
#' @param keepVariants Character vector or \code{NULL}. Per-call override of the
#'   \code{QtlDataset} variant allow-list; \code{NULL} uses the dataset's stored
#'   value.
#' @param L Integer. Maximum number of SuSiE single effects. Default \code{10}.
#' @param Lgreedy Integer or \code{NULL}. Number of greedily-added effects in
#'   the SuSiE-inf refinement (the greedy-L loop). \code{NULL} (default)
#'   disables the greedy loop and fits \code{L} directly.
#' @param twasWeights Optional \code{\link{TwasWeights}} resume cache to reuse
#'   previously fitted weights; \code{NULL} fits fresh.
#' @param dataDrivenPriorWeightsCutoff Numeric or \code{NULL}. Cutoff below
#'   which data-driven prior weights are zeroed; \code{NULL} disables the
#'   cutoff.
#' @param naAction Character. How to handle missing values in the extracted
#'   phenotype/genotype data.
#' @param verbose Verbosity (0 silent, 1 default). Default \code{1}.
#' @param phenotypeCovariatesToResidualize Character vector (or \code{NULL}) of
#'   phenotype-covariate names to residualize against. \code{NULL} (default)
#'   uses every available phenotype covariate. Only meaningful when the input is
#'   a \code{QtlDataset} / \code{MultiStudyQtlDataset} (ignored for sumstat
#'   inputs).
#' @param genotypeCovariatesToResidualize Character vector (or \code{NULL}) of
#'   genotype-covariate column names to residualize against. \code{NULL} uses
#'   every available genotype covariate.
#' @param residualizePhenotypeCovariates Logical (length 1). When \code{TRUE}
#'   (default) residualize against the phenotype-side covariates listed in
#'   \code{phenotypeCovariatesToResidualize}. Set \code{FALSE} to disable
#'   phenotype-covariate residualization entirely. The marginal univariate
#'   effects stored on each \code{FineMappingRow} obey the same
#'   residualization choice as the SuSiE fit itself -- they are computed against
#'   the same residualized \code{X} / \code{Y}.
#' @param residualizeGenotypeCovariates Logical (length 1). When \code{TRUE}
#'   (default) residualize against the genotype-side covariates listed in
#'   \code{genotypeCovariatesToResidualize}. Set \code{FALSE} to disable.
#' @param trim Logical (length 1). When \code{TRUE} (default) the
#'   \code{susieFit} slot on each output \code{FineMappingRow} carries a
#'   trimmed view of the SuSiE fit (the minimal subset needed by downstream
#'   pipelines). When \code{FALSE} the full untrimmed \code{susie()} return is
#'   retained so accessors like \code{getSusieFit()} and non-default-coverage
#'   queries through \code{getCs()} can read the full posterior matrices
#'   (\code{lbf_variable}, \code{mu}, \code{mu2}, \code{V}). The per-variant
#'   \code{topLoci} table is always fully populated regardless of \code{trim}.
#' @param fullFit Logical (length 1, default \code{FALSE}). Master switch for
#'   per-credible-set variant-level export in \code{topLoci}. When \code{FALSE},
#'   only the always-on \code{within_cs_pip} scalar column is added (each
#'   variant's \code{alpha} in its assigned primary-coverage credible set; NA
#'   when not in a set). When \code{TRUE}, wide per-CS columns are added, one
#'   set of columns per credible set (labelled \code{cs<k>}).
#' @param fullFitAlphaOnly Logical (length 1, default \code{TRUE}). No-op when
#'   \code{fullFit=FALSE}. When \code{TRUE}, only \code{alpha} is widened
#'   (\code{within_cs_pip_cs<k>}). When \code{FALSE}, all four per-effect
#'   matrices are widened: \code{within_cs_pip_cs<k>} (alpha),
#'   \code{cs_logbf_cs<k>} (lbf_variable), \code{cs_effect_cs<k>} (mu, unscaled)
#'   and \code{cs_effect_var_cs<k>} (posterior variance, unscaled).
#' @param includeAllCs Logical (length 1, default \code{FALSE}). No-op when
#'   \code{fullFit=FALSE}. When \code{FALSE}, only effects that produced a
#'   passing (purity/coverage-filtered) credible set are widened. When
#'   \code{TRUE}, every effect \code{L} is widened (including filtered-out ones,
#'   labelled \code{L<k>} instead of \code{cs<k>}).
#' @param serFallback Logical (length 1, default \code{FALSE}).
#'   \code{QtlSumStats} / \code{GwasSumStats} (summary-statistics SuSiE-RSS
#'   paths) only. When \code{TRUE}, after each standard multi-effect SuSiE-RSS
#'   fit (\code{susie} / \code{susieAsh}) the pipeline reads susieR's
#'   finite-sample R diagnostics and, if
#'   \code{fit$R_finite_diagnostics$R_reliability_flag} is \code{TRUE}, reports
#'   the single-effect (\code{ser_model}) result for that region instead of the
#'   multi-effect fit. Defaults \code{FALSE} so existing callers are unchanged.
#'   See \code{keepFullFit} for retaining the multi-effect fit.
#' @param rFinite \code{QtlSumStats} / \code{GwasSumStats} only. Finite-sample
#'   size for susieR's \code{R_finite} correction, forwarded to
#'   \code{susieR::susie_rss()}. \code{NULL} (default) uses susieR's default,
#'   except when a finite/EB mode is active (\code{serFallback=TRUE} or
#'   \code{rMismatch != "none"}) and \code{rFinite} is \code{NULL}, in which
#'   case it defaults to the LD-panel sample size \code{getNSamples(ldSketch)}.
#' @param rMismatch \code{QtlSumStats} / \code{GwasSumStats} only. LD-mismatch
#'   correction mode forwarded to \code{susieR::susie_rss()} as
#'   \code{R_mismatch}: \code{"none"} (default), \code{"eb"} (empirical Bayes),
#'   or \code{"eb_mix"} (residual-mixture EB).
#' @param rssControl \code{QtlSumStats} / \code{GwasSumStats} only. Optional
#'   named list of \code{susieR::susie_rss_control()} settings (e.g.
#'   \code{check_prior}, \code{mismatch_estimator}), forwarded as
#'   \code{susie_rss()}'s \code{control} argument. Default \code{NULL} leaves
#'   the \code{susie_rss_control()} defaults in place.
#' @param keepFullFit \code{QtlSumStats} / \code{GwasSumStats} only. Controls
#'   retention of the pre-fallback multi-effect SuSiE-RSS fit when
#'   \code{serFallback=TRUE}: \code{"fallback"} (default) keeps it only for
#'   regions that fell back to SER; \code{"all"} keeps it for every region;
#'   \code{"none"} keeps none. The retained fit and the decision are stored on
#'   the entry's SuSiE fit and read via \code{getSusieFit(res)$multiEffectFit},
#'   \code{getSusieFit(res)$R_reliability_flag}, and
#'   \code{getSusieFit(res)$serFallbackUsed}.
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
    standardGeneric("fineMappingPipeline")
})


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
# Accepts: * character vector c("susie", "susieInf") -> empty kwargs per token *
# named list list(susie = list(L = 1), ...) -> per-token kwargs
#
# Names of the returned `methodArgs` always equal `tokens` (one entry per
# token, empty list when the user supplied none). The fitters then
# `modifyList`-merge each entry into the base arg list before do.call.
#
# Mirrors the convention of .twasNormalizeMethods so the two pipelines
# expose the same shape on the user side.
# @noRd
.fmNormalizeMethods <- function(methods, L = 10L, Lgreedy = NULL) {
    if (is.null(methods) || length(methods) == 0L) {
        msg <- glue(
            "fineMappingPipeline: `methods` must be a non-empty character ",
            "vector or named list of <token> = <kwargs> entries."
        )
        abort(msg)
    }
    parsed <- if (is.character(methods)) {
        .fmMethodsFromTokens(methods)
    } else if (is.list(methods)) {
        .fmMethodsFromKwargs(methods)
    } else {
        cls <- class(methods)[[1L]]
        msg <- glue(
            "fineMappingPipeline: `methods` must be a character vector or ",
            "named list. Got class '{cls}'."
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

# A named list of <token> = <kwargs> entries.
# @noRd
.fmMethodsFromKwargs <- function(methods) {
    if (is.null(names(methods)) || any(names(methods) == "")) {
        msg <- glue(
            "fineMappingPipeline: when `methods` is a list it must be ",
            "named (one entry per method token)."
        )
        abort(msg)
    }
    nonListChild <- !map_lgl(methods, is.list)
    if (any(nonListChild)) {
        badNames <- str_flatten(names(methods)[nonListChild], ", ")
        msg <- glue(
            "fineMappingPipeline: each entry of the `methods` list must ",
            "itself be a list of named kwargs (got non-list value ",
            "for: {badNames})."
        )
        abort(msg)
    }
    tokens <- unique(names(methods))
    list(tokens = tokens, methodArgs = methods[tokens])
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
.fmGwasFitRows <- function(
    gr,
    zn,
    st,
    blockId,
    toRun,
    ldSketch,
    addSusieInf,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    methodArgs,
    verbose,
    fullFit,
    fullFitAlphaOnly,
    includeAllCs,
    serFallback,
    rFiniteResolved,
    rMismatch,
    rssControl,
    keepFullFit
) {
    z <- set_names(zn$z, zn$variantIds)
    ldMat <- .ldFromSketch(
        ldSketch,
        zn$variantIds,
        label = "fineMappingPipeline"
    )
    ents <- .fmFitRssBlock(
        z,
        ldMat,
        zn$n,
        toRun,
        addSusieInf,
        coverage,
        secondaryCoverage,
        signalCutoff,
        minAbsCorr,
        methodArgs,
        verbose,
        label = glue("GWAS (study='{st}', region='{blockId}')"),
        af = .fmAfByVar(gr, zn$variantIds),
        nVar = zn$nVar,
        fullFit = fullFit,
        fullFitAlphaOnly = fullFitAlphaOnly,
        includeAllCs = includeAllCs,
        serFallback = serFallback,
        rFinite = rFiniteResolved,
        rMismatch = rMismatch,
        rssControl = rssControl,
        keepFullFit = keepFullFit
    )
    map(names(ents), .fmGwasRowFor, st = st, blockId = blockId, ents = ents)
}

# All result rows for one GwasSumStats entry: cache hits + freshly-fitted
# tokens, or an empty set when the region was screened out. Returns
# list(rows, skipped) so the caller sums the skip flags functionally.
# @noRd
.fmGwasEntryRows <- function(
    i,
    data,
    studyCol,
    tokens,
    ldSketch,
    verbose,
    fineMappingResult,
    mafCutoff,
    macCutoff,
    imissCutoff,
    addSusieInf,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    methodArgs,
    fullFit,
    fullFitAlphaOnly,
    includeAllCs,
    serFallback,
    rFiniteResolved,
    rMismatch,
    rssControl,
    keepFullFit
) {
    st <- studyCol[[i]]
    gr <- .collectionEntry(data, i)
    skip <- .fmEntrySkipInfo(data, i)
    if (isTRUE(skip$skipped)) {
        if (verbose >= 1) {
            reason <- skip$reason
            msg <- glue(
                "fineMappingPipeline(GwasSumStats): study='{st}' region ",
                "skipped: {reason}"
            )
            inform(msg)
        }
        return(list(rows = list(), skipped = TRUE))
    }
    zn <- .fmExtractZn(
        gr,
        glue("fineMappingPipeline(GwasSumStats): study='{st}'"),
        ldSketch = ldSketch,
        cutoffs = .panelCutoffs(
            mafCutoff = mafCutoff,
            macCutoff = macCutoff,
            imissCutoff = imissCutoff
        )
    )
    blockId <- .fmGwasBlockId(gr)
    lookups <- map(
        tokens,
        .fmGwasLookup,
        fineMappingResult = fineMappingResult,
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
    computed <- .fmGwasFitRows(
        gr,
        zn,
        st,
        blockId,
        toRun,
        ldSketch = ldSketch,
        addSusieInf = addSusieInf,
        coverage = coverage,
        secondaryCoverage = secondaryCoverage,
        signalCutoff = signalCutoff,
        minAbsCorr = minAbsCorr,
        methodArgs = methodArgs,
        verbose = verbose,
        fullFit = fullFit,
        fullFitAlphaOnly = fullFitAlphaOnly,
        includeAllCs = includeAllCs,
        serFallback = serFallback,
        rFiniteResolved = rFiniteResolved,
        rMismatch = rMismatch,
        rssControl = rssControl,
        keepFullFit = keepFullFit
    )
    list(rows = c(cachedRows, computed), skipped = FALSE)
}

# Fit the still-to-run RSS tokens for one QtlSumStats entry and return one
# row-record per fitted token.
# @noRd
.fmRssFitRows <- function(
    i,
    st,
    ctx,
    tr,
    toRun,
    data,
    ldSketch,
    addSusieInf,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    methodArgs,
    verbose,
    fullFit,
    fullFitAlphaOnly,
    includeAllCs,
    serFallback,
    rFiniteResolved,
    rMismatch,
    rssControl,
    keepFullFit,
    mafCutoff,
    macCutoff,
    imissCutoff
) {
    entry <- .collectionEntry(data, i)
    zn <- .fmExtractZn(
        entry,
        glue(
            "fineMappingPipeline(QtlSumStats): entry {i} (study='{st}', ",
            "context='{ctx}', trait='{tr}')"
        ),
        ldSketch = ldSketch,
        cutoffs = .panelCutoffs(
            mafCutoff = mafCutoff,
            macCutoff = macCutoff,
            imissCutoff = imissCutoff
        )
    )
    z <- set_names(zn$z, zn$variantIds)
    ldMat <- .ldFromSketch(
        ldSketch,
        zn$variantIds,
        label = "fineMappingPipeline"
    )
    ents <- .fmFitRssBlock(
        z,
        ldMat,
        zn$n,
        toRun,
        addSusieInf,
        coverage,
        secondaryCoverage,
        signalCutoff,
        minAbsCorr,
        methodArgs,
        verbose,
        label = glue("(study='{st}', context='{ctx}', trait='{tr}')"),
        af = .fmAfByVar(entry, zn$variantIds),
        nVar = zn$nVar,
        fullFit = fullFit,
        fullFitAlphaOnly = fullFitAlphaOnly,
        includeAllCs = includeAllCs,
        serFallback = serFallback,
        rFinite = rFiniteResolved,
        rMismatch = rMismatch,
        rssControl = rssControl,
        keepFullFit = keepFullFit
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
    fineMappingResult,
    data,
    ldSketch,
    addSusieInf,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    methodArgs,
    verbose,
    fullFit,
    fullFitAlphaOnly,
    includeAllCs,
    serFallback,
    rFiniteResolved,
    rMismatch,
    rssControl,
    keepFullFit,
    mafCutoff,
    macCutoff,
    imissCutoff
) {
    st <- studyCol[i]
    ctx <- contextCol[i]
    tr <- traitCol[i]
    lookups <- map(
        univTokens,
        .fmQtlLookup,
        fineMappingResult = fineMappingResult,
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
    skip <- .fmEntrySkipInfo(data, i)
    if (isTRUE(skip$skipped)) {
        if (verbose >= 1) {
            reason <- skip$reason
            msg <- glue(
                "fineMappingPipeline(QtlSumStats): entry {i} ",
                "(study='{st}', context='{ctx}', trait='{tr}') ",
                "skipped: {reason}"
            )
            inform(msg)
        }
        return(list(rows = cachedRows, skipped = TRUE))
    }
    computed <- .fmRssFitRows(
        i,
        st,
        ctx,
        tr,
        toRun,
        data = data,
        ldSketch = ldSketch,
        addSusieInf = addSusieInf,
        coverage = coverage,
        secondaryCoverage = secondaryCoverage,
        signalCutoff = signalCutoff,
        minAbsCorr = minAbsCorr,
        methodArgs = methodArgs,
        verbose = verbose,
        fullFit = fullFit,
        fullFitAlphaOnly = fullFitAlphaOnly,
        includeAllCs = includeAllCs,
        serFallback = serFallback,
        rFiniteResolved = rFiniteResolved,
        rMismatch = rMismatch,
        rssControl = rssControl,
        keepFullFit = keepFullFit,
        mafCutoff = mafCutoff,
        macCutoff = macCutoff,
        imissCutoff = imissCutoff
    )
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


# Wrap one finemapping fit into a FineMappingRow via the surviving
# post-processing helpers (postprocessFinemappingFits +
# formatFinemappingOutput). Returns a bare FineMappingRow payload, ready
# to be inserted into a FineMappingResult.
# @noRd
# Look up residualization flags from the enclosing setMethod frame
# and call `getResidualized{Phenotypes,Genotypes}` with them. Each
# fineMappingPipeline / twasWeightsPipeline method exposes the four
# convenience flags listed in `.resFlagNames`; the wrapper threads them
# through to the accessor so per-call-site changes aren't needed.
.resFlagNames <- c(
    "phenotypeCovariatesToResidualize",
    "genotypeCovariatesToResidualize",
    "residualizePhenotypeCovariates",
    "residualizeGenotypeCovariates"
)

.resPickFlags <- function() {
    # Walk up from the immediate caller; the public setMethod frame is where
    # the user-facing args live. Each flag takes its value from the first
    # frame that supplies a non-NULL one -- assigning NULL into a list is a
    # no-op, so the original loop kept searching past those too.
    frames <- sys.frames()
    compact(set_names(
        map(.resFlagNames, .resFlagFromFrames, frames = frames),
        .resFlagNames
    ))
}

# @noRd
.resFrameValue <- function(fr, nm) {
    if (!exists(nm, envir = fr, inherits = FALSE)) {
        return(NULL)
    }
    get(nm, envir = fr, inherits = FALSE)
}

# @noRd
.resFlagIsSupplied <- function(value) {
    !is.null(value)
}

# @noRd
.resFlagFromFrames <- function(nm, frames) {
    detect(map(frames, .resFrameValue, nm = nm), .resFlagIsSupplied)
}

.fmResidPheno <- function(x, ...) {
    resArgs <- c(list(x = x, ...), .resPickFlags())
    exec(getResidualizedPhenotypes, !!!resArgs)
}

.fmResidGeno <- function(x, ...) {
    resArgs <- c(list(x = x, ...), .resPickFlags())
    exec(getResidualizedGenotypes, !!!resArgs)
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
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    csInput = NULL,
    af = NULL,
    n = NULL,
    region = NULL,
    trim = NULL,
    medianAbsCorr = NULL,
    conditionIdx = NULL,
    fullFit = NULL,
    fullFitAlphaOnly = NULL,
    includeAllCs = NULL
) {
    d <- .fmPostprocessDefaults(
        trim,
        medianAbsCorr,
        fullFit,
        fullFitAlphaOnly,
        includeAllCs
    )
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
        coverage,
        secondaryCoverage,
        signalCutoff,
        minAbsCorr,
        d
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
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    d
) {
    post <- postprocessFinemappingFits(
        fits = set_names(list(fit), method),
        dataX = dataX,
        dataY = dataY,
        af = af,
        n = n,
        coverage = coverage,
        secondaryCoverage = secondaryCoverage,
        signalCutoff = signalCutoff,
        minAbsCorr = minAbsCorr,
        medianAbsCorr = d$medianAbsCorr,
        region = region,
        csInput = csInput,
        conditionIdx = conditionIdx,
        trim = d$trim,
        fullFit = d$fullFit,
        fullFitAlphaOnly = d$fullFitAlphaOnly,
        includeAllCs = d$includeAllCs
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

# isTRUE-normalize the postprocess knobs, applying each knob's default when a
# caller passes NULL. Callers that expose these knobs (e.g. the RSS path) thread
# them explicitly via `p$trim` / `p$medianAbsCorr`; the rest pass NULL and take
# the default (trim = TRUE, medianAbsCorr = NULL).
# @noRd
.fmPostprocessDefaults <- function(
    trim,
    medianAbsCorr,
    fullFit,
    fullFitAlphaOnly,
    includeAllCs
) {
    list(
        trim = isTRUE(trim %||% TRUE),
        medianAbsCorr = medianAbsCorr,
        fullFit = isTRUE(fullFit %||% FALSE),
        fullFitAlphaOnly = isTRUE(fullFitAlphaOnly %||% TRUE),
        includeAllCs = isTRUE(includeAllCs %||% FALSE)
    )
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

# The canonical (non-reweighted) mvSuSiE mixture prior for residual variance
# `V`:
# create_mixture_prior(R) restricted to the group's conditions.
# @noRd
.fmCanonicalPrior <- function(V, conditionNames, R) {
    list(
        priorVariance = mvsusieR::create_mixture_prior(
            R = R,
            include_indices = conditionNames
        ),
        residualVariance = V
    )
}

# Rebuild the mvSuSiE data-driven *reweighted* mixture prior + residual variance
# from a stored mr.mash fit -- the lean payload
# (list(dataDrivenPriorMatrices, w0, V)) that mrmashWeights(retainFit = TRUE)
# attaches and twasWeightsPipeline keeps on the mrmash TwasWeightsRow. Shared
# by the fine-mapping mvsusie consumer and the twas mvsusie_weights consumer.
#
# Reproduces the deleted multivariate_pipeline.R reweighting bit-identically:
# rescaleCovW0(w0) collapses the expanded mr.mash weights onto the original
# data-driven covariance matrices ($U), filters to surviving components, and
# create_mixture_prior() wraps them, restricted to the fit's conditions
# (`conditionNames` = colnames(Y)). `V` becomes mvsusie's residual_variance.
# A NULL fit, NULL matrices, or no surviving component falls back to the
# canonical create_mixture_prior(R), matching the legacy `else` branch.
# Returns list(priorVariance, residualVariance) (residualVariance NULL only
# when no fit was supplied at all).
# @noRd
.buildMvsusieReweightedPrior <- function(
    fitParts,
    conditionNames,
    weightsTol = 1e-10,
    overrideU = NULL
) {
    R <- length(conditionNames)
    if (is.null(fitParts)) {
        return(.fmCanonicalPrior(NULL, conditionNames, R))
    }
    # `overrideU` (mode C / hybrid): reuse this fit's reweighted mixture weights
    # (w0) and residual variance (V) but swap in a different set of data-driven
    # covariance matrices -- the per-fold mash prior U. Components are matched
    # to w0 by name, so the override U must share component names with the fit.
    ddpm <- if (!is.null(overrideU)) {
        overrideU
    } else {
        fitParts$dataDrivenPriorMatrices
    }
    if (is.null(ddpm) || is.null(ddpm$U)) {
        return(.fmCanonicalPrior(fitParts$V, conditionNames, R))
    }
    rescaled <- rescaleCovW0(fitParts$w0)
    w0Updated <- rescaled[is_in(names(rescaled), names(ddpm$U))]
    if (length(w0Updated) == 0L) {
        return(.fmCanonicalPrior(fitParts$V, conditionNames, R))
    }
    mixture <- list(matrices = ddpm$U[names(w0Updated)], weights = w0Updated)
    list(
        priorVariance = mvsusieR::create_mixture_prior(
            mixture_prior = mixture,
            weights_tol = weightsTol,
            include_indices = conditionNames
        ),
        residualVariance = fitParts$V
    )
}

# Rows of a TwasWeights whose method is mr.mash and whose identity matches
# every axis the caller FIXED. A NULL axis means "do not filter on it".
# @noRd
.fmMrmashSelector <- function(twasWeights, study, trait, context) {
    axes <- compact(list(study = study, trait = trait, context = context))
    reduce(
        names(axes),
        .fmMrmashNarrow,
        twasWeights = twasWeights,
        axes = axes,
        .init = as.character(twasWeights$method) == "mrmash"
    )
}

# `.tupleColumn()` not `[[`: on a RangedTupleList `[[` extracts an ELEMENT
# (the variant set), while the identity axes live in mcols.
# @noRd
.fmMrmashNarrow <- function(sel, axis, twasWeights, axes) {
    sel & as.character(.tupleColumn(twasWeights, axis)) == axes[[axis]]
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
    twasWeights,
    study = NULL,
    trait = NULL,
    context = NULL
) {
    if (is.null(twasWeights)) {
        return(NULL)
    }
    # Each per-context mr.mash row of a joint group carries the SHARED joint
    # fit, so the consumer matches the FIXED axes and leaves the jointed axis
    # NULL (match-any). study/trait/context = NULL means "skip that axis".
    sel <- .fmMrmashSelector(twasWeights, study, trait, context)
    for (i in which(sel)) {
        f <- getFits(.twrRowParts(twasWeights, i))
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
    twasWeights,
    study = NULL,
    trait = NULL,
    context = NULL
) {
    if (is.null(twasWeights)) {
        return(NULL)
    }
    sel <- .fmMrmashSelector(twasWeights, study, trait, context)
    for (i in which(sel)) {
        cv <- getCvResult(.twrRowParts(twasWeights, i))
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

# Per-column marginal-association z-scores of y on each column of X (univariate
# regression z = betahat / sebetahat), used by the individual-level absZ screen.
# @noRd
.marginalZ <- function(X, y) {
    ur <- susieR::univariate_regression(X, y)
    ur$betahat / ur$sebetahat
}

# Single-effect (SER) pre-screen, individual-level. Reports whether a
# residualized (X, y) block shows a strong enough signal (by the chosen metric)
# to be worth a full fit. `screen` is a screen spec (see .asScreen): a legacy
# PIP cutoff (numeric scalar, 0 = off) OR a resolved list(metric, cutoff) for
# one of pip / absZ / bf / logBf. Ports the deleted multivariate_pipeline.R
# `skipConditions` / susie_twas `pip_cutoff_to_skip` logic (the individual-level
# analog of the sumstat-path `.applyEntryScreen`):
#   * no screen (NULL / 0 / non-scalar numeric) -> always keep.
#   * pip cutoff < 0 uses the adaptive 3 / nVariants threshold.
#   * absZ needs no susie fit; pip/bf/logBf fit susie L = 1 once (its
#     $lbf_variable gives the per-variant logBF for bf/logBf, $pip for pip).
#   * NA entries of `y` are dropped before fitting.
# The screen is advisory: too few samples/variants or a fit failure returns
# `fallback` -- TRUE (default) keeps the block rather than discard a potentially
# real signal (fineMapping / joint paths); colocboost passes FALSE to drop an
# outcome it cannot screen. This is the single L = 1 SuSiE pre-screen shared by
# .fmSerScreenColumns (joint) and .cbPipSkipOutcomes (colocboost).
# @noRd
.fmSerScreen <- function(X, y, screen, fallback = TRUE) {
    scr <- .asScreen(screen)
    if (is.null(scr)) {
        return(TRUE)
    }
    ok <- !is.na(y)
    if (sum(ok) < 2L || ncol(X) < 1L) {
        return(fallback)
    }
    raw <- X[ok, , drop = FALSE]
    # susieR needs a double X.
    Xs <- if (is.double(raw)) raw else `storage.mode<-`(raw, "double")
    ys <- y[ok]
    metric <- scr$metric
    cutoff <- scr$cutoff
    if (metric == "absZ") {
        z <- try_fetch(.marginalZ(Xs, ys), error = function(cnd) NULL)
        if (is.null(z)) {
            return(fallback)
        }
        return(any(abs(z) > cutoff, na.rm = TRUE))
    }
    fit <- try_fetch(
        suppressMessages(susieR::susie(Xs, ys, L = 1L)),
        error = function(cnd) NULL
    )
    if (is.null(fit)) {
        return(fallback)
    }
    if (metric == "pip") {
        thr <- if (cutoff < 0) 3 / ncol(Xs) else cutoff
        return(any(fit$pip > thr, na.rm = TRUE))
    }
    maxLbf <- suppressWarnings(max(as.numeric(fit$lbf_variable), na.rm = TRUE))
    if (!is.finite(maxLbf)) {
        return(fallback)
    }
    # bf: cutoff on the raw BF scale -> compare in log space; logBf: log scale.
    maxLbf > (if (metric == "bf") log(cutoff) else cutoff)
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
    fm <- names(formals(match.fun(fn)))
    if (!is_in("methodArgs", fm)) {
        return(args)
    }
    isFormal <- is_in(names(args), setdiff(fm, "methodArgs"))
    c(args[isFormal], list(methodArgs = args[!isFormal]))
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

# Per-fold mvsusie weights. Reuses the data-driven reweighted prior + residual
# covariance from the full-data mr.mash fit on every fold -- the prior is over
# conditions, identical across folds (only samples are held out). NULL mvPrior
# -> canonical prior (unchanged behavior).
# @noRd
.fmFoldWeightsMv <- function(Xtr, Ytr, coverage, userArgs, mvPrior) {
    baseArgs <- list(
        X = Xtr,
        Y = Ytr,
        coverage = coverage,
        prior_variance = if (is.null(mvPrior)) {
            mvsusieR::create_mixture_prior(R = ncol(Ytr))
        } else {
            mvPrior$priorVariance
        }
    )
    withPrior <- list_assign(
        baseArgs,
        !!!compact(list(residual_variance = mvPrior$residualVariance))
    )
    mvArgs <- .fmMergeUserArgs(withPrior, "mvsusie", userArgs)
    fit <- exec(fitMvsusie, !!!.splitMethodArgs(fitMvsusie, mvArgs))
    raw <- as.matrix(mvsusieWeights(mvsusieFit = fit))
    W <- `rownames<-`(raw, rownames(raw) %||% colnames(Xtr))
    `attr<-`(W, "fit", .fmLeanFoldFit(fit, "mvsusie"))
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
        retainFits = TRUE
    )
    # Recorded at the producer so a consumer can refuse fits trained on a
    # different split; without it a mismatched partition would leak held-out
    # samples into training and inflate the CV metrics.
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
    cisWindow
) {
    if (is.null(rg)) {
        .fmResidGeno(
            data,
            contexts = ctx,
            traitId = traitId,
            cisWindow = cisWindow,
            samples = samples
        )
    } else {
        .fmResidGeno(data, contexts = ctx, region = rg, samples = samples)
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
.fmFitXBlockP <- function(
    X,
    y,
    tokens,
    addSusieInf,
    ctx,
    label,
    afVec,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    methodArgs,
    verbose,
    cvFolds,
    cvThreads,
    samplePartition,
    seed,
    fullFit,
    fullFitAlphaOnly,
    includeAllCs
) {
    .fmFitXBlock(
        X,
        y,
        tokens,
        addSusieInf,
        coverage,
        secondaryCoverage,
        signalCutoff,
        minAbsCorr,
        methodArgs,
        verbose,
        ctx,
        label,
        cvFolds = cvFolds,
        cvThreads = cvThreads,
        samplePartition = samplePartition,
        af = afVec,
        fullFit = fullFit,
        fullFitAlphaOnly = fullFitAlphaOnly,
        includeAllCs = includeAllCs,
        seed = seed
    )
}

# Fit one univariate X window for (ctx, trait): residualize genotype, align to
# Y, SER pre-screen, then .fmFitXBlock. Errors when too few shared samples
# (a hard data problem), returns list() when the window screens out.
# @noRd
.fmUnivBlockFit <- function(
    rg,
    ctx,
    tid,
    Y,
    toRun,
    data,
    cisWindow,
    screen,
    addSusieInf,
    verbose,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    methodArgs,
    cvFolds,
    cvThreads,
    samplePartition,
    seed,
    fullFit,
    fullFitAlphaOnly,
    includeAllCs
) {
    allX <- .fmResidGenoBlock(
        ctx,
        tid,
        rg,
        rownames(Y),
        data = data,
        cisWindow = cisWindow
    )
    common <- intersect(rownames(allX), rownames(Y))
    if (length(common) < 2L) {
        msg <- glue(
            "fineMappingPipeline: too few shared samples between ",
            "residualized X and Y for (context='{ctx}', trait='{tid}')."
        )
        abort(msg)
    }
    X <- allX[common, , drop = FALSE]
    yBlock <- Y[common, , drop = FALSE]
    y <- if (ncol(yBlock) > 1L) {
        yBlock[, 1L, drop = TRUE]
    } else {
        drop(yBlock)
    }
    if (!.fmSerScreen(X, y, screen)) {
        if (verbose >= 1) {
            msg <- glue(
                "Skipping (context='{ctx}', trait='{tid}'): SER ",
                "pre-screen found no signal above the cutoff."
            )
            inform(msg)
        }
        return(list())
    }
    afVec <- .fmAfForX(
        data,
        X,
        traitId = tid,
        region = rg,
        cisWindow = cisWindow
    )
    .fmFitXBlockP(
        X,
        y,
        toRun,
        addSusieInf,
        ctx,
        tid,
        afVec,
        coverage = coverage,
        secondaryCoverage = secondaryCoverage,
        signalCutoff = signalCutoff,
        minAbsCorr = minAbsCorr,
        methodArgs = methodArgs,
        verbose = verbose,
        cvFolds = cvFolds,
        cvThreads = cvThreads,
        samplePartition = samplePartition,
        seed = seed,
        fullFit = fullFit,
        fullFitAlphaOnly = fullFitAlphaOnly,
        includeAllCs = includeAllCs
    )
}

# All univariate row-records for one (context, trait): cache hits, then (when
# tokens remain) per-window fits merged per token.
# @noRd
.fmUnivTraitRows <- function(
    tid,
    ctx,
    data,
    study,
    univTokens,
    xRegions,
    naAction,
    fineMappingResult,
    cisWindow,
    screen,
    addSusieInf,
    verbose,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    methodArgs,
    cvFolds,
    cvThreads,
    samplePartition,
    seed,
    fullFit,
    fullFitAlphaOnly,
    includeAllCs
) {
    lookups <- map(
        univTokens,
        .fmUnivLookup,
        fineMappingResult = fineMappingResult,
        study = study,
        ctx = ctx,
        tid = tid
    )
    cachedRows <- map(
        keep(lookups, .fmHasCached),
        .fmUnivCachedRow,
        study = study,
        ctx = ctx,
        tid = tid
    )
    toRun <- map_chr(keep(lookups, .fmNotCached), "tk")
    if (length(toRun) == 0L) {
        return(cachedRows)
    }
    Y <- .fmResidPheno(
        data,
        contexts = ctx,
        traitId = tid,
        naAction = naAction
    )
    blockEntries <- map(
        xRegions,
        .fmUnivBlockFit,
        ctx = ctx,
        tid = tid,
        data = data,
        cisWindow = cisWindow,
        screen = screen,
        addSusieInf = addSusieInf,
        verbose = verbose,
        coverage = coverage,
        secondaryCoverage = secondaryCoverage,
        signalCutoff = signalCutoff,
        minAbsCorr = minAbsCorr,
        methodArgs = methodArgs,
        cvFolds = cvFolds,
        cvThreads = cvThreads,
        samplePartition = samplePartition,
        seed = seed,
        fullFit = fullFit,
        fullFitAlphaOnly = fullFitAlphaOnly,
        includeAllCs = includeAllCs,
        Y = Y,
        toRun = toRun
    )
    computed <- .fmMergeTokenRows(study, ctx, tid, toRun, blockEntries)
    c(cachedRows, computed)
}

# Fit one PCA X window for (ctx, pcName): residualize genotype, align to the
# PC scores, SER pre-screen, then univariate-susie .fmFitXBlock. Returns
# list() when too few shared samples or the window screens out (soft skip).
# @noRd
.fmPcaBlockFit <- function(
    rg,
    ctx,
    traits,
    pcName,
    pcY,
    samples,
    data,
    cisWindow,
    screen,
    addSusieInf,
    verbose,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    methodArgs,
    cvFolds,
    cvThreads,
    samplePartition,
    seed,
    fullFit,
    fullFitAlphaOnly,
    includeAllCs
) {
    X <- .fmResidGenoBlock(
        ctx,
        traits,
        rg,
        samples,
        data = data,
        cisWindow = cisWindow
    )
    common <- intersect(rownames(X), names(pcY))
    if (length(common) < 2L) {
        return(list())
    }
    Xb <- X[common, , drop = FALSE]
    if (!.fmSerScreen(Xb, pcY[common], screen)) {
        return(list())
    }
    afVec <- .fmAfForX(
        data,
        Xb,
        traitId = traits,
        region = rg,
        cisWindow = cisWindow
    )
    .fmFitXBlockP(
        Xb,
        pcY[common],
        "susie",
        FALSE,
        ctx,
        pcName,
        afVec,
        coverage = coverage,
        secondaryCoverage = secondaryCoverage,
        signalCutoff = signalCutoff,
        minAbsCorr = minAbsCorr,
        methodArgs = methodArgs,
        verbose = verbose,
        cvFolds = cvFolds,
        cvThreads = cvThreads,
        samplePartition = samplePartition,
        seed = seed,
        fullFit = fullFit,
        fullFitAlphaOnly = fullFitAlphaOnly,
        includeAllCs = includeAllCs
    )
}

# Row-records for one PC pseudo-trait: a cache hit, else per-window fits merged
# into the single "susie" token (trait = the PC name).
# @noRd
.fmPcaScoreRows <- function(
    pcName,
    ctx,
    traits,
    scores,
    fineMappingResult,
    study,
    xRegions,
    data,
    cisWindow,
    screen,
    addSusieInf,
    verbose,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    methodArgs,
    cvFolds,
    cvThreads,
    samplePartition,
    seed,
    fullFit,
    fullFitAlphaOnly,
    includeAllCs
) {
    cached <- .fmCacheLookup(fineMappingResult, study, ctx, pcName, "susie")
    if (!is.null(cached)) {
        return(list(.fmQtlRow(study, ctx, pcName, "susie", cached)))
    }
    pcY <- scores[, pcName]
    samples <- rownames(scores)
    blockEntries <- map(
        xRegions,
        .fmPcaBlockFit,
        ctx = ctx,
        traits = traits,
        pcName = pcName,
        pcY = pcY,
        samples = samples,
        data = data,
        cisWindow = cisWindow,
        screen = screen,
        addSusieInf = addSusieInf,
        verbose = verbose,
        coverage = coverage,
        secondaryCoverage = secondaryCoverage,
        signalCutoff = signalCutoff,
        minAbsCorr = minAbsCorr,
        methodArgs = methodArgs,
        cvFolds = cvFolds,
        cvThreads = cvThreads,
        samplePartition = samplePartition,
        seed = seed,
        fullFit = fullFit,
        fullFitAlphaOnly = fullFitAlphaOnly,
        includeAllCs = includeAllCs
    )
    .fmMergeTokenRows(study, ctx, pcName, "susie", blockEntries)
}

# All usePCA row-records for one context: PCA-reduce the multi-trait phenotype
# and fine-map each top PC. A single-trait context (or one with no usable PC
# scores) contributes nothing.
# @noRd
.fmPcaContextRows <- function(
    ctx,
    data,
    study,
    perCtxTraits,
    nPCs,
    naAction,
    fineMappingResult,
    xRegions,
    cisWindow,
    screen,
    addSusieInf,
    verbose,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    methodArgs,
    cvFolds,
    cvThreads,
    samplePartition,
    seed,
    fullFit,
    fullFitAlphaOnly,
    includeAllCs
) {
    traits <- perCtxTraits[[ctx]]
    if (length(traits) < 2L) {
        return(list())
    }
    Yctx <- .fmResidPheno(
        data,
        contexts = ctx,
        traitId = traits,
        naAction = naAction
    )
    scores <- .fmTopPcScores(Yctx, nPCs)
    if (is.null(scores)) {
        return(list())
    }
    if (verbose >= 1) {
        nPc <- ncol(scores)
        nTr <- length(traits)
        msg <- glue(
            "usePCA: fine-mapping {nPc} top PC(s) of context='{ctx}' ",
            "({nTr} traits) ..."
        )
        inform(msg)
    }
    list_flatten(map(
        colnames(scores),
        .fmPcaScoreRows,
        ctx = ctx,
        traits = traits,
        scores = scores,
        fineMappingResult = fineMappingResult,
        study = study,
        xRegions = xRegions,
        data = data,
        cisWindow = cisWindow,
        screen = screen,
        addSusieInf = addSusieInf,
        verbose = verbose,
        coverage = coverage,
        secondaryCoverage = secondaryCoverage,
        signalCutoff = signalCutoff,
        minAbsCorr = minAbsCorr,
        methodArgs = methodArgs,
        cvFolds = cvFolds,
        cvThreads = cvThreads,
        samplePartition = samplePartition,
        seed = seed,
        fullFit = fullFit,
        fullFitAlphaOnly = fullFitAlphaOnly,
        includeAllCs = includeAllCs
    ))
}

# QtlDataset method
# =============================================================================

# Run the joint engine for a QtlDataset with the given spec + token set. Shared
# by the explicit-jointSpecification path and the auto-detected multivariate
# path; only the spec and token arguments differ between them.
# @noRd
.fmQdsJointDispatch <- function(
    jointSpec,
    tokens,
    methodArgs,
    data,
    contexts,
    traitId,
    cisWindow,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    verbose,
    xRegions,
    twasWeights,
    dataDrivenPriorWeightsCutoff,
    cvFolds,
    cvThreads,
    samplePartition,
    screen,
    fineMappingResult,
    fullFit,
    fullFitAlphaOnly,
    includeAllCs,
    seed
) {
    .fmDispatchJointSpecsQtlDataset(
        jointSpec,
        data,
        tokens,
        contexts,
        traitId,
        cisWindow,
        coverage,
        secondaryCoverage,
        signalCutoff,
        minAbsCorr,
        verbose,
        methodArgs = methodArgs,
        xRegions = xRegions,
        twasWeights = twasWeights,
        dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
        cvFolds = cvFolds,
        cvThreads = cvThreads,
        samplePartition = samplePartition,
        pipCutoffToSkip = screen,
        fineMappingResult = fineMappingResult,
        fullFit = fullFit,
        fullFitAlphaOnly = fullFitAlphaOnly,
        includeAllCs = includeAllCs,
        seed = seed
    )
}

# Resolve the signal screen, apply per-call filter overrides to a validated
# copy of the dataset, and derive the X windows. Rejects the region + cisWindow
# combination. Returns `p` extended with data / screen / xRegions.
# @noRd
.fmQdsResolveInputs <- function(
    data,
    pipCutoffToSkip,
    absZCutoffToSkip,
    bfCutoffToSkip,
    logBfCutoffToSkip,
    mafCutoff,
    macCutoff,
    xvarCutoff,
    imissCutoff,
    keepIndel,
    keepSamples,
    keepVariants,
    region,
    cisWindow,
    jointRegions
) {
    screen <- .resolveScreenMetric(
        pipCutoffToSkip,
        absZCutoffToSkip,
        bfCutoffToSkip,
        logBfCutoffToSkip
    )
    data <- .qtlApplyFilterOverrides(
        data,
        mafCutoff,
        macCutoff,
        xvarCutoff,
        imissCutoff,
        keepIndel,
        keepSamples,
        keepVariants
    )
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

# Normalize methods, capability-check them, and run any EXPLICIT
# jointSpecification up front (removing the joint methods from the per-tuple
# token set). `exhaustedByJoint` flags the case where the explicit spec
# consumed every method.
# @noRd
.fmQdsResolveTokens <- function(
    jointSpecification,
    methods,
    L,
    Lgreedy,
    data,
    contexts,
    traitId,
    cisWindow,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    verbose,
    xRegions,
    twasWeights,
    dataDrivenPriorWeightsCutoff,
    cvFolds,
    cvThreads,
    samplePartition,
    screen,
    fineMappingResult,
    fullFit,
    fullFitAlphaOnly,
    includeAllCs,
    seed
) {
    parsedJointSpec <- parseJointSpecification(jointSpecification, data)
    norm <- .fmNormalizeMethods(methods, L = L, Lgreedy = Lgreedy)
    allTokens <- norm$tokens
    allMethodArgs <- norm$methodArgs
    .fmCheckMethodCapabilities(allTokens, "QtlDataset")
    hadJointSpec <- length(parsedJointSpec) > 0L
    jointResult <- if (!hadJointSpec) {
        NULL
    } else {
        .fmQdsJointDispatch(
            parsedJointSpec,
            intersect(allTokens, c("mvsusie", "fsusie")),
            allMethodArgs,
            data = data,
            contexts = contexts,
            traitId = traitId,
            cisWindow = cisWindow,
            coverage = coverage,
            secondaryCoverage = secondaryCoverage,
            signalCutoff = signalCutoff,
            minAbsCorr = minAbsCorr,
            verbose = verbose,
            xRegions = xRegions,
            twasWeights = twasWeights,
            dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
            cvFolds = cvFolds,
            cvThreads = cvThreads,
            samplePartition = samplePartition,
            screen = screen,
            fineMappingResult = fineMappingResult,
            fullFit = fullFit,
            fullFitAlphaOnly = fullFitAlphaOnly,
            includeAllCs = includeAllCs,
            seed = seed
        )
    }
    # The joint dispatch consumed the multivariate tokens; the univariate
    # phase runs on whatever is left.
    narrowed <- .fmNarrowAfterJoint(hadJointSpec, allTokens, allMethodArgs)
    list(
        tokens = narrowed$tokens,
        methodArgs = narrowed$methodArgs,
        jointResult = jointResult,
        exhaustedByJoint = hadJointSpec && length(narrowed$tokens) == 0L
    )
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
.fmQdsDispatchRows <- function(
    addSusieInf,
    cisWindow,
    coverage,
    cvFolds,
    cvThreads,
    data,
    fineMappingResult,
    fullFit,
    fullFitAlphaOnly,
    includeAllCs,
    methodArgs,
    minAbsCorr,
    nPCs,
    naAction,
    perCtxTraits,
    samplePartition,
    screen,
    secondaryCoverage,
    seed,
    signalCutoff,
    study,
    univTokens,
    useCtx,
    usePCA,
    verbose,
    xRegions
) {
    univRows <- if (length(univTokens) > 0L) {
        list_flatten(map(
            useCtx,
            .fmUnivContextRows,
            data = data,
            study = study,
            univTokens = univTokens,
            xRegions = xRegions,
            naAction = naAction,
            fineMappingResult = fineMappingResult,
            cisWindow = cisWindow,
            screen = screen,
            addSusieInf = addSusieInf,
            verbose = verbose,
            coverage = coverage,
            secondaryCoverage = secondaryCoverage,
            signalCutoff = signalCutoff,
            minAbsCorr = minAbsCorr,
            methodArgs = methodArgs,
            cvFolds = cvFolds,
            cvThreads = cvThreads,
            samplePartition = samplePartition,
            seed = seed,
            fullFit = fullFit,
            fullFitAlphaOnly = fullFitAlphaOnly,
            includeAllCs = includeAllCs,
            perCtxTraits = perCtxTraits
        ))
    } else {
        list()
    }
    pcaRows <- if (isTRUE(usePCA)) {
        list_flatten(map(
            useCtx,
            .fmPcaContextRows,
            data = data,
            study = study,
            xRegions = xRegions,
            naAction = naAction,
            fineMappingResult = fineMappingResult,
            cisWindow = cisWindow,
            screen = screen,
            addSusieInf = addSusieInf,
            verbose = verbose,
            coverage = coverage,
            secondaryCoverage = secondaryCoverage,
            signalCutoff = signalCutoff,
            minAbsCorr = minAbsCorr,
            methodArgs = methodArgs,
            cvFolds = cvFolds,
            cvThreads = cvThreads,
            samplePartition = samplePartition,
            seed = seed,
            fullFit = fullFit,
            fullFitAlphaOnly = fullFitAlphaOnly,
            includeAllCs = includeAllCs,
            perCtxTraits = perCtxTraits,
            nPCs = nPCs
        ))
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
    data,
    contexts,
    traitId,
    cisWindow,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    verbose,
    xRegions,
    twasWeights,
    dataDrivenPriorWeightsCutoff,
    cvFolds,
    cvThreads,
    samplePartition,
    screen,
    fineMappingResult,
    fullFit,
    fullFitAlphaOnly,
    includeAllCs,
    seed,
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
        data = data,
        contexts = contexts,
        traitId = traitId,
        cisWindow = cisWindow,
        coverage = coverage,
        secondaryCoverage = secondaryCoverage,
        signalCutoff = signalCutoff,
        minAbsCorr = minAbsCorr,
        verbose = verbose,
        xRegions = xRegions,
        twasWeights = twasWeights,
        dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
        cvFolds = cvFolds,
        cvThreads = cvThreads,
        samplePartition = samplePartition,
        screen = screen,
        fineMappingResult = fineMappingResult,
        fullFit = fullFit,
        fullFitAlphaOnly = fullFitAlphaOnly,
        includeAllCs = includeAllCs,
        seed = seed
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

# QtlDataset fine-mapping worker. `p` is the setMethod's captured arguments;
# each phase extends the bundle via list_modify so the dispatch helpers read
# everything from `p`. The susieInf chaining is applied downstream inside
# .fmFitXBlock / .fmFitRssBlock (which recompute .fmResolveSusieChain), so no
# chain config is threaded here.
# @noRd
.fmPipelineQtlDataset <- function(
    data,
    methods,
    contexts,
    traitId,
    region,
    cisWindow,
    mafCutoff,
    macCutoff,
    xvarCutoff,
    imissCutoff,
    keepIndel,
    keepSamples,
    keepVariants,
    jointRegions,
    jointSpecification,
    addSusieInf,
    L,
    Lgreedy,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    fineMappingResult,
    twasWeights,
    dataDrivenPriorWeightsCutoff,
    verbose,
    cvFolds,
    cvThreads,
    samplePartition,
    seed,
    naAction,
    usePCA,
    nPCs,
    pipCutoffToSkip,
    absZCutoffToSkip,
    bfCutoffToSkip,
    logBfCutoffToSkip,
    fullFit,
    fullFitAlphaOnly,
    includeAllCs
) {
    naAction <- arg_match(naAction, c("drop", "impute"))
    if (!is.null(seed)) {
        withr::local_seed(as.integer(seed))
    }
    # Each stage returns only the values it derives; nothing is grafted onto a
    # captured environment.
    resolved <- .fmQdsResolveInputs(
        data = data,
        pipCutoffToSkip = pipCutoffToSkip,
        absZCutoffToSkip = absZCutoffToSkip,
        bfCutoffToSkip = bfCutoffToSkip,
        logBfCutoffToSkip = logBfCutoffToSkip,
        mafCutoff = mafCutoff,
        macCutoff = macCutoff,
        xvarCutoff = xvarCutoff,
        imissCutoff = imissCutoff,
        keepIndel = keepIndel,
        keepSamples = keepSamples,
        keepVariants = keepVariants,
        region = region,
        cisWindow = cisWindow,
        jointRegions = jointRegions
    )
    data <- resolved$data
    screen <- resolved$screen
    xRegions <- resolved$xRegions
    rt <- .fmQdsResolveTokens(
        jointSpecification = jointSpecification,
        methods = methods,
        L = L,
        Lgreedy = Lgreedy,
        data = data,
        contexts = contexts,
        traitId = traitId,
        cisWindow = cisWindow,
        coverage = coverage,
        secondaryCoverage = secondaryCoverage,
        signalCutoff = signalCutoff,
        minAbsCorr = minAbsCorr,
        verbose = verbose,
        xRegions = xRegions,
        twasWeights = twasWeights,
        dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
        cvFolds = cvFolds,
        cvThreads = cvThreads,
        samplePartition = samplePartition,
        screen = screen,
        fineMappingResult = fineMappingResult,
        fullFit = fullFit,
        fullFitAlphaOnly = fullFitAlphaOnly,
        includeAllCs = includeAllCs,
        seed = seed
    )
    if (rt$exhaustedByJoint) {
        if (is.null(rt$jointResult)) {
            msg <- glue(
                "fineMappingPipeline(QtlDataset): no joint fits produced. ",
                "Check that the jointSpecification scope intersects the ",
                "available studies / contexts / traits."
            )
            abort(msg)
        }
        return(rt$jointResult)
    }
    tokens <- rt$tokens
    methodArgs <- rt$methodArgs
    jointResult <- rt$jointResult
    ctxInfo <- .fmQdsResolveContexts(
        data = data,
        contexts = contexts,
        traitId = traitId,
        region = region
    )
    study <- ctxInfo$study
    useCtx <- ctxInfo$useCtx
    perCtxTraits <- ctxInfo$perCtxTraits
    nCtx <- ctxInfo$nCtx
    nTraits <- ctxInfo$nTraits
    split <- .fmQdsSplitTokens(tokens, nCtx, nTraits)
    univTokens <- split$univTokens
    mvTokens <- split$mvTokens
    fsTokens <- split$fsTokens
    rows <- .fmQdsDispatchRows(
        addSusieInf = addSusieInf,
        cisWindow = cisWindow,
        coverage = coverage,
        cvFolds = cvFolds,
        cvThreads = cvThreads,
        data = data,
        fineMappingResult = fineMappingResult,
        fullFit = fullFit,
        fullFitAlphaOnly = fullFitAlphaOnly,
        includeAllCs = includeAllCs,
        methodArgs = methodArgs,
        minAbsCorr = minAbsCorr,
        nPCs = nPCs,
        naAction = naAction,
        perCtxTraits = perCtxTraits,
        samplePartition = samplePartition,
        screen = screen,
        secondaryCoverage = secondaryCoverage,
        seed = seed,
        signalCutoff = signalCutoff,
        study = study,
        univTokens = univTokens,
        useCtx = useCtx,
        usePCA = usePCA,
        verbose = verbose,
        xRegions = xRegions
    )
    .fmQdsAssemble(
        data,
        rows,
        .fmQdsAutoJoint(
            data = data,
            contexts = contexts,
            traitId = traitId,
            cisWindow = cisWindow,
            coverage = coverage,
            secondaryCoverage = secondaryCoverage,
            signalCutoff = signalCutoff,
            minAbsCorr = minAbsCorr,
            verbose = verbose,
            xRegions = xRegions,
            twasWeights = twasWeights,
            dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
            cvFolds = cvFolds,
            cvThreads = cvThreads,
            samplePartition = samplePartition,
            screen = screen,
            fineMappingResult = fineMappingResult,
            fullFit = fullFit,
            fullFitAlphaOnly = fullFitAlphaOnly,
            includeAllCs = includeAllCs,
            seed = seed,
            mvTokens = mvTokens,
            fsTokens = fsTokens,
            jointResult = jointResult,
            methodArgs = methodArgs,
            nCtx = nCtx,
            nTraits = nTraits
        )
    )
}

#' @rdname fineMappingPipeline
#' @importFrom purrr list_c list_flatten list_modify list_rbind list_assign
#' @importFrom purrr discard walk2 detect zap
#' @export
setMethod(
    "fineMappingPipeline",
    "QtlDataset",
    function(
        data,
        methods,
        contexts = NULL,
        traitId = NULL,
        region = NULL,
        cisWindow = NULL,
        # Per-call genotype-filter overrides; NULL = use the QtlDataset's
        # construct-time slot value (applied lazily at extraction).
        mafCutoff = NULL,
        macCutoff = NULL,
        xvarCutoff = NULL,
        imissCutoff = NULL,
        keepIndel = NULL,
        keepSamples = NULL,
        keepVariants = NULL,
        jointRegions = FALSE,
        jointSpecification = NULL,
        addSusieInf = TRUE,
        L = 10L,
        Lgreedy = NULL,
        coverage = 0.95,
        secondaryCoverage = c(0.7, 0.5),
        signalCutoff = 0.025,
        minAbsCorr = 0.8,
        medianAbsCorr = NULL,
        fineMappingResult = NULL,
        cvFolds = 0,
        cvThreads = 1,
        samplePartition = NULL,
        pipCutoffToSkip = 0,
        absZCutoffToSkip = 0,
        bfCutoffToSkip = 0,
        logBfCutoffToSkip = 0,
        usePCA = FALSE,
        nPCs = 10L,
        seed = NULL,
        twasWeights = NULL,
        dataDrivenPriorWeightsCutoff = 1e-10,
        naAction = c("drop", "impute"),
        verbose = 1,
        trim = TRUE,
        fullFit = FALSE,
        fullFitAlphaOnly = TRUE,
        includeAllCs = FALSE,
        phenotypeCovariatesToResidualize = NULL,
        genotypeCovariatesToResidualize = NULL,
        residualizePhenotypeCovariates = TRUE,
        residualizeGenotypeCovariates = TRUE,
        ...
    ) {
        .fmPipelineQtlDataset(
            data = data,
            methods = methods,
            contexts = contexts,
            traitId = traitId,
            region = region,
            cisWindow = cisWindow,
            mafCutoff = mafCutoff,
            macCutoff = macCutoff,
            xvarCutoff = xvarCutoff,
            imissCutoff = imissCutoff,
            keepIndel = keepIndel,
            keepSamples = keepSamples,
            keepVariants = keepVariants,
            jointRegions = jointRegions,
            jointSpecification = jointSpecification,
            addSusieInf = addSusieInf,
            L = L,
            Lgreedy = Lgreedy,
            coverage = coverage,
            secondaryCoverage = secondaryCoverage,
            signalCutoff = signalCutoff,
            minAbsCorr = minAbsCorr,
            fineMappingResult = fineMappingResult,
            twasWeights = twasWeights,
            dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
            verbose = verbose,
            cvFolds = cvFolds,
            cvThreads = cvThreads,
            samplePartition = samplePartition,
            seed = seed,
            naAction = naAction,
            usePCA = usePCA,
            nPCs = nPCs,
            pipCutoffToSkip = pipCutoffToSkip,
            absZCutoffToSkip = absZCutoffToSkip,
            bfCutoffToSkip = bfCutoffToSkip,
            logBfCutoffToSkip = logBfCutoffToSkip,
            fullFit = fullFit,
            fullFitAlphaOnly = fullFitAlphaOnly,
            includeAllCs = includeAllCs
        )
    }
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
    fmArgs <- c(
        list(
            data = qd,
            methods = m,
            contexts = cfg$contexts,
            traitId = cfg$traitId,
            region = cfg$region,
            cisWindow = cfg$cisWindow,
            jointRegions = cfg$jointRegions,
            jointSpecification = NULL,
            addSusieInf = cfg$addSusieInf,
            coverage = cfg$coverage,
            secondaryCoverage = cfg$secondaryCoverage,
            signalCutoff = cfg$signalCutoff,
            minAbsCorr = cfg$minAbsCorr,
            fineMappingResult = cfg$fineMappingResult,
            cvFolds = cfg$cvFolds,
            cvThreads = cfg$cvThreads,
            samplePartition = cfg$samplePartition,
            pipCutoffToSkip = cfg$pipCutoffToSkip,
            absZCutoffToSkip = cfg$absZCutoffToSkip,
            bfCutoffToSkip = cfg$bfCutoffToSkip,
            logBfCutoffToSkip = cfg$logBfCutoffToSkip,
            seed = cfg$seed,
            naAction = cfg$naAction,
            verbose = cfg$verbose
        ),
        cfg$dotArgs
    )
    exec(fineMappingPipeline, !!!fmArgs)
}

# Embedded-sumstats fine-mapping worker for .multiStudyPipelineDriver: recurse
# fineMappingPipeline on the QtlSumStats with the sumstat-capable methods.
# @noRd
.fmSumStats <- function(ss, cfg) {
    m <- .fmFilterMethodsForKind(cfg$methods, "sumstatImpl")
    if (length(m) == 0L) {
        return(NULL)
    }
    fmArgs <- c(
        list(
            data = ss,
            methods = m,
            contexts = cfg$contexts,
            traitId = cfg$traitId,
            jointSpecification = NULL,
            addSusieInf = cfg$addSusieInf,
            coverage = cfg$coverage,
            secondaryCoverage = cfg$secondaryCoverage,
            signalCutoff = cfg$signalCutoff,
            minAbsCorr = cfg$minAbsCorr,
            fineMappingResult = cfg$fineMappingResult,
            verbose = cfg$verbose
        ),
        cfg$dotArgs
    )
    exec(fineMappingPipeline, !!!fmArgs)
}

# Resolve method tokens for a MultiStudyQtlDataset run and run any EXPLICIT
# jointSpecification (per-component axis dispatcher), removing the joint methods
# from the per-tuple recursion. Returns the still-pending tokens, the forwarded
# `methods` (kwargs-preserving), the joint result, and `exhaustedByJoint`.
# @noRd
.fmMsResolveTokens <- function(
    data,
    methods,
    contexts,
    traitId,
    cisWindow,
    xRegions,
    jointSpecification,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    twasWeights,
    dataDrivenPriorWeightsCutoff,
    verbose
) {
    parsedJointSpec <- parseJointSpecification(jointSpecification, data)
    norm <- .fmNormalizeMethods(methods)
    allTokens <- norm$tokens
    allMethodArgs <- norm$methodArgs
    .fmCheckMethodCapabilities(allTokens, "MultiStudyQtlDataset")
    hadJointSpec <- length(parsedJointSpec) > 0L
    jointResult <- if (!hadJointSpec) {
        NULL
    } else {
        .fmDispatchJointSpecsMultiStudy(
            parsedJointSpec,
            data,
            intersect(allTokens, c("mvsusie", "fsusie")),
            contexts,
            traitId,
            cisWindow,
            coverage,
            secondaryCoverage,
            signalCutoff,
            minAbsCorr,
            verbose,
            methodArgs = allMethodArgs,
            xRegions = xRegions,
            twasWeights = twasWeights,
            dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff
        )
    }
    # The joint dispatch consumed the multivariate tokens; the univariate
    # phase runs on whatever is left.
    narrowed <- .fmNarrowAfterJoint(hadJointSpec, allTokens, allMethodArgs)
    list(
        tokens = narrowed$tokens,
        methods = if (!hadJointSpec) {
            methods
        } else if (length(narrowed$methodArgs) > 0L) {
            narrowed$methodArgs
        } else {
            narrowed$tokens
        },
        jointResult = jointResult,
        exhaustedByJoint = hadJointSpec && length(narrowed$tokens) == 0L
    )
}

# Assemble the per-component driver config for a MultiStudyQtlDataset run:
# individual-capable methods route to the per-study QtlDatasets, sumstat-capable
# methods (incl. the sumstat-only `ser`) to the embedded QtlSumStats.
# @noRd
.fmMsConfig <- function(
    methods,
    contexts,
    traitId,
    region,
    cisWindow,
    jointRegions,
    addSusieInf,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    fineMappingResult,
    verbose,
    cvFolds,
    cvThreads,
    samplePartition,
    seed,
    naAction,
    pipCutoffToSkip,
    absZCutoffToSkip,
    bfCutoffToSkip,
    logBfCutoffToSkip,
    dotArgs
) {
    list(
        methods = methods,
        contexts = contexts,
        traitId = traitId,
        region = region,
        cisWindow = cisWindow,
        jointRegions = jointRegions,
        addSusieInf = addSusieInf,
        coverage = coverage,
        secondaryCoverage = secondaryCoverage,
        signalCutoff = signalCutoff,
        minAbsCorr = minAbsCorr,
        fineMappingResult = fineMappingResult,
        cvFolds = cvFolds,
        cvThreads = cvThreads,
        samplePartition = samplePartition,
        pipCutoffToSkip = pipCutoffToSkip,
        absZCutoffToSkip = absZCutoffToSkip,
        bfCutoffToSkip = bfCutoffToSkip,
        logBfCutoffToSkip = logBfCutoffToSkip,
        seed = seed,
        naAction = naAction,
        verbose = verbose,
        dotArgs = dotArgs
    )
}

# MultiStudyQtlDataset fine-mapping worker. `p` is the setMethod's captured
# arguments (plus `dotArgs`); after resolving the explicit joint spec it routes
# each remaining method to the components it supports via the shared
# multi-study driver.
# @noRd
.fmPipelineMultiStudy <- function(
    data,
    methods,
    contexts,
    traitId,
    region,
    cisWindow,
    jointRegions,
    jointSpecification,
    addSusieInf,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    fineMappingResult,
    twasWeights,
    dataDrivenPriorWeightsCutoff,
    verbose,
    cvFolds,
    cvThreads,
    samplePartition,
    seed,
    naAction,
    pipCutoffToSkip,
    absZCutoffToSkip,
    bfCutoffToSkip,
    logBfCutoffToSkip,
    dotArgs
) {
    naAction <- arg_match(naAction, c("drop", "impute"))
    if (!is.null(region) && !is.null(cisWindow)) {
        msg <- glue(
            "fineMappingPipeline(MultiStudyQtlDataset): specify either ",
            "`region` or `cisWindow`, not both."
        )
        abort(msg)
    }
    xRegions <- .makeXRegions(region, jointRegions)
    rt <- .fmMsResolveTokens(
        data = data,
        methods = methods,
        contexts = contexts,
        traitId = traitId,
        cisWindow = cisWindow,
        xRegions = xRegions,
        jointSpecification = jointSpecification,
        coverage = coverage,
        secondaryCoverage = secondaryCoverage,
        signalCutoff = signalCutoff,
        minAbsCorr = minAbsCorr,
        twasWeights = twasWeights,
        dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
        verbose = verbose
    )
    if (rt$exhaustedByJoint) {
        if (is.null(rt$jointResult)) {
            msg <- glue(
                "fineMappingPipeline(MultiStudyQtlDataset): no joint fits ",
                "produced. Check that the jointSpecification scope ",
                "intersects the available data."
            )
            abort(msg)
        }
        return(rt$jointResult)
    }
    methods <- rt$methods
    .multiStudyPipelineDriver(
        data,
        rt$jointResult,
        .fmPerStudy,
        .fmSumStats,
        .fmMsConfig(
            methods = methods,
            contexts = contexts,
            traitId = traitId,
            region = region,
            cisWindow = cisWindow,
            jointRegions = jointRegions,
            addSusieInf = addSusieInf,
            coverage = coverage,
            secondaryCoverage = secondaryCoverage,
            signalCutoff = signalCutoff,
            minAbsCorr = minAbsCorr,
            fineMappingResult = fineMappingResult,
            verbose = verbose,
            cvFolds = cvFolds,
            cvThreads = cvThreads,
            samplePartition = samplePartition,
            seed = seed,
            naAction = naAction,
            pipCutoffToSkip = pipCutoffToSkip,
            absZCutoffToSkip = absZCutoffToSkip,
            bfCutoffToSkip = bfCutoffToSkip,
            logBfCutoffToSkip = logBfCutoffToSkip,
            dotArgs = dotArgs
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
    function(
        data,
        methods,
        contexts = NULL,
        traitId = NULL,
        region = NULL,
        cisWindow = NULL,
        jointRegions = FALSE,
        jointSpecification = NULL,
        addSusieInf = TRUE,
        coverage = 0.95,
        secondaryCoverage = c(0.7, 0.5),
        signalCutoff = 0.025,
        minAbsCorr = 0.8,
        medianAbsCorr = NULL,
        fineMappingResult = NULL,
        twasWeights = NULL,
        dataDrivenPriorWeightsCutoff = 1e-10,
        cvFolds = 0,
        cvThreads = 1,
        samplePartition = NULL,
        pipCutoffToSkip = 0,
        absZCutoffToSkip = 0,
        bfCutoffToSkip = 0,
        logBfCutoffToSkip = 0,
        seed = NULL,
        naAction = c("drop", "impute"),
        verbose = 1,
        trim = TRUE,
        fullFit = FALSE,
        fullFitAlphaOnly = TRUE,
        includeAllCs = FALSE,
        phenotypeCovariatesToResidualize = NULL,
        genotypeCovariatesToResidualize = NULL,
        residualizePhenotypeCovariates = TRUE,
        residualizeGenotypeCovariates = TRUE,
        ...
    ) {
        .fmPipelineMultiStudy(
            data = data,
            methods = methods,
            contexts = contexts,
            traitId = traitId,
            region = region,
            cisWindow = cisWindow,
            jointRegions = jointRegions,
            jointSpecification = jointSpecification,
            addSusieInf = addSusieInf,
            coverage = coverage,
            secondaryCoverage = secondaryCoverage,
            signalCutoff = signalCutoff,
            minAbsCorr = minAbsCorr,
            fineMappingResult = fineMappingResult,
            twasWeights = twasWeights,
            dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
            verbose = verbose,
            cvFolds = cvFolds,
            cvThreads = cvThreads,
            samplePartition = samplePartition,
            seed = seed,
            naAction = naAction,
            pipCutoffToSkip = pipCutoffToSkip,
            absZCutoffToSkip = absZCutoffToSkip,
            bfCutoffToSkip = bfCutoffToSkip,
            logBfCutoffToSkip = logBfCutoffToSkip,
            dotArgs = list(...)
        )
    }
)


# =============================================================================
# QtlSumStats method
# =============================================================================

# Resolve the method tokens for a QtlSumStats run and run any EXPLICIT
# jointSpecification up front (mvsusie via the axis dispatcher), removing the
# joint methods from the per-tuple token set. `exhaustedByJoint` flags the case
# where an explicit spec consumed every method, so the caller returns the joint
# result directly (or errors when it produced nothing).
# @noRd
.fmQssResolveTokens <- function(
    data,
    methods,
    contexts,
    traitId,
    jointSpecification,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    fineMappingResult,
    twasWeights,
    dataDrivenPriorWeightsCutoff,
    verbose,
    fullFit,
    fullFitAlphaOnly,
    includeAllCs,
    mafCutoff,
    macCutoff,
    imissCutoff
) {
    parsedJointSpec <- parseJointSpecification(jointSpecification, data)
    norm <- .fmNormalizeMethods(methods)
    allTokens <- norm$tokens
    allMethodArgs <- norm$methodArgs
    .fmCheckMethodCapabilities(allTokens, "QtlSumStats")
    hadJointSpec <- length(parsedJointSpec) > 0L
    jointResult <- if (!hadJointSpec) {
        NULL
    } else {
        .fmDispatchJointSpecsQtlSumStats(
            parsedJointSpec,
            data,
            intersect(allTokens, "mvsusie"),
            contexts,
            traitId,
            coverage,
            secondaryCoverage,
            signalCutoff,
            minAbsCorr,
            verbose,
            methodArgs = allMethodArgs,
            twasWeights = twasWeights,
            dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
            fineMappingResult = fineMappingResult,
            fullFit = fullFit,
            fullFitAlphaOnly = fullFitAlphaOnly,
            includeAllCs = includeAllCs,
            mafCutoff = mafCutoff %||% 0,
            macCutoff = macCutoff %||% 0,
            imissCutoff = imissCutoff %||% 1
        )
    }
    # The joint dispatch consumed the multivariate tokens; the univariate
    # phase runs on whatever is left.
    narrowed <- .fmNarrowAfterJoint(hadJointSpec, allTokens, allMethodArgs)
    list(
        tokens = narrowed$tokens,
        methodArgs = narrowed$methodArgs,
        jointResult = jointResult,
        exhaustedByJoint = hadJointSpec && length(narrowed$tokens) == 0L
    )
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
# explicit-spec result already in `p$jointResult`.
# @noRd
.fmQssAutoJoint <- function(
    data,
    contexts,
    traitId,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    fineMappingResult,
    twasWeights,
    dataDrivenPriorWeightsCutoff,
    verbose,
    fullFit,
    fullFitAlphaOnly,
    includeAllCs,
    mafCutoff,
    macCutoff,
    imissCutoff,
    methodArgs,
    mvTokens,
    jointResult,
    ldSketch
) {
    if (length(mvTokens) == 0L) {
        return(jointResult)
    }
    autoJoint <- .fmDispatchJointSpecsQtlSumStats(
        list(list(axes = "context", scope = NULL)),
        data,
        mvTokens,
        contexts,
        traitId,
        coverage,
        secondaryCoverage,
        signalCutoff,
        minAbsCorr,
        verbose,
        methodArgs = methodArgs,
        twasWeights = twasWeights,
        dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
        fineMappingResult = fineMappingResult,
        fullFit = fullFit,
        fullFitAlphaOnly = fullFitAlphaOnly,
        includeAllCs = includeAllCs,
        mafCutoff = mafCutoff %||% 0,
        macCutoff = macCutoff %||% 0,
        imissCutoff = imissCutoff %||% 1
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

# QtlSumStats fine-mapping worker. `p` is the setMethod's captured arguments;
# it is extended via list_modify with the resolved tokens / row selection / LD
# sketch so the per-entry dispatch helpers read everything from one bundle.
# @noRd
.fmPipelineQtlSumStats <- function(
    data,
    methods,
    contexts,
    traitId,
    jointSpecification,
    addSusieInf,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    medianAbsCorr,
    fineMappingResult,
    twasWeights,
    dataDrivenPriorWeightsCutoff,
    verbose,
    trim,
    fullFit,
    fullFitAlphaOnly,
    includeAllCs,
    serFallback,
    rFinite,
    rMismatch,
    rssControl,
    keepFullFit,
    mafCutoff,
    macCutoff,
    imissCutoff
) {
    .fmAssertQcd(data)
    rt <- .fmQssResolveTokens(
        data = data,
        methods = methods,
        contexts = contexts,
        traitId = traitId,
        jointSpecification = jointSpecification,
        coverage = coverage,
        secondaryCoverage = secondaryCoverage,
        signalCutoff = signalCutoff,
        minAbsCorr = minAbsCorr,
        fineMappingResult = fineMappingResult,
        twasWeights = twasWeights,
        dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
        verbose = verbose,
        fullFit = fullFit,
        fullFitAlphaOnly = fullFitAlphaOnly,
        includeAllCs = includeAllCs,
        mafCutoff = mafCutoff,
        macCutoff = macCutoff,
        imissCutoff = imissCutoff
    )
    if (rt$exhaustedByJoint) {
        if (is.null(rt$jointResult)) {
            msg <- glue(
                "fineMappingPipeline(QtlSumStats): no joint fits produced. ",
                "Check that the jointSpecification scope intersects the ",
                "available data."
            )
            abort(msg)
        }
        return(rt$jointResult)
    }
    sel <- .fmQssSelectRows(data = data, contexts = contexts, traitId = traitId)
    split <- .fmQssSplitTokens(rt$tokens, sel)
    ldSketch <- getLdSketch(data)
    # Derived values are locals, not fields grafted onto a parameter pack.
    tokens <- rt$tokens
    methodArgs <- rt$methodArgs
    jointResult <- rt$jointResult
    studyCol <- sel$studyCol
    contextCol <- sel$contextCol
    traitCol <- sel$traitCol
    selRows <- sel$selRows
    univTokens <- split$univTokens
    mvTokens <- split$mvTokens
    rFiniteResolved <- .fmResolveRFinite(
        rFinite,
        serFallback,
        rMismatch,
        ldSketch
    )
    univOut <- if (length(univTokens) > 0L) {
        map(
            selRows,
            .fmRssEntryRows,
            studyCol = studyCol,
            contextCol = contextCol,
            traitCol = traitCol,
            univTokens = univTokens,
            fineMappingResult = fineMappingResult,
            data = data,
            ldSketch = ldSketch,
            addSusieInf = addSusieInf,
            coverage = coverage,
            secondaryCoverage = secondaryCoverage,
            signalCutoff = signalCutoff,
            minAbsCorr = minAbsCorr,
            methodArgs = methodArgs,
            verbose = verbose,
            fullFit = fullFit,
            fullFitAlphaOnly = fullFitAlphaOnly,
            includeAllCs = includeAllCs,
            serFallback = serFallback,
            rFiniteResolved = rFiniteResolved,
            rMismatch = rMismatch,
            rssControl = rssControl,
            keepFullFit = keepFullFit,
            mafCutoff = mafCutoff,
            macCutoff = macCutoff,
            imissCutoff = imissCutoff
        )
    } else {
        list()
    }
    rows <- list_flatten(map(univOut, "rows"))
    nSkipped <- sum(map_lgl(univOut, "skipped"))
    .fmQssAssemble(
        data = data,
        ldSketch = ldSketch,
        rows = rows,
        nSkipped = nSkipped,
        jointResult = .fmQssAutoJoint(
            data = data,
            contexts = contexts,
            traitId = traitId,
            coverage = coverage,
            secondaryCoverage = secondaryCoverage,
            signalCutoff = signalCutoff,
            minAbsCorr = minAbsCorr,
            fineMappingResult = fineMappingResult,
            twasWeights = twasWeights,
            dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
            verbose = verbose,
            fullFit = fullFit,
            fullFitAlphaOnly = fullFitAlphaOnly,
            includeAllCs = includeAllCs,
            mafCutoff = mafCutoff,
            macCutoff = macCutoff,
            imissCutoff = imissCutoff,
            methodArgs = methodArgs,
            mvTokens = mvTokens,
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
        coverage = 0.95,
        secondaryCoverage = c(0.7, 0.5),
        signalCutoff = 0.025,
        minAbsCorr = 0.8,
        medianAbsCorr = NULL,
        fineMappingResult = NULL,
        twasWeights = NULL,
        dataDrivenPriorWeightsCutoff = 1e-10,
        verbose = 1,
        trim = TRUE,
        fullFit = FALSE,
        fullFitAlphaOnly = TRUE,
        includeAllCs = FALSE,
        serFallback = FALSE,
        rFinite = NULL,
        rMismatch = "none",
        rssControl = NULL,
        keepFullFit = "fallback",
        mafCutoff = 0,
        macCutoff = 0,
        imissCutoff = 1,
        ...
    ) {
        .fmPipelineQtlSumStats(
            data = data,
            methods = methods,
            contexts = contexts,
            traitId = traitId,
            jointSpecification = jointSpecification,
            addSusieInf = addSusieInf,
            coverage = coverage,
            secondaryCoverage = secondaryCoverage,
            signalCutoff = signalCutoff,
            minAbsCorr = minAbsCorr,
            medianAbsCorr = medianAbsCorr,
            fineMappingResult = fineMappingResult,
            twasWeights = twasWeights,
            dataDrivenPriorWeightsCutoff = dataDrivenPriorWeightsCutoff,
            verbose = verbose,
            trim = trim,
            fullFit = fullFit,
            fullFitAlphaOnly = fullFitAlphaOnly,
            includeAllCs = includeAllCs,
            serFallback = serFallback,
            rFinite = rFinite,
            rMismatch = rMismatch,
            rssControl = rssControl,
            keepFullFit = keepFullFit,
            mafCutoff = mafCutoff,
            macCutoff = macCutoff,
            imissCutoff = imissCutoff
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

# GwasSumStats fine-mapping worker. `p` is the setMethod's captured arguments
# resolved tokens / LD sketch / finite-sample size are ordinary locals; the
# tokens / LD sketch / finite-sample size so the per-entry dispatch helpers
# read everything from the single bundle. One GwasSumStats is one LD
# block (the caller builds one collection per block when sweeping the genome);
# we fine-map each (study, method) tuple across the whole entry, no in-pipeline
# block partitioning.
# @noRd
.fmPipelineGwas <- function(
    data,
    methods,
    L,
    Lgreedy,
    addSusieInf,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    fineMappingResult,
    verbose,
    fullFit,
    fullFitAlphaOnly,
    includeAllCs,
    serFallback,
    rFinite,
    rMismatch,
    rssControl,
    keepFullFit,
    mafCutoff,
    macCutoff,
    imissCutoff
) {
    .fmAssertQcd(data)
    norm <- .fmNormalizeMethods(methods, L = L, Lgreedy = Lgreedy)
    .fmCheckMethodCapabilities(norm$tokens, "GwasSumStats")
    ldSketch <- getLdSketch(data)
    # Derived values are ordinary locals now, not fields grafted onto a bundle.
    tokens <- norm$tokens
    methodArgs <- norm$methodArgs
    rFiniteResolved <- .fmResolveRFinite(
        rFinite,
        serFallback,
        rMismatch,
        ldSketch
    )
    studyCol <- as.character(data$study)
    entryOut <- map(
        seq_len(nrow(data)),
        .fmGwasEntryRows,
        data = data,
        studyCol = studyCol,
        tokens = tokens,
        ldSketch = ldSketch,
        verbose = verbose,
        fineMappingResult = fineMappingResult,
        mafCutoff = mafCutoff,
        macCutoff = macCutoff,
        imissCutoff = imissCutoff,
        addSusieInf = addSusieInf,
        coverage = coverage,
        secondaryCoverage = secondaryCoverage,
        signalCutoff = signalCutoff,
        minAbsCorr = minAbsCorr,
        methodArgs = methodArgs,
        fullFit = fullFit,
        fullFitAlphaOnly = fullFitAlphaOnly,
        includeAllCs = includeAllCs,
        serFallback = serFallback,
        rFiniteResolved = rFiniteResolved,
        rMismatch = rMismatch,
        rssControl = rssControl,
        keepFullFit = keepFullFit
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
        L = 10L,
        Lgreedy = NULL,
        coverage = 0.95,
        secondaryCoverage = c(0.7, 0.5),
        signalCutoff = 0.025,
        minAbsCorr = 0.8,
        medianAbsCorr = NULL,
        fineMappingResult = NULL,
        verbose = 1,
        trim = TRUE,
        fullFit = FALSE,
        fullFitAlphaOnly = TRUE,
        includeAllCs = FALSE,
        serFallback = FALSE,
        rFinite = NULL,
        rMismatch = "none",
        rssControl = NULL,
        keepFullFit = "fallback",
        mafCutoff = 0,
        macCutoff = 0,
        imissCutoff = 1,
        ...
    ) {
        .fmPipelineGwas(
            data = data,
            methods = methods,
            L = L,
            Lgreedy = Lgreedy,
            addSusieInf = addSusieInf,
            coverage = coverage,
            secondaryCoverage = secondaryCoverage,
            signalCutoff = signalCutoff,
            minAbsCorr = minAbsCorr,
            fineMappingResult = fineMappingResult,
            verbose = verbose,
            fullFit = fullFit,
            fullFitAlphaOnly = fullFitAlphaOnly,
            includeAllCs = includeAllCs,
            serFallback = serFallback,
            rFinite = rFinite,
            rMismatch = rMismatch,
            rssControl = rssControl,
            keepFullFit = keepFullFit,
            mafCutoff = mafCutoff,
            macCutoff = macCutoff,
            imissCutoff = imissCutoff
        )
    }
)


# =============================================================================
# ANY fallback
# =============================================================================

#' @rdname fineMappingPipeline
#' @export
setMethod("fineMappingPipeline", "ANY", function(data, ...) {
    cls <- class(data)[[1L]]
    msg <- glue(
        "fineMappingPipeline does not accept inputs of class '{cls}'. ",
        "Pass a QtlDataset, MultiStudyQtlDataset, QtlSumStats, or ",
        "GwasSumStats. Use summaryStatsQc() on SumStats inputs first."
    )
    abort(msg)
})

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
.fmUnivContextRows <- function(
    ctx,
    data,
    study,
    univTokens,
    xRegions,
    naAction,
    fineMappingResult,
    cisWindow,
    screen,
    addSusieInf,
    verbose,
    coverage,
    secondaryCoverage,
    signalCutoff,
    minAbsCorr,
    methodArgs,
    cvFolds,
    cvThreads,
    samplePartition,
    seed,
    fullFit,
    fullFitAlphaOnly,
    includeAllCs,
    perCtxTraits
) {
    list_flatten(map(
        perCtxTraits[[ctx]],
        .fmUnivTraitRows,
        ctx = ctx,
        data = data,
        study = study,
        univTokens = univTokens,
        xRegions = xRegions,
        naAction = naAction,
        fineMappingResult = fineMappingResult,
        cisWindow = cisWindow,
        screen = screen,
        addSusieInf = addSusieInf,
        verbose = verbose,
        coverage = coverage,
        secondaryCoverage = secondaryCoverage,
        signalCutoff = signalCutoff,
        minAbsCorr = minAbsCorr,
        methodArgs = methodArgs,
        cvFolds = cvFolds,
        cvThreads = cvThreads,
        samplePartition = samplePartition,
        seed = seed,
        fullFit = fullFit,
        fullFitAlphaOnly = fullFitAlphaOnly,
        includeAllCs = includeAllCs
    ))
}
