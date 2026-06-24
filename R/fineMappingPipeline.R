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
#' @section Chained initialisation:
#' When \code{susieInf} is requested alongside \code{susie} and/or
#' \code{susieAsh} and \code{addSusieInf = TRUE} (default), the
#' SuSiE-inf fit is computed first and used as initialisation for
#' the SuSiE / SuSiE-ash fits, mirroring the legacy
#' \code{univariateAnalysisPipeline} / \code{susieRssPipeline} chained
#' init behaviour. SuSiE-inf is dropped from the final result when
#' the caller did not explicitly request it (only used as init).
#'
#' @section QC contract:
#' All \code{QtlSumStats} and \code{GwasSumStats} inputs must have
#' been QC'd via \code{\link{summaryStatsQc}}; the pipeline errors on
#' inputs where \code{length(getQcInfo(x)) == 0L}.
#' \code{summaryStatsQc} also drops variants absent from the
#' \code{ldSketch}, so by the time per-entry processing runs every
#' variant is guaranteed to be present in the LD panel and a local LD
#' matrix can be built with \code{extractBlockGenotypes} +
#' \code{computeLd("sample")}.
#'
#' @section Optional resume cache:
#' Supplying \code{fineMappingResult} of an existing
#' \code{FineMappingResult} skips re-fitting any
#' \code{(study, context, trait, method)} tuple that already has a
#' matching row; cached entries are merged with the newly-fit entries
#' in the returned collection.
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
#' @param region Optional \code{GRanges} for QtlDataset trait
#'   selection. Mutually exclusive with \code{traitId}.
#' @param cisWindow For QtlDataset: cis-window (bp) around each trait's
#'   genomic position when extracting variants. Required when
#'   \code{traitId} is supplied. Mutually exclusive with \code{region}.
#' @param jointRegions For QtlDataset with a multi-range \code{region}:
#'   \code{FALSE} (default) fits each range independently and merges the
#'   per-range results into one entry per (study, context, trait, method) —
#'   the merged \code{susieFit} is a named list of per-region fits and
#'   credible-set labels are renumbered to stay unique. \code{TRUE}
#'   concatenates the ranges' genotypes into one joint fit. Ignored for a
#'   single-range / cis (\code{traitId} + \code{cisWindow}) request.
#' @param addSusieInf Logical. When \code{susieInf} is in
#'   \code{methods} alongside \code{susie} and/or \code{susieAsh},
#'   controls whether the SuSiE-inf fit initialises the chained
#'   downstream method(s). Default \code{TRUE}.
#' @param coverage Primary credible-set coverage (numeric, length 1).
#'   Default \code{0.95}.
#' @param secondaryCoverage Secondary coverages forwarded to
#'   \code{postprocessFinemappingFits}. Default \code{c(0.7, 0.5)}.
#' @param signalCutoff PIP cutoff for top-loci selection. Default
#'   \code{0.025}.
#' @param minAbsCorr Minimum absolute correlation for credible-set
#'   purity. Default \code{0.8}.
#' @param medianAbsCorr Optional median absolute correlation for
#'   credible-set purity, routed to \code{susieR::susie_get_cs}. A set is
#'   kept if it passes either \code{minAbsCorr} or \code{medianAbsCorr}
#'   (OR-logic). Default \code{NULL} (off).
#' @param fineMappingResult Optional existing \code{FineMappingResult}
#'   to use as a resume cache; tuples already present are not refit.
#' @param cvFolds Integer. Number of cross-validation folds. Default
#'   \code{0} (no CV). When \code{> 1}, each method is refit on the
#'   training samples of every fold and used to predict the held-out
#'   samples; the fold partition plus per-fold out-of-fold predictions and
#'   metrics are stored on each \code{FineMappingEntry}'s \code{cvResult}
#'   slot (see \code{\link{getCvResult}}). \code{twasWeightsPipeline} reuses
#'   this partition and feeds these predictions into the SR-TWAS ensemble.
#'   Individual-level (\code{QtlDataset} / \code{MultiStudyQtlDataset})
#'   input only; ignored for sumstat inputs.
#' @param samplePartition Optional pre-defined CV partition
#'   \code{data.frame} with columns \code{Sample} and \code{Fold}. When
#'   supplied (and \code{cvFolds > 1}), every method reuses this exact
#'   partition; otherwise a fresh partition is generated per
#'   \code{(study, context, trait)}.
#' @param seed Optional integer. When non-NULL, \code{set.seed(seed)} is
#'   called once at the start of the call for reproducible fits. Default
#'   \code{NULL} (no seeding).
#' @param pipCutoffToSkip Numeric (length 1). Individual-level single-effect
#'   (SER) pre-screen applied to each residualized \code{(X, y)} block before a
#'   full fit: a susie model with \code{L = 1} is fit and the block is skipped
#'   when no PIP exceeds the cutoff (no potentially significant variant). The
#'   summary-statistics analog lives in \code{summaryStatsQc()}. \code{0}
#'   (default) disables the screen; a negative value uses the adaptive
#'   \code{3 / nVariants} threshold.
#' @param jointSpecification Optional joint-fit specification (NULL by
#'   default). When NULL, the pipeline runs the implicit multi-context /
#'   multi-trait mvSuSiE / fSuSiE branches as before. When non-NULL, the
#'   argument is parsed and validated via the joint-spec grammar
#'   documented under \code{parseJointSpecification} (a character vector
#'   of axes, or a list of \code{list(axes, scope)} specs); the
#'   per-spec axis dispatcher implementation is in progress and a
#'   non-NULL value currently errors with an informative message.
#'   See the design notes in \code{R/jointSpecification.R} for the
#'   accepted grammar.
#' @param ldBlocks For \code{GwasSumStats} input only: an
#'   \code{LdBlocks} object describing the LD-block partition. The
#'   pipeline performs SuSiE-RSS fine-mapping per (study, ldBlock).
#'   Required for the GwasSumStats method.
#' @param verbose Verbosity (0 silent, 1 default). Default \code{1}.
#' @param phenotypeCovariatesToResidualize Character vector (or
#'   \code{NULL}) of phenotype-covariate names to residualize against.
#'   \code{NULL} (default) uses every available phenotype covariate.
#'   Only meaningful when the input is a \code{QtlDataset} /
#'   \code{MultiStudyQtlDataset} (ignored for sumstat inputs).
#' @param genotypeCovariatesToResidualize Character vector (or
#'   \code{NULL}) of genotype-covariate column names to residualize
#'   against. \code{NULL} uses every available genotype covariate.
#' @param residualizePhenotypeCovariates Logical (length 1). When
#'   \code{TRUE} (default) residualize against the phenotype-side
#'   covariates listed in \code{phenotypeCovariatesToResidualize}. Set
#'   \code{FALSE} to disable phenotype-covariate residualization
#'   entirely. The marginal univariate effects stored on each
#'   \code{FineMappingEntry} obey the same residualization choice as
#'   the SuSiE fit itself — they are computed against the same
#'   residualized \code{X} / \code{Y}.
#' @param residualizeGenotypeCovariates Logical (length 1). When
#'   \code{TRUE} (default) residualize against the genotype-side
#'   covariates listed in \code{genotypeCovariatesToResidualize}. Set
#'   \code{FALSE} to disable.
#' @param trim Logical (length 1). When \code{TRUE} (default) the
#'   \code{susieFit} slot on each output \code{FineMappingEntry} carries
#'   a trimmed view of the SuSiE fit (the minimal subset needed by
#'   downstream pipelines). When \code{FALSE} the full untrimmed
#'   \code{susie()} return is retained so accessors like
#'   \code{getSusieFit()} and non-default-coverage queries through
#'   \code{getCs()} can read the full posterior matrices
#'   (\code{lbf_variable}, \code{mu}, \code{mu2}, \code{V}). The
#'   per-variant \code{topLoci} table is always fully populated
#'   regardless of \code{trim}.
#' @param ... Reserved for future per-method arguments.
#'
#' @return A \code{\link{FineMappingResult}} collection keyed by
#'   \code{(study, context, trait, method)}. The \code{ldSketch} slot
#'   is set automatically: \code{NULL} for individual-level
#'   (QtlDataset / all-individual-level MultiStudyQtlDataset) fits, the
#'   input's \code{ldSketch} for RSS-derived fits.
#' @export
setGeneric("fineMappingPipeline",
  function(data, ...) standardGeneric("fineMappingPipeline"))


# =============================================================================
# Method capability table — unified naming, individual vs sumstat dispatch
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
# `mrmash` is intentionally listed with both impls NULL so the capability
# checker emits a clear rejection ("mr.mash is a TWAS-weight-oriented
# method — use twasWeightsPipeline()").
#
# @noRd
.fineMappingMethodCapabilities <- list(
  susie = list(
    individualImpl    = "susieR::susie",
    sumstatImpl       = "susieR::susie_rss",
    multivariate      = FALSE,
    gwasAllowed       = TRUE,
    unmappableEffects = "none",
    args              = list()),
  susieInf = list(
    individualImpl    = "susieR::susie",
    sumstatImpl       = "susieR::susie_rss",
    multivariate      = FALSE,
    gwasAllowed       = TRUE,
    unmappableEffects = "inf",
    args              = list()),
  susieAsh = list(
    individualImpl    = "susieR::susie",
    sumstatImpl       = "susieR::susie_rss",
    multivariate      = FALSE,
    gwasAllowed       = TRUE,
    unmappableEffects = "ash",
    args              = list()),
  mvsusie = list(
    individualImpl    = "mvsusieR::mvsusie",
    sumstatImpl       = "mvsusieR::mvsusie_rss",
    multivariate      = TRUE,
    gwasAllowed       = FALSE,
    unmappableEffects = NA_character_,
    args              = list()),
  fsusie = list(
    individualImpl    = "fsusieR::susiF",
    sumstatImpl       = NULL,
    multivariate      = TRUE,
    gwasAllowed       = FALSE,
    unmappableEffects = NA_character_,
    args              = list()),
  mrmash = list(
    individualImpl    = NULL,
    sumstatImpl       = NULL,
    multivariate      = TRUE,
    gwasAllowed       = FALSE,
    unmappableEffects = NA_character_,
    args              = list()))


# Normalize a user-supplied `methods` argument into a character vector of
# canonical tokens. Mirrors `.twasNormalizeMethods` but the fine-mapping
# pipeline takes only a character vector (no preset strings, no list form).
# @noRd
# Normalize a user-supplied `methods` argument into `(tokens, methodArgs)`.
#
# Accepts:
#   * character vector  c("susie", "susieInf")            -> empty kwargs per token
#   * named list        list(susie = list(L = 1), ...)    -> per-token kwargs
#
# Names of the returned `methodArgs` always equal `tokens` (one entry per
# token, empty list when the user supplied none). The fitters then
# `modifyList`-merge each entry into the base arg list before do.call.
#
# Mirrors the convention of .twasNormalizeMethods so the two pipelines
# expose the same shape on the user side.
# @noRd
.fmNormalizeMethods <- function(methods) {
  if (is.null(methods) || length(methods) == 0L) {
    stop("fineMappingPipeline: `methods` must be a non-empty character ",
         "vector or named list of <token> = <kwargs> entries.")
  }
  if (is.character(methods)) {
    tokens     <- unique(methods)
    methodArgs <- setNames(rep(list(list()), length(tokens)), tokens)
  } else if (is.list(methods)) {
    if (is.null(names(methods)) || any(names(methods) == "")) {
      stop("fineMappingPipeline: when `methods` is a list it must be ",
           "named (one entry per method token).")
    }
    nonListChild <- vapply(methods, function(x) !is.list(x), logical(1))
    if (any(nonListChild)) {
      stop("fineMappingPipeline: each entry of the `methods` list must ",
           "itself be a list of named kwargs (got non-list value for: ",
           paste(names(methods)[nonListChild], collapse = ", "), ").")
    }
    tokens     <- unique(names(methods))
    methodArgs <- methods[tokens]
  } else {
    stop("fineMappingPipeline: `methods` must be a character vector or ",
         "named list. Got class '", class(methods)[[1L]], "'.")
  }
  list(tokens = tokens, methodArgs = methodArgs)
}


# Enforce input-class / method compatibility against the fine-mapping
# capability table. Hard-rejects `mrmash` (a TWAS-weight-oriented
# method). Routes the input class through individual / sumstat / GWAS
# branches and emits a single error listing every offending token.
# @noRd
.fmCheckMethodCapabilities <- function(tokens, inputKind) {
  if (length(tokens) == 0L) return(invisible(NULL))
  caps <- .fineMappingMethodCapabilities
  unknown <- setdiff(tokens, names(caps))
  if (length(unknown) > 0L) {
    stop(sprintf(
      "fineMappingPipeline: unknown method token(s): %s. Known tokens: %s.",
      paste(unknown, collapse = ", "),
      paste(names(caps), collapse = ", ")))
  }
  hardRejections <- list(
    mrmash = "mr.mash is a TWAS-weight-oriented method; use twasWeightsPipeline()")
  individualKinds <- c("QtlDataset", "MultiStudyQtlDataset")
  bad <- character(0); reason <- character(0)
  for (tk in tokens) {
    info <- caps[[tk]]
    if (tk %in% names(hardRejections)) {
      bad <- c(bad, tk); reason <- c(reason, hardRejections[[tk]])
      next
    }
    if (inputKind %in% individualKinds) {
      if (is.null(info$individualImpl)) {
        bad <- c(bad, tk)
        reason <- c(reason, "is sumstat-only on this pipeline")
      }
    } else if (inputKind == "QtlSumStats") {
      if (is.null(info$sumstatImpl)) {
        bad <- c(bad, tk)
        reason <- c(reason, "is individual-only (use a QtlDataset input)")
      }
    } else if (inputKind == "GwasSumStats") {
      if (!isTRUE(info$gwasAllowed) || is.null(info$sumstatImpl)) {
        bad <- c(bad, tk)
        reason <- c(reason,
          "is not supported on GwasSumStats (only the SuSiE-RSS family is)")
      }
    }
  }
  if (length(bad) > 0L) {
    stop(sprintf(
      "fineMappingPipeline: the following method(s) are not available for input class '%s': %s. %s.",
      inputKind,
      paste(unique(bad), collapse = ", "),
      paste(sprintf("%s %s", bad, reason), collapse = "; ")))
  }
}


# Reject SumStats inputs that have not been QC'd via summaryStatsQc.
# @noRd
.fmAssertQcd <- function(sumstats) {
  if (length(getQcInfo(sumstats)) == 0L) {
    stop("fineMappingPipeline: the supplied ",
         class(sumstats)[[1L]],
         " has no QC record (qcInfo is empty). Call summaryStatsQc() ",
         "first and pass the QC-applied result.")
  }
}


# Given a `methods` vector, decide whether the SuSiE-inf chained-init
# shortcut applies. Returns a list of (chainSusie, chainAsh, runInf,
# keepInf): runInf is TRUE when susieInf must be fitted (either user
# requested it OR a chained init needs it); keepInf is TRUE when the
# user asked for "susieInf" in `methods` directly.
# @noRd
.fmResolveSusieChain <- function(tokens, addSusieInf) {
  hasInf <- "susieInf" %in% tokens
  hasSu  <- "susie"    %in% tokens
  hasAsh <- "susieAsh" %in% tokens
  chainSusie <- isTRUE(addSusieInf) && hasInf && hasSu
  chainAsh   <- isTRUE(addSusieInf) && hasInf && hasAsh
  runInf     <- hasInf || chainSusie || chainAsh
  keepInf    <- hasInf
  list(chainSusie = chainSusie, chainAsh = chainAsh,
       runInf = runInf, keepInf = keepInf)
}


# Optional resume-cache lookup. Returns the matching FineMappingEntry from
# `fineMappingResult` for the tuple (study, context, trait, method), or
# NULL when there is no hit. Returns NULL silently when fineMappingResult
# is NULL or not a QtlFineMappingResult.
# @noRd
.fmCacheLookup <- function(fineMappingResult, study, context, trait, method) {
  if (is.null(fineMappingResult)) return(NULL)
  if (!is(fineMappingResult, "QtlFineMappingResult")) return(NULL)
  idx <- .matchTupleRows(fineMappingResult,
                          list(study = study, context = context,
                               trait = trait, method = method))
  if (length(idx) == 0L) return(NULL)
  fineMappingResult$entry[[idx[[1L]]]]
}

# GwasFineMappingResult cache lookup using the (study, method,
# region_id) 3-tuple. Multi-block FMRs can carry one entry per
# (study, method, region_id) triple, so the cache key must include
# region_id to correctly retrieve the cached fit for a specific block.
# @noRd
.fmCacheLookupGwas <- function(fineMappingResult, study, method, region_id) {
  if (is.null(fineMappingResult)) return(NULL)
  if (!is(fineMappingResult, "GwasFineMappingResult")) return(NULL)
  idx <- .matchTupleRows(fineMappingResult,
                          list(study = study, method = method,
                               region_id = region_id))
  if (length(idx) == 0L) return(NULL)
  fineMappingResult$entry[[idx[[1L]]]]
}


# Build a QtlFineMappingResult collection from per-tuple parallel vectors.
# `jointStudies`, `jointContexts`, `jointTraits` are optional character
# vectors (length matches `studies`) describing semicolon-joined joint
# members for cross-study / cross-context / cross-trait joint fits; pass
# `NULL` (default) to omit the column entirely.
# @noRd
.fmBuildQtlResult <- function(studies, contexts, traits, methods, entries,
                              jointStudies  = NULL,
                              jointContexts = NULL,
                              jointTraits   = NULL,
                              ldSketch      = NULL) {
  if (length(entries) == 0L) {
    stop("fineMappingPipeline: no (study, context, trait, method) tuples ",
         "produced a fine-mapping result.")
  }
  QtlFineMappingResult(
    study         = studies,
    context       = contexts,
    trait         = traits,
    method        = methods,
    entry         = entries,
    jointStudies  = jointStudies,
    jointContexts = jointContexts,
    jointTraits   = jointTraits,
    ldSketch      = ldSketch)
}

# Build a GwasFineMappingResult collection from per-row vectors. When
# `region_ids` is NULL, falls through to the constructor's synthetic
# defaults (region_1, region_2, ...). For callers that have meaningful
# block labels (e.g. derived from a GwasSumStats entry's GRanges),
# pass them explicitly so downstream consumers can join on region.
# @noRd
.fmBuildGwasResult <- function(studies, methods, entries,
                               region_ids = NULL, ldSketch = NULL) {
  if (length(entries) == 0L) {
    stop("fineMappingPipeline: no (study, method, region_id) tuples produced a ",
         "fine-mapping result.")
  }
  GwasFineMappingResult(
    study     = studies,
    method    = methods,
    region_id = region_ids,
    entry     = entries,
    ldSketch  = ldSketch)
}

# Combine an optional joint column across two collections. Returns NULL
# when neither input carries the column (so the rebuilt collection
# omits it too); otherwise pads the missing side with NA_character_ so
# both halves contribute a same-length character vector.
# @noRd
.combineJointCol <- function(a, b, colName) {
  hasA <- colName %in% names(a)
  hasB <- colName %in% names(b)
  if (!hasA && !hasB) return(NULL)
  aVals <- if (hasA) as.character(a[[colName]])
           else rep(NA_character_, nrow(a))
  bVals <- if (hasB) as.character(b[[colName]])
           else rep(NA_character_, nrow(b))
  c(aVals, bVals)
}

# Concatenate two FineMappingResult collections row-wise. Routes to the
# right constructor based on input class. rbind() on a DFrame subclass
# does not reliably preserve the ldSketch slot, so this helper rebuilds
# the collection via the constructor.
# @noRd
.rbindFineMappingResult <- function(a, b, ldSketch = NULL) {
  if (!is(a, "FineMappingResultBase") || !is(b, "FineMappingResultBase")) {
    stop(".rbindFineMappingResult expects two FineMappingResultBase inputs.")
  }
  if (!identical(class(a)[[1L]], class(b)[[1L]])) {
    stop(".rbindFineMappingResult: inputs must be the same concrete class ",
         "(got '", class(a)[[1L]], "' and '", class(b)[[1L]], "').")
  }
  if (is(a, "QtlFineMappingResult")) {
    QtlFineMappingResult(
      study         = c(as.character(a$study),   as.character(b$study)),
      context       = c(as.character(a$context), as.character(b$context)),
      trait         = c(as.character(a$trait),   as.character(b$trait)),
      method        = c(as.character(a$method),  as.character(b$method)),
      entry         = c(as.list(a$entry), as.list(b$entry)),
      jointStudies  = .combineJointCol(a, b, "jointStudies"),
      jointContexts = .combineJointCol(a, b, "jointContexts"),
      jointTraits   = .combineJointCol(a, b, "jointTraits"),
      ldSketch      = ldSketch)
  } else {
    GwasFineMappingResult(
      study     = c(as.character(a$study),     as.character(b$study)),
      method    = c(as.character(a$method),    as.character(b$method)),
      region_id = c(as.character(a$region_id), as.character(b$region_id)),
      entry     = c(as.list(a$entry), as.list(b$entry)),
      ldSketch  = ldSketch)
  }
}


# Build an LD correlation matrix from an LD sketch genotype handle for a
# specific variant subset. Thin wrapper over the shared `.ldFromSketch`
# helper (R/ld.R).
# @noRd
.fmLdFromSketch <- function(ldSketch, variantIds) {
  .ldFromSketch(ldSketch, variantIds, label = ".fmLdFromSketch")
}


# Wrap one finemapping fit into a FineMappingEntry via the surviving
# post-processing helpers (postprocessFinemappingFits +
# formatFinemappingOutput). Returns a bare FineMappingEntry payload, ready
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
  "residualizeGenotypeCovariates")

.resPickFlags <- function() {
  out <- list()
  # Walk up from the immediate caller; the public setMethod frame is
  # where the user-facing args live. sys.frames()[[1]] is the global
  # env so stop before that.
  frames <- sys.frames()
  for (i in seq_along(frames)) {
    fr <- frames[[i]]
    for (nm in .resFlagNames) {
      if (!nm %in% names(out) && exists(nm, envir = fr, inherits = FALSE)) {
        out[[nm]] <- get(nm, envir = fr, inherits = FALSE)
      }
    }
  }
  out
}

.fmResidPheno <- function(x, ...) {
  do.call(getResidualizedPhenotypes,
          c(list(x = x, ...), .resPickFlags()))
}

.fmResidGeno <- function(x, ...) {
  do.call(getResidualizedGenotypes,
          c(list(x = x, ...), .resPickFlags()))
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
.fmAfForX <- function(data, X, traitId = NULL, region = NULL,
                      cisWindow = NULL) {
  if (is.null(X) || ncol(X) == 0L || nrow(X) == 0L) return(NULL)
  if (!is(data, "QtlDataset")) return(NULL)
  afAll <- tryCatch(
    if (is.null(region)) {
      getAf(data, traitId = traitId, cisWindow = cisWindow,
            samples = rownames(X))
    } else {
      getAf(data, region = region, samples = rownames(X))
    },
    error = function(e) NULL)
  if (is.null(afAll) || length(afAll) == 0L) return(NULL)
  unname(afAll[colnames(X)])
}

.fmPostprocessOne <- function(fit, method, dataX, dataY,
                              coverage, secondaryCoverage, signalCutoff,
                              minAbsCorr, csInput = NULL, af = NULL,
                              region = NULL, trim = NULL,
                              medianAbsCorr = NULL) {
  # Inherit `trim` from the calling method's frame if not passed in
  # explicitly. The 10 internal call sites don't currently forward it
  # (they predate the trim knob) so we look it up from the caller. This
  # keeps the patch surface minimal: each public setMethod gains a
  # `trim = TRUE` parameter and that value naturally reaches here.
  if (is.null(trim)) {
    trim <- tryCatch(get("trim", envir = parent.frame()),
                     error = function(e) TRUE)
  }
  # `medianAbsCorr` is inherited the same way (each public setMethod gains a
  # `medianAbsCorr = NULL` parameter); NULL is a no-op (OR-logic purity off).
  if (is.null(medianAbsCorr)) {
    medianAbsCorr <- tryCatch(get("medianAbsCorr", envir = parent.frame()),
                              error = function(e) NULL)
  }
  fits <- setNames(list(fit), method)
  post <- postprocessFinemappingFits(
    fits = fits, dataX = dataX, dataY = dataY,
    af = af, coverage = coverage,
    secondaryCoverage = secondaryCoverage,
    signalCutoff = signalCutoff, minAbsCorr = minAbsCorr,
    medianAbsCorr = medianAbsCorr,
    region = region,
    csInput = csInput, trim = isTRUE(trim))
  out <- formatFinemappingOutput(post, primaryMethod = method)
  # `formatFinemappingOutput` returns a list with $finemappingEntry as a
  # bare FineMappingEntry per the helper's contract.
  if (!is(out$finemappingEntry, "FineMappingEntry")) {
    stop(".fmPostprocessOne: postprocess output did not carry a ",
         "FineMappingEntry payload — check pecotmr internal contract.")
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
  if (is.null(region)) {
    list(NULL)
  } else if (isTRUE(jointRegions)) {
    list(region)
  } else {
    lapply(seq_along(region), function(i) region[i])
  }
}

# Single-effect (SER) pre-screen, individual-level. Fits susie with L = 1 on a
# residualized (X, y) block and reports whether any PIP clears `cutoff` -- i.e.
# whether the block shows any potentially significant variant worth a full fit.
# Ports the deleted multivariate_pipeline.R `skipConditions` / susie_twas
# `pip_cutoff_to_skip` logic (the individual-level analog of the sumstat-path
# `.applyPipScreen`):
#   * `cutoff == 0` (or NULL/non-scalar) disables the screen -> always keep.
#   * `cutoff < 0` uses the adaptive 3 / nVariants threshold.
#   * NA entries of `y` are dropped before fitting.
# The screen is advisory: too few samples/variants or a fit failure keeps the
# block (returns TRUE) rather than discarding a potentially real signal.
# @noRd
.fmSerScreen <- function(X, y, cutoff) {
  if (is.null(cutoff) || length(cutoff) != 1L || is.na(cutoff) || cutoff == 0)
    return(TRUE)
  ok <- !is.na(y)
  if (sum(ok) < 2L || ncol(X) < 1L) return(TRUE)
  Xs <- X[ok, , drop = FALSE]
  ys <- y[ok]
  thr <- if (cutoff < 0) 3 / ncol(Xs) else cutoff
  pip <- tryCatch(suppressMessages(susieR::susie(Xs, ys, L = 1L))$pip,
                  error = function(e) NULL)
  if (is.null(pip)) return(TRUE)
  any(pip > thr)
}

# Is the SER pre-screen enabled? Only a finite, non-zero scalar activates it;
# this gates the extra screening extraction so the default (cutoff 0) costs
# nothing.
# @noRd
.fmScreenActive <- function(cutoff) {
  !is.null(cutoff) && length(cutoff) == 1L && !is.na(cutoff) && cutoff != 0
}

# Per-condition SER pre-screen for a joint (multi-context / multi-trait) fit:
# returns a logical vector over the columns of `Y` (the conditions) marking
# which show single-effect signal. The multivariate analog of `.fmSerScreen`
# and a port of the deleted `skipConditions`: callers drop the FALSE columns
# (null contexts / traits) before the joint mvSuSiE fit.
# @noRd
.fmSerScreenColumns <- function(X, Y, cutoff) {
  vapply(seq_len(ncol(Y)),
         function(j) .fmSerScreen(X, Y[, j], cutoff),
         logical(1L))
}

# Fit every requested univariate token on one residualized (X, y) block,
# returning a named list (token -> FineMappingEntry). Extracted from the
# univariate dispatch so the same logic serves the cis path (one block), the
# jointRegions=TRUE path (one concatenated block) and the jointRegions=FALSE
# path (one block per region, merged afterwards via .fmMergeEntries).
.fmFitXBlock <- function(X, y, toRun, addSusieInf, coverage,
                         secondaryCoverage, signalCutoff, minAbsCorr,
                         methodArgs, verbose, ctx, tid,
                         cvFolds = 0, samplePartition = NULL, af = NULL) {
  chainLocal <- .fmResolveSusieChain(toRun, addSusieInf)
  infFit <- NULL
  if (chainLocal$runInf) {
    if (verbose >= 1)
      message(sprintf("Fitting susieInf for (context='%s', trait='%s') ...",
                      ctx, tid))
    infFit <- .fmFitSusieIndiv(X, y, "susieInf", coverage = coverage,
                               userArgs = methodArgs[["susieInf"]])
  }
  out <- list()
  for (tk in toRun) {
    if (tk == "susieInf") {
      if (!chainLocal$keepInf) next
      fit <- infFit
    } else {
      chainFrom <- if ((tk == "susie"    && chainLocal$chainSusie) ||
                       (tk == "susieAsh" && chainLocal$chainAsh))
                     infFit else NULL
      if (verbose >= 1)
        message(sprintf("Fitting %s for (context='%s', trait='%s') ...",
                        tk, ctx, tid))
      fit <- .fmFitSusieIndiv(X, y, tk, chainFromInf = chainFrom,
                              coverage = coverage, userArgs = methodArgs[[tk]])
    }
    out[[tk]] <- .fmPostprocessOne(
      fit = fit, method = tk, dataX = X, dataY = y, coverage = coverage,
      secondaryCoverage = secondaryCoverage, signalCutoff = signalCutoff,
      minAbsCorr = minAbsCorr, af = af, csInput = "X")
  }
  # Per-fold cross-validation across the fitted univariate methods; attach
  # each method's out-of-fold predictions to its entry.
  if (cvFolds > 1L && length(out) > 0L) {
    if (verbose >= 1)
      message(sprintf("Cross-validating (%d folds) for (context='%s', trait='%s') ...",
                      cvFolds, ctx, tid))
    cv <- .fmCrossValidate(X, y, names(out), methodArgs, cvFolds,
                           samplePartition = samplePartition,
                           coverage = coverage, verbose = verbose)
    for (tk in names(out)) {
      out[[tk]] <- .fmAttachCv(out[[tk]], .fmSliceCv(cv, tk))
    }
  }
  out
}

# Extract integer credible-set indices from a "<method>_<idx>" vector.
.fmCsIdx <- function(csVec) {
  suppressWarnings(as.integer(sub("^.*_([0-9]+)$", "\\1", as.character(csVec))))
}

# Re-number credible-set membership labels by `offset`, preserving the
# "<method>_0" (not-in-any-CS) sentinel.
.fmRelabelCs <- function(csVec, offset) {
  csVec <- as.character(csVec)
  if (offset == 0L) return(csVec)
  parts <- regmatches(csVec, regexec("^(.*)_([0-9]+)$", csVec))
  vapply(seq_along(csVec), function(j) {
    p <- parts[[j]]
    if (length(p) != 3L) return(csVec[[j]])
    idx <- as.integer(p[[3L]])
    if (idx == 0L) csVec[[j]] else paste0(p[[2L]], "_", idx + offset)
  }, character(1))
}

# Merge per-region FineMappingEntry payloads (same study/context/trait/method,
# independent fits) into one entry: concatenate variants and topLoci rows,
# renumber credible sets so per-region indices do not collide, and keep the
# per-region SuSiE fits as a named list in `susieFit` (consumers needing a
# single fit must iterate the list).
.fmMergeEntries <- function(entries) {
  entries <- entries[!vapply(entries, is.null, logical(1))]
  if (length(entries) == 0L) return(NULL)
  if (length(entries) == 1L) return(entries[[1L]])
  variantIds <- unlist(lapply(entries, function(e) e@variantIds),
                       use.names = FALSE)
  tls <- lapply(entries, function(e) e@topLoci)
  csCols <- grep("^cs_[0-9]+$",
                 unique(unlist(lapply(tls, names))), value = TRUE)
  offsets <- setNames(integer(length(csCols)), csCols)
  for (i in seq_along(tls)) {
    tl <- tls[[i]]
    for (cc in csCols) {
      if (!cc %in% names(tl)) next
      idx <- .fmCsIdx(tl[[cc]])
      tl[[cc]] <- .fmRelabelCs(tl[[cc]], offsets[[cc]])
      offsets[[cc]] <- offsets[[cc]] + max(c(0L, idx), na.rm = TRUE)
    }
    tls[[i]] <- tl
  }
  topLoci <- do.call(rbind, tls)
  rownames(topLoci) <- NULL
  susieFit <- setNames(lapply(entries, function(e) e@susieFit),
                       paste0("region", seq_along(entries)))
  # Per-region CV partitions/predictions share the same sample set; keep them
  # per region under region* names so a multi-region entry retains each block's
  # cross-validated predictions (NULL when no region carried CV).
  cvList <- setNames(lapply(entries, function(e) e@cvResult),
                     paste0("region", seq_along(entries)))
  cvResult <- if (all(vapply(cvList, is.null, logical(1)))) NULL else cvList
  FineMappingEntry(variantIds = variantIds, susieFit = susieFit,
                   topLoci = topLoci, cvResult = cvResult)
}

# Run a joint-method fit (mvsusie / fsusie) once per region block via the
# method-specific `fitOneRegion(rg)` closure (returns one FineMappingEntry per
# region), then merge across regions into a single shared entry. A single block
# (cis or jointRegions=TRUE) returns its entry unchanged.
.fmJointBlocks <- function(xRegions, fitOneRegion) {
  ents <- lapply(seq_along(xRegions), function(i) fitOneRegion(xRegions[[i]]))
  ents <- ents[!vapply(ents, is.null, logical(1))]
  if (length(ents) == 0L) return(NULL)
  if (length(ents) == 1L) ents[[1L]] else .fmMergeEntries(ents)
}


# Merge per-method user kwargs onto a base arg list. `userArgs` is the
# per-token kwargs supplied by the caller (e.g. `list(L = 1, refine =
# FALSE)`); the capability table's `args` default fills in any keys the
# user did not set. User-supplied values always win over base, capability
# defaults, and chain-derived args. Returns the merged list.
# @noRd
.fmMergeUserArgs <- function(baseArgs, token, userArgs = NULL) {
  if (is.null(userArgs)) userArgs <- list()
  info <- .fineMappingMethodCapabilities[[token]]
  capDefaults <- if (!is.null(info) && !is.null(info$args)) info$args else list()
  # Order matters: base < capability defaults < user overrides.
  if (length(capDefaults) > 0L) baseArgs <- modifyList(baseArgs, capDefaults)
  if (length(userArgs)   > 0L) baseArgs <- modifyList(baseArgs, userArgs)
  baseArgs
}


# Fit one of the SuSiE-family individual-level methods on (X, y). When
# `chainFromInf` is non-NULL, the susieInf fit it points at is used as
# initialisation (with prepareSusieFromInfArgs); otherwise a plain fit
# with the requested `unmappable_effects` is performed. `userArgs` are
# spliced via .fmMergeUserArgs (user wins over chain/base/capability
# defaults), so the caller can override things like L, max_iter,
# estimate_residual_method, refine, etc.
# @noRd
.fmFitSusieIndiv <- function(X, y, token, chainFromInf = NULL,
                             coverage = 0.95, userArgs = NULL) {
  info <- .fineMappingMethodCapabilities[[token]]
  if (is.null(info) || identical(info$unmappableEffects, NA_character_)) {
    stop(".fmFitSusieIndiv: token '", token, "' is not a SuSiE-family method.")
  }
  baseArgs <- list(X = X, y = y, coverage = coverage,
                   unmappable_effects = info$unmappableEffects)
  if (token == "susieInf") {
    baseArgs$convergence_method <- "pip"
    baseArgs$refine <- FALSE
    baseArgs$model_init <- NULL
  } else if (!is.null(chainFromInf)) {
    chainedArgs <- prepareSusieFromInfArgs(
      list(),
      chainFromInf,
      refineDefault = if (token == "susie") TRUE else NULL,
      unmappableEffects = if (token == "susieAsh") "ash" else "none")
    baseArgs <- modifyList(baseArgs, chainedArgs)
    # chainedArgs already supplies unmappable_effects + model_init; X / y
    # / coverage stay in baseArgs.
    baseArgs$X <- X; baseArgs$y <- y; baseArgs$coverage <- coverage
  } else if (token == "susieAsh") {
    baseArgs$convergence_method <- "pip"
  }
  baseArgs <- .fmMergeUserArgs(baseArgs, token, userArgs)
  fit <- do.call(susieR::susie, baseArgs)
  .setFinemappingFitClass(fit, token)
}


# Sumstat counterpart of .fmFitSusieIndiv. Calls susieR::susie_rss with
# the same unmappable_effects switch, chained init, and userArgs merge.
# @noRd
.fmFitSusieRss <- function(z, R, n, token, chainFromInf = NULL,
                           coverage = 0.95, userArgs = NULL) {
  info <- .fineMappingMethodCapabilities[[token]]
  if (is.null(info) || identical(info$unmappableEffects, NA_character_)) {
    stop(".fmFitSusieRss: token '", token, "' is not a SuSiE-family method.")
  }
  baseArgs <- list(z = z, R = R, n = n, coverage = coverage,
                   unmappable_effects = info$unmappableEffects)
  if (token == "susieInf") {
    baseArgs$convergence_method <- "pip"
    baseArgs$refine <- FALSE
    baseArgs$model_init <- NULL
  } else if (!is.null(chainFromInf)) {
    chainedArgs <- prepareSusieFromInfArgs(
      list(),
      chainFromInf,
      refineDefault = if (token == "susie") TRUE else NULL,
      unmappableEffects = if (token == "susieAsh") "ash" else "none")
    baseArgs <- modifyList(baseArgs, chainedArgs)
    baseArgs$z <- z; baseArgs$R <- R; baseArgs$n <- n
    baseArgs$coverage <- coverage
  } else if (token == "susieAsh") {
    baseArgs$convergence_method <- "pip"
  }
  baseArgs <- .fmMergeUserArgs(baseArgs, token, userArgs)
  fit <- do.call(susieR::susie_rss, baseArgs)
  # All susie_rss fits get the "susieRss" S3 class for post-processing
  # (this drives the Xcorr cs-input mode). Token-level distinction stays
  # in the `method` column of the FineMappingResult.
  .setFinemappingFitClass(fit, "susieRss")
}


# Extract variant ids + Z + (median) N from a single QtlSumStats /
# GwasSumStats entry GRanges. Errors when Z or N is missing. Wraps the
# shared `.entryToSumstatDf` helper (R/sumstatsQc.R).
# @noRd
.fmExtractZN <- function(gr, label) {
  df <- .entryToSumstatDf(gr,
                          require = c("SNP", "Z", "N"),
                          label = label)
  list(variantIds = df$variant_id,
       z = df$z,
       n = stats::median(df$N, na.rm = TRUE))
}


# =============================================================================
# Per-fold cross-validation of fine-mapping methods
# -----------------------------------------------------------------------------
# fineMappingPipeline mirrors twasWeightsPipeline's cross-validation: when
# cvFolds > 1, each fine-mapping method is refit on the training samples of
# every fold, its weights extracted and used to predict the held-out samples,
# yielding out-of-fold predictions + per-outcome metrics. The partition and
# predictions are stored on each FineMappingEntry's cvResult slot so
# twasWeightsPipeline can (a) reuse the identical fold partition and (b) feed
# fine-mapping's own cross-validated predictions straight into the SR-TWAS
# ensemble instead of recomputing them. Output shape mirrors twasWeightsCv()
# (samplePartition + per-method <key>_predicted / <key>_performance), keyed by
# the TWAS snake method name (adapter methodKey) for a drop-in merge.
# =============================================================================

# Generate a Sample/Fold partition over the rows of X, matching the scheme in
# twasWeightsCv() (shuffle samples, then cut into `fold` contiguous blocks).
# @noRd
.fmMakeSamplePartition <- function(sampleNames, fold) {
  idx <- sample(length(sampleNames))
  folds <- cut(seq_along(sampleNames), breaks = fold, labels = FALSE)
  data.frame(Sample = sampleNames[idx], Fold = folds, stringsAsFactors = FALSE)
}

# Snake method key (e.g. "susie_inf") for a fine-mapping token, taken from the
# shared adapter registry so fineMapping CV keys match the TwasWeights `method`
# column and twasWeightsCv()'s prediction keys.
# @noRd
.fmTwasMethodKey <- function(token) {
  adapter <- .twasFineMappingMethodAdapters[[token]]
  if (is.null(adapter)) return(token)
  sub("_weights$", "", adapter$methodKey)
}

# Compact CV metric row (corr, rsq, adj_rsq, pval, RMSE, MAE) for one outcome,
# mirroring the metric block of twasWeightsCv().
# @noRd
.fmCvMetricRow <- function(pred, actual) {
  out <- setNames(rep(NA_real_, 6L),
                  c("corr", "rsq", "adj_rsq", "pval", "RMSE", "MAE"))
  ok <- !is.na(pred) & !is.na(actual)
  pred <- pred[ok]; actual <- actual[ok]
  if (length(pred) < 3L || stats::sd(pred) == 0) return(out)
  lmFit <- stats::lm(actual ~ pred); s <- summary(lmFit)
  out["corr"]    <- stats::cor(actual, pred)
  out["rsq"]     <- s$r.squared
  out["adj_rsq"] <- s$adj.r.squared
  out["pval"]    <- if (nrow(s$coefficients) >= 2L) s$coefficients[2L, 4L] else NA_real_
  res <- actual - pred
  out["RMSE"] <- sqrt(mean(res^2))
  out["MAE"]  <- mean(abs(res))
  out
}

# Fit one fine-mapping method on (Xtr, Ytr) for a CV fold and return a
# variants x outcomes weight matrix (rownames = colnames(Xtr)). susie-family
# tokens are fit independently (no chained init) per fold, matching
# twasWeightsCv's per-fold refit. Returns NULL on failure (caller skips it).
# @noRd
.fmFoldWeights <- function(token, Xtr, Ytr, coverage, userArgs, pos) {
  asMat <- function(w) {
    if (is.matrix(w)) return(w)
    matrix(w, ncol = 1L, dimnames = list(names(w), NULL))
  }
  if (token %in% c("susie", "susieInf", "susieAsh")) {
    y <- if (is.matrix(Ytr)) Ytr[, 1L] else Ytr
    fit <- .fmFitSusieIndiv(Xtr, y, token, coverage = coverage,
                            userArgs = userArgs)
    w <- switch(token,
      susie    = susieWeights(susieFit = fit),
      susieInf = susieInfWeights(susieInfFit = fit),
      susieAsh = susieAshWeights(susieAshFit = fit))
    w <- as.numeric(w)
    names(w) <- colnames(Xtr)
    return(asMat(w))
  }
  if (token == "mvsusie") {
    pv <- mvsusieR::create_mixture_prior(R = ncol(Ytr))
    fit <- do.call(fitMvsusie,
                   .fmMergeUserArgs(list(X = Xtr, Y = Ytr, prior_variance = pv,
                                         coverage = coverage),
                                    "mvsusie", userArgs))
    W <- as.matrix(mvsusieWeights(mvsusieFit = fit))
    if (is.null(rownames(W))) rownames(W) <- colnames(Xtr)
    return(W)
  }
  if (token == "fsusie") {
    fit <- do.call(fitFsusie,
                   .fmMergeUserArgs(list(X = Xtr, Y = Ytr, pos = pos),
                                    "fsusie", userArgs))
    W <- fsusieWeights(fsusieFit = fit, variantIds = colnames(Xtr))
    return(as.matrix(W))
  }
  NULL
}

# Cross-validate a homogeneous set of fine-mapping `tokens` over (X, Y). For
# univariate tokens Y is a single column; for mvsusie/fsusie Y carries one
# column per condition/feature (and fsusie additionally needs `pos`). Returns
# a list(samplePartition, prediction, performance) shaped like twasWeightsCv().
# @noRd
.fmCrossValidate <- function(X, Y, tokens, methodArgs, fold,
                             samplePartition = NULL, coverage = 0.95,
                             pos = NULL, verbose = 1) {
  if (length(tokens) == 0L) return(NULL)
  if (!is.matrix(Y)) {
    Y <- matrix(Y, ncol = 1L,
                dimnames = list(rownames(X), NULL))
  }
  if (is.null(rownames(Y))) rownames(Y) <- rownames(X)
  sampleNames <- rownames(X)
  if (is.null(samplePartition)) {
    samplePartition <- .fmMakeSamplePartition(sampleNames, fold)
  }
  foldIds <- sort(unique(samplePartition$Fold))

  preds <- setNames(
    lapply(tokens, function(tk) {
      matrix(NA_real_, nrow(Y), ncol(Y), dimnames = dimnames(Y))
    }), tokens)

  for (j in foldIds) {
    testIds <- samplePartition$Sample[samplePartition$Fold == j]
    isTest  <- rownames(X) %in% testIds
    if (all(isTest) || !any(isTest)) next
    Xtr <- X[!isTest, , drop = FALSE]
    Xte <- X[isTest, , drop = FALSE]
    Ytr <- Y[!isTest, , drop = FALSE]
    # Drop columns with zero variance in this training fold.
    keepCol <- .nonzeroVarColumns(Xtr)
    XtrK <- Xtr[, keepCol, drop = FALSE]
    for (tk in tokens) {
      W <- tryCatch(
        .fmFoldWeights(tk, XtrK, Ytr, coverage, methodArgs[[tk]], pos),
        error = function(e) {
          if (verbose >= 1)
            message(sprintf("  CV fold %s, method %s failed: %s",
                            j, tk, conditionMessage(e)))
          NULL
        })
      if (is.null(W)) next
      common <- intersect(colnames(Xte), rownames(W))
      if (length(common) == 0L) next
      yhat <- Xte[, common, drop = FALSE] %*% W[common, , drop = FALSE]
      preds[[tk]][rownames(Xte), ] <- yhat
    }
  }

  prediction  <- list()
  performance <- list()
  for (tk in tokens) {
    key <- .fmTwasMethodKey(tk)
    prediction[[paste0(key, "_predicted")]] <- preds[[tk]]
    perf <- t(vapply(seq_len(ncol(Y)), function(r) {
      .fmCvMetricRow(preds[[tk]][, r], Y[, r])
    }, numeric(6L)))
    rownames(perf) <- colnames(Y)
    performance[[paste0(key, "_performance")]] <- perf
  }
  list(samplePartition = samplePartition,
       prediction = prediction, performance = performance)
}

# Slice a full .fmCrossValidate() result down to one method's payload, keeping
# the shared samplePartition. Stored on that method's FineMappingEntry.
# @noRd
.fmSliceCv <- function(cv, token) {
  if (is.null(cv)) return(NULL)
  key <- .fmTwasMethodKey(token)
  pk <- paste0(key, "_predicted")
  mk <- paste0(key, "_performance")
  if (!pk %in% names(cv$prediction)) return(NULL)
  list(samplePartition = cv$samplePartition,
       prediction  = cv$prediction[pk],
       performance = cv$performance[mk])
}

# Rebuild a FineMappingEntry with a cvResult attached (the class is immutable).
# @noRd
.fmAttachCv <- function(entry, cvResult) {
  if (is.null(entry) || is.null(cvResult)) return(entry)
  FineMappingEntry(variantIds = entry@variantIds,
                   susieFit   = entry@susieFit,
                   topLoci    = entry@topLoci,
                   cvResult   = cvResult)
}


# =============================================================================
# QtlDataset method
# =============================================================================

#' @rdname fineMappingPipeline
#' @export
setMethod("fineMappingPipeline", "QtlDataset",
  function(data,
           methods,
           contexts           = NULL,
           traitId            = NULL,
           region             = NULL,
           cisWindow          = NULL,
           jointRegions       = FALSE,
           jointSpecification = NULL,
           addSusieInf        = TRUE,
           coverage           = 0.95,
           secondaryCoverage  = c(0.7, 0.5),
           signalCutoff       = 0.025,
           minAbsCorr         = 0.8,
           medianAbsCorr      = NULL,
           fineMappingResult  = NULL,
           cvFolds            = 0,
           samplePartition    = NULL,
           pipCutoffToSkip    = 0,
           seed               = NULL,
           naAction           = c("drop", "impute"),
           verbose            = 1,
           trim               = TRUE,
           phenotypeCovariatesToResidualize = NULL,
           genotypeCovariatesToResidualize  = NULL,
           residualizePhenotypeCovariates   = TRUE,
           residualizeGenotypeCovariates    = TRUE,
           ...) {
    naAction <- match.arg(naAction)
    if (!is.null(seed)) set.seed(as.integer(seed))
    # `cisWindow` expands a trait's own coordinates; `region` is taken
    # literally. Supplying both signals a misunderstanding -> reject.
    if (!is.null(region) && !is.null(cisWindow)) {
      stop("fineMappingPipeline(QtlDataset): specify either `region` or ",
           "`cisWindow`, not both. `cisWindow` expands each trait's own ",
           "coordinates, whereas `region` is the literal variant window.")
    }
    xRegions <- .makeXRegions(region, jointRegions)
    parsedJointSpec <- parseJointSpecification(jointSpecification, data)
    norm       <- .fmNormalizeMethods(methods)
    tokens     <- norm$tokens
    methodArgs <- norm$methodArgs
    .fmCheckMethodCapabilities(tokens, "QtlDataset")

    # Explicit jointSpecification path: run the per-spec axis dispatcher for
    # the multi-axis methods (mvsusie / fsusie) and remove them from the
    # per-tuple loop's token set so they aren't fitted twice. Non-joint
    # methods continue through the existing per-(context, trait) iteration
    # below.
    jointResult <- NULL
    if (length(parsedJointSpec) > 0L) {
      jointResult <- .fmDispatchJointSpecsQtlDataset(
        parsedJointSpec, data, intersect(tokens, c("mvsusie", "fsusie")),
        contexts, traitId, cisWindow,
        coverage, secondaryCoverage, signalCutoff, minAbsCorr, verbose,
        methodArgs = methodArgs, xRegions = xRegions)
      tokens <- setdiff(tokens, c("mvsusie", "fsusie"))
      methodArgs <- methodArgs[tokens]
      if (length(tokens) == 0L) {
        if (is.null(jointResult))
          stop("fineMappingPipeline(QtlDataset): no joint fits produced. ",
               "Check that the jointSpecification scope intersects the ",
               "available studies / contexts / traits.")
        return(jointResult)
      }
    }

    study <- getStudy(data)
    allCtx <- getContexts(data)
    useCtx <- if (is.null(contexts)) allCtx else {
      bad <- setdiff(contexts, allCtx)
      if (length(bad) > 0L) {
        stop("fineMappingPipeline(QtlDataset): unknown context(s): ",
             paste(bad, collapse = ", "))
      }
      contexts
    }

    # Per-context trait list (intersect requested traitId or region with
    # each context's available rows). Mirrors twasWeightsPipeline.
    perCtxTraits <- vector("list", length(useCtx))
    names(perCtxTraits) <- useCtx
    for (ctx in useCtx) {
      se <- getPhenotypes(data, contexts = ctx)
      ids <- rownames(se)
      if (!is.null(traitId)) {
        ids <- intersect(ids, traitId)
      } else if (!is.null(region)) {
        rr <- SummarizedExperiment::rowRanges(se)
        ids <- ids[IRanges::overlapsAny(rr, region)]
      }
      perCtxTraits[[ctx]] <- ids
    }
    allTraits <- unique(unlist(perCtxTraits, use.names = FALSE))
    if (length(allTraits) == 0L) {
      stop("fineMappingPipeline(QtlDataset): no traits selected.")
    }

    nCtx <- length(useCtx)
    nTraits <- length(allTraits)

    # Partition tokens by univariate / multivariate / fsusie.
    isUniv  <- tokens %in% c("susie", "susieInf", "susieAsh")
    univTokens   <- tokens[isUniv]
    isMv    <- tokens == "mvsusie"
    mvTokens     <- tokens[isMv]
    isFs    <- tokens == "fsusie"
    fsTokens     <- tokens[isFs]

    # Multivariate guard: mvsusie requires multi-trait OR multi-context;
    # fsusie always requires multi-trait per context.
    if (length(mvTokens) > 0L && nCtx < 2L && nTraits < 2L) {
      stop("fineMappingPipeline(QtlDataset): mvsusie requires multi-trait or ",
           "multi-context input (got ", nTraits, " trait(s) x ", nCtx,
           " context(s)).")
    }
    if (length(fsTokens) > 0L && nTraits < 2L) {
      stop("fineMappingPipeline(QtlDataset): fsusie requires multi-trait ",
           "input within a context (got ", nTraits, " trait(s)).")
    }

    chain <- .fmResolveSusieChain(univTokens, addSusieInf)

    rowStudy   <- character(0)
    rowContext <- character(0)
    rowTrait   <- character(0)
    rowMethod  <- character(0)
    rowEntries <- list()

    pushRow <- function(st, ctx, tr, mt, ent) {
      rowStudy   <<- c(rowStudy,   st)
      rowContext <<- c(rowContext, ctx)
      rowTrait   <<- c(rowTrait,   tr)
      rowMethod  <<- c(rowMethod,  mt)
      rowEntries[[length(rowEntries) + 1L]] <<- ent
    }

    # ---- Univariate dispatch: per (context, trait), per method.
    # X is drawn from each window in `xRegions` (cis = one trait-derived block;
    # multi-region = one block per range), fitted independently, then merged
    # per token via .fmMergeEntries so every (study, context, trait, method)
    # produces a single entry.
    if (length(univTokens) > 0L) {
      for (ctx in useCtx) {
        for (tid in perCtxTraits[[ctx]]) {
          # Resume-cache lookup per (ctx, tid, token).
          toRun <- character(0)
          for (tk in univTokens) {
            cached <- .fmCacheLookup(fineMappingResult, study, ctx, tid, tk)
            if (!is.null(cached)) {
              pushRow(study, ctx, tid, tk, cached)
            } else {
              toRun <- c(toRun, tk)
            }
          }
          if (length(toRun) == 0L) next

          Y <- .fmResidPheno(
            data, contexts = ctx, traitId = tid, naAction = naAction)

          blockEntries <- lapply(xRegions, function(rg) {
            X <- if (is.null(rg)) {
              .fmResidGeno(data, contexts = ctx, traitId = tid,
                           cisWindow = cisWindow, samples = rownames(Y))
            } else {
              .fmResidGeno(data, contexts = ctx, region = rg,
                           samples = rownames(Y))
            }
            common <- intersect(rownames(X), rownames(Y))
            if (length(common) < 2L) {
              stop(sprintf(
                "fineMappingPipeline: too few shared samples between residualized X and Y for (context='%s', trait='%s').",
                ctx, tid))
            }
            X <- X[common, , drop = FALSE]
            y <- Y[common, , drop = FALSE]
            if (ncol(y) > 1L) y <- y[, 1L, drop = TRUE] else y <- drop(y)
            # SER pre-screen: skip this block when a single-effect fit finds no
            # PIP above pipCutoffToSkip (no potentially significant variant).
            if (!.fmSerScreen(X, y, pipCutoffToSkip)) {
              if (verbose >= 1)
                message(sprintf(
                  "Skipping (context='%s', trait='%s'): SER pre-screen found no PIP above pipCutoffToSkip.",
                  ctx, tid))
              return(list())
            }
            afVec <- .fmAfForX(data, X, traitId = tid, region = rg,
                               cisWindow = cisWindow)
            .fmFitXBlock(X, y, toRun, addSusieInf, coverage,
                         secondaryCoverage, signalCutoff, minAbsCorr,
                         methodArgs, verbose, ctx, tid,
                         cvFolds = cvFolds, samplePartition = samplePartition,
                         af = afVec)
          })

          for (tk in toRun) {
            ents <- lapply(blockEntries, function(be) be[[tk]])
            if (any(vapply(ents, is.null, logical(1)))) next
            entry <- if (length(ents) == 1L) ents[[1L]] else .fmMergeEntries(ents)
            pushRow(study, ctx, tid, tk, entry)
          }
        }
      }
    }

    # ---- mvsusie dispatch: joint over selected (contexts, traits).
    if (length(mvTokens) > 0L) {
      if (!requireNamespace("mvsusieR", quietly = TRUE)) {
        stop("mvsusie requires the mvsusieR package. Install with: ",
             "devtools::install_github('stephenslab/mvsusieR')")
      }
      # Detection: when multiple contexts AND single trait => multi-context
      # mvsusie (group by trait). When single context AND multiple traits =>
      # multi-trait mvsusie (group by context). When both multi, iterate per
      # context for the multi-trait fit (same convention as the design doc:
      # "sequential per-context multi-trait when both are multi").
      mvJobs <- list()
      if (nCtx >= 2L && nTraits == 1L) {
        # Single trait across many contexts.
        mvJobs[[length(mvJobs) + 1L]] <- list(
          mode = "multiContext", trait = allTraits[[1L]],
          contexts = useCtx)
      } else if (nCtx == 1L && nTraits >= 2L) {
        mvJobs[[length(mvJobs) + 1L]] <- list(
          mode = "multiTrait", context = useCtx[[1L]],
          traits = perCtxTraits[[useCtx[[1L]]]])
      } else {
        # Both multi => sequential per-context multi-trait fit.
        for (ctx in useCtx) {
          tr <- perCtxTraits[[ctx]]
          if (length(tr) < 2L) next
          mvJobs[[length(mvJobs) + 1L]] <- list(
            mode = "multiTrait", context = ctx, traits = tr)
        }
      }

      for (job in mvJobs) {
        if (identical(job$mode, "multiContext")) {
          tid <- job$trait
          # Joint Y across contexts for this single trait. X is drawn from each
          # region block (cis or explicit region) and merged across regions.
          contextsHere <- job$contexts
          Yres <- .fmResidPheno(
            data, contexts = contextsHere, traitId = tid, naAction = naAction)
          if (length(contextsHere) == 1L)
            Yres <- setNames(list(Yres), contextsHere)
          baseSamples <- Reduce(intersect, lapply(Yres, rownames))

          # Resume cache: every (study, ctx, tid, mvsusie) row.
          allCached <- TRUE
          for (ctx in contextsHere) {
            if (is.null(.fmCacheLookup(fineMappingResult, study, ctx, tid, "mvsusie"))) {
              allCached <- FALSE; break
            }
          }
          if (allCached) {
            for (ctx in contextsHere) {
              pushRow(study, ctx, tid, "mvsusie",
                .fmCacheLookup(fineMappingResult, study, ctx, tid, "mvsusie"))
            }
            next
          }

          # SER pre-screen: drop contexts with no single-effect signal before
          # the joint fit (faithful port of skipConditions). Screen the first
          # region block; skip the trait entirely when < 2 contexts survive.
          if (.fmScreenActive(pipCutoffToSkip)) {
            rg0  <- xRegions[[1L]]
            Xscr <- if (is.null(rg0)) {
              .fmResidGeno(data, contexts = contextsHere, traitId = tid,
                           cisWindow = cisWindow, samples = baseSamples)
            } else {
              .fmResidGeno(data, contexts = contextsHere, region = rg0,
                           samples = baseSamples)
            }
            csS <- intersect(baseSamples, rownames(Xscr))
            if (length(csS) >= 2L) {
              Yscr <- do.call(cbind, lapply(contextsHere,
                function(ctx) Yres[[ctx]][csS, 1L]))
              kept <- contextsHere[.fmSerScreenColumns(
                Xscr[csS, , drop = FALSE], Yscr, pipCutoffToSkip)]
              if (length(kept) < 2L) {
                if (verbose >= 1)
                  message(sprintf(
                    "Skipping mvsusie (multi-context) for trait='%s': < 2 contexts pass the SER pre-screen.",
                    tid))
                next
              }
              if (length(kept) < length(contextsHere)) {
                if (verbose >= 1)
                  message(sprintf(
                    "mvsusie (multi-context) trait='%s': SER pre-screen kept %d of %d contexts.",
                    tid, length(kept), length(contextsHere)))
                contextsHere <- kept
              }
            }
          }

          if (verbose >= 1)
            message(sprintf("Fitting mvsusie (multi-context) for trait='%s' ...", tid))
          fitOneRegion <- function(rg) {
            X <- if (is.null(rg)) {
              .fmResidGeno(data, contexts = contextsHere, traitId = tid,
                           cisWindow = cisWindow, samples = baseSamples)
            } else {
              .fmResidGeno(data, contexts = contextsHere, region = rg,
                           samples = baseSamples)
            }
            cs <- intersect(baseSamples, rownames(X))
            if (length(cs) < 2L) {
              stop("fineMappingPipeline(QtlDataset, mvsusie multi-context): ",
                   "insufficient shared samples across selected contexts.")
            }
            Xc <- X[cs, , drop = FALSE]
            afVec <- .fmAfForX(data, Xc, traitId = tid, region = rg,
                               cisWindow = cisWindow)
            Yc <- do.call(cbind, lapply(contextsHere, function(ctx) {
              ym <- Yres[[ctx]][cs, , drop = FALSE]
              colnames(ym) <- ctx
              ym
            }))
            mvBaseArgs <- list(
              X = Xc, Y = Yc,
              prior_variance = mvsusieR::create_mixture_prior(R = ncol(Yc)),
              coverage = coverage)
            fit <- do.call(fitMvsusie,
                           .fmMergeUserArgs(mvBaseArgs, "mvsusie",
                                            methodArgs[["mvsusie"]]))
            fit <- .setFinemappingFitClass(fit, "mvsusie")
            entry <- .fmPostprocessOne(
              fit = fit, method = "mvsusie", dataX = Xc, dataY = NULL,
              coverage = coverage, secondaryCoverage = secondaryCoverage,
              signalCutoff = signalCutoff, minAbsCorr = minAbsCorr,
              af = afVec, csInput = "X")
            if (cvFolds > 1L) {
              cv <- .fmCrossValidate(Xc, Yc, "mvsusie", methodArgs, cvFolds,
                                     samplePartition = samplePartition,
                                     coverage = coverage, verbose = verbose)
              entry <- .fmAttachCv(entry, .fmSliceCv(cv, "mvsusie"))
            }
            entry
          }
          entry <- .fmJointBlocks(xRegions, fitOneRegion)
          # Share the joint (merged) entry across contexts via copy-on-modify.
          for (ctx in contextsHere) {
            pushRow(study, ctx, tid, "mvsusie", entry)
          }

        } else {  # multiTrait
          ctx <- job$context
          traits <- job$traits
          # Resume cache: every (study, ctx, trait, mvsusie) row.
          allCached <- TRUE
          for (tid in traits) {
            if (is.null(.fmCacheLookup(fineMappingResult, study, ctx, tid, "mvsusie"))) {
              allCached <- FALSE; break
            }
          }
          if (allCached) {
            for (tid in traits) {
              pushRow(study, ctx, tid, "mvsusie",
                .fmCacheLookup(fineMappingResult, study, ctx, tid, "mvsusie"))
            }
            next
          }

          Y <- .fmResidPheno(
            data, contexts = ctx, traitId = traits, naAction = naAction)

          # SER pre-screen: drop traits with no single-effect signal before the
          # joint fit (faithful port of skipConditions). Skip the context's
          # mvsusie when < 2 traits survive.
          if (.fmScreenActive(pipCutoffToSkip)) {
            rg0  <- xRegions[[1L]]
            Xscr <- if (is.null(rg0)) {
              .fmResidGeno(data, contexts = ctx, traitId = traits,
                           cisWindow = cisWindow, samples = rownames(Y))
            } else {
              .fmResidGeno(data, contexts = ctx, region = rg0,
                           samples = rownames(Y))
            }
            csS <- intersect(rownames(Xscr), rownames(Y))
            if (length(csS) >= 2L) {
              keep <- .fmSerScreenColumns(
                Xscr[csS, , drop = FALSE], Y[csS, , drop = FALSE],
                pipCutoffToSkip)
              if (sum(keep) < 2L) {
                if (verbose >= 1)
                  message(sprintf(
                    "Skipping mvsusie (multi-trait) for context='%s': < 2 traits pass the SER pre-screen.",
                    ctx))
                next
              }
              if (sum(keep) < length(traits)) {
                if (verbose >= 1)
                  message(sprintf(
                    "mvsusie (multi-trait) context='%s': SER pre-screen kept %d of %d traits.",
                    ctx, sum(keep), length(traits)))
                traits <- traits[keep]
                Y <- Y[, keep, drop = FALSE]
              }
            }
          }

          if (verbose >= 1)
            message(sprintf("Fitting mvsusie (multi-trait) for context='%s' ...", ctx))
          fitOneRegion <- function(rg) {
            X <- if (is.null(rg)) {
              .fmResidGeno(data, contexts = ctx, traitId = traits,
                           cisWindow = cisWindow, samples = rownames(Y))
            } else {
              .fmResidGeno(data, contexts = ctx, region = rg,
                           samples = rownames(Y))
            }
            common <- intersect(rownames(X), rownames(Y))
            if (length(common) < 2L) {
              stop(sprintf(
                "fineMappingPipeline(QtlDataset, mvsusie multi-trait): too few shared samples in context '%s'.",
                ctx))
            }
            Xc <- X[common, , drop = FALSE]
            Yc <- Y[common, , drop = FALSE]
            afVec <- .fmAfForX(data, Xc, traitId = traits, region = rg,
                               cisWindow = cisWindow)
            mvBaseArgs <- list(
              X = Xc, Y = Yc,
              prior_variance = mvsusieR::create_mixture_prior(R = ncol(Yc)),
              coverage = coverage)
            fit <- do.call(fitMvsusie,
                           .fmMergeUserArgs(mvBaseArgs, "mvsusie",
                                            methodArgs[["mvsusie"]]))
            fit <- .setFinemappingFitClass(fit, "mvsusie")
            entry <- .fmPostprocessOne(
              fit = fit, method = "mvsusie", dataX = Xc, dataY = NULL,
              coverage = coverage, secondaryCoverage = secondaryCoverage,
              signalCutoff = signalCutoff, minAbsCorr = minAbsCorr,
              af = afVec, csInput = "X")
            if (cvFolds > 1L) {
              cv <- .fmCrossValidate(Xc, Yc, "mvsusie", methodArgs, cvFolds,
                                     samplePartition = samplePartition,
                                     coverage = coverage, verbose = verbose)
              entry <- .fmAttachCv(entry, .fmSliceCv(cv, "mvsusie"))
            }
            entry
          }
          entry <- .fmJointBlocks(xRegions, fitOneRegion)
          for (tid in traits) {
            pushRow(study, ctx, tid, "mvsusie", entry)
          }
        }
      }
    }

    # ---- fsusie dispatch: joint multi-trait per context.
    if (length(fsTokens) > 0L) {
      if (!requireNamespace("fsusieR", quietly = TRUE)) {
        stop("fsusie requires the fsusieR package. Install with: ",
             "devtools::install_github('stephenslab/fsusieR')")
      }
      for (ctx in useCtx) {
        traits <- perCtxTraits[[ctx]]
        if (length(traits) < 2L) {
          stop(sprintf(
            "fineMappingPipeline(QtlDataset, fsusie): context '%s' has %d trait(s); fsusie needs at least 2 within a context.",
            ctx, length(traits)))
        }
        # Resume cache.
        allCached <- TRUE
        for (tid in traits) {
          if (is.null(.fmCacheLookup(fineMappingResult, study, ctx, tid, "fsusie"))) {
            allCached <- FALSE; break
          }
        }
        if (allCached) {
          for (tid in traits) {
            pushRow(study, ctx, tid, "fsusie",
              .fmCacheLookup(fineMappingResult, study, ctx, tid, "fsusie"))
          }
          next
        }

        Y <- .fmResidPheno(
          data, contexts = ctx, traitId = traits, naAction = naAction)

        # Per-trait genomic positions for the wavelet model. Region-independent
        # (depends on the trait set / Y columns): midpoint of each trait range.
        se <- getPhenotypes(data, contexts = ctx, traitId = traits)
        rrIds <- rownames(se)
        ord <- match(colnames(Y), rrIds)
        if (anyNA(ord)) {
          stop("fineMappingPipeline(QtlDataset, fsusie): unable to align trait positions to Y columns.")
        }
        rr <- SummarizedExperiment::rowRanges(se)[ord]
        pos <- (GenomicRanges::start(rr) + GenomicRanges::end(rr)) / 2

        if (verbose >= 1)
          message(sprintf("Fitting fsusie for context='%s' (multi-trait, %d traits) ...",
                          ctx, length(traits)))
        fitOneRegion <- function(rg) {
          X <- if (is.null(rg)) {
            .fmResidGeno(data, contexts = ctx, traitId = traits,
                         cisWindow = cisWindow, samples = rownames(Y))
          } else {
            .fmResidGeno(data, contexts = ctx, region = rg,
                         samples = rownames(Y))
          }
          common <- intersect(rownames(X), rownames(Y))
          if (length(common) < 2L) {
            stop(sprintf("fineMappingPipeline(QtlDataset, fsusie): too few shared samples in context '%s'.", ctx))
          }
          Xc <- X[common, , drop = FALSE]
          Yc <- Y[common, , drop = FALSE]
          afVec <- .fmAfForX(data, Xc, traitId = traits, region = rg,
                             cisWindow = cisWindow)
          fit <- do.call(fitFsusie,
                         .fmMergeUserArgs(list(X = Xc, Y = Yc, pos = pos),
                                          "fsusie", methodArgs[["fsusie"]]))
          # Collapse the functional fit to a variants x features TWAS weight
          # matrix now, while fitted_wc/csd_X are still present (trimming drops
          # them). Stored on $coef so a trimmed fit can still yield weights.
          fit$coef <- tryCatch(
            fsusieWeights(fsusieFit = fit, variantIds = colnames(Xc)),
            error = function(e) NULL)
          fit <- .setFinemappingFitClass(fit, "fsusie")
          entry <- .fmPostprocessOne(
            fit = fit, method = "fsusie", dataX = Xc, dataY = NULL,
            coverage = coverage, secondaryCoverage = secondaryCoverage,
            signalCutoff = signalCutoff, minAbsCorr = minAbsCorr,
            af = afVec, csInput = "fsusie")
          if (cvFolds > 1L) {
            cv <- .fmCrossValidate(Xc, Yc, "fsusie", methodArgs, cvFolds,
                                   samplePartition = samplePartition,
                                   coverage = coverage, pos = pos,
                                   verbose = verbose)
            entry <- .fmAttachCv(entry, .fmSliceCv(cv, "fsusie"))
          }
          entry
        }
        entry <- .fmJointBlocks(xRegions, fitOneRegion)
        for (tid in traits) {
          pushRow(study, ctx, tid, "fsusie", entry)
        }
      }
    }

    perTupleResult <- if (length(rowEntries) > 0L)
      .fmBuildQtlResult(rowStudy, rowContext, rowTrait, rowMethod, rowEntries,
                        ldSketch = NULL)
      else NULL

    if (is.null(jointResult)) {
      if (is.null(perTupleResult))
        stop("fineMappingPipeline: no (study, context, trait, method) tuples ",
             "produced a fine-mapping result.")
      return(perTupleResult)
    }
    if (is.null(perTupleResult)) return(jointResult)
    .rbindFineMappingResult(perTupleResult, jointResult, ldSketch = NULL)
  })


# =============================================================================
# MultiStudyQtlDataset method
# =============================================================================

#' @rdname fineMappingPipeline
#' @export
setMethod("fineMappingPipeline", "MultiStudyQtlDataset",
  function(data,
           methods,
           contexts           = NULL,
           traitId            = NULL,
           region             = NULL,
           cisWindow          = NULL,
           jointRegions       = FALSE,
           jointSpecification = NULL,
           addSusieInf        = TRUE,
           coverage           = 0.95,
           secondaryCoverage  = c(0.7, 0.5),
           signalCutoff       = 0.025,
           minAbsCorr         = 0.8,
           medianAbsCorr      = NULL,
           fineMappingResult  = NULL,
           cvFolds            = 0,
           samplePartition    = NULL,
           pipCutoffToSkip    = 0,
           seed               = NULL,
           naAction           = c("drop", "impute"),
           verbose            = 1,
           trim               = TRUE,
           phenotypeCovariatesToResidualize = NULL,
           genotypeCovariatesToResidualize  = NULL,
           residualizePhenotypeCovariates   = TRUE,
           residualizeGenotypeCovariates    = TRUE,
           ...) {
    naAction <- match.arg(naAction)
    if (!is.null(region) && !is.null(cisWindow)) {
      stop("fineMappingPipeline(MultiStudyQtlDataset): specify either ",
           "`region` or `cisWindow`, not both.")
    }
    xRegions <- .makeXRegions(region, jointRegions)
    parsedJointSpec <- parseJointSpecification(jointSpecification, data)
    norm       <- .fmNormalizeMethods(methods)
    tokens     <- norm$tokens
    methodArgs <- norm$methodArgs
    .fmCheckMethodCapabilities(tokens, "MultiStudyQtlDataset")

    # Explicit jointSpecification path: run the per-component, per-spec
    # joint dispatcher and remove joint-eligible methods from the
    # per-tuple recursion below.
    jointResult <- NULL
    if (length(parsedJointSpec) > 0L) {
      jointResult <- .fmDispatchJointSpecsMultiStudy(
        parsedJointSpec, data, intersect(tokens, c("mvsusie", "fsusie")),
        contexts, traitId, cisWindow,
        coverage, secondaryCoverage, signalCutoff, minAbsCorr, verbose,
        methodArgs = methodArgs, xRegions = xRegions)
      # Forward the still-pending (non-joint) tokens + their kwargs to the
      # per-QtlDataset recursion below, preserving the list shape so
      # methodArgs land on the right tokens.
      tokens <- setdiff(tokens, c("mvsusie", "fsusie"))
      methodArgs <- methodArgs[tokens]
      methods <- if (length(methodArgs) > 0L) methodArgs else tokens
      if (length(tokens) == 0L) {
        if (is.null(jointResult))
          stop("fineMappingPipeline(MultiStudyQtlDataset): no joint fits produced. ",
               "Check that the jointSpecification scope intersects the available data.")
        return(jointResult)
      }
    }

    qtlDatasets <- getQtlDatasets(data)
    sumStats <- getSumStats(data)

    out <- NULL
    embeddedLd <- NULL
    for (qdName in names(qtlDatasets)) {
      qd <- qtlDatasets[[qdName]]
      res <- fineMappingPipeline(
        data               = qd,
        methods            = methods,
        contexts           = contexts,
        traitId            = traitId,
        region             = region,
        cisWindow          = cisWindow,
        jointRegions       = jointRegions,
        jointSpecification = NULL,
        addSusieInf        = addSusieInf,
        coverage           = coverage,
        secondaryCoverage  = secondaryCoverage,
        signalCutoff       = signalCutoff,
        minAbsCorr         = minAbsCorr,
        fineMappingResult  = fineMappingResult,
        cvFolds            = cvFolds,
        samplePartition    = samplePartition,
        pipCutoffToSkip    = pipCutoffToSkip,
        seed               = seed,
        naAction           = naAction,
        verbose            = verbose,
        ...)
      out <- if (is.null(out)) res else .rbindFineMappingResult(out, res, ldSketch = NULL)
    }

    if (!is.null(sumStats)) {
      ssRes <- fineMappingPipeline(
        data               = sumStats,
        methods            = methods,
        contexts           = contexts,
        traitId            = traitId,
        jointSpecification = NULL,
        addSusieInf        = addSusieInf,
        coverage           = coverage,
        secondaryCoverage  = secondaryCoverage,
        signalCutoff       = signalCutoff,
        minAbsCorr         = minAbsCorr,
        fineMappingResult  = fineMappingResult,
        verbose            = verbose,
        ...)
      embeddedLd <- getLdSketch(ssRes)
      out <- if (is.null(out)) ssRes else .rbindFineMappingResult(out, ssRes,
        ldSketch = embeddedLd)
    }

    perTupleResult <- if (!is.null(out)) {
      # ldSketch: NULL if all studies were individual-level; the embedded
      # sumStats's ldSketch otherwise.
      QtlFineMappingResult(
        study         = as.character(out$study),
        context       = as.character(out$context),
        trait         = as.character(out$trait),
        method        = as.character(out$method),
        entry         = as.list(out$entry),
        jointStudies  = if ("jointStudies"  %in% names(out))
                          as.character(out$jointStudies)  else NULL,
        jointContexts = if ("jointContexts" %in% names(out))
                          as.character(out$jointContexts) else NULL,
        jointTraits   = if ("jointTraits"   %in% names(out))
                          as.character(out$jointTraits)   else NULL,
        ldSketch      = embeddedLd)
    } else NULL

    if (is.null(jointResult)) {
      if (is.null(perTupleResult))
        stop("fineMappingPipeline(MultiStudyQtlDataset): no entries produced a result.")
      return(perTupleResult)
    }
    if (is.null(perTupleResult)) return(jointResult)
    .rbindFineMappingResult(perTupleResult, jointResult, ldSketch = embeddedLd)
  })


# =============================================================================
# QtlSumStats method
# =============================================================================

#' @rdname fineMappingPipeline
#' @export
setMethod("fineMappingPipeline", "QtlSumStats",
  function(data,
           methods,
           contexts           = NULL,
           traitId            = NULL,
           jointSpecification = NULL,
           addSusieInf        = TRUE,
           coverage           = 0.95,
           secondaryCoverage  = c(0.7, 0.5),
           signalCutoff       = 0.025,
           minAbsCorr         = 0.8,
           medianAbsCorr      = NULL,
           fineMappingResult  = NULL,
           verbose            = 1,
           trim               = TRUE,
           ...) {
    .fmAssertQcd(data)
    parsedJointSpec <- parseJointSpecification(jointSpecification, data)
    norm       <- .fmNormalizeMethods(methods)
    tokens     <- norm$tokens
    methodArgs <- norm$methodArgs
    .fmCheckMethodCapabilities(tokens, "QtlSumStats")

    jointResult <- NULL
    if (length(parsedJointSpec) > 0L) {
      jointResult <- .fmDispatchJointSpecsQtlSumStats(
        parsedJointSpec, data, intersect(tokens, "mvsusie"),
        contexts, traitId,
        coverage, secondaryCoverage, signalCutoff, minAbsCorr, verbose,
        methodArgs = methodArgs)
      tokens <- setdiff(tokens, c("mvsusie", "fsusie"))
      methodArgs <- methodArgs[tokens]
      if (length(tokens) == 0L) {
        if (is.null(jointResult))
          stop("fineMappingPipeline(QtlSumStats): no joint fits produced. ",
               "Check that the jointSpecification scope intersects the available data.")
        return(jointResult)
      }
    }

    studyCol   <- as.character(data$study)
    contextCol <- as.character(data$context)
    traitCol   <- as.character(data$trait)

    selRows <- seq_len(nrow(data))
    if (!is.null(contexts)) selRows <- selRows[contextCol[selRows] %in% contexts]
    if (!is.null(traitId))  selRows <- selRows[traitCol[selRows]   %in% traitId]
    if (length(selRows) == 0L) {
      stop("fineMappingPipeline(QtlSumStats): no entries matched the supplied ",
           "contexts / traitId filters.")
    }

    isUniv <- tokens %in% c("susie", "susieInf", "susieAsh")
    univTokens <- tokens[isUniv]
    mvTokens   <- tokens[tokens == "mvsusie"]

    if (length(mvTokens) > 0L) {
      groupKey <- paste(studyCol[selRows], traitCol[selRows], sep = "||")
      perGroupNCtx <- vapply(split(contextCol[selRows], groupKey),
                             length, integer(1))
      if (all(perGroupNCtx < 2L)) {
        stop("fineMappingPipeline(QtlSumStats): mvsusie requires at least two ",
             "contexts per (study, trait); the supplied collection has only ",
             "one context per trait.")
      }
    }

    ldSketch <- getLdSketch(data)

    rowStudy   <- character(0)
    rowContext <- character(0)
    rowTrait   <- character(0)
    rowMethod  <- character(0)
    rowEntries <- list()
    pushRow <- function(st, ctx, tr, mt, ent) {
      rowStudy   <<- c(rowStudy,   st)
      rowContext <<- c(rowContext, ctx)
      rowTrait   <<- c(rowTrait,   tr)
      rowMethod  <<- c(rowMethod,  mt)
      rowEntries[[length(rowEntries) + 1L]] <<- ent
    }

    # ---- Univariate dispatch: per (study, context, trait), per method.
    if (length(univTokens) > 0L) {
      for (i in selRows) {
        st <- studyCol[i]; ctx <- contextCol[i]; tr <- traitCol[i]

        # Cache hits first, then determine which tokens still need fitting.
        toRun <- character(0)
        for (tk in univTokens) {
          cached <- .fmCacheLookup(fineMappingResult, st, ctx, tr, tk)
          if (!is.null(cached)) {
            pushRow(st, ctx, tr, tk, cached)
          } else {
            toRun <- c(toRun, tk)
          }
        }
        if (length(toRun) == 0L) next

        entry <- data$entry[[i]]
        zn <- .fmExtractZN(entry,
          sprintf("fineMappingPipeline(QtlSumStats): entry %d (study='%s', context='%s', trait='%s')", i, st, ctx, tr))
        variantIds <- zn$variantIds
        z <- zn$z
        n <- zn$n
        # Effect-allele frequency for export as `af` (entry MAF mcol, post-QC
        # harmonized/complemented); aligned to variantIds, NULL -> af NA.
        .qmc <- S4Vectors::mcols(entry)
        afByVar <- if ("MAF" %in% colnames(.qmc))
          setNames(as.numeric(.qmc$MAF), as.character(.qmc$SNP))[variantIds] else NULL
        ldMat <- .fmLdFromSketch(ldSketch, variantIds)
        names(z) <- variantIds

        chainLocal <- .fmResolveSusieChain(toRun, addSusieInf)
        infFit <- NULL
        if (chainLocal$runInf) {
          if (verbose >= 1)
            message(sprintf("Fitting susieInf (RSS) for (study='%s', context='%s', trait='%s') ...", st, ctx, tr))
          infFit <- .fmFitSusieRss(z, ldMat, n, "susieInf",
                                   coverage = coverage,
                                   userArgs = methodArgs[["susieInf"]])
        }
        for (tk in toRun) {
          if (tk == "susieInf") {
            if (!chainLocal$keepInf) next
            fit <- infFit
          } else {
            chainFrom <- if ((tk == "susie"    && chainLocal$chainSusie) ||
                             (tk == "susieAsh" && chainLocal$chainAsh))
                          infFit else NULL
            if (verbose >= 1)
              message(sprintf("Fitting %s (RSS) for (study='%s', context='%s', trait='%s') ...",
                              tk, st, ctx, tr))
            fit <- .fmFitSusieRss(z, ldMat, n, tk,
                                  chainFromInf = chainFrom,
                                  coverage = coverage,
                                  userArgs = methodArgs[[tk]])
          }
          ent <- .fmPostprocessOne(
            fit = fit, method = "susieRss",
            dataX = ldMat, dataY = list(z = z),
            coverage = coverage,
            secondaryCoverage = secondaryCoverage,
            signalCutoff = signalCutoff,
            minAbsCorr = minAbsCorr,
            af = afByVar,
            csInput = "Xcorr")
          # The method column on the FineMappingResult carries the bare
          # token (susie / susieInf / susieAsh), independent of which
          # postprocess class was used.
          pushRow(st, ctx, tr, tk, ent)
        }
      }
    }

    # ---- mvsusie dispatch: per (study, trait) across selected contexts.
    if (length(mvTokens) > 0L) {
      if (!requireNamespace("mvsusieR", quietly = TRUE)) {
        stop("mvsusie requires the mvsusieR package. Install with: ",
             "devtools::install_github('stephenslab/mvsusieR')")
      }
      groupKey <- paste(studyCol[selRows], traitCol[selRows], sep = "||")
      groups <- split(selRows, groupKey)
      for (gkey in names(groups)) {
        gIdx <- groups[[gkey]]
        if (length(gIdx) < 2L) next
        st <- studyCol[gIdx[[1L]]]
        tr <- traitCol[gIdx[[1L]]]
        ctxNames <- contextCol[gIdx]

        # Resume cache.
        allCached <- TRUE
        for (ctx in ctxNames) {
          if (is.null(.fmCacheLookup(fineMappingResult, st, ctx, tr, "mvsusie"))) {
            allCached <- FALSE; break
          }
        }
        if (allCached) {
          for (ctx in ctxNames) {
            pushRow(st, ctx, tr, "mvsusie",
              .fmCacheLookup(fineMappingResult, st, ctx, tr, "mvsusie"))
          }
          next
        }

        firstMc <- S4Vectors::mcols(data$entry[[gIdx[[1L]]]])
        if (!"SNP" %in% colnames(firstMc))
          stop("fineMappingPipeline(QtlSumStats, mvsusie): entry has no SNP mcol.")
        variantIds <- as.character(firstMc$SNP)
        Z <- matrix(NA_real_, nrow = length(variantIds), ncol = length(gIdx),
                    dimnames = list(variantIds, ctxNames))
        nVec <- numeric(length(gIdx))
        for (kk in seq_along(gIdx)) {
          mc <- S4Vectors::mcols(data$entry[[gIdx[kk]]])
          if (!identical(as.character(mc$SNP), variantIds)) {
            stop("fineMappingPipeline(QtlSumStats, mvsusie): every entry in ",
                 "the (study='", st, "', trait='", tr, "') group must share an ",
                 "identical SNP order after summaryStatsQc().")
          }
          Z[, kk] <- as.numeric(mc$Z)
          nVec[kk] <- stats::median(as.numeric(mc$N), na.rm = TRUE)
        }
        ldMat <- .fmLdFromSketch(ldSketch, variantIds)

        if (verbose >= 1)
          message(sprintf("Fitting mvsusie (RSS) for (study='%s', trait='%s', %d contexts) ...",
                          st, tr, length(ctxNames)))
        mvBaseArgs <- list(
          Z = Z, R = ldMat, N = as.numeric(stats::median(nVec)),
          prior_variance = mvsusieR::create_mixture_prior(R = ncol(Z)),
          coverage = coverage)
        fit <- do.call(fitMvsusieRss,
                       .fmMergeUserArgs(mvBaseArgs, "mvsusie",
                                        methodArgs[["mvsusie"]]))
        fit <- .setFinemappingFitClass(fit, "mvsusie")
        ent <- .fmPostprocessOne(
          fit = fit, method = "mvsusie",
          dataX = ldMat, dataY = NULL,
          coverage = coverage,
          secondaryCoverage = secondaryCoverage,
          signalCutoff = signalCutoff,
          minAbsCorr = minAbsCorr,
          csInput = "Xcorr")
        for (ctx in ctxNames) {
          pushRow(st, ctx, tr, "mvsusie", ent)
        }
      }
    }

    perTupleResult <- if (length(rowEntries) > 0L)
      .fmBuildQtlResult(rowStudy, rowContext, rowTrait, rowMethod, rowEntries,
                        ldSketch = ldSketch)
      else NULL
    if (is.null(jointResult)) {
      if (is.null(perTupleResult))
        stop("fineMappingPipeline(QtlSumStats): no entries produced a result.")
      return(perTupleResult)
    }
    if (is.null(perTupleResult)) return(jointResult)
    .rbindFineMappingResult(perTupleResult, jointResult, ldSketch = ldSketch)
  })


# =============================================================================
# GwasSumStats method
# =============================================================================

#' @rdname fineMappingPipeline
#' @export
setMethod("fineMappingPipeline", "GwasSumStats",
  function(data,
           methods,
           addSusieInf       = TRUE,
           coverage          = 0.95,
           secondaryCoverage = c(0.7, 0.5),
           signalCutoff      = 0.025,
           minAbsCorr        = 0.8,
           medianAbsCorr     = NULL,
           fineMappingResult = NULL,
           verbose           = 1,
           trim              = TRUE,
           ...) {
    .fmAssertQcd(data)
    norm       <- .fmNormalizeMethods(methods)
    tokens     <- norm$tokens
    methodArgs <- norm$methodArgs
    .fmCheckMethodCapabilities(tokens, "GwasSumStats")

    # Per the design contract, one GwasSumStats represents the GWAS
    # sumstats for a single LD block (the user is responsible for
    # building one collection per block when sweeping the genome). We
    # iterate per study row and fine-map each (study, method) tuple
    # across the entire entry's variant set; no in-pipeline LD-block
    # partitioning.
    ldSketch <- getLdSketch(data)
    studyCol <- as.character(data$study)

    rowStudy   <- character(0)
    rowMethod  <- character(0)
    rowRegion  <- character(0)
    rowEntries <- list()
    pushRow <- function(st, mt, rg, ent) {
      rowStudy   <<- c(rowStudy,   st)
      rowMethod  <<- c(rowMethod,  mt)
      rowRegion  <<- c(rowRegion,  rg)
      rowEntries[[length(rowEntries) + 1L]] <<- ent
    }

    for (i in seq_len(nrow(data))) {
      st <- studyCol[[i]]
      gr <- data$entry[[i]]
      zn <- .fmExtractZN(gr,
        sprintf("fineMappingPipeline(GwasSumStats): study='%s'", st))
      variantIds <- zn$variantIds
      z <- zn$z
      n <- zn$n
      # Effect-allele frequency for export as the `af` column (post-QC the
      # entry carries the harmonized, complemented frequency in its MAF mcol).
      # Aligned to `variantIds`; NULL when absent -> af exported as NA.
      .gmc <- S4Vectors::mcols(gr)
      afByVar <- if ("MAF" %in% colnames(.gmc))
        setNames(as.numeric(.gmc$MAF), as.character(.gmc$SNP))[variantIds] else NULL
      # Derive a region_id from the entry's GRanges so multi-block
      # genome-wide GWAS sweeps can carry one row per block without
      # tripping (study, method, region_id) uniqueness. Format:
      # "{seqname}_{minPos}_{maxPos}" (e.g. "chr22_10516173_17379581").
      region_id <- sprintf("%s_%d_%d",
                           as.character(GenomicRanges::seqnames(gr))[[1L]],
                           min(GenomicRanges::start(gr)),
                           max(GenomicRanges::start(gr)))

      # The .fmCacheLookup helper takes a 4-tuple key for the QTL cache
      # shape. For GWAS resume we look up using the GwasFineMappingResult
      # 3-tuple shape (study, method, region_id).
      toRun <- character(0)
      for (tk in tokens) {
        cached <- if (!is.null(fineMappingResult) &&
                      is(fineMappingResult, "GwasFineMappingResult")) {
          .fmCacheLookupGwas(fineMappingResult, st, tk, region_id)
        } else NULL
        if (!is.null(cached)) {
          pushRow(st, tk, region_id, cached)
        } else {
          toRun <- c(toRun, tk)
        }
      }
      if (length(toRun) == 0L) next

      ldMat <- .fmLdFromSketch(ldSketch, variantIds)
      names(z) <- variantIds
      chainLocal <- .fmResolveSusieChain(toRun, addSusieInf)
      infFit <- NULL
      if (chainLocal$runInf) {
        if (verbose >= 1)
          message(sprintf("Fitting susieInf (RSS) for GWAS (study='%s', region='%s') ...",
                          st, region_id))
        infFit <- .fmFitSusieRss(z, ldMat, n, "susieInf",
                                 coverage = coverage,
                                 userArgs = methodArgs[["susieInf"]])
      }
      for (tk in toRun) {
        if (tk == "susieInf") {
          if (!chainLocal$keepInf) next
          fit <- infFit
        } else {
          chainFrom <- if ((tk == "susie"    && chainLocal$chainSusie) ||
                           (tk == "susieAsh" && chainLocal$chainAsh))
                        infFit else NULL
          if (verbose >= 1)
            message(sprintf("Fitting %s (RSS) for GWAS (study='%s', region='%s') ...",
                            tk, st, region_id))
          fit <- .fmFitSusieRss(z, ldMat, n, tk,
                                chainFromInf = chainFrom,
                                coverage = coverage,
                                userArgs = methodArgs[[tk]])
        }
        ent <- .fmPostprocessOne(
          fit = fit, method = "susieRss",
          dataX = ldMat, dataY = list(z = z),
          coverage = coverage,
          secondaryCoverage = secondaryCoverage,
          signalCutoff = signalCutoff,
          minAbsCorr = minAbsCorr,
          af = afByVar,
          csInput = "Xcorr")
        pushRow(st, tk, region_id, ent)
      }
    }

    .fmBuildGwasResult(rowStudy, rowMethod, rowEntries,
                       region_ids = rowRegion,
                       ldSketch   = ldSketch)
  })


# =============================================================================
# ANY fallback
# =============================================================================

#' @rdname fineMappingPipeline
#' @export
setMethod("fineMappingPipeline", "ANY",
  function(data, ...) {
    stop("fineMappingPipeline does not accept inputs of class '",
         class(data)[[1L]], "'. Pass a QtlDataset, MultiStudyQtlDataset, ",
         "QtlSumStats, or GwasSumStats. Use summaryStatsQc() on SumStats ",
         "inputs first.")
  })
