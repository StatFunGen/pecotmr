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
#'   \item \code{pipCutoffToSkip} pre-screen: this lived in the old
#'         pipelines but is not ported. Callers can run a one-shot
#'         \code{susieR::susie} check externally if needed.
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
#'   \code{traitId} is supplied.
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
           jointSpecification = NULL,
           addSusieInf        = TRUE,
           coverage           = 0.95,
           secondaryCoverage  = c(0.7, 0.5),
           signalCutoff       = 0.025,
           minAbsCorr         = 0.8,
           medianAbsCorr      = NULL,
           fineMappingResult  = NULL,
           naAction           = c("drop", "impute"),
           verbose            = 1,
           trim               = TRUE,
           phenotypeCovariatesToResidualize = NULL,
           genotypeCovariatesToResidualize  = NULL,
           residualizePhenotypeCovariates   = TRUE,
           residualizeGenotypeCovariates    = TRUE,
           ...) {
    naAction <- match.arg(naAction)
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
        methodArgs = methodArgs)
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
          X <- .fmResidGeno(
            data, contexts = ctx, traitId = tid,
            cisWindow = cisWindow, samples = rownames(Y))
          common <- intersect(rownames(X), rownames(Y))
          if (length(common) < 2L) {
            stop(sprintf(
              "fineMappingPipeline: too few shared samples between residualized X and Y for (context='%s', trait='%s').",
              ctx, tid))
          }
          X <- X[common, , drop = FALSE]
          y <- Y[common, , drop = FALSE]
          if (ncol(y) > 1L) y <- y[, 1L, drop = TRUE] else y <- drop(y)

          chainLocal <- .fmResolveSusieChain(toRun, addSusieInf)
          infFit <- NULL
          if (chainLocal$runInf) {
            if (verbose >= 1)
              message(sprintf("Fitting susieInf for (context='%s', trait='%s') ...", ctx, tid))
            infFit <- .fmFitSusieIndiv(X, y, "susieInf",
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
                message(sprintf("Fitting %s for (context='%s', trait='%s') ...",
                                tk, ctx, tid))
              fit <- .fmFitSusieIndiv(X, y, tk,
                                     chainFromInf = chainFrom,
                                     coverage = coverage,
                                     userArgs = methodArgs[[tk]])
            }
            entry <- .fmPostprocessOne(
              fit = fit, method = tk,
              dataX = X, dataY = y,
              coverage = coverage,
              secondaryCoverage = secondaryCoverage,
              signalCutoff = signalCutoff,
              minAbsCorr = minAbsCorr,
              csInput = "X")
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
          # Build Y matrix per context for this single trait. Sample
          # intersection across contexts; phenotypeCovariates differ per
          # context but getResidualizedPhenotypes already residualises.
          contextsHere <- job$contexts
          # Use the union of per-context cis-windows for variant extraction.
          Yres <- .fmResidPheno(
            data, contexts = contextsHere, traitId = tid, naAction = naAction)
          if (length(contextsHere) == 1L)
            Yres <- setNames(list(Yres), contextsHere)
          commonSamples <- Reduce(intersect, lapply(Yres, rownames))
          X <- .fmResidGeno(
            data, contexts = contextsHere, traitId = tid,
            cisWindow = cisWindow, samples = commonSamples)
          commonSamples <- intersect(commonSamples, rownames(X))
          if (length(commonSamples) < 2L) {
            stop("fineMappingPipeline(QtlDataset, mvsusie multi-context): ",
                 "insufficient shared samples across selected contexts.")
          }
          X <- X[commonSamples, , drop = FALSE]
          Y <- do.call(cbind, lapply(contextsHere, function(ctx) {
            ym <- Yres[[ctx]][commonSamples, , drop = FALSE]
            colnames(ym) <- ctx
            ym
          }))

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

          if (verbose >= 1)
            message(sprintf("Fitting mvsusie (multi-context) for trait='%s' ...", tid))
          mvBaseArgs <- list(
            X = X, Y = Y,
            prior_variance = mvsusieR::create_mixture_prior(R = ncol(Y)),
            coverage = coverage)
          fit <- do.call(fitMvsusie,
                         .fmMergeUserArgs(mvBaseArgs, "mvsusie",
                                          methodArgs[["mvsusie"]]))
          fit <- .setFinemappingFitClass(fit, "mvsusie")
          entry <- .fmPostprocessOne(
            fit = fit, method = "mvsusie",
            dataX = X, dataY = NULL,
            coverage = coverage,
            secondaryCoverage = secondaryCoverage,
            signalCutoff = signalCutoff,
            minAbsCorr = minAbsCorr,
            csInput = "X")
          # Share the joint fit across contexts via copy-on-modify.
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
          X <- .fmResidGeno(
            data, contexts = ctx, traitId = traits,
            cisWindow = cisWindow, samples = rownames(Y))
          common <- intersect(rownames(X), rownames(Y))
          if (length(common) < 2L) {
            stop(sprintf(
              "fineMappingPipeline(QtlDataset, mvsusie multi-trait): too few shared samples in context '%s'.",
              ctx))
          }
          X <- X[common, , drop = FALSE]
          Y <- Y[common, , drop = FALSE]

          if (verbose >= 1)
            message(sprintf("Fitting mvsusie (multi-trait) for context='%s' ...", ctx))
          mvBaseArgs <- list(
            X = X, Y = Y,
            prior_variance = mvsusieR::create_mixture_prior(R = ncol(Y)),
            coverage = coverage)
          fit <- do.call(fitMvsusie,
                         .fmMergeUserArgs(mvBaseArgs, "mvsusie",
                                          methodArgs[["mvsusie"]]))
          fit <- .setFinemappingFitClass(fit, "mvsusie")
          entry <- .fmPostprocessOne(
            fit = fit, method = "mvsusie",
            dataX = X, dataY = NULL,
            coverage = coverage,
            secondaryCoverage = secondaryCoverage,
            signalCutoff = signalCutoff,
            minAbsCorr = minAbsCorr,
            csInput = "X")
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
        X <- .fmResidGeno(
          data, contexts = ctx, traitId = traits,
          cisWindow = cisWindow, samples = rownames(Y))
        common <- intersect(rownames(X), rownames(Y))
        if (length(common) < 2L) {
          stop(sprintf("fineMappingPipeline(QtlDataset, fsusie): too few shared samples in context '%s'.", ctx))
        }
        X <- X[common, , drop = FALSE]
        Y <- Y[common, , drop = FALSE]

        # Per-trait genomic positions for the wavelet model. Use the
        # midpoint of each trait's rowRanges in this context.
        se <- getPhenotypes(data, contexts = ctx, traitId = traits)
        rr <- SummarizedExperiment::rowRanges(se)
        # Reorder rr to the column order of Y.
        rrIds <- rownames(se)
        ord <- match(colnames(Y), rrIds)
        if (anyNA(ord)) {
          stop("fineMappingPipeline(QtlDataset, fsusie): unable to align trait positions to Y columns.")
        }
        rr <- rr[ord]
        pos <- (GenomicRanges::start(rr) + GenomicRanges::end(rr)) / 2

        if (verbose >= 1)
          message(sprintf("Fitting fsusie for context='%s' (multi-trait, %d traits) ...",
                          ctx, length(traits)))
        fit <- do.call(fitFsusie,
                       .fmMergeUserArgs(list(X = X, Y = Y, pos = pos),
                                        "fsusie", methodArgs[["fsusie"]]))
        fit <- .setFinemappingFitClass(fit, "fsusie")
        entry <- .fmPostprocessOne(
          fit = fit, method = "fsusie",
          dataX = X, dataY = NULL,
          coverage = coverage,
          secondaryCoverage = secondaryCoverage,
          signalCutoff = signalCutoff,
          minAbsCorr = minAbsCorr,
          csInput = "fsusie")
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
           jointSpecification = NULL,
           addSusieInf        = TRUE,
           coverage           = 0.95,
           secondaryCoverage  = c(0.7, 0.5),
           signalCutoff       = 0.025,
           minAbsCorr         = 0.8,
           medianAbsCorr      = NULL,
           fineMappingResult  = NULL,
           naAction           = c("drop", "impute"),
           verbose            = 1,
           trim               = TRUE,
           phenotypeCovariatesToResidualize = NULL,
           genotypeCovariatesToResidualize  = NULL,
           residualizePhenotypeCovariates   = TRUE,
           residualizeGenotypeCovariates    = TRUE,
           ...) {
    naAction <- match.arg(naAction)
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
        methodArgs = methodArgs)
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
        jointSpecification = NULL,
        addSusieInf        = addSusieInf,
        coverage           = coverage,
        secondaryCoverage  = secondaryCoverage,
        signalCutoff       = signalCutoff,
        minAbsCorr         = minAbsCorr,
        fineMappingResult  = fineMappingResult,
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
