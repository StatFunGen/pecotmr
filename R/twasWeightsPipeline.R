# =============================================================================
# Helpers + S4 dispatch surface for twasWeightsPipeline
# =============================================================================

# Concatenate two TwasWeights collections row-wise. `rbind` on DFrame
# subclasses does not reliably preserve the `ldSketch` slot, so this
# helper rebuilds via the constructor. Optional joint columns
# (`jointStudies`, `jointContexts`, `jointTraits`) are carried through
# via .combineJointCol() so a mixed rbind of joint + non-joint rows
# pads the non-joint side with NA_character_.
# @noRd
.rbindTwasWeights <- function(a, b, ldSketch = NULL) {
  if (!is(a, "TwasWeights") || !is(b, "TwasWeights")) {
    stop(".rbindTwasWeights expects two TwasWeights inputs.")
  }
  TwasWeights(
    study         = c(as.character(a$study),   as.character(b$study)),
    context       = c(as.character(a$context), as.character(b$context)),
    trait         = c(as.character(a$trait),   as.character(b$trait)),
    method        = c(as.character(a$method),  as.character(b$method)),
    entry         = c(as.list(a$entry), as.list(b$entry)),
    jointStudies  = .combineJointCol(a, b, "jointStudies"),
    jointContexts = .combineJointCol(a, b, "jointContexts"),
    jointTraits   = .combineJointCol(a, b, "jointTraits"),
    ldSketch      = ldSketch)
}

# --- Multi-region (jointRegions) helpers for the QtlDataset method ----------

# Label a region block for per-region reporting: the genomic coordinate of a
# single-range window, or "cis" for the trait-derived (region = NULL) block.
.twasRegionLabel <- function(rg) {
  if (is.null(rg)) return("cis")
  paste0(as.character(GenomicRanges::seqnames(rg))[[1L]], ":",
         GenomicRanges::start(rg)[[1L]], "-", GenomicRanges::end(rg)[[1L]])
}

# Select the per-region fine-mapping fits for region block `i`. A
# jointRegions=FALSE multi-region fine-mapping stores its per-region SuSiE fits
# as a named list (region1, region2, ...); pick the matching element. With a
# single block the fits are returned unchanged; a non-region-list fit under
# multiple blocks cannot be aligned and is dropped (the method learns fresh).
.twasFitsForRegion <- function(fits, i, nBlocks) {
  if (length(fits) == 0L || nBlocks == 1L) return(fits)
  out <- lapply(fits, function(f) {
    if (is.list(f) && !is.null(names(f)) && length(f) > 0L &&
        all(nzchar(names(f))) && all(startsWith(names(f), "region"))) {
      if (i <= length(f)) f[[i]] else NULL
    } else {
      NULL
    }
  })
  out[!vapply(out, is.null, logical(1))]
}

# Flat per-region cvPerformance reporting table: one row per region carrying the
# region label plus that region's CV metric columns. Per-sample predictions are
# intentionally omitted — this is a summary-reporting structure.
.twasRegionCvDf <- function(entries, regionLabels) {
  rows <- Map(function(e, lab) {
    cv <- getCvPerformance(e)
    if (is.null(cv) || is.null(cv$metrics)) return(NULL)
    cbind(data.frame(region = lab, stringsAsFactors = FALSE),
          as.data.frame(as.list(cv$metrics), check.names = FALSE))
  }, entries, regionLabels)
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) return(NULL)
  do.call(rbind, rows)
}

# Concatenate one method's per-region TwasWeightsEntry payloads into a single
# entry. Variants/weights are stacked (regions are disjoint), the per-region
# fits are kept as a named list, and cvPerformance becomes the flat per-region
# reporting data.frame.
.twasMergeRegionEntries <- function(entries, regionLabels) {
  keep <- !vapply(entries, is.null, logical(1))
  entries <- entries[keep]; regionLabels <- regionLabels[keep]
  if (length(entries) == 0L) return(NULL)
  if (length(entries) == 1L) return(entries[[1L]])
  wList <- lapply(entries, getWeights)
  weights <- if (is.matrix(wList[[1L]])) do.call(rbind, wList)
             else unlist(wList, use.names = FALSE)
  TwasWeightsEntry(
    variantIds    = unlist(lapply(entries, getVariantIds), use.names = FALSE),
    weights       = weights,
    fits          = setNames(lapply(entries, getFits), regionLabels),
    cvPerformance = .twasRegionCvDf(entries, regionLabels),
    standardized  = getStandardized(entries[[1L]]),
    dataType      = getDataType(entries[[1L]]))
}

# Merge per-region TwasWeights collections (same study/context/trait, same
# methods) into one collection by concatenating each method's entry.
.twasMergeRegions <- function(twList, regionLabels) {
  keep <- !vapply(twList, is.null, logical(1))
  twList <- twList[keep]; regionLabels <- regionLabels[keep]
  if (length(twList) == 0L) return(NULL)
  if (length(twList) == 1L) return(twList[[1L]])
  base <- twList[[1L]]
  mergedEntries <- lapply(seq_along(base$method), function(r) {
    key <- c(as.character(base$study[[r]]),   as.character(base$context[[r]]),
             as.character(base$trait[[r]]),    as.character(base$method[[r]]))
    perRegion <- lapply(twList, function(tw) {
      hit <- which(as.character(tw$study)   == key[[1L]] &
                   as.character(tw$context) == key[[2L]] &
                   as.character(tw$trait)   == key[[3L]] &
                   as.character(tw$method)  == key[[4L]])
      if (length(hit)) tw$entry[[hit[[1L]]]] else NULL
    })
    .twasMergeRegionEntries(perRegion, regionLabels)
  })
  TwasWeights(
    study   = as.character(base$study),
    context = as.character(base$context),
    trait   = as.character(base$trait),
    method  = as.character(base$method),
    entry   = mergedEntries)
}

# Splice per-(method, outcome) cross-validated predictions and the 6-metric
# performance row from a `twasWeightsCv()` result into the matching
# `TwasWeightsEntry$cvPerformance` slot of every row in a TwasWeights
# collection. Rebuilds the collection because TwasWeightsEntry is treated
# as immutable. Rows for which no CV result is available (method not in
# the CV run, or trait not in the CV prediction matrix's columns) are
# emitted unchanged.
#
# The CV result keys carry a method suffix (`<m>_predicted`,
# `<m>_performance` in snake form, or `<m>Predicted`, `<m>Performance` in
# camel form); the TwasWeights `method` column carries the bare token
# (e.g. "lasso"). The trait column carries the outcome name, which must
# match the column name of the CV prediction matrix.
# @noRd
.spliceCvIntoTwasWeights <- function(twasWeights, twasCvResult,
                                      ldSketch = NULL) {
  if (is.null(twasCvResult) || is.null(twasCvResult$prediction) ||
      is.null(twasCvResult$performance)) {
    return(twasWeights)
  }
  predKeyBase <- sub("(_predicted|Predicted)$", "",
                     names(twasCvResult$prediction))
  perfKeyBase <- sub("(_performance|Performance)$", "",
                     names(twasCvResult$performance))

  pickKey <- function(bare, keys, base) {
    hit <- which(base == bare)
    if (length(hit) == 0L) NA_character_ else keys[[hit[[1L]]]]
  }

  studies   <- as.character(twasWeights$study)
  contexts  <- as.character(twasWeights$context)
  traits    <- as.character(twasWeights$trait)
  methodsV  <- as.character(twasWeights$method)
  newEntries <- as.list(twasWeights$entry)

  for (i in seq_along(newEntries)) {
    bare <- methodsV[[i]]
    pKey <- pickKey(bare, names(twasCvResult$prediction), predKeyBase)
    mKey <- pickKey(bare, names(twasCvResult$performance), perfKeyBase)
    if (is.na(pKey) || is.na(mKey)) next
    predMat <- twasCvResult$prediction[[pKey]]
    perfMat <- twasCvResult$performance[[mKey]]
    if (is.null(predMat) || is.null(perfMat)) next

    tr  <- traits[[i]]
    predCols <- colnames(predMat)
    perfRows <- rownames(perfMat)
    colHit <- if (!is.null(predCols) && tr %in% predCols) tr
              else if (ncol(predMat) == 1L) 1L else NA_integer_
    rowHit <- if (!is.null(perfRows) && tr %in% perfRows) tr
              else if (nrow(perfMat) == 1L) 1L else NA_integer_
    if (is.na(colHit) || is.na(rowHit)) next

    predVec <- predMat[, colHit, drop = TRUE]
    metRow  <- perfMat[rowHit, , drop = TRUE]
    cv <- list(
      samplePartition = twasCvResult$samplePartition,
      predictions     = predVec,
      metrics         = metRow)
    entry <- newEntries[[i]]
    newEntries[[i]] <- TwasWeightsEntry(
      variantIds    = getVariantIds(entry),
      weights       = getWeights(entry),
      fits          = getFits(entry),
      cvPerformance = cv,
      standardized  = getStandardized(entry),
      dataType      = getDataType(entry))
  }

  TwasWeights(
    study    = studies,
    context  = contexts,
    trait    = traits,
    method   = methodsV,
    entry    = newEntries,
    ldSketch = ldSketch)
}

# Mapping from short / canonical TWAS weight-method name to dispatch
# capability. Used to reject incompatible (input class, method) pairs.
#
# `allowsIndiv`  : may be invoked on a QtlDataset (individual-level X, Y).
# `allowsRss`    : may be invoked on a QtlSumStats / GwasSumStats (RSS).
# `multivariate` : requires a multi-trait / multi-context Y (mvsusie /
#                  mr.mash family).
#
# Rules from `dev/refactor-design.md` (`twasWeightsPipeline` row):
# - PRS-CS is RSS-only.
# - BGLR / CRAN-stable qgg methods (bayes_a/b/c/l/n/r, b_lasso, dpr_*)
#   are individual-level only.
# - mr.mash / mvsusie follow the multi-trait / multi-context rules of
#   the mvSuSiE fine-mapping family.
# @noRd
# User-facing TWAS method tokens are unified across input classes;
# auto-dispatch picks the individual-level vs sumstat implementation based
# on the QtlDataset / QtlSumStats input. Each entry records:
#   individualImpl  Function name to call on QtlDataset input (NULL = not
#                   supported on individual-level input).
#   sumstatImpl     Function name to call on QtlSumStats input (NULL = not
#                   supported on sumstat input).
#   multivariate    Whether the method requires multi-trait / multi-context
#                   structure (mvsusie / mrmash / mvsusieRss / mrmashRss).
#
# Per the design: BGLR / qgg "Bayes alphabet" methods (bayes_a/b/c/l/n/r,
# b_lasso) are individual-only until the qgg CRAN release adds qBayes
# sumstat support. dpr_gibbs has the SDPR sumstat counterpart;
# dpr_vb / dpr_adaptive_gibbs remain individual-only. enet has no cpp11
# sumstat solver yet (lassosumRssRcpp is pure L1, no alpha mixing) and is
# documented as individual-only for now. prsCs has no individual-level
# counterpart (it is a sumstat-only Bayesian shrinkage method).
.twasMethodCapabilities <- list(
  # NOTE: fine-mapping methods (susie / susieInf / susieAsh / mvsusie /
  # fsusie) are NOT listed here. Their availability is governed by
  # .fineMappingMethodCapabilities (the same registry fineMappingPipeline
  # uses) and gated by .twasCheckFineMappingMethods, which delegates
  # input-class compatibility to .fmCheckMethodCapabilities.
  mrash               = list(individualImpl = "mrashWeights",
                             sumstatImpl    = "mrAshRssWeights",
                             multivariate   = FALSE),
  lasso               = list(individualImpl = "lassoWeights",
                             sumstatImpl    = "lassosumRssWeights",
                             multivariate   = FALSE),
  scad                = list(individualImpl = "scadWeights",
                             sumstatImpl    = "scadRssWeights",
                             multivariate   = FALSE),
  mcp                 = list(individualImpl = "mcpWeights",
                             sumstatImpl    = "mcpRssWeights",
                             multivariate   = FALSE),
  l0learn             = list(individualImpl = "l0learnWeights",
                             sumstatImpl    = "l0learnRssWeights",
                             multivariate   = FALSE),
  mrmash              = list(individualImpl = "mrmashWeights",
                             sumstatImpl    = "mrmashRssWeights",
                             multivariate   = TRUE),
  dpr_gibbs           = list(individualImpl = "dprGibbsWeights",
                             sumstatImpl    = "sdprWeights",
                             multivariate   = FALSE),
  # Individual-only — no cpp11 sumstat solver yet.
  enet                = list(individualImpl = "enetWeights",
                             sumstatImpl    = NULL,
                             multivariate   = FALSE),
  # Individual-only DPR variants (sumstat counterparts not implemented).
  dpr_vb              = list(individualImpl = "dprVbWeights",
                             sumstatImpl    = NULL,
                             multivariate   = FALSE),
  dpr_adaptive_gibbs  = list(individualImpl = "dprAdaptiveGibbsWeights",
                             sumstatImpl    = NULL,
                             multivariate   = FALSE),
  # qgg Bayes alphabet — individual-only until qgg CRAN release.
  bayes_a             = list(individualImpl = "bayesAWeights",
                             sumstatImpl    = NULL,
                             multivariate   = FALSE),
  bayes_b             = list(individualImpl = "bayesBWeights",
                             sumstatImpl    = NULL,
                             multivariate   = FALSE),
  bayes_c             = list(individualImpl = "bayesCWeights",
                             sumstatImpl    = NULL,
                             multivariate   = FALSE),
  bayes_l             = list(individualImpl = "bLassoWeights",
                             sumstatImpl    = NULL,
                             multivariate   = FALSE),
  bayes_n             = list(individualImpl = "bayesNWeights",
                             sumstatImpl    = NULL,
                             multivariate   = FALSE),
  bayes_r             = list(individualImpl = "bayesRWeights",
                             sumstatImpl    = NULL,
                             multivariate   = FALSE),
  b_lasso             = list(individualImpl = "bLassoWeights",
                             sumstatImpl    = NULL,
                             multivariate   = FALSE),
  # Sumstat-only Bayesian shrinkage (no individual-level analogue).
  prsCs               = list(individualImpl = NULL,
                             sumstatImpl    = "prsCsWeights",
                             multivariate   = FALSE))

# Normalize a user-supplied `methods` argument (character vector, preset
# string, or named list per `.twasMethodLookup`) into a (token, args) pair
# suitable for `.twasWeightsPipelineMatrix` / the sumstat sub-pipelines.
# Returns a list with `tokens` (canonical short names, used for capability
# lookup) and `methodList` (the `<token>_weights = args` list passed to
# `learnTwasWeights` / sumstat helpers).
# @noRd
.twasNormalizeMethods <- function(methods) {
  if (is.null(methods)) {
    methodList <- .twasMethodLookup("default")
    tokens <- .twasTokensFromMethodList(methodList)
  } else if (is.character(methods)) {
    # Fine-mapping tokens without a .twasMethodLookup entry (e.g. fsusie)
    # are recognised here so the downstream gate can produce a method-
    # specific error rather than "Unknown TWAS method" from the lookup.
    fmExtra <- setdiff(intersect(methods, .twasFineMappingTokens()),
                       .twasKnownMethodLookupNames())
    regular <- setdiff(methods, fmExtra)
    methodList <- if (length(regular) > 0L) .twasMethodLookup(regular)
                  else list()
    if (length(fmExtra) > 0L) {
      # Append stub entries (empty args) for fine-mapping tokens with no
      # learner counterpart (e.g. fsusie). The gate will reject these.
      for (tk in fmExtra) {
        snake <- paste0(tk, "_weights")
        methodList[[snake]] <- list()
      }
    }
    # Tokens come from the user input (canonical camelCase) — the snake
    # keys in methodList are an internal detail of learnTwasWeights.
    tokens <- unique(methods)
  } else if (is.list(methods)) {
    methodList <- methods
    tokens <- .twasTokensFromMethodList(methodList)
  } else {
    stop("`methods` must be a character vector, preset string, or named list.")
  }
  list(tokens = tokens, methodList = methodList)
}

# Canonical (camelCase) tokens known to .twasMethodLookup, for use by
# .twasNormalizeMethods. Source of truth: the methodMap inside
# .twasMethodLookup.
# @noRd
.twasKnownMethodLookupNames <- function() {
  c("susie", "susieAsh", "susieInf", "mrash", "enet", "lasso",
    "bayes_r", "bayes_l", "bayes_a", "bayes_b", "bayes_c", "bayes_n",
    "b_lasso", "dpr_vb", "dpr_gibbs", "dpr_adaptive_gibbs",
    "scad", "mcp", "l0learn", "mvsusie", "mrmash")
}

# Convert a methodList (snake_case keys like `susie_inf_weights`) back to
# canonical camelCase tokens (susieInf). Falls back to the snake form for
# unknown keys.
# @noRd
.twasTokensFromMethodList <- function(methodList) {
  snake <- sub("(_weights|Weights)$", "", names(methodList))
  snakeToCanonical <- c(
    susie = "susie", susie_ash = "susieAsh", susie_inf = "susieInf",
    susie_ash_inf = "susieAsh",
    mrash = "mrash", enet = "enet", lasso = "lasso",
    bayes_r = "bayes_r", bayes_l = "bayes_l", bayes_a = "bayes_a",
    bayes_b = "bayes_b", bayes_c = "bayes_c", bayes_n = "bayes_n",
    b_lasso = "b_lasso", dpr_vb = "dpr_vb", dpr_gibbs = "dpr_gibbs",
    dpr_adaptive_gibbs = "dpr_adaptive_gibbs",
    scad = "scad", mcp = "mcp", l0learn = "l0learn",
    mvsusie = "mvsusie", mrmash = "mrmash", prsCs = "prsCs",
    fsusie = "fsusie")
  vapply(snake, function(s) {
    if (!is.na(snakeToCanonical[s])) snakeToCanonical[[s]] else s
  }, character(1), USE.NAMES = FALSE)
}

# Enforce input-class / method compatibility against the TWAS
# capability table. Routes the input class through individual /
# sumstat branches; the twasWeightsPipeline has no GwasSumStats input
# path so that branch is omitted. Emits a single error listing every
# offending token.
# @noRd
.twasCheckMethodCapabilities <- function(tokens, inputKind) {
  if (length(tokens) == 0L) return(invisible(NULL))
  caps <- .twasMethodCapabilities
  # Fine-mapping tokens are governed by .twasCheckFineMappingMethods (and
  # delegate input-class compat to .fmCheckMethodCapabilities); skip them
  # here so they aren't reported as "unknown".
  fmTokens <- intersect(tokens, .twasFineMappingTokens())
  tokens   <- setdiff(tokens, fmTokens)
  if (length(tokens) == 0L) return(invisible(NULL))
  unknown <- setdiff(tokens, names(caps))
  if (length(unknown) > 0L) {
    stop(sprintf(
      "twasWeightsPipeline: unknown method token(s): %s. Known tokens: %s.",
      paste(unknown, collapse = ", "),
      paste(c(names(caps), .twasFineMappingTokens()), collapse = ", ")))
  }
  individualKinds <- c("QtlDataset", "MultiStudyQtlDataset")
  bad <- character(0); reason <- character(0)
  for (tk in tokens) {
    info <- caps[[tk]]
    if (inputKind %in% individualKinds) {
      if (is.null(info$individualImpl)) {
        bad <- c(bad, tk)
        reason <- c(reason, "is sumstat-only (use a QtlSumStats input)")
      }
    } else if (inputKind == "QtlSumStats") {
      if (is.null(info$sumstatImpl)) {
        bad <- c(bad, tk)
        reason <- c(reason, "is individual-only (use a QtlDataset input)")
      }
    }
    # twasWeightsPipeline does not support GwasSumStats input.
  }
  if (length(bad) > 0L) {
    stop(sprintf(
      "twasWeightsPipeline: the following method(s) are not available for input class '%s': %s. %s.",
      inputKind,
      paste(unique(bad), collapse = ", "),
      paste(sprintf("%s %s", bad, reason), collapse = "; ")))
  }
}

# Adapter registry mapping each fine-mapping method (whose existence is
# governed by .fineMappingMethodCapabilities) to its TWAS-weight extractor
# wrapper. The wrapper names follow the *Weights / *RssWeights convention,
# and the *Fit argument receives the pre-fitted fine-mapping object.
# fSuSiE is multivariate (it collapses a functional fit to a variants x
# features weight matrix via fsusieWeights) and has no RSS counterpart.
# @noRd
.twasFineMappingMethodAdapters <- list(
  susie    = list(weightFn = "susieWeights",
                  rssWeightFn = "susieRssWeights",
                  fitArg = "susieFit",
                  rssFitArg = "susieRssFit",
                  methodKey = "susie_weights"),
  susieInf = list(weightFn = "susieInfWeights",
                  rssWeightFn = "susieInfRssWeights",
                  fitArg = "susieInfFit",
                  rssFitArg = "susieInfRssFit",
                  methodKey = "susie_inf_weights"),
  susieAsh = list(weightFn = "susieAshWeights",
                  rssWeightFn = "susieAshRssWeights",
                  fitArg = "susieAshFit",
                  rssFitArg = "susieAshRssFit",
                  methodKey = "susie_ash_weights"),
  mvsusie  = list(weightFn = "mvsusieWeights",
                  rssWeightFn = "mvsusieRssWeights",
                  fitArg = "mvsusieFit",
                  rssFitArg = "mvsusieRssFit",
                  methodKey = "mvsusie_weights"),
  fsusie   = list(weightFn = "fsusieWeights",
                  rssWeightFn = NULL,
                  fitArg = "fsusieFit",
                  rssFitArg = NULL,
                  methodKey = "fsusie_weights"))

# Canonical list of fine-mapping tokens recognised by twasWeightsPipeline.
# Sourced from fineMappingPipeline's registry minus mrmash (which
# fineMappingPipeline hard-rejects as a TWAS-only method).
# @noRd
.twasFineMappingTokens <- function() {
  setdiff(names(.fineMappingMethodCapabilities), "mrmash")
}

# Reject fine-mapping methods (susie / susieInf / susieAsh / mvsusie /
# fsusie) when no FineMappingResult is supplied. twasWeightsPipeline is
# not allowed to re-fit fine-mapping models from scratch; users must run
# fineMappingPipeline() first and pass the result via `fineMappingResult`.
# Input-class compatibility (e.g. fsusie has no QtlSumStats path) is
# delegated to .fmCheckMethodCapabilities so the rule set stays in lock-
# step with fineMappingPipeline. Methods with no TWAS-weight extractor
# (fsusie) are rejected with a method-specific message.
# @noRd
.twasCheckFineMappingMethods <- function(tokens, fineMappingResult, inputKind) {
  if (length(tokens) == 0L) return(invisible(NULL))
  fmTokens <- intersect(tokens, .twasFineMappingTokens())
  if (length(fmTokens) == 0L) return(invisible(NULL))

  # Defer input-class compatibility to fineMappingPipeline. e.g. this
  # rejects fsusie on QtlSumStats (fsusie has no RSS impl).
  .fmCheckMethodCapabilities(fmTokens, inputKind)

  # Reject fine-mapping methods that have no TWAS-weight extractor
  # (currently only fsusie).
  noAdapter <- setdiff(fmTokens, names(.twasFineMappingMethodAdapters))
  if (length(noAdapter) > 0L) {
    stop(sprintf(
      "twasWeightsPipeline: method(s) %s have no TWAS-weight extractor. For multi-trait fine-mapping use mvsusie via fineMappingResult.",
      paste(noAdapter, collapse = ", ")))
  }

  if (is.null(fineMappingResult)) {
    stop(sprintf(
      "twasWeightsPipeline: method(s) %s are fine-mapping methods and may not be re-fit by twasWeightsPipeline. Run fineMappingPipeline() first and pass the result via `fineMappingResult = <FineMappingResult>`.",
      paste(unique(fmTokens), collapse = ", ")))
  }
  if (!is(fineMappingResult, "FineMappingResultBase")) {
    stop("`fineMappingResult` must be a FineMappingResult or NULL.")
  }
  invisible(NULL)
}

# Look up the multivariate flag for a token. Checks the TWAS-regression
# capability table first; if absent, falls back to the fine-mapping
# capability table (the source of truth for susie / mvsusie / fsusie /
# etc.). Returns FALSE for unknown tokens.
# @noRd
.twasIsMultivariateToken <- function(token) {
  info <- .twasMethodCapabilities[[token]]
  if (!is.null(info)) return(isTRUE(info$multivariate))
  fmInfo <- .fineMappingMethodCapabilities[[token]]
  if (!is.null(fmInfo)) return(isTRUE(fmInfo$multivariate))
  FALSE
}

# Enforce the multi-trait / multi-context rule for mvsusie / mr.mash
# methods (same family as the fine-mapping mvSuSiE rule in the design
# doc). Multivariate methods need at least 2 traits *or* at least 2
# contexts in the Y matrix passed to learnTwasWeights.
# @noRd
.twasCheckMultivariateY <- function(tokens, nTraits, nContexts) {
  multivariateTokens <- tokens[vapply(tokens, .twasIsMultivariateToken,
                                       logical(1))]
  if (length(multivariateTokens) == 0L) return(invisible(NULL))
  if (nTraits < 2L && nContexts < 2L) {
    stop(sprintf(
      "twasWeightsPipeline: method(s) %s require multi-trait or multi-context input (got %d trait(s) x %d context(s)).",
      paste(multivariateTokens, collapse = ", "),
      nTraits, nContexts))
  }
}

# Reject SumStats inputs that have not been QC'd via summaryStatsQc.
# @noRd
.twasAssertQcd <- function(sumstats) {
  if (length(getQcInfo(sumstats)) == 0L) {
    stop("twasWeightsPipeline: the supplied ",
         class(sumstats)[[1L]],
         " has no QC record (qcInfo is empty). Call summaryStatsQc() ",
         "first and pass the QC-applied result.")
  }
}

# Extract a correlation matrix from a GenotypeHandle (LD sketch) for the
# variant subset given by `variantIds`. Thin wrapper over the shared
# `.ldFromSketch` helper.
# @noRd
.twasLdFromSketch <- function(ldSketch, variantIds) {
  .ldFromSketch(ldSketch, variantIds, label = ".twasLdFromSketch")
}

# Optional resume-cache lookup for twasWeightsPipeline. Returns the
# matching TwasWeightsEntry from `twasWeights` for the tuple (study,
# context, trait, method), or NULL when there is no hit. Returns NULL
# silently when twasWeights is NULL or not a TwasWeights collection.
# Mirrors .fmCacheLookup (R/fineMappingPipeline.R).
# @noRd
.twasCacheLookup <- function(twasWeights, study, context, trait, method) {
  if (is.null(twasWeights)) return(NULL)
  if (!is(twasWeights, "TwasWeights")) return(NULL)
  idx <- .matchTupleRows(twasWeights,
                          list(study = study, context = context,
                               trait = trait, method = method))
  if (length(idx) == 0L) return(NULL)
  twasWeights$entry[[idx[[1L]]]]
}

# Build a TwasWeights collection from a list of cached entries keyed by
# short-method-name (the value of the `method` column). Helper for the
# resume-cache short-circuit path in twasWeightsPipeline.
# @noRd
.twasBuildFromCachedRows <- function(cachedRows, study, context, trait,
                                     ldSketch = NULL) {
  if (length(cachedRows) == 0L) return(NULL)
  TwasWeights(
    study   = rep(study,   length(cachedRows)),
    context = rep(context, length(cachedRows)),
    trait   = rep(trait,   length(cachedRows)),
    method  = names(cachedRows),
    entry   = unname(cachedRows),
    ldSketch = ldSketch)
}

# Convert a FineMappingResult (single-method susie/susie_inf row matched
# to the requested study/context/trait) into a `fittedModels` list
# suitable for `learnTwasWeights`. Pulls the trimmedFit from the matching
# entry. Returns a (possibly empty) list.
# @noRd
.twasFineMappingFits <- function(fineMappingResult, study, context, trait) {
  if (is.null(fineMappingResult)) return(list())
  if (!is(fineMappingResult, "FineMappingResultBase")) {
    stop("`fineMappingResult` must be a FineMappingResult or NULL.")
  }
  out <- list()
  methods <- as.character(fineMappingResult$method)
  for (canonical in c("susie", "susieInf", "susieAsh", "mvsusie", "fsusie")) {
    candidates <- c(canonical,
                    paste0(tolower(substring(canonical, 1L, 1L)),
                           substring(canonical, 2L)),
                    gsub("([A-Z])", "_\\1", canonical))
    candidates <- tolower(candidates)
    idx <- which(tolower(methods) %in% candidates &
                 as.character(fineMappingResult$study)   == study &
                 as.character(fineMappingResult$context) == context &
                 as.character(fineMappingResult$trait)   == trait)
    if (length(idx) > 0L) {
      out[[canonical]] <- getSusieFit(fineMappingResult$entry[[idx[[1L]]]])
    }
  }
  out
}

# Locate a fine-mapping fit for one (study, context, trait, token) tuple.
# Used by the QtlSumStats sumstat dispatcher to pass the precomputed fit
# into susieRssWeights / susieInfRssWeights / susieAshRssWeights /
# mvsusieRssWeights via their respective *Fit arguments.
# @noRd
.twasFineMappingFitFor <- function(fineMappingResult, study, context, trait,
                                    token) {
  if (is.null(fineMappingResult)) return(NULL)
  fits <- .twasFineMappingFits(fineMappingResult,
                                study = study, context = context, trait = trait)
  fits[[token]]
}

# Collect the cross-validation payload that fineMappingPipeline stored on the
# FineMappingResult for one (study, context, trait) tuple. fineMapping records
# one cvResult per (study, context, trait, method) entry (samplePartition +
# per-fold predictions/metrics, keyed by the TWAS snake method name); this
# merges them across the fine-mapping methods of the tuple into a single
# twasWeightsCv()-shaped list so twasWeightsPipeline can reuse the partition and
# feed those out-of-fold predictions into the SR-TWAS ensemble without re-
# fitting the fine-mapping models. A multi-region entry stores cvResult as a
# per-region list; the first region carrying CV is used. Returns NULL when no
# fine-mapping entry for the tuple recorded CV.
# @noRd
.twasCvResultFor <- function(fineMappingResult, study, context, trait) {
  if (is.null(fineMappingResult)) return(NULL)
  if (!is(fineMappingResult, "FineMappingResultBase")) return(NULL)
  idx <- which(as.character(fineMappingResult$study)   == study &
               as.character(fineMappingResult$context) == context &
               as.character(fineMappingResult$trait)   == trait)
  if (length(idx) == 0L) return(NULL)
  samplePartition <- NULL
  prediction  <- list()
  performance <- list()
  for (i in idx) {
    cv <- getCvResult(fineMappingResult$entry[[i]])
    if (is.null(cv)) next
    # Multi-region entries store cvResult as a named per-region list; pick the
    # first region that carries a partition.
    if (is.null(cv$samplePartition)) {
      hit <- Filter(function(z) is.list(z) && !is.null(z$samplePartition), cv)
      if (length(hit) == 0L) next
      cv <- hit[[1L]]
    }
    if (is.null(samplePartition)) samplePartition <- cv$samplePartition
    prediction  <- c(prediction,  cv$prediction)
    performance <- c(performance, cv$performance)
  }
  if (length(prediction) == 0L) return(NULL)
  list(samplePartition = samplePartition,
       prediction = prediction, performance = performance)
}

#' TWAS Weights Pipeline
#'
#' S4-dispatched per-region pipeline for learning TWAS prediction weights.
#' Accepts:
#' \itemize{
#'   \item a \code{\link{QtlDataset}} for individual-level cohort fits;
#'   \item a \code{\link{QtlSumStats}} for per-trait RSS fits;
#'   \item a \code{\link{GwasSumStats}} for per-LD-block PRS-CS-style fits
#'         from GWAS summary statistics.
#' }
#'
#' Method-restriction rules (enforced):
#' \itemize{
#'   \item \code{mr.mash}, \code{mvsusie} follow the multi-trait /
#'         multi-context rules of the fine-mapping \code{mvsusie} family
#'         (require at least two traits OR at least two contexts).
#'   \item RSS-only methods (PRS-CS, \code{lassosumRss}, SDPR, all
#'         \code{*Rss} variants) are rejected on \code{QtlDataset}
#'         input.
#'   \item Individual-level-only methods (BGLR and CRAN-stable qgg:
#'         \code{bayes_a/b/c/l/n/r}, \code{b_lasso}, \code{dpr_*}) are
#'         rejected on \code{QtlSumStats} / \code{GwasSumStats} input.
#' }
#'
#' Both \code{QtlSumStats} and \code{GwasSumStats} inputs must have been
#' QC'd via \code{\link{summaryStatsQc}} first; otherwise an error is
#' raised pointing at that function.
#'
#' The returned \code{\link{TwasWeights}} collection's \code{ldSketch}
#' slot is set automatically: \code{NULL} for individual-level fits,
#' the input's \code{ldSketch} for RSS-derived fits.
#'
#' Optionally a \code{\link{FineMappingResult}} may be supplied as a
#' source of pre-fit SuSiE / SuSiE-inf / SuSiE-ash objects; their
#' \code{trimmedFit} payloads are passed through to \code{learnTwasWeights}
#' / the RSS sub-pipelines via the \code{fittedModels} slot, avoiding
#' a re-fit.
#'
#' When the supplied \code{FineMappingResult} was produced with
#' cross-validation (\code{fineMappingPipeline(..., cvFolds > 1)}), each
#' matching \code{(study, context, trait)} entry's \code{cvResult} is
#' reused: its fold partition becomes the CV partition (unless
#' \code{samplePartition} is given explicitly) and its per-fold out-of-fold
#' predictions/metrics are fed directly into the SR-TWAS ensemble in place
#' of re-fitting those fine-mapping methods here. Non-fine-mapping methods
#' (lasso, enet, ...) are still cross-validated on the same shared
#' partition.
#'
#' @param data A \code{QtlDataset}, \code{MultiStudyQtlDataset}, or
#'   \code{QtlSumStats}. The \code{MultiStudyQtlDataset} method iterates
#'   the embedded individual-level \code{QtlDataset} entries and the
#'   optional embedded \code{QtlSumStats}, then rbinds the results.
#' @param methods A character vector of short method names, a preset
#'   string (\code{"default"} or \code{"fast_default"}), or a named list
#'   of \code{<method>_weights = args} entries. For QtlSumStats / GwasSumStats
#'   inputs the default switches to the RSS preset
#'   (\code{c("susieRss", "susieInfRss", "lassosumRss", "prsCs", "sdpr")}).
#' @param contexts Optional character vector of contexts to restrict
#'   processing to (QtlDataset / QtlSumStats inputs). Default \code{NULL}
#'   (use all contexts).
#' @param traitId Optional character vector of trait identifiers to
#'   restrict processing to (QtlDataset / QtlSumStats inputs). Default
#'   \code{NULL}.
#' @param region Optional \code{GRanges} for QtlDataset trait selection.
#'   Mutually exclusive with \code{traitId}.
#' @param cisWindow For QtlDataset: cis-window (bp) around each trait's
#'   genomic position when extracting variants. Required when
#'   \code{traitId} is supplied. Mutually exclusive with \code{region}.
#' @param minTwasMaf For QtlDataset: optional minimum minor-allele frequency
#'   applied to the variant set used for TWAS weight learning, on top of the
#'   dataset's construct-time \code{mafCutoff} (the effective cutoff is the
#'   larger of the two). Lets the TWAS pass use a stricter MAF threshold than
#'   fine mapping. \code{NULL} (default) leaves the construct-time cutoff in
#'   place.
#' @param minTwasXvar As \code{minTwasMaf} but for the per-variant genotype
#'   variance cutoff (\code{xvarCutoff}). \code{NULL} (default) leaves the
#'   construct-time cutoff in place.
#' @param jointRegions For QtlDataset with a multi-range \code{region}:
#'   \code{FALSE} (default) learns weights for each range independently and
#'   concatenates them into one entry per (study, context, trait, method);
#'   the per-region fits are kept as a named list and per-region CV is
#'   recorded as a flat \code{cvPerformance} data frame (one row per region).
#'   \code{TRUE} concatenates the ranges' genotypes into one joint fit.
#'   Ignored for a single-range / cis request.
#' @param jointSpecification Optional joint-fit specification (NULL by
#'   default). When NULL, the pipeline runs the implicit multi-trait /
#'   multi-context mr.mash branches as before. When non-NULL, the
#'   argument is parsed and validated via the joint-spec grammar
#'   documented under \code{parseJointSpecification}; the per-spec axis
#'   dispatcher implementation is in progress and a non-NULL value
#'   currently errors with an informative message.
#' @param fineMappingResult Optional \code{\link{FineMappingResult}}.
#'   When supplied, its SuSiE / SuSiE-inf / SuSiE-ash trimmed fits for
#'   the matching (study, context, trait) tuples are injected into
#'   \code{learnTwasWeights} via \code{fittedModels} so SuSiE-family
#'   weight methods reuse the prior fit instead of refitting.
#' @param twasWeights Optional \code{\link{TwasWeights}} resume cache.
#'   For each requested \code{(study, context, trait, method)} tuple
#'   already present in this collection, the cached
#'   \code{TwasWeightsEntry} is copied through and the corresponding
#'   weight fit is skipped. Only the un-cached method subset is fit;
#'   the cached and fresh entries are concatenated in the returned
#'   collection. Per-tuple matching mirrors the \code{fineMappingResult}
#'   cache in \code{\link{fineMappingPipeline}}. Multivariate dispatch
#'   (\code{mvsusie}, \code{mr.mash}) is unaffected because those
#'   methods produce one fit jointly across multiple
#'   \code{(context, trait)} columns.
#' @param cvFolds Integer. Cross-validation folds. Default 5. Set to 0
#'   to skip CV (and ensemble).
#' @param samplePartition Optional pre-defined CV partition data.frame.
#' @param maxCvVariants Maximum number of variants for CV. Default -1
#'   (no limit).
#' @param cvThreads Threads for CV parallelism. Default 1.
#' @param cvWeightMethods Optional override of methods used for CV.
#' @param ensemble Logical. Compute SR-TWAS ensemble weights. Default
#'   \code{TRUE}.
#' @param ensembleR2Threshold Minimum CV R-squared for ensemble
#'   inclusion. Default 0.01.
#' @param ensembleSolver Solver for ensemble stacking. Default
#'   \code{"quadprog"}.
#' @param ensembleAlpha Elastic-net mixing parameter (only when
#'   \code{ensembleSolver = "glmnet"}). Default 1.
#' @param estimatePi If TRUE, estimate spike-and-slab sparsity from
#'   mr.ash before BGLR / qgg spike-and-slab methods that consume it.
#' @param phenotypeCovariatesToResidualize,genotypeCovariatesToResidualize
#'   Character vector (or \code{NULL}) of covariate column names to
#'   residualize against. Forwarded to
#'   \code{\link{getResidualizedPhenotypes}} /
#'   \code{\link{getResidualizedGenotypes}} for \code{QtlDataset} /
#'   \code{MultiStudyQtlDataset} input. Default \code{NULL} (use all
#'   available covariates). Ignored for sumstat inputs.
#' @param residualizePhenotypeCovariates Logical (length 1). When
#'   \code{TRUE} (default) residualize against the phenotype-side
#'   covariates listed in \code{phenotypeCovariatesToResidualize}; set
#'   \code{FALSE} to disable.
#' @param residualizeGenotypeCovariates Logical (length 1). When
#'   \code{TRUE} (default) residualize against the genotype-side
#'   covariates listed in \code{genotypeCovariatesToResidualize}; set
#'   \code{FALSE} to disable.
#' @param dataType Optional data-type tag stamped into every
#'   \code{TwasWeightsEntry$dataType} (e.g. \code{"expression"}).
#' @param verbose Verbosity (0 silent, 1 default, 2 includes external
#'   package messages).
#' @param ... Reserved for method-specific arguments.
#'
#' @return A \code{\link{TwasWeights}} collection keyed by
#'   \code{(study, context, trait, method)}. The \code{ldSketch} slot is
#'   \code{NULL} for individual-level fits and equals the input's
#'   \code{ldSketch} for RSS-derived fits.
#' @export
setGeneric("twasWeightsPipeline",
  function(data, ...) standardGeneric("twasWeightsPipeline"))

#' @rdname twasWeightsPipeline
#' @export
setMethod("twasWeightsPipeline", "QtlDataset",
  function(data,
           methods                = "default",
           contexts               = NULL,
           traitId                = NULL,
           region                 = NULL,
           cisWindow              = NULL,
           minTwasMaf             = NULL,
           minTwasXvar            = NULL,
           jointRegions           = FALSE,
           jointSpecification     = NULL,
           fineMappingResult      = NULL,
           twasWeights            = NULL,
           cvFolds                = 5,
           samplePartition        = NULL,
           maxCvVariants          = -1,
           cvThreads              = 1,
           cvWeightMethods        = NULL,
           ensemble               = TRUE,
           ensembleR2Threshold    = 0.01,
           ensembleSolver         = "quadprog",
           ensembleAlpha          = 1,
           estimatePi             = TRUE,
           retainFit              = TRUE,
           retainFitDetail        = c("slim", "full"),
           phenotypeCovariatesToResidualize = NULL,
           genotypeCovariatesToResidualize  = NULL,
           residualizePhenotypeCovariates   = TRUE,
           residualizeGenotypeCovariates    = TRUE,
           dataType               = NULL,
           naAction               = c("drop", "impute"),
           verbose                = 1,
           ...) {
    naAction <- match.arg(naAction)
    retainFitDetail <- match.arg(retainFitDetail)
    # `cisWindow` expands a trait's own coordinates; `region` is literal.
    # Supplying both signals a misunderstanding -> reject.
    if (!is.null(region) && !is.null(cisWindow)) {
      stop("twasWeightsPipeline(QtlDataset): specify either `region` or ",
           "`cisWindow`, not both. `cisWindow` expands each trait's own ",
           "coordinates, whereas `region` is the literal variant window.")
    }
    xRegions <- .makeXRegions(region, jointRegions)
    # TWAS-specific variant filters: tighten the QtlDataset maf/xvar cutoffs for
    # weight learning (distinct from the construct-time / fine-mapping cutoffs).
    # Modifying the local `data` copy elevates them everywhere downstream
    # (runOne, runMultivariate, and the jointSpec dispatcher all extract from it).
    if (!is.null(minTwasMaf))
      data@mafCutoff  <- max(data@mafCutoff,  as.numeric(minTwasMaf))
    if (!is.null(minTwasXvar))
      data@xvarCutoff <- max(data@xvarCutoff, as.numeric(minTwasXvar))
    parsedJointSpec <- parseJointSpecification(jointSpecification, data)
    norm <- .twasNormalizeMethods(methods)
    .twasCheckMethodCapabilities(norm$tokens, "QtlDataset")
    .twasCheckFineMappingMethods(norm$tokens, fineMappingResult, "QtlDataset")

    # Explicit jointSpecification path: run the per-spec axis dispatcher for
    # mr.mash. Other (univariate) methods continue through the existing
    # per-(context, trait) iteration below.
    jointResult <- NULL
    if (length(parsedJointSpec) > 0L) {
      jointResult <- .twasDispatchJointSpecsQtlDataset(
        parsedJointSpec, data, intersect(norm$tokens, "mrmash"),
        contexts, traitId, cisWindow, dataType, verbose, xRegions = xRegions,
        retainFit = retainFit, retainFitDetail = retainFitDetail)
      drop <- intersect(norm$tokens, "mrmash")
      keep <- setdiff(norm$tokens, drop)
      if (length(keep) == 0L) {
        if (is.null(jointResult))
          stop("twasWeightsPipeline(QtlDataset): no joint fits produced. ",
               "Check that the jointSpecification scope intersects the ",
               "available studies / contexts / traits.")
        return(jointResult)
      }
      # Drop joint-eligible tokens from the per-tuple loop
      norm$tokens <- keep
      keepKeys <- which(sub("(_weights|Weights)$", "",
                            names(norm$methodList)) %in% keep)
      norm$methodList <- norm$methodList[keepKeys]
      methods <- norm$methodList
    }

    study <- getStudy(data)
    allCtx <- getContexts(data)
    useCtx <- if (is.null(contexts)) allCtx else {
      bad <- setdiff(contexts, allCtx)
      if (length(bad) > 0L)
        stop("twasWeightsPipeline(QtlDataset): unknown context(s): ",
             paste(bad, collapse = ", "))
      contexts
    }

    # Collect traits to iterate over. When traitId is specified, use it;
    # when region is specified, use the per-context overlap; when neither
    # is supplied, iterate over every trait in every selected context.
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
    allTraits <- unique(unlist(perCtxTraits))
    if (length(allTraits) == 0L) {
      stop("twasWeightsPipeline(QtlDataset): no traits selected.")
    }

    # Multivariate guard: gate on (nTraits, nContexts).
    nCtx <- length(useCtx)
    .twasCheckMultivariateY(norm$tokens, length(allTraits), nCtx)

    # Multivariate if any requested token is multivariate in either the TWAS
    # capability table (mrmash) or the fine-mapping one (mvsusie / fsusie).
    multivariate <- any(vapply(norm$tokens, .twasIsMultivariateToken,
                               logical(1)))

    runOne <- function(ctx, tid) {
      # Resume cache: per-method check against the supplied `twasWeights`
      # collection. Methods present in the cache for (study, ctx, tid)
      # are pulled directly; the remaining methods (if any) are fit via
      # .twasWeightsPipelineMatrix with a subset weightMethods list.
      cachedRows <- list()
      remaining  <- norm$methodList
      for (mName in names(norm$methodList)) {
        shortMethod <- sub("(_weights|Weights)$", "", mName)
        cached <- .twasCacheLookup(twasWeights, study, ctx, tid, shortMethod)
        if (!is.null(cached)) {
          cachedRows[[shortMethod]] <- cached
          remaining[[mName]] <- NULL
        }
      }
      cachedTw <- .twasBuildFromCachedRows(cachedRows, study, ctx, tid)
      if (length(remaining) == 0L) return(cachedTw)

      Y <- .fmResidPheno(
        data, contexts = ctx, traitId = tid,
        phenotypeCovariatesToResidualize = phenotypeCovariatesToResidualize,
        genotypeCovariatesToResidualize  = genotypeCovariatesToResidualize,
        naAction = naAction)

      # Fine-mapping fits for this (study, ctx, trait); a multi-region fit is a
      # per-region list and is selected blockwise inside the loop.
      allFits <- .twasFineMappingFits(fineMappingResult,
                                      study = study, context = ctx, trait = tid)
      # Fine-mapping's own cross-validated predictions (shared fold partition),
      # reused by the ensemble instead of re-fitting the fine-mapping methods.
      fmCv <- .twasCvResultFor(fineMappingResult, study, ctx, tid)
      nBlocks <- length(xRegions)
      perBlockTw <- lapply(seq_len(nBlocks), function(bi) {
        rg <- xRegions[[bi]]
        X <- if (is.null(rg)) {
          .fmResidGeno(
            data, contexts = ctx, traitId = tid, cisWindow = cisWindow,
            phenotypeCovariatesToResidualize = phenotypeCovariatesToResidualize,
            genotypeCovariatesToResidualize  = genotypeCovariatesToResidualize,
            samples = rownames(Y))
        } else {
          .fmResidGeno(
            data, contexts = ctx, region = rg,
            phenotypeCovariatesToResidualize = phenotypeCovariatesToResidualize,
            genotypeCovariatesToResidualize  = genotypeCovariatesToResidualize,
            samples = rownames(Y))
        }
        common <- intersect(rownames(X), rownames(Y))
        if (length(common) < 2L) {
          stop(sprintf(
            "twasWeightsPipeline: too few shared samples between residualized X and Y for (context='%s', trait='%s').",
            ctx, tid))
        }
        .twasWeightsPipelineMatrix(
          X = X[common, , drop = FALSE], y = Y[common, , drop = FALSE],
          study = study, context = ctx, trait = tid,
          fittedModels = .twasFitsForRegion(allFits, bi, nBlocks),
          cvFolds = cvFolds,
          samplePartition = samplePartition,
          fineMappingCv = fmCv,
          weightMethods = remaining,
          maxCvVariants = maxCvVariants,
          cvThreads = cvThreads,
          cvWeightMethods = cvWeightMethods,
          ensemble = ensemble,
          ensembleR2Threshold = ensembleR2Threshold,
          ensembleSolver = ensembleSolver,
          ensembleAlpha = ensembleAlpha,
          estimatePi = estimatePi,
          standardized = FALSE,
          dataType = dataType,
          ldSketch = NULL,
          verbose = verbose)$twasWeights
      })
      # Single block (cis or jointRegions=TRUE) returns unchanged; multiple
      # blocks (jointRegions=FALSE) concatenate per method into one entry.
      freshTw <- .twasMergeRegions(
        perBlockTw, vapply(xRegions, .twasRegionLabel, character(1)))
      if (is.null(cachedTw)) freshTw
      else .rbindTwasWeights(freshTw, cachedTw, ldSketch = NULL)
    }

    runMultivariate <- function(traits) {
      # Joint over selected (contexts, traits): residualize, intersect
      # samples across contexts, drop subjects with any-NA in Y.
      # Sample basis for Y construction (residualized genotypes are
      # region-independent in their sample set): use the cis window when no
      # explicit region is given, otherwise the first range.
      Xlist <- lapply(useCtx, function(ctx) {
        if (is.null(region)) {
          .fmResidGeno(
            data, contexts = ctx, traitId = traits, cisWindow = cisWindow,
            phenotypeCovariatesToResidualize = phenotypeCovariatesToResidualize,
            genotypeCovariatesToResidualize  = genotypeCovariatesToResidualize)
        } else {
          .fmResidGeno(
            data, contexts = ctx, region = region[1L],
            phenotypeCovariatesToResidualize = phenotypeCovariatesToResidualize,
            genotypeCovariatesToResidualize  = genotypeCovariatesToResidualize)
        }
      })
      # Intersect samples across contexts.
      commonSamples <- Reduce(intersect, lapply(Xlist, rownames))
      if (length(commonSamples) < 2L) {
        stop("twasWeightsPipeline(QtlDataset, multivariate): insufficient samples shared across selected contexts.")
      }

      Yres <- .fmResidPheno(
        data, contexts = useCtx, traitId = traits,
        phenotypeCovariatesToResidualize = phenotypeCovariatesToResidualize,
        genotypeCovariatesToResidualize  = genotypeCovariatesToResidualize,
        naAction = naAction)
      if (length(useCtx) == 1L) Yres <- setNames(list(Yres), useCtx)
      # Concatenate per-context residualized phenotypes column-wise,
      # restricting to commonSamples. Column names become
      # "<context>__<trait>".
      Ymats <- list()
      colMeta <- list()
      for (ctx in names(Yres)) {
        rn <- intersect(commonSamples, rownames(Yres[[ctx]]))
        Ym <- Yres[[ctx]][rn, , drop = FALSE]
        # Pad missing rows so columns line up across contexts.
        if (length(rn) < length(commonSamples)) {
          full <- matrix(NA_real_, nrow = length(commonSamples),
                         ncol = ncol(Ym),
                         dimnames = list(commonSamples, colnames(Ym)))
          full[rn, ] <- Ym
          Ym <- full
        } else {
          Ym <- Ym[commonSamples, , drop = FALSE]
        }
        colnames(Ym) <- paste(ctx, colnames(Ym), sep = "__")
        Ymats[[ctx]] <- Ym
        colMeta[[ctx]] <- data.frame(
          context = ctx, trait = colnames(Yres[[ctx]]),
          stringsAsFactors = FALSE)
      }
      Y <- do.call(cbind, Ymats)
      meta <- do.call(rbind, colMeta)
      # Drop subjects with any NA across Y columns (joint over contexts).
      keep <- complete.cases(Y)
      if (sum(keep) < 2L) {
        stop("twasWeightsPipeline(QtlDataset, multivariate): too few subjects with complete Y across selected (context, trait) columns.")
      }
      Y <- Y[keep, , drop = FALSE]

      # Per-region fit + merge. The cis block reuses the already-extracted
      # genotypes; an explicit region re-extracts that window (genotype
      # residualization is context-independent, so one context suffices).
      # mvsusie/mr.mash joint fits are stored once per (context, trait) row in
      # the FineMappingResult; pull via the first (context, trait) of the group
      # and (for a multi-region fit) select the per-region element.
      perBlockTw <- lapply(seq_along(xRegions), function(bi) {
        rg <- xRegions[[bi]]
        Xr <- if (is.null(rg)) {
          Xlist[[1L]]
        } else {
          .fmResidGeno(
            data, contexts = useCtx[[1L]], region = rg,
            phenotypeCovariatesToResidualize = phenotypeCovariatesToResidualize,
            genotypeCovariatesToResidualize  = genotypeCovariatesToResidualize)
        }
        Xr <- Xr[rownames(Y), , drop = FALSE]
        jointFits <- .twasFitsForRegion(
          .twasFineMappingFits(fineMappingResult, study = study,
                               context = meta$context[[1L]],
                               trait   = meta$trait[[1L]]),
          bi, length(xRegions))
        fmCv <- .twasCvResultFor(fineMappingResult, study,
                                 meta$context[[1L]], meta$trait[[1L]])
        .twasWeightsPipelineMatrix(
          X = Xr, y = Y,
          study   = study,
          context = meta$context,
          trait   = meta$trait,
          # Retain the mr.mash fit parts ({dataDrivenPriorMatrices, w0, V}) on
          # the entry's `fits` slot so fineMappingPipeline can rebuild the
          # mvSuSiE reweighted prior + residual variance from this shared fit.
          # `retainFitDetail` selects the slim payload (default) or the full
          # mr.mash fit.
          retainFits = TRUE,
          retainFitDetail = retainFitDetail,
          fittedModels = jointFits,
          cvFolds = cvFolds,
          samplePartition = samplePartition,
          fineMappingCv = fmCv,
          weightMethods = norm$methodList,
          maxCvVariants = maxCvVariants,
          cvThreads = cvThreads,
          cvWeightMethods = cvWeightMethods,
          ensemble = ensemble,
          ensembleR2Threshold = ensembleR2Threshold,
          ensembleSolver = ensembleSolver,
          ensembleAlpha = ensembleAlpha,
          estimatePi = estimatePi,
          standardized = FALSE,
          dataType = dataType,
          ldSketch = NULL,
          verbose = verbose)$twasWeights
      })
      .twasMergeRegions(
        perBlockTw, vapply(xRegions, .twasRegionLabel, character(1)))
    }

    # Top-level dispatch within the QtlDataset method body.
    if (multivariate) {
      # mvsusie / mr.mash: joint fit. If both nCtx == 1 and nTraits == 1
      # we already rejected above via .twasCheckMultivariateY.
      tw <- runMultivariate(allTraits)
    } else {
      # Univariate methods: sequential over (context, trait).
      out <- NULL
      for (ctx in useCtx) {
        for (tid in perCtxTraits[[ctx]]) {
          twi <- runOne(ctx, tid)
          out <- if (is.null(out)) twi else .rbindTwasWeights(out, twi, ldSketch = NULL)
        }
      }
      tw <- out
    }
    if (is.null(tw) && is.null(jointResult)) {
      stop("twasWeightsPipeline(QtlDataset): no (context, trait) pair produced any weights.")
    }
    if (is.null(tw)) return(jointResult)
    if (is.null(jointResult)) return(tw)
    .rbindTwasWeights(tw, jointResult, ldSketch = NULL)
  })

#' @rdname twasWeightsPipeline
#' @export
setMethod("twasWeightsPipeline", "QtlSumStats",
  function(data,
           methods            = NULL,
           contexts           = NULL,
           traitId            = NULL,
           jointSpecification = NULL,
           fineMappingResult  = NULL,
           twasWeights        = NULL,
           retainFit          = TRUE,
           retainFitDetail    = c("slim", "full"),
           dataType           = NULL,
           verbose            = 1L,
           ...) {
    retainFitDetail <- match.arg(retainFitDetail)
    # summaryStatsQc() is mandatory before twasWeightsPipeline for SumStats
    # input; it also drops variants not present in the ldSketch, so by the
    # time we reach this method every entry's SNP set is a subset of the
    # ldSketch panel.
    .twasAssertQcd(data)

    parsedJointSpec <- parseJointSpecification(jointSpecification, data)

    # Normalize the methods argument into (tokens, methodArgs). The default
    # set excludes fine-mapping methods (susie / susieInf / susieAsh /
    # mvsusie); those must be requested explicitly together with a
    # FineMappingResult passed via `fineMappingResult`.
    if (is.null(methods)) {
      tokens <- c("lasso", "prsCs", "dpr_gibbs")
      methodArgs <- setNames(rep(list(list()), length(tokens)), tokens)
    } else if (is.character(methods)) {
      tokens <- methods
      methodArgs <- setNames(rep(list(list()), length(tokens)), tokens)
    } else if (is.list(methods)) {
      tokens <- names(methods)
      methodArgs <- methods
    } else {
      stop("`methods` must be NULL, a character vector, or a named list ",
           "of <token> = <args> entries.")
    }
    .twasCheckMethodCapabilities(tokens, "QtlSumStats")
    .twasCheckFineMappingMethods(tokens, fineMappingResult, "QtlSumStats")

    jointResult <- NULL
    if (length(parsedJointSpec) > 0L) {
      jointResult <- .twasDispatchJointSpecsQtlSumStats(
        parsedJointSpec, data, intersect(tokens, "mrmash"),
        contexts, traitId, dataType, verbose,
        retainFit = retainFit, retainFitDetail = retainFitDetail)
      keep <- setdiff(tokens, "mrmash")
      if (length(keep) == 0L) {
        if (is.null(jointResult))
          stop("twasWeightsPipeline(QtlSumStats): no joint fits produced.")
        return(jointResult)
      }
      tokens <- keep
      methodArgs <- methodArgs[keep]
    }

    studyCol   <- as.character(data$study)
    contextCol <- as.character(data$context)
    traitCol   <- as.character(data$trait)

    selRows <- seq_len(nrow(data))
    if (!is.null(contexts)) selRows <- selRows[contextCol[selRows] %in% contexts]
    if (!is.null(traitId))  selRows <- selRows[traitCol[selRows]   %in% traitId]
    if (length(selRows) == 0L) {
      stop("twasWeightsPipeline(QtlSumStats): no entries matched the ",
           "supplied contexts / traitId filters.")
    }

    # Partition method tokens by univariate vs multivariate dispatch.
    isMv <- vapply(tokens, .twasIsMultivariateToken, logical(1))
    multivariateTokens <- tokens[isMv]
    univariateTokens   <- tokens[!isMv]

    if (length(multivariateTokens) > 0L) {
      groupKey <- paste(studyCol[selRows], traitCol[selRows], sep = "||")
      perGroupNCtx <- vapply(split(contextCol[selRows], groupKey),
                             length, integer(1))
      if (all(perGroupNCtx < 2L)) {
        stop(sprintf(
          "twasWeightsPipeline(QtlSumStats): multivariate method(s) %s require at least two contexts per (study, trait); the supplied collection has only one context per trait.",
          paste(multivariateTokens, collapse = ", ")))
      }
    }

    ldSketch <- getLdSketch(data)

    rowStudy   <- character(0)
    rowContext <- character(0)
    rowTrait   <- character(0)
    rowMethod  <- character(0)
    rowEntries <- list()

    # ---- Univariate dispatch: per (study, context, trait), per method.
    for (i in selRows) {
      st <- studyCol[i]; ctx <- contextCol[i]; tr <- traitCol[i]

      # Resume cache: pull cached entries up front and reduce the
      # per-entry fit work to the un-cached tokens. When every requested
      # token hits the cache the expensive Z/N/varY/ldMat setup is
      # skipped entirely.
      toFitTokens <- character(0)
      for (tk in univariateTokens) {
        cached <- .twasCacheLookup(twasWeights, st, ctx, tr, tk)
        if (!is.null(cached)) {
          rowStudy   <- c(rowStudy,   st)
          rowContext <- c(rowContext, ctx)
          rowTrait   <- c(rowTrait,   tr)
          rowMethod  <- c(rowMethod,  tk)
          rowEntries[[length(rowEntries) + 1L]] <- cached
        } else {
          toFitTokens <- c(toFitTokens, tk)
        }
      }
      if (length(toFitTokens) == 0L) next

      df <- getSumstatDf(data, study = st, context = ctx, trait = tr,
                          require = c("Z", "N"), derive = "zFromBetaSe")
      variantIds <- df$variant_id
      n <- stats::median(df$N, na.rm = TRUE)
      varY <- getVarY(data, study = st, context = ctx, trait = tr)
      if (is.null(varY)) varY <- 1
      stat <- list(z = df$z, n = n, varY = varY,
                   variantNames = variantIds)
      ldMat <- .twasLdFromSketch(ldSketch, variantIds)

      for (tk in toFitTokens) {
        adapter <- .twasFineMappingMethodAdapters[[tk]]
        fn <- if (!is.null(adapter)) adapter$rssWeightFn
              else .twasMethodCapabilities[[tk]]$sumstatImpl
        userArgs <- methodArgs[[tk]]
        if (is.null(userArgs)) userArgs <- list()
        # When the token is a fine-mapping method, pass the precomputed
        # fit into the *Rss weight function via its dedicated *Fit arg
        # (e.g. susieRssFit, susieInfRssFit, susieAshRssFit). The gate
        # above ensures fineMappingResult is non-NULL here.
        if (!is.null(adapter)) {
          fit <- .twasFineMappingFitFor(fineMappingResult,
                                         study = st, context = ctx, trait = tr,
                                         token = tk)
          if (is.null(fit)) {
            warning(sprintf(
              "twasWeightsPipeline: no '%s' fit found in fineMappingResult for (study=%s, context=%s, trait=%s); skipping.",
              tk, st, ctx, tr))
            next
          }
          userArgs[[adapter$rssFitArg]] <- fit
        }
        weights <- tryCatch(
          do.call(get(fn, mode = "function"),
                  c(list(stat = stat, LD = ldMat), userArgs)),
          error = function(e) {
            warning(sprintf(
              "twasWeightsPipeline: method '%s' failed for (study=%s, context=%s, trait=%s): %s",
              tk, st, ctx, tr, conditionMessage(e)))
            NULL
          })
        if (is.null(weights)) next
        fitAttr <- attr(weights, "fit")
        attr(weights, "fit") <- NULL
        rowStudy   <- c(rowStudy,   st)
        rowContext <- c(rowContext, ctx)
        rowTrait   <- c(rowTrait,   tr)
        rowMethod  <- c(rowMethod,  tk)
        rowEntries[[length(rowEntries) + 1L]] <- TwasWeightsEntry(
          variantIds    = variantIds,
          weights       = as.numeric(weights),
          fits          = fitAttr,
          cvPerformance = NULL,        # Q5: no CV on the sumstat path
          standardized  = TRUE,        # Q4: sumstat-derived weights are standardized
          dataType      = dataType)
      }
    }

    # ---- Multivariate dispatch: per (study, trait), all selected contexts.
    if (length(multivariateTokens) > 0L) {
      groupKey <- paste(studyCol[selRows], traitCol[selRows], sep = "||")
      groups   <- split(selRows, groupKey)
      for (gkey in names(groups)) {
        gIdx <- groups[[gkey]]
        if (length(gIdx) < 2L) next
        st <- studyCol[gIdx[[1L]]]
        tr <- traitCol[gIdx[[1L]]]
        ctxNames <- contextCol[gIdx]

        # Build (variants x contexts) Z matrix. All entries in a (study, trait)
        # group must share an identical variant order after summaryStatsQc().
        firstDf <- getSumstatDf(data,
                                 study = st, context = ctxNames[[1L]],
                                 trait = tr,
                                 require = c("Z", "N"),
                                 derive = "zFromBetaSe")
        variantIds <- firstDf$variant_id
        Z <- matrix(NA_real_, nrow = length(variantIds), ncol = length(gIdx),
                    dimnames = list(variantIds, ctxNames))
        nVec <- numeric(length(gIdx))
        for (kk in seq_along(gIdx)) {
          d <- getSumstatDf(data,
                             study = st, context = ctxNames[[kk]],
                             trait = tr,
                             require = c("Z", "N"),
                             derive = "zFromBetaSe")
          if (!identical(d$variant_id, variantIds))
            stop("twasWeightsPipeline(QtlSumStats, multivariate): every ",
                 "entry for (study='", st, "', trait='", tr,
                 "') must share an identical SNP order after ",
                 "summaryStatsQc(). Use the same ldSketch on every entry.")
          Z[, kk] <- d$z
          nVec[kk] <- stats::median(d$N, na.rm = TRUE)
        }
        names(nVec) <- ctxNames
        stat <- list(z = Z, n = nVec, variantNames = variantIds)
        ldMat <- .twasLdFromSketch(ldSketch, variantIds)

        for (tk in multivariateTokens) {
          adapter <- .twasFineMappingMethodAdapters[[tk]]
          fn <- if (!is.null(adapter)) adapter$rssWeightFn
                else .twasMethodCapabilities[[tk]]$sumstatImpl
          userArgs <- methodArgs[[tk]]
          if (is.null(userArgs)) userArgs <- list()
          # mr.mash (no fine-mapping adapter) is the producer of the mvSuSiE
          # data-driven prior: retain its (slim by default) fit so a downstream
          # mvsusie_rss fineMappingPipeline run can rebuild the reweighted prior.
          # Mirrors the individual-level path, which hardcodes retainFits = TRUE.
          # Respect an explicit caller override of either knob.
          if (is.null(adapter) && tk == "mrmash") {
            if (is.null(userArgs$retainFit)) userArgs$retainFit <- TRUE
            if (is.null(userArgs$fitDetail)) userArgs$fitDetail <- retainFitDetail
          }
          # mvsusie is fine-mapping; thread the pre-fit through. mr.mash is
          # not, so this branch only fires for mvsusie.
          if (!is.null(adapter)) {
            fit <- .twasFineMappingFitFor(fineMappingResult,
                                           study = st,
                                           context = ctxNames[[1L]],
                                           trait = tr,
                                           token = tk)
            if (is.null(fit)) {
              warning(sprintf(
                "twasWeightsPipeline: no '%s' fit found in fineMappingResult for (study=%s, trait=%s); skipping.",
                tk, st, tr))
              next
            }
            userArgs[[adapter$rssFitArg]] <- fit
          }
          weights <- tryCatch(
            do.call(get(fn, mode = "function"),
                    c(list(stat = stat, LD = ldMat), userArgs)),
            error = function(e) {
              warning(sprintf(
                "twasWeightsPipeline: multivariate method '%s' failed for (study=%s, trait=%s): %s",
                tk, st, tr, conditionMessage(e)))
              NULL
            })
          if (is.null(weights)) next
          if (!is.matrix(weights)) weights <- as.matrix(weights)
          fitAttr <- attr(weights, "fit")
          attr(weights, "fit") <- NULL
          for (kk in seq_along(ctxNames)) {
            rowStudy   <- c(rowStudy,   st)
            rowContext <- c(rowContext, ctxNames[[kk]])
            rowTrait   <- c(rowTrait,   tr)
            rowMethod  <- c(rowMethod,  tk)
            rowEntries[[length(rowEntries) + 1L]] <- TwasWeightsEntry(
              variantIds    = variantIds,
              weights       = as.numeric(weights[, kk]),
              # Share the underlying joint fit on the first row only;
              # remaining rows reference the same fit by leaving fits NULL.
              fits          = if (kk == 1L) fitAttr else NULL,
              cvPerformance = NULL,
              standardized  = TRUE,
              dataType      = dataType)
          }
        }
      }
    }

    perTupleResult <- if (length(rowEntries) > 0L)
      TwasWeights(
        study    = rowStudy,
        context  = rowContext,
        trait    = rowTrait,
        method   = rowMethod,
        entry    = rowEntries,
        ldSketch = ldSketch)
      else NULL
    if (is.null(jointResult)) {
      if (is.null(perTupleResult))
        stop("twasWeightsPipeline(QtlSumStats): no entries produced weights.")
      return(perTupleResult)
    }
    if (is.null(perTupleResult)) return(jointResult)
    .rbindTwasWeights(perTupleResult, jointResult, ldSketch = ldSketch)
  })


# =============================================================================
# MultiStudyQtlDataset method
# =============================================================================
# Mirrors the fineMappingPipeline(MultiStudyQtlDataset) recursion: iterates
# the embedded individual-level QtlDataset entries, then processes the
# optional embedded QtlSumStats. The result rows from the two phases are
# rbind'd; the joint columns (when populated by either phase) are carried
# through .rbindTwasWeights.

#' @rdname twasWeightsPipeline
#' @export
setMethod("twasWeightsPipeline", "MultiStudyQtlDataset",
  function(data,
           methods            = "default",
           contexts           = NULL,
           traitId            = NULL,
           region             = NULL,
           cisWindow          = NULL,
           jointRegions       = FALSE,
           jointSpecification = NULL,
           fineMappingResult  = NULL,
           twasWeights        = NULL,
           retainFit          = TRUE,
           retainFitDetail    = c("slim", "full"),
           naAction           = c("drop", "impute"),
           verbose            = 1,
           phenotypeCovariatesToResidualize = NULL,
           genotypeCovariatesToResidualize  = NULL,
           residualizePhenotypeCovariates   = TRUE,
           residualizeGenotypeCovariates    = TRUE,
           ...) {
    naAction <- match.arg(naAction)
    retainFitDetail <- match.arg(retainFitDetail)
    if (!is.null(region) && !is.null(cisWindow)) {
      stop("twasWeightsPipeline(MultiStudyQtlDataset): specify either ",
           "`region` or `cisWindow`, not both.")
    }
    xRegions <- .makeXRegions(region, jointRegions)
    parsedJointSpec <- parseJointSpecification(jointSpecification, data)

    # Gate fine-mapping methods early so the recursion into the embedded
    # QtlDataset / QtlSumStats components doesn't re-run fine-mapping.
    {
      gateTokens <- if (is.character(methods)) methods
                    else if (is.list(methods))
                      sub("(_weights|Weights)$", "", names(methods))
                    else character(0)
      .twasCheckFineMappingMethods(gateTokens, fineMappingResult,
                                    "MultiStudyQtlDataset")
    }

    jointResult <- NULL
    if (length(parsedJointSpec) > 0L) {
      jointMethods <- character(0)
      if (is.character(methods)) jointMethods <- intersect(methods, "mrmash")
      else if (is.list(methods))
        jointMethods <- intersect(sub("(_weights|Weights)$", "",
                                       names(methods)), "mrmash")
      jointResult <- .twasDispatchJointSpecsMultiStudy(
        parsedJointSpec, data, jointMethods,
        contexts, traitId, cisWindow, NULL, verbose, xRegions = xRegions,
        retainFit = retainFit, retainFitDetail = retainFitDetail)
      # Strip mrmash from the methods passed to the per-component recursion.
      if (is.character(methods)) methods <- setdiff(methods, "mrmash")
      else if (is.list(methods)) {
        keep <- sub("(_weights|Weights)$", "", names(methods)) != "mrmash"
        methods <- methods[keep]
      }
      if ((is.character(methods) && length(methods) == 0L) ||
          (is.list(methods) && length(methods) == 0L)) {
        if (is.null(jointResult))
          stop("twasWeightsPipeline(MultiStudyQtlDataset): no joint fits produced.")
        return(jointResult)
      }
    }

    qtlDatasets <- getQtlDatasets(data)
    sumStats <- getSumStats(data)

    out <- NULL
    embeddedLd <- NULL
    for (qdName in names(qtlDatasets)) {
      qd <- qtlDatasets[[qdName]]
      res <- twasWeightsPipeline(
        data               = qd,
        methods            = methods,
        contexts           = contexts,
        traitId            = traitId,
        region             = region,
        cisWindow          = cisWindow,
        jointRegions       = jointRegions,
        jointSpecification = NULL,
        fineMappingResult  = fineMappingResult,
        twasWeights        = twasWeights,
        naAction           = naAction,
        verbose            = verbose,
        ...)
      out <- if (is.null(out)) res
             else .rbindTwasWeights(out, res, ldSketch = NULL)
    }

    if (!is.null(sumStats)) {
      ssRes <- twasWeightsPipeline(
        data               = sumStats,
        methods            = methods,
        contexts           = contexts,
        traitId            = traitId,
        jointSpecification = NULL,
        fineMappingResult  = fineMappingResult,
        twasWeights        = twasWeights,
        verbose            = verbose,
        ...)
      embeddedLd <- getLdSketch(ssRes)
      out <- if (is.null(out)) ssRes
             else .rbindTwasWeights(out, ssRes, ldSketch = embeddedLd)
    }

    perTupleResult <- if (!is.null(out)) {
      # ldSketch: NULL when all studies are individual-level; the embedded
      # sumStats's ldSketch otherwise.
      TwasWeights(
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
        stop("twasWeightsPipeline(MultiStudyQtlDataset): no entries produced weights.")
      return(perTupleResult)
    }
    if (is.null(perTupleResult)) return(jointResult)
    .rbindTwasWeights(perTupleResult, jointResult, ldSketch = embeddedLd)
  })


#' @rdname twasWeightsPipeline
#' @export
setMethod("twasWeightsPipeline", "ANY",
  function(data, ...) {
    stop("twasWeightsPipeline does not accept inputs of class '",
         class(data)[[1L]], "'. Pass a QtlDataset, MultiStudyQtlDataset, ",
         "or QtlSumStats. (GwasSumStats inputs are not supported; ",
         "GWAS-side per-LD-block weights are produced inside the new ",
         "ctwasPipeline / qtlEnrichmentPipeline.)")
  })

# =============================================================================
# Internal matrix-driven TWAS weights pipeline
# =============================================================================
#
# This is the legacy matrix-based pipeline retained as an internal worker.
# The exported, S4-dispatched `twasWeightsPipeline` defined above extracts
# (X, Y) blocks from QtlDataset / QtlSumStats / GwasSumStats and calls this
# function per (study, context, trait) tuple. It returns a single-tuple
# `TwasWeights` collection (one row per method, plus an optional ensemble
# row) along with auxiliary state used during stacking.
#
# Method restrictions imposed at the dispatch layer:
# - PRS-CS, lassosumRss, sdpr, susieRss, susieInfRss, susieAshRss,
#   mrAshRss, mrmashRss, mvsusieRss: RSS-only (refuse QtlDataset).
# - bglrWeights / qgg methods (bayesA/B/C/L/N/R, bLasso, dpr*): individual
#   level only (refuse QtlSumStats / GwasSumStats).
# - mr.mash / mvsusie: multi-trait / multi-context (same rule family as
#   the fine-mapping mvSuSiE family in the design doc).
#
# @noRd
.twasWeightsPipelineMatrix <- function(X,
                                y,
                                study = "",
                                context = "",
                                trait = "",
                                susieFit = NULL,
                                fittedModels = NULL,
                                cvFolds = 5,
                                samplePartition = NULL,
                                fineMappingCv = NULL,
                                weightMethods = "default",
                                maxCvVariants = -1,
                                cvThreads = 1,
                                cvWeightMethods = NULL,
                                ensemble = TRUE,
                                ensembleR2Threshold = 0.01,
                                ensembleSolver = "quadprog",
                                ensembleAlpha = 1,
                                estimatePi = TRUE,
                                standardized = FALSE,
                                dataType = NULL,
                                ldSketch = NULL,
                                retainFits = FALSE,
                                retainFitDetail = c("slim", "full"),
                                verbose = 1) {
  retainFitDetail <- match.arg(retainFitDetail)
  if (is.character(weightMethods)) {
    weightMethods <- .twasMethodLookup(weightMethods)
  }
  if (is.null(fittedModels)) fittedModels <- list()
  if (!is.null(susieFit)) fittedModels[["susie"]] <- susieFit

  # Inject precomputed fine-mapping fits into the per-method args so the
  # corresponding *Weights wrapper extracts coefficients from the fit
  # rather than refitting. The adapter table (.twasFineMappingMethodAdapters)
  # gives the snake_case methodList key and the *Fit argument name for
  # each fine-mapping method.
  for (canonical in names(.twasFineMappingMethodAdapters)) {
    adapter <- .twasFineMappingMethodAdapters[[canonical]]
    if (!is.null(fittedModels[[canonical]]) &&
        !is.null(weightMethods[[adapter$methodKey]]) &&
        is.null(weightMethods[[adapter$methodKey]][[adapter$fitArg]])) {
      weightMethods[[adapter$methodKey]][[adapter$fitArg]] <-
        fittedModels[[canonical]]
    }
  }

  res <- list()
  st <- proc.time()
  if (verbose >= 1) {
    message("Performing TWAS weights computation for univariate analysis methods ...")
    tic()
  }

  if (!is.null(fittedModels[["susie"]]) && !is.null(weightMethods$susie_weights)) {
    res$susieWeightsIntermediate <- .susieWeightIntermediate(fittedModels[["susie"]], X)
  }

  # Check if empirical pi estimation is needed for spike-and-slab methods
  bayesCneedsPi <- "bayes_c_weights" %in% names(weightMethods) &&
    !"pi" %in% names(weightMethods$bayes_c_weights)
  bayesBneedsPi <- "bayes_b_weights" %in% names(weightMethods) &&
    !"probIn" %in% names(weightMethods$bayes_b_weights)
  needsPiEstimation <- (bayesCneedsPi || bayesBneedsPi) && estimatePi

  learnArgs <- list(
    study = study, context = context, trait = trait,
    standardized = standardized, dataType = dataType,
    ldSketch = ldSketch, retainFitDetail = retainFitDetail)

  if (needsPiEstimation) {
    # Run mr.ash first to estimate sparsity
    mrashMethods <- list(mrash_weights = weightMethods[["mrash_weights"]] %||% list())

    if (verbose >= 1) message("  Estimating sparsity from mr.ash ...")
    mrashWeights <- do.call(learnTwasWeights, c(
      list(X = X, Y = y, weightMethods = mrashMethods,
           retainFits = TRUE, verbose = verbose),
      learnArgs))

    empiricalPi <- estimateSparsity(mrashWeights)
    if (verbose >= 1) message(sprintf("  Empirical sparsity estimate: %.4f", empiricalPi))
    res$empiricalPi <- empiricalPi

    # Inject into spike-and-slab methods that need it
    if (bayesCneedsPi) weightMethods$bayes_c_weights$pi <- as.numeric(empiricalPi)
    if (bayesBneedsPi) weightMethods$bayes_b_weights$probIn <- as.numeric(empiricalPi)

    # Run remaining methods (those not already computed)
    remainingFnNames <- setdiff(names(weightMethods), "mrash_weights")

    if (length(remainingFnNames) > 0) {
      remainingMethods <- weightMethods[remainingFnNames]
      remainingTw <- do.call(learnTwasWeights, c(
        list(X = X, Y = y, weightMethods = remainingMethods,
             fittedModels = fittedModels, retainFits = retainFits,
             verbose = verbose),
        learnArgs))
      res$twasWeights <- .rbindTwasWeights(mrashWeights, remainingTw,
                                            ldSketch = ldSketch)
    } else {
      res$twasWeights <- mrashWeights
    }

    # Remove mr.ash if it was not in the original weightMethods
    if (!"mrash_weights" %in% names(weightMethods)) {
      tw <- res$twasWeights
      keep <- as.character(tw$method) != "mrash"
      res$twasWeights <- TwasWeights(
        study   = as.character(tw$study)[keep],
        context = as.character(tw$context)[keep],
        trait   = as.character(tw$trait)[keep],
        method  = as.character(tw$method)[keep],
        entry   = as.list(tw$entry)[keep],
        ldSketch = ldSketch)
    }
  } else {
    # Run all methods at once
    res$twasWeights <- do.call(learnTwasWeights, c(
      list(X = X, Y = y, weightMethods = weightMethods,
           fittedModels = fittedModels, retainFits = retainFits,
           verbose = verbose),
      learnArgs))
  }
  if (verbose >= 1) {
    elapsed <- toc(quiet = TRUE)
    message(sprintf("TWAS weights fitting done in %.1fs", elapsed$toc - elapsed$tic))
  }
  res$twasPredictions <- twasPredict(X, res$twasWeights)

  if (cvFolds > 1) {
    # A few cutting corners to run CV faster at the disadvantage of SuSiE and mr.ash:
    # 1. reset SuSiE to not using refine or adaptive L but to use L from previous analysis
    # 2. at most 100 iterations for mr.ash allowed
    # 3. only use a subset of variants randomly selected to avoid bias
    if (!is.null(fittedModels[["susieInf"]]) && !is.null(weightMethods$susie_inf_weights)) {
      weightMethods$susie_inf_weights$L <- length(fittedModels[["susieInf"]]$V)
      weightMethods$susie_inf_weights$refine <- FALSE
    }
    if (!is.null(weightMethods$susie_weights)) {
      susieCvFit <- fittedModels[["susie"]]
      if (is.null(susieCvFit)) susieCvFit <- fittedModels[["susieInf"]]
      if (!is.null(susieCvFit)) {
        weightMethods$susie_weights$L <- length(susieCvFit$V)
        weightMethods$susie_weights$refine <- FALSE
      }
    }
    if (is.null(cvWeightMethods)) {
      cvWeightMethods <- .filterZeroWeightMethods(weightMethods, res$twasWeights)
    }

    # Fine-mapping handoff: when fineMappingPipeline supplied cross-validated
    # predictions for some methods (shared fold partition + per-fold out-of-
    # fold predictions), reuse them rather than refitting those methods here.
    # Drop them from the CV refit set, adopt the shared partition (unless the
    # caller passed one explicitly), and merge their predictions/metrics into
    # the CV result below so the SR-TWAS ensemble consumes fine-mapping's own
    # cross-validation.
    fmCvPrediction <- NULL; fmCvPerformance <- NULL
    if (!is.null(fineMappingCv) && length(fineMappingCv$prediction) > 0L) {
      if (is.null(samplePartition) && !is.null(fineMappingCv$samplePartition)) {
        samplePartition <- fineMappingCv$samplePartition
      }
      fmBase <- sub("(_predicted|Predicted)$", "", names(fineMappingCv$prediction))
      cvWeightMethods <- cvWeightMethods[
        setdiff(names(cvWeightMethods), paste0(fmBase, "_weights"))]
      yMat <- if (is.matrix(y)) y
              else matrix(y, ncol = 1L, dimnames = list(names(y), NULL))
      sampleNames  <- rownames(X)
      outcomeNames <- colnames(yMat)
      alignFmPred <- function(mat) {
        out <- matrix(NA_real_, length(sampleNames),
                      max(1L, length(outcomeNames)),
                      dimnames = list(sampleNames, outcomeNames))
        rs <- intersect(rownames(mat), sampleNames)
        cs <- if (!is.null(colnames(mat)) && !is.null(outcomeNames))
                intersect(colnames(mat), outcomeNames) else character(0)
        if (length(cs) > 0L) {
          out[rs, cs] <- mat[rs, cs, drop = FALSE]
        } else if (ncol(mat) == ncol(out)) {
          out[rs, ] <- mat[rs, , drop = FALSE]
        }
        out
      }
      fmCvPrediction <- setNames(lapply(fineMappingCv$prediction, alignFmPred),
                                 names(fineMappingCv$prediction))
      fmCvPerformance <- fineMappingCv$performance
    }

    variantsForCv <- c()
    if (maxCvVariants <= 0) {
      maxCvVariants <- Inf
    }
    if (ncol(X) > maxCvVariants) {
      variantsForCv <- sample(colnames(X), maxCvVariants, replace = FALSE)
    }

    if (length(cvWeightMethods) > 0L) {
      if (verbose >= 1) {
        message("Performing cross-validation to assess TWAS weights ...")
        tic()
      }
      res$twasCvResult <- twasWeightsCv(
        X,
        y,
        fold = cvFolds,
        samplePartitions = samplePartition,
        weightMethods = cvWeightMethods,
        maxNumVariants = maxCvVariants,
        numThreads = cvThreads,
        verbose = verbose,
        variantsToKeep = if (length(variantsForCv) > 0) variantsForCv else NULL
      )
      if (verbose >= 1) {
        elapsed <- toc(quiet = TRUE)
        message(sprintf("Cross-validation done in %.1fs", elapsed$toc - elapsed$tic))
      }
    } else {
      # Every CV method came from fine-mapping; no refit needed here.
      res$twasCvResult <- list(samplePartition = samplePartition,
                               prediction = list(), performance = list())
    }

    # Merge fine-mapping's cross-validated predictions/metrics into the CV
    # result so downstream splicing + ensemble treat them as first-class.
    if (!is.null(fmCvPrediction)) {
      res$twasCvResult$prediction  <- c(res$twasCvResult$prediction,
                                         fmCvPrediction)
      res$twasCvResult$performance <- c(res$twasCvResult$performance,
                                        fmCvPerformance)
      if (is.null(res$twasCvResult$samplePartition)) {
        res$twasCvResult$samplePartition <- samplePartition
      }
    }

    # Number of methods participating in cross-validation / ensemble (refit
    # here plus those handed over by fine-mapping).
    nCvMethods <- length(res$twasCvResult$prediction)

    # Splice per-(method, outcome) CV predictions + metrics into the
    # corresponding TwasWeightsEntry$cvPerformance slot.
    res$twasWeights <- .spliceCvIntoTwasWeights(res$twasWeights,
                                                 res$twasCvResult,
                                                 ldSketch = ldSketch)

    # Ensemble learning: learn optimal method combination via stacked regression
    if (isTRUE(ensemble) && nCvMethods <= 1) {
      if (verbose >= 1) message("Ensemble model skipped: only ", nCvMethods,
              " weight method provided (need >= 2 for ensemble learning).")
    }
    if (isTRUE(ensemble) && nCvMethods > 1) {
      if (!is.null(res$twasCvResult$performance)) {
        # Extract R-squared for each method from CV performance table
        methodRsq <- vapply(res$twasCvResult$performance, function(perf) {
          perf[1, "rsq"]
        }, numeric(1))
        names(methodRsq) <- sub("(_performance|Performance)$", "", names(methodRsq))

        # NA R-squared already implies the method is unusable for the ensemble: a
        # method whose CV predictions are degenerate (zero variance across all
        # held-out folds) yields cor(predictions, y) = NA and therefore rsq = NA.
        # So !is.na(methodRsq) is sufficient to drop both NA-rsq and degenerate
        # methods - no separate variance check needed.
        passing <- !is.na(methodRsq) & methodRsq >= ensembleR2Threshold
        nPassing <- sum(passing)

        if (nPassing < 2) {
          # Ensemble (stacked regression) requires at least 2 base learners.
          # Build a per-method status line so the user can see which methods
          # dropped out and why (NA R-squared from degenerate CV predictions,
          # or simply R-squared below the cutoff).
          reason <- ifelse(passing, "(passed)",
                    ifelse(is.na(methodRsq),
                           "(dropped: NA R-squared - likely degenerate CV predictions)",
                           "(dropped: R-squared below cutoff)"))
          passedInfo <- paste0("  ", names(methodRsq), ": R-squared = ",
                               round(methodRsq, 4), " ", reason)
          surviving <- if (nPassing == 1) {
            paste0(" Use the surviving method's weights directly: ",
                   names(methodRsq)[passing], ".")
          } else ""
          if (verbose >= 1) message("Ensemble TWAS skipped: ", nPassing, " of ", length(methodRsq),
                  " methods passed the R-squared cutoff of ", ensembleR2Threshold,
                  " (need >= 2).", surviving, "\n",
                  "Method R-squared values:\n",
                  paste(passedInfo, collapse = "\n"))
        } else {
          passingBase <- names(methodRsq)[passing]

          # Subset cvResults predictions to passing methods, matching on the
          # base name regardless of whether the prediction key uses snake
          # ("lasso_predicted") or camel ("lassoPredicted") form.
          filteredCv <- res$twasCvResult
          predBaseNames <- sub("(_predicted|Predicted)$", "", names(filteredCv$prediction))
          filteredCv$prediction <- filteredCv$prediction[match(passingBase, predBaseNames)]

          # Subset twas_weights to passing methods.
          # Method names on a TwasWeights collection are stored as bare
          # tokens (e.g. "lasso") in the `method` column; the ensemble
          # helper wants snake_case "<method>_weights" keys.
          tw <- res$twasWeights
          twMethodNames <- as.character(tw$method)
          filteredWeights <- setNames(
            lapply(passingBase, function(bn) {
              idx <- which(twMethodNames == bn)
              if (length(idx) == 0L) return(NULL)
              w <- getWeights(tw$entry[[idx[[1L]]]])
              if (!is.matrix(w)) w <- matrix(w, ncol = 1)
              w
            }),
            paste0(passingBase, "_weights"))
          filteredWeights <- Filter(Negate(is.null), filteredWeights)

          if (verbose >= 1) {
            message("Computing ensemble TWAS weights via stacked regression ",
                    "using ", nPassing, " methods: ",
                    paste(passingBase, collapse = ", "), " ...")
            tic()
          }
          ensResult <- ensembleWeights(
            cvResults = filteredCv,
            Y = y,
            twasWeightList = filteredWeights,
            solver = ensembleSolver,
            alpha = ensembleAlpha
          )
          if (verbose >= 1) {
            elapsed <- toc(quiet = TRUE)
            message(sprintf("Ensemble learning done in %.1fs", elapsed$toc - elapsed$tic))
          }

          # Add ensemble weights alongside individual method weights as a
          # new row in the TwasWeights collection.
          if (!is.null(ensResult$ensembleTwasWeights)) {
            ensWt <- ensResult$ensembleTwasWeights
            if (!is.matrix(ensWt)) ensWt <- matrix(ensWt, ncol = 1)
            tw <- res$twasWeights
            # Use the first existing row's (study, context, trait) as the
            # identity tuple for the ensemble row.
            existingStudy   <- as.character(tw$study)[1L]
            existingContext <- as.character(tw$context)[1L]
            existingTrait   <- as.character(tw$trait)[1L]
            existingStd     <- getStandardized(tw$entry[[1L]])
            ensWtVec <- if (ncol(ensWt) == 1L) drop(ensWt) else ensWt
            ensVarIds <- if (!is.null(rownames(ensWt))) rownames(ensWt)
                         else colnames(X)
            ensEntry <- TwasWeightsEntry(
              variantIds   = ensVarIds,
              weights      = ensWtVec,
              cvPerformance = list(
                methodCoef        = ensResult$methodCoef,
                methodPerformance = ensResult$methodPerformance),
              standardized = existingStd)
            ensRow <- TwasWeights(
              study   = existingStudy,
              context = existingContext,
              trait   = existingTrait,
              method  = "ensemble",
              entry   = list(ensEntry),
              ldSketch = ldSketch)
            res$twasWeights <- .rbindTwasWeights(tw, ensRow, ldSketch = ldSketch)
            res$twasPredictions$ensemble_predicted <- X %*% ensWt
          }
          res$ensemble <- ensResult
        }
      }
    }
  }
  res$totalTimeElapsed <- proc.time() - st

  return(res)
}

# Solve ensemble stacking via quadprog (constrained QP with sum-to-1 and non-negativity).
# @param Pvalid Matrix of CV predictions for valid methods (n x Kvalid).
# @param yObs Observed outcome vector (n).
# @param Kvalid Number of valid methods.
# @return Normalized coefficient vector of length Kvalid.
# @noRd
.solveEnsembleQuadprog <- function(Pvalid, yObs, Kvalid) {
  if (!requireNamespace("quadprog", quietly = TRUE)) {
    stop("Package 'quadprog' is required for solver='quadprog'. ",
         "Install with: install.packages('quadprog')")
  }

  Dmat <- crossprod(Pvalid)
  dvec <- as.vector(crossprod(Pvalid, yObs))
  # Ridge term for numerical stability (small relative to trace)
  Dmat <- Dmat + 1e-8 * mean(diag(Dmat)) * diag(Kvalid)

  # Constraint matrix: first constraint is equality (sum = 1), then Kvalid
  # non-negativity constraints.
  Amat <- cbind(rep(1, Kvalid), diag(Kvalid))
  bvec <- c(1, rep(0, Kvalid))

  qpSol <- tryCatch(
    solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1),
    error = function(e) {
      warning("QP solver failed: ", conditionMessage(e),
              ". Falling back to equal weights among valid methods.")
      NULL
    }
  )

  if (is.null(qpSol)) {
    return(rep(1 / Kvalid, Kvalid))
  }

  # Numerical cleanup: clamp to non-negative and renormalize
  zetaValid <- pmax(qpSol$solution, 0)
  zetaSum <- sum(zetaValid)
  if (zetaSum <= 0) {
    warning("QP returned all-zero solution. Falling back to equal weights.")
    return(rep(1 / Kvalid, Kvalid))
  }
  zetaValid / zetaSum
}

# Solve ensemble stacking via NNLS (non-negative least squares, then normalize).
# This is the approach used by SuperLearner (Lawson-Hanson algorithm).
# @param Pvalid Matrix of CV predictions for valid methods (n x Kvalid).
# @param yObs Observed outcome vector (n).
# @param Kvalid Number of valid methods.
# @return Normalized coefficient vector of length Kvalid.
# @noRd
.solveEnsembleNnls <- function(Pvalid, yObs, Kvalid) {
  if (!requireNamespace("nnls", quietly = TRUE)) {
    stop("Package 'nnls' is required for solver='nnls'. ",
         "Install with: install.packages('nnls')")
  }

  fit <- tryCatch(
    nnls::nnls(Pvalid, yObs),
    error = function(e) {
      warning("NNLS solver failed: ", conditionMessage(e),
              ". Falling back to equal weights.")
      NULL
    }
  )

  if (is.null(fit)) {
    return(rep(1 / Kvalid, Kvalid))
  }

  zetaValid <- fit$x
  zetaSum <- sum(zetaValid)
  if (zetaSum <= 0) {
    warning("NNLS returned all-zero solution. Falling back to equal weights.")
    return(rep(1 / Kvalid, Kvalid))
  }
  zetaValid / zetaSum
}

# Solve ensemble stacking via L-BFGS-B (box-constrained optimization, then normalize).
# Uses base R optim() with analytical gradient. No extra dependencies.
# @param Pvalid Matrix of CV predictions for valid methods (n x Kvalid).
# @param yObs Observed outcome vector (n).
# @param Kvalid Number of valid methods.
# @return Normalized coefficient vector of length Kvalid.
# @noRd
.solveEnsembleLbfgsb <- function(Pvalid, yObs, Kvalid) {
  PtP <- crossprod(Pvalid)
  Pty <- as.vector(crossprod(Pvalid, yObs))

  fn <- function(z) sum((yObs - Pvalid %*% z)^2)
  gr <- function(z) as.vector(2 * (PtP %*% z - Pty))

  fit <- tryCatch(
    optim(
      par = rep(1 / Kvalid, Kvalid),
      fn = fn, gr = gr,
      method = "L-BFGS-B",
      lower = rep(0, Kvalid)
    ),
    error = function(e) {
      warning("L-BFGS-B solver failed: ", conditionMessage(e),
              ". Falling back to equal weights.")
      NULL
    }
  )

  if (is.null(fit)) {
    return(rep(1 / Kvalid, Kvalid))
  }

  zetaValid <- pmax(fit$par, 0)
  zetaSum <- sum(zetaValid)
  if (zetaSum <= 0) {
    warning("L-BFGS-B returned all-zero solution. Falling back to equal weights.")
    return(rep(1 / Kvalid, Kvalid))
  }
  zetaValid / zetaSum
}

# Solve ensemble stacking via glmnet (penalized regression with non-negativity).
# Uses cv.glmnet for automatic lambda selection. The alpha parameter controls
# the elastic net mixing: alpha=1 is lasso (sparse), alpha=0 is ridge.
# @param Pvalid Matrix of CV predictions for valid methods (n x Kvalid).
# @param yObs Observed outcome vector (n).
# @param Kvalid Number of valid methods.
# @param alpha Elastic net mixing parameter (default 1 = lasso).
# @return Normalized coefficient vector of length Kvalid.
# @noRd
.solveEnsembleGlmnet <- function(Pvalid, yObs, Kvalid, alpha = 1) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for solver='glmnet'. ",
         "Install with: install.packages('glmnet')")
  }

  fit <- tryCatch(
    glmnet::cv.glmnet(
      x = Pvalid, y = yObs,
      lower.limits = 0,
      alpha = alpha,
      intercept = FALSE
    ),
    error = function(e) {
      warning("glmnet solver failed: ", conditionMessage(e),
              ". Falling back to equal weights.")
      NULL
    }
  )

  if (is.null(fit)) {
    return(rep(1 / Kvalid, Kvalid))
  }

  zetaValid <- as.numeric(coef(fit, s = "lambda.min"))[-1]  # drop intercept
  zetaValid <- pmax(zetaValid, 0)
  zetaSum <- sum(zetaValid)
  if (zetaSum <= 0) {
    warning("glmnet returned all-zero solution. Falling back to equal weights.")
    return(rep(1 / Kvalid, Kvalid))
  }
  zetaValid / zetaSum
}


#' Ensemble TWAS Weights via Stacked Regression
#'
#' Given cross-validated predictions from multiple TWAS weight methods, learns
#' non-negative combination coefficients (summing to 1) via constrained least
#' squares. Returns ensemble weights and per-method performance metrics.
#'
#' This implements the stacked regression approach of SR-TWAS (Dai et al.,
#' Nature Communications, 2024, \doi{10.1038/s41467-024-50983-w}). The ensemble
#' provides a principled way to combine predictions from many TWAS weight
#' methods without requiring the user to pick one method a priori or pay a
#' multiple-testing penalty for running several.
#'
#' For single-dataset usage, pass one \code{twasWeightsCv()} result directly.
#' For multi-dataset ensemble (e.g., combining cell types or reference panels
#' such as CUMC1 + MIT), pass a list of \code{twasWeightsCv()} results along
#' with a list of observed Y vectors - this learns a single joint set of
#' coefficients.
#'
#' @param cvResults Output of \code{\link{twasWeightsCv}}, with \code{$prediction}
#'   (named list of method -> out-of-fold prediction matrix, keys like
#'   \code{"susie_predicted"}). For multi-dataset: a list of such objects.
#' @param Y Observed outcome vector or matrix (samples x contexts). For
#'   multi-dataset: a list of vectors/matrices, one per dataset.
#' @param twasWeightList Optional named list of weight matrices from
#'   \code{\link{learnTwasWeights}}, with keys like \code{"susie_weights"}. Used to
#'   construct the final combined TWAS weight vector. For multi-dataset: a list
#'   of such lists (the first is used as the weight template).
#' @param contextIndex Integer indicating which column of Y to use when Y is a
#'   matrix. Default is 1 (univariate).
#' @param solver Character string specifying the optimization backend.
#'   One of \code{"quadprog"} (default), \code{"nnls"}, \code{"lbfgsb"}, or
#'   \code{"glmnet"}.
#'   \code{"quadprog"} solves a constrained QP with sum-to-1 and non-negativity
#'   constraints. \code{"nnls"} uses non-negative least squares (Lawson-Hanson
#'   algorithm, as in SuperLearner) and normalizes post-hoc. \code{"lbfgsb"}
#'   uses \code{optim(method = "L-BFGS-B")} with non-negativity bounds and
#'   normalizes post-hoc. \code{"glmnet"} uses \code{cv.glmnet} with
#'   \code{lower.limits = 0} for penalized non-negative regression, providing
#'   automatic method selection via regularization. All solvers fall back to
#'   equal weights on failure.
#' @param alpha Elastic net mixing parameter, used only when
#'   \code{solver = "glmnet"}. \code{alpha = 1} (default) is lasso (sparse
#'   method selection), \code{alpha = 0} is ridge, and intermediate values
#'   give elastic net.
#'
#' @return A list with components:
#' \describe{
#'   \item{methodCoef}{Named numeric vector of combination coefficients
#'     (\eqn{\zeta_k}), non-negative and summing to 1. Names are method
#'     base names (e.g., \code{"susie"}, \code{"enet"}).}
#'   \item{ensembleTwasWeights}{Final combined weight vector
#'     \eqn{w = \sum_k \zeta_k w_k}, or NULL if \code{twasWeightList}
#'     is not provided. Returned as a vector for univariate Y, matrix otherwise.}
#'   \item{methodPerformance}{Named numeric vector of per-method R-squared
#'     computed from out-of-fold CV predictions. Preserved so users can still
#'     report individual method performance.}
#' }
#'
#' @details
#' The stacked regression solves:
#' \deqn{\min_{\zeta} \|y - P\zeta\|^2 \quad \text{s.t.} \quad \zeta_k \geq 0,\ \sum_k \zeta_k = 1}
#' where P is the \eqn{n \times K} matrix of out-of-fold predictions from K
#' methods. Four solver backends are available: \code{"quadprog"} enforces
#' both constraints during optimization; \code{"nnls"}, \code{"lbfgsb"}, and
#' \code{"glmnet"} enforce non-negativity only, then normalize coefficients
#' to sum to 1. The \code{"glmnet"} solver additionally applies
#' regularization, which can produce sparse solutions (method selection).
#' If any solver fails, the function falls back to equal weights with a
#' warning.
#'
#' Methods whose CV predictions have zero variance (e.g., when all weights are
#' zero) are excluded from the optimization and assigned \eqn{\zeta_k = 0}.
#'
#' Predictions and Y are aligned by sample names (rownames) when available,
#' rather than assuming positional order.
#'
#' @seealso \code{\link{twasWeightsCv}}, \code{\link{learnTwasWeights}},
#'   \code{\link{twasWeightsPipeline}}
#'
#' @examples
#' \dontrun{
#' # After running twasWeightsPipeline with CV:
#' res <- twasWeightsPipeline(X, y, cvFolds = 5, weightMethods = methods)
#'
#' ens <- ensembleWeights(
#'   cvResults = res$twasCvResult,
#'   Y = y,
#'   twasWeightList = res$twasWeights
#' )
#' ens$methodCoef           # combination weights, sum to 1
#'
#' # Multi-dataset ensemble (e.g., CUMC1 + MIT cell types):
#' ens_multi <- ensembleWeights(
#'   cvResults = list(res_cumc$twasCvResult, res_mit$twasCvResult),
#'   Y = list(y_cumc, y_mit),
#'   twasWeightList = list(res_cumc$twasWeights, res_mit$twasWeights)
#' )
#' }
#'
#' @importFrom stats optim coef complete.cases sd cor
#' @export
ensembleWeights <- function(cvResults, Y, twasWeightList = NULL,
                            contextIndex = 1,
                            solver = c("quadprog", "nnls", "lbfgsb", "glmnet"),
                            alpha = 1) {
  # --- Input validation ---
  solver <- match.arg(solver)
  if (is.null(cvResults)) {
    stop("'cvResults' is required.")
  }
  if (is.null(Y)) {
    stop("'Y' is required.")
  }
  if (!is.numeric(contextIndex) || length(contextIndex) != 1 || contextIndex < 1) {
    stop("'contextIndex' must be a positive integer scalar.")
  }

  # --- Normalize single vs multi-dataset input ---
  # Single dataset: cvResults has $prediction directly (is a twasWeightsCv() output).
  # Multi-dataset: cvResults is a list of such outputs.
  isSingle <- !is.null(cvResults$prediction)
  if (isSingle) {
    cvResults <- list(cvResults)
    Y <- list(Y)
    if (!is.null(twasWeightList)) twasWeightList <- list(twasWeightList)
  } else {
    # Multi-dataset: validate list consistency
    if (!is.list(cvResults) || length(cvResults) == 0) {
      stop("For multi-dataset ensemble, 'cvResults' must be a non-empty list of ",
           "twasWeightsCv() outputs.")
    }
    if (!is.list(Y) || length(Y) != length(cvResults)) {
      stop("'Y' must be a list of the same length as 'cvResults' for ",
           "multi-dataset ensemble.")
    }
    if (!is.null(twasWeightList)) {
      if (!is.list(twasWeightList) || length(twasWeightList) != length(cvResults)) {
        stop("'twasWeightList' must be a list of the same length as 'cvResults'.")
      }
    }
    for (d in seq_along(cvResults)) {
      if (is.null(cvResults[[d]]$prediction)) {
        stop("cvResults[[", d, "]] does not contain '$prediction'. ",
             "Expected a twasWeightsCv() output.")
      }
    }
  }

  # --- Extract and validate method names ---
  predNames <- names(cvResults[[1]]$prediction)
  if (is.null(predNames) || any(predNames == "")) {
    stop("cvResults$prediction must be a named list (output of twasWeightsCv).")
  }
  baseNames <- sub("(_predicted|Predicted)$", "", predNames)
  K <- length(baseNames)

  if (K < 2) {
    stop("Ensemble learning requires at least 2 methods. Found: ", K, ".")
  }

  # Consistency: all datasets must report the same methods in the same order
  for (d in seq_along(cvResults)) {
    if (!identical(names(cvResults[[d]]$prediction), predNames)) {
      stop("All cvResults must have the same method names (in $prediction) ",
           "in the same order. Dataset 1 has: ", paste(predNames, collapse = ", "),
           "; dataset ", d, " has: ",
           paste(names(cvResults[[d]]$prediction), collapse = ", "))
    }
  }

  # --- Build stacked prediction matrix P and observed y vector ---
  predList <- list()
  yList <- list()

  for (d in seq_along(cvResults)) {
    predsD <- cvResults[[d]]$prediction
    yRaw <- Y[[d]]

    # Get sample names from predictions and Y for alignment
    predSamples <- rownames(predsD[[predNames[1]]])
    yNames <- if (is.matrix(yRaw) || is.data.frame(yRaw)) {
      rownames(yRaw)
    } else {
      names(yRaw)
    }

    # Determine sample alignment
    if (!is.null(predSamples) && !is.null(yNames)) {
      common <- intersect(predSamples, yNames)
      if (length(common) == 0) {
        stop("No common sample names between predictions and Y in dataset ", d, ".")
      }
      if (length(common) < length(predSamples) || length(common) < length(yNames)) {
        message("Dataset ", d, ": using ", length(common), " common samples ",
                "(predictions: ", length(predSamples), ", Y: ", length(yNames), ").")
      }
      # Extract y aligned to common samples
      yD <- if (is.matrix(yRaw) || is.data.frame(yRaw)) {
        if (contextIndex > ncol(yRaw)) {
          stop("contextIndex (", contextIndex, ") exceeds number of columns in Y[[",
               d, "]] (", ncol(yRaw), ").")
        }
        as.numeric(as.matrix(yRaw)[match(common, yNames), contextIndex])
      } else {
        as.numeric(yRaw[match(common, yNames)])
      }
      predOrder <- match(common, predSamples)
      nD <- length(common)
    } else {
      # No sample names available: fall back to positional alignment
      yD <- if (is.matrix(yRaw) || is.data.frame(yRaw)) {
        if (contextIndex > ncol(yRaw)) {
          stop("contextIndex (", contextIndex, ") exceeds number of columns in Y[[",
               d, "]] (", ncol(yRaw), ").")
        }
        as.numeric(as.matrix(yRaw)[, contextIndex])
      } else {
        as.numeric(yRaw)
      }
      nD <- length(yD)
      predOrder <- seq_len(nD)
    }

    Pd <- matrix(NA_real_, nrow = nD, ncol = K)
    colnames(Pd) <- baseNames
    for (k in seq_along(predNames)) {
      predMat <- predsD[[predNames[k]]]
      pCol <- if (is.matrix(predMat)) predMat[predOrder, contextIndex] else as.numeric(predMat)[predOrder]
      if (length(pCol) != nD) {
        stop("Prediction length for method '", predNames[k], "' in dataset ", d,
             " (", length(pCol), ") does not match number of aligned samples (", nD, ").")
      }
      Pd[, k] <- pCol
    }
    predList[[d]] <- Pd
    yList[[d]] <- yD
  }

  P <- do.call(rbind, predList)   # (nTotal x K)
  yObs <- unlist(yList)           # (nTotal)

  # Remove rows with any NA (in P or y)
  complete <- complete.cases(P, yObs)
  nDropped <- sum(!complete)
  if (nDropped > 0) {
    message("Dropping ", nDropped, " observation(s) with NA predictions or outcomes.")
  }
  if (sum(complete) < K + 1) {
    stop("Too few complete observations (", sum(complete), ") for ", K,
         " methods. Need at least ", K + 1, ".")
  }
  P <- P[complete, , drop = FALSE]
  yObs <- yObs[complete]

  # --- Identify methods with non-zero variance predictions ---
  methodSds <- apply(P, 2, sd)
  validMethods <- methodSds > .Machine$double.eps
  nValid <- sum(validMethods)

  if (nValid < 1) {
    stop("All methods have zero-variance predictions. Cannot compute ensemble. ",
         "This typically means all methods returned zero weights - check that ",
         "the input data has sufficient signal.")
  }

  # --- Solve for combination coefficients ---
  if (nValid == 1) {
    # Only one method has signal: assign it full weight
    zeta <- rep(0, K)
    zeta[validMethods] <- 1
    names(zeta) <- baseNames
    message("Only one method ('", baseNames[validMethods],
            "') has non-zero variance predictions. Assigning it full weight.")
  } else {
    Pvalid <- P[, validMethods, drop = FALSE]
    Kvalid <- ncol(Pvalid)

    zetaValid <- switch(solver,
      quadprog = .solveEnsembleQuadprog(Pvalid, yObs, Kvalid),
      nnls     = .solveEnsembleNnls(Pvalid, yObs, Kvalid),
      lbfgsb   = .solveEnsembleLbfgsb(Pvalid, yObs, Kvalid),
      glmnet   = .solveEnsembleGlmnet(Pvalid, yObs, Kvalid, alpha = alpha)
    )

    zeta <- rep(0, K)
    zeta[validMethods] <- zetaValid
    names(zeta) <- baseNames
  }

  # --- Performance metrics ---
  methodRsq <- setNames(vapply(seq_len(K), function(k) {
    if (methodSds[k] > 0) cor(yObs, P[, k])^2 else NA_real_
  }, numeric(1)), baseNames)

  # --- Build ensemble TWAS weight vector (uses first dataset's weights) ---
  ensembleTwasWt <- NULL
  if (!is.null(twasWeightList)) {
    wtList <- twasWeightList[[1]]
    if (!is.list(wtList) || length(wtList) == 0) {
      warning("twasWeightList[[1]] is empty or not a list; skipping weight combination.")
    } else {
      wtKeys <- paste0(baseNames, "_weights")
      matched <- wtKeys %in% names(wtList)

      if (any(matched)) {
        firstWt <- wtList[[wtKeys[which(matched)[1]]]]
        if (!is.matrix(firstWt)) firstWt <- matrix(firstWt, ncol = 1)
        p <- nrow(firstWt)
        nContexts <- ncol(firstWt)

        ensembleTwasWt <- matrix(0, nrow = p, ncol = nContexts)
        rownames(ensembleTwasWt) <- rownames(firstWt)
        colnames(ensembleTwasWt) <- colnames(firstWt)

        for (i in which(matched)) {
          wMat <- wtList[[wtKeys[i]]]
          if (!is.matrix(wMat)) wMat <- matrix(wMat, ncol = 1)
          if (!identical(dim(wMat), dim(ensembleTwasWt))) {
            warning("Weight matrix for '", wtKeys[i],
                    "' has inconsistent dimensions; skipping.")
            next
          }
          ensembleTwasWt <- ensembleTwasWt + zeta[i] * wMat
        }

        # For univariate case, return as vector
        if (nContexts == 1) {
          ensembleTwasWt <- setNames(
            as.numeric(ensembleTwasWt),
            rownames(ensembleTwasWt)
          )
        }
      } else {
        warning("No matching weight keys found in twasWeightList. ",
                "Expected keys like: ",
                paste(wtKeys[seq_len(min(3, K))], collapse = ", "))
      }
    }
  }

  list(
    methodCoef = zeta,
    ensembleTwasWeights = ensembleTwasWt,
    methodPerformance = methodRsq
  )
}

