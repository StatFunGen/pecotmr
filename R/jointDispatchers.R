# =============================================================================
# Joint-specification dispatchers for fineMappingPipeline and
# twasWeightsPipeline. The two pipelines share the same axis taxonomy
# (cross-context, cross-trait, cross-study, composed multi-axis) and the
# same per-axis row-enumeration logic; they differ only in (a) which
# joint method tokens they accept, (b) which fitter they call, and (c)
# whether the result is a `QtlFineMappingResult` or a `TwasWeights`
# collection.
#
# The shared bits live at the top of this file as `.buildJoint*` /
# `.enumerate*` helpers. Each pipeline then keeps a per-axis worker that
# wires those helpers to its fit + result-row construction.
# =============================================================================


# =============================================================================
# Shared helpers
# =============================================================================

# Resolve which studies / contexts / traits participate in `spec` given
# `data`. Filters data scope through the spec's `scope` and any explicit
# pipeline-level `contexts` / `traitIds` arguments. Returns a list with
# `studies` (character), `contexts` (named list keyed by study), `traits`
# (named list keyed by study).
# @noRd
.fmResolveSpecScope <- function(spec, data, contexts = NULL,
                                traitIds = NULL) {
  scope <- spec$scope
  studies <- .spListStudies(data)
  if (!is.null(scope$study))
    studies <- intersect(studies, scope$study)

  contextsOut <- list()
  traitsOut <- list()
  for (s in studies) {
    ctxAvail <- .spListContexts(data, s)
    if (!is.null(scope$context))
      ctxAvail <- intersect(ctxAvail, scope$context)
    if (!is.null(contexts)) {
      if (is.list(contexts) && s %in% names(contexts))
        ctxAvail <- intersect(ctxAvail, contexts[[s]])
      else if (is.character(contexts))
        ctxAvail <- intersect(ctxAvail, contexts)
    }
    contextsOut[[s]] <- ctxAvail

    trAvail <- .spListTraits(data, study = s)
    if (!is.null(scope$trait))
      trAvail <- intersect(trAvail, scope$trait)
    if (!is.null(traitIds)) {
      if (is.character(traitIds))
        trAvail <- intersect(trAvail, traitIds)
      else if (is.list(traitIds) && s %in% names(traitIds)) {
        tv <- traitIds[[s]]
        if (is.character(tv)) trAvail <- intersect(trAvail, tv)
      }
    }
    traitsOut[[s]] <- trAvail
  }
  list(studies = studies, contexts = contextsOut, traits = traitsOut)
}


# Build a (variants × tupleRows) Z matrix from a QtlSumStats subset,
# requiring all rows to share an identical SNP order (the post-
# summaryStatsQc contract). Returns list(Z, nVec, variantIds).
# `errorLabel` is woven into the SNP-order error to identify the caller.
# @noRd
.buildJointSumstatZMatrix <- function(data, tupleRows, colLabels, errorLabel) {
  studyCol   <- as.character(data$study)
  contextCol <- as.character(data$context)
  traitCol   <- as.character(data$trait)
  firstDf <- getSumstatDf(data,
                           study   = studyCol[[tupleRows[[1L]]]],
                           context = contextCol[[tupleRows[[1L]]]],
                           trait   = traitCol[[tupleRows[[1L]]]],
                           require = c("SNP", "Z", "N"))
  variantIds <- firstDf$variant_id
  Z <- matrix(NA_real_, nrow = length(variantIds), ncol = length(tupleRows),
              dimnames = list(variantIds, colLabels))
  nVec <- numeric(length(tupleRows))
  for (kk in seq_along(tupleRows)) {
    i <- tupleRows[[kk]]
    d <- getSumstatDf(data,
                       study   = studyCol[[i]],
                       context = contextCol[[i]],
                       trait   = traitCol[[i]],
                       require = c("SNP", "Z", "N"))
    if (!identical(d$variant_id, variantIds))
      stop(sprintf("%s: every entry in a joint group must share an identical SNP order after summaryStatsQc().",
                   errorLabel))
    Z[, kk] <- d$z
    nVec[kk] <- stats::median(d$N, na.rm = TRUE)
  }
  list(Z = Z, nVec = nVec, variantIds = variantIds)
}


# Build a multi-context Y matrix for a single (study, trait) from an
# individual-level QtlDataset. Returns list(X, Y, perTraitContexts) or
# NULL when fewer than 2 contexts carry `tid` or the sample / complete-Y
# subset is too small to fit.
# @noRd
.buildIndividualCrossContextXY <- function(data, tid, scopedContexts,
                                           cisWindow, verbose, label) {
  perTraitContexts <- character(0)
  for (cx in scopedContexts) {
    se <- getPhenotypes(data, contexts = cx)
    if (tid %in% rownames(se))
      perTraitContexts <- c(perTraitContexts, cx)
  }
  if (length(perTraitContexts) < 2L) {
    if (verbose >= 1)
      message(sprintf(
        "%s: trait '%s' present in %d scoped context(s); skipping.",
        label, tid, length(perTraitContexts)))
    return(NULL)
  }
  X <- getResidualizedGenotypes(
    data, contexts = perTraitContexts, traitId = tid,
    cisWindow = cisWindow)
  Yres <- getResidualizedPhenotypes(
    data, contexts = perTraitContexts, traitId = tid)
  commonSamples <- Reduce(intersect,
    c(list(rownames(X)), lapply(Yres, rownames)))
  if (length(commonSamples) < 2L) {
    if (verbose >= 1)
      message(sprintf(
        "%s: trait '%s' has too few shared samples across contexts; skipping.",
        label, tid))
    return(NULL)
  }
  X <- X[commonSamples, , drop = FALSE]
  Y <- do.call(cbind, lapply(perTraitContexts, function(cx) {
    ym <- Yres[[cx]][commonSamples, , drop = FALSE]
    colnames(ym) <- cx
    ym
  }))
  keep <- stats::complete.cases(Y)
  if (sum(keep) < 2L) {
    if (verbose >= 1)
      message(sprintf(
        "%s: trait '%s' has too few complete-Y subjects; skipping.",
        label, tid))
    return(NULL)
  }
  list(X = X[keep, , drop = FALSE], Y = Y[keep, , drop = FALSE],
       perTraitContexts = perTraitContexts)
}


# Build a multi-trait Y matrix for a single (study, context) from an
# individual-level QtlDataset. Returns list(X, Y, traitsHere, se) or NULL
# when fewer than 2 traits live in the context or the sample / complete-Y
# subset is too small.
# @noRd
.buildIndividualCrossTraitXY <- function(data, cx, scopedTraits,
                                         cisWindow, verbose, label, study) {
  se <- getPhenotypes(data, contexts = cx)
  traitsHere <- intersect(scopedTraits, rownames(se))
  if (length(traitsHere) < 2L) {
    if (verbose >= 1)
      message(sprintf(
        "%s: context '%s' (study '%s') has %d scoped trait(s); skipping.",
        label, cx, study, length(traitsHere)))
    return(NULL)
  }
  X <- getResidualizedGenotypes(
    data, contexts = cx, traitId = traitsHere, cisWindow = cisWindow)
  Y <- getResidualizedPhenotypes(
    data, contexts = cx, traitId = traitsHere)
  common <- intersect(rownames(X), rownames(Y))
  if (length(common) < 2L) return(NULL)
  X <- X[common, , drop = FALSE]; Y <- Y[common, , drop = FALSE]
  keep <- stats::complete.cases(Y)
  if (sum(keep) < 2L) return(NULL)
  list(X = X[keep, , drop = FALSE], Y = Y[keep, , drop = FALSE],
       traitsHere = traitsHere, se = se)
}


# Build a composed-axes (context, trait) X/Y for individual-level
# QtlDataset. Returns list(X, Y, tuples) or NULL.
# @noRd
.buildComposedIndividualXY <- function(data, scope, study, cisWindow,
                                       verbose, label) {
  scopedContexts <- scope$contexts[[study]]
  scopedTraits   <- scope$traits[[study]]
  tuples <- list()
  for (cx in scopedContexts) {
    se <- getPhenotypes(data, contexts = cx)
    for (tid in intersect(scopedTraits, rownames(se))) {
      tuples[[length(tuples) + 1L]] <- list(context = cx, trait = tid)
    }
  }
  if (length(tuples) < 2L) {
    if (verbose >= 1)
      message(sprintf(
        "%s: study '%s' has %d (context, trait) tuple(s) in scope; skipping.",
        label, study, length(tuples)))
    return(NULL)
  }
  allContexts <- unique(vapply(tuples, function(t) t$context, character(1L)))
  allTraits   <- unique(vapply(tuples, function(t) t$trait,   character(1L)))
  X <- getResidualizedGenotypes(
    data, contexts = allContexts, traitId = allTraits, cisWindow = cisWindow)
  YresList <- getResidualizedPhenotypes(
    data, contexts = allContexts, traitId = allTraits)
  if (length(allContexts) == 1L) YresList <- setNames(list(YresList), allContexts)
  commonSamples <- Reduce(intersect,
    c(list(rownames(X)), lapply(YresList, rownames)))
  if (length(commonSamples) < 2L) return(NULL)
  X <- X[commonSamples, , drop = FALSE]
  yCols <- list(); colLabels <- character(0)
  for (t in tuples) {
    ym <- YresList[[t$context]]
    if (!(t$trait %in% colnames(ym))) next
    col <- ym[commonSamples, t$trait, drop = FALSE]
    colnames(col) <- paste(t$context, t$trait, sep = ":")
    yCols[[length(yCols) + 1L]] <- col
    colLabels <- c(colLabels, paste(t$context, t$trait, sep = ":"))
  }
  if (length(yCols) < 2L) return(NULL)
  Y <- do.call(cbind, yCols)
  keep <- stats::complete.cases(Y)
  if (sum(keep) < 2L) return(NULL)
  list(X = X[keep, , drop = FALSE], Y = Y[keep, , drop = FALSE],
       tuples = tuples)
}


# Enumerate composed-axes row groups for a QtlSumStats input. Returns the
# list of (rowIdx) per group along with the per-axis identity columns
# needed to label the output row. Groups containing fewer than 2 rows
# are returned unfiltered; the caller decides whether to skip.
# @noRd
.enumerateComposedSumstatGroups <- function(spec, data, scope) {
  axes <- spec$axes
  complement <- setdiff(c("study", "context", "trait"), axes)
  studyCol   <- as.character(data$study)
  contextCol <- as.character(data$context)
  traitCol   <- as.character(data$trait)
  inScope <- vapply(seq_len(nrow(data)), function(i) {
    s <- studyCol[i]; cx <- contextCol[i]; tr <- traitCol[i]
    (s %in% scope$studies) &&
    (cx %in% scope$contexts[[s]]) &&
    (tr %in% scope$traits[[s]])
  }, logical(1L))
  rowIdx <- which(inScope)
  if (length(rowIdx) == 0L) return(NULL)
  groupKey <- if (length(complement) == 0L) {
    rep("__all__", length(rowIdx))
  } else {
    do.call(paste, c(lapply(complement, function(a)
      switch(a, study = studyCol[rowIdx],
                context = contextCol[rowIdx],
                trait = traitCol[rowIdx])),
      sep = "||"))
  }
  groups <- split(rowIdx, groupKey)
  list(groups = groups, axes = axes,
       studyCol = studyCol, contextCol = contextCol, traitCol = traitCol)
}


# =============================================================================
# Fine-mapping dispatchers
# =============================================================================

# Cross-context joint dispatcher for QtlDataset. For each trait in scope
# with >= 2 contexts in scope, fits mvsusieR::mvsusie on the multi-column
# Y matrix and emits ONE result row with context = "joint" and
# jointContexts = "ctx1;ctx2;...".
# @noRd
.fmDispatchCrossContextQtlDataset <- function(spec, data, methods,
                                               contexts, traitIds,
                                               cisWindow,
                                               coverage, secondaryCoverage,
                                               signalCutoff, minAbsCorr,
                                               verbose) {
  jointMethods <- intersect(methods, "mvsusie")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  study <- getStudy(data)
  if (!(study %in% scope$studies)) return(NULL)
  scopedContexts <- scope$contexts[[study]]
  scopedTraits   <- scope$traits[[study]]
  if (length(scopedContexts) < 2L) {
    if (verbose >= 1)
      message(sprintf(
        "jointCrossContext: study '%s' has %d context(s) in scope; skipping cross-context fits.",
        study, length(scopedContexts)))
    return(NULL)
  }

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list(); rowJointContexts <- character(0)

  for (tid in scopedTraits) {
    xy <- .buildIndividualCrossContextXY(
      data, tid, scopedContexts, cisWindow, verbose,
      label = "jointCrossContext")
    if (is.null(xy)) next

    if (verbose >= 1)
      message(sprintf(
        "jointCrossContext: fitting mvsusie for (study='%s', trait='%s') across contexts (%s) ...",
        study, tid, paste(xy$perTraitContexts, collapse = ", ")))
    fit <- fitMvsusie(
      X = xy$X, Y = xy$Y,
      prior_variance = mvsusieR::create_mixture_prior(R = ncol(xy$Y)),
      coverage = coverage)
    fit <- .setFinemappingFitClass(fit, "mvsusie")
    entry <- .fmPostprocessOne(
      fit = fit, method = "mvsusie",
      dataX = xy$X, dataY = NULL,
      coverage = coverage,
      secondaryCoverage = secondaryCoverage,
      signalCutoff = signalCutoff,
      minAbsCorr = minAbsCorr,
      csInput = "X")
    rowStudy   <- c(rowStudy,   study)
    rowContext <- c(rowContext, "joint")
    rowTrait   <- c(rowTrait,   tid)
    rowMethod  <- c(rowMethod,  "mvsusie")
    rowEntries[[length(rowEntries) + 1L]] <- entry
    rowJointContexts <- c(rowJointContexts,
                          paste(xy$perTraitContexts, collapse = ";"))
  }

  if (length(rowStudy) == 0L) return(NULL)
  QtlFineMappingResult(
    study         = rowStudy,
    context       = rowContext,
    trait         = rowTrait,
    method        = rowMethod,
    entry         = rowEntries,
    jointContexts = rowJointContexts,
    ldSketch      = NULL)
}


# Cross-trait joint dispatcher for QtlDataset. Per (study, context), fits
# mvsusieR::mvsusie or fsusieR::susiF (when in `methods`) jointly across
# the scoped traits within that context. Emits ONE result row per
# (study, context, method) with trait = "joint" and jointTraits populated.
# @noRd
.fmDispatchCrossTraitQtlDataset <- function(spec, data, methods,
                                             contexts, traitIds,
                                             cisWindow,
                                             coverage, secondaryCoverage,
                                             signalCutoff, minAbsCorr,
                                             verbose) {
  jointMethods <- intersect(methods, c("mvsusie", "fsusie"))
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  study <- getStudy(data)
  if (!(study %in% scope$studies)) return(NULL)
  scopedContexts <- scope$contexts[[study]]
  scopedTraits   <- scope$traits[[study]]

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list(); rowJointTraits <- character(0)

  for (cx in scopedContexts) {
    xy <- .buildIndividualCrossTraitXY(
      data, cx, scopedTraits, cisWindow, verbose,
      label = "jointCrossTrait", study = study)
    if (is.null(xy)) next

    for (mm in jointMethods) {
      if (verbose >= 1)
        message(sprintf(
          "jointCrossTrait: fitting %s for (study='%s', context='%s') across traits (%s) ...",
          mm, study, cx, paste(xy$traitsHere, collapse = ", ")))
      if (mm == "mvsusie") {
        fit <- fitMvsusie(
          X = xy$X, Y = xy$Y,
          prior_variance = mvsusieR::create_mixture_prior(R = ncol(xy$Y)),
          coverage = coverage)
        fit <- .setFinemappingFitClass(fit, "mvsusie")
        entry <- .fmPostprocessOne(
          fit = fit, method = "mvsusie",
          dataX = xy$X, dataY = NULL,
          coverage = coverage,
          secondaryCoverage = secondaryCoverage,
          signalCutoff = signalCutoff,
          minAbsCorr = minAbsCorr,
          csInput = "X")
      } else {
        rr <- SummarizedExperiment::rowRanges(xy$se)
        ord <- match(colnames(xy$Y), rownames(xy$se))
        rr <- rr[ord]
        pos <- (GenomicRanges::start(rr) + GenomicRanges::end(rr)) / 2
        fit <- fitFsusie(X = xy$X, Y = xy$Y, pos = pos)
        fit <- .setFinemappingFitClass(fit, "fsusie")
        entry <- .fmPostprocessOne(
          fit = fit, method = "fsusie",
          dataX = xy$X, dataY = NULL,
          coverage = coverage,
          secondaryCoverage = secondaryCoverage,
          signalCutoff = signalCutoff,
          minAbsCorr = minAbsCorr,
          csInput = "fsusie")
      }
      rowStudy   <- c(rowStudy,   study)
      rowContext <- c(rowContext, cx)
      rowTrait   <- c(rowTrait,   "joint")
      rowMethod  <- c(rowMethod,  mm)
      rowEntries[[length(rowEntries) + 1L]] <- entry
      rowJointTraits <- c(rowJointTraits,
                          paste(xy$traitsHere, collapse = ";"))
    }
  }

  if (length(rowStudy) == 0L) return(NULL)
  QtlFineMappingResult(
    study       = rowStudy,
    context     = rowContext,
    trait       = rowTrait,
    method      = rowMethod,
    entry       = rowEntries,
    jointTraits = rowJointTraits,
    ldSketch    = NULL)
}


# Cross-context joint dispatcher for QtlSumStats input. Groups the
# selected sumstats rows by (study, trait); each group with >= 2 contexts
# in scope produces one mvsusie_rss fit and one result row with context =
# "joint" and jointContexts populated.
# @noRd
.fmDispatchCrossContextQtlSumStats <- function(spec, data, methods,
                                                contexts, traitIds,
                                                coverage, secondaryCoverage,
                                                signalCutoff, minAbsCorr,
                                                verbose) {
  jointMethods <- intersect(methods, "mvsusie")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  ldSketch <- getLdSketch(data)
  studyCol   <- as.character(data$study)
  contextCol <- as.character(data$context)
  traitCol   <- as.character(data$trait)

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list(); rowJointContexts <- character(0)

  for (s in scope$studies) {
    scopedContexts <- scope$contexts[[s]]
    scopedTraits   <- scope$traits[[s]]
    if (length(scopedContexts) < 2L) {
      if (verbose >= 1)
        message(sprintf(
          "jointCrossContext (QtlSumStats): study '%s' has %d context(s) in scope; skipping.",
          s, length(scopedContexts)))
      next
    }
    for (tid in scopedTraits) {
      tupleRows <- which(studyCol == s & traitCol == tid &
                         contextCol %in% scopedContexts)
      if (length(tupleRows) < 2L) {
        if (verbose >= 1)
          message(sprintf(
            "jointCrossContext (QtlSumStats): (study='%s', trait='%s') has %d scoped context(s); skipping.",
            s, tid, length(tupleRows)))
        next
      }
      ctxNames <- contextCol[tupleRows]
      jz <- .buildJointSumstatZMatrix(
        data, tupleRows, ctxNames,
        errorLabel = "jointCrossContext (QtlSumStats)")
      ldMat <- .fmLdFromSketch(ldSketch, jz$variantIds)
      if (verbose >= 1)
        message(sprintf(
          "jointCrossContext (QtlSumStats): fitting mvsusie_rss for (study='%s', trait='%s', %d contexts) ...",
          s, tid, length(ctxNames)))
      fit <- fitMvsusieRss(
        Z = jz$Z, R = ldMat, N = as.numeric(stats::median(jz$nVec)),
        prior_variance = mvsusieR::create_mixture_prior(R = ncol(jz$Z)),
        coverage = coverage)
      fit <- .setFinemappingFitClass(fit, "mvsusie")
      entry <- .fmPostprocessOne(
        fit = fit, method = "mvsusie",
        dataX = ldMat, dataY = NULL,
        coverage = coverage,
        secondaryCoverage = secondaryCoverage,
        signalCutoff = signalCutoff,
        minAbsCorr = minAbsCorr,
        csInput = "Xcorr")
      rowStudy   <- c(rowStudy,   s)
      rowContext <- c(rowContext, "joint")
      rowTrait   <- c(rowTrait,   tid)
      rowMethod  <- c(rowMethod,  "mvsusie")
      rowEntries[[length(rowEntries) + 1L]] <- entry
      rowJointContexts <- c(rowJointContexts,
                            paste(ctxNames, collapse = ";"))
    }
  }

  if (length(rowStudy) == 0L) return(NULL)
  QtlFineMappingResult(
    study         = rowStudy,
    context       = rowContext,
    trait         = rowTrait,
    method        = rowMethod,
    entry         = rowEntries,
    jointContexts = rowJointContexts,
    ldSketch      = ldSketch)
}


# Cross-trait joint dispatcher for QtlSumStats: groups by (study, context),
# requires >= 2 scoped traits per group. mvsusie_rss only -- no RSS fsusie.
# @noRd
.fmDispatchCrossTraitQtlSumStats <- function(spec, data, methods,
                                              contexts, traitIds,
                                              coverage, secondaryCoverage,
                                              signalCutoff, minAbsCorr,
                                              verbose) {
  if ("fsusie" %in% methods)
    stop("jointCrossTrait (QtlSumStats): fsusie has no RSS variant; ",
         "fsusie cannot participate in sumstats-based joint fits.")
  jointMethods <- intersect(methods, "mvsusie")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  ldSketch <- getLdSketch(data)
  studyCol   <- as.character(data$study)
  contextCol <- as.character(data$context)
  traitCol   <- as.character(data$trait)

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list(); rowJointTraits <- character(0)

  for (s in scope$studies) {
    scopedContexts <- scope$contexts[[s]]
    scopedTraits   <- scope$traits[[s]]
    for (cx in scopedContexts) {
      tupleRows <- which(studyCol == s & contextCol == cx &
                         traitCol %in% scopedTraits)
      if (length(tupleRows) < 2L) {
        if (verbose >= 1)
          message(sprintf(
            "jointCrossTrait (QtlSumStats): (study='%s', context='%s') has %d scoped trait(s); skipping.",
            s, cx, length(tupleRows)))
        next
      }
      trNames <- traitCol[tupleRows]
      jz <- .buildJointSumstatZMatrix(
        data, tupleRows, trNames,
        errorLabel = "jointCrossTrait (QtlSumStats)")
      ldMat <- .fmLdFromSketch(ldSketch, jz$variantIds)
      if (verbose >= 1)
        message(sprintf(
          "jointCrossTrait (QtlSumStats): fitting mvsusie_rss for (study='%s', context='%s', %d traits) ...",
          s, cx, length(trNames)))
      fit <- fitMvsusieRss(
        Z = jz$Z, R = ldMat, N = as.numeric(stats::median(jz$nVec)),
        prior_variance = mvsusieR::create_mixture_prior(R = ncol(jz$Z)),
        coverage = coverage)
      fit <- .setFinemappingFitClass(fit, "mvsusie")
      entry <- .fmPostprocessOne(
        fit = fit, method = "mvsusie",
        dataX = ldMat, dataY = NULL,
        coverage = coverage,
        secondaryCoverage = secondaryCoverage,
        signalCutoff = signalCutoff,
        minAbsCorr = minAbsCorr,
        csInput = "Xcorr")
      rowStudy   <- c(rowStudy,   s)
      rowContext <- c(rowContext, cx)
      rowTrait   <- c(rowTrait,   "joint")
      rowMethod  <- c(rowMethod,  "mvsusie")
      rowEntries[[length(rowEntries) + 1L]] <- entry
      rowJointTraits <- c(rowJointTraits,
                          paste(trNames, collapse = ";"))
    }
  }
  if (length(rowStudy) == 0L) return(NULL)
  QtlFineMappingResult(
    study       = rowStudy,
    context     = rowContext,
    trait       = rowTrait,
    method      = rowMethod,
    entry       = rowEntries,
    jointTraits = rowJointTraits,
    ldSketch    = ldSketch)
}


# Cross-study joint dispatcher for QtlSumStats: groups by (context, trait),
# requires >= 2 scoped studies per group. Sumstats-only by definition;
# individual-level studies are excluded with a message at the caller.
# mvsusie_rss only.
# @noRd
.fmDispatchCrossStudyQtlSumStats <- function(spec, data, methods,
                                              contexts, traitIds,
                                              coverage, secondaryCoverage,
                                              signalCutoff, minAbsCorr,
                                              verbose) {
  if ("fsusie" %in% methods)
    stop("jointCrossStudy: fsusie cannot participate (no RSS variant).")
  jointMethods <- intersect(methods, "mvsusie")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  ldSketch <- getLdSketch(data)
  studyCol   <- as.character(data$study)
  contextCol <- as.character(data$context)
  traitCol   <- as.character(data$trait)

  allCtxs <- unique(unlist(scope$contexts, use.names = FALSE))
  allTrs  <- unique(unlist(scope$traits,   use.names = FALSE))

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list(); rowJointStudies <- character(0)

  for (cx in allCtxs) {
    for (tid in allTrs) {
      tupleRows <- which(contextCol == cx & traitCol == tid &
                         studyCol %in% scope$studies)
      keep <- logical(length(tupleRows))
      for (k in seq_along(tupleRows)) {
        s <- studyCol[tupleRows[k]]
        keep[k] <- (cx %in% scope$contexts[[s]]) &&
                   (tid %in% scope$traits[[s]])
      }
      tupleRows <- tupleRows[keep]
      if (length(tupleRows) < 2L) {
        if (length(tupleRows) > 0L && verbose >= 1)
          message(sprintf(
            "jointCrossStudy: (context='%s', trait='%s') has %d study(ies) in scope; skipping.",
            cx, tid, length(tupleRows)))
        next
      }
      stNames <- studyCol[tupleRows]
      jz <- .buildJointSumstatZMatrix(
        data, tupleRows, stNames,
        errorLabel = "jointCrossStudy")
      ldMat <- .fmLdFromSketch(ldSketch, jz$variantIds)
      if (verbose >= 1)
        message(sprintf(
          "jointCrossStudy: fitting mvsusie_rss for (context='%s', trait='%s', %d studies) ...",
          cx, tid, length(stNames)))
      fit <- fitMvsusieRss(
        Z = jz$Z, R = ldMat, N = as.numeric(stats::median(jz$nVec)),
        prior_variance = mvsusieR::create_mixture_prior(R = ncol(jz$Z)),
        coverage = coverage)
      fit <- .setFinemappingFitClass(fit, "mvsusie")
      entry <- .fmPostprocessOne(
        fit = fit, method = "mvsusie",
        dataX = ldMat, dataY = NULL,
        coverage = coverage,
        secondaryCoverage = secondaryCoverage,
        signalCutoff = signalCutoff,
        minAbsCorr = minAbsCorr,
        csInput = "Xcorr")
      rowStudy   <- c(rowStudy,   "joint")
      rowContext <- c(rowContext, cx)
      rowTrait   <- c(rowTrait,   tid)
      rowMethod  <- c(rowMethod,  "mvsusie")
      rowEntries[[length(rowEntries) + 1L]] <- entry
      rowJointStudies <- c(rowJointStudies,
                           paste(stNames, collapse = ";"))
    }
  }
  if (length(rowStudy) == 0L) return(NULL)
  QtlFineMappingResult(
    study        = rowStudy,
    context      = rowContext,
    trait        = rowTrait,
    method       = rowMethod,
    entry        = rowEntries,
    jointStudies = rowJointStudies,
    ldSketch     = ldSketch)
}


# Composed multi-axis joint dispatcher for QtlDataset. Only axes =
# c("context", "trait") is meaningful for a single-study individual-
# level input. Iterates per (study) (just one), enumerates the
# (context, trait) tuples in scope where the trait exists in the
# context, and fits one mvsusie joint over those tuples.
# @noRd
.fmDispatchComposedQtlDataset <- function(spec, data, methods,
                                           contexts, traitIds, cisWindow,
                                           coverage, secondaryCoverage,
                                           signalCutoff, minAbsCorr,
                                           verbose) {
  axes <- spec$axes
  if ("study" %in% axes)
    stop("composed jointSpecification (QtlDataset): axes including 'study' require sumstats input.")
  if (!setequal(axes, c("context", "trait")))
    stop(sprintf(
      "composed jointSpecification (QtlDataset): unsupported axes (%s) for individual-level input.",
      paste(axes, collapse = ", ")))
  jointMethods <- intersect(methods, "mvsusie")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  study <- getStudy(data)
  if (!(study %in% scope$studies)) return(NULL)
  xy <- .buildComposedIndividualXY(data, scope, study, cisWindow,
                                    verbose,
                                    label = "composed joint (QtlDataset)")
  if (is.null(xy)) return(NULL)

  if (verbose >= 1)
    message(sprintf(
      "composed joint (QtlDataset): fitting mvsusie for study='%s' over %d (context, trait) columns ...",
      study, ncol(xy$Y)))
  fit <- fitMvsusie(
    X = xy$X, Y = xy$Y,
    prior_variance = mvsusieR::create_mixture_prior(R = ncol(xy$Y)),
    coverage = coverage)
  fit <- .setFinemappingFitClass(fit, "mvsusie")
  entry <- .fmPostprocessOne(
    fit = fit, method = "mvsusie",
    dataX = xy$X, dataY = NULL,
    coverage = coverage,
    secondaryCoverage = secondaryCoverage,
    signalCutoff = signalCutoff,
    minAbsCorr = minAbsCorr,
    csInput = "X")
  QtlFineMappingResult(
    study         = study,
    context       = "joint",
    trait         = "joint",
    method        = "mvsusie",
    entry         = list(entry),
    jointContexts = paste(vapply(xy$tuples, function(t) t$context,
                                  character(1)), collapse = ";"),
    jointTraits   = paste(vapply(xy$tuples, function(t) t$trait,
                                  character(1)), collapse = ";"),
    ldSketch      = NULL)
}


# Composed multi-axis joint dispatcher for QtlSumStats. Handles any
# `axes` subset of {study, context, trait} of size >= 2 by iterating the
# complement-axis Cartesian product and emitting one joint fit per
# iteration unit.
# @noRd
.fmDispatchComposedQtlSumStats <- function(spec, data, methods,
                                            contexts, traitIds,
                                            coverage, secondaryCoverage,
                                            signalCutoff, minAbsCorr,
                                            verbose) {
  if ("fsusie" %in% methods)
    stop("composed jointSpecification (QtlSumStats): fsusie has no RSS variant.")
  jointMethods <- intersect(methods, "mvsusie")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  ldSketch <- getLdSketch(data)
  groupInfo <- .enumerateComposedSumstatGroups(spec, data, scope)
  if (is.null(groupInfo)) return(NULL)
  axes <- groupInfo$axes
  studyCol <- groupInfo$studyCol
  contextCol <- groupInfo$contextCol
  traitCol <- groupInfo$traitCol

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list()
  rowJointStudies  <- character(0)
  rowJointContexts <- character(0)
  rowJointTraits   <- character(0)

  for (gIdx in groupInfo$groups) {
    if (length(gIdx) < 2L) {
      if (verbose >= 1)
        message(sprintf(
          "composed joint (QtlSumStats): group has %d row(s); skipping.",
          length(gIdx)))
      next
    }
    colLabels <- vapply(gIdx, function(i)
      paste(studyCol[i], contextCol[i], traitCol[i], sep = ":"),
      character(1L))
    jz <- .buildJointSumstatZMatrix(
      data, gIdx, colLabels,
      errorLabel = "composed joint (QtlSumStats)")
    ldMat <- .fmLdFromSketch(ldSketch, jz$variantIds)
    if (verbose >= 1)
      message(sprintf(
        "composed joint (QtlSumStats): fitting mvsusie_rss for axes=(%s), %d columns ...",
        paste(axes, collapse = ", "), length(gIdx)))
    fit <- fitMvsusieRss(
      Z = jz$Z, R = ldMat, N = as.numeric(stats::median(jz$nVec)),
      prior_variance = mvsusieR::create_mixture_prior(R = ncol(jz$Z)),
      coverage = coverage)
    fit <- .setFinemappingFitClass(fit, "mvsusie")
    entry <- .fmPostprocessOne(
      fit = fit, method = "mvsusie",
      dataX = ldMat, dataY = NULL,
      coverage = coverage,
      secondaryCoverage = secondaryCoverage,
      signalCutoff = signalCutoff,
      minAbsCorr = minAbsCorr,
      csInput = "Xcorr")

    repStudy   <- if ("study"   %in% axes) "joint" else studyCol[gIdx[[1L]]]
    repContext <- if ("context" %in% axes) "joint" else contextCol[gIdx[[1L]]]
    repTrait   <- if ("trait"   %in% axes) "joint" else traitCol[gIdx[[1L]]]
    rowStudy   <- c(rowStudy,   repStudy)
    rowContext <- c(rowContext, repContext)
    rowTrait   <- c(rowTrait,   repTrait)
    rowMethod  <- c(rowMethod,  "mvsusie")
    rowEntries[[length(rowEntries) + 1L]] <- entry
    rowJointStudies <- c(rowJointStudies,
      if ("study" %in% axes) paste(studyCol[gIdx], collapse = ";")
      else NA_character_)
    rowJointContexts <- c(rowJointContexts,
      if ("context" %in% axes) paste(contextCol[gIdx], collapse = ";")
      else NA_character_)
    rowJointTraits <- c(rowJointTraits,
      if ("trait" %in% axes) paste(traitCol[gIdx], collapse = ";")
      else NA_character_)
  }

  if (length(rowStudy) == 0L) return(NULL)
  jsArg <- if (all(is.na(rowJointStudies))) NULL else rowJointStudies
  jcArg <- if (all(is.na(rowJointContexts))) NULL else rowJointContexts
  jtArg <- if (all(is.na(rowJointTraits))) NULL else rowJointTraits
  QtlFineMappingResult(
    study         = rowStudy,
    context       = rowContext,
    trait         = rowTrait,
    method        = rowMethod,
    entry         = rowEntries,
    jointStudies  = jsArg,
    jointContexts = jcArg,
    jointTraits   = jtArg,
    ldSketch      = ldSketch)
}


# Top-level joint dispatcher for fineMappingPipeline(QtlDataset).
# @noRd
.fmDispatchJointSpecsQtlDataset <- function(parsedJointSpec, data,
                                             methods, contexts, traitIds,
                                             cisWindow,
                                             coverage, secondaryCoverage,
                                             signalCutoff, minAbsCorr,
                                             verbose) {
  out <- NULL
  for (i in seq_along(parsedJointSpec)) {
    spec <- parsedJointSpec[[i]]
    axes <- spec$axes
    if (length(axes) > 1L) {
      res <- .fmDispatchComposedQtlDataset(
        spec, data, methods, contexts, traitIds, cisWindow,
        coverage, secondaryCoverage, signalCutoff, minAbsCorr, verbose)
      if (!is.null(res))
        out <- if (is.null(out)) res
               else .rbindFineMappingResult(out, res, ldSketch = NULL)
      next
    }
    axis <- axes[[1L]]
    res <- switch(axis,
      context = .fmDispatchCrossContextQtlDataset(
        spec, data, methods, contexts, traitIds, cisWindow,
        coverage, secondaryCoverage, signalCutoff, minAbsCorr, verbose),
      trait = .fmDispatchCrossTraitQtlDataset(
        spec, data, methods, contexts, traitIds, cisWindow,
        coverage, secondaryCoverage, signalCutoff, minAbsCorr, verbose),
      study = stop(
        "fineMappingPipeline(QtlDataset): jointSpecification with axes = 'study' requires sumstats input. ",
        "QtlDataset represents a single individual-level study; cross-study joints operate on the sumstats slot of MultiStudyQtlDataset or on QtlSumStats directly."),
      stop(sprintf("Unsupported axis: %s", axis)))
    if (!is.null(res))
      out <- if (is.null(out)) res
             else .rbindFineMappingResult(out, res, ldSketch = NULL)
  }
  out
}


# Top-level joint dispatcher for fineMappingPipeline(QtlSumStats).
# @noRd
.fmDispatchJointSpecsQtlSumStats <- function(parsedJointSpec, data,
                                              methods, contexts, traitIds,
                                              coverage, secondaryCoverage,
                                              signalCutoff, minAbsCorr,
                                              verbose) {
  out <- NULL
  for (i in seq_along(parsedJointSpec)) {
    spec <- parsedJointSpec[[i]]
    axes <- spec$axes
    if (length(axes) > 1L) {
      res <- .fmDispatchComposedQtlSumStats(
        spec, data, methods, contexts, traitIds,
        coverage, secondaryCoverage, signalCutoff, minAbsCorr, verbose)
      if (!is.null(res))
        out <- if (is.null(out)) res
               else .rbindFineMappingResult(out, res, ldSketch = getLdSketch(data))
      next
    }
    axis <- axes[[1L]]
    res <- switch(axis,
      context = .fmDispatchCrossContextQtlSumStats(
        spec, data, methods, contexts, traitIds,
        coverage, secondaryCoverage, signalCutoff, minAbsCorr, verbose),
      trait = .fmDispatchCrossTraitQtlSumStats(
        spec, data, methods, contexts, traitIds,
        coverage, secondaryCoverage, signalCutoff, minAbsCorr, verbose),
      study = .fmDispatchCrossStudyQtlSumStats(
        spec, data, methods, contexts, traitIds,
        coverage, secondaryCoverage, signalCutoff, minAbsCorr, verbose),
      stop(sprintf("Unsupported axis: %s", axis)))
    if (!is.null(res))
      out <- if (is.null(out)) res
             else .rbindFineMappingResult(out, res,
                                          ldSketch = getLdSketch(data))
  }
  out
}


# Top-level joint dispatcher for fineMappingPipeline(MultiStudyQtlDataset).
# Routes per-component AND per-axis: a spec with `axes = "study"` only
# touches the sumStats slot; `axes = "context"` and `axes = "trait"` run
# on every component.
# @noRd
.fmDispatchJointSpecsMultiStudy <- function(parsedJointSpec, data,
                                             methods, contexts, traitIds,
                                             cisWindow,
                                             coverage, secondaryCoverage,
                                             signalCutoff, minAbsCorr,
                                             verbose) {
  out <- NULL
  embeddedLd <- NULL
  qtlDatasets <- getQtlDatasets(data)
  sumStats <- getSumStats(data)

  studyAxisSpecs <- parsedJointSpec[vapply(parsedJointSpec,
    function(s) "study" %in% s$axes, logical(1L))]
  nonStudyAxisSpecs <- parsedJointSpec[vapply(parsedJointSpec,
    function(s) !("study" %in% s$axes), logical(1L))]

  if (length(studyAxisSpecs) > 0L && length(qtlDatasets) > 0L && verbose >= 1) {
    message(sprintf(
      "jointCrossStudy: excluding individual-level studies (%s) from cross-study fits (no LD sketch available); sumstats studies participate.",
      paste(names(qtlDatasets), collapse = ", ")))
  }

  if (length(nonStudyAxisSpecs) > 0L) {
    for (qdName in names(qtlDatasets)) {
      qd <- qtlDatasets[[qdName]]
      qdRes <- .fmDispatchJointSpecsQtlDataset(
        nonStudyAxisSpecs, qd, methods, contexts, traitIds, cisWindow,
        coverage, secondaryCoverage, signalCutoff, minAbsCorr, verbose)
      if (!is.null(qdRes))
        out <- if (is.null(out)) qdRes
               else .rbindFineMappingResult(out, qdRes, ldSketch = NULL)
    }
  }

  if (!is.null(sumStats)) {
    ssRes <- .fmDispatchJointSpecsQtlSumStats(
      parsedJointSpec, sumStats, methods, contexts, traitIds,
      coverage, secondaryCoverage, signalCutoff, minAbsCorr, verbose)
    if (!is.null(ssRes)) {
      embeddedLd <- getLdSketch(ssRes)
      out <- if (is.null(out)) ssRes
             else .rbindFineMappingResult(out, ssRes,
                                          ldSketch = embeddedLd)
    }
  } else if (length(studyAxisSpecs) > 0L && verbose >= 1) {
    message("jointCrossStudy: no sumStats slot present on this MultiStudyQtlDataset; cross-study specs produce no result.")
  }
  out
}


# =============================================================================
# TWAS-weights dispatchers
# =============================================================================

# Cross-context joint dispatcher for QtlDataset (twas). Mr.mash across
# scoped contexts per (study, trait).
# @noRd
.twasDispatchCrossContextQtlDataset <- function(spec, data, methods,
                                                 contexts, traitIds,
                                                 cisWindow, dataType,
                                                 verbose) {
  jointMethods <- intersect(methods, "mrmash")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  study <- getStudy(data)
  if (!(study %in% scope$studies)) return(NULL)
  scopedContexts <- scope$contexts[[study]]
  scopedTraits   <- scope$traits[[study]]
  if (length(scopedContexts) < 2L) {
    if (verbose >= 1)
      message(sprintf(
        "jointCrossContext (twas QtlDataset): study '%s' has %d context(s) in scope; skipping.",
        study, length(scopedContexts)))
    return(NULL)
  }

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list(); rowJointContexts <- character(0)

  for (tid in scopedTraits) {
    xy <- .buildIndividualCrossContextXY(
      data, tid, scopedContexts, cisWindow, verbose,
      label = "jointCrossContext (twas QtlDataset)")
    if (is.null(xy)) next

    if (verbose >= 1)
      message(sprintf(
        "jointCrossContext (twas QtlDataset): fitting mr.mash for (study='%s', trait='%s') across contexts (%s) ...",
        study, tid, paste(xy$perTraitContexts, collapse = ", ")))
    weights <- mrmashWeights(X = xy$X, Y = xy$Y)
    if (is.null(rownames(weights))) rownames(weights) <- colnames(xy$X)
    entry <- TwasWeightsEntry(
      variantIds   = rownames(weights),
      weights      = weights,
      standardized = FALSE,
      dataType     = dataType)
    rowStudy   <- c(rowStudy,   study)
    rowContext <- c(rowContext, "joint")
    rowTrait   <- c(rowTrait,   tid)
    rowMethod  <- c(rowMethod,  "mrmash")
    rowEntries[[length(rowEntries) + 1L]] <- entry
    rowJointContexts <- c(rowJointContexts,
                          paste(xy$perTraitContexts, collapse = ";"))
  }

  if (length(rowStudy) == 0L) return(NULL)
  TwasWeights(
    study         = rowStudy,
    context       = rowContext,
    trait         = rowTrait,
    method        = rowMethod,
    entry         = rowEntries,
    jointContexts = rowJointContexts,
    ldSketch      = NULL)
}


# Cross-trait joint dispatcher for QtlDataset (twas). Mr.mash per
# (study, context) across scoped traits.
# @noRd
.twasDispatchCrossTraitQtlDataset <- function(spec, data, methods,
                                               contexts, traitIds,
                                               cisWindow, dataType,
                                               verbose) {
  jointMethods <- intersect(methods, "mrmash")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  study <- getStudy(data)
  if (!(study %in% scope$studies)) return(NULL)
  scopedContexts <- scope$contexts[[study]]
  scopedTraits   <- scope$traits[[study]]

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list(); rowJointTraits <- character(0)

  for (cx in scopedContexts) {
    xy <- .buildIndividualCrossTraitXY(
      data, cx, scopedTraits, cisWindow, verbose,
      label = "jointCrossTrait (twas)", study = study)
    if (is.null(xy)) next

    if (verbose >= 1)
      message(sprintf(
        "jointCrossTrait (twas): fitting mr.mash for (study='%s', context='%s') across traits (%s) ...",
        study, cx, paste(xy$traitsHere, collapse = ", ")))
    weights <- mrmashWeights(X = xy$X, Y = xy$Y)
    if (is.null(rownames(weights))) rownames(weights) <- colnames(xy$X)
    entry <- TwasWeightsEntry(
      variantIds   = rownames(weights),
      weights      = weights,
      standardized = FALSE,
      dataType     = dataType)
    rowStudy   <- c(rowStudy,   study)
    rowContext <- c(rowContext, cx)
    rowTrait   <- c(rowTrait,   "joint")
    rowMethod  <- c(rowMethod,  "mrmash")
    rowEntries[[length(rowEntries) + 1L]] <- entry
    rowJointTraits <- c(rowJointTraits,
                        paste(xy$traitsHere, collapse = ";"))
  }
  if (length(rowStudy) == 0L) return(NULL)
  TwasWeights(
    study       = rowStudy,
    context     = rowContext,
    trait       = rowTrait,
    method      = rowMethod,
    entry       = rowEntries,
    jointTraits = rowJointTraits,
    ldSketch    = NULL)
}


# Cross-context joint dispatcher for QtlSumStats (twas). Mr.mash.rss per
# (study, trait).
# @noRd
.twasDispatchCrossContextQtlSumStats <- function(spec, data, methods,
                                                  contexts, traitIds,
                                                  dataType, verbose) {
  jointMethods <- intersect(methods, "mrmash")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  ldSketch <- getLdSketch(data)
  studyCol   <- as.character(data$study)
  contextCol <- as.character(data$context)
  traitCol   <- as.character(data$trait)

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list(); rowJointContexts <- character(0)

  for (s in scope$studies) {
    scopedContexts <- scope$contexts[[s]]
    scopedTraits   <- scope$traits[[s]]
    if (length(scopedContexts) < 2L) {
      if (verbose >= 1)
        message(sprintf(
          "jointCrossContext (twas QtlSumStats): study '%s' has %d context(s) in scope; skipping.",
          s, length(scopedContexts)))
      next
    }
    for (tid in scopedTraits) {
      tupleRows <- which(studyCol == s & traitCol == tid &
                         contextCol %in% scopedContexts)
      if (length(tupleRows) < 2L) {
        if (verbose >= 1)
          message(sprintf(
            "jointCrossContext (twas QtlSumStats): (study='%s', trait='%s') has %d scoped context(s); skipping.",
            s, tid, length(tupleRows)))
        next
      }
      ctxNames <- contextCol[tupleRows]
      jz <- .buildJointSumstatZMatrix(
        data, tupleRows, ctxNames,
        errorLabel = "jointCrossContext (twas QtlSumStats)")
      ldMat <- .fmLdFromSketch(ldSketch, jz$variantIds)
      stat <- list(z = jz$Z, N = jz$nVec)
      if (verbose >= 1)
        message(sprintf(
          "jointCrossContext (twas QtlSumStats): fitting mr.mash.rss for (study='%s', trait='%s', %d contexts) ...",
          s, tid, length(ctxNames)))
      weights <- mrmashRssWeights(stat = stat, LD = ldMat)
      if (is.null(rownames(weights))) rownames(weights) <- jz$variantIds
      entry <- TwasWeightsEntry(
        variantIds   = rownames(weights),
        weights      = weights,
        standardized = TRUE,
        dataType     = dataType)
      rowStudy   <- c(rowStudy,   s)
      rowContext <- c(rowContext, "joint")
      rowTrait   <- c(rowTrait,   tid)
      rowMethod  <- c(rowMethod,  "mrmash")
      rowEntries[[length(rowEntries) + 1L]] <- entry
      rowJointContexts <- c(rowJointContexts,
                            paste(ctxNames, collapse = ";"))
    }
  }
  if (length(rowStudy) == 0L) return(NULL)
  TwasWeights(
    study         = rowStudy,
    context       = rowContext,
    trait         = rowTrait,
    method        = rowMethod,
    entry         = rowEntries,
    jointContexts = rowJointContexts,
    ldSketch      = ldSketch)
}


# Cross-trait joint dispatcher for QtlSumStats (twas). Mr.mash.rss per
# (study, context).
# @noRd
.twasDispatchCrossTraitQtlSumStats <- function(spec, data, methods,
                                                contexts, traitIds,
                                                dataType, verbose) {
  jointMethods <- intersect(methods, "mrmash")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  ldSketch <- getLdSketch(data)
  studyCol   <- as.character(data$study)
  contextCol <- as.character(data$context)
  traitCol   <- as.character(data$trait)

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list(); rowJointTraits <- character(0)

  for (s in scope$studies) {
    scopedContexts <- scope$contexts[[s]]
    scopedTraits   <- scope$traits[[s]]
    for (cx in scopedContexts) {
      tupleRows <- which(studyCol == s & contextCol == cx &
                         traitCol %in% scopedTraits)
      if (length(tupleRows) < 2L) {
        if (verbose >= 1)
          message(sprintf(
            "jointCrossTrait (twas QtlSumStats): (study='%s', context='%s') has %d scoped trait(s); skipping.",
            s, cx, length(tupleRows)))
        next
      }
      trNames <- traitCol[tupleRows]
      jz <- .buildJointSumstatZMatrix(
        data, tupleRows, trNames,
        errorLabel = "jointCrossTrait (twas QtlSumStats)")
      ldMat <- .fmLdFromSketch(ldSketch, jz$variantIds)
      stat <- list(z = jz$Z, N = jz$nVec)
      if (verbose >= 1)
        message(sprintf(
          "jointCrossTrait (twas QtlSumStats): fitting mr.mash.rss for (study='%s', context='%s', %d traits) ...",
          s, cx, length(trNames)))
      weights <- mrmashRssWeights(stat = stat, LD = ldMat)
      if (is.null(rownames(weights))) rownames(weights) <- jz$variantIds
      entry <- TwasWeightsEntry(
        variantIds   = rownames(weights),
        weights      = weights,
        standardized = TRUE,
        dataType     = dataType)
      rowStudy   <- c(rowStudy,   s)
      rowContext <- c(rowContext, cx)
      rowTrait   <- c(rowTrait,   "joint")
      rowMethod  <- c(rowMethod,  "mrmash")
      rowEntries[[length(rowEntries) + 1L]] <- entry
      rowJointTraits <- c(rowJointTraits,
                          paste(trNames, collapse = ";"))
    }
  }
  if (length(rowStudy) == 0L) return(NULL)
  TwasWeights(
    study       = rowStudy,
    context     = rowContext,
    trait       = rowTrait,
    method      = rowMethod,
    entry       = rowEntries,
    jointTraits = rowJointTraits,
    ldSketch    = ldSketch)
}


# Cross-study joint dispatcher for QtlSumStats (twas). Mr.mash.rss per
# (context, trait).
# @noRd
.twasDispatchCrossStudyQtlSumStats <- function(spec, data, methods,
                                                contexts, traitIds,
                                                dataType, verbose) {
  jointMethods <- intersect(methods, "mrmash")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  ldSketch <- getLdSketch(data)
  studyCol   <- as.character(data$study)
  contextCol <- as.character(data$context)
  traitCol   <- as.character(data$trait)

  allCtxs <- unique(unlist(scope$contexts, use.names = FALSE))
  allTrs  <- unique(unlist(scope$traits,   use.names = FALSE))

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list(); rowJointStudies <- character(0)

  for (cx in allCtxs) {
    for (tid in allTrs) {
      tupleRows <- which(contextCol == cx & traitCol == tid &
                         studyCol %in% scope$studies)
      keep <- logical(length(tupleRows))
      for (k in seq_along(tupleRows)) {
        s <- studyCol[tupleRows[k]]
        keep[k] <- (cx %in% scope$contexts[[s]]) &&
                   (tid %in% scope$traits[[s]])
      }
      tupleRows <- tupleRows[keep]
      if (length(tupleRows) < 2L) {
        if (length(tupleRows) > 0L && verbose >= 1)
          message(sprintf(
            "jointCrossStudy (twas): (context='%s', trait='%s') has %d study(ies) in scope; skipping.",
            cx, tid, length(tupleRows)))
        next
      }
      stNames <- studyCol[tupleRows]
      jz <- .buildJointSumstatZMatrix(
        data, tupleRows, stNames,
        errorLabel = "jointCrossStudy (twas)")
      ldMat <- .fmLdFromSketch(ldSketch, jz$variantIds)
      stat <- list(z = jz$Z, N = jz$nVec)
      if (verbose >= 1)
        message(sprintf(
          "jointCrossStudy (twas): fitting mr.mash.rss for (context='%s', trait='%s', %d studies) ...",
          cx, tid, length(stNames)))
      weights <- mrmashRssWeights(stat = stat, LD = ldMat)
      if (is.null(rownames(weights))) rownames(weights) <- jz$variantIds
      entry <- TwasWeightsEntry(
        variantIds   = rownames(weights),
        weights      = weights,
        standardized = TRUE,
        dataType     = dataType)
      rowStudy   <- c(rowStudy,   "joint")
      rowContext <- c(rowContext, cx)
      rowTrait   <- c(rowTrait,   tid)
      rowMethod  <- c(rowMethod,  "mrmash")
      rowEntries[[length(rowEntries) + 1L]] <- entry
      rowJointStudies <- c(rowJointStudies,
                           paste(stNames, collapse = ";"))
    }
  }
  if (length(rowStudy) == 0L) return(NULL)
  TwasWeights(
    study        = rowStudy,
    context      = rowContext,
    trait        = rowTrait,
    method       = rowMethod,
    entry        = rowEntries,
    jointStudies = rowJointStudies,
    ldSketch     = ldSketch)
}


# Composed multi-axis joint dispatcher for QtlSumStats (twas).
# @noRd
.twasDispatchComposedQtlSumStats <- function(spec, data, methods,
                                              contexts, traitIds,
                                              dataType, verbose) {
  jointMethods <- intersect(methods, "mrmash")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  ldSketch <- getLdSketch(data)
  groupInfo <- .enumerateComposedSumstatGroups(spec, data, scope)
  if (is.null(groupInfo)) return(NULL)
  axes <- groupInfo$axes
  studyCol <- groupInfo$studyCol
  contextCol <- groupInfo$contextCol
  traitCol <- groupInfo$traitCol

  rowStudy <- character(0); rowContext <- character(0)
  rowTrait <- character(0); rowMethod <- character(0)
  rowEntries <- list()
  rowJointStudies  <- character(0)
  rowJointContexts <- character(0)
  rowJointTraits   <- character(0)

  for (gIdx in groupInfo$groups) {
    if (length(gIdx) < 2L) {
      if (verbose >= 1)
        message(sprintf(
          "composed joint (twas QtlSumStats): group has %d row(s); skipping.",
          length(gIdx)))
      next
    }
    colLabels <- vapply(gIdx, function(i)
      paste(studyCol[i], contextCol[i], traitCol[i], sep = ":"),
      character(1L))
    jz <- .buildJointSumstatZMatrix(
      data, gIdx, colLabels,
      errorLabel = "composed joint (twas QtlSumStats)")
    ldMat <- .fmLdFromSketch(ldSketch, jz$variantIds)
    stat <- list(z = jz$Z, N = jz$nVec)
    if (verbose >= 1)
      message(sprintf(
        "composed joint (twas QtlSumStats): fitting mr.mash.rss for axes=(%s), %d columns ...",
        paste(axes, collapse = ", "), length(gIdx)))
    weights <- mrmashRssWeights(stat = stat, LD = ldMat)
    if (is.null(rownames(weights))) rownames(weights) <- jz$variantIds
    entry <- TwasWeightsEntry(
      variantIds   = rownames(weights),
      weights      = weights,
      standardized = TRUE,
      dataType     = dataType)

    repStudy   <- if ("study"   %in% axes) "joint" else studyCol[gIdx[[1L]]]
    repContext <- if ("context" %in% axes) "joint" else contextCol[gIdx[[1L]]]
    repTrait   <- if ("trait"   %in% axes) "joint" else traitCol[gIdx[[1L]]]
    rowStudy   <- c(rowStudy,   repStudy)
    rowContext <- c(rowContext, repContext)
    rowTrait   <- c(rowTrait,   repTrait)
    rowMethod  <- c(rowMethod,  "mrmash")
    rowEntries[[length(rowEntries) + 1L]] <- entry
    rowJointStudies <- c(rowJointStudies,
      if ("study" %in% axes) paste(studyCol[gIdx], collapse = ";")
      else NA_character_)
    rowJointContexts <- c(rowJointContexts,
      if ("context" %in% axes) paste(contextCol[gIdx], collapse = ";")
      else NA_character_)
    rowJointTraits <- c(rowJointTraits,
      if ("trait" %in% axes) paste(traitCol[gIdx], collapse = ";")
      else NA_character_)
  }
  if (length(rowStudy) == 0L) return(NULL)
  jsArg <- if (all(is.na(rowJointStudies))) NULL else rowJointStudies
  jcArg <- if (all(is.na(rowJointContexts))) NULL else rowJointContexts
  jtArg <- if (all(is.na(rowJointTraits))) NULL else rowJointTraits
  TwasWeights(
    study         = rowStudy,
    context       = rowContext,
    trait         = rowTrait,
    method        = rowMethod,
    entry         = rowEntries,
    jointStudies  = jsArg,
    jointContexts = jcArg,
    jointTraits   = jtArg,
    ldSketch      = ldSketch)
}


# Composed multi-axis joint dispatcher for QtlDataset (twas). axes =
# c("context", "trait") only.
# @noRd
.twasDispatchComposedQtlDataset <- function(spec, data, methods,
                                             contexts, traitIds, cisWindow,
                                             dataType, verbose) {
  axes <- spec$axes
  if ("study" %in% axes)
    stop("composed jointSpecification (twas QtlDataset): axes including 'study' require sumstats input.")
  if (!setequal(axes, c("context", "trait")))
    stop(sprintf("composed jointSpecification (twas QtlDataset): unsupported axes (%s) for individual-level input.",
                 paste(axes, collapse = ", ")))
  jointMethods <- intersect(methods, "mrmash")
  if (length(jointMethods) == 0L) return(NULL)

  scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                traitIds = traitIds)
  study <- getStudy(data)
  if (!(study %in% scope$studies)) return(NULL)
  xy <- .buildComposedIndividualXY(data, scope, study, cisWindow,
                                    verbose,
                                    label = "composed joint (twas QtlDataset)")
  if (is.null(xy)) return(NULL)

  if (verbose >= 1)
    message(sprintf(
      "composed joint (twas QtlDataset): fitting mr.mash for study='%s' over %d (context, trait) columns ...",
      study, ncol(xy$Y)))
  weights <- mrmashWeights(X = xy$X, Y = xy$Y)
  if (is.null(rownames(weights))) rownames(weights) <- colnames(xy$X)
  entry <- TwasWeightsEntry(
    variantIds   = rownames(weights),
    weights      = weights,
    standardized = FALSE,
    dataType     = dataType)
  TwasWeights(
    study         = study,
    context       = "joint",
    trait         = "joint",
    method        = "mrmash",
    entry         = list(entry),
    jointContexts = paste(vapply(xy$tuples, function(t) t$context,
                                  character(1)), collapse = ";"),
    jointTraits   = paste(vapply(xy$tuples, function(t) t$trait,
                                  character(1)), collapse = ";"),
    ldSketch      = NULL)
}


# Top-level joint dispatcher for twasWeightsPipeline(QtlDataset).
# @noRd
.twasDispatchJointSpecsQtlDataset <- function(parsedJointSpec, data,
                                               methods, contexts, traitIds,
                                               cisWindow, dataType,
                                               verbose) {
  out <- NULL
  for (i in seq_along(parsedJointSpec)) {
    spec <- parsedJointSpec[[i]]
    axes <- spec$axes
    if (length(axes) > 1L) {
      res <- .twasDispatchComposedQtlDataset(
        spec, data, methods, contexts, traitIds, cisWindow, dataType, verbose)
      if (!is.null(res))
        out <- if (is.null(out)) res
               else .rbindTwasWeights(out, res, ldSketch = NULL)
      next
    }
    axis <- axes[[1L]]
    res <- switch(axis,
      context = .twasDispatchCrossContextQtlDataset(
        spec, data, methods, contexts, traitIds, cisWindow, dataType, verbose),
      trait = .twasDispatchCrossTraitQtlDataset(
        spec, data, methods, contexts, traitIds, cisWindow, dataType, verbose),
      study = stop(
        "twasWeightsPipeline(QtlDataset): jointSpecification with axes = 'study' requires sumstats input."),
      stop(sprintf("Unsupported axis: %s", axis)))
    if (!is.null(res))
      out <- if (is.null(out)) res
             else .rbindTwasWeights(out, res, ldSketch = NULL)
  }
  out
}


# Top-level joint dispatcher for twasWeightsPipeline(QtlSumStats).
# @noRd
.twasDispatchJointSpecsQtlSumStats <- function(parsedJointSpec, data,
                                                methods, contexts, traitIds,
                                                dataType, verbose) {
  out <- NULL
  for (i in seq_along(parsedJointSpec)) {
    spec <- parsedJointSpec[[i]]
    axes <- spec$axes
    if (length(axes) > 1L) {
      res <- .twasDispatchComposedQtlSumStats(
        spec, data, methods, contexts, traitIds, dataType, verbose)
      if (!is.null(res))
        out <- if (is.null(out)) res
               else .rbindTwasWeights(out, res, ldSketch = getLdSketch(data))
      next
    }
    axis <- axes[[1L]]
    res <- switch(axis,
      context = .twasDispatchCrossContextQtlSumStats(
        spec, data, methods, contexts, traitIds, dataType, verbose),
      trait = .twasDispatchCrossTraitQtlSumStats(
        spec, data, methods, contexts, traitIds, dataType, verbose),
      study = .twasDispatchCrossStudyQtlSumStats(
        spec, data, methods, contexts, traitIds, dataType, verbose),
      stop(sprintf("Unsupported axis: %s", axis)))
    if (!is.null(res))
      out <- if (is.null(out)) res
             else .rbindTwasWeights(out, res, ldSketch = getLdSketch(data))
  }
  out
}


# Top-level joint dispatcher for twasWeightsPipeline(MultiStudyQtlDataset).
# @noRd
.twasDispatchJointSpecsMultiStudy <- function(parsedJointSpec, data,
                                               methods, contexts, traitIds,
                                               cisWindow, dataType, verbose) {
  out <- NULL
  embeddedLd <- NULL
  qtlDatasets <- getQtlDatasets(data)
  sumStats <- getSumStats(data)

  studyAxisSpecs <- parsedJointSpec[vapply(parsedJointSpec,
    function(s) "study" %in% s$axes, logical(1L))]
  nonStudyAxisSpecs <- parsedJointSpec[vapply(parsedJointSpec,
    function(s) !("study" %in% s$axes), logical(1L))]

  if (length(studyAxisSpecs) > 0L && length(qtlDatasets) > 0L && verbose >= 1) {
    message(sprintf(
      "jointCrossStudy (twas): excluding individual-level studies (%s) from cross-study fits; sumstats studies participate.",
      paste(names(qtlDatasets), collapse = ", ")))
  }

  if (length(nonStudyAxisSpecs) > 0L) {
    for (qdName in names(qtlDatasets)) {
      qd <- qtlDatasets[[qdName]]
      qdRes <- .twasDispatchJointSpecsQtlDataset(
        nonStudyAxisSpecs, qd, methods, contexts, traitIds, cisWindow,
        dataType, verbose)
      if (!is.null(qdRes))
        out <- if (is.null(out)) qdRes
               else .rbindTwasWeights(out, qdRes, ldSketch = NULL)
    }
  }

  if (!is.null(sumStats)) {
    ssRes <- .twasDispatchJointSpecsQtlSumStats(
      parsedJointSpec, sumStats, methods, contexts, traitIds, dataType,
      verbose)
    if (!is.null(ssRes)) {
      embeddedLd <- getLdSketch(ssRes)
      out <- if (is.null(out)) ssRes
             else .rbindTwasWeights(out, ssRes, ldSketch = embeddedLd)
    }
  } else if (length(studyAxisSpecs) > 0L && verbose >= 1) {
    message("jointCrossStudy (twas): no sumStats slot present on this MultiStudyQtlDataset; cross-study specs produce no result.")
  }
  out
}
