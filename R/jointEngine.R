# =============================================================================
# Joint-analysis engine (Phase 2; dev/jointSpecification-s4-refactor.md)
# -----------------------------------------------------------------------------
# Replaces the ~14 hand-written joint-dispatch leaf functions with: the uniform
# JointGroup contract (R/JointGroup.R), per-(dataForm, pipeline) `fitJointGroup`
# methods, one enumerator per (pattern, dataForm), one `.jointDispatchTable`
# wiring row per valid cell, and the `.runJointCell` engine.
#
# Identity model: a group's `conditions` data.frame (one row per Y/Z column)
# carries each fitted condition's (study, context, trait). The output row keying
# is DERIVED -- an axis that varies across conditions collapses to "joint" with
# members in jointStudies/jointContexts/jointTraits; a constant axis keeps its
# value. cross-context / cross-trait / cross-study are the single-varying-axis
# case; composed is >1 varying axis. Fitters are shared across patterns; only
# enumeration differs.
# =============================================================================

#' @include AllGenerics.R JointGroup.R
NULL

# ---- identity derivation ----------------------------------------------------

# The data-driven-prior LOOKUP key for a group's conditions: a varying (jointed)
# axis -> NULL (match-any, because the shared joint mr.mash fit lives on every
# per-context row), a constant axis -> its single value. Used only to find the
# mr.mash fit; the OUTPUT rows carry each condition's REAL (study, context,
# trait).
.jointPriorKey <- function(conditions) {
  axisVal <- function(ax) {
    u <- unique(as.character(conditions[[ax]]))
    if (length(u) > 1L) NULL else u[[1L]]
  }
  list(study = axisVal("study"), context = axisVal("context"),
       trait = axisVal("trait"))
}

# The ";"-joined distinct members of a varying axis (the per-row provenance tag
# jointStudies/Contexts/Traits), or NA when the axis is constant.
.jointAxisMembers <- function(conditions, ax) {
  u <- unique(as.character(conditions[[ax]]))
  if (length(u) > 1L) paste(u, collapse = ";") else NA_character_
}

# Slice a fine-mapping per-method CV payload (.fmSliceCv output:
# list(samplePartition, prediction = list(<m>_predicted = sample x condition),
# performance = list(<m>_performance = condition x 6))) down to one condition r,
# so each per-context FineMappingEntry carries that context's CV.
.fmSliceCvCondition <- function(cv, r) {
  if (is.null(cv)) return(NULL)
  out <- list(samplePartition = cv$samplePartition)
  if (!is.null(cv$prediction))
    out$prediction <- lapply(cv$prediction, function(m) m[, r, drop = FALSE])
  if (!is.null(cv$performance))
    out$performance <- lapply(cv$performance, function(m) m[r, , drop = FALSE])
  out
}

# Slice a twas joint cvResult (.jointTwasCvResult output: list(samplePartition,
# predictions = sample x condition, metrics = condition x 6, foldFits)) to one
# condition r. The per-fold mr.mash fits span all conditions, so foldFits is
# shared unchanged.
.sliceTwasCvResultToCondition <- function(cvRes, r) {
  if (is.null(cvRes)) return(NULL)
  list(samplePartition = cvRes$samplePartition,
       predictions = if (!is.null(cvRes$predictions))
                       cvRes$predictions[, r, drop = TRUE] else NULL,
       metrics     = if (!is.null(cvRes$metrics))
                       cvRes$metrics[r, , drop = TRUE] else NULL,
       foldFits    = cvRes$foldFits)
}

# Mutable accumulator for the joint rows the engine assembles; each add()
# appends one fitted group as one result row, deriving its identity + joint*
# members from the group's conditions.
.jointRows <- function() {
  e <- new.env(parent = emptyenv())
  e$study <- character(0); e$context <- character(0); e$trait <- character(0)
  e$method <- character(0); e$entries <- list()
  e$jointStudies <- character(0); e$jointContexts <- character(0)
  e$jointTraits <- character(0)
  e$add <- function(study, context, trait, method, entry,
                    jointStudies = NA_character_,
                    jointContexts = NA_character_,
                    jointTraits = NA_character_) {
    e$study   <- c(e$study, study)
    e$context <- c(e$context, context)
    e$trait   <- c(e$trait, trait)
    e$method  <- c(e$method, method)
    e$entries[[length(e$entries) + 1L]] <- entry
    e$jointStudies  <- c(e$jointStudies,  jointStudies)
    e$jointContexts <- c(e$jointContexts, jointContexts)
    e$jointTraits   <- c(e$jointTraits,   jointTraits)
  }
  e
}

# ---- fitters (fitJointGroup) ------------------------------------------------

# (individual, fine-mapping) -> mvSuSiE joint fit + honest per-fold CV prior.
setMethod("fitJointGroup", signature("IndividualJointGroup", "FmJointPipeline"),
  function(group, pipeline, token, args) {
    cfg <- pipeline@config
    Xc <- group@X; Yc <- group@Y; nCond <- ncol(Yc)
    # fsusie: functional SuSiE over the trait domain (cross-trait, individual-
    # level only -- no multi-context, no RSS). One per-condition entry per trait
    # (no data-driven prior); uses the enumerator-stored per-trait `pos`.
    if (identical(token, "fsusie")) {
      if (length(group@pos) != nCond)
        stop("fitJointGroup: fsusie requires per-trait positions ('pos'); ",
             "it is cross-trait individual-level only.")
      verbose <- if (is.null(cfg$verbose)) 1 else cfg$verbose
      fit <- do.call(fitFsusie,
                     .fmMergeUserArgs(list(X = Xc, Y = Yc, pos = group@pos),
                                      "fsusie", args$methodArgs[["fsusie"]]))
      # Collapse the functional fit to a variants x features weight matrix now,
      # while fitted_wc/csd_X are still present (trimming drops them); store on
      # $coef so a trimmed fit can still yield TWAS weights (port of mvJobs).
      fit$coef <- tryCatch(
        fsusieWeights(fsusieFit = fit, variantIds = colnames(Xc)),
        error = function(e) NULL)
      fit <- .setFinemappingFitClass(fit, "fsusie")
      cvM <- NULL
      cvFolds <- if (is.null(cfg$cvFolds)) 0L else cfg$cvFolds
      if (cvFolds > 1L) {
        cv <- .fmCrossValidate(Xc, Yc, "fsusie", args$methodArgs, cvFolds,
                               samplePartition = cfg$samplePartition,
                               coverage = cfg$coverage, pos = group@pos,
                               verbose = verbose)
        cvM <- .fmSliceCv(cv, "fsusie")
      }
      return(lapply(seq_len(nCond), function(r) {
        e <- .fmPostprocessOne(
          fit = fit, method = "fsusie", dataX = Xc, dataY = NULL, conditionIdx = r,
          coverage = cfg$coverage, secondaryCoverage = cfg$secondaryCoverage,
          signalCutoff = cfg$signalCutoff, minAbsCorr = cfg$minAbsCorr,
          csInput = "fsusie")
        if (!is.null(cvM)) e <- .fmAttachCv(e, .fmSliceCvCondition(cvM, r))
        e
      }))
    }
    if (!identical(token, "mvsusie"))
      stop("fitJointGroup(IndividualJointGroup, FmJointPipeline): unsupported ",
           "token '", token, "' (expected 'mvsusie' or 'fsusie').")
    ddCut <- if (is.null(cfg$dataDrivenPriorWeightsCutoff)) 1e-10
             else cfg$dataDrivenPriorWeightsCutoff
    verbose <- if (is.null(cfg$verbose)) 1 else cfg$verbose
    # SER pre-screen: drop conditions with no single-effect signal before the
    # joint fit (port of mvJobs' skipConditions). < 2 survivors -> skip the whole
    # joint (an all-NULL list, which .runJointCell turns into zero rows). The fit,
    # prior, and CV all run on the surviving columns only; each survivor's per-
    # condition posterior is later sliced at its position in the fitted set.
    keep <- rep(TRUE, nCond)
    if (.fmScreenActive(args$pipCutoffToSkip)) {
      keep <- as.logical(.fmSerScreenColumns(Xc, Yc, args$pipCutoffToSkip))
      if (sum(keep) < 2L) {
        if (verbose >= 1)
          message(sprintf(
            "Skipping mvsusie joint fit: < 2 of %d conditions pass the SER pre-screen.",
            nCond))
        return(vector("list", nCond))
      }
      if (sum(keep) < nCond && verbose >= 1)
        message(sprintf(
          "mvsusie joint fit: SER pre-screen kept %d of %d conditions.",
          sum(keep), nCond))
    }
    survivors <- which(keep)
    Ys <- Yc[, survivors, drop = FALSE]
    key <- .jointPriorKey(group@conditions)
    mvFitParts <- .fmLookupMrmashFit(args$twasWeights, key$study, key$trait,
                                     context = key$context)
    mvCv <- .fmLookupMrmashCv(args$twasWeights, key$study, key$trait,
                              context = key$context)
    mvPrior <- .buildMvsusieReweightedPrior(mvFitParts, colnames(Ys), ddCut)
    mvBaseArgs <- list(X = Xc, Y = Ys,
                       prior_variance = mvPrior$priorVariance,
                       coverage = cfg$coverage)
    if (!is.null(mvPrior$residualVariance))
      mvBaseArgs$residual_variance <- mvPrior$residualVariance
    fit <- do.call(fitMvsusie,
                   .fmMergeUserArgs(mvBaseArgs, "mvsusie",
                                    args$methodArgs[["mvsusie"]]))
    fit <- .setFinemappingFitClass(fit, "mvsusie")
    cvM <- NULL
    cvFolds <- if (is.null(cfg$cvFolds)) 0L else cfg$cvFolds
    if (cvFolds > 1L) {
      sp <- cfg$samplePartition
      if (is.null(sp) && !is.null(mvCv)) sp <- mvCv$samplePartition
      mvPriorCv <- .fmBuildMvsusiePriorCv(mvCv, mvFitParts, colnames(Ys), ddCut)
      cv <- .fmCrossValidate(Xc, Ys, "mvsusie", args$methodArgs, cvFolds,
                             samplePartition = sp, coverage = cfg$coverage,
                             verbose = verbose, mvPrior = mvPrior,
                             mvPriorCv = mvPriorCv)
      cvM <- .fmSliceCv(cv, "mvsusie")
    }
    # One per-context entry per ORIGINAL condition: NULL for screened-out columns
    # (skipped downstream), else that condition's posterior (sliced at its
    # position in the fitted survivor set) with shared pip/cs + its CV slice.
    lapply(seq_len(nCond), function(i) {
      if (!keep[i]) return(NULL)
      r <- match(i, survivors)
      e <- .fmPostprocessOne(
        fit = fit, method = "mvsusie", dataX = Xc, dataY = NULL, conditionIdx = r,
        coverage = cfg$coverage, secondaryCoverage = cfg$secondaryCoverage,
        signalCutoff = cfg$signalCutoff, minAbsCorr = cfg$minAbsCorr,
        csInput = "X")
      if (!is.null(cvM)) e <- .fmAttachCv(e, .fmSliceCvCondition(cvM, r))
      e
    })
  })

# (sumstats, fine-mapping) -> mvSuSiE-rss joint fit. RSS has no sample folds and
# no fsusie variant.
setMethod("fitJointGroup", signature("SumStatsJointGroup", "FmJointPipeline"),
  function(group, pipeline, token, args) {
    if (identical(token, "fsusie"))
      stop("fsusie has no RSS variant; it requires individual-level input.")
    if (!identical(token, "mvsusie"))
      stop("fitJointGroup(SumStatsJointGroup, FmJointPipeline): unsupported ",
           "token '", token, "' (expected 'mvsusie').")
    cfg <- pipeline@config
    ddCut <- if (is.null(cfg$dataDrivenPriorWeightsCutoff)) 1e-10
             else cfg$dataDrivenPriorWeightsCutoff
    key <- .jointPriorKey(group@conditions)
    mvFitParts <- .fmLookupMrmashFit(args$twasWeights, key$study, key$trait,
                                     context = key$context)
    mvPrior <- .buildMvsusieReweightedPrior(mvFitParts, colnames(group@Z), ddCut)
    mvBaseArgs <- list(Z = group@Z, R = group@R,
                       N = as.numeric(stats::median(group@N)),
                       prior_variance = mvPrior$priorVariance,
                       coverage = cfg$coverage)
    if (!is.null(mvPrior$residualVariance))
      mvBaseArgs$residual_variance <- mvPrior$residualVariance
    fit <- do.call(fitMvsusieRss,
                   .fmMergeUserArgs(mvBaseArgs, "mvsusie",
                                    args$methodArgs[["mvsusie"]]))
    fit <- .setFinemappingFitClass(fit, "mvsusie")
    # One per-condition entry (RSS has no sample folds).
    lapply(seq_len(ncol(group@Z)), function(r) .fmPostprocessOne(
      fit = fit, method = "mvsusie", dataX = group@R, dataY = NULL,
      conditionIdx = r, coverage = cfg$coverage,
      secondaryCoverage = cfg$secondaryCoverage, signalCutoff = cfg$signalCutoff,
      minAbsCorr = cfg$minAbsCorr, csInput = "Xcorr"))
  })

# Reshape a twasWeightsCv() result into the single joint entry's cvResult: the
# out-of-fold prediction matrix, the per-condition metric rows, and the per-fold
# mr.mash fits (named fold_<j>) that fineMappingPipeline's mvSuSiE path consumes.
.jointTwasCvResult <- function(cv, token) {
  if (is.null(cv)) return(NULL)
  pickByBase <- function(lst) {
    if (is.null(lst) || length(lst) == 0L) return(NULL)
    bare <- sub("(_predicted|Predicted|_performance|Performance)$", "", names(lst))
    hit <- which(bare == token)
    if (length(hit) == 0L) NULL else lst[[hit[[1L]]]]
  }
  ffKey <- paste0(token, "_weights")
  foldFits <- if (!is.null(cv$foldFits)) {
    ff <- lapply(cv$foldFits, function(f) f[[ffKey]])
    if (all(vapply(ff, is.null, logical(1)))) NULL else ff
  } else NULL
  list(samplePartition = cv$samplePartition,
       predictions     = pickByBase(cv$prediction),
       metrics         = pickByBase(cv$performance),
       foldFits        = foldFits)
}

# learnTwasWeights key for a bare token (fine-mapping tokens key differently,
# e.g. susieInf -> susie_inf_weights).
.twasMethodKey <- function(token) {
  ad <- .twasFineMappingMethodAdapters[[token]]
  if (!is.null(ad)) ad$methodKey else paste0(token, "_weights")
}

# Fine-mapping CV handoff for one twas method: extract that method's out-of-fold
# predictions + performance from fineMappingPipeline's retained CV (shared fold
# partition), shaped like .jointTwasCvResult so the per-condition slice reuses it
# instead of re-cross-validating an FM-derived method (susie / mvsusie / ...).
.twasFmHandoffCv <- function(fineMappingCv, token) {
  if (is.null(fineMappingCv) || is.null(fineMappingCv$prediction)) return(NULL)
  base <- sub("(_predicted|Predicted)$", "", names(fineMappingCv$prediction))
  hit <- which(base == token)
  if (length(hit) == 0L) return(NULL)
  pBase <- sub("(_performance|Performance)$", "", names(fineMappingCv$performance))
  pHit <- which(pBase == token)
  list(samplePartition = fineMappingCv$samplePartition,
       predictions = fineMappingCv$prediction[[hit[[1L]]]],
       metrics = if (length(pHit)) fineMappingCv$performance[[pHit[[1L]]]] else NULL,
       foldFits = NULL)
}

# (individual, twas) -> ONE weight method fit over the group's conditions, as
# per-condition entries (sliced from the variants x conditions weight matrix),
# each with its full-data weights + retained fit + per-condition CV slice. This
# is the SHARED per-method twas fitting (one method per call, like the FM
# fitters); the SR-TWAS ensemble combines methods in a layer above (see
# .twasEnsembleLayer). Owns the orchestration formerly in .twasWeightsPipelineMatrix:
# FM-fit injection (FM-derived tokens extract from the precomputed fit), the FM
# CV handoff (reuse fine-mapping's own CV), spike-and-slab pi from an internal
# mr.ash fit, CV knobs, and fitFullData = FALSE (CV-only) entries.
setMethod("fitJointGroup", signature("IndividualJointGroup", "TwasJointPipeline"),
  function(group, pipeline, token, args) {
    cfg <- pipeline@config
    Xc <- group@X; Yc <- group@Y; nCond <- ncol(Yc)
    cond <- group@conditions
    fitFullData <- if (is.null(cfg$fitFullData)) TRUE else isTRUE(cfg$fitFullData)
    cvFolds <- if (is.null(cfg$cvFolds)) 0L else cfg$cvFolds
    rfd <- if (is.null(cfg$retainFitDetail)) "slim" else cfg$retainFitDetail
    verbose <- if (is.null(cfg$verbose)) 1 else cfg$verbose
    stdz <- isTRUE(cfg$standardized)
    estimatePi <- isTRUE(cfg$estimatePi)
    methodKey <- .twasMethodKey(token)
    # Method args: prefer the full methodList (the unified pipeline path), fall
    # back to methodArgs (the explicit-jointSpec dispatchers).
    ma <- if (!is.null(args$methodList) && methodKey %in% names(args$methodList))
            args$methodList[[methodKey]]
          else if (!is.null(args$methodArgs)) args$methodArgs[[methodKey]] else NULL
    if (is.null(ma)) ma <- list()
    fittedModels <- if (!is.null(args$fittedModels)) args$fittedModels else list()
    # FM-fit injection: an FM-derived token extracts its weights from the
    # precomputed fine-mapping fit rather than refitting.
    adapter <- .twasFineMappingMethodAdapters[[token]]
    if (!is.null(adapter) && !is.null(fittedModels[[token]]) &&
        is.null(ma[[adapter$fitArg]])) {
      ma[[adapter$fitArg]] <- fittedModels[[token]]
    }
    # Spike-and-slab pi from an internal mr.ash fit (self-contained per method).
    if (estimatePi && token %in% c("bayes_c", "bayes_b")) {
      mrA <- learnTwasWeights(Xc, Yc, weightMethods = list(mrash_weights = list()),
               study = as.character(cond$study[1L]),
               context = as.character(cond$context[1L]),
               trait = as.character(cond$trait[1L]),
               retainFits = TRUE, standardized = stdz,
               dataType = cfg$dataType, verbose = 0)
      piHat <- as.numeric(estimateSparsity(mrA))
      if (token == "bayes_c" && is.null(ma$pi))     ma$pi     <- piHat
      if (token == "bayes_b" && is.null(ma$probIn)) ma$probIn <- piHat
    }
    wm <- setNames(list(ma), methodKey)

    W <- NULL; fitParts <- NULL; vids <- colnames(Xc)
    if (fitFullData) {
      tw <- learnTwasWeights(Xc, Yc, weightMethods = wm,
              study = as.character(cond$study[1L]),
              context = as.character(cond$context[1L]),
              trait = as.character(cond$trait[1L]),
              fittedModels = fittedModels,
              retainFits = TRUE, retainFitDetail = rfd,
              standardized = stdz, dataType = cfg$dataType, verbose = verbose)
      base <- tw$entry[[1L]]
      W <- getWeights(base)
      if (!is.matrix(W))
        W <- matrix(W, ncol = nCond, dimnames = list(getVariantIds(base), NULL))
      fitParts <- getFits(base); vids <- getVariantIds(base)
    }
    cvRes <- NULL
    if (cvFolds > 1L) {
      # FM-derived method: reuse fine-mapping's own CV (shared partition); a
      # method whose full-data weights are all zero is skipped (nothing to CV).
      cvRes <- .twasFmHandoffCv(args$fineMappingCv, token)
      if (is.null(cvRes) && !(!is.null(W) && all(W == 0))) {
        # Leakage guard: a single full-data data-driven mr.mash prior reused
        # across folds means each fold's prior saw its own held-out samples.
        if (is.null(args$dataDrivenPriorMatricesCv) &&
            !is.null(ma$dataDrivenPriorMatrices)) {
          warning("Cross-validating mr.mash with a single data-driven prior ",
                  "computed on the full data: the same prior is reused for ",
                  "every fold, so each fold's prior was informed by its own ",
                  "held-out samples (information leakage). Supply per-fold ",
                  "priors via dataDrivenPriorMatricesCv (--mixture-prior-cv) ",
                  "for honest cross-validation.")
        }
        sp <- if (!is.null(args$samplePartition)) args$samplePartition
              else cfg$samplePartition
        mcv <- if (is.null(cfg$maxCvVariants) || cfg$maxCvVariants <= 0) Inf
               else cfg$maxCvVariants
        cv <- twasWeightsCv(Xc, Yc, fold = cvFolds, samplePartitions = sp,
                weightMethods = wm, retainFits = TRUE, maxNumVariants = mcv,
                numThreads = if (is.null(cfg$cvThreads)) 1 else cfg$cvThreads,
                data_driven_prior_matrices_cv = args$dataDrivenPriorMatricesCv,
                verbose = verbose)
        cvRes <- .jointTwasCvResult(cv, token)
      }
    }
    # One per-condition entry: that condition's weight column + the shared fit +
    # its CV slice. fitFullData = FALSE -> CV-only entry.
    lapply(seq_len(nCond), function(r) {
      cvR <- if (!is.null(cvRes)) .sliceTwasCvResultToCondition(cvRes, r) else NULL
      if (is.null(W)) {
        TwasWeightsEntry(variantIds = character(0), weights = NULL,
                         cvResult = cvR, standardized = stdz,
                         dataType = cfg$dataType)
      } else {
        TwasWeightsEntry(variantIds = vids, weights = W[, r], fits = fitParts,
                         cvResult = cvR, standardized = stdz,
                         dataType = cfg$dataType)
      }
    })
  })

# (sumstats, twas) -> mr.mash-rss joint fit as ONE matrix entry. No sample folds.
setMethod("fitJointGroup", signature("SumStatsJointGroup", "TwasJointPipeline"),
  function(group, pipeline, token, args) {
    cfg <- pipeline@config
    rfd <- if (is.null(cfg$retainFitDetail)) "slim" else cfg$retainFitDetail
    weights <- mrmashRssWeights(stat = list(z = group@Z, N = group@N),
                                LD = group@R, retainFit = TRUE, fitDetail = rfd)
    vids <- rownames(weights); if (is.null(vids)) vids <- rownames(group@Z)
    fitParts <- attr(weights, "fit")
    if (!is.matrix(weights))
      weights <- matrix(weights, ncol = ncol(group@Z), dimnames = list(vids, NULL))
    # One per-condition entry: that condition's weight column + the shared fit.
    lapply(seq_len(ncol(weights)), function(r)
      TwasWeightsEntry(variantIds = vids, weights = weights[, r], fits = fitParts,
                       standardized = TRUE, dataType = cfg$dataType))
  })

# ---- result construction (construct) ----------------------------------------

# Both pipelines assemble identically-shaped joint rows; only the result
# collection differs (the axis-3 divergence the markers encode). Only the joint*
# columns for axes that actually vary are attached.
.constructJointArgs <- function(pipeline, rows) {
  a <- list(study = rows$study, context = rows$context, trait = rows$trait,
            method = rows$method, entry = rows$entries,
            ldSketch = pipeline@config$ldSketch)
  if (any(!is.na(rows$jointStudies)))  a$jointStudies  <- rows$jointStudies
  if (any(!is.na(rows$jointContexts))) a$jointContexts <- rows$jointContexts
  if (any(!is.na(rows$jointTraits)))   a$jointTraits   <- rows$jointTraits
  a
}

setMethod("construct", "FmJointPipeline",
  function(pipeline, rows, ...) {
    if (length(rows$entries) == 0L) return(NULL)
    do.call(QtlFineMappingResult, .constructJointArgs(pipeline, rows))
  })

setMethod("construct", "TwasJointPipeline",
  function(pipeline, rows, ...) {
    if (length(rows$entries) == 0L) return(NULL)
    do.call(TwasWeights, .constructJointArgs(pipeline, rows))
  })

# ---- enumerators (pattern x dataForm -> list<JointGroup>) --------------------

# cross-context / individual: one group per scoped trait present in >= 2 scoped
# contexts (the conditions are those (study, context, trait) rows).
.enumCrossContextIndividual <- function(data, scope, args = list()) {
  study <- getStudy(data)
  if (!(study %in% scope$studies)) return(list())
  scopedContexts <- scope$contexts[[study]]
  scopedTraits   <- scope$traits[[study]]
  if (length(scopedContexts) < 2L) return(list())
  verbose <- if (is.null(args$verbose)) 1 else args$verbose
  groups <- list()
  for (tid in scopedTraits) {
    xy <- .buildIndividualCrossContextXY(
      data, tid, scopedContexts, args$cisWindow, verbose,
      label = "jointCrossContext", region = args$region)
    if (is.null(xy)) next
    groups[[length(groups) + 1L]] <- new("IndividualJointGroup",
      conditions = data.frame(study = study, context = xy$perTraitContexts,
                              trait = tid, stringsAsFactors = FALSE),
      X = xy$X, Y = xy$Y)
  }
  groups
}

# cross-context / sumstats.
.enumCrossContextSumstats <- function(data, scope, args = list()) {
  ldSketch   <- getLdSketch(data)
  studyCol   <- as.character(data$study)
  contextCol <- as.character(data$context)
  traitCol   <- as.character(data$trait)
  groups <- list()
  for (s in scope$studies) {
    scopedContexts <- scope$contexts[[s]]
    scopedTraits   <- scope$traits[[s]]
    if (length(scopedContexts) < 2L) next
    for (tid in scopedTraits) {
      tupleRows <- which(studyCol == s & traitCol == tid &
                         contextCol %in% scopedContexts)
      if (length(tupleRows) < 2L) next
      ctxNames <- contextCol[tupleRows]
      jz <- .buildJointSumstatZMatrix(
        data, tupleRows, ctxNames,
        errorLabel = "jointCrossContext (QtlSumStats)")
      ldMat <- .fmLdFromSketch(ldSketch, jz$variantIds)
      groups[[length(groups) + 1L]] <- new("SumStatsJointGroup",
        conditions = data.frame(study = s, context = ctxNames, trait = tid,
                                stringsAsFactors = FALSE),
        Z = jz$Z, R = ldMat, N = jz$nVec)
    }
  }
  groups
}

# cross-trait / individual: one group per scoped context with >= 2 scoped traits.
.enumCrossTraitIndividual <- function(data, scope, args = list()) {
  study <- getStudy(data)
  if (!(study %in% scope$studies)) return(list())
  scopedContexts <- scope$contexts[[study]]
  scopedTraits   <- scope$traits[[study]]
  verbose <- if (is.null(args$verbose)) 1 else args$verbose
  groups <- list()
  for (cx in scopedContexts) {
    xy <- .buildIndividualCrossTraitXY(
      data, cx, scopedTraits, args$cisWindow, verbose,
      label = "jointCrossTrait", study = study, region = args$region)
    if (is.null(xy)) next
    # Functional positions (one per trait column) for fsusie's domain; mvsusie
    # ignores them. Matches the trait order of Y.
    rr  <- SummarizedExperiment::rowRanges(xy$se)
    rr  <- rr[match(colnames(xy$Y), rownames(xy$se))]
    pos <- (GenomicRanges::start(rr) + GenomicRanges::end(rr)) / 2
    groups[[length(groups) + 1L]] <- new("IndividualJointGroup",
      conditions = data.frame(study = study, context = cx,
                              trait = xy$traitsHere, stringsAsFactors = FALSE),
      X = xy$X, Y = xy$Y, pos = as.numeric(pos))
  }
  groups
}

# cross-trait / sumstats.
.enumCrossTraitSumstats <- function(data, scope, args = list()) {
  ldSketch   <- getLdSketch(data)
  studyCol   <- as.character(data$study)
  contextCol <- as.character(data$context)
  traitCol   <- as.character(data$trait)
  groups <- list()
  for (s in scope$studies) {
    scopedContexts <- scope$contexts[[s]]
    scopedTraits   <- scope$traits[[s]]
    for (cx in scopedContexts) {
      tupleRows <- which(studyCol == s & contextCol == cx &
                         traitCol %in% scopedTraits)
      if (length(tupleRows) < 2L) next
      trNames <- traitCol[tupleRows]
      jz <- .buildJointSumstatZMatrix(
        data, tupleRows, trNames, errorLabel = "jointCrossTrait (QtlSumStats)")
      ldMat <- .fmLdFromSketch(ldSketch, jz$variantIds)
      groups[[length(groups) + 1L]] <- new("SumStatsJointGroup",
        conditions = data.frame(study = s, context = cx, trait = trNames,
                                stringsAsFactors = FALSE),
        Z = jz$Z, R = ldMat, N = jz$nVec)
    }
  }
  groups
}

# cross-study / sumstats (no individual form: individual-level studies have
# disjoint samples). One group per (context, trait) present in >= 2 scoped
# studies; the study axis varies -> "joint" + jointStudies.
.enumCrossStudySumstats <- function(data, scope, args = list()) {
  ldSketch   <- getLdSketch(data)
  studyCol   <- as.character(data$study)
  contextCol <- as.character(data$context)
  traitCol   <- as.character(data$trait)
  allCtxs <- unique(unlist(scope$contexts, use.names = FALSE))
  allTrs  <- unique(unlist(scope$traits,   use.names = FALSE))
  groups <- list()
  for (cx in allCtxs) {
    for (tid in allTrs) {
      tupleRows <- which(contextCol == cx & traitCol == tid &
                         studyCol %in% scope$studies)
      keep <- vapply(tupleRows, function(r) {
        s <- studyCol[r]
        (cx %in% scope$contexts[[s]]) && (tid %in% scope$traits[[s]])
      }, logical(1))
      tupleRows <- tupleRows[keep]
      if (length(tupleRows) < 2L) next
      stNames <- studyCol[tupleRows]
      jz <- .buildJointSumstatZMatrix(
        data, tupleRows, stNames, errorLabel = "jointCrossStudy")
      ldMat <- .fmLdFromSketch(ldSketch, jz$variantIds)
      groups[[length(groups) + 1L]] <- new("SumStatsJointGroup",
        conditions = data.frame(study = stNames, context = cx, trait = tid,
                                stringsAsFactors = FALSE),
        Z = jz$Z, R = ldMat, N = jz$nVec)
    }
  }
  groups
}

# composed / individual: ONE group joining every scoped (context, trait) tuple
# for the study. Both context and trait vary across conditions, so both collapse
# to "joint" (the conditions model handles multi-varying-axis uniformly; if the
# tuples happen to share a context it degrades to cross-trait, and vice versa).
.enumComposedIndividual <- function(data, scope, args = list()) {
  study <- getStudy(data)
  if (!(study %in% scope$studies)) return(list())
  verbose <- if (is.null(args$verbose)) 1 else args$verbose
  xy <- .buildComposedIndividualXY(
    data, scope, study, args$cisWindow, verbose,
    label = "composed", region = args$region)
  if (is.null(xy)) return(list())
  # Conditions follow the fitted Y columns ("context:trait"), so dropped tuples
  # don't desync conditions from Y. Split on the first ":" (contexts are simple
  # labels; trait ids may themselves contain ":").
  labs <- colnames(xy$Y)
  conds <- data.frame(
    study   = study,
    context = sub(":.*$", "", labs),
    trait   = sub("^[^:]*:", "", labs),
    stringsAsFactors = FALSE)
  list(new("IndividualJointGroup", conditions = conds, X = xy$X, Y = xy$Y))
}

# univariate / individual: one 1-condition group per (study, context, trait) in
# scope -- the per-(context, trait) iteration expressed as engine groups, so
# univariate methods (lasso / enet / susie / ...) flow through the SAME per-
# method fitter + ensemble layer as the joint ones (minGroup = 1).
.enumUnivariateIndividual <- function(data, scope, args = list()) {
  study <- getStudy(data)
  if (!(study %in% scope$studies)) return(list())
  naAction <- if (is.null(args$naAction)) "drop" else args$naAction
  groups <- list()
  for (cx in scope$contexts[[study]]) {
    se <- getPhenotypes(data, contexts = cx)
    for (tid in intersect(scope$traits[[study]], rownames(se))) {
      Y <- .fmResidPheno(data, contexts = cx, traitId = tid, naAction = naAction)
      X <- if (is.null(args$region))
             .fmResidGeno(data, contexts = cx, traitId = tid,
                          cisWindow = args$cisWindow)
           else .fmResidGeno(data, contexts = cx, region = args$region)
      common <- intersect(rownames(X), rownames(Y))
      if (length(common) < 2L) next
      groups[[length(groups) + 1L]] <- new("IndividualJointGroup",
        conditions = data.frame(study = study, context = cx, trait = tid,
                                stringsAsFactors = FALSE),
        X = X[common, , drop = FALSE], Y = Y[common, , drop = FALSE])
    }
  }
  groups
}

# composed / sumstats: general N-axis joint. `args$axes` (subset of study /
# context / trait) names the collapsed axes; rows split by the complement
# (fixed) axes form one group each. Reuses .enumerateComposedSumstatGroups.
.enumComposedSumstats <- function(data, scope, args = list()) {
  axes <- args$axes
  if (is.null(axes)) axes <- c("context", "trait")
  ldSketch <- getLdSketch(data)
  gi <- .enumerateComposedSumstatGroups(list(axes = axes), data, scope)
  if (is.null(gi)) return(list())
  groups <- list()
  for (gIdx in gi$groups) {
    if (length(gIdx) < 2L) next
    colLabels <- vapply(gIdx, function(i)
      paste(gi$studyCol[i], gi$contextCol[i], gi$traitCol[i], sep = ":"),
      character(1L))
    jz <- .buildJointSumstatZMatrix(
      data, gIdx, colLabels, errorLabel = "composed (QtlSumStats)")
    ldMat <- .fmLdFromSketch(ldSketch, jz$variantIds)
    groups[[length(groups) + 1L]] <- new("SumStatsJointGroup",
      conditions = data.frame(study = gi$studyCol[gIdx],
                              context = gi$contextCol[gIdx],
                              trait = gi$traitCol[gIdx], stringsAsFactors = FALSE),
      Z = jz$Z, R = ldMat, N = jz$nVec)
  }
  groups
}

# ---- engine -----------------------------------------------------------------

# Twas per-group args: resolve the group's fine-mapping fits + CV (keyed on its
# first condition -- the joint fit is shared across conditions) and fix ONE
# shared fold partition (so every method's out-of-fold CV predictions align for
# the ensemble layer). Returns `args` unchanged for fine-mapping pipelines.
.twasGroupArgs <- function(g, pipeline, args) {
  if (!is(pipeline, "TwasJointPipeline")) return(args)
  cfg <- pipeline@config
  cond <- g@conditions
  out <- args
  fmRes <- args$fineMappingResult
  if (!is.null(fmRes)) {
    s1 <- as.character(cond$study[[1L]])
    c1 <- as.character(cond$context[[1L]])
    t1 <- as.character(cond$trait[[1L]])
    nR <- if (is.null(args$nRegions))   1L else args$nRegions
    bi <- if (is.null(args$regionIndex)) 1L else args$regionIndex
    af <- .twasFineMappingFits(fmRes, study = s1, context = c1, trait = t1)
    out$fittedModels  <- if (is.null(af)) list() else .twasFitsForRegion(af, bi, nR)
    out$fineMappingCv <- .twasCvResultFor(fmRes, s1, c1, t1)
  }
  cvF <- if (is.null(cfg$cvFolds)) 0L else cfg$cvFolds
  if (cvF > 1L && is(g, "IndividualJointGroup")) {
    sp <- args$samplePartition
    if (is.null(sp)) sp <- cfg$samplePartition
    if (is.null(sp) && !is.null(out$fineMappingCv))
      sp <- out$fineMappingCv$samplePartition
    if (is.null(sp))
      sp <- .normalizeCvFolds(cvF, NULL, rownames(g@X))$samplePartition
    out$samplePartition <- sp
  }
  out
}

# Run one dispatch cell: enumerate joint groups, fit each method (S4 dispatch on
# the group x pipeline pair) per group, accumulate per-context rows, build the
# per-pipeline result. The loop is GROUP-outer / token-inner so the twas ensemble
# layer can combine a group's per-method fits in place (FM is unaffected by the
# loop order). Per-method fitting is identical for FM and twas -- one method ->
# per-condition entries; the SR-TWAS ensemble is a layer ON TOP of that.
.runJointCell <- function(cell, pipeline, data, scope, tokens, args = list()) {
  groups <- cell@enumerate(data, scope, args)
  groups <- Filter(function(g) nrow(g@conditions) >= cell@minGroup, groups)
  if (length(groups) == 0L) return(NULL)
  doEnsemble <- is(pipeline, "TwasJointPipeline") &&
                isTRUE(pipeline@config$ensemble)
  rows <- .jointRows()
  for (g in groups) {
    cond <- g@conditions
    # Provenance: the ";"-joined members of each varying axis, identical on every
    # per-context row of this joint group.
    js <- .jointAxisMembers(cond, "study")
    jc <- .jointAxisMembers(cond, "context")
    jt <- .jointAxisMembers(cond, "trait")
    addEntries <- function(entries, method) {
      for (i in seq_len(min(length(entries), nrow(cond)))) {
        e <- entries[[i]]
        if (is.null(e)) next
        rows$add(study = as.character(cond$study[[i]]),
                 context = as.character(cond$context[[i]]),
                 trait = as.character(cond$trait[[i]]),
                 method = method, entry = e,
                 jointStudies = js, jointContexts = jc, jointTraits = jt)
      }
    }
    # Twas: resolve this group's fine-mapping fits + CV (keyed on its first
    # condition; the joint fit is shared across conditions) and fix ONE fold
    # partition up front, so every method's out-of-fold CV predictions are
    # aligned for the ensemble layer. FM leaves args untouched.
    fitArgs <- .twasGroupArgs(g, pipeline, args)
    # Per-method fit -> per-condition entries -> rows (shared FM + twas). Retain
    # each method's entries so the twas ensemble layer can combine them.
    perTokenEntries <- list()
    for (token in tokens) {
      # Resume cache: if every condition of this group is already present in the
      # prior partial result (args$cache), reuse those entries instead of
      # refitting. All-or-nothing per group. FM passes a QtlFineMappingResult;
      # twas passes a TwasWeights -- the lookup matches the pipeline.
      entries <- NULL
      if (!is.null(args$cache)) {
        lookup <- if (is(pipeline, "TwasJointPipeline")) .twasCacheLookup
                  else .fmCacheLookup
        cached <- lapply(seq_len(nrow(cond)), function(i)
          lookup(args$cache, as.character(cond$study[[i]]),
                 as.character(cond$context[[i]]),
                 as.character(cond$trait[[i]]), token))
        if (!any(vapply(cached, is.null, logical(1)))) entries <- cached
      }
      if (is.null(entries)) entries <- fitJointGroup(g, pipeline, token, fitArgs)
      if (is.null(entries) || length(entries) == 0L) next
      perTokenEntries[[token]] <- entries
      addEntries(entries, token)
    }
    # SR-TWAS ensemble layer: combine the group's per-method per-condition fits
    # (CV predictions + weights) into ensemble per-context rows -- built ON TOP
    # of the shared per-method fitting above, never inside it.
    if (doEnsemble && length(perTokenEntries) >= 2L) {
      addEntries(.twasEnsembleLayer(g, perTokenEntries, pipeline@config),
                 "ensemble")
    }
  }
  construct(pipeline, rows)
}

# SR-TWAS ensemble LAYER (twas only): combine a group's per-method per-condition
# fits into ensemble per-condition entries -- built ON TOP of the shared per-
# method fitting, never inside it. For each condition r, gather the methods'
# retained out-of-fold CV predictions + weights + R^2, drop methods below the
# R^2 cutoff (stacking needs >= 2), and combine via the `ensembleWeights`
# primitive PER CONTEXT (the sliced single-condition inputs -> contextIndex = 1).
# Returns a length-nCond list of ensemble TwasWeightsEntry (NULL where < 2
# methods qualify). All methods share the group's fold partition (the runner
# fixes it before fitting), so their out-of-fold predictions are comparable.
.twasEnsembleLayer <- function(group, perTokenEntries, cfg) {
  tokens <- names(perTokenEntries)
  Y <- group@Y
  r2Cut  <- if (is.null(cfg$ensembleR2Threshold)) 0.01 else cfg$ensembleR2Threshold
  solver <- if (is.null(cfg$ensembleSolver)) "quadprog" else cfg$ensembleSolver
  alpha  <- if (is.null(cfg$ensembleAlpha)) 1 else cfg$ensembleAlpha
  stdz   <- isTRUE(cfg$standardized)
  lapply(seq_len(nrow(group@conditions)), function(r) {
    preds <- list(); wts <- list(); rsq <- c()
    for (tk in tokens) {
      e <- perTokenEntries[[tk]][[r]]
      if (is.null(e)) next
      cv <- getCvResult(e); w <- getWeights(e)
      if (is.null(cv) || is.null(cv$predictions) || is.null(w)) next
      pr <- cv$predictions
      preds[[paste0(tk, "_predicted")]] <-
        matrix(as.numeric(pr), ncol = 1L, dimnames = list(names(pr), NULL))
      wts[[paste0(tk, "_weights")]] <-
        matrix(as.numeric(w), ncol = 1L, dimnames = list(getVariantIds(e), NULL))
      mt <- cv$metrics
      rsq[tk] <- if (!is.null(mt) && "rsq" %in% names(mt)) mt[["rsq"]] else NA_real_
    }
    passing <- names(rsq)[!is.na(rsq) & rsq >= r2Cut]
    if (length(passing) < 2L) return(NULL)
    ens <- tryCatch(ensembleWeights(
      cvResults = list(prediction = preds[paste0(passing, "_predicted")]),
      Y = Y[, r], twasWeightList = wts[paste0(passing, "_weights")],
      contextIndex = 1, solver = solver, alpha = alpha),
      error = function(err) NULL)
    if (is.null(ens) || is.null(ens$ensembleTwasWeights)) return(NULL)
    ew <- ens$ensembleTwasWeights
    vids <- if (!is.null(names(ew))) names(ew) else rownames(ew)
    if (is.null(vids)) vids <- getVariantIds(perTokenEntries[[passing[1L]]][[r]])
    TwasWeightsEntry(
      variantIds = vids, weights = as.numeric(ew),
      cvResult = list(methodCoef = ens$methodCoef,
                      methodPerformance = ens$methodPerformance),
      standardized = stdz, dataType = cfg$dataType)
  })
}

# ---- wiring table -----------------------------------------------------------
# Valid cells are rows; invalid cells are absences (a lookup miss is the error).
.jointDispatchTable <- list(
  new("JointDispatchCell", pattern = "context", dataForm = "individual",
      enumerate = .enumCrossContextIndividual, minGroup = 2L),
  new("JointDispatchCell", pattern = "context", dataForm = "sumstats",
      enumerate = .enumCrossContextSumstats, minGroup = 2L),
  new("JointDispatchCell", pattern = "trait", dataForm = "individual",
      enumerate = .enumCrossTraitIndividual, minGroup = 2L),
  new("JointDispatchCell", pattern = "trait", dataForm = "sumstats",
      enumerate = .enumCrossTraitSumstats, minGroup = 2L),
  new("JointDispatchCell", pattern = "study", dataForm = "sumstats",
      enumerate = .enumCrossStudySumstats, minGroup = 2L),
  new("JointDispatchCell", pattern = "composed", dataForm = "individual",
      enumerate = .enumComposedIndividual, minGroup = 2L),
  new("JointDispatchCell", pattern = "composed", dataForm = "sumstats",
      enumerate = .enumComposedSumstats, minGroup = 2L),
  # Univariate: per-(context, trait) 1-condition groups (twas individual only),
  # so univariate methods route through the same engine fitter + ensemble layer.
  new("JointDispatchCell", pattern = "univariate", dataForm = "individual",
      enumerate = .enumUnivariateIndividual, minGroup = 1L)
)

.lookupJointCell <- function(pattern, dataForm) {
  for (cell in .jointDispatchTable)
    if (cell@pattern == pattern && cell@dataForm == dataForm) return(cell)
  stop(sprintf("No joint dispatch cell for pattern='%s', dataForm='%s'.",
               pattern, dataForm))
}

# Run a parsed jointSpecification through the engine: for each spec resolve its
# scope, map its axes to a (pattern, dataForm) cell, and run every requested
# joint method (token) through `.runJointCell`, rbinding the per-spec results.
# Shared by the fm + twas QtlDataset / QtlSumStats / MultiStudy dispatchers --
# the marker (pipeline) selects the result type and the rbind. `args` is the
# per-run engine payload (twasWeights, methodArgs, cisWindow, region, ...).
.runJointSpecs <- function(parsedJointSpec, data, dataForm, pipeline,
                           jointMethods, contexts, traitIds, args = list()) {
  if (length(jointMethods) == 0L || length(parsedJointSpec) == 0L) return(NULL)
  ldSketch <- pipeline@config$ldSketch
  isFm <- is(pipeline, "FmJointPipeline")
  out <- NULL
  for (spec in parsedJointSpec) {
    scope <- .fmResolveSpecScope(spec, data, contexts = contexts,
                                 traitIds = traitIds)
    # Region mode WITHOUT an explicit traitId: restrict scoped traits to the
    # genes overlapping the locus (matches fineMappingPipeline's univariate
    # region trait selection, where traitId -- when given -- already pins the
    # gene set and takes precedence over region). Gene coordinates are context-
    # independent, so the first scoped context's SE provides them.
    if (dataForm == "individual" && is.null(traitIds) && !is.null(args$region)) {
      for (st in names(scope$traits)) {
        ctxs <- scope$contexts[[st]]
        if (length(ctxs) == 0L) next
        se <- getPhenotypes(data, contexts = ctxs[[1L]])
        scope$traits[[st]] <- .fmTraitsInRegion(
          se, intersect(scope$traits[[st]], rownames(se)), args$region)
      }
    }
    pattern <- if (length(spec$axes) > 1L) "composed" else spec$axes[[1L]]
    cell <- .lookupJointCell(pattern, dataForm)
    spArgs <- c(args, list(axes = spec$axes))
    # One call per spec with ALL methods: .runJointCell loops them per group so
    # the twas ensemble layer can combine a group's per-method fits.
    res <- .runJointCell(cell, pipeline, data, scope, jointMethods, spArgs)
    if (is.null(res)) next
    out <- if (is.null(out)) res
           else if (isFm) .rbindFineMappingResult(out, res, ldSketch = ldSketch)
           else .rbindTwasWeights(out, res, ldSketch = ldSketch)
  }
  out
}

# Individual-level (QtlDataset) input cannot joint over study: studies have
# disjoint samples (cross-study joints live on the sumstats slot). Preserve the
# historical axis-specific error messages.
.jointRejectStudyOnIndividual <- function(parsedJointSpec) {
  for (spec in parsedJointSpec) {
    if ("study" %in% spec$axes) {
      if (length(spec$axes) > 1L)
        stop("composed joint axes including 'study' require sumstats input.")
      stop("jointSpecification with axis 'study' requires sumstats input ",
           "(QtlDataset is a single individual-level study).")
    }
  }
}
