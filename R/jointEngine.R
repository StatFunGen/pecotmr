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

# A group's constant value on axis `ax` (study/context/trait), or NULL when the
# axis varies (jointed) so the prior lookup matches any value there.
# @noRd
.jointAxisVal <- function(ax, conditions) {
    u <- unique(as.character(conditions[[ax]]))
    if (length(u) > 1L) NULL else u[[1L]]
}

# The data-driven-prior LOOKUP key for a group's conditions: a varying (jointed)
# axis -> NULL (match-any, because the shared joint mr.mash fit lives on every
# per-context row), a constant axis -> its single value. Used only to find the
# mr.mash fit; the OUTPUT rows carry each condition's REAL (study, context,
# trait).
.jointPriorKey <- function(conditions) {
    list(
        study = .jointAxisVal("study", conditions),
        context = .jointAxisVal("context", conditions),
        trait = .jointAxisVal("trait", conditions)
    )
}

# The ";"-joined distinct members of a varying axis (the per-row provenance tag
# jointStudies/Contexts/Traits), or NA when the axis is constant.
.jointAxisMembers <- function(conditions, ax) {
    u <- unique(as.character(conditions[[ax]]))
    if (length(u) > 1L) str_flatten(u, ";") else NA_character_
}

# Slice a fine-mapping per-method CV payload (.fmSliceCv output:
# list(samplePartition, prediction = list(<m>_predicted = sample x condition),
# performance = list(<m>_performance = condition x 6))) down to one condition r,
# so each per-context FineMappingRow carries that context's CV.
.fmSliceCvCondition <- function(cv, r) {
    if (is.null(cv)) {
        return(NULL)
    }
    out <- list(samplePartition = cv$samplePartition)
    if (!is.null(cv$prediction)) {
        out$prediction <- map(cv$prediction, .fmCvSliceCol, r = r)
    }
    if (!is.null(cv$performance)) {
        out$performance <- map(cv$performance, .fmCvSliceRow, r = r)
    }
    out
}

# Slice a twas joint cvResult (.jointTwasCvResult output: list(samplePartition,
# predictions = sample x condition, metrics = condition x 6, foldFits)) to one
# condition r. The per-fold mr.mash fits span all conditions, so foldFits is
# shared unchanged.
.sliceTwasCvResultToCondition <- function(cvRes, r) {
    if (is.null(cvRes)) {
        return(NULL)
    }
    list(
        samplePartition = cvRes$samplePartition,
        predictions = if (!is.null(cvRes$predictions)) {
            cvRes$predictions[, r, drop = TRUE]
        } else {
            NULL
        },
        metrics = if (!is.null(cvRes$metrics)) {
            cvRes$metrics[r, , drop = TRUE]
        } else {
            NULL
        },
        foldFits = cvRes$foldFits
    )
}

# The engine assembles one per-row RECORD per fitted entry: a named list whose
# names match the target collection constructor's parameters (see
# .jointEntryRecords()). The driver collects these into a plain list and
# construct() folds them into a collection via .buildJointResult(). Nothing here
# enumerates columns, so adding a column needs only that it be put in the
# record.

# --- Trait-position and fine-mapping-region provenance ----------------------
# Two DISTINCT per-row anchors, kept separate because they mean different things
# and diverge by input type:
#   traitPos = the bare trait position, NO cis-window. QtlDataset -> the trait's
#              own phenotype rowRanges; QtlSumStats -> the SUPPLIED `traitPos`
#              column (NULL when absent -- the true position cannot be inferred
#              from summary statistics alone).
#   region   = the fine-mapping window. QtlDataset -> traitPos +/- cisWindow;
#              QtlSumStats -> the entry's variant span (sumstats carry no
#              cis-window, so the fitted span IS the region).
# Each returns a length-1 GRanges, or NULL when the anchor is unavailable (the
# accumulator / builder then records a chrUn sentinel).
.traitPosFor <- function(data, context, trait) {
    if (methods::is(data, "QtlDataset")) {
        se <- tryCatch(
            getPhenotypes(data, contexts = context),
            error = function(e) NULL
        )
        if (is.null(se)) {
            return(NULL)
        }
        rr <- SummarizedExperiment::rowRanges(se)
        if (!is_in(trait, names(rr))) {
            return(NULL)
        }
        return(GenomicRanges::granges(rr[trait])[1L])
    }
    if (methods::is(data, "QtlSumStats")) {
        if (!is_in("traitPos", names(data))) {
            return(NULL)
        }
        idx <- which(
            as.character(data$context) == context &
                as.character(data$trait) == trait
        )
        if (length(idx) == 0L) {
            return(NULL)
        }
        tp <- data$traitPos[idx[[1L]]]
        if (length(tp) == 0L) {
            return(NULL)
        }
        return(GenomicRanges::granges(tp)[1L])
    }
    NULL
}

.fitRegionFor <- function(data, context, trait, cisWindow = NULL) {
    if (methods::is(data, "QtlDataset")) {
        tp <- .traitPosFor(data, context, trait)
        if (is.null(tp)) {
            return(NULL)
        }
        w <- if (is.null(cisWindow)) 0L else as.integer(cisWindow)
        return(GenomicRanges::GRanges(
            as.character(GenomicRanges::seqnames(tp))[[1L]],
            IRanges::IRanges(
                max(1L, GenomicRanges::start(tp) - w),
                GenomicRanges::end(tp) + w
            )
        ))
    }
    if (methods::is(data, "QtlSumStats")) {
        idx <- which(
            as.character(data$context) == context &
                as.character(data$trait) == trait
        )
        if (length(idx) == 0L) {
            return(NULL)
        }
        gr <- .collectionEntry(data, idx[[1L]])
        if (length(gr) == 0L) {
            return(NULL)
        }
        return(GenomicRanges::GRanges(
            as.character(GenomicRanges::seqnames(gr))[[1L]],
            IRanges::IRanges(
                min(GenomicRanges::start(gr)),
                max(GenomicRanges::end(gr))
            )
        ))
    }
    NULL
}

# Vectorised per-row anchors for the pushRow builders: map the scalar helper
# over aligned (context, trait) vectors and assemble ONE GRanges of length n,
# using a chrUn sentinel where the anchor is unavailable. Built from plain
# chr/start/end vectors so mixed seqlevels never trip seqinfo merging.
.anchorVector <- function(
    data,
    contexts,
    traits,
    kind = c("traitPos", "region"),
    cisWindow = NULL
) {
    kind <- arg_match(kind)
    n <- length(traits)
    chrs <- rep("chrUn", n)
    starts <- rep(1L, n)
    ends <- rep(1L, n)
    anyFound <- FALSE
    for (i in seq_len(n)) {
        g <- if (kind == "traitPos") {
            .traitPosFor(data, contexts[[i]], traits[[i]])
        } else {
            .fitRegionFor(data, contexts[[i]], traits[[i]], cisWindow)
        }
        if (!is.null(g)) {
            anyFound <- TRUE
            chrs[i] <- as.character(GenomicRanges::seqnames(g))[[1L]]
            starts[i] <- GenomicRanges::start(g)[[1L]]
            ends[i] <- GenomicRanges::end(g)[[1L]]
        }
    }
    # Nothing resolved (e.g. a QtlSumStats with no supplied traitPos): return
    # NULL
    # so the builder omits the column entirely and getTraitPosition() reports
    # NA,
    # rather than a column full of chrUn sentinels.
    if (!anyFound) {
        return(NULL)
    }
    GenomicRanges::GRanges(
        chrs,
        IRanges::IRanges(start = starts, end = pmax(ends, starts))
    )
}

# ---- fitters (fitJointGroup) ------------------------------------------------

# (individual, fine-mapping) -> mvSuSiE joint fit + honest per-fold CV prior.
setMethod(
    "fitJointGroup",
    signature("IndividualJointGroup", "FmJointPipeline"),
    function(group, pipeline, token, args) {
        cfg <- .jpConfig(pipeline)
        Xc <- .jgX(group)
        Yc <- .jgY(group)
        nCond <- ncol(Yc)
        if (identical(token, "fsusie")) {
            return(.jointFitFsusie(group, Xc, Yc, nCond, cfg, args))
        }
        if (!identical(token, "mvsusie")) {
            msg <- glue(
                "fitJointGroup(IndividualJointGroup, FmJointPipeline): ",
                "unsupported token '{token}' ",
                "(expected 'mvsusie' or 'fsusie')."
            )
            abort(msg)
        }
        .jointFitMvsusie(group, Xc, Yc, nCond, cfg, args)
    }
)

# fsusie joint fit (functional SuSiE over the trait domain; individual-level,
# cross-trait). One per-condition entry per trait, with an optional CV slice.
# @noRd
.jointFitFsusie <- function(group, Xc, Yc, nCond, cfg, args) {
    if (length(.jgPos(group)) != nCond) {
        msg <- glue(
            "fitJointGroup: fsusie requires per-trait positions ('pos'); ",
            "it is cross-trait individual-level only."
        )
        abort(msg)
    }
    verbose <- if (is.null(cfg$verbose)) 1 else cfg$verbose
    fitArgs <- .fmMergeUserArgs(
        list(X = Xc, Y = Yc, pos = .jgPos(group)),
        "fsusie",
        args$methodArgs[["fsusie"]]
    )
    fit <- exec(fitFsusie, !!!fitArgs)
    # Collapse the functional fit to a variants x features weight matrix now
    # (trimming later drops fitted_wc/csd_X); store on $coef so a trimmed fit
    # can still yield TWAS weights.
    fit$coef <- tryCatch(
        fsusieWeights(fsusieFit = fit, variantIds = colnames(Xc)),
        error = function(e) NULL
    )
    fit <- .setFinemappingFitClass(fit, "fsusie")
    cvM <- .jointFsusieCv(Xc, Yc, group, cfg, args, verbose)
    map(
        seq_len(nCond),
        .jointFsusieEntry,
        fit = fit,
        cvM = cvM,
        Xc = Xc,
        cfg = cfg
    )
}

# Per-fold fsusie CV slice, or NULL when CV is disabled.
# @noRd
.jointFsusieCv <- function(Xc, Yc, group, cfg, args, verbose) {
    cvFolds <- if (is.null(cfg$cvFolds)) 0L else cfg$cvFolds
    if (cvFolds <= 1L) {
        return(NULL)
    }
    cv <- .fmWeightsCv(
        Xc,
        Yc,
        "fsusie",
        args$methodArgs,
        cvFolds,
        samplePartition = cfg$samplePartition,
        coverage = cfg$coverage,
        pos = .jgPos(group),
        verbose = verbose,
        numThreads = if (is.null(cfg$cvThreads)) 1L else cfg$cvThreads,
        seed = cfg$seed
    )
    .fmSliceCv(cv, "fsusie")
}

# One fsusie per-condition (trait) FineMappingRow, with its CV slice attached.
# @noRd
.jointFsusieEntry <- function(r, fit, cvM, Xc, cfg) {
    e <- .fmPostprocessOne(
        fit = fit,
        method = "fsusie",
        dataX = Xc,
        dataY = NULL,
        conditionIdx = r,
        coverage = cfg$coverage,
        secondaryCoverage = cfg$secondaryCoverage,
        signalCutoff = cfg$signalCutoff,
        minAbsCorr = cfg$minAbsCorr,
        csInput = "fsusie",
        fullFit = cfg$fullFit,
        fullFitAlphaOnly = cfg$fullFitAlphaOnly,
        includeAllCs = cfg$includeAllCs
    )
    if (!is.null(cvM)) {
        e <- .fmAttachCv(e, .fmSliceCvCondition(cvM, r))
    }
    e
}

# mvsusie joint fit: SER pre-screen the conditions, fit the survivors, then emit
# one per-context entry per ORIGINAL condition (NULL for screened-out columns).
# @noRd
.jointFitMvsusie <- function(group, Xc, Yc, nCond, cfg, args) {
    ddCut <- if (is.null(cfg$dataDrivenPriorWeightsCutoff)) {
        1e-10
    } else {
        cfg$dataDrivenPriorWeightsCutoff
    }
    verbose <- if (is.null(cfg$verbose)) 1 else cfg$verbose
    keep <- .jointMvSerScreen(Xc, Yc, nCond, args, verbose)
    if (is.null(keep)) {
        return(vector("list", nCond))
    }
    survivors <- which(keep)
    fitted <- .jointMvFit(
        group,
        Xc,
        Yc[, survivors, drop = FALSE],
        cfg,
        args,
        ddCut,
        verbose
    )
    map(
        seq_len(nCond),
        .jointMvEntry,
        fitted = fitted,
        Xc = Xc,
        keep = keep,
        survivors = survivors,
        cfg = cfg
    )
}

# SER pre-screen mask over the conditions: TRUE-all when no screen, the survivor
# mask when active, or NULL to signal < 2 survivors (skip the whole joint).
# @noRd
.jointMvSerScreen <- function(Xc, Yc, nCond, args, verbose) {
    keep <- rep(TRUE, nCond)
    if (!.fmScreenActive(args$pipCutoffToSkip)) {
        return(keep)
    }
    keep <- as.logical(.fmSerScreenColumns(Xc, Yc, args$pipCutoffToSkip))
    if (sum(keep) < 2L) {
        if (verbose >= 1) {
            msg <- glue(
                "Skipping mvsusie joint fit: < 2 of {nCond} conditions pass ",
                "the SER pre-screen."
            )
            inform(msg)
        }
        return(NULL)
    }
    if (sum(keep) < nCond && verbose >= 1) {
        msg <- glue(
            "mvsusie joint fit: SER pre-screen kept {sum(keep)} of ",
            "{nCond} conditions."
        )
        inform(msg)
    }
    keep
}

# Fit mvsusie over the surviving conditions with the data-driven reweighted
# prior, returning list(fit, cvM).
# @noRd
.jointMvFit <- function(group, Xc, Ys, cfg, args, ddCut, verbose) {
    key <- .jointPriorKey(.jgConditions(group))
    mvFitParts <- .fmLookupMrmashFit(
        args$twasWeights,
        key$study,
        key$trait,
        context = key$context
    )
    mvCv <- .fmLookupMrmashCv(
        args$twasWeights,
        key$study,
        key$trait,
        context = key$context
    )
    mvPrior <- .buildMvsusieReweightedPrior(mvFitParts, colnames(Ys), ddCut)
    mvBaseArgs <- list(
        X = Xc,
        Y = Ys,
        prior_variance = mvPrior$priorVariance,
        coverage = cfg$coverage
    )
    if (!is.null(mvPrior$residualVariance)) {
        mvBaseArgs$residual_variance <- mvPrior$residualVariance
    }
    fitArgs <- .fmMergeUserArgs(
        mvBaseArgs,
        "mvsusie",
        args$methodArgs[["mvsusie"]]
    )
    fit <- exec(fitMvsusie, !!!fitArgs)
    fit <- .setFinemappingFitClass(fit, "mvsusie")
    cvM <- .jointMvCv(
        Xc,
        Ys,
        cfg,
        args,
        mvPrior,
        mvFitParts,
        mvCv,
        ddCut,
        verbose
    )
    list(fit = fit, cvM = cvM)
}

# Per-fold mvsusie CV slice (reusing the mr.mash prior per fold), or NULL.
# @noRd
.jointMvCv <- function(
    Xc,
    Ys,
    cfg,
    args,
    mvPrior,
    mvFitParts,
    mvCv,
    ddCut,
    verbose
) {
    cvFolds <- if (is.null(cfg$cvFolds)) 0L else cfg$cvFolds
    if (cvFolds <= 1L) {
        return(NULL)
    }
    sp <- cfg$samplePartition
    if (is.null(sp) && !is.null(mvCv)) {
        sp <- mvCv$samplePartition
    }
    mvPriorCv <- .fmBuildMvsusiePriorCv(mvCv, mvFitParts, colnames(Ys), ddCut)
    cv <- .fmWeightsCv(
        Xc,
        Ys,
        "mvsusie",
        args$methodArgs,
        cvFolds,
        samplePartition = sp,
        coverage = cfg$coverage,
        verbose = verbose,
        mvPrior = mvPrior,
        mvPriorCv = mvPriorCv,
        numThreads = if (is.null(cfg$cvThreads)) 1L else cfg$cvThreads,
        seed = cfg$seed
    )
    .fmSliceCv(cv, "mvsusie")
}

# One mvsusie per-condition entry (NULL for a screened-out column), sliced at
# the condition's position in the fitted survivor set + its CV slice.
# @noRd
.jointMvEntry <- function(i, fitted, Xc, keep, survivors, cfg) {
    if (!keep[i]) {
        return(NULL)
    }
    r <- match(i, survivors)
    e <- .fmPostprocessOne(
        fit = fitted$fit,
        method = "mvsusie",
        dataX = Xc,
        dataY = NULL,
        conditionIdx = r,
        coverage = cfg$coverage,
        secondaryCoverage = cfg$secondaryCoverage,
        signalCutoff = cfg$signalCutoff,
        minAbsCorr = cfg$minAbsCorr,
        csInput = "X",
        fullFit = cfg$fullFit,
        fullFitAlphaOnly = cfg$fullFitAlphaOnly,
        includeAllCs = cfg$includeAllCs
    )
    if (!is.null(fitted$cvM)) {
        e <- .fmAttachCv(e, .fmSliceCvCondition(fitted$cvM, r))
    }
    e
}

# (sumstats, fine-mapping) -> mvSuSiE-rss joint fit. RSS has no sample folds and
# no fsusie variant.
setMethod(
    "fitJointGroup",
    signature("SumStatsJointGroup", "FmJointPipeline"),
    function(group, pipeline, token, args) {
        if (identical(token, "fsusie")) {
            abort(
                "fsusie has no RSS variant; it requires individual-level input."
            )
        }
        if (!identical(token, "mvsusie")) {
            msg <- glue(
                "fitJointGroup(SumStatsJointGroup, FmJointPipeline): ",
                "unsupported token '{token}' (expected 'mvsusie')."
            )
            abort(msg)
        }
        cfg <- .jpConfig(pipeline)
        fit <- .jointFitMvsusieRss(group, cfg, args)
        # One per-condition entry (RSS has no sample folds).
        map(
            seq_len(ncol(.jgZ(group))),
            .jointRssEntry,
            fit = fit,
            group = group,
            cfg = cfg
        )
    }
)

# mvSuSiE-RSS joint fit over summary statistics, with the data-driven reweighted
# prior. Returns the class-tagged fit.
# @noRd
.jointFitMvsusieRss <- function(group, cfg, args) {
    ddCut <- if (is.null(cfg$dataDrivenPriorWeightsCutoff)) {
        1e-10
    } else {
        cfg$dataDrivenPriorWeightsCutoff
    }
    key <- .jointPriorKey(.jgConditions(group))
    mvFitParts <- .fmLookupMrmashFit(
        args$twasWeights,
        key$study,
        key$trait,
        context = key$context
    )
    mvPrior <- .buildMvsusieReweightedPrior(
        mvFitParts,
        colnames(.jgZ(group)),
        ddCut
    )
    mvBaseArgs <- list(
        Z = .jgZ(group),
        R = .jgR(group),
        N = as.numeric(stats::median(.jgN(group))),
        prior_variance = mvPrior$priorVariance,
        coverage = cfg$coverage
    )
    if (!is.null(mvPrior$residualVariance)) {
        mvBaseArgs$residual_variance <- mvPrior$residualVariance
    }
    fitArgs <- .fmMergeUserArgs(
        mvBaseArgs,
        "mvsusie",
        args$methodArgs[["mvsusie"]]
    )
    fit <- exec(fitMvsusieRss, !!!fitArgs)
    .setFinemappingFitClass(fit, "mvsusie")
}

# One RSS per-condition FineMappingRow (csInput = "Xcorr").
# @noRd
.jointRssEntry <- function(r, fit, group, cfg) {
    .fmPostprocessOne(
        fit = fit,
        method = "mvsusie",
        dataX = .jgR(group),
        dataY = NULL,
        conditionIdx = r,
        coverage = cfg$coverage,
        secondaryCoverage = cfg$secondaryCoverage,
        signalCutoff = cfg$signalCutoff,
        minAbsCorr = cfg$minAbsCorr,
        csInput = "Xcorr",
        fullFit = cfg$fullFit,
        fullFitAlphaOnly = cfg$fullFitAlphaOnly,
        includeAllCs = cfg$includeAllCs
    )
}

# Select the list element whose name (stripped of a _predicted/_performance
# suffix) matches `token`; NULL if none.
# @noRd
.jointPickByBase <- function(lst, token) {
    if (is.null(lst) || length(lst) == 0L) {
        return(NULL)
    }
    bare <- str_remove(
        names(lst),
        "(_predicted|Predicted|_performance|Performance)$"
    )
    hit <- which(bare == token)
    if (length(hit) == 0L) NULL else lst[[hit[[1L]]]]
}

# Reshape a twasWeightsCv() result into the single joint entry's cvResult: the
# out-of-fold prediction matrix, the per-condition metric rows, and the per-fold
# mr.mash fits (named fold_<j>) that fineMappingPipeline's mvSuSiE path
# consumes.
.jointTwasCvResult <- function(cv, token) {
    if (is.null(cv)) {
        return(NULL)
    }
    ffKey <- str_c(token, "_weights")
    foldFits <- if (!is.null(cv$foldFits)) {
        ff <- map(cv$foldFits, ffKey)
        if (all(map_lgl(ff, is.null))) NULL else ff
    } else {
        NULL
    }
    list(
        samplePartition = cv$samplePartition,
        predictions = .jointPickByBase(cv$prediction, token),
        metrics = .jointPickByBase(cv$performance, token),
        foldFits = foldFits
    )
}

# learnTwasWeights key for a bare token (fine-mapping tokens key differently,
# e.g. susieInf -> susie_inf_weights).
.twasMethodKey <- function(token) {
    ad <- .twasFineMappingMethodAdapters[[token]]
    if (!is.null(ad)) ad$methodKey else str_c(token, "_weights")
}

# Fine-mapping CV handoff for one twas method: extract that method's out-of-fold
# predictions + performance from fineMappingPipeline's retained CV (shared fold
# partition), shaped like .jointTwasCvResult so the per-condition slice reuses
# it instead of re-cross-validating an FM-derived method (susie / mvsusie /
# ...).
.twasFmHandoffCv <- function(fineMappingCv, token) {
    if (is.null(fineMappingCv) || is.null(fineMappingCv$prediction)) {
        return(NULL)
    }
    base <- str_remove(
        names(fineMappingCv$prediction),
        "(_predicted|Predicted)$"
    )
    hit <- which(base == token)
    if (length(hit) == 0L) {
        return(NULL)
    }
    pBase <- str_remove(
        names(fineMappingCv$performance),
        "(_performance|Performance)$"
    )
    pHit <- which(pBase == token)
    list(
        samplePartition = fineMappingCv$samplePartition,
        predictions = fineMappingCv$prediction[[hit[[1L]]]],
        metrics = if (length(pHit)) {
            fineMappingCv$performance[[pHit[[1L]]]]
        } else {
            NULL
        },
        foldFits = NULL
    )
}

# (individual, twas) -> ONE weight method fit over the group's conditions, as
# per-condition entries (sliced from the variants x conditions weight matrix),
# each with its full-data weights + retained fit + per-condition CV slice. This
# is the SHARED per-method twas fitting (one method per call, like the FM
# fitters); the SR-TWAS ensemble combines methods in a layer above (see
# .twasEnsembleLayer). Owns the orchestration formerly in
# .twasWeightsPipelineMatrix: FM-fit injection (FM-derived tokens extract from
# the precomputed fit), the FM CV handoff (reuse fine-mapping's own CV),
# spike-and-slab pi from an internal mr.ash fit, CV knobs, and fitFullData =
# FALSE (CV-only) entries.
setMethod(
    "fitJointGroup",
    signature("IndividualJointGroup", "TwasJointPipeline"),
    function(group, pipeline, token, args) {
        cfg <- .jpConfig(pipeline)
        Xc <- .jgX(group)
        Yc <- .jgY(group)
        nCond <- ncol(Yc)
        cond <- .jgConditions(group)
        methodKey <- .twasMethodKey(token)
        stdz <- cfg$standardized
        fittedModels <- args$fittedModels %||% list()
        ma <- .jointTwasMethodArgs(
            args,
            methodKey,
            token,
            fittedModels,
            Xc,
            Yc,
            cond,
            cfg,
            stdz
        )
        wm <- set_names(list(ma), methodKey)
        full <- .jointTwasFitFull(
            Xc,
            Yc,
            cond,
            wm,
            fittedModels,
            cfg,
            stdz,
            nCond
        )
        cvRes <- .jointTwasCv(Xc, Yc, wm, ma, full$W, args, cfg, token)
        # One per-condition entry: that condition's weight column + the shared
        # fit + its CV slice. fitFullData = FALSE -> CV-only entry.
        map(
            seq_len(nCond),
            .jointTwasEntry,
            full = full,
            cvRes = cvRes,
            stdz = stdz,
            cfg = cfg
        )
    }
)

# Resolve the method args for a TWAS token: prefer the unified methodList, else
# the explicit-jointSpec methodArgs; inject an FM-derived fit; compute a
# spike-and-slab pi for bayes_c / bayes_b when estimatePi is on.
# @noRd
.jointTwasMethodArgs <- function(
    args,
    methodKey,
    token,
    fittedModels,
    Xc,
    Yc,
    cond,
    cfg,
    stdz
) {
    ma <- if (
        !is.null(args$methodList) && is_in(methodKey, names(args$methodList))
    ) {
        args$methodList[[methodKey]]
    } else if (!is.null(args$methodArgs)) {
        args$methodArgs[[methodKey]]
    } else {
        NULL
    }
    if (is.null(ma)) {
        ma <- list()
    }
    # FM-fit injection: an FM-derived token extracts its weights from the
    # precomputed fine-mapping fit rather than refitting.
    adapter <- .twasFineMappingMethodAdapters[[token]]
    if (
        !is.null(adapter) &&
            !is.null(fittedModels[[token]]) &&
            is.null(ma[[adapter$fitArg]])
    ) {
        ma[[adapter$fitArg]] <- fittedModels[[token]]
    }
    if (isTRUE(cfg$estimatePi) && is_in(token, c("bayes_c", "bayes_b"))) {
        ma <- .jointTwasSpikeSlabPi(ma, token, Xc, Yc, cond, cfg, stdz)
    }
    ma
}

# Spike-and-slab pi from an internal mr.ash fit (self-contained per method).
# @noRd
.jointTwasSpikeSlabPi <- function(ma, token, Xc, Yc, cond, cfg, stdz) {
    mrA <- learnTwasWeights(
        Xc,
        Yc,
        weightMethods = list(mrash_weights = list()),
        study = as.character(cond$study[1L]),
        context = as.character(cond$context[1L]),
        trait = as.character(cond$trait[1L]),
        retainFits = TRUE,
        standardized = stdz,
        dataType = cfg$dataType,
        verbose = 0,
        seed = cfg$seed
    )
    piHat <- as.numeric(estimateSparsity(mrA))
    if (token == "bayes_c" && is.null(ma$pi)) {
        ma$pi <- piHat
    }
    if (token == "bayes_b" && is.null(ma$probIn)) {
        ma$probIn <- piHat
    }
    ma
}

# Full-data TWAS weight fit for a joint group. Returns list(W, fitParts, vids);
# W is NULL (a CV-only run) when fitFullData is FALSE.
# @noRd
.jointTwasFitFull <- function(
    Xc,
    Yc,
    cond,
    wm,
    fittedModels,
    cfg,
    stdz,
    nCond
) {
    fitFullData <- if (is.null(cfg$fitFullData)) {
        TRUE
    } else {
        isTRUE(cfg$fitFullData)
    }
    if (!fitFullData) {
        return(list(W = NULL, fitParts = NULL, vids = colnames(Xc)))
    }
    rfd <- if (is.null(cfg$retainFitDetail)) "slim" else cfg$retainFitDetail
    verbose <- if (is.null(cfg$verbose)) 1 else cfg$verbose
    tw <- learnTwasWeights(
        Xc,
        Yc,
        weightMethods = wm,
        study = as.character(cond$study[1L]),
        context = as.character(cond$context[1L]),
        trait = as.character(cond$trait[1L]),
        fittedModels = fittedModels,
        retainFits = TRUE,
        retainFitDetail = rfd,
        standardized = stdz,
        dataType = cfg$dataType,
        verbose = verbose,
        seed = cfg$seed
    )
    base <- .twrRowParts(tw, 1L)
    vids <- .twrPartsVariantIds(base)
    W <- getWeights(base)
    if (!is.matrix(W)) {
        W <- matrix(W, ncol = nCond, dimnames = list(vids, NULL))
    }
    list(W = W, fitParts = getFits(base), vids = vids)
}

# Warn when cross-validating mr.mash with a single full-data data-driven prior
# reused across folds (information leakage; supply per-fold priors instead).
# @noRd
.jointTwasLeakageWarn <- function(args, ma) {
    if (
        is.null(args$dataDrivenPriorMatricesCv) &&
            !is.null(ma$dataDrivenPriorMatrices)
    ) {
        msg <- glue(
            "Cross-validating mr.mash with a single data-driven prior ",
            "computed on the full data: the same prior is reused for every ",
            "fold, so each fold's prior was informed by its own held-out ",
            "samples (information leakage). Supply per-fold priors via ",
            "dataDrivenPriorMatricesCv (--mixture-prior-cv) for honest ",
            "cross-validation."
        )
        warn(msg)
    }
}

# Cross-validated prediction result for a TWAS token: reuse fine-mapping's own
# CV when available, else run twasWeightsCv (skipping all-zero-weight methods).
# @noRd
.jointTwasCv <- function(Xc, Yc, wm, ma, W, args, cfg, token) {
    cvFolds <- if (is.null(cfg$cvFolds)) 0L else cfg$cvFolds
    if (cvFolds <= 1L) {
        return(NULL)
    }
    cvRes <- .twasFmHandoffCv(args$fineMappingCv, token)
    if (!is.null(cvRes) || (!is.null(W) && all(W == 0))) {
        return(cvRes)
    }
    .jointTwasLeakageWarn(args, ma)
    verbose <- if (is.null(cfg$verbose)) 1 else cfg$verbose
    sp <- if (!is.null(args$samplePartition)) {
        args$samplePartition
    } else {
        cfg$samplePartition
    }
    mcv <- if (is.null(cfg$maxCvVariants) || cfg$maxCvVariants <= 0) {
        Inf
    } else {
        cfg$maxCvVariants
    }
    cv <- twasWeightsCv(
        Xc,
        Yc,
        fold = cvFolds,
        samplePartitions = sp,
        weightMethods = wm,
        retainFits = TRUE,
        maxNumVariants = mcv,
        numThreads = if (is.null(cfg$cvThreads)) 1 else cfg$cvThreads,
        data_driven_priorMatricesCv = args$dataDrivenPriorMatricesCv,
        verbose = verbose,
        seed = cfg$seed
    )
    .jointTwasCvResult(cv, token)
}

# One per-condition TwasWeightsRow (empty when there are no full-data weights,
# else that condition's weight column + shared fit + CV slice).
# @noRd
.jointTwasEntry <- function(r, full, cvRes, stdz, cfg) {
    cvR <- if (!is.null(cvRes)) {
        .sliceTwasCvResultToCondition(cvRes, r)
    } else {
        NULL
    }
    if (is.null(full$W)) {
        return(twasWeightsRow(
            variantIds = character(0),
            weights = NULL,
            cvResult = cvR,
            standardized = stdz,
            dataType = cfg$dataType
        ))
    }
    twasWeightsRow(
        variantIds = full$vids,
        weights = full$W[, r],
        fits = full$fitParts,
        cvResult = cvR,
        standardized = stdz,
        dataType = cfg$dataType
    )
}

# (sumstats, twas) -> mr.mash-rss joint fit as ONE matrix entry. No sample
# folds.
setMethod(
    "fitJointGroup",
    signature("SumStatsJointGroup", "TwasJointPipeline"),
    function(group, pipeline, token, args) {
        cfg <- .jpConfig(pipeline)
        rfd <- if (is.null(cfg$retainFitDetail)) "slim" else cfg$retainFitDetail
        weights <- mrmashRssWeights(
            stat = list(z = .jgZ(group), n = .jgN(group)),
            LD = .jgR(group),
            retainFit = TRUE,
            fitDetail = rfd
        )
        vids <- rownames(weights)
        if (is.null(vids)) {
            vids <- rownames(.jgZ(group))
        }
        fitParts <- attr(weights, "fit")
        if (!is.matrix(weights)) {
            weights <- matrix(
                weights,
                ncol = ncol(.jgZ(group)),
                dimnames = list(vids, NULL)
            )
        }
        # One per-condition entry: that condition's weight column + the shared
        # fit.
        map(
            seq_len(ncol(weights)),
            .jointColEntry,
            vids = vids,
            weights = weights,
            fitParts = fitParts,
            cfg = cfg
        )
    }
)

# ---- result construction (construct) ----------------------------------------

# Both pipelines assemble identically-shaped joint rows; only the result
# collection differs (the axis-3 divergence the markers encode). Only the joint*
# columns for axes that actually vary are attached.
# Fold the accumulator's per-row records into one collection: build each record
# into a 1-row collection via `constructor`, then union them with the generic
# .rbindCollections (which aligns optional columns and sets the ldSketch slot).
# The only per-row transforms are: wrap `entry` in a list, and drop an NA joint*
# value so the 1-row collection omits that column (the union re-adds it, padded,
# only when some row is joint). Every other field -- study/context/trait/method,
# region, and any future column -- flows through untouched.
.buildJointResult <- function(constructor, records, ldSketch = NULL) {
    if (length(records) == 0L) {
        return(NULL)
    }
    parts <- map(records, .jointBuildRecordPart, constructor = constructor)
    .rbindCollections(parts, ldSketch = ldSketch)
}

setMethod("construct", "FmJointPipeline", function(pipeline, records, ...) {
    .buildJointResult(
        QtlFineMappingResult,
        records,
        .jpConfig(pipeline)$ldSketch
    )
})

setMethod("construct", "TwasJointPipeline", function(pipeline, records, ...) {
    .buildJointResult(TwasWeights, records, .jpConfig(pipeline)$ldSketch)
})

# ---- enumerators (pattern x dataForm -> list<JointGroup>) --------------------

# cross-context / individual: one group per scoped trait present in >= 2 scoped
# contexts (the conditions are those (study, context, trait) rows).
.enumCrossContextIndividual <- function(data, scope, args = list()) {
    study <- getStudy(data)
    if (!is_in(study, scope$studies)) {
        return(list())
    }
    scopedContexts <- scope$contexts[[study]]
    scopedTraits <- scope$traits[[study]]
    if (length(scopedContexts) < 2L) {
        return(list())
    }
    verbose <- if (is.null(args$verbose)) 1 else args$verbose
    groups <- list()
    for (tid in scopedTraits) {
        xy <- .buildIndividualCrossContextXy(
            data,
            tid,
            scopedContexts,
            args$cisWindow,
            verbose,
            label = "jointCrossContext",
            region = args$region
        )
        if (is.null(xy)) {
            next
        }
        groups[[length(groups) + 1L]] <- new(
            "IndividualJointGroup",
            conditions = tibble(
                study = study,
                context = xy$perTraitContexts,
                trait = tid
            ),
            X = xy$X,
            Y = xy$Y
        )
    }
    groups
}

# cross-context / sumstats.
.enumCrossContextSumstats <- function(data, scope, args = list()) {
    ldSketch <- getLdSketch(data)
    studyCol <- as.character(data$study)
    contextCol <- as.character(data$context)
    traitCol <- as.character(data$trait)
    groups <- list()
    for (s in scope$studies) {
        scopedContexts <- scope$contexts[[s]]
        scopedTraits <- scope$traits[[s]]
        if (length(scopedContexts) < 2L) {
            next
        }
        for (tid in scopedTraits) {
            tupleRows <- which(
                studyCol == s &
                    traitCol == tid &
                    is_in(contextCol, scopedContexts)
            )
            if (length(tupleRows) < 2L) {
                next
            }
            ctxNames <- contextCol[tupleRows]
            jz <- .buildJointSumstatZMatrix(
                data,
                tupleRows,
                ctxNames,
                errorLabel = "jointCrossContext (QtlSumStats)",
                ldSketch = ldSketch,
                cutoffs = args$cutoffs
            )
            ldMat <- .fmLdFromSketch(ldSketch, jz$variantIds)
            groups[[length(groups) + 1L]] <- new(
                "SumStatsJointGroup",
                conditions = tibble(
                    study = s,
                    context = ctxNames,
                    trait = tid
                ),
                Z = jz$Z,
                R = ldMat,
                N = jz$nVec
            )
        }
    }
    groups
}

# cross-trait / individual: one group per scoped context with >= 2 scoped
# traits.
.enumCrossTraitIndividual <- function(data, scope, args = list()) {
    study <- getStudy(data)
    if (!is_in(study, scope$studies)) {
        return(list())
    }
    scopedContexts <- scope$contexts[[study]]
    scopedTraits <- scope$traits[[study]]
    verbose <- if (is.null(args$verbose)) 1 else args$verbose
    groups <- list()
    for (cx in scopedContexts) {
        xy <- .buildIndividualCrossTraitXy(
            data,
            cx,
            scopedTraits,
            args$cisWindow,
            verbose,
            label = "jointCrossTrait",
            study = study,
            region = args$region
        )
        if (is.null(xy)) {
            next
        }
        # Functional positions (one per trait column) for fsusie's domain;
        # mvsusie ignores them. Matches the trait order of Y.
        rr <- SummarizedExperiment::rowRanges(xy$se)
        rr <- rr[match(colnames(xy$Y), rownames(xy$se))]
        pos <- (GenomicRanges::start(rr) + GenomicRanges::end(rr)) / 2
        groups[[length(groups) + 1L]] <- new(
            "IndividualJointGroup",
            conditions = tibble(
                study = study,
                context = cx,
                trait = xy$traitsHere
            ),
            X = xy$X,
            Y = xy$Y,
            pos = as.numeric(pos)
        )
    }
    groups
}

# cross-trait / sumstats.
.enumCrossTraitSumstats <- function(data, scope, args = list()) {
    ldSketch <- getLdSketch(data)
    studyCol <- as.character(data$study)
    contextCol <- as.character(data$context)
    traitCol <- as.character(data$trait)
    groups <- list()
    for (s in scope$studies) {
        scopedContexts <- scope$contexts[[s]]
        scopedTraits <- scope$traits[[s]]
        for (cx in scopedContexts) {
            tupleRows <- which(
                studyCol == s & contextCol == cx & is_in(traitCol, scopedTraits)
            )
            if (length(tupleRows) < 2L) {
                next
            }
            trNames <- traitCol[tupleRows]
            jz <- .buildJointSumstatZMatrix(
                data,
                tupleRows,
                trNames,
                errorLabel = "jointCrossTrait (QtlSumStats)",
                ldSketch = ldSketch,
                cutoffs = args$cutoffs
            )
            ldMat <- .fmLdFromSketch(ldSketch, jz$variantIds)
            groups[[length(groups) + 1L]] <- new(
                "SumStatsJointGroup",
                conditions = tibble(
                    study = s,
                    context = cx,
                    trait = trNames
                ),
                Z = jz$Z,
                R = ldMat,
                N = jz$nVec
            )
        }
    }
    groups
}

# cross-study / sumstats (no individual form: individual-level studies have
# disjoint samples). One group per (context, trait) present in >= 2 scoped
# studies; the study axis varies -> "joint" + jointStudies.
.enumCrossStudySumstats <- function(data, scope, args = list()) {
    ldSketch <- getLdSketch(data)
    cols <- list(
        study = as.character(data$study),
        context = as.character(data$context),
        trait = as.character(data$trait)
    )
    allCtxs <- unique(unlist(scope$contexts, use.names = FALSE))
    allTrs <- unique(unlist(scope$traits, use.names = FALSE))
    groups <- list()
    for (cx in allCtxs) {
        for (tid in allTrs) {
            g <- .enumCrossStudyGroup(
                data,
                scope,
                args,
                cols,
                cx,
                tid,
                ldSketch
            )
            if (!is.null(g)) {
                groups[[length(groups) + 1L]] <- g
            }
        }
    }
    groups
}

# One cross-study joint group for a (context, trait) cell, or NULL when fewer
# than two studies survive the scope filter -- a "joint" fit over one study is
# just that study, so the cell contributes nothing.
# @noRd
.enumCrossStudyGroup <- function(data, scope, args, cols, cx, tid, ldSketch) {
    tupleRows <- which(
        cols$context == cx &
            cols$trait == tid &
            is_in(cols$study, scope$studies)
    )
    keep <- map_lgl(
        tupleRows,
        .jointTupleRowInScope,
        studyCol = cols$study,
        cx = cx,
        tid = tid,
        scope = scope
    )
    tupleRows <- tupleRows[keep]
    if (length(tupleRows) < 2L) {
        return(NULL)
    }
    stNames <- cols$study[tupleRows]
    jz <- .buildJointSumstatZMatrix(
        data,
        tupleRows,
        stNames,
        errorLabel = "jointCrossStudy",
        ldSketch = ldSketch,
        cutoffs = args$cutoffs
    )
    new(
        "SumStatsJointGroup",
        conditions = tibble(
            study = stNames,
            context = cx,
            trait = tid
        ),
        Z = jz$Z,
        R = .fmLdFromSketch(ldSketch, jz$variantIds),
        N = jz$nVec
    )
}

# composed / individual: ONE group joining every scoped (context, trait) tuple
# for the study. Both context and trait vary across conditions, so both collapse
# to "joint" (the conditions model handles multi-varying-axis uniformly; if the
# tuples happen to share a context it degrades to cross-trait, and vice versa).
.enumComposedIndividual <- function(data, scope, args = list()) {
    study <- getStudy(data)
    if (!is_in(study, scope$studies)) {
        return(list())
    }
    verbose <- if (is.null(args$verbose)) 1 else args$verbose
    xy <- .buildComposedIndividualXy(
        data,
        scope,
        study,
        args$cisWindow,
        verbose,
        label = "composed",
        region = args$region
    )
    if (is.null(xy)) {
        return(list())
    }
    # Conditions follow the fitted Y columns ("context:trait"), so dropped
    # tuples don't desync conditions from Y. Split on the first ":" (contexts
    # are simple labels; trait ids may themselves contain ":").
    labs <- colnames(xy$Y)
    conds <- tibble(
        study = study,
        context = str_remove(labs, ":.*$"),
        trait = str_remove(labs, "^[^:]*:")
    )
    list(new("IndividualJointGroup", conditions = conds, X = xy$X, Y = xy$Y))
}

# univariate / individual: one 1-condition group per (study, context, trait) in
# scope -- the per-(context, trait) iteration expressed as engine groups, so
# univariate methods (lasso / enet / susie / ...) flow through the SAME per-
# method fitter + ensemble layer as the joint ones (minGroup = 1).
.enumUnivariateIndividual <- function(data, scope, args = list()) {
    study <- getStudy(data)
    if (!is_in(study, scope$studies)) {
        return(list())
    }
    naAction <- if (is.null(args$naAction)) "drop" else args$naAction
    groups <- list()
    for (cx in scope$contexts[[study]]) {
        se <- getPhenotypes(data, contexts = cx)
        for (tid in intersect(scope$traits[[study]], rownames(se))) {
            Y <- .fmResidPheno(
                data,
                contexts = cx,
                traitId = tid,
                naAction = naAction
            )
            X <- if (is.null(args$region)) {
                .fmResidGeno(
                    data,
                    contexts = cx,
                    traitId = tid,
                    cisWindow = args$cisWindow
                )
            } else {
                .fmResidGeno(data, contexts = cx, region = args$region)
            }
            common <- intersect(rownames(X), rownames(Y))
            if (length(common) < 2L) {
                next
            }
            groups[[length(groups) + 1L]] <- new(
                "IndividualJointGroup",
                conditions = tibble(
                    study = study,
                    context = cx,
                    trait = tid
                ),
                X = X[common, , drop = FALSE],
                Y = Y[common, , drop = FALSE]
            )
        }
    }
    groups
}

# composed / sumstats: general N-axis joint. `args$axes` (subset of study /
# context / trait) names the collapsed axes; rows split by the complement
# (fixed) axes form one group each. Reuses .enumerateComposedSumstatGroups.
.enumComposedSumstats <- function(data, scope, args = list()) {
    axes <- args$axes
    if (is.null(axes)) {
        axes <- c("context", "trait")
    }
    ldSketch <- getLdSketch(data)
    gi <- .enumerateComposedSumstatGroups(list(axes = axes), data, scope)
    if (is.null(gi)) {
        return(list())
    }
    groups <- list()
    for (gIdx in gi$groups) {
        if (length(gIdx) < 2L) {
            next
        }
        colLabels <- map_chr(gIdx, .jointGroupColLabel, gi = gi)
        jz <- .buildJointSumstatZMatrix(
            data,
            gIdx,
            colLabels,
            errorLabel = "composed (QtlSumStats)",
            ldSketch = ldSketch,
            cutoffs = args$cutoffs
        )
        ldMat <- .fmLdFromSketch(ldSketch, jz$variantIds)
        groups[[length(groups) + 1L]] <- new(
            "SumStatsJointGroup",
            conditions = tibble(
                study = gi$studyCol[gIdx],
                context = gi$contextCol[gIdx],
                trait = gi$traitCol[gIdx]
            ),
            Z = jz$Z,
            R = ldMat,
            N = jz$nVec
        )
    }
    groups
}

# ---- engine -----------------------------------------------------------------

# Twas per-group args: resolve the group's fine-mapping fits + CV (keyed on its
# first condition -- the joint fit is shared across conditions) and fix ONE
# shared fold partition (so every method's out-of-fold CV predictions align for
# the ensemble layer). Returns `args` unchanged for fine-mapping pipelines.
.twasGroupArgs <- function(g, pipeline, args) {
    if (!is(pipeline, "TwasJointPipeline")) {
        return(args)
    }
    cfg <- .jpConfig(pipeline)
    cond <- .jgConditions(g)
    out <- args
    fmRes <- args$fineMappingResult
    if (!is.null(fmRes)) {
        s1 <- as.character(cond$study[[1L]])
        c1 <- as.character(cond$context[[1L]])
        t1 <- as.character(cond$trait[[1L]])
        nR <- if (is.null(args$nRegions)) 1L else args$nRegions
        bi <- if (is.null(args$regionIndex)) 1L else args$regionIndex
        af <- .twasFineMappingFits(fmRes, study = s1, context = c1, trait = t1)
        out$fittedModels <- if (is.null(af)) {
            list()
        } else {
            .twasFitsForRegion(af, bi, nR)
        }
        out$fineMappingCv <- .twasCvResultFor(fmRes, s1, c1, t1)
    }
    cvF <- if (is.null(cfg$cvFolds)) 0L else cfg$cvFolds
    if (cvF > 1L && is(g, "IndividualJointGroup")) {
        sp <- args$samplePartition
        if (is.null(sp)) {
            sp <- cfg$samplePartition
        }
        if (is.null(sp) && !is.null(out$fineMappingCv)) {
            sp <- out$fineMappingCv$samplePartition
        }
        if (is.null(sp)) {
            sp <- .normalizeCvFolds(
                cvF,
                NULL,
                rownames(.jgX(g))
            )$samplePartition
        }
        out$samplePartition <- sp
    }
    out
}

# Append one output record per (condition, method) to the joint-rows accumulator
# `rows` (mutated by reference), resolving each condition's fine-mapping region
# + trait position. `grp` bundles the per-group invariants list(cond, js, jc,
# jt).
# @noRd
.jointEntryRecords <- function(entries, method, grp, data, cisWindow) {
    cond <- grp$cond
    recs <- map(
        seq_len(min(length(entries), nrow(cond))),
        .jointEntryRecordAt,
        entries = entries,
        cond = cond,
        method = method,
        grp = grp,
        data = data,
        cisWindow = cisWindow
    )
    compact(recs)
}

# Run one dispatch cell: enumerate joint groups, fit each method (S4 dispatch on
# the group x pipeline pair) per group, accumulate per-context rows, build the
# per-pipeline result. The loop is GROUP-outer / token-inner so the twas
# ensemble layer can combine a group's per-method fits in place (FM is
# unaffected by the loop order). Per-method fitting is identical for FM and twas
# -- one method -> per-condition entries; the SR-TWAS ensemble is a layer ON TOP
# of that.
.runJointCell <- function(cell, pipeline, data, scope, tokens, args = list()) {
    groups <- .jcEnumerate(cell)(data, scope, args)
    groups <- keep(groups, .jointGroupMeetsMin, minGroup = .jcMinGroup(cell))
    if (length(groups) == 0L) {
        return(NULL)
    }
    doEnsemble <- is(pipeline, "TwasJointPipeline") &&
        isTRUE(.jpConfig(pipeline)$ensemble)
    records <- list_flatten(map(
        groups,
        .jointCellGroup,
        pipeline = pipeline,
        data = data,
        tokens = tokens,
        args = args,
        doEnsemble = doEnsemble
    ))
    construct(pipeline, records)
}

# Fit every token for one joint group -> records, plus the optional SR-TWAS
# ensemble layer (built ON TOP of the shared per-method fits, when >= 2 methods
# produced entries).
# @noRd
.jointCellGroup <- function(g, pipeline, data, tokens, args, doEnsemble) {
    cond <- .jgConditions(g)
    # Provenance: the ";"-joined members of each varying axis (identical on
    # every per-context row of this group).
    grp <- list(
        cond = cond,
        js = .jointAxisMembers(cond, "study"),
        jc = .jointAxisMembers(cond, "context"),
        jt = .jointAxisMembers(cond, "trait")
    )
    fitArgs <- .twasGroupArgs(g, pipeline, args)
    records <- list()
    perTokenEntries <- list()
    for (token in tokens) {
        entries <- .jointGroupTokenEntries(
            g,
            pipeline,
            token,
            fitArgs,
            cond,
            args
        )
        if (is.null(entries) || length(entries) == 0L) {
            next
        }
        perTokenEntries[[token]] <- entries
        records <- c(
            records,
            .jointEntryRecords(entries, token, grp, data, args$cisWindow)
        )
    }
    if (doEnsemble && length(perTokenEntries) >= 2L) {
        records <- c(
            records,
            .jointEntryRecords(
                .twasEnsembleLayer(g, perTokenEntries, .jpConfig(pipeline)),
                "ensemble",
                grp,
                data,
                args$cisWindow
            )
        )
    }
    records
}

# Entries for one token on one group: reuse the resume cache when it fully
# covers the group's conditions, else fit.
# @noRd
.jointGroupTokenEntries <- function(g, pipeline, token, fitArgs, cond, args) {
    entries <- .jointTokenCacheLookup(pipeline, token, cond, args$cache)
    if (is.null(entries)) {
        entries <- fitJointGroup(g, pipeline, token, fitArgs)
    }
    entries
}

# All-or-nothing resume-cache lookup for a token across a group's conditions
# (twas uses the TwasWeights cache, FM the QtlFineMappingResult cache). NULL
# unless every condition is cached.
# @noRd
.jointTokenCacheLookup <- function(pipeline, token, cond, cache) {
    if (is.null(cache)) {
        return(NULL)
    }
    lookup <- if (is(pipeline, "TwasJointPipeline")) {
        .twasCacheLookup
    } else {
        .fmCacheLookup
    }
    cached <- map(
        seq_len(nrow(cond)),
        .jointCacheLookupAt,
        lookup = lookup,
        cache = cache,
        cond = cond,
        token = token
    )
    if (any(map_lgl(cached, is.null))) NULL else cached
}

# SR-TWAS ensemble LAYER (twas only): combine a group's per-method per-condition
# fits into ensemble per-condition entries -- built ON TOP of the shared per-
# method fitting, never inside it. For each condition r, gather the methods'
# retained out-of-fold CV predictions + weights + R^2, drop methods below the
# R^2 cutoff (stacking needs >= 2), and combine via the `ensembleWeights`
# primitive PER CONTEXT (the sliced single-condition inputs -> contextIndex =
# 1). Returns a length-nCond list of ensemble TwasWeightsRow (NULL where < 2
# methods qualify). All methods share the group's fold partition (the runner
# fixes it before fitting), so their out-of-fold predictions are comparable.
.twasEnsembleLayer <- function(group, perTokenEntries, cfg) {
    tokens <- names(perTokenEntries)
    Y <- .jgY(group)
    r2Cut <- if (is.null(cfg$ensembleR2Threshold)) {
        0.01
    } else {
        cfg$ensembleR2Threshold
    }
    solver <- if (is.null(cfg$ensembleSolver)) {
        "quadprog"
    } else {
        cfg$ensembleSolver
    }
    alpha <- if (is.null(cfg$ensembleAlpha)) 1 else cfg$ensembleAlpha
    stdz <- cfg$standardized
    map(
        seq_len(nrow(.jgConditions(group))),
        .twasEnsembleCondition,
        perTokenEntries = perTokenEntries,
        tokens = tokens,
        Y = Y,
        r2Cut = r2Cut,
        solver = solver,
        alpha = alpha,
        stdz = stdz,
        cfg = cfg
    )
}

# Per-condition CV predictions + weights + R^2 across the group's methods.
# Returns list(preds, wts, rsq), skipping methods with no CV / weights.
# @noRd
.twasEnsembleCollect <- function(r, perTokenEntries, tokens) {
    preds <- list()
    wts <- list()
    rsq <- c()
    for (tk in tokens) {
        e <- perTokenEntries[[tk]][[r]]
        if (is.null(e)) {
            next
        }
        cv <- .rowCvResult(e)
        w <- .rowWeights(e)
        if (is.null(cv) || is.null(cv$predictions) || is.null(w)) {
            next
        }
        pr <- cv$predictions
        preds[[str_c(tk, "_predicted")]] <- matrix(
            as.numeric(pr),
            ncol = 1L,
            dimnames = list(names(pr), NULL)
        )
        wts[[str_c(tk, "_weights")]] <- matrix(
            as.numeric(w),
            ncol = 1L,
            dimnames = list(.rowVariantIds(e), NULL)
        )
        mt <- cv$metrics
        rsq[tk] <- if (!is.null(mt) && is_in("rsq", names(mt))) {
            mt[["rsq"]]
        } else {
            NA_real_
        }
    }
    list(preds = preds, wts = wts, rsq = rsq)
}

# SR-TWAS ensemble entry for one condition: combine the R^2-passing methods
# (>= 2) via ensembleWeights; NULL when too few pass or the solve fails.
# @noRd
.twasEnsembleCondition <- function(
    r,
    perTokenEntries,
    tokens,
    Y,
    r2Cut,
    solver,
    alpha,
    stdz,
    cfg
) {
    coll <- .twasEnsembleCollect(r, perTokenEntries, tokens)
    passing <- names(coll$rsq)[!is.na(coll$rsq) & coll$rsq >= r2Cut]
    if (length(passing) < 2L) {
        return(NULL)
    }
    ens <- tryCatch(
        ensembleWeights(
            cvResults = list(
                prediction = coll$preds[str_c(passing, "_predicted")]
            ),
            Y = Y[, r],
            twasWeightList = coll$wts[str_c(passing, "_weights")],
            contextIndex = 1,
            solver = solver,
            alpha = alpha
        ),
        error = function(err) NULL
    )
    if (is.null(ens) || is.null(ens$ensembleTwasWeights)) {
        return(NULL)
    }
    ew <- ens$ensembleTwasWeights
    vids <- if (!is.null(names(ew))) names(ew) else rownames(ew)
    if (is.null(vids)) {
        vids <- getVariantIds(perTokenEntries[[passing[1L]]][[r]])
    }
    twasWeightsRow(
        variantIds = vids,
        weights = as.numeric(ew),
        cvResult = list(
            methodCoef = ens$methodCoef,
            methodPerformance = ens$methodPerformance
        ),
        standardized = stdz,
        dataType = cfg$dataType
    )
}

# ---- wiring table -----------------------------------------------------------
# Valid cells are rows; invalid cells are absences (a lookup miss is the error).
.jointDispatchTable <- list(
    new(
        "JointDispatchCell",
        pattern = "context",
        dataForm = "individual",
        enumerate = .enumCrossContextIndividual,
        minGroup = 2L
    ),
    new(
        "JointDispatchCell",
        pattern = "context",
        dataForm = "sumstats",
        enumerate = .enumCrossContextSumstats,
        minGroup = 2L
    ),
    new(
        "JointDispatchCell",
        pattern = "trait",
        dataForm = "individual",
        enumerate = .enumCrossTraitIndividual,
        minGroup = 2L
    ),
    new(
        "JointDispatchCell",
        pattern = "trait",
        dataForm = "sumstats",
        enumerate = .enumCrossTraitSumstats,
        minGroup = 2L
    ),
    new(
        "JointDispatchCell",
        pattern = "study",
        dataForm = "sumstats",
        enumerate = .enumCrossStudySumstats,
        minGroup = 2L
    ),
    new(
        "JointDispatchCell",
        pattern = "composed",
        dataForm = "individual",
        enumerate = .enumComposedIndividual,
        minGroup = 2L
    ),
    new(
        "JointDispatchCell",
        pattern = "composed",
        dataForm = "sumstats",
        enumerate = .enumComposedSumstats,
        minGroup = 2L
    ),
    # Univariate: per-(context, trait) 1-condition groups (twas individual
    # only),
    # so univariate methods route through the same engine fitter + ensemble
    # layer.
    new(
        "JointDispatchCell",
        pattern = "univariate",
        dataForm = "individual",
        enumerate = .enumUnivariateIndividual,
        minGroup = 1L
    )
)

.lookupJointCell <- function(pattern, dataForm) {
    for (cell in .jointDispatchTable) {
        if (.jcPattern(cell) == pattern && .jcDataForm(cell) == dataForm) {
            return(cell)
        }
    }
    msg <- glue(
        "No joint dispatch cell for pattern='{pattern}', dataForm='{dataForm}'."
    )
    abort(msg)
}

# Run a parsed jointSpecification through the engine: for each spec resolve its
# scope, map its axes to a (pattern, dataForm) cell, and run every requested
# joint method (token) through `.runJointCell`, rbinding the per-spec results.
# Shared by the fm + twas QtlDataset / QtlSumStats / MultiStudy dispatchers --
# the marker (pipeline) selects the result type and the rbind. `args` is the
# per-run engine payload (twasWeights, methodArgs, cisWindow, region, ...).
.runJointSpecs <- function(
    parsedJointSpec,
    data,
    dataForm,
    pipeline,
    jointMethods,
    contexts,
    traitIds,
    args = list()
) {
    if (length(jointMethods) == 0L || length(parsedJointSpec) == 0L) {
        return(NULL)
    }
    ldSketch <- .jpConfig(pipeline)$ldSketch
    isFm <- is(pipeline, "FmJointPipeline")
    out <- NULL
    for (spec in parsedJointSpec) {
        res <- .runOneJointSpec(
            spec,
            data,
            dataForm,
            pipeline,
            jointMethods,
            contexts,
            traitIds,
            args
        )
        if (is.null(res)) {
            next
        }
        out <- if (is.null(out)) {
            res
        } else if (isFm) {
            .rbindFineMappingResult(out, res, ldSketch = ldSketch)
        } else {
            .rbindTwasWeights(out, res, ldSketch = ldSketch)
        }
    }
    out
}

# Run all joint methods for ONE spec: resolve its scope (optionally
# region-restricting the traits), pick the joint cell, and dispatch to
# .runJointCell (one call with ALL methods so the twas ensemble layer can
# combine a group's per-method fits).
# @noRd
.runOneJointSpec <- function(
    spec,
    data,
    dataForm,
    pipeline,
    jointMethods,
    contexts,
    traitIds,
    args
) {
    scope <- .fmResolveSpecScope(
        spec,
        data,
        contexts = contexts,
        traitIds = traitIds
    )
    if (
        dataForm == "individual" &&
            is.null(traitIds) &&
            !is.null(args$region)
    ) {
        scope <- .jointRestrictRegionTraits(scope, data, args$region)
    }
    pattern <- if (length(spec$axes) > 1L) "composed" else spec$axes[[1L]]
    cell <- .lookupJointCell(pattern, dataForm)
    spArgs <- c(args, list(axes = spec$axes))
    .runJointCell(cell, pipeline, data, scope, jointMethods, spArgs)
}

# Region mode without an explicit traitId: restrict each study's scoped traits
# to the genes overlapping the locus (matches fineMappingPipeline). Gene coords
# are context-independent, so the first scoped context's SE supplies them.
# @noRd
.jointRestrictRegionTraits <- function(scope, data, region) {
    for (st in names(scope$traits)) {
        ctxs <- scope$contexts[[st]]
        if (length(ctxs) == 0L) {
            next
        }
        se <- getPhenotypes(data, contexts = ctxs[[1L]])
        scope$traits[[st]] <- .fmTraitsInRegion(
            se,
            intersect(scope$traits[[st]], rownames(se)),
            region
        )
    }
    scope
}

# Individual-level (QtlDataset) input cannot joint over study: studies have
# disjoint samples (cross-study joints live on the sumstats slot). Preserve the
# historical axis-specific error messages.
.jointRejectStudyOnIndividual <- function(parsedJointSpec) {
    for (spec in parsedJointSpec) {
        if (is_in("study", spec$axes)) {
            if (length(spec$axes) > 1L) {
                msg <- glue(
                    "composed joint axes including 'study' require sumstats ",
                    "input."
                )
                abort(msg)
            }
            msg <- glue(
                "jointSpecification with axis 'study' requires sumstats input ",
                "(QtlDataset is a single individual-level study)."
            )
            abort(msg)
        }
    }
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# Condition `r`'s column of a per-method CV prediction matrix (kept 2-D).
# @noRd
.fmCvSliceCol <- function(m, r) {
    m[, r, drop = FALSE]
}

# Condition `r`'s row of a per-method CV performance matrix (kept 2-D).
# @noRd
.fmCvSliceRow <- function(m, r) {
    m[r, , drop = FALSE]
}

# One per-condition TwasWeightsRow from column `r` of a joint weight matrix.
# @noRd
.jointColEntry <- function(r, vids, weights, fitParts, cfg) {
    twasWeightsRow(
        variantIds = vids,
        weights = weights[, r],
        fits = fitParts,
        standardized = TRUE,
        dataType = cfg$dataType
    )
}

# Build one record into a 1-row collection: wrap `entry` in a list and drop any
# NA joint* axis (so the union re-adds it, padded, only when some row is joint).
# @noRd
.jointBuildRecordPart <- function(rec, constructor) {
    rec$entry <- list(rec$entry)
    for (jc in c("jointStudies", "jointContexts", "jointTraits")) {
        if (
            !is.null(rec[[jc]]) &&
                length(rec[[jc]]) == 1L &&
                is.na(rec[[jc]])
        ) {
            rec[[jc]] <- NULL
        }
    }
    exec(constructor, !!!rec)
}

# TRUE when tuple row `r`'s study keeps context `cx` and trait `tid` in scope.
# @noRd
.jointTupleRowInScope <- function(r, studyCol, cx, tid, scope) {
    s <- studyCol[r]
    is_in(cx, scope$contexts[[s]]) && is_in(tid, scope$traits[[s]])
}

# The "study:context:trait" column label for group member index `i`.
# @noRd
.jointGroupColLabel <- function(i, gi) {
    str_c(gi$studyCol[i], gi$contextCol[i], gi$traitCol[i], sep = ":")
}

# One joint output record for condition `i`, or NULL when that entry is absent.
# traitPos = the bare trait position (NULL when a QtlSumStats caller supplied
# none -> the column is omitted and getTraitPosition() reports NA). The
# fine-mapping window is NOT recorded: the element's own span is the region
# (section 4.4), and a stored window would go stale under subsetRegion().
# @noRd
.jointEntryRecordAt <- function(
    i,
    entries,
    cond,
    method,
    grp,
    data,
    cisWindow
) {
    e <- entries[[i]]
    if (is.null(e)) {
        return(NULL)
    }
    ctx <- as.character(cond$context[[i]])
    tid <- as.character(cond$trait[[i]])
    tpos <- .traitPosFor(data, ctx, tid)
    list(
        study = as.character(cond$study[[i]]),
        context = ctx,
        trait = tid,
        method = method,
        entry = e,
        jointStudies = grp$js,
        jointContexts = grp$jc,
        jointTraits = grp$jt,
        traitPos = tpos
    )
}

# Resume-cache lookup for condition `i` via the pipeline-specific `lookup` fn.
# @noRd
.jointCacheLookupAt <- function(i, lookup, cache, cond, token) {
    lookup(
        cache,
        as.character(cond$study[[i]]),
        as.character(cond$context[[i]]),
        as.character(cond$trait[[i]]),
        token
    )
}

# TRUE when a joint group has at least `minGroup` conditions.
# @noRd
.jointGroupMeetsMin <- function(g, minGroup) {
    nrow(.jgConditions(g)) >= minGroup
}
