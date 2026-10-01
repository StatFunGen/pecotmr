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
    c(
        list(samplePartition = cv$samplePartition),
        compact(list(
            prediction = if (!is.null(cv$prediction)) {
                map(cv$prediction, .fmCvSliceCol, r = r)
            },
            performance = if (!is.null(cv$performance)) {
                map(cv$performance, .fmCvSliceRow, r = r)
            }
        ))
    )
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
#' @importFrom rlang try_fetch
.traitPosFor <- function(data, context, trait) {
    if (methods::is(data, "QtlDataset")) {
        se <- try_fetch(
            getPhenotypes(data, contexts = context),
            error = function(cnd) NULL
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
        # colnames(), not names(): a collection's names() are its GRangesList
        # ELEMENT names (empty here), while the per-row columns live in mcols.
        # Reading names() made this guard always fire, so the whole branch
        # below was dead and every sumstats-derived row fell back to the
        # chrUn sentinel instead of its real trait position.
        if (!is_in("traitPos", colnames(data))) {
            return(NULL)
        }
        idx <- which(
            as.character(data$context) == context &
                as.character(data$trait) == trait
        )
        if (length(idx) == 0L) {
            return(NULL)
        }
        # traitPos is a GRanges column and `idx` came from which(), so this
        # is always a length-1 range.
        tp <- data$traitPos[idx[[1L]]]
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
    anchors <- map(
        seq_along(traits),
        .jointAnchorAt,
        data = data,
        contexts = contexts,
        traits = traits,
        kind = kind,
        cisWindow = cisWindow
    )
    anyFound <- any(map_lgl(anchors, "found"))
    chrs <- map_chr(anchors, "chr")
    starts <- map_int(anchors, "start")
    ends <- map_int(anchors, "end")
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

# One trait's anchor, or the chrUn sentinel when it does not resolve. `found`
# records which it was, so the caller can tell "nothing resolved" from "every
# anchor really is chrUn:1-1".
# @noRd
.jointAnchorAt <- function(i, data, contexts, traits, kind, cisWindow) {
    g <- if (kind == "traitPos") {
        .traitPosFor(data, contexts[[i]], traits[[i]])
    } else {
        .fitRegionFor(data, contexts[[i]], traits[[i]], cisWindow)
    }
    if (is.null(g)) {
        return(list(found = FALSE, chr = "chrUn", start = 1L, end = 1L))
    }
    list(
        found = TRUE,
        chr = as.character(GenomicRanges::seqnames(g))[[1L]],
        start = GenomicRanges::start(g)[[1L]],
        end = GenomicRanges::end(g)[[1L]]
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
    if (length(.jgTraitPos(group)) != nCond) {
        msg <- glue(
            "fitJointGroup: fsusie requires per-trait positions ('pos'); ",
            "it is cross-trait individual-level only."
        )
        abort(msg)
    }
    verbose <- if (is.null(cfg$verbose)) 1 else cfg$verbose
    fitArgs <- .fmMergeUserArgs(
        list(X = Xc, Y = Yc, pos = .jgTraitPos(group)),
        "fsusie",
        args$methodArgs[["fsusie"]]
    )
    raw <- exec(fitFsusie, !!!.splitMethodArgs(fitFsusie, fitArgs))
    # Collapse the functional fit to a variants x features weight matrix now
    # (trimming later drops fitted_wc/csd_X); store on $coef so a trimmed fit
    # can still yield TWAS weights.
    fit <- list_assign(
        raw,
        coef = try_fetch(
            fsusieWeights(fsusieFit = raw, variantIds = colnames(Xc)),
            error = function(cnd) NULL
        )
    ) |>
        .setFinemappingFitClass("fsusie")
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
    cv <- cfg$crossValidationArgs %||% list()
    cvFolds <- cv$folds %||% 0L
    if (cvFolds <= 1L) {
        return(NULL)
    }
    cv <- .fmWeightsCv(
        Xc,
        Yc,
        "fsusie",
        args$methodArgs,
        cvFolds,
        samplePartition = cfg$crossValidationArgs$samplePartition,
        coverage = cfg$credibleSetArgs$coverage,
        pos = .jgTraitPos(group),
        verbose = verbose,
        numThreads = cv$threads %||% 1L,
        seed = cfg$seed
    )
    .fmSliceCv(cv, "fsusie")
}

# One fsusie per-condition (trait) FineMappingRow, with its CV slice attached.
# @noRd
.jointFsusieEntry <- function(
    r,
    fit,
    cvM,
    Xc,
    cfg,
    credibleSetArgs,
    fitRetention
) {
    bare <- .fmPostprocessOne(
        fit = fit,
        method = "fsusie",
        dataX = Xc,
        dataY = NULL,
        conditionIdx = r,
        csInput = "fsusie",
        credibleSetArgs = cfg$credibleSetArgs,
        fitRetention = cfg$fitRetention
    )
    e <- if (is.null(cvM)) {
        bare
    } else {
        .fmAttachCv(bare, .fmSliceCvCondition(cvM, r))
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
    if (!.fmScreenActive(args$pipCutoffToSkip)) {
        return(rep(TRUE, nCond))
    }
    keep <- as.logical(.fmSerScreenColumns(Xc, Yc, args$pipCutoffToSkip))
    if (sum(keep) < 2L) {
        if (verbose >= 1) {
            inform(glue(
                "Skipping mvsusie joint fit: < 2 of {nCond} conditions pass ",
                "the SER pre-screen."
            ))
        }
        return(NULL)
    }
    if (sum(keep) < nCond && verbose >= 1) {
        inform(glue(
            "mvsusie joint fit: SER pre-screen kept {sum(keep)} of ",
            "{nCond} conditions."
        ))
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
    mvBaseArgs <- c(
        list(
            X = Xc,
            Y = Ys,
            prior_variance = mvPrior$priorVariance,
            coverage = cfg$credibleSetArgs$coverage
        ),
        compact(list(residual_variance = mvPrior$residualVariance))
    )
    fitArgs <- .fmMergeUserArgs(
        mvBaseArgs,
        "mvsusie",
        args$methodArgs[["mvsusie"]]
    )
    fit <- exec(fitMvsusie, !!!.splitMethodArgs(fitMvsusie, fitArgs)) |>
        .setFinemappingFitClass("mvsusie")
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
    cv <- cfg$crossValidationArgs %||% list()
    cvFolds <- cv$folds %||% 0L
    if (cvFolds <= 1L) {
        return(NULL)
    }
    sp <- cfg$crossValidationArgs$samplePartition %||% mvCv$samplePartition
    mvPriorCv <- .fmBuildMvsusiePriorCv(mvCv, mvFitParts, colnames(Ys), ddCut)
    cv <- .fmWeightsCv(
        Xc,
        Ys,
        "mvsusie",
        args$methodArgs,
        cvFolds,
        samplePartition = sp,
        coverage = cfg$credibleSetArgs$coverage,
        verbose = verbose,
        mvPrior = mvPrior,
        mvPriorCv = mvPriorCv,
        numThreads = cv$threads %||% 1L,
        seed = cfg$seed
    )
    .fmSliceCv(cv, "mvsusie")
}

# One mvsusie per-condition entry (NULL for a screened-out column), sliced at
# the condition's position in the fitted survivor set + its CV slice.
# @noRd
.jointMvEntry <- function(
    i,
    fitted,
    Xc,
    keep,
    survivors,
    cfg,
    credibleSetArgs,
    fitRetention
) {
    if (!keep[i]) {
        return(NULL)
    }
    r <- match(i, survivors)
    bare <- .fmPostprocessOne(
        fit = fitted$fit,
        method = "mvsusie",
        dataX = Xc,
        dataY = NULL,
        conditionIdx = r,
        csInput = "X",
        credibleSetArgs = cfg$credibleSetArgs,
        fitRetention = cfg$fitRetention
    )
    e <- if (is.null(fitted$cvM)) {
        bare
    } else {
        .fmAttachCv(bare, .fmSliceCvCondition(fitted$cvM, r))
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
        # Derived ONCE here and threaded down: .jointRssEntry runs per
        # condition, so deriving inside it would rebuild the same matrix for
        # every column of Z.
        ldMat <- .jgLdMatrix(group)
        fit <- .jointFitMvsusieRss(group, ldMat, cfg, args)
        # One per-condition entry (RSS has no sample folds).
        map(
            seq_len(ncol(.jgZ(group))),
            .jointRssEntry,
            fit = fit,
            ldMat = ldMat,
            cfg = cfg
        )
    }
)

# mvSuSiE-RSS joint fit over summary statistics, with the data-driven reweighted
# prior. Returns the class-tagged fit.
# @noRd
.jointFitMvsusieRss <- function(group, ldMat, cfg, args) {
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
        R = ldMat,
        N = as.numeric(stats::median(.jgN(group))),
        prior_variance = mvPrior$priorVariance,
        coverage = cfg$credibleSetArgs$coverage
    ) |>
        c(compact(list(residual_variance = mvPrior$residualVariance)))
    fitArgs <- .fmMergeUserArgs(
        mvBaseArgs,
        "mvsusie",
        args$methodArgs[["mvsusie"]]
    )
    fit <- exec(
        fitMvsusieRss,
        !!!.splitMethodArgs(fitMvsusieRss, fitArgs)
    )
    .setFinemappingFitClass(fit, "mvsusie")
}

# One RSS per-condition FineMappingRow (csInput = "Xcorr").
# @noRd
# `credibleSet` / `fitRetention` come off `cfg`, which is what the caller
# supplies; they were also declared as formals with no defaults, so nothing
# could ever have passed them without erroring. `group` was unused too.
.jointRssEntry <- function(r, fit, ldMat, cfg) {
    .fmPostprocessOne(
        fit = fit,
        method = "mvsusie",
        dataX = ldMat,
        dataY = NULL,
        conditionIdx = r,
        csInput = "Xcorr",
        credibleSetArgs = cfg$credibleSetArgs,
        fitRetention = cfg$fitRetention
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
        # The fine-mapping CV refits the method on each fold's training rows;
        # those fits are what let a SuSiE-family method be cross-validated at
        # all now that its weight wrappers never fit.
        foldFits = fineMappingCv$foldFits
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
        stdz <- cfg$standardized %||% FALSE
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
    supplied <- if (
        !is.null(args$methodList) && is_in(methodKey, names(args$methodList))
    ) {
        args$methodList[[methodKey]]
    } else if (!is.null(args$methodArgs)) {
        args$methodArgs[[methodKey]]
    } else {
        NULL
    }
    # FM-fit injection: an FM-derived token extracts its weights from the
    # precomputed fine-mapping fit rather than refitting.
    withFit <- .jointTwasInjectFit(
        supplied %||% list(),
        token,
        fittedModels
    )
    if (!isTRUE(cfg$estimatePi) || !is_in(token, c("bayesC", "bayesB"))) {
        return(withFit)
    }
    .jointTwasSpikeSlabPi(withFit, token, Xc, Yc, cond, cfg, stdz)
}

# The method args with the precomputed fine-mapping fit injected, when the
# token has an adapter, a fit exists, and the caller did not pass one.
# @noRd
.jointTwasInjectFit <- function(ma, token, fittedModels) {
    adapter <- .twasFineMappingMethodAdapters[[token]]
    if (
        is.null(adapter) ||
            is.null(fittedModels[[token]]) ||
            !is.null(ma[[adapter$fitArg]])
    ) {
        return(ma)
    }
    list_assign(ma, !!!set_names(list(fittedModels[[token]]), adapter$fitArg))
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
        fitRetention = "slim",
        standardized = stdz,
        dataType = cfg$dataType,
        verbose = 0,
        seed = cfg$seed
    )
    piHat <- as.numeric(estimateSparsity(mrA))
    list_assign(
        ma,
        !!!compact(list(
            pi = if (token == "bayesC" && is.null(ma$pi)) piHat,
            probIn = if (token == "bayesB" && is.null(ma$probIn)) piHat
        ))
    )
}

# The retention level a joint fit runs at. The engine always keeps its fit --
# the joint layers downstream read it -- so "none" from the pipeline means
# "keep the least we can", not "keep nothing".
# @noRd
.jointRetentionLevel <- function(cfg) {
    level <- cfg$fitRetention %||% "slim"
    if (identical(level, "none")) "slim" else level
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
    retention <- .jointRetentionLevel(cfg)
    verbose <- if (is.null(cfg$verbose)) 1 else cfg$verbose
    tw <- learnTwasWeights(
        Xc,
        Yc,
        weightMethods = wm,
        study = as.character(cond$study[1L]),
        context = as.character(cond$context[1L]),
        trait = as.character(cond$trait[1L]),
        fittedModels = fittedModels,
        fitRetention = retention,
        standardized = stdz,
        dataType = cfg$dataType,
        verbose = verbose,
        seed = cfg$seed
    )
    base <- .twrRowParts(tw, 1L)
    vids <- .twrPartsVariantIds(base)
    raw <- getWeights(base)
    W <- if (is.matrix(raw)) {
        raw
    } else {
        matrix(raw, ncol = nCond, dimnames = list(vids, NULL))
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
# Whether `token` should be cross-validated. `cvWeightMethods` is the caller's
# explicit override of the CV method set; NULL (the default) means "every
# method that produced non-zero weights", which is what the per-method
# all-zero check below enforces.
# @noRd
.jointTwasCvRequested <- function(cvWeightMethods, token) {
    if (is.null(cvWeightMethods)) {
        return(TRUE)
    }
    requested <- if (is.list(cvWeightMethods)) {
        names(cvWeightMethods)
    } else {
        as.character(cvWeightMethods)
    }
    # Accept the short token, the `<token>_weights` method key, or the
    # camelCase weight function. Suffix-stripping alone is not enough for a
    # multi-word token: `susie_inf_weights` strips to `susie_inf`, which is
    # not the canonical `susieInf`.
    canonical <- map_chr(requested, .twasFmTokenFor)
    bare <- str_remove(requested, "(_weights|Weights)$")
    is_in(token, c(canonical[!is.na(canonical)], bare))
}

# @noRd
# Whether this token is barred from cross-validation: not requested, or its
# weights are all zero so there is nothing to validate. A fine-mapping method
# with no recoverable per-fold fits is an error rather than a silent skip.
# @noRd
.jointTwasCvBlocked <- function(W, cfg, token) {
    if (!.jointTwasCvRequested(cfg$crossValidationArgs$weightMethods, token)) {
        return(TRUE)
    }
    if (!is.null(W) && all(W == 0)) {
        # A method whose weights are all zero contributes nothing to
        # cross-validation, and dropping it silently made an empty ensemble
        # look like a modelling result.
        warn(glue(
            "twasWeightsPipeline: method '{token}' is excluded from ",
            "cross-validation because all of its weights are zero."
        ))
        return(TRUE)
    }
    if (is_in(token, names(.twasFineMappingMethodAdapters))) {
        # No handoff means this tuple's fine-mapping entry carries no CV, and
        # nothing here fine-maps -- so the fold fits cannot be recovered.
        # Refusing beats refitting behind the user's back.
        abort(glue(
            "twasWeightsPipeline: cross-validating method '{token}' needs ",
            "each fold's own fine-mapping fit, and the supplied ",
            "fineMappingResult has no cross-validation for this ",
            "(study, context, trait). Run fineMappingPipeline() with ",
            "cvFolds > 1."
        ))
    }
    FALSE
}

# twasWeightsCv() over this cell, with the CV settings resolved from the
# per-call args first and the pipeline config as the fallback.
# @noRd
.jointTwasRunCv <- function(Xc, Yc, wm, args, cfg, cvFolds) {
    cv <- cfg$crossValidationArgs %||% list()
    maxCv <- if ((cv$maxVariants %||% -1) <= 0) {
        Inf
    } else {
        cv$maxVariants
    }
    twasWeightsCv(
        Xc,
        Yc,
        fold = cvFolds,
        samplePartitions = args$samplePartition %||% cv$samplePartition,
        weightMethods = wm,
        fitRetention = "slim",
        maxNumVariants = maxCv,
        numThreads = cfg$crossValidationArgs$threads %||% 1,
        dataDrivenPriorMatricesCv = args$dataDrivenPriorMatricesCv,
        verbose = cfg$verbose %||% 1,
        seed = cfg$seed
    )
}

.jointTwasCv <- function(Xc, Yc, wm, ma, W, args, cfg, token) {
    cv <- cfg$crossValidationArgs %||% list()
    cvFolds <- cv$folds %||% 0L
    if (cvFolds <= 1L) {
        return(NULL)
    }
    # A fine-mapping handoff already carries this token's per-fold fits.
    cvRes <- .twasFmHandoffCv(args$fineMappingCv, token)
    if (!is.null(cvRes)) {
        return(cvRes)
    }
    if (.jointTwasCvBlocked(W, cfg, token)) {
        return(NULL)
    }
    .jointTwasLeakageWarn(args, ma)
    .jointTwasCvResult(
        .jointTwasRunCv(Xc, Yc, wm, args, cfg, cvFolds),
        token
    )
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
        weights <- mrmashRssWeights(
            stat = list(z = .jgZ(group), n = .jgN(group)),
            LD = .jgLdMatrix(group),
            fitRetention = .jointRetentionLevel(cfg)
        )
        vids <- rownames(weights) %||% rownames(.jgZ(group))
        fitParts <- attr(weights, "fit")
        # A single-condition fit comes back as a bare vector; the per-condition
        # split below reads it by column either way.
        wMatrix <- if (is.matrix(weights)) {
            weights
        } else {
            matrix(
                weights,
                ncol = ncol(.jgZ(group)),
                dimnames = list(vids, NULL)
            )
        }
        # One per-condition entry: that condition's weight column + the shared
        # fit.
        map(
            seq_len(ncol(wMatrix)),
            .jointColEntry,
            vids = vids,
            weights = wMatrix,
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

setMethod("construct", "FmJointPipeline", function(pipeline, records) {
    .buildJointResult(
        QtlFineMappingResult,
        records,
        .jpConfig(pipeline)$ldSketch
    )
})

setMethod("construct", "TwasJointPipeline", function(pipeline, records) {
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
    compact(map(
        scopedTraits,
        .enumCrossContextGroupFor,
        data = data,
        study = study,
        scopedContexts = scopedContexts,
        args = args,
        verbose = verbose
    ))
}

# One trait's cross-context group, or NULL when it has no usable (X, Y).
# @noRd
.enumCrossContextGroupFor <- function(
    tid,
    data,
    study,
    scopedContexts,
    args,
    verbose
) {
    xy <- .buildIndividualCrossContextXy(
        data,
        tid,
        scopedContexts,
        args$cisWindow,
        verbose,
        label = "jointCrossContext",
        region = args$region,
        residualizationArgs = args$residualizationArgs
    )
    if (is.null(xy)) {
        return(NULL)
    }
    new(
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

# cross-context / sumstats.
.enumCrossContextSumstats <- function(data, scope, args = list()) {
    ldSketch <- getLdSketch(data)
    studyCol <- as.character(data$study)
    contextCol <- as.character(data$context)
    traitCol <- as.character(data$trait)
    .jeConcat(map(
        scope$studies,
        .enumCrossContextSumstatsForStudy,
        data = data,
        scope = scope,
        args = args,
        ldSketch = ldSketch,
        studyCol = studyCol,
        contextCol = contextCol,
        traitCol = traitCol
    ))
}

# Concatenate per-item lists, empty-safe.
# @noRd
.jeConcat <- function(pieces) {
    if (length(pieces) == 0L) {
        return(list())
    }
    list_c(pieces)
}

# One study's cross-context groups (none when it has fewer than two contexts).
# @noRd
.enumCrossContextSumstatsForStudy <- function(
    s,
    data,
    scope,
    args,
    ldSketch,
    studyCol,
    contextCol,
    traitCol
) {
    scopedContexts <- scope$contexts[[s]]
    if (length(scopedContexts) < 2L) {
        return(list())
    }
    compact(map(
        scope$traits[[s]],
        .enumCrossContextSumstatsGroup,
        data = data,
        s = s,
        scopedContexts = scopedContexts,
        args = args,
        ldSketch = ldSketch,
        studyCol = studyCol,
        contextCol = contextCol,
        traitCol = traitCol
    ))
}

# One (study, trait) cross-context group, or NULL when fewer than two of its
# contexts carry the trait.
# @noRd
.enumCrossContextSumstatsGroup <- function(
    tid,
    data,
    s,
    scopedContexts,
    args,
    ldSketch,
    studyCol,
    contextCol,
    traitCol
) {
    tupleRows <- which(
        studyCol == s & traitCol == tid & is_in(contextCol, scopedContexts)
    )
    if (length(tupleRows) < 2L) {
        return(NULL)
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
    new(
        "SumStatsJointGroup",
        conditions = tibble(study = s, context = ctxNames, trait = tid),
        Z = jz$Z,
        ldSketch = ldSketch,
        N = jz$nVec
    )
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
    compact(map(
        scopedContexts,
        .enumCrossTraitGroupFor,
        data = data,
        study = study,
        scopedTraits = scopedTraits,
        args = args,
        verbose = verbose
    ))
}

# One context's cross-trait group, or NULL when it has no usable (X, Y).
# @noRd
.enumCrossTraitGroupFor <- function(
    cx,
    data,
    study,
    scopedTraits,
    args,
    verbose
) {
    xy <- .buildIndividualCrossTraitXy(
        data,
        cx,
        scopedTraits,
        args$cisWindow,
        verbose,
        label = "jointCrossTrait",
        study = study,
        region = args$region,
        residualizationArgs = args$residualizationArgs
    )
    if (is.null(xy)) {
        return(NULL)
    }
    # Functional positions (one per trait column) for fsusie's domain; mvsusie
    # ignores them. Reordered to match the trait order of Y.
    rr <- SummarizedExperiment::rowRanges(xy$se)[
        match(colnames(xy$Y), rownames(xy$se))
    ]
    new(
        "IndividualJointGroup",
        conditions = tibble(
            study = study,
            context = cx,
            trait = xy$traitsHere
        ),
        X = xy$X,
        Y = xy$Y,
        traitPos = as.numeric(
            (GenomicRanges::start(rr) + GenomicRanges::end(rr)) / 2
        )
    )
}

# cross-trait / sumstats.
.enumCrossTraitSumstats <- function(data, scope, args = list()) {
    ldSketch <- getLdSketch(data)
    studyCol <- as.character(data$study)
    contextCol <- as.character(data$context)
    traitCol <- as.character(data$trait)
    .jeConcat(map(
        scope$studies,
        .enumCrossTraitSumstatsForStudy,
        data = data,
        scope = scope,
        args = args,
        ldSketch = ldSketch,
        studyCol = studyCol,
        contextCol = contextCol,
        traitCol = traitCol
    ))
}

# @noRd
.enumCrossTraitSumstatsForStudy <- function(
    s,
    data,
    scope,
    args,
    ldSketch,
    studyCol,
    contextCol,
    traitCol
) {
    compact(map(
        scope$contexts[[s]],
        .enumCrossTraitSumstatsGroup,
        data = data,
        s = s,
        scopedTraits = scope$traits[[s]],
        args = args,
        ldSketch = ldSketch,
        studyCol = studyCol,
        contextCol = contextCol,
        traitCol = traitCol
    ))
}

# One (study, context) cross-trait group, or NULL when fewer than two of its
# traits are present.
# @noRd
.enumCrossTraitSumstatsGroup <- function(
    cx,
    data,
    s,
    scopedTraits,
    args,
    ldSketch,
    studyCol,
    contextCol,
    traitCol
) {
    tupleRows <- which(
        studyCol == s & contextCol == cx & is_in(traitCol, scopedTraits)
    )
    if (length(tupleRows) < 2L) {
        return(NULL)
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
    new(
        "SumStatsJointGroup",
        conditions = tibble(study = s, context = cx, trait = trNames),
        Z = jz$Z,
        ldSketch = ldSketch,
        N = jz$nVec
    )
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
    allCtxs <- unique(unname(list_c(scope$contexts)))
    allTrs <- unique(unname(list_c(scope$traits)))
    .jeConcat(map(
        allCtxs,
        .enumCrossStudyForContext,
        data = data,
        scope = scope,
        args = args,
        cols = cols,
        allTrs = allTrs,
        ldSketch = ldSketch
    ))
}

# One context's cross-study groups, one per trait that has enough studies.
# @noRd
.enumCrossStudyForContext <- function(
    cx,
    data,
    scope,
    args,
    cols,
    allTrs,
    ldSketch
) {
    compact(map(
        allTrs,
        .enumCrossStudyGroupFor,
        data = data,
        scope = scope,
        args = args,
        cols = cols,
        cx = cx,
        ldSketch = ldSketch
    ))
}

# `map()` hands the trait first; .enumCrossStudyGroup takes it sixth.
# @noRd
.enumCrossStudyGroupFor <- function(
    tid,
    data,
    scope,
    args,
    cols,
    cx,
    ldSketch
) {
    .enumCrossStudyGroup(data, scope, args, cols, cx, tid, ldSketch)
}

# One cross-study joint group for a (context, trait) cell, or NULL when fewer
# than two studies survive the scope filter -- a "joint" fit over one study is
# just that study, so the cell contributes nothing.
# @noRd
.enumCrossStudyGroup <- function(data, scope, args, cols, cx, tid, ldSketch) {
    candidates <- which(
        cols$context == cx &
            cols$trait == tid &
            is_in(cols$study, scope$studies)
    )
    keep <- map_lgl(
        candidates,
        .jointTupleRowInScope,
        studyCol = cols$study,
        cx = cx,
        tid = tid,
        scope = scope
    )
    tupleRows <- candidates[keep]
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
        ldSketch = ldSketch,
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
        region = args$region,
        residualizationArgs = args$residualizationArgs
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
    .jeConcat(map(
        scope$contexts[[study]],
        .enumUnivariateForContext,
        data = data,
        scope = scope,
        args = args,
        study = study,
        naAction = naAction
    ))
}

# One context's univariate groups, one per trait it carries.
# @noRd
.enumUnivariateForContext <- function(cx, data, scope, args, study, naAction) {
    se <- getPhenotypes(data, contexts = cx)
    compact(map(
        intersect(scope$traits[[study]], rownames(se)),
        .enumUnivariateGroupFor,
        data = data,
        args = args,
        study = study,
        cx = cx,
        naAction = naAction
    ))
}

# One (context, trait) univariate group, or NULL when fewer than two samples
# are shared between its genotypes and its phenotype.
# @noRd
.enumUnivariateGroupFor <- function(tid, data, args, study, cx, naAction) {
    Y <- .fmResidPheno(
        data,
        contexts = cx,
        traitId = tid,
        naAction = naAction,
        residualizationArgs = args$residualizationArgs
    )
    X <- if (is.null(args$region)) {
        .fmResidGeno(
            data,
            contexts = cx,
            traitId = tid,
            cisWindow = args$cisWindow,
            residualizationArgs = args$residualizationArgs
        )
    } else {
        .fmResidGeno(
            data,
            contexts = cx,
            region = args$region,
            residualizationArgs = args$residualizationArgs
        )
    }
    common <- intersect(rownames(X), rownames(Y))
    if (length(common) < 2L) {
        return(NULL)
    }
    new(
        "IndividualJointGroup",
        conditions = tibble(study = study, context = cx, trait = tid),
        X = X[common, , drop = FALSE],
        Y = Y[common, , drop = FALSE]
    )
}

# composed / sumstats: general N-axis joint. `args$axes` (subset of study /
# context / trait) names the collapsed axes; rows split by the complement
# (fixed) axes form one group each. Reuses .enumerateComposedSumstatGroups.
.enumComposedSumstats <- function(data, scope, args = list()) {
    axes <- args$axes %||% c("context", "trait")
    ldSketch <- getLdSketch(data)
    gi <- .enumerateComposedSumstatGroups(list(axes = axes), data, scope)
    if (is.null(gi)) {
        return(list())
    }
    # `gi$groups` comes from split(), so it is NAMED; map() would carry those
    # names onto the groups, which the record assembly downstream does not
    # expect (the loop this replaced appended positionally).
    unname(compact(map(
        gi$groups,
        .enumComposedGroupFor,
        data = data,
        gi = gi,
        args = args,
        ldSketch = ldSketch
    )))
}

# One composed cell's group, or NULL when it spans fewer than two tuples --
# a "joint" fit over one tuple is just that tuple.
# @noRd
.enumComposedGroupFor <- function(gIdx, data, gi, args, ldSketch) {
    if (length(gIdx) < 2L) {
        return(NULL)
    }
    jz <- .buildJointSumstatZMatrix(
        data,
        gIdx,
        map_chr(gIdx, .jointGroupColLabel, gi = gi),
        errorLabel = "composed (QtlSumStats)",
        ldSketch = ldSketch,
        cutoffs = args$cutoffs
    )
    new(
        "SumStatsJointGroup",
        conditions = tibble(
            study = gi$studyCol[gIdx],
            context = gi$contextCol[gIdx],
            trait = gi$traitCol[gIdx]
        ),
        Z = jz$Z,
        ldSketch = ldSketch,
        N = jz$nVec
    )
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
    out <- list_assign(args, !!!.twasGroupFmArgs(args, cond))
    cv <- cfg$crossValidationArgs %||% list()
    cvF <- cv$folds %||% 0L
    if (cvF <= 1L || !is(g, "IndividualJointGroup")) {
        return(out)
    }
    list_assign(
        out,
        samplePartition = .jointCvPartition(
            fmCv = out$fineMappingCv,
            userSp = args$samplePartition %||% cv$samplePartition,
            sampleIds = rownames(.jgX(g)),
            cvFolds = cvF
        )
    )
}

# The group's fine-mapping fits + CV, keyed on its first condition. Empty when
# the caller supplied no fine-mapping result.
# @noRd
.twasGroupFmArgs <- function(args, cond) {
    fmRes <- args$fineMappingResult
    if (is.null(fmRes)) {
        return(list())
    }
    s1 <- as.character(cond$study[[1L]])
    c1 <- as.character(cond$context[[1L]])
    t1 <- as.character(cond$trait[[1L]])
    nR <- if (is.null(args$nRegions)) 1L else args$nRegions
    bi <- if (is.null(args$regionIndex)) 1L else args$regionIndex
    af <- .twasFineMappingFits(fmRes, study = s1, context = c1, trait = t1)
    list(
        fittedModels = if (is.null(af)) {
            list()
        } else {
            .twasFitsForRegion(af, bi, nR)
        },
        fineMappingCv = .twasCvResultFor(fmRes, s1, c1, t1)
    )
}

# The fold partition every method in this group is scored on. A fine-mapping
# CV carries fits tied to its own folds, so when one is present it governs:
# the other weight methods are cross-validated on those same folds, which is
# what makes their scores comparable to the fine-mapping method's.
# @noRd
.jointCvPartition <- function(fmCv, userSp, sampleIds, cvFolds) {
    fmSp <- fmCv$samplePartition
    if (is.null(fmSp)) {
        if (!is.null(userSp)) {
            return(userSp)
        }
        return(.normalizeCvFolds(cvFolds, NULL, sampleIds)$samplePartition)
    }
    .jointCheckFmCvSamples(fmSp, sampleIds)
    .jointCheckPartitionAgreement(fmSp, userSp)
    fmSp
}

# Every sample the fine-mapping folds name must be a sample of the dataset
# being scored. Otherwise a fold's held-out rows cannot be located here, and
# the fits would be scored against samples they were never separated from.
# @noRd
.jointCheckFmCvSamples <- function(fmSp, sampleIds) {
    unmatched <- setdiff(as.character(fmSp$Sample), as.character(sampleIds))
    if (length(unmatched) == 0L) {
        return(invisible(NULL))
    }
    shown <- str_flatten(head(unmatched, 5L), ", ")
    more <- if (length(unmatched) > 5L) {
        str_c(" (and ", length(unmatched) - 5L, " more)")
    } else {
        ""
    }
    msg <- glue(
        "twasWeightsPipeline: {length(unmatched)} sample(s) in the ",
        "fine-mapping cross-validation folds are not in this dataset: ",
        "{shown}{more}. The fold fits must come from a fineMappingPipeline() ",
        "run on these same samples."
    )
    abort(msg)
}

# An explicitly supplied partition that disagrees with the fine-mapping CV's
# is ambiguous -- the fold fits belong to one of them, and silently choosing
# would score some methods on folds the others never saw.
# @noRd
.jointCheckPartitionAgreement <- function(fmSp, userSp) {
    if (is.null(userSp)) {
        return(invisible(NULL))
    }
    if (identical(.cvPartitionKey(userSp), .cvPartitionKey(fmSp))) {
        return(invisible(NULL))
    }
    msg <- glue(
        "twasWeightsPipeline: the supplied `samplePartition` differs from ",
        "the fine-mapping cross-validation's own folds. The fine-mapping ",
        "fold fits belong to that partition, so pass it (or omit ",
        "`samplePartition` and let it be reused)."
    )
    abort(msg)
}

# Append one output record per (condition, method) to the joint-rows accumulator
# `rows` (mutated by reference), resolving each condition's fine-mapping region
# + trait position. `grp` bundles the per-group invariants list(cond, js, jc,
# jt).
# @noRd
.jointEntryRecords <- function(entries, method, grp, data) {
    cond <- grp$cond
    recs <- map(
        seq_len(min(length(entries), nrow(cond))),
        .jointEntryRecordAt,
        entries = entries,
        cond = cond,
        method = method,
        grp = grp,
        data = data
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
    groups <- .jcEnumerate(cell)(data, scope, args) |>
        keep(.jointGroupMeetsMin, minGroup = .jcMinGroup(cell))
    if (length(groups) == 0L) {
        return(NULL)
    }
    doEnsemble <- is(pipeline, "TwasJointPipeline") &&
        isTRUE(.jpConfig(pipeline)$ensembleArgs$enabled)
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
    perTokenEntries <- compact(set_names(
        map(
            tokens,
            .jointTokenEntriesOrNull,
            g = g,
            pipeline = pipeline,
            fitArgs = fitArgs,
            cond = cond,
            args = args
        ),
        tokens
    ))
    tokenRecords <- .jeConcat(map(
        names(perTokenEntries),
        .jointTokenRecords,
        perTokenEntries = perTokenEntries,
        grp = grp,
        data = data
    ))
    if (!doEnsemble || length(perTokenEntries) < 2L) {
        return(tokenRecords)
    }
    c(
        tokenRecords,
        .jointEntryRecords(
            .twasEnsembleLayer(g, perTokenEntries, .jpConfig(pipeline)),
            "ensemble",
            grp,
            data
        )
    )
}

# One token's entries, or NULL when it produced none.
# @noRd
.jointTokenEntriesOrNull <- function(token, g, pipeline, fitArgs, cond, args) {
    entries <- .jointGroupTokenEntries(
        g,
        pipeline,
        token,
        fitArgs,
        cond,
        args
    )
    if (is.null(entries) || length(entries) == 0L) {
        return(NULL)
    }
    entries
}

# @noRd
.jointTokenRecords <- function(token, perTokenEntries, grp, data) {
    .jointEntryRecords(perTokenEntries[[token]], token, grp, data)
}

# Entries for one token on one group: reuse the resume cache when it fully
# covers the group's conditions, else fit.
# @noRd
.jointGroupTokenEntries <- function(g, pipeline, token, fitArgs, cond, args) {
    entries <- .jointTokenCacheLookup(pipeline, token, cond, args$cache) %||%
        fitJointGroup(g, pipeline, token, fitArgs)
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
    r2Cut <- if (is.null(cfg$ensembleArgs$r2Threshold)) {
        0.01
    } else {
        cfg$ensembleArgs$r2Threshold
    }
    solver <- if (is.null(cfg$ensembleArgs$solver)) {
        "quadprog"
    } else {
        cfg$ensembleArgs$solver
    }
    alpha <- if (is.null(cfg$ensembleArgs$alpha)) 1 else cfg$ensembleArgs$alpha
    stdz <- cfg$standardized %||% FALSE
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
# One token's ensemble inputs for condition `r`, or NULL when it has no
# cross-validated predictions or no weights to contribute.
# @noRd
.twasEnsembleTokenPart <- function(tk, r, perTokenEntries) {
    e <- perTokenEntries[[tk]][[r]]
    if (is.null(e)) {
        return(NULL)
    }
    cv <- .rowCvResult(e)
    w <- .rowWeights(e)
    if (is.null(cv) || is.null(cv$predictions) || is.null(w)) {
        return(NULL)
    }
    pr <- cv$predictions
    list(
        token = tk,
        pred = matrix(
            as.numeric(pr),
            ncol = 1L,
            dimnames = list(names(pr), NULL)
        ),
        weights = matrix(
            as.numeric(w),
            ncol = 1L,
            dimnames = list(.rowVariantIds(e), NULL)
        ),
        rsq = if (!is.null(cv$metrics) && is_in("rsq", names(cv$metrics))) {
            cv$metrics[["rsq"]]
        } else {
            NA_real_
        }
    )
}

.twasEnsembleCollect <- function(r, perTokenEntries, tokens) {
    parts <- compact(map(
        tokens,
        .twasEnsembleTokenPart,
        r = r,
        perTokenEntries = perTokenEntries
    ))
    if (length(parts) == 0L) {
        return(list(preds = list(), wts = list(), rsq = c()))
    }
    contributing <- map_chr(parts, "token")
    list(
        preds = set_names(
            map(parts, "pred"),
            str_c(contributing, "_predicted")
        ),
        wts = set_names(
            map(parts, "weights"),
            str_c(contributing, "_weights")
        ),
        rsq = set_names(map_dbl(parts, "rsq"), contributing)
    )
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
    ens <- try_fetch(
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
        error = function(cnd) NULL
    )
    if (is.null(ens) || is.null(ens$ensembleTwasWeights)) {
        return(NULL)
    }
    ew <- ens$ensembleTwasWeights
    vids <- (if (!is.null(names(ew))) names(ew) else rownames(ew)) %||%
        getVariantIds(perTokenEntries[[passing[1L]]][[r]])
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

# Top-PC enumeration: PCA-reduce each context's multi-trait phenotype and make
# ONE group per top principal component, each a single-column Y named topPCk.
# Structurally identical to .enumUnivariateIndividual -- same residualization
# helpers, same IndividualJointGroup -- so everything downstream (CV,
# ensemble, retention) applies unchanged. It is the TWAS peer of
# fineMappingPipeline's .fmPcaContextRows, sharing .fmTopPcScores().
# @noRd
.enumTopPcIndividual <- function(data, scope, args = list()) {
    study <- getStudy(data)
    if (!is_in(study, scope$studies)) {
        return(list())
    }
    naAction <- if (is.null(args$naAction)) "drop" else args$naAction
    .jeConcat(map(
        scope$contexts[[study]],
        .enumTopPcForContext,
        data = data,
        scope = scope,
        args = args,
        study = study,
        naAction = naAction
    ))
}

# @noRd
.enumTopPcForContext <- function(cx, data, scope, args, study, naAction) {
    se <- getPhenotypes(data, contexts = cx)
    traits <- intersect(scope$traits[[study]], rownames(se))
    # PCA is undefined for a single trait, matching the fine-mapping rule.
    if (length(traits) < 2L) {
        return(list())
    }
    Y <- .fmResidPheno(
        data,
        contexts = cx,
        traitId = traits,
        naAction = naAction,
        residualizationArgs = args$residualizationArgs
    )
    scores <- .fmTopPcScores(Y, args$nPCs %||% 10L)
    if (is.null(scores)) {
        return(list())
    }
    X <- .enumTopPcX(data, cx, traits, args)
    common <- intersect(rownames(X), rownames(scores))
    if (length(common) < 2L) {
        return(list())
    }
    compact(map(
        colnames(scores),
        .enumTopPcGroup,
        X = X[common, , drop = FALSE],
        scores = scores[common, , drop = FALSE],
        study = study,
        cx = cx
    ))
}

# @noRd
.enumTopPcX <- function(data, cx, traits, args) {
    if (is.null(args$region)) {
        return(.fmResidGeno(
            data,
            contexts = cx,
            traitId = traits,
            cisWindow = args$cisWindow,
            residualizationArgs = args$residualizationArgs
        ))
    }
    .fmResidGeno(
        data,
        contexts = cx,
        region = args$region,
        residualizationArgs = args$residualizationArgs
    )
}

# @noRd
.enumTopPcGroup <- function(pcName, X, scores, study, cx) {
    new(
        "IndividualJointGroup",
        conditions = tibble(study = study, context = cx, trait = pcName),
        X = X,
        Y = scores[, pcName, drop = FALSE]
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
    ),
    new(
        "JointDispatchCell",
        pattern = "topPc",
        dataForm = "individual",
        enumerate = .enumTopPcIndividual,
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
    results <- compact(map(
        parsedJointSpec,
        .runOneJointSpecFor,
        data = data,
        dataForm = dataForm,
        pipeline = pipeline,
        jointMethods = jointMethods,
        contexts = contexts,
        traitIds = traitIds,
        args = args
    ))
    if (length(results) == 0L) {
        return(NULL)
    }
    rbindFn <- if (isFm) .rbindFineMappingResult else .rbindTwasWeights
    reduce(
        results,
        .jointRbindWithSketch,
        rbindFn = rbindFn,
        ldSketch = ldSketch
    )
}

# `map()` hands the spec first, which is also where .runOneJointSpec wants it.
# @noRd
.runOneJointSpecFor <- function(
    spec,
    data,
    dataForm,
    pipeline,
    jointMethods,
    contexts,
    traitIds,
    args
) {
    .runOneJointSpec(
        spec,
        data,
        dataForm,
        pipeline,
        jointMethods,
        contexts,
        traitIds,
        args
    )
}

# @noRd
.jointRbindWithSketch <- function(acc, res, rbindFn, ldSketch) {
    rbindFn(acc, res, ldSketch = ldSketch)
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
    resolved <- .fmResolveSpecScope(
        spec,
        data,
        contexts = contexts,
        traitIds = traitIds
    )
    restrictToRegion <- dataForm == "individual" &&
        is.null(traitIds) &&
        !is.null(args$region)
    scope <- if (restrictToRegion) {
        .jointRestrictRegionTraits(resolved, data, args$region)
    } else {
        resolved
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
# One study's traits narrowed to those inside `region`; left alone when the
# study has no context to read a phenotype from.
# @noRd
.jointStudyTraitsInRegion <- function(st, scope, data, region) {
    ctxs <- scope$contexts[[st]]
    if (length(ctxs) == 0L) {
        return(scope$traits[[st]])
    }
    se <- getPhenotypes(data, contexts = ctxs[[1L]])
    .fmTraitsInRegion(
        se,
        intersect(scope$traits[[st]], rownames(se)),
        region
    )
}

.jointRestrictRegionTraits <- function(scope, data, region) {
    studies <- names(scope$traits)
    list_assign(
        scope,
        traits = set_names(
            map(
                studies,
                .jointStudyTraitsInRegion,
                scope = scope,
                data = data,
                region = region
            ),
            studies
        )
    )
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
# A joint-axis column holding a lone NA means "not a joint axis", so it is
# dropped rather than carried as a column of NA.
# @noRd
.jointAxisIsAbsent <- function(value) {
    !is.null(value) && length(value) == 1L && is.na(value)
}

.jointBuildRecordPart <- function(rec, constructor) {
    jointCols <- c("jointStudies", "jointContexts", "jointTraits")
    absent <- jointCols[map_lgl(rec[jointCols], .jointAxisIsAbsent)]
    kept <- rec[setdiff(names(rec), absent)]
    # `list_assign()`, not `list_modify()`: the entry payload is itself a
    # list, and list_modify would merge into it rather than wrap it.
    exec(constructor, !!!list_assign(kept, entry = list(rec$entry)))
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
.jointEntryRecordAt <- function(i, entries, cond, method, grp, data) {
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
