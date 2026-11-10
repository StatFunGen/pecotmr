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

# The joint key of a group over conditions `cond` (.JOINT_COLUMNS): what its
# rows record and what a lookup of its fit matches on.
# @noRd
.jointKeyOfConditions <- function(cond) {
    set_names(
        map(
            c("study", "context", "trait"),
            .jointAxisMembers,
            conditions = cond
        ),
        .JOINT_COLUMNS
    )
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
#              QtlSumStats -> the tuple's own `region`: the window it was
#              given, or the span of its variants as first built, before QC
#              (a pipeline's cisWindow sets it from traitPos first).
# Each returns a length-1 GRanges, or NULL when the anchor is unavailable (the
# accumulator / builder then records a chrUn sentinel). A row is named by its
# whole (study, context, trait) tuple: a QtlSumStats can carry the same
# (context, trait) in several studies, each with its own position and window.
# A QtlDataset holds one study, so `study` has nothing to choose there; a
# MultiStudyQtlDataset's individual-level study answers for itself.
#' @importFrom rlang try_fetch
.traitPosFor <- function(data, study, context, trait) {
    if (.isMultiStudyComponent(data, study)) {
        return(.traitPosFor(data[[study]], study, context, trait))
    }
    if (methods::is(data, "QtlDataset")) {
        se <- try_fetch(
            .qtlContextSe(data, context = context),
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
        idx <- .qssTupleRows(data, study, context, trait)
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

.fitRegionFor <- function(data, study, context, trait, cisWindow = NULL) {
    if (.isMultiStudyComponent(data, study)) {
        return(.fitRegionFor(data[[study]], study, context, trait, cisWindow))
    }
    if (methods::is(data, "QtlDataset")) {
        tp <- .traitPosFor(data, study, context, trait)
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
        idx <- .qssTupleRows(data, study, context, trait)
        if (length(idx) == 0L) {
            return(NULL)
        }
        return(.rowRegionAt(data, idx[[1L]]))
    }
    NULL
}

# Whether `study` is one of a MultiStudyQtlDataset's individual-level studies.
# @noRd
.isMultiStudyComponent <- function(data, study) {
    methods::is(data, "MultiStudyQtlDataset") && is_in(study, names(data))
}

# The rows of a QtlSumStats holding one (study, context, trait) tuple.
# @noRd
.qssTupleRows <- function(data, study, context, trait) {
    which(
        as.character(data$study) == study &
            as.character(data$context) == context &
            as.character(data$trait) == trait
    )
}

# Vectorised per-row anchors for the pushRow builders: map the scalar helper
# over aligned (study, context, trait) vectors and assemble ONE GRanges of
# length n, using a chrUn sentinel where the anchor is unavailable. Built from
# plain chr/start/end vectors so mixed seqlevels never trip seqinfo merging.
.anchorVector <- function(
    data,
    studies,
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
        studies = studies,
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
    # NULL so the builder omits the column entirely and traitPosition()
    # reports an empty GRanges, rather than a column full of chrUn sentinels.
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
.jointAnchorAt <- function(
    i,
    data,
    studies,
    contexts,
    traits,
    kind,
    cisWindow
) {
    g <- if (kind == "traitPos") {
        .traitPosFor(data, studies[[i]], contexts[[i]], traits[[i]])
    } else {
        .fitRegionFor(
            data,
            studies[[i]],
            contexts[[i]],
            traits[[i]],
            cisWindow
        )
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

# The window a joint group was analysed over. An explicit region wins; failing
# that the group's traits define it, as the union span of their own cis
# windows -- one joint fit over several traits covers all of them.
# @noRd
.jointGroupWindow <- function(data, cond, args) {
    explicit <- .asAnalysisWindow(args$region)
    if (!is.null(explicit)) {
        return(explicit)
    }
    per <- try_fetch(
        .anchorVector(
            data,
            as.character(cond$study),
            as.character(cond$context),
            as.character(cond$trait),
            "region",
            args$cisWindow
        ),
        error = function(cnd) NULL
    )
    if (is.null(per) || length(per) == 0L) {
        return(NULL)
    }
    range(per)[1L]
}


# Per-variant allele frequency for a joint group's design matrix, or NULL.
#
# The univariate/RSS paths pass `af` into postprocessFinemappingFits(); the
# joint ones did not, so every mvsusie / fsusie row exported af = NA (and
# therefore maf = NA) even though the genotypes were right there. The window is
# taken from the fit's OWN variants rather than a trait cis-window: a joint
# group spans several traits, so no single trait's window describes X.
#
# A cross-study group carries its own, one column per condition, since its
# conditions are different studies' samples and it was read from no one block.
# @noRd
.jointAfForGroup <- function(g, block) {
    if (!is(g, "IndividualJointGroup")) {
        return(NULL)
    }
    if (length(.jgAf(g)) > 0L) {
        return(.jgAf(g))
    }
    if (is.null(block$af)) {
        return(NULL)
    }
    block$af[match(colnames(.jgX(g)), colnames(block$X))]
}


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
# cross-trait). Keep the first post-processed entry for this joint fit.
#
# The credible-set arguments supply what .fmOwnedEngineArgs("fsusieR::susiF")
# says is supplied: these are refused in FsusieOptions, so the pipeline is
# their only source. They are fitFsusie() formals, so .splitMethodArgs()
# routes them past the constructor that refuses them.
#
# The functional fit is collapsed to a variants x features weight matrix now
# (trimming later drops fitted_wc) and stored on $coef, so a trimmed fit can
# still yield TWAS weights.
#
# Returns ONE entry for the whole joint fit, unlike mvsusie's
# one-per-condition. mvsusie's posteriors carry a conditions axis, so
# .fmSliceCond() gives each condition its own 2-D slice and N rows partition
# the fit at no cost. An fsusie fit has no such axis: each trait is a POSITION
# on a shared functional grid, so `fitted_func` / `cred_band` / `outing_grid`
# describe all traits at once and cannot be sliced per trait. One row per
# trait would therefore store N copies of the same curve.
#
# The trait list itself is not lost -- every row carries the full joint
# membership in its `jointTraits` column. Only the grid POSITIONS are
# unrecoverable from that, and they are what maps a curve position back to a
# trait, so they ride on the fit.
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
    fit <- if (.jointFitsFullData(cfg)) {
        .jointFsusieFullFit(group, Xc, Yc, cfg, args)
    } else {
        NULL
    }
    cvM <- .jointFsusieCv(Xc, Yc, group, cfg, args, verbose)
    entry <- .jointFsusieEntry(
        r = 1L,
        fit = fit,
        cvM = cvM,
        Xc = Xc,
        cfg = cfg,
        af = args$af,
        region = args$region
    )
    if (!is.null(fit)) {
        entry@susieFit$trait_positions <- .jgTraitPos(group)
    }
    list(entry)
}

# The full-data fSuSiE fit of one joint group, with its functional fit
# collapsed to a variants x features weight matrix on $coef.
# @noRd
.jointFsusieFullFit <- function(group, Xc, Yc, cfg, args) {
    csp <- cfg$credibleSetParam
    fitArgs <- .fmMergeUserArgs(
        list(
            X = Xc,
            Y = Yc,
            pos = .jgTraitPos(group),
            maxNumSingleEffects = csp$maxNumSingleEffects,
            covLev = csp$coverage,
            minPurity = csp$minAbsCorr
        ),
        args$methodArgs[["fsusie"]]
    )
    fitFsusieArgs <- .splitMethodArgs(fitFsusie, fitArgs)
    raw <- exec(fitFsusie, !!!fitFsusieArgs)
    list_assign(
        raw,
        coef = try_fetch(
            fsusieWeights(fsusieFit = raw, variantIds = colnames(Xc)),
            error = function(cnd) NULL
        )
    ) |>
        .setFinemappingFitClass("fsusie")
}

# Per-fold fsusie CV slice, or NULL when CV is disabled.
# @noRd
.jointFsusieCv <- function(Xc, Yc, group, cfg, args, verbose) {
    cv <- cfg$crossValidationParam %||% list()
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
        samplePartition = cfg$crossValidationParam$samplePartition,
        coverage = cfg$credibleSetParam$coverage,
        pos = .jgTraitPos(group),
        verbose = verbose,
        BPPARAM = cv$BPPARAM %||% SerialParam(),
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
    credibleSetParam,
    fitRetention,
    af = NULL,
    region = NULL
) {
    bare <- if (is.null(fit)) {
        .fmCvOnlyRow()
    } else {
        .fmPostprocessOne(
            fit = fit,
            method = "fsusie",
            dataX = Xc,
            dataY = NULL,
            af = af,
            region = region,
            conditionIdx = r,
            csInput = "fsusie",
            credibleSetParam = cfg$credibleSetParam,
            fitRetention = cfg$fitRetention
        )
    }
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
    xy <- .jointMvFitRows(Xc, Yc[, survivors, drop = FALSE])
    fitted <- .jointMvFit(
        group,
        xy$X,
        xy$Y,
        cfg,
        list_assign(args, methodArgs = .jointMvMethodArgs(args)),
        ddCut,
        verbose
    )
    map(
        seq_len(nCond),
        .jointMvEntry,
        fitted = fitted,
        Xc = Xc,
        Yc = Yc,
        keep = keep,
        survivors = survivors,
        cfg = cfg,
        af = args$af,
        region = args$region
    )
}

# The rows a joint fit is given: those observing at least one of the
# conditions it fits. A cross-study group's rows each observe one study, so a
# study the screen drops leaves rows that observe nothing. Every row of any
# other group observes every condition, and the group is passed on as it is.
# @noRd
.jointMvFitRows <- function(Xc, Ys) {
    observed <- rowSums(!is.na(Ys)) > 0L
    if (all(observed)) {
        return(list(X = Xc, Y = Ys))
    }
    list(
        X = Xc[observed, , drop = FALSE],
        Y = Ys[observed, , drop = FALSE]
    )
}

# The run's per-method options, mvsusie's `standardize` defaulting to FALSE
# when the residualizer scaled X and Y within each study
# (`ResidualizationParam(scaleResiduals = TRUE)`, the default), so the shared
# prior components describe per-SD effects. mvsusieR's own standardization
# scales each outcome's prior by sd(y) / sqrt(n): on studies of unequal size,
# its shared components would mean equal effects on the z scale, letting the
# smaller study's per-SD effect be sqrt(n_large / n_small) times larger. Where
# every condition has the same samples the two agree. A user's own
# `standardize` wins.
# @noRd
.jointMvMethodArgs <- function(args) {
    methodArgs <- args$methodArgs %||% list()
    user <- as.list(methodArgs[["mvsusie"]] %||% list())
    scaled <- .qtlScaleResiduals(args$residualizationParam)
    if (!scaled || is_in("standardize", names(user))) {
        return(args$methodArgs)
    }
    list_assign(methodArgs, mvsusie = list_assign(user, standardize = FALSE))
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
        args$mrmashPrior,
        key$study,
        key$trait,
        context = key$context,
        conditionNames = colnames(Ys)
    )
    mvCv <- .fmLookupMrmashCv(
        args$mrmashPrior,
        key$study,
        key$trait,
        context = key$context,
        conditionNames = colnames(Ys)
    )
    mvPrior <- .buildMvsusieReweightedPrior(mvFitParts, colnames(Ys), ddCut)
    fit <- if (.jointFitsFullData(cfg)) {
        .jointMvFullFit(Xc, Ys, cfg, args, mvPrior)
    } else {
        NULL
    }
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

# The full-data mvSuSiE fit of one joint group, under the reweighted prior.
# @noRd
.jointMvFullFit <- function(Xc, Ys, cfg, args, mvPrior) {
    mvBaseArgs <- c(
        list(
            X = Xc,
            Y = Ys,
            prior_variance = mvPrior$priorVariance,
            coverage = cfg$credibleSetParam$coverage
        ),
        compact(list(residual_variance = mvPrior$residualVariance))
    )
    fitArgs <- .fmMergeUserArgs(mvBaseArgs, args$methodArgs[["mvsusie"]])
    fitMvsusieArgs <- .splitMethodArgs(fitMvsusie, fitArgs)
    exec(fitMvsusie, !!!fitMvsusieArgs) |>
        .setFinemappingFitClass("mvsusie")
}

# Whether a joint run fits the full data as well as cross-validating, from the
# CrossValidationParam's fitFullData (TRUE when unset).
# @noRd
.jointFitsFullData <- function(cfg) {
    isTRUE(cfg$crossValidationParam$fitFullData %||% TRUE)
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
    cv <- cfg$crossValidationParam %||% list()
    cvFolds <- cv$folds %||% 0L
    if (cvFolds <= 1L) {
        return(NULL)
    }
    sp <- cfg$crossValidationParam$samplePartition %||% mvCv$samplePartition
    mvPriorCv <- .fmBuildMvsusiePriorCv(mvCv, mvFitParts, colnames(Ys), ddCut)
    cv <- .fmWeightsCv(
        Xc,
        Ys,
        "mvsusie",
        args$methodArgs,
        cvFolds,
        samplePartition = sp,
        coverage = cfg$credibleSetParam$coverage,
        verbose = verbose,
        mvPrior = mvPrior,
        mvPriorCv = mvPriorCv,
        BPPARAM = cv$BPPARAM %||% SerialParam(),
        seed = cfg$seed
    )
    .fmSliceCv(cv, "mvsusie")
}

# One mvsusie per-condition entry (NULL for a screened-out column), sliced at
# the condition's position in the fitted survivor set + its CV slice. Its
# credible sets' purity and its allele frequencies are the condition's own:
# read on the samples that observe it, which for a cross-study group are its
# study's alone.
# @noRd
.jointMvEntry <- function(
    i,
    fitted,
    Xc,
    Yc,
    keep,
    survivors,
    cfg,
    af = NULL,
    region = NULL
) {
    if (!keep[i]) {
        return(NULL)
    }
    r <- match(i, survivors)
    bare <- if (is.null(fitted$fit)) {
        .fmCvOnlyRow()
    } else {
        .fmPostprocessOne(
            fit = fitted$fit,
            method = "mvsusie",
            dataX = .jointConditionX(Xc, Yc, i),
            dataY = NULL,
            af = if (is.matrix(af)) af[, i] else af,
            region = region,
            conditionIdx = r,
            csInput = "X",
            credibleSetParam = cfg$credibleSetParam,
            fitRetention = cfg$fitRetention
        )
    }
    e <- if (is.null(fitted$cvM)) {
        bare
    } else {
        .fmAttachCv(bare, .fmSliceCvCondition(fitted$cvM, r))
    }
    e
}

# The genotypes of the samples observing condition `i`: all of them, except in
# a cross-study group, where each study's samples observe its own condition.
# @noRd
.jointConditionX <- function(Xc, Yc, i) {
    observed <- !is.na(Yc[, i])
    if (all(observed)) {
        return(Xc)
    }
    Xc[observed, , drop = FALSE]
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
            cfg = cfg,
            region = args$region
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
        args$mrmashPrior,
        key$study,
        key$trait,
        context = key$context,
        conditionNames = colnames(.jgZ(group))
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
        coverage = cfg$credibleSetParam$coverage
    ) |>
        c(compact(list(residual_variance = mvPrior$residualVariance)))
    fitArgs <- .fmMergeUserArgs(mvBaseArgs, args$methodArgs[["mvsusie"]])
    fitMvsusieRssArgs <- .splitMethodArgs(fitMvsusieRss, fitArgs)
    fit <- exec(fitMvsusieRss, !!!fitMvsusieRssArgs)
    .setFinemappingFitClass(fit, "mvsusie")
}

# One RSS per-condition FineMappingRow (csInput = "Xcorr").
# `credibleSetParam` / `fitRetention` come off `cfg`, which is what the caller
# supplies; they were also declared as formals with no defaults, so nothing
# could ever have passed them without erroring. `group` was unused too.
# @noRd
.jointRssEntry <- function(r, fit, ldMat, cfg, region = NULL) {
    .fmPostprocessOne(
        fit = fit,
        method = "mvsusie",
        dataX = ldMat,
        dataY = NULL,
        region = region,
        conditionIdx = r,
        csInput = "Xcorr",
        credibleSetParam = cfg$credibleSetParam,
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
# out-of-fold prediction matrix, the per-condition metric rows, and -- for
# mr.mash only, the one method whose fold fits are read -- the per-fold fits
# (named fold_<j>) fineMappingPipeline builds each fold's mvSuSiE prior from.
.jointTwasCvResult <- function(cv, token) {
    if (is.null(cv)) {
        return(NULL)
    }
    ffKey <- str_c(token, "_weights")
    foldFits <- if (.twasKeepsFoldFits(token) && !is.null(cv$foldFits)) {
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
    # The fine-mapping CV names its tables by the method's key (susie_inf),
    # which differs from the twas token (susieInf) for some methods.
    key <- .fmTwasMethodKey(token)
    base <- str_remove(
        names(fineMappingCv$prediction),
        "(_predicted|Predicted)$"
    )
    hit <- which(base == key)
    if (length(hit) == 0L) {
        return(NULL)
    }
    pBase <- str_remove(
        names(fineMappingCv$performance),
        "(_performance|Performance)$"
    )
    pHit <- which(pBase == key)
    list(
        samplePartition = fineMappingCv$samplePartition,
        predictions = fineMappingCv$prediction[[hit[[1L]]]],
        metrics = if (length(pHit)) {
            fineMappingCv$performance[[pHit[[1L]]]]
        } else {
            NULL
        }
    )
}

# Whether a TWAS method's cross-validation keeps its per-fold fits: only
# mr.mash's are read (fineMappingPipeline's per-fold mvSuSiE priors).
# @noRd
.twasKeepsFoldFits <- function(token) {
    identical(token, "mrmash")
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
        resolved <- .jointTwasMethodArgs(
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
        ma <- .jointTwasMultitrait(
            .jointTwasCoupled(resolved, token, cond, cfg),
            token,
            cond,
            cfg
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
            nCond,
            token
        )
        cvRes <- .jointTwasCv(Xc, Yc, wm, ma, full$W, args, cfg, token)
        # One per-condition entry: that condition's weight column + the shared
        # fit + its CV slice. fitFullData = FALSE -> CV-only entry. A
        # fine-mapping fit's columns are its own conditions', each found by
        # its row (args$fmColumns).
        map2(
            seq_len(nCond),
            args$fmColumns[[token]] %||% seq_len(nCond),
            .jointTwasEntry,
            full = full,
            cvRes = cvRes,
            stdz = stdz,
            cfg = cfg
        )
    }
)

# The first of `keys` that `lst` holds, or NULL.
# @noRd
.jointTwasSuppliedArgs <- function(lst, keys) {
    hit <- intersect(keys, names(lst))
    if (length(hit) == 0L) NULL else lst[[hit[[1L]]]]
}

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
    # A method's arguments sit under its .twasMethodMap list key
    # (bayes_r_weights), which differs from its engine key (bayesR_weights)
    # for every multi-word token, so both are looked up.
    keys <- unique(c(methodKey, .twasMethodMap[[token]]$fn))
    supplied <- .jointTwasSuppliedArgs(args$methodList, keys) %||%
        .jointTwasSuppliedArgs(args$methodArgs, keys)
    # FM-fit injection: an FM-derived token extracts its weights from the
    # precomputed fine-mapping fit rather than refitting.
    withFit <- .twasWithH2(
        .jointTwasInjectFit(supplied %||% list(), token, fittedModels),
        token,
        args$h2,
        "individual"
    )
    if (!isTRUE(cfg$estimatePi) || !is_in(token, c("bayesC", "bayesB"))) {
        return(withFit)
    }
    .jointTwasSpikeSlabPi(withFit, token, Xc, Yc, cond, cfg, stdz)
}

# The weight functions the univariate methods fit a group across studies
# with, under the penalty on the difference between the studies' effects
# (.lassoCoupledWeights() and its peers). mcp and scad have none: ncvreg takes
# no sparse design, and the stacked one is mostly zeros.
# @noRd
.TWAS_COUPLED_IMPL <- c(
    lasso = ".lassoCoupledWeights",
    enet = ".enetCoupledWeights",
    ridge = ".ridgeCoupledWeights",
    l0learn = ".l0learnCoupledWeights"
)

# A univariate method's arguments on a group across studies, in a run that
# couples them (`studyCoupling` given): they name its coupled implementation,
# which .resolveMethodFunction() reads, and carry the run's coupling
# strengths. Any other group's, and each study's own fit for the ensemble
# across studies (.twasStackEntries()), are left as they are.
# @noRd
.jointTwasCoupled <- function(ma, token, cond, cfg) {
    impl <- .TWAS_COUPLED_IMPL[token]
    if (
        is.na(impl) ||
            is.null(cfg$studyCoupling) ||
            n_distinct(cond$study) < 2L
    ) {
        return(ma)
    }
    `attr<-`(
        list_assign(ma, strengths = cfg$studyCoupling),
        "impl",
        unname(impl)
    )
}

# A Bayesian alphabet method's arguments on a group of two or more
# conditions, in a run that fits it jointly there (`multitraitTokens`: bayesC
# and bayesN on BGLR, under a jointSpecification): they name the whole-Y
# implementation that fits the group's conditions at once through BGLR's
# Multitrait() (.bglrJointWeights()), and say whether the group spans
# studies, whose residuals are independent. A one-condition group's, and
# each study's own fit for the ensemble across studies (.twasStackEntries()),
# are left as they are.
# @noRd
.jointTwasMultitrait <- function(ma, token, cond, cfg) {
    if (!is_in(token, cfg$multitraitTokens) || nrow(cond) < 2L) {
        return(ma)
    }
    `attr<-`(
        list_assign(
            ma,
            method = token,
            acrossStudies = n_distinct(cond$study) >= 2L
        ),
        "impl",
        ".bglrJointWeights"
    )
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
    fittedArg <- set_names(list(fittedModels[[token]]), adapter$fitArg)
    list_assign(ma, !!!fittedArg)
}

# Spike-and-slab pi from an internal mr.ash fit (self-contained per method).
# @noRd
.jointTwasSpikeSlabPi <- function(ma, token, Xc, Yc, cond, cfg, stdz) {
    # A pi the caller set wins, and then no mr.ash fit is needed.
    if (!is.null(ma$pi)) {
        return(ma)
    }
    mrA <- learnTwasWeights(
        Xc,
        Yc,
        weightMethods = list(mrash_weights = list()),
        studyName = as.character(cond$study[1L]),
        context = as.character(cond$context[1L]),
        trait = as.character(cond$trait[1L]),
        fitRetention = "slim",
        weightStandardized = stdz,
        verbose = 0,
        seed = cfg$seed
    )
    list_assign(ma, pi = as.numeric(estimateSparsity(mrA)))
}

# The retention level a joint fit runs at. The engine keeps its fit -- the
# joint layers downstream read it -- so "none" from the pipeline means "keep
# the least we can", not "keep nothing". A fine-mapping token's fit is the
# fineMappingResult's own, which its TWAS rows keep none of
# (.twasRowRetention()).
# @noRd
.jointRetentionLevel <- function(cfg, token = NULL) {
    level <- cfg$fitRetention %||% "slim"
    .twasRowRetention(token, if (identical(level, "none")) "slim" else level)
}

# Whether a joint fit keeps its copy of X: asked for, by a method whose fit has
# one. The run was refused at any level but "full", and a fine-mapping token,
# whose rows keep no fit, is not such a method.
# @noRd
.jointRetainX <- function(cfg, token) {
    isTRUE(cfg$retainX) && is_in(token, .twasTokensKeepingX())
}

# The mr.mash RSS fit of one summary-statistics joint group, with the caller's
# per-method options split across the wrapper's own formals and the engine's
# arguments exactly as the per-tuple path splits them (.twasWeightCallArgs()).
# @noRd
.jointMrmashRssWeights <- function(group, cfg, args) {
    callArgs <- .twasWeightCallArgs(
        "mrmashRssWeights",
        list(
            stat = list(
                z = .jgZ(group),
                n = .jgN(group),
                variantNames = .jgVariantIds(group)
            ),
            ldReference = .jgOneLdRef(group, "mrmash"),
            fitRetention = .jointRetentionLevel(cfg)
        ),
        .twasUserArgs(args$methodArgs, "mrmash")
    )
    exec(mrmashRssWeights, !!!callArgs)
}

# Full-data TWAS weight fit for a joint group. Returns list(W, fitParts, vids);
# W is NULL (a CV-only run) when the CrossValidationParam's fitFullData is
# FALSE.
# @noRd
.jointTwasFitFull <- function(
    Xc,
    Yc,
    cond,
    wm,
    fittedModels,
    cfg,
    stdz,
    nCond,
    token
) {
    if (!.jointFitsFullData(cfg)) {
        return(list(W = NULL, fitParts = NULL, vids = colnames(Xc)))
    }
    retention <- .jointRetentionLevel(cfg, token)
    verbose <- if (is.null(cfg$verbose)) 1 else cfg$verbose
    tw <- learnTwasWeights(
        Xc,
        Yc,
        weightMethods = wm,
        studyName = as.character(cond$study[1L]),
        context = as.character(cond$context[1L]),
        trait = as.character(cond$trait[1L]),
        fittedModels = fittedModels,
        fitRetention = retention,
        retainX = .jointRetainX(cfg, token),
        weightStandardized = stdz,
        verbose = verbose,
        seed = cfg$seed
    )
    base <- .twrRowParts(tw, 1L)
    vids <- .twrPartsVariantIds(base)
    raw <- weights(base)
    W <- if (is.matrix(raw)) {
        raw
    } else {
        matrix(raw, ncol = nCond, dimnames = list(vids, NULL))
    }
    list(W = W, fitParts = methodFits(base), vids = vids)
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

# Whether this token is barred from cross-validation: not requested, or its
# weights are all zero so there is nothing to validate. A fine-mapping method
# with no recoverable per-fold fits is an error rather than a silent skip.
# @noRd
.jointTwasCvBlocked <- function(W, cfg, token) {
    if (!.jointTwasCvRequested(cfg$crossValidationParam$weightMethods, token)) {
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
.jointTwasRunCv <- function(Xc, Yc, wm, args, cfg, cvFolds, token) {
    cv <- cfg$crossValidationParam %||% list()
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
        fitRetention = if (.twasKeepsFoldFits(token)) "slim" else "none",
        maxNumVariants = maxCv,
        BPPARAM = cv$BPPARAM %||% SerialParam(),
        dataDrivenPriorMatricesCv = args$dataDrivenPriorMatricesCv,
        verbose = cfg$verbose %||% 1,
        seed = cfg$seed
    )
}

# @noRd
.jointTwasCv <- function(Xc, Yc, wm, ma, W, args, cfg, token) {
    cv <- cfg$crossValidationParam %||% list()
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
        .jointTwasRunCv(Xc, Yc, wm, args, cfg, cvFolds, token),
        token
    )
}

# Condition `r`'s TwasWeightsRow (empty when there are no full-data weights,
# else column `column` of the weights + shared fit + CV slice). NULL for a
# condition the fit has no column of.
# @noRd
.jointTwasEntry <- function(r, column, full, cvRes, stdz, cfg) {
    if (is.na(column)) {
        return(NULL)
    }
    cvR <- if (!is.null(cvRes)) {
        .sliceTwasCvResultToCondition(cvRes, r)
    } else {
        NULL
    }
    if (is.null(full$W)) {
        return(TwasWeightsRow(
            variantIds = character(0),
            weights = NULL,
            cvResult = cvR,
            weightStandardized = stdz
        ))
    }
    TwasWeightsRow(
        variantIds = full$vids,
        weights = full$W[, column],
        methodFits = full$fitParts,
        cvResult = cvR,
        weightStandardized = stdz
    )
}

# (sumstats, twas) -> ONE matrix entry: mr.mash-rss fitted on the group, or a
# fine-mapping method's weights read off its fit of the group
# (.jointTwasRssWeights()). No sample folds.
setMethod(
    "fitJointGroup",
    signature("SumStatsJointGroup", "TwasJointPipeline"),
    function(group, pipeline, token, args) {
        cfg <- .jpConfig(pipeline)
        weights <- .jointTwasRssWeights(group, cfg, args, token)
        if (is.null(weights)) {
            return(NULL)
        }
        vids <- rownames(weights) %||% rownames(.jgZ(group))
        fitParts <- attr(weights, "fit")
        nCond <- ncol(.jgZ(group))
        # A single-condition fit comes back as a bare vector; the per-condition
        # split below reads it by column either way.
        wMatrix <- if (is.matrix(weights)) {
            weights
        } else {
            matrix(weights, ncol = nCond, dimnames = list(vids, NULL))
        }
        # One per-condition entry: that condition's weight column + the shared
        # fit. A fine-mapping fit's columns are its own conditions', each found
        # by its row (args$fmColumns).
        map2(
            seq_len(nCond),
            args$fmColumns[[token]] %||% seq_len(nCond),
            .jointColEntry,
            vids = vids,
            weights = wMatrix,
            fitParts = fitParts,
            cfg = cfg
        )
    }
)

# One summary-statistics group's weights for `token`: mr.mash-rss fitted on
# it, or a fine-mapping method's read off the group's own fine-mapping fit
# (.twasGroupFmArgs()), checked against the group's block. NULL when that
# fit has no part for this region, or, with a warning, when the fit fails:
# the other groups still run.
# @noRd
#' @importFrom rlang try_fetch
.jointTwasRssWeights <- function(group, cfg, args, token) {
    coupled <- .jointTwasIsCoupledRss(group, cfg, token)
    if (!coupled) {
        # Refused before the fit, which only warns: a group of studies on
        # their own references is no group for a fit taking one LD matrix.
        .jgOneLdRef(group, token)
    }
    try_fetch(
        if (coupled) {
            .jointTwasCoupledRss(group, cfg, args, token)
        } else {
            .jointTwasRssFit(group, cfg, args, token)
        },
        error = function(cnd) {
            warn(
                glue(
                    "twasWeightsPipeline: {token} failed for ",
                    "{.twasConditionsLabel(.jgConditions(group))}."
                ),
                parent = cnd
            )
            NULL
        }
    )
}

# @noRd
.jointTwasRssFit <- function(group, cfg, args, token) {
    adapter <- .twasFineMappingMethodAdapters[[token]]
    if (is.null(adapter)) {
        return(.jointMrmashRssWeights(group, cfg, args))
    }
    fit <- args$fittedModels[[token]]
    if (is.null(fit)) {
        return(NULL)
    }
    callArgs <- set_names(
        list(
            list(
                z = .jgZ(group),
                n = .jgN(group),
                variantNames = .jgVariantIds(group)
            ),
            .jgOneLdRef(group, token),
            fit,
            .jointRetentionLevel(cfg, token)
        ),
        c("stat", "ldReference", adapter$rssFitArg, "fitRetention")
    )
    exec(adapter$rssWeightFn, !!!callArgs)
}

# The penalties the summary-statistics methods fit a group across studies
# with, under the penalty on the difference between the studies' effects
# (.rssCoupledWeights()), from each method's own settings.
# @noRd
.TWAS_COUPLED_RSS <- c(
    lasso = "lasso",
    enet = "enet",
    ridge = "ridge",
    mcp = "MCP",
    scad = "SCAD",
    l0learn = "L0"
)

# The samplers fitted across studies, each by the function fitting it:
# PRS-CS as JointPRS fits ancestries, SDPR as SDPRX fits them.
# @noRd
.TWAS_COUPLED_MCMC <- c(
    prsCs = ".prsCsCoupledWeights",
    dprGibbs = ".sdprCoupledWeights"
)

# The summary-statistics methods fitted across studies: the penalized ones,
# under the penalty on the difference between the studies' effects, and the
# samplers.
# @noRd
.TWAS_COUPLED_SUMSTATS <- c(
    names(.TWAS_COUPLED_RSS),
    names(.TWAS_COUPLED_MCMC)
)

# The path options a PenalizedRssOptions holds, which a penalized method's
# settings carry flattened.
# @noRd
.RSS_PATH_OPTIONS <- c("lambda", "thr", "maxiter", "strategy")

# Whether `token` is fitted across `group`'s studies: a method fitted across
# studies, on a group of two or more studies, in a run that fits them
# (`studyCoupling` given).
# @noRd
.jointTwasIsCoupledRss <- function(group, cfg, token) {
    is_in(token, .TWAS_COUPLED_SUMSTATS) &&
        !is.null(cfg$studyCoupling) &&
        n_distinct(.jgConditions(group)$study) >= 2L
}

# A method's weights across `group`'s studies: each study's z-scores over the
# variants it carries, against its own LD reference, a sampler's by its chain
# across studies (.TWAS_COUPLED_MCMC) and a penalized method's under its own
# penalty, shrinkage grid and path options (.rssCoupledWeights()).
# @noRd
.jointTwasCoupledRss <- function(group, cfg, args, token) {
    Z <- .jgZ(group)
    .jointCoupledFit(
        map(seq_len(ncol(Z)), .jointCoupledStat, Z = Z, N = .jgN(group)),
        group,
        cfg,
        args,
        token,
        .jointRetentionLevel(cfg, token)
    )
}

# `token`'s weights across the studies whose statistics are `stats`, each on
# its condition's LD reference in `group`, at retention `retention`.
# @noRd
.jointCoupledFit <- function(stats, group, cfg, args, token, retention) {
    user <- .twasUserArgs(args$methodArgs, token)
    label <- glue("{token} across studies")
    Z <- .jgZ(group)
    refs <- map(seq_len(ncol(Z)), .jgLdRefAt, g = group)
    if (is_in(token, names(.TWAS_COUPLED_MCMC))) {
        return(exec(
            .TWAS_COUPLED_MCMC[[token]],
            stats,
            refs,
            rownames(Z),
            user,
            retention,
            label
        ))
    }
    .coupledRssRefuseSelection(user$selection, label)
    .rssCoupledWeights(
        stats,
        refs,
        rownames(Z),
        .coupledRssSpec(token, user, label),
        user$s %||% c(0.2, 0.5, 0.9, 1.0),
        cfg$studyCoupling,
        user[intersect(names(user), .RSS_PATH_OPTIONS)],
        retention,
        label
    )
}

# Study `k`'s statistics: its z-scores and sample size over the variants it
# carries (those not NA in its column).
# @noRd
.jointCoupledStat <- function(k, Z, N) {
    held <- !is.na(Z[, k])
    list(
        z = unname(Z[held, k]),
        n = rep(N[[k]], sum(held)),
        variantNames = rownames(Z)[held]
    )
}

# The penalty a method fits across studies with: its own, from its settings
# (enet's alpha, mcp's and scad's gamma and alpha). l0learn fits the L0
# penalty alone, as the individual-level coupled fit does.
# @noRd
.coupledRssSpec <- function(token, user, label) {
    penalty <- .TWAS_COUPLED_RSS[[token]]
    if (identical(penalty, "L0") && !is.null(user$penalty)) {
        .coupledRequireL0(user$penalty)
    }
    switch(
        penalty,
        lasso = .rssPenalty("lasso"),
        enet = .rssPenalty("lasso", alpha = user$alpha %||% 0.5),
        ridge = .rssPenalty("lasso", alpha = 0),
        L0 = .rssPenalty("L0"),
        .rssPenalty(penalty, alpha = user$alpha %||% 1, gamma = user$gamma)
    )
}

# A fit across studies picks each study's grid point on its own held-out
# statistics, the rule the per-study fits default to.
# @noRd
.coupledRssRefuseSelection <- function(selection, label) {
    if (is.null(selection) || identical(selection, "heldOut")) {
        return(invisible(NULL))
    }
    abort(c(
        glue(
            "{label}: selection = \"{selection}\" does not apply; each ",
            "study picks its grid point on its own held-out statistics."
        ),
        i = str_c(
            "Leave `selection` unset, or pass studyCoupling = NULL to fit ",
            "each study on its own."
        )
    ))
}

# ---- result construction (construct) ----------------------------------------

# Both pipelines assemble identically-shaped joint rows; only the result
# collection differs (the axis-3 divergence the markers encode). Only the joint*
# columns for axes that actually vary are attached.
# Fold the accumulator's per-row records into one collection: build each record
# into a 1-row collection via `constructor`, then union them with the generic
# .rbindCollections (which aligns optional columns and sets the ldReference
# slot). The only per-row transforms are: wrap `entry` in a list, and drop an NA
# joint* value so the 1-row collection omits that column (the union re-adds it,
# padded, only when some row is joint). Every other field --
# study/context/trait/method, region, and any future column -- flows through
# untouched.
.buildJointResult <- function(constructor, records) {
    if (length(records) == 0L) {
        return(NULL)
    }
    parts <- map(records, .jointBuildRecordPart, constructor = constructor)
    .rbindCollections(parts)
}

setMethod("construct", "FmJointPipeline", function(pipeline, records) {
    .buildJointResult(QtlFineMappingResult, records)
})

setMethod("construct", "TwasJointPipeline", function(pipeline, records) {
    .buildJointResult(TwasWeights, records)
})

# ---- enumerators (pattern x dataForm -> list<JointGroup>) --------------------

# cross-context / individual: one group per scoped trait present in >= 2 scoped
# contexts (the conditions are those (study, context, trait) rows).
.enumCrossContextIndividual <- function(data, scope, args = list()) {
    studyName <- studyName(data)
    if (!is_in(studyName, scope$studies)) {
        return(list())
    }
    scopedContexts <- scope$contexts[[studyName]]
    scopedTraits <- scope$traits[[studyName]]
    if (length(scopedContexts) < 2L) {
        return(list())
    }
    verbose <- if (is.null(args$verbose)) 1 else args$verbose
    compact(map(
        scopedTraits,
        .deferCrossContext,
        data = data,
        studyName = studyName,
        scopedContexts = scopedContexts,
        args = args,
        verbose = verbose
    ))
}

# One trait's deferred cross-context unit, or NULL when it cannot be fit. Its
# phenotypes are built here: the samples a fit is on decide what it can share an
# X with, and whether it needs one at all.
# @noRd
.deferCrossContext <- function(
    tid,
    data,
    studyName,
    scopedContexts,
    args,
    verbose
) {
    planned <- .crossContextPlan(
        data,
        tid,
        scopedContexts,
        verbose,
        "jointCrossContext",
        args$residualizationParam,
        args$phenotypeFilter
    )
    if (is.null(planned)) {
        return(NULL)
    }
    .deferGroups(
        tid,
        .enumCrossContextGroupFor,
        xTraits = tid,
        xContexts = planned$perTraitContexts,
        samples = planned$fitSamples,
        studyName = studyName,
        planned = planned
    )
}

# One trait's cross-context group on `X`, or NULL when it has no usable (X, Y).
# @noRd
.enumCrossContextGroupFor <- function(tid, studyName, planned, X) {
    xy <- .crossContextAssemble(
        X,
        planned$Yres,
        planned$perTraitContexts,
        planned$fitSamples,
        planned$report
    )
    if (is.null(xy)) {
        return(NULL)
    }
    new(
        "IndividualJointGroup",
        conditions = tibble(
            study = studyName,
            context = xy$perTraitContexts,
            trait = tid
        ),
        X = xy$X,
        Y = xy$Y
    )
}

# cross-context / sumstats.
.enumCrossContextSumstats <- function(data, scope, args = list()) {
    studyCol <- as.character(data$study)
    contextCol <- as.character(data$context)
    traitCol <- as.character(data$trait)
    .jeConcat(map(
        scope$studies,
        .enumCrossContextSumstatsForStudy,
        data = data,
        scope = scope,
        args = args,
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

# The one LD reference a summary-statistics joint group's rows share. A joint
# fit reads a single LD matrix for all its conditions, so rows on different
# references cannot be fit together; that is refused, naming them, rather
# than fitting some rows against LD they were never harmonized to.
# @noRd
.jointRowsLdRef <- function(data, rows, label) {
    keys <- .ldRefKeys(data)[rows]
    if (n_distinct(keys) > 1L) {
        tuples <- str_c(
            data$study[rows],
            data$context[rows],
            data$trait[rows],
            sep = ":"
        )
        abort(c(
            glue(
                "{label}: a joint fit uses one LD matrix, but its rows use ",
                "different LD references."
            ),
            i = str_c(
                "Rows and their references: ",
                str_flatten(str_c(tuples, " -> ", keys), "; ")
            ),
            i = "Fit them separately, or give them one LD reference."
        ))
    }
    .rowLdRef(data, rows[[1L]])
}

# One study's cross-context groups (none when it has fewer than two contexts).
# @noRd
.enumCrossContextSumstatsForStudy <- function(
    s,
    data,
    scope,
    args,
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
    label <- "jointCrossContext (QtlSumStats)"
    ldReference <- .jointRowsLdRef(data, tupleRows, label)
    ctxNames <- contextCol[tupleRows]
    jz <- .buildJointSumstatZMatrix(
        data,
        tupleRows,
        ctxNames,
        errorLabel = label,
        ldReference = ldReference,
        cutoffs = args$cutoffs
    )
    new(
        "SumStatsJointGroup",
        conditions = tibble(study = s, context = ctxNames, trait = tid),
        Z = jz$Z,
        ldReference = ldReference,
        N = jz$nVec
    )
}

# cross-trait / individual: one group per scoped context with >= 2 scoped
# traits.
.enumCrossTraitIndividual <- function(data, scope, args = list()) {
    studyName <- studyName(data)
    if (!is_in(studyName, scope$studies)) {
        return(list())
    }
    scopedContexts <- scope$contexts[[studyName]]
    scopedTraits <- scope$traits[[studyName]]
    verbose <- if (is.null(args$verbose)) 1 else args$verbose
    compact(map(
        scopedContexts,
        .deferCrossTrait,
        data = data,
        studyName = studyName,
        scopedTraits = scopedTraits,
        args = args,
        verbose = verbose
    ))
}

# One context's deferred cross-trait unit, or NULL when it cannot be fit. Its
# window spans the scoped traits the context carries.
# @noRd
.deferCrossTrait <- function(
    cx,
    data,
    studyName,
    scopedTraits,
    args,
    verbose
) {
    planned <- .crossTraitPlan(
        data,
        cx,
        scopedTraits,
        studyName,
        verbose,
        "jointCrossTrait",
        args$residualizationParam,
        args$phenotypeFilter
    )
    if (is.null(planned)) {
        return(NULL)
    }
    .deferGroups(
        cx,
        .enumCrossTraitGroupFor,
        xTraits = planned$traitsHere,
        xContexts = cx,
        samples = planned$fitSamples,
        studyName = studyName,
        planned = planned
    )
}

# One context's cross-trait group on `X`, or NULL when it has no usable (X, Y).
# @noRd
.enumCrossTraitGroupFor <- function(cx, studyName, planned, X) {
    xy <- .crossTraitAssemble(X, planned$allY, planned$traitsHere, planned$se)
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
            study = studyName,
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
    studyCol <- as.character(data$study)
    contextCol <- as.character(data$context)
    traitCol <- as.character(data$trait)
    .jeConcat(map(
        scope$studies,
        .enumCrossTraitSumstatsForStudy,
        data = data,
        scope = scope,
        args = args,
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
    label <- "jointCrossTrait (QtlSumStats)"
    ldReference <- .jointRowsLdRef(data, tupleRows, label)
    trNames <- traitCol[tupleRows]
    jz <- .buildJointSumstatZMatrix(
        data,
        tupleRows,
        trNames,
        errorLabel = label,
        ldReference = ldReference,
        cutoffs = args$cutoffs
    )
    new(
        "SumStatsJointGroup",
        conditions = tibble(study = s, context = cx, trait = trNames),
        Z = jz$Z,
        ldReference = ldReference,
        N = jz$nVec
    )
}

# cross-study / sumstats. One group per (context, trait) present in >= 2
# scoped studies; the study axis varies -> "joint" + jointStudies. The
# individual-level form is .enumCrossStudyIndividual().
.enumCrossStudySumstats <- function(data, scope, args = list()) {
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
        allTrs = allTrs
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
    allTrs
) {
    compact(map(
        allTrs,
        .enumCrossStudyGroupFor,
        data = data,
        scope = scope,
        args = args,
        cols = cols,
        cx = cx
    ))
}

# `map()` hands the trait first; .enumCrossStudyGroup takes it sixth.
# @noRd
.enumCrossStudyGroupFor <- function(tid, data, scope, args, cols, cx) {
    .enumCrossStudyGroup(data, scope, args, cols, cx, tid)
}

# One cross-study joint group for a (context, trait) cell, or NULL when fewer
# than two studies survive the scope filter -- a "joint" fit over one study is
# just that study, so the cell contributes nothing. Studies on one LD
# reference share its SNP order; studies each on their own reference (each
# ancestry its own panel) hold their z-scores over the union of their
# variants (.buildJointSumstatZUnion()), and the group one reference per
# condition, which only the fits across studies under the coupling read.
# @noRd
.enumCrossStudyGroup <- function(data, scope, args, cols, cx, tid) {
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
    jz <- .crossStudyGroupZ(data, tupleRows, cols, args$cutoffs)
    new(
        "SumStatsJointGroup",
        conditions = tibble(
            study = cols$study[tupleRows],
            context = cx,
            trait = tid
        ),
        Z = jz$Z,
        ldReference = jz$ldReference,
        N = jz$nVec
    )
}

# The z-scores, sample sizes and LD of a cross-study group's rows: one
# reference and its SNP order when the rows share one, else each row's own
# reference and the union of their variants (.buildJointSumstatZUnion()).
# @noRd
.crossStudyGroupZ <- function(data, tupleRows, cols, cutoffs) {
    if (n_distinct(.ldRefKeys(data)[tupleRows]) > 1L) {
        refs <- map(tupleRows, .rowLdRef, x = data)
        jz <- .buildJointSumstatZUnion(
            tupleRows,
            cols,
            refs,
            list(data = data, cutoffs = cutoffs, label = "jointCrossStudy")
        )
        return(list_assign(jz, ldReference = refs))
    }
    ldReference <- .rowLdRef(data, tupleRows[[1L]])
    jz <- .buildJointSumstatZMatrix(
        data,
        tupleRows,
        cols$study[tupleRows],
        errorLabel = "jointCrossStudy",
        ldReference = ldReference,
        cutoffs = cutoffs
    )
    list_assign(jz, ldReference = ldReference)
}

# cross-study / individual: a MultiStudyQtlDataset's individual-level studies
# fitted together, one group per (context, trait) that two or more of its
# scoped studies carry. Their samples are disjoint, so a group stacks them:
# each study's samples in rows of their own, the union of the studies' variants
# in the columns, and each study's phenotype observed on its own rows alone
# (block-missing Y). Each study keeps its own residualization, filters and
# scaling. The study axis varies -> "joint" + jointStudies.
.enumCrossStudyIndividual <- function(data, scope, args = list()) {
    studies <- intersect(names(data), scope$studies)
    if (length(studies) < 2L) {
        return(list())
    }
    .jeConcat(map(
        unique(list_c(scope$contexts[studies])),
        .enumCrossStudyIndividualForContext,
        data = data,
        scope = scope,
        studies = studies,
        args = args,
        phenotypeFilter = .jointPhenotypeFilter(args)
    ))
}

# One context's cross-study units, one per trait two or more studies carry in
# it.
# @noRd
.enumCrossStudyIndividualForContext <- function(
    cx,
    data,
    scope,
    studies,
    args,
    phenotypeFilter
) {
    carried <- set_names(
        map(studies, .xsStudyTraits, cx = cx, data = data, scope = scope),
        studies
    )
    compact(map(
        unique(list_c(carried)),
        .deferCrossStudy,
        cx = cx,
        carried = carried,
        data = data,
        args = args,
        phenotypeFilter = phenotypeFilter
    ))
}

# The scoped, selected traits study `st` carries in context `cx`; none when the
# context is not one of its scoped ones.
# @noRd
.xsStudyTraits <- function(st, cx, data, scope) {
    if (!is_in(cx, scope$contexts[[st]])) {
        return(character(0))
    }
    .enumContextTraits(cx, data[[st]], scope, st)
}

# One (context, trait) cross-study unit, or NULL when fewer than two studies
# carry the trait with samples to fit. The phenotypes are read here, before any
# genotype, so the studies' samples are known to be disjoint before a block is
# read.
# @noRd
.deferCrossStudy <- function(tid, cx, carried, data, args, phenotypeFilter) {
    plans <- compact(map(
        names(keep(carried, .xsCarries, tid = tid)),
        .xsStudyPlan,
        tid = tid,
        cx = cx,
        data = data,
        args = args,
        phenotypeFilter = phenotypeFilter
    ))
    if (length(plans) < 2L) {
        return(NULL)
    }
    .xsRequireDisjointSamples(plans, cx, tid)
    .deferSelfReading(
        tid,
        .enumCrossStudyIndividualGroupFor,
        cx = cx,
        plans = plans,
        args = args
    )
}

# @noRd
.xsCarries <- function(traits, tid) {
    is_in(tid, traits)
}

# A cross-study unit's conditions, known before any genotype is read: one per
# study planned into it (.deferCrossStudy()).
# @noRd
.xsUnitConditions <- function(unit) {
    tibble(
        study = map_chr(unit$args$plans, "study"),
        context = unit$args$cx,
        trait = unit$args[[1L]]
    )
}

# One study's part of a cross-study unit, before any genotype is read: its
# phenotype, residualized as the study's own fits residualize it, and the
# genotypes it is to be read on. NULL when fewer than two samples have one.
# @noRd
.xsStudyPlan <- function(st, tid, cx, data, args, phenotypeFilter) {
    qd <- data[[st]]
    Y <- .fmResidPheno(
        qd,
        contexts = cx,
        traitId = tid,
        phenotypeFilter = phenotypeFilter,
        residualizationParam = args$residualizationParam
    )
    if (nrow(Y) < 2L) {
        return(NULL)
    }
    list(
        study = st,
        dataset = qd,
        Y = Y,
        spec = .qtlFitSpec(
            .fmGenoWindow(qd, tid, args$region, args$cisWindow),
            cx,
            rownames(Y)
        )
    )
}

# The studies of a cross-study group share no sample: its rows are each one
# study's, with residuals independent across studies, so a sample in two of
# them would be counted twice as two unrelated ones.
# @noRd
.xsRequireDisjointSamples <- function(plans, cx, tid) {
    samples <- list_c(map(plans, .xsPlanSamples))
    shared <- unique(samples[duplicated(samples)])
    if (length(shared) == 0L) {
        return(invisible(NULL))
    }
    abort(c(
        glue(
            "jointCrossStudy: the studies fitted for context '{cx}', trait ",
            "'{tid}' share {length(shared)} sample(s)."
        ),
        i = glue("Studies: {str_flatten(map_chr(plans, 'study'), ', ')}."),
        i = str_c(
            "A cross-study fit stacks its studies' samples as rows of their ",
            "own, so it needs studies on disjoint samples."
        )
    ))
}

# @noRd
.xsPlanSamples <- function(plan) {
    rownames(plan$Y)
}

# One (context, trait) cross-study group, each study's genotypes read on a plan
# of its own (.xsStudyPiece()) and the studies stacked (.xsStack()); NULL when
# fewer than two studies keep samples with genotypes.
# @noRd
.enumCrossStudyIndividualGroupFor <- function(tid, cx, plans, args) {
    pieces <- compact(map(plans, .xsStudyPiece, args = args))
    if (length(pieces) < 2L) {
        return(NULL)
    }
    stacked <- .xsStack(pieces)
    new(
        "IndividualJointGroup",
        conditions = tibble(
            study = map_chr(pieces, "study"),
            context = cx,
            trait = tid
        ),
        X = stacked$X,
        Y = stacked$Y,
        af = stacked$af
    )
}

# One study's residualized genotypes over its window, on its phenotype's
# samples, with their allele frequencies: read on the study's own dataset under
# the run's genotype filter and residualization, as its own fits read them.
# @noRd
.xsStudyPiece <- function(plan, args) {
    fits <- .qtlRunWork(
        .qtlWork(list(plan), list(plan$spec), .xsPieceOf),
        plan$dataset,
        args$genotypeFilter,
        args$residualizationParam,
        withAf = TRUE
    )
    fits[[1L]]
}

# A study's piece of its block: the samples the block shares with the
# phenotype, or NULL when fewer than two, or no variant, remain.
# @noRd
.xsPieceOf <- function(plan, block) {
    X <- block$X
    common <- intersect(rownames(X), rownames(plan$Y))
    if (length(common) < 2L || ncol(X) == 0L) {
        return(NULL)
    }
    list(
        study = plan$study,
        X = .qtlRowsOf(X, common),
        af = block$af %||% rep(NA_real_, ncol(X)),
        Y = plan$Y[common, , drop = FALSE]
    )
}

# The studies' pieces as one design. Rows: each study's samples. Columns: the
# union of the studies' variants, matched by allele (matchVariants()), in
# genomic order. Where a study spells a variant's alleles the other way round,
# its genotypes are negated and its frequencies complemented into the union's
# frame. A variant a study lacks is a column of zeros in its rows -- its
# residualized genotypes have mean zero, so the column carries nothing there --
# and NA among its frequencies. Y holds each study's phenotype on its own rows
# and NA elsewhere; af is variants x studies.
# @noRd
.xsStack <- function(pieces) {
    ids <- .xsVariantUnion(map(pieces, .xsPieceIds))
    aligned <- map(pieces, .xsAlignPiece, ids = ids)
    studies <- map_chr(pieces, "study")
    xs <- map(aligned, "X")
    ys <- map(seq_along(pieces), .xsStudyY, pieces = pieces, studies = studies)
    afs <- map(aligned, "af")
    list(
        X = exec(rbind, !!!xs),
        Y = exec(rbind, !!!ys),
        af = `dimnames<-`(exec(cbind, !!!afs), list(ids, studies))
    )
}

# @noRd
.xsPieceIds <- function(piece) {
    colnames(piece$X)
}

# The union of the studies' variant ids, in genomic order: each study adds the
# variants no earlier one holds, matched by allele, in its own spelling.
# @noRd
.xsVariantUnion <- function(idSets) {
    ids <- reduce(idSets, .xsAddVariants)
    parsed <- parseVariantId(ids)
    ids[order(parsed$chrom, parsed$pos)]
}

# @noRd
.xsAddVariants <- function(held, ids) {
    m <- matchVariants(ids, held, removeStrandAmbiguous = FALSE)
    c(held, ids[!is_in(seq_along(ids), m$idxA)])
}

# A study's piece in the union's columns: its genotypes signed into the union's
# allele frame and zero where it lacks the variant; its frequencies likewise,
# complemented where flipped and NA where it lacks the variant.
# @noRd
.xsAlignPiece <- function(piece, ids) {
    m <- matchVariants(colnames(piece$X), ids, removeStrandAmbiguous = FALSE)
    at <- match(seq_along(ids), m$idxB)
    src <- coalesce(m$idxA[at], 1L)
    sign <- coalesce(m$sign[at], 0)
    af <- piece$af[src]
    X <- sweep(piece$X[, src, drop = FALSE], 2L, sign, "*")
    list(
        X = `dimnames<-`(X, list(rownames(piece$X), ids)),
        af = if_else(is.na(at), NA_real_, if_else(sign < 0, 1 - af, af))
    )
}

# Study `k`'s rows of the stacked phenotype: its own phenotype in its column,
# NA in the other studies'.
# @noRd
.xsStudyY <- function(k, pieces, studies) {
    y <- pieces[[k]]$Y[, 1L]
    cols <- map(seq_along(studies), .xsYColumn, y = y, k = k)
    `dimnames<-`(exec(cbind, !!!cols), list(rownames(pieces[[k]]$Y), studies))
}

# @noRd
.xsYColumn <- function(j, y, k) {
    if (j == k) unname(y) else rep(NA_real_, length(y))
}

# composed / individual: ONE group joining every scoped (context, trait) tuple
# for the study. Both context and trait vary across conditions, so both collapse
# to "joint" (the conditions model handles multi-varying-axis uniformly; if the
# tuples happen to share a context it degrades to cross-trait, and vice versa).
.enumComposedIndividual <- function(data, scope, args = list()) {
    studyName <- studyName(data)
    if (!is_in(studyName, scope$studies)) {
        return(list())
    }
    verbose <- if (is.null(args$verbose)) 1 else args$verbose
    planned <- .composedPlan(
        data,
        scope,
        studyName,
        verbose,
        "composed",
        args$residualizationParam,
        args$phenotypeFilter
    )
    if (is.null(planned)) {
        return(list())
    }
    list(.deferGroups(
        studyName,
        .enumComposedIndividualGroup,
        xTraits = planned$allTraits,
        xContexts = planned$allContexts,
        samples = planned$fitSamples,
        planned = planned
    ))
}

# The composed group on `X`, or NULL when it has no usable (X, Y).
# @noRd
.enumComposedIndividualGroup <- function(studyName, planned, X) {
    xy <- .composedAlignXy(
        list(X = X, YresList = planned$YresList),
        planned$tuples
    )
    if (is.null(xy)) {
        return(NULL)
    }
    # Conditions follow the fitted Y columns ("context:trait"), so dropped
    # tuples don't desync conditions from Y. Split on the first ":" (contexts
    # are simple labels; trait ids may themselves contain ":").
    labs <- colnames(xy$Y)
    conds <- tibble(
        study = studyName,
        context = str_remove(labs, ":.*$"),
        trait = str_remove(labs, "^[^:]*:")
    )
    new("IndividualJointGroup", conditions = conds, X = xy$X, Y = xy$Y)
}

# The run's resolved phenotype QC, or the univariate fits' own default -- drop
# a sample with a missing value, keep outliers -- when the engine record
# carries none.
# @noRd
.jointPhenotypeFilter <- function(args) {
    args$phenotypeFilter %||%
        .qtlResolvePhenotypeFilter(
            PhenotypeFilterParam(),
            naDefault = "drop",
            naAllowed = c("drop", "impute"),
            caller = "the joint engine"
        )
}

# univariate / individual: one 1-condition group per (study, context, trait) in
# scope -- the per-(context, trait) iteration expressed as engine groups, so
# univariate methods (lasso / enet / susie / ...) flow through the SAME per-
# method fitter + ensemble layer as the joint ones (minGroup = 1).
.enumUnivariateIndividual <- function(data, scope, args = list()) {
    studyName <- studyName(data)
    if (!is_in(studyName, scope$studies)) {
        return(list())
    }
    phenotypeFilter <- .jointPhenotypeFilter(args)
    .jeConcat(map(
        scope$contexts[[studyName]],
        .enumUnivariateForContext,
        data = data,
        scope = scope,
        args = args,
        studyName = studyName,
        phenotypeFilter = phenotypeFilter
    ))
}

# The traits of context `cx` that a study's scope selects.
# @noRd
.enumContextTraits <- function(cx, data, scope, studyName) {
    se <- .qtlContextSe(data, context = cx)
    keep(
        intersect(scope$traits[[studyName]], rownames(se)),
        .traitSelectedAt,
        traitId = scope$traitSelection,
        s = studyName,
        cx = cx
    )
}

# One context's univariate units, one per trait it carries. Their phenotypes are
# built here, before any genotype is read: the samples a fit is on are what
# decide which traits can share a residualized X.
# @noRd
.enumUnivariateForContext <- function(
    cx,
    data,
    scope,
    args,
    studyName,
    phenotypeFilter
) {
    traits <- .enumContextTraits(cx, data, scope, studyName)
    map(
        traits,
        .deferUnivariate,
        data = data,
        args = args,
        studyName = studyName,
        cx = cx,
        phenotypeFilter = phenotypeFilter
    )
}

# One trait's deferred univariate unit.
# @noRd
.deferUnivariate <- function(
    tid,
    data,
    args,
    studyName,
    cx,
    phenotypeFilter
) {
    Y <- .fmResidPheno(
        data,
        contexts = cx,
        traitId = tid,
        phenotypeFilter = phenotypeFilter,
        residualizationParam = args$residualizationParam
    )
    .deferGroups(
        tid,
        .enumUnivariateGroupFor,
        xTraits = tid,
        xContexts = cx,
        samples = rownames(Y),
        studyName = studyName,
        cx = cx,
        Y = Y
    )
}

# One (context, trait) univariate group on `X`, or NULL when fewer than two
# samples are shared between its genotypes and its phenotype.
# @noRd
.enumUnivariateGroupFor <- function(tid, studyName, cx, Y, X) {
    common <- intersect(rownames(X), rownames(Y))
    if (length(common) < 2L) {
        return(NULL)
    }
    new(
        "IndividualJointGroup",
        conditions = tibble(study = studyName, context = cx, trait = tid),
        X = .qtlRowsOf(X, common),
        Y = Y[common, , drop = FALSE]
    )
}

# composed / sumstats: general N-axis joint. `args$axes` (subset of study /
# context / trait) names the collapsed axes; rows split by the complement
# (fixed) axes form one group each. Reuses .enumerateComposedSumstatGroups.
.enumComposedSumstats <- function(data, scope, args = list()) {
    axes <- args$axes %||% c("context", "trait")
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
        args = args
    )))
}

# One composed cell's group, or NULL when it spans fewer than two tuples --
# a "joint" fit over one tuple is just that tuple.
# @noRd
.enumComposedGroupFor <- function(gIdx, data, gi, args) {
    if (length(gIdx) < 2L) {
        return(NULL)
    }
    label <- "composed (QtlSumStats)"
    ldReference <- .jointRowsLdRef(data, gIdx, label)
    jz <- .buildJointSumstatZMatrix(
        data,
        gIdx,
        map_chr(gIdx, .jointGroupColLabel, gi = gi),
        errorLabel = label,
        ldReference = ldReference,
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
        ldReference = ldReference,
        N = jz$nVec
    )
}

# ---- engine -----------------------------------------------------------------

# Twas per-group args: resolve the group's fine-mapping fits + CV (keyed on its
# first condition -- the joint fit is shared across conditions) and fix ONE
# shared fold partition (so every method's out-of-fold CV predictions align for
# the ensemble layer). Returns `args` unchanged for fine-mapping pipelines.
.twasGroupArgs <- function(g, pipeline, args, tokens) {
    if (!is(pipeline, "TwasJointPipeline")) {
        return(args)
    }
    cfg <- .jpConfig(pipeline)
    cond <- .jgConditions(g)
    groupFmArgs <- .twasGroupFmArgs(args, cond, tokens)
    out <- list_assign(args, !!!groupFmArgs)
    cv <- cfg$crossValidationParam %||% list()
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

# The group's fine-mapping fits and CV, from the fine-mapping fit over exactly
# its conditions (.twasRequireFmGroup()): the fits read off its first
# condition's row (every condition's row of a joint fit keys the one fit), and
# each condition's column of each fit (.twasGroupFmColumns()). A group of one
# condition reads that condition's own fits. Empty when the caller supplied no
# fine-mapping result.
# @noRd
.twasGroupFmArgs <- function(args, cond, tokens) {
    fmRes <- args$fineMappingResult
    if (is.null(fmRes)) {
        return(list())
    }
    joint <- .jointKeyOfConditions(cond)
    .twasRequireFmGroup(fmRes, cond, joint, tokens)
    nR <- if (is.null(args$nRegions)) 1L else args$nRegions
    bi <- if (is.null(args$regionIndex)) 1L else args$regionIndex
    af <- .twasFineMappingFits(
        fmRes,
        studyName = as.character(cond$study[[1L]]),
        context = as.character(cond$context[[1L]]),
        trait = as.character(cond$trait[[1L]]),
        joint = joint
    )
    list(
        fittedModels = .twasFitsForRegion(af, bi, nR),
        fineMappingCv = .twasGroupFmCv(fmRes, cond, joint),
        fmColumns = .twasGroupFmColumns(fmRes, cond, joint, names(af))
    )
}

# Fine-mapping and TWAS group conditions alike, so a group reads a joint
# fine-mapping method only from that method's fit over exactly its conditions.
# Fits of its conditions over others are a divergence between the two runs,
# not a fit to borrow; with none at all, there is nothing to read.
# @noRd
.twasRequireFmGroup <- function(fmRes, cond, joint, tokens) {
    if (nrow(cond) < 2L) {
        return(invisible(NULL))
    }
    walk(
        keep(
            intersect(tokens, names(.twasFineMappingMethodAdapters)),
            .twasIsMultivariateToken
        ),
        .twasRequireFmGroupFor,
        fmRes = fmRes,
        cond = cond,
        joint = joint
    )
}

# @noRd
.twasRequireFmGroupFor <- function(token, fmRes, cond, joint) {
    own <- .twasFmConditionRows(fmRes, token, cond, joint)
    if (length(own) > 0L) {
        return(invisible(NULL))
    }
    others <- .twasFmConditionRows(fmRes, token, cond, NULL)
    groups <- unique(map_chr(others, .twasFmGroupLabelAt, fmRes = fmRes))
    found <- if (length(groups) == 0L) {
        "it holds none of these conditions"
    } else {
        str_c(
            "its fits of these conditions are over ",
            str_flatten(groups, "; ")
        )
    }
    abort(c(
        glue(
            "twasWeightsPipeline: fineMappingResult has no {token} fit over ",
            "exactly {.twasJointKeyLabel(joint)} ",
            "({.twasConditionsLabel(cond)}), the group this run fits them ",
            "in; {found}."
        ),
        i = str_c(
            "TWAS reads a joint fine-mapping fit only for the group it was ",
            "fitted over: run both pipelines with the same jointSpecification ",
            "and the same contexts / traitId selection."
        )
    ))
}

# The rows `token` holds for the conditions `cond`, those of joint group
# `joint` (any group when NULL).
# @noRd
.twasFmConditionRows <- function(fmRes, token, cond, joint) {
    list_c(
        pmap(
            list(
                studyName = as.character(cond$study),
                context = as.character(cond$context),
                trait = as.character(cond$trait)
            ),
            .twasFmTokenRows,
            x = fmRes,
            token = token,
            joint = joint
        ),
        ptype = integer()
    )
}

# Row `i`'s group as a label (.twasJointKeyLabel()).
# @noRd
.twasFmGroupLabelAt <- function(i, fmRes) {
    .twasJointKeyLabel(.jointKeyAt(fmRes, i))
}

# A joint key as text: its set axes, or "no joint group".
# @noRd
.twasJointKeyLabel <- function(joint) {
    set <- discard(joint, is.na)
    if (length(set) == 0L) {
        return("no joint group")
    }
    str_flatten(str_c(names(set), " = ", unlist(set)), ", ")
}

# The first conditions of `cond` as "(study, context, trait)" text.
# @noRd
.twasConditionsLabel <- function(cond) {
    shown <- head(str_c(cond$study, ", ", cond$context, ", ", cond$trait), 3L)
    more <- if (nrow(cond) > 3L) {
        str_c(" and ", nrow(cond) - 3L, " more")
    } else {
        ""
    }
    str_c(str_flatten(str_c("(", shown, ")"), ", "), more)
}

# Each condition's column of the group's fine-mapping fit, per method among
# `tokens`: its row's position in the fit's condition axis, NA for a
# condition the fit has no row of. NULL for a group of one condition, whose
# fits have one column.
# @noRd
.twasGroupFmColumns <- function(fmRes, cond, joint, tokens) {
    if (nrow(cond) < 2L) {
        return(NULL)
    }
    set_names(
        map(
            tokens,
            .twasFmColumnsFor,
            fmRes = fmRes,
            cond = cond,
            joint = joint
        ),
        tokens
    )
}

# @noRd
.twasFmColumnsFor <- function(token, fmRes, cond, joint) {
    map_int(
        seq_len(nrow(cond)),
        .twasFmColumnAt,
        token = token,
        fmRes = fmRes,
        cond = cond,
        joint = joint
    )
}

# @noRd
.twasFmColumnAt <- function(r, token, fmRes, cond, joint) {
    rows <- .twasFmTokenRows(
        fmRes,
        token,
        as.character(cond$study[[r]]),
        as.character(cond$context[[r]]),
        as.character(cond$trait[[r]]),
        joint
    )
    if (length(rows) == 0L) {
        return(NA_integer_)
    }
    as.integer(.fmrJointIndexAt(fmRes, rows[[1L]]))
}

# The group's fine-mapping CV with one prediction column and one performance
# row per condition, in the group's order -- the shape .jointTwasCvResult()
# gives, which the per-condition slice indexes. A joint fit stores each
# condition's slice on that condition's row, so the first row alone covered
# only the first condition. The rows are those of the group's own fit
# (`joint`). A method is handed off only when every condition carries its
# CV, and only from one fold partition.
# @noRd
.twasGroupFmCv <- function(fmRes, cond, joint) {
    cvs <- pmap(
        list(
            studyName = as.character(cond$study),
            context = as.character(cond$context),
            trait = as.character(cond$trait)
        ),
        .twasCvResultFor,
        fineMappingResult = fmRes,
        joint = joint
    )
    if (length(cvs) == 1L) {
        return(cvs[[1L]])
    }
    if (any(map_lgl(cvs, is.null))) {
        return(NULL)
    }
    .twasCheckGroupPartitions(map(cvs, "samplePartition"), cond)
    list(
        samplePartition = cvs[[1L]]$samplePartition,
        prediction = .twasCvBindKeys(map(cvs, "prediction"), cbind),
        performance = .twasCvBindKeys(map(cvs, "performance"), rbind)
    )
}

# The conditions of one joint group were cross-validated together, so their
# fine-mapping rows share one fold partition. Rows from separate runs would
# pair out-of-fold predictions made on different folds.
# @noRd
.twasCheckGroupPartitions <- function(partitions, cond) {
    if (all(map_lgl(partitions, identical, partitions[[1L]]))) {
        return(invisible(NULL))
    }
    labels <- str_c(cond$study, cond$context, cond$trait, sep = ":")
    abort(glue(
        "twasWeightsPipeline: the fine-mapping cross-validation of ",
        "{str_c(labels, collapse = ', ')} used different fold partitions; ",
        "a joint group needs the CV of one joint fit."
    ))
}

# Per key every condition carries, the conditions' pieces bound into one
# matrix: cbind for sample x condition predictions, rbind for condition x
# metric performance.
# @noRd
.twasCvBindKeys <- function(pieces, bind) {
    keys <- reduce(map(pieces, names), intersect)
    set_names(map(keys, .twasCvBindKey, pieces = pieces, bind = bind), keys)
}

# @noRd
.twasCvBindKey <- function(key, pieces, bind) {
    keyValues <- map(pieces, key)
    exec(bind, !!!keyValues)
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
    n <- min(length(entries), nrow(cond))
    recs <- map(
        seq_len(n),
        .jointEntryRecordAt,
        entries = entries,
        cond = cond,
        method = method,
        grp = grp,
        data = data,
        jointIndex = .jointIndexOf(entries[seq_len(n)], grp)
    )
    compact(recs)
}

# Each entry's position in its joint fit's condition axis: its rank among the
# group's fitted conditions, as a condition the SER screen dropped has no
# entry. NA for a group that is not a joint fit.
# @noRd
.jointIndexOf <- function(entries, grp) {
    out <- rep(NA_integer_, length(entries))
    if (!.jointGroupIsJoint(grp)) {
        return(out)
    }
    fitted <- !map_lgl(entries, is.null)
    replace(out, fitted, seq_len(sum(fitted)))
}

# Whether group `grp` joins conditions (a joint* axis is set).
# @noRd
.jointGroupIsJoint <- function(grp) {
    any(map_lgl(list(grp$js, grp$jc, grp$jt), .jointAxisIsSet))
}

# @noRd
.jointAxisIsSet <- function(value) {
    length(value) == 1L && !is.na(value)
}

# Run one dispatch cell: enumerate joint groups, fit each method (S4 dispatch on
# the group x pipeline pair) per group, accumulate per-context rows, build the
# per-pipeline result. The loop is GROUP-outer / token-inner so the twas
# ensemble layer can combine a group's per-method fits in place (FM is
# unaffected by the loop order). Per-method fitting is identical for FM and twas
# -- one method -> per-condition entries; the SR-TWAS ensemble is a layer ON TOP
# of that.
.runJointCell <- function(cell, pipeline, data, scope, tokens, args = list()) {
    .qtlRunWork(
        .jointCellWork(cell, pipeline, data, scope, tokens, args),
        data,
        args$genotypeFilter,
        args$residualizationParam,
        is(pipeline, "FmJointPipeline")
    )
}

# One joint cell's fits, planned: its units enumerated, none of them built yet.
# @noRd
.jointCellWork <- function(cell, pipeline, data, scope, tokens, args = list()) {
    # PC groups are labelled by component, not by trait: their enumerator
    # applies the selection itself, to the traits the components come from.
    selection <- if (identical(.jcPattern(cell), "topPc")) {
        NULL
    } else {
        scope$traitSelection
    }
    doEnsemble <- is(pipeline, "TwasJointPipeline") &&
        isTRUE(.jpConfig(pipeline)$ensembleParam$enabled)
    units <- .jcEnumerate(cell)(data, scope, args)
    .qtlWork(
        units,
        map(units, .unitSpec, data = data, args = args),
        .jointCellUnit,
        list(
            selection = selection,
            minGroup = .jcMinGroup(cell),
            pipeline = pipeline,
            data = data,
            tokens = tokens,
            args = args,
            doEnsemble = doEnsemble
        ),
        .jointCellResult,
        list(pipeline = pipeline)
    )
}

# A joint cell's result, from its units' fits.
# @noRd
.jointCellResult <- function(fits, pipeline) {
    fitted <- discard(fits, is.null)
    if (length(fitted) == 0L) {
        return(NULL)
    }
    construct(pipeline, .jeConcat(fitted))
}

# One enumerated unit built, filtered and fitted: the records of its groups,
# or NULL when none survives the selection and the minimum size. The unit is
# built here, inside the loop, so a residualized X of its own is let go before
# the next unit's is built.
# @noRd
.jointCellUnit <- function(
    unit,
    selection,
    minGroup,
    pipeline,
    data,
    tokens,
    args,
    doEnsemble,
    block = NULL
) {
    groups <- .unitGroups(unit, block) |>
        map(.jointGroupSelect, traitId = selection) |>
        compact() |>
        keep(.jointGroupMeetsMin, minGroup = minGroup)
    if (length(groups) == 0L) {
        return(NULL)
    }
    list_flatten(map(
        groups,
        .jointCellGroup,
        pipeline = pipeline,
        data = data,
        tokens = tokens,
        args = args,
        doEnsemble = doEnsemble,
        block = block
    ))
}

# A unit of joint groups not built yet: the function that builds it from a
# residualized X, its arguments, `x` first, and the genotypes it will be built
# on (the traits that name their window, the contexts of their covariates and
# the samples it may fit). The individual-level enumerators return these: their
# phenotypes are in hand, so what can share an X is known before any genotype is
# read, and a unit's X exists only while its groups are fitted. A unit is one
# group, or a context's principal components, which share one X.
# @noRd
.deferGroups <- function(x, builder, xTraits, xContexts, samples, ...) {
    list(
        builder = builder,
        args = list(x, ...),
        xTraits = xTraits,
        xContexts = xContexts,
        samples = samples
    )
}

# A unit not built yet that reads its own genotypes when it is: a cross-study
# group needs a block from each study, which no one dataset's read gives, so the
# builder reads them itself and is handed no X.
# @noRd
.deferSelfReading <- function(x, builder, ...) {
    list(builder = builder, args = list(x, ...), readsOwn = TRUE)
}

# The genotypes a unit's fits are on; none for a group already built, which
# holds its own, or for a unit that reads its own.
# @noRd
.unitSpec <- function(unit, data, args) {
    if (is(unit, "JointGroup") || isTRUE(unit$readsOwn)) {
        return(NULL)
    }
    .qtlFitSpec(
        .fmGenoWindow(data, unit$xTraits, args$region, args$cisWindow),
        unit$xContexts,
        unit$samples
    )
}

# The groups of an enumerated unit: built on `block` when it was deferred, as
# given when it was enumerated built (the summary-statistics groups, which hold
# no X).
# @noRd
.unitGroups <- function(unit, block = NULL) {
    built <- if (is(unit, "JointGroup")) unit else .unitBuild(unit, block)
    if (is(built, "JointGroup")) {
        return(list(built))
    }
    compact(as.list(built))
}

# A deferred unit built: on `block`, or on the blocks it reads itself.
# @noRd
.unitBuild <- function(unit, block) {
    unitArgs <- unit$args
    if (isTRUE(unit$readsOwn)) {
        return(exec(unit$builder, !!!unitArgs))
    }
    exec(unit$builder, !!!unitArgs, X = block$X)
}

# The group less any (context, trait) condition a per-context `traitId`
# selection leaves out, or NULL when it leaves out every one. A NULL
# selection keeps the group whole.
# @noRd
.jointGroupSelect <- function(g, traitId) {
    if (is.null(traitId)) {
        return(g)
    }
    cond <- .jgConditions(g)
    .jgKeepConditions(
        g,
        map_lgl(
            seq_len(nrow(cond)),
            .jointConditionSelectedAt,
            cond = cond,
            traitId = traitId
        )
    )
}

# Whether the selection keeps condition `i` of a group's conditions table.
# @noRd
.jointConditionSelectedAt <- function(i, cond, traitId) {
    .traitSelectedAt(
        cond$trait[[i]],
        traitId,
        cond$study[[i]],
        cond$context[[i]]
    )
}

# Fit every token for one joint group -> records, plus the optional SR-TWAS
# ensemble layer (built ON TOP of the shared per-method fits, when >= 2 methods
# produced entries).
# @noRd
.jointCellGroup <- function(
    g,
    pipeline,
    data,
    tokens,
    args,
    doEnsemble,
    block = NULL
) {
    grp <- .jointGroupProvenance(g, data, args)
    rowTokens <- .jointRowTokens(tokens, pipeline, grp$cond)
    fitArgs <- .jointGroupFitArgs(g, pipeline, args, block, rowTokens)
    perTokenEntries <- .jointRowEntries(
        rowTokens,
        g,
        pipeline,
        fitArgs,
        grp,
        args
    )
    tokenRecords <- .jeConcat(map(
        names(perTokenEntries),
        .jointTokenRecords,
        perTokenEntries = perTokenEntries,
        grp = grp,
        data = data
    ))
    stacking <- .twasStacksGroup(pipeline, grp$cond)
    if (!doEnsemble || (length(perTokenEntries) < 2L && !stacking)) {
        return(tokenRecords)
    }
    c(
        tokenRecords,
        .jointEntryRecords(
            .jointEnsembleEntries(
                g,
                perTokenEntries,
                pipeline,
                list_assign(args, samplePartition = fitArgs$samplePartition),
                block,
                stacking
            ),
            "ensemble",
            grp,
            data
        )
    )
}

# Each of `tokens`' entries on group `g`, by token, those producing none left
# out.
# @noRd
.jointRowEntries <- function(tokens, g, pipeline, fitArgs, grp, args) {
    compact(set_names(
        map(
            tokens,
            .jointTokenEntriesOrNull,
            g = g,
            pipeline = pipeline,
            fitArgs = fitArgs,
            cond = grp$cond,
            args = args
        ),
        tokens
    ))
}

# The tokens whose fits of a group become its rows. A run stacking across
# studies (EnsembleParam(acrossStudies = TRUE)) also hands a group across
# studies the methods it only stacks, fitted there one study at a time
# (.twasStackEntries()); those make no rows of their own, as each study's
# own run has them. The methods it fits across the studies, coupled or
# jointly (`multitraitTokens`), do.
# @noRd
.jointRowTokens <- function(tokens, pipeline, cond) {
    if (!.twasStacksGroup(pipeline, cond)) {
        return(tokens)
    }
    cfg <- .jpConfig(pipeline)
    coupled <- if (is.null(cfg$studyCoupling)) {
        character(0)
    } else {
        union(names(.TWAS_COUPLED_IMPL), .TWAS_COUPLED_SUMSTATS)
    }
    intersect(tokens, union(coupled, cfg$multitraitTokens))
}

# Whether a run stacks across the studies of a group: its ensemble is on and
# across studies, and the group spans two or more.
# @noRd
.twasStacksGroup <- function(pipeline, cond) {
    ens <- if (is(pipeline, "TwasJointPipeline")) {
        .jpConfig(pipeline)$ensembleParam
    }
    isTRUE(ens$enabled) &&
        isTRUE(ens$acrossStudies) &&
        n_distinct(cond$study) >= 2L
}

# A group's ensemble entries, per condition: its methods' out-of-fold
# predictions stacked (.twasEnsembleLayer()), with each study's own fits
# (.twasStackEntries()) when it is stacked across studies, cross-validated on
# the fold partition its methods share (`args$samplePartition`). On summary
# statistics only a group across studies is stacked here, on held-out splits
# (.twasSumstatEnsembleLayer()), as each study's own methods are stacked one
# tuple at a time.
# @noRd
.jointEnsembleEntries <- function(
    g,
    perTokenEntries,
    pipeline,
    args,
    block,
    stacking
) {
    cfg <- .jpConfig(pipeline)
    if (is(g, "SumStatsJointGroup")) {
        if (!stacking) {
            return(list())
        }
        return(.twasSumstatEnsembleLayer(g, perTokenEntries, cfg, args))
    }
    stack <- if (stacking) .twasStackEntries(g, pipeline, args, block)
    .twasEnsembleLayer(g, perTokenEntries, cfg, stack %||% list())
}

# Each study's own fit of every method a group across studies stacks, as
# entries by condition: fitted on the study's samples alone (the other
# studies' rows are missing from its condition), neither coupled nor joint,
# and never read from the resume cache, which keeps the group's fits across
# studies under the same tokens.
# @noRd
.twasStackEntries <- function(g, pipeline, args, block) {
    tokens <- .jpConfig(pipeline)$ensembleParam$stackTokens
    uncoupled <- new(
        "TwasJointPipeline",
        config = list_assign(
            .jpConfig(pipeline),
            studyCoupling = NULL,
            multitraitTokens = NULL
        )
    )
    fitArgs <- .jointGroupFitArgs(g, uncoupled, args, block, tokens)
    compact(set_names(
        map(
            tokens,
            .twasStackFit,
            g = g,
            pipeline = uncoupled,
            fitArgs = fitArgs
        ),
        tokens
    ))
}

# @noRd
#' @importFrom rlang try_fetch
.twasStackFit <- function(token, g, pipeline, fitArgs) {
    entries <- try_fetch(
        fitJointGroup(g, pipeline, token, fitArgs),
        error = function(cnd) {
            warn(
                glue(
                    "twasWeightsPipeline: {token} failed on a study of ",
                    "the ensemble across studies."
                ),
                parent = cnd
            )
            NULL
        }
    )
    if (length(entries) == 0L) NULL else entries
}

# A group's provenance: the ";"-joined members of each varying axis (identical
# on every per-context row of this group), and the window this group was
# analysed over -- the caller's region when one was named, else the group's own
# cis window (the union span of its traits +/- cisWindow). The window is
# recorded because the fit only ever sees whichever variants survived QC; it is
# read back by genomicRegion().
# @noRd
.jointGroupProvenance <- function(g, data, args) {
    cond <- .jgConditions(g)
    key <- .jointKeyOfConditions(cond)
    list(
        cond = cond,
        js = key$jointStudies,
        jc = key$jointContexts,
        jt = key$jointTraits,
        region = .jointGroupWindow(data, cond, args),
        ldReference = .jgLdRefOf(g)
    )
}

# The arguments a joint group's fitters take, with what they cannot recover for
# themselves: allele frequencies live on the block, not inside fitJointGroup(),
# and the requested window is only ever known to the caller. Both are read back
# as args$af / args$region where the entries are post-processed.
# @noRd
.jointGroupFitArgs <- function(g, pipeline, args, block, tokens) {
    af <- .jointAfForGroup(g, block)
    list_assign(
        .twasGroupArgs(g, pipeline, args, tokens),
        af = af,
        region = args$region,
        h2 = .jointGroupH2(g, pipeline, tokens)
    )
}

# What a TWAS run's `h2` gives an individual-level group's one trait
# (.twasRowH2()), made once, before the group's methods are fitted, on the
# genotypes and phenotype they fit.
# @noRd
.jointGroupH2 <- function(g, pipeline, tokens) {
    if (!is(pipeline, "TwasJointPipeline") || !is(g, "IndividualJointGroup")) {
        return(NULL)
    }
    cond <- .jgConditions(g)
    spec <- .jpConfig(pipeline)$h2
    if (is.null(spec) || nrow(cond) != 1L) {
        return(NULL)
    }
    .twasRowH2(
        spec,
        .h2TupleKey(cond$study, cond$context, cond$trait),
        glue(
            "twasWeightsPipeline(QtlDataset): study='{cond$study}', ",
            "context='{cond$context}', trait='{cond$trait}'"
        ),
        .h2IndividualFit,
        list(X = .jgX(g), y = .jgY(g)[, 1L]),
        tokens
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
# unless every condition is cached from a fit of this same group.
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
        token = token,
        joint = .jointKeyOfConditions(cond)
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
.twasEnsembleLayer <- function(group, perTokenEntries, cfg, stack = list()) {
    tokens <- names(perTokenEntries)
    Y <- .jgY(group)
    r2Cut <- if (is.null(cfg$ensembleParam$r2Threshold)) {
        0.01
    } else {
        cfg$ensembleParam$r2Threshold
    }
    solver <- if (is.null(cfg$ensembleParam$solver)) {
        "quadprog"
    } else {
        cfg$ensembleParam$solver
    }
    alpha <- if (is.null(cfg$ensembleParam$alpha)) {
        1
    } else {
        cfg$ensembleParam$alpha
    }
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
        cfg = cfg,
        stack = stack,
        group = group
    )
}

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
        },
        partition = cv$samplePartition
    )
}

# Per-condition CV predictions + weights + R^2 across the group's methods, and
# the fold partition their predictions share. Returns list(preds, wts, rsq,
# partition), skipping methods with no CV / weights.
# @noRd
.twasEnsembleCollect <- function(r, perTokenEntries, tokens) {
    .twasEnsembleParts(compact(map(
        tokens,
        .twasEnsembleTokenPart,
        r = r,
        perTokenEntries = perTokenEntries
    )))
}

# The ensemble inputs of the predictors `parts`, each named by its `token`,
# and the fold partition the first one with out-of-fold predictions shares.
# @noRd
.twasEnsembleParts <- function(parts) {
    if (length(parts) == 0L) {
        return(list(preds = list(), wts = list(), rsq = c()))
    }
    contributing <- map_chr(parts, "token")
    partitions <- compact(map(parts, "partition"))
    list(
        preds = set_names(
            map(parts, "pred"),
            str_c(contributing, "_predicted")
        ),
        wts = set_names(
            map(parts, "weights"),
            str_c(contributing, "_weights")
        ),
        rsq = set_names(map_dbl(parts, "rsq"), contributing),
        partition = if (length(partitions) > 0L) partitions[[1L]]
    )
}

# Two predictor sets' ensemble inputs as one.
# @noRd
.twasEnsembleMerge <- function(a, b) {
    list(
        preds = c(a$preds, b$preds),
        wts = c(a$wts, b$wts),
        rsq = c(a$rsq, b$rsq),
        partition = a$partition %||% b$partition
    )
}

# Condition `r`'s predictors from each study's own fits `stack`, named
# "study/method": its own study's by their out-of-fold predictions, every
# other study's as fixed predictors on its samples
# (.twasStackExternalPart()). None without a stack.
# @noRd
.twasStackCollect <- function(r, stack, group, scale) {
    if (length(stack) == 0L) {
        return(list(preds = list(), wts = list(), rsq = c()))
    }
    studies <- .jgConditions(group)$study
    own <- map2(
        map(
            names(stack),
            .twasEnsembleTokenPart,
            r = r,
            perTokenEntries = stack
        ),
        str_c(studies[[r]], "/", names(stack)),
        .twasRenamePart
    )
    external <- list_flatten(map(
        setdiff(seq_along(studies), r),
        .twasStackExternalParts,
        r = r,
        stack = stack,
        group = group,
        scale = scale
    ))
    .twasEnsembleParts(compact(c(own, external)))
}

# @noRd
.twasRenamePart <- function(part, name) {
    if (is.null(part)) NULL else list_assign(part, token = name)
}

# Study `k`'s own fits as fixed predictors of condition `r`.
# @noRd
.twasStackExternalParts <- function(k, r, stack, group, scale) {
    map(
        names(stack),
        .twasStackExternalPart,
        k = k,
        r = r,
        stack = stack,
        group = group,
        scale = scale
    )
}

# Study k's fit of `tk`, carried to study r's genotypes
# (.twasTransferWeights()) over the variants that vary among r's samples,
# and applied to the group's samples, as its own methods' out-of-fold
# predictions cover them all (the other studies' outcomes are missing from
# r's condition). Its R^2 on r's samples is honest, as none of them fitted
# it. NULL when k has no weights for `tk`.
# @noRd
.twasStackExternalPart <- function(tk, k, r, stack, group, scale) {
    e <- stack[[tk]][[k]]
    w <- if (!is.null(e)) .rowWeights(e)
    if (is.null(w)) {
        return(NULL)
    }
    X <- .jgX(group)
    observed <- !is.na(.jgY(group)[, r])
    carried <- .twasTransferWeights(
        set_names(as.numeric(w), .rowVariantIds(e))[colnames(X)],
        .jgAf(group),
        c(k, r),
        scale
    )
    varying <- .nonzeroVarColumns(X[observed, , drop = FALSE])
    used <- replace(carried, !varying | is.na(carried), 0)
    pred <- drop(X %*% used)
    list(
        token = str_c(.jgConditions(group)$study[[k]], "/", tk),
        pred = matrix(pred, ncol = 1L, dimnames = list(rownames(X), NULL)),
        weights = matrix(used, ncol = 1L, dimnames = list(colnames(X), NULL)),
        rsq = .cvMetricRow(pred[observed], .jgY(group)[observed, r])[["rsq"]],
        partition = NULL
    )
}

# Weights `w` of study `pair[1]` carried to study `pair[2]`'s genotypes: as
# they are per SD ("perSd"), or per allele ("perAllele"), each scaled by the
# ratio of the target's genotype SD to the source's, sqrt(2 p (1 - p)) from
# the studies' allele frequencies `af` (variants x studies), and 0 where the
# source's is 0 or unknown.
# @noRd
.twasTransferWeights <- function(w, af, pair, scale) {
    if (!identical(scale, "perAllele")) {
        return(w)
    }
    sdFrom <- sqrt(
        2 * af[names(w), pair[[1L]]] * (1 - af[names(w), pair[[1L]]])
    )
    sdTo <- sqrt(2 * af[names(w), pair[[2L]]] * (1 - af[names(w), pair[[2L]]]))
    w * if_else(is.finite(sdFrom) & sdFrom > 0, sdTo / sdFrom, 0)
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
    cfg,
    stack = list(),
    group = NULL
) {
    coll <- .twasEnsembleMerge(
        .twasEnsembleCollect(r, perTokenEntries, tokens),
        .twasStackCollect(r, stack, group, cfg$ensembleParam$transferScale)
    )
    passing <- names(coll$rsq)[!is.na(coll$rsq) & coll$rsq >= r2Cut]
    if (length(passing) < 2L) {
        return(NULL)
    }
    ens <- try_fetch(
        ensembleWeights(
            cvResults = list(
                prediction = coll$preds[str_c(passing, "_predicted")],
                samplePartition = coll$partition
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
    .twasEnsembleRow(
        ens,
        coll$partition,
        rownames(coll$wts[[str_c(passing[1L], "_weights")]]),
        stdz,
        if (length(stack) > 0L) .twasConditionVarying(group, r)
    )
}

# Which of a group's variants vary among condition `r`'s samples, by name: in
# a group across studies, a variant r's study lacks is zero there, and the
# ensemble gives it no weight.
# @noRd
.twasConditionVarying <- function(group, r) {
    X <- .jgX(group)[!is.na(.jgY(group)[, r]), , drop = FALSE]
    set_names(.nonzeroVarColumns(X), colnames(X))
}

# The ensemble's row: its combined weights, and a cvResult recording its
# cross-fit out-of-fold predictions and metrics over its methods' folds, as
# each method's does, beside the stacking coefficients and the methods' own
# R^2. `fallbackIds` names the weights when the ensemble's do not.
# @noRd
.twasEnsembleRow <- function(
    ens,
    partition,
    fallbackIds,
    stdz,
    varying = NULL
) {
    ew <- ens$ensembleTwasWeights
    vids <- (if (!is.null(names(ew))) names(ew) else rownames(ew)) %||%
        fallbackIds
    kept <- if (is.null(varying)) 1 else as.numeric(varying[vids])
    TwasWeightsRow(
        variantIds = vids,
        weights = as.numeric(ew) * kept,
        cvResult = compact(list(
            samplePartition = partition,
            predictions = ens$predictions,
            metrics = ens$metrics,
            methodCoef = ens$methodCoef,
            methodPerformance = ens$methodPerformance
        )),
        weightStandardized = stdz
    )
}

# ---- The ensemble across summary-statistics studies -------------------------
#
# A summary-statistics group across studies stacked as LEOPARD combines
# ancestries (EnsembleParam(acrossStudies = TRUE)): per condition, the
# study's own predictors -- the group's fits across studies and each
# method's own fit of the study -- refit on held-out splits of its
# statistics, beside every other study's own fits as fixed predictors
# carried to the study, all through .sumstatEnsemble().

# Each condition's ensemble entry, NULL where fewer than two predictors clear
# the threshold.
# @noRd
.twasSumstatEnsembleLayer <- function(g, perTokenEntries, cfg, args) {
    studyFits <- map(
        seq_len(ncol(.jgZ(g))),
        .twasSumstatStudyFits,
        g = g,
        cfg = cfg,
        args = args
    )
    map(
        seq_along(studyFits),
        .twasSumstatStackCondition,
        g = g,
        perTokenEntries = perTokenEntries,
        studyFits = studyFits,
        cfg = cfg,
        args = args
    )
}

# Condition `k`'s own fits of every method the run stacks, made as the
# per-tuple path makes them (.twasQssUnivariateFitOne()) on the study's
# statistics and LD reference, with the fit context their refits read.
# @noRd
.twasSumstatStudyFits <- function(k, g, cfg, args) {
    cond <- .jgConditions(g)
    key <- c(cond$study[[k]], cond$context[[k]], cond$trait[[k]])
    held <- .jointCoupledStat(k, .jgZ(g), .jgN(g))
    rowCtx <- list(
        variantIds = held$variantNames,
        stat = list(
            z = held$z,
            n = .jgN(g)[[k]],
            variantNames = held$variantNames
        ),
        ldReference = .ldRefInMemory(
            .jgLdRefAt(g, k),
            held$variantNames,
            label = "twasWeightsPipeline"
        )
    )
    tokens <- cfg$ensembleParam$stackTokens
    fitCtx <- c(
        rowCtx,
        list(
            h2 = .twasQssRowH2(
                cfg$h2,
                key[[1L]],
                key[[2L]],
                key[[3L]],
                rowCtx,
                tokens
            )
        )
    )
    list(
        fitCtx = fitCtx,
        rows = compact(map(
            tokens,
            .twasQssUnivariateFitOne,
            st = key[[1L]],
            ctx = key[[2L]],
            tr = key[[3L]],
            fitCtx = fitCtx,
            methodArgs = args$methodArgs,
            fitRetention = "none",
            fineMappingResult = NULL
        ))
    )
}

# Condition `r`'s ensemble entry (.sumstatEnsemble()), or NULL.
# @noRd
.twasSumstatStackCondition <- function(
    r,
    g,
    perTokenEntries,
    studyFits,
    cfg,
    args
) {
    fitCtx <- studyFits[[r]]$fitCtx
    own <- .twasSumstatOwnMethods(r, g, perTokenEntries, studyFits, cfg, args)
    external <- .twasSumstatExternal(
        r,
        g,
        studyFits,
        cfg$ensembleParam$transferScale
    )
    if (length(own) + ncol(external) < 2L) {
        return(NULL)
    }
    stat <- fitCtx$stat
    ld <- .rssLdSource(stat, fitCtx$ldReference, "twasWeightsPipeline")
    ens <- .sumstatEnsemble(
        own,
        external,
        .lassosumClampCor(
            .rssMarginalCor(stat, n = stat$n, p = .ldSourceSize(ld))
        ),
        stat$n,
        ld,
        cfg$ensembleParam$r2Threshold
    )
    if (is.null(ens)) {
        return(NULL)
    }
    .twasSumstatEnsembleRow(ens, fitCtx$variantIds)
}

# Condition `r`'s own predictors: the group's fits across studies, by token,
# and each method's own fit of its study, as "study/method".
# @noRd
.twasSumstatOwnMethods <- function(
    r,
    g,
    perTokenEntries,
    studyFits,
    cfg,
    args
) {
    fits <- studyFits[[r]]
    c(
        compact(map(
            set_names(names(perTokenEntries)),
            .twasCoupledStackMethod,
            r = r,
            g = g,
            perTokenEntries = perTokenEntries,
            cfg = cfg,
            args = args
        )),
        compact(set_names(
            map(
                fits$rows,
                .twasQssStackMethod,
                fitCtx = fits$fitCtx,
                methodArgs = args$methodArgs
            ),
            str_c(
                .jgConditions(g)$study[[r]],
                "/",
                map_chr(fits$rows, "method")
            )
        ))
    )
}

# The group's fit across studies of `token` as condition `r`'s predictor:
# its weights over r's variants, and the refit that replaces r's statistics
# by a held-out split's training part and fits across the studies again
# (.twasCoupledRefit()). NULL when the fit has no weights for r.
# @noRd
.twasCoupledStackMethod <- function(token, r, g, perTokenEntries, cfg, args) {
    e <- perTokenEntries[[token]][[r]]
    w <- if (!is.null(e)) .rowWeights(e)
    if (is.null(w)) {
        return(NULL)
    }
    Z <- .jgZ(g)
    stats <- map(seq_len(ncol(Z)), .jointCoupledStat, Z = Z, N = .jgN(g))
    ids <- stats[[r]]$variantNames
    list(
        weights = replace_na(
            unname(set_names(as.numeric(w), .rowVariantIds(e))[ids]),
            0
        ),
        refit = list(
            fn = ".twasCoupledRefit",
            args = list(
                token = token,
                stats = stats,
                r = r,
                g = g,
                cfg = cfg,
                args = args
            )
        )
    )
}

# The fit across studies of `token` on a group whose condition `r` is
# replaced by a held-out split's training part (`trainCor` on `share` of its
# samples): r's weights over its own variants.
# @noRd
.twasCoupledRefit <- function(trainCor, share, token, stats, r, g, cfg, args) {
    n <- share * stats[[r]]$n
    trained <- list_assign(
        stats[[r]],
        z = trainCor * sqrt((n - 2) / (1 - trainCor^2)),
        n = n
    )
    W <- .jointCoupledFit(
        replace(stats, r, list(trained)),
        g,
        cfg,
        args,
        token,
        "none"
    )
    unname(W[stats[[r]]$variantNames, r])
}

# The other studies' own fits as fixed predictors of condition `r`, over its
# variants (variants x predictors, named "study/method"): each carried to r
# (.twasSumstatTransfer()), 0 where its study lacks a variant.
# @noRd
.twasSumstatExternal <- function(r, g, studyFits, scale) {
    ids <- studyFits[[r]]$fitCtx$variantIds
    columns <- list_flatten(map(
        setdiff(seq_along(studyFits), r),
        .twasSumstatExternalStudy,
        r = r,
        g = g,
        studyFits = studyFits,
        ids = ids,
        scale = scale
    ))
    if (length(columns) == 0L) {
        return(matrix(0, length(ids), 0L))
    }
    exec(cbind, !!!columns)
}

# Study `k`'s own fits over condition r's variants `ids`, one named column
# each.
# @noRd
.twasSumstatExternalStudy <- function(k, r, g, studyFits, ids, scale) {
    rows <- studyFits[[k]]$rows
    set_names(
        map(
            rows,
            .twasSumstatExternalColumn,
            ids = ids,
            refs = list(.jgLdRefAt(g, k), .jgLdRefAt(g, r)),
            scale = scale
        ),
        str_c(.jgConditions(g)$study[[k]], "/", map_chr(rows, "method"))
    )
}

# @noRd
.twasSumstatExternalColumn <- function(row, ids, refs, scale) {
    w <- replace_na(
        unname(set_names(weights(row$entry), variantIds(row$entry))[ids]),
        0
    )
    .twasSumstatTransfer(w, ids, refs, scale)
}

# Weights `w` over variants `ids` of the study on LD reference `refs[[1]]`
# carried to the study on `refs[[2]]`: as they are per SD ("perSd"), or per
# allele ("perAllele"), each scaled by the ratio of the target's genotype SD
# to the source's, from each reference's genotype variances
# (.ldRefVariance()), which needs every source variant's allele frequency.
# @noRd
.twasSumstatTransfer <- function(w, ids, refs, scale) {
    if (!identical(scale, "perAllele") || !any(w != 0)) {
        return(w)
    }
    used <- w != 0
    variances <- map(refs, .twasRefVariance, ids = ids[used])
    if (anyNA(variances[[1L]]) || anyNA(variances[[2L]])) {
        abort(c(
            str_c(
                "twasWeightsPipeline: EnsembleParam(transferScale = ",
                "\"perAllele\") needs each study's allele frequencies, and ",
                "an LD reference records none."
            ),
            i = "Give the LdMatrix's variants an `AF` column, or use \"perSd\"."
        ))
    }
    replace(
        w,
        used,
        w[used] *
            if_else(
                variances[[1L]] > 0,
                sqrt(variances[[2L]] / variances[[1L]]),
                0
            )
    )
}

# Each variant's genotype variance on an LD reference (.ldRefVariance()), in
# the order of `ids`.
# @noRd
.twasRefVariance <- function(ref, ids) {
    .ldRefVariance(
        ref,
        .ldFromRefMatch(ref, ids, "twasWeightsPipeline", "error")$idx
    )
}

# Top-PC enumeration: PCA-reduce each context's multi-trait phenotype and make
# ONE group per top principal component, each a single-column Y named topPCk.
# Structurally identical to .enumUnivariateIndividual -- same residualization
# helpers, same IndividualJointGroup -- so everything downstream (CV,
# ensemble, retention) applies unchanged. It is the TWAS peer of
# fineMappingPipeline's .fmPcaContextPlan, sharing .fmTopPcScores().
# @noRd
.enumTopPcIndividual <- function(data, scope, args = list()) {
    studyName <- studyName(data)
    if (!is_in(studyName, scope$studies)) {
        return(list())
    }
    phenotypeFilter <- .jointPhenotypeFilter(args)
    # One unit per context: its components share the context's one X.
    compact(map(
        scope$contexts[[studyName]],
        .deferTopPc,
        data = data,
        scope = scope,
        args = args,
        studyName = studyName,
        phenotypeFilter = phenotypeFilter
    ))
}

# One context's deferred top-PC unit: its window spans the traits it selects.
# NULL when PCA has nothing to work on: it is undefined for a single trait,
# matching the fine-mapping rule.
# @noRd
.deferTopPc <- function(cx, data, scope, args, studyName, phenotypeFilter) {
    traits <- .enumContextTraits(cx, data, scope, studyName)
    if (length(traits) < 2L) {
        return(NULL)
    }
    Y <- .fmResidPheno(
        data,
        contexts = cx,
        traitId = traits,
        phenotypeFilter = phenotypeFilter,
        residualizationParam = args$residualizationParam
    )
    scores <- .fmTopPcScores(Y, args$nPCs %||% 10L)
    if (is.null(scores)) {
        return(NULL)
    }
    .deferGroups(
        cx,
        .enumTopPcForContext,
        xTraits = traits,
        xContexts = cx,
        samples = rownames(scores),
        studyName = studyName,
        scores = scores
    )
}

# One context's top-PC groups on `X`.
# @noRd
.enumTopPcForContext <- function(cx, studyName, scores, X) {
    common <- intersect(rownames(X), rownames(scores))
    if (length(common) < 2L) {
        return(list())
    }
    compact(map(
        colnames(scores),
        .enumTopPcGroup,
        X = X[common, , drop = FALSE],
        scores = scores[common, , drop = FALSE],
        studyName = studyName,
        cx = cx
    ))
}

# @noRd
.enumTopPcGroup <- function(pcName, X, scores, studyName, cx) {
    new(
        "IndividualJointGroup",
        conditions = tibble(study = studyName, context = cx, trait = pcName),
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
    # A MultiStudyQtlDataset's individual-level studies, block-missing Y.
    new(
        "JointDispatchCell",
        pattern = "study",
        dataForm = "individual",
        enumerate = .enumCrossStudyIndividual,
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
# per-run engine payload (mrmashPrior, methodArgs, cisWindow, region, ...).
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
    .qtlRunWork(
        .jointSpecsWork(
            parsedJointSpec,
            data,
            dataForm,
            pipeline,
            jointMethods,
            contexts,
            traitIds,
            args
        ),
        data,
        args$genotypeFilter,
        args$residualizationParam,
        is(pipeline, "FmJointPipeline")
    )
}

# .runJointSpecs()' fits, planned, for a caller with others to run alongside.
# @noRd
.jointSpecsWork <- function(
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
        return(.qtlWorkValue(NULL))
    }
    .qtlWorkOf(
        map(
            parsedJointSpec,
            .jointSpecWork,
            data = data,
            dataForm = dataForm,
            pipeline = pipeline,
            jointMethods = jointMethods,
            contexts = contexts,
            traitIds = traitIds,
            args = args
        ),
        .jointSpecsResult,
        list(pipeline = pipeline)
    )
}

# The specs' results rbound. Every row keeps the LD reference its group was
# fitted against.
# @noRd
.jointSpecsResult <- function(perSpec, pipeline) {
    results <- compact(perSpec)
    if (length(results) == 0L) {
        return(NULL)
    }
    isFm <- is(pipeline, "FmJointPipeline")
    reduce(results, if (isFm) .rbindFineMappingResult else .rbindTwasWeights)
}

# Plan all joint methods for ONE spec: resolve its scope (optionally
# region-restricting the traits), pick the joint cell, and plan it
# (.jointCellWork(), with ALL methods so the twas ensemble layer can combine a
# group's per-method fits).
# @noRd
.jointSpecWork <- function(
    spec,
    data,
    dataForm,
    pipeline,
    jointMethods,
    contexts,
    traitIds,
    args
) {
    scope <- .jointSpecScope(
        spec,
        data,
        dataForm,
        contexts,
        traitIds,
        args$region
    )
    pattern <- if (length(spec$axes) > 1L) "composed" else spec$axes[[1L]]
    cell <- .lookupJointCell(pattern, dataForm)
    spArgs <- c(args, list(axes = spec$axes))
    .jointCellWork(cell, pipeline, data, scope, jointMethods, spArgs)
}

# One spec's scope: its own, narrowed to the run's selection and, on
# individual-level data with no trait named, to the traits inside `region`.
# @noRd
.jointSpecScope <- function(spec, data, dataForm, contexts, traitIds, region) {
    resolved <- .fmResolveSpecScope(
        spec,
        data,
        contexts = contexts,
        traitIds = traitIds
    )
    restrictToRegion <- dataForm == "individual" &&
        is.null(traitIds) &&
        !is.null(region)
    if (restrictToRegion) {
        .jointRestrictRegionTraits(resolved, data, region)
    } else {
        resolved
    }
}

# One study's traits narrowed to those inside `region`; left alone when the
# study has no context to read a phenotype from. A MultiStudyQtlDataset's
# individual-level study reads its own.
# @noRd
.jointStudyTraitsInRegion <- function(st, scope, data, region) {
    ctxs <- scope$contexts[[st]]
    if (length(ctxs) == 0L) {
        return(scope$traits[[st]])
    }
    study <- if (.isMultiStudyComponent(data, st)) data[[st]] else data
    se <- .qtlContextSe(study, context = ctxs[[1L]])
    .fmTraitsInRegion(
        se,
        intersect(scope$traits[[st]], rownames(se)),
        region
    )
}

# Region mode without an explicit traitId: restrict each study's scoped traits
# to the genes overlapping the locus (matches fineMappingPipeline). Gene coords
# are context-independent, so the first scoped context's SE supplies them.
# @noRd
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

# A QtlDataset is one study, so a study axis has nothing to join on it. Studies
# join on summary statistics, and fine-mapping also joins a
# MultiStudyQtlDataset's individual-level studies (.enumCrossStudyIndividual());
# a composed axis including the study needs summary statistics either way.
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
                "jointSpecification with axis 'study' needs more than one ",
                "study (QtlDataset is a single individual-level study)."
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

# Condition `r`'s TwasWeightsRow from column `column` of a joint weight
# matrix; NULL for a condition the fit has no column of.
# @noRd
.jointColEntry <- function(r, column, vids, weights, fitParts, cfg) {
    if (is.na(column)) {
        return(NULL)
    }
    TwasWeightsRow(
        variantIds = vids,
        weights = weights[, column],
        methodFits = fitParts,
        weightStandardized = TRUE
    )
}

# A joint-axis column holding a lone NA means "not a joint axis", so it is
# dropped rather than carried as a column of NA.
# @noRd
.jointAxisIsAbsent <- function(value) {
    !is.null(value) && length(value) == 1L && is.na(value)
}

# Build one record into a 1-row collection: wrap `entry` in a list and drop any
# NA joint* axis (so the union re-adds it, padded, only when some row is joint).
# @noRd
.jointBuildRecordPart <- function(rec, constructor) {
    jointCols <- c("jointStudies", "jointContexts", "jointTraits")
    absent <- jointCols[map_lgl(rec[jointCols], .jointAxisIsAbsent)]
    kept <- rec[setdiff(names(rec), c(absent, "study", "jointIndex"))]
    # `list_assign()`, not `list_modify()`: the entry payload is itself a
    # list, and list_modify would merge into it rather than wrap it. The
    # study goes in under the constructor's own name, not by partial
    # matching.
    constructorArgs <- list_assign(
        kept,
        studyName = rec$study,
        entry = list(rec$entry)
    )
    part <- exec(constructor, !!!constructorArgs)
    .withJointIndex(part, rec$jointIndex)
}

# `x`, one row of one method, recording its position in its joint fit's
# condition axis (`jointIndex`, internal). A row from no joint fit records
# none, so a collection holding no joint rows has no such column.
# @noRd
.withJointIndex <- function(x, jointIndex) {
    if (is.null(jointIndex) || is.na(jointIndex)) {
        return(x)
    }
    value <- if (methods::is(x, "TwasWeights")) {
        IRanges::IntegerList(jointIndex)
    } else {
        jointIndex
    }
    `$<-`(x, "jointIndex", value)
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
# none -> the column is omitted and traitPosition() reports an empty
# GRanges). region = the window the joint group was analysed over
# (.jointGroupWindow(): an explicit region, else the union of its traits' cis
# windows).
# @noRd
.jointEntryRecordAt <- function(
    i,
    entries,
    cond,
    method,
    grp,
    data,
    jointIndex
) {
    e <- entries[[i]]
    if (is.null(e)) {
        return(NULL)
    }
    st <- as.character(cond$study[[i]])
    ctx <- as.character(cond$context[[i]])
    tid <- as.character(cond$trait[[i]])
    tpos <- .traitPosFor(data, st, ctx, tid)
    list(
        study = st,
        context = ctx,
        trait = tid,
        method = method,
        entry = e,
        jointStudies = grp$js,
        jointContexts = grp$jc,
        jointTraits = grp$jt,
        traitPos = tpos,
        region = grp$region,
        # The reference the group was fitted against, so each result row
        # keeps it -- its own, when each study is on its own -- and NULL for
        # an individual-level group.
        ldReference = if (.jgIsPerConditionLd(grp$ldReference)) {
            grp$ldReference[[i]]
        } else {
            grp$ldReference
        },
        jointIndex = jointIndex[[i]]
    )
}

# Resume-cache lookup for condition `i` of joint group `joint` via the
# pipeline-specific `lookup` fn.
# @noRd
.jointCacheLookupAt <- function(i, lookup, cache, cond, token, joint) {
    lookup(
        cache,
        as.character(cond$study[[i]]),
        as.character(cond$context[[i]]),
        as.character(cond$trait[[i]]),
        token,
        joint
    )
}

# TRUE when a joint group has at least `minGroup` conditions.
# @noRd
.jointGroupMeetsMin <- function(g, minGroup) {
    nrow(.jgConditions(g)) >= minGroup
}
