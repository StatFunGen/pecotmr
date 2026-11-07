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
#              QtlSumStats -> the tuple's own `region`: the window it was
#              given, or the span of its variants as first built, before QC
#              (a pipeline's cisWindow sets it from traitPos first).
# Each returns a length-1 GRanges, or NULL when the anchor is unavailable (the
# accumulator / builder then records a chrUn sentinel). A row is named by its
# whole (study, context, trait) tuple: a QtlSumStats can carry the same
# (context, trait) in several studies, each with its own position and window.
# A QtlDataset holds one study, so `study` has nothing to choose there.
#' @importFrom rlang try_fetch
.traitPosFor <- function(data, study, context, trait) {
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
# @noRd
.jointAfForGroup <- function(g, block) {
    if (!is(g, "IndividualJointGroup") || is.null(block$af)) {
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
        cfg = cfg,
        af = args$af,
        region = args$region
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
        args$mrmashPrior,
        key$study,
        key$trait,
        context = key$context
    )
    mvCv <- .fmLookupMrmashCv(
        args$mrmashPrior,
        key$study,
        key$trait,
        context = key$context
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
# the condition's position in the fitted survivor set + its CV slice.
# @noRd
.jointMvEntry <- function(
    i,
    fitted,
    Xc,
    keep,
    survivors,
    cfg,
    credibleSetParam,
    fitRetention,
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
            dataX = Xc,
            dataY = NULL,
            af = af,
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
            nCond,
            token
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
            ldReference = .jgLdRef(group),
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
        return(TwasWeightsRow(
            variantIds = character(0),
            weights = NULL,
            cvResult = cvR,
            weightStandardized = stdz
        ))
    }
    TwasWeightsRow(
        variantIds = full$vids,
        weights = full$W[, r],
        methodFits = full$fitParts,
        cvResult = cvR,
        weightStandardized = stdz
    )
}

# (sumstats, twas) -> mr.mash-rss joint fit as ONE matrix entry. No sample
# folds.
setMethod(
    "fitJointGroup",
    signature("SumStatsJointGroup", "TwasJointPipeline"),
    function(group, pipeline, token, args) {
        cfg <- .jpConfig(pipeline)
        weights <- .jointMrmashRssWeights(group, cfg, args)
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

# cross-study / sumstats (no individual form: individual-level studies have
# disjoint samples). One group per (context, trait) present in >= 2 scoped
# studies; the study axis varies -> "joint" + jointStudies.
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
# just that study, so the cell contributes nothing.
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
    ldReference <- .jointRowsLdRef(data, tupleRows, "jointCrossStudy")
    stNames <- cols$study[tupleRows]
    jz <- .buildJointSumstatZMatrix(
        data,
        tupleRows,
        stNames,
        errorLabel = "jointCrossStudy",
        ldReference = ldReference,
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
        ldReference = ldReference,
        N = jz$nVec
    )
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
.twasGroupArgs <- function(g, pipeline, args) {
    if (!is(pipeline, "TwasJointPipeline")) {
        return(args)
    }
    cfg <- .jpConfig(pipeline)
    cond <- .jgConditions(g)
    groupFmArgs <- .twasGroupFmArgs(args, cond)
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

# The group's fine-mapping fits, read off its first condition's row (every
# condition's row of a joint fit keys the one fit), and its CV. Empty when
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
    af <- .twasFineMappingFits(fmRes, studyName = s1, context = c1, trait = t1)
    list(
        fittedModels = if (is.null(af)) {
            list()
        } else {
            .twasFitsForRegion(af, bi, nR)
        },
        fineMappingCv = .twasGroupFmCv(fmRes, cond)
    )
}

# The group's fine-mapping CV with one prediction column and one performance
# row per condition, in the group's order -- the shape .jointTwasCvResult()
# gives, which the per-condition slice indexes. A joint fit stores each
# condition's slice on that condition's row, so the first row alone covered
# only the first condition. A method is handed off only when every condition
# carries its CV, and only from one fold partition.
# @noRd
.twasGroupFmCv <- function(fmRes, cond) {
    cvs <- pmap(
        list(
            studyName = as.character(cond$study),
            context = as.character(cond$context),
            trait = as.character(cond$trait)
        ),
        .twasCvResultFor,
        fineMappingResult = fmRes
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

# The genotypes a unit's fits are on; none for a group already built, which
# holds its own.
# @noRd
.unitSpec <- function(unit, data, args) {
    if (is(unit, "JointGroup")) {
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
    built <- if (is(unit, "JointGroup")) {
        unit
    } else {
        unitArgs <- unit$args
        exec(unit$builder, !!!unitArgs, X = block$X)
    }
    if (is(built, "JointGroup")) {
        return(list(built))
    }
    compact(as.list(built))
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
    fitArgs <- .jointGroupFitArgs(g, pipeline, args, block, tokens)
    perTokenEntries <- compact(set_names(
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

# A group's provenance: the ";"-joined members of each varying axis (identical
# on every per-context row of this group), and the window this group was
# analysed over -- the caller's region when one was named, else the group's own
# cis window (the union span of its traits +/- cisWindow). The window is
# recorded because the fit only ever sees whichever variants survived QC; it is
# read back by genomicRegion().
# @noRd
.jointGroupProvenance <- function(g, data, args) {
    cond <- .jgConditions(g)
    list(
        cond = cond,
        js = .jointAxisMembers(cond, "study"),
        jc = .jointAxisMembers(cond, "context"),
        jt = .jointAxisMembers(cond, "trait"),
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
        .twasGroupArgs(g, pipeline, args),
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
        cfg = cfg
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
        rsq = set_names(map_dbl(parts, "rsq"), contributing),
        partition = parts[[1L]]$partition
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
        variantIds(perTokenEntries[[passing[1L]]][[r]]),
        stdz
    )
}

# The ensemble's row: its combined weights, and a cvResult recording its
# cross-fit out-of-fold predictions and metrics over its methods' folds, as
# each method's does, beside the stacking coefficients and the methods' own
# R^2. `fallbackIds` names the weights when the ensemble's do not.
# @noRd
.twasEnsembleRow <- function(ens, partition, fallbackIds, stdz) {
    ew <- ens$ensembleTwasWeights
    vids <- (if (!is.null(names(ew))) names(ew) else rownames(ew)) %||%
        fallbackIds
    TwasWeightsRow(
        variantIds = vids,
        weights = as.numeric(ew),
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
    .jointCellWork(cell, pipeline, data, scope, jointMethods, spArgs)
}

# One study's traits narrowed to those inside `region`; left alone when the
# study has no context to read a phenotype from.
# @noRd
.jointStudyTraitsInRegion <- function(st, scope, data, region) {
    ctxs <- scope$contexts[[st]]
    if (length(ctxs) == 0L) {
        return(scope$traits[[st]])
    }
    se <- .qtlContextSe(data, context = ctxs[[1L]])
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
    TwasWeightsRow(
        variantIds = vids,
        weights = weights[, r],
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
        # keeps it; NULL for an individual-level group.
        ldReference = grp$ldReference,
        jointIndex = jointIndex[[i]]
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
