#' @include MethodParam.R
NULL

# =============================================================================
# Cross-validation engine (shared by twasWeightsCv + fine-mapping CV)
# -----------------------------------------------------------------------------
# A single K-fold CV harness backs both the TWAS predictive-weight CV
# (twasWeightsCv, twasWeights.R) and the fine-mapping CV (.fmWeightsCv,
# fineMappingPipeline.R). Callers supply only a per-fold fit; partitioning,
# the parallel fold loop, prediction, aggregation, the metric block, and the
# output key format live here so the two paths cannot drift.
# =============================================================================

# =============================================================================
# Shared cross-validation engine
# =============================================================================
# One CV harness backs both twasWeightsCv() (predictive weight methods) and the
# fine-mapping CV (.fmWeightsCv, susie/mvsusie/fsusie). Callers differ only in
# the per-fold fit; partitioning, the parallel fold loop, prediction, fold
# aggregation, the metric block, and the output key format live here.

# Sample/Fold partition: shuffle samples, then cut into `fold` contiguous
# blocks.
# @noRd
.cvSamplePartition <- function(sampleNames, fold) {
    idx <- sample(length(sampleNames))
    folds <- cut(seq_along(sampleNames), breaks = fold, labels = FALSE)
    tibble(
        Sample = sampleNames[idx],
        Fold = folds
    )
}

# Canonical CV output key: "<methodKey>_predicted" / "<methodKey>_performance".
# @noRd
.cvOutputKey <- function(methodKey, suffix) str_c(methodKey, "_", suffix)

# One CV metric row (corr, rsq, adj_rsq, pval, RMSE, MAE) for a single outcome.
# Drops NA in either vector; needs >= 3 valid points and non-constant
# predictions.
# @noRd
.cvMetricRow <- function(pred, actual) {
    out <- set_names(
        rep(NA_real_, 6L),
        c("corr", "rsq", "adj_rsq", "pval", "RMSE", "MAE")
    )
    ok <- !is.na(pred) & !is.na(actual)
    pred <- pred[ok]
    actual <- actual[ok]
    if (length(pred) < 3L || stats::sd(pred) == 0) {
        return(out)
    }
    lmFit <- stats::lm(actual ~ pred)
    s <- summary(lmFit)
    res <- actual - pred
    c(
        corr = stats::cor(actual, pred),
        rsq = s$r.squared,
        adj_rsq = s$adj.r.squared,
        pval = if (nrow(s$coefficients) >= 2L) {
            s$coefficients[2L, 4L]
        } else {
            NA_real_
        },
        RMSE = sqrt(mean(res^2)),
        MAE = mean(abs(res))
    )
}

# A canonical fingerprint of a Sample/Fold partition, used to prove that
# per-fold fits handed to twasWeightsCv were trained on the same folds it is
# scoring. Sample order is normalised first, so two partitions that assign the
# same samples to the same folds agree regardless of row order.
# @noRd
.cvPartitionKey <- function(samplePartition) {
    if (is.null(samplePartition)) {
        return(NULL)
    }
    tbl <- samplePartition |>
        arrange(.data$Sample, .data$Fold) |>
        mutate(.key = str_c(.data$Sample, "=", .data$Fold)) |>
        pull(".key")
    rlang::hash(tbl)
}

# One CV fold: split train/test by fold `j`, drop zero-variance training
# columns,
# fit via `fitFold(Xtr, Ytr, j, fitFoldCtx)`, and predict the held-out samples.
# @noRd
.cvRunFold <- function(j, cv) {
    X <- cv$X
    Y <- cv$Y
    samplePartition <- cv$samplePartition
    foldIds <- cv$foldIds
    fitFold <- cv$fitFold
    fitFoldCtx <- cv$fitFoldCtx
    fitRetention <- cv$fitRetention
    verbose <- cv$verbose
    if (verbose >= 1) {
        msg <- glue("  CV fold {j}/{length(foldIds)} ...")
        inform(msg)
    }
    testIds <- samplePartition |> filter(.data$Fold == j) |> pull("Sample")
    isTest <- is_in(rownames(X), testIds)
    if (all(isTest) || !any(isTest)) {
        return(list(preds = list(), fits = list()))
    }
    trainAll <- X[!isTest, , drop = FALSE]
    Xte <- X[isTest, , drop = FALSE]
    Ytr <- Y[!isTest, , drop = FALSE]
    keep <- .nonzeroVarColumns(trainAll)
    Xtr <- trainAll[, keep, drop = FALSE]
    ff <- fitFold(Xtr, Ytr, j, fitFoldCtx)
    preds <- map(ff$weights, .cvFoldPrediction, Xte = Xte)
    list(
        preds = preds,
        weights = if (cv$keepFoldWeights) ff$weights,
        fits = if (identical(fitRetention, "none")) list() else ff$fits
    )
}

# No-op fold fitter: used when the caller only wants the fold partition.
# @noRd
.cvNoopFitFold <- function(Xtr, Ytr, j, fitFoldCtx) {
    list(weights = list(), fits = list())
}

# Shared K-fold cross-validation engine.
#
# `fitFold(Xtrain, Ytrain, foldIndex, fitFoldCtx)` must return
#   list(weights = <named list: methodKey -> (variants x outcomes) weight
#   matrix,
#                   rownames indexing colnames(Xtrain)>,
#        fits    = <named list: methodKey -> fitted model or NULL>)
# The engine drops zero-variance training columns before calling fitFold,
# predicts held-out samples (Xtest[, common] %*% W[common, ]), aggregates the
# out-of-fold predictions into full-Y matrices, scores them with .cvMetricRow,
# and keys the output via .cvOutputKey(). `maxNumVariants` (optional) randomly
# subsamples variants up front to bound compute; `BPPARAM` runs the folds
# (SerialParam() one at a time). `keepFoldWeights` adds each fold's weights
# to the output as `foldWeights`, for a caller that stores them.
#' @importFrom BiocParallel SerialParam
#' @importFrom stats sd lm cor
#' @importFrom dplyr n_distinct
#' @noRd
.crossValidateWeights <- function(
    X,
    Y,
    fold = NULL,
    samplePartitions = NULL,
    fitFold,
    fitFoldCtx = NULL,
    BPPARAM = SerialParam(),
    maxNumVariants = NULL,
    variantsToKeep = NULL,
    fitRetention = c("none", "slim", "full"),
    verbose = 1,
    seed = NULL,
    keepFoldWeights = FALSE
) {
    fitRetention <- arg_match(fitRetention)
    .applySeed(seed)
    prep <- .cvPrepareData(X, Y, fold, verbose)
    X <- .cvSubsampleVariants(prep$X, maxNumVariants, variantsToKeep, verbose)
    Y <- prep$Y
    samplePartition <- .cvResolvePartition(
        samplePartitions,
        fold,
        rownames(Y),
        verbose
    )
    foldIds <- sort(unique(samplePartition$Fold))
    st <- proc.time()
    cvState <- list(
        X = X,
        Y = Y,
        samplePartition = samplePartition,
        foldIds = foldIds,
        fitFold = fitFold,
        fitFoldCtx = fitFoldCtx,
        fitRetention = fitRetention,
        keepFoldWeights = keepFoldWeights,
        verbose = verbose,
        rngSeed = seed
    )
    foldResults <- .cvRunFolds(foldIds, cvState, BPPARAM)
    .cvOutput(cvState, foldResults, st)
}

# The engine's result from its folds': the partition, the out-of-fold
# predictions and their scores, the fits kept, and each fold's weights when
# the caller asked for them. `st` is when the folds started.
# @noRd
.cvOutput <- function(cvState, foldResults, st) {
    agg <- .cvAggregate(foldResults, cvState$Y, cvState$verbose)
    foldIds <- cvState$foldIds
    out <- list(
        samplePartition = cvState$samplePartition,
        prediction = agg$prediction,
        performance = agg$performance,
        foldFits = .cvCollectFoldFits(foldResults, foldIds),
        timeElapsed = proc.time() - st
    )
    if (!cvState$keepFoldWeights) {
        return(out)
    }
    list_assign(out, foldWeights = .cvCollectFoldWeights(foldResults, foldIds))
}

# Validate inputs, coerce a vector Y to a one-column matrix, and set stable
# row/column dimnames on X and Y. Returns list(X, Y).
# @noRd
#' @importFrom checkmate assert assertCount assertMatrix
#' @importFrom checkmate checkAtomicVector checkMatrix
.cvPrepareData <- function(X, Y, fold, verbose) {
    assertCount(fold, positive = TRUE, null.ok = TRUE)
    assertMatrix(X)
    assert(checkMatrix(Y), checkAtomicVector(Y), .var.name = "Y")
    if (is.vector(Y)) {
        Y <- matrix(Y, ncol = 1)
        if (verbose >= 1) {
            msg <- glue(
                "Y converted to matrix of {nrow(Y)} rows and {ncol(Y)} ",
                "columns."
            )
            inform(msg)
        }
    }
    assertMatrix(Y, nrows = nrow(X))
    .cvSetDimnames(X, Y)
}

# @noRd
.cvSetDimnames <- function(X, Y) {
    sampleNames <- if (!is.null(rownames(Y))) {
        rownames(Y)
    } else if (!is.null(rownames(X))) {
        rownames(X)
    } else {
        str_c("sample_", seq_len(nrow(X)))
    }
    list(
        X = `dimnames<-`(
            X,
            list(
                rownames(X) %||% sampleNames,
                colnames(X) %||% str_c("variable_", seq_len(ncol(X)))
            )
        ),
        Y = `dimnames<-`(
            Y,
            list(
                rownames(Y) %||% sampleNames,
                colnames(Y) %||% str_c("context_", seq_len(ncol(Y)))
            )
        )
    )
}

# Optional variant subsample (compute saver; e.g. TWAS weight CV).
# @noRd
.cvSubsampleVariants <- function(X, maxNumVariants, variantsToKeep, verbose) {
    if (is.null(maxNumVariants) || ncol(X) <= maxNumVariants) {
        return(X)
    }
    selected <- .cvSelectVariants(X, maxNumVariants, variantsToKeep, verbose)
    X[, selected, drop = FALSE]
}

# Choose the variant subset: honor variantsToKeep (topping up at random when
# below the cap, down-sampling when above), else a plain random draw.
# @noRd
.cvSelectVariants <- function(X, maxNumVariants, variantsToKeep, verbose) {
    if (is.null(variantsToKeep) || length(variantsToKeep) == 0) {
        selected <- sort(sample(ncol(X), maxNumVariants))
        .cvMsgRandom(verbose, length(selected), ncol(X))
        return(selected)
    }
    variantsToKeep <- intersect(variantsToKeep, colnames(X))
    if (length(variantsToKeep) >= maxNumVariants) {
        selected <- sample(variantsToKeep, maxNumVariants)
        .cvMsgRandomKept(verbose, length(selected), length(variantsToKeep))
        return(selected)
    }
    remaining <- setdiff(colnames(X), variantsToKeep)
    additional <- sample(remaining, maxNumVariants - length(variantsToKeep))
    selected <- union(variantsToKeep, additional)
    .cvMsgKeptPlus(
        verbose,
        length(variantsToKeep),
        length(additional),
        length(selected),
        ncol(X)
    )
    selected
}

# @noRd
.cvMsgRandom <- function(verbose, nSelected, nTotal) {
    if (verbose >= 1) {
        msg <- glue(
            "Randomly selecting {nSelected} out of {nTotal} variants for ",
            "cross validation purpose."
        )
        inform(msg)
    }
}

# @noRd
.cvMsgRandomKept <- function(verbose, nSelected, nKept) {
    if (verbose >= 1) {
        msg <- glue(
            "Randomly selecting {nSelected} out of {nKept} input variants ",
            "for cross validation purpose."
        )
        inform(msg)
    }
}

# @noRd
.cvMsgKeptPlus <- function(verbose, nKept, nAdditional, nSelected, nTotal) {
    if (verbose >= 1) {
        msg <- glue(
            "Including {nKept} specified variants and randomly selecting ",
            "{nAdditional} additional variants, for a total of {nSelected} ",
            "variants out of {nTotal} for cross-validation purpose."
        )
        inform(msg)
    }
}

# Reuse a provided fold partition (validated) or build one from `fold`.
# @noRd
.cvResolvePartition <- function(samplePartitions, fold, sampleNames, verbose) {
    if (!is.null(samplePartitions)) {
        return(.cvValidatePartition(
            samplePartitions,
            fold,
            sampleNames,
            verbose
        ))
    }
    if (!is.null(fold)) {
        return(.cvSamplePartition(sampleNames, fold))
    }
    abort("Either 'fold' or 'samplePartitions' must be provided.")
}

# @noRd
.cvValidatePartition <- function(samplePartitions, fold, sampleNames, verbose) {
    if (!all(is_in(samplePartitions$Sample, sampleNames))) {
        abort("Some samples in 'samplePartitions' do not match 'X' and 'Y'.")
    }
    nF <- n_distinct(samplePartitions$Fold)
    if (!is.null(fold) && verbose >= 1 && fold != nF) {
        msg <- glue(
            "fold number provided does not match with sample partition, ",
            "performing {nF} fold cross validation based on provided sample ",
            "partition. "
        )
        inform(msg)
    }
    samplePartitions
}

# ---- BiocParallel ------------------------------------------------------------

# Refuse a `BPPARAM` that is not a BiocParallelParam.
# @noRd
.assertBpparam <- function(BPPARAM, what = "BPPARAM") {
    if (methods::is(BPPARAM, "BiocParallelParam")) {
        return(invisible(BPPARAM))
    }
    abort(glue(
        "`{what}` must be a BiocParallelParam, such as SerialParam() or ",
        "MulticoreParam(workers = 4) (got {class(BPPARAM)[[1L]]})."
    ))
}

# bplapply() on a copy of the caller's BPPARAM seeded by .seededBpparam(). A
# worker's error comes back as BiocParallel's report on every element; the
# first one is raised again with its own message, as a serial map would have
# raised it. The extra arguments reach `.f` through bplapply(), whose own
# formals include X and FUN, so pass an argument of either name by position.
# @noRd
#' @importFrom BiocParallel bplapply
.bpMap <- function(.x, .f, ..., BPPARAM, seed = NULL) {
    try_fetch(
        bplapply(.x, .f, ..., BPPARAM = .seededBpparam(BPPARAM, seed)),
        bplist_error = .bpRaiseFirstError
    )
}

# @noRd
.bpRaiseFirstError <- function(cnd) {
    failed <- keep(attr(cnd, "result"), .isRemoteError)
    if (length(failed) == 0L) {
        cnd_signal(cnd)
    }
    abort(conditionMessage(failed[[1L]]))
}

# @noRd
.isRemoteError <- function(x) {
    methods::is(x, "remote_error")
}

# The caller's BPPARAM, copied and given the RNG seed its workers draw from:
# `seed` when one is given, else the param's own RNGseed, else one draw from
# the session RNG, so a set.seed() before the call still governs the run. A
# copy is seeded because a BiocParallelParam is a reference object: seeding
# it in place would change the caller's own. Each element then draws from its
# own stream, so the numbers do not depend on the backend or worker count.
# @noRd
#' @importFrom BiocParallel bpRNGseed bpRNGseed<-
.seededBpparam <- function(BPPARAM, seed = NULL) {
    copy <- BPPARAM$copy()
    bpRNGseed(copy) <- as.integer(
        seed %||% bpRNGseed(BPPARAM) %||% .sessionSeedDraw()
    )
    copy
}

# One integer from the session RNG.
# @noRd
.sessionSeedDraw <- function() {
    sample.int(.Machine$integer.max, 1L)
}

# The threads an engine with its own parallelism (OpenMP, a C++ thread pool,
# mclapply) may use under BPPARAM: its worker count when the workers are
# processes on this machine (SerialParam, MulticoreParam), else 1. A
# SnowParam's or BatchtoolsParam's workers can be other machines or
# scheduler jobs, so their count is no budget for threads here. Said once a
# session per backend, as a windowed engine asks once per window.
# @noRd
#' @importFrom BiocParallel bpnworkers
.bpLocalThreads <- function(BPPARAM, caller) {
    if (
        methods::is(BPPARAM, "SerialParam") ||
            methods::is(BPPARAM, "MulticoreParam")
    ) {
        return(as.integer(bpnworkers(BPPARAM)))
    }
    backend <- class(BPPARAM)[[1L]]
    inform(
        glue(
            "{caller}: the workers of a {backend} need not run on this ",
            "machine, so a threaded engine runs on 1 thread."
        ),
        .frequency = "once",
        .frequency_id = str_c("pecotmr-bpLocalThreads-", backend)
    )
    1L
}

# Seed the main-process RNG when a seed is supplied. This covers only the serial
# draws (fold partitioning, variant sub-sampling); BiocParallel workers draw
# from the streams .seededBpparam() seeds. A NULL seed leaves the session RNG
# untouched so an outer set.seed() still governs it.
#
# Scoped to the CALLER's frame rather than set.seed()'s permanent mutation, so
# a seeded call does not leave the user's session RNG changed after it returns.
# Bioconductor also asks packages not to call set.seed() directly.
# @noRd
.applySeed <- function(seed) {
    if (!is.null(seed)) {
        withr::local_seed(seed, .local_envir = parent.frame())
    }
}

# Run each fold on BPPARAM's workers.
# @noRd
.cvRunFolds <- function(foldIds, cvState, BPPARAM) {
    .bpMap(
        foldIds,
        .cvRunFold,
        cv = cvState,
        BPPARAM = BPPARAM,
        seed = cvState$rngSeed
    )
}

# Aggregate per-fold predictions + performance across methods. Returns
# list(prediction, performance).
# @noRd
.cvAggregate <- function(foldResults, Y, verbose) {
    metricNames <- c("corr", "rsq", "adj_rsq", "pval", "RMSE", "MAE")
    methodKeys <- unique(list_c(map(foldResults, .cvPredNames)))
    predMats <- map(methodKeys, .cvPredMatrix, foldResults = foldResults, Y = Y)
    list(
        prediction = set_names(
            predMats,
            map_chr(methodKeys, .cvOutputKey, suffix = "predicted")
        ),
        performance = set_names(
            map2(
                predMats,
                methodKeys,
                .cvPerformanceFor,
                Y = Y,
                verbose = verbose,
                metricNames = metricNames
            ),
            map_chr(methodKeys, .cvOutputKey, suffix = "performance")
        )
    )
}

# @noRd
.cvPerformanceFor <- function(predMat, mk, Y, verbose, metricNames) {
    .cvPerformance(predMat, Y, mk, verbose, metricNames)
}

# @noRd
.cvPredNames <- function(fr) {
    names(fr$preds)
}

# Assemble the (samples x conditions) prediction matrix for a method, scattering
# each fold's held-out predictions by row name.
# @noRd
.cvPredMatrix <- function(mk, foldResults, Y) {
    # Deliberate scatter: the folds partition the samples, so each writes its
    # own held-out rows once and samples in no fold stay NA.
    predMat <- matrix(NA_real_, nrow(Y), ncol(Y), dimnames = dimnames(Y))
    for (yh in compact(map(foldResults, .cvFoldPreds, mk = mk))) {
        predMat[rownames(yh), ] <- yh
    }
    predMat
}

# @noRd
.cvFoldPreds <- function(fr, mk) {
    fr$preds[[mk]]
}

# Per-condition performance metrics (conditions x metrics).
# @noRd
.cvPerformance <- function(predMat, Y, mk, verbose, metricNames) {
    metricRows <- map(
        seq_len(ncol(Y)),
        .cvConditionMetrics,
        predMat = predMat,
        Y = Y,
        mk = mk,
        verbose = verbose
    )
    `dimnames<-`(
        exec(rbind, !!!metricRows),
        list(colnames(Y), metricNames)
    )
}

# @noRd
.cvConditionMetrics <- function(r, predMat, Y, mk, verbose) {
    .cvWarnZeroVariance(predMat[, r], mk, r, verbose)
    .cvMetricRow(predMat[, r], Y[, r])
}

# @noRd
.cvWarnZeroVariance <- function(pred, mk, r, verbose) {
    if (verbose < 1) {
        return(invisible(NULL))
    }
    s <- stats::sd(pred, na.rm = TRUE)
    if (is.na(s) || s == 0) {
        msg <- glue(
            "Predicted values for condition {r} using {mk} have zero ",
            "variance. Filling performance metric with NAs"
        )
        inform(msg)
    }
    invisible(NULL)
}

# Non-NULL per-fold fits keyed by fold id.
# @noRd
.cvCollectFoldFits <- function(foldResults, foldIds) {
    set_names(map(foldResults, .cvCompactFits), str_c("fold_", foldIds))
}

# @noRd
.cvCompactFits <- function(fr) {
    compact(fr$fits)
}

# Every fold's weight matrices (method key -> variants x outcomes, rows named
# by the fold's training variants), keyed fold_<j> like the fits.
# @noRd
.cvCollectFoldWeights <- function(foldResults, foldIds) {
    set_names(map(foldResults, "weights"), str_c("fold_", foldIds))
}

# Test-fold prediction from one method's weight matrix W: align to the test
# genotype columns, then Xte %*% W. NULL for a missing/non-overlapping W.
# @noRd
.cvFoldPrediction <- function(W, Xte) {
    if (is.null(W)) {
        return(NULL)
    }
    W <- replace(W, is.na(W), 0)
    common <- intersect(colnames(Xte), rownames(W))
    if (length(common) == 0L) {
        return(NULL)
    }
    `rownames<-`(
        Xte[, common, drop = FALSE] %*% W[common, , drop = FALSE],
        rownames(Xte)
    )
}

# =============================================================================
# Cross-validation settings
# -----------------------------------------------------------------------------
# One bundle for both fineMappingPipeline() and twasWeightsPipeline(). The two
# had drifted on a setting that means the same thing in each: `cvFolds`
# defaulted to 0 in fine-mapping and 5 in TWAS. Folding them into one
# constructor is the rule the joint specification already follows.
#
# `folds` is 0 everywhere. twasWeightsPipeline used to default to 5, so the
# two pipelines disagreed on what an unspecified `cvFolds` meant; converging
# on "off unless asked" makes one bundle mean one thing. The knock-on is that
# the SR-TWAS ensemble, which needs out-of-fold predictions, no longer runs by
# default -- so `EnsembleParam(enabled =)` defaults to FALSE too, and asking
# for an ensemble
# without folds is an error rather than a silently missing row.
#
# `fitFullData` means the same thing in both pipelines: keep only the
# cross-validation results, skipping the full-data fit.
#
# Two of the six fields are honoured by twasWeightsPipeline() only, and are
# IGNORED rather than rejected by fineMappingPipeline() -- one bundle the
# caller can hand to either pipeline without rewriting it:
#
#   weightMethods -- selects among TWAS *weight* methods, matching tokens like
#     `lasso` / `mrmash` and their `<token>_weights` spellings. Fine mapping's
#     `methods=` are fine-mapping methods and its CV refits all of them, so
#     there is no weight-method set to select from. It is also the knob that
#     decides whether .jointTwasCvBlocked() reaches back for a fold's
#     fine-mapping fit -- a TWAS-to-fine-mapping handoff with no counterpart
#     in the other direction.
#   maxVariants -- caps the CV design matrix for twasWeightsCv(). Fine mapping
#     does not use that engine; .fmWeightsCv() has its own loop and only
#     mirrors twasWeightsCv()'s output shape.
#
# Both are documented as TWAS-only so the silence is stated, not discovered.
#
# `seed` is not here either. It seeds the whole call in both pipelines
# (fine-mapping wraps the call in withr::local_seed; TWAS also seeds the
# BiocParallel RNG for method fitting), so it is a call-level knob that CV
# happens to benefit from, not a CV setting.
# =============================================================================

#' @rdname CrossValidationParam
#' @aliases CrossValidationParam-class
#' @importClassesFrom BiocParallel BiocParallelParam
#' @exportClass CrossValidationParam
setClass(
    "CrossValidationParam",
    contains = "MethodParam",
    slots = c(
        folds = "numeric",
        BPPARAM = "BiocParallelParam",
        samplePartition = "list_OR_NULL",
        maxVariants = "numeric_OR_NULL",
        weightMethods = "list_OR_NULL",
        fitFullData = "logical"
    )
)

setValidity("CrossValidationParam", function(object) {
    problems <- .cvFitFullDataProblem(object@fitFullData, object@folds)
    if (length(problems) == 0L) TRUE else problems
})

# `fitFullData = FALSE` asks for cross-validation results alone, which is only
# a request for something when cross-validation runs: without folds there
# would be no fit of any kind. Checked here, at construction, rather than
# separately by each pipeline that reads it.
# @noRd
.cvFitFullDataProblem <- function(fitFullData, folds) {
    if (length(fitFullData) != 1L || is.na(fitFullData)) {
        return("`fitFullData` must be a single TRUE or FALSE.")
    }
    if (!fitFullData && !.cvEnabled(list(folds = folds))) {
        return(str_c(
            "fitFullData = FALSE requires cross-validation (folds >= 2): ",
            "without folds there would be no fit at all."
        ))
    }
    NULL
}

#' @title Cross-Validation Settings
#' @description How cross-validation runs, shared by
#'   \code{\link{fineMappingPipeline}} and \code{\link{twasWeightsPipeline}}.
#'
#'   \code{folds}, \code{BPPARAM}, \code{samplePartition} and
#'   \code{fitFullData} mean the same thing in both. \code{maxVariants} and
#'   \code{weightMethods} are honoured by \code{\link{twasWeightsPipeline}}
#'   only and are \strong{ignored} by \code{\link{fineMappingPipeline}}, so
#'   one bundle can be handed to either pipeline unchanged.
#' @param folds Integer. Number of folds; \code{0} (the default) or
#'   \code{1} skips cross-validation. The same value in both pipelines ---
#'   \code{twasWeightsPipeline} used to default to \code{5}, which made an
#'   unspecified setting mean two different things.
#' @param BPPARAM A \code{BiocParallelParam} the per-fold refits run on;
#'   \code{SerialParam()} (default) runs them one at a time. Its workers draw
#'   from random-number streams seeded by the pipeline's \code{seed}, else by
#'   the param's own \code{RNGseed}, else by one draw from the session RNG,
#'   so a run gives the same numbers however many workers it has. Only
#'   consulted when \code{folds > 1}.
#' @param samplePartition Optional pre-defined partition \code{data.frame}
#'   with columns \code{Sample} and \code{Fold}. When supplied, every method
#'   reuses this exact partition instead of generating a fresh one.
#' @param maxVariants Integer. Cap on the number of variants used for CV;
#'   unset means no limit. \strong{\code{twasWeightsPipeline} only} ---
#'   \code{fineMappingPipeline} ignores it, because it does not use the
#'   \code{twasWeightsCv} engine this configures.
#' @param weightMethods Optional override of which methods are
#'   cross-validated, as a character vector of tokens or a named method list.
#'   Unset cross-validates every method that produced non-zero weights.
#'   \strong{\code{twasWeightsPipeline} only} ---
#'   \code{fineMappingPipeline} ignores it, because its \code{methods=} are
#'   fine-mapping methods and its CV refits all of them.
#' @param fitFullData Logical. \code{TRUE} (default) fits the final model on
#'   the full data as well as cross-validating it; \code{FALSE} keeps only
#'   the cross-validation results, and requires \code{folds >= 2}. A row
#'   from such a run carries its \code{cvResult} and nothing else, so it
#'   cannot be exported and cannot supply full-data weights or fits
#'   downstream.
#' @return A \code{CrossValidationParam} object, a \code{\link{MethodParam}}.
#' @examples
#' CrossValidationParam(folds = 10, BPPARAM = BiocParallel::SerialParam())
#' CrossValidationParam(folds = 5, fitFullData = FALSE)
#' @export
CrossValidationParam <- function(
    folds = 0,
    BPPARAM = SerialParam(),
    samplePartition = NULL,
    maxVariants = NULL,
    weightMethods = NULL,
    fitFullData = TRUE
) {
    new(
        "CrossValidationParam",
        folds = folds,
        BPPARAM = .assertBpparam(BPPARAM),
        samplePartition = samplePartition,
        maxVariants = maxVariants,
        weightMethods = weightMethods,
        fitFullData = fitFullData
    )
}

# The fields a pipeline reads, as a plain list so call sites index it without
# caring that the caller passed a constructor result.
# @noRd
.cvResolve <- function(crossValidationParam) {
    list(
        folds = crossValidationParam$folds %||% 0,
        BPPARAM = crossValidationParam$BPPARAM %||% SerialParam(),
        samplePartition = crossValidationParam$samplePartition,
        # -1 is twasWeightsCv()'s "no cap" sentinel. fineMappingPipeline
        # resolves these two as well but never reads them.
        maxVariants = crossValidationParam$maxVariants %||% -1,
        weightMethods = crossValidationParam$weightMethods,
        fitFullData = crossValidationParam$fitFullData %||% TRUE
    )
}

# TRUE when cross-validation will actually run. Both 0 and 1 mean "no CV",
# and anything that needs out-of-fold predictions has to ask this rather than
# testing `folds` itself.
# @noRd
.cvEnabled <- function(crossValidationParam) {
    folds <- crossValidationParam$folds %||% 0
    !is.null(folds) && length(folds) == 1L && !is.na(folds) && folds >= 2L
}

# Refuse a cross-validation request on an input that has no samples to hold
# out. Summary statistics carry no individual-level data, so CV there is not
# unimplemented but meaningless -- saying so beats running to completion and
# returning results that were never cross-validated.
# @noRd
.cvRefuseOnSumstats <- function(crossValidationParam, pipeline, cls) {
    set <- names(crossValidationParam)[
        map_lgl(as.list(crossValidationParam), .cvFieldIsSet)
    ]
    if (length(set) == 0L) {
        return(invisible(NULL))
    }
    abort(glue(
        "{pipeline}: cross-validation is not possible on {cls} input -- ",
        "it holds out samples, and summary statistics carry none. ",
        "Remove CrossValidationParam({str_flatten(set, ', ')})."
    ))
}

# A field counts as "set" when it differs from the constructor's own default,
# so passing CrossValidationParam() unchanged is not mistaken for a request.
# fitFullData = TRUE and BPPARAM = SerialParam() are those fields' defaults.
# @noRd
.cvFieldIsSet <- function(x) {
    if (is.null(x) || isTRUE(x) || methods::is(x, "SerialParam")) {
        return(FALSE)
    }
    if (methods::is(x, "BiocParallelParam")) {
        return(TRUE)
    }
    # folds = 0 (or 1) is the no-op default.
    !(length(x) == 1L && !is.na(x) && is.numeric(x) && (x == 0 || x == 1))
}

# =============================================================================
# Cross-validation results in a collection
# -----------------------------------------------------------------------------
# A CV run's rows -- one per method and condition -- share its fold partition
# and the order its out-of-fold predictions are held in. A collection keeps
# each such pair once, in its cvPartition dictionary (see the row dictionaries
# in tupleSelectors.R); every row keeps its predictions without sample names
# and a key to the pair it shares. The rows a collection hands out get both
# back (.cvJoinPartition()), so they carry the CV their payload was built
# with.
# =============================================================================

# @noRd
.CV_PARTITION_DICT <- list(
    slot = "cvPartition",
    keyCol = "cvPartition",
    prefix = "cv"
)

# A collection's shared CV partitions; none for one saved before the slot
# existed, whose rows still hold their own.
# @noRd
.cvPartitionDict <- function(x) {
    if (!methods::.hasSlot(x, "cvPartition")) {
        return(list())
    }
    x@cvPartition
}

# What is wrong with a collection's CV partitions and its rows' keys into
# them, or NULL.
# @noRd
.cvPartitionProblems <- function(x) {
    .dictProblems(
        .cvPartitionDict(x),
        .dictKeys(x, .CV_PARTITION_DICT),
        "cvPartition",
        "CV partition(s)"
    )
}

# The partition row `i` of a collection with mcols `md` shares, NULL for a
# row with none (or a collection saved before rows kept a key).
# @noRd
.mdCvPartition <- function(x, md, i) {
    keys <- md[["cvPartition"]]
    if (is.null(keys)) {
        return(NULL)
    }
    .dictEntryFor(keys[[i]], .cvPartitionDict(x))
}

# Every row's shared partition, NULL for a row with none.
# @noRd
.rowCvPartitions <- function(x) {
    map(
        .dictKeys(x, .CV_PARTITION_DICT),
        .dictEntryFor,
        entries = .cvPartitionDict(x)
    )
}

# The CV partitions of `parts`, row-bound: each row keeps the one it had.
# @noRd
.cvPartitionCombine <- function(parts) {
    perRow <- list_flatten(map(parts, .rowCvPartitions))
    .dictBuild(perRow, .CV_PARTITION_DICT$prefix)
}

# A collection's cvResult column for `cvs`, one row's CV each: the rows as
# they store them (.cvSplitPartition()), every row's key into the partitions
# they share, and those partitions once each.
# @noRd
.cvStoredColumn <- function(cvs) {
    split <- map(cvs, .cvSplitPartition)
    dict <- .dictBuild(map(split, "partition"), .CV_PARTITION_DICT$prefix)
    list(
        cvResult = S4Vectors::SimpleList(map(split, "stored")),
        keys = dict$keys,
        entries = dict$entries
    )
}

# One row's CV as a collection stores it, and the partition it shares: a
# merged row's per-region payloads each split, and a payload with no
# partition stays as it is.
# @noRd
.cvSplitPartition <- function(cv) {
    if (!.cvIsPerRegion(cv)) {
        return(.cvSplitPayload(cv))
    }
    split <- map(cv, .cvSplitPayload)
    shared <- map(split, "partition")
    list(
        stored = map(split, "stored"),
        partition = if (all(map_lgl(shared, is.null))) NULL else shared
    )
}

# One CV payload without its partition and with its predictions unnamed, and
# what it shares: the partition and the samples its predictions were held
# for, in their order.
# @noRd
.cvSplitPayload <- function(cv) {
    sp <- if (is.list(cv)) cv[["samplePartition"]]
    if (!is.data.frame(sp)) {
        return(list(stored = cv, partition = NULL))
    }
    ids <- .cvPredictionIds(cv)
    list(
        stored = list_assign(
            .cvUnnamePredictions(cv, ids),
            samplePartition = zap()
        ),
        partition = list(samplePartition = sp, sampleIds = ids)
    )
}

# The samples a payload's out-of-fold predictions are held for, in their
# order: a fine-mapping payload names them as the rows of its `prediction`
# matrices, a TWAS payload as the elements of its `predictions` vector. A
# TWAS matrix of predictions (one column per condition) is stored as it is,
# its rows named.
# @noRd
.cvPredictionIds <- function(cv) {
    if (length(cv[["prediction"]]) > 0L) {
        return(rownames(cv[["prediction"]][[1L]]))
    }
    names(cv[["predictions"]])
}

# @noRd
.cvUnnamePredictions <- function(cv, ids) {
    if (length(cv[["prediction"]]) > 0L) {
        unnamed <- map(cv[["prediction"]], .cvUnnameRows, ids = ids)
        return(list_assign(cv, prediction = unnamed))
    }
    p <- cv[["predictions"]]
    if (is.null(p) || !identical(names(p), ids)) {
        return(cv)
    }
    list_assign(cv, predictions = unname(p))
}

# @noRd
.cvUnnameRows <- function(m, ids) {
    if (!is.matrix(m) || !identical(rownames(m), ids)) {
        return(m)
    }
    `rownames<-`(m, NULL)
}

# A row's stored CV with the partition it shares put back: the partition
# first, as every CV payload holds it, and the predictions named again.
# @noRd
.cvJoinPartition <- function(stored, partition) {
    if (is.null(partition)) {
        return(stored)
    }
    if (is_in("samplePartition", names(partition))) {
        return(.cvJoinPayload(stored, partition))
    }
    map2(stored, partition, .cvJoinPayload)
}

# @noRd
.cvJoinPayload <- function(stored, partition) {
    if (is.null(partition)) {
        return(stored)
    }
    c(
        list(samplePartition = partition[["samplePartition"]]),
        .cvNamePredictions(stored, partition[["sampleIds"]])
    )
}

# @noRd
.cvNamePredictions <- function(cv, ids) {
    if (is.null(ids)) {
        return(cv)
    }
    if (length(cv[["prediction"]]) > 0L) {
        named <- map(cv[["prediction"]], .cvNameRows, ids = ids)
        return(list_assign(cv, prediction = named))
    }
    p <- cv[["predictions"]]
    if (
        is.null(p) ||
            is.matrix(p) ||
            !is.null(names(p)) ||
            length(p) != length(ids)
    ) {
        return(cv)
    }
    list_assign(cv, predictions = set_names(p, ids))
}

# @noRd
.cvNameRows <- function(m, ids) {
    if (!is.matrix(m) || !is.null(rownames(m)) || nrow(m) != length(ids)) {
        return(m)
    }
    `rownames<-`(m, ids)
}

# A collection saved while each row held its own CV partition and named
# predictions, with each partition kept once, as collections now store them.
# A collection whose rows already key into its partitions is left as is.
# @noRd
.cvSharePartitions <- function(object) {
    cols <- .tupleColumnNames(object)
    if (!is_in("cvResult", cols) || is_in("cvPartition", cols)) {
        return(object)
    }
    md <- mcols(object, use.names = FALSE)
    cv <- .cvStoredColumn(as.list(md$cvResult))
    md$cvResult <- cv$cvResult
    md$cvPartition <- cv$keys
    methods::slot(object, "cvPartition", check = FALSE) <- cv$entries
    S4Vectors::`mcols<-`(object, value = md)
}
