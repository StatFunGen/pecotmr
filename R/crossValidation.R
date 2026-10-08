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
# subsamples variants up front to bound compute; `numThreads` parallelises the
# fold loop (-1 = all cores, 0/1 = serial).
#' @importFrom BiocParallel bplapply multicoreWorkers MulticoreParam
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
    numThreads = 1,
    maxNumVariants = NULL,
    variantsToKeep = NULL,
    fitRetention = c("none", "slim", "full"),
    verbose = 1,
    seed = NULL
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
    numCores <- .cvNumCores(numThreads)
    cvState <- list(
        X = X,
        Y = Y,
        samplePartition = samplePartition,
        foldIds = foldIds,
        fitFold = fitFold,
        fitFoldCtx = fitFoldCtx,
        fitRetention = fitRetention,
        verbose = verbose,
        rngSeed = seed
    )
    foldResults <- .cvRunFolds(foldIds, cvState, numCores)
    agg <- .cvAggregate(foldResults, Y, verbose)
    list(
        samplePartition = samplePartition,
        prediction = agg$prediction,
        performance = agg$performance,
        foldFits = .cvCollectFoldFits(foldResults, foldIds),
        timeElapsed = proc.time() - st
    )
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

# @noRd
# multicoreWorkers() rather than bpworkers(MulticoreParam()): the two answer
# the same number, but constructing a MulticoreParam costs ~0.6s (almost all
# of it garbage collection) and this runs on every fit.
.cvNumCores <- function(numThreads) {
    avail <- multicoreWorkers()
    min(if (numThreads == -1) avail else numThreads, avail)
}

# Run each fold (parallel via BiocParallel when >= 2 cores).
# @noRd
# MulticoreParam for CV / weight-fitting parallelism. A NULL seed keeps the
# historical default RNGseed = 1L (reproducible); a caller-supplied seed
# overrides it so the parallel L'Ecuyer streams are reproducible under it. This
# is the one place a user seed reaches BiocParallel workers, which set.seed()
# in the main process cannot.
# @noRd
.bpSeedParam <- function(numCores, seed = NULL) {
    MulticoreParam(workers = numCores, RNGseed = as.integer(seed %||% 1L))
}

# Seed the main-process RNG when a seed is supplied. This covers only the serial
# draws (fold partitioning, variant sub-sampling, single-threaded fitting);
# BiocParallel workers seed their own L'Ecuyer streams via .bpSeedParam. A NULL
# seed leaves the session RNG untouched so an outer set.seed() still governs it.
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

.cvRunFolds <- function(foldIds, cvState, numCores) {
    if (numCores >= 2) {
        return(bplapply(
            foldIds,
            .cvRunFold,
            cv = cvState,
            BPPARAM = .bpSeedParam(numCores, cvState$rngSeed)
        ))
    }
    map(foldIds, .cvRunFold, cv = cvState)
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

#' @rdname CrossValidationParam
#' @aliases CrossValidationParam-class
#' @exportClass CrossValidationParam
setClass(
    "CrossValidationParam",
    contains = "MethodParam",
    slots = c(
        folds = "numeric",
        numThreads = "numeric",
        samplePartition = "list_OR_NULL",
        maxVariants = "numeric_OR_NULL",
        weightMethods = "list_OR_NULL"
    )
)

#' @title Cross-Validation Settings
#' @description How cross-validation runs, shared by
#'   \code{\link{fineMappingPipeline}} and \code{\link{twasWeightsPipeline}}.
#'
#'   \code{folds}, \code{numThreads} and \code{samplePartition} mean the same
#'   thing in both. \code{maxVariants} and \code{weightMethods} are honoured
#'   by \code{\link{twasWeightsPipeline}} only and are \strong{ignored} by
#'   \code{\link{fineMappingPipeline}}, so one bundle can be handed to
#'   either pipeline unchanged.
#' @param folds Integer. Number of folds; \code{0} (the default) or
#'   \code{1} skips cross-validation. The same value in both pipelines ---
#'   \code{twasWeightsPipeline} used to default to \code{5}, which made an
#'   unspecified setting mean two different things.
#' @param numThreads Integer. Parallel workers for the per-fold refits.
#'   \code{1} (default) is serial; \code{-1} uses all cores. Only consulted
#'   when \code{folds > 1}.
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
#' @return A \code{CrossValidationParam} object, a \code{\link{MethodParam}}.
#' @examples
#' CrossValidationParam(folds = 10, numThreads = 4)
#' @export
CrossValidationParam <- function(
    folds = 0,
    numThreads = 1,
    samplePartition = NULL,
    maxVariants = NULL,
    weightMethods = NULL
) {
    new(
        "CrossValidationParam",
        folds = folds,
        numThreads = numThreads,
        samplePartition = samplePartition,
        maxVariants = maxVariants,
        weightMethods = weightMethods
    )
}

# The fields a pipeline reads, as a plain list so call sites index it without
# caring that the caller passed a constructor result.
# @noRd
.cvResolve <- function(crossValidationArgs) {
    list(
        folds = crossValidationArgs$folds %||% 0,
        numThreads = crossValidationArgs$numThreads %||% 1,
        samplePartition = crossValidationArgs$samplePartition,
        # -1 is twasWeightsCv()'s "no cap" sentinel. fineMappingPipeline
        # resolves these two as well but never reads them.
        maxVariants = crossValidationArgs$maxVariants %||% -1,
        weightMethods = crossValidationArgs$weightMethods
    )
}

# TRUE when cross-validation will actually run. Both 0 and 1 mean "no CV",
# and anything that needs out-of-fold predictions has to ask this rather than
# testing `folds` itself.
# @noRd
.cvEnabled <- function(crossValidationArgs) {
    folds <- crossValidationArgs$folds %||% 0
    !is.null(folds) && length(folds) == 1L && !is.na(folds) && folds >= 2L
}

# Refuse a cross-validation request on an input that has no samples to hold
# out. Summary statistics carry no individual-level data, so CV there is not
# unimplemented but meaningless -- saying so beats running to completion and
# returning results that were never cross-validated.
# @noRd
.cvRefuseOnSumstats <- function(crossValidationArgs, pipeline, cls) {
    set <- names(crossValidationArgs)[
        map_lgl(as.list(crossValidationArgs), .cvFieldIsSet)
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
# @noRd
.cvFieldIsSet <- function(x) {
    if (is.null(x)) {
        return(FALSE)
    }
    # numThreads = 1 and folds = 0 are the no-op defaults.
    !(length(x) == 1L && !is.na(x) && is.numeric(x) && (x == 0 || x == 1))
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
# default -- so `ensemble` defaults to FALSE too, and asking for an ensemble
# without folds is an error rather than a silently missing row.
#
# Two of the five fields are honoured by twasWeightsPipeline() only, and are
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
