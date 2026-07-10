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

# Sample/Fold partition: shuffle samples, then cut into `fold` contiguous blocks.
# @noRd
.cvSamplePartition <- function(sampleNames, fold) {
  idx <- sample(length(sampleNames))
  folds <- cut(seq_along(sampleNames), breaks = fold, labels = FALSE)
  data.frame(Sample = sampleNames[idx], Fold = folds, stringsAsFactors = FALSE)
}

# Canonical CV output key: "<methodKey>_predicted" / "<methodKey>_performance".
# @noRd
.cvOutputKey <- function(methodKey, suffix) paste0(methodKey, "_", suffix)

# One CV metric row (corr, rsq, adj_rsq, pval, RMSE, MAE) for a single outcome.
# Drops NA in either vector; needs >= 3 valid points and non-constant predictions.
# @noRd
.cvMetricRow <- function(pred, actual) {
  out <- setNames(rep(NA_real_, 6L),
                  c("corr", "rsq", "adj_rsq", "pval", "RMSE", "MAE"))
  ok <- !is.na(pred) & !is.na(actual)
  pred <- pred[ok]; actual <- actual[ok]
  if (length(pred) < 3L || stats::sd(pred) == 0) return(out)
  lmFit <- stats::lm(actual ~ pred); s <- summary(lmFit)
  out["corr"]    <- stats::cor(actual, pred)
  out["rsq"]     <- s$r.squared
  out["adj_rsq"] <- s$adj.r.squared
  out["pval"]    <- if (nrow(s$coefficients) >= 2L) s$coefficients[2L, 4L] else NA_real_
  res <- actual - pred
  out["RMSE"] <- sqrt(mean(res^2))
  out["MAE"]  <- mean(abs(res))
  out
}

# Shared K-fold cross-validation engine.
#
# `fitFold(Xtrain, Ytrain, foldIndex)` must return
#   list(weights = <named list: methodKey -> (variants x outcomes) weight matrix,
#                   rownames indexing colnames(Xtrain)>,
#        fits    = <named list: methodKey -> fitted model or NULL>)
# The engine drops zero-variance training columns before calling fitFold,
# predicts held-out samples (Xtest[, common] %*% W[common, ]), aggregates the
# out-of-fold predictions into full-Y matrices, scores them with .cvMetricRow,
# and keys the output via .cvOutputKey(). `maxNumVariants` (optional) randomly
# subsamples variants up front to bound compute; `numThreads` parallelises the
# fold loop (-1 = all cores, 0/1 = serial).
#' @importFrom BiocParallel bplapply bpworkers MulticoreParam
#' @importFrom stats sd lm cor
#' @noRd
.crossValidateWeights <- function(X, Y, fold = NULL, samplePartitions = NULL,
                                  fitFold, numThreads = 1, maxNumVariants = NULL,
                                  variantsToKeep = NULL, retainFits = FALSE,
                                  verbose = 1) {
  if (!is.null(fold) && (!is.numeric(fold) || fold <= 0)) {
    stop("Invalid value for 'fold'. It must be a positive integer.")
  }
  if (!is.matrix(X) || (!is.matrix(Y) && !is.vector(Y))) {
    stop("X must be a matrix and Y must be a matrix or a vector.")
  }
  if (is.vector(Y)) {
    Y <- matrix(Y, ncol = 1)
    if (verbose >= 1) message(sprintf(
      "Y converted to matrix of %d rows and %d columns.", nrow(Y), ncol(Y)))
  }
  if (nrow(X) != nrow(Y)) {
    stop("The number of rows in X and Y must be the same.")
  }
  sampleNames <- if (!is.null(rownames(Y))) rownames(Y)
                 else if (!is.null(rownames(X))) rownames(X)
                 else paste0("sample_", seq_len(nrow(X)))
  if (is.null(rownames(X))) rownames(X) <- sampleNames
  if (is.null(rownames(Y))) rownames(Y) <- sampleNames
  if (is.null(colnames(X))) colnames(X) <- paste0("variable_", seq_len(ncol(X)))
  if (is.null(colnames(Y))) colnames(Y) <- paste0("context_", seq_len(ncol(Y)))

  # Optional variant subsample (compute saver; e.g. TWAS weight CV).
  if (!is.null(maxNumVariants) && ncol(X) > maxNumVariants) {
    if (!is.null(variantsToKeep) && length(variantsToKeep) > 0) {
      variantsToKeep <- intersect(variantsToKeep, colnames(X))
      remaining <- setdiff(colnames(X), variantsToKeep)
      if (length(variantsToKeep) < maxNumVariants) {
        additional <- sample(remaining, maxNumVariants - length(variantsToKeep))
        selected <- union(variantsToKeep, additional)
        if (verbose >= 1) message(sprintf(
          "Including %d specified variants and randomly selecting %d additional variants, for a total of %d variants out of %d for cross-validation purpose.",
          length(variantsToKeep), length(additional), length(selected), ncol(X)))
      } else {
        selected <- sample(variantsToKeep, maxNumVariants)
        if (verbose >= 1) message(sprintf(
          "Randomly selecting %d out of %d input variants for cross validation purpose.",
          length(selected), length(variantsToKeep)))
      }
    } else {
      selected <- sort(sample(ncol(X), maxNumVariants))
      if (verbose >= 1) message(sprintf(
        "Randomly selecting %d out of %d variants for cross validation purpose.",
        length(selected), ncol(X)))
    }
    X <- X[, selected, drop = FALSE]
  }

  # Fold partition: reuse a provided one, else build it.
  if (!is.null(samplePartitions)) {
    if (!all(samplePartitions$Sample %in% sampleNames)) {
      stop("Some samples in 'samplePartitions' do not match 'X' and 'Y'.")
    }
    if (!is.null(fold) && verbose >= 1 &&
        fold != length(unique(samplePartitions$Fold))) {
      message(sprintf(paste0(
        "fold number provided does not match with sample partition, performing ",
        "%d fold cross validation based on provided sample partition. "),
        length(unique(samplePartitions$Fold))))
    }
    samplePartition <- samplePartitions
  } else if (!is.null(fold)) {
    samplePartition <- .cvSamplePartition(sampleNames, fold)
  } else {
    stop("Either 'fold' or 'samplePartitions' must be provided.")
  }
  foldIds <- sort(unique(samplePartition$Fold))

  st <- proc.time()
  runFold <- function(j) {
    if (verbose >= 1) message(sprintf("  CV fold %s/%s ...", j, length(foldIds)))
    testIds <- samplePartition$Sample[samplePartition$Fold == j]
    isTest  <- rownames(X) %in% testIds
    if (all(isTest) || !any(isTest)) return(list(preds = list(), fits = list()))
    Xtr <- X[!isTest, , drop = FALSE]
    Xte <- X[isTest, , drop = FALSE]
    Ytr <- Y[!isTest, , drop = FALSE]
    keep <- .nonzeroVarColumns(Xtr)
    Xtr <- Xtr[, keep, drop = FALSE]
    ff <- fitFold(Xtr, Ytr, j)
    preds <- lapply(ff$weights, function(W) {
      if (is.null(W)) return(NULL)
      W[is.na(W)] <- 0
      common <- intersect(colnames(Xte), rownames(W))
      if (length(common) == 0L) return(NULL)
      yhat <- Xte[, common, drop = FALSE] %*% W[common, , drop = FALSE]
      rownames(yhat) <- rownames(Xte)
      yhat
    })
    list(preds = preds, fits = if (isTRUE(retainFits)) ff$fits else list())
  }

  numCores <- if (numThreads == -1) bpworkers(MulticoreParam()) else numThreads
  numCores <- min(numCores, bpworkers(MulticoreParam()))
  foldResults <- if (numCores >= 2) {
    bplapply(foldIds, runFold,
             BPPARAM = MulticoreParam(workers = numCores, RNGseed = 1L))
  } else {
    lapply(foldIds, runFold)
  }

  metricNames <- c("corr", "rsq", "adj_rsq", "pval", "RMSE", "MAE")
  methodKeys <- unique(unlist(lapply(foldResults, function(fr) names(fr$preds))))
  prediction <- list()
  performance <- list()
  for (mk in methodKeys) {
    predMat <- matrix(NA_real_, nrow(Y), ncol(Y), dimnames = dimnames(Y))
    for (fr in foldResults) {
      yh <- fr$preds[[mk]]
      if (!is.null(yh)) predMat[rownames(yh), ] <- yh
    }
    prediction[[.cvOutputKey(mk, "predicted")]] <- predMat
    perf <- t(vapply(seq_len(ncol(Y)), function(r) {
      if (verbose >= 1) {
        s <- suppressWarnings(stats::sd(predMat[, r], na.rm = TRUE))
        if (is.na(s) || s == 0) message(sprintf(paste0(
          "Predicted values for condition %d using %s have zero variance. ",
          "Filling performance metric with NAs"), r, mk))
      }
      .cvMetricRow(predMat[, r], Y[, r])
    }, setNames(numeric(6L), metricNames)))
    rownames(perf) <- colnames(Y)
    performance[[.cvOutputKey(mk, "performance")]] <- perf
  }

  foldFits <- setNames(lapply(foldResults, function(fr) {
    ff <- fr$fits
    ff[!vapply(ff, is.null, logical(1))]
  }), paste0("fold_", foldIds))

  list(samplePartition = samplePartition, prediction = prediction,
       performance = performance, foldFits = foldFits,
       timeElapsed = proc.time() - st)
}
