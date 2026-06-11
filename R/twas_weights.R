# Evaluate an expression while suppressing external package output.
# Catches both message() output (susieR, qgg) and Rprintf/cat stdout (mr.ash.alpha).
# @param expr An expression to evaluate.
# @return The result of evaluating expr.
# @noRd
.quietEval <- function(expr) {
  invisible(capture.output(
    result <- suppressMessages(expr),
    type = "output"
  ))
  result
}

# Rename a "_weights"/"Weights" suffix to the case-matching equivalent of `target`.
# Snake-case inputs get the underscored snake-case form (e.g. "lasso_weights" -> "lasso_predicted")
# and camelCase inputs get the CamelCase form (e.g. "lassoWeights" -> "lassoPredicted").
# Names without a recognized suffix are returned unchanged.
# @param x Character vector of names ending in "_weights" or "Weights".
# @param target A bare token such as "predicted" or "performance".
# @return Character vector with suffixes rewritten.
# @noRd
.renameSuffix <- function(x, target) {
  cap <- paste0(toupper(substr(target, 1, 1)), substr(target, 2, nchar(target)))
  x <- sub("_weights$", paste0("_", target), x)
  x <- sub("Weights$", cap, x)
  x
}

# Map short method names and presets to weightMethods lists.
# @param methods A character vector of short method names, or a preset string
#   ("default" or "fast_default").
# @return A named list suitable for the weightMethods parameter.
# @noRd
.twasMethodLookup <- function(methods) {
  # `fn` is the snake_case key used in weight method lists; `impl` is the
  # actual camelCase function name implemented by the package.
  methodMap <- list(
    susie = list(fn = "susie_weights", impl = "susieWeights", args = list(refine = FALSE, L = 20, L_greedy = 5)),
    susie_ash = list(fn = "susie_ash_weights", impl = "susieAshWeights", args = list()),
    susie_inf = list(fn = "susie_inf_weights", impl = "susieInfWeights", args = list()),
    mrash = list(fn = "mrash_weights", impl = "mrashWeights", args = list(initPriorSd = TRUE, max.iter = 100)),
    enet = list(fn = "enet_weights", impl = "enetWeights", args = list()),
    lasso = list(fn = "lasso_weights", impl = "lassoWeights", args = list()),
    bayes_r = list(fn = "bayes_r_weights", impl = "bayesRWeights", args = list()),
    bayes_l = list(fn = "bayes_l_weights", impl = "bLassoWeights", args = list()),
    bayes_a = list(fn = "bayes_a_weights", impl = "bayesAWeights", args = list()),
    bayes_b = list(fn = "bayes_b_weights", impl = "bayesBWeights", args = list()),
    bayes_c = list(fn = "bayes_c_weights", impl = "bayesCWeights", args = list()),
    bayes_n = list(fn = "bayes_n_weights", impl = "bayesNWeights", args = list()),
    b_lasso = list(fn = "b_lasso_weights", impl = "bLassoWeights", args = list()),
    dpr_vb = list(fn = "dpr_vb_weights", impl = "dprVbWeights", args = list()),
    dpr_gibbs = list(fn = "dpr_gibbs_weights", impl = "dprGibbsWeights", args = list()),
    dpr_adaptive_gibbs = list(fn = "dpr_adaptive_gibbs_weights", impl = "dprAdaptiveGibbsWeights", args = list()),
    scad = list(fn = "scad_weights", impl = "scadWeights", args = list()),
    mcp = list(fn = "mcp_weights", impl = "mcpWeights", args = list()),
    l0learn = list(fn = "l0learn_weights", impl = "l0learnWeights", args = list()),
    mvsusie = list(fn = "mvsusie_weights", impl = "mvsusieWeights", args = list(L = 30, L_greedy = 5)),
    mrmash = list(fn = "mrmash_weights", impl = "mrmashWeights", args = list())
  )

  # Handle presets
  fastDefault <- c("susie", "susie_inf", "mrash", "enet", "lasso", "mcp", "scad", "l0learn")
  if (length(methods) == 1) {
    if (methods == "fast_default") {
      methods <- fastDefault
    } else if (methods == "default") {
      methods <- c(fastDefault, "bayes_r", "bayes_c")
    }
  }

  # Build reverse map: function name -> short name, so full names are accepted too
  fnToShort <- setNames(
    names(methodMap),
    vapply(methodMap, function(x) x$fn, character(1))
  )
  # Normalize any full function names to short names
  methods <- vapply(methods, function(m) {
    if (m %in% names(fnToShort)) fnToShort[[m]] else m
  }, character(1), USE.NAMES = FALSE)

  unknown <- setdiff(methods, names(methodMap))
  if (length(unknown) > 0) {
    stop(
      "Unknown TWAS method(s): ", paste(unknown, collapse = ", "),
      ". Available methods: ", paste(names(methodMap), collapse = ", ")
    )
  }

  result <- list()
  for (m in methods) {
    entry <- methodMap[[m]]
    args <- entry$args
    # Track the actual function implementation name so downstream dispatchers
    # can resolve snake_case keys to the camelCase implementation.
    attr(args, "impl") <- entry$impl
    result[[entry$fn]] <- args
  }
  result
}

# Resolve the actual function name for a method key. Honors an "impl" attribute
# on the per-method args list (set by .twasMethodLookup), and otherwise applies
# a snake_case -> camelCase transformation as a fallback for user-supplied
# weightMethods lists.
.resolveMethodFunction <- function(methodKey, methodArgs = NULL) {
  # Search pecotmr's namespace explicitly so this works equally well when the
  # function is called either from inside the package or from a user session.
  ns <- asNamespace("pecotmr")
  fnExists <- function(name) {
    exists(name, mode = "function") ||
      exists(name, mode = "function", envir = ns, inherits = FALSE)
  }
  impl <- if (!is.null(methodArgs)) attr(methodArgs, "impl") else NULL
  if (!is.null(impl) && nzchar(impl) && fnExists(impl)) {
    return(impl)
  }
  # Direct match (e.g. caller already passed camelCase)
  if (fnExists(methodKey)) return(methodKey)
  # snake_case_weights -> camelCaseWeights
  parts <- strsplit(methodKey, "_", fixed = TRUE)[[1]]
  capRest <- paste0(toupper(substring(parts[-1], 1, 1)),
                    substring(parts[-1], 2))
  candidate <- paste0(parts[1], paste0(capRest, collapse = ""))
  if (fnExists(candidate)) return(candidate)
  methodKey
}

# Identify non-zero-variance columns of X. Returns a logical vector.
#' @importFrom matrixStats colSds
#' @noRd
.nonzeroVarColumns <- function(X) {
  sds <- colSds(X, na.rm = TRUE)
  !is.na(sds) & sds != 0
}

# Embed a smaller weights matrix into a full-sized zero matrix matching X and Y dimensions.
# @param weightsMatrix The fitted weights (nrow = number of valid columns).
# @param validColumns Logical or character vector identifying which columns of X were used.
# @param XColnames Column names of the original X.
# @param YColnames Column names of Y.
# @noRd
.embedWeights <- function(weightsMatrix, validColumns, nColsX, nColsY,
                          XColnames = NULL, YColnames = NULL) {
  full <- matrix(0, nrow = nColsX, ncol = nColsY)
  if (!is.null(XColnames)) rownames(full) <- XColnames
  if (!is.null(YColnames)) colnames(full) <- YColnames
  full[validColumns, ] <- weightsMatrix
  full
}

# Filter weight methods that produced all-zero weights from CV.
# Returns filtered weightMethods list and warns about removed methods.
# @noRd
.filterZeroWeightMethods <- function(weightMethods, twasWeightsRes) {
  wl <- if (is(twasWeightsRes, "TwasWeights")) getWeights(twasWeightsRes) else twasWeightsRes
  isAllZero <- vapply(wl, function(w) all(w == 0, na.rm = TRUE), logical(1))
  removed <- names(weightMethods)[isAllZero]
  if (length(removed) > 0) {
    warning(sprintf(
      "Methods %s are removed from CV because all their weights are zeros.",
      paste(removed, collapse = ", ")
    ))
  }
  weightMethods[!isAllZero]
}

.susieWeightIntermediate <- function(fit, X) {
  keep <- intersect(c("mu", "lbf_variable", "X_column_scale_factors", "pip", "theta"), names(fit))
  intermediate <- fit[keep]
  if (!is.null(fit$sets$cs)) {
    intermediate$cs_variants <- setNames(lapply(fit$sets$cs, function(L) colnames(X)[L]), names(fit$sets$cs))
    intermediate$cs_purity <- fit$sets$purity
  }
  intermediate
}

.prepareSusieWeightMethods <- function(X, Y, weightMethods, fittedModels = NULL) {
  if (is.vector(Y)) Y <- matrix(Y, ncol = 1)
  if (is.null(fittedModels)) fittedModels <- list()
  hasSusie <- !is.null(weightMethods[["susie_weights"]])
  hasSusieInf <- !is.null(weightMethods[["susie_inf_weights"]])
  susieFit <- if (hasSusie) weightMethods[["susie_weights"]][["susieFit"]] else NULL
  susieInfFit <- if (hasSusieInf) weightMethods[["susie_inf_weights"]][["susieInfFit"]] else NULL
  if (is.null(susieFit)) susieFit <- fittedModels[["susie"]]
  if (is.null(susieInfFit)) susieInfFit <- fittedModels[["susie_inf"]]

  if (!is.null(susieFit)) {
    susieFit <- .setFinemappingFitClass(susieFit, "susie")
  }
  if (!is.null(susieInfFit)) {
    susieInfFit <- .setFinemappingFitClass(susieInfFit, "susie_inf")
  }

  if (hasSusie && hasSusieInf && ncol(Y) == 1 &&
      is.null(susieFit) && is.null(susieInfFit)) {
    fitArgNames <- c("susieFit", "susieInfFit", "retainFit")
    fits <- fitSusieInfThenSusie(
      X,
      Y[, 1],
      args = weightMethods[["susie_weights"]][setdiff(names(weightMethods[["susie_weights"]]), fitArgNames)],
      susieInfArgs = modifyList(
        list(convergence_method = "pip"),
        weightMethods[["susie_inf_weights"]][setdiff(names(weightMethods[["susie_inf_weights"]]), fitArgNames)]
      ),
      fittedModels = list(susie = susieFit, susie_inf = susieInfFit)
    )
    susieFit <- fits[["susie"]]
    susieInfFit <- fits[["susie_inf"]]
  }

  if (!is.null(susieInfFit) && hasSusieInf) {
    weightMethods[["susie_inf_weights"]][["susieInfFit"]] <- susieInfFit
  }
  if (!is.null(susieFit) && hasSusie) {
    weightMethods[["susie_weights"]][["susieFit"]] <- susieFit
  }
  if (hasSusie &&
      is.null(weightMethods[["susie_weights"]][["susieFit"]]) &&
      !is.null(susieInfFit)) {
    weightMethods[["susie_weights"]] <- prepareSusieFromInfArgs(weightMethods[["susie_weights"]], susieInfFit)
  }
  weightMethods
}

#' Cross-Validation for weights selection in Transcriptome-Wide Association Studies (TWAS)
#'
#' Performs cross-validation for TWAS, supporting both univariate and multivariate methods.
#' It can either create folds for cross-validation or use pre-defined sample partitions.
#' For multivariate methods, it applies the method to the entire Y matrix for each fold.
#'
#' @param X A matrix of samples by features, where each row represents a sample and each column a feature.
#' @param Y A matrix (or vector, which will be converted to a matrix) of samples by outcomes, where each row corresponds to a sample.
#' @param fold An optional integer specifying the number of folds for cross-validation.
#' If NULL, 'samplePartitions' must be provided.
#' @param samplePartitions An optional dataframe with predefined sample partitions,
#' containing columns 'Sample' (sample names) and 'Fold' (fold number). If NULL, 'fold' must be provided.
#' @param weightMethods A list of methods and their specific arguments, formatted as list(method1 = method1_args, method2 = method2_args), or alternatively a character vector of method names (eg, c("susie_weights", "enet_weights")) in which case default arguments will be used for all methods.
#' methods in the list can be either univariate (applied to each column of Y) or multivariate (applied to the entire Y matrix).
#' @param maxNumVariants An optional integer to set the randomly selected maximum number of variants to use for CV purpose, to save computing time.
#' @param variantsToKeep An optional integer to ensure that the listed variants are kept in the CV when there is a limit on the maxNumVariants to use.
#' @param numThreads The number of threads to use for parallel processing.
#'        If set to -1, the function uses all available cores.
#'        If set to 0 or 1, no parallel processing is performed.
#'        If set to 2 or more, parallel processing is enabled with that many threads.
#' @param verbose Integer controlling verbosity level: 0 = suppress all messages,
#'   1 = suppress external package messages (default),
#'   2 = show all messages including those from external packages.
#' @return A list with the following components:
#' \itemize{
#'   \item `sample_partition`: A dataframe showing the sample partitioning used in the cross-validation.
#'   \item `prediction`: A list of matrices with predicted Y values for each method and fold.
#'   \item `metrics`: A matrix with rows representing methods and columns for various metrics:
#'     \itemize{
#'       \item `corr`: Pearson's correlation between predicated and observed values.
#'       \item `adj_rsq`: Adjusted R-squared value (which indicates the proportion of variance explained by the model) that accounts for the number of predictors in the model.
#'       \item `pval`: P-value assessing the significance of the model's predictions.
#'       \item `RMSE`: Root Mean Squared Error, a measure of the model's prediction error.
#'       \item `MAE`: Mean Absolute Error, a measure of the average magnitude of errors in a set of predictions.
#'     }
#'   \item `time_elapsed`: The time taken to complete the cross-validation process.
#' }
#' @importFrom purrr map
#' @importFrom BiocParallel bplapply bpworkers MulticoreParam
#' @importFrom quadprog solve.QP
#' @export
twasWeightsCv <- function(X, Y, fold = NULL, samplePartitions = NULL, weightMethods = NULL, maxNumVariants = NULL, variantsToKeep = NULL, numThreads = 1, verbose = 1, ...) {
  splitData <- function(X, Y, samplePartition, fold) {
    testIds <- samplePartition[which(samplePartition$Fold == fold), "Sample"]
    Xtrain <- X[!(rownames(X) %in% testIds), , drop = FALSE]
    Ytrain <- Y[!(rownames(Y) %in% testIds), , drop = FALSE]
    Xtest <- X[rownames(X) %in% testIds, , drop = FALSE]
    Ytest <- Y[rownames(Y) %in% testIds, , drop = FALSE]
    if (nrow(Xtrain) == 0 || nrow(Ytrain) == 0 || nrow(Xtest) == 0 || nrow(Ytest) == 0) {
      stop("Error: One of the datasets (train or test) has zero rows.")
    }
    return(list(Xtrain = Xtrain, Ytrain = Ytrain, Xtest = Xtest, Ytest = Ytest))
  }

  # Validation checks
  if (!is.null(fold) && (!is.numeric(fold) || fold <= 0)) {
    stop("Invalid value for 'fold'. It must be a positive integer.")
  }

  if (!is.matrix(X) || (!is.matrix(Y) && !is.vector(Y))) {
    stop("X must be a matrix and Y must be a matrix or a vector.")
  }

  if (is.vector(Y)) {
    Y <- matrix(Y, ncol = 1)
    if (verbose >= 1) message(paste("Y converted to matrix of", nrow(Y), "rows and", ncol(Y), "columns."))
  }

  if (nrow(X) != nrow(Y)) {
    stop("The number of rows in X and Y must be the same.")
  }
  if (!is.null(rownames(X)) && !is.null(rownames(Y))) {
    if (!identical(rownames(X), rownames(Y))) {
      rownames(X) <- rownames(Y)
    }
    sampleNames <- rownames(Y)
  } else if (!is.null(rownames(Y))) {
    sampleNames <- rownames(Y)
  } else if (!is.null(rownames(X))) {
    sampleNames <- rownames(X)
  } else {
    sampleNames <- paste0("sample_", 1:nrow(X))
  }
  if (is.null(rownames(X))) {
    rownames(X) <- sampleNames
  }
  if (is.null(rownames(Y))) {
    rownames(Y) <- sampleNames
  }

  if (is.null(colnames(X))) {
    colnames(X) <- paste0("variable_", 1:ncol(X))
  }
  if (is.null(colnames(Y))) {
    colnames(Y) <- paste0("context_", 1:ncol(Y))
  }

  if (is.character(weightMethods)) {
    weightMethods <- .twasMethodLookup(weightMethods)
  }

  if (!exists(".Random.seed")) {
    if (verbose >= 1) message("! No seed has been set. Please set seed for reproducable result. ")
  }

  # Select variants if necessary
  if (!is.null(maxNumVariants) && ncol(X) > maxNumVariants) {
    if (!is.null(variantsToKeep) && length(variantsToKeep) > 0) {
      variantsToKeep <- intersect(variantsToKeep, colnames(X))
      remainingColumns <- setdiff(colnames(X), variantsToKeep)
      if (length(variantsToKeep) < maxNumVariants) {
        additionalColumns <- sample(remainingColumns, maxNumVariants - length(variantsToKeep), replace = FALSE)
        selectedColumns <- union(variantsToKeep, additionalColumns)
        if (verbose >= 1) message(sprintf(
          "Including %d specified variants and randomly selecting %d additional variants, for a total of %d variants out of %d for cross-validation purpose.",
          length(variantsToKeep), length(additionalColumns), length(selectedColumns), ncol(X)
        ))
      } else {
        selectedColumns <- sample(variantsToKeep, maxNumVariants, replace = FALSE)
        if (verbose >= 1) message(paste("Randomly selecting", length(selectedColumns), "out of", length(variantsToKeep), "input variants for cross validation purpose."))
      }
    } else {
      selectedColumns <- sort(sample(ncol(X), maxNumVariants, replace = FALSE))
      if (verbose >= 1) message(paste("Randomly selecting", length(selectedColumns), "out of", ncol(X), "variants for cross validation purpose."))
    }
    X <- X[, selectedColumns, drop = FALSE]
  }

  # Create or use provided folds
  if (!is.null(fold)) {
    if (!is.null(samplePartitions)) {
      if (fold != length(unique(samplePartitions$Fold))) {
        if (verbose >= 1) message(paste0(
          "fold number provided does not match with sample partition, performing ", length(unique(samplePartitions$Fold)),
          " fold cross validation based on provided sample partition. "
        ))
      }

      folds <- samplePartitions$Fold
      samplePartition <- samplePartitions
    } else {
      sampleIndices <- sample(nrow(X))
      folds <- cut(seq(1, nrow(X)), breaks = fold, labels = FALSE)
      samplePartition <- data.frame(Sample = sampleNames[sampleIndices], Fold = folds, stringsAsFactors = FALSE)
    }
  } else if (!is.null(samplePartitions)) {
    if (!all(samplePartitions$Sample %in% sampleNames)) {
      stop("Some samples in 'samplePartitions' do not match the samples in 'X' and 'Y'.")
    }
    samplePartition <- samplePartitions
    fold <- length(unique(samplePartition$Fold))
  } else {
    stop("Either 'fold' or 'samplePartitions' must be provided.")
  }

  st <- proc.time()
  if (is.null(weightMethods)) {
    return(list(sample_partition = samplePartition))
  } else {
    # Hardcoded vector of multivariate weightMethods (accept both snake and camel)
    multivariateWeightMethods <- c("mrmash_weights", "mvsusie_weights",
                                    "mrmashWeights", "mvsusieWeights")

    # Determine the number of cores to use
    numCores <- ifelse(numThreads == -1,
      bpworkers(MulticoreParam()),
      numThreads)
    numCores <- min(numCores,
      bpworkers(MulticoreParam()))

    cvArgs <- list(...)

    # Perform CV with parallel processing
    computeMethodPredictions <- function(j) {
      if (verbose >= 1) {
        message(sprintf("  CV fold %d/%d ...", j, fold))
        tic()
      }
      datSplit <- splitData(X, Y, samplePartition = samplePartition, fold = j)
      Xtrain <- datSplit$Xtrain
      Ytrain <- datSplit$Ytrain
      Xtest <- datSplit$Xtest
      Ytest <- datSplit$Ytest

      # Remove columns with zero variance
      validColumns <- .nonzeroVarColumns(Xtrain)
      Xtrain <- Xtrain[, validColumns, drop = FALSE]
      Xtrain <- filterXWithY(Xtrain, Ytrain, missingRateThresh = 1, mafThresh = NULL)
      validColumns <- colnames(Xtrain)
      # Xtest <- Xtest[, validColumns, drop=FALSE]
      foldWeightMethods <- .prepareSusieWeightMethods(Xtrain, Ytrain, weightMethods)

      foldPreds <- setNames(lapply(names(foldWeightMethods), function(method) {
        args <- foldWeightMethods[[method]]
        fnName <- .resolveMethodFunction(method, args)

        if (method %in% multivariateWeightMethods) {
          # Apply multivariate method to entire Y for this fold
          if (!is.null(cvArgs$data_driven_prior_matrices_cv)) {
            if (method %in% c("mrmash_weights", "mrmashWeights")) {
              args$data_driven_prior_matrices <- cvArgs$data_driven_prior_matrices_cv[[j]]
            }
            if (method %in% c("mvsusie_weights", "mvsusieWeights")) {
              args$prior_variance <- cvArgs$reweighted_mixture_prior_cv[[j]]
            }
          }
          weightsMatrix <- if (verbose < 2) {
            .quietEval(do.call(fnName, c(list(X = Xtrain, Y = Ytrain), args)))
          } else {
            do.call(fnName, c(list(X = Xtrain, Y = Ytrain), args))
          }
          rownames(weightsMatrix) <- colnames(Xtrain)
          fullWeightsMatrix <- .embedWeights(weightsMatrix[validColumns, , drop = FALSE], validColumns, ncol(X), ncol(Y), colnames(X), colnames(Y))
          Ypred <- Xtest %*% fullWeightsMatrix
          rownames(Ypred) <- rownames(Xtest)
          return(Ypred)
        } else {
          Ypred <- sapply(1:ncol(Ytrain), function(k) {
            weights <- if (verbose < 2) {
              .quietEval(do.call(fnName, c(list(X = Xtrain, y = Ytrain[, k]), args)))
            } else {
              do.call(fnName, c(list(X = Xtrain, y = Ytrain[, k]), args))
            }
            fullWeights <- rep(0, ncol(X))
            names(fullWeights) <- colnames(X)
            fullWeights[validColumns] <- weights
            # Handle NAs in weights
            fullWeights[is.na(fullWeights)] <- 0
            Xtest %*% fullWeights
          })
          rownames(Ypred) <- rownames(Xtest)
          return(Ypred)
        }
      }), names(foldWeightMethods))
      if (verbose >= 1) {
        elapsed <- toc(quiet = TRUE)
        message(sprintf("  CV fold %d/%d done in %.1fs", j, fold, elapsed$toc - elapsed$tic))
      }
      foldPreds
    }

    if (numCores >= 2) {
      bpParam <- MulticoreParam(workers = numCores,
                                RNGseed = 1L)
      foldResults <- bplapply(1:fold,
        computeMethodPredictions, BPPARAM = bpParam)
    } else {
      foldResults <- map(1:fold, computeMethodPredictions)
    }

    # Reorganize into Ypred
    # After cross validation, each sample should have been in
    # test set at some point, and therefore has predicted value.
    # The prediction matrix is therefore exactly the same dimension as input Y
    Ypred <- setNames(lapply(weightMethods, function(x) `dimnames<-`(matrix(NA, nrow(Y), ncol(Y)), dimnames(Y))), names(weightMethods))
    for (j in seq_along(foldResults)) {
      for (method in names(weightMethods)) {
        Ypred[[method]][rownames(foldResults[[j]][[method]]), ] <- foldResults[[j]][[method]]
      }
    }

    names(Ypred) <- .renameSuffix(names(Ypred), "predicted")

    # Compute rsq, adj rsq, p-value, RMSE, and MAE for each method
    metricsTable <- list()

    for (m in names(weightMethods)) {
      metricsTable[[m]] <- matrix(NA, nrow = ncol(Y), ncol = 6)
      colnames(metricsTable[[m]]) <- c("corr", "rsq", "adj_rsq", "pval", "RMSE", "MAE")
      rownames(metricsTable[[m]]) <- colnames(Y)

      for (r in 1:ncol(Y)) {
        methodPredictions <- Ypred[[.renameSuffix(m, "predicted")]][, r]
        actualValues <- Y[, r]
        # Remove missing values in the first place
        naIndx <- which(is.na(actualValues))
        if (length(naIndx) != 0) {
          methodPredictions <- methodPredictions[-naIndx]
          actualValues <- actualValues[-naIndx]
        }
        if (sd(methodPredictions) != 0) {
          lmFit <- lm(actualValues ~ methodPredictions)

          # Calculate raw correlation and and adjusted R-squared
          metricsTable[[m]][r, "corr"] <- cor(actualValues, methodPredictions)

          metricsTable[[m]][r, "rsq"] <- summary(lmFit)$r.squared
          metricsTable[[m]][r, "adj_rsq"] <- summary(lmFit)$adj.r.squared

          # Calculate p-value
          metricsTable[[m]][r, "pval"] <- summary(lmFit)$coefficients[2, 4]

          # Calculate RMSE
          residuals <- actualValues - methodPredictions
          metricsTable[[m]][r, "RMSE"] <- sqrt(mean(residuals^2))

          # Calculate MAE
          metricsTable[[m]][r, "MAE"] <- mean(abs(residuals))
        } else {
          metricsTable[[m]][r, ] <- NA
          if (verbose >= 1) message(paste0(
            "Predicted values for condition ", r, " using ", m,
            " have zero variance. Filling performance metric with NAs"
          ))
        }
      }
    }
    names(metricsTable) <- .renameSuffix(names(metricsTable), "performance")
    return(list(sample_partition = samplePartition, prediction = Ypred, performance = metricsTable, time_elapsed = proc.time() - st))
  }
}

#' Run multiple TWAS weight methods
#'
#' Applies specified weight methods to the datasets X and Y, returning weight matrices for each method.
#' Handles both univariate and multivariate methods, and filters out columns in X with zero standard error.
#' This function utilizes parallel processing to handle multiple methods.
#'
#' @param X A matrix of samples by features, where each row represents a sample and each column a feature.
#' @param Y A matrix (or vector, which will be converted to a matrix) of samples by outcomes, where each row corresponds to a sample.
#' @param weightMethods A list of methods and their specific arguments, formatted as list(method1 = method1_args, method2 = method2_args), or alternatively a character vector of method names (eg, c("susie_weights", "enet_weights")) in which case default arguments will be used for all methods.
#' methods in the list can be either univariate (applied to each column of Y) or multivariate (applied to the entire Y matrix).
#' @param numThreads The number of threads to use for parallel processing.
#'        If set to -1, the function uses all available cores.
#'        If set to 0 or 1, no parallel processing is performed.
#'        If set to 2 or more, parallel processing is enabled with that many threads.
#' @param fittedModels Optional named list of fitted SuSiE-family models.
#' @param retainFits If TRUE, retain fitted model objects as attributes on
#'   returned weight matrices when supported by the weight method.
#' @param verbose Integer controlling verbosity level: 0 = suppress all messages,
#'   1 = suppress external package messages (default),
#'   2 = show all messages including those from external packages.
#' @return A list where each element is named after a method and contains the weight matrix produced by that method.
#'
#' @export
#' @importFrom purrr map exec
#' @importFrom rlang !!!
#' @importFrom tictoc tic toc
twasWeights <- function(X, Y, weightMethods, numThreads = 1,
                        fittedModels = NULL, retainFits = FALSE, verbose = 1) {
  if (!is.matrix(X) || (!is.matrix(Y) && !is.vector(Y))) {
    stop("X must be a matrix and Y must be a matrix or a vector.")
  }

  if (is.vector(Y)) {
    Y <- matrix(Y, ncol = 1)
  }

  if (nrow(X) != nrow(Y)) {
    stop("The number of rows in X and Y must be the same.")
  }

  if (is.character(weightMethods)) {
    weightMethods <- .twasMethodLookup(weightMethods)
  }

  # Determine number of cores to use
  numCores <- ifelse(numThreads == -1,
    bpworkers(MulticoreParam()),
    numThreads)
  numCores <- min(numCores,
    bpworkers(MulticoreParam()))

  validColumns <- .nonzeroVarColumns(X)
  Xfiltered <- as.matrix(X[, validColumns, drop = FALSE])
  weightMethods <- .prepareSusieWeightMethods(
    Xfiltered, Y, weightMethods, fittedModels
  )

  computeMethodWeights <- function(methodName, weightMethods) {
    shortName <- sub("_weights$", "", methodName)
    if (verbose >= 1) {
      message(sprintf("  Fitting %s ...", shortName))
      tic()
    }

    # Hardcoded vector of multivariate methods (accept both snake and camel)
    multivariateWeightMethods <- c("mrmash_weights", "mvsusie_weights",
                                    "mrmashWeights", "mvsusieWeights")
    args <- weightMethods[[methodName]]
    fnName <- .resolveMethodFunction(methodName, args)

    # Only pass retainFit (or its legacy snake_case alias) to functions that accept it
    if (retainFits) {
      fnFormals <- names(formals(fnName))
      if ("retainFit" %in% fnFormals) {
        args$retainFit <- TRUE
      } else if ("retain_fit" %in% fnFormals) {
        args$retain_fit <- TRUE
      }
    }

    methodFit <- NULL
    if (methodName %in% multivariateWeightMethods) {
      # Apply multivariate method
      weightsMatrix <- if (verbose < 2) {
        .quietEval(do.call(fnName, c(list(X = Xfiltered, Y = Y), args)))
      } else {
        do.call(fnName, c(list(X = Xfiltered, Y = Y), args))
      }
      if (retainFits) methodFit <- attr(weightsMatrix, "fit")
      if (nrow(weightsMatrix) != length(validColumns)) weightsMatrix <- weightsMatrix[names(validColumns), , drop = FALSE]
    } else {
      # Apply univariate method to each column of Y
      # Initialize it with zeros to avoid NA
      weightsMatrix <- matrix(0, nrow = ncol(Xfiltered), ncol = ncol(Y))

      for (k in 1:ncol(Y)) {
        weightsVector <- if (verbose < 2) {
          .quietEval(do.call(fnName, c(list(X = Xfiltered, y = Y[, k]), args)))
        } else {
          do.call(fnName, c(list(X = Xfiltered, y = Y[, k]), args))
        }
        if (retainFits && is.null(methodFit)) {
          methodFit <- attr(weightsVector, "fit")
        }
        if (is.matrix(weightsVector)) weightsVector <- weightsVector[, k]
        weightsMatrix[, k] <- weightsVector
      }
    }

    result <- .embedWeights(weightsMatrix, validColumns, ncol(X), ncol(Y), colnames(X), colnames(Y))
    if (!is.null(methodFit)) attr(result, "fit") <- methodFit
    if (verbose >= 1) {
      elapsed <- toc(quiet = TRUE)
      message(sprintf("  Fitting %s done in %.1fs", shortName, elapsed$toc - elapsed$tic))
    }
    return(result)
  }

  if (numCores >= 2) {
    bpParam <- MulticoreParam(workers = numCores,
                              RNGseed = 1L)
    weightsList <- bplapply(names(weightMethods),
      computeMethodWeights, weightMethods, BPPARAM = bpParam)
  } else {
    weightsList <- names(weightMethods) %>% map(computeMethodWeights, weightMethods)
  }
  names(weightsList) <- names(weightMethods)

  if (!is.null(colnames(X))) {
    weightsList <- lapply(weightsList, function(x) {
      fit <- attr(x, "fit")
      rownames(x) <- colnames(X)
      if (!is.null(fit)) attr(x, "fit") <- fit
      return(x)
    })
  }
  # Create TwasWeights S4 object
  variantIds <- if (!is.null(colnames(X))) colnames(X) else paste0("variant_", seq_len(ncol(X)))
  fitsList <- lapply(weightsList, function(w) attr(w, "fit"))
  hasAnyFit <- any(!sapply(fitsList, is.null))

  # Strip fit attributes from weight matrices before storing in S4
  cleanWeights <- lapply(weightsList, function(w) { attr(w, "fit") <- NULL; w })

  TwasWeights(
    weights = cleanWeights,
    variantIds = variantIds,
    fits = if (hasAnyFit) fitsList else NULL,
    cvPerformance = NULL
  )
}

#' Predict outcomes using TWAS weights
#'
#' This function takes a matrix of predictors (\code{X}) and a list of TWAS (transcriptome-wide
#' association studies) weights (\code{weightsList}), and calculates the predicted outcomes by
#' multiplying \code{X} by each set of weights in \code{weightsList}. The names of the elements
#' in the output list are derived from the names in \code{weightsList}, with "_weights" replaced
#' by "_predicted".
#'
#' @param X A matrix or data frame of predictors where each row is an observation and each
#' column is a variable.
#' @param weightsList A list of numeric vectors representing the weights for each predictor.
#' The names of the list elements should follow the pattern \code{[outcome]_weights}, where
#' \code{[outcome]} is the name of the outcome variable that the weights are associated with.
#'
#' @return A named list of numeric vectors, where each vector is the predicted outcome for the
#' corresponding set of weights in \code{weightsList}. The names of the list elements are
#' derived from the names in \code{weightsList} by replacing "_weights" with "_predicted".
#'
#' @export
#' @examples
#' # Assuming `X` is your matrix of predictors and `weightsList` is your list of weights:
#' predicted_outcomes <- twasPredict(X, weightsList)
#' print(predicted_outcomes)
twasPredict <- function(X, weightsList) {
  if (is(weightsList, "TwasWeights")) {
    wl <- getWeights(weightsList)
  } else {
    wl <- weightsList
  }
  setNames(lapply(wl, function(w) X %*% w), .renameSuffix(names(wl), "predicted"))
}

#' Estimate Sparsity from mr.ash Mixture Proportions
#'
#' Computes an empirical estimate of the proportion of non-zero effects
#' (sparsity) from the mr.ash fit. mr.ash fits a mixture model with a
#' point mass at zero (spike) plus continuous components (slab), and
#' learns the mixture proportions via variational EM. The sparsity
#' estimate \code{1 - pi[1]} is the empirical Bayes estimate of the
#' non-null proportion, which can be used as a data-driven prior for
#' the inclusion probability parameters (\code{pi} for bayesC,
#' \code{probIn} for BayesB) of spike-and-slab Bayesian methods.
#'
#' @param weightResults Named list of weight vectors or matrices as
#'   returned by \code{\link{twasWeights}}. The mr.ash element should
#'   have a \code{"fit"} attribute containing the model fit object
#'   (set \code{retainFits = TRUE} in \code{twasWeights} to obtain this).
#'
#' @return A scalar sparsity estimate (proportion of non-zero effects).
#' @export
estimateSparsity <- function(weightResults) {
  if (is(weightResults, "TwasWeights")) {
    fit <- getFits(weightResults, "mrash_weights")
    if (is.null(fit) || is.null(fit$pi)) {
      stop("mr.ash fit object not found. Run twasWeights() with retainFits = TRUE ",
           "and ensure mrash_weights is included.")
    }
  } else {
    w <- weightResults[["mrash_weights"]]
    if (is.null(w)) {
      stop("mr.ash weights ('mrash_weights') not found in weightResults.")
    }
    fit <- attr(w, "fit")
    if (is.null(fit) || is.null(fit$pi)) {
      stop("mr.ash fit object not found. Run twasWeights() with retainFits = TRUE ",
           "and ensure mrash_weights is included.")
    }
  }

  # fit$pi[1] is the weight on the spike (sa2[1] = 0); 1 - pi[1] = non-null proportion
  return(1 - fit$pi[1])
}

#' TWAS Weights Pipeline
#'
#' This function performs weights computation for Transcriptome-Wide Association Study (TWAS)
#' incorporating various steps such as filtering variants by linkage disequilibrium reference panel variants,
#' fitting models using SuSiE and other methods, and calculating TWAS weights and predictions.
#' Optionally, it can perform cross-validation for TWAS weights.
#'
#' @param X A matrix of genotype data where rows represent samples and columns represent genetic variants.
#' @param y A vector of phenotype measurements for each sample.
#' @param susieFit An object returned by the SuSiE function, containing the SuSiE model fit.
#' @param fittedModels Optional named list of fitted fine-mapping models, such
#'   as \code{list(susie = susieFit, susie_inf = susieInfFit)}.
#' @param cvFolds The number of folds to use for cross-validation. Set to 0 to skip cross-validation. Defaults to 5.
#' @param samplePartition Optional data frame with Sample and Fold columns for cross-validation. If NULL, a random partition is generated.
#' @param weightMethods List of methods to use to compute weights for TWAS; along with their parameters.
#' @param maxCvVariants The maximum number of variants to be included in cross-validation. Defaults to -1 which means no limit.
#' @param cvThreads The number of threads to use for parallel computation in cross-validation. Defaults to 1.
#' @param cvWeightMethods List of methods to use for cross-validation. If NULL, uses the same methods as weightMethods.
#' @param ensemble Logical. If TRUE and cvFolds > 1, learn ensemble combination
#'   weights via stacked regression (SR-TWAS). Requires at least two individual
#'   methods to have been run and to pass the R-squared cutoff. Defaults to TRUE.
#' @param ensembleR2Threshold Minimum cross-validated R-squared for an individual method
#'   to be included in the ensemble. Methods below this threshold are excluded.
#'   Defaults to 0.01.
#' @param ensembleSolver Character string specifying the optimization backend
#'   for ensemble learning. One of \code{"quadprog"}, \code{"nnls"},
#'   \code{"lbfgsb"}, or \code{"glmnet"}. Passed to
#'   \code{\link{ensembleWeights}}. Defaults to \code{"quadprog"}.
#' @param ensembleAlpha Elastic net mixing parameter, used only when
#'   \code{ensembleSolver = "glmnet"}. Defaults to 1 (lasso).
#' @param estimatePi If TRUE, estimate spike-and-slab sparsity from mr.ash
#'   before running Bayesian alphabet methods that need inclusion probabilities.
#' @param verbose Integer controlling verbosity level: 0 = suppress all messages,
#'   1 = show pecotmr messages but suppress external package messages (default),
#'   2 = show all messages including those from external packages.
#'
#' @return A list containing results from the TWAS pipeline, including TWAS weights, predictions, and optionally cross-validation results.
#' @export
#'
#' @examples
#' # Example usage (assuming appropriate objects for X, y, and susieFit are available):
#' twas_results <- twasWeightsPipeline(X, y, susieFit)
twasWeightsPipeline <- function(X,
                                y,
                                susieFit = NULL,
                                fittedModels = NULL,
                                cvFolds = 5,
                                samplePartition = NULL,
                                weightMethods = "default",
                                maxCvVariants = -1,
                                cvThreads = 1,
                                cvWeightMethods = NULL,
                                ensemble = TRUE,
                                ensembleR2Threshold = 0.01,
                                ensembleSolver = "quadprog",
                                ensembleAlpha = 1,
                                estimatePi = TRUE,
                                verbose = 1) {
  if (is.character(weightMethods)) {
    weightMethods <- .twasMethodLookup(weightMethods)
  }
  if (is.null(fittedModels)) fittedModels <- list()
  if (!is.null(susieFit)) fittedModels[["susie"]] <- susieFit

  res <- list()
  st <- proc.time()
  if (verbose >= 1) {
    message("Performing TWAS weights computation for univariate analysis methods ...")
    tic()
  }

  if (!is.null(fittedModels[["susie"]]) && !is.null(weightMethods$susie_weights)) {
    res$susie_weights_intermediate <- .susieWeightIntermediate(fittedModels[["susie"]], X)
  }

  # Check if empirical pi estimation is needed for spike-and-slab methods
  bayesCneedsPi <- "bayes_c_weights" %in% names(weightMethods) &&
    !"pi" %in% names(weightMethods$bayes_c_weights)
  bayesBneedsPi <- "bayes_b_weights" %in% names(weightMethods) &&
    !"probIn" %in% names(weightMethods$bayes_b_weights)
  needsPiEstimation <- (bayesCneedsPi || bayesBneedsPi) && estimatePi

  if (needsPiEstimation) {
    # Run mr.ash first to estimate sparsity
    mrashMethods <- list(mrash_weights = weightMethods[["mrash_weights"]] %||% list())

    if (verbose >= 1) message("  Estimating sparsity from mr.ash ...")
    mrashWeights <- twasWeights(X, y, weightMethods = mrashMethods, retainFits = TRUE, verbose = verbose)

    empiricalPi <- estimateSparsity(mrashWeights)
    if (verbose >= 1) message(sprintf("  Empirical sparsity estimate: %.4f", empiricalPi))
    res$empirical_pi <- empiricalPi

    # Inject into spike-and-slab methods that need it
    if (bayesCneedsPi) weightMethods$bayes_c_weights$pi <- as.numeric(empiricalPi)
    if (bayesBneedsPi) weightMethods$bayes_b_weights$probIn <- as.numeric(empiricalPi)

    # Run remaining methods (those not already computed)
    remainingFnNames <- setdiff(names(weightMethods), "mrash_weights")

    if (length(remainingFnNames) > 0) {
      remainingMethods <- weightMethods[remainingFnNames]
      remainingTw <- twasWeights(
        X,
        y,
        weightMethods = remainingMethods,
        fittedModels = fittedModels,
        verbose = verbose
      )
      # Combine two TwasWeights objects
      combinedWeights <- c(getWeights(mrashWeights), getWeights(remainingTw))
      combinedFits <- c(getFits(mrashWeights), getFits(remainingTw))
      res$twas_weights <- TwasWeights(
        weights = combinedWeights,
        variantIds = getVariantIds(mrashWeights),
        fits = combinedFits
      )
    } else {
      res$twas_weights <- mrashWeights
    }

    # Remove mr.ash if it was not in the original weightMethods
    if (!"mrash_weights" %in% names(weightMethods)) {
      wList <- getWeights(res$twas_weights)
      fList <- getFits(res$twas_weights)
      wList[["mrash_weights"]] <- NULL
      if (!is.null(fList)) fList[["mrash_weights"]] <- NULL
      res$twas_weights <- TwasWeights(
        weights = wList,
        variantIds = getVariantIds(res$twas_weights),
        fits = if (length(fList) > 0) fList else NULL
      )
    }
  } else {
    # Run all methods at once
    res$twas_weights <- twasWeights(
      X,
      y,
      weightMethods = weightMethods,
      fittedModels = fittedModels,
      verbose = verbose
    )
  }
  if (verbose >= 1) {
    elapsed <- toc(quiet = TRUE)
    message(sprintf("TWAS weights fitting done in %.1fs", elapsed$toc - elapsed$tic))
  }
  res$twas_predictions <- twasPredict(X, res$twas_weights)

  if (cvFolds > 1) {
    # A few cutting corners to run CV faster at the disadvantage of SuSiE and mr.ash:
    # 1. reset SuSiE to not using refine or adaptive L but to use L from previous analysis
    # 2. at most 100 iterations for mr.ash allowed
    # 3. only use a subset of variants randomly selected to avoid bias
    if (!is.null(fittedModels[["susie_inf"]]) && !is.null(weightMethods$susie_inf_weights)) {
      weightMethods$susie_inf_weights$L <- length(fittedModels[["susie_inf"]]$V)
      weightMethods$susie_inf_weights$refine <- FALSE
    }
    if (!is.null(weightMethods$susie_weights)) {
      susieCvFit <- fittedModels[["susie"]]
      if (is.null(susieCvFit)) susieCvFit <- fittedModels[["susie_inf"]]
      if (!is.null(susieCvFit)) {
        weightMethods$susie_weights$L <- length(susieCvFit$V)
        weightMethods$susie_weights$refine <- FALSE
      }
    }
    if (is.null(cvWeightMethods)) {
      cvWeightMethods <- .filterZeroWeightMethods(weightMethods, res$twas_weights)
    }

    variantsForCv <- c()
    if (maxCvVariants <= 0) {
      maxCvVariants <- Inf
    }
    if (ncol(X) > maxCvVariants) {
      variantsForCv <- sample(colnames(X), maxCvVariants, replace = FALSE)
    }

    if (verbose >= 1) {
      message("Performing cross-validation to assess TWAS weights ...")
      tic()
    }
    res$twas_cv_result <- twasWeightsCv(
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

    # Ensemble learning: learn optimal method combination via stacked regression
    if (isTRUE(ensemble) && length(cvWeightMethods) <= 1) {
      if (verbose >= 1) message("Ensemble model skipped: only ", length(cvWeightMethods),
              " weight method provided (need >= 2 for ensemble learning).")
    }
    if (isTRUE(ensemble) && length(cvWeightMethods) > 1) {
      if (!is.null(res$twas_cv_result$performance)) {
        # Extract R-squared for each method from CV performance table
        methodRsq <- vapply(res$twas_cv_result$performance, function(perf) {
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
          filteredCv <- res$twas_cv_result
          predBaseNames <- sub("(_predicted|Predicted)$", "", names(filteredCv$prediction))
          filteredCv$prediction <- filteredCv$prediction[match(passingBase, predBaseNames)]

          # Subset twas_weights to passing methods.
          # Original weight keys may use either snake_case (lasso_weights) or
          # camelCase (lassoWeights); match either by stripping the suffix.
          if (is(res$twas_weights, "TwasWeights")) {
            wl <- getWeights(res$twas_weights)
          } else {
            wl <- res$twas_weights
          }
          wlBaseNames <- sub("(_weights|Weights)$", "", names(wl))
          filteredWeights <- wl[match(passingBase, wlBaseNames)]
          # Rename filteredWeights keys to canonical _weights form for
          # downstream ensembleWeights lookup.
          names(filteredWeights) <- paste0(passingBase, "_weights")

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

          # Add ensemble weights alongside individual method weights
          if (!is.null(ensResult$ensemble_twas_weights)) {
            ensWt <- ensResult$ensemble_twas_weights
            if (!is.matrix(ensWt)) ensWt <- matrix(ensWt, ncol = 1)
            # Rebuild TwasWeights S4 with ensemble method added
            tw <- res$twas_weights
            newWeights <- c(getWeights(tw), list(ensembleWeights = ensWt))
            res$twas_weights <- new("TwasWeights",
              weights = newWeights,
              variantIds = getVariantIds(tw),
              methods = c(getMethodNames(tw), "ensembleWeights"),
              fits = getFits(tw),
              cvPerformance = getCvPerformance(tw),
              standardized = getStandardized(tw)
            )
            res$twas_predictions$ensemble_predicted <- X %*% ensWt
          }
          res$ensemble <- ensResult
        }
      }
    }
  }
  res$total_time_elapsed <- proc.time() - st

  return(res)
}

#' TWAS Multivariate Weights Pipeline
#'
#' This function performs weights computation for Transcriptome-Wide Association Study (TWAS)
#' in a multivariate setting. It incorporates steps such as fitting models using mvSuSiE and mr.mash,
#' calculating TWAS weights and predictions, and optionally performing cross-validation for TWAS weights.
#'
#' @param X A matrix of genotype data where rows represent samples and columns represent genetic variants.
#' @param Y A matrix of phenotype measurements, where rows represent samples and columns represent conditions.
#' @param mnmFit An object containing the fitted multivariate models (e.g., mvSuSiE and mr.mash fits).
#' @param L Maximum number of components in mvSuSiE. If NULL, the number of
#'   components in the fitted mvSuSiE object is used.
#' @param Lgreedy Initial greedy number of components in mvSuSiE. Defaults to 5.
#' @param cvFolds The number of folds to use for cross-validation. Defaults to 5. Set to 0 to skip cross-validation.
#' @param samplePartition Optional data frame with Sample and Fold columns for cross-validation. If NULL, a random partition is generated.
#' @param dataDrivenPriorMatrices A list of data-driven covariance matrices for mr.mash weights. Defaults to NULL.
#' @param dataDrivenPriorMatricesCv A list of data-driven covariance matrices for mr.mash weights in cross-validation. Defaults to NULL.
#' @param canonicalPriorMatrices If TRUE, computes canonical covariance matrices for mr.mash. Defaults to FALSE.
#' @param mvsusieMaxIter The maximum number of iterations for mvSuSiE. Defaults to 200.
#' @param mrmashMaxIter The maximum number of iterations for mr.mash. Defaults to 5000.
#' @param maxCvVariants The maximum number of variants to be included in cross-validation. Defaults to -1 which means no limit.
#' @param cvThreads The number of threads to use for parallel computation in cross-validation. Defaults to 1.
#' @param verbose Integer controlling verbosity level: 0 = suppress all messages,
#'   1 = show pecotmr messages but suppress external package messages (default),
#'   2 = show all messages including those from external packages.
#'
#' @return A list containing results from the TWAS pipeline, including TWAS weights, predictions, and optionally cross-validation results.
#' @export
#' @examples
#' # Example usage (assuming appropriate objects for X, Y, and mnmFit are available):
#' twas_results <- twasMultivariateWeightsPipeline(X, Y, mnmFit)
twasMultivariateWeightsPipeline <- function(
    X,
    Y,
    mnmFit,
    L = NULL,
    Lgreedy = 5,
    cvFolds = 5,
    samplePartition = NULL,
    dataDrivenPriorMatrices = NULL,
    dataDrivenPriorMatricesCv = NULL,
    canonicalPriorMatrices = FALSE,
    mvsusieMaxIter = 200,
    mrmashMaxIter = 5000,
    maxCvVariants = -1,
    cvThreads = 1,
    verbose = 1) {
  copyTwasResults <- function(contextNames, variantNames, twasWeight, twasPredictions) {
    wl <- if (is(twasWeight, "TwasWeights")) getWeights(twasWeight) else twasWeight
    setNames(lapply(contextNames, function(ctx) {
      if (ctx %in% colnames(wl[[1]])) {
        list(
          twas_weights = lapply(wl, function(wgts) wgts[, ctx]),
          twas_predictions = lapply(twasPredictions, function(pred) pred[, ctx]),
          variant_names = variantNames
        )
      } else {
        NULL
      }
    }), contextNames)
  }

  copyTwasCvResults <- function(twasResult, twasCvResult) {
    for (i in names(twasResult)) {
      if (i %in% colnames(twasCvResult$prediction[[1]])) {
        twasResult[[i]]$twas_cv_result$sample_partition <- twasCvResult$sample_partition
        twasResult[[i]]$twas_cv_result$prediction <- lapply(
          twasCvResult$prediction,
          function(predicted) {
            as.matrix(predicted[, i], ncol = 1)
          }
        )
        twasResult[[i]]$twas_cv_result$performance <- lapply(
          twasCvResult$performance,
          function(perform) {
            t(as.matrix(perform[i, ], ncol = 1))
          }
        )
        twasResult[[i]]$twas_cv_result$time_elapsed <- twasCvResult$time_elapsed
      }
    }
    return(twasResult)
  }

  # TWAS weights and predictions
  weightMethods <- list(
    mrmash_weights = list(
      mrmash_fit = mnmFit$mrmash_fitted
    ),
    mvsusie_weights = list(
      mvsusie_fit = mnmFit$mvsusie_fitted
    )
  )
  st <- proc.time()
  if (verbose >= 1) {
    message("Extracting TWAS weights for multivariate analysis methods ...")
    tic()
  }
  # get TWAS weights
  twasWeightsRes <- twasWeights(X = X, Y = Y, weightMethods = weightMethods, verbose = verbose)
  if (verbose >= 1) {
    elapsed <- toc(quiet = TRUE)
    message(sprintf("Multivariate TWAS weights fitting done in %.1fs", elapsed$toc - elapsed$tic))
  }
  # get TWAS predictions for possible next steps such as computing correlations between predicted expression values
  twasPredictions <- twasPredict(X, twasWeightsRes)

  # copy TWAS results by condition
  res <- copyTwasResults(colnames(Y), mnmFit$variant_names, twasWeightsRes, twasPredictions)

  # Perform cross-validation if specified
  if (cvFolds > 1) {
    if (is.null(L)) L <- length(mnmFit$mvsusie_fitted$V)
    if (!is.null(Lgreedy)) Lgreedy <- min(Lgreedy, L)
    subVerbose <- verbose >= 2
    weightMethods <- list(
      mrmash_weights = list(
        data_driven_prior_matrices = dataDrivenPriorMatrices,
        canonical_prior_matrices = canonicalPriorMatrices,
        max_iter = mrmashMaxIter,
        verbose = subVerbose
      ),
      mvsusie_weights = list(
        prior_variance = mnmFit$reweighted_mixture_prior,
        residual_variance = mnmFit$mrmash_fitted$V,
        L = L,
        L_greedy = Lgreedy,
        max_iter = mvsusieMaxIter,
        verbose = subVerbose
      )
    )

    weightMethods <- .filterZeroWeightMethods(weightMethods, twasWeightsRes)

    variantsForCv <- c()
    if (maxCvVariants <= 0) maxCvVariants <- Inf
    if (ncol(X) > maxCvVariants) {
      variantsForCv <- sample(colnames(X), maxCvVariants, replace = FALSE)
    }
    if (verbose >= 1) {
      message("Performing cross-validation to assess TWAS weights ...")
      tic()
    }
    twasCvResult <- twasWeightsCv(
      X = X, Y = Y, fold = cvFolds,
      weightMethods = weightMethods,
      samplePartitions = samplePartition,
      numThreads = cvThreads,
      maxNumVariants = maxCvVariants,
      verbose = verbose,
      variantsToKeep = if (length(variantsForCv) > 0) variantsForCv else NULL,
      data_driven_prior_matrices_cv = dataDrivenPriorMatricesCv,
      reweighted_mixture_prior_cv = mnmFit$reweighted_mixture_prior_cv
    )
    if (verbose >= 1) {
      elapsed <- toc(quiet = TRUE)
      message(sprintf("Cross-validation done in %.1fs", elapsed$toc - elapsed$tic))
    }
    res <- copyTwasCvResults(res, twasCvResult)
  }
  totalTimeElapsed <- proc.time() - st
  for (i in seq_along(res)) {
    res[[i]]$total_time_elapsed <- totalTimeElapsed
  }
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
#'   \code{\link{twasWeights}}, with keys like \code{"susie_weights"}. Used to
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
#'   \item{method_coef}{Named numeric vector of combination coefficients
#'     (\eqn{\zeta_k}), non-negative and summing to 1. Names are method
#'     base names (e.g., \code{"susie"}, \code{"enet"}).}
#'   \item{ensemble_twas_weights}{Final combined weight vector
#'     \eqn{w = \sum_k \zeta_k w_k}, or NULL if \code{twasWeightList}
#'     is not provided. Returned as a vector for univariate Y, matrix otherwise.}
#'   \item{method_performance}{Named numeric vector of per-method R-squared
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
#' @seealso \code{\link{twasWeightsCv}}, \code{\link{twasWeights}},
#'   \code{\link{twasWeightsPipeline}}
#'
#' @examples
#' \dontrun{
#' # After running twasWeightsPipeline with CV:
#' res <- twasWeightsPipeline(X, y, cvFolds = 5, weightMethods = methods)
#'
#' ens <- ensembleWeights(
#'   cvResults = res$twas_cv_result,
#'   Y = y,
#'   twasWeightList = res$twas_weights
#' )
#' ens$method_coef           # combination weights, sum to 1
#'
#' # Multi-dataset ensemble (e.g., CUMC1 + MIT cell types):
#' ens_multi <- ensembleWeights(
#'   cvResults = list(res_cumc$twas_cv_result, res_mit$twas_cv_result),
#'   Y = list(y_cumc, y_mit),
#'   twasWeightList = list(res_cumc$twas_weights, res_mit$twas_weights)
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
    method_coef = zeta,
    ensemble_twas_weights = ensembleTwasWt,
    method_performance = methodRsq
  )
}

# =============================================================================
# Summary-statistics TWAS weight training pipeline
# =============================================================================

# Internal: RAISS-impute QTL z-scores for LD-panel variants missing from the
# QTL summary statistics. Used by twasWeightsSumstatPipeline() when
# imputeMissing = TRUE. Returns the (possibly widened) sumstats data frame
# with new rows for imputed variants. Imputed variants with R^2 below the
# threshold are dropped by RAISS's internal filter.
imputeMissingSumstatsForLd <- function(sumstats, ldMat, ldData,
                                       imputeOpts, verbose = 1) {
  ldIds <- rownames(ldMat)
  missingIds <- setdiff(ldIds, sumstats$variant_id)
  if (length(missingIds) == 0) return(sumstats)

  # Build ref_panel covering all LD variants
  if (is(ldData, "LdData")) {
    ldRefPanel <- getRefPanel(ldData)
  } else {
    ldRefPanel <- parseVariantId(ldIds)
    ldRefPanel$variant_id <- ldIds
  }
  refCols <- c("chrom", "pos", "variant_id", "A1", "A2")
  if (!all(refCols %in% colnames(ldRefPanel))) {
    warning("imputeMissingSumstatsForLd: LD ref_panel missing required columns; skipping imputation.")
    return(sumstats)
  }
  if (!all(refCols %in% colnames(sumstats)) || !"z" %in% colnames(sumstats)) {
    warning("imputeMissingSumstatsForLd: sumstats missing required columns; skipping imputation.")
    return(sumstats)
  }

  # RAISS requires inputs sorted by position (within each chromosome)
  refSorted <- ldRefPanel[order(ldRefPanel$chrom, ldRefPanel$pos), refCols, drop = FALSE]
  knownSorted <- sumstats[order(sumstats$chrom, sumstats$pos), c(refCols, "z"), drop = FALSE]
  # Translate snake_case imputeOpts keys to camelCase raiss() arguments.
  imputeOptsRenamed <- imputeOpts
  if ("R2_threshold" %in% names(imputeOptsRenamed)) {
    imputeOptsRenamed$r2Threshold <- imputeOptsRenamed$R2_threshold
    imputeOptsRenamed$R2_threshold <- NULL
  }
  if ("minimum_ld" %in% names(imputeOptsRenamed)) {
    imputeOptsRenamed$minimumLd <- imputeOptsRenamed$minimum_ld
    imputeOptsRenamed$minimum_ld <- NULL
  }
  raissArgs <- c(list(
    refPanel = refSorted,
    knownZscores = knownSorted,
    ldMatrix = ldMat,
    verbose = (verbose >= 2)
  ), imputeOptsRenamed)
  raissOut <- tryCatch(do.call(raiss, raissArgs),
                       error = function(e) {
                         warning(sprintf("RAISS missing-sumstat imputation failed: %s", e$message))
                         NULL
                       })
  if (is.null(raissOut) || is.null(raissOut$result_filter)) return(sumstats)

  newRows <- raissOut$result_filter[
    !raissOut$result_filter$variant_id %in% sumstats$variant_id, , drop = FALSE
  ]
  if (nrow(newRows) == 0) return(sumstats)

  added <- newRows[, c("variant_id", "chrom", "pos", "A1", "A2", "z"), drop = FALSE]
  if ("beta" %in% colnames(sumstats)) added$beta <- newRows$z
  if ("se"   %in% colnames(sumstats)) added$se   <- 1
  for (col in setdiff(colnames(sumstats), colnames(added))) {
    added[[col]] <- NA
  }
  added <- added[, colnames(sumstats), drop = FALSE]
  if (verbose >= 1) {
    message(sprintf("RAISS imputed %d missing QTL sumstat variants from LD reference.",
                    nrow(added)))
  }
  rbind(sumstats, added)
}

#' Train TWAS weights from summary statistics and LD reference
#'
#' Replaces the OTTERS pipeline with a properly integrated workflow that:
#' (1) runs RSS QC on eQTL summary statistics, (2) trains weights via multiple
#' RSS methods, and (3) extracts fine-mapping results from the shared SuSiE-RSS
#' fit. Returns a \code{TwasWeights} S4 object with \code{standardized = TRUE}
#' that feeds directly into \code{harmonize_twas} and \code{twas_analysis}.
#'
#' @param sumstats Data.frame with columns: \code{variant_id}, \code{A1},
#'   \code{A2}, \code{chrom}, \code{pos}, and either \code{z} or both
#'   \code{beta} and \code{se}.
#' @param ldData LdData S4 object, or a legacy list with \code{LD_matrix},
#'   \code{LD_variants}, \code{ref_panel}. Can also be a plain correlation
#'   matrix (variant IDs taken from row/colnames).
#' @param n eQTL study sample size (scalar).
#' @param methods Named list of RSS weight methods and their arguments.
#'   Method names correspond to functions named
#'   \code{<method>_weights(stat, LD, ...)}. Defaults include lassosum_rss,
#'   prs_cs, sdpr, susie_rss, and susie_inf_rss.
#' @param pThresholds Numeric vector of p-value thresholds for P+T weights.
#'   Set to NULL to skip.
#' @param checkLdMethod LD matrix repair method: \code{"eigenfix"} (default),
#'   \code{"shrink"}, or NULL to skip.
#' @param qcMethod RSS QC method for eQTL data: \code{"slalom"},
#'   \code{"dentist"}, or NULL/\code{"none"} to skip.
#' @param keepIndel Whether to keep indels during QC. Default TRUE.
#' @param pipCutoffToSkip PIP threshold for early stopping. Default 0 (off).
#' @param impute Whether to run RAISS imputation of LD-inconsistent variants
#'   flagged by QC (the QC re-imputation path). Default FALSE.
#' @param imputeMissing Logical. When \code{TRUE}, RAISS imputes QTL z-scores
#'   for variants present in the LD reference but absent from the QTL
#'   summary statistics, after QC and before LD/sumstats intersection. This
#'   widens the sumstats panel available to the weight-learning methods so a
#'   richer set of weights can later be applied to GWAS. Independent of
#'   \code{impute}; both can be enabled together. Default \code{FALSE}.
#' @param imputeOpts RAISS imputation parameters; shared by the \code{impute}
#'   QC re-imputation and the \code{imputeMissing} missing-variant path.
#'   Imputed variants with \code{R2 < R2_threshold} are dropped.
#' @param varY Phenotype variance. Default 1.
#' @param verbose Verbosity level.
#'
#' @return A list with:
#' \describe{
#'   \item{twas_weights}{A \code{TwasWeights} S4 object with
#'     \code{standardized = TRUE}.}
#'   \item{finemapping_result}{A \code{FineMappingResult} S4 object from the
#'     SuSiE-RSS fit, or NULL if no SuSiE-RSS method was used.}
#'   \item{qc_summary}{List with outlier counts and QC metadata.}
#' }
#'
#' @export
twasWeightsSumstatPipeline <- function(
    sumstats, ldData, n,
    methods = list(
      lassosum_rss = list(),
      prs_cs = list(phi = 1e-4, n_iter = 1000, n_burnin = 500, thin = 5),
      sdpr = list(iter = 1000, burn = 200, thin = 1, verbose = FALSE),
      susie_rss = list(),
      susie_inf_rss = list()
    ),
    pThresholds = c(0.001, 0.05),
    checkLdMethod = "eigenfix",
    qcMethod = NULL,
    keepIndel = TRUE,
    pipCutoffToSkip = 0,
    impute = TRUE,
    imputeMissing = FALSE,
    imputeOpts = list(rcond = 0.01, R2_threshold = 0.6,
                      minimum_ld = 5, lamb = 0.01),
    varY = 1, verbose = 1) {

  # -----------------------------------------------------------------------
  # 1. RSS QC on eQTL summary statistics
  # -----------------------------------------------------------------------
  needsQc <- !is.null(qcMethod) && !identical(qcMethod, "none")
  if (needsQc || impute || pipCutoffToSkip != 0) {
    qcResult <- summaryStatsQc(
      rssInput = list(sumstats = sumstats, n = n, var_y = varY),
      ldData = ldData,
      keepIndel = keepIndel,
      pipCutoffToSkip = pipCutoffToSkip,
      qcMethod = qcMethod,
      impute = impute,
      imputeOpts = imputeOpts,
      returnOnSkip = "null"
    )
    if (is.null(qcResult) || isSkipped(qcResult)) {
      return(list(twas_weights = NULL, finemapping_result = NULL,
                  qc_summary = list(skipped = TRUE)))
    }
    sumstats <- getRssInput(qcResult)$sumstats
    qcLd <- getLdData(qcResult)
    ldMat <- if (is.null(qcLd)) NULL else if (hasGenotypes(qcLd)) getGenotypes(qcLd) else getCorrelation(qcLd)
    outlierNumber <- getOutlierNumber(qcResult)
  } else {
    # No QC requested: extract LD matrix directly
    if (is.matrix(ldData)) {
      ldMat <- ldData
    } else if (is(ldData, "LdData")) {
      ldMat <- getCorrelation(ldData)
    } else {
      stop("ldData must be a matrix or LdData object.")
    }
    outlierNumber <- 0L
  }

  if (nrow(sumstats) < 2) {
    return(list(twas_weights = NULL, finemapping_result = NULL,
                qc_summary = list(skipped = TRUE, reason = "fewer than 2 variants")))
  }

  # -----------------------------------------------------------------------
  # 2. Compute z-scores and build stat object
  # -----------------------------------------------------------------------
  if (is.null(sumstats$z)) {
    if (!is.null(sumstats$beta) && !is.null(sumstats$se)) {
      sumstats$z <- sumstats$beta / sumstats$se
    } else {
      stop("sumstats must have 'z' or ('beta' and 'se') columns.")
    }
  }

  p <- nrow(sumstats)
  z <- sumstats$z
  variantIds <- sumstats$variant_id
  b <- z / sqrt(n)
  stat <- list(b = b, cor = b, z = z, n = rep(n, p))

  # Optional RAISS imputation: fill QTL z-scores for LD-panel variants absent
  # from the QTL summary statistics. Widens the sumstats panel so weight
  # learners have access to a richer variant set; downstream intersection
  # with LD becomes a near-identity after this step.
  if (isTRUE(imputeMissing) && !is.null(ldMat) && !is.null(rownames(ldMat))) {
    sumstats <- imputeMissingSumstatsForLd(
      sumstats = sumstats,
      ldMat = ldMat,
      ldData = ldData,
      imputeOpts = imputeOpts,
      verbose = verbose
    )
    p <- nrow(sumstats)
    z <- sumstats$z
    variantIds <- sumstats$variant_id
    b <- z / sqrt(n)
    stat <- list(b = b, cor = b, z = z, n = rep(n, p))
  }

  # Align LD matrix to sumstats variant order
  if (!is.null(rownames(ldMat)) && !is.null(variantIds)) {
    common <- intersect(variantIds, rownames(ldMat))
    if (length(common) < p) {
      idx <- match(common, variantIds)
      sumstats <- sumstats[idx, , drop = FALSE]
      z <- sumstats$z
      variantIds <- sumstats$variant_id
      b <- z / sqrt(n)
      stat <- list(b = b, cor = b, z = z, n = rep(n, length(z)))
      p <- length(z)
    }
    ldMat <- ldMat[variantIds, variantIds, drop = FALSE]
  }

  # -----------------------------------------------------------------------
  # 3. LD eigenfix (optional)
  # -----------------------------------------------------------------------
  if (!is.null(checkLdMethod)) {
    ldCheck <- check_ld(ldMat, method = checkLdMethod)
    if (ldCheck$method_applied != "none") {
      if (verbose >= 1) {
        message(sprintf("check_ld: repaired LD via '%s' (min eigenvalue was %.2e, %d negative).",
                        ldCheck$method_applied, ldCheck$min_eigenvalue, ldCheck$n_negative))
      }
    }
    ldMat <- ldCheck$R
  }

  # -----------------------------------------------------------------------
  # 4. Two-stage SuSiE-RSS (shared fit for susie_rss + susie_inf_rss)
  # -----------------------------------------------------------------------
  hasSusieRss <- "susie_rss" %in% names(methods)
  hasSusieInfRss <- "susie_inf_rss" %in% names(methods)
  susieFits <- NULL

  if (hasSusieRss && hasSusieInfRss) {
    susieArgs <- methods[["susie_rss"]]
    susieInfArgs <- methods[["susie_inf_rss"]]
    susieFits <- fitSusieInfThenSusieRss(
      z = z, R = ldMat, n = n,
      susieInfArgs = susieInfArgs,
      susieArgs = susieArgs
    )
  }

  # -----------------------------------------------------------------------
  # 5. P+T weights
  # -----------------------------------------------------------------------
  results <- list()
  if (!is.null(pThresholds)) {
    pvals <- pchisq(z^2, df = 1, lower.tail = FALSE)
    for (thr in pThresholds) {
      selected <- pvals < thr
      w <- ifelse(selected, stat$b, 0)
      results[[paste0("PT_", thr)]] <- w
    }
  }

  # -----------------------------------------------------------------------
  # 6. RSS method dispatch
  # -----------------------------------------------------------------------
  susieRssFitForFm <- NULL

  for (methodName in names(methods)) {
    fnName <- .resolveMethodFunction(paste0(methodName, "_weights"))
    if (!exists(fnName, mode = "function")) {
      warning(sprintf("Method '%s' not found (looking for function '%s'). Skipping.",
                      methodName, fnName))
      next
    }

    methodArgs <- methods[[methodName]]

    # Build call arguments: separate pre-fitted objects from methodArgs
    callArgs <- list(stat = stat, LD = ldMat)
    if (methodName == "susie_rss" && !is.null(susieFits)) {
      callArgs[["susie_rss_fit"]] <- susieFits$susie
    } else if (methodName == "susie_inf_rss" && !is.null(susieFits)) {
      callArgs[["susie_inf_rss_fit"]] <- susieFits$susie_inf
    }

    # SuSiE-RSS methods accept the user-facing method options as a single
    # `methodArgs` list; other methods take them spread directly.
    isSusieRssMethod <- methodName %in% c("susie_rss", "susie_inf_rss", "susie_ash_rss")
    if (isSusieRssMethod) {
      callArgs[["methodArgs"]] <- methodArgs
    } else {
      callArgs <- c(callArgs, methodArgs)
    }

    tryCatch({
      w <- do.call(fnName, callArgs)
      # Capture retained fit for fine-mapping post-processing
      if (methodName == "susie_rss" && !is.null(attr(w, "fit"))) {
        susieRssFitForFm <- attr(w, "fit")
      } else if (methodName == "susie_inf_rss" && is.null(susieRssFitForFm) && !is.null(attr(w, "fit"))) {
        susieRssFitForFm <- attr(w, "fit")
      }
      results[[methodName]] <- as.numeric(w)
    }, error = function(e) {
      warning(sprintf("Method '%s' failed: %s", methodName, e$message))
      results[[methodName]] <<- rep(0, p)
    })
  }

  if (length(results) == 0) {
    return(list(twas_weights = NULL, finemapping_result = NULL,
                qc_summary = list(skipped = TRUE, reason = "all methods failed")))
  }

  # -----------------------------------------------------------------------
  # 7. Fine-mapping from SuSiE-RSS fit (reuses the same fit)
  # -----------------------------------------------------------------------
  finemappingResult <- NULL
  if (!is.null(susieRssFitForFm)) {
    fmFits <- list(susie_rss = susieRssFitForFm)
    tryCatch({
      fmOutput <- postprocessFinemappingFits(
        fits = fmFits,
        dataX = ldMat,
        coverage = 0.95,
        signalCutoff = 0.025,
        csInput = "Xcorr"
      )
      if (!is.null(fmOutput$finemapping_results$susie_rss$finemapping_result)) {
        finemappingResult <- fmOutput$finemapping_results$susie_rss$finemapping_result
      }
    }, error = function(e) {
      warning(sprintf("Fine-mapping post-processing failed: %s", e$message))
    })
  }

  # -----------------------------------------------------------------------
  # 8. Package into TwasWeights S4
  # -----------------------------------------------------------------------
  weightsList <- lapply(results, function(w) {
    matrix(w, ncol = 1, dimnames = list(variantIds, NULL))
  })

  twasWt <- TwasWeights(
    weights = weightsList,
    variantIds = variantIds,
    standardized = TRUE,
    cvPerformance = NULL
  )

  list(
    twas_weights = twasWt,
    finemapping_result = finemappingResult,
    qc_summary = list(
      skipped = FALSE,
      n_variants_input = p,
      n_variants_after_qc = nrow(sumstats),
      outlier_number = outlierNumber,
      methods_succeeded = names(results)
    )
  )
}

# =============================================================================
# Multivariate summary-statistics TWAS weight training pipeline
# =============================================================================

#' Train multi-context TWAS weights from per-context summary statistics
#'
#' Multi-context summary-statistics analog of
#' \code{\link{twasMultivariateWeightsPipeline}}. Bundles per-context RSS
#' QC, cross-context variant alignment, optional RAISS missing-variant
#' imputation, data-driven prior construction (reusing
#' \code{\link{build_mrmash_prior_matrices}} and the same FLASH / diagonal
#' covariance helpers as the individual-level pipeline), and multi-context
#' weight training via \code{\link{mrmash_rss_weights}} and/or
#' \code{\link{mvsusie_rss_weights}}.
#'
#' @param sumstatsList Named list of per-context sumstats data.frames. Each
#'   data.frame must contain \code{variant_id}, \code{chrom}, \code{pos},
#'   \code{A1}, \code{A2}, and either \code{z} or \code{beta}/\code{se}.
#'   List names become condition labels.
#' @param ldData Shared \code{LdData} S4 object covering the union of
#'   variants. Each per-context sumstats is QC'd against this same LD panel.
#' @param n Per-context sample sizes; either a named numeric vector matching
#'   \code{names(sumstatsList)} or a single scalar to broadcast.
#' @param methods Named list of multivariate RSS weight methods to fit.
#'   Function names must match \code{<name>_weights}. Defaults to mr.mash-RSS
#'   and mvSuSiE-RSS with default arguments.
#' @param qcMethod Per-context QC method passed to
#'   \code{\link{summary_stats_qc}}; one of \code{"slalom"}, \code{"dentist"},
#'   or NULL/\code{"none"}. Default \code{NULL} (basic harmonization only).
#' @param keepIndel Passed through to QC. Default TRUE.
#' @param impute Logical. If TRUE, run per-context RAISS re-imputation of
#'   LD-mismatch outliers (QC re-imputation). Default FALSE.
#' @param imputeMissing Logical. If TRUE, after per-context QC and before
#'   cross-context alignment, RAISS imputes per-context z-scores for
#'   LD-reference variants absent from each context's sumstats. Widens the
#'   intersection used downstream. Default FALSE.
#' @param imputeOpts Named list of RAISS parameters shared by both
#'   \code{impute} and \code{imputeMissing}.
#' @param dataDrivenPriorMatrices Optional list of pre-computed prior
#'   matrices (with element \code{U}) passed to
#'   \code{\link{build_mrmash_prior_matrices}}. When NULL and
#'   \code{estimatePriorsFromSumstats = TRUE}, the pipeline estimates a
#'   data-driven covariance from the cross-context \code{Bhat} matrix via
#'   \code{\link{compute_cov_flash}}. If FLASH fails the error is allowed
#'   to propagate; supply \code{dataDrivenPriorMatrices} explicitly or
#'   set \code{estimatePriorsFromSumstats = FALSE} to bypass.
#' @param canonicalPriorMatrices Passed to
#'   \code{\link{build_mrmash_prior_matrices}}. Default TRUE.
#' @param estimatePriorsFromSumstats Logical. When TRUE (default) and
#'   \code{dataDrivenPriorMatrices} is NULL, estimate data-driven priors
#'   from the cross-context Bhat matrix using
#'   \code{\link{compute_cov_flash}}. FLASH errors are not swallowed; see
#'   \code{dataDrivenPriorMatrices} for the explicit-prior path.
#' @param verbose Integer verbosity level.
#'
#' @return A list with
#' \describe{
#'   \item{twas_weights}{A \code{TwasWeights} S4 object with per-context
#'     weight matrices (variants x conditions).}
#'   \item{qc_summary}{Per-context QC and alignment counts.}
#'   \item{Z}{The aligned z-score matrix (variants x conditions) fed to the
#'     weight learners.}
#' }
#' @export
twasMultivariateWeightsSumstatPipeline <- function(
    sumstatsList, ldData, n,
    methods = list(mrmash_rss = list(), mvsusie_rss = list()),
    qcMethod = NULL,
    keepIndel = TRUE,
    impute = FALSE,
    imputeMissing = FALSE,
    imputeOpts = list(rcond = 0.01, R2_threshold = 0.6,
                      minimum_ld = 5, lamb = 0.01),
    dataDrivenPriorMatrices = NULL,
    canonicalPriorMatrices = TRUE,
    estimatePriorsFromSumstats = TRUE,
    verbose = 1) {

  # ----- 1. Validate inputs and normalise per-context n -----
  if (!is.list(sumstatsList) || length(sumstatsList) == 0L) {
    stop("sumstatsList must be a non-empty named list of sumstats data.frames.")
  }
  if (is.null(names(sumstatsList)) || any(names(sumstatsList) == "")) {
    stop("sumstatsList must be a named list; names become condition labels.")
  }
  conditions <- names(sumstatsList)
  K <- length(conditions)
  if (length(n) == 1L) n <- setNames(rep(as.numeric(n), K), conditions)
  if (is.null(names(n))) names(n) <- conditions
  missingN <- setdiff(conditions, names(n))
  if (length(missingN) > 0) {
    stop("n vector missing entries for conditions: ", paste(missingN, collapse = ", "))
  }

  # ----- 2. Per-context QC + optional missing-variant imputation -----
  perContextQc <- list()
  for (cond in conditions) {
    ssC <- sumstatsList[[cond]]
    if (!is.null(ssC$z) || (!is.null(ssC$beta) && !is.null(ssC$se))) {
      if (is.null(ssC$z)) ssC$z <- ssC$beta / ssC$se
    } else {
      stop(sprintf("Context %s: sumstats must contain z or (beta, se).", cond))
    }
    qcRecord <- summary_stats_qc(
      rss_input = list(sumstats = ssC, n = n[[cond]], var_y = 1),
      LD_data = ldData,
      keep_indel = keepIndel,
      qc_method = if (is.null(qcMethod)) "none" else qcMethod,
      impute = impute, impute_opts = imputeOpts,
      return_on_skip = "preprocess",
      study = cond
    )
    ssQced <- getRssInput(qcRecord)$sumstats
    if (isTRUE(imputeMissing) && nrow(ssQced) > 0) {
      qcLd <- getLdData(qcRecord)
      ldForImpute <- if (is.null(qcLd)) {
        if (hasGenotypes(ldData)) getGenotypes(ldData) else getCorrelation(ldData)
      } else if (hasGenotypes(qcLd)) getGenotypes(qcLd) else getCorrelation(qcLd)
      if (!is.null(ldForImpute) && !is.null(rownames(ldForImpute))) {
        ssQced <- imputeMissingSumstatsForLd(
          sumstats = ssQced, ldMat = ldForImpute,
          ldData = ldData, imputeOpts = imputeOpts, verbose = verbose
        )
      }
    }
    perContextQc[[cond]] <- ssQced
  }

  # ----- 3. Cross-context alignment (intersection on variant_id) -----
  variantSets <- lapply(perContextQc, function(df) df$variant_id)
  commonVariants <- Reduce(intersect, variantSets)
  if (length(commonVariants) < 2) {
    return(list(twas_weights = NULL,
                qc_summary = list(skipped = TRUE,
                                  reason = "fewer than 2 shared variants across contexts"),
                Z = NULL))
  }

  # Use LD reference order for the common set
  ldIds <- if (is(ldData, "LdData")) getVariantIds(ldData) else rownames(getCorrelation(ldData))
  commonVariants <- intersect(ldIds, commonVariants)

  # ----- 4. Build Z and Bhat/Shat matrices -----
  Z <- matrix(NA_real_, nrow = length(commonVariants), ncol = K,
              dimnames = list(commonVariants, conditions))
  Bhat <- Z; Shat <- Z
  nVec <- as.numeric(n[conditions])
  for (k in seq_len(K)) {
    df <- perContextQc[[conditions[k]]]
    df <- df[df$variant_id %in% commonVariants, , drop = FALSE]
    idx <- match(commonVariants, df$variant_id)
    Z[, k] <- df$z[idx]
    Bhat[, k] <- Z[, k] / sqrt(nVec[k])
    Shat[, k] <- 1 / sqrt(nVec[k])
  }

  # ----- 5. LD subset to common variants -----
  ldFull <- if (is(ldData, "LdData")) getCorrelation(ldData) else ldData
  if (is.null(ldFull)) {
    if (is(ldData, "LdData") && hasGenotypes(ldData)) {
      Xref <- getGenotypes(ldData)
      ldFull <- compute_LD(Xref[, commonVariants, drop = FALSE], method = "sample")
    } else {
      stop("ldData must provide either a correlation matrix or a genotype handle.")
    }
  }
  ldMat <- ldFull[commonVariants, commonVariants, drop = FALSE]

  # ----- 6. Data-driven prior matrices from the aligned Bhat -----
  priorInput <- dataDrivenPriorMatrices
  if (is.null(priorInput) && isTRUE(estimatePriorsFromSumstats) && K >= 2) {
    if (verbose >= 1) message("Estimating data-driven prior matrices from Bhat ...")
    # Let FLASH errors propagate; callers can supply dataDrivenPriorMatrices
    # explicitly or set estimatePriorsFromSumstats = FALSE to bypass.
    priorInput <- list(U = list(flash = compute_cov_flash(Bhat)))
  }

  # ----- 7. Build stat object and dispatch weight methods -----
  stat <- list(z = Z, Bhat = Bhat, Shat = Shat, n = nVec)
  results <- list()
  for (methodName in names(methods)) {
    fnName <- .resolveMethodFunction(paste0(methodName, "_weights"))
    if (!exists(fnName, mode = "function")) {
      warning(sprintf("Method '%s' not found (looking for function '%s'). Skipping.",
                      methodName, fnName))
      next
    }
    methodArgs <- methods[[methodName]] %||% list()
    if (methodName == "mrmash_rss") {
      methodArgs$data_driven_prior_matrices <- methodArgs$data_driven_prior_matrices %||% priorInput
      methodArgs$canonical_prior_matrices  <- methodArgs$canonical_prior_matrices  %||% canonicalPriorMatrices
    }
    fn <- get(fnName, mode = "function")
    if (verbose >= 1) message("Fitting ", methodName, " ...")
    w <- tryCatch(
      do.call(fn, c(list(stat = stat, LD = ldMat), methodArgs)),
      error = function(e) {
        warning(sprintf("Method '%s' failed: %s", methodName, e$message))
        NULL
      }
    )
    if (!is.null(w)) {
      if (!is.matrix(w)) w <- matrix(w, nrow = length(commonVariants), ncol = K,
                                     dimnames = list(commonVariants, conditions))
      results[[paste0(methodName, "_weights")]] <- w
    }
  }

  # ----- 8. Package into TwasWeights S4 -----
  if (length(results) == 0) {
    return(list(twas_weights = NULL,
                qc_summary = list(skipped = TRUE,
                                  reason = "all methods failed"),
                Z = Z))
  }
  twasWt <- TwasWeights(
    weights = results,
    variantIds = commonVariants,
    standardized = TRUE,
    cvPerformance = NULL
  )
  list(
    twas_weights = twasWt,
    qc_summary = list(
      skipped = FALSE,
      n_per_context = vapply(perContextQc, nrow, integer(1)),
      n_common = length(commonVariants),
      conditions = conditions,
      methods_succeeded = names(results)
    ),
    Z = Z
  )
}
