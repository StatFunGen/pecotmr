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
    susieAsh = list(fn = "susie_ash_weights", impl = "susieAshWeights", args = list()),
    susieInf = list(fn = "susie_inf_weights", impl = "susieInfWeights", args = list()),
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
  fastDefault <- c("susie", "susieInf", "mrash", "enet", "lasso", "mcp", "scad", "l0learn")
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
  if (is(twasWeightsRes, "TwasWeights")) {
    methodTokens <- as.character(twasWeightsRes$method)
    perMethodAllZero <- vapply(seq_len(nrow(twasWeightsRes)), function(i) {
      w <- getWeights(twasWeightsRes$entry[[i]])
      all(w == 0, na.rm = TRUE)
    }, logical(1))
    methodToZero <- tapply(perMethodAllZero, methodTokens, all)
    methodKeys <- names(weightMethods)
    methodBase <- sub("(_weights|Weights)$", "", methodKeys)
    isAllZero <- vapply(methodBase, function(mb) {
      if (mb %in% names(methodToZero)) isTRUE(methodToZero[[mb]]) else FALSE
    }, logical(1))
  } else {
    wl <- twasWeightsRes
    isAllZero <- vapply(wl, function(w) all(w == 0, na.rm = TRUE), logical(1))
  }
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
    intermediate$csVariants <- setNames(lapply(fit$sets$cs, function(L) colnames(X)[L]), names(fit$sets$cs))
    intermediate$csPurity <- .translateSusiePurity(fit$sets$purity)
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
  if (is.null(susieInfFit)) susieInfFit <- fittedModels[["susieInf"]]

  if (!is.null(susieFit)) {
    susieFit <- .setFinemappingFitClass(susieFit, "susie")
  }
  if (!is.null(susieInfFit)) {
    susieInfFit <- .setFinemappingFitClass(susieInfFit, "susieInf")
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
      fittedModels = list(susie = susieFit, susieInf = susieInfFit)
    )
    susieFit <- fits[["susie"]]
    susieInfFit <- fits[["susieInf"]]
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
#'   \item `samplePartition`: A dataframe showing the sample partitioning used in the cross-validation.
#'   \item `prediction`: A list of matrices with predicted Y values for each method and fold.
#'   \item `metrics`: A matrix with rows representing methods and columns for various metrics:
#'     \itemize{
#'       \item `corr`: Pearson's correlation between predicated and observed values.
#'       \item `adj_rsq`: Adjusted R-squared value (which indicates the proportion of variance explained by the model) that accounts for the number of predictors in the model.
#'       \item `pval`: P-value assessing the significance of the model's predictions.
#'       \item `RMSE`: Root Mean Squared Error, a measure of the model's prediction error.
#'       \item `MAE`: Mean Absolute Error, a measure of the average magnitude of errors in a set of predictions.
#'     }
#'   \item `timeElapsed`: The time taken to complete the cross-validation process.
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
    return(list(samplePartition = samplePartition))
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
              args$prior_variance <- cvArgs$reweightedMixturePriorCv[[j]]
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
    return(list(samplePartition = samplePartition, prediction = Ypred, performance = metricsTable, timeElapsed = proc.time() - st))
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
learnTwasWeights <- function(X, Y, weightMethods,
                             study = "", context = "", trait = "",
                             numThreads = 1,
                             fittedModels = NULL,
                             retainFits = FALSE,
                             standardized = FALSE,
                             dataType = NULL,
                             ldSketch = NULL,
                             verbose = 1) {
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

  variantIds <- if (!is.null(colnames(X))) colnames(X) else paste0("variant_", seq_len(ncol(X)))
  traitLabels <- if (!is.null(colnames(Y))) colnames(Y) else paste0("outcome_", seq_len(ncol(Y)))

  # Build one TwasWeightsEntry per (method, trait/outcome) row. For multi-
  # outcome (multivariate) methods the per-method weights matrix has one
  # column per outcome, so the same row in the TwasWeights collection
  # carries the matrix for that method across outcomes via a single
  # `trait` value taken from the input `trait` arg (when length 1) or the
  # corresponding Y column name when `trait` matches `colnames(Y)`.
  buildEntries <- function() {
    studies   <- character(0)
    contexts  <- character(0)
    traits    <- character(0)
    methodsV  <- character(0)
    entries   <- list()
    for (m in names(weightsList)) {
      wMat <- weightsList[[m]]
      fitVal <- attr(wMat, "fit")
      attr(wMat, "fit") <- NULL
      shortMethod <- sub("(_weights|Weights)$", "", m)
      # When trait/context were supplied per-row (length == ncol(Y)), emit
      # one row per (method, outcome). Otherwise emit one row per method
      # and carry the (possibly multi-column) weights matrix as-is.
      perOutcome <- length(trait) == ncol(Y) &&
                    length(context) %in% c(1L, ncol(Y))
      if (perOutcome) {
        contextV <- if (length(context) == 1L) rep(context, ncol(Y)) else context
        studyV   <- if (length(study) == 1L) rep(study, ncol(Y)) else study
        for (k in seq_len(ncol(Y))) {
          studies  <- c(studies,  studyV[k])
          contexts <- c(contexts, contextV[k])
          traits   <- c(traits,   trait[k])
          methodsV <- c(methodsV, shortMethod)
          entries[[length(entries) + 1L]] <- TwasWeightsEntry(
            variantIds    = variantIds,
            weights       = wMat[, k],
            fits          = if (retainFits) fitVal else NULL,
            cvPerformance = NULL,
            standardized  = isTRUE(standardized),
            dataType      = dataType)
        }
      } else {
        studies  <- c(studies,  study[1L])
        contexts <- c(contexts, context[1L])
        traits   <- c(traits,   trait[1L])
        methodsV <- c(methodsV, shortMethod)
        wPayload <- if (ncol(wMat) == 1L) drop(wMat) else wMat
        entries[[length(entries) + 1L]] <- TwasWeightsEntry(
          variantIds    = variantIds,
          weights       = wPayload,
          fits          = if (retainFits) fitVal else NULL,
          cvPerformance = NULL,
          standardized  = isTRUE(standardized),
          dataType      = dataType)
      }
    }
    list(study = studies, context = contexts, trait = traits,
         method = methodsV, entry = entries)
  }

  rows <- buildEntries()
  TwasWeights(
    study   = rows$study,
    context = rows$context,
    trait   = rows$trait,
    method  = rows$method,
    entry   = rows$entry,
    ldSketch = ldSketch)
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
    # Per-row weights vector/matrix payloads. Use the method name as key
    # for compatibility with the legacy snake_case "<method>_predicted"
    # convention; ensembleWeights() rebinds the suffix.
    methodNames <- as.character(weightsList$method)
    wl <- setNames(
      lapply(seq_len(nrow(weightsList)), function(i) {
        getWeights(weightsList$entry[[i]])
      }),
      paste0(methodNames, "_weights"))
  } else {
    wl <- weightsList
  }
  setNames(lapply(wl, function(w) {
    if (!is.matrix(w)) w <- matrix(w, ncol = 1)
    X %*% w
  }), .renameSuffix(names(wl), "predicted"))
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
#'   returned by \code{\link{learnTwasWeights}}. The mr.ash element should
#'   have a \code{"fit"} attribute containing the model fit object
#'   (set \code{retainFits = TRUE} in \code{learnTwasWeights} to obtain this).
#'
#' @return A scalar sparsity estimate (proportion of non-zero effects).
#' @export
estimateSparsity <- function(weightResults) {
  if (is(weightResults, "TwasWeights")) {
    # Method names on the new TwasWeights collection are bare tokens
    # ("mrash"), not the snake_case _weights suffix form.
    methods <- as.character(weightResults$method)
    idx <- which(methods == "mrash")
    if (length(idx) == 0L) {
      stop("mr.ash entry not found in TwasWeights. Run learnTwasWeights() ",
           "with retainFits = TRUE and ensure 'mrash' is in the method list.")
    }
    fit <- getFits(weightResults$entry[[idx[[1L]]]])
    if (is.null(fit) || is.null(fit$pi)) {
      stop("mr.ash fit object not found. Run learnTwasWeights() with retainFits = TRUE ",
           "and ensure mrash_weights is included.")
    }
  } else {
    w <- weightResults[["mrash_weights"]]
    if (is.null(w)) {
      stop("mr.ash weights ('mrash_weights') not found in weightResults.")
    }
    fit <- attr(w, "fit")
    if (is.null(fit) || is.null(fit$pi)) {
      stop("mr.ash fit object not found. Run learnTwasWeights() with retainFits = TRUE ",
           "and ensure mrash_weights is included.")
    }
  }

  # fit$pi[1] is the weight on the spike (sa2[1] = 0); 1 - pi[1] = non-null proportion
  return(1 - fit$pi[1])
}

# =============================================================================
# Helpers + S4 dispatch surface for twasWeightsPipeline
# =============================================================================

# Concatenate two TwasWeights collections row-wise. `rbind` on DFrame
# subclasses does not reliably preserve the `ldSketch` slot, so this
# helper rebuilds via the constructor.
# @noRd
.rbindTwasWeights <- function(a, b, ldSketch = NULL) {
  if (!is(a, "TwasWeights") || !is(b, "TwasWeights")) {
    stop(".rbindTwasWeights expects two TwasWeights inputs.")
  }
  TwasWeights(
    study   = c(as.character(a$study),   as.character(b$study)),
    context = c(as.character(a$context), as.character(b$context)),
    trait   = c(as.character(a$trait),   as.character(b$trait)),
    method  = c(as.character(a$method),  as.character(b$method)),
    entry   = c(as.list(a$entry), as.list(b$entry)),
    ldSketch = ldSketch)
}

# Splice per-(method, outcome) cross-validated predictions and the 6-metric
# performance row from a `twasWeightsCv()` result into the matching
# `TwasWeightsEntry$cvPerformance` slot of every row in a TwasWeights
# collection. Rebuilds the collection because TwasWeightsEntry is treated
# as immutable. Rows for which no CV result is available (method not in
# the CV run, or trait not in the CV prediction matrix's columns) are
# emitted unchanged.
#
# The CV result keys carry a method suffix (`<m>_predicted`,
# `<m>_performance` in snake form, or `<m>Predicted`, `<m>Performance` in
# camel form); the TwasWeights `method` column carries the bare token
# (e.g. "lasso"). The trait column carries the outcome name, which must
# match the column name of the CV prediction matrix.
# @noRd
.spliceCvIntoTwasWeights <- function(twasWeights, twasCvResult,
                                      ldSketch = NULL) {
  if (is.null(twasCvResult) || is.null(twasCvResult$prediction) ||
      is.null(twasCvResult$performance)) {
    return(twasWeights)
  }
  predKeyBase <- sub("(_predicted|Predicted)$", "",
                     names(twasCvResult$prediction))
  perfKeyBase <- sub("(_performance|Performance)$", "",
                     names(twasCvResult$performance))

  pickKey <- function(bare, keys, base) {
    hit <- which(base == bare)
    if (length(hit) == 0L) NA_character_ else keys[[hit[[1L]]]]
  }

  studies   <- as.character(twasWeights$study)
  contexts  <- as.character(twasWeights$context)
  traits    <- as.character(twasWeights$trait)
  methodsV  <- as.character(twasWeights$method)
  newEntries <- as.list(twasWeights$entry)

  for (i in seq_along(newEntries)) {
    bare <- methodsV[[i]]
    pKey <- pickKey(bare, names(twasCvResult$prediction), predKeyBase)
    mKey <- pickKey(bare, names(twasCvResult$performance), perfKeyBase)
    if (is.na(pKey) || is.na(mKey)) next
    predMat <- twasCvResult$prediction[[pKey]]
    perfMat <- twasCvResult$performance[[mKey]]
    if (is.null(predMat) || is.null(perfMat)) next

    tr  <- traits[[i]]
    predCols <- colnames(predMat)
    perfRows <- rownames(perfMat)
    colHit <- if (!is.null(predCols) && tr %in% predCols) tr
              else if (ncol(predMat) == 1L) 1L else NA_integer_
    rowHit <- if (!is.null(perfRows) && tr %in% perfRows) tr
              else if (nrow(perfMat) == 1L) 1L else NA_integer_
    if (is.na(colHit) || is.na(rowHit)) next

    predVec <- predMat[, colHit, drop = TRUE]
    metRow  <- perfMat[rowHit, , drop = TRUE]
    cv <- list(
      samplePartition = twasCvResult$samplePartition,
      predictions     = predVec,
      metrics         = metRow)
    entry <- newEntries[[i]]
    newEntries[[i]] <- TwasWeightsEntry(
      variantIds    = getVariantIds(entry),
      weights       = getWeights(entry),
      fits          = getFits(entry),
      cvPerformance = cv,
      standardized  = getStandardized(entry),
      dataType      = getDataType(entry))
  }

  TwasWeights(
    study    = studies,
    context  = contexts,
    trait    = traits,
    method   = methodsV,
    entry    = newEntries,
    ldSketch = ldSketch)
}

# Mapping from short / canonical TWAS weight-method name to dispatch
# capability. Used to reject incompatible (input class, method) pairs.
#
# `allowsIndiv`  : may be invoked on a QtlDataset (individual-level X, Y).
# `allowsRss`    : may be invoked on a QtlSumStats / GwasSumStats (RSS).
# `multivariate` : requires a multi-trait / multi-context Y (mvsusie /
#                  mr.mash family).
#
# Rules from `dev/refactor-design.md` (`twasWeightsPipeline` row):
# - PRS-CS is RSS-only.
# - BGLR / CRAN-stable qgg methods (bayes_a/b/c/l/n/r, b_lasso, dpr_*)
#   are individual-level only.
# - mr.mash / mvsusie follow the multi-trait / multi-context rules of
#   the mvSuSiE fine-mapping family.
# @noRd
# User-facing TWAS method tokens are unified across input classes;
# auto-dispatch picks the individual-level vs sumstat implementation based
# on the QtlDataset / QtlSumStats input. Each entry records:
#   individualImpl  Function name to call on QtlDataset input (NULL = not
#                   supported on individual-level input).
#   sumstatImpl     Function name to call on QtlSumStats input (NULL = not
#                   supported on sumstat input).
#   multivariate    Whether the method requires multi-trait / multi-context
#                   structure (mvsusie / mrmash / mvsusieRss / mrmashRss).
#
# Per the design: BGLR / qgg "Bayes alphabet" methods (bayes_a/b/c/l/n/r,
# b_lasso) are individual-only until the qgg CRAN release adds qBayes
# sumstat support. dpr_gibbs has the SDPR sumstat counterpart;
# dpr_vb / dpr_adaptive_gibbs remain individual-only. enet has no cpp11
# sumstat solver yet (lassosumRssRcpp is pure L1, no alpha mixing) and is
# documented as individual-only for now. prsCs has no individual-level
# counterpart (it is a sumstat-only Bayesian shrinkage method).
.twasMethodCapabilities <- list(
  susie               = list(individualImpl = "susieWeights",
                             sumstatImpl    = "susieRssWeights",
                             multivariate   = FALSE),
  susieInf            = list(individualImpl = "susieInfWeights",
                             sumstatImpl    = "susieInfRssWeights",
                             multivariate   = FALSE),
  susieAsh            = list(individualImpl = "susieAshWeights",
                             sumstatImpl    = "susieAshRssWeights",
                             multivariate   = FALSE),
  mrash               = list(individualImpl = "mrashWeights",
                             sumstatImpl    = "mrAshRssWeights",
                             multivariate   = FALSE),
  lasso               = list(individualImpl = "lassoWeights",
                             sumstatImpl    = "lassosumRssWeights",
                             multivariate   = FALSE),
  scad                = list(individualImpl = "scadWeights",
                             sumstatImpl    = "scadRssWeights",
                             multivariate   = FALSE),
  mcp                 = list(individualImpl = "mcpWeights",
                             sumstatImpl    = "mcpRssWeights",
                             multivariate   = FALSE),
  l0learn             = list(individualImpl = "l0learnWeights",
                             sumstatImpl    = "l0learnRssWeights",
                             multivariate   = FALSE),
  mvsusie             = list(individualImpl = "mvsusieWeights",
                             sumstatImpl    = "mvsusieRssWeights",
                             multivariate   = TRUE),
  mrmash              = list(individualImpl = "mrmashWeights",
                             sumstatImpl    = "mrmashRssWeights",
                             multivariate   = TRUE),
  dpr_gibbs           = list(individualImpl = "dprGibbsWeights",
                             sumstatImpl    = "sdprWeights",
                             multivariate   = FALSE),
  # Individual-only — no cpp11 sumstat solver yet.
  enet                = list(individualImpl = "enetWeights",
                             sumstatImpl    = NULL,
                             multivariate   = FALSE),
  # Individual-only DPR variants (sumstat counterparts not implemented).
  dpr_vb              = list(individualImpl = "dprVbWeights",
                             sumstatImpl    = NULL,
                             multivariate   = FALSE),
  dpr_adaptive_gibbs  = list(individualImpl = "dprAdaptiveGibbsWeights",
                             sumstatImpl    = NULL,
                             multivariate   = FALSE),
  # qgg Bayes alphabet — individual-only until qgg CRAN release.
  bayes_a             = list(individualImpl = "bayesAWeights",
                             sumstatImpl    = NULL,
                             multivariate   = FALSE),
  bayes_b             = list(individualImpl = "bayesBWeights",
                             sumstatImpl    = NULL,
                             multivariate   = FALSE),
  bayes_c             = list(individualImpl = "bayesCWeights",
                             sumstatImpl    = NULL,
                             multivariate   = FALSE),
  bayes_l             = list(individualImpl = "bLassoWeights",
                             sumstatImpl    = NULL,
                             multivariate   = FALSE),
  bayes_n             = list(individualImpl = "bayesNWeights",
                             sumstatImpl    = NULL,
                             multivariate   = FALSE),
  bayes_r             = list(individualImpl = "bayesRWeights",
                             sumstatImpl    = NULL,
                             multivariate   = FALSE),
  b_lasso             = list(individualImpl = "bLassoWeights",
                             sumstatImpl    = NULL,
                             multivariate   = FALSE),
  # Sumstat-only Bayesian shrinkage (no individual-level analogue).
  prsCs               = list(individualImpl = NULL,
                             sumstatImpl    = "prsCsWeights",
                             multivariate   = FALSE))

# Normalize a user-supplied `methods` argument (character vector, preset
# string, or named list per `.twasMethodLookup`) into a (token, args) pair
# suitable for `.twasWeightsPipelineMatrix` / the sumstat sub-pipelines.
# Returns a list with `tokens` (canonical short names, used for capability
# lookup) and `methodList` (the `<token>_weights = args` list passed to
# `learnTwasWeights` / sumstat helpers).
# @noRd
.twasNormalizeMethods <- function(methods) {
  if (is.null(methods)) {
    methodList <- .twasMethodLookup("default")
  } else if (is.character(methods)) {
    methodList <- .twasMethodLookup(methods)
  } else if (is.list(methods)) {
    methodList <- methods
  } else {
    stop("`methods` must be a character vector, preset string, or named list.")
  }
  tokens <- sub("(_weights|Weights)$", "", names(methodList))
  list(tokens = tokens, methodList = methodList)
}

# Enforce input-class / method compatibility. Errors with a clear listing
# of incompatible methods and what input would accept them.
# @noRd
.twasCheckMethodCapabilities <- function(tokens, inputKind) {
  if (length(tokens) == 0L) return(invisible(NULL))
  caps <- .twasMethodCapabilities
  bad <- character(0)
  reason <- character(0)
  unknownTokens <- character(0)
  for (tk in tokens) {
    info <- caps[[tk]]
    if (is.null(info)) {
      unknownTokens <- c(unknownTokens, tk)
      next
    }
    if (inputKind == "QtlDataset" || inputKind == "MultiTaskQtlDataset") {
      if (is.null(info$individualImpl)) {
        bad <- c(bad, tk)
        reason <- c(reason,
                    "is sumstat-only (use a QtlSumStats input)")
      }
    } else if (inputKind == "QtlSumStats") {
      if (is.null(info$sumstatImpl)) {
        bad <- c(bad, tk)
        reason <- c(reason,
                    "is individual-only (use a QtlDataset input)")
      }
    }
  }
  if (length(unknownTokens) > 0L) {
    stop(sprintf(
      "twasWeightsPipeline: unknown method token(s): %s. Known tokens: %s.",
      paste(unknownTokens, collapse = ", "),
      paste(names(caps), collapse = ", ")))
  }
  if (length(bad) > 0L) {
    stop(sprintf(
      "twasWeightsPipeline: the following method(s) are not available for input class '%s': %s. %s.",
      inputKind,
      paste(bad, collapse = ", "),
      paste(sprintf("%s %s", bad, reason), collapse = "; ")))
  }
}

# Enforce the multi-trait / multi-context rule for mvsusie / mr.mash
# methods (same family as the fine-mapping mvSuSiE rule in the design
# doc). Multivariate methods need at least 2 traits *or* at least 2
# contexts in the Y matrix passed to learnTwasWeights.
# @noRd
.twasCheckMultivariateY <- function(tokens, nTraits, nContexts) {
  caps <- .twasMethodCapabilities
  multivariateTokens <- tokens[vapply(tokens, function(tk) {
    info <- caps[[tk]]
    !is.null(info) && isTRUE(info$multivariate)
  }, logical(1))]
  if (length(multivariateTokens) == 0L) return(invisible(NULL))
  if (nTraits < 2L && nContexts < 2L) {
    stop(sprintf(
      "twasWeightsPipeline: method(s) %s require multi-trait or multi-context input (got %d trait(s) x %d context(s)).",
      paste(multivariateTokens, collapse = ", "),
      nTraits, nContexts))
  }
}

# Reject SumStats inputs that have not been QC'd via summaryStatsQc.
# @noRd
.twasAssertQcd <- function(sumstats) {
  if (length(getQcInfo(sumstats)) == 0L) {
    stop("twasWeightsPipeline: the supplied ",
         class(sumstats)[[1L]],
         " has no QC record (qcInfo is empty). Call summaryStatsQc() ",
         "first and pass the QC-applied result.")
  }
}

# Extract a correlation matrix from a GenotypeHandle (LD sketch) for the
# variant subset given by `variantIds`. Reads dosage and computes
# sample-based LD via computeLd().
# @noRd
.twasLdFromSketch <- function(ldSketch, variantIds) {
  if (!is(ldSketch, "GenotypeHandle")) {
    stop(".twasLdFromSketch: ldSketch must be a GenotypeHandle.")
  }
  snpInfo <- getSnpInfo(ldSketch)
  snpAll  <- as.character(snpInfo$SNP)
  idx <- match(variantIds, snpAll)
  if (anyNA(idx)) {
    stop(sprintf(
      ".twasLdFromSketch: %d variant id(s) are not present in the LD sketch panel.",
      sum(is.na(idx))))
  }
  block <- extractBlockGenotypes(ldSketch, idx, meanImpute = TRUE)
  geno <- t(SummarizedExperiment::assay(block, "dosage"))
  # Subset / order to the requested variantIds.
  colnames(geno) <- as.character(snpInfo$SNP[idx])
  ldMat <- computeLd(geno, method = "sample")
  dimnames(ldMat) <- list(variantIds, variantIds)
  ldMat
}

# Helper to convert a single QtlSumStats / GwasSumStats entry GRanges to
# the data.frame shape (variant_id, chrom, pos, A1, A2, z, N) expected by
# the matrix-based sumstat pipelines.
# @noRd
.twasSumstatsEntryToDf <- function(gr) {
  mc <- as.data.frame(S4Vectors::mcols(gr))
  df <- data.frame(
    variant_id = as.character(mc$SNP),
    chrom      = as.character(GenomicRanges::seqnames(gr)),
    pos        = as.integer(GenomicRanges::start(gr)),
    A1         = as.character(mc$A1),
    A2         = as.character(mc$A2),
    stringsAsFactors = FALSE)
  if ("Z"    %in% colnames(mc)) df$z    <- as.numeric(mc$Z)
  if ("BETA" %in% colnames(mc)) df$beta <- as.numeric(mc$BETA)
  if ("SE"   %in% colnames(mc)) df$se   <- as.numeric(mc$SE)
  if ("N"    %in% colnames(mc)) df$N    <- as.numeric(mc$N)
  if ("MAF"  %in% colnames(mc)) df$maf  <- as.numeric(mc$MAF)
  if (is.null(df$z) && !is.null(df$beta) && !is.null(df$se)) {
    df$z <- df$beta / df$se
  }
  df
}

# Convert a FineMappingResult (single-method susie/susie_inf row matched
# to the requested study/context/trait) into a `fittedModels` list
# suitable for `learnTwasWeights`. Pulls the trimmedFit from the matching
# entry. Returns a (possibly empty) list.
# @noRd
.twasFineMappingFits <- function(fineMappingResult, study, context, trait) {
  if (is.null(fineMappingResult)) return(list())
  if (!is(fineMappingResult, "FineMappingResultBase")) {
    stop("`fineMappingResult` must be a FineMappingResult or NULL.")
  }
  out <- list()
  methods <- as.character(fineMappingResult$method)
  for (canonical in c("susie", "susieInf", "susieAsh")) {
    candidates <- c(canonical,
                    paste0(tolower(substring(canonical, 1L, 1L)),
                           substring(canonical, 2L)),
                    gsub("([A-Z])", "_\\1", canonical))
    candidates <- tolower(candidates)
    idx <- which(tolower(methods) %in% candidates &
                 as.character(fineMappingResult$study)   == study &
                 as.character(fineMappingResult$context) == context &
                 as.character(fineMappingResult$trait)   == trait)
    if (length(idx) > 0L) {
      out[[canonical]] <- getTrimmedFit(fineMappingResult$entry[[idx[[1L]]]])
    }
  }
  out
}

#' TWAS Weights Pipeline
#'
#' S4-dispatched per-region pipeline for learning TWAS prediction weights.
#' Accepts:
#' \itemize{
#'   \item a \code{\link{QtlDataset}} for individual-level cohort fits;
#'   \item a \code{\link{QtlSumStats}} for per-trait RSS fits;
#'   \item a \code{\link{GwasSumStats}} for per-LD-block PRS-CS-style fits
#'         from GWAS summary statistics.
#' }
#'
#' Method-restriction rules (enforced):
#' \itemize{
#'   \item \code{mr.mash}, \code{mvsusie} follow the multi-trait /
#'         multi-context rules of the fine-mapping \code{mvsusie} family
#'         (require at least two traits OR at least two contexts).
#'   \item RSS-only methods (PRS-CS, \code{lassosumRss}, SDPR, all
#'         \code{*Rss} variants) are rejected on \code{QtlDataset}
#'         input.
#'   \item Individual-level-only methods (BGLR and CRAN-stable qgg:
#'         \code{bayes_a/b/c/l/n/r}, \code{b_lasso}, \code{dpr_*}) are
#'         rejected on \code{QtlSumStats} / \code{GwasSumStats} input.
#' }
#'
#' Both \code{QtlSumStats} and \code{GwasSumStats} inputs must have been
#' QC'd via \code{\link{summaryStatsQc}} first; otherwise an error is
#' raised pointing at that function.
#'
#' The returned \code{\link{TwasWeights}} collection's \code{ldSketch}
#' slot is set automatically: \code{NULL} for individual-level fits,
#' the input's \code{ldSketch} for RSS-derived fits.
#'
#' Optionally a \code{\link{FineMappingResult}} may be supplied as a
#' source of pre-fit SuSiE / SuSiE-inf / SuSiE-ash objects; their
#' \code{trimmedFit} payloads are passed through to \code{learnTwasWeights}
#' / the RSS sub-pipelines via the \code{fittedModels} slot, avoiding
#' a re-fit.
#'
#' @param data A \code{QtlDataset}, \code{QtlSumStats}, or
#'   \code{GwasSumStats}.
#' @param methods A character vector of short method names, a preset
#'   string (\code{"default"} or \code{"fast_default"}), or a named list
#'   of \code{<method>_weights = args} entries. For QtlSumStats / GwasSumStats
#'   inputs the default switches to the RSS preset
#'   (\code{c("susieRss", "susieInfRss", "lassosumRss", "prsCs", "sdpr")}).
#' @param contexts Optional character vector of contexts to restrict
#'   processing to (QtlDataset / QtlSumStats inputs). Default \code{NULL}
#'   (use all contexts).
#' @param traitId Optional character vector of trait identifiers to
#'   restrict processing to (QtlDataset / QtlSumStats inputs). Default
#'   \code{NULL}.
#' @param region Optional \code{GRanges} for QtlDataset trait selection.
#'   Mutually exclusive with \code{traitId}.
#' @param cisWindow For QtlDataset: cis-window (bp) around each trait's
#'   genomic position when extracting variants. Required when
#'   \code{traitId} is supplied.
#' @param fineMappingResult Optional \code{\link{FineMappingResult}}.
#'   When supplied, its SuSiE / SuSiE-inf / SuSiE-ash trimmed fits for
#'   the matching (study, context, trait) tuples are injected into
#'   \code{learnTwasWeights} via \code{fittedModels} so SuSiE-family
#'   weight methods reuse the prior fit instead of refitting.
#' @param cvFolds Integer. Cross-validation folds. Default 5. Set to 0
#'   to skip CV (and ensemble).
#' @param samplePartition Optional pre-defined CV partition data.frame.
#' @param maxCvVariants Maximum number of variants for CV. Default -1
#'   (no limit).
#' @param cvThreads Threads for CV parallelism. Default 1.
#' @param cvWeightMethods Optional override of methods used for CV.
#' @param ensemble Logical. Compute SR-TWAS ensemble weights. Default
#'   \code{TRUE}.
#' @param ensembleR2Threshold Minimum CV R-squared for ensemble
#'   inclusion. Default 0.01.
#' @param ensembleSolver Solver for ensemble stacking. Default
#'   \code{"quadprog"}.
#' @param ensembleAlpha Elastic-net mixing parameter (only when
#'   \code{ensembleSolver = "glmnet"}). Default 1.
#' @param estimatePi If TRUE, estimate spike-and-slab sparsity from
#'   mr.ash before BGLR / qgg spike-and-slab methods that consume it.
#' @param phenotypeCovariatesToResidualize,genotypeCovariatesToResidualize
#'   Pass-through to \code{\link{getResidualizedPhenotypes}} and
#'   \code{\link{getResidualizedGenotypes}} for QtlDataset input.
#'   Default \code{NULL} (use all covariates).
#' @param dataType Optional data-type tag stamped into every
#'   \code{TwasWeightsEntry$dataType} (e.g. \code{"expression"}).
#' @param verbose Verbosity (0 silent, 1 default, 2 includes external
#'   package messages).
#' @param ... Reserved for method-specific arguments.
#'
#' @return A \code{\link{TwasWeights}} collection keyed by
#'   \code{(study, context, trait, method)}. The \code{ldSketch} slot is
#'   \code{NULL} for individual-level fits and equals the input's
#'   \code{ldSketch} for RSS-derived fits.
#' @export
setGeneric("twasWeightsPipeline",
  function(data, ...) standardGeneric("twasWeightsPipeline"))

#' @rdname twasWeightsPipeline
#' @export
setMethod("twasWeightsPipeline", "QtlDataset",
  function(data,
           methods                = "default",
           contexts               = NULL,
           traitId                = NULL,
           region                 = NULL,
           cisWindow              = NULL,
           fineMappingResult      = NULL,
           cvFolds                = 5,
           samplePartition        = NULL,
           maxCvVariants          = -1,
           cvThreads              = 1,
           cvWeightMethods        = NULL,
           ensemble               = TRUE,
           ensembleR2Threshold    = 0.01,
           ensembleSolver         = "quadprog",
           ensembleAlpha          = 1,
           estimatePi             = TRUE,
           phenotypeCovariatesToResidualize = NULL,
           genotypeCovariatesToResidualize  = NULL,
           dataType               = NULL,
           verbose                = 1,
           ...) {
    norm <- .twasNormalizeMethods(methods)
    .twasCheckMethodCapabilities(norm$tokens, "QtlDataset")

    study <- getStudy(data)
    allCtx <- getContexts(data)
    useCtx <- if (is.null(contexts)) allCtx else {
      bad <- setdiff(contexts, allCtx)
      if (length(bad) > 0L)
        stop("twasWeightsPipeline(QtlDataset): unknown context(s): ",
             paste(bad, collapse = ", "))
      contexts
    }

    # Collect traits to iterate over. When traitId is specified, use it;
    # when region is specified, use the per-context overlap; when neither
    # is supplied, iterate over every trait in every selected context.
    perCtxTraits <- vector("list", length(useCtx))
    names(perCtxTraits) <- useCtx
    for (ctx in useCtx) {
      se <- getPhenotypes(data, contexts = ctx)[[ctx]]
      ids <- rownames(se)
      if (!is.null(traitId)) {
        ids <- intersect(ids, traitId)
      } else if (!is.null(region)) {
        rr <- SummarizedExperiment::rowRanges(se)
        ids <- ids[IRanges::overlapsAny(rr, region)]
      }
      perCtxTraits[[ctx]] <- ids
    }
    allTraits <- unique(unlist(perCtxTraits))
    if (length(allTraits) == 0L) {
      stop("twasWeightsPipeline(QtlDataset): no traits selected.")
    }

    # Multivariate guard: gate on (nTraits, nContexts).
    nCtx <- length(useCtx)
    .twasCheckMultivariateY(norm$tokens, length(allTraits), nCtx)

    multivariate <- any(vapply(norm$tokens, function(tk) {
      info <- .twasMethodCapabilities[[tk]]
      !is.null(info) && isTRUE(info$multivariate)
    }, logical(1)))

    runOne <- function(ctx, tid) {
      X <- getResidualizedGenotypes(
        data, contexts = ctx, traitId = tid,
        cisWindow = cisWindow,
        phenotypeCovariatesToResidualize = phenotypeCovariatesToResidualize,
        genotypeCovariatesToResidualize  = genotypeCovariatesToResidualize)
      Yres <- getResidualizedPhenotypes(
        data, contexts = ctx, traitId = tid,
        phenotypeCovariatesToResidualize = phenotypeCovariatesToResidualize,
        genotypeCovariatesToResidualize  = genotypeCovariatesToResidualize)
      Y <- Yres[[ctx]]
      common <- intersect(rownames(X), rownames(Y))
      if (length(common) < 2L) {
        stop(sprintf(
          "twasWeightsPipeline: too few shared samples between residualized X and Y for (context='%s', trait='%s').",
          ctx, tid))
      }
      X <- X[common, , drop = FALSE]
      Y <- Y[common, , drop = FALSE]

      fittedModels <- .twasFineMappingFits(fineMappingResult,
                                            study = study,
                                            context = ctx,
                                            trait = tid)
      .twasWeightsPipelineMatrix(
        X = X, y = Y,
        study = study, context = ctx, trait = tid,
        fittedModels = fittedModels,
        cvFolds = cvFolds,
        samplePartition = samplePartition,
        weightMethods = norm$methodList,
        maxCvVariants = maxCvVariants,
        cvThreads = cvThreads,
        cvWeightMethods = cvWeightMethods,
        ensemble = ensemble,
        ensembleR2Threshold = ensembleR2Threshold,
        ensembleSolver = ensembleSolver,
        ensembleAlpha = ensembleAlpha,
        estimatePi = estimatePi,
        standardized = FALSE,
        dataType = dataType,
        ldSketch = NULL,
        verbose = verbose)$twasWeights
    }

    runMultivariate <- function(traits) {
      # Joint over selected (contexts, traits): residualize, intersect
      # samples across contexts, drop subjects with any-NA in Y.
      Xlist <- lapply(useCtx, function(ctx) {
        getResidualizedGenotypes(
          data, contexts = ctx, traitId = traits,
          cisWindow = cisWindow,
          phenotypeCovariatesToResidualize = phenotypeCovariatesToResidualize,
          genotypeCovariatesToResidualize  = genotypeCovariatesToResidualize)
      })
      # Intersect samples across contexts.
      commonSamples <- Reduce(intersect, lapply(Xlist, rownames))
      if (length(commonSamples) < 2L) {
        stop("twasWeightsPipeline(QtlDataset, multivariate): insufficient samples shared across selected contexts.")
      }
      X <- Xlist[[1L]][commonSamples, , drop = FALSE]

      Yres <- getResidualizedPhenotypes(
        data, contexts = useCtx, traitId = traits,
        phenotypeCovariatesToResidualize = phenotypeCovariatesToResidualize,
        genotypeCovariatesToResidualize  = genotypeCovariatesToResidualize)
      # Concatenate per-context residualized phenotypes column-wise,
      # restricting to commonSamples. Column names become
      # "<context>__<trait>".
      Ymats <- list()
      colMeta <- list()
      for (ctx in names(Yres)) {
        Ym <- Yres[[ctx]][commonSamples, , drop = FALSE]
        colnames(Ym) <- paste(ctx, colnames(Ym), sep = "__")
        Ymats[[ctx]] <- Ym
        colMeta[[ctx]] <- data.frame(
          context = ctx, trait = colnames(Yres[[ctx]]),
          stringsAsFactors = FALSE)
      }
      Y <- do.call(cbind, Ymats)
      meta <- do.call(rbind, colMeta)
      # Drop subjects with any NA across Y columns.
      keep <- complete.cases(Y)
      if (sum(keep) < 2L) {
        stop("twasWeightsPipeline(QtlDataset, multivariate): too few subjects with complete Y across selected (context, trait) columns.")
      }
      Y <- Y[keep, , drop = FALSE]
      X <- X[rownames(Y), , drop = FALSE]

      # Build per-column identity tuples for learnTwasWeights so multi-
      # outcome methods emit one row per (context, trait).
      .twasWeightsPipelineMatrix(
        X = X, y = Y,
        study   = study,
        context = meta$context,
        trait   = meta$trait,
        cvFolds = cvFolds,
        samplePartition = samplePartition,
        weightMethods = norm$methodList,
        maxCvVariants = maxCvVariants,
        cvThreads = cvThreads,
        cvWeightMethods = cvWeightMethods,
        ensemble = ensemble,
        ensembleR2Threshold = ensembleR2Threshold,
        ensembleSolver = ensembleSolver,
        ensembleAlpha = ensembleAlpha,
        estimatePi = estimatePi,
        standardized = FALSE,
        dataType = dataType,
        ldSketch = NULL,
        verbose = verbose)$twasWeights
    }

    # Top-level dispatch within the QtlDataset method body.
    if (multivariate) {
      # mvsusie / mr.mash: joint fit. If both nCtx == 1 and nTraits == 1
      # we already rejected above via .twasCheckMultivariateY.
      tw <- runMultivariate(allTraits)
    } else {
      # Univariate methods: sequential over (context, trait).
      out <- NULL
      for (ctx in useCtx) {
        for (tid in perCtxTraits[[ctx]]) {
          twi <- runOne(ctx, tid)
          out <- if (is.null(out)) twi else .rbindTwasWeights(out, twi, ldSketch = NULL)
        }
      }
      tw <- out
    }
    if (is.null(tw)) {
      stop("twasWeightsPipeline(QtlDataset): no (context, trait) pair produced any weights.")
    }
    tw
  })

#' @rdname twasWeightsPipeline
#' @export
setMethod("twasWeightsPipeline", "QtlSumStats",
  function(data,
           methods           = NULL,
           contexts          = NULL,
           traitId           = NULL,
           dataType          = NULL,
           verbose           = 1L,
           ...) {
    # summaryStatsQc() is mandatory before twasWeightsPipeline for SumStats
    # input; it also drops variants not present in the ldSketch, so by the
    # time we reach this method every entry's SNP set is a subset of the
    # ldSketch panel.
    .twasAssertQcd(data)

    # Normalize the methods argument into (tokens, methodArgs).
    if (is.null(methods)) {
      tokens <- c("susie", "susieInf", "lasso", "prsCs", "dpr_gibbs")
      methodArgs <- setNames(rep(list(list()), length(tokens)), tokens)
    } else if (is.character(methods)) {
      tokens <- methods
      methodArgs <- setNames(rep(list(list()), length(tokens)), tokens)
    } else if (is.list(methods)) {
      tokens <- names(methods)
      methodArgs <- methods
    } else {
      stop("`methods` must be NULL, a character vector, or a named list ",
           "of <token> = <args> entries.")
    }
    .twasCheckMethodCapabilities(tokens, "QtlSumStats")

    studyCol   <- as.character(data$study)
    contextCol <- as.character(data$context)
    traitCol   <- as.character(data$trait)

    selRows <- seq_len(nrow(data))
    if (!is.null(contexts)) selRows <- selRows[contextCol[selRows] %in% contexts]
    if (!is.null(traitId))  selRows <- selRows[traitCol[selRows]   %in% traitId]
    if (length(selRows) == 0L) {
      stop("twasWeightsPipeline(QtlSumStats): no entries matched the ",
           "supplied contexts / traitId filters.")
    }

    # Partition method tokens by univariate vs multivariate dispatch.
    isMv <- vapply(tokens, function(tk) {
      isTRUE(.twasMethodCapabilities[[tk]]$multivariate)
    }, logical(1))
    multivariateTokens <- tokens[isMv]
    univariateTokens   <- tokens[!isMv]

    if (length(multivariateTokens) > 0L) {
      groupKey <- paste(studyCol[selRows], traitCol[selRows], sep = "||")
      perGroupNCtx <- vapply(split(contextCol[selRows], groupKey),
                             length, integer(1))
      if (all(perGroupNCtx < 2L)) {
        stop(sprintf(
          "twasWeightsPipeline(QtlSumStats): multivariate method(s) %s require at least two contexts per (study, trait); the supplied collection has only one context per trait.",
          paste(multivariateTokens, collapse = ", ")))
      }
    }

    ldSketch <- getLdSketch(data)

    rowStudy   <- character(0)
    rowContext <- character(0)
    rowTrait   <- character(0)
    rowMethod  <- character(0)
    rowEntries <- list()

    # ---- Univariate dispatch: per (study, context, trait), per method.
    for (i in selRows) {
      st <- studyCol[i]; ctx <- contextCol[i]; tr <- traitCol[i]
      entry <- data$entry[[i]]
      mc <- S4Vectors::mcols(entry)
      variantIds <- as.character(mc$SNP)
      if (!"Z" %in% colnames(mc))
        stop(sprintf(
          "twasWeightsPipeline(QtlSumStats): entry %d (study='%s', context='%s', trait='%s') has no Z column.",
          i, st, ctx, tr))
      if (!"N" %in% colnames(mc))
        stop(sprintf(
          "twasWeightsPipeline(QtlSumStats): entry %d (study='%s', context='%s', trait='%s') has no N column.",
          i, st, ctx, tr))
      n <- stats::median(as.numeric(mc$N), na.rm = TRUE)
      varY <- getVarY(data, study = st, context = ctx, trait = tr)
      if (is.null(varY)) varY <- 1
      stat <- list(z = as.numeric(mc$Z), n = n, varY = varY,
                   variantNames = variantIds)
      ldMat <- .twasLdFromSketch(ldSketch, variantIds)

      for (tk in univariateTokens) {
        fn <- .twasMethodCapabilities[[tk]]$sumstatImpl
        userArgs <- methodArgs[[tk]]
        if (is.null(userArgs)) userArgs <- list()
        weights <- tryCatch(
          do.call(get(fn, mode = "function"),
                  c(list(stat = stat, LD = ldMat), userArgs)),
          error = function(e) {
            warning(sprintf(
              "twasWeightsPipeline: method '%s' failed for (study=%s, context=%s, trait=%s): %s",
              tk, st, ctx, tr, conditionMessage(e)))
            NULL
          })
        if (is.null(weights)) next
        fitAttr <- attr(weights, "fit")
        attr(weights, "fit") <- NULL
        rowStudy   <- c(rowStudy,   st)
        rowContext <- c(rowContext, ctx)
        rowTrait   <- c(rowTrait,   tr)
        rowMethod  <- c(rowMethod,  tk)
        rowEntries[[length(rowEntries) + 1L]] <- TwasWeightsEntry(
          variantIds    = variantIds,
          weights       = as.numeric(weights),
          fits          = fitAttr,
          cvPerformance = NULL,        # Q5: no CV on the sumstat path
          standardized  = TRUE,        # Q4: sumstat-derived weights are standardized
          dataType      = dataType)
      }
    }

    # ---- Multivariate dispatch: per (study, trait), all selected contexts.
    if (length(multivariateTokens) > 0L) {
      groupKey <- paste(studyCol[selRows], traitCol[selRows], sep = "||")
      groups   <- split(selRows, groupKey)
      for (gkey in names(groups)) {
        gIdx <- groups[[gkey]]
        if (length(gIdx) < 2L) next
        st <- studyCol[gIdx[[1L]]]
        tr <- traitCol[gIdx[[1L]]]
        ctxNames <- contextCol[gIdx]

        # Build (variants x contexts) Z matrix. All entries in a (study, trait)
        # group must share an identical variant order after summaryStatsQc().
        firstMc <- S4Vectors::mcols(data$entry[[gIdx[[1L]]]])
        variantIds <- as.character(firstMc$SNP)
        Z <- matrix(NA_real_, nrow = length(variantIds), ncol = length(gIdx),
                    dimnames = list(variantIds, ctxNames))
        nVec <- numeric(length(gIdx))
        for (kk in seq_along(gIdx)) {
          mc <- S4Vectors::mcols(data$entry[[gIdx[kk]]])
          if (!identical(as.character(mc$SNP), variantIds))
            stop("twasWeightsPipeline(QtlSumStats, multivariate): every ",
                 "entry for (study='", st, "', trait='", tr,
                 "') must share an identical SNP order after ",
                 "summaryStatsQc(). Use the same ldSketch on every entry.")
          Z[, kk] <- as.numeric(mc$Z)
          nVec[kk] <- stats::median(as.numeric(mc$N), na.rm = TRUE)
        }
        names(nVec) <- ctxNames
        stat <- list(z = Z, n = nVec, variantNames = variantIds)
        ldMat <- .twasLdFromSketch(ldSketch, variantIds)

        for (tk in multivariateTokens) {
          fn <- .twasMethodCapabilities[[tk]]$sumstatImpl
          userArgs <- methodArgs[[tk]]
          if (is.null(userArgs)) userArgs <- list()
          weights <- tryCatch(
            do.call(get(fn, mode = "function"),
                    c(list(stat = stat, LD = ldMat), userArgs)),
            error = function(e) {
              warning(sprintf(
                "twasWeightsPipeline: multivariate method '%s' failed for (study=%s, trait=%s): %s",
                tk, st, tr, conditionMessage(e)))
              NULL
            })
          if (is.null(weights)) next
          if (!is.matrix(weights)) weights <- as.matrix(weights)
          fitAttr <- attr(weights, "fit")
          attr(weights, "fit") <- NULL
          for (kk in seq_along(ctxNames)) {
            rowStudy   <- c(rowStudy,   st)
            rowContext <- c(rowContext, ctxNames[[kk]])
            rowTrait   <- c(rowTrait,   tr)
            rowMethod  <- c(rowMethod,  tk)
            rowEntries[[length(rowEntries) + 1L]] <- TwasWeightsEntry(
              variantIds    = variantIds,
              weights       = as.numeric(weights[, kk]),
              # Share the underlying joint fit on the first row only;
              # remaining rows reference the same fit by leaving fits NULL.
              fits          = if (kk == 1L) fitAttr else NULL,
              cvPerformance = NULL,
              standardized  = TRUE,
              dataType      = dataType)
          }
        }
      }
    }

    if (length(rowEntries) == 0L) {
      stop("twasWeightsPipeline(QtlSumStats): no entries produced weights.")
    }
    TwasWeights(
      study    = rowStudy,
      context  = rowContext,
      trait    = rowTrait,
      method   = rowMethod,
      entry    = rowEntries,
      ldSketch = ldSketch)
  })

#' @rdname twasWeightsPipeline
#' @export
setMethod("twasWeightsPipeline", "ANY",
  function(data, ...) {
    stop("twasWeightsPipeline does not accept inputs of class '",
         class(data)[[1L]], "'. Pass a QtlDataset, MultiTaskQtlDataset, ",
         "or QtlSumStats. (GwasSumStats inputs are not supported; ",
         "GWAS-side per-LD-block weights are produced inside the new ",
         "ctwasPipeline / qtlEnrichmentPipeline.)")
  })

# =============================================================================
# Internal matrix-driven TWAS weights pipeline
# =============================================================================
#
# This is the legacy matrix-based pipeline retained as an internal worker.
# The exported, S4-dispatched `twasWeightsPipeline` defined above extracts
# (X, Y) blocks from QtlDataset / QtlSumStats / GwasSumStats and calls this
# function per (study, context, trait) tuple. It returns a single-tuple
# `TwasWeights` collection (one row per method, plus an optional ensemble
# row) along with auxiliary state used during stacking.
#
# Method restrictions imposed at the dispatch layer:
# - PRS-CS, lassosumRss, sdpr, susieRss, susieInfRss, susieAshRss,
#   mrAshRss, mrmashRss, mvsusieRss: RSS-only (refuse QtlDataset).
# - bglrWeights / qgg methods (bayesA/B/C/L/N/R, bLasso, dpr*): individual
#   level only (refuse QtlSumStats / GwasSumStats).
# - mr.mash / mvsusie: multi-trait / multi-context (same rule family as
#   the fine-mapping mvSuSiE family in the design doc).
#
# @noRd
.twasWeightsPipelineMatrix <- function(X,
                                y,
                                study = "",
                                context = "",
                                trait = "",
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
                                standardized = FALSE,
                                dataType = NULL,
                                ldSketch = NULL,
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
    res$susieWeightsIntermediate <- .susieWeightIntermediate(fittedModels[["susie"]], X)
  }

  # Check if empirical pi estimation is needed for spike-and-slab methods
  bayesCneedsPi <- "bayes_c_weights" %in% names(weightMethods) &&
    !"pi" %in% names(weightMethods$bayes_c_weights)
  bayesBneedsPi <- "bayes_b_weights" %in% names(weightMethods) &&
    !"probIn" %in% names(weightMethods$bayes_b_weights)
  needsPiEstimation <- (bayesCneedsPi || bayesBneedsPi) && estimatePi

  learnArgs <- list(
    study = study, context = context, trait = trait,
    standardized = standardized, dataType = dataType,
    ldSketch = ldSketch)

  if (needsPiEstimation) {
    # Run mr.ash first to estimate sparsity
    mrashMethods <- list(mrash_weights = weightMethods[["mrash_weights"]] %||% list())

    if (verbose >= 1) message("  Estimating sparsity from mr.ash ...")
    mrashWeights <- do.call(learnTwasWeights, c(
      list(X = X, Y = y, weightMethods = mrashMethods,
           retainFits = TRUE, verbose = verbose),
      learnArgs))

    empiricalPi <- estimateSparsity(mrashWeights)
    if (verbose >= 1) message(sprintf("  Empirical sparsity estimate: %.4f", empiricalPi))
    res$empiricalPi <- empiricalPi

    # Inject into spike-and-slab methods that need it
    if (bayesCneedsPi) weightMethods$bayes_c_weights$pi <- as.numeric(empiricalPi)
    if (bayesBneedsPi) weightMethods$bayes_b_weights$probIn <- as.numeric(empiricalPi)

    # Run remaining methods (those not already computed)
    remainingFnNames <- setdiff(names(weightMethods), "mrash_weights")

    if (length(remainingFnNames) > 0) {
      remainingMethods <- weightMethods[remainingFnNames]
      remainingTw <- do.call(learnTwasWeights, c(
        list(X = X, Y = y, weightMethods = remainingMethods,
             fittedModels = fittedModels, verbose = verbose),
        learnArgs))
      res$twasWeights <- .rbindTwasWeights(mrashWeights, remainingTw,
                                            ldSketch = ldSketch)
    } else {
      res$twasWeights <- mrashWeights
    }

    # Remove mr.ash if it was not in the original weightMethods
    if (!"mrash_weights" %in% names(weightMethods)) {
      tw <- res$twasWeights
      keep <- as.character(tw$method) != "mrash"
      res$twasWeights <- TwasWeights(
        study   = as.character(tw$study)[keep],
        context = as.character(tw$context)[keep],
        trait   = as.character(tw$trait)[keep],
        method  = as.character(tw$method)[keep],
        entry   = as.list(tw$entry)[keep],
        ldSketch = ldSketch)
    }
  } else {
    # Run all methods at once
    res$twasWeights <- do.call(learnTwasWeights, c(
      list(X = X, Y = y, weightMethods = weightMethods,
           fittedModels = fittedModels, verbose = verbose),
      learnArgs))
  }
  if (verbose >= 1) {
    elapsed <- toc(quiet = TRUE)
    message(sprintf("TWAS weights fitting done in %.1fs", elapsed$toc - elapsed$tic))
  }
  res$twasPredictions <- twasPredict(X, res$twasWeights)

  if (cvFolds > 1) {
    # A few cutting corners to run CV faster at the disadvantage of SuSiE and mr.ash:
    # 1. reset SuSiE to not using refine or adaptive L but to use L from previous analysis
    # 2. at most 100 iterations for mr.ash allowed
    # 3. only use a subset of variants randomly selected to avoid bias
    if (!is.null(fittedModels[["susieInf"]]) && !is.null(weightMethods$susie_inf_weights)) {
      weightMethods$susie_inf_weights$L <- length(fittedModels[["susieInf"]]$V)
      weightMethods$susie_inf_weights$refine <- FALSE
    }
    if (!is.null(weightMethods$susie_weights)) {
      susieCvFit <- fittedModels[["susie"]]
      if (is.null(susieCvFit)) susieCvFit <- fittedModels[["susieInf"]]
      if (!is.null(susieCvFit)) {
        weightMethods$susie_weights$L <- length(susieCvFit$V)
        weightMethods$susie_weights$refine <- FALSE
      }
    }
    if (is.null(cvWeightMethods)) {
      cvWeightMethods <- .filterZeroWeightMethods(weightMethods, res$twasWeights)
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
    res$twasCvResult <- twasWeightsCv(
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

    # Splice per-(method, outcome) CV predictions + metrics into the
    # corresponding TwasWeightsEntry$cvPerformance slot.
    res$twasWeights <- .spliceCvIntoTwasWeights(res$twasWeights,
                                                 res$twasCvResult,
                                                 ldSketch = ldSketch)

    # Ensemble learning: learn optimal method combination via stacked regression
    if (isTRUE(ensemble) && length(cvWeightMethods) <= 1) {
      if (verbose >= 1) message("Ensemble model skipped: only ", length(cvWeightMethods),
              " weight method provided (need >= 2 for ensemble learning).")
    }
    if (isTRUE(ensemble) && length(cvWeightMethods) > 1) {
      if (!is.null(res$twasCvResult$performance)) {
        # Extract R-squared for each method from CV performance table
        methodRsq <- vapply(res$twasCvResult$performance, function(perf) {
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
          filteredCv <- res$twasCvResult
          predBaseNames <- sub("(_predicted|Predicted)$", "", names(filteredCv$prediction))
          filteredCv$prediction <- filteredCv$prediction[match(passingBase, predBaseNames)]

          # Subset twas_weights to passing methods.
          # Method names on a TwasWeights collection are stored as bare
          # tokens (e.g. "lasso") in the `method` column; the ensemble
          # helper wants snake_case "<method>_weights" keys.
          tw <- res$twasWeights
          twMethodNames <- as.character(tw$method)
          filteredWeights <- setNames(
            lapply(passingBase, function(bn) {
              idx <- which(twMethodNames == bn)
              if (length(idx) == 0L) return(NULL)
              w <- getWeights(tw$entry[[idx[[1L]]]])
              if (!is.matrix(w)) w <- matrix(w, ncol = 1)
              w
            }),
            paste0(passingBase, "_weights"))
          filteredWeights <- Filter(Negate(is.null), filteredWeights)

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

          # Add ensemble weights alongside individual method weights as a
          # new row in the TwasWeights collection.
          if (!is.null(ensResult$ensembleTwasWeights)) {
            ensWt <- ensResult$ensembleTwasWeights
            if (!is.matrix(ensWt)) ensWt <- matrix(ensWt, ncol = 1)
            tw <- res$twasWeights
            # Use the first existing row's (study, context, trait) as the
            # identity tuple for the ensemble row.
            existingStudy   <- as.character(tw$study)[1L]
            existingContext <- as.character(tw$context)[1L]
            existingTrait   <- as.character(tw$trait)[1L]
            existingStd     <- getStandardized(tw$entry[[1L]])
            ensWtVec <- if (ncol(ensWt) == 1L) drop(ensWt) else ensWt
            ensVarIds <- if (!is.null(rownames(ensWt))) rownames(ensWt)
                         else colnames(X)
            ensEntry <- TwasWeightsEntry(
              variantIds   = ensVarIds,
              weights      = ensWtVec,
              cvPerformance = list(
                methodCoef        = ensResult$methodCoef,
                methodPerformance = ensResult$methodPerformance),
              standardized = existingStd)
            ensRow <- TwasWeights(
              study   = existingStudy,
              context = existingContext,
              trait   = existingTrait,
              method  = "ensemble",
              entry   = list(ensEntry),
              ldSketch = ldSketch)
            res$twasWeights <- .rbindTwasWeights(tw, ensRow, ldSketch = ldSketch)
            res$twasPredictions$ensemble_predicted <- X %*% ensWt
          }
          res$ensemble <- ensResult
        }
      }
    }
  }
  res$totalTimeElapsed <- proc.time() - st

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
          twasWeights = lapply(wl, function(wgts) wgts[, ctx]),
          twasPredictions = lapply(twasPredictions, function(pred) pred[, ctx]),
          variantNames = variantNames
        )
      } else {
        NULL
      }
    }), contextNames)
  }

  copyTwasCvResults <- function(twasResult, twasCvResult) {
    for (i in names(twasResult)) {
      if (i %in% colnames(twasCvResult$prediction[[1]])) {
        twasResult[[i]]$twasCvResult$samplePartition <- twasCvResult$samplePartition
        twasResult[[i]]$twasCvResult$prediction <- lapply(
          twasCvResult$prediction,
          function(predicted) {
            as.matrix(predicted[, i], ncol = 1)
          }
        )
        twasResult[[i]]$twasCvResult$performance <- lapply(
          twasCvResult$performance,
          function(perform) {
            t(as.matrix(perform[i, ], ncol = 1))
          }
        )
        twasResult[[i]]$twasCvResult$timeElapsed <- twasCvResult$timeElapsed
      }
    }
    return(twasResult)
  }

  # TWAS weights and predictions
  weightMethods <- list(
    mrmash_weights = list(
      mrmash_fit = mnmFit$mrmashFitted
    ),
    mvsusie_weights = list(
      mvsusie_fit = mnmFit$mvsusieFitted
    )
  )
  st <- proc.time()
  if (verbose >= 1) {
    message("Extracting TWAS weights for multivariate analysis methods ...")
    tic()
  }
  # get TWAS weights
  twasWeightsRes <- learnTwasWeights(X = X, Y = Y, weightMethods = weightMethods, verbose = verbose)
  if (verbose >= 1) {
    elapsed <- toc(quiet = TRUE)
    message(sprintf("Multivariate TWAS weights fitting done in %.1fs", elapsed$toc - elapsed$tic))
  }
  # get TWAS predictions for possible next steps such as computing correlations between predicted expression values
  twasPredictions <- twasPredict(X, twasWeightsRes)

  # copy TWAS results by condition
  res <- copyTwasResults(colnames(Y), mnmFit$variantNames, twasWeightsRes, twasPredictions)

  # Perform cross-validation if specified
  if (cvFolds > 1) {
    if (is.null(L)) L <- length(mnmFit$mvsusieFitted$V)
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
        prior_variance = mnmFit$reweightedMixturePrior,
        residual_variance = mnmFit$mrmashFitted$V,
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
      reweightedMixturePriorCv = mnmFit$reweightedMixturePriorCv
    )
    if (verbose >= 1) {
      elapsed <- toc(quiet = TRUE)
      message(sprintf("Cross-validation done in %.1fs", elapsed$toc - elapsed$tic))
    }
    res <- copyTwasCvResults(res, twasCvResult)
  }
  totalTimeElapsed <- proc.time() - st
  for (i in seq_along(res)) {
    res[[i]]$totalTimeElapsed <- totalTimeElapsed
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
#'   \code{\link{learnTwasWeights}}, with keys like \code{"susie_weights"}. Used to
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
#'   \item{methodCoef}{Named numeric vector of combination coefficients
#'     (\eqn{\zeta_k}), non-negative and summing to 1. Names are method
#'     base names (e.g., \code{"susie"}, \code{"enet"}).}
#'   \item{ensembleTwasWeights}{Final combined weight vector
#'     \eqn{w = \sum_k \zeta_k w_k}, or NULL if \code{twasWeightList}
#'     is not provided. Returned as a vector for univariate Y, matrix otherwise.}
#'   \item{methodPerformance}{Named numeric vector of per-method R-squared
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
#' @seealso \code{\link{twasWeightsCv}}, \code{\link{learnTwasWeights}},
#'   \code{\link{twasWeightsPipeline}}
#'
#' @examples
#' \dontrun{
#' # After running twasWeightsPipeline with CV:
#' res <- twasWeightsPipeline(X, y, cvFolds = 5, weightMethods = methods)
#'
#' ens <- ensembleWeights(
#'   cvResults = res$twasCvResult,
#'   Y = y,
#'   twasWeightList = res$twasWeights
#' )
#' ens$methodCoef           # combination weights, sum to 1
#'
#' # Multi-dataset ensemble (e.g., CUMC1 + MIT cell types):
#' ens_multi <- ensembleWeights(
#'   cvResults = list(res_cumc$twasCvResult, res_mit$twasCvResult),
#'   Y = list(y_cumc, y_mit),
#'   twasWeightList = list(res_cumc$twasWeights, res_mit$twasWeights)
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
    methodCoef = zeta,
    ensembleTwasWeights = ensembleTwasWt,
    methodPerformance = methodRsq
  )
}

