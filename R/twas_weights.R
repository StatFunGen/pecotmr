# Evaluate an expression while suppressing external package output.
# Catches both message() output (susieR, qgg) and Rprintf/cat stdout (mr.ash.alpha).
# @param expr An expression to evaluate.
# @return The result of evaluating expr.
# @noRd
.quiet_eval <- function(expr) {
  invisible(capture.output(
    result <- suppressMessages(expr),
    type = "output"
  ))
  result
}

# Map short method names and presets to weight_methods lists.
# @param methods A character vector of short method names, or a preset string
#   ("default" or "fast_default").
# @return A named list suitable for the weight_methods parameter.
# @noRd
.twas_method_lookup <- function(methods) {
  method_map <- list(
    susie = list(fn = "susie_weights", args = list(refine = FALSE, L = 20, L_greedy = 5)),
    susie_ash = list(fn = "susie_ash_weights", args = list()),
    susie_inf = list(fn = "susie_inf_weights", args = list()),
    mrash = list(fn = "mrash_weights", args = list(init_prior_sd = TRUE, max.iter = 100)),
    enet = list(fn = "enet_weights", args = list()),
    lasso = list(fn = "lasso_weights", args = list()),
    bayes_r = list(fn = "bayes_r_weights", args = list()),
    bayes_l = list(fn = "bayes_l_weights", args = list()),
    bayes_a = list(fn = "bayes_a_weights", args = list()),
    bayes_b = list(fn = "bayes_b_weights", args = list()),
    bayes_c = list(fn = "bayes_c_weights", args = list()),
    bayes_n = list(fn = "bayes_n_weights", args = list()),
    b_lasso = list(fn = "b_lasso_weights", args = list()),
    dpr_vb = list(fn = "dpr_vb_weights", args = list()),
    dpr_gibbs = list(fn = "dpr_gibbs_weights", args = list()),
    dpr_adaptive_gibbs = list(fn = "dpr_adaptive_gibbs_weights", args = list()),
    scad = list(fn = "scad_weights", args = list()),
    mcp = list(fn = "mcp_weights", args = list()),
    l0learn = list(fn = "l0learn_weights", args = list()),
    mvsusie = list(fn = "mvsusie_weights", args = list(L = 30, L_greedy = 5)),
    mrmash = list(fn = "mrmash_weights", args = list())
  )

  # Handle presets
  fast_default <- c("susie", "susie_inf", "mrash", "enet", "lasso", "mcp", "scad", "l0learn")
  if (length(methods) == 1) {
    if (methods == "fast_default") {
      methods <- fast_default
    } else if (methods == "default") {
      methods <- c(fast_default, "bayes_r", "bayes_c")
    }
  }

  # Build reverse map: function name -> short name, so full names are accepted too
  fn_to_short <- setNames(
    names(method_map),
    vapply(method_map, function(x) x$fn, character(1))
  )
  # Normalize any full function names to short names
  methods <- vapply(methods, function(m) {
    if (m %in% names(fn_to_short)) fn_to_short[[m]] else m
  }, character(1), USE.NAMES = FALSE)

  unknown <- setdiff(methods, names(method_map))
  if (length(unknown) > 0) {
    stop(
      "Unknown TWAS method(s): ", paste(unknown, collapse = ", "),
      ". Available methods: ", paste(names(method_map), collapse = ", ")
    )
  }

  result <- list()
  for (m in methods) {
    entry <- method_map[[m]]
    result[[entry$fn]] <- entry$args
  }
  result
}

# Identify non-zero-variance columns of X. Returns a logical vector.
#' @importFrom matrixStats colSds
#' @noRd
.nonzero_var_columns <- function(X) {
  sds <- colSds(X, na.rm = TRUE)
  !is.na(sds) & sds != 0
}

# Embed a smaller weights matrix into a full-sized zero matrix matching X and Y dimensions.
# @param weights_matrix The fitted weights (nrow = number of valid columns).
# @param valid_columns Logical or character vector identifying which columns of X were used.
# @param X_colnames Column names of the original X.
# @param Y_colnames Column names of Y.
# @noRd
.embed_weights <- function(weights_matrix, valid_columns, n_cols_X, n_cols_Y,
                           X_colnames = NULL, Y_colnames = NULL) {
  full <- matrix(0, nrow = n_cols_X, ncol = n_cols_Y)
  if (!is.null(X_colnames)) rownames(full) <- X_colnames
  if (!is.null(Y_colnames)) colnames(full) <- Y_colnames
  full[valid_columns, ] <- weights_matrix
  full
}

# Filter weight methods that produced all-zero weights from CV.
# Returns filtered weight_methods list and warns about removed methods.
# @noRd
.filter_zero_weight_methods <- function(weight_methods, twas_weights_res) {
  wl <- if (is(twas_weights_res, "TWASWeights")) getWeights(twas_weights_res) else twas_weights_res
  is_all_zero <- vapply(wl, function(w) all(w == 0, na.rm = TRUE), logical(1))
  removed <- names(weight_methods)[is_all_zero]
  if (length(removed) > 0) {
    warning(sprintf(
      "Methods %s are removed from CV because all their weights are zeros.",
      paste(removed, collapse = ", ")
    ))
  }
  weight_methods[!is_all_zero]
}

.susie_weight_intermediate <- function(fit, X) {
  keep <- intersect(c("mu", "lbf_variable", "X_column_scale_factors", "pip", "theta"), names(fit))
  intermediate <- fit[keep]
  if (!is.null(fit$sets$cs)) {
    intermediate$cs_variants <- setNames(lapply(fit$sets$cs, function(L) colnames(X)[L]), names(fit$sets$cs))
    intermediate$cs_purity <- fit$sets$purity
  }
  intermediate
}

.prepare_susie_weight_methods <- function(X, Y, weight_methods, fitted_models = NULL) {
  if (is.vector(Y)) Y <- matrix(Y, ncol = 1)
  if (is.null(fitted_models)) fitted_models <- list()
  has_susie <- !is.null(weight_methods[["susie_weights"]])
  has_susie_inf <- !is.null(weight_methods[["susie_inf_weights"]])
  susie_fit <- if (has_susie) weight_methods[["susie_weights"]][["susie_fit"]] else NULL
  susie_inf_fit <- if (has_susie_inf) weight_methods[["susie_inf_weights"]][["susie_inf_fit"]] else NULL
  if (is.null(susie_fit)) susie_fit <- fitted_models[["susie"]]
  if (is.null(susie_inf_fit)) susie_inf_fit <- fitted_models[["susie_inf"]]

  if (!is.null(susie_fit)) {
    susie_fit <- .set_finemapping_fit_class(susie_fit, "susie")
  }
  if (!is.null(susie_inf_fit)) {
    susie_inf_fit <- .set_finemapping_fit_class(susie_inf_fit, "susie_inf")
  }

  if (has_susie && has_susie_inf && ncol(Y) == 1 &&
      is.null(susie_fit) && is.null(susie_inf_fit)) {
    fit_arg_names <- c("susie_fit", "susie_inf_fit", "retain_fit")
    fits <- fit_susie_inf_then_susie(
      X,
      Y[, 1],
      args = weight_methods[["susie_weights"]][setdiff(names(weight_methods[["susie_weights"]]), fit_arg_names)],
      susie_inf_args = modifyList(
        list(convergence_method = "pip"),
        weight_methods[["susie_inf_weights"]][setdiff(names(weight_methods[["susie_inf_weights"]]), fit_arg_names)]
      ),
      fitted_models = list(susie = susie_fit, susie_inf = susie_inf_fit)
    )
    susie_fit <- fits[["susie"]]
    susie_inf_fit <- fits[["susie_inf"]]
  }

  if (!is.null(susie_inf_fit) && has_susie_inf) {
    weight_methods[["susie_inf_weights"]][["susie_inf_fit"]] <- susie_inf_fit
  }
  if (!is.null(susie_fit) && has_susie) {
    weight_methods[["susie_weights"]][["susie_fit"]] <- susie_fit
  }
  if (has_susie &&
      is.null(weight_methods[["susie_weights"]][["susie_fit"]]) &&
      !is.null(susie_inf_fit)) {
    weight_methods[["susie_weights"]] <- prepare_susie_from_inf_args(weight_methods[["susie_weights"]], susie_inf_fit)
  }
  weight_methods
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
#' If NULL, 'sample_partitions' must be provided.
#' @param sample_partitions An optional dataframe with predefined sample partitions,
#' containing columns 'Sample' (sample names) and 'Fold' (fold number). If NULL, 'fold' must be provided.
#' @param weight_methods A list of methods and their specific arguments, formatted as list(method1 = method1_args, method2 = method2_args), or alternatively a character vector of method names (eg, c("susie_weights", "enet_weights")) in which case default arguments will be used for all methods.
#' methods in the list can be either univariate (applied to each column of Y) or multivariate (applied to the entire Y matrix).
#' @param max_num_variants An optional integer to set the randomly selected maximum number of variants to use for CV purpose, to save computing time.
#' @param variants_to_keep An optional integer to ensure that the listed variants are kept in the CV when there is a limit on the max_num_variants to use.
#' @param num_threads The number of threads to use for parallel processing.
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
twas_weights_cv <- function(X, Y, fold = NULL, sample_partitions = NULL, weight_methods = NULL, max_num_variants = NULL, variants_to_keep = NULL, num_threads = 1, verbose = 1, ...) {
  split_data <- function(X, Y, sample_partition, fold) {
    test_ids <- sample_partition[which(sample_partition$Fold == fold), "Sample"]
    Xtrain <- X[!(rownames(X) %in% test_ids), , drop = FALSE]
    Ytrain <- Y[!(rownames(Y) %in% test_ids), , drop = FALSE]
    Xtest <- X[rownames(X) %in% test_ids, , drop = FALSE]
    Ytest <- Y[rownames(Y) %in% test_ids, , drop = FALSE]
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
    sample_names <- rownames(Y)
  } else if (!is.null(rownames(Y))) {
    sample_names <- rownames(Y)
  } else if (!is.null(rownames(X))) {
    sample_names <- rownames(X)
  } else {
    sample_names <- paste0("sample_", 1:nrow(X))
  }
  if (is.null(rownames(X))) {
    rownames(X) <- sample_names
  }
  if (is.null(rownames(Y))) {
    rownames(Y) <- sample_names
  }

  if (is.null(colnames(X))) {
    colnames(X) <- paste0("variable_", 1:ncol(X))
  }
  if (is.null(colnames(Y))) {
    colnames(Y) <- paste0("context_", 1:ncol(Y))
  }

  if (is.character(weight_methods)) {
    weight_methods <- .twas_method_lookup(weight_methods)
  }

  if (!exists(".Random.seed")) {
    if (verbose >= 1) message("! No seed has been set. Please set seed for reproducable result. ")
  }

  # Select variants if necessary
  if (!is.null(max_num_variants) && ncol(X) > max_num_variants) {
    if (!is.null(variants_to_keep) && length(variants_to_keep) > 0) {
      variants_to_keep <- intersect(variants_to_keep, colnames(X))
      remaining_columns <- setdiff(colnames(X), variants_to_keep)
      if (length(variants_to_keep) < max_num_variants) {
        additional_columns <- sample(remaining_columns, max_num_variants - length(variants_to_keep), replace = FALSE)
        selected_columns <- union(variants_to_keep, additional_columns)
        if (verbose >= 1) message(sprintf(
          "Including %d specified variants and randomly selecting %d additional variants, for a total of %d variants out of %d for cross-validation purpose.",
          length(variants_to_keep), length(additional_columns), length(selected_columns), ncol(X)
        ))
      } else {
        selected_columns <- sample(variants_to_keep, max_num_variants, replace = FALSE)
        if (verbose >= 1) message(paste("Randomly selecting", length(selected_columns), "out of", length(variants_to_keep), "input variants for cross validation purpose."))
      }
    } else {
      selected_columns <- sort(sample(ncol(X), max_num_variants, replace = FALSE))
      if (verbose >= 1) message(paste("Randomly selecting", length(selected_columns), "out of", ncol(X), "variants for cross validation purpose."))
    }
    X <- X[, selected_columns, drop = FALSE]
  }

  # Create or use provided folds
  if (!is.null(fold)) {
    if (!is.null(sample_partitions)) {
      if (fold != length(unique(sample_partitions$Fold))) {
        if (verbose >= 1) message(paste0(
          "fold number provided does not match with sample partition, performing ", length(unique(sample_partitions$Fold)),
          " fold cross validation based on provided sample partition. "
        ))
      }

      folds <- sample_partitions$Fold
      sample_partition <- sample_partitions
    } else {
      sample_indices <- sample(nrow(X))
      folds <- cut(seq(1, nrow(X)), breaks = fold, labels = FALSE)
      sample_partition <- data.frame(Sample = sample_names[sample_indices], Fold = folds, stringsAsFactors = FALSE)
    }
  } else if (!is.null(sample_partitions)) {
    if (!all(sample_partitions$Sample %in% sample_names)) {
      stop("Some samples in 'sample_partitions' do not match the samples in 'X' and 'Y'.")
    }
    sample_partition <- sample_partitions
    fold <- length(unique(sample_partition$Fold))
  } else {
    stop("Either 'fold' or 'sample_partitions' must be provided.")
  }

  st <- proc.time()
  if (is.null(weight_methods)) {
    return(list(sample_partition = sample_partition))
  } else {
    # Hardcoded vector of multivariate weight_methods
    multivariate_weight_methods <- c("mrmash_weights", "mvsusie_weights")

    # Determine the number of cores to use
    num_cores <- ifelse(num_threads == -1,
      bpworkers(MulticoreParam()),
      num_threads)
    num_cores <- min(num_cores,
      bpworkers(MulticoreParam()))

    cv_args <- list(...)

    # Perform CV with parallel processing
    compute_method_predictions <- function(j) {
      if (verbose >= 1) {
        message(sprintf("  CV fold %d/%d ...", j, fold))
        tic()
      }
      dat_split <- split_data(X, Y, sample_partition = sample_partition, fold = j)
      X_train <- dat_split$Xtrain
      Y_train <- dat_split$Ytrain
      X_test <- dat_split$Xtest
      Y_test <- dat_split$Ytest

      # Remove columns with zero variance
      valid_columns <- .nonzero_var_columns(X_train)
      X_train <- X_train[, valid_columns, drop = FALSE]
      X_train <- filter_X_with_Y(X_train, Y_train, missing_rate_thresh = 1, maf_thresh = NULL)
      valid_columns <- colnames(X_train)
      # X_test <- X_test[, valid_columns, drop=FALSE]
      fold_weight_methods <- .prepare_susie_weight_methods(X_train, Y_train, weight_methods)

      fold_preds <- setNames(lapply(names(fold_weight_methods), function(method) {
        args <- fold_weight_methods[[method]]

        if (method %in% multivariate_weight_methods) {
          # Apply multivariate method to entire Y for this fold
          if (!is.null(cv_args$data_driven_prior_matrices_cv)) {
            if (method == "mrmash_weights") {
              args$data_driven_prior_matrices <- cv_args$data_driven_prior_matrices_cv[[j]]
            }
            if (method == "mvsusie_weights") {
              args$prior_variance <- cv_args$reweighted_mixture_prior_cv[[j]]
            }
          }
          weights_matrix <- if (verbose < 2) {
            .quiet_eval(do.call(method, c(list(X = X_train, Y = Y_train), args)))
          } else {
            do.call(method, c(list(X = X_train, Y = Y_train), args))
          }
          rownames(weights_matrix) <- colnames(X_train)
          full_weights_matrix <- .embed_weights(weights_matrix[valid_columns, , drop = FALSE], valid_columns, ncol(X), ncol(Y), colnames(X), colnames(Y))
          Y_pred <- X_test %*% full_weights_matrix
          rownames(Y_pred) <- rownames(X_test)
          return(Y_pred)
        } else {
          Y_pred <- sapply(1:ncol(Y_train), function(k) {
            weights <- if (verbose < 2) {
              .quiet_eval(do.call(method, c(list(X = X_train, y = Y_train[, k]), args)))
            } else {
              do.call(method, c(list(X = X_train, y = Y_train[, k]), args))
            }
            full_weights <- rep(0, ncol(X))
            names(full_weights) <- colnames(X)
            full_weights[valid_columns] <- weights
            # Handle NAs in weights
            full_weights[is.na(full_weights)] <- 0
            X_test %*% full_weights
          })
          rownames(Y_pred) <- rownames(X_test)
          return(Y_pred)
        }
      }), names(fold_weight_methods))
      if (verbose >= 1) {
        elapsed <- toc(quiet = TRUE)
        message(sprintf("  CV fold %d/%d done in %.1fs", j, fold, elapsed$toc - elapsed$tic))
      }
      fold_preds
    }

    if (num_cores >= 2) {
      bp_param <- MulticoreParam(workers = num_cores,
                                  RNGseed = 1L)
      fold_results <- bplapply(1:fold,
        compute_method_predictions, BPPARAM = bp_param)
    } else {
      fold_results <- map(1:fold, compute_method_predictions)
    }

    # Reorganize into Y_pred
    # After cross validation, each sample should have been in
    # test set at some point, and therefore has predicted value.
    # The prediction matrix is therefore exactly the same dimension as input Y
    Y_pred <- setNames(lapply(weight_methods, function(x) `dimnames<-`(matrix(NA, nrow(Y), ncol(Y)), dimnames(Y))), names(weight_methods))
    for (j in seq_along(fold_results)) {
      for (method in names(weight_methods)) {
        Y_pred[[method]][rownames(fold_results[[j]][[method]]), ] <- fold_results[[j]][[method]]
      }
    }

    names(Y_pred) <- gsub("_weights", "_predicted", names(Y_pred))

    # Compute rsq, adj rsq, p-value, RMSE, and MAE for each method
    metrics_table <- list()

    for (m in names(weight_methods)) {
      metrics_table[[m]] <- matrix(NA, nrow = ncol(Y), ncol = 6)
      colnames(metrics_table[[m]]) <- c("corr", "rsq", "adj_rsq", "pval", "RMSE", "MAE")
      rownames(metrics_table[[m]]) <- colnames(Y)

      for (r in 1:ncol(Y)) {
        method_predictions <- Y_pred[[gsub("_weights", "_predicted", m)]][, r]
        actual_values <- Y[, r]
        # Remove missing values in the first place
        na_indx <- which(is.na(actual_values))
        if (length(na_indx) != 0) {
          method_predictions <- method_predictions[-na_indx]
          actual_values <- actual_values[-na_indx]
        }
        if (sd(method_predictions) != 0) {
          lm_fit <- lm(actual_values ~ method_predictions)

          # Calculate raw correlation and and adjusted R-squared
          metrics_table[[m]][r, "corr"] <- cor(actual_values, method_predictions)

          metrics_table[[m]][r, "rsq"] <- summary(lm_fit)$r.squared
          metrics_table[[m]][r, "adj_rsq"] <- summary(lm_fit)$adj.r.squared

          # Calculate p-value
          metrics_table[[m]][r, "pval"] <- summary(lm_fit)$coefficients[2, 4]

          # Calculate RMSE
          residuals <- actual_values - method_predictions
          metrics_table[[m]][r, "RMSE"] <- sqrt(mean(residuals^2))

          # Calculate MAE
          metrics_table[[m]][r, "MAE"] <- mean(abs(residuals))
        } else {
          metrics_table[[m]][r, ] <- NA
          if (verbose >= 1) message(paste0(
            "Predicted values for condition ", r, " using ", m,
            " have zero variance. Filling performance metric with NAs"
          ))
        }
      }
    }
    names(metrics_table) <- gsub("_weights", "_performance", names(metrics_table))
    return(list(sample_partition = sample_partition, prediction = Y_pred, performance = metrics_table, time_elapsed = proc.time() - st))
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
#' @param weight_methods A list of methods and their specific arguments, formatted as list(method1 = method1_args, method2 = method2_args), or alternatively a character vector of method names (eg, c("susie_weights", "enet_weights")) in which case default arguments will be used for all methods.
#' methods in the list can be either univariate (applied to each column of Y) or multivariate (applied to the entire Y matrix).
#' @param num_threads The number of threads to use for parallel processing.
#'        If set to -1, the function uses all available cores.
#'        If set to 0 or 1, no parallel processing is performed.
#'        If set to 2 or more, parallel processing is enabled with that many threads.
#' @param fitted_models Optional named list of fitted SuSiE-family models.
#' @param retain_fits If TRUE, retain fitted model objects as attributes on
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
twas_weights <- function(X, Y, weight_methods, num_threads = 1,
                         fitted_models = NULL, retain_fits = FALSE, verbose = 1) {
  if (!is.matrix(X) || (!is.matrix(Y) && !is.vector(Y))) {
    stop("X must be a matrix and Y must be a matrix or a vector.")
  }

  if (is.vector(Y)) {
    Y <- matrix(Y, ncol = 1)
  }

  if (nrow(X) != nrow(Y)) {
    stop("The number of rows in X and Y must be the same.")
  }

  if (is.character(weight_methods)) {
    weight_methods <- .twas_method_lookup(weight_methods)
  }

  # Determine number of cores to use
  num_cores <- ifelse(num_threads == -1,
    bpworkers(MulticoreParam()),
    num_threads)
  num_cores <- min(num_cores,
    bpworkers(MulticoreParam()))

  valid_columns <- .nonzero_var_columns(X)
  X_filtered <- as.matrix(X[, valid_columns, drop = FALSE])
  weight_methods <- .prepare_susie_weight_methods(
    X_filtered, Y, weight_methods, fitted_models
  )

  compute_method_weights <- function(method_name, weight_methods) {
    short_name <- sub("_weights$", "", method_name)
    if (verbose >= 1) {
      message(sprintf("  Fitting %s ...", short_name))
      tic()
    }

    # Hardcoded vector of multivariate methods
    multivariate_weight_methods <- c("mrmash_weights", "mvsusie_weights")
    args <- weight_methods[[method_name]]

    # Only pass retain_fit to functions that accept it
    if (retain_fits && "retain_fit" %in% names(formals(method_name))) {
      args$retain_fit <- TRUE
    }

    method_fit <- NULL
    if (method_name %in% multivariate_weight_methods) {
      # Apply multivariate method
      weights_matrix <- if (verbose < 2) {
        .quiet_eval(do.call(method_name, c(list(X = X_filtered, Y = Y), args)))
      } else {
        do.call(method_name, c(list(X = X_filtered, Y = Y), args))
      }
      if (retain_fits) method_fit <- attr(weights_matrix, "fit")
      if (nrow(weights_matrix) != length(valid_columns)) weights_matrix <- weights_matrix[names(valid_columns), , drop = FALSE]
    } else {
      # Apply univariate method to each column of Y
      # Initialize it with zeros to avoid NA
      weights_matrix <- matrix(0, nrow = ncol(X_filtered), ncol = ncol(Y))

      for (k in 1:ncol(Y)) {
        weights_vector <- if (verbose < 2) {
          .quiet_eval(do.call(method_name, c(list(X = X_filtered, y = Y[, k]), args)))
        } else {
          do.call(method_name, c(list(X = X_filtered, y = Y[, k]), args))
        }
        if (retain_fits && is.null(method_fit)) {
          method_fit <- attr(weights_vector, "fit")
        }
        if (is.matrix(weights_vector)) weights_vector <- weights_vector[, k]
        weights_matrix[, k] <- weights_vector
      }
    }

    result <- .embed_weights(weights_matrix, valid_columns, ncol(X), ncol(Y), colnames(X), colnames(Y))
    if (!is.null(method_fit)) attr(result, "fit") <- method_fit
    if (verbose >= 1) {
      elapsed <- toc(quiet = TRUE)
      message(sprintf("  Fitting %s done in %.1fs", short_name, elapsed$toc - elapsed$tic))
    }
    return(result)
  }

  if (num_cores >= 2) {
    bp_param <- MulticoreParam(workers = num_cores,
                                RNGseed = 1L)
    weights_list <- bplapply(names(weight_methods),
      compute_method_weights, weight_methods, BPPARAM = bp_param)
  } else {
    weights_list <- names(weight_methods) %>% map(compute_method_weights, weight_methods)
  }
  names(weights_list) <- names(weight_methods)

  if (!is.null(colnames(X))) {
    weights_list <- lapply(weights_list, function(x) {
      fit <- attr(x, "fit")
      rownames(x) <- colnames(X)
      if (!is.null(fit)) attr(x, "fit") <- fit
      return(x)
    })
  }
  # Create TWASWeights S4 object
  variant_ids <- if (!is.null(colnames(X))) colnames(X) else paste0("variant_", seq_len(ncol(X)))
  fits_list <- lapply(weights_list, function(w) attr(w, "fit"))
  has_any_fit <- any(!sapply(fits_list, is.null))

  # Strip fit attributes from weight matrices before storing in S4
  clean_weights <- lapply(weights_list, function(w) { attr(w, "fit") <- NULL; w })

  TWASWeights(
    weights = clean_weights,
    variant_ids = variant_ids,
    fits = if (has_any_fit) fits_list else NULL,
    cv_performance = NULL
  )
}

#' Predict outcomes using TWAS weights
#'
#' This function takes a matrix of predictors (\code{X}) and a list of TWAS (transcriptome-wide
#' association studies) weights (\code{weights_list}), and calculates the predicted outcomes by
#' multiplying \code{X} by each set of weights in \code{weights_list}. The names of the elements
#' in the output list are derived from the names in \code{weights_list}, with "_weights" replaced
#' by "_predicted".
#'
#' @param X A matrix or data frame of predictors where each row is an observation and each
#' column is a variable.
#' @param weights_list A list of numeric vectors representing the weights for each predictor.
#' The names of the list elements should follow the pattern \code{[outcome]_weights}, where
#' \code{[outcome]} is the name of the outcome variable that the weights are associated with.
#'
#' @return A named list of numeric vectors, where each vector is the predicted outcome for the
#' corresponding set of weights in \code{weights_list}. The names of the list elements are
#' derived from the names in \code{weights_list} by replacing "_weights" with "_predicted".
#'
#' @export
#' @examples
#' # Assuming `X` is your matrix of predictors and `weights_list` is your list of weights:
#' predicted_outcomes <- twas_predict(X, weights_list)
#' print(predicted_outcomes)
twas_predict <- function(X, weights_list) {
  if (is(weights_list, "TWASWeights")) {
    wl <- getWeights(weights_list)
  } else {
    wl <- weights_list
  }
  setNames(lapply(wl, function(w) X %*% w), gsub("_weights", "_predicted", names(wl)))
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
#' @param weight_results Named list of weight vectors or matrices as
#'   returned by \code{\link{twas_weights}}. The mr.ash element should
#'   have a \code{"fit"} attribute containing the model fit object
#'   (set \code{retain_fits = TRUE} in \code{twas_weights} to obtain this).
#'
#' @return A scalar sparsity estimate (proportion of non-zero effects).
#' @export
estimate_sparsity <- function(weight_results) {
  if (is(weight_results, "TWASWeights")) {
    fit <- getFits(weight_results, "mrash_weights")
    if (is.null(fit) || is.null(fit$pi)) {
      stop("mr.ash fit object not found. Run twas_weights() with retain_fits = TRUE ",
           "and ensure mrash_weights is included.")
    }
  } else {
    w <- weight_results[["mrash_weights"]]
    if (is.null(w)) {
      stop("mr.ash weights ('mrash_weights') not found in weight_results.")
    }
    fit <- attr(w, "fit")
    if (is.null(fit) || is.null(fit$pi)) {
      stop("mr.ash fit object not found. Run twas_weights() with retain_fits = TRUE ",
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
#' @param susie_fit An object returned by the SuSiE function, containing the SuSiE model fit.
#' @param fitted_models Optional named list of fitted fine-mapping models, such
#'   as \code{list(susie = susie_fit, susie_inf = susie_inf_fit)}.
#' @param cv_folds The number of folds to use for cross-validation. Set to 0 to skip cross-validation. Defaults to 5.
#' @param sample_partition Optional data frame with Sample and Fold columns for cross-validation. If NULL, a random partition is generated.
#' @param weight_methods List of methods to use to compute weights for TWAS; along with their parameters.
#' @param max_cv_variants The maximum number of variants to be included in cross-validation. Defaults to -1 which means no limit.
#' @param cv_threads The number of threads to use for parallel computation in cross-validation. Defaults to 1.
#' @param cv_weight_methods List of methods to use for cross-validation. If NULL, uses the same methods as weight_methods.
#' @param ensemble Logical. If TRUE and cv_folds > 1, learn ensemble combination
#'   weights via stacked regression (SR-TWAS). Requires at least two individual
#'   methods to have been run and to pass the R-squared cutoff. Defaults to TRUE.
#' @param ensemble_r2_threshold Minimum cross-validated R-squared for an individual method
#'   to be included in the ensemble. Methods below this threshold are excluded.
#'   Defaults to 0.01.
#' @param ensemble_solver Character string specifying the optimization backend
#'   for ensemble learning. One of \code{"quadprog"}, \code{"nnls"},
#'   \code{"lbfgsb"}, or \code{"glmnet"}. Passed to
#'   \code{\link{ensemble_weights}}. Defaults to \code{"quadprog"}.
#' @param ensemble_alpha Elastic net mixing parameter, used only when
#'   \code{ensemble_solver = "glmnet"}. Defaults to 1 (lasso).
#' @param estimate_pi If TRUE, estimate spike-and-slab sparsity from mr.ash
#'   before running Bayesian alphabet methods that need inclusion probabilities.
#' @param verbose Integer controlling verbosity level: 0 = suppress all messages,
#'   1 = show pecotmr messages but suppress external package messages (default),
#'   2 = show all messages including those from external packages.
#'
#' @return A list containing results from the TWAS pipeline, including TWAS weights, predictions, and optionally cross-validation results.
#' @export
#'
#' @examples
#' # Example usage (assuming appropriate objects for X, y, and susie_fit are available):
#' twas_results <- twas_weights_pipeline(X, y, susie_fit)
twas_weights_pipeline <- function(X,
                                  y,
                                  susie_fit = NULL,
                                  fitted_models = NULL,
                                  cv_folds = 5,
                                  sample_partition = NULL,
                                  weight_methods = "default",
                                  max_cv_variants = -1,
                                  cv_threads = 1,
                                  cv_weight_methods = NULL,
                                  ensemble = TRUE,
                                  ensemble_r2_threshold = 0.01,
                                  ensemble_solver = "quadprog",
                                  ensemble_alpha = 1,
                                  estimate_pi = TRUE,
                                  verbose = 1) {
  if (is.character(weight_methods)) {
    weight_methods <- .twas_method_lookup(weight_methods)
  }
  if (is.null(fitted_models)) fitted_models <- list()
  if (!is.null(susie_fit)) fitted_models[["susie"]] <- susie_fit

  res <- list()
  st <- proc.time()
  if (verbose >= 1) {
    message("Performing TWAS weights computation for univariate analysis methods ...")
    tic()
  }

  if (!is.null(fitted_models[["susie"]]) && !is.null(weight_methods$susie_weights)) {
    res$susie_weights_intermediate <- .susie_weight_intermediate(fitted_models[["susie"]], X)
  }

  # Check if empirical pi estimation is needed for spike-and-slab methods
  bayes_c_needs_pi <- "bayes_c_weights" %in% names(weight_methods) &&
    !"pi" %in% names(weight_methods$bayes_c_weights)
  bayes_b_needs_pi <- "bayes_b_weights" %in% names(weight_methods) &&
    !"probIn" %in% names(weight_methods$bayes_b_weights)
  needs_pi_estimation <- (bayes_c_needs_pi || bayes_b_needs_pi) && estimate_pi

  if (needs_pi_estimation) {
    # Run mr.ash first to estimate sparsity
    mrash_methods <- list(mrash_weights = weight_methods[["mrash_weights"]] %||% list())

    if (verbose >= 1) message("  Estimating sparsity from mr.ash ...")
    mrash_weights <- twas_weights(X, y, weight_methods = mrash_methods, retain_fits = TRUE, verbose = verbose)

    empirical_pi <- estimate_sparsity(mrash_weights)
    if (verbose >= 1) message(sprintf("  Empirical sparsity estimate: %.4f", empirical_pi))
    res$empirical_pi <- empirical_pi

    # Inject into spike-and-slab methods that need it
    if (bayes_c_needs_pi) weight_methods$bayes_c_weights$pi <- as.numeric(empirical_pi)
    if (bayes_b_needs_pi) weight_methods$bayes_b_weights$probIn <- as.numeric(empirical_pi)

    # Run remaining methods (those not already computed)
    remaining_fn_names <- setdiff(names(weight_methods), "mrash_weights")

    if (length(remaining_fn_names) > 0) {
      remaining_methods <- weight_methods[remaining_fn_names]
      remaining_tw <- twas_weights(
        X,
        y,
        weight_methods = remaining_methods,
        fitted_models = fitted_models,
        verbose = verbose
      )
      # Combine two TWASWeights objects
      combined_weights <- c(getWeights(mrash_weights), getWeights(remaining_tw))
      combined_fits <- c(getFits(mrash_weights), getFits(remaining_tw))
      res$twas_weights <- TWASWeights(
        weights = combined_weights,
        variant_ids = getVariantIds(mrash_weights),
        fits = combined_fits
      )
    } else {
      res$twas_weights <- mrash_weights
    }

    # Remove mr.ash if it was not in the original weight_methods
    if (!"mrash_weights" %in% names(weight_methods)) {
      w_list <- getWeights(res$twas_weights)
      f_list <- getFits(res$twas_weights)
      w_list[["mrash_weights"]] <- NULL
      if (!is.null(f_list)) f_list[["mrash_weights"]] <- NULL
      res$twas_weights <- TWASWeights(
        weights = w_list,
        variant_ids = getVariantIds(res$twas_weights),
        fits = if (length(f_list) > 0) f_list else NULL
      )
    }
  } else {
    # Run all methods at once
    res$twas_weights <- twas_weights(
      X,
      y,
      weight_methods = weight_methods,
      fitted_models = fitted_models,
      verbose = verbose
    )
  }
  if (verbose >= 1) {
    elapsed <- toc(quiet = TRUE)
    message(sprintf("TWAS weights fitting done in %.1fs", elapsed$toc - elapsed$tic))
  }
  res$twas_predictions <- twas_predict(X, res$twas_weights)

  if (cv_folds > 1) {
    # A few cutting corners to run CV faster at the disadvantage of SuSiE and mr.ash:
    # 1. reset SuSiE to not using refine or adaptive L but to use L from previous analysis
    # 2. at most 100 iterations for mr.ash allowed
    # 3. only use a subset of variants randomly selected to avoid bias
    if (!is.null(fitted_models[["susie_inf"]]) && !is.null(weight_methods$susie_inf_weights)) {
      weight_methods$susie_inf_weights$L <- length(fitted_models[["susie_inf"]]$V)
      weight_methods$susie_inf_weights$refine <- FALSE
    }
    if (!is.null(weight_methods$susie_weights)) {
      susie_cv_fit <- fitted_models[["susie"]]
      if (is.null(susie_cv_fit)) susie_cv_fit <- fitted_models[["susie_inf"]]
      if (!is.null(susie_cv_fit)) {
        weight_methods$susie_weights$L <- length(susie_cv_fit$V)
        weight_methods$susie_weights$refine <- FALSE
      }
    }
    if (is.null(cv_weight_methods)) {
      cv_weight_methods <- .filter_zero_weight_methods(weight_methods, res$twas_weights)
    }

    variants_for_cv <- c()
    if (max_cv_variants <= 0) {
      max_cv_variants <- Inf
    }
    if (ncol(X) > max_cv_variants) {
      variants_for_cv <- sample(colnames(X), max_cv_variants, replace = FALSE)
    }

    if (verbose >= 1) {
      message("Performing cross-validation to assess TWAS weights ...")
      tic()
    }
    res$twas_cv_result <- twas_weights_cv(
      X,
      y,
      fold = cv_folds,
      sample_partitions = sample_partition,
      weight_methods = cv_weight_methods,
      max_num_variants = max_cv_variants,
      num_threads = cv_threads,
      verbose = verbose,
      variants_to_keep = if (length(variants_for_cv) > 0) variants_for_cv else NULL
    )
    if (verbose >= 1) {
      elapsed <- toc(quiet = TRUE)
      message(sprintf("Cross-validation done in %.1fs", elapsed$toc - elapsed$tic))
    }

    # Ensemble learning: learn optimal method combination via stacked regression
    if (isTRUE(ensemble) && length(cv_weight_methods) <= 1) {
      if (verbose >= 1) message("Ensemble model skipped: only ", length(cv_weight_methods),
              " weight method provided (need >= 2 for ensemble learning).")
    }
    if (isTRUE(ensemble) && length(cv_weight_methods) > 1) {
      if (!is.null(res$twas_cv_result$performance)) {
        # Extract R-squared for each method from CV performance table
        method_rsq <- vapply(res$twas_cv_result$performance, function(perf) {
          perf[1, "rsq"]
        }, numeric(1))
        names(method_rsq) <- gsub("_performance$", "", names(method_rsq))

        # NA R-squared already implies the method is unusable for the ensemble: a
        # method whose CV predictions are degenerate (zero variance across all
        # held-out folds) yields cor(predictions, y) = NA and therefore rsq = NA.
        # So !is.na(method_rsq) is sufficient to drop both NA-rsq and degenerate
        # methods - no separate variance check needed.
        passing <- !is.na(method_rsq) & method_rsq >= ensemble_r2_threshold
        n_passing <- sum(passing)

        if (n_passing < 2) {
          # Ensemble (stacked regression) requires at least 2 base learners.
          # Build a per-method status line so the user can see which methods
          # dropped out and why (NA R-squared from degenerate CV predictions,
          # or simply R-squared below the cutoff).
          reason <- ifelse(passing, "(passed)",
                    ifelse(is.na(method_rsq),
                           "(dropped: NA R-squared - likely degenerate CV predictions)",
                           "(dropped: R-squared below cutoff)"))
          passed_info <- paste0("  ", names(method_rsq), ": R-squared = ",
                                round(method_rsq, 4), " ", reason)
          surviving <- if (n_passing == 1) {
            paste0(" Use the surviving method's weights directly: ",
                   names(method_rsq)[passing], ".")
          } else ""
          if (verbose >= 1) message("Ensemble TWAS skipped: ", n_passing, " of ", length(method_rsq),
                  " methods passed the R-squared cutoff of ", ensemble_r2_threshold,
                  " (need >= 2).", surviving, "\n",
                  "Method R-squared values:\n",
                  paste(passed_info, collapse = "\n"))
        } else {
          passing_base <- names(method_rsq)[passing]
          passing_pred_names <- paste0(passing_base, "_predicted")
          passing_weight_names <- paste0(passing_base, "_weights")

          # Subset cv_results predictions to passing methods
          filtered_cv <- res$twas_cv_result
          filtered_cv$prediction <- filtered_cv$prediction[passing_pred_names]

          # Subset twas_weights to passing methods
          if (is(res$twas_weights, "TWASWeights")) {
            wl <- getWeights(res$twas_weights)
            filtered_weights <- wl[passing_weight_names]
          } else {
            filtered_weights <- res$twas_weights[passing_weight_names]
          }

          if (verbose >= 1) {
            message("Computing ensemble TWAS weights via stacked regression ",
                    "using ", n_passing, " methods: ",
                    paste(passing_base, collapse = ", "), " ...")
            tic()
          }
          ens_result <- ensemble_weights(
            cv_results = filtered_cv,
            Y = y,
            twas_weight_list = filtered_weights,
            solver = ensemble_solver,
            alpha = ensemble_alpha
          )
          if (verbose >= 1) {
            elapsed <- toc(quiet = TRUE)
            message(sprintf("Ensemble learning done in %.1fs", elapsed$toc - elapsed$tic))
          }

          # Add ensemble weights alongside individual method weights
          if (!is.null(ens_result$ensemble_twas_weights)) {
            ens_wt <- ens_result$ensemble_twas_weights
            if (!is.matrix(ens_wt)) ens_wt <- matrix(ens_wt, ncol = 1)
            # Rebuild TWASWeights S4 with ensemble method added
            tw <- res$twas_weights
            new_weights <- c(getWeights(tw), list(ensemble_weights = ens_wt))
            res$twas_weights <- new("TWASWeights",
              weights = new_weights,
              variant_ids = getVariantIds(tw),
              methods = c(getMethodNames(tw), "ensemble_weights"),
              fits = getFits(tw),
              cv_performance = getCVPerformance(tw),
              standardized = getStandardized(tw)
            )
            res$twas_predictions$ensemble_predicted <- X %*% ens_wt
          }
          res$ensemble <- ens_result
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
#' @param mnm_fit An object containing the fitted multivariate models (e.g., mvSuSiE and mr.mash fits).
#' @param L Maximum number of components in mvSuSiE. If NULL, the number of
#'   components in the fitted mvSuSiE object is used.
#' @param L_greedy Initial greedy number of components in mvSuSiE. Defaults to 5.
#' @param cv_folds The number of folds to use for cross-validation. Defaults to 5. Set to 0 to skip cross-validation.
#' @param sample_partition Optional data frame with Sample and Fold columns for cross-validation. If NULL, a random partition is generated.
#' @param data_driven_prior_matrices A list of data-driven covariance matrices for mr.mash weights. Defaults to NULL.
#' @param data_driven_prior_matrices_cv A list of data-driven covariance matrices for mr.mash weights in cross-validation. Defaults to NULL.
#' @param canonical_prior_matrices If TRUE, computes canonical covariance matrices for mr.mash. Defaults to FALSE.
#' @param mvsusie_max_iter The maximum number of iterations for mvSuSiE. Defaults to 200.
#' @param mrmash_max_iter The maximum number of iterations for mr.mash. Defaults to 5000.
#' @param max_cv_variants The maximum number of variants to be included in cross-validation. Defaults to -1 which means no limit.
#' @param cv_threads The number of threads to use for parallel computation in cross-validation. Defaults to 1.
#' @param verbose Integer controlling verbosity level: 0 = suppress all messages,
#'   1 = show pecotmr messages but suppress external package messages (default),
#'   2 = show all messages including those from external packages.
#'
#' @return A list containing results from the TWAS pipeline, including TWAS weights, predictions, and optionally cross-validation results.
#' @export
#' @examples
#' # Example usage (assuming appropriate objects for X, Y, and mnm_fit are available):
#' twas_results <- twas_multivariate_weights_pipeline(X, Y, mnm_fit)
twas_multivariate_weights_pipeline <- function(
    X,
    Y,
    mnm_fit,
    L = NULL,
    L_greedy = 5,
    cv_folds = 5,
    sample_partition = NULL,
    data_driven_prior_matrices = NULL,
    data_driven_prior_matrices_cv = NULL,
    canonical_prior_matrices = FALSE,
    mvsusie_max_iter = 200,
    mrmash_max_iter = 5000,
    max_cv_variants = -1,
    cv_threads = 1,
    verbose = 1) {
  copy_twas_results <- function(context_names, variant_names, twas_weight, twas_predictions) {
    wl <- if (is(twas_weight, "TWASWeights")) getWeights(twas_weight) else twas_weight
    setNames(lapply(context_names, function(ctx) {
      if (ctx %in% colnames(wl[[1]])) {
        list(
          twas_weights = lapply(wl, function(wgts) wgts[, ctx]),
          twas_predictions = lapply(twas_predictions, function(pred) pred[, ctx]),
          variant_names = variant_names
        )
      } else {
        NULL
      }
    }), context_names)
  }

  copy_twas_cv_results <- function(twas_result, twas_cv_result) {
    for (i in names(twas_result)) {
      if (i %in% colnames(twas_cv_result$prediction[[1]])) {
        twas_result[[i]]$twas_cv_result$sample_partition <- twas_cv_result$sample_partition
        twas_result[[i]]$twas_cv_result$prediction <- lapply(
          twas_cv_result$prediction,
          function(predicted) {
            as.matrix(predicted[, i], ncol = 1)
          }
        )
        twas_result[[i]]$twas_cv_result$performance <- lapply(
          twas_cv_result$performance,
          function(perform) {
            t(as.matrix(perform[i, ], ncol = 1))
          }
        )
        twas_result[[i]]$twas_cv_result$time_elapsed <- twas_cv_result$time_elapsed
      }
    }
    return(twas_result)
  }

  # TWAS weights and predictions
  weight_methods <- list(
    mrmash_weights = list(
      mrmash_fit = mnm_fit$mrmash_fitted
    ),
    mvsusie_weights = list(
      mvsusie_fit = mnm_fit$mvsusie_fitted
    )
  )
  st <- proc.time()
  if (verbose >= 1) {
    message("Extracting TWAS weights for multivariate analysis methods ...")
    tic()
  }
  # get TWAS weights
  twas_weights_res <- twas_weights(X = X, Y = Y, weight_methods = weight_methods, verbose = verbose)
  if (verbose >= 1) {
    elapsed <- toc(quiet = TRUE)
    message(sprintf("Multivariate TWAS weights fitting done in %.1fs", elapsed$toc - elapsed$tic))
  }
  # get TWAS predictions for possible next steps such as computing correlations between predicted expression values
  twas_predictions <- twas_predict(X, twas_weights_res)

  # copy TWAS results by condition
  res <- copy_twas_results(colnames(Y), mnm_fit$variant_names, twas_weights_res, twas_predictions)

  # Perform cross-validation if specified
  if (cv_folds > 1) {
    if (is.null(L)) L <- length(mnm_fit$mvsusie_fitted$V)
    if (!is.null(L_greedy)) L_greedy <- min(L_greedy, L)
    sub_verbose <- verbose >= 2
    weight_methods <- list(
      mrmash_weights = list(
        data_driven_prior_matrices = data_driven_prior_matrices,
        canonical_prior_matrices = canonical_prior_matrices,
        max_iter = mrmash_max_iter,
        verbose = sub_verbose
      ),
      mvsusie_weights = list(
        prior_variance = mnm_fit$reweighted_mixture_prior,
        residual_variance = mnm_fit$mrmash_fitted$V,
        L = L,
        L_greedy = L_greedy,
        max_iter = mvsusie_max_iter,
        verbose = sub_verbose
      )
    )

    weight_methods <- .filter_zero_weight_methods(weight_methods, twas_weights_res)

    variants_for_cv <- c()
    if (max_cv_variants <= 0) max_cv_variants <- Inf
    if (ncol(X) > max_cv_variants) {
      variants_for_cv <- sample(colnames(X), max_cv_variants, replace = FALSE)
    }
    if (verbose >= 1) {
      message("Performing cross-validation to assess TWAS weights ...")
      tic()
    }
    twas_cv_result <- twas_weights_cv(
      X = X, Y = Y, fold = cv_folds,
      weight_methods = weight_methods,
      sample_partitions = sample_partition,
      num_threads = cv_threads,
      max_num_variants = max_cv_variants,
      verbose = verbose,
      variants_to_keep = if (length(variants_for_cv) > 0) variants_for_cv else NULL,
      data_driven_prior_matrices_cv = data_driven_prior_matrices_cv,
      reweighted_mixture_prior_cv = mnm_fit$reweighted_mixture_prior_cv
    )
    if (verbose >= 1) {
      elapsed <- toc(quiet = TRUE)
      message(sprintf("Cross-validation done in %.1fs", elapsed$toc - elapsed$tic))
    }
    res <- copy_twas_cv_results(res, twas_cv_result)
  }
  total_time_elapsed <- proc.time() - st
  for (i in seq_along(res)) {
    res[[i]]$total_time_elapsed <- total_time_elapsed
  }
  return(res)
}


# Solve ensemble stacking via quadprog (constrained QP with sum-to-1 and non-negativity).
# @param P_valid Matrix of CV predictions for valid methods (n x K_valid).
# @param y_obs Observed outcome vector (n).
# @param K_valid Number of valid methods.
# @return Normalized coefficient vector of length K_valid.
# @noRd
.solve_ensemble_quadprog <- function(P_valid, y_obs, K_valid) {
  if (!requireNamespace("quadprog", quietly = TRUE)) {
    stop("Package 'quadprog' is required for solver='quadprog'. ",
         "Install with: install.packages('quadprog')")
  }

  Dmat <- crossprod(P_valid)
  dvec <- as.vector(crossprod(P_valid, y_obs))
  # Ridge term for numerical stability (small relative to trace)
  Dmat <- Dmat + 1e-8 * mean(diag(Dmat)) * diag(K_valid)

  # Constraint matrix: first constraint is equality (sum = 1), then K_valid
  # non-negativity constraints.
  Amat <- cbind(rep(1, K_valid), diag(K_valid))
  bvec <- c(1, rep(0, K_valid))

  qp_sol <- tryCatch(
    solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1),
    error = function(e) {
      warning("QP solver failed: ", conditionMessage(e),
              ". Falling back to equal weights among valid methods.")
      NULL
    }
  )

  if (is.null(qp_sol)) {
    return(rep(1 / K_valid, K_valid))
  }

  # Numerical cleanup: clamp to non-negative and renormalize
  zeta_valid <- pmax(qp_sol$solution, 0)
  zeta_sum <- sum(zeta_valid)
  if (zeta_sum <= 0) {
    warning("QP returned all-zero solution. Falling back to equal weights.")
    return(rep(1 / K_valid, K_valid))
  }
  zeta_valid / zeta_sum
}

# Solve ensemble stacking via NNLS (non-negative least squares, then normalize).
# This is the approach used by SuperLearner (Lawson-Hanson algorithm).
# @param P_valid Matrix of CV predictions for valid methods (n x K_valid).
# @param y_obs Observed outcome vector (n).
# @param K_valid Number of valid methods.
# @return Normalized coefficient vector of length K_valid.
# @noRd
.solve_ensemble_nnls <- function(P_valid, y_obs, K_valid) {
  if (!requireNamespace("nnls", quietly = TRUE)) {
    stop("Package 'nnls' is required for solver='nnls'. ",
         "Install with: install.packages('nnls')")
  }

  fit <- tryCatch(
    nnls::nnls(P_valid, y_obs),
    error = function(e) {
      warning("NNLS solver failed: ", conditionMessage(e),
              ". Falling back to equal weights.")
      NULL
    }
  )

  if (is.null(fit)) {
    return(rep(1 / K_valid, K_valid))
  }

  zeta_valid <- fit$x
  zeta_sum <- sum(zeta_valid)
  if (zeta_sum <= 0) {
    warning("NNLS returned all-zero solution. Falling back to equal weights.")
    return(rep(1 / K_valid, K_valid))
  }
  zeta_valid / zeta_sum
}

# Solve ensemble stacking via L-BFGS-B (box-constrained optimization, then normalize).
# Uses base R optim() with analytical gradient. No extra dependencies.
# @param P_valid Matrix of CV predictions for valid methods (n x K_valid).
# @param y_obs Observed outcome vector (n).
# @param K_valid Number of valid methods.
# @return Normalized coefficient vector of length K_valid.
# @noRd
.solve_ensemble_lbfgsb <- function(P_valid, y_obs, K_valid) {
  PtP <- crossprod(P_valid)
  Pty <- as.vector(crossprod(P_valid, y_obs))

  fn <- function(z) sum((y_obs - P_valid %*% z)^2)
  gr <- function(z) as.vector(2 * (PtP %*% z - Pty))

  fit <- tryCatch(
    optim(
      par = rep(1 / K_valid, K_valid),
      fn = fn, gr = gr,
      method = "L-BFGS-B",
      lower = rep(0, K_valid)
    ),
    error = function(e) {
      warning("L-BFGS-B solver failed: ", conditionMessage(e),
              ". Falling back to equal weights.")
      NULL
    }
  )

  if (is.null(fit)) {
    return(rep(1 / K_valid, K_valid))
  }

  zeta_valid <- pmax(fit$par, 0)
  zeta_sum <- sum(zeta_valid)
  if (zeta_sum <= 0) {
    warning("L-BFGS-B returned all-zero solution. Falling back to equal weights.")
    return(rep(1 / K_valid, K_valid))
  }
  zeta_valid / zeta_sum
}

# Solve ensemble stacking via glmnet (penalized regression with non-negativity).
# Uses cv.glmnet for automatic lambda selection. The alpha parameter controls
# the elastic net mixing: alpha=1 is lasso (sparse), alpha=0 is ridge.
# @param P_valid Matrix of CV predictions for valid methods (n x K_valid).
# @param y_obs Observed outcome vector (n).
# @param K_valid Number of valid methods.
# @param alpha Elastic net mixing parameter (default 1 = lasso).
# @return Normalized coefficient vector of length K_valid.
# @noRd
.solve_ensemble_glmnet <- function(P_valid, y_obs, K_valid, alpha = 1) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for solver='glmnet'. ",
         "Install with: install.packages('glmnet')")
  }

  fit <- tryCatch(
    glmnet::cv.glmnet(
      x = P_valid, y = y_obs,
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
    return(rep(1 / K_valid, K_valid))
  }

  zeta_valid <- as.numeric(coef(fit, s = "lambda.min"))[-1]  # drop intercept
  zeta_valid <- pmax(zeta_valid, 0)
  zeta_sum <- sum(zeta_valid)
  if (zeta_sum <= 0) {
    warning("glmnet returned all-zero solution. Falling back to equal weights.")
    return(rep(1 / K_valid, K_valid))
  }
  zeta_valid / zeta_sum
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
#' For single-dataset usage, pass one \code{twas_weights_cv()} result directly.
#' For multi-dataset ensemble (e.g., combining cell types or reference panels
#' such as CUMC1 + MIT), pass a list of \code{twas_weights_cv()} results along
#' with a list of observed Y vectors - this learns a single joint set of
#' coefficients.
#'
#' @param cv_results Output of \code{\link{twas_weights_cv}}, with \code{$prediction}
#'   (named list of method -> out-of-fold prediction matrix, keys like
#'   \code{"susie_predicted"}). For multi-dataset: a list of such objects.
#' @param Y Observed outcome vector or matrix (samples x contexts). For
#'   multi-dataset: a list of vectors/matrices, one per dataset.
#' @param twas_weight_list Optional named list of weight matrices from
#'   \code{\link{twas_weights}}, with keys like \code{"susie_weights"}. Used to
#'   construct the final combined TWAS weight vector. For multi-dataset: a list
#'   of such lists (the first is used as the weight template).
#' @param context_index Integer indicating which column of Y to use when Y is a
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
#'     \eqn{w = \sum_k \zeta_k w_k}, or NULL if \code{twas_weight_list}
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
#' @seealso \code{\link{twas_weights_cv}}, \code{\link{twas_weights}},
#'   \code{\link{twas_weights_pipeline}}
#'
#' @examples
#' \dontrun{
#' # After running twas_weights_pipeline with CV:
#' res <- twas_weights_pipeline(X, y, cv_folds = 5, weight_methods = methods)
#'
#' ens <- ensemble_weights(
#'   cv_results = res$twas_cv_result,
#'   Y = y,
#'   twas_weight_list = res$twas_weights
#' )
#' ens$method_coef           # combination weights, sum to 1
#'
#' # Multi-dataset ensemble (e.g., CUMC1 + MIT cell types):
#' ens_multi <- ensemble_weights(
#'   cv_results = list(res_cumc$twas_cv_result, res_mit$twas_cv_result),
#'   Y = list(y_cumc, y_mit),
#'   twas_weight_list = list(res_cumc$twas_weights, res_mit$twas_weights)
#' )
#' }
#'
#' @importFrom stats optim coef complete.cases sd cor
#' @export
ensemble_weights <- function(cv_results, Y, twas_weight_list = NULL,
                             context_index = 1,
                             solver = c("quadprog", "nnls", "lbfgsb", "glmnet"),
                             alpha = 1) {
  # --- Input validation ---
  solver <- match.arg(solver)
  if (is.null(cv_results)) {
    stop("'cv_results' is required.")
  }
  if (is.null(Y)) {
    stop("'Y' is required.")
  }
  if (!is.numeric(context_index) || length(context_index) != 1 || context_index < 1) {
    stop("'context_index' must be a positive integer scalar.")
  }

  # --- Normalize single vs multi-dataset input ---
  # Single dataset: cv_results has $prediction directly (is a twas_weights_cv() output).
  # Multi-dataset: cv_results is a list of such outputs.
  is_single <- !is.null(cv_results$prediction)
  if (is_single) {
    cv_results <- list(cv_results)
    Y <- list(Y)
    if (!is.null(twas_weight_list)) twas_weight_list <- list(twas_weight_list)
  } else {
    # Multi-dataset: validate list consistency
    if (!is.list(cv_results) || length(cv_results) == 0) {
      stop("For multi-dataset ensemble, 'cv_results' must be a non-empty list of ",
           "twas_weights_cv() outputs.")
    }
    if (!is.list(Y) || length(Y) != length(cv_results)) {
      stop("'Y' must be a list of the same length as 'cv_results' for ",
           "multi-dataset ensemble.")
    }
    if (!is.null(twas_weight_list)) {
      if (!is.list(twas_weight_list) || length(twas_weight_list) != length(cv_results)) {
        stop("'twas_weight_list' must be a list of the same length as 'cv_results'.")
      }
    }
    for (d in seq_along(cv_results)) {
      if (is.null(cv_results[[d]]$prediction)) {
        stop("cv_results[[", d, "]] does not contain '$prediction'. ",
             "Expected a twas_weights_cv() output.")
      }
    }
  }

  # --- Extract and validate method names ---
  pred_names <- names(cv_results[[1]]$prediction)
  if (is.null(pred_names) || any(pred_names == "")) {
    stop("cv_results$prediction must be a named list (output of twas_weights_cv).")
  }
  base_names <- gsub("_predicted$", "", pred_names)
  K <- length(base_names)

  if (K < 2) {
    stop("Ensemble learning requires at least 2 methods. Found: ", K, ".")
  }

  # Consistency: all datasets must report the same methods in the same order
  for (d in seq_along(cv_results)) {
    if (!identical(names(cv_results[[d]]$prediction), pred_names)) {
      stop("All cv_results must have the same method names (in $prediction) ",
           "in the same order. Dataset 1 has: ", paste(pred_names, collapse = ", "),
           "; dataset ", d, " has: ",
           paste(names(cv_results[[d]]$prediction), collapse = ", "))
    }
  }

  # --- Build stacked prediction matrix P and observed y vector ---
  pred_list <- list()
  y_list <- list()

  for (d in seq_along(cv_results)) {
    preds_d <- cv_results[[d]]$prediction
    y_raw <- Y[[d]]

    # Get sample names from predictions and Y for alignment
    pred_samples <- rownames(preds_d[[pred_names[1]]])
    y_names <- if (is.matrix(y_raw) || is.data.frame(y_raw)) {
      rownames(y_raw)
    } else {
      names(y_raw)
    }

    # Determine sample alignment
    if (!is.null(pred_samples) && !is.null(y_names)) {
      common <- intersect(pred_samples, y_names)
      if (length(common) == 0) {
        stop("No common sample names between predictions and Y in dataset ", d, ".")
      }
      if (length(common) < length(pred_samples) || length(common) < length(y_names)) {
        message("Dataset ", d, ": using ", length(common), " common samples ",
                "(predictions: ", length(pred_samples), ", Y: ", length(y_names), ").")
      }
      # Extract y aligned to common samples
      y_d <- if (is.matrix(y_raw) || is.data.frame(y_raw)) {
        if (context_index > ncol(y_raw)) {
          stop("context_index (", context_index, ") exceeds number of columns in Y[[",
               d, "]] (", ncol(y_raw), ").")
        }
        as.numeric(as.matrix(y_raw)[match(common, y_names), context_index])
      } else {
        as.numeric(y_raw[match(common, y_names)])
      }
      pred_order <- match(common, pred_samples)
      n_d <- length(common)
    } else {
      # No sample names available: fall back to positional alignment
      y_d <- if (is.matrix(y_raw) || is.data.frame(y_raw)) {
        if (context_index > ncol(y_raw)) {
          stop("context_index (", context_index, ") exceeds number of columns in Y[[",
               d, "]] (", ncol(y_raw), ").")
        }
        as.numeric(as.matrix(y_raw)[, context_index])
      } else {
        as.numeric(y_raw)
      }
      n_d <- length(y_d)
      pred_order <- seq_len(n_d)
    }

    P_d <- matrix(NA_real_, nrow = n_d, ncol = K)
    colnames(P_d) <- base_names
    for (k in seq_along(pred_names)) {
      pred_mat <- preds_d[[pred_names[k]]]
      p_col <- if (is.matrix(pred_mat)) pred_mat[pred_order, context_index] else as.numeric(pred_mat)[pred_order]
      if (length(p_col) != n_d) {
        stop("Prediction length for method '", pred_names[k], "' in dataset ", d,
             " (", length(p_col), ") does not match number of aligned samples (", n_d, ").")
      }
      P_d[, k] <- p_col
    }
    pred_list[[d]] <- P_d
    y_list[[d]] <- y_d
  }

  P <- do.call(rbind, pred_list)   # (n_total x K)
  y_obs <- unlist(y_list)           # (n_total)

  # Remove rows with any NA (in P or y)
  complete <- complete.cases(P, y_obs)
  n_dropped <- sum(!complete)
  if (n_dropped > 0) {
    message("Dropping ", n_dropped, " observation(s) with NA predictions or outcomes.")
  }
  if (sum(complete) < K + 1) {
    stop("Too few complete observations (", sum(complete), ") for ", K,
         " methods. Need at least ", K + 1, ".")
  }
  P <- P[complete, , drop = FALSE]
  y_obs <- y_obs[complete]

  # --- Identify methods with non-zero variance predictions ---
  method_sds <- apply(P, 2, sd)
  valid_methods <- method_sds > .Machine$double.eps
  n_valid <- sum(valid_methods)

  if (n_valid < 1) {
    stop("All methods have zero-variance predictions. Cannot compute ensemble. ",
         "This typically means all methods returned zero weights - check that ",
         "the input data has sufficient signal.")
  }

  # --- Solve for combination coefficients ---
  if (n_valid == 1) {
    # Only one method has signal: assign it full weight
    zeta <- rep(0, K)
    zeta[valid_methods] <- 1
    names(zeta) <- base_names
    message("Only one method ('", base_names[valid_methods],
            "') has non-zero variance predictions. Assigning it full weight.")
  } else {
    P_valid <- P[, valid_methods, drop = FALSE]
    K_valid <- ncol(P_valid)

    zeta_valid <- switch(solver,
      quadprog = .solve_ensemble_quadprog(P_valid, y_obs, K_valid),
      nnls     = .solve_ensemble_nnls(P_valid, y_obs, K_valid),
      lbfgsb   = .solve_ensemble_lbfgsb(P_valid, y_obs, K_valid),
      glmnet   = .solve_ensemble_glmnet(P_valid, y_obs, K_valid, alpha = alpha)
    )

    zeta <- rep(0, K)
    zeta[valid_methods] <- zeta_valid
    names(zeta) <- base_names
  }

  # --- Performance metrics ---
  method_rsq <- setNames(vapply(seq_len(K), function(k) {
    if (method_sds[k] > 0) cor(y_obs, P[, k])^2 else NA_real_
  }, numeric(1)), base_names)

  # --- Build ensemble TWAS weight vector (uses first dataset's weights) ---
  ensemble_twas_wt <- NULL
  if (!is.null(twas_weight_list)) {
    wt_list <- twas_weight_list[[1]]
    if (!is.list(wt_list) || length(wt_list) == 0) {
      warning("twas_weight_list[[1]] is empty or not a list; skipping weight combination.")
    } else {
      wt_keys <- paste0(base_names, "_weights")
      matched <- wt_keys %in% names(wt_list)

      if (any(matched)) {
        first_wt <- wt_list[[wt_keys[which(matched)[1]]]]
        if (!is.matrix(first_wt)) first_wt <- matrix(first_wt, ncol = 1)
        p <- nrow(first_wt)
        n_contexts <- ncol(first_wt)

        ensemble_twas_wt <- matrix(0, nrow = p, ncol = n_contexts)
        rownames(ensemble_twas_wt) <- rownames(first_wt)
        colnames(ensemble_twas_wt) <- colnames(first_wt)

        for (i in which(matched)) {
          w_mat <- wt_list[[wt_keys[i]]]
          if (!is.matrix(w_mat)) w_mat <- matrix(w_mat, ncol = 1)
          if (!identical(dim(w_mat), dim(ensemble_twas_wt))) {
            warning("Weight matrix for '", wt_keys[i],
                    "' has inconsistent dimensions; skipping.")
            next
          }
          ensemble_twas_wt <- ensemble_twas_wt + zeta[i] * w_mat
        }

        # For univariate case, return as vector
        if (n_contexts == 1) {
          ensemble_twas_wt <- setNames(
            as.numeric(ensemble_twas_wt),
            rownames(ensemble_twas_wt)
          )
        }
      } else {
        warning("No matching weight keys found in twas_weight_list. ",
                "Expected keys like: ",
                paste(wt_keys[seq_len(min(3, K))], collapse = ", "))
      }
    }
  }

  list(
    method_coef = zeta,
    ensemble_twas_weights = ensemble_twas_wt,
    method_performance = method_rsq
  )
}

# =============================================================================
# Summary-statistics TWAS weight training pipeline
# =============================================================================

# Internal: RAISS-impute QTL z-scores for LD-panel variants missing from the
# QTL summary statistics. Used by twas_weights_sumstat_pipeline() when
# impute_missing = TRUE. Returns the (possibly widened) sumstats data frame
# with new rows for imputed variants. Imputed variants with R^2 below the
# threshold are dropped by RAISS's internal filter.
impute_missing_sumstats_for_LD <- function(sumstats, LD_mat, LD_data,
                                           impute_opts, verbose = 1) {
  ld_ids <- rownames(LD_mat)
  missing_ids <- setdiff(ld_ids, sumstats$variant_id)
  if (length(missing_ids) == 0) return(sumstats)

  # Build ref_panel covering all LD variants
  if (is(LD_data, "LDData")) {
    ld_ref_panel <- getRefPanel(LD_data)
  } else {
    ld_ref_panel <- parse_variant_id(ld_ids)
    ld_ref_panel$variant_id <- ld_ids
  }
  ref_cols <- c("chrom", "pos", "variant_id", "A1", "A2")
  if (!all(ref_cols %in% colnames(ld_ref_panel))) {
    warning("impute_missing_sumstats_for_LD: LD ref_panel missing required columns; skipping imputation.")
    return(sumstats)
  }
  if (!all(ref_cols %in% colnames(sumstats)) || !"z" %in% colnames(sumstats)) {
    warning("impute_missing_sumstats_for_LD: sumstats missing required columns; skipping imputation.")
    return(sumstats)
  }

  # RAISS requires inputs sorted by position (within each chromosome)
  ref_sorted <- ld_ref_panel[order(ld_ref_panel$chrom, ld_ref_panel$pos), ref_cols, drop = FALSE]
  known_sorted <- sumstats[order(sumstats$chrom, sumstats$pos), c(ref_cols, "z"), drop = FALSE]
  raiss_args <- c(list(
    ref_panel = ref_sorted,
    known_zscores = known_sorted,
    LD_matrix = LD_mat,
    verbose = (verbose >= 2)
  ), impute_opts)
  raiss_out <- tryCatch(do.call(raiss, raiss_args),
                        error = function(e) {
                          warning(sprintf("RAISS missing-sumstat imputation failed: %s", e$message))
                          NULL
                        })
  if (is.null(raiss_out) || is.null(raiss_out$result_filter)) return(sumstats)

  new_rows <- raiss_out$result_filter[
    !raiss_out$result_filter$variant_id %in% sumstats$variant_id, , drop = FALSE
  ]
  if (nrow(new_rows) == 0) return(sumstats)

  added <- new_rows[, c("variant_id", "chrom", "pos", "A1", "A2", "z"), drop = FALSE]
  if ("beta" %in% colnames(sumstats)) added$beta <- new_rows$z
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
#' fit. Returns a \code{TWASWeights} S4 object with \code{standardized = TRUE}
#' that feeds directly into \code{harmonize_twas} and \code{twas_analysis}.
#'
#' @param sumstats Data.frame with columns: \code{variant_id}, \code{A1},
#'   \code{A2}, \code{chrom}, \code{pos}, and either \code{z} or both
#'   \code{beta} and \code{se}.
#' @param LD_data LDData S4 object, or a legacy list with \code{LD_matrix},
#'   \code{LD_variants}, \code{ref_panel}. Can also be a plain correlation
#'   matrix (variant IDs taken from row/colnames).
#' @param n eQTL study sample size (scalar).
#' @param methods Named list of RSS weight methods and their arguments.
#'   Method names correspond to functions named
#'   \code{<method>_weights(stat, LD, ...)}. Defaults include lassosum_rss,
#'   prs_cs, sdpr, susie_rss, and susie_inf_rss.
#' @param p_thresholds Numeric vector of p-value thresholds for P+T weights.
#'   Set to NULL to skip.
#' @param check_ld_method LD matrix repair method: \code{"eigenfix"} (default),
#'   \code{"shrink"}, or NULL to skip.
#' @param qc_method RSS QC method for eQTL data: \code{"slalom"},
#'   \code{"dentist"}, or NULL/\code{"none"} to skip.
#' @param keep_indel Whether to keep indels during QC. Default TRUE.
#' @param pip_cutoff_to_skip PIP threshold for early stopping. Default 0 (off).
#' @param impute Whether to run RAISS imputation of LD-inconsistent variants
#'   flagged by QC (the QC re-imputation path). Default FALSE.
#' @param impute_missing Logical. When \code{TRUE}, RAISS imputes QTL z-scores
#'   for variants present in the LD reference but absent from the QTL
#'   summary statistics, after QC and before LD/sumstats intersection. This
#'   widens the sumstats panel available to the weight-learning methods so a
#'   richer set of weights can later be applied to GWAS. Independent of
#'   \code{impute}; both can be enabled together. Default \code{FALSE}.
#' @param impute_opts RAISS imputation parameters; shared by the \code{impute}
#'   QC re-imputation and the \code{impute_missing} missing-variant path.
#'   Imputed variants with \code{R2 < R2_threshold} are dropped.
#' @param var_y Phenotype variance. Default 1.
#' @param verbose Verbosity level.
#'
#' @return A list with:
#' \describe{
#'   \item{twas_weights}{A \code{TWASWeights} S4 object with
#'     \code{standardized = TRUE}.}
#'   \item{finemapping_result}{A \code{FineMappingResult} S4 object from the
#'     SuSiE-RSS fit, or NULL if no SuSiE-RSS method was used.}
#'   \item{qc_summary}{List with outlier counts and QC metadata.}
#' }
#'
#' @export
twas_weights_sumstat_pipeline <- function(
    sumstats, LD_data, n,
    methods = list(
      lassosum_rss = list(),
      prs_cs = list(phi = 1e-4, n_iter = 1000, n_burnin = 500, thin = 5),
      sdpr = list(iter = 1000, burn = 200, thin = 1, verbose = FALSE),
      susie_rss = list(),
      susie_inf_rss = list()
    ),
    p_thresholds = c(0.001, 0.05),
    check_ld_method = "eigenfix",
    qc_method = NULL,
    keep_indel = TRUE,
    pip_cutoff_to_skip = 0,
    impute = TRUE,
    impute_missing = FALSE,
    impute_opts = list(rcond = 0.01, R2_threshold = 0.6,
                       minimum_ld = 5, lamb = 0.01),
    var_y = 1, verbose = 1) {

  # -----------------------------------------------------------------------
  # 1. RSS QC on eQTL summary statistics
  # -----------------------------------------------------------------------
  needs_qc <- !is.null(qc_method) && !identical(qc_method, "none")
  if (needs_qc || impute || pip_cutoff_to_skip != 0) {
    qc_result <- summary_stats_qc(
      rss_input = list(sumstats = sumstats, n = n, var_y = var_y),
      LD_data = LD_data,
      keep_indel = keep_indel,
      pip_cutoff_to_skip = pip_cutoff_to_skip,
      qc_method = qc_method,
      impute = impute,
      impute_opts = impute_opts,
      return_on_skip = "null"
    )
    if (is.null(qc_result) || isSkipped(qc_result)) {
      return(list(twas_weights = NULL, finemapping_result = NULL,
                  qc_summary = list(skipped = TRUE)))
    }
    sumstats <- getRSSInput(qc_result)$sumstats
    qc_ld <- getLDData(qc_result)
    LD_mat <- if (is.null(qc_ld)) NULL else if (hasGenotypes(qc_ld)) getGenotypes(qc_ld) else getCorrelation(qc_ld)
    outlier_number <- getOutlierNumber(qc_result)
  } else {
    # No QC requested: extract LD matrix directly
    if (is.matrix(LD_data)) {
      LD_mat <- LD_data
    } else if (is(LD_data, "LDData")) {
      LD_mat <- getCorrelation(LD_data)
    } else {
      stop("LD_data must be a matrix or LDData object.")
    }
    outlier_number <- 0L
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
  variant_ids <- sumstats$variant_id
  b <- z / sqrt(n)
  stat <- list(b = b, cor = b, z = z, n = rep(n, p))

  # Optional RAISS imputation: fill QTL z-scores for LD-panel variants absent
  # from the QTL summary statistics. Widens the sumstats panel so weight
  # learners have access to a richer variant set; downstream intersection
  # with LD becomes a near-identity after this step.
  if (isTRUE(impute_missing) && !is.null(LD_mat) && !is.null(rownames(LD_mat))) {
    sumstats <- impute_missing_sumstats_for_LD(
      sumstats = sumstats,
      LD_mat = LD_mat,
      LD_data = LD_data,
      impute_opts = impute_opts,
      verbose = verbose
    )
    p <- nrow(sumstats)
    z <- sumstats$z
    variant_ids <- sumstats$variant_id
    b <- z / sqrt(n)
    stat <- list(b = b, cor = b, z = z, n = rep(n, p))
  }

  # Align LD matrix to sumstats variant order
  if (!is.null(rownames(LD_mat)) && !is.null(variant_ids)) {
    common <- intersect(variant_ids, rownames(LD_mat))
    if (length(common) < p) {
      idx <- match(common, variant_ids)
      sumstats <- sumstats[idx, , drop = FALSE]
      z <- sumstats$z
      variant_ids <- sumstats$variant_id
      b <- z / sqrt(n)
      stat <- list(b = b, cor = b, z = z, n = rep(n, length(z)))
      p <- length(z)
    }
    LD_mat <- LD_mat[variant_ids, variant_ids, drop = FALSE]
  }

  # -----------------------------------------------------------------------
  # 3. LD eigenfix (optional)
  # -----------------------------------------------------------------------
  if (!is.null(check_ld_method)) {
    ld_check <- check_ld(LD_mat, method = check_ld_method)
    if (ld_check$method_applied != "none") {
      if (verbose >= 1) {
        message(sprintf("check_ld: repaired LD via '%s' (min eigenvalue was %.2e, %d negative).",
                        ld_check$method_applied, ld_check$min_eigenvalue, ld_check$n_negative))
      }
    }
    LD_mat <- ld_check$R
  }

  # -----------------------------------------------------------------------
  # 4. Two-stage SuSiE-RSS (shared fit for susie_rss + susie_inf_rss)
  # -----------------------------------------------------------------------
  has_susie_rss <- "susie_rss" %in% names(methods)
  has_susie_inf_rss <- "susie_inf_rss" %in% names(methods)
  susie_fits <- NULL

  if (has_susie_rss && has_susie_inf_rss) {
    susie_args <- methods[["susie_rss"]]
    susie_inf_args <- methods[["susie_inf_rss"]]
    susie_fits <- fit_susie_inf_then_susie_rss(
      z = z, R = LD_mat, n = n,
      susie_inf_args = susie_inf_args,
      susie_args = susie_args
    )
  }

  # -----------------------------------------------------------------------
  # 5. P+T weights
  # -----------------------------------------------------------------------
  results <- list()
  if (!is.null(p_thresholds)) {
    pvals <- pchisq(z^2, df = 1, lower.tail = FALSE)
    for (thr in p_thresholds) {
      selected <- pvals < thr
      w <- ifelse(selected, stat$b, 0)
      results[[paste0("PT_", thr)]] <- w
    }
  }

  # -----------------------------------------------------------------------
  # 6. RSS method dispatch
  # -----------------------------------------------------------------------
  susie_rss_fit_for_fm <- NULL

  for (method_name in names(methods)) {
    fn_name <- paste0(method_name, "_weights")
    if (!exists(fn_name, mode = "function")) {
      warning(sprintf("Method '%s' not found (looking for function '%s'). Skipping.",
                      method_name, fn_name))
      next
    }

    method_args <- methods[[method_name]]

    # Build call arguments: separate pre-fitted objects from method_args
    call_args <- list(stat = stat, LD = LD_mat)
    if (method_name == "susie_rss" && !is.null(susie_fits)) {
      call_args[["susie_rss_fit"]] <- susie_fits$susie
    } else if (method_name == "susie_inf_rss" && !is.null(susie_fits)) {
      call_args[["susie_inf_rss_fit"]] <- susie_fits$susie_inf
    }

    # SuSiE-RSS methods use method_args; others spread args directly
    is_susie_rss_method <- method_name %in% c("susie_rss", "susie_inf_rss", "susie_ash_rss")
    if (is_susie_rss_method) {
      call_args[["method_args"]] <- method_args
    } else {
      call_args <- c(call_args, method_args)
    }

    tryCatch({
      w <- do.call(fn_name, call_args)
      # Capture retained fit for fine-mapping post-processing
      if (method_name == "susie_rss" && !is.null(attr(w, "fit"))) {
        susie_rss_fit_for_fm <- attr(w, "fit")
      } else if (method_name == "susie_inf_rss" && is.null(susie_rss_fit_for_fm) && !is.null(attr(w, "fit"))) {
        susie_rss_fit_for_fm <- attr(w, "fit")
      }
      results[[method_name]] <- as.numeric(w)
    }, error = function(e) {
      warning(sprintf("Method '%s' failed: %s", method_name, e$message))
      results[[method_name]] <<- rep(0, p)
    })
  }

  if (length(results) == 0) {
    return(list(twas_weights = NULL, finemapping_result = NULL,
                qc_summary = list(skipped = TRUE, reason = "all methods failed")))
  }

  # -----------------------------------------------------------------------
  # 7. Fine-mapping from SuSiE-RSS fit (reuses the same fit)
  # -----------------------------------------------------------------------
  finemapping_result <- NULL
  if (!is.null(susie_rss_fit_for_fm)) {
    fm_fits <- list(susie_rss = susie_rss_fit_for_fm)
    tryCatch({
      fm_output <- postprocess_finemapping_fits(
        fits = fm_fits,
        data_x = LD_mat,
        coverage = 0.95,
        signal_cutoff = 0.025,
        cs_input = "Xcorr"
      )
      if (!is.null(fm_output$finemapping_results$susie_rss$finemapping_result)) {
        finemapping_result <- fm_output$finemapping_results$susie_rss$finemapping_result
      }
    }, error = function(e) {
      warning(sprintf("Fine-mapping post-processing failed: %s", e$message))
    })
  }

  # -----------------------------------------------------------------------
  # 8. Package into TWASWeights S4
  # -----------------------------------------------------------------------
  weights_list <- lapply(results, function(w) {
    matrix(w, ncol = 1, dimnames = list(variant_ids, NULL))
  })

  twas_wt <- TWASWeights(
    weights = weights_list,
    variant_ids = variant_ids,
    standardized = TRUE,
    cv_performance = NULL
  )

  list(
    twas_weights = twas_wt,
    finemapping_result = finemapping_result,
    qc_summary = list(
      skipped = FALSE,
      n_variants_input = p,
      n_variants_after_qc = nrow(sumstats),
      outlier_number = outlier_number,
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
#' \code{\link{twas_multivariate_weights_pipeline}}. Bundles per-context RSS
#' QC, cross-context variant alignment, optional RAISS missing-variant
#' imputation, data-driven prior construction (reusing
#' \code{\link{build_mrmash_prior_matrices}} and the same FLASH / diagonal
#' covariance helpers as the individual-level pipeline), and multi-context
#' weight training via \code{\link{mrmash_rss_weights}} and/or
#' \code{\link{mvsusie_rss_weights}}.
#'
#' @param sumstats_list Named list of per-context sumstats data.frames. Each
#'   data.frame must contain \code{variant_id}, \code{chrom}, \code{pos},
#'   \code{A1}, \code{A2}, and either \code{z} or \code{beta}/\code{se}.
#'   List names become condition labels.
#' @param LD_data Shared \code{LDData} S4 object covering the union of
#'   variants. Each per-context sumstats is QC'd against this same LD panel.
#' @param n Per-context sample sizes; either a named numeric vector matching
#'   \code{names(sumstats_list)} or a single scalar to broadcast.
#' @param methods Named list of multivariate RSS weight methods to fit.
#'   Function names must match \code{<name>_weights}. Defaults to mr.mash-RSS
#'   and mvSuSiE-RSS with default arguments.
#' @param qc_method Per-context QC method passed to
#'   \code{\link{summary_stats_qc}}; one of \code{"slalom"}, \code{"dentist"},
#'   or NULL/\code{"none"}. Default \code{NULL} (basic harmonization only).
#' @param keep_indel Passed through to QC. Default TRUE.
#' @param impute Logical. If TRUE, run per-context RAISS re-imputation of
#'   LD-mismatch outliers (QC re-imputation). Default FALSE.
#' @param impute_missing Logical. If TRUE, after per-context QC and before
#'   cross-context alignment, RAISS imputes per-context z-scores for
#'   LD-reference variants absent from each context's sumstats. Widens the
#'   intersection used downstream. Default FALSE.
#' @param impute_opts Named list of RAISS parameters shared by both
#'   \code{impute} and \code{impute_missing}.
#' @param data_driven_prior_matrices Optional list of pre-computed prior
#'   matrices (with element \code{U}) passed to
#'   \code{\link{build_mrmash_prior_matrices}}. When NULL and
#'   \code{estimate_priors_from_sumstats = TRUE}, the pipeline estimates a
#'   data-driven covariance from the cross-context \code{Bhat} matrix via
#'   \code{\link{compute_cov_flash}}. If FLASH fails the error is allowed
#'   to propagate; supply \code{data_driven_prior_matrices} explicitly or
#'   set \code{estimate_priors_from_sumstats = FALSE} to bypass.
#' @param canonical_prior_matrices Passed to
#'   \code{\link{build_mrmash_prior_matrices}}. Default TRUE.
#' @param estimate_priors_from_sumstats Logical. When TRUE (default) and
#'   \code{data_driven_prior_matrices} is NULL, estimate data-driven priors
#'   from the cross-context Bhat matrix using
#'   \code{\link{compute_cov_flash}}. FLASH errors are not swallowed; see
#'   \code{data_driven_prior_matrices} for the explicit-prior path.
#' @param verbose Integer verbosity level.
#'
#' @return A list with
#' \describe{
#'   \item{twas_weights}{A \code{TWASWeights} S4 object with per-context
#'     weight matrices (variants x conditions).}
#'   \item{qc_summary}{Per-context QC and alignment counts.}
#'   \item{Z}{The aligned z-score matrix (variants x conditions) fed to the
#'     weight learners.}
#' }
#' @export
twas_multivariate_weights_sumstat_pipeline <- function(
    sumstats_list, LD_data, n,
    methods = list(mrmash_rss = list(), mvsusie_rss = list()),
    qc_method = NULL,
    keep_indel = TRUE,
    impute = FALSE,
    impute_missing = FALSE,
    impute_opts = list(rcond = 0.01, R2_threshold = 0.6,
                       minimum_ld = 5, lamb = 0.01),
    data_driven_prior_matrices = NULL,
    canonical_prior_matrices = TRUE,
    estimate_priors_from_sumstats = TRUE,
    verbose = 1) {

  # ----- 1. Validate inputs and normalise per-context n -----
  if (!is.list(sumstats_list) || length(sumstats_list) == 0L) {
    stop("sumstats_list must be a non-empty named list of sumstats data.frames.")
  }
  if (is.null(names(sumstats_list)) || any(names(sumstats_list) == "")) {
    stop("sumstats_list must be a named list; names become condition labels.")
  }
  conditions <- names(sumstats_list)
  K <- length(conditions)
  if (length(n) == 1L) n <- setNames(rep(as.numeric(n), K), conditions)
  if (is.null(names(n))) names(n) <- conditions
  missing_n <- setdiff(conditions, names(n))
  if (length(missing_n) > 0) {
    stop("n vector missing entries for conditions: ", paste(missing_n, collapse = ", "))
  }

  # ----- 2. Per-context QC + optional missing-variant imputation -----
  per_context_qc <- list()
  for (cond in conditions) {
    ss_c <- sumstats_list[[cond]]
    if (!is.null(ss_c$z) || (!is.null(ss_c$beta) && !is.null(ss_c$se))) {
      if (is.null(ss_c$z)) ss_c$z <- ss_c$beta / ss_c$se
    } else {
      stop(sprintf("Context %s: sumstats must contain z or (beta, se).", cond))
    }
    qc_record <- summary_stats_qc(
      rss_input = list(sumstats = ss_c, n = n[[cond]], var_y = 1),
      LD_data = LD_data,
      keep_indel = keep_indel,
      qc_method = if (is.null(qc_method)) "none" else qc_method,
      impute = impute, impute_opts = impute_opts,
      return_on_skip = "preprocess",
      study = cond
    )
    ss_qced <- getRSSInput(qc_record)$sumstats
    if (isTRUE(impute_missing) && nrow(ss_qced) > 0) {
      qc_ld <- getLDData(qc_record)
      LD_for_impute <- if (is.null(qc_ld)) {
        if (hasGenotypes(LD_data)) getGenotypes(LD_data) else getCorrelation(LD_data)
      } else if (hasGenotypes(qc_ld)) getGenotypes(qc_ld) else getCorrelation(qc_ld)
      if (!is.null(LD_for_impute) && !is.null(rownames(LD_for_impute))) {
        ss_qced <- impute_missing_sumstats_for_LD(
          sumstats = ss_qced, LD_mat = LD_for_impute,
          LD_data = LD_data, impute_opts = impute_opts, verbose = verbose
        )
      }
    }
    per_context_qc[[cond]] <- ss_qced
  }

  # ----- 3. Cross-context alignment (intersection on variant_id) -----
  variant_sets <- lapply(per_context_qc, function(df) df$variant_id)
  common_variants <- Reduce(intersect, variant_sets)
  if (length(common_variants) < 2) {
    return(list(twas_weights = NULL,
                qc_summary = list(skipped = TRUE,
                                  reason = "fewer than 2 shared variants across contexts"),
                Z = NULL))
  }

  # Use LD reference order for the common set
  ld_ids <- if (is(LD_data, "LDData")) getVariantIds(LD_data) else rownames(getCorrelation(LD_data))
  common_variants <- intersect(ld_ids, common_variants)

  # ----- 4. Build Z and Bhat/Shat matrices -----
  Z <- matrix(NA_real_, nrow = length(common_variants), ncol = K,
              dimnames = list(common_variants, conditions))
  Bhat <- Z; Shat <- Z
  n_vec <- as.numeric(n[conditions])
  for (k in seq_len(K)) {
    df <- per_context_qc[[conditions[k]]]
    df <- df[df$variant_id %in% common_variants, , drop = FALSE]
    idx <- match(common_variants, df$variant_id)
    Z[, k] <- df$z[idx]
    Bhat[, k] <- Z[, k] / sqrt(n_vec[k])
    Shat[, k] <- 1 / sqrt(n_vec[k])
  }

  # ----- 5. LD subset to common variants -----
  LD_full <- if (is(LD_data, "LDData")) getCorrelation(LD_data) else LD_data
  if (is.null(LD_full)) {
    if (is(LD_data, "LDData") && hasGenotypes(LD_data)) {
      X_ref <- getGenotypes(LD_data)
      LD_full <- compute_LD(X_ref[, common_variants, drop = FALSE], method = "sample")
    } else {
      stop("LD_data must provide either a correlation matrix or a genotype handle.")
    }
  }
  LD_mat <- LD_full[common_variants, common_variants, drop = FALSE]

  # ----- 6. Data-driven prior matrices from the aligned Bhat -----
  prior_input <- data_driven_prior_matrices
  if (is.null(prior_input) && isTRUE(estimate_priors_from_sumstats) && K >= 2) {
    if (verbose >= 1) message("Estimating data-driven prior matrices from Bhat ...")
    # Let FLASH errors propagate; callers can supply data_driven_prior_matrices
    # explicitly or set estimate_priors_from_sumstats = FALSE to bypass.
    prior_input <- list(U = list(flash = compute_cov_flash(Bhat)))
  }

  # ----- 7. Build stat object and dispatch weight methods -----
  stat <- list(z = Z, Bhat = Bhat, Shat = Shat, n = n_vec)
  results <- list()
  for (method_name in names(methods)) {
    fn_name <- paste0(method_name, "_weights")
    if (!exists(fn_name, mode = "function")) {
      warning(sprintf("Method '%s' not found (looking for function '%s'). Skipping.",
                      method_name, fn_name))
      next
    }
    method_args <- methods[[method_name]] %||% list()
    if (method_name == "mrmash_rss") {
      method_args$data_driven_prior_matrices <- method_args$data_driven_prior_matrices %||% prior_input
      method_args$canonical_prior_matrices  <- method_args$canonical_prior_matrices  %||% canonical_prior_matrices
    }
    fn <- get(fn_name, mode = "function")
    if (verbose >= 1) message("Fitting ", method_name, " ...")
    w <- tryCatch(
      do.call(fn, c(list(stat = stat, LD = LD_mat), method_args)),
      error = function(e) {
        warning(sprintf("Method '%s' failed: %s", method_name, e$message))
        NULL
      }
    )
    if (!is.null(w)) {
      if (!is.matrix(w)) w <- matrix(w, nrow = length(common_variants), ncol = K,
                                     dimnames = list(common_variants, conditions))
      results[[paste0(method_name, "_weights")]] <- w
    }
  }

  # ----- 8. Package into TWASWeights S4 -----
  if (length(results) == 0) {
    return(list(twas_weights = NULL,
                qc_summary = list(skipped = TRUE,
                                  reason = "all methods failed"),
                Z = Z))
  }
  twas_wt <- TWASWeights(
    weights = results,
    variant_ids = common_variants,
    standardized = TRUE,
    cv_performance = NULL
  )
  list(
    twas_weights = twas_wt,
    qc_summary = list(
      skipped = FALSE,
      n_per_context = vapply(per_context_qc, nrow, integer(1)),
      n_common = length(common_variants),
      conditions = conditions,
      methods_succeeded = names(results)
    ),
    Z = Z
  )
}
