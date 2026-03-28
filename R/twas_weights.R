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
#' @importFrom future plan multicore availableCores
#' @importFrom furrr future_map furrr_options
#' @importFrom purrr map
#' @export
twas_weights_cv <- function(X, Y, fold = NULL, sample_partitions = NULL, weight_methods = NULL, max_num_variants = NULL, variants_to_keep = NULL, num_threads = 1, ...) {
  split_data <- function(X, Y, sample_partition, fold) {
    if (is.null(rownames(X))) {
      warning("Row names in X are missing. Using row indices.")
      rownames(X) <- 1:nrow(X)
    }
    if (is.null(rownames(Y))) {
      warning("Row names in Y are missing. Using row indices.")
      rownames(Y) <- 1:nrow(Y)
    }
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
    message(paste("Y converted to matrix of", nrow(Y), "rows and", ncol(Y), "columns."))
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
    weight_methods <- lapply(setNames(nm = weight_methods), function(x) list())
  }

  if (!exists(".Random.seed")) {
    message("! No seed has been set. Please set seed for reproducable result. ")
  }

  # Select variants if necessary
  if (!is.null(max_num_variants) && ncol(X) > max_num_variants) {
    if (!is.null(variants_to_keep) && length(variants_to_keep) > 0) {
      variants_to_keep <- intersect(variants_to_keep, colnames(X))
      remaining_columns <- setdiff(colnames(X), variants_to_keep)
      if (length(variants_to_keep) < max_num_variants) {
        additional_columns <- sample(remaining_columns, max_num_variants - length(variants_to_keep), replace = FALSE)
        selected_columns <- union(variants_to_keep, additional_columns)
        message(sprintf(
          "Including %d specified variants and randomly selecting %d additional variants, for a total of %d variants out of %d for cross-validation purpose.",
          length(variants_to_keep), length(additional_columns), length(selected_columns), ncol(X)
        ))
      } else {
        selected_columns <- sample(variants_to_keep, max_num_variants, replace = FALSE)
        message(paste("Randomly selecting", length(selected_columns), "out of", length(variants_to_keep), "input variants for cross validation purpose."))
      }
    } else {
      selected_columns <- sort(sample(ncol(X), max_num_variants, replace = FALSE))
      message(paste("Randomly selecting", length(selected_columns), "out of", ncol(X), "variants for cross validation purpose."))
    }
    X <- X[, selected_columns, drop = FALSE]
  }

  # Create or use provided folds
  if (!is.null(fold)) {
    if (!is.null(sample_partitions)) {
      if (fold != length(unique(sample_partitions$Fold))) {
        message(paste0(
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
    num_cores <- ifelse(num_threads == -1, availableCores(), num_threads)
    num_cores <- min(num_cores, availableCores())

    cv_args <- list(...)

    # Perform CV with parallel processing
    compute_method_predictions <- function(j) {
      dat_split <- split_data(X, Y, sample_partition = sample_partition, fold = j)
      X_train <- dat_split$Xtrain
      Y_train <- dat_split$Ytrain
      X_test <- dat_split$Xtest
      Y_test <- dat_split$Ytest

      # Remove columns with zero standard error
      valid_columns <- apply(X_train, 2, function(col) sd(col) != 0)
      X_train <- X_train[, valid_columns, drop = FALSE]
      X_train <- filter_X_with_Y(X_train, Y_train, missing_rate_thresh = 1, maf_thresh = NULL)
      valid_columns <- colnames(X_train)
      # X_test <- X_test[, valid_columns, drop=FALSE]

      setNames(lapply(names(weight_methods), function(method) {
        args <- weight_methods[[method]]

        if (method %in% multivariate_weight_methods) {
          # Apply multivariate method to entire Y for this fold
          if (!is.null(cv_args$data_driven_prior_matrices_cv)) {
            if (method == "mrmash_weights") {
              args$data_driven_prior_matrices <- cv_args$data_driven_prior_matrices_cv[[j]]
            }
            if (method == "mvsusie_weights") {
              args$prior_variance <- cv_args$reweighted_data_driven_prior_matrices_cv[[j]]
            }
          }
          weights_matrix <- do.call(method, c(list(X = X_train, Y = Y_train), args))
          rownames(weights_matrix) <- colnames(X_train)
          # Adjust the weights matrix to include zeros for invalid columns
          full_weights_matrix <- matrix(0, nrow = ncol(X), ncol = ncol(Y))
          rownames(full_weights_matrix) <- colnames(X)
          colnames(full_weights_matrix) <- colnames(Y)
          full_weights_matrix[valid_columns, ] <- weights_matrix[valid_columns, ]
          Y_pred <- X_test %*% full_weights_matrix
          rownames(Y_pred) <- rownames(X_test)
          return(Y_pred)
        } else {
          Y_pred <- sapply(1:ncol(Y_train), function(k) {
            weights <- do.call(method, c(list(X = X_train, y = Y_train[, k]), args))
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
      }), names(weight_methods))
    }

    if (num_cores >= 2) {
      plan(multicore, workers = num_cores)
      fold_results <- future_map(1:fold, compute_method_predictions, .options = furrr_options(seed = TRUE, globals = c("sample_partition", "weight_methods", "args", "cv_args")))
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
          message(paste0(
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
#' @return A list where each element is named after a method and contains the weight matrix produced by that method.
#'
#' @export
#' @importFrom future plan multicore availableCores
#' @importFrom furrr future_map furrr_options
#' @importFrom purrr map exec
#' @importFrom rlang !!!
twas_weights <- function(X, Y, weight_methods, num_threads = 1) {
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
    weight_methods <- lapply(setNames(nm = weight_methods), function(x) list())
  }

  # Determine number of cores to use
  num_cores <- ifelse(num_threads == -1, availableCores(), num_threads)
  num_cores <- min(num_cores, availableCores())

  compute_method_weights <- function(method_name, weight_methods) {
    # Hardcoded vector of multivariate methods
    multivariate_weight_methods <- c("mrmash_weights", "mvsusie_weights")
    args <- weight_methods[[method_name]]

    # Remove columns with zero standard error
    valid_columns <- apply(X, 2, function(col) sd(col) != 0)
    X_filtered <- as.matrix(X[, valid_columns, drop = FALSE])

    if (method_name %in% multivariate_weight_methods) {
      # Apply multivariate method
      weights_matrix <- do.call(method_name, c(list(X = X_filtered, Y = Y), args))
      if (nrow(weights_matrix) != length(valid_columns)) weights_matrix <- weights_matrix[names(valid_columns), , drop = FALSE]
    } else {
      # Apply univariate method to each column of Y
      # Initialize it with zeros to avoid NA
      weights_matrix <- matrix(0, nrow = ncol(X_filtered), ncol = ncol(Y))

      for (k in 1:ncol(Y)) {
        weights_vector <- do.call(method_name, c(list(X = X_filtered, y = Y[, k]), args))
        if (is.matrix(weights_vector)) weights_vector <- weights_vector[, k]
        weights_matrix[, k] <- weights_vector
      }
    }

    # Adjust the weights matrix to include zeros for invalid columns
    full_weights_matrix <- matrix(0, nrow = ncol(X), ncol = ncol(Y))
    rownames(full_weights_matrix) <- colnames(X)
    colnames(full_weights_matrix) <- colnames(Y)
    full_weights_matrix[valid_columns, ] <- weights_matrix

    return(full_weights_matrix)
  }

  if (num_cores >= 2) {
    # Set up parallel backend to use multiple cores
    plan(multicore, workers = num_cores)
    weights_list <- names(weight_methods) %>% future_map(compute_method_weights, weight_methods, .options = furrr_options(seed = TRUE))
  } else {
    weights_list <- names(weight_methods) %>% map(compute_method_weights, weight_methods)
  }
  names(weights_list) <- names(weight_methods)

  if (!is.null(colnames(X))) {
    weights_list <- lapply(weights_list, function(x) {
      rownames(x) <- colnames(X)
      return(x)
    })
  }
  return(weights_list)
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
  setNames(lapply(weights_list, function(w) X %*% w), gsub("_weights", "_predicted", names(weights_list)))
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
#' @param cv_folds The number of folds to use for cross-validation. Set to 0 to skip cross-validation. Defaults to 5.
#' @param weight_methods List of methods to use to compute weights for TWAS; along with their parameters.
#' @param max_cv_variants The maximum number of variants to be included in cross-validation. Defaults to -1 which means no limit.
#' @param cv_threads The number of threads to use for parallel computation in cross-validation. Defaults to 1.
#' @param cv_weight_methods List of methods to use for cross-validation. If NULL, uses the same methods as weight_methods.
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
                                  cv_folds = 5,
                                  sample_partition = NULL,
                                  weight_methods = list(
                                    enet_weights = list(),
                                    lasso_weights = list(),
                                    bayes_r_weights = list(),
                                    bayes_l_weights = list(),
                                    mrash_weights = list(init_prior_sd = TRUE, max.iter = 100),
                                    susie_weights = list(refine = FALSE, init_L = 5, max_L = 20)
                                  ),
                                  max_cv_variants = -1,
                                  cv_threads = 1,
                                  cv_weight_methods = NULL) {
  res <- list()
  st <- proc.time()
  message("Performing TWAS weights computation for univariate analysis methods ...")

  # TWAS weights and predictions
  if (!is.null(susie_fit) && !is.null(weight_methods$susie_weights)) {
    weight_methods$susie_weights <- list(susie_fit = susie_fit)
    res$susie_weights_intermediate <- susie_fit[c("mu", "lbf_variable", "X_column_scale_factors", "pip")]
    if (!is.null(susie_fit$sets$cs)) {
      res$susie_weights_intermediate$cs_variants <- setNames(lapply(susie_fit$sets$cs, function(L) colnames(X)[L]), names(susie_fit$sets$cs))
      res$susie_weights_intermediate$cs_purity <- susie_fit$sets$purity
    }
  }
  res$twas_weights <- twas_weights(X, y, weight_methods = weight_methods)
  res$twas_predictions <- twas_predict(X, res$twas_weights)

  if (cv_folds > 1) {
    # A few cutting corners to run CV faster at the disadvantage of SuSiE and mr.ash:
    # 1. reset SuSiE to not using refine or adaptive L but to use L from previous analysis
    # 2. at most 100 iterations for mr.ash allowed
    # 3. only use a subset of variants randomly selected to avoid bias
    if (!is.null(susie_fit) && !is.null(weight_methods$susie_weights)) {
      max_L <- length(susie_fit$V)
      weight_methods$susie_weights <- list(refine = FALSE, init_L = max_L, max_L = max_L)
    }
    if (is.null(cv_weight_methods)) {
      # Filter function to exclude methods with all zero weights and track removed methods
      is_all_zero <- sapply(res$twas_weights, function(weights) all(weights == 0, na.rm = TRUE))

      # Identify removed methods
      removed_methods <- names(weight_methods)[is_all_zero]

      # Issue a warning if any methods have been removed
      if (length(removed_methods) > 0) {
        warning(sprintf(
          "Methods %s are removed from CV because all their weights are zeros.",
          paste(removed_methods, collapse = ", ")
        ))
      }
      # Extract the filtered methods retaining their specific configurations
      cv_weight_methods <- names(weight_methods)[!is_all_zero]
    }

    variants_for_cv <- c()
    if (max_cv_variants <= 0) {
      max_cv_variants <- Inf
    }
    if (ncol(X) > max_cv_variants) {
      variants_for_cv <- sample(colnames(X), max_cv_variants, replace = FALSE)
    }

    message("Performing cross-validation to assess TWAS weights ...")
    res$twas_cv_result <- twas_weights_cv(
      X,
      y,
      fold = cv_folds,
      sample_partition = sample_partition,
      weight_methods = cv_weight_methods,
      max_num_variants = max_cv_variants,
      num_threads = cv_threads,
      variants_to_keep = if (length(variants_for_cv) > 0) variants_for_cv else NULL
    )
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
#' @param cv_folds The number of folds to use for cross-validation. Defaults to 5. Set to 0 to skip cross-validation.
#' @param sample_partition An optional vector specifying the partition of samples for cross-validation. If NULL, a random partition is generated.
#' @param data_driven_prior_matrices A list of data-driven covariance matrices for mr.mash weights. Defaults to NULL.
#' @param data_driven_prior_matrices_cv A list of data-driven covariance matrices for mr.mash weights in cross-validation. Defaults to NULL.
#' @param canonical_prior_matrices If TRUE, computes canonical covariance matrices for mr.mash. Defaults to FALSE.
#' @param mvsusie_max_iter The maximum number of iterations for mvSuSiE. Defaults to 200.
#' @param mrmash_max_iter The maximum number of iterations for mr.mash. Defaults to 5000.
#' @param max_cv_variants The maximum number of variants to be included in cross-validation. Defaults to -1 which means no limit.
#' @param cv_threads The number of threads to use for parallel computation in cross-validation. Defaults to 1.
#' @param verbose If TRUE, provides more detailed output during execution. Defaults to FALSE.
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
    cv_folds = 5,
    sample_partition = NULL,
    data_driven_prior_matrices = NULL,
    data_driven_prior_matrices_cv = NULL,
    canonical_prior_matrices = FALSE,
    mvsusie_max_iter = 200,
    mrmash_max_iter = 5000,
    max_cv_variants = -1,
    cv_threads = 1,
    verbose = FALSE) {
  copy_twas_results <- function(context_names, variant_names, twas_weight, twas_predictions) {
    res <- setNames(vector("list", length(context_names)), context_names)
    for (i in names(res)) {
      if (i %in% colnames(twas_weights_res[[1]])) {
        res[[i]]$twas_weights <- lapply(twas_weight, function(wgts) {
          wgts[, i]
        })
        res[[i]]$twas_predictions <- lapply(twas_predictions, function(pred) {
          pred[, i]
        })
        res[[i]]$variant_names <- variant_names
      }
    }
    return(res)
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
  message("Extracting TWAS weights for multivariate analysis methods ...")
  # get TWAS weights
  twas_weights_res <- twas_weights(X = X, Y = Y, weight_methods = weight_methods)
  # get TWAS predictions for possible next steps such as computing correlations between predicted expression values
  twas_predictions <- twas_predict(X, twas_weights_res)

  # copy TWAS results by condition
  res <- copy_twas_results(colnames(Y), mnm_fit$variant_names, twas_weights_res, twas_predictions)

  # Perform cross-validation if specified
  if (cv_folds > 1) {
    # max_L <- length(which(mnm_fit$mvsusie_fitted$V > 1E-9)) + 2
    # To be fair in comparion with other methods mvSuSiE should have the same input max_L as when the weights were originally computed
    max_L <- length(mnm_fit$mvsusie_fitted$V)
    weight_methods <- list(
      mrmash_weights = list(
        data_driven_prior_matrices = data_driven_prior_matrices,
        canonical_prior_matrices = canonical_prior_matrices,
        max_iter = mrmash_max_iter,
        verbose = verbose
      ),
      mvsusie_weights = list(
        prior_variance = mnm_fit$reweighted_data_driven_prior_matrices,
        residual_variance = mnm_fit$mrmash_fitted$V,
        L = max_L,
        max_iter = mvsusie_max_iter,
        verbosity = verbose
      )
    )

    # Filter function to exclude methods with all zero weights and track removed methods
    is_all_zero <- sapply(twas_weights_res, function(weights) all(weights == 0, na.rm = TRUE))

    # Identify removed methods
    removed_methods <- names(weight_methods)[is_all_zero]

    # Issue a warning if any methods have been removed
    if (length(removed_methods) > 0) {
      warning(sprintf(
        "Methods %s are removed from CV because all their weights are zeros.",
        paste(removed_methods, collapse = ", ")
      ))
    }

    # Extract the filtered methods retaining their specific configurations
    weight_methods <- weight_methods[!is_all_zero]

    variants_for_cv <- c()
    if (max_cv_variants <= 0) max_cv_variants <- Inf
    if (ncol(X) > max_cv_variants) {
      variants_for_cv <- sample(colnames(X), max_cv_variants, replace = FALSE)
    }
    message("Performing cross-validation to assess TWAS weights ...")
    twas_cv_result <- twas_weights_cv(
      X = X, Y = Y, fold = cv_folds,
      weight_methods = weight_methods,
      sample_partition = sample_partition,
      num_threads = cv_threads,
      max_num_variants = max_cv_variants,
      variants_to_keep = if (length(variants_for_cv) > 0) variants_for_cv else NULL,
      data_driven_prior_matrices_cv = data_driven_prior_matrices_cv,
      reweighted_data_driven_prior_matrices_cv = mnm_fit$reweighted_data_driven_prior_matrices_cv
    )
    res <- copy_twas_cv_results(res, twas_cv_result)
  }
  total_time_elapsed <- proc.time() - st
  for (i in seq_along(res)) {
    res[[i]]$total_time_elapsed <- total_time_elapsed
  }
  return(res)
}
