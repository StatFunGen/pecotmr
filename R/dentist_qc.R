#' Detect Outliers Using Dentist Algorithm
#'
#' DENTIST (Detecting Errors iN analyses of summary staTISTics) is a quality control
#' tool for GWAS summary data. It uses linkage disequilibrium (LD) information from a reference
#' panel to identify and correct problematic variants by comparing observed GWAS statistics to
#' predicted values. It can detect errors in genotyping/imputation, allelic errors, and
#' heterogeneity between GWAS and LD reference samples.
#'
#' @param sum_stat A data frame containing summary statistics, including 'pos' or 'position' and 'z' or 'zscore' columns.
#' @param LD_mat A matrix containing LD (linkage disequilibrium) information.
#' @param nSample The number of samples.
#' @param window_size The size of the window for dividing the genomic region. Default is 2000000.
#' @param pValueThreshold The p-value threshold for significance. Default is 5e-8.
#' @param propSVD The proportion of singular value decomposition (SVD) to use. Default is 0.4.
#' @param gcControl Logical indicating whether genomic control should be applied. Default is FALSE.
#' @param nIter The number of iterations for the Dentist algorithm. Default is 10.
#' @param gPvalueThreshold The genomic p-value threshold for significance. Default is 0.05.
#' @param duprThreshold The absolute correlation r value threshold to be considered duplicate. Default is 0.99.
#' @param ncpus The number of CPU cores to use for parallel processing. Default is 1.
#' @param seed The random seed for reproducibility. Default is 999.
#' @param correct_chen_et_al_bug Logical indicating whether to correct the Chen et al. bug. Default is TRUE.
#'
#' @return A data frame containing the imputed result and detected outliers.
#'
#' The returned data frame includes the following columns:
#'
#' \describe{
#'   \item{\code{original_z}}{The original z-score values from the input \code{sum_stat}.}
#'   \item{\code{imputed_z}}{The imputed z-score values computed by the Dentist algorithm.}
#'   \item{\code{rsq}}{The coefficient of determination (R-squared) between original and imputed z-scores.}
#'   \item{\code{iter_to_correct}}{The number of iterations required to correct the z-scores, if applicable.}
#'   \item{\code{index_within_window}}{The index of the observation within the window.}
#'   \item{\code{index_global}}{The global index of the observation.}
#'   \item{\code{outlier_stat}}{The computed statistical value based on the original and imputed z-scores and R-squared.}
#'   \item{\code{outlier}}{A logical indicator specifying whether the observation is identified as an outlier based on the statistical test.}
#' }
#'
#'
#' @examples
#' # Example usage of dentist
#' dentist(sum_stat, LD_mat, nSample)
#'
#' @details
#' # correct_chen_et_al_bug may affect the result in three parts compared with the original DENTIST method:
#' # 1. the way that window is divided
#'    when there is only one window and window size >= position range, it will run two windows
#'    where the first window is the whole range (e.g., 0-4122) and second is the subset of it, (e.g., 742-4122), which adds no information
#'    so if we set correct_chen_et_al_bug=TRUE and this is the scenario for window_size and positions, then it only runs on the whole window
#'    A more general description for this is:
#'        - if correct_chen_et_al_bug = FALSE, then there can never be a single window
#'        - if window_size is invalid (<=0 or >=range):
#'            - if correct_chen_et_al_bug==TRUE, run dentist single window, so just window (0-4122)
#'            - if correct_chen_et_al_bug==FALSE, run in the dentist original way, so window (0-4122) and (742-4122)
#'        - if window_size is valid:
#'            - if correct_chen_et_al_bug==TRUE, divide the window but don't add additional window if the first window covers all, so just (0-4122)
#'            - if correct_chen_et_al_bug==FALSE, divide the window in the dentist original way, so window (0-4122) and (742-4122)
#' 2. comparison between iteration index t and nIter (explained in the source code)
#' 3. !grouping_tmp (explained in the source code)
#'
#' @export
dentist <- function(sum_stat, LD_mat, nSample,
                    window_size = 2000000, pValueThreshold = 5.0369e-8, propSVD = 0.4, gcControl = FALSE,
                    nIter = 10, gPvalueThreshold = 0.05, duprThreshold = 0.99, ncpus = 1, seed = 999, correct_chen_et_al_bug = TRUE) {
  # detect for column names and order by pos
  if (!any(tolower(c("pos", "position")) %in% tolower(colnames(sum_stat))) ||
    !any(tolower(c("z", "zscore")) %in% tolower(colnames(sum_stat)))) {
    stop("Input sum_stat is missing either 'pos'/'position' or 'z'/'zscore' column.")
  }
  # rename to common column name
  if (!tolower("pos") %in% tolower(colnames(sum_stat))) {
    colnames(sum_stat)[which(tolower(colnames(sum_stat)) %in% tolower(c("position")))] <- "pos"
  }

  if (!tolower("z") %in% tolower(colnames(sum_stat))) {
    colnames(sum_stat)[which(tolower(colnames(sum_stat)) %in% tolower(c("zscore")))] <- "z"
  }

  sum_stat <- sum_stat %>% arrange(pos)
  if (window_size <= 0 | ((window_size >= max(sum_stat$pos) - min(sum_stat$pos) | is.na(window_size)) & (correct_chen_et_al_bug == TRUE))) {
    dentist_result <- dentist_single_window(
      sum_stat$z, LD_mat, nSample,
      pValueThreshold, propSVD, gcControl,
      nIter, gPvalueThreshold, duprThreshold,
      ncpus, seed, correct_chen_et_al_bug
    )
  } else {
    # divide windows
    window_divided_res <- divide_into_windows(sum_stat$pos, window_size = window_size, correct_chen_et_al_bug = TRUE)
    # compute dentist result for each window
    dentist_result_by_window <- list()
    for (k in 1:nrow(window_divided_res)) {
      zScore_k <- sum_stat$z[window_divided_res$windowStartIdx[k]:window_divided_res$windowEndIdx[k]]
      LD_mat_k <- LD_mat[
        window_divided_res$windowStartIdx[k]:window_divided_res$windowEndIdx[k],
        window_divided_res$windowStartIdx[k]:window_divided_res$windowEndIdx[k]
      ]
      dentist_result_by_window[[k]] <- dentist_single_window(
        zScore_k, LD_mat_k, nSample,
        pValueThreshold, propSVD, gcControl,
        nIter, gPvalueThreshold, duprThreshold,
        ncpus, seed, correct_chen_et_al_bug
      )
    }
    # merge single window result and generate a final dentist_result (similar to dentist_result above)
    dentist_result <- merge_windows(dentist_result_by_window, window_divided_res)
  }
  return(dentist_result)
}

#' Perform DENTIST on a single window
#'
#' This function performs imputation of summary statistics for a single genomic window
#' using the Dentist algorithm.
#'
#' @param zScore A numeric vector containing the z-score values for variants within the window.
#' @param LD_mat A square matrix containing linkage disequilibrium (LD) information for variants within the window.
#' @param nSample The total number of samples.
#' @param pValueThreshold The p-value threshold for significance. Default is 5e-8.
#' @param propSVD The proportion of singular value decomposition (SVD) to use. Default is 0.4.
#' @param gcControl Logical indicating whether genomic control should be applied. Default is FALSE.
#' @param nIter The number of iterations for the Dentist algorithm. Default is 10.
#' @param gPvalueThreshold The genomic p-value threshold for significance. Default is 0.05.
#' @param duprThreshold The absolute correlation r value threshold to be considered duplicate. Default is 0.99.
#' @param ncpus The number of CPU cores to use for parallel processing. Default is 1.
#' @param seed The random seed for reproducibility. Default is 999.
#' @param correct_chen_et_al_bug Logical indicating whether to correct the Chen et al. bug. Default is TRUE.
#'
#' @return data frame includes columns representing the imputed summary statistics and outlier detected.
#'
#' @examples
#' # Example usage of dentist_impute_single_window
#' library(MASS)
#' library(corpcor)
#' set.seed(999)
#' # Set the number of SNPs, sample size, and number of outliers
#' n_snps <- 1000
#' sample_size <- 10000
#' n_outliers <- 5
#'
#' # Generate a correlation matrix with more off-diagonal correlation
#' cor_matrix <- matrix(0, nrow = n_snps, ncol = n_snps)
#' for (i in 1:(n_snps - 1)) {
#'   for (j in (i + 1):n_snps) {
#'     cor_matrix[i, j] <- runif(1, 0.2, 0.8) # Generate random correlations between 0.2 and 0.8
#'     cor_matrix[j, i] <- cor_matrix[i, j]
#'   }
#' }
#' diag(cor_matrix) <- 1
#'
#' # Convert the correlation matrix to a positive definite matrix
#' ld_matrix <- cov2cor(make.positive.definite(cor_matrix))
#'
#' # Simulate Z-scores based on the LD matrix
#' z_scores <- mvrnorm(n = 1, mu = rep(0, n_snps), Sigma = ld_matrix)
#'
#' # Introduce outliers
#' outlier_indices <- sample(1:n_snps, n_outliers)
#' z_scores[outlier_indices] <- rnorm(n_outliers, mean = 0, sd = 5)
#' dentist_single_window(zScore, LD_mat, nSample)
#'
#' @seealso
#' \code{\link{dentist}} for detecting outliers using the Dentist algorithm.
#'
#' @references
#' https://github.com/Yves-CHEN/DENTIST
#' @export
dentist_single_window <- function(zScore, LD_mat, nSample,
                                  pValueThreshold = 5e-8, propSVD = 0.4, gcControl = FALSE,
                                  nIter = 10, gPvalueThreshold = 0.05, duprThreshold = 0.99,
                                  ncpus = 1, seed = 999, correct_chen_et_al_bug = TRUE) {
  calculate_stat <- function(impOp_zScores, impOp_imputed, impOp_rsq) {
    (impOp_zScores - impOp_imputed)^2 / (1 - impOp_rsq)
  }

  outlier_test <- function(stat, lambda, alpha = 5e-8) {
    minusLogPvalueChisq <- function(stat) {
      p <- pchisq(stat, df = 1, lower.tail = FALSE)
      return(-log10(p))
    }
    ifelse(minusLogPvalueChisq(stat / lambda) > -log10(alpha), TRUE, FALSE)
  }
  # Check that number of variants cannot be below 2000
  if (length(zScore) < 2000) {
    warning("The number of variants is below 2000. The algorithm may not work as expected, as suggested by the original DENTIST.")
  }
  # Check that LD_mat dimensions match the length of zScore
  if (!is.matrix(LD_mat) || nrow(LD_mat) != ncol(LD_mat) || nrow(LD_mat) != length(zScore)) {
    stop("LD_mat must be a square matrix with dimensions equal to the length of zScore.")
  }
  # Remove dups
  org_Zscore <- zScore
  dedup_res <- NULL
  if (duprThreshold < 1.0) {
    dedup_res <- find_duplicate_variants(zScore, LD_mat, duprThreshold)
    num_dup <- sum(dedup_res$dupBearer != -1)
    if (num_dup > 0) {
      message(paste(num_dup, "duplicated variants out of a total of", length(zScore), "were found at r threshold of", duprThreshold))
    }
    zScore <- dedup_res$filteredZ
    LD_mat <- dedup_res$filteredLD
  }

  # Define a custom condition to capture warnings
  warning_handler <- function(w) {
    # Check if the warning message matches the specified pattern
    if (grepl("Adjusted rsq_eigen value exceeding 1", w$message)) {
      # Convert the warning to an error
      stop(w$message)
    }
    # Otherwise, invoke the default warning handler
    invokeRestart("muffleWarning")
  }
  res <- tryCatch(
    {
      dentist_iterative_impute(
        LD_mat, nSample, zScore,
        pValueThreshold, propSVD, gcControl, nIter,
        gPvalueThreshold, ncpus, seed, correct_chen_et_al_bug
      )
    },
    warning = warning_handler
  )
  
  res <- as.data.frame(res)
  # Recover dups
  if (duprThreshold < 1.0) {
    res <- add_dups_back_dentist(org_Zscore, res, dedup_res)
  }
  # detect outlier
  lambda_original <- 1
  res %>%
    mutate(
      outlier_stat = z_diff^2,
      outlier = outlier_test(outlier_stat, lambda_original)
    ) %>%
    select(-z_diff) %>%
    filter(!(imputed_z == 0 & rsq == 0))
}

#' Add duplicates back to DENTIST output
#'
#' This function takes the output from the DENTIST algorithm and adds back the duplicated variants
#' based on the output from the `find_duplicate_variants` function.
#' @param zScore The original zScore
#' @param dentist_output A data frame containing the output from the DENTIST algorithm.
#' @param find_dup_output A list containing the output from the `find_duplicate_variants` function.
#'
#' @return A data frame with duplicated variants added back and an additional column indicating duplicates.
#'
#' @noRd
add_dups_back_dentist <- function(zScore, dentist_output, find_dup_output) {
  # Extract relevant columns from the DENTIST output
  original_z <- dentist_output$original_z
  imputed_z <- dentist_output$imputed_z
  iter_to_correct <- dentist_output$iter_to_correct
  rsq <- dentist_output$rsq
  z_diff <- dentist_output$z_diff

  # Extract output from find_duplicate_variants
  dupBearer <- find_dup_output$dupBearer
  sign <- find_dup_output$sign

  # Get the number of rows in dupBearer
  nrows_dup <- length(dupBearer)

  if (nrow(dentist_output) != sum(dupBearer == -1)) {
    stop("The number of rows in the input data does not match the occurrences of -1 in dupBearer.")
  }

  if (length(zScore) != nrows_dup) {
    stop("Input zScore and find_dup_output have inconsistent dimension")
  }

  # Initialize assignIdx vector
  count <- 1
  assignIdx <- rep(0, nrows_dup)

  for (i in seq_along(dupBearer)) {
    if (dupBearer[i] == -1) {
      assignIdx[i] <- count
      count <- count + 1
    } else {
      assignIdx[i] <- dupBearer[i]
    }
  }

  # Create a new data frame to store the updated values
  updated_data <- data.frame(
    original_z = numeric(nrows_dup),
    imputed_z = numeric(nrows_dup),
    iter_to_correct = numeric(nrows_dup),
    rsq = numeric(nrows_dup),
    z_diff = numeric(nrows_dup),
    is_duplicate = logical(nrows_dup)
  )

  for (i in seq_len(nrows_dup)) {
    updated_data$original_z[i] <- zScore[i]
    updated_data$iter_to_correct[i] <- iter_to_correct[assignIdx[i]]
    updated_data$rsq[i] <- rsq[assignIdx[i]]
    # This may be wrong in sign but it does not really matter because we will take sqrt of it anyways
    updated_data$z_diff[i] <- z_diff[assignIdx[i]]
    if (dupBearer[i] == -1) {
      updated_data$imputed_z[i] <- imputed_z[assignIdx[i]]
      updated_data$is_duplicate[i] <- FALSE
    } else {
      updated_data$imputed_z[i] <- imputed_z[assignIdx[i]] * sign[i]
      updated_data$is_duplicate[i] <- TRUE
    }
  }

  return(updated_data)
}

#' Divide Genomic Region into Windows
#'
#' This function divides a genomic region into windows based on the specified window size and other parameters.
#'
#' @param pos A numeric vector containing the positions of variants.
#' @param window_size The size of the window for dividing the genomic region.
#' @param correct_chen_et_al_bug Logical indicating whether to correct the Chen et al. bug.
#'
#' @return A data frame containing information about the divided windows, including start and end indices for both windows and fillers.
#'
#' @details
#' - If \code{correct_chen_et_al_bug = FALSE}, then there can never be a single window.
#' - If the \code{window_size} is invalid (<=0 or >=range):
#'   - If \code{correct_chen_et_al_bug==TRUE}, the function runs with a single window, covering the entire genomic region.
#'   - If \code{correct_chen_et_al_bug==FALSE}, the function runs in the original Dentist way, dividing the genomic region into multiple windows.
#' - If the \code{window_size} is valid:
#'   - If \code{correct_chen_et_al_bug==TRUE}, the function divides the window but doesn't add additional window if the first window covers all.
#'   - If \code{correct_chen_et_al_bug==FALSE}, the function divides the window in the original Dentist way.
#'
#' @seealso
#' \code{\link{dentist_single_window}} for detecting outlier for summary statistics for a single window using the Dentist algorithm.
#'
#' @noRd
divide_into_windows <- function(pos, window_size, correct_chen_et_al_bug) {
  windowStartIdx <- c()
  windowEndIdx <- c()
  fillStartIdx <- c()
  fillEndIdx <- c()
  input_pos_start <- min(pos)
  input_pos_end <- max(pos)
  pos_range <- input_pos_end - input_pos_start
  if (window_size <= 0 | (window_size >= pos_range & correct_chen_et_al_bug == TRUE)) {
    windowStartIdx <- 1
    windowEndIdx <- length(pos)
    fillStartIdx <- 1
    fillEndIdx <- length(pos)
  } else {
    for (i in 1:(2 * ceiling(pos_range / window_size))) {
      start_pointer <- ((i - 1) * window_size) * 0.5 + 1
      end_pointer <- start_pointer + window_size
      start_idx <- which.min(abs(pos - start_pointer - input_pos_start))
      windowStartIdx <- c(windowStartIdx, start_idx)
      end_idx <- which.min(abs(pos - end_pointer - input_pos_start))
      windowEndIdx <- c(windowEndIdx, end_idx)
    }
    # avoid the situation where the last window is too small
    cutoff <- which.min(abs(pos - input_pos_end + window_size))
    closest_idx <- which.min(abs(windowStartIdx - cutoff))
    windowStartIdx <- windowStartIdx[windowStartIdx <= cutoff | seq_along(windowStartIdx) == closest_idx]
    windowEndIdx <- windowEndIdx[1:length(windowStartIdx)]

    if (length(windowStartIdx) == 1 & (windowEndIdx[1] < 4122 | correct_chen_et_al_bug == FALSE)) {
      # to avoid only 1 window
      windowEndIdx <- c(windowEndIdx, length(pos))
      start_pointer <- (pos[windowEndIdx] - pos[windowStartIdx]) * 0.25 + pos[windowStartIdx] # position
      start_idx <- which.min(abs(pos - start_pointer))
      windowStartIdx <- c(windowStartIdx, start_idx)
    }
    # decide the fillStartIdx and fillEndIdx
    if (length(windowStartIdx) == 1) { # for single window (correct_chen_et_al_bug must be TRUE)
      fillStartIdx <- 1
      fillEndIdx <- length(pos)
    } else { # for multiple windows
      for (i in 1:length(windowStartIdx)) {
        if (i == 1 & i != length(windowStartIdx)) {
          # for the first window and not the last window, fill Start from 1 and fill End is 0.75 quantile
          fill_start_idx <- 1
          fillStartIdx <- c(fillStartIdx, fill_start_idx)
          fill_end_pointer <- 0.75 * window_size + input_pos_start
          fill_end_idx <- which.max(pos[pos <= fill_end_pointer])
          fillEndIdx <- c(fillEndIdx, fill_end_idx)
        } else if (i != length(windowStartIdx)) {
          # for the windows in the middle, the start is the next variant of the previous fillendIdx, and the end is the 0.75 quantile
          fillStartIdx <- c(fillStartIdx, fillEndIdx[i - 1] + 1)
          fill_end_pointer <- pos[windowStartIdx[i]] + 0.75 * window_size
          fill_end_idx <- which.max(pos[pos <= fill_end_pointer])
          fillEndIdx <- c(fillEndIdx, fill_end_idx)
        }
      }
      # for the last window, the start is the next variant of the previous fillEndIdx, and the end is the last variant
      fillStartIdx <- c(fillStartIdx, fillEndIdx[length(fillEndIdx)] + 1)
      fillEndIdx <- c(fillEndIdx, length(pos))
    }
  }
  # combine window information into a data frame
  window_index <- 1:length(windowStartIdx)
  window_divided_res <- data.frame(
    windowIdx = window_index,
    windowStartIdx = windowStartIdx,
    windowEndIdx = windowEndIdx,
    fillStartIdx = fillStartIdx,
    fillEndIdx = fillEndIdx
  )
  # check if the divided windows are valid
  invalid_ends <- window_divided_res$fillStartIdx[1] != 1 | window_divided_res$fillEndIdx[nrow(window_divided_res)] != length(pos)
  if (nrow(window_divided_res) == 1 & invalid_ends) {
    stop("Invalid window divided!")
  } else if (nrow(window_divided_res) > 1) {
    invalid_middle <- 0
    for (i in 1:(nrow(window_divided_res) - 1)) {
      if (fillStartIdx[i + 1] - fillEndIdx[i] != 1) {
        invalid_middle <- invalid_middle + 1
      }
    }
    if (invalid_ends | invalid_middle != 0) {
      stop("Invalid window divided!")
    }
  }
  return(window_divided_res)
}

#' Merge dentist Results by Window
#'
#' This function merges DENTIST results by window into a single data frame.
#'
#' @param dentist_result_by_window A list containing imputed results for each window.
#' @param window_divided_res A data frame containing information about the divided windows.
#'
#' @return A data frame containing merged results.
#'
#' @details
#' The function checks if the number of imputed results matches the number of windows.
#' It then merges the results by window, adding an index within the window and a global index.
#' Finally, it extracts the results within the fillers and combines them into a single data frame.
#'
#' @noRd
merge_windows <- function(dentist_result_by_window, window_divided_res) {
  if (length(dentist_result_by_window) != nrow(window_divided_res)) {
    stop("Different number of windows and imputed results!")
  }
  merged_results <- c()
  for (k in 1:nrow(window_divided_res)) {
    imputed_k <- dentist_result_by_window[[k]]
    imputed_k$index_within_window <- seq(1:nrow(imputed_k))
    imputed_k <- imputed_k %>%
      mutate(index_global = index_within_window + window_divided_res$windowStartIdx[k] - 1)
    extracted_results <- imputed_k %>%
      filter(index_global >= window_divided_res$fillStartIdx[k] & index_global <= window_divided_res$fillEndIdx[k])
    merged_results <- rbind(merged_results, extracted_results)
  }
  return(merged_results)
}
