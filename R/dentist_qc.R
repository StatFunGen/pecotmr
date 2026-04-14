#' Resolve LD Input: Accept Either R (LD matrix) or X (Genotype Matrix)
#'
#' Internal helper that validates and resolves the LD input for QC functions.
#' Exactly one of \code{R} or \code{X} must be provided. When \code{X} is
#' provided, LD is computed via \code{compute_LD(X)} and \code{nSample}
#' defaults to \code{nrow(X)}.
#'
#' @param R Square LD correlation matrix, or NULL.
#' @param X Genotype matrix (samples x SNPs), or NULL.
#' @param nSample Sample size. Required when \code{R} is provided and
#'   \code{need_nSample} is TRUE; inferred from \code{X} when \code{X} is provided.
#' @param need_nSample Logical; if TRUE, \code{nSample} must be available
#'   (either provided or inferred from \code{X}).
#'
#' @return A list with components \code{R} (LD correlation matrix) and
#'   \code{nSample} (integer or NULL).
#'
#' @noRd
resolve_LD_input <- function(R = NULL, X = NULL, nSample = NULL, need_nSample = FALSE,
                             ld_method = "sample") {
  if (is.null(R) && is.null(X)) {
    stop("Either R (LD matrix) or X (genotype matrix) must be provided.")
  }
  if (!is.null(R) && !is.null(X)) {
    stop("Provide either R or X, not both.")
  }
  if (!is.null(X)) {
    if (!is.matrix(X)) X <- as.matrix(X)
    if (is.null(nSample)) nSample <- nrow(X)
    R <- compute_LD(X, method = ld_method)
  }
  if (need_nSample && is.null(nSample)) {
    stop("nSample is required when providing an LD matrix R.")
  }
  list(R = R, nSample = nSample)
}

#' Detect Outliers Using Dentist Algorithm
#'
#' DENTIST (Detecting Errors iN analyses of summary staTISTics) is a quality control
#' tool for GWAS summary data. It uses linkage disequilibrium (LD) information from a reference
#' panel to identify and correct problematic variants by comparing observed GWAS statistics to
#' predicted values. It can detect errors in genotyping/imputation, allelic errors, and
#' heterogeneity between GWAS and LD reference samples.
#'
#' @param sum_stat A data frame containing summary statistics, including 'pos' or 'position' and 'z' or 'zscore' columns.
#' @param R Square LD correlation matrix. Provide either \code{R} or \code{X}.
#' @param X Genotype matrix (samples x SNPs). If provided, LD is computed via
#'   \code{compute_LD(X)} and \code{nSample} defaults to \code{nrow(X)}.
#' @param nSample The number of samples in the LD reference panel (NOT the GWAS sample
#'   size). This controls the SVD truncation rank K = min(idx_size, nSample) * propSVD.
#'   Required when \code{R} is provided; inferred from \code{X} when \code{X} is provided.
#' @param window_size The size of the window for dividing the genomic region
#'   in distance mode (base pairs). Default is 2000000 (2 Mb). Only used when
#'   \code{window_mode = "distance"}.
#' @param window_mode Character string specifying the windowing strategy:
#'   \code{"distance"} (default) creates windows by physical distance using
#'   \code{\link{segment_by_dist}} (C++ \code{--wind-dist}), and
#'   \code{"count"} creates windows by variant count using
#'   \code{\link{segment_by_count}} (C++ \code{--wind}).
#' @param pValueThreshold The p-value threshold for significance. Default is 5e-8.
#' @param propSVD The proportion of singular value decomposition (SVD) to use. Default is 0.4.
#' @param gcControl Logical indicating whether genomic control should be applied. Default is FALSE.
#' @param nIter The number of iterations for the Dentist algorithm. Default is 10.
#' @param gPvalueThreshold The genomic p-value threshold for significance. Default is 0.05.
#' @param duprThreshold The absolute correlation r value threshold to be considered duplicate. Default is 0.99.
#' @param ncpus The number of CPU cores to use for parallel processing. Default is 1.
#' @param correct_chen_et_al_bug Logical indicating whether to correct the Chen et al. bug. Default is TRUE.
#' @param min_dim In distance mode: minimum number of SNPs per block (default 2000).
#'   In count mode: the number of variants per window (i.e., the window size).
#' @param ld_method Character string specifying the LD computation method when
#'   \code{X} is provided. Passed to \code{\link{compute_LD}}. One of
#'   \code{"sample"} (default), \code{"population"}, or \code{"gcta"}.
#'   Ignored when \code{R} is provided directly.
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
#' @examples
#' # Example usage of dentist
#' dentist(sum_stat, R = LD_mat, nSample = nSample)
#'
#' @details
#' Windowing supports two modes matching the original DENTIST C++ binary:
#' \itemize{
#'   \item \code{"distance"} (default): Uses the \code{segmentingByDist} algorithm
#'     (C++ \code{--wind-dist}), implemented in \code{\link{segment_by_dist}}.
#'     Windows span a fixed physical distance (\code{window_size} bp).
#'   \item \code{"count"}: Uses the \code{segmentedQCed} algorithm
#'     (C++ \code{--wind}), implemented in \code{\link{segment_by_count}}.
#'     Windows contain a fixed number of variants (\code{min_dim}).
#'     Useful when regions have sparse variants where distance-based windows
#'     would create windows with too few variants.
#' }
#' The \code{correct_chen_et_al_bug} parameter affects the iterative filtering
#' in two ways:
#' \enumerate{
#'   \item Comparison between iteration index \code{t} and \code{nIter} (explained in source code)
#'   \item The \code{!grouping_tmp} operator bug (explained in source code)
#' }
#'
#' @export
dentist <- function(sum_stat, R = NULL, X = NULL, nSample = NULL,
                    window_size = 2000000, window_mode = c("distance", "count"),
                    pValueThreshold = 5.0369e-8, propSVD = 0.4, gcControl = FALSE,
                    nIter = 10, gPvalueThreshold = 0.05, duprThreshold = 0.99, ncpus = 1,
                    correct_chen_et_al_bug = TRUE, min_dim = 2000,
                    ld_method = "sample") {
  # Resolve LD matrix and sample size from R or X
  resolved <- resolve_LD_input(R = R, X = X, nSample = nSample, need_nSample = TRUE,
                               ld_method = ld_method)
  LD_mat <- resolved$R
  nSample <- resolved$nSample

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

  window_mode <- match.arg(window_mode)

  # If the data has fewer SNPs than min_dim, run as a single window directly.
  n_snps <- nrow(sum_stat)
  if (n_snps < min_dim) {
    dentist_result <- dentist_single_window(
      sum_stat$z, R = LD_mat, nSample = nSample,
      pValueThreshold = pValueThreshold, propSVD = propSVD, gcControl = gcControl,
      nIter = nIter, gPvalueThreshold = gPvalueThreshold, duprThreshold = duprThreshold,
      ncpus = ncpus, correct_chen_et_al_bug = correct_chen_et_al_bug
    )
  } else {
    # Windowing: dispatch by mode (C++ --wind-dist vs --wind)
    if (window_mode == "distance") {
      window_divided_res <- segment_by_dist(sum_stat$pos, max_dist = window_size, min_dim = min_dim)
    } else {
      window_divided_res <- segment_by_count(sum_stat$pos, max_count = min_dim)
    }
    dentist_result_by_window <- list()
    for (k in 1:nrow(window_divided_res)) {
      # windowEndIdx is 1-based exclusive (one past last element), so convert to
      # inclusive range by subtracting 1.
      idx_range <- window_divided_res$windowStartIdx[k]:(window_divided_res$windowEndIdx[k] - 1L)
      zScore_k <- sum_stat$z[idx_range]
      LD_mat_k <- LD_mat[idx_range, idx_range]
      dentist_result_by_window[[k]] <- dentist_single_window(
        zScore_k, R = LD_mat_k, nSample = nSample,
        pValueThreshold = pValueThreshold, propSVD = propSVD, gcControl = gcControl,
        nIter = nIter, gPvalueThreshold = gPvalueThreshold, duprThreshold = duprThreshold,
        ncpus = ncpus, correct_chen_et_al_bug = correct_chen_et_al_bug
      )
    }
    dentist_result <- merge_windows(dentist_result_by_window, window_divided_res)
  }
  return(dentist_result)
}

#' Perform DENTIST on a single window
#'
#' Detect outliers in GWAS summary statistics using LD-based iterative imputation.
#' Provide either an LD correlation matrix \code{R} or a genotype matrix \code{X}
#' (from which LD and sample size are derived automatically).
#'
#' @param zScore Numeric vector of z-scores.
#' @param R Square LD correlation matrix. Provide either \code{R} or \code{X}.
#' @param X Genotype matrix (samples x SNPs). If provided, LD is computed via
#'   \code{compute_LD(X)} and \code{nSample} defaults to \code{nrow(X)}.
#' @param nSample Number of samples in the LD reference panel (NOT the GWAS sample
#'   size). Controls the SVD truncation rank. Required when \code{R} is provided;
#'   inferred from \code{X} when \code{X} is provided.
#' @param pValueThreshold P-value threshold for outlier detection. Default is 5e-8.
#' @param propSVD SVD truncation proportion. Default is 0.4.
#' @param gcControl Logical; apply genomic control. Default is FALSE.
#' @param nIter Number of iterations. Default is 10.
#' @param gPvalueThreshold Grouping p-value threshold. Default is 0.05.
#' @param duprThreshold Duplicate r-squared threshold. Default is 0.99.
#' @param ncpus Number of CPU cores. Default is 1.
#' @param correct_chen_et_al_bug Correct the original DENTIST operator! bug. Default is TRUE.
#' @param ld_method Character string specifying the LD computation method when
#'   \code{X} is provided. Passed to \code{\link{compute_LD}}. One of
#'   \code{"sample"} (default), \code{"population"}, or \code{"gcta"}.
#'   Ignored when \code{R} is provided directly.
#'
#' @return Data frame with columns: original_z, imputed_z, iter_to_correct, rsq,
#'   is_duplicate, outlier_stat, outlier.
#'
#' @seealso \code{\link{dentist}}, \code{\link{slalom}}
#' @references \url{https://github.com/Yves-CHEN/DENTIST}
#' @export
dentist_single_window <- function(zScore, R = NULL, X = NULL, nSample = NULL,
                                  pValueThreshold = 5e-8, propSVD = 0.4, gcControl = FALSE,
                                  nIter = 10, gPvalueThreshold = 0.05, duprThreshold = 0.99,
                                  ncpus = 1, correct_chen_et_al_bug = TRUE,
                                  ld_method = "sample") {
  # Resolve LD matrix and sample size from R or X
  LD_mat <- resolve_LD_input(R = R, X = X, nSample = nSample, need_nSample = TRUE,
                             ld_method = ld_method)
  nSample <- LD_mat$nSample
  LD_mat <- LD_mat$R

  if (length(zScore) < 2000) {
    warning(sprintf(
      "The number of variants (%d) is below 2000. The algorithm may not work as expected, as suggested by the original DENTIST. Consider using window_mode = 'count' with an appropriate min_dim to control window sizes by variant count.",
      length(zScore)
    ))
  }
  if (!is.matrix(LD_mat) || nrow(LD_mat) != ncol(LD_mat) || nrow(LD_mat) != length(zScore)) {
    stop("LD_mat must be a square matrix with dimensions equal to the length of zScore.")
  }

  # Deduplicate variants
  org_Zscore <- zScore
  dedup_res <- NULL
  rThreshold <- round(sqrt(duprThreshold) * 1000) / 1000
  if (duprThreshold < 1.0) {
    dedup_res <- find_duplicate_variants(zScore, LD_mat, rThreshold)
    num_dup <- sum(dedup_res$dupBearer != -1)
    if (num_dup > 0) {
      message(paste(num_dup, "duplicated variants out of a total of", length(zScore), "were found at r threshold of", rThreshold))
    }
    zScore <- dedup_res$filteredZ
    LD_mat <- dedup_res$filteredLD
  }

  # Run C++ iterative imputation (collect rsq warnings)
  rsq_warnings <- character(0)
  warning_handler <- function(w) {
    if (grepl("Adjusted rsq_eigen value exceeding 1", w$message)) {
      rsq_warnings <<- c(rsq_warnings, w$message)
      invokeRestart("muffleWarning")
    }
  }
  verbose_iter <- getOption("pecotmr.dentist.verbose", FALSE)
  res <- withCallingHandlers(
    dentist_iterative_impute(
      LD_mat, nSample, zScore,
      pValueThreshold, propSVD, gcControl, nIter,
      gPvalueThreshold, ncpus, correct_chen_et_al_bug,
      verbose_iter
    ),
    warning = warning_handler
  )
  if (length(rsq_warnings) > 0) {
    warning(sprintf("%d rsq_eigen values exceeded 1 (capped at 1.0). Max reported: %s",
                    length(rsq_warnings), rsq_warnings[length(rsq_warnings)]))
  }
  res <- as.data.frame(res)

  # Recover duplicates
  if (duprThreshold < 1.0) {
    res <- add_dups_back_dentist(org_Zscore, res, dedup_res)
  }

  # Compute outlier stat: (z - imputed)^2 / (1 - rsq), matching binary formula
  res %>%
    mutate(
      outlier_stat = (original_z - imputed_z)^2 / pmax(1 - rsq, 1e-8),
      outlier = -log10(pchisq(outlier_stat, df = 1, lower.tail = FALSE)) > -log10(pValueThreshold)
    ) %>%
    select(-z_diff)
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
    if (dupBearer[i] == -1) {
      # Non-duplicate: copy values directly from de-duplicated output
      updated_data$imputed_z[i] <- imputed_z[assignIdx[i]]
      updated_data$z_diff[i] <- z_diff[assignIdx[i]]
      updated_data$is_duplicate[i] <- FALSE
    } else {
      # Duplicate: sign-flip imputed_z and recompute z_diff from this SNP's own z-score.
      # The original binary computes output stat as (z - imputed)^2 / (1 - rsq) using each
      # SNP's own z-score (DENTIST.h line 706), not zScore_e^2 from the bearer. We must
      # recompute z_diff here so that z_diff^2 matches the binary's stat.
      updated_data$imputed_z[i] <- imputed_z[assignIdx[i]] * sign[i]
      denom <- sqrt(max(1 - updated_data$rsq[i], 1e-8))
      updated_data$z_diff[i] <- (zScore[i] - updated_data$imputed_z[i]) / denom
      updated_data$is_duplicate[i] <- TRUE
    }
  }

  return(updated_data)
}

# ---- Segmentation helpers ----
# detect_gaps(), build_segment_result(), and sliding_window_loop() are shared
# by both segment_by_dist() and segment_by_count() to avoid code duplication.
# The core overlapping-window loop lives in sliding_window_loop(); each mode
# only supplies mode-specific callbacks for fill, step, and block-skip logic.

#' Detect Gaps in Genomic Positions
#'
#' Finds positions where the inter-SNP distance exceeds a threshold,
#' e.g., centromeric regions. Returns a vector of 1-based block boundaries.
#'
#' @param pos Sorted numeric vector of base pair positions.
#' @param gap_threshold Numeric distance threshold for gap detection.
#' @param verbose Logical; print gap info. Default is FALSE.
#'
#' @return Integer vector of 1-based block boundaries, including
#'   \code{1} (start) and \code{length(pos) + 1} (end sentinel).
#'
#' @noRd
detect_gaps <- function(pos, gap_threshold, verbose = FALSE) {
  n <- length(pos)
  diffs <- diff(pos)
  allGaps <- c(1L)
  for (i in seq_along(diffs)) {
    if (diffs[i] > gap_threshold) {
      allGaps <- c(allGaps, i + 1L)
    }
  }
  allGaps <- c(allGaps, n + 1L)

  if (verbose && length(allGaps) - 2 > 0) {
    message(sprintf("No. of gaps found: %d", length(allGaps) - 2))
    for (i in 2:(length(allGaps) - 1)) {
      message(sprintf("  Gap %d: %d - %d", i - 1, pos[allGaps[i] - 1], pos[allGaps[i]]))
    }
  }
  allGaps
}

#' Build Segment Result Data Frame
#'
#' Validates, caps indices, optionally prints verbose info, and returns the
#' standardized segmentation result data frame.
#'
#' @param startList Integer vector of window start indices.
#' @param endList Integer vector of window end indices (exclusive).
#' @param fillStartList Integer vector of fill start indices.
#' @param fillEndList Integer vector of fill end indices (exclusive).
#' @param n Total number of positions.
#' @param verbose Logical; print interval info. Default is FALSE.
#'
#' @return A data frame with columns: windowIdx, windowStartIdx, windowEndIdx,
#'   fillStartIdx, fillEndIdx.
#'
#' @noRd
build_segment_result <- function(startList, endList, fillStartList, fillEndList, n, verbose = FALSE) {
  if (length(startList) == 0) stop("No intervals created by segmentation")

  # Cap end indices at n+1 (one past the last valid 1-based index)
  endList <- pmin(endList, n + 1L)
  fillEndList <- pmin(fillEndList, n + 1L)

  if (verbose) {
    message("Intervals:")
    for (i in seq_along(startList)) {
      message(sprintf("  %d: SNPs %d-%d (fill %d-%d)",
                      i, startList[i], endList[i], fillStartList[i], fillEndList[i]))
    }
  }

  data.frame(
    windowIdx = seq_along(startList),
    windowStartIdx = startList,
    windowEndIdx = endList,
    fillStartIdx = fillStartList,
    fillEndIdx = fillEndList
  )
}

#' Sliding Window Loop for Genomic Segmentation
#'
#' Core overlapping-window loop shared by both distance-based and count-based
#' segmentation strategies. Iterates over contiguous blocks (separated by gaps),
#' creates overlapping windows within each block using mode-specific callbacks,
#' and assembles the result.
#'
#' @param allGaps Integer vector of 1-based block boundaries from
#'   \code{\link{detect_gaps}}.
#' @param n Total number of positions.
#' @param min_block_fn Function(blockSize) -> logical; returns TRUE if the block
#'   is large enough to process.
#' @param init_end_fn Function(startIdx, blockEnd) -> integer; computes the
#'   initial window end index for the first window in a block.
#' @param fill_fn Function(startIdx, endIdx, notStartInterval, notLastInterval)
#'   -> list(start, end); computes fill boundaries for each window.
#' @param step_fn Function(startIdx, blockEnd) -> list(startIdx, endIdx);
#'   advances to the next window.
#' @param adjust_last_fn Optional function(startIdx, old_startIdx, endIdx, blockEnd)
#'   -> integer; adjusts startIdx when the last interval is detected.
#'   Used by distance mode for small-last-interval correction. Default is NULL (no adjustment).
#' @param verbose Logical; print interval info. Default is FALSE.
#'
#' @return A data frame with columns: windowIdx, windowStartIdx, windowEndIdx,
#'   fillStartIdx, fillEndIdx.
#'
#' @noRd
sliding_window_loop <- function(allGaps, n,
                                min_block_fn,
                                init_end_fn,
                                fill_fn,
                                step_fn,
                                adjust_last_fn = NULL,
                                verbose = FALSE) {
  startList <- integer(0)
  endList <- integer(0)
  fillStartList <- integer(0)
  fillEndList <- integer(0)

  for (k in seq_len(length(allGaps) - 1)) {
    firstSegIdx <- length(startList) + 1
    blockStart <- allGaps[k]
    blockEnd <- allGaps[k + 1]
    blockSize <- blockEnd - blockStart

    if (!min_block_fn(blockSize)) next

    startIdx <- blockStart
    endIdx <- init_end_fn(startIdx, blockEnd)

    old_startIdx <- startIdx
    notStartInterval <- FALSE
    notLastInterval <- TRUE
    times <- 0

    repeat {
      times <- times + 1
      if (times > 400) stop("Windowing iteration limit exceeded")

      # Compute fill boundaries BEFORE any startIdx adjustment.
      # In the original C++ code, fill is recorded using the pre-adjustment
      # startIdx, then startIdx is optionally moved backward for the window.
      # This ensures fill boundaries remain non-overlapping between windows.
      fill_startIdx <- startIdx

      # Check if this is the last window
      if (blockEnd <= endIdx) {
        notLastInterval <- FALSE
        # Optional: adjust startIdx for the last window (distance mode only)
        if (!is.null(adjust_last_fn)) {
          startIdx <- adjust_last_fn(startIdx, old_startIdx, endIdx, blockEnd)
        }
      }

      # Compute fill boundaries using the pre-adjustment startIdx
      fills <- fill_fn(fill_startIdx, endIdx, notStartInterval, notLastInterval)

      startList <- c(startList, startIdx)
      endList <- c(endList, min(endIdx, blockEnd))
      fillStartList <- c(fillStartList, fills$start)
      fillEndList <- c(fillEndList, fills$end)

      if (!notLastInterval) break

      # Step to next window (mode-specific)
      old_startIdx <- startIdx
      stepped <- step_fn(startIdx, blockEnd)
      startIdx <- stepped$startIdx
      endIdx <- stepped$endIdx
      notStartInterval <- TRUE
    }

    # Fix first and last fill boundaries for this block:
    # first window's fill starts at window start, last window's fill ends at window end
    if (length(startList) >= firstSegIdx) {
      fillStartList[firstSegIdx] <- startList[firstSegIdx]
      fillEndList[length(fillEndList)] <- endList[length(endList)]
    }
  }

  build_segment_result(startList, endList, fillStartList, fillEndList, n, verbose)
}

#' Segment Genomic Region by Distance (Original DENTIST Algorithm)
#'
#' Implements the same windowing/segmentation algorithm as the original DENTIST C++ binary's
#' \code{segmentingByDist} function. Windows are created using quarter-distance SNP index
#' lookups, with gap detection for centromeres and large gaps.
#'
#' @param pos Integer vector of base pair positions (must be sorted).
#' @param max_dist Maximum distance (bp) between SNPs for windowing. Default is 2000000.
#' @param min_dim Minimum number of SNPs per window. Default is 2000.
#' @param verbose Logical; print segmentation info. Default is FALSE.
#'
#' @return A data frame with columns: windowIdx, windowStartIdx, windowEndIdx,
#'   fillStartIdx, fillEndIdx. Start indices are 1-based inclusive;
#'   end indices (windowEndIdx, fillEndIdx) are 1-based exclusive (one past last element),
#'   matching the C++ convention. Use \code{startIdx:(endIdx - 1)} for R inclusive ranges.
#'
#' @details
#' This is a faithful R translation of the C++ \code{segmentingByDist} function.
#' The algorithm:
#' \enumerate{
#'   \item Precomputes for each SNP: the index of the farthest SNP within \code{max_dist},
#'         and the index of the SNP at \code{max_dist/4} distance.
#'   \item Detects gaps > \code{max_dist/4} in the position vector (e.g., centromeres).
#'   \item Creates overlapping windows that slide by half the distance cutoff, with fill
#'         regions covering the inner three-quarters of each window.
#'   \item The first window's fill starts at the window start; the last window's fill
#'         ends at the window end.
#' }
#'
#' @seealso \code{\link{dentist_single_window}}, \code{\link{dentist}}
#'
#' @noRd
segment_by_dist <- function(pos, max_dist = 2000000, min_dim = 2000, verbose = FALSE) {
  n <- length(pos)
  if (n == 0) stop("No positions provided")

  cutoff <- max_dist
  minBlockSize <- min_dim

  # Precompute nextIdx: for each SNP i, the farthest SNP index within cutoff distance.
  # C++ uses 0-based; we translate to 1-based. Key: loop boundaries must allow
  # j to reach n+1 (one past end) so that j-1 = n (last valid 1-based index).
  nextIdx <- integer(n)
  for (i in 1:n) {
    if (i == 1) {
      j <- 2
      while (j <= n && pos[j] - pos[1] < cutoff) j <- j + 1
      nextIdx[1] <- min(j, n)
    } else {
      j <- nextIdx[i - 1]
      while (j <= n && pos[j] - pos[i] < cutoff) j <- j + 1
      nextIdx[i] <- min(j, n)
    }
  }

  # Precompute quaterIdx: for each SNP i, the last SNP index within cutoff/4 distance.
  # C++ logic: starting from the previous quaterIdx value, advance j while
  # pos[j] < cutoff/4 + pos[i], then store j-1.
  quaterIdx <- integer(n)
  # First element: find largest index where pos < cutoff/4 + pos[1]
  j <- 1
  while (j <= n && pos[j] < cutoff / 4 + as.numeric(pos[1])) j <- j + 1
  quaterIdx[1] <- max(j - 1, 1L)
  # Rest: advance from previous value
  for (i in 2:n) {
    j <- quaterIdx[i - 1]
    while (j <= n && pos[j] < cutoff / 4 + as.numeric(pos[i])) j <- j + 1
    quaterIdx[i] <- max(j - 1, 1L)
  }
  # Clamp to valid range [1, n]
  quaterIdx <- pmin(quaterIdx, n)
  quaterIdx <- pmax(quaterIdx, 1L)

  # Helper to chain quaterIdx lookups (equivalent to quaterIdx[quaterIdx[x]] in C++)
  q1 <- function(x) quaterIdx[x]
  q2 <- function(x) quaterIdx[quaterIdx[x]]
  q3 <- function(x) quaterIdx[quaterIdx[quaterIdx[x]]]
  q4 <- function(x) quaterIdx[quaterIdx[quaterIdx[quaterIdx[x]]]]

  # Find gaps > cutoff/4
  allGaps <- detect_gaps(pos, gap_threshold = cutoff / 4, verbose = verbose)

  sliding_window_loop(
    allGaps, n,
    min_block_fn = function(blockSize) {
      blockSize >= minBlockSize / 2 && (blockSize - min_dim) >= 0
    },
    init_end_fn = function(startIdx, blockEnd) {
      min(q4(startIdx) + 1, blockEnd)
    },
    fill_fn = function(startIdx, endIdx, notStartInterval, notLastInterval) {
      # Distance mode: fill is always q1 to q3 (inner 50% by distance);
      # first/last corrections are handled by fix_block_fills in the loop
      list(start = q1(startIdx), end = q3(startIdx))
    },
    step_fn = function(startIdx, blockEnd) {
      next_start <- q2(startIdx)
      list(startIdx = next_start, endIdx = min(q4(next_start) + 1, blockEnd))
    },
    adjust_last_fn = function(startIdx, old_startIdx, endIdx, blockEnd) {
      # If last interval is small, go back one step
      if (as.numeric(pos[min(endIdx - 1, n)]) - as.numeric(pos[q1(old_startIdx)]) < cutoff) {
        q1(old_startIdx)
      } else {
        startIdx
      }
    },
    verbose = verbose
  )
}

#' Segment Genomic Region by Variant Count
#'
#' Implements the windowing algorithm from the original DENTIST C++ binary's
#' \code{segmentedQCed} function. Windows contain a fixed number of variants
#' rather than spanning a fixed physical distance.
#'
#' @param pos Integer vector of base pair positions (must be sorted).
#' @param max_count Maximum number of variants per window.
#' @param gap_dist Physical distance threshold for centromeric gap detection.
#'   Default is 1e6 (matching the C++ hardcoded value).
#' @param verbose Logical; print segmentation info. Default is FALSE.
#'
#' @return A data frame with the same structure as \code{\link{segment_by_dist}}:
#'   windowIdx, windowStartIdx, windowEndIdx, fillStartIdx, fillEndIdx.
#'   End indices are 1-based exclusive (one past last element).
#'
#' @details
#' This is a faithful R translation of the C++ \code{segmentedQCed} windowing
#' algorithm. Key differences from \code{segment_by_dist}:
#' \itemize{
#'   \item Windows are sized by variant count, not physical distance.
#'   \item Uses simple index arithmetic (step = max_count/2) instead of
#'         distance-based quarter-index lookups.
#'   \item Gap detection uses a fixed 1 Mb threshold (centromeres) instead of
#'         distance/4.
#'   \item Adaptive tail absorption: if fewer than \code{max_count/2} variants
#'         remain after a window, the window extends to cover the rest.
#' }
#'
#' @seealso \code{\link{segment_by_dist}}, \code{\link{dentist}}
#'
#' @noRd
segment_by_count <- function(pos, max_count, gap_dist = 1e6, verbose = FALSE) {
  n <- length(pos)
  if (n == 0) stop("No positions provided")

  cutoff <- as.integer(max_count)
  quarter <- cutoff %/% 4L
  half <- cutoff %/% 2L

  # Detect centromeric gaps (C++ line 784: diff > 1e6)
  allGaps <- detect_gaps(pos, gap_threshold = gap_dist, verbose = verbose)

  sliding_window_loop(
    allGaps, n,
    min_block_fn = function(blockSize) blockSize >= half,
    init_end_fn = function(startIdx, blockEnd) {
      if (blockEnd - half > startIdx + cutoff) startIdx + cutoff else blockEnd
    },
    fill_fn = function(startIdx, endIdx, notStartInterval, notLastInterval) {
      # Count mode: fill based on index arithmetic (inner 50%)
      list(
        start = if (notStartInterval) startIdx + quarter else startIdx,
        end = if (notLastInterval) endIdx - quarter else endIdx
      )
    },
    step_fn = function(startIdx, blockEnd) {
      next_start <- startIdx + half
      endIdx <- if (blockEnd - half > next_start + cutoff) next_start + cutoff else blockEnd
      list(startIdx = next_start, endIdx = endIdx)
    },
    verbose = verbose
  )
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
      filter(index_global >= window_divided_res$fillStartIdx[k] & index_global < window_divided_res$fillEndIdx[k])
    merged_results <- rbind(merged_results, extracted_results)
  }
  return(merged_results)
}

### File-I/O functions (dentist_from_files, read_dentist_sumstat, parse_dentist_output)
### have been removed. Use the standard pipeline: load genotypes via
### load_genotype_region(), compute LD via compute_LD(), then call dentist()
### or ld_mismatch_qc() directly.
