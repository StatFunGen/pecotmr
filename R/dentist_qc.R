#' Resolve LD Input: Accept Either R (LD matrix) or X (Genotype Matrix)
#'
#' Internal helper that validates and resolves the LD input for QC functions.
#' Exactly one of \code{R} or \code{X} must be provided. When \code{X} is
#' provided, LD is computed via \code{computeLd(X)} and \code{nSample}
#' defaults to \code{nrow(X)}.
#'
#' @param R Square LD correlation matrix, or NULL.
#' @param X Genotype matrix (samples x SNPs), or NULL.
#' @param nSample Sample size. Required when \code{R} is provided and
#'   \code{needNSample} is TRUE; inferred from \code{X} when \code{X} is provided.
#' @param needNSample Logical; if TRUE, \code{nSample} must be available
#'   (either provided or inferred from \code{X}).
#'
#' @return A list with components \code{R} (LD correlation matrix) and
#'   \code{nSample} (integer or NULL).
#'
#' @noRd
resolveLdInput <- function(R = NULL, X = NULL, nSample = NULL, needNSample = FALSE,
                           ldMethod = "sample") {
  if (is.null(R) && is.null(X)) {
    stop("Either R (LD matrix) or X (genotype matrix) must be provided.")
  }
  if (!is.null(R) && !is.null(X)) {
    stop("Provide either R or X, not both.")
  }
  if (!is.null(X)) {
    if (!is.matrix(X)) X <- as.matrix(X)
    if (is.null(nSample)) nSample <- nrow(X)
    R <- computeLd(X, method = ldMethod)
  }
  if (needNSample && is.null(nSample)) {
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
#' @param sumStat A data frame containing summary statistics, including 'pos' or 'position' and 'z' or 'zscore' columns.
#' @param R Square LD correlation matrix. Provide either \code{R} or \code{X}.
#' @param X Genotype matrix (samples x SNPs). If provided, LD is computed via
#'   \code{computeLd(X)} and \code{nSample} defaults to \code{nrow(X)}.
#' @param nSample The number of samples in the LD reference panel (NOT the GWAS sample
#'   size). This controls the SVD truncation rank K = min(idx_size, nSample) * propSVD.
#'   Required when \code{R} is provided; inferred from \code{X} when \code{X} is provided.
#' @param windowSize The size of the window for dividing the genomic region
#'   in distance mode (base pairs). Default is 2000000 (2 Mb). Only used when
#'   \code{windowMode = "distance"}.
#' @param windowMode Character string specifying the windowing strategy:
#'   \code{"distance"} (default) creates windows by physical distance using
#'   \code{\link{segmentByDist}} (C++ \code{--wind-dist}), and
#'   \code{"count"} creates windows by variant count using
#'   \code{\link{segmentByCount}} (C++ \code{--wind}).
#' @param pValueThreshold The p-value threshold for significance. Default is 5e-8.
#' @param propSVD The proportion of singular value decomposition (SVD) to use. Default is 0.4.
#' @param gcControl Logical indicating whether genomic control should be applied. Default is FALSE.
#' @param nIter The number of iterations for the Dentist algorithm. Default is 10.
#' @param gPvalueThreshold The genomic p-value threshold for significance. Default is 0.05.
#' @param duprThreshold The absolute correlation r value threshold to be considered duplicate. Default is 0.99.
#' @param ncpus The number of CPU cores to use for parallel processing. Default is 1.
#' @param correctChenEtAlBug Logical indicating whether to correct the Chen et al. bug. Default is TRUE.
#' @param minDim In distance mode: minimum number of SNPs per block (default 2000).
#'   In count mode: the number of variants per window (i.e., the window size).
#' @param ldMethod Character string specifying the LD computation method when
#'   \code{X} is provided. Passed to \code{\link{computeLd}}. One of
#'   \code{"sample"} (default), \code{"population"}, or \code{"gcta"}.
#'   Ignored when \code{R} is provided directly.
#'
#' @return A data frame containing the imputed result and detected outliers.
#'
#' The returned data frame includes the following columns:
#'
#' \describe{
#'   \item{\code{original_z}}{The original z-score values from the input \code{sumStat}.}
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
#' dentist(sumStat, R = ldMat, nSample = nSample)
#'
#' @details
#' Windowing supports two modes matching the original DENTIST C++ binary:
#' \itemize{
#'   \item \code{"distance"} (default): Uses the \code{segmentingByDist} algorithm
#'     (C++ \code{--wind-dist}), implemented in \code{\link{segmentByDist}}.
#'     Windows span a fixed physical distance (\code{windowSize} bp).
#'   \item \code{"count"}: Uses the \code{segmentedQCed} algorithm
#'     (C++ \code{--wind}), implemented in \code{\link{segmentByCount}}.
#'     Windows contain a fixed number of variants (\code{minDim}).
#'     Useful when regions have sparse variants where distance-based windows
#'     would create windows with too few variants.
#' }
#' The \code{correctChenEtAlBug} parameter affects the iterative filtering
#' in two ways:
#' \enumerate{
#'   \item Comparison between iteration index \code{t} and \code{nIter} (explained in source code)
#'   \item The \code{!grouping_tmp} operator bug (explained in source code)
#' }
#'
#' @export
dentist <- function(sumStat, R = NULL, X = NULL, nSample = NULL,
                    windowSize = 2000000, windowMode = c("distance", "count"),
                    pValueThreshold = 5.0369e-8, propSVD = 0.4, gcControl = FALSE,
                    nIter = 10, gPvalueThreshold = 0.05, duprThreshold = 0.99, ncpus = 1,
                    correctChenEtAlBug = TRUE, minDim = 2000,
                    ldMethod = "sample") {
  # Resolve LD matrix and sample size from R or X
  resolved <- resolveLdInput(R = R, X = X, nSample = nSample, needNSample = TRUE,
                             ldMethod = ldMethod)
  ldMat <- resolved$R
  nSample <- resolved$nSample

  # detect for column names and order by pos
  if (!any(tolower(c("pos", "position")) %in% tolower(colnames(sumStat))) ||
    !any(tolower(c("z", "zscore")) %in% tolower(colnames(sumStat)))) {
    stop("Input sumStat is missing either 'pos'/'position' or 'z'/'zscore' column.")
  }
  # rename to common column name
  if (!tolower("pos") %in% tolower(colnames(sumStat))) {
    colnames(sumStat)[which(tolower(colnames(sumStat)) %in% tolower(c("position")))] <- "pos"
  }

  if (!tolower("z") %in% tolower(colnames(sumStat))) {
    colnames(sumStat)[which(tolower(colnames(sumStat)) %in% tolower(c("zscore")))] <- "z"
  }

  sumStat <- sumStat %>% arrange(pos)

  windowMode <- match.arg(windowMode)

  # If the data has fewer SNPs than minDim, run as a single window directly.
  nSnps <- nrow(sumStat)
  if (nSnps < minDim) {
    dentistResult <- dentistSingleWindow(
      sumStat$z, R = ldMat, nSample = nSample,
      pValueThreshold = pValueThreshold, propSVD = propSVD, gcControl = gcControl,
      nIter = nIter, gPvalueThreshold = gPvalueThreshold, duprThreshold = duprThreshold,
      ncpus = ncpus, correctChenEtAlBug = correctChenEtAlBug
    )
  } else {
    # Windowing: dispatch by mode (C++ --wind-dist vs --wind)
    if (windowMode == "distance") {
      windowDividedRes <- segmentByDist(sumStat$pos, maxDist = windowSize, minDim = minDim)
    } else {
      windowDividedRes <- segmentByCount(sumStat$pos, maxCount = minDim)
    }
    dentistResultByWindow <- list()
    for (k in 1:nrow(windowDividedRes)) {
      # windowEndIdx is 1-based exclusive (one past last element), so convert to
      # inclusive range by subtracting 1.
      idxRange <- windowDividedRes$windowStartIdx[k]:(windowDividedRes$windowEndIdx[k] - 1L)
      zScoreK <- sumStat$z[idxRange]
      ldMatK <- ldMat[idxRange, idxRange]
      dentistResultByWindow[[k]] <- dentistSingleWindow(
        zScoreK, R = ldMatK, nSample = nSample,
        pValueThreshold = pValueThreshold, propSVD = propSVD, gcControl = gcControl,
        nIter = nIter, gPvalueThreshold = gPvalueThreshold, duprThreshold = duprThreshold,
        ncpus = ncpus, correctChenEtAlBug = correctChenEtAlBug
      )
    }
    dentistResult <- mergeWindows(dentistResultByWindow, windowDividedRes)
  }
  return(dentistResult)
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
#'   \code{computeLd(X)} and \code{nSample} defaults to \code{nrow(X)}.
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
#' @param correctChenEtAlBug Correct the original DENTIST operator! bug. Default is TRUE.
#' @param ldMethod Character string specifying the LD computation method when
#'   \code{X} is provided. Passed to \code{\link{computeLd}}. One of
#'   \code{"sample"} (default), \code{"population"}, or \code{"gcta"}.
#'   Ignored when \code{R} is provided directly.
#'
#' @return Data frame with columns: original_z, imputed_z, iter_to_correct, rsq,
#'   is_duplicate, outlier_stat, outlier.
#'
#' @seealso \code{\link{dentist}}, \code{\link{slalom}}
#' @references \url{https://github.com/Yves-CHEN/DENTIST}
#' @export
dentistSingleWindow <- function(zScore, R = NULL, X = NULL, nSample = NULL,
                                pValueThreshold = 5e-8, propSVD = 0.4, gcControl = FALSE,
                                nIter = 10, gPvalueThreshold = 0.05, duprThreshold = 0.99,
                                ncpus = 1, correctChenEtAlBug = TRUE,
                                ldMethod = "sample") {
  # Resolve LD matrix and sample size from R or X
  ldMat <- resolveLdInput(R = R, X = X, nSample = nSample, needNSample = TRUE,
                          ldMethod = ldMethod)
  nSample <- ldMat$nSample
  ldMat <- ldMat$R

  if (length(zScore) < 2000) {
    warning(sprintf(
      "The number of variants (%d) is below 2000. The algorithm may not work as expected, as suggested by the original DENTIST. Consider using windowMode = 'count' with an appropriate minDim to control window sizes by variant count.",
      length(zScore)
    ))
  }
  if (!is.matrix(ldMat) || nrow(ldMat) != ncol(ldMat) || nrow(ldMat) != length(zScore)) {
    stop("ldMat must be a square matrix with dimensions equal to the length of zScore.")
  }

  # Deduplicate variants
  orgZscore <- zScore
  dedupRes <- NULL
  rThreshold <- round(sqrt(duprThreshold) * 1000) / 1000
  if (duprThreshold < 1.0) {
    dedupRes <- findDuplicateVariants(zScore, ldMat, rThreshold)
    numDup <- sum(dedupRes$dupBearer != -1)
    if (numDup > 0) {
      message(paste(numDup, "duplicated variants out of a total of", length(zScore), "were found at r threshold of", rThreshold))
    }
    zScore <- dedupRes$filteredZ
    ldMat <- dedupRes$filteredLD
  }

  # Run C++ iterative imputation (collect rsq warnings)
  rsqWarnings <- character(0)
  warningHandler <- function(w) {
    if (grepl("Adjusted rsq_eigen value exceeding 1", w$message)) {
      rsqWarnings <<- c(rsqWarnings, w$message)
      invokeRestart("muffleWarning")
    }
  }
  verboseIter <- getOption("pecotmr.dentist.verbose", FALSE)
  res <- withCallingHandlers(
    # cpp11 requires exact integer types for int parameters
    dentistIterativeImpute(
      ldMat, as.integer(nSample), zScore,
      pValueThreshold, propSVD, gcControl, as.integer(nIter),
      gPvalueThreshold, as.integer(ncpus), correctChenEtAlBug,
      verboseIter
    ),
    warning = warningHandler
  )
  if (length(rsqWarnings) > 0) {
    warning(sprintf("%d rsq_eigen values exceeded 1 (capped at 1.0). Max reported: %s",
                    length(rsqWarnings), rsqWarnings[length(rsqWarnings)]))
  }
  res <- as.data.frame(res)
  # cpp11 wrapper returns camelCase keys; convert to documented snake_case columns
  names(res)[names(res) == "originalZ"] <- "original_z"
  names(res)[names(res) == "imputedZ"] <- "imputed_z"
  names(res)[names(res) == "zDiff"] <- "z_diff"
  names(res)[names(res) == "iterToCorrect"] <- "iter_to_correct"

  # Recover duplicates
  if (duprThreshold < 1.0) {
    res <- addDupsBackDentist(orgZscore, res, dedupRes)
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
#' based on the output from the `findDuplicateVariants` function.
#' @param zScore The original zScore
#' @param dentistOutput A data frame containing the output from the DENTIST algorithm.
#' @param findDupOutput A list containing the output from the `findDuplicateVariants` function.
#'
#' @return A data frame with duplicated variants added back and an additional column indicating duplicates.
#'
#' @noRd
addDupsBackDentist <- function(zScore, dentistOutput, findDupOutput) {
  # Extract relevant columns from the DENTIST output
  originalZ <- dentistOutput$original_z
  imputedZ <- dentistOutput$imputed_z
  iterToCorrect <- dentistOutput$iter_to_correct
  rsq <- dentistOutput$rsq
  zDiff <- dentistOutput$z_diff

  # Extract output from findDuplicateVariants
  dupBearer <- findDupOutput$dupBearer
  sign <- findDupOutput$sign

  # Get the number of rows in dupBearer
  nrowsDup <- length(dupBearer)

  if (nrow(dentistOutput) != sum(dupBearer == -1)) {
    stop("The number of rows in the input data does not match the occurrences of -1 in dupBearer.")
  }

  if (length(zScore) != nrowsDup) {
    stop("Input zScore and findDupOutput have inconsistent dimension")
  }

  # Initialize assignIdx vector
  count <- 1
  assignIdx <- rep(0, nrowsDup)

  for (i in seq_along(dupBearer)) {
    if (dupBearer[i] == -1) {
      assignIdx[i] <- count
      count <- count + 1
    } else {
      assignIdx[i] <- dupBearer[i]
    }
  }

  # Create a new data frame to store the updated values
  updatedData <- data.frame(
    original_z = numeric(nrowsDup),
    imputed_z = numeric(nrowsDup),
    iter_to_correct = numeric(nrowsDup),
    rsq = numeric(nrowsDup),
    z_diff = numeric(nrowsDup),
    is_duplicate = logical(nrowsDup)
  )

  for (i in seq_len(nrowsDup)) {
    updatedData$original_z[i] <- zScore[i]
    updatedData$iter_to_correct[i] <- iterToCorrect[assignIdx[i]]
    updatedData$rsq[i] <- rsq[assignIdx[i]]
    if (dupBearer[i] == -1) {
      # Non-duplicate: copy values directly from de-duplicated output
      updatedData$imputed_z[i] <- imputedZ[assignIdx[i]]
      updatedData$z_diff[i] <- zDiff[assignIdx[i]]
      updatedData$is_duplicate[i] <- FALSE
    } else {
      # Duplicate: sign-flip imputed_z and recompute z_diff from this SNP's own z-score.
      # The original binary computes output stat as (z - imputed)^2 / (1 - rsq) using each
      # SNP's own z-score (DENTIST.h line 706), not zScore_e^2 from the bearer. We must
      # recompute z_diff here so that z_diff^2 matches the binary's stat.
      updatedData$imputed_z[i] <- imputedZ[assignIdx[i]] * sign[i]
      denom <- sqrt(max(1 - updatedData$rsq[i], 1e-8))
      updatedData$z_diff[i] <- (zScore[i] - updatedData$imputed_z[i]) / denom
      updatedData$is_duplicate[i] <- TRUE
    }
  }

  return(updatedData)
}

# ---- Segmentation helpers ----
# detectGaps(), buildSegmentResult(), and slidingWindowLoop() are shared
# by both segmentByDist() and segmentByCount() to avoid code duplication.
# The core overlapping-window loop lives in slidingWindowLoop(); each mode
# only supplies mode-specific callbacks for fill, step, and block-skip logic.

#' Detect Gaps in Genomic Positions
#'
#' Finds positions where the inter-SNP distance exceeds a threshold,
#' e.g., centromeric regions. Returns a vector of 1-based block boundaries.
#'
#' @param pos Sorted numeric vector of base pair positions.
#' @param gapThreshold Numeric distance threshold for gap detection.
#' @param verbose Logical; print gap info. Default is FALSE.
#'
#' @return Integer vector of 1-based block boundaries, including
#'   \code{1} (start) and \code{length(pos) + 1} (end sentinel).
#'
#' @noRd
detectGaps <- function(pos, gapThreshold, verbose = FALSE) {
  n <- length(pos)
  diffs <- diff(pos)
  allGaps <- c(1L)
  for (i in seq_along(diffs)) {
    if (diffs[i] > gapThreshold) {
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
buildSegmentResult <- function(startList, endList, fillStartList, fillEndList, n, verbose = FALSE) {
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
#'   \code{\link{detectGaps}}.
#' @param n Total number of positions.
#' @param minBlockFn Function(blockSize) -> logical; returns TRUE if the block
#'   is large enough to process.
#' @param initEndFn Function(startIdx, blockEnd) -> integer; computes the
#'   initial window end index for the first window in a block.
#' @param fillFn Function(startIdx, endIdx, notStartInterval, notLastInterval)
#'   -> list(start, end); computes fill boundaries for each window.
#' @param stepFn Function(startIdx, blockEnd) -> list(startIdx, endIdx);
#'   advances to the next window.
#' @param adjustLastFn Optional function(startIdx, oldStartIdx, endIdx, blockEnd)
#'   -> integer; adjusts startIdx when the last interval is detected.
#'   Used by distance mode for small-last-interval correction. Default is NULL (no adjustment).
#' @param verbose Logical; print interval info. Default is FALSE.
#'
#' @return A data frame with columns: windowIdx, windowStartIdx, windowEndIdx,
#'   fillStartIdx, fillEndIdx.
#'
#' @noRd
slidingWindowLoop <- function(allGaps, n,
                              minBlockFn,
                              initEndFn,
                              fillFn,
                              stepFn,
                              adjustLastFn = NULL,
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

    if (!minBlockFn(blockSize)) next

    startIdx <- blockStart
    endIdx <- initEndFn(startIdx, blockEnd)

    oldStartIdx <- startIdx
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
      fillStartIdx <- startIdx

      # Check if this is the last window
      if (blockEnd <= endIdx) {
        notLastInterval <- FALSE
        # Optional: adjust startIdx for the last window (distance mode only)
        if (!is.null(adjustLastFn)) {
          startIdx <- adjustLastFn(startIdx, oldStartIdx, endIdx, blockEnd)
        }
      }

      # Compute fill boundaries using the pre-adjustment startIdx
      fills <- fillFn(fillStartIdx, endIdx, notStartInterval, notLastInterval)

      startList <- c(startList, startIdx)
      endList <- c(endList, min(endIdx, blockEnd))
      fillStartList <- c(fillStartList, fills$start)
      fillEndList <- c(fillEndList, fills$end)

      if (!notLastInterval) break

      # Step to next window (mode-specific)
      oldStartIdx <- startIdx
      stepped <- stepFn(startIdx, blockEnd)
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

  buildSegmentResult(startList, endList, fillStartList, fillEndList, n, verbose)
}

#' Segment Genomic Region by Distance (Original DENTIST Algorithm)
#'
#' Implements the same windowing/segmentation algorithm as the original DENTIST C++ binary's
#' \code{segmentingByDist} function. Windows are created using quarter-distance SNP index
#' lookups, with gap detection for centromeres and large gaps.
#'
#' @param pos Integer vector of base pair positions (must be sorted).
#' @param maxDist Maximum distance (bp) between SNPs for windowing. Default is 2000000.
#' @param minDim Minimum number of SNPs per window. Default is 2000.
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
#'   \item Precomputes for each SNP: the index of the farthest SNP within \code{maxDist},
#'         and the index of the SNP at \code{maxDist/4} distance.
#'   \item Detects gaps > \code{maxDist/4} in the position vector (e.g., centromeres).
#'   \item Creates overlapping windows that slide by half the distance cutoff, with fill
#'         regions covering the inner three-quarters of each window.
#'   \item The first window's fill starts at the window start; the last window's fill
#'         ends at the window end.
#' }
#'
#' @seealso \code{\link{dentistSingleWindow}}, \code{\link{dentist}}
#'
#' @noRd
segmentByDist <- function(pos, maxDist = 2000000, minDim = 2000, verbose = FALSE) {
  n <- length(pos)
  if (n == 0) stop("No positions provided")

  cutoff <- maxDist
  minBlockSize <- minDim

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
  allGaps <- detectGaps(pos, gapThreshold = cutoff / 4, verbose = verbose)

  slidingWindowLoop(
    allGaps, n,
    minBlockFn = function(blockSize) {
      blockSize >= minBlockSize / 2 && (blockSize - minDim) >= 0
    },
    initEndFn = function(startIdx, blockEnd) {
      min(q4(startIdx) + 1, blockEnd)
    },
    fillFn = function(startIdx, endIdx, notStartInterval, notLastInterval) {
      # Distance mode: fill is always q1 to q3 (inner 50% by distance);
      # first/last corrections are handled by fix_block_fills in the loop
      list(start = q1(startIdx), end = q3(startIdx))
    },
    stepFn = function(startIdx, blockEnd) {
      nextStart <- q2(startIdx)
      list(startIdx = nextStart, endIdx = min(q4(nextStart) + 1, blockEnd))
    },
    adjustLastFn = function(startIdx, oldStartIdx, endIdx, blockEnd) {
      # If last interval is small, go back one step
      if (as.numeric(pos[min(endIdx - 1, n)]) - as.numeric(pos[q1(oldStartIdx)]) < cutoff) {
        q1(oldStartIdx)
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
#' @param maxCount Maximum number of variants per window.
#' @param gapDist Physical distance threshold for centromeric gap detection.
#'   Default is 1e6 (matching the C++ hardcoded value).
#' @param verbose Logical; print segmentation info. Default is FALSE.
#'
#' @return A data frame with the same structure as \code{\link{segmentByDist}}:
#'   windowIdx, windowStartIdx, windowEndIdx, fillStartIdx, fillEndIdx.
#'   End indices are 1-based exclusive (one past last element).
#'
#' @details
#' This is a faithful R translation of the C++ \code{segmentedQCed} windowing
#' algorithm. Key differences from \code{segmentByDist}:
#' \itemize{
#'   \item Windows are sized by variant count, not physical distance.
#'   \item Uses simple index arithmetic (step = maxCount/2) instead of
#'         distance-based quarter-index lookups.
#'   \item Gap detection uses a fixed 1 Mb threshold (centromeres) instead of
#'         distance/4.
#'   \item Adaptive tail absorption: if fewer than \code{maxCount/2} variants
#'         remain after a window, the window extends to cover the rest.
#' }
#'
#' @seealso \code{\link{segmentByDist}}, \code{\link{dentist}}
#'
#' @noRd
segmentByCount <- function(pos, maxCount, gapDist = 1e6, verbose = FALSE) {
  n <- length(pos)
  if (n == 0) stop("No positions provided")

  cutoff <- as.integer(maxCount)
  quarter <- cutoff %/% 4L
  half <- cutoff %/% 2L

  # Detect centromeric gaps (C++ line 784: diff > 1e6)
  allGaps <- detectGaps(pos, gapThreshold = gapDist, verbose = verbose)

  slidingWindowLoop(
    allGaps, n,
    minBlockFn = function(blockSize) blockSize >= half,
    initEndFn = function(startIdx, blockEnd) {
      if (blockEnd - half > startIdx + cutoff) startIdx + cutoff else blockEnd
    },
    fillFn = function(startIdx, endIdx, notStartInterval, notLastInterval) {
      # Count mode: fill based on index arithmetic (inner 50%)
      list(
        start = if (notStartInterval) startIdx + quarter else startIdx,
        end = if (notLastInterval) endIdx - quarter else endIdx
      )
    },
    stepFn = function(startIdx, blockEnd) {
      nextStart <- startIdx + half
      endIdx <- if (blockEnd - half > nextStart + cutoff) nextStart + cutoff else blockEnd
      list(startIdx = nextStart, endIdx = endIdx)
    },
    verbose = verbose
  )
}

#' Merge dentist Results by Window
#'
#' This function merges DENTIST results by window into a single data frame.
#'
#' @param dentistResultByWindow A list containing imputed results for each window.
#' @param windowDividedRes A data frame containing information about the divided windows.
#'
#' @return A data frame containing merged results.
#'
#' @details
#' The function checks if the number of imputed results matches the number of windows.
#' It then merges the results by window, adding an index within the window and a global index.
#' Finally, it extracts the results within the fillers and combines them into a single data frame.
#'
#' @noRd
mergeWindows <- function(dentistResultByWindow, windowDividedRes) {
  if (length(dentistResultByWindow) != nrow(windowDividedRes)) {
    stop("Different number of windows and imputed results!")
  }
  mergedResults <- c()
  for (k in 1:nrow(windowDividedRes)) {
    imputedK <- dentistResultByWindow[[k]]
    imputedK$index_within_window <- seq(1:nrow(imputedK))
    imputedK <- imputedK %>%
      mutate(index_global = index_within_window + windowDividedRes$windowStartIdx[k] - 1)
    extractedResults <- imputedK %>%
      filter(index_global >= windowDividedRes$fillStartIdx[k] & index_global < windowDividedRes$fillEndIdx[k])
    mergedResults <- rbind(mergedResults, extractedResults)
  }
  return(mergedResults)
}

### File-I/O functions (dentist_from_files, read_dentist_sumstat, parse_dentist_output)
### have been removed. Use the standard pipeline: load genotypes via
### loadGenotypeRegion(), compute LD via computeLd(), then call dentist()
### or ldMismatchQc() directly.
