# =============================================================================
# Summary-statistic QC pipeline
# -----------------------------------------------------------------------------
# Consolidated QC suite for GwasSumStats / QtlSumStats objects. The
# top-level entry point is `summaryStatsQc()` (below), which orchestrates
# the individual passes:
#
#   * Allele harmonization (harmonizeAlleles, in R/variantId.R)
#   * RAISS sumstats imputation (fills variants present on the LD panel
#     but missing from sumstats)
#   * SLALoM (single-causal-variant ABF outlier detection)
#   * DENTIST (test-based LD-mismatch detection)
#   * Univariate RSS diagnostics (post-finemap mismatch diagnostics)
#
# Each pass occupies its own section below; the orchestrator lives at the
# bottom. Pure individual-level sample QC (relatedness etc.) is in
# R/relatednessQc.R, not here.
# =============================================================================

#' @importFrom GenomicRanges seqnames
#' @importFrom S4Vectors mcols
NULL

# Allele harmonization (harmonizeAlleles) now lives in R/variantId.R, alongside
# the other variant-matching primitives (parseVariantId, matchVariants). It is
# package-internal and called from here (.matchAgainstSketch), ctwasPipeline,
# and the pipeline join sites.

# Standardize a variant set (GRanges or data.frame) to a chrom/pos/alt/ref frame.
# @noRd
.variantsToDf <- function(x) {
  if (is(x, "GRanges")) {
    mc <- as.data.frame(mcols(x))
    mc$chrom <- as.character(seqnames(x))
    mc$pos <- start(x)
    mc[, c("chrom", "pos", "alt", "ref")]
  } else {
    as.data.frame(x)[, c("chrom", "pos", "alt", "ref")]
  }
}

# Canonical per-variant key from sorted alleles, so strand/allele flips collide.
# @noRd
.canonicalAlleleKey <- function(df) {
  aMin <- pmin(df$alt, df$ref)
  aMax <- pmax(df$alt, df$ref)
  paste(df$chrom, df$pos, aMin, aMax)
}

#' Merge variant info from two sources with allele-flip-aware matching
#'
#' Merges variant metadata (chromosome, position, ref, alt) from two sources,
#' detecting and correcting allele flips (where alt/ref are swapped). Creates
#' a canonical key from sorted alleles to match across datasets.
#'
#' @param variants1 A data.frame with columns \code{chrom}, \code{pos},
#'   \code{alt}, \code{ref}, or a \code{GRanges} with corresponding metadata
#'   columns.
#' @param variants2 A data.frame or \code{GRanges} with the same columns.
#' @param all Logical. If TRUE (default), returns the union of both sets.
#'   If FALSE, returns only variants from \code{variants2} (flipped to match
#'   \code{variants1}'s allele orientation).
#' @return A data.frame with columns \code{chrom}, \code{pos}, \code{alt},
#'   \code{ref}, deduplicated by position and alleles.
#' @export
mergeVariantInfo <- function(variants1, variants2, all = TRUE) {
  df1 <- .variantsToDf(variants1)
  df2 <- .variantsToDf(variants2)

  key1 <- .canonicalAlleleKey(df1)
  key2 <- .canonicalAlleleKey(df2)

  # Detect flips: where df2's alt matches df1's ref at the same key
  matchIdx <- match(key2, key1)
  hasMatch <- !is.na(matchIdx)

  flip <- rep(FALSE, nrow(df2))
  mi <- matchIdx[hasMatch]
  flip[hasMatch] <- df2$alt[hasMatch] == df1$ref[mi] &
                    df2$ref[hasMatch] == df1$alt[mi]

  # Apply flips to df2
  flipRows <- which(hasMatch)[flip[hasMatch]]
  if (length(flipRows) > 0) {
    tmp <- df2$alt[flipRows]
    df2$alt[flipRows] <- df2$ref[flipRows]
    df2$ref[flipRows] <- tmp
  }

  if (all) {
    combined <- rbind(
      df1[, c("chrom", "pos", "alt", "ref")],
      df2[, c("chrom", "pos", "alt", "ref")])
    combined[!duplicated(paste(combined$chrom, combined$pos,
                               combined$alt, combined$ref)), ]
  } else {
    df2[, c("chrom", "pos", "alt", "ref")]
  }
}



# =============================================================================
# DENTIST: deterministic test-based LD-mismatch detection
# =============================================================================

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

# Module-scoped accumulator for the rsq_eigen warnings raised (per call) by the
# C++ dentistIterativeImpute(); dentistSingleWindow resets $msgs before the call
# and reads it back after. dentist windows are processed serially, so the shared
# env carries no cross-call state.
.dentistRsqAcc <- new.env(parent = emptyenv())

# withCallingHandlers warning handler: collect the "rsq_eigen exceeding 1"
# warnings into .dentistRsqAcc and muffle them (they are summarized afterwards).
# @noRd
.dentistRsqHandler <- function(w) {
  if (grepl("Adjusted rsq_eigen value exceeding 1", w$message)) {
    .dentistRsqAcc$msgs <- c(.dentistRsqAcc$msgs, w$message)
    invokeRestart("muffleWarning")
  }
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
    dedupRes <- .findDuplicateVariants(zScore, ldMat, rThreshold)
    numDup <- sum(dedupRes$dupBearer != -1)
    if (numDup > 0) {
      message(paste(numDup, "duplicated variants out of a total of", length(zScore), "were found at r threshold of", rThreshold))
    }
    zScore <- dedupRes$filteredZ
    ldMat <- dedupRes$filteredLD
  }

  # Run C++ iterative imputation; rsq_eigen warnings are collected into the
  # module-level .dentistRsqAcc (reset here, read back after the call).
  .dentistRsqAcc$msgs <- character(0)
  verboseIter <- getOption("pecotmr.dentist.verbose", FALSE)
  res <- withCallingHandlers(
    # cpp11 requires exact integer types for int parameters
    dentistIterativeImpute(
      ldMat, as.integer(nSample), zScore,
      pValueThreshold, propSVD, gcControl, as.integer(nIter),
      gPvalueThreshold, as.integer(ncpus), correctChenEtAlBug,
      verboseIter
    ),
    warning = .dentistRsqHandler
  )
  rsqWarnings <- .dentistRsqAcc$msgs
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

# Apply the quarter-distance index map `quaterIdx` `n` times to `x` (n = 1..4
# gives the 1st..4th quarter boundary from x). Used by segmentByDist.
# @noRd
.nthQuaterIdx <- function(x, n, quaterIdx) {
  for (i in seq_len(n)) x <- quaterIdx[x]
  x
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

  # Find gaps > cutoff/4
  allGaps <- detectGaps(pos, gapThreshold = cutoff / 4, verbose = verbose)

  slidingWindowLoop(
    allGaps, n,
    minBlockFn = function(blockSize) {
      blockSize >= minBlockSize / 2 && (blockSize - minDim) >= 0
    },
    initEndFn = function(startIdx, blockEnd) {
      min(.nthQuaterIdx(startIdx, 4, quaterIdx) + 1, blockEnd)
    },
    fillFn = function(startIdx, endIdx, notStartInterval, notLastInterval) {
      # Distance mode: fill is always q1 to q3 (inner 50% by distance);
      # first/last corrections are handled by fix_block_fills in the loop
      list(start = .nthQuaterIdx(startIdx, 1, quaterIdx), end = .nthQuaterIdx(startIdx, 3, quaterIdx))
    },
    stepFn = function(startIdx, blockEnd) {
      nextStart <- .nthQuaterIdx(startIdx, 2, quaterIdx)
      list(startIdx = nextStart, endIdx = min(.nthQuaterIdx(nextStart, 4, quaterIdx) + 1, blockEnd))
    },
    adjustLastFn = function(startIdx, oldStartIdx, endIdx, blockEnd) {
      # If last interval is small, go back one step
      if (as.numeric(pos[min(endIdx - 1, n)]) - as.numeric(pos[.nthQuaterIdx(oldStartIdx, 1, quaterIdx)]) < cutoff) {
        .nthQuaterIdx(oldStartIdx, 1, quaterIdx)
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


# =============================================================================
# SLALoM: Approximate Bayes Factor single-causal-variant outlier detection
# =============================================================================

# Numerically stable log-sum-exp.
# @noRd
.logSumExp <- function(x) {
  maxX <- max(x, na.rm = TRUE)
  sumExp <- sum(exp(x - maxX), na.rm = TRUE)
  return(maxX + log(sumExp))
}

# Approximate (Wakefield) Bayes factors from z-scores and standard errors:
# per-variant log-BF and normalized posterior probabilities.
# @noRd
.approxBayesFactor <- function(z, se, W = 0.04) {
  V <- se^2
  r <- W / (W + V)
  lbf <- 0.5 * (log(1 - r) + (r * z^2))
  denom <- .logSumExp(lbf)
  prob <- exp(lbf - denom)
  return(list(lbf = lbf, prob = prob))
}

# Greedy credible set: indices (ordered by decreasing prob) whose cumulative
# posterior first exceeds `coverage`.
# @noRd
.slalomCredibleSet <- function(prob, coverage = 0.95) {
  ordering <- order(prob, decreasing = TRUE)
  cumprob <- cumsum(prob[ordering])
  idx <- which(cumprob > coverage)[1]
  cs <- ordering[1:idx]
  return(cs)
}

#' Slalom Function for Summary Statistics QC for Fine-Mapping Analysis
#'
#' Performs Approximate Bayesian Factor (ABF) analysis, identifies credible sets,
#' and annotates lead variants based on fine-mapping results. It computes p-values
#' from z-scores assuming a two-sided standard normal distribution.
#'
#' Provide either an LD correlation matrix \code{R} or a genotype matrix \code{X}
#' (from which LD is derived automatically via \code{computeLd}).
#'
#' @param zScore Numeric vector of z-scores corresponding to each variant.
#' @param R Square LD correlation matrix. Provide either \code{R} or \code{X}.
#' @param X Genotype matrix (samples x SNPs). If provided, LD is computed via
#'   \code{computeLd(X)}.
#' @param standardError Optional numeric vector of standard errors corresponding
#'   to each z-score. If not provided, a default value of 1 is assumed for all variants.
#' @param abfPriorVariance Numeric, the prior effect size variance for ABF calculations.
#'   Default is 0.04.
#' @param nlog10pDentistSThreshold Numeric, the -log10 DENTIST-S P value threshold
#'   for identifying outlier variants for prediction. Default is 4.0.
#' @param r2Threshold Numeric, the r2 threshold for DENTIST-S outlier variants
#'   for prediction. Default is 0.6.
#' @param leadVariantChoice Character, method to choose the lead variant, either
#'   "pvalue" or "abf", with default "pvalue".
#' @param ldMethod Character string specifying the LD computation method when
#'   \code{X} is provided. Passed to \code{\link{computeLd}}. One of
#'   \code{"sample"} (default), \code{"population"}, or \code{"gcta"}.
#'   Ignored when \code{R} is provided directly.
#' @return A list containing the annotated LD matrix with ABF results, credible sets,
#'   lead variant, and DENTIST-S statistics; and a summary dataframe with aggregate statistics.
#' @examples
#' results <- slalom(zScore, R = R, standardError = standardError)
#' @seealso \code{\link{dentistSingleWindow}}, \code{\link{resolveLdInput}}
#' @export
#'
slalom <- function(zScore, R = NULL, X = NULL, standardError = rep(1, length(zScore)),
                   abfPriorVariance = 0.04, nlog10pDentistSThreshold = 4.0,
                   r2Threshold = 0.6, leadVariantChoice = "pvalue",
                   ldMethod = "sample") {
  if (is.null(R) && is.null(X)) {
    stop("Either R (LD matrix) or X (genotype matrix) must be provided.")
  }
  if (!is.null(R) && !is.null(X)) {
    stop("Provide either R or X, not both.")
  }

  # One-sided p-value matching the original Python implementation (stats.norm.cdf).
  # This selects the most negative z-score as lead when leadVariantChoice == "pvalue".
  pvalue <- pnorm(zScore)

  abfResults <- .approxBayesFactor(zScore, standardError, W = abfPriorVariance)
  lbf <- abfResults$lbf
  prob <- abfResults$prob

  cs <- .slalomCredibleSet(prob, coverage = 0.95)
  cs99 <- .slalomCredibleSet(prob, coverage = 0.99)

  leadIdx <- if (leadVariantChoice == "pvalue") {
    which.min(pvalue)
  } else {
    which.max(prob)
  }

  # Only the lead column of R is needed for DENTIST-S.
  # When X is provided, compute just that column instead of the full p x p matrix.
  if (!is.null(X)) {
    if (!is.matrix(X)) X <- as.matrix(X)
    rLead <- as.numeric(cor(X, X[, leadIdx]))
  } else {
    if (!is.matrix(R) || nrow(R) != ncol(R) || nrow(R) != length(zScore)) {
      stop("R must be a square matrix matching the length of zScore.")
    }
    rLead <- R[, leadIdx]
  }

  r2Lead <- rLead^2
  tDentistS <- (zScore - rLead * zScore[leadIdx])^2 / (1 - r2Lead)
  tDentistS[tDentistS < 0] <- Inf
  nlog10pDentistS <- -log10(pchisq(tDentistS, df = 1, lower.tail = FALSE))
  outliers <- (r2Lead > r2Threshold) & (nlog10pDentistS > nlog10pDentistSThreshold)

  nR2 <- sum(r2Lead > r2Threshold)
  nDentistSOutlier <- sum(outliers, na.rm = TRUE)
  maxPip <- max(prob)

  summary <- list(
    leadPipVariant = leadIdx,
    nTotal = length(zScore),
    nR2 = nR2,
    nDentistSOutlier = nDentistSOutlier,
    fraction = ifelse(nR2 > 0, nDentistSOutlier / nR2, 0),
    maxPip = maxPip,
    cs95 = cs,
    cs99 = cs99
  )
  result <- as.data.frame(list(original_z = zScore, prob = prob, pvalue = pvalue, outliers = outliers, nlog10p_dentist_s = nlog10pDentistS))

  return(list(data = result, summary = summary))
}


# =============================================================================
# Univariate RSS diagnostics (post-finemap)
# =============================================================================

#' Extract the trimmed SuSiE fit from a finemapping pipeline result
#'
#' Returns the trimmed model fit underlying \code{con_data$finemappingEntry}
#' (a \code{FineMappingEntry} S4 object), or NULL if no fine-mapping entry
#' is attached.
#'
#' @param conData List. The method-layer entry from a finemapping pipeline
#'   result, expected to carry \code{$finemappingEntry} as a
#'   \code{FineMappingEntry} object.
#' @return The trimmed fit (a list with \code{pip}, \code{sets}, etc.) or NULL.
#' @export
getSusieResult <- function(conData) {
  if (length(conData) == 0) return(NULL)
  fm <- conData$finemappingEntry
  if (is.null(fm) || !is(fm, "FineMappingEntry")) return(NULL)
  trimmed <- getSusieFit(fm)
  if (length(trimmed) == 0) return(NULL)
  trimmed
}

#' Process Credible Sets (CS) from Finemapping Results
#'
#' This function extracts and processes information for each Credible Set (CS) 
#' from finemapping results, typically obtained from a finemapping RDS file.
#'
#' @param conData List. The method layer data from a finemapping RDS file that is not empty.
#' @param csNames Character vector. Names of the Credible Sets, usually in the format "L_<number>".
#' @param topLociTable Data frame. The $top_loci layer data from the finemapping results.
#'
#' @return A data frame with one row per CS, containing the following columns:
#'   \item{cs_name}{Name of the Credible Set}
#'   \item{variants_per_cs}{Number of variants in the CS}
#'   \item{top_variant}{ID of the variant with the highest PIP in the CS}
#'   \item{top_variant_index}{Global index of the top variant}
#'   \item{top_pip}{Highest Posterior Inclusion Probability (PIP) in the CS}
#'   \item{top_z}{Z-score of the top variant}
#'   \item{p_value}{P-value calculated from the top Z-score}
#'   \item{cs_corr}{Pairwise correlations of other CSs in this RDS with the CS of 
#'     the current row, delimited by '|', if there is more than one CS in this RDS file}
#'
#' @details
#' This function is designed to be used only when there is at least one Credible Set 
#' in the finemapping results usually for a given study and block. It processes each CS, 
#' extracting key information such as the top variant, its statistics, and 
#' correlation information between multiple CS if available.
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#'
#' @export
extractCsInfo <- function(conData, csNames, topLociTable) {
  fm <- conData$finemappingEntry
  trimmed <- getSusieFit(fm)
  variantNames <- getVariantIds(fm)
  results <- map(seq_along(csNames), function(i) {
    csName <- csNames[i]
    indices <- trimmed$sets$cs[[csName]]

    # Get variants for this CS using the full variant names list
    csVariants <- variantNames[indices]
    csData <- topLociTable[topLociTable$variant_id %in% csVariants, ]
    topRow <- which.max(csData$pip)

    topVariant <- csData$variant_id[topRow]
    # Find the global index of the top variant
    topVariantGlobalIndex <- which(variantNames == topVariant)
    topPip <- csData$pip[topRow]
    topZ <- csData$z[topRow]
    pValue <- .zToPvalue(topZ)

    # Extract cs_corr
    csCorr <- if (length(csNames) > 1) {
      trimmed$cs_corr[i,]
    } else {
      NA  # Use NA for the second CS or when there's only one CS
    }

    # Return results for this CS as a one-row data.frame
    result <- tibble(
      cs_name = csName,
      variants_per_cs = length(csVariants),
      top_variant = topVariant,
      top_variant_index = topVariantGlobalIndex,
      top_pip = topPip,
      top_z = topZ,
      p_value = pValue,
      cs_corr = list(paste(csCorr, collapse = ","))  # list column if csCorr is a vector
    )
    return(result)
  })
  # Combine all tibbles into one data frame
  finalResult <- bind_rows(results)
  return(finalResult)
}

#' Extract Information for Top Variant from Finemapping Results
#'
#' This function extracts information about the variant with the highest Posterior 
#' Inclusion Probability (PIP) from finemapping results, typically used when no 
#' Credible Sets (CS) are identified in the analysis.
#'
#' @param conData List. The method layer data from a finemapping RDS file.
#'
#' @return A data frame with one row containing the following columns:
#'   \item{cs_name}{NA (as no CS is identified)}
#'   \item{variants_per_cs}{NA (as no CS is identified)}
#'   \item{top_variant}{ID of the variant with the highest PIP}
#'   \item{top_variant_index}{Index of the top variant in the original data}
#'   \item{top_pip}{Highest Posterior Inclusion Probability (PIP)}
#'   \item{top_z}{Z-score of the top variant}
#'   \item{p_value}{P-value calculated from the top Z-score}
#'   \item{cs_corr}{NA (as no CS correlation is available)}
#'
#' @details
#' This function is designed to be used when no Credible Sets are identified in 
#' the finemapping results, but information about the most significant variant 
#' is still desired. It identifies the variant with the highest PIP and extracts 
#' relevant statistical information.
#'
#' @note
#' This function is particularly useful for capturing information about potentially 
#' important variants that might be included in Credible Sets under different 
#' analysis parameters or lower coverage. It maintains a structure similar to 
#' the output of `extract_cs_info()` for consistency in downstream analyses.
#'
#' @seealso
#' \code{\link{extractCsInfo}} for processing when Credible Sets are present.
#'
#' @export
extractTopPipInfo <- function(conData) {
  fm <- conData$finemappingEntry
  trimmed <- getSusieFit(fm)
  variantNames <- getVariantIds(fm)
  # Find the variant with the highest PIP
  topPipIndex <- which.max(trimmed$pip)
  topPip <- trimmed$pip[topPipIndex]
  topVariant <- variantNames[topPipIndex]
  topZ <- conData$sumstats$z[topPipIndex]
  pValue <- .zToPvalue(topZ)

  list(
    cs_name = NA,
    variants_per_cs = NA,
    top_variant = topVariant,
    top_variant_index = topPipIndex,
    top_pip = topPip,
    top_z = topZ,
    p_value = pValue,
    cs_corr = NA  # or NULL
  )
}

# Parse one comma-joined cs_corr string into its values + max/min |corr| (self-
# correlations equal to 1 excluded); invalid/empty input yields empties + NA.
# @noRd
.extractCorrelations <- function(x) {
  # Early return if x is invalid
  if(is.na(x) || x == "" || is.null(x) || !grepl(",", as.character(x))) {
    return(list(values = numeric(0), max_corr = NA_real_, min_corr = NA_real_))
  }

  # Convert and filter values
  values <- as.numeric(unlist(strsplit(x, ",")))
  valuesFiltered <- abs(values[values != 1])

  # Return list with NA if no valid correlations
  list(
    values = values,
    max_corr = if(length(valuesFiltered) > 0) max(abs(valuesFiltered), na.rm = TRUE) else NA_real_,
    min_corr = if(length(valuesFiltered) > 0) min(abs(valuesFiltered), na.rm = TRUE) else NA_real_
  )
}

#' Parse Credible Set Correlations from extractCsInfo() Output
#'
#' This function takes the output from `extractCsInfo()` and expands the `cs_corr` column
#' into multiple columns, preserving the original order of correlations. It also
#' calculates maximum and minimum correlation values for each Credible Set.
#'
#' @param df Data frame. The output from `extractCsInfo()` function,
#'           containing a `cs_corr` column with correlation information.
#'
#' @return A data frame with the original columns from the input, plus:
#'   \item{cs_corr_1, cs_corr_2, ...}{Individual correlation values, with column names
#'         based on their position in the original string}
#'   \item{cs_corr_max}{Maximum absolute correlation value (excluding 1)}
#'   \item{cs_corr_min}{Minimum absolute correlation value}
#'
#' @details
#' The function splits the `cs_corr` column, which typically contains correlation
#' values separated by '|', into individual columns. It preserves the order of
#' these correlations, allowing for easy interpretation in a matrix-like format.
#'
#' @note
#' - This function converts the input to a data frame if it isn't already one.
#' - It handles cases where correlation values might be missing or not in the expected format.
#' - The function assumes that correlation values of 1 represent self-correlations and excludes
#'   these when calculating max and min correlations.
#'
#' @export
parseCsCorr <- function(df) {
  # Ensure we work with a data frame
  df <- as.data.frame(df)

  # Process correlations
  processedResults <- lapply(df$cs_corr, .extractCorrelations)
  # If no valid results, add NA columns and return
  if(all(sapply(processedResults, function(x) length(x$values) == 0))) {
    df$cs_corr_max <- NA_real_
    df$cs_corr_min <- NA_real_
    return(df)
  }

  # Determine max number of correlations
  maxCorrCount <- max(sapply(processedResults, function(x) length(x$values)))

  # Create and add correlation columns
  colNames <- paste0("cs_corr_", 1:maxCorrCount)

  for(i in seq_along(colNames)) {
    df[[colNames[i]]] <- sapply(processedResults, function(x) {
      if(length(x$values) >= i) x$values[i] else NA_real_
    })
  }

  # Add max and min correlation columns
  df$cs_corr_max <- sapply(processedResults, `[[`, "max_corr")
  df$cs_corr_min <- sapply(processedResults, `[[`, "min_corr")

  return(df)
}

#' Process Credible Set Information and Determine Updating Strategy
#'
#' This function categorizes Credible Sets (CS) within a study block into different 
#' updating strategies based on their statistical properties and correlations.
#'
#' @param df Data frame. Contains information about Credible Sets for a specific study and block.
#' @param highCorrCols Character vector. Names of columns in df that represent high correlations.
#'
#' @return A modified data frame with additional columns attached to the diagnostic table:
#'   \item{top_cs}{Logical. TRUE for the CS with the highest absolute Z-score.}
#'   \item{tagged_cs}{Logical. TRUE for CS that are considered "tagged" based on p-value and correlation criteria.}
#'   \item{method}{Character. The determined updating strategy ("BVSR", "SER", or "BCR").}
#'
#' @details
#' This function performs the following steps:
#' 1. Identifies the top CS based on the highest absolute Z-score.
#' 2. Identifies tagged CS based on high p-value and high correlations.
#' 3. Counts total, tagged, and remaining CS.
#' 4. Determines the appropriate updating method based on these counts.
#'
#' The updating methods are:
#' - BVSR (Bayesian Variable Selection Regression): Used when there's only one CS or all CS are accounted for.
#' - SER (Single Effect Regression): Used when there are tagged CS but no remaining untagged CS.
#' - BCR (Bayesian Conditional Regression): Used when there are remaining untagged CS.
#'
#' @note
#' This function is part of a developing methodology for automatically handling 
#' finemapping results. The thresholds and criteria used (e.g., p-value > 1e-4 for tagging) 
#' are subject to refinement and may change in future versions.
#'
#' @importFrom dplyr case_when
#'
#' @export
autoDecision <- function(df, highCorrCols) {
  # Identify top_cs
  topCsIndex <- which.max(abs(df$top_z))
  df$top_cs <- FALSE
  df$top_cs[topCsIndex] <- TRUE

  # Identify tagged_cs
  df$tagged_cs <- sapply(1:nrow(df), function(i) {
    if (df$top_cs[i]) return(FALSE)
    if (df$p_value[i] > 1e-4) return(TRUE)
    if (length(highCorrCols) == 0) return(FALSE)
    any(sapply(highCorrCols, function(col) df[i, ..col] == 1))
  })

  # Count total and remaining CS
  totalCs <- nrow(df)
  print("total_cs")
  print(totalCs)
  taggedCsCount <- sum(df$tagged_cs)
  if (totalCs > 0) {
    remainingCs <- totalCs - 1 - taggedCsCount
  } else {
    remainingCs <- 0
  }
  # Determine method
  df$method <- case_when(
  taggedCsCount == 0 & totalCs > 1 ~ "BVSR",
  (remainingCs == 0 & totalCs > 1) | (totalCs == 1) ~ "SER",
  remainingCs > 0 ~ "BCR",
  TRUE ~ NA_character_
)


  return(df)
}



# =============================================================================
# RAISS: regression-based sumstats imputation
# =============================================================================

#' Core RAISS implementation for a single LD matrix
#'
#' @param refPanel A data frame containing 'chrom', 'pos', 'variant_id', 'A1', and 'A2'.
#' @param knownZscores A data frame containing 'chrom', 'pos', 'variant_id', 'A1', 'A2', and 'z' values.
#' @param ldMatrix A square matrix of dimension equal to the number of rows in refPanel.
#' @param lamb Regularization term added to the diagonal of the ldMatrix.
#' @param rcond Threshold for filtering eigenvalues in the pseudo-inverse computation.
#' @param r2Threshold R square threshold below which SNPs are filtered from the output.
#' @param minimumLd Minimum LD score threshold for SNP filtering.
#' @param verbose Logical indicating whether to print progress information.
#'
#' @return A list containing filtered and unfiltered results, and filtered LD matrix.
#' @importFrom MASS ginv
#' @importFrom dplyr arrange
#' @noRd
raissSingleMatrix <- function(refPanel, knownZscores, ldMatrix, lamb = 0.01, rcond = 0.01,
                              r2Threshold = 0.6, minimumLd = 5, verbose = TRUE) {
  # Check that refPanel and knownZscores are both increasing in terms of pos
  if (is.unsorted(refPanel$pos) || is.unsorted(knownZscores$pos)) {
    stop("refPanel and knownZscores must be in increasing order of pos.")
  }

  # Convert ldMatrix to matrix if it's a data frame
  if (is.data.frame(ldMatrix)) {
    ldMatrix <- as.matrix(ldMatrix)
  }

  # Define knowns and unknowns
  knownsId <- intersect(knownZscores$variant_id, refPanel$variant_id)
  knowns <- which(refPanel$variant_id %in% knownsId)
  unknowns <- which(!refPanel$variant_id %in% knownsId)

  # Handle edge cases
  if (length(knowns) == 0) {
    if (verbose) message("No known variants found, cannot perform imputation.")
    return(NULL)
  }

  if (length(unknowns) == 0) {
    if (verbose) message("No unknown variants to impute, returning known variants.")
    return(list(
      resultNofilter = knownZscores,
      resultFilter = knownZscores,
      ldMat = ldMatrix
    ))
  }

  # Extract zt, sigT, and sigIT
  zt <- knownZscores$z
  sigT <- ldMatrix[knowns, knowns, drop = FALSE]
  sigIT <- ldMatrix[unknowns, knowns, drop = FALSE]

  # Call raissModel
  results <- raissModel(zt, sigT, sigIT, lamb, rcond)
  # Format the results
  results <- formatRaissDf(results, refPanel, unknowns)
  # Filter output
  results <- filterRaissOutput(results, r2Threshold, minimumLd, verbose)

  # Merge with known z-scores
  resultNofilter <- mergeRaissDf(results$zscoresNofilter, knownZscores) %>% arrange(pos)
  resultFilter <- mergeRaissDf(results$zscores, knownZscores) %>% arrange(pos)

  # Filter out variants not included in the imputation result
  filteredOutVariant <- setdiff(refPanel$variant_id, resultFilter$variant_id)

  # Update the LD matrix excluding filtered variants
  ldExtractFiltered <- if (length(filteredOutVariant) > 0) {
    filteredOutId <- match(filteredOutVariant, refPanel$variant_id)
    as.matrix(ldMatrix)[-filteredOutId, -filteredOutId]
  } else {
    as.matrix(ldMatrix)
  }
  # Return results
  return(list(
    resultNofilter = resultNofilter,
    resultFilter = resultFilter,
    ldMat = ldExtractFiltered
  ))
}

#' Core RAISS implementation from a genotype matrix X (SVD-based)
#'
#' Performs the same imputation as \code{raissSingleMatrix} but works directly
#' with the genotype matrix X instead of the LD correlation matrix R. This avoids
#' forming the p x p LD matrix, saving O(p^2) memory and O(np^2) compute.
#'
#' The reformulation is mathematically exact: using the thin SVD of Xt (the known
#' variant columns), all RAISS quantities (mu, var, ld_score) are computed in the
#' SVD basis without ever forming R = X'X/(n-1).
#'
#' @param refPanel A data frame containing 'chrom', 'pos', 'variant_id', 'A1', and 'A2'.
#' @param knownZscores A data frame containing 'chrom', 'pos', 'variant_id', 'A1', 'A2', and 'z' values.
#' @param X Centered and scaled genotype matrix (nSamples x pVariants). Column order must
#'   match the variant order in refPanel.
#' @param lamb Regularization term (same role as in the LD-based path).
#' @param svdTol Relative tolerance for filtering small singular values in the SVD of Xt.
#' @param r2Threshold R square threshold below which SNPs are filtered from the output.
#' @param minimumLd Minimum LD score threshold for SNP filtering.
#' @param verbose Logical indicating whether to print progress information.
#'
#' @return A list containing filtered and unfiltered results, and ldMat = NULL.
#' @importFrom dplyr arrange
#' @noRd
raissSingleMatrixFromX <- function(refPanel, knownZscores, X, lamb = 0.01,
                                   svdTol = 1e-8, r2Threshold = 0.6,
                                   minimumLd = 5, verbose = TRUE) {
  # Check that refPanel and knownZscores are both increasing in terms of pos
  if (is.unsorted(refPanel$pos) || is.unsorted(knownZscores$pos)) {
    stop("refPanel and knownZscores must be in increasing order of pos.")
  }

  nSamples <- nrow(X)

  # Define knowns and unknowns (same logic as raissSingleMatrix)
  knownsId <- intersect(knownZscores$variant_id, refPanel$variant_id)
  knowns <- which(refPanel$variant_id %in% knownsId)
  unknowns <- which(!refPanel$variant_id %in% knownsId)

  # Handle edge cases
  if (length(knowns) == 0) {
    if (verbose) message("No known variants found, cannot perform imputation.")
    return(NULL)
  }

  if (length(unknowns) == 0) {
    if (verbose) message("No unknown variants to impute, returning known variants.")
    return(list(
      resultNofilter = knownZscores,
      resultFilter = knownZscores,
      ldMat = NULL
    ))
  }

  # Extract known columns for SVD (unavoidable copy for LAPACK).
  # We do NOT copy X_i - instead we compute X' %*% [w|U] on the full X
  # and index the unknown rows, saving O(n*m) memory.
  Xt <- X[, knowns, drop = FALSE]
  zt <- knownZscores$z

  # Compute thin SVD of Xt (n x k -> U: n x r, d: r, V: k x r)
  svdResult <- .safeSvd(Xt, tol = svdTol)
  U <- svdResult$u
  d <- svdResult$d
  V <- svdResult$v
  rm(Xt)  # free n*k memory; no longer needed

  # Precompute regularization and weight vectors (length r, cheap)
  cReg <- lamb * (nSamples - 1)
  d2 <- d^2
  d2PlusC <- d2 + cReg

  # --- Build w (n x 1): the projection of zt through the regularized SVD ---
  # w = U %*% diag(d / (d^2 + c)) %*% V' zt
  VtZt <- crossprod(V, zt)                       # r x 1
  w <- U %*% (d / d2PlusC * VtZt)                # n x 1

  # --- Single BLAS call: X' %*% [w | U] -> p x (1+r) ---
  # This avoids copying X_i (n x m) entirely.
  # Row unknowns of column 1 gives mu; rows unknowns of columns 2:(r+1) gives A.
  XtWU <- crossprod(X, cbind(w, U))               # p x (1+r), one dgemm call
  mu <- as.numeric(XtWU[unknowns, 1])             # m x 1
  A <- XtWU[unknowns, -1, drop = FALSE]           # m x r (subset, not copy of X)
  rm(XtWU)                                         # free p*(1+r)

  # --- Variance and LD score in one pass over A^2 ---
  # var = (1+lamb) - (1/(n-1)) * A^2 %*% (d^2/(d^2+c))
  # ld_score = (1/(n-1))^2 * A^2 %*% d^2
  # Compute A^2 once, multiply by [d_weights_var | d^2] in one dgemm.
  ASq <- A^2                                       # m x r (one allocation)
  rm(A)                                            # free m*r
  dWeights <- cbind(d2 / d2PlusC, d2)              # r x 2
  scores <- ASq %*% dWeights                       # m x 2 (one dgemm)
  rm(ASq)                                          # free m*r

  nm1 <- nSamples - 1
  varRaw <- (1 + lamb) - scores[, 1] / nm1
  raissLdScore <- scores[, 2] / nm1^2
  rm(scores)

  # --- Condition number (scalar, expanded to vector by formatRaissDf) ---
  conditionNumber <- rep(d[1] / d[length(d)], length(unknowns))
  correctInversion <- rep(TRUE, length(unknowns))

  # --- R2 correction (same as raissModel) ---
  varNorm <- varInBoundaries(varRaw, lamb)
  R2 <- (1 + lamb) - varNorm
  mu <- mu / sqrt(R2)

  # Package results in the same format as raissModel output
  imp <- list(
    var = varNorm,
    mu = mu,
    raissLdScore = raissLdScore,
    conditionNumber = conditionNumber,
    correctInversion = correctInversion
  )

  # Reuse existing formatting and filtering functions
  results <- formatRaissDf(imp, refPanel, unknowns)
  results <- filterRaissOutput(results, r2Threshold, minimumLd, verbose)

  # Merge with known z-scores
  resultNofilter <- mergeRaissDf(results$zscoresNofilter, knownZscores) %>% arrange(pos)
  resultFilter <- mergeRaissDf(results$zscores, knownZscores) %>% arrange(pos)

  return(list(
    resultNofilter = resultNofilter,
    resultFilter = resultFilter,
    ldMat = NULL
  ))
}

# Sequentially append RAISS block results, resolving the shared boundary variant
# (duplicated across adjacent blocks) by keeping the higher-R2 imputation.
# @noRd
.combineWithBoundaryCheck <- function(combinedResult, newResult) {
  # If either is empty, simply return the non-empty one or empty data frame
  if (is.null(combinedResult)) {
    return(newResult)
  }
  if (is.null(newResult)) {
    return(combinedResult)
  }

  # Check if the last variant of combined matches the first of new
  lastVar <- combinedResult$variant_id[nrow(combinedResult)]
  firstVar <- newResult$variant_id[1]

  if (lastVar == firstVar) {
    newR2 <- newResult$raissR2[1]
    oldR2 <- combinedResult$raissR2[nrow(combinedResult)]
    if (is.na(newR2) && is.na(oldR2)) {
      # Both are NA - keep the existing one
    } else if (is.na(oldR2)) {
      # Old is NA but new is not - use new
      combinedResult[nrow(combinedResult), ] <- newResult[1, ]
    } else if (is.na(newR2)) {
      # New is NA but old is not - keep old
    } else if (newR2 > oldR2) {
      # Both are non-NA and new is better - use new
      combinedResult[nrow(combinedResult), ] <- newResult[1, ]
    }

    # Add remaining rows from new (excluding first)
    if (nrow(newResult) > 1) {
      combinedResult <- bind_rows(combinedResult, newResult[-1, ])
    }
  } else {
    # No overlap - combine all rows
    combinedResult <- bind_rows(combinedResult, newResult)
  }

  return(combinedResult)
}

#' Impute Summary Statistics Using LD (RAISS)
#'
#' This function is a part of the statistical library for SNP imputation from:
#' https://gitlab.pasteur.fr/statistical-genetics/raiss/-/blob/master/raiss/stat_models.py
#' It is R implementation of the imputation model described in the paper by Bogdan Pasaniuc,
#' Noah Zaitlen, et al., titled "Fast and accurate imputation of summary
#' statistics enhances evidence of functional enrichment", published in
#' Bioinformatics in 2014.
#'
#' This function can process either a single LD matrix or a list of LD matrices for different blocks.
#' For a list of matrices, it processes each block separately and combines the results.
#' Alternatively, it can accept a genotype matrix X directly, avoiding the need to form
#' the p x p LD matrix (memory and compute savings when n << p).
#'
#' @param refPanel A data frame containing 'chrom', 'pos', 'variant_id', 'A1', and 'A2'.
#' @param knownZscores A data frame containing 'chrom', 'pos', 'variant_id', 'A1', 'A2', and 'z' values.
#' @param ldMatrix Either a square matrix or a list of matrices for LD blocks.
#'   Provide either \code{ldMatrix} or \code{genotypeMatrix}, not both.
#' @param genotypeMatrix A centered and scaled genotype matrix (n x p) as an alternative
#'   to \code{ldMatrix}. Column order must match the variant order in \code{refPanel}.
#'   When provided, the imputation uses an SVD-based approach that avoids forming the
#'   p x p LD matrix.
#' @param lamb Regularization term added to the diagonal of the ldMatrix.
#' @param rcond Threshold for filtering eigenvalues in the pseudo-inverse computation
#'   (only used with ldMatrix path).
#' @param svdTol Relative tolerance for filtering small singular values
#'   (only used with genotypeMatrix path).
#' @param r2Threshold R square threshold below which SNPs are filtered from the output.
#' @param minimumLd Minimum LD score threshold for SNP filtering.
#' @param verbose Logical indicating whether to print progress information.
#'
#' @return A list containing filtered and unfiltered results, and filtered LD matrix
#'   (ldMat is NULL when using genotypeMatrix path).
#' @importFrom dplyr arrange bind_rows
#' @export
raiss <- function(refPanel, knownZscores, ldMatrix = NULL,
                  genotypeMatrix = NULL, lamb = 0.01, rcond = 0.01,
                  svdTol = 1e-8, r2Threshold = 0.6, minimumLd = 5,
                  verbose = TRUE) {
  # --- Genotype matrix path (SVD-based, avoids forming R) ---
  if (!is.null(genotypeMatrix)) {
    if (!is.null(ldMatrix)) {
      stop("Provide either ldMatrix or genotypeMatrix, not both.")
    }
    if (is.matrix(genotypeMatrix)) {
      if (verbose) message("Processing genotype matrix via SVD-based imputation...")
      return(raissSingleMatrixFromX(
        refPanel, knownZscores, genotypeMatrix,
        lamb, svdTol, r2Threshold, minimumLd, verbose
      ))
    }
    if (is.list(genotypeMatrix)) {
      # List of genotype matrices (block processing)
      if (verbose) message("Processing multiple genotype matrix blocks via SVD-based imputation...")
      resultsList <- list()
      for (i in seq_along(genotypeMatrix)) {
        if (verbose) message(paste("Processing block", i, "of", length(genotypeMatrix)))
        blockResult <- raissSingleMatrixFromX(
          refPanel, knownZscores, genotypeMatrix[[i]],
          lamb, svdTol, r2Threshold, minimumLd,
          verbose = FALSE
        )
        if (!is.null(blockResult)) {
          resultsList[[length(resultsList) + 1]] <- blockResult
        }
      }
      if (length(resultsList) == 0) {
        if (verbose) message("No blocks could be processed.")
        return(NULL)
      }
      combinedNofilter <- do.call(bind_rows, lapply(resultsList, `[[`, "resultNofilter"))
      combinedFilter <- do.call(bind_rows, lapply(resultsList, `[[`, "resultFilter"))
      return(list(
        resultNofilter = combinedNofilter %>% arrange(pos),
        resultFilter = combinedFilter %>% arrange(pos),
        ldMat = NULL
      ))
    }
    stop("genotypeMatrix must be a matrix or a list of matrices.")
  }

  # --- LD matrix path (original implementation) ---
  if (is.null(ldMatrix)) {
    stop("Provide either ldMatrix or genotypeMatrix.")
  }
  # Determine if we can process as a single matrix
  isSingleMatrixCase <- is.matrix(ldMatrix) ||
    (is.list(ldMatrix) && !is.null(ldMatrix$ldMatrices) &&
      length(ldMatrix$ldMatrices) == 1)

  if (isSingleMatrixCase) {
    if (verbose) message("Processing single LD matrix", if (!is.matrix(ldMatrix)) " from list", "...")

    # Extract the matrix if it's in a list
    if (!is.matrix(ldMatrix)) {
      ldMatrix <- ldMatrix$ldMatrices[[1]]
    }

    return(raissSingleMatrix(
      refPanel, knownZscores, ldMatrix,
      lamb, rcond, r2Threshold, minimumLd, verbose
    ))
  }

  # For list of matrices, process each block
  if (verbose) message("Processing multiple LD blocks...")

  resultsList <- list()
  variantIndices <- ldMatrix$variantIndices
  blockIds <- unique(variantIndices$blockId)

  for (blockId in blockIds) {
    if (verbose) message(paste("Processing block", blockId, "of", length(blockIds)))

    blockVariantIds <- variantIndices$variant_id[variantIndices$blockId == blockId]

    # Subset refPanel and ldMatrix for this block
    blockIndices <- match(blockVariantIds, refPanel$variant_id)
    blockRefPanel <- refPanel[blockIndices, ]
    blockLdMatrix <- ldMatrix$ldMatrices[[blockId]]
    blockKnownZscores <- knownZscores %>% filter(variant_id %in% blockVariantIds)
    if (nrow(blockLdMatrix) != nrow(blockRefPanel)) {
      stop(paste("Block", blockId, ": LD matrix dimension does not match number of variants in reference panel"))
    }

    # Process the block using the core function
    blockResult <- raissSingleMatrix(
      blockRefPanel, blockKnownZscores, blockLdMatrix,
      lamb, rcond, r2Threshold, minimumLd,
      verbose = FALSE
    )
    # Skip if block returned NULL (no known variants)
    if (!is.null(blockResult)) {
      resultsList[[blockId]] <- blockResult
    }
  }

  if (length(resultsList) == 0) {
    if (verbose) message("No blocks could be processed. Check that knownZscores overlap with variants in the blocks.")
    return(NULL)
  }

  # Combine results sequentially to handle boundary duplicates
  combinedNofilter <- resultsList[[1]]$resultNofilter
  combinedFilter <- resultsList[[1]]$resultFilter

  if (length(resultsList) > 1) {
    for (i in 2:length(resultsList)) {
      combinedNofilter <- .combineWithBoundaryCheck(
        combinedNofilter,
        resultsList[[i]]$resultNofilter
      )

      combinedFilter <- .combineWithBoundaryCheck(
        combinedFilter,
        resultsList[[i]]$resultFilter
      )
    }
  }

  ldFilteredList <- lapply(resultsList, function(x) x$ldMat)
  variantList <- lapply(ldFilteredList, function(ld) data.frame(variants = colnames(ld)))
  ldMatrix <- createLdMatrix(
    ldMatrices = ldFilteredList,
    variants = variantList
  )

  return(list(
    resultNofilter = combinedNofilter,
    resultFilter = combinedFilter,
    ldMat = ldMatrix
  ))
}

#' @param zt Vector of known z scores.
#' @param sigT Matrix of known linkage disequilibrium (LD) correlation.
#' @param sigIT Correlation matrix with rows corresponding to unknown SNPs (to impute)
#'               and columns to known SNPs.
#' @param lamb Regularization term added to the diagonal of the sigT matrix.
#' @param rcond Threshold for filtering eigenvalues in the pseudo-inverse computation.
#' @param batch Boolean indicating whether batch processing is used.
#'
#' @return A list containing the variance 'var', estimation 'mu', LD score 'raissLdScore',
#'         condition number 'conditionNumber', and correctness of inversion
#'         'correctInversion'.
#' @noRd
raissModel <- function(zt, sigT, sigIT, lamb = 0.01, rcond = 0.01, batch = TRUE, reportConditionNumber = FALSE) {
  sigTInv <- invertMatRecursive(sigT, lamb, rcond)
  if (!is.numeric(zt) || !is.numeric(sigT) || !is.numeric(sigIT)) {
    stop("zt, sigT, and sigIT must be numeric.")
  }
  if (batch) {
    conditionNumber <- if (reportConditionNumber) rep(kappa(sigT, exact = TRUE, norm = "2"), nrow(sigIT)) else NA
    correctInversion <- rep(checkInversion(sigT, sigTInv), nrow(sigIT))
  } else {
    conditionNumber <- if (reportConditionNumber) kappa(sigT, exact = TRUE, norm = "2") else NA
    correctInversion <- checkInversion(sigT, sigTInv)
  }

  varRaissLdScore <- computeVar(sigIT, sigTInv, lamb, batch)
  var <- varRaissLdScore$var
  raissLdScore <- varRaissLdScore$raissLdScore

  mu <- computeMu(sigIT, sigTInv, zt)
  varNorm <- varInBoundaries(var, lamb)

  R2 <- ((1 + lamb) - varNorm)
  mu <- mu / sqrt(R2)

  return(list(var = varNorm, mu = mu, raissLdScore = raissLdScore, conditionNumber = conditionNumber, correctInversion = correctInversion))
}

#' @param imp is the output of raissModel()
#' @param refPanel is a data frame with columns 'chrom', 'pos', 'variant_id', 'ref', and 'alt'.
#' @noRd
formatRaissDf <- function(imp, refPanel, unknowns) {
  resultDf <- data.frame(
    chrom = refPanel[unknowns, "chrom"],
    pos = refPanel[unknowns, "pos"],
    variant_id = refPanel[unknowns, "variant_id"],
    A1 = refPanel[unknowns, "A1"],
    A2 = refPanel[unknowns, "A2"],
    z = imp$mu,
    Var = imp$var,
    raissLdScore = imp$raissLdScore,
    conditionNumber = imp$conditionNumber,
    correctInversion = imp$correctInversion
  )

  # Specify the column order
  columnOrder <- c(
    "chrom", "pos", "variant_id", "A1", "A2", "z", "Var", "raissLdScore", "conditionNumber",
    "correctInversion"
  )

  # Reorder the columns
  resultDf <- resultDf[, columnOrder]
  return(resultDf)
}

mergeRaissDf <- function(raissDf, knownZscores) {
  # Merge the data frames
  mergedDf <- merge(raissDf, knownZscores, by = c("chrom", "pos", "variant_id", "A1", "A2"), all = TRUE)

  # Identify rows that came from knownZscores
  fromKnown <- !is.na(mergedDf$z.y) & is.na(mergedDf$z.x)

  # Set Var to -1 and raissLdScore to Inf for these rows
  mergedDf$Var[fromKnown] <- -1
  mergedDf$raissLdScore[fromKnown] <- Inf

  # If there are overlapping columns (e.g., z.x and z.y), resolve them
  # For example, use z from knownZscores where available, otherwise use z from raissDf
  mergedDf$z <- ifelse(fromKnown, mergedDf$z.y, mergedDf$z.x)

  # Remove the extra columns resulted from the merge (e.g., z.x, z.y)
  mergedDf <- mergedDf[, !colnames(mergedDf) %in% c("z.x", "z.y")]
  mergedDf <- arrange(mergedDf, pos)
  # assign imputed variants beta, se as NA to avoid confusion, since they are not imputed
  mergedDf$beta[mergedDf$Var == -1] <- NA
  mergedDf$se[mergedDf$Var == -1] <- NA
  return(mergedDf)
}

# Format one aligned "label: value" report line (label left-padded to
# maxLabelLength).
# @noRd
.formatRaissLine <- function(label, value, maxLabelLength) {
  sprintf("%-*s %d", maxLabelLength, paste0(label, ":"), value)
}

filterRaissOutput <- function(zscores, r2Threshold = 0.6, minimumLd = 5, verbose = TRUE) {
  # Reset the index and subset the data frame
  zscores <- zscores[, c("chrom", "pos", "variant_id", "A1", "A2", "z", "Var", "raissLdScore")]
  zscores$raissR2 <- 1 - zscores$Var

  # Count statistics before filtering
  nSnpsBfFilt <- nrow(zscores)
  nSnpsInitial <- sum(zscores$raissR2 == 2.0, na.rm = TRUE)
  nSnpsImputed <- sum(zscores$raissR2 != 2.0, na.rm = TRUE)
  nSnpsLdFilt <- sum(zscores$raissLdScore < minimumLd, na.rm = TRUE)
  nSnpsR2Filt <- sum(zscores$raissR2 < r2Threshold, na.rm = TRUE)

  # Apply filters
  zscoresNofilter <- zscores
  zscores <- zscores[zscores$raissR2 > r2Threshold & zscores$raissLdScore >= minimumLd, ]
  nSnpsAfFilt <- nrow(zscores)

  # Print report
  if (verbose) {
    maxLabelLength <- max(nchar(c(
      "Variants before filter:",
      "Non-imputed variants:",
      "Imputed variants:",
      "Variants filtered because of low LD score:",
      "Variants filtered because of low R2:",
      "Remaining variants after filter:"
    )))

    message("IMPUTATION REPORT\n")
    message(.formatRaissLine("Variants before filter", nSnpsBfFilt, maxLabelLength))
    message(.formatRaissLine("Non-imputed variants", nSnpsInitial, maxLabelLength))
    message(.formatRaissLine("Imputed variants", nSnpsImputed, maxLabelLength))
    message(.formatRaissLine("Variants filtered because of low LD score", nSnpsLdFilt, maxLabelLength))
    message(.formatRaissLine("Variants filtered because of low R2", nSnpsR2Filt, maxLabelLength))
    message(.formatRaissLine("Remaining variants after filter", nSnpsAfFilt, maxLabelLength))
  }
  return(zscore_list = list(zscoresNofilter = zscoresNofilter, zscores = zscores))
}

computeMu <- function(sigIT, sigTInv, zt) {
  return(sigIT %*% (sigTInv %*% zt))
}

computeVar <- function(sigIT, sigTInv, lamb, batch = TRUE) {
  if (batch) {
    var <- (1 + lamb) - rowSums((sigIT %*% sigTInv) * sigIT)
    raissLdScore <- rowSums(sigIT^2)
  } else {
    var <- (1 + lamb) - (sigIT %*% (sigTInv %*% t(sigIT)))
    raissLdScore <- sum(sigIT^2)
  }
  return(list(var = var, raissLdScore = raissLdScore))
}

checkInversion <- function(sigT, sigTInv) {
  return(all.equal(sigT, sigT %*% (sigTInv %*% sigT), tolerance = 1e-5))
}

varInBoundaries <- function(var, lamb) {
  var[var < 0] <- 0
  var[var > (0.99999 + lamb)] <- 1
  return(var)
}

invertMat <- function(mat, lamb, rcond) {
  tryCatch(
    {
      # Modify the diagonal elements of mat
      diag(mat) <- 1 + lamb
      # Compute the pseudo-inverse
      matInv <- ginv(mat, tol = rcond)
      return(matInv)
    },
    error = function(e) {
      # Second attempt with updated lamb and rcond in case of an error
      diag(mat) <- 1 + lamb * 1.1
      matInv <- ginv(mat, tol = rcond * 1.1)
      return(matInv)
    }
  )
}

invertMatRecursive <- function(mat, lamb, rcond) {
  tryCatch(
    {
      # Modify the diagonal elements of mat
      diag(mat) <- 1 + lamb
      # Compute the pseudo-inverse
      matInv <- ginv(mat, tol = rcond)
      return(matInv)
    },
    error = function(e) {
      # Recursive call with updated lamb and rcond in case of an error
      invertMat(mat, lamb * 1.1, rcond * 1.1)
    }
  )
}

invertMatEigen <- function(mat, tol = 1e-3) {
  eigenMat <- eigen(mat)
  L <- which(cumsum(eigenMat$values) / sum(eigenMat$values) > 1 - tol)[1]
  if (is.na(L)) {
    # all eigen values are extremely small
    stop("Cannot invert the input matrix because all its eigen values are negative or close to zero")
  }
  matInv <- eigenMat$vectors[, 1:L] %*%
    diag(1 / eigenMat$values[1:L]) %*%
    t(eigenMat$vectors[, 1:L])

  return(matInv)
}


# =============================================================================
# Top-level summaryStatsQc() pipeline + helpers
# =============================================================================



#' Detect LD-Summary Statistic Mismatches
#'
#' Unified wrapper for detecting outlier variants due to LD-summary statistic
#' mismatches. Dispatches to either \code{\link{dentistSingleWindow}} or
#' \code{\link{slalom}} based on the \code{method} argument.
#'
#' @param zScore Numeric vector of z-scores.
#' @param R Square LD correlation matrix. Provide either \code{R} or \code{X}.
#' @param X Genotype matrix (samples x SNPs). If provided, LD is computed via
#'   \code{\link{computeLd}} and \code{nSample} defaults to \code{nrow(X)}.
#' @param nSample Number of samples in the LD reference panel. Required when
#'   \code{R} is provided and \code{method = "dentist"}; inferred from \code{X}
#'   when \code{X} is provided.
#' @param method Character string specifying the QC method: \code{"slalom"}
#'   (default) or \code{"dentist"}.
#' @param ldMethod Character string specifying the LD computation method when
#'   \code{X} is provided. One of \code{"sample"} (default), \code{"population"},
#'   or \code{"gcta"}. Ignored when \code{R} is provided directly.
#' @param ... Additional arguments passed to the underlying QC method
#'   (\code{\link{dentistSingleWindow}} or \code{\link{slalom}}).
#'
#' @return A data frame with at least a logical \code{outlier} column indicating
#'   which variants are identified as outliers. The remaining columns depend on
#'   the method used.
#'
#' @seealso \code{\link{dentistSingleWindow}}, \code{\link{slalom}},
#'   \code{\link{summaryStatsQc}}
#' @importFrom dplyr mutate row_number filter pull
#' @export
ldMismatchQc <- function(zScore, R = NULL, X = NULL, nSample = NULL,
                         method = c("slalom", "dentist"),
                         ldMethod = "sample", ...) {
  method <- match.arg(method)
  if (method == "dentist") {
    qcResults <- dentistSingleWindow(zScore, R = R, X = X, nSample = nSample,
                                     ldMethod = ldMethod, ...)
    return(qcResults)
  } else {
    qcResults <- slalom(zScore, R = R, X = X, ldMethod = ldMethod, ...)
    # Standardize output: slalom uses "outliers", rename to "outlier" for consistency
    result <- qcResults$data
    if ("outliers" %in% colnames(result) && !"outlier" %in% colnames(result)) {
      colnames(result)[colnames(result) == "outliers"] <- "outlier"
    }
    return(result)
  }
}

.resolveZMismatchQc <- function(zMismatchQc) {
  if (is.null(zMismatchQc)) return("none")
  match.arg(zMismatchQc, c("none", "slalom", "dentist"))
}

#' Effective sample size for a case/control study
#'
#' Computes the effective sample size
#' \code{N_eff = 4 / (1/nCase + 1/nControl) = 4 * nCase * nControl /
#' (nCase + nControl)} for case/control GWAS. Balanced studies
#' (\code{nCase == nControl}) recover the total \code{nCase + nControl};
#' imbalanced studies give a smaller value, which is the statistically
#' correct sample size for the RSS likelihood, residual-variance estimation,
#' kriging, and the N-cutoff filter. Vectorized over \code{nCase} /
#' \code{nControl}; entries where either count is \code{NA} or \code{<= 0}
#' return \code{NA_real_}.
#'
#' @param nCase Numeric vector of case counts.
#' @param nControl Numeric vector of control counts.
#' @return Numeric vector of effective sample sizes (\code{NA_real_} where a
#'   count is missing or non-positive).
#' @references Privé et al., "Identifying and correcting for misspecifications
#'   in GWAS summary statistics and polygenic scores", HGG Advances 2022.
#' @export
effectiveN <- function(nCase, nControl) {
  nCase <- as.numeric(nCase)
  nControl <- as.numeric(nControl)
  out <- 4 / (1 / nCase + 1 / nControl)
  out[is.na(nCase) | is.na(nControl) | nCase <= 0 | nControl <= 0] <- NA_real_
  out
}

#' Kriging-style LD-consistency outlier QC
#'
#' Flags variants whose observed z-score is inconsistent with the value
#' predicted from its LD neighbours, using susieR's kriging diagnostic.
#' \code{susieR::kriging_rss()} computes the leave-one-out conditional
#' distribution of each \code{z_i} given the rest (with the LD-mismatch scale
#' \code{s} defaulting to \code{susieR::estimate_s_rss()}) and a per-variant
#' \code{logLR} for the allele-switch hypothesis. This helper reuses susieR's
#' own allele-switch rule --- \code{logLR > logLRThreshold & abs(z) > zThreshold}
#' (the same \code{logLR > 2 & |z| > 2} used in \code{susie_rss_utils}) --- to
#' flag variants whose sign should be flipped. RSS-only helper, opt-in via
#' \code{alleleFlipKriging}; never wired into \code{alleleQc()} /
#' \code{matchRefPanel()}. Requires a susieR that provides \code{kriging_rss()}
#' and \code{estimate_s_rss()}.
#'
#' @param zScore Numeric vector of harmonized z-scores.
#' @param R Square LD correlation matrix aligned to \code{zScore}.
#' @param n Sample size, forwarded to \code{susieR::kriging_rss()} (whose
#'   default \code{s} is \code{susieR::estimate_s_rss()}).
#' @param variantIds Optional variant IDs for the diagnostics table.
#' @param zThreshold Absolute-z cutoff for the allele-switch rule
#'   (default \code{2}, matching susieR).
#' @param logLRThreshold Log-likelihood-ratio cutoff for the allele-switch rule
#'   (default \code{2}, matching susieR).
#' @return A list with \code{flip} (logical vector; \code{TRUE} = allele switch,
#'   z-score should be sign-flipped) and \code{diagnostics} (data frame of
#'   per-variant \code{z}, \code{condmean}, \code{z_std_diff}, \code{logLR}, and
#'   the \code{flipped} flag).
#' @importFrom stats pnorm
#' @export
krigingOutlierQc <- function(zScore, R, n, variantIds = NULL,
                             zThreshold = 2, logLRThreshold = 2) {
  zScore <- as.numeric(zScore)
  m <- length(zScore)
  if (is.null(R) || !is.matrix(R) || nrow(R) != m || ncol(R) != m) {
    stop("krigingOutlierQc requires a square LD matrix aligned to zScore.")
  }
  if (missing(n) || length(n) != 1L || is.na(n) || !is.finite(n) || n <= 0) {
    stop("krigingOutlierQc requires a single positive sample size 'n'.")
  }
  if (!requireNamespace("susieR", quietly = TRUE) ||
      !all(c("estimate_s_rss", "kriging_rss") %in% getNamespaceExports("susieR"))) {
    # nocov start
    stop("krigingOutlierQc requires a susieR that provides estimate_s_rss() and ",
         "kriging_rss(); the installed susieR does not. Install a susieR with the ",
         "kriging RSS diagnostic, or disable alleleFlipKriging.")
    # nocov end
  }
  if (is.null(variantIds)) variantIds <- rownames(R)
  # susieR's kriging RSS diagnostic. Called exactly as the single-panel/notebook
  # call (`kriging_rss(z, R, n)`), so `s` defaults to estimate_s_rss() and the
  # returned `logLR` matches susieR's own allele-switch selection.
  cd <- susieR::kriging_rss(z = zScore, R = R, n = n)$conditional_dist
  condMean <- as.numeric(cd$condmean)
  zStdDiff <- as.numeric(cd$z_std_diff)
  logLR    <- as.numeric(cd$logLR)
  # susieR's allele-switch rule (susie_rss_utils.R): logLR > 2 & |z| > 2.
  flip <- !is.na(logLR) & !is.na(zScore) &
    logLR > logLRThreshold & abs(zScore) > zThreshold
  list(
    flip = flip,
    diagnostics = data.frame(
      variant_id = if (is.null(variantIds)) seq_len(m) else variantIds,
      z = zScore, condmean = condMean, z_std_diff = zStdDiff, logLR = logLR,
      flipped = flip,
      stringsAsFactors = FALSE
    )
  )
}

# =============================================================================
# summaryStatsQc — SumStats-input QC pipeline (replaces the previous
# data.frame/LdData/QcResult-based summaryStatsQc and rssBasicQc).
# =============================================================================

# Convert one entry's GRanges into a flat data.frame with the column shape
# harmonizeAlleles expects (lower-case chrom/pos plus the CapsCase mcols).
.entryGrangesToDf <- function(gr) {
  mc <- as.data.frame(S4Vectors::mcols(gr), stringsAsFactors = FALSE)
  out <- data.frame(
    chrom = sub("^chr", "", as.character(GenomicRanges::seqnames(gr)),
                ignore.case = TRUE),
    pos   = GenomicRanges::start(gr),
    stringsAsFactors = FALSE)
  cbind(out, mc)
}

# Build a refVariants data.frame (chrom, pos, A1, A2, variant_id) from the
# ldSketch GenotypeHandle's snpInfo so harmonizeAlleles can join by (chrom, pos).
.refVariantsFromSketch <- function(handle) {
  si <- getSnpInfo(handle)
  chr <- sub("^chr", "", as.character(si$CHR), ignore.case = TRUE)
  data.frame(
    chrom      = chr,
    pos        = as.integer(si$BP),
    A1         = as.character(si$A1),
    A2         = as.character(si$A2),
    variant_id = as.character(si$SNP),
    stringsAsFactors = FALSE)
}

# Reassemble a harmonized data.frame into a GRanges with the SumStats mcol
# shape (SNP, A1, A2, Z, N, ... optional MAF/INFO/BETA/SE/P kept if present).
.dfToEntryGranges <- function(df) {
  # Short-circuit on empty input: `paste0("chr", character(0))` returns
  # "chr" (a length-1 vector), not character(0), so we cannot rely on the
  # paste/IRanges constructors to handle the zero-row case cleanly.
  chrRaw <- as.character(df$chrom)
  if (length(chrRaw) == 0L) {
    gr <- GenomicRanges::GRanges()
    return(gr)
  }
  chr <- withChrPrefix(chrRaw)
  gr <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges   = IRanges::IRanges(start = as.integer(df$pos), width = 1L))
  if (!is.null(df$variant_id) && is.null(df$SNP)) df$SNP <- df$variant_id
  baseCols <- c("SNP", "A1", "A2", "Z", "N")
  optCols  <- c("MAF", "INFO", "BETA", "SE", "P", "N_CASE", "N_CONTROL")
  use <- intersect(c(baseCols, optCols), colnames(df))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(df[, use, drop = FALSE])
  gr
}

# -----------------------------------------------------------------------------
# Shared entry-to-sumstat data.frame converter
# -----------------------------------------------------------------------------

# Internal: convert one sumstat-entry GRanges into a flat data.frame with
# (variant_id, chrom, pos, A1, A2) base columns plus a configurable set
# of stats columns (z, beta, se, N, maf). Shared by the four pipelines
# that walk sumstats GRanges (fineMappingPipeline, twasWeights,
# ctwasPipeline, colocboostPipeline).
#
# `require`   character vector of mcol names that MUST be present; errors
#             when any is missing. Use this for the strict callers
#             (e.g. `.fmExtractZn` needs SNP + Z + N to proceed).
# `derive`    when "zFromBetaSe" and `z` is absent but BETA + SE are
#             present, set z := BETA/SE. Default "none".
# `label`     error-message prefix for missing-`require` errors.
# `keepChrPrefix`  when TRUE keep the seqname as-is ("chr1"); when FALSE,
#                  strip any leading "chr" so callers that expect numeric
#                  chrom (ctwas, colocboost) see "1".
.entryToSumstatDf <- function(gr,
                              require = character(0),
                              derive = c("none", "zFromBetaSe"),
                              label = "entry",
                              keepChrPrefix = TRUE) {
  derive <- match.arg(derive)
  mc <- S4Vectors::mcols(gr)
  for (col in require) {
    if (!(col %in% colnames(mc))) {
      stop(sprintf("%s: entry has no %s mcol.", label, col))
    }
  }
  chr <- as.character(GenomicRanges::seqnames(gr))
  if (!keepChrPrefix) chr <- sub("^chr", "", chr, ignore.case = TRUE)
  df <- data.frame(
    variant_id = if ("SNP" %in% colnames(mc)) as.character(mc$SNP)
                 else rep(NA_character_, length(gr)),
    chrom      = chr,
    pos        = as.integer(GenomicRanges::start(gr)),
    A1         = if ("A1" %in% colnames(mc)) as.character(mc$A1)
                 else rep(NA_character_, length(gr)),
    A2         = if ("A2" %in% colnames(mc)) as.character(mc$A2)
                 else rep(NA_character_, length(gr)),
    stringsAsFactors = FALSE)
  if ("Z"    %in% colnames(mc)) df$z    <- as.numeric(mc$Z)
  if ("BETA" %in% colnames(mc)) df$beta <- as.numeric(mc$BETA)
  if ("SE"   %in% colnames(mc)) df$se   <- as.numeric(mc$SE)
  if ("N"    %in% colnames(mc)) df$N    <- as.numeric(mc$N)
  if ("MAF"  %in% colnames(mc)) df$maf  <- as.numeric(mc$MAF)
  if (derive == "zFromBetaSe" && is.null(df$z) &&
      !is.null(df$beta) && !is.null(df$se)) {
    df$z <- df$beta / df$se
  }
  df
}

# Derive BETA and SE columns from signed Z when the entry has only Z.
# Formula (Zhu et al. 2016 / RAISS):
#   se   = 1 / sqrt(2 * maf * (1 - maf) * (N + z^2))
#   beta = z * se
# Requires Z, MAF, and N to all be present in `df`. No-op if BETA and SE
# are already there, or if any required column is missing. Returns:
#   list(df = <data.frame>, audit = NULL | list(nDerived = <int>))
# .zToPvalue (two-tailed normal p-value from a signed Z) is defined once in
# pvalCombine.R and shared package-wide.

# Internal: thin SVD with numerical-stability filtering. Drops singular
# values below `tol * max(d)` and caps the retained rank at `maxRank`.
# Used by RAISS imputation (raissSingleMatrixFromX) to invert a panel
# genotype matrix safely under rank deficiency.
.safeSvd <- function(mat, tol = 1e-8, maxRank = NULL) {
  if (max(abs(mat)) == 0)
    stop("Cannot compute SVD of an all-zero matrix.")
  s <- svd(mat)
  d <- s$d
  keep <- if (tol > 0 && length(d) > 0) {
    out <- d / d[1] > tol
    if (!any(out))
      stop("All singular values are below the tolerance threshold.")
    out
  } else rep(TRUE, length(d))
  if (!is.null(maxRank) && maxRank > 0) {
    nKeep <- min(sum(keep), maxRank)
    keepIdx <- which(keep)
    if (length(keepIdx) > nKeep)
      keep[keepIdx[(nKeep + 1):length(keepIdx)]] <- FALSE
  }
  list(u = s$u[, keep, drop = FALSE],
       d = d[keep],
       v = s$v[, keep, drop = FALSE])
}

# Internal: identify LD-correlated duplicate variants. Walks the LD
# matrix left-to-right, marking each variant as either a unique anchor
# (dupBearer == -1) or a duplicate of an earlier anchor (dupBearer == k,
# the anchor's index). Returns filtered z / LD plus per-variant
# bookkeeping that DENTIST's addDupsBackDentist uses to splice the
# dropped variants back into the output.
.findDuplicateVariants <- function(z, ld, rThreshold) {
  p <- length(z)
  dupBearer <- rep(-1, p)
  corABS <- rep(0, p)
  sign <- rep(1, p)
  count <- 1L
  minValue <- 1
  for (i in seq_len(p - 1L)) {
    if (dupBearer[i] != -1) next
    idx <- (i + 1L):p
    corVec <- abs(ld[i, idx])
    dupIdx <- which(dupBearer[idx] == -1 & corVec > rThreshold)
    if (length(dupIdx) > 0) {
      j <- idx[dupIdx]
      sign[j] <- ifelse(ld[i, j] < 0, -1, sign[j])
      corABS[j] <- corVec[dupIdx]
      dupBearer[j] <- count
    }
    minValue <- min(minValue, min(corVec))
    count <- count + 1L
  }
  filteredZ <- z[dupBearer == -1]
  filteredLD <- ld[dupBearer == -1, dupBearer == -1, drop = FALSE]
  list(filteredZ = filteredZ, filteredLD = filteredLD,
       dupBearer = dupBearer, corABS = corABS, sign = sign,
       minValue = minValue)
}

.deriveBetaSeFromZ <- function(df) {
  hasBeta <- "BETA" %in% colnames(df)
  hasSe   <- "SE"   %in% colnames(df)
  if (hasBeta && hasSe) return(list(df = df, audit = NULL))
  hasZ   <- "Z"   %in% colnames(df)
  hasMaf <- "MAF" %in% colnames(df)
  hasN   <- "N"   %in% colnames(df)
  if (!(hasZ && hasMaf && hasN)) return(list(df = df, audit = NULL))
  z   <- as.numeric(df$Z)
  maf <- as.numeric(df$MAF)
  n   <- as.numeric(df$N)
  bs   <- .zToBetaSe(z, maf, n)
  se   <- bs$se
  beta <- bs$beta
  if (!hasBeta) df$BETA <- beta
  if (!hasSe)   df$SE   <- se
  list(df = df, audit = list(nDerived = sum(!is.na(se))))
}

# Drop variants whose (chrom, pos) overlaps any user-supplied skipRegion.
# skipRegion may be a character vector of "chr:start-end" strings or a GRanges.
.applySkipRegion <- function(df, skipRegion) {
  if (is.null(skipRegion) || length(skipRegion) == 0L) return(df)
  if (is.character(skipRegion)) {
    parsed <- do.call(rbind, lapply(skipRegion, function(s) {
      m <- regmatches(s, regexec("^([^:]+):([0-9]+)-([0-9]+)$", s))[[1]]
      if (length(m) != 4L)
        stop("skipRegion entry must be 'chr:start-end'; got '", s, "'")
      data.frame(chrom = sub("^chr", "", m[2], ignore.case = TRUE),
                 start = as.integer(m[3]),
                 end   = as.integer(m[4]),
                 stringsAsFactors = FALSE)
    }))
  } else if (methods::is(skipRegion, "GRanges")) {
    parsed <- data.frame(
      chrom = sub("^chr", "", as.character(GenomicRanges::seqnames(skipRegion)),
                  ignore.case = TRUE),
      start = GenomicRanges::start(skipRegion),
      end   = GenomicRanges::end(skipRegion),
      stringsAsFactors = FALSE)
  } else {
    stop("skipRegion must be a character vector of 'chr:start-end' ",
         "strings or a GRanges.")
  }
  dropMask <- rep(FALSE, nrow(df))
  dfChr <- sub("^chr", "", as.character(df$chrom), ignore.case = TRUE)
  for (i in seq_len(nrow(parsed))) {
    dropMask <- dropMask |
      (dfChr == parsed$chrom[i] &
       df$pos >= parsed$start[i] &
       df$pos <= parsed$end[i])
  }
  df[!dropMask, , drop = FALSE]
}

# Apply the panel-vs-sumstats allele harmonization using the slim
# harmonizeAlleles against the ldSketch's variant info. Threads the
# variant-level filters (indels, strand-ambiguous, duplicates) through
# so the LD-panel-anchored pass handles them in a single sweep.
.matchAgainstSketch <- function(df, ldSketch, matchMinProp,
                                removeIndels = FALSE,
                                removeStrandAmbiguous = TRUE,
                                removeDups = TRUE) {
  refVariants <- .refVariantsFromSketch(ldSketch)
  flipCandidates <- c("Z", "BETA")
  colToFlip <- intersect(flipCandidates, colnames(df))
  if (length(colToFlip) == 0L)
    stop("summaryStatsQc: input entry must contain at least one of Z or BETA ",
         "before panel harmonization.")
  colToComplement <- intersect("MAF", colnames(df))
  if (!"A1" %in% colnames(df) || !"A2" %in% colnames(df))
    stop("summaryStatsQc: input entry must contain A1 and A2 columns.")
  res <- harmonizeAlleles(
    targetData            = df,
    refVariants           = refVariants,
    colToFlip             = colToFlip,
    colToComplement       = colToComplement,
    matchMinProp          = matchMinProp,
    removeUnmatched       = TRUE,
    removeIndels          = removeIndels,
    removeStrandAmbiguous = removeStrandAmbiguous,
    removeDups            = removeDups)
  out <- res$harmonizedData
  if (!"chrom" %in% colnames(out) && "chr" %in% colnames(out))
    colnames(out)[colnames(out) == "chr"] <- "chrom"
  attr(out, "qcCounts") <- attr(res, "qcCounts")
  out
}

# Variant-content filters (MAF / INFO / N). Pure data-frame column
# filters; no Bioconductor genome packages needed.
#
# mafCutoff:  drop rows where MAF (or FRQ) < mafCutoff. Requires either
#             column when mafCutoff > 0; errors if neither is present.
# infoCutoff: drop rows where INFO < infoCutoff. Requires INFO column
#             when infoCutoff > 0.
# nCutoff:    drop rows whose N is more than nCutoff median-absolute-
#             deviations from the median (a 5-MAD-from-median cap on
#             per-variant N). Set nCutoff = 0 to disable. Rows with NA N
#             are always dropped.
.applyContentFilters <- function(df, mafCutoff = 0, infoCutoff = 0,
                                 nCutoff = 5) {
  audit <- list()
  if (mafCutoff > 0) {
    mafCol <- intersect(c("MAF", "FRQ"), colnames(df))[1L]
    if (is.na(mafCol))
      stop(".applyContentFilters: mafCutoff > 0 requires a MAF or FRQ column.")
    before <- nrow(df)
    mafVals <- as.numeric(df[[mafCol]])
    # Normalise effect-allele frequency to MAF: take min(af, 1-af).
    mafVals <- pmin(mafVals, 1 - mafVals, na.rm = FALSE)
    df <- df[!is.na(mafVals) & mafVals >= mafCutoff, , drop = FALSE]
    audit$mafDropped <- before - nrow(df)
  }
  if (infoCutoff > 0) {
    if (!"INFO" %in% colnames(df))
      stop(".applyContentFilters: infoCutoff > 0 requires an INFO column.")
    before <- nrow(df)
    infoVals <- as.numeric(df$INFO)
    df <- df[!is.na(infoVals) & infoVals >= infoCutoff, , drop = FALSE]
    audit$infoDropped <- before - nrow(df)
  }
  if (nCutoff > 0 && "N" %in% colnames(df) && nrow(df) > 0L) {
    nVals <- as.numeric(df$N)
    before <- nrow(df)
    if (any(is.na(nVals))) {
      df <- df[!is.na(nVals), , drop = FALSE]
      nVals <- nVals[!is.na(nVals)]
    }
    if (length(nVals) > 0L) {
      medN <- stats::median(nVals)
      madN <- stats::mad(nVals, constant = 1)
      if (madN > 0) {
        zN <- abs(nVals - medN) / madN
        df <- df[zN <= nCutoff, , drop = FALSE]
      }
    }
    audit$nDropped <- before - nrow(df)
  }
  list(df = df, audit = audit)
}

# Per-row variant sanity / hygiene checks ported from MungeSumstats's
# check_*.R series but rewritten as pure data.frame operations with no
# genome / dbSNP dependency. Each step is gated by its own flag so a
# caller can disable any single check.
#
# Steps (in order; each contributes a count to audit):
#   - coerceNumeric: cast signed columns to numeric (catches stray "0.5"
#       strings). NA-introducing coercions are counted.
#   - normalizeChr:   strip "chr"/"ch" prefix, uppercase X/Y/MT, map
#       23->X, 24->Y, M->MT. Optional dropNonstandardChr removes rows
#       whose CHR is outside 1..22, X, Y, MT after normalization.
#   - dropMissData:   drop rows with NA in any vital column (chrom, pos,
#       A1, A2, and at least one of Z / BETA).
#   - dropPOutOfRange: drop rows where P < 0 or P > 1 (corrupt p-values).
#       Only fires when a P column is present.
#   - clampSmallP:    floor 0 <= P <= smallPFloor to smallPFloor so
#       -log10(P) stays finite downstream.
#   - dropZeroEffect: drop rows where any effect column is exactly 0
#       (BETA / LOG_ODDS / SIGNED_SUMSTAT) or OR is exactly 1. MungeSumstats
#       treats these as degenerate / artefactual.
#   - dropNonpositiveSe: drop rows where SE <= 0.
.applySanityChecks <- function(df,
                                coerceNumeric        = TRUE,
                                normalizeChr         = TRUE,
                                dropNonstandardChr   = TRUE,
                                dropMissData         = TRUE,
                                dropPOutOfRange      = TRUE,
                                clampSmallP          = TRUE,
                                smallPFloor          = 5e-324,
                                dropZeroEffect       = TRUE,
                                dropNonpositiveSe    = TRUE) {
  audit <- list()
  if (nrow(df) == 0L) return(list(df = df, audit = audit))

  if (coerceNumeric) {
    numericCols <- intersect(
      c("Z", "BETA", "SE", "OR", "LOG_ODDS", "SIGNED_SUMSTAT",
        "P", "MAF", "FRQ", "INFO", "N"),
      colnames(df))
    naIntroduced <- 0L
    for (col in numericCols) {
      orig <- df[[col]]
      if (is.numeric(orig)) next
      coerced <- suppressWarnings(as.numeric(orig))
      naIntroduced <- naIntroduced +
        sum(is.na(coerced) & !is.na(orig))
      df[[col]] <- coerced
    }
    if (naIntroduced > 0L) audit$nonNumericCoerced <- naIntroduced
  }

  if (normalizeChr && "chrom" %in% colnames(df)) {
    chr <- as.character(df$chrom)
    chr <- sub("^chr", "", chr, ignore.case = TRUE)
    chr <- sub("^ch",  "", chr, ignore.case = TRUE)
    chr <- toupper(chr)
    chr[chr == "23"] <- "X"
    chr[chr == "24"] <- "Y"
    chr[chr == "M"]  <- "MT"
    df$chrom <- chr
    if (dropNonstandardChr) {
      before <- nrow(df)
      standardChrs <- c(as.character(1:22), "X", "Y", "MT")
      df <- df[chr %in% standardChrs, , drop = FALSE]
      dropped <- before - nrow(df)
      if (dropped > 0L) audit$nonstandardChrDropped <- dropped
    }
  }

  if (dropMissData && nrow(df) > 0L) {
    vital <- intersect(c("chrom", "pos", "A1", "A2"), colnames(df))
    signedCol <- intersect(c("Z", "BETA"), colnames(df))[1L]
    if (!is.na(signedCol)) vital <- c(vital, signedCol)
    if (length(vital) > 0L) {
      before <- nrow(df)
      bad <- Reduce(`|`, lapply(vital, function(c) is.na(df[[c]])))
      if (any(bad)) df <- df[!bad, , drop = FALSE]
      dropped <- before - nrow(df)
      if (dropped > 0L) audit$missDataDropped <- dropped
    }
  }

  if (dropPOutOfRange && "P" %in% colnames(df) && nrow(df) > 0L) {
    before <- nrow(df)
    p <- as.numeric(df$P)
    bad <- !is.na(p) & (p < 0 | p > 1)
    if (any(bad)) df <- df[!bad, , drop = FALSE]
    dropped <- before - nrow(df)
    if (dropped > 0L) audit$pOutOfRangeDropped <- dropped
  }

  if (clampSmallP && "P" %in% colnames(df) && nrow(df) > 0L) {
    p <- as.numeric(df$P)
    smallMask <- !is.na(p) & p >= 0 & p < smallPFloor
    nClamped <- sum(smallMask)
    if (nClamped > 0L) {
      df$P[smallMask] <- smallPFloor
      audit$smallPClamped <- nClamped
    }
  }

  if (dropZeroEffect && nrow(df) > 0L) {
    effectCols <- intersect(
      c("BETA", "LOG_ODDS", "SIGNED_SUMSTAT", "OR"), colnames(df))
    if (length(effectCols) > 0L) {
      before <- nrow(df)
      badMask <- rep(FALSE, nrow(df))
      for (col in effectCols) {
        vals <- as.numeric(df[[col]])
        sentinel <- if (col == "OR") 1 else 0
        badMask <- badMask | (!is.na(vals) & vals == sentinel)
      }
      if (any(badMask)) df <- df[!badMask, , drop = FALSE]
      dropped <- before - nrow(df)
      if (dropped > 0L) audit$zeroEffectDropped <- dropped
    }
  }

  if (dropNonpositiveSe && "SE" %in% colnames(df) && nrow(df) > 0L) {
    before <- nrow(df)
    se <- as.numeric(df$SE)
    bad <- !is.na(se) & se <= 0
    if (any(bad)) df <- df[!bad, , drop = FALSE]
    dropped <- before - nrow(df)
    if (dropped > 0L) audit$nonpositiveSeDropped <- dropped
  }

  list(df = df, audit = audit)
}

# Apply ldMismatchQc (SLALOM/DENTIST) against the LD sketch. Returns the
# filtered df, outlier count, and the full per-variant diagnostics table
# (the data.frame returned by ldMismatchQc(), prepended with a
# variant_id column for downstream joins). Callers stamp `diagnostics`
# into the entry's qcInfo audit so the per-variant detail is available
# for plotting / postprocessing instead of just the outlier count.
.applyLdMismatchQcToEntry <- function(df, ldSketch, method) {
  variantIds <- df$SNP
  if (is.null(variantIds) || any(is.na(variantIds)))
    stop("summaryStatsQc: ldMismatchQc requires SNP column on the entry.")
  # Panel LD for the entry variants via the shared LD-from-sketch helper (tuple
  # match with chr-prefix tolerance, strand-ambiguous variants kept; errors if
  # any variant is absent from the panel).
  R <- .ldFromSketch(ldSketch, variantIds, label = "summaryStatsQc: zMismatchQc")
  qc <- ldMismatchQc(zScore = df$Z, R = R, nSample = getNSamples(ldSketch),
                     method = method)
  # slalom / dentist can leave NA in the outlier column when their
  # per-variant statistic is undefined (e.g. a degenerate dentist
  # chisq for variants effectively orthogonal to the lead). Treat NA as
  # "no evidence of being an outlier" (conservative: keep the variant)
  # so the downstream df / sum() / IRanges construction stay finite.
  outlierFlags <- qc$outlier
  outlierFlags[is.na(outlierFlags)] <- FALSE
  # Attach the variant_id column so the diagnostics data.frame stays
  # self-describing once it's separated from the input df.
  diagnostics <- if (is.data.frame(qc)) {
    cbind(variant_id = as.character(variantIds), qc,
          stringsAsFactors = FALSE)
  } else NULL
  list(df = df[!outlierFlags, , drop = FALSE],
       outliers = sum(outlierFlags),
       diagnostics = diagnostics)
}

# -----------------------------------------------------------------------------
# Signal screen: skip an entry/block with no strong signal by a chosen metric.
# -----------------------------------------------------------------------------
# The screen is driven by ONE metric at a time (enforced by .resolveScreenMetric):
#   pip   : max single-effect PIP                (susie_ser $pip); cutoff<0 => 3/nVar
#   absZ  : max |Z|                              (no model fit)
#   logBf : max per-variant single-effect logBF  (susie_ser $lbf_variable)
#   bf    : same evidence on the raw BF scale    (compare maxlogBF > log(cutoff))
# A "screen spec" flowing through the pipelines is EITHER a legacy PIP cutoff
# (numeric scalar, 0 = off -- the historical `pipCutoffToSkip`) OR a resolved
# screen object `list(metric, cutoff)`. .asScreen() canonicalizes either into
# `list(metric, cutoff)` or NULL (no screen); the pipeline channels stay
# untyped so only the screen primitives need to interpret the spec.

# Canonicalize a screen spec into list(metric, cutoff) or NULL (no screen).
.asScreen <- function(spec) {
  if (is.null(spec)) return(NULL)
  if (is.list(spec)) {                                   # already a screen object
    if (is.null(spec$metric) || is.null(spec$cutoff) ||
        is.na(spec$cutoff) || spec$cutoff == 0) return(NULL)
    return(spec)
  }
  # Legacy numeric: a PIP cutoff. Only a non-zero scalar activates it (a
  # non-scalar / NA / 0 means no screen), matching the historical behaviour.
  if (length(spec) != 1L || is.na(spec) || spec == 0) return(NULL)
  list(metric = "pip", cutoff = as.numeric(spec))
}

# Turn the four public cutoff arguments into a single screen object (or NULL).
# Enforces one-metric-at-a-time and rejects meaningless negative cutoffs. The
# pip metric keeps its `< 0 => 3 / nVariants` adaptive convention (resolved at
# screen time); absZ / bf must be > 0; logBf may be any non-zero value.
.resolveScreenMetric <- function(pipCutoffToSkip = 0, absZCutoffToSkip = 0,
                                 bfCutoffToSkip = 0, logBfCutoffToSkip = 0) {
  scalar <- function(x) if (is.null(x) || length(x) != 1L || is.na(x)) 0 else as.numeric(x)
  cuts <- c(pip = scalar(pipCutoffToSkip), absZ = scalar(absZCutoffToSkip),
            bf = scalar(bfCutoffToSkip), logBf = scalar(logBfCutoffToSkip))
  on <- cuts[cuts != 0]
  if (length(on) == 0L) return(NULL)
  if (length(on) > 1L)
    stop("Only one signal screen may be enabled at a time, but these are ",
         "non-zero: ", paste(sprintf("%s=%g", names(on), on), collapse = ", "),
         ". Set all but one of pipCutoffToSkip / absZCutoffToSkip / ",
         "bfCutoffToSkip / logBfCutoffToSkip to 0.")
  metric <- names(on); cutoff <- unname(on[[1L]])
  if (metric == "absZ" && cutoff < 0)
    stop("absZCutoffToSkip must be > 0 (it screens on max|Z|).")
  if (metric == "bf" && cutoff < 0)
    stop("bfCutoffToSkip must be > 0 (Bayes factors are positive).")
  list(metric = metric, cutoff = cutoff)
}

# Decide whether the z-scores of one entry clear the chosen screen. Returns
# list(ok = logical, reason = character). susie_ser is fit at most once, and
# only for the metrics that need it (absZ stays model-free).
.entryScreenPass <- function(z, n, nVar, scr) {
  metric <- scr$metric; cutoff <- scr$cutoff
  if (metric == "absZ") {
    m <- suppressWarnings(max(abs(as.numeric(z)), na.rm = TRUE))
    return(list(ok = is.finite(m) && m > cutoff,
                reason = sprintf("no variant with |Z| above %g (max |Z| = %g)",
                                 cutoff, m)))
  }
  ser <- susieR::susie_ser(z = z, n = n, coverage = NULL)
  if (metric == "pip") {
    eff <- if (cutoff < 0) 3 / nVar else cutoff
    return(list(ok = any(ser$pip > eff),
                reason = sprintf("no signals above PIP threshold %g", eff)))
  }
  maxLbf <- suppressWarnings(max(as.numeric(ser$lbf_variable), na.rm = TRUE))
  if (metric == "logBf")
    return(list(ok = is.finite(maxLbf) && maxLbf > cutoff,
                reason = sprintf("no variant with logBF above %g (max logBF = %g)",
                                 cutoff, maxLbf)))
  # metric == "bf": compare in log space to avoid overflow of exp(maxLbf).
  list(ok = is.finite(maxLbf) && maxLbf > log(cutoff),
       reason = sprintf("no variant with BF above %g (max BF = %g)",
                        cutoff, exp(maxLbf)))
}

# Per-entry signal screen. `screen` is a screen spec (see .asScreen). Skips
# (empties) the entry when the chosen metric shows no signal above its cutoff.
.applyEntryScreen <- function(df, n, screen) {
  scr <- .asScreen(screen)
  if (is.null(scr)) return(list(df = df, skipped = FALSE))
  res <- .entryScreenPass(df$Z, n = n, nVar = nrow(df), scr = scr)
  if (!res$ok) {
    return(list(df = df[FALSE, , drop = FALSE], skipped = TRUE,
                reason = res$reason))
  }
  list(df = df, skipped = FALSE)
}

# Back-compat thin wrapper for the original PIP-only screen (a bare numeric
# cutoff, 0 = off). Delegates to the generalized screen via .asScreen.
.applyPipScreen <- function(df, n, cutoff) {
  .applyEntryScreen(df, n = n, screen = cutoff)
}

# Prefix QC-track log lines with the entry label `lbl` (as `[lbl] ...`), or emit
# them bare when `lbl` is NA.
# @noRd
.qcEmit <- function(lbl, ...) {
  if (is.na(lbl)) message(...) else message("[", lbl, "] ", ...)
}

# Internal: canonicalize the working per-variant `N`. The N source is resolved
# by a four-level priority: (1) per-variant N_CASE / N_CONTROL columns, (2) study
# -level nCase / nControl scalars, (3) a per-variant N column, (4) a study-level
# nSample scalar (total N). Levels 1-2 give the effective sample size (default)
# or the raw total (escape hatch). Returns list(df=, nSource=), where nSource is
# "effective" | "column" | "total" | "study-n" | NA_character_ (no source).
# The entry label `lbl` prefixes the counts-win override log. A study-level
# scalar fills `df$N` even when the entry has no per-variant N column.
.resolveEffectiveN <- function(df, opts, lbl) {
  hasCols <- all(c("N_CASE", "N_CONTROL") %in% colnames(df))
  hasScalar <- !is.null(opts$nCase) && !is.null(opts$nControl) &&
    length(opts$nCase) == 1L && length(opts$nControl) == 1L &&
    is.finite(opts$nCase) && is.finite(opts$nControl) &&
    opts$nCase > 0 && opts$nControl > 0
  hasNSample <- !is.null(opts$nSample) && length(opts$nSample) == 1L &&
    is.finite(opts$nSample) && opts$nSample > 0
  hasN <- "N" %in% colnames(df)
  nRow <- nrow(df)

  if (!isTRUE(opts$effectiveN)) {
    # Escape hatch: prefer the raw N column; else raw total from counts; else the
    # study-level total N. No effective-N conversion, no override.
    if (hasN)       return(list(df = df, nSource = "column"))
    if (hasCols) {
      df$N <- as.numeric(df$N_CASE) + as.numeric(df$N_CONTROL)
      return(list(df = df, nSource = "total"))
    }
    if (hasScalar) {
      df$N <- rep(opts$nCase + opts$nControl, nRow)
      return(list(df = df, nSource = "total"))
    }
    if (hasNSample) {
      df$N <- rep(opts$nSample, nRow)
      return(list(df = df, nSource = "study-n"))
    }
    return(list(df = df, nSource = NA_character_))
  }

  # Default: per-variant c/c -> study c/c -> per-variant N -> study nSample.
  if (hasCols) {
    if (hasN) .qcEmit(lbl, "QC track: N overridden by effective N from per-variant ",
                   "n_case/n_control.")
    df$N <- effectiveN(df$N_CASE, df$N_CONTROL)
    return(list(df = df, nSource = "effective"))
  }
  if (hasScalar) {
    if (hasN) .qcEmit(lbl, "QC track: N overridden by effective N from study ",
                   "nCase/nControl.")
    df$N <- rep(effectiveN(opts$nCase, opts$nControl), nRow)
    return(list(df = df, nSource = "effective"))
  }
  if (hasN) return(list(df = df, nSource = "column"))
  if (hasNSample) {
    df$N <- rep(opts$nSample, nRow)
    return(list(df = df, nSource = "study-n"))
  }
  list(df = df, nSource = NA_character_)
}

# Internal: per-entry pipeline. Returns the cleaned GRanges and an audit list.
.runEntrySummaryStatsQc <- function(gr, ldSketch, refGenome, opts,
                                    entryLabel = NULL) {
  entryAudit <- list()
  # Counter capture for per-step "kept N of M" messages + per-entry rollup.
  # Each step records what it removed / added so we can summarise without
  # re-scanning the data. Skipped steps are left as NA / 0 and omitted
  # from the rollup.
  qcCount <- list(
    harmCorrSign = 0L, harmCorrStrand = 0L, harmDropped = 0L,
    krigingFlipped = 0L, mismatchRemoved = 0L,
    imputeBefore = NA_integer_, imputeAfter = NA_integer_)
  lbl <- if (!is.null(entryLabel) && nzchar(entryLabel)) entryLabel
         else NA_character_
  df <- .entryGrangesToDf(gr)
  entryAudit$variantsIn <- nrow(df)
  nStudyIn <- nrow(df)

  # 1. Per-row sanity checks (drop bad P / zero effect / non-positive SE,
  # clamp tiny P, coerce numeric, normalize CHR, drop missing-data rows).
  # Runs before any other filtering so downstream steps see clean values.
  nSanIn <- nrow(df)
  sanity <- .applySanityChecks(
    df,
    coerceNumeric      = opts$coerceNumeric,
    normalizeChr       = opts$normalizeChr,
    dropNonstandardChr = opts$dropNonstandardChr,
    dropMissData       = opts$dropMissData,
    dropPOutOfRange    = opts$dropPOutOfRange,
    clampSmallP        = opts$clampSmallP,
    smallPFloor        = opts$smallPFloor,
    dropZeroEffect     = opts$dropZeroEffect,
    dropNonpositiveSe  = opts$dropNonpositiveSe)
  df <- sanity$df
  if (length(sanity$audit) > 0L) entryAudit$sanityChecks <- sanity$audit
  if (nSanIn > 0L && nrow(df) != nSanIn) {
    .qcEmit(lbl, "QC track: sanity checks kept ", nrow(df), " of ", nSanIn,
         " variant(s).")
  }

  # Canonicalize N to the effective sample size for case/control input,
  # BEFORE the N-cutoff filter so the filter, kriging (median df$N), and the
  # rebuilt entry all consume N_eff. No-op for quantitative traits (no counts).
  nRes <- .resolveEffectiveN(df, opts, lbl)
  df <- nRes$df
  entryAudit$nSource <- nRes$nSource
  # Keep the PIP-screen / kriging sample size consistent with the canonicalized
  # N: when effective (or total) N was applied, nForPip (set upstream from the
  # entry's raw N) must follow the same N the N-filter/fit use.
  if (isTRUE(nRes$nSource %in% c("effective", "total")) && "N" %in% colnames(df))
    opts$nForPip <- stats::median(as.numeric(df$N), na.rm = TRUE)

  # 2. Variant-content filters (MAF / INFO / N). Pure column-numeric
  # filters; the indel / strand-ambiguous variant-allele filtering happens
  # inside .matchAgainstSketch via harmonizeAlleles against the LD panel.
  nFiltIn <- nrow(df)
  contentFiltered <- .applyContentFilters(
    df,
    mafCutoff  = opts$mafCutoff,
    infoCutoff = opts$infoCutoff,
    nCutoff    = opts$nCutoff)
  df <- contentFiltered$df
  if (length(contentFiltered$audit) > 0L)
    entryAudit$contentFilters <- contentFiltered$audit
  if (nFiltIn > 0L && nrow(df) != nFiltIn) {
    .qcEmit(lbl, "QC track: MAF/INFO/N filters kept ", nrow(df), " of ", nFiltIn,
         " variant(s).")
  }

  # 3. Derive BETA / SE from signed Z when the input only carries Z.
  # Formula (Zhu et al. 2016 / RAISS): se = 1/sqrt(2 * maf * (1-maf) *
  # (N + z^2)); beta = z * se. Requires Z, MAF, and N to all be present.
  derived <- .deriveBetaSeFromZ(df)
  df <- derived$df
  if (!is.null(derived$audit)) entryAudit$betaSeFromZ <- derived$audit

  # 4. Derive P-values from Z when the entry carries Z but not P.
  # Standard two-tailed normal: p = 2 * pnorm(-|z|). Very large |Z| can
  # underflow to P = 0, so re-apply the small-P clamp afterwards.
  if ("Z" %in% colnames(df) && !"P" %in% colnames(df)) {
    df$P <- .zToPvalue(df$Z)
    entryAudit$pValueFromZ <- sum(!is.na(df$P))
    if (isTRUE(opts$clampSmallP) && nrow(df) > 0L) {
      smallMask <- !is.na(df$P) & df$P >= 0 & df$P < opts$smallPFloor
      nClamped <- sum(smallMask)
      if (nClamped > 0L) {
        df$P[smallMask] <- opts$smallPFloor
        prev <- entryAudit$sanityChecks$smallPClamped %||% 0L
        if (is.null(entryAudit$sanityChecks))
          entryAudit$sanityChecks <- list()
        entryAudit$sanityChecks$smallPClamped <- prev + nClamped
      }
    }
  }

  # 2. keepVariants subset.
  if (length(opts$keepVariants) > 0L) {
    before <- nrow(df)
    df <- df[df$SNP %in% opts$keepVariants, , drop = FALSE]
    entryAudit$keepVariantsDropped <- before - nrow(df)
  }

  # 3. skipRegion drop.
  if (!is.null(opts$skipRegion) && length(opts$skipRegion) > 0L) {
    before <- nrow(df)
    df <- .applySkipRegion(df, opts$skipRegion)
    entryAudit$skipRegionDropped <- before - nrow(df)
  }

  if (nrow(df) < 2L) {
    entryAudit$earlyExit <- "fewer than two variants after pre-harmonization QC"
    return(list(gr = .dfToEntryGranges(df), audit = entryAudit))
  }

  # 4. Panel-vs-sumstats allele harmonization.
  nHarmIn <- nrow(df)
  df <- .matchAgainstSketch(
    df, ldSketch,
    matchMinProp          = opts$matchMinProp,
    removeIndels          = opts$removeIndels,
    removeStrandAmbiguous = opts$removeStrandAmbiguous,
    removeDups            = TRUE)
  harmCounts <- attr(df, "qcCounts")
  attr(df, "qcCounts") <- NULL
  # Re-key SNP to the harmonized id. .matchAgainstSketch rewrites variant_id to
  # the panel's allele orientation (chr:pos:A2:A1) and sign-flips Z/BETA for
  # swapped variants, but leaves the original SNP untouched. Every subsequent
  # panel lookup keys on SNP (kriging, z-mismatch QC) as does the downstream
  # fine-mapping LD lookup, so a stale SNP makes sign-flipped variants miss the
  # panel. Sync it here so harmonized alleles and id stay consistent.
  if (!is.null(df$variant_id)) df$SNP <- df$variant_id
  entryAudit$matchedAgainstSketch <- nrow(df)
  if (!is.null(harmCounts)) {
    qcCount$harmCorrSign   <- harmCounts$signFlip
    qcCount$harmCorrStrand <- harmCounts$strandFlip
    qcCount$harmDropped    <- nHarmIn - nrow(df)
    .qcEmit(lbl, "QC track: harmonization kept ", nrow(df), " of ", nHarmIn,
         " variant(s) (corrected: sign-flipped ", harmCounts$signFlip,
         ", strand-flipped ", harmCounts$strandFlip,
         "; dropped ", qcCount$harmDropped, ").")
  } else {
    qcCount$harmDropped <- nHarmIn - nrow(df)
    .qcEmit(lbl, "QC track: harmonization kept ", nrow(df), " of ", nHarmIn,
         " variant(s).")
  }

  # 5. Optional signal screen (after harmonization, on panel-aligned variants).
  #    One metric at a time: PIP / max|Z| / max BF / max logBF (opts$screen).
  if (!is.null(opts$screen)) {
    scr <- .applyEntryScreen(df, n = opts$nForPip, screen = opts$screen)
    df <- scr$df
    entryAudit$pipScreenSkipped <- isTRUE(scr$skipped)
    if (isTRUE(scr$skipped)) entryAudit$pipScreenReason <- scr$reason
  }

  # 6. Optional kriging allele-flip QC. Uses susieR's allele-switch rule
  #    (logLR > 2 & |z| > 2) to sign-flip inconsistent z-scores in place --
  #    variants are corrected and RETAINED, not dropped.
  if (isTRUE(opts$alleleFlipKriging) && nrow(df) >= 2L) {
    nKrIn <- nrow(df)
    R <- .ldFromSketch(ldSketch, df$SNP,
                       label = "summaryStatsQc: kriging prefilter")
    nKrig <- if (!is.null(opts$nForPip) && is.finite(opts$nForPip)) opts$nForPip
             else stats::median(as.numeric(df$N), na.rm = TRUE)
    kr <- krigingOutlierQc(df$Z, R, n = nKrig, variantIds = df$SNP)
    nKr <- sum(kr$flip)
    if (nKr > 0L) {
      df$Z[kr$flip] <- -df$Z[kr$flip]
      if ("BETA" %in% colnames(df)) df$BETA[kr$flip] <- -df$BETA[kr$flip]
    }
    entryAudit$krigingFlipped <- nKr
    entryAudit$krigingDiagnostics <- kr$diagnostics
    qcCount$krigingFlipped <- nKr
    .qcEmit(lbl, "QC track: kriging sign-flipped ", nKr, " of ", nKrIn,
         " LD-inconsistent variant(s).")
  }

  # 7. Optional LD-mismatch QC.
  if (!identical(opts$zMismatchQc, "none") && nrow(df) >= 2L) {
    nMmIn <- nrow(df)
    ldQc <- .applyLdMismatchQcToEntry(df, ldSketch, opts$zMismatchQc)
    df <- ldQc$df
    entryAudit$ldMismatchOutliersDropped <- ldQc$outliers
    entryAudit$ldMismatchMethod          <- opts$zMismatchQc
    # Preserve the full per-variant SLALOM/DENTIST diagnostics table for
    # downstream plotting / postprocessing (the outlier-only summary kept
    # historically dropped the per-variant detail). Stored as a data.frame
    # on the entry's audit; absent when zMismatchQc = "none".
    if (!is.null(ldQc$diagnostics))
      entryAudit$ldMismatchDiagnostics <- ldQc$diagnostics
    qcCount$mismatchRemoved <- ldQc$outliers
    .qcEmit(lbl, "QC track: ", opts$zMismatchQc, " removed ", ldQc$outliers, " of ",
         nMmIn, " LD-mismatch outlier(s).")
  }

  # 8. Optional RAISS imputation against the ldSketch.
  if (isTRUE(opts$impute) && nrow(df) >= 1L) {
    qcCount$imputeBefore <- nrow(df)
    # Scope the reference panel + dosage to the analysis-region window. RAISS
    # imputation is local (each missing variant is filled from its LD
    # neighbours), so a per-chromosome / genome-wide sketch must NOT materialize
    # its full dosage here -- restrict to [min(pos) - flank, max(pos) + flank]
    # computed PER CHROMOSOME (a multi-chromosome entry must not impute across a
    # cross-chromosome span), mirroring the region-scoping every other QC step
    # (harmonization, kriging, LD-mismatch) already does via df$SNP.
    flank    <- if (is.null(opts$imputeOpts$flank)) 0L
                else as.integer(opts$imputeOpts$flank)
    dfChrom  <- sub("^chr", "", as.character(df$chrom), ignore.case = TRUE)
    dfPos    <- as.integer(df$pos)
    loByChr  <- tapply(dfPos, dfChrom, min, na.rm = TRUE) - flank
    hiByChr  <- tapply(dfPos, dfChrom, max, na.rm = TRUE) + flank
    sketchSnpInfo <- getSnpInfo(ldSketch)
    skChrom   <- sub("^chr", "", as.character(sketchSnpInfo$CHR), ignore.case = TRUE)
    skBp      <- as.integer(sketchSnpInfo$BP)
    windowIdx <- which(skChrom %in% names(loByChr) &
                       skBp >= loByChr[skChrom] & skBp <= hiByChr[skChrom])

    if (length(windowIdx) == 0L) {
      .qcEmit(lbl, "QC track: RAISS imputation skipped (no LD-panel variants in the ",
           "region window).")
      entryAudit$raissImputedVariants <- 0L
      qcCount$imputeAfter <- nrow(df)
    } else {
      # Canonicalize the (windowed) panel ids so they align with the QC'd entry
      # ids (df$SNP is already canonical) across chr-prefix / separator
      # differences -- raiss() matches knownZscores$variant_id against
      # refPanel$variant_id by string.
      refPanel <- .refVariantsFromSketch(ldSketch)[windowIdx, , drop = FALSE]
      refPanel$variant_id <- normalizeVariantId(refPanel$variant_id)
      refPanel <- refPanel[order(refPanel$pos), , drop = FALSE]

      knownVariantIds <- if (!is.null(df$SNP)) as.character(df$SNP)
                         else as.character(df$variant_id)
      knownZ <- data.frame(
        chrom      = as.character(df$chrom),
        pos        = as.integer(df$pos),
        variant_id = knownVariantIds,
        A1         = as.character(df$A1),
        A2         = as.character(df$A2),
        z          = as.numeric(df$Z),
        stringsAsFactors = FALSE)
      if ("N"    %in% colnames(df)) knownZ$n    <- as.numeric(df$N)
      if ("BETA" %in% colnames(df)) knownZ$beta <- as.numeric(df$BETA)
      if ("SE"   %in% colnames(df)) knownZ$se   <- as.numeric(df$SE)
      knownZ <- knownZ[order(knownZ$pos), , drop = FALSE]

      # Materialize dosage for the region window only, in panel-order matching
      # refPanel.
      dosage <- .dosageMatrix(ldSketch, windowIdx, meanImpute = TRUE)
      colnames(dosage) <- normalizeVariantId(as.character(sketchSnpInfo$SNP[windowIdx]))
      dosage <- dosage[, refPanel$variant_id, drop = FALSE]
      scaledDosage <- scale(dosage)
      scaledDosage[is.na(scaledDosage)] <- 0

      imputed <- raiss(
        refPanel       = refPanel,
        knownZscores   = knownZ,
        genotypeMatrix = scaledDosage,
        svdTol         = if (is.null(opts$imputeOpts$svdTol)) 1e-12
                         else opts$imputeOpts$svdTol,
        lamb           = if (is.null(opts$imputeOpts$lamb)) 0.01
                         else opts$imputeOpts$lamb,
        r2Threshold    = if (is.null(opts$imputeOpts$r2Threshold)) 0.6
                         else opts$imputeOpts$r2Threshold,
        minimumLd      = if (is.null(opts$imputeOpts$minimumLd)) 5
                         else opts$imputeOpts$minimumLd,
        verbose        = FALSE)
      if (!is.null(imputed) && !is.null(imputed$resultFilter)) {
        impDf <- imputed$resultFilter
        out <- data.frame(
          chrom = impDf$chrom,
          pos   = impDf$pos,
          SNP   = impDf$variant_id,
          A1    = impDf$A1,
          A2    = impDf$A2,
          Z     = impDf$z,
          stringsAsFactors = FALSE)
        if ("n"    %in% colnames(impDf)) out$N    <- impDf$n
        if ("beta" %in% colnames(impDf)) out$BETA <- impDf$beta
        if ("se"   %in% colnames(impDf)) out$SE   <- impDf$se
        if ("N" %in% colnames(out) && any(is.na(out$N)))
          out$N[is.na(out$N)] <- stats::median(out$N, na.rm = TRUE)
        entryAudit$raissTotalVariants    <- nrow(out)
        entryAudit$raissImputedVariants  <- nrow(out) - nrow(knownZ)
        df <- out
      } else {
        entryAudit$raissImputedVariants <- 0L
      }
      qcCount$imputeAfter <- nrow(df)
      .qcEmit(lbl, "QC track: RAISS imputation ", qcCount$imputeBefore, " -> ",
           qcCount$imputeAfter, " variant(s) (net ",
           sprintf("%+d", qcCount$imputeAfter - qcCount$imputeBefore), ").")
    }
  }

  # Per-entry QC rollup: corrected (sign/strand flip, retained), removed
  # (drops at each step), imputed (added). Kept as distinct categories
  # because imputation adds variants back, so "in -> out" is not monotonic.
  # Skipped steps omitted.
  removedSegs <- character(0)
  sc <- entryAudit$sanityChecks
  if (!is.null(sc)) {
    if (!is.null(sc$nonstandardChrDropped) && sc$nonstandardChrDropped > 0L)
      removedSegs <- c(removedSegs,
                       paste0("nonstdChr ", sc$nonstandardChrDropped))
    if (!is.null(sc$missDataDropped) && sc$missDataDropped > 0L)
      removedSegs <- c(removedSegs,
                       paste0("missData ", sc$missDataDropped))
    if (!is.null(sc$pOutOfRangeDropped) && sc$pOutOfRangeDropped > 0L)
      removedSegs <- c(removedSegs,
                       paste0("badP ", sc$pOutOfRangeDropped))
    if (!is.null(sc$zeroEffectDropped) && sc$zeroEffectDropped > 0L)
      removedSegs <- c(removedSegs,
                       paste0("zeroEffect ", sc$zeroEffectDropped))
    if (!is.null(sc$nonpositiveSeDropped) && sc$nonpositiveSeDropped > 0L)
      removedSegs <- c(removedSegs,
                       paste0("badSE ", sc$nonpositiveSeDropped))
  }
  cf <- entryAudit$contentFilters
  if (!is.null(cf)) {
    if (!is.null(cf$mafDropped) && cf$mafDropped > 0L)
      removedSegs <- c(removedSegs, paste0("maf ", cf$mafDropped))
    if (!is.null(cf$infoDropped) && cf$infoDropped > 0L)
      removedSegs <- c(removedSegs, paste0("info ", cf$infoDropped))
    if (!is.null(cf$nDropped) && cf$nDropped > 0L)
      removedSegs <- c(removedSegs, paste0("nCutoff ", cf$nDropped))
  }
  if (qcCount$harmDropped > 0L)
    removedSegs <- c(removedSegs,
                     paste0("harmonization ", qcCount$harmDropped))
  if (!identical(opts$zMismatchQc, "none"))
    removedSegs <- c(removedSegs,
                     paste0("mismatch ", qcCount$mismatchRemoved))
  correctedSeg <- paste0("sign-flip ", qcCount$harmCorrSign,
                          ", strand-flip ", qcCount$harmCorrStrand)
  if (isTRUE(opts$alleleFlipKriging))
    correctedSeg <- paste0(correctedSeg,
                           ", kriging-flip ", qcCount$krigingFlipped)
  impSeg <- if (isTRUE(opts$impute) && !is.na(qcCount$imputeAfter)) {
    paste0(" | imputed ",
           sprintf("%+d", qcCount$imputeAfter - qcCount$imputeBefore))
  } else ""
  .qcEmit(lbl, "QC summary: ", nStudyIn, " in -> ", nrow(df), " out",
       " | corrected: ", correctedSeg,
       if (length(removedSegs) > 0L)
         paste0(" | removed: ", paste(removedSegs, collapse = ", "))
       else "",
       impSeg)

  entryAudit$variantsOut <- nrow(df)
  list(gr = .dfToEntryGranges(df), audit = entryAudit)
}

# Shrink an LD-sketch GenotypeHandle to the panel variants inside the summary
# statistics' per-chromosome position span. `entries` is a list/SimpleList of
# per-study (or per-tuple) GRanges. A genome-wide sketch otherwise carries a
# full-genome snpInfo; only variants inside [min,max] BP of each represented
# chromosome are reachable by harmonization or within-range imputation, so the
# rest is dropped at load time. NULL-safe; a no-op when the span already covers
# the panel. See [[.subsetGenotypeHandle]] for why this is read-safe.
# @noRd
.subsetSketchToRange <- function(ldSketch, entries) {
  if (is.null(ldSketch)) return(NULL)
  chrom <- unlist(lapply(entries, function(gr)
    canonChrom(as.character(GenomicRanges::seqnames(gr)))), use.names = FALSE)
  pos <- unlist(lapply(entries, function(gr)
    as.integer(GenomicRanges::start(gr))), use.names = FALSE)
  ok <- !is.na(chrom) & !is.na(pos)
  chrom <- chrom[ok]; pos <- pos[ok]
  if (length(pos) == 0L) return(ldSketch)
  lo <- tapply(pos, chrom, min)
  hi <- tapply(pos, chrom, max)
  si <- getSnpInfo(ldSketch)
  siChrom <- canonChrom(as.character(si$CHR))
  siBp    <- as.integer(si$BP)
  keep <- siChrom %in% names(lo) &
          siBp >= lo[siChrom] & siBp <= hi[siChrom]
  keep[is.na(keep)] <- FALSE
  .subsetGenotypeHandle(ldSketch, keep)
}

# Shrink an LD-sketch GenotypeHandle to EXACTLY the variants present across the
# QC'd `entries` (imputation may have added variants; QC may have dropped some),
# matched by canonical variant id. Applied at the end of summaryStatsQc so the
# retained sketch mirrors the object's final variant set. NULL-safe.
# @noRd
.subsetSketchToIds <- function(ldSketch, entries) {
  if (is.null(ldSketch)) return(NULL)
  ids <- unlist(lapply(entries, function(gr) {
    if (is.null(gr)) return(character(0))
    snp <- S4Vectors::mcols(gr)$SNP
    if (is.null(snp)) character(0) else as.character(snp)
  }), use.names = FALSE)
  if (length(ids) == 0L) return(ldSketch)
  si <- getSnpInfo(ldSketch)
  keep <- normalizeVariantId(as.character(si$SNP)) %in%
          normalizeVariantId(unique(ids))
  .subsetGenotypeHandle(ldSketch, keep)
}

#' Run QC on a SumStats Collection
#'
#' Applies a single QC pass to a \code{QtlSumStats} or \code{GwasSumStats}
#' collection: per-row sanity checks via \code{.applySanityChecks} (drop
#' rows with out-of-range / zero P, BETA == 0, SE <= 0, NA in vital
#' columns; clamp tiny P; normalize CHR; coerce signed columns to
#' numeric), variant-content filters (MAF / INFO / N) via
#' \code{.applyContentFilters}, optional \code{skipRegion} drop, optional
#' PIP screen, panel-vs-sumstats allele harmonization against the
#' \code{ldSketch} via \code{harmonizeAlleles} (which handles indels,
#' strand-ambiguous variants, sign / strand flips, and duplicate drops in
#' a single sweep), optional SLALOM/DENTIST LD-mismatch QC, and optional
#' RAISS imputation. No Bioconductor genome / dbSNP packages required.
#'
#' The returned collection has its \code{qcInfo} slot populated with a
#' per-entry audit record (variant counts, drop counts at each step,
#' which filters fired, etc.). Fine-mapping and TWAS-weights pipelines
#' reject SumStats inputs where \code{length(getQcInfo(x)) == 0L}.
#'
#' Column-availability error contract: a non-zero \code{mafCutoff}
#' requires every entry to carry a \code{MAF} column; non-zero
#' \code{infoCutoff} requires \code{INFO}; non-zero \code{nCutoff}
#' requires \code{N}. Missing column with a non-zero cutoff is a hard
#' error.
#'
#' @param sumstats A \code{QtlSumStats} or \code{GwasSumStats}
#'   collection.
#' @param removeIndels Logical (length 1). When \code{TRUE}, drop
#'   indels during panel harmonization. Default \code{FALSE}.
#' @param removeStrandAmbiguous Logical (length 1). When \code{TRUE},
#'   drop A/T and C/G strand-ambiguous variants. Default \code{TRUE}.
#' @param mafCutoff Numeric (length 1). MAF threshold (variants with
#'   \code{MAF < mafCutoff} are dropped). Default 0. Requires \code{MAF}
#'   or \code{FRQ} column when non-zero.
#' @param infoCutoff Numeric (length 1). INFO score threshold. Default
#'   0. Requires \code{INFO} column when non-zero.
#' @param nCutoff Numeric (length 1). Sample-size deviation threshold:
#'   drop variants whose \code{N} is more than \code{nCutoff}
#'   median-absolute-deviations from the median. Set to 0 to disable.
#'   Default 5.
#' @param keepVariants Optional character vector of variant IDs (SNP
#'   column) to retain prior to harmonization.
#' @param skipRegion Optional character vector of \code{"chr:start-end"}
#'   strings, or a \code{GRanges}, of regions to drop.
#' @param pipCutoffToSkip Numeric (length 1). When \code{!= 0}, run an
#'   LD-independent single-effect SER screen and skip the entry if no
#'   PIP exceeds the cutoff. \code{< 0} resolves to \code{3 / nVariants}.
#'   Default 0 (no screen).
#' @param absZCutoffToSkip Numeric (length 1). Alternative signal screen:
#'   skip the entry when \code{max(abs(Z))} does not exceed the cutoff. No
#'   model fit. Default 0 (off).
#' @param bfCutoffToSkip Numeric (length 1). Alternative signal screen: skip
#'   the entry when the largest per-variant single-effect Bayes factor (from
#'   the same \code{susie_ser} fit as the PIP screen) does not exceed the
#'   cutoff. Compared in log space (\code{maxlogBF > log(cutoff)}). Must be
#'   \code{> 0}. Default 0 (off).
#' @param logBfCutoffToSkip Numeric (length 1). As \code{bfCutoffToSkip} but
#'   the cutoff is on the log Bayes factor scale (\code{maxlogBF > cutoff}).
#'   Default 0 (off).
#'   Exactly one of \code{pipCutoffToSkip} / \code{absZCutoffToSkip} /
#'   \code{bfCutoffToSkip} / \code{logBfCutoffToSkip} may be non-zero (one
#'   screening metric at a time); enabling more than one is an error.
#' @param zMismatchQc One of \code{"none"} (default), \code{"slalom"},
#'   \code{"dentist"}.
#' @param alleleFlipKriging Logical (length 1). Opt-in kriging
#'   LD-consistency prefilter run before SLALOM/DENTIST. Default
#'   \code{FALSE}.
#' @param effectiveN Logical (length 1). When \code{TRUE} (default) and the
#'   input carries case/control counts --- per-variant \code{N_CASE} /
#'   \code{N_CONTROL} mcols, else the study-level \code{nCase} /
#'   \code{nControl} scalars --- the working per-variant \code{N} is set to
#'   the effective sample size \code{effectiveN(nCase, nControl)} BEFORE the
#'   N-cutoff filter, so the filter, kriging, and the downstream fit all use
#'   \code{N_eff}. When both counts and an \code{N} column are present the
#'   counts win: \code{N} is overridden and the override is logged. Inputs
#'   with no counts (quantitative traits) are unchanged. The escape hatch
#'   \code{effectiveN = FALSE} restores the raw \code{N} column (or, when
#'   there is no \code{N}, the raw total \code{nCase + nControl}) with no
#'   override. \code{qcInfo$options$effectiveN} records the setting and each
#'   entry's \code{nSource} is one of \code{"effective"}, \code{"column"},
#'   \code{"total"}, or \code{NA}.
#' @param impute Logical (length 1). Run RAISS imputation against the
#'   \code{ldSketch}. Default \code{FALSE}. (Note: RAISS against the
#'   sketch is not yet fully wired for the new path; the option is
#'   accepted but currently emits a warning and is skipped.)
#' @param imputeOpts Named list of RAISS parameters. RAISS imputation scopes its
#'   reference panel to the analysis-region window (so a per-chromosome /
#'   genome-wide \code{ldSketch} does not materialize its full dosage);
#'   \code{flank} (default 0) widens that window by the given number of base
#'   pairs on each side to retain LD context for edge variants.
#' @param matchMinProp Minimum proportion of LD panel variants that must
#'   be matched by the sumstats; default 0.
#' @param coerceNumeric Logical. Coerce signed columns
#'   (Z/BETA/SE/OR/LOG_ODDS/SIGNED_SUMSTAT/P/MAF/FRQ/INFO/N) to numeric.
#'   Default \code{TRUE}.
#' @param normalizeChr Logical. Strip \code{"chr"} prefix, uppercase the
#'   chromosome label, and map 23->X, 24->Y, M->MT. Default \code{TRUE}.
#' @param dropNonstandardChr Logical. Drop variants whose CHR (after
#'   normalization) is outside 1..22, X, Y, MT. Default \code{TRUE}.
#' @param dropMissData Logical. Drop rows with NA in any vital column
#'   (chrom, pos, A1, A2, and at least one of Z / BETA). Default
#'   \code{TRUE}.
#' @param dropPOutOfRange Logical. Drop rows where \code{P < 0} or
#'   \code{P > 1}. Default \code{TRUE}.
#' @param clampSmallP Logical. Floor non-negative P values below
#'   \code{smallPFloor} to \code{smallPFloor} so \code{-log10(P)} stays
#'   finite. Applied to both input and Z-derived P values. Default
#'   \code{TRUE}.
#' @param smallPFloor Numeric (length 1). Floor for \code{clampSmallP}.
#'   Default \code{5e-324} (R's smallest positive double).
#' @param dropZeroEffect Logical. Drop rows where any effect column is
#'   exactly 0 (\code{BETA}, \code{LOG_ODDS}, \code{SIGNED_SUMSTAT}) or
#'   \code{OR} is exactly 1. Default \code{TRUE}.
#' @param dropNonpositiveSe Logical. Drop rows where \code{SE <= 0}.
#'   Default \code{TRUE}.
#' @return A new \code{QtlSumStats} / \code{GwasSumStats} with cleaned
#'   entries and \code{qcInfo} populated.
#' @export
summaryStatsQc <- function(sumstats,
                           removeIndels           = FALSE,
                           removeStrandAmbiguous  = TRUE,
                           mafCutoff              = 0,
                           infoCutoff             = 0,
                           nCutoff                = 5,
                           keepVariants           = NULL,
                           skipRegion             = NULL,
                           pipCutoffToSkip        = 0,
                           absZCutoffToSkip       = 0,
                           bfCutoffToSkip         = 0,
                           logBfCutoffToSkip      = 0,
                           zMismatchQc            = c("none", "slalom",
                                                     "dentist"),
                           alleleFlipKriging      = FALSE,
                           effectiveN             = TRUE,
                           impute                 = FALSE,
                           imputeOpts             = list(rcond = 0.01,
                                                        r2Threshold = 0.6,
                                                        minimumLd = 5,
                                                        lamb = 0.01),
                           matchMinProp           = 0,
                           coerceNumeric          = TRUE,
                           normalizeChr           = TRUE,
                           dropNonstandardChr     = TRUE,
                           dropMissData           = TRUE,
                           dropPOutOfRange        = TRUE,
                           clampSmallP            = TRUE,
                           smallPFloor            = 5e-324,
                           dropZeroEffect         = TRUE,
                           dropNonpositiveSe      = TRUE) {
  if (!methods::is(sumstats, "QtlSumStats") &&
      !methods::is(sumstats, "GwasSumStats")) {
    stop("summaryStatsQc requires a QtlSumStats or GwasSumStats input.")
  }
  zMismatchQc <- match.arg(zMismatchQc)

  # Column-availability checks across all entries.
  for (i in seq_len(nrow(sumstats))) {
    mc <- S4Vectors::mcols(sumstats$entry[[i]])
    cols <- colnames(mc)
    if (mafCutoff > 0 && !any(c("MAF", "FRQ") %in% cols))
      stop("summaryStatsQc: mafCutoff > 0 requires every entry to carry a ",
           "MAF or FRQ column; entry ", i, " does not.")
    if (infoCutoff > 0 && !"INFO" %in% cols)
      stop("summaryStatsQc: infoCutoff > 0 requires every entry to carry an ",
           "INFO column; entry ", i, " does not.")
  }

  opts <- list(
    removeIndels           = removeIndels,
    removeStrandAmbiguous  = removeStrandAmbiguous,
    mafCutoff              = mafCutoff,
    infoCutoff             = infoCutoff,
    nCutoff                = nCutoff,
    keepVariants           = as.character(keepVariants),
    skipRegion             = skipRegion,
    screen                 = .resolveScreenMetric(pipCutoffToSkip,
                                                  absZCutoffToSkip,
                                                  bfCutoffToSkip,
                                                  logBfCutoffToSkip),
    zMismatchQc            = zMismatchQc,
    alleleFlipKriging      = alleleFlipKriging,
    effectiveN             = effectiveN,
    nCase                  = NULL,
    nControl               = NULL,
    impute                 = impute,
    imputeOpts             = imputeOpts,
    matchMinProp           = matchMinProp,
    coerceNumeric          = coerceNumeric,
    normalizeChr           = normalizeChr,
    dropNonstandardChr     = dropNonstandardChr,
    dropMissData           = dropMissData,
    dropPOutOfRange        = dropPOutOfRange,
    clampSmallP            = clampSmallP,
    smallPFloor            = smallPFloor,
    dropZeroEffect         = dropZeroEffect,
    dropNonpositiveSe      = dropNonpositiveSe,
    nForPip                = NULL)

  newEntries <- vector("list", nrow(sumstats))
  entryAudits <- vector("list", nrow(sumstats))
  isQtl <- methods::is(sumstats, "QtlSumStats")
  for (i in seq_len(nrow(sumstats))) {
    opts$nForPip <- if ("N" %in% colnames(S4Vectors::mcols(sumstats$entry[[i]])))
      stats::median(S4Vectors::mcols(sumstats$entry[[i]])$N, na.rm = TRUE)
    else NULL
    # Per-entry study-level case/control counts (used only when the entry has
    # no per-variant N_CASE / N_CONTROL columns). NULL for quantitative traits.
    opts$nCase    <- if ("nCase"    %in% names(sumstats)) as.numeric(sumstats$nCase)[[i]]    else NULL
    opts$nControl <- if ("nControl" %in% names(sumstats)) as.numeric(sumstats$nControl)[[i]] else NULL
    # Study-level total N: the level-4 fallback used when the entry has no
    # per-variant N_CASE/N_CONTROL, no study nCase/nControl, and no N column.
    opts$nSample  <- if ("nSample"  %in% names(sumstats)) as.numeric(sumstats$nSample)[[i]]  else NULL
    # Per-entry label woven into QC log messages and the rollup. For
    # QtlSumStats it's (study/context/trait); for GwasSumStats it's the
    # study identifier.
    entryLabel <- if (isQtl) {
      paste(as.character(sumstats$study)[[i]],
            as.character(sumstats$context)[[i]],
            as.character(sumstats$trait)[[i]], sep = "/")
    } else {
      as.character(sumstats$study)[[i]]
    }
    result <- .runEntrySummaryStatsQc(
      gr         = sumstats$entry[[i]],
      ldSketch   = getLdSketch(sumstats),
      refGenome  = getGenome(sumstats),
      opts       = opts,
      entryLabel = entryLabel)
    newEntries[[i]] <- result$gr
    entryAudits[[i]] <- result$audit
  }

  qcInfo <- list(
    timestamp        = NA_character_,
    options          = list(
      removeIndels          = removeIndels,
      removeStrandAmbiguous = removeStrandAmbiguous,
      mafCutoff             = mafCutoff,
      infoCutoff            = infoCutoff,
      nCutoff               = nCutoff,
      pipCutoffToSkip       = pipCutoffToSkip,
      absZCutoffToSkip      = absZCutoffToSkip,
      bfCutoffToSkip        = bfCutoffToSkip,
      logBfCutoffToSkip     = logBfCutoffToSkip,
      zMismatchQc           = zMismatchQc,
      alleleFlipKriging     = alleleFlipKriging,
      effectiveN            = effectiveN,
      impute                = impute,
      coerceNumeric         = coerceNumeric,
      normalizeChr          = normalizeChr,
      dropNonstandardChr    = dropNonstandardChr,
      dropMissData          = dropMissData,
      dropPOutOfRange       = dropPOutOfRange,
      clampSmallP           = clampSmallP,
      smallPFloor           = smallPFloor,
      dropZeroEffect        = dropZeroEffect,
      dropNonpositiveSe     = dropNonpositiveSe),
    entryAudit       = entryAudits)

  # Trim the LD sketch to the variants actually present after QC (imputation may
  # have added some; filters removed others) so it mirrors the final object.
  newLdSketch <- .subsetSketchToIds(getLdSketch(sumstats), newEntries)

  # Rebuild the SumStats with new entries and qcInfo.
  if (methods::is(sumstats, "GwasSumStats")) {
    GwasSumStats(
      study    = as.character(sumstats$study),
      entry    = newEntries,
      genome   = getGenome(sumstats),
      ldSketch = newLdSketch,
      varY     = as.numeric(sumstats$varY),
      # Preserve the optional per-study case/control counts through QC.
      nCase    = if ("nCase" %in% names(sumstats))
                   as.numeric(sumstats$nCase) else NULL,
      nControl = if ("nControl" %in% names(sumstats))
                   as.numeric(sumstats$nControl) else NULL,
      nSample  = if ("nSample" %in% names(sumstats))
                   as.numeric(sumstats$nSample) else NULL,
      qcInfo   = qcInfo)
  } else {
    QtlSumStats(
      study    = as.character(sumstats$study),
      context  = as.character(sumstats$context),
      trait    = as.character(sumstats$trait),
      entry    = newEntries,
      genome   = getGenome(sumstats),
      ldSketch = newLdSketch,
      varY     = as.numeric(sumstats$varY),
      # Preserve the optional per-tuple study-level total N through QC.
      nSample  = if ("nSample" %in% names(sumstats))
                   as.numeric(sumstats$nSample) else NULL,
      qcInfo   = qcInfo)
  }
}
