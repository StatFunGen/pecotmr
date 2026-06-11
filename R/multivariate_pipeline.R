#' Convert QC'ed regional summary-statistic data to mvSuSiE RSS input
#'
#' @param sumstatData The \code{sumstat_data} component from a QC'ed regional
#'   object, typically \code{qced_regional_data$sumstat_data}.
#' @param ldName Optional name of the LD reference to use. If \code{NULL}, a
#'   unique \code{LD_match} entry is used. If \code{LD_match} points to multiple
#'   references, the first one is used with a message. When no \code{LD_match}
#'   is available, the first LD reference containing all shared variants is used.
#' @return A list with \code{mvsusie_rss_input}, ready for
#'   \code{mvsusieR::mvsusie_rss()}, and \code{source_info}.
#' @export
regionDataToMvsusieRssInput <- function(sumstatData, ldName = NULL) {
  if (is.null(sumstatData) || is.null(sumstatData$sumstats) ||
      is.null(sumstatData$LD_data)) {
    stop("sumstatData must contain post-QC sumstats and LD_data entries.")
  }

  sumstats <- sumstatData$sumstats
  if (length(sumstats) < 2L) {
    stop("mvSuSiE RSS input requires at least two summary-statistic studies.")
  }

  studyNames <- names(sumstats)
  if (is.null(studyNames) || any(!nzchar(studyNames))) {
    studyNames <- paste0("study", seq_along(sumstats))
    names(sumstats) <- studyNames
  }

  sumstats <- lapply(sumstats, function(record) {
    record$sumstats <- .rss_sumstats_with_variant_id(record$sumstats)
    record
  })
  overlap <- Reduce(intersect, lapply(sumstats, function(record) record$sumstats$variant_id))
  if (length(overlap) < 2L) {
    stop("mvSuSiE RSS input requires at least two shared variants across studies.")
  }

  Z <- do.call(cbind, lapply(sumstats, function(record) {
    record$sumstats$z[match(overlap, record$sumstats$variant_id)]
  }))
  rownames(Z) <- overlap
  colnames(Z) <- studyNames

  nVec <- vapply(sumstats, function(record) {
    n <- record$n
    if (is.null(n) || length(n) == 0L) return(NA_real_)
    stats::median(as.numeric(n), na.rm = TRUE)
  }, numeric(1))
  names(nVec) <- studyNames

  LD_data <- sumstatData$LD_data
  ldNames <- names(LD_data)
  if (is.null(ldNames) || any(!nzchar(ldNames))) {
    ldNames <- paste0("LD", seq_along(LD_data))
    names(LD_data) <- ldNames
  }

  LD_match <- sumstatData$LD_match
  matchedLdNames <- character()
  if (!is.null(LD_match)) {
    LD_match <- as.character(LD_match)
    if (!is.null(names(LD_match)) && all(studyNames %in% names(LD_match))) {
      matchedLdNames <- LD_match[studyNames]
    } else if (length(LD_match) >= length(studyNames)) {
      matchedLdNames <- LD_match[seq_along(studyNames)]
    }
    matchedLdNames <- unique(matchedLdNames[!is.na(matchedLdNames) & nzchar(matchedLdNames)])
  }

  selectedLdName <- ldName
  if (is.null(selectedLdName) && length(matchedLdNames) == 1L) {
    if (!matchedLdNames %in% ldNames) {
      stop("LD_match points to an LD reference that is not present in sumstatData$LD_data.")
    }
    selectedLdName <- matchedLdNames
  }
  if (is.null(selectedLdName) && length(matchedLdNames) > 1L) {
    selectedLdName <- matchedLdNames[[1]]
    message("mvSuSiE RSS input: multiple LD_match references were found; using the first reference '",
            selectedLdName, "'. Provide ldName to choose a different reference.")
  }
  if (is.null(selectedLdName)) {
    containsOverlap <- vapply(LD_data, function(ld) {
      is(ld, "LdData") && all(overlap %in% getVariantIds(ld))
    }, logical(1))
    if (!any(containsOverlap)) {
      stop("No LD_data entry contains all shared variants.")
    }
    selectedLdName <- ldNames[which(containsOverlap)[1]]
  }
  if (!selectedLdName %in% ldNames) {
    stop("ldName is not present in sumstatData$LD_data.")
  }

  ld <- LD_data[[selectedLdName]]
  if (!is(ld, "LdData")) stop("Selected LD_data entry must be an LdData object.")

  referenceIds <- getVariantIds(ld)
  if (hasGenotypes(ld)) {
    X <- getGenotypes(ld)
    if (is.list(X) && !is.matrix(X)) {
      stop("regionDataToMvsusieRssInput requires a single genotype reference, not mixture panels.")
    }
    Xoverlap <- .subset_rss_matrix_columns(
      X, overlap, referenceIds, "genotype reference panel"
    )
    R <- computeLd(Xoverlap, method = "sample")
    rownames(R) <- colnames(R) <- overlap
  } else {
    R <- ld@correlation
    if (is.null(R) || (is.list(R) && !is.matrix(R))) {
      stop("regionDataToMvsusieRssInput requires one correlation matrix or one genotype matrix.")
    }
    R <- .subset_rss_ld_matrix(R, overlap, referenceIds)
  }

  list(
    mvsusie_rss_input = list(
      Z = Z,
      R = R,
      N = if (all(is.na(nVec))) NA_real_ else max(nVec, na.rm = TRUE)
    ),
    source_info = list(
      studies = studyNames,
      variants = overlap,
      n = nVec,
      ld_name = selectedLdName
    )
  )
}

#' Multivariate Analysis Pipeline
#'
#' This function performs weights computation for Transcriptome-Wide Association Study (TWAS) with fitting
#' models using mvSuSiE and mr.mash with the option of using a limited number of variants selected from
#' mvSuSiE fine-mapping for computing TWAS weights with cross-validation.
#'
#' @param X A matrix of genotype data where rows represent samples and columns represent genetic variants.
#' @param Y A matrix of phenotype measurements, representing samples and columns represent conditions.
#' @param maf Optional vector of minor allele frequencies for each variant in X,
#'   used ONLY for \code{mafCutoff} filtering and never exported. When \code{af}
#'   is supplied the filtering MAF is derived from it (\code{min(af, 1 - af)}) and
#'   a supplied \code{maf} is ignored (with a warning if they disagree). Default
#'   NULL; if neither \code{maf} nor \code{af} is supplied and \code{mafCutoff}
#'   is set, the call errors.
#' @param af Optional vector of directional effect-allele frequencies (frequency
#'   of \code{a1}) aligned to the columns of X. When supplied it is exported as
#'   the \code{top_loci$af} column; when NULL, \code{af} is \code{NA_real_}.
#'   Default NULL.
#' @param L Maximum number of components in mvSuSiE. Default is 30.
#' @param L_greedy Initial greedy number of components in mvSuSiE. Default is 5.
#' @param ld_reference_meta_file An optional path to a file containing linkage disequilibrium reference data. If provided, variants in X are filtered based on this reference.
#' @param pip_cutoff_to_skip Cutoff value for skipping conditions based on PIP values. Default is 0.
#' @param signal_cutoff Cutoff value for signal identification in PIP values. Default is 0.025.
#' @param coverage A vector of coverage probabilities, with the first element being the primary coverage and the rest being secondary coverage probabilities for credible set refinement. Defaults to c(0.95, 0.7, 0.5).
#' @param min_abs_corr Minimum absolute correlation for credible set purity filtering. Default is 0.8,
#'   which is stricter than the susieR default of 0.5.
#' @param data_driven_prior_matrices A list of data-driven covariance matrices for mr.mash weights.
#' @param data_driven_prior_matrices_cv A list of data-driven covariance matrices for mr.mash weights in cross-validation.
#' @param canonical_prior_matrices If set to TRUE, will compute canonical covariance matrices and add them into the prior covariance matrix list in mrmash_wrapper. Default is TRUE.
#' @param sample_partition Optional data frame with Sample and Fold columns for cross-validation.
#' @param mrmash_max_iter The maximum number of iterations for mr.mash. Default is 5000.
#' @param mvsusie_max_iter The maximum number of iterations for mvSuSiE. Default is 200.
#' @param estimate_residual_variance Passed to \code{mvsusieR::mvsusie()}. Default is TRUE.
#' @param min_cv_maf The minimum minor allele frequency for variants to be included in cross-validation. Default is 0.05.
#' @param max_cv_variants The maximum number of variants to be included in cross-validation. Defaults to -1 which means no limit.
#' @param cv_folds The number of folds to use for cross-validation. Set to 0 to skip cross-validation. Default is 5.
#' @param cv_threads The number of threads to use for parallel computation in cross-validation. Defaults to 1.
#' @param data_driven_prior_weights_cutoff The minimum weight for prior covariance matrices. Default is 1e-4.
#' @param verbose Verbosity level. Default is 0.
#'
#' @return A list containing the multivariate analysis results.
#' @examples
#' library(pecotmr)
#'
#' data(multitrait_data)
#' attach(multitrait_data)
#'
#' data_driven_prior_matrices <- list(
#'   U = prior_matrices,
#'   w = rep(1 / length(prior_matrices), length(prior_matrices))
#' )
#'
#' data_driven_prior_matrices_cv <- lapply(prior_matrices_cv, function(x) {
#'   list(U = x, w = rep(1 / length(x), length(x)))
#' })
#'
#' result <- multivariateAnalysisPipeline(
#'   X = multitrait_data$X,
#'   Y = multitrait_data$Y,
#'   maf = colMeans(multitrait_data$X),
#'   L = 10,
#'   lGreedy = 5,
#'   ldReferenceMetaFile = NULL,
#'   maxCvVariants = -1,
#'   pipCutoffToSkip = 0,
#'   signalCutoff = 0.025,
#'   dataDrivenPriorMatrices = dataDrivenPriorMatrices,
#'   dataDrivenPriorMatricesCv = dataDrivenPriorMatricesCv,
#'   canonicalPriorMatrices = TRUE,
#'   samplePartition = NULL,
#'   cvFolds = 5,
#'   cvThreads = 2,
#'   dataDrivenPriorWeightsCutoff = 1e-4
#' )
#' @export
multivariateAnalysisPipeline <- function(
    # input data
    X,
    Y,
    maf = NULL,
    af = NULL,
    xVariance = NULL,
    otherQuantities = list(),
    region = NULL,
    # filters
    imissCutoff = 1.0,
    mafCutoff = 0.01,
    xvarCutoff = 0.01,
    ldReferenceMetaFile = NULL,
    pipCutoffToSkip = 0,
    # methods parameter configuration
    L = 30,
    lGreedy = 5,
    dataDrivenPriorMatrices = NULL,
    dataDrivenPriorMatricesCv = NULL,
    dataDrivenPriorWeightsCutoff = 1e-4,
    canonicalPriorMatrices = TRUE,
    mrmashMaxIter = 5000,
    mvsusieMaxIter = 200,
    estimateResidualVariance = TRUE,
    # fine-mapping results summary
    signalCutoff = 0.025,
    coverage = c(0.95, 0.7, 0.5),
    minAbsCorr = 0.8,
    # TWAS weights and CV for TWAS weights
    twasWeights = TRUE,
    samplePartition = NULL,
    maxCvVariants = -1,
    cvFolds = 5,
    cvThreads = 1,
    verbose = 0) {
  # Make sure mvsusieR is installed
  if (!requireNamespace("mvsusieR", quietly = TRUE)) {
    stop("To use this function, please install mvsusieR: https://github.com/stephenslab/mvsusieR")
  }
  # Skip conditions based on univariate PIP values
  skipConditions <- function(X, Y, pipCutoffToSkip) {
    if (length(pipCutoffToSkip) == 1 && is.numeric(pipCutoffToSkip)) {
      pipCutoffToSkip <- rep(pipCutoffToSkip, ncol(Y))
    } else if (length(pipCutoffToSkip) != ncol(Y)) {
      stop("pipCutoffToSkip must be a single number or a vector of the same length as ncol(Y).")
    }
    colsToKeep <- logical(ncol(Y))
    for (r in 1:ncol(Y)) {
      if (pipCutoffToSkip[r] != 0) {
        nonMissingIndices <- which(!is.na(Y[, r]))
        XNonMissing <- X[match(names(Y[, r])[nonMissingIndices], rownames(X)), ]
        YNonMissing <- Y[nonMissingIndices, r]
        if (pipCutoffToSkip[r] < 0) {
          # automatically determine the cutoff to use
          pipCutoffToSkip[r] <- 3 * 1 / ncol(XNonMissing)
        }
        topModelPip <- susie(XNonMissing, YNonMissing, L = 1)$pip

        if (any(topModelPip > pipCutoffToSkip[r])) {
          colsToKeep[r] <- TRUE
        } else {
          message(paste0(
            "Skipping condition ", colnames(Y)[r], ", because all top_model_pip < pipCutoffToSkip = ",
            pipCutoffToSkip[r], ". Top loci model does not show any potentially significant variants."
          ))
        }
      } else {
        colsToKeep[r] <- TRUE
      }
    }

    YFiltered <- Y[, colsToKeep, drop = FALSE]

    if (ncol(YFiltered) <= 1) {
      warning("After filtering by potential association signals, Y has ", ncol(YFiltered), " context left. Returning NULL.")
      return(NULL)
    } else {
      message("After filtering by potential association signals, Y has ", ncol(YFiltered), " contexts left.")
      return(YFiltered)
    }
  }

  initializeMvsusiePrior <- function(conditionNames, dataDrivenPriorMatrices,
                                     dataDrivenPriorMatricesCv, cvFolds, priorWeights, dataDrivenPriorWeightsCutoff) {
    if (!is.null(dataDrivenPriorMatrices)) {
      # update w based on mrmash prior weights
      message("Updating prior weights based on mrmash_fitted. ")
      dataDrivenPriorMatrices$w <- priorWeights
      dataDrivenPriorMatrices$U <- dataDrivenPriorMatrices$U[names(priorWeights)]
      dataDrivenPriorMatrices <- list(matrices = dataDrivenPriorMatrices$U, weights = dataDrivenPriorMatrices$w)
      dataDrivenPriorMatrices <- mvsusieR::create_mixture_prior(mixture_prior = dataDrivenPriorMatrices, weights_tol = dataDrivenPriorWeightsCutoff, include_indices = conditionNames)
    } else {
      dataDrivenPriorMatrices <- mvsusieR::create_mixture_prior(R = length(conditionNames), include_indices = conditionNames)
    }

    if (!is.null(dataDrivenPriorMatricesCv)) {
      dataDrivenPriorMatricesCv <- lapply(
        dataDrivenPriorMatricesCv,
        function(x) {
          x$U <- x$U[names(priorWeights)]
          x <- list(matrices = x$U, weights = priorWeights)
          mvsusieR::create_mixture_prior(mixture_prior = x, weights_tol = dataDrivenPriorWeightsCutoff, include_indices = conditionNames)
        }
      )
    } else {
      if (!is.null(dataDrivenPriorMatrices)) {
        dataDrivenPriorMatricesCv <- lapply(1:cvFolds, function(x) {
          return(dataDrivenPriorMatrices)
        })
      }
    }
    return(list(
      dataDrivenPriorMatrices = dataDrivenPriorMatrices, dataDrivenPriorMatricesCv = dataDrivenPriorMatricesCv
    ))
  }

  # filter X and Y missing, specific to multivariate analysis where some conditions are skipped we have to updated X matrix
  filterXYMissing <- function(X, Y) {
    YRowsWithMissing <- apply(Y, 1, function(row) all(is.na(row)))
    if (any(YRowsWithMissing)) {
      YFiltered <- Y[-which(YRowsWithMissing), , drop = FALSE]
    } else {
      YFiltered <- Y
    }
    XFiltered <- X[match(rownames(YFiltered), rownames(X)), ]
    XColumnsWithMissing <- apply(XFiltered, 2, function(column) all(is.na(column)))
    if (any(XColumnsWithMissing)) {
      columnsToRemove <- which(XColumnsWithMissing)
      XFiltered <- XFiltered[, -columnsToRemove, drop = FALSE]
    }
    return(list(XFiltered = XFiltered, YFiltered = YFiltered))
  }

  # Input validation
  if (!is.matrix(X) || !is.numeric(X)) stop("X must be a numeric matrix")
  if (!is.matrix(Y) || !is.numeric(Y)) stop("Y must be a numeric matrix")
  if (nrow(X) != nrow(Y)) stop("X and Y must have the same number of rows")
  if (!is.null(maf)) {
    if (!is.numeric(maf) || length(maf) != ncol(X)) stop("maf must be NULL or a numeric vector with length equal to the number of columns in X")
    if (any(maf < 0 | maf > 1, na.rm = TRUE)) stop("maf values must be between 0 and 1")
  }
  if (!is.null(af)) {
    if (!is.numeric(af) || length(af) != ncol(X)) stop("af must be NULL or a numeric vector with length equal to the number of columns in X")
    if (any(af < 0 | af > 1, na.rm = TRUE)) stop("af values must be between 0 and 1")
  }
  # Single source of truth = af. When af is available, derive the filtering MAF
  # from it (min(af, 1 - af)); a supplied directionless maf is only a fallback
  # and, if it disagrees with the af-derived value, af wins (with a warning).
  if (!is.null(af)) {
    afDerivedMaf <- pmin(af, 1 - af)
    if (!is.null(maf) && any(abs(maf - afDerivedMaf) > 1e-6, na.rm = TRUE)) {
      warning("Both 'maf' and 'af' were supplied and disagree; using the ",
              "af-derived MAF for filtering (af is the single source of truth).")
    }
    maf <- afDerivedMaf
  }
  if (is.null(maf) && !is.null(mafCutoff) && is.numeric(mafCutoff) && mafCutoff > 0) {
    stop("mafCutoff is set but neither 'af' nor 'maf' was supplied; provide ",
         "one so MAF can be derived for filtering.")
  }
  if (!is.numeric(L) || L <= 0) stop("L must be a positive integer")
  if (!is.null(lGreedy) && (!is.numeric(lGreedy) || lGreedy <= 0)) stop("lGreedy must be NULL or a positive integer")
  if (!is.null(lGreedy)) lGreedy <- min(lGreedy, L)

  # main analysis codes
  Y <- skipConditions(X, Y, pipCutoffToSkip)
  if (is.null(Y)) {
    return(list())
  }

  # filter X and Y missing data
  XYFiltered <- filterXYMissing(X, Y)
  X <- XYFiltered$XFiltered
  Y <- XYFiltered$YFiltered
  if (nrow(Y) == 0 || is.null(Y)) {
    return(list())
  }

  # filter variants by ld reference panel
  if (!is.null(ldReferenceMetaFile)) {
    variantsKept <- filterVariantsByLdReference(colnames(X), ldReferenceMetaFile)
    X <- X[, variantsKept$data, drop = FALSE]
    if (!is.null(maf)) maf <- maf[variantsKept$idx]
    if (!is.null(af)) af <- af[variantsKept$idx]
  }

  # filter X based on Y subjects
  if (!is.null(imissCutoff) || !is.null(mafCutoff)) {
    X <- filterXWithY(X, Y, imissCutoff, mafCutoff, varThresh = xvarCutoff, maf = maf, xVariance = xVariance)
    if (!is.null(maf)) maf <- maf[colnames(X)]
    if (!is.null(af)) af <- af[colnames(X)]
  }

  # filter data driven prior matrices
  if (!is.null(dataDrivenPriorMatrices)) {
    dataDrivenPriorMatrices <- filterMixtureComponents(
      colnames(Y),
      dataDrivenPriorMatrices$U, dataDrivenPriorMatrices$w,
      dataDrivenPriorWeightsCutoff
    )
  }

  st <- proc.time()
  res <- list()
  message("Fitting mr.mash model on input data ...")
  res$mrmash_fitted <- mrmashWrapper(
    X = X, Y = Y, dataDrivenPriorMatrices = dataDrivenPriorMatrices,
    canonicalPriorMatrices = canonicalPriorMatrices, maxIter = mrmashMaxIter
  )

  # For input into mvSuSiE
  residY <- res$mrmash_fitted$V
  w0Updated <- rescaleCovW0(res$mrmash_fitted$w0)
  if (length(w0Updated) == 0) {
    return(list())
  }
  w0Updated <- w0Updated[names(w0Updated) %in% names(dataDrivenPriorMatrices$U)]
  dataDrivenPriorMatrices$U <- dataDrivenPriorMatrices$U[names(w0Updated)]
  dataDrivenPriorMatrices$w <- dataDrivenPriorMatrices$w[names(w0Updated)]

  if (!is.null(dataDrivenPriorMatricesCv)) {
    for (fold in seq_along(dataDrivenPriorMatricesCv)) {
      dataDrivenPriorMatricesCv[[fold]] <- filterMixtureComponents(
        colnames(Y), dataDrivenPriorMatricesCv[[fold]]$U,
        dataDrivenPriorMatricesCv[[fold]]$w, dataDrivenPriorWeightsCutoff
      )
      dataDrivenPriorMatricesCv[[fold]]$w <- dataDrivenPriorMatricesCv[[fold]]$w[names(dataDrivenPriorMatricesCv[[fold]]$w) %in% names(w0Updated)]
      dataDrivenPriorMatricesCv[[fold]]$w <- w0Updated[names(dataDrivenPriorMatricesCv[[fold]]$w)]
      dataDrivenPriorMatricesCv[[fold]]$U <- dataDrivenPriorMatricesCv[[fold]]$U[names(dataDrivenPriorMatricesCv[[fold]]$U) %in% names(w0Updated)]
    }
  } else if (is.null(dataDrivenPriorMatricesCv) && !is.null(dataDrivenPriorMatrices)) {
    dataDrivenPriorMatricesCv <- lapply(1:cvFolds, function(fold) dataDrivenPriorMatrices)
    names(dataDrivenPriorMatricesCv) <- paste0("fold_", 1:cvFolds)
  }

  mvsusieReweightedMixturePrior <- initializeMvsusiePrior(
    colnames(Y), dataDrivenPriorMatrices,
    dataDrivenPriorMatricesCv, cvFolds, w0Updated, dataDrivenPriorWeightsCutoff
  )
  res$reweighted_mixture_prior <- mvsusieReweightedMixturePrior$dataDrivenPriorMatrices
  res$reweighted_mixture_prior_cv <- mvsusieReweightedMixturePrior$dataDrivenPriorMatricesCv

  # Fit mvSuSiE
  message("Fitting mvSuSiE model on input data ...")
  res$mvsusie_fitted <- mvsusieR::mvsusie(X, Y,
    L = L, L_greedy = lGreedy,
    prior_variance = mvsusieReweightedMixturePrior$dataDrivenPriorMatrices,
    residual_variance = residY, estimate_residual_variance = estimateResidualVariance,
    max_iter = mvsusieMaxIter,
    verbose = verbose, coverage = coverage[1]
  )

  # Process mvSuSiE results
  secCoverage <- if (length(coverage) > 1) coverage[-1] else NULL
  mvsusiePost <- postprocessFinemappingFits(
    fits = list(mvsusie = .setFinemappingFitClass(res$mvsusie_fitted, "mvsusie")),
    dataX = X,
    dataY = NULL,
    xScalar = 1,
    yScalar = 1,
    af = af,
    coverage = coverage[1],
    secondaryCoverage = secCoverage,
    signalCutoff = signalCutoff,
    minAbsCorr = minAbsCorr,
    otherQuantities = otherQuantities,
    region = region
  )
  res <- c(res, formatFinemappingOutput(mvsusiePost, primaryMethod = "mvsusie"))
  res$total_time_elapsed <- proc.time() - st

  # Run TWAS weights and optionally CV
  if (twasWeights) {
    res$twas_weights_result <- twasMultivariateWeightsPipeline(X, Y, res,
      cvFolds = cvFolds, samplePartition = samplePartition,
      maxCvVariants = maxCvVariants,
      mvsusieMaxIter = mvsusieMaxIter, mrmashMaxIter = mrmashMaxIter,
      canonicalPriorMatrices = canonicalPriorMatrices,
      dataDrivenPriorMatrices = dataDrivenPriorMatrices,
      dataDrivenPriorMatricesCv = dataDrivenPriorMatricesCv,
      L = L, Lgreedy = lGreedy,
      cvThreads = cvThreads, verbose = verbose
    )
  }
  return(res)
}

