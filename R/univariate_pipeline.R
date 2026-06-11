#' Univariate Analysis Pipeline
#'
#' This function performs univariate analysis for fine-mapping and Transcriptome-Wide Association Study (TWAS)
#' with optional cross-validation. By default, fine-mapping fits SuSiE-inf first
#' and then fits SuSiE initialized from the SuSiE-inf result.
#'
#' @param X A matrix of genotype data where rows represent samples and columns represent genetic variants.
#' @param Y A vector of phenotype measurements.
#' @param X_scalar A scalar or vector to rescale X to its original scale.
#' @param Y_scalar A scalar to rescale Y to its original scale.
#' @param maf Optional vector of minor allele frequencies for each variant in X,
#'   used ONLY for \code{maf_cutoff} filtering and never exported. \code{af} is
#'   the single source of truth: when \code{af} is supplied the filtering MAF is
#'   derived from it (\code{min(af, 1 - af)}) and a supplied \code{maf} is
#'   ignored (with a warning if they disagree). Default NULL; if neither
#'   \code{maf} nor \code{af} is supplied and \code{maf_cutoff} is set, the call
#'   errors.
#' @param af Optional vector of directional effect-allele frequencies (frequency
#'   of \code{a1}) aligned to the columns of X. When supplied it is exported as
#'   the \code{top_loci$af} column; when NULL, \code{af} is \code{NA_real_}.
#'   Default NULL.
#' @param X_variance Optional variance of X. Default is NULL.
#' @param other_quantities A list of other quantities to be carried into fine-mapping post-processing. Default is an empty list.
#' @param region Optional \code{"chr:start-end"} string for the analysis region. Default is NULL.
#' @param imiss_cutoff Individual missingness cutoff. Default is 1.0.
#' @param maf_cutoff Minor allele frequency cutoff. Default is NULL.
#' @param xvar_cutoff Variance cutoff for X. Default is 0.05.
#' @param ld_reference_meta_file An optional path to a file containing linkage disequilibrium reference data. Default is NULL.
#' @param pip_cutoff_to_skip Cutoff value for skipping analysis based on PIP values. Default is 0.
#' @param L Maximum number of components in SuSiE. Default is 20.
#' @param L_greedy Initial greedy number of components in SuSiE. Default is 5.
#' @param signal_cutoff Cutoff value for signal identification in PIP values. Default is 0.025.
#' @param coverage A vector of coverage probabilities for credible sets. Default is c(0.95, 0.7, 0.5).
#' @param min_abs_corr Minimum absolute correlation for credible set purity filtering. Default is 0.8,
#'   which is stricter than the susieR default of 0.5.
#' @param finemapping_extra_opts Additional options passed to \code{susieR::susie()}.
#'   SuSiE-inf is always fitted with \code{refine = FALSE}; the ordinary SuSiE
#'   fit keeps these options and is initialized with \code{model_init}.
#' @param estimate_residual_variance Passed to \code{susieR::susie()}. Default is TRUE.
#' @param methods Optional character vector selecting which SuSiE variants to
#'   fit. Any subset of \code{c("susie", "susie_inf", "susie_ash")}. Default
#'   \code{NULL} falls back to the legacy \code{add_susie_inf} behavior:
#'   \code{add_susie_inf = TRUE} (default) maps to
#'   \code{methods = c("susie_inf", "susie")} with SuSiE-inf chained into the
#'   SuSiE fit as initialization; \code{add_susie_inf = FALSE} maps to
#'   \code{methods = "susie"} (plain SuSiE alone). When \code{methods} is
#'   passed explicitly, each requested method is fitted; if
#'   \code{"susie_inf"} is paired with \code{"susie"} or \code{"susie_ash"}
#'   (or both) and \code{add_susie_inf = TRUE}, the SuSiE-inf fit
#'   initialises each chained downstream method. This gives five distinct
#'   fitting modes: SuSiE alone, SuSiE with SuSiE-inf init, SuSiE-inf alone,
#'   SuSiE-ash alone, and SuSiE-ash with SuSiE-inf init.
#' @param add_susie_inf When \code{methods} is \code{NULL}, controls whether
#'   SuSiE-inf is fitted and chained into SuSiE. When \code{methods} is set
#'   explicitly, controls whether the chained-init shortcut is applied to
#'   any \code{"susie"} or \code{"susie_ash"} method present alongside
#'   \code{"susie_inf"}. Default \code{TRUE}.
#' @param twas_weights Whether to compute TWAS weights. Default is TRUE.
#' @param sample_partition Optional data frame with Sample and Fold columns for cross-validation. Default is NULL.
#' @param max_cv_variants The maximum number of variants to be included in cross-validation. Default is -1 (no limit).
#' @param cv_folds The number of folds to use for cross-validation. Default is 5.
#' @param cv_threads The number of threads to use for parallel computation in cross-validation. Default is 1.
#' @param verbose Verbosity level. Default is 0.
#'
#' @return A list containing the univariate analysis results.
#' @importFrom susieR susie
#' @export
univariateAnalysisPipeline <- function(
    # input data
    X,
    Y,
    maf = NULL,
    af = NULL,
    xScalar = 1,
    yScalar = 1,
    xVariance = NULL,
    otherQuantities = list(),
    region = NULL,
    # filters
    imissCutoff = 1.0,
    mafCutoff = NULL,
    xvarCutoff = 0,
    ldReferenceMetaFile = NULL,
    pipCutoffToSkip = 0,
    # methods parameter configuration
    L = 20,
    lGreedy = 5,
    # fine-mapping results summary
    signalCutoff = 0.025,
    coverage = c(0.95, 0.7, 0.5),
    minAbsCorr = 0.8,
    finemappingExtraOpts = list(refine = TRUE),
    estimateResidualVariance = TRUE,
    methods = NULL,
    addSusieInf = TRUE,
    # TWAS weights and CV for TWAS weights
    twasWeights = TRUE,
    samplePartition = NULL,
    maxCvVariants = -1,
    cvFolds = 5,
    cvThreads = 1,
    verbose = 0) {
  # Input validation
  if (!is.matrix(X) || !is.numeric(X)) stop("X must be a numeric matrix")
  if (!is.vector(Y) && !(is.matrix(Y) && ncol(Y) == 1) || !is.numeric(Y)) stop("Y must be a numeric vector or a single column matrix")
  if (nrow(X) != length(Y)) stop("X and Y must have the same number of rows/length")
  # maf is optional (directionless, used ONLY for mafCutoff filtering, never
  # exported). af (directional effect-allele frequency) is optional and is the
  # single source of truth: when supplied it is exported as top_loci$af and the
  # filtering MAF is derived from it.
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
  # If a MAF cutoff is requested, a frequency must be available to derive it.
  if (is.null(maf) && !is.null(mafCutoff) && is.numeric(mafCutoff) && mafCutoff > 0) {
    stop("mafCutoff is set but neither 'af' nor 'maf' was supplied; provide ",
         "one so MAF can be derived for filtering.")
  }
  if (!is.numeric(xScalar) || (length(xScalar) != 1 && length(xScalar) != ncol(X))) stop("xScalar must be a numeric scalar or vector with length equal to the number of columns in X")
  if (!is.numeric(yScalar) || length(yScalar) != 1) stop("yScalar must be a numeric scalar")
  if (!is.numeric(L) || L <= 0) stop("L must be a positive integer")
  if (!is.null(lGreedy) && (!is.numeric(lGreedy) || lGreedy <= 0)) stop("lGreedy must be NULL or a positive integer")
  if (!is.logical(addSusieInf) || length(addSusieInf) != 1 || is.na(addSusieInf)) {
    stop("addSusieInf must be TRUE or FALSE")
  }

  # Resolve effective methods. NULL => backward-compat via addSusieInf.
  validMethods <- c("susie", "susie_inf", "susie_ash")
  if (is.null(methods)) {
    methods <- if (isTRUE(addSusieInf)) c("susie_inf", "susie") else "susie"
  } else {
    if (!is.character(methods) || length(methods) == 0L) {
      stop("methods must be a non-empty character vector of method names.")
    }
    bad <- setdiff(methods, validMethods)
    if (length(bad) > 0) {
      stop("Unknown method(s): ", paste(bad, collapse = ", "),
           ". Valid options: ", paste(validMethods, collapse = ", "))
    }
    methods <- unique(methods)
  }
  # SuSiE-inf initialisation chains into SuSiE and/or SuSiE-ash whenever
  # either of them is requested alongside SuSiE-inf and addSusieInf is TRUE.
  chainInfToSusie     <- isTRUE(addSusieInf) &&
    all(c("susie_inf", "susie") %in% methods)
  chainInfToSusieAsh  <- isTRUE(addSusieInf) &&
    all(c("susie_inf", "susie_ash") %in% methods)
  anyChainedInit <- chainInfToSusie || chainInfToSusieAsh
  if (isTRUE(twasWeights) && !("susie" %in% methods)) {
    stop("twasWeights = TRUE requires \"susie\" to be in methods.")
  }
  if (isTRUE(twasWeights) && !chainInfToSusie) {
    stop("twasWeights = TRUE requires SuSiE to be initialised from SuSiE-inf; ",
         "set methods = c(\"susie_inf\", \"susie\") and addSusieInf = TRUE.")
  }

  # Initial PIP check
  if (pipCutoffToSkip != 0) {
    if (pipCutoffToSkip < 0) {
      # automatically determine the cutoff to use
      pipCutoffToSkip <- 3 * 1 / ncol(X)
    }
    topModelPip <- susie(X, Y, L = 1)$pip
    if (!any(topModelPip > pipCutoffToSkip)) {
      message(paste("Skipping follow-up analysis: No signals above PIP threshold", pipCutoffToSkip, "in initial model screening."))
      return(list())
    } else {
      message(paste("Follow-up on region because signals above PIP threshold", pipCutoffToSkip, "were detected in initial model screening."))
    }
  }

  # Filter variants if LD reference is provided
  if (!is.null(ldReferenceMetaFile)) {
    variantsKept <- filterVariantsByLdReference(colnames(X), ldReferenceMetaFile)
    X <- X[, variantsKept$data, drop = FALSE]
    if (!is.null(maf)) maf <- maf[variantsKept$idx]
    if (!is.null(af)) af <- af[variantsKept$idx]
    if (length(xScalar) > 1) xScalar <- xScalar[variantsKept$idx]
  }

  # Filter X based on missingness, MAF, and variance
  if (!is.null(imissCutoff) || !is.null(mafCutoff)) {
    XFiltered <- filterX(X, imissCutoff, mafCutoff, varThresh = xvarCutoff, maf = maf, xVariance = xVariance)
    keptIndices <- match(colnames(XFiltered), colnames(X))
    if (!is.null(maf)) maf <- maf[keptIndices]
    if (!is.null(af)) af <- af[keptIndices]
    if (length(xScalar) > 1) xScalar <- xScalar[keptIndices]
    X <- XFiltered
  }

  # Main analysis
  st <- proc.time()
  res <- list()

  susieArgs <- modifyList(
    finemappingExtraOpts,
    list(L = L, L_greedy = lGreedy, coverage = coverage[1],
         estimate_residual_variance = estimateResidualVariance)
  )
  fittedModels <- list()

  if ("susie_inf" %in% methods || anyChainedInit) {
    message("Fitting SuSiE-inf model on input data ...")
    infArgs <- modifyList(susieArgs, list(
      X = X, y = Y,
      unmappable_effects = "inf",
      convergence_method = "pip",
      refine = FALSE, model_init = NULL
    ))
    infFit <- do.call(susie, infArgs)
    fittedModels[["susie_inf"]] <- .setFinemappingFitClass(infFit, "susie_inf")
  }

  if ("susie" %in% methods) {
    if (chainInfToSusie) {
      message("Fitting SuSiE model initialized by SuSiE-inf ...")
      suArgs <- prepareSusieFromInfArgs(susieArgs,
                                        fittedModels[["susie_inf"]],
                                        refineDefault = TRUE,
                                        unmappableEffects = "none")
      suFit <- do.call(susie, c(list(X = X, y = Y), suArgs))
    } else {
      message("Fitting SuSiE model on input data ...")
      suFit <- do.call(susie, c(list(X = X, y = Y), susieArgs))
    }
    fittedModels[["susie"]] <- .setFinemappingFitClass(suFit, "susie")
  }

  if ("susie_ash" %in% methods) {
    if (chainInfToSusieAsh) {
      message("Fitting SuSiE-ash model initialized by SuSiE-inf ...")
      ashArgs <- prepareSusieFromInfArgs(susieArgs,
                                         fittedModels[["susie_inf"]],
                                         refineDefault = NULL,
                                         unmappableEffects = "ash")
      ashFit <- do.call(susie, c(list(X = X, y = Y), ashArgs))
    } else {
      message("Fitting SuSiE-ash model on input data ...")
      ashArgs <- modifyList(susieArgs, list(
        X = X, y = Y,
        unmappable_effects = "ash",
        convergence_method = "pip"
      ))
      ashFit <- do.call(susie, ashArgs)
    }
    fittedModels[["susie_ash"]] <- .setFinemappingFitClass(ashFit, "susie_ash")
  }

  # Drop susie_inf from post-processing if it was only fit to provide init for
  # SuSiE / SuSiE-ash (i.e. caller did not request "susie_inf" in methods).
  if (anyChainedInit && !("susie_inf" %in% methods)) {
    fittedModels[["susie_inf"]] <- NULL
  }

  # Back-compat slots for the most common methods
  res$susie_inf_fitted <- fittedModels[["susie_inf"]]
  res$susie_fitted     <- fittedModels[["susie"]]
  res$susie_ash_fitted <- fittedModels[["susie_ash"]]

  # Process SuSiE results
  susiePost <- postprocessFinemappingFits(
    fits = fittedModels,
    dataX = X,
    dataY = Y,
    xScalar = xScalar,
    yScalar = yScalar,
    af = af,
    coverage = coverage[1],
    secondaryCoverage = if (length(coverage) > 1) coverage[-1] else NULL,
    signalCutoff = signalCutoff,
    minAbsCorr = minAbsCorr,
    otherQuantities = otherQuantities,
    region = region
  )
  # Primary method drives root-level finemapping_result / sumstats / etc.
  # Preference order favors "susie" for backward compatibility, then
  # falls back to the first requested method actually fitted.
  primaryMethod <- if ("susie" %in% names(fittedModels)) "susie" else names(fittedModels)[1]
  res <- c(res, formatFinemappingOutput(susiePost, primaryMethod = primaryMethod))
  susieInfFm <- susiePost$finemapping_results$susie_inf$finemapping_result
  res$susie_inf_result_trimmed <- if (!is.null(susieInfFm)) getTrimmedFit(susieInfFm) else NULL
  susieAshFm <- susiePost$finemapping_results$susie_ash$finemapping_result
  res$susie_ash_result_trimmed <- if (!is.null(susieAshFm)) getTrimmedFit(susieAshFm) else NULL
  res$total_time_elapsed <- proc.time() - st

  # TWAS weights and cross-validation
  if (twasWeights) {
    res$twas_weights_result <- twasWeightsPipeline(
      X, Y, fittedModels = fittedModels,
      cvFolds = cvFolds,
      maxCvVariants = maxCvVariants,
      cvThreads = cvThreads,
      samplePartition = samplePartition
    )
    if ("top_loci" %in% names(res) && !is.null(res$twas_weights_result$susie_weights_intermediate)) {
      res$twas_weights_result$susie_weights_intermediate$top_loci <- res$top_loci
    }
  }

  return(res)
}

#' Load LD for a study, supporting single or mixture panels.
#'
#' @param ldPath A single LD metadata TSV path, or comma-separated paths for
#'   mixture panels (e.g., "ld_EUR.tsv,ld_AFR.tsv").
#' @param region Region string "chr:start-end".
#' @return An \code{LdData} S4 object. For single panels, returns the result of
#'   \code{loadLdMatrix()} unchanged. For mixture panels, \code{genotypeHandle}
#'   is a list of per-panel genotype handles sharing the first panel's variants.
#' @export
loadStudyLd <- function(ldPath, region) {
  paths <- strsplit(ldPath, ",")[[1]]
  if (length(paths) == 1) {
    return(loadLdMatrix(paths, region, returnGenotype = "auto"))
  }
  # Mixture: load each panel; combine handles into a list
  base <- loadLdMatrix(paths[1], region, returnGenotype = TRUE)
  otherHandles <- lapply(paths[-1], function(p) {
    ld <- loadLdMatrix(p, region, returnGenotype = TRUE)
    ld@genotypeHandle
  })
  allHandles <- c(list(base@genotypeHandle), otherHandles)
  LdData(
    correlation = NULL,
    genotypeHandle = allHandles,
    snpIdx = base@snpIdx,
    variants = base@variants,
    blockMetadata = base@blockMetadata,
    nRef = base@nRef
  )
}

.rss_variant_ids <- function(sumstats) {
  if ("variant_id" %in% names(sumstats)) return(as.character(sumstats$variant_id))
  if ("variant" %in% names(sumstats)) return(as.character(sumstats$variant))
  rn <- rownames(sumstats)
  if (!is.null(rn) && length(rn) == nrow(sumstats) &&
      !all(grepl("^[0-9]+$", rn))) {
    return(rn)
  }
  stop("RSS sumstats must contain a variant_id or variant column.")
}

.rss_sumstats_with_variant_id <- function(sumstats) {
  if (!"variant_id" %in% names(sumstats)) {
    sumstats$variant_id <- .rss_variant_ids(sumstats)
  }
  sumstats$variant_id <- as.character(sumstats$variant_id)
  sumstats
}

.match_rss_variants <- function(variants, referenceIds, referenceName) {
  idx <- match(variants, referenceIds)
  if (anyNA(idx)) {
    idx <- match(stripChrPrefix(stripBuildSuffix(variants)),
                 stripChrPrefix(stripBuildSuffix(referenceIds)))
  }
  if (anyNA(idx)) {
    missingVar <- variants[is.na(idx)]
    stop(referenceName, " is missing ", length(missingVar),
         " variant(s): ", paste(utils::head(missingVar, 3), collapse = ", "))
  }
  idx
}

.subset_rss_matrix_columns <- function(X, variants, referenceIds, referenceName) {
  if (is.null(colnames(X)) && length(referenceIds) == ncol(X)) {
    colnames(X) <- referenceIds
  }
  idx <- .match_rss_variants(variants, colnames(X), referenceName)
  Xout <- X[, idx, drop = FALSE]
  colnames(Xout) <- variants
  Xout
}

.subset_rss_ld_matrix <- function(R, variants, referenceIds) {
  if (is.null(rownames(R)) && length(referenceIds) == nrow(R)) {
    rownames(R) <- referenceIds
  }
  if (is.null(colnames(R)) && length(referenceIds) == ncol(R)) {
    colnames(R) <- referenceIds
  }
  idx <- .match_rss_variants(variants, rownames(R), "LD matrix")
  Rout <- R[idx, idx, drop = FALSE]
  rownames(Rout) <- colnames(Rout) <- variants
  Rout
}

#' Convert one loaded RSS record to direct SuSiE RSS input
#'
#' @param rssInput A single loaded RSS record, usually one element of
#'   \code{qced_regional_data$sumstat_data$sumstats}. It must contain
#'   \code{sumstats}, \code{n}, and \code{var_y}.
#' @param ldData A matching \code{LdData} object for the same study.
#' @return A list with \code{susie_rss_input}, ready to pass to
#'   \code{\link{susieRssPipeline}}, and \code{source_info}.
#' @export
regionDataToSusieRssInput <- function(rssInput, ldData) {
  if (!is.list(rssInput) || is.null(rssInput$sumstats)) {
    stop("rssInput must be a single RSS record with a sumstats element.")
  }
  if (is.null(ldData) || !is(ldData, "LdData")) {
    stop("ldData must be an LdData object.")
  }

  sumstats <- .rss_sumstats_with_variant_id(rssInput$sumstats)
  variants <- sumstats$variant_id
  if (length(variants) == 0L) {
    stop("rssInput$sumstats contains no variants.")
  }

  referenceIds <- getVariantIds(ldData)
  if (hasGenotypes(ldData)) {
    X <- getGenotypes(ldData)
    xMat <- if (is.list(X) && !is.matrix(X)) {
      lapply(X, .subset_rss_matrix_columns, variants = variants,
             referenceIds = referenceIds,
             referenceName = "genotype reference panel")
    } else {
      .subset_rss_matrix_columns(X, variants, referenceIds,
                                 referenceName = "genotype reference panel")
    }
    ldMat <- NULL
  } else {
    R <- ldData@correlation
    if (is.null(R) || (is.list(R) && !is.matrix(R))) {
      stop("ldData must contain one correlation matrix or genotype data.")
    }
    ldMat <- .subset_rss_ld_matrix(R, variants, referenceIds)
    xMat <- NULL
  }

  varY <- rssInput$varY %||% rssInput$var_y

  list(
    susie_rss_input = list(
      sumstats = sumstats,
      ldMat = ldMat,
      xMat = xMat,
      n = rssInput$n,
      var_y = varY
    ),
    source_info = list(
      n_variants = length(variants),
      variants = variants,
      uses_X_ref = !is.null(xMat),
      has_LD = !is.null(ldMat)
    )
  )
}

#' RSS Analysis Pipeline
#'
#' End-to-end pipeline for summary statistics fine-mapping via SuSiE RSS.
#' Supports both z+R (correlation matrix) and z+X (genotype matrix) interfaces.
#'
#' @param sumstatPath File path to the summary statistics.
#' @param columnFilePath File path to the column mapping file.
#' @param ldData An \code{LdData} S4 object from \code{loadLdMatrix()}. When
#'   \code{hasGenotypes(ldData)} is TRUE (from \code{returnGenotype=TRUE}),
#'   susie_rss uses the z+X interface via \code{getGenotypes()}. Local R is
#'   computed only for QC stages that require a correlation matrix.
#' @param nSample Sample size. If 0, retrieved from the sumstat file.
#' @param nCase Number of cases (for case-control studies).
#' @param nControl Number of controls (for case-control studies).
#' @param binaryTraitModel How to handle case-control summary statistics.
#'   The default \code{"rss"} uses the z-score RSS interface and does not pass
#'   a phenotype variance to \code{susieR::susie_rss()}. Use \code{"ols"} only
#'   when \code{beta} and \code{se} are from OLS on a centered 0/1 phenotype;
#'   then \code{varY} is computed from \code{nCase/n} and passed through to
#'   select the \code{bhat/shat/var_y} sufficient-statistic interface.
#' @param region Region string "chr:start-end" for tabix subsetting.
#' @param skipRegion Character vector of regions to skip (format "chrom:start-end").
#' @param extractRegionName Gene/phenotype name to subset.
#' @param regionNameCol Column to filter for extractRegionName.
#' @param qcMethod Summary-statistic QC method. \code{"slalom"} and
#'   \code{"dentist"} run basic allele harmonization plus LD-mismatch QC;
#'   \code{"none"} runs basic allele harmonization only.
#' @param finemappingMethod Iteration mode for the SuSiE-RSS fit (when
#'   \code{"susie_rss"} is among \code{methods}). One of \code{"susie_rss"}
#'   (default normal IBSS), \code{"single_effect"} (L=1, single iteration),
#'   or \code{"bayesian_conditional_regression"} (full L, single iteration).
#' @param methods Optional character vector selecting which SuSiE-RSS
#'   variants to fit. Any subset of \code{c("susie_rss", "susie_inf_rss",
#'   "susie_ash_rss")}. Default \code{NULL} preserves legacy single-method
#'   behavior via \code{finemappingMethod}. When set explicitly, every
#'   requested method contributes rows to the unified \code{top_loci}; when
#'   \code{"susie_inf_rss"} is paired with \code{"susie_rss"} or
#'   \code{"susie_ash_rss"} (or both) and \code{addSusieInf = TRUE}, the
#'   SuSiE-inf-RSS fit initialises the chained downstream method(s).
#' @param addSusieInf Logical controlling chained init when
#'   \code{"susie_inf_rss"} is in \code{methods} alongside
#'   \code{"susie_rss"} and/or \code{"susie_ash_rss"}. Default \code{TRUE}.
#' @param finemappingOpts List of fine-mapping options (L, L_greedy, coverage,
#'   signal_cutoff, min_abs_corr).
#' @param impute Whether to impute missing variants via RAISS (default TRUE).
#' @param imputeOpts List of imputation options (rcond, R2_threshold, minimum_ld, lamb).
#' @param pipCutoffToSkip PIP threshold for early stopping (default 0, no skip).
#' @param rFinite Controls variance inflation to account for finite reference LD.
#'   Passed to \code{susieR::susie_rss()}.
#' @param rMismatch LD mismatch correction method passed to \code{susieR::susie_rss()}.
#'   Default NULL disables mismatch correction.
#' @param keepIndel Whether to keep indel variants (default TRUE).
#' @param commentString Comment character for sumstat file (default "#").
#' @param diagnostics Whether to include diagnostic info (default FALSE).
#'
#' @return A list with fine-mapping results and analyzed summary statistics.
#' @importFrom magrittr %>%
#' @importFrom susieR susie_rss
#' @export
rssAnalysisPipeline <- function(
    sumstatPath, columnFilePath, ldData,
    nSample = 0, nCase = 0, nControl = 0, region = NULL, skipRegion = NULL,
    extractRegionName = NULL, regionNameCol = NULL,
    qcMethod = c("slalom", "dentist", "none"),
    finemappingMethod = c("susie_rss", "single_effect", "bayesian_conditional_regression"),
    methods = NULL,
    addSusieInf = TRUE,
    finemappingOpts = list(
      L = 20, L_greedy = 5,
      coverage = c(0.95, 0.7, 0.5), signal_cutoff = 0.025,
      min_abs_corr = 0.8
    ),
    impute = TRUE, imputeOpts = list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5, lamb = 0.01),
    pipCutoffToSkip = 0, rFinite = NULL, rMismatch = NULL,
    keepIndel = TRUE, commentString = "#", diagnostics = FALSE,
    binaryTraitModel = c("rss", "ols")) {
  binaryTraitModel <- match.arg(binaryTraitModel)
  if (!is(ldData, "LdData")) {
    stop("ldData must be an LdData object")
  }
  res <- list()
  rssInput <- loadRssData(
    sumstatPath = sumstatPath, columnFilePath = columnFilePath,
    nSample = nSample, nCase = nCase, nControl = nControl,
    extractRegionName = extractRegionName, region = region,
    regionNameCol = regionNameCol, commentString = commentString,
    binaryTraitModel = binaryTraitModel
  )

  sumstats <- rssInput$sumstats
  n <- rssInput$n
  varY <- rssInput$var_y

  if (nrow(sumstats) == 0) {
    return(list(rss_data_analyzed = sumstats))
  }

  qcMethodArg <- if (is.null(qcMethod)) NULL else match.arg(qcMethod)
  qcMethod <- qcMethodArg
  qcRecord <- summaryStatsQc(
    rssInput = rssInput,
    ldData = ldData,
    keepIndel = keepIndel,
    skipRegion = skipRegion,
    pipCutoffToSkip = pipCutoffToSkip,
    qcMethod = if (is.null(qcMethodArg)) "none" else qcMethodArg,
    impute = impute,
    imputeOpts = imputeOpts,
    returnOnSkip = "preprocess",
    R_finite = rFinite,
    R_mismatch = rMismatch
  )
  if (!is(qcRecord, "QcResult")) {
    stop("summaryStatsQc must return a QcResult object.")
  }
  rssRecord <- getRssInput(qcRecord)
  sumstats <- rssRecord$sumstats
  n <- rssRecord$n
  varY <- rssRecord$var_y
  preprocessSnapshot <- getPreprocess(qcRecord)
  preprocessLd <- preprocessSnapshot$ldData %||% preprocessSnapshot$ld_data
  preprocessResults <- list(
    sumstats = preprocessSnapshot$sumstats
  )
  qcResults <- list(outlier_number = getOutlierNumber(qcRecord))

  if (nrow(sumstats) == 0) {
    message("No variants left after preprocessing. Returning empty results.")
    return(list(rss_data_analyzed = sumstats))
  }
  if (isSkipped(qcRecord)) {
    return(list(rss_data_analyzed = sumstats))
  }

  qcLd <- getLdData(qcRecord)
  susieReady <- regionDataToSusieRssInput(rssRecord, qcLd)$susie_rss_input

  # Fine-mapping: use xMat if available, otherwise R
  if (!is.null(finemappingMethod)) {
    priCoverage <- finemappingOpts$coverage[1]
    secCoverage <- if (length(finemappingOpts$coverage) > 1) finemappingOpts$coverage[-1] else NULL

    finemappingOptsL_greedy <- finemappingOpts$L_greedy %||% finemappingOpts$lGreedy
    finemappingOptsSignalCutoff <- finemappingOpts$signal_cutoff %||% finemappingOpts$signalCutoff
    finemappingOptsMinAbsCorr <- finemappingOpts$min_abs_corr %||% finemappingOpts$minAbsCorr
    res <- do.call(susieRssPipeline, c(susieReady, list(
      L = finemappingOpts$L, lGreedy = finemappingOptsL_greedy,
      analysisMethod = finemappingMethod,
      methods = methods,
      addSusieInf = addSusieInf,
      coverage = priCoverage,
      secondaryCoverage = secCoverage,
      signalCutoff = finemappingOptsSignalCutoff,
      minAbsCorr = finemappingOptsMinAbsCorr,
      rFinite = rFinite,
      rMismatch = rMismatch
    )))
    if (!is.null(qcMethod)) {
      res$outlier_number <- qcResults$outlier_number
    }
  }
  .makeMethodName <- function(method, qcMethod, impute) {
    suffix <- if (!is.null(qcMethod) && impute) {
      paste0(toupper(qcMethod), "_RAISS_imputed")
    } else if (!is.null(qcMethod)) {
      toupper(qcMethod)
    } else {
      "NO_QC"
    }
    paste0(method, "_", suffix)
  }

  .runReanalysis <- function(sumstats, ldData, method, finemappingOpts, priCoverage, secCoverage) {
    reanalysisInput <- regionDataToSusieRssInput(
      list(sumstats = sumstats, n = n, var_y = varY),
      ldData
    )$susie_rss_input
    do.call(susieRssPipeline, c(reanalysisInput, list(
      L = finemappingOpts$L, lGreedy = finemappingOpts$L_greedy %||% finemappingOpts$lGreedy,
      analysisMethod = method,
      coverage = priCoverage,
      secondaryCoverage = secCoverage,
      signalCutoff = finemappingOptsSignalCutoff,
      minAbsCorr = finemappingOptsMinAbsCorr,
      rFinite = rFinite,
      rMismatch = rMismatch
    )))
  }

  methodName <- .makeMethodName(finemappingMethod, qcMethod, impute)
  resultList <- list()
  resultList[[methodName]] <- res
  resultList[["rss_data_analyzed"]] <- sumstats

  blockCsMetrics <- list()
  if (diagnostics) {
    if (length(res) > 0) {
        bvsrRes <- getSusieResult(res)
        bvsrCsNum <- if(!is.null(bvsrRes)) length(bvsrRes$sets$cs) else NULL
        if (isTRUE(bvsrCsNum > 0)) { # have CS
            csNamesBvsr <- names(bvsrRes$sets$cs)
            blockCsMetrics <- extractCsInfo(conData = res, csNames = csNamesBvsr, topLociTable = res$top_loci)
        } else { # no CS
            if (sum(bvsrRes$pip > finemappingOpts$signal_cutoff) > 0) {
                blockCsMetrics <- extractTopPipInfo(res)
            }
        }
    }
    # sensitive check for additional analyses
    if (!is.null(blockCsMetrics) && length(blockCsMetrics) > 0) {
      blockCsMetrics <- parseCsCorr(blockCsMetrics)
      csRow <- blockCsMetrics %>% filter(!is.na(blockCsMetrics$variants_per_cs))
      if (nrow(csRow) > 1) {# CS > 1
        blockCsMetrics <- blockCsMetrics %>%
          mutate(max_cs_corr_study_block = if(all(is.na(cs_corr_max))) {
            NA_real_
          } else {
            max(cs_corr_max, na.rm = TRUE)
          })
        if (any(blockCsMetrics$p_value > 1e-4 | blockCsMetrics$max_cs_corr_study_block > 0.5)) {
          bcr <- .runReanalysis(sumstats, qcLd, "bayesian_conditional_regression",
            finemappingOpts, priCoverage, secCoverage)
          if (!is.null(qcMethod)) {
            bcr$outlier_number <- qcResults$outlier_number
          }
          resultList[[.makeMethodName("bayesian_conditional_regression", qcMethod, impute)]] <- bcr
          ser <- .runReanalysis(preprocessResults$sumstats, preprocessLd,
            "single_effect", finemappingOpts, priCoverage, secCoverage)
          resultList[["single_effect_NO_QC"]] <- ser
        }
      } else { # CS = 1 or NA
        ser <- .runReanalysis(preprocessResults$sumstats, preprocessLd,
          "single_effect", finemappingOpts, priCoverage, secCoverage)
        resultList[["single_effect_NO_QC"]] <- ser
      }
    resultList[["diagnostics"]] <- blockCsMetrics
    }
  }
  return(resultList)
}

