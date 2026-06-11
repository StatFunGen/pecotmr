#' Convert Log Bayes Factors to Single Effects PIP
#'
#' This function converts log Bayes factors (LBF) to alpha values, optionally
#' using prior weights. It handles numerical stability by adjusting with the
#' maximum LBF value.
#'
#' @param lbf Numeric vector of log Bayes factors.
#' @param priorWeights Optional numeric vector of prior weights for each element in lbf.
#' @return A named numeric vector of alpha values corresponding to the input LBF.
#' @examples
#' lbf <- c(-0.5, 1.2, 0.3)
#' alpha <- lbfToAlphaVector(lbf)
#' print(alpha)
#' @noRd
lbfToAlphaVector <- function(lbf, priorWeights = NULL) {
  if (is.null(priorWeights)) priorWeights <- rep(1 / length(lbf), length(lbf))
  maxlbf <- max(lbf)

  # If maxlbf is 0, return a vector of zeros
  if (maxlbf == 0) {
    return(setNames(rep(0, length(lbf)), names(lbf)))
  }

  # w is proportional to BF, subtract max for numerical stability
  w <- exp(lbf - maxlbf)

  # Posterior prob for each SNP
  wWeighted <- w * priorWeights
  weightedSumW <- sum(wWeighted)
  alpha <- wWeighted / weightedSumW

  return(alpha)
}

#' Applies the 'lbfToAlphaVector' function row-wise to a matrix of log Bayes factors
#' to convert them to Single Effect PIP values.
#'
#' @param lbf Matrix of log Bayes factors.
#' @return A matrix of alpha values with the same dimensions as the input LBF matrix.
#' @examples
#' lbfMatrix <- matrix(c(-0.5, 1.2, 0.3, 0.7, -1.1, 0.4), nrow = 2)
#' alphaMatrix <- lbfToAlpha(lbfMatrix)
#' print(alphaMatrix)
#' @export
lbfToAlpha <- function(lbf) {
  alphaMatrix <- t(apply(as.matrix(lbf), 1, lbfToAlphaVector))
  if (ncol(lbf) == 1) alphaMatrix <- matrix(alphaMatrix, ncol = 1, dimnames = list(NULL, colnames(lbf)))
  return(alphaMatrix)
}

formatPipColumn <- function(method) {
  paste0("pip_", method)
}

resolvePipColumn <- function(topLoci, method = NULL) {
  if (is.null(topLoci) || nrow(topLoci) == 0) return(NULL)
  if (!is.null(method)) {
    pipCol <- formatPipColumn(method)
    if (pipCol %in% names(topLoci)) return(pipCol)
  }
  if ("pip" %in% names(topLoci)) return("pip")
  pipCols <- grep("^pip_", names(topLoci), value = TRUE)
  if (length(pipCols) == 1) return(pipCols)
  NULL
}

formatCsColumn <- function(coverage, method) {
  pct <- as.numeric(coverage) * 100
  if (is.na(pct)) stop("coverage must be numeric.")
  label <- if (abs(pct - round(pct)) < 1e-8) {
    as.character(as.integer(round(pct)))
  } else {
    gsub("\\.", "_", format(pct, scientific = FALSE, trim = TRUE))
  }
  paste0("CS_", label, "_", method)
}

.translateLegacyCsColumnName <- function(coverage) {
  if (is.null(coverage)) return(NULL)
  vapply(coverage, function(x) {
    x <- as.character(x)
    oldMatch <- regexec("^cs_coverage_([0-9.]+)$", x, ignore.case = TRUE)
    oldParts <- regmatches(x, oldMatch)[[1]]
    if (length(oldParts) == 2) return(formatCsColumn(as.numeric(oldParts[[2]]), "susie"))
    x
  }, character(1), USE.NAMES = FALSE)
}

.translateLegacyTopLociCsColumns <- function(topLoci) {
  if (!is.data.frame(topLoci)) return(topLoci)
  names(topLoci) <- .translateLegacyCsColumnName(names(topLoci))
  if ("pip_susie" %in% names(topLoci) && !"pip" %in% names(topLoci)) {
    names(topLoci)[names(topLoci) == "pip_susie"] <- "pip"
  }
  topLoci
}

.setFinemappingFitClass <- function(fit, method) {
  if (is.null(fit)) return(NULL)
  methodClass <- switch(method,
    susie = "susie",
    susie_inf = "susie_inf",
    susie_rss = "susie_rss",
    single_effect = "susie_rss",
    bayesian_conditional_regression = "susie_rss",
    fsusie = "susiF",
    mvsusie = "mvsusie",
    NULL
  )
  if (!is.null(methodClass)) class(fit) <- unique(c(methodClass, class(fit)))
  fit
}

# Build the argument list for a SuSiE / SuSiE-ash fit initialised from a
# prior SuSiE-inf fit. `unmappableEffects` controls which branch the
# downstream fit takes: "none" yields the standard SuSiE-inf-initialised
# SuSiE; "ash" yields SuSiE-ash with the SuSiE-inf warm start.
prepareSusieFromInfArgs <- function(args, susieInfFit, refineDefault = NULL,
                                    unmappableEffects = c("none", "ash")) {
  unmappableEffects <- match.arg(unmappableEffects)
  L <- args[["L"]]
  if (is.null(L)) L <- length(susieInfFit$V)
  if (is.null(args[["refine"]]) && !is.null(refineDefault)) args[["refine"]] <- refineDefault
  args[["unmappable_effects"]] <- unmappableEffects
  args[["model_init"]] <- susieInfFit
  if (unmappableEffects == "ash") {
    args[["convergence_method"]] <- args[["convergence_method"]] %||% "pip"
  }
  if (!is.null(args[["L_greedy"]])) args[["L_greedy"]] <- min(length(susieInfFit$V), L)
  args
}

fitSusieInfThenSusie <- function(X, y, args = list(),
                                 susieInfArgs = list(),
                                 susieArgs = list(),
                                 fittedModels = NULL) {
  if (is.null(fittedModels)) fittedModels <- list()
  susieInfFit <- fittedModels[["susie_inf"]]
  susieFit <- fittedModels[["susie"]]

  if (is.null(susieInfFit)) {
    fitArgs <- modifyList(args, susieInfArgs)
    fitArgs <- modifyList(fitArgs, list(
      X = X, y = y, unmappable_effects = "inf",
      convergence_method = "pip", refine = FALSE, model_init = NULL
    ))
    susieInfFit <- do.call(susie, fitArgs)
  }
  susieInfFit <- .setFinemappingFitClass(susieInfFit, "susie_inf")

  if (is.null(susieFit)) {
    fitArgs <- prepareSusieFromInfArgs(modifyList(args, susieArgs), susieInfFit, refineDefault = TRUE)
    susieFit <- do.call(susie, c(list(X = X, y = y), fitArgs))
  }
  susieFit <- .setFinemappingFitClass(susieFit, "susie")

  list(susie = susieFit, susie_inf = susieInfFit)
}

#' Two-stage SuSiE-RSS Fine-mapping
#'
#' RSS analog of \code{fitSusieInfThenSusie}. Fits SuSiE-inf via
#' \code{susie_rss} first, then initialises standard SuSiE-RSS from
#' the SuSiE-inf result. The single pair of fits can be used both for
#' fine-mapping post-processing and TWAS weight extraction.
#'
#' @param z Numeric vector of z-scores.
#' @param R LD correlation matrix.
#' @param n Sample size (scalar).
#' @param args Default arguments forwarded to both fits.
#' @param susieInfArgs SuSiE-inf-specific overrides.
#' @param susieArgs Standard SuSiE-RSS-specific overrides.
#' @param fittedModels Optional list with pre-fitted \code{$susie} and/or
#'   \code{$susie_inf} objects to skip re-fitting.
#' @return A list with \code{susie} and \code{susie_inf} fit objects.
#' @importFrom susieR susie_rss
#' @export
fitSusieInfThenSusieRss <- function(z, R, n, args = list(),
                                    susieInfArgs = list(),
                                    susieArgs = list(),
                                    fittedModels = NULL) {
  if (is.null(fittedModels)) fittedModels <- list()
  susieInfFit <- fittedModels[["susie_inf"]]
  susieFit <- fittedModels[["susie"]]

  if (is.null(susieInfFit)) {
    fitArgs <- modifyList(args, susieInfArgs)
    fitArgs <- modifyList(fitArgs, list(
      z = z, R = R, n = n, unmappable_effects = "inf",
      convergence_method = "pip", refine = FALSE, model_init = NULL
    ))
    susieInfFit <- do.call(susie_rss, fitArgs)
  }
  susieInfFit <- .setFinemappingFitClass(susieInfFit, "susie_inf")

  if (is.null(susieFit)) {
    fitArgs <- prepareSusieFromInfArgs(modifyList(args, susieArgs), susieInfFit, refineDefault = TRUE)
    susieFit <- do.call(susie_rss, c(list(z = z, R = R, n = n), fitArgs))
  }
  susieFit <- .setFinemappingFitClass(susieFit, "susie_rss")

  list(susie = susieFit, susie_inf = susieInfFit)
}

#' Post-process Fine-mapping Fits
#'
#' Applies method-aware post-processing to one or more SuSiE-family fits and
#' builds both a method-specific result list and shared top-loci tables.
#'
#' @param fits Named list of fine-mapping fits. Names define method identity,
#'   for example \code{susie}, \code{susie_inf}, \code{susie_rss},
#'   \code{mvsusie}, or \code{fsusie}.
#' @param dataX Genotype matrix, LD/correlation matrix, or other method-specific
#'   input used for credible-set purity and correlations.
#' @param dataY Phenotype vector/matrix or summary statistics. Default NULL.
#' @param xScalar Scaling factor for genotype effects. Default 1.
#' @param yScalar Scaling factor for phenotype effects. Default 1.
#' @param af Effect-allele frequencies (exported as the \code{af} column; never
#'   MAF). Default NULL.
#' @param coverage Primary credible-set coverage.
#' @param secondaryCoverage Additional credible-set coverages.
#' @param signalCutoff PIP cutoff for including non-CS variants in top loci.
#' @param otherQuantities Optional list carried into each method result.
#' @param priorEffTol Tolerance for retaining effects by prior variance.
#' @param minAbsCorr Minimum absolute correlation for credible-set purity.
#' @return A list with \code{finemapping_results} (per-method post-processed
#'   objects, each carrying a trimmed fit and method-specific intermediates)
#'   and a single unified \code{top_loci} table in the fixed 22-column shape
#'   (see \code{\link{buildTopLoci}}). Per-method contributions are
#'   row-bound into \code{top_loci} by an outer method for-loop.
#' @export
postprocessFinemappingFits <- function(fits, dataX, dataY = NULL,
                                       xScalar = 1, yScalar = 1,
                                       af = NULL, coverage = NULL,
                                       secondaryCoverage = c(0.7, 0.5),
                                       signalCutoff = 0.1,
                                       otherQuantities = NULL,
                                       region = NULL,
                                       priorEffTol = 1e-9,
                                       minAbsCorr = 0.8,
                                       csInput = NULL) {
  fits <- fits[!vapply(fits, is.null, logical(1))]
  if (length(fits) == 0) stop("At least one fine-mapping fit must be supplied.")
  if (is.null(names(fits)) || any(names(fits) == "")) {
    stop("fits must be a named list; names define method identity.")
  }

  # One method for-loop: each method calls buildTopLoci() once per fit; the
  # per-method 22-column contributions are row-bound below into the single
  # final `top_loci` table. There is no separately exposed long or wide table.
  posts <- lapply(names(fits), function(method) {
    fit <- .setFinemappingFitClass(fits[[method]], method)
    postprocessFinemappingFit(
      fit, method = method, dataX = dataX, dataY = dataY,
      xScalar = xScalar, yScalar = yScalar, af = af,
      coverage = coverage, secondaryCoverage = secondaryCoverage,
      signalCutoff = signalCutoff, otherQuantities = otherQuantities,
      region = region,
      priorEffTol = priorEffTol, minAbsCorr = minAbsCorr,
      csInput = csInput
    )
  })
  names(posts) <- names(fits)

  perMethod <- lapply(posts, function(x) x$top_loci)
  perMethod <- perMethod[!vapply(perMethod, is.null, logical(1))]
  topLoci <- if (length(perMethod) == 0L) {
    .emptyTopLoci()
  } else {
    do.call(rbind, perMethod)
  }
  rownames(topLoci) <- NULL
  posts <- lapply(posts, function(x) {
    x$top_loci <- NULL
    x
  })

  list(
    finemapping_results = posts,
    top_loci = topLoci
  )
}

postprocessFinemappingFit <- function(fit, ...) {
  UseMethod("postprocessFinemappingFit")
}

#' @exportS3Method
postprocessFinemappingFit.susie <- function(fit, method = "susie", csInput = NULL, ...) {
  if (is.null(csInput)) csInput <- "X"
  .postprocessFinemappingFitCommon(fit, method = method, csInput = csInput, ...)
}

#' @exportS3Method
postprocessFinemappingFit.susie_inf <- function(fit, method = "susie_inf", csInput = NULL, ...) {
  if (is.null(csInput)) csInput <- "X"
  .postprocessFinemappingFitCommon(fit, method = method, csInput = csInput, ...)
}

#' @exportS3Method
postprocessFinemappingFit.susie_rss <- function(fit, method = "susie_rss", csInput = NULL, ...) {
  if (is.null(csInput)) csInput <- "Xcorr"
  .postprocessFinemappingFitCommon(fit, method = method, csInput = csInput, ...)
}

#' @exportS3Method
postprocessFinemappingFit.mvsusie <- function(fit, method = "mvsusie", csInput = NULL, ...) {
  if (is.null(csInput)) csInput <- "X"
  .postprocessFinemappingFitCommon(fit, method = method, csInput = csInput, ...)
}

#' @exportS3Method
postprocessFinemappingFit.susiF <- function(fit, method = "fsusie", csInput = NULL, ...) {
  if (is.null(csInput)) csInput <- "fsusie"
  .postprocessFinemappingFitCommon(fit, method = method, csInput = csInput, ...)
}

.postprocessFinemappingFitCommon <- function(fit, method, dataX, dataY = NULL,
                                             xScalar = 1, yScalar = 1,
                                             af = NULL, coverage = NULL,
                                             secondaryCoverage = c(0.7, 0.5),
                                             signalCutoff = 0.1,
                                             otherQuantities = NULL,
                                             region = NULL,
                                             priorEffTol = 1e-9,
                                             minAbsCorr = 0.8,
                                             csInput = c("X", "Xcorr", "fsusie")) {
  csInput <- match.arg(csInput)
  variantNames <- extractVariantNames(fit)
  sumstats <- extractSumstats(fit, dataX, dataY, xScalar, yScalar, method)
  effectIdx <- selectEffects(fit, priorEffTol)
  csTables <- computeCsTables(
    fit, dataX = dataX, coverage = coverage,
    secondaryCoverage = secondaryCoverage, method = method,
    csInput = csInput, minAbsCorr = minAbsCorr
  )
  topLoci <- buildTopLoci(
    fit, csTables, variantNames = variantNames, sumstats = sumstats,
    af = af, method = method, signalCutoff = signalCutoff,
    dataX = dataX, dataY = dataY, otherQuantities = otherQuantities,
    region = region
  )

  trimmed <- trimFinemappingFit(fit, effectIdx, method, csTables)

  # Build FineMappingResult S4 object. The S4 contract (validity check,
  # vcf_writer, getPIP, getCS) still expects `variant_id`, `pip`, and an
  # integer `cs` column on the slot. To avoid rippling renames into
  # AllClasses / AllMethods / vcf_writer for this change, we project the
  # new 22-column `top_loci` into the legacy slot shape here, in
  # susie_wrapper only. The wrapper-facing `top_loci` returned to callers
  # is unchanged.
  s4TopLoci <- .topLociForS4Slot(topLoci)
  fmResult <- FineMappingResult(
    variantNames = variantNames,
    trimmedFit = trimmed,
    topLoci = s4TopLoci,
    method = method,
    sumstats = sumstats
  )

  res <- list(
    top_loci = topLoci,
    finemapping_result = fmResult
  )
  if (!is.null(sumstats)) res$sumstats <- sumstats
  sampleNames <- .sampleNamesFromDataY(dataY)
  if (!is.null(sampleNames)) res$sample_names <- sampleNames
  if (method == "mvsusie" && !is.null(fit$outcome_names)) res$context_names <- fit$outcome_names
  analysisScript <- loadScript()
  if (analysisScript != "") res$analysis_script <- analysisScript
  if (!is.null(otherQuantities)) res$other_quantities <- otherQuantities
  res
}

extractVariantNames <- function(fit) {
  variantNames <- names(fit$pip)
  if (is.null(variantNames)) variantNames <- colnames(fit$alpha)
  if (is.null(variantNames)) variantNames <- paste0("variant_", seq_along(fit$pip))
  tryCatch(normalizeVariantId(variantNames), error = function(e) variantNames)
}

extractSumstats <- function(fit, dataX, dataY, xScalar = 1, yScalar = 1, method = "susie") {
  if (is.null(dataY)) return(NULL)
  if (method == "susie_rss") return(dataY)
  if (is.list(dataY) && !is.data.frame(dataY) &&
      any(c("betahat", "sebetahat", "z") %in% names(dataY))) {
    return(dataY)
  }
  if (is.null(dataX)) return(NULL)
  if (is.matrix(dataY) || is.data.frame(dataY)) {
    if (ncol(as.matrix(dataY)) != 1) return(NULL)
  }
  sumstats <- univariate_regression(dataX, dataY)
  yScalar <- if (is.null(yScalar) || all(yScalar == 1)) 1 else yScalar
  xScalar <- if (is.null(xScalar) || all(xScalar == 1)) 1 else xScalar
  sumstats$betahat <- sumstats$betahat * yScalar / xScalar
  sumstats$sebetahat <- sumstats$sebetahat * yScalar / xScalar
  sumstats
}

.sampleNamesFromDataY <- function(dataY) {
  if (is.null(dataY) || is.list(dataY)) return(NULL)
  rownames(as.matrix(dataY))
}

selectEffects <- function(fit, priorEffTol = 1e-9) {
  alpha <- .asEffectMatrix(fit$alpha)
  nEffects <- nrow(alpha)
  if (nEffects == 0) return(integer(0))
  if (!is.null(fit$V)) {
    which(fit$V > priorEffTol)
  } else {
    seq_len(nEffects)
  }
}

.asEffectMatrix <- function(x) {
  if (is.null(x)) return(matrix(numeric(0), nrow = 0))
  if (is.list(x) && !is.data.frame(x)) return(do.call(rbind, x))
  as.matrix(x)
}

.asLbfMatrix <- function(fit) {
  if (!is.null(fit$lbf_variable)) return(.asEffectMatrix(fit$lbf_variable))
  if (!is.null(fit$lBF)) return(.asEffectMatrix(fit$lBF))
  NULL
}

#' @importFrom susieR get_cs_correlation
#' @noRd
computeCsTables <- function(fit, dataX, coverage = NULL,
                            secondaryCoverage = c(0.7, 0.5),
                            method = "susie", csInput = c("X", "Xcorr", "fsusie"),
                            minAbsCorr = 0.8) {
  csInput <- match.arg(csInput)
  primaryCoverage <- coverage
  if (is.null(primaryCoverage)) primaryCoverage <- fit$sets$requested_coverage
  if (is.null(primaryCoverage)) primaryCoverage <- 0.95
  coverages <- unique(c(primaryCoverage, secondaryCoverage))
  coverages <- coverages[!is.na(coverages)]

  tables <- lapply(coverages, function(cov) {
    computeCsTable(fit, dataX, coverage = cov, csInput = csInput, minAbsCorr = minAbsCorr)
  })
  names(tables) <- vapply(coverages, formatCsColumn, character(1), method = method)
  attr(tables, "coverage") <- coverages
  tables
}

computeCsTable <- function(fit, dataX, coverage, csInput = c("X", "Xcorr", "fsusie"),
                           minAbsCorr = 0.8) {
  csInput <- match.arg(csInput)
  if (csInput == "fsusie") {
    sets <- tryCatch(
      fsusieGetCs(fit, dataX, requestedCoverage = coverage),
      error = function(e) list(cs = list(), requested_coverage = coverage)
    )
    if (is.null(sets$cs) || length(sets$cs) == 0 || all(vapply(sets$cs, is.null, logical(1)))) {
      sets$cs <- list()
      return(list(sets = sets, cs_corr = NULL, pip = fit$pip))
    }
    tmp <- fit
    tmp$sets <- sets
    csCorr <- if (requireNamespace("fsusieR", quietly = TRUE)) {
      tryCatch(fsusieR::cal_cor_cs(tmp, dataX), error = function(e) NULL)
    } else {
      NULL
    }
    return(list(sets = sets, cs_corr = csCorr, pip = fit$pip))
  }

  if (csInput == "X") {
    sets <- susie_get_cs(fit, X = dataX, coverage = coverage, min_abs_corr = minAbsCorr)
    out <- list(sets = sets, pip = fit$pip)
    out$cs_corr <- get_cs_correlation(out, X = dataX)
  } else {
    sets <- susie_get_cs(fit, Xcorr = dataX, coverage = coverage, min_abs_corr = minAbsCorr)
    out <- list(sets = sets, pip = fit$pip)
    out$cs_corr <- get_cs_correlation(out, Xcorr = dataX)
  }
  out
}

#' Build the unified top-loci table for one fit and one method.
#'
#' Returns the per-fit, per-method contribution to the unified \code{top_loci}
#' table in the fixed 22-column shape. \code{postprocessFinemappingFits()}
#' calls this once per method per fit and row-binds the results into the
#' single \code{top_loci} returned by \code{formatFinemappingOutput()}.
#'
#' Output columns, in order: \code{#chr}, \code{start}, \code{end}, \code{a1},
#' \code{a2}, \code{variant}, \code{gene}, \code{event}, \code{n}, \code{af},
#' \code{beta}, \code{se}, \code{pip}, \code{posterior_effect_mean},
#' \code{posterior_effect_se}, \code{cs_95}, \code{cs_70}, \code{cs_50},
#' \code{cs_95_purity}, \code{method}, \code{grange_start}, \code{grange_end}.
#'
#' \code{cs_95} / \code{cs_70} / \code{cs_50} are character strings of the
#' form \code{"<method>_<cs_index>"} where each method numbers credible sets
#' independently from 1. Variants retained by the PIP cutoff but not in any
#' credible set at a coverage carry \code{"<method>_0"}. \code{cs_95_purity}
#' is the 0.95-coverage purity for the row's \code{(method, cs_95)}; rows
#' whose \code{cs_95} is \code{"<method>_0"} carry \code{0}.
#'
#' Row uniqueness is \code{(variant, gene, cs_membership)} at the given
#' \code{method}; overlapping CS within the same method produces one row per
#' CS.
#'
#' @param fit Fitted SuSiE-family object (must expose \code{alpha},
#'   \code{mu}, \code{mu2}, \code{pip}).
#' @param csTables List of CS tables (one per coverage) from
#'   \code{computeCsTables()}.
#' @param variantNames Character vector of variant IDs
#'   (\code{chr:pos:A2:A1}).
#' @param sumstats Optional marginal-association summary (\code{betahat},
#'   \code{sebetahat}) filling \code{beta} / \code{se}.
#' @param af Optional numeric vector of effect-allele frequencies (frequency of
#'   the final effect allele / \code{a1} after allele harmonization against the
#'   LD/reference variants). Exported directly as the \code{af} column. MAF is
#'   never exported; derive it from \code{af} at filter time. Default NULL ->
#'   \code{af = NA_real_}.
#' @param method Method name (e.g. \code{"susie"}, \code{"susie_inf"}). Required.
#' @param signalCutoff PIP cutoff for retaining PIP-only (non-CS) variants.
#' @param dataX Optional regional genotype matrix.
#' @param dataY Optional regional phenotype matrix; \code{nrow(dataY)} fills
#'   \code{n}, \code{colnames(dataY)[1]} fills \code{gene}.
#' @param otherQuantities Optional list. Default is NULL.
#' @param region Optional \code{"chr:start-end"} string. Default is NULL.
#' @return A data frame in the fixed 22-column shape for this fit and method,
#'   or an empty data frame if nothing is retained.
#' @export
buildTopLoci <- function(fit, csTables, variantNames, sumstats = NULL,
                         af = NULL, method, signalCutoff = 0.1,
                         dataX = NULL, dataY = NULL,
                         otherQuantities = NULL,
                         region = NULL) {
  if (missing(method) || is.null(method) ||
      length(method) != 1L || is.na(method) || !nzchar(method)) {
    stop("buildTopLoci: `method` is required (e.g. \"susie\", \"susie_inf\").")
  }
  if (length(csTables) == 0) return(.emptyTopLoci())
  coverageValues <- attr(csTables, "coverage")
  if (is.null(coverageValues)) coverageValues <- rep(NA_real_, length(csTables))

  # Per-fit constants.
  dataYMat <- if (!is.null(dataY)) as.matrix(dataY) else NULL
  fitN    <- if (is.null(dataYMat)) NA_integer_ else as.integer(nrow(dataYMat))
  fitGene <- if (!is.null(dataYMat) && !is.null(colnames(dataYMat))) {
    colnames(dataYMat)[1]
  } else NA_character_
  fitEvent <- if (!is.null(otherQuantities$condition_id) &&
                   !is.na(fitGene) && nzchar(fitGene)) {
    paste(otherQuantities$condition_id, fitGene, sep = "_")
  } else NA_character_
  grange <- .parseGrange(region)

  # Per-variant posterior effect / SE, computed once across all variants.
  alpha <- as.matrix(fit$alpha)
  mu    <- if (!is.null(fit$mu))  as.matrix(fit$mu)  else NULL
  mu2   <- if (!is.null(fit$mu2)) as.matrix(fit$mu2) else NULL
  postMean <- if (!is.null(mu) && all(dim(alpha) == dim(mu))) {
    colSums(alpha * mu)
  } else rep(NA_real_, length(variantNames))
  postSe <- if (!is.null(mu2) && all(dim(alpha) == dim(mu2))) {
    sqrt(pmax(colSums(alpha * mu2) - postMean^2, 0))
  } else rep(NA_real_, length(variantNames))

  # Collect CS-membership records (variant_idx, cs_idx, coverage) across all
  # requested coverages. This is the only intermediate; the 22-column shape
  # is projected from it below.
  csRecords <- do.call(rbind, lapply(seq_along(csTables), function(i) {
    ct <- csTables[[i]]
    info <- getCsInfo(ct$sets$cs, getTopVariantsIdx(ct, signalCutoff))
    if (is.null(info) || nrow(info) == 0) return(NULL)
    data.frame(variant_idx = as.integer(info$variant_idx),
               cs_idx      = as.integer(info$cs_idx),
               coverage    = as.numeric(coverageValues[[i]]),
               stringsAsFactors = FALSE)
  }))
  if (is.null(csRecords) || nrow(csRecords) == 0) return(.emptyTopLoci())

  # Key grid: one row per (variant_idx, cs_idx). Overlapping CS membership
  # within this method is preserved as separate keys.
  keyGrid <- unique(csRecords[, c("variant_idx", "cs_idx"), drop = FALSE])
  rownames(keyGrid) <- NULL
  nKeys  <- nrow(keyGrid)
  keyStr <- paste(keyGrid$variant_idx, keyGrid$cs_idx, sep = ":")

  # For each requested coverage, which keys appear in csRecords at that
  # coverage? Returns the key's cs_idx if present, else 0L.
  idxAt <- function(cov) {
    at <- csRecords[abs(csRecords$coverage - cov) < 1e-12, , drop = FALSE]
    hits <- paste(at$variant_idx, at$cs_idx, sep = ":")
    ifelse(keyStr %in% hits, keyGrid$cs_idx, 0L)
  }
  idx95 <- idxAt(0.95); idx70 <- idxAt(0.70); idx50 <- idxAt(0.50)

  # Per-coverage CS purity vectors (indexed by 1-based CS index). Only the
  # 0.95-coverage purity is currently exported (as cs_95_purity); per-CS
  # purities for the other coverages are kept here for downstream / future
  # use even though they are not part of the 22-column output.
  purityPerCov <- lapply(csTables, .csPurityVec)
  cov95        <- which(abs(coverageValues - 0.95) < 1e-12)
  purity95     <- if (length(cov95) > 0L) purityPerCov[[cov95[1]]] else numeric()
  cs95Purity   <- vapply(idx95, function(i) {
    if (i <= 0L || i > length(purity95)) return(0)
    v <- purity95[i]; if (is.na(v)) 0 else as.numeric(v)
  }, numeric(1))

  vIdx          <- keyGrid$variant_idx
  variantIdVec <- variantNames[vIdx]
  parsed <- tryCatch(
    suppressWarnings(parseVariantId(variantIdVec)),
    error = function(e) stop("buildTopLoci: parseVariantId failed: ",
                             conditionMessage(e))
  )
  if (is.null(parsed) || nrow(parsed) != length(variantIdVec)) {
    stop("buildTopLoci: parseVariantId did not return one row per variant.")
  }
  invalid <- is.na(parsed$chrom) | is.na(parsed$pos) |
    is.na(parsed$A1) | !nzchar(parsed$A1) |
    is.na(parsed$A2) | !nzchar(parsed$A2)
  if (any(invalid)) {
    stop("buildTopLoci: parseVariantId produced invalid coordinates ",
         "for variant_id: ", variantIdVec[which(invalid)[[1]]])
  }
  pick <- function(x) if (is.null(x)) rep(NA_real_, nKeys) else x[vIdx]

  out <- data.frame(
    "#chr"                = parsed$chrom,
    start                 = as.integer(parsed$pos) - 1L,
    end                   = as.integer(parsed$pos),
    a1                    = parsed$A1,
    a2                    = parsed$A2,
    variant               = variantIdVec,
    gene                  = rep(fitGene, nKeys),
    event                 = rep(fitEvent, nKeys),
    n                     = rep(fitN, nKeys),
    af                    = pick(af),
    beta                  = pick(sumstats$betahat),
    se                    = pick(sumstats$sebetahat),
    pip                   = as.numeric(fit$pip[vIdx]),
    posterior_effect_mean = postMean[vIdx],
    posterior_effect_se   = postSe[vIdx],
    cs_95                 = paste0(method, "_", idx95),
    cs_70                 = paste0(method, "_", idx70),
    cs_50                 = paste0(method, "_", idx50),
    cs_95_purity          = cs95Purity,
    method                = rep(method, nKeys),
    grange_start          = rep(grange[["start"]], nKeys),
    grange_end            = rep(grange[["end"]],   nKeys),
    stringsAsFactors      = FALSE,
    check.names           = FALSE
  )
  rownames(out) <- NULL
  out
}

# Per-CS purity from one cs_table: prefer susieR's sets$purity$min.abs.corr;
# fall back to cs_corr when purity is unavailable.
.csPurityVec <- function(ct) {
  sp <- ct$sets$purity
  if (!is.null(sp) && "min.abs.corr" %in% names(sp)) {
    return(as.numeric(sp$min.abs.corr))
  }
  if (!is.null(ct$cs_corr)) {
    return(vapply(ct$cs_corr, function(m) {
      if (is.null(m)) return(NA_real_)
      if (!is.matrix(m) || nrow(m) <= 1) return(1)
      min(abs(m[upper.tri(m)]))
    }, numeric(1)))
  }
  rep(NA_real_, length(ct$sets$cs))
}

.emptyTopLoci <- function() {
  data.frame(
    "#chr"                = character(),
    start                 = integer(),
    end                   = integer(),
    a1                    = character(),
    a2                    = character(),
    variant               = character(),
    gene                  = character(),
    event                 = character(),
    n                     = integer(),
    af                    = numeric(),
    beta                  = numeric(),
    se                    = numeric(),
    pip                   = numeric(),
    posterior_effect_mean = numeric(),
    posterior_effect_se   = numeric(),
    cs_95                 = character(),
    cs_70                 = character(),
    cs_50                 = character(),
    cs_95_purity          = numeric(),
    method                = character(),
    grange_start          = integer(),
    grange_end            = integer(),
    stringsAsFactors      = FALSE,
    check.names           = FALSE
  )
}

.parseGrange <- function(regionStr) {
  if (is.null(regionStr) || length(regionStr) == 0L ||
      is.na(regionStr) || !nzchar(as.character(regionStr))) {
    return(c(start = NA_integer_, end = NA_integer_))
  }
  pr <- tryCatch(parseRegion(as.character(regionStr)),
                 error = function(e) NULL)
  if (is.null(pr) || !is.data.frame(pr)) {
    return(c(start = NA_integer_, end = NA_integer_))
  }
  c(start = as.integer(pr$start), end = as.integer(pr$end))
}

# Project the new 22-column `top_loci` into the legacy shape expected by the
# FineMappingResult S4 slot, vcf_writer, getPIP, and getCS. We add backward-
# compatible aliases without renaming any column in the wrapper-facing
# `top_loci`:
#
#   * `variant_id` — copy of `variant`
#   * `cs`         — integer credible-set index derived from `cs_95` strings of
#                    the form `<method>_<idx>` (PIP-only `<method>_0` -> 0L)
#
# This isolates the schema change to susie_wrapper.R so AllClasses.R,
# AllMethods.R, and vcf_writer.R do not have to change.
.topLociForS4Slot <- function(topLoci) {
  if (is.null(topLoci) || nrow(topLoci) == 0) {
    return(data.frame(variant_id = character(0),
                      method     = character(0),
                      stringsAsFactors = FALSE))
  }
  out <- topLoci
  if ("variant" %in% names(out) && !"variant_id" %in% names(out)) {
    out$variant_id <- out$variant
  }
  if ("cs_95" %in% names(out) && !"cs" %in% names(out)) {
    out$cs <- vapply(out$cs_95, function(s) {
      if (is.na(s) || !nzchar(s)) return(0L)
      tailStr <- sub("^.*_", "", s)
      suppressWarnings(as.integer(tailStr))
    }, integer(1))
    out$cs[is.na(out$cs)] <- 0L
  }
  out
}

trimFinemappingFit <- function(fit, effectIdx, method, csTables) {
  alpha <- .asEffectMatrix(fit$alpha)
  lbfVariable <- .asLbfMatrix(fit)
  primary <- csTables[[1]]
  secondary <- if (length(csTables) > 1) {
    lapply(csTables[-1], function(x) x[names(x) != "pip"])
  } else {
    NULL
  }

  trimmed <- list(
    pip = as.numeric(fit$pip),
    sets = primary$sets,
    cs_corr = primary$cs_corr,
    sets_secondary = secondary,
    alpha = alpha[effectIdx, , drop = FALSE],
    lbf_variable = if (!is.null(lbfVariable)) lbfVariable[effectIdx, , drop = FALSE] else NULL,
    V = if (!is.null(fit$V)) fit$V[effectIdx] else NULL,
    niter = fit$niter,
    max_L = nrow(alpha),
    n_effects = nrow(alpha)
  )

  if (!is.null(fit$X_column_scale_factors)) trimmed$X_column_scale_factors <- fit$X_column_scale_factors
  if (!is.null(fit$mu)) {
    trimmed$mu <- if (length(dim(fit$mu)) == 3) fit$mu[effectIdx, , , drop = FALSE] else fit$mu[effectIdx, , drop = FALSE]
  }
  if (!is.null(fit$mu2)) {
    # mu2 is L x p for univariate susie and L x p x R for multivariate (mvsusie).
    # Match the shape handling used for mu just above.
    trimmed$mu2 <- if (length(dim(fit$mu2)) == 3) fit$mu2[effectIdx, , , drop = FALSE] else fit$mu2[effectIdx, , drop = FALSE]
  }
  if (!is.null(fit$theta)) trimmed$theta <- fit$theta
  if (!is.null(fit$omega_weights)) trimmed$omega_weights <- fit$omega_weights

  if (method == "mvsusie") {
    if (!is.null(fit$mu2_diag)) trimmed$mu2_diag <- fit$mu2_diag[effectIdx, , , drop = FALSE]
    if (requireNamespace("mvsusieR", quietly = TRUE)) {
      trimmed$coef <- mvsusieR::coef.mvsusie(fit)[-1, , drop = FALSE]
    }
    if (!is.null(fit$conditional_lfsr)) trimmed$clfsr <- fit$conditional_lfsr[effectIdx, , , drop = FALSE]
  }

  class(trimmed) <- unique(c(method, "susie"))
  trimmed
}

#' Format Fine-mapping Post-processing for Protocol Output
#'
#' Converts method-aware fine-mapping post-processing output into the root-level
#' fields consumed by protocol RDS files. The primary method's
#' \code{FineMappingResult} S4 object is promoted to the \code{finemapping_result}
#' field; use its accessors (\code{getTrimmedFit}, \code{getVariantNames},
#' \code{getTopLoci}, etc.) instead of legacy list keys.
#'
#' @param post Output from \code{\link{postprocessFinemappingFits}}.
#' @param primaryMethod Method whose result should populate root-level fields.
#' @return A list with root-level fields including \code{finemapping_result}
#'   and \code{top_loci}.
#' @export
formatFinemappingOutput <- function(post, primaryMethod) {
  methodPost <- post$finemapping_results[[primaryMethod]]
  if (is.null(methodPost)) {
    stop("primaryMethod was not found in finemapping_results: ", primaryMethod)
  }
  c(
    methodPost,
    list(
      top_loci = post$top_loci
    )
  )
}

#' Adjust SuSiE Weights
#'
#' Adjusts SuSiE TWAS weights by subsetting to intersected variants and
#' optionally running allele QC against LD reference variants.
#'
#' @param twasWeightsResults A list containing TWAS weight data (nested structure).
#' @param keepVariants Vector of variant names to keep.
#' @param runAlleleQc Whether to run allele_qc to align alleles. Default TRUE.
#' @param variableNameObj Path to variant names in the nested list.
#' @param susieObj Path to susie result in the nested list.
#' @param twasWeightsTable Path to weights table in the nested list.
#' @param ldVariants Vector of LD reference variant IDs for allele QC.
#' @param matchMinProp Minimum proportion of matched variants. Default 0.2.
#' @return A list with adjusted_susie_weights and remained_variants_ids.
#' @export
adjustSusieWeights <- function(twasWeightsResults, keepVariants, runAlleleQc = TRUE,
                               variableNameObj = c("susie_results", context, "variant_names"),
                               susieObj = c("susie_results", context, "susie_result_trimmed"),
                               twasWeightsTable = c("weights", context), ldVariants, matchMinProp = 0.2) {
  # Intersect the rownames of weights with keepVariants
  twasWeightsVariants <- getNestedElement(twasWeightsResults, variableNameObj)
  # Normalize to canonical format (with chr prefix)
  twasWeightsVariants <- normalizeVariantId(twasWeightsVariants)
  # allele flip twas weights matrix variants name
  if (runAlleleQc) {
    weightsMatrix <- getNestedElement(twasWeightsResults, twasWeightsTable)
    if (!all(c("chrom", "pos", "A2", "A1") %in% colnames(weightsMatrix))) {
      weightsMatrix <- cbind(parseVariantId(twasWeightsVariants), weightsMatrix)
    }
    weightsMatrixQced <- matchRefPanel(weightsMatrix, ldVariants, colnames(weightsMatrix)[!colnames(weightsMatrix) %in% c(
      "chrom",
      "pos", "A2", "A1"
    )], matchMinProp = matchMinProp)
    # matchRefPanel outputs canonical variant_ids (with chr prefix)
    qcSummaryDf <- getQcSummary(weightsMatrixQced)
    originalIdx <- match(qcSummaryDf$variants_id_original, twasWeightsVariants)
    intersectedIndices <- originalIdx[qcSummaryDf$keep == TRUE]
  } else {
    # Normalize keepVariants to canonical format for matching
    keepVariantsNormalized <- normalizeVariantId(keepVariants)
    intersectedVariants <- intersect(twasWeightsVariants, keepVariantsNormalized)
    intersectedIndices <- match(intersectedVariants, twasWeightsVariants)
  }
  if (length(intersectedIndices) == 0) {
    stop("Error: No intersected variants found. Please check 'twas_weights' and 'keep_variants' inputs to make sure there are variants left to use.")
  }
  # Subset lbf_matrix, mu, and x_column_scale_factors
  lbfMatrix <- getNestedElement(twasWeightsResults, c(susieObj, "lbf_variable"))
  mu <- getNestedElement(twasWeightsResults, c(susieObj, "mu"))
  xColumnScalFactors <- getNestedElement(twasWeightsResults, c(susieObj, "X_column_scale_factors"))

  lbfMatrixSubset <- lbfMatrix[, intersectedIndices, drop = FALSE]
  muSubset <- mu[, intersectedIndices, drop = FALSE]
  xColumnScalFactorsSubset <- xColumnScalFactors[intersectedIndices]

  # Convert lbf_matrix to alpha and calculate adjusted xQTL coefficients
  adjustedXqtlAlpha <- lbfToAlpha(lbfMatrixSubset)
  adjustedXqtlCoef <- colSums(adjustedXqtlAlpha * muSubset) / xColumnScalFactorsSubset
  # alleleQc now outputs canonical variant_ids (with chr prefix) -- no need to add chr
  remainedVariantsIds <- if (runAlleleQc) {
    getHarmonizedData(weightsMatrixQced)$variant_id
  } else {
    intersectedVariants
  }
  return(list(adjusted_susie_weights = adjustedXqtlCoef, remained_variants_ids = remainedVariantsIds))
}

#' Run the SuSiE RSS pipeline
#'
#' Runs SuSiE RSS analysis with one or more SuSiE-family variants. Supports
#' both z+R (correlation matrix) and z+X (genotype matrix) interfaces.
#'
#' @param sumstats Data frame with 'z' or ('beta' and 'se') columns.
#' @param ldMat LD correlation matrix. Mutually exclusive with xMat.
#' @param xMat Genotype matrix (samples x variants). Mutually exclusive with ldMat.
#' @param n Sample size.
#' @param L Maximum number of causal configurations (default: 30).
#' @param lGreedy Initial greedy number of causal configurations (default: 5).
#' @param analysisMethod Iteration mode for the \code{"susie_rss"} fit:
#'   \code{"susie_rss"} (default, normal IBSS), \code{"single_effect"} (L=1,
#'   single iteration), or \code{"bayesian_conditional_regression"}
#'   (full L, single iteration). Only affects the \code{"susie_rss"}
#'   method; ignored for \code{"susie_inf_rss"} and \code{"susie_ash_rss"}.
#' @param methods Optional character vector selecting which RSS variants to
#'   fit. Any subset of \code{c("susie_rss", "susie_inf_rss",
#'   "susie_ash_rss")}. Default \code{NULL} falls back to a single-method fit
#'   driven by \code{analysisMethod} (backward-compatible behavior). When
#'   \code{methods} is passed explicitly, each requested method is fitted;
#'   if \code{"susie_inf_rss"} is paired with \code{"susie_rss"} or
#'   \code{"susie_ash_rss"} (or both) and \code{addSusieInf = TRUE}, the
#'   SuSiE-inf-RSS fit initialises the downstream method. This exposes five
#'   distinct fitting modes mirroring the individual-level pipeline.
#' @param addSusieInf Logical. When \code{methods} contains
#'   \code{"susie_inf_rss"} alongside \code{"susie_rss"} and/or
#'   \code{"susie_ash_rss"}, controls whether SuSiE-inf-RSS is chained into
#'   the downstream method(s) as initialisation. Default \code{TRUE}.
#' @param coverage Coverage level (default: 0.95).
#' @param secondaryCoverage Secondary coverage levels (default: c(0.7, 0.5)).
#' @param signalCutoff PIP cutoff for selecting top loci (default: 0.1).
#' @param minAbsCorr Minimum absolute correlation for CS purity (default: 0.8).
#' @param rFinite Controls variance inflation to account for estimating
#'   the R matrix from a finite reference panel. NULL (default): no
#'   variance inflation. Passed directly to susie_rss.
#' @param rMismatch LD mismatch correction method passed directly to susie_rss.
#'   Default NULL disables mismatch correction.
#' @param ... Additional parameters passed to susie_rss. Supplying
#'   \code{var_y} here, together with \code{beta} and \code{se} columns in
#'   \code{sumstats}, selects the \code{bhat/shat/var_y} sufficient-statistic
#'   interface. Without \code{var_y}, this wrapper uses the z-score RSS
#'   interface. For binary traits, \code{rss_analysis_pipeline()} can compute
#'   the observed-scale OLS \code{var_y} automatically via
#'   \code{binary_trait_model = "ols"}; ordinary RSS mode leaves it absent.
#' @return A list with post-processed SuSiE RSS results. The unified
#'   \code{top_loci} table contains rows from every requested method,
#'   distinguished by the \code{method} column.
#' @importFrom susieR susie_rss
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange select
#' @export
susieRssPipeline <- function(sumstats, ldMat = NULL, xMat = NULL, n = NULL,
                             L = 30, lGreedy = 5,
                             analysisMethod = c("susie_rss", "single_effect", "bayesian_conditional_regression"),
                             methods = NULL,
                             addSusieInf = TRUE,
                             coverage = 0.95,
                             secondaryCoverage = c(0.7, 0.5),
                             signalCutoff = 0.1,
                             minAbsCorr = 0.8,
                             rFinite = NULL, rMismatch = NULL, ...) {
  analysisMethod <- match.arg(analysisMethod)
  if (is.null(ldMat) && is.null(xMat)) stop("Either ldMat or xMat must be provided.")
  if (!is.null(ldMat) && !is.null(xMat)) stop("Only one of ldMat or xMat should be provided, not both.")
  if (!is.null(lGreedy)) lGreedy <- min(lGreedy, L)

  # Resolve effective methods. NULL => legacy single-method via analysisMethod.
  validRssMethods <- c("susie_rss", "susie_inf_rss", "susie_ash_rss")
  if (is.null(methods)) {
    # Backward-compatible: single fit using analysisMethod, labeled accordingly.
    fitMethods <- analysisMethod
  } else {
    if (!is.character(methods) || length(methods) == 0L) {
      stop("methods must be a non-empty character vector of method names.")
    }
    bad <- setdiff(methods, validRssMethods)
    if (length(bad) > 0) {
      stop("Unknown RSS method(s): ", paste(bad, collapse = ", "),
           ". Valid options: ", paste(validRssMethods, collapse = ", "))
    }
    fitMethods <- unique(methods)
  }
  chainInfToSusieRss     <- isTRUE(addSusieInf) &&
    all(c("susie_inf_rss", "susie_rss") %in% fitMethods)
  chainInfToSusieAshRss <- isTRUE(addSusieInf) &&
    all(c("susie_inf_rss", "susie_ash_rss") %in% fitMethods)
  anyChainedInitRss <- chainInfToSusieRss || chainInfToSusieAshRss

  if (!is.null(sumstats$z)) {
    z <- sumstats$z
  } else if (!is.null(sumstats$beta) && !is.null(sumstats$se)) {
    z <- sumstats$beta / sumstats$se
  } else {
    stop("sumstats must have 'z' or ('beta' and 'se') columns.")
  }
  if (is.null(names(z)) && !is.null(sumstats$variant_id) && length(sumstats$variant_id) == length(z)) {
    names(z) <- sumstats$variant_id
  }
  if (is.null(names(z)) && !is.null(rownames(sumstats)) && length(rownames(sumstats)) == length(z)) {
    names(z) <- rownames(sumstats)
  }

  dots <- list(...)
  varY <- dots$varY
  dots$varY <- NULL
  if (!is.null(dots$bhat) || !is.null(dots$shat)) {
    stop("Pass summary effects as 'beta' and 'se' columns in sumstats; ",
         "susieRssPipeline constructs bhat and shat internally.")
  }
  if (!is.null(varY)) {
    if (is.null(sumstats$beta) || is.null(sumstats$se)) {
      stop("Supplying varY requires sumstats columns 'beta' and 'se'.")
    }
    if (isTRUE(attr(sumstats, "pecotmr_beta_se_from_z"))) {
      stop("Supplying varY requires observed beta and se columns; this ",
           "sumstats object has beta/se placeholders derived from z-scores.")
    }
    if (length(varY) != 1 || is.na(varY) || !is.finite(varY) ||
        varY <= 0) {
      stop("varY must be a positive finite scalar.")
    }
    if (is.list(xMat) && !is.matrix(xMat)) {
      stop("varY is not supported with multi-panel or list-backed xMat. ",
           "Use the z-score RSS interface instead.")
    }
    bhat <- sumstats$beta
    shat <- sumstats$se
    names(bhat) <- names(shat) <- names(z)
    common <- c(list(bhat = bhat, shat = shat, var_y = varY, n = n,
                     coverage = coverage, R_finite = rFinite,
                     R_mismatch = rMismatch), dots)
  } else {
    common <- c(list(z = z, n = n, coverage = coverage,
                     R_finite = rFinite, R_mismatch = rMismatch), dots)
  }
  if (!is.null(xMat)) common$X <- xMat else common$R <- ldMat

  fitOneSusieRss <- function() {
    if (analysisMethod == "single_effect") {
      do.call(susie_rss, c(common, list(L = 1, L_greedy = NULL, max_iter = 1)))
    } else if (analysisMethod == "bayesian_conditional_regression") {
      do.call(susie_rss, c(common, list(L = L, L_greedy = lGreedy, max_iter = 1)))
    } else {
      do.call(susie_rss, c(common, list(L = L, L_greedy = lGreedy)))
    }
  }
  fitOneSusieInfRss <- function() {
    do.call(susie_rss, c(common, list(L = L, L_greedy = lGreedy,
                                       unmappable_effects = "inf",
                                       convergence_method = "pip",
                                       refine = FALSE, model_init = NULL)))
  }
  fitOneSusieAshRss <- function() {
    do.call(susie_rss, c(common, list(L = L, L_greedy = lGreedy,
                                       unmappable_effects = "ash",
                                       convergence_method = "pip")))
  }

  fittedModels <- list()
  if ("susie_inf_rss" %in% fitMethods || anyChainedInitRss) {
    infFit <- fitOneSusieInfRss()
    fittedModels[["susie_inf_rss"]] <- .setFinemappingFitClass(infFit, "susie_inf_rss")
  }
  if ("susie_rss" %in% fitMethods ||
      identical(fitMethods, "single_effect") ||
      identical(fitMethods, "bayesian_conditional_regression")) {
    if (chainInfToSusieRss) {
      chainedArgs <- prepareSusieFromInfArgs(
        list(L = L, L_greedy = lGreedy),
        fittedModels[["susie_inf_rss"]], refineDefault = TRUE,
        unmappableEffects = "none"
      )
      rssFit <- do.call(susie_rss, c(common, chainedArgs))
    } else {
      rssFit <- fitOneSusieRss()
    }
    # Label by analysisMethod when in legacy single-method mode, else "susie_rss"
    rssLabel <- if (is.null(methods)) analysisMethod else "susie_rss"
    fittedModels[[rssLabel]] <- .setFinemappingFitClass(rssFit, rssLabel)
  }
  if ("susie_ash_rss" %in% fitMethods) {
    if (chainInfToSusieAshRss) {
      chainedArgs <- prepareSusieFromInfArgs(
        list(L = L, L_greedy = lGreedy),
        fittedModels[["susie_inf_rss"]], refineDefault = NULL,
        unmappableEffects = "ash"
      )
      ashFit <- do.call(susie_rss, c(common, chainedArgs))
    } else {
      ashFit <- fitOneSusieAshRss()
    }
    fittedModels[["susie_ash_rss"]] <- .setFinemappingFitClass(ashFit, "susie_ash_rss")
  }

  # Drop SuSiE-inf-RSS from post-processing if it was only fit for init
  if (anyChainedInitRss && !("susie_inf_rss" %in% fitMethods)) {
    fittedModels[["susie_inf_rss"]] <- NULL
  }

  # For post-processing, pass genotype matrix X directly when available.
  if (!is.null(ldMat)) {
    dataX <- ldMat
    ppCsInput <- "Xcorr"
  } else if (is.list(xMat) && !is.matrix(xMat)) {
    dataX <- do.call(rbind, xMat)[, seq_along(z), drop = FALSE]
    ppCsInput <- "X"
  } else {
    dataX <- xMat[, seq_along(z), drop = FALSE]
    ppCsInput <- "X"
  }

  # Effect-allele frequency for top_loci$af. Carried only when the harmonized
  # sumstats declares it (effect-allele AF); aligned to the z / variant order.
  # MAF is never exported here; it is an internal QC quantity derived from af.
  af <- if (!is.null(sumstats$af)) as.numeric(sumstats$af) else NULL

  post <- postprocessFinemappingFits(
    fits = fittedModels,
    dataX = dataX,
    dataY = list(z = z),
    af = af,
    coverage = coverage,
    secondaryCoverage = secondaryCoverage,
    signalCutoff = signalCutoff,
    minAbsCorr = minAbsCorr,
    csInput = ppCsInput
  )
  # Primary method preference: "susie_rss" > other names > first fit
  primary <- if ("susie_rss" %in% names(fittedModels)) "susie_rss" else names(fittedModels)[1]
  formatFinemappingOutput(post, primaryMethod = primary)
}

#' @noRd
getCsIndex <- function(snpsIdx, susieCs) {
  # Return ALL CS indices that contain this variant (not just one)
  idx <- which(vapply(susieCs, function(x) snpsIdx %in% x, logical(1)))
  if (length(idx) == 0) return(NA_integer_)
  return(idx)
}
#' @noRd
getTopVariantsIdx <- function(susieOutput, signalCutoff) {
  c(which(susieOutput$pip >= signalCutoff), unlist(susieOutput$sets$cs)) %>%
    unique() %>%
    sort()
}
# Returns a data.frame(variant_idx, cs_idx) with one row per (variant, CS) pair.
# Variants in multiple CSs get multiple rows.
#' @importFrom stringr str_replace
#' @noRd
getCsInfo <- function(susieOutputSetsCs, topVariantsIdx) {
  csNames <- names(susieOutputSetsCs)
  rows <- lapply(topVariantsIdx, function(vi) {
    idx <- getCsIndex(vi, susieOutputSetsCs)
    if (length(idx) == 1 && is.na(idx)) {
      data.frame(variant_idx = vi, cs_idx = 0L, stringsAsFactors = FALSE)
    } else {
      csNums <- as.integer(str_replace(csNames[idx], "L", ""))
      data.frame(variant_idx = rep(vi, length(csNums)), cs_idx = csNums, stringsAsFactors = FALSE)
    }
  })
  do.call(rbind, rows)
}
