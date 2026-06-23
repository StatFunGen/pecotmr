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

# Translate a camelCase pecotmr method identifier (e.g. "susieInfRss") into the
# snake_case form (e.g. "susie_inf_rss") used in the documented top_loci schema.
# Single-word identifiers (e.g. "susie", "mvsusie", "fsusie") pass through.
.camelToSnakeMethod <- function(method) {
  if (is.null(method) || length(method) == 0L) return(method)
  lookup <- c(
    susieInf                      = "susie_inf",
    susieAsh                      = "susie_ash",
    susieRss                      = "susie_rss",
    susieInfRss                   = "susie_inf_rss",
    susieAshRss                   = "susie_ash_rss",
    singleEffect                  = "single_effect",
    bayesianConditionalRegression = "bayesian_conditional_regression"
  )
  vapply(method, function(m) {
    if (m %in% names(lookup)) lookup[[m]] else m
  }, character(1), USE.NAMES = FALSE)
}

.setFinemappingFitClass <- function(fit, method) {
  if (is.null(fit)) return(NULL)
  methodClass <- switch(method,
    susie = "susie",
    susieInf = "susieInf",
    susieRss = "susieRss",
    singleEffect = "susieRss",
    bayesianConditionalRegression = "susieRss",
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
  susieInfFit <- fittedModels[["susieInf"]]
  susieFit <- fittedModels[["susie"]]

  if (is.null(susieInfFit)) {
    fitArgs <- modifyList(args, susieInfArgs)
    fitArgs <- modifyList(fitArgs, list(
      X = X, y = y, unmappable_effects = "inf",
      convergence_method = "pip", refine = FALSE, model_init = NULL
    ))
    susieInfFit <- do.call(susie, fitArgs)
  }
  susieInfFit <- .setFinemappingFitClass(susieInfFit, "susieInf")

  if (is.null(susieFit)) {
    fitArgs <- prepareSusieFromInfArgs(modifyList(args, susieArgs), susieInfFit, refineDefault = TRUE)
    susieFit <- do.call(susie, c(list(X = X, y = y), fitArgs))
  }
  susieFit <- .setFinemappingFitClass(susieFit, "susie")

  list(susie = susieFit, susieInf = susieInfFit)
}

#' Two-stage SuSiE-RSS Fine-mapping
#'
#' RSS analog of \code{fitSusieInfThenSusie}. Fits SuSiE-inf via
#' \code{susieRss} first, then initialises standard SuSiE-RSS from
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
#'   \code{$susieInf} objects to skip re-fitting.
#' @return A list with \code{susie} and \code{susieInf} fit objects.
#' @importFrom susieR susie_rss
#' @export
fitSusieInfThenSusieRss <- function(z, R, n, args = list(),
                                    susieInfArgs = list(),
                                    susieArgs = list(),
                                    fittedModels = NULL) {
  if (is.null(fittedModels)) fittedModels <- list()
  susieInfFit <- fittedModels[["susieInf"]]
  susieFit <- fittedModels[["susie"]]

  if (is.null(susieInfFit)) {
    fitArgs <- modifyList(args, susieInfArgs)
    fitArgs <- modifyList(fitArgs, list(
      z = z, R = R, n = n, unmappable_effects = "inf",
      convergence_method = "pip", refine = FALSE, model_init = NULL
    ))
    susieInfFit <- do.call(susie_rss, fitArgs)
  }
  susieInfFit <- .setFinemappingFitClass(susieInfFit, "susieInf")

  if (is.null(susieFit)) {
    fitArgs <- prepareSusieFromInfArgs(modifyList(args, susieArgs), susieInfFit, refineDefault = TRUE)
    susieFit <- do.call(susie_rss, c(list(z = z, R = R, n = n), fitArgs))
  }
  susieFit <- .setFinemappingFitClass(susieFit, "susieRss")

  list(susie = susieFit, susieInf = susieInfFit)
}

#' Post-process Fine-mapping Fits
#'
#' Applies method-aware post-processing to one or more SuSiE-family fits and
#' builds both a method-specific result list and shared top-loci tables.
#'
#' @param fits Named list of fine-mapping fits. Names define method identity,
#'   for example \code{susie}, \code{susieInf}, \code{susieRss},
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
#' @return A list with \code{finemappingResults} (per-method post-processed
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
                                       medianAbsCorr = NULL,
                                       csInput = NULL,
                                       trim = TRUE) {
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
      medianAbsCorr = medianAbsCorr,
      csInput = csInput,
      trim = trim
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
    finemappingResults = posts,
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
postprocessFinemappingFit.susieInf <- function(fit, method = "susieInf", csInput = NULL, ...) {
  if (is.null(csInput)) csInput <- "X"
  .postprocessFinemappingFitCommon(fit, method = method, csInput = csInput, ...)
}

#' @exportS3Method
postprocessFinemappingFit.susieRss <- function(fit, method = "susieRss", csInput = NULL, ...) {
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
                                             trim = TRUE,
                                             minAbsCorr = 0.8,
                                             medianAbsCorr = NULL,
                                             csInput = c("X", "Xcorr", "fsusie")) {
  csInput <- match.arg(csInput)
  variantNames <- extractVariantNames(fit)
  sumstats <- extractSumstats(fit, dataX, dataY, xScalar, yScalar, method)
  effectIdx <- selectEffects(fit, priorEffTol)
  csTables <- computeCsTables(
    fit, dataX = dataX, coverage = coverage,
    secondaryCoverage = secondaryCoverage, method = method,
    csInput = csInput, minAbsCorr = minAbsCorr, medianAbsCorr = medianAbsCorr
  )
  # Always build the canonical unfiltered table; the FineMappingEntry
  # slot stores it as-is so accessors can filter by PIP at query time.
  # The wrapper-facing `top_loci` (in `res` below) preserves the legacy
  # `signalCutoff` behaviour for non-S4 callers.
  topLociFull <- buildTopLoci(
    fit, csTables, variantNames = variantNames, sumstats = sumstats,
    af = af, method = method, signalCutoff = 0,
    dataX = dataX, dataY = dataY, otherQuantities = otherQuantities,
    region = region
  )

  # When `trim = TRUE` we store a minimal subset of the fit on the
  # entry; when `trim = FALSE` we keep the full untrimmed susie return so
  # downstream code can access `mu` / `mu2` / `lbf_variable` / `V` / etc.
  storedFit <- if (isTRUE(trim)) {
    trimFinemappingFit(fit, effectIdx, method, csTables)
  } else {
    fit
  }

  fmEntry <- FineMappingEntry(
    variantIds = variantNames,
    susieFit   = storedFit,
    topLoci    = topLociFull)

  topLociWrapper <- topLociFull
  if (!is.null(signalCutoff) && signalCutoff > 0 && nrow(topLociWrapper) > 0L) {
    keep <- !is.na(topLociWrapper$pip) & topLociWrapper$pip > signalCutoff
    topLociWrapper <- topLociWrapper[keep, , drop = FALSE]
  }

  res <- list(
    top_loci = topLociWrapper,
    finemappingEntry = fmEntry,
    method = method
  )
  if (!is.null(sumstats)) res$sumstats <- sumstats
  sampleNames <- .sampleNamesFromDataY(dataY)
  if (!is.null(sampleNames)) res$sampleNames <- sampleNames
  if (method == "mvsusie" && !is.null(fit$outcome_names)) res$contextNames <- fit$outcome_names
  if (!is.null(otherQuantities)) res$otherQuantities <- otherQuantities
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
  if (method == "susieRss") return(dataY)
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
                            minAbsCorr = 0.8, medianAbsCorr = NULL) {
  csInput <- match.arg(csInput)
  primaryCoverage <- coverage
  if (is.null(primaryCoverage)) primaryCoverage <- fit$sets$requested_coverage
  if (is.null(primaryCoverage)) primaryCoverage <- 0.95
  coverages <- unique(c(primaryCoverage, secondaryCoverage))
  coverages <- coverages[!is.na(coverages)]

  tables <- lapply(coverages, function(cov) {
    computeCsTable(fit, dataX, coverage = cov, csInput = csInput,
                   minAbsCorr = minAbsCorr, medianAbsCorr = medianAbsCorr)
  })
  names(tables) <- vapply(coverages, formatCsColumn, character(1), method = method)
  attr(tables, "coverage") <- coverages
  tables
}

computeCsTable <- function(fit, dataX, coverage, csInput = c("X", "Xcorr", "fsusie"),
                           minAbsCorr = 0.8, medianAbsCorr = NULL) {
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

  # Purity thresholds for credible-set extraction. min_abs_corr / median_abs_corr
  # are isolated from finemappingOpts upstream and routed here; pass each only
  # when set. `fit` is passed positionally as `res`.
  csArgs <- list(coverage = coverage)
  if (!is.null(minAbsCorr)) csArgs$min_abs_corr <- minAbsCorr
  if (!is.null(medianAbsCorr)) csArgs$median_abs_corr <- medianAbsCorr
  if (csInput == "X") {
    sets <- do.call(susie_get_cs, c(list(fit), csArgs, list(X = dataX)))
    out <- list(sets = sets, pip = fit$pip)
    out$cs_corr <- get_cs_correlation(out, X = dataX)
  } else {
    sets <- do.call(susie_get_cs, c(list(fit), csArgs, list(Xcorr = dataX)))
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
#' @param method Method name (e.g. \code{"susie"}, \code{"susieInf"}). Required.
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
                         af = NULL, method, signalCutoff = 0,
                         dataX = NULL, dataY = NULL,
                         otherQuantities = NULL,
                         region = NULL) {
  if (missing(method) || is.null(method) ||
      length(method) != 1L || is.na(method) || !nzchar(method)) {
    stop("buildTopLoci: `method` is required (e.g. \"susie\", \"susieInf\").")
  }
  if (length(variantNames) == 0L) return(.emptyTopLoci())
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
  postSd <- if (!is.null(mu2) && all(dim(alpha) == dim(mu2))) {
    sqrt(pmax(colSums(alpha * mu2) - postMean^2, 0))
  } else rep(NA_real_, length(variantNames))

  # Parse variant IDs into chrom/pos/A1/A2 (one row per variant).
  parsed <- tryCatch(
    suppressWarnings(parseVariantId(variantNames)),
    error = function(e) stop("buildTopLoci: parseVariantId failed: ",
                             conditionMessage(e)))
  if (is.null(parsed) || nrow(parsed) != length(variantNames)) {
    stop("buildTopLoci: parseVariantId did not return one row per variant.")
  }
  invalid <- is.na(parsed$chrom) | is.na(parsed$pos) |
    is.na(parsed$A1) | !nzchar(parsed$A1) |
    is.na(parsed$A2) | !nzchar(parsed$A2)
  if (any(invalid)) {
    stop("buildTopLoci: parseVariantId produced invalid coordinates ",
         "for variant_id: ", variantNames[which(invalid)[[1]]])
  }

  # Marginal univariate effects (β, SE, Z, p). Per-variant; populated
  # uniformly across individual-level and RSS paths (the caller computes
  # the underlying sumstats list).
  nV <- length(variantNames)
  marginalBeta <- if (!is.null(sumstats$betahat))   as.numeric(sumstats$betahat)
                  else rep(NA_real_, nV)
  marginalSe   <- if (!is.null(sumstats$sebetahat)) as.numeric(sumstats$sebetahat)
                  else rep(NA_real_, nV)
  marginalZ    <- if (!is.null(sumstats$z))         as.numeric(sumstats$z)
                  else if (any(!is.na(marginalBeta)) && any(!is.na(marginalSe)))
                    marginalBeta / marginalSe
                  else rep(NA_real_, nV)
  marginalP    <- if (!is.null(sumstats$p))         as.numeric(sumstats$p)
                  else if (any(!is.na(marginalZ)))  2 * stats::pnorm(-abs(marginalZ))
                  else rep(NA_real_, nV)

  # Per-coverage CS membership: for each variant, which CS at each
  # coverage level (cs_idx, or 0 if not in any). If a variant belongs
  # to multiple CSs at a given coverage, the smallest cs_idx wins.
  csIdxAtCoverage <- function(targetCov) {
    out <- integer(nV)
    hit <- which(abs(coverageValues - targetCov) < 1e-12)
    if (length(hit) == 0L) return(out)
    sets <- csTables[[hit[1L]]]$sets$cs
    if (is.null(sets) || length(sets) == 0L) return(out)
    for (csIdx in seq_along(sets)) {
      vi <- as.integer(sets[[csIdx]])
      vi <- vi[vi >= 1L & vi <= nV & out[vi] == 0L]
      out[vi] <- csIdx
    }
    out
  }
  idx95 <- csIdxAtCoverage(0.95)
  idx70 <- csIdxAtCoverage(0.70)
  idx50 <- csIdxAtCoverage(0.50)

  # 0.95-coverage CS purity, per-variant (0 for non-CS variants).
  purityPerCs <- {
    h <- which(abs(coverageValues - 0.95) < 1e-12)
    if (length(h) > 0L) .csPurityVec(csTables[[h[1L]]]) else numeric()
  }
  cs95Purity <- vapply(idx95, function(i) {
    if (i <= 0L || i > length(purityPerCs)) return(0)
    v <- purityPerCs[i]; if (is.na(v)) 0 else as.numeric(v)
  }, numeric(1))

  methodTag <- .camelToSnakeMethod(method)
  out <- data.frame(
    variant_id     = as.character(variantNames),
    chrom          = parsed$chrom,
    pos            = as.integer(parsed$pos),
    A1             = parsed$A1,
    A2             = parsed$A2,
    N              = rep(fitN, nV),
    af             = if (is.null(af)) rep(NA_real_, nV) else as.numeric(af),
    marginal_beta  = marginalBeta,
    marginal_se    = marginalSe,
    marginal_z     = marginalZ,
    marginal_p     = marginalP,
    pip            = as.numeric(fit$pip),
    posterior_mean = postMean,
    posterior_sd   = postSd,
    cs_95          = paste0(methodTag, "_", idx95),
    cs_70          = paste0(methodTag, "_", idx70),
    cs_50          = paste0(methodTag, "_", idx50),
    cs_95_purity   = cs95Purity,
    method         = rep(method, nV),
    gene           = rep(fitGene, nV),
    event          = rep(fitEvent, nV),
    grange_start   = rep(grange[["start"]], nV),
    grange_end     = rep(grange[["end"]],   nV),
    stringsAsFactors = FALSE)
  if (!is.null(signalCutoff) && signalCutoff > 0) {
    keep <- !is.na(out$pip) & out$pip > signalCutoff
    out <- out[keep, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

# Translate susieR's snake-case `sets$purity` columns into pecotmr camelCase.
# Accepts a data.frame, matrix, or NULL; preserves type and column order.
.translateSusiePurity <- function(p) {
  if (is.null(p)) return(p)
  lookup <- c("min.abs.corr"    = "minAbsCorr",
              "mean.abs.corr"   = "meanAbsCorr",
              "median.abs.corr" = "medianAbsCorr")
  if (is.data.frame(p)) {
    nm <- names(p)
    names(p) <- ifelse(nm %in% names(lookup), lookup[nm], nm)
  } else if (is.matrix(p)) {
    cn <- colnames(p)
    if (!is.null(cn)) {
      colnames(p) <- ifelse(cn %in% names(lookup), lookup[cn], cn)
    }
  }
  p
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
    variant_id     = character(),
    chrom          = character(),
    pos            = integer(),
    A1             = character(),
    A2             = character(),
    N              = numeric(),
    af             = numeric(),
    marginal_beta  = numeric(),
    marginal_se    = numeric(),
    marginal_z     = numeric(),
    marginal_p     = numeric(),
    pip            = numeric(),
    posterior_mean = numeric(),
    posterior_sd   = numeric(),
    cs_95          = character(),
    cs_70          = character(),
    cs_50          = character(),
    cs_95_purity   = numeric(),
    method         = character(),
    gene           = character(),
    event          = character(),
    grange_start   = integer(),
    grange_end     = integer(),
    stringsAsFactors = FALSE
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
# This isolates the schema change to susieWrapper.R so allClasses.R,
# allMethods.R, and vcfWriter.R do not have to change.
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
#' Promotes the primary method's per-method post-processing payload to the
#' root level and attaches the unified \code{top_loci} table. The primary
#' method's bare \code{FineMappingEntry} appears at \code{$finemappingEntry};
#' wrap it into a \code{FineMappingResult} collection at the pipeline level
#' once (study, context, trait, method) identity tags are known.
#'
#' @param post Output from \code{\link{postprocessFinemappingFits}}.
#' @param primaryMethod Method whose result should populate root-level fields.
#' @return A list with root-level fields including \code{finemappingEntry}
#'   (a bare \code{FineMappingEntry} S4 payload) and \code{top_loci}.
#' @export
formatFinemappingOutput <- function(post, primaryMethod) {
  methodPost <- post$finemappingResults[[primaryMethod]]
  if (is.null(methodPost)) {
    stop("primaryMethod was not found in finemappingResults: ", primaryMethod)
  }
  c(
    methodPost,
    list(
      top_loci = post$top_loci
    )
  )
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
#' @title  Calculate Purity Measures for Credible Sets
#'
#' @description As an extension of the internal cal_purity function. This function computes purity metrics (minimum, mean, and median absolute correlations)
#' for each credible set in a list of credible set indices, based on the provided X matrix.
#' The output Purity depends on the method specified: for the 'min' method,
#' it returns a single value for single-element sets or the minimum absolute correlation for others.
#' For other methods, it returns a vector of three values (min, mean, median) for each set.
#'
#' @param lCs A list of credible set indices, where each element is a vector of indices
#'             corresponding to variables in a credible set.
#' @param X The data matrix used to compute correlations between variables in each credible set.
#' @param method A character string specifying the method to use for calculating purity.
#'               Defaults to 'min'. Other methods return a vector of min, mean, and median
#'               absolute correlations for each credible set.
#' @return A list where each element corresponds to a credible set and contains either a single
#'         purity value (for 'min' method and single-element sets) or a vector of purity metrics
#'         (for other methods and multi-element sets).
#' @noRd

calPurity <- function(lCs, X, method = "min") {
  tt <- list()

  for (k in seq_along(lCs)) {
    csIndices <- unlist(lCs[[k]])
    if (method == "min") {
      if (length(csIndices) == 1) {
        tt[[k]] <- 1
      } else {
        x <- abs(computeLd(X[, csIndices, drop = FALSE], method = "sample"))
        x[col(x) == row(x)] <- NA
        tt[[k]] <- min(x, na.rm = TRUE)
      }
    } else {
      if (length(csIndices) == 1) {
        tt[[k]] <- c(1, 1, 1)
      } else {
        x <- abs(computeLd(X[, csIndices, drop = FALSE], method = "sample"))
        x[col(x) == row(x)] <- NA
        tt[[k]] <- c(
          min(x, na.rm = TRUE),
          mean(x, na.rm = TRUE),
          median(x, na.rm = TRUE)
        )
      }
    }
  }

  return(tt)
}


#'  @title Create Sets Similar to SuSiE Output from fSuSiE Object
#'
#' @description This function constructs a list that mimics the structure of SuSiE output sets
#' from a fSuSiE object. It includes credible sets (cs) with their names, a purity
#' dataframe, coverage information, and the requested coverage level.
#'
#' @param fsusieObj A fSuSiE object containing the results from a fSuSiE analysis.
#' expected to at least have 'cs' and 'alpha' components.
#' @param requestedCoverage A numeric value specifying the desired coverage level for the
#'  credible sets. This is purely for record purpose so should be
#'  manually ensured that it correctly reflect the actual coverage used. Defaults to 0.95.
#' @return A list containing named credible sets (cs), a dataframe of purity metrics
#'         (minAbsCorr, meanAbsCorr, medianAbsCorr), an index of credible sets (cs_index),
#'         coverage values for each set, and the requested coverage level. Similar to the SuSiE set output
#' @export
fsusieGetCs <- function(fsusieObj, X, requestedCoverage = 0.95) {
  # Create 'cs' set with names
  csNamed <- setNames(object = fsusieObj$cs, nm = paste0("L", seq_along(fsusieObj$cs)))

  # Create 'purity' data frame
  purityDf <- do.call(rbind, lapply(calPurity(fsusieObj$cs, X = X, method = "susie"), function(x) as.data.frame(t(x))))
  rownames(purityDf) <- names(csNamed)
  colnames(purityDf) <- c("minAbsCorr", "meanAbsCorr", "medianAbsCorr")

  # Create 'coverage' without
  coverageVector <- numeric(length(fsusieObj$alpha))
  for (i in seq_along(fsusieObj$alpha)) {
    alphaI <- fsusieObj$alpha[[i]]
    csI <- fsusieObj$cs[[i]]
    coverageVector[i] <- sum(alphaI[csI])
  }

  # Combine all elements into a list
  sets <- list(
    cs = csNamed,
    purity = purityDf,
    cs_index = seq_along(fsusieObj$cs),
    coverage = coverageVector,
    requested_coverage = requestedCoverage
  )

  return(sets)
}

#' @title Wrapper for fsusie Function with Automatic Post-Processing
#'
#' @description This function serves as a wrapper for the fsusie function, facilitating
#' automatic post-processing such as removing dummy credible sets (cs) that don't meet
#' the minimum purity threshold and calculating correlations for the remaining cs.
#' The function parameters are identical to those of the fSuSiE function.
#'
#' @param X Residual genotype matrix.
#' @param Y Response phenotype matrix.
#' @param pos Genomics position of phenotypes, used for specifying the wavelet model.
#' @param L The maximum number of the credible set.
#' @param prior method to generate the prior.
#' @param maxSnpEm maximum number of SNP used for learning the prior.
#' @param covLev Coverage level for the credible sets.
#' @param maxScale numeric, define the maximum of wavelet coefficients used in the analysis (2^maxScale).
#'        Set 10 true by default.
#' @param minPurity Minimum purity threshold for credible sets to be retained.
#' @param ... Additional arguments passed to the fsusie function.
#' @return A modified fsusie object with the susie sets list, correlations for cs, alpha as df like susie,
#'         and without the dummy cs that do not meet the minimum purity requirement.
#' @export

fsusieWrapper <- function(X, Y, pos, L, prior, maxSnpEm, covLev, minPurity, maxScale, ...) {
  # Make sure fsusieR installed
  if (!requireNamespace("fsusieR", quietly = TRUE)) {
    stop("To use this function, please install fsusieR: https://github.com/stephenslab/fsusieR")
  }
  # Run fsusie
  fsusieObj <- fsusieR::susiF(
    X = X, Y = Y, pos = pos, L = L, prior = prior,
    max_SNP_EM = maxSnpEm, cov_lev = covLev,
    min_purity = minPurity, max_scale = maxScale, ...
  )

  # Remove dummy cs based on purity threshold
  if (all(abs(as.numeric(fsusieObj$purity)) < minPurity)) {
    fsusieObj$cs <- list(NULL)
    fsusieObj$sets <- list(cs = list(NULL), requested_coverage = covLev)
    fsusieObj$cs_corr <- NULL # Set cs correlations to NULL if no credible sets meet purity criteria
  } else {
    # Create sets and add correlation for CS if purity criteria are met
    fsusieObj$sets <- fsusieGetCs(fsusieObj, X, requestedCoverage = covLev)
    fsusieObj$cs_corr <- fsusieR::cal_cor_cs(fsusieObj, X)
  }
  # Put alpha into df
  fsusieObj$alpha <- do.call(rbind, lapply(fsusieObj$alpha, function(x) as.data.frame(t(x))))
  return(fsusieObj)
}



# =============================================================================
# Uniform fit wrappers for mvSuSiE (individual + RSS)
# -----------------------------------------------------------------------------
# Thin wrappers around mvsusieR::mvsusie and mvsusieR::mvsusie_rss. Every
# inline call across the package routes through these so the indirection
# is testable in one place and so future changes to the underlying mvsusieR
# API only need updating here.
# =============================================================================

#' Fit mvSuSiE on individual-level (X, Y) data
#'
#' Wrapper around \code{mvsusieR::mvsusie} with the canonical argument
#' set used inside fine-mapping and TWAS-weight pipelines.
#'
#' @param X Numeric matrix of genotypes (samples x variants).
#' @param Y Numeric matrix of multi-trait / multi-context outcomes
#'   (samples x conditions).
#' @param prior_variance Prior variance matrix; pass the output of
#'   \code{mvsusieR::create_mixture_prior(R = ncol(Y))} unless you have
#'   a domain-specific prior.
#' @param coverage Credible set coverage (default 0.95).
#' @param ... Additional arguments forwarded to
#'   \code{mvsusieR::mvsusie}.
#' @return The fit object returned by \code{mvsusieR::mvsusie}.
#' @export
fitMvsusie <- function(X, Y, prior_variance, coverage = 0.95, ...) {
  mvsusieR::mvsusie(X = X, Y = Y,
                    prior_variance = prior_variance,
                    coverage = coverage, ...)
}

#' Fit mvSuSiE-RSS on summary-statistic (Z, R, N) data
#'
#' Wrapper around \code{mvsusieR::mvsusie_rss}. The underlying function
#' was renamed from \code{mvsusieRss} to \code{mvsusie_rss} upstream;
#' this wrapper insulates pecotmr from that naming.
#'
#' @param Z Numeric matrix of Z-scores (variants x conditions).
#' @param R Variant-by-variant LD correlation matrix.
#' @param N Scalar sample size (median across conditions when N varies).
#' @param prior_variance Prior variance matrix.
#' @param coverage Credible set coverage (default 0.95).
#' @param ... Additional arguments forwarded to
#'   \code{mvsusieR::mvsusie_rss}.
#' @return The fit object returned by \code{mvsusieR::mvsusie_rss}.
#' @export
fitMvsusieRss <- function(Z, R, N, prior_variance, coverage = 0.95, ...) {
  mvsusieR::mvsusie_rss(Z = Z, R = R, N = N,
                         prior_variance = prior_variance,
                         coverage = coverage, ...)
}

#' Fit fSuSiE on individual-level (X, Y, pos) data
#'
#' Thin wrapper around \code{fsusieR::susiF}.
#'
#' @param X Numeric matrix of genotypes (samples x variants).
#' @param Y Numeric matrix of multi-trait outcomes (samples x traits).
#' @param pos Numeric vector of trait positions (length \code{ncol(Y)}).
#' @param ... Additional arguments forwarded to \code{fsusieR::susiF}.
#' @return The fit object returned by \code{fsusieR::susiF}.
#' @export
fitFsusie <- function(X, Y, pos, ...) {
  fsusieR::susiF(X = X, Y = Y, pos = pos, ...)
}
