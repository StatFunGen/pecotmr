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

  logSumExp <- function(x) {
    maxX <- max(x, na.rm = TRUE)
    sumExp <- sum(exp(x - maxX), na.rm = TRUE)
    return(maxX + log(sumExp))
  }

  abf <- function(z, se, W = 0.04) {
    V <- se^2
    r <- W / (W + V)
    lbf <- 0.5 * (log(1 - r) + (r * z^2))
    denom <- logSumExp(lbf)
    prob <- exp(lbf - denom)
    return(list(lbf = lbf, prob = prob))
  }

  abfResults <- abf(zScore, standardError, W = abfPriorVariance)
  lbf <- abfResults$lbf
  prob <- abfResults$prob

  getCs <- function(prob, coverage = 0.95) {
    ordering <- order(prob, decreasing = TRUE)
    cumprob <- cumsum(prob[ordering])
    idx <- which(cumprob > coverage)[1]
    cs <- ordering[1:idx]
    return(cs)
  }

  cs <- getCs(prob, coverage = 0.95)
  cs99 <- getCs(prob, coverage = 0.99)

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
    lead_pip_variant = leadIdx,
    n_total = length(zScore),
    n_r2 = nR2,
    n_dentist_s_outlier = nDentistSOutlier,
    fraction = ifelse(nR2 > 0, nDentistSOutlier / nR2, 0),
    max_pip = maxPip,
    cs_95 = cs,
    cs_99 = cs99
  )
  result <- as.data.frame(list(original_z = zScore, prob = prob, pvalue = pvalue, outliers = outliers, nlog10p_dentist_s = nlog10pDentistS))

  return(list(data = result, summary = summary))
}
