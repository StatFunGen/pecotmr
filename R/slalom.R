#' Slalom Function for Summary Statistics QC for Fine-Mapping Analysis
#'
#' Performs Approximate Bayesian Factor (ABF) analysis, identifies credible sets,
#' and annotates lead variants based on fine-mapping results. It computes p-values
#' from z-scores assuming a two-sided standard normal distribution.
#'
#' Provide either an LD correlation matrix \code{R} or a genotype matrix \code{X}
#' (from which LD is derived automatically via \code{compute_LD}).
#'
#' @param zScore Numeric vector of z-scores corresponding to each variant.
#' @param R Square LD correlation matrix. Provide either \code{R} or \code{X}.
#' @param X Genotype matrix (samples x SNPs). If provided, LD is computed via
#'   \code{compute_LD(X)}.
#' @param standard_error Optional numeric vector of standard errors corresponding
#'   to each z-score. If not provided, a default value of 1 is assumed for all variants.
#' @param abf_prior_variance Numeric, the prior effect size variance for ABF calculations.
#'   Default is 0.04.
#' @param nlog10p_dentist_s_threshold Numeric, the -log10 DENTIST-S P value threshold
#'   for identifying outlier variants for prediction. Default is 4.0.
#' @param r2_threshold Numeric, the r2 threshold for DENTIST-S outlier variants
#'   for prediction. Default is 0.6.
#' @param lead_variant_choice Character, method to choose the lead variant, either
#'   "pvalue" or "abf", with default "pvalue".
#' @param ld_method Character string specifying the LD computation method when
#'   \code{X} is provided. Passed to \code{\link{compute_LD}}. One of
#'   \code{"sample"} (default), \code{"population"}, or \code{"gcta"}.
#'   Ignored when \code{R} is provided directly.
#' @return A list containing the annotated LD matrix with ABF results, credible sets,
#'   lead variant, and DENTIST-S statistics; and a summary dataframe with aggregate statistics.
#' @examples
#' results <- slalom(zScore, R = R, standard_error = standard_error)
#' @seealso \code{\link{dentist_single_window}}, \code{\link{resolve_LD_input}}
#' @export
#'
slalom <- function(zScore, R = NULL, X = NULL, standard_error = rep(1, length(zScore)),
                   abf_prior_variance = 0.04, nlog10p_dentist_s_threshold = 4.0,
                   r2_threshold = 0.6, lead_variant_choice = "pvalue",
                   ld_method = "sample") {
  # Resolve LD matrix from R or X
  ld_resolved <- resolve_LD_input(R = R, X = X, need_nSample = FALSE,
                                  ld_method = ld_method)
  LD_mat <- ld_resolved$R

  if (!is.matrix(LD_mat) || nrow(LD_mat) != ncol(LD_mat) || nrow(LD_mat) != length(zScore)) {
    stop("LD_mat must be a square matrix matching the length of zScore.")
  }

  # One-sided p-value matching the original Python implementation (stats.norm.cdf).
  # This selects the most negative z-score as lead when lead_variant_choice == "pvalue".
  pvalue <- pnorm(zScore)

  logSumExp <- function(x) {
    max_x <- max(x, na.rm = TRUE)
    sum_exp <- sum(exp(x - max_x), na.rm = TRUE)
    return(max_x + log(sum_exp))
  }

  abf <- function(z, se, W = 0.04) {
    V <- se^2
    r <- W / (W + V)
    lbf <- 0.5 * (log(1 - r) + (r * z^2))
    denom <- logSumExp(lbf)
    prob <- exp(lbf - denom)
    return(list(lbf = lbf, prob = prob))
  }

  abf_results <- abf(zScore, standard_error, W = abf_prior_variance)
  lbf <- abf_results$lbf
  prob <- abf_results$prob

  get_cs <- function(prob, coverage = 0.95) {
    ordering <- order(prob, decreasing = TRUE)
    cumprob <- cumsum(prob[ordering])
    idx <- which(cumprob > coverage)[1]
    cs <- ordering[1:idx]
    return(cs)
  }

  cs <- get_cs(prob, coverage = 0.95)
  cs_99 <- get_cs(prob, coverage = 0.99)

  lead_idx <- if (lead_variant_choice == "pvalue") {
    which.min(pvalue)
  } else {
    which.max(prob)
  }

  r2 <- LD_mat^2
  t_dentist_s <- (zScore - LD_mat[, lead_idx] * zScore[lead_idx])^2 / (1 - r2[, lead_idx])
  t_dentist_s[t_dentist_s < 0] <- Inf
  nlog10p_dentist_s <- -log10(pchisq(t_dentist_s, df = 1, lower.tail = FALSE))
  outliers <- (r2[, lead_idx] > r2_threshold) & (nlog10p_dentist_s > nlog10p_dentist_s_threshold)

  n_r2 <- sum(r2[, lead_idx] > r2_threshold)
  n_dentist_s_outlier <- sum(outliers, na.rm = TRUE)
  max_pip <- max(prob)

  summary <- list(
    lead_pip_variant = lead_idx,
    n_total = length(zScore),
    n_r2 = n_r2,
    n_dentist_s_outlier = n_dentist_s_outlier,
    fraction = ifelse(n_r2 > 0, n_dentist_s_outlier / n_r2, 0),
    max_pip = max_pip,
    cs_95 = cs,
    cs_99 = cs_99
  )
  result <- as.data.frame(list(original_z = zScore, prob = prob, pvalue = pvalue, outliers = outliers, nlog10p_dentist_s = nlog10p_dentist_s))

  return(list(data = result, summary = summary))
}
