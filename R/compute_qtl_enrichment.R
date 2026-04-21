#' @title Implementation of enrichment analysis described in https://doi.org/10.1371/journal.pgen.1006646
#'
#' @description Largely follows from fastenloc https://github.com/xqwen/fastenloc
#' but uses `susieR` fitted objects as input to estimate prior for use with `coloc` package (coloc v5, aka SuSiE-coloc).
#' The main differences are 1) now enrichment is based on all QTL variants whether or not they are inside signal clusters;
#' 2) Causal QTL are sampled from SuSiE single effects, not signal clusters;
#' 3) Allow a variant to be QTL for not only multiple conditions (eg cell types) but also multiple regions (eg genes).
#' Other minor improvements include 1) Make GSL RNG thread-safe; 2) Release memory from QTL binary annotation samples immediately after they are used.
#' @details Uses output of \code{\link[susieR]{susie}} from the
#'   \code{susieR} package.
#'
#' @param gwas_pip This is a vector of GWAS PIP, genome-wide.
#' @param susie_qtl_regions This is a list of SuSiE fitted objects per QTL unit analyzed
#' @param num_gwas This parameter is highly important if GWAS input does not contain all SNPs interrogated (e.g., in some cases, only fine-mapped geomic regions are included).
#' Then users must pick a value of total_variants and estimate pi_gwas beforehand by: sum(gwas_pip$pip)/num_gwas. If num_gwas is null, pi_gwas would be sum(gwas_pip$pip)/total_variants.
#' @param pi_qtl This parameter can be safely left to default if your input QTL data has enough regions to estimate it.
#' @param lambda Similar to the shrinkage parameter used in ridge regression. It takes any non-negative value and shrinks the enrichment estimate towards 0.
#' When it is set to 0, no shrinkage will be applied. A large value indicates strong shrinkage. The default value is set to 1.0.
#' @param ImpN Rounds of multiple imputation to draw QTL from, default is 25.
#' @param num_threads Number of Simultaneous running CPU threads for multiple imputation, default is 1.
#' @return A list of enrichment parameter estimates
#'
#' @examples
#'
#' # Simulate fake data for gwas_pip
#' n_gwas_pip <- 1000
#' gwas_pip <- runif(n_gwas_pip)
#' names(gwas_pip) <- paste0("snp", 1:n_gwas_pip)
#' gwas_fit <- list(pip = gwas_pip)
#' # Simulate fake data for a single SuSiEFit object
#' simulate_susiefit <- function(n, p) {
#'   pip <- runif(n)
#'   names(pip) <- paste0("snp", 1:n)
#'   alpha <- t(matrix(runif(n * p), nrow = n))
#'   alpha <- t(apply(alpha, 1, function(row) row / sum(row)))
#'   list(
#'     pip = pip,
#'     alpha = alpha,
#'     prior_variance = runif(p)
#'   )
#' }
#' # Simulate multiple SuSiEFit objects
#' n_susie_fits <- 2
#' susie_fits <- replicate(n_susie_fits, simulate_susiefit(n_gwas_pip, 10), simplify = FALSE)
#' # Add these fits to a list, providing names to each element
#' names(susie_fits) <- paste0("fit", 1:length(susie_fits))
#' # Set other parameters
#' ImpN <- 10
#' lambda <- 1
#' num_threads <- 1
#' library(pecotmr)
#' en <- compute_qtl_enrichment(gwas_fit, susie_fits, lambda = lambda, ImpN = ImpN, num_threads = num_threads)
#'
#' @seealso \code{\link[susieR]{susie}}
#' @useDynLib pecotmr, .registration = TRUE
#' @export
#'
compute_qtl_enrichment <- function(gwas_pip, susie_qtl_regions,
                                   num_gwas = NULL, pi_qtl = NULL,
                                   lambda = 1.0, ImpN = 25,
                                   double_shrinkage = FALSE,
                                   bessel_correction = TRUE,
                                   num_threads = 1, verbose = TRUE) {
  if (is.null(num_gwas)) {
    warning("num_gwas is not provided. Estimating pi_gwas from the data. Note that this estimate may be biased if the input gwas_pip does not contain genome-wide variants.")
    pi_gwas <- sum(gwas_pip) / length(gwas_pip)
    if (verbose) {
      message(paste("Estimated pi_gwas: ", round(pi_gwas, 5), "\n"))
    }
  } else {
    pi_gwas <- sum(gwas_pip) / num_gwas
  }

  if (is.null(pi_qtl)) {
    warning("pi_qtl is not provided. Estimating pi_qtl from the data. Note that this estimate may be biased if either 1) the input susie_qtl_regions does not have enough data, or 2) the single effects only include variables inside of credible sets or signal clusters.")
    num_signal <- 0
    num_test <- 0
    for (d in susie_qtl_regions) {
      num_signal <- num_signal + sum(d$pip)
      num_test <- num_test + length(d$pip)
    }
    pi_qtl <- num_signal / num_test
    if (verbose) {
      message(paste("Estimated pi_qtl: ", round(pi_qtl, 5), "\n"))
    }
  }

  if (pi_gwas == 0) stop("Cannot perform enrichment analysis. No association signal found in GWAS data.")
  if (pi_qtl == 0) stop("Cannot perform enrichment analysis. No QTL associated with the molecular phenotype.")

  # Check if names of gwas_pip and susie_qtl_regions$pip are both available
  if (is.null(names(gwas_pip))) {
    stop("Variant names are missing in gwas_pip. Please provide named gwas_pip data.")
  }
  if (!all(sapply(susie_qtl_regions, function(x) !is.null(names(x$pip))))) {
    stop("Variant names are missing in susie_qtl_regions$pip. Please provide susie_qtl_regions with named pip data.")
  }

  # Align the names of susie_qtl_regions$pip to gwas_pip names and document unmatched variants
  aligned_susie_qtl_regions <- lapply(susie_qtl_regions, function(x) {
    alignment_result <- align_variant_names(names(x$pip), names(gwas_pip))
    names(x$pip) <- alignment_result$aligned_variants
    if (length(alignment_result$unmatched_indices) > 0) {
      x$unmatched_variants <- names(x$pip)[alignment_result$unmatched_indices]
    }
    x
  })
  unmatched_variants <- lapply(aligned_susie_qtl_regions, function(x) x$unmatched_variants)

  # Update susie_qtl_regions with the aligned variant names
  susie_qtl_regions <- lapply(aligned_susie_qtl_regions, function(x) {
    x$unmatched_variants <- NULL
    x
  })

  # cpp11 requires exact integer types for int parameters
  en <- qtl_enrichment_rcpp(
    r_gwas_pip = gwas_pip,
    r_qtl_susie_fit = susie_qtl_regions,
    pi_gwas = pi_gwas,
    pi_qtl = pi_qtl,
    ImpN = as.integer(ImpN),
    shrinkage_lambda = lambda,
    double_shrinkage = double_shrinkage,
    bessel_correction = bessel_correction,
    num_threads = as.integer(num_threads)
  )

  # Add the unmatched variants to the output
  en <- list(en)
  en$unused_xqtl_variants <- unmatched_variants

  return(en)
}
