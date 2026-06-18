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
#' @param gwasPip This is a vector of GWAS PIP, genome-wide.
#' @param susieQtlRegions This is a list of SuSiE fitted objects per QTL unit analyzed
#' @param numGwas This parameter is highly important if GWAS input does not contain all SNPs interrogated (e.g., in some cases, only fine-mapped geomic regions are included).
#' Then users must pick a value of total_variants and estimate piGwas beforehand by: sum(gwasPip$pip)/numGwas. If numGwas is null, piGwas would be sum(gwasPip$pip)/total_variants.
#' @param piQtl This parameter can be safely left to default if your input QTL data has enough regions to estimate it.
#' @param lambda Similar to the shrinkage parameter used in ridge regression. It takes any non-negative value and shrinks the enrichment estimate towards 0.
#' When it is set to 0, no shrinkage will be applied. A large value indicates strong shrinkage. The default value is set to 1.0.
#' @param impN Rounds of multiple imputation to draw QTL from, default is 25.
#' @param numThreads Number of Simultaneous running CPU threads for multiple imputation, default is 1.
#' @return A list of enrichment parameter estimates
#'
#' @examples
#'
#' # Simulate fake data for gwasPip
#' nGwasPip <- 1000
#' gwasPip <- runif(nGwasPip)
#' names(gwasPip) <- paste0("snp", 1:nGwasPip)
#' gwasFit <- list(pip = gwasPip)
#' # Simulate fake data for a single SuSiEFit object
#' simulateSusiefit <- function(n, p) {
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
#' nSusieFits <- 2
#' susieFits <- replicate(nSusieFits, simulateSusiefit(nGwasPip, 10), simplify = FALSE)
#' # Add these fits to a list, providing names to each element
#' names(susieFits) <- paste0("fit", 1:length(susieFits))
#' # Set other parameters
#' impN <- 10
#' lambda <- 1
#' numThreads <- 1
#' library(pecotmr)
#' en <- computeQtlEnrichment(gwasFit, susieFits, lambda = lambda, impN = impN, numThreads = numThreads)
#'
#' @seealso \code{\link[susieR]{susie}}
#' @useDynLib pecotmr, .registration = TRUE
#' @export
#'
computeQtlEnrichment <- function(gwasPip, susieQtlRegions,
                                 numGwas = NULL, piQtl = NULL,
                                 lambda = 1.0, impN = 25,
                                 doubleShrinkage = FALSE,
                                 besselCorrection = TRUE,
                                 numThreads = 1, verbose = TRUE) {
  if (is.null(numGwas)) {
    warning("numGwas is not provided. Estimating piGwas from the data. Note that this estimate may be biased if the input gwasPip does not contain genome-wide variants.")
    piGwas <- sum(gwasPip) / length(gwasPip)
    if (verbose) {
      message(paste("Estimated piGwas: ", round(piGwas, 5), "\n"))
    }
  } else {
    piGwas <- sum(gwasPip) / numGwas
  }

  if (is.null(piQtl)) {
    warning("piQtl is not provided. Estimating piQtl from the data. Note that this estimate may be biased if either 1) the input susieQtlRegions does not have enough data, or 2) the single effects only include variables inside of credible sets or signal clusters.")
    numSignal <- 0
    numTest <- 0
    for (d in susieQtlRegions) {
      numSignal <- numSignal + sum(d$pip)
      numTest <- numTest + length(d$pip)
    }
    piQtl <- numSignal / numTest
    if (verbose) {
      message(paste("Estimated piQtl: ", round(piQtl, 5), "\n"))
    }
  }

  if (piGwas == 0) stop("Cannot perform enrichment analysis. No association signal found in GWAS data.")
  if (piQtl == 0) stop("Cannot perform enrichment analysis. No QTL associated with the molecular phenotype.")

  # Check if names of gwasPip and susieQtlRegions$pip are both available
  if (is.null(names(gwasPip))) {
    stop("Variant names are missing in gwasPip. Please provide named gwasPip data.")
  }
  if (!all(sapply(susieQtlRegions, function(x) !is.null(names(x$pip))))) {
    stop("Variant names are missing in susieQtlRegions$pip. Please provide susieQtlRegions with named pip data.")
  }

  # Align the names of susieQtlRegions$pip to gwasPip names and document unmatched variants
  alignedSusieQtlRegions <- lapply(susieQtlRegions, function(x) {
    alignmentResult <- alignVariantNames(names(x$pip), names(gwasPip))
    names(x$pip) <- alignmentResult$alignedVariants
    if (length(alignmentResult$unmatchedIndices) > 0) {
      x$unmatched_variants <- names(x$pip)[alignmentResult$unmatchedIndices]
    }
    x
  })
  unmatchedVariants <- lapply(alignedSusieQtlRegions, function(x) x$unmatched_variants)

  # Update susieQtlRegions with the aligned variant names
  susieQtlRegions <- lapply(alignedSusieQtlRegions, function(x) {
    x$unmatched_variants <- NULL
    x
  })

  # cpp11 requires exact integer types for int parameters
  en <- qtlEnrichmentRcpp(
    rGwasPip = gwasPip,
    rQtlSusieFit = susieQtlRegions,
    piGwas = piGwas,
    piQtl = piQtl,
    ImpN = as.integer(impN),
    shrinkageLambda = lambda,
    doubleShrinkage = doubleShrinkage,
    besselCorrection = besselCorrection,
    numThreads = as.integer(numThreads)
  )

  # Add the unmatched variants to the output
  en <- list(en)
  en$unused_xqtl_variants <- unmatchedVariants

  return(en)
}

