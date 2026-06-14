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
    csIndices <- unlist(lCs[[k]]) # Extract indices for the current credible set
    # Calculate purity based on the specified method
    if (method == "min") {
      if (length(csIndices) == 1) {
        tt[[k]] <- 1 # Set purity to 1 for non-"min" methods
      } else {
        x <- abs(cor(X[, csIndices])) # Compute the absolute correlation matrix
        x[col(x) == row(x)] <- NA # Set diagonal elements to NA to exclude them
        tt[[k]] <- min(x, na.rm = TRUE) # Calculate minimum off-diagonal correlation for "min" method
        # Check if the credible set has only one element and the method is not "min"
      }
    } else {
      if (length(csIndices) == 1) {
        tt[[k]] <- c(1, 1, 1) # Set purity to 1 for non-"min" methods
      } else {
        x <- abs(cor(X[, csIndices])) # Compute the absolute correlation matrix
        x[col(x) == row(x)] <- NA # Set diagonal elements to NA to exclude them
        # Calculate min, mean, and median of off-diagonal correlations for other methods
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

