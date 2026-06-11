computeQvalues <- function(pvalues) {
  # Make sure qvalue is installed
  if (!requireNamespace("qvalue", quietly = TRUE)) {
    stop("To use this function, please install qvalue: https://www.bioconductor.org/packages/release/bioc/html/qvalue.html")
  }
  if (all(is.na(pvalues))) {
    message("All p-values are NA. Returning NA vector.")
    return(rep(NA_real_, length(pvalues)))
  }      
  tryCatch(
    {
      if (length(pvalues) < 2) {
        return(pvalues)
      } else {
        return(qvalue::qvalue(pvalues)$qvalues)
      }
    },
    error = function(e) {
      message("Too few p-values to calculate qvalue, fall back to BH")
      qvalue::qvalue(pvalues, pi0 = 1)$qvalues
    }
  )
}

matxMax <- function(mtx) {
  return(arrayInd(which.max(mtx), dim(mtx)))
}

computeMaf <- function(geno) {
  f <- mean(geno, na.rm = TRUE) / 2
  return(min(f, 1 - f))
}

#' Derive minor-allele frequency from effect-allele frequency
#'
#' MAF is an internal QC/filtering quantity only; it is never exported. Use this
#' helper wherever a MAF is needed from a (directional) effect-allele frequency
#' \code{af}, instead of carrying a separate \code{maf} column. NA in -> NA out.
#'
#' @param af Numeric vector of effect-allele frequencies in \code{[0, 1]}.
#' @return Numeric vector \code{pmin(af, 1 - af)}, preserving NA.
#' @noRd
mafFromAf <- function(af) {
  af <- as.numeric(af)
  pmin(af, 1 - af)
}

computeMissing <- function(geno) {
  miss <- sum(is.na(geno)) / length(geno)
  return(miss)
}

computeNonMissingY <- function(y) {
  nonmiss <- sum(!is.na(y))
  return(nonmiss)
}

computeAllMissingY <- function(y) {
  allmiss <- all(is.na(y))
  return(allmiss)
}

meanImpute <- function(geno) {
  f <- apply(geno, 2, function(x) mean(x, na.rm = TRUE))
  for (i in seq_along(f)) geno[, i][which(is.na(geno[, i]))] <- f[i]
  return(geno)
}

isZeroVariance <- function(x) length(unique(x)) == 1

#' Safe truncated SVD with numerical stability
#'
#' Computes a thin SVD and optionally truncates small singular values.
#' Useful for avoiding numerical issues when working with rank-deficient
#' or near-singular matrices.
#'
#' @param mat Input matrix (n x p).
#' @param tol Relative tolerance for filtering singular values.
#'   Singular values smaller than \code{tol * max(d)} are discarded.
#'   Set to 0 to keep all singular values.
#' @param maxRank Optional maximum number of singular values to retain.
#'   If NULL, all singular values passing the tolerance filter are kept.
#' @return A list with components:
#'   \describe{
#'     \item{u}{Left singular vectors (n x r matrix).}
#'     \item{d}{Singular values (length-r numeric vector).}
#'     \item{v}{Right singular vectors (p x r matrix).}
#'   }
#'   where r is the number of retained singular values.
#' @noRd
safeSvd <- function(mat, tol = 1e-8, maxRank = NULL) {
  if (max(abs(mat)) == 0) {
    stop("Cannot compute SVD of an all-zero matrix.")
  }
  # Compute thin SVD
  s <- svd(mat)
  d <- s$d
  # Filter by relative tolerance
  if (tol > 0 && length(d) > 0) {
    keep <- d / d[1] > tol
    if (!any(keep)) {
      stop("All singular values are below the tolerance threshold.")
    }
  } else {
    keep <- rep(TRUE, length(d))
  }
  # Apply maxRank cap
  if (!is.null(maxRank) && maxRank > 0) {
    nKeep <- min(sum(keep), maxRank)
    keepIdx <- which(keep)
    if (length(keepIdx) > nKeep) {
      keep[keepIdx[(nKeep + 1):length(keepIdx)]] <- FALSE
    }
  }
  r <- sum(keep)
  list(
    u = s$u[, keep, drop = FALSE],
    d = d[keep],
    v = s$v[, keep, drop = FALSE]
  )
}

#' Compute LD (Linkage Disequilibrium) Correlation Matrix from Genotypes
#'
#' Computes a pairwise Pearson correlation matrix from a genotype matrix.
#' Supports three variance conventions:
#' \describe{
#'   \item{\code{"sample"}}{Standard sample variance with N-1 denominator (default).
#'     Uses mean imputation for missing genotypes, then \code{Rfast::cora} (if available)
#'     or base \code{cor()}.}
#'   \item{\code{"population"}}{Population variance with N denominator, matching
#'     GCTA-style tools (e.g. DENTIST, GCTA --make-grm). Per-SNP means are computed
#'     from non-missing values; missing entries are set to zero after centering so they
#'     do not contribute to cross-products. Cross-products are normalized by the total
#'     sample count N, not by pairwise non-missing counts.}
#'   \item{\code{"gcta"}}{GCTA per-pair missing data correction. Like \code{"population"}
#'     but applies a correction term for each SNP pair based on the number of jointly
#'     non-missing samples. Matches the exact formula from the DENTIST C++ binary's
#'     \code{calcLDFromBfile_gcta}. Use this when missingness varies substantially
#'     across SNPs and accuracy of individual LD entries matters.}
#' }
#'
#' @param X Numeric genotype matrix (samples x SNPs). May contain \code{NA}
#'   for missing genotypes.
#' @param method Character, one of \code{"sample"} (default, N-1 denominator),
#'   \code{"population"} (N denominator, GCTA-style), or \code{"gcta"} (per-pair
#'   missing data correction). Partial matching is supported.
#' @param backend Character, one of \code{"internal"} (default), \code{"snprelate"},
#'   or \code{"snpstats"}. Controls which library computes the correlation matrix
#'   when \code{method = "sample"}:
#'   \describe{
#'     \item{\code{"internal"}}{Uses \code{Rfast::cora} if available, otherwise
#'       base \code{cor()}.}
#'     \item{\code{"snprelate"}}{Requires a temporary GDS file; uses
#'       \code{SNPRelate::snpgdsLDMat(method = "corr")}.}
#'     \item{\code{"snpstats"}}{Converts to \code{SnpMatrix}; uses
#'       \code{snpStats::ld(, stat = "R")}.}
#'   }
#'   The \code{"snprelate"} and \code{"snpstats"} backends are only supported
#'   with \code{method = "sample"}; combining them with other methods will
#'   raise an error.
#' @param trimSamples Logical. If \code{TRUE} and \code{method} is
#'   \code{"population"} or \code{"gcta"}, drops trailing samples so that
#'   \code{nrow(X)} is a multiple of 4, matching PLINK .bed file chunk processing.
#'   Ignored when \code{method = "sample"}. Default is \code{FALSE}.
#' @param shrinkage Numeric in (0, 1]. Shrink the LD matrix toward the identity:
#'   \code{R_s = (1 - shrinkage) * R + shrinkage * I}. Useful for regularizing
#'   LD for summary-statistics-based methods such as lassosum (Mak et al 2017).
#'   Default is 0 (no shrinkage).
#'
#' @return A symmetric correlation matrix with row and column names taken from
#'   \code{colnames(X)}.
#'
#' @details
#' \strong{Missing data handling.}
#' With \code{method = "sample"}, missing values are mean-imputed per SNP
#' before computing the full Pearson correlation matrix.
#' With \code{method = "population"}, per-SNP means are computed from
#' non-missing values, the matrix is centered, then \code{NA}s are set to 0
#' so that missing pairs contribute nothing to the cross-product.
#' The denominator is always the total sample count \code{N}
#' (after optional trimming), matching the original GCTA formula:
#' \deqn{\text{Var}(X_i) = E[X_i^2] - E[X_i]^2}
#' \deqn{\text{Cor}(X_i, X_j) = \frac{\text{Cov}(X_i, X_j)}{\sqrt{\text{Var}(X_i)\,\text{Var}(X_j)}}}
#'
#' \strong{Zero-variance SNPs.}
#' Any monomorphic SNP will have zero variance, producing \code{NaN}
#' correlations. These are set to 0 in the returned matrix; the diagonal
#' is forced to 1.
#'
#' @examples
#' \dontrun{
#' X <- matrix(sample(0:2, 500, replace = TRUE), nrow = 50)
#' colnames(X) <- paste0("rs", 1:10)
#'
#' # Standard sample correlation (default)
#' R1 <- computeLd(X)
#'
#' # GCTA-style population variance
#' R2 <- computeLd(X, method = "population")
#'
#' # GCTA-style with per-pair missing data correction
#' R3 <- computeLd(X, method = "gcta")
#' }
#'
#' @export
computeLd <- function(X, method = c("sample", "population", "gcta"),
                      backend = c("internal", "snprelate", "snpstats"),
                      trimSamples = FALSE, shrinkage = 0) {
  if (is.null(X)) {
    stop("X must be provided.")
  }
  method <- match.arg(method)
  backend <- match.arg(backend)
  nms <- colnames(X)

  if (method == "sample") {
    # ---- Standard sample correlation (N-1 denominator) ----
    if (backend == "snprelate") {
      R <- .computeLdSnprelate(X)
    } else if (backend == "snpstats") {
      R <- .computeLdSnpstats(X)
    } else {
      # internal backend: Rfast::cora if available, else base cor()
      # Mean impute only if NAs exist (PLINK2 data typically has none)
      X_imp <- X
      if (anyNA(X_imp)) {
        colMeansX <- colMeans(X_imp, na.rm = TRUE)
        naPos <- which(is.na(X_imp), arr.ind = TRUE)
        X_imp[naPos] <- colMeansX[naPos[, 2]]
      }
      if (requireNamespace("Rfast", quietly = TRUE)) {
        # large=FALSE uses tcrossprod internally, ~40x faster than large=TRUE
        R <- Rfast::cora(X_imp, large = FALSE)
      } else {
        R <- cor(X_imp)
      }
    }
  } else if (method == "population") {
    if (backend != "internal") {
      stop("backend '", backend, "' is only supported with method='sample'.")
    }
    # ---- Population variance (N denominator, GCTA-style) ----
    # Optionally trim trailing samples to a multiple of 4 (matches .bed processing)
    if (trimSamples) {
      N_kept <- (nrow(X) %/% 4L) * 4L
      if (N_kept < nrow(X)) X <- X[seq_len(N_kept), , drop = FALSE]
    }
    N <- nrow(X)
    # Per-SNP means from non-missing values
    colMeansX <- colMeans(X, na.rm = TRUE)
    # Population variance: E[X^2] - E[X]^2
    colVarsX <- colMeans(X^2, na.rm = TRUE) - colMeansX^2
    # Center; set NA -> 0 so missing pairs don't contribute to cross-products.
    # NOTE: the covariance divides by total N (not pairwise non-missing counts),
    # which is an approximation that assumes uniform missingness across SNPs.
    # With heterogeneous missingness, correlations between high-missing and
    # low-missing columns will be slightly deflated. This matches the GCTA
    # convention and is standard for PLINK-style LD computation.
    if (anyNA(X)) {
      naRates <- colMeans(is.na(X))
      if (max(naRates) - min(naRates) > 0.1) {
        warning("Population LD method with heterogeneous missingness ",
                "(max NA rate ", round(max(naRates), 3),
                ", min ", round(min(naRates), 3),
                "): correlations may be biased. Consider using method='sample' ",
                "which handles missingness via mean imputation.")
      }
    }
    X_c <- sweep(X, 2, colMeansX)
    X_c[is.na(X_c)] <- 0
    # Covariance with N denominator
    covMat <- crossprod(X_c) / N
    # Correlation
    sdVec <- sqrt(colVarsX)
    R <- covMat / outer(sdVec, sdVec)
  } else {
    if (backend != "internal") {
      stop("backend '", backend, "' is only supported with method='sample'.")
    }
    # ---- GCTA per-pair missing data correction ----
    # Matches the DENTIST binary's calcLDFromBfile_gcta formula exactly.
    # Unlike "population" which divides by total N, this method tracks
    # per-pair missing counts and applies a correction term.
    if (trimSamples) {
      N_kept <- (nrow(X) %/% 4L) * 4L
      if (N_kept < nrow(X)) X <- X[seq_len(N_kept), , drop = FALSE]
    }
    N <- nrow(X)
    p <- ncol(X)

    # Marginal statistics from non-missing values
    colMeansX <- colMeans(X, na.rm = TRUE)
    colMeanSq <- colMeans(X^2, na.rm = TRUE)
    colVarsX <- colMeanSq - colMeansX^2

    # Build indicator matrix for non-missing values
    notNa <- !is.na(X)
    # Replace NA with 0 for cross-product computation
    X_zero <- X
    X_zero[is.na(X_zero)] <- 0

    # Per-pair non-missing counts: notNa'notNa gives count of jointly observed
    pairCounts <- crossprod(notNa * 1.0)
    nMissing <- N - pairCounts

    # Per-pair sums: sum of X_i over samples where both i and j are observed
    # For the correction term we need E_i2 = sum_i_pair / N (pair-specific mean)
    # X_zero' %*% notNa gives, for each (i,j), sum of X_i where j is not missing
    pairSums <- crossprod(X_zero, notNa * 1.0)

    # Cross-product sum: sum(X_i * X_j) over jointly non-missing samples
    sum_XY <- crossprod(X_zero)

    # GCTA correction formula:
    # E_i2[i,j] = pairSums[i,j] / N  (mean of SNP i restricted to non-missing-j samples, divided by N)
    # cov = sum_XY/N + E[i]*E[j]*(N-m)/N - E[i]*E_j2 - E_i2*E[j]
    E_i2 <- pairSums / N  # p x p: row i, col j = sum of X_i where j non-missing, / N
    E_j2 <- t(E_i2)        # transposed version

    covMat <- sum_XY / N +
      outer(colMeansX, colMeansX) * (pairCounts / N) -
      colMeansX * E_j2 -
      E_i2 * rep(colMeansX, each = p)

    # Correlation
    sdVec <- sqrt(colVarsX)
    sdOuter <- outer(sdVec, sdVec)
    R <- matrix(0.001, p, p)
    valid <- sdOuter > 0
    R[valid] <- covMat[valid] / sdOuter[valid]
  }

  # Ensure clean output
  diag(R) <- 1.0
  R[is.na(R) | is.nan(R)] <- 0

  # Optional shrinkage toward identity: R_s = (1 - shrinkage) * R + shrinkage * I
  # Used e.g. by lassosum (Mak et al 2017) to regularize LD for RSS methods.
  if (shrinkage > 0 && shrinkage <= 1) {
    R <- (1 - shrinkage) * R + shrinkage * diag(nrow(R))
  }

  colnames(R) <- rownames(R) <- nms
  R
}

#' Compute LD via SNPRelate (creates a temporary GDS file from the dosage matrix).
#' @param X Numeric genotype matrix (samples x SNPs).
#' @return Correlation matrix.
#' @noRd
.computeLdSnprelate <- function(X) {
  if (!requireNamespace("SNPRelate", quietly = TRUE))
    stop("Package 'SNPRelate' is required for backend='snprelate'")
  if (!requireNamespace("gdsfmt", quietly = TRUE))
    stop("Package 'gdsfmt' is required for backend='snprelate'")

  tmpGds <- tempfile(fileext = ".gds")
  on.exit(unlink(tmpGds), add = TRUE)

  # Round to integer dosage for GDS (0/1/2)
  X_int <- round(X)
  storage.mode(X_int) <- "integer"
  X_int[is.na(X_int)] <- 3L  # GDS missing code

  snpIds <- colnames(X) %||% seq_len(ncol(X))
  sampleIds <- rownames(X) %||% seq_len(nrow(X))

  SNPRelate::snpgdsCreateGeno(tmpGds,
    genmat = X_int,
    sample.id = sampleIds,
    snp.id = snpIds,
    snp.chromosome = rep(1L, ncol(X)),
    snp.position = seq_len(ncol(X)),
    snpfirstdim = FALSE
  )

  gds <- SNPRelate::snpgdsOpen(tmpGds, readonly = TRUE)
  on.exit(SNPRelate::snpgdsClose(gds), add = TRUE)

  ldObj <- SNPRelate::snpgdsLDMat(gds, method = "corr",
                                  slide = -1, verbose = FALSE)
  ldObj$LD
}

#' Compute LD via snpStats (converts dosage matrix to SnpMatrix).
#' @param X Numeric genotype matrix (samples x SNPs).
#' @return Correlation matrix (r, not r²).
#' @noRd
.computeLdSnpstats <- function(X) {
  if (!requireNamespace("snpStats", quietly = TRUE))
    stop("Package 'snpStats' is required for backend='snpstats'")

  # snpStats expects counts of the B allele as raw codes: 1=AA, 2=AB, 3=BB, 0=NA
  # pecotmr dosage is ALT count (0/1/2), so map: 0->1, 1->2, 2->3, NA->0
  X_raw <- round(X) + 1L
  X_raw[is.na(X) | X_raw < 1L] <- 0L
  X_raw[X_raw > 3L] <- 3L
  storage.mode(X_raw) <- "raw"
  sm <- new("SnpMatrix", X_raw)

  R <- as.matrix(snpStats::ld(sm, stats = "R", depth = ncol(X) - 1L))
  # snpStats::ld returns a sparse-like matrix; ensure full dense
  R[is.na(R)] <- 0
  diag(R) <- 1
  R
}

#' @importFrom matrixStats colVars
filterX <- function(X, missingRateThresh, mafThresh, varThresh = 0, maf = NULL, xVariance = NULL) {
  totalVariants <- ncol(X)
  if (!is.null(missingRateThresh) && missingRateThresh < 1.0) {
    rmCol <- which(apply(X, 2, computeMissing) > missingRateThresh)
    if (length(rmCol)) X <- X[, -rmCol, drop = FALSE]
  }

  # Check if non-NA values are valid genotypes before MAF filtering
  if (!is.null(mafThresh) && mafThresh > 0.0) {
    validGenotypes <- all(sapply(1:ncol(X), function(i) {
      x <- X[!is.na(X[, i]), i]
      all(x %in% c(0, 1, 2))
    }))

    if (validGenotypes || !is.null(maf)) {
      rmCol <- if (!is.null(maf)) which(maf <= mafThresh) else which(apply(X, 2, computeMaf) <= mafThresh)
      if (length(rmCol)) X <- X[, -rmCol, drop = FALSE]
    } else {
      message("Skipping MAF filtering as X does not appear to be 0/1/2 matrix, and no external MAF information is provided")
    }
  }

  rmCol <- which(apply(X, 2, isZeroVariance))
  if (length(rmCol)) X <- X[, -rmCol, drop = FALSE]
  X <- meanImpute(X)
  if (varThresh > 0) {
    rmCol <- if (!is.null(xVariance)) which(xVariance < varThresh) else which(colVars(X) < varThresh)
    if (length(rmCol)) X <- X[, -rmCol, drop = FALSE]
  }
  nDropped <- totalVariants - ncol(X)
  if (nDropped > 0) {
    message(paste0(nDropped, " out of ", totalVariants, " total variants dropped due to quality control on X matrix."))
  }
  return(X)
}

#' This function performing filters on X variants based on Y subjects for TWAS analysis. This function checks
#' whether the absence (NA) of certain subjects would lead to monomorphic in some variants in X after removing
#' of these subjects data from X.
#' @param missingRateThresh Maximum individual missingness cutoff.
#' @param mafThresh Minimum minor allele frequency (MAF) cutoff.
#' @param varThresh Minimum variance cutoff for a variant. Default is 0.
#' @param xVariance A vector of variance for X variants.
filterXWithY <- function(X, Y, missingRateThresh, mafThresh, varThresh = 0, maf = NULL, xVariance = NULL) {
  totalVariants <- ncol(X)
  X <- filterX(X, missingRateThresh, mafThresh, varThresh = varThresh, maf = maf, xVariance = xVariance)
  dropIdx <- do.call(c, lapply(colnames(Y), function(context) {
    subjectsWithNaY <- rownames(Y)[is.na(Y[, context])]
    X_temp <- X
    X_temp[subjectsWithNaY, ] <- NA
    rmCol <- which(apply(X_temp, 2, function(x) isZeroVariance(na.omit(x))))
    return(unique(rmCol))
  }))
  dropIdx <- unique(sort(dropIdx))
  if (length(dropIdx)) X <- X[, -dropIdx, drop = FALSE]
  if (length(dropIdx) > 0) {
    message(paste0("Additional ", length(dropIdx), " variants dropped after considering missing data in Y matrix, with ", ncol(X), " variants left."))
  }
  return(X)
}

filterY <- function(Y, nNonmiss) {
  rmCol <- which(apply(Y, 2, computeNonMissingY) < nNonmiss)
  if (length(rmCol)) Y <- Y[, -rmCol]
  rmRows <- NULL
  if (is.matrix(Y)) {
    rmRows <- which(apply(Y, 1, computeAllMissingY))
    if (length(rmRows)) Y <- Y[-rmRows, ]
  } else {
    Y <- Y[which(!is.na(Y))]
  }
  return(list(Y = Y, rmRows = rmRows))
}

# Retrieve a nested element from a list structure
#' @export
getNestedElement <- function(nestedList, nameVector) {
  if (is.null(nameVector)) {
    return(NULL)
  }
  nameVector <- nameVector[nameVector!='']
  currentElement <- nestedList
  for (name in nameVector) {
    if (is.null(currentElement[[name]])) {
      stop("Element not found in the list")
    }
    currentElement <- currentElement[[name]]
  }
  return(currentElement)
}



#' Utility function to specify the path to access the target list item in a nested list, especially when some list layers
#' in between are dynamic or uncertain.
#' @export
findData <- function(x, depthObj, showPath = FALSE, rmNull = TRUE, rmDup = FALSE, docall = c, lastObj = NULL) {
  depth <- as.integer(depthObj[1])
  listName <- if (length(depthObj) > 1) depthObj[2:length(depthObj)] else NULL
  if (depth == 1 || depth == 0) {
    if (!is.null(listName)) {
      if (listName[1] %in% names(x)) {
        if (any(grepl("^[0-9]+$", listName))) { # list names, indx name, list names
          secondDepth <- which(grepl("^[0-9]+$", listName))[1]
          data <- getNestedElement(x, listName[1:secondDepth[1] - 1])
          remainingPath <- listName[secondDepth:length(listName)]
          return(findData(data, remainingPath,
            showPath = showPath,
            rmNull = rmNull, rmDup = rmDup, lastObj = names(data)
          ))
        }
        return(getNestedElement(x, listName))
      }
    } else {
      return(x)
    }
  } else if (is.list(x)) {
    result <- lapply(x, findData,
      depthObj = c(depth - 1, listName), showPath = showPath,
      rmNull = rmNull, rmDup = rmDup, lastObj = names(x)
    )
    sharedListNames <- list()
    if (isTRUE(rmNull)) {
      result <- result[!sapply(result, is.null)]
      result <- result[!sapply(result, function(x) length(x) == 0)]
    }
    if (isTRUE(rmDup)) {
      uniqueResult <- list()
      uniqueCounter <- 1
      for (i in seq_along(result)) {
        duplicateFound <- FALSE
        for (j in seq_along(uniqueResult)) {
          if (identical(result[[i]], uniqueResult[[j]])) {
            duplicateFound <- TRUE
            sharedListNames[[paste0("unique_list_", j)]] <- c(sharedListNames[[paste0("unique_list_", j)]], names(result)[i])
            break
          }
        }
        if (!duplicateFound) {
          uniqueName <- paste0("unique_list_", uniqueCounter)
          uniqueResult[[names(result)[i]]] <- result[[i]]
          sharedListNames[[uniqueName]] <- names(result)[i]
          uniqueCounter <- uniqueCounter + 1
        }
      }
      result <- uniqueResult
    }

    if (isTRUE(showPath)) {
      if (length(sharedListNames) > 0 & depth == 2) result$shared_list_names <- sharedListNames
      return(result) # Carry original list structure
    } else {
      flatResult <- do.call(docall, unname(result))
      if (length(sharedListNames) > 0 & depth == 2) {
        names(result) <- paste0("unique_list_", seq_along(result))
        result$shared_list_names <- sharedListNames
        return(result)
      } else {
        return(flatResult) # Only return values
      }
    }
  } else {
    message(paste0("list ", depthObj[length(depthObj)], " is not found in ", lastObj, ".  \n"))
  }
}


thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    ## Rscript
    path <- cmdArgs[match]
    path <- gsub("\\~\\+\\~", " ", path)
    return(normalizePath(sub(needle, "", path)))
  } else {
    ## 'source'd via R console
    return(sys.frames()[[1]]$ofile)
  }
}

loadScript <- function() {
  fileName <- thisFile()
  return(ifelse(!is.null(fileName) && file.exists(fileName),
    readChar(fileName, file.info(fileName)$size), ""
  ))
}

#' Find Valid File Path
findValidFilePath <- function(referenceFilePath, targetFilePath) {
  # Check if the reference file path exits
  tryReference <- function() {
    if (file.exists(referenceFilePath)) {
      return(referenceFilePath)
    } else {
      return(NULL)
    }
  }
  # Check if the target file path exists
  tryTarget <- function() {
    if (file.exists(targetFilePath)) {
      return(targetFilePath)
    } else {
      # If not, construct a new target path by combining the directory of the reference file path with the target file path
      targetFullPath <- file.path(dirname(referenceFilePath), targetFilePath)
      if (file.exists(targetFullPath)) {
        return(targetFullPath)
      } else {
        return(NULL)
      }
    }
  }

  targetResult <- tryTarget()
  if (!is.null(targetResult)) {
    return(targetResult)
  }

  referenceResult <- tryReference()
  if (!is.null(referenceResult)) {
    return(referenceResult)
  }

  stop(sprintf(
    "Both reference and target file paths do not work. Tried paths: '%s' and '%s'",
    referenceFilePath, file.path(dirname(referenceFilePath), targetFilePath)
  ))
}

findValidFilePaths <- function(referenceFilePath, targetFilePaths) sapply(targetFilePaths, function(x) findValidFilePath(referenceFilePath, x))

#' Filter a vector based on a correlation matrix
#'
#' This function filters a vector `z` based on a correlation matrix `ld` and a correlation threshold `rThreshold`.
#' It keeps only one element among those having an absolute correlation value greater than the threshold.
#'
#' @param z A numeric vector to be filtered.
#' @param ld A square correlation matrix with dimensions equal to the length of `z`.
#' @param rThreshold The correlation threshold for filtering.
#'
#' @return A list containing the following elements:
#'   \describe{
#'     \item{filteredZ}{The filtered vector `z` based on the correlation threshold.}
#'     \item{filteredLD}{The filtered matrix `ld` based on the correlation threshold.}
#'     \item{dupBearer}{A vector indicating the duplicate status of each element in `z`.}
#'     \item{corABS}{A vector storing the absolute correlation values of duplicates.}
#'     \item{sign}{A vector storing the sign of the correlation values (-1 for negative, 1 for positive).}
#'     \item{minValue}{The minimum absolute correlation value encountered.}
#'   }
#'
#' @examples
#' z <- c(1, 2, 3, 4, 5)
#' ld <- matrix(c(
#'   1.0, 0.8, 0.2, 0.1, 0.3,
#'   0.8, 1.0, 0.4, 0.2, 0.5,
#'   0.2, 0.4, 1.0, 0.6, 0.1,
#'   0.1, 0.2, 0.6, 1.0, 0.3,
#'   0.3, 0.5, 0.1, 0.3, 1.0
#' ), nrow = 5, ncol = 5)
#' rThreshold <- 0.5
#'
#' result <- findDuplicateVariants(z, ld, rThreshold)
#' print(result)
#'
#' @export
findDuplicateVariants <- function(z, ld, rThreshold) {
  p <- length(z)
  dupBearer <- rep(-1, p)
  corABS <- rep(0, p)
  sign <- rep(1, p)
  count <- 1
  minValue <- 1

  for (i in 1:(p - 1)) {
    if (dupBearer[i] != -1) next

    idx <- (i + 1):p
    corVec <- abs(ld[i, idx])
    dupIdx <- which(dupBearer[idx] == -1 & corVec > rThreshold)

    if (length(dupIdx) > 0) {
      j <- idx[dupIdx]
      sign[j] <- ifelse(ld[i, j] < 0, -1, sign[j])
      corABS[j] <- corVec[dupIdx]
      dupBearer[j] <- count
    }

    minValue <- min(minValue, min(corVec))
    count <- count + 1
  }

  # Filter z based on dupBearer
  filteredZ <- z[dupBearer == -1]
  filteredLD <- ld[dupBearer == -1, dupBearer == -1, drop = FALSE]

  return(list(filteredZ = filteredZ, filteredLD = filteredLD, dupBearer = dupBearer, corABS = corABS, sign = sign, minValue = minValue))
}

#' Convert Z-scores to Beta and Standard Error
#'
#' This function estimates the effect sizes (beta) and standard errors (SE) from
#' given z-scores, minor allele frequencies (MAF), and a sample size (n) in genetic studies.
#' It supports vector inputs for z-scores and MAFs to process multiple variants simultaneously.
#'
#' @param z Numeric vector. The z-scores of the genetic variants.
#' @param maf Numeric vector. The minor allele frequencies of the genetic variants (0 < maf <= 0.5).
#' @param n Integer. The sample size of the study (assumed to be the same for all variants).
#'
#' @return A data frame containing three columns:
#' \describe{
#'   \item{beta}{The estimated effect sizes.}
#'   \item{se}{The estimated standard errors.}
#'   \item{maf}{The input minor allele frequencies (possibly adjusted if > 0.5).}
#' }
#'
#' @details
#' The function uses the following formulas to estimate beta and SE:
#' Beta = z / sqrt(2p(1-p)(n + z^2))
#' SE = 1 / sqrt(2p(1-p)(n + z^2))
#' Where p is the minor allele frequency.
#'
#' @examples
#' z <- c(2.5, -1.8, 3.2, 0.7)
#' maf <- c(0.3, 0.1, 0.4, 0.05)
#' n <- 10000
#' result <- zToBetaSe(z, maf, n)
#' print(result)
#' test_data_with_results <- cbind(test_data, results)
#' print(test_data_with_results)
#'
#' @note
#' This function assumes that the input z-scores are normally distributed and
#' that the genetic model is additive. It may not be accurate for rare variants
#' or in cases of imperfect imputation. The function automatically adjusts MAF > 0.5
#' to ensure it's always working with the minor allele.
#' @noRd
zToBetaSe <- function(z, maf, n) {
  if (length(z) != length(maf)) {
    stop("z and maf must be vectors of the same length")
  }
  # Ensure MAF is the minor allele frequency
  p <- pmin(maf, 1 - maf)
  denominator <- sqrt(2 * p * (1 - p) * (n + z^2))
  beta <- z / denominator
  se <- 1 / denominator
  return(data.frame(beta = beta, se = se, maf = p))
}

#' Convert Z-scores to P-values
#'
#' This function calculates p-values from given z-scores using a two-tailed normal distribution.
#' It supports vector input to process multiple z-scores simultaneously.
#'
#' @param z Numeric vector. The z-scores to be converted to p-values.
#'
#' @return A numeric vector of p-values corresponding to the input z-scores.
#'
#' @details
#' The function uses the following formula to calculate p-values:
#' p-value = 2 * Phi(-|z|)
#' Where Phi is the cumulative distribution function of the standard normal distribution.
#'
#' @examples
#' z <- c(2.5, -1.8, 3.2, 0.7)
#' pvalues <- zToPvalue(z)
#' print(pvalues)
#'
#' @note
#' This function assumes that the input z-scores are from a two-tailed test and
#' are normally distributed. It calculates two-sided p-values.
#' For extremely large absolute z-scores, the resulting p-values may be computed as zero
#' due to floating-point limitations in R. This occurs when the absolute z-score > 37.
#'
#' @export
zToPvalue <- function(z) {
  2 * pnorm(-abs(z))
}
                                                                                 
#' Filter events based on provided context name pattern
#'
#' @param events A character vector of event names
#' @param filters A data frame with character column of type_pattern, valid_pattern, and exclude_pattern.
#' @param condition Optional label context name
#' @param removeAllGroup Logical if \code{TRUE}, removes all events from the same group and character-defined context.
filterMolecularEvents <- function(events, filters, condition = NULL, removeAllGroup = FALSE) {
  # filters is a list of filter specifications
  # Each filter spec must have:
  #   type_pattern: pattern to identify event type
  #   And at least ONE of:
  #   valid_pattern: pattern that must exist in group
  #   exclude_pattern: pattern to exclude

  filteredEvents <- events
  for (filter in filters) {
    if (is.null(filter$type_pattern) ||
      (is.null(filter$valid_pattern) && is.null(filter$exclude_pattern))) {
      stop("Each filter must specify type_pattern and at least one of valid_pattern or exclude_pattern")
    }
    # Get events of this type
    typeEvents <- filteredEvents[grepl(filter$type_pattern, filteredEvents)]
    typeEventsAll <- typeEvents
    if (length(typeEvents) == 0) next
    # Apply valid pattern if specified
    if (!is.null(filter$valid_pattern)) {
      filter$valid_pattern <- strsplit(filter$valid_pattern, ",")[[1]]
      validGroups <- unique(gsub(
        filter$type_pattern, "\\1",
        typeEvents[grepl(paste(filter$valid_pattern, collapse = "|"), typeEvents)]
      ))
      if (length(validGroups) > 0) {
        typeEvents <- typeEvents[grepl(paste(filter$valid_pattern, collapse = "|"), typeEvents)] # filter for valid pattern in type events
      } else {
        typeEvents <- character(0)
      }
    }
    # Apply exclusions if specified
    if (!is.null(filter$exclude_pattern)) {
      filter$exclude_pattern <- strsplit(filter$exclude_pattern, ",")[[1]]
      typeEvents <- typeEvents[!grepl(paste(filter$exclude_pattern, collapse = "|"), typeEvents)]
    }
    if (is.null(condition)) condition <- events
    if (length(typeEvents) == length(events)) {
      message(paste("All events matching", filter$type_pattern, "in", condition, "included in following analysis."))
    } else if (length(typeEvents) == 0) {
      message(paste("No events matching", filter$type_pattern, "in", condition, "pass the filtering."))
      return(NULL)
    } else {
      excludeEvents <- paste0(setdiff(typeEventsAll, typeEvents), collapse = ";")
      message(paste("Some events,", excludeEvents, "in", condition, "are removed. \n"))
      if (removeAllGroup) {
        excludeEvents <- setdiff(typeEventsAll, typeEvents)
        excludeGroups <- gsub(filter$type_pattern, "\\1",
                              excludeEvents[grepl(paste(filter$exclude_pattern, collapse = "|"), excludeEvents)]
        )
        for (i in seq_along(excludeEvents)) {
            #if (!any(grepl(excludeGroups[i], typeEvents))) next  # skip the event if the corresponding group is all removed
            for (x in filter$exclude_pattern) excludeEvents[i] <- gsub(x, ".*", excludeEvents[i]) # remove exclude pattern from the context
            contextKey <- gsub("\\b\\d+\\b", "", excludeEvents[i]) # remove stand alone numbers (strings such as "lf2" or "chr8" will be kept)
            # General pattern to match all events of same group ID and similar character structure
            patternToRemove <- paste0(".*", excludeGroups[i], ".*")
            # Identify all events that match both the context structure and group ID
            sameGroupEvents <- typeEvents[grepl(patternToRemove, typeEvents) & grepl(gsub("\\d+", "", contextKey), gsub("\\d+", "", typeEvents))]
            typeEvents <- setdiff(typeEvents, sameGroupEvents)
        }
      }
    }
    # Update events list
    filteredEvents <- unique(c(
      filteredEvents[!grepl(filter$type_pattern, filteredEvents)],
      typeEvents
    ))
  }

  return(filteredEvents)
}


#' Robust Mahalanobis Distance
#'
#' Drop-in replacement for \code{\link[stats]{mahalanobis}} that handles
#' singular (rank-deficient) covariance matrices by falling back to the
#' Moore–Penrose pseudoinverse via \code{MASS::ginv}.
#'
#' @param x Numeric matrix (samples x features) or vector.
#' @param center Numeric vector of column means (length = number of features).
#'   If \code{NULL}, computed from \code{x}.
#' @param cov Covariance matrix. If \code{NULL}, computed from \code{x}.
#' @param inverted Logical; if \code{TRUE}, \code{cov} is already inverted.
#' @return Named numeric vector of Mahalanobis distances.
#' @importFrom MASS ginv
#' @importFrom stats cov quantile
#' @export
robustMahalanobis <- function(x, center = NULL, cov = NULL,
                              inverted = FALSE) {
  x <- if (is.vector(x)) matrix(x, ncol = length(x)) else as.matrix(x)
  if (is.null(center)) center <- colMeans(x)
  if (is.null(cov)) cov <- cov(x)
  x <- sweep(x, 2L, center)
  if (!inverted) {
    cov <- tryCatch(solve(cov), error = function(cond) {
      ginv(cov)
    })
  }
  setNames(rowSums(x %*% cov * x), rownames(x))
}

#' Detect Outliers via Mahalanobis Distance
#'
#' Identifies outlier samples in a numeric matrix (e.g., PCA scores) using
#' Mahalanobis distance with chi-squared-based p-values. Useful for QC
#' in genotype PCA or expression PCA workflows.
#'
#' @param x Numeric matrix (samples x features). Rownames are used as
#'   sample IDs in the output.
#' @param prob Numeric in (0, 1); quantile threshold for the Mahalanobis
#'   distance cutoff (default 0.99).
#' @param pvalThreshold P-value threshold for outlier classification
#'   (default 0.05). A sample is flagged only if its distance exceeds
#'   the quantile cutoff \emph{and} its p-value is below this threshold.
#' @return A data.frame with columns:
#'   \describe{
#'     \item{sample_id}{Row names from \code{x}, or row indices if unnamed.}
#'     \item{mahal}{Mahalanobis distance.}
#'     \item{pvalue}{Chi-squared p-value (df = number of features).}
#'     \item{is_outlier}{Logical; TRUE if distance > quantile cutoff and
#'       p-value < \code{pvalThreshold}.}
#'   }
#' @export
detectOutliersMahalanobis <- function(x, prob = 0.99,
                                      pvalThreshold = 0.05) {
  x <- as.matrix(x)
  sampleIds <- rownames(x) %||% as.character(seq_len(nrow(x)))
  center <- colMeans(x)
  covMat <- cov(x)
  d <- robustMahalanobis(x, center, covMat)
  p <- ncol(x)
  pvals <- pchisq(d, df = p, lower.tail = FALSE)
  cutoff <- quantile(d, probs = prob)
  data.frame(
    sample_id = sampleIds,
    mahal = as.numeric(d),
    pvalue = pvals,
    is_outlier = (d > cutoff) & (pvals < pvalThreshold),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
