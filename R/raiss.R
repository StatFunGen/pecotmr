#' Core RAISS implementation for a single LD matrix
#'
#' @param refPanel A data frame containing 'chrom', 'pos', 'variant_id', 'A1', and 'A2'.
#' @param knownZscores A data frame containing 'chrom', 'pos', 'variant_id', 'A1', 'A2', and 'z' values.
#' @param ldMatrix A square matrix of dimension equal to the number of rows in refPanel.
#' @param lamb Regularization term added to the diagonal of the ldMatrix.
#' @param rcond Threshold for filtering eigenvalues in the pseudo-inverse computation.
#' @param r2Threshold R square threshold below which SNPs are filtered from the output.
#' @param minimumLd Minimum LD score threshold for SNP filtering.
#' @param verbose Logical indicating whether to print progress information.
#'
#' @return A list containing filtered and unfiltered results, and filtered LD matrix.
#' @importFrom MASS ginv
#' @importFrom dplyr arrange
#' @noRd
raissSingleMatrix <- function(refPanel, knownZscores, ldMatrix, lamb = 0.01, rcond = 0.01,
                              r2Threshold = 0.6, minimumLd = 5, verbose = TRUE) {
  # Check that refPanel and knownZscores are both increasing in terms of pos
  if (is.unsorted(refPanel$pos) || is.unsorted(knownZscores$pos)) {
    stop("refPanel and knownZscores must be in increasing order of pos.")
  }

  # Convert ldMatrix to matrix if it's a data frame
  if (is.data.frame(ldMatrix)) {
    ldMatrix <- as.matrix(ldMatrix)
  }

  # Define knowns and unknowns
  knownsId <- intersect(knownZscores$variant_id, refPanel$variant_id)
  knowns <- which(refPanel$variant_id %in% knownsId)
  unknowns <- which(!refPanel$variant_id %in% knownsId)

  # Handle edge cases
  if (length(knowns) == 0) {
    if (verbose) message("No known variants found, cannot perform imputation.")
    return(NULL)
  }

  if (length(unknowns) == 0) {
    if (verbose) message("No unknown variants to impute, returning known variants.")
    return(list(
      resultNofilter = knownZscores,
      resultFilter = knownZscores,
      ldMat = ldMatrix
    ))
  }

  # Extract zt, sigT, and sigIT
  zt <- knownZscores$z
  sigT <- ldMatrix[knowns, knowns, drop = FALSE]
  sigIT <- ldMatrix[unknowns, knowns, drop = FALSE]

  # Call raissModel
  results <- raissModel(zt, sigT, sigIT, lamb, rcond)
  # Format the results
  results <- formatRaissDf(results, refPanel, unknowns)
  # Filter output
  results <- filterRaissOutput(results, r2Threshold, minimumLd, verbose)

  # Merge with known z-scores
  resultNofilter <- mergeRaissDf(results$zscoresNofilter, knownZscores) %>% arrange(pos)
  resultFilter <- mergeRaissDf(results$zscores, knownZscores) %>% arrange(pos)

  # Filter out variants not included in the imputation result
  filteredOutVariant <- setdiff(refPanel$variant_id, resultFilter$variant_id)

  # Update the LD matrix excluding filtered variants
  ldExtractFiltered <- if (length(filteredOutVariant) > 0) {
    filteredOutId <- match(filteredOutVariant, refPanel$variant_id)
    as.matrix(ldMatrix)[-filteredOutId, -filteredOutId]
  } else {
    as.matrix(ldMatrix)
  }
  # Return results
  return(list(
    resultNofilter = resultNofilter,
    resultFilter = resultFilter,
    ldMat = ldExtractFiltered
  ))
}

#' Core RAISS implementation from a genotype matrix X (SVD-based)
#'
#' Performs the same imputation as \code{raissSingleMatrix} but works directly
#' with the genotype matrix X instead of the LD correlation matrix R. This avoids
#' forming the p x p LD matrix, saving O(p^2) memory and O(np^2) compute.
#'
#' The reformulation is mathematically exact: using the thin SVD of Xt (the known
#' variant columns), all RAISS quantities (mu, var, ld_score) are computed in the
#' SVD basis without ever forming R = X'X/(n-1).
#'
#' @param refPanel A data frame containing 'chrom', 'pos', 'variant_id', 'A1', and 'A2'.
#' @param knownZscores A data frame containing 'chrom', 'pos', 'variant_id', 'A1', 'A2', and 'z' values.
#' @param X Centered and scaled genotype matrix (nSamples x pVariants). Column order must
#'   match the variant order in refPanel.
#' @param lamb Regularization term (same role as in the LD-based path).
#' @param svdTol Relative tolerance for filtering small singular values in the SVD of Xt.
#' @param r2Threshold R square threshold below which SNPs are filtered from the output.
#' @param minimumLd Minimum LD score threshold for SNP filtering.
#' @param verbose Logical indicating whether to print progress information.
#'
#' @return A list containing filtered and unfiltered results, and ldMat = NULL.
#' @importFrom dplyr arrange
#' @noRd
raissSingleMatrixFromX <- function(refPanel, knownZscores, X, lamb = 0.01,
                                   svdTol = 1e-8, r2Threshold = 0.6,
                                   minimumLd = 5, verbose = TRUE) {
  # Check that refPanel and knownZscores are both increasing in terms of pos
  if (is.unsorted(refPanel$pos) || is.unsorted(knownZscores$pos)) {
    stop("refPanel and knownZscores must be in increasing order of pos.")
  }

  nSamples <- nrow(X)

  # Define knowns and unknowns (same logic as raissSingleMatrix)
  knownsId <- intersect(knownZscores$variant_id, refPanel$variant_id)
  knowns <- which(refPanel$variant_id %in% knownsId)
  unknowns <- which(!refPanel$variant_id %in% knownsId)

  # Handle edge cases
  if (length(knowns) == 0) {
    if (verbose) message("No known variants found, cannot perform imputation.")
    return(NULL)
  }

  if (length(unknowns) == 0) {
    if (verbose) message("No unknown variants to impute, returning known variants.")
    return(list(
      resultNofilter = knownZscores,
      resultFilter = knownZscores,
      ldMat = NULL
    ))
  }

  # Extract known columns for SVD (unavoidable copy for LAPACK).
  # We do NOT copy X_i - instead we compute X' %*% [w|U] on the full X
  # and index the unknown rows, saving O(n*m) memory.
  Xt <- X[, knowns, drop = FALSE]
  zt <- knownZscores$z

  # Compute thin SVD of Xt (n x k -> U: n x r, d: r, V: k x r)
  svdResult <- safeSvd(Xt, tol = svdTol)
  U <- svdResult$u
  d <- svdResult$d
  V <- svdResult$v
  rm(Xt)  # free n*k memory; no longer needed

  # Precompute regularization and weight vectors (length r, cheap)
  cReg <- lamb * (nSamples - 1)
  d2 <- d^2
  d2PlusC <- d2 + cReg

  # --- Build w (n x 1): the projection of zt through the regularized SVD ---
  # w = U %*% diag(d / (d^2 + c)) %*% V' zt
  VtZt <- crossprod(V, zt)                       # r x 1
  w <- U %*% (d / d2PlusC * VtZt)                # n x 1

  # --- Single BLAS call: X' %*% [w | U] -> p x (1+r) ---
  # This avoids copying X_i (n x m) entirely.
  # Row unknowns of column 1 gives mu; rows unknowns of columns 2:(r+1) gives A.
  XtWU <- crossprod(X, cbind(w, U))               # p x (1+r), one dgemm call
  mu <- as.numeric(XtWU[unknowns, 1])             # m x 1
  A <- XtWU[unknowns, -1, drop = FALSE]           # m x r (subset, not copy of X)
  rm(XtWU)                                         # free p*(1+r)

  # --- Variance and LD score in one pass over A^2 ---
  # var = (1+lamb) - (1/(n-1)) * A^2 %*% (d^2/(d^2+c))
  # ld_score = (1/(n-1))^2 * A^2 %*% d^2
  # Compute A^2 once, multiply by [d_weights_var | d^2] in one dgemm.
  ASq <- A^2                                       # m x r (one allocation)
  rm(A)                                            # free m*r
  dWeights <- cbind(d2 / d2PlusC, d2)              # r x 2
  scores <- ASq %*% dWeights                       # m x 2 (one dgemm)
  rm(ASq)                                          # free m*r

  nm1 <- nSamples - 1
  varRaw <- (1 + lamb) - scores[, 1] / nm1
  raissLdScore <- scores[, 2] / nm1^2
  rm(scores)

  # --- Condition number (scalar, expanded to vector by formatRaissDf) ---
  conditionNumber <- rep(d[1] / d[length(d)], length(unknowns))
  correctInversion <- rep(TRUE, length(unknowns))

  # --- R2 correction (same as raissModel) ---
  varNorm <- varInBoundaries(varRaw, lamb)
  R2 <- (1 + lamb) - varNorm
  mu <- mu / sqrt(R2)

  # Package results in the same format as raissModel output
  imp <- list(
    var = varNorm,
    mu = mu,
    raissLdScore = raissLdScore,
    conditionNumber = conditionNumber,
    correctInversion = correctInversion
  )

  # Reuse existing formatting and filtering functions
  results <- formatRaissDf(imp, refPanel, unknowns)
  results <- filterRaissOutput(results, r2Threshold, minimumLd, verbose)

  # Merge with known z-scores
  resultNofilter <- mergeRaissDf(results$zscoresNofilter, knownZscores) %>% arrange(pos)
  resultFilter <- mergeRaissDf(results$zscores, knownZscores) %>% arrange(pos)

  return(list(
    resultNofilter = resultNofilter,
    resultFilter = resultFilter,
    ldMat = NULL
  ))
}

#' Impute Summary Statistics Using LD (RAISS)
#'
#' This function is a part of the statistical library for SNP imputation from:
#' https://gitlab.pasteur.fr/statistical-genetics/raiss/-/blob/master/raiss/stat_models.py
#' It is R implementation of the imputation model described in the paper by Bogdan Pasaniuc,
#' Noah Zaitlen, et al., titled "Fast and accurate imputation of summary
#' statistics enhances evidence of functional enrichment", published in
#' Bioinformatics in 2014.
#'
#' This function can process either a single LD matrix or a list of LD matrices for different blocks.
#' For a list of matrices, it processes each block separately and combines the results.
#' Alternatively, it can accept a genotype matrix X directly, avoiding the need to form
#' the p x p LD matrix (memory and compute savings when n << p).
#'
#' @param refPanel A data frame containing 'chrom', 'pos', 'variant_id', 'A1', and 'A2'.
#' @param knownZscores A data frame containing 'chrom', 'pos', 'variant_id', 'A1', 'A2', and 'z' values.
#' @param ldMatrix Either a square matrix or a list of matrices for LD blocks.
#'   Provide either \code{ldMatrix} or \code{genotypeMatrix}, not both.
#' @param genotypeMatrix A centered and scaled genotype matrix (n x p) as an alternative
#'   to \code{ldMatrix}. Column order must match the variant order in \code{refPanel}.
#'   When provided, the imputation uses an SVD-based approach that avoids forming the
#'   p x p LD matrix.
#' @param lamb Regularization term added to the diagonal of the ldMatrix.
#' @param rcond Threshold for filtering eigenvalues in the pseudo-inverse computation
#'   (only used with ldMatrix path).
#' @param svdTol Relative tolerance for filtering small singular values
#'   (only used with genotypeMatrix path).
#' @param r2Threshold R square threshold below which SNPs are filtered from the output.
#' @param minimumLd Minimum LD score threshold for SNP filtering.
#' @param verbose Logical indicating whether to print progress information.
#'
#' @return A list containing filtered and unfiltered results, and filtered LD matrix
#'   (ldMat is NULL when using genotypeMatrix path).
#' @importFrom dplyr arrange bind_rows
#' @export
raiss <- function(refPanel, knownZscores, ldMatrix = NULL,
                  genotypeMatrix = NULL, lamb = 0.01, rcond = 0.01,
                  svdTol = 1e-8, r2Threshold = 0.6, minimumLd = 5,
                  verbose = TRUE) {
  # --- Genotype matrix path (SVD-based, avoids forming R) ---
  if (!is.null(genotypeMatrix)) {
    if (!is.null(ldMatrix)) {
      stop("Provide either ldMatrix or genotypeMatrix, not both.")
    }
    if (is.matrix(genotypeMatrix)) {
      if (verbose) message("Processing genotype matrix via SVD-based imputation...")
      return(raissSingleMatrixFromX(
        refPanel, knownZscores, genotypeMatrix,
        lamb, svdTol, r2Threshold, minimumLd, verbose
      ))
    }
    if (is.list(genotypeMatrix)) {
      # List of genotype matrices (block processing)
      if (verbose) message("Processing multiple genotype matrix blocks via SVD-based imputation...")
      resultsList <- list()
      for (i in seq_along(genotypeMatrix)) {
        if (verbose) message(paste("Processing block", i, "of", length(genotypeMatrix)))
        blockResult <- raissSingleMatrixFromX(
          refPanel, knownZscores, genotypeMatrix[[i]],
          lamb, svdTol, r2Threshold, minimumLd,
          verbose = FALSE
        )
        if (!is.null(blockResult)) {
          resultsList[[length(resultsList) + 1]] <- blockResult
        }
      }
      if (length(resultsList) == 0) {
        if (verbose) message("No blocks could be processed.")
        return(NULL)
      }
      combinedNofilter <- do.call(bind_rows, lapply(resultsList, `[[`, "resultNofilter"))
      combinedFilter <- do.call(bind_rows, lapply(resultsList, `[[`, "resultFilter"))
      return(list(
        resultNofilter = combinedNofilter %>% arrange(pos),
        resultFilter = combinedFilter %>% arrange(pos),
        ldMat = NULL
      ))
    }
    stop("genotypeMatrix must be a matrix or a list of matrices.")
  }

  # --- LD matrix path (original implementation) ---
  if (is.null(ldMatrix)) {
    stop("Provide either ldMatrix or genotypeMatrix.")
  }
  # Determine if we can process as a single matrix
  isSingleMatrixCase <- is.matrix(ldMatrix) ||
    (is.list(ldMatrix) && !is.null(ldMatrix$ldMatrices) &&
      length(ldMatrix$ldMatrices) == 1)

  if (isSingleMatrixCase) {
    if (verbose) message("Processing single LD matrix", if (!is.matrix(ldMatrix)) " from list", "...")

    # Extract the matrix if it's in a list
    if (!is.matrix(ldMatrix)) {
      ldMatrix <- ldMatrix$ldMatrices[[1]]
    }

    return(raissSingleMatrix(
      refPanel, knownZscores, ldMatrix,
      lamb, rcond, r2Threshold, minimumLd, verbose
    ))
  }

  # For list of matrices, process each block
  if (verbose) message("Processing multiple LD blocks...")

  combineWithBoundaryCheck <- function(combinedResult, newResult) {
    # If either is empty, simply return the non-empty one or empty data frame
    if (is.null(combinedResult)) {
      return(newResult)
    }
    if (is.null(newResult)) {
      return(combinedResult)
    }

    # Check if the last variant of combined matches the first of new
    lastVar <- combinedResult$variant_id[nrow(combinedResult)]
    firstVar <- newResult$variant_id[1]

    if (lastVar == firstVar) {
      newR2 <- newResult$raissR2[1]
      oldR2 <- combinedResult$raissR2[nrow(combinedResult)]
      if (is.na(newR2) && is.na(oldR2)) {
        # Both are NA - keep the existing one
      } else if (is.na(oldR2)) {
        # Old is NA but new is not - use new
        combinedResult[nrow(combinedResult), ] <- newResult[1, ]
      } else if (is.na(newR2)) {
        # New is NA but old is not - keep old
      } else if (newR2 > oldR2) {
        # Both are non-NA and new is better - use new
        combinedResult[nrow(combinedResult), ] <- newResult[1, ]
      }

      # Add remaining rows from new (excluding first)
      if (nrow(newResult) > 1) {
        combinedResult <- bind_rows(combinedResult, newResult[-1, ])
      }
    } else {
      # No overlap - combine all rows
      combinedResult <- bind_rows(combinedResult, newResult)
    }

    return(combinedResult)
  }

  resultsList <- list()
  variantIndices <- ldMatrix$variantIndices
  blockIds <- unique(variantIndices$blockId)

  for (blockId in blockIds) {
    if (verbose) message(paste("Processing block", blockId, "of", length(blockIds)))

    blockVariantIds <- variantIndices$variant_id[variantIndices$blockId == blockId]

    # Subset refPanel and ldMatrix for this block
    blockIndices <- match(blockVariantIds, refPanel$variant_id)
    blockRefPanel <- refPanel[blockIndices, ]
    blockLdMatrix <- ldMatrix$ldMatrices[[blockId]]
    blockKnownZscores <- knownZscores %>% filter(variant_id %in% blockVariantIds)
    if (nrow(blockLdMatrix) != nrow(blockRefPanel)) {
      stop(paste("Block", blockId, ": LD matrix dimension does not match number of variants in reference panel"))
    }

    # Process the block using the core function
    blockResult <- raissSingleMatrix(
      blockRefPanel, blockKnownZscores, blockLdMatrix,
      lamb, rcond, r2Threshold, minimumLd,
      verbose = FALSE
    )
    # Skip if block returned NULL (no known variants)
    if (!is.null(blockResult)) {
      resultsList[[blockId]] <- blockResult
    }
  }

  if (length(resultsList) == 0) {
    if (verbose) message("No blocks could be processed. Check that knownZscores overlap with variants in the blocks.")
    return(NULL)
  }

  # Combine results sequentially to handle boundary duplicates
  combinedNofilter <- resultsList[[1]]$resultNofilter
  combinedFilter <- resultsList[[1]]$resultFilter

  if (length(resultsList) > 1) {
    for (i in 2:length(resultsList)) {
      combinedNofilter <- combineWithBoundaryCheck(
        combinedNofilter,
        resultsList[[i]]$resultNofilter
      )

      combinedFilter <- combineWithBoundaryCheck(
        combinedFilter,
        resultsList[[i]]$resultFilter
      )
    }
  }

  ldFilteredList <- lapply(resultsList, function(x) x$ldMat)
  variantList <- lapply(ldFilteredList, function(ld) data.frame(variants = colnames(ld)))
  ldMatrix <- createLdMatrix(
    ldMatrices = ldFilteredList,
    variants = variantList
  )

  return(list(
    resultNofilter = combinedNofilter,
    resultFilter = combinedFilter,
    ldMat = ldMatrix
  ))
}

#' @param zt Vector of known z scores.
#' @param sigT Matrix of known linkage disequilibrium (LD) correlation.
#' @param sigIT Correlation matrix with rows corresponding to unknown SNPs (to impute)
#'               and columns to known SNPs.
#' @param lamb Regularization term added to the diagonal of the sigT matrix.
#' @param rcond Threshold for filtering eigenvalues in the pseudo-inverse computation.
#' @param batch Boolean indicating whether batch processing is used.
#'
#' @return A list containing the variance 'var', estimation 'mu', LD score 'raissLdScore',
#'         condition number 'conditionNumber', and correctness of inversion
#'         'correctInversion'.
#' @noRd
raissModel <- function(zt, sigT, sigIT, lamb = 0.01, rcond = 0.01, batch = TRUE, reportConditionNumber = FALSE) {
  sigTInv <- invertMatRecursive(sigT, lamb, rcond)
  if (!is.numeric(zt) || !is.numeric(sigT) || !is.numeric(sigIT)) {
    stop("zt, sigT, and sigIT must be numeric.")
  }
  if (batch) {
    conditionNumber <- if (reportConditionNumber) rep(kappa(sigT, exact = TRUE, norm = "2"), nrow(sigIT)) else NA
    correctInversion <- rep(checkInversion(sigT, sigTInv), nrow(sigIT))
  } else {
    conditionNumber <- if (reportConditionNumber) kappa(sigT, exact = TRUE, norm = "2") else NA
    correctInversion <- checkInversion(sigT, sigTInv)
  }

  varRaissLdScore <- computeVar(sigIT, sigTInv, lamb, batch)
  var <- varRaissLdScore$var
  raissLdScore <- varRaissLdScore$raissLdScore

  mu <- computeMu(sigIT, sigTInv, zt)
  varNorm <- varInBoundaries(var, lamb)

  R2 <- ((1 + lamb) - varNorm)
  mu <- mu / sqrt(R2)

  return(list(var = varNorm, mu = mu, raissLdScore = raissLdScore, conditionNumber = conditionNumber, correctInversion = correctInversion))
}

#' @param imp is the output of raissModel()
#' @param refPanel is a data frame with columns 'chrom', 'pos', 'variant_id', 'ref', and 'alt'.
#' @noRd
formatRaissDf <- function(imp, refPanel, unknowns) {
  resultDf <- data.frame(
    chrom = refPanel[unknowns, "chrom"],
    pos = refPanel[unknowns, "pos"],
    variant_id = refPanel[unknowns, "variant_id"],
    A1 = refPanel[unknowns, "A1"],
    A2 = refPanel[unknowns, "A2"],
    z = imp$mu,
    Var = imp$var,
    raissLdScore = imp$raissLdScore,
    conditionNumber = imp$conditionNumber,
    correctInversion = imp$correctInversion
  )

  # Specify the column order
  columnOrder <- c(
    "chrom", "pos", "variant_id", "A1", "A2", "z", "Var", "raissLdScore", "conditionNumber",
    "correctInversion"
  )

  # Reorder the columns
  resultDf <- resultDf[, columnOrder]
  return(resultDf)
}

mergeRaissDf <- function(raissDf, knownZscores) {
  # Merge the data frames
  mergedDf <- merge(raissDf, knownZscores, by = c("chrom", "pos", "variant_id", "A1", "A2"), all = TRUE)

  # Identify rows that came from knownZscores
  fromKnown <- !is.na(mergedDf$z.y) & is.na(mergedDf$z.x)

  # Set Var to -1 and raissLdScore to Inf for these rows
  mergedDf$Var[fromKnown] <- -1
  mergedDf$raissLdScore[fromKnown] <- Inf

  # If there are overlapping columns (e.g., z.x and z.y), resolve them
  # For example, use z from knownZscores where available, otherwise use z from raissDf
  mergedDf$z <- ifelse(fromKnown, mergedDf$z.y, mergedDf$z.x)

  # Remove the extra columns resulted from the merge (e.g., z.x, z.y)
  mergedDf <- mergedDf[, !colnames(mergedDf) %in% c("z.x", "z.y")]
  mergedDf <- arrange(mergedDf, pos)
  # assign imputed variants beta, se as NA to avoid confusion, since they are not imputed
  mergedDf$beta[mergedDf$Var == -1] <- NA
  mergedDf$se[mergedDf$Var == -1] <- NA
  return(mergedDf)
}

filterRaissOutput <- function(zscores, r2Threshold = 0.6, minimumLd = 5, verbose = TRUE) {
  # Reset the index and subset the data frame
  zscores <- zscores[, c("chrom", "pos", "variant_id", "A1", "A2", "z", "Var", "raissLdScore")]
  zscores$raissR2 <- 1 - zscores$Var

  # Count statistics before filtering
  nSnpsBfFilt <- nrow(zscores)
  nSnpsInitial <- sum(zscores$raissR2 == 2.0, na.rm = TRUE)
  nSnpsImputed <- sum(zscores$raissR2 != 2.0, na.rm = TRUE)
  nSnpsLdFilt <- sum(zscores$raissLdScore < minimumLd, na.rm = TRUE)
  nSnpsR2Filt <- sum(zscores$raissR2 < r2Threshold, na.rm = TRUE)

  # Apply filters
  zscoresNofilter <- zscores
  zscores <- zscores[zscores$raissR2 > r2Threshold & zscores$raissLdScore >= minimumLd, ]
  nSnpsAfFilt <- nrow(zscores)

  # Print report
  if (verbose) {
    maxLabelLength <- max(nchar(c(
      "Variants before filter:",
      "Non-imputed variants:",
      "Imputed variants:",
      "Variants filtered because of low LD score:",
      "Variants filtered because of low R2:",
      "Remaining variants after filter:"
    )))

    formatLine <- function(label, value) {
      sprintf("%-*s %d", maxLabelLength, paste0(label, ":"), value)
    }

    message("IMPUTATION REPORT\n")
    message(formatLine("Variants before filter", nSnpsBfFilt))
    message(formatLine("Non-imputed variants", nSnpsInitial))
    message(formatLine("Imputed variants", nSnpsImputed))
    message(formatLine("Variants filtered because of low LD score", nSnpsLdFilt))
    message(formatLine("Variants filtered because of low R2", nSnpsR2Filt))
    message(formatLine("Remaining variants after filter", nSnpsAfFilt))
  }
  return(zscore_list = list(zscoresNofilter = zscoresNofilter, zscores = zscores))
}

computeMu <- function(sigIT, sigTInv, zt) {
  return(sigIT %*% (sigTInv %*% zt))
}

computeVar <- function(sigIT, sigTInv, lamb, batch = TRUE) {
  if (batch) {
    var <- (1 + lamb) - rowSums((sigIT %*% sigTInv) * sigIT)
    raissLdScore <- rowSums(sigIT^2)
  } else {
    var <- (1 + lamb) - (sigIT %*% (sigTInv %*% t(sigIT)))
    raissLdScore <- sum(sigIT^2)
  }
  return(list(var = var, raissLdScore = raissLdScore))
}

checkInversion <- function(sigT, sigTInv) {
  return(all.equal(sigT, sigT %*% (sigTInv %*% sigT), tolerance = 1e-5))
}

varInBoundaries <- function(var, lamb) {
  var[var < 0] <- 0
  var[var > (0.99999 + lamb)] <- 1
  return(var)
}

invertMat <- function(mat, lamb, rcond) {
  tryCatch(
    {
      # Modify the diagonal elements of mat
      diag(mat) <- 1 + lamb
      # Compute the pseudo-inverse
      matInv <- ginv(mat, tol = rcond)
      return(matInv)
    },
    error = function(e) {
      # Second attempt with updated lamb and rcond in case of an error
      diag(mat) <- 1 + lamb * 1.1
      matInv <- ginv(mat, tol = rcond * 1.1)
      return(matInv)
    }
  )
}

invertMatRecursive <- function(mat, lamb, rcond) {
  tryCatch(
    {
      # Modify the diagonal elements of mat
      diag(mat) <- 1 + lamb
      # Compute the pseudo-inverse
      matInv <- ginv(mat, tol = rcond)
      return(matInv)
    },
    error = function(e) {
      # Recursive call with updated lamb and rcond in case of an error
      invertMat(mat, lamb * 1.1, rcond * 1.1)
    }
  )
}

invertMatEigen <- function(mat, tol = 1e-3) {
  eigenMat <- eigen(mat)
  L <- which(cumsum(eigenMat$values) / sum(eigenMat$values) > 1 - tol)[1]
  if (is.na(L)) {
    # all eigen values are extremely small
    stop("Cannot invert the input matrix because all its eigen values are negative or close to zero")
  }
  matInv <- eigenMat$vectors[, 1:L] %*%
    diag(1 / eigenMat$values[1:L]) %*%
    t(eigenMat$vectors[, 1:L])

  return(matInv)
}
