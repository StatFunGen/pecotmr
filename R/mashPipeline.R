#' @title Run mashr Across Multi-Context QTL or GWAS Summary Statistics
#' @description End-to-end driver: from `{strong, random, null}` sumstats
#'   collections, builds the variant × context Bhat / Shat matrices,
#'   estimates the residual correlation (\code{Vhat}), and fits the mash
#'   model with canonical + PCA + flash + ED covariance components,
#'   returning the fitted covariance list and the estimated mixture
#'   weights.
#' @param sumStatsList Named list (or \code{SimpleList}) of
#'   \code{\link{QtlSumStats}} or \code{\link{GwasSumStats}} objects.
#'   Required names: \code{"strong"} (discovery variants), \code{"random"}
#'   (random background). Optional: \code{"null"} (for residual
#'   correlation estimation).
#' @param alpha Numeric (length 1). Variance-stabilising-transform
#'   exponent forwarded to \code{mashr::mash_set_data()}. Use
#'   \code{alpha = 0} on the BETA scale, \code{alpha = 1} on the Z scale.
#' @param residualCorrelation Optional pre-computed residual correlation
#'   matrix (\code{Vhat}). Only consulted when \code{sumStatsList$null}
#'   is absent; otherwise \code{Vhat} is estimated from the null slice
#'   via \code{mashr::estimate_null_correlation_simple()}.
#' @param nPcs Optional integer; number of principal components seeded
#'   into \code{mashr::cov_pca()}. Defaults to \code{ncol(Bhat) - 1}.
#' @param inputScale One of \code{"auto"} (default), \code{"beta"},
#'   \code{"z"}. Controls which (Bhat, Shat) pair is extracted from each
#'   sumstats entry:
#'   \describe{
#'     \item{\code{"beta"}}{Bhat = BETA, Shat = SE — the standard
#'       effect-size scale mashr was designed around. Requires every
#'       entry to carry BETA + SE mcols.}
#'     \item{\code{"z"}}{Bhat = Z, Shat = 1 — z-score scale. Requires Z.}
#'     \item{\code{"auto"}}{Use BETA + SE when every entry carries both;
#'       otherwise fall back to (Z, 1) when every entry carries Z. Mixed
#'       inputs (some entries missing BETA, others missing Z) are a
#'       hard error.}
#'   }
#'   \code{alpha} should be chosen consistently with the resolved scale:
#'   typically \code{alpha = 0} for beta, \code{alpha = 1} for z.
#' @param setSeed Integer. RNG seed for reproducibility of
#'   \code{mashr::cov_flash} and \code{mashr::cov_ed}. Default 999.
#' @return A list with elements \code{U} (the combined covariance list:
#'   canonical + PCA + flash + ED) and \code{w} (the estimated mixture
#'   weights).
#' @export
mashPipeline <- function(sumStatsList, alpha,
                          residualCorrelation = NULL,
                          nPcs = NULL,
                          inputScale = c("auto", "beta", "z"),
                          setSeed = 999) {
  inputScale <- match.arg(inputScale)
  if (!requireNamespace("mashr", quietly = TRUE)) {
    stop("To use this function, please install mashr: ",
         "https://cran.r-project.org/web/packages/mashr/index.html")
  }
  if (!requireNamespace("flashier", quietly = TRUE)) {
    stop("To use this function, please install flashier: ",
         "https://github.com/willwerscheid/flashier")
  }

  # Accept either a base list or a S4Vectors::SimpleList.
  if (methods::is(sumStatsList, "SimpleList")) {
    sumStatsList <- as.list(sumStatsList)
  }
  if (!is.list(sumStatsList) || is.null(names(sumStatsList))) {
    stop("mashPipeline: `sumStatsList` must be a named list (or SimpleList) ",
         "of QtlSumStats / GwasSumStats objects, named with at least ",
         "'strong' and 'random' (optionally 'null').")
  }
  required <- c("strong", "random")
  missingNames <- setdiff(required, names(sumStatsList))
  if (length(missingNames) > 0L) {
    stop("mashPipeline: `sumStatsList` is missing required entr",
         if (length(missingNames) == 1L) "y: " else "ies: ",
         paste(shQuote(missingNames), collapse = ", "), ".")
  }
  extraNames <- setdiff(names(sumStatsList), c("strong", "random", "null"))
  if (length(extraNames) > 0L) {
    stop("mashPipeline: `sumStatsList` has unrecognised entries: ",
         paste(shQuote(extraNames), collapse = ", "),
         ". Only 'strong', 'random', and 'null' are accepted.")
  }

  set.seed(setSeed)

  strongMats <- .mashSumStatsToMatrices(sumStatsList$strong, "strong",
                                         inputScale = inputScale)
  randomMats <- .mashSumStatsToMatrices(sumStatsList$random, "random",
                                         inputScale = inputScale)

  hasNull <- "null" %in% names(sumStatsList) && !is.null(sumStatsList$null)
  if (hasNull) {
    nullMats <- .mashSumStatsToMatrices(sumStatsList$null, "null",
                                         inputScale = inputScale)
  }

  if (!hasNull) {
    if (!is.null(residualCorrelation)) {
      vhat <- residualCorrelation
    } else {
      conditionNum <- ncol(randomMats$b)
      vhat <- diag(rep(1, conditionNum))
    }
  } else {
    vhat <- mashr::estimate_null_correlation_simple(
      mashr::mash_set_data(nullMats$b,
                           Shat = nullMats$s,
                           alpha, zero_Bhat_Shat_reset = 1000))
  }

  mashData <- mashr::mash_set_data(strongMats$b,
                                   Shat = strongMats$s,
                                   V = vhat,
                                   alpha, zero_Bhat_Shat_reset = 1000)

  # Canonical covariance matrices
  U.can <- mashr::cov_canonical(mashData)
  # PCA-based covariance matrices
  if (is.null(nPcs)) {
    nPcs <- ncol(mashData$Bhat) - 1
  }
  U.pca <- mashr::cov_pca(mashData, npc = nPcs)
  # Flash-based covariance matrices (factor analysis)
  U.flash <- mashr::cov_flash(mashData)
  # ED-based covariance matrices (initialized from all others)
  U.ed <- mashr::cov_ed(mashData, Ulist_init = c(U.can, U.pca, U.flash))
  # Combine all covariance matrices
  U.all <- c(U.can, U.pca, U.flash, U.ed)

  # Fit mash to estimate mixture weights
  m <- mashr::mash(mashData, Ulist = U.all, outputlevel = 1)
  w <- mashr::get_estimated_pi(m)

  list(U = U.all, w = w)
}

                 

            

# =============================================================================
# Mash pairwise contrast functions
# =============================================================================

#' Create a pairwise contrast column
#'
#' Sets +1 for the first condition and -1 for the second in a zero vector.
#' Used as a building block for contrast design matrices.
#'
#' @param pair A length-2 character vector naming the two conditions to contrast.
#' @param template A named numeric vector of zeros with names matching all
#'   conditions.
#' @return The template vector with +1 at \code{pair[1]} and -1 at
#'   \code{pair[2]}.
#' @export
makePairwiseContrastCol <- function(pair, template) {
  template[pair[1]] <- 1
  template[pair[2]] <- -1
  template
}

#' Compute pairwise contrasts from mash posterior
#'
#' For a single variant (row index), computes deviation contrasts (each
#' condition vs grand mean) and all pairwise contrasts from the mash posterior
#' mean and covariance. Supports condition grouping for weighted contrasts.
#'
#' @param index Integer row index of the variant in the posterior matrices.
#' @param origMean Matrix of original effect sizes (variants x conditions).
#'   Used to determine which conditions are "tested" (non-zero).
#' @param posteriorMean Matrix of mash posterior means (variants x conditions).
#' @param posteriorVcov 3D array of posterior covariance matrices
#'   (conditions x conditions x variants).
#' @param grouping Named integer vector mapping condition names to group IDs.
#'   Conditions with the same positive group ID are treated as replicates
#'   (e.g., multiple datasets for the same cell type). Use 0 for ungrouped.
#'   If NULL (default), all conditions are treated independently.
#' @return A single-row data.frame with columns
#'   \code{mean_contrast_*}, \code{se_contrast_*}, \code{p_contrast_*} for
#'   both deviation and pairwise contrasts. Returns NULL if fewer than 2
#'   tested conditions.
#' @importFrom stringr str_remove_all
#' @importFrom utils combn
#' @export
fitMashContrast <- function(index, origMean, posteriorMean, posteriorVcov,
                               grouping = NULL) {
  populationNames <- colnames(posteriorMean)
  if (!is.null(populationNames))
    populationNames <- str_remove_all(populationNames, "BETA_")

  origMeanVector <- origMean[index, ]
  names(origMeanVector) <- populationNames
  tested <- names(origMeanVector[origMeanVector != 0])

  if (length(tested) < 2) return(NULL)

  nPop <- length(tested)
  pairwiseVector <- setNames(rep(0, nPop), tested)

  # Default grouping: all independent

  if (is.null(grouping)) {
    grouping <- setNames(rep(0L, nPop), tested)
  } else {
    grouping <- grouping[tested]
  }

  if (nPop > 2) {
    # 1. Deviation contrasts
    dev <- matrix(-1, nPop, nPop, dimnames = list(tested, tested))
    diag(dev) <- nPop - 1

    # Adjust for grouped conditions
    uniqueGroups <- unique(grouping)
    for (grp in uniqueGroups[uniqueGroups > 0]) {
      grpMask <- grouping == grp
      grpSize <- sum(grpMask)
      diag(dev)[grpMask] <- (nPop - 1) / grpSize
      dev[grpMask, grpMask] <- (nPop - 1) / grpSize
    }
    colnames(dev) <- paste0(tested, "_deviation")

    # 2. Pairwise contrasts
    twoCombn <- combn(tested, 2)
    pwNames <- apply(twoCombn, 2, paste, collapse = "_vs_")
    pw <- apply(twoCombn, 2, makePairwiseContrastCol, pairwiseVector)
    colnames(pw) <- pwNames

    # Adjust pairwise contrasts for grouped conditions
    pwAdj <- pw
    for (col in colnames(pw)) {
      groups <- strsplit(col, "_vs_")[[1]]
      groupValues <- grouping[names(grouping) %in% groups]
      relevant <- names(groupValues[groupValues > 0])
      if (length(unique(groupValues)) > 1 && length(relevant) > 0) {
        for (dg in unique(groupValues[groupValues > 0])) {
          rowsInGroup <- names(grouping[grouping == dg])
          matchedRow <- rowsInGroup[rowsInGroup %in% groups]
          if (length(matchedRow) > 0)
            pwAdj[rowsInGroup, col] <- pw[matchedRow, col] / length(rowsInGroup)
        }
      }
    }

    contrastDesign <- cbind(dev / (nPop - 1), pwAdj)
  } else {
    pairwiseVector[tested[1]] <- 1
    pairwiseVector[tested[2]] <- -1
    contrastDesign <- matrix(pairwiseVector, ncol = 1,
                              dimnames = list(tested, paste0(tested[1], "_vs_", tested[2])))
  }

  # Subset posterior to tested conditions
  pm <- posteriorMean[index, tested]
  pv <- posteriorVcov[tested, tested, index]

  # Compute contrasts
  contrastDiff <- drop(t(contrastDesign) %*% pm)
  contrastVcov <- t(contrastDesign) %*% pv %*% contrastDesign
  contrastSe <- sqrt(diag(contrastVcov))
  contrastP <- 2 * (1 - pnorm(abs(contrastDiff) / contrastSe))

  # Build output data.frame
  cnames <- colnames(contrastDesign)
  df <- data.frame(
    row.names = rownames(posteriorMean)[index],
    stringsAsFactors = FALSE)
  for (i in seq_along(cnames)) {
    df[[paste0("mean_contrast_", cnames[i])]] <- contrastDiff[i]
    df[[paste0("se_contrast_", cnames[i])]] <- contrastSe[i]
    df[[paste0("p_contrast_", cnames[i])]] <- contrastP[i]
  }
  df
}

# =============================================================================
# Mash model subsetting functions
# =============================================================================

#' Subset a fitted mash model to a subset of conditions
#'
#' Updates the prior covariance matrices (\code{Ulist}) and mixture weights
#' (\code{pi}) in a fitted \code{mashr} model to match a reduced set of
#' conditions. Handles condition-specific, identity, and data-driven
#' covariance components.
#'
#' @param mashModel A fitted mash model object (from \code{mashr::mash}).
#' @param allSamples Character vector of all original condition names.
#' @param samples Character vector of the conditions to retain.
#' @return The updated mash model with resized covariance matrices and
#'   pruned mixture weights.
#' @export
updateMashModelCov <- function(mashModel, allSamples, samples) {
  cov <- mashModel$fitted_g$Ulist

  # Remove matrices for dropped conditions
  unwanted <- setdiff(allSamples, samples)
  for (d in names(cov)) {
    if (d %in% unwanted || d %in% paste0("ED_", unwanted))
      cov[[d]] <- NULL
  }

  # Resize remaining matrices to match retained conditions
  for (d in names(cov)) {
    if (d %in% samples) {
      # Condition-specific: single 1 on diagonal
      m <- matrix(0, length(samples), length(samples))
      m[which(samples == d), which(samples == d)] <- 1
      cov[[d]] <- m
    } else if (d == "identity") {
      m <- matrix(0, length(samples), length(samples))
      m[1, 1] <- 1
      cov[[d]] <- m
    } else if (is.null(colnames(cov[[d]]))) {
      cov[[d]] <- cov[[d]][seq_len(length(samples)), seq_len(length(samples))]
    } else {
      cov[[d]] <- cov[[d]][samples, samples]
    }
    cov[[d]] <- as.matrix(cov[[d]])
  }

  mashModel$fitted_g$Ulist <- cov

  # Prune mixture weights for removed conditions
  for (s in unwanted) {
    dropIdx <- grep(s, names(mashModel$fitted_g$pi), fixed = TRUE)
    if (length(dropIdx) > 0)
      mashModel$fitted_g$pi <- mashModel$fitted_g$pi[-dropIdx]
  }

  mashModel
}

#' Subset mash data matrices to specific SNPs and conditions
#'
#' Slices the \code{bhat}, \code{sbhat}, and \code{Z} matrices by row (SNPs)
#' and column (samples/conditions), and correspondingly subsets the \code{vhat}
#' covariance matrix.
#'
#' @param data A mash data list with elements \code{bhat}, \code{sbhat},
#'   \code{Z} (matrices), and \code{snp} (character vector).
#' @param vhat A square covariance matrix (conditions x conditions).
#' @param snps Character vector of SNP IDs to retain (row names).
#' @param samples Character vector of condition names to retain (column names).
#' @return A list with \code{data} (sliced data list) and \code{vhat}
#'   (sliced covariance matrix).
#' @export
sliceMashData <- function(data, vhat, snps, samples) {
  data$bhat <- as.matrix(data$bhat[snps, samples])
  data$sbhat <- as.matrix(data$sbhat[snps, samples])
  data$Z <- as.matrix(data$Z[snps, samples])
  vhat <- as.matrix(vhat[samples, samples])
  data$snp <- data$snp[data$snp %in% snps]
  colnames(data$bhat) <- colnames(data$sbhat) <- colnames(data$Z) <- colnames(vhat) <- samples
  list(data = data, vhat = vhat)
}

#' Sanitize NaN/Inf values in mash data
#'
#' Replaces NaN in \code{bhat} with 0 and NaN/Inf in \code{sbhat} with 1e3
#' (indicating high uncertainty).
#'
#' @param data A mash data list with \code{bhat} and \code{sbhat} matrices.
#' @return The data list with sanitized values.
#' @export
sanitizeMashData <- function(data) {
  data$bhat[is.nan(data$bhat)] <- 0
  data$sbhat[is.nan(data$sbhat) | is.infinite(data$sbhat)] <- 1e3
  data
}

#' Random-Effects Meta-Analysis of Mash Pairwise Contrasts
#'
#' For each cell type (condition), gathers all pairwise contrast effect
#' sizes and standard errors involving that cell, then runs a
#' DerSimonian–Laird random-effects meta-analysis per condition.
#' Intended to be run on the output of \code{\link{fitMashContrast}}.
#'
#' @param effectSizes Numeric matrix (features x conditions) of contrast
#'   effect sizes. Column names must follow the pattern
#'   \code{mean_contrast_<cellA>_vs_<cellB>}.
#' @param seValues Numeric matrix (features x conditions) of contrast
#'   standard errors. Must have the same dimensions and column names as
#'   \code{effectSizes}.
#' @param seCutoff Numeric; minimum SE below which a condition is excluded
#'   from the meta-analysis for a given feature (default 0).
#' @return A tibble with columns:
#'   \describe{
#'     \item{cell}{Cell type name.}
#'     \item{condition}{Original pairwise contrast name (without prefix).}
#'     \item{meta_pvalue}{P-value from the random-effects meta-analysis.}
#'     \item{meta_effect}{Pooled absolute effect size estimate.}
#'     \item{meta_se}{Standard error of the pooled estimate.}
#'     \item{tau2}{Between-study variance estimate.}
#'     \item{I2}{Heterogeneity measure (proportion of variance due to
#'       between-study variance), in [0, 1].}
#'   }
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
#' @export
metaAnalysisPerCell <- function(effectSizes, seValues,
                                   seCutoff = 0) {
  stopifnot(identical(dim(effectSizes), dim(seValues)))
  stopifnot(identical(colnames(effectSizes), colnames(seValues)))

  conditions <- sub("^mean_contrast_", "", colnames(effectSizes))
  cells <- unique(c(sub("_vs_.*", "", conditions),
                     sub(".*_vs_", "", conditions)))

  results <- list()
  for (cell in cells) {
    # Columns involving this cell
    cellIdx <- grep(cell, colnames(effectSizes))
    if (length(cellIdx) == 0) next

    cellEffects <- effectSizes[, cellIdx, drop = FALSE]
    cellSes <- seValues[, cellIdx, drop = FALSE]
    cellConditions <- conditions[cellIdx]

    for (i in seq_along(cellConditions)) {
      es <- abs(as.numeric(cellEffects[, i]))
      se <- as.numeric(cellSes[, i])

      # Filter by SE cutoff
      keep <- se > seCutoff & is.finite(es) & is.finite(se)
      es <- es[keep]
      se <- se[keep]

      if (length(es) < 2) {
        results[[length(results) + 1]] <- tibble(
          cell = cell,
          condition = cellConditions[i],
          meta_pvalue = if (length(es) == 1) {
            2 * pnorm(abs(es / se), lower.tail = FALSE)
          } else NA_real_,
          meta_effect = if (length(es) == 1) es else NA_real_,
          meta_se = if (length(es) == 1) se else NA_real_,
          tau2 = NA_real_,
          I2 = NA_real_
        )
        next
      }

      ma <- metaRandomEffects(es, se)
      z <- ma$mean / ma$se
      results[[length(results) + 1]] <- tibble(
        cell = cell,
        condition = cellConditions[i],
        meta_pvalue = 2 * pnorm(abs(z), lower.tail = FALSE),
        meta_effect = ma$mean,
        meta_se = ma$se,
        tau2 = ma$tau2,
        I2 = ma$I2
      )
    }
  }

  bind_rows(results)
}
