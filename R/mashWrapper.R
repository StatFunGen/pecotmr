# Filter rows of a z-score matrix by significance p-value cutoff.
# Returns integer indices of rows where any |z| exceeds the threshold.
# @noRd
filterBySignificance <- function(zMatrix, sigPCutoff) {
  zThreshold <- sqrt(qchisq(sigPCutoff, df = 1, lower.tail = FALSE))
  which(apply(zMatrix, 1, function(row) any(abs(row) >= zThreshold)))
}

#' @importFrom vroom vroom
#' @export
filterInvalidSummaryStat <- function(datList, bhat = NULL, sbhat = NULL, z = NULL, btoz = FALSE, sigPCutoff = 1E-6, filterByMissingRate = 0.2) {
  replaceValues <- function(df, replaceWith) {
    df <- df %>%
      mutate(across(everything(), as.numeric)) %>%
      mutate(across(everything(), ~ replace(., is.nan(.) | is.infinite(.) | is.na(.), replaceWith)))
  }
  # Function to process bhat, sbhat
  if (!is.null(bhat) && !is.null(sbhat) && all(c(bhat, sbhat) %in% names(datList))) {
    # If the element is a list with 'bhat' and 'sbhat'
    if (!is.null(datList[[bhat]]) && !is.null(datList[[sbhat]])) {
      datList[[bhat]] <- as.matrix(replaceValues(datList[[bhat]], 0))
      datList[[sbhat]] <- as.matrix(replaceValues(datList[[sbhat]], 1000))
      if (("null.b" %in% names(datList)) || ("random.b" %in% names(datList))) {
        if (!is.null(filterByMissingRate)) {
          proportionNonzero <- apply(datList[[bhat]], 1, function(row) {
            mean(row != 0)
          })
          datList[[bhat]] <- datList[[bhat]][proportionNonzero >= filterByMissingRate, ]
          datList[[sbhat]] <- datList[[sbhat]][proportionNonzero >= filterByMissingRate, ]
        }
      }
    }
  }
  # Function to filter strong signal using z score
  if (btoz) {
    if (any(grepl("\\.b$", bhat)) | any(grepl("\\.s$", sbhat))) {
      condition <- sub("\\.b$", "", bhat)
      if (!is.null(datList[[bhat]]) && !is.null(datList[[sbhat]])) {
        datList[[paste0(condition, ".z")]] <- as.matrix(datList[[bhat]] / datList[[sbhat]])
      } else {
        datList[paste0(condition, ".z")] <- list(NULL)
      }
    } else {
      if (!is.null(datList[[bhat]]) && !is.null(datList[[sbhat]])) {
        datList[["z"]] <- as.matrix(datList[[bhat]] / datList[[sbhat]])
      } else {
        datList["z"] <- list(NULL)
      }
    }
    if ("strong.z" %in% names(datList)) {
      if (!is.null(sigPCutoff)) {
        keepIndex <- filterBySignificance(datList$strong.z, sigPCutoff)
        datList[["strong.z"]] <- datList$strong.z[keepIndex, ]
        datList[["strong.b"]] <- datList$strong.b[keepIndex, ]
        datList[["strong.s"]] <- datList$strong.s[keepIndex, ]
      }
    }
  }
  # Function to process z-scores and filter directly
  if (!is.null(z)) {
    processZ <- function(zData) {
      zData <- as.matrix(replaceValues(zData, 0))

      if (!is.null(filterByMissingRate)) {
        proportionNonzero <- apply(zData, 1, function(row) mean(row != 0))
        zData <- zData[proportionNonzero >= filterByMissingRate, , drop = FALSE]
      }

      return(zData)
    }

    # Process each component if it exists
    for (comp in c("strong", "random", "null")) {
      if (!is.null(datList[[comp]]) && !is.null(datList[[comp]]$z)) {
        datList[[comp]]$z <- processZ(datList[[comp]]$z)
      }
    }

    # Apply significance cutoff to strong signals if applicable
    if (!is.null(datList$strong) && !is.null(datList$strong$z) && !is.null(sigPCutoff)) {
      keepIndex <- filterBySignificance(datList$strong$z, sigPCutoff)
      datList$strong$z <- datList$strong$z[keepIndex, , drop = FALSE]
    }
  }

  return(datList)
}

#' @importFrom purrr keep
#' @export
filterMixtureComponents <- function(conditionsToKeep, U, w = NULL, wCutoff = 1e-04) {
  # Identify conditions not to keep (to be removed)
  conditionsToFilter <- setdiff(colnames(U[[1]]), conditionsToKeep)
  sumW <- sum(w) # Original total sum of weights

  # Filter U by removing unwanted phenotypes (conditions)
  U <- lapply(U, function(mat, toKeep) {
    missingConditions <- setdiff(toKeep, colnames(mat))
    if (length(missingConditions) > 0) {
      stop(paste("Condition(s)", paste(missingConditions,
        collapse = ", "
      ), "not found in matrix"))
    }
    mat[toKeep, toKeep] # Keep only relevant conditions
  }, conditionsToKeep)

  # Remove matrices where all values are zero or weight is below cutoff
  keepNames <- names(keep(U, function(mat) !all(mat == 0)))
  if (!is.null(w)) {
    keepNames <- intersect(keepNames, names(w[w >= wCutoff]))
  }
  U <- U[keepNames]
  if (!is.null(w)) {
    w <- w[keepNames]
  }

  # Note: Matrices in U may contain very small values on the diagonal
  # even when contexts are not present, due to EM algorithm adjustments.
  # This makes the matrix not exactly zero, so it won't be removed even though it may not
  # have strong context relevance. This behavior arises because the algorithm attempts to
  # ensure matrices are full-rank, slightly changing initial values.

  # We cannot simply remove diagonal matrices as signals on the diagonal can be strong and relevant.
  # So we manually remove the U components that are driven by non-relevant contexts.
  U[conditionsToFilter] <- NULL
  w <- w[!names(w) %in% conditionsToFilter]

  # Recalculate the sum of remaining weights
  sumWnew <- sum(w)

  # Adjust weights to maintain the original sumW
  w <- (w / sumWnew) * sumW

  message(paste(length(U), "components of matrices remained after filtering."))

  return(list(U = U, w = w))
}

#' @importFrom purrr map_dfr
#' @importFrom dplyr bind_rows
mergeSusieCs <- function(susieFit, coverage = "CS_95_susie", method = NULL) {
  if (is.null(coverage)) coverage <- "CS_95_susie"
  coverage <- .translateLegacyCsColumnName(coverage)
  # Identify variant IDs that are associated with more than one credible set
  identifyOverlapSets <- function(variantsSetsAndPipsList) {
    overlapSets <- list()
    for (variantId in names(variantsSetsAndPipsList)) {
      sets <- variantsSetsAndPipsList[[variantId]][["sets"]]
      if (length(sets) > 1) {
        overlapSets[[variantId]] <- sets
      }
    }
    return(overlapSets)
  }
  # Merge overlapping credible sets using connected components.
  mergeAndUpdateOverlapSets <- function(variantsSetsAndPipsList, overlapSets) {
    allSets <- unique(unlist(overlapSets))
    if (length(allSets) == 0) return(list())

    parent <- setNames(allSets, allSets)
    findRoot <- function(x) {
      while (!identical(parent[[x]], x)) x <- parent[[x]]
      x
    }
    unionSets <- function(a, b) {
      rootA <- findRoot(a)
      rootB <- findRoot(b)
      if (!identical(rootA, rootB)) parent[[rootB]] <<- rootA
    }

    for (sets in overlapSets) {
      if (length(sets) > 1) {
        for (s in sets[-1]) unionSets(sets[[1]], s)
      }
    }

    components <- split(names(parent), vapply(names(parent), findRoot, character(1)))
    setNameMap <- list()
    for (members in components) {
      label <- paste(sort(members), collapse = ",")
      for (s in members) {
        setNameMap[[s]] <- label
      }
    }

    # Update each variant's credible set names
    updatedCredibleSets <- lapply(
      setNames(names(variantsSetsAndPipsList), names(variantsSetsAndPipsList)),
      function(variantId) {
        currentSets <- variantsSetsAndPipsList[[variantId]][["sets"]]
        mapped <- intersect(currentSets, names(setNameMap))
        if (length(mapped) > 0) {
          setNameMap[[mapped[1]]]
        } else {
          paste(sort(unique(currentSets)), collapse = ",")
        }
      }
    )
    return(updatedCredibleSets)
  }
  # Loop through each condition and their credible sets
  extractTopLoci <- function(susieFit, coverage) {
    # Build a flat data frame of (variant_id, pip, set_name) across all conditions
    condNames <- names(susieFit[[1]])
    rows <- map_dfr(seq_along(condNames), function(i) {
      condData <- susieFit[[1]][[i]]
      topLoci <- .translateLegacyTopLociCsColumns(condData[["top_loci"]])
      if (is.null(topLoci) || nrow(topLoci) == 0) return(NULL)
      pipCol <- resolvePipColumn(topLoci, method)
      if (is.null(pipCol)) return(NULL)

      setNum <- unique(topLoci[[coverage]])
      setNum <- setNum[!is.na(setNum) & setNum != 0]
      if (length(setNum) == 0) return(NULL)

      map_dfr(setNum, function(sn) {
        rows <- topLoci[topLoci[[coverage]] == sn & !is.na(topLoci[[coverage]]),
                         c("variant_id", pipCol), drop = FALSE]
        names(rows)[names(rows) == pipCol] <- "pip"
        rows$set_name <- paste0("cs_", i, "_", sn)
        rows
      })
    })

    if (is.null(rows) || nrow(rows) == 0) return(list())

    # Aggregate by variant_id preserving first-seen order
    seenOrder <- unique(rows$variant_id)
    splitRows <- split(rows, factor(rows$variant_id, levels = seenOrder))
    lapply(splitRows, function(df) {
      list(sets = df$set_name, pips = df$pip)
    })
  }

  combineTopLoci <- function(extractedResult) {
    if (length(extractedResult) == 0) return(NULL)

    # Compute overlap sets once, outside the per-variant loop
    overlapSets <- identifyOverlapSets(extractedResult)
    hasOverlaps <- length(overlapSets) != 0
    mergedSets <- if (hasOverlaps) {
      mergeAndUpdateOverlapSets(extractedResult, overlapSets = overlapSets)
    } else {
      NULL
    }

    topLociDf <- do.call(rbind, lapply(names(extractedResult), function(variantId) {
      maxPip <- max(unlist(extractedResult[[variantId]]$pips))
      medianPip <- median(unlist(extractedResult[[variantId]]$pips))
      credibleSetNames <- if (hasOverlaps) {
        mergedSets[[variantId]]
      } else {
        paste(sort(unique(unlist(extractedResult[[variantId]]$sets))), collapse = ",")
      }
      data.frame(
        variant_id = variantId, credibleSetNames = credibleSetNames,
        maxPip = maxPip, medianPip = medianPip, stringsAsFactors = FALSE
      )
    }))
    return(topLociDf)
  }

  extractedTopLoci <- extractTopLoci(susieFit, coverage = coverage)
  if (length(extractedTopLoci) == 0) return(NULL)
  combinedTopLociDf <- combineTopLoci(extractedTopLoci)
  if (is.null(combinedTopLociDf) || nrow(combinedTopLociDf) == 0) return(NULL)
  # Clean up row names and make sure variant_id is unique
  combinedTopLociDf <- combinedTopLociDf[!duplicated(combinedTopLociDf$variant_id), ]
  rownames(combinedTopLociDf) <- NULL # Clean up row names
  return(combinedTopLociDf)
}


#' @export
mashRandNullSample <- function(dat, nRandom, nNull, excludeCondition, seed = NULL) {
  # Function to extract one data set
  extractOneData <- function(dat, nRandom, nNull) {
    if (is.null(dat)) {
      return(NULL)
    }

    if ("z" %in% names(dat)) {
      absZ <- abs(dat$z)
      zData <- dat$z
    } else {
      absZ <- abs(dat$bhat / dat$sbhat)
      zData <- NULL
    }

    sampleIdx <- 1:nrow(absZ)
    randomIdx <- sample(sampleIdx, min(nRandom, length(sampleIdx)), replace = FALSE)

    if (!is.null(zData)) {
      random <- list(z = zData[randomIdx, , drop = FALSE])
    } else {
      random <- list(
        bhat = dat$bhat[randomIdx, , drop = FALSE],
        sbhat = dat$sbhat[randomIdx, , drop = FALSE]
      )
    }

    null.id <- which(apply(absZ, 1, max) < 2)
    if (length(null.id) == 0) {
      warning(paste("no variants are included in the null dataset because absZ > 2 for all variants in", dat$region))
      null <- list()
    } else {
      if (length(null.id) < ncol(absZ)) {
        warning(paste("not enough null data to estimate null correlation in", dat$region))
        null <- list()
      } else {
        nullIdx <- sample(null.id, min(nNull, length(null.id)), replace = FALSE)
        if (!is.null(zData)) {
          null <- list(z = zData[nullIdx, , drop = FALSE])
        } else {
          null <- list(
            bhat = dat$bhat[nullIdx, , drop = FALSE],
            sbhat = dat$sbhat[nullIdx, , drop = FALSE]
          )
        }
      }
    }
    dat <- list(random = random, null = null)
    return(dat)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (length(excludeCondition) > 0) {
    colsToCheck <- if ("z" %in% names(dat)) "z" else "bhat"
    if (!all(excludeCondition %in% colnames(dat[[colsToCheck]]))) {
      stop(paste("Error: excludeCondition are not present in", dat$region))
    }
    for (key in intersect(names(dat), c("z", "bhat", "sbhat"))) {
      dat[[key]] <- dat[[key]][, -excludeCondition, drop = FALSE]
    }
  }

  result <- extractOneData(dat, nRandom, nNull)
  return(result)
}

#' @export
mergeMashData <- function(resData, oneData) {
  if (length(resData) == 0 || is.null(resData)) return(oneData)
  if (length(oneData) == 0 || is.null(oneData)) return(resData)

  combinedData <- lapply(names(oneData), function(d) {
    od <- oneData[[d]]
    rd <- resData[[d]]
    if (length(od) == 0 || is.null(od)) return(rd)
    if (is.null(rd) || length(rd) == 0) return(od)

    # bind_rows auto-aligns columns, filling missing with NA; replace with NaN
    rnRes <- rownames(as.data.frame(rd))
    rnOne <- rownames(as.data.frame(od))
    combined <- bind_rows(as.data.frame(rd), as.data.frame(od))
    combined[is.na(combined)] <- NaN
    rnAll <- make.names(c(rnRes, rnOne), unique = TRUE)
    rownames(combined) <- rnAll
    combined
  })
  names(combinedData) <- names(oneData)
  return(combinedData)
}

# Internal: convert a single SumStats object (post-QC) into a (Bhat, Shat)
# pair of matrices keyed by context. Each row of the matrix corresponds to
# one (variantId × (study, trait)) cell from the SumStats entries; each
# column corresponds to a context (from QtlSumStats $context; from
# GwasSumStats $study, which is the per-study mash column).
#
# For QtlSumStats:
#   * Pivots entries on (study, trait) so each (study, trait) becomes a
#     block of rows and each context becomes a column. Missing
#     (study, trait, context) cells are filled with NA.
# For GwasSumStats:
#   * Each row of the collection is one study; we treat each study as a
#     mash "context" (single block of rows per study, columns = studies).
#     This is rarely used on its own but lets a flat GwasSumStats pass
#     through alongside (or instead of) a QtlSumStats without special
#     casing further upstream.
#
# Variant alignment within a (study, trait) block uses the entry's
# variant order; missing variants in any one context are filled with NA.
# NA in Bhat is mapped to 0 and NA in Shat is mapped to a large value
# (1000) inside mashr::mash_set_data via its `zero_Bhat_Shat_reset`
# pathway, matching the prior pipeline's handling of incomplete cells.
# @noRd
.mashSumStatsToMatrices <- function(x, role,
                                    inputScale = c("auto", "beta", "z")) {
  inputScale <- match.arg(inputScale)
  if (!methods::is(x, "QtlSumStats") && !methods::is(x, "GwasSumStats")) {
    stop(sprintf(
      "mashPipeline: '%s' input must be a QtlSumStats or GwasSumStats; got %s.",
      role, paste(class(x), collapse = "/")))
  }
  if (length(getQcInfo(x)) == 0L) {
    stop(sprintf(
      "mashPipeline: '%s' SumStats has no QC info (length(getQcInfo(x)) == 0L). ",
      role),
      "Run summaryStatsQc() on the SumStats before passing it to mashPipeline().")
  }
  if (nrow(x) == 0L) {
    stop(sprintf(
      "mashPipeline: '%s' SumStats has no entries (nrow == 0).", role))
  }

  isQtl <- methods::is(x, "QtlSumStats")
  if (isQtl) {
    studyCol <- as.character(x$study)
    traitCol <- as.character(x$trait)
    contextCol <- as.character(x$context)
    blockKeys <- paste(studyCol, traitCol, sep = "::")
    columnLabels <- unique(contextCol)
  } else {
    studyCol <- as.character(x$study)
    blockKeys <- studyCol
    contextCol <- studyCol
    columnLabels <- unique(studyCol)
  }

  # Per-entry GRanges (avoid `@`; use the public list-column accessor).
  entries <- x$entry

  # Resolve per-entry scale: which (Bhat, Shat) source to pull. mashr
  # expects one coherent convention per call:
  #   "beta" → Bhat = BETA, Shat = SE   (effect-size scale; standard)
  #   "z"    → Bhat = Z,    Shat = 1    (z-score scale)
  # "auto" picks "beta" when every entry has BETA + SE, else "z" if
  # every entry has Z. Mixed inputs (some entries missing BETA, others
  # missing Z) are a hard error.
  entryCaps <- lapply(seq_len(nrow(x)), function(i) {
    mc <- S4Vectors::mcols(entries[[i]])
    list(hasBetaSe = all(c("BETA", "SE") %in% colnames(mc)),
         hasZ      = "Z" %in% colnames(mc))
  })
  allHaveBetaSe <- all(vapply(entryCaps, `[[`, logical(1), "hasBetaSe"))
  allHaveZ      <- all(vapply(entryCaps, `[[`, logical(1), "hasZ"))
  resolvedScale <- switch(inputScale,
    beta = {
      if (!allHaveBetaSe)
        stop(sprintf(
          "mashPipeline: inputScale = 'beta' requires every '%s' entry to ",
          role),
          "carry both BETA and SE mcols.")
      "beta"
    },
    z = {
      if (!allHaveZ)
        stop(sprintf(
          "mashPipeline: inputScale = 'z' requires every '%s' entry to ",
          role), "carry a Z mcol.")
      "z"
    },
    auto = {
      if (allHaveBetaSe) "beta"
      else if (allHaveZ) "z"
      else stop(sprintf(
        "mashPipeline: '%s' SumStats has no usable scale — every entry ",
        role),
        "must carry (BETA, SE) or Z mcols.")
    })

  # Group rows of x by (study, trait) block; within each block, build a
  # variant × context matrix for Bhat and Shat.
  uniqBlocks <- unique(blockKeys)
  bhatBlocks <- vector("list", length(uniqBlocks))
  shatBlocks <- vector("list", length(uniqBlocks))
  variantRowNames <- vector("list", length(uniqBlocks))

  for (bi in seq_along(uniqBlocks)) {
    bkey <- uniqBlocks[[bi]]
    rowsInBlock <- which(blockKeys == bkey)

    # Variant universe for this block = union of variant IDs (SNP)
    # across the contexts in this block, preserving first-seen order.
    variantOrder <- character()
    perContextB  <- list()
    perContextSe <- list()
    requireCols <- if (resolvedScale == "beta") c("SNP", "BETA", "SE")
                   else                          c("SNP", "Z")
    for (rIdx in rowsInBlock) {
      df <- if (isQtl) {
        getSumstatDf(x,
                     study   = studyCol[[rIdx]],
                     context = contextCol[[rIdx]],
                     trait   = traitCol[[rIdx]],
                     require = requireCols)
      } else {
        getSumstatDf(x, study = studyCol[[rIdx]], require = requireCols)
      }
      snps <- df$variant_id
      newSnps <- setdiff(snps, variantOrder)
      variantOrder <- c(variantOrder, newSnps)
      ctx <- contextCol[[rIdx]]
      if (resolvedScale == "beta") {
        perContextB[[ctx]]  <- setNames(df$beta, snps)
        perContextSe[[ctx]] <- setNames(df$se,   snps)
      } else {
        perContextB[[ctx]]  <- setNames(df$z, snps)
        perContextSe[[ctx]] <- setNames(rep(1, length(snps)), snps)
      }
    }

    nVar <- length(variantOrder)
    bMat <- matrix(NA_real_, nrow = nVar, ncol = length(columnLabels),
                   dimnames = list(variantOrder, columnLabels))
    sMat <- matrix(NA_real_, nrow = nVar, ncol = length(columnLabels),
                   dimnames = list(variantOrder, columnLabels))
    for (ctx in names(perContextB)) {
      bMat[names(perContextB[[ctx]]), ctx] <- perContextB[[ctx]]
      sMat[names(perContextSe[[ctx]]), ctx] <- perContextSe[[ctx]]
    }
    # Disambiguate rownames across blocks to avoid silent dedup.
    rownames(bMat) <- paste(bkey, variantOrder, sep = "::")
    rownames(sMat) <- rownames(bMat)
    bhatBlocks[[bi]] <- bMat
    shatBlocks[[bi]] <- sMat
    variantRowNames[[bi]] <- rownames(bMat)
  }

  bhat <- do.call(rbind, bhatBlocks)
  shat <- do.call(rbind, shatBlocks)

  # Replace NAs: bhat NA -> 0, shat NA -> 1000 (the same convention as
  # filterInvalidSummaryStat() and mash_set_data()'s
  # zero_Bhat_Shat_reset; ensures missing-cell variants do not drive the
  # fit).
  bhat[is.na(bhat)] <- 0
  shat[is.na(shat) | shat <= 0] <- 1000

  list(b = bhat, s = shat)
}

#' Estimate mash covariance matrices and mixture weights from SumStats.
#'
#' Genome-wide pipeline that estimates the \pkg{mashr} canonical, PCA,
#' flash, and ED covariance matrices plus mixture weights between
#' contexts. The pipeline is memory-intensive and not gene-parallelizable
#' (see \code{dev/refactor-design.md}).
#'
#' @param sumStatsList A named \code{list} (or \code{S4Vectors::SimpleList})
#'   of \code{\link{QtlSumStats}} and/or \code{\link{GwasSumStats}}
#'   objects. The list MUST be named with at least \code{"strong"} and
#'   \code{"random"}; \code{"null"} is optional. Each element is one
#'   SumStats collection whose entries are pivoted internally into a
#'   variants \eqn{\times} contexts \eqn{Bhat} / \eqn{Shat} matrix pair.
#'   For \code{QtlSumStats} the columns of the resulting matrix are the
#'   \code{context} values; for \code{GwasSumStats} the columns are the
#'   \code{study} values. Every SumStats must have been processed by
#'   \code{\link{summaryStatsQc}} (the pipeline rejects inputs where
#'   \code{length(getQcInfo(x)) == 0L}).
#' @param alpha mash \code{alpha} parameter (passed to
#'   \code{mashr::mash_set_data()}).
#' @param residualCorrelation Optional residual correlation matrix. Used
#'   in place of the null-data-derived V matrix when no \code{"null"}
#'   entry is supplied (matches the prior pipeline).
#' @param nPcs Number of principal components for PCA-based covariance
#'   matrices. Defaults to \code{ncol(Bhat) - 1}.
#' @param setSeed Integer seed for reproducibility (default 999).
#'
#' @return A list with two elements: \code{U} (the combined list of mash
#'   covariance matrices: canonical + PCA + flash + ED) and \code{w}
#'   (the estimated mixture weights from \code{mashr::get_estimated_pi}).
#'
#' @section Behavioural notes (changes vs. the legacy contract):
#' The legacy \code{mashInput} list of six pre-built matrices
#' (\code{strong.b}/\code{strong.s}/\code{random.b}/\code{random.s}/
#' \code{null.b}/\code{null.s}) is replaced by a \code{list} of SumStats
#' objects. Construction of the per-partition matrices now happens
#' inside \code{mashPipeline} via \code{getSumStats()} /
#' \code{S4Vectors::mcols()} accessors — no \code{@@} slot access. The
#' mashr algorithm itself is unchanged: same \code{cov_canonical},
#' \code{cov_pca}, \code{cov_flash}, \code{cov_ed}, and
#' \code{mash(..., outputlevel = 1)} sequence as before.
#'
#' @export
