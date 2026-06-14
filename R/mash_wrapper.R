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
# This function extracts tensorQTL results for given region for multiple
# summary statistics files
#' @export
loadMultitraitTensorqtlSumstat <- function(
    sumstatsPaths, region, gene = NULL,
    traitNames = NULL, topLoci = FALSE, filterFile = NULL, removeAnyMissing = FALSE,
    maxRowsSelected = 300, nanRemove = TRUE) {
  if (!is.vector(sumstatsPaths) || !all(file.exists(sumstatsPaths))) {
    stop("sumstatsPaths must be a vector of existing file paths.")
  }
  if (!is.character(region) || length(region) != 1) {
    stop("region must be a single character string.")
  }
  if (!is.character(traitNames)) {
    stop("traitNames must be a vector of character strings.")
  }


  extractTensorqtlData <- function(path, region) {
    tabixRegion(path, region) %>%
      # first four columns are 'chrom','pos','alt','ref'
      mutate(variants = paste(.[[1]], .[[2]], .[[3]], .[[4]], sep = ":")) %>%
      distinct(variants, molecular_trait_id, .keep_all = TRUE)
  }


  mergeMatrices <- function(matrixList, valueColumn, idColumn = "variants",
                             removeAnyMissing = FALSE) {
    # Convert matrices to data frames
    dfList <- lapply(seq_along(matrixList), function(i) {
      df <- as.data.frame(matrixList[[i]])
      df2 <- df[, c(idColumn, valueColumn)]
      # Rename columns to avoid duplication
      colnames(df2) <- c(idColumn, paste0(valueColumn, "_", i))
      return(df2)
    })

    # Iteratively merge the data frames
    mergedDf <- Reduce(
      function(x, y) merge(x, y, by = idColumn, all = TRUE),
      dfList
    )

    # Optionally, remove rows with any missing values
    if (removeAnyMissing) {
      mergedDf <- mergedDf[complete.cases(mergedDf), ]
    }
    return(mergedDf)
  }

  splitVariantsAndMatch <- function(variant, filterFile, maxRowsSelected) {
    if (!file.exists(filterFile)) {
      stop("Filter file does not exist.")
    }

    # Split the variant vector into components
    variantSplit <- strsplit(variant, ":")
    variantDf <- data.frame(chr = sapply(variantSplit, `[`, 1), pos = sapply(
      variantSplit,
      `[`, 2
    ), stringsAsFactors = FALSE)
    variantDf$pos <- as.numeric(variantDf$pos)

    # get the region of interest
    minPos <- min(variantDf$pos)
    maxPos <- max(variantDf$pos)
    chrom <- unique(variantDf$chr)
    if (length(chrom) != 1) {
      stop("Variants are from multiple chromosomes. Cannot create a single range string.")
    }
    region <- paste0(chrom, ":", minPos, "-", maxPos)
    refTable <- tabixRegion(filterFile, region)
    if (is.null(refTable)) {
      stop("No variants in the region.")
    }
    colnames(refTable)[1:2] <- c("#CHROM", "POS")
    if (!all(c("#CHROM", "POS") %in% colnames(refTable))) {
      stop("Filter file must contain columns: #CHROM, POS.")
    }
    matchedIndices <- which(variantDf$chr %in% refTable$`#CHROM` & variantDf$pos %in%
      refTable$POS)
    if (!is.null(maxRowsSelected) && maxRowsSelected > 0 && maxRowsSelected <
      length(matchedIndices)) {
      selectedRows <- sample(length(matchedIndices), maxRowsSelected)
      matchedIndices <- matchedIndices[selectedRows]
    }
    return(matchedIndices)
  }

  Y <- lapply(sumstatsPaths, function(x, region, gene) {
    out <- extractTensorqtlData(x, region)
    if (!is.null(gene)) {
      out <- out[which(out$molecular_trait_id %in% gene), ]
    }
    sorted <- out[order(-abs(out$beta / out$se)), c("variants", "molecular_trait_id")]
    topV <- apply(sorted[1:2, ], 1, function(row) paste(row, collapse = "_")) # paste the variant (chr:pos:alt:ref) with gene_id with '_'
    out <- as.list(out)
    out$topVariants <- topV
    return(out)
  }, region = region, gene = gene)

  ## Y is list of data frames where colnames(Y[[1]])=
  ## c('chrom','pos','alt','ref','variant_id','molecular_trait_id','start_distance',
  ## 'end_distance','af','ma_samples','ma_count','pvalue','beta','se','molecular_trait_object_id',
  ## 'n','variant''topVariants')

  ### The step below assigns condition names to the Y; in case the filename
  ### itself does not contain any condition names, users can input the
  ### condition names via assigning traitNames if traitNames left blank, then
  ### the file names will be assigned as traitNames, where if the text before
  ### '.' (extension) can differentiate the conditions, we will use the shorter
  ### names as trait_name
  if (is.null(traitNames)) {
    traitNames <- gsub("\\..*", "", basename(sumstatsPaths)) # extract condition name that is listed before the first appearance of '.'

    if (length(traitNames[duplicated(traitNames)]) >= 1) {
      traitNames <- basename(sumstatsPaths)
    }
  }
  names(Y) <- traitNames

  bhat <- mergeMatrices(Y,
    valueColumn = "beta", idColumn = c("variants", "molecular_trait_id"),
    removeAnyMissing
  )
  sbhat <- mergeMatrices(Y,
    valueColumn = "se", idColumn = c("variants", "molecular_trait_id"),
    removeAnyMissing
  )
  out <- list(bhat = bhat, sbhat = sbhat)

  # Check if variants are the same in both bhat and sbhat
  if (!identical(out$bhat$variants, out$sbhat$variants)) {
    stop("Error: Variants in bhat and sbhat are not the same.")
  }

  varIdx <- 1:nrow(out$bhat)

  # match with filterFile
  if (!is.null(filterFile)) {
    if (!file.exists(filterFile)) {
      stop("Filter file does not exist.")
    }
    variants <- paste0(out$bhat$variants, "_", out$bhat$molecular_trait_id)
    varIdx <- splitVariantsAndMatch(out$bhat$variants, filterFile, maxRowsSelected)
  }

  if (topLoci) {
    unionTopLoci <- unique(unlist(lapply(Y, function(item) item$topVariants)))
    varIdx <- which(variants %in% unionTopLoci) # varIdx may end up empty if maxRowsSelected number too small
  }

  # Extract only subset of data
  variants <- paste0(out$bhat$variants[varIdx], "_", out$bhat$molecular_trait_id[varIdx])
  out$bhat <- out$bhat[varIdx, ]
  out$sbhat <- out$sbhat[varIdx, ]

  if (nanRemove) {
    out <- filterInvalidSummaryStat(out, bhat = "beta", sbhat = "se", nanRemove)
  }

  rownames(out$bhat) <- rownames(out$sbhat) <- variants
  colnames(out$bhat)[which(startsWith(colnames(out$bhat), "beta"))] <- colnames(out$sbhat)[which(startsWith(
    colnames(out$sbhat),
    "se"
  ))] <- traitNames
  out$region <- region
  out$topVariants <- lapply(Y, function(x) x$topVariants)
  return(out)
}

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
loadMultitraitRSumstat <- function(susieFit, sumstatsDb, coverage = NULL, extractInf = "z", topLoci = FALSE, filterFile = NULL, excludeCondition = NULL, ldMetaFile = NULL, removeAnyMissing = TRUE, maxRowsSelected = 300, nanRemove = FALSE, conditionFilter = FALSE) {
  # Internal recursive filtering function
  filterNestedList <- function(inputList, validConditions) {
    if (is.null(validConditions) || length(validConditions) == 0) {
      return(inputList)
    }

    filteredList <- list()

    # Recursively process list
    for (name in names(inputList)) {
      # Check if current name matches any valid condition
      if (name %in% validConditions) {
        filteredList[[name]] <- inputList[[name]]
      } else {
        # Recursively check nested lists
        if (is.list(inputList[[name]])) {
          nestedResult <- filterNestedList(inputList[[name]], validConditions)
          if (length(nestedResult) > 0) {
            filteredList[[name]] <- nestedResult
          }
        }
      }
    }

    return(filteredList)
  }

  # Apply condition filtering before processing
  if (!is.null(conditionFilter)) {
    # Convert conditionFilter to character vector if it's not already
    if (!is.character(conditionFilter)) {
      conditionFilter <- as.character(conditionFilter)
    }

    # Split if comma-separated string
    if (length(conditionFilter) == 1 && grepl(",", conditionFilter)) {
      conditionFilter <- strsplit(conditionFilter, ",")[[1]]
    }

    # Trim whitespace
    conditionFilter <- trimws(conditionFilter)

    # Filter susieFit and sumstatsDb
    susieFit <- filterNestedList(susieFit, conditionFilter)
    sumstatsDb <- filterNestedList(sumstatsDb, conditionFilter)
  }


  splitVariantsAndMatch <- function(variant, filterFile, maxRowsSelected) {
    if (!file.exists(filterFile)) {
      stop("Filter file does not exist.")
    }

    # Split the variant vector into components
    variantDf <- parse_variant_id(variant)
    conv <- attr(variantDf, "convention")

    # get the region of interest
    minPos <- min(variantDf$pos)
    maxPos <- max(variantDf$pos)
    chrom <- unique(variantDf$chrom)
    if (length(chrom) != 1) {
      stop("Variants are from multiple chromosomes. Cannot create a single range string.")
    }
    # Reconstruct chrom string with original convention for tabix
    chromStr <- if (conv$hasChr) paste0("chr", chrom) else as.character(chrom)
    region <- paste0(chromStr, ":", minPos, "-", maxPos)
    refTable <- tabixRegion(filterFile, region)
    if (is.null(refTable)) {
      stop("No variants in the region.")
    }
    colnames(refTable)[1:2] <- c("#CHROM", "POS")
    if (!all(c("#CHROM", "POS") %in% colnames(refTable))) {
      stop("Filter file must contain columns: #CHROM, POS.")
    }
    refChrom <- as.integer(stripChrPrefix(refTable$`#CHROM`))
    matchedIndices <- which(variantDf$chrom %in% refChrom & variantDf$pos %in% refTable$POS)
    if (!is.null(maxRowsSelected) && maxRowsSelected > 0 && maxRowsSelected < length(matchedIndices)) {
      selectedRows <- sample(length(matchedIndices), maxRowsSelected)
      matchedIndices <- matchedIndices[selectedRows]
    }
    return(matchedIndices)
  }


  results <- lapply(sumstatsDb[[1]], function(data) extractFlattenSumstatsFromNested(data, extractInf))
  traitNames <- names(results)
  zScores <- mergeMatrices(results, valueColumn = extractInf, ldMetaFile, idColumn = "variants", removeAnyMissing)
  out <- list(z = zScores)
  varIdx <- 1:nrow(out[[1]])
  if (!is.null(filterFile)) {
    variants <- out[[1]]$variants
    varIdx <- splitVariantsAndMatch(variants, filterFile, maxRowsSelected)
  }

  if (topLoci) {
    unionTopLoci <- mergeSusieCs(susieFit, coverage)
    if (!is.null(unionTopLoci)) {
      strongSignalDf <- unionTopLoci %>%
        group_by(credibleSetNames) %>%
        filter(medianPip == max(medianPip)) %>%
        slice(1) %>%
        ungroup()
      varIdx <- which(out[[1]]$variants %in% strongSignalDf$variant_id)
    } else {
      varIdx <- NULL
    }
  }

  # Extract only subset of data
  variants <- out[[1]]$variants[varIdx]
  for (key in names(out)) {
    out[[key]] <- out[[key]][varIdx, , drop = FALSE]
    rownames(out[[key]]) <- variants
    colnames(out[[key]])[2:ncol(out[[key]])] <- traitNames
    out[[key]] <- out[[key]][, -which(names(out[[key]]) == "variants"), drop = FALSE]
  }
  out$region <- names(susieFit)

  if (!is.null(excludeCondition) && length(excludeCondition) != 0) {
    for (key in setdiff(names(out), "region")) {
      if (all(excludeCondition %in% colnames(out[[key]]))) {
        out[[key]] <- out[[key]][, -excludeCondition]
      } else {
        stop(paste("Error: excludeCondition are not present in", out$region))
      }
    }
  }

  return(out)
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

#' @export
mashPipeline <- function(mashInput, alpha, residualCorrelation = NULL, nPcs = NULL, setSeed = 999) {
  if (!requireNamespace("mashr", quietly = TRUE)) {
    stop("To use this function, please install mashr: https://cran.r-project.org/web/packages/mashr/index.html")
  }
  if (!requireNamespace("flashier", quietly = TRUE)) {
    stop("To use this function, please install flashier: https://github.com/willwerscheid/flashier")
  }
  set.seed(setSeed)
  if (length(mashInput$null.b) == 0 && length(mashInput$null.s) == 0) {
    if (!is.null(residualCorrelation)) {
      vhat <- residualCorrelation
    } else {
      conditionNum <- ncol(mashInput$random.b)
      vhat <- diag(rep(1, conditionNum))
    }
  } else {
    vhat <- mashr::estimate_null_correlation_simple(mashr::mash_set_data(mashInput$null.b,
      Shat = mashInput$null.s,
      alpha, zero_Bhat_Shat_reset = 1000
    ))
  }

  mashData <- mashr::mash_set_data(mashInput$strong.b,
    Shat = mashInput$strong.s, V = vhat,
    alpha, zero_Bhat_Shat_reset = 1000
  )

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

  return(list(U = U.all, w = w))
}

#' Merge a List of Matrices or Data Frames with Optional Allele Flipping
#'
#' @description
#' This function merges a list of matrices or data frames by a shared identifier column,
#' optionally aligning to a reference panel using allele QC procedures.
#'
#' @param matrixList A named or unnamed list of data frames or matrices.
#' @param valueColumn Character string. The name of the column containing values to extract (e.g., z-scores or betas).
#' @param refPanel Optional data frame. A reference panel for allele QC (must be compatible with `allele_qc`).
#' @param idColumn Character string. The name of the column identifying variant IDs. Default is `"variants"`.
#' @param removeAnyMissing Logical. If `TRUE`, rows with any missing values will be removed after merging.
#'
#' @return A data frame containing merged values, one column per dataset with suffix `_i`.
#' @examples
#' \dontrun{
#' merged <- mergeSumstatsMatrices(list(df1, df2), valueColumn = "variants", refPanel = ref_df)
#' }
#' @import dplyr
#' @export

mergeSumstatsMatrices <- function(matrixList, valueColumn, refPanel = NULL, ldMetaFile = NULL, idColumn = "variants",
                             removeAnyMissing = FALSE) {
    # Input validation
    if (!is.list(matrixList) || length(matrixList) == 0) {
      stop("matrixList must be a non-empty list")
    }
    if (!is.character(valueColumn) || length(valueColumn) != 1) {
      stop("valueColumn must be a single string")
    }
    if (!is.character(idColumn) || length(idColumn) != 1) {
      stop("idColumn must be a single string")
    }

    dfList <- lapply(seq_along(matrixList), function(i) {
      tryCatch(
        {
           # Step 1: Convert matrix to data frame and extract relevant columns
           df <- as.data.frame(matrixList[[i]])
           if (!(idColumn %in% colnames(df)) || !(valueColumn %in% colnames(df))) {
            stop(paste("Required columns", idColumn, "or", valueColumn, "not found in dataset", i))
           }
           df2 <- df[, c(idColumn, valueColumn)]
             if (!is.null(ldMetaFile)) {
            # Step 2: Split 'variants' to extract chromosomal info
            cohortVariantsDf <- parse_variant_id(df2[, c(idColumn)])
            # Step 3: Combine extracted chromosomal info with value column
            cohortDf <- cbind(cohortVariantsDf, value = df2[, valueColumn, drop = FALSE])

            # Step 4: Merge with LD reference and filter
            # Normalize ldMetaFile chrom to integer to match parse_variant_id output
            ldMetaFile$chrom <- as.integer(stripChrPrefix(as.character(ldMetaFile$chrom)))
            variantsLdBlockMatch <- merge(cohortDf, ldMetaFile, by = "chrom", allow.cartesian = TRUE) %>%
              filter(pos > start & pos < end) %>%
              select(-path)

            # Function to process each group
            processGroup <- function(data) {
              # Construct file path
              bimFilePath <- unique(data$bim_path)
              ldBimFile <- vroom(bimFilePath)

              # Perform allele quality control
              flippedData <- getHarmonizedData(matchRefPanel(data, ldBimFile$V2,
                colToFlip = c(valueColumn),
                matchMinProp = 0, removeDups = FALSE,
                removeIndels = FALSE, removeStrandAmbiguous = FALSE,
                flipStrand = FALSE, removeUnmatched = TRUE
              ))
              return(flippedData)
            }

            finalDf <- variantsLdBlockMatch %>%
              group_by(start, end) %>%
              group_map(~ processGroup(.x)) %>%
              bind_rows() %>%
              select(c("variant_id", valueColumn)) %>%
              rename("variants" = "variant_id")
            # Rename columns to avoid duplication
            colnames(finalDf) <- c(idColumn, paste0(valueColumn, "_", i))
          } else if (!is.null(refPanel)) {
            # Step 2: Split 'variants' to extract chromosomal info
            cohortVariantsDf <- parse_variant_id(df2[, c(idColumn)])
            # Step 3: Combine extracted chromosomal info with value column
            cohortDf <- cbind(cohortVariantsDf, value = df2[, valueColumn, drop = FALSE])

            flippedData <- getHarmonizedData(matchRefPanel(cohortDf, refPanel, colToFlip = c(valueColumn),
                matchMinProp = 0, removeDups = FALSE,
                removeIndels = FALSE, removeStrandAmbiguous = FALSE,
                flipStrand = FALSE, removeUnmatched = TRUE, removeSameVars = FALSE))

            finalDf <- flippedData %>%
                select(c("variant_id", valueColumn))
            colnames(finalDf) <- c(idColumn, paste0(valueColumn, "_", i))
          } else {
            finalDf <- df2
            colnames(finalDf) <- c(idColumn, paste0(valueColumn, "_", i))
          }
          return(finalDf)
        },
        error = function(e) {
          message(paste("Error processing dataset", i, ":", e$message))
          return(NULL)
        }
      )
    })

    # Remove any NULL results from errors
    dfList <- dfList[!sapply(dfList, is.null)]
    if (length(dfList) == 0) {
        message("No valid datasets after processing")
        return(NULL)
    }

    # Iteratively merge the data frames
    mergedDf <- Reduce(
      function(x, y) merge(x, y, by = idColumn, all = TRUE),
      dfList
    )
    # Optionally, remove rows with any missing values
    if (removeAnyMissing) {
      mergedDf <- mergedDf[complete.cases(mergedDf), ]
    }
    return(mergedDf)
  }
                 
#' Load and Align Summary Statistics for a Given Gene and Condition
#'
#' @description
#' This function processes summary statistics matrices for a target gene across contexts,
#' optionally aligning with a reference panel and updating an existing result list.
#'
#' @param datList A named list of matrices or data.frames, each element corresponding to a summary statistics type (e.g., z, beta).
#' @param signalDf A data.frame containing signal information including `variant_ID`, `gene_ID`, and `event_ID`.
#' @param cond Character. Condition type: "strong", "null", or "random".
#' @param region Character. Target gene ID.
#' @param extractInfs Character vector. Names of summary statistics to extract (e.g., `"z"`, `"beta"`).
#' @param tagPatterns Optional named pattern list used to classify context.
#' @param resultListFormat A nested list used as a running result container.
#'
#' @importFrom stringr str_detect str_remove_all
#' @importFrom rlang .data sym
#' @importFrom purrr keep map_dfr map_chr
#' @importFrom utils combn
#' @import dplyr tidyr tibble
#' @return The updated `resultListFormat` with processed results for the specified gene and condition.
#' @export
loadMulticontextSumstats <- function(datList, signalDf, cond, region, extractInfs = "z", tagPatterns = NULL, resultListFormat) {
  # Initialize output list
  out <- list()
  traitNames <- names(datList[[1]])
    if (cond == "strong" && region %in% signalDf$gene_ID){
  events <- signalDf %>% filter(gene_ID == region) %>% pull(event_ID) %>% unique()
  for (j in seq_along(events)){
        refDfFiltered <- signalDf %>% filter(gene_ID == region, event_ID == events[j]) %>%
            filter(!str_detect(context_classify, "NE"))
        if(dim(refDfFiltered)[1] == 0) next
        ## generate the reference panel for allele flipping
        refPanel <- parse_variant_id(refDfFiltered$variant_ID%>%unique())

        varIdx <- c()
        variants <- c()
        sumstatsDf <- list()
        eventIDextracted <- c()

        # Flatten the nested list
        for (extractInf in extractInfs) {
        extractedMatrix <- mergeSumstatsMatrices(datList[[extractInf]], valueColumn = extractInf, refPanel = refPanel, idColumn = "variants", removeAnyMissing = FALSE)
        if(is.null(extractedMatrix)||dim(extractedMatrix)[1]==0) return(resultListFormat)
        out[[extractInf]] <- extractedMatrix
        # Set variant order on first iteration
        if (is.null(varIdx)&& is.null(variants)) {
            varIdx <- 1:nrow(out[[extractInf]])
            variants <- out[[extractInf]]$variants[varIdx]
        }
        numberIndex <- str_extract(colnames(out[[extractInf]]), "\\d+")[-1]
        out[[extractInf]] <- out[[extractInf]][varIdx, , drop = FALSE]
        rownames(out[[extractInf]]) <- variants
        colnames(out[[extractInf]])[2:ncol(out[[extractInf]])] <- traitNames[as.integer(numberIndex)]
        out[[extractInf]] <- out[[extractInf]][, -which(names(out[[extractInf]]) == "variants"), drop = FALSE]

        df <- as.data.frame(t(out[[extractInf]]))
        df <- rownames_to_column(df, var = "context")

            # Match context to tag
        df <- df %>%
                  mutate(context_classify = if (is.null(tagPatterns) || length(tagPatterns) == 0) {
                    context
                  } else {
                    map_chr(context, function(ctx) {
                      matched <- names(tagPatterns)[str_detect(ctx, tagPatterns)]
                      if (length(matched) == 0) NA_character_ else matched[1]
                    })
                  })

        numericCol <- colnames(df)[2]

         if (extractInf == "z"){
                # Make a copy to store added rows
              addedDf <- data.frame()

                # Ensure the column name of the numeric column
                if (any(grepl("sQTL|pQTL|gpQTL", df$context_classify))) {
                  if (any(grepl("sQTL|pQTL|gpQTL", refDfFiltered$context_classify))) {

                    # Extract sQTL contexts to loop over
                    xQTLspecificContexts <- unique(str_subset(refDfFiltered$context_classify, "sQTL|pQTL|gpQTL"))

                    for (cont in xQTLspecificContexts) {
                          eventIDsExtracted <- refDfFiltered %>%
                                    filter(context_classify == cont) %>%
                                    pull(event_IDs)

                    # Filter matching rows in df
                          contextRows <- df %>%
                                filter(context_classify == cont, str_detect(context, paste(eventIDsExtracted, collapse = "|")))

                      if (nrow(contextRows) > 0) {
                        # Get the row with median absolute value
                        absValues <- abs(contextRows[[numericCol]])
                        medianVal <- median(absValues, na.rm = TRUE)
                        medianIdx <- which.min(abs(absValues - medianVal))  # Closest to median
                        selectedDf <- contextRows[medianIdx, , drop = FALSE]

                        addedDf <- bind_rows(addedDf, selectedDf)
                        df <- df %>% filter(context_classify !=cont)
                      }
                    }
                    # Combine updated sQTL-specific rows back into df
                    df <- bind_rows(df, addedDf)
                  }
                }
                sumstatsDf[[extractInf]] <- df %>%
                  filter(!str_detect(context_classify, "NE") & context_classify != 'NA')%>%
                  group_by(context_classify) %>%
                  slice_min(order_by = abs(.data[[numericCol]] - median(abs(.data[[numericCol]]), na.rm = TRUE)), n = 1, with_ties = FALSE) %>%
                  ungroup()%>%
                  rename(!!numericCol := !!sym(numericCol))
                eventIDextracted <- sumstatsDf[[extractInf]]%>%pull(context)
             } else if (is.null(eventIDextracted)){
                    warning("Please provide 'z-score'")
             } else {
                 sumstatsDf[[extractInf]] <- df %>% filter(context%in%eventIDextracted)%>%
                                        rename(!!numericCol := !!sym(numericCol))
             }
             resultDf <- sumstatsDf[[extractInf]] %>%
                  select(-context) %>%
                  rename(value = !!sym(numericCol)) %>%
                  pivot_wider(names_from = context_classify, values_from = value) %>%
                  mutate(
                    variant_ID = numericCol,
                    gene_ID = region
                  ) %>%
                 select(variant_ID, gene_ID, everything())
                 resultListFormat[[cond]][[extractInf]]  <- resultListFormat[[cond]][[extractInf]]%>% rows_update(resultDf, by = c("variant_ID", "gene_ID"))
     }
  }
}
  # Handle "null" condition
  if (cond%in%c("null","random") && region %in% signalDf$gene_ID) {
    refDfFiltered <- signalDf %>% filter(gene_ID == region)
    refPanel <- parse_variant_id(refDfFiltered$variant_ID %>% unique())

    varIdx <- c()
    variants <- c()
    sumstatsDf <- list()
    eventIDextracted <- list()
    for (extractInf in extractInfs){
         # Flatten the nested list
         extractedMatrix <- mergeSumstatsMatrices(datList[[extractInf]], valueColumn = extractInf, refPanel = refPanel, idColumn = "variants", removeAnyMissing = FALSE)
          if (is.null(extractedMatrix)||dim(extractedMatrix)[1]==0) return(resultListFormat)
         out[[extractInf]] <- extractedMatrix
          # Set variant order on first iteration
          if (is.null(varIdx)&& is.null(variants)) {
                varIdx <- 1:nrow(out[[extractInf]])
                variants <- out[[extractInf]]$variants[varIdx]
          }
          numberIndex <- str_extract(colnames(out[[extractInf]]), "\\d+")[-1]
          out[[extractInf]] <- out[[extractInf]][varIdx, , drop = FALSE]
          rownames(out[[extractInf]]) <- variants
          colnames(out[[extractInf]])[2:ncol(out[[extractInf]])] <- traitNames[as.integer(numberIndex)]
          out[[extractInf]] <- out[[extractInf]][, -which(names(out[[extractInf]]) == "variants"), drop = FALSE]

          for (k in 1: dim(out[[extractInf]])[1]){
               df <- as.data.frame(t(out[[extractInf]][k,]))
               df <- rownames_to_column(df, var = "context")

              # Match context to tag
               df <- df %>%
                   mutate(context_classify = if (is.null(tagPatterns) || length(tagPatterns) == 0) {
                     context
                   } else {
                     map_chr(context, function(ctx) {
                       matched <- names(tagPatterns)[str_detect(ctx, tagPatterns)]
                       if (length(matched) == 0) NA_character_ else matched[1]
                     })
                   })

              numericCol <- colnames(df)[2]
            if (extractInf == "z"){
                sumstatsDf[[extractInf]] <- df %>%
                          filter(!str_detect(context_classify, "NE") & context_classify != 'NA')
                if(cond == "null"){
                  sumstatsDf[[extractInf]] <- sumstatsDf[[extractInf]] %>%
                        group_by(context_classify) %>%
                        filter(
                            !is.na(.data[[numericCol]]),
                            if (any(str_detect(context_classify, "sQTL|pQTL|gpQTL"))) {
                              abs(.data[[numericCol]]) < 2
                            } else {
                              TRUE
                            }
                          )%>%
                         slice_min(
                            order_by = abs(.data[[numericCol]] - median(abs(.data[[numericCol]]), na.rm = TRUE)),
                            n = 1,
                            with_ties = FALSE
                          ) %>%
                          ungroup() %>%
                          rename(!!numericCol := !!sym(numericCol))
                } else if (cond == "random") {
                    sumstatsDf[[extractInf]] <-  sumstatsDf[[extractInf]] %>%
                          group_by(context_classify) %>%
                          slice_min(
                            order_by = abs(.data[[numericCol]] - median(abs(.data[[numericCol]]), na.rm = TRUE)),
                            n = 1,
                            with_ties = FALSE
                          ) %>%
                          ungroup() %>%
                          rename(!!numericCol := !!sym(numericCol))
                }
                eventIDextracted[[k]] <- sumstatsDf[[extractInf]]%>%pull(context)
            }  else if (is.null(eventIDextracted)){
                    warning("Please provide 'z-score'")
            } else {
                sumstatsDf[[extractInf]] <- df %>% filter(context%in%eventIDextracted[[k]])%>%
                                        rename(!!numericCol := !!sym(numericCol))
            }
            resultDf <- sumstatsDf[[extractInf]] %>%
                  select(-context) %>%
                  rename(value = !!sym(numericCol)) %>%
                  pivot_wider(names_from = context_classify, values_from = value) %>%
                  mutate(
                    variant_ID = numericCol,
                    gene_ID = region
                  ) %>%
                 select(variant_ID, gene_ID, everything())
            resultListFormat[[cond]][[extractInf]]  <- resultListFormat[[cond]][[extractInf]] %>% rows_update(resultDf, by = c("variant_ID", "gene_ID"))
            }
          }
       }
     return(resultListFormat)
   }

            
#' Extract Summary Statistics from Nested Data Structure
#'
#' @description
#' Recursively searches a nested list to extract summary statistics (z, beta, or se)
#' using `variantNames` and `sumstats`. Computes `z` if needed from `betahat` and `sebetahat`.
#'
#' @param data A nested list structure potentially containing `variantNames` and `sumstats`.
#' @param extractInf Character. One of `"z"`, `"beta"`, or `"se"`.
#' @param maxDepth Integer. Maximum depth to search within the list. Default is 3.
#'
#' @return A data.frame with columns `variants` and the requested summary statistic.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- extractFlattenSumstatsFromNested(nestedListObject, extractInf = "z")
#' }

extractFlattenSumstatsFromNested <- function(data, extractInf = "z", maxDepth = 3) {
  # Validate input
  if (!extractInf %in% c("z", "beta", "se")) {
    stop("extractInf must be one of: 'z', 'beta', or 'se'")
  }

  # Internal recursive function
  findNested <- function(element, currentDepth = 0) {
    if (currentDepth >= maxDepth) {
      message("Maximum search depth reached. Could not find 'variantNames' and 'sumstats' together.")
      return(NULL)
    }

    if (is.list(element)) {
      hasFm <- !is.null(element$finemappingResult) && is(element$finemappingResult, "FineMappingResult")
      hasSumstats <- "sumstats" %in% names(element)
      if (hasSumstats && hasFm) {
        variantNames <- getVariantNames(element$finemappingResult)
        sumstats <- element$sumstats

        # Extract based on type
        resultColumn <- switch(
          extractInf,
          "z" = {
            if (all(c("betahat", "sebetahat") %in% names(sumstats))) {
              sumstats$betahat / sumstats$sebetahat
            } else if ("z" %in% names(sumstats)) {
              sumstats$z
            } else {
              message("Cannot compute z: missing 'betahat' and 'sebetahat', and 'z' not available.")
              return(NULL)
            }
          },
          "beta" = {
            if ("betahat" %in% names(sumstats)) {
              sumstats$betahat
            } else {
              message("Missing 'betahat' for beta extraction.")
              return(NULL)
            }
          },
          "se" = {
            if ("sebetahat" %in% names(sumstats)) {
              sumstats$sebetahat
            } else {
              message("Missing 'sebetahat' for se extraction.")
              return(NULL)
            }
          }
        )

        result <- data.frame(variants = variantNames)
        result[[extractInf]] <- resultColumn

        # Normalize variants to canonical format (with chr prefix)
        result$variants <- normalizeVariantId(result$variants)

        return(result)
      }

      # Recurse into nested elements
      for (name in names(element)) {
        result <- findNested(element[[name]], currentDepth + 1)
        if (!is.null(result)) {
          result$variants <- normalizeVariantId(result$variants)
          return(result)
        }
      }
    }

    return(NULL)
  }

  # Start search
  return(findNested(data))
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
