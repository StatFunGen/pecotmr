mashPipeline <- function(sumStatsList, alpha,
                          residualCorrelation = NULL,
                          nPcs = NULL,
                          setSeed = 999) {
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

  strongMats <- .mashSumStatsToMatrices(sumStatsList$strong, "strong")
  randomMats <- .mashSumStatsToMatrices(sumStatsList$random, "random")

  hasNull <- "null" %in% names(sumStatsList) && !is.null(sumStatsList$null)
  if (hasNull) {
    nullMats <- .mashSumStatsToMatrices(sumStatsList$null, "null")
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
              flippedData <- .matchRefPanel(data, ldBimFile$V2,
                colToFlip = c(valueColumn),
                matchMinProp = 0,
                flipStrand = FALSE, removeUnmatched = TRUE
              )$harmonizedData
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

            flippedData <- .matchRefPanel(cohortDf, refPanel, colToFlip = c(valueColumn),
                matchMinProp = 0,
                flipStrand = FALSE, removeUnmatched = TRUE)$harmonizedData

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
      hasFm <- !is.null(element$finemappingEntry) && is(element$finemappingEntry, "FineMappingEntry")
      hasSumstats <- "sumstats" %in% names(element)
      if (hasSumstats && hasFm) {
        variantNames <- getVariantIds(element$finemappingEntry)
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
