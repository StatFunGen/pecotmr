#' Extract the trimmed SuSiE fit from a finemapping pipeline result
#'
#' Returns the trimmed model fit underlying \code{con_data$finemappingResult}
#' (a \code{FineMappingResult} S4 object), or NULL if no fine-mapping result
#' is attached.
#'
#' @param conData List. The method-layer entry from a finemapping pipeline
#'   result, expected to carry \code{$finemappingResult} as a
#'   \code{FineMappingResult} object.
#' @return The trimmed fit (a list with \code{pip}, \code{sets}, etc.) or NULL.
#' @export
getSusieResult <- function(conData) {
  if (length(conData) == 0) return(NULL)
  fm <- conData$finemappingResult
  if (is.null(fm) || !is(fm, "FineMappingResult")) return(NULL)
  trimmed <- getTrimmedFit(fm)
  if (length(trimmed) == 0) return(NULL)
  trimmed
}

#' Process Credible Sets (CS) from Finemapping Results
#'
#' This function extracts and processes information for each Credible Set (CS) 
#' from finemapping results, typically obtained from a finemapping RDS file.
#'
#' @param conData List. The method layer data from a finemapping RDS file that is not empty.
#' @param csNames Character vector. Names of the Credible Sets, usually in the format "L_<number>".
#' @param topLociTable Data frame. The $top_loci layer data from the finemapping results.
#'
#' @return A data frame with one row per CS, containing the following columns:
#'   \item{cs_name}{Name of the Credible Set}
#'   \item{variants_per_cs}{Number of variants in the CS}
#'   \item{top_variant}{ID of the variant with the highest PIP in the CS}
#'   \item{top_variant_index}{Global index of the top variant}
#'   \item{top_pip}{Highest Posterior Inclusion Probability (PIP) in the CS}
#'   \item{top_z}{Z-score of the top variant}
#'   \item{p_value}{P-value calculated from the top Z-score}
#'   \item{cs_corr}{Pairwise correlations of other CSs in this RDS with the CS of 
#'     the current row, delimited by '|', if there is more than one CS in this RDS file}
#'
#' @details
#' This function is designed to be used only when there is at least one Credible Set 
#' in the finemapping results usually for a given study and block. It processes each CS, 
#' extracting key information such as the top variant, its statistics, and 
#' correlation information between multiple CS if available.
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#'
#' @export
extractCsInfo <- function(conData, csNames, topLociTable) {
  fm <- conData$finemappingResult
  trimmed <- getTrimmedFit(fm)
  variantNames <- getVariantNames(fm)
  results <- map(seq_along(csNames), function(i) {
    csName <- csNames[i]
    indices <- trimmed$sets$cs[[csName]]

    # Get variants for this CS using the full variant names list
    csVariants <- variantNames[indices]
    csData <- topLociTable[topLociTable$variant_id %in% csVariants, ]
    topRow <- which.max(csData$pip)

    topVariant <- csData$variant_id[topRow]
    # Find the global index of the top variant
    topVariantGlobalIndex <- which(variantNames == topVariant)
    topPip <- csData$pip[topRow]
    topZ <- csData$z[topRow]
    pValue <- zToPvalue(topZ)

    # Extract cs_corr
    csCorr <- if (length(csNames) > 1) {
      trimmed$cs_corr[i,]
    } else {
      NA  # Use NA for the second CS or when there's only one CS
    }

    # Return results for this CS as a one-row data.frame
    result <- tibble(
      cs_name = csName,
      variants_per_cs = length(csVariants),
      top_variant = topVariant,
      top_variant_index = topVariantGlobalIndex,
      top_pip = topPip,
      top_z = topZ,
      p_value = pValue,
      cs_corr = list(paste(csCorr, collapse = ","))  # list column if csCorr is a vector
    )
    return(result)
  })
  # Combine all tibbles into one data frame
  finalResult <- bind_rows(results)
  return(finalResult)
}

#' Extract Information for Top Variant from Finemapping Results
#'
#' This function extracts information about the variant with the highest Posterior 
#' Inclusion Probability (PIP) from finemapping results, typically used when no 
#' Credible Sets (CS) are identified in the analysis.
#'
#' @param conData List. The method layer data from a finemapping RDS file.
#'
#' @return A data frame with one row containing the following columns:
#'   \item{cs_name}{NA (as no CS is identified)}
#'   \item{variants_per_cs}{NA (as no CS is identified)}
#'   \item{top_variant}{ID of the variant with the highest PIP}
#'   \item{top_variant_index}{Index of the top variant in the original data}
#'   \item{top_pip}{Highest Posterior Inclusion Probability (PIP)}
#'   \item{top_z}{Z-score of the top variant}
#'   \item{p_value}{P-value calculated from the top Z-score}
#'   \item{cs_corr}{NA (as no CS correlation is available)}
#'
#' @details
#' This function is designed to be used when no Credible Sets are identified in 
#' the finemapping results, but information about the most significant variant 
#' is still desired. It identifies the variant with the highest PIP and extracts 
#' relevant statistical information.
#'
#' @note
#' This function is particularly useful for capturing information about potentially 
#' important variants that might be included in Credible Sets under different 
#' analysis parameters or lower coverage. It maintains a structure similar to 
#' the output of `extract_cs_info()` for consistency in downstream analyses.
#'
#' @seealso
#' \code{\link{extractCsInfo}} for processing when Credible Sets are present.
#'
#' @export
extractTopPipInfo <- function(conData) {
  fm <- conData$finemappingResult
  trimmed <- getTrimmedFit(fm)
  variantNames <- getVariantNames(fm)
  # Find the variant with the highest PIP
  topPipIndex <- which.max(trimmed$pip)
  topPip <- trimmed$pip[topPipIndex]
  topVariant <- variantNames[topPipIndex]
  topZ <- conData$sumstats$z[topPipIndex]
  pValue <- zToPvalue(topZ)

  list(
    cs_name = NA,
    variants_per_cs = NA,
    top_variant = topVariant,
    top_variant_index = topPipIndex,
    top_pip = topPip,
    top_z = topZ,
    p_value = pValue,
    cs_corr = NA  # or NULL
  )
}

#' Parse Credible Set Correlations from extractCsInfo() Output
#'
#' This function takes the output from `extractCsInfo()` and expands the `cs_corr` column
#' into multiple columns, preserving the original order of correlations. It also
#' calculates maximum and minimum correlation values for each Credible Set.
#'
#' @param df Data frame. The output from `extractCsInfo()` function,
#'           containing a `cs_corr` column with correlation information.
#'
#' @return A data frame with the original columns from the input, plus:
#'   \item{cs_corr_1, cs_corr_2, ...}{Individual correlation values, with column names
#'         based on their position in the original string}
#'   \item{cs_corr_max}{Maximum absolute correlation value (excluding 1)}
#'   \item{cs_corr_min}{Minimum absolute correlation value}
#'
#' @details
#' The function splits the `cs_corr` column, which typically contains correlation
#' values separated by '|', into individual columns. It preserves the order of
#' these correlations, allowing for easy interpretation in a matrix-like format.
#'
#' @note
#' - This function converts the input to a data frame if it isn't already one.
#' - It handles cases where correlation values might be missing or not in the expected format.
#' - The function assumes that correlation values of 1 represent self-correlations and excludes
#'   these when calculating max and min correlations.
#'
#' @export
parseCsCorr <- function(df) {
  # Ensure we work with a data frame
  df <- as.data.frame(df)

  extractCorrelations <- function(x) {
    # Early return if x is invalid
    if(is.na(x) || x == "" || is.null(x) || !grepl(",", as.character(x))) {
      return(list(values = numeric(0), max_corr = NA_real_, min_corr = NA_real_))
    }

    # Convert and filter values
    values <- as.numeric(unlist(strsplit(x, ",")))
    valuesFiltered <- abs(values[values != 1])

    # Return list with NA if no valid correlations
    list(
      values = values,
      max_corr = if(length(valuesFiltered) > 0) max(abs(valuesFiltered), na.rm = TRUE) else NA_real_,
      min_corr = if(length(valuesFiltered) > 0) min(abs(valuesFiltered), na.rm = TRUE) else NA_real_
    )
  }
  # Process correlations
  processedResults <- lapply(df$cs_corr, extractCorrelations)
  # If no valid results, add NA columns and return
  if(all(sapply(processedResults, function(x) length(x$values) == 0))) {
    df$cs_corr_max <- NA_real_
    df$cs_corr_min <- NA_real_
    return(df)
  }

  # Determine max number of correlations
  maxCorrCount <- max(sapply(processedResults, function(x) length(x$values)))

  # Create and add correlation columns
  colNames <- paste0("cs_corr_", 1:maxCorrCount)

  for(i in seq_along(colNames)) {
    df[[colNames[i]]] <- sapply(processedResults, function(x) {
      if(length(x$values) >= i) x$values[i] else NA_real_
    })
  }

  # Add max and min correlation columns
  df$cs_corr_max <- sapply(processedResults, `[[`, "max_corr")
  df$cs_corr_min <- sapply(processedResults, `[[`, "min_corr")

  return(df)
}

#' Process Credible Set Information and Determine Updating Strategy
#'
#' This function categorizes Credible Sets (CS) within a study block into different 
#' updating strategies based on their statistical properties and correlations.
#'
#' @param df Data frame. Contains information about Credible Sets for a specific study and block.
#' @param highCorrCols Character vector. Names of columns in df that represent high correlations.
#'
#' @return A modified data frame with additional columns attached to the diagnostic table:
#'   \item{top_cs}{Logical. TRUE for the CS with the highest absolute Z-score.}
#'   \item{tagged_cs}{Logical. TRUE for CS that are considered "tagged" based on p-value and correlation criteria.}
#'   \item{method}{Character. The determined updating strategy ("BVSR", "SER", or "BCR").}
#'
#' @details
#' This function performs the following steps:
#' 1. Identifies the top CS based on the highest absolute Z-score.
#' 2. Identifies tagged CS based on high p-value and high correlations.
#' 3. Counts total, tagged, and remaining CS.
#' 4. Determines the appropriate updating method based on these counts.
#'
#' The updating methods are:
#' - BVSR (Bayesian Variable Selection Regression): Used when there's only one CS or all CS are accounted for.
#' - SER (Single Effect Regression): Used when there are tagged CS but no remaining untagged CS.
#' - BCR (Bayesian Conditional Regression): Used when there are remaining untagged CS.
#'
#' @note
#' This function is part of a developing methodology for automatically handling 
#' finemapping results. The thresholds and criteria used (e.g., p-value > 1e-4 for tagging) 
#' are subject to refinement and may change in future versions.
#'
#' @importFrom dplyr case_when
#'
#' @export
autoDecision <- function(df, highCorrCols) {
  # Identify top_cs
  topCsIndex <- which.max(abs(df$top_z))
  df$top_cs <- FALSE
  df$top_cs[topCsIndex] <- TRUE

  # Identify tagged_cs
  df$tagged_cs <- sapply(1:nrow(df), function(i) {
    if (df$top_cs[i]) return(FALSE)
    if (df$p_value[i] > 1e-4) return(TRUE)
    if (length(highCorrCols) == 0) return(FALSE)
    any(sapply(highCorrCols, function(col) df[i, ..col] == 1))
  })

  # Count total and remaining CS
  totalCs <- nrow(df)
  print("total_cs")
  print(totalCs)
  taggedCsCount <- sum(df$tagged_cs)
  if (totalCs > 0) {
    remainingCs <- totalCs - 1 - taggedCsCount
  } else {
    remainingCs <- 0
  }
  # Determine method
  df$method <- case_when(
  taggedCsCount == 0 & totalCs > 1 ~ "BVSR",
  (remainingCs == 0 & totalCs > 1) | (totalCs == 1) ~ "SER",
  remainingCs > 0 ~ "BCR",
  TRUE ~ NA_character_
)


  return(df)
}

