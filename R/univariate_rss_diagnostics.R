#' Extract SuSiE Results from Finemapping Data
#'
#' This function extracts the trimmed SuSiE results from a finemapping data object,
#' typically obtained from a finemapping RDS file. It's designed to work with
#' the method layer of these files, often named as 'method_RAISS_imputed', 'method',
#' or 'method_NO_QC'. This layer is right under the study layer.
#'
#' @param con_data List. The method layer data from a finemapping RDS file.
#'
#' @return The trimmed SuSiE results (`$susie_result_trimmed`) if available,
#' otherwise NULL.
#'
#' @details
#' The function checks if the input data is empty or if the `$susie_result_trimmed`
#' element is missing. It returns NULL in these cases. If `$susie_result_trimmed`
#' exists and is not empty, it returns this element.
#'
#' @note
#' This function is particularly useful when working with large datasets
#' where not all method layers may contain valid SuSiE results or method layer.
#'
#' @export
get_susie_result <- function(con_data) {
    if (length(con_data) == 0) return(NULL)
    if (length(con_data$susie_result_trimmed) == 0) {
        return(NULL)
        print(paste("$susie_result_trimmed is null for", con_data))
    } else {
        return(con_data$susie_result_trimmed)
    }
}

#' Process Credible Sets (CS) from Finemapping Results
#'
#' This function extracts and processes information for each Credible Set (CS) 
#' from finemapping results, typically obtained from a finemapping RDS file.
#'
#' @param con_data List. The method layer data from a finemapping RDS file that is not empty.
#' @param cs_names Character vector. Names of the Credible Sets, usually in the format "L_<number>".
#' @param top_loci_table Data frame. The $top_loci layer data from the finemapping results.
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
extract_cs_info <- function(con_data, cs_names, top_loci_table) {  
  results <- map(seq_along(cs_names), function(i) {
    cs_name <- cs_names[i]
    indices <- con_data$susie_result_trimmed$sets$cs[[cs_name]]
    
    # Get variants for this CS using the full variant_names list
    cs_variants <- con_data$variant_names[indices]
    cs_data <- top_loci_table[top_loci_table$variant_id %in% cs_variants, ]
    top_row <- which.max(cs_data$pip)
    
    top_variant <- cs_data$variant_id[top_row]
    # Find the global index of the top variant
    top_variant_global_index = which(con_data$variant_names == top_variant)
    top_pip <- cs_data$pip[top_row]
    top_z <- cs_data$z[top_row]
    p_value <- z_to_pvalue(top_z)
    
    # Extract cs_corr
    cs_corr <- if (length(cs_names) > 1) {
      con_data$susie_result_trimmed$cs_corr[i,]
    } else {
      NA  # Use NA for the second CS or when there's only one CS
    }
    
    # Return results for this CS as a one-row data.frame
    result = tibble(
      cs_name = cs_name,
      variants_per_cs = length(cs_variants),
      top_variant = top_variant,
      top_variant_index = top_variant_global_index,
      top_pip = top_pip,
      top_z = top_z,
      p_value = p_value,
      cs_corr = list(paste(cs_corr, collapse = ","))  # list column if cs_corr is a vector
    )
    return(result)
  })
  # Combine all tibbles into one data frame
  final_result <- bind_rows(results)
  return(final_result)
}

#' Extract Information for Top Variant from Finemapping Results
#'
#' This function extracts information about the variant with the highest Posterior 
#' Inclusion Probability (PIP) from finemapping results, typically used when no 
#' Credible Sets (CS) are identified in the analysis.
#'
#' @param con_data List. The method layer data from a finemapping RDS file.
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
#' \code{\link{extract_cs_info}} for processing when Credible Sets are present.
#'
#' @export
extract_top_pip_info <- function(con_data) {
  # Find the variant with the highest PIP
  top_pip_index <- which.max(con_data$susie_result_trimmed$pip)
  top_pip <- con_data$susie_result_trimmed$pip[top_pip_index]
  top_variant <- con_data$variant_names[top_pip_index]
  top_z <- con_data$sumstats$z[top_pip_index]
  p_value <- z_to_pvalue(top_z)
  
  list(
    cs_name = NA,
    variants_per_cs = NA,
    top_variant = top_variant,
    top_variant_index = top_pip_index,
    top_pip = top_pip,
    top_z = top_z,
    p_value = p_value,
    cs_corr = NA  # or NULL
  )
}

#' Parse Credible Set Correlations from extract_cs_info() Output
#'
#' This function takes the output from `extract_cs_info()` and expands the `cs_corr` column
#' into multiple columns, preserving the original order of correlations. It also
#' calculates maximum and minimum correlation values for each Credible Set.
#'
#' @param df Data frame. The output from `extract_cs_info()` function,
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
parse_cs_corr <- function(df) {
  # Ensure we work with a data frame
  df <- as.data.frame(df)

  extract_correlations <- function(x) {
    # Early return if x is invalid
    if(is.na(x) || x == "" || is.null(x) || !grepl(",", as.character(x))) {
      return(list(values = numeric(0), max_corr = NA_real_, min_corr = NA_real_))
    }

    # Convert and filter values
    values <- as.numeric(unlist(strsplit(x, ",")))
    values_filtered <- abs(values[values != 1])

    # Return list with NA if no valid correlations
    list(
      values = values,
      max_corr = if(length(values_filtered) > 0) max(abs(values_filtered), na.rm = TRUE) else NA_real_,
      min_corr = if(length(values_filtered) > 0) min(abs(values_filtered), na.rm = TRUE) else NA_real_
    )
  }
  # Process correlations
  processed_results <- lapply(df$cs_corr, extract_correlations)
  # If no valid results, add NA columns and return
  if(all(sapply(processed_results, function(x) length(x$values) == 0))) {
    df$cs_corr_max <- NA_real_
    df$cs_corr_min <- NA_real_
    return(df)
  }

  # Determine max number of correlations
  max_corr_count <- max(sapply(processed_results, function(x) length(x$values)))

  # Create and add correlation columns
  col_names <- paste0("cs_corr_", 1:max_corr_count)

  for(i in seq_along(col_names)) {
    df[[col_names[i]]] <- sapply(processed_results, function(x) {
      if(length(x$values) >= i) x$values[i] else NA_real_
    })
  }

  # Add max and min correlation columns
  df$cs_corr_max <- sapply(processed_results, `[[`, "max_corr")
  df$cs_corr_min <- sapply(processed_results, `[[`, "min_corr")

  return(df)
}

#' Process Credible Set Information and Determine Updating Strategy
#'
#' This function categorizes Credible Sets (CS) within a study block into different 
#' updating strategies based on their statistical properties and correlations.
#'
#' @param df Data frame. Contains information about Credible Sets for a specific study and block.
#' @param high_corr_cols Character vector. Names of columns in df that represent high correlations.
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
auto_decision <- function(df, high_corr_cols) {
  # Identify top_cs
  top_cs_index <- which.max(abs(df$top_z))
  df$top_cs <- FALSE
  df$top_cs[top_cs_index] <- TRUE

  # Identify tagged_cs
  df$tagged_cs <- sapply(1:nrow(df), function(i) {
    if (df$top_cs[i]) return(FALSE)
    if (df$p_value[i] > 1e-4) return(TRUE)
    if (length(high_corr_cols) == 0) return(FALSE)
    any(sapply(high_corr_cols, function(col) df[i, ..col] == 1))
  })

  # Count total and remaining CS
  total_cs <- nrow(df)
  print("total_cs")
  print(total_cs)
  tagged_cs_count <- sum(df$tagged_cs)
  if (total_cs > 0) {
    remaining_cs <- total_cs - 1 - tagged_cs_count
  } else {
    remaining_cs <- 0
  }
  # Determine method
  df$method <- case_when(
  tagged_cs_count == 0 & total_cs > 1 ~ "BVSR",
  (remaining_cs == 0 & total_cs > 1) | (total_cs == 1) ~ "SER",
  remaining_cs > 0 ~ "BCR",
  TRUE ~ NA_character_
)


  return(df)
}
