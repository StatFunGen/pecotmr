#' xQTL GWAS Enrichment Analysis
#'
#' This function processes GWAS and xQTL finemapped data files and then computes QTL enrichment.
#' For details on the parameters `pi_gwas`, `pi_qtl`, `lambda`, `ImpN`, and `num_threads`,
#' refer to the documentation of the `compute_qtl_enrichment` function.
#'
#' @param xqtl_files Vector of xQTL RDS file paths.
#' @param gwas_files Vector of GWAS RDS file paths.
#' @param xqtl_finemapping_obj Optional table name in xQTL RDS files (default 'susie_fit').
#' @param gwas_finemapping_obj Optional table name in GWAS RDS files (default 'susie_fit').
#' @param xqtl_varname_obj Optional table name in xQTL RDS files (default 'susie_fit').
#' @param gwas_varname_obj Optional table name in GWAS RDS files (default 'susie_fit').
#' @param pi_gwas Optional parameter for GWAS enrichment estimation (see `compute_qtl_enrichment`).
#' @param pi_qtl Optional parameter for xQTL enrichment estimation (see `compute_qtl_enrichment`).
#' @param lambda Shrinkage parameter for enrichment computation (see `compute_qtl_enrichment`).
#' @param ImpN Importance parameter for enrichment computation (see `compute_qtl_enrichment`).
#' @param num_threads Number of threads for parallel processing (see `compute_qtl_enrichment`).
#' @return The output from the compute_qtl_enrichment function.
#' @examples
#' gwas_files <- c("gwas_file1.rds", "gwas_file2.rds")
#' xqtl_files <- c("xqtl_file1.rds", "xqtl_file2.rds")
#' result <- xqtl_enrichment_wrapper(gwas_files, xqtl_files)
#' @export
xqtl_enrichment_wrapper <- function(xqtl_files, gwas_files,
                                    xqtl_finemapping_obj = NULL, gwas_finemapping_obj = NULL,
                                    xqtl_varname_obj = NULL, gwas_varname_obj = NULL,
                                    num_gwas = NULL, pi_qtl = NULL,
                                    lambda = 1.0, ImpN = 25,
                                    double_shrinkage = FALSE,
                                    bessel_correction = TRUE,
                                    num_threads = 1) {
  process_finemapped_data <- function(xqtl_files, gwas_files,
                                      xqtl_finemapping_obj = NULL, gwas_finemapping_obj = NULL,
                                      xqtl_varname_obj = NULL, gwas_varname_obj = NULL) {
    # Load and process GWAS data
    gwas_pip_list <- map(gwas_files, function(file) {
      raw_data <- readRDS(file)[[1]]
      gwas_data <- if (!is.null(gwas_finemapping_obj)) get_nested_element(raw_data, gwas_finemapping_obj) else raw_data
      pip <- gwas_data$pip
      if (!is.null(gwas_varname_obj)) names(pip) <- get_nested_element(raw_data, gwas_varname_obj)
      pip
    })

    # Check for unique variant names in GWAS pip vectors
    all_variant_names <- unique(unlist(map(gwas_pip_list, names)))
    if (length(unique(all_variant_names)) != length(all_variant_names)) {
      stop("Non-unique variant names found in GWAS data with different pip values.")
    }
    gwas_pip <- unlist(gwas_pip_list)

    # Process xQTL data
    xqtl_data <- lapply(xqtl_files, function(file) {
      raw_data <- readRDS(file)[[1]]
      xqtl_data <- tryCatch(
        {
          if (!is.null(xqtl_finemapping_obj)) get_nested_element(raw_data, xqtl_finemapping_obj) else raw_data
        },
        error = function(e) {
          return(NULL)
        }
      )
      if (!is.null(xqtl_data)) {
        list(
          alpha = xqtl_data$alpha,
          pip = setNames(xqtl_data$pip, get_nested_element(raw_data, xqtl_varname_obj)),
          prior_variance = xqtl_data$V
        )
      } else {
        NULL
      }
    })

    # Return results as a list
    return(list(gwas_pip = gwas_pip, xqtl_data = xqtl_data))
  }

  # Load data
  dat <- process_finemapped_data(xqtl_files, gwas_files, xqtl_finemapping_obj, gwas_finemapping_obj, xqtl_varname_obj, gwas_varname_obj)
  # Compute QTL enrichment
  return(compute_qtl_enrichment(
    gwas_pip = dat$gwas_pip, susie_qtl_regions = dat$xqtl_data,
    num_gwas = num_gwas, pi_qtl = pi_qtl,
    lambda = lambda, ImpN = ImpN,
    double_shrinkage = double_shrinkage,
    bessel_correction = bessel_correction,
    num_threads = num_threads
  ))
}

#' Function to filter and order colocalization results
#' @noRd
filter_and_order_coloc_results <- function(coloc_results_fil) {
  # Ensure the input has more than one column
  if (ncol(coloc_results_fil) <= 1) {
    stop("Insufficient number of columns in colocalization results")
  }

  cs_num <- ncol(coloc_results_fil) - 1
  map(seq_len(cs_num), function(n) {
    coloc_results_fil[, c(1, n + 1)] %>% .[order(.[, 2], decreasing = TRUE), ]
  })
}

#' Function to calculate cumulative sum
#' @noRd
calculate_cumsum <- function(coloc_results) {
  cumsum(coloc_results[, 2])
}

#' Load LD matrix for a set of variants, narrowing the region and aligning names.
#' @importFrom stringr str_split
#' @noRd
extract_ld_for_variants <- function(ld_meta_file_path, analysis_region, variants) {
  var_pos <- as.numeric(str_split(variants, ":", simplify = TRUE)[, 2])
  chr <- str_split(analysis_region, ":", simplify = TRUE)[, 1]
  region_narrow <- paste0(chr, ":", min(var_pos), "-", max(var_pos))
  ld_data <- load_LD_matrix(ld_meta_file_path, region = region_narrow,
                            return_genotype = "auto")
  if (!is(ld_data, "LDData")) {
    stop("load_LD_matrix must return an LDData object")
  }
  ld_variants <- getVariantIds(ld_data)
  has_geno <- hasGenotypes(ld_data)
  aligned <- align_variant_names(ld_variants, variants)
  # When genotypes available, compute R only for the needed variant subset
  if (has_geno) {
    X <- getGenotypes(ld_data)
    colnames(X) <- aligned$aligned_variants
    X_sub <- X[, variants, drop = FALSE]
    ld_matrix <- compute_LD(X_sub, method = "sample")
  } else {
    ld_matrix <- getCorrelation(ld_data)
    colnames(ld_matrix) <- rownames(ld_matrix) <- aligned$aligned_variants
    ld_matrix <- ld_matrix[variants, variants]
  }
  ld_matrix
}

#' Function to calculate purity
#' @noRd
calculate_purity <- function(variants, ext_ld, squared = FALSE) {
  # This is a placeholder for calculating purity, adjust as per your actual function
  purity <- matrix(susieR:::get_purity(variants, Xcorr = ext_ld, squared), 1, 3)
  purity
}

#' Main processing function
#' This function is designed to summarize coloc results based on the following criteria:
#' 1. Among the colocalized variant pairs, PPH4 has the highest value compared to PPH0-PPH3.
#' 2. PPH4 exceeds threshold, default as 0 since we advocate not using PPH4 concept but rather use CoS
#' 3. We aggregate variants and cumulatively sum their PPH4 values to form a credible set until the threshold, default as 0.95.
#' 4. The cs's purity is computed with the `get_purity` function from the `gaow/susieR` package, and the same purity criteria are employed to filter the credibility set.
#' @noRd
process_coloc_results <- function(coloc_result, LD_meta_file_path, analysis_region, PPH4_thres = 0, coverage = 0.95, min_abs_corr = 0.8, null_index = 0, coloc_index = "PP.H4.abf") {
  # Extract PIP values from coloc_result summary
  coloc_summary <- as.data.frame(coloc_result$summary)
  coloc_pip <- coloc_summary[, grepl("PP", colnames(coloc_summary))]

  # Filter and extract relevant columns from coloc_result results
  # PP.H4 is highest and > 0.8
  coloc_results_df <- as.data.frame(coloc_result$results)
  coloc_filter <- apply(coloc_pip, 1, function(row) {
    max_index <- which.max(row)
    max_value <- row[max_index]
    return(max_value > PPH4_thres && colnames(coloc_pip)[max_index] == coloc_index)
  })

  coloc_res <- list()

  if (sum(coloc_filter) > 0) {
    coloc_results_fil <- coloc_results_df[, c(1, which(coloc_filter) + 1), drop = FALSE]
    coloc_summary_fil <- coloc_summary[which(coloc_filter), , drop = FALSE]

    # prepare to calculate purity
    ordered_results <- filter_and_order_coloc_results(coloc_results_fil)
    cs <- map(ordered_results, function(res) {
      csm <- calculate_cumsum(res)
      res[, 1][1:min(which(csm > coverage))]
    })

    purity <- map_dfr(seq_along(cs), function(n) {
      variants <- normalize_variant_id(cs[[n]])
      if (null_index > 0 && null_index %in% variants) {
        data.frame(min.abs.corr = -9, mean.abs.corr = -9, median.abs.corr = -9)
      } else {
        ext_ld <- extract_ld_for_variants(LD_meta_file_path, analysis_region, variants)
        p <- calculate_purity(variants, ext_ld)
        data.frame(min.abs.corr = p[1, 1], mean.abs.corr = p[1, 2], median.abs.corr = p[1, 3])
      }
    })
    is_pure <- which(purity[, 1] >= min_abs_corr)

    # Finalize the result
    if (length(is_pure) > 0) {
      cs <- cs[is_pure]
      purity <- purity[is_pure, ]
      true_summary <- coloc_summary_fil[is_pure, ]
      coloc_res$sets <- list(cs = cs, purity = purity, true_summary = true_summary)
    }
  } else {
    message("Coloc results did not find any variants that satisfy the condition of PP.H4 being the highest value and > ", PPH4_thres)
    coloc_res$sets <- list(cs = NULL)
  }

  return(coloc_res)
}

# Extract and filter an LBF matrix from a finemapped data object.
# @noRd
.extract_lbf_matrix <- function(raw_data, finemapping_obj, varname_obj,
                                filter_lbf_cs, filter_lbf_cs_secondary, prior_tol) {
  fm_data <- if (!is.null(finemapping_obj)) {
    tryCatch(get_nested_element(raw_data, finemapping_obj),
      error = function(e) {
        message(paste("no", finemapping_obj[2], "in", finemapping_obj[1]))
        NULL
      }
    )
  } else {
    raw_data
  }
  if (is.null(fm_data)) return(NULL)

  lbf_matrix <- as.data.frame(fm_data$lbf_variable)
  # fSuSiE has a different structure
  if (is.null(lbf_matrix) || nrow(lbf_matrix) == 0) {
    lbf_matrix <- do.call(rbind, raw_data[[1]]$fsusie_result$lBF) %>% as.data.frame()
    if (nrow(lbf_matrix) > 0) message("This is a fSuSiE case")
  }

  # Filter rows
  if (filter_lbf_cs && is.null(filter_lbf_cs_secondary)) {
    lbf_matrix <- lbf_matrix[fm_data$sets$cs_index, , drop = FALSE]
  } else if (!is.null(filter_lbf_cs_secondary)) {
    lbf_matrix <- lbf_matrix[get_filter_lbf_index(fm_data, coverage = filter_lbf_cs_secondary), , drop = FALSE]
  } else {
    if ("V" %in% names(fm_data)) {
      lbf_matrix <- lbf_matrix[fm_data$V > prior_tol, , drop = FALSE]
    } else {
      message("No V found in original data.")
    }
  }

  # Set variant names and remove NA columns
  if (!is.null(varname_obj)) colnames(lbf_matrix) <- get_nested_element(raw_data, varname_obj)
  lbf_matrix <- lbf_matrix[, !is.na(colnames(lbf_matrix))]

  list(lbf_matrix = lbf_matrix, fm_data = fm_data)
}

# Extract LBF matrix from an rss_analysis_pipeline result object.
# Unlike .extract_lbf_matrix which navigates RDS-loaded nested lists,
# this works directly with the in-memory pipeline output structure.
# @noRd
.extract_lbf_from_pipeline_result <- function(pipeline_result,
                                               filter_lbf_cs, filter_lbf_cs_secondary,
                                               prior_tol) {
  method_names <- setdiff(names(pipeline_result), "rss_data_analyzed")
  if (length(method_names) == 0) return(NULL)

  method_result <- pipeline_result[[method_names[1]]]
  fm_result <- method_result$finemapping_result
  if (is.null(fm_result) || !is(fm_result, "FineMappingResult")) return(NULL)
  fm_data <- getTrimmedFit(fm_result)
  variant_names <- getVariantNames(fm_result)
  if (is.null(fm_data) || is.null(fm_data$lbf_variable)) return(NULL)

  lbf_matrix <- as.data.frame(fm_data$lbf_variable)

  # Row filtering — same logic as .extract_lbf_matrix
  if (filter_lbf_cs && is.null(filter_lbf_cs_secondary)) {
    lbf_matrix <- lbf_matrix[fm_data$sets$cs_index, , drop = FALSE]
  } else if (!is.null(filter_lbf_cs_secondary)) {
    lbf_matrix <- lbf_matrix[get_filter_lbf_index(fm_data, coverage = filter_lbf_cs_secondary), , drop = FALSE]
  } else if ("V" %in% names(fm_data)) {
    lbf_matrix <- lbf_matrix[fm_data$V > prior_tol, , drop = FALSE]
  }

  if (!is.null(variant_names) && length(variant_names) == ncol(lbf_matrix)) {
    colnames(lbf_matrix) <- variant_names
  }
  lbf_matrix <- lbf_matrix[, !is.na(colnames(lbf_matrix))]
  list(lbf_matrix = lbf_matrix, fm_data = fm_data)
}

# Save inline fine-mapping result to disk in a format compatible with the
# file-based reading path (readRDS(file)[[1]] + gwas_finemapping_obj/gwas_varname_obj).
# @noRd
.save_finemapping_result <- function(pipeline_result, save_path) {
  if (is.null(save_path) || is.null(pipeline_result)) return(invisible(NULL))
  method_names <- setdiff(names(pipeline_result), "rss_data_analyzed")
  if (length(method_names) == 0) return(invisible(NULL))
  method_result <- pipeline_result[[method_names[1]]]
  fm_result <- method_result$finemapping_result
  if (is.null(fm_result) || !is(fm_result, "FineMappingResult")) return(invisible(NULL))
  save_data <- list(
    susie_fit = getTrimmedFit(fm_result),
    variant_names = getVariantNames(fm_result)
  )
  saveRDS(list(save_data), save_path)
  message("Fine-mapping result saved to: ", save_path,
          "\n  Reuse with: gwas_files = '", save_path,
          "', gwas_finemapping_obj = 'susie_fit', gwas_varname_obj = 'variant_names'")
  invisible(save_path)
}

#' Colocalization Analysis Wrapper
#'
#' Processes xQTL and GWAS finemapped data for colocalization analysis.
#' GWAS data can come from pre-computed RDS files or from inline fine-mapping
#' via \code{\link{rss_analysis_pipeline}}.
#'
#' @param xqtl_file Path to the xQTL RDS file.
#' @param gwas_files Vector of paths to GWAS RDS files. Required when
#'   \code{run_finemapping = FALSE}. Ignored when \code{run_finemapping = TRUE}.
#' @param xqtl_finemapping_obj Optional path in xQTL RDS to the finemapping object.
#' @param gwas_finemapping_obj Optional path in GWAS RDS to the finemapping object.
#' @param xqtl_varname_obj Optional path in xQTL RDS to variant names.
#' @param gwas_varname_obj Optional path in GWAS RDS to variant names.
#' @param xqtl_region_obj Optional path in xQTL RDS to region info.
#' @param gwas_region_obj Optional path in GWAS RDS to region info.
#' @param filter_lbf_cs Logical. Filter LBF rows by credible set index.
#' @param filter_lbf_cs_secondary Coverage for secondary LBF filtering.
#' @param prior_tol Minimum prior variance to retain an effect (default 1e-9).
#' @param p1 Prior probability a SNP is associated with trait 1 (default 1e-4).
#' @param p2 Prior probability a SNP is associated with trait 2 (default 1e-4).
#' @param p12 Prior probability a SNP is associated with both traits (default 5e-6).
#' @param run_finemapping Logical. If TRUE, run GWAS fine-mapping inline via
#'   \code{\link{rss_analysis_pipeline}}. Default FALSE.
#' @param sumstat_path Path to GWAS summary statistics file. Required when
#'   \code{run_finemapping = TRUE}.
#' @param column_file_path Path to column mapping file for summary statistics.
#' @param LD_data LD reference data (LDData object or list). Required when
#'   \code{run_finemapping = TRUE}.
#' @param n_sample Sample size for GWAS.
#' @param n_case Number of cases for binary traits.
#' @param n_control Number of controls for binary traits.
#' @param region Genomic region string (e.g., "chr1:1000-2000").
#' @param qc_method QC method: "slalom", "dentist", or "none". Default "slalom".
#' @param finemapping_method Fine-mapping method. Default "susie_rss".
#' @param finemapping_opts List of fine-mapping options passed to
#'   \code{\link{rss_analysis_pipeline}}.
#' @param impute Logical. Run RAISS imputation. Default TRUE.
#' @param impute_opts List of imputation options.
#' @param save_finemapping_path Path to save fine-mapping result as RDS. The
#'   saved file can be reused via \code{gwas_files} with
#'   \code{gwas_finemapping_obj = "susie_fit"} and
#'   \code{gwas_varname_obj = "variant_names"}.
#' @param return_finemapping Logical. If TRUE and \code{run_finemapping = TRUE},
#'   include full fine-mapping result under \code{$gwas_finemapping}.
#' @param ... Additional arguments (currently unused).
#' @return A list containing the coloc results and the summarized sets.
#' @seealso \code{\link{rss_analysis_pipeline}}, \code{\link{coloc_post_processor}}
#' @importFrom dplyr bind_rows mutate across
#' @importFrom tidyr replace_na
#' @importFrom coloc coloc.bf_bf
#' @importFrom purrr map map_dfr
#' @export
coloc_wrapper <- function(xqtl_file, gwas_files = NULL,
                          xqtl_finemapping_obj = NULL, xqtl_varname_obj = NULL, xqtl_region_obj = NULL,
                          gwas_finemapping_obj = NULL, gwas_varname_obj = NULL, gwas_region_obj = NULL,
                          filter_lbf_cs = FALSE, filter_lbf_cs_secondary = NULL,
                          prior_tol = 1e-9, p1 = 1e-4, p2 = 1e-4, p12 = 5e-6,
                          run_finemapping = FALSE,
                          sumstat_path = NULL, column_file_path = NULL,
                          LD_data = NULL,
                          n_sample = 0, n_case = 0, n_control = 0,
                          region = NULL,
                          qc_method = "slalom",
                          finemapping_method = "susie_rss",
                          finemapping_opts = list(
                            L = 20, L_greedy = 5,
                            coverage = c(0.95, 0.7, 0.5),
                            signal_cutoff = 0.025,
                            min_abs_corr = 0.8
                          ),
                          impute = TRUE,
                          impute_opts = list(rcond = 0.01, R2_threshold = 0.6,
                                             minimum_ld = 5, lamb = 0.01),
                          save_finemapping_path = NULL,
                          return_finemapping = FALSE,
                          ...) {
  # --- Input validation ---
  if (!run_finemapping && is.null(gwas_files)) {
    stop("Either set run_finemapping = TRUE with GWAS sumstat inputs, or provide gwas_files paths to pre-computed results.")
  }
  if (run_finemapping && !is.null(gwas_files)) {
    warning("Both run_finemapping = TRUE and gwas_files provided. Inline fine-mapping will be used; gwas_files ignored.")
    gwas_files <- NULL
  }
  if (run_finemapping) {
    if (is.null(sumstat_path)) stop("sumstat_path is required when run_finemapping = TRUE.")
    if (is.null(LD_data)) stop("LD_data is required when run_finemapping = TRUE.")
  }

  gwas_pipeline_result <- NULL

  if (run_finemapping) {
    # --- Inline fine-mapping path: QC runs inside rss_analysis_pipeline ---
    gwas_pipeline_result <- rss_analysis_pipeline(
      sumstat_path = sumstat_path, column_file_path = column_file_path,
      LD_data = LD_data,
      n_sample = n_sample, n_case = n_case, n_control = n_control,
      region = region,
      qc_method = qc_method, finemapping_method = finemapping_method,
      finemapping_opts = finemapping_opts,
      impute = impute, impute_opts = impute_opts
    )

    # Save to disk before extraction (useful even if extraction fails)
    .save_finemapping_result(gwas_pipeline_result, save_finemapping_path)

    gwas_extracted <- .extract_lbf_from_pipeline_result(
      gwas_pipeline_result, filter_lbf_cs, filter_lbf_cs_secondary, prior_tol
    )
    if (is.null(gwas_extracted)) {
      coloc_res <- list("No GWAS fine-mapping results produced by inline pipeline.")
      result <- c(coloc_res, analysis_region = region)
      if (return_finemapping) result$gwas_finemapping <- gwas_pipeline_result
      return(result)
    }
    combined_gwas_lbf_matrix <- gwas_extracted$lbf_matrix %>%
      as.data.frame() %>% mutate(across(everything(), ~ replace_na(., 0)))
  } else {
    # --- File-based path (unchanged) ---
    gwas_lbf_matrices <- map(gwas_files, function(file) {
      raw_data <- readRDS(file)[[1]]
      .extract_lbf_matrix(raw_data, gwas_finemapping_obj, gwas_varname_obj,
                          filter_lbf_cs, filter_lbf_cs_secondary, prior_tol)$lbf_matrix
    })
    combined_gwas_lbf_matrix <- bind_rows(gwas_lbf_matrices) %>%
      mutate(across(everything(), ~ replace_na(., 0)))
  }

  # Process xQTL data
  xqtl_raw_data <- readRDS(xqtl_file)[[1]]
  xqtl_extracted <- .extract_lbf_matrix(xqtl_raw_data, xqtl_finemapping_obj, xqtl_varname_obj,
                                        filter_lbf_cs, filter_lbf_cs_secondary, prior_tol)

  if (!is.null(xqtl_extracted)) {
    xqtl_lbf_matrix <- xqtl_extracted$lbf_matrix
    if (nrow(combined_gwas_lbf_matrix) > 0 && nrow(xqtl_lbf_matrix) > 0) {
      colnames(xqtl_lbf_matrix) <- align_variant_names(colnames(xqtl_lbf_matrix), colnames(combined_gwas_lbf_matrix))$aligned_variants
      common_colnames <- intersect(colnames(xqtl_lbf_matrix), colnames(combined_gwas_lbf_matrix))

      num_dropped_cols <- ncol(xqtl_lbf_matrix) - length(common_colnames)
      if (num_dropped_cols > 0) {
        message("Number of columns dropped from xQTL matrix: ", num_dropped_cols)
      }

      xqtl_lbf_matrix <- xqtl_lbf_matrix[, common_colnames, drop = FALSE] %>% as.matrix()
      combined_gwas_lbf_matrix <- combined_gwas_lbf_matrix[, common_colnames, drop = FALSE] %>% as.matrix()

      convert_to_string <- function(df) paste0("chr", df$chrom, ":", df$start, "-", df$end)
      analysis_region_out <- if (!is.null(xqtl_region_obj)) {
        get_nested_element(xqtl_raw_data, xqtl_region_obj) %>% convert_to_string()
      } else {
        region
      }

      coloc_res <- coloc.bf_bf(xqtl_lbf_matrix, combined_gwas_lbf_matrix, p1 = p1, p2 = p2, p12 = p12)
    } else {
      coloc_res <- list("No coloc results due to the absence of a GWAS log Bayes factor matrix filtered by prior tolerance.")
      analysis_region_out <- region
    }
  } else {
    coloc_res <- list(paste("no", xqtl_finemapping_obj[2], "in", xqtl_finemapping_obj[1]))
    analysis_region_out <- region
  }

  result <- c(coloc_res, analysis_region = analysis_region_out)
  if (return_finemapping && !is.null(gwas_pipeline_result)) {
    result$gwas_finemapping <- gwas_pipeline_result
  }
  return(result)
}

#' coloc_post_processor function
#' @param coloc_res coloc results from coloc.susie.
#' @param LD_meta_file_path Path to the metadata of LD reference.
#' @param analysis_region Path to the analysis region of coloc result.
#' @return A list containing the coloc results and post processed coloc sets.
#' @export
coloc_post_processor <- function(coloc_res, LD_meta_file_path = NULL, analysis_region = NULL, ...) {
  if (!is.null(LD_meta_file_path)) {
    if (is.null(analysis_region)) {
      stop("LD_meta_file_path is provided but analysis_region is not provided. Please provide analysis_region for purity filter.")
    }
    # Perform purity filter using LD_meta_file_path and analysis_region
    coloc_res <- c(coloc_res, process_coloc_results(coloc_res, LD_meta_file_path, analysis_region = analysis_region))
  } else {
    if (!is.null(analysis_region)) {
      warning("Analysis_region is provided but will not be used as LD_meta_file_path is not provided.")
    }
    warning("LD_meta_file_path not provided. Purity filter cannot be applied.")
  }
  return(coloc_res)
}

# In practice, analysis will contain two lines:
# res <- coloc_wrapper(...)
# post_processed_res <- coloc_post_processor
