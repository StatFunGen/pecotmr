#' Function to perform allele flip QC and harmonization on the weights and GWAS against LD for a region.
#' FIXME: GWAS loading function from Haochen for both tabix & column-mapping yml application
#'
#' Function Conditions:
#' - processes data in the format of either the output from load_twas_weights/generate_twas_db or
#'   refined_twas_weights_data from twas pipeline.
#' - For the first format, we expect there is only one gene/events's information, that can be accessed through `region_info_obj`
#'   and refined_twas_weights_data contains per region multiple gene/event's refined weights data.
#'
#' Main Steps:
#' 1. allele QC for TWAS weights against the LD meta
#' 2. allele QC for GWA summary stats against the LD meta
#' 3. adjust susie/mvsusie weights based on the overlap variants
#'
#' @param twas_weights_data List of list of twas weights output from from generate_twas_db function.
#' @param gwas_meta_file A file path for a dataframe table with column of "study_id", "chrom" (integer), "file_path",
#' "column_mapping_file". Each file in "file_path" column is tab-delimited dataframe of GWAS summary statistics with column name
#' "chrom" (or #chrom" if tabix-indexed), "pos", "A2", "A1".
#' @param ld_meta_file_path Path to LD reference: either a PLINK2/PLINK1 prefix, or a tab-delimited
#'   metadata file with columns "#chrom", "start", "end", "path" (auto-detected).
#' @param ld_reference_sample_size Sample size of the LD reference panel (integer). Required.
#'   Used to compute per-variant variance as 2*p*(1-p)*n/(n-1). For ADSP R4, use 17000.
#' @return A list of list for harmonized weights and dataframe of gwas summary statistics that is add to the original input of
#' twas_weights_data under each context.
#' @importFrom vroom vroom
#' @importFrom readr parse_number
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom IRanges IRanges findOverlaps start end reduce
#' @export
harmonize_twas <- function(twas_weights_data, ld_meta_file_path, gwas_meta_file,
                           ld_reference_sample_size, column_file_path = NULL, comment_string = "#") {
  # Function to group contexts based on start and end positions
  group_contexts_by_region <- function(twas_weights_data, molecular_id, chrom, tolerance = 5000) {
    region_info_df <- do.call(rbind, lapply(names(twas_weights_data$weights), function(context) {
      wgt_range <- parse_variant_id(rownames(twas_weights_data[["weights"]][[context]]))$pos
      data.frame(context = context, start = min(wgt_range), end = max(wgt_range))
    }))
    if (nrow(region_info_df) == 1) {
      # Handle case with only one context
      single_context_group <- list(
        context_group_1 = list(
          contexts = region_info_df$context,
          query_region = paste0(chrom, ":", region_info_df$start, "-", region_info_df$end),
          all_variants = unique(rownames(twas_weights_data[["weights"]][[region_info_df$context]]))
        )
      )
      return(single_context_group)
    }
    # Calculate distance matrix and perform hierarchical clustering
    clusters <- cutree(hclust(dist(region_info_df[, c("start", "end")])), h = tolerance)
    # Group contexts and determine query regions
    region_groups <- split(region_info_df, clusters) %>%
      lapply(function(group) {
        list(contexts = group$context, query_region = paste0(
          chrom, ":", min(group$start),
          "-", max(group$end)
        ))
      })
    # Create IRanges objects and merge overlapping intervals
    intervals <- IRanges(start = unlist(lapply(region_groups, function(context_group) {
      as.numeric(gsub(
        "^.*:\\s*|\\s*-.*$", "",
        context_group$query_region
      ))
    })), end = unlist(lapply(
      region_groups,
      function(context_group) as.numeric(sub("^.*?\\-", "", context_group$query_region))
    )))
    reduced_intervals <- reduce(intervals)

    # Find which original groups are merged, and update region_groups lists
    overlaps <- findOverlaps(intervals, reduced_intervals)
    # Create merged groups based on overlap mapping
    merged_groups <- lapply(seq_along(reduced_intervals), function(i) {
      context_indices <- queryHits(overlaps)[subjectHits(overlaps) == i]
      merged_contexts <- unlist(lapply(context_indices, function(idx) region_groups[[idx]]$contexts))
      list(contexts = merged_contexts, query_region = paste0(chrom, ":", start(reduced_intervals[i]), "-", end(reduced_intervals[i])))
    })
    names(merged_groups) <- paste0("context_group_", seq_along(merged_groups))
    # add variant names for coordinate extraction
    for (group in names(merged_groups)) {
      contexts <- merged_groups[[group]]$contexts
      merged_groups[[group]]$all_variants <- unique(do.call(c, lapply(
        contexts,
        function(context) {
          rownames(twas_weights_data[["weights"]][[context]])
        }
      )))
    }
    return(merged_groups)
  }

  # Step 1: load TWAS weights data
  molecular_ids <- names(twas_weights_data)
  chrom <- as.integer(parse_number(gsub(":.*$", "", rownames(twas_weights_data[[1]]$weights[[1]])[1])))
  gwas_meta_df <- as.data.frame(vroom(gwas_meta_file))
  gwas_files <- unique(gwas_meta_df$file_path[gwas_meta_df$chrom == chrom])
  names(gwas_files) <- unique(gwas_meta_df$study_id[gwas_meta_df$chrom == chrom])
  results <- list()

  # Step 2: Load LD for all events/genes by clustered context region
  for (molecular_id in molecular_ids) {
    twas_weights_data[[molecular_id]][["variant_names"]] <- lapply(twas_weights_data[[molecular_id]]$weights, function(x) rownames(x))
  }
  region_variants <- variant_id_to_df(unique(do.call(c, find_data(twas_weights_data, c(2, "variant_names")))))
  region_of_interest <- data.frame(chrom = chrom, start = min(region_variants$pos), end = max(region_variants$pos))
  LD_list <- load_LD_matrix(ld_meta_file_path, region_of_interest, region_variants,
                            n_sample = ld_reference_sample_size)
  # remove duplicate variants
  dup_idx <- which(duplicated(LD_list$LD_variants))
  if (length(dup_idx) >= 1) {
    LD_list$LD_variants <- LD_list$LD_variants[-dup_idx]
    LD_list$LD_matrix <- LD_list$LD_matrix[-dup_idx, -dup_idx]
    LD_list$ref_panel <- LD_list$ref_panel[-dup_idx, ]
  }

  # loop through genes/events:
  for (molecular_id in molecular_ids) {
    mol_data <- twas_weights_data[[molecular_id]]
    mol_res <- list(chrom = chrom, variant_names = list())
    mol_res[["data_type"]] <- if ("data_type" %in% names(mol_data)) mol_data$data_type

    # group contexts based on the variant position
    context_clusters <- group_contexts_by_region(mol_data, molecular_id, chrom, tolerance = 5000)

    # loop through contexts: grouping contexts can be useful during TWAS data harmonization to stratify variants for LD loading
    for (context_group in names(context_clusters)) {
      cluster <- context_clusters[[context_group]]
      contexts <- cluster$contexts
      query_region <- cluster$query_region
      region_of_interest <- region_to_df(query_region)
      all_variants <- variant_id_to_df(cluster$all_variants)

      # Step 3: load GWAS data for clustered context groups
      for (study in names(gwas_files)) {
        gwas_file <- gwas_files[study]
        gwas_data_sumstats <- harmonize_gwas(gwas_file, query_region=query_region,
                                              LD_list$LD_variants, c("beta", "z"),
                                              match_min_prop = 0, column_file_path = column_file_path, comment_string = comment_string)
        if(is.null(gwas_data_sumstats)) next
        # loop through context within the context group:
        for (context in contexts) {
          weights_matrix <- mol_data[["weights"]][[context]]

          # Step 4: harmonize weights, flip allele
          weights_matrix <- cbind(variant_id_to_df(rownames(weights_matrix)), weights_matrix)
          weights_matrix_qced <- match_ref_panel(weights_matrix, LD_list$LD_variants,
            colnames(weights_matrix)[!colnames(weights_matrix) %in% c("chrom", "pos", "A2", "A1")],
            match_min_prop = 0
          )
          qced_data <- weights_matrix_qced$target_data_qced
          weights_matrix_subset <- as.matrix(qced_data[, !colnames(qced_data) %in% c(
            "chrom", "pos", "A2", "A1", "variant_id", "variants_id_original"
          ), drop = FALSE])
          rownames(weights_matrix_subset) <- qced_data$variant_id

          # intersect post-qc gwas and post-qc weight variants (all now in canonical chr-prefix format)
          gwas_LD_variants <- intersect(gwas_data_sumstats$variant_id, LD_list$LD_variants)
          weights_matrix_subset <- weights_matrix_subset[rownames(weights_matrix_subset) %in% gwas_data_sumstats$variant_id, , drop = FALSE]
          if (nrow(weights_matrix_subset) == 0) next
          postqc_weight_variants <- rownames(weights_matrix_subset)

          # Step 5: adjust SuSiE weights based on available variants
          if ("susie_weights" %in% colnames(mol_data[["weights"]][[context]])) {
            adjusted_susie_weights <- adjust_susie_weights(mol_data,
              keep_variants = postqc_weight_variants, run_allele_qc = TRUE,
              variable_name_obj = c("variant_names", context),
              susie_obj = c("susie_results", context),
              twas_weights_table = c("weights", context), postqc_weight_variants, match_min_prop = 0
            )
            weights_matrix_subset <- cbind(
              susie_weights = setNames(adjusted_susie_weights$adjusted_susie_weights, adjusted_susie_weights$remained_variants_ids),
              weights_matrix_subset[adjusted_susie_weights$remained_variants_ids, !colnames(weights_matrix_subset) %in% "susie_weights", drop = FALSE]
            )
            susie_intermediate <- mol_data$susie_results[[context]][c("pip", "cs_variants", "cs_purity")]
            names(susie_intermediate[["pip"]]) <- rownames(weights_matrix) # original variants that is not qced yet
            pip <- susie_intermediate[["pip"]]
            pip_qced <- match_ref_panel(cbind(parse_variant_id(names(pip)), pip), LD_list$LD_variants, "pip", match_min_prop = 0)
            susie_intermediate[["pip"]] <- abs(pip_qced$target_data_qced$pip)
            names(susie_intermediate[["pip"]]) <- pip_qced$target_data_qced$variant_id
            susie_intermediate[["cs_variants"]] <- lapply(susie_intermediate[["cs_variants"]], function(x) {
              variant_qc <- match_ref_panel(x, LD_list$LD_variants, match_min_prop = 0)
              variant_qc$target_data_qced$variant_id[variant_qc$target_data_qced$variant_id %in% postqc_weight_variants]
            })
            mol_res[["susie_weights_intermediate_qced"]][[context]] <- susie_intermediate
          }
          rm(weights_matrix) # context specific original weight matrix
          gc()

          if (nrow(weights_matrix_subset) == 0) {
              warning("weights_matrix_subset is empty. Skipping this context.")
              next
          }
          mol_res[["variant_names"]][[context]][[study]] <- rownames(weights_matrix_subset)

          # Step 6: scale weights by variance (from ref_panel, populated by load_LD_matrix)
          variance <- LD_list$ref_panel$variance[match(rownames(weights_matrix_subset), LD_list$ref_panel$variant_id)]
          mol_res[["weights_qced"]][[context]][[study]] <- list(scaled_weights = weights_matrix_subset * sqrt(variance), weights = weights_matrix_subset)
        }
        # Combine gwas sumstat across different context for a single context group (all variant_ids now in canonical format)
        gwas_data_sumstats <- gwas_data_sumstats[gwas_data_sumstats$variant_id %in% unique(find_data(mol_res[["variant_names"]], c(2, study))), , drop = FALSE]
        mol_res[["gwas_qced"]][[study]] <- rbind(mol_res[["gwas_qced"]][[study]], gwas_data_sumstats)
        gwas_qced <- mol_res[["gwas_qced"]][[study]]
        mol_res[["gwas_qced"]][[study]] <- gwas_qced[!duplicated(gwas_qced[, c("variant_id", "z")]), ]
      }
    }
    twas_weights_data[[molecular_id]] <- NULL
    # extract LD matrix for variants intersect with gwas and twas weights at molecular_id level
    all_molecular_variants <- unique(find_data(mol_res[["gwas_qced"]], c(2, "variant_id")))
    if (is.null(all_molecular_variants)) {
      results[[molecular_id]] <- NULL
    } else {
      # All variant IDs are now in canonical chr-prefix format
      var_indx <- match(all_molecular_variants, LD_list$LD_variants)
      mol_res[["LD"]] <- as.matrix(LD_list$LD_matrix[var_indx, var_indx])
      results[[molecular_id]] <- mol_res
    }
  }
  # return results
  return(list(twas_data_qced = results, ref_panel = LD_list$ref_panel))
}

#' Harmonize GWAS Summary Statistics 
#' perform harmonization on gwas summary statistics for a chromosome data or specific queried region
#' @param gwas_file A string for the file path of gwas summary statistics file that is already tabix indexed 
#' @param query_region A string for region of query for tabix-indexed gwas summary statistics file in the format of chr:start-end
#' @noRd
#' @export
harmonize_gwas <- function(gwas_file, query_region, ld_variants, col_to_flip=NULL, match_min_prop=0, column_file_path=NULL, comment_string="#"){
    if(is.null(gwas_file)| is.na(gwas_file)) stop("No GWAS file path provided. ")
    if (!is.null(column_file_path)) {
      rss_result <- load_rss_data(
        sumstat_path = gwas_file,
        column_file_path = column_file_path,
        region = query_region,
        comment_string = comment_string
      )
      gwas_data_sumstats <- rss_result$sumstats
    } else {
      gwas_data_sumstats <- as.data.frame(tabix_region(gwas_file, query_region))
      if (nrow(gwas_data_sumstats) > 0) {
        gwas_data_sumstats <- standardise_sumstats_columns(gwas_data_sumstats)
      }
    }
    if (nrow(gwas_data_sumstats) == 0) {
        if (length(names(gwas_file))==0) names(gwas_file) <- gwas_file
        warning(paste0("No GWAS summary statistics found for the region of ", query_region, " in ", names(gwas_file), ". "))
        return(NULL)
    }
    # Check if sumstats has z-scores or (beta and se)
    if (!is.null(gwas_data_sumstats$z)) {
      # z-scores already present, nothing to do
    } else if (!is.null(gwas_data_sumstats$beta) && !is.null(gwas_data_sumstats$se)) {
      gwas_data_sumstats$z <- gwas_data_sumstats$beta / gwas_data_sumstats$se
    } else {
      stop("gwas_data_sumstats should have 'z' or ('beta' and 'se') columns")
    }
    # check for overlapping variants
    if (!any(gwas_data_sumstats$pos %in% gsub("\\:.*$", "", sub("^.*?\\:", "", ld_variants)))) return(NULL)
    gwas_allele_flip <- match_ref_panel(gwas_data_sumstats, ld_variants, col_to_flip=col_to_flip, match_min_prop = match_min_prop)
    gwas_data_sumstats <- gwas_allele_flip$target_data_qced # post-qc gwas data that is flipped and corrected - gwas study level
    gwas_data_sumstats <- gwas_data_sumstats[!is.na(gwas_data_sumstats$z) & !is.infinite(gwas_data_sumstats$z), ]
    return(gwas_data_sumstats)
}

#' Function to perform TWAS analysis for across multiple contexts.
#' This function peforms TWAS analysis for multiple contexts for imputable genes within an LD region and summarize the twas results.
#' @param twas_weights_data List of list of twas weights output from generate_twas_db function.
#' @param region_block A string with LD region informaiton of chromosome number, star and end position of LD block conneced with "_".
#' @return A list of list containing twas result table and formatted TWAS data compatible with ctwas_sumstats() function.
#' \itemize{
#'   \item{twas_table}{ A dataframe of twas results summary is generated for each gene-contexts-method pair of all methods for imputable genes.}
#'   \item{twas_data}{ A list of list containing formatted TWAS data.}
#' }
# Shared shape for twas_analysis() result rows. Internal.
build_twas_score_row <- function(twas_rs, weight_db, context, study) {
  if (is.null(twas_rs)) return(data.frame())
  data.frame(
    gwas_study   = study,
    method       = sub("_[^_]+$", "", names(twas_rs)),
    twas_z       = find_data(twas_rs, c(2, "z")),
    twas_pval    = find_data(twas_rs, c(2, "pval")),
    context      = context,
    molecular_id = weight_db
  )
}

#' @importFrom stringr str_remove
#' @export
twas_pipeline <- function(twas_weights_data,
                          ld_meta_file_path,
                          gwas_meta_file,
                          region_block,
                          ld_reference_sample_size,
                          rsq_cutoff = 0.01,
                          rsq_pval_cutoff = 0.05,
                          rsq_option = c("rsq", "adj_rsq"),
                          rsq_pval_option = c("pval", "adj_rsq_pval"),
                          mr_pval_cutoff = 0.05,
                          mr_coverage_column = "cs_coverage_0.95",
                          output_twas_data = FALSE,
                          event_filters=NULL,
                          column_file_path = NULL,
                          comment_string="#") {
  # internal function to format TWAS output
  format_twas_data <- function(post_qc_twas_data, twas_table) {
    weights_list <- purrr::map(names(post_qc_twas_data), function(molecular_id) {
      mol <- post_qc_twas_data[[molecular_id]]
      contexts <- names(mol[["weights_qced"]])
      mol_chrom <- mol[["chrom"]]
      model_sel <- mol[["model_selection"]]

      purrr::map(contexts, function(context) {
        data_type <- mol[["data_type"]][[context]]
        if (!is.null(model_sel) && is.list(model_sel) && length(model_sel) > 0) {
          is_imputable <- model_sel[[context]]$is_imputable
          model_selected <- if (isTRUE(is_imputable)) model_sel[[context]]$selected_model else NA
        } else {
          model_selected <- NA
          is_imputable <- NA
        }
        if (is.null(model_selected) || !isTRUE(is_imputable)) return(NULL)

        gwas_studies <- names(mol[["weights_qced"]][[context]])
        weight_key <- paste0(molecular_id, "|", data_type, "_", context)
        study_entries <- purrr::map(gwas_studies, function(study) {
          ctx_weights <- mol[["weights_qced"]][[context]][[study]]
          scaled_wgt <- ctx_weights[["scaled_weights"]][, paste0(model_selected, "_weights"), drop = FALSE]
          colnames(scaled_wgt) <- "weight"
          context_variants <- rownames(ctx_weights[["scaled_weights"]])
          context_range <- parse_variant_id(context_variants)$pos
          entry <- list(list(
            chrom = mol_chrom, p0 = min(context_range), p1 = max(context_range),
            wgt = scaled_wgt, molecular_id = molecular_id,
            weight_name = paste0(data_type, "_", context), type = data_type,
            context = context, n_wgt = length(context_variants)
          ))
          names(entry) <- study
          result <- list(entry)
          names(result) <- weight_key
          result
        }) %>% purrr::list_flatten()
        study_entries
      }) %>% purrr::compact() %>% purrr::list_flatten()
    }) %>% purrr::list_flatten()
    weights <- purrr::compact(weights_list)
    # Optional susie_weights_intermediate_qced processing
    if ("susie_weights_intermediate_qced" %in% names(post_qc_twas_data[[1]])) {
      susie_weights_intermediate_qced <- setNames(lapply(
        names(post_qc_twas_data),
        function(x) post_qc_twas_data[[x]]$susie_weights_intermediate_qced
      ), names(post_qc_twas_data))
    } else {
      susie_weights_intermediate_qced <- NULL
    }

    # gene_z table
    if ("is_selected_method" %in% colnames(twas_table)) {
      twas_table <- twas_table[na.omit(twas_table$is_selected_method), , drop = FALSE]
    }
    if (nrow(twas_table) > 0) {
      twas_table$id <- paste0(twas_table$molecular_id, "|", twas_table$type, "_", twas_table$context)
      twas_table$group <- paste0(twas_table$context, "|", twas_table$type)
      
      twas_table$z <- twas_table$twas_z
      
      output_columns <- c("id", "z", "type", "context", "group", "gwas_study")
      twas_table <- twas_table[, intersect(output_columns, colnames(twas_table)), drop = FALSE]
      studies <- unique(twas_table$gwas_study)
      z_gene_list <- list()
      z_snp <- list()
      for (study in studies) {
        z_gene_list[[study]] <- twas_table[twas_table$gwas_study == study, , drop = FALSE]
      }
      result <- list(weights = weights, z_gene = z_gene_list)
      if (!is.null(susie_weights_intermediate_qced)) {
        result$susie_weights_intermediate_qced <- susie_weights_intermediate_qced
      }
      return(result)
    } else {
      return(NULL)
    }
  }
  pick_best_model <- function(twas_data_combined, rsq_cutoff, rsq_pval_cutoff, rsq_option, rsq_pval_option) {
    best_rsq <- rsq_cutoff
    # Determine if a gene/region is imputable and select the best model
    model_selection <- lapply(names(twas_data_combined$weights), function(context) {
      selected_model <- NULL
      available_models <- do.call(c, lapply(names(twas_data_combined$twas_cv_performance[[context]]), function(model) {
        if (!is.na(twas_data_combined$twas_cv_performance[[context]][[model]][, rsq_option])) {
          return(model)
        }
      }))
      if (length(available_models) <= 0) {
        message(paste0("No model provided TWAS cross validation performance metrics information at context ", context, ". "))
        return(NULL)
      }
      for (model in available_models) {
        model_data <- twas_data_combined$twas_cv_performance[[context]][[model]]
        if (model_data[, rsq_option] >= best_rsq & model_data[, colnames(model_data)[which(colnames(model_data) %in% rsq_pval_option)]] < rsq_pval_cutoff) {
          best_rsq <- model_data[, rsq_option]
          selected_model <- model
        }
      }
      if (is.null(selected_model)) {
        message(paste0(
          "No model has p-value < ", rsq_pval_cutoff, " and r2 >= ", rsq_cutoff, ", skipping context ", context,
          " at region ", unique(twas_data_combined$molecular_id), ". "
        ))
        return(list(selected_model = c("context_non_imputable"), is_imputable = FALSE)) # No significant model found
      } else {
        selected_model <- unlist(strsplit(selected_model, "_performance"))
        message(paste0("The selected best performing model for context ", context, " at region ", twas_data_combined$molecular_id, " is ", selected_model, ". "))
        return(list(selected_model = selected_model, is_imputable = TRUE))
      }
    })
    names(model_selection) <- names(twas_data_combined$weights)
    return(model_selection)
  }

  # Step 1: TWAS and MR analysis for all methods for imputable gene
  rsq_option <- match.arg(rsq_option)
  # filter events 
  if (!is.null(event_filters)){
    for (weight_db in names(twas_weights_data)){
      contexts <- names(twas_weights_data[[weight_db]]$weights)
      filtered_events <- filter_molecular_events(contexts, event_filters, remove_all_group=TRUE)
      if (length(filtered_events!=0)){
          for (db in names(twas_weights_data[[weight_db]])){
             twas_weights_data[[weight_db]][[db]] <- twas_weights_data[[weight_db]][[db]][filtered_events]
          }   
      } else {
          twas_weights_data[[weight_db]] <- NULL
      }
    }
  }
  if (length(twas_weights_data)==0) {
    return(list(NULL))
  }

  # harmonize twas weights and gwas sumstats against LD
  twas_data_qced_result <- harmonize_twas(twas_weights_data, ld_meta_file_path, gwas_meta_file,
                                          ld_reference_sample_size = ld_reference_sample_size,
                                          column_file_path = column_file_path, comment_string = comment_string)
  twas_results_db <- lapply(names(twas_weights_data), function(weight_db) {
    twas_weights_data[[weight_db]][["molecular_id"]] <- weight_db
    twas_data_qced <- twas_data_qced_result$twas_data_qced
    if (length(twas_data_qced[[weight_db]]) == 0 | is.null(twas_data_qced[[weight_db]])) {
      warning(paste0("No data harmonized for ", weight_db, ". Returning NULL for TWAS result for this region."))
      return(NULL)
    }
    if (rsq_cutoff > 0) {
      message("Selecting the best model based on criteria...")
      best_model_selection <- pick_best_model(
        twas_weights_data[[weight_db]],
        rsq_cutoff = rsq_cutoff,
        rsq_pval_cutoff = rsq_pval_cutoff,
        rsq_option = rsq_option,
        rsq_pval_option = rsq_pval_option
      )
      twas_data_qced[[weight_db]][["model_selection"]] <- setNames(best_model_selection, names(twas_weights_data[[weight_db]]$weights))
    } else {
      message("Skipping best model selection. Assigning NA of model_selection to all weights.")
      twas_data_qced[[weight_db]][["model_selection"]] <- setNames(
        rep(NA, length(names(twas_weights_data[[weight_db]]$weights))),
        names(twas_weights_data[[weight_db]]$weights)
      )
    }
    if (!"data_type" %in% names(twas_weights_data[[weight_db]])) {
      twas_data_qced[[weight_db]][["data_type"]] <- setNames(rep(
        list(NA),
        length(names(twas_weights_data[[weight_db]]$weights))
      ), names(twas_weights_data[[weight_db]]$weights))
    }
    if (length(weight_db) < 1) stop(paste0("No data harmonized for ", weight_db, ". "))
    contexts <- names(twas_data_qced[[weight_db]][["weights_qced"]])
    gwas_studies <- names(twas_data_qced[[weight_db]][["gwas_qced"]])

    # Combined loop for TWAS and MR analysis
    mr_cols <- c("gene_name", "num_CS", "num_IV", "cpip", "meta_eff", "se_meta_eff", "meta_pval", "Q", "Q_pval", "I2")

    # Nested lapply for contexts and gwas studies
    twas_gene_results <- lapply(contexts, function(context) {
      study_results <- lapply(gwas_studies, function(study) {
        twas_variants <- Reduce(intersect, list(rownames(twas_data_qced[[weight_db]][["weights_qced"]][[context]][[study]][["weights"]]), 
          twas_data_qced[[weight_db]][["variant_names"]][[context]][[study]],
          twas_data_qced[[weight_db]][["gwas_qced"]][[study]]$variant_id)
        )
        if (length(twas_variants) == 0) {
          return(list(twas_rs_df = data.frame(), mr_rs_df = data.frame()))
        }
        # twas analysis
        twas_rs <- twas_analysis(
          twas_data_qced[[weight_db]][["weights_qced"]][[context]][[study]][["weights"]], twas_data_qced[[weight_db]][["gwas_qced"]][[study]],
          twas_data_qced[[weight_db]][["LD"]], twas_variants
        )
        if (is.null(twas_rs)) {
          return(list(twas_rs_df = data.frame(), mr_rs_df = data.frame()))
        }
        twas_rs_df <- build_twas_score_row(twas_rs, weight_db, context, study)
        # MR analysis
        if (!is.null(twas_weights_data[[weight_db]]$susie_results) &&
          any(na.omit(twas_rs_df$twas_pval) < mr_pval_cutoff) &&
          "top_loci" %in% names(twas_weights_data[[weight_db]]$susie_results[[context]])) {
          if (!"effect_allele_frequency" %in% colnames(twas_data_qced[[weight_db]][["gwas_qced"]][[study]])) {
            warning(paste0("skip MR for ", weight_db, " for ", study, ", the effect_allele_frequency information is not available."))
            return(list(twas_rs_df = twas_rs_df, mr_rs_df = data.frame()))
          }
          combined_ld_meta_df <- twas_data_qced_result$ref_panel
          mr_formatted_input <- mr_format(twas_weights_data[[weight_db]], context, twas_data_qced[[weight_db]][["gwas_qced"]][[study]],
            coverage = mr_coverage_column, run_allele_qc = TRUE, molecular_name_obj = c("molecular_id"), ld_meta_df = combined_ld_meta_df
          )
          if (all(is.na(mr_formatted_input$bhat_y))) {
            # FIXME: after updating gwas beta and se NA problem, mr analysis will be restored
            mr_rs_df <- as.data.frame(matrix(rep(NA, length(mr_cols)), nrow = 1))
            colnames(mr_rs_df) <- mr_cols
          } else {
            mr_rs_df <- as.data.frame(mr_analysis(mr_formatted_input, cpip_cutoff = 0.1))
          }
        } else {
          mr_rs_df <- as.data.frame(matrix(rep(NA, length(mr_cols)), nrow = 1))
          colnames(mr_rs_df) <- mr_cols
        }
        mr_rs_df$context <- context
        mr_rs_df$gwas_study <- study
        mr_rs_df$gene_name <- weight_db
        return(list(twas_rs_df = twas_rs_df, mr_rs_df = mr_rs_df))
      })
      twas_context_table <- do.call(rbind, lapply(study_results, function(x) x$twas_rs_df))
      mr_context_table <- do.call(rbind, lapply(study_results, function(x) x$mr_rs_df))
      return(list(twas_context_table = twas_context_table, mr_context_table = mr_context_table))
    })
    twas_data_qced[[weight_db]][["LD"]] <- NULL
    twas_weights_data[[weight_db]] <- NULL
    twas_gene_table <- do.call(rbind, lapply(twas_gene_results, function(x) x$twas_context_table))
    mr_gene_table <- do.call(rbind, lapply(twas_gene_results, function(x) x$mr_context_table))
    return(list(twas_table = twas_gene_table, twas_data_qced = twas_data_qced[weight_db], mr_result = mr_gene_table))
  })
  rm(twas_data_qced_result)
  gc()
  twas_results_db <- twas_results_db[!sapply(twas_results_db, function(x) is.null(x) || (is.list(x) && all(sapply(x, is.null))))]
  if (length(twas_results_db) == 0) {
    return(list(NULL))
  }
  twas_results_table <- do.call(rbind, lapply(twas_results_db, function(x) x$twas_table))
  mr_results <- do.call(rbind, lapply(twas_results_db, function(x) x$mr_result))
  twas_data <- do.call(c, lapply(twas_results_db, function(x) x$twas_data_qced))
  # snp_info <- do.call(c, lapply(twas_results_db, function(x) x$snp_info))
  rm(twas_results_db)
  gc()

  # Step 2: Summarize and merge twas cv results and region information for all methods for all contexts for imputable genes.
  twas_table <- do.call(rbind, lapply(names(twas_data), function(molecular_id) {
    contexts <- names(twas_weights_data[[molecular_id]]$weights)
    # merge twas_cv information for same gene across all weight db files, loop through each context for all methods
    gene_table <- do.call(rbind, lapply(contexts, function(context) {
      methods <- sub("_[^_]+$", "", names(twas_weights_data[[molecular_id]]$twas_cv_performance[[context]]))
      is_imputable <- twas_data[[molecular_id]][["model_selection"]][[context]]$is_imputable
      selected_method <- twas_data[[molecular_id]][["model_selection"]][[context]]$selected_model
      if (is.null(selected_method)) selected_method <- NA
      is_selected_method <- ifelse(methods == selected_method, TRUE, FALSE)

      cv_rsqs <- sapply(twas_weights_data[[molecular_id]]$twas_cv_performance[[context]], function(x) x[, rsq_option])
      cv_pvals <- sapply(twas_weights_data[[molecular_id]]$twas_cv_performance[[context]], function(x) x[, colnames(x)[which(colnames(x) %in% rsq_pval_option)]])

      context_table <- data.frame(
        context = context, method = methods,
        is_imputable = is_imputable,
        is_selected_method = is_selected_method,
        rsq_cv = cv_rsqs, pval_cv = cv_pvals,
        type = twas_weights_data[[molecular_id]][["data_type"]][[context]]
      )
      return(context_table)
    }))
    gene_table$molecular_id <- molecular_id
    return(gene_table)
  }))
  twas_table$chr <- as.integer(strip_chr_prefix(gsub("\\_.*", "", region_block)))
  twas_table$block <- region_block

  # Step 3. merge twas result table and twas input into twas_data to output
  colname_ordered <- c("chr", "molecular_id", "context", "gwas_study", "method", "is_imputable", "is_selected_method", "rsq_cv", "pval_cv", "twas_z", "twas_pval", "type", "block")
  if (nrow(twas_results_table) == 0) {
    return(list(twas_result = NULL, twas_data = NULL, mr_result = NULL))
  }
  twas_table <- merge(twas_table, twas_results_table, by = c("molecular_id", "context", "method"))
  twas_table <- twas_table[twas_table$is_imputable, , drop = FALSE]
  if (output_twas_data & nrow(twas_table) > 0) {
    twas_data_subset <- format_twas_data(twas_data, twas_table)
    # if (!is.null(twas_data_subset)) twas_data_subset$snp_info <- snp_info
  } else {
    twas_data_subset <- NULL
  }
  return(list(twas_result = twas_table[, colname_ordered], twas_data = twas_data_subset, mr_result = mr_results))
}

#' Calculate TWAS z-score and p-value
#'
#' This function calculates the TWAS z-score and p-value given the weights, z-scores,
#' and optionally the correlation matrix (R) or the genotype matrix (X).
#'
#' @param weights A numeric vector of weights.
#' @param z A numeric vector of z-scores.
#' @param R An optional correlation matrix. If not provided, it will be calculated from the genotype matrix X.
#' @param X An optional genotype matrix. If R is not provided, X must be supplied to calculate the correlation matrix.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item z: The TWAS z-score.
#'   \item pval: The corresponding p-value.
#' }
#'
#' @importFrom stats cor pchisq
#'
#' @export
twas_z <- function(weights, z, R = NULL, X = NULL) {
  # Check that weights and z-scores have the same length
  if (length(weights) != length(z)) {
    stop("Weights and z-scores must have the same length.")
  }

  if (is.null(R)) R <- compute_LD(X)

  stat <- t(weights) %*% z
  denom <- t(weights) %*% R %*% weights
  zscore <- stat / sqrt(denom)
  pval <- pchisq(zscore * zscore, 1, lower.tail = FALSE)

  return(list(z = zscore, pval = pval))
}

#' Multi-condition TWAS joint test
#'
#' This function performs a multi-condition TWAS joint test using the GBJ method.
#' It assumes that the input genotype matrix (X) is standardized.
#'
#' @param R An optional correlation matrix. If not provided, it will be calculated from the genotype matrix X.
#' @param X An optional genotype matrix. If R is not provided, X must be supplied to calculate the correlation matrix.
#' @param weights A matrix of weights, where each column corresponds to a different condition.
#' @param z A vector of GWAS z-scores.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item Z: A matrix of TWAS z-scores and p-values for each condition.
#'   \item GBJ: The result of the GBJ test.
#' }
#'
#' @importFrom stats cor pnorm
#' @export
twas_joint_z <- function(weights, z, R = NULL, X = NULL) {
  # Make sure GBJ is installed
  if (!requireNamespace("GBJ", quietly = TRUE)) {
    stop("To use this function, please install GBJ: https://cran.r-project.org/web/packages/GBJ/index.html")
  }
  # Check that weights and z-scores have the same number of rows
  if (nrow(weights) != length(z)) {
    stop("Number of rows in weights must match the length of z-scores.")
  }
  if (is.null(R)) R <- compute_LD(X)
  idx <- which(rownames(R) %in% rownames(weights))
  D <- R[idx, idx]
  cov_y <- crossprod(weights, D) %*% weights
  y_sd <- sqrt(diag(cov_y))
  x_sd <- rep(1, nrow(weights)) # Assuming X is standardized

  # Get gamma matrix MxM (snp x snp)
  g <- lapply(colnames(weights), function(x) {
    gm <- diag(x_sd / y_sd[x], length(x_sd), length(x_sd))
    return(gm)
  })
  names(g) <- colnames(weights)

  ######### Get TWAS - Z statistics & P-value, GBJ test ########
  z_matrix <- do.call(rbind, lapply(colnames(weights), function(x) {
    Zi <- crossprod(weights[, x], g[[x]]) %*% as.numeric(z)
    pval <- 2 * pnorm(abs(Zi), lower.tail = FALSE)
    Zp <- c(Zi, pval)
    names(Zp) <- c("Z", "pval")
    return(Zp)
  }))
  rownames(z_matrix) <- colnames(weights)

  # GBJ test
  lam <- matrix(rep(NA, ncol(weights) * nrow(weights)), nrow = ncol(weights))
  rownames(lam) <- colnames(weights)
  for (p in colnames(weights)) {
    la <- as.matrix(weights[, p] %*% g[[p]])
    lam[p, ] <- la
  }
  sig <- tcrossprod((lam %*% D), lam)
  gbj <- GBJ::GBJ(test_stats = z_matrix[, 1], cor_mat = sig)

  rs <- list("Z" = z_matrix, "GBJ" = gbj)
  return(rs)
}

#' TWAS Analysis
#'
#' Performs TWAS analysis using the provided weights matrix, GWAS summary statistics database,
#' and LD matrix. It extracts the necessary GWAS summary statistics and LD matrix based on the
#' specified variants and computes the z-score and p-value for each gene.
#'
#' @param weights_matrix A matrix containing weights for all methods.
#' @param gwas_sumstats_db A data frame containing the GWAS summary statistics.
#' @param LD_matrix A matrix representing linkage disequilibrium between variants.
#' @param extract_variants_objs A vector of variant identifiers to extract from the GWAS and LD matrix.
#'
#' @return A list with TWAS z-scores and p-values across four methods for each gene.
#' @export
twas_analysis <- function(weights_matrix, gwas_sumstats_db, LD_matrix, extract_variants_objs) {
  # Extract gwas_sumstats
  gwas_sumstats_subset <- gwas_sumstats_db[match(extract_variants_objs, gwas_sumstats_db$variant_id), ]
  # Validate that the GWAS subset is not empty
  if (nrow(gwas_sumstats_subset) == 0 | all(is.na(gwas_sumstats_subset))) {
    warning("No GWAS summary statistics found for the specified variants.")
    return(NULL)
  }      
  # Check if extract_variants_objs are in the rownames of LD_matrix
  valid_indices <- extract_variants_objs %in% rownames(LD_matrix)
  if (!any(valid_indices)) {
    warning("None of the specified variants are present in the LD matrix. Skipping this context.")
    return(NULL)
  }    
  # Extract only the valid indices from extract_variants_objs
  valid_variants_objs <- extract_variants_objs[valid_indices]
  # Extract LD_matrix subset using valid indices
  LD_matrix_subset <- LD_matrix[valid_variants_objs, valid_variants_objs]
  # Extract weight matrix subset using valid indices
  weights_matrix <- weights_matrix[valid_variants_objs, , drop = FALSE]
  # Caculate the z score and pvalue of each gene
  twas_z_pval <- apply(
    as.matrix(weights_matrix), 2,
    function(x) twas_z(x, gwas_sumstats_subset$z, R = LD_matrix_subset)
  )
  return(twas_z_pval)
}
