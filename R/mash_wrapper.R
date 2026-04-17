# Filter rows of a z-score matrix by significance p-value cutoff.
# Returns integer indices of rows where any |z| exceeds the threshold.
# @noRd
filter_by_significance <- function(z_matrix, sig_p_cutoff) {
  z_threshold <- sqrt(qchisq(sig_p_cutoff, df = 1, lower.tail = FALSE))
  which(apply(z_matrix, 1, function(row) any(abs(row) >= z_threshold)))
}

#' @importFrom vroom vroom
#' @export
filter_invalid_summary_stat <- function(dat_list, bhat = NULL, sbhat = NULL, z = NULL, btoz = FALSE, sig_p_cutoff = 1E-6, filter_by_missing_rate = 0.2) {
  replace_values <- function(df, replace_with) {
    df <- df %>%
      mutate(across(everything(), as.numeric)) %>%
      mutate(across(everything(), ~ replace(., is.nan(.) | is.infinite(.) | is.na(.), replace_with)))
  }
  # Function to process bhat, sbhat
  if (!is.null(bhat) && !is.null(sbhat) && all(c(bhat, sbhat) %in% names(dat_list))) {
    # If the element is a list with 'bhat' and 'sbhat'
    if (!is.null(dat_list[[bhat]]) && !is.null(dat_list[[sbhat]])) {
      dat_list[[bhat]] <- as.matrix(replace_values(dat_list[[bhat]], 0))
      dat_list[[sbhat]] <- as.matrix(replace_values(dat_list[[sbhat]], 1000))
      if (("null.b" %in% names(dat_list)) || ("random.b" %in% names(dat_list))) {
        if (!is.null(filter_by_missing_rate)) {
          proportion_nonzero <- apply(dat_list[[bhat]], 1, function(row) {
            mean(row != 0)
          })
          dat_list[[bhat]] <- dat_list[[bhat]][proportion_nonzero >= filter_by_missing_rate, ]
          dat_list[[sbhat]] <- dat_list[[sbhat]][proportion_nonzero >= filter_by_missing_rate, ]
        }
      }
    }
  }
  # Function to filter strong signal using z score
  if (btoz) {
    if (any(grepl("\\.b$", bhat)) | any(grepl("\\.s$", sbhat))) {
      condition <- sub("\\.b$", "", bhat)
      if (!is.null(dat_list[[bhat]]) && !is.null(dat_list[[sbhat]])) {
        dat_list[[paste0(condition, ".z")]] <- as.matrix(dat_list[[bhat]] / dat_list[[sbhat]])
      } else {
        dat_list[paste0(condition, ".z")] <- list(NULL)
      }
    } else {
      if (!is.null(dat_list[[bhat]]) && !is.null(dat_list[[sbhat]])) {
        dat_list[["z"]] <- as.matrix(dat_list[[bhat]] / dat_list[[sbhat]])
      } else {
        dat_list["z"] <- list(NULL)
      }
    }
    if ("strong.z" %in% names(dat_list)) {
      if (!is.null(sig_p_cutoff)) {
        keep_index <- filter_by_significance(dat_list$strong.z, sig_p_cutoff)
        dat_list[["strong.z"]] <- dat_list$strong.z[keep_index, ]
        dat_list[["strong.b"]] <- dat_list$strong.b[keep_index, ]
        dat_list[["strong.s"]] <- dat_list$strong.s[keep_index, ]
      }
    }
  }
  # Function to process z-scores and filter directly
  if (!is.null(z)) {
    process_z <- function(z_data) {
      z_data <- as.matrix(replace_values(z_data, 0))

      if (!is.null(filter_by_missing_rate)) {
        proportion_nonzero <- apply(z_data, 1, function(row) mean(row != 0))
        z_data <- z_data[proportion_nonzero >= filter_by_missing_rate, , drop = FALSE]
      }

      return(z_data)
    }

    # Process each component if it exists
    for (comp in c("strong", "random", "null")) {
      if (!is.null(dat_list[[comp]]) && !is.null(dat_list[[comp]]$z)) {
        dat_list[[comp]]$z <- process_z(dat_list[[comp]]$z)
      }
    }

    # Apply significance cutoff to strong signals if applicable
    if (!is.null(dat_list$strong) && !is.null(dat_list$strong$z) && !is.null(sig_p_cutoff)) {
      keep_index <- filter_by_significance(dat_list$strong$z, sig_p_cutoff)
      dat_list$strong$z <- dat_list$strong$z[keep_index, , drop = FALSE]
    }
  }

  return(dat_list)
}

#' @export
filter_mixture_components <- function(conditions_to_keep, U, w = NULL, w_cutoff = 1e-04) {
  # Identify conditions not to keep (to be removed)
  conditions_to_filter <- setdiff(colnames(U[[1]]), conditions_to_keep)
  sum_w <- sum(w) # Original total sum of weights

  # Filter U by removing unwanted phenotypes (conditions)
  U <- lapply(U, function(mat, to_keep) {
    missing_conditions <- setdiff(to_keep, colnames(mat))
    if (length(missing_conditions) > 0) {
      stop(paste("Condition(s)", paste(missing_conditions,
        collapse = ", "
      ), "not found in matrix"))
    }
    mat[to_keep, to_keep] # Keep only relevant conditions
  }, conditions_to_keep)

  # Remove matrices where all values are zero or weight is below cutoff
  keep_names <- names(purrr::keep(U, function(mat) !all(mat == 0)))
  if (!is.null(w)) {
    keep_names <- intersect(keep_names, names(w[w >= w_cutoff]))
  }
  U <- U[keep_names]
  if (!is.null(w)) {
    w <- w[keep_names]
  }

  # Note: Matrices in U may contain very small values on the diagonal
  # even when contexts are not present, due to EM algorithm adjustments.
  # This makes the matrix not exactly zero, so it won't be removed even though it may not
  # have strong context relevance. This behavior arises because the algorithm attempts to
  # ensure matrices are full-rank, slightly changing initial values.

  # We cannot simply remove diagonal matrices as signals on the diagonal can be strong and relevant.
  # So we manually remove the U components that are driven by non-relevant contexts.
  U[conditions_to_filter] <- NULL
  w <- w[!names(w) %in% conditions_to_filter]

  # Recalculate the sum of remaining weights
  sum_w_new <- sum(w)

  # Adjust weights to maintain the original sum_w
  w <- (w / sum_w_new) * sum_w

  message(paste(length(U), "components of matrices remained after filtering."))

  return(list(U = U, w = w))
}
# This function extracts tensorQTL results for given region for multiple
# summary statistics files
#' @export
load_multitrait_tensorqtl_sumstat <- function(
    sumstats_paths, region, gene = NULL,
    trait_names = NULL, top_loci = FALSE, filter_file = NULL, remove_any_missing = FALSE,
    max_rows_selected = 300, nan_remove = TRUE) {
  if (!is.vector(sumstats_paths) || !all(file.exists(sumstats_paths))) {
    stop("sumstats_paths must be a vector of existing file paths.")
  }
  if (!is.character(region) || length(region) != 1) {
    stop("region must be a single character string.")
  }
  if (!is.character(trait_names)) {
    stop("trait_names must be a vector of character strings.")
  }


  extract_tensorqtl_data <- function(path, region) {
    tabix_region(path, region) %>%
      # first four columns are 'chrom','pos','alt','ref'
      mutate(variants = paste(.[[1]], .[[2]], .[[3]], .[[4]], sep = ":")) %>%
      distinct(variants, molecular_trait_id, .keep_all = TRUE)
  }


  merge_matrices <- function(matrix_list, value_column, id_column = "variants",
                             remove_any_missing = FALSE) {
    # Convert matrices to data frames
    df_list <- lapply(seq_along(matrix_list), function(i) {
      df <- as.data.frame(matrix_list[[i]])
      df2 <- df[, c(id_column, value_column)]
      # Rename columns to avoid duplication
      colnames(df2) <- c(id_column, paste0(value_column, "_", i))
      return(df2)
    })

    # Iteratively merge the data frames
    merged_df <- Reduce(
      function(x, y) merge(x, y, by = id_column, all = TRUE),
      df_list
    )

    # Optionally, remove rows with any missing values
    if (remove_any_missing) {
      merged_df <- merged_df[complete.cases(merged_df), ]
    }
    return(merged_df)
  }

  split_variants_and_match <- function(variant, filter_file, max_rows_selected) {
    if (!file.exists(filter_file)) {
      stop("Filter file does not exist.")
    }

    # Split the variant vector into components
    variant_split <- strsplit(variant, ":")
    variant_df <- data.frame(chr = sapply(variant_split, `[`, 1), pos = sapply(
      variant_split,
      `[`, 2
    ), stringsAsFactors = FALSE)
    variant_df$pos <- as.numeric(variant_df$pos)

    # get the region of interest
    min_pos <- min(variant_df$pos)
    max_pos <- max(variant_df$pos)
    chrom <- unique(variant_df$chr)
    if (length(chrom) != 1) {
      stop("Variants are from multiple chromosomes. Cannot create a single range string.")
    }
    region <- paste0(chrom, ":", min_pos, "-", max_pos)
    ref_table <- tabix_region(filter_file, region)
    if (is.null(ref_table)) {
      stop("No variants in the region.")
    }
    colnames(ref_table)[1:2] <- c("#CHROM", "POS")
    if (!all(c("#CHROM", "POS") %in% colnames(ref_table))) {
      stop("Filter file must contain columns: #CHROM, POS.")
    }
    matched_indices <- which(variant_df$chr %in% ref_table$`#CHROM` & variant_df$pos %in%
      ref_table$POS)
    if (!is.null(max_rows_selected) && max_rows_selected > 0 && max_rows_selected <
      length(matched_indices)) {
      selected_rows <- sample(length(matched_indices), max_rows_selected)
      matched_indices <- matched_indices[selected_rows]
    }
    return(matched_indices)
  }

  Y <- lapply(sumstats_paths, function(x, region, gene) {
    out <- extract_tensorqtl_data(x, region)
    if (!is.null(gene)) {
      out <- out[which(out$molecular_trait_id %in% gene), ]
    }
    sorted <- out[order(-abs(out$beta / out$se)), c("variants", "molecular_trait_id")]
    top_v <- apply(sorted[1:2, ], 1, function(row) paste(row, collapse = "_")) # paste the variant (chr:pos:alt:ref) with gene_id with '_'
    out <- as.list(out)
    out$top_variants <- top_v
    return(out)
  }, region = region, gene = gene)

  ## Y is list of data frames where colnames(Y[[1]])=
  ## c('chrom','pos','alt','ref','variant_id','molecular_trait_id','start_distance',
  ## 'end_distance','af','ma_samples','ma_count','pvalue','beta','se','molecular_trait_object_id',
  ## 'n','variant''top_variants')

  ### The step below assigns condition names to the Y; in case the filename
  ### itself does not contain any condition names, users can input the
  ### condition names via assigning trait_names if trait_names left blank, then
  ### the file names will be assigned as trait_names, where if the text before
  ### '.' (extension) can differentiate the conditions, we will use the shorter
  ### names as trait_name
  if (is.null(trait_names)) {
    trait_names <- gsub("\\..*", "", basename(sumstats_paths)) # extract condition name that is listed before the first appearance of '.'

    if (length(trait_names[duplicated(trait_names)]) >= 1) {
      trait_names <- basename(sumstats_paths)
    }
  }
  names(Y) <- trait_names

  bhat <- merge_matrices(Y,
    value_column = "beta", id_column = c("variants", "molecular_trait_id"),
    remove_any_missing
  )
  sbhat <- merge_matrices(Y,
    value_column = "se", id_column = c("variants", "molecular_trait_id"),
    remove_any_missing
  )
  out <- list(bhat = bhat, sbhat = sbhat)

  # Check if variants are the same in both bhat and sbhat
  if (!identical(out$bhat$variants, out$sbhat$variants)) {
    stop("Error: Variants in bhat and sbhat are not the same.")
  }

  var_idx <- 1:nrow(out$bhat)

  # match with filter_file
  if (!is.null(filter_file)) {
    if (!file.exists(filter_file)) {
      stop("Filter file does not exist.")
    }
    variants <- paste0(out$bhat$variants, "_", out$bhat$molecular_trait_id)
    var_idx <- split_variants_and_match(out$bhat$variants, filter_file, max_rows_selected)
  }

  if (top_loci) {
    union_top_loci <- unique(unlist(lapply(Y, function(item) item$top_variants)))
    var_idx <- which(variants %in% union_top_loci) # var_idx may end up empty if max_rows_selected number too small
  }

  # Extract only subset of data
  variants <- paste0(out$bhat$variants[var_idx], "_", out$bhat$molecular_trait_id[var_idx])
  out$bhat <- out$bhat[var_idx, ]
  out$sbhat <- out$sbhat[var_idx, ]

  if (nan_remove) {
    out <- filter_invalid_summary_stat(out, bhat = "beta", sbhat = "se", nan_remove)
  }

  rownames(out$bhat) <- rownames(out$sbhat) <- variants
  colnames(out$bhat)[which(startsWith(colnames(out$bhat), "beta"))] <- colnames(out$sbhat)[which(startsWith(
    colnames(out$sbhat),
    "se"
  ))] <- trait_names
  out$region <- region
  out$top_variants <- lapply(Y, function(x) x$top_variants)
  return(out)
}

merge_susie_cs <- function(susie_fit, coverage = "cs_coverage_0.95", complementary = FALSE) {
  # Initialize an empty list for the combined_sets
  combined_sets <- list()

  # Identify variant IDs that are associated with more than one credible set
  identify_overlap_sets <- function(variants_sets_and_pips_list) {
    overlap_sets <- list()
    for (variant_id in names(variants_sets_and_pips_list)) {
      sets <- variants_sets_and_pips_list[[variant_id]][["sets"]]
      if (length(sets) > 1) {
        overlap_sets[[variant_id]] <- sets
      }
    }
    return(overlap_sets)
  }
  # Merge overlapping credible sets using graph connected components
  merge_and_update_overlap_sets <- function(variants_sets_and_pips_list, overlap_sets) {
    # Collect all unique set names from overlap_sets
    all_sets <- unique(unlist(overlap_sets))
    if (length(all_sets) == 0) return(list())

    # Build edges: each overlap_sets entry links its sets together
    edges <- unlist(lapply(overlap_sets, function(sets) {
      if (length(sets) < 2) return(character(0))
      pairs <- combn(sets, 2)
      as.character(pairs)
    }))

    g <- igraph::make_graph(edges, directed = FALSE)
    comp <- igraph::components(g)

    # Build mapping: each set name -> comma-separated merged component label
    set_name_map <- list()
    for (i in seq_len(comp$no)) {
      members <- names(comp$membership[comp$membership == i])
      label <- paste(sort(members), collapse = ",")
      for (s in members) {
        set_name_map[[s]] <- label
      }
    }

    # Update each variant's credible set names
    updated_credible_sets <- lapply(
      stats::setNames(names(variants_sets_and_pips_list), names(variants_sets_and_pips_list)),
      function(variant_id) {
        current_sets <- variants_sets_and_pips_list[[variant_id]][["sets"]]
        mapped <- intersect(current_sets, names(set_name_map))
        if (length(mapped) > 0) {
          set_name_map[[mapped[1]]]
        } else {
          paste(sort(unique(current_sets)), collapse = ",")
        }
      }
    )
    return(updated_credible_sets)
  }
  # Loop through each condition and their credible sets
  extract_top_loci <- function(susie_fit, complementary, coverage) {
    # Build a flat data frame of (variant_id, pip, set_name) across all conditions
    cond_names <- names(susie_fit[[1]])
    rows <- purrr::map_dfr(seq_along(cond_names), function(i) {
      cond_data <- susie_fit[[1]][[i]]
      top_loci <- cond_data[["top_loci"]]
      if (is.null(top_loci) || nrow(top_loci) == 0) return(NULL)

      if (!complementary) {
        set_num <- unique(get_nested_element(cond_data, c("top_loci", coverage)))
        set_num <- set_num[set_num != 0]
      } else {
        set_num <- 0
      }
      if (length(set_num) == 0) return(NULL)

      purrr::map_dfr(set_num, function(sn) {
        top_loci %>%
          filter(!!sym(coverage) == sn) %>%
          select(variant_id, pip) %>%
          mutate(set_name = paste0("cs_", i, "_", sn))
      })
    })

    if (is.null(rows) || nrow(rows) == 0) return(list())

    # Aggregate by variant_id preserving first-seen order
    seen_order <- unique(rows$variant_id)
    split_rows <- split(rows, factor(rows$variant_id, levels = seen_order))
    lapply(split_rows, function(df) {
      list(sets = df$set_name, pips = df$pip)
    })
  }

  combine_top_loci <- function(extracted_result) {
    # Compute overlap sets once, outside the per-variant loop
    overlap_sets <- identify_overlap_sets(extracted_result)
    has_overlaps <- length(overlap_sets) != 0
    merged_sets <- if (has_overlaps) {
      merge_and_update_overlap_sets(extracted_result, overlap_sets = overlap_sets)
    } else {
      NULL
    }

    top_loci_df <- do.call(rbind, lapply(names(extracted_result), function(variant_id) {
      max_pip <- max(unlist(extracted_result[[variant_id]]$pips))
      median_pip <- median(unlist(extracted_result[[variant_id]]$pips))
      credible_set_names <- if (has_overlaps) {
        merged_sets[[variant_id]]
      } else {
        paste(sort(unique(unlist(extracted_result[[variant_id]]$sets))), collapse = ",")
      }
      data.frame(
        variant_id = variant_id, credible_set_names = credible_set_names,
        max_pip = max_pip, median_pip = median_pip, stringsAsFactors = FALSE
      )
    }))
    return(top_loci_df)
  }

  extracted_top_loci <- extract_top_loci(susie_fit, complementary, coverage = coverage)
  combined_top_loci_df <- combine_top_loci(extracted_top_loci)
  # Clean up row names and make sure variant_id is unique
  combined_top_loci_df <- combined_top_loci_df[!duplicated(combined_top_loci_df$variant_id), ]
  rownames(combined_top_loci_df) <- NULL # Clean up row names
  return(combined_top_loci_df)
}

#' @export
load_multitrait_R_sumstat <- function(susie_fit, sumstats_db, coverage = NULL, extract_inf = "z", top_loci = FALSE, filter_file = NULL, exclude_condition = NULL, ld_meta_file = NULL, remove_any_missing = TRUE, max_rows_selected = 300, nan_remove = FALSE, condition_filter = FALSE) {
  # Internal recursive filtering function
  filter_nested_list <- function(input_list, valid_conditions) {
    if (is.null(valid_conditions) || length(valid_conditions) == 0) {
      return(input_list)
    }

    filtered_list <- list()

    # Recursively process list
    for (name in names(input_list)) {
      # Check if current name matches any valid condition
      if (name %in% valid_conditions) {
        filtered_list[[name]] <- input_list[[name]]
      } else {
        # Recursively check nested lists
        if (is.list(input_list[[name]])) {
          nested_result <- filter_nested_list(input_list[[name]], valid_conditions)
          if (length(nested_result) > 0) {
            filtered_list[[name]] <- nested_result
          }
        }
      }
    }

    return(filtered_list)
  }

  # Apply condition filtering before processing
  if (!is.null(condition_filter)) {
    # Convert condition_filter to character vector if it's not already
    if (!is.character(condition_filter)) {
      condition_filter <- as.character(condition_filter)
    }

    # Split if comma-separated string
    if (length(condition_filter) == 1 && grepl(",", condition_filter)) {
      condition_filter <- strsplit(condition_filter, ",")[[1]]
    }

    # Trim whitespace
    condition_filter <- trimws(condition_filter)

    # Filter susie_fit and sumstats_db
    susie_fit <- filter_nested_list(susie_fit, condition_filter)
    sumstats_db <- filter_nested_list(sumstats_db, condition_filter)
  }


  split_variants_and_match <- function(variant, filter_file, max_rows_selected) {
    if (!file.exists(filter_file)) {
      stop("Filter file does not exist.")
    }

    # Split the variant vector into components
    variant_df <- parse_variant_id(variant)
    conv <- attr(variant_df, "convention")

    # get the region of interest
    min_pos <- min(variant_df$pos)
    max_pos <- max(variant_df$pos)
    chrom <- unique(variant_df$chrom)
    if (length(chrom) != 1) {
      stop("Variants are from multiple chromosomes. Cannot create a single range string.")
    }
    # Reconstruct chrom string with original convention for tabix
    chrom_str <- if (conv$has_chr) paste0("chr", chrom) else as.character(chrom)
    region <- paste0(chrom_str, ":", min_pos, "-", max_pos)
    ref_table <- tabix_region(filter_file, region)
    if (is.null(ref_table)) {
      stop("No variants in the region.")
    }
    colnames(ref_table)[1:2] <- c("#CHROM", "POS")
    if (!all(c("#CHROM", "POS") %in% colnames(ref_table))) {
      stop("Filter file must contain columns: #CHROM, POS.")
    }
    ref_chrom <- as.integer(strip_chr_prefix(ref_table$`#CHROM`))
    matched_indices <- which(variant_df$chrom %in% ref_chrom & variant_df$pos %in% ref_table$POS)
    if (!is.null(max_rows_selected) && max_rows_selected > 0 && max_rows_selected < length(matched_indices)) {
      selected_rows <- sample(length(matched_indices), max_rows_selected)
      matched_indices <- matched_indices[selected_rows]
    }
    return(matched_indices)
  }


  results <- lapply(sumstats_db[[1]], function(data) extract_flatten_sumstats_from_nested(data,extract_inf))
  trait_names <- names(results)
  z_scores <- merge_matrices(results, value_column = extract_inf, ld_meta_file, id_column = "variants", remove_any_missing)
  out <- list(z = z_scores)
  var_idx <- 1:nrow(out[[1]])
  if (!is.null(filter_file)) {
    variants <- out[[1]]$variants
    var_idx <- split_variants_and_match(variants, filter_file, max_rows_selected)
  }

  if (top_loci) {
    union_top_loci <- merge_susie_cs(susie_fit, coverage)
    if (!is.null(union_top_loci)) {
      strong_signal_df <- union_top_loci %>%
        group_by(credible_set_names) %>%
        filter(median_pip == max(median_pip)) %>%
        slice(1) %>%
        ungroup()
      var_idx <- which(out[[1]]$variants %in% strong_signal_df$variant_id)
    } else {
      var_idx <- NULL
    }
  }

  # Extract only subset of data
  variants <- out[[1]]$variants[var_idx]
  for (key in names(out)) {
    out[[key]] <- out[[key]][var_idx, , drop = FALSE]
    rownames(out[[key]]) <- variants
    colnames(out[[key]])[2:ncol(out[[key]])] <- trait_names
    out[[key]] <- out[[key]][, -which(names(out[[key]]) == "variants"), drop = FALSE]
  }
  out$region <- names(susie_fit)

  if (!is.null(exclude_condition) && length(exclude_condition) != 0) {
    for (key in setdiff(names(out), "region")) {
      if (all(exclude_condition %in% colnames(out[[key]]))) {
        out[[key]] <- out[[key]][, -exclude_condition]
      } else {
        stop(paste("Error: exclude_condition are not present in", out$region))
      }
    }
  }

  return(out)
}

#' @export
mash_rand_null_sample <- function(dat, n_random, n_null, exclude_condition, seed = NULL) {
  # Function to extract one data set
  extract_one_data <- function(dat, n_random, n_null) {
    if (is.null(dat)) {
      return(NULL)
    }

    if ("z" %in% names(dat)) {
      abs_z <- abs(dat$z)
      z_data <- dat$z
    } else {
      abs_z <- abs(dat$bhat / dat$sbhat)
      z_data <- NULL
    }

    sample_idx <- 1:nrow(abs_z)
    random_idx <- sample(sample_idx, min(n_random, length(sample_idx)), replace = FALSE)

    if (!is.null(z_data)) {
      random <- list(z = z_data[random_idx, , drop = FALSE])
    } else {
      random <- list(
        bhat = dat$bhat[random_idx, , drop = FALSE],
        sbhat = dat$sbhat[random_idx, , drop = FALSE]
      )
    }

    null.id <- which(apply(abs_z, 1, max) < 2)
    if (length(null.id) == 0) {
      warning(paste("no variants are included in the null dataset because abs_z > 2 for all variants in", dat$region))
      null <- list()
    } else {
      if (length(null.id) < ncol(abs_z)) {
        warning(paste("not enough null data to estimate null correlation in", dat$region))
        null <- list()
      } else {
        null_idx <- sample(null.id, min(n_null, length(null.id)), replace = FALSE)
        if (!is.null(z_data)) {
          null <- list(z = z_data[null_idx, , drop = FALSE])
        } else {
          null <- list(
            bhat = dat$bhat[null_idx, , drop = FALSE],
            sbhat = dat$sbhat[null_idx, , drop = FALSE]
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

  if (length(exclude_condition) > 0) {
    cols_to_check <- if ("z" %in% names(dat)) "z" else "bhat"
    if (!all(exclude_condition %in% colnames(dat[[cols_to_check]]))) {
      stop(paste("Error: exclude_condition are not present in", dat$region))
    }
    for (key in intersect(names(dat), c("z", "bhat", "sbhat"))) {
      dat[[key]] <- dat[[key]][, -exclude_condition, drop = FALSE]
    }
  }

  result <- extract_one_data(dat, n_random, n_null)
  return(result)
}

#' @export
merge_mash_data <- function(res_data, one_data) {
  if (length(res_data) == 0 || is.null(res_data)) return(one_data)
  if (length(one_data) == 0 || is.null(one_data)) return(res_data)

  combined_data <- lapply(names(one_data), function(d) {
    od <- one_data[[d]]
    rd <- res_data[[d]]
    if (length(od) == 0 || is.null(od)) return(rd)
    if (is.null(rd) || length(rd) == 0) return(od)

    # bind_rows auto-aligns columns, filling missing with NA; replace with NaN
    rn_res <- rownames(as.data.frame(rd))
    rn_one <- rownames(as.data.frame(od))
    combined <- dplyr::bind_rows(as.data.frame(rd), as.data.frame(od))
    combined[is.na(combined)] <- NaN
    rn_all <- make.names(c(rn_res, rn_one), unique = TRUE)
    rownames(combined) <- rn_all
    combined
  })
  names(combined_data) <- names(one_data)
  return(combined_data)
}

#' @export
mash_pipeline <- function(mash_input, alpha, residual_correlation = NULL, n_pcs = NULL, set_seed = 999) {
  if (!requireNamespace("mashr", quietly = TRUE)) {
    stop("To use this function, please install mashr: https://cran.r-project.org/web/packages/mashr/index.html")
  }
  if (!requireNamespace("flashier", quietly = TRUE)) {
    stop("To use this function, please install flashier: https://github.com/willwerscheid/flashier")
  }
  set.seed(set_seed)
  if (length(mash_input$null.b) == 0 && length(mash_input$null.s) == 0) {
    if (!is.null(residual_correlation)) {
      vhat <- residual_correlation
    } else {
      condition_num <- ncol(mash_input$random.b)
      vhat <- diag(rep(1, condition_num))
    }
  } else {
    vhat <- mashr::estimate_null_correlation_simple(mashr::mash_set_data(mash_input$null.b,
      Shat = mash_input$null.s,
      alpha, zero_Bhat_Shat_reset = 1000
    ))
  }

  mash_data <- mashr::mash_set_data(mash_input$strong.b,
    Shat = mash_input$strong.s, V = vhat,
    alpha, zero_Bhat_Shat_reset = 1000
  )

  # Canonical covariance matrices
  U.can <- mashr::cov_canonical(mash_data)
  # PCA-based covariance matrices
  if (is.null(n_pcs)) {
    n_pcs <- ncol(mash_data$Bhat) - 1
  }
  U.pca <- mashr::cov_pca(mash_data, npc = n_pcs)
  # Flash-based covariance matrices (factor analysis)
  U.flash <- mashr::cov_flash(mash_data)
  # ED-based covariance matrices (initialized from all others)
  U.ed <- mashr::cov_ed(mash_data, Ulist_init = c(U.can, U.pca, U.flash))
  # Combine all covariance matrices
  U.all <- c(U.can, U.pca, U.flash, U.ed)

  # Fit mash to estimate mixture weights
  m <- mashr::mash(mash_data, Ulist = U.all, outputlevel = 1)
  w <- mashr::get_estimated_pi(m)

  return(list(U = U.all, w = w))
}

#' Merge a List of Matrices or Data Frames with Optional Allele Flipping
#'
#' @description
#' This function merges a list of matrices or data frames by a shared identifier column,
#' optionally aligning to a reference panel using allele QC procedures.
#'
#' @param matrix_list A named or unnamed list of data frames or matrices.
#' @param value_column Character string. The name of the column containing values to extract (e.g., z-scores or betas).
#' @param ref_panel Optional data frame. A reference panel for allele QC (must be compatible with `allele_qc`).
#' @param id_column Character string. The name of the column identifying variant IDs. Default is `"variants"`.
#' @param remove_any_missing Logical. If `TRUE`, rows with any missing values will be removed after merging.
#'
#' @return A data frame containing merged values, one column per dataset with suffix `_i`.
#' @examples
#' \dontrun{
#' merged <- merge_matrices(list(df1, df2), value_column = "variants", ref_panel = ref_df)
#' }
#' @import dplyr  
#' @export 

merge_sumstats_matrices <- function(matrix_list, value_column, ref_panel = NULL, ld_meta_file = NULL, id_column = "variants",
                             remove_any_missing = FALSE) {
    # Input validation
    if (!is.list(matrix_list) || length(matrix_list) == 0) {
      stop("matrix_list must be a non-empty list")
    }
    if (!is.character(value_column) || length(value_column) != 1) {
      stop("value_column must be a single string")
    }
    if (!is.character(id_column) || length(id_column) != 1) {
      stop("id_column must be a single string")
    }

    df_list <- lapply(seq_along(matrix_list), function(i) {
      tryCatch(
        {  
           # Step 1: Convert matrix to data frame and extract relevant columns
           df <- as.data.frame(matrix_list[[i]])
           if (!(id_column %in% colnames(df)) || !(value_column %in% colnames(df))) {
            stop(paste("Required columns", id_column, "or", value_column, "not found in dataset", i))
           }
           df2 <- df[, c(id_column, value_column)]
             if (!is.null(ld_meta_file)) {
            # Step 2: Split 'variants' to extract chromosomal info
            cohort_variants_df <- parse_variant_id(df2[, c(id_column)])
            # Step 3: Combine extracted chromosomal info with value column
            cohort_df <- cbind(cohort_variants_df, value = df2[, value_column, drop = FALSE])

            # Step 4: Merge with LD reference and filter
            variants_ld_block_match <- merge(cohort_df, ld_meta_file, by = "chrom", allow.cartesian = TRUE) %>%
              filter(pos > start & pos < end) %>%
              select(-path)

            # Function to process each group
            process_group <- function(data) {
              # Construct file path
              bim_file_path <- unique(data$bim_path)
              ld_bim_file <- vroom(bim_file_path)

              # Perform allele quality control
              flipped_data <- allele_qc(data, ld_bim_file$V2,
                col_to_flip = c(value_column),
                match_min_prop = 0, remove_dups = FALSE,
                remove_indels = FALSE, remove_strand_ambiguous = FALSE,
                flip_strand = FALSE, remove_unmatched = TRUE
              )$target_data_qced
              return(flipped_data)
            }

            final_df <- variants_ld_block_match %>%
              group_by(start, end) %>%
              group_map(~ process_group(.x)) %>%
              bind_rows() %>%
              select(c("variant_id", value_column)) %>%
              rename("variants" = "variant_id")
            # Rename columns to avoid duplication
            colnames(final_df) <- c(id_column, paste0(value_column, "_", i))
          } else if (!is.null(ref_panel)) {
            # Step 2: Split 'variants' to extract chromosomal info
            cohort_variants_df <- parse_variant_id(df2[, c(id_column)])
            # Step 3: Combine extracted chromosomal info with value column
            cohort_df <- cbind(cohort_variants_df, value = df2[, value_column, drop = FALSE])
          
            flipped_data <- allele_qc(cohort_df, ref_panel, col_to_flip = c(value_column),
                match_min_prop = 0, remove_dups = FALSE,
                remove_indels = FALSE, remove_strand_ambiguous = FALSE,
                flip_strand = FALSE, remove_unmatched = TRUE, remove_same_vars = FALSE)$target_data_qced
              
            final_df <- flipped_data %>%
                select(c("variant_id", value_column))
            colnames(final_df) <- c(id_column, paste0(value_column, "_", i))
          } else {
            final_df <- df2
            colnames(final_df) <- c(id_column, paste0(value_column, "_", i))
          }
          return(final_df)
        },
        error = function(e) {
          message(paste("Error processing dataset", i, ":", e$message))
          return(NULL)
        }
      )
    })

    # Remove any NULL results from errors
    df_list <- df_list[!sapply(df_list, is.null)]
    if (length(df_list) == 0) {
        message("No valid datasets after processing")
        return(NULL)
    }

    # Iteratively merge the data frames
    merged_df <- Reduce(
      function(x, y) merge(x, y, by = id_column, all = TRUE),
      df_list
    )
    # Optionally, remove rows with any missing values
    if (remove_any_missing) {
      merged_df <- merged_df[complete.cases(merged_df), ]
    }
    return(merged_df)
  }
                 
#' Load and Align Summary Statistics for a Given Gene and Condition
#'
#' @description
#' This function processes summary statistics matrices for a target gene across contexts, 
#' optionally aligning with a reference panel and updating an existing result list.
#'
#' @param dat_list A named list of matrices or data.frames, each element corresponding to a summary statistics type (e.g., z, beta).
#' @param signal_df A data.frame containing signal information including `variant_ID`, `gene_ID`, and `event_ID`.
#' @param cond Character. Condition type: "strong", "null", or "random".
#' @param region Character. Target gene ID.
#' @param extract_infs Character vector. Names of summary statistics to extract (e.g., `"z"`, `"beta"`).
#' @param tag_patterns Optional named pattern list used to classify context.
#' @param result_list_format A nested list used as a running result container.
#'
#' @importFrom stringr str_detect
#' @import dplyr tidyr tibble                
#' @return The updated `result_list_format` with processed results for the specified gene and condition.
#' @export
load_multicontext_sumstats <- function(dat_list, signal_df, cond, region, extract_infs = "z", tag_patterns = NULL, result_list_format) {
  # Initialize output list
  out <- list()
  trait_names <- names(dat_list[[1]])
    if (cond == "strong" && region %in% signal_df$gene_ID){
  events <- signal_df %>% filter(gene_ID == region) %>% pull(event_ID) %>% unique()                    
  for (j in seq_along(events)){                           
        ref_df_filtered <- signal_df %>% filter(gene_ID == region, event_ID == events[j]) %>% 
            filter(!str_detect(context_classify, "NE"))
        if(dim(ref_df_filtered)[1] == 0) next
        ## generate the reference panel for allele flipping
        ref_panel <- parse_variant_id(ref_df_filtered$variant_ID%>%unique())
        
        var_idx <- c()
        variants <- c()                             
        sumstats_df <- list()
        event_ID_extracted <- c()

        # Flatten the nested list
        for (extract_inf in extract_infs) {
        extracted_matrix <- merge_sumstats_matrices(dat_list[[extract_inf]], value_column = extract_inf, ref_panel = ref_panel, id_column = "variants", remove_any_missing = FALSE)
        if(is.null(extracted_matrix)||dim(extracted_matrix)[1]==0) return(result_list_format)
        out[[extract_inf]] <- extracted_matrix
        # Set variant order on first iteration
        if (is.null(var_idx)&& is.null(variants)) {
            var_idx <- 1:nrow(out[[extract_inf]])
            variants <- out[[extract_inf]]$variants[var_idx]
        }
        number_index <- str_extract(colnames(out[[extract_inf]]), "\\d+")[-1]
        out[[extract_inf]] <- out[[extract_inf]][var_idx, , drop = FALSE]
        rownames(out[[extract_inf]]) <- variants
        colnames(out[[extract_inf]])[2:ncol(out[[extract_inf]])] <- trait_names[as.integer(number_index)]
        out[[extract_inf]] <- out[[extract_inf]][, -which(names(out[[extract_inf]]) == "variants"), drop = FALSE]

        df <- as.data.frame(t(out[[extract_inf]]))
        df <- rownames_to_column(df, var = "context")

            # Match context to tag
        df <- df %>%
                  mutate(context_classify = if (is.null(tag_patterns) || length(tag_patterns) == 0) {
                    context
                  } else {
                    purrr::map_chr(context, function(ctx) {
                      matched <- names(tag_patterns)[str_detect(ctx, tag_patterns)]
                      if (length(matched) == 0) NA_character_ else matched[1]
                    })
                  })

        numeric_col <- colnames(df)[2]

         if (extract_inf == "z"){
                # Make a copy to store added rows
              added_df <- data.frame()

                # Ensure the column name of the numeric column
                if (any(grepl("sQTL|pQTL|gpQTL", df$context_classify))) {
                  if (any(grepl("sQTL|pQTL|gpQTL", ref_df_filtered$context_classify))) {
    
                    # Extract sQTL contexts to loop over
                    xQTL_specific_contexts <- unique(str_subset(ref_df_filtered$context_classify, "sQTL|pQTL|gpQTL"))
    
                    for (cont in xQTL_specific_contexts) {
                          event_IDs_extracted <- ref_df_filtered %>%
                                    filter(context_classify == cont) %>%
                                    pull(event_IDs)
      
                    # Filter matching rows in df
                          context_rows <- df %>%
                                filter(context_classify == cont, str_detect(context, paste(event_IDs_extracted, collapse = "|")))
      
                      if (nrow(context_rows) > 0) {
                        # Get the row with median absolute value
                        abs_values <- abs(context_rows[[numeric_col]])
                        median_val <- median(abs_values, na.rm = TRUE)
                        median_idx <- which.min(abs(abs_values - median_val))  # Closest to median
                        selected_df <- context_rows[median_idx, , drop = FALSE]
        
                        added_df <- bind_rows(added_df, selected_df)
                        df <- df %>% filter(context_classify !=cont)
                      }
                    }
                    # Combine updated sQTL-specific rows back into df
                    df <- bind_rows(df, added_df)
                  }
                }               
                sumstats_df[[extract_inf]] <- df %>%
                  filter(!str_detect(context_classify, "NE") & context_classify != 'NA')%>%
                  group_by(context_classify) %>%
                  slice_min(order_by = abs(.data[[numeric_col]] - median(abs(.data[[numeric_col]]), na.rm = TRUE)), n = 1, with_ties = FALSE) %>%
                  ungroup()%>%
                  rename(!!numeric_col := !!sym(numeric_col))
                event_ID_extracted <- sumstats_df[[extract_inf]]%>%pull(context)
             } else if (is.null(event_ID_extracted)){
                    warning("Please provide 'z-score'")
             } else {
                 sumstats_df[[extract_inf]] <- df %>% filter(context%in%event_ID_extracted)%>%
                                        rename(!!numeric_col := !!sym(numeric_col))
             }
             result_df <- sumstats_df[[extract_inf]] %>%
                  select(-context) %>%
                  rename(value = !!sym(numeric_col)) %>%
                  pivot_wider(names_from = context_classify, values_from = value) %>%
                  mutate(
                    variant_ID = numeric_col,
                    gene_ID = region
                  ) %>%
                 select(variant_ID, gene_ID, everything())
                 result_list_format[[cond]][[extract_inf]]  <- result_list_format[[cond]][[extract_inf]]%>% rows_update(result_df, by = c("variant_ID", "gene_ID"))
     }
  }
}
  # Handle "null" condition
  if (cond%in%c("null","random") && region %in% signal_df$gene_ID) {
    ref_df_filtered <- signal_df %>% filter(gene_ID == region)
    ref_panel <- parse_variant_id(ref_df_filtered$variant_ID %>% unique())

    var_idx <- c()
    variants <- c()                             
    sumstats_df <- list()
    event_ID_extracted <- list()
    for (extract_inf in extract_infs){
         # Flatten the nested list
         extracted_matrix <- merge_sumstats_matrices(dat_list[[extract_inf]], value_column = extract_inf, ref_panel = ref_panel, id_column = "variants", remove_any_missing = FALSE)
          if (is.null(extracted_matrix)||dim(extracted_matrix)[1]==0) return(result_list_format)
         out[[extract_inf]] <- extracted_matrix
          # Set variant order on first iteration
          if (is.null(var_idx)&& is.null(variants)) {
                var_idx <- 1:nrow(out[[extract_inf]])
                variants <- out[[extract_inf]]$variants[var_idx]
          }
          number_index <- str_extract(colnames(out[[extract_inf]]), "\\d+")[-1]
          out[[extract_inf]] <- out[[extract_inf]][var_idx, , drop = FALSE]
          rownames(out[[extract_inf]]) <- variants
          colnames(out[[extract_inf]])[2:ncol(out[[extract_inf]])] <- trait_names[as.integer(number_index)]
          out[[extract_inf]] <- out[[extract_inf]][, -which(names(out[[extract_inf]]) == "variants"), drop = FALSE]

          for (k in 1: dim(out[[extract_inf]])[1]){
               df <- as.data.frame(t(out[[extract_inf]][k,]))
               df <- rownames_to_column(df, var = "context")

              # Match context to tag
               df <- df %>%
                   mutate(context_classify = if (is.null(tag_patterns) || length(tag_patterns) == 0) {
                     context
                   } else {
                     purrr::map_chr(context, function(ctx) {
                       matched <- names(tag_patterns)[str_detect(ctx, tag_patterns)]
                       if (length(matched) == 0) NA_character_ else matched[1]
                     })
                   })

              numeric_col <- colnames(df)[2]
            if (extract_inf == "z"){
                sumstats_df[[extract_inf]] <- df %>%
                          filter(!str_detect(context_classify, "NE") & context_classify != 'NA')
                if(cond == "null"){
                  sumstats_df[[extract_inf]] <- sumstats_df[[extract_inf]] %>% 
                        group_by(context_classify) %>%
                        filter(
                            !is.na(.data[[numeric_col]]),
                            if (any(str_detect(context_classify, "sQTL|pQTL|gpQTL"))) {
                              abs(.data[[numeric_col]]) < 2
                            } else {
                              TRUE
                            }
                          )%>%
                         slice_min(
                            order_by = abs(.data[[numeric_col]] - median(abs(.data[[numeric_col]]), na.rm = TRUE)),
                            n = 1,
                            with_ties = FALSE
                          ) %>%
                          ungroup() %>%
                          rename(!!numeric_col := !!sym(numeric_col))
                } else if (cond == "random") {
                    sumstats_df[[extract_inf]] <-  sumstats_df[[extract_inf]] %>%
                          group_by(context_classify) %>%
                          slice_min(
                            order_by = abs(.data[[numeric_col]] - median(abs(.data[[numeric_col]]), na.rm = TRUE)),
                            n = 1,
                            with_ties = FALSE
                          ) %>%
                          ungroup() %>%
                          rename(!!numeric_col := !!sym(numeric_col))
                }
                event_ID_extracted[[k]] <- sumstats_df[[extract_inf]]%>%pull(context)
            }  else if (is.null(event_ID_extracted)){
                    warning("Please provide 'z-score'")
            } else {
                sumstats_df[[extract_inf]] <- df %>% filter(context%in%event_ID_extracted[[k]])%>%
                                        rename(!!numeric_col := !!sym(numeric_col))
            }
            result_df <- sumstats_df[[extract_inf]] %>%
                  select(-context) %>%
                  rename(value = !!sym(numeric_col)) %>%
                  pivot_wider(names_from = context_classify, values_from = value) %>%
                  mutate(
                    variant_ID = numeric_col,
                    gene_ID = region
                  ) %>%
                 select(variant_ID, gene_ID, everything())
            result_list_format[[cond]][[extract_inf]]  <- result_list_format[[cond]][[extract_inf]] %>% rows_update(result_df, by = c("variant_ID", "gene_ID"))
            }
          }
       }
     return(result_list_format)
   }

            
#' Extract Summary Statistics from Nested Data Structure
#'
#' @description
#' Recursively searches a nested list to extract summary statistics (z, beta, or se) 
#' using `variant_names` and `sumstats`. Computes `z` if needed from `betahat` and `sebetahat`.
#'
#' @param data A nested list structure potentially containing `variant_names` and `sumstats`.
#' @param extract_inf Character. One of `"z"`, `"beta"`, or `"se"`.
#' @param max_depth Integer. Maximum depth to search within the list. Default is 3.
#'
#' @return A data.frame with columns `variants` and the requested summary statistic.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- extract_data(nested_list_object, extract_inf = "z")
#' }

extract_flatten_sumstats_from_nested <- function(data, extract_inf = "z", max_depth = 3) {
  # Validate input
  if (!extract_inf %in% c("z", "beta", "se")) {
    stop("extract_inf must be one of: 'z', 'beta', or 'se'")
  }

  # Internal recursive function
  find_nested <- function(element, current_depth = 0) {
    if (current_depth >= max_depth) {
      message("Maximum search depth reached. Could not find 'variant_names' and 'sumstats' together.")
      return(NULL)
    }

    if (is.list(element)) {
      if (all(c("variant_names", "sumstats") %in% names(element))) {
        variant_names <- element$variant_names
        sumstats <- element$sumstats

        # Extract based on type
        result_column <- switch(
          extract_inf,
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

        result <- data.frame(variants = variant_names)
        result[[extract_inf]] <- result_column

        # Normalize variants to canonical format (with chr prefix)
        result$variants <- normalize_variant_id(result$variants)

        return(result)
      }

      # Recurse into nested elements
      for (name in names(element)) {
        result <- find_nested(element[[name]], current_depth + 1)
        if (!is.null(result)) {
          result$variants <- normalize_variant_id(result$variants)
          return(result)
        }
      }
    }

    return(NULL)
  }

  # Start search
  return(find_nested(data))
}