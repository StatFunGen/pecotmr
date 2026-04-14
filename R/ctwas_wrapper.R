### File-I/O functions (ctwas_bimfile_loader, get_ctwas_meta_data) have been
### removed. Use ld_loader() and read_bim() from the standard I/O path instead.

#' Function to select variants for ctwas weights input
#' @param region_data A list of list containing weights list and snp_info list data for multiple genes/events within a single LD block region.
#' @param export_twas_weight_db A list of list of fine-mapping result data formatted by generate_twas_db function.
#' @param region_block A string for region information for region_weights, consisted of chromosome number, star and end position of LD block conneced with "_".
#' @export
trim_ctwas_variants <- function(region_data, twas_weight_cutoff = 1e-5, cs_min_cor = 0.8,
                                min_pip_cutoff = 0.1, max_num_variants = 1000) {
  # internal functions to select variants for a gene-context pair weight list
  select_variants <- function(group_name, study_name, region_data, cs_min_cor, min_pip_cutoff, max_num_variants) {
    weight_list <- region_data$weights[[group_name]][[study_name]]
    context <- weight_list$context
    selected_variants_by_context <- c()
    molecular_id <- gsub("\\|.*", "", group_name)

    if ("cs_variants" %in% names(region_data$susie_weights_intermediate[[molecular_id]][[context]]) & length(region_data$susie_weights_intermediate[[molecular_id]][[context]][["cs_variants"]]) != 0) {
      cs_min_abs_cor <- region_data$susie_weights_intermediate[[molecular_id]][[context]]$cs_purity$min.abs.corr
      for (L in seq_along(region_data$susie_weights_intermediate[[molecular_id]][[context]]$cs_variants)) {
        # we includ all variants in $cs_variant if min_abs_corr > cs_min_cor for the set
        if (cs_min_abs_cor[L] >= cs_min_cor) {
          cs_variants <- region_data$susie_weights_intermediate[[molecular_id]][[context]]$cs_variants[[L]]
          selected_variants_by_context <- cs_variants[cs_variants %in% rownames(weight_list$wgt)]
        }
      }
    }
    context_pip <- region_data$susie_weights_intermediate[[molecular_id]][[context]]$pip
    # variant IDs are in canonical chr-prefix format from allele_qc pipeline
    high_pip_variants <- names(context_pip[context_pip > min_pip_cutoff])[names(context_pip[context_pip > min_pip_cutoff]) %in% rownames(weight_list$wgt)]
    selected_variants_by_context <- unique(c(selected_variants_by_context, high_pip_variants))

    # prioritize SNPs based on PIP if max_num_variants different from Inf 
    available_variants <- intersect(rownames(weight_list$wgt), names(context_pip))
    prioritized <- unique(c(selected_variants_by_context, setdiff(available_variants, selected_variants_by_context)))
    prioritized <- prioritized[order(-context_pip[prioritized])]
    selected_variants_by_context <- head(prioritized, max_num_variants) 
    weight_list$wgt <- weight_list$wgt[selected_variants_by_context, , drop = FALSE]
    return(weight_list)
  }
  merge_by_study <- function(weights) {
    weight_list <- list()
    for (group in names(weights)) {
      for (study in names(weights[[group]])) {
        weight_list[[study]][[group]] <- weights[[group]][[study]]
      }
    }
    return(weight_list)
  }

  weights <- setNames(lapply(names(region_data$weights), function(group) {
    for (study in names(region_data$weights[[group]])) {
      region_data$weights[[group]][[study]]$wgt <- region_data$weights[[group]][[study]]$wgt[abs(region_data$weights[[group]][[study]]$wgt[, 1]) >= twas_weight_cutoff, , drop = FALSE]
      if (nrow(region_data$weights[[group]][[study]]$wgt) < 1) {
        region_data$weights[[group]][[study]] <- NULL
        next
      }
      if (all(is.na(region_data$weights[[group]][[study]]$wgt[, 1])) || all(is.nan(region_data$weights[[group]][[study]]$wgt[, 1]))) {
        region_data$weights[[group]][[study]] <- NULL
        next
      }
      if (nrow(region_data$weights[[group]][[study]]$wgt) < max_num_variants) {
        region_data$weights[[group]][[study]]$n_wgt <- nrow(region_data$weights[[group]][[study]]$wgt)
      } else {
        region_data$weights[[group]][[study]] <- select_variants(group, study, region_data, cs_min_cor = cs_min_cor, min_pip_cutoff = min_pip_cutoff, max_num_variants = max_num_variants)
        region_data$weights[[group]][[study]]$n_wgt <- nrow(region_data$weights[[group]][[study]]$wgt)
      }
      region_data$weights[[group]] <- Filter(Negate(is.null), region_data$weights[[group]])
      context_range <- as.integer(sapply(rownames(region_data$weights[[group]][[study]]$wgt), function(variant) strsplit(variant, "\\:")[[1]][2]))
      if(twas_weight_cutoff!=0 | cs_min_cor!=0 | min_pip_cutoff!=0 | max_num_variants!=Inf){
        region_data$weights[[group]][[study]][["p0"]] = min(context_range)# update min max position
        region_data$weights[[group]][[study]][["p1"]] = max(context_range)
      }
    }
    return(region_data$weights[[group]])
  }), names(region_data$weights))
  weights <- Filter(Negate(is.null), weights)
  weights <- merge_by_study(weights)
  return(weights)
}
