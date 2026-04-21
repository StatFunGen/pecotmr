#' heterogeneity:  calculate I2 statistics based on the Cochran's Q statistic
#' @noRd
calc_I2 <- function(Q, Est) {
  Q <- Q[[1]]
  Est <- length(unique(Est))
  I2 <- if (Q > 1e-3) (Q - Est + 1) / Q else 0
  return(if (I2 < 0) 0 else I2)
}

# Create a null data frame with gene_name and NA columns for MR pipeline outputs.
# @noRd
.create_null_mr_df <- function(gene_name, col_spec) {
  n <- length(gene_name)
  cols <- purrr::map(col_spec, function(type) {
    switch(type,
      character = as.character(rep(NA, n)),
      integer = as.integer(rep(NA, n)),
      numeric = as.numeric(rep(NA, n))
    )
  })
  do.call(data.frame, c(list(gene_name = gene_name), cols, list(stringsAsFactors = FALSE)))
}

#' MR Format Function
#'
#' Description of what the function does.
#'
#' @param susie_result A list containing the results of SuSiE analysis. This list should include nested elements such as 'susie_results', 'susie_result_trimmed', and 'top_loci', containing details about the statistical analysis of genetic variants.
#' @param condition A character string specifying the conditions. This is used to select the corresponding subset of results within 'susie_result'.
#' @param gwas_sumstats_db A data frame containing summary statistics from GWAS studies. It should include columns for variant id and their associated statistics such as beta coefficients and standard errors.
#' @param sets A character string indicating the method used to define sets of genetic variants. Defaults to "sets". This parameter is used to specify the type of sets to extract from the 'susie_result' object.
#' @param coverage A character string specifying the coverage threshold for credible sets, used when 'sets' is not "sets". Defaults to "coverage_0.95", indicating a 95% coverage credible set.
#' @param run_allele_qc Whether to run allele QC on variants. Default TRUE.
#' @return A data frame formatted for MR analysis or NULL if cs_list is empty.
#' @importFrom stringr str_remove
#' @export
mr_format <- function(susie_result, condition, gwas_sumstats_db, coverage = "cs_coverage_0.95", run_allele_qc = TRUE,
                      molecular_name_obj = c("susie_results", condition, "region_info", "region_name"), ld_meta_df) {
  mr_format_spec <- c(variant_id = "character", bhat_x = "numeric", sbhat_x = "numeric",
                      cs = "numeric", pip = "numeric", bhat_y = "numeric", sbhat_y = "numeric")
  gene_name <- unique(get_nested_element(susie_result, molecular_name_obj))

  # Attempt to retrieve top_loci; if not found, return null
  top_loci <- tryCatch(
    get_nested_element(susie_result, c("susie_results", condition, "top_loci")),
    error = function(e) {
      message("top_loci does not exist for the specified condition in susie_result.")
      NULL
    }
  )
  if (!is.data.frame(top_loci)) return(.create_null_mr_df(gene_name, mr_format_spec))
  if (!any(unique(get_nested_element(top_loci, coverage)) != 0)) return(.create_null_mr_df(gene_name, mr_format_spec))

  susie_cs_result_formatted <- top_loci %>%
    mutate(gene_name = gene_name) %>%
    filter(coverage >= 1) %>%
    mutate(variant = strip_chr_prefix(variant_id)) %>%
    select(gene_name, variant, betahat, sebetahat, all_of(coverage), pip) %>%
    rename("bhat_x" = "betahat", "sbhat_x" = "sebetahat", "cs" = all_of(coverage))

  susie_pos <- stringr::str_split_i(susie_cs_result_formatted$variant, ":", 2)
  gwas_pos <- stringr::str_split_i(gwas_sumstats_db$variant_id, ":", 2)
  if (!any(susie_pos %in% gwas_pos)) return(.create_null_mr_df(gene_name, mr_format_spec))

  gwas_sumstats_db_extracted <- gwas_sumstats_db %>%
    filter(pos %in% susie_pos) %>%
    mutate(n_sample = if ("n_sample" %in% colnames(.)) n_sample else (n_case + n_control))
  mean_n_sample <- round(mean(gwas_sumstats_db_extracted$n_sample, na.rm = TRUE))
  # Impute `n_sample` and `maf`
  if (any(is.na(gwas_sumstats_db_extracted$effect_allele_frequency))) {
    gwas_sumstats_db_extracted <- gwas_sumstats_db_extracted %>%
      left_join(ld_meta_df %>% select(pos, allele_freq), by = "pos") %>%
      mutate(effect_allele_frequency = ifelse(is.na(effect_allele_frequency), allele_freq, effect_allele_frequency)) %>%
      mutate(n_sample = ifelse(is.na(n_sample), mean_n_sample, n_sample)) %>%
      select(-allele_freq)
  }
  gwas_beta_se <- z_to_beta_se(gwas_sumstats_db_extracted$z, gwas_sumstats_db_extracted$effect_allele_frequency, gwas_sumstats_db_extracted$n_sample)
  gwas_sumstats_db_extracted <- gwas_sumstats_db_extracted %>% mutate(beta = gwas_beta_se$beta, se = gwas_beta_se$se)
  if (run_allele_qc) {
    susie_cs_result_formatted <- match_ref_panel(cbind(variant_id_to_df(susie_cs_result_formatted$variant), susie_cs_result_formatted),
      gwas_sumstats_db_extracted$variant_id, c("bhat_x"),
      match_min_prop = 0
    )
    susie_cs_result_formatted <- susie_cs_result_formatted$target_data_qced[, c("gene_name", "variant_id", "bhat_x", "sbhat_x", "cs", "pip")]
  }
  # Normalize variant IDs to canonical format for matching
  gwas_sumstats_db_extracted$variant_id <- normalize_variant_id(gwas_sumstats_db_extracted$variant_id)
  common_variants <- intersect(susie_cs_result_formatted$variant_id, gwas_sumstats_db_extracted$variant_id)

  susie_cs_result_formatted[match(common_variants, susie_cs_result_formatted$variant), ] %>%
    cbind(., gwas_sumstats_db_extracted[match(common_variants, gwas_sumstats_db_extracted$variant_id), ] %>%
      select(beta, se) %>%
      rename("bhat_y" = "beta", "sbhat_y" = "se"))
}

#' Mendelian Randomization (MR)
#'
#' @param mr_formatted_input the output of twas_mr_format_input function
#' @param cpip_cutoff the threshold of cumulative posterior inclusion probability, default is 0.5
#' @return A single data frame of output with columns "gene_name", "num_CS", "num_IV",
#' "meta_eff", "se_meta_eff", "meta_pval", "Q", "Q_pval" and "I2". "gene_name" is ensemble ID. "num_CS" is the number of credible sets
#' contained in each gene, "num_IV" is the number of variants contained in each gene. "meta_eff", "se_meta_eff" and "meta_pval" are the MR estimate, standard error and pvalue.
#' "Q" is Cochran's Q statistic, "I2" quantifies the heterogeneity, range from 0 to 1.
#' @importFrom dplyr mutate group_by filter ungroup distinct arrange select summarise n left_join
#' @importFrom magrittr %>%
#' @importFrom stats pnorm pchisq
#' @export
mr_analysis <- function(mr_formatted_input, cpip_cutoff = 0.5) {
  mr_output_spec <- c(num_CS = "integer", num_IV = "integer", cpip = "numeric",
                      meta_eff = "numeric", se_meta_eff = "numeric", meta_pval = "numeric",
                      Q = "numeric", Q_pval = "numeric", I2 = "numeric")
  if (all(is.na(mr_formatted_input[, -1]))) {
    return(.create_null_mr_df(unique(mr_formatted_input$gene_name), mr_output_spec))
  }
  output <- mr_formatted_input %>%
    mutate(
      bhat_x = bhat_x / sbhat_x,
      sbhat_x = 1
    ) %>%
    group_by(gene_name, cs) %>%
    mutate(cpip = sum(pip)) %>%
    filter(cpip >= cpip_cutoff)

  if (nrow(output) == 0) {
    return(.create_null_mr_df(unique(mr_formatted_input$gene_name), mr_output_spec))
  }

  # Compute per-CS composite estimates
  cs_summary <- output %>%
    group_by(gene_name, cs) %>%
    summarise(
      cpip = first(cpip),
      composite_bhat = sum((bhat_y / bhat_x * pip) / cpip),
      composite_sbhat = sqrt(
        sum(((bhat_y / bhat_x)^2 + (sbhat_y^2 / bhat_x^2) + ((bhat_y^2 * sbhat_x^2) / bhat_x^4)) * pip / cpip) -
          sum((bhat_y / bhat_x * pip) / cpip)^2
      ),
      .groups = "drop"
    ) %>%
    mutate(wv = composite_sbhat^-2)

  # Compute gene-level meta-analysis
  gene_summary <- cs_summary %>%
    group_by(gene_name) %>%
    summarise(
      num_CS = n(),
      cpip = first(cpip),
      sum_w = sum(wv),
      meta_eff = sum(wv * composite_bhat) / sum(wv),
      se_meta_eff = sqrt(1 / sum(wv)),
      Q = sum(wv * (composite_bhat - sum(wv * composite_bhat) / sum(wv))^2),
      I2 = calc_I2(Q, composite_bhat),
      Q_pval = pchisq(Q, df = n() - 1, lower = FALSE),
      .groups = "drop"
    ) %>%
    mutate(
      meta_pval = 2 * pnorm(abs(meta_eff) / se_meta_eff, lower.tail = FALSE)
    )

  # Add num_IV from original output
  iv_counts <- output %>%
    group_by(gene_name) %>%
    summarise(num_IV = n(), .groups = "drop")

  gene_summary %>%
    left_join(iv_counts, by = "gene_name") %>%
    mutate(
      cpip = round(cpip, 3),
      meta_eff = round(meta_eff, 3),
      se_meta_eff = round(se_meta_eff, 3),
      meta_pval = round(meta_pval, 3),
      Q = round(Q, 3),
      Q_pval = round(Q_pval, 3),
      I2 = round(I2, 3)
    ) %>%
    arrange(meta_pval) %>%
    select(gene_name, num_CS, num_IV, cpip, meta_eff, se_meta_eff, meta_pval, Q, Q_pval, I2)
}
