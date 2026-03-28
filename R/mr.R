#' heterogeneity:  calculate I2 statistics based on the Cochran's Q statistic
#' @noRd
calc_I2 <- function(Q, Est) {
  Q <- Q[[1]]
  Est <- length(unique(Est))
  I2 <- if (Q > 1e-3) (Q - Est + 1) / Q else 0
  return(if (I2 < 0) 0 else I2)
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
#' @param allele_qc Optional. A logical value indicating whether allele qc should be performed on the variants. When TRUE, allele qc are applied to the variants based on the GWAS summary statistics database ('gwas_sumstats_db').
#' @return A data frame formatted for MR analysis or NULL if cs_list is empty.
#' @importFrom stringr str_remove
#' @export
mr_format <- function(susie_result, condition, gwas_sumstats_db, coverage = "cs_coverage_0.95", allele_qc = TRUE,
                      molecular_name_obj = c("susie_results", condition, "region_info", "region_name"), ld_meta_df) {
  # Create null mr_format_input
  create_null_mr_input <- function(gene_name) {
    mr_format_input <- data.frame(
      gene_name = gene_name,
      variant_id = as.character(rep(NA, length(gene_name))),
      bhat_x = as.numeric(rep(NA, length(gene_name))),
      sbhat_x = as.numeric(rep(NA, length(gene_name))),
      cs = as.numeric(rep(NA, length(gene_name))),
      pip = as.numeric(rep(NA, length(gene_name))),
      bhat_y = as.numeric(rep(NA, length(gene_name))),
      sbhat_y = as.numeric(rep(NA, length(gene_name))),
      stringsAsFactors = FALSE # Optional, to prevent factors
    )
  }
  gene_name <- unique(get_nested_element(susie_result, molecular_name_obj))

  # Attempt to retrieve top_loci; if not found, return NULL
  top_loci <- tryCatch(
    {
      get_nested_element(susie_result, c("susie_results", condition, "top_loci"))
    },
    error = function(e) {
      message("top_loci does not exist for the specified condition in susie_result.")
      return(NULL)
    }
  )
  if (is.data.frame(top_loci)) {
    if (any(unique(get_nested_element(top_loci, coverage)) != 0)) {
      susie_cs_result_formatted <- top_loci %>%
        mutate(gene_name = gene_name) %>%
        filter(coverage >= 1) %>%
        mutate(variant = strip_chr_prefix(variant_id)) %>%
        select(gene_name, variant, betahat, sebetahat, all_of(coverage), pip) %>%
        rename("bhat_x" = "betahat", "sbhat_x" = "sebetahat", "cs" = all_of(coverage))
      susie_pos <- sapply(susie_cs_result_formatted$variant, function(variant_id) strsplit(variant_id, "\\:")[[1]][2])
      gwas_pos <- sapply(gwas_sumstats_db$variant_id, function(variant_id) strsplit(variant_id, "\\:")[[1]][2])
      if (any(susie_pos %in% gwas_pos)) {
        gwas_sumstats_db_extracted <- gwas_sumstats_db %>%
          filter(pos %in% susie_pos) %>%
          mutate(n_sample = if ("n_sample" %in% colnames(.)) n_sample else (n_case + n_control))
        mean_n_sample <- round(mean(gwas_sumstats_db_extracted$n_sample, na.rm = TRUE))
        # Impute `n_sample` and `maf`
        if (any(is.na(gwas_sumstats_db_extracted$effect_allele_frequency))) {
          gwas_sumstats_db_extracted_imputed <- gwas_sumstats_db_extracted %>%
            left_join(ld_meta_df %>% select(pos, allele_freq), by = "pos") %>%
            mutate(effect_allele_frequency = ifelse(is.na(effect_allele_frequency), allele_freq, effect_allele_frequency)) %>%
            mutate(n_sample = ifelse(is.na(n_sample), mean_n_sample, n_sample)) %>%
            select(-allele_freq)
        } else {
          gwas_sumstats_db_extracted_imputed <- gwas_sumstats_db_extracted
        }
        gwas_sumstats_db_beta_se <- z_to_beta_se(gwas_sumstats_db_extracted_imputed$z, gwas_sumstats_db_extracted_imputed$effect_allele_frequency, gwas_sumstats_db_extracted_imputed$n_sample)
        gwas_sumstats_db_extracted_imputed <- gwas_sumstats_db_extracted_imputed %>% mutate(beta = gwas_sumstats_db_beta_se$beta, se = gwas_sumstats_db_beta_se$se)
        if (allele_qc) {
          susie_cs_result_formatted <- allele_qc(cbind(variant_id_to_df(susie_cs_result_formatted$variant), susie_cs_result_formatted),
            gwas_sumstats_db_extracted_imputed$variant_id, c("bhat_x"),
            match_min_prop = 0
          )
          susie_cs_result_formatted <- susie_cs_result_formatted$target_data_qced[, c("gene_name", "variant_id", "bhat_x", "sbhat_x", "cs", "pip")]
        }
        # Normalize variant IDs to canonical format for matching
        gwas_sumstats_db_extracted_imputed$variant_id <- normalize_variant_id(gwas_sumstats_db_extracted_imputed$variant_id)
        susie_cs_gwas_variants_merge <- intersect(susie_cs_result_formatted$variant_id, gwas_sumstats_db_extracted_imputed$variant_id)

        mr_format_input <- susie_cs_result_formatted[match(susie_cs_gwas_variants_merge, susie_cs_result_formatted$variant), ] %>%
          cbind(., gwas_sumstats_db_extracted_imputed[match(susie_cs_gwas_variants_merge, gwas_sumstats_db_extracted_imputed$variant_id), ] %>%
            select(beta, se) %>%
            rename("bhat_y" = "beta", "sbhat_y" = "se"))
      } else {
        mr_format_input <- create_null_mr_input(gene_name)
      }
    } else {
      mr_format_input <- create_null_mr_input(gene_name)
    }
  } else {
    mr_format_input <- create_null_mr_input(gene_name)
  }
  return(mr_format_input)
}

#' Mendelian Randomization (MR)
#'
#' @param mr_formatted_input the output of twas_mr_format_input function
#' @param cpip_cutoff the threshold of cumulative posterior inclusion probability, default is 0.5
#' @return A single data frame of output with columns "gene_name", "num_CS", "num_IV",
#' "meta_eff", "se_meta_eff", "meta_pval", "Q", "Q_pval" and "I2". "gene_name" is ensemble ID. "num_CS" is the number of credible sets
#' contained in each gene, "num_IV" is the number of variants contained in each gene. "meta_eff", "se_meta_eff" and "meta_pval" are the MR estimate, standard error and pvalue.
#' "Q" is Cochran's Q statistic, "I2" quantifies the heterogeneity, range from 0 to 1.
#' @importFrom dplyr mutate group_by filter ungroup distinct arrange select
#' @importFrom magrittr %>%
#' @importFrom stats pnorm pchisq
#' @export
mr_analysis <- function(mr_formatted_input, cpip_cutoff = 0.5) {
  create_null_output <- function(gene_name) {
    data.frame(
      gene_name = gene_name,
      num_CS = as.integer(rep(NA, length(gene_name))),
      num_IV = as.integer(rep(NA, length(gene_name))),
      cpip = as.numeric(rep(NA, length(gene_name))),
      meta_eff = as.numeric(rep(NA, length(gene_name))),
      se_meta_eff = as.numeric(rep(NA, length(gene_name))),
      meta_pval = as.numeric(rep(NA, length(gene_name))),
      Q = as.numeric(rep(NA, length(gene_name))),
      Q_pval = as.numeric(rep(NA, length(gene_name))),
      I2 = as.numeric(rep(NA, length(gene_name))),
      stringsAsFactors = FALSE
    )
  }
  if (all(is.na(mr_formatted_input[, -1]))) {
    return(create_null_output(unique(mr_formatted_input$gene_name)))
  }
  output <- mr_formatted_input %>%
    mutate(
      bhat_x = bhat_x / sbhat_x,
      sbhat_x = 1
    ) %>%
    group_by(gene_name, cs) %>%
    mutate(cpip = sum(pip)) %>%
    filter(cpip >= cpip_cutoff) # Cumulative pip greater than a user defined cumulative pip threshold

  if (dim(output)[1] != 0) {
    output <- output %>%
      group_by(gene_name, cs) %>%
      mutate(
        beta_yx = bhat_y / bhat_x,
        se_yx = sqrt((sbhat_y^2 / bhat_x^2) + ((bhat_y^2 * sbhat_x^2) / bhat_x^4)),
        composite_bhat = sum((beta_yx * pip) / cpip),
        composite_sbhat = sum((beta_yx^2 + se_yx^2) * pip / cpip)
      ) %>%
      mutate(
        composite_sbhat = sqrt(composite_sbhat - composite_bhat^2),
        wv = composite_sbhat^-2
      ) %>%
      ungroup() %>%
      mutate(
        meta_eff = sum(unique(wv) * unique(composite_bhat)),
        sum_w = sum(unique(wv)),
        se_meta_eff = sqrt(sum_w^-1),
        num_CS = length(unique(cs))
      ) %>%
      mutate(
        num_IV = length(variant_id),
        meta_eff = meta_eff / sum_w,
        meta_pval = 2 * pnorm(abs(meta_eff) / se_meta_eff, lower.tail = FALSE),
        Q = sum(unique(wv) * (unique(composite_bhat) - unique(meta_eff))^2),
        I2 = calc_I2(Q, composite_bhat),
        Q_pval = pchisq(Q, df = length(unique(composite_bhat)) - 1, lower = FALSE)
      ) %>%
      ungroup() %>%
      distinct(gene_name, .keep_all = TRUE) %>%
      mutate(
        cpip = round(cpip, 3),
        meta_pval = round(meta_pval, 3),
        meta_eff = round(meta_eff, 3),
        se_meta_eff = round(se_meta_eff, 3),
        Q = round(Q, 3),
        Q_pval = round(Q_pval, 3),
        I2 = round(I2, 3)
      ) %>%
      arrange(meta_pval) %>%
      select(gene_name, num_CS, num_IV, cpip, meta_eff, se_meta_eff, meta_pval, Q, Q_pval, I2)
  } else {
    return(create_null_output(unique(mr_formatted_input$gene_name)))
  }
  output
}
