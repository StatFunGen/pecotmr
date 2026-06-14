#' heterogeneity:  calculate I2 statistics based on the Cochran's Q statistic
#' @noRd
calcI2 <- function(Q, Est) {
  Q <- Q[[1]]
  Est <- length(unique(Est))
  I2 <- if (Q > 1e-3) (Q - Est + 1) / Q else 0
  return(if (I2 < 0) 0 else I2)
}

# Create a null data frame with geneName and NA columns for MR pipeline outputs.
# @noRd
.createNullMrDf <- function(geneName, colSpec) {
  n <- length(geneName)
  cols <- map(colSpec, function(type) {
    switch(type,
      character = as.character(rep(NA, n)),
      integer = as.integer(rep(NA, n)),
      numeric = as.numeric(rep(NA, n))
    )
  })
  do.call(data.frame, c(list(geneName = geneName), cols, list(stringsAsFactors = FALSE)))
}

#' MR Format Function
#'
#' Description of what the function does.
#'
#' @param susieResult A list containing the results of SuSiE analysis. This list should include nested elements such as 'susie_results', 'finemappingResult' (a FineMappingResult S4 object), and 'top_loci', containing details about the statistical analysis of genetic variants.
#' @param condition A character string specifying the conditions. This is used to select the corresponding subset of results within 'susieResult'.
#' @param gwasSumstatsDb A data frame containing summary statistics from GWAS studies. It should include columns for variant id and their associated statistics such as beta coefficients and standard errors.
#' @param coverage A character string specifying the credible set column. If
#'   NULL, it is derived from \code{coverageLevel} and \code{method}.
#' @param runAlleleQc Whether to run allele QC on variants. Default TRUE.
#' @param method Fine-mapping method suffix used for method-specific columns.
#' @param coverageLevel Numeric credible set coverage used when \code{coverage}
#'   is NULL.
#' @return A data frame formatted for MR analysis or NULL if cs_list is empty.
#' @importFrom stringr str_remove str_split_i
#' @export
mrFormat <- function(susieResult, condition, gwasSumstatsDb, coverage = NULL,
                     runAlleleQc = TRUE, method = "susie", coverageLevel = 0.95,
                     molecularNameObj = c("susie_results", condition, "region_info", "region_name"), ldMetaDf) {
  mrFormatSpec <- c(variant_id = "character", bhat_x = "numeric", sbhat_x = "numeric",
                    cs = "numeric", pip = "numeric", bhat_y = "numeric", sbhat_y = "numeric")
  geneName <- unique(getNestedElement(susieResult, molecularNameObj))

  # Attempt to retrieve top_loci; if not found, return null
  topLoci <- tryCatch(
    getNestedElement(susieResult, c("susie_results", condition, "top_loci")),
    error = function(e) {
      message("top_loci does not exist for the specified condition in susieResult.")
      NULL
    }
  )
  topLoci <- .translateLegacyTopLociCsColumns(topLoci)
  if (!is.data.frame(topLoci)) return(.createNullMrDf(geneName, mrFormatSpec))
  if (is.null(coverage)) coverage <- formatCsColumn(coverageLevel, method)
  coverage <- .translateLegacyCsColumnName(coverage)
  csValues <- topLoci[[coverage]]
  if (is.null(csValues) || !any(!is.na(csValues) & csValues != 0)) {
    return(.createNullMrDf(geneName, mrFormatSpec))
  }
  pipCol <- resolvePipColumn(topLoci, method)
  if (is.null(pipCol)) return(.createNullMrDf(geneName, mrFormatSpec))

  susieCsResultFormatted <- topLoci[!is.na(csValues) & csValues >= 1, , drop = FALSE] %>%
    mutate(geneName = geneName) %>%
    mutate(
      variant_id = normalizeVariantId(variant_id),
      variant = stripChrPrefix(variant_id)
    ) %>%
    select(geneName, variant_id, variant, betahat, sebetahat, all_of(coverage), all_of(pipCol)) %>%
    rename("bhat_x" = "betahat", "sbhat_x" = "sebetahat", "cs" = all_of(coverage), "pip" = all_of(pipCol))

  susiePos <- str_split_i(susieCsResultFormatted$variant, ":", 2)
  gwasPos <- str_split_i(gwasSumstatsDb$variant_id, ":", 2)
  if (!any(susiePos %in% gwasPos)) return(.createNullMrDf(geneName, mrFormatSpec))

  gwasSumstatsDbExtracted <- gwasSumstatsDb %>%
    filter(pos %in% susiePos) %>%
    mutate(n_sample = if ("n_sample" %in% colnames(.)) n_sample else (n_case + n_control))
  meanNSample <- round(mean(gwasSumstatsDbExtracted$n_sample, na.rm = TRUE))
  # Impute `n_sample` and `maf`
  if (any(is.na(gwasSumstatsDbExtracted$effect_allele_frequency))) {
    gwasSumstatsDbExtracted <- gwasSumstatsDbExtracted %>%
      left_join(ldMetaDf %>% select(pos, allele_freq), by = "pos") %>%
      mutate(effect_allele_frequency = ifelse(is.na(effect_allele_frequency), allele_freq, effect_allele_frequency)) %>%
      mutate(n_sample = ifelse(is.na(n_sample), meanNSample, n_sample)) %>%
      select(-allele_freq)
  }
  gwasBetaSe <- zToBetaSe(gwasSumstatsDbExtracted$z, gwasSumstatsDbExtracted$effect_allele_frequency, gwasSumstatsDbExtracted$n_sample)
  gwasSumstatsDbExtracted <- gwasSumstatsDbExtracted %>% mutate(beta = gwasBetaSe$beta, se = gwasBetaSe$se)
  if (runAlleleQc) {
    susieCsResultFormatted <- matchRefPanel(cbind(variantIdToDf(susieCsResultFormatted$variant), susieCsResultFormatted),
      gwasSumstatsDbExtracted$variant_id, c("bhat_x"),
      matchMinProp = 0
    )
    susieCsResultFormatted <- getHarmonizedData(susieCsResultFormatted)[, c("geneName", "variant_id", "bhat_x", "sbhat_x", "cs", "pip")]
  }
  # Ensure consistent chr prefix convention before intersecting
  if (nrow(susieCsResultFormatted) == 0) return(.createNullMrDf(geneName, mrFormatSpec))
  if (!is.null(susieCsResultFormatted$variant_id) && !is.null(gwasSumstatsDbExtracted$variant_id)) {
    chrMatched <- ensureChrMatch(susieCsResultFormatted$variant_id, gwasSumstatsDbExtracted$variant_id)
    susieCsResultFormatted$variant_id <- chrMatched$idsA
    gwasSumstatsDbExtracted$variant_id <- chrMatched$idsB
  }
  commonVariants <- intersect(susieCsResultFormatted$variant_id, gwasSumstatsDbExtracted$variant_id)
  if (length(commonVariants) == 0) return(.createNullMrDf(geneName, mrFormatSpec))

  susieCsResultFormatted[match(commonVariants, susieCsResultFormatted$variant_id), ] %>%
    cbind(., gwasSumstatsDbExtracted[match(commonVariants, gwasSumstatsDbExtracted$variant_id), ] %>%
      select(beta, se) %>%
      rename("bhat_y" = "beta", "sbhat_y" = "se"))
}

#' Mendelian Randomization (MR)
#'
#' @param mrFormattedInput the output of twas_mr_format_input function
#' @param cpipCutoff the threshold of cumulative posterior inclusion probability, default is 0.5
#' @return A single data frame of output with columns "geneName", "num_CS", "num_IV",
#' "meta_eff", "se_meta_eff", "meta_pval", "Q", "Q_pval" and "I2". "geneName" is ensemble ID. "num_CS" is the number of credible sets
#' contained in each gene, "num_IV" is the number of variants contained in each gene. "meta_eff", "se_meta_eff" and "meta_pval" are the MR estimate, standard error and pvalue.
#' "Q" is Cochran's Q statistic, "I2" quantifies the heterogeneity, range from 0 to 1.
#' @importFrom dplyr mutate group_by filter ungroup distinct arrange select summarise n left_join
#' @importFrom magrittr %>%
#' @importFrom stats pnorm pchisq
#' @export
mrAnalysis <- function(mrFormattedInput, cpipCutoff = 0.5) {
  mrOutputSpec <- c(num_CS = "integer", num_IV = "integer", cpip = "numeric",
                    meta_eff = "numeric", se_meta_eff = "numeric", meta_pval = "numeric",
                    Q = "numeric", Q_pval = "numeric", I2 = "numeric")
  if (all(is.na(mrFormattedInput[, -1]))) {
    return(.createNullMrDf(unique(mrFormattedInput$geneName), mrOutputSpec))
  }
  output <- mrFormattedInput %>%
    mutate(
      bhat_x = bhat_x / sbhat_x,
      sbhat_x = 1
    ) %>%
    group_by(geneName, cs) %>%
    mutate(cpip = sum(pip)) %>%
    filter(cpip >= cpipCutoff)

  if (nrow(output) == 0) {
    return(.createNullMrDf(unique(mrFormattedInput$geneName), mrOutputSpec))
  }

  # Compute per-CS composite estimates
  csSummary <- output %>%
    group_by(geneName, cs) %>%
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
  geneSummary <- csSummary %>%
    group_by(geneName) %>%
    summarise(
      num_CS = n(),
      cpip = first(cpip),
      sum_w = sum(wv),
      meta_eff = sum(wv * composite_bhat) / sum(wv),
      se_meta_eff = sqrt(1 / sum(wv)),
      Q = sum(wv * (composite_bhat - sum(wv * composite_bhat) / sum(wv))^2),
      I2 = calcI2(Q, composite_bhat),
      Q_pval = pchisq(Q, df = n() - 1, lower = FALSE),
      .groups = "drop"
    ) %>%
    mutate(
      meta_pval = 2 * pnorm(abs(meta_eff) / se_meta_eff, lower.tail = FALSE)
    )

  # Add num_IV from original output
  ivCounts <- output %>%
    group_by(geneName) %>%
    summarise(num_IV = n(), .groups = "drop")

  geneSummary %>%
    left_join(ivCounts, by = "geneName") %>%
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
    select(geneName, num_CS, num_IV, cpip, meta_eff, se_meta_eff, meta_pval, Q, Q_pval, I2)
}

#' Fine-mapping-based Mendelian Randomization
#'
#' Performs MR using fine-mapping credible sets. For each exposure gene,
#' computes PIP-weighted composite Wald ratio estimates within each credible
#' set, then meta-analyzes across credible sets with inverse-variance weighting.
#' Reports Cochran's Q and I-squared heterogeneity statistics.
#'
#' @param formattedInput A data.frame/tibble with columns: X_ID (gene ID),
#'   cs (credible set index), pip (posterior inclusion probability),
#'   bhat_x (exposure effect), sbhat_x (exposure SE), bhat_y (outcome effect),
#'   sbhat_y (outcome SE), snp (variant ID).
#' @param cpipCutoff Minimum cumulative PIP to retain a credible set
#'   (default 0.5).
#' @return A tibble with columns: X_ID, num_CS, num_IV, cpip,
#'   composite_bhat, composite_sbhat, meta_eff, se_meta_eff, Q, I2.
#' @export
fineMr <- function(formattedInput, cpipCutoff = 0.5) {
  resultCols <- c("X_ID", "num_CS", "num_IV", "cpip", "composite_bhat",
                  "composite_sbhat", "meta_eff", "se_meta_eff", "Q", "I2")

  filtered <- formattedInput %>%
    mutate(
      bhat_x = bhat_x / sbhat_x,
      sbhat_x = 1) %>%
    group_by(X_ID, cs) %>%
    mutate(cpip = sum(pip)) %>%
    filter(cpip >= cpipCutoff)

  if (nrow(filtered) == 0) {
    return(tibble(!!!setNames(
      rep(list(logical(0)), length(resultCols)), resultCols)))
  }

  filtered %>%
    group_by(X_ID, cs) %>%
    mutate(
      beta_yx = bhat_y / bhat_x,
      se_yx = sqrt(
        (sbhat_y^2 / bhat_x^2) + ((bhat_y^2 * sbhat_x^2) / bhat_x^4)),
      composite_bhat = sum((beta_yx * pip) / cpip),
      composite_sbhat = sum((beta_yx^2 + se_yx^2) * pip / cpip)) %>%
    mutate(
      composite_sbhat = sqrt(composite_sbhat - composite_bhat^2),
      wv = composite_sbhat^-2) %>%
    ungroup() %>%
    group_by(X_ID) %>%
    mutate(
      meta_eff = sum(unique(wv) * unique(composite_bhat)),
      sum_w = sum(unique(wv)),
      se_meta_eff = sqrt(sum_w^-1),
      num_CS = length(unique(cs))) %>%
    mutate(
      num_IV = length(snp),
      meta_eff = meta_eff / sum_w,
      Q = sum(unique(wv) * (unique(composite_bhat) - unique(meta_eff))^2),
      I2 = calcI2(Q, composite_bhat)) %>%
    ungroup() %>%
    distinct(X_ID, .keep_all = TRUE) %>%
    mutate(
      cpip = round(cpip, 3),
      composite_bhat = round(composite_bhat, 3),
      meta_eff = round(meta_eff, 3),
      se_meta_eff = round(se_meta_eff, 3),
      Q = round(Q, 3),
      I2 = round(I2, 3)) %>%
    select(X_ID, num_CS, num_IV, cpip, composite_bhat, composite_sbhat,
                  meta_eff, se_meta_eff, Q, I2)
}

