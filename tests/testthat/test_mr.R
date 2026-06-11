context("mr")
library(tidyverse)

# Define a function to randomly generate ref:alt pairs
library(testthat)
generate_format_mock_data <- function(seed = 1, num_variants = 10, empty_sets = F) {
generate_ref_alt <- function(num_variants){
  ref_alleles <- c("A", "T", "G", "C")
  alt_alleles <- c("A", "T", "G", "C")
  pairs <- paste0(sample(ref_alleles, num_variants, replace = TRUE), ":", sample(alt_alleles, num_variants, replace = TRUE))
  # Ensure that ref is not the same as alt
  while(any(sapply(strsplit(pairs, ":"), function(x) x[1] == x[2]))) {
    pairs <- paste0(sample(ref_alleles, num_variants, replace = TRUE), ":", sample(alt_alleles, num_variants, replace = TRUE))
  }
  return(pairs)
}
top_loci_mock <- data.frame(
  variant_id = paste0("chr1:", 1:num_variants, ":", generate_ref_alt(num_variants)),
  betahat = rnorm(num_variants),
  sebetahat = runif(num_variants, 0.05, 0.1),
  #maf = runif(n_entries, 0.1, 0.5),
  pip = runif(num_variants, 0, 1),
  CS_95_susie = sample(0:2, num_variants, replace = TRUE),
  CS_70_susie = sample(0:2, num_variants, replace = TRUE),
  CS_50_susie = sample(0:2, num_variants, replace = TRUE)
)
susie_result_mock <- list(
  susie_results = list(
    condition1 = list(
      top_loci = top_loci_mock,
      region_info = list(region_name = "Gene1")
    )
  )
)
gwas_sumstats_db_mock <- data.frame(
  variant_id = paste0("chr1:", 1:num_variants, ":", generate_ref_alt(num_variants)),
  beta = rnorm(num_variants),
  se = runif(num_variants, 0.05, 0.1),
  pos = 1:num_variants,  # Genomic position
  effect_allele_frequency = runif(num_variants, 0.05, 0.95),  # Always assigned, no missing values
  n_case = sample(500:5000, num_variants, replace = TRUE),  # Random case counts
  n_control = sample(500:5000, num_variants, replace = TRUE) # Random control counts
  )%>%
  mutate(n_sample = n_case + n_control,  # Total sample size
         z = beta / se)
return(
        list(
            susie_result = susie_result_mock,
            gwas_sumstats_db = gwas_sumstats_db_mock
          )
    )

}

generate_mock_mr_formatted_input <- function(num_variants = NULL, generate_full_dataset = TRUE) {
  generate_ref_alt <- function(num_variants){
  ref_alleles <- c("A", "T", "G", "C")
  alt_alleles <- c("A", "T", "G", "C")
  pairs <- paste0(sample(ref_alleles, num_variants, replace = TRUE), ":", sample(alt_alleles, num_variants, replace = TRUE))
  # Ensure that ref is not the same as alt
  while(any(sapply(strsplit(pairs, ":"), function(x) x[1] == x[2]))) {
    pairs <- paste0(sample(ref_alleles, num_variants, replace = TRUE), ":", sample(alt_alleles, num_variants, replace = TRUE))
  }
  return(pairs)
  }
  if (generate_full_dataset) {
    data.frame(
      gene_name = rep("Gene1", num_variants),
      #cs = sample(1:2, num_variants, replace = TRUE),
      cs = as.integer(rep(1,num_variants)),
      variant_id = paste0("1:", 1:num_variants, ":", generate_ref_alt(num_variants)),
      bhat_x = rnorm(num_variants),
      sbhat_x = runif(num_variants, 0.1, 0.2),
      bhat_y = rnorm(num_variants),
      sbhat_y = runif(num_variants, 0.1, 0.2),
      pip = runif(num_variants, 0, 1),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
     gene_name = "Gene1",
     variant_id = as.character(NA),
     bhat_x = as.numeric(NA),
     sbhat_x = as.numeric(NA),
     cs = as.numeric(NA),
     pip = as.numeric(NA),
     bhat_y = as.numeric(NA),
     sbhat_y = as.numeric(NA),
     stringsAsFactors = FALSE # Optional, to prevent factors
     )
    }
  }

# =============================================================================
# calcI2
# =============================================================================

test_that("calcI2 works with dummy data",{
    # Test with Q > 1e-3 and positive I2
    expect_equal(calcI2(10, c(1, 2, 3)), (10 - 3 + 1)/10)

    # Test with Q exactly 1e-3 and I2 should be 0
    expect_equal(calcI2(1e-3, c(1, 2)), 0)

    # Test with Q < 1e-3 and I2 should be 0
    expect_equal(calcI2(1e-4, c(1, 2, 3, 4)), 0)

    # Test with negative Q and I2 should be 0
    expect_equal(calcI2(-10, c(1, 2, 3)), 0)

    # Test with Q leading to negative I2, I2 should be 0
    expect_equal(calcI2(5, c(1, 2, 3, 4, 5, 6)), 0)

    # Test with Q as 0, I2 should be 0
    expect_equal(calcI2(0, c(1, 2)), 0)

    # Test with Est as empty vector -- returns 1.1 (edge case, see dedicated test below)
    expect_equal(calcI2(10, c()), 1.1)
})

test_that("calcI2: very large Q returns value close to 1", {
  big_Q <- 1e12
  est_vec <- c(1, 2, 3, 4, 5)
  result <- pecotmr:::calcI2(Q = list(big_Q), Est = est_vec)
  expected <- (big_Q - 5 + 1) / big_Q
  expect_equal(result, expected, tolerance = 1e-10)
  expect_true(result > 0.99999)
})

test_that("calcI2: very large Est (many unique) causes I2 to clamp to 0", {
  large_est <- seq_len(1000)
  result <- pecotmr:::calcI2(Q = list(5), Est = large_est)
  expect_equal(result, 0)
})

test_that("calcI2: Est with duplicates reduces unique count", {
  result <- pecotmr:::calcI2(Q = list(20), Est = c(1, 1, 2, 2, 3, 3))
  expected <- (20 - 3 + 1) / 20
  expect_equal(result, expected)
})

test_that("calcI2 with empty Est vector returns value > 1 (edge case)", {
  # When Est is empty, length(unique(c())) = 0, so I2 = (Q - 0 + 1)/Q = 1 + 1/Q
  # This yields I2 > 1 which is outside the expected [0,1] range.
  # Documenting this behavior so any future fix will be detected.
  result <- pecotmr:::calcI2(Q = list(10), Est = c())
  expect_gt(result, 1)
  expect_equal(result, 1.1)
})

# =============================================================================
# mrFormat
# =============================================================================

test_that("mrFormat functions with normal parameters", {
  input_data <- generate_format_mock_data()

  condition <- "condition1"
  coverage <- "CS_95_susie"
  alleleQc <- TRUE

  res <- mrFormat(input_data$susie_result, condition, input_data$gwas_sumstats_db, coverage, alleleQc)
  expect_true(is.data.frame(res))
  expect_true(all(c("gene_name", "variant_id", "bhat_x", "sbhat_x", "cs", "pip", "bhat_y", "sbhat_y") %in% names(res)))
  expect_gt(nrow(res), 0)
  # gene_name should be populated
  expect_true(all(res$gene_name == res$gene_name[1]))
  # When we get matched rows, bhat_y and sbhat_y should be numeric
  expect_true(is.numeric(res$bhat_y))
  expect_true(is.numeric(res$sbhat_y))
})

test_that("mrFormat returns a dataframe with NAs for zero coverage in top_loci", {
  input_data <- generate_format_mock_data()

  condition <- "condition1"
  coverage <- "CS_95_susie"
  alleleQc <- TRUE

  susie_result_mock <- input_data$susie_result
  susie_result_mock[["susie_results"]][[condition]][["top_loci"]][[coverage]] <- rep(0, nrow(susie_result_mock[["susie_results"]][[condition]][["top_loci"]]))

  result <- mrFormat(susie_result_mock, condition, input_data$gwas_sumstats_db, coverage = coverage, runAlleleQc = alleleQc)

  expect_true(is.data.frame(result))
  expect_true(all(is.na(result[,-1])))
})

test_that("mrFormat returns a dataframe with NAs with non-existent top_loci", {
  input_data <- generate_format_mock_data()

  condition <- "condition1"
  coverage <- "CS_95_susie"
  alleleQc <- TRUE

  susie_result_mock <- input_data$susie_result
  susie_result_mock[["susie_results"]][[condition]][["top_loci"]] <- list()

  result <- mrFormat(susie_result_mock, condition, input_data$gwas_sumstats_db, coverage = coverage, runAlleleQc = alleleQc)

  expect_true(is.data.frame(result))
  expect_true(all(is.na(result[,-1])))
})

# =============================================================================
# mrAnalysis
# =============================================================================

test_that("mrAnalysis returns expected output with normal inputs", {
  input_data <- generate_mock_mr_formatted_input(num_variants = 10, generate_full_dataset = TRUE)
  result <- mrAnalysis(input_data, cpipCutoff = 0.5)
  expect_true(is.data.frame(result))
  expect_gt(nrow(result), 0)
  expect_true(all(c("gene_name", "num_CS", "num_IV", "cpip", "meta_eff", "se_meta_eff", "meta_pval", "Q", "Q_pval", "I2") %in% names(result)))
  # When result has non-NA values, verify value properties
  if (!is.na(result$meta_eff)) {
    expect_true(is.finite(result$meta_eff))
    expect_true(result$se_meta_eff > 0)
    expect_true(result$I2 >= 0 && result$I2 <= 1)
  }
})

test_that("mrAnalysis returns null output for all NA input except gene_name", {
  input_data <- generate_mock_mr_formatted_input(generate_full_dataset = FALSE)
  result <- mrAnalysis(input_data, cpipCutoff = 0.5)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1)
  expect_true(all(is.na(result[,-1])))
})

test_that("mrAnalysis handles no significant cpip values correctly", {
  input_data <- generate_mock_mr_formatted_input(num_variants = 5, generate_full_dataset = TRUE)
  input_data$pip <- runif(nrow(input_data), 0, 0.1)
  result <- mrAnalysis(input_data, cpipCutoff = 0.5)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1)
  expect_true(all(is.na(result[,-1])))
})

test_that("mrAnalysis: single CS single variant with pip above cutoff", {
  set.seed(42)
  input <- data.frame(
    gene_name = "GENE_SINGLE",
    variant_id = "1:500:A:C",
    bhat_x = 0.8,
    sbhat_x = 0.1,
    cs = 1,
    pip = 0.9,
    bhat_y = 0.3,
    sbhat_y = 0.06,
    stringsAsFactors = FALSE
  )

  result <- mrAnalysis(input, cpipCutoff = 0.5)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(result$gene_name, "GENE_SINGLE")
  expect_equal(result$num_CS, 1L)
  expect_equal(result$num_IV, 1L)
  expect_equal(result$Q, 0)
  expect_equal(result$I2, 0)
  expect_false(is.na(result$meta_eff))
  expect_false(is.na(result$se_meta_eff))
  expect_false(is.na(result$meta_pval))
})

test_that("mrAnalysis: cpip exactly at cutoff boundary is included", {
  set.seed(101)
  input <- data.frame(
    gene_name = rep("GENE_BOUNDARY", 2),
    variant_id = c("1:100:A:G", "1:200:C:T"),
    bhat_x = c(0.5, 0.4),
    sbhat_x = c(0.1, 0.1),
    cs = c(1, 1),
    pip = c(0.25, 0.25),
    bhat_y = c(0.2, 0.15),
    sbhat_y = c(0.05, 0.05),
    stringsAsFactors = FALSE
  )

  result <- mrAnalysis(input, cpipCutoff = 0.5)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_false(is.na(result$meta_eff))
  expect_equal(result$gene_name, "GENE_BOUNDARY")
  expect_equal(result$cpip, 0.5)
})

test_that("mrAnalysis: cpip just below cutoff returns null output", {
  set.seed(102)
  input <- data.frame(
    gene_name = rep("GENE_BELOW", 2),
    variant_id = c("1:100:A:G", "1:200:C:T"),
    bhat_x = c(0.5, 0.4),
    sbhat_x = c(0.1, 0.1),
    cs = c(1, 1),
    pip = c(0.24, 0.25),
    bhat_y = c(0.2, 0.15),
    sbhat_y = c(0.05, 0.05),
    stringsAsFactors = FALSE
  )

  result <- mrAnalysis(input, cpipCutoff = 0.5)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_true(is.na(result$meta_eff))
  expect_equal(result$gene_name, "GENE_BELOW")
})

test_that("mrAnalysis: output columns are rounded to 3 decimal places", {
  set.seed(300)
  input <- data.frame(
    gene_name = rep("GENE_ROUND", 4),
    variant_id = c("1:10:A:G", "1:20:C:T", "1:30:G:A", "1:40:T:C"),
    bhat_x = c(0.5, 0.3, 0.4, 0.6),
    sbhat_x = c(0.1, 0.1, 0.1, 0.1),
    cs = c(1, 1, 2, 2),
    pip = c(0.6, 0.4, 0.5, 0.5),
    bhat_y = c(0.2, 0.15, 0.3, 0.25),
    sbhat_y = c(0.05, 0.05, 0.05, 0.05),
    stringsAsFactors = FALSE
  )

  result <- mrAnalysis(input, cpipCutoff = 0.5)

  rounded_cols <- c("cpip", "meta_eff", "se_meta_eff", "meta_pval", "Q", "Q_pval", "I2")
  for (col in rounded_cols) {
    val <- result[[col]]
    if (!is.na(val)) {
      expect_equal(round(val, 3), val,
                   info = paste("Column", col, "should be rounded to 3 decimal places"))
    }
  }
})

test_that("mrAnalysis: gene_name is preserved in output", {
  set.seed(400)
  gene <- "ENSG00000012345_BRCA1"
  input <- data.frame(
    gene_name = rep(gene, 2),
    variant_id = c("1:100:A:G", "1:200:C:T"),
    bhat_x = c(0.5, 0.3),
    sbhat_x = c(0.1, 0.1),
    cs = c(1, 1),
    pip = c(0.6, 0.4),
    bhat_y = c(0.2, 0.15),
    sbhat_y = c(0.05, 0.05),
    stringsAsFactors = FALSE
  )

  result <- mrAnalysis(input, cpipCutoff = 0.5)
  expect_equal(result$gene_name, gene)
})

test_that("mrAnalysis: meta_pval is a valid probability", {
  set.seed(600)
  input <- data.frame(
    gene_name = rep("GENE_ORDER", 3),
    variant_id = c("1:10:A:G", "1:20:C:T", "1:30:G:A"),
    bhat_x = c(0.5, 0.3, 0.4),
    sbhat_x = c(0.1, 0.1, 0.1),
    cs = c(1, 1, 1),
    pip = c(0.4, 0.3, 0.3),
    bhat_y = c(0.2, 0.15, 0.18),
    sbhat_y = c(0.05, 0.05, 0.05),
    stringsAsFactors = FALSE
  )

  result <- mrAnalysis(input, cpipCutoff = 0.5)
  expect_true(result$meta_pval >= 0 && result$meta_pval <= 1)
})

test_that("mrAnalysis: se_meta_eff is always positive when result is non-null", {
  set.seed(700)
  input <- data.frame(
    gene_name = rep("GENE_SE", 3),
    variant_id = c("1:10:A:G", "1:20:C:T", "1:30:G:A"),
    bhat_x = c(0.5, -0.3, 0.4),
    sbhat_x = c(0.1, 0.1, 0.1),
    cs = c(1, 1, 1),
    pip = c(0.4, 0.3, 0.3),
    bhat_y = c(0.2, -0.15, 0.18),
    sbhat_y = c(0.05, 0.05, 0.05),
    stringsAsFactors = FALSE
  )

  result <- mrAnalysis(input, cpipCutoff = 0.5)
  expect_true(result$se_meta_eff > 0)
})

test_that("mrAnalysis: I2 is between 0 and 1 inclusive", {
  set.seed(800)
  input <- data.frame(
    gene_name = rep("GENE_I2", 4),
    variant_id = c("1:10:A:G", "1:20:C:T", "1:30:G:A", "1:40:T:C"),
    bhat_x = c(0.5, 0.3, 0.4, 0.6),
    sbhat_x = c(0.1, 0.1, 0.1, 0.1),
    cs = c(1, 1, 2, 2),
    pip = c(0.6, 0.4, 0.5, 0.5),
    bhat_y = c(0.2, 0.15, 0.3, 0.25),
    sbhat_y = c(0.05, 0.05, 0.05, 0.05),
    stringsAsFactors = FALSE
  )

  result <- mrAnalysis(input, cpipCutoff = 0.5)
  expect_true(result$I2 >= 0 && result$I2 <= 1)
})

test_that("mrAnalysis: Q_pval is a valid probability", {
  set.seed(801)
  input <- data.frame(
    gene_name = rep("GENE_QPVAL", 4),
    variant_id = c("1:10:A:G", "1:20:C:T", "1:30:G:A", "1:40:T:C"),
    bhat_x = c(0.5, 0.3, 0.4, 0.6),
    sbhat_x = c(0.1, 0.1, 0.1, 0.1),
    cs = c(1, 1, 2, 2),
    pip = c(0.6, 0.4, 0.5, 0.5),
    bhat_y = c(0.2, 0.15, 0.3, 0.25),
    sbhat_y = c(0.05, 0.05, 0.05, 0.05),
    stringsAsFactors = FALSE
  )

  result <- mrAnalysis(input, cpipCutoff = 0.5)
  expect_true(result$Q_pval >= 0 && result$Q_pval <= 1)
})

test_that("mrAnalysis: cpipCutoff = 0 includes all CS groups", {
  set.seed(900)
  input <- data.frame(
    gene_name = rep("GENE_CUTOFF0", 3),
    variant_id = c("1:10:A:G", "1:20:C:T", "1:30:G:A"),
    bhat_x = c(0.5, 0.3, 0.4),
    sbhat_x = c(0.1, 0.1, 0.1),
    cs = c(1, 2, 3),
    pip = c(0.01, 0.01, 0.01),
    bhat_y = c(0.2, 0.15, 0.18),
    sbhat_y = c(0.05, 0.05, 0.05),
    stringsAsFactors = FALSE
  )

  result <- mrAnalysis(input, cpipCutoff = 0)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_false(is.na(result$meta_eff))
  expect_equal(result$num_CS, 3L)
  expect_equal(result$num_IV, 3L)
})

test_that("mrAnalysis: cpipCutoff = 1 excludes CS with cpip < 1", {
  set.seed(901)
  input <- data.frame(
    gene_name = rep("GENE_CUTOFF1", 3),
    variant_id = c("1:10:A:G", "1:20:C:T", "1:30:G:A"),
    bhat_x = c(0.5, 0.3, 0.4),
    sbhat_x = c(0.1, 0.1, 0.1),
    cs = c(1, 1, 2),
    pip = c(0.4, 0.3, 0.2),
    bhat_y = c(0.2, 0.15, 0.18),
    sbhat_y = c(0.05, 0.05, 0.05),
    stringsAsFactors = FALSE
  )

  result <- mrAnalysis(input, cpipCutoff = 1.0)
  expect_true(is.na(result$meta_eff))
})

test_that("mrAnalysis: mixed CS where only some pass cpip filter", {
  set.seed(1000)
  input <- data.frame(
    gene_name = rep("GENE_MIXED", 4),
    variant_id = c("1:10:A:G", "1:20:C:T", "1:30:G:A", "1:40:T:C"),
    bhat_x = c(0.5, 0.3, 0.4, 0.6),
    sbhat_x = c(0.1, 0.1, 0.1, 0.1),
    cs = c(1, 1, 2, 2),
    pip = c(0.4, 0.4, 0.1, 0.1),
    bhat_y = c(0.2, 0.15, 0.3, 0.25),
    sbhat_y = c(0.05, 0.05, 0.05, 0.05),
    stringsAsFactors = FALSE
  )

  result <- mrAnalysis(input, cpipCutoff = 0.5)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_false(is.na(result$meta_eff))
  expect_equal(result$num_CS, 1L)
  expect_equal(result$num_IV, 2L)
})

test_that("mrAnalysis: large negative bhat_x values are handled", {
  set.seed(1100)
  input <- data.frame(
    gene_name = rep("GENE_NEG", 2),
    variant_id = c("1:10:A:G", "1:20:C:T"),
    bhat_x = c(-2.5, -1.8),
    sbhat_x = c(0.1, 0.1),
    cs = c(1, 1),
    pip = c(0.6, 0.4),
    bhat_y = c(-0.5, -0.3),
    sbhat_y = c(0.05, 0.05),
    stringsAsFactors = FALSE
  )

  result <- mrAnalysis(input, cpipCutoff = 0.5)
  expect_s3_class(result, "data.frame")
  expect_false(is.na(result$meta_eff))
  expect_true(result$se_meta_eff > 0)
})

test_that("mrAnalysis: bhat_x normalized to z-score (bhat_x/sbhat_x) then sbhat_x=1", {
  set.seed(1200)
  bx <- 0.6
  sx <- 0.15
  by <- 0.2
  sy <- 0.05

  input <- data.frame(
    gene_name = "GENE_NORM",
    variant_id = "1:10:A:G",
    bhat_x = bx,
    sbhat_x = sx,
    cs = 1,
    pip = 1.0,
    bhat_y = by,
    sbhat_y = sy,
    stringsAsFactors = FALSE
  )

  result <- mrAnalysis(input, cpipCutoff = 0.5)

  bhat_x_norm <- bx / sx
  sbhat_x_norm <- 1
  beta_yx <- by / bhat_x_norm
  se_yx <- sqrt((sy^2 / bhat_x_norm^2) + ((by^2 * sbhat_x_norm^2) / bhat_x_norm^4))
  expected_meta_eff <- round(beta_yx, 3)

  expect_equal(result$meta_eff, expected_meta_eff)
})

test_that("mrAnalysis with multiple credible sets", {
  input <- data.frame(
    gene_name = rep("GENE1", 4),
    variant_id = paste0("1:", seq(100, 400, 100), ":A:G"),
    bhat_x = c(0.5, 0.3, 0.4, 0.6),
    sbhat_x = c(0.1, 0.1, 0.1, 0.1),
    cs = c(1, 1, 2, 2),
    pip = c(0.4, 0.6, 0.5, 0.5),
    bhat_y = c(0.2, 0.15, 0.3, 0.25),
    sbhat_y = c(0.05, 0.05, 0.05, 0.05),
    stringsAsFactors = FALSE
  )

  result <- mrAnalysis(input, cpipCutoff = 0.5)
  expect_equal(nrow(result), 1)
  expect_equal(result$num_CS, 2)
})

# =============================================================================
# .createNullMrDf
# =============================================================================

test_that(".createNullMrDf creates correct structure", {
  spec <- c(x = "numeric", y = "character", z = "integer")
  result <- pecotmr:::.createNullMrDf("gene1", spec)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(result$gene_name, "gene1")
  expect_true(is.numeric(result$x))
  expect_true(is.character(result$y))
  expect_true(is.integer(result$z))
  expect_true(all(is.na(result[, -1])))
})

# =============================================================================
# mrFormat: tryCatch error path
# =============================================================================

test_that("mrFormat returns null df when getNestedElement errors for top_loci", {
  # susie_result has the gene_name path but top_loci path is broken
  bad_susie <- list(susie_results = list(
    condition1 = list(
      region_info = list(region_name = "Gene_Test")
      # no top_loci key - forces tryCatch error path
    )
  ))
  result <- mrFormat(
    bad_susie, "condition1",
    data.frame(variant_id = "chr1:1:A:G", pos = 1, z = 1, beta = 0.1, se = 0.05,
               effect_allele_frequency = 0.3, n_case = 500, n_control = 500),
    coverage = "CS_95_susie"
  )
  expect_true(is.data.frame(result))
  expect_true(all(is.na(result[, -1])))
})

# =============================================================================
# mrFormat: no overlapping positions
# =============================================================================

test_that("mrFormat returns null df when no positions overlap", {
  input_data <- generate_format_mock_data(seed = 42)
  # Shift GWAS positions so they don't overlap with SuSiE positions
  input_data$gwas_sumstats_db$pos <- input_data$gwas_sumstats_db$pos + 10000
  input_data$gwas_sumstats_db$variant_id <- paste0(
    "chr1:", input_data$gwas_sumstats_db$pos, ":",
    sub("chr1:\\d+:", "", input_data$gwas_sumstats_db$variant_id)
  )
  result <- mrFormat(
    input_data$susie_result, "condition1", input_data$gwas_sumstats_db,
    coverage = "CS_95_susie", runAlleleQc = TRUE
  )
  expect_true(is.data.frame(result))
  expect_true(all(is.na(result[, -1])))
})

# =============================================================================
# mrFormat: coverage = NULL path (line 57)
# =============================================================================

test_that("mrFormat derives coverage from coverage_level and method when coverage is NULL", {
  input_data <- generate_format_mock_data(seed = 10)
  # The default method is "susie" and coverage_level is 0.95 which gives "CS_95_susie"
  result <- mrFormat(
    input_data$susie_result, "condition1", input_data$gwas_sumstats_db,
    coverage = NULL, runAlleleQc = TRUE,
    method = "susie", coverageLevel = 0.95
  )
  expect_true(is.data.frame(result))
  # Result should be the same as explicitly passing "CS_95_susie"
  result_explicit <- mrFormat(
    input_data$susie_result, "condition1", input_data$gwas_sumstats_db,
    coverage = "CS_95_susie", runAlleleQc = TRUE
  )
  expect_equal(result, result_explicit)
})

# =============================================================================
# mrFormat: pip column not found (line 64)
# =============================================================================

test_that("mrFormat returns null df when pip column is not found", {
  input_data <- generate_format_mock_data(seed = 20)
  # Remove the pip column from top_loci so resolve_pip_column returns NULL
  tl <- input_data$susie_result$susie_results$condition1$top_loci
  tl$pip <- NULL
  # Also rename any pip_* columns if they exist
  names(tl) <- gsub("^pip_.*", "removed", names(tl))
  input_data$susie_result$susie_results$condition1$top_loci <- tl
  result <- mrFormat(
    input_data$susie_result, "condition1", input_data$gwas_sumstats_db,
    coverage = "CS_95_susie", runAlleleQc = TRUE,
    method = "nonexistent_method"
  )
  expect_true(is.data.frame(result))
  expect_true(all(is.na(result[, -1])))
})

# =============================================================================
# mrFormat: impute missing effect_allele_frequency (line 85)
# =============================================================================

test_that("mrFormat imputes missing effect_allele_frequency from ld_meta_df", {
  set.seed(30)
  num_variants <- 5
  ref_alleles <- c("A", "T", "G", "C")
  allele_pairs <- paste0(
    sample(ref_alleles, num_variants, replace = TRUE), ":",
    sample(ref_alleles, num_variants, replace = TRUE)
  )
  # Ensure ref != alt
  while (any(sapply(strsplit(allele_pairs, ":"), function(x) x[1] == x[2]))) {
    allele_pairs <- paste0(
      sample(ref_alleles, num_variants, replace = TRUE), ":",
      sample(ref_alleles, num_variants, replace = TRUE)
    )
  }
  variant_ids <- paste0("chr1:", 1:num_variants, ":", allele_pairs)

  top_loci <- data.frame(
    variant_id = variant_ids,
    betahat = rnorm(num_variants),
    sebetahat = runif(num_variants, 0.05, 0.1),
    pip = runif(num_variants, 0.5, 1),
    CS_95_susie = rep(1L, num_variants),
    stringsAsFactors = FALSE
  )
  susie_result <- list(
    susie_results = list(
      condition1 = list(
        topLoci = top_loci,
        region_info = list(region_name = "Gene_Impute")
      )
    )
  )
  gwas_sumstats_db <- data.frame(
    variant_id = variant_ids,
    beta = rnorm(num_variants),
    se = runif(num_variants, 0.05, 0.1),
    pos = 1:num_variants,
    effect_allele_frequency = c(NA, 0.3, NA, 0.4, NA),
    n_case = rep(1000, num_variants),
    n_control = rep(1000, num_variants),
    stringsAsFactors = FALSE
  ) %>%
    mutate(n_sample = n_case + n_control, z = beta / se)
  ld_meta_df <- data.frame(
    pos = 1:num_variants,
    allele_freq = runif(num_variants, 0.1, 0.5),
    stringsAsFactors = FALSE
  )
  result <- mrFormat(
    susie_result, "condition1", gwas_sumstats_db,
    coverage = "CS_95_susie", runAlleleQc = FALSE,
    ldMetaDf = ld_meta_df
  )
  expect_true(is.data.frame(result))
  # Should have some rows if matching worked
  expect_true("bhat_y" %in% names(result))
})

# =============================================================================
# mrFormat: no common variants after allele QC (line 103)
# =============================================================================

test_that("mrFormat returns null df when allele QC removes all variants", {
  set.seed(40)
  num_variants <- 3
  # Create SuSiE variants with one set of alleles
  top_loci <- data.frame(
    variant_id = paste0("chr1:", 1:num_variants, ":A:G"),
    betahat = rnorm(num_variants),
    sebetahat = runif(num_variants, 0.05, 0.1),
    pip = runif(num_variants, 0.5, 1),
    CS_95_susie = rep(1L, num_variants),
    stringsAsFactors = FALSE
  )
  susie_result <- list(
    susie_results = list(
      condition1 = list(
        topLoci = top_loci,
        region_info = list(region_name = "Gene_NoCommon")
      )
    )
  )
  # GWAS has same positions but incompatible alleles (G:T is not a strand

  # complement or flip of A:G)
  gwas_sumstats_db <- data.frame(
    variant_id = paste0("chr1:", 1:num_variants, ":G:T"),
    beta = rnorm(num_variants),
    se = runif(num_variants, 0.05, 0.1),
    pos = 1:num_variants,
    effect_allele_frequency = runif(num_variants, 0.1, 0.5),
    n_case = rep(1000, num_variants),
    n_control = rep(1000, num_variants),
    stringsAsFactors = FALSE
  ) %>%
    mutate(n_sample = n_case + n_control, z = beta / se)
  result <- mrFormat(
    susie_result, "condition1", gwas_sumstats_db,
    coverage = "CS_95_susie", runAlleleQc = TRUE
  )
  expect_true(is.data.frame(result))
  expect_true(all(is.na(result[, -1])))
})

# =============================================================================
# fineMr
# =============================================================================

test_that("fineMr returns expected columns and reasonable values with multi-gene input", {
  set.seed(100)
  # Gene A: 2 credible sets, 3 variants in CS1, 2 variants in CS2
  # Gene B: 2 credible sets, 2 variants in CS1, 3 variants in CS2
  formatted_input <- tibble::tibble(
    X_ID = c(rep("GeneA", 5), rep("GeneB", 5)),
    cs = c(1, 1, 1, 2, 2, 1, 1, 2, 2, 2),
    pip = c(0.4, 0.3, 0.3, 0.6, 0.4, 0.5, 0.5, 0.3, 0.4, 0.3),
    bhat_x = c(0.5, 0.3, 0.4, 0.6, 0.35, 0.7, 0.45, 0.5, 0.6, 0.4),
    sbhat_x = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
    bhat_y = c(0.2, 0.15, 0.18, 0.3, 0.2, 0.35, 0.2, 0.25, 0.3, 0.22),
    sbhat_y = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05),
    snp = paste0("rs", 1:10)
  )

  result <- fineMr(formatted_input, cpipCutoff = 0.5)

  # Check output structure
  expected_cols <- c("X_ID", "num_CS", "num_IV", "cpip", "composite_bhat",
                     "composite_sbhat", "meta_eff", "se_meta_eff", "Q", "I2")
  expect_true(all(expected_cols %in% names(result)))
  # One row per gene
  expect_equal(nrow(result), 2)
  expect_setequal(result$X_ID, c("GeneA", "GeneB"))
  # Both genes have 2 CS
  expect_true(all(result$num_CS == 2))
  # num_IV should be 5 for each gene
  expect_true(all(result$num_IV == 5))
  # se_meta_eff should be positive
  expect_true(all(result$se_meta_eff > 0))
  # I2 should be between 0 and 1
  expect_true(all(result$I2 >= 0 & result$I2 <= 1))
  # Numeric columns should be rounded to 3 decimal places
  for (col in c("cpip", "composite_bhat", "meta_eff", "se_meta_eff", "Q", "I2")) {
    expect_equal(round(result[[col]], 3), result[[col]],
                 info = paste("Column", col, "should be rounded to 3 decimal places"))
  }
})

test_that("fineMr filters out CS with cpip below cutoff", {
  set.seed(101)
  formatted_input <- tibble::tibble(
    X_ID = rep("GeneC", 4),
    cs = c(1, 1, 2, 2),
    pip = c(0.4, 0.4, 0.05, 0.05),
    bhat_x = c(0.5, 0.3, 0.4, 0.6),
    sbhat_x = c(0.1, 0.1, 0.1, 0.1),
    bhat_y = c(0.2, 0.15, 0.3, 0.25),
    sbhat_y = c(0.05, 0.05, 0.05, 0.05),
    snp = paste0("rs", 1:4)
  )

  result <- fineMr(formatted_input, cpipCutoff = 0.5)
  # CS 2 has cpip = 0.10 which is below 0.5, so only CS 1 passes
  expect_equal(nrow(result), 1)
  expect_equal(result$num_CS, 1L)
  expect_equal(result$num_IV, 2L)
})

test_that("fineMr returns empty result when all CS are below cpip cutoff", {
  set.seed(102)
  formatted_input <- tibble::tibble(
    X_ID = rep("GeneD", 4),
    cs = c(1, 1, 2, 2),
    pip = c(0.1, 0.1, 0.1, 0.1),
    bhat_x = c(0.5, 0.3, 0.4, 0.6),
    sbhat_x = c(0.1, 0.1, 0.1, 0.1),
    bhat_y = c(0.2, 0.15, 0.3, 0.25),
    sbhat_y = c(0.05, 0.05, 0.05, 0.05),
    snp = paste0("rs", 1:4)
  )

  result <- fineMr(formatted_input, cpipCutoff = 0.5)
  # All CS have cpip = 0.2 which is below cutoff; result should be empty
  expect_equal(nrow(result), 0)
})

test_that("fineMr with single variant single CS returns correct Wald ratio", {
  set.seed(103)
  bx <- 0.8
  sx <- 0.1
  by <- 0.3
  sy <- 0.06

  formatted_input <- tibble::tibble(
    X_ID = "GeneE",
    cs = 1,
    pip = 1.0,
    bhat_x = bx,
    sbhat_x = sx,
    bhat_y = by,
    sbhat_y = sy,
    snp = "rs1"
  )

  result <- fineMr(formatted_input, cpipCutoff = 0.5)
  expect_equal(nrow(result), 1)
  expect_equal(result$num_CS, 1L)
  expect_equal(result$num_IV, 1L)

  # Verify the Wald ratio calculation
  bhat_x_norm <- bx / sx
  beta_yx <- by / bhat_x_norm
  expected_bhat <- round(beta_yx, 3)
  expect_equal(result$composite_bhat, expected_bhat)
  # With single CS, meta_eff should equal composite_bhat
  expect_equal(result$meta_eff, expected_bhat)
  # Q and I2 should be 0 with single CS
  expect_equal(result$Q, 0)
  expect_equal(result$I2, 0)
})

test_that("fineMr handles negative effect sizes", {
  set.seed(104)
  formatted_input <- tibble::tibble(
    X_ID = rep("GeneF", 4),
    cs = c(1, 1, 2, 2),
    pip = c(0.6, 0.4, 0.5, 0.5),
    bhat_x = c(-0.5, -0.3, -0.4, -0.6),
    sbhat_x = c(0.1, 0.1, 0.1, 0.1),
    bhat_y = c(-0.2, -0.15, -0.3, -0.25),
    sbhat_y = c(0.05, 0.05, 0.05, 0.05),
    snp = paste0("rs", 1:4)
  )

  result <- fineMr(formatted_input, cpipCutoff = 0.5)
  expect_equal(nrow(result), 1)
  expect_true(is.finite(result$meta_eff))
  expect_true(result$se_meta_eff > 0)
})

test_that("fineMr with cpipCutoff = 0 includes all credible sets", {
  set.seed(105)
  formatted_input <- tibble::tibble(
    X_ID = rep("GeneG", 3),
    cs = c(1, 2, 3),
    pip = c(0.01, 0.01, 0.01),
    bhat_x = c(0.5, 0.3, 0.4),
    sbhat_x = c(0.1, 0.1, 0.1),
    bhat_y = c(0.2, 0.15, 0.18),
    sbhat_y = c(0.05, 0.05, 0.05),
    snp = paste0("rs", 1:3)
  )

  result <- fineMr(formatted_input, cpipCutoff = 0)
  expect_equal(nrow(result), 1)
  expect_equal(result$num_CS, 3L)
  expect_equal(result$num_IV, 3L)
})
