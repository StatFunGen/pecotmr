context("univariate_rss_diagnostics")

.test_fm_result <- function(variant_names, trimmed_fit = list(),
                            top_loci = data.frame(variant_id = character(0),
                                                  method = character(0),
                                                  stringsAsFactors = FALSE)) {
  FineMappingResult(
    variant_names = variant_names,
    trimmed_fit = trimmed_fit,
    top_loci = top_loci,
    method = "susie"
  )
}

# ===========================================================================
# get_susie_result
# ===========================================================================

test_that("get_susie_result returns NULL for empty input", {
  result <- get_susie_result(list())
  expect_null(result)
})

test_that("get_susie_result returns NULL when finemapping_result missing", {
  result <- get_susie_result(list(some_data = 42))
  expect_null(result)
})

test_that("get_susie_result returns trimmed result when present", {
  mock_result <- list(pip = c(0.1, 0.5, 0.3), sets = list(cs = list()))
  con_data <- list(finemapping_result = .test_fm_result(
    variant_names = character(0),
    trimmed_fit = mock_result
  ))
  result <- get_susie_result(con_data)
  expect_equal(result, mock_result)
})

# ===========================================================================
# extract_top_pip_info
# ===========================================================================

test_that("extract_top_pip_info finds top PIP variant", {
  con_data <- list(
    finemapping_result = .test_fm_result(
      variant_names = c("1:100:A:G", "1:200:C:T", "1:300:G:A"),
      trimmed_fit = list(pip = c(0.1, 0.7, 0.2))
    ),
    sumstats = list(z = c(1.0, 3.5, -0.5))
  )
  result <- extract_top_pip_info(con_data)
  expect_equal(result$top_variant, "1:200:C:T")
  expect_equal(result$top_pip, 0.7)
  expect_equal(result$top_z, 3.5)
  expect_equal(result$top_variant_index, 2)
  expect_true(is.na(result$cs_name))
  expect_true(is.na(result$variants_per_cs))
})

test_that("extract_top_pip_info computes p_value from z", {
  con_data <- list(
    finemapping_result = .test_fm_result(
      variant_names = c("1:100:A:G", "1:200:C:T", "1:300:G:A"),
      trimmed_fit = list(pip = c(0.9, 0.05, 0.05))
    ),
    sumstats = list(z = c(5.0, 0.5, -0.3))
  )
  result <- extract_top_pip_info(con_data)
  expected_pval <- z_to_pvalue(5.0)
  expect_equal(result$p_value, expected_pval)
})

test_that("extract_top_pip_info handles ties by taking first max", {
  con_data <- list(
    finemapping_result = .test_fm_result(
      variant_names = c("1:100:A:G", "1:200:C:T", "1:300:G:A"),
      trimmed_fit = list(pip = c(0.5, 0.5, 0.5))
    ),
    sumstats = list(z = c(1.0, 2.0, 3.0))
  )
  result <- extract_top_pip_info(con_data)
  expect_equal(result$top_variant_index, 1)
  expect_equal(result$top_pip, 0.5)
})

# ===========================================================================
# extract_cs_info
# ===========================================================================

test_that("extract_cs_info extracts single CS correctly", {
  con_data <- list(
    finemapping_result = .test_fm_result(
      variant_names = c("1:100:A:G", "1:200:C:T", "1:300:G:A"),
      trimmed_fit = list(
        sets = list(cs = list(L_1 = c(1, 2))),
        cs_corr = NULL
      )
    )
  )
  top_loci_table <- data.frame(
    variant_id = c("1:100:A:G", "1:200:C:T"),
    pip = c(0.3, 0.8),
    z = c(2.0, 4.5),
    stringsAsFactors = FALSE
  )
  result <- extract_cs_info(con_data, cs_names = "L_1", top_loci_table = top_loci_table)
  expect_equal(nrow(result), 1)
  expect_equal(result$cs_name, "L_1")
  expect_equal(result$top_variant, "1:200:C:T")
  expect_equal(result$top_pip, 0.8)
  expect_equal(result$variants_per_cs, 2)
  expect_true(grepl("NA", result$cs_corr[[1]]))
})

test_that("extract_cs_info extracts multiple CSs with cs_corr", {
  con_data <- list(
    finemapping_result = .test_fm_result(
      variant_names = c("1:100:A:G", "1:200:C:T", "1:300:G:A", "1:400:T:C"),
      trimmed_fit = list(
        sets = list(
          cs = list(L_1 = c(1, 2), L_2 = c(3, 4))
        ),
        cs_corr = matrix(c(1, 0.3, 0.3, 1), nrow = 2)
      )
    )
  )
  top_loci_table <- data.frame(
    variant_id = c("1:100:A:G", "1:200:C:T", "1:300:G:A", "1:400:T:C"),
    pip = c(0.3, 0.8, 0.6, 0.1),
    z = c(2.0, 4.5, 3.0, 0.5),
    stringsAsFactors = FALSE
  )
  result <- extract_cs_info(con_data, cs_names = c("L_1", "L_2"), top_loci_table = top_loci_table)
  expect_equal(nrow(result), 2)
  expect_equal(result$cs_name[1], "L_1")
  expect_equal(result$cs_name[2], "L_2")
  expect_equal(result$top_variant[1], "1:200:C:T")
  expect_equal(result$top_variant[2], "1:300:G:A")
  expect_true(is.character(result$cs_corr[[1]]))
})

test_that("extract_cs_info computes p_value from z-score", {
  con_data <- list(
    finemapping_result = .test_fm_result(
      variant_names = c("1:100:A:G", "1:200:C:T"),
      trimmed_fit = list(
        sets = list(cs = list(L_1 = c(1, 2))),
        cs_corr = NULL
      )
    )
  )
  top_loci_table <- data.frame(
    variant_id = c("1:100:A:G", "1:200:C:T"),
    pip = c(0.9, 0.1),
    z = c(5.0, 0.5),
    stringsAsFactors = FALSE
  )
  result <- extract_cs_info(con_data, cs_names = "L_1", top_loci_table = top_loci_table)
  expected_pval <- z_to_pvalue(5.0)
  expect_equal(result$p_value, expected_pval, tolerance = 1e-10)
})

# ===========================================================================
# parse_cs_corr
# ===========================================================================

test_that("parse_cs_corr handles NA correlations", {
  df <- data.frame(
    cs_name = "L1",
    top_pip = 0.9,
    cs_corr = NA_character_,
    stringsAsFactors = FALSE
  )
  result <- parse_cs_corr(df)
  expect_true("cs_corr_max" %in% colnames(result))
  expect_true("cs_corr_min" %in% colnames(result))
  expect_true(is.na(result$cs_corr_max))
  expect_true(is.na(result$cs_corr_min))
})

test_that("parse_cs_corr splits comma-separated correlations", {
  df <- data.frame(
    cs_name = c("L1", "L2"),
    top_pip = c(0.9, 0.3),
    cs_corr = c("1,0.3", "0.3,1"),
    stringsAsFactors = FALSE
  )
  result <- parse_cs_corr(df)
  expect_true("cs_corr_1" %in% colnames(result))
  expect_true("cs_corr_2" %in% colnames(result))
  expect_equal(result$cs_corr_max[1], 0.3)
  expect_equal(result$cs_corr_min[1], 0.3)
})

test_that("parse_cs_corr handles empty string", {
  df <- data.frame(
    cs_name = "L1",
    cs_corr = "",
    stringsAsFactors = FALSE
  )
  result <- parse_cs_corr(df)
  expect_true(is.na(result$cs_corr_max))
})

test_that("parse_cs_corr handles multiple correlations", {
  df <- data.frame(
    cs_name = "L1",
    cs_corr = "1,0.5,0.2",
    stringsAsFactors = FALSE
  )
  result <- parse_cs_corr(df)
  expect_equal(result$cs_corr_max, 0.5)
  expect_equal(result$cs_corr_min, 0.2)
  expect_equal(result$cs_corr_1, 1)
  expect_equal(result$cs_corr_2, 0.5)
  expect_equal(result$cs_corr_3, 0.2)
})

test_that("parse_cs_corr handles NULL value", {
  df <- data.frame(
    cs_name = "L1",
    cs_corr = NA,
    stringsAsFactors = FALSE
  )
  result <- parse_cs_corr(df)
  expect_true(is.na(result$cs_corr_max))
  expect_true(is.na(result$cs_corr_min))
})

test_that("parse_cs_corr handles single value without comma", {
  df <- data.frame(
    cs_name = "L1",
    cs_corr = "0.5",
    stringsAsFactors = FALSE
  )
  result <- parse_cs_corr(df)
  expect_true(is.na(result$cs_corr_max))
  expect_true(is.na(result$cs_corr_min))
})

test_that("parse_cs_corr handles all-1 correlations (self-corr only)", {
  df <- data.frame(
    cs_name = c("L1", "L2"),
    cs_corr = c("1,1", "1,1"),
    stringsAsFactors = FALSE
  )
  result <- parse_cs_corr(df)
  expect_true(is.na(result$cs_corr_max[1]))
  expect_true(is.na(result$cs_corr_min[1]))
})

test_that("parse_cs_corr handles mixed valid and NA rows", {
  df <- data.frame(
    cs_name = c("L1", "L2"),
    cs_corr = c("1,0.3,0.7", NA_character_),
    stringsAsFactors = FALSE
  )
  result <- parse_cs_corr(df)
  expect_equal(result$cs_corr_max[1], 0.7)
  expect_equal(result$cs_corr_min[1], 0.3)
  expect_true(is.na(result$cs_corr_max[2]))
})

test_that("parse_cs_corr expands columns for different lengths", {
  df <- data.frame(
    cs_name = c("L1", "L2", "L3"),
    cs_corr = c("1,0.3,0.5", "1,0.2", "1,0.8,0.4,0.1"),
    stringsAsFactors = FALSE
  )
  result <- parse_cs_corr(df)
  expect_true("cs_corr_4" %in% colnames(result))
  expect_true(is.na(result$cs_corr_4[1]))
  expect_equal(result$cs_corr_4[3], 0.1)
})

# ===========================================================================
# auto_decision
# ===========================================================================

test_that("auto_decision assigns BVSR when no CS is tagged", {
  df <- data.frame(
    cs_name = c("L1", "L2"),
    top_z = c(5.0, 3.5),
    p_value = c(1e-10, 1e-6),
    stringsAsFactors = FALSE
  )
  result <- auto_decision(df, high_corr_cols = character(0))
  expect_true("top_cs" %in% colnames(result))
  expect_true("tagged_cs" %in% colnames(result))
  expect_true("method" %in% colnames(result))
  expect_true(all(result$method == "BVSR"))
})

test_that("auto_decision assigns SER when all non-top CSs are tagged", {
  df <- data.frame(
    cs_name = c("L1", "L2"),
    top_z = c(5.0, 0.1),
    p_value = c(1e-10, 0.5),
    stringsAsFactors = FALSE
  )
  result <- auto_decision(df, high_corr_cols = character(0))
  expect_true(result$top_cs[1])
  expect_false(result$top_cs[2])
  expect_true(result$tagged_cs[2])
  expect_true(all(result$method == "SER"))
})

test_that("auto_decision assigns BCR when untagged CS remain", {
  df <- data.frame(
    cs_name = c("L1", "L2", "L3"),
    top_z = c(5.0, 3.5, 0.1),
    p_value = c(1e-10, 1e-6, 0.5),
    stringsAsFactors = FALSE
  )
  result <- auto_decision(df, high_corr_cols = character(0))
  expect_true(all(result$method == "BCR"))
})

test_that("auto_decision assigns SER for single CS", {
  df <- data.frame(
    cs_name = "L1",
    top_z = 5.0,
    p_value = 1e-10,
    stringsAsFactors = FALSE
  )
  result <- auto_decision(df, high_corr_cols = character(0))
  expect_true(result$top_cs[1])
  expect_false(result$tagged_cs[1])
  expect_equal(result$method, "SER")
})
