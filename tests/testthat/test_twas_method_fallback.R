context("TWAS method fallback")

# Helper to build a minimal twas_table for testing apply_method_fallback
make_twas_table <- function(methods, twas_z_values, rsq_values, is_selected, gwas_study = "study1") {
  data.frame(
    molecular_id = "gene1",
    context = "ctx1",
    gwas_study = gwas_study,
    method = methods,
    is_selected_method = is_selected,
    is_imputable = TRUE,
    rsq_cv = rsq_values,
    pval_cv = 0.01,
    twas_z = twas_z_values,
    twas_pval = ifelse(is.na(twas_z_values) | !is.finite(twas_z_values), NA, 0.05),
    type = "eQTL",
    chr = 1,
    block = "chr1_100_200",
    stringsAsFactors = FALSE
  )
}

# Access the internal function
apply_fallback <- pecotmr:::applyMethodFallback

test_that("no fallback when selected method has valid z", {
  df <- make_twas_table(
    methods = c("susie", "enet", "lasso"),
    twas_z_values = c(2.5, 1.8, 1.2),
    rsq_values = c(0.3, 0.2, 0.1),
    is_selected = c(TRUE, FALSE, FALSE)
  )
  result <- apply_fallback(df)
  expect_equal(result$method[result$is_selected_method], "susie")
})

test_that("fallback to next best method when selected has NA z", {
  df <- make_twas_table(
    methods = c("susie", "enet", "lasso"),
    twas_z_values = c(NA, 1.8, 1.2),
    rsq_values = c(0.3, 0.2, 0.1),
    is_selected = c(TRUE, FALSE, FALSE)
  )
  result <- apply_fallback(df)
  expect_equal(result$method[result$is_selected_method], "enet")
  expect_false(result$is_selected_method[result$method == "susie"])
})

test_that("fallback to next best method when selected has Inf z", {
  df <- make_twas_table(
    methods = c("susie", "enet", "lasso"),
    twas_z_values = c(Inf, 1.8, 1.2),
    rsq_values = c(0.3, 0.2, 0.1),
    is_selected = c(TRUE, FALSE, FALSE)
  )
  result <- apply_fallback(df)
  expect_equal(result$method[result$is_selected_method], "enet")
})

test_that("fallback picks highest rsq among valid candidates", {
  df <- make_twas_table(
    methods = c("susie", "enet", "lasso"),
    twas_z_values = c(NA, 1.8, 2.0),
    rsq_values = c(0.3, 0.1, 0.25),
    is_selected = c(TRUE, FALSE, FALSE)
  )
  result <- apply_fallback(df)
  # lasso has higher rsq (0.25) than enet (0.1)
  expect_equal(result$method[result$is_selected_method], "lasso")
})

test_that("all methods NA sets is_imputable to FALSE", {
  df <- make_twas_table(
    methods = c("susie", "enet", "lasso"),
    twas_z_values = c(NA, NA, NA),
    rsq_values = c(0.3, 0.2, 0.1),
    is_selected = c(TRUE, FALSE, FALSE)
  )
  result <- apply_fallback(df)
  expect_true(all(!result$is_imputable))
})

test_that("fallback is per-study: one study needs fallback, another does not", {
  df1 <- make_twas_table(
    methods = c("susie", "enet"),
    twas_z_values = c(NA, 1.5),
    rsq_values = c(0.3, 0.2),
    is_selected = c(TRUE, FALSE),
    gwas_study = "study1"
  )
  df2 <- make_twas_table(
    methods = c("susie", "enet"),
    twas_z_values = c(2.5, 1.8),
    rsq_values = c(0.3, 0.2),
    is_selected = c(TRUE, FALSE),
    gwas_study = "study2"
  )
  df <- rbind(df1, df2)
  result <- apply_fallback(df)
  # study1: fallback to enet
  s1 <- result[result$gwas_study == "study1", ]
  expect_equal(s1$method[s1$is_selected_method], "enet")
  # study2: no change
  s2 <- result[result$gwas_study == "study2", ]
  expect_equal(s2$method[s2$is_selected_method], "susie")
})

test_that("fallback handles empty data frame", {
  df <- data.frame(
    molecularId = character(), context = character(), gwas_study = character(),
    method = character(), is_selected_method = logical(), is_imputable = logical(),
    rsq_cv = numeric(), pval_cv = numeric(), twasZ = numeric(), twas_pval = numeric(),
    stringsAsFactors = FALSE
  )
  result <- apply_fallback(df)
  expect_equal(nrow(result), 0)
})

test_that("fallback handles -Inf z", {
  df <- make_twas_table(
    methods = c("susie", "enet"),
    twas_z_values = c(-Inf, 1.5),
    rsq_values = c(0.3, 0.2),
    is_selected = c(TRUE, FALSE)
  )
  result <- apply_fallback(df)
  expect_equal(result$method[result$is_selected_method], "enet")
})
