context("TWAS method fallback")

# Helper to build a minimal twasTable for testing apply_method_fallback
make_twas_table <- function(methods, twas_z_values, rsq_values, is_selected, gwasStudy = "study1") {
  data.frame(
    molecularId = "gene1",
    context = "ctx1",
    gwasStudy = gwasStudy,
    method = methods,
    isSelectedMethod = is_selected,
    isImputable = TRUE,
    rsqCv = rsq_values,
    pvalCv = 0.01,
    twasZ = twas_z_values,
    twasPval = ifelse(is.na(twas_z_values) | !is.finite(twas_z_values), NA, 0.05),
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
  expect_equal(result$method[result$isSelectedMethod], "susie")
})

test_that("fallback to next best method when selected has NA z", {
  df <- make_twas_table(
    methods = c("susie", "enet", "lasso"),
    twas_z_values = c(NA, 1.8, 1.2),
    rsq_values = c(0.3, 0.2, 0.1),
    is_selected = c(TRUE, FALSE, FALSE)
  )
  result <- apply_fallback(df)
  expect_equal(result$method[result$isSelectedMethod], "enet")
  expect_false(result$isSelectedMethod[result$method == "susie"])
})

test_that("fallback to next best method when selected has Inf z", {
  df <- make_twas_table(
    methods = c("susie", "enet", "lasso"),
    twas_z_values = c(Inf, 1.8, 1.2),
    rsq_values = c(0.3, 0.2, 0.1),
    is_selected = c(TRUE, FALSE, FALSE)
  )
  result <- apply_fallback(df)
  expect_equal(result$method[result$isSelectedMethod], "enet")
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
  expect_equal(result$method[result$isSelectedMethod], "lasso")
})

test_that("all methods NA sets isImputable to FALSE", {
  df <- make_twas_table(
    methods = c("susie", "enet", "lasso"),
    twas_z_values = c(NA, NA, NA),
    rsq_values = c(0.3, 0.2, 0.1),
    is_selected = c(TRUE, FALSE, FALSE)
  )
  result <- apply_fallback(df)
  expect_true(all(!result$isImputable))
})

test_that("fallback is per-study: one study needs fallback, another does not", {
  df1 <- make_twas_table(
    methods = c("susie", "enet"),
    twas_z_values = c(NA, 1.5),
    rsq_values = c(0.3, 0.2),
    is_selected = c(TRUE, FALSE),
    gwasStudy = "study1"
  )
  df2 <- make_twas_table(
    methods = c("susie", "enet"),
    twas_z_values = c(2.5, 1.8),
    rsq_values = c(0.3, 0.2),
    is_selected = c(TRUE, FALSE),
    gwasStudy = "study2"
  )
  df <- rbind(df1, df2)
  result <- apply_fallback(df)
  # study1: fallback to enet
  s1 <- result[result$gwasStudy == "study1", ]
  expect_equal(s1$method[s1$isSelectedMethod], "enet")
  # study2: no change
  s2 <- result[result$gwasStudy == "study2", ]
  expect_equal(s2$method[s2$isSelectedMethod], "susie")
})

test_that("fallback handles empty data frame", {
  df <- data.frame(
    molecularId = character(), context = character(), gwasStudy = character(),
    method = character(), isSelectedMethod = logical(), isImputable = logical(),
    rsqCv = numeric(), pvalCv = numeric(), twasZ = numeric(), twasPval = numeric(),
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
  expect_equal(result$method[result$isSelectedMethod], "enet")
})
