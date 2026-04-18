library(testthat)

# build_twas_score_row is an unexported contract helper shared between
# twas_pipeline (pecotmr) and quantile_twas_pipeline (qQTLR). Tests cover the
# shapes returned by twas_analysis() so both pipelines stay aligned.

test_that("build_twas_score_row returns empty data.frame for NULL input", {
  out <- pecotmr:::build_twas_score_row(NULL,
                                        weight_db = "ENSG001",
                                        context   = "Cortex",
                                        study     = "AD")
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0)
})

test_that("build_twas_score_row packs a single-method twas_rs", {
  # Shape mimics twas_analysis() output: apply(weights, 2, twas_z) returns
  # a named list where each element is list(z=, pval=). The name's final
  # "_<suffix>" is the method marker (see the sub() regex).
  # find_data(twas_rs, c(2, "z")) descends one level and pulls out "z".
  twas_rs <- list(
    enet_weights = list(z = 2.5, pval = 0.012)
  )
  out <- pecotmr:::build_twas_score_row(twas_rs,
                                        weight_db = "ENSG001",
                                        context   = "Cortex",
                                        study     = "AD")
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1)
  expect_setequal(colnames(out),
                  c("gwas_study", "method", "twas_z", "twas_pval",
                    "context", "molecular_id"))
  expect_equal(out$gwas_study,   "AD")
  expect_equal(out$method,       "enet")
  expect_equal(out$twas_z,       2.5)
  expect_equal(out$twas_pval,    0.012)
  expect_equal(out$context,      "Cortex")
  expect_equal(out$molecular_id, "ENSG001")
})

test_that("build_twas_score_row packs a multi-method twas_rs", {
  twas_rs <- list(
    enet_weights  = list(z =  2.5, pval = 0.012),
    lasso_weights = list(z = -1.8, pval = 0.072),
    top_weights   = list(z =  0.4, pval = 0.689)
  )
  out <- pecotmr:::build_twas_score_row(twas_rs,
                                        weight_db = "ENSG002",
                                        context   = "Liver",
                                        study     = "T2D")
  expect_equal(nrow(out), 3)
  expect_equal(out$method,       c("enet", "lasso", "top"))
  expect_equal(out$twas_z,       c(2.5, -1.8, 0.4))
  expect_equal(out$twas_pval,    c(0.012, 0.072, 0.689))
  expect_true(all(out$gwas_study   == "T2D"))
  expect_true(all(out$context      == "Liver"))
  expect_true(all(out$molecular_id == "ENSG002"))
})
