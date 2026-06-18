library(testthat)

# build_twas_score_row is an unexported contract helper shared between
# twasPipeline (pecotmr) and quantile_twas_pipeline (qQTLR). Tests cover the
# shapes returned by twasAnalysis() so both pipelines stay aligned.

test_that("build_twas_score_row returns empty data.frame for NULL input", {
  out <- pecotmr:::buildTwasScoreRow(NULL,
                                        weightDb ="ENSG001",
                                        context   = "Cortex",
                                        study     = "AD")
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0)
})

test_that("build_twas_score_row packs a single-method twas_rs", {
  # Shape mimics twasAnalysis() output: apply(weights, 2, twasZ) returns
  # a named list where each element is list(z=, pval=). The name's final
  # "_<suffix>" is the method marker (see the sub() regex).
  # findData(twas_rs, c(2, "z")) descends one level and pulls out "z".
  twas_rs <- list(
    enetWeights = list(z = 2.5, pval = 0.012)
  )
  out <- pecotmr:::buildTwasScoreRow(twas_rs,
                                        weightDb ="ENSG001",
                                        context   = "Cortex",
                                        study     = "AD")
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1)
  expect_setequal(colnames(out),
                  c("gwasStudy", "method", "twasZ", "twasPval",
                    "context", "molecularId"))
  expect_equal(out$gwasStudy,   "AD")
  expect_equal(out$method,       "enet")
  expect_equal(out$twasZ,       2.5)
  expect_equal(out$twasPval,    0.012)
  expect_equal(out$context,      "Cortex")
  expect_equal(out$molecularId, "ENSG001")
})

test_that("build_twas_score_row packs a multi-method twas_rs", {
  twas_rs <- list(
    enetWeights  = list(z =  2.5, pval = 0.012),
    lassoWeights = list(z = -1.8, pval = 0.072),
    top_weights   = list(z =  0.4, pval = 0.689)
  )
  out <- pecotmr:::buildTwasScoreRow(twas_rs,
                                        weightDb ="ENSG002",
                                        context   = "Liver",
                                        study     = "T2D")
  expect_equal(nrow(out), 3)
  expect_equal(out$method,       c("enet", "lasso", "top"))
  expect_equal(out$twasZ,       c(2.5, -1.8, 0.4))
  expect_equal(out$twasPval,    c(0.012, 0.072, 0.689))
  expect_true(all(out$gwasStudy   == "T2D"))
  expect_true(all(out$context      == "Liver"))
  expect_true(all(out$molecularId == "ENSG002"))
})
