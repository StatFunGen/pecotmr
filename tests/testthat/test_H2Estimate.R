# Tests migrated from test_h2ClassesSumstats.R

# === Tests migrated from test_h2ClassesSumstats.R (H2Estimate) ===

test_that("H2Estimate constructs with all slots", {
  obj <- new("H2Estimate",
    h2 = 0.3,
    h2Se = 0.05,
    intercept = 1.01,
    interceptSe = 0.02,
    local = NULL,
    enrichment = NULL,
    tauBlocks = NULL,
    scoreStats = NULL,
    method = "lder",
    nSnps = 10000L,
    traitName = "height"
  )
  expect_s4_class(obj, "H2Estimate")
  expect_equal(obj@h2, 0.3)
  expect_equal(obj@method, "lder")
})


