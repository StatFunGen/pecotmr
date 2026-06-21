# Tests migrated from test_h2ClassesSumstats.R

# === Tests migrated from test_h2ClassesSumstats.R (LdScore) ===

test_that("LdScore constructs and validates correctly", {
  ldblocks <- make_test_ldblocks()
  n <- 10
  snp_info <- make_test_snp_info(n)

  obj <- new("LdScore",
    ldBlocks = ldblocks,
    snpInfo = snp_info,
    nRef = 500L,
    inSample = FALSE,
    genome = "hg19",
    ldScores = matrix(runif(n), nrow = n, ncol = 1),
    ldScoreWeights = runif(n),
    ldMatrixList = list()
  )
  expect_s4_class(obj, "LdScore")
  expect_true(methods::validObject(obj))
})


test_that("LdScore rejects ld_scores row mismatch with snp_info", {
  ldblocks <- make_test_ldblocks()
  snp_info <- make_test_snp_info(10)

  expect_error(
    methods::validObject(
      new("LdScore",
        ldBlocks = ldblocks,
        snpInfo = snp_info,
        nRef = 500L,
        inSample = FALSE,
        genome = "hg19",
        ldScores = matrix(0, nrow = 5, ncol = 1),  # wrong rows
        ldScoreWeights = runif(10),
        ldMatrixList = list()
      )
    ),
    "ldScores.*must match"
  )
})


