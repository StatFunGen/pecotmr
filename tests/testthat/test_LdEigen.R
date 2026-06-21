# Tests migrated from test_h2ClassesSumstats.R

# === Tests migrated from test_h2ClassesSumstats.R (LdEigen) ===

test_that("LdEigen constructs and validates correctly", {
  ldblocks <- make_test_ldblocks()
  snp_info <- make_test_snp_info()
  eigen_list <- list(
    list(values = c(1, 0.5), vectors = matrix(rnorm(20), 10, 2),
         snpIdx = 1:10),
    list(values = c(0.8), vectors = matrix(rnorm(10), 10, 1),
         snpIdx = 1:10)
  )

  obj <- new("LdEigen",
    ldBlocks = ldblocks,
    snpInfo = snp_info,
    nRef = 500L,
    inSample = FALSE,
    genome = "hg19",
    eigenList = eigen_list,
    eigenvalueTruncation = 0.9
  )
  expect_s4_class(obj, "LdEigen")
  expect_true(methods::validObject(obj))
})


test_that("LdEigen rejects eigen_list length mismatch", {
  ldblocks <- make_test_ldblocks()  # 2 blocks
  # Only 1 element in eigen_list
  expect_error(
    methods::validObject(
      new("LdEigen",
        ldBlocks = ldblocks,
        snpInfo = make_test_snp_info(),
        nRef = 500L,
        inSample = FALSE,
        genome = "hg19",
        eigenList = list(list(values = 1)),
        eigenvalueTruncation = 0.9
      )
    ),
    "eigenList.*must match"
  )
})


test_that("LdEigen rejects invalid eigenvalue_truncation", {
  ldblocks <- make_test_ldblocks()
  expect_error(
    methods::validObject(
      new("LdEigen",
        ldBlocks = ldblocks,
        snpInfo = make_test_snp_info(),
        nRef = 500L,
        inSample = FALSE,
        genome = "hg19",
        eigenList = list(list(), list()),
        eigenvalueTruncation = 0
      )
    ),
    "eigenvalueTruncation"
  )
})


