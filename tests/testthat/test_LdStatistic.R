# Tests for R/LdStatistic.R (virtual base class)
# getGenome() is defined on the virtual LdStatistic and inherited by its
# concrete subclasses (LdEigen / LdScore); exercise it through a concrete
# LdScore instance. Fixtures (make_test_ldblocks / make_test_snp_info) come
# from helper-h2Classes.R.

test_that("getGenome returns the genome build string (via an LdScore subclass)", {
  n <- 10
  obj <- new("LdScore",
    ldBlocks = make_test_ldblocks(),
    snpInfo = make_test_snp_info(n),
    nRef = 500L,
    inSample = FALSE,
    genome = "hg19",
    ldScores = matrix(runif(n), nrow = n, ncol = 1),
    ldScoreWeights = runif(n),
    ldMatrixList = list())
  expect_equal(getGenome(obj), "hg19")
})
