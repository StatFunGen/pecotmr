# Cross-class smoke test: confirms every S4 class with a show() method
# emits some expected output. Per-class show() behaviour tests live in
# the corresponding per-class test files (test_<Class>.R); this file is
# the single deliberate exception to the strict test_X.R <-> R/X.R naming
# rule because it exercises many classes in one pass.

test_that("show() methods do not error", {
  # LdBlocks
  expect_output(show(make_test_ldblocks()), "LdBlocks")

  # GenotypeHandle
  gh <- new("GenotypeHandle",
    path = "/tmp/test.gds",
    format = "gds",
    snpInfo = make_test_snp_info(),
    nSamples = 100L,
    sampleIds = paste0("s", 1:100),
    pgenPtr = NULL
  )
  expect_output(show(gh), "GenotypeHandle")

  # GwasSumStats (via constructor)
  ss <- makeGwasSumStatsFromDf(make_test_sumstats_df(10))
  expect_output(show(ss), "GwasSumStats")

  # AnnotationMatrix
  am <- AnnotationMatrix(
    matrix(0, nrow = 10, ncol = 1),
    make_test_granges(10),
    data.frame(name = "base", tier = "baseline", type = "binary",
               stringsAsFactors = FALSE)
  )
  expect_output(show(am), "AnnotationMatrix")

  # LdEigen
  ldblocks <- make_test_ldblocks()
  eig <- new("LdEigen",
    ldBlocks = ldblocks,
    snpInfo = make_test_snp_info(),
    nRef = 500L,
    inSample = FALSE,
    genome = "hg19",
    eigenList = list(list(), list()),
    eigenvalueTruncation = 0.9
  )
  expect_output(show(eig), "LdEigen")

  # LdScore
  n <- 10
  lsr <- new("LdScore",
    ldBlocks = ldblocks,
    snpInfo = make_test_snp_info(n),
    nRef = 500L,
    inSample = FALSE,
    genome = "hg19",
    ldScores = matrix(1, nrow = n, ncol = 1),
    ldScoreWeights = rep(1, n),
    ldMatrixList = list()
  )
  expect_output(show(lsr), "LdScore")

  # H2Estimate
  h2 <- new("H2Estimate",
    h2 = 0.3, h2Se = 0.05,
    intercept = 1.0, interceptSe = 0.01,
    local = NULL, enrichment = NULL,
    tauBlocks = NULL, scoreStats = NULL,
    method = "lder", nSnps = 1000L, traitName = "test"
  )
  expect_output(show(h2), "H2Estimate")
})
