# Tests migrated from test_h2ClassesSumstats.R

# === Tests migrated from test_h2ClassesSumstats.R (LdBlocks) ===

test_that("LdBlocks constructs and validates correctly", {
  obj <- make_test_ldblocks()
  expect_s4_class(obj, "LdBlocks")
  expect_equal(length(obj@blocks), 2)
  expect_equal(obj@genome, "hg19")
  expect_true(methods::validObject(obj))
})


test_that("LdBlocks rejects genome of length != 1", {
  blocks_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1, end = 5000)
  )
  expect_error(
    methods::validObject(
      new("LdBlocks", blocks = blocks_gr, genome = c("hg19", "hg38"))
    ),
    "genome.*single"
  )
})


