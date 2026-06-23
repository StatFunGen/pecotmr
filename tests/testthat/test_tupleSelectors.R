context("tupleSelectors (internal row-selector helpers)")

# These helpers (.matchTupleRows, .tupleSelectRow, .tupleSelectRowGwasFmr)
# work on anything with `nrow(x)` and `x[[col]]` semantics, so the tests
# below use plain base-R data.frames rather than building S4 collections.

# ===========================================================================
# .matchTupleRows
# ===========================================================================

test_that(".matchTupleRows: empty keys returns every row index", {
  df <- data.frame(study = c("s1", "s2"), method = c("susie", "lasso"))
  expect_equal(pecotmr:::.matchTupleRows(df, list()),
               c(1L, 2L))
})

test_that(".matchTupleRows: AND-matches across multiple (column, value) pairs", {
  df <- data.frame(study   = c("s1", "s1", "s2"),
                   context = c("c1", "c2", "c1"),
                   stringsAsFactors = FALSE)
  expect_equal(pecotmr:::.matchTupleRows(df, list(study = "s1")),
               c(1L, 2L))
  expect_equal(
    pecotmr:::.matchTupleRows(df, list(study = "s1", context = "c2")),
    2L)
  expect_equal(
    pecotmr:::.matchTupleRows(df, list(study = "ghost", context = "c1")),
    integer(0))
})

# ===========================================================================
# .tupleSelectRow (QtlFineMappingResult / TwasWeights shape)
# ===========================================================================

test_that(".tupleSelectRow: zero-row input errors with the class label", {
  empty <- data.frame(study = character(0), context = character(0),
                      trait = character(0), method = character(0),
                      stringsAsFactors = FALSE)
  expect_error(
    pecotmr:::.tupleSelectRow(empty,
      study = "s1", context = "c1", trait = "t1", method = "susie",
      cls = "TwasWeights"),
    "TwasWeights has no rows")
})

test_that(".tupleSelectRow: single-row collection returns 1L without selectors", {
  one <- data.frame(study = "s1", context = "c1", trait = "t1",
                    method = "susie", stringsAsFactors = FALSE)
  expect_equal(pecotmr:::.tupleSelectRow(one), 1L)
})

test_that(".tupleSelectRow: multi-row + missing selectors errors with row count", {
  multi <- data.frame(study   = c("s1", "s1"),
                      context = c("c1", "c2"),
                      trait   = c("t1", "t1"),
                      method  = c("susie", "susie"),
                      stringsAsFactors = FALSE)
  expect_error(
    pecotmr:::.tupleSelectRow(multi, cls = "QtlFineMappingResult"),
    "QtlFineMappingResult has 2 entries")
})

test_that(".tupleSelectRow: non-scalar selectors error", {
  multi <- data.frame(study   = c("s1", "s2"),
                      context = c("c1", "c2"),
                      trait   = c("t1", "t2"),
                      method  = c("susie", "susie"),
                      stringsAsFactors = FALSE)
  expect_error(
    pecotmr:::.tupleSelectRow(multi,
      study = c("s1", "s2"), context = "c1", trait = "t1", method = "susie"),
    "must each be length 1")
})

test_that(".tupleSelectRow: matching tuple returns first row index", {
  multi <- data.frame(study   = c("s1", "s1"),
                      context = c("c1", "c2"),
                      trait   = c("t1", "t1"),
                      method  = c("susie", "susie"),
                      stringsAsFactors = FALSE)
  expect_equal(
    pecotmr:::.tupleSelectRow(multi,
      study = "s1", context = "c2", trait = "t1", method = "susie"),
    2L)
})

test_that(".tupleSelectRow: missing tuple errors with the 4-tuple in the message", {
  multi <- data.frame(study   = c("s1", "s1"),
                      context = c("c1", "c2"),
                      trait   = c("t1", "t1"),
                      method  = c("susie", "susie"),
                      stringsAsFactors = FALSE)
  expect_error(
    pecotmr:::.tupleSelectRow(multi,
      study = "ghost", context = "c1", trait = "t1", method = "susie"),
    "No entry for")
})

# ===========================================================================
# .tupleSelectRowGwasFmr
# ===========================================================================

test_that(".tupleSelectRowGwasFmr: zero-row input errors", {
  empty <- data.frame(study = character(0), method = character(0),
                      region_id = character(0), stringsAsFactors = FALSE)
  expect_error(
    pecotmr:::.tupleSelectRowGwasFmr(empty, study = "g1", method = "susie"),
    "GwasFineMappingResult has no rows")
})

test_that(".tupleSelectRowGwasFmr: single-row collection returns 1L", {
  one <- data.frame(study = "g1", method = "susie", region_id = "region_1",
                    stringsAsFactors = FALSE)
  expect_equal(pecotmr:::.tupleSelectRowGwasFmr(one), 1L)
})

test_that(".tupleSelectRowGwasFmr: missing selectors on multi-row errors", {
  multi <- data.frame(study     = c("g1", "g2"),
                      method    = c("susie", "susie"),
                      region_id = c("region_1", "region_1"),
                      stringsAsFactors = FALSE)
  expect_error(
    pecotmr:::.tupleSelectRowGwasFmr(multi),
    "Pass `study` and `method`")
})

test_that(".tupleSelectRowGwasFmr: non-scalar region errors", {
  multi <- data.frame(study     = c("g1", "g1"),
                      method    = c("susie", "susie"),
                      region_id = c("r1", "r2"),
                      stringsAsFactors = FALSE)
  expect_error(
    pecotmr:::.tupleSelectRowGwasFmr(multi,
      study = "g1", method = "susie", region = c("r1", "r2")),
    "`region` must be length 1")
})

test_that(".tupleSelectRowGwasFmr: region disambiguates per-block rows", {
  # Same (study, method) across two regions; region picks the right row.
  multi <- data.frame(study     = c("g1", "g1"),
                      method    = c("susie", "susie"),
                      region_id = c("chr22_1_100", "chr22_500_600"),
                      stringsAsFactors = FALSE)
  expect_equal(
    pecotmr:::.tupleSelectRowGwasFmr(multi,
      study = "g1", method = "susie", region = "chr22_500_600"),
    2L)
})

test_that(".tupleSelectRowGwasFmr: missing tuple errors and includes region in message", {
  one <- data.frame(study = "g1", method = "susie", region_id = "r1",
                    stringsAsFactors = FALSE)
  expect_error(
    pecotmr:::.tupleSelectRowGwasFmr(one,
      study = "g1", method = "susie", region = "ghost"),
    "region='ghost'")
})

test_that(".tupleSelectRowGwasFmr: ambiguous multi-match (no region) lists candidates", {
  # Two rows share (study, method); .tupleSelectRowGwasFmr should error
  # listing the available region_ids since the caller didn't disambiguate.
  multi <- data.frame(study     = c("g1", "g1"),
                      method    = c("susie", "susie"),
                      region_id = c("region_A", "region_B"),
                      stringsAsFactors = FALSE)
  expect_error(
    pecotmr:::.tupleSelectRowGwasFmr(multi, study = "g1", method = "susie"),
    "pass `region` to disambiguate")
})
