# Tests for R/SldscData.R (S4 container + accessors)

test_that("SldscData constructs from in-memory objects", {
  sd <- .sldscMkData()
  expect_s4_class(sd, "SldscData")
  expect_true(methods::validObject(sd))
})

test_that("SldscData defaults frq to a 0-row data.frame when NULL", {
  sd <- .sldscMkData(withFrq = FALSE)
  expect_s4_class(sd, "SldscData")
  expect_equal(nrow(getFrqData(sd)), 0L)
})

test_that("SldscData errors when `annot` is missing", {
  expect_error(SldscData(), "`annot` is required")
})

# ---- validity ----

test_that("validity rejects annot without CHR/SNP", {
  expect_error(
    SldscData(annot = data.frame(SNP = "rs1", annot_A = 1)),
    "must have columns CHR and SNP")
})

test_that("validity rejects annot with no annotation column", {
  expect_error(
    SldscData(annot = data.frame(CHR = 1, SNP = "rs1")),
    "at least one annotation column")
})

test_that("validity rejects a non-empty frq without SNP/MAF", {
  expect_error(
    SldscData(annot = data.frame(CHR = 1, SNP = "rs1", a = 1),
              frq = data.frame(CHR = 1, foo = 0.2)),
    "non-empty `frq` must have columns SNP and MAF")
})

test_that("validity rejects unnamed traits", {
  run <- .sldscMkRun(c("annot_A_0"))
  expect_error(
    SldscData(annot = data.frame(CHR = 1, SNP = "rs1", a = 1),
              traits = list(list(single = list(run)))),
    "must be a named list")
})

test_that("validity rejects a trait without a `single` element", {
  run <- .sldscMkRun(c("annot_A_0"))
  expect_error(
    SldscData(annot = data.frame(CHR = 1, SNP = "rs1", a = 1),
              traits = list(traitX = list(joint = run))),
    "must be a list with a `single` element")
})

# ---- accessors ----

test_that("getAnnotData / getFrqData return the stored frames", {
  sd <- .sldscMkData()
  expect_s3_class(getAnnotData(sd), "data.frame")
  expect_equal(nrow(getAnnotData(sd)), 6L)
  expect_equal(nrow(getFrqData(sd)), 6L)
})

test_that("getAnnotCols returns the annotation columns only", {
  sd <- .sldscMkData()
  expect_equal(getAnnotCols(sd), c("annot_A", "annot_B"))
})

test_that("getTraitRuns / getTraitNames expose the traits list", {
  sd <- .sldscMkData()
  expect_equal(getTraitNames(sd), c("traitX", "traitY"))
  expect_named(getTraitRuns(sd), c("traitX", "traitY"))
})

test_that("getTraitRun retrieves single (by idx), joint, and NULL cases", {
  sd <- .sldscMkData()
  single1 <- getTraitRun(sd, "traitX", "single", 1L)
  expect_equal(single1$categories, c("annot_A_0", "baselineLD_0"))
  expect_equal(getTraitRun(sd, "traitX", "joint")$categories,
               c("annot_A_0", "annot_B_0", "baselineLD_0"))
  # whole single list when idx is NULL
  expect_length(getTraitRun(sd, "traitX", "single"), 2L)
  # out-of-range idx -> NULL
  expect_null(getTraitRun(sd, "traitX", "single", 5L))
  # unknown trait -> NULL
  expect_null(getTraitRun(sd, "no_such_trait", "joint"))
})

test_that("getTraitRun returns NULL joint when a trait has none", {
  sd <- .sldscMkData(withJoint = FALSE)
  expect_null(getTraitRun(sd, "traitX", "joint"))
})

test_that("show prints a compact summary", {
  sd <- .sldscMkData()
  out <- capture.output(show(sd))
  expect_true(any(grepl("SldscData", out)))
  expect_true(any(grepl("annot_A, annot_B", out)))
  expect_true(any(grepl("traitX, traitY", out)))
})
