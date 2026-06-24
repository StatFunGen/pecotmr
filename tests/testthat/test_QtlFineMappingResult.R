# Tests for R/QtlFineMappingResult.R

# === Tests migrated from test_s4Constructors.R (QtlFineMappingResult) ===

test_that("QtlFineMappingResult: builds a collection keyed by 4-tuple", {
  e1 <- .sc_makeFineMappingEntry(3)
  e2 <- .sc_makeFineMappingEntry(3)
  res <- QtlFineMappingResult(
    study   = c("s1", "s1"),
    context = c("c1", "c2"),
    trait   = c("t1", "t1"),
    method  = c("susie", "susie"),
    entry   = list(e1, e2))
  expect_s4_class(res, "QtlFineMappingResult")
  expect_equal(nrow(res), 2L)
  expect_null(res@ldSketch)
})


test_that("QtlFineMappingResult: stores an LD sketch when supplied", {
  e <- .sc_makeFineMappingEntry(3)
  gh <- .sc_makeGenotypeHandle()
  res <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(e), ldSketch = gh)
  expect_identical(getLdSketch(res), gh)
})


test_that("QtlFineMappingResult: errors on length mismatch", {
  e <- .sc_makeFineMappingEntry(3)
  expect_error(
    QtlFineMappingResult(
      study   = c("s1", "s2"),
      context = c("c1"),
      trait   = c("t1"),
      method  = c("susie"),
      entry   = list(e)),
    "same length"
  )
})


test_that("QtlFineMappingResult: validity rejects duplicate 4-tuples", {
  e1 <- .sc_makeFineMappingEntry(3)
  e2 <- .sc_makeFineMappingEntry(3)
  expect_error(
    QtlFineMappingResult(
      study   = c("s1", "s1"),
      context = c("c1", "c1"),
      trait   = c("t1", "t1"),
      method  = c("susie", "susie"),
      entry   = list(e1, e2)),
    "uniqueness violated"
  )
})


test_that("QtlFineMappingResult: getFineMappingResult returns selected entry", {
  e1 <- .sc_makeFineMappingEntry(3)
  e2 <- .sc_makeFineMappingEntry(3)
  res <- QtlFineMappingResult(
    study   = c("s1", "s1"),
    context = c("c1", "c2"),
    trait   = c("t1", "t1"),
    method  = c("susie", "susie"),
    entry   = list(e1, e2))
  picked <- getFineMappingResult(res,
                                 study = "s1", context = "c2",
                                 trait = "t1", method = "susie")
  expect_identical(picked, e2)
})


test_that("QtlFineMappingResult: getFineMappingResult errors on missing tuple", {
  e <- .sc_makeFineMappingEntry(3)
  res <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(e))
  expect_error(
    getFineMappingResult(res, study = "ghost", context = "c1",
                         trait = "t1", method = "susie"),
    "No entry for"
  )
})


test_that("QtlFineMappingResult: single-row collection allows omitting selectors", {
  e <- .sc_makeFineMappingEntry(3)
  res <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(e))
  expect_identical(getFineMappingResult(res), e)
})


test_that("QtlFineMappingResult: show prints summary", {
  e <- .sc_makeFineMappingEntry(3)
  res <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(e))
  expect_output(show(res), "QtlFineMappingResult")
})


test_that("QtlFineMappingResult: joint columns absent by default", {
  e <- .sc_makeFineMappingEntry(3)
  res <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(e))
  expect_false("jointStudies"  %in% names(res))
  expect_false("jointContexts" %in% names(res))
  expect_false("jointTraits"   %in% names(res))
})


test_that("QtlFineMappingResult: accepts jointContexts column", {
  e <- .sc_makeFineMappingEntry(3)
  res <- QtlFineMappingResult(
    study   = c("s1", "s1"),
    context = c("c1", "joint"),
    trait   = c("t1", "t1"),
    method  = c("susie", "mvsusie"),
    entry   = list(e, e),
    jointContexts = c(NA_character_, "c1;c2"))
  expect_true("jointContexts" %in% names(res))
  expect_identical(res$jointContexts, c(NA_character_, "c1;c2"))
})


test_that("QtlFineMappingResult: jointStudies + jointTraits combine cleanly", {
  e <- .sc_makeFineMappingEntry(3)
  res <- QtlFineMappingResult(
    study   = c("s1", "joint", "joint"),
    context = c("c1", "c1", "joint"),
    trait   = c("t1", "joint", "t1"),
    method  = c("susie", "mvsusieRss", "mvsusie"),
    entry   = list(e, e, e),
    jointStudies  = c(NA_character_, "s1;s2", NA_character_),
    jointContexts = c(NA_character_, NA_character_, "c1;c2"),
    jointTraits   = c(NA_character_, "t1;t2", NA_character_))
  expect_equal(nrow(res), 3L)
  expect_identical(res$jointStudies,
                   c(NA_character_, "s1;s2", NA_character_))
  expect_identical(res$jointTraits,
                   c(NA_character_, "t1;t2", NA_character_))
})


test_that("QtlFineMappingResult: uniqueness distinguishes joint members", {
  e <- .sc_makeFineMappingEntry(3)
  # same 4-tuple but different jointContexts -> distinct
  res <- QtlFineMappingResult(
    study   = c("s1", "s1"),
    context = c("joint", "joint"),
    trait   = c("t1", "t1"),
    method  = c("mvsusie", "mvsusie"),
    entry   = list(e, e),
    jointContexts = c("c1;c2", "c1;c3"))
  expect_equal(nrow(res), 2L)
  # same 4-tuple AND same jointContexts -> duplicate
  expect_error(
    QtlFineMappingResult(
      study   = c("s1", "s1"),
      context = c("joint", "joint"),
      trait   = c("t1", "t1"),
      method  = c("mvsusie", "mvsusie"),
      entry   = list(e, e),
      jointContexts = c("c1;c2", "c1;c2")),
    "uniqueness violated")
})


test_that("QtlFineMappingResult: length-mismatched joint vector errors", {
  e <- .sc_makeFineMappingEntry(3)
  expect_error(
    QtlFineMappingResult(
      study = "s1", context = "c1", trait = "t1", method = "susie",
      entry = list(e),
      jointContexts = c("c1;c2", "c1;c3")),
    "same length"
  )
})

# ===========================================================================
# GwasFineMappingResult
# ===========================================================================



# === Tests migrated from test_showMethods.R (QtlFineMappingResult) ===

test_that("show.QtlFineMappingResult prints entry/study/context/trait/method counts", {
  res <- QtlFineMappingResult(
    study   = c("s1", "s1"),
    context = c("c1", "c2"),
    trait   = c("t1", "t1"),
    method  = c("susie", "susieRss"),
    entry   = list(.sh_makeFmEntry(), .sh_makeFmEntry()))
  out <- capture.output(show(res))
  expect_true(any(grepl("QtlFineMappingResult: 2 entries", out)))
  expect_true(any(grepl("1 studies.*2 contexts.*1 traits.*2 methods", out)))
  expect_true(any(grepl("LD sketch: NULL", out)))
})


test_that("show.QtlFineMappingResult reports the ldSketch source when present", {
  res <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(.sh_makeFmEntry()),
    ldSketch = .sh_makeGenotypeHandle())
  out <- capture.output(show(res))
  expect_true(any(grepl("LD sketch: gds @ /tmp/test.gds", out)))
})



# === Tests migrated from test_collectionAccessors.R (QtlFineMappingResult) ===

test_that("QtlFineMappingResult: getPip returns named pip vector for selected tuple", {
  e1 <- .ca_makeFmEntry(3)
  e2 <- .ca_makeFmEntry(4)
  res <- QtlFineMappingResult(
    study   = c("s1", "s1"),
    context = c("c1", "c2"),
    trait   = c("t1", "t1"),
    method  = c("susie", "susie"),
    entry   = list(e1, e2))
  pip <- getPip(res, study = "s1", context = "c2", trait = "t1", method = "susie")
  expect_equal(length(pip), 4L)
  expect_equal(names(pip)[1L], "chr1:100:A:G")
})


test_that("QtlFineMappingResult: getPip(returnList = TRUE) wraps in pipe-keyed list", {
  e <- .ca_makeFmEntry(3)
  res <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(e))
  out <- getPip(res, returnList = TRUE)
  expect_true(is.list(out))
  expect_equal(names(out), "s1|c1|t1|susie")
})


test_that("QtlFineMappingResult: getCs filters to credible sets", {
  e <- .ca_makeFmEntry(3)  # cs_95 = c("susie_1", "susie_1", "susie_0")
  res <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(e))
  cs <- getCs(res)
  expect_equal(nrow(cs), 2L)
})


test_that("QtlFineMappingResult: getTopLoci returns the entry's topLoci (projected)", {
  e <- .ca_makeFmEntry(3)
  res <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(e))
  tl <- getTopLoci(res, signalCutoff = 0)
  expect_equal(nrow(tl), 3L)
  expect_equal(tl$variant_id, .ca_makeTopLoci(3)$variant_id)
})


test_that("QtlFineMappingResult: getSusieFit reads the entry's trimmedFit", {
  e <- .ca_makeFmEntry(3)
  res <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(e))
  expect_equal(getSusieFit(res), list(payload = "fit_n=3"))
})


test_that("QtlFineMappingResult: getVariantIds reads the entry's variantIds", {
  e <- .ca_makeFmEntry(3)
  res <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(e))
  expect_equal(length(getVariantIds(res)), 3L)
})


test_that("QtlFineMappingResult: getStudy/getContexts/getTraits/getMethodNames are unique", {
  e <- .ca_makeFmEntry(3)
  res <- QtlFineMappingResult(
    study   = c("s1", "s1", "s2"),
    context = c("c1", "c2", "c1"),
    trait   = c("t1", "t1", "t1"),
    method  = c("susie", "susieRss", "susie"),
    entry   = list(e, e, e))
  expect_setequal(getStudy(res), c("s1", "s2"))
  expect_setequal(getContexts(res), c("c1", "c2"))
  expect_equal(getTraits(res), "t1")
  expect_setequal(getMethodNames(res), c("susie", "susieRss"))
})


test_that("getCvResult works at the QtlFineMappingResult collection level", {
  tl <- data.frame(variant_id = "v1", pip = 0.5, stringsAsFactors = FALSE)
  cv <- list(samplePartition = data.frame(Sample = "s1", Fold = 1L),
             prediction = list(susie_predicted = matrix(0, 1, 1)),
             performance = list(susie_performance = matrix(0, 1, 6)))
  e <- FineMappingEntry("v1", list(), tl, cvResult = cv)
  fmr <- QtlFineMappingResult(study = "S", context = "C", trait = "T",
                              method = "susie", entry = list(e))
  expect_identical(
    getCvResult(fmr, study = "S", context = "C", trait = "T", method = "susie"),
    cv)
})

# ===========================================================================
# GwasFineMappingResult collection accessors
# ===========================================================================


