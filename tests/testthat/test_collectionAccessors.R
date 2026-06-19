context("Collection-level accessors")

# ===========================================================================
# Shared fixture
# ===========================================================================

.ca_makeTopLoci <- function(n = 3, withCs = TRUE) {
  tl <- data.frame(
    variant_id = paste0("chr1:", 100 * seq_len(n), ":A:G"),
    pip        = seq(0.9, by = -0.1, length.out = n),
    stringsAsFactors = FALSE)
  if (withCs) tl$cs <- c(1L, 1L, 0L)[seq_len(n)]
  tl
}

.ca_makeFmEntry <- function(n = 3) {
  FineMappingEntry(
    variantIds = paste0("chr1:", 100 * seq_len(n), ":A:G"),
    trimmedFit = list(payload = sprintf("fit_n=%d", n)),
    topLoci    = .ca_makeTopLoci(n))
}

# ===========================================================================
# QtlFineMappingResult collection accessors
# ===========================================================================

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
  e <- .ca_makeFmEntry(3)  # cs = c(1, 1, 0)
  res <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(e))
  cs <- getCs(res)
  expect_equal(nrow(cs), 2L)
  expect_true(all(cs$cs > 0))
})

test_that("QtlFineMappingResult: getTopLoci returns the entry's topLoci", {
  e <- .ca_makeFmEntry(3)
  res <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(e))
  expect_equal(getTopLoci(res), .ca_makeTopLoci(3))
})

test_that("QtlFineMappingResult: getTrimmedFit reads the entry's trimmedFit", {
  e <- .ca_makeFmEntry(3)
  res <- QtlFineMappingResult(
    study = "s1", context = "c1", trait = "t1", method = "susie",
    entry = list(e))
  expect_equal(getTrimmedFit(res), list(payload = "fit_n=3"))
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

# ===========================================================================
# GwasFineMappingResult collection accessors
# ===========================================================================

test_that("GwasFineMappingResult: getPip with study/method selectors", {
  e1 <- .ca_makeFmEntry(3)
  e2 <- .ca_makeFmEntry(4)
  res <- GwasFineMappingResult(
    study  = c("g1", "g2"),
    method = c("susie", "susie"),
    entry  = list(e1, e2))
  pip <- getPip(res, study = "g2", method = "susie")
  expect_equal(length(pip), 4L)
})

test_that("GwasFineMappingResult: getContexts/getTraits return NULL", {
  e <- .ca_makeFmEntry(3)
  res <- GwasFineMappingResult(study = "g1", method = "susie",
                                entry = list(e))
  expect_null(getContexts(res))
  expect_null(getTraits(res))
})

test_that("GwasFineMappingResult: getCs/getTopLoci/getTrimmedFit/getVariantIds dispatch", {
  e <- .ca_makeFmEntry(3)
  res <- GwasFineMappingResult(study = "g1", method = "susie",
                                entry = list(e))
  expect_equal(nrow(getCs(res)), 2L)
  expect_equal(getTopLoci(res), .ca_makeTopLoci(3))
  expect_equal(getTrimmedFit(res), list(payload = "fit_n=3"))
  expect_equal(length(getVariantIds(res)), 3L)
})

test_that("GwasFineMappingResult: .tupleSelectRowGwasFmr requires both selectors for multi-row", {
  e <- .ca_makeFmEntry(3)
  res <- GwasFineMappingResult(
    study  = c("g1", "g2"),
    method = c("susie", "susie"),
    entry  = list(e, e))
  expect_error(getPip(res),
               "Pass `study` and `method`")
  expect_error(getPip(res, study = c("g1", "g2"), method = "susie"),
               "must each be length 1")
  expect_error(getPip(res, study = "ghost", method = "susie"),
               "No entry for")
})

test_that("GwasFineMappingResult: getStudy/getMethodNames inherit from base", {
  e <- .ca_makeFmEntry(3)
  res <- GwasFineMappingResult(
    study  = c("g1", "g2"),
    method = c("susie", "susieRss"),
    entry  = list(e, e))
  expect_setequal(getStudy(res), c("g1", "g2"))
  expect_setequal(getMethodNames(res), c("susie", "susieRss"))
})
