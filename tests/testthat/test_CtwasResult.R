context("CtwasResult")

# CtwasResultEntry (the per-row payload class) is tested in
# test_CtwasResultEntry.R; this file covers the CtwasResult collection.

.cr_entry <- function(ids = c("g1", "g2"), pip = c(0.9, 0.1), prior = NULL) {
  CtwasResultEntry(
    finemap    = data.frame(id = ids, susie_pip = pip, stringsAsFactors = FALSE),
    susieAlpha = data.frame(id = ids, susie_alpha = pip, stringsAsFactors = FALSE),
    param      = prior,
    regionInfo = data.frame(region_id = "r1", stringsAsFactors = FALSE))
}

test_that("CtwasResult: constructor keyed by (gwasStudy, study, context, method)", {
  cr <- CtwasResult(
    gwasStudy = c("D1", "D1"), study = c("Q1", "Q1"),
    context = c("brain", "liver"), method = c("susie", "susie"),
    entry = list(.cr_entry(), .cr_entry(ids = "g3", pip = 0.5)))
  expect_s4_class(cr, "CtwasResult")
  expect_equal(nrow(cr), 2L)
  expect_equal(getMethodNames(cr), "susie")
  expect_setequal(getContexts(cr), c("brain", "liver"))
  expect_equal(getStudy(cr), "Q1")
})

test_that("CtwasResult: getFinemap aggregates rows tagged with run identity", {
  cr <- CtwasResult(
    gwasStudy = c("D1", "D1"), study = c("Q1", "Q1"),
    context = c("brain", "liver"), method = c("susie", "susie"),
    entry = list(.cr_entry(ids = c("g1", "g2")), .cr_entry(ids = "g3", pip = 0.5)))
  fm <- getFinemap(cr)
  expect_true(all(c("gwasStudy", "study", "context", "method", "id", "susie_pip")
                  %in% names(fm)))
  expect_equal(nrow(fm), 3L)                 # 2 + 1 gene rows
  expect_equal(fm$context, c("brain", "brain", "liver"))
})

test_that("CtwasResult: getSusieAlpha aggregates the per-effect table with identity", {
  cr <- CtwasResult(
    gwasStudy = c("D1", "D1"), study = c("Q1", "Q1"),
    context = c("brain", "liver"), method = c("susie", "susie"),
    entry = list(.cr_entry(ids = c("g1", "g2")), .cr_entry(ids = "g3", pip = 0.5)))
  sa <- getSusieAlpha(cr)
  expect_equal(nrow(sa), 3L)
  expect_true(all(c("gwasStudy", "study", "context", "method", "susie_alpha")
                  %in% names(sa)))
  expect_equal(sa$context, c("brain", "brain", "liver"))
})

test_that("CtwasResult: uniqueness includes joint columns", {
  # same (gwasStudy, study, context, method) collides without joint columns
  expect_error(
    CtwasResult(gwasStudy = c("D1", "D1"), study = c("Q1", "Q1"),
                context = c("brain", "brain"), method = c("susie", "susie"),
                entry = list(.cr_entry(), .cr_entry())),
    "uniqueness violated")
  # a single-context row and a multi-context (jointContexts) row over the same
  # context are distinct once jointContexts joins the key
  cr <- CtwasResult(gwasStudy = c("D1", "D1"), study = c("Q1", "Q1"),
                    context = c("brain", "brain"), method = c("susie", "susie"),
                    entry = list(.cr_entry(), .cr_entry()),
                    jointContexts = c(NA, "brain,liver"))
  expect_equal(nrow(cr), 2L)
})

test_that("CtwasResult: rejects a non-CtwasResultEntry payload", {
  expect_error(
    CtwasResult(gwasStudy = "D1", study = "Q1", context = "brain",
                method = "susie", entry = list("not an entry")),
    "must be a CtwasResultEntry")
})

test_that("CtwasResult: constructor rejects mismatched core-vector lengths", {
  expect_error(
    CtwasResult(gwasStudy = c("D1", "D1"), study = "Q1", context = "brain",
                method = "susie", entry = list(.cr_entry())),
    "same length")
})

test_that("CtwasResult: a joint column of the wrong length is rejected", {
  expect_error(
    CtwasResult(gwasStudy = "D1", study = "Q1", context = "brain",
                method = "susie", entry = list(.cr_entry()),
                jointContexts = c("a", "b")),
    "same length as")
})

test_that("CtwasResult: getFinemap/getSusieAlpha are NULL when empty or all-NULL", {
  empty <- CtwasResult(gwasStudy = character(0), study = character(0),
                       context = character(0), method = character(0),
                       entry = list())
  expect_null(getFinemap(empty))
  expect_null(getSusieAlpha(empty))
  nullE <- CtwasResult(gwasStudy = "D1", study = "Q1", context = "brain",
                       method = "susie",
                       entry = list(CtwasResultEntry(finemap = NULL,
                                                     susieAlpha = NULL)))
  expect_null(getFinemap(nullE))
  expect_null(getSusieAlpha(nullE))
})

test_that("CtwasResult: show() prints a one-line-per-run summary", {
  cr <- CtwasResult(gwasStudy = "D1", study = "Q1", context = "brain",
                    method = "susie", entry = list(.cr_entry()))
  expect_output(show(cr), "CtwasResult: 1 run")
  expect_output(show(cr), "susie")
})
