context("CtwasResultEntry")

# Migrated from test_CtwasResult.R: the per-run cTWAS payload class (constructor
# + ctwasPosteriors / susieAlpha / ctwasGroupPriors accessors). The CtwasResult
# collection tests live in test_CtwasResult.R.

test_that("CtwasResultEntry: constructor + accessors round-trip", {
    e <- CtwasResultEntry(
        posteriors = data.frame(
            id = c("g1", "g2"),
            susie_pip = c(0.9, 0.1),
            stringsAsFactors = FALSE
        ),
        susieAlpha = data.frame(
            id = c("g1", "g2"),
            susie_alpha = c(0.9, 0.1),
            stringsAsFactors = FALSE
        ),
        groupPriors = list(group_prior = c(brain = 0.01)),
        regionInfo = data.frame(region_id = "r1", stringsAsFactors = FALSE)
    )
    expect_s4_class(e, "CtwasResultEntry")
    expect_equal(nrow(ctwasPosteriors(e)), 2L)
    expect_equal(nrow(susieAlpha(e)), 2L)
    expect_equal(ctwasGroupPriors(e)$group_prior, c(brain = 0.01))
})

test_that("CtwasResultEntry: an empty entry is valid and its accessors return NULL", {
    e <- CtwasResultEntry()
    expect_s4_class(e, "CtwasResultEntry")
    expect_null(ctwasPosteriors(e))
    expect_null(susieAlpha(e))
    expect_null(ctwasGroupPriors(e))
})
