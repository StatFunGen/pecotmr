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

test_that("as(entry, 'GRanges') gives placed posteriors only", {
    placed <- GenomicRanges::GRanges("chr1:100-200", susie_pip = 0.4)
    expect_identical(
        as(CtwasResultEntry(posteriors = placed), "GRanges"),
        placed
    )
    expect_length(as(CtwasResultEntry(), "GRanges"), 0L)
    expect_error(
        as(CtwasResultEntry(posteriors = data.frame(id = "g1")), "GRanges"),
        "not placed on the genome"
    )
})

test_that("CtwasResultEntry: an empty entry is valid and its accessors return NULL", {
    e <- CtwasResultEntry()
    expect_s4_class(e, "CtwasResultEntry")
    expect_null(ctwasPosteriors(e))
    expect_null(susieAlpha(e))
    expect_null(ctwasGroupPriors(e))
})

test_that("show() counts the genes and regions an entry scores", {
    e <- CtwasResultEntry(
        posteriors = data.frame(
            id = c("g1", "g2", "rs1", "g1"),
            type = c("gene", "gene", "SNP", "gene"),
            region_id = c("r1", "r1", "r1", "r2"),
            susie_pip = c(0.9, 0.1, 0.2, 0.3)
        ),
        groupPriors = list(group_prior = c(brain = 0.01))
    )
    expect_output(
        show(e),
        "CtwasResultEntry: 2 gene\\(s\\) in 2 region\\(s\\); group priors set"
    )
})

test_that("show() reads a ranged entry, and an empty one", {
    placed <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(100L, 200L), width = 1L),
        id = c("g1", "g2"),
        region_id = c("r1", "r1"),
        susie_pip = c(0.9, 0.1)
    )
    expect_output(
        show(CtwasResultEntry(posteriors = placed)),
        "2 gene\\(s\\) in 1 region\\(s\\); group priors not set"
    )
    # No id or region column to count, and no posteriors at all.
    expect_output(
        show(CtwasResultEntry(posteriors = data.frame(susie_pip = 0.5))),
        "0 gene\\(s\\) in 0 region\\(s\\)"
    )
    expect_output(show(CtwasResultEntry()), "0 gene\\(s\\) in 0 region\\(s\\)")
})
