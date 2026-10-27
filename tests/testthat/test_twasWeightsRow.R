# Tests for R/twasWeightsRow.R

# === Tests migrated from test_s4Constructors.R (TwasWeightsRow) ===

test_that("as(row, 'GRanges') carries the weights", {
    row <- TwasWeightsRow(
        variantIds = c("chr1:100:A:G", "chr1:200:C:T"),
        weights = c(0.4, -0.2)
    )
    gr <- as(row, "GRanges")
    expect_s4_class(gr, "GRanges")
    expect_identical(gr$weight, c(0.4, -0.2))
})

test_that("TwasWeightsRow: constructor and accessors round-trip", {
    e <- TwasWeightsRow(
        variantIds = c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G"),
        weights = c(0.1, -0.2, 0.05),
        methodFits = list(model = "lasso"),
        cvResult = list(rsq = 0.4),
        weightStandardized = TRUE,
        weightsDataType = "expression"
    )
    expect_s4_class(e, "TwasWeightsRow")
    expect_equal(
        .twrPartsVariantIds(e),
        c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G")
    )
    expect_equal(weights(e), c(0.1, -0.2, 0.05))
    expect_equal(methodFits(e), list(model = "lasso"))
    expect_equal(cvResult(e), list(rsq = 0.4))
    expect_true(isTRUE(weightStandardized(e)))
    expect_equal(weightsDataType(e), "expression")
})

test_that("TwasWeightsRow: resolveWeights returns the aligned (variantIds, weights) pair", {
    e <- TwasWeightsRow(
        variantIds = c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G"),
        weights = c(0.1, -0.2, 0.05)
    )
    wr <- .twrRowResolveWeights(e)
    expect_equal(
        wr$variantIds,
        c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G")
    )
    expect_equal(wr$weights, c(0.1, -0.2, 0.05))
    # length mismatch (defensive) -> empty
    e2 <- TwasWeightsRow(
        variantIds = c("chr1:100:A:G", "chr1:200:A:G"),
        weights = c(0.1, -0.2)
    )
    e2@weights <- matrix(0, 2, 2) # flattens to length 4 != 2 ids
    expect_length(.twrRowResolveWeights(e2)$variantIds, 0L)
})


test_that("TwasWeightsRow: standardized is coerced via isTRUE() semantics", {
    # isTRUE() only returns TRUE for a length-1 logical TRUE. Non-TRUE
    # input lands as FALSE (the safe default for the standardized flag).
    e_logical <- TwasWeightsRow(
        variantIds = "chr1:100:A:G",
        weights = 0.1,
        weightStandardized = TRUE
    )
    expect_true(isTRUE(weightStandardized(e_logical)))

    e_default <- TwasWeightsRow(
        variantIds = "chr1:100:A:G",
        weights = 0.1,
        weightStandardized = "yes-please"
    )
    expect_false(isTRUE(weightStandardized(e_default)))
})


test_that("TwasWeightsRow: validity rejects matrix weights with wrong nrow", {
    expect_error(
        TwasWeightsRow(
            variantIds = c("chr1:100:A:G", "chr1:200:A:G"),
            weights = matrix(0, nrow = 5, ncol = 1)
        ),
        "nrow\\(weights\\) must equal length\\(variantIds\\)"
    )
})

test_that("TwasWeightsRow: rejects mismatched variantIds/weights lengths", {
    # The weights are per-variant mcols on the entry's element now, so a
    # length mismatch cannot be represented. It is refused at construction
    # rather than producing a gene that gets silently dropped downstream.
    expect_error(
        TwasWeightsRow(
            variantIds = sprintf("chr1:%d:G:A", 100L * 1:5),
            weights = c(0.1, 0.2, 0.3)
        ),
        "length\\(weights\\) is 3 but 5 variants were supplied"
    )
})

# ===========================================================================
# TwasWeights
# ===========================================================================

# === Tests migrated from test_showMethods.R (TwasWeightsRow) ===

test_that("show.TwasWeightsRow reports standardized flag and CV availability", {
    e <- .sh_makeTwEntry(p = 5, weightStandardized = TRUE)
    out <- capture.output(show(e))
    expect_true(any(grepl(
        "TwasWeightsRow: 5 variants.*standardized=TRUE",
        out
    )))
    expect_true(any(grepl("CV performance: TRUE", out)))

    e_no_cv <- TwasWeightsRow(
        variantIds = c("chr1:100:A:G", "chr1:200:A:G"),
        weights = c(0.1, 0.2)
    )
    out2 <- capture.output(show(e_no_cv))
    expect_true(any(grepl("CV performance: FALSE", out2)))
})

# ===========================================================================
# Variant-id identity check
#
# Every fixture above uses well-formed chrom:pos:ref:alt ids, so the branch
# that rejects ids carrying no coordinates was never executed.
# ===========================================================================

test_that("variantIds that encode no coordinates are rejected by name", {
    # An id is a rendering of the variant's range and alleles; "v1" names
    # nothing, and the failure should bite here rather than later when the
    # entry is turned into a ranged element.
    err <- tryCatch(
        TwasWeightsRow(
            variantIds = c("chr1:100:A:G", "v1", "v2"),
            weights = c(0.1, -0.2, 0.05),
            methodFits = list(model = "lasso"),
            cvResult = list(rsq = 0.4),
            weightStandardized = TRUE,
            weightsDataType = "expression"
        ),
        error = function(e) conditionMessage(e)
    )
    expect_match(err, "do not encode coordinates")
    # The message counts and names the offenders.
    expect_match(err, "2 of 3")
    expect_match(err, "v1, v2")
})

test_that("an empty variant set is allowed", {
    # Emptiness is a separate concern from malformedness: there is nothing to
    # parse, so the coordinate rule has nothing to say.
    e <- TwasWeightsRow(
        variantIds = character(0),
        weights = numeric(0),
        methodFits = list(),
        cvResult = list(),
        weightStandardized = TRUE,
        weightsDataType = "expression"
    )
    expect_s4_class(e, "TwasWeightsRow")
    expect_equal(length(.twrPartsVariantIds(e)), 0L)
})
