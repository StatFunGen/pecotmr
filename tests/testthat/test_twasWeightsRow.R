# Tests for R/TwasWeightsRow.R

# === Tests migrated from test_s4Constructors.R (TwasWeightsRow) ===

test_that("TwasWeightsRow: constructor and accessors round-trip", {
    e <- twasWeightsRow(
        variantIds = c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G"),
        weights = c(0.1, -0.2, 0.05),
        fits = list(model = "lasso"),
        cvResult = list(rsq = 0.4),
        standardized = TRUE,
        dataType = "expression"
    )
    expect_s4_class(e, "TwasWeightsRow")
    expect_equal(
        .twrPartsVariantIds(e),
        c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G")
    )
    expect_equal(getWeights(e), c(0.1, -0.2, 0.05))
    expect_equal(getFits(e), list(model = "lasso"))
    expect_equal(getCvResult(e), list(rsq = 0.4))
    expect_true(isTRUE(getStandardized(e)))
    expect_equal(getDataType(e), "expression")
})

test_that("TwasWeightsRow: resolveWeights returns the aligned (variantIds, weights) pair", {
    e <- twasWeightsRow(
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
    e2 <- twasWeightsRow(
        variantIds = c("chr1:100:A:G", "chr1:200:A:G"),
        weights = c(0.1, -0.2)
    )
    e2@weights <- matrix(0, 2, 2) # flattens to length 4 != 2 ids
    expect_length(.twrRowResolveWeights(e2)$variantIds, 0L)
})


test_that("TwasWeightsRow: standardized is coerced via isTRUE() semantics", {
    # isTRUE() only returns TRUE for a length-1 logical TRUE. Non-TRUE
    # input lands as FALSE (the safe default for the standardized flag).
    e_logical <- twasWeightsRow(
        variantIds = "chr1:100:A:G",
        weights = 0.1,
        standardized = TRUE
    )
    expect_true(isTRUE(getStandardized(e_logical)))

    e_default <- twasWeightsRow(
        variantIds = "chr1:100:A:G",
        weights = 0.1,
        standardized = "yes-please"
    )
    expect_false(isTRUE(getStandardized(e_default)))
})


test_that("TwasWeightsRow: validity rejects matrix weights with wrong nrow", {
    expect_error(
        twasWeightsRow(
            variantIds = c("chr1:100:A:G", "chr1:200:A:G"),
            weights = matrix(0, nrow = 5, ncol = 1)
        ),
        "nrow\\(weights\\) must equal length\\(variantIds\\)"
    )
})

# ===========================================================================
# QtlFineMappingResult
# ===========================================================================

test_that("TwasWeights: rejects non-TwasWeightsRow rows", {
    expect_error(
        TwasWeights(
            study = "s1",
            context = "c1",
            trait = "t1",
            method = "lasso",
            entry = list("not_an_entry")
        ),
        "must be a TWAS-weight row"
    )
})


test_that("TwasWeights: validity does not recurse on key subset (#546)", {
    # Building the tuple-uniqueness key via `object[, keyCols]` preserves the
    # TwasWeights class while dropping the required `entry` column; older
    # S4Vectors revalidates that intermediate and fails with "missing columns".
    e <- .sc_makeTwasWeightsRow(p = 5L)
    res <- TwasWeights(
        study = "s1",
        context = "c1",
        trait = "t1",
        method = "lasso",
        entry = list(e)
    )
    expect_s4_class(res, "TwasWeights")
    expect_true(validObject(res))
})


# === Tests migrated from test_showMethods.R (TwasWeightsRow) ===

test_that("show.TwasWeightsRow reports standardized flag and CV availability", {
    e <- .sh_makeTwEntry(p = 5, standardized = TRUE)
    out <- capture.output(show(e))
    expect_true(any(grepl(
        "TwasWeightsRow: 5 variants.*standardized=TRUE",
        out
    )))
    expect_true(any(grepl("CV performance: TRUE", out)))

    e_no_cv <- twasWeightsRow(
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

test_that("variantIds that encode no coordinates are rejected", {
    # An id is a rendering of the variant's range and alleles; "v1" names
    # nothing, and the failure should bite here rather than later when the
    # entry is turned into a ranged element.
    expect_error(
        twasWeightsRow(
            variantIds = c("chr1:100:A:G", "v1", "v2"),
            weights = c(0.1, -0.2, 0.05),
            fits = list(model = "lasso"),
            cvResult = list(rsq = 0.4),
            standardized = TRUE,
            dataType = "expression"
        ),
        "do not encode coordinates"
    )
})

test_that("the identity message names the offenders", {
    err <- tryCatch(
        twasWeightsRow(
            variantIds = c("rsX", "rsY"),
            weights = c(0.1, -0.2),
            fits = list(),
            cvResult = list(),
            standardized = TRUE,
            dataType = "expression"
        ),
        error = function(e) conditionMessage(e)
    )
    expect_match(err, "2 of 2")
    expect_match(err, "rsX, rsY")
})

test_that("an empty variant set is allowed", {
    # Emptiness is a separate concern from malformedness: there is nothing to
    # parse, so the coordinate rule has nothing to say.
    e <- twasWeightsRow(
        variantIds = character(0),
        weights = numeric(0),
        fits = list(),
        cvResult = list(),
        standardized = TRUE,
        dataType = "expression"
    )
    expect_s4_class(e, "TwasWeightsRow")
    expect_equal(length(.twrPartsVariantIds(e)), 0L)
})
