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
    # Named by the row's variants, whose ids the row holds once; the stored
    # weights carry no copy of them.
    expect_equal(
        weights(e),
        set_names(c(0.1, -0.2, 0.05), .twrPartsVariantIds(e))
    )
    expect_null(names(S4Vectors::mcols(e@variants)$weight))
    expect_equal(methodFits(e), list(model = "lasso"))
    expect_equal(cvResult(e), list(rsq = 0.4))
    expect_true(isTRUE(weightStandardized(e)))
    expect_equal(dataType(e), "expression")
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

# ===========================================================================
# The variants a stored fit was made on
# ===========================================================================

.twr_ids <- c("chr1:100:A:G", "chr1:200:C:T")

test_that("a fit's variants are read where its engine names them", {
    f <- pecotmr:::.twasFitVariantIds
    named <- set_names(c(0.1, 0.2), .twr_ids)
    path <- matrix(0, 2L, 3L, dimnames = list(.twr_ids, NULL))
    geno <- matrix(0, 3L, 2L, dimnames = list(NULL, .twr_ids))
    expect_equal(f(list(coef = c(0.1, 0.2), variantIds = .twr_ids)), .twr_ids)
    # glmnet's path, ncvreg's past its intercept, qgg's and BGLR's effects.
    expect_equal(f(list(glmnet.fit = list(beta = path))), .twr_ids)
    intercept <- matrix(0, 1L, 3L, dimnames = list("(Intercept)", NULL))
    expect_equal(f(list(fit = list(beta = rbind(intercept, path)))), .twr_ids)
    expect_equal(f(list(bm = named)), .twr_ids)
    expect_equal(f(list(ETA = list(list(b = named)))), .twr_ids)
    # mr.ash's and DPR's copies of the genotypes, a SuSiE fit's names.
    expect_equal(f(list(data = list(X = geno))), .twr_ids)
    expect_equal(f(list(x = geno)), .twr_ids)
    expect_equal(f(list(pip = named)), .twr_ids)
    # A row joined from regions: each region's variants in turn.
    regions <- list(
        r1 = list(variantIds = .twr_ids[[1L]]),
        r2 = list(variantIds = .twr_ids[[2L]])
    )
    expect_equal(f(regions), .twr_ids)
    # Nothing names them, or there is no fit.
    expect_equal(f(list(coef = c(0.1, 0.2))), character(0))
    expect_equal(f(NULL), character(0))
})

test_that("a kept fit takes the variant ids only when it names none", {
    g <- pecotmr:::.withFitVariantIds
    expect_equal(g(list(coef = c(0.1, 0.2)), .twr_ids)$variantIds, .twr_ids)
    # One that already names them keeps its own names alone.
    named <- list(pip = set_names(c(0.1, 0.9), .twr_ids))
    expect_identical(g(named, rev(.twr_ids)), named)
    expect_null(g(NULL, .twr_ids))
    expect_identical(g(list(coef = 1), NULL), list(coef = 1))
    # An engine's class survives the addition.
    classed <- structure(list(coef = c(0.1, 0.2)), class = "someEngineFit")
    expect_s3_class(g(classed, .twr_ids), "someEngineFit")
    # On weights, only a fit they carry is touched.
    w <- `attr<-`(c(0.1, 0.2), "fit", list(coef = c(0.1, 0.2)))
    expect_equal(
        attr(pecotmr:::.withWeightsFitIds(w, .twr_ids), "fit")$variantIds,
        .twr_ids
    )
    expect_identical(
        pecotmr:::.withWeightsFitIds(c(0.1, 0.2), .twr_ids),
        c(0.1, 0.2)
    )
})

test_that("updateObject() drops a row's saved fitVariants slot", {
    row <- TwasWeightsRow(.twr_ids, c(0.1, 0.2))
    saved <- `attr<-`(row, "fitVariants", .twr_ids)
    up <- updateObject(saved)
    expect_null(attr(up, "fitVariants", exact = TRUE))
    expect_identical(up, row)
    expect_identical(updateObject(row), row)
    expect_message(updateObject(row, verbose = TRUE), "TwasWeightsRow")
})
