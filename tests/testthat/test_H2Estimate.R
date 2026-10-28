test_that("H2Estimate constructs with all slots", {
    obj <- new(
        "H2Estimate",
        estimate = 0.3,
        estimateSe = 0.05,
        intercept = 1.01,
        interceptSe = 0.02,
        localBlocks = NULL,
        enrichment = NULL,
        annotationJackknifeCoefs = NULL,
        scoreStats = NULL,
        method = "lder",
        nSnps = 10000L,
        traitName = "height"
    )
    expect_s4_class(obj, "H2Estimate")
    expect_equal(heritabilityEstimate(obj), 0.3)
    expect_equal(methodNames(obj), "lder")
})

test_that("show(H2Estimate) does not error", {
    h2 <- new(
        "H2Estimate",
        estimate = 0.3,
        estimateSe = 0.05,
        intercept = 1.0,
        interceptSe = 0.01,
        localBlocks = NULL,
        enrichment = NULL,
        annotationJackknifeCoefs = NULL,
        scoreStats = NULL,
        method = "lder",
        nSnps = 1000L,
        traitName = "test"
    )
    expect_output(show(h2), "H2Estimate")
})

test_that("H2Estimate types its optional slots", {
    # Each optional slot is a type-or-NULL union, so a tibble fills the table
    # slots while a value of the wrong type is refused when the object is
    # built.
    h2 <- function(...) {
        new(
            "H2Estimate",
            estimate = 0.3,
            estimateSe = 0.05,
            intercept = 1.0,
            interceptSe = 0.01,
            method = "lder",
            nSnps = 1000L,
            traitName = "test",
            ...
        )
    }
    local <- tibble::tibble(blockId = 1L, h2Local = 0.1, h2LocalSe = 0.01)
    expect_identical(localH2Blocks(h2(localBlocks = local)), local)
    expect_null(localH2Blocks(h2()))
    expect_error(
        h2(localBlocks = list(blockId = 1L)),
        'invalid object for slot "localBlocks"'
    )
    expect_error(
        h2(enrichment = "none"),
        'invalid object for slot "enrichment"'
    )
    expect_error(
        h2(annotationJackknifeCoefs = c(0.1, 0.2)),
        'invalid object for slot "annotationJackknifeCoefs"'
    )
    expect_error(
        h2(scoreStats = c(z = 1.5)),
        'invalid object for slot "scoreStats"'
    )
})

# ===========================================================================
# H2Estimate accessors
# ===========================================================================

test_that("stratifiedHeritabilityEnrichment returns the enrichment slot", {
    h2_obj <- makeTestH2Estimate(with_enrichment = TRUE)
    result <- stratifiedHeritabilityEnrichment(h2_obj)
    expect_true(is.data.frame(result))
    expect_equal(result$annotation, c("annot1", "annot2"))
})

test_that("scoreStats returns the scoreStats slot", {
    scoreStats <- list(
        z = c(1.5, -0.8),
        R = diag(2),
        annotationNames = c("a1", "a2")
    )
    h2_obj <- new(
        "H2Estimate",
        estimate = 0.3,
        estimateSe = 0.05,
        intercept = 1.0,
        interceptSe = 0.01,
        localBlocks = NULL,
        enrichment = NULL,
        annotationJackknifeCoefs = NULL,
        scoreStats = scoreStats,
        method = "lder",
        nSnps = 100L,
        traitName = "test"
    )
    expect_identical(scoreStats(h2_obj), scoreStats)
})

test_that("accessors return NULL when slots are NULL", {
    h2_obj <- new(
        "H2Estimate",
        estimate = 0.3,
        estimateSe = 0.05,
        intercept = 1.0,
        interceptSe = 0.01,
        localBlocks = NULL,
        enrichment = NULL,
        annotationJackknifeCoefs = NULL,
        scoreStats = NULL,
        method = "lder",
        nSnps = 100L,
        traitName = "test"
    )
    expect_null(localH2Blocks(h2_obj))
    expect_null(stratifiedHeritabilityEnrichment(h2_obj))
    expect_null(scoreStats(h2_obj))
})

test_that("updateObject() moves an H2Estimate's estimates off the old names", {
    h2Obj <- new(
        "H2Estimate",
        estimate = 0.3,
        estimateSe = 0.05,
        intercept = 1.0,
        interceptSe = 0.01,
        localBlocks = data.frame(blockId = "b1", h2Local = 0.1, h2LocalSe = 0),
        enrichment = NULL,
        annotationJackknifeCoefs = matrix(0.2, nrow = 2L, ncol = 1L),
        scoreStats = NULL,
        method = "sldsc",
        nSnps = 100L,
        traitName = "test"
    )
    saved <- .savedUnder(
        h2Obj,
        c(
            h2 = "estimate",
            h2Se = "estimateSe",
            local = "localBlocks",
            tauBlocks = "annotationJackknifeCoefs"
        )
    )
    # The inherited update alone drops slots the class no longer has, with a
    # warning, and leaves the estimate empty.
    expect_no_warning(up <- updateObject(saved))
    expect_identical(up, h2Obj)
    expect_identical(updateObject(h2Obj), h2Obj)
    expect_message(updateObject(h2Obj, verbose = TRUE), "H2Estimate")
})
