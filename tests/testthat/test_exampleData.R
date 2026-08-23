# Tests for R/exampleData.R, which is pure roxygen documenting the datasets
# shipped in data/. There is no code to unit-test, so what is checked here is
# the contract that file asserts: every documented dataset exists, loads, and
# still holds an object of the class its documentation claims.
#
# This is not a formality. Twice during the S4 collection refactor a `.rda`
# went stale -- it still loaded, but carried the OLD internal layout, so the
# object was an instance of a class whose definition had moved on. Nothing
# caught it until an unrelated test failed several steps later. A load-and-
# validate sweep turns that into an immediate, obvious failure.

# The datasets R/exampleData.R documents, paired with the class each one is
# documented to be. `NA` means the documentation does not pin a class (a plain
# list or data.frame fixture), so only the load is checked.
.ed_expected <- list(
    colocboostResultExample = "ColocBoostResult",
    ctwasEstExample = NA_character_,
    ctwasFinemapExample = NA_character_,
    ctwasInputsExample = NA_character_,
    ctwasWeightsExample = "TwasWeights",
    eqtlRegionExample = NA_character_,
    fsusieFineMappingExample = "QtlFineMappingResult",
    gwasFineMappingExample = "GwasFineMappingResult",
    gwasFineMappingLbfExample = "GwasFineMappingResult",
    gwasSumStatsExample = NA_character_,
    gwasSumStatsS4Example = "GwasSumStats",
    h2EstimateExample = "H2Estimate",
    ldEigenExample = "LdEigen",
    ldScoreExample = "LdScore",
    mashInputExample = NA_character_,
    mashPosteriorExample = NA_character_,
    multiStudyQtlDatasetExample = "MultiStudyQtlDataset",
    multiTraitData = NA_character_,
    mvsusieFineMappingExample = "QtlFineMappingResult",
    qtlDatasetExample = "QtlDataset",
    qtlFineMappingExample = "QtlFineMappingResult",
    qtlFineMappingLbfExample = "QtlFineMappingResult",
    qtlFineMappingPairedExample = NA_character_,
    qtlSumStatsExample = "QtlSumStats",
    qtlSumStatsMulticontextExample = "QtlSumStats",
    twasWeightsExample = "TwasWeights"
)

test_that("every documented example dataset ships in data/", {
    shipped <- sub(
        "\\.rda$",
        "",
        basename(
            list.files(
                system.file(
                    "..",
                    "data",
                    package = "pecotmr",
                    mustWork = FALSE
                ),
                pattern = "\\.rda$"
            )
        )
    )
    skip_if(length(shipped) == 0L, "data/ not visible from this test run")
    expect_setequal(shipped, names(.ed_expected))
})

test_that("every documented example dataset loads", {
    for (nm in names(.ed_expected)) {
        env <- new.env()
        loaded <- tryCatch(
            {
                data(list = nm, envir = env)
                nm %in% ls(env)
            },
            error = function(e) FALSE,
            warning = function(w) FALSE
        )
        expect_true(loaded, label = str_c("dataset '", nm, "' loads"))
    }
})

test_that("example datasets still hold the class they document", {
    # A stale .rda loads happily but carries the old layout; validObject()
    # is what distinguishes "loaded" from "still valid".
    for (nm in names(.ed_expected)) {
        cls <- .ed_expected[[nm]]
        if (is.na(cls)) {
            next
        }
        env <- new.env()
        data(list = nm, envir = env)
        obj <- get(nm, envir = env)
        expect_s4_class(obj, cls)
        expect_true(
            validObject(obj, test = TRUE) == TRUE,
            label = str_c("dataset '", nm, "' is a valid ", cls)
        )
    }
})

test_that("S4 example collections are non-empty and self-consistent", {
    # nrow() is the collection's element count; an empty example is almost
    # always a build accident rather than an intended fixture.
    collections <- c(
        "qtlFineMappingExample",
        "gwasFineMappingExample",
        "qtlSumStatsExample",
        "gwasSumStatsS4Example",
        "twasWeightsExample"
    )
    for (nm in collections) {
        env <- new.env()
        data(list = nm, envir = env)
        obj <- get(nm, envir = env)
        expect_gt(nrow(obj), 0L)
        expect_equal(length(colnames(obj)), ncol(obj))
    }
})
