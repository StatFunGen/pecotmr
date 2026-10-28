# Tests for R/AllGenerics.R, which holds only setGeneric() declarations -- no
# behaviour of its own. What is checked here is the property that file is
# responsible for: that the generic surface and the method surface agree.
#
# A generic with no methods is dead API. It exports a name, generates a help
# page, and answers every call with "unable to find an inherited method", which
# reads to a user like a bug in their input rather than a gap in the package.
# The refactor produced one of these (`molecularTraits` briefly had no method at
# all), so this is a live failure mode, not a hypothetical.

# @noRd
.agHasMethod <- function(nm) {
    tryCatch(
        length(methods::findMethods(nm, where = asNamespace("pecotmr"))) > 0L,
        error = function(e) FALSE
    )
}

test_that("every generic declared in AllGenerics.R has at least one method", {
    src <- readLines(
        system.file(
            "..",
            "R",
            "AllGenerics.R",
            package = "pecotmr",
            mustWork = FALSE
        ),
        warn = FALSE
    )
    skip_if(length(src) == 0L, "R/AllGenerics.R not visible from this run")
    declared <- str_match(src, '^\\s*setGeneric\\(\\s*"([^"]+)"')[, 2]
    declared <- declared[!is.na(declared)]
    # Multi-line form: setGeneric(\n    "name", ...
    contd <- which(str_detect(src, "^\\s*setGeneric\\($"))
    extra <- str_match(src[contd + 1L], '^\\s*"([^"]+)"')[, 2]
    declared <- unique(c(declared, extra[!is.na(extra)]))
    expect_gt(length(declared), 50L)

    orphans <- declared[!map_lgl(declared, .agHasMethod)]
    expect_equal(
        orphans,
        character(0),
        label = "generics with no method"
    )
})
