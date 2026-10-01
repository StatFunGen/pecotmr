context("filterRelatedness")

# Helper: check no remaining pairs in the kept set exceed threshold
no_related_pairs_remain <- function(
    relatedness,
    excluded,
    threshold,
    iid1 = "IID1",
    iid2 = "IID2",
    value = "PI_HAT"
) {
    kept <- relatedness[
        !(relatedness[[iid1]] %in% excluded) &
            !(relatedness[[iid2]] %in% excluded),
    ]
    all(kept[[value]] < threshold)
}

test_that("maximizeUnrelated removes related individuals and leaves clean set", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("plinkQC")

    # 10 individuals, several pairs above threshold 0.125
    rel <- data.frame(
        IID1 = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "A"),
        IID2 = c("B", "C", "D", "E", "F", "G", "H", "I", "J", "J"),
        PI_HAT = c(0.25, 0.15, 0.30, 0.08, 0.20, 0.05, 0.18, 0.03, 0.22, 0.14),
        stringsAsFactors = FALSE
    )
    threshold <- 0.125

    result <- filterRelatedness(
        relatedness = rel,
        relatednessThreshold = threshold,
        analysisType = "maximizeUnrelated"
    )

    expect_type(result, "character")
    expect_true(length(result) > 0)
    expect_true(no_related_pairs_remain(rel, result, threshold))
})

test_that("no related pairs returns empty exclusion vector", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("plinkQC")

    rel <- data.frame(
        IID1 = c("A", "B", "C", "D"),
        IID2 = c("B", "C", "D", "E"),
        PI_HAT = c(0.01, 0.02, 0.03, 0.04),
        stringsAsFactors = FALSE
    )
    threshold <- 0.125

    result <- filterRelatedness(
        relatedness = rel,
        relatednessThreshold = threshold,
        analysisType = "maximizeUnrelated"
    )

    expect_type(result, "character")
    expect_equal(length(result), 0)
})

test_that("large component pre-pruning removes individuals", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("plinkQC")

    # Build a chain of 30 individuals: 1-2, 2-3, ..., 29-30, all above threshold
    n <- 30
    ids <- paste0("IND", seq_len(n))
    rel <- data.frame(
        IID1 = ids[seq_len(n - 1)],
        IID2 = ids[2:n],
        PI_HAT = rep(0.20, n - 1),
        stringsAsFactors = FALSE
    )
    threshold <- 0.125

    # verbose = TRUE exercises the graph pre-pruning progress messages and the
    # final exclusion-count message.
    result <- filterRelatedness(
        relatedness = rel,
        relatednessThreshold = threshold,
        analysisType = "maximizeUnrelated",
        pruning = relatednessPruning(maxComponentSize = 10),
        verbose = TRUE
    )

    expect_type(result, "character")
    expect_true(length(result) > 0)
    # The remaining kept individuals should have no related pairs
    expect_true(no_related_pairs_remain(rel, result, threshold))
})

test_that("maximizeCases preferentially retains cases", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("plinkQC")

    # Cases: C1, C2, C3; Controls: X1, X2, X3
    # Related pairs above threshold:
    #   C1-X1 (case-control), C2-X2 (case-control), C1-C2 (case-case),
    #   X1-X3 (control-control)
    rel <- data.frame(
        IID1 = c("C1", "C2", "C1", "X1", "C3", "X2"),
        IID2 = c("X1", "X2", "C2", "X3", "X3", "X3"),
        PI_HAT = c(0.25, 0.20, 0.15, 0.18, 0.05, 0.04),
        stringsAsFactors = FALSE
    )

    pheno <- data.frame(
        IID = c("C1", "C2", "C3", "X1", "X2", "X3"),
        pheno = c(1, 1, 1, 0, 0, 0),
        stringsAsFactors = FALSE
    )

    threshold <- 0.125
    result <- filterRelatedness(
        relatedness = rel,
        relatednessThreshold = threshold,
        analysisType = "maximizeCases",
        phenoData = pheno,
        phenoCol = "pheno"
    )

    expect_type(result, "character")
    # Controls related to kept cases should be excluded preferentially
    retained <- setdiff(pheno$IID, result)
    retained_cases <- intersect(retained, pheno$IID[pheno$pheno == 1])
    retained_controls <- intersect(retained, pheno$IID[pheno$pheno == 0])
    # We expect at least 2 of 3 cases kept, and controls sacrificed
    expect_true(length(retained_cases) >= 2)
    expect_true(no_related_pairs_remain(rel, result, threshold))
})

test_that("maximizeCases errors without phenoData", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("plinkQC")

    rel <- data.frame(
        IID1 = c("A", "B"),
        IID2 = c("B", "C"),
        PI_HAT = c(0.25, 0.20),
        stringsAsFactors = FALSE
    )

    expect_error(
        filterRelatedness(
            relatedness = rel,
            relatednessThreshold = 0.125,
            analysisType = "maximizeCases"
        ),
        "Must provide phenoData"
    )
})

test_that("maximizeCases excludes a control listed as IID1 paired with a kept case", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("plinkQC")

    # X1-C1 is a control(IID1)-case(IID2) pair, exercising the mirror branch of
    # the case/control exclusion loop. C1-C2 keeps both cases (sub-threshold);
    # X2-X3 gives the control-control step a non-empty input.
    rel <- data.frame(
        IID1 = c("X1", "C1", "X2"),
        IID2 = c("C1", "C2", "X3"),
        PI_HAT = c(0.25, 0.05, 0.20),
        stringsAsFactors = FALSE
    )
    pheno <- data.frame(
        IID = c("C1", "C2", "X1", "X2", "X3"),
        pheno = c(1, 1, 0, 0, 0),
        stringsAsFactors = FALSE
    )

    result <- filterRelatedness(
        relatedness = rel,
        relatednessThreshold = 0.125,
        analysisType = "maximizeCases",
        phenoData = pheno
    )

    expect_type(result, "character")
    # X1 is a control related to the retained case C1, so it must be excluded.
    expect_true("X1" %in% result)
})

test_that("iterative cleanup loops and warns when related pairs persist", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("plinkQC")

    # Force plinkQC to exclude nobody, so related pairs survive Phase 2 and the
    # Phase-3 iterative cleanup loop runs to exhaustion, emitting the warning.
    local_mocked_bindings(
        relatednessFilter = function(...) {
            list(
                failIDs = data.frame(
                    IID = character(0),
                    stringsAsFactors = FALSE
                )
            )
        },
        .package = "plinkQC"
    )

    rel <- data.frame(
        IID1 = c("A", "B", "C"),
        IID2 = c("B", "C", "D"),
        PI_HAT = c(0.30, 0.30, 0.30),
        stringsAsFactors = FALSE
    )

    expect_warning(
        result <- filterRelatedness(
            relatedness = rel,
            relatednessThreshold = 0.125,
            analysisType = "maximizeUnrelated",
            pruning = relatednessPruning(maxIterations = 2L),
            verbose = TRUE
        ),
        "related pairs remain"
    )
    # Nobody is excluded because the (mocked) filter never fails anyone.
    expect_type(result, "character")
    expect_equal(length(result), 0L)
})

test_that(".relatednessLargestComponent guards an empty component list", {
    # max(integer(0)) warns and returns -Inf; with no related pairs the
    # pruning loop must simply not run, without emitting that warning.
    expect_identical(
        pecotmr:::.relatednessLargestComponent(list(csize = integer(0))),
        0L
    )
    expect_identical(
        pecotmr:::.relatednessLargestComponent(list(csize = c(3L, 7L))),
        7L
    )
})

test_that("relatednessColumns carries the defaults and rejects a typo", {
    cols <- relatednessColumns()
    expect_s4_class(cols, "MethodConfig")
    expect_equal(cols$iid1, "IID1")
    expect_equal(cols$value, "PI_HAT")
    # No `...`, so R's own argument matching is the check.
    expect_error(relatednessColumns(vlaue = "X"), "unused argument")
})

test_that("relatednessColumns drops NULL family-ID columns", {
    cols <- relatednessColumns()
    # A NULL default means "the table has none", so it is not carried as an
    # explicit NULL into plinkQC's argument list.
    expect_false(is_in("fid1", names(cols)))
    expect_null(cols$fid1)
    expect_true(is_in("fid1", names(relatednessColumns(fid1 = "FID1"))))
})

test_that("relatednessPruning carries the defaults and rejects a typo", {
    pr <- relatednessPruning(maxComponentSize = 50L)
    expect_s4_class(pr, "MethodConfig")
    expect_equal(pr$maxComponentSize, 50L)
    expect_equal(pr$maxIterations, 20L)
    expect_error(relatednessPruning(maxComponents = 5L), "unused argument")
})

test_that("plinkQcConfig validates extras against plinkQC's live formals", {
    skip_if_not_installed("plinkQC")
    pq <- plinkQcConfig(otherCriterionThDirection = "le")
    expect_s4_class(pq, "MethodConfig")
    expect_equal(pq$otherCriterionThDirection, "le")
    # The check comes from names(formals(plinkQC::relatednessFilter)), so it
    # cannot drift from the upstream signature.
    expect_error(
        plinkQcConfig(otherCriterionThreshold = 1),
        "unknown argument"
    )
})

test_that("plinkQcConfig refuses the arguments filterRelatedness derives", {
    skip_if_not_installed("plinkQC")
    # relatednessTh is legal for plinkQC, so the constructor accepts it, but
    # filterRelatedness injects its own and the duplicate is caught there.
    expect_error(
        filterRelatedness(
            data.frame(IID1 = "a", IID2 = "b", PI_HAT = 0.5),
            plinkQcArgs = plinkQcConfig(relatednessTh = 0.5)
        ),
        "derived from `relatednessThreshold` and `columns`"
    )
})

test_that("filterRelatedness refuses a bare list where a constructor is due", {
    rel <- data.frame(IID1 = "a", IID2 = "b", PI_HAT = 0.5)
    expect_error(
        filterRelatedness(rel, columns = list(iid1 = "IID1")),
        "must be built with relatednessColumns\\(\\)"
    )
    expect_error(
        filterRelatedness(rel, pruning = list(maxIterations = 2L)),
        "must be built with relatednessPruning\\(\\)"
    )
})
