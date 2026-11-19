# One estimate covering `window` (none for a genome-wide estimate), under the
# settings `param`, with `...` replacing any of its columns.
.h2Est <- function(
    window = GenomicRanges::GRanges(),
    param = .h2GenomeWideSettings("lder", FALSE, inSample = FALSE),
    ...
) {
    cols <- list_assign(
        list(
            study = "s1",
            context = NA_character_,
            trait = NA_character_,
            h2 = 0.3,
            h2_se = 0.05,
            h2_lower = NA_real_,
            h2_upper = NA_real_,
            h2_liability = NA_real_,
            h2_liability_se = NA_real_,
            h2_liability_lower = NA_real_,
            h2_liability_upper = NA_real_,
            intercept = 1,
            intercept_se = 0.01,
            p_value = NA_real_,
            n_snps = 100L,
            n_components = NA_integer_,
            scale = NA_real_,
            reference_overlap = NA_real_,
            enrichment = list(NULL),
            jackknife = list(NULL),
            score_stats = list(NULL)
        ),
        ...
    )
    .h2EstimateNew(
        list(window),
        cols[setdiff(names(cols), .H2_OBJECT_COLS)],
        cols[.H2_OBJECT_COLS],
        GenomicRanges::GRanges(),
        param
    )
}

test_that("H2Estimate holds each estimate's values in mcols", {
    x <- .h2Est()
    expect_s4_class(x, "H2Estimate")
    expect_identical(length(x), 1L)
    # A genome-wide estimate's element holds no ranges.
    expect_length(x[[1]], 0L)
    expect_equal(heritabilityEstimate(x), 0.3)
    expect_equal(heritabilityEstimateSe(x), 0.05)
    expect_equal(heritabilityIntercept(x), 1)
    expect_equal(heritabilityInterceptSe(x), 0.01)
    expect_identical(nSnps(x), 100L)
    expect_identical(studyName(x), "s1")
    # A GWAS estimate has no trait.
    expect_identical(traitNames(x), character(0))
    expect_output(show(x), "H2Estimate: 1 estimate\\(s\\)")
    expect_output(show(x[0]), "0 estimate")
    # A cis estimate's element is the window it covers.
    window <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100L, 900L))
    cis <- .h2Est(window, H2EstimateParam(), trait = "ENSG1")
    expect_identical(IRanges::ranges(cis[[1]]), IRanges::ranges(window))
    expect_identical(traitNames(cis), "ENSG1")
})

test_that("H2Estimate records the settings its estimates were fitted under", {
    x <- .h2Est()
    expect_identical(
        h2EstimateParam(x),
        .h2GenomeWideSettings("lder", FALSE, inSample = FALSE)
    )
    expect_output(
        show(x),
        "settings: eigen / regression / floored noise / reference weightsOnly"
    )
    # A subset keeps them, and so does an empty collection.
    expect_identical(h2EstimateParam(x[0]), h2EstimateParam(x))
    expect_identical(h2EstimateParam(new("H2Estimate")), H2EstimateParam())
})

test_that("H2Estimate refuses a missing or mistyped column", {
    x <- .h2Est()
    md <- S4Vectors::mcols(x)
    expect_error(
        validObject(S4Vectors::`mcols<-`(
            x,
            value = md[-match("scale", colnames(md))]
        )),
        "scale"
    )
    expect_error(.h2Est(h2 = "high"), "`h2` column must be double")
    expect_error(.h2Est(n_snps = 100), "`n_snps` column must be integer")
    bad <- md
    bad$enrichment <- "none"
    expect_error(
        validObject(S4Vectors::`mcols<-`(x, value = bad)),
        "`enrichment` column must be a list"
    )
})

test_that("H2Estimate holds one row per estimate, under one set of settings", {
    x <- .h2Est()
    expect_error(c(x, x), "tuple uniqueness violated")
    # Another study under the same settings is a second estimate.
    both <- c(x, .h2Est(study = "s2", h2 = 0.2))
    expect_identical(length(both), 2L)
    expect_equal(heritabilityEstimate(both), c(0.3, 0.2))
    expect_identical(h2EstimateParam(both), h2EstimateParam(x))
    expect_output(show(both), "2 estimate\\(s\\)")
    # One made under other settings cannot join it.
    hdl <- .h2Est(
        param = .h2GenomeWideSettings("hdl", FALSE, inSample = FALSE),
        study = "s2"
    )
    expect_error(c(x, hdl), "fitted under different settings")
    expect_error(append(x, hdl), "fitted under different settings")
})

test_that("the per-estimate accessors follow the rows", {
    enrich <- tibble(annotation = "a1", tau = 1e-7, enrichment = 2)
    jk <- matrix(0.1, 3L, 1L, dimnames = list(NULL, "a1"))
    ss <- list(z = 1.5, R = diag(1), annotationNames = "c1")
    x <- c(
        .h2Est(),
        .h2Est(
            study = "s2",
            h2 = 0.2,
            enrichment = list(enrich),
            jackknife = list(jk),
            score_stats = list(ss)
        )
    )
    expect_equal(heritabilityEstimate(x), c(0.3, 0.2))
    expect_identical(studyName(x), c("s1", "s2"))
    expect_identical(annotationJackknifeCoefs(x), list(NULL, jk))
    expect_identical(scoreStats(x), list(NULL, ss))
    # Each stratified estimate's table, behind its identity.
    tbl <- stratifiedHeritabilityEnrichment(x)
    expect_identical(tbl$study, "s2")
    expect_identical(tbl$annotation, "a1")
    expect_identical(colnames(tbl)[1:3], .H2_ID_COLS)
    # None stratified: no rows, the identity columns still there.
    none <- stratifiedHeritabilityEnrichment(.h2Est())
    expect_identical(nrow(none), 0L)
    expect_named(none, .H2_ID_COLS)
})

# ===========================================================================
# updateObject: estimates saved in slots, and with a `method` column
# ===========================================================================

# H2Estimates saved before it held its estimates as rows: the shipped example
# (stratified) and an LDER estimate with per-block local estimates, as
# estimateH2() returned them.
.h2Saved <- function() {
    readRDS(test_path("test_data", "h2EstimateSaved.rds"))
}

test_that("updateObject() rebuilds an H2Estimate saved in slots as a row", {
    saved <- .h2Saved()$example
    expect_no_warning(up <- updateObject(saved))
    expect_s4_class(up, "H2Estimate")
    expect_true(validObject(up))
    expect_identical(length(up), 1L)
    expect_length(up[[1]], 0L)
    expect_identical(studyName(up), "test")
    # Its method's settings, stratified as its enrichment table says.
    expect_identical(
        h2EstimateParam(up),
        .h2GenomeWideSettings("lder", TRUE, inSample = FALSE)
    )
    expect_equal(heritabilityEstimate(up), 0.3)
    expect_equal(heritabilityIntercept(up), 1.01)
    expect_identical(nSnps(up), 100L)
    expect_identical(
        annotationJackknifeCoefs(up),
        list(methods::slot(saved, "annotationJackknifeCoefs"))
    )
    expect_equal(
        as.data.frame(stratifiedHeritabilityEnrichment(up)[-(1:3)]),
        methods::slot(saved, "enrichment")
    )
    expect_identical(updateObject(up), up)
    expect_message(updateObject(up, verbose = TRUE), "H2Estimate")
    # ...and one saved under the slots' earlier names, the same way.
    older <- .savedUnder(
        saved,
        c(
            h2 = "estimate",
            h2Se = "estimateSe",
            local = "localBlocks",
            tauBlocks = "annotationJackknifeCoefs"
        )
    )
    expect_identical(updateObject(older), up)
})

test_that("updateObject() gives an estimate saved without an overlap NA", {
    current <- .h2Est()
    saved <- S4Vectors::`mcols<-`(
        current,
        value = S4Vectors::mcols(current, use.names = FALSE)[
            setdiff(
                colnames(S4Vectors::mcols(current)),
                "reference_overlap"
            )
        ]
    )
    expect_error(validObject(saved), "reference_overlap")
    up <- updateObject(saved)
    expect_true(validObject(up))
    expect_true(is.na(S4Vectors::mcols(up)$reference_overlap))
    expect_identical(updateObject(current), current)
})

test_that("updateObject() drops the per-block local estimates", {
    saved <- .h2Saved()$local
    expect_warning(
        up <- updateObject(saved),
        "per-block local estimates are dropped"
    )
    expect_identical(length(up), 1L)
    expect_equal(
        heritabilityEstimate(up),
        methods::slot(saved, "estimate")
    )
})

# The shipped example as saved with a `method` column, before the settings
# were recorded.
.h2SavedWithMethod <- function() {
    readRDS(test_path("test_data", "h2EstimateMethodColumn.rds"))
}

# `x` as it was saved with `method` naming each row's estimator: the column
# back, the settings slot gone.
test_that("updateObject() gives an estimate saved before intervals none", {
    x <- .h2Est()
    old <- x
    md <- S4Vectors::mcols(old, use.names = FALSE)
    attr(old, "elementMetadata") <- md[,
        setdiff(colnames(md), c("h2_lower", "h2_upper"))
    ]
    up <- updateObject(old)
    expect_true(validObject(up))
    expect_identical(S4Vectors::mcols(up)$h2_lower, NA_real_)
    expect_identical(S4Vectors::mcols(up)$h2_upper, NA_real_)
    expect_identical(heritabilityEstimate(up), heritabilityEstimate(x))
})

test_that("updateObject() gives an estimate saved before liability none", {
    x <- .h2Est()
    liability <- c(
        "h2_liability",
        "h2_liability_se",
        "h2_liability_lower",
        "h2_liability_upper"
    )
    old <- x
    md <- S4Vectors::mcols(old, use.names = FALSE)
    attr(old, "elementMetadata") <- md[, setdiff(colnames(md), liability)]
    expect_error(validObject(old), "h2_liability")
    up <- updateObject(old)
    expect_true(validObject(up))
    expect_true(all(is.na(unlist(as.list(S4Vectors::mcols(up)[liability])))))
    expect_identical(heritabilityEstimate(up), heritabilityEstimate(x))
    expect_identical(updateObject(x), x)
})

.savedWithMethod <- function(x, method) {
    old <- x
    md <- S4Vectors::mcols(old, use.names = FALSE)
    attr(old, "elementMetadata") <- cbind(
        md,
        S4Vectors::DataFrame(method = method)
    )
    attr(old, "param") <- NULL
    old
}

test_that("updateObject() turns a saved method column into its settings", {
    saved <- .h2SavedWithMethod()
    expect_true(is.null(attr(saved, "param")))
    up <- updateObject(saved)
    expect_true(validObject(up))
    expect_false(is_in("method", colnames(S4Vectors::mcols(up))))
    expect_identical(
        h2EstimateParam(up),
        .h2GenomeWideSettings("gldsc", TRUE, inSample = FALSE)
    )
    data(h2EstimateExample, envir = environment())
    expect_equal(
        heritabilityEstimate(up),
        heritabilityEstimate(h2EstimateExample)
    )
    expect_identical(h2EstimateParam(up), h2EstimateParam(h2EstimateExample))
    # A per-trait estimator's own settings.
    cis <- .h2Est(param = H2EstimateParam(), trait = "ENSG1")
    expect_identical(
        h2EstimateParam(updateObject(.savedWithMethod(cis, "hdlL"))),
        .h2LocalSettings("hdlL")
    )
})

test_that("updateObject() refuses estimates by several methods in one", {
    two <- c(.h2Est(), .h2Est(study = "s2"))
    expect_error(
        updateObject(.savedWithMethod(two, c("lder", "hdl"))),
        "holds estimates by lder, hdl"
    )
    expect_error(
        updateObject(.savedWithMethod(.h2Est(), "unknown")),
        "Rebuild the estimates"
    )
})
