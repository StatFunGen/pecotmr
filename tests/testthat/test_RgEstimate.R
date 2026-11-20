context("RgEstimate")

# A two-pair RgEstimate built directly: context a against b and c of one gene,
# under `param`, with `values` overriding the defaults.
.rgTestEstimate <- function(
    param = H2EstimateParam(),
    ...,
    conditional = FALSE
) {
    defaults <- list(
        study_1 = c("s1", "s1"),
        context_1 = c("a", "a"),
        trait_1 = c("g1", "g1"),
        study_2 = c("s1", "s1"),
        context_2 = c("b", "c"),
        trait_2 = c("g1", "g1"),
        rg = c(0.5, -0.2),
        rg_se = c(NA_real_, NA_real_),
        rg_lower = c(0.1, -0.7),
        rg_upper = c(0.9, 0.3),
        covariance = c(0.05, -0.01),
        covariance_se = c(NA_real_, NA_real_),
        h2_1 = c(0.1, 0.1),
        h2_2 = c(0.1, 0.05),
        covariance_liability = c(NA_real_, NA_real_),
        covariance_liability_se = c(NA_real_, NA_real_),
        h2_1_liability = c(NA_real_, NA_real_),
        h2_2_liability = c(NA_real_, NA_real_),
        noise_correlation = c(0.2, 0),
        noise_source = c("window", "window"),
        p_value = c(0.01, 0.4),
        r2_1 = c(NA_real_, NA_real_),
        r2_1_lower = c(NA_real_, NA_real_),
        r2_1_upper = c(NA_real_, NA_real_),
        r2_2 = c(NA_real_, NA_real_),
        r2_2_lower = c(NA_real_, NA_real_),
        r2_2_upper = c(NA_real_, NA_real_),
        n_snps = c(100L, 100L),
        n_components = c(40L, 40L),
        rg_loglik = list(NULL, NULL)
    )
    values <- utils::modifyList(defaults, list(...))
    window <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 1000))
    pecotmr:::.rgEstimateNew(
        list(window, window),
        values,
        window,
        param,
        conditional
    )
}

test_that("an RgEstimate holds one row per pair, keyed by both conditions", {
    x <- .rgTestEstimate()
    expect_s4_class(x, "RgEstimate")
    expect_length(x, 2L)
    expect_identical(names(x)[[1L]], "s1|a|g1|s1|b|g1")
    expect_identical(h2EstimateParam(x), H2EstimateParam())
    expect_identical(x$rg, c(0.5, -0.2))
    out <- capture.output(show(x))
    expect_match(out[[1L]], "RgEstimate: 2 pair\\(s\\)")
    expect_match(out[[2L]], "cross noise window")
    # A subset keeps the settings.
    expect_identical(h2EstimateParam(x[1L]), h2EstimateParam(x))
    expect_length(x[1L], 1L)
})

test_that("an empty RgEstimate shows its settings alone", {
    x <- .rgTestEstimate()[0L]
    out <- capture.output(show(x))
    expect_match(out[[1L]], "RgEstimate: 0 pair\\(s\\)")
    expect_length(out, 2L)
})

test_that("RgEstimate validity refuses a malformed table", {
    expect_error(
        .rgTestEstimate(rg = c(1L, 0L)),
        "the `rg` column must be double"
    )
    expect_error(
        .rgTestEstimate(noise_source = c("window", "mash")),
        "the `noise_source` column holds mash"
    )
    expect_error(
        .rgTestEstimate(context_2 = c("b", "b")),
        "pair uniqueness violated"
    )
    expect_error(
        .rgTestEstimate(study_1 = c(1, 1)),
        "the `study_1` column must be character"
    )
    window <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 1000))
    expect_error(
        pecotmr:::.rgEstimateNew(
            list(window),
            list(study_1 = "s1", rg = 0.5),
            window,
            H2EstimateParam()
        ),
        "columns"
    )
})

test_that("c() keeps one set of settings and refuses mixed ones", {
    x <- .rgTestEstimate()
    other <- .rgTestEstimate(context_1 = c("d", "d"))
    both <- c(x, other)
    expect_s4_class(both, "RgEstimate")
    expect_length(both, 4L)
    expect_identical(h2EstimateParam(both), h2EstimateParam(x))
    zero <- .rgTestEstimate(
        param = methods::initialize(
            H2EstimateParam(),
            crossNoise = "zero"
        ),
        context_1 = c("d", "d")
    )
    expect_error(
        c(x, zero),
        "the RgEstimate objects were fitted under different settings"
    )
})

test_that("an RgEstimate records whether its pairs are given their group", {
    x <- .rgTestEstimate()
    given <- .rgTestEstimate(conditional = TRUE, context_1 = c("d", "d"))
    expect_false(isConditional(x))
    expect_true(isConditional(given))
    expect_true(isConditional(given[1L]))
    expect_match(
        capture.output(show(given))[[3L]],
        "each pair given its group's other conditions"
    )
    expect_length(capture.output(show(x[0L])), 2L)
    expect_true(isConditional(c(given, given[0L])))
    expect_error(
        c(x, given),
        "some of the RgEstimate objects report each pair given its group"
    )
    expect_error(
        .rgTestEstimate(conditional = NA),
        "conditional"
    )
})

test_that("an RgEstimate's curves and noise sources are checked", {
    x <- .rgTestEstimate()
    curve <- x
    expect_error(
        {
            S4Vectors::mcols(curve)$rg_loglik <- list(1:3, NULL)
            validObject(curve)
        },
        "rg_loglik"
    )
    source <- x
    expect_error(
        {
            S4Vectors::mcols(source)$noise_source <- c(1, 2)
            validObject(source)
        },
        "noise_source"
    )
})
