test_that("H2EstimateParam's defaults are HDL-L's settings", {
    p <- H2EstimateParam()
    expect_s4_class(p, "H2EstimateParam")
    expect_identical(p$statistic, "eigen")
    expect_identical(p$criterion, "likelihood")
    expect_identical(p$noise, "free")
    expect_identical(p$crossNoise, "window")
    expect_identical(p$effects, "random")
    expect_identical(p$alpha, -1)
    expect_identical(p$referenceNoise, "eigenTerm")
    expect_identical(p$sampleSize, "median")
    expect_identical(p$inference, "lrt")
    expect_identical(p$procedures, H2LikelihoodParam())
    expect_identical(
        p$truncation,
        H2TruncationParam(
            fraction = 0.99,
            of = "variants",
            maxPerSample = NULL,
            minimum = 1
        )
    )
    expect_identical(.h2LocalSettings("hdlL"), p)
    # LAVA's are given in full, whatever the defaults.
    lava <- .h2LocalSettings("lava")
    expect_identical(lava$criterion, "moments")
    expect_identical(lava$effects, "fixed")
    expect_identical(lava$noise, "tied")
    expect_identical(lava$crossNoise, "supplied")
    expect_identical(lava$referenceNoise, "overlapFactor")
    expect_identical(lava$sampleSize, "mean")
    expect_identical(lava$inference, "beta")
    expect_null(lava$procedures)
    expect_identical(
        lava$truncation,
        H2TruncationParam(
            fraction = 0.99,
            of = "eigenvalues",
            maxPerSample = 0.75,
            minimum = 2
        )
    )
})

test_that("each criterion's procedures default to the setting paired with it", {
    lik <- H2EstimateParam(
        criterion = "likelihood",
        noise = "free",
        crossNoise = "window",
        referenceNoise = "eigenTerm",
        inference = "lrt"
    )
    expect_identical(lik$procedures, H2LikelihoodParam())
    # HDL-L's procedures are the likelihood's defaults.
    expect_identical(.h2LocalSettings("hdlL")$procedures, H2LikelihoodParam())
    reg <- H2EstimateParam(
        criterion = "regression",
        noise = "free",
        referenceNoise = "none",
        inference = "information"
    )
    expect_identical(reg$procedures, H2RegressionParam())
    gls <- H2EstimateParam(
        statistic = "variant",
        criterion = "gls",
        noise = "free",
        referenceNoise = "none",
        inference = "information"
    )
    expect_identical(gls$procedures, H2GlsParam())
    reml <- .h2LocalSettings("greml")
    expect_s4_class(reml$procedures, "MethodOptions")
    expect_identical(S4Vectors::metadata(reml$procedures)$engine, "greml")
})

test_that("the weights and steps default to the local settings' own", {
    reg <- H2RegressionParam()
    expect_true(reg$fitWeightIntercept)
    expect_identical(reg$initialWeights, "flat")
    expect_identical(H2LikelihoodParam()$maxit, 1000)
    # A genetic correlation's pairs are fitted one by one, as HDL-L's are.
    expect_false(H2LikelihoodParam()$jointFit)
})

test_that("only the eigen statistic keeps a truncation record", {
    gls <- H2EstimateParam(
        statistic = "variant",
        criterion = "gls",
        noise = "free",
        referenceNoise = "none",
        inference = "information"
    )
    expect_null(gls$truncation)
    expect_null(
        H2EstimateParam(
            statistic = "variant",
            criterion = "gls",
            noise = "free",
            referenceNoise = "none",
            inference = "information",
            truncation = NULL
        )$truncation
    )
    expect_error(
        H2EstimateParam(
            statistic = "variant",
            criterion = "gls",
            noise = "free",
            referenceNoise = "none",
            inference = "information",
            truncation = H2TruncationParam()
        ),
        "the \"variant\" statistic has none"
    )
    # NULL on the eigen statistic is its default record.
    expect_identical(
        H2EstimateParam(truncation = NULL)$truncation,
        H2TruncationParam()
    )
})

test_that("settings that do not belong together are refused", {
    refused <- function(..., pattern) {
        expect_error(H2EstimateParam(...), pattern)
    }
    refused(statistic = "phenotype", pattern = "goes with the reml criterion")
    refused(criterion = "reml", inference = "information", pattern = "reml")
    refused(
        statistic = "variant",
        criterion = "moments",
        referenceNoise = "none",
        pattern = "moments criterion needs the eigen statistic"
    )
    refused(
        criterion = "moments",
        noise = "free",
        pattern = "moments criterion needs tied noise"
    )
    refused(
        criterion = "moments",
        noise = "tied",
        crossNoise = "window",
        pattern = "fit the shared noise"
    )
    refused(effects = "mixed", pattern = "effects")
    expect_error(
        methods::initialize(H2EstimateParam(), effects = "mixed"),
        "effects must be one of random, fixed"
    )
    refused(alpha = Inf, pattern = "alpha must be one finite number")
    refused(alpha = c(0, -1), pattern = "alpha must be one finite number")
    refused(
        statistic = "variant",
        criterion = "likelihood",
        referenceNoise = "none",
        inference = "lrt",
        pattern = "likelihood criterion needs the eigen statistic"
    )
    refused(
        criterion = "gls",
        referenceNoise = "none",
        inference = "information",
        pattern = "gls criterion needs the variant statistic"
    )
    # Local and rg fits alike: the GLS has no tied or floored noise.
    refused(
        statistic = "variant",
        criterion = "gls",
        noise = "tied",
        referenceNoise = "none",
        inference = "information",
        pattern = "gls criterion fits a free noise term only"
    )
    refused(
        criterion = "regression",
        noise = "free",
        referenceNoise = "none",
        inference = "beta",
        pattern = "beta test goes with the moments"
    )
    lava <- .h2LocalSettings("lava")
    expect_error(
        methods::initialize(lava, inference = "lrt"),
        "goes with the likelihood criterion"
    )
    expect_error(
        methods::initialize(lava, inference = "information"),
        "reports the beta test"
    )
    expect_error(
        methods::initialize(lava, inference = "jackknife"),
        "jackknife goes with"
    )
    refused(
        criterion = "likelihood",
        noise = "free",
        referenceNoise = "overlapFactor",
        inference = "lrt",
        pattern = "overlapFactor reference correction goes with the moments"
    )
    refused(
        statistic = "variant",
        criterion = "regression",
        noise = "free",
        referenceNoise = "eigenTerm",
        inference = "information",
        pattern = "eigenTerm reference correction needs the eigen statistic"
    )
    refused(
        criterion = "likelihood",
        noise = "free",
        referenceNoise = "unbiasedR2",
        inference = "lrt",
        pattern = "unbiasedR2 reference correction needs the variant"
    )
    refused(
        criterion = "likelihood",
        noise = "free",
        referenceNoise = "weightsOnly",
        inference = "lrt",
        pattern = "weightsOnly reference goes with the eigen regression"
    )
    refused(
        criterion = "regression",
        noise = "floored",
        referenceNoise = "none",
        inference = "information",
        pattern = "floored noise needs a regression with an interceptFloor"
    )
})

test_that("the procedures must be the criterion's own, and fit its noise", {
    expect_error(
        methods::initialize(
            .h2LocalSettings("lava"),
            procedures = H2GlsParam()
        ),
        "moments criterion takes no procedures"
    )
    expect_error(
        H2EstimateParam(
            criterion = "likelihood",
            noise = "free",
            referenceNoise = "eigenTerm",
            inference = "lrt",
            procedures = H2RegressionParam()
        ),
        "likelihood criterion takes H2LikelihoodParam"
    )
    expect_error(
        H2EstimateParam(
            statistic = "phenotype",
            criterion = "reml",
            referenceNoise = "none",
            inference = "information",
            procedures = SusieOptions()
        ),
        "takes GremlOptions\\(\\)"
    )
    expect_error(
        H2EstimateParam(
            criterion = "likelihood",
            noise = "free",
            referenceNoise = "none",
            inference = "information",
            procedures = H2LikelihoodParam(profileScale = TRUE)
        ),
        "need tied noise"
    )
    expect_error(
        H2EstimateParam(
            criterion = "likelihood",
            referenceNoise = "eigenTerm",
            inference = "information",
            procedures = H2LikelihoodParam(nullSpace = TRUE)
        ),
        "use referenceNoise = \"none\""
    )
    expect_error(
        H2EstimateParam(
            criterion = "regression",
            noise = "tied",
            referenceNoise = "none",
            inference = "information",
            procedures = H2RegressionParam(twoStage = TRUE)
        ),
        "two stages fit a free or floored noise term"
    )
    # S-LDSC's correction weighs variants; LDER's taper, eigen directions.
    expect_error(
        H2EstimateParam(
            criterion = "regression",
            noise = "free",
            referenceNoise = "none",
            inference = "information",
            procedures = H2RegressionParam(ldRedundancy = TRUE)
        ),
        "the LD-redundancy correction weighs variants"
    )
    expect_error(
        H2EstimateParam(
            statistic = "variant",
            criterion = "regression",
            noise = "free",
            referenceNoise = "none",
            inference = "jackknife",
            procedures = H2RegressionParam(eigenvalueTaper = "linear")
        ),
        "the eigenvalue taper weighs eigen directions"
    )
})

test_that("the nested records check their own fields", {
    expect_error(H2TruncationParam(fraction = 0), "fraction must be")
    expect_error(H2TruncationParam(fraction = 1.5), "fraction must be")
    expect_error(H2TruncationParam(of = "trace"), "must be one of")
    expect_error(H2TruncationParam(maxPerSample = 2), "maxPerSample")
    expect_error(H2TruncationParam(minimum = 0.5), "minimum")
    expect_error(H2TruncationParam(keepTies = NA), "keepTies")
    expect_error(H2LikelihoodParam(h2Bounds = c(1, 0)), "h2Bounds")
    expect_error(H2LikelihoodParam(interceptBounds = 1), "interceptBounds")
    expect_error(H2LikelihoodParam(varianceFloor = -1), "varianceFloor")
    expect_error(H2LikelihoodParam(h2Starts = NA_real_), "h2Starts")
    expect_error(H2LikelihoodParam(nullSpace = NA), "nullSpace")
    expect_error(H2LikelihoodParam(jointFit = c(TRUE, FALSE)), "jointFit")
    expect_error(H2LikelihoodParam(maxit = 0), "maxit")
    expect_error(H2RegressionParam(weights = "flat"), "must be one of")
    expect_error(
        H2RegressionParam(fitWeightIntercept = NA),
        "fitWeightIntercept"
    )
    expect_error(H2RegressionParam(initialWeights = "unit"), "must be one of")
    expect_error(
        H2RegressionParam(weights = "fixed", iterations = 2),
        "iterations must be NULL"
    )
    expect_error(H2RegressionParam(iterations = 1.5), "iterations")
    expect_error(H2RegressionParam(chi2Threshold = -1), "chi2Threshold")
    expect_error(H2RegressionParam(chi2Filter = 80), "chi2Filter")
    expect_error(H2RegressionParam(jackknifeBlocks = 1), "jackknifeBlocks")
    expect_error(H2RegressionParam(eigenvalueTaper = "cubic"), "must be one")
    expect_error(H2RegressionParam(interceptFloor = NA_real_), "interceptFloor")
    expect_error(H2GlsParam(ridge = -1), "ridge")
    expect_error(H2GlsParam(omegaFraction = 0), "omegaFraction")
    expect_error(H2GlsParam(plugIn = "iterated"), "plugIn")
})

test_that("GremlOptions refuses what pecotmr supplies itself", {
    skip_if_not_installed("qgg")
    expect_s4_class(GremlOptions(maxit = 50), "MethodOptions")
    expect_error(GremlOptions(y = 1), "supplied by pecotmr")
    expect_error(GremlOptions(GRM = list()), "supplied by pecotmr")
    expect_error(GremlOptions(notAnArgument = 1), "notAnArgument")
})

test_that("the published estimators' settings are valid and as documented", {
    hdl <- .h2LocalSettings("hdlL")
    expect_identical(hdl$truncation$of, "variants")
    expect_null(hdl$truncation$maxPerSample)
    expect_identical(hdl$inference, "lrt")
    fastHer <- .h2LocalSettings("fastHer")
    expect_true(fastHer$procedures$nullSpace && fastHer$procedures$profileScale)
    expect_true(fastHer$truncation$keepTies)
    expect_identical(fastHer$sampleSize, "harmonic")
    expect_identical(.h2LocalSettings("greml")$statistic, "phenotype")
    for (method in c("sldsc", "gldsc", "lder", "hdl")) {
        for (stratified in c(FALSE, TRUE)) {
            expect_true(
                validObject(.h2GenomeWideSettings(method, stratified, FALSE)),
                info = method
            )
        }
    }
    # The kind of reference sets LDER's taper; the stratification, S-LDSC's
    # weights and filter and HDL's optimizer.
    expect_identical(
        .h2GenomeWideSettings("lder", FALSE, TRUE)$procedures$eigenvalueTaper,
        "linear"
    )
    expect_identical(
        .h2GenomeWideSettings("lder", FALSE, FALSE)$procedures$eigenvalueTaper,
        "logistic"
    )
    # LDER's weights hold the intercept at 1 and start flat; S-LDSC's start
    # at its moments; HDL steps the intercept as it steps h2.
    lder <- .h2GenomeWideSettings("lder", FALSE, FALSE)$procedures
    expect_false(lder$fitWeightIntercept)
    expect_identical(lder$initialWeights, "flat")
    expect_identical(
        .h2GenomeWideSettings("sldsc", FALSE, FALSE)$procedures$initialWeights,
        "moments"
    )
    expect_true(
        .h2GenomeWideSettings("hdl", FALSE, FALSE)$procedures$regressionStart
    )
    stratified <- .h2GenomeWideSettings("sldsc", TRUE, FALSE)$procedures
    expect_identical(stratified$weights, "fixed")
    expect_identical(stratified$initialWeights, "moments")
    expect_identical(stratified$chi2Filter, c(0.001, 80))
    expect_false(stratified$twoStage)
    # HDL's are the same stratified or not, every likelihood fitted alike.
    expect_identical(
        .h2GenomeWideSettings("hdl", TRUE, FALSE),
        .h2GenomeWideSettings("hdl", FALSE, FALSE)
    )
})
