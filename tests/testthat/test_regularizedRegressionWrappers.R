context("regularized_regression - wrappers")

# Tests of R/regularizedRegressionWrappers.R: the mr.mash helpers, the weight
# wrappers' engine arguments and fit retention, and their real-fit coverage.
# Each method's own tests live in the matching test_rr*.R topic file.

# =============================================================================
# mrmashWeights fit retention
# =============================================================================

test_that("mrmashWeights retains nothing, the prior parts, or the full fit", {
    skip_if_not_installed("mr.mashr")
    # The slim payload is exactly what fineMappingPipeline needs to rebuild the
    # mvSuSiE reweighted mixture prior (w0 -> rescaleCovW0, original $U) and
    # the residual variance (V); the heavy mu1 coefficient matrix is not
    # retained.
    fakeFit <- list(w0 = c(a_1 = 0.5, a_2 = 0.5), V = diag(2))
    ddpm <- list(U = list(a = diag(2)))
    # Mock coef extraction so we exercise only the retain payload logic, not a
    # real mr.mash fit. coef.mr.mash(fit)[-1, ] -> drop the intercept row.
    local_mocked_bindings(
        coef.mr.mash = function(object, ...) rbind(c(0, 0), c(0.1, 0.2)),
        .package = "mr.mashr"
    )
    fitSlim <- attr(
        mrmashWeights(
            mrmashFit = fakeFit,
            fitRetention = "slim",
            dataDrivenPriorMatrices = ddpm
        ),
        "fit"
    )
    expect_true(is.list(fitSlim))
    expect_setequal(names(fitSlim), c("dataDrivenPriorMatrices", "w0", "V"))
    expect_null(fitSlim$fit) # slim: no full fit
    expect_identical(fitSlim$dataDrivenPriorMatrices, ddpm)
    expect_identical(fitSlim$w0, fakeFit$w0)
    expect_identical(fitSlim$V, fakeFit$V)
    # The default ("none") leaves the weights free of the fit attribute.
    expect_null(attr(
        mrmashWeights(mrmashFit = fakeFit, dataDrivenPriorMatrices = ddpm),
        "fit"
    ))

    fitFull <- attr(
        mrmashWeights(
            mrmashFit = fakeFit,
            fitRetention = "full",
            dataDrivenPriorMatrices = ddpm
        ),
        "fit"
    )
    expect_true("fit" %in% names(fitFull)) # full: the whole fit retained
    expect_identical(fitFull$fit, fakeFit)
    expect_identical(fitFull$w0, fakeFit$w0) # slim fields still present
})

# =============================================================================
# mr.mash helpers: computeW0
# =============================================================================

test_that("computeW0 returns uniform weights when ncomps == 1", {
    Bhat <- matrix(c(1, 0, 0, 2, 0, 0), nrow = 3, ncol = 2)
    result <- pecotmr:::computeW0(Bhat, ncomps = 1)
    expect_equal(result, 1)
})


test_that("computeW0 handles all-zero Bhat by returning uniform weights", {
    # When Bhat is all zero, prop_nonzero = 0
    # w0 = c(1, 0, ..., 0) => sum(w0 != 0) < 2 => fallback to uniform
    Bhat <- matrix(0, nrow = 5, ncol = 3)
    result <- pecotmr:::computeW0(Bhat, ncomps = 4)
    expect_equal(result, rep(1 / 4, 4))
    expect_equal(sum(result), 1)
})


test_that("computeW0 weights by nonzero rows when ncomps > 1", {
    # 2 out of 4 rows have nonzero entries
    Bhat <- matrix(0, nrow = 4, ncol = 2)
    Bhat[1, 1] <- 1
    Bhat[3, 2] <- 2
    result <- pecotmr:::computeW0(Bhat, ncomps = 3)
    expect_equal(length(result), 3)
    expect_equal(sum(result), 1, tolerance = 1e-10)
    # First element should be (1 - prop_nonzero) = 0.5
    expect_equal(result[1], 0.5)
})

# =========================================================================
# rescaleCovW0
# =========================================================================

test_that("rescaleCovW0 removes null component and renormalizes", {
    w0 <- c(
        null = 0.3,
        XtX_1 = 0.2,
        XtX_2 = 0.1,
        FLASH_1 = 0.15,
        FLASH_2 = 0.25
    )
    result <- pecotmr:::rescaleCovW0(w0)
    expect_false("null" %in% names(result))
    expect_equal(sum(result), 1, tolerance = 1e-10)
})


test_that("rescaleCovW0 handles all-zero non-null weights", {
    w0 <- c(null = 1.0, XtX_1 = 0, XtX_2 = 0, FLASH_1 = 0)
    result <- pecotmr:::rescaleCovW0(w0)
    # All non-null weights are zero -> equal weights
    expect_equal(sum(result), 1, tolerance = 1e-10)
    expect_true(all(result == result[1])) # all equal
})


test_that("rescaleCovW0 groups correctly by prior group prefix", {
    w0 <- c(
        null = 0.5,
        PCA_1 = 0.1,
        PCA_2 = 0.2,
        tFLASH_1 = 0.1,
        tFLASH_2 = 0.1
    )
    result <- pecotmr:::rescaleCovW0(w0)
    expect_true("PCA" %in% names(result))
    expect_true("tFLASH" %in% names(result))
    expect_equal(sum(result), 1, tolerance = 1e-10)
})

# =========================================================================
# computeGrid (mr.mash's own exported grid builder)
# =========================================================================

test_that("computeGrid is mr.mash's grid, not a reimplementation of it", {
    skip_if_not_installed("mr.mashr")
    # gridMin/gridMax/autoselectMixsd used to be ported copies here, and
    # gridMin had drifted: upstream grid_min() is min(Shat)/10, the port
    # dropped the /10, so every grid started ten times too high.
    set.seed(3)
    bhat <- matrix(rnorm(20, sd = 2), nrow = 10, ncol = 2)
    sbhat <- matrix(abs(rnorm(20, mean = 0.5, sd = 0.1)), nrow = 10, ncol = 2)
    grid <- pecotmr:::computeGrid(bhat, sbhat)
    expect_equal(
        grid,
        mr.mashr::autoselect.mixsd(
            list(Bhat = bhat, Shat = sbhat),
            mult = sqrt(2)
        )^2
    )
    expect_true(all(grid > 0))
    # The ported floor would have been 10^2 times this one, since the grid
    # scales variances.
    ported <- min(sbhat)
    expect_lt(min(grid), ported^2)
})


test_that("compute_grid handles NA and zero sbhat values", {
    bhat <- c(1, 2, NA, 4, 5)
    sbhat <- c(0.5, 0, NA, 0.3, 0.8)
    result <- pecotmr:::computeGrid(bhat, sbhat)
    expect_true(is.numeric(result))
    expect_true(length(result) > 0)
})

# =========================================================================
# mrmashWrapper input validation
# =========================================================================

# Note: Cannot mock requireNamespace via local_mocked_bindings because it is a
# base R function, not in pecotmr's namespace. We skip these tests and instead
# test the downstream validation that we CAN exercise.

test_that("mrmashWrapper errors when X and Y are not matrices", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mr.mashr")
    expect_error(
        mrmashWrapper(data.frame(x = 1:3), matrix(1:6, nrow = 3, ncol = 2)),
        "matrices"
    )
})


test_that("mrmashWrapper errors when X and Y row counts differ", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mr.mashr")
    expect_error(
        mrmashWrapper(
            matrix(1:6, nrow = 3, ncol = 2),
            matrix(1:8, nrow = 4, ncol = 2)
        ),
        "Assertion on 'Y'.*Must have exactly 3 rows"
    )
})


test_that("MrmashPriorParam refuses a priorGrid that is not a vector", {
    # The typed slot catches this when the record is built, so the wrapper's
    # own is.vector() guard became unreachable and was deleted.
    expect_error(
        MrmashPriorParam(priorGrid = matrix(1:4, nrow = 2)),
        'invalid object for slot "priorGrid"'
    )
    expect_equal(
        MrmashPriorParam(priorGrid = c(0.1, 0.2))$priorGrid,
        c(0.1, 0.2)
    )
})


test_that("mrmashWrapper errors with no data-driven or canonical prior", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mr.mashr")
    X <- matrix(rnorm(12), nrow = 3, ncol = 4)
    Y <- matrix(rnorm(6), nrow = 3, ncol = 2)
    expect_error(
        mrmashWrapper(
            X,
            Y,
            dataDrivenPriorMatrices = NULL,
            prior = MrmashPriorParam(canonicalPriorMatrices = FALSE)
        ),
        "dataDrivenPriorMatrices"
    )
})


test_that("mrmashWrapper warns on missing Y with glasso initialisation", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mr.mashr")
    set.seed(42)
    n <- 20
    p <- 5
    r <- 2
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
    Y[1, 1] <- NA # introduce missing values
    colnames(Y) <- c("cond1", "cond2")

    # Should produce warning about glasso and NAs, then likely fail on the
    # downstream mr.mashr call, but the warning is what we test
    expect_warning(
        tryCatch(
            mrmashWrapper(
                X,
                Y,
                bInitMethod = "glasso",
                dataDrivenPriorMatrices = list(U = list(matrix(1, 2, 2))),
                prior = MrmashPriorParam(canonicalPriorMatrices = FALSE)
            ),
            error = function(e) NULL
        ),
        "glasso"
    )
})

# =========================================================================
# computeCoefficientsGlasso
# =========================================================================

test_that("computeCoefficientsGlasso runs without Xnew", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 50
    p <- 5
    r <- 3
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
    colnames(Y) <- paste0("cond", seq_len(r))
    result <- pecotmr:::computeCoefficientsGlasso(
        X,
        Y,
        standardize = FALSE,
        Xnew = NULL
    )
    expect_true("Bhat" %in% names(result))
    expect_true("Ytrain" %in% names(result))
    expect_equal(nrow(result$Bhat), p)
    expect_equal(ncol(result$Bhat), r)
    expect_null(result$Yhat_new)
})


test_that("computeCoefficientsGlasso runs with Xnew", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 50
    p <- 5
    r <- 3
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
    colnames(Y) <- paste0("cond", seq_len(r))
    Xnew <- matrix(rnorm(10 * p), nrow = 10, ncol = p)
    result <- pecotmr:::computeCoefficientsGlasso(
        X,
        Y,
        standardize = FALSE,
        Xnew = Xnew
    )
    expect_true("Yhat_new" %in% names(result))
    expect_equal(nrow(result$Yhat_new), 10)
    expect_equal(ncol(result$Yhat_new), r)
    expect_equal(colnames(result$Yhat_new), colnames(Y))
})

# =========================================================================
# computeCoefficientsUnivGlmnet
# =========================================================================

test_that("computeCoefficientsUnivGlmnet gives the same fits on two workers", {
    skip_if_not_installed("glmnet")
    skip_on_os("windows")
    set.seed(5)
    X <- matrix(rnorm(60 * 5), 60, 5)
    Y <- matrix(rnorm(60 * 2), 60, 2, dimnames = list(NULL, c("a", "b")))
    fit <- function(bp) {
        pecotmr:::computeCoefficientsUnivGlmnet(
            X,
            Y,
            alpha = 0.5,
            standardize = FALSE,
            BPPARAM = bp
        )
    }
    expect_equal(
        fit(BiocParallel::SerialParam(RNGseed = 9L)),
        fit(BiocParallel::MulticoreParam(2, RNGseed = 9L))
    )
})

test_that("computeCoefficientsUnivGlmnet runs without Xnew", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 60
    p <- 5
    r <- 2
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
    colnames(Y) <- paste0("cond", seq_len(r))
    result <- pecotmr:::computeCoefficientsUnivGlmnet(
        X,
        Y,
        alpha = 0.5,
        standardize = FALSE,
        Xnew = NULL
    )
    expect_true("Bhat" %in% names(result))
    expect_true("intercept" %in% names(result))
    expect_equal(nrow(result$Bhat), p)
    expect_equal(ncol(result$Bhat), r)
    expect_null(result$Yhat_new)
})


test_that("computeCoefficientsUnivGlmnet runs with Xnew", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 60
    p <- 5
    r <- 2
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
    colnames(Y) <- paste0("cond", seq_len(r))
    Xnew <- matrix(rnorm(8 * p), nrow = 8, ncol = p)
    result <- pecotmr:::computeCoefficientsUnivGlmnet(
        X,
        Y,
        alpha = 0.5,
        standardize = FALSE,
        Xnew = Xnew
    )
    expect_true("Yhat_new" %in% names(result))
    expect_equal(nrow(result$Yhat_new), 8)
    expect_equal(ncol(result$Yhat_new), r)
    expect_equal(colnames(result$Yhat_new), colnames(Y))
})


test_that("computeCoefficientsUnivGlmnet handles NA in Y", {
    skip_if_not_installed("glmnet")
    set.seed(42)
    n <- 60
    p <- 5
    r <- 2
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
    colnames(Y) <- paste0("cond", seq_len(r))
    Y[1:5, 1] <- NA # introduce missing values in one condition
    result <- pecotmr:::computeCoefficientsUnivGlmnet(
        X,
        Y,
        alpha = 0.5,
        standardize = FALSE,
        Xnew = NULL
    )
    expect_true("Bhat" %in% names(result))
    expect_equal(nrow(result$Bhat), p)
})

# =========================================================================
# mrmashWrapper seed message
# =========================================================================

# Note: Cannot mock base::exists() via local_mocked_bindings.
# The seed-check message would require removing .Random.seed from the global
# environment, which is not safe to do in tests.

# =============================================================================
# Real-fit coverage for the solver wrappers in R/regularizedRegressionWrappers.R
# -----------------------------------------------------------------------------
# The fine-mapping / TWAS pipelines MOCK these wrappers. Here we drive the
# wrappers on a SMALL real fixture and assert the return shape (weight length
# == #variants; matrix for multivariate; attr(.,"fit") unless fitRetention is
# "none"). MCMC iterations kept tiny. Each method's shape and signal-recovery
# tests live in its test_rr*.R topic file; fsusieWeights is tested in
# test_fineMappingWrappers.R.
# =============================================================================

# (Shared .rrwXy / .rrwStatLd / .rrwMulti fixtures live in
#  helper-rrwFixtures.R so test_fineMappingWrappers.R can reuse them for
#  the SuSiE weight extractors.)

# -------------------------------- individual --------------------------------

# Capture what reaches cv.glmnet, which is all fitGlmnetCv() is responsible
# for -- the engine itself is glmnet's. The fake fit is not a usable
# cv.glmnet object, so an extractor downstream of the call errors; `try`
# swallows that, since only the captured arguments are under test.
.rrw_captureGlmnet <- function(expr) {
    captured <- NULL
    local_mocked_bindings(
        cv.glmnet = function(...) {
            captured <<- list(...)
            structure(list(), class = "cv.glmnet")
        },
        .package = "glmnet"
    )
    try(force(expr), silent = TRUE)
    captured
}

test_that("fitGlmnetCv: NULL formals leave glmnet's own defaults", {
    skip_if_not_installed("glmnet")
    f <- .rrwXy()
    args <- .rrw_captureGlmnet(fitGlmnetCv(f$X, f$y, alpha = 0.5))
    expect_equal(args$alpha, 0.5)
    # An explicit NULL would reach cv.glmnet and override its default with
    # nothing, so an unset formal must be ABSENT from the call.
    for (nm in c("family", "standardize", "intercept", "nfolds")) {
        expect_false(is_in(nm, names(args)))
    }
})

test_that("fitGlmnetCv: set formals use glmnet's own argument names", {
    skip_if_not_installed("glmnet")
    f <- .rrwXy()
    args <- .rrw_captureGlmnet(fitGlmnetCv(
        f$X,
        f$y,
        alpha = 1,
        family = "mgaussian",
        standardize = TRUE,
        intercept = FALSE,
        nfolds = 3
    ))
    expect_equal(args$family, "mgaussian")
    expect_true(args$standardize)
    expect_false(args$intercept)
    expect_equal(args$nfolds, 3)
    expect_false(is_in("nfold", names(args)))
})

test_that("fitGlmnetCv: methodArgs win over the caller's defaults", {
    skip_if_not_installed("glmnet")
    f <- .rrwXy()
    args <- .rrw_captureGlmnet(fitGlmnetCv(
        f$X,
        f$y,
        alpha = 0.5,
        nfolds = 5,
        methodArgs = GlmnetOptions(nfolds = 3)
    ))
    expect_equal(args$nfolds, 3)
    expect_error(
        fitGlmnetCv(f$X, f$y, alpha = 0.5, methodArgs = list(nfolds = 3)),
        "GlmnetOptions"
    )
})

test_that("GlmnetOptions: nfold is refused and names nfolds instead", {
    expect_error(GlmnetOptions(nfold = 5), "fold count is `nfolds`")
})

test_that("GlmnetOptions: parallel is refused outright", {
    # glmnet's `parallel` needs a foreach backend pecotmr does not declare,
    # and weight fitting already runs inside BiocParallel workers over the
    # same CV folds -- so there is no value a caller could usefully set.
    # cv.glmnet takes `...`, so nothing else would report it.
    expect_error(GlmnetOptions(parallel = TRUE), "cannot be set here")
    expect_error(GlmnetOptions(parallel = FALSE), "foreach backend")
    # A refusal must not swallow the names that ARE settable.
    expect_equal(GlmnetOptions(nfolds = 3)$nfolds, 3)
})

test_that("fitGlmnetCv: the engine is never asked to parallelise", {
    skip_if_not_installed("glmnet")
    f <- .rrwXy()
    args <- .rrw_captureGlmnet(fitGlmnetCv(f$X, f$y, alpha = 0.5))
    expect_false(is_in("parallel", names(args)))
})

test_that("glmnetWeights: the fold count reaches the engine as nfolds", {
    skip_if_not_installed("glmnet")
    f <- .rrwXy()
    args <- .rrw_captureGlmnet(glmnetWeights(f$X, f$y, alpha = 0.5))
    expect_equal(args$nfolds, 5)
    expect_false(args$standardize)
    expect_true(args$intercept)
    # Regression: pecotmr passed `nfold = 5`, which reached cv.glmnet only by
    # partial matching, so a user who spelled the engine's own `nfolds`
    # supplied BOTH and glmnet failed with a message about `family`.
    w <- glmnetWeights(f$X, f$y, 0.5, methodArgs = GlmnetOptions(nfolds = 3))
    expect_length(as.numeric(w), f$p)
})

test_that("mr.mash glmnet initialisation routes through fitGlmnetCv", {
    skip_if_not_installed("glmnet")
    f <- .rrwXy()
    Y <- cbind(f$y, f$y + rnorm(length(f$y)))
    gl <- .rrw_captureGlmnet(computeCoefficientsGlasso(
        f$X,
        Y,
        standardize = TRUE
    ))
    expect_equal(gl$family, "mgaussian")
    expect_equal(gl$alpha, 1)
    expect_true(gl$standardize)
    # `parallel = FALSE` was a restatement of cv.glmnet's own default.
    expect_false(is_in("parallel", names(gl)))
    uv <- .rrw_captureGlmnet(pecotmr:::.linreg(
        1L,
        f$X,
        Y,
        alpha = 0.5,
        standardize = FALSE,
        Xnew = NULL
    ))
    expect_equal(uv$family, "gaussian")
    expect_equal(uv$alpha, 0.5)
    expect_false(uv$standardize)
    expect_false(is_in("parallel", names(uv)))
})

test_that("mrashWeights returns length-p weights and can retain the fit", {
    skip_if_not_installed("susieR")
    skip_if_not_installed("glmnet")
    f <- .rrwXy()
    # mr.ash's own default grid is too narrow for this fixture's two strong
    # effects, and mr.ash says so.
    expect_warning(
        w <- mrashWeights(f$X, f$y, fitRetention = "slim"),
        "larger setting of the prior variance"
    )
    expect_length(w, f$p)
    expect_false(is.null(attr(w, "fit")))
})

test_that("qgg Bayes-alphabet weights (N/L/A/C/R) return length-p weights", {
    skip_if_not_installed("qgg")
    f <- .rrwXy()
    mc <- list(nIter = 200, nBurnin = 20, thin = 1)
    fitArgs <- c(list(f$X, f$y), mc)
    wN <- exec(bayesNWeights, !!!fitArgs)
    expect_length(wN, f$p)
    expect_true(is.numeric(wN))
    expect_true(all(is.finite(wN)))
    expect_length(exec(bayesLWeights, !!!fitArgs), f$p)
    expect_length(exec(bayesAWeights, !!!fitArgs), f$p)
    expect_length(exec(bayesCWeights, !!!fitArgs), f$p)
    expect_length(exec(bayesRWeights, !!!fitArgs), f$p)
})

test_that("buildMrmashPriorMatrices exposes expand_covs' zeromat", {
    skip_if_not_installed("mr.mashr")
    # compute_canonical_covs() was already fully exposed (singletons,
    # hetgrid) but expand_covs()'s `zeromat` was hardcoded TRUE, so the
    # null component could not be dropped.
    seen <- NULL
    # The default must match the real expand_covs(): pecotmr no longer
    # passes zeromat unless the caller set it, so an omitted argument has to
    # fall through to the engine's own default here too.
    local_mocked_bindings(
        expand_covs = function(mats, grid, zeromat = TRUE) {
            seen <<- zeromat
            mats
        },
        .package = "mr.mashr"
    )
    set.seed(5)
    Bhat <- matrix(rnorm(12), 6, 2)
    Shat <- matrix(abs(rnorm(12, 0.5, 0.1)), 6, 2)
    invisible(buildMrmashPriorMatrices(
        Bhat,
        Shat,
        expandCovs = MrmashExpandCovsOptions(zeromat = FALSE)
    ))
    expect_false(seen)
    invisible(buildMrmashPriorMatrices(Bhat, Shat))
    expect_true(seen)
})

test_that("bayesAlphabetWeights forwards every MCMC control to gbayes", {
    # `nthin` was declared and documented (default 5) but never put in
    # callArgs, so every fit silently ran at gbayes's own nthin = 1 while
    # its siblings nit and nburn were forwarded. gbayes runs nit + nburn
    # iterations, so the total nIter reaches it net of the burn-in.
    seen <- NULL
    local_mocked_bindings(
        gbayes = function(...) {
            seen <<- list(...)
            list(bm = rep(0, 5))
        },
        .package = "qgg"
    )
    f <- .rrwXy(p = 5)
    invisible(bayesAlphabetWeights(
        f$X,
        f$y,
        method = "bayesN",
        nIter = 300,
        nBurnin = 30,
        thin = 7
    ))
    expect_equal(seen$nit, 270)
    expect_equal(seen$nburn, 30)
    expect_equal(seen$nthin, 7)
})

test_that("dprVbWeights returns length-p weights and retains the fit", {
    skip_if_not_installed("RcppDPR")
    f <- .rrwXy()
    w <- dprVbWeights(f$X, f$y, fitRetention = "slim")
    expect_length(w, f$p)
    expect_false(is.null(attr(w, "fit")))
})

test_that("dprGibbsWeights returns length-p weights", {
    skip_if_not_installed("RcppDPR")
    f <- .rrwXy()
    # Capturing the DPR engine's progress bar is what is slow, not sampling.
    invisible(capture.output(
        w <- dprGibbsWeights(
            f$X,
            f$y,
            sStep = 200,
            methodArgs = DprOptions(show_progress = FALSE)
        )
    ))
    expect_length(w, f$p)
})

test_that("dprAdaptiveGibbsWeights returns length-p weights", {
    skip_if_not_installed("RcppDPR")
    f <- .rrwXy()
    invisible(capture.output(
        w <- dprAdaptiveGibbsWeights(
            f$X,
            f$y,
            methodArgs = DprOptions(s_step = 100, show_progress = FALSE)
        )
    ))
    expect_length(w, f$p)
})

test_that("mrmashWeights fits from (X, Y) and returns p x K weights", {
    skip_if_not_installed("mr.mashr")
    skip_if_not_installed("glmnet")
    set.seed(3)
    m <- .rrwMulti(n = 60, p = 6, K = 3)
    w <- suppressMessages(mrmashWeights(
        X = m$X,
        Y = m$Y,
        methodArgs = MrmashOptions(canonicalPriorMatrices = TRUE)
    ))
    expect_equal(dim(w), c(m$p, m$K))
    expect_true(all(is.finite(w)))
})


# --------------------------- RSS weight wrappers ----------------------------

test_that("lassosumRssWeights returns p weights and records the selection", {
    f <- .rrwStatLd()
    set.seed(1)
    w <- lassosumRssWeights(f$stat, f$ldReference)
    expect_length(w, f$p)
    expect_equal(unname(attr(w, "lassosum_selection")["mode"]), "heldOut")
    # The held-out splits come from R's random number generator.
    set.seed(1)
    expect_identical(lassosumRssWeights(f$stat, f$ldReference), w)
    expect_length(
        lassosumRssWeights(f$stat, f$ldReference, selection = "minFbeta"),
        f$p
    )
})

test_that("mrashRssWeights attaches the mr.ash.rss fit on request", {
    f <- .rrwStatLd()
    # At mr.ash.rss's own tol (1e-8) this fit takes some 10,000 iterations;
    # what is retained does not depend on how tightly it converged.
    opts <- MrashRssOptions(tol = 1e-4)
    none <- mrashRssWeights(f$stat, f$ldReference, methodArgs = opts)
    expect_null(attr(none, "fit"))
    slim <- mrashRssWeights(
        f$stat,
        f$ldReference,
        fitRetention = "slim",
        methodArgs = opts
    )
    expect_equal(as.numeric(slim), as.numeric(none))
    expect_equal(as.numeric(attr(slim, "fit")$mu1), as.numeric(none))
})

test_that("the raw RSS solvers are internal behind their weight adapters", {
    # The adapters are the public route; the Options constructors still read
    # the now-internal solvers' formals, so they keep rejecting a stray name.
    exported <- getNamespaceExports("pecotmr")
    solvers <- c("prsCs", "sdpr", "lassosumRss", "penalizedRss")
    adapters <- c(
        "prsCsWeights",
        "sdprWeights",
        "lassosumRssWeights",
        "enetRssWeights",
        "ridgeRssWeights",
        "scadRssWeights",
        "mcpRssWeights",
        "l0learnRssWeights"
    )
    expect_false(any(is_in(solvers, exported)))
    expect_true(all(is_in(adapters, exported)))
    expect_error(PrsCsOptions(notAnOption = 1), "unknown argument")
    expect_error(SdprOptions(notAnOption = 1), "unknown argument")
    expect_error(LassosumOptions(notAnOption = 1), "unknown argument")
    expect_error(PenalizedRssOptions(notAnOption = 1), "unknown argument")
})

test_that("PenalizedRssOptions refuses the settings its wrappers own", {
    for (nm in c(
        "penalty",
        "alpha",
        "gamma",
        "lambda0",
        "lambda2",
        "maxSwaps"
    )) {
        ownedArg <- set_names(list(1), nm)
        expect_error(
            exec(PenalizedRssOptions, !!!ownedArg),
            str_c(nm, " is supplied by pecotmr"),
            label = nm
        )
    }
    # The rest stay tunable: `lambda` is the only route to the path for the
    # wrappers that take none.
    expect_no_error(
        PenalizedRssOptions(lambda = c(0.1, 0.01), thr = 1e-6, maxiter = 500)
    )
})

test_that("GlmnetOptions points alpha at the wrapper that takes it", {
    expect_error(
        GlmnetOptions(alpha = 0.3),
        "alpha is supplied by pecotmr.*enetWeights' own `alpha`"
    )
})

test_that("prsCsWeights and sdprWeights attach the sampler's estimates", {
    # What the samplers estimate besides the weights reaches the caller as
    # the "fit" attribute; asking for it leaves the weights unchanged.
    f <- .rrwStatLd()
    prsOpts <- PrsCsOptions(nIter = 100, nBurnin = 20, thin = 1, seed = 42L)
    prsNone <- prsCsWeights(f$stat, f$ldReference, methodArgs = prsOpts)
    expect_null(attr(prsNone, "fit"))
    prsSlim <- prsCsWeights(
        f$stat,
        f$ldReference,
        fitRetention = "slim",
        methodArgs = prsOpts
    )
    prsFit <- attr(prsSlim, "fit")
    expect_named(
        prsFit,
        c("betaEst", "psiEst", "sigmaEst", "phiEst", "variantIds")
    )
    expect_equal(prsFit$variantIds, f$stat$variantNames)
    expect_equal(as.numeric(prsSlim), as.numeric(prsNone))
    expect_equal(prsFit$betaEst, as.numeric(prsNone))
    expect_length(prsFit$psiEst, f$p)
    sdprOpts <- SdprOptions(
        iter = 100,
        burn = 20,
        thin = 1,
        verbose = FALSE,
        seed = 42L
    )
    sdprNone <- sdprWeights(f$stat, f$ldReference, methodArgs = sdprOpts)
    expect_null(attr(sdprNone, "fit"))
    sdprFull <- sdprWeights(
        f$stat,
        f$ldReference,
        fitRetention = "full",
        methodArgs = sdprOpts
    )
    expect_named(attr(sdprFull, "fit"), c("betaEst", "h2", "variantIds"))
    expect_equal(attr(sdprFull, "fit")$variantIds, f$stat$variantNames)
    expect_length(attr(sdprFull, "fit")$h2, 1L)
    expect_equal(as.numeric(sdprFull), as.numeric(sdprNone))
    expect_error(
        prsCsWeights(f$stat, f$ldReference, fitRetention = "all"),
        "fitRetention"
    )
})

test_that("lassosumRssWeights keeps its search: path at slim, beta at full", {
    f <- .rrwStatLd()
    s <- c(0.5, 1)
    # One seed per call, so the three draw the same held-out splits.
    set.seed(1)
    none <- lassosumRssWeights(f$stat, f$ldReference, s = s)
    expect_null(attr(none, "fit"))
    set.seed(1)
    slim <- lassosumRssWeights(
        f$stat,
        f$ldReference,
        s = s,
        fitRetention = "slim"
    )
    expect_named(attr(slim, "fit"), c("path", "variantIds"))
    expect_equal(attr(slim, "fit")$variantIds, f$stat$variantNames)
    path <- attr(slim, "fit")$path
    expect_named(
        path,
        c("s", "lambda", "fbeta", "loss", "conv", "nparams", "selected")
    )
    # One row per (s, lambda): lassosumRss() traces 20 lambdas by default.
    expect_equal(nrow(path), length(s) * 20L)
    expect_equal(sum(path$selected), 1L)
    expect_equal(as.numeric(slim), as.numeric(none))
    # The selected row is the one the selection attribute names.
    sel <- attr(slim, "lassosum_selection")
    expect_equal(path$s[path$selected], as.numeric(sel[["s"]]))
    expect_equal(path$lambda[path$selected], as.numeric(sel[["lambda"]]))
    set.seed(1)
    full <- lassosumRssWeights(
        f$stat,
        f$ldReference,
        s = s,
        fitRetention = "full"
    )
    beta <- attr(full, "fit")$beta
    expect_equal(dim(beta), c(f$p, nrow(path)))
    expect_equal(unname(beta[, path$selected]), as.numeric(full))
    expect_equal(colSums(beta != 0), path$nparams)
})

test_that("scad, mcp and l0learn RSS weights pass fitRetention through", {
    f <- .rrwStatLd()
    cols <- c("s", "lambda", "fbeta", "loss", "conv", "nparams", "selected")
    scad <- scadRssWeights(f$stat, f$ldReference, s = 1, fitRetention = "slim")
    mcp <- mcpRssWeights(f$stat, f$ldReference, s = 1, fitRetention = "slim")
    expect_named(attr(scad, "fit")$path, cols)
    expect_named(attr(mcp, "fit")$path, cols)
    l0 <- l0learnRssWeights(
        f$stat,
        f$ldReference,
        s = 1,
        lambda0 = c(0.01, 0.1),
        fitRetention = "full"
    )
    path <- attr(l0, "fit")$path
    expect_named(path, c("s", "lambda0", cols[-1]))
    # Plain L0 traces the single lambda = 0, so one row per lambda0.
    expect_equal(path$lambda0, c(0.01, 0.1))
    expect_equal(ncol(attr(l0, "fit")$beta), 2L)
})

test_that("held-out refits solve at the full data's path values", {
    # The default paths depend on the data -- lambda_max on the
    # correlations, the L0 family's adaptive path on every fit along it --
    # so a refit left to its own would line its candidates up with the full
    # data's by position alone. It is handed the full data's values instead.
    f <- .rrwStatLd()
    real <- pecotmr:::.penalizedRssPath
    calls <- new.env(parent = emptyenv())
    local_mocked_bindings(
        .penalizedRssPath = function(bhat, ld, n, spec, lambda, ...) {
            fit <- real(bhat, ld, n, spec, lambda, ...)
            record <- list(given = lambda, used = fit$lambda, s = ld$shrink)
            calls$log <- c(calls$log, list(record))
            fit
        }
    )
    for (fn in list(lassosumRssWeights, l0learnRssWeights)) {
        calls$log <- list()
        set.seed(1)
        fn(f$stat, f$ldReference, s = c(0.5, 1))
        full <- calls$log[1:2]
        refits <- calls$log[-(1:2)]
        expect_true(all(map_lgl(full, function(x) is.null(x$given))))
        expect_length(refits, 8L)
        for (refit in refits) {
            match <- full[[which(map_dbl(full, "s") == refit$s)]]
            expect_equal(refit$given, match$used)
        }
    }
})

test_that("a misspelled strategy fails before the grid starts", {
    f <- .rrwStatLd()
    cnd <- tryCatch(
        lassosumRssWeights(
            f$stat,
            f$ldReference,
            methodArgs = LassosumOptions(strategy = "sideways")
        ),
        error = identity
    )
    expect_match(conditionMessage(cnd), "strategy")
    expect_null(cnd$parent)
})

test_that("the RSS solvers run on a factor of the LD as on the LD itself", {
    # The factor route keeps F beta rather than R_s beta, so R_s is never
    # formed; the coordinate updates, the L0 swaps and the path are the same.
    # The dense route applies the shrinkage itself too, as fitting the
    # pre-shrunk matrix would.
    set.seed(8)
    X <- matrix(rbinom(80 * 120, 2, 0.3), 80, 120)
    ldFactor <- pecotmr:::.ldFactorFromDosage(X)
    R <- crossprod(ldFactor) / (nrow(ldFactor) - 1)
    beta <- replace(numeric(120), c(5, 40, 90), c(0.3, -0.2, 0.25))
    z <- as.vector(R %*% beta) + rnorm(120, sd = 0.05)
    lambda <- exp(seq(log(0.1), log(0.001), length.out = 8))
    lambda0 <- exp(seq(log(0.05), log(5e-4), length.out = 8))
    # lassosum is the lasso at alpha = 1; below it, the elastic net and ridge.
    # L0L2 adds its fixed ridge weight lambda2, which both routes apply in
    # the shared thresholding, swaps and objective. The L0 family's path is
    # in lambda0; the last two rows take the solver's own path.
    cases <- data.frame(
        penalty = c(
            "lasso",
            "lasso",
            "lasso",
            "MCP",
            "MCP",
            "SCAD",
            "L0",
            "L0L2",
            "L0L2",
            "L0",
            "L0L2"
        ),
        alpha = c(1, 0.5, 0, 1, 0.5, 1, 1, 1, 1, 1, 1),
        lambda2 = c(0, 0, 0, 0, 0, 0, 0, 0.05, 0.5, 0, 0.05),
        ownPath = c(rep(FALSE, 9), TRUE, TRUE)
    )
    # A small screen, so gradient ranking screens at this size.
    control <- modifyList(pecotmr:::.rssPathControl(), list(screenSize = 5L))
    for (s in c(0.2, 0.9)) {
        Rs <- (1 - s) * R + s * diag(120)
        for (strategy in c("strongRule", "gradientRanked")) {
            for (k in seq_len(nrow(cases))) {
                l0 <- startsWith(cases$penalty[[k]], "L0")
                path <- if (cases$ownPath[[k]]) {
                    numeric(0)
                } else if (l0) {
                    lambda0
                } else {
                    lambda
                }
                args <- list(
                    pathR = path,
                    penaltyStr = cases$penalty[[k]],
                    strategyStr = strategy,
                    gamma = 3.7,
                    alpha = cases$alpha[[k]],
                    lambda1 = 0,
                    lambda2 = cases$lambda2[[k]],
                    thr = 1e-6,
                    maxiter = 1000L,
                    maxSwaps = 10L,
                    control = control
                )
                label <- str_c(
                    cases$penalty[[k]],
                    ", alpha = ",
                    cases$alpha[[k]],
                    ", lambda2 = ",
                    cases$lambda2[[k]],
                    ", own path = ",
                    cases$ownPath[[k]],
                    ", ",
                    strategy
                )
                dense <- exec(pecotmr:::penalizedRssCpp, z, R, s, !!!args)
                expect_equal(
                    dense,
                    exec(pecotmr:::penalizedRssCpp, z, Rs, 0, !!!args),
                    tolerance = 1e-10,
                    label = label
                )
                expect_equal(
                    exec(
                        pecotmr:::penalizedRssFactorCpp,
                        z,
                        ldFactor,
                        s,
                        !!!args
                    ),
                    dense,
                    tolerance = 1e-8,
                    label = label
                )
            }
        }
    }
    # The two routes agreeing would not show lambda2 is applied at all; the
    # fit with it must differ from the one without.
    l0Fit <- function(lambda2) {
        pecotmr:::penalizedRssFactorCpp(
            z,
            ldFactor,
            0.2,
            pathR = lambda0,
            penaltyStr = "L0L2",
            strategyStr = "gradientRanked",
            gamma = 3.7,
            alpha = 1,
            lambda1 = 0,
            lambda2 = lambda2,
            thr = 1e-6,
            maxiter = 1000L,
            maxSwaps = 10L,
            control = control
        )$beta
    }
    expect_false(isTRUE(all.equal(l0Fit(0.5), l0Fit(0))))
})

test_that("SDPR on a factor of the LD draws as it does on the matrix", {
    # The factor route holds SDPR's matrices in low-rank form, B = W W', and
    # runs SDPR's own sampler on them, so a seeded chain matches the dense one
    # to rounding. optLlk = 2 reaches B through the Woodbury identity while
    # its low-rank part is narrower than p (one sample size, three sizes, or
    # arrays), and otherwise forms it as the dense route does (three sizes
    # with arrays).
    set.seed(7)
    X <- matrix(rbinom(40 * 120, 2, 0.3), 40, 120)
    ldFactor <- pecotmr:::.ldFactorFromDosage(X)
    R <- crossprod(ldFactor) / (nrow(ldFactor) - 1)
    p <- ncol(ldFactor)
    effects <- replace(numeric(p), c(5, 40, 90), c(0.3, -0.25, 0.2))
    bhat <- as.vector(R %*% effects) + rnorm(p, sd = 0.02)
    sizes <- sample(c(400, 500, 600), p, replace = TRUE)
    arrays <- sample(0:2, p, replace = TRUE)
    base <- list(
        bhatR = bhat,
        n = 500L,
        perVariantSampleSize = NULL,
        array = NULL,
        a = 0.1,
        c = 1,
        M = 50L,
        a0k = 0.5,
        b0k = 0.5,
        # Every sampler step runs on each iteration, so a short chain already
        # exercises burn-in, thinning and the draws the routes must share.
        iter = 100L,
        burn = 25L,
        thin = 5L,
        optLlk = 1L,
        verbose = FALSE,
        seed = 11L,
        raoBlackwell = FALSE,
        track = FALSE
    )
    cases <- list(
        llk1 = list(),
        llk2 = list(optLlk = 2L),
        sizes = list(optLlk = 2L, perVariantSampleSize = sizes),
        arrays = list(optLlk = 2L, array = arrays),
        sizesArrays = list(
            optLlk = 2L,
            perVariantSampleSize = sizes,
            array = arrays
        )
    )
    for (nm in names(cases)) {
        caseSettings <- cases[[nm]]
        args <- list_modify(base, !!!caseSettings)
        expect_equal(
            exec(pecotmr:::sdprFactorCpp, factorR = ldFactor, !!!args),
            exec(pecotmr:::sdprCpp, LD = list(R), !!!args),
            tolerance = 1e-8,
            label = nm
        )
    }
    expect_error(
        exec(pecotmr:::sdprFactorCpp, factorR = ldFactor[, -1], !!!base),
        "The LD factor has 119 columns but there are 120 effects"
    )
})

test_that("prsCsWeights and sdprWeights hand the engine correlations", {
    # Both upstream tools fit the standardized marginal effects, so a z-only
    # stat (the QtlSumStats shape) and b with its standard error both reach
    # the engine on that scale.
    f <- .rrwStatLd()
    z <- as.numeric(f$stat$z)
    want <- z / sqrt(z^2 + f$n - 2)
    local_mocked_bindings(
        .prsCsRun = function(bhat, ld, n, ...) list(betaEst = bhat),
        .sdprRun = function(bhat, ld, n, ...) list(betaEst = bhat)
    )
    ids <- f$stat$variantNames
    zOnly <- list(z = f$stat$z, n = f$n, variantNames = ids)
    withSe <- list(
        b = 3 * f$stat$z,
        seb = rep(3, f$p),
        n = rep(f$n, f$p),
        variantNames = ids
    )
    for (stat in list(zOnly, withSe)) {
        expect_equal(prsCsWeights(stat, f$ldReference), want)
        expect_equal(sdprWeights(stat, f$ldReference), want)
    }
})

test_that("a QTL past r^2 = 0.5 reaches the engines as a correlation", {
    # z / sqrt(n) passes 1 there, and the standard error the matrix engines
    # read, sqrt((1 - bhat^2) / (n - 2)), used to come out NaN.
    data(eqtlRegionExample, envir = environment())
    X <- eqtlRegionExample$X[, 1:30]
    X <- X[, apply(X, 2L, sd) > 0]
    set.seed(17)
    Xs <- scale(X)
    y <- sqrt(0.7) * Xs[, 1L] + rnorm(nrow(Xs), sd = sqrt(0.3))
    r <- as.vector(cor(Xs, y))
    n <- nrow(Xs)
    stat <- list(
        z = r * sqrt((n - 2) / (1 - r^2)),
        n = rep(n, length(r)),
        variantNames = colnames(X)
    )
    expect_gt(max(abs(stat$z)) / sqrt(n), 1)
    ld <- LdMatrix(cor(X), nSamples = n)
    # As given and through SDPR's projection alike.
    for (a in list(NULL, 0.1)) {
        dense <- pecotmr:::.alphabetRssDense(
            stat,
            ld,
            list(shrinkCorr = 0.95, ldProjection = a),
            "qgg"
        )
        expect_equal(dense$bhat, r)
        expect_true(all(is.finite(dense$seb)))
    }
    skip_if_not_installed("qgg")
    w <- bayesCRssWeights(stat, ld, "qgg", nIter = 200, nBurnin = 50)
    expect_true(all(is.finite(w)))
})

# Statistics over the first p variants of the test genotype panel, with the
# same LD as an LdMatrix: a weight method must read the panel's LD (through
# its standardized genotypes) as it reads the matrix.
test_that("RSS weights read a genotype panel's LD as they read the matrix", {
    # The panel route gives the solvers the panel's factor and never forms
    # the LD matrix; it must land on the matrix route's weights. Held-out
    # selection draws its splits through the factor on one route and the
    # matrix's eigenvectors on the other, so a deterministic selection
    # compares the fits.
    d <- .rrwPanelStat()
    fns <- list(
        lassosum = lassosumRssWeights,
        scad = scadRssWeights,
        mcp = mcpRssWeights,
        l0learn = l0learnRssWeights
    )
    for (nm in names(fns)) {
        expect_equal(
            fns[[nm]](d$stat, d$panel, selection = "minFbeta"),
            fns[[nm]](d$stat, d$ldMatrix, selection = "minFbeta"),
            tolerance = 1e-6,
            label = nm
        )
    }
    # mr.ash.rss takes only a matrix, formed from the panel. A looser tol
    # than its own 1e-8 (some 7,000 iterations here) converges in about 100;
    # the two routes fit the same matrix whatever the tolerance.
    mrashOpts <- MrashRssOptions(tol = 1e-4)
    expect_equal(
        mrashRssWeights(d$stat, d$panel, methodArgs = mrashOpts),
        mrashRssWeights(d$stat, d$ldMatrix, methodArgs = mrashOpts),
        tolerance = 1e-8
    )
    # SDPR runs its own sampler on the low-rank form of its matrices, so a
    # seeded chain lands on the matrix route's. A 50-component mixture
    # rather than the default 1000 keeps the chains quick, and statistics
    # consistent with the LD keep them from diverging.
    d <- .rrwPanelStatLd()
    for (optLlk in 1:2) {
        opts <- SdprOptions(
            iter = 100,
            burn = 20,
            thin = 1,
            M = 50,
            optLlk = optLlk,
            verbose = FALSE,
            seed = 1L
        )
        expect_equal(
            sdprWeights(d$stat, d$panel, methodArgs = opts),
            sdprWeights(d$stat, d$ldMatrix, methodArgs = opts),
            tolerance = 1e-6,
            label = str_c("sdpr optLlk ", optLlk)
        )
    }
})

test_that("an RSS weight method matches its LD reference by variant id", {
    f <- .rrwStatLd()
    unnamed <- f$stat
    unnamed$variantNames <- NULL
    expect_error(
        lassosumRssWeights(unnamed, f$ldReference),
        "must carry `variantNames`"
    )
    expect_error(
        mrashRssWeights(unnamed, f$ldReference),
        "must carry `variantNames`"
    )
    # The extractors need the ids only to check against a reference.
    expect_equal(pecotmr:::.rssVariantCount(unnamed, NULL, "t"), f$p)
    expect_equal(pecotmr:::.rssVariantCount(f$stat, f$ldReference, "t"), f$p)
    expect_true(is.na(pecotmr:::.rssVariantCount(list(), NULL, "t")))
})

test_that(".formalDefaults evaluates a function's own defaults", {
    expect_equal(
        pecotmr:::.formalDefaults(prsCs, c("a", "phi", "nIter")),
        list(a = 1, phi = NULL, nIter = 1000)
    )
})

test_that("mrashRssWeights defaults to mr.ash's own prior grid", {
    # On the standardized scale mr.ash.rss fits, its sufficient statistics
    # are those of standardized genotypes, whose grid mr.ash builds itself.
    set.seed(1)
    X <- scale(matrix(rnorm(200 * 10), 200))
    fit <- susieR::mr.ash(X, rnorm(200), verbose = FALSE, max.iter = 5)
    expect_equal(
        pecotmr:::.mrashRssDefaultGrid(200, rnorm(10), numeric(0), NULL),
        fit$data$sa2
    )
    f <- .rrwStatLd()
    # A looser tol than mr.ash.rss's 1e-8, which takes some 10,000 iterations.
    expect_length(
        mrashRssWeights(
            f$stat,
            f$ldReference,
            methodArgs = MrashRssOptions(tol = 1e-4)
        ),
        f$p
    )
    # With b, seb and varY the unstandardized statistics set the scale.
    z <- f$stat$z
    shat <- f$stat$seb
    w <- 2 * (f$n - 1) / (z^2 + f$n - 2) / shat^2
    expect_equal(
        pecotmr:::.mrashRssDefaultGrid(f$n, z, shat, 2),
        (2^((0:24) / 25) - 1)^2 / median(w) * f$n
    )
})


test_that("mrmashRssWeights fits mr.mash.rss and returns p x K weights", {
    .skipIfMrmashCannotFitStudies()
    m <- .rrwMulti(n = 60, p = 6, K = 3)
    # Six variants are too few for mr.mash.rss to estimate the residual
    # correlation from.
    expect_warning(
        w <- mrmashRssWeights(m$stat, m$ldReference),
        "near-null variables"
    )
    expect_equal(dim(w), c(m$p, m$K))
    expect_true(all(is.finite(w)))
})

test_that("mrmashRssWeights reads a genotype panel's LD as a factor", {
    skip_if_not_installed("mr.mashr")
    skip_if_not(
        is_in("X", names(formals(mr.mashr::mr.mash.rss))),
        "mr.mashr predates an LD factor"
    )
    d <- .rrwPanelStat()
    withr::local_seed(4L)
    stat <- list_assign(
        d$stat,
        z = cbind(a = d$stat$z, b = d$stat$z + stats::rnorm(40L))
    )
    # The panel's genotypes as the factor fit as the matrix they give. These
    # statistics take mr.mash.rss past its 5000 iterations; the two forms of
    # the LD agree at every iteration, so a few dozen check it.
    opts <- MrmashRssOptions(verbose = FALSE, max_iter = 50L)
    expect_equal(
        mrmashRssWeights(
            stat,
            d$panel,
            V = diag(2),
            covY = diag(2),
            methodArgs = opts
        ),
        mrmashRssWeights(
            stat,
            d$ldMatrix,
            V = diag(2),
            covY = diag(2),
            methodArgs = opts
        ),
        tolerance = 1e-6
    )
    ld <- pecotmr:::.rssConditionLd(stat, d$panel, "mrmashRssWeights")
    expect_equal(ncol(ld$X), 40L)
    expect_null(ld$R)
    # One reference per condition, the second lacking the last variant: a
    # zero column of its factor, which the statistics' carrying reads.
    partial <- d$panel[seq_len(39L), ]
    perCondition <- pecotmr:::.rssConditionLd(
        stat,
        list(d$panel, partial),
        "mrmashRssWeights"
    )
    expect_equal(unname(perCondition$X[[2L]][, 40L]), rep(0, ncol(d$panel)))
    lacked <- pecotmr:::.rssLackedVariants(perCondition)
    expect_equal(which(lacked), 80L)
    expect_equal(
        unname(pecotmr:::.rssCarriedStat(stat, perCondition)$z[40L, 2L]),
        0
    )
})

test_that("mrmashRssWeights errors on single-context stat$z", {
    skip_if_not_installed("mr.mashr")
    f <- .rrwStatLd()
    oneCol <- list(
        z = matrix(f$stat$z, ncol = 1),
        n = f$n,
        variantNames = f$stat$variantNames
    )
    expect_error(mrmashRssWeights(oneCol, f$ldReference), ">= 2 columns")
})

test_that("mrmashRssWeights takes one LD reference per condition", {
    skip_if_not_installed("mr.mashr")
    skip_if_not(
        is_in("missing_y_method", names(formals(mr.mashr::mr.mash))),
        "mr.mashr predates missing_y_method"
    )
    m <- .rrwMulti(n = 60, p = 6, K = 3)
    opts <- MrmashRssOptions(verbose = FALSE)
    # Six variants are too few to estimate the residual correlation from, so
    # it is given.
    fitOn <- function(stat, ldReference) {
        mrmashRssWeights(
            stat,
            ldReference,
            V = diag(3),
            covY = diag(3),
            methodArgs = opts
        )
    }
    # The same reference for every condition is the one-reference fit.
    expect_equal(
        fitOn(m$stat, rep(list(m$ldReference), 3L)),
        fitOn(m$stat, m$ldReference),
        tolerance = 1e-8
    )
    # Two studies on their own references: the third condition's study lacks
    # the last variant, which carries no information there.
    other <- LdMatrix(`dimnames<-`(
        0.5^abs(outer(1:6, 1:6, "-")),
        list(m$stat$variantNames, m$stat$variantNames)
    ))
    stat <- m$stat
    stat$z[6L, 3L] <- NA
    refs <- list(m$ldReference, m$ldReference, other)
    w <- fitOn(stat, refs)
    expect_equal(dim(w), c(m$p, m$K))
    expect_true(all(is.finite(w)))
})

test_that("one LD reference per condition is laid on the variants", {
    m <- .rrwMulti(n = 60, p = 6, K = 3)
    other <- LdMatrix(`dimnames<-`(
        0.5^abs(outer(1:6, 1:6, "-")),
        list(m$stat$variantNames, m$stat$variantNames)
    ))
    # The third condition's study lacks the last variant.
    stat <- m$stat
    stat$z[6L, 3L] <- NA
    refs <- list(m$ldReference, m$ldReference, other)
    ld <- pecotmr:::.rssConditionLd(stat, refs, "mrmashRssWeights")
    expect_identical(ld$R[[1L]], ld$R[[2L]])
    expect_equal(unname(ld$R[[3L]][6L, ]), rep(0, 6L))
    expect_equal(unname(ld$R[[3L]][1:5, 1:5]), 0.5^abs(outer(1:5, 1:5, "-")))
    expect_equal(ld$n, rep(60, 3L))
    expect_error(
        pecotmr:::.rssConditionLd(stat, refs[1:2], "mrmashRssWeights"),
        "one LD reference per column of stat\\$z \\(3\\); got 2"
    )
    # A statistic a condition lacks carries no information: z and Bhat 0,
    # Shat 1, where its LD is zero.
    expect_equal(
        unname(pecotmr:::.rssCarriedStat(stat, ld)$z[6L, 3L]),
        0
    )
    withSe <- pecotmr:::.rssCarriedStat(
        list(Bhat = matrix(c(NA, 1), 2L), Shat = matrix(c(NA, 2), 2L)),
        list(R = list(diag(2)))
    )
    expect_equal(withSe$Bhat[, 1L], c(0, 1))
    expect_equal(withSe$Shat[, 1L], c(1, 2))
    # On one reference, a missing statistic carries no information either;
    # the z-scores the fit reads keep it missing (.rssFitZ()).
    oneRef <- list(R = diag(6))
    expect_equal(
        unname(pecotmr:::.rssCarriedStat(stat, oneRef)$z[6L, 3L]),
        0
    )
    expect_true(is.na(pecotmr:::.rssFitZ(stat, oneRef)[6L, 3L]))
})

test_that("a variant its reference lacks carries no information there", {
    m <- .rrwMulti(n = 60, p = 6, K = 3)
    ids <- m$stat$variantNames
    # The third condition's reference holds only the first five variants,
    # though that condition has a z-score at the sixth.
    partial <- LdMatrix(`dimnames<-`(
        0.5^abs(outer(1:5, 1:5, "-")),
        list(ids[1:5], ids[1:5])
    ))
    refs <- list(m$ldReference, m$ldReference, partial)
    ld <- pecotmr:::.rssConditionLd(m$stat, refs, "mrmashRssWeights")
    expect_equal(unname(diag(ld$R[[3L]])), c(rep(1, 5L), 0))
    carried <- pecotmr:::.rssCarriedStat(m$stat, ld)
    expect_equal(unname(carried$z[6L, 3L]), 0)
    expect_equal(unname(carried$z[6L, 1:2]), unname(m$stat$z[6L, 1:2]))
    expect_equal(
        pecotmr:::.rssLackedVariants(ld),
        cbind(rep(FALSE, 6L), rep(FALSE, 6L), c(rep(FALSE, 5L), TRUE))
    )
})

test_that("mr.mash.rss takes the mash prior's residual correlation", {
    .skipIfMrmashCannotFitStudies()
    m <- .rrwMulti(n = 60, p = 6, K = 3)
    R <- 0.3 + diag(0.7, 3)
    prior <- list(U = list(shared = matrix(1, 3L, 3L)), residualCorrelation = R)
    oneRef <- list(R = diag(6))
    out <- pecotmr:::.mrmashRssResidual(m$stat, oneRef, prior, NULL, NULL)
    # covY, whose correlations V takes; without h2, mr.mash.rss starts V
    # from covY itself.
    expect_equal(out$covY, R)
    expect_null(out$V)
    expect_false(out$originalScale)
    # Zero between conditions on different references, fitted as different
    # samples.
    perRef <- list(R = rep(list(diag(6)), 3L), group = c(1L, 1L, 2L))
    masked <- pecotmr:::.mrmashRssResidual(m$stat, perRef, prior, NULL, NULL)
    expect_equal(masked$covY[1:2, 3L], c(0, 0))
    expect_equal(masked$covY[1L, 2L], 0.3)
    # Each supplied one is used as it is.
    own <- pecotmr:::.mrmashRssResidual(m$stat, oneRef, prior, diag(3), NULL)
    expect_equal(own$V, diag(3))
    expect_equal(own$covY, R)
    expect_error(
        pecotmr:::.mrmashRssResidual(
            m$stat,
            oneRef,
            list(residualCorrelation = diag(2)),
            NULL,
            NULL
        ),
        "residualCorrelation"
    )
})

test_that("without a mash prior, mr.mash.rss estimates the correlation", {
    .skipIfMrmashCannotFitStudies()
    set.seed(7)
    n <- 300L
    p <- 120L
    X <- matrix(
        rnorm(n * p),
        n,
        p,
        dimnames = list(NULL, sprintf("chr1:%d:A:G", 100L * seq_len(p)))
    )
    E <- matrix(rnorm(n * 2L), n, 2L) %*% chol(0.5 + diag(0.5, 2L))
    Y <- cbind(a = 0.5 * X[, 1L] + E[, 1L], b = 0.4 * X[, 1L] + E[, 2L])
    r <- cor(X, Y)
    stat <- list(
        z = sqrt(n - 2) * r / sqrt(1 - r^2),
        n = n,
        variantNames = colnames(X)
    )
    # Nothing given: mr.mash.rss estimates both from the near-null z-scores
    out <- pecotmr:::.mrmashRssResidual(
        stat,
        list(R = cor(X)),
        NULL,
        NULL,
        NULL
    )
    expect_null(out$covY)
    expect_null(out$V)
    fitV <- function(stat) {
        w <- mrmashRssWeights(
            stat,
            LdMatrix(cor(X)),
            fitRetention = "slim",
            methodArgs = MrmashRssOptions(verbose = FALSE)
        )
        unname(attr(w, "fit")$V)
    }
    zAdj <- sqrt((n - 1) / (stat$z^2 + n - 2)) * stat$z
    nearNull <- unname(mashr::estimate_null_correlation_simple(
        mashr::mash_set_data(zAdj)
    ))
    expect_equal(fitV(stat), nearNull)
    # The conditions' heritabilities give V their 1 - h2, with those
    # correlations
    s <- sqrt(c(0.7, 1))
    expect_equal(fitV(c(stat, list(h2 = c(0.3, NA)))), nearNull * outer(s, s))
    # A z-score missing in one condition is left out of the estimate (a 0
    # would pull the correlation toward 0), without a warning
    gappy <- stat
    gappy$z[3:12, "b"] <- NA
    expect_silent(gappyV <- fitV(gappy))
    expect_equal(
        gappyV,
        unname(mashr::estimate_null_correlation_simple(
            mashr::mash_set_data(zAdj[-(3:12), ])
        ))
    )
})

test_that("the RSS fits read a missing z-score as missing, a lacked one as 0", {
    z <- cbind(a = c(1, NA, 2), b = c(NA, 1, 3))
    stat <- list(z = z)
    # One LD for both conditions: missing stays missing
    expect_identical(pecotmr:::.rssFitZ(stat, list(R = diag(3))), z)
    # One LD per condition: a variant the reference of b lacks is 0
    perCondition <- list(R = list(diag(3), diag(c(0, 1, 1))))
    expect_equal(
        pecotmr:::.rssFitZ(stat, perCondition),
        cbind(a = c(1, NA, 2), b = c(0, 1, 3))
    )
    # The statistics the prior grid and h2 read carry no missing values
    expect_equal(
        pecotmr:::.rssCarriedStat(stat, list(R = diag(3)))$z,
        replace(z, is.na(z), 0)
    )
})

test_that(".mrmashRssStats reads sample sizes per variant, condition or both", {
    z <- matrix(c(1, -2, 0.5, 3, 1, -1), 3L, 2L)
    perVariant <- pecotmr:::.mrmashRssStats(list(z = z, n = c(100, 400, 900)))
    expect_equal(perVariant$Bhat, z / c(10, 20, 30))
    expect_equal(perVariant$Shat, matrix(1 / c(10, 20, 30), 3L, 2L))
    expect_equal(perVariant$nVec, c(400, 400))
    perCondition <- pecotmr:::.mrmashRssStats(list(z = z, n = c(100, 400)))
    expect_equal(
        perCondition$Shat,
        matrix(1 / c(10, 20), 3L, 2L, byrow = TRUE)
    )
    expect_equal(perCondition$nVec, c(100, 400))
    both <- matrix(c(100, 200, 300, 400, 500, 600), 3L, 2L)
    expect_equal(
        pecotmr:::.mrmashRssStats(list(z = z, n = both))$nVec,
        c(200, 500)
    )
    expect_error(
        pecotmr:::.mrmashRssStats(list(z = z, n = 1:4)),
        "one per variant"
    )
})

test_that("mr.mash.rss fits z-scores on the standardized scale", {
    .skipIfNoMrmashRelativeS0()
    m <- .rrwMulti(n = 60, p = 6, K = 3)
    prior <- pecotmr:::.mrmashRssDefaults(
        pecotmr:::.mrmashRssStats(m$stat),
        NULL,
        TRUE,
        NULL,
        NULL
    )
    w <- mrmashRssWeights(
        m$stat,
        m$ldReference,
        S0 = prior$S0,
        w0 = prior$w0,
        covY = diag(3),
        methodArgs = MrmashRssOptions(verbose = FALSE)
    )
    direct <- mr.mashr::mr.mash.rss(
        Z = m$stat$z,
        R = m$LD,
        n = 60,
        S0 = prior$S0,
        w0 = prior$w0,
        covY = diag(3),
        relative_S0 = TRUE,
        verbose = FALSE,
        nthreads = 1L
    )
    expect_equal(unname(w[,]), unname(mr.mashr::coef.mr.mash.rss(direct)))
    # With Bhat, Shat and covY, the fit is on their scale, where h2 leaves
    # 1 - h2 of each condition's variance
    scaled <- pecotmr:::.mrmashRssResidual(
        c(
            m$stat,
            list(Bhat = m$stat$z, Shat = 1 + 0 * m$stat$z, h2 = c(0.5, NA, 0.2))
        ),
        list(R = m$LD),
        NULL,
        NULL,
        diag(c(2, 3, 4))
    )
    expect_true(scaled$originalScale)
    expect_equal(scaled$V, c(1, 3, 3.2))
})

test_that("the mash prior's residual correlation reaches mr.mash.rss", {
    .skipIfMrmashCannotFitStudies()
    m <- .rrwMulti(n = 60, p = 6, K = 3)
    R <- 0.3 + diag(0.7, 3)
    w <- mrmashRssWeights(
        m$stat,
        m$ldReference,
        dataDrivenPriorMatrices = list(
            U = list(shared = matrix(1, 3L, 3L)),
            residualCorrelation = R
        ),
        fitRetention = "slim",
        methodArgs = MrmashRssOptions(verbose = FALSE)
    )
    # mr.mash.rss does not update V, so the fit's V is the one it was given,
    # and the payload keeps the correlation with the prior.
    expect_equal(unname(attr(w, "fit")$V), R)
    expect_equal(attr(w, "fit")$dataDrivenPriorMatrices$residualCorrelation, R)
})

test_that("conditions share a reference group only with an identical one", {
    a <- .rrwLdMatrix(diag(2))
    b <- .rrwLdMatrix(0.5 + diag(0.5, 2))
    expect_equal(pecotmr:::.rssRefGroups(list(a, b, a, b)), c(1, 2, 1, 2))
    expect_equal(
        pecotmr:::.rssMaskBetweenRefs(matrix(0.5, 3L, 3L), c(1L, 2L, 1L)),
        matrix(c(0.5, 0, 0.5, 0, 0.5, 0, 0.5, 0, 0.5), 3L, 3L)
    )
    expect_identical(pecotmr:::.rssMaskBetweenRefs(diag(2), NULL), diag(2))
    # A reference holding none of the variants contributes zero LD.
    expect_equal(
        pecotmr:::.rssEmbedLd(NULL, c("v1", "v2")),
        matrix(0, 2L, 2L, dimnames = list(c("v1", "v2"), c("v1", "v2")))
    )
})


# ------------------------------- pure helpers -------------------------------

test_that(".rssMarginalCor reads cor / z / b and validates length", {
    f <- .rrwStatLd()
    expect_length(
        pecotmr:::.rssMarginalCor(
            list(cor = f$stat$cor),
            n = f$n,
            p = f$p
        ),
        f$p
    )
    expect_equal(
        pecotmr:::.rssMarginalCor(list(z = f$stat$z), n = f$n, p = f$p),
        as.numeric(f$stat$z) / sqrt(as.numeric(f$stat$z)^2 + f$n - 2)
    )
    expect_equal(
        pecotmr:::.rssMarginalCor(list(b = f$stat$b), n = f$n, p = f$p),
        as.numeric(f$stat$b)
    )
    # With its standard error, b is a z-score on any scale.
    expect_equal(
        pecotmr:::.rssMarginalCor(
            list(b = 2 * f$stat$z, seb = rep(2, f$p)),
            n = f$n,
            p = f$p
        ),
        as.numeric(f$stat$z) / sqrt(as.numeric(f$stat$z)^2 + f$n - 2)
    )
    # The correlation behind a t statistic, exactly, however strong: reading
    # z / sqrt(n) instead would put it past 1.
    r <- c(-0.9, 0.1, 0.8)
    n <- 50
    expect_equal(
        pecotmr:::.rssMarginalCor(
            list(z = r * sqrt((n - 2) / (1 - r^2))),
            n = n,
            p = 3L
        ),
        r
    )
    expect_error(
        pecotmr:::.rssMarginalCor(list(), n = f$n, p = f$p),
        "one of"
    )
    expect_error(
        pecotmr:::.rssMarginalCor(
            list(z = f$stat$z[-1]),
            n = f$n,
            p = f$p
        ),
        "must equal"
    )
})

test_that(".lassosumClampCor scales values with |cor| >= 1 below 1", {
    expect_equal(pecotmr:::.lassosumClampCor(c(0.1, 0.5)), c(0.1, 0.5))
    expect_lt(max(abs(pecotmr:::.lassosumClampCor(c(0.5, 1.5, -2)))), 1)
})

test_that(".lassosumFirstMax returns the first index of the maximum", {
    expect_equal(pecotmr:::.lassosumFirstMax(c(1, 3, 3, 2)), 2L)
    expect_equal(pecotmr:::.lassosumFirstMax(c(5, 1, 2)), 1L)
    # A tie to rounding is a tie; a real difference is not.
    expect_equal(pecotmr:::.lassosumFirstMax(c(1, 3, 3 + 1e-14, 2)), 2L)
    expect_equal(pecotmr:::.lassosumFirstMax(c(1, 3, 3 + 1e-6, 2)), 3L)
    expect_equal(pecotmr:::.lassosumFirstMax(c(-Inf, -Inf)), 1L)
})

test_that("held-out selection picks on variance, then on correlation", {
    pick <- pecotmr:::.rssHeldOutPick
    # A real difference in the variance wins whatever the correlation.
    expect_equal(pick(c(0.9, 0.5, 0.6), c(-3, 0.2, 0.4)), 3L)
    expect_equal(pick(c(0.9, 0.5, 0.6), c(0.5, 0.2, 0.4)), 1L)
    # Candidates 2 to 4 tie on the variance, to rounding; of them, 3 has the
    # best correlation.
    expect_equal(pick(c(0.9, 0.5, 0.7, 0.6), c(-1, 0.3, 0.3 + 1e-15, 0.3)), 3L)
    # Tied on both, the first in the grid.
    expect_equal(pick(c(0.4, 0.5, 0.5), c(0.1, 0.2, 0.2)), 2L)
})

test_that("held-out selection passes over an empty fit while another exists", {
    pick <- pecotmr:::.rssHeldOutPick
    # An empty fit explains nothing, which beats every fit that overshoots.
    expect_equal(pick(c(-Inf, 0.2, 0.3), c(0, -0.1, -0.2)), 2L)
    # With nothing else left, the first in the grid.
    expect_equal(pick(c(-Inf, -Inf), c(0, 0)), 1L)
    # Empty on the full data but the best in every split: the fit returned
    # would be empty, so the next one is. Candidate 1's splits fit the
    # statistics themselves; candidate 2's are the same, heavily shrunk.
    corInput <- c(0.3, 0.1, -0.1, 0.2)
    refit <- list(
        fn = function(trainCor, share) cbind(trainCor, 0.01 * trainCor),
        args = list()
    )
    sel <- function(candidates) {
        set.seed(1)
        pecotmr:::.rssSelectHeldOut(
            candidates,
            corInput,
            list(R = diag(4), nSamples = 1000L),
            100,
            refit
        )$index
    }
    expect_equal(sel(cbind(corInput, 0.01 * corInput)), 1L)
    expect_equal(sel(cbind(numeric(4), 0.01 * corInput)), 2L)
})

test_that("held-out selection refuses a fit far larger than the trait", {
    # Real genotypes, a window of strong LD and one causal variant with many
    # tags. The fit at s = 1 scores best on correlation, and a pick on
    # correlation takes it, but its predicted variance is about 10,000.
    data(eqtlRegionExample, envir = environment())
    G <- eqtlRegionExample$X
    af <- colMeans(G) / 2
    G <- G[, pmin(af, 1 - af) > 0.05][, 1301:1500]
    Gs <- scale(G)
    h2 <- 0.6
    set.seed(1)
    g <- as.numeric(Gs %*% replace(numeric(200), 1, 1))
    y <- g * sqrt(h2 / var(g)) + rnorm(nrow(G), sd = sqrt(1 - h2))
    study <- sample(nrow(G), 300L)
    R <- cor(G[study, ])
    stat <- list(
        cor = as.numeric(cor(G[study, ], y[study])),
        n = rep(300, 200),
        variantNames = colnames(G)
    )
    w <- lassosumRssWeights(
        stat,
        LdMatrix(R, nSamples = 300L),
        s = c(0.2, 1),
        fitRetention = "full"
    )
    predicted <- function(b) sum(b * (R %*% b))
    fit <- attr(w, "fit")
    # The least penalized fit at s = 1, chosen by value, not grid position.
    sOne <- which(fit$path$s == 1)
    atOne <- fit$beta[, sOne[which.min(fit$path$lambda[sOne])]]
    expect_gt(predicted(atOne), 100 * h2)
    expect_lt(predicted(as.numeric(w)), 2 * h2)
    expect_equal(
        as.numeric(attr(w, "lassosum_selection")[["s"]]),
        0.2
    )
})

test_that(".lassosumSelectMinFbeta picks the minimum-fbeta candidate", {
    set.seed(7)
    cb <- matrix(rnorm(6 * 4), 6, 4)
    r <- pecotmr:::.lassosumSelectMinFbeta(
        cb,
        data.frame(fbeta = c(3, 1, 2, 4))
    )
    expect_equal(r$index, 2)
    expect_equal(r$mode, "minFbeta")
    expect_equal(r$beta, cb[, 2])
})

test_that(".lassosumSelectLdQuadratic shrinks c by its lfdr like lassosum", {
    skip_if_not_installed("fdrtool")
    # Five strong correlations among 300; candidate 1 fits those five,
    # candidate 2 is the least-squares direction c itself (R = I), which
    # the raw in-sample score always prefers.
    set.seed(8)
    p <- 300L
    corInput <- c(rnorm(p - 5L, sd = 0.03), 0.3, -0.25, 0.2, 0.28, -0.3)
    signal <- replace(numeric(p), (p - 4L):p, sign(corInput[(p - 4L):p]))
    cb <- unname(cbind(signal, corInput))
    lfdr <- fdrtool::fdrtool(
        corInput,
        statistic = "correlation",
        plot = FALSE,
        verbose = FALSE
    )$lfdr
    score <- function(cVec) colSums(cVec * cb) / sqrt(colSums(cb^2))
    expect_equal(which.max(score(corInput)), 2L)
    want <- which.max(score(corInput * (1 - lfdr)))
    expect_equal(want, 1L)
    r <- pecotmr:::.lassosumSelectLdQuadratic(cb, corInput, list(R = diag(p)))
    expect_equal(r$mode, "ldQuadratic")
    expect_equal(r$index, want)
    expect_equal(r$beta, cb[, want])
})

test_that(".rssShrunkLdScore scores c'b against LD shrunk toward I", {
    R <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
    cb <- cbind(c(1, 1), c(1, -1), c(0, 0))
    corInput <- c(0.2, 0.1)
    w <- 0.3
    want <- colSums(corInput * cb[, 1:2]) /
        sqrt((1 - w) * colSums(cb[, 1:2] * (R %*% cb[, 1:2])) + w * 2)
    expect_equal(
        pecotmr:::.rssShrunkLdScore(cb, corInput, list(R = R), w),
        c(want, -Inf)
    )
})

test_that("held-out selection needs the LD reference's sample size", {
    f <- .rrwStatLd()
    unsized <- LdMatrix(f$LD)
    expect_error(
        lassosumRssWeights(f$stat, unsized),
        "lassosumRssWeights: selection = \"heldOut\" needs the LD reference's"
    )
    expect_error(
        mcpRssWeights(f$stat, unsized),
        "mcpRssWeights: selection = \"heldOut\" needs the LD reference's"
    )
    expect_length(
        lassosumRssWeights(f$stat, unsized, selection = "minFbeta"),
        f$p
    )
})

test_that("held-out selection draws its splits through a panel's factor", {
    d <- .rrwPanelStat()
    set.seed(2)
    w <- lassosumRssWeights(d$stat, d$panel)
    expect_equal(unname(attr(w, "lassosum_selection")["mode"]), "heldOut")
    expect_length(w, length(d$stat$z))
    expect_true(all(is.finite(w)))
})

test_that("held-out selection avoids the fit a small reference cannot judge", {
    # Real genotypes: statistics from 250 people, LD from 100 others, over
    # 200 variants. The least-penalized fit loads on directions of the LD
    # the 100-person reference cannot see, and lassosum's in-sample score
    # keeps choosing it. Held-out selection, with the LD shrunk for the
    # reference's size, must pick a sparser fit that predicts the remaining
    # people better.
    data(eqtlRegionExample, envir = environment())
    G <- eqtlRegionExample$X
    af <- colMeans(G) / 2
    G <- G[, pmin(af, 1 - af) > 0.05][, 1001:1200]
    Gs <- scale(G)
    set.seed(11)
    g <- as.numeric(Gs %*% replace(numeric(200), c(40, 150), c(1, -0.8)))
    y <- g * sqrt(0.15 / var(g)) + rnorm(nrow(G), sd = sqrt(0.85))
    perm <- sample(nrow(G))
    study <- perm[1:250]
    test <- perm[351:415]
    R <- cor(G[perm[251:350], ])
    stat <- list(
        cor = as.numeric(cor(G[study, ], y[study])),
        n = rep(250, 200),
        variantNames = colnames(G)
    )
    w <- lassosumRssWeights(
        stat,
        LdMatrix(R, nSamples = 100L),
        fitRetention = "full"
    )
    fit <- attr(w, "fit")
    loosest <- which(
        fit$path$s == 0.2 & fit$path$lambda == min(fit$path$lambda)
    )
    accuracy <- function(b) cor(as.numeric(Gs[test, ] %*% b), g[test])
    expect_lt(sum(w != 0), fit$path$nparams[loosest])
    expect_gt(accuracy(as.numeric(w)), accuracy(fit$beta[, loosest]))
})

test_that("computeCovDiag returns a diagonal condition covariance", {
    m <- .rrwMulti(n = 60, p = 6, K = 3)
    cv <- computeCovDiag(m$Y)
    expect_equal(dim(cv), c(m$K, m$K))
    expect_equal(cv[upper.tri(cv)], rep(0, sum(upper.tri(cv))))
    expect_equal(unname(diag(cv)), unname(apply(m$Y, 2, var)))
})

test_that("buildMrmashPriorMatrices builds an expanded S0 and a prior grid", {
    skip_if_not_installed("mr.mashr")
    set.seed(9)
    res <- buildMrmashPriorMatrices(
        Bhat = matrix(rnorm(18), 6, 3),
        Shat = matrix(0.2, 6, 3),
        K = 3
    )
    expect_true(is.list(res$S0))
    expect_gt(length(res$S0), 1)
    expect_true(is.numeric(res$priorGrid))
    expect_true(all(map_lgl(
        res$S0,
        function(s) all(dim(s) == c(3, 3))
    )))
})

test_that("buildMrmashPriorMatrices needs a canonical or data-driven prior", {
    skip_if_not_installed("mr.mashr")
    expect_error(
        buildMrmashPriorMatrices(
            Bhat = matrix(rnorm(6), 3, 2),
            Shat = matrix(0.2, 3, 2),
            K = 2,
            canonicalPriorMatrices = FALSE
        ),
        "dataDrivenPriorMatrices"
    )
})


# ===========================================================================
# mrmashWeights and the mr.mash engine arguments
# ===========================================================================

# ---- mrmashWeights ----
test_that("mrmashWeights errors when mr.mashr package is not available", {
    skip_if(
        requireNamespace("mr.mashr", quietly = TRUE),
        "mr.mashr is installed; skipping missing-package test"
    )

    expect_error(
        mrmashWeights(
            mrmashFit = NULL,
            X = matrix(1, 10, 5),
            Y = matrix(1, 10, 3)
        ),
        "mr\\.mash\\.alpha"
    )
})

test_that("mrmashWeights errors when X and Y are NULL and fit is NULL", {
    skip_if_not(
        requireNamespace("mr.mashr", quietly = TRUE),
        "mr.mashr not installed"
    )
    expect_error(
        mrmashWeights(mrmashFit = NULL, X = NULL, Y = NULL),
        "Both X and Y must be provided"
    )
})

test_that("MCMC / optimiser control arguments are guarded", {
    expect_error(
        bayesBWeights(NULL, NULL, nIter = 0),
        "nIter.*Must be >= 1"
    )
    expect_error(
        bayesBWeights(NULL, NULL, thin = 0),
        "thin.*Must be >= 1"
    )
    expect_error(
        bayesBWeights(NULL, NULL, pi = 2),
        "pi.*is not <= 1"
    )
    expect_error(
        bayesLWeights(NULL, NULL, nBurnin = -1),
        "nBurnin.*Must be >= 0"
    )
    expect_error(
        dprGibbsWeights(NULL, NULL, sStep = 0),
        "sStep.*Must be >= 1"
    )
    expect_error(
        dprAdaptiveGibbsWeights(NULL, NULL, fitRetention = "sometimes"),
        "must be one of"
    )
    expect_error(
        ncvregWeights(NULL, NULL, penalty = "SCAD", nfolds = 0),
        "nfolds.*Must be >= 1"
    )
    expect_error(
        computeCoefficientsGlasso(NULL, NULL, standardize = NA),
        "standardize.*May not be NA"
    )
})

test_that("RSS solver control arguments are guarded", {
    expect_error(
        lassosumRss(NULL, NULL, NULL, thr = -1),
        "thr.*is not >= 0"
    )
    expect_error(
        lassosumRss(NULL, NULL, NULL, maxiter = 0),
        "maxiter.*Must be >= 1"
    )
    expect_error(
        scadRssWeights(stat = "nope", ldReference = NULL),
        "stat.*Must be of type 'list'"
    )
    expect_error(
        scadRssWeights(stat = list(), ldReference = NULL, s = -1),
        "s.*is not >= 0"
    )
    expect_error(
        mcpRssWeights(stat = "nope", ldReference = NULL),
        "stat.*Must be of type 'list'"
    )
    expect_error(
        mrashRssWeights(stat = "nope", ldReference = NULL),
        "stat.*Must be of type 'list'"
    )
    expect_error(
        mrmashRssWeights(stat = "nope", ldReference = NULL),
        "stat.*Must be of type 'list'"
    )
})

test_that("MrmashPriorParam carries prior options and rejects a typo", {
    pa <- MrmashPriorParam(canonicalPriorMatrices = TRUE)
    expect_s4_class(pa, "MrmashPriorParam")
    expect_true(pa$canonicalPriorMatrices)
    # pecotmr's own builder, so no `...`: R rejects an unknown name. Note
    # `canonicalPrior` would be ACCEPTED, by R's partial matching.
    expect_error(
        MrmashPriorParam(canonicalPriorMatrix = TRUE),
        "unused argument"
    )
    expect_true(MrmashPriorParam(canonicalPrior = TRUE)$canonicalPriorMatrices)
})

test_that("MrmashPriorParam nests each engine function's own arguments", {
    skip_if_not_installed("mr.mashr")
    # hetgrid and singletons belong to compute_canonical_covs, zeromat to
    # expand_covs. Nesting them says where each goes AND gets each checked
    # against that function's live formals rather than a transcribed list.
    pa <- MrmashPriorParam(
        canonicalCovs = MrmashCanonicalCovsOptions(
            hetgrid = c(0, 0.5),
            singletons = FALSE
        ),
        expandCovs = MrmashExpandCovsOptions(zeromat = FALSE)
    )
    expect_equal(pa$canonicalCovs$hetgrid, c(0, 0.5))
    expect_false(pa$canonicalCovs$singletons)
    expect_false(pa$expandCovs$zeromat)
    # Each bundle rejects the other's name, which a single flat bag could not.
    expect_error(
        MrmashCanonicalCovsOptions(zeromat = FALSE),
        "unknown argument"
    )
    expect_error(
        MrmashExpandCovsOptions(singletons = FALSE),
        "unknown argument"
    )
    # Neither carries a default: pecotmr used to copy the engine's own.
    expect_length(MrmashCanonicalCovsOptions(), 0L)
    expect_length(MrmashExpandCovsOptions(), 0L)
})

test_that("MrmashOptions takes mr.mash's own argument names", {
    skip_if_not_installed("mr.mashr")
    ma <- MrmashOptions(max_iter = 10, tol = 0.5)
    expect_equal(ma$max_iter, 10)
    # The renamed pecotmr spellings are gone; mr.mash's names are the contract.
    expect_error(MrmashOptions(maxIter = 10), "unknown argument")
})

test_that("MrmashOptions and MrmashRssOptions refuse nthreads", {
    skip_if_not_installed("mr.mashr")
    expect_error(MrmashOptions(nthreads = 2L), "`nthreads` cannot be set")
    expect_error(MrmashRssOptions(nthreads = 2L), "`nthreads` cannot be set")
})

test_that("mr.mash and mr.mash.rss fit on one thread", {
    # mr.mash names its thread count only when it runs more than one, so its
    # verbose header pins single-threading on both paths. The header prints
    # only with verbose = TRUE, so it also shows that mrmashWeights() hands
    # mr.mash's own options on to the engine (it once spliced them into
    # mrmashWrapper's arguments, which have none of those names).
    skip_if_not_installed("mr.mashr")
    skip_if_not_installed("glmnet")
    m <- .rrwMulti(n = 60, p = 6, K = 3)
    individual <- capture.output(suppressWarnings(suppressMessages(
        mrmashWeights(
            X = m$X,
            Y = m$Y,
            methodArgs = MrmashOptions(
                canonicalPriorMatrices = TRUE,
                verbose = TRUE,
                max_iter = 2
            )
        )
    )))
    rss <- capture.output(suppressWarnings(mrmashRssWeights(
        m$stat,
        m$ldReference,
        methodArgs = MrmashRssOptions(verbose = TRUE, max_iter = 2)
    )))
    for (out in list(individual, rss)) {
        expect_match(out, "Fitting the optimization algorithm", all = FALSE)
        expect_no_match(out, "RcppParallel threads")
    }
})

test_that("mrmashWrapper refuses to let methodArgs override derived values", {
    skip_if_not_installed("mr.mashr")
    skip_if_not_installed("glmnet")
    set.seed(1)
    X <- matrix(rnorm(60 * 4), 60, 4)
    Y <- matrix(rnorm(60 * 2), 60, 2)
    expect_error(
        mrmashWrapper(
            X,
            Y,
            dataDrivenPriorMatrices = list(U = list(diag(2))),
            methodArgs = MrmashOptions(S0 = list(diag(2)))
        ),
        "supplied by pecotmr"
    )
})

test_that("near-duplicate contexts get a floor on V's eigenvalues", {
    skip_if_not_installed("mr.mashr")
    set.seed(1)
    y <- rnorm(60)
    nearDuplicate <- cbind(y, y + 1e-5 * rnorm(60))
    initV <- function(Y, floorValue = 1e-3) {
        pecotmr:::.mrmashInitV(Y, NULL, TRUE, "full", floorValue)
    }
    eigenRatio <- function(V) {
        values <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
        min(values) / max(values)
    }
    expect_warning(floored <- initV(nearDuplicate), "close to singular")
    expect_false(floored$updateV)
    expect_equal(eigenRatio(floored$V), 1e-3, tolerance = 1e-3)
    # The floor is relative: rescaling Y rescales V, not the decision.
    expect_warning(scaled <- initV(1e3 * nearDuplicate), "close to singular")
    expect_equal(eigenRatio(scaled$V), 1e-3, tolerance = 1e-3)
    # Distinct contexts keep V as it was estimated, and its update.
    distinct <- cbind(y, rnorm(60))
    expect_no_warning(kept <- initV(distinct))
    expect_true(kept$updateV)
    expect_equal(
        kept$V,
        pecotmr:::.mrmashComputeVInit(distinct)
    )
    # A user's floor decides whether they count as near-duplicate.
    expect_warning(initV(distinct, 0.99), "close to singular")
})

test_that("V starts at the pairwise covariance, or at its diagonal", {
    set.seed(3)
    Y <- matrix(rnorm(90), 30L, 3L)
    expect_equal(pecotmr:::.mrmashComputeVInit(Y), cov(Y))
    # Two studies observing one condition each: zero between them.
    blocks <- Y[, 1:2]
    blocks[16:30, 1] <- NA
    blocks[1:15, 2] <- NA
    V <- pecotmr:::.mrmashComputeVInit(blocks)
    expect_equal(V[1, 2], 0)
    expect_equal(diag(V), apply(blocks, 2, var, na.rm = TRUE))
    # Pairs observed on different rows can disagree past any covariance:
    # the first two move together, the second and third too, and the first
    # and third against each other. Then only the variances are kept.
    x <- rnorm(30)
    clash <- cbind(
        c(x[1:10], rep(NA, 10), x[21:30]),
        c(x[1:10], x[11:20], rep(NA, 10)),
        c(rep(NA, 10), x[11:20], -x[21:30])
    )
    pairwise <- cov(clash, use = "pairwise.complete.obs")
    expect_lt(min(eigen(pairwise, symmetric = TRUE)$values), 0)
    expect_equal(
        pecotmr:::.mrmashComputeVInit(clash),
        diag(diag(pairwise))
    )
})

test_that("V can start at the mash prior's residual correlation", {
    set.seed(5)
    Y <- cbind(a = rnorm(40, sd = 2), b = rnorm(40, sd = 3), c = rnorm(40))
    R <- matrix(c(1, 0.5, 0.2, 0.5, 1, 0.3, 0.2, 0.3, 1), 3L, 3L)
    initV <- function(Y, ..., floorValue = 1e-3) {
        pecotmr:::.mrmashInitV(Y, NULL, TRUE, "full", floorValue, ...)
    }
    # The correlation on Y's scale, updated as the covariance of Y would be.
    start <- initV(Y, "mash", R)
    sds <- apply(Y, 2, sd)
    expect_equal(
        start$V,
        `dimnames<-`(R * tcrossprod(sds), list(colnames(Y), colnames(Y)))
    )
    expect_true(start$updateV)
    # The default is the covariance of Y.
    expect_equal(initV(Y)$V, cov(Y))
    expect_equal(initV(Y, "data", R)$V, cov(Y))
    # Two studies never observed together: zero between them, as the
    # covariance of Y is.
    blocks <- Y
    blocks[21:40, 1] <- NA
    blocks[1:20, 2:3] <- NA
    V <- initV(blocks, "mash", R)$V
    expect_equal(unname(V[1L, 2:3]), c(0, 0))
    expect_equal(
        unname(V[2L, 3L]),
        R[2L, 3L] *
            sd(blocks[, 2L], na.rm = TRUE) *
            sd(blocks[, 3L], na.rm = TRUE)
    )
    # A near-singular correlation is floored like a near-singular covariance.
    nearDuplicate <- matrix(c(1, 1 - 1e-9, 1 - 1e-9, 1), 2L, 2L)
    expect_warning(
        floored <- initV(Y[, 1:2], "mash", nearDuplicate),
        "close to singular"
    )
    expect_false(floored$updateV)
})

test_that("a mash start needs the correlation; a supplied V is its own", {
    Y <- matrix(rnorm(60), 20L, 3L)
    expect_error(
        pecotmr:::.mrmashInitV(Y, NULL, TRUE, "full", 1e-3, "mash", NULL),
        "carries no `residualCorrelation`"
    )
    expect_error(
        pecotmr:::.mrmashInitV(Y, NULL, TRUE, "full", 1e-3, "mash", diag(2)),
        "residualCorrelation"
    )
    expect_error(
        pecotmr:::.mrmashInitV(Y, diag(3), TRUE, "full", 1e-3, "data", NULL),
        "cannot also be given"
    )
    expect_identical(
        pecotmr:::.mrmashInitV(Y, diag(3), FALSE, "full", 1e-3)$V,
        diag(3)
    )
    # The diagonal update method keeps only the start's variances.
    R <- 0.5 + diag(0.5, 3)
    expect_equal(
        unname(pecotmr:::.mrmashInitV(Y, NULL, TRUE, "diagonal", 1e-3)$V),
        diag(apply(Y, 2, var))
    )
    expect_equal(
        unname(
            pecotmr:::.mrmashInitV(Y, NULL, TRUE, "diagonal", 1e-3, "mash", R)$V
        ),
        diag(apply(Y, 2, var))
    )
})

test_that("vInitMethod reaches the wrapper and stops short of mr.mash", {
    skip_if_not_installed("mr.mashr")
    skip_if_not_installed("glmnet")
    m <- .rrwMulti(n = 60, p = 6, K = 3)
    R <- 0.3 + diag(0.7, 3)
    # With the update off, the V mr.mash fits with is the start.
    fit <- suppressMessages(mrmashWrapper(
        m$X,
        m$Y,
        dataDrivenPriorMatrices = list(
            U = list(shared = matrix(1, 3L, 3L)),
            residualCorrelation = R
        ),
        prior = MrmashPriorParam(canonicalPriorMatrices = TRUE),
        updateV = FALSE,
        methodArgs = MrmashOptions(vInitMethod = "mash", max_iter = 20)
    ))
    expect_equal(unname(fit$V), R * tcrossprod(apply(m$Y, 2, sd)))
})

test_that("vInitMethod must be data or mash, on the individual path only", {
    skip_if_not_installed("mr.mashr")
    for (bad in list("Mash", c("data", "mash"), NA_character_, 1)) {
        expect_error(
            MrmashOptions(vInitMethod = bad),
            "must be \"data\" or \"mash\""
        )
    }
    expect_identical(
        as.list(MrmashOptions(vInitMethod = "mash"))$vInitMethod,
        "mash"
    )
    expect_error(MrmashRssOptions(vInitMethod = "mash"))
})

test_that("mrmashWrapper fits studies stacked with block-missing Y", {
    skip_if_not_installed("mr.mashr")
    skip_if_not_installed("glmnet")
    skip_if_not(
        is_in("missing_y_method", names(formals(mr.mashr::mr.mash))),
        "mr.mashr predates missing_y_method"
    )
    m <- .rrwMulti(n = 80, p = 6, K = 2)
    # Two studies on different individuals, each observing one condition;
    # a variant the second study lacks is zero on its rows.
    Y <- m$Y
    Y[41:80, 1] <- NA
    Y[1:40, 2] <- NA
    X <- m$X
    X[41:80, 6] <- 0
    fit <- suppressMessages(mrmashWrapper(
        X,
        Y,
        prior = MrmashPriorParam(canonicalPriorMatrices = TRUE),
        methodArgs = MrmashOptions(max_iter = 50)
    ))
    expect_equal(dim(fit$mu1), c(6L, 2L))
    expect_null(fit$Y)
    expect_equal(unname(fit$V[1, 2]), 0)
})

test_that("vEigenFloor reaches the wrapper and stops short of mr.mash", {
    skip_if_not_installed("mr.mashr")
    skip_if_not_installed("glmnet")
    m <- .rrwMulti(n = 60, p = 6, K = 3)
    fitWith <- function(floorValue) {
        mrmashWeights(
            X = m$X,
            Y = m$Y,
            methodArgs = MrmashOptions(
                canonicalPriorMatrices = TRUE,
                vEigenFloor = floorValue
            )
        )
    }
    expect_warning(
        suppressMessages(fitWith(0.99)),
        "vEigenFloor = 0.99"
    )
    expect_no_warning(
        suppressMessages(fitWith(0)),
        message = "close to singular"
    )
})

test_that("vEigenFloor must be one fraction, on the individual path only", {
    skip_if_not_installed("mr.mashr")
    for (bad in list(-0.1, 1, NA_real_, c(1e-6, 1e-5), "small")) {
        expect_error(MrmashOptions(vEigenFloor = bad), "must be one number")
    }
    expect_equal(as.list(MrmashOptions(vEigenFloor = 0))$vEigenFloor, 0)
    expect_error(MrmashRssOptions(vEigenFloor = 1e-6))
})

test_that("a flat prior name is routed into the prior group", {
    # The pipeline spells per-token kwargs flat, so mrmashWeights routes
    # prior-construction names into mrmashWrapper's `prior` record.
    split <- pecotmr:::.mrmashSplitPriorArgs(
        list(canonicalPriorMatrices = TRUE, max_iter = 5)
    )
    expect_equal(names(split$prior), "canonicalPriorMatrices")
    expect_equal(names(split$rest), "max_iter")
})

test_that(".assertRetainX takes TRUE only for a full fit", {
    check <- pecotmr:::.assertRetainX
    expect_silent(check(TRUE, "full"))
    expect_silent(check(FALSE, "full"))
    for (level in c("none", "slim")) {
        expect_silent(check(FALSE, level))
        expect_error(
            check(TRUE, level),
            sprintf("requires `fitRetention = \"full\"`, not \"%s\"", level),
            fixed = TRUE
        )
    }
    expect_error(check(NA, "full"), "retainX")
    expect_error(check("yes", "full"), "retainX")
})

test_that("each weight function holding X rejects retainX unless full", {
    X <- matrix(0, 4L, 2L)
    y <- numeric(4L)
    rejects <- function(fn, ...) {
        expect_error(
            fn(X, y, ..., retainX = TRUE),
            "requires `fitRetention",
            fixed = TRUE
        )
        expect_error(
            fn(X, y, ..., fitRetention = "slim", retainX = TRUE),
            "requires `fitRetention",
            fixed = TRUE
        )
    }
    rejects(glmnetWeights, alpha = 1)
    rejects(lassoWeights)
    rejects(enetWeights)
    rejects(ridgeWeights)
    rejects(pecotmr:::ncvregWeights, penalty = "MCP")
    rejects(scadWeights)
    rejects(mcpWeights)
    rejects(mrashWeights)
    rejects(dprWeights)
    rejects(dprVbWeights)
    rejects(dprGibbsWeights)
    rejects(dprAdaptiveGibbsWeights)
})

# ---------------------------------------------------------------------------
# Fits across studies (difference penalty)
# ---------------------------------------------------------------------------

# Two studies on disjoint halves of eqtlRegionExample's people, stacked as
# .xsStack() stacks them: X over both studies' rows, Y one column per study,
# missing outside its own rows. The variants are kept at |r| < 0.95 to one
# another, as perfect-LD copies leave the lasso's optimum unique only up to
# how it splits an effect between them. The last variant is absent from the
# second study, a column of zeros in its rows.
.rrw_coupledStack <- function(nVariants = 40L) {
    data(eqtlRegionExample, envir = environment())
    X0 <- eqtlRegionExample$X
    p <- colMeans(X0) / 2
    cand <- which(pmin(p, 1 - p) >= 0.05)[301:500]
    R <- abs(stats::cor(X0[, cand]))
    kept <- Reduce(
        function(held, j) if (all(R[j, held] < 0.95)) c(held, j) else held,
        seq_along(cand),
        integer(0)
    )
    cols <- cand[kept][seq_len(nVariants)]
    withr::local_seed(7L)
    Xs <- map(list(1:250, 251:415), function(r) scale(X0[r, cols]))
    ys <- map(Xs, function(X) {
        g <- drop(X[, c(5, 20)] %*% c(0.4, -0.3))
        as.numeric(scale(g + stats::rnorm(nrow(X))))
    })
    X <- `colnames<-`(
        rbind(Xs[[1L]], Xs[[2L]]),
        sprintf("chr1:%d:A:G", 100L * seq_len(nVariants))
    )
    X[251:415, nVariants] <- 0
    list(
        X = X,
        Y = cbind(
            A = c(ys[[1L]], rep(NA, 165L)),
            B = c(rep(NA, 250L), ys[[2L]])
        )
    )
}

# glmnet's fit of one study alone, at each of `lambda`.
.rrw_separateGlmnet <- function(X, y, lambda, alpha) {
    map(lambda, function(l) {
        fit <- glmnet::glmnet(
            X,
            y,
            alpha = alpha,
            lambda = l,
            standardize = FALSE,
            thresh = 1e-14
        )
        as.numeric(stats::coef(fit))[-1L]
    })
}

# Three contexts of one study on the same samples, the second missing on
# some: each context's own phenotype from two shared effects.
.rrw_coupledContexts <- function() {
    s <- .rrw_coupledStack(30L)
    X <- s$X
    withr::local_seed(3L)
    g <- drop(X[, c(5, 20)] %*% c(0.4, -0.3))
    Y <- sapply(1:3, function(i) g + stats::rnorm(nrow(X)))
    Y[sample(nrow(X), 40L), 2L] <- NA
    list(X = X, Y = `colnames<-`(Y, c("c1", "c2", "c3")))
}

# The design of `conds` at the strengths of its cell `i`.
.rrw_coupledDesign <- function(conds, i, layout = "conditions") {
    pecotmr:::.coupledDesign(
        conds,
        conds$rows,
        pecotmr:::.coupledStrength(i, conds),
        layout
    )
}

test_that("a fit across studies at c = 0 is the studies' separate fits", {
    skip_if_not_installed("glmnet")
    s <- .rrw_coupledStack()
    conds <- pecotmr:::.coupledConditions(s$X, s$Y, 0)
    design <- .rrw_coupledDesign(conds, 1L)
    expect_equal(
        map_int(design$blocks, function(b) length(b$carried)),
        c(40L, 39L)
    )
    for (alpha in c(1, 0.5, 0)) {
        engine <- list(
            alpha = alpha,
            grid = c(0.2, 0.05),
            args = list(thresh = 1e-14)
        )
        beta <- suppressWarnings(pecotmr:::.coupledGlmnetPath(design, engine))
        for (k in 1:2) {
            block <- design$blocks[[k]]
            rows <- conds$rows[[k]]
            # glmnet divides only its L2 term by the phenotype's SD, so the
            # separate fit matching lambda is the one at lambda * sigma_k.
            sep <- suppressWarnings(.rrw_separateGlmnet(
                s$X[rows, block$carried],
                s$Y[rows, k],
                engine$grid * block$sigma,
                alpha
            ))
            joint <- block$sigma * beta[block$cols, , drop = FALSE]
            expect_equal(
                unname(joint),
                exec(cbind, !!!sep),
                tolerance = 1e-5,
                info = paste("alpha", alpha, "study", k)
            )
        }
    }
})

test_that("a strong coupling gives the studies one pooled fit", {
    skip_if_not_installed("glmnet")
    s <- .rrw_coupledStack()
    keep <- 1:39 # the variants both studies carry
    stack <- list(X = s$X[, keep], Y = s$Y)
    conds <- pecotmr:::.coupledConditions(stack$X, stack$Y, 1000)
    design <- .rrw_coupledDesign(conds, 1L)
    # glmnet's own tolerance: the coupled pairs' steps shrink as c grows, so a
    # far tighter one takes more passes than glmnet allows.
    engine <- list(alpha = 1, grid = 0.01, args = list())
    beta <- pecotmr:::.coupledGlmnetPath(design, engine)
    # The penalty couples each variant's per-SD effect, its coefficient
    # times its SD in the study's own rows.
    perSd <- map2(design$blocks, conds$rows, function(b, rows) {
        Xk <- stack$X[rows, b$carried]
        sqrt(colMeans(sweep(Xk, 2L, colMeans(Xk))^2)) * beta[b$cols, 1L]
    })
    expect_lt(max(abs(perSd[[1L]] - perSd[[2L]])), 1e-3)
    expect_gt(sum(perSd[[1L]] != 0), 0L)
})

test_that("each study takes its own grid point, its absent variants at 0", {
    skip_if_not_installed("glmnet")
    s <- .rrw_coupledStack()
    withr::local_seed(1L)
    W <- suppressWarnings(pecotmr:::.lassoCoupledWeights(
        s$X,
        s$Y,
        c(0.1, 1, 10),
        fitRetention = "slim"
    ))
    expect_equal(dimnames(W), list(colnames(s$X), c("A", "B")))
    expect_equal(W[40L, "B"], 0)
    fit <- attr(W, "fit")
    expect_named(fit, c("strength", "lambda", "cvm", "cvsd", "variantIds"))
    expect_equal(dimnames(fit$strength), list(c("A", "B"), "strength"))
    expect_true(all(fit$strength %in% c(0.1, 1, 10)))
    expect_equal(fit$variantIds, colnames(s$X))
    none <- suppressWarnings(pecotmr:::.lassoCoupledWeights(s$X, s$Y, 1))
    expect_null(attr(none, "fit"))
    full <- attr(
        suppressWarnings(pecotmr:::.ridgeCoupledWeights(
            s$X,
            s$Y,
            c(0.1, 1),
            fitRetention = "full"
        )),
        "fit"
    )
    expect_length(full$paths, 2L)
    expect_equal(nrow(full$cells), 2L)
    expect_length(full$lambdaPath, 100L)
    expect_equal(dim(full$cvmPath$A), c(100L, 2L))
})

test_that("L0Learn fits across studies with the L0 penalty only", {
    skip_if_not_installed("L0Learn")
    s <- .rrw_coupledStack()
    withr::local_seed(1L)
    W <- pecotmr:::.l0learnCoupledWeights(s$X, s$Y, c(0.1, 1))
    expect_equal(dim(W), c(40L, 2L))
    expect_equal(W[40L, "B"], 0)
    expect_error(
        pecotmr:::.l0learnCoupledWeights(s$X, s$Y, 1, penalty = "L0L2"),
        "fits penalty = \"L0\" only, not \"L0L2\""
    )
    expect_error(
        pecotmr:::.l0learnCoupledWeights(
            s$X,
            s$Y,
            1,
            methodArgs = L0learnOptions(excludeFirstK = 1L)
        ),
        "L0learnOptions\\(excludeFirstK\\) does not apply"
    )
})

test_that("a coupled fit refuses what it cannot honour", {
    skip_if_not_installed("glmnet")
    s <- .rrw_coupledStack()
    expect_error(
        pecotmr:::.lassoCoupledWeights(
            s$X,
            s$Y,
            1,
            methodArgs = GlmnetOptions(penalty.factor = rep(1, 40L))
        ),
        "GlmnetOptions\\(penalty.factor\\) does not apply"
    )
    expect_error(
        pecotmr:::.lassoCoupledWeights(
            s$X,
            s$Y,
            1,
            methodArgs = GlmnetOptions(type.measure = "mae")
        ),
        "chooses its grid point by mean squared error"
    )
    expect_error(
        pecotmr:::.lassoCoupledWeights(s$X, s$Y, -1),
        "strengths"
    )
    expect_error(
        pecotmr:::.mcpCoupledWeights(
            s$X,
            s$Y,
            1,
            methodArgs = NcvregOptions(alpha = 0.5)
        ),
        "NcvregOptions\\(alpha\\) does not apply"
    )
    expect_error(
        pecotmr:::.lassoCoupledWeights(s$X, s$Y, 1, backend = "ncvreg"),
        "must be one of glmnet, grpreg"
    )
})

test_that("a coupled fit holds a sample out of every condition at once", {
    skip_if_not_installed("glmnet")
    d <- .rrw_coupledContexts()
    folds <- pecotmr:::.groupedFolds(d$Y, 5L)
    # A sample's rows in every context share its fold, so no context's
    # held-out samples train another's fit.
    expect_length(folds, nrow(d$X))
    expect_setequal(folds, 1:5)
    withr::local_seed(1L)
    W <- pecotmr:::.lassoCoupledWeights(
        d$X,
        d$Y,
        c(0.1, 1),
        fitRetention = "slim"
    )
    expect_equal(dim(W), c(30L, 3L))
    expect_equal(rownames(attr(W, "fit")$strength), c("c1", "c2", "c3"))
})

test_that("ncvreg's coupled MCP is its separate fits at c = 0", {
    skip_if_not_installed("ncvreg")
    d <- .rrw_coupledContexts()
    conds <- pecotmr:::.coupledConditions(d$X, d$Y[, c(1L, 3L)], 0)
    design <- .rrw_coupledDesign(conds, 1L)
    engine <- list(
        penalty = "MCP",
        gamma = 3,
        grid = c(0.2, 0.08),
        args = list(eps = 1e-10)
    )
    beta <- pecotmr:::.coupledNcvregPath(design, engine)
    for (k in 1:2) {
        y <- d$Y[, c(1L, 3L)[[k]]]
        scaled <- (y - mean(y)) / sqrt(mean((y - mean(y))^2))
        own <- ncvreg::ncvreg(
            d$X,
            scaled,
            penalty = "MCP",
            gamma = 3,
            lambda = engine$grid,
            eps = 1e-10
        )
        expect_equal(
            unname(beta[design$blocks[[k]]$cols, ]),
            unname(own$beta[-1L, ]),
            tolerance = 1e-8,
            info = k
        )
    }
})

test_that("ncvreg's coupled MCP solves the coupled objective", {
    skip_if_not_installed("ncvreg")
    d <- .rrw_coupledContexts()
    Y <- d$Y[, c(1L, 3L)]
    strength <- 0.7
    lambda <- 0.08
    conds <- pecotmr:::.coupledConditions(d$X, Y, strength)
    design <- .rrw_coupledDesign(conds, 1L)
    engine <- list(
        penalty = "MCP",
        gamma = 3,
        grid = c(0.3, lambda),
        args = list(eps = 1e-12, max.iter = 100000L)
    )
    beta <- pecotmr:::.coupledNcvregPath(design, engine)[, 2L]
    copies <- pecotmr:::.coupledPlaced(
        pecotmr:::.coupledCopies(conds, conds$rows, "conditions"),
        "conditions"
    )
    perSd <- map(copies, function(cp) beta[cp$cols] * cp$s)
    # sum_k RSS_k / (2 n) + MCP(s b; lambda, 3) + (c / 2) (s b - s' b')^2
    for (k in 1:2) {
        cp <- copies[[k]]
        loss <- -drop(crossprod(cp$X, cp$Y[, 1L] - cp$X %*% beta[cp$cols]))
        coupling <- strength * cp$s * (perSd[[k]] - perSd[[3L - k]])
        t <- perSd[[k]]
        mcp <- cp$s *
            ifelse(abs(t) <= 3 * lambda, sign(t) * (lambda - abs(t) / 3), 0)
        on <- beta[cp$cols] != 0
        expect_gt(sum(on), 0L)
        expect_lt(max(abs((loss / cp$n + coupling + mcp)[on])), 1e-8)
        expect_true(all(
            abs((loss / cp$n + coupling)[!on]) <= lambda * cp$s[!on]
        ))
    }
})

test_that("glmnet couples studies and groups each one's contexts at once", {
    skip_if_not_installed("glmnet")
    d <- .rrw_coupledContexts()
    rows <- list(A = 1:250, B = 251:415)
    Y <- cbind(
        A_c1 = replace(d$Y[, 1L], rows$B, NA),
        A_c3 = replace(d$Y[, 3L], rows$B, NA),
        B_c1 = replace(d$Y[, 1L], rows$A, NA),
        B_c3 = replace(d$Y[, 3L], rows$A, NA)
    )
    cond <- data.frame(
        study = c("A", "A", "B", "B"),
        context = c("c1", "c3", "c1", "c3"),
        trait = "t"
    )
    spec <- pecotmr:::.jointCouplingSpec(cond, list(study = 2))
    expect_equal(spec$unit, c(1L, 1L, 2L, 2L))
    expect_equal(spec$column, c(1L, 2L, 1L, 2L))
    conds <- pecotmr:::.coupledConditions(d$X, Y, spec)
    expect_true(pecotmr:::.coupledUnitsHold(conds))
    design <- .rrw_coupledDesign(conds, 1L, "units")
    lambda <- 0.05
    beta <- suppressWarnings(pecotmr:::.coupledGlmnetPath(
        design,
        list(
            alpha = 1,
            grid = c(0.5, lambda),
            args = list(thresh = 1e-16, maxit = 1e7)
        )
    ))[, 2L]
    copies <- pecotmr:::.coupledPlaced(
        pecotmr:::.coupledCopies(conds, conds$rows, "units"),
        "units"
    )
    width <- ncol(design$X)
    B <- map(copies, function(cp) {
        cbind(beta[cp$cols], beta[width + cp$cols])
    })
    # sum RSS / (2 n) + lambda sum_j ||B_j|| + (c / 2) sum (s b - s' b')^2,
    # a variant's contexts one group in each study.
    for (u in 1:2) {
        cp <- copies[[u]]
        other <- copies[[3L - u]]
        at <- match(cp$carried, other$carried)
        loss <- -crossprod(cp$X, cp$Y - cp$X %*% B[[u]]) / cp$n
        # A variant the other study lacks is not coupled.
        pulled <- 2 * cp$s * (cp$s * B[[u]] - other$s[at] * B[[3L - u]][at, ])
        coupling <- replace(pulled, is.na(pulled), 0)
        norms <- sqrt(rowSums(B[[u]]^2))
        on <- norms > 0
        expect_gt(sum(on), 0L)
        kkt <- loss + coupling + lambda * B[[u]] / ifelse(on, norms, 1)
        expect_lt(max(abs(kkt[on, ])), 1e-6)
    }
    # Coupled along the contexts instead, each condition is its own unit.
    byContext <- pecotmr:::.jointCouplingSpec(cond, list(context = 2))
    expect_false(pecotmr:::.coupledUnitsHold(
        pecotmr:::.coupledConditions(d$X, Y, byContext)
    ))
})

test_that("a pair is coupled along the one coupled axis it differs along", {
    cond <- data.frame(
        study = c("A", "A", "B", "B"),
        context = c("c1", "c3", "c1", "c3"),
        trait = "t"
    )
    spec <- pecotmr:::.jointCouplingSpec(
        cond,
        list(context = c(0.1, 1), study = 5)
    )
    expect_equal(
        spec$pairs,
        matrix(
            c(
                NA,
                "context",
                "study",
                NA,
                "context",
                NA,
                NA,
                "study",
                "study",
                NA,
                NA,
                "context",
                NA,
                "study",
                "context",
                NA
            ),
            4L,
            byrow = TRUE
        )
    )
    expect_named(spec$grids, c("study", "context"))
    conds <- pecotmr:::.coupledConditions(
        matrix(0, 2L, 1L),
        matrix(0, 2L, 4L),
        spec
    )
    expect_equal(nrow(conds$cells), 2L)
    W <- pecotmr:::.coupledStrength(2L, conds)
    expect_equal(W[1L, 2L], 1)
    expect_equal(W[1L, 3L], 5)
    expect_equal(W[1L, 4L], 0)
})

test_that("grpreg and abess couple the conditions with their group penalty", {
    skip_if_not_installed("grpreg")
    skip_if_not_installed("abess")
    d <- .rrw_coupledContexts()
    rlang::local_options(rlib_warning_verbosity = "verbose")
    withr::local_seed(1L)
    expect_warning(
        W <- pecotmr:::.mcpCoupledWeights(
            d$X,
            d$Y,
            c(0.1, 1),
            backend = "grpreg"
        ),
        "the penalty is approximate on this group"
    )
    expect_equal(dim(W), c(30L, 3L))
    # A variant is in every condition's weights or none, at the chosen grid
    # points' shared support where the conditions pick the same one.
    conds <- pecotmr:::.coupledConditions(d$X, d$Y, 1)
    design <- .rrw_coupledDesign(conds, 1L, "groups")
    expect_equal(sort(unique(design$group)), seq_len(30L))
    fit <- pecotmr:::.coupledAbessPath(
        design,
        list(grid = 0:4, args = list())
    )
    # abess's best subset at each size is least squares on its groups, the
    # coupling's rows among them.
    picked <- which(fit[, 4L] != 0)
    ls <- qr.coef(qr(as.matrix(design$X[, picked])), design$y)
    expect_equal(unname(fit[picked, 4L]), unname(ls), tolerance = 1e-6)
    expect_equal(length(picked), 3L * length(unique(design$group[picked])))
    withr::local_seed(1L)
    W0 <- pecotmr:::.l0learnCoupledWeights(
        d$X,
        d$Y,
        c(0.1, 1),
        backend = "abess",
        fitRetention = "slim"
    )
    expect_named(
        attr(W0, "fit"),
        c("strength", "supportSize", "cvm", "cvsd", "variantIds")
    )
})

test_that("ridge fits a group's conditions at one lambda", {
    skip_if_not_installed("glmnet")
    d <- .rrw_coupledContexts()
    withr::local_seed(1L)
    W <- pecotmr:::.ridgeGroupedWeights(d$X, d$Y, fitRetention = "slim")
    fit <- attr(W, "fit")
    expect_null(fit$strength)
    expect_length(unique(fit$lambda), 1L)
    expect_true(all(W != 0))
})

test_that("ncvreg's MCP warns when its columns' couplings differ", {
    skip_if_not_installed("ncvreg")
    s <- .rrw_coupledStack(20L)
    rlang::local_options(rlib_warning_verbosity = "verbose")
    withr::local_seed(1L)
    expect_warning(
        pecotmr:::.mcpCoupledWeights(s$X[, 1:19], s$Y, c(0.1, 1)),
        NA
    )
    # The last variant is absent from study B, so its column in A has no
    # partner while every other has one.
    expect_warning(
        pecotmr:::.mcpCoupledWeights(s$X, s$Y, c(0.1, 1)),
        "mcp on ncvreg across coupled conditions"
    )
    expect_warning(
        pecotmr:::.scadCoupledWeights(s$X[, 1:19], s$Y, c(0.1, 1)),
        "scad on ncvreg across coupled conditions"
    )
})

# ---------------------------------------------------------------------------
# Penalized RSS across studies (difference penalty)
# ---------------------------------------------------------------------------

# The two studies of .rrw_coupledStack() as summary statistics: each one's
# marginal correlations and in-sample LD over the variants it carries, on
# population-scaled genotypes and phenotype, so that the RSS objective is
# exactly twice the individual-level one.
.rrw_coupledRss <- function() {
    s <- .rrw_coupledStack()
    popScale <- function(v) {
        centred <- v - mean(v)
        centred / sqrt(mean(centred^2))
    }
    map(1:2, function(k) {
        rows <- which(!is.na(s$Y[, k]))
        X <- apply(s$X[rows, ], 2L, popScale)
        X <- X[, apply(X, 2L, function(v) all(is.finite(v))), drop = FALSE]
        y <- popScale(s$Y[rows, k])
        list(
            X = X,
            y = y,
            z = as.numeric(crossprod(X, y)) / nrow(X),
            R = crossprod(X) / nrow(X)
        )
    })
}

# The strengths of two conditions coupled at `strength`.
.rrwPairs <- function(strength) {
    strength * (1 - diag(2L))
}

.rrw_coupledRssPath <- function(fit, strength, lambda, penalty = "lasso") {
    penalizedRssCoupledCpp(
        map(fit, "z"),
        map(fit, "R"),
        map(fit, function(f) match(colnames(f$X), colnames(fit[[1L]]$X))),
        ncol(fit[[1L]]$X),
        .rrwPairs(strength),
        FALSE,
        0,
        lambda,
        penalty,
        "strongRule",
        if (penalty == "SCAD") 3.7 else 3,
        1,
        0,
        0,
        1e-10,
        100000L,
        0L,
        pecotmr:::.rssPathControl()
    )
}

test_that("the coupled RSS path at c = 0 is each study's own", {
    fit <- .rrw_coupledRss()
    lambda <- c(0.2, 0.1, 0.05)
    for (penalty in c("lasso", "MCP", "SCAD")) {
        joint <- .rrw_coupledRssPath(fit, 0, lambda, penalty)
        for (k in 1:2) {
            alone <- penalizedRssCpp(
                fit[[k]]$z,
                fit[[k]]$R,
                0,
                lambda,
                penalty,
                "strongRule",
                if (penalty == "SCAD") 3.7 else 3,
                1,
                0,
                0,
                1e-10,
                100000L,
                0L,
                pecotmr:::.rssPathControl()
            )
            expect_equal(
                joint$studies[[k]]$beta,
                alone$beta,
                tolerance = 1e-8,
                info = paste(penalty, k)
            )
        }
    }
})

test_that("coupled RSS on in-sample LD is the coupled individual-level fit", {
    skip_if_not_installed("glmnet")
    # On standardized data the RSS objective is twice the individual-level
    # one, so the same (lambda, c) gives the same coefficients.
    fit <- .rrw_coupledRss()
    s <- .rrw_coupledStack()
    keep <- colnames(fit[[2L]]$X)
    stack <- list(
        X = rbind(fit[[1L]]$X[, keep], fit[[2L]]$X),
        Y = cbind(
            A = c(fit[[1L]]$y, rep(NA, nrow(fit[[2L]]$X))),
            B = c(rep(NA, nrow(fit[[1L]]$X)), fit[[2L]]$y)
        )
    )
    shared <- map(fit, function(f) {
        list(X = f$X[, keep], z = f$z[match(keep, colnames(f$X))])
    })
    rss <- map(shared, function(f) {
        list(X = f$X, z = f$z, R = crossprod(f$X) / nrow(f$X))
    })
    lambda <- c(0.2, 0.1, 0.05)
    for (strength in c(0.1, 1)) {
        joint <- .rrw_coupledRssPath(rss, strength, lambda)
        conds <- pecotmr:::.coupledConditions(stack$X, stack$Y, strength)
        design <- .rrw_coupledDesign(conds, 1L)
        beta <- suppressWarnings(pecotmr:::.coupledGlmnetPath(
            design,
            list(alpha = 1, grid = lambda, args = list(thresh = 1e-14))
        ))
        for (k in 1:2) {
            expect_equal(
                joint$studies[[k]]$beta,
                beta[design$blocks[[k]]$cols, , drop = FALSE],
                tolerance = 1e-5,
                check.attributes = FALSE,
                info = paste(strength, k)
            )
        }
    }
})

test_that("held-out splits of conditions sharing samples correlate", {
    withr::local_seed(1L)
    unit <- map(1:3, function(k) matrix(stats::rnorm(5L * 4000L), 5L))
    residual <- matrix(c(1, 0.6, 0.3, 0.6, 1, 0.2, 0.3, 0.2, 1), 3L)
    # Conditions 1 and 2 on one reference and variants, 3 on its own.
    run <- list(residual = residual, sources = c(1L, 1L, 3L))
    mixed <- pecotmr:::.rssCoupledCorrelate(unit, run)
    expect_equal(
        stats::cor(c(mixed[[1L]]), c(mixed[[2L]])),
        0.6,
        tolerance = 0.03
    )
    expect_equal(stats::var(c(mixed[[2L]])), 1, tolerance = 0.03)
    expect_identical(mixed[[3L]], unit[[3L]])
    # Without a correlation every condition's splits stay apart.
    expect_identical(
        pecotmr:::.rssCoupledCorrelate(unit, list(sources = c(1L, 1L, 3L))),
        unit
    )
    refs <- list("A", "A", "B")
    stats <- list(
        list(variantNames = c("v1", "v2")),
        list(variantNames = c("v1", "v2")),
        list(variantNames = c("v1", "v2"))
    )
    expect_equal(pecotmr:::.rssSharedSources(refs, stats), c(1L, 1L, 3L))
    stats[[2L]]$variantNames <- "v1"
    expect_equal(pecotmr:::.rssSharedSources(refs, stats), 1:3)
})

test_that("the group penalty across studies solves its objective", {
    fit <- .rrw_coupledRss()
    union <- map(fit, function(f) match(colnames(f$X), colnames(fit[[1L]]$X)))
    nUnion <- ncol(fit[[1L]]$X)
    # sum_k [b_k'R_k b_k - 2 b_k'z_k] + 2 P(||b_j||; lambda sqrt(K_j))
    #     + c sum_j (b_j1 - b_j2)^2
    check <- function(penalty, strength, lambda, gamma = 3) {
        path <- penalizedRssCoupledCpp(
            map(fit, "z"),
            map(fit, "R"),
            union,
            nUnion,
            .rrwPairs(strength),
            TRUE,
            0,
            c(0.5, lambda),
            penalty,
            "strongRule",
            gamma,
            1,
            0,
            0,
            1e-12,
            100000L,
            0L,
            pecotmr:::.rssPathControl()
        )
        B <- matrix(NA_real_, nUnion, 2L)
        G <- matrix(NA_real_, nUnion, 2L)
        for (k in 1:2) {
            b <- path$studies[[k]]$beta[, 2L]
            B[union[[k]], k] <- b
            G[union[[k]], k] <- fit[[k]]$z - drop(fit[[k]]$R %*% b)
        }
        pull <- strength * cbind(B[, 1L] - B[, 2L], B[, 2L] - B[, 1L])
        g <- G - replace(pull, is.na(pull), 0)
        for (j in seq_len(nUnion)) {
            m <- which(!is.na(B[j, ]))
            w <- sqrt(length(m))
            r <- sqrt(sum(B[j, m]^2))
            if (r == 0) {
                expect_lte(sqrt(sum(g[j, m]^2)), lambda * w + 1e-10)
                next
            }
            d <- if (penalty == "lasso") {
                lambda * w
            } else {
                max(lambda * w - r / gamma, 0)
            }
            expect_lt(max(abs(g[j, m] - d * B[j, m] / r)), 1e-8)
        }
        # A variant both studies carry is in both or in neither.
        both <- rowSums(is.na(B)) == 0
        expect_identical(B[both, 1L] != 0, B[both, 2L] != 0)
    }
    for (penalty in c("lasso", "MCP")) {
        for (strength in c(0, 1)) {
            check(penalty, strength, 0.05)
        }
    }
})

test_that("the group L0 penalty swaps whole groups", {
    fit <- .rrw_coupledRss()
    union <- map(fit, function(f) match(colnames(f$X), colnames(fit[[1L]]$X)))
    nUnion <- ncol(fit[[1L]]$X)
    run <- function(maxSwaps, lambda0) {
        penalizedRssCoupledCpp(
            map(fit, "z"),
            map(fit, "R"),
            union,
            nUnion,
            .rrwPairs(0),
            TRUE,
            0,
            lambda0,
            "L0",
            "gradientRanked",
            0,
            1,
            0,
            0,
            1e-10,
            100000L,
            maxSwaps,
            pecotmr:::.rssPathControl()
        )
    }
    # sum_k [b_k'R_k b_k - 2 b_k'z_k] + 2 lambda0 (groups in the model)
    objective <- function(path, at) {
        betas <- map(path$studies, function(s) s$beta[, at])
        loss <- sum(map2_dbl(fit, betas, function(f, b) {
            drop(t(b) %*% f$R %*% b) - 2 * sum(f$z * b)
        }))
        B <- matrix(0, nUnion, 2L)
        for (k in 1:2) {
            B[union[[k]], k] <- betas[[k]]
        }
        loss + 2 * path$lambda[[at]] * sum(rowSums(B != 0) > 0)
    }
    # From zero at one lambda0, a swap is kept only where it lowers the
    # objective, so the swaps never end higher.
    for (lambda0 in c(0.004, 0.002, 0.001)) {
        expect_lte(
            objective(run(5L, lambda0), 1L),
            objective(run(0L, lambda0), 1L) + 1e-10
        )
    }
})

test_that("the coupled RSS path reads a factor of each study's LD as it", {
    fit <- .rrw_coupledRss()
    union <- map(fit, function(f) match(colnames(f$X), colnames(fit[[1L]]$X)))
    rest <- list(
        unionIdxList = union,
        nUnion = ncol(fit[[1L]]$X),
        strengthR = .rrwPairs(1),
        grouped = FALSE,
        shrink = 0,
        pathR = c(0.2, 0.1, 0.05),
        penaltyStr = "lasso",
        strategyStr = "strongRule",
        gamma = 3,
        alpha = 1,
        lambda1 = 0,
        lambda2 = 0,
        thr = 1e-10,
        maxiter = 100000L,
        maxSwaps = 0L,
        control = pecotmr:::.rssPathControl()
    )
    dense <- exec(
        penalizedRssCoupledCpp,
        zList = map(fit, "z"),
        ldList = map(fit, "R"),
        !!!rest
    )
    # R = F'F / (m - 1) for F = X sqrt((m - 1) / m), X population-scaled.
    factors <- map(fit, function(f) f$X * sqrt((nrow(f$X) - 1) / nrow(f$X)))
    viaFactor <- exec(
        penalizedRssCoupledFactorCpp,
        zList = map(fit, "z"),
        factorList = factors,
        !!!rest
    )
    for (k in 1:2) {
        expect_equal(
            viaFactor$studies[[k]]$beta,
            dense$studies[[k]]$beta,
            tolerance = 1e-8,
            info = paste("study", k)
        )
    }
})

test_that("the coupled L0 path swaps as each study's own does at c = 0", {
    fit <- .rrw_coupledRss()
    union <- map(fit, function(f) match(colnames(f$X), colnames(fit[[1L]]$X)))
    path <- c(0.02, 0.01, 0.005)
    joint <- penalizedRssCoupledCpp(
        map(fit, "z"),
        map(fit, "R"),
        union,
        ncol(fit[[1L]]$X),
        .rrwPairs(0),
        FALSE,
        0,
        path,
        "L0",
        "gradientRanked",
        0,
        1,
        0,
        0,
        1e-10,
        100000L,
        3L,
        pecotmr:::.rssPathControl()
    )
    for (k in 1:2) {
        alone <- penalizedRssCpp(
            fit[[k]]$z,
            fit[[k]]$R,
            0,
            path,
            "L0",
            "gradientRanked",
            0,
            1,
            0,
            0,
            1e-10,
            100000L,
            3L,
            pecotmr:::.rssPathControl()
        )
        expect_equal(
            joint$studies[[k]]$beta,
            alone$beta,
            tolerance = 1e-8,
            info = paste("study", k)
        )
    }
    # No swap qualifies on that data. On real LD both studies swap along
    # study 1's path, each lowering its objective, and each round re-solves
    # every study as its own swaps would. Nothing may rest on a tie, which
    # rounding or the joint solve's extra sweeps can break either way: the
    # variants keep |r| < 0.98 (a swap test between duplicates compares equal
    # numbers), and the path skips study 1's entry value, which its own fit
    # records as zero but the joint fit, at study 2's larger pathMax, solves.
    studies <- .rrwTwoStudies(401L, maxR = 0.98)
    ownFit <- function(k, lambda, maxSwaps) {
        penalizedRssCpp(
            studies[[k]]$z,
            studies[[k]]$R,
            0,
            lambda,
            "L0",
            "gradientRanked",
            0,
            1,
            0,
            0,
            1e-10,
            100000L,
            maxSwaps,
            pecotmr:::.rssPathControl()
        )
    }
    sharedPath <- ownFit(1L, numeric(0), 0L)$lambda[2:7]
    coupledFit <- function(maxSwaps) {
        penalizedRssCoupledCpp(
            map(studies, "z"),
            map(studies, "R"),
            list(1:30, 1:30),
            30L,
            .rrwPairs(0),
            FALSE,
            0,
            sharedPath,
            "L0",
            "gradientRanked",
            0,
            1,
            0,
            0,
            1e-10,
            100000L,
            maxSwaps,
            pecotmr:::.rssPathControl()
        )
    }
    swapped <- coupledFit(3L)
    unswapped <- coupledFit(0L)
    for (k in 1:2) {
        gain <- unswapped$studies[[k]]$fbeta - swapped$studies[[k]]$fbeta
        expect_gt(max(gain), 0.01, label = paste("study", k, "swap gain"))
        expect_equal(
            swapped$studies[[k]]$beta,
            ownFit(k, sharedPath, 3L)$beta,
            tolerance = 1e-8,
            info = paste("study", k)
        )
    }
})

test_that("the coupled path flags a point it ran out of passes on", {
    fit <- .rrw_coupledRss()
    union <- map(fit, function(f) match(colnames(f$X), colnames(fit[[1L]]$X)))
    for (strategy in c("strongRule", "gradientRanked")) {
        joint <- penalizedRssCoupledCpp(
            map(fit, "z"),
            map(fit, "R"),
            union,
            ncol(fit[[1L]]$X),
            .rrwPairs(1),
            FALSE,
            0,
            c(0.05, 0.01),
            "lasso",
            strategy,
            3,
            1,
            0,
            0,
            1e-14,
            1L,
            0L,
            pecotmr:::.rssPathControl()
        )
        expect_true(any(joint$studies[[1L]]$conv == 0L), info = strategy)
    }
})

test_that("the coupled path refuses LD that does not fit the studies", {
    studies <- .rrwTwoStudies()
    rest <- list(
        zList = map(studies, "z"),
        unionIdxList = list(1:30, 1:30),
        nUnion = 30L,
        strengthR = .rrwPairs(1),
        grouped = FALSE,
        shrink = 0,
        pathR = c(0.05, 0.01),
        penaltyStr = "lasso",
        strategyStr = "strongRule",
        gamma = 3,
        alpha = 1,
        lambda1 = 0,
        lambda2 = 0,
        thr = 1e-10,
        maxiter = 1000L,
        maxSwaps = 0L,
        control = pecotmr:::.rssPathControl()
    )
    expect_error(
        exec(penalizedRssCoupledCpp, ldList = list(studies[[1L]]$R), !!!rest),
        "There are 1 LD sources for 2 studies"
    )
    expect_error(
        exec(
            penalizedRssCoupledCpp,
            ldList = list(studies[[1L]]$R, studies[[2L]]$R[1:29, 1:29]),
            !!!rest
        ),
        "Study 2's LD is 29 x 29 but it has 30 z-scores"
    )
    expect_error(
        exec(
            penalizedRssCoupledFactorCpp,
            factorList = list(studies[[1L]]$X, studies[[2L]]$X[, 1:29]),
            !!!rest
        ),
        "Study 2's LD is 165 x 29 but it has 30 z-scores"
    )
})

test_that("the coupled path's helpers pad, skip and pair as asked", {
    B <- matrix(1:4, 2L)
    # An engine that stopped early stands for the grid's rest.
    expect_equal(pecotmr:::.coupledPadPath(B, 4L), cbind(B, B[, c(2L, 2L)]))
    expect_equal(pecotmr:::.coupledPadPath(B, 1L), B[, 1L, drop = FALSE])
    # A study with no held-out rows in a fold adds nothing there.
    expect_equal(
        pecotmr:::.coupledHeldOutSse(
            1L,
            integer(0),
            NULL,
            matrix(0, 3, 5),
            NULL
        ),
        numeric(5)
    )
    # Conditions sharing no variant, or not coupled, have no pseudo-rows.
    copies <- list(list(carried = c(1L, 2L)), list(carried = c(3L, 4L)))
    W <- matrix(c(0, 1, 1, 0), 2L)
    expect_null(pecotmr:::.coupledPair(c(1L, 2L), copies, W))
    shared <- list(list(carried = 1:2), list(carried = 1:2))
    expect_null(pecotmr:::.coupledPair(c(1L, 2L), shared, 0 * W))
})

test_that("fits across studies read genotype panels through their factors", {
    d <- .rrwPanelStatLd()
    ids <- d$stat$variantNames
    stats <- list(d$stat, list_modify(d$stat, z = 0.8 * d$stat$z))
    refs <- list(d$panel, d$panel)
    withr::local_seed(1L)
    W <- pecotmr:::.rssCoupledWeights(
        stats,
        refs,
        ids,
        pecotmr:::.rssPenalty("lasso"),
        c(0.5, 0.9),
        c(0.1, 1),
        list(),
        "full",
        "test"
    )
    expect_equal(dim(W), c(length(ids), 2L))
    expect_false(is.null(attr(W, "fit")$beta))
    prs <- pecotmr:::.prsCsCoupledWeights(
        stats,
        refs,
        ids,
        list(nIter = 30L, nBurnin = 10L, thin = 2L, seed = 2L),
        "none",
        "test"
    )
    expect_equal(dim(prs), c(length(ids), 2L))
    expect_null(attr(prs, "fit"))
    sdpr <- function(rho) {
        pecotmr:::.sdprCoupledWeights(
            stats,
            refs,
            ids,
            list(
                iter = 30L,
                burn = 10L,
                thin = 2L,
                M = 10L,
                verbose = FALSE,
                seed = 3L,
                rho = rho
            ),
            "slim",
            "test"
        )
    }
    learned <- attr(sdpr(NULL), "fit")$rhoEst[1L, 2L]
    expect_true(learned >= 0 && learned <= 0.99)
    given <- matrix(c(1, 0.5, 0.5, 1), 2L)
    expect_equal(attr(sdpr(given), "fit")$rhoEst, given)
})

test_that("the coupled RSS solver checks each study's inputs", {
    fit <- .rrw_coupledRss()
    expect_error(
        penalizedRssCoupledCpp(
            map(fit, "z"),
            map(fit, "R"),
            list(1:3, 1:3),
            40L,
            .rrwPairs(1),
            FALSE,
            0,
            0.1,
            "lasso",
            "strongRule",
            3,
            1,
            0,
            0,
            1e-8,
            1000L,
            0L,
            pecotmr:::.rssPathControl()
        ),
        "union indices but"
    )
    expect_error(
        penalizedRssCoupledCpp(
            map(fit, "z"),
            map(fit, "R"),
            map(fit, function(f) rep(1L, ncol(f$X))),
            40L,
            .rrwPairs(1),
            FALSE,
            0,
            0.1,
            "lasso",
            "strongRule",
            3,
            1,
            0,
            0,
            1e-8,
            1000L,
            0L,
            pecotmr:::.rssPathControl()
        ),
        "must be distinct"
    )
})

test_that("each study picks its own (s, c, lambda) from summary statistics", {
    fit <- .rrw_coupledRss()
    ids <- colnames(fit[[1L]]$X)
    stats <- map(fit, function(f) {
        list(
            z = f$z * sqrt(nrow(f$X)),
            n = rep(nrow(f$X), ncol(f$X)),
            variantNames = colnames(f$X)
        )
    })
    refs <- map(fit, function(f) LdMatrix(f$R, nSamples = nrow(f$X)))
    withr::local_seed(1L)
    W <- pecotmr:::.rssCoupledWeights(
        stats,
        refs,
        ids,
        pecotmr:::.rssPenalty("lasso"),
        c(0.2, 0.9),
        c(0.1, 1),
        list(),
        "slim",
        "test"
    )
    expect_equal(dim(W), c(length(ids), 2L))
    expect_equal(rownames(W), ids)
    expect_equal(unname(W[setdiff(ids, colnames(fit[[2L]]$X)), 2L]), 0)
    chosen <- attr(W, "fit")$selection
    expect_equal(nrow(chosen), 2L)
    expect_true(all(chosen$strength %in% c(0.1, 1)))
    expect_true(all(chosen$s %in% c(0.2, 0.9)))
    l0 <- pecotmr:::.rssCoupledWeights(
        stats,
        refs,
        ids,
        pecotmr:::.rssPenalty("L0"),
        0.9,
        1,
        list(),
        "none",
        "test"
    )
    expect_null(attr(l0, "fit"))
    expect_equal(dim(l0), c(length(ids), 2L))
})

test_that(".rssSampleGroups splits a reference's conditions by sample group", {
    expect_equal(.rssSampleGroups(c(1, 1, 2), NULL), c(1, 1, 2))
    expect_equal(.rssSampleGroups(c(1, 1, 1), c("s1", "s2", "s2")), c(1, 2, 2))
    expect_equal(.rssSampleGroups(c(1, 2, 1), c("s1", "s1", "s1")), c(1, 2, 1))
})

test_that(".rssResidualVariances gives 1 - h2 of each condition's variance", {
    expect_equal(.rssResidualVariances(c(0.2, NA)), c(0.8, 1))
    expect_equal(.rssResidualVariances(c(0.2, 0.5), c(2, 4)), c(1.6, 2))
    expect_null(.rssResidualVariances(NULL))
})

# ---------------------------------------------------------------------------
# Data-driven mr.mash priors on the data's conditions
# ---------------------------------------------------------------------------

# A data-driven prior over four conditions, c1-c4, named by them, with a
# residual correlation.
.alignTestPrior <- function() {
    ctx <- c("c1", "c2", "c3", "c4")
    m1 <- matrix(1:16, 4, 4, dimnames = list(ctx, ctx))
    m2 <- diag(4) * 2
    dimnames(m2) <- list(ctx, ctx)
    list(
        U = list(a = m1, b = m2),
        residualCorrelation = m1 / 16,
        w = c(0.5, 0.5)
    )
}

# A prior over .rrwMulti()'s conditions, ctx1-ctx3, and a fourth, in another
# order. Its residual correlations differ, so a mix-up of the conditions
# would show.
.alignWiderPrior <- function() {
    ids <- c("ctx4", "ctx3", "ctx1", "ctx2")
    # Correlations 0.1, ..., 0.6 below the diagonal, column by column.
    lower <- replace(matrix(0, 4L, 4L), lower.tri(diag(4L)), (1:6) / 10)
    list(
        U = list(shared = `dimnames<-`(matrix(1, 4L, 4L), list(ids, ids))),
        residualCorrelation = `dimnames<-`(
            lower + t(lower) + diag(4L),
            list(ids, ids)
        )
    )
}

test_that("a named mr.mash prior is taken on the data's conditions by name", {
    prior <- .alignTestPrior()
    on <- c("c3", "c1")
    out <- pecotmr:::.mrmashAlignPriorToContexts(prior, on, 2)
    expect_equal(out$U$a, prior$U$a[on, on])
    expect_equal(out$U$b, prior$U$b[on, on])
    expect_equal(out$residualCorrelation, prior$residualCorrelation[on, on])
    expect_equal(out$w, prior$w)
    # One condition is still a matrix.
    one <- pecotmr:::.mrmashAlignPriorToContexts(prior, "c2", 1)
    expect_equal(one$U$a, prior$U$a["c2", "c2", drop = FALSE])
})

test_that("a mr.mash prior lacking one of the data's conditions is refused", {
    expect_error(
        pecotmr:::.mrmashAlignPriorToContexts(
            .alignTestPrior(),
            c("c1", "zz"),
            2
        ),
        "Condition(s) zz not found",
        fixed = TRUE
    )
})

test_that("a mr.mash prior not named on both sides must be K x K", {
    prior <- list(U = list(diag(3)))
    expect_equal(
        pecotmr:::.mrmashAlignPriorToContexts(prior, NULL, 3)$U[[1]],
        diag(3)
    )
    expect_error(
        pecotmr:::.mrmashAlignPriorToContexts(prior, NULL, 2),
        "is 3 x 3 .*: name its rows and columns by condition"
    )
    # A named prior on data whose conditions are unnamed asks for their
    # names, and is used as it is when it already has their number.
    expect_error(
        pecotmr:::.mrmashAlignPriorToContexts(.alignTestPrior(), NULL, 2),
        "name the data's conditions"
    )
    expect_equal(
        pecotmr:::.mrmashAlignPriorToContexts(.alignTestPrior(), NULL, 4),
        .alignTestPrior()
    )
})

test_that("a mr.mash prior with no U or residual correlation is refused", {
    # A bare list of matrices, as multiTraitData's priors once were: nothing
    # would read it, leaving the canonical prior alone without a word.
    bare <- list(XtX = diag(2), PCA_1 = diag(2))
    expect_error(
        pecotmr:::.mrmashAlignPriorToContexts(bare, c("c1", "c2"), 2),
        "list of covariance matrices goes under `U`",
        fixed = TRUE
    )
})

test_that("buildMrmashPriorMatrices takes a prior over more conditions", {
    skip_if_not_installed("mr.mashr")
    set.seed(1)
    ctx <- c("c2", "c4")
    Bhat <- matrix(rnorm(20), 10, 2, dimnames = list(NULL, ctx))
    Shat <- matrix(1, 10, 2, dimnames = list(NULL, ctx))
    res <- buildMrmashPriorMatrices(
        Bhat,
        Shat,
        dataDrivenPriorMatrices = .alignTestPrior(),
        canonicalPriorMatrices = FALSE
    )
    expect_equal(unique(map(res$S0, dim)), list(c(2L, 2L)))
})

test_that("mrmashWrapper takes a wider prior on Y's conditions by name", {
    skip_if_not_installed("mr.mashr")
    skip_if_not_installed("glmnet")
    m <- .rrwMulti(n = 60, p = 6, K = 3)
    prior <- .alignWiderPrior()
    # With the update off, the V mr.mash fits with is the start: the prior's
    # residual correlation on Y's conditions, scaled by their SDs.
    fit <- suppressMessages(mrmashWrapper(
        m$X,
        m$Y,
        dataDrivenPriorMatrices = prior,
        prior = MrmashPriorParam(canonicalPriorMatrices = TRUE),
        updateV = FALSE,
        methodArgs = MrmashOptions(vInitMethod = "mash", max_iter = 20)
    ))
    R <- unname(prior$residualCorrelation[colnames(m$Y), colnames(m$Y)])
    expect_equal(unname(fit$V), R * tcrossprod(apply(m$Y, 2, sd)))
})

test_that("mrmashRssWeights takes a wider prior on the z-scores' conditions", {
    .skipIfMrmashCannotFitStudies()
    m <- .rrwMulti(n = 60, p = 6, K = 3)
    prior <- .alignWiderPrior()
    w <- mrmashRssWeights(
        m$stat,
        m$ldReference,
        dataDrivenPriorMatrices = prior,
        fitRetention = "slim",
        methodArgs = MrmashRssOptions(verbose = FALSE)
    )
    # mr.mash.rss starts V from covY, the prior's residual correlation on the
    # z-scores' conditions, and does not update it.
    on <- colnames(m$stat$z)
    expect_equal(
        unname(attr(w, "fit")$V),
        unname(prior$residualCorrelation[on, on])
    )
})
