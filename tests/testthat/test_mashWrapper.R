context("mash_wrapper")

# ===========================================================================
# mashModelFit / mashPosterior and the mashr / udr engine-argument
# constructors (R/mashWrapper.R)
# ===========================================================================

delayedAssign(".mwModel", .mashTestModel(.mwSs))

test_that("mashModelFit returns a fitted mash model", {
    skip_if_not_installed("mashr")
    expect_s3_class(.mwModel, "mash")
    expect_false(is.null(.mwModel$fitted_g))
})

test_that("mashModelFit validates the prior and the fitOn entry", {
    skip_if_not_installed("mashr")
    ss <- .mwSs
    expect_error(
        mashModelFit(list(random = ss), alpha = 0, priorCovariances = list()),
        "non-empty named list"
    )
    expect_error(
        mashModelFit(
            list(strong = ss),
            alpha = 0,
            priorCovariances = list(identity = diag(3))
        ),
        "no 'random' entry"
    )
})

test_that("mashModelFit accepts either prior shape and fits identically", {
    skip_if_not_installed("mashr")
    # .mwModel is the fit of this same prior given as a bare list.
    ulist <- list(identity = diag(3), effectA = diag(c(1, 0, 0)))
    wrapped <- list(U = ulist, w = NULL, loglik = NULL)
    fit <- suppressMessages(suppressWarnings(mashModelFit(
        list(random = .mwSs),
        alpha = 0,
        priorCovariances = wrapped,
        vhat = diag(3),
        setSeed = 1L
    )))
    expect_equal(fit$loglik, .mwModel$loglik)
})

test_that("mashPosterior returns posterior matrices with covariance", {
    skip_if_not_installed("mashr")
    post <- suppressMessages(suppressWarnings(
        mashPosterior(.mwModel, .mwSs, alpha = 0, vhat = diag(3))
    ))
    expect_true(all(
        c("PosteriorMean", "PosteriorSD", "lfsr", "PosteriorCov") %in%
            names(post)
    ))
    expect_equal(ncol(post$PosteriorMean), 3L)
    expect_equal(nrow(post$PosteriorMean), nrow(post$PosteriorSD))
})

test_that("mashPosterior outputPosteriorCov = FALSE omits PosteriorCov", {
    skip_if_not_installed("mashr")
    post <- suppressMessages(suppressWarnings(
        mashPosterior(
            .mwModel,
            .mwSs,
            alpha = 0,
            vhat = diag(3),
            outputPosteriorCov = FALSE
        )
    ))
    expect_false("PosteriorCov" %in% names(post))
    expect_equal(ncol(post$PosteriorMean), 3L)
})

test_that("mashPosterior(excludeCondition) drops the condition from model + output", {
    skip_if_not_installed("mashr")
    ss <- .mwSs
    conds <- colnames(.mashSumStatsToMatrices(ss, "strong")$b)
    post <- suppressMessages(suppressWarnings(
        mashPosterior(
            .mwModel,
            ss,
            alpha = 0,
            vhat = diag(3),
            excludeCondition = conds[3]
        )
    ))
    expect_equal(ncol(post$PosteriorMean), 2L)
    expect_false(conds[3] %in% colnames(post$PosteriorMean))
})

test_that("mashPosterior errors on an unknown excludeCondition", {
    skip_if_not_installed("mashr")
    expect_error(
        mashPosterior(
            .mwModel,
            .mwSs,
            alpha = 0,
            vhat = diag(3),
            excludeCondition = "not_a_condition"
        ),
        "not found in the target"
    )
})

test_that("mashModelFit: SimpleList input + default (NULL) vhat", {
    skip_if_not_installed("mashr")
    ss <- .mwSs
    m <- suppressMessages(suppressWarnings(
        mashModelFit(
            S4Vectors::SimpleList(random = ss),
            alpha = 0,
            priorCovariances = list(identity = diag(3)),
            setSeed = 1L
        )
    ))
    expect_s3_class(m, "mash")
})

test_that("mashPosterior: default (NULL) vhat and excludeCondition dropping every condition", {
    skip_if_not_installed("mashr")
    ss <- .mwSs
    model <- .mwModel
    conds <- colnames(.mashSumStatsToMatrices(ss, "strong")$b)
    # NULL vhat path: omit vhat so mashPosterior fills the identity default.
    post <- suppressMessages(suppressWarnings(mashPosterior(
        model,
        ss,
        alpha = 0
    )))
    expect_equal(ncol(post$PosteriorMean), 3L)
    # excludeCondition removing every condition errors.
    expect_error(
        suppressMessages(suppressWarnings(
            mashPosterior(model, ss, alpha = 0, excludeCondition = conds)
        )),
        "drops every condition"
    )
})

test_that(".mashUdFit re-raises an unrelated udr failure unchanged", {
    skip_if_not_installed("udr")
    local_mocked_bindings(
        .mashUdControl = function(...) list(),
        .package = "pecotmr"
    )
    # Only the TED i.i.d. incompatibility is rewrapped; anything else must
    # surface as itself rather than being swallowed into a NULL fit.
    expect_error(
        with_mocked_bindings(
            pecotmr:::.mashUdFit(
                NULL,
                UdFitOptions(unconstrained.update = "ted"),
                2L
            ),
            ud_fit = function(...) stop("totally unrelated failure"),
            .package = "udr"
        ),
        "totally unrelated failure"
    )
})

test_that("mash covariance constructors validate against their mashr callee", {
    skip_if_not_installed("mashr")
    expect_s4_class(CovPcaOptions(subset = 1:5), "MethodOptions")
    expect_equal(
        CovFlashOptions(remove_singleton = TRUE)$remove_singleton,
        TRUE
    )
    expect_error(CovPcaOptions(subsett = 1:5), "unknown argument")
    # A genuinely unmatched name; note `cov_method` would be accepted by R's
    # partial matching, resolving to `cov_methods`.
    expect_error(
        CovCanonicalOptions(covMethods = "identity"),
        "unknown argument"
    )
    # mashr::cov_ed takes `...`, so its names cannot be checked; the record
    # says so rather than pretending otherwise.
    expect_output(show(CovEdOptions()), "NOT checked")
})

test_that("CovUdrOptions splits udr's two functions into their own bundles", {
    skip_if_not_installed("udr")
    a <- CovUdrOptions(
        init = UdInitOptions(n_unconstrained = 20L),
        fit = UdFitOptions(unconstrained.update = "ted", maxiter = 500)
    )
    expect_equal(a$init$n_unconstrained, 20L)
    expect_equal(a$fit$unconstrained.update, "ted")
    expect_equal(a$fit$maxiter, 500)
    # Each half is checked against its own live source: ud_init's formals
    # for one, ud_fit_control_default()'s names for the other. udr's tunables
    # live in a control list, not in ud_fit's formals, which is why the two
    # sources differ.
    expect_error(UdFitOptions(maxitr = 5), "unknown argument")
    expect_error(UdInitOptions(n_unconstraind = 5), "unknown argument")
    # And neither accepts the other's names, which one flat bag could not
    # distinguish.
    expect_error(UdInitOptions(maxiter = 5), "unknown argument")
    expect_error(UdFitOptions(n_unconstrained = 5), "unknown argument")
    # udr's own spellings throughout: the pecotmr relabelling is gone.
    expect_error(
        CovUdrOptions(unconstrainedUpdate = "ted"),
        "unused argument"
    )
})

test_that("the mash constructors refuse the arguments pecotmr owns", {
    expect_error(MashOptions(seed = 1L), "`seed`")
    expect_error(MashOptions(data = 1L), "the sumstats pecotmr assembles")
    expect_error(MashOptions(Ulist = list()), "the prior covariances")
    expect_error(MashOptions(outputlevel = 2L), "fixed by the calling entry")
    expect_error(MashDataOptions(Bhat = 1), "the effect-size matrix")
    expect_error(MashDataOptions(alpha = 0), "the caller's `alpha`")
    expect_error(MashDataOptions(V = diag(2)), "the caller's `vhat`")
    expect_error(MashPosteriorOptions(seed = 1L), "`seed`")
    expect_error(MashPosteriorOptions(g = 1L), "the fitted mash model")
    expect_error(
        MashPosteriorOptions(output_posterior_cov = TRUE),
        "the caller's `outputPosteriorCov`"
    )
    expect_error(MashCorSimpleOptions(data = 1L), "the null-partition data")
    expect_error(MashCorEmOptions(max_iter = 2L), "the caller's `maxIter`")
    expect_error(CorShrinkOptions(data = 1L), "the null z-matrix")
})

test_that("the mash constructors validate against their engine's formals", {
    expect_error(MashOptions(nosuch = 1), "unknown argument")
    expect_error(MashDataOptions(nosuch = 1), "unknown argument")
    expect_error(MashPosteriorOptions(nosuch = 1), "unknown argument")
    expect_error(MashCorSimpleOptions(nosuch = 1), "unknown argument")
    expect_error(CorShrinkOptions(nosuch = 1), "unknown argument")
    # mashr::mash_estimate_corr_em takes dots, so nothing can be checked.
    expect_s4_class(MashCorEmOptions(nosuch = 1), "MethodOptions")
})

test_that("the mash constructors carry their pecotmr defaults", {
    expect_equal(MashDataOptions()$zero_Bhat_Shat_reset, 1000)
    expect_setequal(names(CorShrinkOptions()), c("ash.control", "image"))
    expect_equal(CorShrinkOptions()$image, "null")
    # An explicit value wins over the default.
    expect_equal(
        MashDataOptions(zero_Bhat_Shat_reset = 1)$zero_Bhat_Shat_reset,
        1
    )
    expect_length(MashOptions(), 0L)
})

test_that("mashDataArgs and mashArgs reach mashr::mash_set_data and mash", {
    skip_if_not_installed("mashr")
    # One recorded fit serves both bundles: each mash() call is costly.
    dataSeen <- NULL
    mashSeen <- NULL
    realData <- mashr::mash_set_data
    realMash <- mashr::mash
    ss <- mashFixture(20L)
    suppressMessages(suppressWarnings(with_mocked_bindings(
        mashModelFit(
            list(random = ss),
            alpha = 0,
            priorCovariances = mashTinyPrior(),
            vhat = diag(3),
            mashDataArgs = MashDataOptions(zero_Shat_reset = 0.5),
            mashArgs = MashOptions(nullweight = 7, verbose = FALSE),
            setSeed = 1L
        ),
        mash_set_data = function(...) {
            dataSeen <<- list(...)
            realData(...)
        },
        mash = function(...) {
            mashSeen <<- list(...)
            realMash(...)
        },
        .package = "mashr"
    )))
    expect_equal(dataSeen$zero_Shat_reset, 0.5)
    # The pecotmr default rides along unless the caller overrides it.
    expect_equal(dataSeen$zero_Bhat_Shat_reset, 1000)
    expect_equal(mashSeen$nullweight, 7)
    expect_false(mashSeen$verbose)
    # The entry point keeps owning outputlevel and the data/prior pair.
    expect_equal(mashSeen$outputlevel, 4L)
})

test_that("mashPosterior forwards posteriorArgs", {
    skip_if_not_installed("mashr")
    seen <- NULL
    real <- mashr::mash_compute_posterior_matrices
    # Built before the mock is in place, so only mashPosterior is recorded.
    model <- .mwModel
    suppressMessages(suppressWarnings(with_mocked_bindings(
        mashPosterior(
            model,
            .mwSs,
            alpha = 0,
            vhat = diag(3),
            posteriorArgs = MashPosteriorOptions(pi_thresh = 1e-8)
        ),
        mash_compute_posterior_matrices = function(...) {
            seen <<- list(...)
            real(...)
        },
        .package = "mashr"
    )))
    expect_equal(seen$pi_thresh, 1e-8)
})

test_that("UdFitOptions checks unconstrained.update's value, not just its name", {
    skip_if_not_installed("udr")
    # udr validates the NAME but not the value: compute_penalty() is two
    # `if (update.type == ...)` branches with no else, so an unrecognised
    # value assigns nothing and it dies on `object 'log_penalty' not found`
    # -- a variable the caller never set. It runs whenever lambda is
    # non-zero, which .mashUdControl always makes it.
    for (v in c("ed", "ted", "none")) {
        expect_s4_class(
            UdFitOptions(unconstrained.update = v),
            "MethodOptions"
        )
    }
    # NA is udr's own default and means "choose from the data".
    expect_s4_class(UdFitOptions(unconstrained.update = NA), "MethodOptions")
    expect_equal(udr::ud_fit_control_default()$unconstrained.update, NA)

    expect_error(
        UdFitOptions(unconstrained.update = "TED"),
        "must be 'ed', 'ted', 'none' or NA"
    )
    expect_error(UdFitOptions(unconstrained.update = 1), "must be 'ed'")

    # The sibling updates have their OWN vocabularies (scaled.update is
    # fa/none), which this does not police -- guessing at them would risk
    # rejecting something udr accepts.
    expect_s4_class(UdFitOptions(scaled.update = "fa"), "MethodOptions")
    expect_s4_class(UdFitOptions(maxiter = 500), "MethodOptions")
    # Name checking is unchanged.
    expect_error(UdFitOptions(maxitr = 5), "unknown argument")
})
