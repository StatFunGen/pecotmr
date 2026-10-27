context("pvalCombine")

# ===========================================================================
# waldTestPval
# ===========================================================================

test_that("waldTestPval: returns two-sided p-values matching stats::pt", {
    beta <- c(0.5, -1.2, 0)
    se <- c(0.1, 0.4, 0.5)
    n <- 100
    res <- waldTestPval(beta, se, n)
    expected <- 2 * pt(-abs(beta / se), df = n - 2, lower.tail = TRUE)
    expect_equal(res, expected)
})

# ===========================================================================
# pvalAcat (internal)
# ===========================================================================

test_that("pvalAcat: single p-value passes through", {
    expect_equal(pecotmr:::pvalAcat(0.04), 0.04)
})

test_that("pvalAcat: returns NA when all input is NA", {
    expect_true(is.na(pecotmr:::pvalAcat(c(NA_real_, NA_real_))))
})

test_that("pvalAcat: drops NA by default", {
    with_na <- pecotmr:::pvalAcat(c(0.1, NA_real_, 0.3))
    no_na <- pecotmr:::pvalAcat(c(0.1, 0.3))
    expect_equal(with_na, no_na)
})

test_that("pvalAcat: clips very-near-1 p-values to 0.99", {
    # All p-values get pmin'd to 0.99 — a vector of 0.999s behaves like 0.99s.
    expect_equal(
        pecotmr:::pvalAcat(rep(0.999, 4)),
        pecotmr:::pvalAcat(rep(0.99, 4))
    )
})

test_that("pvalAcat: small p-values produce small combined p", {
    combined <- pecotmr:::pvalAcat(rep(1e-6, 5))
    expect_gt(combined, 0)
    expect_lt(combined, 1e-5)
})

test_that("pvalAcat: very tiny p-values use the asymptotic branch", {
    # Below 1e-15 the small-p approximation tan(pi*(0.5 - p)) ~ 1/(pi*p)
    # kicks in; the result must still be finite, in (0, 1], and smaller
    # than the input.
    combined <- pecotmr:::pvalAcat(rep(1e-20, 3))
    expect_true(is.finite(combined))
    expect_gt(combined, 0)
    expect_lt(combined, 1e-15)
})

# ===========================================================================
# combinePValues: dispatcher
# ===========================================================================

test_that("combinePValues: errors when methods argument is missing or empty", {
    expect_error(combinePValues(pvals = c(0.1, 0.2)), "methods.*is required")
    expect_error(
        combinePValues(pvals = c(0.1, 0.2), methods = character()),
        "methods.*is required"
    )
})

test_that("combinePValues: errors on unknown method", {
    expect_error(
        combinePValues(pvals = c(0.1, 0.2), methods = "bogus"),
        "Unknown method"
    )
})

test_that("combinePValues: errors when correlation-method R is missing", {
    expect_error(
        combinePValues(pvals = c(0.1, 0.2), methods = "fisher"),
        "require an `R` correlation matrix"
    )
})

test_that("combinePValues: errors when signed-z method has no zScores", {
    expect_error(
        combinePValues(pvals = c(0.1, 0.2), methods = "gbj", R = diag(2)),
        "require `zScores`"
    )
})

test_that("combinePValues: errors when neither pvals nor zScores supplied", {
    expect_error(
        combinePValues(methods = "acat"),
        "`pvals` or `zScores` must be supplied"
    )
})

test_that("combinePValues: errors when pvals and zScores lengths disagree", {
    expect_error(
        combinePValues(pvals = c(0.1, 0.2), zScores = 1, methods = "acat"),
        "must have the same length"
    )
})

test_that("combinePValues: derives p-values from z-scores when only z is given", {
    z <- c(-2.5, 1.0, 3.0)
    res <- combinePValues(zScores = z, methods = "acat")
    expected_p <- 2 * pnorm(-abs(z))
    # Internally the dispatcher derives pvals via 2 * pnorm(-|z|) and runs ACAT.
    expect_equal(res$results$acat$pval, pecotmr:::pvalAcat(expected_p))
    expect_equal(res$input$nPvalsIn, 0L)
    expect_equal(res$input$nZScoresIn, 3L)
})

test_that("combinePValues: ACAT result matches the internal pvalAcat", {
    p <- c(0.01, 0.1, 0.4)
    res <- combinePValues(pvals = p, methods = "acat")
    expect_equal(res$results$acat$method, "acat")
    expect_equal(res$results$acat$pval, pecotmr:::pvalAcat(p))
})

test_that("combinePValues: Bonferroni returns min(L * minP, 1)", {
    p <- c(0.04, 0.5, 0.8)
    res <- combinePValues(pvals = p, methods = "bonferroni")
    expect_equal(res$results$bonferroni$pval, min(length(p) * min(p), 1.0))
})

test_that("combinePValues: Bonferroni is capped at 1", {
    res <- combinePValues(pvals = c(0.9, 0.95), methods = "bonferroni")
    expect_equal(res$results$bonferroni$pval, 1.0)
})

test_that("combinePValues: runs multiple methods at once", {
    p <- c(0.01, 0.1, 0.4)
    res <- combinePValues(pvals = p, methods = c("acat", "bonferroni"))
    expect_equal(names(res$results), c("acat", "bonferroni"))
    expect_true(all(map_lgl(
        res$results,
        function(r) is.finite(r$pval)
    )))
})

test_that("combinePValues: drops NA p-values when naRm is TRUE (default)", {
    expect_warning(
        res <- combinePValues(pvals = c(0.1, NA, 0.3), methods = "acat"),
        "dropped"
    )
    expect_equal(res$input$nValid, 2L)
})

test_that("combinePValues: drops invalid (<=0, >=1, non-finite) p-values", {
    expect_warning(
        res <- combinePValues(pvals = c(0.1, 0, 1, Inf, 0.3), methods = "acat"),
        "dropped"
    )
    # 0, 1, and Inf are all invalid -> only 2 valid entries remain.
    expect_equal(res$input$nValid, 2L)
    expect_true(is.finite(res$results$acat$pval))
})

test_that("combinePValues: returns NA when no valid entries remain", {
    expect_warning(
        res <- combinePValues(pvals = c(NA_real_, NA_real_), methods = "acat"),
        "dropped"
    )
    expect_equal(res$input$nValid, 0L)
    expect_true(is.na(res$results$acat$pval))
})

test_that("combinePValues: per-method failure surfaces as NA + warning", {
    # Force an error inside the per-method dispatcher.
    local_mocked_bindings(
        pvalAcat = function(...) stop("synthetic test failure"),
        .package = "pecotmr"
    )
    expect_warning(
        res <- combinePValues(pvals = c(0.1, 0.2), methods = "acat"),
        "method 'acat' failed"
    )
    expect_true(is.na(res$results$acat$pval))
})

# ===========================================================================
# combinePValues: R-matrix alignment (.combinePvalAlignR via dispatcher)
# ===========================================================================

test_that("combinePValues: rejects non-square R", {
    R_bad <- matrix(1, nrow = 3, ncol = 2)
    expect_error(
        combinePValues(pvals = c(0.1, 0.2, 0.3), methods = "fisher", R = R_bad),
        "must be square"
    )
})

test_that("combinePValues: errors when unnamed R has wrong dimension", {
    R_bad <- diag(2)
    expect_error(
        combinePValues(pvals = c(0.1, 0.2, 0.3), methods = "fisher", R = R_bad),
        "Unnamed `R` must have nrow"
    )
})

test_that("combinePValues: errors when named R is missing entries", {
    R_named <- diag(2)
    dimnames(R_named) <- list(c("x", "y"), c("x", "y"))
    p <- c(a = 0.1, b = 0.2)
    expect_error(
        combinePValues(pvals = p, methods = "fisher", R = R_named),
        "missing entries"
    )
})

# ===========================================================================
# Dispatch branches of .combinePvalSingle, pvalAcat non-finite guard, and
# .combinePvalAlignR matrix guards
# ===========================================================================

test_that("pvalAcat returns NA when the Cauchy statistic is non-finite", {
    # naRm = FALSE keeps the NaN, which propagates to a non-finite mean and
    # trips the `!is.finite(stat)` guard.
    expect_true(is.na(pecotmr:::pvalAcat(c(0.1, 0.2, NaN), naRm = FALSE)))
})

test_that(".combinePvalAlignR rejects a non-matrix R", {
    expect_error(
        pecotmr:::.combinePvalAlignR(5, c("a", "b")),
        "Must be of type 'matrix'"
    )
})

test_that(".combinePvalAlignR rejects R whose rownames and colnames differ", {
    m <- matrix(1, 2, 2)
    rownames(m) <- c("a", "b")
    colnames(m) <- c("x", "y")
    expect_error(
        pecotmr:::.combinePvalAlignR(m, c("a", "b")),
        "rownames and colnames must be identical"
    )
})

test_that(".combinePvalSingle errors on an unknown method", {
    # Unreachable through combinePValues (the known-method guard fires first),
    # so the default switch arm is exercised directly here.
    expect_error(
        pecotmr:::.combinePvalSingle(
            "bogus",
            pvals = c(0.1, 0.2),
            zScores = NULL,
            R = NULL
        ),
        "Unknown combination method"
    )
})

test_that("combinePValues dispatches HMP through the unified menu", {
    skip_if_not_installed("harmonicmeanp")
    res <- combinePValues(pvals = c(0.01, 0.05, 0.1), methods = "hmp")
    expect_equal(res$results$hmp$method, "hmp")
    expect_true(is.finite(res$results$hmp$pval))
    expect_gt(res$results$hmp$pval, 0)
    expect_lte(res$results$hmp$pval, 1)
})

test_that("combinePValues dispatches poolr methods (fisher/stouffer/invchisq)", {
    skip_if_not_installed("poolr")
    res <- combinePValues(
        pvals = c(0.01, 0.05, 0.1),
        methods = c("fisher", "stouffer", "invchisq"),
        R = diag(3)
    )
    for (m in c("fisher", "stouffer", "invchisq")) {
        expect_true(is.finite(res$results[[m]]$pval))
        expect_gt(res$results[[m]]$pval, 0)
        expect_lt(res$results[[m]]$pval, 1)
    }
})

test_that("combinePValues dispatches the GBJ-family methods", {
    skip_if_not_installed("GBJ")
    z <- c(2.5, 1.8, 3.0)
    methods <- c("gbj", "bj", "hc", "ghc", "minp", "gbj_omni")
    res <- combinePValues(zScores = z, methods = methods, R = diag(3))
    for (m in methods) {
        expect_true(is.numeric(res$results[[m]]$pval))
        expect_gte(res$results[[m]]$pval, 0)
        expect_lte(res$results[[m]]$pval, 1)
    }
})

test_that("combinePValues dispatches aSPU and GATES", {
    skip_if_not_installed("aSPU")
    set.seed(42)
    res <- combinePValues(
        zScores = c(2.5, 1.8, 3.0),
        pvals = c(0.01, 0.05, 0.1),
        methods = c("aspu", "gates"),
        R = diag(3)
    )
    expect_gte(res$results$aspu$pval, 0)
    expect_lte(res$results$aspu$pval, 1)
    expect_gte(res$results$gates$pval, 0)
    expect_lte(res$results$gates$pval, 1)
})

# ===========================================================================
# pvalHmp
# ===========================================================================

test_that("pvalHmp returns valid p-value", {
    skip_if_not_installed("harmonicmeanp")
    pvals <- c(0.01, 0.05, 0.1)
    result <- pecotmr:::pvalHmp(pvals)
    expect_true(is.numeric(result))
    expect_true(result >= 0 && result <= 1)
    # pLandau's tail at 1 / HMP, HMP = L / sum(1 / p): for these inputs it
    # sits below both their arithmetic mean and the largest of them.
    expect_true(result < mean(pvals))
    expect_true(result < 0.1)
})


test_that("pvalHmp errors when package not available", {
    skip_if(
        requireNamespace("harmonicmeanp", quietly = TRUE),
        "harmonicmeanp is installed, cannot test missing-package path"
    )
    expect_error(pecotmr:::pvalHmp(c(0.01, 0.05)), "harmonicmeanp")
})

# =============================================================================
# pvalAcat
# =============================================================================

test_that("pvalAcat with all large p-values returns large combined p-value", {
    result <- pecotmr:::pvalAcat(c(0.8, 0.9, 0.95))
    expect_true(is.numeric(result))
    expect_true(result > 0.5 && result <= 1)
})


# =============================================================================
# pvalPoolr
# =============================================================================

test_that("generalized poolr tests reduce to the plain ones when independent", {
    # With R = I the generalized adjustment must give the unadjusted test.
    # poolr reads R as the covariance of each test's transformed statistic,
    # so the raw correlations understated Var(-2 log p) as 1 rather than 4:
    # Fisher on these four p-values gave 1.9e-6 instead of 0.0076.
    skip_if_not_installed("poolr")
    p <- c(0.01, 0.2, 0.5, 0.03)
    R <- diag(4)
    expect_equal(pecotmr:::pvalPoolr(p, "fisher", R), poolr::fisher(p)$p)
    expect_equal(pecotmr:::pvalPoolr(p, "stouffer", R), poolr::stouffer(p)$p)
    expect_equal(pecotmr:::pvalPoolr(p, "invchisq", R), poolr::invchisq(p)$p)
})

test_that("pvalPoolr runs poolr's empirical adjustment for every test", {
    # poolr names the running test from its own call; with the function
    # object inlined into the call it could not, and every empirical
    # adjustment failed inside poolr.
    skip_if_not_installed("poolr")
    p <- c(0.01, 0.2, 0.5)
    R <- diag(3)
    R[1, 2] <- R[2, 1] <- 0.3
    direct <- list(
        fisher = function() poolr::fisher(p, adjust = "empirical", R = R),
        stouffer = function() poolr::stouffer(p, adjust = "empirical", R = R),
        invchisq = function() poolr::invchisq(p, adjust = "empirical", R = R)
    )
    for (method in names(direct)) {
        set.seed(1)
        got <- pecotmr:::pvalPoolr(
            p,
            method,
            R,
            PoolrOptions(adjust = "empirical")
        )
        set.seed(1)
        expect_equal(got, direct[[method]]()$p, info = method)
    }
})

test_that(".poolrDependence converts R for the adjustment poolr runs", {
    skip_if_not_installed("poolr")
    R <- 0.6^abs(outer(1:4, 1:4, "-"))
    gen <- list(adjust = "generalized")
    expect_equal(
        pecotmr:::.poolrDependence(R, "fisher", gen),
        poolr::mvnconv(R, target = "m2lp")
    )
    expect_equal(
        pecotmr:::.poolrDependence(R, "stouffer", gen),
        poolr::mvnconv(R, target = "z")
    )
    expect_equal(
        pecotmr:::.poolrDependence(R, "invchisq", gen),
        poolr::mvnconv(R, target = "chisq1")
    )
    # The effective-number adjustments read the p-values' correlations, on
    # the sidedness the caller declared.
    expect_equal(
        pecotmr:::.poolrDependence(
            R,
            "fisher",
            list(adjust = "nyholt", side = 1)
        ),
        poolr::mvnconv(R, side = 1, target = "p", cov2cor = TRUE)
    )
    # "empirical" simulates from the statistics' own correlations.
    expect_identical(
        pecotmr:::.poolrDependence(R, "fisher", list(adjust = "empirical")),
        R
    )
    expect_identical(pecotmr:::.poolrDependence(R, "fisher", list()), R)
})

test_that("pvalPoolr errors on unknown method", {
    skip_if_not_installed("poolr")
    expect_error(
        pecotmr:::pvalPoolr(c(0.01, 0.05), method = "bogus", R = diag(2)),
        "Unknown poolr method"
    )
})

# =============================================================================
# pvalGbj
# =============================================================================

test_that("pvalGbj errors on unknown method", {
    skip_if_not_installed("GBJ")
    expect_error(
        pecotmr:::pvalGbj(c(2.5, 1.8), diag(2), method = "bogus"),
        "Unknown GBJ method"
    )
})

# =============================================================================
# pvalAspu
# =============================================================================

test_that("pvalAspu errors on unknown method", {
    skip_if_not_installed("aSPU")
    expect_error(
        pecotmr:::pvalAspu(zScores = c(1, 2), R = diag(2), method = "bogus"),
        "Unknown aSPU method"
    )
})


# =============================================================================
# Argument guards
# =============================================================================

test_that("pvalCombine: argument guards fire", {
    expect_error(waldTestPval("a", 1, 10), "beta.*Must be of type 'numeric'")
    expect_error(waldTestPval(1, "a", 10), "se.*Must be of type 'numeric'")
    expect_error(waldTestPval(1, 1, "a"), "n.*Must be of type 'numeric'")
    expect_error(
        combinePValues(pvals = 0.5, methods = "fisher", naRm = NA),
        "naRm.*May not be NA"
    )
})

test_that("the p-value constructors refuse the arguments pecotmr owns", {
    expect_error(GbjOptions(test_stats = 1), "the caller's `zScores`")
    expect_error(GbjOptions(cor_mat = diag(2)), "the caller's `R`")
    expect_error(GbjOmniOptions(test_stats = 1), "the caller's `zScores`")
    expect_error(AspuOptions(Zs = 1), "the caller's `zScores`")
    expect_error(AspuOptions(corSNP = diag(2)), "the caller's `R`")
    expect_error(PoolrOptions(p = 0.1), "the caller's `pvals`")
    expect_error(PoolrOptions(R = diag(2)), "the caller's `R`")
})

test_that("the p-value constructors check what their engines accept", {
    expect_error(GbjOptions(nosuch = 1), "unknown argument")
    expect_error(GbjOmniOptions(pairwise_cors = 1), "unknown argument")
    expect_error(AspuOptions(nosuch = 1), "unknown argument")
    # poolr's tests all end in `...`, so no name can be rejected.
    expect_s4_class(PoolrOptions(nosuch = 1), "MethodOptions")
    # The adjustment pecotmr has always used stays the default.
    expect_equal(PoolrOptions()$adjust, "generalized")
    expect_equal(PoolrOptions(adjust = "empirical")$adjust, "empirical")
})

test_that("PvalMethodsParam pairs each method with its own constructor", {
    rec <- PvalMethodsParam(
        aspu = AspuOptions(n.perm = 100),
        gbj = GbjOptions(pairwise_cors = 0.1),
        fisher = list(side = 1)
    )
    # The entries live in one named-list slot rather than ten nullable
    # ones: a record configuring aspu used to show nine NULL slots.
    expect_setequal(names(rec$methods), c("aspu", "gbj", "fisher"))
    expect_setequal(names(rec), "methods")
    # A plain list is routed through that method's constructor, so it picks
    # up the pecotmr default the constructor carries.
    expect_equal(rec$methods$fisher$adjust, "generalized")
    expect_error(
        PvalMethodsParam(gbj = AspuOptions()),
        "was built with the constructor for 'aspu'"
    )
    # The five GBJ tests share one engine, so one constructor serves them all.
    expect_s4_class(PvalMethodsParam(minp = GbjOptions()), "PvalMethodsParam")
    expect_s4_class(PvalMethodsParam(minp = GbjOptions()), "MethodParam")
    # Each entry is still an engine argument bag; only the container changed.
    expect_s4_class(rec$methods$aspu, "MethodOptions")
    expect_null(rec$methods$minp)
    # acat / hmp / bonferroni / gates are real combinePValues methods but
    # take no options, so they are not formals. The ten that do are, which
    # is why R rejects these before the constructor body runs -- the
    # hand-rolled "unknown method" check is no longer reachable.
    expect_error(PvalMethodsParam(acat = list()), "unused argument")
    expect_error(PvalMethodsParam(gates = list()), "unused argument")
    expect_setequal(
        names(formals(PvalMethodsParam)),
        names(pecotmr:::.pvalMethodCtors())
    )
})

test_that("combinePValues finds a configured method by name, not position", {
    # .combinePvalMethodArgs reads the entries out of the `methods` slot
    # and matches by name, so every configured method must resolve to its
    # own entry -- and a plain list() must still work as "no settings".
    rec <- PvalMethodsParam(
        aspu = AspuOptions(n.perm = 100),
        fisher = PoolrOptions(side = 1)
    )
    for (m in names(rec$methods)) {
        expect_identical(
            pecotmr:::.combinePvalMethodArgs(rec, m),
            rec$methods[[m]],
            info = m
        )
    }
    expect_equal(
        pecotmr:::.combinePvalMethodArgs(list(), "fisher")$adjust,
        "generalized"
    )
    # An unconfigured method falls back to its own constructor's defaults,
    # not to an empty list.
    expect_equal(
        pecotmr:::.combinePvalMethodArgs(rec, "stouffer")$adjust,
        "generalized"
    )
})

test_that("combinePValues refuses options for a method it will not run", {
    expect_error(
        combinePValues(
            pvals = c(0.01, 0.2, 0.5),
            methods = "fisher",
            R = diag(3),
            pvalMethodsParam = PvalMethodsParam(
                stouffer = PoolrOptions(side = 1)
            )
        ),
        "configures stouffer, which `methods` does not request"
    )
    expect_error(
        combinePValues(
            pvals = c(0.01, 0.2),
            methods = "acat",
            pvalMethodsParam = list(fisher = list())
        ),
        "PvalMethodsParam"
    )
})

test_that("combinePValues forwards each method's options to its engine", {
    skip_if_not_installed("poolr")
    seen <- NULL
    real <- poolr::stouffer
    with_mocked_bindings(
        combinePValues(
            pvals = c(0.01, 0.2, 0.5),
            methods = "stouffer",
            R = diag(3),
            pvalMethodsParam = PvalMethodsParam(
                stouffer = PoolrOptions(adjust = "empirical", size = 1000)
            )
        ),
        stouffer = function(...) {
            seen <<- list(...)
            # Called by its own name: poolr's empirical adjustment reads
            # which test it is from its call.
            stouffer <- real
            stouffer(...)
        },
        .package = "poolr"
    )
    expect_equal(seen$adjust, "empirical")
    expect_equal(seen$size, 1000)
})

test_that("an unconfigured method still gets its constructor's defaults", {
    skip_if_not_installed("poolr")
    seen <- NULL
    real <- poolr::fisher
    with_mocked_bindings(
        combinePValues(
            pvals = c(0.01, 0.2, 0.5),
            methods = "fisher",
            R = diag(3)
        ),
        fisher = function(...) {
            seen <<- list(...)
            real(...)
        },
        .package = "poolr"
    )
    # The adjustment pecotmr hard-coded before per-method options existed.
    expect_equal(seen$adjust, "generalized")
})

test_that("combinePValues forwards GBJ and aSPU options", {
    skip_if_not_installed("GBJ")
    skip_if_not_installed("aSPU")
    gbjSeen <- NULL
    aspuSeen <- NULL
    realGbj <- GBJ::GBJ
    realAspu <- aSPU::aSPUs
    z <- c(2.5, -0.3, 1.1)
    suppressWarnings(with_mocked_bindings(
        with_mocked_bindings(
            combinePValues(
                zScores = z,
                methods = c("gbj", "aspu"),
                R = diag(3),
                pvalMethodsParam = PvalMethodsParam(
                    gbj = GbjOptions(pairwise_cors = rep(0, 3)),
                    aspu = AspuOptions(n.perm = 200)
                )
            ),
            aSPUs = function(...) {
                aspuSeen <<- list(...)
                realAspu(...)
            },
            .package = "aSPU"
        ),
        GBJ = function(...) {
            gbjSeen <<- list(...)
            realGbj(...)
        },
        .package = "GBJ"
    ))
    expect_equal(gbjSeen$pairwise_cors, rep(0, 3))
    expect_equal(aspuSeen$n.perm, 200)
})
