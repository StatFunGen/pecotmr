context("pvalCombine")

# ===========================================================================
# waldTestPval
# ===========================================================================

test_that("waldTestPval: returns two-sided p-values matching stats::pt", {
  beta <- c(0.5, -1.2, 0)
  se   <- c(0.1, 0.4, 0.5)
  n    <- 100
  res  <- waldTestPval(beta, se, n)
  expected <- 2 * pt(-abs(beta / se), df = n - 2, lower.tail = TRUE)
  expect_equal(res, expected)
})

test_that("waldTestPval: returns 1 when beta is zero", {
  expect_equal(waldTestPval(0, 0.1, 50), 1, tolerance = 1e-12)
})

test_that("waldTestPval: is vectorised over beta/se", {
  res <- waldTestPval(c(1, 2, 3), c(0.5, 0.5, 0.5), 100)
  expect_length(res, 3)
  expect_true(all(res > 0 & res <= 1))
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
  no_na   <- pecotmr:::pvalAcat(c(0.1, 0.3))
  expect_equal(with_na, no_na)
})

test_that("pvalAcat: clips very-near-1 p-values to 0.99", {
  # All p-values get pmin'd to 0.99 — a vector of 0.999s behaves like 0.99s.
  expect_equal(pecotmr:::pvalAcat(rep(0.999, 4)),
               pecotmr:::pvalAcat(rep(0.99, 4)))
})

test_that("pvalAcat: small p-values produce small combined p", {
  combined <- pecotmr:::pvalAcat(rep(1e-6, 5))
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
  expect_error(combinePValues(pvals = c(0.1, 0.2)),
               "methods.*is required")
  expect_error(combinePValues(pvals = c(0.1, 0.2), methods = character()),
               "methods.*is required")
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
    combinePValues(pvals = c(0.1, 0.2), methods = "gbj",
                   R = diag(2)),
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
  expect_equal(res$results$acat$pval,
               pecotmr:::pvalAcat(expected_p))
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
  expect_true(all(vapply(res$results, function(r) is.finite(r$pval),
                         logical(1))))
})

test_that("combinePValues: drops NA p-values when naRm is TRUE (default)", {
  expect_warning(
    res <- combinePValues(pvals = c(0.1, NA, 0.3),
                          methods = "acat"),
    "dropped"
  )
  expect_equal(res$input$nValid, 2L)
})

test_that("combinePValues: drops invalid (<=0, >=1, non-finite) p-values", {
  expect_warning(
    res <- combinePValues(pvals = c(0.1, 0, 1, Inf, 0.3),
                          methods = "acat"),
    "dropped"
  )
  # 0, 1, and Inf are all invalid -> only 2 valid entries remain.
  expect_equal(res$input$nValid, 2L)
  expect_true(is.finite(res$results$acat$pval))
})

test_that("combinePValues: returns NA when no valid entries remain", {
  expect_warning(
    res <- combinePValues(pvals = c(NA_real_, NA_real_),
                          methods = "acat"),
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
# Tests migrated from test_misc.R (p-value combiners + waldTestPval)
# ===========================================================================

test_that("pvalHmp returns valid p-value", {
  skip_if_not_installed("harmonicmeanp")
  pvals <- c(0.01, 0.05, 0.1)
  result <- pecotmr:::pvalHmp(pvals)
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
  # The harmonic mean is L/sum(1/p) where L = length(unique(pvals))
  L <- length(unique(pvals))
  HMP <- L / sum(1 / pvals)
  # Result should be based on pLandau(1/HMP, ...) and be smaller than the arithmetic mean
  expect_true(result < mean(pvals))
  # Verify the result is less than the smallest individual p-value is not required,
  # but it should be in a reasonable range relative to the harmonic mean
  expect_true(result < 0.1)
})


test_that("pvalHmp uses unique p-values only", {
  skip_if_not_installed("harmonicmeanp")
  pvals <- c(0.01, 0.01, 0.05, 0.05, 0.3)
  result <- pecotmr:::pvalHmp(pvals)
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})


test_that("pvalHmp errors when package not available", {
  skip_if(requireNamespace("harmonicmeanp", quietly = TRUE),
          "harmonicmeanp is installed, cannot test missing-package path")
  expect_error(pecotmr:::pvalHmp(c(0.01, 0.05)), "harmonicmeanp")
})

# =============================================================================
# pvalAcat
# =============================================================================


test_that("pvalAcat returns single p-value unchanged", {
  expect_equal(pecotmr:::pvalAcat(0.05), 0.05)
})


test_that("pvalAcat combines multiple p-values", {
  pvals <- c(0.001, 0.01, 0.1)
  combined <- pecotmr:::pvalAcat(pvals)
  expect_true(combined > 0 && combined < 1)
})


test_that("pvalAcat with very small p-values does not return NA", {
  result <- pecotmr:::pvalAcat(c(1e-10, 1e-8, 1e-6))
  expect_true(is.numeric(result))
  expect_true(!is.na(result))
  expect_true(result > 0 && result <= 1)
})


test_that("pvalAcat with all large p-values returns large combined p-value", {
  result <- pecotmr:::pvalAcat(c(0.8, 0.9, 0.95))
  expect_true(is.numeric(result))
  expect_true(result > 0.5 && result <= 1)
})



test_that("pvalAcat uses asymptotic approximation for p < 1e-15", {
  # p-values below 1e-15 use 1/(p*pi) instead of tan()
  result <- pecotmr:::pvalAcat(c(1e-20, 1e-18, 0.01))
  expect_true(is.numeric(result))
  expect_true(result > 0 && result < 1)
  expect_false(is.na(result))
})

# =============================================================================
# pvalPoolr
# =============================================================================


test_that("pvalPoolr fisher method returns valid p-value", {
  skip_if_not_installed("poolr")
  pvals <- c(0.01, 0.05, 0.1)
  R <- diag(3)
  result <- pecotmr:::pvalPoolr(pvals, method = "fisher", R = R)
  expect_true(is.numeric(result))
  expect_true(result > 0 && result < 1)
})


test_that("pvalPoolr stouffer method returns valid p-value", {
  skip_if_not_installed("poolr")
  pvals <- c(0.01, 0.05, 0.1)
  R <- diag(3)
  result <- pecotmr:::pvalPoolr(pvals, method = "stouffer", R = R)
  expect_true(is.numeric(result))
  expect_true(result > 0 && result < 1)
})


test_that("pvalPoolr invchisq method returns valid p-value", {
  skip_if_not_installed("poolr")
  pvals <- c(0.01, 0.05, 0.1)
  R <- diag(3)
  result <- pecotmr:::pvalPoolr(pvals, method = "invchisq", R = R)
  expect_true(is.numeric(result))
  expect_true(result > 0 && result < 1)
})


test_that("pvalPoolr errors on unknown method", {
  skip_if_not_installed("poolr")
  expect_error(pecotmr:::pvalPoolr(c(0.01, 0.05), method = "bogus", R = diag(2)),
               "Unknown poolr method")
})

# =============================================================================
# pvalGbj
# =============================================================================


test_that("pvalGbj gbj method returns valid p-value", {
  skip_if_not_installed("GBJ")
  z <- c(2.5, 1.8, 3.0)
  R <- diag(3)
  result <- pecotmr:::pvalGbj(z, R, method = "gbj")
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})


test_that("pvalGbj hc method returns valid p-value", {
  skip_if_not_installed("GBJ")
  z <- c(2.5, 1.8, 3.0)
  R <- diag(3)
  result <- pecotmr:::pvalGbj(z, R, method = "hc")
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})


test_that("pvalGbj minp method returns valid p-value", {
  skip_if_not_installed("GBJ")
  z <- c(2.5, 1.8, 3.0)
  R <- diag(3)
  result <- pecotmr:::pvalGbj(z, R, method = "minp")
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})


test_that("pvalGbj bj method returns valid p-value", {
  skip_if_not_installed("GBJ")
  z <- c(2.5, 1.8, 3.0)
  R <- diag(3)
  result <- pecotmr:::pvalGbj(z, R, method = "bj")
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})


test_that("pvalGbj ghc method returns valid p-value", {
  skip_if_not_installed("GBJ")
  z <- c(2.5, 1.8, 3.0)
  R <- diag(3)
  result <- pecotmr:::pvalGbj(z, R, method = "ghc")
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})


test_that("pvalGbj gbj_omni method returns valid p-value", {
  skip_if_not_installed("GBJ")
  z <- c(2.5, 1.8, 3.0)
  R <- diag(3)
  result <- pecotmr:::pvalGbj(z, R, method = "gbj_omni")
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})


test_that("pvalGbj errors on unknown method", {
  skip_if_not_installed("GBJ")
  expect_error(pecotmr:::pvalGbj(c(2.5, 1.8), diag(2), method = "bogus"),
               "Unknown GBJ method")
})

# =============================================================================
# pvalAspu
# =============================================================================


test_that("pvalAspu aspu method returns valid p-value", {
  skip_if_not_installed("aSPU")
  set.seed(42)
  z <- c(2.5, 1.8, 3.0)
  R <- diag(3)
  result <- pecotmr:::pvalAspu(zScores = z, R = R, method = "aspu")
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})


test_that("pvalAspu gates method returns valid p-value", {
  skip_if_not_installed("aSPU")
  pvals <- c(0.01, 0.05, 0.1)
  R <- diag(3)
  result <- pecotmr:::pvalAspu(pvals = pvals, R = R, method = "gates")
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})


test_that("pvalAspu errors on unknown method", {
  skip_if_not_installed("aSPU")
  expect_error(pecotmr:::pvalAspu(zScores = c(1, 2), R = diag(2), method = "bogus"),
               "Unknown aSPU method")
})

# =============================================================================
# =============================================================================
# findValidFilePath and findValidFilePaths
# =============================================================================


test_that("waldTestPval computes correct p-values", {
  pval <- waldTestPval(beta = 5, se = 1, n = 100)
  expect_true(pval < 0.001)

  pval_zero <- waldTestPval(beta = 0, se = 1, n = 100)
  expect_equal(pval_zero, 1.0, tolerance = 1e-10)
})


test_that("waldTestPval handles vector inputs", {
  betas <- c(0, 1, 2, 5)
  ses <- c(1, 1, 1, 1)
  pvals <- waldTestPval(betas, ses, n = 100)
  expect_length(pvals, 4)
  expect_true(pvals[1] > pvals[4])
})


test_that("waldTestPval is symmetric in beta sign", {
  pval_pos <- waldTestPval(beta = 3, se = 1, n = 50)
  pval_neg <- waldTestPval(beta = -3, se = 1, n = 50)
  expect_equal(pval_pos, pval_neg, tolerance = 1e-10)
})


test_that("waldTestPval with very large beta gives p near 0", {
  pval <- waldTestPval(beta = 100, se = 1, n = 1000)
  expect_true(pval < 1e-10)
})


test_that("waldTestPval with very large se gives p near 1", {
  pval <- waldTestPval(beta = 1, se = 1000, n = 100)
  expect_true(pval > 0.99)
})

# =============================================================================
# parseRegion
# =============================================================================


