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
