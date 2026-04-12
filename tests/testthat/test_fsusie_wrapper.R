context("fsusie_wrapper")

# ---- cal_purity ----
test_that("cal_purity with min method and single element CS", {
  set.seed(42)
  X <- matrix(rnorm(100), nrow = 10, ncol = 10)
  l_cs <- list(c(1))

  result <- pecotmr:::cal_purity(l_cs, X, method = "min")
  expect_equal(result[[1]], 1)
})

test_that("cal_purity with min method and multi-element CS", {
  set.seed(42)
  X <- matrix(rnorm(200), nrow = 20, ncol = 10)
  l_cs <- list(c(1, 2, 3))

  result <- pecotmr:::cal_purity(l_cs, X, method = "min")
  expect_length(result, 1)
  # Manually compute expected: min off-diagonal |cor|
  cormat <- abs(cor(X[, c(1, 2, 3)]))
  diag(cormat) <- NA
  expect_equal(result[[1]], min(cormat, na.rm = TRUE))
})

test_that("cal_purity with non-min method returns three values", {
  set.seed(42)
  X <- matrix(rnorm(200), nrow = 20, ncol = 10)
  l_cs <- list(c(1, 2, 3))

  result <- pecotmr:::cal_purity(l_cs, X, method = "susie")
  expect_length(result[[1]], 3)  # min, mean, median
  # Manually compute expected values
  cormat <- abs(cor(X[, c(1, 2, 3)]))
  diag(cormat) <- NA
  vals <- cormat[!is.na(cormat)]
  expect_equal(result[[1]][1], min(vals))
  expect_equal(result[[1]][2], mean(vals))
  expect_equal(result[[1]][3], median(vals))
  # min <= mean and min <= median by definition
  expect_true(result[[1]][1] <= result[[1]][2])
  expect_true(result[[1]][1] <= result[[1]][3])
})

test_that("cal_purity with non-min method single element returns (1,1,1)", {
  X <- matrix(rnorm(100), nrow = 10, ncol = 10)
  l_cs <- list(c(1))

  result <- pecotmr:::cal_purity(l_cs, X, method = "susie")
  expect_equal(result[[1]], c(1, 1, 1))
})

test_that("cal_purity with multiple credible sets", {
  set.seed(42)
  X <- matrix(rnorm(200), nrow = 20, ncol = 10)
  l_cs <- list(c(1, 2), c(5, 6, 7))

  result <- pecotmr:::cal_purity(l_cs, X, method = "min")
  expect_length(result, 2)
})

# ---- fsusie_get_cs ----
# ---- fsusie_wrapper ----
test_that("fsusie_wrapper errors when fsusieR is not installed", {
  skip_if(requireNamespace("fsusieR", quietly = TRUE),
          "fsusieR is installed, skipping not-installed test")
  set.seed(1)
  X <- matrix(rnorm(50), nrow = 10, ncol = 5)
  Y <- matrix(rnorm(40), nrow = 10, ncol = 4)
  expect_error(
    fsusie_wrapper(
      X = X, Y = Y, pos = seq_len(4), L = 3, prior = "mixture_normal",
      max_SNP_EM = 100, cov_lev = 0.95, min_purity = 0.5, max_scale = 5
    ),
    "fsusieR"
  )
})

test_that("fsusie_wrapper low-purity branch sets cs to list(NULL) and cs_corr to NULL", {
  skip_if_not_installed("fsusieR")
  fake_fit <- list(
    cs = list(c(1, 2), c(3)),
    purity = c(0.1, 0.05),  # all < min_purity = 0.5
    pip = c(0.1, 0.2, 0.3, 0.05, 0.05),
    alpha = list(matrix(0.1, nrow = 2, ncol = 5), matrix(0.1, nrow = 2, ncol = 5))
  )
  local_mocked_bindings(
    susiF = function(...) fake_fit,
    .package = "fsusieR"
  )
  set.seed(1)
  X <- matrix(rnorm(50), nrow = 10, ncol = 5)
  Y <- matrix(rnorm(40), nrow = 10, ncol = 4)
  out <- fsusie_wrapper(
    X = X, Y = Y, pos = seq_len(4), L = 3, prior = "mixture_normal",
    max_SNP_EM = 100, cov_lev = 0.95, min_purity = 0.5, max_scale = 5
  )
  expect_equal(out$cs, list(NULL))
  expect_equal(out$sets$cs, list(NULL))
  expect_null(out$cs_corr)
})

test_that("fsusie_wrapper high-purity branch builds sets and computes cs_corr", {
  skip_if_not_installed("fsusieR")
  set.seed(2)
  p <- 5
  fake_fit <- list(
    cs = list(c(1, 2), c(3, 4)),
    purity = c(0.95, 0.9),  # all > min_purity = 0.5
    pip = c(0.4, 0.4, 0.6, 0.6, 0.1),
    alpha = list(
      matrix(rep(c(0.4, 0.4, 0.05, 0.05, 0.1), each = 2), nrow = 2, byrow = FALSE),
      matrix(rep(c(0.05, 0.05, 0.45, 0.4, 0.05), each = 2), nrow = 2, byrow = FALSE)
    )
  )
  local_mocked_bindings(
    susiF = function(...) fake_fit,
    cal_cor_cs = function(obj, X) matrix(c(1, 0.9, 0.9, 1), nrow = 2),
    .package = "fsusieR"
  )
  X <- matrix(rnorm(10 * p), nrow = 10, ncol = p)
  Y <- matrix(rnorm(40), nrow = 10, ncol = 4)
  out <- fsusie_wrapper(
    X = X, Y = Y, pos = seq_len(4), L = 3, prior = "mixture_normal",
    max_SNP_EM = 100, cov_lev = 0.95, min_purity = 0.5, max_scale = 5
  )
  expect_length(out$sets$cs, 2)
  expect_equal(names(out$sets$cs), c("L1", "L2"))
  expect_equal(dim(out$cs_corr), c(2, 2))
  expect_equal(out$sets$requested_coverage, 0.95)
})

test_that("fsusie_get_cs creates susie-like sets", {
  set.seed(42)
  X <- matrix(rnorm(200), nrow = 20, ncol = 10)

  fSuSiE_obj <- list(
    cs = list(c(1, 2, 3), c(5, 6)),
    alpha = list(
      c(0.4, 0.3, 0.2, 0.05, 0.02, 0.01, 0.01, 0.005, 0.003, 0.002),
      c(0.01, 0.02, 0.02, 0.05, 0.45, 0.35, 0.05, 0.02, 0.02, 0.01)
    )
  )

  result <- fsusie_get_cs(fSuSiE_obj, X, requested_coverage = 0.95)

  expect_type(result, "list")
  expect_true("cs" %in% names(result))
  expect_true("purity" %in% names(result))
  expect_true("cs_index" %in% names(result))
  expect_true("coverage" %in% names(result))
  expect_true("requested_coverage" %in% names(result))
  expect_equal(result$requested_coverage, 0.95)
  expect_equal(length(result$cs), 2)
  expect_equal(names(result$cs), c("L1", "L2"))
  # Purity should be a data.frame with min/mean/median columns
  expect_true(is.data.frame(result$purity))
  expect_equal(nrow(result$purity), 2)
  # Coverage should be numeric and positive, one per CS
  expect_length(result$coverage, 2)
  expect_true(all(result$coverage > 0 & result$coverage <= 1))
  # cs_index should identify which effects had credible sets
  expect_length(result$cs_index, 2)
})
