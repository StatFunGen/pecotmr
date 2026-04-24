context("otters")

# ---- otters_weights ----
test_that("otters_weights returns named list of weight vectors", {
  set.seed(42)
  p <- 20
  n <- 500
  z <- rnorm(p, sd = 2)
  R <- diag(p)
  sumstats <- data.frame(z = z)
  result <- otters_weights(sumstats, R, n,
    methods = list(lassosum_rss = list()),
    p_thresholds = c(0.05)
  )
  expect_type(result, "list")
  expect_true("PT_0.05" %in% names(result))
  expect_true("lassosum_rss" %in% names(result))
  expect_equal(length(result$PT_0.05), p)
  expect_equal(length(result$lassosum_rss), p)
})

test_that("otters_weights computes z from beta/se if z missing", {
  set.seed(42)
  p <- 10
  n <- 100
  sumstats <- data.frame(beta = rnorm(p, sd = 0.1), se = rep(0.05, p))
  R <- diag(p)
  result <- otters_weights(sumstats, R, n,
    methods = list(lassosum_rss = list()),
    p_thresholds = c(0.05)
  )
  expect_true("lassosum_rss" %in% names(result))
  expect_equal(length(result$lassosum_rss), p)
})

test_that("otters_weights errors when no z or beta/se", {
  sumstats <- data.frame(x = 1:5)
  expect_error(otters_weights(sumstats, diag(5), 100), "z.*beta.*se")
})

test_that("otters_weights P+T selects correct SNPs", {
  set.seed(42)
  p <- 20
  n <- 500
  # Large z-scores for first 3 SNPs (should pass threshold)
  z <- c(rep(5, 3), rep(0.1, 17))
  sumstats <- data.frame(z = z)
  R <- diag(p)
  result <- otters_weights(sumstats, R, n,
    methods = list(), p_thresholds = c(0.001)
  )
  w <- result$PT_0.001
  # First 3 should be non-zero, rest should be zero
  expect_true(all(w[1:3] != 0))
  expect_true(all(w[4:20] == 0))
})

test_that("otters_weights warns on unknown method", {
  sumstats <- data.frame(z = rnorm(5))
  expect_warning(
    otters_weights(sumstats, diag(5), 100,
      methods = list(nonexistent_method = list()),
      p_thresholds = NULL),
    "not found"
  )
})

test_that("otters_weights with multiple methods returns all", {
  set.seed(42)
  p <- 15
  n <- 500
  z <- rnorm(p, sd = 2)
  R <- diag(p)
  for (i in 1:(p - 1)) { R[i, i + 1] <- 0.3; R[i + 1, i] <- 0.3 }
  sumstats <- data.frame(z = z)
  result <- otters_weights(sumstats, R, n,
    methods = list(
      lassosum_rss = list(),
      prs_cs = list(n_iter = 50, n_burnin = 10, thin = 2, seed = 42)
    ),
    p_thresholds = c(0.001, 0.05)
  )
  expect_true(all(c("PT_0.001", "PT_0.05", "lassosum_rss", "prs_cs") %in% names(result)))
  for (nm in names(result)) {
    expect_equal(length(result[[nm]]), p)
    expect_true(all(is.finite(result[[nm]])))
  }
})

test_that("otters_weights passes correlation-scale stat fields to lassosum", {
  p <- 5
  n <- 100
  z <- rnorm(p)
  R <- diag(p)
  sumstats <- data.frame(z = z)
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    lassosum_rss_weights = function(stat, LD, ...) {
      captured$lassosum_stat <- stat
      captured$lassosum_dots <- list(...)
      rep(0.1, nrow(LD))
    },
    prs_cs_weights = function(stat, LD, ...) {
      captured$prs_cs_dots <- list(...)
      rep(0.2, nrow(LD))
    }
  )

  result <- otters_weights(
    sumstats, R, n,
    methods = list(
      lassosum_rss = list(),
      prs_cs = list(phi = 1e-4)
    ),
    p_thresholds = NULL,
    check_ld_method = NULL
  )

  expect_equal(result$lassosum_rss, rep(0.1, p))
  expect_equal(result$prs_cs, rep(0.2, p))
  expect_equal(captured$lassosum_stat$cor, z / sqrt(n))
  expect_equal(captured$lassosum_stat$z, z)
  expect_equal(captured$lassosum_stat$b, z / sqrt(n))
})

# ---- otters_association ----
test_that("otters_association returns correct structure", {
  set.seed(42)
  p <- 20
  gwas_z <- rnorm(p)
  R <- diag(p)
  weights <- list(
    method1 = rnorm(p, sd = 0.01),
    method2 = rnorm(p, sd = 0.01)
  )
  result <- otters_association(weights, gwas_z, R)
  expect_true(is.data.frame(result))
  expect_true(all(c("method", "twas_z", "twas_pval", "n_snps") %in% colnames(result)))
  # Two methods + ACAT combined
  expect_equal(nrow(result), 3)
  expect_true("ACAT_combined" %in% result$method)
})

test_that("otters_association handles all-zero weights gracefully", {
  p <- 10
  gwas_z <- rnorm(p)
  R <- diag(p)
  weights <- list(zero_method = rep(0, p), nonzero = rnorm(p, sd = 0.01))
  result <- otters_association(weights, gwas_z, R)
  zero_row <- result[result$method == "zero_method", ]
  expect_true(is.na(zero_row$twas_z))
  expect_equal(zero_row$n_snps, 0)
})

test_that("otters_association with single method has no combined row", {
  p <- 10
  gwas_z <- rnorm(p)
  R <- diag(p)
  weights <- list(only_method = rnorm(p, sd = 0.01))
  result <- otters_association(weights, gwas_z, R)
  # Only one valid p-value, so no ACAT combination
  expect_false("ACAT_combined" %in% result$method)
})

test_that("otters_association uses HMP when specified", {
  skip_if_not_installed("harmonicmeanp")
  set.seed(42)
  p <- 20
  gwas_z <- rnorm(p)
  R <- diag(p)
  weights <- list(m1 = rnorm(p, sd = 0.01), m2 = rnorm(p, sd = 0.01))
  result <- otters_association(weights, gwas_z, R, combine_method = "hmp")
  expect_true("HMP_combined" %in% result$method)
})

# ---- end-to-end integration ----
test_that("otters_weights + otters_association end-to-end on simulated data", {
  set.seed(2024)
  n_eqtl <- 500
  n_gwas <- 10000
  p <- 20

  # Simulate genotypes and eQTL
  X <- matrix(rbinom(n_eqtl * p, 2, 0.3), nrow = n_eqtl)
  beta_eqtl <- rep(0, p)
  beta_eqtl[c(3, 10)] <- c(0.3, -0.2)
  expr <- X %*% beta_eqtl + rnorm(n_eqtl)
  eqtl_z <- as.vector(cor(expr, X)) * sqrt(n_eqtl)
  R <- cor(X)

  # Simulate GWAS (gene affects trait)
  beta_gwas_gene <- 0.1
  gwas_z <- R %*% (beta_eqtl * beta_gwas_gene * sqrt(n_gwas)) + rnorm(p)

  # Stage I: train weights
  sumstats <- data.frame(z = eqtl_z)
  weights <- otters_weights(sumstats, R, n_eqtl,
    methods = list(lassosum_rss = list()),
    p_thresholds = c(0.05)
  )
  expect_true(length(weights) >= 2)

  # Stage II: test association
  result <- otters_association(weights, as.numeric(gwas_z), R)
  expect_true(is.data.frame(result))
  expect_true(nrow(result) >= 2)
  # At least one method should have a small-ish p-value (gene is truly associated)
  min_pval <- min(result$twas_pval, na.rm = TRUE)
  expect_true(min_pval < 0.5)
})
