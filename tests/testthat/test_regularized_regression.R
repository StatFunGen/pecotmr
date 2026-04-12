context("regularized_regression")

# ---- prs_cs ----
test_that("prs_cs errors on invalid LD input", {
  expect_error(prs_cs(bhat = rnorm(5), LD = "not_a_list", n = 100),
               "valid list of LD blocks")
})

test_that("prs_cs errors on non-positive sample size", {
  expect_error(prs_cs(bhat = rnorm(5), LD = list(blk1 = diag(5)), n = -1),
               "valid sample size")
})

test_that("prs_cs errors on mismatched maf length", {
  expect_error(
    prs_cs(bhat = rnorm(5), LD = list(blk1 = diag(5)), n = 100, maf = rep(0.3, 3)),
    "same as 'maf'"
  )
})

test_that("prs_cs errors on mismatched bhat and LD dimensions", {
  expect_error(
    prs_cs(bhat = rnorm(10), LD = list(blk1 = diag(5)), n = 100),
    "same as the sum"
  )
})

test_that("prs_cs runs successfully with valid input", {
  set.seed(42)
  p <- 10
  n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  for (i in 1:(p - 1)) {
    R[i, i + 1] <- 0.3
    R[i + 1, i] <- 0.3
  }
  result <- prs_cs(bhat = bhat, LD = list(blk1 = R), n = n,
                   maf = rep(0.3, p), n_iter = 50, n_burnin = 10, thin = 2)
  expect_type(result, "list")
  expect_true("beta_est" %in% names(result))
  expect_equal(length(result$beta_est), p)
  expect_true(all(is.finite(result$beta_est)))
})

test_that("prs_cs accepts multiple LD blocks whose dimensions sum to length of bhat", {
  set.seed(42)
  p1 <- 5
  p2 <- 5
  p <- p1 + p2
  n <- 100

  bhat <- rnorm(p, sd = 0.1)
  R1 <- diag(p1)
  R2 <- diag(p2)

  result <- prs_cs(
    bhat = bhat,
    LD = list(blk1 = R1, blk2 = R2),
    n = n,
    maf = rep(0.3, p),
    n_iter = 50, n_burnin = 10, thin = 2
  )

  expect_type(result, "list")
  expect_true("beta_est" %in% names(result))
  expect_equal(length(result$beta_est), p)
  expect_true(all(is.finite(result$beta_est)))
})

test_that("prs_cs with phi = NULL estimates phi automatically", {
  set.seed(42)
  p <- 10; n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  result <- prs_cs(bhat = bhat, LD = list(blk1 = R), n = n,
                   phi = NULL, maf = rep(0.3, p),
                   n_iter = 50, n_burnin = 10, thin = 2)
  expect_true("phi_est" %in% names(result))
})

test_that("prs_cs with explicit phi value", {
  set.seed(42)
  p <- 10; n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  result <- prs_cs(bhat = bhat, LD = list(blk1 = R), n = n,
                   phi = 0.01, maf = rep(0.3, p),
                   n_iter = 50, n_burnin = 10, thin = 2)
  expect_true("phi_est" %in% names(result))
  expect_true("sigma_est" %in% names(result))
  expect_true("psi_est" %in% names(result))
})

test_that("prs_cs works without maf (maf = NULL)", {
  set.seed(42)
  p <- 10; n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  result <- prs_cs(bhat = bhat, LD = list(blk1 = R), n = n,
                   maf = NULL, n_iter = 50, n_burnin = 10, thin = 2)
  expect_equal(length(result$beta_est), p)
})

# ---- prs_cs signal recovery ----
test_that("prs_cs recovers signal direction on simulated genotype data", {
  set.seed(42)
  n <- 500
  p <- 20
  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n)
  beta_true <- rep(0, p)
  beta_true[c(3, 10, 15)] <- c(0.4, -0.3, 0.2)
  y <- X %*% beta_true + rnorm(n)
  bhat <- as.vector(cor(y, X))
  R <- cor(X)
  result <- prs_cs(bhat = bhat, LD = list(blk1 = R), n = n,
                   n_iter = 1000, n_burnin = 500, thin = 5, seed = 42)
  expect_true("beta_est" %in% names(result))
  expect_equal(length(result$beta_est), p)
  expect_true(all(is.finite(result$beta_est)))
  # Sigma should be reasonable (near 1 for standardized data)
  expect_true(result$sigma_est > 0.1 && result$sigma_est < 10)
  # Correlation with truth should be positive (signal recovery)
  expect_gt(cor(result$beta_est, beta_true), 0.5)
})

# ---- prs_cs_weights (wrapper) ----
test_that("prs_cs_weights calls prs_cs and returns beta_est", {
  set.seed(42)
  p <- 10
  n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  for (i in 1:(p - 1)) {
    R[i, i + 1] <- 0.3
    R[i + 1, i] <- 0.3
  }
  stat <- list(b = bhat, n = rep(n, p))
  result <- prs_cs_weights(stat = stat, LD = R,
                           maf = rep(0.3, p), n_iter = 50, n_burnin = 10, thin = 2)
  expect_equal(length(result), p)
  expect_true(is.numeric(result))
})

# ---- sdpr ----
test_that("sdpr errors on mismatched bhat and LD dimensions", {
  expect_error(
    sdpr(bhat = rnorm(10), LD = list(blk1 = diag(5)), n = 100),
    "same as the length of bhat"
  )
})

test_that("sdpr errors on non-positive sample size", {
  expect_error(
    sdpr(bhat = rnorm(5), LD = list(blk1 = diag(5)), n = -1),
    "positive integer"
  )
})

test_that("sdpr errors when M is less than 4", {
  expect_error(
    sdpr(bhat = rnorm(5), LD = list(blk1 = diag(5)), n = 100, M = 3),
    "'M' must be at least 4"
  )
})

test_that("sdpr errors on invalid per_variant_sample_size", {
  expect_error(
    sdpr(bhat = rnorm(5), LD = list(blk1 = diag(5)), n = 100,
         per_variant_sample_size = c(100, -1, 100, 100, 100)),
    "positive values"
  )
})

test_that("sdpr errors on invalid array values", {
  expect_error(
    sdpr(bhat = rnorm(5), LD = list(blk1 = diag(5)), n = 100,
         array = c(0, 1, 3, 1, 0)),
    "0, 1, or 2"
  )
})

test_that("sdpr runs successfully", {
  set.seed(42)
  p <- 10
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  result <- sdpr(bhat = bhat, LD = list(blk1 = R), n = 100,
                 iter = 50, burn = 10, thin = 2, verbose = FALSE)
  expect_type(result, "list")
  expect_true("beta_est" %in% names(result))
  expect_equal(length(result$beta_est), p)
  expect_true(all(is.finite(result$beta_est)))
})

test_that("sdpr with per_variant_sample_size", {
  set.seed(42)
  p <- 10
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  result <- sdpr(bhat = bhat, LD = list(blk1 = R), n = 100,
                 per_variant_sample_size = rep(100, p),
                 iter = 50, burn = 10, thin = 2, verbose = FALSE)
  expect_equal(length(result$beta_est), p)
  expect_true(all(is.finite(result$beta_est)))
})

test_that("sdpr with valid array parameter", {
  set.seed(42)
  p <- 10
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  result <- sdpr(bhat = bhat, LD = list(blk1 = R), n = 100,
                 array = rep(1, p),
                 iter = 50, burn = 10, thin = 2, verbose = FALSE)
  expect_equal(length(result$beta_est), p)
  expect_true(all(is.finite(result$beta_est)))
})

# ---- sdpr signal recovery ----
test_that("sdpr recovers signal direction on simulated genotype data", {
  set.seed(2024)
  n <- 500
  p <- 20
  # Realistic genotype matrix with non-trivial LD (binomial, MAF=0.3)
  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n)
  beta_true <- rep(0, p)
  beta_true[c(3, 10, 15)] <- c(0.4, -0.3, 0.2)
  y <- X %*% beta_true + rnorm(n)
  bhat <- as.vector(cor(y, X))
  R <- cor(X)
  result <- sdpr(bhat = bhat, LD = list(blk1 = R), n = n,
                 iter = 500, burn = 200, thin = 5, verbose = FALSE)
  expect_true("beta_est" %in% names(result))
  expect_equal(length(result$beta_est), p)
  expect_true(all(is.finite(result$beta_est)))
  # Correlation with truth should be positive (signal recovery)
  expect_gt(cor(result$beta_est, beta_true), 0.5)
})

test_that("sdpr accepts multiple LD blocks with realistic genotype data", {
  set.seed(2024)
  n <- 500
  p <- 20
  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n)
  beta_true <- rep(0, p)
  beta_true[c(3, 15)] <- c(0.4, 0.2)
  y <- X %*% beta_true + rnorm(n)
  bhat <- as.vector(cor(y, X))
  R <- cor(X)
  R1 <- R[1:10, 1:10]
  R2 <- R[11:20, 11:20]
  result <- sdpr(bhat = bhat, LD = list(blk1 = R1, blk2 = R2), n = n,
                 iter = 500, burn = 200, thin = 5, verbose = FALSE)
  expect_equal(length(result$beta_est), p)
  expect_true(all(is.finite(result$beta_est)))
})

# ---- sdpr_weights (wrapper) ----
test_that("sdpr_weights calls sdpr and returns beta_est", {
  set.seed(42)
  p <- 10
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  stat <- list(b = bhat, n = rep(100, p))
  result <- sdpr_weights(stat = stat, LD = R,
                         iter = 50, burn = 10, thin = 2, verbose = FALSE)
  expect_equal(length(result), p)
  expect_true(is.numeric(result))
})

# ---- init_prior_sd ----
test_that("init_prior_sd returns correct number of values with expected properties", {
  set.seed(123)
  n <- 50
  p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 2 + rnorm(n)

  result <- init_prior_sd(X, y)
  expect_length(result, 30)
  expect_equal(result[1], 0)
  expect_true(all(diff(result) > 0))

  result_15 <- init_prior_sd(X, y, n = 15)
  expect_length(result_15, 15)
  expect_equal(result_15[1], 0)
  expect_true(all(diff(result_15) > 0))
})

# ---- susie_weights ----
test_that("susie_weights returns zeros when no alpha/mu/X_column_scale_factors", {
  mock_fit <- list(pip = c(0.1, 0.2, 0.3))
  result <- susie_weights(susie_fit = mock_fit)
  expect_equal(result, rep(0, 3))
})

test_that("susie_weights errors on dimension mismatch", {
  mock_fit <- list(pip = c(0.1, 0.2))
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  expect_error(susie_weights(X = X, susie_fit = mock_fit), "Dimension mismatch")
})

test_that("susie_weights with alpha/mu/X_column_scale_factors calls coef.susie", {
  p <- 5
  L <- 2
  mock_fit <- list(
    alpha = matrix(c(
      0.9, 0.05, 0.02, 0.02, 0.01,
      0.1, 0.1, 0.6, 0.1, 0.1
    ), nrow = L, ncol = p, byrow = TRUE),
    mu = matrix(c(
      2.0, 0.1, 0.05, 0.03, 0.01,
      0.5, 0.2, 1.5, 0.1, 0.05
    ), nrow = L, ncol = p, byrow = TRUE),
    X_column_scale_factors = rep(1.0, p),
    pip = runif(p),
    intercept = 0
  )
  result <- susie_weights(susie_fit = mock_fit)
  expect_length(result, p)
  expect_true(is.numeric(result))
  expect_true(any(result != 0))
})

test_that("susie_weights calls susie_wrapper when susie_fit is NULL", {
  set.seed(42)
  p <- 5
  n <- 50
  X <- matrix(rnorm(n * p), nrow = n)
  y <- rnorm(n)
  local_mocked_bindings(
    susie_wrapper = function(...) {
      list(pip = rep(0.1, p))
    }
  )
  result <- susie_weights(X = X, y = y)
  expect_equal(result, rep(0, p))
})

# ---- susie_ash_weights ----
test_that("susie_ash_weights returns zeros when no expected fields", {
  mock_fit <- list(pip = c(0.1, 0.2, 0.3))
  result <- susie_ash_weights(susie_ash_fit = mock_fit)
  expect_equal(result, rep(0, 3))
})

test_that("susie_ash_weights errors on dimension mismatch", {
  mock_fit <- list(pip = c(0.1, 0.2))
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  expect_error(susie_ash_weights(X = X, susie_ash_fit = mock_fit), "Dimension mismatch")
})

test_that("susie_ash_weights with proper fields calls coef.susie", {
  p <- 4
  L <- 2
  mock_fit <- list(
    alpha = matrix(runif(L * p), nrow = L, ncol = p),
    mu = matrix(rnorm(L * p), nrow = L, ncol = p),
    theta = matrix(rnorm(L * p), nrow = L, ncol = p),
    X_column_scale_factors = rep(1.0, p),
    pip = runif(p),
    intercept = 0
  )
  result <- susie_ash_weights(susie_ash_fit = mock_fit)
  expect_true(is.numeric(result))
  expect_true(length(result) >= p)
})

test_that("susie_ash_weights calls susie_wrapper when fit is NULL", {
  set.seed(42)
  p <- 4
  n <- 30
  X <- matrix(rnorm(n * p), nrow = n)
  y <- rnorm(n)
  local_mocked_bindings(
    susie_wrapper = function(...) {
      list(pip = rep(0.1, p))
    }
  )
  result <- susie_ash_weights(X = X, y = y)
  expect_equal(result, rep(0, p))
})

# ---- susie_inf_weights ----
test_that("susie_inf_weights returns zeros when no expected fields", {
  mock_fit <- list(pip = c(0.4, 0.5))
  result <- susie_inf_weights(susie_inf_fit = mock_fit)
  expect_equal(result, rep(0, 2))
})

test_that("susie_inf_weights errors on dimension mismatch", {
  mock_fit <- list(pip = c(0.1, 0.2))
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  expect_error(susie_inf_weights(X = X, susie_inf_fit = mock_fit), "Dimension mismatch")
})

test_that("susie_inf_weights with proper fields calls coef.susie", {
  p <- 4
  L <- 2
  mock_fit <- list(
    alpha = matrix(runif(L * p), nrow = L, ncol = p),
    mu = matrix(rnorm(L * p), nrow = L, ncol = p),
    theta = matrix(rnorm(L * p), nrow = L, ncol = p),
    X_column_scale_factors = rep(1.0, p),
    pip = runif(p),
    intercept = 0
  )
  result <- susie_inf_weights(susie_inf_fit = mock_fit)
  expect_true(is.numeric(result))
  expect_true(length(result) >= p)
})

test_that("susie_inf_weights calls susie_wrapper when fit is NULL", {
  set.seed(42)
  p <- 4
  n <- 30
  X <- matrix(rnorm(n * p), nrow = n)
  y <- rnorm(n)
  local_mocked_bindings(
    susie_wrapper = function(...) {
      list(pip = rep(0.1, p))
    }
  )
  result <- susie_inf_weights(X = X, y = y)
  expect_equal(result, rep(0, p))
})

# ---- mrmash_weights ----
test_that("mrmash_weights errors when mr.mashr package is not available", {
  skip_if(requireNamespace("mr.mashr", quietly = TRUE),
          "mr.mashr is installed; skipping missing-package test")

  expect_error(
    mrmash_weights(mrmash_fit = NULL, X = matrix(1, 10, 5), Y = matrix(1, 10, 3)),
    "mr\\.mash\\.alpha"
  )
})

test_that("mrmash_weights errors when X and Y are NULL and fit is NULL", {
  skip_if_not(requireNamespace("mr.mashr", quietly = TRUE),
              "mr.mashr not installed")
  expect_error(mrmash_weights(mrmash_fit = NULL, X = NULL, Y = NULL),
               "Both X and Y must be provided")
})

# ---- mvsusie_weights ----
test_that("mvsusie_weights errors when mvsusieR package is not available", {
  skip_if(requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR is installed; skipping missing-package test")

  expect_error(
    mvsusie_weights(mvsusie_fit = NULL, X = matrix(1, 10, 5), Y = matrix(1, 10, 3)),
    "mvsusieR"
  )
})

test_that("mvsusie_weights errors when X and Y are NULL and fit is NULL", {
  skip_if_not(requireNamespace("mvsusieR", quietly = TRUE),
              "mvsusieR not installed")
  expect_error(mvsusie_weights(mvsusie_fit = NULL, X = NULL, Y = NULL),
               "Both X and Y must be provided")
})

# ---- glmnet_weights ----
test_that("glmnet_weights computes LASSO weights", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- glmnet_weights(X, y, alpha = 1)
  expect_equal(nrow(result), p)
  expect_equal(ncol(result), 1)
})

test_that("enet_weights computes elastic net weights", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- enet_weights(X, y)
  expect_equal(nrow(result), p)
})

test_that("lasso_weights computes LASSO weights", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- lasso_weights(X, y)
  expect_equal(nrow(result), p)
})

test_that("glmnet_weights handles zero-variance columns", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  X[, 3] <- 5
  y <- X[, 1] * 0.5 + rnorm(n)
  expect_warning(result <- glmnet_weights(X, y, alpha = 1),
                 "glmnet_weights: dropping 1 zero-variance column")
  expect_equal(result[3, 1], 0)
})

test_that("glmnet_weights errors when all columns are constant", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50
  p <- 5
  X <- matrix(rep(1:p, each = n), nrow = n, ncol = p)
  y <- rnorm(n)

  expect_error(glmnet_weights(X, y, alpha = 1), "matrix with 2 or more columns")
})

test_that("glmnet_weights handles NA-producing columns gracefully", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50
  p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)

  X[1, 4] <- NA

  expect_warning(result <- glmnet_weights(X, y, alpha = 1),
                 "glmnet_weights: dropping 1 zero-variance column")
  expect_equal(nrow(result), p)
  expect_equal(ncol(result), 1)
  expect_equal(result[4, 1], 0)
})

test_that("glmnet_weights with alpha = 0 (ridge regression) works", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- glmnet_weights(X, y, alpha = 0)
  expect_equal(nrow(result), p)
})

# ---- bayes_alphabet_weights ----
test_that("bayes_alphabet_weights errors on dimension mismatch", {
  skip_if_not_installed("qgg")
  X <- matrix(rnorm(100), nrow = 10)
  y <- rnorm(5)
  expect_error(bayes_alphabet_weights(X, y, method = "bayesN"), "same number of rows")
})

test_that("bayes_alphabet_weights errors on covariate dimension mismatch", {
  skip_if_not_installed("qgg")
  X <- matrix(rnorm(100), nrow = 10)
  y <- rnorm(10)
  Z <- matrix(rnorm(15), nrow = 5)
  expect_error(bayes_alphabet_weights(X, y, method = "bayesN", Z = Z),
               "same number of rows")
})

test_that("bayes_alphabet_weights does not error on valid Z dimensions", {
  skip_if_not_installed("qgg")
  set.seed(42)
  n <- 50
  p <- 10
  q <- 3
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  Z <- matrix(round(runif(n * q, 0, 0.8), 0), nrow = n)

  result <- bayes_alphabet_weights(X, y, method = "bayesN", Z = Z, nit = 50, nburn = 10)
  expect_equal(length(result), p)
  expect_true(is.numeric(result))
})

test_that("bayes_alphabet_weights with Z = NULL does not error on dimension check", {
  skip_if_not_installed("qgg")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- bayes_alphabet_weights(X, y, method = "bayesN", Z = NULL, nit = 50, nburn = 10)
  expect_equal(length(result), p)
  expect_true(is.numeric(result))
})

# ---- bayes_n/l/a/c/r_weights (wrapper dispatchers) ----
test_that("bayes_n_weights dispatches to bayes_alphabet_weights with bayesN", {
  skip_if_not_installed("qgg")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- bayes_n_weights(X, y, nit = 50, nburn = 10)
  expect_equal(length(result), p)
  expect_true(is.numeric(result))
  expect_true(all(is.finite(result)))
})

# # ---- gbayes_rss ----
# test_that("gbayes_rss errors without LD matrix", {
#   skip_if_not_installed("qgg")
#   sumstats <- data.frame(beta = 0.1, se = 0.05, n = 1000)
#   expect_error(gbayes_rss(sumstats = sumstats), "Must provide LD")
# })
# 
# test_that("gbayes_rss errors on invalid method", {
#   skip_if_not_installed("qgg")
#   sumstats <- data.frame(beta = 0.1, se = 0.05, n = 1000)
#   LD <- diag(1)
#   expect_error(gbayes_rss(sumstats = sumstats, LD = LD, method = "invalidMethod"),
#                "not valid")
# })
# 
# test_that("gbayes_rss errors on non-dataframe input", {
#   skip_if_not_installed("qgg")
#   # A list input fails before the is.data.frame check because nrow(list) is NULL,
#   # causing the dimension comparison to error with "argument is of length zero"
#   expect_error(gbayes_rss(sumstats = list(beta = 0.1), LD = diag(1)),
#                "argument is of length zero")
# })
# 
# test_that("gbayes_rss errors on dimension mismatch", {
#   skip_if_not_installed("qgg")
#   sumstats <- data.frame(beta = c(0.1, 0.2), se = c(0.05, 0.1), n = c(1000, 1000))
#   LD <- diag(3)
#   expect_error(gbayes_rss(sumstats = sumstats, LD = LD), "must correspond")
# })
# 
# test_that("gbayes_rss rejects invalid method name", {
#   skip_if_not_installed("qgg")
#   sumstats <- data.frame(
#     beta = c(0.1, 0.2, 0.3),
#     se = c(0.05, 0.05, 0.05),
#     n = c(1000, 1000, 1000)
#   )
#   LD <- diag(3)
#   expect_error(
#     gbayes_rss(sumstats = sumstats, LD = LD, method = "bayesZ"),
#     "not valid"
#   )
# })
# 
# test_that("gbayes_rss runs successfully with valid method", {
#   skip_if_not_installed("qgg")
#   skip_if_not(exists("sbayes_spa", where = asNamespace("qgg")),
#               "qgg:::sbayes_spa not available in installed qgg version")
#   sumstats <- data.frame(
#     beta = c(0.1, 0.2, 0.3),
#     se = c(0.05, 0.05, 0.05),
#     n = c(1000, 1000, 1000)
#   )
#   LD <- diag(3)
#   result <- gbayes_rss(sumstats = sumstats, LD = LD, method = "bayesN", nit = 10, nburn = 0)
#   expect_true("method" %in% names(result))
# })
# 
# test_that("gbayes_rss errors when sumstats has NA in wy", {
#   skip_if_not_installed("qgg")
#   sumstats <- data.frame(
#     beta = c(0.1, NA),
#     se = c(0.05, 0.05),
#     n = c(1000, 1000)
#   )
#   LD <- diag(2)
#   expect_error(
#     gbayes_rss(sumstats = sumstats, LD = LD, method = "bayesN"),
#     "Missing values"
#   )
# })
# 
# test_that("gbayes_rss sets variant_ids from data when not provided", {
#   skip_if_not_installed("qgg")
#   skip_if_not_installed("coda")
#   skip_if_not(exists("sbayes_spa", where = asNamespace("qgg")),
#               "qgg:::sbayes_spa not available in installed qgg version")
#   sumstats <- data.frame(
#     beta = c(0.1, 0.2, 0.3),
#     se = c(0.05, 0.05, 0.05),
#     n = c(1000, 1000, 1000)
#   )
#   LD <- diag(3)
#   result <- gbayes_rss(sumstats = sumstats, LD = LD, method = "bayesN",
#                        nit = 20, nburn = 5)
#   expect_true("sumstats" %in% names(result))
#   expect_equal(result$sumstats$variant_id[1], "snp1")
# })
# 
# test_that("gbayes_rss uses rsids column when available", {
#   skip_if_not_installed("qgg")
#   skip_if_not_installed("coda")
#   skip_if_not(exists("sbayes_spa", where = asNamespace("qgg")),
#               "qgg:::sbayes_spa not available in installed qgg version")
#   sumstats <- data.frame(
#     rsids = c("rs1", "rs2", "rs3"),
#     beta = c(0.1, 0.2, 0.3),
#     se = c(0.05, 0.05, 0.05),
#     n = c(1000, 1000, 1000)
#   )
#   LD <- diag(3)
#   result <- gbayes_rss(sumstats = sumstats, LD = LD, method = "bayesN",
#                        nit = 20, nburn = 5)
#   expect_true("sumstats" %in% names(result))
# })
# 
# test_that("gbayes_rss uses explicit variant_ids parameter", {
#   skip_if_not_installed("qgg")
#   skip_if_not_installed("coda")
#   skip_if_not(exists("sbayes_spa", where = asNamespace("qgg")),
#               "qgg:::sbayes_spa not available in installed qgg version")
#   sumstats <- data.frame(
#     beta = c(0.1, 0.2, 0.3),
#     se = c(0.05, 0.05, 0.05),
#     n = c(1000, 1000, 1000)
#   )
#   LD <- diag(3)
#   result <- gbayes_rss(sumstats = sumstats, LD = LD, method = "bayesN",
#                        variant_ids = c("my_snp1", "my_snp2", "my_snp3"),
#                        nit = 20, nburn = 5)
#   expect_true("sumstats" %in% names(result))
# })
# 
# test_that("gbayes_rss em-mcmc algorithm path sets algo = 2", {
#   skip_if_not_installed("qgg")
#   skip_if_not_installed("coda")
#   skip_if_not(exists("sbayes_spa", where = asNamespace("qgg")),
#               "qgg:::sbayes_spa not available in installed qgg version")
#   sumstats <- data.frame(
#     beta = c(0.1, 0.2, 0.3),
#     se = c(0.05, 0.05, 0.05),
#     n = c(1000, 1000, 1000)
#   )
#   LD <- diag(3)
#   result <- gbayes_rss(sumstats = sumstats, LD = LD, method = "bayesN",
#                        algorithm = "em-mcmc", nit = 20, nburn = 5)
#   expect_true("method" %in% names(result))
# })
# 
# test_that("gbayes_rss handles sumstats with maf column", {
#   skip_if_not_installed("qgg")
#   skip_if_not_installed("coda")
#   skip_if_not(exists("sbayes_spa", where = asNamespace("qgg")),
#               "qgg:::sbayes_spa not available in installed qgg version")
#   sumstats <- data.frame(
#     variant_id = c("rs1", "rs2", "rs3"),
#     beta = c(0.1, 0.2, 0.3),
#     se = c(0.05, 0.05, 0.05),
#     n = c(1000, 1000, 1000),
#     maf = c(0.1, 0.2, 0.3),
#     pos = c(100, 200, 300),
#     A1 = c("A", "C", "G"),
#     A2 = c("T", "G", "A")
#   )
#   LD <- diag(3)
#   result <- gbayes_rss(sumstats = sumstats, LD = LD, method = "bayesR",
#                        nit = 30, nburn = 10)
#   expect_true("sumstats" %in% names(result))
#   expect_true("vm" %in% colnames(result$sumstats))
# })
# 
# test_that("gbayes_rss handles sumstats without A1/A2/maf columns", {
#   skip_if_not_installed("qgg")
#   skip_if_not_installed("coda")
#   skip_if_not(exists("sbayes_spa", where = asNamespace("qgg")),
#               "qgg:::sbayes_spa not available in installed qgg version")
#   sumstats <- data.frame(
#     beta = c(0.1, 0.2, 0.3),
#     se = c(0.05, 0.05, 0.05),
#     n = c(1000, 1000, 1000)
#   )
#   LD <- diag(3)
#   result <- gbayes_rss(sumstats = sumstats, LD = LD, method = "bayesN",
#                        nit = 20, nburn = 5)
#   expect_true("sumstats" %in% names(result))
#   expect_equal(result$sumstats$A1[1], "Unknown")
#   expect_equal(result$sumstats$A2[1], "Unknown")
# })
# 
# test_that("gbayes_rss bayesR method uses 4-component mixture pi", {
#   skip_if_not_installed("qgg")
#   skip_if_not_installed("coda")
#   skip_if_not(exists("sbayes_spa", where = asNamespace("qgg")),
#               "qgg:::sbayes_spa not available in installed qgg version")
#   sumstats <- data.frame(
#     beta = c(0.1, 0.2, 0.3),
#     se = c(0.05, 0.05, 0.05),
#     n = c(1000, 1000, 1000)
#   )
#   LD <- diag(3)
#   result <- gbayes_rss(sumstats = sumstats, LD = LD, method = "bayesR",
#                        nit = 20, nburn = 5)
#   expect_true("method" %in% names(result))
#   expect_equal(result$method, "bayesR")
# })
# 
# # ---- bayes_alphabet_rss_weights wrappers ----
# test_that("bayes_alphabet_rss_weights dispatches to gbayes_rss", {
#   skip_if_not_installed("qgg")
#   skip_if_not_installed("coda")
#   skip_if_not(exists("sbayes_spa", where = asNamespace("qgg")),
#               "qgg:::sbayes_spa not available in installed qgg version")
#   sumstats <- data.frame(
#     beta = c(0.1, 0.2, 0.3),
#     se = c(0.05, 0.05, 0.05),
#     n = c(1000, 1000, 1000)
#   )
#   LD <- diag(3)
#   result <- bayes_alphabet_rss_weights(sumstats, LD, method = "bayesN", nit = 20, nburn = 5)
#   expect_equal(length(result), 3)
#   expect_true(is.numeric(result))
# })
# 
# test_that("bayes_n_rss_weights dispatches correctly", {
#   skip_if_not_installed("qgg")
#   skip_if_not_installed("coda")
#   skip_if_not(exists("sbayes_spa", where = asNamespace("qgg")),
#               "qgg:::sbayes_spa not available in installed qgg version")
#   sumstats <- data.frame(
#     beta = c(0.1, 0.2), se = c(0.05, 0.05), n = c(1000, 1000)
#   )
#   LD <- diag(2)
#   result <- bayes_n_rss_weights(sumstats, LD, nit = 20, nburn = 5)
#   expect_true(is.numeric(result))
#   expect_equal(length(result), nrow(sumstats))
# })

# ---- mrash_weights ----
test_that("mrash_weights calls lasso_weights as default beta.init", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  local_mocked_bindings(
    lasso_weights = function(X, y) rep(0.01, ncol(X)),
    init_prior_sd = function(X, y, n = 30) seq(0, 3, length.out = n)
  )
  result <- mrash_weights(X, y)
  expect_true(is.numeric(result))
  expect_equal(length(result), p)
})

test_that("mrash_weights with init_prior_sd = FALSE passes NULL sa2", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  local_mocked_bindings(
    lasso_weights = function(X, y) rep(0.01, ncol(X))
  )
  result <- mrash_weights(X, y, init_prior_sd = FALSE)
  expect_true(is.numeric(result))
  expect_equal(length(result), p)
})

# ---- lassosum_rss ----
test_that("lassosum_rss errors on invalid LD input", {
  expect_error(lassosum_rss(bhat = rnorm(5), LD = "not_a_list", n = 100),
               "valid list of LD blocks")
})

test_that("lassosum_rss errors on non-positive sample size", {
  expect_error(lassosum_rss(bhat = rnorm(5), LD = list(blk1 = diag(5)), n = -1),
               "valid sample size")
})

test_that("lassosum_rss errors on mismatched bhat and LD dimensions", {
  expect_error(
    lassosum_rss(bhat = rnorm(10), LD = list(blk1 = diag(5)), n = 100),
    "same as the sum"
  )
})

test_that("lassosum_rss runs successfully with valid input", {
  set.seed(42)
  p <- 10
  n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  for (i in 1:(p - 1)) {
    R[i, i + 1] <- 0.3
    R[i + 1, i] <- 0.3
  }
  result <- lassosum_rss(bhat = bhat, LD = list(blk1 = R), n = n)
  expect_type(result, "list")
  expect_true("beta_est" %in% names(result))
  expect_equal(length(result$beta_est), p)
  expect_true(all(is.finite(result$beta_est)))
  expect_equal(nrow(result$beta), p)
  expect_equal(ncol(result$beta), 20)
})

test_that("lassosum_rss accepts multiple LD blocks", {
  set.seed(42)
  p1 <- 5
  p2 <- 5
  p <- p1 + p2
  n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R1 <- diag(p1)
  R2 <- diag(p2)
  result <- lassosum_rss(bhat = bhat, LD = list(blk1 = R1, blk2 = R2), n = n)
  expect_type(result, "list")
  expect_true("beta_est" %in% names(result))
  expect_equal(length(result$beta_est), p)
  expect_true(all(is.finite(result$beta_est)))
})

test_that("lassosum_rss with large lambda gives all-zero betas", {
  set.seed(42)
  p <- 10
  n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  result <- lassosum_rss(bhat = bhat, LD = list(blk1 = R), n = n,
                         lambda = c(100))
  expect_true(all(result$beta_est == 0))
})

# ---- lassosum_rss_weights (wrapper) ----
test_that("lassosum_rss_weights calls lassosum_rss and returns beta_est", {
  set.seed(42)
  p <- 10
  n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  for (i in 1:(p - 1)) {
    R[i, i + 1] <- 0.3
    R[i + 1, i] <- 0.3
  }
  stat <- list(b = bhat, n = rep(n, p))
  result <- lassosum_rss_weights(stat = stat, LD = R)
  expect_equal(length(result), p)
  expect_true(is.numeric(result))
})

test_that("lassosum_rss_weights rescales bhat when max(abs(bhat)) >= 1", {
  set.seed(42)
  p <- 5
  n <- 100
  # Construct bhat with max abs >= 1 to trigger the clamp branch
  bhat <- c(1.5, -0.2, 0.1, 0.05, -0.3)
  R <- diag(p)
  stat <- list(b = bhat, n = rep(n, p))
  captured <- NULL
  local_mocked_bindings(
    lassosum_rss = function(bhat, LD, n, ...) {
      captured <<- bhat
      list(beta_est = rep(0, length(bhat)), fbeta = 1)
    }
  )
  result <- lassosum_rss_weights(stat = stat, LD = R, s = 0.5)
  # Captured bhat should have been rescaled so max abs is just under 1
  expect_true(max(abs(captured)) < 1)
  expect_equal(max(abs(captured)), 0.9999, tolerance = 1e-9)
  # Sign-preserving rescale
  expect_true(all(sign(captured) == sign(bhat)))
})

test_that("mrash_weights subsets a user-supplied beta.init of length ncol(X) by keep", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  X[, 5] <- 7  # zero-variance column to force subset by keep
  y <- X[, 1] * 0.5 + rnorm(n)
  user_beta_init <- seq(0.01, by = 0.01, length.out = p)
  expect_warning(
    result <- mrash_weights(X, y, beta.init = user_beta_init),
    "mrash_weights: dropping 1 zero-variance column"
  )
  expect_equal(length(result), p)
  expect_equal(result[5], 0)
})

test_that("mrash_weights warns and zero-pads when X has zero-variance columns", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  X[, 3] <- 7
  y <- X[, 1] * 0.5 + rnorm(n)
  local_mocked_bindings(
    lasso_weights = function(X, y) rep(0.01, ncol(X))
  )
  expect_warning(result <- mrash_weights(X, y),
                 "mrash_weights: dropping 1 zero-variance column")
  expect_equal(length(result), p)
  expect_equal(result[3], 0)
})

# ---- bayes_alphabet_weights zero-variance handling ----
test_that("bayes_alphabet_weights warns and zero-pads when X has zero-variance columns", {
  skip_if_not_installed("qgg")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  X[, 4] <- 2
  y <- X[, 1] * 0.5 + rnorm(n)
  expect_warning(result <- bayes_alphabet_weights(X, y, method = "bayesN", nit = 50, nburn = 10),
                 "bayes_alphabet_weights: dropping 1 zero-variance column")
  expect_equal(length(result), p)
  expect_equal(result[4], 0)
  expect_true(all(is.finite(result)))
})

# ---- ncvreg_weights / scad_weights / mcp_weights ----
test_that("ncvreg_weights computes weights with SCAD penalty", {
  skip_if_not_installed("ncvreg")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- ncvreg_weights(X, y, penalty = "SCAD")
  expect_equal(nrow(result), p)
  expect_equal(ncol(result), 1)
  expect_true(all(is.finite(result)))
})

test_that("ncvreg_weights warns and zero-pads when X has zero-variance columns", {
  skip_if_not_installed("ncvreg")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  X[, 5] <- 3
  y <- X[, 1] * 0.5 + rnorm(n)
  expect_warning(result <- ncvreg_weights(X, y, penalty = "SCAD"),
                 "ncvreg_weights: dropping 1 zero-variance column")
  expect_equal(nrow(result), p)
  expect_equal(result[5, 1], 0)
})

test_that("scad_weights computes weights and dispatches to ncvreg_weights", {
  skip_if_not_installed("ncvreg")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- scad_weights(X, y)
  expect_equal(nrow(result), p)
  expect_equal(ncol(result), 1)
  expect_true(all(is.finite(result)))
})

test_that("mcp_weights computes weights and dispatches to ncvreg_weights", {
  skip_if_not_installed("ncvreg")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- mcp_weights(X, y)
  expect_equal(nrow(result), p)
  expect_equal(ncol(result), 1)
  expect_true(all(is.finite(result)))
})

test_that("scad_weights passes nfolds through to cv.ncvreg", {
  skip_if_not_installed("ncvreg")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  # cv.ncvreg requires nfolds >= 2; passing 1 should error from inside ncvreg,
  # which proves nfolds reached it.
  expect_error(scad_weights(X, y, nfolds = 1))
})

# ---- l0learn_weights ----
test_that("l0learn_weights computes weights with default L0 penalty", {
  skip_if_not_installed("L0Learn")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- l0learn_weights(X, y)
  expect_equal(nrow(result), p)
  expect_equal(ncol(result), 1)
  expect_true(all(is.finite(result)))
})

test_that("l0learn_weights computes weights with L0L2 penalty", {
  skip_if_not_installed("L0Learn")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- l0learn_weights(X, y, penalty = "L0L2")
  expect_equal(nrow(result), p)
  expect_equal(ncol(result), 1)
  expect_true(all(is.finite(result)))
})

test_that("l0learn_weights warns and zero-pads when X has zero-variance columns", {
  skip_if_not_installed("L0Learn")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  X[, 6] <- 4
  y <- X[, 1] * 0.5 + rnorm(n)
  expect_warning(result <- l0learn_weights(X, y),
                 "l0learn_weights: dropping 1 zero-variance column")
  expect_equal(nrow(result), p)
  expect_equal(result[6, 1], 0)
})

test_that("l0learn_weights passes nGamma through to L0Learn.cvfit", {
  skip_if_not_installed("L0Learn")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  # nGamma = 0 is invalid; the underlying function should reject it,
  # proving the argument reached L0Learn.cvfit.
  expect_error(l0learn_weights(X, y, penalty = "L0L2", nGamma = 0))
})

# ---- bglr_weights / bayes_b_weights / b_lasso_weights ----
test_that("bglr_weights computes weights with BayesB model", {
  skip_if_not_installed("BGLR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- bglr_weights(X, y, model = "BayesB", nIter = 100, burnIn = 20, thin = 2,
                         eta_args = list(probIn = 0.05))
  expect_equal(length(result), p)
  expect_true(is.numeric(result))
  expect_true(all(is.finite(result)))
})

test_that("bayes_b_weights computes weights and dispatches to bglr_weights", {
  skip_if_not_installed("BGLR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- bayes_b_weights(X, y, nIter = 100, burnIn = 20, thin = 2)
  expect_equal(length(result), p)
  expect_true(all(is.finite(result)))
})

test_that("bayes_b_weights passes probIn through to BGLR ETA", {
  skip_if_not_installed("BGLR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)

  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    BGLR = function(y, ETA, ...) {
      captured$eta <- ETA
      list(ETA = list(list(b = rep(0, ncol(ETA[[1]]$X)))))
    },
    .package = "BGLR"
  )

  result <- bayes_b_weights(X, y, nIter = 100, burnIn = 20, thin = 2, probIn = 0.42)
  expect_equal(length(result), p)
  expect_equal(captured$eta[[1]]$model, "BayesB")
  expect_equal(captured$eta[[1]]$probIn, 0.42)
})

test_that("b_lasso_weights computes weights with BL model", {
  skip_if_not_installed("BGLR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- b_lasso_weights(X, y, nIter = 100, burnIn = 20, thin = 2)
  expect_equal(length(result), p)
  expect_true(all(is.finite(result)))
})

test_that("bglr_weights warns and zero-pads when X has zero-variance columns", {
  skip_if_not_installed("BGLR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  X[, 7] <- 9
  y <- X[, 1] * 0.5 + rnorm(n)
  expect_warning(
    result <- bayes_b_weights(X, y, nIter = 100, burnIn = 20, thin = 2),
    "bglr_weights: dropping 1 zero-variance column"
  )
  expect_equal(length(result), p)
  expect_equal(result[7], 0)
})

test_that("bglr_weights cleans up its tempdir on exit", {
  skip_if_not_installed("BGLR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  before <- list.files(tempdir(), pattern = "^bglr_")
  bayes_b_weights(X, y, nIter = 100, burnIn = 20, thin = 2)
  after <- list.files(tempdir(), pattern = "^bglr_")
  expect_setequal(before, after)
})

# ---- dpr_weights ----
test_that("dpr_weights computes weights with VB fitting method", {
  skip_if_not_installed("RcppDPR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- dpr_weights(X, y, fitting_method = "VB")
  expect_equal(length(result), p)
  expect_true(is.numeric(result))
  expect_true(all(is.finite(result)))
})

test_that("dpr_weights computes weights with Gibbs fitting method", {
  skip_if_not_installed("RcppDPR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- dpr_weights(X, y, fitting_method = "Gibbs")
  expect_equal(length(result), p)
  expect_true(is.numeric(result))
  expect_true(all(is.finite(result)))
})

test_that("dpr_weights warns and zero-pads when X has zero-variance columns", {
  skip_if_not_installed("RcppDPR")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  X[, 8] <- 6
  y <- X[, 1] * 0.5 + rnorm(n)
  expect_warning(result <- dpr_weights(X, y),
                 "dpr_weights: dropping 1 zero-variance column")
  expect_equal(length(result), p)
  expect_equal(result[8], 0)
})

# ---- dispatch verification ----
# These tests use local_mocked_bindings to replace the inner function with a
# stub that captures its arguments, then assert the wrapper forwarded the
# correct values. They catch silent dispatch bugs (wrong method, wrong penalty,
# dropped argument) that the shape-only tests above would not.

test_that("prs_cs_weights dispatches to prs_cs with correct arguments", {
  set.seed(42)
  p <- 10
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  # Heterogeneous n with median != mean: this lets the test catch a regression
  # that swapped median(stat$n) for mean(stat$n) or sum(stat$n).
  # sorted: 10, 20, 30, 40, 50, 60, 70, 80, 90, 1000 -> median = 55, mean = 145.
  stat <- list(b = bhat, n = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 1000))
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    prs_cs = function(bhat, LD, n, ...) {
      captured$bhat <- bhat
      captured$LD <- LD
      captured$n <- n
      captured$dots <- list(...)
      list(beta_est = seq_len(length(bhat)) * 0.01)
    }
  )
  result <- prs_cs_weights(stat = stat, LD = R, maf = rep(0.3, p), n_iter = 17)
  expect_equal(captured$bhat, bhat)
  expect_equal(captured$LD, list(blk1 = R))
  expect_equal(captured$n, 55) # median of stat$n, NOT mean (145)
  expect_equal(captured$dots$maf, rep(0.3, p))
  expect_equal(captured$dots$n_iter, 17)
  expect_equal(result, seq_len(p) * 0.01)
})

test_that("sdpr_weights dispatches to sdpr with correct arguments", {
  set.seed(42)
  p <- 10
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  stat <- list(b = bhat, n = rep(456, p))
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    sdpr = function(bhat, LD, n, ...) {
      captured$bhat <- bhat
      captured$LD <- LD
      captured$n <- n
      captured$dots <- list(...)
      list(beta_est = seq_len(length(bhat)) * 0.02)
    }
  )
  result <- sdpr_weights(stat = stat, LD = R, iter = 19, burn = 3)
  expect_equal(captured$bhat, bhat)
  expect_equal(captured$LD, list(blk1 = R))
  expect_equal(captured$n, 456)
  expect_equal(captured$dots$iter, 19)
  expect_equal(captured$dots$burn, 3)
  expect_equal(result, seq_len(p) * 0.02)
})

test_that("lassosum_rss_weights dispatches to lassosum_rss once per s value", {
  set.seed(42)
  p <- 10
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  for (i in 1:(p - 1)) {
    R[i, i + 1] <- 0.4
    R[i + 1, i] <- 0.4
  }
  stat <- list(b = bhat, n = rep(100, p))
  call_log <- new.env(parent = emptyenv())
  call_log$calls <- list()
  local_mocked_bindings(
    lassosum_rss = function(bhat, LD, n, ...) {
      call_log$calls <- c(call_log$calls,
                          list(list(bhat = bhat, LD = LD, n = n)))
      list(beta_est = rep(0.05, length(bhat)), fbeta = c(1.0, 0.5))
    }
  )
  lassosum_rss_weights(stat = stat, LD = R, s = c(0.2, 0.9))
  expect_equal(length(call_log$calls), 2L)
  expect_equal(call_log$calls[[1]]$n, 100)
  expect_equal(length(call_log$calls[[1]]$bhat), p)
  # LD should differ between calls because s is different
  expect_false(identical(call_log$calls[[1]]$LD, call_log$calls[[2]]$LD))
})

test_that("bayes_{n,l,a,c,r}_weights each dispatch to bayes_alphabet_weights with correct method", {
  set.seed(42)
  X <- matrix(rnorm(50), nrow = 10)
  y <- rnorm(10)
  dispatchers <- list(
    list(fn = bayes_n_weights, expected = "bayesN"),
    list(fn = bayes_l_weights, expected = "bayesL"),
    list(fn = bayes_a_weights, expected = "bayesA"),
    list(fn = bayes_c_weights, expected = "bayesC"),
    list(fn = bayes_r_weights, expected = "bayesR")
  )
  for (d in dispatchers) {
    captured <- new.env(parent = emptyenv())
    local_mocked_bindings(
      bayes_alphabet_weights = function(X, y, method, ...) {
        captured$method <- method
        rep(0, ncol(X))
      }
    )
    d$fn(X, y)
    expect_equal(captured$method, d$expected,
                 label = paste("dispatch for", d$expected))
  }
})

test_that("scad_weights and mcp_weights dispatch to ncvreg_weights with correct penalty", {
  set.seed(42)
  X <- matrix(rnorm(50), nrow = 10)
  y <- rnorm(10)
  dispatchers <- list(
    list(fn = scad_weights, expected_penalty = "SCAD", nfolds = 7),
    list(fn = mcp_weights,  expected_penalty = "MCP",  nfolds = 9)
  )
  for (d in dispatchers) {
    captured <- new.env(parent = emptyenv())
    local_mocked_bindings(
      ncvreg_weights = function(X, y, penalty, nfolds = 5, ...) {
        captured$penalty <- penalty
        captured$nfolds <- nfolds
        matrix(0, nrow = ncol(X), ncol = 1)
      }
    )
    d$fn(X, y, nfolds = d$nfolds)
    expect_equal(captured$penalty, d$expected_penalty,
                 label = paste("penalty for", d$expected_penalty))
    expect_equal(captured$nfolds, d$nfolds,
                 label = paste("nfolds for", d$expected_penalty))
  }
})

test_that("b_lasso_weights dispatches to bglr_weights with model = 'BL'", {
  set.seed(42)
  X <- matrix(rnorm(50), nrow = 10)
  y <- rnorm(10)
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    bglr_weights = function(X, y, model, nIter, burnIn, thin, ...) {
      captured$model <- model
      captured$nIter <- nIter
      captured$burnIn <- burnIn
      captured$thin <- thin
      rep(0, ncol(X))
    }
  )
  b_lasso_weights(X, y, nIter = 77, burnIn = 11, thin = 3)
  expect_equal(captured$model, "BL")
  expect_equal(captured$nIter, 77)
  expect_equal(captured$burnIn, 11)
  expect_equal(captured$thin, 3)
})

test_that("bayes_b_weights dispatches to bglr_weights with model = 'BayesB' and probIn", {
  set.seed(42)
  X <- matrix(rnorm(50), nrow = 10)
  y <- rnorm(10)
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    bglr_weights = function(X, y, model, nIter, burnIn, thin, eta_args = list(), ...) {
      captured$model <- model
      captured$eta_args <- eta_args
      rep(0, ncol(X))
    }
  )
  bayes_b_weights(X, y, nIter = 100, burnIn = 20, thin = 2, probIn = 0.42)
  expect_equal(captured$model, "BayesB")
  expect_equal(captured$eta_args, list(probIn = 0.42))
})

test_that("lasso_weights and enet_weights dispatch to glmnet_weights with correct alpha", {
  set.seed(42)
  X <- matrix(rnorm(50), nrow = 10)
  y <- rnorm(10)
  dispatchers <- list(
    list(fn = lasso_weights, expected_alpha = 1),
    list(fn = enet_weights,  expected_alpha = 0.5)
  )
  for (d in dispatchers) {
    captured <- new.env(parent = emptyenv())
    local_mocked_bindings(
      glmnet_weights = function(X, y, alpha) {
        captured$alpha <- alpha
        matrix(0, nrow = ncol(X), ncol = 1)
      }
    )
    d$fn(X, y)
    expect_equal(captured$alpha, d$expected_alpha,
                 label = paste("alpha for", d$expected_alpha))
  }
})

test_that("susie_weights actually calls susie_wrapper when fit is NULL", {
  set.seed(42)
  p <- 5
  X <- matrix(rnorm(10 * p), nrow = 10)
  y <- rnorm(10)
  captured <- new.env(parent = emptyenv())
  captured$called <- FALSE
  local_mocked_bindings(
    susie_wrapper = function(X, y, ...) {
      captured$called <- TRUE
      captured$X <- X
      captured$y <- y
      list(pip = rep(0.1, ncol(X)))
    }
  )
  susie_weights(X = X, y = y)
  expect_true(captured$called)
  expect_identical(captured$X, X)
  expect_identical(captured$y, y)
})

test_that("susie_ash_weights calls susie_wrapper with ash dispatch arguments", {
  set.seed(42)
  p <- 5
  X <- matrix(rnorm(10 * p), nrow = 10)
  y <- rnorm(10)
  captured <- new.env(parent = emptyenv())
  captured$called <- FALSE
  local_mocked_bindings(
    susie_wrapper = function(X, y, unmappable_effects = NULL, convergence_method = NULL, ...) {
      captured$called <- TRUE
      captured$unmappable_effects <- unmappable_effects
      captured$convergence_method <- convergence_method
      list(pip = rep(0.1, ncol(X)))
    }
  )
  susie_ash_weights(X = X, y = y)
  expect_true(captured$called)
  expect_equal(captured$unmappable_effects, "ash")
  expect_equal(captured$convergence_method, "pip")
})

test_that("susie_inf_weights calls susie_wrapper with inf dispatch arguments", {
  set.seed(42)
  p <- 5
  X <- matrix(rnorm(10 * p), nrow = 10)
  y <- rnorm(10)
  captured <- new.env(parent = emptyenv())
  captured$called <- FALSE
  local_mocked_bindings(
    susie_wrapper = function(X, y, unmappable_effects = NULL, convergence_method = NULL, ...) {
      captured$called <- TRUE
      captured$unmappable_effects <- unmappable_effects
      captured$convergence_method <- convergence_method
      list(pip = rep(0.1, ncol(X)))
    }
  )
  susie_inf_weights(X = X, y = y)
  expect_true(captured$called)
  expect_equal(captured$unmappable_effects, "inf")
  expect_equal(captured$convergence_method, "pip")
})

test_that("mrash_weights actually calls lasso_weights for default beta.init", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  called <- new.env(parent = emptyenv())
  called$lasso <- FALSE
  local_mocked_bindings(
    lasso_weights = function(X, y) {
      called$lasso <- TRUE
      rep(0.01, ncol(X))
    },
    init_prior_sd = function(X, y, n = 30) seq(0, 3, length.out = n)
  )
  mrash_weights(X, y)
  expect_true(called$lasso)
})

test_that("mrash_weights calls init_prior_sd only when init_prior_sd = TRUE", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  called <- new.env(parent = emptyenv())

  # init_prior_sd = TRUE: init_prior_sd should be called
  called$init <- FALSE
  local_mocked_bindings(
    lasso_weights = function(X, y) rep(0.01, ncol(X)),
    init_prior_sd = function(X, y, n = 30) {
      called$init <- TRUE
      seq(0, 3, length.out = n)
    }
  )
  mrash_weights(X, y, init_prior_sd = TRUE)
  expect_true(called$init)

  # init_prior_sd = FALSE: init_prior_sd should NOT be called
  called$init <- FALSE
  mrash_weights(X, y, init_prior_sd = FALSE)
  expect_false(called$init)
})

# ---- end-to-end solver verification ----
# These tests mock the *external-package* solver and assert that the
# user-facing wrapper forwards the correct dispatch parameter all the way
# through the helper hop. This is the same pattern as the BGLR test at L1100,
# applied to every method. Together with the dispatch verification tests
# above, this gives full coverage of the parameter-passing chain
# wrapper -> helper -> solver.
#
# Two patterns are used:
#  - For solvers whose results are consumed via direct field access (BGLR,
#    qgg::gbayes, RcppDPR::fit_model), the mock returns a minimal fake list.
#  - For solvers whose results flow through S3 coef() dispatch (glmnet,
#    ncvreg, L0Learn), faking the result is awkward, so the mock captures
#    the parameter and then raises a sentinel error that the test catches.

test_that("glmnet_weights forwards alpha to glmnet::cv.glmnet", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  X <- matrix(rnorm(50), nrow = 10)
  y <- rnorm(10)
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    cv.glmnet = function(x, y, alpha, ...) {
      captured$alpha <- alpha
      stop("STOP_AFTER_CAPTURE")
    },
    .package = "glmnet"
  )
  expect_error(lasso_weights(X, y), "STOP_AFTER_CAPTURE")
  expect_equal(captured$alpha, 1)

  expect_error(enet_weights(X, y), "STOP_AFTER_CAPTURE")
  expect_equal(captured$alpha, 0.5)
})

test_that("ncvreg_weights forwards penalty to ncvreg::cv.ncvreg", {
  skip_if_not_installed("ncvreg")
  set.seed(42)
  X <- matrix(rnorm(50), nrow = 10)
  y <- rnorm(10)
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    cv.ncvreg = function(X, y, penalty, nfolds = 5, ...) {
      captured$penalty <- penalty
      captured$nfolds <- nfolds
      stop("STOP_AFTER_CAPTURE")
    },
    .package = "ncvreg"
  )
  expect_error(scad_weights(X, y, nfolds = 7), "STOP_AFTER_CAPTURE")
  expect_equal(captured$penalty, "SCAD")
  expect_equal(captured$nfolds, 7)

  expect_error(mcp_weights(X, y, nfolds = 9), "STOP_AFTER_CAPTURE")
  expect_equal(captured$penalty, "MCP")
  expect_equal(captured$nfolds, 9)
})

test_that("l0learn_weights forwards penalty to L0Learn::L0Learn.cvfit", {
  skip_if_not_installed("L0Learn")
  set.seed(42)
  X <- matrix(rnorm(50), nrow = 10)
  y <- rnorm(10)
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    L0Learn.cvfit = function(x, y, penalty, nFolds = 5, ...) {
      captured$penalty <- penalty
      captured$nFolds <- nFolds
      stop("STOP_AFTER_CAPTURE")
    },
    .package = "L0Learn"
  )
  expect_error(l0learn_weights(X, y), "STOP_AFTER_CAPTURE")
  expect_equal(captured$penalty, "L0")

  expect_error(l0learn_weights(X, y, penalty = "L0L2", nFolds = 8),
               "STOP_AFTER_CAPTURE")
  expect_equal(captured$penalty, "L0L2")
  expect_equal(captured$nFolds, 8)
})

test_that("b_lasso_weights forwards model = 'BL' all the way to BGLR::BGLR", {
  skip_if_not_installed("BGLR")
  set.seed(42)
  X <- matrix(rnorm(50), nrow = 10)
  y <- rnorm(10)
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    BGLR = function(y, ETA, nIter, burnIn, thin, ...) {
      captured$model <- ETA[[1]]$model
      captured$nIter <- nIter
      captured$burnIn <- burnIn
      captured$thin <- thin
      list(ETA = list(list(b = rep(0, ncol(ETA[[1]]$X)))))
    },
    .package = "BGLR"
  )
  b_lasso_weights(X, y, nIter = 77, burnIn = 11, thin = 3)
  expect_equal(captured$model, "BL")
  expect_equal(captured$nIter, 77)
  expect_equal(captured$burnIn, 11)
  expect_equal(captured$thin, 3)
})

test_that("bayes_alphabet_weights forwards method to qgg::gbayes for all alphabet variants", {
  skip_if_not_installed("qgg")
  set.seed(42)
  X <- matrix(rnorm(50), nrow = 10)
  y <- rnorm(10)
  Z <- matrix(rnorm(20), nrow = 10)
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    gbayes = function(y, W, X = NULL, method, nit, nburn, ...) {
      captured$method <- method
      captured$nit <- nit
      captured$nburn <- nburn
      captured$X <- X
      list(bm = rep(0, ncol(W)))
    },
    .package = "qgg"
  )
  for (m in c("bayesN", "bayesL", "bayesA", "bayesC", "bayesR")) {
    bayes_alphabet_weights(X, y, method = m, Z = Z, nit = 17, nburn = 4)
    expect_equal(captured$method, m, info = paste("method =", m))
    expect_equal(captured$nit, 17)
    expect_equal(captured$nburn, 4)
    # Z is forwarded to qgg::gbayes via the X argument.
    expect_equal(captured$X, Z, info = paste("Z forwarding for method =", m))
  }
})

test_that("dpr_weights forwards fitting_method to RcppDPR::fit_model", {
  skip_if_not_installed("RcppDPR")
  set.seed(42)
  X <- matrix(rnorm(50), nrow = 10)
  y <- rnorm(10)
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    fit_model = function(y, w, x, rotate_variables, fitting_method, ...) {
      captured$fitting_method <- fitting_method
      captured$rotate_variables <- rotate_variables
      list(beta = rep(0, ncol(x)), alpha = rep(0, ncol(x)))
    },
    .package = "RcppDPR"
  )
  dpr_weights(X, y, fitting_method = "VB")
  expect_equal(captured$fitting_method, "VB")
  expect_false(captured$rotate_variables)

  dpr_weights(X, y, fitting_method = "Gibbs")
  expect_equal(captured$fitting_method, "Gibbs")
})

# ============================================================================
# Signal recovery tests for X-y wrappers
#
# These tests catch breaking regressions in the post-solver coefficient
# extraction code (which the dispatch/end-to-end mocks short-circuit). They
# simulate a small sparse linear model and assert that the returned weights
# correlate with truth (or that the top-K nonzero indices overlap truth).
# ============================================================================

# Helper: simulate a small sparse linear model with binomial genotypes.
.simulate_sparse_xy <- function(n = 200, p = 20,
                                signal_idx = c(3, 10, 15),
                                signal_beta = c(0.8, -0.6, 0.5),
                                seed = 2024) {
  set.seed(seed)
  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n)
  beta_true <- rep(0, p)
  beta_true[signal_idx] <- signal_beta
  y <- as.numeric(X %*% beta_true) + rnorm(n, sd = 0.5)
  list(X = X, y = y, beta_true = beta_true, signal_idx = signal_idx)
}

test_that("lasso_weights recovers signal direction on simulated data", {
  skip_if_not_installed("glmnet")
  sim <- .simulate_sparse_xy()
  w <- lasso_weights(sim$X, sim$y)
  expect_equal(length(w), ncol(sim$X))
  expect_true(all(is.finite(w)))
  expect_gt(cor(w, sim$beta_true), 0.5)
})

test_that("scad_weights recovers signal indices on simulated data", {
  skip_if_not_installed("ncvreg")
  sim <- .simulate_sparse_xy()
  w <- scad_weights(sim$X, sim$y, nfolds = 5)
  expect_equal(length(w), ncol(sim$X))
  expect_true(all(is.finite(w)))
  # SCAD is sparse: top-3 by absolute weight should heavily overlap truth.
  top3 <- order(abs(w), decreasing = TRUE)[1:3]
  expect_gte(length(intersect(top3, sim$signal_idx)), 2L)
})

test_that("l0learn_weights with default L0 penalty recovers signal indices", {
  skip_if_not_installed("L0Learn")
  sim <- .simulate_sparse_xy()
  w <- l0learn_weights(sim$X, sim$y, penalty = "L0", nFolds = 5)
  expect_equal(length(w), ncol(sim$X))
  expect_true(all(is.finite(w)))
  top3 <- order(abs(w), decreasing = TRUE)[1:3]
  expect_gte(length(intersect(top3, sim$signal_idx)), 2L)
})

test_that("l0learn_weights with L0L2 penalty recovers signal direction", {
  skip_if_not_installed("L0Learn")
  sim <- .simulate_sparse_xy()
  w <- l0learn_weights(sim$X, sim$y, penalty = "L0L2", nFolds = 5)
  expect_equal(length(w), ncol(sim$X))
  expect_true(all(is.finite(w)))
  # L0L2 selects (gamma, lambda) over a 2D grid; correlation is the safer bet.
  expect_gt(cor(w, sim$beta_true), 0.4)
})

test_that("dpr_weights with VB fitting recovers signal direction", {
  skip_if_not_installed("RcppDPR")
  sim <- .simulate_sparse_xy()
  w <- dpr_weights(sim$X, sim$y, fitting_method = "VB")
  expect_equal(length(w), ncol(sim$X))
  expect_true(all(is.finite(w)))
  expect_gt(cor(w, sim$beta_true), 0.4)
})

test_that("lassosum_rss_weights recovers signal direction on simulated data", {
  skip_if_not_installed("Rcpp")
  set.seed(2024)
  n <- 500
  p <- 20
  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n)
  beta_true <- rep(0, p)
  beta_true[c(3, 10, 15)] <- c(0.4, -0.3, 0.2)
  y <- as.numeric(X %*% beta_true) + rnorm(n)
  bhat <- as.vector(cor(y, X))
  R <- cor(X)
  # Heterogeneous n forces the median(stat$n) path; mean would be wrong.
  stat <- list(b = bhat, n = c(rep(n, p - 1), 10 * n))
  w <- lassosum_rss_weights(stat = stat, LD = R, s = c(0.2, 0.5, 0.9))
  expect_equal(length(w), p)
  expect_true(all(is.finite(w)))
  expect_gt(cor(w, beta_true), 0.5)
})

# ============================================================================
# mr_ash_rss_weights — dispatch + smoke test
# ============================================================================

test_that("mr_ash_rss_weights forwards arguments to susieR::mr.ash.rss", {
  set.seed(42)
  p <- 10
  bhat <- rnorm(p, sd = 0.1)
  shat <- rep(0.05, p)
  R <- diag(p)
  # Heterogeneous n with median != mean: median = 55, mean = 145.
  stat <- list(
    b = bhat,
    seb = shat,
    n = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 1000)
  )
  z_vec <- bhat / shat
  captured <- new.env(parent = emptyenv())
  # mr.ash.rss is `@importFrom`'d into pecotmr's namespace, so mock the
  # pecotmr binding (no .package argument), not the susieR binding.
  local_mocked_bindings(
    mr.ash.rss = function(bhat, shat, z, R, var_y, n, sigma2_e, s0, w0, ...) {
      captured$bhat <- bhat
      captured$shat <- shat
      captured$z <- z
      captured$R <- R
      captured$var_y <- var_y
      captured$n <- n
      captured$sigma2_e <- sigma2_e
      captured$s0 <- s0
      captured$w0 <- w0
      captured$dots <- list(...)
      list(mu1 = seq_len(length(bhat)) * 0.01)
    }
  )
  result <- mr_ash_rss_weights(
    stat = stat, LD = R, var_y = 1.5, sigma2_e = 0.4,
    s0 = c(0, 0.1, 0.2), w0 = c(0.5, 0.3, 0.2), z = z_vec,
    tol = 1e-6
  )
  expect_equal(captured$bhat, bhat)
  expect_equal(captured$shat, shat)
  expect_equal(captured$z, z_vec)
  expect_equal(captured$R, R)
  expect_equal(captured$var_y, 1.5)
  expect_equal(captured$n, 55) # median of stat$n, NOT mean (145)
  expect_equal(captured$sigma2_e, 0.4)
  expect_equal(captured$s0, c(0, 0.1, 0.2))
  expect_equal(captured$w0, c(0.5, 0.3, 0.2))
  expect_equal(captured$dots$tol, 1e-6)
  expect_equal(result, seq_len(p) * 0.01)
})

test_that("mr_ash_rss_weights returns a numeric vector of expected length", {
  skip_if_not_installed("susieR")
  set.seed(2024)
  n <- 500
  p <- 10
  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n)
  beta_true <- rep(0, p)
  beta_true[c(3, 7)] <- c(0.4, -0.3)
  y <- as.numeric(X %*% beta_true) + rnorm(n)
  bhat <- as.vector(cor(y, X))
  shat <- rep(1 / sqrt(n - 1), p)
  R <- cor(X)
  stat <- list(b = bhat, seb = shat, n = rep(n, p))
  w <- mr_ash_rss_weights(
    stat = stat, LD = R, var_y = var(y),
    sigma2_e = NULL,
    s0 = c(0, 0.01, 0.1, 0.5),
    w0 = rep(1 / 4, 4)
  )
  expect_equal(length(w), p)
  expect_true(all(is.finite(w)))
})
