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
