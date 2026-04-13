context("regularized_regression — dispatch verification")

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

gc()
