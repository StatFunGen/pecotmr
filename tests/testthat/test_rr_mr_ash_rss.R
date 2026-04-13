context("regularized_regression — mr_ash_rss")

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
