context("regularized_regression - lassosumRss")

# ---- lassosumRss ----
test_that("lassosumRss errors on invalid LD input", {
  expect_error(lassosumRss(bhat = rnorm(5), LD = "not_a_list", n = 100),
               "valid list of LD blocks")
})

test_that("lassosumRss errors on non-positive sample size", {
  expect_error(lassosumRss(bhat = rnorm(5), LD = list(blk1 = diag(5)), n = -1),
               "valid sample size")
})

test_that("lassosumRss errors on mismatched bhat and LD dimensions", {
  expect_error(
    lassosumRss(bhat = rnorm(10), LD = list(blk1 = diag(5)), n = 100),
    "same as the sum"
  )
})

test_that("lassosumRss runs successfully with valid input", {
  set.seed(42)
  p <- 10
  n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  for (i in 1:(p - 1)) {
    R[i, i + 1] <- 0.3
    R[i + 1, i] <- 0.3
  }
  result <- lassosumRss(bhat = bhat, LD = list(blk1 = R), n = n)
  expect_type(result, "list")
  expect_true("betaEst" %in% names(result))
  expect_equal(length(result$betaEst), p)
  expect_true(all(is.finite(result$betaEst)))
  expect_equal(nrow(result$beta), p)
  expect_equal(ncol(result$beta), 20)
})

test_that("lassosumRss accepts multiple LD blocks", {
  set.seed(42)
  p1 <- 5
  p2 <- 5
  p <- p1 + p2
  n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R1 <- diag(p1)
  R2 <- diag(p2)
  result <- lassosumRss(bhat = bhat, LD = list(blk1 = R1, blk2 = R2), n = n)
  expect_type(result, "list")
  expect_true("betaEst" %in% names(result))
  expect_equal(length(result$betaEst), p)
  expect_true(all(is.finite(result$betaEst)))
})

test_that("lassosumRss with large lambda gives all-zero betas", {
  set.seed(42)
  p <- 10
  n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  result <- lassosumRss(bhat = bhat, LD = list(blk1 = R), n = n,
                         lambda = c(100))
  expect_true(all(result$betaEst == 0))
})

# ---- lassosumRssWeights (wrapper) ----
test_that("lassosumRssWeights calls lassosumRss and returns betaEst", {
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
  expected <- seq_len(p) * 0.02
  local_mocked_bindings(
    lassosumRss = function(bhat, LD, n, ...) {
      list(
        beta = cbind(rep(0, length(bhat)), expected),
        lambda = c(0.05, 0.01),
        fbeta = c(1, 0.1)
      )
    }
  )
  result <- lassosumRssWeights(stat = stat, LD = R, s = 0.5)
  expect_equal(length(result), p)
  expect_true(is.numeric(result))
  expect_equal(c(result), expected)
  expect_equal(unname(attr(result, "lassosum_selection")["mode"]), "ld_quadratic")
})

test_that("lassosumRssWeights clamps correlation input before scaling", {
  set.seed(42)
  p <- 5
  n <- 100
  # Construct bhat with max abs >= 1 to trigger the clamp branch
  bhat <- c(1.5, -0.2, 0.1, 0.05, -0.3)
  R <- diag(p)
  stat <- list(b = bhat, n = rep(n, p))
  captured <- NULL
  local_mocked_bindings(
    lassosumRss = function(bhat, LD, n, ...) {
      captured <<- bhat
      list(
        beta = matrix(0, nrow = length(bhat), ncol = 2),
        lambda = c(0.05, 0.01),
        fbeta = c(1, 0.5)
      )
    }
  )
  result <- lassosumRssWeights(stat = stat, LD = R, s = 0.5, selection = "min_fbeta")
  scaled_cor <- captured / sqrt(n)
  expect_true(max(abs(scaled_cor)) < 1)
  expect_equal(max(abs(scaled_cor)), 0.9999, tolerance = 1e-9)
  # Sign-preserving rescale
  expect_true(all(sign(scaled_cor) == sign(bhat)))
})
