context("regularized_regression — end-to-end solver verification")

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
