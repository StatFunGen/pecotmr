context("Molecular trait imputation pipeline")

# Helper: low-rank synthetic matrix with NA mask.
.make_imputation_inputs <- function(n = 60, p = 15, rank = 3, na_frac = 0.1,
                                     seed = 1) {
  set.seed(seed)
  U <- matrix(rnorm(n * rank), n, rank)
  V <- matrix(rnorm(p * rank), p, rank)
  truth <- U %*% t(V) + matrix(rnorm(n * p, sd = 0.1), n, p)
  X <- truth
  na_mask <- matrix(runif(n * p) < na_frac, n, p)
  X[na_mask] <- NA
  rownames(X) <- paste0("S", seq_len(n))
  colnames(X) <- paste0("F", seq_len(p))
  list(X = X, truth = truth, na_mask = na_mask)
}

# =============================================================================
# softimpute_imputation
# =============================================================================

test_that("softimpute_imputation: returns matrix with same dims, no NAs", {
  skip_if_not_installed("softImpute")
  inp <- .make_imputation_inputs()
  Xi <- softimpute_imputation(inp$X)
  expect_equal(dim(Xi), dim(inp$X))
  expect_false(any(is.na(Xi)))
  expect_equal(rownames(Xi), rownames(inp$X))
})

test_that("softimpute_imputation: handles all-NA columns by removing them", {
  skip_if_not_installed("softImpute")
  inp <- .make_imputation_inputs()
  X <- inp$X
  X[, 5] <- NA
  expect_message(
    Xi <- softimpute_imputation(X, verbose = TRUE),
    "Removed 1 column"
  )
  expect_equal(ncol(Xi), ncol(X) - 1L)
})

test_that("softimpute_imputation: default tunes lambda via CV and attaches grid", {
  skip_if_not_installed("softImpute")
  inp <- .make_imputation_inputs()
  Xi <- softimpute_imputation(inp$X)  # default lambda = NULL -> tune
  expect_false(is.null(attr(Xi, "softimpute_lambda")))
  grid <- attr(Xi, "softimpute_lambda_cv")
  expect_s3_class(grid, "data.frame")
  expect_named(grid, c("lambda", "rmse"))
  expect_true(all(diff(grid$lambda) <= 0))  # grid sorted descending
  expect_true(any(is.finite(grid$rmse)))
  # The chosen lambda must be one of the grid points
  expect_true(attr(Xi, "softimpute_lambda") %in% grid$lambda)
})

test_that("softimpute_imputation: scalar lambda bypasses tuning", {
  skip_if_not_installed("softImpute")
  inp <- .make_imputation_inputs()
  l0 <- softImpute::lambda0(inp$X)
  Xi <- softimpute_imputation(inp$X, lambda = 0.3 * l0)
  expect_equal(attr(Xi, "softimpute_lambda"), 0.3 * l0)
  expect_null(attr(Xi, "softimpute_lambda_cv"))
})

test_that("softimpute_imputation: vector lambda restricts tuning to that grid", {
  skip_if_not_installed("softImpute")
  inp <- .make_imputation_inputs()
  l0 <- softImpute::lambda0(inp$X)
  custom_grid <- l0 * c(0.05, 0.1, 0.25, 0.5)
  Xi <- softimpute_imputation(inp$X, lambda = custom_grid)
  grid <- attr(Xi, "softimpute_lambda_cv")
  expect_setequal(grid$lambda, custom_grid)
  expect_true(attr(Xi, "softimpute_lambda") %in% custom_grid)
})

test_that("softimpute_imputation: chosen lambda minimises CV RMSE on the grid", {
  skip_if_not_installed("softImpute")
  inp <- .make_imputation_inputs(n = 80, p = 20, rank = 3, seed = 11)
  Xi <- softimpute_imputation(inp$X)
  grid <- attr(Xi, "softimpute_lambda_cv")
  chosen <- attr(Xi, "softimpute_lambda")
  ok <- which(is.finite(grid$rmse))
  best_idx <- ok[which.min(grid$rmse[ok])]
  expect_equal(chosen, grid$lambda[best_idx])
})

# =============================================================================
# flashier_imputation
# =============================================================================

test_that("flashier_imputation: returns matrix with same dims, no NAs", {
  skip_if_not_installed("flashier")
  skip_if_not_installed("ebnm")
  inp <- .make_imputation_inputs()
  Xi <- flashier_imputation(inp$X, verbose = 0)
  expect_equal(dim(Xi), dim(inp$X))
  expect_false(any(is.na(Xi)))
})

test_that("flashier_imputation: recovers low-rank truth better than column means", {
  skip_if_not_installed("flashier")
  skip_if_not_installed("ebnm")
  inp <- .make_imputation_inputs()
  Xi <- flashier_imputation(inp$X, verbose = 0)
  # Baseline: column-mean imputation
  Xbase <- inp$X
  cm <- colMeans(Xbase, na.rm = TRUE)
  for (j in seq_len(ncol(Xbase))) Xbase[is.na(Xbase[, j]), j] <- cm[j]
  rmse_flash <- sqrt(mean((Xi[inp$na_mask] - inp$truth[inp$na_mask])^2))
  rmse_mean  <- sqrt(mean((Xbase[inp$na_mask] - inp$truth[inp$na_mask])^2))
  expect_lt(rmse_flash, rmse_mean)
})

# =============================================================================
# molecular_trait_imputation_pipeline
# =============================================================================

test_that("pipeline: single-method dispatch returns named list keyed by method", {
  skip_if_not_installed("softImpute")
  inp <- .make_imputation_inputs()
  res <- molecular_trait_imputation_pipeline(inp$X, methods = "softimpute")
  expect_named(res$imputed, "softimpute")
  expect_equal(dim(res$imputed[["softimpute"]]), dim(inp$X))
  expect_null(res$evaluation)
})

test_that("pipeline: multi-method runs all three and returns one matrix per method", {
  skip_if_not_installed("softImpute")
  skip_if_not_installed("flashier")
  skip_if_not_installed("xgboost")
  inp <- .make_imputation_inputs()
  res <- molecular_trait_imputation_pipeline(
    inp$X, methods = c("xgboost", "softimpute", "flashier"),
    method_args = list(
      xgboost = list(maxiter = 2, nrounds = 10, verbose = FALSE)
    )
  )
  expect_named(res$imputed, c("xgboost", "softimpute", "flashier"))
  for (m in names(res$imputed)) {
    expect_equal(dim(res$imputed[[m]]), dim(inp$X))
    expect_false(any(is.na(res$imputed[[m]])))
  }
})

test_that("pipeline: evaluate = TRUE reports RMSE on held-out cells", {
  skip_if_not_installed("softImpute")
  skip_if_not_installed("flashier")
  inp <- .make_imputation_inputs()
  res <- molecular_trait_imputation_pipeline(
    inp$X, methods = c("softimpute", "flashier"),
    evaluate = TRUE, eval_fraction = 0.1, eval_seed = 7
  )
  expect_s3_class(res$evaluation, "data.frame")
  expect_equal(nrow(res$evaluation), 2L)
  expect_named(res$evaluation,
               c("method", "rmse", "n_held_out", "time_seconds"))
  expect_true(all(res$evaluation$rmse > 0))
  expect_true(all(res$evaluation$n_held_out > 0))
})

test_that("pipeline: missingness summary matches input NA count", {
  skip_if_not_installed("softImpute")
  inp <- .make_imputation_inputs()
  res <- molecular_trait_imputation_pipeline(inp$X, methods = "softimpute")
  expect_equal(res$missingness$n_missing, sum(is.na(inp$X)))
  expect_equal(res$missingness$fraction,
               sum(is.na(inp$X)) / prod(dim(inp$X)))
})

test_that("pipeline: rejects unknown method", {
  inp <- .make_imputation_inputs()
  expect_error(
    molecular_trait_imputation_pipeline(inp$X, methods = "knnimpute"),
    "Unknown imputation method"
  )
})

test_that("pipeline: rejects empty methods vector", {
  inp <- .make_imputation_inputs()
  expect_error(
    molecular_trait_imputation_pipeline(inp$X, methods = character(0)),
    "non-empty character vector"
  )
})
