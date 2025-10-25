context("twas")
library(tidyverse)
library(snpStats)

generate_X_Y <- function(seed=1, num_samples=10, num_features=10, X_rownames=TRUE, y_rownames=TRUE) {
  set.seed(seed)
  X <- scale(
    matrix(rnorm(num_samples * num_features), nrow = num_samples),
    center = TRUE, scale = TRUE)
  
  if (X_rownames) {
    rownames(X) <- paste0("sample", 1:num_samples)
  } else {
    rownames(X) <- NULL
  }
  
  beta = rep(0, num_features)
  beta[1:4] = 1
  y <- X %*% beta + rnorm(num_samples)
  y <- matrix(y, nrow = num_samples, ncol = 1)
  if (y_rownames) {
    rownames(y) <- paste0("sample", 1:num_samples)
  } else {
    rownames(y) <- NULL
  }
  colnames(y) <- c("Outcome")
  
  return(list(X=X, Y=y))
}

generate_susie_obj <- function(X, y) {
    return(susie(X,y, L=10,
        max_iter=500,
        estimate_residual_variance=TRUE,
        estimate_prior_variance=TRUE,
        refine=TRUE,
        compute_univariate_zscore=FALSE,
        min_abs_corr=0.5,
        coverage=0.95))
}

test_that("Check twas_z works", {
    X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3)
    weights <- c(0.5, 0.5)
    z <- c(2, 3)
    R <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
    result <- twas_z(weights, z, R)
    expect_true(is.list(result))
    expect_true(all(names(result) %in% c("z", "pval")))
    # Add checks for expected values of z and pval
})

test_that("Check twas_z weights and z-scores length match", {
    expect_error(twas_z(c(1, 2), c(1, 2, 3)))
})

test_that("Check twas_z handle NA in X", {
    X <- matrix(c(1, NA, 3, 4, 5, 6), nrow = 3)
    weights <- c(0.5, 0.5)
    z <- c(2, 3)
    result <- twas_z(weights, z, X = X)
    expect_true(is.list(result))
})

test_that("Check twas_z R is NULL and X is provided", {
    X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3)
    weights <- c(0.5, 0.5)
    z <- c(2, 3)
    result <- twas_z(weights, z, X = X)
    expect_true(is.list(result))
})

# Use fold with made up partition dataframe

test_that("Check twas_weights_cv works with minimum data", {
    sim <- generate_X_Y()
    X <- sim$X
    y = sim$Y
    local_mocked_bindings(
        susie_weights = function(X, y, ...) rnorm(ncol(X)),
        glmnet_weights = function(X, y, ...) runif(ncol(X))
    )
    weight_methods_test <- list(susie_weights = c("init_prior_sd"), glmnet_weights = c("init_prior_sd"))
    result_min <- twas_weights_cv(X, y, fold = 2, weight_methods = weight_methods_test)
    expect_is(result_min, "list")
})

test_that("Check twas_weights_cv works with no weight methods", {
    sim <- generate_X_Y()
    X <- sim$X
    y = sim$Y
    result_no_methods <- twas_weights_cv(X, y, fold = 5)
    expect_is(result_no_methods, "list")
    expect_false("prediction" %in% names(result_no_methods))
})

test_that("twas_weights_cv is reproducible with seed", {
    sim <- generate_X_Y(seed=1)
    X <- sim$X
    y = sim$Y
    local_mocked_bindings(
        susie_weights = function(X, y, ...) rnorm(ncol(X)),
        glmnet_weights = function(X, y, ...) runif(ncol(X))
    )
    weight_methods_test <- list(susie_weights = c("init_prior_sd"), glmnet_weights = c("init_prior_sd"))
    set.seed(1)
    result_seed1 <- twas_weights_cv(X, y, fold = 2, weight_methods = weight_methods_test)
    set.seed(1)
    result_seed2 <- twas_weights_cv(X, y, fold = 2, weight_methods = weight_methods_test)
    expect_equal(result_seed1$sample_partition, result_seed2$sample_partition)
})

test_that("twas_weights_cv handles errors appropriately", {
    sim <- generate_X_Y(seed=1)
    X <- sim$X
    y = sim$Y
    local_mocked_bindings(
        susie_weights = function(X, y, ...) rnorm(ncol(X)),
        glmnet_weights = function(X, y, ...) runif(ncol(X))
    )
    weight_methods_test <- list(susie_weights = c("init_prior_sd"), glmnet_weights = c("init_prior_sd"))
    expect_error(twas_weights_cv(X, y, fold = NULL))
    expect_error(twas_weights_cv(X, y, fold = "invalid"))
    expect_error(twas_weights_cv(X, y, fold = -1))
    expect_error(twas_weights_cv(2, y, fold = 2))
    expect_error(twas_weights_cv(X, 2, fold = 2))
    expect_error(twas_weights_cv(matrix(rnorm(4, nrow=2)), matrix(rnorm(2, nrow=1)), fold = 2))
    expect_error(twas_weights_cv(X, y))
    #expect_error(twas_weights_cv(X, y, sample_partitions = data.frame(Sample = c("sample1", "sample2", "sample3"), Fold = c(1, 2, 3))))
})

# test_that("twas_weights_cv handles parallel processing", {
#     RNGkind("L'Ecuyer-CMRG")
#     sim <- generate_X_Y(seed=1, num_samples=30)
#     X <- sim$X
#     y = sim$Y
#     weight_methods_test <- list(
#         glmnet_weights = list(alpha = 0.5))
#     set.seed(1)
#     result_parallel <- twas_weights_cv(X, y, fold = 2, weight_methods = weight_methods_test, num_threads = 2)
#     set.seed(1)
#     result_single <- twas_weights_cv(X, y, fold = 2, weight_methods = weight_methods_test, num_threads = 1)
#     expect_is(result_parallel, "list")
#     expect_is(result_single, "list")
#     expect_equal(result_parallel$sample_partition, result_single$sample_partition)
#     expect_equal(result_parallel$prediction$glmnet_predicted, result_single$prediction$glmnet_predicted)
#     RNGkind("default")
# })
#
test_that("Check twas_weights works with minimum data", {
    sim <- generate_X_Y(seed=1)
    X <- sim$X
    y = sim$Y
    local_mocked_bindings(
        susie_weights = function(X, y, ...) rnorm(ncol(X)),
        glmnet_weights = function(X, y, ...) runif(ncol(X))
    )
    weight_methods_test <- list(susie_weights = c("init_prior_sd"), glmnet_weights = c("init_prior_sd"))
    result_min <- twas_weights(X, y, weight_methods = weight_methods_test)
    expect_is(result_min, "list")
})

test_that("twas_weights handles errors appropriately", {
    sim <- generate_X_Y(seed=1)
    X <- sim$X
    y = sim$Y
    local_mocked_bindings(
        susie_weights = function(X, y, ...) rnorm(ncol(X)),
        glmnet_weights = function(X, y, ...) runif(ncol(X))
    )
    weight_methods_test <- list(susie_weights = c("init_prior_sd"), glmnet_weights = c("init_prior_sd"))
    expect_error(twas_weights(matrix(rnorm(4, nrow=2)), matrix(rnorm(2, nrow=1))))
    expect_error(twas_weights(X, y))
})

# test_that("twas_weights handles parallel processing", {
#     RNGkind("L'Ecuyer-CMRG")
#     sim <- generate_X_Y(seed=1, num_samples=30)
#     X <- sim$X
#     y = sim$Y
#     weight_methods_test <- list(
#         glmnet_weights = list(alpha = 0.5))
#     set.seed(1)
#     result_parallel <- twas_weights(X, y, weight_methods = weight_methods_test, num_threads = 2)
#     set.seed(1)
#     result_single <- twas_weights(X, y, weight_methods = weight_methods_test, num_threads = 1)
#     expect_equal(result_parallel, result_single)
#     RNGkind("default")
# })
#
test_that("Check pval_acat works", {
    set.seed(1)
    expect_equal(pval_acat(c(0.05)), 0.05)
    expect_true(!is.null(pval_acat(runif(10))))
})

test_that("Check pval_hmp works", {
    set.seed(1)
    expect_true(!is.null(pval_hmp(runif(10))))
})

test_that("Check pval_global works", {
    set.seed(1)
    pvals <- runif(10)
    expect_true(!is.null(pval_global(runif(10)))) # naive returned
    expect_true(!is.null(pval_global(runif(10), naive=TRUE))) # acat returned
    expect_true(!is.null(pval_global(runif(10), naive=TRUE))) # hmp returned
})

generate_twas_joint_z_data <- function(num_samples=10, num_snps=10, num_conditions = 5) {
  X <- matrix(sample(0:2, num_samples * num_snps, replace = TRUE), nrow = num_snps, ncol = num_samples)
  rownames(X) <- paste0("Sample", 1:num_samples)
  colnames(X) <- paste0("SNP", 1:num_snps)

  weights <- matrix(rnorm(num_snps * num_conditions), nrow = num_snps, ncol = num_conditions)
  rownames(weights) <- paste0("SNP", 1:num_snps)
  colnames(weights) <- paste0("Cond", 1:num_conditions)

  z <- rnorm(num_snps)
  names(z) <- paste0("SNP", 1:num_snps)

  R <- cor(X)
  rownames(R) <- colnames(X) <- paste0("SNP", 1:ncol(X))

  return(list(X=X, weights=weights, z=z, R=R))
}

test_that("twas_joint_z fails when z noteq weights", {
    expect_error(twas_joint_z(weights = matrix(rnorm(10), nrow=10), z = rnorm(11)))
})

test_that("twas_join_z works without R", {
    data <- generate_twas_joint_z_data()
    result <- twas_joint_z(data$weights, data$z, X=data$X)
    expect_is(result, "list")
    expect_true(all(names(result) %in% c("Z", "GBJ")))
})

test_that("Check twas_joint_z works with R", {
    data <- generate_twas_joint_z_data()
    result <- twas_joint_z(data$weights, data$z, R=data$R)
    expect_is(result, "list")
    expect_true(all(names(result) %in% c("Z", "GBJ")))
})
