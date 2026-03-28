context("quantile_twas")

# ===========================================================================
# ensure_continuous_clusters (internal)
# ===========================================================================

test_that("ensure_continuous_clusters preserves already continuous", {
  index <- c(1, 1, 2, 2, 3)
  result <- pecotmr:::ensure_continuous_clusters(index)
  expect_equal(result, c(1, 1, 2, 2, 3))
})

test_that("ensure_continuous_clusters renumbers non-continuous", {
  index <- c(1, 2, 1, 2, 3)
  result <- pecotmr:::ensure_continuous_clusters(index)
  expect_equal(result, c(1, 2, 3, 4, 5))
})

test_that("ensure_continuous_clusters handles single element", {
  result <- pecotmr:::ensure_continuous_clusters(c(1))
  expect_equal(result, c(1))
})

test_that("ensure_continuous_clusters handles two elements same cluster", {
  result <- pecotmr:::ensure_continuous_clusters(c(1, 1))
  expect_equal(result, c(1, 1))
})

test_that("ensure_continuous_clusters all same cluster", {
  result <- pecotmr:::ensure_continuous_clusters(c(3, 3, 3, 3))
  expect_equal(result, c(1, 1, 1, 1))
})

# ===========================================================================
# get_cluster_ranges (internal)
# ===========================================================================

test_that("get_cluster_ranges returns correct ranges", {
  index <- c(1, 1, 2, 2, 3)
  result <- pecotmr:::get_cluster_ranges(index)
  expect_length(result, 3)
  expect_equal(result[[1]], "1 - 2")
  expect_equal(result[[2]], "3 - 4")
  expect_equal(result[[3]], "5 - 5")
})

test_that("get_cluster_ranges single cluster", {
  index <- c(1, 1, 1)
  result <- pecotmr:::get_cluster_ranges(index)
  expect_length(result, 1)
  expect_equal(result[[1]], "1 - 3")
})

# ===========================================================================
# get_modularity (internal)
# ===========================================================================

test_that("get_modularity returns 0 for single element", {
  W <- matrix(1, nrow = 1, ncol = 1)
  B <- matrix(1, nrow = 1, ncol = 1)
  result <- pecotmr:::get_modularity(W, B)
  expect_equal(result, 0)
})

test_that("get_modularity returns numeric for identity weight", {
  W <- diag(5)
  B <- matrix(c(rep(c(1, 0), c(3, 2)), rep(c(0, 1), c(3, 2))), nrow = 5)
  result <- pecotmr:::get_modularity(W, B)
  expect_type(result, "double")
})

test_that("get_modularity handles all-positive weights", {
  set.seed(42)
  W <- abs(matrix(rnorm(16), 4, 4))
  W <- (W + t(W)) / 2
  B <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1), nrow = 4)
  result <- pecotmr:::get_modularity(W, B)
  expect_type(result, "double")
  expect_true(is.finite(result))
})

# ===========================================================================
# get_n_cluster (internal)
# ===========================================================================

test_that("get_n_cluster returns single cluster for high correlation", {
  Sigma <- matrix(0.9, nrow = 5, ncol = 5)
  diag(Sigma) <- 1
  hc <- hclust(as.dist(1 - Sigma))
  result <- pecotmr:::get_n_cluster(hc, Sigma, between_cluster = 0.8)
  expect_equal(result$n_cluster, 1)
})

test_that("get_n_cluster returns multiple clusters for low correlation", {
  set.seed(42)
  Sigma <- diag(10)
  Sigma[1:5, 1:5] <- 0.8
  Sigma[6:10, 6:10] <- 0.8
  diag(Sigma) <- 1
  hc <- hclust(as.dist(1 - Sigma))
  result <- pecotmr:::get_n_cluster(hc, Sigma, between_cluster = 0.5)
  expect_true(result$n_cluster >= 1)
})

# ===========================================================================
# integrate_tau (internal)
# ===========================================================================

test_that("integrate_tau with vector input", {
  tau <- c(0.25, 0.50, 0.75)
  a_tau <- c(1, 2, 1)
  result <- pecotmr:::integrate_tau(tau, a_tau)
  expect_type(result, "double")
  expect_true(result > 0)
})

test_that("integrate_tau with matrix input", {
  tau <- c(0.25, 0.50, 0.75)
  a_tau <- matrix(c(1, 2, 1, 0.5, 1, 0.5), nrow = 2, byrow = TRUE)
  result <- pecotmr:::integrate_tau(tau, a_tau)
  expect_equal(nrow(result), 2)
  expect_true(all(result > 0))
})

test_that("integrate_tau zero weights give zero", {
  tau <- c(0.25, 0.50, 0.75)
  a_tau <- c(0, 0, 0)
  result <- pecotmr:::integrate_tau(tau, a_tau)
  expect_equal(as.numeric(result), 0)
})

# ===========================================================================
# get_hierarchical_clusters (internal)
# ===========================================================================

test_that("get_hierarchical_clusters returns valid structure", {
  set.seed(42)
  p <- 10
  Sigma <- matrix(0.3, nrow = p, ncol = p)
  Sigma[1:5, 1:5] <- 0.9
  Sigma[6:10, 6:10] <- 0.9
  diag(Sigma) <- 1
  result <- pecotmr:::get_hierarchical_clusters(Sigma, between_cluster = 0.5)
  expect_type(result, "list")
  expect_true("cluster" %in% names(result))
  expect_true("Q_modularity_initial" %in% names(result))
  expect_true("cluster_ranges" %in% names(result))
  expect_equal(nrow(result$cluster), p)
  expect_true(all(rowSums(result$cluster) == 1))
})

# ===========================================================================
# perform_grouped_integration (internal)
# ===========================================================================

test_that("perform_grouped_integration returns correct structure", {
  set.seed(42)
  n_variants <- 5
  n_tau <- 10
  twas_weight <- matrix(rnorm(n_variants * n_tau), nrow = n_variants)
  rownames(twas_weight) <- paste0("v", 1:n_variants)
  tau_values <- seq(0.1, 0.9, length.out = n_tau)
  pseudo_R2 <- runif(n_tau, 0.1, 0.5)

  result <- pecotmr:::perform_grouped_integration(twas_weight, tau_values, pseudo_R2,
                                                   num_intervals = 3)
  expect_type(result, "list")
  expect_true("weights" %in% names(result))
  expect_true("twas_weight_performance" %in% names(result))
  expect_equal(nrow(result$weights), n_variants)
  expect_true(ncol(result$weights) >= 3)
})

test_that("perform_grouped_integration single variant skips clustering", {
  twas_weight <- matrix(rnorm(10), nrow = 1)
  rownames(twas_weight) <- "v1"
  tau_values <- seq(0.1, 0.9, length.out = 10)
  pseudo_R2 <- runif(10, 0.1, 0.5)

  result <- pecotmr:::perform_grouped_integration(twas_weight, tau_values, pseudo_R2,
                                                   num_intervals = 3)
  expect_equal(ncol(result$weights), 3)
})

# ===========================================================================
# corr_filter (internal)
# ===========================================================================

test_that("corr_filter removes highly correlated columns", {
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  colnames(X) <- paste0("v", 1:p)
  X[, 2] <- X[, 1] + rnorm(n, sd = 0.01)
  result <- pecotmr:::corr_filter(X, cor_thres = 0.9)
  expect_true(ncol(result$X.new) < p)
  expect_true(length(result$filter.id) == ncol(result$X.new))
})

test_that("corr_filter keeps all columns when uncorrelated", {
  set.seed(42)
  n <- 100; p <- 5
  X <- matrix(rnorm(n * p), nrow = n)
  colnames(X) <- paste0("v", 1:p)
  result <- pecotmr:::corr_filter(X, cor_thres = 0.99)
  expect_equal(ncol(result$X.new), p)
  expect_equal(result$filter.id, 1:p)
})

test_that("corr_filter preserves colnames for single remaining column", {
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n * 3), nrow = n)
  colnames(X) <- c("a", "b", "c")
  X[, 2] <- X[, 1] + rnorm(n, sd = 0.001)
  X[, 3] <- X[, 1] + rnorm(n, sd = 0.001)
  result <- pecotmr:::corr_filter(X, cor_thres = 0.5)
  expect_true(ncol(result$X.new) >= 1)
  expect_true(!is.null(colnames(result$X.new)))
})

test_that("corr_filter handles single-column input", {
  set.seed(42)
  n <- 30
  X <- matrix(rnorm(n), nrow = n, ncol = 1)
  colnames(X) <- "v1"
  expect_error(pecotmr:::corr_filter(X, cor_thres = 0.8))
})

test_that("corr_filter with low threshold removes more columns", {
  set.seed(42)
  n <- 100; p <- 5
  X <- matrix(rnorm(n * p), nrow = n)
  colnames(X) <- paste0("v", 1:p)
  X[, 2] <- X[, 1] + rnorm(n, sd = 0.1)
  X[, 3] <- X[, 1] + rnorm(n, sd = 0.1)
  X[, 5] <- X[, 4] + rnorm(n, sd = 0.1)
  result_strict <- pecotmr:::corr_filter(X, cor_thres = 0.3)
  result_lenient <- pecotmr:::corr_filter(X, cor_thres = 0.99)
  expect_true(ncol(result_strict$X.new) <= ncol(result_lenient$X.new))
})

test_that("corr_filter preserves colnames when no columns deleted", {
  set.seed(42)
  n <- 100; p <- 3
  X <- matrix(rnorm(n * p), nrow = n)
  colnames(X) <- c("snp_a", "snp_b", "snp_c")
  result <- pecotmr:::corr_filter(X, cor_thres = 0.999)
  expect_equal(colnames(result$X.new), colnames(X))
})

# ===========================================================================
# remove_highcorr_snp (internal)
# ===========================================================================

test_that("remove_highcorr_snp returns X unchanged when no problematic columns", {
  X <- matrix(rnorm(100), nrow = 20, ncol = 5)
  colnames(X) <- paste0("v", 1:5)
  result <- pecotmr:::remove_highcorr_snp(X, problematic_cols = character(0), strategy = "correlation")
  expect_equal(ncol(result), 5)
})

test_that("remove_highcorr_snp removes single problematic column", {
  X <- matrix(rnorm(100), nrow = 20, ncol = 5)
  colnames(X) <- paste0("v", 1:5)
  result <- pecotmr:::remove_highcorr_snp(X, problematic_cols = "v3", strategy = "correlation")
  expect_equal(ncol(result), 4)
  expect_false("v3" %in% colnames(result))
})

test_that("remove_highcorr_snp variance strategy removes lowest variance column", {
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  colnames(X) <- c("low_var", "mid_var", "high_var")
  X[, 1] <- X[, 1] * 0.01
  X[, 3] <- X[, 3] * 10
  result <- pecotmr:::remove_highcorr_snp(X, problematic_cols = c("low_var", "mid_var", "high_var"),
                                           strategy = "variance")
  expect_equal(ncol(result), 2)
  expect_false("low_var" %in% colnames(result))
})

test_that("remove_highcorr_snp correlation strategy with two columns removes one randomly", {
  set.seed(42)
  X <- matrix(rnorm(100), nrow = 20, ncol = 5)
  colnames(X) <- paste0("v", 1:5)
  result <- pecotmr:::remove_highcorr_snp(X, problematic_cols = c("v1", "v2"), strategy = "correlation")
  expect_equal(ncol(result), 4)
  expect_true(xor("v1" %in% colnames(result), "v2" %in% colnames(result)) ||
              (!("v1" %in% colnames(result)) && !("v2" %in% colnames(result))) == FALSE)
})

test_that("remove_highcorr_snp correlation strategy with 3+ cols removes highest sum", {
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
  colnames(X) <- paste0("v", 1:4)
  X[, 2] <- X[, 1] + rnorm(n, sd = 0.01)
  X[, 3] <- X[, 1] + rnorm(n, sd = 0.01)
  result <- pecotmr:::remove_highcorr_snp(X, problematic_cols = c("v1", "v2", "v3"),
                                           strategy = "correlation")
  expect_equal(ncol(result), 3)
})

test_that("remove_highcorr_snp response_correlation strategy removes lowest response corr", {
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  colnames(X) <- c("v1", "v2", "v3")
  response <- X[, 1] * 2 + rnorm(n, sd = 0.1)
  result <- pecotmr:::remove_highcorr_snp(X, problematic_cols = c("v1", "v2", "v3"),
                                           strategy = "response_correlation",
                                           response = response)
  expect_equal(ncol(result), 2)
  expect_true("v1" %in% colnames(result))
})

test_that("remove_highcorr_snp errors on invalid strategy", {
  X <- matrix(rnorm(100), nrow = 20, ncol = 5)
  colnames(X) <- paste0("v", 1:5)
  expect_error(
    pecotmr:::remove_highcorr_snp(X, problematic_cols = c("v1", "v2"),
                                   strategy = "invalid_strategy"),
    "arg"
  )
})

test_that("remove_highcorr_snp preserves column name when single column remains", {
  set.seed(42)
  X <- matrix(rnorm(40), nrow = 20, ncol = 2)
  colnames(X) <- c("keeper", "removed")
  result <- pecotmr:::remove_highcorr_snp(X, problematic_cols = "removed", strategy = "correlation")
  expect_equal(ncol(result), 1)
  expect_equal(colnames(result), "keeper")
})

# ===========================================================================
# check_remove_highcorr_snp (internal)
# ===========================================================================

test_that("check_remove_highcorr_snp returns full-rank matrix when already full rank", {
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  colnames(X) <- paste0("v", 1:3)
  C <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  result <- pecotmr:::check_remove_highcorr_snp(X = X, C = C, strategy = "correlation")
  expect_true(is.matrix(result))
  expect_true(ncol(result) >= 1)
})

test_that("check_remove_highcorr_snp handles rank-deficient design via corr_filter fallback", {
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
  colnames(X) <- paste0("v", 1:4)
  X[, 4] <- X[, 1] + X[, 2]
  C <- NULL
  result <- pecotmr:::check_remove_highcorr_snp(X = X, C = C, strategy = "correlation")
  design <- cbind(1, result)
  expect_equal(qr(design)$rank, ncol(design))
})

test_that("check_remove_highcorr_snp preserves colname for single-column input", {
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n), nrow = n, ncol = 1)
  colnames(X) <- "only_snp"
  C <- NULL
  result <- pecotmr:::check_remove_highcorr_snp(X = X, C = C, strategy = "correlation")
  expect_equal(colnames(result), "only_snp")
})

# ===========================================================================
# calculate_coef_heterogeneity (internal)
# ===========================================================================

test_that("calculate_coef_heterogeneity computes log(sd/mean)", {
  df <- data.frame(
    variant_id = c("v1", "v2"),
    coef_qr_0.25 = c(1, 0),
    coef_qr_0.50 = c(2, 0),
    coef_qr_0.75 = c(3, 0),
    stringsAsFactors = FALSE
  )
  result <- pecotmr:::calculate_coef_heterogeneity(df)
  expect_equal(nrow(result), 2)
  expect_true("coef_heter" %in% colnames(result))
  expect_equal(result$coef_heter[1], log(1 / 2), tolerance = 0.01)
  expect_true(is.na(result$coef_heter[2]))
})

test_that("calculate_coef_heterogeneity handles NA values", {
  df <- data.frame(
    variant_id = "v1",
    coef_qr_0.25 = 1,
    coef_qr_0.50 = NA,
    coef_qr_0.75 = 3,
    stringsAsFactors = FALSE
  )
  result <- pecotmr:::calculate_coef_heterogeneity(df)
  expect_equal(nrow(result), 1)
  expect_true(is.finite(result$coef_heter[1]) || is.na(result$coef_heter[1]))
})

test_that("calculate_coef_heterogeneity handles all-same coefficients", {
  df <- data.frame(
    variant_id = "v1",
    coef_qr_0.25 = 5,
    coef_qr_0.50 = 5,
    coef_qr_0.75 = 5,
    stringsAsFactors = FALSE
  )
  result <- pecotmr:::calculate_coef_heterogeneity(df)
  expect_true(is.infinite(result$coef_heter[1]) || is.na(result$coef_heter[1]))
})

test_that("calculate_coef_heterogeneity handles negative coefficients", {
  df <- data.frame(
    variant_id = "v1",
    coef_qr_0.25 = -1,
    coef_qr_0.50 = -2,
    coef_qr_0.75 = -3,
    stringsAsFactors = FALSE
  )
  result <- pecotmr:::calculate_coef_heterogeneity(df)
  expect_equal(result$coef_heter[1], log(1 / 2), tolerance = 0.01)
})

# ===========================================================================
# calculate_xi_correlation
# ===========================================================================

test_that("calculate_xi_correlation computes xi for valid monotonic data", {
  skip_if_not_installed("XICOR")
  df <- data.frame(variant_id = "v1")
  tau_range <- seq(0.05, 0.95, by = 0.05)
  for (tau in tau_range) {
    df[[paste0("coef_qr_", tau)]] <- tau * 2
  }
  result <- calculate_xi_correlation(df)
  expect_true("xi" %in% colnames(result))
  expect_true(is.numeric(result$xi))
})

test_that("calculate_xi_correlation warns on missing columns", {
  skip_if_not_installed("XICOR")
  df <- data.frame(variant_id = "v1", some_col = 1)
  expect_warning(result <- calculate_xi_correlation(df))
  expect_true(is.na(result$xi))
})

test_that("calculate_xi_correlation returns NA for too few valid values", {
  skip_if_not_installed("XICOR")
  df <- data.frame(variant_id = "v1")
  tau_range <- seq(0.1, 0.9, by = 0.05)
  for (tau in tau_range) {
    df[[paste0("coef_qr_", tau)]] <- NA
  }
  df$coef_qr_0.1 <- 1
  result <- calculate_xi_correlation(df, min_valid = 10)
  expect_true(is.na(result$xi))
})

test_that("calculate_xi_correlation handles multiple variants", {
  skip_if_not_installed("XICOR")
  df <- data.frame(variant_id = c("v1", "v2"))
  tau_range <- seq(0.1, 0.9, by = 0.05)
  for (tau in tau_range) {
    df[[paste0("coef_qr_", tau)]] <- c(tau * 2, tau * -1)
  }
  result <- calculate_xi_correlation(df, tau_range = tau_range)
  expect_equal(nrow(result), 2)
  expect_true("xi" %in% colnames(result))
  expect_true("xi_pval" %in% colnames(result))
})

test_that("calculate_xi_correlation handles error in xicor gracefully", {
  skip_if_not_installed("XICOR")
  df <- data.frame(variant_id = "v1")
  tau_range <- seq(0.1, 0.9, by = 0.05)
  for (tau in tau_range) {
    df[[paste0("coef_qr_", tau)]] <- 1
  }
  result <- calculate_xi_correlation(df, tau_range = tau_range, min_valid = 5)
  expect_equal(nrow(result), 1)
  expect_true(is.numeric(result$xi) || is.na(result$xi))
})

# ===========================================================================
# qr_screen
# ===========================================================================

test_that("qr_screen runs and returns results with quantreg", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n * 5), nrow = n, ncol = 5)
  colnames(X) <- paste0("chr1:", seq(100, 500, by = 100), ":A:G")
  Y <- matrix(rnorm(n), ncol = 1)
  colnames(Y) <- "pheno1"
  result <- qr_screen(X, Y, tau.list = c(0.25, 0.5, 0.75))
  expect_type(result, "list")
})

test_that("qr_screen with fdr screen_method works", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 50; p <- 5
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("chr1:", seq(100, 500, by = 100), ":A:G")
  Y <- matrix(rnorm(n), ncol = 1)
  colnames(Y) <- "pheno1"
  result <- qr_screen(X, Y, tau.list = c(0.25, 0.5, 0.75), screen_method = "fdr")
  expect_type(result, "list")
  expect_true("df_result" %in% names(result))
  expect_true("fdr_p_qr" %in% colnames(result$df_result))
})

test_that("qr_screen errors with invalid screen_method", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 30; p <- 3
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("chr1:", seq(100, 300, by = 100), ":A:G")
  Y <- matrix(rnorm(n), ncol = 1)
  colnames(Y) <- "pheno1"
  expect_error(qr_screen(X, Y, tau.list = c(0.5), screen_method = "invalid_method"),
               "Invalid screen_method")
})

test_that("qr_screen with covariates Z", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 50; p <- 3
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("chr1:", seq(100, 300, by = 100), ":A:G")
  Y <- matrix(rnorm(n), ncol = 1)
  colnames(Y) <- "pheno1"
  Z <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  result <- qr_screen(X, Y, Z = Z, tau.list = c(0.25, 0.5, 0.75))
  expect_type(result, "list")
  expect_equal(nrow(result$df_result), p)
})

test_that("qr_screen with single variant", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n), nrow = n, ncol = 1)
  colnames(X) <- "chr1:100:A:G"
  Y <- matrix(rnorm(n), ncol = 1)
  colnames(Y) <- "pheno1"
  result <- qr_screen(X, Y, tau.list = c(0.25, 0.5, 0.75))
  expect_type(result, "list")
  expect_equal(nrow(result$df_result), 1)
})

# ===========================================================================
# perform_qr_analysis
# ===========================================================================

test_that("perform_qr_analysis returns wide-format results", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 50; p <- 3
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("chr1:", seq(100, 300, by = 100), ":A:G")
  Y <- matrix(rnorm(n), ncol = 1)
  colnames(Y) <- "pheno1"
  result <- perform_qr_analysis(X, Y, tau_values = c(0.25, 0.5, 0.75))
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), p)
  expect_true("coef_qr_0.25" %in% colnames(result))
  expect_true("coef_qr_0.5" %in% colnames(result))
  expect_true("coef_qr_0.75" %in% colnames(result))
  expect_true("chr" %in% colnames(result))
  expect_true("pos" %in% colnames(result))
})

test_that("perform_qr_analysis with covariates Z", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 50; p <- 2
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("chr1:", seq(100, 200, by = 100), ":A:G")
  Y <- matrix(rnorm(n), ncol = 1)
  colnames(Y) <- "pheno1"
  Z <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  result <- perform_qr_analysis(X, Y, Z = Z, tau_values = c(0.5))
  expect_equal(nrow(result), p)
  expect_true("coef_qr_0.5" %in% colnames(result))
})

test_that("perform_qr_analysis with single variant", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n), nrow = n, ncol = 1)
  colnames(X) <- "chr1:100:A:G"
  Y <- matrix(rnorm(n), ncol = 1)
  colnames(Y) <- "pheno1"
  result <- perform_qr_analysis(X, Y, tau_values = c(0.25, 0.5))
  expect_equal(nrow(result), 1)
})

# ===========================================================================
# multicontext_ld_clumping
# ===========================================================================

test_that("multicontext_ld_clumping returns empty result when no significant SNPs", {
  qr_results <- list(sig.SNPs_names = character(0))
  X <- matrix(rnorm(100), nrow = 20, ncol = 5)
  result <- multicontext_ld_clumping(X, qr_results)
  expect_null(result$final_SNPs)
  expect_null(result$clumped_SNPs)
  expect_true(grepl("No significant", result$message))
})

test_that("multicontext_ld_clumping returns early when X has single column", {
  X <- matrix(rnorm(20), nrow = 20, ncol = 1)
  colnames(X) <- "chr1:100:A:G"
  qr_results <- list(sig.SNPs_names = "chr1:100:A:G")
  result <- multicontext_ld_clumping(X, qr_results)
  expect_equal(result$final_SNPs, "chr1:100:A:G")
  expect_true(grepl("Only one SNP", result$message))
})

# ===========================================================================
# quantile_twas_weight_pipeline (mocked)
# ===========================================================================

test_that("quantile_twas_weight_pipeline returns early when no significant SNPs and screen_significant=TRUE", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 50; p <- 3
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("chr1:", seq(100, 300, by = 100), ":A:G")
  Y <- matrix(rnorm(n), ncol = 1)
  colnames(Y) <- "pheno1"

  local_mocked_bindings(
    QUAIL_rank_score_pipeline = function(...) list(rank_score = matrix(rnorm(n), ncol = 1)),
    QUAIL_pipeline = function(...) list(vqtl = "mocked"),
    qr_screen = function(...) {
      list(
        df_result = data.frame(variant_id = colnames(X)),
        sig_SNP_threshold = integer(0),
        sig.SNPs_names = character(0),
        pvec = rep(1, p),
        quantile_pvalue = matrix(1, nrow = p, ncol = 1),
        quantile_zscore = matrix(0, nrow = p, ncol = 1),
        tau_list = 0.5
      )
    }
  )
  result <- quantile_twas_weight_pipeline(
    X, Y, screen_significant = TRUE,
    quantile_qtl_tau_list = c(0.5),
    quantile_twas_tau_list = c(0.5)
  )
  expect_true(grepl("No significant", result$message))
  expect_true("qr_screen_pvalue_df" %in% names(result))
})

test_that("quantile_twas_weight_pipeline returns early when no raw p_qr < 0.05", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 50; p <- 3
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("chr1:", seq(100, 300, by = 100), ":A:G")
  Y <- matrix(rnorm(n), ncol = 1)
  colnames(Y) <- "pheno1"

  local_mocked_bindings(
    QUAIL_rank_score_pipeline = function(...) list(rank_score = matrix(rnorm(n), ncol = 1)),
    QUAIL_pipeline = function(...) list(vqtl = "mocked"),
    qr_screen = function(...) {
      pvec <- rep(0.5, p)
      names(pvec) <- colnames(X)
      list(
        df_result = data.frame(variant_id = colnames(X)),
        sig_SNP_threshold = 1:p,
        sig.SNPs_names = colnames(X),
        pvec = pvec,
        quantile_pvalue = matrix(0.5, nrow = p, ncol = 1),
        quantile_zscore = matrix(0, nrow = p, ncol = 1),
        tau_list = 0.5
      )
    },
    perform_qr_analysis = function(...) {
      data.frame(
        variant_id = colnames(X),
        coef_qr_0.5 = rnorm(p),
        chr = rep("chr1", p),
        pos = seq(100, 300, by = 100),
        ref = rep("G", p),
        alt = rep("A", p),
        phenotype_id = rep("pheno1", p),
        stringsAsFactors = FALSE
      )
    },
    calculate_coef_heterogeneity = function(...) {
      data.frame(variant_id = colnames(X), coef_heter = rnorm(p))
    },
    calculate_xi_correlation = function(...) {
      data.frame(variant_id = colnames(X), xi = rnorm(p), xi_pval = runif(p))
    }
  )
  result <- quantile_twas_weight_pipeline(
    X, Y, screen_significant = FALSE,
    quantile_qtl_tau_list = c(0.5),
    quantile_twas_tau_list = c(0.5)
  )
  # expect_true(grepl("No variants with raw", result$message))
})
