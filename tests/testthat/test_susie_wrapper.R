context("susie_wrapper")

# =============================================================================
# lbf_to_alpha_vector (internal)
# =============================================================================

test_that("lbf_to_alpha_vector converts correctly", {
  lbf <- c(a = -0.5, b = 1.2, c = 0.3)
  alpha <- pecotmr:::lbf_to_alpha_vector(lbf)
  expect_length(alpha, 3)
  expect_equal(sum(alpha), 1, tolerance = 1e-10)
  expect_true(alpha["b"] > alpha["a"])
})

test_that("lbf_to_alpha_vector with prior weights", {
  lbf <- c(a = 1, b = 1, c = 1)  # Equal LBFs
  pw <- c(0.5, 0.25, 0.25)
  alpha <- pecotmr:::lbf_to_alpha_vector(lbf, prior_weights = pw)
  expect_true(alpha[1] > alpha[2])
  expect_equal(sum(alpha), 1, tolerance = 1e-10)
})

test_that("lbf_to_alpha_vector returns zeros for all-zero lbf", {
  lbf <- c(a = 0, b = 0, c = 0)
  alpha <- pecotmr:::lbf_to_alpha_vector(lbf)
  expect_true(all(alpha == 0))
})

test_that("lbf_to_alpha_vector handles single element", {
  lbf <- c(a = 2.0)
  alpha <- pecotmr:::lbf_to_alpha_vector(lbf)
  expect_length(alpha, 1)
  expect_equal(alpha[["a"]], 1.0)
})

test_that("lbf_to_alpha_vector handles very large LBFs without overflow", {
  lbf <- c(a = 500, b = 500.1, c = 499)
  alpha <- pecotmr:::lbf_to_alpha_vector(lbf)
  expect_true(all(is.finite(alpha)))
  expect_equal(sum(alpha), 1, tolerance = 1e-10)
  expect_true(alpha["b"] > alpha["a"])
})

test_that("lbf_to_alpha_vector handles very negative LBFs", {
  lbf <- c(a = -1000, b = -999, c = -1001)
  alpha <- pecotmr:::lbf_to_alpha_vector(lbf)
  expect_true(all(is.finite(alpha)))
  expect_equal(sum(alpha), 1, tolerance = 1e-10)
  expect_true(alpha["b"] > alpha["a"])
})

test_that("lbf_to_alpha_vector with unequal prior weights", {
  lbf <- c(a = 0.5, b = 0.5, c = 0.5)
  pw <- c(0.8, 0.1, 0.1)
  alpha <- pecotmr:::lbf_to_alpha_vector(lbf, prior_weights = pw)
  expect_equal(sum(alpha), 1, tolerance = 1e-10)
  expect_true(alpha[1] > 0.7)
})

# =============================================================================
# lbf_to_alpha (matrix version)
# =============================================================================

test_that("lbf_to_alpha converts log BFs to posteriors", {
  lbf <- matrix(c(0, 3, 2, 1, 4, 0), nrow = 2, ncol = 3)
  alpha <- pecotmr:::lbf_to_alpha(lbf)
  expect_equal(dim(alpha), c(2, 3))
  expect_equal(rowSums(alpha), c(1, 1), tolerance = 1e-10)
  expect_true(alpha[1, 3] > alpha[1, 1])
  expect_true(alpha[2, 1] > alpha[2, 3])
})

test_that("lbf_to_alpha handles uniform lbf", {
  lbf <- matrix(1, nrow = 1, ncol = 5)
  alpha <- pecotmr:::lbf_to_alpha(lbf)
  expect_equal(as.numeric(alpha), rep(0.2, 5), tolerance = 1e-10)
})

test_that("lbf_to_alpha handles single-row matrix", {
  lbf <- matrix(c(1.0, 2.0, 0.5), nrow = 1)
  colnames(lbf) <- paste0("v", 1:3)
  result <- lbf_to_alpha(lbf)
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 3)
  expect_equal(sum(result), 1, tolerance = 1e-10)
})

test_that("lbf_to_alpha handles large matrix", {
  set.seed(42)
  lbf <- matrix(rnorm(100), nrow = 10, ncol = 10)
  colnames(lbf) <- paste0("v", 1:10)
  result <- lbf_to_alpha(lbf)
  expect_equal(dim(result), c(10, 10))
  expect_equal(rowSums(result), rep(1, 10), tolerance = 1e-10)
})

test_that("lbf_to_alpha with mixed zero and nonzero rows", {
  lbf <- matrix(c(0, 0, 0, 1, 2, 3), nrow = 2, byrow = TRUE)
  colnames(lbf) <- paste0("v", 1:3)
  result <- lbf_to_alpha(lbf)
  expect_true(all(result[1, ] == 0))
  expect_equal(sum(result[2, ]), 1, tolerance = 1e-10)
})

# =============================================================================
# get_cs_index (internal)
# =============================================================================

test_that("get_cs_index finds variant in credible set", {
  susie_cs <- list(L1 = c(1, 2, 3), L2 = c(4, 5))
  idx <- pecotmr:::get_cs_index(2, susie_cs)
  expect_equal(unname(idx), 1)
})

test_that("get_cs_index returns NA for variant not in any CS", {
  susie_cs <- list(L1 = c(1, 2), L2 = c(4, 5))
  idx <- pecotmr:::get_cs_index(3, susie_cs)
  expect_true(is.na(idx))
})

test_that("get_cs_index returns all CS indices when variant in multiple", {
  susie_cs <- list(L1 = c(1, 2, 3), L2 = c(2, 4, 5))
  idx <- pecotmr:::get_cs_index(2, susie_cs)
  expect_equal(unname(idx), c(1, 2))
})

test_that("get_cs_index returns all matching CS regardless of size", {
  susie_cs <- list(L1 = c(1, 2, 3, 4, 5), L2 = c(2, 3))
  result <- pecotmr:::get_cs_index(2, susie_cs)
  expect_equal(unname(result), c(1, 2))
})

test_that("get_cs_index handles empty CS list", {
  susie_cs <- list()
  result <- pecotmr:::get_cs_index(1, susie_cs)
  expect_true(is.na(result))
})

test_that("get_cs_index returns correct CS assignment with real susie fit", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 200
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(2, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  y <- X %*% beta + rnorm(n, sd = 0.5)
  fit <- susieR::susie(X, y, L = 5)
  # With beta[1]=2 and sd=0.5, susie should find a CS containing variant 1
  expect_false(is.null(fit$sets$cs))
  idx <- pecotmr:::get_cs_index(1, fit$sets$cs)
  expect_true(is.numeric(unname(idx)))
  expect_true(all(idx >= 1))
})

# =============================================================================
# get_top_variants_idx (internal)
# =============================================================================

test_that("get_top_variants_idx returns combined PIP and CS variants", {
  susie_output <- list(
    pip = c(0.01, 0.15, 0.02, 0.5, 0.01),
    sets = list(cs = list(L1 = c(1, 2)))
  )
  result <- pecotmr:::get_top_variants_idx(susie_output, signal_cutoff = 0.1)
  expect_true(1 %in% result)
  expect_true(2 %in% result)
  expect_true(4 %in% result)
  expect_true(all(result == sort(result)))
})

test_that("get_top_variants_idx with no CS", {
  susie_output <- list(
    pip = c(0.01, 0.5, 0.02, 0.8, 0.01),
    sets = list(cs = NULL)
  )
  result <- pecotmr:::get_top_variants_idx(susie_output, signal_cutoff = 0.1)
  expect_equal(result, c(2, 4))
})

test_that("get_top_variants_idx with all low PIPs", {
  susie_output <- list(
    pip = c(0.01, 0.02, 0.03),
    sets = list(cs = list(L1 = c(1, 2)))
  )
  result <- pecotmr:::get_top_variants_idx(susie_output, signal_cutoff = 0.5)
  expect_equal(result, c(1, 2))
})

test_that("get_top_variants_idx with high cutoff and no CS", {
  susie_output <- list(
    pip = c(0.01, 0.02, 0.03),
    sets = list(cs = NULL)
  )
  result <- pecotmr:::get_top_variants_idx(susie_output, signal_cutoff = 0.5)
  expect_length(result, 0)
})

# =============================================================================
# get_cs_info (internal)
# =============================================================================

test_that("get_cs_info maps variants to CS numbers", {
  susie_cs <- list(L1 = c(1, 2), L3 = c(4, 5, 6))
  top_idx <- c(1, 3, 5)
  result <- pecotmr:::get_cs_info(susie_cs, top_idx)
  # Now returns data.frame(variant_idx, cs_idx) with one row per (variant, CS) pair
  expect_true(is.data.frame(result))
  expect_equal(result$variant_idx, c(1, 3, 5))
  expect_equal(result$cs_idx, c(1L, 0L, 3L))
})

test_that("get_cs_info handles all variants outside CS", {
  susie_cs <- list(L1 = c(1, 2))
  top_idx <- c(5, 6, 7)
  result <- pecotmr:::get_cs_info(susie_cs, top_idx)
  expect_true(is.data.frame(result))
  expect_true(all(result$cs_idx == 0))
})

test_that("get_cs_info reports variant in multiple CSs as multiple rows", {
  susie_cs <- list(L1 = c(1, 2, 3), L3 = c(2, 3, 4))
  top_idx <- c(1, 2, 4)
  result <- pecotmr:::get_cs_info(susie_cs, top_idx)
  expect_true(is.data.frame(result))
  # variant 2 is in both L1 and L3, so it gets two rows
  expect_equal(nrow(result), 4)
  expect_equal(sum(result$variant_idx == 2), 2)
  expect_equal(sort(result$cs_idx[result$variant_idx == 2]), c(1L, 3L))
})

# =============================================================================
# susie_rss_pipeline
# =============================================================================

test_that("susie_rss_pipeline errors on missing z and beta/se", {
  sumstats <- data.frame(x = 1)
  LD_mat <- matrix(1)
  expect_error(susie_rss_pipeline(sumstats, LD_mat), "must have 'z'")
})

test_that("susie_rss_pipeline errors on invalid method", {
  sumstats <- list(z = rnorm(5))
  LD_mat <- diag(5)
  expect_error(susie_rss_pipeline(sumstats, LD_mat, analysis_method = "invalid"))
})

test_that("susie_rss_pipeline runs with single_effect method", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 20
  z <- rnorm(n)
  names(z) <- paste0("chr1:", seq_len(n), ":A:G")
  R <- diag(n)
  colnames(R) <- rownames(R) <- names(z)
  sumstats <- list(z = z)

  result <- susie_rss_pipeline(sumstats, R, analysis_method = "single_effect")
  expect_true(is.list(result))
  expect_true("variant_names" %in% names(result))
  expect_true("susie_result_trimmed" %in% names(result))
  # PIPs should be numeric, in [0,1], and sum to at most 1 (L=1)
  pip <- result$susie_result_trimmed$pip
  expect_true(is.numeric(pip))
  expect_length(pip, n)
  expect_true(all(pip >= 0 & pip <= 1))
  expect_true(sum(pip) <= 1 + 1e-6)
  # Credible sets, if any, should contain valid indices
  cs_list <- result$susie_result_trimmed$sets$cs
  if (!is.null(cs_list)) {
    for (cs in cs_list) {
      expect_true(all(cs >= 1 & cs <= n))
    }
  }
})

test_that("susie_rss_pipeline runs with bayesian_conditional_regression", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 20
  z <- rnorm(n)
  names(z) <- paste0("chr1:", seq_len(n), ":A:G")
  R <- diag(n)
  colnames(R) <- rownames(R) <- names(z)
  sumstats <- list(z = z)

  result <- susie_rss_pipeline(sumstats, R,
    analysis_method = "bayesian_conditional_regression",
    L = 5, max_L = 5
  )
  expect_true(is.list(result))
  expect_true("susie_result_trimmed" %in% names(result))
  pip <- result$susie_result_trimmed$pip
  expect_true(is.numeric(pip))
  expect_length(pip, n)
  expect_true(all(pip >= 0 & pip <= 1))
  # With L=5, sum of PIPs can be up to L
  expect_true(sum(pip) <= 5 + 1e-6)
  cs_list <- result$susie_result_trimmed$sets$cs
  if (!is.null(cs_list)) {
    for (cs in cs_list) {
      expect_true(all(cs >= 1 & cs <= n))
    }
  }
})

test_that("susie_rss_pipeline uses beta/se when z not provided", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 15
  beta <- rnorm(n, sd = 0.1)
  se <- rep(0.1, n)
  names(beta) <- paste0("chr1:", seq_len(n), ":A:G")
  R <- diag(n)
  colnames(R) <- rownames(R) <- names(beta)
  sumstats <- list(beta = beta, se = se)

  result <- susie_rss_pipeline(sumstats, R,
    analysis_method = "susie_rss",
    L = 5, max_L = 5
  )
  expect_true(is.list(result))
  expect_true("susie_result_trimmed" %in% names(result))
  pip <- result$susie_result_trimmed$pip
  expect_true(is.numeric(pip))
  expect_length(pip, n)
  expect_true(all(pip >= 0 & pip <= 1))
  expect_true(sum(pip) <= 5 + 1e-6)
  cs_list <- result$susie_result_trimmed$sets$cs
  if (!is.null(cs_list)) {
    for (cs in cs_list) {
      expect_true(all(cs >= 1 & cs <= n))
    }
  }
})

# =============================================================================
# susie_rss_wrapper
# =============================================================================

test_that("susie_rss_wrapper with L=1 runs single effect", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 10
  R <- diag(p)
  z <- rnorm(p)
  result <- susie_rss_wrapper(z = z, R = R, L = 1)
  expect_true("pip" %in% names(result))
  expect_length(result$pip, p)
  expect_true(is.numeric(result$pip))
  expect_true(all(result$pip >= 0 & result$pip <= 1))
  # L=1 so PIPs sum to at most 1
  expect_true(sum(result$pip) <= 1 + 1e-6)
  if (!is.null(result$sets$cs)) {
    for (cs in result$sets$cs) {
      expect_true(all(cs >= 1 & cs <= p))
    }
  }
})

test_that("susie_rss_wrapper with L equal to max_L", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 15
  z <- rnorm(n)
  R <- diag(n)
  result <- susie_rss_wrapper(z = z, R = R, L = 5, max_L = 5)
  expect_true("pip" %in% names(result))
})

test_that("susie_rss_wrapper dynamic L with no CS found", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 10
  R <- diag(p)
  z <- rep(0.1, p)
  result <- susie_rss_wrapper(z = z, R = R, L = 2, max_L = 10, l_step = 2)
  expect_true("pip" %in% names(result))
})

# =============================================================================
# susie_wrapper
# =============================================================================

test_that("susie_wrapper runs with init_L equal to max_L", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 50
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] * 2 + rnorm(n)

  result <- susie_wrapper(X, y, init_L = 5, max_L = 5)
  expect_true("pip" %in% names(result))
  expect_length(result$pip, p)
})

test_that("susie_wrapper dynamically adjusts L", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 200
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] + rnorm(n)
  result <- susie_wrapper(X, y, init_L = 1, max_L = 10, l_step = 2)
  expect_true("pip" %in% names(result))
  expect_true(!is.null(result$sets))
})

# =============================================================================
# susie_weights
# =============================================================================

test_that("susie_weights returns zeros when fit lacks alpha/mu", {
  fake_fit <- list(pip = rep(0.01, 5))
  result <- susie_weights(susie_fit = fake_fit)
  expect_equal(result, rep(0, 5))
})

test_that("susie_weights checks dimension mismatch", {
  set.seed(42)
  X <- matrix(rnorm(100), 20, 5)
  fake_fit <- list(pip = rep(0.01, 10))
  expect_error(susie_weights(X = X, susie_fit = fake_fit), "Dimension mismatch")
})

# =============================================================================
# susie_ash_weights
# =============================================================================

test_that("susie_ash_weights returns zeros without proper fit structure", {
  fake_fit <- list(pip = rep(0.01, 5))
  result <- susie_ash_weights(susie_ash_fit = fake_fit)
  expect_equal(result, rep(0, 5))
})

# =============================================================================
# susie_inf_weights
# =============================================================================

test_that("susie_inf_weights returns zeros without proper fit structure", {
  fake_fit <- list(pip = rep(0.01, 5))
  result <- susie_inf_weights(susie_inf_fit = fake_fit)
  expect_equal(result, rep(0, 5))
})

# =============================================================================
# glmnet_weights
# =============================================================================

test_that("glmnet_weights produces non-zero weights for correlated data", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 100
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  beta_true <- c(3, -2, rep(0, p - 2))
  y <- X %*% beta_true + rnorm(n)

  w <- glmnet_weights(X, y, alpha = 0.5)
  expect_length(w, p)
  expect_true(any(w != 0))
})

test_that("glmnet_weights handles zero-variance columns", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 100
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  X[, 3] <- 1  # zero variance column
  y <- X[, 1] * 2 + rnorm(n)

  w <- glmnet_weights(X, y, alpha = 1)
  expect_length(w, p)
  expect_equal(w[3], 0)
})

# =============================================================================
# init_prior_sd
# =============================================================================

test_that("init_prior_sd returns n standard deviations", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n_samples <- 50
  p <- 10
  X <- matrix(rnorm(n_samples * p), n_samples, p)
  y <- X[, 1] * 2 + rnorm(n_samples)

  sds <- pecotmr:::init_prior_sd(X, y, n = 15)
  expect_length(sds, 15)
  expect_equal(sds[1], 0)
  expect_true(all(diff(sds) >= 0))
})

# =============================================================================
# adjust_susie_weights
# =============================================================================

# Helper: build a minimal twas_weights_results object with the nested structure
# adjust_susie_weights expects (susie_results / weights paths).
make_adjust_obj <- function(variant_ids, L = 3, ctx = "ctx") {
  set.seed(123)
  p <- length(variant_ids)
  weights_df <- data.frame(
    susie = rnorm(p), enet = rnorm(p),
    row.names = variant_ids, stringsAsFactors = FALSE
  )
  list(
    susie_results = setNames(list(list(
      variant_names = variant_ids,
      susie_result_trimmed = list(
        lbf_variable = matrix(rnorm(L * p), nrow = L, ncol = p),
        mu = matrix(rnorm(L * p), nrow = L, ncol = p),
        X_column_scale_factors = rep(1, p)
      )
    )), ctx),
    weights = setNames(list(weights_df), ctx)
  )
}

# Use non-strand-ambiguous alleles (A2="A", A1="G") so allele_qc keeps them.
adjust_vids <- function(positions = 1:6) {
  paste0("chr1:", positions, ":A:G")
}

# =============================================================================
# susie_rss_wrapper validation (Tier 1)
# =============================================================================

test_that("susie_rss_wrapper errors when neither R nor X is provided", {
  expect_error(
    susie_rss_wrapper(z = rnorm(5)),
    "Either R or X must be provided"
  )
})

test_that("susie_rss_wrapper errors when both R and X are provided", {
  expect_error(
    susie_rss_wrapper(z = rnorm(5), R = diag(5), X = matrix(1, 10, 5)),
    "Only one of R or X should be provided, not both"
  )
})

# =============================================================================
# susie_post_processor: analysis_script and fSuSiE V=NULL branches (Tier 1)
# =============================================================================

# Helper: build a minimal synthetic susie_output for susie_post_processor
make_fake_susie_output <- function(p = 5, L = 3, has_V = TRUE) {
  vnames <- paste0("chr1:", 1:p, ":A:G")
  out <- list(
    pip = setNames(rep(0.01, p), vnames),
    alpha = matrix(1 / p, nrow = L, ncol = p),
    lbf_variable = matrix(0, nrow = L, ncol = p),
    sets = list(
      cs = NULL,
      requested_coverage = 0.95
    ),
    niter = 10
  )
  if (has_V) {
    out$V <- rep(1, L)
  }
  out
}

test_that("susie_post_processor stores analysis_script when load_script returns non-empty", {
  skip_if_not_installed("susieR")
  p <- 5
  fake_output <- make_fake_susie_output(p)
  R <- diag(p)
  colnames(R) <- rownames(R) <- names(fake_output$pip)
  local_mocked_bindings(
    load_script = function() "fake_script_content"
  )
  result <- susie_post_processor(
    fake_output,
    data_x = R,
    data_y = list(z = rnorm(p)),
    mode = "susie_rss"
  )
  expect_equal(result$analysis_script, "fake_script_content")
})

test_that("susie_post_processor uses 1:max_L for eff_idx when V is NULL (fSuSiE)", {
  skip_if_not_installed("susieR")
  p <- 5
  L <- 3
  fake_output <- make_fake_susie_output(p, L = L, has_V = FALSE)
  R <- diag(p)
  colnames(R) <- rownames(R) <- names(fake_output$pip)
  result <- susie_post_processor(
    fake_output,
    data_x = R,
    data_y = list(z = rnorm(p)),
    mode = "susie_rss"
  )
  # With V=NULL, eff_idx = 1:L, so trimmed alpha should keep all L rows
  expect_equal(nrow(result$susie_result_trimmed$alpha), L)
  # V should be NULL in trimmed output
  expect_null(result$susie_result_trimmed$V)
})

# =============================================================================
# susie_wrapper: dynamic-L break when cs is NULL (Tier 2)
# =============================================================================

test_that("susie_wrapper breaks out of dynamic-L loop when cs is NULL", {
  call_count <- 0
  local_mocked_bindings(
    susie = function(X, y, L, ...) {
      call_count <<- call_count + 1
      list(
        sets = list(cs = NULL),
        pip = rep(0.01, ncol(X)),
        time_elapsed = 0
      )
    }
  )
  X <- matrix(rnorm(50), 10, 5)
  y <- rnorm(10)
  expect_message(
    result <- susie_wrapper(X, y, init_L = 5, max_L = 20),
    "Total time elapsed"
  )
  # susie should be called exactly once: cs=NULL triggers immediate break
  expect_equal(call_count, 1)
})

# =============================================================================
# susie_rss_wrapper: dynamic-L increment then break (Tier 2)
# =============================================================================

test_that("susie_rss_wrapper increments L in while-loop then breaks", {
  call_count <- 0
  p <- 10
  local_mocked_bindings(
    susie_rss = function(z, L, ...) {
      call_count <<- call_count + 1
      if (call_count == 1) {
        # First call: length(cs) >= L, so L should increment
        list(
          sets = list(cs = lapply(1:L, function(i) 1:3)),
          pip = rep(0.1, length(z))
        )
      } else {
        # Second call: length(cs) < L, so loop breaks
        list(
          sets = list(cs = list(L1 = 1:3)),
          pip = rep(0.1, length(z))
        )
      }
    }
  )
  z <- rnorm(p)
  result <- susie_rss_wrapper(z = z, R = diag(p), L = 2, max_L = 20, l_step = 3)
  # Called twice: once at L=2 (incremented to 5), once at L=5 (broke)
  expect_equal(call_count, 2)
  expect_true("pip" %in% names(result))
})

# =============================================================================
# susie_rss_wrapper and susie_rss_pipeline X-mode branches (Tier 2)
# =============================================================================

test_that("susie_rss_wrapper passes X (not R) when X is provided", {
  captured_args <- NULL
  p <- 5
  n <- 20
  local_mocked_bindings(
    susie_rss = function(...) {
      args <- list(...)
      captured_args <<- args
      list(
        sets = list(cs = NULL),
        pip = rep(0.1, length(args$z))
      )
    }
  )
  z <- rnorm(p)
  X <- matrix(rnorm(n * p), n, p)
  result <- susie_rss_wrapper(z = z, X = X, L = 5, max_L = 5)
  expect_true("X" %in% names(captured_args))
  expect_null(captured_args$R)
})

test_that("susie_rss_pipeline X-mode passes X to wrapper and computes LD from X for post-processor", {
  skip_if_not_installed("susieR")
  p <- 5
  n <- 20
  z <- rnorm(p)
  vnames <- paste0("chr1:", 1:p, ":A:G")
  names(z) <- vnames
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- vnames

  captured_wrapper_args <- NULL
  captured_pp_data_x <- NULL
  local_mocked_bindings(
    susie_rss_wrapper = function(...) {
      captured_wrapper_args <<- list(...)
      list(
        pip = setNames(rep(0.01, p), vnames),
        alpha = matrix(1 / p, nrow = 5, ncol = p),
        lbf_variable = matrix(0, nrow = 5, ncol = p),
        V = rep(1, 5),
        sets = list(cs = NULL, requested_coverage = 0.95),
        niter = 10
      )
    },
    susie_post_processor = function(susie_output, data_x, ...) {
      captured_pp_data_x <<- data_x
      list(variant_names = vnames)
    }
  )
  result <- susie_rss_pipeline(list(z = z), X_mat = X)
  # Wrapper should have received X, not R
  expect_true("X" %in% names(captured_wrapper_args))
  expect_null(captured_wrapper_args$R)
  # Post-processor should have received a p x p matrix (LD computed from X)
  expect_equal(dim(captured_pp_data_x), c(p, p))
})

# =============================================================================
# susie_rss_pipeline: mixture-panel (list of X) branch (Tier 2)
# =============================================================================

test_that("susie_rss_pipeline computes LD from first panel when X_mat is a list", {
  skip_if_not_installed("susieR")
  p <- 5
  n1 <- 20
  n2 <- 15
  z <- rnorm(p)
  vnames <- paste0("chr1:", 1:p, ":A:G")
  names(z) <- vnames
  X1 <- matrix(rnorm(n1 * p), n1, p)
  X2 <- matrix(rnorm(n2 * p), n2, p)
  colnames(X1) <- colnames(X2) <- vnames
  X_list <- list(panel1 = X1, panel2 = X2)

  captured_pp_data_x <- NULL
  local_mocked_bindings(
    susie_rss_wrapper = function(...) {
      list(
        pip = setNames(rep(0.01, p), vnames),
        alpha = matrix(1 / p, nrow = 5, ncol = p),
        lbf_variable = matrix(0, nrow = 5, ncol = p),
        V = rep(1, 5),
        sets = list(cs = NULL, requested_coverage = 0.95),
        niter = 10
      )
    },
    susie_post_processor = function(susie_output, data_x, ...) {
      captured_pp_data_x <<- data_x
      list(variant_names = vnames)
    }
  )
  result <- susie_rss_pipeline(list(z = z), X_mat = X_list)
  # data_x should be a p x p correlation matrix computed from X1 (first panel)
  expect_equal(dim(captured_pp_data_x), c(p, p))
  # It should be a symmetric matrix (correlation/LD)
  expect_equal(captured_pp_data_x, t(captured_pp_data_x), tolerance = 1e-10)
})

# =============================================================================
# adjust_susie_weights
# =============================================================================

test_that("adjust_susie_weights errors when no variants intersect", {
  vids <- adjust_vids(1:5)
  obj <- make_adjust_obj(vids)
  expect_error(
    adjust_susie_weights(
      obj,
      keep_variants = paste0("chr2:", 1:5, ":A:G"),
      run_allele_qc = FALSE,
      variable_name_obj = c("susie_results", "ctx", "variant_names"),
      susie_obj = c("susie_results", "ctx", "susie_result_trimmed"),
      twas_weights_table = c("weights", "ctx"),
      LD_variants = NULL
    ),
    "No intersected variants"
  )
})

test_that("adjust_susie_weights run_allele_qc=FALSE returns intersect coefs", {
  vids <- adjust_vids(1:6)
  obj <- make_adjust_obj(vids)
  keep <- vids[2:5]
  out <- adjust_susie_weights(
    obj,
    keep_variants = keep, run_allele_qc = FALSE,
    variable_name_obj = c("susie_results", "ctx", "variant_names"),
    susie_obj = c("susie_results", "ctx", "susie_result_trimmed"),
    twas_weights_table = c("weights", "ctx"),
    LD_variants = NULL
  )
  expect_length(out$adjusted_susie_weights, 4)
  expect_equal(out$remained_variants_ids, normalize_variant_id(keep))
  expect_true(all(is.finite(out$adjusted_susie_weights)))
})

test_that("adjust_susie_weights run_allele_qc=FALSE normalizes variant ids before matching", {
  # Object has non-canonical (no chr prefix) variant ids
  vids_raw <- c("1:1:A:G", "1:2:A:G", "1:3:A:G", "1:4:A:G")
  obj <- make_adjust_obj(vids_raw)
  # keep_variants supplied with chr prefix
  keep <- c("chr1:2:A:G", "chr1:3:A:G")
  out <- adjust_susie_weights(
    obj,
    keep_variants = keep, run_allele_qc = FALSE,
    variable_name_obj = c("susie_results", "ctx", "variant_names"),
    susie_obj = c("susie_results", "ctx", "susie_result_trimmed"),
    twas_weights_table = c("weights", "ctx"),
    LD_variants = NULL
  )
  expect_length(out$adjusted_susie_weights, 2)
  expect_equal(out$remained_variants_ids, c("chr1:2:A:G", "chr1:3:A:G"))
})

test_that("adjust_susie_weights run_allele_qc=TRUE returns adjusted xQTL coefs", {
  vids <- adjust_vids(1:5)
  obj <- make_adjust_obj(vids)
  out <- adjust_susie_weights(
    obj,
    keep_variants = vids, run_allele_qc = TRUE, LD_variants = vids,
    variable_name_obj = c("susie_results", "ctx", "variant_names"),
    susie_obj = c("susie_results", "ctx", "susie_result_trimmed"),
    twas_weights_table = c("weights", "ctx"),
    match_min_prop = 0.1
  )
  expect_true(length(out$adjusted_susie_weights) > 0)
  expect_true(all(is.finite(out$adjusted_susie_weights)))
  expect_true(all(grepl("^chr1:", out$remained_variants_ids)))
})

test_that("adjust_susie_weights run_allele_qc=TRUE auto-prepends chrom/pos/A2/A1", {
  vids <- adjust_vids(1:5)
  obj <- make_adjust_obj(vids)
  # Confirm the helper produced a weights matrix WITHOUT chrom/pos/A2/A1 cols
  expect_false(any(c("chrom", "pos", "A2", "A1") %in% colnames(obj$weights$ctx)))
  out <- adjust_susie_weights(
    obj,
    keep_variants = vids, run_allele_qc = TRUE, LD_variants = vids,
    variable_name_obj = c("susie_results", "ctx", "variant_names"),
    susie_obj = c("susie_results", "ctx", "susie_result_trimmed"),
    twas_weights_table = c("weights", "ctx"),
    match_min_prop = 0.1
  )
  expect_true(length(out$adjusted_susie_weights) > 0)
})
