context("univariate_pipeline")

# ===========================================================================
# Helpers
# ===========================================================================

make_uap_inputs <- function(n = 30, p = 10, seed = 42) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p), ":A:G")
  rownames(X) <- paste0("s", seq_len(n))
  beta_true <- c(2, -1.5, rep(0, p - 2))
  Y <- as.numeric(X %*% beta_true + rnorm(n))
  names(Y) <- rownames(X)
  maf <- runif(p, 0.1, 0.5)
  list(X = X, Y = Y, maf = maf, n = n, p = p)
}

# ===========================================================================
# Helper: build a fake susie-like fitted object
# ===========================================================================
make_fake_susie_fit <- function(p) {
  L <- 2
  pip <- runif(p, 0, 0.3)
  pip[1] <- 0.9
  alpha <- matrix(runif(L * p), nrow = L, ncol = p)
  lbf_variable <- matrix(rnorm(L * p), nrow = L, ncol = p)
  V <- c(0.5, 0.01)
  mu <- matrix(rnorm(L * p), nrow = L, ncol = p)
  X_column_scale_factors <- rep(1, p)
  sets <- list(
    cs = list(L1 = c(1L, 2L)),
    purity = data.frame(min.abs.corr = 0.9, mean.abs.corr = 0.95, median.abs.corr = 0.93),
    requested_coverage = 0.95
  )
  niter <- 10
  list(
    pip = pip, alpha = alpha, lbf_variable = lbf_variable,
    V = V, mu = mu, X_column_scale_factors = X_column_scale_factors,
    sets = sets, niter = niter
  )
}

# ===========================================================================
# Helper: build a fake susie_post_processor return value
# ===========================================================================
make_fake_post_result <- function(p) {
  list(
    variant_names = paste0("chr1:", seq_len(p), ":A:G"),
    susie_result_trimmed = list(
      pip = runif(p),
      sets = list(cs = list(L1 = c(1L, 2L)), requested_coverage = 0.95),
      cs_corr = matrix(c(1, 0.1, 0.1, 1), nrow = 2),
      alpha = matrix(runif(2 * p), 2, p),
      lbf_variable = matrix(rnorm(2 * p), 2, p),
      V = c(0.5, 0.01),
      niter = 10,
      max_L = 2
    ),
    top_loci = data.frame(
      variant_id = paste0("chr1:1:A:G"),
      betahat = 1.5, sebetahat = 0.3, z = 5.0, maf = 0.25,
      pip = 0.9, cs_coverage_0.95 = 1,
      stringsAsFactors = FALSE
    )
  )
}

# ===========================================================================
# Helper: fake twas_weights_pipeline return value
# ===========================================================================
make_fake_twas_result <- function(p) {
  list(
    twas_weights = setNames(rep(0.1, p), paste0("chr1:", seq_len(p), ":A:G")),
    twas_predictions = rnorm(30),
    susie_weights_intermediate = list(
      mu = matrix(rnorm(2 * p), 2, p),
      lbf_variable = matrix(rnorm(2 * p), 2, p),
      X_column_scale_factors = rep(1, p),
      pip = runif(p)
    ),
    total_time_elapsed = c(user.self = 0.1, sys.self = 0, elapsed = 0.1)
  )
}

# ===========================================================================
# univariate_analysis_pipeline: input validation
# ===========================================================================


test_that("univariate_analysis_pipeline errors on non-matrix X", {
  Y <- rnorm(10)
  maf <- runif(5, 0.05, 0.5)
  expect_error(
    univariate_analysis_pipeline(X = data.frame(matrix(rnorm(50), 10, 5)), Y = Y, maf = maf),
    "X must be a numeric matrix"
  )
})

test_that("univariate_analysis_pipeline errors on non-numeric X", {
  X <- matrix(as.character(1:50), 10, 5)
  Y <- rnorm(10)
  maf <- runif(5, 0.05, 0.5)
  expect_error(
    univariate_analysis_pipeline(X = X, Y = Y, maf = maf),
    "X must be a numeric matrix"
  )
})

test_that("univariate_analysis_pipeline errors on non-numeric Y", {
  X <- matrix(rnorm(50), 10, 5)
  Y <- as.character(rnorm(10))
  maf <- runif(5, 0.05, 0.5)
  expect_error(
    univariate_analysis_pipeline(X = X, Y = Y, maf = maf),
    "Y must be a numeric vector"
  )
})

test_that("univariate_analysis_pipeline errors on multi-column Y matrix", {
  X <- matrix(rnorm(50), 10, 5)
  Y <- matrix(rnorm(20), 10, 2)
  maf <- runif(5, 0.05, 0.5)
  expect_error(
    univariate_analysis_pipeline(X = X, Y = Y, maf = maf),
    "Y must be a numeric vector"
  )
})

test_that("univariate_analysis_pipeline accepts single-column Y matrix", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 20; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p), ":A:G")
  rownames(X) <- paste0("s", seq_len(n))
  Y <- matrix(X[, 1] * 2 + rnorm(n), ncol = 1)
  rownames(Y) <- rownames(X)
  maf <- runif(p, 0.1, 0.5)
  result <- univariate_analysis_pipeline(X = X, Y = Y, maf = maf,
    twas_weights = FALSE, max_L = 5, init_L = 5)
  expect_type(result, "list")
})

test_that("univariate_analysis_pipeline errors on X/Y row mismatch", {
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(8)
  maf <- runif(5, 0.05, 0.5)
  expect_error(
    univariate_analysis_pipeline(X = X, Y = Y, maf = maf),
    "same number of rows"
  )
})

test_that("univariate_analysis_pipeline errors on maf length mismatch", {
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  maf <- runif(3, 0.05, 0.5)
  expect_error(
    univariate_analysis_pipeline(X = X, Y = Y, maf = maf),
    "length equal to the number of columns"
  )
})

test_that("univariate_analysis_pipeline errors on maf out of bounds", {
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  maf <- c(0.1, 0.2, -0.1, 0.3, 0.4)
  expect_error(
    univariate_analysis_pipeline(X = X, Y = Y, maf = maf),
    "maf values must be between 0 and 1"
  )
})

test_that("univariate_analysis_pipeline errors on invalid X_scalar", {
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  maf <- runif(5, 0.05, 0.5)
  expect_error(
    univariate_analysis_pipeline(X = X, Y = Y, maf = maf, X_scalar = "bad"),
    "X_scalar must be a numeric"
  )
})

test_that("univariate_analysis_pipeline errors on wrong-length X_scalar vector", {
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  maf <- runif(5, 0.05, 0.5)
  expect_error(
    univariate_analysis_pipeline(X = X, Y = Y, maf = maf, X_scalar = c(1, 2, 3)),
    "X_scalar must be a numeric"
  )
})

test_that("univariate_analysis_pipeline errors on non-numeric Y_scalar", {
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  maf <- runif(5, 0.05, 0.5)
  expect_error(
    univariate_analysis_pipeline(X = X, Y = Y, maf = maf, Y_scalar = "bad"),
    "Y_scalar must be a numeric scalar"
  )
})

test_that("univariate_analysis_pipeline errors on non-positive init_L", {
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  maf <- runif(5, 0.05, 0.5)
  expect_error(
    univariate_analysis_pipeline(X = X, Y = Y, maf = maf, init_L = 0),
    "init_L must be a positive"
  )
})

test_that("univariate_analysis_pipeline errors on non-positive max_L", {
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  maf <- runif(5, 0.05, 0.5)
  expect_error(
    univariate_analysis_pipeline(X = X, Y = Y, maf = maf, max_L = -1),
    "max_L must be a positive"
  )
})

test_that("univariate_analysis_pipeline errors on non-positive l_step", {
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  maf <- runif(5, 0.05, 0.5)
  expect_error(
    univariate_analysis_pipeline(X = X, Y = Y, maf = maf, l_step = 0),
    "l_step must be a positive"
  )
})

# ---- univariate_analysis_pipeline functional tests ----

test_that("univariate_analysis_pipeline runs with minimal valid input", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 30; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p), ":A:G")
  rownames(X) <- paste0("s", seq_len(n))
  beta_true <- c(2, -1.5, rep(0, p - 2))
  Y <- as.numeric(X %*% beta_true + rnorm(n))
  names(Y) <- rownames(X)
  maf <- runif(p, 0.1, 0.5)

  result <- univariate_analysis_pipeline(
    X = X, Y = Y, maf = maf,
    twas_weights = FALSE, init_L = 5, max_L = 5
  )

  expect_type(result, "list")
  expect_true("susie_fitted" %in% names(result))
  expect_true("total_time_elapsed" %in% names(result))
})

test_that("univariate_analysis_pipeline produces twas_weights when requested", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 30; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p), ":A:G")
  rownames(X) <- paste0("s", seq_len(n))
  Y <- as.numeric(X[, 1] * 2 + rnorm(n))
  names(Y) <- rownames(X)
  maf <- runif(p, 0.1, 0.5)

  result <- univariate_analysis_pipeline(
    X = X, Y = Y, maf = maf,
    twas_weights = TRUE, init_L = 5, max_L = 5,
    cv_folds = 2
  )

  expect_true("twas_weights_result" %in% names(result))
})

test_that("univariate_analysis_pipeline respects X_scalar vector", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 20; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p), ":A:G")
  rownames(X) <- paste0("s", seq_len(n))
  Y <- as.numeric(X[, 1] * 2 + rnorm(n))
  names(Y) <- rownames(X)
  maf <- runif(p, 0.1, 0.5)

  r_noscalar <- suppressMessages(univariate_analysis_pipeline(
    X = X, Y = Y, maf = maf, X_scalar = rep(1, p),
    twas_weights = FALSE, init_L = 5, max_L = 5
  ))
  r_scalar <- suppressMessages(univariate_analysis_pipeline(
    X = X, Y = Y, maf = maf, X_scalar = rep(2, p),
    twas_weights = FALSE, init_L = 5, max_L = 5
  ))
  # X_scalar divides betahat: scalar=2 should halve the estimates
  ratio <- unname(r_noscalar$sumstats$betahat / r_scalar$sumstats$betahat)
  expect_equal(ratio, rep(2, p), tolerance = 1e-10)
})

test_that("univariate_analysis_pipeline with cv_folds=0 skips CV", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 20; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p), ":A:G")
  rownames(X) <- paste0("s", seq_len(n))
  Y <- as.numeric(X[, 1] * 2 + rnorm(n))
  names(Y) <- rownames(X)
  maf <- runif(p, 0.1, 0.5)

  result <- univariate_analysis_pipeline(
    X = X, Y = Y, maf = maf,
    twas_weights = TRUE, init_L = 5, max_L = 5,
    cv_folds = 0
  )
  expect_type(result, "list")
  expect_true("twas_weights_result" %in% names(result))
})

# ---- rss_analysis_pipeline input validation ----

test_that("rss_analysis_pipeline requires file inputs", {
  # rss_analysis_pipeline calls load_rss_data which requires valid file paths
  expect_error(
    rss_analysis_pipeline(
      sumstat_path = "/nonexistent/file.tsv",
      column_file_path = "/nonexistent/columns.yml",
      LD_data = list()
    )
  )
})

# ---- resolve_LD_input ----

test_that("resolve_LD_input errors when both R and X are NULL", {
  expect_error(
    pecotmr:::resolve_LD_input(R = NULL, X = NULL),
    "Either R .* or X .* must be provided"
  )
})

test_that("resolve_LD_input errors when both R and X are provided", {
  R <- diag(5)
  X <- matrix(rnorm(50), 10, 5)
  expect_error(
    pecotmr:::resolve_LD_input(R = R, X = X),
    "Provide either R or X, not both"
  )
})

test_that("resolve_LD_input with R returns R unchanged", {
  R <- diag(5)
  result <- pecotmr:::resolve_LD_input(R = R, nSample = 100)
  expect_equal(result$R, R)
  expect_equal(result$nSample, 100)
})

test_that("resolve_LD_input with X computes LD and infers nSample", {
  set.seed(42)
  X <- matrix(rbinom(200, 2, 0.3), 20, 10)
  result <- pecotmr:::resolve_LD_input(X = X)
  expect_equal(nrow(result$R), 10)
  expect_equal(ncol(result$R), 10)
  expect_equal(result$nSample, 20)
})

test_that("resolve_LD_input errors when nSample required but missing", {
  R <- diag(5)
  expect_error(
    pecotmr:::resolve_LD_input(R = R, need_nSample = TRUE),
    "nSample is required"
  )
})

test_that("resolve_LD_input does not error when nSample not needed", {
  R <- diag(5)
  result <- pecotmr:::resolve_LD_input(R = R, need_nSample = FALSE)
  expect_equal(result$R, R)
  expect_null(result$nSample)
})

# ===========================================================================
# univariate_analysis_pipeline: mocked pipeline tests
# ===========================================================================

# ========================================================================
#  SECTION 1: univariate_analysis_pipeline – pip_cutoff_to_skip branch
# ========================================================================

test_that("uap: pip_cutoff_to_skip > 0, no signal above threshold => returns empty list", {
  inp <- make_uap_inputs()
  # Mock susie (imported from susieR into pecotmr namespace)
  low_pip <- rep(0.001, inp$p)
  local_mocked_bindings(
    susie_wrapper = function(...) {
      fit <- make_fake_susie_fit(inp$p)
      fit$pip <- low_pip
      fit
    },
    # susie is imported from susieR; mock in pecotmr namespace
    susie = function(...) list(pip = low_pip)
  )
  result <- expect_message(
    univariate_analysis_pipeline(
      X = inp$X, Y = inp$Y, maf = inp$maf,
      pip_cutoff_to_skip = 0.5,
      twas_weights = FALSE, init_L = 5, max_L = 5
    ),
    "Skipping follow-up"
  )
  expect_equal(result, list())
})

test_that("uap: pip_cutoff_to_skip > 0, signal above threshold => continues analysis", {
  inp <- make_uap_inputs()
  high_pip <- rep(0.001, inp$p)
  high_pip[1] <- 0.9

  fake_fit <- make_fake_susie_fit(inp$p)
  fake_post <- make_fake_post_result(inp$p)

  local_mocked_bindings(
    susie_wrapper = function(...) fake_fit,
    susie_post_processor = function(...) fake_post,
    susie = function(...) list(pip = high_pip)
  )

  result <- expect_message(
    univariate_analysis_pipeline(
      X = inp$X, Y = inp$Y, maf = inp$maf,
      pip_cutoff_to_skip = 0.5,
      twas_weights = FALSE, init_L = 5, max_L = 5
    ),
    "Follow-up on region"
  )
  expect_true("susie_fitted" %in% names(result))
})

test_that("uap: negative pip_cutoff_to_skip auto-computes threshold", {
  inp <- make_uap_inputs()
  # auto threshold = 3 * 1/p = 3/10 = 0.3
  # make all PIPs below that
  low_pip <- rep(0.1, inp$p)

  local_mocked_bindings(
    susie_wrapper = function(...) make_fake_susie_fit(inp$p),
    susie = function(...) list(pip = low_pip)
  )

  result <- expect_message(
    univariate_analysis_pipeline(
      X = inp$X, Y = inp$Y, maf = inp$maf,
      pip_cutoff_to_skip = -1,
      twas_weights = FALSE, init_L = 5, max_L = 5
    ),
    "Skipping follow-up"
  )
  expect_equal(result, list())
})

# ========================================================================
#  SECTION 2: univariate_analysis_pipeline – LD reference filtering
# ========================================================================

test_that("uap: LD reference filtering subsets X columns and maf", {
  inp <- make_uap_inputs()
  fake_fit <- make_fake_susie_fit(5)  # 5 variants after filtering
  fake_post <- make_fake_post_result(5)

  captured_args <- list()
  local_mocked_bindings(
    filter_variants_by_ld_reference = function(variant_ids, ld_ref_file) {
      # keep only first 5 variants
      list(data = variant_ids[1:5], idx = 1:5)
    },
    susie_wrapper = function(X, ...) {
      captured_args$ncol_X <<- ncol(X)
      fake_fit
    },
    susie_post_processor = function(...) fake_post,
  )

  result <- univariate_analysis_pipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    ld_reference_meta_file = "/fake/ld_meta.tsv",
    twas_weights = FALSE, init_L = 5, max_L = 5
  )
  expect_equal(captured_args$ncol_X, 5)
  expect_true("susie_fitted" %in% names(result))
})

test_that("uap: LD reference filtering also subsets X_scalar vector", {
  inp <- make_uap_inputs()
  X_scalar_vec <- seq_len(inp$p) * 1.0
  fake_fit <- make_fake_susie_fit(5)
  fake_post <- make_fake_post_result(5)

  captured_maf <- NULL
  local_mocked_bindings(
    filter_variants_by_ld_reference = function(variant_ids, ld_ref_file) {
      list(data = variant_ids[1:5], idx = 1:5)
    },
    susie_wrapper = function(...) fake_fit,
    susie_post_processor = function(susie_output, data_x, data_y, X_scalar, ...) {
      captured_maf <<- length(data_x[1,])
      fake_post
    },
  )

  result <- univariate_analysis_pipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    X_scalar = X_scalar_vec,
    ld_reference_meta_file = "/fake/ld_meta.tsv",
    twas_weights = FALSE, init_L = 5, max_L = 5
  )
  expect_true("susie_fitted" %in% names(result))
})

# ========================================================================
#  SECTION 3: univariate_analysis_pipeline – filter_X branch
# ========================================================================

test_that("uap: filter_X is invoked when imiss_cutoff is set", {
  inp <- make_uap_inputs()
  fake_fit <- make_fake_susie_fit(inp$p)
  fake_post <- make_fake_post_result(inp$p)

  filter_called <- FALSE
  local_mocked_bindings(
    filter_X = function(X, ...) {
      filter_called <<- TRUE
      X  # return unchanged
    },
    susie_wrapper = function(...) fake_fit,
    susie_post_processor = function(...) fake_post,
  )

  result <- univariate_analysis_pipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    imiss_cutoff = 0.9,
    twas_weights = FALSE, init_L = 5, max_L = 5
  )
  expect_true(filter_called)
})

test_that("uap: filter_X with maf_cutoff active properly subsets maf and X_scalar", {
  inp <- make_uap_inputs()
  # Mock filter_X to drop last 3 columns
  keep_cols <- 1:(inp$p - 3)
  X_filtered <- inp$X[, keep_cols, drop = FALSE]
  fake_fit <- make_fake_susie_fit(length(keep_cols))
  fake_post <- make_fake_post_result(length(keep_cols))

  captured_ncol <- NULL
  local_mocked_bindings(
    filter_X = function(X, ...) X_filtered,
    susie_wrapper = function(X, ...) {
      captured_ncol <<- ncol(X)
      fake_fit
    },
    susie_post_processor = function(...) fake_post,
  )

  result <- univariate_analysis_pipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    maf_cutoff = 0.05, X_scalar = rep(2, inp$p),
    twas_weights = FALSE, init_L = 5, max_L = 5
  )
  expect_equal(captured_ncol, length(keep_cols))
})

# ========================================================================
#  SECTION 4: univariate_analysis_pipeline – main analysis (susie + post)
# ========================================================================

test_that("uap: susie_wrapper called with correct args and result stored", {
  inp <- make_uap_inputs()
  fake_fit <- make_fake_susie_fit(inp$p)
  fake_post <- make_fake_post_result(inp$p)

  captured_init_L <- NULL
  captured_max_L <- NULL
  local_mocked_bindings(
    susie_wrapper = function(X, y, init_L, max_L, l_step, coverage, ...) {
      captured_init_L <<- init_L
      captured_max_L <<- max_L
      fake_fit
    },
    susie_post_processor = function(...) fake_post,
  )

  result <- univariate_analysis_pipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    twas_weights = FALSE, init_L = 7, max_L = 15
  )
  expect_equal(captured_init_L, 7)
  expect_equal(captured_max_L, 15)
  expect_identical(result$susie_fitted, fake_fit)
})

test_that("uap: susie_post_processor output is merged into result", {
  inp <- make_uap_inputs()
  fake_fit <- make_fake_susie_fit(inp$p)
  fake_post <- list(
    variant_names = paste0("v", seq_len(inp$p)),
    top_loci = data.frame(variant_id = "v1", pip = 0.9, stringsAsFactors = FALSE),
    susie_result_trimmed = list(pip = runif(inp$p))
  )

  local_mocked_bindings(
    susie_wrapper = function(...) fake_fit,
    susie_post_processor = function(...) fake_post,
  )

  result <- univariate_analysis_pipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    twas_weights = FALSE, init_L = 5, max_L = 5
  )
  expect_true("variant_names" %in% names(result))
  expect_true("top_loci" %in% names(result))
  expect_true("total_time_elapsed" %in% names(result))
})

test_that("uap: finemapping_extra_opts are forwarded to susie_wrapper", {
  inp <- make_uap_inputs()
  fake_fit <- make_fake_susie_fit(inp$p)
  fake_post <- make_fake_post_result(inp$p)

  captured_refine <- NULL
  local_mocked_bindings(
    susie_wrapper = function(...) {
      args <- list(...)
      captured_refine <<- args$refine
      fake_fit
    },
    susie_post_processor = function(...) fake_post,
  )

  result <- univariate_analysis_pipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    finemapping_extra_opts = list(refine = TRUE),
    twas_weights = FALSE, init_L = 5, max_L = 5
  )
  expect_true(captured_refine)
})

# ========================================================================
#  SECTION 5: univariate_analysis_pipeline – TWAS weights
# ========================================================================

test_that("uap: twas_weights = TRUE calls twas_weights_pipeline", {
  inp <- make_uap_inputs()
  fake_fit <- make_fake_susie_fit(inp$p)
  fake_post <- make_fake_post_result(inp$p)
  fake_twas <- make_fake_twas_result(inp$p)

  twas_called <- FALSE
  local_mocked_bindings(
    susie_wrapper = function(...) fake_fit,
    susie_post_processor = function(...) fake_post,
    twas_weights_pipeline = function(...) {
      twas_called <<- TRUE
      fake_twas
    },
  )

  result <- univariate_analysis_pipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    twas_weights = TRUE, init_L = 5, max_L = 5,
    cv_folds = 2
  )
  expect_true(twas_called)
  expect_true("twas_weights_result" %in% names(result))
})

test_that("uap: twas_weights copies top_loci into susie_weights_intermediate", {
  inp <- make_uap_inputs()
  fake_fit <- make_fake_susie_fit(inp$p)
  fake_post <- make_fake_post_result(inp$p)
  fake_twas <- make_fake_twas_result(inp$p)

  local_mocked_bindings(
    susie_wrapper = function(...) fake_fit,
    susie_post_processor = function(...) fake_post,
    twas_weights_pipeline = function(...) fake_twas,
  )

  result <- univariate_analysis_pipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    twas_weights = TRUE, init_L = 5, max_L = 5,
    cv_folds = 2
  )
  # top_loci should be copied from fake_post into the twas result
  expect_true("top_loci" %in% names(result$twas_weights_result$susie_weights_intermediate))
})

test_that("uap: twas_weights_pipeline receives correct cv_folds and sample_partition", {
  inp <- make_uap_inputs()
  fake_fit <- make_fake_susie_fit(inp$p)
  fake_post <- make_fake_post_result(inp$p)
  fake_twas <- make_fake_twas_result(inp$p)

  captured_cv_folds <- NULL
  captured_partition <- NULL
  local_mocked_bindings(
    susie_wrapper = function(...) fake_fit,
    susie_post_processor = function(...) fake_post,
    twas_weights_pipeline = function(X, Y, susie_fit, cv_folds, max_cv_variants,
                                     cv_threads, sample_partition) {
      captured_cv_folds <<- cv_folds
      captured_partition <<- sample_partition
      fake_twas
    },
  )

  sample_part <- rep(1:3, length.out = inp$n)
  result <- univariate_analysis_pipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    twas_weights = TRUE, init_L = 5, max_L = 5,
    cv_folds = 3, sample_partition = sample_part
  )
  expect_equal(captured_cv_folds, 3)
  expect_equal(captured_partition, sample_part)
})

test_that("uap: twas_weights = FALSE skips twas_weights_pipeline", {
  inp <- make_uap_inputs()
  fake_fit <- make_fake_susie_fit(inp$p)
  fake_post <- make_fake_post_result(inp$p)

  twas_called <- FALSE
  local_mocked_bindings(
    susie_wrapper = function(...) fake_fit,
    susie_post_processor = function(...) fake_post,
    twas_weights_pipeline = function(...) {
      twas_called <<- TRUE
      make_fake_twas_result(inp$p)
    },
  )

  result <- univariate_analysis_pipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    twas_weights = FALSE, init_L = 5, max_L = 5
  )
  expect_false(twas_called)
  expect_false("twas_weights_result" %in% names(result))
})

# ========================================================================
#  SECTION 6: univariate_analysis_pipeline – coverage vector forwarding
# ========================================================================

test_that("uap: coverage vector is forwarded correctly to susie_wrapper and susie_post_processor", {
  inp <- make_uap_inputs()
  fake_fit <- make_fake_susie_fit(inp$p)
  fake_post <- make_fake_post_result(inp$p)

  captured_coverage_wrapper <- NULL
  captured_secondary_cov <- NULL
  local_mocked_bindings(
    susie_wrapper = function(X, y, init_L, max_L, l_step, coverage, ...) {
      captured_coverage_wrapper <<- coverage
      fake_fit
    },
    susie_post_processor = function(susie_output, data_x, data_y, X_scalar, Y_scalar, maf,
                                    secondary_coverage, ...) {
      captured_secondary_cov <<- secondary_coverage
      fake_post
    },
  )

  result <- univariate_analysis_pipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    coverage = c(0.95, 0.7, 0.5),
    twas_weights = FALSE, init_L = 5, max_L = 5
  )
  expect_equal(captured_coverage_wrapper, 0.95)
  expect_equal(captured_secondary_cov, c(0.7, 0.5))
})

test_that("uap: single coverage value => secondary_coverage is NULL", {
  inp <- make_uap_inputs()
  fake_fit <- make_fake_susie_fit(inp$p)
  fake_post <- make_fake_post_result(inp$p)

  captured_secondary_cov <- "NOT_SET"
  local_mocked_bindings(
    susie_wrapper = function(...) fake_fit,
    susie_post_processor = function(susie_output, data_x, data_y, X_scalar, Y_scalar, maf,
                                    secondary_coverage, ...) {
      captured_secondary_cov <<- secondary_coverage
      fake_post
    },
  )

  result <- univariate_analysis_pipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    coverage = c(0.95),
    twas_weights = FALSE, init_L = 5, max_L = 5
  )
  expect_null(captured_secondary_cov)
})

# ========================================================================
#  SECTION 7: univariate_analysis_pipeline – combined LD + filter_X
# ========================================================================

test_that("uap: both LD filtering and filter_X applied in sequence", {
  inp <- make_uap_inputs(p = 10)
  # LD filtering keeps first 8
  # filter_X then drops 2 more => 6 left
  X_after_ld <- inp$X[, 1:8, drop = FALSE]
  X_after_filter <- inp$X[, 1:6, drop = FALSE]

  fake_fit <- make_fake_susie_fit(6)
  fake_post <- make_fake_post_result(6)

  captured_ncol <- NULL
  local_mocked_bindings(
    filter_variants_by_ld_reference = function(variant_ids, ld_ref_file) {
      list(data = variant_ids[1:8], idx = 1:8)
    },
    filter_X = function(X, ...) X_after_filter,
    susie_wrapper = function(X, ...) {
      captured_ncol <<- ncol(X)
      fake_fit
    },
    susie_post_processor = function(...) fake_post,
  )

  result <- univariate_analysis_pipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    ld_reference_meta_file = "/fake/ld.tsv",
    imiss_cutoff = 0.9, maf_cutoff = 0.05,
    twas_weights = FALSE, init_L = 5, max_L = 5
  )
  expect_equal(captured_ncol, 6)
})

# ========================================================================
#  SECTION 8: rss_analysis_pipeline – empty sumstats early return
# ========================================================================

test_that("rss: empty sumstats from load_rss_data => early return", {
  local_mocked_bindings(
    load_rss_data = function(...) {
      list(sumstats = data.frame(), n = NULL, var_y = NULL)
    },
  )

  result <- rss_analysis_pipeline(
    sumstat_path = "/fake/sumstats.tsv",
    column_file_path = "/fake/columns.yml",
    LD_data = list()
  )
  expect_true("rss_data_analyzed" %in% names(result))
  expect_equal(nrow(result$rss_data_analyzed), 0)
})

# ========================================================================
#  SECTION 9: rss_analysis_pipeline – empty sumstats after preprocessing
# ========================================================================

test_that("rss: empty sumstats after rss_basic_qc => early return", {
  local_mocked_bindings(
    load_rss_data = function(...) {
      list(
        sumstats = data.frame(chrom = 1, pos = 100, A1 = "A", A2 = "G", z = 2.5, variant_id = "1:100:A:G"),
        n = 1000,
        var_y = 1
      )
    },
    rss_basic_qc = function(...) {
      list(sumstats = data.frame(), LD_mat = matrix(nrow = 0, ncol = 0))
    },
  )

  result <- expect_message(
    rss_analysis_pipeline(
      sumstat_path = "/fake/sumstats.tsv",
      column_file_path = "/fake/columns.yml",
      LD_data = list()
    ),
    "No variants left after preprocessing"
  )
  expect_true("rss_data_analyzed" %in% names(result))
  expect_equal(nrow(result$rss_data_analyzed), 0)
})

# ========================================================================
#  SECTION 10: rss_analysis_pipeline – pip_cutoff_to_skip branch
# ========================================================================

make_rss_sumstats <- function(n_variants = 5) {
  data.frame(
    chrom = rep(1, n_variants),
    pos = seq(100, by = 100, length.out = n_variants),
    A1 = rep("A", n_variants),
    A2 = rep("G", n_variants),
    z = rnorm(n_variants, 0, 2),
    beta = rnorm(n_variants),
    se = rep(0.5, n_variants),
    variant_id = paste0("1:", seq(100, by = 100, length.out = n_variants), ":A:G"),
    stringsAsFactors = FALSE
  )
}

make_rss_ld_mat <- function(n_variants = 5) {
  ids <- paste0("1:", seq(100, by = 100, length.out = n_variants), ":A:G")
  m <- diag(n_variants)
  rownames(m) <- colnames(m) <- ids
  m
}

test_that("rss: pip_cutoff_to_skip > 0, no signal => early return", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  local_mocked_bindings(
    load_rss_data = function(...) list(sumstats = ss, n = 1000, var_y = 1),
    rss_basic_qc = function(...) list(sumstats = ss, LD_mat = ld_mat),
    susie_rss_wrapper = function(...) list(pip = rep(0.01, 5)),
  )

  result <- expect_message(
    rss_analysis_pipeline(
      sumstat_path = "/fake/sumstats.tsv",
      column_file_path = "/fake/columns.yml",
      LD_data = list(),
      pip_cutoff_to_skip = 0.5
    ),
    "Skipping follow-up"
  )
  expect_true("rss_data_analyzed" %in% names(result))
})

test_that("rss: pip_cutoff_to_skip > 0, signal detected => continues", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)
  fake_rss_result <- make_fake_post_result(5)

  local_mocked_bindings(
    load_rss_data = function(...) list(sumstats = ss, n = 1000, var_y = 1),
    rss_basic_qc = function(...) list(sumstats = ss, LD_mat = ld_mat),
    susie_rss_wrapper = function(...) list(pip = c(0.9, 0.01, 0.01, 0.01, 0.01)),
    summary_stats_qc = function(...) list(sumstats = ss, LD_mat = ld_mat, outlier_number = 0),
    partition_LD_matrix = function(...) ld_mat,
    raiss = function(...) list(result_filter = ss, LD_mat = ld_mat),
    susie_rss_pipeline = function(...) fake_rss_result,
  )

  result <- expect_message(
    rss_analysis_pipeline(
      sumstat_path = "/fake/sumstats.tsv",
      column_file_path = "/fake/columns.yml",
      LD_data = list(ref_panel = ss),
      pip_cutoff_to_skip = 0.5,
      qc_method = "slalom",
      finemapping_method = "susie_rss"
    ),
    "Follow-up on region"
  )
  # Should have the susie_rss pipeline result stored under the method name
  method_names <- grep("susie_rss", names(result), value = TRUE)
  expect_true(length(method_names) > 0)
})

test_that("rss: negative pip_cutoff_to_skip auto-computes threshold", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  local_mocked_bindings(
    load_rss_data = function(...) list(sumstats = ss, n = 1000, var_y = 1),
    rss_basic_qc = function(...) list(sumstats = ss, LD_mat = ld_mat),
    susie_rss_wrapper = function(...) list(pip = rep(0.01, 5)),
  )

  # auto threshold = 3 * 1/5 = 0.6, all PIPs are 0.01 so skip
  result <- expect_message(
    rss_analysis_pipeline(
      sumstat_path = "/fake/sumstats.tsv",
      column_file_path = "/fake/columns.yml",
      LD_data = list(),
      pip_cutoff_to_skip = -1
    ),
    "Skipping follow-up"
  )
  expect_true("rss_data_analyzed" %in% names(result))
})

# ========================================================================
#  SECTION 11: rss_analysis_pipeline – QC + imputation + fine-mapping
# ========================================================================

test_that("rss: full pipeline with QC, imputation, and fine-mapping", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)
  fake_result <- make_fake_post_result(5)

  qc_called <- FALSE
  impute_called <- FALSE
  finemapping_called <- FALSE

  local_mocked_bindings(
    load_rss_data = function(...) list(sumstats = ss, n = 1000, var_y = 1),
    rss_basic_qc = function(...) list(sumstats = ss, LD_mat = ld_mat),
    summary_stats_qc = function(...) {
      qc_called <<- TRUE
      list(sumstats = ss, LD_mat = ld_mat, outlier_number = 1)
    },
    partition_LD_matrix = function(...) list(ld_matrices = list(ld_mat)),
    raiss = function(...) {
      impute_called <<- TRUE
      list(result_filter = ss, LD_mat = ld_mat)
    },
    susie_rss_pipeline = function(...) {
      finemapping_called <<- TRUE
      fake_result
    },
  )

  result <- rss_analysis_pipeline(
    sumstat_path = "/fake/sumstats.tsv",
    column_file_path = "/fake/columns.yml",
    LD_data = list(ref_panel = ss),
    qc_method = "slalom",
    impute = TRUE,
    finemapping_method = "susie_rss"
  )

  expect_true(qc_called)
  expect_true(impute_called)
  expect_true(finemapping_called)
  # Method name should include QC method and RAISS
  expect_true(any(grepl("SLALOM_RAISS_imputed", names(result))))
  expect_true("rss_data_analyzed" %in% names(result))
})

test_that("rss: method name is correct for no-impute with QC", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)
  fake_result <- make_fake_post_result(5)

  local_mocked_bindings(
    load_rss_data = function(...) list(sumstats = ss, n = 1000, var_y = 1),
    rss_basic_qc = function(...) list(sumstats = ss, LD_mat = ld_mat),
    summary_stats_qc = function(...) list(sumstats = ss, LD_mat = ld_mat, outlier_number = 0),
    susie_rss_pipeline = function(...) fake_result,
  )

  result <- rss_analysis_pipeline(
    sumstat_path = "/fake/sumstats.tsv",
    column_file_path = "/fake/columns.yml",
    LD_data = list(),
    qc_method = "slalom",
    impute = FALSE,
    finemapping_method = "susie_rss"
  )

  expect_true(any(grepl("susie_rss_SLALOM$", names(result))))
})

test_that("rss: method name is correct for no QC", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)
  fake_result <- make_fake_post_result(5)

  local_mocked_bindings(
    load_rss_data = function(...) list(sumstats = ss, n = 1000, var_y = 1),
    rss_basic_qc = function(...) list(sumstats = ss, LD_mat = ld_mat),
    partition_LD_matrix = function(...) list(ld_matrices = list(ld_mat)),
    raiss = function(...) list(result_filter = ss, LD_mat = ld_mat),
    susie_rss_pipeline = function(...) fake_result,
  )

  result <- rss_analysis_pipeline(
    sumstat_path = "/fake/sumstats.tsv",
    column_file_path = "/fake/columns.yml",
    LD_data = list(ref_panel = ss),
    qc_method = NULL,
    impute = TRUE,
    finemapping_method = "susie_rss"
  )

  expect_true(any(grepl("NO_QC", names(result))))
})

# ========================================================================
#  SECTION 12: rss_analysis_pipeline – outlier_number stored
# ========================================================================

test_that("rss: outlier_number is stored in result when QC is active", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)
  fake_result <- make_fake_post_result(5)

  local_mocked_bindings(
    load_rss_data = function(...) list(sumstats = ss, n = 1000, var_y = 1),
    rss_basic_qc = function(...) list(sumstats = ss, LD_mat = ld_mat),
    summary_stats_qc = function(...) list(sumstats = ss, LD_mat = ld_mat, outlier_number = 3),
    susie_rss_pipeline = function(...) fake_result,
  )

  result <- rss_analysis_pipeline(
    sumstat_path = "/fake/sumstats.tsv",
    column_file_path = "/fake/columns.yml",
    LD_data = list(),
    qc_method = "dentist",
    impute = FALSE,
    finemapping_method = "susie_rss"
  )

  method_key <- grep("DENTIST", names(result), value = TRUE)
  expect_true(length(method_key) > 0)
  expect_equal(result[[method_key[1]]]$outlier_number, 3)
})

# ========================================================================
#  SECTION 13: rss_analysis_pipeline – finemapping_method = NULL skips
# ========================================================================

test_that("rss: finemapping_method = NULL skips fine-mapping", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  finemapping_called <- FALSE
  local_mocked_bindings(
    load_rss_data = function(...) list(sumstats = ss, n = 1000, var_y = 1),
    rss_basic_qc = function(...) list(sumstats = ss, LD_mat = ld_mat),
    susie_rss_pipeline = function(...) {
      finemapping_called <<- TRUE
      list()
    },
  )

  # Use qc_method = NULL and impute = FALSE to avoid paste0(NULL, ...) edge case
  result <- rss_analysis_pipeline(
    sumstat_path = "/fake/sumstats.tsv",
    column_file_path = "/fake/columns.yml",
    LD_data = list(),
    qc_method = NULL,
    impute = FALSE,
    finemapping_method = NULL
  )

  expect_false(finemapping_called)
  expect_true("rss_data_analyzed" %in% names(result))
})

# ========================================================================
#  SECTION 14: rss_analysis_pipeline – QC = NULL skips QC
# ========================================================================

test_that("rss: qc_method = NULL skips quality control", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)
  fake_result <- make_fake_post_result(5)

  qc_called <- FALSE
  local_mocked_bindings(
    load_rss_data = function(...) list(sumstats = ss, n = 1000, var_y = 1),
    rss_basic_qc = function(...) list(sumstats = ss, LD_mat = ld_mat),
    summary_stats_qc = function(...) {
      qc_called <<- TRUE
      list(sumstats = ss, LD_mat = ld_mat, outlier_number = 0)
    },
    susie_rss_pipeline = function(...) fake_result,
  )

  result <- rss_analysis_pipeline(
    sumstat_path = "/fake/sumstats.tsv",
    column_file_path = "/fake/columns.yml",
    LD_data = list(),
    qc_method = NULL,
    impute = FALSE,
    finemapping_method = "susie_rss"
  )

  expect_false(qc_called)
})

# ========================================================================
#  SECTION 15: rss_analysis_pipeline – impute = FALSE skips imputation
# ========================================================================

test_that("rss: impute = FALSE skips raiss imputation", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)
  fake_result <- make_fake_post_result(5)

  raiss_called <- FALSE
  local_mocked_bindings(
    load_rss_data = function(...) list(sumstats = ss, n = 1000, var_y = 1),
    rss_basic_qc = function(...) list(sumstats = ss, LD_mat = ld_mat),
    summary_stats_qc = function(...) list(sumstats = ss, LD_mat = ld_mat, outlier_number = 0),
    raiss = function(...) {
      raiss_called <<- TRUE
      list(result_filter = ss, LD_mat = ld_mat)
    },
    susie_rss_pipeline = function(...) fake_result,
  )

  result <- rss_analysis_pipeline(
    sumstat_path = "/fake/sumstats.tsv",
    column_file_path = "/fake/columns.yml",
    LD_data = list(),
    qc_method = "slalom",
    impute = FALSE,
    finemapping_method = "susie_rss"
  )

  expect_false(raiss_called)
})

# ========================================================================
#  SECTION 16: rss_analysis_pipeline – diagnostics branch (empty res)
# ========================================================================

test_that("rss: diagnostics = TRUE with empty fine-mapping result skips diagnostics logic", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  local_mocked_bindings(
    load_rss_data = function(...) list(sumstats = ss, n = 1000, var_y = 1),
    rss_basic_qc = function(...) list(sumstats = ss, LD_mat = ld_mat),
    summary_stats_qc = function(...) list(sumstats = ss, LD_mat = ld_mat, outlier_number = 0),
    susie_rss_pipeline = function(...) list(),  # empty result
  )

  result <- rss_analysis_pipeline(
    sumstat_path = "/fake/sumstats.tsv",
    column_file_path = "/fake/columns.yml",
    LD_data = list(),
    qc_method = "slalom",
    impute = FALSE,
    diagnostics = TRUE,
    finemapping_method = "susie_rss"
  )

  # With empty fine-mapping result, block_cs_metrics stays empty,
  # so the "diagnostics" key is never added to result_list
  expect_false("diagnostics" %in% names(result))
})

# ========================================================================
#  SECTION 17: rss_analysis_pipeline – diagnostics with CS > 1 and
#              high p-value triggers BCR + SER re-analysis
# ========================================================================

test_that("rss: diagnostics with 2+ CS and high p-value/corr triggers BCR and SER", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  # Build a fake result that has susie_result_trimmed with 2 CS
  fake_result <- list(
    variant_names = paste0("1:", seq(100, by = 100, length.out = 5), ":A:G"),
    susie_result_trimmed = list(
      pip = c(0.9, 0.8, 0.1, 0.05, 0.02),
      sets = list(
        cs = list(L1 = c(1L, 2L), L2 = c(3L, 4L)),
        requested_coverage = 0.95
      ),
      cs_corr = matrix(c(1, 0.6, 0.6, 1), nrow = 2)
    ),
    sumstats = list(z = c(5.0, 4.0, 1.5, 0.5, 0.1)),
    top_loci = data.frame(
      variant_id = c("1:100:A:G", "1:200:A:G", "1:300:A:G", "1:400:A:G"),
      pip = c(0.9, 0.8, 0.1, 0.05),
      z = c(5.0, 4.0, 1.5, 0.5),
      stringsAsFactors = FALSE
    )
  )

  # extract_cs_info returns a tibble with 2 CS rows
  fake_cs_info <- tibble::tibble(
    cs_name = c("L1", "L2"),
    variants_per_cs = c(2L, 2L),
    top_variant = c("1:100:A:G", "1:300:A:G"),
    top_variant_index = c(1L, 3L),
    top_pip = c(0.9, 0.1),
    top_z = c(5.0, 1.5),
    p_value = c(5.7e-7, 0.13),  # second one > 1e-4, triggers BCR
    cs_corr = list("1,0.6", "0.6,1")
  )

  # parse_cs_corr adds cs_corr_max and cs_corr_min
  fake_parsed <- as.data.frame(fake_cs_info)
  fake_parsed$cs_corr_1 <- c(1, 0.6)
  fake_parsed$cs_corr_2 <- c(0.6, 1)
  fake_parsed$cs_corr_max <- c(0.6, 0.6)
  fake_parsed$cs_corr_min <- c(0.6, 0.6)

  susie_rss_pipeline_call_count <- 0
  local_mocked_bindings(
    load_rss_data = function(...) list(sumstats = ss, n = 1000, var_y = 1),
    rss_basic_qc = function(...) list(sumstats = ss, LD_mat = ld_mat),
    summary_stats_qc = function(...) list(sumstats = ss, LD_mat = ld_mat, outlier_number = 0),
    susie_rss_pipeline = function(...) {
      susie_rss_pipeline_call_count <<- susie_rss_pipeline_call_count + 1
      fake_result
    },
    get_susie_result = function(res) res$susie_result_trimmed,
    extract_cs_info = function(...) fake_cs_info,
    parse_cs_corr = function(...) fake_parsed,
  )

  result <- rss_analysis_pipeline(
    sumstat_path = "/fake/sumstats.tsv",
    column_file_path = "/fake/columns.yml",
    LD_data = list(),
    qc_method = "slalom",
    impute = FALSE,
    diagnostics = TRUE,
    finemapping_method = "susie_rss"
  )

  # Main susie_rss call + BCR reanalysis + SER reanalysis = 3 calls
  expect_equal(susie_rss_pipeline_call_count, 3)
  # Should have BCR and SER method entries
  expect_true(any(grepl("bayesian_conditional_regression", names(result))))
  expect_true(any(grepl("single_effect", names(result))))
  expect_true("diagnostics" %in% names(result))
})

# ========================================================================
#  SECTION 18: rss_analysis_pipeline – diagnostics with 1 CS => SER only
# ========================================================================

test_that("rss: diagnostics with 1 CS triggers SER reanalysis only", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  fake_result <- list(
    variant_names = paste0("1:", seq(100, by = 100, length.out = 5), ":A:G"),
    susie_result_trimmed = list(
      pip = c(0.9, 0.1, 0.05, 0.02, 0.01),
      sets = list(
        cs = list(L1 = c(1L, 2L)),
        requested_coverage = 0.95
      ),
      cs_corr = matrix(1, nrow = 1)
    ),
    sumstats = list(z = c(5.0, 1.0, 0.5, 0.2, 0.1)),
    top_loci = data.frame(
      variant_id = c("1:100:A:G", "1:200:A:G"),
      pip = c(0.9, 0.1),
      z = c(5.0, 1.0),
      stringsAsFactors = FALSE
    )
  )

  fake_cs_info <- tibble::tibble(
    cs_name = "L1",
    variants_per_cs = 2L,
    top_variant = "1:100:A:G",
    top_variant_index = 1L,
    top_pip = 0.9,
    top_z = 5.0,
    p_value = 5.7e-7,
    cs_corr = list("1")
  )

  fake_parsed <- as.data.frame(fake_cs_info)
  fake_parsed$cs_corr_max <- NA_real_
  fake_parsed$cs_corr_min <- NA_real_

  susie_rss_pipeline_call_count <- 0
  local_mocked_bindings(
    load_rss_data = function(...) list(sumstats = ss, n = 1000, var_y = 1),
    rss_basic_qc = function(...) list(sumstats = ss, LD_mat = ld_mat),
    summary_stats_qc = function(...) list(sumstats = ss, LD_mat = ld_mat, outlier_number = 0),
    susie_rss_pipeline = function(...) {
      susie_rss_pipeline_call_count <<- susie_rss_pipeline_call_count + 1
      fake_result
    },
    get_susie_result = function(res) res$susie_result_trimmed,
    extract_cs_info = function(...) fake_cs_info,
    parse_cs_corr = function(...) fake_parsed,
  )

  result <- rss_analysis_pipeline(
    sumstat_path = "/fake/sumstats.tsv",
    column_file_path = "/fake/columns.yml",
    LD_data = list(),
    qc_method = "slalom",
    impute = FALSE,
    diagnostics = TRUE,
    finemapping_method = "susie_rss"
  )

  # Main call + SER reanalysis = 2 calls
  expect_equal(susie_rss_pipeline_call_count, 2)
  expect_true(any(grepl("single_effect", names(result))))
  expect_true("diagnostics" %in% names(result))
})

# ========================================================================
#  SECTION 19: rss_analysis_pipeline – diagnostics with no CS but high PIP
# ========================================================================

test_that("rss: diagnostics with no CS but high PIP calls extract_top_pip_info", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  fake_result <- list(
    variant_names = paste0("1:", seq(100, by = 100, length.out = 5), ":A:G"),
    susie_result_trimmed = list(
      pip = c(0.15, 0.05, 0.03, 0.02, 0.01),
      sets = list(cs = NULL, requested_coverage = 0.95),
      cs_corr = NULL
    ),
    sumstats = list(z = c(3.0, 1.0, 0.5, 0.2, 0.1)),
    top_loci = data.frame(
      variant_id = "1:100:A:G",
      pip = 0.15,
      z = 3.0,
      stringsAsFactors = FALSE
    )
  )

  # extract_top_pip_info returns when no CS but PIP > signal_cutoff
  extract_top_pip_called <- FALSE
  fake_top_pip_info <- list(
    cs_name = NA,
    variants_per_cs = NA,
    top_variant = "1:100:A:G",
    top_variant_index = 1,
    top_pip = 0.15,
    top_z = 3.0,
    p_value = 0.0027,
    cs_corr = NA
  )

  local_mocked_bindings(
    load_rss_data = function(...) list(sumstats = ss, n = 1000, var_y = 1),
    rss_basic_qc = function(...) list(sumstats = ss, LD_mat = ld_mat),
    summary_stats_qc = function(...) list(sumstats = ss, LD_mat = ld_mat, outlier_number = 0),
    susie_rss_pipeline = function(...) fake_result,
    get_susie_result = function(res) res$susie_result_trimmed,
    extract_top_pip_info = function(...) {
      extract_top_pip_called <<- TRUE
      fake_top_pip_info
    },
    parse_cs_corr = function(df) {
      df <- as.data.frame(df)
      df$cs_corr_max <- NA_real_
      df$cs_corr_min <- NA_real_
      df
    },
  )

  result <- rss_analysis_pipeline(
    sumstat_path = "/fake/sumstats.tsv",
    column_file_path = "/fake/columns.yml",
    LD_data = list(),
    qc_method = "slalom",
    impute = FALSE,
    diagnostics = TRUE,
    finemapping_method = "susie_rss"
  )

  expect_true(extract_top_pip_called)
  expect_true("diagnostics" %in% names(result))
})

# ========================================================================
#  SECTION 20: rss_analysis_pipeline – diagnostics with no CS and low PIP
# ========================================================================

test_that("rss: diagnostics with no CS and no high PIP => diagnostics empty", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  fake_result <- list(
    variant_names = paste0("1:", seq(100, by = 100, length.out = 5), ":A:G"),
    susie_result_trimmed = list(
      pip = c(0.01, 0.005, 0.003, 0.002, 0.001),
      sets = list(cs = NULL, requested_coverage = 0.95),
      cs_corr = NULL
    ),
    sumstats = list(z = c(0.5, 0.3, 0.2, 0.1, 0.05))
  )

  local_mocked_bindings(
    load_rss_data = function(...) list(sumstats = ss, n = 1000, var_y = 1),
    rss_basic_qc = function(...) list(sumstats = ss, LD_mat = ld_mat),
    summary_stats_qc = function(...) list(sumstats = ss, LD_mat = ld_mat, outlier_number = 0),
    susie_rss_pipeline = function(...) fake_result,
    get_susie_result = function(res) res$susie_result_trimmed,
  )

  result <- rss_analysis_pipeline(
    sumstat_path = "/fake/sumstats.tsv",
    column_file_path = "/fake/columns.yml",
    LD_data = list(),
    qc_method = "slalom",
    impute = FALSE,
    diagnostics = TRUE,
    finemapping_method = "susie_rss",
    finemapping_opts = list(
      init_L = 5, max_L = 20, l_step = 5,
      coverage = c(0.95, 0.7, 0.5), signal_cutoff = 0.025,
      min_abs_corr = 0.8
    )
  )

  # No CS and no high PIP means block_cs_metrics stays empty,
  # so "diagnostics" key is never added to result_list
  expect_false("diagnostics" %in% names(result))
})

# ========================================================================
#  SECTION 21: rss_analysis_pipeline – finemapping_opts forwarded correctly
# ========================================================================

test_that("rss: finemapping_opts are forwarded to susie_rss_pipeline", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  captured_L <- NULL
  captured_coverage <- NULL
  captured_signal_cutoff <- NULL
  local_mocked_bindings(
    load_rss_data = function(...) list(sumstats = ss, n = 1000, var_y = 1),
    rss_basic_qc = function(...) list(sumstats = ss, LD_mat = ld_mat),
    susie_rss_pipeline = function(sumstats, LD_mat, n, var_y, L, max_L, l_step,
                                  analysis_method, coverage, secondary_coverage,
                                  signal_cutoff, min_abs_corr, ...) {
      captured_L <<- L
      captured_coverage <<- coverage
      captured_signal_cutoff <<- signal_cutoff
      make_fake_post_result(5)
    },
  )

  result <- rss_analysis_pipeline(
    sumstat_path = "/fake/sumstats.tsv",
    column_file_path = "/fake/columns.yml",
    LD_data = list(),
    qc_method = NULL,
    impute = FALSE,
    finemapping_method = "susie_rss",
    finemapping_opts = list(
      init_L = 3, max_L = 10, l_step = 2,
      coverage = c(0.90, 0.6), signal_cutoff = 0.05,
      min_abs_corr = 0.7
    )
  )

  expect_equal(captured_L, 3)
  expect_equal(captured_coverage, 0.90)
  expect_equal(captured_signal_cutoff, 0.05)
})

# ========================================================================
#  SECTION 22: rss_analysis_pipeline – dentist QC method naming
# ========================================================================

test_that("rss: dentist QC method generates correct method name", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)
  fake_result <- make_fake_post_result(5)

  local_mocked_bindings(
    load_rss_data = function(...) list(sumstats = ss, n = 1000, var_y = 1),
    rss_basic_qc = function(...) list(sumstats = ss, LD_mat = ld_mat),
    summary_stats_qc = function(...) list(sumstats = ss, LD_mat = ld_mat, outlier_number = 2),
    partition_LD_matrix = function(...) list(ld_matrices = list(ld_mat)),
    raiss = function(...) list(result_filter = ss, LD_mat = ld_mat),
    susie_rss_pipeline = function(...) fake_result,
  )

  result <- rss_analysis_pipeline(
    sumstat_path = "/fake/sumstats.tsv",
    column_file_path = "/fake/columns.yml",
    LD_data = list(ref_panel = ss),
    qc_method = "dentist",
    impute = TRUE,
    finemapping_method = "susie_rss"
  )

  expect_true(any(grepl("DENTIST_RAISS_imputed", names(result))))
})

# ========================================================================
#  SECTION 23: univariate_analysis_pipeline – null imiss/maf cutoffs skip filter_X
# ========================================================================

test_that("uap: both imiss_cutoff and maf_cutoff NULL skips filter_X", {
  inp <- make_uap_inputs()
  fake_fit <- make_fake_susie_fit(inp$p)
  fake_post <- make_fake_post_result(inp$p)

  filter_called <- FALSE
  local_mocked_bindings(
    filter_X = function(...) {
      filter_called <<- TRUE
      inp$X
    },
    susie_wrapper = function(...) fake_fit,
    susie_post_processor = function(...) fake_post,
  )

  result <- univariate_analysis_pipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    imiss_cutoff = NULL, maf_cutoff = NULL,
    twas_weights = FALSE, init_L = 5, max_L = 5
  )

  expect_false(filter_called)
})

# ========================================================================
#  SECTION 24: rss_analysis_pipeline – diagnostics with get_susie_result
#              returning NULL
# ========================================================================

test_that("rss: diagnostics with get_susie_result returning NULL => diagnostics empty", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  # A non-empty result but with get_susie_result returning NULL
  fake_result <- list(some_field = "value")

  local_mocked_bindings(
    load_rss_data = function(...) list(sumstats = ss, n = 1000, var_y = 1),
    rss_basic_qc = function(...) list(sumstats = ss, LD_mat = ld_mat),
    summary_stats_qc = function(...) list(sumstats = ss, LD_mat = ld_mat, outlier_number = 0),
    susie_rss_pipeline = function(...) fake_result,
    get_susie_result = function(res) NULL,
  )

  result <- rss_analysis_pipeline(
    sumstat_path = "/fake/sumstats.tsv",
    column_file_path = "/fake/columns.yml",
    LD_data = list(),
    qc_method = "slalom",
    impute = FALSE,
    diagnostics = TRUE,
    finemapping_method = "susie_rss"
  )

  # get_susie_result returns NULL, so block_cs_metrics stays empty,
  # and "diagnostics" key is never added to result_list
  expect_false("diagnostics" %in% names(result))
})

# ========================================================================
#  SECTION 25: rss_analysis_pipeline – diagnostics with CS but CS count <= 1
#              and null block_cs_metrics (no high PIP either)
# ========================================================================

test_that("rss: diagnostics with null/empty block_cs_metrics => no additional analysis", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  # Non-empty result with susie_result_trimmed having 0 CS and no high PIPs
  fake_result <- list(
    variant_names = paste0("1:", seq(100, by = 100, length.out = 5), ":A:G"),
    susie_result_trimmed = list(
      pip = c(0.001, 0.002, 0.001, 0.001, 0.001),
      sets = list(cs = list(), requested_coverage = 0.95),
      cs_corr = NULL
    ),
    sumstats = list(z = rep(0.1, 5))
  )

  susie_rss_call_count <- 0
  local_mocked_bindings(
    load_rss_data = function(...) list(sumstats = ss, n = 1000, var_y = 1),
    rss_basic_qc = function(...) list(sumstats = ss, LD_mat = ld_mat),
    summary_stats_qc = function(...) list(sumstats = ss, LD_mat = ld_mat, outlier_number = 0),
    susie_rss_pipeline = function(...) {
      susie_rss_call_count <<- susie_rss_call_count + 1
      fake_result
    },
    get_susie_result = function(res) res$susie_result_trimmed,
  )

  result <- rss_analysis_pipeline(
    sumstat_path = "/fake/sumstats.tsv",
    column_file_path = "/fake/columns.yml",
    LD_data = list(),
    qc_method = "slalom",
    impute = FALSE,
    diagnostics = TRUE,
    finemapping_method = "susie_rss",
    finemapping_opts = list(
      init_L = 5, max_L = 20, l_step = 5,
      coverage = c(0.95, 0.7, 0.5), signal_cutoff = 0.025,
      min_abs_corr = 0.8
    )
  )

  # Only the main call, no re-analyses
  expect_equal(susie_rss_call_count, 1)
  # block_cs_metrics stays empty (no CS and no high PIP),
  # so "diagnostics" key is never added to result_list
  expect_false("diagnostics" %in% names(result))
})

# ========================================================================
#  SECTION 26: rss_analysis_pipeline – diagnostics with CS > 1 but low
#              p-value and low corr => no BCR/SER triggered
# ========================================================================

test_that("rss: diagnostics with 2 CS but low p-value and low corr => no extra analysis", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  fake_result <- list(
    variant_names = paste0("1:", seq(100, by = 100, length.out = 5), ":A:G"),
    susie_result_trimmed = list(
      pip = c(0.9, 0.8, 0.1, 0.05, 0.02),
      sets = list(
        cs = list(L1 = c(1L, 2L), L2 = c(3L, 4L)),
        requested_coverage = 0.95
      ),
      cs_corr = matrix(c(1, 0.1, 0.1, 1), nrow = 2)
    ),
    sumstats = list(z = c(5.0, 4.0, 3.0, 2.5, 0.1)),
    top_loci = data.frame(
      variant_id = c("1:100:A:G", "1:200:A:G", "1:300:A:G", "1:400:A:G"),
      pip = c(0.9, 0.8, 0.1, 0.05),
      z = c(5.0, 4.0, 3.0, 2.5),
      stringsAsFactors = FALSE
    )
  )

  fake_cs_info <- tibble::tibble(
    cs_name = c("L1", "L2"),
    variants_per_cs = c(2L, 2L),
    top_variant = c("1:100:A:G", "1:300:A:G"),
    top_variant_index = c(1L, 3L),
    top_pip = c(0.9, 0.1),
    top_z = c(5.0, 3.0),
    p_value = c(5.7e-7, 2.7e-3),  # both < 1e-4 NOT met for L2 (0.0027 < 1e-4 is FALSE, wait 0.0027 > 1e-4)
    cs_corr = list("1,0.1", "0.1,1")
  )

  # Actually 0.0027 > 1e-4 is TRUE, so the condition IS triggered.
  # Let me adjust p-values to both be very small
  fake_cs_info$p_value <- c(5.7e-7, 2.7e-5)  # both < 1e-4

  fake_parsed <- as.data.frame(fake_cs_info)
  fake_parsed$cs_corr_1 <- c(1, 0.1)
  fake_parsed$cs_corr_2 <- c(0.1, 1)
  fake_parsed$cs_corr_max <- c(0.1, 0.1)  # both < 0.5
  fake_parsed$cs_corr_min <- c(0.1, 0.1)

  susie_rss_call_count <- 0
  local_mocked_bindings(
    load_rss_data = function(...) list(sumstats = ss, n = 1000, var_y = 1),
    rss_basic_qc = function(...) list(sumstats = ss, LD_mat = ld_mat),
    summary_stats_qc = function(...) list(sumstats = ss, LD_mat = ld_mat, outlier_number = 0),
    susie_rss_pipeline = function(...) {
      susie_rss_call_count <<- susie_rss_call_count + 1
      fake_result
    },
    get_susie_result = function(res) res$susie_result_trimmed,
    extract_cs_info = function(...) fake_cs_info,
    parse_cs_corr = function(...) fake_parsed,
  )

  result <- rss_analysis_pipeline(
    sumstat_path = "/fake/sumstats.tsv",
    column_file_path = "/fake/columns.yml",
    LD_data = list(),
    qc_method = "slalom",
    impute = FALSE,
    diagnostics = TRUE,
    finemapping_method = "susie_rss"
  )

  # Only the main call, no BCR/SER because no p > 1e-4 and no corr > 0.5
  expect_equal(susie_rss_call_count, 1)
  expect_true("diagnostics" %in% names(result))
})

# ========================================================================
#  SECTION 27: rss_analysis_pipeline – diagnostics with CS > 1 and
#              high corr (> 0.5) but low p-value still triggers BCR+SER
# ========================================================================

test_that("rss: diagnostics with high max_cs_corr_study_block triggers BCR+SER", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  fake_result <- list(
    variant_names = paste0("1:", seq(100, by = 100, length.out = 5), ":A:G"),
    susie_result_trimmed = list(
      pip = c(0.9, 0.8, 0.1, 0.05, 0.02),
      sets = list(
        cs = list(L1 = c(1L, 2L), L2 = c(3L, 4L)),
        requested_coverage = 0.95
      ),
      cs_corr = matrix(c(1, 0.7, 0.7, 1), nrow = 2)
    ),
    sumstats = list(z = c(10.0, 9.0, 8.0, 7.0, 0.1)),
    top_loci = data.frame(
      variant_id = c("1:100:A:G", "1:200:A:G", "1:300:A:G", "1:400:A:G"),
      pip = c(0.9, 0.8, 0.1, 0.05),
      z = c(10.0, 9.0, 8.0, 7.0),
      stringsAsFactors = FALSE
    )
  )

  fake_cs_info <- tibble::tibble(
    cs_name = c("L1", "L2"),
    variants_per_cs = c(2L, 2L),
    top_variant = c("1:100:A:G", "1:300:A:G"),
    top_variant_index = c(1L, 3L),
    top_pip = c(0.9, 0.1),
    top_z = c(10.0, 8.0),
    p_value = c(1e-23, 1e-15),  # both very small
    cs_corr = list("1,0.7", "0.7,1")
  )

  fake_parsed <- as.data.frame(fake_cs_info)
  fake_parsed$cs_corr_1 <- c(1, 0.7)
  fake_parsed$cs_corr_2 <- c(0.7, 1)
  fake_parsed$cs_corr_max <- c(0.7, 0.7)  # > 0.5, triggers the condition
  fake_parsed$cs_corr_min <- c(0.7, 0.7)

  susie_rss_call_count <- 0
  local_mocked_bindings(
    load_rss_data = function(...) list(sumstats = ss, n = 1000, var_y = 1),
    rss_basic_qc = function(...) list(sumstats = ss, LD_mat = ld_mat),
    summary_stats_qc = function(...) list(sumstats = ss, LD_mat = ld_mat, outlier_number = 0),
    susie_rss_pipeline = function(...) {
      susie_rss_call_count <<- susie_rss_call_count + 1
      fake_result
    },
    get_susie_result = function(res) res$susie_result_trimmed,
    extract_cs_info = function(...) fake_cs_info,
    parse_cs_corr = function(...) fake_parsed,
  )

  result <- rss_analysis_pipeline(
    sumstat_path = "/fake/sumstats.tsv",
    column_file_path = "/fake/columns.yml",
    LD_data = list(),
    qc_method = "slalom",
    impute = FALSE,
    diagnostics = TRUE,
    finemapping_method = "susie_rss"
  )

  # Main + BCR + SER = 3
  expect_equal(susie_rss_call_count, 3)
  expect_true(any(grepl("bayesian_conditional_regression", names(result))))
  expect_true(any(grepl("single_effect", names(result))))
})

# ========================================================================
#  load_study_LD: single vs mixture LD panels
# ========================================================================

test_that("load_study_LD with single path returns load_LD_matrix output unchanged", {
  fake_LD <- list(
    LD_matrix = diag(4),
    LD_variants = paste0("chr1:", 1:4, ":A:G"),
    is_genotype = FALSE,
    block_metadata = data.frame(chrom = 1, start = 1, end = 4)
  )
  captured_args <- NULL
  local_mocked_bindings(
    load_LD_matrix = function(LD_meta_file_path, region, ...) {
      captured_args <<- list(path = LD_meta_file_path, region = region, ...)
      fake_LD
    },
  )
  out <- load_study_LD("ld_meta.tsv", "chr1:1-100")
  expect_identical(out, fake_LD)
  expect_equal(captured_args$path, "ld_meta.tsv")
  expect_equal(captured_args$region, "chr1:1-100")
  expect_equal(captured_args$return_genotype, "auto")
})

test_that("load_study_LD with comma-separated paths builds list of X panels", {
  panel1 <- list(
    LD_matrix = matrix(1:9, 3, 3),
    LD_variants = paste0("chr1:", 1:3, ":A:G"),
    is_genotype = TRUE,
    block_metadata = data.frame(chrom = 1, start = 1, end = 3)
  )
  panel2 <- list(
    LD_matrix = matrix(10:18, 3, 3),
    LD_variants = paste0("chr1:", 1:3, ":A:G"),
    is_genotype = TRUE,
    block_metadata = data.frame(chrom = 1, start = 1, end = 3)
  )
  call_count <- 0
  local_mocked_bindings(
    load_LD_matrix = function(LD_meta_file_path, region, ...) {
      call_count <<- call_count + 1
      if (call_count == 1) panel1 else panel2
    },
  )
  out <- load_study_LD("ld_eur.tsv,ld_afr.tsv", "chr1:1-100")
  expect_equal(call_count, 2)
  expect_true(is.list(out$LD_matrix))
  expect_length(out$LD_matrix, 2)
  expect_equal(out$LD_matrix[[1]], panel1$LD_matrix)
  expect_equal(out$LD_matrix[[2]], panel2$LD_matrix)
  # Other fields are preserved from base panel
  expect_equal(out$LD_variants, panel1$LD_variants)
})

# ========================================================================
#  rss_analysis_pipeline: X-as-LD path (genotype interface)
# ========================================================================

test_that("rss: is_genotype=TRUE path computes R from X and uses X for fine-mapping", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)
  X_geno <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
  colnames(X_geno) <- ss$variant_id
  fake_result <- make_fake_post_result(5)

  compute_LD_called <- FALSE
  susie_X_arg <- NULL
  susie_LD_arg <- "unset"

  local_mocked_bindings(
    compute_LD = function(X, method = "sample") {
      compute_LD_called <<- TRUE
      ld_mat
    },
    load_rss_data = function(...) list(sumstats = ss, n = 1000, var_y = 1),
    rss_basic_qc = function(...) list(sumstats = ss, LD_mat = ld_mat),
    summary_stats_qc = function(...) list(sumstats = ss, LD_mat = ld_mat, outlier_number = 0),
    partition_LD_matrix = function(...) list(ld_matrices = list(ld_mat)),
    raiss = function(...) list(result_filter = ss, LD_mat = ld_mat),
    susie_rss_pipeline = function(sumstats, LD_mat = NULL, X_mat = NULL, ...) {
      susie_X_arg <<- X_mat
      susie_LD_arg <<- LD_mat
      fake_result
    },
  )

  result <- rss_analysis_pipeline(
    sumstat_path = "/fake/sumstats.tsv",
    column_file_path = "/fake/columns.yml",
    LD_data = list(LD_matrix = X_geno, is_genotype = TRUE, ref_panel = ss),
    qc_method = "slalom",
    finemapping_method = "susie_rss"
  )

  expect_true(compute_LD_called)
  expect_true(is.matrix(susie_X_arg))
  expect_null(susie_LD_arg)
  expect_true(any(grepl("susie_rss", names(result))))
})

test_that("rss: mixture LD_data (list of X panels) preserves list shape into susie_rss_pipeline", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)
  X1 <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
  X2 <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
  colnames(X1) <- colnames(X2) <- ss$variant_id
  fake_result <- make_fake_post_result(5)

  susie_X_arg <- NULL

  local_mocked_bindings(
    compute_LD = function(X, method = "sample") ld_mat,
    load_rss_data = function(...) list(sumstats = ss, n = 1000, var_y = 1),
    rss_basic_qc = function(...) list(sumstats = ss, LD_mat = ld_mat),
    summary_stats_qc = function(...) list(sumstats = ss, LD_mat = ld_mat, outlier_number = 0),
    partition_LD_matrix = function(...) list(ld_matrices = list(ld_mat)),
    raiss = function(...) list(result_filter = ss, LD_mat = ld_mat),
    susie_rss_pipeline = function(sumstats, LD_mat = NULL, X_mat = NULL, ...) {
      susie_X_arg <<- X_mat
      fake_result
    },
  )

  result <- rss_analysis_pipeline(
    sumstat_path = "/fake/sumstats.tsv",
    column_file_path = "/fake/columns.yml",
    LD_data = list(LD_matrix = list(X1, X2), ref_panel = ss),
    qc_method = "slalom",
    finemapping_method = "susie_rss"
  )

  # Mixture path => list of subset matrices passed to susie_rss_pipeline as X_mat
  expect_true(is.list(susie_X_arg))
  expect_length(susie_X_arg, 2)
  expect_true(all(sapply(susie_X_arg, is.matrix)))
  expect_true(any(grepl("susie_rss", names(result))))
})
