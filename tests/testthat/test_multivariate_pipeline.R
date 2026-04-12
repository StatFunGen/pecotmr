context("multivariate_pipeline")

# ===========================================================================
# Helpers: small synthetic multivariate data
# ===========================================================================

make_mv_data <- function(n = 30, p = 10, r = 3, seed = 42) {
  set.seed(seed)
  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n, ncol = p)
  rownames(X) <- paste0("sample_", 1:n)
  colnames(X) <- paste0("chr1:", seq(100, by = 100, length.out = p), ":A:G")
  Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
  rownames(Y) <- rownames(X)
  colnames(Y) <- paste0("cond_", 1:r)
  maf <- colMeans(X) / 2
  list(X = X, Y = Y, maf = maf)
}

# NOTE: In multivariate_analysis_pipeline, the mvsusieR check (requireNamespace)
# runs BEFORE input validation. So when mvsusieR is not installed, the function
# stops with "please install mvsusieR" before ever reaching the X/Y/maf checks.
# Input validation tests therefore require mvsusieR to be installed.

# ===========================================================================
# Input validation tests (require mvsusieR since validation comes after the
# mvsusieR availability check in the pipeline)
# ===========================================================================

test_that("multivariate_analysis_pipeline errors when X is not a numeric matrix", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, input validation unreachable")
  d <- make_mv_data()
  expect_error(
    multivariate_analysis_pipeline(
      X = as.data.frame(d$X), Y = d$Y, maf = d$maf, pip_cutoff_to_skip = 0
    ),
    "X must be a numeric matrix"
  )
})

test_that("multivariate_analysis_pipeline errors when Y is not a numeric matrix", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, input validation unreachable")
  d <- make_mv_data()
  expect_error(
    multivariate_analysis_pipeline(
      X = d$X, Y = as.data.frame(d$Y), maf = d$maf, pip_cutoff_to_skip = 0
    ),
    "Y must be a numeric matrix"
  )
})

test_that("multivariate_analysis_pipeline errors when nrow(X) != nrow(Y)", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, input validation unreachable")
  d <- make_mv_data()
  Y_short <- d$Y[1:10, , drop = FALSE]
  expect_error(
    multivariate_analysis_pipeline(
      X = d$X, Y = Y_short, maf = d$maf, pip_cutoff_to_skip = 0
    ),
    "X and Y must have the same number of rows"
  )
})

# ---------- maf validation -------------------------------------------------

test_that("multivariate_analysis_pipeline errors when maf has wrong length", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, input validation unreachable")
  d <- make_mv_data()
  expect_error(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = c(0.1, 0.2), pip_cutoff_to_skip = 0
    ),
    "maf must be a numeric vector with length equal to the number of columns in X"
  )
})

test_that("multivariate_analysis_pipeline errors when maf is not numeric", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, input validation unreachable")
  d <- make_mv_data()
  expect_error(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = rep("0.1", ncol(d$X)), pip_cutoff_to_skip = 0
    ),
    "maf must be a numeric vector"
  )
})

test_that("multivariate_analysis_pipeline errors when maf values are below 0", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, input validation unreachable")
  d <- make_mv_data()
  bad_maf <- d$maf
  bad_maf[1] <- -0.1
  expect_error(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = bad_maf, pip_cutoff_to_skip = 0
    ),
    "maf values must be between 0 and 1"
  )
})

test_that("multivariate_analysis_pipeline errors when maf values exceed 1", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, input validation unreachable")
  d <- make_mv_data()
  bad_maf <- d$maf
  bad_maf[2] <- 1.5
  expect_error(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = bad_maf, pip_cutoff_to_skip = 0
    ),
    "maf values must be between 0 and 1"
  )
})

test_that("multivariate_analysis_pipeline accepts maf of all zeros past validation", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, input validation unreachable")
  d <- make_mv_data()
  zero_maf <- rep(0, ncol(d$X))
  # maf = 0 is between 0 and 1 so should pass maf validation;
  # the pipeline will proceed past validation into skip_conditions
  result <- tryCatch(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = zero_maf, pip_cutoff_to_skip = 0
    ),
    error = function(e) e
  )
  # If it errors, the error should NOT be about maf validation
  if (inherits(result, "error")) {
    expect_false(grepl("maf values must be between", result$message))
  }
})

test_that("multivariate_analysis_pipeline accepts maf of all ones past validation", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, input validation unreachable")
  d <- make_mv_data()
  one_maf <- rep(1, ncol(d$X))
  # maf = 1 is between 0 and 1 so passes maf validation
  result <- tryCatch(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = one_maf, pip_cutoff_to_skip = 0
    ),
    error = function(e) e
  )
  # If it errors, the error should NOT be about maf validation
  if (inherits(result, "error")) {
    expect_false(grepl("maf values must be between", result$message))
  }
})

# ===========================================================================
# mvsusieR dependency check
# ===========================================================================

test_that("multivariate_analysis_pipeline errors when mvsusieR is not installed", {
  skip_if(requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR is installed, skipping not-installed test")
  d <- make_mv_data()
  expect_error(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = d$maf, pip_cutoff_to_skip = 0
    ),
    "mvsusieR"
  )
})

test_that("multivariate_analysis_pipeline error message includes install URL", {
  skip_if(requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR is installed, skipping not-installed test")
  d <- make_mv_data()
  expect_error(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = d$maf, pip_cutoff_to_skip = 0
    ),
    "github.com/stephenslab/mvsusieR"
  )
})

# ===========================================================================
# pip_cutoff_to_skip validation (tested via pipeline)
#
# skip_conditions is a local function inside the pipeline, so it can only
# be reached through the pipeline. When mvsusieR is not installed, we can
# still test cases where skip_conditions runs before mvsusieR is needed
# (when pip_cutoff_to_skip = 0, skip_conditions keeps all columns and
# the pipeline proceeds to mvsusieR check).
# ===========================================================================

test_that("pipeline with pip_cutoff_to_skip=0 passes skip_conditions and reaches mvsusieR check", {
  skip_if(requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR is installed, skipping not-installed test")
  d <- make_mv_data()
  # pip_cutoff_to_skip = 0 means keep all columns (no susie call needed)
  # Should reach mvsusieR check
  expect_error(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = d$maf, pip_cutoff_to_skip = 0
    ),
    "mvsusieR"
  )
})

test_that("pipeline with pip_cutoff_to_skip wrong length vector errors", {
  d <- make_mv_data() # Y has 3 columns
  # Provide a vector of length 2 (not 1 and not ncol(Y)=3)
  # This should fail in skip_conditions before mvsusieR is called
  # But skip_conditions is inside the pipeline after the mvsusieR check,
  # so if mvsusieR is not installed, it will error on mvsusieR first.
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, cannot reach skip_conditions")
  expect_error(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = d$maf, pip_cutoff_to_skip = c(0.1, 0.2)
    ),
    "pip_cutoff_to_skip must be a single number or a vector of the same length"
  )
})

test_that("pipeline with pip_cutoff_to_skip vector matching ncol(Y) is accepted", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, cannot reach skip_conditions")
  d <- make_mv_data() # Y has 3 columns
  # Provide a vector of length 3 with all zeros to bypass susie calls
  # This should pass skip_conditions validation
  # It may still fail later in the pipeline, but not on pip_cutoff_to_skip
  result <- tryCatch(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = d$maf,
      pip_cutoff_to_skip = rep(0, 3)
    ),
    error = function(e) e
  )
  # If it errors, the error should not be about pip_cutoff_to_skip
  if (inherits(result, "error")) {
    expect_false(grepl("pip_cutoff_to_skip", result$message))
  }
})

# ===========================================================================
# filter_X_Y_missing behavior (tested indirectly through pipeline)
#
# Since filter_X_Y_missing is a local function, we test its logic by
# providing data that triggers its filtering paths.
# ===========================================================================

test_that("pipeline returns empty list when Y has all-NA rows leaving no data", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, cannot reach filter_X_Y_missing")
  n <- 10
  p <- 5
  r <- 2
  set.seed(55)
  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n, ncol = p)
  rownames(X) <- paste0("s", 1:n)
  colnames(X) <- paste0("chr1:", 1:p * 100, ":A:G")
  # All rows of Y are NA
  Y <- matrix(NA_real_, nrow = n, ncol = r)
  rownames(Y) <- rownames(X)
  colnames(Y) <- paste0("c", 1:r)
  maf <- rep(0.3, p)
  result <- multivariate_analysis_pipeline(
    X = X, Y = Y, maf = maf, pip_cutoff_to_skip = 0
  )
  expect_true(is.list(result))
  expect_equal(length(result), 0)
})

# ===========================================================================
# Edge case: single-column Y (should warn and return empty list)
# ===========================================================================

test_that("pipeline warns and returns empty list when Y has single column", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, cannot test single-column behavior")
  d <- make_mv_data(r = 1)
  # With pip_cutoff_to_skip = 0, skip_conditions keeps the single column,
  # then checks ncol(Y_filtered) <= 1 which triggers warning + NULL return
  expect_warning(
    result <- multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = d$maf, pip_cutoff_to_skip = 0
    ),
    "After filtering.*1 context left"
  )
  expect_true(is.list(result))
  expect_equal(length(result), 0)
})

test_that("multivariate_analysis_pipeline errors when maf is NULL", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, input validation unreachable")
  d <- make_mv_data()
  expect_error(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = NULL, pip_cutoff_to_skip = 0
    ),
    "maf must be a numeric vector"
  )
})

# ===========================================================================
# skip_conditions auto-cutoff (negative pip_cutoff_to_skip)
# ===========================================================================

test_that("skip_conditions with negative pip_cutoff_to_skip auto-computes threshold and skips all", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, cannot reach skip_conditions")
  skip_if(!requireNamespace("susieR", quietly = TRUE),
          "susieR not installed")
  d <- make_mv_data()
  # Mock susieR::susie to return all-low PIPs so the auto-cutoff branch
  # filters out every condition. Pipeline should warn + return empty list.
  # multivariate_pipeline.R imports susieR::susie via NAMESPACE, so mock the
  # binding in the pecotmr namespace.
  local_mocked_bindings(
    susie = function(X, Y, ...) list(pip = rep(1e-6, ncol(X))),
  )
  result <- expect_warning(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = d$maf, pip_cutoff_to_skip = -1
    ),
    "After filtering"
  )
  expect_true(is.list(result))
  expect_equal(length(result), 0)
})

# ===========================================================================
# LD reference filtering
# ===========================================================================

test_that("multivariate_pipeline filters X by LD reference variants and short-circuits on empty w0", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed")
  d <- make_mv_data()
  filter_called <- FALSE
  mrmash_called <- FALSE
  # filter_variants_by_ld_reference: keep first half of variants
  local_mocked_bindings(
    filter_variants_by_ld_reference = function(variant_ids, ld_reference_meta_file, ...) {
      filter_called <<- TRUE
      kept <- seq_len(length(variant_ids) %/% 2)
      list(data = variant_ids[kept], idx = kept)
    },
    # Short-circuit before mvSuSiE/mr.mash internals: return w0 of all "null"
    # so rescale_cov_w0 yields length 0 -> pipeline returns list().
    mrmash_wrapper = function(X, Y, ...) {
      mrmash_called <<- TRUE
      list(V = diag(ncol(Y)), w0 = c(null = 1.0), w1 = matrix(1, nrow = ncol(X), ncol = 1))
    },
  )
  result <- multivariate_analysis_pipeline(
    X = d$X, Y = d$Y, maf = d$maf,
    ld_reference_meta_file = "fake_ld_meta.tsv",
    pip_cutoff_to_skip = 0
  )
  expect_true(filter_called)
  expect_true(mrmash_called)
  expect_true(is.list(result))
  expect_equal(length(result), 0)
})

# ===========================================================================
# empty w0_updated returns empty list
# ===========================================================================

test_that("pipeline returns empty list when rescale_cov_w0 yields length 0", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed")
  d <- make_mv_data()
  local_mocked_bindings(
    mrmash_wrapper = function(X, Y, ...) {
      # Only a "null" component -> rescale_cov_w0 strips it, length = 0.
      list(V = diag(ncol(Y)), w0 = c(null = 1.0), w1 = matrix(1, nrow = ncol(X), ncol = 1))
    },
  )
  result <- multivariate_analysis_pipeline(
    X = d$X, Y = d$Y, maf = d$maf, pip_cutoff_to_skip = 0
  )
  expect_true(is.list(result))
  expect_equal(length(result), 0)
})

# ===========================================================================
# initialize_mvsusie_prior with provided data_driven_prior_matrices
# ===========================================================================

# ===========================================================================
# skip_conditions keeps columns when PIP exceeds cutoff (Tier 1)
# ===========================================================================

test_that("skip_conditions keeps columns when top_model_pip exceeds cutoff", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, cannot reach skip_conditions")
  d <- make_mv_data()
  p <- ncol(d$X)
  # Mock susie so that at least one PIP exceeds the cutoff for every condition
  local_mocked_bindings(
    susie = function(X, Y, ...) list(pip = c(0.9, rep(0.01, ncol(X) - 1)))
  )
  # Use tryCatch to catch the downstream mrmash_wrapper error
  # The message from skip_conditions tells us both columns were kept
  expect_message(
    tryCatch(
      suppressWarnings(multivariate_analysis_pipeline(
        X = d$X, Y = d$Y, maf = d$maf, pip_cutoff_to_skip = 0.5
      )),
      error = function(e) NULL
    ),
    "After filtering by potential association signals, Y has 3 contexts left"
  )
})

# ===========================================================================
# initialize_mvsusie_prior with provided data_driven_prior_matrices
# ===========================================================================

test_that("initialize_mvsusie_prior runs with provided data_driven_prior_matrices", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed")
  d <- make_mv_data()
  r <- ncol(d$Y)
  # Build a non-NULL prior matrices spec keyed under a non-"null" group name
  # so rescale_cov_w0 produces non-empty output.
  prior_U <- list(
    udd_1 = matrix(0.5, r, r) + 0.5 * diag(r),
    udd_2 = diag(r)
  )
  for (k in seq_along(prior_U)) {
    rownames(prior_U[[k]]) <- colnames(prior_U[[k]]) <- colnames(d$Y)
  }
  prior_mats <- list(U = prior_U, w = c(udd_1 = 0.5, udd_2 = 0.5))

  # mrmash returns w0 with names matching "udd_1_<...>" so rescale_cov_w0
  # group prefix "udd_1"/"udd_2" survives.
  fake_w0 <- c("udd_1_a" = 0.4, "udd_2_a" = 0.4, "null" = 0.2)
  local_mocked_bindings(
    mrmash_wrapper = function(X, Y, ...) {
      list(V = diag(ncol(Y)), w0 = fake_w0,
           w1 = matrix(0.1, nrow = ncol(X), ncol = 1))
    },
    susie_post_processor = function(...) list(),
  )
  # Mock mvsusieR::mvsusie + create_mixture_prior to avoid heavy fits
  local_mocked_bindings(
    mvsusie = function(...) list(pip = rep(0.1, ncol(d$X)),
                                 sets = list(cs = NULL)),
    create_mixture_prior = function(...) list(matrices = prior_U, weights = c(0.5, 0.5)),
    .package = "mvsusieR"
  )

  result <- multivariate_analysis_pipeline(
    X = d$X, Y = d$Y, maf = d$maf,
    pip_cutoff_to_skip = 0,
    data_driven_prior_matrices = prior_mats,
    twas_weights = FALSE
  )
  expect_true(is.list(result))
  # initialize_mvsusie_prior path stores reweighted_mixture_prior
  expect_true("reweighted_mixture_prior" %in% names(result))
  expect_true("reweighted_mixture_prior_cv" %in% names(result))
  expect_true("mvsusie_fitted" %in% names(result))
})
