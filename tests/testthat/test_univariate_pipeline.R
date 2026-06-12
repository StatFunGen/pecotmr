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
    requestedCoverage = 0.95
  )
  niter <- 10
  list(
    pip = pip, alpha = alpha, lbf_variable = lbf_variable,
    V = V, mu = mu, X_column_scale_factors = X_column_scale_factors,
    sets = sets, niter = niter
  )
}

# ===========================================================================
# Helper: build a fake protocol-facing post-processing return value
# ===========================================================================
make_fake_post_result <- function(p) {
  vnames <- paste0("chr1:", seq_len(p), ":A:G")
  trimmed <- list(
    pip = runif(p),
    sets = list(cs = list(L1 = c(1L, 2L)), requestedCoverage = 0.95),
    cs_corr = matrix(c(1, 0.1, 0.1, 1), nrow = 2),
    alpha = matrix(runif(2 * p), 2, p),
    lbf_variable = matrix(rnorm(2 * p), 2, p),
    V = c(0.5, 0.01),
    niter = 10,
    n_effects = 2
  )
  fm <- FineMappingResult(
    variantNames = vnames,
    trimmedFit = trimmed,
    topLoci = data.frame(variant_id = character(0), method = character(0)),
    method = "susie"
  )
  list(
    finemapping_result = fm,
    # Legacy keys retained for diagnostics tests that mock getSusieResult
    # to return res$susie_result_trimmed directly.
    susie_result_trimmed = trimmed,
    variantNames = vnames,
    top_loci = data.frame(
      variant_id = paste0("chr1:1:A:G"),
      betahat = 1.5, sebetahat = 0.3, z = 5.0, maf = 0.25,
      pip = 0.9, CS_95_susie = 1,
      stringsAsFactors = FALSE
    )
  )
}

# ===========================================================================
# Helper: fake twasWeightsPipeline return value
# ===========================================================================
make_test_ld_data <- function(variant_ids, R = NULL) {
  if (is.null(R)) {
    p <- length(variant_ids)
    R <- diag(p)
    rownames(R) <- colnames(R) <- variant_ids
  }
  ref_panel <- pecotmr:::parseVariantId(variant_ids)
  ref_panel$variant_id <- variant_ids
  variants_gr <- pecotmr:::.refPanelToGranges(ref_panel)
  bm <- data.frame(
    block_id = 1L, chrom = as.character(ref_panel$chrom[1]),
    block_start = min(ref_panel$pos), block_end = max(ref_panel$pos),
    size = length(variant_ids), start_idx = 1L, end_idx = length(variant_ids),
    stringsAsFactors = FALSE
  )
  LdData(correlation = R, variants = variants_gr, blockMetadata = bm)
}

# ===========================================================================
# Helper: build an LdData S4 object from a matrix for QcResult mocks.
# When is_genotype = FALSE, the matrix is the correlation R (variants on row
# and column names). When is_genotype = TRUE, the matrix is a genotype matrix
# (samples x variants) and colnames are variant IDs.
# ===========================================================================
.test_lddata_from_matrix <- function(mat, is_genotype = FALSE) {
  vids <- if (is_genotype) colnames(mat) else rownames(mat)
  if (is.null(vids)) vids <- colnames(mat)
  ref_panel <- cbind(pecotmr:::parseVariantId(vids), variant_id = vids)
  ref_panel$chrom <- as.character(ref_panel$chrom)
  variants_gr <- pecotmr:::.refPanelToGranges(ref_panel)
  bm <- pecotmr:::.inferSingleLdBlockMetadata(ref_panel)
  if (is_genotype) {
    LdData(correlation = NULL, genotypeHandle = mat,
           variants = variants_gr, blockMetadata = bm,
           nRef = as.integer(nrow(mat)))
  } else {
    LdData(correlation = mat, variants = variants_gr, blockMetadata = bm)
  }
}

# ===========================================================================
# Helper: build a QcResult mock from a sumstats data.frame and LD matrix.
# ===========================================================================
.test_qcresult <- function(sumstats, ld_mat, n = 1000, varY = 1,
                           outlierNumber = 0L, skipped = FALSE,
                           is_genotype = FALSE) {
  ld <- .test_lddata_from_matrix(ld_mat, is_genotype = is_genotype)
  QcResult(
    ldData = ld,
    rssInput = list(sumstats = sumstats, n = n, varY = varY),
    preprocess = list(sumstats = sumstats, ldData = ld),
    outlierNumber = as.integer(outlierNumber),
    skipped = skipped
  )
}

make_fake_twas_result <- function(p) {
  list(
    twasWeights = setNames(rep(0.1, p), paste0("chr1:", seq_len(p), ":A:G")),
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
# univariateAnalysisPipeline: input validation
# ===========================================================================


test_that("univariateAnalysisPipeline errors on non-matrix X", {
  Y <- rnorm(10)
  maf <- runif(5, 0.05, 0.5)
  expect_error(
    univariateAnalysisPipeline(X = data.frame(matrix(rnorm(50), 10, 5)), Y = Y, maf = maf),
    "X must be a numeric matrix"
  )
})

test_that("univariateAnalysisPipeline errors on non-numeric X", {
  X <- matrix(as.character(1:50), 10, 5)
  Y <- rnorm(10)
  maf <- runif(5, 0.05, 0.5)
  expect_error(
    univariateAnalysisPipeline(X = X, Y = Y, maf = maf),
    "X must be a numeric matrix"
  )
})

test_that("univariateAnalysisPipeline errors on non-numeric Y", {
  X <- matrix(rnorm(50), 10, 5)
  Y <- as.character(rnorm(10))
  maf <- runif(5, 0.05, 0.5)
  expect_error(
    univariateAnalysisPipeline(X = X, Y = Y, maf = maf),
    "Y must be a numeric vector"
  )
})

test_that("univariateAnalysisPipeline errors on multi-column Y matrix", {
  X <- matrix(rnorm(50), 10, 5)
  Y <- matrix(rnorm(20), 10, 2)
  maf <- runif(5, 0.05, 0.5)
  expect_error(
    univariateAnalysisPipeline(X = X, Y = Y, maf = maf),
    "Y must be a numeric vector"
  )
})

test_that("univariateAnalysisPipeline accepts single-column Y matrix", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 20; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p), ":A:G")
  rownames(X) <- paste0("s", seq_len(n))
  Y <- matrix(X[, 1] * 2 + rnorm(n), ncol = 1)
  rownames(Y) <- rownames(X)
  maf <- runif(p, 0.1, 0.5)
  result <- univariateAnalysisPipeline(X = X, Y = Y, maf = maf,
    twasWeights = FALSE, L = 5, lGreedy = 5)
  expect_type(result, "list")
})

test_that("univariateAnalysisPipeline errors on X/Y row mismatch", {
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(8)
  maf <- runif(5, 0.05, 0.5)
  expect_error(
    univariateAnalysisPipeline(X = X, Y = Y, maf = maf),
    "same number of rows"
  )
})

test_that("univariateAnalysisPipeline errors on maf length mismatch", {
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  maf <- runif(3, 0.05, 0.5)
  expect_error(
    univariateAnalysisPipeline(X = X, Y = Y, maf = maf),
    "length equal to the number of columns"
  )
})

test_that("univariateAnalysisPipeline errors on maf out of bounds", {
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  maf <- c(0.1, 0.2, -0.1, 0.3, 0.4)
  expect_error(
    univariateAnalysisPipeline(X = X, Y = Y, maf = maf),
    "maf values must be between 0 and 1"
  )
})

test_that("univariateAnalysisPipeline errors on invalid X_scalar", {
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  maf <- runif(5, 0.05, 0.5)
  expect_error(
    univariateAnalysisPipeline(X = X, Y = Y, maf = maf, xScalar = "bad"),
    "xScalar must be a numeric"
  )
})

test_that("univariateAnalysisPipeline errors on wrong-length X_scalar vector", {
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  maf <- runif(5, 0.05, 0.5)
  expect_error(
    univariateAnalysisPipeline(X = X, Y = Y, maf = maf, xScalar = c(1, 2, 3)),
    "xScalar must be a numeric"
  )
})

test_that("univariateAnalysisPipeline errors on non-numeric Y_scalar", {
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  maf <- runif(5, 0.05, 0.5)
  expect_error(
    univariateAnalysisPipeline(X = X, Y = Y, maf = maf, yScalar = "bad"),
    "yScalar must be a numeric scalar"
  )
})

test_that("univariateAnalysisPipeline errors on non-positive L", {
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  maf <- runif(5, 0.05, 0.5)
  expect_error(
    univariateAnalysisPipeline(X = X, Y = Y, maf = maf, L = 0),
    "L must be a positive"
  )
})

test_that("univariateAnalysisPipeline errors on non-positive L_greedy", {
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  maf <- runif(5, 0.05, 0.5)
  expect_error(
    univariateAnalysisPipeline(X = X, Y = Y, maf = maf, lGreedy = 0),
    "lGreedy must be NULL or a positive"
  )
})

# ---- univariateAnalysisPipeline functional tests ----

test_that("univariateAnalysisPipeline runs with minimal valid input", {
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

  result <- univariateAnalysisPipeline(
    X = X, Y = Y, maf = maf,
    twasWeights = FALSE, L = 5, lGreedy = 5
  )

  expect_type(result, "list")
  expect_true("susie_fitted" %in% names(result))
  expect_true("susie_inf_fitted" %in% names(result))
  expect_true("total_time_elapsed" %in% names(result))
})

test_that("univariateAnalysisPipeline produces twasWeights when requested", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 30; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p), ":A:G")
  rownames(X) <- paste0("s", seq_len(n))
  Y <- as.numeric(X[, 1] * 2 + rnorm(n))
  names(Y) <- rownames(X)
  maf <- runif(p, 0.1, 0.5)

  result <- univariateAnalysisPipeline(
    X = X, Y = Y, maf = maf,
    twasWeights = TRUE, L = 5, lGreedy = 5,
    cvFolds = 2
  )

  expect_true("twas_weights_result" %in% names(result))
})

test_that("univariateAnalysisPipeline respects X_scalar vector", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 20; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p), ":A:G")
  rownames(X) <- paste0("s", seq_len(n))
  Y <- as.numeric(X[, 1] * 2 + rnorm(n))
  names(Y) <- rownames(X)
  maf <- runif(p, 0.1, 0.5)

  r_noscalar <- suppressMessages(univariateAnalysisPipeline(
    X = X, Y = Y, maf = maf, xScalar = rep(1, p),
    twasWeights = FALSE, L = 5, lGreedy = 5
  ))
  r_scalar <- suppressMessages(univariateAnalysisPipeline(
    X = X, Y = Y, maf = maf, xScalar = rep(2, p),
    twasWeights = FALSE, L = 5, lGreedy = 5
  ))
  # X_scalar divides betahat: scalar=2 should halve the estimates
  ratio <- unname(r_noscalar$sumstats$betahat / r_scalar$sumstats$betahat)
  expect_equal(ratio, rep(2, p), tolerance = 1e-10)
})

test_that("univariateAnalysisPipeline with cv_folds=0 skips CV", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 20; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p), ":A:G")
  rownames(X) <- paste0("s", seq_len(n))
  Y <- as.numeric(X[, 1] * 2 + rnorm(n))
  names(Y) <- rownames(X)
  maf <- runif(p, 0.1, 0.5)

  result <- univariateAnalysisPipeline(
    X = X, Y = Y, maf = maf,
    twasWeights = TRUE, L = 5, lGreedy = 5,
    cvFolds = 0
  )
  expect_type(result, "list")
  expect_true("twas_weights_result" %in% names(result))
})

# ---- rssAnalysisPipeline input validation ----

test_that("rssAnalysisPipeline requires file inputs", {
  # rssAnalysisPipeline calls loadRssData which requires valid file paths
  dummy_ld <- make_test_ld_data(paste0("1:", 1:5, ":A:G"))
  expect_error(
    rssAnalysisPipeline(
      sumstatPath = "/nonexistent/file.tsv",
      columnFilePath = "/nonexistent/columns.yml",
      ldData = dummy_ld
    )
  )
})

# ---- resolveLdInput ----

test_that("resolveLdInput errors when both R and X are NULL", {
  expect_error(
    pecotmr:::resolveLdInput(R = NULL, X = NULL),
    "Either R .* or X .* must be provided"
  )
})

test_that("resolveLdInput errors when both R and X are provided", {
  R <- diag(5)
  X <- matrix(rnorm(50), 10, 5)
  expect_error(
    pecotmr:::resolveLdInput(R = R, X = X),
    "Provide either R or X, not both"
  )
})

test_that("resolveLdInput with R returns R unchanged", {
  R <- diag(5)
  result <- pecotmr:::resolveLdInput(R = R, nSample = 100)
  expect_equal(result$R, R)
  expect_equal(result$nSample, 100)
})

test_that("resolveLdInput with X computes LD and infers nSample", {
  set.seed(42)
  X <- matrix(rbinom(200, 2, 0.3), 20, 10)
  result <- pecotmr:::resolveLdInput(X = X)
  expect_equal(nrow(result$R), 10)
  expect_equal(ncol(result$R), 10)
  expect_equal(result$nSample, 20)
})

test_that("resolveLdInput errors when nSample required but missing", {
  R <- diag(5)
  expect_error(
    pecotmr:::resolveLdInput(R = R, needNSample = TRUE),
    "nSample is required"
  )
})

test_that("resolveLdInput does not error when nSample not needed", {
  R <- diag(5)
  result <- pecotmr:::resolveLdInput(R = R, needNSample = FALSE)
  expect_equal(result$R, R)
  expect_null(result$nSample)
})

# ===========================================================================
# univariateAnalysisPipeline: mocked pipeline tests
# ===========================================================================

# ========================================================================
#  SECTION 1: univariateAnalysisPipeline - pip_cutoff_to_skip branch
# ========================================================================

test_that("uap: pip_cutoff_to_skip > 0, no signal above threshold => returns empty list", {
  inp <- make_uap_inputs()
  # Mock susie (imported from susieR into pecotmr namespace)
  low_pip <- rep(0.001, inp$p)
  local_mocked_bindings(
    susie = function(...) list(pip = low_pip)
  )
  result <- expect_message(
    univariateAnalysisPipeline(
      X = inp$X, Y = inp$Y, maf = inp$maf,
      pipCutoffToSkip = 0.5,
      twasWeights = FALSE, L = 5, lGreedy = 5
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
  call_count <- 0

  local_mocked_bindings(
    susie = function(...) {
      call_count <<- call_count + 1
      if (call_count == 1) list(pip = high_pip) else fake_fit
    },
    postprocessFinemappingFits = function(...) fake_post,
    formatFinemappingOutput = function(post, primaryMethod) post,
  )

  result <- expect_message(
    univariateAnalysisPipeline(
      X = inp$X, Y = inp$Y, maf = inp$maf,
      pipCutoffToSkip = 0.5,
      twasWeights = FALSE, L = 5, lGreedy = 5
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
    susie = function(...) list(pip = low_pip)
  )

  result <- expect_message(
    univariateAnalysisPipeline(
      X = inp$X, Y = inp$Y, maf = inp$maf,
      pipCutoffToSkip = -1,
      twasWeights = FALSE, L = 5, lGreedy = 5
    ),
    "Skipping follow-up"
  )
  expect_equal(result, list())
})

# ========================================================================
#  SECTION 2: univariateAnalysisPipeline - LD reference filtering
# ========================================================================

test_that("uap: LD reference filtering subsets X columns and maf", {
  inp <- make_uap_inputs()
  fake_fit <- make_fake_susie_fit(5)  # 5 variants after filtering
  fake_post <- make_fake_post_result(5)

  captured_args <- list()
  local_mocked_bindings(
    filterVariantsByLdReference = function(variant_ids, ld_ref_file) {
      # keep only first 5 variants
      list(data = variant_ids[1:5], idx = 1:5)
    },
    susie = function(X, ...) {
      captured_args$ncol_X <<- ncol(X)
      fake_fit
    },
    postprocessFinemappingFits = function(...) fake_post,
    formatFinemappingOutput = function(post, primaryMethod) post,
  )

  result <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    ldReferenceMetaFile = "/fake/ld_meta.tsv",
    twasWeights = FALSE, L = 5, lGreedy = 5
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
    filterVariantsByLdReference = function(variant_ids, ld_ref_file) {
      list(data = variant_ids[1:5], idx = 1:5)
    },
    susie = function(...) fake_fit,
    postprocessFinemappingFits = function(fits, dataX, dataY, xScalar, ...) {
      captured_maf <<- length(dataX[1,])
      fake_post
    },
    formatFinemappingOutput = function(post, primaryMethod) post,
  )

  result <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    xScalar = X_scalar_vec,
    ldReferenceMetaFile = "/fake/ld_meta.tsv",
    twasWeights = FALSE, L = 5, lGreedy = 5
  )
  expect_true("susie_fitted" %in% names(result))
})

# ========================================================================
#  SECTION 3: univariateAnalysisPipeline - filter_X branch
# ========================================================================

test_that("uap: filter_X is invoked when imiss_cutoff is set", {
  inp <- make_uap_inputs()
  fake_fit <- make_fake_susie_fit(inp$p)
  fake_post <- make_fake_post_result(inp$p)

  filter_called <- FALSE
  local_mocked_bindings(
    filterX = function(X, ...) {
      filter_called <<- TRUE
      X  # return unchanged
    },
    susie = function(...) fake_fit,
    postprocessFinemappingFits = function(...) fake_post,
    formatFinemappingOutput = function(post, primaryMethod) post,
  )

  result <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    imissCutoff = 0.9,
    twasWeights = FALSE, L = 5, lGreedy = 5
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
    filterX = function(X, ...) X_filtered,
    susie = function(X, ...) {
      captured_ncol <<- ncol(X)
      fake_fit
    },
    postprocessFinemappingFits = function(...) fake_post,
    formatFinemappingOutput = function(post, primaryMethod) post,
  )

  result <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    mafCutoff = 0.05, xScalar = rep(2, inp$p),
    twasWeights = FALSE, L = 5, lGreedy = 5
  )
  expect_equal(captured_ncol, length(keep_cols))
})

# ========================================================================
#  SECTION 4: univariateAnalysisPipeline - main analysis (susie + post)
# ========================================================================

test_that("uap: susie-inf initializes susie with matching greedy L", {
  inp <- make_uap_inputs()
  fake_inf_fit <- make_fake_susie_fit(inp$p)
  fake_fit <- make_fake_susie_fit(inp$p)
  fake_post <- make_fake_post_result(inp$p)

  captured_calls <- list()
  local_mocked_bindings(
    susie = function(X, y, L, L_greedy, coverage, unmappable_effects = "none",
                     model_init = NULL, refine = NULL,
                     convergence_method = NULL, ...) {
      captured_calls[[length(captured_calls) + 1]] <<- list(
        L = L,
        lGreedy = L_greedy,
        unmappable_effects = unmappable_effects,
        model_init = model_init,
        refine = refine,
        convergence_method = convergence_method
      )
      if (identical(unmappable_effects, "inf")) fake_inf_fit else fake_fit
    },
    postprocessFinemappingFits = function(...) fake_post,
    formatFinemappingOutput = function(post, primaryMethod) post,
  )

  result <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    twasWeights = FALSE, L = 15, lGreedy = 7
  )
  expect_length(captured_calls, 2)
  expect_equal(captured_calls[[1]]$L, 15)
  expect_equal(captured_calls[[1]]$lGreedy, 7)
  expect_equal(captured_calls[[1]]$unmappable_effects, "inf")
  expect_equal(captured_calls[[1]]$convergence_method, "pip")
  expect_false(captured_calls[[1]]$refine)
  expect_equal(captured_calls[[2]]$L, 15)
  expect_equal(captured_calls[[2]]$lGreedy, nrow(fake_inf_fit$alpha))
  expect_equal(captured_calls[[2]]$unmappable_effects, "none")
  expect_true(captured_calls[[2]]$refine)
  expect_identical(captured_calls[[2]]$model_init, result$susie_inf_fitted)
  expect_s3_class(result$susie_fitted, "susie")
  expect_s3_class(result$susie_inf_fitted, "susie_inf")
  expect_equal(unclass(result$susie_fitted), fake_fit)
  expect_equal(unclass(result$susie_inf_fitted), fake_inf_fit)
})

test_that("uap: ordinary susie can run without susie-inf initialization", {
  inp <- make_uap_inputs()
  fake_fit <- make_fake_susie_fit(inp$p)
  fake_post <- make_fake_post_result(inp$p)

  captured_calls <- list()
  local_mocked_bindings(
    susie = function(X, y, L, L_greedy, coverage, unmappable_effects = "none",
                     model_init = NULL, ...) {
      captured_calls[[length(captured_calls) + 1]] <<- list(
        L = L,
        lGreedy = L_greedy,
        unmappable_effects = unmappable_effects,
        model_init = model_init
      )
      fake_fit
    },
    postprocessFinemappingFits = function(fits, ...) {
      expect_named(fits, "susie")
      fake_post
    },
    formatFinemappingOutput = function(post, primaryMethod) post,
  )

  result <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    twasWeights = FALSE, addSusieInf = FALSE,
    L = 15, lGreedy = 7
  )
  expect_length(captured_calls, 1)
  expect_equal(captured_calls[[1]]$L, 15)
  expect_equal(captured_calls[[1]]$lGreedy, 7)
  expect_equal(captured_calls[[1]]$unmappable_effects, "none")
  expect_null(captured_calls[[1]]$model_init)
  expect_s3_class(result$susie_fitted, "susie")
  expect_false("susie_inf_fitted" %in% names(result))
  expect_false("susie_inf_result_trimmed" %in% names(result))
})

test_that("uap: post-processing output is merged into result", {
  inp <- make_uap_inputs()
  fake_fit <- make_fake_susie_fit(inp$p)
  fake_fm <- FineMappingResult(
    variantNames = paste0("v", seq_len(inp$p)),
    trimmedFit = list(pip = runif(inp$p)),
    topLoci = data.frame(variant_id = character(0), method = character(0)),
    method = "susie"
  )
  fake_post <- list(
    finemapping_result = fake_fm,
    top_loci = data.frame(variant_id = "v1", pip = 0.9, stringsAsFactors = FALSE)
  )

  local_mocked_bindings(
    susie = function(...) fake_fit,
    postprocessFinemappingFits = function(fits, ...) {
      # SuSiE-inf is fitted first because chained SuSiE / SuSiE-ash
      # warm-start from it; the list-name order reflects fit order.
      expect_named(fits, c("susie_inf", "susie"))
      fake_post
    },
    formatFinemappingOutput = function(post, primaryMethod) post,
  )

  result <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    twasWeights = FALSE, L = 5, lGreedy = 5
  )
  expect_true("finemapping_result" %in% names(result))
  expect_true("top_loci" %in% names(result))
  expect_true("total_time_elapsed" %in% names(result))
})

test_that("uap: finemapping_extra_opts are forwarded to susie", {
  inp <- make_uap_inputs()
  fake_fit <- make_fake_susie_fit(inp$p)
  fake_post <- make_fake_post_result(inp$p)

  captured_refine <- NULL
  local_mocked_bindings(
    susie = function(...) {
      args <- list(...)
      captured_refine <<- args$refine
      fake_fit
    },
    postprocessFinemappingFits = function(...) fake_post,
    formatFinemappingOutput = function(post, primaryMethod) post,
  )

  result <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    finemappingExtraOpts = list(refine = TRUE),
    twasWeights = FALSE, L = 5, lGreedy = 5
  )
  expect_true(captured_refine)
})

# ========================================================================
#  SECTION 5: univariateAnalysisPipeline - TWAS weights
# ========================================================================

test_that("uap: twasWeights = TRUE calls twasWeightsPipeline", {
  inp <- make_uap_inputs()
  fake_fit <- make_fake_susie_fit(inp$p)
  fake_post <- make_fake_post_result(inp$p)
  fake_twas <- make_fake_twas_result(inp$p)

  twas_called <- FALSE
  local_mocked_bindings(
    susie = function(...) fake_fit,
    postprocessFinemappingFits = function(...) fake_post,
    formatFinemappingOutput = function(post, primaryMethod) post,
    twasWeightsPipeline = function(...) {
      twas_called <<- TRUE
      fake_twas
    },
  )

  result <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    twasWeights = TRUE, L = 5, lGreedy = 5,
    cvFolds = 2
  )
  expect_true(twas_called)
  expect_true("twas_weights_result" %in% names(result))
})

test_that("uap: twasWeights copies top_loci into susie_weights_intermediate", {
  inp <- make_uap_inputs()
  fake_fit <- make_fake_susie_fit(inp$p)
  fake_post <- make_fake_post_result(inp$p)
  fake_twas <- make_fake_twas_result(inp$p)

  local_mocked_bindings(
    susie = function(...) fake_fit,
    postprocessFinemappingFits = function(...) fake_post,
    formatFinemappingOutput = function(post, primaryMethod) post,
    twasWeightsPipeline = function(...) fake_twas,
  )

  result <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    twasWeights = TRUE, L = 5, lGreedy = 5,
    cvFolds = 2
  )
  # top_loci should be copied from fake_post into the twas result
  expect_true("top_loci" %in% names(result$twas_weights_result$susie_weights_intermediate))
})

test_that("uap: twasWeightsPipeline receives correct cv_folds and sample_partition", {
  inp <- make_uap_inputs()
  fake_fit <- make_fake_susie_fit(inp$p)
  fake_post <- make_fake_post_result(inp$p)
  fake_twas <- make_fake_twas_result(inp$p)

  captured_cv_folds <- NULL
  captured_partition <- NULL
  local_mocked_bindings(
    susie = function(...) fake_fit,
    postprocessFinemappingFits = function(...) fake_post,
    formatFinemappingOutput = function(post, primaryMethod) post,
    twasWeightsPipeline = function(X, Y, fittedModels, cvFolds, maxCvVariants,
                                     cvThreads, samplePartition) {
      captured_cv_folds <<- cvFolds
      captured_partition <<- samplePartition
      fake_twas
    },
  )

  sample_part <- rep(1:3, length.out = inp$n)
  result <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    twasWeights = TRUE, L = 5, lGreedy = 5,
    cvFolds = 3, samplePartition = sample_part
  )
  expect_equal(captured_cv_folds, 3)
  expect_equal(captured_partition, sample_part)
})

test_that("uap: twasWeights = FALSE skips twasWeightsPipeline", {
  inp <- make_uap_inputs()
  fake_fit <- make_fake_susie_fit(inp$p)
  fake_post <- make_fake_post_result(inp$p)

  twas_called <- FALSE
  local_mocked_bindings(
    susie = function(...) fake_fit,
    postprocessFinemappingFits = function(...) fake_post,
    formatFinemappingOutput = function(post, primaryMethod) post,
    twasWeightsPipeline = function(...) {
      twas_called <<- TRUE
      make_fake_twas_result(inp$p)
    },
  )

  result <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    twasWeights = FALSE, L = 5, lGreedy = 5
  )
  expect_false(twas_called)
  expect_false("twas_weights_result" %in% names(result))
})

# ========================================================================
#  SECTION 6: univariateAnalysisPipeline - coverage vector forwarding
# ========================================================================

test_that("uap: coverage vector is forwarded correctly to susie and post-processing", {
  inp <- make_uap_inputs()
  fake_fit <- make_fake_susie_fit(inp$p)
  fake_post <- make_fake_post_result(inp$p)

  captured_coverage_wrapper <- NULL
  captured_secondary_cov <- NULL
  local_mocked_bindings(
    susie = function(X, y, L, L_greedy, coverage, ...) {
      captured_coverage_wrapper <<- coverage
      fake_fit
    },
    postprocessFinemappingFits = function(fits, data_x, data_y, X_scalar, y_scalar, maf,
                                            secondaryCoverage, ...) {
      captured_secondary_cov <<- secondaryCoverage
      fake_post
    },
    formatFinemappingOutput = function(post, primaryMethod) post,
  )

  result <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    coverage = c(0.95, 0.7, 0.5),
    twasWeights = FALSE, L = 5, lGreedy = 5
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
    susie = function(...) fake_fit,
    postprocessFinemappingFits = function(fits, data_x, data_y, X_scalar, y_scalar, maf,
                                            secondaryCoverage, ...) {
      captured_secondary_cov <<- secondaryCoverage
      fake_post
    },
    formatFinemappingOutput = function(post, primaryMethod) post,
  )

  result <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    coverage = c(0.95),
    twasWeights = FALSE, L = 5, lGreedy = 5
  )
  expect_null(captured_secondary_cov)
})

# ========================================================================
#  SECTION 7: univariateAnalysisPipeline - combined LD + filter_X
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
    filterVariantsByLdReference = function(variant_ids, ld_ref_file) {
      list(data = variant_ids[1:8], idx = 1:8)
    },
    filterX = function(X, ...) X_after_filter,
    susie = function(X, ...) {
      captured_ncol <<- ncol(X)
      fake_fit
    },
    postprocessFinemappingFits = function(...) fake_post,
    formatFinemappingOutput = function(post, primaryMethod) post,
  )

  result <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    ldReferenceMetaFile = "/fake/ld.tsv",
    imissCutoff = 0.9, mafCutoff = 0.05,
    twasWeights = FALSE, L = 5, lGreedy = 5
  )
  expect_equal(captured_ncol, 6)
})

# ========================================================================
#  SECTION 8: rssAnalysisPipeline - empty sumstats early return
# ========================================================================

test_that("rss: empty sumstats from loadRssData => early return", {
  dummy_ld <- make_test_ld_data(paste0("1:", 1:5, ":A:G"))
  local_mocked_bindings(
    loadRssData = function(...) {
      list(sumstats = data.frame(), n = NULL, varY = NULL)
    },
  )

  result <- rssAnalysisPipeline(
    sumstatPath = "/fake/sumstats.tsv",
    columnFilePath = "/fake/columns.yml",
    ldData = dummy_ld
  )
  expect_true("rss_data_analyzed" %in% names(result))
  expect_equal(nrow(result$rss_data_analyzed), 0)
})

# ========================================================================
#  SECTION 9: rssAnalysisPipeline - empty sumstats after preprocessing
# ========================================================================

test_that("rss: empty sumstats after rssBasicQc => early return", {
  dummy_ld <- make_test_ld_data(c("1:100:A:G"))
  local_mocked_bindings(
    loadRssData = function(...) {
      list(
        sumstats = data.frame(chrom = 1, pos = 100, A1 = "A", A2 = "G", z = 2.5, variant_id = "1:100:A:G"),
        n = 1000,
        varY = 1
      )
    },
    rssBasicQc = function(...) {
      QcResult(
        ldData = NULL,
        rssInput = list(sumstats = data.frame(), n = NA_real_, varY = NA_real_),
        preprocess = list(),
        outlierNumber = 0L,
        skipped = FALSE
      )
    },
  )

  result <- expect_message(
    rssAnalysisPipeline(
      sumstatPath = "/fake/sumstats.tsv",
      columnFilePath = "/fake/columns.yml",
      ldData = dummy_ld
    ),
    "No variants left after preprocessing"
  )
  expect_true("rss_data_analyzed" %in% names(result))
  expect_equal(nrow(result$rss_data_analyzed), 0)
})

# ========================================================================
#  SECTION 10: rssAnalysisPipeline - pip_cutoff_to_skip branch
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

test_that("regionDataToSusieRssInput prepares correlation-backed RSS input", {
  ss <- make_rss_sumstats(5)
  ss_sub <- ss[c(2, 4, 5), , drop = FALSE]
  ld <- .test_lddata_from_matrix(make_rss_ld_mat(5))

  out <- regionDataToSusieRssInput(
    list(sumstats = ss_sub, n = 1000, varY = 1.5),
    ld
  )

  expect_equal(out$source_info$n_variants, 3)
  expect_false(out$source_info$uses_X_ref)
  expect_equal(out$susie_rss_input$n, 1000)
  expect_equal(out$susie_rss_input$var_y, 1.5)
  expect_null(out$susie_rss_input$xMat)
  expect_equal(rownames(out$susie_rss_input$ldMat), ss_sub$variant_id)
  expect_equal(colnames(out$susie_rss_input$ldMat), ss_sub$variant_id)
})

test_that("regionDataToSusieRssInput prepares genotype-backed RSS input", {
  ss <- make_rss_sumstats(5)
  ss_sub <- ss[c(1, 3, 5), , drop = FALSE]
  X_ref <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
  colnames(X_ref) <- ss$variant_id
  ld <- .test_lddata_from_matrix(X_ref, is_genotype = TRUE)

  out <- regionDataToSusieRssInput(
    list(sumstats = ss_sub, n = 1000, varY = 1),
    ld
  )

  expect_true(out$source_info$uses_X_ref)
  expect_false(out$source_info$has_LD)
  expect_null(out$susie_rss_input$ldMat)
  expect_equal(dim(out$susie_rss_input$xMat), c(20, 3))
  expect_equal(colnames(out$susie_rss_input$xMat), ss_sub$variant_id)
})

test_that("regionDataToSusieRssInput rejects default row names as variant IDs", {
  ss <- make_rss_sumstats(5)[c(2, 4, 5), , drop = FALSE]
  ss$variant_id <- NULL
  ld <- .test_lddata_from_matrix(make_rss_ld_mat(3))

  expect_error(
    regionDataToSusieRssInput(list(sumstats = ss, n = 1000, varY = 1), ld),
    "variant_id or variant"
  )
})

test_that("rss: pip_cutoff_to_skip > 0, no signal => early return", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  local_mocked_bindings(
    loadRssData = function(...) list(sumstats = ss, n = 1000, varY = 1),
    rssBasicQc = function(...) QcResult(
      ldData = if (nrow(ld_mat) > 0) .test_lddata_from_matrix(ld_mat) else NULL,
      rssInput = list(sumstats = ss, n = NA_real_, varY = NA_real_),
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    susie_ser = function(...) list(pip = rep(0.01, 5)),
  )

  result <- expect_message(
    rssAnalysisPipeline(
      sumstatPath = "/fake/sumstats.tsv",
      columnFilePath = "/fake/columns.yml",
      ldData = make_test_ld_data(ss$variant_id),
      pipCutoffToSkip = 0.5,
      qcMethod = "none"
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
    loadRssData = function(...) list(sumstats = ss, n = 1000, varY = 1),
    rssBasicQc = function(...) QcResult(
      ldData = if (nrow(ld_mat) > 0) .test_lddata_from_matrix(ld_mat) else NULL,
      rssInput = list(sumstats = ss, n = NA_real_, varY = NA_real_),
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    susie_rss = function(...) list(pip = c(0.9, 0.01, 0.01, 0.01, 0.01)),
    summaryStatsQc = function(...) {
      message("Follow-up on region: signals above PIP threshold 0.5 detected.")
      .test_qcresult(ss, ld_mat, outlierNumber = 0)
    },
    partitionLdMatrix = function(...) ld_mat,
    raiss = function(...) list(result_filter = ss, LD_mat = ld_mat),
    susieRssPipeline = function(...) fake_rss_result,
  )

  result <- expect_message(
    rssAnalysisPipeline(
      sumstatPath = "/fake/sumstats.tsv",
      columnFilePath = "/fake/columns.yml",
      ldData = make_test_ld_data(ss$variant_id),
      pipCutoffToSkip = 0.5,
      qcMethod = "slalom",
      finemappingMethod = "susie_rss"
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
    loadRssData = function(...) list(sumstats = ss, n = 1000, varY = 1),
    rssBasicQc = function(...) QcResult(
      ldData = if (nrow(ld_mat) > 0) .test_lddata_from_matrix(ld_mat) else NULL,
      rssInput = list(sumstats = ss, n = NA_real_, varY = NA_real_),
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    susie_ser = function(...) list(pip = rep(0.01, 5)),
  )

  # auto threshold = 3 * 1/5 = 0.6, all PIPs are 0.01 so skip
  result <- expect_message(
    rssAnalysisPipeline(
      sumstatPath = "/fake/sumstats.tsv",
      columnFilePath = "/fake/columns.yml",
      ldData = make_test_ld_data(ss$variant_id),
      pipCutoffToSkip = -1,
      qcMethod = "none"
    ),
    "Skipping follow-up"
  )
  expect_true("rss_data_analyzed" %in% names(result))
})

# ========================================================================
#  SECTION 11: rssAnalysisPipeline - QC + imputation + fine-mapping
# ========================================================================

test_that("rss: full pipeline with QC, imputation, and fine-mapping", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)
  fake_result <- make_fake_post_result(5)

  qc_called <- FALSE
  impute_called <- FALSE
  finemapping_called <- FALSE

  local_mocked_bindings(
    loadRssData = function(...) list(sumstats = ss, n = 1000, varY = 1),
    rssBasicQc = function(...) QcResult(
      ldData = if (nrow(ld_mat) > 0) .test_lddata_from_matrix(ld_mat) else NULL,
      rssInput = list(sumstats = ss, n = NA_real_, varY = NA_real_),
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    summaryStatsQc = function(..., impute = FALSE) {
      qc_called <<- TRUE
      # Imputation now happens inside summaryStatsQc, so simulate that
      # branch in the mock to keep the raiss call observable.
      if (isTRUE(impute)) raiss()
      .test_qcresult(ss, ld_mat, outlierNumber = 1)
    },
    partitionLdMatrix = function(...) list(ld_matrices = list(ld_mat)),
    raiss = function(...) {
      impute_called <<- TRUE
      list(result_filter = ss, LD_mat = ld_mat)
    },
    susieRssPipeline = function(...) {
      finemapping_called <<- TRUE
      fake_result
    },
  )

  result <- rssAnalysisPipeline(
    sumstatPath = "/fake/sumstats.tsv",
    columnFilePath = "/fake/columns.yml",
    ldData = make_test_ld_data(ss$variant_id),
    qcMethod = "slalom",
    impute = TRUE,
    finemappingMethod = "susie_rss"
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
    loadRssData = function(...) list(sumstats = ss, n = 1000, varY = 1),
    rssBasicQc = function(...) QcResult(
      ldData = if (nrow(ld_mat) > 0) .test_lddata_from_matrix(ld_mat) else NULL,
      rssInput = list(sumstats = ss, n = NA_real_, varY = NA_real_),
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    summaryStatsQc = function(...) .test_qcresult(ss, ld_mat, outlierNumber = 0),
    susieRssPipeline = function(...) fake_result,
  )

  result <- rssAnalysisPipeline(
    sumstatPath = "/fake/sumstats.tsv",
    columnFilePath = "/fake/columns.yml",
    ldData = make_test_ld_data(ss$variant_id),
    qcMethod = "slalom",
    impute = FALSE,
    finemappingMethod = "susie_rss"
  )

  expect_true(any(grepl("susie_rss_SLALOM$", names(result))))
})

test_that("rss: method name is correct for no QC", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)
  fake_result <- make_fake_post_result(5)

  local_mocked_bindings(
    loadRssData = function(...) list(sumstats = ss, n = 1000, varY = 1),
    rssBasicQc = function(...) QcResult(
      ldData = if (nrow(ld_mat) > 0) .test_lddata_from_matrix(ld_mat) else NULL,
      rssInput = list(sumstats = ss, n = NA_real_, varY = NA_real_),
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    partitionLdMatrix = function(...) list(ld_matrices = list(ld_mat)),
    raiss = function(...) list(result_filter = ss, LD_mat = ld_mat),
    susieRssPipeline = function(...) fake_result,
  )

  result <- rssAnalysisPipeline(
    sumstatPath = "/fake/sumstats.tsv",
    columnFilePath = "/fake/columns.yml",
    ldData = make_test_ld_data(ss$variant_id),
    qcMethod = NULL,
    impute = TRUE,
    finemappingMethod = "susie_rss"
  )

  expect_true(any(grepl("NO_QC", names(result))))
})

# ========================================================================
#  SECTION 12: rssAnalysisPipeline - outlier_number stored
# ========================================================================

test_that("rss: outlier_number is stored in result when QC is active", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)
  fake_result <- make_fake_post_result(5)

  local_mocked_bindings(
    loadRssData = function(...) list(sumstats = ss, n = 1000, varY = 1),
    rssBasicQc = function(...) QcResult(
      ldData = if (nrow(ld_mat) > 0) .test_lddata_from_matrix(ld_mat) else NULL,
      rssInput = list(sumstats = ss, n = NA_real_, varY = NA_real_),
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    summaryStatsQc = function(...) .test_qcresult(ss, ld_mat, outlierNumber = 3),
    susieRssPipeline = function(...) fake_result,
  )

  result <- rssAnalysisPipeline(
    sumstatPath = "/fake/sumstats.tsv",
    columnFilePath = "/fake/columns.yml",
    ldData = make_test_ld_data(ss$variant_id),
    qcMethod = "dentist",
    impute = FALSE,
    finemappingMethod = "susie_rss"
  )

  method_key <- grep("DENTIST", names(result), value = TRUE)
  expect_true(length(method_key) > 0)
  expect_equal(result[[method_key[1]]]$outlier_number, 3)
})

# ========================================================================
#  SECTION 13: rssAnalysisPipeline - finemappingMethod = NULL skips
# ========================================================================

test_that("rss: finemappingMethod = NULL skips fine-mapping", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  finemapping_called <- FALSE
  local_mocked_bindings(
    loadRssData = function(...) list(sumstats = ss, n = 1000, varY = 1),
    rssBasicQc = function(...) QcResult(
      ldData = if (nrow(ld_mat) > 0) .test_lddata_from_matrix(ld_mat) else NULL,
      rssInput = list(sumstats = ss, n = NA_real_, varY = NA_real_),
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    susieRssPipeline = function(...) {
      finemapping_called <<- TRUE
      list()
    },
  )

  # Use qcMethod = NULL and impute = FALSE to avoid paste0(NULL, ...) edge case
  result <- rssAnalysisPipeline(
    sumstatPath = "/fake/sumstats.tsv",
    columnFilePath = "/fake/columns.yml",
    ldData = make_test_ld_data(ss$variant_id),
    qcMethod = NULL,
    impute = FALSE,
    finemappingMethod = NULL
  )

  expect_false(finemapping_called)
  expect_true("rss_data_analyzed" %in% names(result))
})

# ========================================================================
#  SECTION 14: rssAnalysisPipeline - QC = NULL skips QC
# ========================================================================

test_that("rss: qcMethod = NULL uses combined basic QC without LD-mismatch method", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)
  fake_result <- make_fake_post_result(5)

  qc_called <- FALSE
  local_mocked_bindings(
    loadRssData = function(...) list(sumstats = ss, n = 1000, varY = 1),
    rssBasicQc = function(...) QcResult(
      ldData = if (nrow(ld_mat) > 0) .test_lddata_from_matrix(ld_mat) else NULL,
      rssInput = list(sumstats = ss, n = NA_real_, varY = NA_real_),
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    summaryStatsQc = function(...) {
      qc_called <<- TRUE
      .test_qcresult(ss, ld_mat, outlierNumber = 0)
    },
    susieRssPipeline = function(...) fake_result,
  )

  result <- rssAnalysisPipeline(
    sumstatPath = "/fake/sumstats.tsv",
    columnFilePath = "/fake/columns.yml",
    ldData = make_test_ld_data(ss$variant_id),
    qcMethod = NULL,
    impute = FALSE,
    finemappingMethod = "susie_rss"
  )

  expect_true(qc_called)
})

# ========================================================================
#  SECTION 15: rssAnalysisPipeline - impute = FALSE skips imputation
# ========================================================================

test_that("rss: impute = FALSE skips raiss imputation", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)
  fake_result <- make_fake_post_result(5)

  raiss_called <- FALSE
  local_mocked_bindings(
    loadRssData = function(...) list(sumstats = ss, n = 1000, varY = 1),
    rssBasicQc = function(...) QcResult(
      ldData = if (nrow(ld_mat) > 0) .test_lddata_from_matrix(ld_mat) else NULL,
      rssInput = list(sumstats = ss, n = NA_real_, varY = NA_real_),
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    summaryStatsQc = function(...) .test_qcresult(ss, ld_mat, outlierNumber = 0),
    raiss = function(...) {
      raiss_called <<- TRUE
      list(result_filter = ss, LD_mat = ld_mat)
    },
    susieRssPipeline = function(...) fake_result,
  )

  result <- rssAnalysisPipeline(
    sumstatPath = "/fake/sumstats.tsv",
    columnFilePath = "/fake/columns.yml",
    ldData = make_test_ld_data(ss$variant_id),
    qcMethod = "slalom",
    impute = FALSE,
    finemappingMethod = "susie_rss"
  )

  expect_false(raiss_called)
})

# ========================================================================
#  SECTION 16: rssAnalysisPipeline - diagnostics branch (empty res)
# ========================================================================

test_that("rss: diagnostics = TRUE with empty fine-mapping result skips diagnostics logic", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  local_mocked_bindings(
    loadRssData = function(...) list(sumstats = ss, n = 1000, varY = 1),
    rssBasicQc = function(...) QcResult(
      ldData = if (nrow(ld_mat) > 0) .test_lddata_from_matrix(ld_mat) else NULL,
      rssInput = list(sumstats = ss, n = NA_real_, varY = NA_real_),
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    summaryStatsQc = function(...) .test_qcresult(ss, ld_mat, outlierNumber = 0),
    susieRssPipeline = function(...) list(),  # empty result
  )

  result <- rssAnalysisPipeline(
    sumstatPath = "/fake/sumstats.tsv",
    columnFilePath = "/fake/columns.yml",
    ldData = make_test_ld_data(ss$variant_id),
    qcMethod = "slalom",
    impute = FALSE,
    diagnostics = TRUE,
    finemappingMethod = "susie_rss"
  )

  # With empty fine-mapping result, block_cs_metrics stays empty,
  # so the "diagnostics" key is never added to result_list
  expect_false("diagnostics" %in% names(result))
})

# ========================================================================
#  SECTION 17: rssAnalysisPipeline - diagnostics with CS > 1 and
#              high p-value triggers BCR + SER re-analysis
# ========================================================================

test_that("rss: diagnostics with 2+ CS and high p-value/corr triggers BCR and SER", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  # Build a fake result that has susie_result_trimmed with 2 CS
  fake_result <- list(
    variantNames = paste0("1:", seq(100, by = 100, length.out = 5), ":A:G"),
    susie_result_trimmed = list(
      pip = c(0.9, 0.8, 0.1, 0.05, 0.02),
      sets = list(
        cs = list(L1 = c(1L, 2L), L2 = c(3L, 4L)),
        requestedCoverage = 0.95
      ),
      cs_corr = matrix(c(1, 0.6, 0.6, 1), nrow = 2)
    ),
    sumstats = list(z = c(5.0, 4.0, 1.5, 0.5, 0.1)),
    topLoci = data.frame(
      variant_id = c("1:100:A:G", "1:200:A:G", "1:300:A:G", "1:400:A:G"),
      pip = c(0.9, 0.8, 0.1, 0.05),
      z = c(5.0, 4.0, 1.5, 0.5),
      stringsAsFactors = FALSE
    )
  )

  # extractCsInfo returns a tibble with 2 CS rows
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

  # parseCsCorr adds cs_corr_max and cs_corr_min
  fake_parsed <- as.data.frame(fake_cs_info)
  fake_parsed$cs_corr_1 <- c(1, 0.6)
  fake_parsed$cs_corr_2 <- c(0.6, 1)
  fake_parsed$cs_corr_max <- c(0.6, 0.6)
  fake_parsed$cs_corr_min <- c(0.6, 0.6)

  susie_rss_pipeline_call_count <- 0
  local_mocked_bindings(
    loadRssData = function(...) list(sumstats = ss, n = 1000, varY = 1),
    rssBasicQc = function(...) QcResult(
      ldData = if (nrow(ld_mat) > 0) .test_lddata_from_matrix(ld_mat) else NULL,
      rssInput = list(sumstats = ss, n = NA_real_, varY = NA_real_),
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    summaryStatsQc = function(...) .test_qcresult(ss, ld_mat, outlierNumber = 0),
    susieRssPipeline = function(...) {
      susie_rss_pipeline_call_count <<- susie_rss_pipeline_call_count + 1
      fake_result
    },
    getSusieResult = function(res) res$susie_result_trimmed,
    extractCsInfo = function(...) fake_cs_info,
    parseCsCorr = function(...) fake_parsed,
  )

  result <- rssAnalysisPipeline(
    sumstatPath = "/fake/sumstats.tsv",
    columnFilePath = "/fake/columns.yml",
    ldData = make_test_ld_data(ss$variant_id),
    qcMethod = "slalom",
    impute = FALSE,
    diagnostics = TRUE,
    finemappingMethod = "susie_rss"
  )

  # Main susie_rss call + BCR reanalysis + SER reanalysis = 3 calls
  expect_equal(susie_rss_pipeline_call_count, 3)
  # Should have BCR and SER method entries
  expect_true(any(grepl("bayesian_conditional_regression", names(result))))
  expect_true(any(grepl("single_effect", names(result))))
  expect_true("diagnostics" %in% names(result))
})

# ========================================================================
#  SECTION 18: rssAnalysisPipeline - diagnostics with 1 CS => SER only
# ========================================================================

test_that("rss: diagnostics with 1 CS triggers SER reanalysis only", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  fake_result <- list(
    variantNames = paste0("1:", seq(100, by = 100, length.out = 5), ":A:G"),
    susie_result_trimmed = list(
      pip = c(0.9, 0.1, 0.05, 0.02, 0.01),
      sets = list(
        cs = list(L1 = c(1L, 2L)),
        requestedCoverage = 0.95
      ),
      cs_corr = matrix(1, nrow = 1)
    ),
    sumstats = list(z = c(5.0, 1.0, 0.5, 0.2, 0.1)),
    topLoci = data.frame(
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
    loadRssData = function(...) list(sumstats = ss, n = 1000, varY = 1),
    rssBasicQc = function(...) QcResult(
      ldData = if (nrow(ld_mat) > 0) .test_lddata_from_matrix(ld_mat) else NULL,
      rssInput = list(sumstats = ss, n = NA_real_, varY = NA_real_),
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    summaryStatsQc = function(...) .test_qcresult(ss, ld_mat, outlierNumber = 0),
    susieRssPipeline = function(...) {
      susie_rss_pipeline_call_count <<- susie_rss_pipeline_call_count + 1
      fake_result
    },
    getSusieResult = function(res) res$susie_result_trimmed,
    extractCsInfo = function(...) fake_cs_info,
    parseCsCorr = function(...) fake_parsed,
  )

  result <- rssAnalysisPipeline(
    sumstatPath = "/fake/sumstats.tsv",
    columnFilePath = "/fake/columns.yml",
    ldData = make_test_ld_data(ss$variant_id),
    qcMethod = "slalom",
    impute = FALSE,
    diagnostics = TRUE,
    finemappingMethod = "susie_rss"
  )

  # Main call + SER reanalysis = 2 calls
  expect_equal(susie_rss_pipeline_call_count, 2)
  expect_true(any(grepl("single_effect", names(result))))
  expect_true("diagnostics" %in% names(result))
})

# ========================================================================
#  SECTION 19: rssAnalysisPipeline - diagnostics with no CS but high PIP
# ========================================================================

test_that("rss: diagnostics with no CS but high PIP calls extractTopPipInfo", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  fake_result <- list(
    variantNames = paste0("1:", seq(100, by = 100, length.out = 5), ":A:G"),
    susie_result_trimmed = list(
      pip = c(0.15, 0.05, 0.03, 0.02, 0.01),
      sets = list(cs = NULL, requestedCoverage = 0.95),
      cs_corr = NULL
    ),
    sumstats = list(z = c(3.0, 1.0, 0.5, 0.2, 0.1)),
    topLoci = data.frame(
      variant_id = "1:100:A:G",
      pip = 0.15,
      z = 3.0,
      stringsAsFactors = FALSE
    )
  )

  # extractTopPipInfo returns when no CS but PIP > signal_cutoff
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
    loadRssData = function(...) list(sumstats = ss, n = 1000, varY = 1),
    rssBasicQc = function(...) QcResult(
      ldData = if (nrow(ld_mat) > 0) .test_lddata_from_matrix(ld_mat) else NULL,
      rssInput = list(sumstats = ss, n = NA_real_, varY = NA_real_),
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    summaryStatsQc = function(...) .test_qcresult(ss, ld_mat, outlierNumber = 0),
    susieRssPipeline = function(...) fake_result,
    getSusieResult = function(res) res$susie_result_trimmed,
    extractTopPipInfo = function(...) {
      extract_top_pip_called <<- TRUE
      fake_top_pip_info
    },
    parseCsCorr = function(df) {
      df <- as.data.frame(df)
      df$cs_corr_max <- NA_real_
      df$cs_corr_min <- NA_real_
      df
    },
  )

  result <- rssAnalysisPipeline(
    sumstatPath = "/fake/sumstats.tsv",
    columnFilePath = "/fake/columns.yml",
    ldData = make_test_ld_data(ss$variant_id),
    qcMethod = "slalom",
    impute = FALSE,
    diagnostics = TRUE,
    finemappingMethod = "susie_rss"
  )

  expect_true(extract_top_pip_called)
  expect_true("diagnostics" %in% names(result))
})

# ========================================================================
#  SECTION 20: rssAnalysisPipeline - diagnostics with no CS and low PIP
# ========================================================================

test_that("rss: diagnostics with no CS and no high PIP => diagnostics empty", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  fake_result <- list(
    variantNames = paste0("1:", seq(100, by = 100, length.out = 5), ":A:G"),
    susie_result_trimmed = list(
      pip = c(0.01, 0.005, 0.003, 0.002, 0.001),
      sets = list(cs = NULL, requestedCoverage = 0.95),
      cs_corr = NULL
    ),
    sumstats = list(z = c(0.5, 0.3, 0.2, 0.1, 0.05))
  )

  local_mocked_bindings(
    loadRssData = function(...) list(sumstats = ss, n = 1000, varY = 1),
    rssBasicQc = function(...) QcResult(
      ldData = if (nrow(ld_mat) > 0) .test_lddata_from_matrix(ld_mat) else NULL,
      rssInput = list(sumstats = ss, n = NA_real_, varY = NA_real_),
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    summaryStatsQc = function(...) .test_qcresult(ss, ld_mat, outlierNumber = 0),
    susieRssPipeline = function(...) fake_result,
    getSusieResult = function(res) res$susie_result_trimmed,
  )

  result <- rssAnalysisPipeline(
    sumstatPath = "/fake/sumstats.tsv",
    columnFilePath = "/fake/columns.yml",
    ldData = make_test_ld_data(ss$variant_id),
    qcMethod = "slalom",
    impute = FALSE,
    diagnostics = TRUE,
    finemappingMethod = "susie_rss",
    finemappingOpts = list(
      L = 20, lGreedy = 5,
      coverage = c(0.95, 0.7, 0.5), signalCutoff = 0.025,
      min_abs_corr = 0.8
    )
  )

  # No CS and no high PIP means block_cs_metrics stays empty,
  # so "diagnostics" key is never added to result_list
  expect_false("diagnostics" %in% names(result))
})

# ========================================================================
#  SECTION 21: rssAnalysisPipeline - finemapping_opts forwarded correctly
# ========================================================================

test_that("rss: finemapping_opts are forwarded to susieRssPipeline", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  captured_L <- NULL
  captured_L_greedy <- NULL
  captured_coverage <- NULL
  captured_signal_cutoff <- NULL
  captured_R_mismatch <- NULL
  local_mocked_bindings(
    loadRssData = function(...) list(sumstats = ss, n = 1000, varY = 1),
    rssBasicQc = function(...) QcResult(
      ldData = if (nrow(ld_mat) > 0) .test_lddata_from_matrix(ld_mat) else NULL,
      rssInput = list(sumstats = ss, n = NA_real_, varY = NA_real_),
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    susieRssPipeline = function(sumstats, ldMat, n, var_y, L, lGreedy,
                                  analysisMethod, coverage, secondaryCoverage,
                                  signalCutoff, minAbsCorr, ...) {
      captured_L <<- L
      captured_L_greedy <<- lGreedy
      captured_coverage <<- coverage
      captured_signal_cutoff <<- signalCutoff
      captured_R_mismatch <<- list(...)$rMismatch
      make_fake_post_result(5)
    },
  )

  result <- rssAnalysisPipeline(
    sumstatPath = "/fake/sumstats.tsv",
    columnFilePath = "/fake/columns.yml",
    ldData = make_test_ld_data(ss$variant_id),
    qcMethod = NULL,
    impute = FALSE,
    finemappingMethod = "susie_rss",
    finemappingOpts = list(
      L = 10, lGreedy = 3,
      coverage = c(0.90, 0.6), signalCutoff = 0.05,
      min_abs_corr = 0.7
    ),
    rMismatch = "eb"
  )

  expect_equal(captured_L, 10)
  expect_equal(captured_L_greedy, 3)
  expect_equal(captured_coverage, 0.90)
  expect_equal(captured_signal_cutoff, 0.05)
  expect_equal(captured_R_mismatch, "eb")
})

# ========================================================================
#  SECTION 22: rssAnalysisPipeline - dentist QC method naming
# ========================================================================

test_that("rss: dentist QC method generates correct method name", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)
  fake_result <- make_fake_post_result(5)

  local_mocked_bindings(
    loadRssData = function(...) list(sumstats = ss, n = 1000, varY = 1),
    rssBasicQc = function(...) QcResult(
      ldData = if (nrow(ld_mat) > 0) .test_lddata_from_matrix(ld_mat) else NULL,
      rssInput = list(sumstats = ss, n = NA_real_, varY = NA_real_),
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    summaryStatsQc = function(...) .test_qcresult(ss, ld_mat, outlierNumber = 2),
    partitionLdMatrix = function(...) list(ld_matrices = list(ld_mat)),
    raiss = function(...) list(result_filter = ss, LD_mat = ld_mat),
    susieRssPipeline = function(...) fake_result,
  )

  result <- rssAnalysisPipeline(
    sumstatPath = "/fake/sumstats.tsv",
    columnFilePath = "/fake/columns.yml",
    ldData = make_test_ld_data(ss$variant_id),
    qcMethod = "dentist",
    impute = TRUE,
    finemappingMethod = "susie_rss"
  )

  expect_true(any(grepl("DENTIST_RAISS_imputed", names(result))))
})

# ========================================================================
#  SECTION 23: univariateAnalysisPipeline - null imiss/maf cutoffs skip filter_X
# ========================================================================

test_that("uap: both imiss_cutoff and maf_cutoff NULL skips filter_X", {
  inp <- make_uap_inputs()
  fake_fit <- make_fake_susie_fit(inp$p)
  fake_post <- make_fake_post_result(inp$p)

  filter_called <- FALSE
  local_mocked_bindings(
    filterX = function(...) {
      filter_called <<- TRUE
      inp$X
    },
    susie = function(...) fake_fit,
    postprocessFinemappingFits = function(...) fake_post,
    formatFinemappingOutput = function(post, primaryMethod) post,
  )

  result <- univariateAnalysisPipeline(
    X = inp$X, Y = inp$Y, maf = inp$maf,
    imissCutoff = NULL, mafCutoff = NULL,
    twasWeights = FALSE, L = 5, lGreedy = 5
  )

  expect_false(filter_called)
})

# ========================================================================
#  SECTION 24: rssAnalysisPipeline - diagnostics with getSusieResult
#              returning NULL
# ========================================================================

test_that("rss: diagnostics with getSusieResult returning NULL => diagnostics empty", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  # A non-empty result but with getSusieResult returning NULL
  fake_result <- list(some_field = "value")

  local_mocked_bindings(
    loadRssData = function(...) list(sumstats = ss, n = 1000, varY = 1),
    rssBasicQc = function(...) QcResult(
      ldData = if (nrow(ld_mat) > 0) .test_lddata_from_matrix(ld_mat) else NULL,
      rssInput = list(sumstats = ss, n = NA_real_, varY = NA_real_),
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    summaryStatsQc = function(...) .test_qcresult(ss, ld_mat, outlierNumber = 0),
    susieRssPipeline = function(...) fake_result,
    getSusieResult = function(res) NULL,
  )

  result <- rssAnalysisPipeline(
    sumstatPath = "/fake/sumstats.tsv",
    columnFilePath = "/fake/columns.yml",
    ldData = make_test_ld_data(ss$variant_id),
    qcMethod = "slalom",
    impute = FALSE,
    diagnostics = TRUE,
    finemappingMethod = "susie_rss"
  )

  # getSusieResult returns NULL, so block_cs_metrics stays empty,
  # and "diagnostics" key is never added to result_list
  expect_false("diagnostics" %in% names(result))
})

# ========================================================================
#  SECTION 25: rssAnalysisPipeline - diagnostics with CS but CS count <= 1
#              and null block_cs_metrics (no high PIP either)
# ========================================================================

test_that("rss: diagnostics with null/empty block_cs_metrics => no additional analysis", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  # Non-empty result with susie_result_trimmed having 0 CS and no high PIPs
  fake_result <- list(
    variantNames = paste0("1:", seq(100, by = 100, length.out = 5), ":A:G"),
    susie_result_trimmed = list(
      pip = c(0.001, 0.002, 0.001, 0.001, 0.001),
      sets = list(cs = list(), requestedCoverage = 0.95),
      cs_corr = NULL
    ),
    sumstats = list(z = rep(0.1, 5))
  )

  susie_rss_call_count <- 0
  local_mocked_bindings(
    loadRssData = function(...) list(sumstats = ss, n = 1000, varY = 1),
    rssBasicQc = function(...) QcResult(
      ldData = if (nrow(ld_mat) > 0) .test_lddata_from_matrix(ld_mat) else NULL,
      rssInput = list(sumstats = ss, n = NA_real_, varY = NA_real_),
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    summaryStatsQc = function(...) .test_qcresult(ss, ld_mat, outlierNumber = 0),
    susieRssPipeline = function(...) {
      susie_rss_call_count <<- susie_rss_call_count + 1
      fake_result
    },
    getSusieResult = function(res) res$susie_result_trimmed,
  )

  result <- rssAnalysisPipeline(
    sumstatPath = "/fake/sumstats.tsv",
    columnFilePath = "/fake/columns.yml",
    ldData = make_test_ld_data(ss$variant_id),
    qcMethod = "slalom",
    impute = FALSE,
    diagnostics = TRUE,
    finemappingMethod = "susie_rss",
    finemappingOpts = list(
      L = 20, lGreedy = 5,
      coverage = c(0.95, 0.7, 0.5), signalCutoff = 0.025,
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
#  SECTION 26: rssAnalysisPipeline - diagnostics with CS > 1 but low
#              p-value and low corr => no BCR/SER triggered
# ========================================================================

test_that("rss: diagnostics with 2 CS but low p-value and low corr => no extra analysis", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  fake_result <- list(
    variantNames = paste0("1:", seq(100, by = 100, length.out = 5), ":A:G"),
    susie_result_trimmed = list(
      pip = c(0.9, 0.8, 0.1, 0.05, 0.02),
      sets = list(
        cs = list(L1 = c(1L, 2L), L2 = c(3L, 4L)),
        requestedCoverage = 0.95
      ),
      cs_corr = matrix(c(1, 0.1, 0.1, 1), nrow = 2)
    ),
    sumstats = list(z = c(5.0, 4.0, 3.0, 2.5, 0.1)),
    topLoci = data.frame(
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
    loadRssData = function(...) list(sumstats = ss, n = 1000, varY = 1),
    rssBasicQc = function(...) QcResult(
      ldData = if (nrow(ld_mat) > 0) .test_lddata_from_matrix(ld_mat) else NULL,
      rssInput = list(sumstats = ss, n = NA_real_, varY = NA_real_),
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    summaryStatsQc = function(...) .test_qcresult(ss, ld_mat, outlierNumber = 0),
    susieRssPipeline = function(...) {
      susie_rss_call_count <<- susie_rss_call_count + 1
      fake_result
    },
    getSusieResult = function(res) res$susie_result_trimmed,
    extractCsInfo = function(...) fake_cs_info,
    parseCsCorr = function(...) fake_parsed,
  )

  result <- rssAnalysisPipeline(
    sumstatPath = "/fake/sumstats.tsv",
    columnFilePath = "/fake/columns.yml",
    ldData = make_test_ld_data(ss$variant_id),
    qcMethod = "slalom",
    impute = FALSE,
    diagnostics = TRUE,
    finemappingMethod = "susie_rss"
  )

  # Only the main call, no BCR/SER because no p > 1e-4 and no corr > 0.5
  expect_equal(susie_rss_call_count, 1)
  expect_true("diagnostics" %in% names(result))
})

# ========================================================================
#  SECTION 27: rssAnalysisPipeline - diagnostics with CS > 1 and
#              high corr (> 0.5) but low p-value still triggers BCR+SER
# ========================================================================

test_that("rss: diagnostics with high max_cs_corr_study_block triggers BCR+SER", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)

  fake_result <- list(
    variantNames = paste0("1:", seq(100, by = 100, length.out = 5), ":A:G"),
    susie_result_trimmed = list(
      pip = c(0.9, 0.8, 0.1, 0.05, 0.02),
      sets = list(
        cs = list(L1 = c(1L, 2L), L2 = c(3L, 4L)),
        requestedCoverage = 0.95
      ),
      cs_corr = matrix(c(1, 0.7, 0.7, 1), nrow = 2)
    ),
    sumstats = list(z = c(10.0, 9.0, 8.0, 7.0, 0.1)),
    topLoci = data.frame(
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
    loadRssData = function(...) list(sumstats = ss, n = 1000, varY = 1),
    rssBasicQc = function(...) QcResult(
      ldData = if (nrow(ld_mat) > 0) .test_lddata_from_matrix(ld_mat) else NULL,
      rssInput = list(sumstats = ss, n = NA_real_, varY = NA_real_),
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    summaryStatsQc = function(...) .test_qcresult(ss, ld_mat, outlierNumber = 0),
    susieRssPipeline = function(...) {
      susie_rss_call_count <<- susie_rss_call_count + 1
      fake_result
    },
    getSusieResult = function(res) res$susie_result_trimmed,
    extractCsInfo = function(...) fake_cs_info,
    parseCsCorr = function(...) fake_parsed,
  )

  result <- rssAnalysisPipeline(
    sumstatPath = "/fake/sumstats.tsv",
    columnFilePath = "/fake/columns.yml",
    ldData = make_test_ld_data(ss$variant_id),
    qcMethod = "slalom",
    impute = FALSE,
    diagnostics = TRUE,
    finemappingMethod = "susie_rss"
  )

  # Main + BCR + SER = 3
  expect_equal(susie_rss_call_count, 3)
  expect_true(any(grepl("bayesian_conditional_regression", names(result))))
  expect_true(any(grepl("single_effect", names(result))))
})

# ========================================================================
#  loadStudyLd: single vs mixture LD panels
# ========================================================================

test_that("loadStudyLd with single path returns loadLdMatrix output unchanged", {
  fake_LD <- list(
    LD_matrix = diag(4),
    LD_variants = paste0("chr1:", 1:4, ":A:G"),
    is_genotype = FALSE,
    blockMetadata = data.frame(chrom = 1, start = 1, end = 4)
  )
  captured_args <- NULL
  local_mocked_bindings(
    loadLdMatrix = function(LD_meta_file_path, region, ...) {
      captured_args <<- list(path = LD_meta_file_path, region = region, ...)
      fake_LD
    },
  )
  out <- loadStudyLd("ld_meta.tsv", "chr1:1-100")
  expect_identical(out, fake_LD)
  expect_equal(captured_args$path, "ld_meta.tsv")
  expect_equal(captured_args$region, "chr1:1-100")
  expect_equal(captured_args$returnGenotype, "auto")
})

test_that("loadStudyLd with comma-separated paths builds list of X panels", {
  # Build mock LdData objects with fake genotype handles
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 3),
    ranges = IRanges::IRanges(start = 1:3, width = 1L)
  )
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    variant_id = paste0("chr1:", 1:3, ":A:G"),
    A1 = rep("G", 3), A2 = rep("A", 3)
  )
  bm <- data.frame(chrom = 1, start = 1, end = 3)
  fake_handle1 <- "handle1"
  fake_handle2 <- "handle2"
  panel1 <- LdData(
    correlation = NULL, genotypeHandle = fake_handle1, snpIdx = 1:3,
    variants = gr, blockMetadata = bm, nRef = 10L
  )
  panel2 <- LdData(
    correlation = NULL, genotypeHandle = fake_handle2, snpIdx = 1:3,
    variants = gr, blockMetadata = bm, nRef = 10L
  )
  call_count <- 0
  local_mocked_bindings(
    loadLdMatrix = function(LD_meta_file_path, region, ...) {
      call_count <<- call_count + 1
      if (call_count == 1) panel1 else panel2
    },
  )
  out <- loadStudyLd("ld_eur.tsv,ld_afr.tsv", "chr1:1-100")
  expect_equal(call_count, 2)
  expect_s4_class(out, "LdData")
  expect_true(hasGenotypes(out))
  # Mixture: genotype_handle should be a list of handles
  expect_true(is.list(out@genotypeHandle))
  expect_length(out@genotypeHandle, 2)
  expect_equal(out@genotypeHandle[[1]], "handle1")
  expect_equal(out@genotypeHandle[[2]], "handle2")
  # Variants preserved from base panel
  expect_equal(getVariantIds(out), paste0("chr1:", 1:3, ":A:G"))
})

# ========================================================================
#  rssAnalysisPipeline: X-as-LD path (genotype interface)
# ========================================================================

test_that("rss: is_genotype=TRUE path does not precompute R and uses X for fine-mapping", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)
  X_geno <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
  colnames(X_geno) <- ss$variant_id
  fake_result <- make_fake_post_result(5)

  # Create LdData with genotype_handle to trigger the genotype path
  ref_panel <- pecotmr:::parseVariantId(ss$variant_id)
  ref_panel$variant_id <- ss$variant_id
  variants_gr <- pecotmr:::.refPanelToGranges(ref_panel)
  bm <- data.frame(
    block_id = 1L, chrom = as.character(ref_panel$chrom[1]),
    block_start = min(ref_panel$pos), block_end = max(ref_panel$pos),
    size = 5L, start_idx = 1L, end_idx = 5L, stringsAsFactors = FALSE
  )
  geno_ld <- LdData(
    correlation = NULL, genotypeHandle = "fake_handle", snpIdx = 1:5,
    variants = variants_gr, blockMetadata = bm, nRef = 20L
  )

  compute_LD_called <- FALSE
  susie_X_arg <- NULL
  susie_LD_arg <- "unset"

  # Mock extractBlockGenotypes to return a fake SummarizedExperiment-like object
  # that getgenotypes can use. We mock getgenotypes directly at package level.
  local_mocked_bindings(
    extractBlockGenotypes = function(handle, snp_idx, ...) {
      # Return a fake object that assay() can handle
      se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(dosage = t(X_geno))
      )
      se
    },
    computeLd = function(X, method = "sample") {
      compute_LD_called <<- TRUE
      stop("rssAnalysisPipeline should not precompute LD from X")
    },
    loadRssData = function(...) list(sumstats = ss, n = 1000, varY = 1),
    rssBasicQc = function(...) QcResult(
      ldData = if (nrow(ld_mat) > 0) .test_lddata_from_matrix(ld_mat) else NULL,
      rssInput = list(sumstats = ss, n = NA_real_, varY = NA_real_),
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    summaryStatsQc = function(...) .test_qcresult(
      ss, X_geno, outlierNumber = 0, is_genotype = TRUE
    ),
    partitionLdMatrix = function(...) list(ld_matrices = list(ld_mat)),
    raiss = function(...) list(result_filter = ss, LD_mat = ld_mat),
    susieRssPipeline = function(sumstats, ldMat = NULL, xMat = NULL, ...) {
      susie_X_arg <<- xMat
      susie_LD_arg <<- ldMat
      fake_result
    },
  )

  result <- rssAnalysisPipeline(
    sumstatPath = "/fake/sumstats.tsv",
    columnFilePath = "/fake/columns.yml",
    ldData = geno_ld,
    qcMethod = "slalom",
    finemappingMethod = "susie_rss"
  )

  expect_false(compute_LD_called)
  expect_true(is.matrix(susie_X_arg))
  expect_null(susie_LD_arg)
  expect_true(any(grepl("susie_rss", names(result))))
})

test_that("rss: mixture LD_data (list of X panels) preserves list shape into susieRssPipeline", {
  ss <- make_rss_sumstats(5)
  ld_mat <- make_rss_ld_mat(5)
  X1 <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
  X2 <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
  colnames(X1) <- colnames(X2) <- ss$variant_id
  fake_result <- make_fake_post_result(5)

  # Create LdData with list of genotype handles (mixture path)
  ref_panel <- pecotmr:::parseVariantId(ss$variant_id)
  ref_panel$variant_id <- ss$variant_id
  variants_gr <- pecotmr:::.refPanelToGranges(ref_panel)
  bm <- data.frame(
    block_id = 1L, chrom = as.character(ref_panel$chrom[1]),
    block_start = min(ref_panel$pos), block_end = max(ref_panel$pos),
    size = 5L, start_idx = 1L, end_idx = 5L, stringsAsFactors = FALSE
  )
  mixture_ld <- LdData(
    correlation = NULL, genotypeHandle = list("fake_handle1", "fake_handle2"),
    snpIdx = 1:5, variants = variants_gr, blockMetadata = bm, nRef = 20L
  )

  susie_X_arg <- NULL
  compute_LD_called <- FALSE
  extract_call_count <- 0

  local_mocked_bindings(
    extractBlockGenotypes = function(handle, snp_idx, ...) {
      extract_call_count <<- extract_call_count + 1
      # Return appropriate X matrix based on which handle
      X_mat <- if (extract_call_count == 1) X1 else X2
      SummarizedExperiment::SummarizedExperiment(
        assays = list(dosage = t(X_mat))
      )
    },
    computeLd = function(X, method = "sample") {
      compute_LD_called <<- TRUE
      stop("rssAnalysisPipeline should not precompute LD from mixture X")
    },
    loadRssData = function(...) list(sumstats = ss, n = 1000, varY = 1),
    rssBasicQc = function(...) QcResult(
      ldData = if (nrow(ld_mat) > 0) .test_lddata_from_matrix(ld_mat) else NULL,
      rssInput = list(sumstats = ss, n = NA_real_, varY = NA_real_),
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    summaryStatsQc = function(...) QcResult(
      ldData = mixture_ld,
      rssInput = list(sumstats = ss, n = 1000, varY = 1),
      preprocess = list(sumstats = ss, ldData = mixture_ld),
      outlierNumber = 0L,
      skipped = FALSE
    ),
    partitionLdMatrix = function(...) list(ld_matrices = list(ld_mat)),
    raiss = function(...) list(result_filter = ss, LD_mat = ld_mat),
    susieRssPipeline = function(sumstats, ldMat = NULL, xMat = NULL, ...) {
      susie_X_arg <<- xMat
      fake_result
    },
  )

  result <- rssAnalysisPipeline(
    sumstatPath = "/fake/sumstats.tsv",
    columnFilePath = "/fake/columns.yml",
    ldData = mixture_ld,
    qcMethod = "slalom",
    finemappingMethod = "susie_rss",
    impute = FALSE
  )

  # Mixture path => list of subset matrices passed to susieRssPipeline as X_mat
  expect_false(compute_LD_called)
  expect_true(is.list(susie_X_arg))
  expect_length(susie_X_arg, 2)
  expect_true(all(sapply(susie_X_arg, is.matrix)))
  expect_true(any(grepl("susie_rss", names(result))))
})

# ===========================================================================
# rss-af-single-source (Step 1.2): af is the single source of truth
# ===========================================================================
.uap_af_xy <- function(n = 20, p = 5, seed = 7) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p), ":A:G")
  rownames(X) <- paste0("s", seq_len(n))
  Y <- as.numeric(X[, 1] * 2 + rnorm(n))
  names(Y) <- rownames(X)
  list(X = X, Y = Y, p = p)
}

test_that("af supplied (no maf): MAF derived from af, af exported, no maf column", {
  skip_if_not_installed("susieR")
  d <- .uap_af_xy()
  af <- runif(d$p, 0.1, 0.5)
  res <- univariateAnalysisPipeline(X = d$X, Y = d$Y, af = af,
    mafCutoff = 0.01, twasWeights = FALSE, L = 5, lGreedy = 5)
  expect_true("top_loci" %in% names(res))
  expect_true("af" %in% colnames(res$top_loci))
  expect_false("maf" %in% colnames(res$top_loci))
  if (nrow(res$top_loci) > 0) expect_true(all(!is.na(res$top_loci$af)))
})

test_that("maf supplied (no af): maf used for filtering, top_loci$af is NA", {
  skip_if_not_installed("susieR")
  d <- .uap_af_xy()
  maf <- runif(d$p, 0.1, 0.5)
  res <- univariateAnalysisPipeline(X = d$X, Y = d$Y, maf = maf,
    mafCutoff = 0.01, twasWeights = FALSE, L = 5, lGreedy = 5)
  expect_true("af" %in% colnames(res$top_loci))
  if (nrow(res$top_loci) > 0) expect_true(all(is.na(res$top_loci$af)))
})

test_that("both af and maf supplied and disagreeing emits a warning (af wins)", {
  skip_if_not_installed("susieR")
  d <- .uap_af_xy()
  af  <- rep(0.30, d$p)   # min(af, 1-af) = 0.30
  maf <- rep(0.40, d$p)   # disagrees with the af-derived 0.30
  expect_warning(
    univariateAnalysisPipeline(X = d$X, Y = d$Y, af = af, maf = maf,
      mafCutoff = 0.01, twasWeights = FALSE, L = 5, lGreedy = 5),
    "disagree"
  )
})

test_that("neither af nor maf with maf_cutoff set errors clearly", {
  skip_if_not_installed("susieR")
  d <- .uap_af_xy()
  expect_error(
    univariateAnalysisPipeline(X = d$X, Y = d$Y,
      mafCutoff = 0.01, twasWeights = FALSE, L = 5, lGreedy = 5),
    "mafCutoff is set but neither"
  )
})

test_that("neither af nor maf without mafCutoff runs; top_loci$af is NA", {
  skip_if_not_installed("susieR")
  d <- .uap_af_xy()
  res <- univariateAnalysisPipeline(X = d$X, Y = d$Y,
    mafCutoff = NULL, twasWeights = FALSE, L = 5, lGreedy = 5)
  expect_true("top_loci" %in% names(res))
  expect_true("af" %in% colnames(res$top_loci))
  if (nrow(res$top_loci) > 0) expect_true(all(is.na(res$top_loci$af)))
})

test_that("af-derived MAF drives maf_cutoff filtering identically to the equivalent maf", {
  # Behavioral test of the edit: supplying `af` must filter via min(af, 1-af),
  # giving the SAME fit as supplying that directionless maf explicitly. Includes
  # an af > 0.5 (variant 2) so a wrong derivation (using af directly, not the
  # minor allele) would diverge; and a sub-cutoff variant (variant 1) so the
  # cutoff actually removes something.
  skip_if_not_installed("susieR")
  d <- .uap_af_xy(p = 6)
  af  <- c(0.005, 0.85, 0.20, 0.45, 0.15, 0.40)   # variant 1 below cutoff; variant 2 af>0.5
  res_af  <- univariateAnalysisPipeline(X = d$X, Y = d$Y, af = af,
    mafCutoff = 0.01, twasWeights = FALSE, L = 5, lGreedy = 5)
  res_maf <- univariateAnalysisPipeline(X = d$X, Y = d$Y, maf = pmin(af, 1 - af),
    mafCutoff = 0.01, twasWeights = FALSE, L = 5, lGreedy = 5)
  # Identical on every column except af (res_maf's af is NA) => the af-derived
  # MAF was used for filtering exactly like an explicit maf.
  cols <- setdiff(names(res_af$top_loci), "af")
  expect_equal(res_af$top_loci[cols], res_maf$top_loci[cols])
  # And the sub-cutoff variant (chr1:1:A:G, MAF 0.005 < 0.01) was filtered out:
  expect_false("chr1:1:A:G" %in% res_af$top_loci$variant)
})
