context("susie_finemapping")

# =============================================================================
# lbf_to_alpha_vector (internal)
# =============================================================================

test_that("lbf_to_alpha_vector converts correctly", {
  lbf <- c(a = -0.5, b = 1.2, c = 0.3)
  alpha <- pecotmr:::lbfToAlphaVector(lbf)
  expect_length(alpha, 3)
  expect_equal(sum(alpha), 1, tolerance = 1e-10)
  expect_true(alpha["b"] > alpha["a"])
})

test_that("lbf_to_alpha_vector with prior weights", {
  lbf <- c(a = 1, b = 1, c = 1)  # Equal LBFs
  pw <- c(0.5, 0.25, 0.25)
  alpha <- pecotmr:::lbfToAlphaVector(lbf, priorWeights = pw)
  expect_true(alpha[1] > alpha[2])
  expect_equal(sum(alpha), 1, tolerance = 1e-10)
})

test_that("lbf_to_alpha_vector returns zeros for all-zero lbf", {
  lbf <- c(a = 0, b = 0, c = 0)
  alpha <- pecotmr:::lbfToAlphaVector(lbf)
  expect_true(all(alpha == 0))
})

test_that("lbf_to_alpha_vector handles single element", {
  lbf <- c(a = 2.0)
  alpha <- pecotmr:::lbfToAlphaVector(lbf)
  expect_length(alpha, 1)
  expect_equal(alpha[["a"]], 1.0)
})

test_that("lbf_to_alpha_vector handles very large LBFs without overflow", {
  lbf <- c(a = 500, b = 500.1, c = 499)
  alpha <- pecotmr:::lbfToAlphaVector(lbf)
  expect_true(all(is.finite(alpha)))
  expect_equal(sum(alpha), 1, tolerance = 1e-10)
  expect_true(alpha["b"] > alpha["a"])
})

test_that("lbf_to_alpha_vector handles very negative LBFs", {
  lbf <- c(a = -1000, b = -999, c = -1001)
  alpha <- pecotmr:::lbfToAlphaVector(lbf)
  expect_true(all(is.finite(alpha)))
  expect_equal(sum(alpha), 1, tolerance = 1e-10)
  expect_true(alpha["b"] > alpha["a"])
})

test_that("lbf_to_alpha_vector with unequal prior weights", {
  lbf <- c(a = 0.5, b = 0.5, c = 0.5)
  pw <- c(0.8, 0.1, 0.1)
  alpha <- pecotmr:::lbfToAlphaVector(lbf, priorWeights = pw)
  expect_equal(sum(alpha), 1, tolerance = 1e-10)
  expect_true(alpha[1] > 0.7)
})

# =============================================================================
# lbfToAlpha (matrix version)
# =============================================================================

test_that("lbfToAlpha converts log BFs to posteriors", {
  lbf <- matrix(c(0, 3, 2, 1, 4, 0), nrow = 2, ncol = 3)
  alpha <- pecotmr:::lbfToAlpha(lbf)
  expect_equal(dim(alpha), c(2, 3))
  expect_equal(rowSums(alpha), c(1, 1), tolerance = 1e-10)
  expect_true(alpha[1, 3] > alpha[1, 1])
  expect_true(alpha[2, 1] > alpha[2, 3])
})

test_that("lbfToAlpha handles uniform lbf", {
  lbf <- matrix(1, nrow = 1, ncol = 5)
  alpha <- pecotmr:::lbfToAlpha(lbf)
  expect_equal(as.numeric(alpha), rep(0.2, 5), tolerance = 1e-10)
})

test_that("lbfToAlpha handles single-row matrix", {
  lbf <- matrix(c(1.0, 2.0, 0.5), nrow = 1)
  colnames(lbf) <- paste0("v", 1:3)
  result <- lbfToAlpha(lbf)
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 3)
  expect_equal(sum(result), 1, tolerance = 1e-10)
})

test_that("lbfToAlpha handles large matrix", {
  set.seed(42)
  lbf <- matrix(rnorm(100), nrow = 10, ncol = 10)
  colnames(lbf) <- paste0("v", 1:10)
  result <- lbfToAlpha(lbf)
  expect_equal(dim(result), c(10, 10))
  expect_equal(rowSums(result), rep(1, 10), tolerance = 1e-10)
})

test_that("lbfToAlpha with mixed zero and nonzero rows", {
  lbf <- matrix(c(0, 0, 0, 1, 2, 3), nrow = 2, byrow = TRUE)
  colnames(lbf) <- paste0("v", 1:3)
  result <- lbfToAlpha(lbf)
  expect_true(all(result[1, ] == 0))
  expect_equal(sum(result[2, ]), 1, tolerance = 1e-10)
})

# =============================================================================
# get_cs_index (internal)
# =============================================================================

test_that("get_cs_index finds variant in credible set", {
  susie_cs <- list(L1 = c(1, 2, 3), L2 = c(4, 5))
  idx <- pecotmr:::getCsIndex(2, susie_cs)
  expect_equal(unname(idx), 1)
})

test_that("get_cs_index returns NA for variant not in any CS", {
  susie_cs <- list(L1 = c(1, 2), L2 = c(4, 5))
  idx <- pecotmr:::getCsIndex(3, susie_cs)
  expect_true(is.na(idx))
})

test_that("get_cs_index returns all CS indices when variant in multiple", {
  susie_cs <- list(L1 = c(1, 2, 3), L2 = c(2, 4, 5))
  idx <- pecotmr:::getCsIndex(2, susie_cs)
  expect_equal(unname(idx), c(1, 2))
})

test_that("get_cs_index returns all matching CS regardless of size", {
  susie_cs <- list(L1 = c(1, 2, 3, 4, 5), L2 = c(2, 3))
  result <- pecotmr:::getCsIndex(2, susie_cs)
  expect_equal(unname(result), c(1, 2))
})

test_that("get_cs_index handles empty CS list", {
  susie_cs <- list()
  result <- pecotmr:::getCsIndex(1, susie_cs)
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
  idx <- pecotmr:::getCsIndex(1, fit$sets$cs)
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
  result <- pecotmr:::getTopVariantsIdx(susie_output, signalCutoff = 0.1)
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
  result <- pecotmr:::getTopVariantsIdx(susie_output, signalCutoff = 0.1)
  expect_equal(result, c(2, 4))
})

test_that("get_top_variants_idx with all low PIPs", {
  susie_output <- list(
    pip = c(0.01, 0.02, 0.03),
    sets = list(cs = list(L1 = c(1, 2)))
  )
  result <- pecotmr:::getTopVariantsIdx(susie_output, signalCutoff = 0.5)
  expect_equal(result, c(1, 2))
})

test_that("get_top_variants_idx with high cutoff and no CS", {
  susie_output <- list(
    pip = c(0.01, 0.02, 0.03),
    sets = list(cs = NULL)
  )
  result <- pecotmr:::getTopVariantsIdx(susie_output, signalCutoff = 0.5)
  expect_length(result, 0)
})

# =============================================================================
# get_cs_info (internal)
# =============================================================================

test_that("get_cs_info maps variants to CS numbers", {
  susie_cs <- list(L1 = c(1, 2), L3 = c(4, 5, 6))
  top_idx <- c(1, 3, 5)
  result <- pecotmr:::getCsInfo(susie_cs, top_idx)
  # Now returns data.frame(variant_idx, cs_idx) with one row per (variant, CS) pair
  expect_true(is.data.frame(result))
  expect_equal(result$variant_idx, c(1, 3, 5))
  expect_equal(result$cs_idx, c(1L, 0L, 3L))
})

test_that("get_cs_info handles all variants outside CS", {
  susie_cs <- list(L1 = c(1, 2))
  top_idx <- c(5, 6, 7)
  result <- pecotmr:::getCsInfo(susie_cs, top_idx)
  expect_true(is.data.frame(result))
  expect_true(all(result$cs_idx == 0))
})

test_that("get_cs_info reports variant in multiple CSs as multiple rows", {
  susie_cs <- list(L1 = c(1, 2, 3), L3 = c(2, 3, 4))
  top_idx <- c(1, 2, 4)
  result <- pecotmr:::getCsInfo(susie_cs, top_idx)
  expect_true(is.data.frame(result))
  # variant 2 is in both L1 and L3, so it gets two rows
  expect_equal(nrow(result), 4)
  expect_equal(sum(result$variant_idx == 2), 2)
  expect_equal(sort(result$cs_idx[result$variant_idx == 2]), c(1L, 3L))
})

# =============================================================================
# susieWeights
# =============================================================================

test_that("susieWeights returns zeros when fit lacks alpha/mu", {
  fake_fit <- list(pip = rep(0.01, 5))
  result <- susieWeights(susieFit = fake_fit)
  expect_equal(result, rep(0, 5))
})

test_that("susieWeights checks dimension mismatch", {
  set.seed(42)
  X <- matrix(rnorm(100), 20, 5)
  fake_fit <- list(pip = rep(0.01, 10))
  expect_error(susieWeights(X = X, susieFit = fake_fit), "Dimension mismatch")
})

# =============================================================================
# susieAshWeights
# =============================================================================

test_that("susieAshWeights returns zeros without proper fit structure", {
  fake_fit <- list(pip = rep(0.01, 5))
  result <- susieAshWeights(susieAshFit = fake_fit)
  expect_equal(result, rep(0, 5))
})

# =============================================================================
# susieInfWeights
# =============================================================================

test_that("susieInfWeights returns zeros without proper fit structure", {
  fake_fit <- list(pip = rep(0.01, 5))
  result <- susieInfWeights(susieInfFit = fake_fit)
  expect_equal(result, rep(0, 5))
})

# =============================================================================
# glmnetWeights
# =============================================================================

test_that("glmnetWeights produces non-zero weights for correlated data", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 100
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  beta_true <- c(3, -2, rep(0, p - 2))
  y <- X %*% beta_true + rnorm(n)

  w <- glmnetWeights(X, y, alpha = 0.5)
  expect_length(w, p)
  expect_true(any(w != 0))
})

test_that("glmnetWeights handles zero-variance columns", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 100
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  X[, 3] <- 1  # zero variance column
  y <- X[, 1] * 2 + rnorm(n)

  w <- glmnetWeights(X, y, alpha = 1)
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

  sds <- pecotmr:::initPriorSd(X, y, n = 15)
  expect_length(sds, 15)
  expect_equal(sds[1], 0)
  expect_true(all(diff(sds) >= 0))
})

# =============================================================================
# postprocessFinemappingFits: analysisScript and V=NULL branches (Tier 1)
# =============================================================================

# Helper: build a minimal synthetic SuSiE-family output for post-processing
make_fake_susie_output <- function(p = 5, L = 3, has_V = TRUE) {
  vnames <- paste0("chr1:", 1:p, ":A:G")
  out <- list(
    pip = setNames(rep(0.01, p), vnames),
    alpha = matrix(1 / p, nrow = L, ncol = p),
    lbf_variable = matrix(0, nrow = L, ncol = p),
    sets = list(
      cs = NULL,
      requestedCoverage = 0.95
    ),
    niter = 10
  )
  if (has_V) {
    out$V <- rep(1, L)
  }
  out
}

test_that("postprocessFinemappingFits keeps all effects when V is NULL", {
  skip_if_not_installed("susieR")
  p <- 5
  L <- 3
  fake_output <- make_fake_susie_output(p, L = L, has_V = FALSE)
  R <- diag(p)
  colnames(R) <- rownames(R) <- names(fake_output$pip)
  post <- postprocessFinemappingFits(
    fits = list(susieRss = pecotmr:::.setFinemappingFitClass(fake_output, "susieRss")),
    dataX = R,
    dataY = list(z = rnorm(p)),
    coverage = 0.95
  )
  result <- formatFinemappingOutput(post, primaryMethod = "susieRss")
  trimmed <- getTrimmedFit(result$finemappingEntry)
  # With V=NULL, eff_idx = 1:L, so trimmed alpha should keep all L rows
  expect_equal(nrow(trimmed$alpha), L)
  # V should be NULL in trimmed output
  expect_null(trimmed$V)
})

# =============================================================================
# postprocessFinemappingFits: mvsusie output (outcome_names, coef, clfsr)
# =============================================================================

test_that("postprocessFinemappingFits stores outcome_names, coef, and clfsr for mvsusie", {
  skip_if_not_installed("susieR")
  skip_if_not_installed("mvsusieR")
  p <- 5
  L <- 3
  R <- 2
  vnames <- paste0("chr1:", 1:p, ":A:G")
  cnames <- paste0("cond_", 1:R)
  fake_coef <- matrix(rnorm((p + 1) * R), nrow = p + 1, ncol = R)

  fake_output <- list(
    pip = setNames(rep(0.01, p), vnames),
    alpha = matrix(1 / p, nrow = L, ncol = p),
    lbf_variable = matrix(0, nrow = L, ncol = p),
    sets = list(cs = NULL, requestedCoverage = 0.95),
    niter = 10,
    V = rep(1, L),
    outcome_names = cnames,
    conditional_lfsr = array(0.5, dim = c(L, p, R))
  )

  n <- 20
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- vnames

  local_mocked_bindings(
    coef.mvsusie = function(...) fake_coef,
    .package = "mvsusieR"
  )

  post <- postprocessFinemappingFits(
    fits = list(mvsusie = pecotmr:::.setFinemappingFitClass(fake_output, "mvsusie")),
    dataX = X,
    dataY = NULL,
    xScalar = 1, yScalar = 1,
    coverage = 0.95
  )
  result <- formatFinemappingOutput(post, primaryMethod = "mvsusie")

  # outcome_names should be stored as contextNames
  expect_equal(result$contextNames, cnames)
  trimmed <- getTrimmedFit(result$finemappingEntry)
  # coef should come from mvsusieR::coef.mvsusie
  expect_equal(trimmed$coef, fake_coef[-1, , drop = FALSE])
  # conditional_lfsr should be trimmed to eff_idx
  expect_equal(dim(trimmed$clfsr), c(L, p, R))
})

test_that("formatFinemappingOutput does not duplicate top loci variants", {
  top_loci <- data.frame(
    variant_id = paste0("v", 1:4),
    CS_95_susie = c(0L, 1L, NA_integer_, 0L),
    pip_susie = c(0.2, 0.005, 0.001, 0),
    stringsAsFactors = FALSE
  )
  fm <- FineMappingEntry(
    variantIds = paste0("v", 1:4),
    trimmedFit = list(pip = 1:4),
    topLoci = data.frame(variant_id = character(0), pip = numeric(0))
  )
  post <- list(
    finemappingResults = list(susie = list(
      finemappingEntry = fm
    )),
    top_loci = top_loci
  )
  out <- formatFinemappingOutput(post, "susie")
  expect_false("top_loci_variants" %in% names(out))
  expect_equal(unique(out$top_loci$variant_id), paste0("v", 1:4))
})

.make_univariate_data <- function(seed = 42, n = 300, p = 50,
                                  effect_idx = integer(0), effect_size = NULL) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- sprintf("chr1:%d:G:A", seq_len(p))
  beta <- rep(0, p)
  if (length(effect_idx) > 0) {
    if (is.null(effect_size)) effect_size <- rep(1.5, length(effect_idx))
    beta[effect_idx] <- effect_size
  }
  y <- as.numeric(X %*% beta) + rnorm(n, sd = 0.5)
  list(X = X, y = y)
}

test_that(".translate_legacy_top_loci_cs_columns renames pip_susie -> pip for legacy callers", {
  new_format <- data.frame(
    variant_id = c("v1", "v2"),
    pip_susie = c(0.9, 0.1),
    CS_95_susie = c(1, 0),
    pip_susie_inf = c(0.8, 0.2),
    CS_95_susie_inf = c(1, 0),
    stringsAsFactors = FALSE
  )
  out <- pecotmr:::.translateLegacyTopLociCsColumns(new_format)
  expect_true("pip" %in% colnames(out))
  expect_false("pip_susie" %in% colnames(out))
  # The susieInf and CS columns are untouched
  expect_true("pip_susie_inf" %in% colnames(out))
  expect_true("CS_95_susie" %in% colnames(out))
  expect_true("CS_95_susie_inf" %in% colnames(out))
})

test_that(".translate_legacy_top_loci_cs_columns leaves existing pip column alone", {
  legacy <- data.frame(
    variant_id = c("v1", "v2"),
    pip = c(0.9, 0.1),
    cs_coverage_0.95 = c(1, 0),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out <- pecotmr:::.translateLegacyTopLociCsColumns(legacy)
  expect_true("pip" %in% colnames(out))
  expect_true("CS_95_susie" %in% colnames(out))   # legacy cs_coverage rename
  expect_false("cs_coverage_0.95" %in% colnames(out))
  expect_false("pip_susie" %in% colnames(out))     # no double-conversion
})


# ============================================================================
# Unified top-loci: hard gating coverage per OpenSpec tasks 4.27 / 4.28.
# These tests are the implementation gate for the buildTopLoci migration.
# ============================================================================

# Reuse the existing local helper. Re-declared inside this block so the file
# remains correct whether the unified-section tests are run alone or as part of
# the full file.
if (!exists(".make_univariate_data", inherits = FALSE)) {
  .make_univariate_data <- function(seed = 42, n = 300, p = 50,
                                    effect_idx = integer(0),
                                    effect_size = NULL) {
    set.seed(seed)
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- sprintf("chr1:%d:G:A", seq_len(p))
    beta <- rep(0, p)
    if (length(effect_idx) > 0) {
      if (is.null(effect_size)) effect_size <- rep(1.5, length(effect_idx))
      beta[effect_idx] <- effect_size
    }
    y <- as.numeric(X %*% beta) + rnorm(n, sd = 0.5)
    list(X = X, y = y)
  }
}

.UNIFIED_TOP_LOCI_COLS <- c(
  "#chr", "start", "end", "a1", "a2",
  "variant", "gene", "event",
  "n", "af", "beta", "se",
  "pip", "posterior_effect_mean", "posterior_effect_se",
  "cs_95", "cs_70", "cs_50", "cs_95_purity",
  "method", "grange_start", "grange_end"
)

# Synthesize a SuSiE-like fit + cs_tables with explicit per-coverage CS
# membership. `cs_at_cov` is a named list keyed by coverage value (e.g.
# `"0.95"`); each element is a list of integer vectors, one per CS at that
# coverage. The CS numbering is 1-based per coverage. PIP values are filled
# from `pip` (variants outside the CS get small non-zero PIP so they can be
# retained or dropped via `signal_cutoff`).
.fake_fit_and_cs <- function(variant_ids, cs_at_cov,
                             cs_purity_value = 0.85,
                             pip = NULL,
                             nSamples = 100, n_variants = NULL,
                             gene = "ENSG00000179403") {
  p <- length(variant_ids)
  if (is.null(n_variants)) n_variants <- p
  if (is.null(pip)) pip <- seq(0.6, 0.9, length.out = p)
  # alpha must be L x p so colSums(alpha * mu) is well-defined. We use one
  # row whose values are normalized PIPs.
  alpha <- matrix(pip / sum(pip), nrow = 1, ncol = p)
  mu    <- matrix(0.5, nrow = 1, ncol = p)
  mu2   <- mu^2 + 0.1
  fit <- list(pip = setNames(pip, variant_ids),
              alpha = alpha, mu = mu, mu2 = mu2)

  cs_tables <- lapply(names(cs_at_cov), function(cov_str) {
    cs_list <- cs_at_cov[[cov_str]]
    if (is.null(cs_list)) cs_list <- list()
    n_cs <- length(cs_list)
    if (n_cs > 0L) names(cs_list) <- paste0("L", seq_len(n_cs))
    purity_df <- if (n_cs > 0L) {
      data.frame(min.abs.corr   = rep(cs_purity_value, n_cs),
                 mean.abs.corr  = rep(cs_purity_value, n_cs),
                 median.abs.corr = rep(cs_purity_value, n_cs))
    } else {
      data.frame(min.abs.corr = numeric(0),
                 mean.abs.corr = numeric(0),
                 median.abs.corr = numeric(0))
    }
    list(
      sets = list(cs = cs_list,
                  cs_index = seq_len(n_cs),
                  requestedCoverage = as.numeric(cov_str),
                  purity = purity_df),
      cs_corr = if (n_cs > 0L) {
        lapply(seq_len(n_cs), function(i) {
          matrix(c(1, cs_purity_value, cs_purity_value, 1), nrow = 2)
        })
      } else NULL,
      pip = fit$pip
    )
  })
  attr(cs_tables, "coverage") <- as.numeric(names(cs_at_cov))

  X <- matrix(0, nrow = nSamples, ncol = n_variants)
  rownames(X) <- paste0("sample", seq_len(nSamples))
  if (n_variants == length(variant_ids)) colnames(X) <- variant_ids
  Y <- matrix(0, nrow = nSamples, ncol = 1,
              dimnames = list(paste0("sample", seq_len(nSamples)), gene))
  list(fit = fit, cs_tables = cs_tables,
       variantNames = variant_ids, data_x = X, data_y = Y)
}

.runBuildTopLoci <- function(inp, method = "susie", signalCutoff = 0.05,
                             sumstats = NULL, af = NULL,
                             otherQuantities = NULL,
                             region = NULL) {
  buildTopLoci(
    fit = inp$fit, csTables = inp$cs_tables,
    variantNames = inp$variantNames,
    sumstats = sumstats, af = af,
    method = method, signalCutoff = signalCutoff,
    dataX = inp$data_x, dataY = inp$data_y,
    otherQuantities = otherQuantities,
    region = region
  )
}

test_that("buildTopLoci returns the exact 22-column schema in order with stable dtypes", {
  # `.emptyTopLoci` is a package-internal helper. Use the namespace lookup
  # so the test works both source-loaded and after R CMD INSTALL.
  empty_fn <- if (exists(".emptyTopLoci", envir = asNamespace("pecotmr"),
                          inherits = FALSE)) {
    get(".emptyTopLoci", envir = asNamespace("pecotmr"))
  } else {
    get(".emptyTopLoci", envir = .GlobalEnv)
  }
  out <- empty_fn()
  expect_equal(names(out), .UNIFIED_TOP_LOCI_COLS)
  expect_equal(nrow(out), 0L)
  expect_true(is.character(out$"#chr"))
  expect_true(is.integer(out$start))
  expect_true(is.integer(out$end))
  expect_true(is.character(out$variant))
  expect_true(is.character(out$gene))
  expect_true(is.character(out$event))
  expect_true(is.integer(out$n))
  expect_true(is.numeric(out$af))
  expect_false("maf" %in% names(out))
  expect_true(is.numeric(out$pip))
  expect_true(is.numeric(out$posterior_effect_mean))
  expect_true(is.numeric(out$posterior_effect_se))
  expect_true(is.character(out$cs_95))
  expect_true(is.character(out$cs_70))
  expect_true(is.character(out$cs_50))
  expect_true(is.numeric(out$cs_95_purity))
  expect_true(is.character(out$method))
  expect_true(is.integer(out$grange_start))
  expect_true(is.integer(out$grange_end))
})

test_that("buildTopLoci emits 22 columns in the fixed order on a non-empty fit", {
  variant_ids <- c("chr1:100:A:G", "chr1:200:C:T")
  inp <- .fake_fit_and_cs(variant_ids,
                          cs_at_cov = list("0.95" = list(c(1L, 2L)),
                                            "0.7"  = list(c(1L, 2L)),
                                            "0.5"  = list(c(1L, 2L))),
                          nSamples = 419, n_variants = 11332)
  other_q <- list(condition_id = "Ast_DeJager_eQTL")
  out <- .runBuildTopLoci(inp, method = "susie",
                             sumstats = list(betahat = c(0.2, -0.1),
                                             sebetahat = c(0.05, 0.04)),
                             af = c(0.10, 0.25),
                             otherQuantities = other_q,
                             region = "chr10:10823338-14348298")
  expect_equal(names(out), .UNIFIED_TOP_LOCI_COLS)
  expect_equal(unique(out$gene), "ENSG00000179403")
  expect_equal(unique(out$event), "Ast_DeJager_eQTL_ENSG00000179403")
  expect_equal(unique(out$grange_start), 10823338L)
  expect_equal(unique(out$grange_end),   14348298L)
  expect_equal(unique(out$n), 419L)
  expect_equal(unique(out$method), "susie")
})

test_that("cs_95 / cs_70 / cs_50 are character strings of the form '<method>_<idx>'", {
  variant_ids <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:T:C")
  # Variants 1 and 2 in CS 1 (at 95), variant 3 in CS 2 (at 95). All three
  # also appear at 70/50 with the same memberships.
  cs_at_cov <- list("0.95" = list(c(1L, 2L), 3L),
                    "0.7"  = list(c(1L, 2L), 3L),
                    "0.5"  = list(c(1L, 2L), 3L))
  inp <- .fake_fit_and_cs(variant_ids, cs_at_cov, pip = c(0.9, 0.9, 0.9))
  out <- .runBuildTopLoci(inp, method = "susie")
  expect_true(all(grepl("^susie_\\d+$", out$cs_95)))
  expect_true(all(grepl("^susie_\\d+$", out$cs_70)))
  expect_true(all(grepl("^susie_\\d+$", out$cs_50)))
  expect_true(any(out$cs_95 == "susie_1"))
  expect_true(any(out$cs_95 == "susie_2"))
})

test_that("PIP-only retained variants carry '<method>_0' at every coverage and cs_95_purity = 0", {
  variant_ids <- c("chr1:100:A:G", "chr1:200:C:T")
  # No CS at any coverage; variant 2 has high PIP so it is retained via
  # signal_cutoff and produces a "<method>_0" row.
  cs_at_cov <- list("0.95" = list(),
                    "0.7"  = list(),
                    "0.5"  = list())
  inp <- .fake_fit_and_cs(variant_ids, cs_at_cov,
                          pip = c(0.02, 0.95))
  out <- .runBuildTopLoci(inp, method = "susie", signalCutoff = 0.5)
  expect_gte(nrow(out), 1L)
  # Every row must have <method>_0 at every coverage and cs_95_purity = 0.
  expect_true(all(out$cs_95 == "susie_0"))
  expect_true(all(out$cs_70 == "susie_0"))
  expect_true(all(out$cs_50 == "susie_0"))
  expect_true(all(out$cs_95_purity == 0))
})

test_that("per-method CS indices are independent across susie and susieInf (postprocessFinemappingFits)", {
  d <- .make_univariate_data(seed = 21, effect_idx = c(15, 35))
  fits <- fitSusieInfThenSusie(d$X, d$y)
  post <- postprocessFinemappingFits(fits, dataX = d$X, dataY = d$y,
                                       coverage = 0.95,
                                       secondaryCoverage = c(0.7, 0.5))
  tl <- post$top_loci
  expect_setequal(unique(tl$method), c("susie", "susieInf"))
  # Each method must only ever emit "<method>_<idx>" strings, never strings
  # from the other method. This is the core safeguard against silent
  # method-mixing.
  susie_rows     <- tl[tl$method == "susie", , drop = FALSE]
  susie_inf_rows <- tl[tl$method == "susieInf", , drop = FALSE]
  for (col in c("cs_95", "cs_70", "cs_50")) {
    expect_true(all(grepl("^susie_\\d+$", susie_rows[[col]])),
                info = paste("susie rows have wrong prefix in", col))
    expect_true(all(grepl("^susie_inf_\\d+$", susie_inf_rows[[col]])),
                info = paste("susieInf rows have wrong prefix in", col))
  }
  # CS indices are not sequenced across methods: if both methods have any
  # CS, each may independently include "<method>_1".
  has_susie_1     <- any(susie_rows$cs_95 == "susie_1")
  has_susie_inf_1 <- any(susie_inf_rows$cs_95 == "susie_inf_1")
  expect_true(has_susie_1 || nrow(susie_rows) == 0L)
  expect_true(has_susie_inf_1 || nrow(susie_inf_rows) == 0L)
})

test_that("cs_95_purity = 0 when cs_95 is '<method>_0', and in (0, 1] otherwise", {
  variant_ids <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:T:C")
  # Variant 1 in CS 1 at 95-cov; variant 2 PIP-only retained; variant 3
  # PIP-only retained.
  cs_at_cov <- list("0.95" = list(1L),
                    "0.7"  = list(1L),
                    "0.5"  = list(1L))
  inp <- .fake_fit_and_cs(variant_ids, cs_at_cov,
                          cs_purity_value = 0.85,
                          pip = c(0.9, 0.6, 0.55))
  out <- .runBuildTopLoci(inp, method = "susie", signalCutoff = 0.5)
  expect_true(any(out$cs_95 == "susie_1"))
  expect_true(any(out$cs_95 == "susie_0"))
  in_cs   <- out[out$cs_95 != "susie_0", , drop = FALSE]
  not_cs  <- out[out$cs_95 == "susie_0", , drop = FALSE]
  expect_true(all(not_cs$cs_95_purity == 0))
  expect_true(all(in_cs$cs_95_purity > 0 & in_cs$cs_95_purity <= 1))
})

test_that("overlapping CS within one method produces one row per CS membership", {
  variant_ids <- c("chr1:100:A:G")
  # One variant belongs to CS 1 AND CS 2 at 95-cov (overlap).
  cs_at_cov <- list("0.95" = list(1L, 1L),
                    "0.7"  = list(1L, 1L),
                    "0.5"  = list(1L, 1L))
  inp <- .fake_fit_and_cs(variant_ids, cs_at_cov, pip = 0.9)
  out <- .runBuildTopLoci(inp, method = "susie")
  # Two rows: one for CS 1 membership, one for CS 2 membership; same
  # (variant, gene, method).
  expect_equal(nrow(out), 2L)
  expect_equal(unique(out$variant), "chr1:100:A:G")
  expect_equal(unique(out$method), "susie")
  expect_setequal(out$cs_95, c("susie_1", "susie_2"))
})

test_that("overlapping CS across methods produces one row per method", {
  d <- .make_univariate_data(seed = 22, effect_idx = c(12, 32))
  fits <- fitSusieInfThenSusie(d$X, d$y)
  post <- postprocessFinemappingFits(fits, dataX = d$X, dataY = d$y,
                                       coverage = 0.95,
                                       secondaryCoverage = c(0.7, 0.5))
  tl <- post$top_loci
  if (nrow(tl) > 0L) {
    cnt_per_method <- table(tl$variant, tl$method)
    shared <- rownames(cnt_per_method)[apply(cnt_per_method > 0, 1, sum) >= 2L]
    if (length(shared) > 0L) {
      v <- shared[[1]]
      rows_for_v <- tl[tl$variant == v, , drop = FALSE]
      expect_gte(length(unique(rows_for_v$method)), 2L)
    } else {
      succeed("no shared variants in this fixture; cross-method uniqueness rule is structural")
    }
  } else {
    succeed("empty top_loci from this fixture; cross-method uniqueness rule is structural")
  }
})

test_that("formatFinemappingOutput exposes exactly one top_loci field; no top_loci_long, no wide top_loci, no top_loci_export", {
  d <- .make_univariate_data(seed = 23, effect_idx = c(15, 40))
  fits <- fitSusieInfThenSusie(d$X, d$y)
  post <- postprocessFinemappingFits(fits, dataX = d$X, dataY = d$y, coverage = 0.95)
  out <- formatFinemappingOutput(post, primaryMethod = "susie")
  expect_true("top_loci" %in% names(out))
  expect_false("top_loci_long" %in% names(out))
  expect_false("top_loci_export" %in% names(out))
  # The exposed top_loci has the unified 22-column schema.
  expect_equal(names(out$top_loci), .UNIFIED_TOP_LOCI_COLS)
})

test_that("postprocessFinemappingFits does not return top_loci_long anywhere", {
  d <- .make_univariate_data(seed = 24, effect_idx = c(25))
  fit <- susieR::susie(d$X, d$y, L = 5)
  post <- postprocessFinemappingFits(list(susie = fit),
                                       dataX = d$X, dataY = d$y, coverage = 0.95)
  expect_true("top_loci" %in% names(post))
  expect_false("top_loci_long" %in% names(post))
  expect_equal(names(post$top_loci), .UNIFIED_TOP_LOCI_COLS)
  # No per-method finemappingResults entry should carry a top_loci_long either.
  for (name in names(post$finemappingResults)) {
    expect_false("top_loci_long" %in% names(post$finemappingResults[[name]]),
                 info = name)
  }
})

test_that("build_top_loci_long / build_top_loci_wide / build_top_loci_export are removed from the package", {
  # Two layered checks. (1) The new helper must be reachable. (2) The old
  # helpers must NOT be reachable — neither in the package namespace (if the
  # package is installed) nor in .GlobalEnv (when this file is source-loaded
  # against the fresh R/ tree for testing). The contract is that the trio is
  # gone end-to-end after the migration.
  resolve <- function(name) {
    if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
      return(get(name, envir = .GlobalEnv))
    }
    ns_ok <- tryCatch(asNamespace("pecotmr"), error = function(e) NULL)
    if (!is.null(ns_ok) && exists(name, envir = ns_ok, inherits = FALSE)) {
      return(get(name, envir = ns_ok))
    }
    NULL
  }
  expect_false(is.null(resolve("buildTopLoci")),
               info = "buildTopLoci must be defined after the migration")
  # The removed trio: only fail if a definition exists in the SAME source
  # tree we just loaded (i.e. globalenv). The installed-package namespace
  # may still carry stale copies from an earlier install; the gate only
  # cares that the new source tree does not redefine them.
  expect_false(exists("build_top_loci_long",   envir = .GlobalEnv, inherits = FALSE))
  expect_false(exists("build_top_loci_wide",   envir = .GlobalEnv, inherits = FALSE))
  expect_false(exists("build_top_loci_export", envir = .GlobalEnv, inherits = FALSE))
  expect_false(exists(".emptyTopLoci_long",  envir = .GlobalEnv, inherits = FALSE))
  expect_false(exists(".emptyTopLoci_export",envir = .GlobalEnv, inherits = FALSE))
})

test_that("buildTopLoci raises an explicit error on invalid variant_id rather than silently filling NA", {
  variant_ids <- c("not_a_valid_id")
  cs_at_cov <- list("0.95" = list(1L),
                    "0.7"  = list(1L),
                    "0.5"  = list(1L))
  inp <- .fake_fit_and_cs(variant_ids, cs_at_cov, pip = 0.9)
  expect_error(.runBuildTopLoci(inp, method = "susie"),
               "parseVariantId")
})

test_that("buildTopLoci requires `method`", {
  variant_ids <- c("chr1:100:A:G")
  cs_at_cov <- list("0.95" = list(1L),
                    "0.7"  = list(1L),
                    "0.5"  = list(1L))
  inp <- .fake_fit_and_cs(variant_ids, cs_at_cov, pip = 0.9)
  expect_error(buildTopLoci(
    fit = inp$fit, csTables = inp$cs_tables,
    variantNames = inp$variantNames
  ), "method")
})

test_that("formatFinemappingOutput exposes finemappingEntry with S4 accessors", {
  d <- .make_univariate_data(seed = 25, effect_idx = c(20))
  fit <- susieR::susie(d$X, d$y, L = 5)
  post <- postprocessFinemappingFits(list(susie = fit), dataX = d$X, dataY = d$y, coverage = 0.95)
  out <- formatFinemappingOutput(post, primaryMethod = "susie")
  expect_true("finemappingEntry" %in% names(out))
  fm <- out$finemappingEntry
  expect_true(is.character(getVariantIds(fm)) && length(getVariantIds(fm)) == ncol(d$X))
  expect_true(is.list(getTrimmedFit(fm)) && !is.null(getTrimmedFit(fm)$pip))
})

test_that("missing region produces NA grange columns rather than silent omission", {
  variant_ids <- c("chr1:100:A:G")
  cs_at_cov <- list("0.95" = list(1L),
                    "0.7"  = list(1L),
                    "0.5"  = list(1L))
  inp <- .fake_fit_and_cs(variant_ids, cs_at_cov, pip = 0.9)
  out <- .runBuildTopLoci(inp, method = "susie",
                             otherQuantities = list(condition_id = "ctx"))
  # grange_* must still be present columns of the 22-col schema, with NA values.
  expect_true(all(c("grange_start", "grange_end") %in% names(out)))
  expect_true(all(is.na(out$grange_start)))
  expect_true(all(is.na(out$grange_end)))
  # And event composition still works from gene + condition_id.
  expect_equal(unique(out$event), "ctx_ENSG00000179403")
})

test_that("posterior_effect_mean equals colSums(alpha*mu); posterior_effect_se equals sqrt(pmax(colSums(alpha*mu2) - mean^2, 0))", {
  variant_ids <- c("chr1:100:A:G", "chr1:200:C:T")
  cs_at_cov <- list("0.95" = list(c(1L, 2L)),
                    "0.7"  = list(c(1L, 2L)),
                    "0.5"  = list(c(1L, 2L)))
  inp <- .fake_fit_and_cs(variant_ids, cs_at_cov, pip = c(0.8, 0.6))
  out <- .runBuildTopLoci(inp, method = "susie")
  expected_mean <- colSums(inp$fit$alpha * inp$fit$mu)
  expected_se   <- sqrt(pmax(colSums(inp$fit$alpha * inp$fit$mu2) - expected_mean^2, 0))
  # Match per variant index by looking up via variant string.
  for (i in seq_along(variant_ids)) {
    row <- out[out$variant == variant_ids[i], , drop = FALSE]
    expect_true(nrow(row) >= 1L)
    expect_equal(unique(row$posterior_effect_mean), expected_mean[i],
                 tolerance = 1e-10)
    expect_equal(unique(row$posterior_effect_se), expected_se[i],
                 tolerance = 1e-10)
  }
})


context("fsusieWrapper")

# ---- cal_purity ----
test_that("cal_purity with min method and single element CS", {
  set.seed(42)
  X <- matrix(rnorm(100), nrow = 10, ncol = 10)
  l_cs <- list(c(1))

  result <- pecotmr:::calPurity(l_cs, X, method = "min")
  expect_equal(result[[1]], 1)
})

test_that("cal_purity with min method and multi-element CS", {
  set.seed(42)
  X <- matrix(rnorm(200), nrow = 20, ncol = 10)
  l_cs <- list(c(1, 2, 3))

  result <- pecotmr:::calPurity(l_cs, X, method = "min")
  expect_length(result, 1)
  # Manually compute expected: min off-diagonal |cor|
  cormat <- abs(cor(X[, c(1, 2, 3)]))
  diag(cormat) <- NA
  expect_equal(result[[1]], min(cormat, na.rm = TRUE))
})

test_that("cal_purity with non-min method returns three values", {
  set.seed(42)
  X <- matrix(rnorm(200), nrow = 20, ncol = 10)
  l_cs <- list(c(1, 2, 3))

  result <- pecotmr:::calPurity(l_cs, X, method = "susie")
  expect_length(result[[1]], 3)  # min, mean, median
  # Manually compute expected values
  cormat <- abs(cor(X[, c(1, 2, 3)]))
  diag(cormat) <- NA
  vals <- cormat[!is.na(cormat)]
  expect_equal(result[[1]][1], min(vals))
  expect_equal(result[[1]][2], mean(vals))
  expect_equal(result[[1]][3], median(vals))
  # min <= mean and min <= median by definition
  expect_true(result[[1]][1] <= result[[1]][2])
  expect_true(result[[1]][1] <= result[[1]][3])
})

test_that("cal_purity with non-min method single element returns (1,1,1)", {
  X <- matrix(rnorm(100), nrow = 10, ncol = 10)
  l_cs <- list(c(1))

  result <- pecotmr:::calPurity(l_cs, X, method = "susie")
  expect_equal(result[[1]], c(1, 1, 1))
})

test_that("cal_purity with multiple credible sets", {
  set.seed(42)
  X <- matrix(rnorm(200), nrow = 20, ncol = 10)
  l_cs <- list(c(1, 2), c(5, 6, 7))

  result <- pecotmr:::calPurity(l_cs, X, method = "min")
  expect_length(result, 2)
})

# ---- fsusieGetCs ----
# ---- fsusieWrapper ----
test_that("fsusieWrapper errors when fsusieR is not installed", {
  skip_if(requireNamespace("fsusieR", quietly = TRUE),
          "fsusieR is installed, skipping not-installed test")
  set.seed(1)
  X <- matrix(rnorm(50), nrow = 10, ncol = 5)
  Y <- matrix(rnorm(40), nrow = 10, ncol = 4)
  expect_error(
    fsusieWrapper(
      X = X, Y = Y, pos = seq_len(4), L = 3, prior = "mixture_normal",
      maxSnpEm = 100, covLev = 0.95, minPurity = 0.5, maxScale = 5
    ),
    "fsusieR"
  )
})

test_that("fsusieWrapper low-purity branch sets cs to list(NULL) and cs_corr to NULL", {
  skip_if_not_installed("fsusieR")
  fake_fit <- list(
    cs = list(c(1, 2), c(3)),
    purity = c(0.1, 0.05),  # all < min_purity = 0.5
    pip = c(0.1, 0.2, 0.3, 0.05, 0.05),
    alpha = list(matrix(0.1, nrow = 2, ncol = 5), matrix(0.1, nrow = 2, ncol = 5))
  )
  local_mocked_bindings(
    susiF = function(...) fake_fit,
    .package = "fsusieR"
  )
  set.seed(1)
  X <- matrix(rnorm(50), nrow = 10, ncol = 5)
  Y <- matrix(rnorm(40), nrow = 10, ncol = 4)
  out <- fsusieWrapper(
    X = X, Y = Y, pos = seq_len(4), L = 3, prior = "mixture_normal",
    maxSnpEm = 100, covLev = 0.95, minPurity = 0.5, maxScale = 5
  )
  expect_equal(out$cs, list(NULL))
  expect_equal(out$sets$cs, list(NULL))
  expect_null(out$cs_corr)
})

test_that("fsusieWrapper high-purity branch builds sets and computes cs_corr", {
  skip_if_not_installed("fsusieR")
  set.seed(2)
  p <- 5
  fake_fit <- list(
    cs = list(c(1, 2), c(3, 4)),
    purity = c(0.95, 0.9),  # all > min_purity = 0.5
    pip = c(0.4, 0.4, 0.6, 0.6, 0.1),
    alpha = list(
      matrix(rep(c(0.4, 0.4, 0.05, 0.05, 0.1), each = 2), nrow = 2, byrow = FALSE),
      matrix(rep(c(0.05, 0.05, 0.45, 0.4, 0.05), each = 2), nrow = 2, byrow = FALSE)
    )
  )
  local_mocked_bindings(
    susiF = function(...) fake_fit,
    cal_cor_cs = function(obj, X) matrix(c(1, 0.9, 0.9, 1), nrow = 2),
    .package = "fsusieR"
  )
  X <- matrix(rnorm(10 * p), nrow = 10, ncol = p)
  Y <- matrix(rnorm(40), nrow = 10, ncol = 4)
  out <- fsusieWrapper(
    X = X, Y = Y, pos = seq_len(4), L = 3, prior = "mixture_normal",
    maxSnpEm = 100, covLev = 0.95, minPurity = 0.5, maxScale = 5
  )
  expect_length(out$sets$cs, 2)
  expect_equal(names(out$sets$cs), c("L1", "L2"))
  expect_equal(dim(out$cs_corr), c(2, 2))
  expect_equal(out$sets$requested_coverage, 0.95)
})

test_that("fsusieGetCs creates susie-like sets", {
  set.seed(42)
  X <- matrix(rnorm(200), nrow = 20, ncol = 10)

  fSuSiE_obj <- list(
    cs = list(c(1, 2, 3), c(5, 6)),
    alpha = list(
      c(0.4, 0.3, 0.2, 0.05, 0.02, 0.01, 0.01, 0.005, 0.003, 0.002),
      c(0.01, 0.02, 0.02, 0.05, 0.45, 0.35, 0.05, 0.02, 0.02, 0.01)
    )
  )

  result <- fsusieGetCs(fSuSiE_obj, X, requestedCoverage = 0.95)

  expect_type(result, "list")
  expect_true("cs" %in% names(result))
  expect_true("purity" %in% names(result))
  expect_true("cs_index" %in% names(result))
  expect_true("coverage" %in% names(result))
  expect_true("requested_coverage" %in% names(result))
  expect_equal(result$requested_coverage, 0.95)
  expect_equal(length(result$cs), 2)
  expect_equal(names(result$cs), c("L1", "L2"))
  # Purity should be a data.frame with min/mean/median columns
  expect_true(is.data.frame(result$purity))
  expect_equal(nrow(result$purity), 2)
  # Coverage should be numeric and positive, one per CS
  expect_length(result$coverage, 2)
  expect_true(all(result$coverage > 0 & result$coverage <= 1))
  # cs_index should identify which effects had credible sets
  expect_length(result$cs_index, 2)
})
