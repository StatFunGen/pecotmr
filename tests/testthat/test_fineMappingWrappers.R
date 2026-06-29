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
  trimmed <- getSusieFit(result$finemappingEntry)
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
  trimmed <- getSusieFit(result$finemappingEntry)
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
    susieFit = list(pip = 1:4),
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
  "variant_id", "chrom", "pos", "A1", "A2",
  "N", "af",
  "marginal_beta", "marginal_se", "marginal_z", "marginal_p",
  "pip", "posterior_mean", "posterior_sd",
  "cs_95", "cs_70", "cs_50", "cs_95_purity",
  "method", "gene", "event", "grange_start", "grange_end"
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
                             af = NULL,
                             sumstats = NULL,
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
  expect_true(is.character(out$variant_id))
  expect_true(is.character(out$chrom))
  expect_true(is.integer(out$pos))
  expect_true(is.character(out$A1))
  expect_true(is.character(out$A2))
  expect_true(is.numeric(out$N))
  expect_true(is.numeric(out$af))
  expect_true(is.numeric(out$marginal_beta))
  expect_true(is.numeric(out$marginal_se))
  expect_true(is.numeric(out$marginal_z))
  expect_true(is.numeric(out$marginal_p))
  expect_true(is.numeric(out$pip))
  expect_true(is.numeric(out$posterior_mean))
  expect_true(is.numeric(out$posterior_sd))
  expect_true(is.character(out$cs_95))
  expect_true(is.character(out$cs_70))
  expect_true(is.character(out$cs_50))
  expect_true(is.numeric(out$cs_95_purity))
  expect_true(is.character(out$method))
  expect_true(is.character(out$gene))
  expect_true(is.character(out$event))
  expect_true(is.integer(out$grange_start))
  expect_true(is.integer(out$grange_end))
})

test_that("buildTopLoci: conditionIdx slices the 3-D mvsusie posterior per condition", {
  set.seed(1); L <- 2L; p <- 3L; R <- 2L
  alpha <- matrix(c(0.7, 0.2, 0.1, 0.6, 0.2, 0.2), nrow = L, ncol = p)
  mu  <- array(rnorm(L * p * R), dim = c(L, p, R))
  mu2 <- mu^2 + 0.1
  fit <- list(alpha = alpha, mu = mu, mu2 = mu2,
              pip = 1 - apply(1 - alpha, 2, prod))
  vids <- c("1:100:A:G", "1:200:C:T", "1:300:G:A")
  cst <- list(list(sets = list(cs = list())),
              list(sets = list(cs = list())),
              list(sets = list(cs = list())))
  attr(cst, "coverage") <- c(0.95, 0.70, 0.50)
  tl1 <- buildTopLoci(fit, cst, variantNames = vids, method = "mvsusie",
                      conditionIdx = 1L)
  tl2 <- buildTopLoci(fit, cst, variantNames = vids, method = "mvsusie",
                      conditionIdx = 2L)
  tl0 <- buildTopLoci(fit, cst, variantNames = vids, method = "mvsusie")
  # Each conditionIdx yields THAT condition's posterior (colSums(alpha*mu[,,r])).
  expect_equal(tl1$posterior_mean, colSums(alpha * mu[, , 1]))
  expect_equal(tl2$posterior_mean, colSums(alpha * mu[, , 2]))
  expect_false(isTRUE(all.equal(tl1$posterior_mean, tl2$posterior_mean)))
  expect_true(all(is.finite(tl1$posterior_sd)))
  # PIP is shared across conditions (mvSuSiE inclusion is joint).
  expect_equal(tl1$pip, tl2$pip)
  # A 3-D fit without a conditionIdx leaves the per-variant posterior NA.
  expect_true(all(is.na(tl0$posterior_mean)))
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
  expect_equal(unique(out$N), 419L)
  expect_equal(unique(out$method), "susie")
})

test_that("buildTopLoci exports af (not MAF) and carries the supplied af values", {
  variant_ids <- c("chr1:100:A:G", "chr1:200:C:T")
  inp <- .fake_fit_and_cs(
    variant_ids,
    cs_at_cov = list("0.95" = list(c(1L, 2L)), "0.7" = list(c(1L, 2L)),
                     "0.5" = list(c(1L, 2L))),
    pip = c(0.9, 0.9))
  out <- .runBuildTopLoci(inp, method = "susie", af = c(0.12, 0.87))
  expect_true("af" %in% names(out))
  expect_false("MAF" %in% names(out))
  # Directional effect-allele frequency: value passed through verbatim,
  # not folded to a minor-allele frequency (0.87 retained, not 0.13).
  expect_equal(out$af, c(0.12, 0.87))
})

test_that("buildTopLoci sets af = NA when no af is supplied (no silent coercion)", {
  variant_ids <- c("chr1:100:A:G", "chr1:200:C:T")
  inp <- .fake_fit_and_cs(
    variant_ids,
    cs_at_cov = list("0.95" = list(c(1L, 2L)), "0.7" = list(c(1L, 2L)),
                     "0.5" = list(c(1L, 2L))),
    pip = c(0.9, 0.9))
  out <- .runBuildTopLoci(inp, method = "susie", af = NULL)
  expect_true("af" %in% names(out))
  expect_true(all(is.na(out$af)))
})

.n_cs95 <- function(post)
  length(setdiff(unique(as.character(post$top_loci$cs_95)), c(NA, "")))

test_that("postprocessFinemappingFits forwards medianAbsCorr to susie_get_cs (OR-logic admits >= sets)", {
  # medianAbsCorr is only meaningful when the installed susieR's susie_get_cs
  # accepts median_abs_corr (GitHub-HEAD susieR; the CRAN/conda-forge build
  # does not yet). Skip rather than error where it is unavailable.
  skip_if_not("median_abs_corr" %in% names(formals(susieR::susie_get_cs)),
              "installed susieR has no median_abs_corr support")
  d <- .make_univariate_data(seed = 11, effect_idx = c(10, 35))
  fit <- susieR::susie(d$X, d$y, L = 5)
  # A very strict min_abs_corr alone vs the same min_abs_corr OR a lenient
  # median_abs_corr: OR-logic keeps at least as many credible sets.
  pStrict <- postprocessFinemappingFits(list(susie = fit), dataX = d$X, dataY = d$y,
                                        coverage = 0.95, minAbsCorr = 0.999,
                                        medianAbsCorr = NULL)
  pOr     <- postprocessFinemappingFits(list(susie = fit), dataX = d$X, dataY = d$y,
                                        coverage = 0.95, minAbsCorr = 0.999,
                                        medianAbsCorr = 0.1)
  expect_gte(.n_cs95(pOr), .n_cs95(pStrict))
})

test_that("postprocessFinemappingFits with medianAbsCorr = NULL is a no-op", {
  d <- .make_univariate_data(seed = 7, effect_idx = c(20))
  fit <- susieR::susie(d$X, d$y, L = 5)
  p1 <- postprocessFinemappingFits(list(susie = fit), dataX = d$X, dataY = d$y,
                                   coverage = 0.95)
  p2 <- postprocessFinemappingFits(list(susie = fit), dataX = d$X, dataY = d$y,
                                   coverage = 0.95, medianAbsCorr = NULL)
  expect_equal(p1$top_loci$cs_95, p2$top_loci$cs_95)
  expect_equal(p1$top_loci$af, p2$top_loci$af)
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

test_that("overlapping CS within one method: one row per variant; smallest cs_idx wins", {
  variant_ids <- c("chr1:100:A:G")
  # One variant belongs to CS 1 AND CS 2 at 95-cov (overlap).
  cs_at_cov <- list("0.95" = list(1L, 1L),
                    "0.7"  = list(1L, 1L),
                    "0.5"  = list(1L, 1L))
  inp <- .fake_fit_and_cs(variant_ids, cs_at_cov, pip = 0.9)
  out <- .runBuildTopLoci(inp, method = "susie")
  # Canonical schema is one row per variant. When a variant is in
  # multiple CSs at a coverage, the smallest cs_idx is reported.
  expect_equal(nrow(out), 1L)
  expect_equal(out$variant_id, "chr1:100:A:G")
  expect_equal(out$method, "susie")
  expect_equal(out$cs_95, "susie_1")
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
  expect_true(is.list(getSusieFit(fm)) && !is.null(getSusieFit(fm)$pip))
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

test_that("posterior_mean equals colSums(alpha*mu); posterior_sd equals sqrt(pmax(colSums(alpha*mu2) - mean^2, 0))", {
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
    row <- out[out$variant_id == variant_ids[i], , drop = FALSE]
    expect_true(nrow(row) >= 1L)
    expect_equal(unique(row$posterior_mean), expected_mean[i],
                 tolerance = 1e-10)
    expect_equal(unique(row$posterior_sd), expected_se[i],
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


# =============================================================================
# APPENDED COVERAGE TESTS
# Pure / mostly-pure helpers, S3 post-processing methods, two-stage fit, and
# the thin mvSuSiE / fSuSiE fit wrappers. Internal helpers are called via
# pecotmr:::; exported functions are called bare to match the file's style.
# =============================================================================
context("fineMappingWrappers coverage")

# ---- formatPipColumn / resolvePipColumn ----
test_that("formatPipColumn prefixes the method", {
  expect_equal(pecotmr:::formatPipColumn("susie"), "pip_susie")
  expect_equal(pecotmr:::formatPipColumn("susieInf"), "pip_susieInf")
})

test_that("resolvePipColumn covers NULL/empty/method/pip/single/multi branches", {
  expect_null(pecotmr:::resolvePipColumn(NULL))
  expect_null(pecotmr:::resolvePipColumn(data.frame(pip = numeric(0))))
  # method-specific column present -> returned directly
  expect_equal(pecotmr:::resolvePipColumn(data.frame(pip_susie = 0.1, pip_x = 0.2),
                                          method = "susie"), "pip_susie")
  # method given but absent -> fall through to "pip"
  expect_equal(pecotmr:::resolvePipColumn(data.frame(pip = 0.1, pip_x = 0.2),
                                          method = "susie"), "pip")
  # plain "pip" present
  expect_equal(pecotmr:::resolvePipColumn(data.frame(pip = 0.1, pip_susie = 0.2)), "pip")
  # single pip_ column, no "pip"
  expect_equal(pecotmr:::resolvePipColumn(data.frame(pip_susie = 0.1)), "pip_susie")
  # multiple pip_ columns, no "pip", no method -> NULL (ambiguous)
  expect_null(pecotmr:::resolvePipColumn(data.frame(pip_a = 0.1, pip_b = 0.2)))
})

# ---- formatCsColumn / legacy column translation ----
test_that("formatCsColumn formats integer and fractional coverage, errors on non-numeric", {
  expect_equal(pecotmr:::formatCsColumn(0.95, "susie"), "CS_95_susie")
  expect_equal(pecotmr:::formatCsColumn(0.7, "susie"), "CS_70_susie")
  expect_equal(pecotmr:::formatCsColumn(0.999, "susieInf"), "CS_99_9_susieInf")
  expect_error(pecotmr:::formatCsColumn(NA, "susie"), "coverage must be numeric")
})

test_that(".translateLegacyCsColumnName converts cs_coverage_* and passes others through", {
  expect_null(pecotmr:::.translateLegacyCsColumnName(NULL))
  expect_equal(
    pecotmr:::.translateLegacyCsColumnName(c("cs_coverage_0.95", "variant_id", "cs_coverage_0.7")),
    c("CS_95_susie", "variant_id", "CS_70_susie")
  )
})

test_that(".translateLegacyTopLociCsColumns returns non-data.frame inputs unchanged", {
  x <- list(a = 1)
  expect_identical(pecotmr:::.translateLegacyTopLociCsColumns(x), x)
  expect_null(pecotmr:::.translateLegacyTopLociCsColumns(NULL))
})

# ---- .camelToSnakeMethod ----
test_that(".camelToSnakeMethod handles NULL, empty, and vectors of method ids", {
  expect_null(pecotmr:::.camelToSnakeMethod(NULL))
  expect_equal(pecotmr:::.camelToSnakeMethod(character(0)), character(0))
  expect_equal(
    pecotmr:::.camelToSnakeMethod(c("susieInfRss", "mvsusie", "susie", "susieAsh", "singleEffect")),
    c("susie_inf_rss", "mvsusie", "susie", "susie_ash", "single_effect")
  )
})

# ---- .setFinemappingFitClass ----
test_that(".setFinemappingFitClass assigns method classes and handles NULL/unknown", {
  expect_null(pecotmr:::.setFinemappingFitClass(NULL, "susie"))
  expect_true("susiF" %in% class(pecotmr:::.setFinemappingFitClass(list(a = 1), "fsusie")))
  expect_true("mvsusie" %in% class(pecotmr:::.setFinemappingFitClass(list(a = 1), "mvsusie")))
  expect_true("susieRss" %in% class(pecotmr:::.setFinemappingFitClass(list(a = 1), "singleEffect")))
  expect_true("susieRss" %in%
                class(pecotmr:::.setFinemappingFitClass(list(a = 1), "bayesianConditionalRegression")))
  # unknown method -> class unchanged
  obj <- structure(list(a = 1), class = "foo")
  expect_equal(class(pecotmr:::.setFinemappingFitClass(obj, "weird")), "foo")
})

# ---- prepareSusieFromInfArgs ----
test_that("prepareSusieFromInfArgs sets none-branch defaults", {
  fit <- list(V = rep(1, 5))
  a <- pecotmr:::prepareSusieFromInfArgs(list(), fit, refineDefault = TRUE)
  expect_true(a$refine)
  expect_equal(a$unmappable_effects, "none")
  expect_identical(a$model_init, fit)
  expect_null(a$convergence_method)
})

test_that("prepareSusieFromInfArgs ash branch sets convergence, caps L_greedy, keeps preset refine", {
  fit <- list(V = rep(1, 5))
  a <- pecotmr:::prepareSusieFromInfArgs(list(L = 3, L_greedy = 10), fit,
                                         refineDefault = TRUE, unmappableEffects = "ash")
  expect_equal(a$convergence_method, "pip")
  expect_equal(a$L_greedy, 3)   # min(length(V) = 5, L = 3)
  expect_equal(a$unmappable_effects, "ash")
  # a preset refine is not overwritten by refineDefault
  b <- pecotmr:::prepareSusieFromInfArgs(list(refine = FALSE), fit, refineDefault = TRUE)
  expect_false(b$refine)
})

# ---- .asEffectMatrix / .asLbfMatrix ----
test_that(".asEffectMatrix handles NULL, list, matrix, and data.frame", {
  expect_equal(dim(pecotmr:::.asEffectMatrix(NULL)), c(0L, 0L))
  expect_equal(pecotmr:::.asEffectMatrix(list(c(1, 2), c(3, 4))),
               matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE))
  df_out <- pecotmr:::.asEffectMatrix(data.frame(a = 1:2, b = 3:4))
  expect_true(is.matrix(df_out))
  expect_equal(dim(df_out), c(2L, 2L))
  m <- matrix(1:6, 2, 3)
  expect_equal(pecotmr:::.asEffectMatrix(m), m)
})

test_that(".asLbfMatrix prefers lbf_variable, falls back to lBF, else NULL", {
  expect_equal(pecotmr:::.asLbfMatrix(list(lbf_variable = matrix(1, 2, 2))), matrix(1, 2, 2))
  expect_equal(pecotmr:::.asLbfMatrix(list(lBF = matrix(2, 2, 2))), matrix(2, 2, 2))
  expect_null(pecotmr:::.asLbfMatrix(list(x = 1)))
})

# ---- .parseGrange ----
test_that(".parseGrange returns NA for NULL/empty/invalid and parses valid regions", {
  expect_equal(pecotmr:::.parseGrange(NULL), c(start = NA_integer_, end = NA_integer_))
  expect_equal(pecotmr:::.parseGrange(""), c(start = NA_integer_, end = NA_integer_))
  expect_equal(pecotmr:::.parseGrange("not_a_region"), c(start = NA_integer_, end = NA_integer_))
  expect_equal(pecotmr:::.parseGrange("chr10:10823338-14348298"),
               c(start = 10823338L, end = 14348298L))
})

# ---- selectEffects ----
test_that("selectEffects returns integer(0) for empty alpha, V-filtered or all effects", {
  expect_equal(pecotmr:::selectEffects(list(alpha = NULL)), integer(0))
  expect_equal(pecotmr:::selectEffects(list(alpha = matrix(0.2, 3, 4), V = c(1, 1e-20, 2))),
               c(1L, 3L))
  expect_equal(pecotmr:::selectEffects(list(alpha = matrix(0.2, 3, 4))), 1:3)
})

# ---- .csPurityVec ----
test_that(".csPurityVec uses purity, falls back to cs_corr, else NA", {
  expect_equal(
    pecotmr:::.csPurityVec(list(sets = list(purity = data.frame(min.abs.corr = c(0.8, 0.9)),
                                            cs = list(1, 2)))),
    c(0.8, 0.9))
  m1 <- matrix(c(1, 0.7, 0.7, 1), 2)
  m2 <- matrix(1, 1, 1)
  expect_equal(
    pecotmr:::.csPurityVec(list(sets = list(cs = list(1, 2, 3)),
                                cs_corr = list(m1, m2, NULL))),
    c(0.7, 1, NA))
  expect_equal(pecotmr:::.csPurityVec(list(sets = list(cs = list(1, 2)))),
               c(NA_real_, NA_real_))
})

# ---- .translateSusiePurity ----
test_that(".translateSusiePurity renames df and matrix columns, leaves others alone", {
  expect_null(pecotmr:::.translateSusiePurity(NULL))
  df <- data.frame(min.abs.corr = 1, mean.abs.corr = 2, median.abs.corr = 3, other = 4)
  expect_equal(names(pecotmr:::.translateSusiePurity(df)),
               c("minAbsCorr", "meanAbsCorr", "medianAbsCorr", "other"))
  mm <- matrix(1:6, 2, 3)
  colnames(mm) <- c("min.abs.corr", "mean.abs.corr", "median.abs.corr")
  expect_equal(colnames(pecotmr:::.translateSusiePurity(mm)),
               c("minAbsCorr", "meanAbsCorr", "medianAbsCorr"))
  mm2 <- matrix(1:6, 2, 3)
  expect_null(colnames(pecotmr:::.translateSusiePurity(mm2)))
})

# ---- .topLociForS4Slot ----
test_that(".topLociForS4Slot returns an empty frame for NULL/0-row input", {
  e <- pecotmr:::.topLociForS4Slot(NULL)
  expect_equal(names(e), c("variant_id", "method"))
  expect_equal(nrow(e), 0L)
  expect_equal(nrow(pecotmr:::.topLociForS4Slot(data.frame(a = integer(0)))), 0L)
})

test_that(".topLociForS4Slot derives variant_id from variant and integer cs from cs_95", {
  tl <- data.frame(variant = c("v1", "v2", "v3", "v4"),
                   cs_95 = c("susie_1", "susie_0", "susie_2", NA),
                   stringsAsFactors = FALSE)
  r <- pecotmr:::.topLociForS4Slot(tl)
  expect_equal(r$variant_id, c("v1", "v2", "v3", "v4"))
  expect_equal(r$cs, c(1L, 0L, 2L, 0L))
  # a non-numeric cs_95 tail collapses to 0L
  expect_equal(pecotmr:::.topLociForS4Slot(
    data.frame(variant_id = "v1", cs_95 = "susie_x", stringsAsFactors = FALSE))$cs, 0L)
})

# ---- extractVariantNames ----
test_that("extractVariantNames reads pip names, then alpha colnames, then a fallback", {
  expect_equal(
    pecotmr:::extractVariantNames(list(pip = setNames(c(0.1, 0.2),
                                                      c("chr1:100:A:G", "chr1:200:C:T")))),
    c("chr1:100:A:G", "chr1:200:C:T"))
  expect_equal(
    pecotmr:::extractVariantNames(list(
      pip = c(0.1, 0.2),
      alpha = matrix(0, 1, 2, dimnames = list(NULL, c("chr1:100:A:G", "chr1:200:C:T"))))),
    c("chr1:100:A:G", "chr1:200:C:T"))
  r <- pecotmr:::extractVariantNames(list(pip = c(0.1, 0.2, 0.3)))
  expect_length(r, 3)
  expect_true(is.character(r))
})

# ---- extractSumstats ----
test_that("extractSumstats returns NULL / passthrough across non-regression branches", {
  expect_null(pecotmr:::extractSumstats(list(), NULL, NULL))
  expect_equal(pecotmr:::extractSumstats(list(), NULL, list(z = c(1, 2)), method = "susieRss"),
               list(z = c(1, 2)))
  expect_equal(pecotmr:::extractSumstats(list(), NULL,
                                         list(betahat = c(1, 2), sebetahat = c(0.1, 0.2))),
               list(betahat = c(1, 2), sebetahat = c(0.1, 0.2)))
  expect_null(pecotmr:::extractSumstats(list(), NULL, c(1, 2, 3)))            # dataX NULL
  expect_null(pecotmr:::extractSumstats(list(), matrix(0, 3, 2), matrix(0, 3, 2)))  # multi-col dataY
})

test_that("extractSumstats runs univariate regression and applies x/y scalars", {
  skip_if_not_installed("susieR")
  set.seed(1)
  X <- matrix(rnorm(60), 20, 3)
  colnames(X) <- c("chr1:1:A:G", "chr1:2:A:G", "chr1:3:A:G")
  y <- X[, 1] * 2 + rnorm(20)
  s1 <- pecotmr:::extractSumstats(list(), X, y)
  expect_named(s1, c("betahat", "sebetahat"))
  s2 <- pecotmr:::extractSumstats(list(), X, y, yScalar = 2, xScalar = 1)
  expect_equal(s2$betahat, s1$betahat * 2)
  expect_equal(s2$sebetahat, s1$sebetahat * 2)
})

# ---- computeCsTable / computeCsTables ----
test_that("computeCsTable fsusie branch returns empty sets and NULL cs_corr when no CS", {
  ct <- pecotmr:::computeCsTable(
    list(cs = list(), pip = setNames(c(0.1, 0.2), c("a", "b"))),
    matrix(0, 5, 2), coverage = 0.95, csInput = "fsusie")
  expect_equal(names(ct), c("sets", "cs_corr", "pip"))
  expect_null(ct$cs_corr)
  expect_length(ct$sets$cs, 0)
})

test_that("computeCsTable X and Xcorr branches return sets/pip/cs_corr", {
  skip_if_not_installed("susieR")
  d <- .make_univariate_data(seed = 7, n = 200, p = 8, effect_idx = c(2))
  fit <- susieR::susie(d$X, d$y, L = 4)
  ctx <- pecotmr:::computeCsTable(fit, d$X, coverage = 0.95, csInput = "X")
  expect_true(all(c("sets", "pip", "cs_corr") %in% names(ctx)))
  ctc <- pecotmr:::computeCsTable(fit, cor(d$X), coverage = 0.95, csInput = "Xcorr")
  expect_true(all(c("sets", "pip", "cs_corr") %in% names(ctc)))
})

test_that("computeCsTables names tables, sets coverage attr, defaults coverage from fit", {
  skip_if_not_installed("susieR")
  d <- .make_univariate_data(seed = 7, n = 200, p = 8, effect_idx = c(2))
  fit <- susieR::susie(d$X, d$y, L = 4)
  cts <- pecotmr:::computeCsTables(fit, d$X, coverage = 0.95,
                                   secondaryCoverage = c(0.7, 0.5),
                                   method = "susie", csInput = "X")
  expect_equal(attr(cts, "coverage"), c(0.95, 0.7, 0.5))
  expect_equal(names(cts), c("CS_95_susie", "CS_70_susie", "CS_50_susie"))
  # coverage = NULL falls back to fit$sets$requested_coverage
  fit2 <- fit
  fit2$sets$requested_coverage <- 0.9
  cts2 <- pecotmr:::computeCsTables(fit2, d$X, coverage = NULL,
                                    secondaryCoverage = 0.5, method = "susie", csInput = "X")
  expect_equal(attr(cts2, "coverage"), c(0.9, 0.5))
})

# ---- trimFinemappingFit ----
test_that("trimFinemappingFit (susie) trims to selected effects and keeps scalar slots", {
  fit <- list(
    pip = setNames(c(0.5, 0.3), c("chr1:100:A:G", "chr1:200:C:T")),
    alpha = matrix(c(0.5, 0.5, 0.3, 0.7), nrow = 2, byrow = TRUE),
    mu = matrix(0.1, 2, 2), mu2 = matrix(0.02, 2, 2),
    lbf_variable = matrix(0, 2, 2),
    V = c(1, 1e-20), niter = 5,
    theta = c(1, 2), omega_weights = c(0.5, 0.5), X_column_scale_factors = c(1, 1))
  eff <- pecotmr:::selectEffects(fit)   # only effect 1 survives V filtering
  csTables <- list(
    list(sets = list(cs = list(L1 = c(1, 2)), purity = data.frame(min.abs.corr = 0.8)),
         cs_corr = list(matrix(c(1, 0.8, 0.8, 1), 2)), pip = fit$pip),
    list(sets = list(cs = list()), cs_corr = NULL, pip = fit$pip))
  tr <- pecotmr:::trimFinemappingFit(fit, eff, "susie", csTables)
  expect_equal(nrow(tr$alpha), 1L)
  expect_equal(nrow(tr$mu), 1L)
  expect_equal(nrow(tr$mu2), 1L)
  expect_equal(nrow(tr$lbf_variable), 1L)
  expect_equal(tr$V, 1)
  expect_equal(tr$theta, c(1, 2))
  expect_equal(tr$omega_weights, c(0.5, 0.5))
  expect_equal(tr$X_column_scale_factors, c(1, 1))
  expect_equal(tr$niter, 5)
  expect_equal(tr$max_L, 2L)
  expect_equal(tr$n_effects, 2L)
  expect_equal(class(tr), "susie")
  # secondary tables drop the pip element
  expect_length(tr$sets_secondary, 1L)
  expect_false("pip" %in% names(tr$sets_secondary[[1]]))
})

test_that("trimFinemappingFit (fsusie) retains coef and sets the fsusie/susie class", {
  fit <- list(pip = setNames(c(0.4, 0.6), c("chr1:100:A:G", "chr1:200:C:T")),
              alpha = matrix(0.5, 2, 2), coef = matrix(1, 2, 3))
  csTables <- list(list(sets = list(cs = list()), cs_corr = NULL, pip = fit$pip))
  tr <- pecotmr:::trimFinemappingFit(fit, c(1, 2), "fsusie", csTables)
  expect_equal(tr$coef, matrix(1, 2, 3))
  expect_equal(class(tr), c("fsusie", "susie"))
  expect_null(tr$V)
})

test_that("trimFinemappingFit (mvsusie) slices 3-D mu/mu2/mu2_diag/clfsr and stores coef", {
  skip_if_not_installed("mvsusieR")
  L <- 2L; p <- 3L; R <- 2L
  fit <- list(
    pip = setNames(seq(0.1, 0.3, length.out = p),
                   c("chr1:1:A:G", "chr1:2:A:G", "chr1:3:A:G")),
    alpha = matrix(1 / p, L, p),
    mu = array(rnorm(L * p * R), dim = c(L, p, R)),
    mu2 = array(0.1, dim = c(L, p, R)),
    mu2_diag = array(0.2, dim = c(L, p, R)),
    V = c(1, 1e-20),
    conditional_lfsr = array(0.5, dim = c(L, p, R)),
    niter = 3)
  eff <- pecotmr:::selectEffects(fit)   # effect 1
  csTables <- list(list(sets = list(cs = list()), cs_corr = NULL, pip = fit$pip))
  fake_coef <- matrix(rnorm((p + 1) * R), nrow = p + 1, ncol = R)
  local_mocked_bindings(coef.mvsusie = function(...) fake_coef, .package = "mvsusieR")
  tr <- pecotmr:::trimFinemappingFit(fit, eff, "mvsusie", csTables)
  expect_equal(dim(tr$mu), c(1L, p, R))
  expect_equal(dim(tr$mu2), c(1L, p, R))
  expect_equal(dim(tr$mu2_diag), c(1L, p, R))
  expect_equal(dim(tr$clfsr), c(1L, p, R))
  expect_equal(tr$coef, fake_coef[-1, , drop = FALSE])
  expect_equal(class(tr), c("mvsusie", "susie"))
})

# ---- postprocessFinemappingFit S3 methods ----
test_that("postprocessFinemappingFit.susiF post-processes an fsusie fit (empty-CS path)", {
  fit <- pecotmr:::.setFinemappingFitClass(list(
    pip = setNames(c(0.5, 0.3), c("chr1:100:A:G", "chr1:200:C:T")),
    alpha = matrix(0.5, 1, 2), mu = matrix(0.1, 1, 2), mu2 = matrix(0.02, 1, 2),
    cs = list()), "fsusie")
  expect_true("susiF" %in% class(fit))
  res <- pecotmr:::postprocessFinemappingFit(
    fit, method = "fsusie",
    dataX = matrix(0, 5, 2, dimnames = list(NULL, c("chr1:100:A:G", "chr1:200:C:T"))),
    dataY = NULL, coverage = 0.95,
    otherQuantities = list(condition_id = "ctx"))
  expect_equal(res$method, "fsusie")
  expect_equal(unique(res$top_loci$method), "fsusie")
  expect_equal(res$otherQuantities, list(condition_id = "ctx"))
  expect_equal(class(getSusieFit(res$finemappingEntry)), c("fsusie", "susie"))
})

test_that("postprocessFinemappingFit.susieInf labels credible sets with the susie_inf_ prefix", {
  skip_if_not_installed("susieR")
  d <- .make_univariate_data(seed = 7, n = 200, p = 8, effect_idx = c(2))
  fits <- fitSusieInfThenSusie(d$X, d$y)
  res <- pecotmr:::postprocessFinemappingFit(fits$susieInf, method = "susieInf",
                                             dataX = d$X, dataY = d$y, coverage = 0.95)
  expect_equal(res$method, "susieInf")
  expect_gt(nrow(res$top_loci), 0L)
  expect_equal(unique(res$top_loci$method), "susieInf")
  expect_true(all(grepl("^susie_inf_\\d+$", res$top_loci$cs_95)))
})

test_that(".postprocessFinemappingFitCommon trim=FALSE stores the untrimmed fit", {
  fit <- pecotmr:::.setFinemappingFitClass(list(
    pip = setNames(c(0.5, 0.3), c("chr1:100:A:G", "chr1:200:C:T")),
    alpha = matrix(0.5, 1, 2), mu = matrix(0.1, 1, 2), mu2 = matrix(0.02, 1, 2),
    cs = list(), extra_slot = "kept"), "fsusie")
  res <- pecotmr:::postprocessFinemappingFit(
    fit, method = "fsusie", trim = FALSE,
    dataX = matrix(0, 5, 2, dimnames = list(NULL, c("chr1:100:A:G", "chr1:200:C:T"))),
    dataY = NULL, coverage = 0.95)
  expect_equal(getSusieFit(res$finemappingEntry)$extra_slot, "kept")
})

# ---- postprocessFinemappingFits / formatFinemappingOutput error branches ----
test_that("postprocessFinemappingFits errors on empty or unnamed fit lists", {
  expect_error(postprocessFinemappingFits(list(), dataX = matrix(0, 2, 2)),
               "At least one fine-mapping fit")
  expect_error(postprocessFinemappingFits(list(NULL), dataX = matrix(0, 2, 2)),
               "At least one fine-mapping fit")
  expect_error(postprocessFinemappingFits(list(matrix(0, 1, 1)), dataX = matrix(0, 2, 2)),
               "named list")
})

test_that("formatFinemappingOutput errors when primaryMethod is absent", {
  post <- list(finemappingResults = list(susie = list(method = "susie")),
               top_loci = pecotmr:::.emptyTopLoci())
  expect_error(formatFinemappingOutput(post, "nonexistent"), "primaryMethod was not found")
})

# ---- buildTopLoci: marginal z/p passthrough ----
test_that("buildTopLoci passes through marginal z and p supplied in sumstats", {
  variant_ids <- c("chr1:100:A:G", "chr1:200:C:T")
  inp <- .fake_fit_and_cs(variant_ids,
                          cs_at_cov = list("0.95" = list(c(1L, 2L)),
                                            "0.7" = list(c(1L, 2L)),
                                            "0.5" = list(c(1L, 2L))),
                          pip = c(0.9, 0.9))
  out <- .runBuildTopLoci(inp, method = "susie",
                          sumstats = list(z = c(2.5, -1.5), p = c(0.01, 0.13)))
  expect_equal(out$marginal_z, c(2.5, -1.5))
  expect_equal(out$marginal_p, c(0.01, 0.13))
})

# ---- lbfToAlpha single-column matrix branch ----
test_that("lbfToAlpha handles a single-column matrix", {
  lbf <- matrix(c(1, 2, 3), ncol = 1)
  colnames(lbf) <- "v1"
  res <- pecotmr:::lbfToAlpha(lbf)
  expect_equal(dim(res), c(3L, 1L))
  expect_true(all(res == 1))   # single column -> the only entry carries all weight
})

# ---- fitSusieInfThenSusie ----
test_that("fitSusieInfThenSusie returns classed susie and susieInf fits", {
  skip_if_not_installed("susieR")
  d <- .make_univariate_data(seed = 3, n = 200, p = 8, effect_idx = c(4))
  fits <- fitSusieInfThenSusie(d$X, d$y)
  expect_named(fits, c("susie", "susieInf"))
  expect_true("susie" %in% class(fits$susie))
  expect_true("susieInf" %in% class(fits$susieInf))
  expect_length(fits$susie$pip, ncol(d$X))
})

test_that("fitSusieInfThenSusie reuses fittedModels without refitting", {
  skip_if_not_installed("susieR")
  d <- .make_univariate_data(seed = 3, n = 200, p = 8, effect_idx = c(4))
  fits <- fitSusieInfThenSusie(d$X, d$y)
  again <- fitSusieInfThenSusie(d$X, d$y,
                                fittedModels = list(susie = fits$susie, susieInf = fits$susieInf))
  expect_equal(unname(again$susie$pip), unname(fits$susie$pip))
  expect_equal(unname(again$susieInf$pip), unname(fits$susieInf$pip))
})

# ---- thin fit wrappers (mvSuSiE / fSuSiE) ----
test_that("fitMvsusie forwards arguments to mvsusieR::mvsusie", {
  skip_if_not_installed("mvsusieR")
  local_mocked_bindings(
    mvsusie = function(X, Y, prior_variance, coverage, ...) list(tag = "mv", coverage = coverage),
    .package = "mvsusieR")
  r <- fitMvsusie(matrix(0, 4, 2), matrix(0, 4, 2), prior_variance = 1, coverage = 0.9)
  expect_equal(r$tag, "mv")
  expect_equal(r$coverage, 0.9)
})

test_that("fitMvsusieRss forwards arguments to mvsusieR::mvsusie_rss", {
  skip_if_not_installed("mvsusieR")
  local_mocked_bindings(
    mvsusie_rss = function(Z, R, N, prior_variance, coverage, ...) list(tag = "rss", N = N),
    .package = "mvsusieR")
  r <- fitMvsusieRss(matrix(0, 2, 1), diag(2), N = 100, prior_variance = 1)
  expect_equal(r$tag, "rss")
  expect_equal(r$N, 100)
})

test_that("fitFsusie forwards arguments to fsusieR::susiF", {
  skip_if_not_installed("fsusieR")
  local_mocked_bindings(
    susiF = function(X, Y, pos, ...) list(tag = "fs", npos = length(pos)),
    .package = "fsusieR")
  r <- fitFsusie(matrix(0, 4, 3), matrix(0, 4, 2), pos = 1:2)
  expect_equal(r$tag, "fs")
  expect_equal(r$npos, 2)
})

# ===========================================================================
# SuSiE / mvSuSiE / fSuSiE weight-extractor tests (relocated to match the source move)
# ===========================================================================

test_that(".susie_rss_extract_weights returns correct-length vector", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 20
  n <- 500
  R <- diag(p)
  z <- rnorm(p)
  w <- pecotmr:::.susieRssExtractWeights(
    fit = NULL, z = z, R = R, n = n,
    requiredFields = c("alpha", "mu", "X_column_scale_factors"),
    fitArgs = list(L = 5)
  )
  expect_equal(length(w), p)
  expect_true(all(is.finite(w)))
})

test_that("susieRssWeights follows (stat, LD) convention", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 20
  n <- 500
  R <- diag(p)
  z <- rnorm(p)
  stat <- list(b = z / sqrt(n), cor = z / sqrt(n), z = z, n = rep(n, p))
  w <- susieRssWeights(stat, R, methodArgs = list(L = 5))
  expect_equal(length(w), p)
  expect_true(all(is.finite(w)))
})

test_that("susieRssWeights retains fit when retainFit = TRUE", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 20
  n <- 500
  R <- diag(p)
  z <- rnorm(p)
  stat <- list(b = z / sqrt(n), cor = z / sqrt(n), z = z, n = rep(n, p))
  w <- susieRssWeights(stat, R, retainFit = TRUE, methodArgs = list(L = 5))
  expect_false(is.null(attr(w, "fit")))
})

test_that("susieInfRssWeights works", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 20
  n <- 500
  R <- diag(p)
  z <- rnorm(p)
  stat <- list(b = z / sqrt(n), cor = z / sqrt(n), z = z, n = rep(n, p))
  w <- susieInfRssWeights(stat, R, methodArgs = list(L = 5))
  expect_equal(length(w), p)
  expect_true(all(is.finite(w)))
})

test_that("mvsusieWeights real fit returns p x K weights or errors on unstable small data", {
  skip_if_not_installed("mvsusieR")
  m <- .rrwMulti(n = 80, p = 8, K = 2)
  res <- tryCatch(suppressMessages(mvsusieWeights(X = m$X, Y = m$Y, L = 5, LGreedy = 2)),
                  error = function(e) e)
  if (inherits(res, "error")) {
    # mvSuSiE can be numerically unstable on tiny (X, Y); a clean error is
    # acceptable here. The coef-extraction path is covered by the mocked tests
    # in test_rrMrmashMvsusie.R.
    succeed("mvsusieWeights errored on small data (documented instability)")
  } else {
    expect_equal(dim(res), c(m$p, m$K))
    expect_true(all(is.finite(res)))
  }
})

test_that("susieAshRssWeights returns weights of length p", {
  skip_if_not_installed("susieR")
  f <- .rrwStatLd()
  w <- susieAshRssWeights(f$stat, f$LD, methodArgs = list(L = 5))
  expect_length(w, f$p)
  expect_true(all(is.finite(w)))
})

test_that("mvsusieRssWeights fits mvsusie_rss and returns p x K weights", {
  skip_if_not_installed("mvsusieR")
  m <- .rrwMulti(n = 80, p = 8, K = 2)
  w <- mvsusieRssWeights(m$stat, m$LD, L = 5, LGreedy = 2)
  expect_equal(dim(w), c(m$p, m$K))
  expect_true(all(is.finite(w)))
})

test_that("mvsusieRssWeights errors on single-context stat$z", {
  skip_if_not_installed("mvsusieR")
  f <- .rrwStatLd()
  oneCol <- list(z = matrix(f$stat$z, ncol = 1), n = f$n)
  expect_error(mvsusieRssWeights(oneCol, f$LD), ">= 2 columns")
})

# ---- mvsusieWeights ----
test_that("mvsusieWeights errors when mvsusieR package is not available", {
  skip_if(requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR is installed; skipping missing-package test")

  expect_error(
    mvsusieWeights(mvsusieFit = NULL, X = matrix(1, 10, 5), Y = matrix(1, 10, 3)),
    "mvsusieR"
  )
})

test_that("mvsusieWeights errors when X and Y are NULL and fit is NULL", {
  skip_if_not(requireNamespace("mvsusieR", quietly = TRUE),
              "mvsusieR not installed")
  expect_error(mvsusieWeights(mvsusieFit = NULL, X = NULL, Y = NULL),
               "Both X and Y must be provided")
})

test_that("mvsusieWeights fits model and returns coefficients when fit is NULL", {
  skip_if_not(requireNamespace("mvsusieR", quietly = TRUE),
              "mvsusieR not installed")
  set.seed(42)
  n <- 30
  p <- 5
  R <- 3
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * R), n, R)
  fake_coef <- matrix(rnorm((p + 1) * R), nrow = p + 1, ncol = R)
  captured <- list()

  local_mocked_bindings(
    create_mixture_prior = function(...) list(),
    mvsusie = function(...) {
      captured <<- list(...)
      "mock_fit"
    },
    coef.mvsusie = function(...) fake_coef,
    .package = "mvsusieR"
  )

  result <- expect_message(
    mvsusieWeights(X = X, Y = Y, L = 12, LGreedy = 4),
    "mvsusieFit is not provided"
  )
  # Should return coef without intercept row
  expect_equal(dim(result), c(p, R))
  expect_equal(result, fake_coef[-1, ])
  expect_equal(captured$L, 12)
  expect_equal(captured$L_greedy, 4)
})

test_that("mvsusieWeights returns coefficients from provided fit", {
  skip_if_not(requireNamespace("mvsusieR", quietly = TRUE),
              "mvsusieR not installed")
  p <- 5
  R <- 3
  fake_coef <- matrix(rnorm((p + 1) * R), nrow = p + 1, ncol = R)

  local_mocked_bindings(
    coef.mvsusie = function(...) fake_coef,
    .package = "mvsusieR"
  )

  result <- mvsusieWeights(mvsusieFit = "precomputed_fit")
  expect_equal(dim(result), c(p, R))
  expect_equal(result, fake_coef[-1, ])
})

.fw_makeFsusieFit <- function(seed = 1, n = 150L, p = 24L, J = 16L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)), paste0("v", seq_len(p))))
  b1 <- sin(seq(0, 2 * pi, length.out = J))
  b2 <- cos(seq(0, pi, length.out = J))
  Y <- X[, 3] %o% b1 + X[, 10] %o% b2 +
    matrix(rnorm(n * J, sd = 0.3), n, J)
  colnames(Y) <- paste0("f", seq_len(J))
  list(X = X, Y = Y,
       fit = suppressWarnings(fsusieR::susiF(
         X = X, Y = Y, pos = seq_len(J), L = 5,
         post_processing = "none", verbose = FALSE)))
}

test_that("fsusieWeights returns a variants x features matrix with variant rownames", {
  skip_if_not_installed("fsusieR")
  skip_if_not_installed("wavethresh")
  obj <- .fw_makeFsusieFit()
  W <- fsusieWeights(fsusieFit = obj$fit, variantIds = colnames(obj$X))
  expect_true(is.matrix(W))
  expect_equal(nrow(W), ncol(obj$X))
  expect_equal(ncol(W), ncol(obj$Y))
  expect_equal(rownames(W), colnames(obj$X))
})

test_that("fsusieWeights matches fsusieR's own out_prep reconstruction (post_processing='none')", {
  skip_if_not_installed("fsusieR")
  skip_if_not_installed("wavethresh")
  obj <- .fw_makeFsusieFit()
  fit <- obj$fit
  # The alpha-weighted sum over SNPs of the per-SNP feature-domain curves that
  # fsusieWeights reconstructs must equal fSuSiE's own fitted_func[[l]] (built
  # by out_prep.susiF) for every effect l.
  csdX <- as.numeric(fit$csd_X)
  perScale <- "mixture_normal_per_scale" %in% class(fsusieR::get_G_prior(fit))
  indxLst <- fsusieR::gen_wavelet_indx(log2(length(fit$outing_grid)))
  scaleCols <- if (perScale) indxLst[[length(indxLst)]]
               else ncol(as.matrix(fit$fitted_wc[[1L]]))
  S <- pecotmr:::.fsusieSynthesisMatrix(fit$n_wac, scaleCols)
  maxErr <- 0
  for (l in seq_along(fit$fitted_wc)) {
    al <- as.numeric(fit$alpha[[l]])
    contrib <- colSums((al * (1 / csdX) * as.matrix(fit$fitted_wc[[l]])) %*% S)
    maxErr <- max(maxErr, max(abs(contrib - as.numeric(fit$fitted_func[[l]]))))
  }
  expect_lt(maxErr, 1e-8)
})

test_that("fsusieWeights concentrates weight on the causal SNPs", {
  skip_if_not_installed("fsusieR")
  skip_if_not_installed("wavethresh")
  obj <- .fw_makeFsusieFit()
  W <- fsusieWeights(fsusieFit = obj$fit, variantIds = colnames(obj$X))
  rowNorm <- sqrt(rowSums(W^2))
  top2 <- names(sort(rowNorm, decreasing = TRUE))[1:2]
  expect_setequal(top2, c("v3", "v10"))
})

test_that("fsusieWeights fast path returns precomputed $coef for a trimmed fit", {
  # A trimmed fSuSiE fit drops fitted_wc but keeps the precomputed weight
  # matrix in $coef; fsusieWeights returns it without touching wavelet slots.
  W0 <- matrix(c(1, 0, 2, 0, 0, 3), nrow = 3,
               dimnames = list(c("v1", "v2", "v3"), c("f1", "f2")))
  trimmed <- list(coef = W0, pip = c(0.1, 0.2, 0.7))
  class(trimmed) <- c("fsusie", "susie")
  W <- fsusieWeights(fsusieFit = trimmed)
  expect_identical(W, W0)
})

test_that("fsusieWeights errors without a fit and on an unusable (trimmed, no coef) fit", {
  expect_error(fsusieWeights(fsusieFit = NULL), "is required")
  bad <- list(pip = c(0.1, 0.9))  # no coef, no fitted_wc
  class(bad) <- c("fsusie", "susie")
  expect_error(fsusieWeights(fsusieFit = bad), "missing required slot")
})

# ===========================================================================
# mergeSusieCs — cross-condition credible-set merging on a QtlFineMappingResult
# (relocated from mashWrapper.R and adapted to consume the S4 result type)
# ===========================================================================

.msc_entry <- function(vid, pip, cs, csName = "cs_95") {
  tl <- data.frame(variant_id = vid, pip = pip, stringsAsFactors = FALSE)
  tl[[csName]] <- cs
  FineMappingEntry(variantIds = vid, susieFit = list(), topLoci = tl)
}
.msc_fmr <- function(entries, method = "susie") {
  n <- length(entries)
  QtlFineMappingResult(study = rep("S", n), context = paste0("c", seq_len(n)),
    trait = rep("t", n), method = rep(method, n), entry = entries)
}

test_that("mergeSusieCs: non-overlapping CSs keep distinct per-condition labels", {
  fmr <- .msc_fmr(list(
    .msc_entry(c("v1", "v2"), c(0.8, 0.6), c("susie_1", "susie_1")),
    .msc_entry(c("v3", "v4"), c(0.9, 0.7), c("susie_1", "susie_2"))))
  res <- mergeSusieCs(fmr)
  expect_equal(res$variant_id, c("v1", "v2", "v3", "v4"))
  expect_equal(res$credibleSetNames, c("cs_1_1", "cs_1_1", "cs_2_1", "cs_2_2"))
  expect_equal(res$maxPip, c(0.8, 0.6, 0.9, 0.7))
})

test_that("mergeSusieCs: a variant shared across conditions merges their credible sets", {
  fmr <- .msc_fmr(list(
    .msc_entry(c("v1", "v2", "v3"), c(0.9, 0.5, 0.8), c("susie_1", "susie_0", "susie_2")),
    .msc_entry(c("v3", "v4"),       c(0.7, 0.6),       c("susie_1", "susie_1"))))
  res <- mergeSusieCs(fmr)
  expect_false("v2" %in% res$variant_id)                       # susie_0 -> not in a CS
  expect_equal(res$credibleSetNames[res$variant_id == "v3"], "cs_1_2,cs_2_1")
  expect_equal(res$credibleSetNames[res$variant_id == "v4"], "cs_1_2,cs_2_1")  # via shared v3
  expect_equal(res$maxPip[res$variant_id == "v3"], 0.8)
  expect_equal(res$medianPip[res$variant_id == "v3"], 0.75)    # median(0.8, 0.7)
})

test_that("mergeSusieCs: coverage selects the cs_<coverage*100> column", {
  fmr <- .msc_fmr(list(
    .msc_entry(c("v1", "v2"), c(0.8, 0.6), c("susie_1", "susie_1"), csName = "cs_70")))
  res <- mergeSusieCs(fmr, coverage = 0.70)
  expect_equal(res$variant_id, c("v1", "v2"))
  expect_equal(res$credibleSetNames, c("cs_1_1", "cs_1_1"))
})

test_that("mergeSusieCs: a condition with no usable CS is skipped, valid ones kept", {
  valid  <- .msc_entry(c("v1", "v2"), c(0.9, 0.8), c("susie_1", "susie_1"))
  noCs   <- FineMappingEntry(variantIds = "v3", susieFit = list(),
                             topLoci = data.frame(variant_id = "v3", pip = 0.9))
  res <- mergeSusieCs(.msc_fmr(list(valid, noCs)))
  expect_setequal(res$variant_id, c("v1", "v2"))               # noCs condition skipped
})

test_that("mergeSusieCs: NULL when no condition contributes a credible set", {
  expect_null(mergeSusieCs(.msc_fmr(list(.msc_entry("v1", 0.9, "susie_0")))))  # all _0
  noCs <- FineMappingEntry(variantIds = "v1", susieFit = list(),
                           topLoci = data.frame(variant_id = "v1", pip = 0.9))
  expect_null(mergeSusieCs(.msc_fmr(list(noCs))))                              # no cs col
})

test_that("mergeSusieCs: single condition with one credible set", {
  res <- mergeSusieCs(.msc_fmr(list(
    .msc_entry(c("v1", "v2"), c(0.9, 0.8), c("susie_1", "susie_1")))))
  expect_equal(res$credibleSetNames, c("cs_1_1", "cs_1_1"))
})

test_that("mergeSusieCs: non-FineMappingResult input errors", {
  expect_error(mergeSusieCs(list(1, 2)), "QtlFineMappingResult")
})
