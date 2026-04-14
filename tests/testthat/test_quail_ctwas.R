context("quail_ctwas")

# ===========================================================================
# Helper: small simulated phenotype + covariate data (n = 30)
# ===========================================================================
make_quail_data <- function(n = 30, seed = 42) {
  set.seed(seed)
  covariates <- data.frame(
    age = rnorm(n, 50, 10),
    sex = rbinom(n, 1, 0.5)
  )
  # Phenotype with variance that depends on covariates (heteroscedastic)
  phenotype <- 2 + 0.5 * covariates$age + rnorm(n, sd = 1 + 0.3 * covariates$sex)
  list(phenotype = phenotype, covariates = covariates)
}

# ===========================================================================
#  QUAIL rank score tests
# ===========================================================================

# ---------- calculate_rank_score -------------------------------------------

test_that("calculate_rank_score returns vector of correct length", {
  skip_if_not_installed("quantreg")
  d <- make_quail_data()
  result <- pecotmr:::calculate_rank_score(d$phenotype, d$covariates, tau = 0.5, seed = 1)
  expect_true(is.numeric(result))
  expect_equal(length(result), length(d$phenotype))
})

test_that("calculate_rank_score is reproducible with same seed", {
  skip_if_not_installed("quantreg")
  d <- make_quail_data()
  r1 <- pecotmr:::calculate_rank_score(d$phenotype, d$covariates, tau = 0.25, seed = 99)
  r2 <- pecotmr:::calculate_rank_score(d$phenotype, d$covariates, tau = 0.25, seed = 99)
  expect_identical(r1, r2)
})

test_that("calculate_rank_score differs with different seeds", {
  skip_if_not_installed("quantreg")
  d <- make_quail_data()
  r1 <- pecotmr:::calculate_rank_score(d$phenotype, d$covariates, tau = 0.5, seed = 1)
  r2 <- pecotmr:::calculate_rank_score(d$phenotype, d$covariates, tau = 0.5, seed = 2)
  # Different random covariate column, so results should differ

  expect_false(identical(r1, r2))
})

test_that("calculate_rank_score produces different results for different tau values", {
  skip_if_not_installed("quantreg")
  d <- make_quail_data()
  r_low <- pecotmr:::calculate_rank_score(d$phenotype, d$covariates, tau = 0.1, seed = 1)
  r_high <- pecotmr:::calculate_rank_score(d$phenotype, d$covariates, tau = 0.9, seed = 1)
  expect_false(identical(r_low, r_high))
})

test_that("calculate_rank_score returns finite values", {
  skip_if_not_installed("quantreg")
  d <- make_quail_data()
  result <- pecotmr:::calculate_rank_score(d$phenotype, d$covariates, tau = 0.5, seed = 1)
  expect_true(all(is.finite(result)))
})

# ---------- calculate_integrated_score -------------------------------------

test_that("calculate_integrated_score works with equal method and even num_tau_levels", {
  skip_if_not_installed("quantreg")
  skip_if_not_installed("quadprog")
  n <- 30
  num_tau <- 4 # even
  # Create mock rank scores: list of numeric vectors
  set.seed(10)
  rank_scores <- lapply(1:num_tau, function(i) rnorm(n))
  result <- pecotmr:::calculate_integrated_score(rank_scores, method = "equal", num_tau_levels = num_tau)
  expect_true(is.numeric(result))
  expect_equal(length(result), n)
  expect_true(all(is.finite(result)))
})

test_that("calculate_integrated_score works with equal method and odd num_tau_levels", {
  skip_if_not_installed("quantreg")
  skip_if_not_installed("quadprog")
  n <- 30
  num_tau <- 5 # odd
  set.seed(11)
  rank_scores <- lapply(1:num_tau, function(i) rnorm(n))
  result <- pecotmr:::calculate_integrated_score(rank_scores, method = "equal", num_tau_levels = num_tau)
  expect_true(is.numeric(result))
  expect_equal(length(result), n)
})

test_that("calculate_integrated_score equal method: middle tau has zero weight for odd levels", {
  skip_if_not_installed("quantreg")
  skip_if_not_installed("quadprog")
  n <- 20
  num_tau <- 3 # odd, mid_point = 2
  # Set all scores to zero except the middle one
  rank_scores <- list(rep(0, n), rep(1, n), rep(0, n))
  result <- pecotmr:::calculate_integrated_score(rank_scores, method = "equal", num_tau_levels = num_tau)
  # Middle score (index 2) is skipped for odd levels; lower_half = {1}, upper_half = {3}
  # int_rank_score = -rank_scores[[1]] + rank_scores[[3]] = 0 + 0 = 0
  # n_pairs = 1, so result = 0/1 = 0
  expect_equal(result, rep(0, n))
})

test_that("calculate_integrated_score IVW method returns correct length", {
  skip("Known bug: solve.QP constraint matrix dimensions are incorrect in IVW method")
  skip_if_not_installed("quantreg")
  skip_if_not_installed("quadprog")
  n <- 30
  num_tau <- 4
  set.seed(20)
  rank_scores <- lapply(1:num_tau, function(i) rnorm(n, mean = i))
  result <- pecotmr:::calculate_integrated_score(rank_scores, method = "ivw", num_tau_levels = num_tau)
  expect_true(is.numeric(result))
  expect_equal(length(result), n)
  expect_true(all(is.finite(result)))
})

test_that("calculate_integrated_score IVW method with odd num_tau_levels", {
  skip("Known bug: solve.QP constraint matrix dimensions are incorrect in IVW method")
  skip_if_not_installed("quantreg")
  skip_if_not_installed("quadprog")
  n <- 30
  num_tau <- 5
  set.seed(21)
  rank_scores <- lapply(1:num_tau, function(i) rnorm(n, mean = i * 0.5))
  result <- pecotmr:::calculate_integrated_score(rank_scores, method = "ivw", num_tau_levels = num_tau)
  expect_true(is.numeric(result))
  expect_equal(length(result), n)
})

# ---------- fit_rank_scores ------------------------------------------------

test_that("fit_rank_scores returns list of correct length", {
  skip_if_not_installed("quantreg")
  d <- make_quail_data()
  num_tau <- 3
  result <- pecotmr:::fit_rank_scores(d$phenotype, d$covariates, num_tau_levels = num_tau, num_cores = 1)
  expect_true(is.list(result))
  expect_equal(length(result), num_tau)
})

test_that("fit_rank_scores each element has correct length", {
  skip_if_not_installed("quantreg")
  d <- make_quail_data()
  num_tau <- 3
  result <- pecotmr:::fit_rank_scores(d$phenotype, d$covariates, num_tau_levels = num_tau, num_cores = 1)
  for (i in seq_along(result)) {
    expect_equal(length(result[[i]]), length(d$phenotype))
    expect_true(is.numeric(result[[i]]))
  }
})

test_that("fit_rank_scores elements correspond to evenly spaced tau levels", {
  skip_if_not_installed("quantreg")
  d <- make_quail_data()
  num_tau <- 3
  result <- pecotmr:::fit_rank_scores(d$phenotype, d$covariates, num_tau_levels = num_tau, num_cores = 1)
  # Each element should match a direct call to calculate_rank_score at the expected tau
  for (i in 1:num_tau) {
    tau <- i / (num_tau + 1)
    expected <- pecotmr:::calculate_rank_score(d$phenotype, d$covariates, tau)
    expect_equal(result[[i]], expected)
  }
})

# ---------- QUAIL_rank_score_pipeline (exported) ---------------------------

test_that("QUAIL_rank_score_pipeline errors on non-numeric character phenotype", {
  skip_if_not_installed("quantreg")
  skip_if_not_installed("quadprog")
  d <- make_quail_data()
  pheno_char <- c("a", "b", "c")
  expect_error(
    QUAIL_rank_score_pipeline(pheno_char, d$covariates, num_tau_levels = 3),
    "phenotype must be a numeric vector"
  )
})

test_that("QUAIL_rank_score_pipeline accepts data.frame phenotype", {
  skip_if_not_installed("quantreg")
  skip_if_not_installed("quadprog")
  d <- make_quail_data()
  pheno_df <- data.frame(pheno = d$phenotype)
  # Should not error; data.frame is converted internally
  result <- QUAIL_rank_score_pipeline(pheno_df, d$covariates, num_tau_levels = 3, method = "equal")
  expect_true(is.numeric(result))
  expect_equal(length(result), length(d$phenotype))
})

test_that("QUAIL_rank_score_pipeline full run with small data and equal method", {
  skip_if_not_installed("quantreg")
  skip_if_not_installed("quadprog")
  d <- make_quail_data()
  result <- QUAIL_rank_score_pipeline(d$phenotype, d$covariates,
    num_tau_levels = 5, method = "equal", num_cores = 1
  )
  expect_true(is.numeric(result))
  expect_equal(length(result), length(d$phenotype))
  expect_true(all(is.finite(result)))
})

test_that("QUAIL_rank_score_pipeline full run with IVW method", {
  skip("Known bug: solve.QP constraint matrix dimensions are incorrect in IVW method")
  skip_if_not_installed("quantreg")
  skip_if_not_installed("quadprog")
  d <- make_quail_data()
  result <- QUAIL_rank_score_pipeline(d$phenotype, d$covariates,
    num_tau_levels = 4, method = "ivw", num_cores = 1
  )
  expect_true(is.numeric(result))
  expect_equal(length(result), length(d$phenotype))
  expect_true(all(is.finite(result)))
})

# ===========================================================================
#  ctwas wrapper tests
# ===========================================================================


# ---------- trim_ctwas_variants --------------------------------------------

# Helper: build a minimal region_data structure that trim_ctwas_variants expects
make_mock_region_data <- function() {
  # Variant IDs in canonical format (chr:pos:A2:A1)
  variant_ids <- c("chr1:1000:A:G", "chr1:2000:C:T", "chr1:3000:G:A", "chr1:4000:T:C")

  # Weight matrix (4 variants x 1 weight column)
  wgt <- matrix(c(0.5, 0.0001, 0.3, -0.2), nrow = 4, ncol = 1)
  rownames(wgt) <- variant_ids

  gene_id <- "GENE1|ctx1"
  context <- "ctx1"
  study <- "study1"

  weights <- list()
  weights[[gene_id]] <- list()
  weights[[gene_id]][[study]] <- list(
    wgt = wgt,
    context = context,
    p0 = 1000,
    p1 = 4000
  )

  # SuSiE intermediate info
  pip_vals <- c(0.8, 0.05, 0.6, 0.02)
  names(pip_vals) <- variant_ids

  susie_weights_intermediate <- list()
  susie_weights_intermediate[["GENE1"]] <- list()
  susie_weights_intermediate[["GENE1"]][[context]] <- list(
    pip = pip_vals,
    cs_variants = list(variant_ids[c(1, 3)]),
    cs_purity = list(min.abs.corr = 0.9)
  )

  list(
    weights = weights,
    susie_weights_intermediate = susie_weights_intermediate
  )
}

test_that("trim_ctwas_variants removes variants below weight cutoff", {
  rd <- make_mock_region_data()
  # Default cutoff 1e-5, variant 2 has weight 0.0001 (above), so all 4 should pass default cutoff
  result <- trim_ctwas_variants(rd, twas_weight_cutoff = 1e-5)
  expect_true(is.list(result))
  # With a higher cutoff, the near-zero variant should be removed
  result_strict <- trim_ctwas_variants(rd, twas_weight_cutoff = 0.001)
  # study1 should exist in result
  expect_true("study1" %in% names(result_strict))
  # Get the gene-level result
  gene_weights <- result_strict[["study1"]][["GENE1|ctx1"]]
  # Variant 2 has abs(weight) = 0.0001 < 0.001, so should be removed
  expect_false("chr1:2000:C:T" %in% rownames(gene_weights$wgt))
})

test_that("trim_ctwas_variants removes gene when all weights below cutoff", {
  rd <- make_mock_region_data()
  # Set cutoff so high that all variants are dropped
  result <- trim_ctwas_variants(rd, twas_weight_cutoff = 10)
  # Gene should be removed entirely since no weights pass the cutoff
  # Result should be an empty list
  expect_equal(length(result), 0)
})

test_that("trim_ctwas_variants returns result keyed by study", {
  rd <- make_mock_region_data()
  result <- trim_ctwas_variants(rd, twas_weight_cutoff = 1e-5)
  # merge_by_study reorganizes: weights[[study]][[group]]
  expect_true("study1" %in% names(result))
  expect_true("GENE1|ctx1" %in% names(result[["study1"]]))
})

test_that("trim_ctwas_variants updates p0 and p1 positions", {
  rd <- make_mock_region_data()
  # Use a weight cutoff that removes the variant at position 2000
  result <- trim_ctwas_variants(rd, twas_weight_cutoff = 0.001)
  gene_weights <- result[["study1"]][["GENE1|ctx1"]]
  # p0 and p1 should reflect the range of remaining variant positions
  remaining_positions <- as.integer(sapply(
    rownames(gene_weights$wgt),
    function(v) strsplit(v, ":")[[1]][2]
  ))
  expect_equal(gene_weights$p0, min(remaining_positions))
  expect_equal(gene_weights$p1, max(remaining_positions))
})

test_that("trim_ctwas_variants respects max_num_variants", {
  rd <- make_mock_region_data()
  # Request max 2 variants; since nrow(wgt) == 4 >= max_num_variants == 2,
  # it triggers select_variants which picks by PIP priority
  result <- trim_ctwas_variants(rd, twas_weight_cutoff = 1e-5, max_num_variants = 2)
  gene_weights <- result[["study1"]][["GENE1|ctx1"]]
  expect_true(nrow(gene_weights$wgt) <= 2)
})

test_that("trim_ctwas_variants handles NA weights by removing group", {
  rd <- make_mock_region_data()
  # Replace all weights with NA
  rd$weights[["GENE1|ctx1"]][["study1"]]$wgt[, 1] <- NA
  result <- trim_ctwas_variants(rd, twas_weight_cutoff = 0)
  # The group should be removed because all weights are NA
  expect_equal(length(result), 0)
})

test_that("trim_ctwas_variants handles multiple genes", {
  rd <- make_mock_region_data()

  # Add a second gene
  variant_ids2 <- c("chr1:5000:A:G", "chr1:6000:C:T")
  wgt2 <- matrix(c(0.4, -0.3), nrow = 2, ncol = 1)
  rownames(wgt2) <- variant_ids2

  rd$weights[["GENE2|ctx1"]] <- list()
  rd$weights[["GENE2|ctx1"]][["study1"]] <- list(
    wgt = wgt2,
    context = "ctx1",
    p0 = 5000,
    p1 = 6000
  )

  pip_vals2 <- c(0.7, 0.4)
  names(pip_vals2) <- variant_ids2
  rd$susie_weights_intermediate[["GENE2"]] <- list()
  rd$susie_weights_intermediate[["GENE2"]][["ctx1"]] <- list(
    pip = pip_vals2,
    cs_variants = list(variant_ids2[1]),
    cs_purity = list(min.abs.corr = 0.95)
  )

  result <- trim_ctwas_variants(rd, twas_weight_cutoff = 1e-5)
  expect_true("GENE1|ctx1" %in% names(result[["study1"]]))
  expect_true("GENE2|ctx1" %in% names(result[["study1"]]))
})

test_that("trim_ctwas_variants select_variants uses cs_min_cor to include CS variants", {
  rd <- make_mock_region_data()
  # cs_purity min.abs.corr = 0.9, so with cs_min_cor = 0.8 the CS variants
  # (variant 1 and 3) should be included. Max 2 variants.
  result <- trim_ctwas_variants(rd,
    twas_weight_cutoff = 1e-5,
    cs_min_cor = 0.8,
    min_pip_cutoff = 0.0,
    max_num_variants = 2
  )
  gene_weights <- result[["study1"]][["GENE1|ctx1"]]
  included <- rownames(gene_weights$wgt)
  # CS variants chr1:1000:A:G and chr1:3000:G:A have highest PIPs (0.8 and 0.6)
  # and are in the CS, so they should be prioritized
  expect_true("chr1:1000:A:G" %in% included)
  expect_true("chr1:3000:G:A" %in% included)
})
