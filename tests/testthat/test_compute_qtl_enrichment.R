context("compute_qtl_enrichment")

generate_mock_data <- function(seed=1, num_pips = 1000, num_susie_fits = 2) {
  # Simulate fake data for gwas_pip
  n_gwas_pip <- num_pips
  gwas_pip <- runif(n_gwas_pip)
  names(gwas_pip) <- paste0("snp", 1:n_gwas_pip)
  gwas_fit <- list(pip=gwas_pip)

  # Simulate fake data for a single SuSiEFit object
  simulate_susiefit <- function(n, p) {
    pip <- runif(n)
    names(pip) <- paste0("snp", 1:n)
    alpha <- t(matrix(runif(n * p), nrow = n))
    alpha <- t(apply(alpha, 1, function(row) row / sum(row)))
    list(
      pip = pip,
      alpha = alpha,
      prior_variance = runif(p)
    )
  }
  
  # Simulate multiple SuSiEFit objects
  n_susie_fits <- num_susie_fits
  susie_fits <- replicate(n_susie_fits, simulate_susiefit(n_gwas_pip, 10), simplify = FALSE)
  # Add these fits to a list, providing names to each element
  names(susie_fits) <- paste0("fit", 1:length(susie_fits))
  return(list(gwas_fit=gwas_fit, susie_fits=susie_fits))
}

test_that("compute_qtl_enrichment dummy data single-threaded works",{
  local_mocked_bindings(
      qtl_enrichment_rcpp = function(...) TRUE)
  input_data <- generate_mock_data(seed=1, num_pips=10)
  expect_warning(
    compute_qtl_enrichment(input_data$gwas_fit$pip, input_data$susie_fits, lambda = 1, ImpN = 10, num_threads = 1),
    "num_gwas is not provided. Estimating pi_gwas from the data. Note that this estimate may be biased if the input gwas_pip does not contain genome-wide variants.")
  expect_warning(
    compute_qtl_enrichment(input_data$gwas_fit$pip, input_data$susie_fits, lambda = 1, ImpN = 10, num_threads = 1),
    "pi_qtl is not provided. Estimating pi_qtl from the data. Note that this estimate may be biased if either 1) the input susie_qtl_regions does not have enough data, or 2) the single effects only include variables inside of credible sets or signal clusters.")
  res <- expect_warning(compute_qtl_enrichment(input_data$gwas_fit$pip, input_data$susie_fits, num_gwas=5000, pi_qtl=0.49819, lambda = 1, ImpN = 10, num_threads = 1))
  expect_true(length(res) > 0)
})

test_that("compute_qtl_enrichment dummy data single thread and multi-threaded are equivalent",{
  local_mocked_bindings(
      qtl_enrichment_rcpp = function(...) TRUE)
  input_data <- generate_mock_data(seed=1, num_pips=10)
  res_single <- expect_warning(compute_qtl_enrichment(input_data$gwas_fit$pip, input_data$susie_fits, num_gwas=5000, pi_qtl=0.49819, lambda = 1, ImpN = 10, num_threads = 1))
  res_multi <- expect_warning(compute_qtl_enrichment(input_data$gwas_fit$pip, input_data$susie_fits, num_gwas=5000, pi_qtl=0.49819, lambda = 1, ImpN = 10, num_threads = 2))
  expect_equal(res_single, res_multi)
})

# ---- error paths (compute_qtl_enrichment.R lines 86, 87, 91) ----
test_that("compute_qtl_enrichment errors when pi_gwas is zero", {
  gwas_pip <- rep(0, 10)
  names(gwas_pip) <- paste0("snp", 1:10)
  susie_fits <- list(fit1 = list(pip = setNames(runif(10), paste0("snp", 1:10)),
                                  alpha = matrix(1, 1, 10),
                                  prior_variance = 1))
  expect_error(
    compute_qtl_enrichment(gwas_pip, susie_fits, pi_qtl = 0.5),
    "No association signal found in GWAS data"
  )
})

test_that("compute_qtl_enrichment errors when pi_qtl is zero", {
  gwas_pip <- runif(10)
  names(gwas_pip) <- paste0("snp", 1:10)
  susie_fits <- list(fit1 = list(pip = setNames(rep(0, 10), paste0("snp", 1:10)),
                                  alpha = matrix(1, 1, 10),
                                  prior_variance = 1))
  expect_error(
    suppressWarnings(compute_qtl_enrichment(gwas_pip, susie_fits, num_gwas = 1000, pi_qtl = 0)),
    "No QTL associated"
  )
})

test_that("compute_qtl_enrichment errors when gwas_pip has no names", {
  gwas_pip <- runif(10)  # no names
  susie_fits <- list(fit1 = list(pip = setNames(runif(10), paste0("snp", 1:10)),
                                  alpha = matrix(1, 1, 10),
                                  prior_variance = 1))
  expect_error(
    suppressWarnings(compute_qtl_enrichment(gwas_pip, susie_fits, num_gwas = 1000, pi_qtl = 0.5)),
    "Variant names are missing in gwas_pip"
  )
})

# ---- real C++ qtl_enrichment_rcpp integration test ----
test_that("compute_qtl_enrichment calls real C++ enrichment code and returns expected keys", {
  set.seed(42)
  n_snps <- 50
  variant_names <- paste0("1:", 1:n_snps, ":A:G")

  # GWAS PIPs: sparse signal
  gwas_pip <- rep(0.01, n_snps)
  gwas_pip[c(5, 20, 35)] <- c(0.8, 0.6, 0.9)
  names(gwas_pip) <- variant_names

  # SuSiE fit with 2 single effects over same variants
  L <- 2
  alpha <- matrix(1 / n_snps, nrow = L, ncol = n_snps)
  # Concentrate probability on causal variants
  alpha[1, ] <- 0.001; alpha[1, 5] <- 0.95; alpha[1, ] <- alpha[1, ] / sum(alpha[1, ])
  alpha[2, ] <- 0.001; alpha[2, 20] <- 0.95; alpha[2, ] <- alpha[2, ] / sum(alpha[2, ])
  pip <- colSums(alpha)
  names(pip) <- variant_names

  susie_fits <- list(
    fit1 = list(pip = pip, alpha = alpha, prior_variance = c(0.5, 0.3))
  )

  # Call without mocking — exercises the real C++ code
  res <- suppressWarnings(
    compute_qtl_enrichment(gwas_pip, susie_fits,
                           num_gwas = 5000, pi_qtl = 0.5,
                           lambda = 1, ImpN = 5, num_threads = 1)
  )
  expect_type(res, "list")
  # The enrichment results are in res[[1]] (the C++ output list)
  en <- res[[1]]
  expected_keys <- c("Intercept", "Enrichment (no shrinkage)", "Enrichment (w/ shrinkage)",
                     "sd (no shrinkage)", "sd (w/ shrinkage)",
                     "Alternative (coloc) p1", "Alternative (coloc) p2", "Alternative (coloc) p12")
  for (key in expected_keys) {
    expect_true(key %in% names(en), info = paste("Missing key:", key))
  }
  # All numeric and finite
  numeric_vals <- unlist(en[expected_keys])
  expect_true(all(is.finite(numeric_vals)))
})

# ---- unmatched variants tracking (compute_qtl_enrichment.R line 102) ----
test_that("compute_qtl_enrichment tracks unmatched QTL variants", {
  local_mocked_bindings(
    qtl_enrichment_rcpp = function(...) TRUE
  )
  gwas_pip <- runif(10)
  names(gwas_pip) <- paste0("1:", 1:10, ":A:G")
  # QTL has some variants not in GWAS
  qtl_pip <- runif(5)
  names(qtl_pip) <- c(paste0("1:", 1:3, ":A:G"), "1:999:A:G", "1:998:A:G")
  susie_fits <- list(fit1 = list(pip = qtl_pip,
                                  alpha = matrix(runif(5), 1, 5),
                                  prior_variance = 1))
  res <- suppressWarnings(
    compute_qtl_enrichment(gwas_pip, susie_fits, num_gwas = 1000, pi_qtl = 0.5)
  )
  expect_true("unused_xqtl_variants" %in% names(res))
})