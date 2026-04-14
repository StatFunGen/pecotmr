context("raiss")
library(tidyverse)
library(MASS)

generate_dummy_data <- function(seed=1, ref_panel_ordered=TRUE, known_zscores_ordered=TRUE) {
    set.seed(seed)

    n_variants <- 100
    ref_panel <- data.frame(
        chrom = rep(1, n_variants),
        pos = seq(1, n_variants * 10, 10),
        variant_id = paste0("rs", seq_len(n_variants)),
        A1 = sample(c("A", "T", "G", "C"), n_variants, replace = TRUE),
        A2 = sample(c("A", "T", "G", "C"), n_variants, replace = TRUE)
    )

    n_known <- 50
    known_zscores <- data.frame(
        chrom = rep(1, n_known),
        pos = sample(ref_panel$pos, n_known),
        variant_id = sample(ref_panel$variant_id, n_known),
        A1 = sample(c("A", "T", "G", "C"), n_known, replace = TRUE),
        A2 = sample(c("A", "T", "G", "C"), n_known, replace = TRUE),
        z = rnorm(n_known)
    )

    LD_matrix <- matrix(rnorm(n_variants^2), nrow = n_variants, ncol = n_variants)
    diag(LD_matrix) <- 1 
    known_zscores <- if (known_zscores_ordered) known_zscores[order(known_zscores$pos),] else known_zscores
    ref_panel <- if (ref_panel_ordered) ref_panel else ref_panel[order(ref_panel$pos, decreasing = TRUE),]
    return(list(ref_panel=ref_panel, known_zscores=known_zscores, LD_matrix=LD_matrix))
}

test_that("Input validation for raiss works correctly", {
    input_data <- generate_dummy_data()
    input_data_ref_panel_unordered <- generate_dummy_data(ref_panel_ordered=FALSE)
    input_data_zscores_unordered <- generate_dummy_data(known_zscores_ordered=FALSE)
    expect_error(raiss(input_data_ref_panel_unordered$ref_panel, input_data$known_zscores, input_data$LD_matrix))
    expect_error(raiss(input_data$ref_panel, input_data_zscores_unordered$known_zscores, input_data$LD_matrix))
})

test_that("Default parameters for raiss work correctly", {
    input_data <- generate_dummy_data()
    result <- raiss(input_data$ref_panel, input_data$known_zscores, input_data$LD_matrix)
    expect_true(is.list(result))
    # Expected list elements
    expect_true(all(c("result_nofilter", "result_filter", "LD_mat") %in% names(result)))
    # result_nofilter should be a data frame with expected columns
    expect_true(is.data.frame(result$result_nofilter))
    expect_true(all(c("variant_id", "z", "Var", "raiss_ld_score") %in% names(result$result_nofilter)))
    # Imputed z-scores should be numeric and finite
    expect_true(is.numeric(result$result_nofilter$z))
    expect_true(all(is.finite(result$result_nofilter$z)))
    # Output should cover all ref_panel variants (known + imputed)
    expect_equal(nrow(result$result_nofilter), nrow(input_data$ref_panel))
    # Filtered result should be a subset of unfiltered
    expect_true(nrow(result$result_filter) <= nrow(result$result_nofilter))
    # LD_mat should be a matrix
    expect_true(is.matrix(result$LD_mat))
})

test_that("Test Default Parameters for raiss_model", {
  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

  result <- raiss_model(zt, sig_t, sig_i_t)

  expect_true(is.list(result))
  expect_true(all(c("var", "mu", "raiss_ld_score", "condition_number", "correct_inversion") %in% names(result)))
  # mu (imputed z-scores) should be numeric and finite
  expect_true(is.numeric(result$mu))
  expect_true(all(is.finite(result$mu)))
  # var should be numeric
  expect_true(is.numeric(result$var))
  # raiss_ld_score should be numeric and non-negative
  expect_true(is.numeric(result$raiss_ld_score))
  expect_true(all(result$raiss_ld_score >= 0))
})

test_that("Test with Different lamb Values for raiss_model", {
  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

  lamb_values <- c(0.01, 0.05, 0.1)
  for (lamb in lamb_values) {
    result <- raiss_model(zt, sig_t, sig_i_t, lamb)
    expect_true(is.list(result))
    expect_true(all(c("var", "mu", "raiss_ld_score") %in% names(result)))
    expect_true(is.numeric(result$mu))
    expect_true(all(is.finite(result$mu)))
  }
})

test_that("Report Condition Number in raiss_model", {
  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

  result_with_cn <- raiss_model(zt, sig_t, sig_i_t, report_condition_number = TRUE)
  result_without_cn <- raiss_model(zt, sig_t, sig_i_t, report_condition_number = FALSE)

  expect_true(is.list(result_with_cn))
  expect_true(is.list(result_without_cn))
  # With condition number reporting, condition_number should be populated
  expect_true(is.numeric(result_with_cn$condition_number))
  expect_true(all(is.finite(result_with_cn$mu)))
  expect_true(all(is.finite(result_without_cn$mu)))
})

test_that("Input Validation of raiss_model", {

  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)
  zt_invalid <- "not a numeric vector"
  sig_t_invalid <- "not a matrix"
  sig_i_t_invalid <- "not a matrix"

  expect_error(raiss_model(zt_invalid, sig_t, sig_i_t))
  expect_error(raiss_model(zt, sig_t_invalid, sig_i_t))
  expect_error(raiss_model(zt, sig_t, sig_i_t_invalid))
})

test_that("Boundary Conditions of raiss_model", {
  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

  zt_empty <- numeric(0)
  sig_t_empty <- matrix(numeric(0), nrow = 0)
  sig_i_t_empty <- matrix(numeric(0), nrow = 0)

  expect_error(raiss_model(zt_empty, sig_t, sig_i_t))
  expect_error(raiss_model(zt, sig_t_empty, sig_i_t))
  expect_error(raiss_model(zt, sig_t, sig_i_t_empty))
})

test_that("Test with Different rcond Values for raiss_model", {
  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

  rcond_values <- c(0.01, 0.05, 0.1)
  for (rcond in rcond_values) {
    result <- raiss_model(zt, sig_t, sig_i_t, lamb = 0.01, rcond = rcond)
    expect_true(is.list(result))
    expect_true(all(c("var", "mu", "raiss_ld_score", "condition_number", "correct_inversion") %in% names(result)))
    expect_true(is.numeric(result$mu))
    expect_true(all(is.finite(result$mu)))
  }
})

test_that("format_raiss_df returns correctly formatted data frame", {
  imp <- list(
    mu = rnorm(5),
    var = runif(5),
    raiss_ld_score = rnorm(5),
    condition_number = runif(5),
    correct_inversion = sample(c(TRUE, FALSE), 5, replace = TRUE)
  )
  
  ref_panel <- data.frame(
    chrom = sample(1:22, 10, replace = TRUE),
    pos = sample(1:10000, 10),
    variant_id = paste0("rs", 1:10),
    A1 = sample(c("A", "T", "G", "C"), 10, replace = TRUE),
    A2 = sample(c("A", "T", "G", "C"), 10, replace = TRUE)
  )

  unknowns <- sample(1:nrow(ref_panel), 5)

  result <- format_raiss_df(imp, ref_panel, unknowns)

  expect_true(is.data.frame(result))
  expect_equal(ncol(result), 10)
  expect_equal(colnames(result), c('chrom', 'pos', 'variant_id', 'A1', 'A2', 'z', 'Var', 'raiss_ld_score', 'condition_number', 'correct_inversion'))

  for (col in c('chrom', 'pos', 'variant_id', 'A1', 'A2')) {
    expect_equal(setNames(unlist(result[col]), NULL), unlist(ref_panel[unknowns, col, drop = TRUE]))
  }
  for (col in c('z', 'Var', 'raiss_ld_score', 'condition_number', 'correct_inversion')) {
    expected_col <- if (col == "z") "mu" else if (col == "Var") "var" else col
    expect_equal(setNames(unlist(result[col]), NULL), setNames(unlist(imp[expected_col]), NULL))
  }
})

test_that("Merge operation is correct for merge_raiss_df", {
    raiss_df_example <- data.frame(
        chrom = c("chr21", "chr22"),
        pos = c(123, 456),
        variant_id = c("var1", "var2"),
        A1 = c("A", "T"),
        A2 = c("T", "A"),
        z = c(0.5, 1.5),
        Var = c(0.2, 0.3),
        raiss_ld_score = c(10, 20),
        raiss_R2 = c(0.8, 0.7))

    known_zscores_example <- data.frame(
        chrom = c("chr21", "chr22"),
        pos = c(123, 456),
        variant_id = c("var1", "var2"),
        A1 = c("A", "T"),
        A2 = c("T", "A"),
        z = c(0.5, 1.5))

    merged_df <- merge_raiss_df(raiss_df_example, known_zscores_example)
    expect_equal(nrow(merged_df), 2)
    expect_true(all(c("chr21", "chr22") %in% merged_df$chrom))
})

generate_fro_test_data <- function(seed=1) {
    set.seed(seed)
    return(data.frame(
        chrom = paste0("chr", rep(22, 10)),
        pos = seq(1, 100, 10),
        variant_id = 1:10,
        A1 = rep("A", 10),
        A2 = rep("T", 10),
        z = rnorm(10),
        Var = runif(10, 0, 1),
        raiss_ld_score = rnorm(10, 5, 2)
    ))
}

test_that("Correct columns are selected in filter_raiss_output", {
    test_data <- generate_fro_test_data()
    output <- filter_raiss_output(test_data)$zscores
    expect_true(all(c('variant_id', 'A1', 'A2', 'z', 'Var', 'raiss_ld_score') %in% names(output)))
})

test_that("raiss_R2 is calculated correctly in filter_raiss_output", {
    test_data <- generate_fro_test_data()
    output <- filter_raiss_output(test_data)$zscores
    expected_R2 <- 1 - test_data[which(test_data$raiss_ld_score >= 5),]$Var
    expect_equal(output$raiss_R2, expected_R2[which(expected_R2 > 0.6)])
})

test_that("Filtering is applied correctly in filter_raiss_output", {
    test_data <- generate_fro_test_data()
    R2_threshold <- 0.6
    minimum_ld <- 5
    output <- filter_raiss_output(test_data, R2_threshold, minimum_ld)$zscores

    expect_true(all(output$raiss_R2 > R2_threshold))
    expect_true(all(output$raiss_ld_score >= minimum_ld))
})

test_that("Function returns the correct subset in filter_raiss_output", {
    test_data <- generate_fro_test_data()
    test_data$raiss_R2 <- 1 - test_data$Var
    output <- filter_raiss_output(test_data)$zscores

    manual_filter <- test_data[test_data$raiss_R2 > 0.6 & test_data$raiss_ld_score >= 5, ]

    expect_equal(nrow(output), nrow(manual_filter))
    expect_equal(sum(output$variant_id != manual_filter$variant_id), 0)
})

test_that("compute_mu basic functionality", {
    sig_i_t <- matrix(c(1, 2, 3, 4), nrow = 2)
    sig_t_inv <- matrix(c(5, 6, 7, 8), nrow = 2)
    zt <- matrix(c(9, 10, 11, 12), nrow = 2)

    expected_result <- matrix(c(517, 766, 625, 926), nrow = 2)
    result <- compute_mu(sig_i_t, sig_t_inv, zt)
    expect_equal(result, expected_result)
})

generate_mock_data_for_compute_var <- function(seed=1) {
    return(
        list(
            sig_i_t_1 = matrix(c(1, 2, 3, 4), nrow = 2),
            sig_t_inv_1 = matrix(c(5, 6, 7, 8), nrow = 2),
            lamb_1 = 0.5))
}

test_that("compute_var returns correct output for batch = TRUE", {
    input_data <- generate_mock_data_for_compute_var()
    result <- compute_var(input_data$sig_i_t_1, input_data$sig_t_inv_1, input_data$lamb_1, batch = TRUE)
    expect_true(is.list(result))
    expect_length(result, 2)
    expect_true(all(c("var", "raiss_ld_score") %in% names(result)))
    expect_true(is.numeric(result$var))
    expect_true(is.numeric(result$raiss_ld_score))
})

test_that("compute_var returns correct output for batch = FALSE", {
    input_data <- generate_mock_data_for_compute_var()
    result <- compute_var(input_data$sig_i_t_1, input_data$sig_t_inv_1, input_data$lamb_1, batch = FALSE)
    expect_true(is.list(result))
    expect_length(result, 2)
    expect_true(all(c("var", "raiss_ld_score") %in% names(result)))
    expect_true(is.numeric(result$var))
    expect_true(is.numeric(result$raiss_ld_score))
})

test_that("check_inversion correctly identifies inverse matrices in", {
  sig_t <- matrix(c(1, 2, 3, 4), nrow=2, ncol=2)
  sig_t_inv <- solve(sig_t)  
  expect_true(check_inversion(sig_t, sig_t_inv))
})

test_that("var_in_boundaries sets boundaries correctly", {
  lamb_test <- 0.05
  var <- c(-1, 0, 0.5, 1.04, 1.05)  

  result <- var_in_boundaries(var, lamb_test)

  expect_equal(result[1], 0)                   # Value less than 0 should be set to 0
  expect_equal(result[2], 0)                   # Value within lower boundary should remain unchanged
  expect_equal(result[3], 0.5)                 # Value within boundaries should remain unchanged
  expect_equal(result[4], 1.04)                   # Value greater than 0.99999 + lamb should be set to 1
  expect_equal(result[5], 1)                   # Value greater than 0.99999 + lamb should be set to 1
})

test_that("invert_mat computes correct pseudo-inverse", {
  mat <- matrix(c(1, 2, 3, 4), nrow = 2)
  lamb <- 0.5
  rcond <- 1e-7
  result <- invert_mat(mat, lamb, rcond)
  expect_true(is.matrix(result))
})

test_that("invert_mat handles errors and retries", {
  mat <- matrix(c(0, 0, 0, 0), nrow = 2) 
  lamb <- 0.1
  rcond <- 1e-7
  result <- invert_mat(mat, lamb, rcond)
  expect_true(is.matrix(result))
})

test_that("invert_mat_recursive correctly inverts a valid square matrix", {
  mat <- matrix(c(2, -1, -1, 2), nrow = 2)
  lamb <- 0.5
  rcond <- 0.01
  result <- invert_mat_recursive(mat, lamb, rcond)
  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(mat))
})

test_that("invert_mat_recursive handles non-square matrices appropriately", {
  mat <- matrix(1:6, nrow = 2)
  lamb <- 0.5
  rcond <- 0.01
  expect_silent(invert_mat_recursive(mat, lamb, rcond))
})

test_that("invert_mat_recursive handles errors and performs recursive call correctly", {
  mat <- "not a matrix"
  lamb <- 0.5
  rcond <- 0.01
  expect_error(invert_mat_recursive(mat, lamb, rcond))
})

# Test with Different Tolerance Levels
test_that("invert_mat_eigen behaves differently with varying tolerance levels", {
  mat <- matrix(c(1, 0, 0, 1e-4), nrow = 2)
  tol_high <- 1e-2
  tol_low <- 1e-6
  result_high_tol <- invert_mat_eigen(mat, tol_high)
  result_low_tol <- invert_mat_eigen(mat, tol_low)
  expect_true(!is.logical(all.equal(result_high_tol, result_low_tol)))
})

test_that("invert_mat_eigen handles non-square matrices", {
  mat <- matrix(1:6, nrow = 2)
  expect_error(invert_mat_eigen(mat))
})

test_that("invert_mat_eigen returns the same matrix for an identity matrix", {
    mat <- diag(2)
    expected <- mat
    actual <- invert_mat_eigen(mat)
    expect_equal(actual, expected)
})

test_that("invert_mat_eigen returns a zero matrix for a zero matrix input", {
    mat <- matrix(0, nrow = 2, ncol = 2)
    expected <- mat
    expect_error(invert_mat_eigen(mat),
      "Cannot invert the input matrix because all its eigen values are negative or close to zero")
})

test_that("invert_mat_eigen handles matrices with negative eigenvalues", {
    mat <- matrix(c(-2, 0, 0, -3), nrow = 2)
    expect_silent(invert_mat_eigen(mat))
})

# ===========================================================================
# raiss_single_matrix edge cases
# ===========================================================================

test_that("raiss_single_matrix returns NULL when no known variants overlap", {
  set.seed(42)
  ref_panel <- data.frame(
    chrom = rep(1, 10), pos = seq(10, 100, 10),
    variant_id = paste0("rs", 1:10),
    A1 = rep("A", 10), A2 = rep("G", 10),
    stringsAsFactors = FALSE
  )
  # known_zscores has variant IDs that don't match ref_panel at all
  known_zscores <- data.frame(
    chrom = rep(1, 3), pos = c(200, 300, 400),
    variant_id = paste0("other", 1:3),
    A1 = rep("A", 3), A2 = rep("G", 3),
    z = rnorm(3), stringsAsFactors = FALSE
  )
  LD_matrix <- diag(10)
  result <- raiss_single_matrix(ref_panel, known_zscores, LD_matrix, verbose = FALSE)
  expect_null(result)
})

test_that("raiss_single_matrix returns known zscores when no unknowns to impute", {
  set.seed(42)
  ref_panel <- data.frame(
    chrom = rep(1, 5), pos = seq(10, 50, 10),
    variant_id = paste0("rs", 1:5),
    A1 = rep("A", 5), A2 = rep("G", 5),
    stringsAsFactors = FALSE
  )
  # All ref_panel variants are known — nothing to impute
  known_zscores <- data.frame(
    chrom = rep(1, 5), pos = seq(10, 50, 10),
    variant_id = paste0("rs", 1:5),
    A1 = rep("A", 5), A2 = rep("G", 5),
    z = rnorm(5), stringsAsFactors = FALSE
  )
  LD_matrix <- diag(5)
  result <- raiss_single_matrix(ref_panel, known_zscores, LD_matrix, verbose = FALSE)
  expect_true(is.list(result))
  expect_equal(result$result_nofilter, known_zscores)
  expect_equal(result$result_filter, known_zscores)
  expect_equal(result$LD_mat, LD_matrix)
})

# ===========================================================================
# raiss_single_matrix_from_X edge cases
# ===========================================================================

test_that("raiss_single_matrix_from_X returns NULL when no known variants overlap", {
  set.seed(42)
  n <- 50
  p <- 10
  ref_panel <- data.frame(
    chrom = rep(1, p), pos = seq(10, p * 10, 10),
    variant_id = paste0("rs", 1:p),
    A1 = rep("A", p), A2 = rep("G", p),
    stringsAsFactors = FALSE
  )
  known_zscores <- data.frame(
    chrom = rep(1, 3), pos = c(200, 300, 400),
    variant_id = paste0("other", 1:3),
    A1 = rep("A", 3), A2 = rep("G", 3),
    z = rnorm(3), stringsAsFactors = FALSE
  )
  X <- scale(matrix(sample(0:2, n * p, replace = TRUE), nrow = n))
  X[is.na(X)] <- 0
  colnames(X) <- ref_panel$variant_id
  result <- raiss_single_matrix_from_X(ref_panel, known_zscores, X, verbose = FALSE)
  expect_null(result)
})

test_that("raiss_single_matrix_from_X returns known zscores when no unknowns to impute", {
  set.seed(42)
  n <- 50
  p <- 5
  ref_panel <- data.frame(
    chrom = rep(1, p), pos = seq(10, p * 10, 10),
    variant_id = paste0("rs", 1:p),
    A1 = rep("A", p), A2 = rep("G", p),
    stringsAsFactors = FALSE
  )
  known_zscores <- data.frame(
    chrom = rep(1, p), pos = seq(10, p * 10, 10),
    variant_id = paste0("rs", 1:p),
    A1 = rep("A", p), A2 = rep("G", p),
    z = rnorm(p), stringsAsFactors = FALSE
  )
  X <- scale(matrix(sample(0:2, n * p, replace = TRUE), nrow = n))
  X[is.na(X)] <- 0
  colnames(X) <- ref_panel$variant_id
  result <- raiss_single_matrix_from_X(ref_panel, known_zscores, X, verbose = FALSE)
  expect_true(is.list(result))
  expect_equal(result$result_nofilter, known_zscores)
  expect_null(result$LD_mat)
})

# ===========================================================================
# raiss() dispatch paths: single-matrix list and genotype_matrix list
# ===========================================================================

test_that("raiss with single-matrix LD list dispatches to single matrix path", {
  set.seed(42)
  n_variants <- 20
  ref_panel <- data.frame(
    chrom = rep(1, n_variants), pos = seq(10, n_variants * 10, 10),
    variant_id = paste0("rs", 1:n_variants),
    A1 = rep("A", n_variants), A2 = rep("G", n_variants),
    stringsAsFactors = FALSE
  )
  n_known <- 10
  known_idx <- sort(sample(seq_len(n_variants), n_known))
  known_zscores <- data.frame(
    chrom = rep(1, n_known), pos = ref_panel$pos[known_idx],
    variant_id = ref_panel$variant_id[known_idx],
    A1 = ref_panel$A1[known_idx], A2 = ref_panel$A2[known_idx],
    z = rnorm(n_known), stringsAsFactors = FALSE
  )
  R <- diag(n_variants)
  colnames(R) <- rownames(R) <- ref_panel$variant_id

  # Wrap in list structure with ld_matrices
  LD_list <- list(ld_matrices = list(R))

  result <- raiss(ref_panel, known_zscores, LD_matrix = LD_list,
                  R2_threshold = 0, minimum_ld = 0, verbose = FALSE)
  expect_true(is.list(result))
  expect_true("result_nofilter" %in% names(result))
  expect_equal(nrow(result$result_nofilter), n_variants)
})

test_that("raiss with genotype_matrix list processes multiple blocks", {
  set.seed(42)
  n <- 50
  p <- 20
  # Each block has its own ref_panel subset matching the X columns
  ref_panel <- data.frame(
    chrom = rep(1, p), pos = seq(10, p * 10, 10),
    variant_id = paste0("rs", 1:p),
    A1 = rep("A", p), A2 = rep("G", p),
    stringsAsFactors = FALSE
  )
  # Use only first block's variants as known (so second block has unknowns)
  known_idx <- sort(sample(1:10, 5))
  known_zscores <- data.frame(
    chrom = rep(1, length(known_idx)), pos = ref_panel$pos[known_idx],
    variant_id = ref_panel$variant_id[known_idx],
    A1 = ref_panel$A1[known_idx], A2 = ref_panel$A2[known_idx],
    z = rnorm(length(known_idx)), stringsAsFactors = FALSE
  )
  X <- scale(matrix(sample(0:2, n * p, replace = TRUE), nrow = n))
  X[is.na(X)] <- 0
  colnames(X) <- ref_panel$variant_id

  # Use the full matrix as a single-element list — the simplest valid list input
  X_list <- list(X)

  result <- raiss(ref_panel, known_zscores, genotype_matrix = X_list,
                  R2_threshold = 0, minimum_ld = 0, verbose = FALSE)
  expect_true(is.list(result))
  expect_true("result_nofilter" %in% names(result))
  expect_true(nrow(result$result_nofilter) > 0)
  expect_null(result$LD_mat)
})

test_that("raiss with genotype_matrix list returns NULL when all blocks fail", {
  set.seed(42)
  n <- 50
  p <- 10
  ref_panel <- data.frame(
    chrom = rep(1, p), pos = seq(10, p * 10, 10),
    variant_id = paste0("rs", 1:p),
    A1 = rep("A", p), A2 = rep("G", p),
    stringsAsFactors = FALSE
  )
  # known_zscores has no overlap with ref_panel
  known_zscores <- data.frame(
    chrom = rep(1, 3), pos = c(200, 300, 400),
    variant_id = paste0("other", 1:3),
    A1 = rep("A", 3), A2 = rep("G", 3),
    z = rnorm(3), stringsAsFactors = FALSE
  )
  X <- scale(matrix(sample(0:2, n * p, replace = TRUE), nrow = n))
  X[is.na(X)] <- 0
  colnames(X) <- ref_panel$variant_id
  X_list <- list(X[, 1:5, drop = FALSE], X[, 6:10, drop = FALSE])

  result <- raiss(ref_panel, known_zscores, genotype_matrix = X_list,
                  verbose = FALSE)
  expect_null(result)
})

# Block-Diagonal LD data generator for RAISS testing
# Corrected function to generate proper block-diagonal test data
generate_block_diagonal_test_data <- function(seed = 123, block_structure = "overlapping", n_variants = 30) {
  set.seed(seed)
  
  # Create reference panel with variants
  ref_panel <- data.frame(
    chrom = rep(1, n_variants),
    pos = seq(1, n_variants * 10, 10),
    variant_id = paste0("var", seq_len(n_variants)),
    A1 = sample(c("A", "T", "G", "C"), n_variants, replace = TRUE),
    A2 = sample(c("A", "T", "G", "C"), n_variants, replace = TRUE),
    stringsAsFactors = FALSE
  )
  
  # Create known z-scores for every other variant
  known_indices <- seq(1, n_variants, by = 2)
  known_zscores <- data.frame(
    chrom = rep(1, length(known_indices)),
    pos = ref_panel$pos[known_indices],
    variant_id = ref_panel$variant_id[known_indices],
    A1 = ref_panel$A1[known_indices],
    A2 = ref_panel$A2[known_indices],
    z = rnorm(length(known_indices)),
    stringsAsFactors = FALSE
  )
  
  # Define block boundaries based on requested structure
  if (block_structure == "overlapping") {
    block_boundaries <- list(
      c(1, 11),    # Block 1: variants 1-11
      c(11, 21),   # Block 2: variants 11-21 (overlap at var11)
      c(21, n_variants)  # Block 3: variants 21-30 (overlap at var21)
    )
  } else if (block_structure == "non_overlapping") {
    block_boundaries <- list(
      c(1, 10),
      c(11, 20),
      c(21, n_variants)
    )
  } else if (block_structure == "uneven") {
    block_boundaries <- list(
      c(1, 5),
      c(6, 20),
      c(21, n_variants)
    )
  } else if (block_structure == "many_small") {
    block_size <- 5
    n_blocks <- ceiling(n_variants / block_size)
    block_boundaries <- list()
    for (i in 1:n_blocks) {
      start_idx <- (i-1) * block_size + 1
      end_idx <- min(i * block_size, n_variants)
      block_boundaries[[i]] <- c(start_idx, end_idx)
    }
  } else if (block_structure == "single_block") {
    block_boundaries <- list(c(1, n_variants))
  }
  
  # First, create independent block matrices
  block_matrices <- list()
  for (i in seq_along(block_boundaries)) {
    start_idx <- block_boundaries[[i]][1]
    end_idx <- block_boundaries[[i]][2]
    block_variant_ids <- ref_panel$variant_id[start_idx:end_idx]
    n_block <- length(block_variant_ids)
    
    # Create the block matrix with correlations ONLY within the block
    block_matrix <- matrix(0, nrow = n_block, ncol = n_block)
    for (a in 1:n_block) {
      for (b in 1:n_block) {
        if (a == b) {
          block_matrix[a, b] <- 1
        } else {
          # Use positions within the block, not absolute positions
          block_matrix[a, b] <- 0.95^abs(a - b)
        }
      }
    }
    rownames(block_matrix) <- block_variant_ids
    colnames(block_matrix) <- block_variant_ids
    
    block_matrices[[i]] <- block_matrix
  }
  
  # Create variant indices data frame
  variant_indices <- data.frame(
    variant_id = character(),
    block_id = integer(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(block_boundaries)) {
    start_idx <- block_boundaries[[i]][1]
    end_idx <- block_boundaries[[i]][2]
    block_variant_ids <- ref_panel$variant_id[start_idx:end_idx]
    
    block_indices <- data.frame(
      variant_id = block_variant_ids,
      block_id = i,
      stringsAsFactors = FALSE
    )
    variant_indices <- rbind(variant_indices, block_indices)
  }
  
  # Create block metadata
  block_sizes <- sapply(block_boundaries, function(b) b[2] - b[1] + 1)
  block_metadata <- data.frame(
    block_id = seq_along(block_boundaries),
    chrom = rep(1, length(block_boundaries)),
    size = block_sizes,
    start_idx = sapply(seq_along(block_boundaries), function(i) {
      # Adjust for 1-based indexing in R
      if (i == 1) return(1)
      # Count unique variants before this block
      sum(sapply(1:(i-1), function(j) {
        # If there's an overlap with the next block, count one less
        if (j < length(block_boundaries) && 
            block_boundaries[[j]][2] == block_boundaries[[j+1]][1]) {
          return(block_boundaries[[j]][2] - block_boundaries[[j]][1])
        } else {
          return(block_boundaries[[j]][2] - block_boundaries[[j]][1] + 1)
        }
      })) + 1
    }),
    end_idx = sapply(seq_along(block_boundaries), function(i) {
      # Count all unique variants up to and including this block
      sum(sapply(1:i, function(j) {
        # If there's an overlap with the next block, count one less
        if (j < i && block_boundaries[[j]][2] == block_boundaries[[j+1]][1]) {
          return(block_boundaries[[j]][2] - block_boundaries[[j]][1])
        } else {
          return(block_boundaries[[j]][2] - block_boundaries[[j]][1] + 1)
        }
      }))
    }),
    stringsAsFactors = FALSE
  )
  
  # Build the full matrix correctly ensuring proper block structure
  # IMPORTANT: Initialize a matrix with zeros - ensure no correlations between blocks
  all_variant_ids <- unique(variant_indices$variant_id)
  LD_matrix_full <- matrix(0, nrow = length(all_variant_ids), ncol = length(all_variant_ids))
  rownames(LD_matrix_full) <- all_variant_ids
  colnames(LD_matrix_full) <- all_variant_ids
  
  # For each block, fill in only the relevant section of the full matrix
  for (i in seq_along(block_matrices)) {
    block_matrix <- block_matrices[[i]]
    block_vars <- rownames(block_matrix)
    
    for (var_a in block_vars) {
      for (var_b in block_vars) {
        LD_matrix_full[var_a, var_b] <- block_matrix[var_a, var_b]
      }
    }
  }
  
  # Create the block structure for RAISS
  LD_matrix_blocks <- list(
    ld_matrices = block_matrices,
    variant_indices = variant_indices,
    block_metadata = block_metadata,
    LD_variants = all_variant_ids
  )
  
  return(list(
    ref_panel = ref_panel,
    known_zscores = known_zscores,
    LD_matrix_full = LD_matrix_full,
    LD_matrix_blocks = LD_matrix_blocks,
    variant_indices = variant_indices,
    block_boundaries = block_boundaries,
    block_metadata = block_metadata
  ))
}

test_that("full matrix and block processing produce identical results", {
  # Only test non-overlapping structures for exact z-score matching
  block_structures <- c("non_overlapping", "single_block")
  
  for (structure in block_structures) {
    test_data <- generate_block_diagonal_test_data(seed = 123, block_structure = structure)
    
    # Prepare ld_data for partition_LD_matrix
    ld_data <- list(
      LD_matrix = test_data$LD_matrix_full,
      LD_variants = test_data$ref_panel$variant_id,
      block_metadata = test_data$block_metadata
    )
    
    # For non-overlapping structures, use partition_LD_matrix
    partitioned <- partition_LD_matrix(
      ld_data,
      merge_small_blocks = FALSE
    )
    
    # Run RAISS with full matrix
    result_full <- raiss(
      ref_panel = test_data$ref_panel,
      known_zscores = test_data$known_zscores,
      LD_matrix = test_data$LD_matrix_full,
      lamb = 0.01,
      rcond = 0.01,
      R2_threshold = 0.3,
      minimum_ld = 1,
      verbose = FALSE
    )
    
    # Run RAISS with partitioned blocks
    result_blocks <- raiss(
      ref_panel = test_data$ref_panel,
      known_zscores = test_data$known_zscores,
      LD_matrix = partitioned,
      lamb = 0.01,
      rcond = 0.01,
      R2_threshold = 0.3,
      minimum_ld = 1,
      verbose = FALSE
    )
    
    # For non-overlapping blocks, we compare all variants
    result_full_sorted <- result_full$result_nofilter %>% arrange(variant_id)
    result_blocks_sorted <- result_blocks$result_nofilter %>% arrange(variant_id)
    
    # Compare variant IDs
    expect_equal(
      sort(result_full$result_nofilter$variant_id),
      sort(result_blocks$result_nofilter$variant_id),
      info = paste("Variant IDs should match for", structure)
    )
    
    # Compare z-scores with appropriate tolerance
    expect_equal(
      result_full_sorted$z,
      result_blocks_sorted$z,
      tolerance = 0.01,
      info = paste("Z-scores should match for", structure)
    )
    
    # Compare filtered results if present
    if (!is.null(result_full$result_filter) && !is.null(result_blocks$result_filter) &&
        nrow(result_full$result_filter) > 0 && nrow(result_blocks$result_filter) > 0) {
      expect_equal(
        sort(result_full$result_filter$variant_id),
        sort(result_blocks$result_filter$variant_id),
        info = paste("Filtered variant IDs should match for", structure)
      )
      
      result_full_filter_sorted <- result_full$result_filter %>% arrange(variant_id)
      result_blocks_filter_sorted <- result_blocks$result_filter %>% arrange(variant_id)
      
      expect_equal(
        result_full_filter_sorted$z,
        result_blocks_filter_sorted$z,
        tolerance = 0.01,
        info = paste("Filtered Z-scores should match for", structure)
      )
    }
  }
})

test_that("overlapping blocks preserve variant IDs but may have different z-scores", {
  # Test only overlapping structure
  test_data <- generate_block_diagonal_test_data(seed = 123, block_structure = "overlapping")
  
  # Run RAISS with full matrix
  result_full <- raiss(
    ref_panel = test_data$ref_panel,
    known_zscores = test_data$known_zscores,
    LD_matrix = test_data$LD_matrix_full,
    lamb = 0.01,
    rcond = 0.01,
    R2_threshold = 0.3,
    minimum_ld = 1,
    verbose = FALSE
  )
  
  # Run RAISS with block processing
  result_blocks <- raiss(
    ref_panel = test_data$ref_panel,
    known_zscores = test_data$known_zscores,
    LD_matrix = test_data$LD_matrix_blocks,
    lamb = 0.01,
    rcond = 0.01,
    R2_threshold = 0.3,
    minimum_ld = 1,
    verbose = FALSE
  )
  
  # Test 1: Verify all variants are present in both results
  expect_equal(
    sort(result_full$result_nofilter$variant_id),
    sort(result_blocks$result_nofilter$variant_id),
    info = "Both methods should have the same set of variant IDs"
  )
  
  # Test 2: For overlapping blocks, verify boundary variants exist and have valid values
  # Identify boundary variants
  boundary_variants <- character(0)
  for (i in 1:(length(test_data$block_boundaries) - 1)) {
    overlap_pos <- test_data$block_boundaries[[i]][2]
    boundary_variants <- c(boundary_variants, paste0("var", overlap_pos))
  }
  
  # Verify boundary variants exist in results
  expect_true(
    all(boundary_variants %in% result_blocks$result_nofilter$variant_id),
    info = "All boundary variants should be present in block results"
  )
  
  # Verify boundary variants have valid z-scores
  boundary_results <- result_blocks$result_nofilter %>%
    filter(variant_id %in% boundary_variants)
  
  expect_true(
    all(!is.na(boundary_results$z)),
    info = "Boundary variants should have valid z-scores in block results"
  )
  
  # Test 3: Verify non-boundary variants have z-scores with reasonable range
  non_boundary_results <- result_blocks$result_nofilter %>%
    filter(!variant_id %in% boundary_variants)
  
  expect_true(
    all(!is.na(non_boundary_results$z)),
    info = "Non-boundary variants should have valid z-scores"
  )
  
  expect_true(
    all(abs(non_boundary_results$z) < 10),
    info = "Non-boundary variant z-scores should be in reasonable range"
  )
  
  # We deliberately do NOT compare z-score values between full matrix and block processing
  # for overlapping blocks, as differences are expected and valid
})

test_that("raiss handles block boundaries correctly", {
  # Generate test data with overlapping blocks
  test_data <- generate_block_diagonal_test_data(seed = 456, block_structure = "overlapping")
  
  # Define the thresholds explicitly
  test_R2_threshold <- 0.3
  test_minimum_ld <- 1
  
  # Run RAISS with block processing
  result <- raiss(
    ref_panel = test_data$ref_panel,
    known_zscores = test_data$known_zscores,
    LD_matrix = test_data$LD_matrix_blocks,
    lamb = 0.01,
    rcond = 0.01,
    R2_threshold = test_R2_threshold,
    minimum_ld = test_minimum_ld,
    verbose = FALSE
  )
  
  # First verify that the required columns exist in the results
  expect_true(
    "variant_id" %in% names(result$result_nofilter),
    info = "result_nofilter should contain a variant_id column"
  )
  
  expect_true(
    "raiss_R2" %in% names(result$result_nofilter),
    info = "result_nofilter should contain a raiss_R2 column"
  )
  
  expect_true(
    "raiss_ld_score" %in% names(result$result_nofilter),
    info = "result_nofilter should contain a raiss_ld_score column"
  )
  
  # Check that we have only one entry per variant ID (no duplicates)
  expect_equal(
    length(unique(result$result_nofilter$variant_id)),
    length(result$result_nofilter$variant_id),
    info = "Result should have no duplicate variant IDs"
  )
  
  # Check that boundary variants have reasonable values
  boundary_variants <- character(0)
  for (i in 1:(length(test_data$block_boundaries) - 1)) {
    overlap_pos <- test_data$block_boundaries[[i]][2]
    boundary_variants <- c(boundary_variants, paste0("var", overlap_pos))
  }
  
  # Verify that boundary variants exist in the results
  expect_true(
    all(boundary_variants %in% result$result_nofilter$variant_id),
    info = "All boundary variants should be present in the results"
  )
  
  # Get the boundary variant results
  boundary_results <- result$result_nofilter %>%
    filter(variant_id %in% boundary_variants)
  
  # Check R² values for non-NA boundary variants
  non_na_r2 <- boundary_results$raiss_R2[!is.na(boundary_results$raiss_R2)]
  if (length(non_na_r2) > 0) {
    expect_true(
      all(non_na_r2 >= 0 & non_na_r2 <= 1),
      info = "Non-NA boundary variant R² values should be between 0 and 1"
    )
  }
  
  # Check LD scores for non-NA boundary variants
  non_na_ld <- boundary_results$raiss_ld_score[!is.na(boundary_results$raiss_ld_score)]
  if (length(non_na_ld) > 0) {
    expect_true(
      all(non_na_ld >= 0),
      info = "Non-NA boundary variant LD scores should be non-negative"
    )
  }
  
  # Verify that pre-filtering and post-filtering steps handle boundary variants correctly
  if (!is.null(result$result_filter) && nrow(result$result_filter) > 0) {
    # First check if filtered results have the required columns
    expect_true(
      "variant_id" %in% names(result$result_filter),
      info = "result_filter should contain a variant_id column"
    )
    
    expect_true(
      "raiss_R2" %in% names(result$result_filter),
      info = "result_filter should contain a raiss_R2 column"
    )
    
    expect_true(
      "raiss_ld_score" %in% names(result$result_filter),
      info = "result_filter should contain a raiss_ld_score column"
    )
    
    # Check which boundary variants passed the filtering
    boundary_in_filtered <- boundary_variants %in% result$result_filter$variant_id
    
    if (any(boundary_in_filtered)) {
      # Get the filtered boundary variants
      boundary_filtered <- result$result_filter %>%
        filter(variant_id %in% boundary_variants)
      
      # Check that non-NA R² values meet the threshold
      non_na_r2_filtered <- boundary_filtered$raiss_R2[!is.na(boundary_filtered$raiss_R2)]
      if (length(non_na_r2_filtered) > 0) {
        expect_true(
          all(non_na_r2_filtered >= test_R2_threshold),
          info = paste("Non-NA filtered boundary variant R² values should meet the threshold of", test_R2_threshold)
        )
      }
      
      # Check that non-NA LD scores meet the threshold
      non_na_ld_filtered <- boundary_filtered$raiss_ld_score[!is.na(boundary_filtered$raiss_ld_score)]
      if (length(non_na_ld_filtered) > 0) {
        expect_true(
          all(non_na_ld_filtered >= test_minimum_ld),
          info = paste("Non-NA filtered boundary variant LD scores should meet the threshold of", test_minimum_ld)
        )
      }
    }
  }
})

test_that("partition_LD_matrix integrates correctly with RAISS", {
  test_data <- generate_block_diagonal_test_data(seed = 456, block_structure = "non_overlapping")
  
  ld_data <- list(
    LD_matrix = test_data$LD_matrix_full,
    LD_variants = test_data$ref_panel$variant_id,
    block_metadata = test_data$block_metadata
  )
  
  partitioned <- partition_LD_matrix(
    ld_data,
    merge_small_blocks = FALSE
  )
  
  result_full <- raiss(
    ref_panel = test_data$ref_panel,
    known_zscores = test_data$known_zscores,
    LD_matrix = test_data$LD_matrix_full,
    lamb = 0.01,
    rcond = 0.01,
    R2_threshold = 0.3,
    minimum_ld = 1,
    verbose = FALSE
  )
  
  result_partitioned <- raiss(
    ref_panel = test_data$ref_panel,
    known_zscores = test_data$known_zscores,
    LD_matrix = partitioned,
    lamb = 0.01,
    rcond = 0.01,
    R2_threshold = 0.3,
    minimum_ld = 1,
    verbose = FALSE
  )
  
  result_full_sorted <- result_full$result_nofilter %>% arrange(variant_id)
  result_partitioned_sorted <- result_partitioned$result_nofilter %>% arrange(variant_id)
  
  expect_equal(
    result_full_sorted$variant_id,
    result_partitioned_sorted$variant_id,
    info = "Variant IDs should match"
  )
  
  expect_equal(
    result_full_sorted$z,
    result_partitioned_sorted$z,
    tolerance = 1e-4,
    info = "Z-scores should match"
  )
})

# Test 3: Boundary overlap handling
test_that("boundary overlaps are handled correctly", {
  test_data <- generate_block_diagonal_test_data(seed = 789, block_structure = "overlapping")
  
  result_blocks <- raiss(
    ref_panel = test_data$ref_panel,
    known_zscores = test_data$known_zscores,
    LD_matrix = test_data$LD_matrix_blocks,
    lamb = 0.01,
    rcond = 0.01,
    R2_threshold = 0.1,
    minimum_ld = 1,
    verbose = FALSE
  )
  
  variant_counts <- table(test_data$variant_indices$variant_id)
  boundary_vars <- names(variant_counts[variant_counts > 1])
  
  for (var in boundary_vars) {
    expect_equal(
      sum(result_blocks$result_nofilter$variant_id == var),
      1,
      info = paste("Boundary variant", var, "should appear once")
    )
  }
  
  expect_equal(
    nrow(result_blocks$result_nofilter),
    length(unique(result_blocks$result_nofilter$variant_id)),
    info = "No duplicate variants in results"
  )
})

# Test 4: Single-block case
test_that("RAISS handles single-block list correctly", {
  test_data <- generate_block_diagonal_test_data(seed = 202, block_structure = "single_block")
  
  result_full <- raiss(
    ref_panel = test_data$ref_panel,
    known_zscores = test_data$known_zscores,
    LD_matrix = test_data$LD_matrix_full,
    lamb = 0.01,
    rcond = 0.01,
    R2_threshold = 0.3,
    minimum_ld = 1,
    verbose = FALSE
  )
  
  result_single_block <- raiss(
    ref_panel = test_data$ref_panel,
    known_zscores = test_data$known_zscores,
    LD_matrix = test_data$LD_matrix_blocks,
    lamb = 0.01,
    rcond = 0.01,
    R2_threshold = 0.3,
    minimum_ld = 1,
    verbose = FALSE
  )
  
  result_full_sorted <- result_full$result_nofilter %>% arrange(variant_id)
  result_single_block_sorted <- result_single_block$result_nofilter %>% arrange(variant_id)
  
  expect_equal(
    result_full_sorted$z,
    result_single_block_sorted$z,
    tolerance = 1e-6,
    info = "Z-scores should match for single block"
  )
})

# ============================================================================
# Tests for SVD-based genotype matrix path (raiss_single_matrix_from_X)
# ============================================================================

#' Helper: generate a genotype matrix X with corresponding ref_panel,
#' known_zscores, and LD matrix R for equivalence testing.
generate_X_test_data <- function(n = 200, p = 100, n_known = 50, seed = 42) {
  set.seed(seed)
  # Generate genotype-like matrix (dosages 0/1/2)
  X_raw <- matrix(sample(0:2, n * p, replace = TRUE, prob = c(0.25, 0.5, 0.25)),
                  nrow = n, ncol = p)
  # Center and scale
  X <- scale(X_raw)
  X[is.na(X)] <- 0  # zero-variance columns become 0

  # ref_panel for all p variants
  ref_panel <- data.frame(
    chrom = rep(1, p),
    pos = seq(1, p * 10, 10),
    variant_id = paste0("rs", seq_len(p)),
    A1 = sample(c("A", "T", "G", "C"), p, replace = TRUE),
    A2 = sample(c("A", "T", "G", "C"), p, replace = TRUE),
    stringsAsFactors = FALSE
  )
  colnames(X) <- ref_panel$variant_id

  # Select known variants
  known_idx <- sort(sample(seq_len(p), n_known))
  known_zscores <- data.frame(
    chrom = rep(1, n_known),
    pos = ref_panel$pos[known_idx],
    variant_id = ref_panel$variant_id[known_idx],
    A1 = ref_panel$A1[known_idx],
    A2 = ref_panel$A2[known_idx],
    z = rnorm(n_known),
    stringsAsFactors = FALSE
  )

  # LD matrix from X
  R <- cor(X_raw)
  R[is.na(R)] <- 0
  colnames(R) <- rownames(R) <- ref_panel$variant_id

  list(X = X, R = R, ref_panel = ref_panel, known_zscores = known_zscores,
       n = n, p = p, n_known = n_known)
}

test_that("safe_svd basic functionality", {
  set.seed(1)
  mat <- matrix(rnorm(20), nrow = 5, ncol = 4)
  s <- safe_svd(mat)
  expect_equal(length(s$d), min(5, 4))
  expect_true(all(s$d > 0))
  # Reconstruct
  reconstructed <- s$u %*% diag(s$d) %*% t(s$v)
  expect_equal(mat, reconstructed, tolerance = 1e-10)
})

test_that("safe_svd filters small singular values", {
  set.seed(2)
  # Create rank-2 matrix
  u <- matrix(rnorm(10), nrow = 5, ncol = 2)
  v <- matrix(rnorm(8), nrow = 4, ncol = 2)
  mat <- u %*% t(v) + matrix(rnorm(20) * 1e-12, nrow = 5, ncol = 4)
  s <- safe_svd(mat, tol = 1e-6)
  expect_equal(length(s$d), 2)
})

test_that("safe_svd max_rank works", {
  set.seed(3)
  mat <- matrix(rnorm(50), nrow = 10, ncol = 5)
  s <- safe_svd(mat, max_rank = 2)
  expect_equal(length(s$d), 2)
  expect_equal(ncol(s$u), 2)
  expect_equal(ncol(s$v), 2)
})

test_that("safe_svd rejects all-zero matrix", {
  mat <- matrix(0, nrow = 5, ncol = 3)
  expect_error(safe_svd(mat), "all-zero")
})

test_that("X path matches R path: basic equivalence (n > p)", {
  data <- generate_X_test_data(n = 200, p = 100, n_known = 50, seed = 42)

  result_R <- raiss(data$ref_panel, data$known_zscores, LD_matrix = data$R,
                    lamb = 0.01, rcond = 0.01, R2_threshold = 0, minimum_ld = 0,
                    verbose = FALSE)
  result_X <- raiss(data$ref_panel, data$known_zscores, genotype_matrix = data$X,
                    lamb = 0.01, svd_tol = 1e-12, R2_threshold = 0, minimum_ld = 0,
                    verbose = FALSE)

  # Compare imputed z-scores (sort by variant_id for alignment)
  r_sorted <- result_R$result_nofilter %>% arrange(variant_id)
  x_sorted <- result_X$result_nofilter %>% arrange(variant_id)

  expect_equal(r_sorted$z, x_sorted$z, tolerance = 1e-4,
               info = "Imputed z-scores should match between X and R paths")
  expect_equal(r_sorted$Var, x_sorted$Var, tolerance = 1e-4,
               info = "Variance should match between X and R paths")
  expect_equal(r_sorted$raiss_ld_score, x_sorted$raiss_ld_score, tolerance = 1e-4,
               info = "LD scores should match between X and R paths")
})

test_that("X path matches R path: n < p regime", {
  data <- generate_X_test_data(n = 50, p = 200, n_known = 100, seed = 123)

  result_R <- raiss(data$ref_panel, data$known_zscores, LD_matrix = data$R,
                    lamb = 0.01, rcond = 0.01, R2_threshold = 0, minimum_ld = 0,
                    verbose = FALSE)
  result_X <- raiss(data$ref_panel, data$known_zscores, genotype_matrix = data$X,
                    lamb = 0.01, svd_tol = 1e-12, R2_threshold = 0, minimum_ld = 0,
                    verbose = FALSE)

  r_sorted <- result_R$result_nofilter %>% arrange(variant_id)
  x_sorted <- result_X$result_nofilter %>% arrange(variant_id)

  expect_equal(r_sorted$z, x_sorted$z, tolerance = 1e-4,
               info = "z-scores should match in n < p regime")
  expect_equal(r_sorted$Var, x_sorted$Var, tolerance = 1e-4,
               info = "Variance should match in n < p regime")
})

test_that("X path matches R path: n >> p regime", {
  data <- generate_X_test_data(n = 500, p = 50, n_known = 25, seed = 99)

  result_R <- raiss(data$ref_panel, data$known_zscores, LD_matrix = data$R,
                    lamb = 0.01, rcond = 0.01, R2_threshold = 0, minimum_ld = 0,
                    verbose = FALSE)
  result_X <- raiss(data$ref_panel, data$known_zscores, genotype_matrix = data$X,
                    lamb = 0.01, svd_tol = 1e-12, R2_threshold = 0, minimum_ld = 0,
                    verbose = FALSE)

  r_sorted <- result_R$result_nofilter %>% arrange(variant_id)
  x_sorted <- result_X$result_nofilter %>% arrange(variant_id)

  expect_equal(r_sorted$z, x_sorted$z, tolerance = 1e-4)
  expect_equal(r_sorted$Var, x_sorted$Var, tolerance = 1e-4)
})

test_that("X path matches R path: varying lambda", {
  data <- generate_X_test_data(n = 150, p = 80, n_known = 40, seed = 7)

  for (lamb in c(0.001, 0.01, 0.1)) {
    result_R <- raiss(data$ref_panel, data$known_zscores, LD_matrix = data$R,
                      lamb = lamb, rcond = 0.01, R2_threshold = 0, minimum_ld = 0,
                      verbose = FALSE)
    result_X <- raiss(data$ref_panel, data$known_zscores, genotype_matrix = data$X,
                      lamb = lamb, svd_tol = 1e-12, R2_threshold = 0, minimum_ld = 0,
                      verbose = FALSE)

    r_sorted <- result_R$result_nofilter %>% arrange(variant_id)
    x_sorted <- result_X$result_nofilter %>% arrange(variant_id)

    expect_equal(r_sorted$z, x_sorted$z, tolerance = 1e-4,
                 info = paste("z-scores should match for lamb =", lamb))
    expect_equal(r_sorted$Var, x_sorted$Var, tolerance = 1e-4,
                 info = paste("Variance should match for lamb =", lamb))
  }
})

test_that("X path handles all-known edge case", {
  data <- generate_X_test_data(n = 100, p = 50, n_known = 50, seed = 10)
  # Make all variants known
  all_known <- data.frame(
    chrom = data$ref_panel$chrom,
    pos = data$ref_panel$pos,
    variant_id = data$ref_panel$variant_id,
    A1 = data$ref_panel$A1,
    A2 = data$ref_panel$A2,
    z = rnorm(nrow(data$ref_panel)),
    stringsAsFactors = FALSE
  )
  result <- raiss(data$ref_panel, all_known, genotype_matrix = data$X,
                  verbose = FALSE)
  expect_equal(nrow(result$result_nofilter), nrow(data$ref_panel))
})

test_that("X path handles single unknown variant", {
  data <- generate_X_test_data(n = 100, p = 50, n_known = 49, seed = 15)

  result_R <- raiss(data$ref_panel, data$known_zscores, LD_matrix = data$R,
                    lamb = 0.01, rcond = 0.01, R2_threshold = 0, minimum_ld = 0,
                    verbose = FALSE)
  result_X <- raiss(data$ref_panel, data$known_zscores, genotype_matrix = data$X,
                    lamb = 0.01, svd_tol = 1e-12, R2_threshold = 0, minimum_ld = 0,
                    verbose = FALSE)

  r_sorted <- result_R$result_nofilter %>% arrange(variant_id)
  x_sorted <- result_X$result_nofilter %>% arrange(variant_id)

  expect_equal(r_sorted$z, x_sorted$z, tolerance = 1e-4)
})

test_that("X path handles single known variant", {
  data <- generate_X_test_data(n = 100, p = 50, n_known = 1, seed = 20)

  result_X <- raiss(data$ref_panel, data$known_zscores, genotype_matrix = data$X,
                    lamb = 0.01, svd_tol = 1e-12, R2_threshold = 0, minimum_ld = 0,
                    verbose = FALSE)
  expect_true(is.data.frame(result_X$result_nofilter))
  expect_equal(nrow(result_X$result_nofilter), nrow(data$ref_panel))
})

test_that("X path R2 filtering matches R path", {
  data <- generate_X_test_data(n = 200, p = 100, n_known = 50, seed = 42)

  result_R <- raiss(data$ref_panel, data$known_zscores, LD_matrix = data$R,
                    lamb = 0.01, rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5,
                    verbose = FALSE)
  result_X <- raiss(data$ref_panel, data$known_zscores, genotype_matrix = data$X,
                    lamb = 0.01, svd_tol = 1e-12, R2_threshold = 0.6, minimum_ld = 5,
                    verbose = FALSE)

  # Same variants should pass filtering
  r_filtered_ids <- sort(result_R$result_filter$variant_id)
  x_filtered_ids <- sort(result_X$result_filter$variant_id)
  expect_equal(r_filtered_ids, x_filtered_ids,
               info = "Same variants should pass R2/LD filtering")
})

test_that("raiss rejects both LD_matrix and genotype_matrix", {
  data <- generate_X_test_data(n = 50, p = 20, n_known = 10, seed = 1)
  expect_error(
    raiss(data$ref_panel, data$known_zscores,
          LD_matrix = data$R, genotype_matrix = data$X),
    "not both"
  )
})

test_that("raiss rejects neither LD_matrix nor genotype_matrix", {
  data <- generate_X_test_data(n = 50, p = 20, n_known = 10, seed = 1)
  expect_error(
    raiss(data$ref_panel, data$known_zscores),
    "Provide either"
  )
})

test_that("X path with collinear variants matches R path", {
  set.seed(55)
  n <- 150
  p <- 60
  # Create X with some near-duplicate columns
  X_raw <- matrix(sample(0:2, n * p, replace = TRUE, prob = c(0.25, 0.5, 0.25)),
                  nrow = n, ncol = p)
  # Make columns 5 and 6 nearly identical
  X_raw[, 6] <- X_raw[, 5] + sample(c(0, 0, 0, 0, 1), n, replace = TRUE)
  X_raw[X_raw > 2] <- 2

  X <- scale(X_raw)
  X[is.na(X)] <- 0

  ref_panel <- data.frame(
    chrom = rep(1, p), pos = seq(1, p * 10, 10),
    variant_id = paste0("rs", seq_len(p)),
    A1 = rep("A", p), A2 = rep("G", p),
    stringsAsFactors = FALSE
  )
  colnames(X) <- ref_panel$variant_id

  n_known <- 30
  known_idx <- sort(sample(seq_len(p), n_known))
  known_zscores <- data.frame(
    chrom = rep(1, n_known), pos = ref_panel$pos[known_idx],
    variant_id = ref_panel$variant_id[known_idx],
    A1 = ref_panel$A1[known_idx], A2 = ref_panel$A2[known_idx],
    z = rnorm(n_known), stringsAsFactors = FALSE
  )

  R <- cor(X_raw)
  R[is.na(R)] <- 0
  colnames(R) <- rownames(R) <- ref_panel$variant_id

  result_R <- raiss(ref_panel, known_zscores, LD_matrix = R,
                    lamb = 0.01, rcond = 0.01, R2_threshold = 0, minimum_ld = 0,
                    verbose = FALSE)
  result_X <- raiss(ref_panel, known_zscores, genotype_matrix = X,
                    lamb = 0.01, svd_tol = 1e-12, R2_threshold = 0, minimum_ld = 0,
                    verbose = FALSE)

  r_sorted <- result_R$result_nofilter %>% arrange(variant_id)
  x_sorted <- result_X$result_nofilter %>% arrange(variant_id)

  expect_equal(r_sorted$z, x_sorted$z, tolerance = 1e-3,
               info = "Collinear case: z-scores should be close")
})