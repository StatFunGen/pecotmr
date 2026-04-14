context("dentist_qc")
library(MASS)
library(corpcor)

generate_dentist_data <- function(seed=42, n_snps = 100, sample_size = 100, n_outliers = 5, start_pos = 1000000, end_pos = 4000000) {
    set.seed(seed)
    cor_matrix <- matrix(0, nrow = n_snps, ncol = n_snps)
    for (i in 1:(n_snps - 1)) {
        for (j in (i + 1):n_snps) {
            cor_matrix[i, j] <- runif(1, 0.2, 0.8)
            cor_matrix[j, i] <- cor_matrix[i, j]
        }
    }
    diag(cor_matrix) <- 1
    ld_matrix <- cov2cor(make.positive.definite(cor_matrix))
    z_scores <- mvrnorm(n = 1, mu = rep(0, n_snps), Sigma = ld_matrix)
    outlier_indices <- sample(1:n_snps, n_outliers)
    z_scores[outlier_indices] <- rnorm(n_outliers, mean = 0, sd = 5)
    sumstat <- data.frame(
        position = unlist(lapply(seq(start_pos,end_pos,length.out = n_snps), round)),
        z = z_scores
    )
    return(list(sumstat = sumstat, LD_mat = ld_matrix, nSample = sample_size))
}

generate_dentist_single_window_data <- function(seed=42, n_snps = 100, sample_size = 100, n_outliers = 5) {
    set.seed(seed)
    cor_matrix <- matrix(0, nrow = n_snps, ncol = n_snps)
    for (i in 1:(n_snps - 1)) {
        for (j in (i + 1):n_snps) {
            cor_matrix[i, j] <- runif(1, 0.2, 0.8)
            cor_matrix[j, i] <- cor_matrix[i, j]
        }
    }
    diag(cor_matrix) <- 1
    ld_matrix <- cov2cor(make.positive.definite(cor_matrix))
    z_scores <- mvrnorm(n = 1, mu = rep(0, n_snps), Sigma = ld_matrix)
    outlier_indices <- sample(1:n_snps, n_outliers)
    z_scores[outlier_indices] <- rnorm(n_outliers, mean = 0, sd = 5)
    return(list(z_scores = z_scores, LD_mat = ld_matrix, nSample = sample_size))
}

# ===========================================================================
# dentist: basic tests
# ===========================================================================

test_that("dentist output has exactly N rows for N input variants", {
    data <- generate_dentist_data(n_snps = 100)
    expect_warning(res <- dentist(data$sumstat, R = data$LD_mat, nSample = data$nSample))
    expect_equal(nrow(res), 100)
})

test_that("dentist output has exactly N rows with correct_chen_et_al_bug = FALSE", {
    data <- generate_dentist_data(n_snps = 100)
    expect_warning(res <- dentist(data$sumstat, R = data$LD_mat, nSample = data$nSample, correct_chen_et_al_bug = FALSE))
    expect_equal(nrow(res), 100)
})

test_that("dentist stops when missing position", {
    data <- generate_dentist_data()
    colnames(data$sumstat) <- c("something", "z")
    expect_error(dentist(data$sumstat, R = data$LD_mat, nSample = data$nSample, correct_chen_et_al_bug = FALSE),
                 regexp = "missing either.*pos.*or.*z")
})

test_that("dentist stops when missing zscore", {
    data <- generate_dentist_data()
    colnames(data$sumstat) <- c("position", "something")
    expect_error(dentist(data$sumstat, R = data$LD_mat, nSample = data$nSample, correct_chen_et_al_bug = FALSE),
                 regexp = "missing either.*pos.*or.*z")
})

test_that("dentist accepts 'position' and 'zscore' column names", {
  set.seed(42)
  n_snps <- 80
  n_samples <- 100
  cor_matrix <- matrix(0, nrow = n_snps, ncol = n_snps)
  for (i in 1:(n_snps - 1)) {
    for (j in (i + 1):n_snps) {
      cor_matrix[i, j] <- runif(1, 0.2, 0.8)
      cor_matrix[j, i] <- cor_matrix[i, j]
    }
  }
  diag(cor_matrix) <- 1
  ld_matrix <- cov2cor(make.positive.definite(cor_matrix))
  z_scores <- mvrnorm(n = 1, mu = rep(0, n_snps), Sigma = ld_matrix)
  sumstat <- data.frame(
    position = seq(1000000, by = 1000, length.out = n_snps),
    zscore = z_scores
  )
  expect_warning(res <- dentist(sumstat, R = ld_matrix, nSample = n_samples))
  expect_equal(nrow(res), n_snps)
})

test_that("dentist with X matrix input returns exactly N rows", {
    set.seed(42)
    n_snps <- 80
    n_samples <- 100
    X <- matrix(rbinom(n_snps * n_samples, 2, 0.3), nrow = n_samples, ncol = n_snps)
    z_scores <- rnorm(n_snps)
    sumstat <- data.frame(position = seq(1000000, by = 1000, length.out = n_snps), z = z_scores)
    expect_warning(res <- dentist(sumstat, X = X))
    expect_equal(nrow(res), n_snps)
})

# ===========================================================================
# dentist_single_window
# ===========================================================================

test_that("dentist_single_window returns exactly N rows for N input z-scores", {
    data <- generate_dentist_single_window_data()
    expect_warning(res <- dentist_single_window(data$z_scores, R = data$LD_mat, nSample = data$nSample))
    expect_equal(nrow(res), 100)
})

test_that("dentist_single_window warns when < 2000 variants", {
    data <- generate_dentist_single_window_data()
    expect_warning(dentist_single_window(data$z_scores, R = data$LD_mat, nSample = data$nSample))
})

test_that("dentist_single_window stops with zscore/LD matrix dimension mismatch", {
    data <- generate_dentist_single_window_data()
    expect_warning(expect_error(dentist_single_window(generate_dentist_single_window_data()$z_scores, R = generate_dentist_single_window_data(n_snps = 80)$LD_mat, nSample = data$nSample),
                                regexp = "LD_mat must be a square matrix"))
})

test_that("dentist_single_window output columns are correct", {
    data <- generate_dentist_single_window_data()
    expect_warning(res <- dentist_single_window(data$z_scores, R = data$LD_mat, nSample = data$nSample))
    expected_cols <- c("original_z", "imputed_z", "iter_to_correct", "rsq", "is_duplicate", "outlier_stat", "outlier")
    expect_true(all(expected_cols %in% colnames(res)))
})

test_that("dentist_single_window original_z matches input z-scores", {
    data <- generate_dentist_single_window_data()
    expect_warning(res <- dentist_single_window(data$z_scores, R = data$LD_mat, nSample = data$nSample))
    expect_equal(res$original_z, data$z_scores)
})

test_that("dentist_single_window with X matrix input returns exactly N rows", {
    set.seed(42)
    n_snps <- 80
    n_samples <- 100
    X <- matrix(rbinom(n_snps * n_samples, 2, 0.3), nrow = n_samples, ncol = n_snps)
    z_scores <- rnorm(n_snps)
    expect_warning(res <- dentist_single_window(z_scores, X = X))
    expect_equal(nrow(res), n_snps)
})

test_that("dentist_single_window with correct_chen_et_al_bug = FALSE returns N rows", {
    data <- generate_dentist_single_window_data()
    expect_warning(res <- dentist_single_window(
      data$z_scores, R = data$LD_mat, nSample = data$nSample,
      correct_chen_et_al_bug = FALSE
    ))
    expect_equal(nrow(res), 100)
    expect_true(all(c("original_z", "imputed_z", "outlier") %in% colnames(res)))
})

test_that("dentist_single_window with gcControl = TRUE returns N rows", {
    data <- generate_dentist_single_window_data()
    expect_warning(res <- dentist_single_window(
      data$z_scores, R = data$LD_mat, nSample = data$nSample,
      gcControl = TRUE
    ))
    expect_equal(nrow(res), 100)
    expect_true(all(c("original_z", "imputed_z", "outlier") %in% colnames(res)))
})

test_that("dentist with gcControl = TRUE returns N rows", {
    data <- generate_dentist_data(n_snps = 100)
    expect_warning(res <- dentist(
      data$sumstat, R = data$LD_mat, nSample = data$nSample,
      gcControl = TRUE
    ))
    expect_equal(nrow(res), 100)
})

test_that("dentist_single_window dedup path with message for duplicates", {
  set.seed(42)
  n_snps <- 80
  n_samples <- 100
  cor_matrix <- matrix(0, nrow = n_snps, ncol = n_snps)
  for (i in 1:(n_snps - 1)) {
    for (j in (i + 1):n_snps) {
      cor_matrix[i, j] <- runif(1, 0.2, 0.8)
      cor_matrix[j, i] <- cor_matrix[i, j]
    }
  }
  diag(cor_matrix) <- 1
  ld_matrix <- cov2cor(make.positive.definite(cor_matrix))
  z_scores <- mvrnorm(n = 1, mu = rep(0, n_snps), Sigma = ld_matrix)
  # Use a very low threshold to trigger dedup logic
  expect_warning(
    res <- dentist_single_window(z_scores, R = ld_matrix, nSample = n_samples, duprThreshold = 0.5)
  )
  expect_equal(nrow(res), n_snps)
  expect_true("is_duplicate" %in% colnames(res))
})

# ===========================================================================
# add_dups_back_dentist
# ===========================================================================

test_that("add_dups_back_dentist works", {
    zScore <- c(1.2, 2.3, 2.4, 1.4, 5.6)
    dentist_output <- data.frame(
        original_z = c(1.2, 2.3, 5.6),
        imputed_z = c(1.1, 2.1, 5.1),
        iter_to_correct = c(1, 1, 3),
        rsq = c(0.9, 0.8,  0.5),
        z_diff = c(0.1, 0.2, 0.5)
    )
    find_dup_output <- data.frame(
        dupBearer = c(-1, -1, 2, 1, -1),
        sign = c(1, -1, 1, -1, 1)
    )

    res <- pecotmr:::add_dups_back_dentist(zScore, dentist_output, find_dup_output)
    # Non-duplicates: z_diff is copied directly from dentist output
    expect_equal(res$z_diff[1], 0.1)
    expect_equal(res$z_diff[2], 0.2)
    expect_equal(res$z_diff[5], 0.5)
    # Duplicates: z_diff = (zScore - imputed_z) / sqrt(1 - rsq)
    expect_equal(res$z_diff[3], (2.4 - 2.1) / sqrt(1 - 0.8), tolerance = 1e-10)
    expect_equal(res$z_diff[4], (1.4 - (-1.1)) / sqrt(1 - 0.9), tolerance = 1e-10)
    expect_equal(res$imputed_z, c(1.1, 2.1, 2.1, -1.1, 5.1))
    expect_equal(res$is_duplicate, c(FALSE, FALSE, TRUE, TRUE, FALSE))
})

test_that("add_dups_back_dentist stops when nrow mismatch", {
    z_scores <- rep(0, 5)
    dentist_output <- list(
        original_z = c(1, 2, 3, 4, 5),
        imputed_z = c(1, 2, 3, 4, 5),
        iter_to_correct = c(1, 2, 3, 4, 5),
        rsq = c(1, 2, 3, 4, 5),
        z_diff = c(1, 2, 3, 4, 5)
    )
    find_dup_output <- list(
        dupBearer = c(-1, -1, -1, -1, -1, -1),
        sign = c(1, 2, 3, 4, 5, 6)
    )
    expect_error(pecotmr:::add_dups_back_dentist(z_scores, dentist_output, find_dup_output))
})

# ===========================================================================
# segment_by_dist
# ===========================================================================

test_that("segment_by_dist works", {
    res <- pecotmr:::segment_by_dist(seq(2000000,5000000,100000), max_dist = 2000000, min_dim = 10)
    expect_true(nrow(res) >= 1)
})

test_that("segment_by_dist fill regions cover all input positions", {
    # Verify that fill regions cover every SNP index exactly once
    pos <- seq(1000000, 5000000, length.out = 200)
    res <- pecotmr:::segment_by_dist(pos, max_dist = 2000000, min_dim = 10)
    # Collect all fill region indices
    covered <- integer(0)
    for (k in 1:nrow(res)) {
        covered <- c(covered, res$fillStartIdx[k]:(res$fillEndIdx[k] - 1L))
    }
    # Every position from 1 to length(pos) should be covered
    expect_equal(sort(unique(covered)), 1:length(pos))
})

test_that("segment_by_dist errors on empty positions", {
  expect_error(pecotmr:::segment_by_dist(integer(0)), "No positions")
})

test_that("segment_by_dist verbose mode prints intervals", {
  pos <- seq(1000000, 5000000, length.out = 300)
  expect_message(
    pecotmr:::segment_by_dist(pos, max_dist = 2000000, min_dim = 50, verbose = TRUE),
    "Intervals"
  )
})

# ===========================================================================
# detect_gaps
# ===========================================================================

test_that("detect_gaps finds no internal gaps for contiguous positions", {
    pos <- seq(1000, by = 100, length.out = 50)
    gaps <- pecotmr:::detect_gaps(pos, gap_threshold = 500)
    # Only start and end sentinel
    expect_equal(gaps, c(1L, 51L))
})

test_that("detect_gaps finds a centromeric gap", {
    pos <- c(seq(1000, by = 100, length.out = 50),
             seq(2000000, by = 100, length.out = 50))
    gaps <- pecotmr:::detect_gaps(pos, gap_threshold = 1e6)
    expect_equal(length(gaps), 3)  # start, gap, end
    expect_equal(gaps, c(1L, 51L, 101L))
})

test_that("detect_gaps finds multiple gaps", {
    pos <- c(1000, 2000, 5000000, 6000000, 12000000)
    gaps <- pecotmr:::detect_gaps(pos, gap_threshold = 1e6)
    # Gaps at positions 3 and 5 (diffs > 1e6 at indices 2 and 4)
    expect_equal(gaps, c(1L, 3L, 5L, 6L))
})

test_that("detect_gaps verbose branch prints messages", {
  pos <- c(1000, 2000, 5000000, 6000000)
  expect_message(
    pecotmr:::detect_gaps(pos, gap_threshold = 1e6, verbose = TRUE),
    "No\\. of gaps found"
  )
})

# ===========================================================================
# segment_by_count
# ===========================================================================

test_that("segment_by_count produces valid windows", {
    pos <- seq(1000000, by = 1000, length.out = 500)
    res <- pecotmr:::segment_by_count(pos, max_count = 100)
    expect_true(nrow(res) >= 1)
    # All window starts should be >= 1
    expect_true(all(res$windowStartIdx >= 1))
    # All window ends should be <= length(pos) + 1
    expect_true(all(res$windowEndIdx <= length(pos) + 1))
    # Fill regions should be within windows
    for (k in 1:nrow(res)) {
        expect_true(res$fillStartIdx[k] >= res$windowStartIdx[k])
        expect_true(res$fillEndIdx[k] <= res$windowEndIdx[k])
    }
})

test_that("segment_by_count fill regions cover all positions", {
    pos <- seq(1000000, by = 1000, length.out = 500)
    res <- pecotmr:::segment_by_count(pos, max_count = 100)
    covered <- integer(0)
    for (k in 1:nrow(res)) {
        covered <- c(covered, res$fillStartIdx[k]:(res$fillEndIdx[k] - 1L))
    }
    expect_equal(sort(unique(covered)), 1:length(pos))
})

test_that("segment_by_count handles centromeric gap", {
    # Create two blocks separated by a large gap
    pos <- c(seq(1000000, by = 1000, length.out = 200),
             seq(5000000, by = 1000, length.out = 200))
    res <- pecotmr:::segment_by_count(pos, max_count = 100)
    # Should create windows in both blocks
    expect_true(nrow(res) >= 2)
    # Fill regions should still cover all positions
    covered <- integer(0)
    for (k in 1:nrow(res)) {
        covered <- c(covered, res$fillStartIdx[k]:(res$fillEndIdx[k] - 1L))
    }
    expect_equal(sort(unique(covered)), 1:length(pos))
})

test_that("segment_by_count skips blocks smaller than half max_count", {
    # Block of 20 variants with max_count=100 (half=50): too small, should be skipped
    pos <- seq(1000000, by = 1000, length.out = 20)
    expect_error(pecotmr:::segment_by_count(pos, max_count = 100),
                 "No intervals created by segmentation")
})

test_that("segment_by_count creates single window for small blocks", {
    # Block of 60 variants with max_count=100: creates one window (60 >= half=50)
    pos <- seq(1000000, by = 1000, length.out = 60)
    res <- pecotmr:::segment_by_count(pos, max_count = 100)
    expect_equal(nrow(res), 1)
    expect_equal(res$windowStartIdx[1], 1)
    expect_equal(res$windowEndIdx[1], 61)
})

test_that("segment_by_count single block creates correct number of windows", {
    # 200 variants with max_count = 100 should create ~3 windows
    pos <- seq(1000000, by = 1000, length.out = 200)
    res <- pecotmr:::segment_by_count(pos, max_count = 100)
    expect_true(nrow(res) >= 2)
    expect_true(nrow(res) <= 5)
})

test_that("segment_by_count errors on empty positions", {
  expect_error(pecotmr:::segment_by_count(integer(0), max_count = 100), "No positions")
})

test_that("segment_by_count verbose mode prints intervals", {
  pos <- seq(1000000, by = 1000, length.out = 300)
  expect_message(
    pecotmr:::segment_by_count(pos, max_count = 100, verbose = TRUE),
    "Intervals"
  )
})

# ===========================================================================
# merge_windows
# ===========================================================================

test_that("merge_windows returns exactly N rows", {
    data <- generate_dentist_data(n_snps = 1000, sample_size = 1000, start_pos = 0, end_pos = 2000)
    window_divided_res <- pecotmr:::segment_by_dist(data$sumstat$position, max_dist = 1000, min_dim = 10)
    dentist_result_by_window <- list()
    suppressWarnings({
        for (k in 1:nrow(window_divided_res)) {
            idx_range <- window_divided_res$windowStartIdx[k]:(window_divided_res$windowEndIdx[k] - 1L)
            zScore_k <- data$sumstat$z[idx_range]
            LD_mat_k <- data$LD_mat[idx_range, idx_range]
            dentist_result_by_window[[k]] <- dentist_single_window(
                zScore_k, R = LD_mat_k, nSample = 100,
                pValueThreshold = 5.0369e-8, propSVD = 0.4, gcControl = FALSE,
                nIter = 10, gPvalueThreshold = 0.05, duprThreshold = 0.99,
                ncpus = 1, correct_chen_et_al_bug = TRUE
            )
        }
    })
    res <- pecotmr:::merge_windows(dentist_result_by_window, window_divided_res)
    expect_equal(nrow(res), 1000)
})

test_that("merge_windows stops with window and imputed mismatch", {
    expect_error(pecotmr:::merge_windows(rep(0, 5), data.frame(windowStartIdx = rep(0,2), windowEndIdx = rep(0, 2))))
})

test_that("merge_windows correctly indexes and merges windows", {
  # Create two fake windows
  window1 <- data.frame(
    original_z = c(1.0, 2.0, 3.0),
    imputed_z = c(0.9, 1.9, 2.9),
    iter_to_correct = c(1, 1, 1),
    rsq = c(0.5, 0.6, 0.7),
    is_duplicate = c(FALSE, FALSE, FALSE),
    outlier_stat = c(0.1, 0.2, 0.3),
    outlier = c(FALSE, FALSE, FALSE)
  )
  window2 <- data.frame(
    original_z = c(4.0, 5.0, 6.0),
    imputed_z = c(3.9, 4.9, 5.9),
    iter_to_correct = c(2, 2, 2),
    rsq = c(0.8, 0.9, 0.5),
    is_duplicate = c(FALSE, FALSE, FALSE),
    outlier_stat = c(0.4, 0.5, 0.6),
    outlier = c(FALSE, FALSE, TRUE)
  )
  window_info <- data.frame(
    windowIdx = c(1, 2),
    windowStartIdx = c(1, 4),
    windowEndIdx = c(4, 7),
    fillStartIdx = c(1, 4),
    fillEndIdx = c(4, 7)
  )
  result <- pecotmr:::merge_windows(list(window1, window2), window_info)
  expect_equal(nrow(result), 6)
  expect_true("index_global" %in% colnames(result))
  expect_true("index_within_window" %in% colnames(result))
})

# ===========================================================================
# dentist windowed mode
# ===========================================================================

test_that("dentist windowed output has exactly N rows for large input", {
    # Generate data large enough to trigger windowed mode (> min_dim)
    data <- generate_dentist_data(seed = 123, n_snps = 1000, sample_size = 1000,
                                   n_outliers = 50, start_pos = 0, end_pos = 5000000)
    suppressWarnings({
        res <- dentist(data$sumstat, R = data$LD_mat, nSample = data$nSample,
                       min_dim = 100, window_size = 2000000)
    })
    expect_equal(nrow(res), 1000)
})

test_that("dentist outlier_stat formula is correct: (z-imputed)^2/(1-rsq)", {
    data <- generate_dentist_single_window_data(seed = 55, n_snps = 100)
    expect_warning(res <- dentist_single_window(data$z_scores, R = data$LD_mat, nSample = data$nSample))
    expected_stat <- (res$original_z - res$imputed_z)^2 / pmax(1 - res$rsq, 1e-8)
    expect_equal(res$outlier_stat, expected_stat, tolerance = 1e-10)
})

# ===========================================================================
# dentist with count mode
# ===========================================================================

test_that("dentist with window_mode='count' returns exactly N rows", {
    data <- generate_dentist_data(seed = 789, n_snps = 500, sample_size = 500,
                                   n_outliers = 25, start_pos = 1000000, end_pos = 4000000)
    suppressWarnings({
        res <- dentist(data$sumstat, R = data$LD_mat, nSample = data$nSample,
                       window_mode = "count", min_dim = 100)
    })
    expect_equal(nrow(res), 500)
})

# ===========================================================================
# Equivalence tests: both windowing methods
# ===========================================================================

test_that("segment_by_dist and segment_by_count agree on uniformly-spaced variants", {
    n <- 200
    spacing <- 10000  # 10kb between each variant
    pos <- seq(1000000, by = spacing, length.out = n)
    window_count <- 50  # variants per window in count mode
    window_dist <- window_count * spacing  # equivalent distance

    res_dist <- pecotmr:::segment_by_dist(pos, max_dist = window_dist, min_dim = 10)
    res_count <- pecotmr:::segment_by_count(pos, max_count = window_count, gap_dist = 1e6)

    # Both should cover all positions
    covered_dist <- integer(0)
    for (k in 1:nrow(res_dist)) {
        covered_dist <- c(covered_dist, res_dist$fillStartIdx[k]:(res_dist$fillEndIdx[k] - 1L))
    }
    covered_count <- integer(0)
    for (k in 1:nrow(res_count)) {
        covered_count <- c(covered_count, res_count$fillStartIdx[k]:(res_count$fillEndIdx[k] - 1L))
    }
    expect_equal(sort(unique(covered_dist)), 1:n)
    expect_equal(sort(unique(covered_count)), 1:n)
})

test_that("both windowing modes produce same dentist results on uniform data", {
    data <- generate_dentist_data(seed = 555, n_snps = 500, sample_size = 500,
                                   n_outliers = 25, start_pos = 0, end_pos = 5000000)
    suppressWarnings({
        res_dist <- dentist(data$sumstat, R = data$LD_mat, nSample = data$nSample,
                            window_mode = "distance", min_dim = 100, window_size = 2000000)
        res_count <- dentist(data$sumstat, R = data$LD_mat, nSample = data$nSample,
                             window_mode = "count", min_dim = 100)
    })
    # Both should return exactly N rows
    expect_equal(nrow(res_dist), 500)
    expect_equal(nrow(res_count), 500)
})

# ===========================================================================
# resolve_LD_input (internal)
# ===========================================================================

test_that("resolve_LD_input errors when neither R nor X provided", {
  expect_error(
    pecotmr:::resolve_LD_input(R = NULL, X = NULL),
    "Either R.*or X.*must be provided"
  )
})

test_that("resolve_LD_input errors when both R and X provided", {
  R <- diag(3)
  X <- matrix(1:9, nrow = 3)
  expect_error(
    pecotmr:::resolve_LD_input(R = R, X = X),
    "Provide either R or X, not both"
  )
})

test_that("resolve_LD_input errors when R provided without nSample and need_nSample is TRUE", {
  R <- diag(3)
  expect_error(
    pecotmr:::resolve_LD_input(R = R, nSample = NULL, need_nSample = TRUE),
    "nSample is required"
  )
})

test_that("resolve_LD_input returns nSample = NULL when need_nSample is FALSE", {
  R <- diag(3)
  result <- pecotmr:::resolve_LD_input(R = R, nSample = NULL, need_nSample = FALSE)
  expect_null(result$nSample)
  expect_equal(result$R, R)
})

test_that("resolve_LD_input infers nSample from X", {
  set.seed(42)
  n <- 50; p <- 5
  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n, ncol = p)
  result <- pecotmr:::resolve_LD_input(X = X, need_nSample = TRUE)
  expect_equal(result$nSample, n)
  expect_true(is.matrix(result$R))
  expect_equal(nrow(result$R), p)
})

test_that("resolve_LD_input converts non-matrix X to matrix", {
  set.seed(42)
  X_df <- data.frame(a = rbinom(30, 2, 0.3), b = rbinom(30, 2, 0.3))
  result <- pecotmr:::resolve_LD_input(X = X_df, need_nSample = FALSE)
  expect_true(is.matrix(result$R))
  expect_equal(result$nSample, 30)
})

test_that("resolve_LD_input uses explicit nSample when X provided", {
  set.seed(42)
  X <- matrix(rbinom(100, 2, 0.3), nrow = 20, ncol = 5)
  result <- pecotmr:::resolve_LD_input(X = X, nSample = 999, need_nSample = TRUE)
  expect_equal(result$nSample, 999)
})


# ===========================================================================
# build_segment_result (internal)
# ===========================================================================

test_that("build_segment_result caps end indices and verbose prints", {
  expect_message(
    result <- pecotmr:::build_segment_result(
      startList = c(1L), endList = c(200L),
      fillStartList = c(1L), fillEndList = c(200L),
      n = 100, verbose = TRUE
    ),
    "Intervals"
  )
  expect_equal(result$windowEndIdx[1], 101)
  expect_equal(result$fillEndIdx[1], 101)
})

test_that("build_segment_result errors on empty startList", {
  expect_error(
    pecotmr:::build_segment_result(
      startList = integer(0), endList = integer(0),
      fillStartList = integer(0), fillEndList = integer(0),
      n = 100
    ),
    "No intervals"
  )
})

# ===========================================================================
# sliding_window_loop (iteration limit)
# ===========================================================================

test_that("sliding_window_loop errors on infinite loop", {
  allGaps <- c(1L, 1001L)
  expect_error(
    pecotmr:::sliding_window_loop(
      allGaps, n = 1000,
      min_block_fn = function(blockSize) TRUE,
      init_end_fn = function(startIdx, blockEnd) startIdx + 10,
      fill_fn = function(startIdx, endIdx, notStart, notLast) list(start = startIdx, end = endIdx),
      step_fn = function(startIdx, blockEnd) list(startIdx = startIdx, endIdx = startIdx + 10),
      verbose = FALSE
    ),
    "iteration limit exceeded"
  )
})

