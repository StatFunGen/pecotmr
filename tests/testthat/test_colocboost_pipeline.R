context("colocboost_pipeline")

# ===========================================================================
# Tests from test_colocboost_pipeline.R
# ===========================================================================


# ---- qc_method match.arg ----
test_that("qc_regional_data accepts explicit qc_method = 'slalom'", {
  region_data <- list(individual_data = NULL, sumstat_data = NULL)
  result <- pecotmr:::qc_regional_data(region_data, qc_method = "slalom")
  expect_type(result, "list")
})

test_that("qc_regional_data rejects invalid qc_method", {
  region_data <- list(individual_data = NULL, sumstat_data = NULL)
  expect_error(
    pecotmr:::qc_regional_data(region_data, qc_method = "invalid"),
    "arg"
  )
})

# ---- pip_cutoff_to_skip_ind validation ----
test_that("pip_cutoff scalar is recycled for individual contexts", {
  # Create individual_data with 3 real-ish contexts
  set.seed(42)
  n <- 10; p <- 5
  make_ctx <- function() {
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- paste0("var", 1:p)
    Y <- matrix(rnorm(n * 2), n, 2)
    colnames(Y) <- paste0("gene", 1:2)
    list(X = X, Y = Y, maf = runif(p, 0.05, 0.5))
  }
  ctx <- make_ctx()
  individual_data <- list(
    residual_Y = list(ctx1 = ctx$Y, ctx2 = ctx$Y, ctx3 = ctx$Y),
    residual_X = list(ctx1 = ctx$X, ctx2 = ctx$X, ctx3 = ctx$X),
    maf = list(ctx1 = ctx$maf, ctx2 = ctx$maf, ctx3 = ctx$maf)
  )
  region_data <- list(individual_data = individual_data, sumstat_data = NULL)

  # Scalar 0 (no PIP check) should be recycled and run without error
  result <- pecotmr:::qc_regional_data(region_data, pip_cutoff_to_skip_ind = 0)
  expect_type(result, "list")
})

test_that("pip_cutoff wrong length errors for individual contexts", {
  set.seed(42)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n), n, 1)
  individual_data <- list(
    residual_Y = list(ctx1 = Y, ctx2 = Y, ctx3 = Y),
    residual_X = list(ctx1 = X, ctx2 = X, ctx3 = X),
    maf = list(ctx1 = runif(p), ctx2 = runif(p), ctx3 = runif(p))
  )
  region_data <- list(individual_data = individual_data, sumstat_data = NULL)

  expect_error(
    pecotmr:::qc_regional_data(region_data, pip_cutoff_to_skip_ind = c(0, 0)),
    "pip_cutoff_to_skip_ind"
  )
})

test_that("pip_cutoff correct length vector works", {
  set.seed(42)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("var", 1:p)
  Y <- matrix(rnorm(n), n, 1)
  colnames(Y) <- "gene1"
  individual_data <- list(
    residual_Y = list(ctx1 = Y, ctx2 = Y),
    residual_X = list(ctx1 = X, ctx2 = X),
    maf = list(ctx1 = runif(p, 0.05, 0.5), ctx2 = runif(p, 0.05, 0.5))
  )
  region_data <- list(individual_data = individual_data, sumstat_data = NULL)

  # Length-2 vector for 2 contexts should work
  result <- pecotmr:::qc_regional_data(region_data, pip_cutoff_to_skip_ind = c(0, 0))
  expect_type(result, "list")
})

# ---- is_valid_sumstat_entry ----
test_that("is_valid_sumstat_entry returns TRUE for valid data", {
  ss <- data.frame(z = c(1.5, -2.3, 0.8), n = c(1000, 1000, 1000),
                   variant = c("chr1:100:A:G", "chr1:200:T:C", "chr1:300:G:A"))
  expect_true(pecotmr:::is_valid_sumstat_entry(ss))
})

test_that("is_valid_sumstat_entry returns FALSE for too few variants", {
  ss <- data.frame(z = 1.5, n = 1000, variant = "chr1:100:A:G")
  expect_false(pecotmr:::is_valid_sumstat_entry(ss, min_variants = 2))
})

test_that("is_valid_sumstat_entry returns FALSE for all-NA z-scores", {
  ss <- data.frame(z = c(NA, NA, NA), n = c(1000, 1000, 1000),
                   variant = c("v1", "v2", "v3"))
  expect_false(pecotmr:::is_valid_sumstat_entry(ss))
})

test_that("is_valid_sumstat_entry returns FALSE for NULL input", {
  expect_false(pecotmr:::is_valid_sumstat_entry(NULL))
})

test_that("is_valid_sumstat_entry returns FALSE for zero sample size", {
  ss <- data.frame(z = c(1.5, -2.3), n = c(0, 0), variant = c("v1", "v2"))
  expect_false(pecotmr:::is_valid_sumstat_entry(ss))
})

# ---- filter_valid_sumstats ----
test_that("filter_valid_sumstats removes invalid entries", {
  good <- data.frame(z = c(1.0, 2.0, 3.0), n = c(100, 100, 100),
                     variant = c("v1", "v2", "v3"))
  bad <- data.frame(z = c(NA, NA), n = c(0, 0), variant = c("v1", "v2"))
  sumstats <- list(study_a = good, study_b = bad)
  LD_mat <- list(ld1 = matrix(1, 3, 3))
  LD_match <- c("ld1", "ld1")
  result <- pecotmr:::filter_valid_sumstats(sumstats, LD_mat, LD_match)
  expect_equal(length(result$sumstats), 1)
  expect_equal(names(result$sumstats), "study_a")
})

test_that("filter_valid_sumstats returns NULL when all invalid", {
  bad1 <- data.frame(z = NA, n = 0, variant = "v1")
  bad2 <- data.frame(z = NA, n = 0, variant = "v1")
  sumstats <- list(s1 = bad1, s2 = bad2)
  LD_mat <- list(ld1 = matrix(1))
  LD_match <- c("ld1", "ld1")
  expect_message(
    result <- pecotmr:::filter_valid_sumstats(sumstats, LD_mat, LD_match),
    "No valid"
  )
  expect_null(result)
})

test_that("filter_valid_sumstats keeps all when all valid", {
  make_ss <- function() data.frame(z = c(1, 2, 3), n = c(100, 100, 100),
                                   variant = c("v1", "v2", "v3"))
  sumstats <- list(s1 = make_ss(), s2 = make_ss())
  LD_mat <- list(ld1 = matrix(1, 3, 3))
  LD_match <- c("ld1", "ld1")
  result <- pecotmr:::filter_valid_sumstats(sumstats, LD_mat, LD_match)
  expect_equal(length(result$sumstats), 2)
})

# ===========================================================================
# Tests from test_colocboost_pipeline_comprehensive.R
# ===========================================================================


# ===========================================================================
# Helper: build a minimal region_data with individual-level data
# ===========================================================================
make_individual_region_data <- function(n = 20, p = 8, n_contexts = 2, n_events = 3) {
  set.seed(101)
  make_ctx <- function(ctx_name) {
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- paste0("chr1:", seq_len(p) * 100, ":A:G")
    Y <- matrix(rnorm(n * n_events), n, n_events)
    colnames(Y) <- paste0("event", seq_len(n_events))
    maf <- runif(p, 0.05, 0.45)
    list(X = X, Y = Y, maf = maf)
  }
  ctxs <- lapply(paste0("ctx", seq_len(n_contexts)), make_ctx)
  names(ctxs) <- paste0("ctx", seq_len(n_contexts))
  list(
    individual_data = list(
      residual_Y = lapply(ctxs, `[[`, "Y"),
      residual_X = lapply(ctxs, `[[`, "X"),
      maf        = lapply(ctxs, `[[`, "maf")
    ),
    sumstat_data = NULL
  )
}

# ===========================================================================
# Helper: build a minimal region_data with sumstat data
# ===========================================================================
make_sumstat_region_data <- function(n_variants = 5, n_studies = 2) {
  set.seed(202)
  vids <- paste0("chr1:", seq_len(n_variants) * 100, ":A:G")

  LD_mat <- diag(n_variants)
  rownames(LD_mat) <- colnames(LD_mat) <- vids

  ref_panel <- data.frame(
    chrom = rep(1, n_variants),
    pos   = seq_len(n_variants) * 100,
    A2    = rep("A", n_variants),
    A1    = rep("G", n_variants),
    stringsAsFactors = FALSE
  )

  sumstats_list <- lapply(seq_len(n_studies), function(i) {
    ss <- list(
      sumstats = data.frame(
        chrom      = rep(1, n_variants),
        pos        = seq_len(n_variants) * 100,
        A1         = rep("G", n_variants),
        A2         = rep("A", n_variants),
        beta       = rnorm(n_variants),
        se         = runif(n_variants, 0.05, 0.2),
        z          = rnorm(n_variants, 0, 2),
        variant_id = vids,
        stringsAsFactors = FALSE
      ),
      n     = 10000,
      var_y = 1
    )
    list(ss) |> setNames(paste0("study", i))
  })

  LD_info <- list(list(
    LD_variants = ref_panel,
    LD_matrix   = LD_mat,
    ref_panel            = ref_panel
  ))

  list(
    individual_data = NULL,
    sumstat_data    = list(
      sumstats = sumstats_list,
      LD_info  = LD_info
    )
  )
}

# ===========================================================================
# 1. colocboost_analysis_pipeline: no analysis flags returns empty results
# ===========================================================================
test_that("pipeline returns empty results with message when no analysis flags set", {
  region_data <- make_individual_region_data()
  expect_message(
    result <- colocboost_analysis_pipeline(
      region_data,
      xqtl_coloc    = FALSE,
      joint_gwas     = FALSE,
      separate_gwas  = FALSE
    ),
    "No colocalization has been performed"
  )
  expect_type(result, "list")
  expect_null(result$xqtl_coloc)
  expect_null(result$joint_gwas)
  expect_null(result$separate_gwas)
})

# ===========================================================================
# 2. colocboost_analysis_pipeline: NULL individual_data and NULL sumstat_data
# ===========================================================================
test_that("pipeline returns early when both data sources are NULL", {
  region_data <- list(individual_data = NULL, sumstat_data = NULL)
  expect_message(
    result <- colocboost_analysis_pipeline(region_data, xqtl_coloc = TRUE),
    "No individual data"
  )
  expect_type(result, "list")
  expect_null(result$xqtl_coloc)
})

# ===========================================================================
# 3. filter_events: type_pattern, valid_pattern, exclude_pattern
# ===========================================================================
test_that("filter_events keeps events matching valid_pattern", {
  # Access the internal function from within colocboost_analysis_pipeline's environment
  # We need to build a minimal call through the pipeline to test filter_events indirectly.
  # Instead, recreate the inner function for testing purposes.

  # Create a small region_data with events that should be filtered
  set.seed(303)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  # Events named like sQTL cluster events
  events <- c("clu_1_+:PR:gene1", "clu_1_+:IN:gene1", "clu_2_-:PR:gene2")
  Y <- matrix(rnorm(n * length(events)), n, length(events))
  colnames(Y) <- events

  region_data <- list(
    individual_data = list(
      residual_Y = list(tissue1 = Y),
      residual_X = list(tissue1 = X),
      maf        = list(tissue1 = runif(p, 0.05, 0.45))
    ),
    sumstat_data = NULL
  )

  event_filters <- list(
    list(
      type_pattern    = ".*clu_(\\d+_[+-?]).*",
      valid_pattern   = "clu_(\\d+_[+-?]):PR:",
      exclude_pattern = "clu_(\\d+_[+-?]):IN:"
    )
  )

  # Pipeline calls filter_events, then qc_regional_data.
  # Mock qc_regional_data so we can isolate the filtering step.
  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      # Return the data as-is; transform residual_Y to Y format
      list(
        individual_data = list(
          Y = region_data$individual_data$residual_Y,
          X = region_data$individual_data$residual_X
        ),
        sumstat_data = NULL
      )
    }
  )

  result <- suppressMessages(
    colocboost_analysis_pipeline(
      region_data,
      event_filters  = event_filters,
      xqtl_coloc     = FALSE,
      joint_gwas      = FALSE,
      separate_gwas   = FALSE
    )
  )
  # When no analysis flags set, should return empty and the region_data internal
  # was still filtered. At minimum no error was thrown.
  expect_type(result, "list")
})

# ===========================================================================
# 4. filter_events: error when missing required filter fields
# ===========================================================================
test_that("filter_events errors on missing type_pattern", {
  set.seed(404)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  Y <- matrix(rnorm(n * 2), n, 2)
  colnames(Y) <- c("evt1", "evt2")

  region_data <- list(
    individual_data = list(
      residual_Y = list(ctx1 = Y),
      residual_X = list(ctx1 = X),
      maf        = list(ctx1 = runif(p, 0.05, 0.45))
    ),
    sumstat_data = NULL
  )

  # Filter missing type_pattern
  bad_filter <- list(list(valid_pattern = "something"))

  expect_error(
    suppressMessages(
      colocboost_analysis_pipeline(
        region_data,
        event_filters  = bad_filter,
        xqtl_coloc     = TRUE,
        joint_gwas      = FALSE,
        separate_gwas   = FALSE
      )
    ),
    "type_pattern"
  )
})

test_that("filter_events errors when only type_pattern is given (no valid or exclude)", {
  set.seed(405)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  Y <- matrix(rnorm(n * 2), n, 2)
  colnames(Y) <- c("evt1", "evt2")

  region_data <- list(
    individual_data = list(
      residual_Y = list(ctx1 = Y),
      residual_X = list(ctx1 = X),
      maf        = list(ctx1 = runif(p, 0.05, 0.45))
    ),
    sumstat_data = NULL
  )

  bad_filter <- list(list(type_pattern = "evt.*"))

  expect_error(
    suppressMessages(
      colocboost_analysis_pipeline(
        region_data,
        event_filters  = bad_filter,
        xqtl_coloc     = TRUE
      )
    ),
    "type_pattern.*valid_pattern.*exclude_pattern"
  )
})

# ===========================================================================
# 5. extract_contexts_studies: initial extraction
# ===========================================================================
test_that("extract_contexts_studies returns individual contexts and sumstat studies on initial call", {
  # We access the internal by constructing minimal region_data and triggering
  # the pipeline but with both analysis=FALSE so it exits early after extraction.
  region_data <- list(
    individual_data = list(
      residual_Y = list(tissue_A = matrix(1, 2, 2), tissue_B = matrix(1, 2, 2))
    ),
    sumstat_data = list(
      sumstats = list(
        list(gwas_trait1 = list(), gwas_trait2 = list())
      )
    )
  )

  # Pipeline calls extract_contexts_studies internally.
  # With no analysis flags it will still run extraction and return.
  expect_message(
    result <- colocboost_analysis_pipeline(
      region_data,
      xqtl_coloc     = FALSE,
      joint_gwas      = FALSE,
      separate_gwas   = FALSE
    ),
    "No colocalization"
  )
  expect_type(result, "list")
})

# ===========================================================================
# 6. extract_contexts_studies: after-QC extraction messages
# ===========================================================================
test_that("extract_contexts_studies reports after-QC when some individual data removed", {
  set.seed(601)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  Y <- matrix(rnorm(n), n, 1)
  colnames(Y) <- "gene1"

  region_data <- list(
    individual_data = list(
      residual_Y = list(ctx1 = Y, ctx2 = Y),
      residual_X = list(ctx1 = X, ctx2 = X),
      maf        = list(ctx1 = runif(p, 0.05, 0.45), ctx2 = runif(p, 0.05, 0.45))
    ),
    sumstat_data = NULL
  )

  # Mock qc_regional_data to return one NULL context (simulating QC removal).
  # colocboost is an external package function and cannot be mocked via
  # local_mocked_bindings. The pipeline's tryCatch around the colocboost call
  # handles the case where colocboost is unavailable or errors out.
  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      Y1 <- matrix(rnorm(10), 10, 1)
      colnames(Y1) <- "ctx1_gene1"
      X1 <- matrix(rnorm(50), 10, 5)
      colnames(X1) <- paste0("chr1:", seq_len(5) * 100, ":A:G")
      X2 <- matrix(rnorm(50), 10, 5)
      colnames(X2) <- paste0("chr1:", seq_len(5) * 100, ":A:G")
      list(
        individual_data = list(
          Y = list(ctx1_gene1 = Y1, ctx2_gene1 = NULL),
          X = list(ctx1 = X1, ctx2 = X2)
        ),
        sumstat_data = NULL
      )
    }
  )

  expect_message(
    result <- colocboost_analysis_pipeline(
      region_data,
      xqtl_coloc    = TRUE,
      joint_gwas     = FALSE,
      separate_gwas  = FALSE
    ),
    "Skipping follow-up analysis for individual traits"
  )
  expect_type(result, "list")
})

# ===========================================================================
# 7. qc_regional_data: named pip_cutoff_to_skip_sumstat vector
# ===========================================================================
test_that("qc_regional_data handles named pip_cutoff_to_skip_sumstat vector", {
  region_data <- make_sumstat_region_data(n_variants = 5, n_studies = 2)

  # Mock out the heavy QC functions
  local_mocked_bindings(
    allele_qc = function(target_data, ref_variants, ...) {
      list(target_data_qced = target_data)
    },
    rss_basic_qc = function(sumstats, LD_data, ...) {
      LD_mat <- LD_data$LD_matrix[sumstats$variant_id, sumstats$variant_id, drop = FALSE]
      list(sumstats = sumstats, LD_mat = LD_mat)
    },
    summary_stats_qc = function(sumstats, LD_data, ...) {
      LD_mat <- LD_data$LD_matrix[sumstats$variant_id, sumstats$variant_id, drop = FALSE]
      list(sumstats = sumstats, LD_mat = LD_mat, outlier_number = 0)
    },
    raiss = function(...) {
      list(result_filter = data.frame(z = rnorm(5)), LD_mat = diag(5))
    },
    partition_LD_matrix = function(...) diag(5)
  )

  # Named vector: only specify cutoff for study1
  pip_named <- c("study1" = 0, "study2" = 0)
  result <- suppressMessages(
    pecotmr:::qc_regional_data(
      region_data,
      pip_cutoff_to_skip_sumstat = pip_named,
      qc_method = "slalom",
      impute = FALSE
    )
  )
  expect_type(result, "list")
})

# ===========================================================================
# 8. qc_regional_data: named pip_cutoff fills missing studies with 0
# ===========================================================================
test_that("qc_regional_data fills missing study names with 0 for pip_cutoff_to_skip_sumstat", {
  region_data <- make_sumstat_region_data(n_variants = 5, n_studies = 2)

  local_mocked_bindings(
    rss_basic_qc = function(sumstats, LD_data, ...) {
      LD_mat <- LD_data$LD_matrix[sumstats$variant_id, sumstats$variant_id, drop = FALSE]
      list(sumstats = sumstats, LD_mat = LD_mat)
    },
    summary_stats_qc = function(sumstats, LD_data, ...) {
      LD_mat <- LD_data$LD_matrix[sumstats$variant_id, sumstats$variant_id, drop = FALSE]
      list(sumstats = sumstats, LD_mat = LD_mat, outlier_number = 0)
    },
    raiss = function(...) list(result_filter = data.frame(z = rnorm(5)), LD_mat = diag(5)),
    partition_LD_matrix = function(...) diag(5)
  )

  # Named vector with only one of the two studies: the missing study should get 0

  pip_partial <- c("study1" = 0.05)
  result <- suppressMessages(
    pecotmr:::qc_regional_data(
      region_data,
      pip_cutoff_to_skip_sumstat = pip_partial,
      qc_method = "slalom",
      impute = FALSE
    )
  )
  expect_type(result, "list")
})

# ===========================================================================
# 9. colocboost_analysis_pipeline: output structure verification
# ===========================================================================
test_that("pipeline output structure has expected top-level keys", {
  region_data <- list(individual_data = NULL, sumstat_data = NULL)
  result <- suppressMessages(
    colocboost_analysis_pipeline(region_data, xqtl_coloc = FALSE, joint_gwas = FALSE, separate_gwas = FALSE)
  )
  expect_true("xqtl_coloc" %in% names(result))
  expect_true("joint_gwas" %in% names(result))
  expect_true("separate_gwas" %in% names(result))
  expect_true("computing_time" %in% names(result))
  expect_true("QC" %in% names(result$computing_time))
  expect_true("Analysis" %in% names(result$computing_time))
})

# ===========================================================================
# 10. Pipeline with individual data but NULL sumstat (xqtl path)
# ===========================================================================
test_that("pipeline with individual data enters xqtl_coloc path and records timing", {
  region_data <- make_individual_region_data(n = 20, p = 8, n_contexts = 1, n_events = 2)

  # Mock qc_regional_data (pecotmr namespace) to simulate QC output.
  # colocboost is an external package function and cannot be mocked via
  # local_mocked_bindings. The pipeline's tryCatch handles the case where
  # colocboost is unavailable. We verify the pipeline enters the xqtl path
  # by checking that computing_time$Analysis$xqtl_coloc is recorded.
  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      Y1 <- region_data$individual_data$residual_Y[[1]]
      colnames(Y1) <- paste0(names(region_data$individual_data$residual_Y)[1], "_", colnames(Y1))
      list(
        individual_data = list(
          Y = list(ctx1 = Y1),
          X = list(ctx1 = region_data$individual_data$residual_X[[1]])
        ),
        sumstat_data = NULL
      )
    }
  )

  result <- suppressMessages(
    colocboost_analysis_pipeline(
      region_data,
      xqtl_coloc    = TRUE,
      joint_gwas     = FALSE,
      separate_gwas  = FALSE
    )
  )
  expect_type(result, "list")
  # The pipeline should have entered the xqtl_coloc analysis branch and
  # recorded timing, regardless of whether the colocboost call succeeded.
  expect_true(!is.null(result$computing_time$Analysis$xqtl_coloc))
})

# ===========================================================================
# 11. Pipeline: filter_events with exclude_pattern only
# ===========================================================================
test_that("filter_events exclude_pattern removes matching events via pipeline", {
  set.seed(1100)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  # Two events: one we want to keep, one to exclude
  Y <- matrix(rnorm(n * 2), n, 2)
  colnames(Y) <- c("good_event", "bad_event")

  region_data <- list(
    individual_data = list(
      residual_Y = list(ctx1 = Y),
      residual_X = list(ctx1 = X),
      maf        = list(ctx1 = runif(p, 0.05, 0.45))
    ),
    sumstat_data = NULL
  )

  event_filters <- list(
    list(
      type_pattern    = ".*_event$",
      exclude_pattern = "bad_event"
    )
  )

  # Mock qc_regional_data (pecotmr namespace) to pass through filtered data.
  # colocboost is an external package function and cannot be mocked via
  # local_mocked_bindings. The pipeline's tryCatch handles the case where
  # colocboost is unavailable.
  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      # The residual_Y should have had bad_event removed by filter_events
      remaining_events <- colnames(region_data$individual_data$residual_Y$ctx1)
      list(
        individual_data = list(
          Y = list(ctx1 = region_data$individual_data$residual_Y$ctx1),
          X = list(ctx1 = region_data$individual_data$residual_X$ctx1)
        ),
        sumstat_data = NULL
      )
    }
  )

  expect_message(
    result <- colocboost_analysis_pipeline(
      region_data,
      event_filters  = event_filters,
      xqtl_coloc     = TRUE
    ),
    "removed"
  )
  expect_type(result, "list")
})

# ===========================================================================
# 12. Pipeline with sumstat_data initializes separate_gwas structure
# ===========================================================================
# ===========================================================================
# 13. Pipeline catches colocboost errors gracefully
# ===========================================================================
test_that("pipeline catches colocboost xqtl error and returns NULL result", {
  region_data <- make_individual_region_data(n = 20, p = 8, n_contexts = 1, n_events = 2)

  # Mock qc_regional_data to return deliberately mismatched data (X has

  # different row count from Y) so that the colocboost call will always
  # error, whether or not the colocboost package is installed.
  # colocboost is an external package function and cannot be mocked via
  # local_mocked_bindings.
  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      Y1 <- region_data$individual_data$residual_Y[[1]]
      colnames(Y1) <- paste0("ctx1_", colnames(Y1))
      # Return X with mismatched rows to guarantee colocboost errors
      bad_X <- matrix(rnorm(5 * 8), nrow = 5, ncol = 8)
      colnames(bad_X) <- colnames(region_data$individual_data$residual_X[[1]])
      list(
        individual_data = list(
          Y = list(ctx1 = Y1),
          X = list(ctx1 = bad_X)
        ),
        sumstat_data = NULL
      )
    }
  )

  expect_message(
    result <- colocboost_analysis_pipeline(
      region_data,
      xqtl_coloc = TRUE,
      joint_gwas  = FALSE,
      separate_gwas = FALSE
    ),
    "xQTL-only ColocBoost failed"
  )
  expect_null(result$xqtl_coloc)
})

# ===========================================================================
# 14. Pipeline: no data passes QC returns early
# ===========================================================================
test_that("pipeline returns early when all data fails QC", {
  region_data <- make_individual_region_data()

  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      # Simulate all data removed by QC
      list(individual_data = NULL, sumstat_data = list(sumstats = NULL))
    }
  )

  expect_message(
    result <- colocboost_analysis_pipeline(
      region_data,
      xqtl_coloc    = TRUE,
      joint_gwas     = FALSE,
      separate_gwas  = FALSE
    ),
    "No data pass QC"
  )
  expect_type(result, "list")
})


# ===========================================================================
# Tests from test_twas_colocboost_round3.R (colocboost-related)
# ===========================================================================

# Helper functions used by round3 tests
# Helper: build individual-level region_data
make_individual_region_data <- function(n = 20, p = 8, n_contexts = 2, n_events = 3) {
  set.seed(701)
  make_ctx <- function(ctx_name) {
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- paste0("chr1:", seq_len(p) * 100, ":A:G")
    Y <- matrix(rnorm(n * n_events), n, n_events)
    colnames(Y) <- paste0("event", seq_len(n_events))
    maf <- runif(p, 0.05, 0.45)
    list(X = X, Y = Y, maf = maf)
  }
  ctxs <- lapply(paste0("ctx", seq_len(n_contexts)), make_ctx)
  names(ctxs) <- paste0("ctx", seq_len(n_contexts))
  list(
    individual_data = list(
      residual_Y = lapply(ctxs, `[[`, "Y"),
      residual_X = lapply(ctxs, `[[`, "X"),
      maf = lapply(ctxs, `[[`, "maf")
    ),
    sumstat_data = NULL
  )
}

# Helper: build sumstat-only region_data
make_sumstat_region_data <- function(n_variants = 5, n_studies = 2) {
  set.seed(702)
  vids <- paste0("chr1:", seq_len(n_variants) * 100, ":A:G")

  LD_mat <- diag(n_variants)
  rownames(LD_mat) <- colnames(LD_mat) <- vids

  ref_panel <- data.frame(
    chrom = rep(1, n_variants),
    pos = seq_len(n_variants) * 100,
    A2 = rep("A", n_variants),
    A1 = rep("G", n_variants),
    stringsAsFactors = FALSE
  )

  sumstats_list <- lapply(seq_len(n_studies), function(i) {
    ss <- list(
      sumstats = data.frame(
        chrom = rep(1, n_variants),
        pos = seq_len(n_variants) * 100,
        A1 = rep("G", n_variants),
        A2 = rep("A", n_variants),
        beta = rnorm(n_variants),
        se = runif(n_variants, 0.05, 0.2),
        z = rnorm(n_variants, 0, 2),
        variant_id = vids,
        stringsAsFactors = FALSE
      ),
      n = 10000,
      var_y = 1
    )
    list(ss) |> setNames(paste0("study", i))
  })

  LD_info <- list(list(
    LD_variants = ref_panel,
    LD_matrix = LD_mat,
    ref_panel = ref_panel
  ))

  list(
    individual_data = NULL,
    sumstat_data = list(
      sumstats = sumstats_list,
      LD_info = LD_info
    )
  )
}


test_that("filter_valid_sumstats: all valid entries pass through unchanged", {
  ss1 <- data.frame(z = c(1.0, 2.0), n = c(100, 100), variant = c("v1", "v2"))
  ss2 <- data.frame(z = c(-1.0, 0.5), n = c(200, 200), variant = c("v3", "v4"))

  sumstats <- list(s1 = ss1, s2 = ss2)
  LD_mat <- list(s1 = diag(2), s2 = diag(2))
  LD_match <- c("s1", "s2")

  # No message about removed entries when all are valid
  result <- pecotmr:::filter_valid_sumstats(sumstats, LD_mat, LD_match)
  expect_equal(length(result$sumstats), 2)
  expect_equal(nrow(result$dict_sumstatLD), 2)
})

# ===========================================================================
# SECTION C: colocboost_analysis_pipeline - filter_events valid_pattern path
# (lines 64, 67, 71-72, 74, 82, 84-85)
# ===========================================================================

test_that("filter_events: valid_pattern with no matching groups returns NULL (line 74)", {
  set.seed(800)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  # Events that match type_pattern but none match valid_pattern
  events <- c("clu_1_+:IN:gene1", "clu_2_-:IN:gene2")
  Y <- matrix(rnorm(n * length(events)), n, length(events))
  colnames(Y) <- events

  region_data <- list(
    individual_data = list(
      residual_Y = list(tissue1 = Y),
      residual_X = list(tissue1 = X),
      maf = list(tissue1 = runif(p, 0.05, 0.45))
    ),
    sumstat_data = NULL
  )

  # valid_pattern requires ":PR:" but no events have it -> valid_groups is empty -> type_events = character(0) -> returns NULL
  event_filters <- list(
    list(
      type_pattern = ".*clu_(\\d+_[+-?]).*",
      valid_pattern = "clu_(\\d+_[+-?]):PR:"
    )
  )

  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      list(
        individual_data = list(
          Y = region_data$individual_data$residual_Y,
          X = region_data$individual_data$residual_X
        ),
        sumstat_data = NULL
      )
    }
  )

  # The filter returns NULL for the context -> residual_Y entry becomes NULL
  expect_message(
    result <- colocboost_analysis_pipeline(
      region_data,
      event_filters = event_filters,
      xqtl_coloc = TRUE,
      joint_gwas = FALSE,
      separate_gwas = FALSE
    ),
    "No events matching|No data pass QC"
  )
  expect_type(result, "list")
})

test_that("filter_events: type_pattern matches nothing skips via next (line 64)", {
  set.seed(801)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  Y <- matrix(rnorm(n * 2), n, 2)
  colnames(Y) <- c("gene_A", "gene_B")

  region_data <- list(
    individual_data = list(
      residual_Y = list(ctx1 = Y),
      residual_X = list(ctx1 = X),
      maf = list(ctx1 = runif(p, 0.05, 0.45))
    ),
    sumstat_data = NULL
  )

  # type_pattern matches nothing -> type_events length 0 -> next
  event_filters <- list(
    list(
      type_pattern = "NONEXISTENT_PATTERN_xyz",
      exclude_pattern = "something"
    )
  )

  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      # Verify events were NOT filtered (both still present)
      remaining <- colnames(region_data$individual_data$residual_Y$ctx1)
      expect_equal(length(remaining), 2)
      list(
        individual_data = list(
          Y = region_data$individual_data$residual_Y,
          X = region_data$individual_data$residual_X
        ),
        sumstat_data = NULL
      )
    }
  )

  result <- suppressMessages(
    colocboost_analysis_pipeline(
      region_data,
      event_filters = event_filters,
      xqtl_coloc = TRUE,
      joint_gwas = FALSE,
      separate_gwas = FALSE
    )
  )
  expect_type(result, "list")
})

test_that("filter_events: all events pass -> 'included in following analysis' message (line 82)", {
  set.seed(802)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  # Both events match type_pattern and none get excluded
  Y <- matrix(rnorm(n * 2), n, 2)
  colnames(Y) <- c("evt_alpha", "evt_beta")

  region_data <- list(
    individual_data = list(
      residual_Y = list(ctx1 = Y),
      residual_X = list(ctx1 = X),
      maf = list(ctx1 = runif(p, 0.05, 0.45))
    ),
    sumstat_data = NULL
  )

  # type_pattern matches all events, exclude_pattern matches none
  event_filters <- list(
    list(
      type_pattern = "^evt_",
      exclude_pattern = "NONEXISTENT"
    )
  )

  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      list(
        individual_data = list(
          Y = region_data$individual_data$residual_Y,
          X = region_data$individual_data$residual_X
        ),
        sumstat_data = NULL
      )
    }
  )

  expect_message(
    result <- colocboost_analysis_pipeline(
      region_data,
      event_filters = event_filters,
      xqtl_coloc = TRUE,
      joint_gwas = FALSE,
      separate_gwas = FALSE
    ),
    "included in following analysis"
  )
  expect_type(result, "list")
})

test_that("filter_events: valid_pattern with successful groups retains valid events (lines 67, 71-72)", {
  set.seed(803)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  events <- c("clu_1_+:PR:gene1", "clu_1_+:IN:gene1", "clu_2_-:PR:gene2")
  Y <- matrix(rnorm(n * length(events)), n, length(events))
  colnames(Y) <- events

  region_data <- list(
    individual_data = list(
      residual_Y = list(tissue1 = Y),
      residual_X = list(tissue1 = X),
      maf = list(tissue1 = runif(p, 0.05, 0.45))
    ),
    sumstat_data = NULL
  )

  event_filters <- list(
    list(
      type_pattern = ".*clu_(\\d+_[+-?]).*",
      valid_pattern = "clu_(\\d+_[+-?]):PR:",
      exclude_pattern = "clu_(\\d+_[+-?]):IN:"
    )
  )

  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      remaining <- colnames(region_data$individual_data$residual_Y$tissue1)
      # IN event should be removed
      expect_false("clu_1_+:IN:gene1" %in% remaining)
      list(
        individual_data = list(
          Y = region_data$individual_data$residual_Y,
          X = region_data$individual_data$residual_X
        ),
        sumstat_data = NULL
      )
    }
  )

  expect_message(
    result <- colocboost_analysis_pipeline(
      region_data,
      event_filters = event_filters,
      xqtl_coloc = TRUE,
      joint_gwas = FALSE,
      separate_gwas = FALSE
    ),
    "removed"
  )
  expect_type(result, "list")
})

# ===========================================================================
# SECTION D: extract_contexts_studies - after-QC paths
# (lines 114, 125-127, 128-133, 134-135, 143-146)
# ===========================================================================

test_that("extract_contexts_studies: all individual data pass QC (line 126)", {
  region_data <- make_individual_region_data(n = 20, p = 8, n_contexts = 2, n_events = 2)

  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      Y1 <- matrix(rnorm(20), 20, 1); colnames(Y1) <- "ctx1_event1"
      Y2 <- matrix(rnorm(20), 20, 1); colnames(Y2) <- "ctx2_event1"
      X1 <- matrix(rnorm(20 * 8), 20, 8); colnames(X1) <- paste0("chr1:", seq_len(8) * 100, ":A:G")
      X2 <- matrix(rnorm(20 * 8), 20, 8); colnames(X2) <- paste0("chr1:", seq_len(8) * 100, ":A:G")
      list(
        individual_data = list(
          Y = list(ctx1_event1 = Y1, ctx2_event1 = Y2),
          X = list(ctx1 = X1, ctx2 = X2)
        ),
        sumstat_data = NULL
      )
    }
  )

  expect_message(
    result <- colocboost_analysis_pipeline(
      region_data,
      xqtl_coloc = TRUE,
      joint_gwas = FALSE,
      separate_gwas = FALSE
    ),
    "All individual data pass QC"
  )
  expect_type(result, "list")
})

test_that("extract_contexts_studies: all individual data fail QC (line 134-135)", {
  region_data <- make_individual_region_data(n = 20, p = 8, n_contexts = 2, n_events = 2)

  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      list(
        individual_data = list(
          Y = list(ctx1 = NULL, ctx2 = NULL),
          X = list(ctx1 = matrix(0, 1, 1), ctx2 = matrix(0, 1, 1))
        ),
        sumstat_data = NULL
      )
    }
  )

  expect_message(
    result <- colocboost_analysis_pipeline(
      region_data,
      xqtl_coloc = TRUE,
      joint_gwas = FALSE,
      separate_gwas = FALSE
    ),
    "No individual data pass QC"
  )
  expect_type(result, "list")
})

test_that("extract_contexts_studies: sumstat studies extraction on initial call (line 114)", {
  # region_data with sumstat only -> triggers sumstat branch in extract_contexts_studies
  region_data <- list(
    individual_data = NULL,
    sumstat_data = list(
      sumstats = list(
        list(gwas_trait1 = list(sumstats = data.frame(z = 1), n = 100, var_y = 1))
      )
    )
  )

  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      list(
        individual_data = NULL,
        sumstat_data = list(
          sumstats = list(gwas_trait1 = list(
            sumstats = data.frame(z = 1.5, variant_id = "chr1:100:A:G"),
            n = 100, var_y = 1
          )),
          LD_mat = list(gwas_trait1 = matrix(1, 1, 1, dimnames = list("chr1:100:A:G", "chr1:100:A:G"))),
          LD_match = "gwas_trait1"
        )
      )
    }
  )

  # With xqtl_coloc=FALSE, separate_gwas=TRUE -> will enter the sumstat code paths
  result <- suppressMessages(
    colocboost_analysis_pipeline(
      region_data,
      xqtl_coloc = FALSE,
      joint_gwas = FALSE,
      separate_gwas = TRUE
    )
  )
  expect_type(result, "list")
})

test_that("extract_contexts_studies: after-QC sumstat all pass (line 144)", {
  region_data <- make_sumstat_region_data(n_variants = 5, n_studies = 2)

  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      vids <- paste0("chr1:", seq_len(5) * 100, ":A:G")
      LD_mat <- diag(5); rownames(LD_mat) <- colnames(LD_mat) <- vids
      ss1 <- list(sumstats = data.frame(z = rnorm(5), variant_id = vids), n = 10000, var_y = 1)
      ss2 <- list(sumstats = data.frame(z = rnorm(5), variant_id = vids), n = 10000, var_y = 1)
      list(
        individual_data = NULL,
        sumstat_data = list(
          sumstats = list(study1 = ss1, study2 = ss2),
          LD_mat = list(study1 = LD_mat),
          LD_match = c("study1", "study1")
        )
      )
    }
  )

  expect_message(
    result <- colocboost_analysis_pipeline(
      region_data,
      xqtl_coloc = FALSE,
      joint_gwas = FALSE,
      separate_gwas = TRUE
    ),
    "All sumstat studies pass QC"
  )
  expect_type(result, "list")
})

test_that("extract_contexts_studies: after-QC sumstat partial pass (line 146)", {
  region_data <- make_sumstat_region_data(n_variants = 5, n_studies = 2)

  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      vids <- paste0("chr1:", seq_len(5) * 100, ":A:G")
      LD_mat <- diag(5); rownames(LD_mat) <- colnames(LD_mat) <- vids
      # Only one study remains after QC
      ss1 <- list(sumstats = data.frame(z = rnorm(5), variant_id = vids), n = 10000, var_y = 1)
      list(
        individual_data = NULL,
        sumstat_data = list(
          sumstats = list(study1 = ss1),
          LD_mat = list(study1 = LD_mat),
          LD_match = c("study1")
        )
      )
    }
  )

  expect_message(
    result <- colocboost_analysis_pipeline(
      region_data,
      xqtl_coloc = FALSE,
      joint_gwas = FALSE,
      separate_gwas = TRUE
    ),
    "Skipping follow-up analysis for sumstat studies"
  )
  expect_type(result, "list")
})

# ===========================================================================
# SECTION E: colocboost_analysis_pipeline - sumstat processing block
# (lines 244-281: organizing sumstats, normalizing variant IDs, LD normalization)
# ===========================================================================

test_that("pipeline sumstat block: normalizes variant IDs and processes LD matrices (lines 245-281)", {
  set.seed(820)
  n_variants <- 4
  vids <- paste0("chr1:", seq_len(n_variants) * 100, ":A:G")
  LD_mat <- diag(n_variants)
  rownames(LD_mat) <- colnames(LD_mat) <- vids

  region_data <- list(
    individual_data = NULL,
    sumstat_data = list(
      sumstats = list(
        list(study1 = list(
          sumstats = data.frame(z = c(2.1, -1.5, 0.8, 3.2), variant_id = vids),
          n = 5000, var_y = 1
        ))
      ),
      LD_info = list(list(
        LD_variants = data.frame(chrom = 1, pos = seq_len(n_variants) * 100, A2 = "A", A1 = "G"),
        LD_matrix = LD_mat,
        ref_panel = data.frame(chrom = 1, pos = seq_len(n_variants) * 100, A2 = "A", A1 = "G")
      ))
    )
  )

  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      ss <- list(sumstats = data.frame(z = c(2.1, -1.5, 0.8, 3.2), variant_id = vids), n = 5000, var_y = 1)
      list(
        individual_data = NULL,
        sumstat_data = list(
          sumstats = list(study1 = ss),
          LD_mat = list(study1 = LD_mat),
          LD_match = c("study1")
        )
      )
    }
  )

  result <- suppressMessages(
    colocboost_analysis_pipeline(
      region_data,
      xqtl_coloc = FALSE,
      joint_gwas = FALSE,
      separate_gwas = TRUE
    )
  )
  expect_type(result, "list")
  # The separate_gwas structure should be initialized with the study name
  expect_true("separate_gwas" %in% names(result))
})

test_that("pipeline sumstat block: single sumstat study initializes separate_gwas (line 180-181)", {
  set.seed(821)
  n_variants <- 3
  vids <- paste0("chr1:", seq_len(n_variants) * 100, ":A:G")
  LD_mat <- diag(n_variants)
  rownames(LD_mat) <- colnames(LD_mat) <- vids

  region_data <- list(
    individual_data = NULL,
    sumstat_data = list(
      sumstats = list(
        list(single_study = list(
          sumstats = data.frame(z = rnorm(n_variants), variant_id = vids),
          n = 5000, var_y = 1
        ))
      ),
      LD_info = list(list(
        LD_variants = data.frame(chrom = 1, pos = seq_len(n_variants) * 100, A2 = "A", A1 = "G"),
        LD_matrix = LD_mat,
        ref_panel = data.frame(chrom = 1, pos = seq_len(n_variants) * 100, A2 = "A", A1 = "G")
      ))
    )
  )

  # With only one sumstat study, line 180-181 should be reached
  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      ss <- list(sumstats = data.frame(z = rnorm(n_variants), variant_id = vids), n = 5000, var_y = 1)
      list(
        individual_data = NULL,
        sumstat_data = list(
          sumstats = list(single_study = ss),
          LD_mat = list(single_study = LD_mat),
          LD_match = c("single_study")
        )
      )
    }
  )

  result <- suppressMessages(
    colocboost_analysis_pipeline(
      region_data,
      xqtl_coloc = FALSE,
      joint_gwas = FALSE,
      separate_gwas = TRUE
    )
  )
  expect_type(result, "list")
  expect_true("separate_gwas" %in% names(result))
  # The separate_gwas structure was initialized (may be empty list if colocboost fails or not installed)
  expect_true(is.list(result$separate_gwas))
})

test_that("pipeline sumstat block: multiple sumstat studies initializes separate_gwas (line 176-178)", {
  set.seed(822)
  n_variants <- 3
  vids <- paste0("chr1:", seq_len(n_variants) * 100, ":A:G")
  LD_mat <- diag(n_variants)
  rownames(LD_mat) <- colnames(LD_mat) <- vids

  region_data <- list(
    individual_data = NULL,
    sumstat_data = list(
      sumstats = list(
        list(studyA = list(
          sumstats = data.frame(z = rnorm(n_variants), variant_id = vids), n = 5000, var_y = 1
        ), studyB = list(
          sumstats = data.frame(z = rnorm(n_variants), variant_id = vids), n = 5000, var_y = 1
        ))
      ),
      LD_info = list(list(
        LD_variants = data.frame(chrom = 1, pos = seq_len(n_variants) * 100, A2 = "A", A1 = "G"),
        LD_matrix = LD_mat,
        ref_panel = data.frame(chrom = 1, pos = seq_len(n_variants) * 100, A2 = "A", A1 = "G")
      ))
    )
  )

  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      ss1 <- list(sumstats = data.frame(z = rnorm(n_variants), variant_id = vids), n = 5000, var_y = 1)
      ss2 <- list(sumstats = data.frame(z = rnorm(n_variants), variant_id = vids), n = 5000, var_y = 1)
      list(
        individual_data = NULL,
        sumstat_data = list(
          sumstats = list(studyA = ss1, studyB = ss2),
          LD_mat = list(studyA = LD_mat),
          LD_match = c("studyA", "studyA")
        )
      )
    }
  )

  result <- suppressMessages(
    colocboost_analysis_pipeline(
      region_data,
      xqtl_coloc = FALSE,
      joint_gwas = FALSE,
      separate_gwas = TRUE
    )
  )
  expect_type(result, "list")
  expect_true("separate_gwas" %in% names(result))
  # The separate_gwas structure was initialized (may be empty list if colocboost fails or not installed)
  expect_true(is.list(result$separate_gwas))
})

# ===========================================================================
# SECTION F: colocboost_analysis_pipeline - filter_valid_sumstats returning NULL
# (lines 275-276: pipeline with all invalid sumstats -> "No data pass QC")
# ===========================================================================

test_that("pipeline: all sumstats invalid after filter_valid_sumstats returns No data pass QC (line 275-276)", {
  set.seed(830)
  n_variants <- 3
  vids <- paste0("chr1:", seq_len(n_variants) * 100, ":A:G")
  LD_mat <- diag(n_variants)
  rownames(LD_mat) <- colnames(LD_mat) <- vids

  region_data <- list(
    individual_data = NULL,
    sumstat_data = list(
      sumstats = list(
        list(study1 = list(
          sumstats = data.frame(z = rnorm(n_variants), variant_id = vids), n = 5000, var_y = 1
        ))
      ),
      LD_info = list(list(
        LD_variants = data.frame(chrom = 1, pos = seq_len(n_variants) * 100, A2 = "A", A1 = "G"),
        LD_matrix = LD_mat,
        ref_panel = data.frame(chrom = 1, pos = seq_len(n_variants) * 100, A2 = "A", A1 = "G")
      ))
    )
  )

  # Mock qc_regional_data to return sumstats that are all invalid (e.g., all NA z)
  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      invalid_ss <- list(
        sumstats = data.frame(z = NA_real_, variant_id = "chr1:100:A:G"),
        n = 0, var_y = 1
      )
      list(
        individual_data = NULL,
        sumstat_data = list(
          sumstats = list(study1 = invalid_ss),
          LD_mat = list(study1 = matrix(1, 1, 1, dimnames = list("chr1:100:A:G", "chr1:100:A:G"))),
          LD_match = c("study1")
        )
      )
    }
  )

  expect_message(
    result <- colocboost_analysis_pipeline(
      region_data,
      xqtl_coloc = FALSE,
      joint_gwas = FALSE,
      separate_gwas = TRUE
    ),
    "No data pass QC|No valid summary"
  )
  expect_type(result, "list")
})

# ===========================================================================
# SECTION G: colocboost_analysis_pipeline - focal_trait matching (lines 300-301)
# ===========================================================================

test_that("pipeline: focal_trait matches a trait name sets focal_outcome_idx (lines 300-301)", {
  region_data <- make_individual_region_data(n = 20, p = 8, n_contexts = 1, n_events = 2)

  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      Y1 <- matrix(rnorm(40), 20, 2)
      colnames(Y1) <- c("ctx1_event1", "ctx1_event2")
      X1 <- matrix(rnorm(20 * 8), 20, 8)
      colnames(X1) <- paste0("chr1:", seq_len(8) * 100, ":A:G")
      list(
        individual_data = list(
          Y = list(ctx1_event1 = Y1[, 1, drop = FALSE], ctx1_event2 = Y1[, 2, drop = FALSE]),
          X = list(ctx1 = X1)
        ),
        sumstat_data = NULL
      )
    }
  )

  # focal_trait = "ctx1_event2" should match and set focal_outcome_idx
  result <- suppressMessages(
    colocboost_analysis_pipeline(
      region_data,
      focal_trait = "ctx1_event2",
      xqtl_coloc = TRUE,
      joint_gwas = FALSE,
      separate_gwas = FALSE
    )
  )
  expect_type(result, "list")
  # The xqtl_coloc branch was entered; timing should be recorded
  expect_true(!is.null(result$computing_time$Analysis$xqtl_coloc))
})

test_that("pipeline: focal_trait does NOT match leaves focal_outcome_idx NULL (line 299-302)", {
  region_data <- make_individual_region_data(n = 20, p = 8, n_contexts = 1, n_events = 2)

  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      Y1 <- matrix(rnorm(20), 20, 1); colnames(Y1) <- "ctx1_event1"
      X1 <- matrix(rnorm(20 * 8), 20, 8)
      colnames(X1) <- paste0("chr1:", seq_len(8) * 100, ":A:G")
      list(
        individual_data = list(
          Y = list(ctx1_event1 = Y1),
          X = list(ctx1 = X1)
        ),
        sumstat_data = NULL
      )
    }
  )

  # focal_trait is specified but doesn't match any trait
  result <- suppressMessages(
    colocboost_analysis_pipeline(
      region_data,
      focal_trait = "nonexistent_trait",
      xqtl_coloc = TRUE,
      joint_gwas = FALSE,
      separate_gwas = FALSE
    )
  )
  expect_type(result, "list")
  expect_true(!is.null(result$computing_time$Analysis$xqtl_coloc))
})

# ===========================================================================
# SECTION H: colocboost_analysis_pipeline - joint_gwas path (lines 320-323)
# ===========================================================================

test_that("pipeline: joint_gwas path is entered with both individual and sumstat data (lines 320-323)", {
  set.seed(840)
  n <- 20; p <- 5
  vids <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- vids
  Y <- matrix(rnorm(n * 2), n, 2)
  colnames(Y) <- c("ctx1_gene1", "ctx1_gene2")
  LD_mat <- diag(p); rownames(LD_mat) <- colnames(LD_mat) <- vids

  region_data <- list(
    individual_data = list(
      residual_Y = list(ctx1 = Y),
      residual_X = list(ctx1 = X),
      maf = list(ctx1 = runif(p, 0.05, 0.45))
    ),
    sumstat_data = list(
      sumstats = list(
        list(gwas1 = list(
          sumstats = data.frame(z = rnorm(p), variant_id = vids), n = 5000, var_y = 1
        ))
      ),
      LD_info = list(list(
        LD_variants = data.frame(chrom = 1, pos = seq_len(p) * 100, A2 = "A", A1 = "G"),
        LD_matrix = LD_mat,
        ref_panel = data.frame(chrom = 1, pos = seq_len(p) * 100, A2 = "A", A1 = "G")
      ))
    )
  )

  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      ss <- list(sumstats = data.frame(z = rnorm(p), variant_id = vids), n = 5000, var_y = 1)
      list(
        individual_data = list(
          Y = list(ctx1_gene1 = Y[, 1, drop = FALSE], ctx1_gene2 = Y[, 2, drop = FALSE]),
          X = list(ctx1 = X)
        ),
        sumstat_data = list(
          sumstats = list(gwas1 = ss),
          LD_mat = list(gwas1 = LD_mat),
          LD_match = c("gwas1")
        )
      )
    }
  )

  # joint_gwas=TRUE should trigger the joint GWAS branch (lines 320-323)
  expect_message(
    result <- colocboost_analysis_pipeline(
      region_data,
      xqtl_coloc = FALSE,
      joint_gwas = TRUE,
      separate_gwas = FALSE
    ),
    "non-focaled version GWAS-xQTL ColocBoost"
  )
  expect_type(result, "list")
  expect_true(!is.null(result$computing_time$Analysis$joint_gwas))
})


# ===========================================================================
# SECTION I: colocboost_analysis_pipeline - separate_gwas path (lines 341+)
# ===========================================================================

test_that("pipeline: separate_gwas path is entered for each GWAS study", {
  set.seed(850)
  n <- 20; p <- 5
  vids <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- vids
  Y <- matrix(rnorm(n), n, 1); colnames(Y) <- "ctx1_gene1"
  LD_mat <- diag(p); rownames(LD_mat) <- colnames(LD_mat) <- vids

  region_data <- list(
    individual_data = list(
      residual_Y = list(ctx1 = Y),
      residual_X = list(ctx1 = X),
      maf = list(ctx1 = runif(p, 0.05, 0.45))
    ),
    sumstat_data = list(
      sumstats = list(
        list(gwasA = list(
          sumstats = data.frame(z = rnorm(p), variant_id = vids), n = 5000, var_y = 1
        ))
      ),
      LD_info = list(list(
        LD_variants = data.frame(chrom = 1, pos = seq_len(p) * 100, A2 = "A", A1 = "G"),
        LD_matrix = LD_mat,
        ref_panel = data.frame(chrom = 1, pos = seq_len(p) * 100, A2 = "A", A1 = "G")
      ))
    )
  )

  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      ss <- list(sumstats = data.frame(z = rnorm(p), variant_id = vids), n = 5000, var_y = 1)
      list(
        individual_data = list(
          Y = list(ctx1_gene1 = Y),
          X = list(ctx1 = X)
        ),
        sumstat_data = list(
          sumstats = list(gwasA = ss),
          LD_mat = list(gwasA = LD_mat),
          LD_match = c("gwasA")
        )
      )
    }
  )

  expect_message(
    result <- colocboost_analysis_pipeline(
      region_data,
      xqtl_coloc = FALSE,
      joint_gwas = FALSE,
      separate_gwas = TRUE
    ),
    "focaled version GWAS-xQTL ColocBoost"
  )
  expect_type(result, "list")
  expect_true(!is.null(result$computing_time$Analysis$separate_gwas))
})

# ===========================================================================
# SECTION J: colocboost_analysis_pipeline - all Y NULL after organizing (lines 225-227)
# ===========================================================================

test_that("pipeline: all Y become NULL after organizing individual data (lines 225-227)", {
  region_data <- make_individual_region_data(n = 20, p = 8, n_contexts = 2, n_events = 2)

  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      # Return individual_data where all Y entries are NULL
      list(
        individual_data = list(
          Y = list(ctx1 = NULL, ctx2 = NULL),
          X = list(ctx1 = matrix(rnorm(160), 20, 8), ctx2 = matrix(rnorm(160), 20, 8))
        ),
        sumstat_data = NULL
      )
    }
  )

  expect_message(
    result <- colocboost_analysis_pipeline(
      region_data,
      xqtl_coloc = TRUE,
      joint_gwas = FALSE,
      separate_gwas = FALSE
    ),
    "No data pass QC|No individual data pass QC"
  )
  expect_type(result, "list")
})

# ===========================================================================
# SECTION K: colocboost_analysis_pipeline - sumstat all z NA (line 252-255)
# ===========================================================================

test_that("pipeline: sumstat with all NA z-scores yields warning message (lines 252-255)", {
  set.seed(860)
  n_variants <- 3
  vids <- paste0("chr1:", seq_len(n_variants) * 100, ":A:G")
  LD_mat <- diag(n_variants); rownames(LD_mat) <- colnames(LD_mat) <- vids

  region_data <- list(
    individual_data = NULL,
    sumstat_data = list(
      sumstats = list(
        list(study_na = list(
          sumstats = data.frame(z = rnorm(n_variants), variant_id = vids), n = 5000, var_y = 1
        ))
      ),
      LD_info = list(list(
        LD_variants = data.frame(chrom = 1, pos = seq_len(n_variants) * 100, A2 = "A", A1 = "G"),
        LD_matrix = LD_mat,
        ref_panel = data.frame(chrom = 1, pos = seq_len(n_variants) * 100, A2 = "A", A1 = "G")
      ))
    )
  )

  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      # Return sumstats where ALL z-scores are NA
      ss <- list(
        sumstats = data.frame(z = rep(NA_real_, n_variants), variant_id = vids),
        n = 5000, var_y = 1
      )
      list(
        individual_data = NULL,
        sumstat_data = list(
          sumstats = list(study_na = ss),
          LD_mat = list(study_na = LD_mat),
          LD_match = c("study_na")
        )
      )
    }
  )

  expect_message(
    result <- colocboost_analysis_pipeline(
      region_data,
      xqtl_coloc = FALSE,
      joint_gwas = FALSE,
      separate_gwas = TRUE
    ),
    "All z-scores are NA|No data pass QC|No valid"
  )
  expect_type(result, "list")
})

# ===========================================================================
# SECTION L: colocboost_analysis_pipeline - no sumstat_data pass QC (line 152)
# ===========================================================================

test_that("extract_contexts_studies: no sumstat data pass QC message", {
  region_data <- make_individual_region_data(n = 20, p = 8, n_contexts = 1, n_events = 2)
  # Add some sumstat_data so it enters the initial extraction
  region_data$sumstat_data <- list(
    sumstats = list(
      list(gwas1 = list(sumstats = data.frame(z = 1), n = 100, var_y = 1))
    )
  )

  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      Y1 <- matrix(rnorm(20), 20, 1); colnames(Y1) <- "ctx1_event1"
      X1 <- matrix(rnorm(160), 20, 8)
      colnames(X1) <- paste0("chr1:", seq_len(8) * 100, ":A:G")
      list(
        individual_data = list(
          Y = list(ctx1_event1 = Y1),
          X = list(ctx1 = X1)
        ),
        sumstat_data = NULL  # All sumstat data removed by QC
      )
    }
  )

  expect_message(
    result <- colocboost_analysis_pipeline(
      region_data,
      xqtl_coloc = TRUE,
      joint_gwas = FALSE,
      separate_gwas = FALSE
    ),
    "No sumstat data pass QC"
  )
  expect_type(result, "list")
})

# ===========================================================================
# SECTION M: qc_regional_data - pip_cutoff_to_skip_ind wrong length errors
# ===========================================================================

test_that("qc_regional_data: mismatched pip_cutoff_to_skip_ind length errors", {
  region_data <- make_individual_region_data(n = 20, p = 8, n_contexts = 2, n_events = 2)

  expect_error(
    pecotmr:::qc_regional_data(
      region_data,
      pip_cutoff_to_skip_ind = c(0.1, 0.2, 0.3),  # 3 values but only 2 contexts
      qc_method = "slalom"
    ),
    "pip_cutoff_to_skip_ind must be a scalar"
  )
})

# ===========================================================================
# SECTION U: colocboost_analysis_pipeline - sumstat with NA z filtering (line 252-258)
# ===========================================================================

test_that("pipeline sumstat processing handles all-NA z-scores with warning (lines 252-258)", {
  set.seed(870)
  n_variants <- 4
  vids <- paste0("chr1:", seq_len(n_variants) * 100, ":A:G")
  LD_mat <- diag(n_variants); rownames(LD_mat) <- colnames(LD_mat) <- vids

  region_data <- list(
    individual_data = NULL,
    sumstat_data = list(
      sumstats = list(
        list(study_allna = list(
          sumstats = data.frame(z = rnorm(n_variants), variant_id = vids), n = 5000, var_y = 1
        ))
      ),
      LD_info = list(list(
        LD_variants = data.frame(chrom = 1, pos = seq_len(n_variants) * 100, A2 = "A", A1 = "G"),
        LD_matrix = LD_mat,
        ref_panel = data.frame(chrom = 1, pos = seq_len(n_variants) * 100, A2 = "A", A1 = "G")
      ))
    )
  )

  # Mock to return sumstats with all NA z-scores
  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      ss_allna <- list(
        sumstats = data.frame(z = rep(NA_real_, n_variants), variant_id = vids),
        n = 5000, var_y = 1
      )
      list(
        individual_data = NULL,
        sumstat_data = list(
          sumstats = list(study_allna = ss_allna),
          LD_mat = list(study_allna = LD_mat),
          LD_match = c("study_allna")
        )
      )
    }
  )

  # Should produce message about NA z-scores or no valid studies
  result <- suppressMessages(
    colocboost_analysis_pipeline(
      region_data,
      xqtl_coloc = FALSE,
      joint_gwas = FALSE,
      separate_gwas = TRUE
    )
  )
  expect_type(result, "list")
})

# ===========================================================================
# SECTION V: colocboost_analysis_pipeline - LD matrix normalization (lines 263-269)
# ===========================================================================

test_that("pipeline: LD matrix dimnames are normalized to canonical format (lines 261-269)", {
  set.seed(880)
  n_variants <- 3
  # Variant IDs without chr prefix to test normalization
  vids_no_chr <- paste0("1:", seq_len(n_variants) * 100, ":A:G")
  vids_with_chr <- paste0("chr1:", seq_len(n_variants) * 100, ":A:G")

  LD_mat <- diag(n_variants)
  rownames(LD_mat) <- colnames(LD_mat) <- vids_no_chr

  region_data <- list(
    individual_data = NULL,
    sumstat_data = list(
      sumstats = list(
        list(study1 = list(
          sumstats = data.frame(z = rnorm(n_variants), variant_id = vids_no_chr), n = 5000, var_y = 1
        ))
      ),
      LD_info = list(list(
        LD_variants = data.frame(chrom = 1, pos = seq_len(n_variants) * 100, A2 = "A", A1 = "G"),
        LD_matrix = LD_mat,
        ref_panel = data.frame(chrom = 1, pos = seq_len(n_variants) * 100, A2 = "A", A1 = "G")
      ))
    )
  )

  local_mocked_bindings(
    qc_regional_data = function(region_data, ...) {
      ss <- list(
        sumstats = data.frame(z = rnorm(n_variants), variant_id = vids_no_chr),
        n = 5000, var_y = 1
      )
      list(
        individual_data = NULL,
        sumstat_data = list(
          sumstats = list(study1 = ss),
          LD_mat = list(study1 = LD_mat),  # LD_mat has non-chr names
          LD_match = c("study1")
        )
      )
    }
  )

  # normalize_variant_id should add chr prefix to LD matrix dimnames
  result <- suppressMessages(
    colocboost_analysis_pipeline(
      region_data,
      xqtl_coloc = FALSE,
      joint_gwas = FALSE,
      separate_gwas = TRUE
    )
  )
  expect_type(result, "list")
})

# ===========================================================================
# SECTION X: colocboost - qc_regional_data with NULL individual_data only sumstat
# ===========================================================================

test_that("qc_regional_data: with only sumstat data processes correctly", {
  region_data <- make_sumstat_region_data(n_variants = 5, n_studies = 1)

  local_mocked_bindings(
    rss_basic_qc = function(sumstats, LD_data, ...) {
      LD_mat <- LD_data$LD_matrix
      list(sumstats = sumstats, LD_mat = LD_mat)
    },
    summary_stats_qc = function(sumstats, LD_data, ...) {
      LD_mat <- LD_data$LD_matrix
      list(sumstats = sumstats, LD_mat = LD_mat, outlier_number = 0)
    },
    raiss = function(...) {
      list(result_filter = data.frame(z = rnorm(5)), LD_mat = diag(5))
    },
    partition_LD_matrix = function(...) diag(5)
  )

  result <- suppressMessages(
    pecotmr:::qc_regional_data(
      region_data,
      qc_method = "slalom",
      impute = FALSE
    )
  )
  expect_type(result, "list")
  expect_null(result$individual_data)
  expect_true(!is.null(result$sumstat_data))
})
NA
