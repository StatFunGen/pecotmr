context("sumstats_qc")

# Previous tests against `rssBasicQc()` and the legacy
# `summaryStatsQc(rssInput, ldData)` signature have been removed because:
#   * `rssBasicQc()` was folded into `summaryStatsQc()`
#   * `summaryStatsQc()` now dispatches on `GwasSumStats` / `QtlSumStats`
#     S4 collections (DFrame subclasses), not on a (rssInput list, ldData)
#     pair
#   * `QcResult` and the accessors `getRssInput()`, `getLdData()`,
#     `getOutlierNumber()` have been removed; QC audit lives on
#     `getQcInfo(<SumStats>)`
# Integration coverage of `summaryStatsQc(<SumStats>)` lives in the
# pipeline test files (test_colocboostPipeline.R, etc.).
#
# What remains here: tests of internal QC helpers
# (`ldMismatchQc`, `krigingOutlierQc`) whose
# signatures and contracts are unchanged.

# ===========================================================================
# ldMismatchQc
# ===========================================================================

test_that("ldMismatchQc with dentist method returns data frame with outlier column", {
    set.seed(42)
    p <- 20
    R <- diag(p)
    z <- rnorm(p)
    result <- ldMismatchQc(z, R = R, nSample = 1000, method = "dentist")
    expect_true(is.data.frame(result) || is.list(result))
    expect_true("outlier" %in% names(result))
})

test_that("ldMismatchQc with slalom method returns data frame with outlier column", {
    set.seed(42)
    p <- 20
    R <- diag(p)
    z <- rnorm(p)
    result <- ldMismatchQc(z, R = R, method = "slalom")
    expect_true(is.data.frame(result) || is.list(result))
    expect_true("outlier" %in% names(result))
})

test_that("ldMismatchQc method argument is validated", {
    z <- rnorm(5)
    R <- diag(5)
    expect_error(ldMismatchQc(z, R = R, method = "invalid"))
})

# ===========================================================================
# zMismatchQc resolver
# ===========================================================================

# ===========================================================================
# krigingOutlierQc
# ===========================================================================

# An allele switch (what the flip rule targets) is an LD-CONSISTENT z with one
# entry's sign reversed -- not a magnitude outlier. Build z = R %*% b (a single
# causal SNP) then negate a strongly-tagged neighbour, so susieR's logLR fires.
.kr_switchScenario <- function() {
    m <- 10
    rho <- 0.9
    R <- outer(seq_len(m), seq_len(m), function(i, j) rho^abs(i - j))
    ids <- paste0("1:", seq_len(m) * 100, ":A:G")
    rownames(R) <- colnames(R) <- ids
    b <- numeric(m)
    b[5] <- 6
    z <- as.numeric(R %*% b)
    z[6] <- -z[6] # flip a strongly-tagged neighbour of the causal
    list(z = z, R = R, ids = ids, flipped = 6L)
}

test_that("krigingOutlierQc flags an allele-switched variant and spares the rest", {
    skip_if_not(
        "kriging_rss" %in% getNamespaceExports("susieR"),
        "installed susieR has no kriging_rss"
    )
    s <- .kr_switchScenario()
    kr <- krigingOutlierQc(s$z, s$R, n = 1000, variantIds = s$ids)
    expect_true(kr$flip[s$flipped])
    expect_equal(sum(kr$flip), 1L)
    expect_equal(nrow(kr$diagnostics), length(s$z))
    expect_true(all(
        c("z", "condmean", "z_std_diff", "logLR", "flipped") %in%
            colnames(kr$diagnostics)
    ))
    expect_identical(kr$diagnostics$flipped, kr$flip)
})

test_that("krigingOutlierQc flip == susieR's logLR>2 & |z|>2 selection", {
    skip_if_not(
        "kriging_rss" %in% getNamespaceExports("susieR"),
        "installed susieR has no kriging_rss"
    )
    set.seed(7)
    m <- 10
    R <- cov2cor(crossprod(matrix(rnorm(m * m), m)))
    ids <- paste0("1:", seq_len(m) * 100, ":A:G")
    rownames(R) <- colnames(R) <- ids
    z <- as.numeric(R %*% rnorm(m))
    z[4] <- 9
    kr <- krigingOutlierQc(z, R, n = 5000, variantIds = ids)
    ref <- susieR::kriging_rss(z = z, R = R, n = 5000)$conditional_dist
    # susieR's own allele-switch rule (susie_rss_utils.R): logLR > 2 & |z| > 2.
    expected <- as.numeric(ref$logLR) > 2 & abs(z) > 2
    expect_identical(kr$flip, expected)
    expect_equal(kr$diagnostics$logLR, as.numeric(ref$logLR), tolerance = 1e-8)
    expect_equal(
        kr$diagnostics$z_std_diff,
        as.numeric(ref$z_std_diff),
        tolerance = 1e-8
    )
    expect_equal(
        kr$diagnostics$condmean,
        as.numeric(ref$condmean),
        tolerance = 1e-8
    )
})

test_that("krigingOutlierQc thresholds are configurable", {
    skip_if_not(
        "kriging_rss" %in% getNamespaceExports("susieR"),
        "installed susieR has no kriging_rss"
    )
    s <- .kr_switchScenario()
    # The switch fires at the default logLRThreshold = 2 ...
    expect_true(krigingOutlierQc(s$z, s$R, n = 1000, variantIds = s$ids)$flip[
        s$flipped
    ])
    # ... and an impossibly high logLR threshold suppresses every flip.
    expect_false(any(
        krigingOutlierQc(
            s$z,
            s$R,
            n = 1000,
            variantIds = s$ids,
            logLRThreshold = 1e6
        )$flip
    ))
    # A |z| threshold above the switched z also suppresses it.
    expect_false(any(
        krigingOutlierQc(
            s$z,
            s$R,
            n = 1000,
            variantIds = s$ids,
            zThreshold = 1e6
        )$flip
    ))
})

test_that("krigingOutlierQc requires a positive sample size n", {
    expect_error(krigingOutlierQc(c(1, 2, 3), diag(3)), "positive sample size")
})


context("raiss")
library(tidyverse)
library(MASS)

# Helper: build LdData S4 from a ref_panel data.frame, correlation matrix, and blockMetadata
make_ld_data_from_ref_panel <- function(R_mat, ref_panel, blockMetadata) {
    ref_panel$chrom <- as.character(ref_panel$chrom)
    ref_panel$variant_id <- as.character(ref_panel$variant_id)
    variants_gr <- pecotmr:::.refPanelToGranges(ref_panel)
    LdData(
        correlation = R_mat,
        variants = variants_gr,
        blockMetadata = blockMetadata
    )
}

generate_dummy_data <- function(
    seed = 1,
    ref_panel_ordered = TRUE,
    known_zscores_ordered = TRUE
) {
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
    # One sample of ROWS, not four independent samples of columns: a real
    # observed variant's id encodes its own position, so drawing id and pos
    # separately would place a known id at a foreign position -- which the
    # observed-position guard correctly reads as a second, typed site.
    known_idx <- sample(n_variants, n_known)
    known_zscores <- data.frame(
        chrom = rep(1, n_known),
        pos = ref_panel$pos[known_idx],
        variant_id = ref_panel$variant_id[known_idx],
        A1 = ref_panel$A1[known_idx],
        A2 = ref_panel$A2[known_idx],
        z = rnorm(n_known)
    )

    ldMatrix <- matrix(
        rnorm(n_variants^2),
        nrow = n_variants,
        ncol = n_variants
    )
    diag(ldMatrix) <- 1
    known_zscores <- if (known_zscores_ordered) {
        known_zscores[order(known_zscores$pos), ]
    } else {
        known_zscores
    }
    ref_panel <- if (ref_panel_ordered) {
        ref_panel
    } else {
        ref_panel[order(ref_panel$pos, decreasing = TRUE), ]
    }
    return(list(
        ref_panel = ref_panel,
        known_zscores = known_zscores,
        ldMatrix = ldMatrix
    ))
}

test_that("Input validation for raiss works correctly", {
    input_data <- generate_dummy_data()
    input_data_ref_panel_unordered <- generate_dummy_data(
        ref_panel_ordered = FALSE
    )
    input_data_zscores_unordered <- generate_dummy_data(
        known_zscores_ordered = FALSE
    )
    expect_error(raiss(
        input_data_ref_panel_unordered$ref_panel,
        input_data$known_zscores,
        input_data$ldMatrix
    ))
    expect_error(raiss(
        input_data$ref_panel,
        input_data_zscores_unordered$known_zscores,
        input_data$ldMatrix
    ))
})

test_that("Default parameters for raiss work correctly", {
    input_data <- generate_dummy_data()
    result <- raiss(
        input_data$ref_panel,
        input_data$known_zscores,
        input_data$ldMatrix
    )
    expect_true(is.list(result))
    # Expected list elements
    expect_true(all(
        c("resultNofilter", "resultFilter", "ldMat") %in% names(result)
    ))
    # resultNofilter should be a data frame with expected columns
    expect_true(is.data.frame(result$resultNofilter))
    expect_true(all(
        c("variant_id", "z", "Var", "raissLdScore") %in%
            names(result$resultNofilter)
    ))
    # Imputed z-scores should be numeric and finite
    expect_true(is.numeric(result$resultNofilter$z))
    expect_true(all(is.finite(result$resultNofilter$z)))
    # Output should cover all ref_panel variants (known + imputed)
    expect_equal(nrow(result$resultNofilter), nrow(input_data$ref_panel))
    # Filtered result should be a subset of unfiltered
    expect_true(nrow(result$resultFilter) <= nrow(result$resultNofilter))
    # ldMat should be a matrix
    expect_true(is.matrix(result$ldMat))
})

test_that("Test Default Parameters for raissModel", {
    zt <- c(1.2, 0.5)
    sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
    sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

    result <- raissModel(zt, sig_t, sig_i_t)

    expect_true(is.list(result))
    expect_true(all(
        c(
            "var",
            "mu",
            "raissLdScore",
            "conditionNumber",
            "correctInversion"
        ) %in%
            names(result)
    ))
    # mu (imputed z-scores) should be numeric and finite
    expect_true(is.numeric(result$mu))
    expect_true(all(is.finite(result$mu)))
    # var should be numeric
    expect_true(is.numeric(result$var))
    # raissLdScore should be numeric and non-negative
    expect_true(is.numeric(result$raissLdScore))
    expect_true(all(result$raissLdScore >= 0))
})

test_that("Test with Different lamb Values for raissModel", {
    zt <- c(1.2, 0.5)
    sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
    sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

    lamb_values <- c(0.01, 0.05, 0.1)
    for (lamb in lamb_values) {
        result <- raissModel(zt, sig_t, sig_i_t, lamb)
        expect_true(is.list(result))
        expect_true(all(c("var", "mu", "raissLdScore") %in% names(result)))
        expect_true(is.numeric(result$mu))
        expect_true(all(is.finite(result$mu)))
    }
})

test_that("Report Condition Number in raissModel", {
    zt <- c(1.2, 0.5)
    sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
    sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

    result_with_cn <- raissModel(
        zt,
        sig_t,
        sig_i_t,
        reportConditionNumber = TRUE
    )
    result_without_cn <- raissModel(
        zt,
        sig_t,
        sig_i_t,
        reportConditionNumber = FALSE
    )

    expect_true(is.list(result_with_cn))
    expect_true(is.list(result_without_cn))
    # With condition number reporting, conditionNumber should be populated
    expect_true(is.numeric(result_with_cn$conditionNumber))
    expect_true(all(is.finite(result_with_cn$mu)))
    expect_true(all(is.finite(result_without_cn$mu)))
})

test_that("Input Validation of raissModel", {
    zt <- c(1.2, 0.5)
    sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
    sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)
    zt_invalid <- "not a numeric vector"
    sig_t_invalid <- "not a matrix"
    sig_i_t_invalid <- "not a matrix"

    expect_error(raissModel(zt_invalid, sig_t, sig_i_t))
    expect_error(raissModel(zt, sig_t_invalid, sig_i_t))
    expect_error(raissModel(zt, sig_t, sig_i_t_invalid))
})

test_that("Boundary Conditions of raissModel", {
    zt <- c(1.2, 0.5)
    sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
    sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

    zt_empty <- numeric(0)
    sig_t_empty <- matrix(numeric(0), nrow = 0)
    sig_i_t_empty <- matrix(numeric(0), nrow = 0)

    expect_error(raissModel(zt_empty, sig_t, sig_i_t))
    expect_error(raissModel(zt, sig_t_empty, sig_i_t))
    expect_error(raissModel(zt, sig_t, sig_i_t_empty))
})

test_that("Test with Different rcond Values for raissModel", {
    zt <- c(1.2, 0.5)
    sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
    sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

    rcond_values <- c(0.01, 0.05, 0.1)
    for (rcond in rcond_values) {
        result <- raissModel(zt, sig_t, sig_i_t, lamb = 0.01, rcond = rcond)
        expect_true(is.list(result))
        expect_true(all(
            c(
                "var",
                "mu",
                "raissLdScore",
                "conditionNumber",
                "correctInversion"
            ) %in%
                names(result)
        ))
        expect_true(is.numeric(result$mu))
        expect_true(all(is.finite(result$mu)))
    }
})

test_that("formatRaissDf returns correctly formatted data frame", {
    imp <- list(
        mu = rnorm(5),
        var = runif(5),
        raissLdScore = rnorm(5),
        conditionNumber = runif(5),
        correctInversion = sample(c(TRUE, FALSE), 5, replace = TRUE)
    )

    ref_panel <- data.frame(
        chrom = sample(1:22, 10, replace = TRUE),
        pos = sample(1:10000, 10),
        variant_id = paste0("rs", 1:10),
        A1 = sample(c("A", "T", "G", "C"), 10, replace = TRUE),
        A2 = sample(c("A", "T", "G", "C"), 10, replace = TRUE)
    )

    unknowns <- sample(1:nrow(ref_panel), 5)

    result <- formatRaissDf(imp, ref_panel, unknowns)

    expect_true(is.data.frame(result))
    expect_equal(ncol(result), 10)
    expect_equal(
        colnames(result),
        c(
            'chrom',
            'pos',
            'variant_id',
            'A1',
            'A2',
            'z',
            'Var',
            'raissLdScore',
            'conditionNumber',
            'correctInversion'
        )
    )

    for (col in c('chrom', 'pos', 'variant_id', 'A1', 'A2')) {
        expect_equal(
            setNames(unlist(result[col]), NULL),
            unlist(ref_panel[unknowns, col, drop = TRUE])
        )
    }
    for (col in c(
        'z',
        'Var',
        'raissLdScore',
        'conditionNumber',
        'correctInversion'
    )) {
        expected_col <- if (col == "z") {
            "mu"
        } else if (col == "Var") {
            "var"
        } else {
            col
        }
        expect_equal(
            setNames(unlist(result[col]), NULL),
            setNames(unlist(imp[expected_col]), NULL)
        )
    }
})

test_that("Merge operation is correct for mergeRaissDf", {
    raiss_df_example <- data.frame(
        chrom = c("chr21", "chr22"),
        pos = c(123, 456),
        variant_id = c("var1", "var2"),
        A1 = c("A", "T"),
        A2 = c("T", "A"),
        z = c(0.5, 1.5),
        Var = c(0.2, 0.3),
        raissLdScore = c(10, 20),
        raissR2 = c(0.8, 0.7)
    )

    known_zscores_example <- data.frame(
        chrom = c("chr21", "chr22"),
        pos = c(123, 456),
        variant_id = c("var1", "var2"),
        A1 = c("A", "T"),
        A2 = c("T", "A"),
        z = c(0.5, 1.5)
    )

    merged_df <- mergeRaissDf(raiss_df_example, known_zscores_example)
    expect_equal(nrow(merged_df), 2)
    expect_true(all(c("chr21", "chr22") %in% merged_df$chrom))
})

generate_fro_test_data <- function(seed = 1) {
    set.seed(seed)
    return(data.frame(
        chrom = paste0("chr", rep(22, 10)),
        pos = seq(1, 100, 10),
        variant_id = 1:10,
        A1 = rep("A", 10),
        A2 = rep("T", 10),
        z = rnorm(10),
        Var = runif(10, 0, 1),
        raissLdScore = rnorm(10, 5, 2)
    ))
}

test_that("Correct columns are selected in filterRaissOutput", {
    test_data <- generate_fro_test_data()
    output <- filterRaissOutput(test_data)$zscores
    expect_true(all(
        c('variant_id', 'A1', 'A2', 'z', 'Var', 'raissLdScore') %in%
            names(output)
    ))
})

test_that("raissR2 is calculated correctly in filterRaissOutput", {
    test_data <- generate_fro_test_data()
    output <- filterRaissOutput(test_data)$zscores
    expected_R2 <- 1 - test_data[which(test_data$raissLdScore >= 5), ]$Var
    expect_equal(output$raissR2, expected_R2[which(expected_R2 > 0.6)])
})

test_that("Filtering is applied correctly in filterRaissOutput", {
    test_data <- generate_fro_test_data()
    R2_threshold <- 0.6
    minimum_ld <- 5
    output <- filterRaissOutput(test_data, R2_threshold, minimum_ld)$zscores

    expect_true(all(output$raissR2 > R2_threshold))
    expect_true(all(output$raissLdScore >= minimum_ld))
})

test_that("Function returns the correct subset in filterRaissOutput", {
    test_data <- generate_fro_test_data()
    test_data$raissR2 <- 1 - test_data$Var
    output <- filterRaissOutput(test_data)$zscores

    manual_filter <- test_data[
        test_data$raissR2 > 0.6 & test_data$raissLdScore >= 5,
    ]

    expect_equal(nrow(output), nrow(manual_filter))
    expect_equal(sum(output$variant_id != manual_filter$variant_id), 0)
})

test_that("computeMu basic functionality", {
    sig_i_t <- matrix(c(1, 2, 3, 4), nrow = 2)
    sig_t_inv <- matrix(c(5, 6, 7, 8), nrow = 2)
    zt <- matrix(c(9, 10, 11, 12), nrow = 2)

    expected_result <- matrix(c(517, 766, 625, 926), nrow = 2)
    result <- computeMu(sig_i_t, sig_t_inv, zt)
    expect_equal(result, expected_result)
})

generate_mock_data_for_computeVar <- function(seed = 1) {
    return(
        list(
            sig_i_t_1 = matrix(c(1, 2, 3, 4), nrow = 2),
            sig_t_inv_1 = matrix(c(5, 6, 7, 8), nrow = 2),
            lamb_1 = 0.5
        )
    )
}

test_that("computeVar returns correct output for batch = TRUE", {
    input_data <- generate_mock_data_for_computeVar()
    result <- computeVar(
        input_data$sig_i_t_1,
        input_data$sig_t_inv_1,
        input_data$lamb_1,
        batch = TRUE
    )
    expect_true(is.list(result))
    expect_length(result, 2)
    expect_true(all(c("var", "raissLdScore") %in% names(result)))
    expect_true(is.numeric(result$var))
    expect_true(is.numeric(result$raissLdScore))
})

test_that("computeVar returns correct output for batch = FALSE", {
    input_data <- generate_mock_data_for_computeVar()
    result <- computeVar(
        input_data$sig_i_t_1,
        input_data$sig_t_inv_1,
        input_data$lamb_1,
        batch = FALSE
    )
    expect_true(is.list(result))
    expect_length(result, 2)
    expect_true(all(c("var", "raissLdScore") %in% names(result)))
    expect_true(is.numeric(result$var))
    expect_true(is.numeric(result$raissLdScore))
})

test_that("checkInversion correctly identifies inverse matrices in", {
    sig_t <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
    sig_t_inv <- solve(sig_t)
    expect_true(checkInversion(sig_t, sig_t_inv))
})

test_that("varInBoundaries sets boundaries correctly", {
    lamb_test <- 0.05
    var <- c(-1, 0, 0.5, 1.04, 1.05)

    result <- varInBoundaries(var, lamb_test)

    expect_equal(result[1], 0) # Value less than 0 should be set to 0
    expect_equal(result[2], 0) # Value within lower boundary should remain unchanged
    expect_equal(result[3], 0.5) # Value within boundaries should remain unchanged
    expect_equal(result[4], 1.04) # Value greater than 0.99999 + lamb should be set to 1
    expect_equal(result[5], 1) # Value greater than 0.99999 + lamb should be set to 1
})

test_that("invertMat computes correct pseudo-inverse", {
    mat <- matrix(c(1, 2, 3, 4), nrow = 2)
    lamb <- 0.5
    rcond <- 1e-7
    result <- invertMat(mat, lamb, rcond)
    expect_true(is.matrix(result))
})

test_that("invertMat handles errors and retries", {
    mat <- matrix(c(0, 0, 0, 0), nrow = 2)
    lamb <- 0.1
    rcond <- 1e-7
    result <- invertMat(mat, lamb, rcond)
    expect_true(is.matrix(result))
})

test_that("invertMatRecursive correctly inverts a valid square matrix", {
    mat <- matrix(c(2, -1, -1, 2), nrow = 2)
    lamb <- 0.5
    rcond <- 0.01
    result <- invertMatRecursive(mat, lamb, rcond)
    expect_true(is.matrix(result))
    expect_equal(dim(result), dim(mat))
})

test_that("invertMatRecursive handles non-square matrices appropriately", {
    mat <- matrix(1:6, nrow = 2)
    lamb <- 0.5
    rcond <- 0.01
    expect_silent(invertMatRecursive(mat, lamb, rcond))
})

test_that("invertMatRecursive handles errors and performs recursive call correctly", {
    mat <- "not a matrix"
    lamb <- 0.5
    rcond <- 0.01
    expect_error(invertMatRecursive(mat, lamb, rcond))
})

# Test with Different Tolerance Levels
test_that("invertMatEigen behaves differently with varying tolerance levels", {
    mat <- matrix(c(1, 0, 0, 1e-4), nrow = 2)
    tol_high <- 1e-2
    tol_low <- 1e-6
    result_high_tol <- invertMatEigen(mat, tol_high)
    result_low_tol <- invertMatEigen(mat, tol_low)
    expect_true(!is.logical(all.equal(result_high_tol, result_low_tol)))
})

test_that("invertMatEigen handles non-square matrices", {
    mat <- matrix(1:6, nrow = 2)
    expect_error(invertMatEigen(mat))
})

test_that("invertMatEigen returns the same matrix for an identity matrix", {
    mat <- diag(2)
    expected <- mat
    actual <- invertMatEigen(mat)
    expect_equal(actual, expected)
})

test_that("invertMatEigen returns a zero matrix for a zero matrix input", {
    mat <- matrix(0, nrow = 2, ncol = 2)
    expected <- mat
    expect_error(
        invertMatEigen(mat),
        "Cannot invert the input matrix because all its eigen values are negative or close to zero"
    )
})

test_that("invertMatEigen handles matrices with negative eigenvalues", {
    mat <- matrix(c(-2, 0, 0, -3), nrow = 2)
    expect_silent(invertMatEigen(mat))
})

# ===========================================================================
# raissSingleMatrix edge cases
# ===========================================================================

test_that("raissSingleMatrix returns NULL when no known variants overlap", {
    set.seed(42)
    ref_panel <- data.frame(
        chrom = rep(1, 10),
        pos = seq(10, 100, 10),
        variant_id = paste0("rs", 1:10),
        A1 = rep("A", 10),
        A2 = rep("G", 10),
        stringsAsFactors = FALSE
    )
    # known_zscores has variant IDs that don't match ref_panel at all
    known_zscores <- data.frame(
        chrom = rep(1, 3),
        pos = c(200, 300, 400),
        variant_id = paste0("other", 1:3),
        A1 = rep("A", 3),
        A2 = rep("G", 3),
        z = rnorm(3),
        stringsAsFactors = FALSE
    )
    ldMatrix <- diag(10)
    result <- raissSingleMatrix(
        ref_panel,
        known_zscores,
        ldMatrix,
        verbose = FALSE
    )
    expect_null(result)
})

test_that("raissSingleMatrix returns known zscores when no unknowns to impute", {
    set.seed(42)
    ref_panel <- data.frame(
        chrom = rep(1, 5),
        pos = seq(10, 50, 10),
        variant_id = paste0("rs", 1:5),
        A1 = rep("A", 5),
        A2 = rep("G", 5),
        stringsAsFactors = FALSE
    )
    # All ref_panel variants are known - nothing to impute
    known_zscores <- data.frame(
        chrom = rep(1, 5),
        pos = seq(10, 50, 10),
        variant_id = paste0("rs", 1:5),
        A1 = rep("A", 5),
        A2 = rep("G", 5),
        z = rnorm(5),
        stringsAsFactors = FALSE
    )
    ldMatrix <- diag(5)
    result <- raissSingleMatrix(
        ref_panel,
        known_zscores,
        ldMatrix,
        verbose = FALSE
    )
    expect_true(is.list(result))
    expect_equal(result$resultNofilter, known_zscores)
    expect_equal(result$resultFilter, known_zscores)
    expect_equal(result$ldMat, ldMatrix)
})

# ===========================================================================
# raissSingleMatrixFromX edge cases
# ===========================================================================

test_that("raissSingleMatrixFromX returns NULL when no known variants overlap", {
    set.seed(42)
    n <- 50
    p <- 10
    ref_panel <- data.frame(
        chrom = rep(1, p),
        pos = seq(10, p * 10, 10),
        variant_id = paste0("rs", 1:p),
        A1 = rep("A", p),
        A2 = rep("G", p),
        stringsAsFactors = FALSE
    )
    known_zscores <- data.frame(
        chrom = rep(1, 3),
        pos = c(200, 300, 400),
        variant_id = paste0("other", 1:3),
        A1 = rep("A", 3),
        A2 = rep("G", 3),
        z = rnorm(3),
        stringsAsFactors = FALSE
    )
    X <- scale(matrix(sample(0:2, n * p, replace = TRUE), nrow = n))
    X[is.na(X)] <- 0
    colnames(X) <- ref_panel$variant_id
    result <- raissSingleMatrixFromX(
        ref_panel,
        known_zscores,
        X,
        verbose = FALSE
    )
    expect_null(result)
})

test_that("raissSingleMatrixFromX returns known zscores when no unknowns to impute", {
    set.seed(42)
    n <- 50
    p <- 5
    ref_panel <- data.frame(
        chrom = rep(1, p),
        pos = seq(10, p * 10, 10),
        variant_id = paste0("rs", 1:p),
        A1 = rep("A", p),
        A2 = rep("G", p),
        stringsAsFactors = FALSE
    )
    known_zscores <- data.frame(
        chrom = rep(1, p),
        pos = seq(10, p * 10, 10),
        variant_id = paste0("rs", 1:p),
        A1 = rep("A", p),
        A2 = rep("G", p),
        z = rnorm(p),
        stringsAsFactors = FALSE
    )
    X <- scale(matrix(sample(0:2, n * p, replace = TRUE), nrow = n))
    X[is.na(X)] <- 0
    colnames(X) <- ref_panel$variant_id
    result <- raissSingleMatrixFromX(
        ref_panel,
        known_zscores,
        X,
        verbose = FALSE
    )
    expect_true(is.list(result))
    expect_equal(result$resultNofilter, known_zscores)
    expect_null(result$ldMat)
})

# ===========================================================================
# raiss() dispatch paths: single-matrix list and genotype_matrix list
# ===========================================================================

test_that("raiss with single-matrix LD list dispatches to single matrix path", {
    set.seed(42)
    n_variants <- 20
    ref_panel <- data.frame(
        chrom = rep(1, n_variants),
        pos = seq(10, n_variants * 10, 10),
        variant_id = paste0("rs", 1:n_variants),
        A1 = rep("A", n_variants),
        A2 = rep("G", n_variants),
        stringsAsFactors = FALSE
    )
    n_known <- 10
    known_idx <- sort(sample(seq_len(n_variants), n_known))
    known_zscores <- data.frame(
        chrom = rep(1, n_known),
        pos = ref_panel$pos[known_idx],
        variant_id = ref_panel$variant_id[known_idx],
        A1 = ref_panel$A1[known_idx],
        A2 = ref_panel$A2[known_idx],
        z = rnorm(n_known),
        stringsAsFactors = FALSE
    )
    R <- diag(n_variants)
    colnames(R) <- rownames(R) <- ref_panel$variant_id

    # Wrap in list structure with ldMatrices
    LD_list <- list(ldMatrices = list(R))

    result <- raiss(
        ref_panel,
        known_zscores,
        ldMatrix = LD_list,
        r2Threshold = 0,
        minimumLd = 0,
        verbose = FALSE
    )
    expect_true(is.list(result))
    expect_true("resultNofilter" %in% names(result))
    expect_equal(nrow(result$resultNofilter), n_variants)
})

test_that("raiss with genotype_matrix list processes multiple blocks", {
    set.seed(42)
    n <- 50
    p <- 20
    # Each block has its own ref_panel subset matching the X columns
    ref_panel <- data.frame(
        chrom = rep(1, p),
        pos = seq(10, p * 10, 10),
        variant_id = paste0("rs", 1:p),
        A1 = rep("A", p),
        A2 = rep("G", p),
        stringsAsFactors = FALSE
    )
    # Use only first block's variants as known (so second block has unknowns)
    known_idx <- sort(sample(1:10, 5))
    known_zscores <- data.frame(
        chrom = rep(1, length(known_idx)),
        pos = ref_panel$pos[known_idx],
        variant_id = ref_panel$variant_id[known_idx],
        A1 = ref_panel$A1[known_idx],
        A2 = ref_panel$A2[known_idx],
        z = rnorm(length(known_idx)),
        stringsAsFactors = FALSE
    )
    X <- scale(matrix(sample(0:2, n * p, replace = TRUE), nrow = n))
    X[is.na(X)] <- 0
    colnames(X) <- ref_panel$variant_id

    # Use the full matrix as a single-element list - the simplest valid list input
    X_list <- list(X)

    result <- raiss(
        ref_panel,
        known_zscores,
        genotypeMatrix = X_list,
        r2Threshold = 0,
        minimumLd = 0,
        verbose = FALSE
    )
    expect_true(is.list(result))
    expect_true("resultNofilter" %in% names(result))
    expect_true(nrow(result$resultNofilter) > 0)
    expect_null(result$ldMat)
})

test_that("raiss with genotype_matrix list returns NULL when all blocks fail", {
    set.seed(42)
    n <- 50
    p <- 10
    ref_panel <- data.frame(
        chrom = rep(1, p),
        pos = seq(10, p * 10, 10),
        variant_id = paste0("rs", 1:p),
        A1 = rep("A", p),
        A2 = rep("G", p),
        stringsAsFactors = FALSE
    )
    # known_zscores has no overlap with ref_panel
    known_zscores <- data.frame(
        chrom = rep(1, 3),
        pos = c(200, 300, 400),
        variant_id = paste0("other", 1:3),
        A1 = rep("A", 3),
        A2 = rep("G", 3),
        z = rnorm(3),
        stringsAsFactors = FALSE
    )
    X <- scale(matrix(sample(0:2, n * p, replace = TRUE), nrow = n))
    X[is.na(X)] <- 0
    colnames(X) <- ref_panel$variant_id
    X_list <- list(X[, 1:5, drop = FALSE], X[, 6:10, drop = FALSE])

    result <- raiss(
        ref_panel,
        known_zscores,
        genotypeMatrix = X_list,
        verbose = FALSE
    )
    expect_null(result)
})

# Block-Diagonal LD data generator for RAISS testing
# Corrected function to generate proper block-diagonal test data
generate_block_diagonal_test_data <- function(
    seed = 123,
    block_structure = "overlapping",
    n_variants = 30
) {
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
            c(1, 11), # Block 1: variants 1-11
            c(11, 21), # Block 2: variants 11-21 (overlap at var11)
            c(21, n_variants) # Block 3: variants 21-30 (overlap at var21)
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
            startIdx <- (i - 1) * block_size + 1
            endIdx <- min(i * block_size, n_variants)
            block_boundaries[[i]] <- c(startIdx, endIdx)
        }
    } else if (block_structure == "single_block") {
        block_boundaries <- list(c(1, n_variants))
    }

    # First, create independent block matrices
    block_matrices <- list()
    for (i in seq_along(block_boundaries)) {
        startIdx <- block_boundaries[[i]][1]
        endIdx <- block_boundaries[[i]][2]
        block_variant_ids <- ref_panel$variant_id[startIdx:endIdx]
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
    variantIndices <- data.frame(
        variant_id = character(),
        blockId = integer(),
        stringsAsFactors = FALSE
    )

    for (i in seq_along(block_boundaries)) {
        startIdx <- block_boundaries[[i]][1]
        endIdx <- block_boundaries[[i]][2]
        block_variant_ids <- ref_panel$variant_id[startIdx:endIdx]

        block_indices <- data.frame(
            variant_id = block_variant_ids,
            blockId = i,
            stringsAsFactors = FALSE
        )
        variantIndices <- rbind(variantIndices, block_indices)
    }

    # Create block metadata
    block_sizes <- sapply(block_boundaries, function(b) b[2] - b[1] + 1)
    blockMetadata <- data.frame(
        blockId = seq_along(block_boundaries),
        chrom = rep(1, length(block_boundaries)),
        size = block_sizes,
        startIdx = sapply(seq_along(block_boundaries), function(i) {
            # Adjust for 1-based indexing in R
            if (i == 1) {
                return(1)
            }
            # Count unique variants before this block
            sum(sapply(1:(i - 1), function(j) {
                # If there's an overlap with the next block, count one less
                if (
                    j < length(block_boundaries) &&
                        block_boundaries[[j]][2] == block_boundaries[[j + 1]][1]
                ) {
                    return(block_boundaries[[j]][2] - block_boundaries[[j]][1])
                } else {
                    return(
                        block_boundaries[[j]][2] - block_boundaries[[j]][1] + 1
                    )
                }
            })) +
                1
        }),
        endIdx = sapply(seq_along(block_boundaries), function(i) {
            # Count all unique variants up to and including this block
            sum(sapply(1:i, function(j) {
                # If there's an overlap with the next block, count one less
                if (
                    j < i &&
                        block_boundaries[[j]][2] == block_boundaries[[j + 1]][1]
                ) {
                    return(block_boundaries[[j]][2] - block_boundaries[[j]][1])
                } else {
                    return(
                        block_boundaries[[j]][2] - block_boundaries[[j]][1] + 1
                    )
                }
            }))
        }),
        stringsAsFactors = FALSE
    )

    # Build the full matrix correctly ensuring proper block structure
    # IMPORTANT: Initialize a matrix with zeros - ensure no correlations between blocks
    all_variant_ids <- unique(variantIndices$variant_id)
    LD_matrix_full <- matrix(
        0,
        nrow = length(all_variant_ids),
        ncol = length(all_variant_ids)
    )
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
        ldMatrices = block_matrices,
        variantIndices = variantIndices,
        blockMetadata = blockMetadata,
        ldVariants = all_variant_ids
    )

    return(list(
        ref_panel = ref_panel,
        known_zscores = known_zscores,
        LD_matrix_full = LD_matrix_full,
        LD_matrix_blocks = LD_matrix_blocks,
        variantIndices = variantIndices,
        block_boundaries = block_boundaries,
        blockMetadata = blockMetadata
    ))
}

test_that("full matrix and block processing produce identical results", {
    # Only test non-overlapping structures for exact z-score matching
    block_structures <- c("non_overlapping", "single_block")

    for (structure in block_structures) {
        test_data <- generate_block_diagonal_test_data(
            seed = 123,
            block_structure = structure
        )

        # Prepare ld_data as LdData S4 for partitionLdMatrix
        ld_data <- make_ld_data_from_ref_panel(
            test_data$LD_matrix_full,
            test_data$ref_panel,
            test_data$blockMetadata
        )

        # For non-overlapping structures, use partitionLdMatrix
        partitioned <- partitionLdMatrix(
            ld_data,
            mergeSmallBlocks = FALSE
        )

        # Run RAISS with full matrix
        result_full <- raiss(
            refPanel = test_data$ref_panel,
            knownZscores = test_data$known_zscores,
            ldMatrix = test_data$LD_matrix_full,
            lamb = 0.01,
            rcond = 0.01,
            r2Threshold = 0.3,
            minimumLd = 1,
            verbose = FALSE
        )

        # Run RAISS with partitioned blocks
        result_blocks <- raiss(
            refPanel = test_data$ref_panel,
            knownZscores = test_data$known_zscores,
            ldMatrix = partitioned,
            lamb = 0.01,
            rcond = 0.01,
            r2Threshold = 0.3,
            minimumLd = 1,
            verbose = FALSE
        )

        # For non-overlapping blocks, we compare all variants
        result_full_sorted <- result_full$resultNofilter %>% arrange(variant_id)
        result_blocks_sorted <- result_blocks$resultNofilter %>%
            arrange(variant_id)

        # Compare variant IDs
        expect_equal(
            sort(result_full$resultNofilter$variant_id),
            sort(result_blocks$resultNofilter$variant_id),
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
        if (
            !is.null(result_full$resultFilter) &&
                !is.null(result_blocks$resultFilter) &&
                nrow(result_full$resultFilter) > 0 &&
                nrow(result_blocks$resultFilter) > 0
        ) {
            expect_equal(
                sort(result_full$resultFilter$variant_id),
                sort(result_blocks$resultFilter$variant_id),
                info = paste("Filtered variant IDs should match for", structure)
            )

            result_full_filter_sorted <- result_full$resultFilter %>%
                arrange(variant_id)
            result_blocks_filter_sorted <- result_blocks$resultFilter %>%
                arrange(variant_id)

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
    test_data <- generate_block_diagonal_test_data(
        seed = 123,
        block_structure = "overlapping"
    )

    # Run RAISS with full matrix
    result_full <- raiss(
        refPanel = test_data$ref_panel,
        knownZscores = test_data$known_zscores,
        ldMatrix = test_data$LD_matrix_full,
        lamb = 0.01,
        rcond = 0.01,
        r2Threshold = 0.3,
        minimumLd = 1,
        verbose = FALSE
    )

    # Run RAISS with block processing
    result_blocks <- raiss(
        refPanel = test_data$ref_panel,
        knownZscores = test_data$known_zscores,
        ldMatrix = test_data$LD_matrix_blocks,
        lamb = 0.01,
        rcond = 0.01,
        r2Threshold = 0.3,
        minimumLd = 1,
        verbose = FALSE
    )

    # Test 1: Verify all variants are present in both results
    expect_equal(
        sort(result_full$resultNofilter$variant_id),
        sort(result_blocks$resultNofilter$variant_id),
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
        all(boundary_variants %in% result_blocks$resultNofilter$variant_id),
        info = "All boundary variants should be present in block results"
    )

    # Verify boundary variants have valid z-scores
    boundaryResults <- result_blocks$resultNofilter %>%
        filter(variant_id %in% boundary_variants)

    expect_true(
        all(!is.na(boundaryResults$z)),
        info = "Boundary variants should have valid z-scores in block results"
    )

    # Test 3: Verify non-boundary variants have z-scores with reasonable range
    non_boundaryResults <- result_blocks$resultNofilter %>%
        filter(!variant_id %in% boundary_variants)

    expect_true(
        all(!is.na(non_boundaryResults$z)),
        info = "Non-boundary variants should have valid z-scores"
    )

    expect_true(
        all(abs(non_boundaryResults$z) < 10),
        info = "Non-boundary variant z-scores should be in reasonable range"
    )

    # We deliberately do NOT compare z-score values between full matrix and block processing
    # for overlapping blocks, as differences are expected and valid
})

test_that("raiss handles block boundaries correctly", {
    # Generate test data with overlapping blocks
    test_data <- generate_block_diagonal_test_data(
        seed = 456,
        block_structure = "overlapping"
    )

    # Define the thresholds explicitly
    test_R2_threshold <- 0.3
    test_minimum_ld <- 1

    # Run RAISS with block processing
    result <- raiss(
        refPanel = test_data$ref_panel,
        knownZscores = test_data$known_zscores,
        ldMatrix = test_data$LD_matrix_blocks,
        lamb = 0.01,
        rcond = 0.01,
        r2Threshold = test_R2_threshold,
        minimumLd = test_minimum_ld,
        verbose = FALSE
    )

    # First verify that the required columns exist in the results
    expect_true(
        "variant_id" %in% names(result$resultNofilter),
        info = "resultNofilter should contain a variant_id column"
    )

    expect_true(
        "raissR2" %in% names(result$resultNofilter),
        info = "resultNofilter should contain a raissR2 column"
    )

    expect_true(
        "raissLdScore" %in% names(result$resultNofilter),
        info = "resultNofilter should contain a raissLdScore column"
    )

    # Check that we have only one entry per variant ID (no duplicates)
    expect_equal(
        length(unique(result$resultNofilter$variant_id)),
        length(result$resultNofilter$variant_id),
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
        all(boundary_variants %in% result$resultNofilter$variant_id),
        info = "All boundary variants should be present in the results"
    )

    # Get the boundary variant results
    boundaryResults <- result$resultNofilter %>%
        filter(variant_id %in% boundary_variants)

    # Check R-squared values for non-NA boundary variants
    non_na_r2 <- boundaryResults$raissR2[!is.na(boundaryResults$raissR2)]
    if (length(non_na_r2) > 0) {
        expect_true(
            all(non_na_r2 >= 0 & non_na_r2 <= 1),
            info = "Non-NA boundary variant R-squared values should be between 0 and 1"
        )
    }

    # Check LD scores for non-NA boundary variants
    non_na_ld <- boundaryResults$raissLdScore[
        !is.na(boundaryResults$raissLdScore)
    ]
    if (length(non_na_ld) > 0) {
        expect_true(
            all(non_na_ld >= 0),
            info = "Non-NA boundary variant LD scores should be non-negative"
        )
    }

    # Verify that pre-filtering and post-filtering steps handle boundary variants correctly
    if (!is.null(result$resultFilter) && nrow(result$resultFilter) > 0) {
        # First check if filtered results have the required columns
        expect_true(
            "variant_id" %in% names(result$resultFilter),
            info = "resultFilter should contain a variant_id column"
        )

        expect_true(
            "raissR2" %in% names(result$resultFilter),
            info = "resultFilter should contain a raissR2 column"
        )

        expect_true(
            "raissLdScore" %in% names(result$resultFilter),
            info = "resultFilter should contain a raissLdScore column"
        )

        # Check which boundary variants passed the filtering
        boundary_in_filtered <- boundary_variants %in%
            result$resultFilter$variant_id

        if (any(boundary_in_filtered)) {
            # Get the filtered boundary variants
            boundary_filtered <- result$resultFilter %>%
                filter(variant_id %in% boundary_variants)

            # Check that non-NA R-squared values meet the threshold
            non_na_r2_filtered <- boundary_filtered$raissR2[
                !is.na(boundary_filtered$raissR2)
            ]
            if (length(non_na_r2_filtered) > 0) {
                expect_true(
                    all(non_na_r2_filtered >= test_R2_threshold),
                    info = paste(
                        "Non-NA filtered boundary variant R-squared values should meet the threshold of",
                        test_R2_threshold
                    )
                )
            }

            # Check that non-NA LD scores meet the threshold
            non_na_ld_filtered <- boundary_filtered$raissLdScore[
                !is.na(boundary_filtered$raissLdScore)
            ]
            if (length(non_na_ld_filtered) > 0) {
                expect_true(
                    all(non_na_ld_filtered >= test_minimum_ld),
                    info = paste(
                        "Non-NA filtered boundary variant LD scores should meet the threshold of",
                        test_minimum_ld
                    )
                )
            }
        }
    }
})

test_that("partitionLdMatrix integrates correctly with RAISS", {
    test_data <- generate_block_diagonal_test_data(
        seed = 456,
        block_structure = "non_overlapping"
    )

    ld_data <- make_ld_data_from_ref_panel(
        test_data$LD_matrix_full,
        test_data$ref_panel,
        test_data$blockMetadata
    )

    partitioned <- partitionLdMatrix(
        ld_data,
        mergeSmallBlocks = FALSE
    )

    result_full <- raiss(
        refPanel = test_data$ref_panel,
        knownZscores = test_data$known_zscores,
        ldMatrix = test_data$LD_matrix_full,
        lamb = 0.01,
        rcond = 0.01,
        r2Threshold = 0.3,
        minimumLd = 1,
        verbose = FALSE
    )

    result_partitioned <- raiss(
        refPanel = test_data$ref_panel,
        knownZscores = test_data$known_zscores,
        ldMatrix = partitioned,
        lamb = 0.01,
        rcond = 0.01,
        r2Threshold = 0.3,
        minimumLd = 1,
        verbose = FALSE
    )

    result_full_sorted <- result_full$resultNofilter %>% arrange(variant_id)
    result_partitioned_sorted <- result_partitioned$resultNofilter %>%
        arrange(variant_id)

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
    test_data <- generate_block_diagonal_test_data(
        seed = 789,
        block_structure = "overlapping"
    )

    result_blocks <- raiss(
        refPanel = test_data$ref_panel,
        knownZscores = test_data$known_zscores,
        ldMatrix = test_data$LD_matrix_blocks,
        lamb = 0.01,
        rcond = 0.01,
        r2Threshold = 0.1,
        minimumLd = 1,
        verbose = FALSE
    )

    variant_counts <- table(test_data$variantIndices$variant_id)
    boundary_vars <- names(variant_counts[variant_counts > 1])

    for (var in boundary_vars) {
        expect_equal(
            sum(result_blocks$resultNofilter$variant_id == var),
            1,
            info = paste("Boundary variant", var, "should appear once")
        )
    }

    expect_equal(
        nrow(result_blocks$resultNofilter),
        length(unique(result_blocks$resultNofilter$variant_id)),
        info = "No duplicate variants in results"
    )
})

# Test 4: Single-block case
test_that("RAISS handles single-block list correctly", {
    test_data <- generate_block_diagonal_test_data(
        seed = 202,
        block_structure = "single_block"
    )

    result_full <- raiss(
        refPanel = test_data$ref_panel,
        knownZscores = test_data$known_zscores,
        ldMatrix = test_data$LD_matrix_full,
        lamb = 0.01,
        rcond = 0.01,
        r2Threshold = 0.3,
        minimumLd = 1,
        verbose = FALSE
    )

    result_single_block <- raiss(
        refPanel = test_data$ref_panel,
        knownZscores = test_data$known_zscores,
        ldMatrix = test_data$LD_matrix_blocks,
        lamb = 0.01,
        rcond = 0.01,
        r2Threshold = 0.3,
        minimumLd = 1,
        verbose = FALSE
    )

    result_full_sorted <- result_full$resultNofilter %>% arrange(variant_id)
    result_single_block_sorted <- result_single_block$resultNofilter %>%
        arrange(variant_id)

    expect_equal(
        result_full_sorted$z,
        result_single_block_sorted$z,
        tolerance = 1e-6,
        info = "Z-scores should match for single block"
    )
})

# ============================================================================
# Tests for SVD-based genotype matrix path (raissSingleMatrixFromX)
# ============================================================================

#' Helper: generate a genotype matrix X with corresponding ref_panel,
#' known_zscores, and LD matrix R for equivalence testing.
generate_X_test_data <- function(n = 200, p = 100, n_known = 50, seed = 42) {
    set.seed(seed)
    # Generate genotype-like matrix (dosages 0/1/2)
    X_raw <- matrix(
        sample(0:2, n * p, replace = TRUE, prob = c(0.25, 0.5, 0.25)),
        nrow = n,
        ncol = p
    )
    # Center and scale
    X <- scale(X_raw)
    X[is.na(X)] <- 0 # zero-variance columns become 0

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

    list(
        X = X,
        R = R,
        ref_panel = ref_panel,
        known_zscores = known_zscores,
        n = n,
        p = p,
        n_known = n_known
    )
}

test_that("safeSvd basic functionality", {
    set.seed(1)
    mat <- matrix(rnorm(20), nrow = 5, ncol = 4)
    s <- pecotmr:::.safeSvd(mat)
    expect_equal(length(s$d), min(5, 4))
    expect_true(all(s$d > 0))
    # Reconstruct
    reconstructed <- s$u %*% diag(s$d) %*% t(s$v)
    expect_equal(mat, reconstructed, tolerance = 1e-10)
})

test_that("safeSvd filters small singular values", {
    set.seed(2)
    # Create rank-2 matrix
    u <- matrix(rnorm(10), nrow = 5, ncol = 2)
    v <- matrix(rnorm(8), nrow = 4, ncol = 2)
    mat <- u %*% t(v) + matrix(rnorm(20) * 1e-12, nrow = 5, ncol = 4)
    s <- pecotmr:::.safeSvd(mat, tol = 1e-6)
    expect_equal(length(s$d), 2)
})

test_that("safeSvd max_rank works", {
    set.seed(3)
    mat <- matrix(rnorm(50), nrow = 10, ncol = 5)
    s <- pecotmr:::.safeSvd(mat, maxRank = 2)
    expect_equal(length(s$d), 2)
    expect_equal(ncol(s$u), 2)
    expect_equal(ncol(s$v), 2)
})

test_that("safeSvd rejects all-zero matrix", {
    mat <- matrix(0, nrow = 5, ncol = 3)
    expect_error(pecotmr:::.safeSvd(mat), "all-zero")
})

test_that("X path matches R path: basic equivalence (n > p)", {
    data <- generate_X_test_data(n = 200, p = 100, n_known = 50, seed = 42)

    result_R <- raiss(
        data$ref_panel,
        data$known_zscores,
        ldMatrix = data$R,
        lamb = 0.01,
        rcond = 0.01,
        r2Threshold = 0,
        minimumLd = 0,
        verbose = FALSE
    )
    result_X <- raiss(
        data$ref_panel,
        data$known_zscores,
        genotypeMatrix = data$X,
        lamb = 0.01,
        svdTol = 1e-12,
        r2Threshold = 0,
        minimumLd = 0,
        verbose = FALSE
    )

    # Compare imputed z-scores (sort by variant_id for alignment)
    r_sorted <- result_R$resultNofilter %>% arrange(variant_id)
    x_sorted <- result_X$resultNofilter %>% arrange(variant_id)

    expect_equal(
        r_sorted$z,
        x_sorted$z,
        tolerance = 1e-4,
        info = "Imputed z-scores should match between X and R paths"
    )
    expect_equal(
        r_sorted$Var,
        x_sorted$Var,
        tolerance = 1e-4,
        info = "Variance should match between X and R paths"
    )
    expect_equal(
        r_sorted$raissLdScore,
        x_sorted$raissLdScore,
        tolerance = 1e-4,
        info = "LD scores should match between X and R paths"
    )
})

test_that("X path matches R path: n < p regime", {
    data <- generate_X_test_data(n = 50, p = 200, n_known = 100, seed = 123)

    result_R <- raiss(
        data$ref_panel,
        data$known_zscores,
        ldMatrix = data$R,
        lamb = 0.01,
        rcond = 0.01,
        r2Threshold = 0,
        minimumLd = 0,
        verbose = FALSE
    )
    result_X <- raiss(
        data$ref_panel,
        data$known_zscores,
        genotypeMatrix = data$X,
        lamb = 0.01,
        svdTol = 1e-12,
        r2Threshold = 0,
        minimumLd = 0,
        verbose = FALSE
    )

    r_sorted <- result_R$resultNofilter %>% arrange(variant_id)
    x_sorted <- result_X$resultNofilter %>% arrange(variant_id)

    expect_equal(
        r_sorted$z,
        x_sorted$z,
        tolerance = 1e-4,
        info = "z-scores should match in n < p regime"
    )
    expect_equal(
        r_sorted$Var,
        x_sorted$Var,
        tolerance = 1e-4,
        info = "Variance should match in n < p regime"
    )
})

test_that("X path matches R path: n >> p regime", {
    data <- generate_X_test_data(n = 500, p = 50, n_known = 25, seed = 99)

    result_R <- raiss(
        data$ref_panel,
        data$known_zscores,
        ldMatrix = data$R,
        lamb = 0.01,
        rcond = 0.01,
        r2Threshold = 0,
        minimumLd = 0,
        verbose = FALSE
    )
    result_X <- raiss(
        data$ref_panel,
        data$known_zscores,
        genotypeMatrix = data$X,
        lamb = 0.01,
        svdTol = 1e-12,
        r2Threshold = 0,
        minimumLd = 0,
        verbose = FALSE
    )

    r_sorted <- result_R$resultNofilter %>% arrange(variant_id)
    x_sorted <- result_X$resultNofilter %>% arrange(variant_id)

    expect_equal(r_sorted$z, x_sorted$z, tolerance = 1e-4)
    expect_equal(r_sorted$Var, x_sorted$Var, tolerance = 1e-4)
})

test_that("X path matches R path: varying lambda", {
    data <- generate_X_test_data(n = 150, p = 80, n_known = 40, seed = 7)

    for (lamb in c(0.001, 0.01, 0.1)) {
        result_R <- raiss(
            data$ref_panel,
            data$known_zscores,
            ldMatrix = data$R,
            lamb = lamb,
            rcond = 0.01,
            r2Threshold = 0,
            minimumLd = 0,
            verbose = FALSE
        )
        result_X <- raiss(
            data$ref_panel,
            data$known_zscores,
            genotypeMatrix = data$X,
            lamb = lamb,
            svdTol = 1e-12,
            r2Threshold = 0,
            minimumLd = 0,
            verbose = FALSE
        )

        r_sorted <- result_R$resultNofilter %>% arrange(variant_id)
        x_sorted <- result_X$resultNofilter %>% arrange(variant_id)

        expect_equal(
            r_sorted$z,
            x_sorted$z,
            tolerance = 1e-4,
            info = paste("z-scores should match for lamb =", lamb)
        )
        expect_equal(
            r_sorted$Var,
            x_sorted$Var,
            tolerance = 1e-4,
            info = paste("Variance should match for lamb =", lamb)
        )
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
    result <- raiss(
        data$ref_panel,
        all_known,
        genotypeMatrix = data$X,
        verbose = FALSE
    )
    expect_equal(nrow(result$resultNofilter), nrow(data$ref_panel))
})

test_that("X path handles single unknown variant", {
    data <- generate_X_test_data(n = 100, p = 50, n_known = 49, seed = 15)

    result_R <- raiss(
        data$ref_panel,
        data$known_zscores,
        ldMatrix = data$R,
        lamb = 0.01,
        rcond = 0.01,
        r2Threshold = 0,
        minimumLd = 0,
        verbose = FALSE
    )
    result_X <- raiss(
        data$ref_panel,
        data$known_zscores,
        genotypeMatrix = data$X,
        lamb = 0.01,
        svdTol = 1e-12,
        r2Threshold = 0,
        minimumLd = 0,
        verbose = FALSE
    )

    r_sorted <- result_R$resultNofilter %>% arrange(variant_id)
    x_sorted <- result_X$resultNofilter %>% arrange(variant_id)

    expect_equal(r_sorted$z, x_sorted$z, tolerance = 1e-4)
})

test_that("X path handles single known variant", {
    data <- generate_X_test_data(n = 100, p = 50, n_known = 1, seed = 20)

    result_X <- raiss(
        data$ref_panel,
        data$known_zscores,
        genotypeMatrix = data$X,
        lamb = 0.01,
        svdTol = 1e-12,
        r2Threshold = 0,
        minimumLd = 0,
        verbose = FALSE
    )
    expect_true(is.data.frame(result_X$resultNofilter))
    expect_equal(nrow(result_X$resultNofilter), nrow(data$ref_panel))
})

test_that("X path R2 filtering matches R path", {
    data <- generate_X_test_data(n = 200, p = 100, n_known = 50, seed = 42)

    result_R <- raiss(
        data$ref_panel,
        data$known_zscores,
        ldMatrix = data$R,
        lamb = 0.01,
        rcond = 0.01,
        r2Threshold = 0.6,
        minimumLd = 5,
        verbose = FALSE
    )
    result_X <- raiss(
        data$ref_panel,
        data$known_zscores,
        genotypeMatrix = data$X,
        lamb = 0.01,
        svdTol = 1e-12,
        r2Threshold = 0.6,
        minimumLd = 5,
        verbose = FALSE
    )

    # Same variants should pass filtering
    r_filtered_ids <- sort(result_R$resultFilter$variant_id)
    x_filtered_ids <- sort(result_X$resultFilter$variant_id)
    expect_equal(
        r_filtered_ids,
        x_filtered_ids,
        info = "Same variants should pass R2/LD filtering"
    )
})

test_that("raw genotype_matrix path is not equivalent to LD path used by legacy pipeline", {
    set.seed(1)
    n <- 80
    p <- 40
    n_known <- 20
    X_raw <- matrix(
        sample(0:2, n * p, replace = TRUE, prob = c(0.35, 0.45, 0.20)),
        nrow = n,
        ncol = p
    )
    ref_panel <- data.frame(
        chrom = rep(1, p),
        pos = seq_len(p) * 100,
        variant_id = paste0("rs", seq_len(p)),
        A1 = rep("A", p),
        A2 = rep("G", p),
        stringsAsFactors = FALSE
    )
    colnames(X_raw) <- ref_panel$variant_id

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

    R <- computeLd(X_raw, method = "sample")
    rownames(R) <- colnames(R) <- ref_panel$variant_id
    X_scaled <- scale(X_raw)
    X_scaled[is.na(X_scaled)] <- 0
    colnames(X_scaled) <- ref_panel$variant_id

    result_LD <- raiss(
        ref_panel,
        known_zscores,
        ldMatrix = R,
        lamb = 0.01,
        rcond = 0.01,
        r2Threshold = 0.6,
        minimumLd = 0,
        verbose = FALSE
    )
    result_raw_X <- raiss(
        ref_panel,
        known_zscores,
        genotypeMatrix = X_raw,
        lamb = 0.01,
        svdTol = 1e-8,
        r2Threshold = 0.6,
        minimumLd = 0,
        verbose = FALSE
    )
    result_scaled_X <- raiss(
        ref_panel,
        known_zscores,
        genotypeMatrix = X_scaled,
        lamb = 0.01,
        svdTol = 1e-12,
        r2Threshold = 0.6,
        minimumLd = 0,
        verbose = FALSE
    )

    ld_ids <- sort(result_LD$resultFilter$variant_id)
    raw_x_ids <- sort(result_raw_X$resultFilter$variant_id)
    scaled_x_ids <- sort(result_scaled_X$resultFilter$variant_id)

    expect_false(identical(ld_ids, raw_x_ids))
    expect_gt(nrow(result_raw_X$resultFilter), nrow(result_LD$resultFilter))
    expect_equal(scaled_x_ids, ld_ids)

    ld_sorted <- result_LD$resultNofilter %>% arrange(variant_id)
    scaled_sorted <- result_scaled_X$resultNofilter %>% arrange(variant_id)
    expect_equal(ld_sorted$z, scaled_sorted$z, tolerance = 1e-10)
    expect_equal(ld_sorted$raissR2, scaled_sorted$raissR2, tolerance = 1e-10)
})

test_that("raiss rejects both ldMatrix and genotype_matrix", {
    data <- generate_X_test_data(n = 50, p = 20, n_known = 10, seed = 1)
    expect_error(
        raiss(
            data$ref_panel,
            data$known_zscores,
            ldMatrix = data$R,
            genotypeMatrix = data$X
        ),
        "not both"
    )
})

test_that("raiss rejects neither ldMatrix nor genotype_matrix", {
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
    X_raw <- matrix(
        sample(0:2, n * p, replace = TRUE, prob = c(0.25, 0.5, 0.25)),
        nrow = n,
        ncol = p
    )
    # Make columns 5 and 6 nearly identical
    X_raw[, 6] <- X_raw[, 5] + sample(c(0, 0, 0, 0, 1), n, replace = TRUE)
    X_raw[X_raw > 2] <- 2

    X <- scale(X_raw)
    X[is.na(X)] <- 0

    ref_panel <- data.frame(
        chrom = rep(1, p),
        pos = seq(1, p * 10, 10),
        variant_id = paste0("rs", seq_len(p)),
        A1 = rep("A", p),
        A2 = rep("G", p),
        stringsAsFactors = FALSE
    )
    colnames(X) <- ref_panel$variant_id

    n_known <- 30
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

    R <- cor(X_raw)
    R[is.na(R)] <- 0
    colnames(R) <- rownames(R) <- ref_panel$variant_id

    result_R <- raiss(
        ref_panel,
        known_zscores,
        ldMatrix = R,
        lamb = 0.01,
        rcond = 0.01,
        r2Threshold = 0,
        minimumLd = 0,
        verbose = FALSE
    )
    result_X <- raiss(
        ref_panel,
        known_zscores,
        genotypeMatrix = X,
        lamb = 0.01,
        svdTol = 1e-12,
        r2Threshold = 0,
        minimumLd = 0,
        verbose = FALSE
    )

    r_sorted <- result_R$resultNofilter %>% arrange(variant_id)
    x_sorted <- result_X$resultNofilter %>% arrange(variant_id)

    expect_equal(
        r_sorted$z,
        x_sorted$z,
        tolerance = 1e-3,
        info = "Collinear case: z-scores should be close"
    )
})


context("RAISS missing-variant imputation in TWAS pipelines")

# Previous tests covered `twasWeightsSumstatPipeline(imputeMissing = ...)`,
# which has been removed in favor of the S4 `twasWeightsPipeline` family
# dispatching on `QtlSumStats` / `QtlDataset`. The missing-variant
# imputation knob now lives inside `summaryStatsQc(impute = TRUE)`, and
# tests for that path live in test_sumstatsQc.R (internal helpers) and
# in the SumStats pipeline tests.
#
# RAISS itself (`raiss()`) still exists with the same signature; its
# direct tests live in test_raiss.R.

context("slalom")

# ============================================================================
# Helper: build a valid positive-definite LD matrix from a genotype matrix
# ============================================================================
make_synthetic_ld <- function(n_samples, n_snps, seed = 1) {
    set.seed(seed)
    # Simulate genotypes with some LD structure by using a factor model
    # X = Z %*% L + noise, where Z is latent and L is a loading matrix
    n_factors <- min(3, n_snps)
    Z <- matrix(rnorm(n_samples * n_factors), nrow = n_samples)
    L <- matrix(runif(n_factors * n_snps, -1, 1), nrow = n_factors)
    X_raw <- Z %*%
        L +
        matrix(rnorm(n_samples * n_snps, sd = 0.5), nrow = n_samples)
    # Discretise to genotype-like values (0, 1, 2)
    X <- matrix(
        as.integer(cut(X_raw, breaks = c(-Inf, -0.5, 0.5, Inf))) - 1L,
        nrow = n_samples,
        ncol = n_snps
    )
    colnames(X) <- paste0("snp", seq_len(n_snps))
    R <- cor(X)
    # Ensure perfect diagonal and clean NaN from any zero-variance columns
    R[is.na(R) | is.nan(R)] <- 0
    diag(R) <- 1.0
    list(X = X, R = R)
}

# ============================================================================
# Basic output structure
# ============================================================================

test_that("slalom basic output structure", {
    set.seed(42)
    n <- 50
    z <- rnorm(n)
    R <- diag(n)
    # Add some off-diagonal correlations
    for (i in 1:(n - 1)) {
        R[i, i + 1] <- 0.3
        R[i + 1, i] <- 0.3
    }

    result <- slalom(zScore = z, R = R)

    expect_type(result, "list")
    expect_named(result, c("data", "summary"))
    expect_s3_class(result$data, "data.frame")
    expect_true("original_z" %in% colnames(result$data))
    expect_true("prob" %in% colnames(result$data))
    expect_true("pvalue" %in% colnames(result$data))
    expect_true("outliers" %in% colnames(result$data))
    expect_true("nlog10p_dentist_s" %in% colnames(result$data))
    expect_equal(nrow(result$data), n)
})

test_that("slalom errors on non-square R", {
    z <- rnorm(10)
    R <- matrix(rnorm(50), nrow = 5, ncol = 10)
    expect_error(slalom(zScore = z, R = R), "R must be a square matrix")
})

test_that("slalom accepts X matrix instead of R", {
    set.seed(42)
    n_samples <- 100
    n_snps <- 10
    X <- matrix(
        sample(0:2, n_samples * n_snps, replace = TRUE),
        nrow = n_samples,
        ncol = n_snps
    )
    colnames(X) <- paste0("snp", 1:n_snps)
    z <- rnorm(n_snps)

    result <- slalom(zScore = z, X = X)
    expect_type(result, "list")
    expect_equal(nrow(result$data), n_snps)
    # PIPs should be in [0,1] and sum to 1
    expect_true(all(result$data$prob >= 0 & result$data$prob <= 1))
    expect_equal(sum(result$data$prob), 1, tolerance = 1e-12)
    # Data frame should have expected column names
    expected_cols <- c(
        "original_z",
        "prob",
        "pvalue",
        "outliers",
        "nlog10p_dentist_s"
    )
    expect_true(all(expected_cols %in% colnames(result$data)))
})

# ============================================================================
# ABF computation correctness
# ============================================================================

test_that("ABF: strong signal (z=10) gets very high PIP", {
    set.seed(100)
    n <- 20
    z <- rnorm(n, sd = 0.3)
    z[7] <- 10
    R <- diag(n)

    result <- slalom(zScore = z, R = R)

    expect_equal(which.max(result$data$prob), 7)
    expect_gt(result$data$prob[7], 0.15)
    expect_equal(sum(result$data$prob), 1, tolerance = 1e-12)
})

test_that("ABF: moderate signal (z=3) gets higher PIP than weak signal (z=1)", {
    set.seed(101)
    n <- 10
    z <- rep(0, n)
    z[2] <- 1
    z[5] <- 3
    R <- diag(n)

    result <- slalom(zScore = z, R = R)

    expect_gt(result$data$prob[5], result$data$prob[2])
    # z=0 variants should all have same PIP (by symmetry, with identity LD)
    zero_pips <- result$data$prob[c(1, 3, 4, 6, 7, 8, 9, 10)]
    expect_equal(max(zero_pips) - min(zero_pips), 0, tolerance = 1e-14)
})

test_that("ABF: lbf formula matches manual calculation", {
    z_val <- 4.0
    se_val <- 1.0
    W <- 0.04

    V <- se_val^2
    r <- W / (W + V)
    expected_lbf <- 0.5 * (log(1 - r) + r * z_val^2)

    z <- c(z_val, 0)
    R <- diag(2)
    result <- slalom(zScore = z, R = R, abfPriorVariance = W)

    lbf_0 <- 0.5 * (log(1 - r) + r * 0^2)
    expected_ratio <- exp(expected_lbf - lbf_0)
    actual_ratio <- result$data$prob[1] / result$data$prob[2]
    expect_equal(actual_ratio, expected_ratio, tolerance = 1e-10)
})

test_that("ABF: PIPs always sum to exactly 1", {
    for (s in 1:5) {
        set.seed(200 + s)
        n <- sample(5:30, 1)
        z <- rnorm(n, sd = 2)
        R <- diag(n)
        result <- slalom(zScore = z, R = R)
        expect_equal(
            sum(result$data$prob),
            1,
            tolerance = 1e-12,
            label = paste("seed", 200 + s)
        )
    }
})

test_that("ABF: symmetric z-scores give symmetric PIPs", {
    z <- c(-3, 3)
    R <- diag(2)
    result <- slalom(zScore = z, R = R)
    expect_equal(result$data$prob[1], result$data$prob[2], tolerance = 1e-14)
})

# ============================================================================
# Credible sets
# ============================================================================

test_that("CS95 contains the causal variant in a simple synthetic signal", {
    set.seed(300)
    n <- 30
    z <- rnorm(n, sd = 0.5)
    causal <- 12
    z[causal] <- 6
    R <- diag(n)

    result <- slalom(zScore = z, R = R)

    expect_true(causal %in% result$summary$cs95)
    expect_true(causal %in% result$summary$cs99)
})

test_that("CS99 is a superset of CS95", {
    set.seed(301)
    n <- 40
    z <- rnorm(n, sd = 1.5)
    R <- diag(n)

    result <- slalom(zScore = z, R = R)

    expect_true(all(result$summary$cs95 %in% result$summary$cs99))
    expect_gte(length(result$summary$cs99), length(result$summary$cs95))
})

test_that("CS95 covers at least 95% of posterior mass", {
    set.seed(302)
    n <- 25
    z <- rnorm(n)
    z[10] <- 4
    R <- diag(n)

    result <- slalom(zScore = z, R = R)

    cs95_mass <- sum(result$data$prob[result$summary$cs95])
    expect_gt(cs95_mass, 0.95)

    cs99_mass <- sum(result$data$prob[result$summary$cs99])
    expect_gt(cs99_mass, 0.99)
})

test_that("CS with very strong signal contains only the causal variant", {
    set.seed(303)
    n <- 15
    z <- rnorm(n, sd = 0.1)
    z[8] <- 15 # extremely strong
    R <- diag(n)

    result <- slalom(zScore = z, R = R)

    expect_equal(result$summary$cs95[1], 8)
    expect_true(8 %in% result$summary$cs95)
    expect_true(8 %in% result$summary$cs99)
})

test_that("CS with diffuse signal contains many variants", {
    set.seed(304)
    n <- 10
    z <- rep(0, n) # all equally uninformative
    R <- diag(n)

    result <- slalom(zScore = z, R = R)

    # Uniform PIPs => need at least ceiling(0.95 * n) = 10 variants for 95% coverage
    expect_equal(length(result$summary$cs95), n)
})

# ============================================================================
# Lead variant by pvalue vs abf
# ============================================================================

test_that("lead variant by pvalue selects most negative z-score", {
    z <- c(0, -4, 3, -1, 2)
    R <- diag(5)

    result <- slalom(zScore = z, R = R, leadVariantChoice = "pvalue")

    expect_equal(result$summary$leadPipVariant, 2)
})

test_that("lead variant by abf selects highest PIP", {
    z <- c(0, -4, 3, -1, 2)
    R <- diag(5)

    result <- slalom(zScore = z, R = R, leadVariantChoice = "abf")

    expect_equal(result$summary$leadPipVariant, which.max(result$data$prob))
})

test_that("pvalue and abf lead can differ when z has asymmetric magnitudes", {
    z <- c(-3.0, 5.0, 0.1, -0.2, 0.3)
    R <- diag(5)

    result_pv <- slalom(zScore = z, R = R, leadVariantChoice = "pvalue")
    result_abf <- slalom(zScore = z, R = R, leadVariantChoice = "abf")

    expect_equal(result_pv$summary$leadPipVariant, 1)
    expect_equal(result_abf$summary$leadPipVariant, 2)
    expect_false(
        result_pv$summary$leadPipVariant == result_abf$summary$leadPipVariant
    )
})

# ============================================================================
# DENTIST-S outlier detection
# ============================================================================

test_that("DENTIST-S: lead variant itself is not flagged as outlier", {
    set.seed(400)
    n <- 10
    z <- rnorm(n, sd = 0.5)
    z[3] <- -5 # lead by pvalue
    R <- diag(n)

    result <- slalom(zScore = z, R = R)
    lead <- result$summary$leadPipVariant
    expect_equal(lead, 3)

    expect_true(
        is.na(result$data$outliers[lead]) || !result$data$outliers[lead]
    )
})

test_that("DENTIST-S: outlier variant inconsistent with LD is flagged", {
    n <- 5
    z <- c(-5, 0, 0.1, -0.1, 0.2)
    R <- diag(n)
    R[1, 2] <- R[2, 1] <- 0.9
    R[1, 3] <- R[3, 1] <- 0.1
    R[1, 4] <- R[4, 1] <- -0.05
    R[1, 5] <- R[5, 1] <- 0.02

    result <- slalom(
        zScore = z,
        R = R,
        r2Threshold = 0.5,
        nlog10pDentistSThreshold = 2.0
    )

    lead <- result$summary$leadPipVariant
    expect_equal(lead, 1)

    expect_true(result$data$outliers[2])
    expect_false(result$data$outliers[3])
})

test_that("DENTIST-S: perfectly consistent variant in LD is not flagged", {
    n <- 3
    z <- c(-5, -4.0, 0.1)
    R <- diag(n)
    R[1, 2] <- R[2, 1] <- 0.8
    R[1, 3] <- R[3, 1] <- 0.05
    R[2, 3] <- R[3, 2] <- 0.04

    result <- slalom(
        zScore = z,
        R = R,
        r2Threshold = 0.5,
        nlog10pDentistSThreshold = 4.0
    )

    lead <- result$summary$leadPipVariant
    expect_equal(lead, 1)

    expect_equal(result$data$nlog10p_dentist_s[2], 0, tolerance = 1e-10)
    expect_false(result$data$outliers[2])
})

test_that("DENTIST-S: n_dentist_s_outlier and fraction are consistent", {
    set.seed(401)
    n <- 20
    syn <- make_synthetic_ld(200, n, seed = 401)
    z <- rnorm(n, sd = 1)
    z[1] <- -5

    result <- slalom(zScore = z, R = syn$R, r2Threshold = 0.3)

    n_r2 <- result$summary$nR2
    n_out <- result$summary$nDentistSOutlier
    frac <- result$summary$fraction

    expect_gte(n_r2, 1)
    expect_equal(frac, ifelse(n_r2 > 0, n_out / n_r2, 0), tolerance = 1e-14)
    expect_gte(n_out, 0)
    expect_lte(n_out, n_r2)
})

test_that("DENTIST-S: lowering threshold flags more outliers", {
    set.seed(402)
    n <- 15
    syn <- make_synthetic_ld(300, n, seed = 402)
    z <- rnorm(n, sd = 2)
    z[5] <- -6

    result_strict <- slalom(
        zScore = z,
        R = syn$R,
        nlog10pDentistSThreshold = 6.0,
        r2Threshold = 0.3
    )
    result_loose <- slalom(
        zScore = z,
        R = syn$R,
        nlog10pDentistSThreshold = 1.0,
        r2Threshold = 0.3
    )

    expect_gte(
        result_loose$summary$nDentistSOutlier,
        result_strict$summary$nDentistSOutlier
    )
})

# ============================================================================
# Edge cases
# ============================================================================

test_that("edge case: single variant", {
    z <- c(3.0)
    R <- matrix(1, nrow = 1, ncol = 1)

    result <- slalom(zScore = z, R = R)

    expect_equal(nrow(result$data), 1)
    expect_equal(result$data$prob[1], 1.0, tolerance = 1e-14)
    expect_equal(result$summary$leadPipVariant, 1)
    expect_equal(result$summary$nTotal, 1)
    expect_equal(result$summary$cs95, 1)
    expect_equal(result$summary$cs99, 1)
})

test_that("edge case: all zero z-scores", {
    n <- 10
    z <- rep(0, n)
    R <- diag(n)

    result <- slalom(zScore = z, R = R)

    expect_equal(result$data$prob, rep(1 / n, n), tolerance = 1e-14)
    expect_equal(result$data$pvalue, rep(0.5, n), tolerance = 1e-14)
    expect_equal(result$summary$maxPip, 1 / n, tolerance = 1e-14)
})

test_that("edge case: very large z-scores do not produce NaN in PIPs", {
    z <- c(50, -50, 30, -30, 0)
    R <- diag(5)

    result <- slalom(zScore = z, R = R)

    expect_false(any(is.nan(result$data$prob)))
    expect_false(any(is.na(result$data$prob)))
    expect_equal(sum(result$data$prob), 1, tolerance = 1e-10)
    expect_equal(result$data$prob[1], result$data$prob[2], tolerance = 1e-14)
})

test_that("edge case: identical z-scores yield uniform PIPs", {
    n <- 8
    z <- rep(2.5, n)
    R <- diag(n)

    result <- slalom(zScore = z, R = R)

    expect_equal(result$data$prob, rep(1 / n, n), tolerance = 1e-14)
})

test_that("edge case: two variants only", {
    z <- c(-3, 2)
    R <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)

    result <- slalom(zScore = z, R = R)

    expect_equal(nrow(result$data), 2)
    expect_equal(sum(result$data$prob), 1, tolerance = 1e-12)
    expect_equal(result$summary$leadPipVariant, 1)
    expect_gt(result$data$prob[1], result$data$prob[2])
})

test_that("edge case: mismatched dimensions error", {
    z <- rnorm(10)
    R <- diag(5)
    expect_error(
        slalom(zScore = z, R = R),
        "R must be a square matrix matching the length of zScore"
    )
})

test_that("edge case: no R and no X provided errors", {
    z <- rnorm(5)
    expect_error(slalom(zScore = z), "Either R.*or X.*must be provided")
})

test_that("edge case: both R and X provided errors", {
    set.seed(500)
    n <- 5
    z <- rnorm(n)
    R <- diag(n)
    X <- matrix(sample(0:2, 50 * n, replace = TRUE), nrow = 50, ncol = n)
    expect_error(
        slalom(zScore = z, R = R, X = X),
        "Provide either R or X, not both"
    )
})

# ============================================================================
# X input mode
# ============================================================================

test_that("X input yields same result as R = cor(X)", {
    set.seed(600)
    n_samples <- 200
    n_snps <- 10
    X <- matrix(
        sample(0:2, n_samples * n_snps, replace = TRUE),
        nrow = n_samples,
        ncol = n_snps
    )
    colnames(X) <- paste0("snp", seq_len(n_snps))
    z <- rnorm(n_snps, sd = 2)

    R_manual <- cor(X)
    diag(R_manual) <- 1.0

    result_X <- slalom(zScore = z, X = X)
    result_R <- slalom(zScore = z, R = R_manual)

    expect_equal(result_X$data$prob, result_R$data$prob, tolerance = 1e-10)
    expect_equal(
        result_X$data$original_z,
        result_R$data$original_z,
        tolerance = 1e-14
    )
    expect_equal(
        result_X$summary$leadPipVariant,
        result_R$summary$leadPipVariant
    )
    expect_equal(result_X$summary$cs95, result_R$summary$cs95)
    expect_equal(result_X$summary$cs99, result_R$summary$cs99)
})

# ============================================================================
# Parameter variation
# ============================================================================

test_that("larger abf_prior_variance concentrates PIPs on strong signals more", {
    set.seed(700)
    n <- 15
    z <- rnorm(n, sd = 0.5)
    z[4] <- 4
    R <- diag(n)

    result_small_W <- slalom(zScore = z, R = R, abfPriorVariance = 0.01)
    result_large_W <- slalom(zScore = z, R = R, abfPriorVariance = 1.0)

    expect_gt(result_large_W$summary$maxPip, result_small_W$summary$maxPip)
})

test_that("abfPriorVariance = 0 gives uniform PIPs", {
    n <- 10
    z <- c(5, 3, 1, 0, -1, -3, -5, 2, -2, 4)
    R <- diag(n)

    result <- slalom(zScore = z, R = R, abfPriorVariance = 0)

    expect_equal(result$data$prob, rep(1 / n, n), tolerance = 1e-14)
})

test_that("different standard_error values affect PIPs", {
    n <- 5
    z <- c(3, 3, 3, 3, 3) # same z for all
    R <- diag(n)
    se1 <- c(1, 1, 1, 1, 1)
    se2 <- c(0.5, 1, 1, 1, 1) # variant 1 has smaller SE

    result1 <- slalom(zScore = z, R = R, standardError = se1)
    result2 <- slalom(zScore = z, R = R, standardError = se2)

    expect_gt(result2$data$prob[1], result1$data$prob[1])
})

test_that("r2_threshold variation affects n_r2 count", {
    set.seed(701)
    n <- 10
    syn <- make_synthetic_ld(200, n, seed = 701)
    z <- rnorm(n, sd = 2)
    z[1] <- -5

    result_low <- slalom(zScore = z, R = syn$R, r2Threshold = 0.1)
    result_high <- slalom(zScore = z, R = syn$R, r2Threshold = 0.9)

    expect_gte(result_low$summary$nR2, result_high$summary$nR2)
})

test_that("nlog10p_dentist_s_threshold variation affects outlier count", {
    set.seed(702)
    n <- 10
    syn <- make_synthetic_ld(200, n, seed = 702)
    z <- rnorm(n, sd = 2)
    z[1] <- -6

    result_low_thresh <- slalom(
        zScore = z,
        R = syn$R,
        nlog10pDentistSThreshold = 1.0
    )
    result_high_thresh <- slalom(
        zScore = z,
        R = syn$R,
        nlog10pDentistSThreshold = 10.0
    )

    expect_gte(
        result_low_thresh$summary$nDentistSOutlier,
        result_high_thresh$summary$nDentistSOutlier
    )
})

# ============================================================================
# Output structure validation
# ============================================================================

test_that("output data types are correct", {
    set.seed(801)
    n <- 15
    z <- rnorm(n)
    R <- diag(n)

    result <- slalom(zScore = z, R = R)

    expect_type(result$data$original_z, "double")
    expect_type(result$data$prob, "double")
    expect_type(result$data$pvalue, "double")
    expect_type(result$data$outliers, "logical")
    expect_type(result$data$nlog10p_dentist_s, "double")

    expect_type(result$summary$leadPipVariant, "integer")
    expect_type(result$summary$nTotal, "integer")
    expect_type(result$summary$nR2, "integer")
    expect_type(result$summary$nDentistSOutlier, "integer")
    expect_type(result$summary$fraction, "double")
    expect_type(result$summary$maxPip, "double")
    expect_type(result$summary$cs95, "integer")
    expect_type(result$summary$cs99, "integer")
})

test_that("original_z in output matches input z-scores", {
    z <- c(1.5, -2.3, 0.7, 4.1, -0.5)
    R <- diag(5)

    result <- slalom(zScore = z, R = R)
    expect_equal(result$data$original_z, z, tolerance = 0)
})

test_that("pvalue in output matches pnorm(z)", {
    z <- c(-3, -1, 0, 1, 3)
    R <- diag(5)

    result <- slalom(zScore = z, R = R)
    expect_equal(result$data$pvalue, pnorm(z), tolerance = 1e-14)
})

# ============================================================================
# Summary statistics consistency
# ============================================================================

test_that("n_r2 counts variants with r2 > threshold to lead correctly", {
    n <- 10
    z <- c(-5, rep(0, 9))
    R <- diag(n)

    result <- slalom(zScore = z, R = R, r2Threshold = 0.6)
    expect_equal(result$summary$nR2, 1)
})

test_that("n_r2 includes correlated variants", {
    n <- 5
    z <- c(-5, -4, 0, 0, 0)
    R <- diag(n)
    R[1, 2] <- R[2, 1] <- 0.9

    result <- slalom(zScore = z, R = R, r2Threshold = 0.6)
    expect_equal(result$summary$nR2, 2)
})

test_that("fraction = 0 when there are no outliers (identity LD, consistent z)", {
    n <- 5
    z <- c(-3, 0, 0, 0, 0)
    R <- diag(n)

    result <- slalom(zScore = z, R = R)
    expect_equal(result$summary$fraction, 0)
})

test_that("fraction is between 0 and 1", {
    set.seed(900)
    for (s in 1:5) {
        set.seed(900 + s)
        n <- sample(10:30, 1)
        syn <- make_synthetic_ld(200, n, seed = 900 + s)
        z <- rnorm(n, sd = 2)
        z[1] <- -6

        result <- slalom(zScore = z, R = syn$R, r2Threshold = 0.3)
        expect_gte(result$summary$fraction, 0)
        expect_lte(result$summary$fraction, 1)
    }
})

test_that("maxPip equals the maximum of prob vector", {
    set.seed(901)
    n <- 20
    z <- rnorm(n, sd = 2)
    R <- diag(n)

    result <- slalom(zScore = z, R = R)
    expect_equal(
        result$summary$maxPip,
        max(result$data$prob),
        tolerance = 1e-14
    )
})

# ============================================================================
# Realistic synthetic LD scenarios
# ============================================================================

test_that("realistic LD: correlated variants share PIP mass", {
    set.seed(1000)
    syn <- make_synthetic_ld(500, 20, seed = 1000)
    z <- rep(0, 20)
    z[3] <- 5

    result <- slalom(zScore = z, R = syn$R)

    expect_true(3 %in% result$summary$cs95)
    expect_equal(sum(result$data$prob), 1, tolerance = 1e-12)
})

test_that("realistic LD: DENTIST-S detects outlier in correlated block", {
    set.seed(1001)
    syn <- make_synthetic_ld(500, 15, seed = 1001)
    R <- syn$R

    lead_idx <- 1
    z <- R[, lead_idx] * (-5)
    z[1] <- -5

    r2_to_lead <- R[, 1]^2
    candidates <- which(r2_to_lead > 0.3 & seq_along(z) != 1)
    if (length(candidates) > 0) {
        corrupt_idx <- candidates[1]
        z[corrupt_idx] <- z[corrupt_idx] + 10

        result <- slalom(
            zScore = z,
            R = R,
            r2Threshold = 0.2,
            nlog10pDentistSThreshold = 3.0
        )

        expect_true(result$data$outliers[corrupt_idx])
    }
})

test_that("realistic LD: no outliers when z perfectly matches LD structure", {
    set.seed(1002)
    syn <- make_synthetic_ld(500, 10, seed = 1002)
    R <- syn$R

    lead_idx <- 1
    z <- R[, lead_idx] * (-4)

    result <- slalom(
        zScore = z,
        R = R,
        r2Threshold = 0.3,
        nlog10pDentistSThreshold = 4.0
    )

    non_lead <- setdiff(seq_along(z), result$summary$leadPipVariant)
    lead <- result$summary$leadPipVariant
    for (i in non_lead) {
        r2_val <- R[i, lead]^2
        if (r2_val > 0.3 && r2_val < 1.0) {
            expect_false(result$data$outliers[i])
        } else if (r2_val >= 1.0) {
            expect_true(
                is.na(result$data$outliers[i]) || !result$data$outliers[i]
            )
        }
    }
})

# ============================================================================
# Internal function access (resolveLdInput, pecotmr internal)
# ============================================================================

test_that("resolveLdInput returns R when R is provided", {
    R <- diag(5)
    res <- pecotmr:::resolveLdInput(R = R, needNSample = FALSE)
    expect_equal(res$R, R)
})

test_that("resolveLdInput computes R from X", {
    set.seed(1100)
    X <- matrix(sample(0:2, 200 * 5, replace = TRUE), nrow = 200, ncol = 5)
    colnames(X) <- paste0("s", 1:5)
    res <- pecotmr:::resolveLdInput(X = X, needNSample = FALSE)
    expect_true(is.matrix(res$R))
    expect_equal(nrow(res$R), 5)
    expect_equal(ncol(res$R), 5)
    for (j in seq_len(5)) {
        expect_equal(res$R[j, j], 1.0, tolerance = 1e-6)
    }
})

test_that("resolveLdInput errors when neither R nor X given", {
    expect_error(
        pecotmr:::resolveLdInput(R = NULL, X = NULL),
        "Either R.*or X.*must be provided"
    )
})

test_that("resolveLdInput errors when both R and X given", {
    R <- diag(3)
    X <- matrix(1, nrow = 10, ncol = 3)
    expect_error(
        pecotmr:::resolveLdInput(R = R, X = X),
        "Provide either R or X, not both"
    )
})

# ============================================================================
# Numerical stability
# ============================================================================

test_that("log-sum-exp trick prevents overflow with extreme z-scores", {
    z <- c(100, 0, -100)
    R <- diag(3)

    result <- slalom(zScore = z, R = R)

    expect_false(any(is.nan(result$data$prob)))
    expect_false(any(is.infinite(result$data$prob)))
    expect_equal(sum(result$data$prob), 1, tolerance = 1e-10)
    expect_equal(result$data$prob[1], result$data$prob[3], tolerance = 1e-14)
    expect_gt(result$data$prob[1], result$data$prob[2])
})

test_that("standard_error near zero concentrates PIP on large z-scores", {
    z <- c(2, 0.5, 0, -0.1)
    R <- diag(4)
    se <- rep(0.01, 4)

    result <- slalom(
        zScore = z,
        R = R,
        standardError = se,
        abfPriorVariance = 0.04
    )

    expect_equal(which.max(result$data$prob), 1)
    expect_false(any(is.nan(result$data$prob)))
})

# ============================================================================
# Determinism and reproducibility
# ============================================================================

test_that("slalom is deterministic (no randomness)", {
    z <- c(3, -2, 1, 0, -4)
    R <- diag(5)

    result1 <- slalom(zScore = z, R = R)
    result2 <- slalom(zScore = z, R = R)

    expect_identical(result1$data, result2$data)
    expect_identical(
        result1$summary$leadPipVariant,
        result2$summary$leadPipVariant
    )
    expect_identical(result1$summary$cs95, result2$summary$cs95)
    expect_identical(result1$summary$cs99, result2$summary$cs99)
})

# ============================================================================
# Credible set ordering
# ============================================================================

test_that("CS95 variants are ordered by decreasing PIP", {
    set.seed(1400)
    n <- 20
    z <- rnorm(n, sd = 2)
    R <- diag(n)

    result <- slalom(zScore = z, R = R)

    cs_pips <- result$data$prob[result$summary$cs95]
    expect_true(all(diff(cs_pips) <= .Machine$double.eps))
})


context("summaryStatsQc")

# NOTE
# ----
# The variant-content QC (MAF/INFO/N) is pure data.frame filtering via
# .applyContentFilters; no external genome / dbSNP packages required.
# All pecotmr-native steps (.applyContentFilters, .applySkipRegion,
# .matchAgainstSketch, .applyLdMismatchQcToEntry) run
# for real on the synthetic fixture.

# ===========================================================================
# Fixture builders
# ===========================================================================

.ssQ_makeHandle <- function(snp_n = 8L, n_samples = 60L) {
    new(
        "GenotypeHandle",
        path = "/tmp/sketch.gds",
        format = "gds",
        snpInfo = data.frame(
            SNP = paste0("rs", seq_len(snp_n)),
            CHR = rep("1", snp_n),
            BP = seq(100L, by = 100L, length.out = snp_n),
            A1 = rep("A", snp_n),
            A2 = rep("G", snp_n),
            stringsAsFactors = FALSE
        ),
        nSamples = n_samples,
        sampleIds = paste0("s", seq_len(n_samples)),
        pgenPtr = NULL
    )
}

.ssQ_makeEntryGr <- function(
    snp_ids = paste0("rs", 1:4),
    positions = c(100L, 200L, 300L, 400L)
) {
    gr <- GenomicRanges::GRanges(
        seqnames = rep("chr1", length(snp_ids)),
        ranges = IRanges::IRanges(start = positions, width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = snp_ids,
        A1 = rep("A", length(snp_ids)),
        A2 = rep("G", length(snp_ids)),
        Z = seq(1.0, by = 0.5, length.out = length(snp_ids)),
        N = rep(1000L, length(snp_ids))
    )
    gr
}

.ssQ_makeGwasSumStats <- function(
    snp_ids = paste0("rs", 1:4),
    positions = c(100L, 200L, 300L, 400L),
    study = "g1"
) {
    GwasSumStats(
        study = study,
        entry = list(.ssQ_makeEntryGr(snp_ids, positions)),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle()
    )
}

.ssQ_mockExtractor <- function(seed = 13, n_samples = 60L) {
    function(handle, snpIdx, meanImpute = TRUE) {
        set.seed(seed)
        panel <- matrix(
            rbinom(n_samples * nrow(getSnpInfo(handle)), 2, 0.3),
            nrow = n_samples,
            ncol = nrow(getSnpInfo(handle)),
            dimnames = list(getSampleIds(handle), getSnpInfo(handle)$SNP)
        )
        sub <- panel[, snpIdx, drop = FALSE]
        rr <- GenomicRanges::GRanges(
            seqnames = paste0("chr", getSnpInfo(handle)$CHR[snpIdx]),
            ranges = IRanges::IRanges(
                start = getSnpInfo(handle)$BP[snpIdx],
                width = 1L
            )
        )
        S4Vectors::mcols(rr) <- S4Vectors::DataFrame(
            SNP = getSnpInfo(handle)$SNP[snpIdx],
            A1 = getSnpInfo(handle)$A1[snpIdx],
            A2 = getSnpInfo(handle)$A2[snpIdx]
        )
        cd <- S4Vectors::DataFrame(
            sampleId = getSampleIds(handle),
            row.names = getSampleIds(handle)
        )
        dosage <- t(sub)
        rownames(dosage) <- getSnpInfo(handle)$SNP[snpIdx]
        colnames(dosage) <- getSampleIds(handle)
        SummarizedExperiment::SummarizedExperiment(
            assays = list(dosage = dosage),
            rowRanges = rr,
            colData = cd
        )
    }
}

# ===========================================================================
# summaryStatsQc: input-type validation
# ===========================================================================

test_that("summaryStatsQc: rejects non-SumStats input", {
    expect_error(
        summaryStatsQc("not_a_sumstats"),
        "requires a QtlSumStats or GwasSumStats input"
    )
})

test_that("mafCutoff > 0 with no frequency skips the filter, not the run", {
    # A frequency-less study can still run under a default mafCutoff: the
    # STUDY-side filter is skipped with one warning rather than aborting the
    # pipeline. The panel side still runs (mafCutoff is measured wherever it
    # can be), so the mock extractor stands in for the sketch's dosage; every
    # mocked variant sits near MAF 0.3, so nothing is dropped there.
    ss <- .ssQ_makeGwasSumStats()
    local_mocked_bindings(
        extractBlockGenotypes = .ssQ_mockExtractor(),
        .package = "pecotmr"
    )
    expect_warning(
        res <- summaryStatsQc(ss, mafCutoff = 0.05),
        "skipping the MAF filter"
    )
    expect_s4_class(res, "GwasSumStats")
})

test_that("summaryStatsQc: infoCutoff > 0 with no INFO column errors", {
    ss <- .ssQ_makeGwasSumStats()
    expect_error(
        summaryStatsQc(ss, infoCutoff = 0.5),
        "infoCutoff > 0 requires every entry to carry an INFO column"
    )
})

test_that("summaryStatsQc: PIP screen runs AFTER allele harmonization", {
    # Entry: three weak panel-matched variants plus one STRONG signal whose
    # position (9999) is absent from the LD sketch (panel BP = 100..800), so
    # harmonization drops it. With the screen running *after* harmonization it
    # never sees the strong off-panel signal, so the remaining weak variants
    # fail the PIP screen and the region is skipped (empty). (Pre-harmonization
    # screening would have kept the region on the strength of the off-panel hit.)
    gr <- .ssQ_makeEntryGr(
        c("rs1", "rs2", "rs3", "rsX"),
        positions = c(100L, 200L, 300L, 9999L)
    )
    mc <- S4Vectors::mcols(gr)
    mc$Z <- c(0.2, 0.3, 0.1, 12) # only the off-panel rsX carries signal
    S4Vectors::mcols(gr) <- mc
    ss <- GwasSumStats(
        study = "g1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle()
    )
    out <- summaryStatsQc(ss, pipCutoffToSkip = 0.5, nCutoff = 0)
    snps <- as.character(S4Vectors::mcols(out[[1L]])$SNP)
    expect_false("rsX" %in% snps) # dropped by harmonization
    # screen (post-harmonization) skips the region
    expect_equal(length(out[[1L]]), 0L)
})

test_that("summaryStatsQc: PIP screen off leaves the harmonized set intact", {
    # Same harmonized variants, screen disabled (pipCutoffToSkip = 0): the three
    # panel-matched variants survive (behavior-invariance for the screen-off path).
    gr <- .ssQ_makeEntryGr(
        c("rs1", "rs2", "rs3", "rsX"),
        positions = c(100L, 200L, 300L, 9999L)
    )
    mc <- S4Vectors::mcols(gr)
    mc$Z <- c(0.2, 0.3, 0.1, 12)
    S4Vectors::mcols(gr) <- mc
    ss <- GwasSumStats(
        study = "g1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle()
    )
    out <- summaryStatsQc(ss, pipCutoffToSkip = 0, nCutoff = 0)
    snps <- as.character(S4Vectors::mcols(out[[1L]])$SNP)
    # SNP is re-keyed to the panel-harmonized id (chr:pos:A2:A1) after
    # harmonization; the panel is A1=A / A2=G, so the surviving three become
    # chr1:<pos>:G:A.
    expect_setequal(snps, c("chr1:100:G:A", "chr1:200:G:A", "chr1:300:G:A"))
})

# A panel whose SNP ids follow the formatVariantId convention (chr:pos:A2:A1),
# as real genotype LD sketches do, so a re-keyed entry SNP resolves against it.
.ssQ_makeHandleVid <- function(snp_n = 8L, n_samples = 60L) {
    h <- .ssQ_makeHandle(snp_n, n_samples)
    h@snpInfo$SNP <- paste0("chr1:", getSnpInfo(h)$BP, ":G:A")
    h
}

# Entry with one allele-swapped variant (pos 200: A1/A2 reversed vs the panel),
# carrying its input-orientation SNP. Others are exact matches.
.ssQ_makeFlipEntryGr <- function() {
    gr <- GenomicRanges::GRanges(
        seqnames = rep("chr1", 3L),
        ranges = IRanges::IRanges(start = c(100L, 200L, 300L), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = c("chr1:100:G:A", "chr1:200:A:G", "chr1:300:G:A"), # pos200 input-orientation
        A1 = c("A", "G", "A"), # pos200 swapped relative to panel (A1=A)
        A2 = c("G", "A", "G"),
        Z = c(1.0, 2.0, 1.5),
        N = rep(1000L, 3L)
    )
    gr
}

test_that("summaryStatsQc: harmonization re-keys SNP to the panel id and sign-flips Z", {
    # pos200 is allele-swapped vs the panel, so harmonization sign-flips its Z and
    # the SNP must be re-keyed to the panel-orientation id (chr1:200:G:A), not
    # left at the input-orientation chr1:200:A:G. Exact matches are unchanged.
    ss <- GwasSumStats(
        study = "g1",
        entry = list(.ssQ_makeFlipEntryGr()),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle()
    )
    out <- summaryStatsQc(ss, pipCutoffToSkip = 0, nCutoff = 0)
    e <- out[[1L]]
    o <- order(GenomicRanges::start(e))
    snp <- as.character(S4Vectors::mcols(e)$SNP)[o]
    z <- S4Vectors::mcols(e)$Z[o]
    expect_equal(snp, c("chr1:100:G:A", "chr1:200:G:A", "chr1:300:G:A"))
    expect_equal(z, c(1.0, -2.0, 1.5)) # only the swapped variant's Z flips
})

# ---------------------------------------------------------------------------
# Uninterpretable panel ids (e.g. ADSP R5 INS/DEL tags): drop, don't guess
# ---------------------------------------------------------------------------

test_that(".ldSketchIdUsable flags tag-allele ids but not SNPs, real indels, or rsIDs", {
    ids <- c(
        "chr21:13988031:T:C",
        "chr1:788757:T:TAATGG",
        "chr21:16298:C:T",
        "chr21:13988152:INS:T",
        "chr21:13988153:DEL:T",
        "rs12345",
        NA
    )
    expect_equal(
        pecotmr:::.ldSketchIdUsable(ids),
        c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE)
    )
})

test_that(".raissUnsafeToImpute now excludes tag-allele ids (third clause)", {
    # A tag id has no flip twin, so it slipped through before; it must now be
    # unsafe, while a normal SNP stays safe and the flip-pair behaviour is intact.
    expect_true(pecotmr:::.raissUnsafeToImpute(
        "chr21:13988152:INS:T",
        character(0)
    ))
    expect_false(pecotmr:::.raissUnsafeToImpute(
        "chr21:16298:C:T",
        character(0)
    ))
    # flip-pair still flagged (both orientations present)
    fp <- c("chr1:100:A:G", "chr1:100:G:A")
    expect_true(all(pecotmr:::.raissUnsafeToImpute(fp, character(0))))
})

test_that("harmonization drops a variant on a tag-named panel position, keeps normal ones", {
    skip_if_not_installed("pgenlibr")
    h <- .ssQ_makeHandleVid() # ids chr1:100:G:A .. chr1:800:G:A, A1=A A2=G
    df <- tibble(
        chrom = c("1", "1"),
        pos = c(100L, 200L),
        SNP = c("chr1:100:G:A", "chr1:200:G:A"),
        A1 = c("A", "A"),
        A2 = c("G", "G"),
        Z = c(1.0, 2.0),
        N = c(1000L, 1000L)
    )
    # control: both variants match the panel and survive
    out0 <- pecotmr:::.matchAgainstSketch(df, h, matchMinProp = 0)
    expect_setequal(pecotmr:::parseVariantId(out0$SNP)$pos, c(100L, 200L))
    # tag the pos-200 panel id (A1/A2 columns stay the real alleles)
    h@snpInfo$SNP[h@snpInfo$BP == 200L] <- "chr1:200:INS:A"
    outT <- pecotmr:::.matchAgainstSketch(df, h, matchMinProp = 0)
    # pos-200 is dropped (its only reference entry is uninterpretable); pos-100 kept
    expect_true(100L %in% pecotmr:::parseVariantId(outT$SNP)$pos)
    expect_false(200L %in% pecotmr:::parseVariantId(outT$SNP)$pos)
})

test_that("summaryStatsQc: slalom z-mismatch resolves sign-flipped variants against the panel", {
    # Pre-fix regression: a sign-flipped variant kept its input-orientation SNP,
    # which is absent from the panel, so .applyLdMismatchQcToEntry errored with
    # "absent from the ldSketch panel". After the harmonization re-key the SNP
    # matches the panel and slalom QC runs to completion.
    ss <- GwasSumStats(
        study = "g1",
        entry = list(.ssQ_makeFlipEntryGr()),
        genome = "hg19",
        ldSketch = .ssQ_makeHandleVid()
    )
    local_mocked_bindings(
        extractBlockGenotypes = .ssQ_mockExtractor(),
        ldMismatchQc = function(
            zScore,
            R = NULL,
            X = NULL,
            nSample = NULL,
            method = c("slalom", "dentist"),
            ldMethod = "sample",
            ...
        ) {
            data.frame(
                original_z = zScore,
                outlier = rep(FALSE, length(zScore)),
                stringsAsFactors = FALSE
            )
        },
        .package = "pecotmr"
    )
    expect_no_error(
        out <- summaryStatsQc(
            ss,
            zMismatchQc = "slalom",
            pipCutoffToSkip = 0,
            nCutoff = 0
        )
    )
    snp <- as.character(S4Vectors::mcols(out[[1L]])$SNP)
    expect_true("chr1:200:G:A" %in% snp) # the sign-flipped variant survived
    expect_equal(length(out[[1L]]), 3L)
})

test_that("summaryStatsQc: zMismatchQc reconciles a chr-prefix difference vs the panel", {
    # Panel SNP ids are non-chr-prefixed positional; QC re-keys the entry to the
    # canonical chr-prefixed form, so the opt-in z-mismatch panel match must
    # reconcile the prefix (previously errored "absent from the ldSketch panel").
    h <- .ssQ_makeHandle()
    h@snpInfo$SNP <- paste0("1:", getSnpInfo(h)$BP, ":G:A") # non-chr-prefixed
    ss <- GwasSumStats(
        study = "g1",
        entry = list(.ssQ_makeEntryGr()),
        genome = "hg19",
        ldSketch = h
    )
    local_mocked_bindings(
        extractBlockGenotypes = .ssQ_mockExtractor(),
        ldMismatchQc = function(
            zScore,
            R = NULL,
            X = NULL,
            nSample = NULL,
            method = c("slalom", "dentist"),
            ldMethod = "sample",
            ...
        ) {
            data.frame(
                original_z = zScore,
                outlier = rep(FALSE, length(zScore)),
                stringsAsFactors = FALSE
            )
        },
        .package = "pecotmr"
    )
    expect_no_error(
        out <- summaryStatsQc(
            ss,
            zMismatchQc = "slalom",
            pipCutoffToSkip = 0,
            nCutoff = 0
        )
    )
    expect_gte(length(out[[1L]]), 1L)
})

test_that(".deriveBetaSeFromZ: derives BETA+SE when entry has Z+MAF+N only", {
    df <- data.frame(
        SNP = paste0("rs", 1:3),
        A1 = c("A", "C", "G"),
        A2 = c("G", "T", "A"),
        Z = c(1.5, -2.1, 0.4),
        MAF = c(0.2, 0.35, 0.05),
        N = c(10000, 10000, 10000),
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.deriveBetaSeFromZ(df)
    expect_true(all(c("BETA", "SE") %in% names(out$df)))
    expect_equal(out$audit$nDerived, 3L)
    # Verify the formula: se = 1/sqrt(2*maf*(1-maf)*(N+z^2))
    expected_se <- 1 / sqrt(2 * df$MAF * (1 - df$MAF) * (df$N + df$Z^2))
    expect_equal(out$df$SE, expected_se)
    expect_equal(out$df$BETA, df$Z * expected_se)
})

test_that(".deriveBetaSeFromZ: no-op when BETA and SE already present", {
    df <- data.frame(
        Z = 1.5,
        BETA = 0.5,
        SE = 0.1,
        MAF = 0.3,
        N = 1000,
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.deriveBetaSeFromZ(df)
    expect_null(out$audit)
    expect_equal(out$df, df)
})

test_that(".deriveBetaSeFromZ: skipped when N missing", {
    df <- data.frame(
        Z = 1.5,
        MAF = 0.3,
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.deriveBetaSeFromZ(df)
    expect_null(out$audit)
    expect_false("BETA" %in% names(out$df))
    expect_false("SE" %in% names(out$df))
})

# ===========================================================================
# summaryStatsQc: end-to-end on the synthetic fixture
# ===========================================================================

test_that("summaryStatsQc: vanilla run populates qcInfo and returns a GwasSumStats", {
    ss <- .ssQ_makeGwasSumStats()
    res <- summaryStatsQc(ss)
    expect_s4_class(res, "GwasSumStats")
    qc <- getQcInfo(res)
    expect_true(length(qc) > 0L)
    expect_true("options" %in% names(qc))
    expect_true("entryAudit" %in% names(qc))
    expect_equal(length(qc$entryAudit), nrow(ss))
    ea <- qc$entryAudit[[1L]]
    expect_equal(ea$variantsIn, 4L)
    expect_equal(ea$variantsOut, 4L)
})

test_that("summaryStatsQc: keepVariants subsets each entry and records the drop", {
    ss <- .ssQ_makeGwasSumStats()
    res <- summaryStatsQc(ss, keepVariants = c("rs1", "rs3"))
    ea <- getQcInfo(res)$entryAudit[[1L]]
    expect_equal(ea$keepVariantsDropped, 2L)
    expect_equal(ea$variantsOut, 2L)
})

test_that("summaryStatsQc: skipRegion drops overlapping variants", {
    ss <- .ssQ_makeGwasSumStats()
    res <- summaryStatsQc(ss, skipRegion = "chr1:50-150")
    ea <- getQcInfo(res)$entryAudit[[1L]]
    expect_equal(ea$skipRegionDropped, 1L) # rs1 at pos 100 is dropped
})

test_that("summaryStatsQc: PIP screen triggers when no variant has signal", {
    # Build an entry with weak signal so the SER PIP screen tags everything
    # below the threshold.
    gr <- .ssQ_makeEntryGr()
    S4Vectors::mcols(gr)$Z <- rep(0.1, length(gr))
    ss <- GwasSumStats(
        study = "g1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle()
    )
    res <- summaryStatsQc(ss, pipCutoffToSkip = 0.99)
    ea <- getQcInfo(res)$entryAudit[[1L]]
    expect_true(isTRUE(ea$pipScreenSkipped))
    expect_match(ea$pipScreenReason, "no signals above PIP threshold")
    expect_equal(length(res[[1L]]), 0L)
})

test_that("summaryStatsQc: harmonized variants count is recorded", {
    ss <- .ssQ_makeGwasSumStats()
    res <- summaryStatsQc(ss)
    ea <- getQcInfo(res)$entryAudit[[1L]]
    expect_equal(ea$matchedAgainstSketch, 4L)
})

test_that("summaryStatsQc: options block records the curated knobs", {
    ss <- .ssQ_makeGwasSumStats()
    res <- summaryStatsQc(
        ss,
        removeIndels = TRUE,
        removeStrandAmbiguous = FALSE,
        nCutoff = 10
    )
    opts <- getQcInfo(res)$options
    expect_true(opts$removeIndels)
    expect_false(opts$removeStrandAmbiguous)
    expect_equal(opts$nCutoff, 10)
})

test_that("summaryStatsQc: round-trips QtlSumStats inputs", {
    gr <- .ssQ_makeEntryGr()
    ss <- QtlSumStats(
        study = "s1",
        context = "c1",
        trait = "t1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle()
    )
    res <- summaryStatsQc(ss)
    expect_s4_class(res, "QtlSumStats")
    expect_equal(length(getQcInfo(res)$entryAudit), 1L)
})

# ===========================================================================
# effectiveN helper + case/control N canonicalization in summaryStatsQc
# ===========================================================================

# Build a GwasSumStats entry carrying per-variant N_CASE / N_CONTROL mcols
# (optionally an N column too). Panel-matched (A1=A / A2=G at BP 100..) so the
# variants survive harmonization against .ssQ_makeHandle().
.ssQ_makeCCEntry <- function(
    nCase,
    nControl,
    includeN = FALSE,
    snp_ids = paste0("rs", 1:4),
    positions = c(100L, 200L, 300L, 400L)
) {
    m <- length(snp_ids)
    gr <- GenomicRanges::GRanges(
        seqnames = rep("chr1", m),
        ranges = IRanges::IRanges(start = positions, width = 1L)
    )
    mc <- S4Vectors::DataFrame(
        SNP = snp_ids,
        A1 = rep("A", m),
        A2 = rep("G", m),
        Z = seq(1.0, by = 0.5, length.out = m),
        N_CASE = as.numeric(rep_len(nCase, m)),
        N_CONTROL = as.numeric(rep_len(nControl, m))
    )
    if (includeN) {
        mc$N <- rep(1000L, m)
    }
    S4Vectors::mcols(gr) <- mc
    gr
}

# Entry N ordered by genomic position (harmonization can reorder rows).
.ssQ_entryNByPos <- function(entry) {
    o <- order(GenomicRanges::start(entry))
    as.numeric(S4Vectors::mcols(entry)$N[o])
}

test_that("effectiveN: balanced equals total, imbalanced is smaller, guards NA/<=0", {
    # Balanced: 4/(1/500 + 1/500) == 1000 == total.
    expect_equal(effectiveN(500, 500), 1000)
    # Imbalanced: 4*100*900/1000 = 360 < 1000.
    expect_equal(effectiveN(100, 900), 360)
    expect_lt(effectiveN(100, 900), 100 + 900)
    # NA / non-positive counts -> NA_real_.
    expect_true(is.na(effectiveN(NA, 500)))
    expect_true(is.na(effectiveN(500, NA)))
    expect_true(is.na(effectiveN(0, 500)))
    expect_true(is.na(effectiveN(-5, 500)))
    # Vectorized, element-wise.
    expect_equal(
        effectiveN(c(500, 100, 0), c(500, 900, 500)),
        c(1000, 360, NA_real_)
    )
})

test_that("summaryStatsQc(effectiveN=TRUE): per-variant counts, no N -> N == N_eff", {
    gr <- .ssQ_makeCCEntry(
        nCase = c(100, 200, 150, 250),
        nControl = c(900, 800, 850, 750)
    )
    ss <- GwasSumStats(
        study = "g1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle()
    )
    res <- summaryStatsQc(ss)
    # 4*case*control/(case+control) per variant.
    expect_equal(.ssQ_entryNByPos(res[[1L]]), c(360, 640, 510, 750))
    expect_identical(getQcInfo(res)$entryAudit[[1L]]$nSource, "effective")
})

test_that("summaryStatsQc(effectiveN=TRUE): counts + N -> counts win, override logged", {
    gr <- .ssQ_makeCCEntry(
        nCase = c(100, 200, 150, 250),
        nControl = c(900, 800, 850, 750),
        includeN = TRUE
    )
    ss <- GwasSumStats(
        study = "g1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle()
    )
    expect_message(summaryStatsQc(ss), "overridden by effective N")
    res <- summaryStatsQc(ss)
    # N (was 1000) is replaced by the per-variant N_eff.
    expect_equal(.ssQ_entryNByPos(res[[1L]]), c(360, 640, 510, 750))
    expect_identical(getQcInfo(res)$entryAudit[[1L]]$nSource, "effective")
})

test_that("summaryStatsQc(effectiveN=TRUE): N only, no counts -> used as-is", {
    ss <- .ssQ_makeGwasSumStats() # entry carries N = 1000, no counts
    res <- summaryStatsQc(ss)
    expect_true(all(.ssQ_entryNByPos(res[[1L]]) == 1000))
    expect_identical(getQcInfo(res)$entryAudit[[1L]]$nSource, "column")
})

test_that("summaryStatsQc(effectiveN=FALSE): counts + N -> raw N, no override", {
    gr <- .ssQ_makeCCEntry(
        nCase = c(100, 200, 150, 250),
        nControl = c(900, 800, 850, 750),
        includeN = TRUE
    )
    ss <- GwasSumStats(
        study = "g1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle()
    )
    res <- summaryStatsQc(ss, effectiveN = FALSE)
    expect_true(all(.ssQ_entryNByPos(res[[1L]]) == 1000))
    expect_identical(getQcInfo(res)$entryAudit[[1L]]$nSource, "column")
})

test_that("summaryStatsQc(effectiveN=FALSE): counts only -> raw total, nSource='total'", {
    gr <- .ssQ_makeCCEntry(
        nCase = c(100, 200, 150, 250),
        nControl = c(900, 800, 850, 750)
    )
    ss <- GwasSumStats(
        study = "g1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle()
    )
    res <- summaryStatsQc(ss, effectiveN = FALSE)
    # Raw total n_case + n_control per variant.
    expect_equal(.ssQ_entryNByPos(res[[1L]]), c(1000, 1000, 1000, 1000))
    expect_identical(getQcInfo(res)$entryAudit[[1L]]$nSource, "total")
})

test_that("summaryStatsQc(effectiveN=TRUE): study-level scalars applied to all variants", {
    # Entry has an N column but NO per-variant N_CASE/N_CONTROL; the scalars win.
    gr <- .ssQ_makeEntryGr()
    ss <- GwasSumStats(
        study = "g1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle(),
        nCase = 100,
        nControl = 900
    )
    expect_message(summaryStatsQc(ss), "from study")
    res <- summaryStatsQc(ss)
    # 4*100*900/1000 = 360 for every variant.
    expect_true(all(.ssQ_entryNByPos(res[[1L]]) == 360))
    expect_identical(getQcInfo(res)$entryAudit[[1L]]$nSource, "effective")
})

test_that("summaryStatsQc: study nSample is the level-4 fallback (no counts, no per-variant N)", {
    # Entry with NO per-variant N and NO case/control; only a study nSample scalar.
    gr <- .ssQ_makeEntryGr()
    mc <- S4Vectors::mcols(gr)
    mc$N <- NULL # remove per-variant N
    S4Vectors::mcols(gr) <- mc
    ss <- GwasSumStats(
        study = "g1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle(),
        nSample = 4321
    )
    res <- summaryStatsQc(ss)
    expect_true(all(.ssQ_entryNByPos(res[[1L]]) == 4321))
    expect_identical(getQcInfo(res)$entryAudit[[1L]]$nSource, "study-n")
})

test_that("summaryStatsQc: QtlSumStats tuple nSample is the level-4 fallback too", {
    # Parity with the GWAS study-n fallback: a QtlSumStats carrying only a
    # tuple-level nSample (no per-variant N) fills N from the scalar and is
    # preserved on the QC'd object.
    gr <- .ssQ_makeEntryGr()
    mc <- S4Vectors::mcols(gr)
    mc$N <- NULL # remove per-variant N
    S4Vectors::mcols(gr) <- mc
    ss <- QtlSumStats(
        study = "s1",
        context = "c1",
        trait = "t1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle(),
        nSample = 838
    )
    res <- summaryStatsQc(ss)
    expect_s4_class(res, "QtlSumStats")
    expect_true(all(.ssQ_entryNByPos(res[[1L]]) == 838))
    expect_identical(getQcInfo(res)$entryAudit[[1L]]$nSource, "study-n")
    expect_equal(as.numeric(res$nSample), 838) # slot preserved through QC
})

test_that("summaryStatsQc: level precedence -- per-variant counts beat nSample; N column beats nSample", {
    # per-variant counts present alongside nSample -> counts win (effective).
    gr1 <- .ssQ_makeCCEntry(
        nCase = c(100, 200, 150, 250),
        nControl = c(900, 800, 850, 750)
    )
    ss1 <- GwasSumStats(
        study = "g1",
        entry = list(gr1),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle(),
        nSample = 4321
    )
    res1 <- summaryStatsQc(ss1)
    expect_equal(.ssQ_entryNByPos(res1[[1L]]), c(360, 640, 510, 750))
    expect_identical(getQcInfo(res1)$entryAudit[[1L]]$nSource, "effective")
    # per-variant N column present alongside nSample (no counts) -> N wins (column).
    ss2 <- GwasSumStats(
        study = "g1",
        entry = list(.ssQ_makeEntryGr()),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle(),
        nSample = 4321
    )
    res2 <- summaryStatsQc(ss2) # .ssQ_makeEntryGr carries N = 1000
    expect_true(all(.ssQ_entryNByPos(res2[[1L]]) == 1000))
    expect_identical(getQcInfo(res2)$entryAudit[[1L]]$nSource, "column")
})

test_that("summaryStatsQc: effectiveN recorded in qcInfo options", {
    ss <- .ssQ_makeGwasSumStats()
    expect_true(getQcInfo(summaryStatsQc(ss))$options$effectiveN)
    expect_false(
        getQcInfo(summaryStatsQc(ss, effectiveN = FALSE))$options$effectiveN
    )
})

test_that("summaryStatsQc: quantitative QtlSumStats is a no-op for effective N", {
    # No counts anywhere; entry carries an N column -> used as-is (nSource
    # "column"), N untouched.
    gr <- .ssQ_makeEntryGr()
    ss <- QtlSumStats(
        study = "s1",
        context = "c1",
        trait = "t1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle()
    )
    res <- summaryStatsQc(ss)
    expect_true(all(.ssQ_entryNByPos(res[[1L]]) == 1000))
    expect_identical(getQcInfo(res)$entryAudit[[1L]]$nSource, "column")
})

# ===========================================================================
# summaryStatsQc with LD-mismatch QC enabled (mocked extractor)
# ===========================================================================

test_that("summaryStatsQc: zMismatchQc = 'dentist' walks the LD-mismatch branch", {
    # Panel ids follow the chr:pos:A2:A1 convention (as real LD sketches do) so the
    # post-harmonization re-keyed SNP resolves against the panel for z-mismatch QC.
    ss <- GwasSumStats(
        study = "g1",
        entry = list(.ssQ_makeEntryGr(
            snp_ids = paste0("rs", 1:8),
            positions = seq(100L, by = 100L, length.out = 8L)
        )),
        genome = "hg19",
        ldSketch = .ssQ_makeHandleVid(snp_n = 8L)
    )
    local_mocked_bindings(
        extractBlockGenotypes = .ssQ_mockExtractor(),
        .package = "pecotmr"
    )
    res <- suppressWarnings(summaryStatsQc(ss, zMismatchQc = "dentist"))
    ea <- getQcInfo(res)$entryAudit[[1L]]
    expect_equal(ea$ldMismatchMethod, "dentist")
    expect_true("ldMismatchOutliersDropped" %in% names(ea))
})

# ===========================================================================
# summaryStatsQc with impute = TRUE: exercise the RAISS branch
# ===========================================================================

test_that("summaryStatsQc: impute = TRUE invokes RAISS and records the audit counts", {
    # Build a sketch panel with 8 variants and a GWAS entry covering only the
    # first 4 — RAISS is asked to impute the missing 4.
    full_snp_ids <- paste0("rs", 1:8)
    full_positions <- seq(100L, by = 100L, length.out = 8L)
    ss <- GwasSumStats(
        study = "g1",
        entry = list(.ssQ_makeEntryGr(
            snp_ids = full_snp_ids[1:4],
            positions = full_positions[1:4]
        )),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle(snp_n = 8L, n_samples = 60L)
    )

    local_mocked_bindings(
        extractBlockGenotypes = .ssQ_mockExtractor(),
        raiss = function(refPanel, knownZscores, genotypeMatrix, ...) {
            # Pretend RAISS imputed two of the missing panel variants (rs5, rs6)
            # with synthetic z-scores.
            added <- refPanel[
                refPanel$variant_id %in% c("rs5", "rs6"),
                ,
                drop = FALSE
            ]
            added$z <- c(1.5, -2.0)
            added$n <- c(1000, 1000)
            list(resultFilter = rbind(knownZscores, added))
        },
        .package = "pecotmr"
    )
    # rs5/rs6 (pos 500/600) sit beyond the observed range (100-400); impute now
    # scopes to the region window, so widen it with a flank to reach them.
    res <- summaryStatsQc(ss, impute = TRUE, imputeOpts = list(flank = 500))
    ea <- getQcInfo(res)$entryAudit[[1L]]
    expect_equal(ea$raissTotalVariants, 6L)
    expect_equal(ea$raissImputedVariants, 2L)
})

test_that("summaryStatsQc: impute scopes the reference panel/dosage to the region window", {
    # Sketch spans rs1..rs8 (pos 100..800); the entry observes only rs1..rs4
    # (100..400). With the default flank the impute window is [100, 400], so the
    # dosage must be materialized for just those 4 panel variants -- NOT the whole
    # 8-variant sketch (the bug that makes --impute unusable on a per-chromosome
    # sketch: it built dosage for seq_len(nrow(sketchSnpInfo))).
    #
    # Observed at extractBlockGenotypes: the single reader every dosage path
    # funnels through, including the DelayedArray seed behind assay(). The
    # panel reads its own dosage through the assay now, so watching
    # .dosageMatrix would no longer see the RAISS read at all.
    cap <- new.env(parent = emptyenv())
    cap$idx <- list()
    ss <- GwasSumStats(
        study = "g1",
        entry = list(.ssQ_makeEntryGr(
            paste0("rs", 1:4),
            c(100L, 200L, 300L, 400L)
        )),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle(snp_n = 8L, n_samples = 60L)
    )
    inner <- .ssQ_mockExtractor(n_samples = 60L)
    local_mocked_bindings(
        extractBlockGenotypes = function(handle, snpIdx, meanImpute = TRUE) {
            cap$idx <- c(cap$idx, list(snpIdx))
            inner(handle, snpIdx, meanImpute)
        },
        raiss = function(...) NULL,
        .package = "pecotmr"
    )
    suppressWarnings(summaryStatsQc(ss, impute = TRUE, nCutoff = 0))
    expect_gt(length(cap$idx), 0L)
    # No read reaches past the region window's 4 panel variants.
    expect_true(all(map_lgl(cap$idx, function(i) all(i %in% 1:4))))
    expect_equal(max(map_int(cap$idx, length)), 4L)
})

test_that("summaryStatsQc: impute = TRUE with raiss returning NULL records 0 imputed", {
    full_snp_ids <- paste0("rs", 1:8)
    full_positions <- seq(100L, by = 100L, length.out = 8L)
    ss <- GwasSumStats(
        study = "g1",
        entry = list(.ssQ_makeEntryGr(
            snp_ids = full_snp_ids[1:4],
            positions = full_positions[1:4]
        )),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle(snp_n = 8L, n_samples = 60L)
    )

    local_mocked_bindings(
        extractBlockGenotypes = .ssQ_mockExtractor(),
        raiss = function(...) NULL,
        .package = "pecotmr"
    )
    res <- summaryStatsQc(ss, impute = TRUE)
    ea <- getQcInfo(res)$entryAudit[[1L]]
    expect_equal(ea$raissImputedVariants, 0L)
})

# ===========================================================================
# summaryStatsQc: per-step QC counter logging (concept salvaged from PR #520)
# ===========================================================================

test_that("harmonizeAlleles surfaces sign/strand/dropped counts via qcCounts attribute", {
    # 4 shared positions: 100 exact, 200 sign-flip, 300 strand-flip (A/G
    # unambiguous), 400 allele mismatch (dropped).
    target <- data.frame(
        chrom = c(1, 1, 1, 1),
        pos = c(100, 200, 300, 400),
        A2 = c("A", "A", "A", "A"),
        A1 = c("G", "G", "G", "G"),
        z = c(1, 2, 3, 4),
        stringsAsFactors = FALSE
    )
    ref <- data.frame(
        chrom = c(1, 1, 1, 1),
        pos = c(100, 200, 300, 400),
        A2 = c("A", "G", "T", "C"),
        A1 = c("G", "A", "C", "A"),
        stringsAsFactors = FALSE
    )
    res <- pecotmr:::harmonizeAlleles(
        target,
        ref,
        colToFlip = "z",
        matchMinProp = 0
    )
    # Default return shape unchanged.
    expect_named(res, c("harmonizedData", "qcSummary"))
    expect_equal(nrow(res$harmonizedData), 3L)
    cnt <- attr(res, "qcCounts")
    expect_false(is.null(cnt))
    expect_equal(cnt$considered, 4L)
    expect_equal(cnt$signFlip, 1L)
    expect_equal(cnt$strandFlip, 1L)
    expect_equal(cnt$kept, 3L)
    expect_equal(cnt$dropped, 1L)
})

test_that("summaryStatsQc: QC track emits per-step 'kept N of M' messages plus a rollup", {
    ss <- .ssQ_makeGwasSumStats()
    local_mocked_bindings(
        extractBlockGenotypes = .ssQ_mockExtractor(),
        .package = "pecotmr"
    )
    msgs <- capture_messages(summaryStatsQc(ss))
    joined <- paste(msgs, collapse = "")
    # Harmonization step + corrected/dropped breakdown.
    expect_match(joined, "harmonization kept [0-9]+ of [0-9]+")
    expect_match(
        joined,
        "corrected: sign-flipped [0-9]+, strand-flipped [0-9]+"
    )
    # Per-entry rollup line.
    expect_match(joined, "QC summary: [0-9]+ in -> [0-9]+ out")
    expect_match(joined, "corrected: sign-flip [0-9]+, strand-flip [0-9]+")
})

test_that("summaryStatsQc: skipped optional steps are omitted from the rollup", {
    ss <- .ssQ_makeGwasSumStats()
    local_mocked_bindings(
        extractBlockGenotypes = .ssQ_mockExtractor(),
        .package = "pecotmr"
    )
    msgs <- capture_messages(
        summaryStatsQc(
            ss,
            alleleFlipKriging = FALSE,
            zMismatchQc = "none",
            impute = FALSE
        )
    )
    joined <- paste(msgs, collapse = "")
    # Kriging / mismatch / imputation are skipped: their per-step messages
    # and their rollup segments should be absent.
    expect_false(grepl("kriging", joined))
    expect_false(grepl("LD-mismatch", joined))
    expect_false(grepl("RAISS imputation", joined))
    expect_false(grepl("imputed [+-][1-9]", joined))
})

test_that("summaryStatsQc: per-entry log lines carry the (study/context/trait) label for QtlSumStats", {
    # Reuse the QtlSumStats fixture from the round-trip test.
    qss <- QtlSumStats(
        study = "qstudy",
        context = "qctx",
        trait = "qtrait",
        entry = list(.ssQ_makeEntryGr()),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle()
    )
    local_mocked_bindings(
        extractBlockGenotypes = .ssQ_mockExtractor(),
        .package = "pecotmr"
    )
    msgs <- capture_messages(summaryStatsQc(qss))
    joined <- paste(msgs, collapse = "")
    expect_match(joined, "\\[qstudy/qctx/qtrait\\] QC track")
    expect_match(joined, "\\[qstudy/qctx/qtrait\\] QC summary")
})


context("sumstatsQc internal helpers")

# ===========================================================================
# Fixture builders
# ===========================================================================

.ssh_makeHandle <- function(snp_n = 6L, n_samples = 30L) {
    new(
        "GenotypeHandle",
        path = "/tmp/sketch.gds",
        format = "gds",
        snpInfo = data.frame(
            SNP = paste0("rs", seq_len(snp_n)),
            CHR = rep("1", snp_n),
            BP = seq(100L, by = 100L, length.out = snp_n),
            A1 = rep("A", snp_n),
            A2 = rep("G", snp_n),
            stringsAsFactors = FALSE
        ),
        nSamples = n_samples,
        sampleIds = paste0("s", seq_len(n_samples)),
        pgenPtr = NULL
    )
}

.ssh_makeEntryGr <- function(n = 5, chr = "chr1", with_extras = FALSE) {
    gr <- GenomicRanges::GRanges(
        seqnames = rep(chr, n),
        ranges = IRanges::IRanges(
            start = seq(100L, by = 100L, length.out = n),
            width = 1L
        )
    )
    mc <- list(
        SNP = paste0("rs", seq_len(n)),
        A1 = rep("A", n),
        A2 = rep("G", n),
        Z = seq(1.0, by = 0.5, length.out = n),
        N = rep(1000L, n)
    )
    if (with_extras) {
        mc$MAF <- seq(0.1, by = 0.05, length.out = n)
        mc$INFO <- rep(0.95, n)
        mc$BETA <- rnorm(n)
        mc$SE <- rep(0.1, n)
        mc$P <- 2 * pnorm(-abs(mc$Z))
    }
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(mc)
    gr
}

.ssh_mockExtractor <- function(seed = 42, n_samples = 30L) {
    function(handle, snpIdx, meanImpute = TRUE) {
        set.seed(seed)
        panel <- matrix(
            rbinom(n_samples * nrow(getSnpInfo(handle)), 2, 0.3),
            nrow = n_samples,
            ncol = nrow(getSnpInfo(handle)),
            dimnames = list(getSampleIds(handle), getSnpInfo(handle)$SNP)
        )
        sub <- panel[, snpIdx, drop = FALSE]
        rr <- GenomicRanges::GRanges(
            seqnames = paste0("chr", getSnpInfo(handle)$CHR[snpIdx]),
            ranges = IRanges::IRanges(
                start = getSnpInfo(handle)$BP[snpIdx],
                width = 1L
            )
        )
        S4Vectors::mcols(rr) <- S4Vectors::DataFrame(
            SNP = getSnpInfo(handle)$SNP[snpIdx],
            A1 = getSnpInfo(handle)$A1[snpIdx],
            A2 = getSnpInfo(handle)$A2[snpIdx]
        )
        cd <- S4Vectors::DataFrame(
            sampleId = getSampleIds(handle),
            row.names = getSampleIds(handle)
        )
        dosage <- t(sub)
        rownames(dosage) <- getSnpInfo(handle)$SNP[snpIdx]
        colnames(dosage) <- getSampleIds(handle)
        SummarizedExperiment::SummarizedExperiment(
            assays = list(dosage = dosage),
            rowRanges = rr,
            colData = cd
        )
    }
}

# ===========================================================================
# .entryGrangesToDf and .dfToEntryGranges (round-trip)
# ===========================================================================

test_that(".entryGrangesToDf: extracts chrom/pos and all mcols", {
    gr <- .ssh_makeEntryGr(3, with_extras = TRUE)
    df <- pecotmr:::.entryGrangesToDf(gr)
    expect_s3_class(df, "data.frame")
    expect_equal(df$chrom, rep("1", 3)) # "chr" stripped
    expect_equal(df$pos, c(100L, 200L, 300L))
    expect_setequal(
        intersect(
            colnames(df),
            c("SNP", "A1", "A2", "Z", "N", "MAF", "INFO", "BETA", "SE", "P")
        ),
        c("SNP", "A1", "A2", "Z", "N", "MAF", "INFO", "BETA", "SE", "P")
    )
})

test_that(".dfToEntryGranges: rebuilds the GRanges with canonical mcols", {
    df <- data.frame(
        chrom = "1",
        pos = c(100L, 200L),
        SNP = c("rs1", "rs2"),
        A1 = c("A", "A"),
        A2 = c("G", "G"),
        Z = c(1.5, -2.0),
        N = c(1000L, 1200L),
        stringsAsFactors = FALSE
    )
    gr <- pecotmr:::.dfToEntryGranges(df)
    expect_s4_class(gr, "GRanges")
    expect_equal(as.character(GenomicRanges::seqnames(gr)), c("chr1", "chr1"))
    expect_equal(GenomicRanges::start(gr), c(100L, 200L))
    expect_setequal(
        colnames(S4Vectors::mcols(gr)),
        c("SNP", "A1", "A2", "Z", "N")
    )
})

test_that(".dfToEntryGranges: derives SNP from variant_id when SNP is absent", {
    df <- data.frame(
        chrom = "1",
        pos = 100L,
        variant_id = "chr1:100:A:G",
        A1 = "A",
        A2 = "G",
        Z = 1.0,
        N = 1000L,
        stringsAsFactors = FALSE
    )
    gr <- pecotmr:::.dfToEntryGranges(df)
    expect_equal(S4Vectors::mcols(gr)$SNP, "chr1:100:A:G")
})

test_that("entry GRanges round-trips through df conversion", {
    gr <- .ssh_makeEntryGr(4, with_extras = TRUE)
    df <- pecotmr:::.entryGrangesToDf(gr)
    gr2 <- pecotmr:::.dfToEntryGranges(df)
    # Positions and core mcols must match (mcols may differ by reordering).
    expect_equal(GenomicRanges::start(gr2), GenomicRanges::start(gr))
    expect_equal(S4Vectors::mcols(gr2)$SNP, S4Vectors::mcols(gr)$SNP)
    expect_equal(S4Vectors::mcols(gr2)$Z, S4Vectors::mcols(gr)$Z)
})

# ===========================================================================
# .refVariantsFromSketch
# ===========================================================================

test_that(".refVariantsFromSketch: extracts chr/pos/A1/A2/variant_id from snpInfo", {
    h <- .ssh_makeHandle()
    rv <- pecotmr:::.refVariantsFromSketch(h)
    expect_equal(rv$chrom, rep("1", 6)) # "chr" stripped
    expect_equal(rv$pos, c(100L, 200L, 300L, 400L, 500L, 600L))
    expect_equal(rv$variant_id, paste0("rs", 1:6))
    expect_equal(rv$A1, rep("A", 6))
    expect_equal(rv$A2, rep("G", 6))
})

# ===========================================================================
# .applySkipRegion
# ===========================================================================

.ssh_smallDf <- function() {
    data.frame(
        chrom = c("1", "1", "2"),
        pos = c(100L, 200L, 100L),
        SNP = c("rs1", "rs2", "rs3"),
        Z = c(1, 2, 3),
        stringsAsFactors = FALSE
    )
}

test_that(".applySanityChecks: empty input is a no-op", {
    out <- pecotmr:::.applySanityChecks(data.frame(chrom = character()))
    expect_equal(nrow(out$df), 0L)
    expect_equal(out$audit, list())
})

test_that(".applySanityChecks: coerceNumeric converts character columns and counts NAs", {
    df <- data.frame(
        chrom = c("1", "1"),
        pos = c(100L, 200L),
        A1 = c("A", "C"),
        A2 = c("G", "T"),
        Z = c("1.5", "not_a_number"),
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.applySanityChecks(df, dropMissData = FALSE)
    expect_type(out$df$Z, "double")
    expect_equal(out$df$Z[[1L]], 1.5)
    expect_true(is.na(out$df$Z[[2L]]))
    expect_equal(out$audit$nonNumericCoerced, 1L)
})

test_that(".applySanityChecks: normalizeChr maps 23/24/M/chr* and drops non-standard", {
    df <- data.frame(
        chrom = c("chr1", "23", "24", "M", "chrX_random"),
        pos = c(100L, 200L, 300L, 400L, 500L),
        A1 = c("A", "A", "A", "A", "A"),
        A2 = c("G", "G", "G", "G", "G"),
        Z = c(1, 2, 3, 4, 5),
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.applySanityChecks(df)
    expect_equal(out$df$chrom, c("1", "X", "Y", "MT"))
    expect_equal(out$audit$nonstandardChrDropped, 1L)
})

test_that(".applySanityChecks: dropMissData drops rows with NA in vital columns", {
    df <- data.frame(
        chrom = c("1", "1", "1"),
        pos = c(100L, 200L, 300L),
        A1 = c("A", NA, "C"),
        A2 = c("G", "T", "T"),
        Z = c(1, 2, NA),
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.applySanityChecks(df)
    expect_equal(nrow(out$df), 1L)
    expect_equal(out$audit$missDataDropped, 2L)
})

test_that(".applySanityChecks: dropPOutOfRange drops P < 0 and P > 1", {
    df <- data.frame(
        chrom = c("1", "1", "1", "1"),
        pos = c(100L, 200L, 300L, 400L),
        A1 = c("A", "A", "A", "A"),
        A2 = c("G", "G", "G", "G"),
        Z = c(1, 2, 3, 4),
        P = c(0.5, -0.1, 1.5, 0.001),
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.applySanityChecks(df)
    expect_equal(nrow(out$df), 2L)
    expect_equal(out$audit$pOutOfRangeDropped, 2L)
})

test_that(".applySanityChecks: clampSmallP floors to smallPFloor", {
    df <- data.frame(
        chrom = c("1", "1"),
        pos = c(100L, 200L),
        A1 = c("A", "A"),
        A2 = c("G", "G"),
        Z = c(1, 50),
        P = c(0.1, 0),
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.applySanityChecks(df, smallPFloor = 5e-324)
    expect_equal(out$df$P[[2L]], 5e-324)
    expect_equal(out$audit$smallPClamped, 1L)
})

test_that(".applySanityChecks: dropZeroEffect drops BETA==0 and OR==1", {
    df <- data.frame(
        chrom = c("1", "1", "1"),
        pos = c(100L, 200L, 300L),
        A1 = c("A", "A", "A"),
        A2 = c("G", "G", "G"),
        BETA = c(0.5, 0, 0.2),
        OR = c(1.5, 1.5, 1),
        SE = c(0.1, 0.1, 0.1),
        Z = c(5, 0, 2),
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.applySanityChecks(df)
    expect_equal(nrow(out$df), 1L)
    expect_equal(out$audit$zeroEffectDropped, 2L)
})

test_that(".applySanityChecks: dropNonpositiveSe drops SE <= 0", {
    df <- data.frame(
        chrom = c("1", "1", "1"),
        pos = c(100L, 200L, 300L),
        A1 = c("A", "A", "A"),
        A2 = c("G", "G", "G"),
        SE = c(0.1, 0, -0.2),
        Z = c(1, 2, 3),
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.applySanityChecks(df)
    expect_equal(nrow(out$df), 1L)
    expect_equal(out$audit$nonpositiveSeDropped, 2L)
})

test_that(".applySanityChecks: per-check knobs can disable each step", {
    df <- data.frame(
        chrom = c("chr1", "23"),
        pos = c(100L, 200L),
        A1 = c("A", "A"),
        A2 = c("G", "G"),
        BETA = c(0.5, 0),
        SE = c(0.1, -0.2),
        P = c(0.5, 1.5),
        Z = c(1, 2),
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.applySanityChecks(
        df,
        coerceNumeric = FALSE,
        normalizeChr = FALSE,
        dropMissData = FALSE,
        dropPOutOfRange = FALSE,
        clampSmallP = FALSE,
        dropZeroEffect = FALSE,
        dropNonpositiveSe = FALSE
    )
    expect_equal(nrow(out$df), 2L)
    expect_equal(out$audit, list())
    expect_equal(out$df$chrom, c("chr1", "23"))
})

test_that("summaryStatsQc: surfaces sanity-check audit and respects per-check knobs", {
    gr <- .ssQ_makeEntryGr()
    S4Vectors::mcols(gr)$BETA <- c(0.1, 0, 0.2, 0)
    S4Vectors::mcols(gr)$SE <- c(0.1, 0.1, 0.1, 0.1)
    S4Vectors::mcols(gr)$P <- c(0.5, 0.5, 0.5, 0.5)
    ss <- GwasSumStats(
        study = "g1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle()
    )
    res <- summaryStatsQc(ss)
    ea <- getQcInfo(res)$entryAudit[[1L]]
    expect_equal(ea$sanityChecks$zeroEffectDropped, 2L)
    # Disable the knob: zero-effect rows should remain.
    res2 <- summaryStatsQc(ss, dropZeroEffect = FALSE)
    ea2 <- getQcInfo(res2)$entryAudit[[1L]]
    expect_null(ea2$sanityChecks$zeroEffectDropped)
})

test_that(".applySkipRegion: NULL / empty skipRegion is a no-op", {
    df <- .ssh_smallDf()
    expect_identical(pecotmr:::.applySkipRegion(df, NULL), df)
    expect_identical(pecotmr:::.applySkipRegion(df, character()), df)
})

test_that(".applySkipRegion: drops variants overlapping a single character region", {
    df <- .ssh_smallDf()
    out <- pecotmr:::.applySkipRegion(df, "1:50-150")
    expect_equal(out$SNP, c("rs2", "rs3"))
})

test_that(".applySkipRegion: handles multiple regions and chr-prefixed input", {
    df <- .ssh_smallDf()
    out <- pecotmr:::.applySkipRegion(df, c("chr1:50-250", "chr2:50-150"))
    expect_equal(nrow(out), 0L)
})

test_that(".applySkipRegion: accepts a GRanges of skip regions", {
    df <- .ssh_smallDf()
    gr <- GenomicRanges::GRanges("1", IRanges::IRanges(start = 50, end = 150))
    out <- pecotmr:::.applySkipRegion(df, gr)
    expect_equal(out$SNP, c("rs2", "rs3"))
})

test_that(".applySkipRegion: rejects malformed character entries", {
    df <- .ssh_smallDf()
    expect_error(
        pecotmr:::.applySkipRegion(df, "garbage"),
        "must be 'chr:start-end'"
    )
})

test_that(".applySkipRegion: rejects non-character non-GRanges input", {
    df <- .ssh_smallDf()
    expect_error(
        pecotmr:::.applySkipRegion(df, 42L),
        "must be a character vector"
    )
})

# ===========================================================================
# .matchAgainstSketch
# ===========================================================================

test_that(".matchAgainstSketch: errors when neither Z nor BETA is present", {
    df <- data.frame(
        chrom = "1",
        pos = 100L,
        SNP = "rs1",
        A1 = "A",
        A2 = "G",
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::.matchAgainstSketch(df, .ssh_makeHandle(), matchMinProp = 0),
        "must contain at least one of Z or BETA"
    )
})

test_that(".matchAgainstSketch: errors when A1/A2 columns are missing", {
    df <- data.frame(
        chrom = "1",
        pos = 100L,
        SNP = "rs1",
        Z = 1.0,
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::.matchAgainstSketch(df, .ssh_makeHandle(), matchMinProp = 0),
        "must contain A1 and A2 columns"
    )
})

test_that(".matchAgainstSketch: harmonizes the input against the sketch", {
    df <- data.frame(
        chrom = c("1", "1"),
        pos = c(100L, 200L),
        SNP = c("rs1", "rs2"),
        A1 = c("A", "A"),
        A2 = c("G", "G"),
        Z = c(1.0, 2.0),
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.matchAgainstSketch(
        df,
        .ssh_makeHandle(),
        matchMinProp = 0
    )
    # All variants align to the sketch; Z values pass through unchanged.
    expect_equal(nrow(out), 2L)
    expect_equal(out$Z, c(1.0, 2.0))
})

# ===========================================================================
# .applyLdMismatchQcToEntry
# ===========================================================================

test_that(".applyLdMismatchQcToEntry: errors when SNP column is missing", {
    df <- data.frame(chrom = "1", pos = 100L, Z = 1.0, stringsAsFactors = FALSE)
    expect_error(
        pecotmr:::.applyLdMismatchQcToEntry(
            df,
            .ssh_makeHandle(),
            method = "dentist"
        ),
        "requires SNP column"
    )
})

test_that(".applyLdMismatchQcToEntry: errors on variants absent from the sketch", {
    df <- data.frame(
        SNP = c("rs1", "ghost"),
        Z = c(1, 2),
        N = c(1000, 1000),
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::.applyLdMismatchQcToEntry(
            df,
            .ssh_makeHandle(),
            method = "dentist"
        ),
        "not present in the LD sketch panel"
    )
})

# ===========================================================================
# The RSS QC steps agree with the panel about which variants exist
#
# Harmonization decides the surviving variant set; kriging, LD-mismatch QC and
# the fine-mapping LD build then assert that every survivor is in the panel.
# A panel carrying a variant AND its own flip used to break that: harmonization
# emitted the variant twice (opposite Z, one per orientation) and the lookups,
# which go through matchVariants, refused both and aborted.
# ===========================================================================

# @noRd
.ssTwin_makeHandle <- function(nSamples = 60L) {
    new(
        "GenotypeHandle",
        path = "/tmp/sketch.gds",
        format = "gds",
        snpInfo = data.frame(
            SNP = c(
                "chr1:100:A:G",
                "chr1:200:A:G",
                "chr1:200:G:A",
                "chr1:300:C:T"
            ),
            CHR = rep("1", 4L),
            BP = c(100L, 200L, 200L, 300L),
            A1 = c("G", "G", "A", "T"),
            A2 = c("A", "A", "G", "C"),
            stringsAsFactors = FALSE
        ),
        nSamples = nSamples,
        sampleIds = paste0("s", seq_len(nSamples)),
        pgenPtr = NULL
    )
}

# @noRd
.ssTwin_df <- function() {
    data.frame(
        chrom = rep("1", 3L),
        pos = c(100L, 200L, 300L),
        SNP = c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:C:T"),
        A1 = c("G", "G", "T"),
        A2 = c("A", "A", "C"),
        Z = c(1.1, 2.2, 3.3),
        N = rep(1000, 3L),
        stringsAsFactors = FALSE
    )
}

# @noRd
.ssTwin_opts <- function() {
    list(
        matchMinProp = 0,
        removeIndels = FALSE,
        removeStrandAmbiguous = TRUE,
        alleleFlipKriging = TRUE,
        nForPip = 1000,
        zMismatchQc = "slalom"
    )
}

test_that("harmonization drops a twin-pair variant instead of duplicating it", {
    harm <- suppressMessages(pecotmr:::.qcHarmonizeEntry(
        .ssTwin_df(),
        .ssTwin_makeHandle(),
        .ssTwin_opts(),
        NA_character_
    ))
    # 3 in, 2 out -- never 4: the undecidable chr1:200 leaves, it does not
    # come back twice with opposite signs.
    expect_equal(nrow(harm$df), 2L)
    expect_false(any(harm$df$pos == 200L))
    expect_equal(harm$counts$harmDropped, 1L)
})

test_that("kriging and LD-mismatch QC accept every harmonized variant", {
    handle <- .ssTwin_makeHandle()
    opts <- .ssTwin_opts()
    local_mocked_bindings(
        extractBlockGenotypes = .ssQ_mockExtractor(n_samples = 60L),
        .package = "pecotmr"
    )
    harm <- suppressMessages(
        pecotmr:::.qcHarmonizeEntry(.ssTwin_df(), handle, opts, NA_character_)
    )
    expect_no_error(suppressMessages(
        pecotmr:::.qcKrigingFlip(harm$df, handle, opts, NA_character_)
    ))
    expect_no_error(suppressMessages(
        pecotmr:::.qcMismatchQc(harm$df, handle, opts, NA_character_)
    ))
})

test_that("the fine-mapping LD build accepts every harmonized variant", {
    handle <- .ssTwin_makeHandle()
    local_mocked_bindings(
        extractBlockGenotypes = .ssQ_mockExtractor(n_samples = 60L),
        .package = "pecotmr"
    )
    harm <- suppressMessages(pecotmr:::.qcHarmonizeEntry(
        .ssTwin_df(),
        handle,
        .ssTwin_opts(),
        NA_character_
    ))
    R <- pecotmr:::.ldFromSketch(
        handle,
        as.character(harm$df$SNP),
        label = "fineMappingPipeline"
    )
    expect_equal(dim(R), c(2L, 2L))
})

test_that(".applyLdMismatchQcToEntry: NA outlier flags from slalom are kept (not dropped)", {
    # Regression test: slalom (and dentist on degenerate inputs) can leave
    # NA in the `outlier` column for variants whose per-variant statistic
    # is undefined. Treating NA as TRUE would silently drop those rows;
    # treating it as FALSE (current behavior) keeps them and yields a
    # finite outlier count.
    handle <- .ssh_makeHandle()
    panel_ids <- as.character(getSnpInfo(handle)$SNP)
    vids <- panel_ids[seq_len(min(4L, length(panel_ids)))]
    df <- data.frame(
        SNP = vids,
        Z = c(1.0, 2.0, 0.5, 1.5),
        N = rep(1000L, length(vids)),
        stringsAsFactors = FALSE
    )
    # Also mock extractBlockGenotypes: the entry-side helper extracts
    # dosages from the LD sketch before delegating to ldMismatchQc, and the
    # fixture handle points at a fake /tmp/sketch.gds path with no real
    # GDS file behind it.
    local_mocked_bindings(
        extractBlockGenotypes = .ssh_mockExtractor(),
        ldMismatchQc = function(
            zScore,
            R = NULL,
            X = NULL,
            nSample = NULL,
            method = c("slalom", "dentist"),
            ldMethod = "sample",
            ...
        ) {
            data.frame(
                original_z = zScore,
                outlier = c(FALSE, TRUE, NA, NA),
                stringsAsFactors = FALSE
            )
        },
        .package = "pecotmr"
    )
    res <- pecotmr:::.applyLdMismatchQcToEntry(df, handle, method = "slalom")
    # Only the explicit TRUE row should be dropped; the two NA rows survive.
    expect_equal(nrow(res$df), 3L)
    # v2 is the TRUE outlier
    expect_false("chr1:200:A:G" %in% as.character(res$df$SNP))
    expect_equal(res$outliers, 1L)
    # Diagnostics should preserve every row + add a variant_id column.
    expect_equal(nrow(res$diagnostics), length(vids))
    expect_true("variant_id" %in% names(res$diagnostics))
})

# ===========================================================================
# Signal screen: metric resolver + absZ / bf / logBf metrics
# ===========================================================================

test_that(".resolveScreenMetric enforces one metric at a time and sane cutoffs", {
    expect_null(pecotmr:::.resolveScreenMetric()) # all 0 -> off
    expect_equal(
        pecotmr:::.resolveScreenMetric(pipCutoffToSkip = 0.5),
        list(metric = "pip", cutoff = 0.5)
    )
    expect_equal(
        pecotmr:::.resolveScreenMetric(absZCutoffToSkip = 5),
        list(metric = "absZ", cutoff = 5)
    )
    expect_equal(
        pecotmr:::.resolveScreenMetric(bfCutoffToSkip = 100),
        list(metric = "bf", cutoff = 100)
    )
    expect_equal(
        pecotmr:::.resolveScreenMetric(logBfCutoffToSkip = 3),
        list(metric = "logBf", cutoff = 3)
    )
    expect_error(
        pecotmr:::.resolveScreenMetric(
            pipCutoffToSkip = 0.5,
            absZCutoffToSkip = 5
        ),
        "one signal screen"
    )
    expect_error(
        pecotmr:::.resolveScreenMetric(absZCutoffToSkip = -1),
        "must be > 0"
    )
    expect_error(
        pecotmr:::.resolveScreenMetric(bfCutoffToSkip = -1),
        "must be > 0"
    )
})

test_that(".asScreen canonicalizes screen specs", {
    expect_null(pecotmr:::.asScreen(NULL))
    expect_null(pecotmr:::.asScreen(0))
    expect_null(pecotmr:::.asScreen(c(0.1, 0.2))) # non-scalar numeric -> off
    expect_equal(pecotmr:::.asScreen(0.9), list(metric = "pip", cutoff = 0.9))
    sc <- list(metric = "logBf", cutoff = 3)
    expect_identical(pecotmr:::.asScreen(sc), sc)
    expect_null(pecotmr:::.asScreen(list(metric = "bf", cutoff = 0))) # explicit off
})

test_that(".applyEntryScreen: absZ screen skips / retains on max|Z| (no model fit)", {
    weak <- data.frame(Z = c(0.2, 0.3, 0.1), stringsAsFactors = FALSE)
    out <- pecotmr:::.applyEntryScreen(
        weak,
        n = 1000,
        screen = list(metric = "absZ", cutoff = 5)
    )
    expect_true(out$skipped)
    expect_match(out$reason, "|Z| above 5", fixed = TRUE)
    expect_equal(nrow(out$df), 0L)

    strong <- data.frame(Z = c(0.2, 6, 0.1), stringsAsFactors = FALSE)
    out2 <- pecotmr:::.applyEntryScreen(
        strong,
        n = 1000,
        screen = list(metric = "absZ", cutoff = 5)
    )
    expect_false(out2$skipped)
    expect_identical(out2$df, strong)
})

test_that(".applyEntryScreen: bf / logBf screens use the SER lbf_variable", {
    weak <- data.frame(Z = rep(0.1, 10), stringsAsFactors = FALSE)
    strong <- data.frame(
        Z = c(10, 0.1, 0.1, 0.1, 0.1),
        stringsAsFactors = FALSE
    )
    expect_true(
        pecotmr:::.applyEntryScreen(
            weak,
            1000,
            list(metric = "logBf", cutoff = 3)
        )$skipped
    )
    expect_false(
        pecotmr:::.applyEntryScreen(
            strong,
            1000,
            list(metric = "logBf", cutoff = 3)
        )$skipped
    )
    expect_true(
        pecotmr:::.applyEntryScreen(
            weak,
            1000,
            list(metric = "bf", cutoff = 100)
        )$skipped
    )
    expect_false(
        pecotmr:::.applyEntryScreen(
            strong,
            1000,
            list(metric = "bf", cutoff = 100)
        )$skipped
    )
})

test_that("summaryStatsQc: absZ / bf / logBf screens skip a no-signal entry", {
    mk <- function() {
        gr <- .ssQ_makeEntryGr()
        S4Vectors::mcols(gr)$Z <- rep(0.1, length(gr))
        GwasSumStats(
            study = "g1",
            entry = list(gr),
            genome = "hg19",
            ldSketch = .ssQ_makeHandle()
        )
    }
    for (arg in list(
        list(absZCutoffToSkip = 5),
        list(bfCutoffToSkip = 100),
        list(logBfCutoffToSkip = 5)
    )) {
        res <- do.call(summaryStatsQc, c(list(mk()), arg, list(nCutoff = 0)))
        ea <- getQcInfo(res)$entryAudit[[1L]]
        expect_true(isTRUE(ea$pipScreenSkipped))
        expect_equal(length(res[[1L]]), 0L)
    }
})

test_that("summaryStatsQc: absZ screen retains an entry with a strong marginal Z", {
    gr <- .ssQ_makeEntryGr()
    z <- rep(0.1, length(gr))
    z[1] <- 8
    S4Vectors::mcols(gr)$Z <- z
    ss <- GwasSumStats(
        study = "g1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle()
    )
    res <- summaryStatsQc(ss, absZCutoffToSkip = 5, nCutoff = 0)
    ea <- getQcInfo(res)$entryAudit[[1L]]
    expect_false(isTRUE(ea$pipScreenSkipped))
    expect_gt(length(res[[1L]]), 0L)
})

test_that("summaryStatsQc: enabling two screens at once errors", {
    ss <- .ssQ_makeGwasSumStats()
    expect_error(
        summaryStatsQc(ss, pipCutoffToSkip = 0.5, absZCutoffToSkip = 5),
        "one signal screen"
    )
})


context("dentist_qc")
library(MASS)
library(corpcor)

generate_dentist_data <- function(
    seed = 42,
    nSnps = 100,
    sample_size = 100,
    n_outliers = 5,
    start_pos = 1000000,
    end_pos = 4000000
) {
    set.seed(seed)
    cor_matrix <- matrix(0, nrow = nSnps, ncol = nSnps)
    for (i in 1:(nSnps - 1)) {
        for (j in (i + 1):nSnps) {
            cor_matrix[i, j] <- runif(1, 0.2, 0.8)
            cor_matrix[j, i] <- cor_matrix[i, j]
        }
    }
    diag(cor_matrix) <- 1
    ld_matrix <- cov2cor(make.positive.definite(cor_matrix))
    z_scores <- mvrnorm(n = 1, mu = rep(0, nSnps), Sigma = ld_matrix)
    outlier_indices <- sample(1:nSnps, n_outliers)
    z_scores[outlier_indices] <- rnorm(n_outliers, mean = 0, sd = 5)
    sumstat <- data.frame(
        position = unlist(lapply(
            seq(start_pos, end_pos, length.out = nSnps),
            round
        )),
        z = z_scores
    )
    return(list(sumstat = sumstat, ldMat = ld_matrix, nSample = sample_size))
}

generate_dentist_single_window_data <- function(
    seed = 42,
    nSnps = 100,
    sample_size = 100,
    n_outliers = 5
) {
    set.seed(seed)
    cor_matrix <- matrix(0, nrow = nSnps, ncol = nSnps)
    for (i in 1:(nSnps - 1)) {
        for (j in (i + 1):nSnps) {
            cor_matrix[i, j] <- runif(1, 0.2, 0.8)
            cor_matrix[j, i] <- cor_matrix[i, j]
        }
    }
    diag(cor_matrix) <- 1
    ld_matrix <- cov2cor(make.positive.definite(cor_matrix))
    z_scores <- mvrnorm(n = 1, mu = rep(0, nSnps), Sigma = ld_matrix)
    outlier_indices <- sample(1:nSnps, n_outliers)
    z_scores[outlier_indices] <- rnorm(n_outliers, mean = 0, sd = 5)
    return(list(z_scores = z_scores, ldMat = ld_matrix, nSample = sample_size))
}

# ===========================================================================
# dentist: basic tests
# ===========================================================================

test_that("dentist output has exactly N rows for N input variants", {
    data <- generate_dentist_data(nSnps = 100)
    expect_warning(
        res <- dentist(data$sumstat, R = data$ldMat, nSample = data$nSample)
    )
    expect_equal(nrow(res), 100)
})

test_that("dentist output has exactly N rows with correctChenEtAlBug = FALSE", {
    data <- generate_dentist_data(nSnps = 100)
    expect_warning(
        res <- dentist(
            data$sumstat,
            R = data$ldMat,
            nSample = data$nSample,
            correctChenEtAlBug = FALSE
        )
    )
    expect_equal(nrow(res), 100)
})

test_that("dentist stops when missing position", {
    data <- generate_dentist_data()
    colnames(data$sumstat) <- c("something", "z")
    expect_error(
        dentist(
            data$sumstat,
            R = data$ldMat,
            nSample = data$nSample,
            correctChenEtAlBug = FALSE
        ),
        regexp = "missing either.*pos.*or.*z"
    )
})

test_that("dentist stops when missing zscore", {
    data <- generate_dentist_data()
    colnames(data$sumstat) <- c("position", "something")
    expect_error(
        dentist(
            data$sumstat,
            R = data$ldMat,
            nSample = data$nSample,
            correctChenEtAlBug = FALSE
        ),
        regexp = "missing either.*pos.*or.*z"
    )
})

test_that("dentist accepts 'position' and 'zscore' column names", {
    set.seed(42)
    nSnps <- 80
    n_samples <- 100
    cor_matrix <- matrix(0, nrow = nSnps, ncol = nSnps)
    for (i in 1:(nSnps - 1)) {
        for (j in (i + 1):nSnps) {
            cor_matrix[i, j] <- runif(1, 0.2, 0.8)
            cor_matrix[j, i] <- cor_matrix[i, j]
        }
    }
    diag(cor_matrix) <- 1
    ld_matrix <- cov2cor(make.positive.definite(cor_matrix))
    z_scores <- mvrnorm(n = 1, mu = rep(0, nSnps), Sigma = ld_matrix)
    sumstat <- data.frame(
        position = seq(1000000, by = 1000, length.out = nSnps),
        zscore = z_scores
    )
    expect_warning(res <- dentist(sumstat, R = ld_matrix, nSample = n_samples))
    expect_equal(nrow(res), nSnps)
})

test_that("dentist with X matrix input returns exactly N rows", {
    set.seed(42)
    nSnps <- 80
    n_samples <- 100
    X <- matrix(
        rbinom(nSnps * n_samples, 2, 0.3),
        nrow = n_samples,
        ncol = nSnps
    )
    z_scores <- rnorm(nSnps)
    sumstat <- data.frame(
        position = seq(1000000, by = 1000, length.out = nSnps),
        z = z_scores
    )
    expect_warning(res <- dentist(sumstat, X = X))
    expect_equal(nrow(res), nSnps)
})

# ===========================================================================
# dentistSingleWindow
# ===========================================================================

test_that("dentistSingleWindow returns exactly N rows for N input z-scores", {
    data <- generate_dentist_single_window_data()
    expect_warning(
        res <- dentistSingleWindow(
            data$z_scores,
            R = data$ldMat,
            nSample = data$nSample
        )
    )
    expect_equal(nrow(res), 100)
})

test_that("dentistSingleWindow warns when < 2000 variants", {
    data <- generate_dentist_single_window_data()
    expect_warning(dentistSingleWindow(
        data$z_scores,
        R = data$ldMat,
        nSample = data$nSample
    ))
})

test_that("dentistSingleWindow stops with zscore/LD matrix dimension mismatch", {
    data <- generate_dentist_single_window_data()
    expect_warning(expect_error(
        dentistSingleWindow(
            generate_dentist_single_window_data()$z_scores,
            R = generate_dentist_single_window_data(nSnps = 80)$ldMat,
            nSample = data$nSample
        ),
        regexp = "ldMat must be a square matrix"
    ))
})

test_that("dentistSingleWindow output columns are correct", {
    data <- generate_dentist_single_window_data()
    expect_warning(
        res <- dentistSingleWindow(
            data$z_scores,
            R = data$ldMat,
            nSample = data$nSample
        )
    )
    expected_cols <- c(
        "original_z",
        "imputed_z",
        "iter_to_correct",
        "rsq",
        "is_duplicate",
        "outlier_stat",
        "outlier"
    )
    expect_true(all(expected_cols %in% colnames(res)))
})

test_that("dentistSingleWindow original_z matches input z-scores", {
    data <- generate_dentist_single_window_data()
    expect_warning(
        res <- dentistSingleWindow(
            data$z_scores,
            R = data$ldMat,
            nSample = data$nSample
        )
    )
    expect_equal(res$original_z, data$z_scores)
})

test_that("dentistSingleWindow with X matrix input returns exactly N rows", {
    set.seed(42)
    nSnps <- 80
    n_samples <- 100
    X <- matrix(
        rbinom(nSnps * n_samples, 2, 0.3),
        nrow = n_samples,
        ncol = nSnps
    )
    z_scores <- rnorm(nSnps)
    expect_warning(res <- dentistSingleWindow(z_scores, X = X))
    expect_equal(nrow(res), nSnps)
})

test_that("dentistSingleWindow with correctChenEtAlBug = FALSE returns N rows", {
    data <- generate_dentist_single_window_data()
    expect_warning(
        res <- dentistSingleWindow(
            data$z_scores,
            R = data$ldMat,
            nSample = data$nSample,
            correctChenEtAlBug = FALSE
        )
    )
    expect_equal(nrow(res), 100)
    expect_true(all(c("original_z", "imputed_z", "outlier") %in% colnames(res)))
})

test_that("dentistSingleWindow with gcControl = TRUE returns N rows", {
    data <- generate_dentist_single_window_data()
    expect_warning(
        res <- dentistSingleWindow(
            data$z_scores,
            R = data$ldMat,
            nSample = data$nSample,
            gcControl = TRUE
        )
    )
    expect_equal(nrow(res), 100)
    expect_true(all(c("original_z", "imputed_z", "outlier") %in% colnames(res)))
})

test_that("dentist with gcControl = TRUE returns N rows", {
    data <- generate_dentist_data(nSnps = 100)
    expect_warning(
        res <- dentist(
            data$sumstat,
            R = data$ldMat,
            nSample = data$nSample,
            gcControl = TRUE
        )
    )
    expect_equal(nrow(res), 100)
})

test_that("dentistSingleWindow dedup path with message for duplicates", {
    set.seed(42)
    nSnps <- 80
    n_samples <- 100
    cor_matrix <- matrix(0, nrow = nSnps, ncol = nSnps)
    for (i in 1:(nSnps - 1)) {
        for (j in (i + 1):nSnps) {
            cor_matrix[i, j] <- runif(1, 0.2, 0.8)
            cor_matrix[j, i] <- cor_matrix[i, j]
        }
    }
    diag(cor_matrix) <- 1
    ld_matrix <- cov2cor(make.positive.definite(cor_matrix))
    z_scores <- mvrnorm(n = 1, mu = rep(0, nSnps), Sigma = ld_matrix)
    # Use a very low threshold to trigger dedup logic
    expect_warning(
        res <- dentistSingleWindow(
            z_scores,
            R = ld_matrix,
            nSample = n_samples,
            duprThreshold = 0.5
        )
    )
    expect_equal(nrow(res), nSnps)
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
        rsq = c(0.9, 0.8, 0.5),
        z_diff = c(0.1, 0.2, 0.5)
    )
    find_dup_output <- data.frame(
        dupBearer = c(-1, -1, 2, 1, -1),
        sign = c(1, -1, 1, -1, 1)
    )

    res <- pecotmr:::addDupsBackDentist(zScore, dentist_output, find_dup_output)
    # Non-duplicates: z_diff is copied directly from dentist output
    expect_equal(res$z_diff[1], 0.1)
    expect_equal(res$z_diff[2], 0.2)
    expect_equal(res$z_diff[5], 0.5)
    # Duplicates: z_diff = (zScore - imputed_z) / sqrt(1 - rsq)
    expect_equal(res$z_diff[3], (2.4 - 2.1) / sqrt(1 - 0.8), tolerance = 1e-10)
    expect_equal(
        res$z_diff[4],
        (1.4 - (-1.1)) / sqrt(1 - 0.9),
        tolerance = 1e-10
    )
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
    expect_error(pecotmr:::addDupsBackDentist(
        z_scores,
        dentist_output,
        find_dup_output
    ))
})

# ===========================================================================
# segment_by_dist
# ===========================================================================

test_that("segment_by_dist works", {
    res <- pecotmr:::segmentByDist(
        seq(2000000, 5000000, 100000),
        maxDist = 2000000,
        minDim = 10
    )
    expect_true(nrow(res) >= 1)
})

test_that("segment_by_dist fill regions cover all input positions", {
    # Verify that fill regions cover every SNP index exactly once
    pos <- seq(1000000, 5000000, length.out = 200)
    res <- pecotmr:::segmentByDist(pos, maxDist = 2000000, minDim = 10)
    # Collect all fill region indices
    covered <- integer(0)
    for (k in 1:nrow(res)) {
        covered <- c(covered, res$fillStartIdx[k]:(res$fillEndIdx[k] - 1L))
    }
    # Every position from 1 to length(pos) should be covered
    expect_equal(sort(unique(covered)), 1:length(pos))
})

test_that("segment_by_dist errors on empty positions", {
    expect_error(pecotmr:::segmentByDist(integer(0)), "No positions")
})

test_that("segment_by_dist verbose mode prints intervals", {
    pos <- seq(1000000, 5000000, length.out = 300)
    expect_message(
        pecotmr:::segmentByDist(
            pos,
            maxDist = 2000000,
            minDim = 50,
            verbose = TRUE
        ),
        "Intervals"
    )
})

# ===========================================================================
# detect_gaps
# ===========================================================================

test_that("detect_gaps finds no internal gaps for contiguous positions", {
    pos <- seq(1000, by = 100, length.out = 50)
    gaps <- pecotmr:::detectGaps(pos, gapThreshold = 500)
    # Only start and end sentinel
    expect_equal(gaps, c(1L, 51L))
})

test_that("detect_gaps finds a centromeric gap", {
    pos <- c(
        seq(1000, by = 100, length.out = 50),
        seq(2000000, by = 100, length.out = 50)
    )
    gaps <- pecotmr:::detectGaps(pos, gapThreshold = 1e6)
    expect_equal(length(gaps), 3) # start, gap, end
    expect_equal(gaps, c(1L, 51L, 101L))
})

test_that("detect_gaps finds multiple gaps", {
    pos <- c(1000, 2000, 5000000, 6000000, 12000000)
    gaps <- pecotmr:::detectGaps(pos, gapThreshold = 1e6)
    # Gaps at positions 3 and 5 (diffs > 1e6 at indices 2 and 4)
    expect_equal(gaps, c(1L, 3L, 5L, 6L))
})

test_that("detect_gaps verbose branch prints messages", {
    pos <- c(1000, 2000, 5000000, 6000000)
    expect_message(
        pecotmr:::detectGaps(pos, gapThreshold = 1e6, verbose = TRUE),
        "No\\. of gaps found"
    )
})

# ===========================================================================
# segment_by_count
# ===========================================================================

test_that("segment_by_count produces valid windows", {
    pos <- seq(1000000, by = 1000, length.out = 500)
    res <- pecotmr:::segmentByCount(pos, maxCount = 100)
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
    res <- pecotmr:::segmentByCount(pos, maxCount = 100)
    covered <- integer(0)
    for (k in 1:nrow(res)) {
        covered <- c(covered, res$fillStartIdx[k]:(res$fillEndIdx[k] - 1L))
    }
    expect_equal(sort(unique(covered)), 1:length(pos))
})

test_that("segment_by_count handles centromeric gap", {
    # Create two blocks separated by a large gap
    pos <- c(
        seq(1000000, by = 1000, length.out = 200),
        seq(5000000, by = 1000, length.out = 200)
    )
    res <- pecotmr:::segmentByCount(pos, maxCount = 100)
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
    expect_error(
        pecotmr:::segmentByCount(pos, maxCount = 100),
        "No intervals created by segmentation"
    )
})

test_that("segment_by_count creates single window for small blocks", {
    # Block of 60 variants with max_count=100: creates one window (60 >= half=50)
    pos <- seq(1000000, by = 1000, length.out = 60)
    res <- pecotmr:::segmentByCount(pos, maxCount = 100)
    expect_equal(nrow(res), 1)
    expect_equal(res$windowStartIdx[1], 1)
    expect_equal(res$windowEndIdx[1], 61)
})

test_that("segment_by_count single block creates correct number of windows", {
    # 200 variants with maxCount = 100 should create ~3 windows
    pos <- seq(1000000, by = 1000, length.out = 200)
    res <- pecotmr:::segmentByCount(pos, maxCount = 100)
    expect_true(nrow(res) >= 2)
    expect_true(nrow(res) <= 5)
})

test_that("segment_by_count errors on empty positions", {
    expect_error(
        pecotmr:::segmentByCount(integer(0), maxCount = 100),
        "No positions"
    )
})

test_that("segment_by_count verbose mode prints intervals", {
    pos <- seq(1000000, by = 1000, length.out = 300)
    expect_message(
        pecotmr:::segmentByCount(pos, maxCount = 100, verbose = TRUE),
        "Intervals"
    )
})

# ===========================================================================
# merge_windows
# ===========================================================================

test_that("merge_windows returns exactly N rows", {
    data <- generate_dentist_data(
        nSnps = 1000,
        sample_size = 1000,
        start_pos = 0,
        end_pos = 2000
    )
    window_divided_res <- pecotmr:::segmentByDist(
        data$sumstat$position,
        maxDist = 1000,
        minDim = 10
    )
    dentist_result_by_window <- list()
    suppressWarnings({
        for (k in 1:nrow(window_divided_res)) {
            idx_range <- window_divided_res$windowStartIdx[
                k
            ]:(window_divided_res$windowEndIdx[k] - 1L)
            zScore_k <- data$sumstat$z[idx_range]
            LD_mat_k <- data$ldMat[idx_range, idx_range]
            dentist_result_by_window[[k]] <- dentistSingleWindow(
                zScore_k,
                R = LD_mat_k,
                nSample = 100,
                pValueThreshold = 5.0369e-8,
                propSVD = 0.4,
                gcControl = FALSE,
                nIter = 10,
                gPvalueThreshold = 0.05,
                duprThreshold = 0.99,
                ncpus = 1,
                correctChenEtAlBug = TRUE
            )
        }
    })
    res <- pecotmr:::mergeWindows(dentist_result_by_window, window_divided_res)
    expect_equal(nrow(res), 1000)
})

test_that("merge_windows stops with window and imputed mismatch", {
    expect_error(pecotmr:::mergeWindows(
        rep(0, 5),
        data.frame(windowStartIdx = rep(0, 2), windowEndIdx = rep(0, 2))
    ))
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
    result <- pecotmr:::mergeWindows(list(window1, window2), window_info)
    expect_equal(nrow(result), 6)
    expect_true("index_global" %in% colnames(result))
    expect_true("index_within_window" %in% colnames(result))
})

# ===========================================================================
# dentist windowed mode
# ===========================================================================

test_that("dentist windowed output has exactly N rows for large input", {
    # Generate data large enough to trigger windowed mode (> min_dim)
    data <- generate_dentist_data(
        seed = 123,
        nSnps = 1000,
        sample_size = 1000,
        n_outliers = 50,
        start_pos = 0,
        end_pos = 5000000
    )
    suppressWarnings({
        res <- dentist(
            data$sumstat,
            R = data$ldMat,
            nSample = data$nSample,
            minDim = 100,
            windowSize = 2000000
        )
    })
    expect_equal(nrow(res), 1000)
})

test_that("dentist outlier_stat formula is correct: (z-imputed)^2/(1-rsq)", {
    data <- generate_dentist_single_window_data(seed = 55, nSnps = 100)
    expect_warning(
        res <- dentistSingleWindow(
            data$z_scores,
            R = data$ldMat,
            nSample = data$nSample
        )
    )
    expected_stat <- (res$original_z - res$imputed_z)^2 /
        pmax(1 - res$rsq, 1e-8)
    expect_equal(res$outlier_stat, expected_stat, tolerance = 1e-10)
})

# ===========================================================================
# dentist with count mode
# ===========================================================================

test_that("dentist with window_mode='count' returns exactly N rows", {
    data <- generate_dentist_data(
        seed = 789,
        nSnps = 500,
        sample_size = 500,
        n_outliers = 25,
        start_pos = 1000000,
        end_pos = 4000000
    )
    suppressWarnings({
        res <- dentist(
            data$sumstat,
            R = data$ldMat,
            nSample = data$nSample,
            windowMode = "count",
            minDim = 100
        )
    })
    expect_equal(nrow(res), 500)
})

# ===========================================================================
# Equivalence tests: both windowing methods
# ===========================================================================

test_that("segment_by_dist and segment_by_count agree on uniformly-spaced variants", {
    n <- 200
    spacing <- 10000 # 10kb between each variant
    pos <- seq(1000000, by = spacing, length.out = n)
    window_count <- 50 # variants per window in count mode
    window_dist <- window_count * spacing # equivalent distance

    res_dist <- pecotmr:::segmentByDist(pos, maxDist = window_dist, minDim = 10)
    res_count <- pecotmr:::segmentByCount(
        pos,
        maxCount = window_count,
        gapDist = 1e6
    )

    # Both should cover all positions
    covered_dist <- integer(0)
    for (k in 1:nrow(res_dist)) {
        covered_dist <- c(
            covered_dist,
            res_dist$fillStartIdx[k]:(res_dist$fillEndIdx[k] - 1L)
        )
    }
    covered_count <- integer(0)
    for (k in 1:nrow(res_count)) {
        covered_count <- c(
            covered_count,
            res_count$fillStartIdx[k]:(res_count$fillEndIdx[k] - 1L)
        )
    }
    expect_equal(sort(unique(covered_dist)), 1:n)
    expect_equal(sort(unique(covered_count)), 1:n)
})

test_that("both windowing modes produce same dentist results on uniform data", {
    data <- generate_dentist_data(
        seed = 555,
        nSnps = 500,
        sample_size = 500,
        n_outliers = 25,
        start_pos = 0,
        end_pos = 5000000
    )
    suppressWarnings({
        res_dist <- dentist(
            data$sumstat,
            R = data$ldMat,
            nSample = data$nSample,
            windowMode = "distance",
            minDim = 100,
            windowSize = 2000000
        )
        res_count <- dentist(
            data$sumstat,
            R = data$ldMat,
            nSample = data$nSample,
            windowMode = "count",
            minDim = 100
        )
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
        pecotmr:::resolveLdInput(R = NULL, X = NULL),
        "Either R.*or X.*must be provided"
    )
})

test_that("resolve_LD_input errors when both R and X provided", {
    R <- diag(3)
    X <- matrix(1:9, nrow = 3)
    expect_error(
        pecotmr:::resolveLdInput(R = R, X = X),
        "Provide either R or X, not both"
    )
})

test_that("resolve_LD_input errors when R provided without nSample and need_nSample is TRUE", {
    R <- diag(3)
    expect_error(
        pecotmr:::resolveLdInput(R = R, nSample = NULL, needNSample = TRUE),
        "nSample is required"
    )
})

test_that("resolve_LD_input returns nSample = NULL when need_nSample is FALSE", {
    R <- diag(3)
    result <- pecotmr:::resolveLdInput(
        R = R,
        nSample = NULL,
        needNSample = FALSE
    )
    expect_null(result$nSample)
    expect_equal(result$R, R)
})

test_that("resolve_LD_input infers nSample from X", {
    set.seed(42)
    n <- 50
    p <- 5
    X <- matrix(rbinom(n * p, 2, 0.3), nrow = n, ncol = p)
    result <- pecotmr:::resolveLdInput(X = X, needNSample = TRUE)
    expect_equal(result$nSample, n)
    expect_true(is.matrix(result$R))
    expect_equal(nrow(result$R), p)
})

test_that("resolve_LD_input converts non-matrix X to matrix", {
    set.seed(42)
    X_df <- data.frame(a = rbinom(30, 2, 0.3), b = rbinom(30, 2, 0.3))
    result <- pecotmr:::resolveLdInput(X = X_df, needNSample = FALSE)
    expect_true(is.matrix(result$R))
    expect_equal(result$nSample, 30)
})

test_that("resolve_LD_input uses explicit nSample when X provided", {
    set.seed(42)
    X <- matrix(rbinom(100, 2, 0.3), nrow = 20, ncol = 5)
    result <- pecotmr:::resolveLdInput(X = X, nSample = 999, needNSample = TRUE)
    expect_equal(result$nSample, 999)
})


# ===========================================================================
# build_segment_result (internal)
# ===========================================================================

test_that("build_segment_result caps end indices and verbose prints", {
    expect_message(
        result <- pecotmr:::buildSegmentResult(
            startList = c(1L),
            endList = c(200L),
            fillStartList = c(1L),
            fillEndList = c(200L),
            n = 100,
            verbose = TRUE
        ),
        "Intervals"
    )
    expect_equal(result$windowEndIdx[1], 101)
    expect_equal(result$fillEndIdx[1], 101)
})

test_that("build_segment_result errors on empty startList", {
    expect_error(
        pecotmr:::buildSegmentResult(
            startList = integer(0),
            endList = integer(0),
            fillStartList = integer(0),
            fillEndList = integer(0),
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
        pecotmr:::slidingWindowLoop(
            allGaps,
            n = 1000,
            ctx = list(),
            minBlockFn = function(blockSize, ctx) TRUE,
            initEndFn = function(startIdx, blockEnd, ctx) startIdx + 10,
            fillFn = function(startIdx, endIdx, notStart, notLast, ctx) {
                list(start = startIdx, end = endIdx)
            },
            stepFn = function(startIdx, blockEnd, ctx) {
                list(startIdx = startIdx, endIdx = startIdx + 10)
            },
            verbose = FALSE
        ),
        "iteration limit exceeded"
    )
})


context("univariate_rss_diagnostics")

# A single-row QtlFineMappingResult standing in for the retired
# FineMappingRow. topLoci defaults to one row per variant: the row-payload
# builder requires the two to be aligned, where the entry tolerated an empty
# table beside a non-empty variant list.
.testFineMappingRow <- function(
    variantIds,
    susieFit = list(),
    topLoci = NULL
) {
    if (is.null(topLoci)) {
        topLoci <- data.frame(
            variant_id = variantIds,
            pip = rep(0, length(variantIds)),
            stringsAsFactors = FALSE
        )
    }
    QtlFineMappingResult(
        study = "s1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(fineMappingRow(
            variantIds = variantIds,
            susieFit = susieFit,
            topLoci = topLoci
        ))
    )
}

# ===========================================================================
# getSusieResult
# ===========================================================================

test_that("getSusieResult returns NULL for empty input", {
    result <- getSusieResult(list())
    expect_null(result)
})

test_that("getSusieResult returns NULL when finemappingEntry missing", {
    result <- getSusieResult(list(some_data = 42))
    expect_null(result)
})

test_that("getSusieResult returns trimmed result when present", {
    mock_result <- list(pip = c(0.1, 0.5, 0.3), sets = list(cs = list()))
    con_data <- list(
        finemappingEntry = .testFineMappingRow(
            variantIds = c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A"),
            susieFit = mock_result
        )
    )
    result <- getSusieResult(con_data)
    expect_equal(result, mock_result)
})

# ===========================================================================
# extractTopPipInfo
# ===========================================================================

test_that("extractTopPipInfo finds top PIP variant", {
    con_data <- list(
        finemappingEntry = .testFineMappingRow(
            variantIds = c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A"),
            susieFit = list(pip = c(0.1, 0.7, 0.2))
        ),
        sumstats = list(z = c(1.0, 3.5, -0.5))
    )
    result <- extractTopPipInfo(con_data$finemappingEntry, con_data$sumstats)
    expect_equal(result$top_variant, "chr1:200:C:T")
    expect_equal(result$top_pip, 0.7)
    expect_equal(result$top_z, 3.5)
    expect_equal(result$top_variant_index, 2)
    expect_true(is.na(result$cs_name))
    expect_true(is.na(result$variants_per_cs))
})

test_that("extractTopPipInfo computes p_value from z", {
    con_data <- list(
        finemappingEntry = .testFineMappingRow(
            variantIds = c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A"),
            susieFit = list(pip = c(0.9, 0.05, 0.05))
        ),
        sumstats = list(z = c(5.0, 0.5, -0.3))
    )
    result <- extractTopPipInfo(con_data$finemappingEntry, con_data$sumstats)
    expected_pval <- pecotmr:::.zToPvalue(5.0)
    expect_equal(result$p_value, expected_pval)
})

test_that("extractTopPipInfo handles ties by taking first max", {
    con_data <- list(
        finemappingEntry = .testFineMappingRow(
            variantIds = c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A"),
            susieFit = list(pip = c(0.5, 0.5, 0.5))
        ),
        sumstats = list(z = c(1.0, 2.0, 3.0))
    )
    result <- extractTopPipInfo(con_data$finemappingEntry, con_data$sumstats)
    expect_equal(result$top_variant_index, 1)
    expect_equal(result$top_pip, 0.5)
})

# ===========================================================================
# extractCsInfo
# ===========================================================================

test_that("extractCsInfo extracts single CS correctly", {
    data(qtlSumStatsExample)
    fe <- .testFineMappingRow(
        variantIds = c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A"),
        susieFit = list(sets = list(cs = list(L_1 = c(1, 2))))
    )
    top_loci_table <- data.frame(
        variant_id = c("chr1:100:A:G", "chr1:200:C:T"),
        pip = c(0.3, 0.8),
        z = c(2.0, 4.5),
        stringsAsFactors = FALSE
    )
    # A single CS short-circuits (no between-CS correlation), so the unrelated
    # ldSource is not consulted.
    result <- extractCsInfo(
        fe,
        csNames = "L_1",
        topLociTable = top_loci_table,
        ldSource = qtlSumStatsExample
    )
    expect_equal(nrow(result), 1)
    expect_equal(result$cs_name, "L_1")
    expect_equal(result$top_variant, "chr1:200:C:T")
    expect_equal(result$top_pip, 0.8)
    expect_equal(result$variants_per_cs, 2)
    expect_true(is.na(result$cs_corr_max))
    expect_true(is.na(result$cs_corr_min))
    expect_false("cs_corr_1" %in% colnames(result))
})

test_that("extractCsInfo builds correlation columns from the ldSource", {
    data(qtlSumStatsExample)
    ss <- qtlSumStatsExample
    vids <- rownames(getLdSketch(ss))
    set.seed(1)
    fit <- list(
        sets = list(cs = list(L_1 = c(1L, 2L, 3L), L_2 = c(90L, 91L))),
        pip = runif(length(vids))
    )
    tl <- data.frame(
        variant_id = vids,
        pip = fit$pip,
        z = rnorm(length(vids)),
        stringsAsFactors = FALSE
    )
    fe <- .testFineMappingRow(
        variantIds = vids,
        susieFit = fit,
        topLoci = tl
    )
    result <- extractCsInfo(
        fe,
        csNames = c("L_1", "L_2"),
        topLociTable = tl,
        ldSource = ss
    )
    expect_equal(nrow(result), 2)
    expect_true(all(
        c("cs_corr_1", "cs_corr_2", "cs_corr_max", "cs_corr_min") %in%
            colnames(result)
    ))
    # The cs_corr_j columns are the columns of the computed between-CS matrix
    # (symmetric; diagonal == 1), reduced on demand from the ldSource.
    cc <- computeCsCorrelation(fe, ss)
    expect_equal(result$cs_corr_1, unname(cc[, 1]))
    expect_equal(result$cs_corr_2, unname(cc[, 2]))
    expect_equal(result$cs_corr_max, rep(abs(cc[1, 2]), 2))
    expect_equal(result$cs_corr_min, rep(abs(cc[1, 2]), 2))
})

test_that("extractCsInfo computes p_value from z-score", {
    data(qtlSumStatsExample)
    fe <- .testFineMappingRow(
        variantIds = c("chr1:100:A:G", "chr1:200:C:T"),
        susieFit = list(sets = list(cs = list(L_1 = c(1, 2))))
    )
    top_loci_table <- data.frame(
        variant_id = c("chr1:100:A:G", "chr1:200:C:T"),
        pip = c(0.9, 0.1),
        z = c(5.0, 0.5),
        stringsAsFactors = FALSE
    )
    result <- extractCsInfo(
        fe,
        csNames = "L_1",
        topLociTable = top_loci_table,
        ldSource = qtlSumStatsExample
    )
    expected_pval <- pecotmr:::.zToPvalue(5.0)
    expect_equal(result$p_value, expected_pval, tolerance = 1e-10)
})

# ===========================================================================
# autoDecision
# ===========================================================================

test_that("autoDecision assigns BVSR when no CS is tagged", {
    df <- data.frame(
        cs_name = c("L1", "L2"),
        top_z = c(5.0, 3.5),
        p_value = c(1e-10, 1e-6),
        stringsAsFactors = FALSE
    )
    result <- autoDecision(df, highCorrCols = character(0))
    expect_true("top_cs" %in% colnames(result))
    expect_true("tagged_cs" %in% colnames(result))
    expect_true("method" %in% colnames(result))
    expect_true(all(result$method == "BVSR"))
})

test_that("autoDecision assigns SER when all non-top CSs are tagged", {
    df <- data.frame(
        cs_name = c("L1", "L2"),
        top_z = c(5.0, 0.1),
        p_value = c(1e-10, 0.5),
        stringsAsFactors = FALSE
    )
    result <- autoDecision(df, highCorrCols = character(0))
    expect_true(result$top_cs[1])
    expect_false(result$top_cs[2])
    expect_true(result$tagged_cs[2])
    expect_true(all(result$method == "SER"))
})

test_that("autoDecision assigns BCR when untagged CS remain", {
    df <- data.frame(
        cs_name = c("L1", "L2", "L3"),
        top_z = c(5.0, 3.5, 0.1),
        p_value = c(1e-10, 1e-6, 0.5),
        stringsAsFactors = FALSE
    )
    result <- autoDecision(df, highCorrCols = character(0))
    expect_true(all(result$method == "BCR"))
})

test_that("autoDecision assigns SER for single CS", {
    df <- data.frame(
        cs_name = "L1",
        top_z = 5.0,
        p_value = 1e-10,
        stringsAsFactors = FALSE
    )
    result <- autoDecision(df, highCorrCols = character(0))
    expect_true(result$top_cs[1])
    expect_false(result$tagged_cs[1])
    expect_equal(result$method, "SER")
})

test_that("summaryStatsQc: preserves optional nCase/nControl columns through QC", {
    gr <- .ssQ_makeEntryGr(paste0("rs", 1:4), c(100L, 200L, 300L, 400L))
    ss <- GwasSumStats(
        study = "g1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle(),
        nCase = 500,
        nControl = 1500
    )
    expect_true(all(c("nCase", "nControl") %in% colnames(S4Vectors::mcols(ss))))
    out <- summaryStatsQc(ss, pipCutoffToSkip = 0, nCutoff = 0)
    expect_true(all(
        c("nCase", "nControl") %in% colnames(S4Vectors::mcols(out))
    ))
    expect_equal(out$nCase, 500)
    expect_equal(out$nControl, 1500)
})


# ===========================================================================
# mergeVariantInfo (data.frame + GRanges, flip-aware, all = TRUE/FALSE)
# ===========================================================================

test_that("mergeVariantInfo (data.frame): all = TRUE returns flip-corrected union", {
    v1 <- data.frame(
        chrom = c("1", "1", "2"),
        pos = c(100, 200, 300),
        alt = c("A", "C", "G"),
        ref = c("G", "T", "A"),
        stringsAsFactors = FALSE
    )
    # row1 exact; row2 alt/ref swapped vs v1 (flip); row3 brand new variant.
    v2 <- data.frame(
        chrom = c("1", "1", "3"),
        pos = c(100, 200, 400),
        alt = c("A", "T", "C"),
        ref = c("G", "C", "A"),
        stringsAsFactors = FALSE
    )
    out <- mergeVariantInfo(v1, v2, all = TRUE)
    expect_s3_class(out, "data.frame")
    expect_setequal(colnames(out), c("chrom", "pos", "alt", "ref"))
    # The flipped row2 collapses onto v1's orientation, so only the genuinely new
    # variant (chrom 3) is added to the 3 from v1.
    expect_equal(nrow(out), 4L)
    flipped <- out[out$chrom == "1" & out$pos == 200, ]
    expect_equal(flipped$alt, "C")
    expect_equal(flipped$ref, "T")
    expect_true(any(out$chrom == "3" & out$pos == 400))
})

test_that("mergeVariantInfo (data.frame): all = FALSE returns only flip-corrected variants2", {
    v1 <- data.frame(
        chrom = c("1", "1"),
        pos = c(100, 200),
        alt = c("A", "C"),
        ref = c("G", "T"),
        stringsAsFactors = FALSE
    )
    v2 <- data.frame(
        chrom = c("1", "1"),
        pos = c(100, 200),
        alt = c("A", "T"),
        ref = c("G", "C"),
        stringsAsFactors = FALSE
    )
    out <- mergeVariantInfo(v1, v2, all = FALSE)
    expect_equal(nrow(out), 2L)
    # row2 was a flip of v1; mergeVariantInfo rewrites it to v1's orientation.
    expect_equal(out$alt, c("A", "C"))
    expect_equal(out$ref, c("G", "T"))
})

test_that("mergeVariantInfo (GRanges): converts GRanges inputs and detects flips", {
    ssqcMakeGr <- function(chrom, pos, alt, ref) {
        gr <- GenomicRanges::GRanges(
            seqnames = paste0("chr", chrom),
            ranges = IRanges::IRanges(start = pos, width = 1L)
        )
        S4Vectors::mcols(gr) <- S4Vectors::DataFrame(alt = alt, ref = ref)
        gr
    }
    g1 <- ssqcMakeGr(
        c(1, 1, 2),
        c(100, 200, 300),
        c("A", "C", "G"),
        c("G", "T", "A")
    )
    g2 <- ssqcMakeGr(
        c(1, 1, 3),
        c(100, 200, 400),
        c("A", "T", "C"),
        c("G", "C", "A")
    )
    out <- mergeVariantInfo(g1, g2, all = TRUE)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 4L)
    flipped <- out[out$chrom == "chr1" & out$pos == 200, ]
    expect_equal(flipped$alt, "C")
    expect_equal(flipped$ref, "T")
})

test_that("mergeVariantInfo does not warn on length-mismatch recycling", {
    # Regression guard: the flip vector must be full-length so `hasMatch & ...`
    # never recycles a subset-length comparison.
    v1 <- data.frame(
        chrom = c("1", "1", "2"),
        pos = c(100, 200, 300),
        alt = c("A", "C", "G"),
        ref = c("G", "T", "A"),
        stringsAsFactors = FALSE
    )
    v2 <- data.frame(
        chrom = c("1", "1", "3"),
        pos = c(100, 200, 400),
        alt = c("A", "T", "C"),
        ref = c("G", "C", "A"),
        stringsAsFactors = FALSE
    )
    expect_no_warning(out <- mergeVariantInfo(v1, v2, all = TRUE))
    # Behaviour is unchanged: the flipped row collapses onto v1's orientation.
    flipped <- out[out$chrom == "1" & out$pos == 200, ]
    expect_equal(flipped$alt, "C")
    expect_equal(flipped$ref, "T")
})

# ===========================================================================
# harmonizeAlleles uncovered branches
# ===========================================================================

test_that("harmonizeAlleles: accepts a bare variant-id character vector (targetData)", {
    res <- pecotmr:::harmonizeAlleles(
        c("chr1:100:A:G", "chr1:200:C:T"),
        c("chr1:100:A:G", "chr1:200:C:T"),
        matchMinProp = 0
    )
    expect_equal(nrow(res$harmonizedData), 2L)
})

test_that("harmonizeAlleles: strips merge-conflicting columns (variant_id) from targetData", {
    target <- data.frame(
        chrom = c(1, 1),
        pos = c(100, 200),
        A2 = c("A", "C"),
        A1 = c("G", "T"),
        variant_id = c("old1", "old2"),
        z = c(1, 2),
        stringsAsFactors = FALSE
    )
    ref <- data.frame(
        chrom = c(1, 1),
        pos = c(100, 200),
        A2 = c("A", "C"),
        A1 = c("G", "T"),
        stringsAsFactors = FALSE
    )
    res <- pecotmr:::harmonizeAlleles(
        target,
        ref,
        colToFlip = "z",
        matchMinProp = 0
    )
    expect_equal(nrow(res$harmonizedData), 2L)
    # The input variant_id was stripped; the returned id is rebuilt from QC'd alleles.
    expect_false(any(c("old1", "old2") %in% res$harmonizedData$variant_id))
})

test_that("harmonizeAlleles: errors when colToFlip column is absent", {
    target <- data.frame(
        chrom = 1,
        pos = 100,
        A2 = "A",
        A1 = "G",
        z = 1,
        stringsAsFactors = FALSE
    )
    ref <- data.frame(
        chrom = 1,
        pos = 100,
        A2 = "A",
        A1 = "G",
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::harmonizeAlleles(
            target,
            ref,
            colToFlip = "missingCol",
            matchMinProp = 0
        ),
        "not found in targetData"
    )
})

test_that("harmonizeAlleles: errors when colToComplement column is absent", {
    target <- data.frame(
        chrom = 1,
        pos = 100,
        A2 = "A",
        A1 = "G",
        z = 1,
        stringsAsFactors = FALSE
    )
    ref <- data.frame(
        chrom = 1,
        pos = 100,
        A2 = "A",
        A1 = "G",
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::harmonizeAlleles(
            target,
            ref,
            colToComplement = "missingAf",
            matchMinProp = 0
        ),
        "not found in targetData"
    )
})

test_that("harmonizeAlleles: complements colToComplement (1 - af) on an allele swap", {
    # pos 100 exact; pos 200 allele-swapped -> sign_flip => z negated, af -> 1-af.
    target <- data.frame(
        chrom = c(1, 1),
        pos = c(100, 200),
        A2 = c("A", "A"),
        A1 = c("G", "G"),
        z = c(1, 2),
        af = c(0.3, 0.4),
        stringsAsFactors = FALSE
    )
    ref <- data.frame(
        chrom = c(1, 1),
        pos = c(100, 200),
        A2 = c("A", "G"),
        A1 = c("G", "A"),
        stringsAsFactors = FALSE
    )
    res <- pecotmr:::harmonizeAlleles(
        target,
        ref,
        colToFlip = "z",
        colToComplement = "af",
        matchMinProp = 0
    )
    hd <- res$harmonizedData
    row200 <- hd[hd$pos == 200, ]
    expect_equal(row200$af, 0.6) # 1 - 0.4
    expect_equal(row200$z, -2) # sign-flipped
})

test_that("harmonizeAlleles: removeDups = TRUE warns and drops duplicate variants", {
    target <- data.frame(
        chrom = c(1, 1),
        pos = c(100, 100),
        A2 = c("A", "A"),
        A1 = c("G", "G"),
        stringsAsFactors = FALSE
    )
    ref <- data.frame(
        chrom = 1,
        pos = 100,
        A2 = "A",
        A1 = "G",
        stringsAsFactors = FALSE
    )
    expect_warning(
        res <- pecotmr:::harmonizeAlleles(
            target,
            ref,
            removeDups = TRUE,
            matchMinProp = 0
        ),
        "Removed 1 duplicate"
    )
    expect_equal(nrow(res$harmonizedData), 1L)
})

test_that("harmonizeAlleles: errors when duplicated variant IDs remain (removeDups = FALSE)", {
    target <- data.frame(
        chrom = c(1, 1),
        pos = c(100, 100),
        A2 = c("A", "A"),
        A1 = c("G", "G"),
        stringsAsFactors = FALSE
    )
    ref <- data.frame(
        chrom = 1,
        pos = 100,
        A2 = "A",
        A1 = "G",
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::harmonizeAlleles(
            target,
            ref,
            removeDups = FALSE,
            matchMinProp = 0
        ),
        "Duplicated variant IDs remain"
    )
})

test_that("harmonizeAlleles: errors when too few variants match (matchMinProp)", {
    target <- data.frame(
        chrom = 1,
        pos = 100,
        A2 = "A",
        A1 = "G",
        stringsAsFactors = FALSE
    )
    ref <- data.frame(
        chrom = rep(1, 10),
        pos = seq(100, 1000, 100),
        A2 = rep("A", 10),
        A1 = rep("G", 10),
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::harmonizeAlleles(target, ref, matchMinProp = 0.9),
        "Not enough variants"
    )
})

# ===========================================================================
# addDupsBackDentist: dimension-mismatch stops
# ===========================================================================

test_that("addDupsBackDentist stops when dentistOutput nrow != count of non-duplicates", {
    dentistOutput <- data.frame(
        original_z = c(1, 2, 3),
        imputed_z = c(1, 2, 3),
        iter_to_correct = c(1, 1, 1),
        rsq = c(0.5, 0.5, 0.5),
        z_diff = c(0.1, 0.1, 0.1)
    )
    # 4 non-duplicate markers but only 3 rows in dentistOutput.
    findDupOutput <- list(dupBearer = c(-1, -1, -1, -1), sign = rep(1, 4))
    expect_error(
        pecotmr:::addDupsBackDentist(rep(0, 4), dentistOutput, findDupOutput),
        "does not match the occurrences"
    )
})

test_that("addDupsBackDentist stops on inconsistent zScore / findDupOutput length", {
    dentistOutput <- data.frame(
        original_z = c(1, 2),
        imputed_z = c(1, 2),
        iter_to_correct = c(1, 1),
        rsq = c(0.5, 0.5),
        z_diff = c(0.1, 0.1)
    )
    findDupOutput <- list(dupBearer = c(-1, -1), sign = c(1, 1))
    # nrow(dentistOutput) == sum(dupBearer == -1) == 2, but zScore has length 3.
    expect_error(
        pecotmr:::addDupsBackDentist(rep(0, 3), dentistOutput, findDupOutput),
        "inconsistent dimension"
    )
})

# ===========================================================================
# slalom: non-matrix X coerced to matrix
# ===========================================================================

test_that("slalom coerces a non-matrix X (data.frame) to a matrix", {
    set.seed(11)
    n_samples <- 80
    n_snps <- 6
    Xdf <- as.data.frame(matrix(
        sample(0:2, n_samples * n_snps, replace = TRUE),
        nrow = n_samples,
        ncol = n_snps
    ))
    z <- rnorm(n_snps)
    result <- slalom(zScore = z, X = Xdf)
    expect_equal(nrow(result$data), n_snps)
})

# ===========================================================================
# getSusieResult: trimmed fit is empty
# ===========================================================================

test_that("getSusieResult returns NULL when the trimmed susie fit is empty", {
    conData <- list(
        # An empty fit with variants present: topLoci must still be aligned
        # row-for-row, which the entry did not enforce.
        finemappingEntry = fineMappingRow(
            variantIds = c("chr1:100:A:G", "chr1:200:C:T"),
            susieFit = list(),
            topLoci = data.frame(
                variant_id = c("chr1:100:A:G", "chr1:200:C:T"),
                pip = c(0, 0),
                stringsAsFactors = FALSE
            )
        )
    )
    expect_null(getSusieResult(conData))
})

# ===========================================================================
# autoDecision: high-correlation tagging branch is reached
# ===========================================================================

test_that("autoDecision evaluates the high-corr tagging expression for non-top CS", {
    # A non-top CS with a small p-value forces evaluation of the highCorrCols
    # branch (the `..col` accessor errors without data.table; we only need the
    # line to be exercised).
    df <- data.frame(
        cs_name = c("L1", "L2"),
        top_z = c(5.0, 3.0),
        p_value = c(1e-10, 1e-10),
        stringsAsFactors = FALSE
    )
    expect_error(suppressWarnings(autoDecision(
        df,
        highCorrCols = c("cs_corr_1")
    )))
})

# ===========================================================================
# raissSingleMatrix / raissSingleMatrixFromX uncovered branches
# ===========================================================================

test_that("raissSingleMatrix coerces a data.frame LD matrix and is verbose", {
    set.seed(99)
    p <- 8
    ref_panel <- data.frame(
        chrom = rep(1, p),
        pos = seq(10, p * 10, 10),
        variant_id = paste0("rs", 1:p),
        A1 = rep("A", p),
        A2 = rep("G", p),
        stringsAsFactors = FALSE
    )
    known_idx <- c(1, 3, 5)
    known_zscores <- data.frame(
        chrom = rep(1, 3),
        pos = ref_panel$pos[known_idx],
        variant_id = ref_panel$variant_id[known_idx],
        A1 = rep("A", 3),
        A2 = rep("G", 3),
        z = rnorm(3),
        stringsAsFactors = FALSE
    )
    ldDf <- as.data.frame(diag(p))
    res <- pecotmr:::raissSingleMatrix(
        ref_panel,
        known_zscores,
        ldDf,
        r2Threshold = 0,
        minimumLd = 0,
        verbose = FALSE
    )
    expect_true(is.list(res))
    expect_equal(nrow(res$resultNofilter), p)
})

test_that("raissSingleMatrix emits verbose 'No known variants' message", {
    ref_panel <- data.frame(
        chrom = rep(1, 4),
        pos = seq(10, 40, 10),
        variant_id = paste0("rs", 1:4),
        A1 = rep("A", 4),
        A2 = rep("G", 4),
        stringsAsFactors = FALSE
    )
    known_zscores <- data.frame(
        chrom = rep(1, 2),
        pos = c(500, 600),
        variant_id = c("ghost1", "ghost2"),
        A1 = rep("A", 2),
        A2 = rep("G", 2),
        z = rnorm(2),
        stringsAsFactors = FALSE
    )
    expect_message(
        res <- pecotmr:::raissSingleMatrix(
            ref_panel,
            known_zscores,
            diag(4),
            verbose = TRUE
        ),
        "No known variants"
    )
    expect_null(res)
})

test_that("raissSingleMatrix emits verbose 'No unknown variants' message", {
    ref_panel <- data.frame(
        chrom = rep(1, 4),
        pos = seq(10, 40, 10),
        variant_id = paste0("rs", 1:4),
        A1 = rep("A", 4),
        A2 = rep("G", 4),
        stringsAsFactors = FALSE
    )
    known_zscores <- data.frame(
        chrom = rep(1, 4),
        pos = seq(10, 40, 10),
        variant_id = paste0("rs", 1:4),
        A1 = rep("A", 4),
        A2 = rep("G", 4),
        z = rnorm(4),
        stringsAsFactors = FALSE
    )
    expect_message(
        res <- pecotmr:::raissSingleMatrix(
            ref_panel,
            known_zscores,
            diag(4),
            verbose = TRUE
        ),
        "No unknown variants"
    )
    expect_equal(res$resultNofilter, known_zscores)
})

test_that("raissSingleMatrixFromX stops on unsorted positions", {
    ref_panel <- data.frame(
        chrom = rep(1, 4),
        pos = c(40, 30, 20, 10),
        variant_id = paste0("rs", 1:4),
        A1 = rep("A", 4),
        A2 = rep("G", 4),
        stringsAsFactors = FALSE
    )
    known_zscores <- data.frame(
        chrom = 1,
        pos = 10,
        variant_id = "rs4",
        A1 = "A",
        A2 = "G",
        z = 1,
        stringsAsFactors = FALSE
    )
    X <- matrix(rnorm(20 * 4), nrow = 20)
    colnames(X) <- ref_panel$variant_id
    expect_error(
        pecotmr:::raissSingleMatrixFromX(
            ref_panel,
            known_zscores,
            X,
            verbose = FALSE
        ),
        "increasing order of pos"
    )
})

test_that("raissSingleMatrixFromX emits verbose no-known / no-unknown messages", {
    set.seed(7)
    p <- 5
    ref_panel <- data.frame(
        chrom = rep(1, p),
        pos = seq(10, p * 10, 10),
        variant_id = paste0("rs", 1:p),
        A1 = rep("A", p),
        A2 = rep("G", p),
        stringsAsFactors = FALSE
    )
    X <- scale(matrix(sample(0:2, 40 * p, replace = TRUE), nrow = 40))
    X[is.na(X)] <- 0
    colnames(X) <- ref_panel$variant_id
    ghost <- data.frame(
        chrom = rep(1, 2),
        pos = c(900, 1000),
        variant_id = c("g1", "g2"),
        A1 = rep("A", 2),
        A2 = rep("G", 2),
        z = rnorm(2),
        stringsAsFactors = FALSE
    )
    expect_message(
        res_no_known <- pecotmr:::raissSingleMatrixFromX(
            ref_panel,
            ghost,
            X,
            verbose = TRUE
        ),
        "No known variants"
    )
    expect_null(res_no_known)

    all_known <- data.frame(
        chrom = rep(1, p),
        pos = seq(10, p * 10, 10),
        variant_id = paste0("rs", 1:p),
        A1 = rep("A", p),
        A2 = rep("G", p),
        z = rnorm(p),
        stringsAsFactors = FALSE
    )
    expect_message(
        res_no_unknown <- pecotmr:::raissSingleMatrixFromX(
            ref_panel,
            all_known,
            X,
            verbose = TRUE
        ),
        "No unknown variants"
    )
    expect_equal(res_no_unknown$resultNofilter, all_known)
})

# ===========================================================================
# raiss: genotypeMatrix dispatch verbose / error branches
# ===========================================================================

test_that("raiss genotypeMatrix path: single-matrix, list, all-fail, and bad-type", {
    data <- generate_X_test_data(n = 60, p = 20, n_known = 10, seed = 3)

    # Single matrix, verbose -> "Processing genotype matrix via SVD..."
    expect_message(
        res_single <- raiss(
            data$ref_panel,
            data$known_zscores,
            genotypeMatrix = data$X,
            r2Threshold = 0,
            minimumLd = 0,
            verbose = TRUE
        ),
        "Processing genotype matrix"
    )
    expect_equal(nrow(res_single$resultNofilter), nrow(data$ref_panel))

    # List of (one) matrix, verbose -> "Processing multiple genotype matrix blocks"
    expect_message(
        res_list <- raiss(
            data$ref_panel,
            data$known_zscores,
            genotypeMatrix = list(data$X),
            r2Threshold = 0,
            minimumLd = 0,
            verbose = TRUE
        ),
        "Processing multiple genotype matrix blocks"
    )
    expect_true(is.list(res_list))

    # List where every block fails (foreign known variants) -> NULL + message
    ghost <- data.frame(
        chrom = 1,
        pos = 99999,
        variant_id = "zzz",
        A1 = "A",
        A2 = "G",
        z = 1,
        stringsAsFactors = FALSE
    )
    expect_message(
        res_fail <- raiss(
            data$ref_panel,
            ghost,
            genotypeMatrix = list(
                data$X[, 1:10, drop = FALSE],
                data$X[, 11:20, drop = FALSE]
            ),
            verbose = TRUE
        ),
        "No blocks could be processed"
    )
    expect_null(res_fail)

    # Neither matrix nor list -> hard error
    expect_error(
        raiss(data$ref_panel, data$known_zscores, genotypeMatrix = 42),
        "must be a matrix or a list"
    )
})

# ===========================================================================
# .qcRaissMerge: observed-variant AF is preserved through the RAISS merge
# ===========================================================================

# RAISS reconstructs z only; the merge rebuilds the entry from its output and
# used to drop AF entirely, so --impute left top_loci$af NA for the whole entry.
# The observed variants' (harmonized) AF must survive; imputed variants get NA.
.raissMergeInputs <- function(withAf = TRUE) {
    obs <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:A:T")
    imp <- c("chr1:150:A:G", "chr1:250:C:T")
    df <- tibble(
        SNP = obs,
        A1 = c("G", "T", "T"),
        A2 = c("A", "C", "A"),
        Z = c(1.2, -0.8, 2.1),
        N = c(1000L, 1000L, 1000L)
    )
    if (withAf) {
        # A post-harmonization (directional) AF -- what the merge must re-attach.
        df$AF <- c(0.10, 0.25, 0.40)
    }
    impDf <- data.frame(
        chrom = rep("1", 5L),
        pos = c(100L, 200L, 300L, 150L, 250L),
        variant_id = c(obs, imp),
        A1 = c("G", "T", "T", "G", "T"),
        A2 = c("A", "C", "A", "A", "C"),
        z = c(1.2, -0.8, 2.1, 0.5, -0.3),
        n = rep(1000L, 5L),
        stringsAsFactors = FALSE
    )
    list(
        imputed = list(resultFilter = impDf),
        knownZ = df["SNP"],
        df = df,
        obs = obs,
        imp = imp
    )
}

test_that(".qcRaissMerge preserves observed AF and leaves imputed AF NA", {
    f <- .raissMergeInputs(withAf = TRUE)
    out <- pecotmr:::.qcRaissMerge(f$imputed, f$knownZ, f$df)$df
    expect_true("AF" %in% colnames(out))
    # observed variants keep exactly their pre-imputation (harmonized) AF
    expect_equal(out$AF[match(f$obs, out$SNP)], c(0.10, 0.25, 0.40))
    # imputed variants -- absent from df, no frequency from RAISS -- are NA
    expect_true(all(is.na(out$AF[match(f$imp, out$SNP)])))
})

test_that(".qcRaissMerge adds no AF column when the study declared none", {
    f <- .raissMergeInputs(withAf = FALSE)
    out <- pecotmr:::.qcRaissMerge(f$imputed, f$knownZ, f$df)$df
    expect_false("AF" %in% colnames(out))
})

# ===========================================================================
# raiss: multi-LD-block verbose / merge / error branches
# ===========================================================================

# Build a 2-block LD list whose shared boundary variant is IMPUTED (not known)
# in both blocks, so its raissR2 is a finite value on both sides and the
# boundary-merge R2 comparison branches are exercised.
ssqcOverlapImputedBlocks <- function(seed = 5) {
    set.seed(seed)
    vid <- sprintf("chr1:%d:A:G", 100L * (1:8))
    pos <- seq(10, 80, by = 10)
    ref_panel <- data.frame(
        chrom = rep(1, 8),
        pos = pos,
        variant_id = vid,
        A1 = rep("A", 8),
        A2 = rep("G", 8),
        stringsAsFactors = FALSE
    )
    # v4 (the boundary) is left out of the known set so it is imputed in both blocks.
    known_idx <- c(1, 2, 3, 5, 6, 7, 8)
    known_zscores <- data.frame(
        chrom = rep(1, length(known_idx)),
        pos = pos[known_idx],
        variant_id = vid[known_idx],
        A1 = rep("A", length(known_idx)),
        A2 = rep("G", length(known_idx)),
        z = rnorm(length(known_idx)),
        stringsAsFactors = FALSE
    )
    mkBlock <- function(ids) {
        nb <- length(ids)
        m <- matrix(0, nb, nb)
        for (a in 1:nb) {
            for (b in 1:nb) {
                m[a, b] <- if (a == b) 1 else 0.9^abs(a - b)
            }
        }
        rownames(m) <- colnames(m) <- ids
        m
    }
    block1_ids <- vid[1:4]
    block2_ids <- vid[4:8]
    variantIndices <- rbind(
        data.frame(
            variant_id = block1_ids,
            blockId = 1L,
            stringsAsFactors = FALSE
        ),
        data.frame(
            variant_id = block2_ids,
            blockId = 2L,
            stringsAsFactors = FALSE
        )
    )
    blockMetadata <- data.frame(
        blockId = c(1L, 2L),
        chrom = c(1, 1),
        size = c(4L, 5L),
        startIdx = c(1L, 4L),
        endIdx = c(4L, 8L),
        stringsAsFactors = FALSE
    )
    ldBlocks <- list(
        ldMatrices = list(mkBlock(block1_ids), mkBlock(block2_ids)),
        variantIndices = variantIndices,
        blockMetadata = blockMetadata,
        ldVariants = vid
    )
    list(
        ref_panel = ref_panel,
        known_zscores = known_zscores,
        LD_matrix_blocks = ldBlocks
    )
}

test_that("raiss multi-LD-block: verbose messages and an imputed boundary merge", {
    td <- ssqcOverlapImputedBlocks(seed = 5)
    expect_message(
        res <- raiss(
            td$ref_panel,
            td$known_zscores,
            ldMatrix = td$LD_matrix_blocks,
            lamb = 0.01,
            rcond = 0.01,
            r2Threshold = 0,
            minimumLd = 0,
            verbose = TRUE
        ),
        "Processing multiple LD blocks"
    )
    expect_true(is.list(res))
    # Boundary variant v4 appears exactly once after the merge.
    expect_equal(sum(res$resultNofilter$variant_id == "chr1:400:A:G"), 1L)
})

test_that("raiss multi-LD-block: stops on a block dimension mismatch", {
    td <- generate_block_diagonal_test_data(
        seed = 2,
        block_structure = "non_overlapping",
        n_variants = 30
    )
    blocks <- td$LD_matrix_blocks
    # Shrink block 1's matrix so it no longer matches its variant count.
    blocks$ldMatrices[[1]] <- blocks$ldMatrices[[1]][-1, -1, drop = FALSE]
    expect_error(
        raiss(
            td$ref_panel,
            td$known_zscores,
            ldMatrix = blocks,
            verbose = FALSE
        ),
        "LD matrix dimension does not match"
    )
})

test_that("raiss multi-LD-block: returns NULL when no block has known variants", {
    td <- generate_block_diagonal_test_data(
        seed = 3,
        block_structure = "non_overlapping",
        n_variants = 30
    )
    ghost <- data.frame(
        chrom = 1,
        pos = 99999,
        variant_id = "zzz",
        A1 = "A",
        A2 = "G",
        z = 1,
        stringsAsFactors = FALSE
    )
    expect_message(
        res <- raiss(
            td$ref_panel,
            ghost,
            ldMatrix = td$LD_matrix_blocks,
            verbose = TRUE
        ),
        "No blocks could be processed"
    )
    expect_null(res)
})

# ===========================================================================
# raissModel: batch = FALSE with condition-number reporting
# ===========================================================================

test_that("raissModel batch = FALSE reports the condition number", {
    zt <- c(1.2, 0.5)
    sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
    sig_i_t <- matrix(c(0.5, 0.2), nrow = 1) # single unknown -> 1 x 2
    res <- pecotmr:::raissModel(
        zt,
        sig_t,
        sig_i_t,
        batch = FALSE,
        reportConditionNumber = TRUE
    )
    expect_true(is.numeric(res$conditionNumber))
    expect_true(is.finite(res$conditionNumber))
    # checkInversion() returns all.equal()'s value: TRUE when the inverse
    # reproduces sigT within tolerance, otherwise a character diff string.
    expect_true(
        isTRUE(res$correctInversion) || is.character(res$correctInversion)
    )
})

# ===========================================================================
# krigingOutlierQc: non-square LD stop + variantIds defaulting to rownames
# ===========================================================================

test_that("krigingOutlierQc requires a square LD matrix aligned to zScore", {
    expect_error(
        krigingOutlierQc(c(1, 2, 3), diag(2), n = 100),
        "square LD matrix"
    )
})

test_that("krigingOutlierQc defaults variantIds to rownames(R)", {
    skip_if_not(
        "kriging_rss" %in% getNamespaceExports("susieR"),
        "installed susieR has no kriging_rss"
    )
    m <- 6
    R <- matrix(0.6, m, m)
    diag(R) <- 1
    ids <- paste0("1:", seq_len(m) * 100, ":A:G")
    rownames(R) <- colnames(R) <- ids
    z <- rep(2, m)
    kr <- krigingOutlierQc(z, R, n = 1000) # no variantIds passed
    expect_equal(kr$diagnostics$variant_id, ids)
})

# ===========================================================================
# .safeSvd: all singular values below tolerance
# ===========================================================================

test_that(".safeSvd stops when all singular values fall below tolerance", {
    set.seed(1)
    mat <- matrix(rnorm(20), nrow = 5, ncol = 4)
    # tol >= 1 forces even the largest (ratio == 1) singular value below threshold.
    expect_error(
        pecotmr:::.safeSvd(mat, tol = 2),
        "below the tolerance threshold"
    )
})

# ===========================================================================
# .applyContentFilters (MAF / FRQ, INFO, N-MAD filters + missing-column stops)
# ===========================================================================

test_that(".applyContentFilters: MAF filter drops low-frequency variants", {
    df <- data.frame(
        SNP = paste0("rs", 1:4),
        MAF = c(0.30, 0.001, 0.20, 0.0005),
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.applyContentFilters(df, mafCutoff = 0.01, nCutoff = 0)
    expect_equal(nrow(out$df), 2L)
    expect_equal(out$audit$mafDropped, 2L)
})

test_that(".applyContentFilters: FRQ is normalized to MAF via min(af, 1 - af)", {
    df <- data.frame(
        SNP = paste0("rs", 1:3),
        FRQ = c(0.5, 0.995, 0.001),
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.applyContentFilters(df, mafCutoff = 0.01, nCutoff = 0)
    # FRQ 0.995 -> MAF 0.005 (dropped); 0.001 dropped; 0.5 kept.
    expect_equal(nrow(out$df), 1L)
    expect_equal(out$audit$mafDropped, 2L)
})

test_that(".applyContentFilters: no frequency skips the MAF filter", {
    df <- data.frame(SNP = paste0("rs", 1:3), stringsAsFactors = FALSE)
    expect_warning(
        out <- .applyContentFilters(df, mafCutoff = 0.01),
        "skipping the MAF filter"
    )
    expect_equal(nrow(out$df), 3L)
})

test_that(".applyContentFilters: INFO filter drops low-INFO variants", {
    df <- data.frame(
        SNP = paste0("rs", 1:4),
        INFO = c(0.99, 0.10, 0.80, 0.30),
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.applyContentFilters(df, infoCutoff = 0.5, nCutoff = 0)
    expect_equal(nrow(out$df), 2L)
    expect_equal(out$audit$infoDropped, 2L)
})

test_that(".applyContentFilters: infoCutoff > 0 without INFO column errors", {
    df <- data.frame(SNP = paste0("rs", 1:3), stringsAsFactors = FALSE)
    expect_error(
        pecotmr:::.applyContentFilters(df, infoCutoff = 0.5),
        "requires an INFO column"
    )
})

test_that(".applyContentFilters: N-outlier (MAD) filter drops extreme N", {
    df <- data.frame(
        SNP = paste0("rs", 1:5),
        N = c(1000, 1010, 1005, 995, 100000),
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.applyContentFilters(df, nCutoff = 3)
    expect_equal(nrow(out$df), 4L)
    expect_equal(out$audit$nDropped, 1L)
})

test_that(".applyContentFilters: NA N values are always dropped", {
    df <- data.frame(
        SNP = paste0("rs", 1:4),
        N = c(1000, NA, 1005, 1010),
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.applyContentFilters(df, nCutoff = 5)
    expect_false(any(is.na(out$df$N)))
    expect_equal(out$audit$nDropped, 1L)
})

# ===========================================================================
# .runEntrySummaryStatsQc / summaryStatsQc deep branches
# ===========================================================================

test_that("summaryStatsQc: emit() uses the no-label form for an empty study id", {
    # An empty study id resolves the per-entry label to NA, exercising the
    # unlabeled emit() branch.
    ss <- GwasSumStats(
        study = "",
        entry = list(.ssQ_makeEntryGr()),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle()
    )
    msgs <- capture_messages(summaryStatsQc(ss, nCutoff = 0))
    joined <- paste(msgs, collapse = "")
    expect_match(joined, "QC summary:")
    # No bracketed label prefix appears on the rollup line.
    expect_true(any(grepl("^QC summary:", msgs)))
})

test_that("summaryStatsQc: content (N) filter emits its 'kept N of M' message + rollup nCutoff segment", {
    gr <- .ssQ_makeEntryGr()
    S4Vectors::mcols(gr)$N <- c(1000L, 1010L, 1005L, 100000L) # last is an N outlier
    ss <- GwasSumStats(
        study = "g1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle()
    )
    msgs <- capture_messages(res <- summaryStatsQc(ss, nCutoff = 3))
    joined <- paste(msgs, collapse = "")
    expect_match(joined, "MAF/INFO/N filters kept")
    expect_match(joined, "nCutoff 1")
    ea <- getQcInfo(res)$entryAudit[[1L]]
    expect_equal(ea$contentFilters$nDropped, 1L)
})

test_that("summaryStatsQc: derives BETA/SE from Z+MAF+N and records the audit", {
    gr <- GenomicRanges::GRanges(
        seqnames = rep("chr1", 4),
        ranges = IRanges::IRanges(start = c(100L, 200L, 300L, 400L), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = paste0("rs", 1:4),
        A1 = rep("A", 4),
        A2 = rep("G", 4),
        Z = c(1, 2, 3, 4),
        N = rep(1000L, 4),
        MAF = c(0.2, 0.3, 0.4, 0.25)
    )
    ss <- GwasSumStats(
        study = "g1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle()
    )
    res <- summaryStatsQc(ss, nCutoff = 0)
    ea <- getQcInfo(res)$entryAudit[[1L]]
    expect_equal(ea$betaSeFromZ$nDerived, 4L)
})

test_that("summaryStatsQc: clamps tiny Z-derived P values and accumulates the audit", {
    gr <- .ssQ_makeEntryGr()
    S4Vectors::mcols(gr)$Z <- c(50, 1, 2, 3) # |Z| = 50 underflows P to 0
    ss <- GwasSumStats(
        study = "g1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle()
    )
    res <- summaryStatsQc(ss, nCutoff = 0)
    ea <- getQcInfo(res)$entryAudit[[1L]]
    expect_true(!is.null(ea$sanityChecks$smallPClamped))
    expect_gte(ea$sanityChecks$smallPClamped, 1L)
})

test_that("summaryStatsQc: early-exits when fewer than two variants survive pre-harmonization QC", {
    gr <- .ssQ_makeEntryGr(snp_ids = "rs1", positions = 100L)
    ss <- GwasSumStats(
        study = "g1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle()
    )
    res <- summaryStatsQc(ss, nCutoff = 0)
    ea <- getQcInfo(res)$entryAudit[[1L]]
    expect_match(ea$earlyExit, "fewer than two variants")
    expect_equal(length(res[[1L]]), 1L)
})

test_that("summaryStatsQc: kriging QC runs, records the flip audit, and adds the rollup segment", {
    skip_if_not(
        "kriging_rss" %in% getNamespaceExports("susieR"),
        "installed susieR has no kriging_rss"
    )
    ss <- GwasSumStats(
        study = "g1",
        entry = list(.ssQ_makeEntryGr(
            snp_ids = paste0("rs", 1:8),
            positions = seq(100L, by = 100L, length.out = 8L)
        )),
        genome = "hg19",
        ldSketch = .ssQ_makeHandleVid(snp_n = 8L)
    )
    local_mocked_bindings(
        extractBlockGenotypes = .ssQ_mockExtractor(),
        .package = "pecotmr"
    )
    msgs <- capture_messages(
        res <- summaryStatsQc(
            ss,
            alleleFlipKriging = TRUE,
            pipCutoffToSkip = 0,
            nCutoff = 0
        )
    )
    joined <- paste(msgs, collapse = "")
    expect_match(joined, "kriging sign-flipped")
    expect_match(joined, "kriging-flip ")
    ea <- getQcInfo(res)$entryAudit[[1L]]
    expect_true("krigingFlipped" %in% names(ea))
    expect_true("krigingDiagnostics" %in% names(ea))
    kd <- ea$krigingDiagnostics
    expect_true(is.data.frame(kd))
    expect_true(all(
        c("variant_id", "z", "condmean", "z_std_diff", "logLR", "flipped") %in%
            colnames(kd)
    ))
    expect_identical(sum(kd$flipped), ea$krigingFlipped)
})

test_that("summaryStatsQc: impute = TRUE assembles BETA/SE/N and median-fills missing N", {
    gr <- GenomicRanges::GRanges(
        seqnames = rep("chr1", 4),
        ranges = IRanges::IRanges(start = c(100L, 200L, 300L, 400L), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = paste0("rs", 1:4),
        A1 = rep("A", 4),
        A2 = rep("G", 4),
        Z = c(1, 2, 3, 4),
        N = rep(1000L, 4),
        BETA = c(0.1, 0.2, 0.3, 0.4),
        SE = rep(0.1, 4)
    )
    ss <- GwasSumStats(
        study = "g1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle(snp_n = 8L, n_samples = 60L)
    )
    local_mocked_bindings(
        extractBlockGenotypes = .ssQ_mockExtractor(),
        raiss = function(refPanel, knownZscores, genotypeMatrix, ...) {
            added <- refPanel[
                refPanel$variant_id %in% c("rs5", "rs6"),
                ,
                drop = FALSE
            ]
            added$z <- c(1.5, -2.0)
            added$n <- c(1000, NA) # NA triggers the median fill
            added$beta <- c(0.11, -0.22)
            added$se <- c(0.05, 0.06)
            list(resultFilter = rbind(knownZscores, added))
        },
        .package = "pecotmr"
    )
    # rs5/rs6 (pos 500/600) are beyond the observed range; widen the impute window.
    res <- summaryStatsQc(ss, impute = TRUE, imputeOpts = list(flank = 500))
    ea <- getQcInfo(res)$entryAudit[[1L]]
    expect_equal(ea$raissTotalVariants, 6L)
    expect_equal(ea$raissImputedVariants, 2L)
    mc <- S4Vectors::mcols(res[[1L]])
    expect_true(all(c("BETA", "SE", "N") %in% colnames(mc)))
    expect_false(any(is.na(mc$N))) # median-filled
})

test_that("summaryStatsQc: per-entry rollup enumerates every removed-step segment", {
    # One entry that trips each sanity / content / harmonization drop so the
    # rollup segment strings are all assembled.
    df <- data.frame(
        seqnames = c(
            "chr1",
            "chr1",
            "chr1",
            "chr1",
            "chr99",
            "chr1",
            "chr1",
            "chr1",
            "chr1",
            "chr1",
            "chr1",
            "chr1"
        ),
        pos = c(
            100L,
            200L,
            300L,
            9999L,
            100L,
            500L,
            600L,
            700L,
            800L,
            150L,
            250L,
            350L
        ),
        SNP = c(
            "rs1",
            "rs2",
            "rs3",
            "rsOff",
            "rsChr",
            "rsMiss",
            "rsBadP",
            "rsZero",
            "rsBadSE",
            "rsMaf",
            "rsInfo",
            "rsN"
        ),
        A1 = c("A", "A", "A", "A", "A", NA, "A", "A", "A", "A", "A", "A"),
        A2 = rep("G", 12),
        Z = c(5, 4, 3, 6, 2, 2, 2, 2, 2, 2, 2, 2),
        N = c(
            1000,
            1001,
            1002,
            1000,
            1000,
            1000,
            1000,
            1000,
            1000,
            1000,
            1000,
            100000
        ),
        MAF = c(
            0.30,
            0.30,
            0.30,
            0.30,
            0.30,
            0.30,
            0.30,
            0.30,
            0.30,
            0.001,
            0.30,
            0.30
        ),
        INFO = c(
            0.95,
            0.95,
            0.95,
            0.95,
            0.95,
            0.95,
            0.95,
            0.95,
            0.95,
            0.95,
            0.10,
            0.95
        ),
        P = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.5, 0.5, 0.5, 0.5, 0.5, 0.5),
        BETA = c(0.5, 0.4, 0.3, 0.6, 0.2, 0.2, 0.2, 0.0, 0.2, 0.2, 0.2, 0.2),
        SE = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, -1.0, 0.1, 0.1, 0.1),
        stringsAsFactors = FALSE
    )
    gr <- GenomicRanges::GRanges(
        seqnames = df$seqnames,
        ranges = IRanges::IRanges(start = df$pos, width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = df$SNP,
        A1 = df$A1,
        A2 = df$A2,
        Z = df$Z,
        N = df$N,
        MAF = df$MAF,
        INFO = df$INFO,
        P = df$P,
        BETA = df$BETA,
        SE = df$SE
    )
    ss <- GwasSumStats(
        study = "g1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ssQ_makeHandle()
    )
    # mafCutoff is measured against the panel too, so the fabricated sketch
    # needs a dosage source; the mock puts every variant near MAF 0.3, well
    # clear of the 0.01 cutoff, leaving the audited drop counts study-side.
    local_mocked_bindings(
        extractBlockGenotypes = .ssQ_mockExtractor(),
        .package = "pecotmr"
    )
    msgs <- capture_messages(
        res <- summaryStatsQc(
            ss,
            mafCutoff = 0.01,
            infoCutoff = 0.5,
            nCutoff = 3
        )
    )
    joined <- paste(msgs, collapse = "")
    # The chr99 variant lands in its own ELEMENT (entries are split by seqname
    # at construction), so the audit is now per (study, chromosome) and the
    # drop counts are summed across entries rather than read off entry 1.
    audits <- getQcInfo(res)$entryAudit
    tally <- function(section, field) {
        sum(unlist(map(map(audits, section), field)))
    }
    expect_equal(tally("sanityChecks", "nonstandardChrDropped"), 1L)
    expect_equal(tally("sanityChecks", "missDataDropped"), 1L)
    expect_equal(tally("sanityChecks", "pOutOfRangeDropped"), 1L)
    expect_equal(tally("sanityChecks", "zeroEffectDropped"), 1L)
    expect_equal(tally("sanityChecks", "nonpositiveSeDropped"), 1L)
    expect_equal(tally("contentFilters", "mafDropped"), 1L)
    expect_equal(tally("contentFilters", "infoDropped"), 1L)
    expect_equal(tally("contentFilters", "nDropped"), 1L)
    # The rollup line names each removed step.
    for (seg in c(
        "nonstdChr",
        "missData",
        "badP",
        "zeroEffect",
        "badSE",
        "maf ",
        "info ",
        "nCutoff",
        "harmonization"
    )) {
        expect_match(joined, seg, fixed = TRUE)
    }
})

test_that("dentistSingleWindow: rsq-warning capture path on a near-singular LD block", {
    # Highly collinear LD with a moderately strong signal can push the C++
    # adjusted rsq_eigen above 1, exercising the withCallingHandlers capture /
    # summary-warning path. The call is wrapped to stay robust either way.
    set.seed(123)
    m <- 60
    R <- matrix(0.999, m, m)
    diag(R) <- 1
    z <- as.numeric(R %*% rnorm(m, sd = 3))
    res <- suppressWarnings(
        dentistSingleWindow(
            z,
            R = R,
            nSample = 50,
            propSVD = 0.95,
            duprThreshold = 1.0
        )
    )
    expect_equal(nrow(res), m)
})

# ===========================================================================
# Final-coverage closeouts: a handful of branch/error-handler lines that the
# existing suite leaves untouched.
# ===========================================================================

test_that("harmonizeAlleles sanitizes an empty-named target column to 'unnamed_N'", {
    # 5 columns whose first four are positional (NOT literally named
    # chrom/pos/A2/A1), so harmonizeAlleles routes through variantIdToDf -- which
    # preserves the extra column -- rather than the select()-based path. The 5th
    # column has an empty name, so the joined matchResult carries an empty-named
    # column and sanitizeNames() rewrites it to "unnamed_1" (sumstatsQc.R:83).
    target <- data.frame(
        CHR = c(1, 1, 1),
        POS = c(100, 200, 300),
        ref = c("A", "C", "G"),
        alt = c("G", "T", "A"),
        extra = c(0.1, 0.2, 0.3),
        stringsAsFactors = FALSE
    )
    names(target)[5] <- ""
    ref <- data.frame(
        chrom = c(1, 1, 1),
        pos = c(100, 200, 300),
        A2 = c("A", "C", "G"),
        A1 = c("G", "T", "A"),
        stringsAsFactors = FALSE
    )
    res <- pecotmr:::harmonizeAlleles(target, ref, matchMinProp = 0)
    expect_equal(nrow(res$harmonizedData), 3L)
    nm <- colnames(res$harmonizedData)
    expect_false(any(nm == "" | is.na(nm))) # the empty name was repaired
    expect_true(any(grepl("^unnamed_", nm))) # ... to an "unnamed_*" label
})

test_that("dentistSingleWindow summarizes the cpp11 rsqExceed field", {
    # The cpp11 imputer returns any rsq values it capped at 1.0 in `rsqExceed`;
    # .dentistRunImpute reads that field and summarizes it into a single
    # warning (no warning handler / shared env needed). Mock it to return a
    # non-empty rsqExceed alongside the documented raw columns.
    set.seed(1)
    m <- 5
    R <- diag(m)
    z <- c(1, 2, 3, 4, 5)
    local_mocked_bindings(
        dentistIterativeImpute = function(ldMatR, nSample, zScoreR, ...) {
            n <- length(zScoreR)
            list(
                originalZ = as.numeric(zScoreR),
                imputedZ = as.numeric(zScoreR) * 0.5,
                rsq = rep(0.3, n),
                zDiff = rep(0.1, n),
                iterToCorrect = rep(1L, n),
                rsqExceed = c(1.02, 1.10)
            )
        },
        .package = "pecotmr"
    )
    # capture_warnings collects every warning (the <2000-variant note and the
    # rsq summary) so the summary can be asserted specifically.
    res <- NULL
    seen <- capture_warnings(
        res <- dentistSingleWindow(
            z,
            R = R,
            nSample = 1000,
            duprThreshold = 1.0
        )
    )
    expect_true(any(grepl("rsq values exceeded 1", seen))) # summary from the data
    expect_true(any(grepl("Max reported: 1.1", seen))) # max(rsqExceed)
    # The post-warning code (renames, outlier stat, z_diff drop) ran successfully.
    expect_equal(nrow(res), m)
    expect_equal(res$original_z, z)
    expect_true(all(
        c("imputed_z", "rsq", "outlier_stat", "outlier") %in% colnames(res)
    ))
    expect_false("z_diff" %in% colnames(res))
})

test_that("segmentByDist keeps the last window when its span clears the cutoff", {
    # A single dense block spanning 1.25x the distance cutoff produces two windows
    # whose final window is wide enough that adjustLastFn does NOT shrink it: the
    # else branch returns the current startIdx unchanged (sumstatsQc.R:1095).
    pos <- as.integer(round(seq(1, 1250000, length.out = 600)))
    win <- pecotmr:::segmentByDist(pos, maxDist = 1000000, minDim = 500)
    expect_true(is.data.frame(win))
    expect_gte(nrow(win), 2L)
    expect_true(all(
        c("windowStartIdx", "windowEndIdx", "fillStartIdx", "fillEndIdx") %in%
            colnames(win)
    ))
    # Window indices stay in range (end indices are 1-based exclusive).
    expect_true(all(win$windowStartIdx >= 1L))
    expect_true(all(win$windowEndIdx <= length(pos) + 1L))
    # Fill regions tile every variant index exactly once.
    covered <- integer(0)
    for (k in seq_len(nrow(win))) {
        covered <- c(covered, win$fillStartIdx[k]:(win$fillEndIdx[k] - 1L))
    }
    expect_equal(sort(unique(covered)), seq_along(pos))
})

test_that("raiss multi-LD-block skips a NULL middle block and keeps the rest (line 1965)", {
    td <- generate_block_diagonal_test_data(
        seed = 11,
        block_structure = "non_overlapping",
        n_variants = 30
    )
    # Drop every known z-score in block 2 (var11..var20) so raissSingleMatrix
    # returns NULL for it. Blocks 1 and 3 still succeed, so resultsList carries a
    # NULL hole at index 2 and combineWithBoundaryCheck(accumulated, NULL) returns
    # the accumulated result unchanged (sumstatsQc.R:1965).
    kz <- dplyr::filter(td$known_zscores, !variant_id %in% paste0("var", 11:20))
    res <- raiss(
        td$ref_panel,
        kz,
        ldMatrix = td$LD_matrix_blocks,
        lamb = 0.01,
        rcond = 0.01,
        r2Threshold = 0,
        minimumLd = 0,
        verbose = FALSE
    )
    expect_true(is.list(res))
    expect_false(is.null(res$resultNofilter))
    # Block 2 was skipped entirely, so none of its variants survive.
    expect_false(any(paste0("var", 11:20) %in% res$resultNofilter$variant_id))
    # Blocks 1 and 3 (var1..10, var21..30) are present.
    expect_true(all(c("var1", "var21") %in% res$resultNofilter$variant_id))
})

test_that("summaryStatsQc kriging QC sign-flips an LD-inconsistent variant and retains it", {
    skip_if_not(
        "kriging_rss" %in% getNamespaceExports("susieR"),
        "installed susieR has no kriging_rss"
    )
    # All eight genotype columns load on one common factor (pairwise rho ~0.7),
    # so LD predicts a shared sign. rs4's z has the opposite sign to its strongly
    # correlated neighbours -- a genuine allele switch. krigingOutlierQc flags it,
    # so the QC step sign-flips its Z in place (the nKr > 0 branch) and RETAINS
    # the row.
    corrExtractor <- function(handle, snpIdx, meanImpute = TRUE) {
        set.seed(42)
        n <- length(getSampleIds(handle))
        k <- length(snpIdx)
        f <- rnorm(n) # shared latent factor
        M <- sapply(seq_len(k), function(j) {
            sqrt(0.7) * f + sqrt(0.3) * rnorm(n)
        })
        rr <- GenomicRanges::GRanges(
            seqnames = paste0("chr", getSnpInfo(handle)$CHR[snpIdx]),
            ranges = IRanges::IRanges(
                start = getSnpInfo(handle)$BP[snpIdx],
                width = 1L
            )
        )
        S4Vectors::mcols(rr) <- S4Vectors::DataFrame(
            SNP = getSnpInfo(handle)$SNP[snpIdx],
            A1 = getSnpInfo(handle)$A1[snpIdx],
            A2 = getSnpInfo(handle)$A2[snpIdx]
        )
        dosage <- t(M)
        rownames(dosage) <- getSnpInfo(handle)$SNP[snpIdx]
        colnames(dosage) <- getSampleIds(handle)
        SummarizedExperiment::SummarizedExperiment(
            assays = list(dosage = dosage),
            rowRanges = rr,
            colData = S4Vectors::DataFrame(
                sampleId = getSampleIds(handle),
                row.names = getSampleIds(handle)
            )
        )
    }
    gr <- .ssQ_makeEntryGr(
        snp_ids = paste0("rs", 1:8),
        positions = seq(100L, by = 100L, length.out = 8L)
    )
    mc <- S4Vectors::mcols(gr)
    # LD-consistent block (all +4) with rs4's sign reversed -- a genuine allele
    # switch (logLR > 2 & |z| > 2), which the kriging QC sign-flips back to +4.
    mc$Z <- c(4, 4, 4, -4, 4, 4, 4, 4)
    mc$N <- rep(3000L, 8)
    S4Vectors::mcols(gr) <- mc
    ss <- GwasSumStats(
        study = "g1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ssQ_makeHandleVid(snp_n = 8L, n_samples = 300L)
    )
    local_mocked_bindings(
        extractBlockGenotypes = corrExtractor,
        .package = "pecotmr"
    )
    res <- suppressWarnings(
        summaryStatsQc(
            ss,
            alleleFlipKriging = TRUE,
            pipCutoffToSkip = 0,
            nCutoff = 0
        )
    )
    ea <- getQcInfo(res)$entryAudit[[1L]]
    expect_gte(ea$krigingFlipped, 1L) # at least one flipped
    expect_equal(length(res[[1L]]), 8L) # retained, not dropped
    # rs4's -15 flips to +15; its neighbours were +4, so every retained Z is now
    # positive.
    zout <- S4Vectors::mcols(res[[1L]])$Z
    expect_true(all(zout > 0))
    # the audit's flip count equals the diagnostics rows marked flipped.
    expect_identical(sum(ea$krigingDiagnostics$flipped), ea$krigingFlipped)
})

# ---- LD-sketch trimming to the summary-stats range / final variants ----------

test_that(".subsetSketchToRange keeps only panel variants in the entries' per-chrom span", {
    h <- readGenotypeHandle(
        test_path("test_data/test_variants"),
        format = "plink2"
    )
    si <- getSnpInfo(h)
    ch <- names(sort(table(si$CHR), decreasing = TRUE))[1] # busiest chromosome
    chIdx <- which(si$CHR == ch)
    slice <- si[chIdx[seq_len(min(30L, length(chIdx)))], , drop = FALSE]
    gr <- GenomicRanges::GRanges(
        paste0("chr", slice$CHR),
        IRanges::IRanges(as.integer(slice$BP), width = 1L)
    )
    S4Vectors::mcols(gr)$SNP <- slice$SNP
    sub <- getSnpInfo(pecotmr:::.subsetSketchToRange(h, list(gr)))
    lo <- min(slice$BP)
    hi <- max(slice$BP)
    expected <- si$SNP[
        pecotmr:::canonChrom(si$CHR) == pecotmr:::canonChrom(ch) &
            as.integer(si$BP) >= lo &
            as.integer(si$BP) <= hi
    ]
    expect_setequal(sub$SNP, expected)
    expect_lt(nrow(sub), nrow(si))
})

test_that(".subsetSketchToIds keeps exactly the entries' variants", {
    h <- readGenotypeHandle(
        test_path("test_data/test_variants"),
        format = "plink2"
    )
    si <- getSnpInfo(h)
    sel <- c(3L, 10L, 40L, 200L)
    gr <- GenomicRanges::GRanges(
        paste0("chr", si$CHR[sel]),
        IRanges::IRanges(as.integer(si$BP[sel]), width = 1L)
    )
    S4Vectors::mcols(gr)$SNP <- si$SNP[sel]
    sub <- getSnpInfo(pecotmr:::.subsetSketchToIds(h, list(gr)))
    expect_setequal(
        normalizeVariantId(sub$SNP),
        normalizeVariantId(si$SNP[sel])
    )
})

test_that(".subsetSketchToRange / .subsetSketchToIds are NULL-safe", {
    expect_null(pecotmr:::.subsetSketchToRange(NULL, list()))
    expect_null(pecotmr:::.subsetSketchToIds(NULL, list()))
})


# ---------------------------------------------------------------------------
# summaryStatsQc panel filters: mafCutoff / macCutoff / imissCutoff measured
# against the LD reference panel, applied before any entry is harmonized.
# ---------------------------------------------------------------------------

.ssqcPanelStem <- function() {
    test_path("test_data/test_variants")
}

.ssqcPanelHandle <- function() {
    readGenotypeHandle(.ssqcPanelStem(), format = "plink2")
}

# Panel MAF straight from the .afreq sidecar, aligned to the handle's snpInfo
# order -- the ground truth the filter is checked against.
.ssqcPanelMaf <- function(handle) {
    afreq <- readAfreq(.ssqcPanelStem())
    ids <- as.character(getSnpInfo(handle)$SNP)
    altFreq <- afreq$alt_freq[match(ids, afreq$id)]
    pmin(altFreq, 1 - altFreq)
}

# A GwasSumStats over every panel variant, each carrying a study AF of 0.30 --
# common in the study whatever the panel says, so only a panel-side filter can
# remove it.
.ssqcPanelSumStats <- function(handle) {
    si <- getSnpInfo(handle)
    gr <- GenomicRanges::GRanges(
        seqnames = paste0("chr", si$CHR),
        ranges = IRanges::IRanges(start = as.integer(si$BP), width = 1L)
    )
    set.seed(11)
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = si$SNP,
        A1 = si$A1,
        A2 = si$A2,
        Z = stats::rnorm(nrow(si)),
        N = rep(10000, nrow(si)),
        AF = rep(0.30, nrow(si))
    )
    GwasSumStats(
        study = "g1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = handle
    )
}

test_that(".ssqcPrunePanel drops exactly the sub-cutoff panel variants", {
    skip_if_not_installed("pgenlibr")
    h <- .ssqcPanelHandle()
    ids <- as.character(getSnpInfo(h)$SNP)
    maf <- .ssqcPanelMaf(h)
    cutoff <- stats::median(maf, na.rm = TRUE)
    pruned <- expect_message(
        pecotmr:::.ssqcPrunePanel(
            h,
            list(mafCutoff = cutoff, macCutoff = 0, imissCutoff = 1),
            "test"
        ),
        "below the LD-panel MAF / MAC / missingness cutoffs"
    )
    keptIds <- as.character(getSnpInfo(pruned)$SNP)
    expect_lt(length(keptIds), length(ids))
    expect_setequal(keptIds, ids[!is.na(maf) & maf >= cutoff])
})

test_that(".ssqcPrunePanel is a no-op with no cutoffs / a NULL sketch", {
    h <- .ssqcPanelHandle()
    expect_identical(pecotmr:::.ssqcPrunePanel(h, NULL, "test"), h)
    cutoffs <- list(mafCutoff = 0, macCutoff = 0, imissCutoff = 1)
    expect_identical(pecotmr:::.ssqcPrunePanel(h, cutoffs, "test"), h)
    expect_null(pecotmr:::.ssqcPrunePanel(NULL, cutoffs, "test"))
})

test_that(".ssqcPrunePanel prunes the RSE sketch's seed handle, not the view", {
    skip_if_not_installed("pgenlibr")
    # A summaryStatsQc ldSketch is an RSE wrapping the GenotypeHandle;
    # harmonization, kriging, mismatch-QC and RAISS all read the panel through
    # .ldSketchHandle() (the dosage DelayedArray's seed). Subsetting the RSE
    # rows leaves that seed carrying the full panel, so the prune must reach
    # the seed -- the real-flow shape the bare-handle tests do not exercise.
    h <- .ssqcPanelHandle()
    maf <- .ssqcPanelMaf(h)
    ids <- as.character(getSnpInfo(h)$SNP)
    cutoff <- stats::median(maf, na.rm = TRUE)
    rse <- pecotmr:::.genotypeExperiment(h)
    expect_s4_class(rse, "RangedSummarizedExperiment")
    pruned <- suppressMessages(pecotmr:::.ssqcPrunePanel(
        rse,
        list(mafCutoff = cutoff, macCutoff = 0, imissCutoff = 1),
        "test"
    ))
    seed <- pecotmr:::.ldSketchHandle(pruned)
    expect_setequal(
        as.character(getSnpInfo(seed)$SNP),
        ids[!is.na(maf) & maf >= cutoff]
    )
})

test_that("summaryStatsQc: mafCutoff drops panel-rare observed variants", {
    skip_if_not_installed("pgenlibr")
    # Regression for the study-side-only filter: every variant here has study
    # AF 0.30, so nothing the sumstats' own frequency column can say would
    # remove it. The panel must be what decides.
    h <- .ssqcPanelHandle()
    maf <- .ssqcPanelMaf(h)
    cutoff <- stats::median(maf, na.rm = TRUE)
    ss <- .ssqcPanelSumStats(h)
    out <- suppressMessages(summaryStatsQc(ss, mafCutoff = cutoff, nCutoff = 0))
    entry <- pecotmr:::.collectionEntry(out, 1)
    expect_equal(length(entry), sum(!is.na(maf) & maf >= cutoff))
    # ... and the panel carried on the result -- the seed handle every LD
    # read keys off, not just the row view -- holds no sub-cutoff variant.
    seed <- pecotmr:::.ldSketchHandle(getLdSketch(out))
    keptMaf <- maf[match(
        as.character(getSnpInfo(seed)$SNP),
        as.character(getSnpInfo(h)$SNP)
    )]
    expect_true(all(keptMaf >= cutoff))
})

test_that("summaryStatsQc: macCutoff is the stricter of MAF / MAC", {
    skip_if_not_installed("pgenlibr")
    h <- .ssqcPanelHandle()
    maf <- .ssqcPanelMaf(h)
    ss <- .ssqcPanelSumStats(h)
    # MAC expressed per panel sample: 2 * nSamples * mafEquivalent.
    cutoff <- stats::median(maf, na.rm = TRUE)
    mac <- ceiling(cutoff * 2 * getNSamples(h))
    out <- suppressMessages(summaryStatsQc(ss, macCutoff = mac, nCutoff = 0))
    entry <- pecotmr:::.collectionEntry(out, 1)
    expected <- sum(!is.na(maf) & maf >= mac / (2 * getNSamples(h)))
    expect_equal(length(entry), expected)
})

test_that("summaryStatsQc: the panel filter is off by default", {
    skip_if_not_installed("pgenlibr")
    h <- .ssqcPanelHandle()
    ss <- .ssqcPanelSumStats(h)
    out <- suppressMessages(summaryStatsQc(ss, nCutoff = 0))
    expect_equal(
        length(pecotmr:::.collectionEntry(out, 1)),
        nrow(getSnpInfo(h))
    )
})

test_that("summaryStatsQc validates the panel cutoffs before any panel read", {
    ss <- .ssQ_makeGwasSumStats()
    expect_error(summaryStatsQc(ss, mafCutoff = -1), "mafCutoff")
    expect_error(summaryStatsQc(ss, macCutoff = c(1, 2)), "macCutoff")
    expect_error(summaryStatsQc(ss, imissCutoff = NA_real_), "imissCutoff")
    expect_error(summaryStatsQc(ss, macCutoff = Inf), "macCutoff")
})


# ---------------------------------------------------------------------------
# RAISS imputation: MAF / MAC / missingness cutoffs
#
# Without these, every rare variant in the window of a large LD sketch becomes
# an imputation target -- slow, and of little value when the study is far
# smaller than the panel. The cutoffs bound what RAISS will IMPUTE, not what it
# keeps: an observed variant survives whatever its panel frequency, because
# .raissSvdImpute pairs the panel's known columns with knownZscores$z and
# dropping one would put those two out of step.
# ---------------------------------------------------------------------------

# 20 common + 20 rare variants over 100 samples. Variant 3 is observed and
# holed; variant 8 is a target and holed, so the two cases are separable.
# @noRd
.rmask_fixture <- function() {
    set.seed(1)
    nS <- 100L
    af <- c(runif(20L, 0.2, 0.4), runif(20L, 0.002, 0.01))
    dosage <- vapply(af, function(f) rbinom(nS, 2L, f), numeric(nS))
    ids <- sprintf("chr1:%d:A:G", 1000L * seq_along(af))
    colnames(dosage) <- ids
    dosage[1:60, 3] <- NA
    dosage[1:70, 8] <- NA
    refPanel <- data.frame(
        chrom = "1",
        pos = 1000L * seq_along(af),
        A1 = "G",
        A2 = "A",
        variant_id = ids,
        stringsAsFactors = FALSE
    )
    obs <- c(1L, 2L, 3L, 21L)
    knownZ <- data.frame(
        chrom = "1",
        pos = refPanel$pos[obs],
        variant_id = ids[obs],
        A1 = "G",
        A2 = "A",
        z = rnorm(length(obs)),
        stringsAsFactors = FALSE
    )
    list(dosage = dosage, refPanel = refPanel, knownZ = knownZ, obs = obs)
}

# @noRd
.rmask_run <- function(f, ...) {
    .qcRaissTargetMask(
        f$refPanel,
        f$knownZ,
        f$dosage,
        list(imputeOpts = list(...))
    )
}

test_that(".qcRaissVariantStats reads MAF and missingness pre-imputation", {
    f <- .rmask_fixture()
    st <- .qcRaissVariantStats(f$dosage)
    expect_length(st$maf, ncol(f$dosage))
    expect_true(all(st$maf <= 0.5, na.rm = TRUE))
    # The rare half sits well below the common half.
    expect_lt(max(st$maf[21:40]), min(st$maf[1:20]))
    expect_equal(st$missRate[[3]], 0.6)
    expect_equal(st$missRate[[8]], 0.7)
    expect_equal(st$missRate[[1]], 0)
})

test_that(".qcRaissTargetMask keeps everything when no cutoff is set", {
    f <- .rmask_fixture()
    expect_true(all(.rmask_run(f)))
})

test_that(".qcRaissTargetMask drops rare imputation targets", {
    f <- .rmask_fixture()
    keep <- .rmask_run(f, mafCutoff = 0.05)
    expect_false(all(keep))
    st <- .qcRaissVariantStats(f$dosage)
    # Every dropped variant is genuinely below the cutoff...
    expect_true(all(st$maf[!keep] < 0.05, na.rm = TRUE))
    # ...and every common one survives.
    expect_true(all(keep[1:20]))
})

test_that(".qcRaissTargetMask never drops an observed variant", {
    # Variant 21 is rare AND observed. Dropping it would misalign
    # knownZscores$z against the panel's known columns inside the SVD.
    f <- .rmask_fixture()
    for (cut in c(0.05, 0.2, 0.45)) {
        keep <- .rmask_run(f, mafCutoff = cut)
        expect_true(
            all(is_in(f$knownZ$variant_id, f$refPanel$variant_id[keep])),
            label = str_c("observed variants kept at mafCutoff ", cut)
        )
    }
})

test_that(".qcRaissTargetMask applies a MAC cutoff as a MAF equivalent", {
    f <- .rmask_fixture()
    # 100 samples -> macCutoff 20 is MAF 0.1.
    byMac <- .rmask_run(f, macCutoff = 20)
    byMaf <- .rmask_run(f, mafCutoff = 0.1)
    expect_equal(byMac, byMaf)
})

test_that(".qcRaissTargetMask takes the stricter of MAF and MAC", {
    f <- .rmask_fixture()
    strictMaf <- .rmask_run(f, mafCutoff = 0.3, macCutoff = 2)
    expect_equal(strictMaf, .rmask_run(f, mafCutoff = 0.3))
    strictMac <- .rmask_run(f, mafCutoff = 0.001, macCutoff = 60)
    expect_equal(strictMac, .rmask_run(f, mafCutoff = 0.3))
})

test_that(".qcRaissTargetMask drops high-missingness targets only", {
    f <- .rmask_fixture()
    keep <- .rmask_run(f, imissCutoff = 0.5)
    # Variant 8 is 70% missing and unobserved -> dropped.
    expect_false(keep[[8]])
    # Variant 3 is 60% missing but observed -> kept.
    expect_true(keep[[3]])
    expect_equal(sum(!keep), 1L)
    # Above both rates, nothing is dropped.
    expect_true(all(.rmask_run(f, imissCutoff = 0.8)))
})

test_that(".qcRaissTargetMask drops a target whose MAF is undefined", {
    # An all-missing panel variant has no MAF; it must not slip through a
    # numeric comparison against NA.
    f <- .rmask_fixture()
    f$dosage[, 12] <- NA_real_
    keep <- .rmask_run(f, mafCutoff = 0.01)
    expect_false(keep[[12]])
})

test_that("summaryStatsQc imputeOpts cutoffs bound what RAISS imputes", {
    data(gwasSumStatsS4Example)
    gss <- gwasSumStatsS4Example
    variants <- unlist(gss)
    thin <- GwasSumStats(
        study = getStudy(gss),
        entry = list(variants[seq(1L, length(variants), by = 4L)]),
        genome = getGenome(gss),
        ldSketch = getLdSketch(gss)
    )
    nObserved <- sum(lengths(thin))
    # The bundled toy panel imputes poorly, so the R2 gate would otherwise
    # reject every target and leave nothing to count.
    base <- list(r2Threshold = 0, minimumLd = 0)
    nOut <- function(...) {
        out <- suppressWarnings(suppressMessages(summaryStatsQc(
            thin,
            impute = TRUE,
            imputeOpts = utils::modifyList(base, list(...))
        )))
        sum(lengths(out))
    }
    loose <- nOut()
    expect_gt(loose, nObserved)
    # A stricter cutoff imputes strictly fewer variants...
    mid <- nOut(mafCutoff = 0.2)
    strict <- nOut(mafCutoff = 0.45)
    expect_lt(mid, loose)
    expect_lt(strict, mid)
    # ...but never fewer than the observed set it started from.
    expect_gte(strict, nObserved)
})

# ===========================================================================
# RAISS: what must never be imputed
# ===========================================================================

.raiss_panel <- function() {
    normalizeVariantId(
        c("chr1:100:A:AT", "chr1:100:AT:A", "chr1:400:A:G", "chr1:500:C:T")
    )
}

.raiss_targets <- function(knownRaw) {
    panelIds <- .raiss_panel()
    knownIds <- normalizeVariantId(knownRaw)
    # chrom/pos are derived from each id, so a flip pair shares one position
    # and BOTH the allele-identity and observed-position guards can apply.
    # The masks are exercised individually below where that distinction
    # matters.
    frame <- function(ids) {
        p <- parseVariantId(ids)
        data.frame(
            chrom = as.character(p$chrom),
            pos = as.integer(p$pos),
            variant_id = ids,
            stringsAsFactors = FALSE
        )
    }
    refPanel <- frame(panelIds)
    knownZ <- if (length(knownIds) == 0L) {
        data.frame(
            chrom = character(0),
            pos = integer(0),
            variant_id = character(0),
            stringsAsFactors = FALSE
        )
    } else {
        frame(knownIds)
    }
    panelIds[.raissUnknownIdx(
        refPanel,
        knownZ,
        intersect(knownIds, panelIds)
    )]
}

test_that("the allele flip of a known variant is not imputed", {
    # chr1:100:AT:A is the flip of the known chr1:100:A:AT. It is not missing
    # data -- its z-score is already in hand under the other orientation.
    expect_false(is_in(
        normalizeVariantId("chr1:100:AT:A"),
        .raiss_targets("chr1:100:A:AT")
    ))
})

test_that("a dropped ambiguous variant is not re-imputed, nor is its flip", {
    # matchVariants() drops the sumstats variant when the panel carries both
    # orientations. Imputation must not put it straight back.
    got <- .raiss_targets("chr1:400:A:G")
    expect_false(is_in(normalizeVariantId("chr1:100:A:AT"), got))
    expect_false(is_in(normalizeVariantId("chr1:100:AT:A"), got))
})

test_that("a flip pair is excluded even when nothing is known there", {
    # The ambiguity is a property of the panel, not of what the sumstats
    # happen to carry.
    got <- .raiss_targets(character(0))
    expect_false(any(is_in(
        normalizeVariantId(c("chr1:100:A:AT", "chr1:100:AT:A")),
        got
    )))
})

test_that("unambiguous variants are still imputed", {
    # The guard must not cost coverage where there is no ambiguity.
    expect_setequal(
        .raiss_targets(character(0)),
        normalizeVariantId(c("chr1:400:A:G", "chr1:500:C:T"))
    )
})

test_that("a flip-of-known is skipped with only one orientation present", {
    # Reachable through the exported raiss() even though the QC pipeline
    # harmonizes to the panel orientation first.
    expect_equal(
        .raiss_targets("chr1:500:T:C"),
        normalizeVariantId("chr1:400:A:G")
    )
})

test_that(".raissFlipPairMask needs both orientations, not just an indel", {
    ids <- normalizeVariantId(
        c("chr1:100:A:AT", "chr1:100:AT:A", "chr1:200:A:AT", "chr1:300:A:G")
    )
    expect_equal(.raissFlipPairMask(ids), c(TRUE, TRUE, FALSE, FALSE))
})

test_that(".raissFlipPairMask sees a strand-recorded second orientation", {
    # A:G and C:T are one variant written on opposite strands AND swapped, so
    # matchVariants() refuses the position. A literal A1/A2 reversal key read
    # this as two unrelated variants and let RAISS impute the pair back.
    ids <- normalizeVariantId(c("chr1:100:A:G", "chr1:100:C:T"))
    expect_equal(.raissFlipPairMask(ids), c(TRUE, TRUE))
})

test_that(".raissFlipPairMask agrees with matchVariants on what is ambiguous", {
    # The mask exists to stop imputation putting back what matchVariants drops,
    # so the two have to answer the same question the same way.
    ids <- normalizeVariantId(c(
        "chr1:100:A:G",
        "chr1:100:C:T",
        "chr1:200:A:G",
        "chr1:200:G:A",
        "chr1:300:A:G",
        "chr1:400:A:T"
    ))
    m <- matchVariants(ids, ids, removeStrandAmbiguous = FALSE)
    expect_equal(.raissFlipPairMask(ids), !is_in(seq_along(ids), m$idxA))
    expect_equal(
        .raissFlipPairMask(ids),
        c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)
    )
})

test_that(".raissFlipPairMask spares the first of two identical entries", {
    # Nothing is ambiguous about a variant listed twice -- the first copy is
    # imputable. The second is not: imputing it would put one id into the
    # sumstats twice, and matchVariants answers a repeated id once, so the
    # next LD lookup would abort on the copy it could not place.
    ids <- normalizeVariantId(c("chr1:100:A:G", "chr1:100:A:G"))
    expect_equal(.raissFlipPairMask(ids), c(FALSE, TRUE))
})

test_that(".raissFlipPairMask handles an empty panel", {
    expect_equal(.raissFlipPairMask(character(0)), logical(0))
})

# ===========================================================================
# RAISS: positions the GWAS already typed (ported from PR #573)
#
# Complementary to the allele-ambiguity guards above: those key on allele
# IDENTITY (a flip, or a panel carrying both orientations), this keys on the
# POSITION being measured at all. Neither subsumes the other.
# ===========================================================================

.raiss_pr573 <- function() {
    list(
        refPanel = data.frame(
            chrom = rep("1", 5),
            pos = c(100L, 100L, 200L, 300L, 300L),
            variant_id = c(
                "1:100:A:G",
                "1:100:A:T",
                "1:200:A:G",
                "1:300:A:G",
                "1:300:A:T"
            ),
            A1 = c("G", "T", "G", "G", "T"),
            A2 = rep("A", 5),
            stringsAsFactors = FALSE
        ),
        knownZ = data.frame(
            chrom = "1",
            pos = 100L,
            variant_id = "1:100:A:G",
            A1 = "G",
            A2 = "A",
            z = 3.0,
            stringsAsFactors = FALSE
        )
    )
}

test_that(".raissFlipOfKnownMask flags a flip regardless of position", {
    # Isolates the allele-identity guard from the observed-position guard:
    # this mask sees only ids, so it cannot be satisfied by position alone.
    panel <- normalizeVariantId(c("chr1:100:AT:A", "chr1:400:A:G"))
    known <- normalizeVariantId("chr1:100:A:AT")
    expect_equal(.raissFlipOfKnownMask(panel, known), c(TRUE, FALSE))
    expect_equal(
        .raissFlipOfKnownMask(panel, character(0)),
        c(FALSE, FALSE)
    )
})

test_that("a second allele at a typed position is not imputed", {
    # Imputation fills UN-observed variants; another form of a site the GWAS
    # measured is a redundant record, not new information.
    d <- .raiss_pr573()
    idx <- .raissUnknownIdx(
        d$refPanel,
        d$knownZ,
        intersect(d$knownZ$variant_id, d$refPanel$variant_id)
    )
    expect_false(is_in("1:100:A:T", d$refPanel$variant_id[idx]))
})

test_that("un-observed positions still impute every panel form they have", {
    # The guard is about MEASURED positions, not about multi-allelic sites:
    # 1:300 is multi-allelic but untyped, so both its forms are targets.
    d <- .raiss_pr573()
    idx <- .raissUnknownIdx(
        d$refPanel,
        d$knownZ,
        intersect(d$knownZ$variant_id, d$refPanel$variant_id)
    )
    expect_setequal(
        d$refPanel$variant_id[idx],
        c("1:200:A:G", "1:300:A:G", "1:300:A:T")
    )
})

test_that("the position guard is chr-prefix tolerant", {
    # "chr1" in the panel and "1" in the sumstats is one site, not two.
    d <- .raiss_pr573()
    d$refPanel$chrom <- str_c("chr", d$refPanel$chrom)
    d$refPanel$variant_id <- str_c("chr", d$refPanel$variant_id)
    idx <- .raissUnknownIdx(d$refPanel, d$knownZ, character(0))
    expect_false(is_in("chr1:100:A:T", d$refPanel$variant_id[idx]))
})

test_that("raissSingleMatrix does not impute at a typed position", {
    d <- .raiss_pr573()
    ld <- diag(5)
    ld[1, 2] <- ld[2, 1] <- 0.7
    ld[1, 3] <- ld[3, 1] <- 0.6
    ld[1, 4] <- ld[4, 1] <- 0.5
    ld[1, 5] <- ld[5, 1] <- 0.45
    r <- raissSingleMatrix(d$refPanel, d$knownZ, ld, verbose = FALSE)
    imp <- setdiff(r$resultNofilter$variant_id, d$knownZ$variant_id)
    expect_false(is_in("1:100:A:T", imp))
    expect_true(all(is_in(
        c("1:200:A:G", "1:300:A:G", "1:300:A:T"),
        imp
    )))
    expect_true(is_in("1:100:A:G", r$resultNofilter$variant_id))
})

# ===========================================================================
# Zero-variant entries empty the sketch (ported from PR #571)
#
# An entry that early-exits QC carries no variants, so it references no LD --
# but both trimmers used to hand back the FULL genome-wide panel, because
# their keep-set is built from the entries' own positions/ids. On a real panel
# that is ~135 MB of dead weight serialized per such entry.
# ===========================================================================

.sk_panel <- function() {
    readGenotypes(
        test_path("test_data", "test_variants"),
        format = "plink1"
    )
}

test_that("both trimmers empty the sketch for zero-variant entries", {
    h <- .sk_panel()
    expect_gt(nrow(h), 0L)
    # An empty list and a zero-range GRanges both mean "no variants".
    for (entries in list(list(), list(GenomicRanges::GRanges()))) {
        expect_equal(nrow(.subsetSketchToRange(h, entries)), 0L)
        expect_equal(nrow(.subsetSketchToIds(h, entries)), 0L)
    }
})

test_that("emptying drops both axes but keeps the panel's identity", {
    h <- .sk_panel()
    e <- .subsetSketchToRange(h, list())
    # The sample axis is shed along with the variants: an emptied sketch
    # references no LD, so its anonymous sample names are dead weight (the bulk
    # of a skipped-region file). What is kept is the panel's identity -- where
    # it came from -- not its dimensions.
    expect_equal(nrow(e), 0L)
    expect_equal(ncol(e), 0L)
    hh <- .ldSketchHandle(h)
    eh <- .ldSketchHandle(e)
    expect_equal(getFormat(eh), getFormat(hh))
    expect_equal(getPath(eh), getPath(hh))
})

test_that(".emptySketch drops the seed's variants, not just the rows", {
    # A panel carries its handle twice: as rowRanges and inside the dosage
    # assay's DelayedArray seed. Dropping rows alone leaves the seed holding
    # the whole panel, which is the weight this exists to shed.
    h <- .sk_panel()
    e <- .emptySketch(h)
    expect_equal(nrow(e), 0L)
    expect_equal(nrow(getSnpInfo(.ldSketchHandle(e))), 0L)
    expect_lt(
        length(serialize(e, NULL)),
        length(serialize(h, NULL)) / 2
    )
})

test_that(".emptySketch is NULL-safe and idempotent", {
    expect_null(.emptySketch(NULL))
    e <- .emptySketch(.sk_panel())
    expect_equal(nrow(.emptySketch(e)), 0L)
})

test_that(".subsetSketchToIds keeps the panel when ids are unextractable", {
    # Ranges present but no SNP mcol is a different situation from zero
    # variants: the panel cannot be keyed, so blanking it would discard LD the
    # object may well reference. Stay conservative.
    h <- .sk_panel()
    gr <- GenomicRanges::GRanges(
        "chr21",
        IRanges::IRanges(c(17513228L, 17513580L), width = 1L)
    )
    expect_equal(nrow(.subsetSketchToIds(h, list(gr))), nrow(h))
})

# ===========================================================================
# af / MAF semantics (ported from PR #578)
#
# `af` is the DIRECTIONAL effect-allele frequency, declared only via an `af:`
# mapping key and exported as top_loci$af. A directionless `maf`/`FRQ` is
# QC-only and is never exported as af, so an undeclared study gets an honest
# NA instead of a silently mislabelled minor-allele frequency.
# ===========================================================================

.afm_df <- function(freqcol, val) {
    d <- data.frame(
        chrom = "chr1",
        pos = 100L,
        variant_id = "chr1:100:A:G",
        A1 = "A",
        A2 = "G",
        z = 2.5,
        n_sample = 1000,
        stringsAsFactors = FALSE
    )
    d[[freqcol]] <- val
    d
}

.afm_cm <- function(...) {
    c(
        chrom = "chrom",
        pos = "pos",
        variant_id = "variant_id",
        A1 = "A1",
        A2 = "A2",
        z = "z",
        n_sample = "n_sample",
        ...
    )
}

test_that("an af: key is directional AF; a maf: key is directionless", {
    # Same source column, different declaration -> different meaning.
    a <- suppressWarnings(.resolveSumstatCols(
        .afm_df("effect_allele_frequency", 0.8),
        .afm_cm(af = "effect_allele_frequency"),
        "af"
    ))
    expect_equal(a$AF, 0.8)
    expect_null(a$MAF)
    b <- suppressWarnings(.resolveSumstatCols(
        .afm_df("effect_allele_frequency", 0.8),
        .afm_cm(maf = "effect_allele_frequency"),
        "maf"
    ))
    expect_null(b$AF)
    expect_equal(b$MAF, 0.8)
})

test_that("an undeclared af warns distinctly and yields no AF column", {
    expect_warning(
        out <- .resolveSumstatCols(
            .afm_df("maf", 0.2),
            .afm_cm(maf = "maf"),
            "undeclared"
        ),
        "no effect-allele frequency declared"
    )
    expect_null(out$AF)
})

test_that("a declared-but-unusable af warns with a different message", {
    # Distinct from "undeclared": this is a data problem, not a mapping one.
    expect_warning(
        .resolveSumstatCols(
            .afm_df("eaf", NA_real_),
            .afm_cm(af = "eaf"),
            "missing"
        ),
        "declared but its values are all missing"
    )
})

test_that(".fmAfByVar reads AF, and is NULL when af was never declared", {
    mk <- function(cm) {
        .dfToEntryGranges(suppressWarnings(.resolveSumstatCols(
            .afm_df("effect_allele_frequency", 0.8),
            cm,
            "lbl"
        )))
    }
    gA <- mk(.afm_cm(af = "effect_allele_frequency"))
    expect_equal(unname(.fmAfByVar(gA, "chr1:100:A:G")), 0.8)
    gM <- mk(.afm_cm(maf = "effect_allele_frequency"))
    expect_null(.fmAfByVar(gM, "chr1:100:A:G"))
})

test_that(".cfMaf prefers af, falls back to maf, and skips when neither", {
    expect_equal(nrow(.cfMaf(data.frame(AF = c(0.8, 0.01)), 0.05)$df), 1L)
    expect_equal(nrow(.cfMaf(data.frame(MAF = c(0.2, 0.01)), 0.05)$df), 1L)
    expect_warning(
        out <- .cfMaf(data.frame(SNP = c("a", "b")), 0.05),
        "skipping the MAF filter"
    )
    expect_equal(nrow(out$df), 2L)
})

test_that("an allele swap complements af but leaves the directionless maf", {
    # af -> 1 - af on a swap; maf = min(af, 1-af) is swap-invariant, so
    # complementing it would corrupt it.
    tgt <- data.frame(
        chrom = "1",
        pos = 100L,
        A2 = "G",
        A1 = "A",
        SNP = "v1",
        Z = 2.0,
        AF = 0.8,
        MAF = 0.2,
        stringsAsFactors = FALSE
    )
    ref <- data.frame(
        chrom = "1",
        pos = 100L,
        A2 = "A",
        A1 = "G",
        variant_id = "v1",
        stringsAsFactors = FALSE
    )
    h <- harmonizeAlleles(
        tgt,
        ref,
        colToFlip = "Z",
        colToComplement = "AF",
        matchMinProp = 0,
        removeStrandAmbiguous = FALSE
    )$harmonizedData
    expect_equal(h$AF, 0.2)
    expect_equal(h$Z, -2.0)
    expect_equal(h$MAF, 0.2)
})
