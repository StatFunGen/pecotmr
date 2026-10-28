context("sumstats_qc")

.gsd_makeHandle <- function(snp_n = 2L) {
    new(
        "GenotypeHandle",
        path = "/tmp/sketch",
        format = "plink2",
        snpInfo = data.frame(
            SNP = paste0("rs", seq_len(snp_n)),
            CHR = rep("chr1", snp_n),
            BP = seq(100L, by = 100L, length.out = snp_n),
            A1 = rep("A", snp_n),
            A2 = rep("G", snp_n),
            stringsAsFactors = FALSE
        ),
        nSamples = 10L,
        sampleIds = paste0("s", seq_len(10L)),
        pgenPtr = NULL
    )
}

.gsd_makeGwasSumStats <- function(mcolsDf, snp_n = nrow(mcolsDf)) {
    gr <- GenomicRanges::GRanges(
        seqnames = rep("chr1", snp_n),
        ranges = IRanges::IRanges(
            start = seq(100L, by = 100L, length.out = snp_n),
            width = 1L
        )
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(mcolsDf)
    GwasSumStats(
        studyName = "G1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .gsd_makeHandle(snp_n),
        qcOptions = list(prebuilt = "synthetic")
    )
}

# ===========================================================================
# ldMismatchQc
# ===========================================================================

# Captures the options slalom was called with, so a test can assert that a
# constructor's settings actually reached the engine.
.ssq_captureSlalom <- function(capture) {
    function(
        zScore,
        R = NULL,
        X = NULL,
        standardError = rep(1, length(zScore)),
        abfPriorVariance = 0.04,
        nlog10pDentistSThreshold = 4.0,
        r2Threshold = 0.6,
        leadVariantChoice = "pvalue",
        ldMethod = "sample"
    ) {
        capture$r2Threshold <- r2Threshold
        capture$abfPriorVariance <- abfPriorVariance
        list(data = data.frame(outlier = rep(FALSE, length(zScore))))
    }
}

test_that(".entryToSumstatDf: derives maf from a directional af", {
    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(100L * (1:3), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = sprintf("chr1:%d:A:G", 100L * (1:3)),
        A1 = "A",
        A2 = "G",
        Z = c(2, -1.5, 1.8),
        N = 1000L,
        AF = c(0.3, 0.8, 0.5)
    )
    df <- pecotmr:::.entryToSumstatDf(gr)
    # af is directional, maf is not: 0.8 folds to 0.2.
    expect_equal(df$af, c(0.3, 0.8, 0.5))
    expect_equal(df$maf, c(0.3, 0.2, 0.5))
    # ... which is what takes the MR ratio off the z-scale fallback.
    expect_true(pecotmr:::.cipGwasHasScale(df$maf, df$N))
})

test_that(".entryToSumstatDf: a declared MAF is not overwritten by af", {
    gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100L, width = 1L))
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = "chr1:100:A:G",
        A1 = "A",
        A2 = "G",
        Z = 2,
        N = 1000L,
        AF = 0.8,
        MAF = 0.45
    )
    df <- pecotmr:::.entryToSumstatDf(gr)
    # The study asserted 0.45; af is only another route to the same thing.
    expect_equal(df$maf, 0.45)
    expect_equal(df$af, 0.8)
})

test_that(".sumStatsDf flattens the tuple sumStats() selects", {
    data(qtlSumStatsExample, package = "pecotmr", envir = environment())
    qss <- qtlSumStatsExample
    st <- as.character(qss$study[[1L]])
    ctx <- as.character(qss$context[[1L]])
    tr <- as.character(qss$trait[[1L]])
    expect_identical(
        .sumStatsDf(qss, studyName = st, context = ctx, trait = tr),
        .entryToSumstatDf(
            sumStats(qss, studyName = st, context = ctx, trait = tr)
        )
    )
    expect_error(
        .sumStatsDf(
            qss,
            studyName = st,
            context = ctx,
            trait = tr,
            require = "NOPE"
        ),
        str_c("QtlSumStats[", st, "/", ctx, "/", tr, "]: entry has no NOPE"),
        fixed = TRUE
    )
})

test_that(".sumStatsDf names a GWAS by study, an unset selector <auto>", {
    data(
        gwasSumStatsS4Example,
        qtlSumStatsExample,
        package = "pecotmr",
        envir = environment()
    )
    gss <- gwasSumStatsS4Example
    st <- as.character(gss$study[[1L]])
    expect_error(
        .sumStatsDf(gss, studyName = st, require = "NOPE"),
        str_c("GwasSumStats[", st, "]: entry has no NOPE"),
        fixed = TRUE
    )
    expect_error(
        .sumStatsDf(qtlSumStatsExample, require = "NOPE"),
        "QtlSumStats[<auto>/<auto>/<auto>]: entry has no NOPE",
        fixed = TRUE
    )
})

test_that("ldMismatchQc dentist returns a data frame with an outlier column", {
    set.seed(42)
    p <- 20
    R <- diag(p)
    z <- rnorm(p)
    expect_warning(
        result <- ldMismatchQc(z, R = R, nSample = 1000, method = "dentist"),
        "below 2000"
    )
    expect_true(is.data.frame(result) || is.list(result))
    expect_true("outlier" %in% names(result))
})

test_that("ldMismatchQc slalom returns a data frame with an outlier column", {
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
    expect_error(
        ldMismatchQc(z, R = R, method = "invalid"),
        "`method`: unknown engine 'invalid'. Known: slalom, dentist.",
        fixed = TRUE
    )
})

# ===========================================================================
# raiss and its helpers
# ===========================================================================

context("raiss")
library(tidyverse)
library(MASS)

generateDummyData <- function(
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
    input_data <- generateDummyData()
    input_data_ref_panel_unordered <- generateDummyData(
        ref_panel_ordered = FALSE
    )
    input_data_zscores_unordered <- generateDummyData(
        known_zscores_ordered = FALSE
    )
    expect_error(
        raiss(
            input_data_ref_panel_unordered$ref_panel,
            input_data$known_zscores,
            input_data$ldMatrix
        ),
        "refPanel and knownZscores must be in increasing order of pos"
    )
    expect_error(
        raiss(
            input_data$ref_panel,
            input_data_zscores_unordered$known_zscores,
            input_data$ldMatrix
        ),
        "refPanel and knownZscores must be in increasing order of pos"
    )
})

test_that("raiss refuses a plain list of block matrices", {
    # Such a list does not say which variants each block holds; it used to
    # come back as NULL with nothing said.
    input_data <- generateDummyData()
    R <- input_data$ldMatrix
    expect_error(
        raiss(input_data$ref_panel, input_data$known_zscores, list(R, R)),
        "square LD matrix or an LdMatrix"
    )
})

test_that("raiss needs refPanel with a plain LD matrix", {
    input_data <- generateDummyData()
    expect_error(
        raiss(
            knownZscores = input_data$known_zscores,
            ldMatrix = input_data$ldMatrix
        ),
        "`refPanel` is required"
    )
})

test_that("raiss takes an LdMatrix and uses its variants", {
    set.seed(3)
    p <- 30L
    R <- 0.6^abs(outer(seq_len(p), seq_len(p), "-"))
    pos <- seq(100L, by = 100L, length.out = p)
    ref_panel <- data.frame(
        chrom = "1",
        pos = pos,
        variant_id = sprintf("chr1:%d:A:G", pos),
        A1 = "G",
        A2 = "A"
    )
    known <- ref_panel[seq_len(p) %% 3L != 0L, ]
    known$z <- rnorm(nrow(known))
    full <- raiss(
        ref_panel,
        known,
        ldMatrix = R,
        r2Threshold = 0.3,
        minimumLd = 1,
        verbose = FALSE
    )
    ld <- LdMatrix(
        R,
        variants = data.frame(
            chrom = "chr1",
            pos = pos,
            SNP = ref_panel$variant_id,
            A1 = "G",
            A2 = "A"
        )
    )
    # No refPanel: the LdMatrix's own variants supply it. The z-scores'
    # chromosomes are read in UCSC style, so "1" and "chr1" both meet them;
    # z-scores on a chromosome the reference lacks have nothing to impute
    # from.
    expect_error(
        raiss(
            knownZscores = transform(known, chrom = "2"),
            ldMatrix = ld,
            verbose = FALSE
        ),
        "shares no chromosome with the LD reference"
    )
    fromLd <- raiss(
        knownZscores = known,
        ldMatrix = ld,
        r2Threshold = 0.3,
        minimumLd = 1,
        verbose = FALSE
    )
    expect_equal(
        arrange(fromLd$resultNofilter, variant_id)$z,
        arrange(full$resultNofilter, variant_id)$z
    )
    fromChr <- raiss(
        knownZscores = transform(known, chrom = "chr1"),
        ldMatrix = ld,
        r2Threshold = 0.3,
        minimumLd = 1,
        verbose = FALSE
    )
    expect_equal(fromChr, fromLd)
    expect_error(
        raiss(knownZscores = known, genotypes = ld),
        "must be a genotype panel or a GenotypeMixture"
    )
})

test_that("raiss takes a summary-statistics GRanges as knownZscores", {
    set.seed(3)
    p <- 30L
    R <- 0.6^abs(outer(seq_len(p), seq_len(p), "-"))
    pos <- seq(100L, by = 100L, length.out = p)
    ref_panel <- data.frame(
        chrom = "1",
        pos = pos,
        variant_id = sprintf("chr1:%d:A:G", pos),
        A1 = "G",
        A2 = "A"
    )
    known <- ref_panel[seq_len(p) %% 3L != 0L, ]
    known$z <- rnorm(nrow(known))
    # One entry as sumStats() returns it: the variants as ranges, the
    # statistics as SNP / A1 / A2 / Z columns.
    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(known$pos, width = 1L),
        SNP = known$variant_id,
        A1 = known$A1,
        A2 = known$A2,
        Z = known$z
    )
    fromTable <- raiss(
        ref_panel,
        known,
        ldMatrix = R,
        r2Threshold = 0.3,
        minimumLd = 1,
        verbose = FALSE
    )
    fromRanges <- raiss(
        ref_panel,
        gr,
        ldMatrix = R,
        r2Threshold = 0.3,
        minimumLd = 1,
        verbose = FALSE
    )
    expect_equal(
        arrange(fromRanges$resultNofilter, variant_id)$z,
        arrange(fromTable$resultNofilter, variant_id)$z
    )
    expect_equal(
        arrange(fromRanges$resultFilter, variant_id)$variant_id,
        arrange(fromTable$resultFilter, variant_id)$variant_id
    )
    expect_error(
        raiss(ref_panel, gr[, "Z"], ldMatrix = R),
        "lacks column\\(s\\) SNP, A1, A2"
    )
})

test_that("Default parameters for raiss work correctly", {
    input_data <- generateDummyData()
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

    expect_error(
        raissModel(zt_invalid, sig_t, sig_i_t),
        "Assertion on 'zt' failed: Must be of type 'numeric'"
    )
    expect_error(
        raissModel(zt, sig_t_invalid, sig_i_t),
        "Assertion on 'sigT' failed: Must be of type 'numeric'"
    )
    expect_error(
        raissModel(zt, sig_t, sig_i_t_invalid),
        "Assertion on 'sigIT' failed: Must be of type 'numeric'"
    )
})

test_that("Boundary Conditions of raissModel", {
    zt <- c(1.2, 0.5)
    sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
    sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

    zt_empty <- numeric(0)
    sig_t_empty <- matrix(numeric(0), nrow = 0)
    sig_i_t_empty <- matrix(numeric(0), nrow = 0)

    expect_error(
        raissModel(zt_empty, sig_t, sig_i_t),
        "Assertion on 'zt' failed: Must have length >= 1"
    )
    expect_error(
        raissModel(zt, sig_t_empty, sig_i_t),
        "Assertion on 'sigT' failed: Must have exactly 2 rows"
    )
    expect_error(
        raissModel(zt, sig_t, sig_i_t_empty),
        "Assertion on 'sigIT' failed: Must have at least 1 rows"
    )
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

    unknowns <- sample(seq_len(nrow(ref_panel)), 5)

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

generateFroTestData <- function(seed = 1) {
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
    test_data <- generateFroTestData()
    output <- filterRaissOutput(test_data)$zscores
    expect_true(all(
        c('variant_id', 'A1', 'A2', 'z', 'Var', 'raissLdScore') %in%
            names(output)
    ))
})

test_that("raissR2 is calculated correctly in filterRaissOutput", {
    test_data <- generateFroTestData()
    output <- filterRaissOutput(test_data)$zscores
    expected_R2 <- 1 - test_data[which(test_data$raissLdScore >= 5), ]$Var
    expect_equal(output$raissR2, expected_R2[which(expected_R2 > 0.6)])
})

test_that("Function returns the correct subset in filterRaissOutput", {
    test_data <- generateFroTestData()
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

generateMockDataForComputeVar <- function(seed = 1) {
    return(
        list(
            sig_i_t_1 = matrix(c(1, 2, 3, 4), nrow = 2),
            sig_t_inv_1 = matrix(c(5, 6, 7, 8), nrow = 2),
            lamb_1 = 0.5
        )
    )
}

test_that("computeVar returns var and raissLdScore with and without batch", {
    input_data <- generateMockDataForComputeVar()
    for (batch in c(TRUE, FALSE)) {
        result <- computeVar(
            input_data$sig_i_t_1,
            input_data$sig_t_inv_1,
            input_data$lamb_1,
            batch = batch
        )
        lbl <- str_c("batch = ", batch)
        expect_true(is.list(result), info = lbl)
        expect_length(result, 2)
        expect_true(
            all(c("var", "raissLdScore") %in% names(result)),
            info = lbl
        )
        expect_true(is.numeric(result$var), info = lbl)
        expect_true(is.numeric(result$raissLdScore), info = lbl)
    }
})

test_that("checkInversion correctly identifies inverse matrices in", {
    sig_t <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
    sig_t_inv <- solve(sig_t)
    expect_true(checkInversion(sig_t, sig_t_inv))
})

test_that("checkInversion reports a failed inversion as FALSE, not text", {
    # all.equal() describes a mismatch as a string; passed through, that put
    # text in the logical correctInversion column exactly when it failed.
    sig_t <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
    expect_false(checkInversion(sig_t, diag(2)))
})

test_that("varInBoundaries sets boundaries correctly", {
    lamb_test <- 0.05
    var <- c(-1, 0, 0.5, 1.04, 1.05)

    result <- varInBoundaries(var, lamb_test)

    # Below 0 is floored to 0; above 0.99999 + lamb is capped to 1; anything
    # in between, the bounds included, is unchanged.
    expect_equal(result[1], 0)
    expect_equal(result[2], 0)
    expect_equal(result[3], 0.5)
    expect_equal(result[4], 1.04)
    expect_equal(result[5], 1)
})

test_that("invertMat computes correct pseudo-inverse", {
    # The diagonal is replaced by 1 + lamb, not added to: what is inverted
    # is the ridged matrix.
    mat <- matrix(c(1, 2, 3, 4), nrow = 2)
    result <- invertMat(mat, lamb = 0.5, rcond = 1e-7)
    ridged <- matrix(c(1.5, 2, 3, 1.5), nrow = 2)
    expect_equal(ridged %*% result, diag(2))
    expect_equal(result %*% ridged, diag(2))
    # A singular ridged matrix has no inverse; the result is its
    # Moore-Penrose pseudo-inverse.
    ones <- matrix(1, nrow = 2, ncol = 2)
    pinv <- invertMat(ones, lamb = 0, rcond = 1e-7)
    expect_equal(pinv, matrix(0.25, nrow = 2, ncol = 2))
    expect_equal(ones %*% pinv %*% ones, ones)
    expect_equal(pinv %*% ones %*% pinv, pinv)
})

test_that("invertMat handles errors and retries", {
    # ginv() can fail, as when LAPACK's SVD does not converge. invertMat()
    # then tries once more with the ridge and the cutoff each raised 10%.
    ridges <- numeric(0)
    tols <- numeric(0)
    local_mocked_bindings(
        ginv = function(X, tol = sqrt(.Machine$double.eps)) {
            ridges <<- c(ridges, X[1, 1])
            tols <<- c(tols, tol)
            if (length(tols) == 1L) {
                stop("SVD did not converge")
            }
            MASS::ginv(X, tol = tol)
        },
        .package = "pecotmr"
    )
    mat <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
    result <- invertMat(mat, lamb = 0.1, rcond = 1e-7)
    expect_equal(ridges, c(1.1, 1.11))
    expect_equal(tols, c(1e-7, 1.1e-7))
    # What comes back is the retry's inverse, of the matrix ridged to 1.11.
    retried <- matrix(c(1.11, 0.5, 0.5, 1.11), nrow = 2)
    expect_equal(retried %*% result, diag(2))
})

test_that("invertMatRecursive correctly inverts a valid square matrix", {
    mat <- matrix(c(2, -1, -1, 2), nrow = 2)
    result <- invertMatRecursive(mat, lamb = 0.5, rcond = 0.01)
    # What is inverted is mat with its diagonal set to 1 + lamb.
    ridged <- matrix(c(1.5, -1, -1, 1.5), nrow = 2)
    expect_equal(ridged %*% result, diag(2))
    expect_equal(result %*% ridged, diag(2))
})

test_that("invertMatRecursive gives a non-square matrix its pseudo-inverse", {
    mat <- matrix(1:6, nrow = 2)
    result <- invertMatRecursive(mat, lamb = 0.5, rcond = 0.01)
    # The ridge sets the leading diagonal, [1, 1] and [2, 2]. The ridged
    # 2 x 3 matrix A has full row rank, so its Moore-Penrose inverse is the
    # right inverse t(A) (A t(A))^-1.
    ridged <- matrix(c(1.5, 2, 3, 1.5, 5, 6), nrow = 2)
    expect_equal(dim(result), c(3L, 2L))
    expect_equal(ridged %*% result, diag(2))
    expect_equal(result, t(ridged) %*% solve(ridged %*% t(ridged)))
})

test_that("invertMatRecursive retries via invertMat, three attempts in all", {
    # A failed inversion falls back to invertMat() with lamb and rcond
    # raised 10%, which retries once more itself: three attempts, each 10%
    # above the last, before the error reaches the caller.
    failures <- 2L
    ridges <- numeric(0)
    tols <- numeric(0)
    local_mocked_bindings(
        ginv = function(X, tol = sqrt(.Machine$double.eps)) {
            ridges <<- c(ridges, X[1, 1])
            tols <<- c(tols, tol)
            if (length(tols) <= failures) {
                stop("SVD did not converge")
            }
            MASS::ginv(X, tol = tol)
        },
        .package = "pecotmr"
    )
    mat <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
    result <- invertMatRecursive(mat, lamb = 0.1, rcond = 0.01)
    expect_equal(ridges, 1 + 0.1 * c(1, 1.1, 1.21))
    expect_equal(tols, 0.01 * c(1, 1.1, 1.21))
    lastTry <- matrix(c(1.121, 0.5, 0.5, 1.121), nrow = 2)
    expect_equal(lastTry %*% result, diag(2))
    # A third failure is not retried: it reaches the caller.
    failures <- 3L
    tols <- numeric(0)
    expect_error(
        invertMatRecursive(mat, lamb = 0.1, rcond = 0.01),
        "SVD did not converge"
    )
    expect_length(tols, 3L)
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

test_that("raissSingleMatrix returns known z when nothing to impute", {
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

test_that("raissSingleMatrixFromX returns NULL with no known overlap", {
    set.seed(42)
    n <- 50
    p <- 10
    ref_panel <- data.frame(
        chrom = rep(1, p),
        pos = seq(10, p * 10, 10),
        variant_id = paste0("rs", seq_len(p)),
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

test_that("raissSingleMatrixFromX returns known z when nothing to impute", {
    set.seed(42)
    n <- 50
    p <- 5
    ref_panel <- data.frame(
        chrom = rep(1, p),
        pos = seq(10, p * 10, 10),
        variant_id = paste0("rs", seq_len(p)),
        A1 = rep("A", p),
        A2 = rep("G", p),
        stringsAsFactors = FALSE
    )
    known_zscores <- data.frame(
        chrom = rep(1, p),
        pos = seq(10, p * 10, 10),
        variant_id = paste0("rs", seq_len(p)),
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

# ============================================================================
# Tests for SVD-based genotype matrix path (raissSingleMatrixFromX)
# ============================================================================

#' Helper: generate a genotype matrix X with corresponding ref_panel,
#' known_zscores, and LD matrix R for equivalence testing.
generateXTestData <- function(n = 200, p = 100, n_known = 50, seed = 42) {
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

# The SVD route on a standardized genotype matrix. raiss() takes genotypes
# only as a panel, which it standardizes itself, so these numeric checks call
# the route it hands that matrix to.
.raissX <- function(
    refPanel,
    knownZscores,
    X,
    lamb = 0.01,
    rcond = 0.01,
    svdTol = 1e-8,
    r2Threshold = 0.6,
    minimumLd = 5,
    verbose = FALSE
) {
    pecotmr:::.raissGenotypePath(
        refPanel,
        knownZscores,
        X,
        list(
            lamb = lamb,
            rcond = rcond,
            svdTol = svdTol,
            r2Threshold = r2Threshold,
            minimumLd = minimumLd,
            verbose = verbose
        )
    )
}

test_that(".safeSvd basic functionality", {
    set.seed(1)
    mat <- matrix(rnorm(20), nrow = 5, ncol = 4)
    s <- pecotmr:::.safeSvd(mat)
    expect_equal(length(s$d), min(5, 4))
    expect_true(all(s$d > 0))
    # Reconstruct
    reconstructed <- s$u %*% diag(s$d) %*% t(s$v)
    expect_equal(mat, reconstructed, tolerance = 1e-10)
})

test_that(".safeSvd filters small singular values", {
    set.seed(2)
    # Create rank-2 matrix
    u <- matrix(rnorm(10), nrow = 5, ncol = 2)
    v <- matrix(rnorm(8), nrow = 4, ncol = 2)
    mat <- u %*% t(v) + matrix(rnorm(20) * 1e-12, nrow = 5, ncol = 4)
    s <- pecotmr:::.safeSvd(mat, tol = 1e-6)
    expect_equal(length(s$d), 2)
})

test_that(".safeSvd maxRank works", {
    set.seed(3)
    mat <- matrix(rnorm(50), nrow = 10, ncol = 5)
    s <- pecotmr:::.safeSvd(mat, maxRank = 2)
    expect_equal(length(s$d), 2)
    expect_equal(ncol(s$u), 2)
    expect_equal(ncol(s$v), 2)
})

test_that(".safeSvd rejects all-zero matrix", {
    mat <- matrix(0, nrow = 5, ncol = 3)
    expect_error(pecotmr:::.safeSvd(mat), "all-zero")
})

test_that("X path matches R path: basic equivalence (n > p)", {
    data <- generateXTestData(n = 200, p = 100, n_known = 50, seed = 42)

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
    result_X <- .raissX(
        data$ref_panel,
        data$known_zscores,
        data$X,
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
    data <- generateXTestData(n = 50, p = 200, n_known = 100, seed = 123)

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
    result_X <- .raissX(
        data$ref_panel,
        data$known_zscores,
        data$X,
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


test_that("the genotype and LD routes drop the same weak directions", {
    # Blocks of near-identical columns: a few large eigenvalues and many tiny
    # ones, so rcond's cutoff removes directions the ridge alone would keep.
    set.seed(17)
    n <- 150
    p <- 40
    base <- matrix(sample(0:2, n * 5, replace = TRUE), nrow = n)
    X_raw <- base[, rep(1:5, each = 8)]
    flips <- cbind(sample(n, 6 * p, replace = TRUE), rep(seq_len(p), each = 6))
    X_raw[flips] <- sample(0:2, nrow(flips), replace = TRUE)
    ref_panel <- data.frame(
        chrom = 1,
        pos = seq_len(p),
        variant_id = paste0("rs", seq_len(p)),
        A1 = "A",
        A2 = "G"
    )
    X <- `colnames<-`(scale(X_raw), ref_panel$variant_id)
    R <- cor(X_raw)
    dimnames(R) <- list(ref_panel$variant_id, ref_panel$variant_id)
    known_idx <- which(seq_len(p) %% 4L != 0L)
    known <- ref_panel[known_idx, ]
    known$z <- rnorm(nrow(known))
    ridged <- `diag<-`(R[known_idx, known_idx], 1.01)
    ev <- eigen(ridged, symmetric = TRUE, only.values = TRUE)$values
    expect_true(any(ev < 0.01 * ev[[1L]]))

    viaLd <- raiss(
        ref_panel,
        known,
        ldMatrix = R,
        r2Threshold = 0,
        minimumLd = 0,
        verbose = FALSE
    )$resultNofilter
    viaX <- .raissX(
        ref_panel,
        known,
        X,
        svdTol = 1e-12,
        r2Threshold = 0,
        minimumLd = 0
    )$resultNofilter
    expect_equal(viaX$z, viaLd$z, tolerance = 1e-8)
    expect_equal(viaX$Var, viaLd$Var, tolerance = 1e-8)
    expect_equal(viaX$raissLdScore, viaLd$raissLdScore, tolerance = 1e-8)
    # Without the cutoff the genotype route keeps those directions and its
    # imputation differs.
    uncut <- .raissX(
        ref_panel,
        known,
        X,
        rcond = 0,
        svdTol = 1e-12,
        r2Threshold = 0,
        minimumLd = 0
    )$resultNofilter
    expect_false(isTRUE(all.equal(uncut$z, viaLd$z, tolerance = 1e-6)))
})

test_that("X path matches R path: varying lambda", {
    data <- generateXTestData(n = 150, p = 80, n_known = 40, seed = 7)

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
        result_X <- .raissX(
            data$ref_panel,
            data$known_zscores,
            data$X,
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

test_that("X path handles single unknown variant", {
    data <- generateXTestData(n = 100, p = 50, n_known = 49, seed = 15)

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
    result_X <- .raissX(
        data$ref_panel,
        data$known_zscores,
        data$X,
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
    data <- generateXTestData(n = 100, p = 50, n_known = 1, seed = 20)

    result_X <- .raissX(
        data$ref_panel,
        data$known_zscores,
        data$X,
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
    data <- generateXTestData(n = 200, p = 100, n_known = 50, seed = 42)

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
    result_X <- .raissX(
        data$ref_panel,
        data$known_zscores,
        data$X,
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

test_that("the SVD route needs standardized genotypes", {
    # Why raiss() standardizes a panel's dosages itself: raw dosages through
    # the same route give different answers from the LD path.
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
    result_raw_X <- .raissX(
        ref_panel,
        known_zscores,
        X_raw,
        lamb = 0.01,
        svdTol = 1e-8,
        r2Threshold = 0.6,
        minimumLd = 0,
        verbose = FALSE
    )
    result_scaled_X <- .raissX(
        ref_panel,
        known_zscores,
        X_scaled,
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

test_that("raiss rejects both ldMatrix and genotypes", {
    data <- generateXTestData(n = 50, p = 20, n_known = 10, seed = 1)
    expect_error(
        raiss(
            data$ref_panel,
            data$known_zscores,
            ldMatrix = data$R,
            genotypes = data$X
        ),
        "not both"
    )
})

test_that("raiss refuses genotypes that are not a panel", {
    # A bare matrix leaves centring, scaling, column order and allele coding
    # to the caller, and gets every one of them wrong silently.
    data <- generateXTestData(n = 50, p = 20, n_known = 10, seed = 1)
    expect_error(
        raiss(data$ref_panel, data$known_zscores, genotypes = data$X),
        "`genotypes` must be a genotype panel or a GenotypeMixture"
    )
})

test_that("raiss rejects neither ldMatrix nor genotypes", {
    data <- generateXTestData(n = 50, p = 20, n_known = 10, seed = 1)
    expect_error(
        raiss(data$ref_panel, data$known_zscores),
        "Provide either"
    )
})

test_that("raiss imputes from a genotype panel through the SVD route", {
    # The panel's dosages are standardized here and the p x p matrix is never
    # formed: the fit is the SVD route's on the same standardized dosages.
    set.seed(11)
    panel <- readGenotypes(test_path("test_data", "test_variants.bed"))[1:30, ]
    ref_panel <- pecotmr:::.refVariantsFromLdRef(panel)
    known <- ref_panel[seq_len(30L) %% 3L != 0L, ]
    known$z <- rnorm(nrow(known))
    expect_message(
        fromPanel <- raiss(knownZscores = known, genotypes = panel),
        "SVD-based"
    )
    direct <- .raissX(
        ref_panel,
        known,
        pecotmr:::.ldRefFactor(panel, seq_len(30L)),
        verbose = FALSE
    )
    expect_equal(fromPanel$resultNofilter, direct$resultNofilter)
    expect_error(
        raiss(knownZscores = known, ldMatrix = panel),
        "Pass genotypes as `genotypes`"
    )
})

test_that("the SVD route reports its progress when verbose", {
    data <- generateXTestData(n = 60, p = 20, n_known = 10, seed = 3)
    expect_message(
        res <- .raissX(
            data$ref_panel,
            data$known_zscores,
            data$X,
            r2Threshold = 0,
            minimumLd = 0,
            verbose = TRUE
        ),
        "Processing genotype matrix"
    )
    expect_equal(nrow(res$resultNofilter), nrow(data$ref_panel))
})

context("slalom")

# ============================================================================
# Helper: build a valid positive-definite LD matrix from a genotype matrix
# ============================================================================
makeSyntheticLd <- function(n_samples, n_snps, seed = 1) {
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
    for (i in seq_len(n - 1)) {
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
    expect_error(slalom(zScore = z, R = R), "R.*Must have exactly 10 rows")
})

test_that("slalom honours ldMethod when given genotypes", {
    # `ldMethod` is documented as "Passed to computeLd", but slalom used to
    # correlate the lead column of X with cor() directly, so "population"
    # and "gcta" were accepted and silently given sample LD. It now resolves
    # LD the same way dentist does.
    set.seed(11)
    X <- matrix(sample(0:2, 100 * 6, replace = TRUE), nrow = 100, ncol = 6)
    colnames(X) <- paste0("snp", seq_len(6))
    z <- rnorm(6)
    seen <- NULL
    local_mocked_bindings(
        computeLd = function(X, method = "sample", ...) {
            seen <<- method
            stats::cor(X)
        },
        .package = "pecotmr"
    )
    invisible(slalom(zScore = z, X = X, ldMethod = "gcta"))
    expect_equal(seen, "gcta")
    # The sample estimator reads the lead column off X's LD factor instead,
    # never forming the matrix.
    seen <- NULL
    invisible(slalom(zScore = z, X = X))
    expect_null(seen)
})

test_that("slalom from genotypes matches slalom from their LD", {
    set.seed(12)
    X <- matrix(sample(0:2, 80 * 7, replace = TRUE), nrow = 80, ncol = 7)
    X[, 4] <- 2
    X[c(3, 17), 2] <- NA
    colnames(X) <- paste0("snp", seq_len(7))
    z <- rnorm(7)
    expect_equal(
        slalom(zScore = z, X = X),
        slalom(zScore = z, R = computeLd(X))
    )
})

test_that("slalom still refuses both or neither of R and X", {
    # resolveLdInput now carries the exactly-one-of rule that .slalomValidate
    # used to duplicate; the errors must not have gone with it.
    z <- rnorm(4)
    R <- diag(4)
    X <- matrix(rnorm(40), nrow = 10, ncol = 4)
    expect_error(slalom(zScore = z), "Either R .* or X .* must be provided")
    expect_error(
        slalom(zScore = z, R = R, X = X),
        "Provide either R or X, not both"
    )
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
    colnames(X) <- paste0("snp", seq_len(n_snps))
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
    # ... and, as the top-PIP variant, it heads the CS.
    expect_equal(result$summary$cs95[1], causal)
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

test_that("pvalue lead is the most negative z, abf lead the highest PIP", {
    # Asymmetric magnitudes, so the two choices disagree.
    z <- c(-3.0, 5.0, 0.1, -0.2, 0.3)
    R <- diag(5)

    result_pv <- slalom(zScore = z, R = R, leadVariantChoice = "pvalue")
    result_abf <- slalom(zScore = z, R = R, leadVariantChoice = "abf")

    expect_equal(result_pv$summary$leadPipVariant, 1)
    expect_equal(result_abf$summary$leadPipVariant, 2)
    expect_equal(
        result_abf$summary$leadPipVariant,
        which.max(result_abf$data$prob)
    )
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
    syn <- makeSyntheticLd(200, n, seed = 401)
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
    syn <- makeSyntheticLd(300, n, seed = 402)
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

test_that("edge case: two variants only", {
    z <- c(-3, 2)
    R <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)

    result <- slalom(zScore = z, R = R)

    expect_equal(nrow(result$data), 2)
    expect_equal(sum(result$data$prob), 1, tolerance = 1e-12)
    expect_equal(result$summary$leadPipVariant, 1)
    expect_gt(result$data$prob[1], result$data$prob[2])
})

# ============================================================================
# Parameter variation
# ============================================================================

test_that("larger abf_prior_variance concentrates PIPs on strong signals", {
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
    syn <- makeSyntheticLd(200, n, seed = 701)
    z <- rnorm(n, sd = 2)
    z[1] <- -5

    result_low <- slalom(zScore = z, R = syn$R, r2Threshold = 0.1)
    result_high <- slalom(zScore = z, R = syn$R, r2Threshold = 0.9)

    expect_gte(result_low$summary$nR2, result_high$summary$nR2)
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

test_that("fraction = 0 with no outliers (identity LD, consistent z)", {
    n <- 5
    z <- c(-3, 0, 0, 0, 0)
    R <- diag(n)

    result <- slalom(zScore = z, R = R)
    expect_equal(result$summary$fraction, 0)
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
    syn <- makeSyntheticLd(500, 20, seed = 1000)
    z <- rep(0, 20)
    z[3] <- 5

    result <- slalom(zScore = z, R = syn$R)

    expect_true(3 %in% result$summary$cs95)
    expect_equal(sum(result$data$prob), 1, tolerance = 1e-12)
})

test_that("realistic LD: DENTIST-S detects outlier in correlated block", {
    # LD decaying with distance from the lead, r = 0.9^d, so its neighbour
    # is in LD with it by construction (r2 = 0.81, above the 0.2 cutoff).
    p <- 15L
    R <- 0.9^abs(outer(seq_len(p), seq_len(p), "-"))
    # Every z-score as the lead's LD predicts it, then the neighbour's
    # corrupted by +10.
    z <- R[, 1] * (-5)
    z[2] <- z[2] + 10

    result <- slalom(
        zScore = z,
        R = R,
        r2Threshold = 0.2,
        nlog10pDentistSThreshold = 3.0
    )

    expect_equal(result$summary$leadPipVariant, 1L)
    # The corrupted variant, and only it, is flagged.
    expect_equal(which(result$data$outliers), 2L)
})

test_that("realistic LD: no outliers when z perfectly matches LD structure", {
    # LD decaying with distance from the lead, r = 0.9^d: the lead and
    # variants 2-6 clear the 0.3 r2 cutoff, so only their DENTIST-S
    # statistic keeps them from being flagged.
    p <- 10L
    R <- 0.9^abs(outer(seq_len(p), seq_len(p), "-"))
    z <- R[, 1] * (-4)

    result <- slalom(
        zScore = z,
        R = R,
        r2Threshold = 0.3,
        nlog10pDentistSThreshold = 4.0
    )

    expect_equal(result$summary$leadPipVariant, 1L)
    expect_equal(result$summary$nR2, 6L)
    # Each z is exactly what the lead's LD predicts, so every statistic is 0
    # and nothing is flagged.
    expect_equal(result$data$nlog10p_dentist_s[-1], rep(0, p - 1L))
    expect_false(any(result$data$outliers[-1]))
    expect_equal(result$summary$nDentistSOutlier, 0L)
})

# ============================================================================
# Numerical stability
# ============================================================================

test_that("log-sum-exp trick prevents overflow with extreme z-scores", {
    z <- c(100, 0, -100)
    R <- diag(3)

    result <- slalom(zScore = z, R = R)

    expect_false(any(is.na(result$data$prob)))
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
# .matchAgainstLdRef, .applyLdMismatchQcToEntry) run
# for real on the synthetic fixture.

# ===========================================================================
# Fixture builders
# ===========================================================================

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
    studyName = "g1"
) {
    GwasSumStats(
        studyName = studyName,
        entry = list(.ssQ_makeEntryGr(snp_ids, positions)),
        genome = "hg19",
        ldReference = .ssQ_makeHandle()
    )
}

.ssQ_mockExtractor <- function(seed = 13, n_samples = 60L) {
    function(handle, snpIdx, meanImpute = TRUE) {
        set.seed(seed)
        panel <- matrix(
            rbinom(n_samples * nrow(snpInfo(handle)), 2, 0.3),
            nrow = n_samples,
            ncol = nrow(snpInfo(handle)),
            dimnames = list(sampleIds(handle), snpInfo(handle)$SNP)
        )
        sub <- panel[, snpIdx, drop = FALSE]
        rr <- GenomicRanges::GRanges(
            seqnames = paste0("chr", snpInfo(handle)$CHR[snpIdx]),
            ranges = IRanges::IRanges(
                start = snpInfo(handle)$BP[snpIdx],
                width = 1L
            )
        )
        S4Vectors::mcols(rr) <- S4Vectors::DataFrame(
            SNP = snpInfo(handle)$SNP[snpIdx],
            A1 = snpInfo(handle)$A1[snpIdx],
            A2 = snpInfo(handle)$A2[snpIdx]
        )
        cd <- S4Vectors::DataFrame(
            sampleId = sampleIds(handle),
            row.names = sampleIds(handle)
        )
        dosage <- t(sub)
        rownames(dosage) <- snpInfo(handle)$SNP[snpIdx]
        colnames(dosage) <- sampleIds(handle)
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
        "Must inherit from class 'QtlSumStats'/'GwasSumStats'"
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
        res <- summaryStatsQc(
            ss,
            panelFilterParam = PanelFilterParam(mafCutoff = 0.05)
        ),
        "skipping the MAF filter"
    )
    expect_s4_class(res, "GwasSumStats")
})

test_that("summaryStatsQc: infoCutoff > 0 with no INFO column errors", {
    ss <- .ssQ_makeGwasSumStats()
    expect_error(
        summaryStatsQc(
            ss,
            sumstatsFilterParam = SumstatsFilterParam(infoCutoff = 0.5)
        ),
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
        studyName = "g1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ssQ_makeHandle()
    )
    out <- summaryStatsQc(
        ss,
        signalScreenParam = SignalScreenParam(pip = 0.5),
        sumstatsFilterParam = SumstatsFilterParam(nCutoff = 0)
    )
    snps <- as.character(S4Vectors::mcols(out[[1L]])$SNP)
    expect_false("rsX" %in% snps) # dropped by harmonization
    # screen (post-harmonization) skips the region
    expect_equal(length(out[[1L]]), 0L)
})

test_that("summaryStatsQc: PIP screen off leaves the harmonized set intact", {
    # Same harmonized variants, screen disabled (no metric set): the three
    # panel-matched variants survive (behavior-invariance for the screen-off path).
    gr <- .ssQ_makeEntryGr(
        c("rs1", "rs2", "rs3", "rsX"),
        positions = c(100L, 200L, 300L, 9999L)
    )
    mc <- S4Vectors::mcols(gr)
    mc$Z <- c(0.2, 0.3, 0.1, 12)
    S4Vectors::mcols(gr) <- mc
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ssQ_makeHandle()
    )
    out <- summaryStatsQc(
        ss,
        signalScreenParam = SignalScreenParam(pip = 0),
        sumstatsFilterParam = SumstatsFilterParam(nCutoff = 0)
    )
    snps <- as.character(S4Vectors::mcols(out[[1L]])$SNP)
    # SNP is re-keyed to the panel-harmonized id (chr:pos:A2:A1) after
    # harmonization; the panel is A1=A / A2=G, so the surviving three become
    # chr1:<pos>:G:A.
    expect_setequal(snps, c("chr1:100:G:A", "chr1:200:G:A", "chr1:300:G:A"))
})

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

test_that("summaryStatsQc: harmonization re-keys SNP and sign-flips Z", {
    # pos200 is allele-swapped vs the panel, so harmonization sign-flips its Z and
    # the SNP must be re-keyed to the panel-orientation id (chr1:200:G:A), not
    # left at the input-orientation chr1:200:A:G. Exact matches are unchanged.
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(.ssQ_makeFlipEntryGr()),
        genome = "hg19",
        ldReference = .ssQ_makeHandle()
    )
    out <- summaryStatsQc(
        ss,
        signalScreenParam = SignalScreenParam(pip = 0),
        sumstatsFilterParam = SumstatsFilterParam(nCutoff = 0)
    )
    e <- out[[1L]]
    o <- order(GenomicRanges::start(e))
    snp <- as.character(S4Vectors::mcols(e)$SNP)[o]
    z <- S4Vectors::mcols(e)$Z[o]
    expect_equal(snp, c("chr1:100:G:A", "chr1:200:G:A", "chr1:300:G:A"))
    expect_equal(z, c(1.0, -2.0, 1.5)) # only the swapped variant's Z flips
})

# ---------------------------------------------------------------------------
# Tag-allele panel ids (e.g. ADSP R5 INS/DEL): repair from the allele columns,
# never drop the variant
# ---------------------------------------------------------------------------

test_that(".raissUnsafeToImpute does not refuse a tag-allele id", {
    # A tag id is repaired at the panel boundary, so RAISS sees a normal indel
    # id; only the flip-pair / flip-of-known ambiguities are unsafe.
    expect_false(pecotmr:::.raissUnsafeToImpute(
        "chr21:13988152:A:AT",
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

test_that("harmonization keeps a variant on a tag-named panel position", {
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
    out0 <- pecotmr:::.matchAgainstLdRef(df, h, matchMinProp = 0)
    expect_setequal(pecotmr:::parseVariantId(out0$SNP)$pos, c(100L, 200L))
    # tag the pos-200 panel id (A1/A2 columns stay the real alleles)
    h@snpInfo$SNP[h@snpInfo$BP == 200L] <- "chr1:200:INS:A"
    outT <- pecotmr:::.matchAgainstLdRef(df, h, matchMinProp = 0)
    # the tag costs nothing: both variants survive, with unchanged Z
    expect_setequal(pecotmr:::parseVariantId(outT$SNP)$pos, c(100L, 200L))
    expect_equal(outT$SNP, out0$SNP)
    expect_equal(outT$Z, out0$Z)
})

test_that("a tag-named panel entry survives the end-to-end sketch subset", {
    h <- .ssQ_makeHandleVid()
    h@snpInfo$SNP[h@snpInfo$BP == 200L] <- "chr1:200:INS:A"
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(.ssQ_makeEntryGr()),
        genome = "hg19",
        ldReference = h
    )
    out <- summaryStatsQc(
        ss,
        signalScreenParam = SignalScreenParam(pip = 0),
        sumstatsFilterParam = SumstatsFilterParam(nCutoff = 0)
    )
    snp <- as.character(S4Vectors::mcols(out[[1L]])$SNP)
    expect_true("chr1:200:G:A" %in% snp)
    # the retained sketch keeps the panel row the QC'd entry still refers to
    expect_true(
        "chr1:200:G:A" %in% pecotmr:::.ldRefMatchIds(ldReference(out)[[1L]])
    )
})

test_that("the sketch survives a chr:pos:A1:A2 panel (PLINK .bim order)", {
    # Pre-fix regression: harmonization canonicalizes the entry ids to
    # chr:pos:A2:A1, while panel ids are passed through verbatim whenever their
    # allele fields are already valid DNA. A panel keyed the other way round --
    # which is how a PLINK .bim writes ids, whereas a .pvar writes REF:ALT and
    # is already canonical -- matched nothing under string equality, so
    # .subsetLdRefToIds() silently emptied the sketch and RSS fine-mapping /
    # TWAS weights aborted later with "not present in the LD reference".
    h <- .ssQ_makeHandleVid()
    si <- snpInfo(h)
    h@snpInfo$SNP <- paste0("chr1:", si$BP, ":", si$A1, ":", si$A2)
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(.ssQ_makeEntryGr(
            snp_ids = paste0("chr1:", c(100L, 200L, 300L, 400L), ":G:A")
        )),
        genome = "hg19",
        ldReference = h
    )
    out <- summaryStatsQc(
        ss,
        signalScreenParam = SignalScreenParam(pip = 0),
        sumstatsFilterParam = SumstatsFilterParam(nCutoff = 0)
    )
    kept <- pecotmr:::.ldRefMatchIds(ldReference(out)[[1L]])
    expect_length(kept, 4L)
    expect_setequal(
        pecotmr:::parseVariantId(kept)$pos,
        c(100L, 200L, 300L, 400L)
    )
})

test_that("summaryStatsQc: slalom z-mismatch resolves sign-flipped variants", {
    # Pre-fix regression: a sign-flipped variant kept its input-orientation SNP,
    # which is absent from the panel, so .applyLdMismatchQcToEntry errored with
    # "absent from the ldReference panel". After the harmonization re-key the
    # SNP matches the panel and slalom QC runs to completion.
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(.ssQ_makeFlipEntryGr()),
        genome = "hg19",
        ldReference = .ssQ_makeHandleVid()
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
            ldMismatchQcMethod = "slalom",
            signalScreenParam = SignalScreenParam(pip = 0),
            sumstatsFilterParam = SumstatsFilterParam(nCutoff = 0)
        )
    )
    snp <- as.character(S4Vectors::mcols(out[[1L]])$SNP)
    expect_true("chr1:200:G:A" %in% snp) # the sign-flipped variant survived
    expect_equal(length(out[[1L]]), 3L)
})

test_that("summaryStatsQc: ldMismatchQc reconciles a chr-prefix difference", {
    # Panel SNP ids are non-chr-prefixed positional; QC re-keys the entry to the
    # canonical chr-prefixed form, so the opt-in z-mismatch panel match must
    # reconcile the prefix (previously errored "absent from the ldReference
    # panel").
    h <- .ssQ_makeHandle()
    h@snpInfo$SNP <- paste0("1:", snpInfo(h)$BP, ":G:A") # non-chr-prefixed
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(.ssQ_makeEntryGr()),
        genome = "hg19",
        ldReference = h
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
            ldMismatchQcMethod = "slalom",
            signalScreenParam = SignalScreenParam(pip = 0),
            sumstatsFilterParam = SumstatsFilterParam(nCutoff = 0)
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

test_that("summaryStatsQc: vanilla run records its options and audits", {
    ss <- .ssQ_makeGwasSumStats()
    res <- summaryStatsQc(ss)
    expect_s4_class(res, "GwasSumStats")
    expect_true("options" %in% names(qcOptions(res)))
    expect_false("entryAudit" %in% names(qcOptions(res)))
    expect_equal(length(res$qcAudit), nrow(ss))
    ea <- res$qcAudit[[1L]]
    expect_equal(ea$variantsIn, 4L)
    expect_equal(ea$variantsOut, 4L)
    expect_equal(ea$matchedAgainstLdReference, 4L)
    # The per-row table: identity columns, then the audit's scalar counts.
    tbl <- qcInfo(res)
    expect_s4_class(tbl, "DataFrame")
    expect_equal(nrow(tbl), nrow(ss))
    expect_equal(tbl$study, as.character(res$study))
    expect_equal(tbl$variantsIn, 4L)
    expect_equal(tbl$variantsOut, 4L)
})

test_that("summaryStatsQc: keepVariants subsets entries, records the drop", {
    ss <- .ssQ_makeGwasSumStats()
    res <- summaryStatsQc(ss, keepVariants = c("rs1", "rs3"))
    ea <- res$qcAudit[[1L]]
    expect_equal(ea$keepVariantsDropped, 2L)
    expect_equal(ea$variantsOut, 2L)
})

test_that("summaryStatsQc: skipRegion drops overlapping variants", {
    ss <- .ssQ_makeGwasSumStats()
    res <- summaryStatsQc(
        ss,
        skipRegion = GenomicRanges::GRanges("chr1:50-150")
    )
    ea <- res$qcAudit[[1L]]
    expect_equal(ea$skipRegionDropped, 1L) # rs1 at pos 100 is dropped
})

test_that("summaryStatsQc: PIP screen triggers when no variant has signal", {
    # Build an entry with weak signal so the SER PIP screen tags everything
    # below the threshold.
    gr <- .ssQ_makeEntryGr()
    S4Vectors::mcols(gr)$Z <- rep(0.1, length(gr))
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ssQ_makeHandle()
    )
    res <- summaryStatsQc(ss, signalScreenParam = SignalScreenParam(pip = 0.99))
    ea <- res$qcAudit[[1L]]
    expect_true(isTRUE(ea$pipScreenSkipped))
    expect_match(ea$pipScreenReason, "no signals above PIP threshold")
    expect_equal(length(res[[1L]]), 0L)
})

test_that("summaryStatsQc: options block records the curated knobs", {
    ss <- .ssQ_makeGwasSumStats()
    res <- summaryStatsQc(
        ss,
        sumstatsFilterParam = SumstatsFilterParam(
            removeIndels = TRUE,
            removeStrandAmbiguous = FALSE,
            nCutoff = 10
        )
    )
    opts <- qcOptions(res)$options
    expect_true(opts$removeIndels)
    expect_false(opts$removeStrandAmbiguous)
    expect_equal(opts$nCutoff, 10)
})

test_that("summaryStatsQc: round-trips QtlSumStats inputs", {
    gr <- .ssQ_makeEntryGr()
    ss <- QtlSumStats(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ssQ_makeHandle()
    )
    res <- summaryStatsQc(ss)
    expect_s4_class(res, "QtlSumStats")
    expect_equal(length(res$qcAudit), 1L)
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

test_that("effectiveN: balanced == total, imbalanced smaller, guards NA", {
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

test_that("summaryStatsQc(effectiveN=TRUE): counts, no N -> N == N_eff", {
    gr <- .ssQ_makeCCEntry(
        nCase = c(100, 200, 150, 250),
        nControl = c(900, 800, 850, 750)
    )
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ssQ_makeHandle()
    )
    res <- summaryStatsQc(ss)
    # 4*case*control/(case+control) per variant.
    expect_equal(.ssQ_entryNByPos(res[[1L]]), c(360, 640, 510, 750))
    expect_identical(res$qcAudit[[1L]]$nSource, "effective")
})

test_that("summaryStatsQc(effectiveN=TRUE): counts + N -> counts win", {
    gr <- .ssQ_makeCCEntry(
        nCase = c(100, 200, 150, 250),
        nControl = c(900, 800, 850, 750),
        includeN = TRUE
    )
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ssQ_makeHandle()
    )
    expect_message(res <- summaryStatsQc(ss), "overridden by effective N")
    # N (was 1000) is replaced by the per-variant N_eff.
    expect_equal(.ssQ_entryNByPos(res[[1L]]), c(360, 640, 510, 750))
    expect_identical(res$qcAudit[[1L]]$nSource, "effective")
})

test_that("summaryStatsQc(effectiveN=TRUE): N only, no counts -> used as-is", {
    ss <- .ssQ_makeGwasSumStats() # entry carries N = 1000, no counts
    res <- summaryStatsQc(ss)
    expect_true(all(.ssQ_entryNByPos(res[[1L]]) == 1000))
    expect_identical(res$qcAudit[[1L]]$nSource, "column")
})

test_that("summaryStatsQc(effectiveN=FALSE): counts + N -> raw N", {
    gr <- .ssQ_makeCCEntry(
        nCase = c(100, 200, 150, 250),
        nControl = c(900, 800, 850, 750),
        includeN = TRUE
    )
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ssQ_makeHandle()
    )
    res <- summaryStatsQc(ss, effectiveN = FALSE)
    expect_true(all(.ssQ_entryNByPos(res[[1L]]) == 1000))
    expect_identical(res$qcAudit[[1L]]$nSource, "column")
})

test_that("summaryStatsQc(effectiveN=FALSE): counts -> raw total", {
    gr <- .ssQ_makeCCEntry(
        nCase = c(100, 200, 150, 250),
        nControl = c(900, 800, 850, 750)
    )
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ssQ_makeHandle()
    )
    res <- summaryStatsQc(ss, effectiveN = FALSE)
    # Raw total n_case + n_control per variant.
    expect_equal(.ssQ_entryNByPos(res[[1L]]), c(1000, 1000, 1000, 1000))
    expect_identical(res$qcAudit[[1L]]$nSource, "total")
})

test_that("summaryStatsQc(effectiveN=TRUE): study scalars fill all variants", {
    # Entry has an N column but NO per-variant N_CASE/N_CONTROL; the scalars win.
    gr <- .ssQ_makeEntryGr()
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ssQ_makeHandle(),
        nCase = 100,
        nControl = 900
    )
    expect_message(res <- summaryStatsQc(ss), "from study")
    # 4*100*900/1000 = 360 for every variant.
    expect_true(all(.ssQ_entryNByPos(res[[1L]]) == 360))
    expect_identical(res$qcAudit[[1L]]$nSource, "effective")
})

test_that("summaryStatsQc: study nSample is the level-4 fallback", {
    # Entry with NO per-variant N and NO case/control; only a study nSample scalar.
    gr <- .ssQ_makeEntryGr()
    mc <- S4Vectors::mcols(gr)
    mc$N <- NULL # remove per-variant N
    S4Vectors::mcols(gr) <- mc
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ssQ_makeHandle(),
        nSample = 4321
    )
    res <- summaryStatsQc(ss)
    expect_true(all(.ssQ_entryNByPos(res[[1L]]) == 4321))
    expect_identical(res$qcAudit[[1L]]$nSource, "study-n")
})

test_that("summaryStatsQc: QtlSumStats tuple nSample falls back too", {
    # Parity with the GWAS study-n fallback: a QtlSumStats carrying only a
    # tuple-level nSample (no per-variant N) fills N from the scalar and is
    # preserved on the QC'd object.
    gr <- .ssQ_makeEntryGr()
    mc <- S4Vectors::mcols(gr)
    mc$N <- NULL # remove per-variant N
    S4Vectors::mcols(gr) <- mc
    ss <- QtlSumStats(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ssQ_makeHandle(),
        nSample = 838
    )
    res <- summaryStatsQc(ss)
    expect_s4_class(res, "QtlSumStats")
    expect_true(all(.ssQ_entryNByPos(res[[1L]]) == 838))
    expect_identical(res$qcAudit[[1L]]$nSource, "study-n")
    expect_equal(as.numeric(res$nSample), 838) # slot preserved through QC
})

test_that("summaryStatsQc: counts beat nSample, and an N column beats it", {
    # per-variant counts present alongside nSample -> counts win (effective).
    gr1 <- .ssQ_makeCCEntry(
        nCase = c(100, 200, 150, 250),
        nControl = c(900, 800, 850, 750)
    )
    ss1 <- GwasSumStats(
        studyName = "g1",
        entry = list(gr1),
        genome = "hg19",
        ldReference = .ssQ_makeHandle(),
        nSample = 4321
    )
    res1 <- summaryStatsQc(ss1)
    expect_equal(.ssQ_entryNByPos(res1[[1L]]), c(360, 640, 510, 750))
    expect_identical(res1$qcAudit[[1L]]$nSource, "effective")
    # per-variant N column present alongside nSample (no counts) -> N wins (column).
    ss2 <- GwasSumStats(
        studyName = "g1",
        entry = list(.ssQ_makeEntryGr()),
        genome = "hg19",
        ldReference = .ssQ_makeHandle(),
        nSample = 4321
    )
    res2 <- summaryStatsQc(ss2) # .ssQ_makeEntryGr carries N = 1000
    expect_true(all(.ssQ_entryNByPos(res2[[1L]]) == 1000))
    expect_identical(res2$qcAudit[[1L]]$nSource, "column")
})

test_that("summaryStatsQc: effectiveN recorded in qcOptions", {
    ss <- .ssQ_makeGwasSumStats()
    expect_true(qcOptions(summaryStatsQc(ss))$options$effectiveN)
    expect_false(
        qcOptions(summaryStatsQc(ss, effectiveN = FALSE))$options$effectiveN
    )
})

test_that("summaryStatsQc: quantitative QtlSumStats: no effective N", {
    # No counts anywhere; entry carries an N column -> used as-is (nSource
    # "column"), N untouched.
    gr <- .ssQ_makeEntryGr()
    ss <- QtlSumStats(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ssQ_makeHandle()
    )
    res <- summaryStatsQc(ss)
    expect_true(all(.ssQ_entryNByPos(res[[1L]]) == 1000))
    expect_identical(res$qcAudit[[1L]]$nSource, "column")
})

# ===========================================================================
# summaryStatsQc with LD-mismatch QC enabled (mocked extractor)
# ===========================================================================

test_that("summaryStatsQc: ldMismatchQcMethod 'dentist' walks the LD branch", {
    # Panel ids follow the chr:pos:A2:A1 convention (as real LD sketches do) so the
    # post-harmonization re-keyed SNP resolves against the panel for z-mismatch QC.
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(.ssQ_makeEntryGr(
            snp_ids = paste0("rs", 1:8),
            positions = seq(100L, by = 100L, length.out = 8L)
        )),
        genome = "hg19",
        ldReference = .ssQ_makeHandleVid(snp_n = 8L)
    )
    local_mocked_bindings(
        extractBlockGenotypes = .ssQ_mockExtractor(),
        .package = "pecotmr"
    )
    res <- suppressWarnings(summaryStatsQc(ss, ldMismatchQcMethod = "dentist"))
    ea <- res$qcAudit[[1L]]
    expect_equal(ea$ldMismatchMethod, "dentist")
    expect_true("ldMismatchOutliersDropped" %in% names(ea))
})

# ===========================================================================
# summaryStatsQc with impute = TRUE: exercise the RAISS branch
# ===========================================================================

test_that("summaryStatsQc: impute = TRUE invokes RAISS, records counts", {
    # Build a sketch panel with 8 variants and a GWAS entry covering only the
    # first 4 — RAISS is asked to impute the missing 4.
    full_snp_ids <- paste0("rs", 1:8)
    full_positions <- seq(100L, by = 100L, length.out = 8L)
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(.ssQ_makeEntryGr(
            snp_ids = full_snp_ids[1:4],
            positions = full_positions[1:4]
        )),
        genome = "hg19",
        ldReference = .ssQ_makeHandle(snp_n = 8L, n_samples = 60L)
    )

    local_mocked_bindings(
        extractBlockGenotypes = .ssQ_mockExtractor(),
        .raissGenotypePath = function(
            refPanel,
            knownZscores,
            genotypeMatrix,
            p
        ) {
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
    res <- summaryStatsQc(
        ss,
        impute = TRUE,
        raissParam = RaissParam(flank = 500)
    )
    ea <- res$qcAudit[[1L]]
    expect_equal(ea$raissTotalVariants, 6L)
    expect_equal(ea$raissImputedVariants, 2L)
})

test_that("summaryStatsQc: impute scopes panel/dosage to the region", {
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
        studyName = "g1",
        entry = list(.ssQ_makeEntryGr(
            paste0("rs", 1:4),
            c(100L, 200L, 300L, 400L)
        )),
        genome = "hg19",
        ldReference = .ssQ_makeHandle(snp_n = 8L, n_samples = 60L)
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
    suppressWarnings(summaryStatsQc(
        ss,
        impute = TRUE,
        sumstatsFilterParam = SumstatsFilterParam(nCutoff = 0)
    ))
    expect_gt(length(cap$idx), 0L)
    # No read reaches past the region window's 4 panel variants.
    expect_true(all(map_lgl(cap$idx, function(i) all(i %in% 1:4))))
    expect_equal(max(map_int(cap$idx, length)), 4L)
})

test_that("summaryStatsQc: impute with a NULL raiss records 0 imputed", {
    full_snp_ids <- paste0("rs", 1:8)
    full_positions <- seq(100L, by = 100L, length.out = 8L)
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(.ssQ_makeEntryGr(
            snp_ids = full_snp_ids[1:4],
            positions = full_positions[1:4]
        )),
        genome = "hg19",
        ldReference = .ssQ_makeHandle(snp_n = 8L, n_samples = 60L)
    )

    local_mocked_bindings(
        extractBlockGenotypes = .ssQ_mockExtractor(),
        raiss = function(...) NULL,
        .package = "pecotmr"
    )
    res <- summaryStatsQc(ss, impute = TRUE)
    ea <- res$qcAudit[[1L]]
    expect_equal(ea$raissImputedVariants, 0L)
})

# ===========================================================================
# summaryStatsQc: per-step QC counter logging (concept salvaged from PR #520)
# ===========================================================================

test_that("summaryStatsQc: QC track emits per-step messages and a rollup", {
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

test_that("summaryStatsQc: skipped steps are omitted from the rollup", {
    ss <- .ssQ_makeGwasSumStats()
    local_mocked_bindings(
        extractBlockGenotypes = .ssQ_mockExtractor(),
        .package = "pecotmr"
    )
    msgs <- capture_messages(
        summaryStatsQc(
            ss,
            alleleFlipKriging = FALSE,
            ldMismatchQcMethod = "none",
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

test_that("summaryStatsQc: QtlSumStats log lines carry the tuple label", {
    # Reuse the QtlSumStats fixture from the round-trip test.
    qss <- QtlSumStats(
        studyName = "qstudy",
        context = "qctx",
        trait = "qtrait",
        entry = list(.ssQ_makeEntryGr()),
        genome = "hg19",
        ldReference = .ssQ_makeHandle()
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
        path = "/tmp/sketch",
        format = "plink2",
        snpInfo = data.frame(
            SNP = paste0("rs", seq_len(snp_n)),
            CHR = rep("chr1", snp_n),
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
            rbinom(n_samples * nrow(snpInfo(handle)), 2, 0.3),
            nrow = n_samples,
            ncol = nrow(snpInfo(handle)),
            dimnames = list(sampleIds(handle), snpInfo(handle)$SNP)
        )
        sub <- panel[, snpIdx, drop = FALSE]
        rr <- GenomicRanges::GRanges(
            seqnames = paste0("chr", snpInfo(handle)$CHR[snpIdx]),
            ranges = IRanges::IRanges(
                start = snpInfo(handle)$BP[snpIdx],
                width = 1L
            )
        )
        S4Vectors::mcols(rr) <- S4Vectors::DataFrame(
            SNP = snpInfo(handle)$SNP[snpIdx],
            A1 = snpInfo(handle)$A1[snpIdx],
            A2 = snpInfo(handle)$A2[snpIdx]
        )
        cd <- S4Vectors::DataFrame(
            sampleId = sampleIds(handle),
            row.names = sampleIds(handle)
        )
        dosage <- t(sub)
        rownames(dosage) <- snpInfo(handle)$SNP[snpIdx]
        colnames(dosage) <- sampleIds(handle)
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
    expect_equal(df$chrom, rep("chr1", 3))
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
# .refVariantsFromLdRef
# ===========================================================================

test_that(".refVariantsFromLdRef: extracts chr/pos/alleles from snpInfo", {
    h <- .ssh_makeHandle()
    rv <- pecotmr:::.refVariantsFromLdRef(h)
    expect_equal(rv$chrom, rep("chr1", 6))
    expect_equal(rv$pos, c(100L, 200L, 300L, 400L, 500L, 600L))
    expect_equal(rv$variant_id, paste0("rs", 1:6))
    expect_equal(rv$A1, rep("A", 6))
    expect_equal(rv$A2, rep("G", 6))
})

# ===========================================================================
# .applySanityChecks
# ===========================================================================

test_that(".applySanityChecks: empty input is a no-op", {
    out <- pecotmr:::.applySanityChecks(data.frame(chrom = character()))
    expect_equal(nrow(out$df), 0L)
    expect_equal(out$audit, list())
})

test_that(".applySanityChecks: coerceNumeric converts and counts NAs", {
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

test_that(".applySanityChecks: normalizeChr maps 23/24/M, drops others", {
    df <- data.frame(
        chrom = c("chr1", "23", "24", "M", "chrX_random"),
        pos = c(100L, 200L, 300L, 400L, 500L),
        A1 = c("A", "A", "A", "A", "A"),
        A2 = c("G", "G", "G", "G", "G"),
        Z = c(1, 2, 3, 4, 5),
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.applySanityChecks(df)
    expect_equal(out$df$chrom, c("chr1", "chrX", "chrY", "chrM"))
    expect_equal(out$audit$nonstandardChrDropped, 1L)
})

test_that(".applySanityChecks: dropMissData drops NA in vital columns", {
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

test_that("summaryStatsQc: surfaces the sanity audit, honours the knobs", {
    gr <- .ssQ_makeEntryGr()
    S4Vectors::mcols(gr)$BETA <- c(0.1, 0, 0.2, 0)
    S4Vectors::mcols(gr)$SE <- c(0.1, 0.1, 0.1, 0.1)
    S4Vectors::mcols(gr)$P <- c(0.5, 0.5, 0.5, 0.5)
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ssQ_makeHandle()
    )
    res <- summaryStatsQc(ss)
    ea <- res$qcAudit[[1L]]
    expect_equal(ea$sanityChecks$zeroEffectDropped, 2L)
    # Disable the knob: zero-effect rows should remain.
    res2 <- summaryStatsQc(
        ss,
        sumstatsCleaningParam = SumstatsCleaningParam(dropZeroEffect = FALSE)
    )
    ea2 <- res2$qcAudit[[1L]]
    expect_null(ea2$sanityChecks$zeroEffectDropped)
})

# ===========================================================================
# .applySkipRegion
# ===========================================================================

.ssh_smallDf <- function() {
    data.frame(
        chrom = c("chr1", "chr1", "chr2"),
        pos = c(100L, 200L, 100L),
        SNP = c("rs1", "rs2", "rs3"),
        Z = c(1, 2, 3),
        stringsAsFactors = FALSE
    )
}

test_that(".applySkipRegion: NULL / empty skipRegion is a no-op", {
    df <- .ssh_smallDf()
    expect_identical(pecotmr:::.applySkipRegion(df, NULL), df)
    expect_identical(
        pecotmr:::.applySkipRegion(df, GenomicRanges::GRanges()),
        df
    )
})

test_that(".applySkipRegion: handles multiple regions", {
    df <- .ssh_smallDf()
    out <- pecotmr:::.applySkipRegion(
        df,
        GenomicRanges::GRanges(c("chr1:50-250", "chr2:50-150"))
    )
    expect_equal(nrow(out), 0L)
})

test_that(".applySkipRegion: accepts a GRanges of skip regions", {
    df <- .ssh_smallDf()
    # A window spelled "1" is read as chr1, as every range argument is.
    gr <- GenomicRanges::GRanges("1", IRanges::IRanges(start = 50, end = 150))
    out <- pecotmr:::.applySkipRegion(df, gr)
    expect_equal(out$SNP, c("rs2", "rs3"))
})

test_that(".applySkipRegion: refuses anything but a GRanges", {
    df <- .ssh_smallDf()
    # A window given in R is a GRanges; a string is refused, not read.
    expect_error(
        pecotmr:::.applySkipRegion(df, "1:50-150"),
        "`skipRegion` must be a GRanges \\(got character\\)"
    )
    expect_error(
        pecotmr:::.applySkipRegion(df, 42L),
        "`skipRegion` must be a GRanges"
    )
})

# ===========================================================================
# .matchAgainstLdRef
# ===========================================================================

test_that(".matchAgainstLdRef: errors when neither Z nor BETA is present", {
    df <- data.frame(
        chrom = "1",
        pos = 100L,
        SNP = "rs1",
        A1 = "A",
        A2 = "G",
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::.matchAgainstLdRef(df, .ssh_makeHandle(), matchMinProp = 0),
        "must contain at least one of Z or BETA"
    )
})

test_that(".matchAgainstLdRef: errors when A1/A2 columns are missing", {
    df <- data.frame(
        chrom = "1",
        pos = 100L,
        SNP = "rs1",
        Z = 1.0,
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::.matchAgainstLdRef(df, .ssh_makeHandle(), matchMinProp = 0),
        "must contain A1 and A2 columns"
    )
})

test_that(".matchAgainstLdRef: harmonizes the input against the sketch", {
    df <- data.frame(
        chrom = c("1", "1"),
        pos = c(100L, 200L),
        SNP = c("rs1", "rs2"),
        A1 = c("A", "A"),
        A2 = c("G", "G"),
        Z = c(1.0, 2.0),
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.matchAgainstLdRef(
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

test_that(".applyLdMismatchQcToEntry: drops variants absent from the sketch", {
    # A variant can survive study-side QC while its LD-panel partner is removed
    # by the panel MAF/MAC/missingness filter, so it is absent from the sketch
    # the QC lookup sees. .panelVariantFilter passes such a variant through and
    # leaves the drop to .ldFromRef's onMissing, so the z-mismatch lookup
    # drops it (recording the count) instead of aborting the run.
    local_mocked_bindings(
        extractBlockGenotypes = .ssh_mockExtractor(),
        .package = "pecotmr"
    )
    df <- data.frame(
        SNP = c(paste0("rs", 1:6), "ghost"),
        Z = c(seq(1.0, by = 0.5, length.out = 6L), 2.0),
        N = rep(1000L, 7L),
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.applyLdMismatchQcToEntry(
        df,
        .ssh_makeHandle(snp_n = 6L),
        method = "slalom"
    )
    expect_false("ghost" %in% out$df$SNP)
    expect_equal(out$panelUnsupportedDropped, 1L)
    expect_true(all(out$df$SNP %in% paste0("rs", 1:6)))
})

test_that(".qcKrigingFlip: drops a panel-unsupported variant, no abort", {
    # Same contract at the kriging prefilter: the orphan (no sketch entry) is
    # dropped and counted, and the run continues on the panel-supported subset.
    local_mocked_bindings(
        extractBlockGenotypes = .ssh_mockExtractor(),
        .package = "pecotmr"
    )
    df <- data.frame(
        SNP = c(paste0("rs", 1:6), "ghost"),
        Z = c(seq(1.0, by = 0.5, length.out = 6L), 2.0),
        N = rep(1000L, 7L),
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.qcKrigingFlip(
        df,
        .ssh_makeHandle(snp_n = 6L),
        opts = list(alleleFlipKriging = TRUE, nForPip = NULL),
        lbl = NA_character_
    )
    expect_false("ghost" %in% out$df$SNP)
    expect_equal(out$audit$panelUnsupportedDropped, 1L)
    expect_equal(nrow(out$df), 6L)
})

test_that(".qcKrigingFlip: no-orphan entry is unchanged (drop is a no-op)", {
    # When every variant has a sketch entry, nothing is dropped and the result
    # matches the pre-change behaviour (panelUnsupportedDropped == 0).
    local_mocked_bindings(
        extractBlockGenotypes = .ssh_mockExtractor(),
        .package = "pecotmr"
    )
    df <- data.frame(
        SNP = paste0("rs", 1:6),
        Z = seq(1.0, by = 0.5, length.out = 6L),
        N = rep(1000L, 6L),
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.qcKrigingFlip(
        df,
        .ssh_makeHandle(snp_n = 6L),
        opts = list(alleleFlipKriging = TRUE, nForPip = NULL),
        lbl = NA_character_
    )
    expect_equal(out$audit$panelUnsupportedDropped, 0L)
    expect_equal(nrow(out$df), 6L)
    expect_setequal(out$df$SNP, paste0("rs", 1:6))
})

test_that("the kriging and mismatch checks read an entry's panel once", {
    # Both run on the same panel-supported variants, so the reference is read
    # once for the two, and the panel's drop is recorded once.
    local_mocked_bindings(
        extractBlockGenotypes = .ssh_mockExtractor(),
        .package = "pecotmr"
    )
    reads <- 0L
    real <- pecotmr:::.qcPanelSupportedLd
    local_mocked_bindings(
        .qcPanelSupportedLd = function(...) {
            reads <<- reads + 1L
            real(...)
        },
        .package = "pecotmr"
    )
    df <- data.frame(
        SNP = c(paste0("rs", 1:6), "ghost"),
        Z = c(seq(1.0, by = 0.5, length.out = 6L), 2.0),
        N = rep(1000L, 7L),
        stringsAsFactors = FALSE
    )
    state <- list(
        df = df,
        ldReference = .ssh_makeHandle(snp_n = 6L),
        opts = list(
            alleleFlipKriging = TRUE,
            nForPip = NULL,
            ldMismatchQcMethod = "slalom"
        ),
        lbl = NA_character_,
        entryAudit = list(),
        qcCount = list()
    )
    out <- pecotmr:::.qcPhaseLdChecks(state)
    expect_identical(reads, 1L)
    expect_equal(out$entryAudit$panelUnsupportedDropped, 1L)
    expect_false("ghost" %in% out$df$SNP)
})

test_that("fine-mapping keeps onMissing='error' for an absent variant", {
    # The post-QC invariant: QC now drops panel-unsupported variants, so every
    # variant reaching fine-mapping has a panel entry. .ldFromRef's default
    # onMissing='error' (used by the fine-mapping LD build) MUST still abort on
    # a genuinely absent variant, guaranteeing that invariant.
    expect_error(
        pecotmr:::.ldFromRef(
            .ssh_makeHandle(snp_n = 6L),
            c("rs1", "ghost"),
            label = "fine-mapping"
        ),
        "not present in the LD reference"
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
        path = "/tmp/sketch",
        format = "plink2",
        snpInfo = data.frame(
            SNP = c(
                "chr1:100:A:G",
                "chr1:200:A:G",
                "chr1:200:G:A",
                "chr1:300:C:T"
            ),
            CHR = rep("chr1", 4L),
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
        sumstatsFilterParam = SumstatsFilterParam(
            removeIndels = FALSE,
            removeStrandAmbiguous = TRUE
        ),
        alleleFlipKriging = TRUE,
        nForPip = 1000,
        ldMismatchQcMethod = "slalom"
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
    R <- pecotmr:::.ldFromRef(
        handle,
        as.character(harm$df$SNP),
        label = "fineMappingPipeline"
    )
    expect_equal(dim(R), c(2L, 2L))
})

test_that(".applyLdMismatchQcToEntry: NA slalom outlier flags are kept", {
    # Regression test: slalom (and dentist on degenerate inputs) can leave
    # NA in the `outlier` column for variants whose per-variant statistic
    # is undefined. Treating NA as TRUE would silently drop those rows;
    # treating it as FALSE (current behavior) keeps them and yields a
    # finite outlier count.
    handle <- .ssh_makeHandle()
    panel_ids <- as.character(snpInfo(handle)$SNP)
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

test_that(".asScreen canonicalizes screen specs", {
    expect_null(pecotmr:::.asScreen(NULL))
    expect_null(pecotmr:::.asScreen(0))
    expect_null(pecotmr:::.asScreen(c(0.1, 0.2))) # non-scalar numeric -> off
    expect_equal(pecotmr:::.asScreen(0.9), list(metric = "pip", cutoff = 0.9))
    sc <- list(metric = "logBf", cutoff = 3)
    expect_identical(pecotmr:::.asScreen(sc), sc)
    expect_null(pecotmr:::.asScreen(list(metric = "bf", cutoff = 0))) # explicit off
})

test_that(".applyEntryScreen: absZ skips / retains on max|Z|, no fit", {
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
    gr <- .ssQ_makeEntryGr()
    S4Vectors::mcols(gr)$Z <- rep(0.1, length(gr))
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ssQ_makeHandle()
    )
    screens <- list(
        absZ = SignalScreenParam(absZ = 5),
        bf = SignalScreenParam(bf = 100),
        logBf = SignalScreenParam(logBf = 5)
    )
    for (metric in names(screens)) {
        res <- summaryStatsQc(
            ss,
            signalScreenParam = screens[[metric]],
            sumstatsFilterParam = SumstatsFilterParam(nCutoff = 0)
        )
        ea <- res$qcAudit[[1L]]
        expect_true(isTRUE(ea$pipScreenSkipped), info = metric)
        expect_equal(length(res[[1L]]), 0L, info = metric)
    }
})

test_that("summaryStatsQc: absZ screen keeps a strong marginal Z", {
    gr <- .ssQ_makeEntryGr()
    z <- rep(0.1, length(gr))
    z[1] <- 8
    S4Vectors::mcols(gr)$Z <- z
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ssQ_makeHandle()
    )
    res <- summaryStatsQc(
        ss,
        signalScreenParam = SignalScreenParam(absZ = 5),
        sumstatsFilterParam = SumstatsFilterParam(nCutoff = 0)
    )
    ea <- res$qcAudit[[1L]]
    expect_false(isTRUE(ea$pipScreenSkipped))
    expect_gt(length(res[[1L]]), 0L)
})

test_that("summaryStatsQc: enabling two screens at once errors", {
    ss <- .ssQ_makeGwasSumStats()
    expect_error(
        summaryStatsQc(
            ss,
            signalScreenParam = SignalScreenParam(pip = 0.5, absZ = 5)
        ),
        "only one screening metric"
    )
})

context("dentist_qc")
library(MASS)
library(corpcor)

generateDentistData <- function(
    seed = 42,
    nSnps = 100,
    sample_size = 100,
    n_outliers = 5,
    start_pos = 1000000,
    end_pos = 4000000
) {
    set.seed(seed)
    cor_matrix <- matrix(0, nrow = nSnps, ncol = nSnps)
    for (i in seq_len(nSnps - 1)) {
        for (j in (i + 1):nSnps) {
            cor_matrix[i, j] <- runif(1, 0.2, 0.8)
            cor_matrix[j, i] <- cor_matrix[i, j]
        }
    }
    diag(cor_matrix) <- 1
    ld_matrix <- cov2cor(make.positive.definite(cor_matrix))
    z_scores <- mvrnorm(n = 1, mu = rep(0, nSnps), Sigma = ld_matrix)
    outlier_indices <- sample(seq_len(nSnps), n_outliers)
    z_scores[outlier_indices] <- rnorm(n_outliers, mean = 0, sd = 5)
    sumstat <- data.frame(
        position = map_dbl(seq(start_pos, end_pos, length.out = nSnps), round),
        z = z_scores
    )
    return(list(sumstat = sumstat, ldMat = ld_matrix, nSample = sample_size))
}

generateDentistSingleWindowData <- function(
    seed = 42,
    nSnps = 100,
    sample_size = 100,
    n_outliers = 5
) {
    set.seed(seed)
    cor_matrix <- matrix(0, nrow = nSnps, ncol = nSnps)
    for (i in seq_len(nSnps - 1)) {
        for (j in (i + 1):nSnps) {
            cor_matrix[i, j] <- runif(1, 0.2, 0.8)
            cor_matrix[j, i] <- cor_matrix[i, j]
        }
    }
    diag(cor_matrix) <- 1
    ld_matrix <- cov2cor(make.positive.definite(cor_matrix))
    z_scores <- mvrnorm(n = 1, mu = rep(0, nSnps), Sigma = ld_matrix)
    outlier_indices <- sample(seq_len(nSnps), n_outliers)
    z_scores[outlier_indices] <- rnorm(n_outliers, mean = 0, sd = 5)
    return(list(z_scores = z_scores, ldMat = ld_matrix, nSample = sample_size))
}

# ===========================================================================
# dentist: basic tests
# ===========================================================================

test_that("dentist returns one row per variant under each DentistParam", {
    data <- generateDentistData(nSnps = 100)
    params <- list(
        default = DentistParam(),
        correctChenEtAlBug = DentistParam(correctChenEtAlBug = FALSE),
        gcControl = DentistParam(gcControl = TRUE)
    )
    for (nm in names(params)) {
        expect_warning(
            res <- dentist(
                data$sumstat,
                R = data$ldMat,
                nSample = data$nSample,
                dentistParam = params[[nm]]
            ),
            info = nm
        )
        expect_equal(nrow(res), 100, info = nm)
    }
})

test_that("dentist stops when the position or z column is missing", {
    data <- generateDentistData()
    for (cols in list(c("something", "z"), c("position", "something"))) {
        sumstat <- `colnames<-`(data$sumstat, cols)
        expect_error(
            dentist(
                sumstat,
                R = data$ldMat,
                nSample = data$nSample,
                dentistParam = DentistParam(correctChenEtAlBug = FALSE)
            ),
            regexp = "missing either.*pos.*or.*z",
            info = str_flatten(cols, " / ")
        )
    }
})

test_that("dentist accepts 'position' and 'zscore' column names", {
    set.seed(42)
    nSnps <- 80
    n_samples <- 100
    cor_matrix <- matrix(0, nrow = nSnps, ncol = nSnps)
    for (i in seq_len(nSnps - 1)) {
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

# Genotypes with strong LD clusters, a near-duplicate pair and a monomorphic
# variant, plus z-scores with two planted sign errors.
.dentistFactorData <- function(n = 250, p = 300, seed = 21) {
    set.seed(seed)
    base <- matrix(rnorm(n * 30), n)
    X <- base[, sample(30, p, replace = TRUE)] +
        matrix(rnorm(n * p, sd = 0.8), n)
    X[, 7] <- X[, 6] + rnorm(n, sd = 0.01)
    X[, 40] <- 1
    R <- computeLd(X)
    z <- as.vector(R %*% c(rnorm(6) * 3, rep(0, p - 6)) + rnorm(p))
    z[c(90, 200)] <- -3 * z[c(90, 200)]
    list(X = X, R = R, z = z, n = n)
}

test_that("dentistSingleWindow from genotypes matches DENTIST on their LD", {
    # Genotypes go to DENTIST's factor kernel, which never forms the LD; the
    # duplicate scan reads the factor a block of rows at a time.
    d <- .dentistFactorData()
    viaX <- suppressWarnings(dentistSingleWindow(d$z, X = d$X))
    viaR <- suppressWarnings(dentistSingleWindow(d$z, R = d$R, nSample = d$n))
    expect_equal(nrow(viaX), length(d$z))
    expect_equal(viaX, viaR, tolerance = 1e-8)
})

test_that("dentist from genotypes matches DENTIST on their LD per window", {
    d <- .dentistFactorData()
    sumStat <- data.frame(pos = seq_along(d$z) * 1000, z = d$z)
    viaX <- suppressWarnings(
        dentist(sumStat, X = d$X, windowMode = "count", minDim = 120)
    )
    viaR <- suppressWarnings(dentist(
        sumStat,
        R = d$R,
        nSample = d$n,
        windowMode = "count",
        minDim = 120
    ))
    expect_equal(viaX, viaR, tolerance = 1e-8)
})

test_that("the duplicate scan reads a factor as it reads the LD", {
    d <- .dentistFactorData()
    f <- pecotmr:::.ldFactorFromDosage(d$X)
    viaR <- pecotmr:::.findDuplicateVariants(d$z, d$R, 0.995)
    viaF <- pecotmr:::.findDuplicateVariantsFactor(d$z, f, 0.995)
    expect_gt(sum(viaR$dupBearer != -1), 0L)
    for (k in c("filteredZ", "dupBearer", "corABS", "sign", "minValue")) {
        expect_equal(viaF[[k]], viaR[[k]], tolerance = 1e-10)
    }
    kept <- viaF$filteredFactor
    expect_equal(
        unname(crossprod(kept) / (nrow(kept) - 1)),
        unname(viaR$filteredLD)
    )
})

# ===========================================================================
# dentistSingleWindow
# ===========================================================================

test_that("dentistSingleWindow returns one row per z-score under each option", {
    # Every option set keeps the input z, reports the documented columns and
    # an outlier_stat of (z - imputed)^2 / (1 - rsq); fewer than 2000
    # variants also draw DENTIST's size warning.
    data <- generateDentistSingleWindowData()
    options <- list(
        default = list(),
        correctChenEtAlBug = list(correctChenEtAlBug = FALSE),
        gcControl = list(gcControl = TRUE)
    )
    expectedCols <- c(
        "original_z",
        "imputed_z",
        "iter_to_correct",
        "rsq",
        "is_duplicate",
        "outlier_stat",
        "outlier"
    )
    for (nm in names(options)) {
        expect_warning(
            res <- exec(
                dentistSingleWindow,
                data$z_scores,
                R = data$ldMat,
                nSample = data$nSample,
                !!!options[[nm]]
            ),
            "below 2000",
            info = nm
        )
        expect_equal(nrow(res), 100, info = nm)
        expect_true(all(expectedCols %in% colnames(res)), info = nm)
        expect_type(res$outlier, "logical")
        expect_equal(res$original_z, data$z_scores, info = nm)
        expect_equal(
            res$outlier_stat,
            (res$original_z - res$imputed_z)^2 / pmax(1 - res$rsq, 1e-8),
            tolerance = 1e-10,
            info = nm
        )
    }
})

test_that("dentistSingleWindow stops on a zscore/LD dimension mismatch", {
    data <- generateDentistSingleWindowData()
    expect_warning(expect_error(
        dentistSingleWindow(
            generateDentistSingleWindowData()$z_scores,
            R = generateDentistSingleWindowData(nSnps = 80)$ldMat,
            nSample = data$nSample
        ),
        regexp = "ldMat must be a square matrix"
    ))
})

# ===========================================================================
# addDupsBackDentist
# ===========================================================================

test_that("addDupsBackDentist works", {
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

# ===========================================================================
# segmentByDist
# ===========================================================================

test_that("segmentByDist fill regions cover all input positions", {
    # Verify that fill regions cover every SNP index exactly once
    pos <- seq(1000000, 5000000, length.out = 200)
    res <- pecotmr:::segmentByDist(pos, maxDist = 2000000, minDim = 10)
    # Collect all fill region indices
    covered <- integer(0)
    for (k in seq_len(nrow(res))) {
        covered <- c(covered, res$fillStartIdx[k]:(res$fillEndIdx[k] - 1L))
    }
    # Every position from 1 to length(pos) should be covered
    expect_equal(sort(unique(covered)), seq_along(pos))
})

test_that("segmentByDist errors on empty positions", {
    expect_error(pecotmr:::segmentByDist(integer(0)), "No positions")
})

test_that("segmentByDist verbose mode prints intervals", {
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
# detectGaps
# ===========================================================================

test_that("detectGaps finds no internal gaps for contiguous positions", {
    pos <- seq(1000, by = 100, length.out = 50)
    gaps <- pecotmr:::detectGaps(pos, gapThreshold = 500)
    # Only start and end sentinel
    expect_equal(gaps, c(1L, 51L))
})

test_that("detectGaps finds multiple gaps", {
    pos <- c(1000, 2000, 5000000, 6000000, 12000000)
    gaps <- pecotmr:::detectGaps(pos, gapThreshold = 1e6)
    # Gaps at positions 3 and 5 (diffs > 1e6 at indices 2 and 4)
    expect_equal(gaps, c(1L, 3L, 5L, 6L))
})

test_that("detectGaps verbose branch prints messages", {
    pos <- c(1000, 2000, 5000000, 6000000)
    expect_message(
        pecotmr:::detectGaps(pos, gapThreshold = 1e6, verbose = TRUE),
        "No\\. of gaps found"
    )
})

# ===========================================================================
# segmentByCount
# ===========================================================================

test_that("segmentByCount tiles one block with valid, covering windows", {
    # 200 variants with maxCount = 100 make a few overlapping windows.
    pos <- seq(1000000, by = 1000, length.out = 200)
    res <- pecotmr:::segmentByCount(pos, maxCount = 100)
    expect_true(nrow(res) >= 2)
    expect_true(nrow(res) <= 5)
    # Window indices stay in range (end indices are 1-based exclusive) and
    # each fill region sits inside its window.
    expect_true(all(res$windowStartIdx >= 1))
    expect_true(all(res$windowEndIdx <= length(pos) + 1))
    expect_true(all(res$fillStartIdx >= res$windowStartIdx))
    expect_true(all(res$fillEndIdx <= res$windowEndIdx))
    # The fill regions tile every position.
    covered <- integer(0)
    for (k in seq_len(nrow(res))) {
        covered <- c(covered, res$fillStartIdx[k]:(res$fillEndIdx[k] - 1L))
    }
    expect_equal(sort(unique(covered)), seq_along(pos))
})

test_that("segmentByCount handles centromeric gap", {
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
    for (k in seq_len(nrow(res))) {
        covered <- c(covered, res$fillStartIdx[k]:(res$fillEndIdx[k] - 1L))
    }
    expect_equal(sort(unique(covered)), seq_along(pos))
})

test_that("segmentByCount skips blocks smaller than half maxCount", {
    # 20 variants with maxCount = 100 (half = 50): too small, so skipped.
    pos <- seq(1000000, by = 1000, length.out = 20)
    expect_error(
        pecotmr:::segmentByCount(pos, maxCount = 100),
        "No intervals created by segmentation"
    )
})

test_that("segmentByCount creates single window for small blocks", {
    # 60 variants with maxCount = 100 make one window (60 >= half = 50).
    pos <- seq(1000000, by = 1000, length.out = 60)
    res <- pecotmr:::segmentByCount(pos, maxCount = 100)
    expect_equal(nrow(res), 1)
    expect_equal(res$windowStartIdx[1], 1)
    expect_equal(res$windowEndIdx[1], 61)
})

test_that("segmentByCount errors on empty positions", {
    expect_error(
        pecotmr:::segmentByCount(integer(0), maxCount = 100),
        "No positions"
    )
})

test_that("segmentByCount verbose mode prints intervals", {
    pos <- seq(1000000, by = 1000, length.out = 300)
    expect_message(
        pecotmr:::segmentByCount(pos, maxCount = 100, verbose = TRUE),
        "Intervals"
    )
})

# ===========================================================================
# mergeWindows
# ===========================================================================

test_that("mergeWindows returns exactly N rows", {
    # 300 variants over 2 kb with maxDist = 1000 make four overlapping
    # windows; the LD build is cubic in the variant count, so no more.
    nSnps <- 300L
    data <- generateDentistData(
        nSnps = nSnps,
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
        for (k in seq_len(nrow(window_divided_res))) {
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
                BPPARAM = BiocParallel::SerialParam(),
                correctChenEtAlBug = TRUE
            )
        }
    })
    res <- pecotmr:::mergeWindows(dentist_result_by_window, window_divided_res)
    expect_equal(nrow(res), nSnps)
})

test_that("mergeWindows stops with window and imputed mismatch", {
    expect_error(
        pecotmr:::mergeWindows(
            rep(0, 5),
            data.frame(windowStartIdx = rep(0, 2), windowEndIdx = rep(0, 2))
        ),
        "Different number of windows and imputed results"
    )
})

test_that("mergeWindows correctly indexes and merges windows", {
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
    # 300 variants over 5 Mb with windowSize = 2e6 make five overlapping
    # windows of more than minDim variants each; the LD build is cubic in
    # the variant count, so no more.
    nSnps <- 300L
    data <- generateDentistData(
        seed = 123,
        nSnps = nSnps,
        sample_size = 1000,
        n_outliers = 15,
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
    expect_equal(nrow(res), nSnps)
})

# ===========================================================================
# dentist with count mode
# ===========================================================================

test_that("dentist with window_mode='count' returns exactly N rows", {
    data <- generateDentistData(
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
# resolveLdInput (internal)
# ===========================================================================

test_that("resolveLdInput passes R through, with or without its nSample", {
    R <- diag(3)
    result <- pecotmr:::resolveLdInput(R = R, needNSample = FALSE)
    expect_equal(result$R, R)
    expect_null(result$nSample)
    result <- pecotmr:::resolveLdInput(
        R = R,
        nSample = 100L,
        needNSample = TRUE
    )
    expect_equal(result$R, R)
    expect_equal(result$nSample, 100L)
})

test_that("resolveLdInput computes the LD and nSample from genotypes", {
    set.seed(42)
    n <- 50
    p <- 5
    X <- matrix(rbinom(n * p, 2, 0.3), nrow = n, ncol = p)
    result <- pecotmr:::resolveLdInput(X = X, needNSample = TRUE)
    expect_equal(result$nSample, n)
    expect_true(is.matrix(result$R))
    expect_equal(dim(result$R), c(p, p))
    expect_true(isSymmetric(result$R))
    expect_equal(diag(result$R), rep(1, p), tolerance = 1e-6)
})

test_that("resolveLdInput converts a non-matrix X to a matrix", {
    set.seed(42)
    X_df <- data.frame(a = rbinom(30, 2, 0.3), b = rbinom(30, 2, 0.3))
    result <- pecotmr:::resolveLdInput(X = X_df, needNSample = FALSE)
    expect_true(is.matrix(result$R))
    expect_equal(result$nSample, 30)
})

test_that("resolveLdInput keeps an explicit nSample given with X", {
    set.seed(42)
    X <- matrix(rbinom(100, 2, 0.3), nrow = 20, ncol = 5)
    result <- pecotmr:::resolveLdInput(X = X, nSample = 999, needNSample = TRUE)
    expect_equal(result$nSample, 999)
})

test_that("resolveLdInput needs exactly one of R and X", {
    expect_error(
        pecotmr:::resolveLdInput(),
        "Either R.*or X.*must be provided"
    )
    expect_error(
        pecotmr:::resolveLdInput(R = diag(3), X = matrix(1:9, nrow = 3)),
        "Provide either R or X, not both"
    )
})

test_that("resolveLdInput needs nSample with R when one is required", {
    expect_error(
        pecotmr:::resolveLdInput(R = diag(3), needNSample = TRUE),
        "nSample is required"
    )
})

# ===========================================================================
# buildSegmentResult (internal)
# ===========================================================================

test_that("buildSegmentResult caps end indices and verbose prints", {
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

test_that("buildSegmentResult errors on empty startList", {
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
# slidingWindowLoop (iteration limit)
# ===========================================================================

test_that("slidingWindowLoop errors on infinite loop", {
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

# ===========================================================================
# summaryStatsQc: study-level case/control counts
# ===========================================================================

test_that("summaryStatsQc: preserves nCase/nControl columns through QC", {
    gr <- .ssQ_makeEntryGr(paste0("rs", 1:4), c(100L, 200L, 300L, 400L))
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ssQ_makeHandle(),
        nCase = 500,
        nControl = 1500
    )
    expect_true(all(c("nCase", "nControl") %in% colnames(S4Vectors::mcols(ss))))
    out <- summaryStatsQc(
        ss,
        signalScreenParam = SignalScreenParam(pip = 0),
        sumstatsFilterParam = SumstatsFilterParam(nCutoff = 0)
    )
    expect_true(all(
        c("nCase", "nControl") %in% colnames(S4Vectors::mcols(out))
    ))
    expect_equal(out$nCase, 500)
    expect_equal(out$nControl, 1500)
})

# ===========================================================================
# mergeVariantInfo (data.frame + GRanges, flip-aware, all = TRUE/FALSE)
# ===========================================================================

test_that("mergeVariantInfo (data.frame): all = TRUE returns the union", {
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

test_that("mergeVariantInfo (data.frame): all = FALSE returns variants2", {
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

test_that("mergeVariantInfo (GRanges): converts inputs and detects flips", {
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
# addDupsBackDentist: dimension-mismatch stops
# ===========================================================================

test_that("addDupsBackDentist stops on a dentistOutput nrow mismatch", {
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

test_that("addDupsBackDentist stops on inconsistent input lengths", {
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
# autoDecision: high-correlation tagging branch
# ===========================================================================

test_that("autoDecision evaluates the high-corr tag for a non-top CS", {
    # A non-top CS with a small p-value can only be tagged through
    # highCorrCols: a 1 in its correlation column tags it, leaving SER.
    # cs_corr_1 is each CS's correlation with L1, as extractCsInfo() lays
    # it out, so L1's own entry is 1; being the top CS, L1 is never tagged.
    df <- data.frame(
        cs_name = c("L1", "L2"),
        top_z = c(5.0, 3.0),
        p_value = c(1e-10, 1e-10),
        cs_corr_1 = c(1, 1),
        stringsAsFactors = FALSE
    )
    tagged <- autoDecision(df, highCorrCols = "cs_corr_1")
    expect_equal(tagged$tagged_cs, c(FALSE, TRUE))
    expect_equal(tagged$method, c("SER", "SER"))
    # Correlated with L1 short of 1, L2 stays untagged.
    df$cs_corr_1 <- c(1, 0.5)
    untagged <- autoDecision(df, highCorrCols = "cs_corr_1")
    expect_equal(untagged$tagged_cs, c(FALSE, FALSE))
    expect_equal(untagged$method, c("BVSR", "BVSR"))
})

test_that("autoDecision never tags a CS by its own correlation column", {
    # Every cs_corr_* column passed: L2's diagonal entry (its self-
    # correlation, always 1) must not tag it; only its correlation with
    # L1 can.
    df <- data.frame(
        cs_name = c("L1", "L2"),
        top_z = c(5.0, 3.0),
        p_value = c(1e-10, 1e-10),
        cs_corr_1 = c(1, 0.2),
        cs_corr_2 = c(0.2, 1),
        stringsAsFactors = FALSE
    )
    cols <- c("cs_corr_1", "cs_corr_2")
    untagged <- autoDecision(df, highCorrCols = cols)
    expect_equal(untagged$tagged_cs, c(FALSE, FALSE))
    expect_equal(untagged$method, c("BVSR", "BVSR"))
    df$cs_corr_1[2] <- 1
    df$cs_corr_2[1] <- 1
    expect_equal(
        autoDecision(df, highCorrCols = cols)$tagged_cs,
        c(FALSE, TRUE)
    )
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
        variant_id = paste0("rs", seq_len(p)),
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

test_that("raissSingleMatrixFromX emits no-known / no-unknown messages", {
    set.seed(7)
    p <- 5
    ref_panel <- data.frame(
        chrom = rep(1, p),
        pos = seq(10, p * 10, 10),
        variant_id = paste0("rs", seq_len(p)),
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
        variant_id = paste0("rs", seq_len(p)),
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
    # The inverse is checked against the ridged matrix it was computed from,
    # so a well-conditioned sigT passes.
    expect_true(res$correctInversion)
})

test_that("raissModel flags an inversion the rcond cutoff truncated", {
    # Ridged, this sigT has eigenvalues 2.01 and 0.01; their ratio is under
    # rcond = 0.01, so the pseudo-inverse drops one and cannot reproduce it.
    sig_t <- matrix(1, nrow = 2, ncol = 2)
    sig_i_t <- matrix(c(0.5, 0.2, 0.1, 0.3), nrow = 2)
    res <- pecotmr:::raissModel(c(1.2, 0.5), sig_t, sig_i_t)
    expect_identical(res$correctInversion, c(FALSE, FALSE))
    ok <- pecotmr:::raissModel(c(1.2, 0.5), diag(2), sig_i_t)
    expect_identical(ok$correctInversion, c(TRUE, TRUE))
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

test_that(".applyContentFilters: FRQ normalizes to MAF via min(af, 1-af)", {
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

test_that("summaryStatsQc: emit() drops the label for an empty study id", {
    # An empty study id resolves the per-entry label to NA, exercising the
    # unlabeled emit() branch.
    ss <- GwasSumStats(
        studyName = "",
        entry = list(.ssQ_makeEntryGr()),
        genome = "hg19",
        ldReference = .ssQ_makeHandle()
    )
    msgs <- capture_messages(summaryStatsQc(
        ss,
        sumstatsFilterParam = SumstatsFilterParam(nCutoff = 0)
    ))
    joined <- paste(msgs, collapse = "")
    expect_match(joined, "QC summary:")
    # No bracketed label prefix appears on the rollup line.
    expect_true(any(grepl("^QC summary:", msgs)))
})

test_that("summaryStatsQc: the N filter emits its message and rollup", {
    gr <- .ssQ_makeEntryGr()
    S4Vectors::mcols(gr)$N <- c(1000L, 1010L, 1005L, 100000L) # last is an N outlier
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ssQ_makeHandle()
    )
    msgs <- capture_messages(
        res <- summaryStatsQc(
            ss,
            sumstatsFilterParam = SumstatsFilterParam(nCutoff = 3)
        )
    )
    joined <- paste(msgs, collapse = "")
    expect_match(joined, "MAF/INFO/N filters kept")
    expect_match(joined, "nCutoff 1")
    ea <- res$qcAudit[[1L]]
    expect_equal(ea$contentFilters$nDropped, 1L)
})

test_that("summaryStatsQc: derives BETA/SE from Z+MAF+N, records it", {
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
        studyName = "g1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ssQ_makeHandle()
    )
    res <- summaryStatsQc(
        ss,
        sumstatsFilterParam = SumstatsFilterParam(nCutoff = 0)
    )
    ea <- res$qcAudit[[1L]]
    expect_equal(ea$betaSeFromZ$nDerived, 4L)
})

test_that("summaryStatsQc: clamps tiny Z-derived P and audits it", {
    gr <- .ssQ_makeEntryGr()
    S4Vectors::mcols(gr)$Z <- c(50, 1, 2, 3) # |Z| = 50 underflows P to 0
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ssQ_makeHandle()
    )
    res <- summaryStatsQc(
        ss,
        sumstatsFilterParam = SumstatsFilterParam(nCutoff = 0)
    )
    ea <- res$qcAudit[[1L]]
    expect_true(!is.null(ea$sanityChecks$smallPClamped))
    expect_gte(ea$sanityChecks$smallPClamped, 1L)
})

test_that("summaryStatsQc: early-exits below two pre-harmonization variants", {
    gr <- .ssQ_makeEntryGr(snp_ids = "rs1", positions = 100L)
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ssQ_makeHandle()
    )
    res <- summaryStatsQc(
        ss,
        sumstatsFilterParam = SumstatsFilterParam(nCutoff = 0)
    )
    ea <- res$qcAudit[[1L]]
    expect_match(ea$earlyExit, "fewer than two variants")
    expect_equal(length(res[[1L]]), 1L)
})

test_that("summaryStatsQc: kriging QC records its audit and rollup", {
    skip_if_not(
        "kriging_rss" %in% getNamespaceExports("susieR"),
        "installed susieR has no kriging_rss"
    )
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(.ssQ_makeEntryGr(
            snp_ids = paste0("rs", 1:8),
            positions = seq(100L, by = 100L, length.out = 8L)
        )),
        genome = "hg19",
        ldReference = .ssQ_makeHandleVid(snp_n = 8L)
    )
    local_mocked_bindings(
        extractBlockGenotypes = .ssQ_mockExtractor(),
        .package = "pecotmr"
    )
    msgs <- capture_messages(
        res <- summaryStatsQc(
            ss,
            alleleFlipKriging = TRUE,
            signalScreenParam = SignalScreenParam(pip = 0),
            sumstatsFilterParam = SumstatsFilterParam(nCutoff = 0)
        )
    )
    joined <- paste(msgs, collapse = "")
    expect_match(joined, "kriging sign-flipped")
    expect_match(joined, "kriging-flip ")
    ea <- res$qcAudit[[1L]]
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

test_that("summaryStatsQc: impute assembles BETA/SE/N, median-fills N", {
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
        studyName = "g1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ssQ_makeHandle(snp_n = 8L, n_samples = 60L)
    )
    local_mocked_bindings(
        extractBlockGenotypes = .ssQ_mockExtractor(),
        .raissGenotypePath = function(
            refPanel,
            knownZscores,
            genotypeMatrix,
            p
        ) {
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
    res <- summaryStatsQc(
        ss,
        impute = TRUE,
        raissParam = RaissParam(flank = 500)
    )
    ea <- res$qcAudit[[1L]]
    expect_equal(ea$raissTotalVariants, 6L)
    expect_equal(ea$raissImputedVariants, 2L)
    mc <- S4Vectors::mcols(res[[1L]])
    expect_true(all(c("BETA", "SE", "N") %in% colnames(mc)))
    expect_false(any(is.na(mc$N))) # median-filled
})

test_that("summaryStatsQc: the rollup enumerates every removed step", {
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
        studyName = "g1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ssQ_makeHandle()
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
            sumstatsFilterParam = SumstatsFilterParam(
                infoCutoff = 0.5,
                nCutoff = 3
            ),
            panelFilterParam = PanelFilterParam(mafCutoff = 0.01)
        )
    )
    joined <- paste(msgs, collapse = "")
    # The chr99 variant lands in its own ELEMENT (entries are split by seqname
    # at construction), so the audit is now per (study, chromosome) and the
    # drop counts are summed across entries rather than read off entry 1.
    audits <- res$qcAudit
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

# ===========================================================================
# Final-coverage closeouts: a handful of branch/error-handler lines that the
# existing suite leaves untouched.
# ===========================================================================

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

test_that("segmentByDist keeps a last window that clears the cutoff", {
    # A single dense block spanning 1.25x the distance cutoff produces two
    # windows whose final window is wide enough that adjustLastFn does NOT
    # shrink it: it keeps the current startIdx unchanged.
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

test_that("summaryStatsQc kriging QC sign-flips and keeps a bad variant", {
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
        n <- length(sampleIds(handle))
        k <- length(snpIdx)
        f <- rnorm(n) # shared latent factor
        M <- exec(
            cbind,
            !!!map(seq_len(k), function(j) {
                sqrt(0.7) * f + sqrt(0.3) * rnorm(n)
            })
        )
        rr <- GenomicRanges::GRanges(
            seqnames = paste0("chr", snpInfo(handle)$CHR[snpIdx]),
            ranges = IRanges::IRanges(
                start = snpInfo(handle)$BP[snpIdx],
                width = 1L
            )
        )
        S4Vectors::mcols(rr) <- S4Vectors::DataFrame(
            SNP = snpInfo(handle)$SNP[snpIdx],
            A1 = snpInfo(handle)$A1[snpIdx],
            A2 = snpInfo(handle)$A2[snpIdx]
        )
        dosage <- t(M)
        rownames(dosage) <- snpInfo(handle)$SNP[snpIdx]
        colnames(dosage) <- sampleIds(handle)
        SummarizedExperiment::SummarizedExperiment(
            assays = list(dosage = dosage),
            rowRanges = rr,
            colData = S4Vectors::DataFrame(
                sampleId = sampleIds(handle),
                row.names = sampleIds(handle)
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
        studyName = "g1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ssQ_makeHandleVid(snp_n = 8L, n_samples = 300L)
    )
    local_mocked_bindings(
        extractBlockGenotypes = corrExtractor,
        .package = "pecotmr"
    )
    res <- suppressWarnings(
        summaryStatsQc(
            ss,
            alleleFlipKriging = TRUE,
            signalScreenParam = SignalScreenParam(pip = 0),
            sumstatsFilterParam = SumstatsFilterParam(nCutoff = 0)
        )
    )
    ea <- res$qcAudit[[1L]]
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

test_that(".subsetLdRefToRange keeps panel variants in the entry span", {
    h <- readGenotypeHandle(
        test_path("test_data/test_variants"),
        format = "plink2"
    )
    si <- snpInfo(h)
    ch <- names(sort(table(si$CHR), decreasing = TRUE))[1] # busiest chromosome
    chIdx <- which(si$CHR == ch)
    slice <- si[chIdx[seq_len(min(30L, length(chIdx)))], , drop = FALSE]
    gr <- GenomicRanges::GRanges(
        slice$CHR,
        IRanges::IRanges(as.integer(slice$BP), width = 1L)
    )
    S4Vectors::mcols(gr)$SNP <- slice$SNP
    sub <- snpInfo(pecotmr:::.subsetLdRefToRange(h, list(gr)))
    lo <- min(slice$BP)
    hi <- max(slice$BP)
    expected <- si$SNP[
        si$CHR == ch &
            as.integer(si$BP) >= lo &
            as.integer(si$BP) <= hi
    ]
    expect_setequal(sub$SNP, expected)
    expect_lt(nrow(sub), nrow(si))
})

test_that(".subsetLdRefToIds keeps exactly the entries' variants", {
    h <- readGenotypeHandle(
        test_path("test_data/test_variants"),
        format = "plink2"
    )
    si <- snpInfo(h)
    sel <- c(3L, 10L, 40L, 200L)
    gr <- GenomicRanges::GRanges(
        si$CHR[sel],
        IRanges::IRanges(as.integer(si$BP[sel]), width = 1L)
    )
    S4Vectors::mcols(gr)$SNP <- si$SNP[sel]
    sub <- snpInfo(pecotmr:::.subsetLdRefToIds(h, list(gr)))
    expect_setequal(
        normalizeVariantId(sub$SNP),
        normalizeVariantId(si$SNP[sel])
    )
})

test_that(".subsetLdRefToRange / .subsetLdRefToIds are NULL-safe", {
    expect_null(pecotmr:::.subsetLdRefToRange(NULL, list()))
    expect_null(pecotmr:::.subsetLdRefToIds(NULL, list()))
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
    ids <- as.character(snpInfo(handle)$SNP)
    altFreq <- afreq$alt_freq[match(ids, afreq$id)]
    pmin(altFreq, 1 - altFreq)
}

# A GwasSumStats over every panel variant, each carrying a study AF of 0.30 --
# common in the study whatever the panel says, so only a panel-side filter can
# remove it.
.ssqcPanelSumStats <- function(handle) {
    si <- snpInfo(handle)
    gr <- GenomicRanges::GRanges(
        seqnames = si$CHR,
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
        studyName = "g1",
        entry = list(gr),
        genome = "hg19",
        ldReference = handle
    )
}

test_that(".ssqcPrunePanel drops exactly the sub-cutoff panel variants", {
    skip_if_not_installed("pgenlibr")
    h <- .ssqcPanelHandle()
    ids <- as.character(snpInfo(h)$SNP)
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
    keptIds <- as.character(snpInfo(pruned)$SNP)
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
    # A summaryStatsQc ldReference is an RSE wrapping the GenotypeHandle;
    # harmonization, kriging, mismatch-QC and RAISS all read the panel through
    # .ldRefHandle() (the dosage DelayedArray's seed). Subsetting the RSE
    # rows leaves that seed carrying the full panel, so the prune must reach
    # the seed -- the real-flow shape the bare-handle tests do not exercise.
    h <- .ssqcPanelHandle()
    maf <- .ssqcPanelMaf(h)
    ids <- as.character(snpInfo(h)$SNP)
    cutoff <- stats::median(maf, na.rm = TRUE)
    rse <- pecotmr:::.genotypeExperiment(h)
    expect_s4_class(rse, "RangedSummarizedExperiment")
    pruned <- suppressMessages(pecotmr:::.ssqcPrunePanel(
        rse,
        list(mafCutoff = cutoff, macCutoff = 0, imissCutoff = 1),
        "test"
    ))
    seed <- pecotmr:::.ldRefHandle(pruned)
    expect_setequal(
        as.character(snpInfo(seed)$SNP),
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
    out <- suppressMessages(summaryStatsQc(
        ss,
        sumstatsFilterParam = SumstatsFilterParam(nCutoff = 0),
        panelFilterParam = PanelFilterParam(mafCutoff = cutoff)
    ))
    entry <- pecotmr:::.collectionEntry(out, 1)
    expect_equal(length(entry), sum(!is.na(maf) & maf >= cutoff))
    # ... and the panel carried on the result -- the seed handle every LD
    # read keys off, not just the row view -- holds no sub-cutoff variant.
    seed <- pecotmr:::.ldRefHandle(ldReference(out)[[1L]])
    keptMaf <- maf[match(
        as.character(snpInfo(seed)$SNP),
        as.character(snpInfo(h)$SNP)
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
    mac <- ceiling(cutoff * 2 * nSamples(h))
    out <- suppressMessages(summaryStatsQc(
        ss,
        sumstatsFilterParam = SumstatsFilterParam(nCutoff = 0),
        panelFilterParam = PanelFilterParam(macCutoff = mac)
    ))
    entry <- pecotmr:::.collectionEntry(out, 1)
    expected <- sum(!is.na(maf) & maf >= mac / (2 * nSamples(h)))
    expect_equal(length(entry), expected)
})

test_that("summaryStatsQc: the panel filter is off by default", {
    skip_if_not_installed("pgenlibr")
    h <- .ssqcPanelHandle()
    ss <- .ssqcPanelSumStats(h)
    out <- suppressMessages(summaryStatsQc(
        ss,
        sumstatsFilterParam = SumstatsFilterParam(nCutoff = 0)
    ))
    expect_equal(
        length(pecotmr:::.collectionEntry(out, 1)),
        nrow(snpInfo(h))
    )
})

test_that("summaryStatsQc validates the panel cutoffs before any panel read", {
    ss <- .ssQ_makeGwasSumStats()
    expect_error(
        summaryStatsQc(
            ss,
            panelFilterParam = PanelFilterParam(mafCutoff = -1)
        ),
        "mafCutoff"
    )
    expect_error(
        summaryStatsQc(
            ss,
            panelFilterParam = PanelFilterParam(macCutoff = c(1, 2))
        ),
        "macCutoff"
    )
    expect_error(
        summaryStatsQc(
            ss,
            panelFilterParam = PanelFilterParam(imissCutoff = NA_real_)
        ),
        "imissCutoff"
    )
    expect_error(
        summaryStatsQc(
            ss,
            panelFilterParam = PanelFilterParam(macCutoff = Inf)
        ),
        "macCutoff"
    )
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
    dosage <- exec(cbind, !!!map(af, function(f) rbinom(nS, 2L, f)))
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

# The mask as .qcRaissBuildInputs() computes it for a panel of real
# genotypes: statistics off the window's dosages, only when a cutoff is set,
# and a MAC cutoff counting the panel's samples.
# @noRd
.rmask_run <- function(f, ...) {
    rp <- as.list(RaissParam(...))
    .qcRaissTargetMask(
        f$refPanel,
        f$knownZ,
        .qcRaissTargetStats(NULL, NULL, f$dosage, rp),
        list(raissParam = rp),
        nSamples = nrow(f$dosage)
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

test_that("summaryStatsQc raissParam cutoffs bound what RAISS imputes", {
    data(gwasSumStatsS4Example)
    gss <- gwasSumStatsS4Example
    variants <- unlist(gss)
    thin <- GwasSumStats(
        studyName = studyName(gss),
        entry = list(variants[seq(1L, length(variants), by = 4L)]),
        genome = unique(unname(GenomeInfoDb::genome(gss))),
        ldReference = ldReference(gss)
    )
    nObserved <- sum(lengths(thin))
    # The bundled toy panel imputes poorly, so the R2 gate would otherwise
    # reject every target and leave nothing to count.
    base <- list(r2Threshold = 0, minimumLd = 0)
    nOut <- function(...) {
        out <- suppressWarnings(suppressMessages(summaryStatsQc(
            thin,
            impute = TRUE,
            raissParam = do.call(
                RaissParam,
                utils::modifyList(base, list(...))
            )
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

test_that("with nothing known, only the flip pair escapes imputation", {
    # The ambiguity is a property of the panel, not of what the sumstats
    # happen to carry -- and the guard must not cost coverage where there is
    # no ambiguity: exactly the unambiguous variants are imputed.
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

test_that(".raissFlipPairMask agrees with matchVariants on what is ambiguous", {
    # The mask exists to stop imputation putting back what matchVariants drops,
    # so the two have to answer the same question the same way. chr1:100's
    # A:G and C:T are one variant written on opposite strands AND swapped, so
    # matchVariants() refuses the position; a literal A1/A2 reversal key read
    # this as two unrelated variants and let RAISS impute the pair back.
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

test_that("a typed position imputes no second allele, untyped ones all forms", {
    # Imputation fills UN-observed variants; another form of a site the GWAS
    # measured is a redundant record, not new information, so 1:100:A:T is
    # no target. The guard is about MEASURED positions, not about
    # multi-allelic sites: 1:300 is multi-allelic but untyped, so both its
    # forms are targets.
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

test_that("raiss() reads either chromosome spelling as one site", {
    # "chr1" in the panel and "1" in the sumstats is one site, not two:
    # raiss() spells both tables' chromosomes in UCSC style as they come in.
    d <- .raiss_pr573()
    d$refPanel$chrom <- str_c("chr", d$refPanel$chrom)
    ld <- diag(5)
    ld[1, 2] <- ld[2, 1] <- 0.7
    ld[1, 3] <- ld[3, 1] <- 0.6
    ld[1, 4] <- ld[4, 1] <- 0.5
    ld[1, 5] <- ld[5, 1] <- 0.45
    r <- raiss(d$refPanel, d$knownZ, ldMatrix = ld, verbose = FALSE)
    imp <- setdiff(r$resultNofilter$variant_id, d$knownZ$variant_id)
    expect_false(is_in("1:100:A:T", imp))
    expect_true(all(r$resultNofilter$chrom == "chr1"))
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

test_that("both trimmers empty the sketch for zero-variant entries", {
    h <- .sk_panel()
    expect_gt(nrow(h), 0L)
    # An empty list and a zero-range GRanges both mean "no variants".
    for (entries in list(list(), list(GenomicRanges::GRanges()))) {
        expect_equal(nrow(.subsetLdRefToRange(h, entries)), 0L)
        expect_equal(nrow(.subsetLdRefToIds(h, entries)), 0L)
    }
})

test_that("emptying drops both axes but keeps the panel's identity", {
    h <- .sk_panel()
    e <- .subsetLdRefToRange(h, list())
    # The sample axis is shed along with the variants: an emptied sketch
    # references no LD, so its anonymous sample names are dead weight (the bulk
    # of a skipped-region file). What is kept is the panel's identity -- where
    # it came from -- not its dimensions.
    expect_equal(nrow(e), 0L)
    expect_equal(ncol(e), 0L)
    expect_equal(.ldRefFormat(e), .ldRefFormat(h))
    expect_equal(.ldRefSource(e)$path, .ldRefSource(h)$path)
})

test_that(".subsetLdRefToIds keeps the panel when ids are unextractable", {
    # Ranges present but no SNP mcol is a different situation from zero
    # variants: the panel cannot be keyed, so blanking it would discard LD the
    # object may well reference. Stay conservative.
    h <- .sk_panel()
    gr <- GenomicRanges::GRanges(
        "chr21",
        IRanges::IRanges(c(17513228L, 17513580L), width = 1L)
    )
    expect_equal(nrow(.subsetLdRefToIds(h, list(gr))), nrow(h))
})

# ===========================================================================
# af / MAF semantics (ported from PR #578)
#
# `af` is the DIRECTIONAL effect-allele frequency, declared only via an `af:`
# mapping key and exported as top_loci$af. A directionless `maf`/`FRQ` is
# QC-only and is never exported as af, so an undeclared study gets an honest
# NA instead of a silently mislabelled minor-allele frequency.
# ===========================================================================

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

test_that("sumstatsQc: argument guards fire", {
    expect_error(
        mergeVariantInfo(data.frame(), data.frame(), all = NA),
        "all.*May not be NA"
    )
    expect_error(
        autoDecision(data.frame(), highCorrCols = 1L),
        "highCorrCols.*Must be of type 'character'"
    )
})

test_that(".entrySnpIds answers character(0) for a NULL entry", {
    expect_identical(pecotmr:::.entrySnpIds(NULL), character(0))
})

test_that("SlalomParam / DentistParam carry defaults and reject a typo", {
    sa <- SlalomParam(r2Threshold = 0.8)
    expect_s4_class(sa, "SlalomParam")
    expect_equal(sa$r2Threshold, 0.8)
    expect_equal(sa$abfPriorVariance, 0.04)
    da <- DentistParam(propSVD = 0.5)
    expect_equal(da$propSVD, 0.5)
    expect_equal(da$nIter, 10)
    # Our own engines, so no `...`: R itself rejects an unknown name.
    expect_error(SlalomParam(zzz = 1), "unused argument")
    expect_error(DentistParam(zzz = 1), "unused argument")
})

test_that("SlalomParam leaves the data-dependent standardError to slalom", {
    # slalom defaults it to rep(1, length(zScore)), which no static default
    # can express, so a NULL default is dropped rather than forwarded.
    expect_false(is_in("standardError", names(SlalomParam())))
    expect_true(is_in("standardError", names(SlalomParam(standardError = 1))))
})

test_that("ldMismatchQc takes a name or a constructor for method", {
    set.seed(1)
    n <- 60L
    p <- 8L
    X <- matrix(rnorm(n * p), n, p)
    z <- rnorm(p)
    R <- cor(X)
    byName <- ldMismatchQc(z, R = R, nSample = n, method = "slalom")
    byCtor <- ldMismatchQc(z, R = R, nSample = n, method = SlalomParam())
    # A bare name means the engine with no options set, which is exactly what
    # an unconfigured constructor produces.
    expect_identical(byName, byCtor)
})

test_that("ldMismatchQc forwards the constructor's options to the engine", {
    set.seed(1)
    n <- 60L
    p <- 8L
    z <- rnorm(p)
    R <- cor(matrix(rnorm(n * p), n, p))
    cap <- new.env(parent = emptyenv())
    local_mocked_bindings(
        slalom = .ssq_captureSlalom(cap),
        .package = "pecotmr"
    )
    ldMismatchQc(
        z,
        R = R,
        nSample = n,
        method = SlalomParam(r2Threshold = 0.9)
    )
    expect_equal(cap$r2Threshold, 0.9)
})

test_that("ldMismatchQc rejects a constructor for another engine", {
    set.seed(1)
    z <- rnorm(4L)
    R <- diag(4L)
    expect_error(
        ldMismatchQc(z, R = R, nSample = 50L, method = PlinkQcOptions()),
        "unknown engine"
    )
})

test_that("summaryStatsQc takes a constructor for ldMismatchQcMethod", {
    expect_error(
        pecotmr:::.resolveLdMismatchChoice("nope"),
        "unknown engine 'nope'"
    )
    expect_equal(pecotmr:::.resolveLdMismatchChoice("none"), "none")
    expect_s4_class(
        pecotmr:::.resolveLdMismatchChoice(DentistParam()),
        "DentistParam"
    )
})

# =============================================================================
# RaissParam()
# =============================================================================

test_that("RaissParam carries every field the QC imputation path reads", {
    a <- RaissParam()
    expect_s4_class(a, "RaissParam")
    expect_setequal(
        names(a),
        c(
            "lamb",
            "rcond",
            "svdTol",
            "r2Threshold",
            "minimumLd",
            "mafCutoff",
            "macCutoff",
            "imissCutoff",
            "flank"
        )
    )
})

test_that("RaissParam carries rcond for an LdMatrix reference", {
    # A genotype panel reaches RAISS's SVD route (svdTol); an LdMatrix
    # reaches its LD-matrix route, which takes rcond. Both are forwarded, and
    # each route reads its own.
    expect_equal(RaissParam(rcond = 0.05)$rcond, 0.05)
    expect_equal(RaissParam()$rcond, formals(raiss)$rcond)
})

test_that("summaryStatsQc refuses a bare list for raissParam", {
    data(gwasSumStatsS4Example)
    expect_error(
        summaryStatsQc(gwasSumStatsS4Example, raissParam = list(flank = 500)),
        "must be built with RaissParam"
    )
})

test_that("an empty raissParam resolves to the full default bundle", {
    resolved <- pecotmr:::.ssqcResolveImputeArgs(list())
    expect_equal(resolved, as.list(RaissParam()))
    # Every read site indexes this without a fallback, so a missing field
    # would surface as NULL rather than a default.
    expect_false(any(map_lgl(resolved, is.null)))
})

test_that("the forwarded half of RaissParam is exactly raiss's own formals", {
    fwd <- pecotmr:::.raissForwardedNames()
    expect_setequal(
        fwd,
        c("lamb", "rcond", "svdTol", "r2Threshold", "minimumLd")
    )
    expect_true(all(fwd %in% names(formals(raiss))))
    # The other half scopes what RAISS is asked to impute and is consumed
    # before the call; splicing it into raiss would be an error.
    rest <- setdiff(names(RaissParam()), fwd)
    expect_false(any(rest %in% names(formals(raiss))))
    expect_setequal(
        names(pecotmr:::.raissForwardedArgs(RaissParam())),
        fwd
    )
})

test_that("SumstatsCleaningParam is separate from SumstatsFilterParam", {
    # Cleaning decides whether a row is a well-formed record; filtering
    # applies quality thresholds to rows that already are. No field overlaps.
    expect_length(
        intersect(
            names(formals(SumstatsCleaningParam)),
            names(formals(SumstatsFilterParam))
        ),
        0L
    )
})

test_that("every SumstatsCleaningParam field is one the resolver hands on", {
    expect_setequal(
        names(formals(SumstatsCleaningParam)),
        names(pecotmr:::.sumstatsCleaningResolve(SumstatsCleaningParam()))
    )
    # And every one of them is an .applySanityChecks argument, since the
    # resolved record is spliced straight into it.
    expect_true(all(
        names(pecotmr:::.sumstatsCleaningResolve(SumstatsCleaningParam())) %in%
            names(formals(pecotmr:::.applySanityChecks))
    ))
})

test_that(".sumstatsCleaningResolve fills defaults an empty bundle omits", {
    # .applySanityChecks tests these with `if (!flag)`, so a missing field
    # must not arrive as NULL.
    res <- pecotmr:::.sumstatsCleaningResolve(list())
    expect_equal(
        res,
        pecotmr:::.sumstatsCleaningResolve(SumstatsCleaningParam())
    )
    expect_true(res$coerceNumeric)
    expect_equal(res$smallPFloor, 5e-324)
    expect_error(
        pecotmr:::.sumstatsCleaningResolve(list(coerceNumeric = FALSE)),
        "must be built with SumstatsCleaningParam"
    )
})


test_that("summaryStatsQc gives one result from an LdMatrix and its panel", {
    data(gwasSumStatsS4Example)
    gss <- gwasSumStatsS4Example
    panel <- ldReference(gss)[[1L]]
    ld <- LdMatrix(
        unname(computeLd(panel)),
        variants = SummarizedExperiment::rowRanges(panel),
        nSamples = ncol(panel)
    )
    withLd <- GwasSumStats(
        studyName = studyName(gss),
        entry = list(.collectionEntry(gss, 1L)),
        genome = "hg19",
        ldReference = ld,
        qcOptions = qcOptions(gss)
    )
    viaPanel <- suppressMessages(
        summaryStatsQc(gss, ldMismatchQcMethod = "slalom")
    )
    viaMatrix <- suppressMessages(
        summaryStatsQc(withLd, ldMismatchQcMethod = "slalom")
    )
    expect_equal(
        as.data.frame(unlist(viaMatrix)),
        as.data.frame(unlist(viaPanel))
    )
    # The reference is trimmed to the surviving variants, still a matrix.
    expect_s4_class(ldReference(viaMatrix)[[1L]], "LdMatrix")
    expect_equal(
        length(ldReference(viaMatrix)[[1L]]),
        length(unlist(viaMatrix))
    )
})

test_that("DENTIST needs an LdMatrix's size; SLALOM does not", {
    ld <- LdMatrix(
        diag(2),
        variants = data.frame(chrom = "1", pos = 1:2, A1 = "A", A2 = "G")
    )
    expect_error(
        .qcMismatchNSample(ld, "dentist"),
        "DENTIST needs the LD reference's sample size"
    )
    expect_null(.qcMismatchNSample(ld, "slalom"))
})

test_that("RAISS imputation in QC runs on an LdMatrix reference", {
    data(gwasSumStatsS4Example)
    gss <- gwasSumStatsS4Example
    panel <- ldReference(gss)[[1L]]
    ld <- LdMatrix(
        unname(computeLd(panel)),
        variants = SummarizedExperiment::rowRanges(panel),
        nSamples = ncol(panel)
    )
    entry <- .collectionEntry(gss, 1L)
    thin <- GwasSumStats(
        studyName = studyName(gss),
        entry = list(entry[seq(1L, length(entry), by = 2L)]),
        genome = "hg19",
        ldReference = ld
    )
    res <- suppressMessages(summaryStatsQc(thin, impute = TRUE))
    ea <- res$qcAudit[[1L]]
    expect_gt(ea$raissImputedVariants, 0L)
    # Target filters read the LdMatrix's AF, which this one does not record.
    expect_error(
        suppressMessages(summaryStatsQc(
            thin,
            impute = TRUE,
            raissParam = RaissParam(mafCutoff = 0.05)
        )),
        "needs allele frequencies"
    )
})

test_that("summaryStatsQc QCs each row against its own LD reference", {
    # One row on the panel, another on the same LD as a matrix: each is
    # harmonized against its own and keeps it, cut to its own variants.
    data(qtlSumStatsMulticontextExample)
    mc <- qtlSumStatsMulticontextExample
    panel <- ldReference(mc)[[1L]]
    ldm <- LdMatrix(
        unname(computeLd(panel)),
        variants = SummarizedExperiment::rowRanges(panel)
    )
    refs <- rep(list(panel), nrow(mc))
    refs[[2L]] <- ldm
    out <- suppressMessages(summaryStatsQc(pecotmr:::.withRowLdRefs(mc, refs)))
    expect_length(ldReference(out), 2L)
    expect_s4_class(
        pecotmr:::.rowLdRef(out, 1L),
        "RangedSummarizedExperiment"
    )
    expect_s4_class(pecotmr:::.rowLdRef(out, 2L), "LdMatrix")
    expect_equal(length(pecotmr:::.rowLdRef(out, 2L)), length(out[[2L]]))
})

# ===========================================================================
# RAISS's flank stops at a row's own LD block
# ===========================================================================

# Every fourth variant of the bundled example, so RAISS has targets, split
# into two LD blocks of one study. `blockId` names the blocks, or leaves
# their coordinates as their keys.
.qcb_thinBlocked <- function(blockId = NULL) {
    data(gwasSumStatsS4Example)
    gss <- gwasSumStatsS4Example
    variants <- unlist(gss)
    blocks <- GenomicRanges::GRanges(
        "chr22",
        IRanges::IRanges(c(1L, 15200001L), c(15200000L, 51304566L))
    )
    if (!is.null(blockId)) {
        S4Vectors::mcols(blocks)$blockId <- blockId
    }
    GwasSumStats(
        studyName = studyName(gss),
        entry = list(variants[seq(1L, length(variants), by = 4L)]),
        genome = unique(unname(GenomeInfoDb::genome(gss))),
        ldReference = ldReference(gss),
        qcOptions = list(step1 = "ok"),
        region = blocks
    )
}

test_that(".ssqcRaissBounds keeps each row of a split study to its block", {
    x <- .qcb_thinBlocked()
    on <- list(impute = TRUE, raissParam = list(flank = 1e5))
    bounds <- pecotmr:::.ssqcRaissBounds(x, on)
    # The blocks' ranges were given, so the blocks are their own bounds,
    # whether the manifest keys them by their coordinates or by name.
    expect_identical(
        map_chr(bounds, as.character),
        c("chr22:1-15200000", "chr22:15200001-51304566")
    )
    named <- pecotmr:::.ssqcRaissBounds(.qcb_thinBlocked(c("a", "b")), on)
    expect_identical(
        map_chr(named, as.character),
        c("chr22:1-15200000", "chr22:15200001-51304566")
    )
    # Rows whose ranges are only their spans fall back to halfway between.
    computed <- x
    computed$regionGiven <- c(FALSE, FALSE)
    labelled <- pecotmr:::.ssqcRaissBounds(computed, on)
    spans <- range(x)
    gap <- (GenomicRanges::end(spans[[1L]]) +
        GenomicRanges::start(spans[[2L]])) %/%
        2L
    expect_identical(GenomicRanges::end(labelled[[1L]]), gap)
    expect_identical(GenomicRanges::start(labelled[[2L]]), gap + 1L)
    # No flank, or no imputation, needs no bound; nor does a lone row.
    expect_true(all(map_lgl(
        pecotmr:::.ssqcRaissBounds(
            x,
            list(impute = TRUE, raissParam = list(flank = 0))
        ),
        is.null
    )))
    expect_true(all(map_lgl(
        pecotmr:::.ssqcRaissBounds(x, list(impute = FALSE)),
        is.null
    )))
    data(gwasSumStatsS4Example)
    expect_null(
        pecotmr:::.ssqcRaissBounds(gwasSumStatsS4Example, on)[[1L]]
    )
    # A lone row given its block keeps to it all the same: what lies past
    # its range is another block, loaded or not.
    lone <- x[1L]
    expect_identical(
        as.character(pecotmr:::.ssqcRaissBounds(lone, on)[[1L]]),
        "chr22:1-15200000"
    )
})

test_that(".ssqcMidpointBounds splits gaps; an overlap keeps its own span", {
    own <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100L, 200L))
    sib <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(10L, 300L), c(50L, 400L))
    )
    mid <- pecotmr:::.ssqcMidpointBounds(own, sib)
    expect_identical(GenomicRanges::start(mid), 76L)
    expect_identical(GenomicRanges::end(mid), 250L)
    # Neighbours on one side only leave the other side open.
    below <- pecotmr:::.ssqcMidpointBounds(own, sib[1L])
    expect_identical(GenomicRanges::end(below), .Machine$integer.max)
    above <- pecotmr:::.ssqcMidpointBounds(own, sib[2L])
    expect_identical(GenomicRanges::start(above), 1L)
    over <- GenomicRanges::GRanges("chr1", IRanges::IRanges(150L, 400L))
    expect_identical(pecotmr:::.ssqcMidpointBounds(own, over), own)
})

test_that(".qcRaissClampBounds narrows only the clamp's chromosome", {
    bounds <- tibble::tibble(
        chrom = c("chr1", "chr2"),
        lo = c(50L, 50L),
        hi = c(500L, 500L)
    )
    clamp <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100L, 300L))
    out <- pecotmr:::.qcRaissClampBounds(bounds, clamp)
    expect_identical(out$lo, c(100L, 50L))
    expect_identical(out$hi, c(300L, 500L))
    expect_identical(pecotmr:::.qcRaissClampBounds(bounds, NULL), bounds)
})

test_that("RAISS with a flank imputes each variant of a split study once", {
    x <- .qcb_thinBlocked()
    out <- suppressWarnings(suppressMessages(summaryStatsQc(
        x,
        impute = TRUE,
        raissParam = RaissParam(flank = 2e5, r2Threshold = 0, minimumLd = 0)
    )))
    expect_gt(sum(lengths(out)), sum(lengths(x)))
    ids <- map(as.list(out), function(el) S4Vectors::mcols(el)$SNP)
    expect_length(intersect(ids[[1L]], ids[[2L]]), 0L)
    # Each row stays inside its own block.
    expect_lte(max(GenomicRanges::end(out[[1L]])), 15200000L)
    expect_gte(min(GenomicRanges::start(out[[2L]])), 15200001L)
})

test_that("a signal screen may empty every block of a study", {
    x <- .qcb_thinBlocked()
    out <- suppressMessages(summaryStatsQc(
        x,
        signalScreenParam = SignalScreenParam(absZ = 1e6)
    ))
    expect_equal(unname(lengths(out)), c(0L, 0L))
    expect_identical(out$blockId, x$blockId)
    # The blocks were given their ranges, which they keep with no variant
    # left in them.
    expect_identical(out$regionGiven, c(TRUE, TRUE))
    expect_identical(out$region, x$region)
})

test_that("a QC pass carries each row's window unchanged", {
    x <- .qcb_thinBlocked(c("a", "b"))
    out <- suppressMessages(summaryStatsQc(x))
    expect_identical(out$blockId, c("a", "b"))
    expect_identical(out$regionGiven, c(TRUE, TRUE))
    expect_identical(out$region, x$region)
    # A computed window is the span the row was first built with: QC dropping
    # an edge variant does not move it, and range(x) gives what is left.
    data(gwasSumStatsS4Example)
    v <- sumStats(gwasSumStatsS4Example)
    S4Vectors::mcols(v)$Z[[1L]] <- NA
    built <- GwasSumStats(
        studyName = "s",
        entry = list(v),
        genome = "hg19",
        ldReference = ldReference(gwasSumStatsS4Example)[[1L]],
        qcOptions = list(step1 = "ok")
    )
    qc <- suppressWarnings(suppressMessages(summaryStatsQc(built)))
    expect_false(qc$regionGiven)
    expect_identical(qc$region, built$region)
    expect_gt(
        GenomicRanges::start(range(qc[[1L]])),
        GenomicRanges::start(qc$region)
    )
})

test_that("a QC pass keeps a QtlSumStats's trait position and columns", {
    data(qtlSumStatsExample)
    q <- qtlSumStatsExample
    g <- q[[1L]]
    tp <- GenomicRanges::GRanges(
        as.character(GenomicRanges::seqnames(g))[[1L]],
        IRanges::IRanges(GenomicRanges::start(g)[[50L]], width = 1L)
    )
    y <- QtlSumStats(
        studyName = q$study,
        context = q$context,
        trait = q$trait,
        entry = list(g),
        genome = unique(GenomeInfoDb::genome(q)),
        ldReference = ldReference(q)[[1L]],
        traitPos = tp,
        extraCols = list(tag = "kept")
    )
    out <- suppressMessages(summaryStatsQc(y))
    expect_identical(as.character(out$traitPos), as.character(tp))
    expect_identical(out$tag, "kept")
    expect_identical(out$region, y$region)
})

test_that(".ssqcRaissBounds leaves an emptied row unbounded", {
    x <- .qcb_thinBlocked()
    emptied <- GwasSumStats(
        studyName = "trait1",
        entry = list(sumStats(x)[0L]),
        genome = "hg19",
        ldReference = ldReference(x)[[1L]],
        qcOptions = list(step1 = "ok"),
        blockId = "screened"
    )
    both <- combineGwasSumStats(x, emptied)
    expect_equal(unname(lengths(both)), c(unname(lengths(x)), 0L))
    bounds <- pecotmr:::.ssqcRaissBounds(
        both,
        list(impute = TRUE, raissParam = list(flank = 1e5))
    )
    # An empty row has no chromosome to share, so it is never bounded; the
    # two blocks still bound each other.
    expect_null(bounds[[3L]])
    expect_false(is.null(bounds[[1L]]))
    expect_false(is.null(bounds[[2L]]))
})

# ===========================================================================
# On a real genotype panel: dentistSingleWindow
# ===========================================================================

test_that("dentistSingleWindow detects injected outliers", {
    skip_if_not_installed("pgenlibr")
    # 100 samples x 349 variants on chr21
    X <- loadGenotypeRegion(test_path("test_data", "test_variants"))
    R <- computeLd(X, method = "sample")
    set.seed(42)
    z <- rnorm(ncol(X))
    # Inject extreme outliers
    z[1] <- 50
    z[2] <- -50
    result <- suppressWarnings(dentistSingleWindow(z, R = R, nSample = nrow(X)))
    # At least one of the injected values should be flagged
    expect_true(any(result$outlier))
})

# ===========================================================================
# .entryToSumstatDf
# ===========================================================================

test_that(".entryToSumstatDf: returns the canonical column layout", {
    ss <- .gsd_makeGwasSumStats(data.frame(
        SNP = c("rs1", "rs2"),
        A1 = c("A", "A"),
        A2 = c("G", "G"),
        Z = c(1.0, -2.5),
        N = c(1000L, 1500L),
        MAF = c(0.1, 0.3)
    ))
    df <- pecotmr:::.entryToSumstatDf(sumStats(ss))
    expect_s3_class(df, "data.frame")
    expect_equal(df$variant_id, c("rs1", "rs2"))
    expect_equal(df$chrom, c("chr1", "chr1"))
    expect_equal(df$pos, c(100L, 200L))
    expect_equal(df$z, c(1.0, -2.5))
    expect_equal(df$N, c(1000, 1500))
    expect_equal(df$maf, c(0.1, 0.3))
})

test_that(".entryToSumstatDf: derives z from beta/se on request", {
    ss <- .gsd_makeGwasSumStats(
        data.frame(
            SNP = "rs1",
            A1 = "A",
            A2 = "G",
            BETA = 0.5,
            SE = 0.1,
            N = 1000L
        ),
        snp_n = 1L
    )
    df <- pecotmr:::.entryToSumstatDf(sumStats(ss), derive = "zFromBetaSe")
    expect_equal(df$beta, 0.5)
    expect_equal(df$se, 0.1)
    expect_equal(df$z, 0.5 / 0.1)
})

test_that(".entryToSumstatDf: omits optional columns when absent", {
    ss <- .gsd_makeGwasSumStats(
        data.frame(
            SNP = "rs1",
            A1 = "A",
            A2 = "G"
        ),
        snp_n = 1L
    )
    df <- pecotmr:::.entryToSumstatDf(sumStats(ss))
    expect_false("z" %in% names(df))
    expect_false("beta" %in% names(df))
    expect_false("se" %in% names(df))
    expect_false("N" %in% names(df))
    expect_false("maf" %in% names(df))
})

test_that(".entryToSumstatDf: require = c('Z', 'N') errors when missing", {
    ss <- .gsd_makeGwasSumStats(
        data.frame(
            SNP = "rs1",
            A1 = "A",
            A2 = "G",
            Z = 1.5
        ),
        snp_n = 1L
    )
    expect_error(
        pecotmr:::.entryToSumstatDf(sumStats(ss), require = c("Z", "N")),
        "no N mcol"
    )
})
