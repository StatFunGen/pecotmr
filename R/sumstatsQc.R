# =============================================================================
# Summary-statistic QC pipeline
# -----------------------------------------------------------------------------
# Consolidated QC suite for GwasSumStats / QtlSumStats objects. The
# top-level entry point is `summaryStatsQc()` (below), which orchestrates
# the individual passes:
#
#   * Allele harmonization (harmonizeAlleles, in R/variantId.R)
#   * RAISS sumstats imputation (fills variants present on the LD panel
#     but missing from sumstats)
#   * SLALoM (single-causal-variant ABF outlier detection)
#   * DENTIST (test-based LD-mismatch detection)
#   * Univariate RSS diagnostics (post-finemap mismatch diagnostics)
#
# Each pass occupies its own section below; the orchestrator lives at the
# bottom. Pure individual-level sample QC (relatedness etc.) is in
# R/relatednessQc.R, not here.
# =============================================================================

#' @importFrom GenomicRanges seqnames
#' @importFrom S4Vectors mcols
NULL

# Allele harmonization (harmonizeAlleles) now lives in R/variantId.R, alongside
# the other variant-matching primitives (parseVariantId, matchVariants). It is
# package-internal and called from here (.matchAgainstSketch), ctwasPipeline,
# and the pipeline join sites.

# Standardize a variant set (GRanges or data.frame) to a chrom/pos/alt/ref
# frame.
# @noRd
.variantsToDf <- function(x) {
    df <- if (is(x, "GRanges")) {
        mc <- as.data.frame(mcols(x))
        mc$chrom <- as.character(seqnames(x))
        mc$pos <- start(x)
        as_tibble(mc)
    } else {
        as_tibble(x)
    }
    select(df, all_of(c("chrom", "pos", "alt", "ref")))
}

# Canonical per-variant key from sorted alleles, so strand/allele flips collide.
# @noRd
.canonicalAlleleKey <- function(df) {
    aMin <- pmin(df$alt, df$ref)
    aMax <- pmax(df$alt, df$ref)
    str_c(df$chrom, df$pos, aMin, aMax, sep = " ")
}

#' Merge variant info from two sources with allele-flip-aware matching
#'
#' Merges variant metadata (chromosome, position, ref, alt) from two sources,
#' detecting and correcting allele flips (where alt/ref are swapped). Creates a
#' canonical key from sorted alleles to match across datasets.
#'
#' @param variants1 A data.frame with columns \code{chrom}, \code{pos},
#'   \code{alt}, \code{ref}, or a \code{GRanges} with corresponding metadata
#'   columns.
#' @param variants2 A data.frame or \code{GRanges} with the same columns.
#' @param all Logical. If TRUE (default), returns the union of both sets. If
#'   FALSE, returns only variants from \code{variants2} (flipped to match
#'   \code{variants1}'s allele orientation).
#' @return A data.frame with columns \code{chrom}, \code{pos}, \code{alt},
#'   \code{ref}, deduplicated by position and alleles.
#' @examples
#' v1 <- data.frame(chrom = "1", pos = 1:3, alt = "A", ref = "G")
#' v2 <- data.frame(chrom = "1", pos = 2:4, alt = "A", ref = "G")
#' mergeVariantInfo(v1, v2, all = TRUE)
#' @export
mergeVariantInfo <- function(variants1, variants2, all = TRUE) {
    df1 <- .variantsToDf(variants1)
    df2 <- .variantsToDf(variants2)

    key1 <- .canonicalAlleleKey(df1)
    key2 <- .canonicalAlleleKey(df2)

    # Detect flips: where df2's alt matches df1's ref at the same key
    matchIdx <- match(key2, key1)
    hasMatch <- !is.na(matchIdx)

    flip <- rep(FALSE, nrow(df2))
    mi <- matchIdx[hasMatch]
    flip[hasMatch] <- df2$alt[hasMatch] == df1$ref[mi] &
        df2$ref[hasMatch] == df1$alt[mi]

    # Apply flips to df2
    flipRows <- which(hasMatch)[flip[hasMatch]]
    if (length(flipRows) > 0) {
        tmp <- df2$alt[flipRows]
        df2$alt[flipRows] <- df2$ref[flipRows]
        df2$ref[flipRows] <- tmp
    }

    if (all) {
        distinct(bind_rows(df1, df2))
    } else {
        df2
    }
}


# =============================================================================
# DENTIST: deterministic test-based LD-mismatch detection
# =============================================================================

#' Resolve LD Input: Accept Either R (LD matrix) or X (Genotype Matrix)
#'
#' Internal helper that validates and resolves the LD input for QC functions.
#' Exactly one of \code{R} or \code{X} must be provided. When \code{X} is
#' provided, LD is computed via \code{computeLd(X)} and \code{nSample} defaults
#' to \code{nrow(X)}.
#'
#' @param R Square LD correlation matrix, or NULL.
#' @param X Genotype matrix (samples x SNPs), or NULL.
#' @param nSample Sample size. Required when \code{R} is provided and
#'   \code{needNSample} is TRUE; inferred from \code{X} when \code{X} is
#'   provided.
#' @param needNSample Logical; if TRUE, \code{nSample} must be available (either
#'   provided or inferred from \code{X}).
#'
#' @return A list with components \code{R} (LD correlation matrix) and
#'   \code{nSample} (integer or NULL).
#'
#' @noRd
resolveLdInput <- function(
    R = NULL,
    X = NULL,
    nSample = NULL,
    needNSample = FALSE,
    ldMethod = "sample"
) {
    if (is.null(R) && is.null(X)) {
        abort("Either R (LD matrix) or X (genotype matrix) must be provided.")
    }
    if (!is.null(R) && !is.null(X)) {
        abort("Provide either R or X, not both.")
    }
    if (!is.null(X)) {
        if (!is.matrix(X)) {
            X <- as.matrix(X)
        }
        if (is.null(nSample)) {
            nSample <- nrow(X)
        }
        R <- computeLd(X, method = ldMethod)
    }
    if (needNSample && is.null(nSample)) {
        abort("nSample is required when providing an LD matrix R.")
    }
    list(R = R, nSample = nSample)
}

# --- dentist helpers --------------------------------------------------------

# Validate/rename the pos and z columns (accepting position/zscore) + sort.
.dentistResolveColumns <- function(sumStat) {
    lc <- str_to_lower(colnames(sumStat))
    if (
        !any(is_in(c("pos", "position"), lc)) ||
            !any(is_in(c("z", "zscore"), lc))
    ) {
        msg <- glue(
            "Input sumStat is missing either 'pos'/'position' or ",
            "'z'/'zscore' column."
        )
        abort(msg)
    }
    if (!is_in("pos", lc)) {
        colnames(sumStat)[which(is_in(lc, "position"))] <- "pos"
    }
    if (!is_in("z", lc)) {
        colnames(sumStat)[which(is_in(lc, "zscore"))] <- "z"
    }
    arrange(sumStat, .data$pos)
}

# Run DENTIST on a single window, unpacking the shared tuning parameters.
.dentistCallSingle <- function(zScore, ldMat, nSample, p) {
    dentistSingleWindow(
        zScore,
        R = ldMat,
        nSample = nSample,
        pValueThreshold = p$pValueThreshold,
        propSVD = p$propSVD,
        gcControl = p$gcControl,
        nIter = p$nIter,
        gPvalueThreshold = p$gPvalueThreshold,
        duprThreshold = p$duprThreshold,
        ncpus = p$ncpus,
        correctChenEtAlBug = p$correctChenEtAlBug,
        seed = p$seed
    )
}

# Segment into windows, run DENTIST per window, and merge the results.
.dentistWindows <- function(
    sumStat,
    ldMat,
    nSample,
    windowMode,
    windowSize,
    minDim,
    p
) {
    if (windowMode == "distance") {
        windowDividedRes <- segmentByDist(
            sumStat$pos,
            maxDist = windowSize,
            minDim = minDim
        )
    } else {
        windowDividedRes <- segmentByCount(sumStat$pos, maxCount = minDim)
    }
    dentistResultByWindow <- list()
    for (k in seq_len(nrow(windowDividedRes))) {
        # windowEndIdx is 1-based exclusive; convert to an inclusive range.
        idxRange <- windowDividedRes$windowStartIdx[
            k
        ]:(windowDividedRes$windowEndIdx[k] - 1L)
        zScoreK <- sumStat$z[idxRange]
        ldMatK <- ldMat[idxRange, idxRange]
        dentistResultByWindow[[k]] <- .dentistCallSingle(
            zScoreK,
            ldMatK,
            nSample,
            p
        )
    }
    mergeWindows(dentistResultByWindow, windowDividedRes)
}

#' Detect Outliers Using Dentist Algorithm
#'
#' DENTIST (Detecting Errors iN analyses of summary staTISTics) is a quality
#' control tool for GWAS summary data. It uses linkage disequilibrium (LD)
#' information from a reference panel to identify and correct problematic
#' variants by comparing observed GWAS statistics to predicted values. It can
#' detect errors in genotyping/imputation, allelic errors, and heterogeneity
#' between GWAS and LD reference samples.
#'
#' @param sumStat A data frame containing summary statistics, including 'pos' or
#'   'position' and 'z' or 'zscore' columns.
#' @param R Square LD correlation matrix. Provide either \code{R} or \code{X}.
#' @param X Genotype matrix (samples x SNPs). If provided, LD is computed via
#'   \code{computeLd(X)} and \code{nSample} defaults to \code{nrow(X)}.
#' @param nSample The number of samples in the LD reference panel (NOT the GWAS
#'   sample size). This controls the SVD truncation rank K = min(idx_size,
#'   nSample) * propSVD. Required when \code{R} is provided; inferred from
#'   \code{X} when \code{X} is provided.
#' @param windowSize The size of the window for dividing the genomic region in
#'   distance mode (base pairs). Default is 2000000 (2 Mb). Only used when
#'   \code{windowMode = "distance"}.
#' @param windowMode Character string specifying the windowing strategy:
#'   \code{"distance"} (default) creates windows by physical distance using
#'   \code{segmentByDist} (C++ \code{--wind-dist}), and \code{"count"} creates
#'   windows by variant count using \code{segmentByCount} (C++ \code{--wind}).
#' @param pValueThreshold The p-value threshold for significance. Default is
#'   5e-8.
#' @param propSVD The proportion of singular value decomposition (SVD) to use.
#'   Default is 0.4.
#' @param gcControl Logical indicating whether genomic control should be
#'   applied. Default is FALSE.
#' @param nIter The number of iterations for the Dentist algorithm. Default is
#'   10.
#' @param gPvalueThreshold The genomic p-value threshold for significance.
#'   Default is 0.05.
#' @param duprThreshold The absolute correlation r value threshold to be
#'   considered duplicate. Default is 0.99.
#' @param ncpus The number of CPU cores to use for parallel processing. Default
#'   is 1.
#' @param correctChenEtAlBug Logical indicating whether to correct the Chen et
#'   al. bug. Default is TRUE.
#' @param minDim In distance mode: minimum number of SNPs per block (default
#'   2000). In count mode: the number of variants per window (i.e., the window
#'   size).
#' @param ldMethod Character string specifying the LD computation method when
#'   \code{X} is provided. Passed to \code{computeLd}. One of \code{"sample"}
#'   (default), \code{"population"}, or \code{"gcta"}. Ignored when \code{R} is
#'   provided directly.
#' @param seed Integer or \code{NULL}. Random seed for the iterative
#'   variant-partitioning RNG. \code{NULL} (default) preserves the original
#'   DENTIST hard-coded seeds (\code{10} for the initial partition,
#'   \code{20000 + t * 20000} per iteration) for exact fidelity with the
#'   reference binary; a provided seed makes the partitioning reproducible under
#'   a value of your choosing.
#'
#' @return A data frame containing the imputed result and detected outliers.
#'
#' The returned data frame includes the following columns:
#'
#' \describe{
#'   \item{\code{original_z}}{The original z-score values from the input
#'   \code{sumStat}.}
#'   \item{\code{imputed_z}}{The imputed z-score values computed by the Dentist
#'   algorithm.}
#'   \item{\code{rsq}}{The coefficient of determination (R-squared) between
#'   original and imputed z-scores.}
#'   \item{\code{iter_to_correct}}{The number of iterations required to correct
#'   the z-scores, if applicable.}
#'   \item{\code{index_within_window}}{The index of the observation within the
#'   window.}
#'   \item{\code{index_global}}{The global index of the observation.}
#'   \item{\code{outlier_stat}}{The computed statistical value based on the
#'   original and imputed z-scores and R-squared.}
#'   \item{\code{outlier}}{A logical indicator specifying whether the
#'   observation is identified as an outlier based on the statistical test.}
#' }
#'
#' @examples
#' # Simulate summary statistics for 100 variants and an LD matrix estimated
#' # from a 500-sample reference panel, then screen for LD-outlier variants.
#' set.seed(1)
#' nSample <- 500
#' geno <- matrix(rbinom(nSample * 100, 2, 0.3), nrow = nSample, ncol = 100)
#' ldMat <- cor(geno)
#' sumStat <- data.frame(pos = seq_len(100), z = rnorm(100))
#' result <- dentist(sumStat, R = ldMat, nSample = nSample)
#' head(result)
#'
#' @details
#' Windowing supports two modes matching the original DENTIST C++ binary:
#' \itemize{
#'   \item \code{"distance"} (default): Uses the \code{segmentingByDist}
#'   algorithm
#'     (C++ \code{--wind-dist}), implemented in \code{segmentByDist}.
#'     Windows span a fixed physical distance (\code{windowSize} bp).
#'   \item \code{"count"}: Uses the \code{segmentedQCed} algorithm
#'     (C++ \code{--wind}), implemented in \code{segmentByCount}.
#'     Windows contain a fixed number of variants (\code{minDim}).
#'     Useful when regions have sparse variants where distance-based windows
#'     would create windows with too few variants.
#' }
#' The \code{correctChenEtAlBug} parameter affects the iterative filtering
#' in two ways:
#' \enumerate{
#'   \item Comparison between iteration index \code{t} and \code{nIter}
#'   (explained in source code)
#'   \item The \code{!grouping_tmp} operator bug (explained in source code)
#' }
#'
#' @export
dentist <- function(
    sumStat,
    R = NULL,
    X = NULL,
    nSample = NULL,
    windowSize = 2000000,
    windowMode = c("distance", "count"),
    pValueThreshold = 5.0369e-8,
    propSVD = 0.4,
    gcControl = FALSE,
    nIter = 10,
    gPvalueThreshold = 0.05,
    duprThreshold = 0.99,
    ncpus = 1,
    correctChenEtAlBug = TRUE,
    minDim = 2000,
    ldMethod = "sample",
    seed = NULL
) {
    resolved <- resolveLdInput(
        R = R,
        X = X,
        nSample = nSample,
        needNSample = TRUE,
        ldMethod = ldMethod
    )
    ldMat <- resolved$R
    nSample <- resolved$nSample
    sumStat <- .dentistResolveColumns(sumStat)
    windowMode <- arg_match(windowMode)
    p <- as.list(environment())
    if (nrow(sumStat) < minDim) {
        return(.dentistCallSingle(sumStat$z, ldMat, nSample, p))
    }
    .dentistWindows(
        sumStat,
        ldMat,
        nSample,
        windowMode,
        windowSize,
        minDim,
        p
    )
}

# --- dentistSingleWindow helpers -------------------------------------------

# Warn on small windows; validate the LD matrix shape against zScore.
.dentistValidateInput <- function(zScore, ldMat) {
    if (length(zScore) < 2000) {
        nZ <- length(zScore)
        msg <- glue(
            "The number of variants ({nZ}) is below 2000. The algorithm ",
            "may not work as expected, as suggested by the original ",
            "DENTIST. Consider using windowMode = 'count' with an ",
            "appropriate minDim to control window sizes by variant ",
            "count."
        )
        warn(msg)
    }
    if (
        !is.matrix(ldMat) ||
            nrow(ldMat) != ncol(ldMat) ||
            nrow(ldMat) != length(zScore)
    ) {
        msg <- glue(
            "ldMat must be a square matrix with dimensions equal to ",
            "the length of zScore."
        )
        abort(msg)
    }
}

# Optionally deduplicate near-perfectly-correlated variants before imputation.
.dentistDedup <- function(zScore, ldMat, duprThreshold) {
    dedupRes <- NULL
    rThreshold <- round(sqrt(duprThreshold) * 1000) / 1000
    if (duprThreshold < 1.0) {
        dedupRes <- .findDuplicateVariants(zScore, ldMat, rThreshold)
        numDup <- sum(dedupRes$dupBearer != -1)
        if (numDup > 0) {
            nZ <- length(zScore)
            msg <- glue(
                "{numDup} duplicated variants out of a total of {nZ} ",
                "were found at r threshold of {rThreshold}"
            )
            inform(msg)
        }
        zScore <- dedupRes$filteredZ
        ldMat <- dedupRes$filteredLD
    }
    list(zScore = zScore, ldMat = ldMat, dedupRes = dedupRes)
}

# Run the C++ iterative imputation and snake_case the output. The C++ returns
# any rsq values it capped at 1.0 in `rsqExceed`; summarize them into a single
# warning here (no warning handler / shared env needed).
.dentistRunImpute <- function(ldMat, nSample, zScore, p) {
    verboseIter <- getOption("pecotmr.dentist.verbose", FALSE)
    res <- dentistIterativeImpute(
        # cpp11 requires exact integer types for int parameters
        ldMat,
        as.integer(nSample),
        zScore,
        p$pValueThreshold,
        p$propSVD,
        p$gcControl,
        as.integer(p$nIter),
        p$gPvalueThreshold,
        as.integer(p$ncpus),
        p$correctChenEtAlBug,
        verboseIter,
        if (is.null(p$seed)) NULL else as.integer(p$seed)
    )
    rsqExceed <- res$rsqExceed
    res$rsqExceed <- NULL
    if (length(rsqExceed) > 0) {
        nExceed <- length(rsqExceed)
        maxExceed <- max(rsqExceed)
        msg <- glue(
            "{nExceed} rsq values exceeded 1 (capped at 1.0). ",
            "Max reported: {maxExceed}"
        )
        warn(msg)
    }
    # cpp11 wrapper returns camelCase keys; convert to snake_case columns
    as_tibble(res) |>
        rename(
            original_z = "originalZ",
            imputed_z = "imputedZ",
            z_diff = "zDiff",
            iter_to_correct = "iterToCorrect"
        )
}

# Outlier statistic: (z - imputed)^2 / (1 - rsq), thresholded on the p-value.
.dentistOutlierStat <- function(res, pValueThreshold) {
    res |>
        mutate(
            outlier_stat = (.data$original_z - .data$imputed_z)^2 /
                pmax(1 - .data$rsq, 1e-8),
            outlier = -log10(pchisq(
                .data$outlier_stat,
                df = 1,
                lower.tail = FALSE
            )) >
                -log10(pValueThreshold)
        ) |>
        select(-any_of("z_diff"))
}

#' Perform DENTIST on a single window
#'
#' Detect outliers in GWAS summary statistics using LD-based iterative
#' imputation. Provide either an LD correlation matrix \code{R} or a genotype
#' matrix \code{X} (from which LD and sample size are derived automatically).
#'
#' @param zScore Numeric vector of z-scores.
#' @param R Square LD correlation matrix. Provide either \code{R} or \code{X}.
#' @param X Genotype matrix (samples x SNPs). If provided, LD is computed via
#'   \code{computeLd(X)} and \code{nSample} defaults to \code{nrow(X)}.
#' @param nSample Number of samples in the LD reference panel (NOT the GWAS
#'   sample size). Controls the SVD truncation rank. Required when \code{R} is
#'   provided; inferred from \code{X} when \code{X} is provided.
#' @param pValueThreshold P-value threshold for outlier detection. Default is
#'   5e-8.
#' @param propSVD SVD truncation proportion. Default is 0.4.
#' @param gcControl Logical; apply genomic control. Default is FALSE.
#' @param nIter Number of iterations. Default is 10.
#' @param gPvalueThreshold Grouping p-value threshold. Default is 0.05.
#' @param duprThreshold Duplicate r-squared threshold. Default is 0.99.
#' @param ncpus Number of CPU cores. Default is 1.
#' @param correctChenEtAlBug Correct the original DENTIST operator! bug. Default
#'   is TRUE.
#' @param ldMethod Character string specifying the LD computation method when
#'   \code{X} is provided. Passed to \code{computeLd}. One of \code{"sample"}
#'   (default), \code{"population"}, or \code{"gcta"}. Ignored when \code{R} is
#'   provided directly.
#' @param seed Integer or \code{NULL}. Random seed for the iterative
#'   variant-partitioning RNG. \code{NULL} (default) preserves the original
#'   DENTIST hard-coded seeds (\code{10} for the initial partition,
#'   \code{20000 + t * 20000} per iteration); a provided seed overrides them.
#'
#' @return Data frame with columns: original_z, imputed_z, iter_to_correct, rsq,
#'   is_duplicate, outlier_stat, outlier.
#'
#' @seealso \code{\link{dentist}}, \code{\link{slalom}}
#' @references \url{https://github.com/Yves-CHEN/DENTIST}
#' @examples
#' data(eqtlRegionExample)
#' R <- cor(eqtlRegionExample$X[, 1:20])
#' dentistSingleWindow(zScore = rnorm(20), R = R, nSample = 415)
#' @export
dentistSingleWindow <- function(
    zScore,
    R = NULL,
    X = NULL,
    nSample = NULL,
    pValueThreshold = 5e-8,
    propSVD = 0.4,
    gcControl = FALSE,
    nIter = 10,
    gPvalueThreshold = 0.05,
    duprThreshold = 0.99,
    ncpus = 1,
    correctChenEtAlBug = TRUE,
    ldMethod = "sample",
    seed = NULL
) {
    ld <- resolveLdInput(
        R = R,
        X = X,
        nSample = nSample,
        needNSample = TRUE,
        ldMethod = ldMethod
    )
    nSample <- ld$nSample
    ldMat <- ld$R
    .dentistValidateInput(zScore, ldMat)
    p <- list(
        pValueThreshold = pValueThreshold,
        propSVD = propSVD,
        gcControl = gcControl,
        nIter = nIter,
        gPvalueThreshold = gPvalueThreshold,
        ncpus = ncpus,
        correctChenEtAlBug = correctChenEtAlBug,
        seed = seed
    )
    orgZscore <- zScore
    dedup <- .dentistDedup(zScore, ldMat, duprThreshold)
    res <- .dentistRunImpute(dedup$ldMat, nSample, dedup$zScore, p)
    if (duprThreshold < 1.0) {
        res <- addDupsBackDentist(orgZscore, res, dedup$dedupRes)
    }
    .dentistOutlierStat(res, pValueThreshold)
}

#' Add duplicates back to DENTIST output
#'
#' This function takes the output from the DENTIST algorithm and adds back the
#' duplicated variants based on the output from the `findDuplicateVariants`
#' function.
#' @param zScore The original zScore
#' @param dentistOutput A data frame containing the output from the DENTIST
#'   algorithm.
#' @param findDupOutput A list containing the output from the
#'   `findDuplicateVariants` function.
#'
#' @return A data frame with duplicated variants added back and an additional
#'   column indicating duplicates.
#'
#' @noRd
# --- addDupsBackDentist helpers ---------------------------------------------

# Validate DENTIST output vs the duplicate-bearer bookkeeping.
.dentistValidateDups <- function(zScore, dentistOutput, dupBearer, nrowsDup) {
    if (nrow(dentistOutput) != sum(dupBearer == -1)) {
        msg <- glue(
            "The number of rows in the input data does not match the ",
            "occurrences of -1 in dupBearer."
        )
        abort(msg)
    }
    if (length(zScore) != nrowsDup) {
        abort("Input zScore and findDupOutput have inconsistent dimension")
    }
}

# Map each variant to its row in the de-duplicated DENTIST output.
.dentistBuildAssignIdx <- function(dupBearer, nrowsDup) {
    count <- 1
    assignIdx <- rep(0, nrowsDup)
    for (i in seq_along(dupBearer)) {
        if (dupBearer[i] == -1) {
            assignIdx[i] <- count
            count <- count + 1
        } else {
            assignIdx[i] <- dupBearer[i]
        }
    }
    assignIdx
}

# Rebuild the full per-variant table, recovering duplicates (sign-flipped).
.dentistFillDups <- function(zScore, dentistOutput, findDupOutput, assignIdx) {
    dupBearer <- findDupOutput$dupBearer
    sign <- findDupOutput$sign
    nrowsDup <- length(dupBearer)
    imputedZ <- dentistOutput$imputed_z
    iterToCorrect <- dentistOutput$iter_to_correct
    rsq <- dentistOutput$rsq
    zDiff <- dentistOutput$z_diff
    updatedData <- tibble(
        original_z = numeric(nrowsDup),
        imputed_z = numeric(nrowsDup),
        iter_to_correct = numeric(nrowsDup),
        rsq = numeric(nrowsDup),
        z_diff = numeric(nrowsDup),
        is_duplicate = logical(nrowsDup)
    )
    for (i in seq_len(nrowsDup)) {
        updatedData$original_z[i] <- zScore[i]
        updatedData$iter_to_correct[i] <- iterToCorrect[assignIdx[i]]
        updatedData$rsq[i] <- rsq[assignIdx[i]]
        if (dupBearer[i] == -1) {
            updatedData$imputed_z[i] <- imputedZ[assignIdx[i]]
            updatedData$z_diff[i] <- zDiff[assignIdx[i]]
            updatedData$is_duplicate[i] <- FALSE
        } else {
            # Duplicate: sign-flip imputed_z, recompute z_diff from its own
            # z-score so z_diff^2 matches the binary stat (DENTIST.h l706).
            updatedData$imputed_z[i] <- imputedZ[assignIdx[i]] * sign[i]
            denom <- sqrt(max(1 - updatedData$rsq[i], 1e-8))
            updatedData$z_diff[i] <-
                (zScore[i] - updatedData$imputed_z[i]) / denom
            updatedData$is_duplicate[i] <- TRUE
        }
    }
    updatedData
}

addDupsBackDentist <- function(zScore, dentistOutput, findDupOutput) {
    dupBearer <- findDupOutput$dupBearer
    nrowsDup <- length(dupBearer)
    .dentistValidateDups(zScore, dentistOutput, dupBearer, nrowsDup)
    assignIdx <- .dentistBuildAssignIdx(dupBearer, nrowsDup)
    .dentistFillDups(zScore, dentistOutput, findDupOutput, assignIdx)
}

# ---- Segmentation helpers ----
# detectGaps(), buildSegmentResult(), and slidingWindowLoop() are shared
# by both segmentByDist() and segmentByCount() to avoid code duplication.
# The core overlapping-window loop lives in slidingWindowLoop(); each mode
# only supplies mode-specific callbacks for fill, step, and block-skip logic.

#' Detect Gaps in Genomic Positions
#'
#' Finds positions where the inter-SNP distance exceeds a threshold, e.g.,
#' centromeric regions. Returns a vector of 1-based block boundaries.
#'
#' @param pos Sorted numeric vector of base pair positions.
#' @param gapThreshold Numeric distance threshold for gap detection.
#' @param verbose Logical; print gap info. Default is FALSE.
#'
#' @return Integer vector of 1-based block boundaries, including \code{1}
#'   (start) and \code{length(pos) + 1} (end sentinel).
#'
#' @noRd
detectGaps <- function(pos, gapThreshold, verbose = FALSE) {
    n <- length(pos)
    diffs <- diff(pos)
    allGaps <- c(1L)
    for (i in seq_along(diffs)) {
        if (diffs[i] > gapThreshold) {
            allGaps <- c(allGaps, i + 1L)
        }
    }
    allGaps <- c(allGaps, n + 1L)

    if (verbose && length(allGaps) - 2 > 0) {
        nGaps <- length(allGaps) - 2
        msg <- glue("No. of gaps found: {nGaps}")
        inform(msg)
        for (i in 2:(length(allGaps) - 1)) {
            gapNo <- i - 1
            startPos <- pos[allGaps[i] - 1]
            endPos <- pos[allGaps[i]]
            msg <- glue("  Gap {gapNo}: {startPos} - {endPos}", .trim = FALSE)
            inform(msg)
        }
    }
    allGaps
}

#' Build Segment Result Data Frame
#'
#' Validates, caps indices, optionally prints verbose info, and returns the
#' standardized segmentation result data frame.
#'
#' @param startList Integer vector of window start indices.
#' @param endList Integer vector of window end indices (exclusive).
#' @param fillStartList Integer vector of fill start indices.
#' @param fillEndList Integer vector of fill end indices (exclusive).
#' @param n Total number of positions.
#' @param verbose Logical; print interval info. Default is FALSE.
#'
#' @return A data frame with columns: windowIdx, windowStartIdx, windowEndIdx,
#'   fillStartIdx, fillEndIdx.
#'
#' @noRd
buildSegmentResult <- function(
    startList,
    endList,
    fillStartList,
    fillEndList,
    n,
    verbose = FALSE
) {
    if (length(startList) == 0) {
        abort("No intervals created by segmentation")
    }

    # Cap end indices at n+1 (one past the last valid 1-based index)
    endList <- pmin(endList, n + 1L)
    fillEndList <- pmin(fillEndList, n + 1L)

    if (verbose) {
        inform("Intervals:")
        for (i in seq_along(startList)) {
            s <- startList[i]
            e <- endList[i]
            fs <- fillStartList[i]
            fe <- fillEndList[i]
            msg <- glue("  {i}: SNPs {s}-{e} (fill {fs}-{fe})", .trim = FALSE)
            inform(msg)
        }
    }

    tibble(
        windowIdx = seq_along(startList),
        windowStartIdx = startList,
        windowEndIdx = endList,
        fillStartIdx = fillStartList,
        fillEndIdx = fillEndList
    )
}

#' Sliding Window Loop for Genomic Segmentation
#'
#' Core overlapping-window loop shared by both distance-based and count-based
#' segmentation strategies. Iterates over contiguous blocks (separated by gaps),
#' creates overlapping windows within each block using mode-specific callbacks,
#' and assembles the result.
#'
#' @param allGaps Integer vector of 1-based block boundaries from
#'   \code{\link{detectGaps}}.
#' @param n Total number of positions.
#' @param ctx Named list bundling the caller's mode-specific segmentation state;
#'   passed as the final argument to every callback below (so they can be
#'   top-level functions rather than closures).
#' @param minBlockFn Function(blockSize, ctx) -> logical; returns TRUE if the
#'   block is large enough to process.
#' @param initEndFn Function(startIdx, blockEnd, ctx) -> integer; computes the
#'   initial window end index for the first window in a block.
#' @param fillFn Function(startIdx, endIdx, notStartInterval, notLastInterval,
#'   ctx) -> list(start, end); computes fill boundaries for each window.
#' @param stepFn Function(startIdx, blockEnd, ctx) -> list(startIdx, endIdx);
#'   advances to the next window.
#' @param adjustLastFn Optional function(startIdx, oldStartIdx, endIdx,
#'   blockEnd, ctx) -> integer; adjusts startIdx when the last interval is
#'   detected. Used by distance mode for small-last-interval correction. Default
#'   is NULL (no adjustment).
#' @param verbose Logical; print interval info. Default is FALSE.
#'
#' @return A data frame with columns: windowIdx, windowStartIdx, windowEndIdx,
#'   fillStartIdx, fillEndIdx.
#'
#' @noRd
# --- slidingWindowLoop helpers ----------------------------------------------

# Generate the window tuples for one block via the mode-specific callbacks in
# `fns` (minBlockFn/initEndFn/fillFn/stepFn/adjustLastFn). Returns the per-block
# start/end and fill-start/fill-end vectors (pre first/last fill correction).
.swlLastCheck <- function(
    startIdx,
    oldStartIdx,
    endIdx,
    blockEnd,
    adjustLastFn,
    ctx
) {
    isLast <- blockEnd <= endIdx
    if (isLast && !is.null(adjustLastFn)) {
        startIdx <- adjustLastFn(startIdx, oldStartIdx, endIdx, blockEnd, ctx)
    }
    list(startIdx = startIdx, notLastInterval = !isLast)
}

# One window's last-interval adjustment + fill boundaries. `fillFn` uses the
# PRE-adjustment startIdx (C++ parity: non-overlapping); the returned startIdx
# is
# the POST-adjustment value recorded for the window.
# @noRd
.swlWindow <- function(
    startIdx,
    oldStartIdx,
    endIdx,
    blockEnd,
    notStartInterval,
    fns,
    ctx
) {
    lc <- .swlLastCheck(
        startIdx,
        oldStartIdx,
        endIdx,
        blockEnd,
        fns$adjustLastFn,
        ctx
    )
    fills <- fns$fillFn(
        startIdx,
        endIdx,
        notStartInterval,
        lc$notLastInterval,
        ctx
    )
    list(
        startIdx = lc$startIdx,
        notLastInterval = lc$notLastInterval,
        fills = fills
    )
}

.swlBlockWindows <- function(blockStart, blockEnd, fns, ctx) {
    startIdx <- blockStart
    endIdx <- fns$initEndFn(startIdx, blockEnd, ctx)
    oldStartIdx <- startIdx
    notStartInterval <- FALSE
    notLastInterval <- TRUE
    times <- 0
    starts <- ends <- fillStarts <- fillEnds <- integer(0)
    repeat {
        times <- times + 1
        if (times > 400) {
            abort("Windowing iteration limit exceeded")
        }
        win <- .swlWindow(
            startIdx,
            oldStartIdx,
            endIdx,
            blockEnd,
            notStartInterval,
            fns,
            ctx
        )
        startIdx <- win$startIdx
        notLastInterval <- win$notLastInterval
        starts <- c(starts, startIdx)
        ends <- c(ends, min(endIdx, blockEnd))
        fillStarts <- c(fillStarts, win$fills$start)
        fillEnds <- c(fillEnds, win$fills$end)
        if (!notLastInterval) {
            break
        }
        oldStartIdx <- startIdx
        stepped <- fns$stepFn(startIdx, blockEnd, ctx)
        startIdx <- stepped$startIdx
        endIdx <- stepped$endIdx
        notStartInterval <- TRUE
    }
    list(
        starts = starts,
        ends = ends,
        fillStarts = fillStarts,
        fillEnds = fillEnds
    )
}

slidingWindowLoop <- function(
    allGaps,
    n,
    ctx,
    minBlockFn,
    initEndFn,
    fillFn,
    stepFn,
    adjustLastFn = NULL,
    verbose = FALSE
) {
    startList <- integer(0)
    endList <- integer(0)
    fillStartList <- integer(0)
    fillEndList <- integer(0)
    fns <- list(
        minBlockFn = minBlockFn,
        initEndFn = initEndFn,
        fillFn = fillFn,
        stepFn = stepFn,
        adjustLastFn = adjustLastFn
    )
    for (k in seq_len(length(allGaps) - 1)) {
        blockStart <- allGaps[k]
        blockEnd <- allGaps[k + 1]
        if (!minBlockFn(blockEnd - blockStart, ctx)) {
            next
        }
        w <- .swlBlockWindows(blockStart, blockEnd, fns, ctx)
        # First window's fill starts at the window start; last window's fill
        # ends at the window end.
        w$fillStarts[1] <- w$starts[1]
        w$fillEnds[length(w$fillEnds)] <- w$ends[length(w$ends)]
        startList <- c(startList, w$starts)
        endList <- c(endList, w$ends)
        fillStartList <- c(fillStartList, w$fillStarts)
        fillEndList <- c(fillEndList, w$fillEnds)
    }
    buildSegmentResult(
        startList,
        endList,
        fillStartList,
        fillEndList,
        n,
        verbose
    )
}

# Apply the quarter-distance index map `quaterIdx` `n` times to `x` (n = 1..4
# gives the 1st..4th quarter boundary from x). Used by segmentByDist.
# @noRd
.nthQuaterIdx <- function(x, n, quaterIdx) {
    for (i in seq_len(n)) {
        x <- quaterIdx[x]
    }
    x
}

#' Segment Genomic Region by Distance (Original DENTIST Algorithm)
#'
#' Implements the same windowing/segmentation algorithm as the original DENTIST
#' C++ binary's \code{segmentingByDist} function. Windows are created using
#' quarter-distance SNP index lookups, with gap detection for centromeres and
#' large gaps.
#'
#' @param pos Integer vector of base pair positions (must be sorted).
#' @param maxDist Maximum distance (bp) between SNPs for windowing. Default is
#'   2000000.
#' @param minDim Minimum number of SNPs per window. Default is 2000.
#' @param verbose Logical; print segmentation info. Default is FALSE.
#'
#' @return A data frame with columns: windowIdx, windowStartIdx, windowEndIdx,
#'   fillStartIdx, fillEndIdx. Start indices are 1-based inclusive; end indices
#'   (windowEndIdx, fillEndIdx) are 1-based exclusive (one past last element),
#'   matching the C++ convention. Use \code{startIdx:(endIdx - 1)} for R
#'   inclusive ranges.
#'
#' @details
#' This is a faithful R translation of the C++ \code{segmentingByDist} function.
#' The algorithm:
#' \enumerate{
#'   \item Precomputes for each SNP: the index of the farthest SNP within
#'   \code{maxDist},
#'         and the index of the SNP at \code{maxDist/4} distance.
#'   \item Detects gaps > \code{maxDist/4} in the position vector (e.g.,
#'   centromeres).
#'   \item Creates overlapping windows that slide by half the distance cutoff,
#'   with fill
#'         regions covering the inner three-quarters of each window.
#'   \item The first window's fill starts at the window start; the last
#'   window's fill
#'         ends at the window end.
#' }
#'
#' @seealso \code{\link{dentistSingleWindow}}, \code{\link{dentist}}
#'
#' @noRd
# --- segmentByDist helpers --------------------------------------------------

# For each SNP, the last SNP index within cutoff/4 distance (clamped to [1, n]).
.segByDistQuaterIdx <- function(pos, n, cutoff) {
    quaterIdx <- integer(n)
    j <- 1
    while (j <= n && pos[j] < cutoff / 4 + as.numeric(pos[1])) {
        j <- j + 1
    }
    quaterIdx[1] <- max(j - 1, 1L)
    for (i in 2:n) {
        j <- quaterIdx[i - 1]
        while (j <= n && pos[j] < cutoff / 4 + as.numeric(pos[i])) {
            j <- j + 1
        }
        quaterIdx[i] <- max(j - 1, 1L)
    }
    quaterIdx <- pmin(quaterIdx, n)
    pmax(quaterIdx, 1L)
}

# Advance the window start by one quarter-step and recompute its end.
.segByDistStep <- function(startIdx, blockEnd, quaterIdx) {
    nextStart <- .nthQuaterIdx(startIdx, 2, quaterIdx)
    list(
        startIdx = nextStart,
        endIdx = min(.nthQuaterIdx(nextStart, 4, quaterIdx) + 1, blockEnd)
    )
}

# Distance-mode sliding-window segmentation (fill = inner 50% by distance).
.segByDistRun <- function(
    allGaps,
    n,
    pos,
    cutoff,
    minDim,
    minBlockSize,
    quaterIdx,
    verbose
) {
    ctx <- list(
        minBlockSize = minBlockSize,
        minDim = minDim,
        quaterIdx = quaterIdx,
        pos = pos,
        cutoff = cutoff,
        n = n
    )
    slidingWindowLoop(
        allGaps,
        n,
        ctx = ctx,
        minBlockFn = .segDistMinBlock,
        initEndFn = .segDistInitEnd,
        fillFn = .segDistFill,
        stepFn = .segDistStep,
        adjustLastFn = .segDistAdjustLast,
        verbose = verbose
    )
}

segmentByDist <- function(
    pos,
    maxDist = 2000000,
    minDim = 2000,
    verbose = FALSE
) {
    n <- length(pos)
    if (n == 0) {
        abort("No positions provided")
    }
    cutoff <- maxDist
    minBlockSize <- minDim
    quaterIdx <- .segByDistQuaterIdx(pos, n, cutoff)
    allGaps <- detectGaps(pos, gapThreshold = cutoff / 4, verbose = verbose)
    .segByDistRun(
        allGaps,
        n,
        pos,
        cutoff,
        minDim,
        minBlockSize,
        quaterIdx,
        verbose
    )
}

#' Segment Genomic Region by Variant Count
#'
#' Implements the windowing algorithm from the original DENTIST C++ binary's
#' \code{segmentedQCed} function. Windows contain a fixed number of variants
#' rather than spanning a fixed physical distance.
#'
#' @param pos Integer vector of base pair positions (must be sorted).
#' @param maxCount Maximum number of variants per window.
#' @param gapDist Physical distance threshold for centromeric gap detection.
#'   Default is 1e6 (matching the C++ hardcoded value).
#' @param verbose Logical; print segmentation info. Default is FALSE.
#'
#' @return A data frame with the same structure as \code{segmentByDist}:
#'   windowIdx, windowStartIdx, windowEndIdx, fillStartIdx, fillEndIdx. End
#'   indices are 1-based exclusive (one past last element).
#'
#' @details
#' This is a faithful R translation of the C++ \code{segmentedQCed} windowing
#' algorithm. Key differences from \code{segmentByDist}:
#' \itemize{
#'   \item Windows are sized by variant count, not physical distance.
#'   \item Uses simple index arithmetic (step = maxCount/2) instead of
#'         distance-based quarter-index lookups.
#'   \item Gap detection uses a fixed 1 Mb threshold (centromeres) instead of
#'         distance/4.
#'   \item Adaptive tail absorption: if fewer than \code{maxCount/2} variants
#'         remain after a window, the window extends to cover the rest.
#' }
#'
#' @seealso \code{segmentByDist}, \code{\link{dentist}}
#'
#' @noRd
segmentByCount <- function(pos, maxCount, gapDist = 1e6, verbose = FALSE) {
    n <- length(pos)
    if (n == 0) {
        abort("No positions provided")
    }

    cutoff <- as.integer(maxCount)
    quarter <- cutoff %/% 4L
    half <- cutoff %/% 2L

    # Detect centromeric gaps (C++ line 784: diff > 1e6)
    allGaps <- detectGaps(pos, gapThreshold = gapDist, verbose = verbose)

    ctx <- list(half = half, cutoff = cutoff, quarter = quarter)
    slidingWindowLoop(
        allGaps,
        n,
        ctx = ctx,
        minBlockFn = .segCountMinBlock,
        initEndFn = .segCountInitEnd,
        fillFn = .segCountFill,
        stepFn = .segCountStep,
        verbose = verbose
    )
}

#' Merge dentist Results by Window
#'
#' This function merges DENTIST results by window into a single data frame.
#'
#' @param dentistResultByWindow A list containing imputed results for each
#'   window.
#' @param windowDividedRes A data frame containing information about the divided
#'   windows.
#'
#' @return A data frame containing merged results.
#'
#' @details The function checks if the number of imputed results matches the
#'   number of windows. It then merges the results by window, adding an index
#'   within the window and a global index. Finally, it extracts the results
#'   within the fillers and combines them into a single data frame.
#'
#' @noRd
mergeWindows <- function(dentistResultByWindow, windowDividedRes) {
    if (length(dentistResultByWindow) != nrow(windowDividedRes)) {
        abort("Different number of windows and imputed results!")
    }
    mergedResults <- c()
    for (k in seq_len(nrow(windowDividedRes))) {
        imputedK <- dentistResultByWindow[[k]]
        imputedK$index_within_window <- seq_len(nrow(imputedK))
        imputedK <- imputedK |>
            mutate(
                index_global = .data$index_within_window +
                    windowDividedRes$windowStartIdx[k] -
                    1
            )
        extractedResults <- imputedK |>
            filter(
                .data$index_global >= windowDividedRes$fillStartIdx[k] &
                    .data$index_global < windowDividedRes$fillEndIdx[k]
            )
        mergedResults <- bind_rows(mergedResults, extractedResults)
    }
    return(mergedResults)
}

# ## File-I/O functions (dentist_from_files, read_dentist_sumstat,
# parse_dentist_output)
### have been removed. Use the standard pipeline: load genotypes via
### loadGenotypeRegion(), compute LD via computeLd(), then call dentist()
### or ldMismatchQc() directly.

# =============================================================================
# SLALoM: Approximate Bayes Factor single-causal-variant outlier detection
# =============================================================================

# Numerically stable log-sum-exp.
# @noRd
.logSumExp <- function(x) {
    maxX <- max(x, na.rm = TRUE)
    sumExp <- sum(exp(x - maxX), na.rm = TRUE)
    return(maxX + log(sumExp))
}

# Approximate (Wakefield) Bayes factors from z-scores and standard errors:
# per-variant log-BF and normalized posterior probabilities.
# @noRd
.approxBayesFactor <- function(z, se, W = 0.04) {
    V <- se^2
    r <- W / (W + V)
    lbf <- 0.5 * (log(1 - r) + (r * z^2))
    denom <- .logSumExp(lbf)
    prob <- exp(lbf - denom)
    return(list(lbf = lbf, prob = prob))
}

# Greedy credible set: indices (ordered by decreasing prob) whose cumulative
# posterior first exceeds `coverage`.
# @noRd
.slalomCredibleSet <- function(prob, coverage = 0.95) {
    ordering <- order(prob, decreasing = TRUE)
    cumprob <- cumsum(prob[ordering])
    idx <- which(cumprob > coverage)[1]
    cs <- ordering[seq_len(idx)]
    return(cs)
}

# --- slalom helpers ---------------------------------------------------------

# Require exactly one of R (LD matrix) or X (genotype matrix).
.slalomValidate <- function(R, X) {
    if (is.null(R) && is.null(X)) {
        abort("Either R (LD matrix) or X (genotype matrix) must be provided.")
    }
    if (!is.null(R) && !is.null(X)) {
        abort("Provide either R or X, not both.")
    }
}

# Lead-variant LD column. From X, compute just that column (avoid full p x p).
.slalomLeadR <- function(zScore, R, X, leadIdx) {
    if (!is.null(X)) {
        if (!is.matrix(X)) {
            X <- as.matrix(X)
        }
        return(as.numeric(cor(X, X[, leadIdx])))
    }
    if (!is.matrix(R) || nrow(R) != ncol(R) || nrow(R) != length(zScore)) {
        abort("R must be a square matrix matching the length of zScore.")
    }
    R[, leadIdx]
}

# DENTIST-S outlier statistic against the lead variant.
.slalomDentistS <- function(
    zScore,
    rLead,
    leadIdx,
    r2Threshold,
    nlog10pDentistSThreshold
) {
    r2Lead <- rLead^2
    tDentistS <- (zScore - rLead * zScore[leadIdx])^2 / (1 - r2Lead)
    tDentistS[tDentistS < 0] <- Inf
    nlog10pDentistS <- -log10(pchisq(tDentistS, df = 1, lower.tail = FALSE))
    outliers <- (r2Lead > r2Threshold) &
        (nlog10pDentistS > nlog10pDentistSThreshold)
    list(
        r2Lead = r2Lead,
        nlog10pDentistS = nlog10pDentistS,
        outliers = outliers
    )
}

# Assemble the SLALOM per-variant table and summary.
.slalomResult <- function(
    zScore,
    prob,
    pvalue,
    ds,
    leadIdx,
    cs,
    cs99,
    r2Threshold
) {
    nR2 <- sum(ds$r2Lead > r2Threshold)
    nDentistSOutlier <- sum(ds$outliers, na.rm = TRUE)
    summary <- list(
        leadPipVariant = leadIdx,
        nTotal = length(zScore),
        nR2 = nR2,
        nDentistSOutlier = nDentistSOutlier,
        fraction = if_else(nR2 > 0, nDentistSOutlier / nR2, 0),
        maxPip = max(prob),
        cs95 = cs,
        cs99 = cs99
    )
    result <- tibble(
        original_z = zScore,
        prob = prob,
        pvalue = pvalue,
        outliers = ds$outliers,
        nlog10p_dentist_s = ds$nlog10pDentistS
    )
    list(data = result, summary = summary)
}

#' Slalom Function for Summary Statistics QC for Fine-Mapping Analysis
#'
#' Performs Approximate Bayesian Factor (ABF) analysis, identifies credible
#' sets, and annotates lead variants based on fine-mapping results. It computes
#' p-values from z-scores assuming a two-sided standard normal distribution.
#'
#' Provide either an LD correlation matrix \code{R} or a genotype matrix
#' \code{X} (from which LD is derived automatically via \code{computeLd}).
#'
#' @param zScore Numeric vector of z-scores corresponding to each variant.
#' @param R Square LD correlation matrix. Provide either \code{R} or \code{X}.
#' @param X Genotype matrix (samples x SNPs). If provided, LD is computed via
#'   \code{computeLd(X)}.
#' @param standardError Optional numeric vector of standard errors corresponding
#'   to each z-score. If not provided, a default value of 1 is assumed for all
#'   variants.
#' @param abfPriorVariance Numeric, the prior effect size variance for ABF
#'   calculations. Default is 0.04.
#' @param nlog10pDentistSThreshold Numeric, the -log10 DENTIST-S P value
#'   threshold for identifying outlier variants for prediction. Default is 4.0.
#' @param r2Threshold Numeric, the r2 threshold for DENTIST-S outlier variants
#'   for prediction. Default is 0.6.
#' @param leadVariantChoice Character, method to choose the lead variant, either
#'   "pvalue" or "abf", with default "pvalue".
#' @param ldMethod Character string specifying the LD computation method when
#'   \code{X} is provided. Passed to \code{computeLd}. One of \code{"sample"}
#'   (default), \code{"population"}, or \code{"gcta"}. Ignored when \code{R} is
#'   provided directly.
#' @return A list containing the annotated LD matrix with ABF results, credible
#'   sets, lead variant, and DENTIST-S statistics; and a summary dataframe with
#'   aggregate statistics.
#' @examples
#' # Simulate z-scores and an LD matrix from a reference panel, then screen
#' # for LD-outlier variants.
#' set.seed(1)
#' geno <- matrix(rbinom(500 * 100, 2, 0.3), nrow = 500, ncol = 100)
#' ldMat <- cor(geno)
#' zScore <- rnorm(100)
#' result <- slalom(zScore, R = ldMat)
#' @seealso \code{\link{dentistSingleWindow}}, \code{resolveLdInput}
#' @importFrom stats pchisq
#' @export
#'
slalom <- function(
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
    .slalomValidate(R, X)
    # One-sided p-value matching the original Python (stats.norm.cdf): lead is
    # the most negative z-score when leadVariantChoice == "pvalue".
    pvalue <- pnorm(zScore)
    abfResults <- .approxBayesFactor(
        zScore,
        standardError,
        W = abfPriorVariance
    )
    prob <- abfResults$prob
    cs <- .slalomCredibleSet(prob, coverage = 0.95)
    cs99 <- .slalomCredibleSet(prob, coverage = 0.99)
    leadIdx <- if (leadVariantChoice == "pvalue") {
        which.min(pvalue)
    } else {
        which.max(prob)
    }
    rLead <- .slalomLeadR(zScore, R, X, leadIdx)
    ds <- .slalomDentistS(
        zScore,
        rLead,
        leadIdx,
        r2Threshold,
        nlog10pDentistSThreshold
    )
    .slalomResult(zScore, prob, pvalue, ds, leadIdx, cs, cs99, r2Threshold)
}


# =============================================================================
# Univariate RSS diagnostics (post-finemap)
# =============================================================================

#' Extract the trimmed SuSiE fit from a finemapping pipeline result
#'
#' Returns the trimmed model fit underlying \code{con_data$finemappingEntry} (a
#' \code{FineMappingRow} S4 object), or NULL if no fine-mapping entry is
#' attached.
#'
#' @param conData List. The method-layer entry from a finemapping pipeline
#'   result, expected to carry \code{$finemappingEntry} as a
#'   \code{FineMappingRow} object.
#' @return The trimmed fit (a list with \code{pip}, \code{sets}, etc.) or NULL.
#' @examples
#' data(qtlSumStatsExample)
#' getSusieResult(qtlSumStatsExample)
#' @export
getSusieResult <- function(conData) {
    if (length(conData) == 0) {
        return(NULL)
    }
    fm <- conData$finemappingEntry
    if (is.null(fm) || !is(fm, "FineMappingResultBase")) {
        return(NULL)
    }
    trimmed <- .fmrPartsSusieFit(fm)
    if (length(trimmed) == 0) {
        return(NULL)
    }
    trimmed
}

#' Process Credible Sets (CS) from Finemapping Results
#'
#' This function extracts and processes information for each Credible Set (CS)
#' from finemapping results, typically obtained from a finemapping RDS file.
#'
#' @param fmRow A \code{\link{fineMappingRow}} carrying the SuSiE
#'   fit and variant ids (e.g. from \code{\link{getFineMappingResult}}).
#' @param csNames Character vector. Names of the Credible Sets, usually in the
#'   format "L_<number>".
#' @param topLociTable Data frame. The top-loci table (e.g. from
#'   \code{\link{getTopLoci}}) carrying \code{variant_id}, \code{pip}, and
#'   \code{z} columns.
#' @param ldSource The LD source from which the between-credible-set correlation
#'   is derived on demand: a \code{QtlDataset} (individual-level) or a
#'   \code{QtlSumStats} / \code{GwasSumStats} (summary statistics). See
#'   \code{\link{computeCsCorrelation}}.
#'
#' @return A data frame with one row per CS, containing the following columns:
#'   \item{cs_name}{Name of the Credible Set}
#'   \item{variants_per_cs}{Number of variants in the CS}
#'   \item{top_variant}{ID of the variant with the highest PIP in the CS}
#'   \item{top_variant_index}{Global index of the top variant}
#'   \item{top_pip}{Highest Posterior Inclusion Probability (PIP) in the CS}
#'   \item{top_z}{Z-score of the top variant}
#'   \item{p_value}{P-value calculated from the top Z-score}
#'   \item{cs_corr_1, cs_corr_2, ...}{Each CS's pairwise correlation with every
#'     CS (its row of the between-CS matrix, self-correlation on the diagonal),
#'     computed on demand from \code{ldSource}. Absent when there are fewer than
#'     two credible sets.}
#'   \item{cs_corr_max}{Maximum absolute between-CS correlation (excluding the
#'     self == 1); \code{NA} for a single CS.}
#'   \item{cs_corr_min}{Minimum absolute between-CS correlation; \code{NA} for a
#'     single CS.}
#'
#' @details This function is designed to be used only when there is at least one
#'   Credible Set in the finemapping results usually for a given study and
#'   block. It processes each CS, extracting key information such as the top
#'   variant, its statistics, and correlation information between multiple CS if
#'   available.
#'
#' @importFrom purrr map map_dbl map_int
#' @importFrom dplyr bind_rows
#'
#' @examples
#' data(qtlSumStatsExample)
#' vids <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")
#' fit <- list(pip = c(0.1, 0.7, 0.2), sets = list(cs = list(L_1 = c(1, 2))))
#' tl <- data.frame(variant_id = vids, pip = c(0.1, 0.7, 0.2),
#'   z = c(1.0, 3.5, -0.5))
#' fe <- fineMappingRow(variantIds = vids, susieFit = fit, topLoci = tl)
#' # A single credible set has no between-CS correlation (cs_corr_* are NA), so
#' # the ldSource is not consulted here.
#' extractCsInfo(fe, csNames = "L_1", topLociTable = tl,
#'   ldSource = qtlSumStatsExample)
#'
#' @export
extractCsInfo <- function(fmRow, csNames, topLociTable, ldSource) {
    fm <- fmRow
    trimmed <- .fmrPartsSusieFit(fm)
    variantNames <- .fmrPartsVariantIds(fm)
    csCorr <- .rowCsCorrelation(fm, ldSource)
    rows <- map(
        seq_along(csNames),
        .extractCsInfoRow,
        csNames = csNames,
        trimmed = trimmed,
        variantNames = variantNames,
        topLociTable = topLociTable
    )
    .csAppendCorrelationCols(bind_rows(rows), csCorr)
}

#' Extract Information for Top Variant from Finemapping Results
#'
#' This function extracts information about the variant with the highest
#' Posterior Inclusion Probability (PIP) from finemapping results, typically
#' used when no Credible Sets (CS) are identified in the analysis.
#'
#' @param fmRow A \code{\link{fineMappingRow}} carrying the SuSiE
#'   fit and variant ids (e.g. from \code{\link{getFineMappingResult}}).
#' @param sumstats A list or data frame carrying a \code{z} element aligned to
#'   the fit's variants (\code{sumstats$z}).
#'
#' @return A data frame with one row containing the following columns:
#'   \item{cs_name}{NA (as no CS is identified)}
#'   \item{variants_per_cs}{NA (as no CS is identified)}
#'   \item{top_variant}{ID of the variant with the highest PIP}
#'   \item{top_variant_index}{Index of the top variant in the original data}
#'   \item{top_pip}{Highest Posterior Inclusion Probability (PIP)}
#'   \item{top_z}{Z-score of the top variant}
#'   \item{p_value}{P-value calculated from the top Z-score}
#'   \item{cs_corr_max}{NA (no between-CS correlation without a CS)}
#'   \item{cs_corr_min}{NA (no between-CS correlation without a CS)}
#'
#' @details This function is designed to be used when no Credible Sets are
#'   identified in the finemapping results, but information about the most
#'   significant variant is still desired. It identifies the variant with the
#'   highest PIP and extracts relevant statistical information.
#'
#' @note This function is particularly useful for capturing information about
#'   potentially important variants that might be included in Credible Sets
#'   under different analysis parameters or lower coverage. It maintains a
#'   structure similar to the output of `extract_cs_info()` for consistency in
#'   downstream analyses.
#'
#' @seealso \code{\link{extractCsInfo}} for processing when Credible Sets are
#'   present.
#'
#' @examples
#' vids <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")
#' fit <- list(pip = c(0.1, 0.7, 0.2))
#' tl <- data.frame(variant_id = vids, pip = c(0.1, 0.7, 0.2))
#' fe <- fineMappingRow(variantIds = vids, susieFit = fit, topLoci = tl)
#' extractTopPipInfo(fe, sumstats = list(z = c(1.0, 3.5, -0.5)))
#'
#' @export
extractTopPipInfo <- function(fmRow, sumstats) {
    fm <- fmRow
    trimmed <- .fmrPartsSusieFit(fm)
    variantNames <- .fmrPartsVariantIds(fm)
    # Find the variant with the highest PIP
    topPipIndex <- which.max(trimmed$pip)
    topPip <- trimmed$pip[topPipIndex]
    topVariant <- variantNames[topPipIndex]
    topZ <- sumstats$z[topPipIndex]
    pValue <- .zToPvalue(topZ)

    list(
        cs_name = NA,
        variants_per_cs = NA,
        top_variant = topVariant,
        top_variant_index = topPipIndex,
        top_pip = topPip,
        top_z = topZ,
        p_value = pValue,
        cs_corr_max = NA_real_,
        cs_corr_min = NA_real_
    )
}

# Reduce one credible set's correlation vector (a row of the between-CS matrix,
# the self-correlation == 1 on the diagonal) to its |corr| max/min, excluding
# every self / perfect correlation (== 1). An empty result yields NA.
# @noRd
.extractCorrelations <- function(x) {
    filtered <- abs(x[x != 1])
    if (length(filtered) == 0L) {
        return(list(max_corr = NA_real_, min_corr = NA_real_))
    }
    list(
        max_corr = max(filtered, na.rm = TRUE),
        min_corr = min(filtered, na.rm = TRUE)
    )
}

# Append the between-CS correlation columns to the per-CS summary `base` from
# the m x m matrix `csCorr` (whose rows are aligned to `base`): the expanded
# cs_corr_1..m (each CS's row, self-correlation on the diagonal) plus
# cs_corr_max / cs_corr_min (|corr| excluding the self == 1). A NULL matrix
# (fewer than two credible sets) yields NA max/min and no expanded columns.
# @noRd
.csAppendCorrelationCols <- function(base, csCorr) {
    if (is.null(csCorr)) {
        return(mutate(base, cs_corr_max = NA_real_, cs_corr_min = NA_real_))
    }
    perRow <- apply(csCorr, 1L, .extractCorrelations, simplify = FALSE)
    expanded <- as_tibble(csCorr, .name_repair = "minimal")
    names(expanded) <- str_c("cs_corr_", seq_len(ncol(csCorr)))
    # unname(): apply() names its result by the matrix rownames, which map_dbl
    # then carries into the column (tibbles preserve element names).
    base |>
        bind_cols(expanded) |>
        mutate(
            cs_corr_max = unname(map_dbl(perRow, "max_corr")),
            cs_corr_min = unname(map_dbl(perRow, "min_corr"))
        )
}

#' Process Credible Set Information and Determine Updating Strategy
#'
#' This function categorizes Credible Sets (CS) within a study block into
#' different updating strategies based on their statistical properties and
#' correlations.
#'
#' @param df Data frame. Contains information about Credible Sets for a specific
#'   study and block.
#' @param highCorrCols Character vector. Names of columns in df that represent
#'   high correlations.
#'
#' @return A modified data frame with additional columns attached to the
#' diagnostic table:
#'   \item{top_cs}{Logical. TRUE for the CS with the highest absolute Z-score.}
#'   \item{tagged_cs}{Logical. TRUE for CS that are considered "tagged" based
#'   on p-value and correlation criteria.}
#'   \item{method}{Character. The determined updating strategy ("BVSR", "SER",
#'   or "BCR").}
#'
#' @details This function performs the following steps: 1. Identifies the top CS
#'   based on the highest absolute Z-score. 2. Identifies tagged CS based on
#'   high p-value and high correlations. 3. Counts total, tagged, and remaining
#'   CS. 4. Determines the appropriate updating method based on these counts.
#'
#' The updating methods are: - BVSR (Bayesian Variable Selection Regression):
#' Used when there's only one CS or all CS are accounted for. - SER (Single
#' Effect Regression): Used when there are tagged CS but no remaining untagged
#' CS. - BCR (Bayesian Conditional Regression): Used when there are remaining
#' untagged CS.
#'
#' @note This function is part of a developing methodology for automatically
#'   handling finemapping results. The thresholds and criteria used (e.g.,
#'   p-value > 1e-4 for tagging) are subject to refinement and may change in
#'   future versions.
#'
#' @importFrom dplyr case_when filter select all_of row_number
#'
#' @examples
#' df <- data.frame(cs_name = c("L1", "L2"), top_z = c(5, 3.5),
#'   p_value = c(1e-10, 1e-6))
#' autoDecision(df, highCorrCols = character(0))
#' @export
autoDecision <- function(df, highCorrCols) {
    # Identify top_cs
    topCsIndex <- which.max(abs(df$top_z))
    df$top_cs <- FALSE
    df$top_cs[topCsIndex] <- TRUE

    # Identify tagged_cs
    df$tagged_cs <- map_lgl(
        seq_len(nrow(df)),
        .autoDecisionTagged,
        df = df,
        highCorrCols = highCorrCols
    )

    # Count total and remaining CS
    totalCs <- nrow(df)
    taggedCsCount <- sum(df$tagged_cs)
    if (totalCs > 0) {
        remainingCs <- totalCs - 1 - taggedCsCount
    } else {
        remainingCs <- 0
    }
    # Determine method
    df$method <- case_when(
        taggedCsCount == 0 & totalCs > 1 ~ "BVSR",
        (remainingCs == 0 & totalCs > 1) | (totalCs == 1) ~ "SER",
        remainingCs > 0 ~ "BCR",
        TRUE ~ NA_character_
    )

    return(df)
}


# =============================================================================
# RAISS: regression-based sumstats imputation
# =============================================================================

#' Core RAISS implementation for a single LD matrix
#'
#' @param refPanel A data frame containing 'chrom', 'pos', 'variant_id', 'A1',
#'   and 'A2'.
#' @param knownZscores A data frame containing 'chrom', 'pos', 'variant_id',
#'   'A1', 'A2', and 'z' values.
#' @param ldMatrix A square matrix of dimension equal to the number of rows in
#'   refPanel.
#' @param lamb Regularization term added to the diagonal of the ldMatrix.
#' @param rcond Threshold for filtering eigenvalues in the pseudo-inverse
#'   computation.
#' @param r2Threshold R square threshold below which SNPs are filtered from the
#'   output.
#' @param minimumLd Minimum LD score threshold for SNP filtering.
#' @param verbose Logical indicating whether to print progress information.
#'
#' @return A list containing filtered and unfiltered results, and filtered LD
#'   matrix.
#' @importFrom MASS ginv
#' @importFrom dplyr arrange
#' @noRd
# Drop LD-matrix rows/cols for variants filtered out of the imputation result.
.raissFilterLd <- function(ldMatrix, refPanel, resultFilter) {
    filteredOutVariant <- setdiff(refPanel$variant_id, resultFilter$variant_id)
    if (length(filteredOutVariant) > 0) {
        filteredOutId <- match(filteredOutVariant, refPanel$variant_id)
        as.matrix(ldMatrix)[-filteredOutId, -filteredOutId]
    } else {
        as.matrix(ldMatrix)
    }
}

raissSingleMatrix <- function(
    refPanel,
    knownZscores,
    ldMatrix,
    lamb = 0.01,
    rcond = 0.01,
    r2Threshold = 0.6,
    minimumLd = 5,
    verbose = TRUE
) {
    if (is.unsorted(refPanel$pos) || is.unsorted(knownZscores$pos)) {
        abort("refPanel and knownZscores must be in increasing order of pos.")
    }
    if (is.data.frame(ldMatrix)) {
        ldMatrix <- as.matrix(ldMatrix)
    }
    knownsId <- intersect(knownZscores$variant_id, refPanel$variant_id)
    knowns <- which(is_in(refPanel$variant_id, knownsId))
    unknowns <- which(!is_in(refPanel$variant_id, knownsId))
    edge <- .raissEdgeCase(knowns, unknowns, knownZscores, verbose, ldMatrix)
    if (!is.null(edge)) {
        return(edge$value)
    }
    zt <- knownZscores$z
    sigT <- ldMatrix[knowns, knowns, drop = FALSE]
    sigIT <- ldMatrix[unknowns, knowns, drop = FALSE]
    results <- raissModel(zt, sigT, sigIT, lamb, rcond)
    results <- formatRaissDf(results, refPanel, unknowns)
    results <- filterRaissOutput(results, r2Threshold, minimumLd, verbose)
    resultNofilter <- mergeRaissDf(results$zscoresNofilter, knownZscores) |>
        arrange(.data$pos)
    resultFilter <- mergeRaissDf(results$zscores, knownZscores) |>
        arrange(.data$pos)
    list(
        resultNofilter = resultNofilter,
        resultFilter = resultFilter,
        ldMat = .raissFilterLd(ldMatrix, refPanel, resultFilter)
    )
}

#' Core RAISS implementation from a genotype matrix X (SVD-based)
#'
#' Performs the same imputation as \code{raissSingleMatrix} but works directly
#' with the genotype matrix X instead of the LD correlation matrix R. This
#' avoids forming the p x p LD matrix, saving O(p^2) memory and O(np^2) compute.
#'
#' The reformulation is mathematically exact: using the thin SVD of Xt (the
#' known variant columns), all RAISS quantities (mu, var, ld_score) are computed
#' in the SVD basis without ever forming R = X'X/(n-1).
#'
#' @param refPanel A data frame containing 'chrom', 'pos', 'variant_id', 'A1',
#'   and 'A2'.
#' @param knownZscores A data frame containing 'chrom', 'pos', 'variant_id',
#'   'A1', 'A2', and 'z' values.
#' @param X Centered and scaled genotype matrix (nSamples x pVariants). Column
#'   order must match the variant order in refPanel.
#' @param lamb Regularization term (same role as in the LD-based path).
#' @param svdTol Relative tolerance for filtering small singular values in the
#'   SVD of Xt.
#' @param r2Threshold R square threshold below which SNPs are filtered from the
#'   output.
#' @param minimumLd Minimum LD score threshold for SNP filtering.
#' @param verbose Logical indicating whether to print progress information.
#'
#' @return A list containing filtered and unfiltered results, and ldMat = NULL.
#' @importFrom dplyr arrange
#' @noRd
# --- raissSingleMatrixFromX helpers ----------------------------------------

# Early-return cases: no known variants (NULL) or nothing to impute (knowns).
# Returns NULL to continue, or list(value = <return value>) to short-circuit.
.raissEdgeCase <- function(
    knowns,
    unknowns,
    knownZscores,
    verbose,
    ldMat = NULL
) {
    if (length(knowns) == 0) {
        if (verbose) {
            inform("No known variants found, cannot perform imputation.")
        }
        return(list(value = NULL))
    }
    if (length(unknowns) == 0) {
        if (verbose) {
            inform("No unknown variants to impute, returning known variants.")
        }
        return(list(
            value = list(
                resultNofilter = knownZscores,
                resultFilter = knownZscores,
                ldMat = ldMat
            )
        ))
    }
    NULL
}

# SVD-based imputation from the genotype matrix (avoids forming R). Computes
# X' %*% [w | U] in a single BLAS call to skip an O(n*m) copy of X_unknown.
.raissSvdImpute <- function(X, knowns, unknowns, zt, lamb, svdTol, nSamples) {
    # Thin SVD of the known columns (n x k -> U: n x r, d: r, V: k x r).
    Xt <- X[, knowns, drop = FALSE]
    svdResult <- .safeSvd(Xt, tol = svdTol)
    U <- svdResult$u
    d <- svdResult$d
    V <- svdResult$v
    rm(Xt)
    cReg <- lamb * (nSamples - 1)
    d2 <- d^2
    d2PlusC <- d2 + cReg
    # w = U %*% diag(d / (d^2 + c)) %*% V' zt  (projection of zt through SVD).
    VtZt <- crossprod(V, zt)
    w <- U %*% (d / d2PlusC * VtZt)
    # Single dgemm X' %*% [w | U]: col 1 (rows unknowns) -> mu; rest -> A.
    XtWU <- crossprod(X, cbind(w, U))
    mu <- as.numeric(XtWU[unknowns, 1])
    A <- XtWU[unknowns, -1, drop = FALSE]
    rm(XtWU)
    # Variance and LD score in one pass over A^2.
    ASq <- A^2
    rm(A)
    dWeights <- cbind(d2 / d2PlusC, d2)
    scores <- ASq %*% dWeights
    rm(ASq)
    nm1 <- nSamples - 1
    varRaw <- (1 + lamb) - scores[, 1] / nm1
    raissLdScore <- scores[, 2] / nm1^2
    rm(scores)
    conditionNumber <- rep(d[1] / d[length(d)], length(unknowns))
    correctInversion <- rep(TRUE, length(unknowns))
    # R2 correction (same as raissModel).
    varNorm <- varInBoundaries(varRaw, lamb)
    R2 <- (1 + lamb) - varNorm
    mu <- mu / sqrt(R2)
    list(
        var = varNorm,
        mu = mu,
        raissLdScore = raissLdScore,
        conditionNumber = conditionNumber,
        correctInversion = correctInversion
    )
}

raissSingleMatrixFromX <- function(
    refPanel,
    knownZscores,
    X,
    lamb = 0.01,
    svdTol = 1e-8,
    r2Threshold = 0.6,
    minimumLd = 5,
    verbose = TRUE
) {
    if (is.unsorted(refPanel$pos) || is.unsorted(knownZscores$pos)) {
        abort("refPanel and knownZscores must be in increasing order of pos.")
    }
    knownsId <- intersect(knownZscores$variant_id, refPanel$variant_id)
    knowns <- which(is_in(refPanel$variant_id, knownsId))
    unknowns <- which(!is_in(refPanel$variant_id, knownsId))
    edge <- .raissEdgeCase(knowns, unknowns, knownZscores, verbose)
    if (!is.null(edge)) {
        return(edge$value)
    }
    imp <- .raissSvdImpute(
        X,
        knowns,
        unknowns,
        knownZscores$z,
        lamb,
        svdTol,
        nrow(X)
    )
    results <- formatRaissDf(imp, refPanel, unknowns)
    results <- filterRaissOutput(results, r2Threshold, minimumLd, verbose)
    resultNofilter <- mergeRaissDf(results$zscoresNofilter, knownZscores) |>
        arrange(.data$pos)
    resultFilter <- mergeRaissDf(results$zscores, knownZscores) |>
        arrange(.data$pos)
    list(
        resultNofilter = resultNofilter,
        resultFilter = resultFilter,
        ldMat = NULL
    )
}

# Sequentially append RAISS block results, resolving the shared boundary variant
# (duplicated across adjacent blocks) by keeping the higher-R2 imputation.
# @noRd
.combineWithBoundaryCheck <- function(combinedResult, newResult) {
    # If either is empty, simply return the non-empty one or empty data frame
    if (is.null(combinedResult)) {
        return(newResult)
    }
    if (is.null(newResult)) {
        return(combinedResult)
    }

    # Check if the last variant of combined matches the first of new
    lastVar <- combinedResult$variant_id[nrow(combinedResult)]
    firstVar <- newResult$variant_id[1]

    if (lastVar == firstVar) {
        newR2 <- newResult$raissR2[1]
        oldR2 <- combinedResult$raissR2[nrow(combinedResult)]
        if (is.na(newR2) && is.na(oldR2)) {
            # Both are NA - keep the existing one
        } else if (is.na(oldR2)) {
            # Old is NA but new is not - use new
            combinedResult[nrow(combinedResult), ] <- newResult[1, ]
        } else if (is.na(newR2)) {
            # New is NA but old is not - keep old
        } else if (newR2 > oldR2) {
            # Both are non-NA and new is better - use new
            combinedResult[nrow(combinedResult), ] <- newResult[1, ]
        }

        # Add remaining rows from new (excluding first)
        if (nrow(newResult) > 1) {
            combinedResult <- bind_rows(combinedResult, newResult[-1, ])
        }
    } else {
        # No overlap - combine all rows
        combinedResult <- bind_rows(combinedResult, newResult)
    }

    return(combinedResult)
}

# --- raiss: genotype-matrix and LD-matrix path helpers ---------------------

# List of genotype-matrix blocks: impute each via SVD, then row-bind.
.raissGenotypeBlocks <- function(refPanel, knownZscores, genotypeMatrix, p) {
    if (p$verbose) {
        msg <- glue(
            "Processing multiple genotype matrix blocks via SVD-based ",
            "imputation..."
        )
        inform(msg)
    }
    resultsList <- list()
    for (i in seq_along(genotypeMatrix)) {
        if (p$verbose) {
            nBlocks <- length(genotypeMatrix)
            msg <- glue("Processing block {i} of {nBlocks}")
            inform(msg)
        }
        blockResult <- raissSingleMatrixFromX(
            refPanel,
            knownZscores,
            genotypeMatrix[[i]],
            p$lamb,
            p$svdTol,
            p$r2Threshold,
            p$minimumLd,
            verbose = FALSE
        )
        if (!is.null(blockResult)) {
            resultsList[[length(resultsList) + 1]] <- blockResult
        }
    }
    if (length(resultsList) == 0) {
        if (p$verbose) {
            inform("No blocks could be processed.")
        }
        return(NULL)
    }
    combinedNofilter <- bind_rows(map(resultsList, "resultNofilter"))
    combinedFilter <- bind_rows(map(resultsList, "resultFilter"))
    list(
        resultNofilter = combinedNofilter |> arrange(.data$pos),
        resultFilter = combinedFilter |> arrange(.data$pos),
        ldMat = NULL
    )
}

# Genotype-matrix imputation path: single matrix, list of blocks, or error.
.raissGenotypePath <- function(refPanel, knownZscores, genotypeMatrix, p) {
    if (is.matrix(genotypeMatrix)) {
        if (p$verbose) {
            inform("Processing genotype matrix via SVD-based imputation...")
        }
        return(raissSingleMatrixFromX(
            refPanel,
            knownZscores,
            genotypeMatrix,
            p$lamb,
            p$svdTol,
            p$r2Threshold,
            p$minimumLd,
            p$verbose
        ))
    }
    if (is.list(genotypeMatrix)) {
        return(.raissGenotypeBlocks(refPanel, knownZscores, genotypeMatrix, p))
    }
    abort("genotypeMatrix must be a matrix or a list of matrices.")
}

# Impute one LD block against its reference-panel subset.
.raissLdBlockOne <- function(
    refPanel,
    knownZscores,
    ldMatrix,
    variantIndices,
    blockId,
    p
) {
    blockVariantIds <- variantIndices$variant_id[
        variantIndices$blockId == blockId
    ]
    blockIndices <- match(blockVariantIds, refPanel$variant_id)
    blockRefPanel <- refPanel[blockIndices, ]
    blockLdMatrix <- ldMatrix$ldMatrices[[blockId]]
    blockKnownZscores <- knownZscores |>
        filter(is_in(.data$variant_id, blockVariantIds))
    if (nrow(blockLdMatrix) != nrow(blockRefPanel)) {
        msg <- glue(
            "Block {blockId} : LD matrix dimension does not match number ",
            "of variants in reference panel"
        )
        abort(msg)
    }
    raissSingleMatrix(
        blockRefPanel,
        blockKnownZscores,
        blockLdMatrix,
        p$lamb,
        p$rcond,
        p$r2Threshold,
        p$minimumLd,
        verbose = FALSE
    )
}

# Combine per-block imputation results (boundary-dedup) + rebuild LD matrix.
.raissLdBlocksCombine <- function(resultsList) {
    combinedNofilter <- resultsList[[1]]$resultNofilter
    combinedFilter <- resultsList[[1]]$resultFilter
    if (length(resultsList) > 1) {
        for (i in 2:length(resultsList)) {
            combinedNofilter <- .combineWithBoundaryCheck(
                combinedNofilter,
                resultsList[[i]]$resultNofilter
            )
            combinedFilter <- .combineWithBoundaryCheck(
                combinedFilter,
                resultsList[[i]]$resultFilter
            )
        }
    }
    ldFilteredList <- map(resultsList, "ldMat")
    variantList <- map(ldFilteredList, .ldVariantsDf)
    ldMatrix <- createLdMatrix(
        ldMatrices = ldFilteredList,
        variants = variantList
    )
    list(
        resultNofilter = combinedNofilter,
        resultFilter = combinedFilter,
        ldMat = ldMatrix
    )
}

# LD-block imputation path: impute each block then combine.
.raissLdBlocksPath <- function(refPanel, knownZscores, ldMatrix, p) {
    if (p$verbose) {
        inform("Processing multiple LD blocks...")
    }
    variantIndices <- ldMatrix$variantIndices
    blockIds <- unique(variantIndices$blockId)
    resultsList <- list()
    for (blockId in blockIds) {
        if (p$verbose) {
            nBlocks <- length(blockIds)
            msg <- glue("Processing block {blockId} of {nBlocks}")
            inform(msg)
        }
        blockResult <- .raissLdBlockOne(
            refPanel,
            knownZscores,
            ldMatrix,
            variantIndices,
            blockId,
            p
        )
        if (!is.null(blockResult)) resultsList[[blockId]] <- blockResult
    }
    if (length(resultsList) == 0) {
        if (p$verbose) {
            msg <- glue(
                "No blocks could be processed. Check that knownZscores ",
                "overlap with variants in the blocks."
            )
            inform(msg)
        }
        return(NULL)
    }
    .raissLdBlocksCombine(resultsList)
}

# Single LD-matrix imputation path (extracts matrix from a 1-block list).
.raissSingleLdPath <- function(refPanel, knownZscores, ldMatrix, p) {
    if (p$verbose) {
        fromList <- if (!is.matrix(ldMatrix)) " from list" else ""
        msg <- glue("Processing single LD matrix{fromList}...")
        inform(msg)
    }
    if (!is.matrix(ldMatrix)) {
        ldMatrix <- ldMatrix$ldMatrices[[1]]
    }
    raissSingleMatrix(
        refPanel,
        knownZscores,
        ldMatrix,
        p$lamb,
        p$rcond,
        p$r2Threshold,
        p$minimumLd,
        p$verbose
    )
}

#' Impute Summary Statistics Using LD (RAISS)
#'
#' This function is a part of the statistical library for SNP imputation from:
#' https://gitlab.pasteur.fr/statistical-genetics/raiss/-/blob/master/raiss/stat_models.py
#' It is R implementation of the imputation model described in the paper by
#' Bogdan Pasaniuc, Noah Zaitlen, et al., titled "Fast and accurate imputation
#' of summary statistics enhances evidence of functional enrichment", published
#' in Bioinformatics in 2014.
#'
#' This function can process either a single LD matrix or a list of LD matrices
#' for different blocks. For a list of matrices, it processes each block
#' separately and combines the results. Alternatively, it can accept a genotype
#' matrix X directly, avoiding the need to form the p x p LD matrix (memory and
#' compute savings when n << p).
#'
#' @param refPanel A data frame containing 'chrom', 'pos', 'variant_id', 'A1',
#'   and 'A2'.
#' @param knownZscores A data frame containing 'chrom', 'pos', 'variant_id',
#'   'A1', 'A2', and 'z' values.
#' @param ldMatrix Either a square matrix or a list of matrices for LD blocks.
#'   Provide either \code{ldMatrix} or \code{genotypeMatrix}, not both.
#' @param genotypeMatrix A centered and scaled genotype matrix (n x p) as an
#'   alternative to \code{ldMatrix}. Column order must match the variant order
#'   in \code{refPanel}. When provided, the imputation uses an SVD-based
#'   approach that avoids forming the p x p LD matrix.
#' @param lamb Regularization term added to the diagonal of the ldMatrix.
#' @param rcond Threshold for filtering eigenvalues in the pseudo-inverse
#'   computation (only used with ldMatrix path).
#' @param svdTol Relative tolerance for filtering small singular values (only
#'   used with genotypeMatrix path).
#' @param r2Threshold R square threshold below which SNPs are filtered from the
#'   output.
#' @param minimumLd Minimum LD score threshold for SNP filtering.
#' @param verbose Logical indicating whether to print progress information.
#'
#' @return A list containing filtered and unfiltered results, and filtered LD
#'   matrix (ldMat is NULL when using genotypeMatrix path).
#' @importFrom dplyr arrange bind_rows
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:30]
#' refPanel <- data.frame(chrom = 22, pos = 1:30,
#'   variant_id = paste0("22:", 1:30, ":A:G"), A1 = "A", A2 = "G")
#' knownZscores <- data.frame(chrom = 22, pos = 1:20,
#'   variant_id = paste0("22:", 1:20, ":A:G"), A1 = "A", A2 = "G",
#'   z = rnorm(20))
#' raiss(refPanel = refPanel, knownZscores = knownZscores,
#'   ldMatrix = cor(X))
#' @export
raiss <- function(
    refPanel,
    knownZscores,
    ldMatrix = NULL,
    genotypeMatrix = NULL,
    lamb = 0.01,
    rcond = 0.01,
    svdTol = 1e-8,
    r2Threshold = 0.6,
    minimumLd = 5,
    verbose = TRUE
) {
    p <- list(
        lamb = lamb,
        rcond = rcond,
        svdTol = svdTol,
        r2Threshold = r2Threshold,
        minimumLd = minimumLd,
        verbose = verbose
    )
    if (!is.null(genotypeMatrix)) {
        if (!is.null(ldMatrix)) {
            abort("Provide either ldMatrix or genotypeMatrix, not both.")
        }
        return(.raissGenotypePath(refPanel, knownZscores, genotypeMatrix, p))
    }
    if (is.null(ldMatrix)) {
        abort("Provide either ldMatrix or genotypeMatrix.")
    }
    isSingleMatrixCase <- is.matrix(ldMatrix) ||
        (is.list(ldMatrix) &&
            !is.null(ldMatrix$ldMatrices) &&
            length(ldMatrix$ldMatrices) == 1)
    if (isSingleMatrixCase) {
        return(.raissSingleLdPath(refPanel, knownZscores, ldMatrix, p))
    }
    .raissLdBlocksPath(refPanel, knownZscores, ldMatrix, p)
}

#' @param zt Vector of known z scores.
#' @param sigT Matrix of known linkage disequilibrium (LD) correlation.
#' @param sigIT Correlation matrix with rows corresponding to unknown SNPs (to
#'   impute) and columns to known SNPs.
#' @param lamb Regularization term added to the diagonal of the sigT matrix.
#' @param rcond Threshold for filtering eigenvalues in the pseudo-inverse
#'   computation.
#' @param batch Boolean indicating whether batch processing is used.
#'
#' @return A list containing the variance 'var', estimation 'mu', LD score
#'   'raissLdScore', condition number 'conditionNumber', and correctness of
#'   inversion 'correctInversion'.
#' @noRd
raissModel <- function(
    zt,
    sigT,
    sigIT,
    lamb = 0.01,
    rcond = 0.01,
    batch = TRUE,
    reportConditionNumber = FALSE
) {
    sigTInv <- invertMatRecursive(sigT, lamb, rcond)
    if (!is.numeric(zt) || !is.numeric(sigT) || !is.numeric(sigIT)) {
        abort("zt, sigT, and sigIT must be numeric.")
    }
    if (batch) {
        conditionNumber <- if (reportConditionNumber) {
            rep(kappa(sigT, exact = TRUE, norm = "2"), nrow(sigIT))
        } else {
            NA
        }
        correctInversion <- rep(checkInversion(sigT, sigTInv), nrow(sigIT))
    } else {
        conditionNumber <- if (reportConditionNumber) {
            kappa(sigT, exact = TRUE, norm = "2")
        } else {
            NA
        }
        correctInversion <- checkInversion(sigT, sigTInv)
    }

    varRaissLdScore <- computeVar(sigIT, sigTInv, lamb, batch)
    var <- varRaissLdScore$var
    raissLdScore <- varRaissLdScore$raissLdScore

    mu <- computeMu(sigIT, sigTInv, zt)
    varNorm <- varInBoundaries(var, lamb)

    R2 <- ((1 + lamb) - varNorm)
    mu <- mu / sqrt(R2)

    return(list(
        var = varNorm,
        mu = mu,
        raissLdScore = raissLdScore,
        conditionNumber = conditionNumber,
        correctInversion = correctInversion
    ))
}

#' @param imp is the output of raissModel()
#' @param refPanel is a data frame with columns 'chrom', 'pos', 'variant_id',
#'   'ref', and 'alt'.
#' @noRd
formatRaissDf <- function(imp, refPanel, unknowns) {
    # z / Var / raissLdScore come back from the RAISS matrix algebra as
    # 1-column matrices (e.g. imp$mu <- sigIT %*% ...). Flatten them to plain
    # numeric vectors EXPLICITLY with as.numeric() rather than leaning on
    # data.frame()'s implicit matrix-to-vector coercion (tibble would keep them
    # as matrix columns and break the downstream mergeRaissDf if_else). Extract
    # refPanel columns via [[ ]] + row index so we get a vector whatever class
    # refPanel is, instead of relying on the [, "col"] single-column drop.
    resultDf <- tibble(
        chrom = refPanel[["chrom"]][unknowns],
        pos = refPanel[["pos"]][unknowns],
        variant_id = refPanel[["variant_id"]][unknowns],
        A1 = refPanel[["A1"]][unknowns],
        A2 = refPanel[["A2"]][unknowns],
        z = as.numeric(imp$mu),
        Var = as.numeric(imp$var),
        raissLdScore = as.numeric(imp$raissLdScore),
        conditionNumber = imp$conditionNumber,
        correctInversion = imp$correctInversion
    )

    # Specify the column order
    columnOrder <- c(
        "chrom",
        "pos",
        "variant_id",
        "A1",
        "A2",
        "z",
        "Var",
        "raissLdScore",
        "conditionNumber",
        "correctInversion"
    )

    # Reorder the columns
    resultDf <- select(resultDf, all_of(columnOrder))
    return(resultDf)
}

#' @importFrom dplyr full_join
#' @noRd
mergeRaissDf <- function(raissDf, knownZscores) {
    # Full outer join keeps every variant from both frames (imputed + known).
    mergedDf <- full_join(
        raissDf,
        knownZscores,
        by = c("chrom", "pos", "variant_id", "A1", "A2")
    )

    # Identify rows that came from knownZscores
    fromKnown <- !is.na(mergedDf$z.y) & is.na(mergedDf$z.x)

    # Set Var to -1 and raissLdScore to Inf for these rows
    mergedDf$Var[fromKnown] <- -1
    mergedDf$raissLdScore[fromKnown] <- Inf

    # If there are overlapping columns (e.g., z.x and z.y), resolve them For
    # example, use z from knownZscores where available, otherwise use z from
    # raissDf
    mergedDf$z <- if_else(fromKnown, mergedDf$z.y, mergedDf$z.x)

    # Remove the extra columns produced by the join (z.x, z.y).
    mergedDf <- select(mergedDf, -all_of(c("z.x", "z.y")))
    mergedDf <- arrange(mergedDf, .data$pos)
    # assign imputed variants beta, se as NA to avoid confusion, since they are
    # not imputed. beta/se are optional (knownZscores may omit them), so guard
    # on column presence explicitly rather than relying on a data.frame
    # silently creating an all-NA column on `$col[mask] <- NA`.
    if (is_in("beta", colnames(mergedDf))) {
        mergedDf$beta[mergedDf$Var == -1] <- NA
    }
    if (is_in("se", colnames(mergedDf))) {
        mergedDf$se[mergedDf$Var == -1] <- NA
    }
    return(mergedDf)
}

# Format one aligned "label: value" report line (label left-padded to
# maxLabelLength).
# @noRd
.formatRaissLine <- function(label, value, maxLabelLength) {
    sprintf("%-*s %d", maxLabelLength, str_c(label, ":"), value)
}

# Print the RAISS imputation filter report (counts from pre/post frames).
.filterRaissReport <- function(
    zscoresNofilter,
    zscores,
    r2Threshold,
    minimumLd
) {
    r2 <- zscoresNofilter$raissR2
    counts <- c(
        nrow(zscoresNofilter),
        sum(r2 == 2.0, na.rm = TRUE),
        sum(r2 != 2.0, na.rm = TRUE),
        sum(zscoresNofilter$raissLdScore < minimumLd, na.rm = TRUE),
        sum(r2 < r2Threshold, na.rm = TRUE),
        nrow(zscores)
    )
    labels <- c(
        "Variants before filter",
        "Non-imputed variants",
        "Imputed variants",
        "Variants filtered because of low LD score",
        "Variants filtered because of low R2",
        "Remaining variants after filter"
    )
    maxLabelLength <- max(str_length(str_c(labels, ":")))
    inform("IMPUTATION REPORT\n")
    for (i in seq_along(labels)) {
        inform(.formatRaissLine(labels[i], counts[i], maxLabelLength))
    }
}

filterRaissOutput <- function(
    zscores,
    r2Threshold = 0.6,
    minimumLd = 5,
    verbose = TRUE
) {
    zscores <- select(
        zscores,
        all_of(c(
            "chrom",
            "pos",
            "variant_id",
            "A1",
            "A2",
            "z",
            "Var",
            "raissLdScore"
        ))
    )
    zscores$raissR2 <- 1 - zscores$Var
    zscoresNofilter <- zscores
    zscores <- filter(
        zscores,
        .data$raissR2 > r2Threshold & .data$raissLdScore >= minimumLd
    )
    if (verbose) {
        .filterRaissReport(zscoresNofilter, zscores, r2Threshold, minimumLd)
    }
    list(zscoresNofilter = zscoresNofilter, zscores = zscores)
}

computeMu <- function(sigIT, sigTInv, zt) {
    return(sigIT %*% (sigTInv %*% zt))
}

computeVar <- function(sigIT, sigTInv, lamb, batch = TRUE) {
    if (batch) {
        var <- (1 + lamb) - rowSums((sigIT %*% sigTInv) * sigIT)
        raissLdScore <- rowSums(sigIT^2)
    } else {
        var <- (1 + lamb) - (sigIT %*% (sigTInv %*% t(sigIT)))
        raissLdScore <- sum(sigIT^2)
    }
    return(list(var = var, raissLdScore = raissLdScore))
}

checkInversion <- function(sigT, sigTInv) {
    return(all.equal(sigT, sigT %*% (sigTInv %*% sigT), tolerance = 1e-5))
}

varInBoundaries <- function(var, lamb) {
    var[var < 0] <- 0
    var[var > (0.99999 + lamb)] <- 1
    return(var)
}

invertMat <- function(mat, lamb, rcond) {
    tryCatch(
        {
            # Modify the diagonal elements of mat
            diag(mat) <- 1 + lamb
            # Compute the pseudo-inverse
            matInv <- ginv(mat, tol = rcond)
            return(matInv)
        },
        error = function(e) {
            # Second attempt with updated lamb and rcond in case of an error
            diag(mat) <- 1 + lamb * 1.1
            matInv <- ginv(mat, tol = rcond * 1.1)
            return(matInv)
        }
    )
}

invertMatRecursive <- function(mat, lamb, rcond) {
    tryCatch(
        {
            # Modify the diagonal elements of mat
            diag(mat) <- 1 + lamb
            # Compute the pseudo-inverse
            matInv <- ginv(mat, tol = rcond)
            return(matInv)
        },
        error = function(e) {
            # Recursive call with updated lamb and rcond in case of an error
            invertMat(mat, lamb * 1.1, rcond * 1.1)
        }
    )
}

invertMatEigen <- function(mat, tol = 1e-3) {
    eigenMat <- eigen(mat)
    L <- which(cumsum(eigenMat$values) / sum(eigenMat$values) > 1 - tol)[1]
    if (is.na(L)) {
        # all eigen values are extremely small
        msg <- glue(
            "Cannot invert the input matrix because all its eigen ",
            "values are negative or close to zero"
        )
        abort(msg)
    }
    matInv <- eigenMat$vectors[, seq_len(L)] %*%
        diag(1 / eigenMat$values[seq_len(L)]) %*%
        t(eigenMat$vectors[, seq_len(L)])

    return(matInv)
}


# =============================================================================
# Top-level summaryStatsQc() pipeline + helpers
# =============================================================================

#' Detect LD-Summary Statistic Mismatches
#'
#' Unified wrapper for detecting outlier variants due to LD-summary statistic
#' mismatches. Dispatches to either \code{\link{dentistSingleWindow}} or
#' \code{\link{slalom}} based on the \code{method} argument.
#'
#' @param zScore Numeric vector of z-scores.
#' @param R Square LD correlation matrix. Provide either \code{R} or \code{X}.
#' @param X Genotype matrix (samples x SNPs). If provided, LD is computed via
#'   \code{computeLd} and \code{nSample} defaults to \code{nrow(X)}.
#' @param nSample Number of samples in the LD reference panel. Required when
#'   \code{R} is provided and \code{method = "dentist"}; inferred from \code{X}
#'   when \code{X} is provided.
#' @param method Character string specifying the QC method: \code{"slalom"}
#'   (default) or \code{"dentist"}.
#' @param ldMethod Character string specifying the LD computation method when
#'   \code{X} is provided. One of \code{"sample"} (default),
#'   \code{"population"}, or \code{"gcta"}. Ignored when \code{R} is provided
#'   directly.
#' @param ... Additional arguments passed to the underlying QC method
#'   (\code{\link{dentistSingleWindow}} or \code{\link{slalom}}).
#'
#' @return A data frame with at least a logical \code{outlier} column indicating
#'   which variants are identified as outliers. The remaining columns depend on
#'   the method used.
#'
#' @seealso \code{\link{dentistSingleWindow}}, \code{\link{slalom}},
#'   \code{\link{summaryStatsQc}}
#' @importFrom dplyr mutate row_number filter pull
#' @examples
#' data(eqtlRegionExample)
#' R <- cor(eqtlRegionExample$X[, 1:20])
#' ldMismatchQc(zScore = rnorm(20), R = R, nSample = 415)
#' @export
ldMismatchQc <- function(
    zScore,
    R = NULL,
    X = NULL,
    nSample = NULL,
    method = c("slalom", "dentist"),
    ldMethod = "sample",
    ...
) {
    method <- arg_match(method)
    if (method == "dentist") {
        qcResults <- dentistSingleWindow(
            zScore,
            R = R,
            X = X,
            nSample = nSample,
            ldMethod = ldMethod,
            ...
        )
        return(qcResults)
    } else {
        qcResults <- slalom(zScore, R = R, X = X, ldMethod = ldMethod, ...)
        # Standardize output: slalom uses "outliers", rename to "outlier" for
        # consistency
        result <- qcResults$data
        if (
            is_in("outliers", colnames(result)) &&
                !is_in("outlier", colnames(result))
        ) {
            result <- rename(result, outlier = "outliers")
        }
        return(result)
    }
}

.resolveZMismatchQc <- function(zMismatchQc) {
    if (is.null(zMismatchQc)) {
        return("none")
    }
    arg_match(zMismatchQc, c("none", "slalom", "dentist"))
}

#' Effective sample size for a case/control study
#'
#' Computes the effective sample size \code{N_eff = 4 / (1/nCase + 1/nControl) =
#' 4 * nCase * nControl / (nCase + nControl)} for case/control GWAS. Balanced
#' studies (\code{nCase == nControl}) recover the total \code{nCase + nControl};
#' imbalanced studies give a smaller value, which is the statistically correct
#' sample size for the RSS likelihood, residual-variance estimation, kriging,
#' and the N-cutoff filter. Vectorized over \code{nCase} / \code{nControl};
#' entries where either count is \code{NA} or \code{<= 0} return
#' \code{NA_real_}.
#'
#' @param nCase Numeric vector of case counts.
#' @param nControl Numeric vector of control counts.
#' @return Numeric vector of effective sample sizes (\code{NA_real_} where a
#'   count is missing or non-positive).
#' @references Prive et al., "Identifying and correcting for misspecifications
#'   in GWAS summary statistics and polygenic scores", HGG Advances 2022.
#' @examples
#' nCase <- 1000
#' nControl <- 2000
#' effectiveN(nCase = nCase, nControl = nControl)
#' @export
effectiveN <- function(nCase, nControl) {
    nCase <- as.numeric(nCase)
    nControl <- as.numeric(nControl)
    out <- 4 / (1 / nCase + 1 / nControl)
    out[is.na(nCase) | is.na(nControl) | nCase <= 0 | nControl <= 0] <- NA_real_
    out
}

# Require a susieR that provides the kriging RSS diagnostic.
.krigingCheckSusie <- function() {
    if (
        !requireNamespace("susieR", quietly = TRUE) ||
            !all(
                is_in(
                    c("estimate_s_rss", "kriging_rss"),
                    getNamespaceExports("susieR")
                )
            )
    ) {
        # nocov start
        msg <- glue(
            "krigingOutlierQc requires a susieR that provides ",
            "estimate_s_rss() and kriging_rss(); the installed susieR does ",
            "not. Install a susieR with the kriging RSS diagnostic, or ",
            "disable alleleFlipKriging."
        )
        abort(msg)
        # nocov end
    }
}

#' Kriging-style LD-consistency outlier QC
#'
#' Flags variants whose observed z-score is inconsistent with the value
#' predicted from its LD neighbours, using susieR's kriging diagnostic.
#' \code{susieR::kriging_rss()} computes the leave-one-out conditional
#' distribution of each \code{z_i} given the rest (with the LD-mismatch scale
#' \code{s} defaulting to \code{susieR::estimate_s_rss()}) and a per-variant
#' \code{logLR} for the allele-switch hypothesis. This helper reuses susieR's
#' own allele-switch rule --- \code{logLR > logLRThreshold & abs(z) >
#' zThreshold} (the same \code{logLR > 2 & |z| > 2} used in
#' \code{susie_rss_utils}) --- to flag variants whose sign should be flipped.
#' RSS-only helper, opt-in via \code{alleleFlipKriging}; never wired into
#' \code{alleleQc()} / \code{matchRefPanel()}. Requires a susieR that provides
#' \code{kriging_rss()} and \code{estimate_s_rss()}.
#'
#' @param zScore Numeric vector of harmonized z-scores.
#' @param R Square LD correlation matrix aligned to \code{zScore}.
#' @param n Sample size, forwarded to \code{susieR::kriging_rss()} (whose
#'   default \code{s} is \code{susieR::estimate_s_rss()}).
#' @param variantIds Optional variant IDs for the diagnostics table.
#' @param zThreshold Absolute-z cutoff for the allele-switch rule (default
#'   \code{2}, matching susieR).
#' @param logLRThreshold Log-likelihood-ratio cutoff for the allele-switch rule
#'   (default \code{2}, matching susieR).
#' @return A list with \code{flip} (logical vector; \code{TRUE} = allele switch,
#'   z-score should be sign-flipped) and \code{diagnostics} (data frame of
#'   per-variant \code{z}, \code{condmean}, \code{z_std_diff}, \code{logLR}, and
#'   the \code{flipped} flag).
#' @importFrom stats pnorm
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:20]
#' R <- cor(X)
#' krigingOutlierQc(
#'   zScore = rnorm(20), R = R, n = 415, variantIds = colnames(X))
#' @export
krigingOutlierQc <- function(
    zScore,
    R,
    n,
    variantIds = NULL,
    zThreshold = 2,
    logLRThreshold = 2
) {
    zScore <- as.numeric(zScore)
    m <- length(zScore)
    if (is.null(R) || !is.matrix(R) || nrow(R) != m || ncol(R) != m) {
        abort("krigingOutlierQc requires a square LD matrix aligned to zScore.")
    }
    if (missing(n) || length(n) != 1L || is.na(n) || !is.finite(n) || n <= 0) {
        abort("krigingOutlierQc requires a single positive sample size 'n'.")
    }
    .krigingCheckSusie()
    if (is.null(variantIds)) {
        variantIds <- rownames(R)
    }
    # susieR kriging RSS diagnostic (kriging_rss(z, R, n)): s defaults to
    # estimate_s_rss(); logLR matches susieR's allele-switch selection.
    cd <- susieR::kriging_rss(z = zScore, R = R, n = n)$conditional_dist
    condMean <- as.numeric(cd$condmean)
    zStdDiff <- as.numeric(cd$z_std_diff)
    logLR <- as.numeric(cd$logLR)
    # susieR's allele-switch rule (susie_rss_utils.R): logLR > 2 & |z| > 2.
    flip <- !is.na(logLR) &
        !is.na(zScore) &
        logLR > logLRThreshold &
        abs(zScore) > zThreshold
    list(
        flip = flip,
        diagnostics = tibble(
            variant_id = if (is.null(variantIds)) seq_len(m) else variantIds,
            z = zScore,
            condmean = condMean,
            z_std_diff = zStdDiff,
            logLR = logLR,
            flipped = flip
        )
    )
}

# =============================================================================
# summaryStatsQc -- SumStats-input QC pipeline (replaces the previous
# data.frame/LdData/QcResult-based summaryStatsQc and rssBasicQc).
# =============================================================================

# Convert one entry's GRanges into a flat tibble with the column shape
# harmonizeAlleles expects (lower-case chrom/pos plus the CapsCase mcols).
.entryGrangesToDf <- function(gr) {
    mc <- as.data.frame(S4Vectors::mcols(gr), stringsAsFactors = FALSE)
    out <- tibble(
        chrom = str_remove(
            as.character(GenomicRanges::seqnames(gr)),
            regex("^chr", ignore_case = TRUE)
        ),
        pos = GenomicRanges::start(gr)
    )
    bind_cols(out, mc)
}

# Build a refVariants data.frame (chrom, pos, A1, A2, variant_id) from the
# ldSketch GenotypeHandle's snpInfo so harmonizeAlleles can join by (chrom,
# pos).
.refVariantsFromSketch <- function(handle) {
    si <- getSnpInfo(handle)
    chr <- str_remove(as.character(si$CHR), regex("^chr", ignore_case = TRUE))
    data.frame(
        chrom = chr,
        pos = as.integer(si$BP),
        A1 = as.character(si$A1),
        A2 = as.character(si$A2),
        variant_id = as.character(si$SNP),
        stringsAsFactors = FALSE
    )
}

# Reassemble a harmonized data.frame into a GRanges with the SumStats mcol
# shape (SNP, A1, A2, Z, N, ... optional MAF/INFO/BETA/SE/P kept if present).
.dfToEntryGranges <- function(df) {
    # Short-circuit on empty input: `paste0("chr", character(0))` returns
    # "chr" (a length-1 vector), not character(0), so we cannot rely on the
    # paste/IRanges constructors to handle the zero-row case cleanly.
    chrRaw <- as.character(df$chrom)
    if (length(chrRaw) == 0L) {
        gr <- GenomicRanges::GRanges()
        return(gr)
    }
    chr <- withChrPrefix(chrRaw)
    gr <- GenomicRanges::GRanges(
        seqnames = chr,
        ranges = IRanges::IRanges(start = as.integer(df$pos), width = 1L)
    )
    if (is_in("variant_id", colnames(df)) && !is_in("SNP", colnames(df))) {
        df$SNP <- df$variant_id
    }
    baseCols <- c("SNP", "A1", "A2", "Z", "N")
    optCols <- c("MAF", "INFO", "BETA", "SE", "P", "N_CASE", "N_CONTROL")
    use <- intersect(c(baseCols, optCols), colnames(df))
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(select(df, all_of(use)))
    gr
}

# -----------------------------------------------------------------------------
# Shared entry-to-sumstat data.frame converter
# -----------------------------------------------------------------------------

# Internal: convert one sumstat-entry GRanges into a flat data.frame with
# (variant_id, chrom, pos, A1, A2) base columns plus a configurable set
# of stats columns (z, beta, se, N, maf). Shared by the four pipelines
# that walk sumstats GRanges (fineMappingPipeline, twasWeights,
# ctwasPipeline, colocboostPipeline).
#
# `require`   character vector of mcol names that MUST be present; errors
#             when any is missing. Use this for the strict callers
#             (e.g. `.fmExtractZn` needs SNP + Z + N to proceed).
# `derive`    when "zFromBetaSe" and `z` is absent but BETA + SE are
#             present, set z := BETA/SE. Default "none".
# `label`     error-message prefix for missing-`require` errors.
# `keepChrPrefix`  when TRUE keep the seqname as-is ("chr1"); when FALSE,
#                  strip any leading "chr" so callers that expect numeric
#                  chrom (ctwas, colocboost) see "1".
# Base (variant_id/chrom/pos/A1/A2) frame for one sumstat-entry GRanges.
.entryDfBase <- function(gr, mc, chr) {
    colOr <- function(name) {
        if (is_in(name, colnames(mc))) {
            as.character(mc[[name]])
        } else {
            rep(NA_character_, length(gr))
        }
    }
    tibble(
        variant_id = colOr("SNP"),
        chrom = chr,
        pos = unname(as.integer(GenomicRanges::start(gr))),
        A1 = colOr("A1"),
        A2 = colOr("A2")
    )
}

# Append the optional numeric stat columns present on the entry's mcols.
.entryDfAddStats <- function(df, mc) {
    statMap <- c(z = "Z", beta = "BETA", se = "SE", N = "N", maf = "MAF")
    for (out in names(statMap)) {
        src <- statMap[[out]]
        if (is_in(src, colnames(mc))) {
            df[[out]] <- as.numeric(mc[[src]])
        }
    }
    df
}

.entryToSumstatDf <- function(
    gr,
    require = character(0),
    derive = c("none", "zFromBetaSe"),
    label = "entry",
    keepChrPrefix = TRUE
) {
    derive <- arg_match(derive)
    mc <- S4Vectors::mcols(gr)
    for (col in require) {
        if (!is_in(col, colnames(mc))) {
            msg <- glue("{label}: entry has no {col} mcol.")
            abort(msg)
        }
    }
    chr <- as.character(GenomicRanges::seqnames(gr))
    if (!keepChrPrefix) {
        chr <- str_remove(chr, regex("^chr", ignore_case = TRUE))
    }
    df <- .entryDfBase(gr, mc, chr)
    df <- .entryDfAddStats(df, mc)
    if (
        derive == "zFromBetaSe" &&
            is.null(df[["z"]]) &&
            !is.null(df[["beta"]]) &&
            !is.null(df[["se"]])
    ) {
        df[["z"]] <- df[["beta"]] / df[["se"]]
    }
    df
}

# Derive BETA and SE columns from signed Z when the entry has only Z.
# Formula (Zhu et al. 2016 / RAISS):
#   se   = 1 / sqrt(2 * maf * (1 - maf) * (N + z^2))
#   beta = z * se
# Requires Z, MAF, and N to all be present in `df`. No-op if BETA and SE
# are already there, or if any required column is missing. Returns:
#   list(df = <data.frame>, audit = NULL | list(nDerived = <int>))
# .zToPvalue (two-tailed normal p-value from a signed Z) is defined once in
# pvalCombine.R and shared package-wide.

# Internal: thin SVD with numerical-stability filtering. Drops singular
# values below `tol * max(d)` and caps the retained rank at `maxRank`.
# Used by RAISS imputation (raissSingleMatrixFromX) to invert a panel
# genotype matrix safely under rank deficiency.
.safeSvd <- function(mat, tol = 1e-8, maxRank = NULL) {
    if (max(abs(mat)) == 0) {
        abort("Cannot compute SVD of an all-zero matrix.")
    }
    s <- svd(mat)
    d <- s$d
    keep <- if (tol > 0 && length(d) > 0) {
        out <- d / d[1] > tol
        if (!any(out)) {
            abort("All singular values are below the tolerance threshold.")
        }
        out
    } else {
        rep(TRUE, length(d))
    }
    if (!is.null(maxRank) && maxRank > 0) {
        nKeep <- min(sum(keep), maxRank)
        keepIdx <- which(keep)
        if (length(keepIdx) > nKeep) {
            keep[keepIdx[(nKeep + 1):length(keepIdx)]] <- FALSE
        }
    }
    list(
        u = s$u[, keep, drop = FALSE],
        d = d[keep],
        v = s$v[, keep, drop = FALSE]
    )
}

# Internal: identify LD-correlated duplicate variants. Walks the LD
# matrix left-to-right, marking each variant as either a unique anchor
# (dupBearer == -1) or a duplicate of an earlier anchor (dupBearer == k,
# the anchor's index). Returns filtered z / LD plus per-variant
# bookkeeping that DENTIST's addDupsBackDentist uses to splice the
# dropped variants back into the output.
.findDuplicateVariants <- function(z, ld, rThreshold) {
    p <- length(z)
    dupBearer <- rep(-1, p)
    corABS <- rep(0, p)
    sign <- rep(1, p)
    count <- 1L
    minValue <- 1
    for (i in seq_len(p - 1L)) {
        if (dupBearer[i] != -1) {
            next
        }
        idx <- (i + 1L):p
        corVec <- abs(ld[i, idx])
        dupIdx <- which(dupBearer[idx] == -1 & corVec > rThreshold)
        if (length(dupIdx) > 0) {
            j <- idx[dupIdx]
            sign[j] <- if_else(ld[i, j] < 0, -1, sign[j])
            corABS[j] <- corVec[dupIdx]
            dupBearer[j] <- count
        }
        minValue <- min(minValue, min(corVec))
        count <- count + 1L
    }
    filteredZ <- z[dupBearer == -1]
    filteredLD <- ld[dupBearer == -1, dupBearer == -1, drop = FALSE]
    list(
        filteredZ = filteredZ,
        filteredLD = filteredLD,
        dupBearer = dupBearer,
        corABS = corABS,
        sign = sign,
        minValue = minValue
    )
}

.deriveBetaSeFromZ <- function(df) {
    hasBeta <- is_in("BETA", colnames(df))
    hasSe <- is_in("SE", colnames(df))
    if (hasBeta && hasSe) {
        return(list(df = df, audit = NULL))
    }
    hasZ <- is_in("Z", colnames(df))
    hasMaf <- is_in("MAF", colnames(df))
    hasN <- is_in("N", colnames(df))
    if (!(hasZ && hasMaf && hasN)) {
        return(list(df = df, audit = NULL))
    }
    z <- as.numeric(df$Z)
    maf <- as.numeric(df$MAF)
    n <- as.numeric(df$N)
    bs <- .zToBetaSe(z, maf, n)
    se <- bs$se
    beta <- bs$beta
    if (!hasBeta) {
        df$BETA <- beta
    }
    if (!hasSe) {
        df$SE <- se
    }
    list(df = df, audit = list(nDerived = sum(!is.na(se))))
}

# Drop variants whose (chrom, pos) overlaps any user-supplied skipRegion.
# skipRegion may be a character vector of "chr:start-end" strings or a GRanges.
# Parse one 'chr:start-end' skip-region string into a 1-row data.frame.
.parseSkipRegionEntry <- function(s) {
    m <- str_match(s, "^([^:]+):([0-9]+)-([0-9]+)$")[1L, ]
    if (is.na(m[[1L]])) {
        msg <- glue("skipRegion entry must be 'chr:start-end'; got '{s}'")
        abort(msg)
    }
    tibble(
        chrom = str_remove(m[[2L]], regex("^chr", ignore_case = TRUE)),
        start = as.integer(m[[3L]]),
        end = as.integer(m[[4L]])
    )
}

# Normalise skipRegion (character vector or GRanges) to a chrom/start/end frame.
.parseSkipRegion <- function(skipRegion) {
    if (is.character(skipRegion)) {
        return(bind_rows(map(skipRegion, .parseSkipRegionEntry)))
    }
    if (methods::is(skipRegion, "GRanges")) {
        return(tibble(
            chrom = str_remove(
                as.character(GenomicRanges::seqnames(skipRegion)),
                regex("^chr", ignore_case = TRUE)
            ),
            start = GenomicRanges::start(skipRegion),
            end = GenomicRanges::end(skipRegion)
        ))
    }
    msg <- glue(
        "skipRegion must be a character vector of 'chr:start-end' ",
        "strings or a GRanges."
    )
    abort(msg)
}

.applySkipRegion <- function(df, skipRegion) {
    if (is.null(skipRegion) || length(skipRegion) == 0L) {
        return(df)
    }
    parsed <- .parseSkipRegion(skipRegion)
    dropMask <- rep(FALSE, nrow(df))
    dfChr <- str_remove(
        as.character(df$chrom),
        regex("^chr", ignore_case = TRUE)
    )
    for (i in seq_len(nrow(parsed))) {
        dropMask <- dropMask |
            (dfChr == parsed$chrom[i] &
                df$pos >= parsed$start[i] &
                df$pos <= parsed$end[i])
    }
    filter(df, !dropMask)
}

# Apply the panel-vs-sumstats allele harmonization using the slim
# harmonizeAlleles against the ldSketch's variant info. Threads the
# variant-level filters (indels, strand-ambiguous, duplicates) through
# so the LD-panel-anchored pass handles them in a single sweep.
.matchAgainstSketch <- function(
    df,
    ldSketch,
    matchMinProp,
    removeIndels = FALSE,
    removeStrandAmbiguous = TRUE,
    removeDups = TRUE
) {
    refVariants <- .refVariantsFromSketch(ldSketch)
    flipCandidates <- c("Z", "BETA")
    colToFlip <- intersect(flipCandidates, colnames(df))
    if (length(colToFlip) == 0L) {
        msg <- glue(
            "summaryStatsQc: input entry must contain at least one of ",
            "Z or BETA before panel harmonization."
        )
        abort(msg)
    }
    colToComplement <- intersect("MAF", colnames(df))
    if (!is_in("A1", colnames(df)) || !is_in("A2", colnames(df))) {
        abort("summaryStatsQc: input entry must contain A1 and A2 columns.")
    }
    res <- harmonizeAlleles(
        targetData = df,
        refVariants = refVariants,
        colToFlip = colToFlip,
        colToComplement = colToComplement,
        matchMinProp = matchMinProp,
        removeUnmatched = TRUE,
        removeIndels = removeIndels,
        removeStrandAmbiguous = removeStrandAmbiguous,
        removeDups = removeDups
    )
    out <- res$harmonizedData
    if (!is_in("chrom", colnames(out)) && is_in("chr", colnames(out))) {
        colnames(out)[colnames(out) == "chr"] <- "chrom"
    }
    attr(out, "qcCounts") <- attr(res, "qcCounts")
    out
}

# Variant-content filters (MAF / INFO / N). Pure data-frame column
# filters; no Bioconductor genome packages needed.
#
# mafCutoff:  drop rows where MAF (or FRQ) < mafCutoff. Requires either
#             column when mafCutoff > 0; errors if neither is present.
# infoCutoff: drop rows where INFO < infoCutoff. Requires INFO column
#             when infoCutoff > 0.
# nCutoff:    drop rows whose N is more than nCutoff median-absolute-
#             deviations from the median (a 5-MAD-from-median cap on
#             per-variant N). Set nCutoff = 0 to disable. Rows with NA N
#             are always dropped.
# MAF/FRQ frequency filter (frequency normalised to minor-allele frequency).
.cfMaf <- function(df, mafCutoff) {
    if (mafCutoff <= 0) {
        return(list(df = df, dropped = NULL))
    }
    mafCol <- intersect(c("MAF", "FRQ"), colnames(df))[1L]
    if (is.na(mafCol)) {
        msg <- glue(
            ".applyContentFilters: mafCutoff > 0 requires a MAF or FRQ ",
            "column."
        )
        abort(msg)
    }
    before <- nrow(df)
    mafVals <- as.numeric(df[[mafCol]])
    # Normalise effect-allele frequency to MAF: take min(af, 1-af).
    mafVals <- pmin(mafVals, 1 - mafVals, na.rm = FALSE)
    df <- filter(df, !is.na(mafVals) & mafVals >= mafCutoff)
    list(df = df, dropped = before - nrow(df))
}

# INFO (imputation-quality) filter.
.cfInfo <- function(df, infoCutoff) {
    if (infoCutoff <= 0) {
        return(list(df = df, dropped = NULL))
    }
    if (!is_in("INFO", colnames(df))) {
        abort(".applyContentFilters: infoCutoff > 0 requires an INFO column.")
    }
    before <- nrow(df)
    infoVals <- as.numeric(df$INFO)
    df <- filter(df, !is.na(infoVals) & infoVals >= infoCutoff)
    list(df = df, dropped = before - nrow(df))
}

# Per-variant N outlier filter (MAD-z on the sample-size column).
.cfN <- function(df, nCutoff) {
    if (!(nCutoff > 0 && is_in("N", colnames(df)) && nrow(df) > 0L)) {
        return(list(df = df, dropped = NULL))
    }
    nVals <- as.numeric(df$N)
    before <- nrow(df)
    if (any(is.na(nVals))) {
        df <- filter(df, !is.na(nVals))
        nVals <- nVals[!is.na(nVals)]
    }
    if (length(nVals) > 0L) {
        medN <- stats::median(nVals)
        madN <- stats::mad(nVals, constant = 1)
        if (madN > 0) {
            zN <- abs(nVals - medN) / madN
            df <- filter(df, zN <= nCutoff)
        }
    }
    list(df = df, dropped = before - nrow(df))
}

.applyContentFilters <- function(
    df,
    mafCutoff = 0,
    infoCutoff = 0,
    nCutoff = 5
) {
    audit <- list()
    r <- .cfMaf(df, mafCutoff)
    df <- r$df
    if (!is.null(r$dropped)) {
        audit$mafDropped <- r$dropped
    }
    r <- .cfInfo(df, infoCutoff)
    df <- r$df
    if (!is.null(r$dropped)) {
        audit$infoDropped <- r$dropped
    }
    r <- .cfN(df, nCutoff)
    df <- r$df
    if (!is.null(r$dropped)) {
        audit$nDropped <- r$dropped
    }
    list(df = df, audit = audit)
}

# Per-row variant sanity / hygiene checks ported from MungeSumstats's
# check_*.R series but rewritten as pure data.frame operations with no
# genome / dbSNP dependency. Each step is gated by its own flag so a
# caller can disable any single check.
#
# Steps (in order; each contributes a count to audit):
#   - coerceNumeric: cast signed columns to numeric (catches stray "0.5"
#       strings). NA-introducing coercions are counted.
#   - normalizeChr:   strip "chr"/"ch" prefix, uppercase X/Y/MT, map
#       23->X, 24->Y, M->MT. Optional dropNonstandardChr removes rows
#       whose CHR is outside 1..22, X, Y, MT after normalization.
#   - dropMissData:   drop rows with NA in any vital column (chrom, pos,
#       A1, A2, and at least one of Z / BETA).
#   - dropPOutOfRange: drop rows where P < 0 or P > 1 (corrupt p-values).
#       Only fires when a P column is present.
#   - clampSmallP:    floor 0 <= P <= smallPFloor to smallPFloor so
#       -log10(P) stays finite downstream.
#   - dropZeroEffect: drop rows where any effect column is exactly 0
#       (BETA / LOG_ODDS / SIGNED_SUMSTAT) or OR is exactly 1. MungeSumstats
#       treats these as degenerate / artefactual.
#   - dropNonpositiveSe: drop rows where SE <= 0.
# --- .applySanityChecks per-check helpers (each guards its own flag) --------

# Coerce known numeric columns; count NAs newly introduced by coercion.
.scCoerceNumeric <- function(df, coerceNumeric) {
    if (!coerceNumeric) {
        return(list(df = df, audit = list()))
    }
    numericCols <- intersect(
        c(
            "Z",
            "BETA",
            "SE",
            "OR",
            "LOG_ODDS",
            "SIGNED_SUMSTAT",
            "P",
            "MAF",
            "FRQ",
            "INFO",
            "N"
        ),
        colnames(df)
    )
    naIntroduced <- 0L
    for (col in numericCols) {
        orig <- df[[col]]
        if (is.numeric(orig)) {
            next
        }
        coerced <- suppressWarnings(as.numeric(orig))
        naIntroduced <- naIntroduced + sum(is.na(coerced) & !is.na(orig))
        df[[col]] <- coerced
    }
    audit <- list()
    if (naIntroduced > 0L) {
        audit$nonNumericCoerced <- naIntroduced
    }
    list(df = df, audit = audit)
}

# Normalize chromosome labels; optionally drop non-standard chromosomes.
.scNormalizeChr <- function(df, normalizeChr, dropNonstandardChr) {
    if (!normalizeChr || !is_in("chrom", colnames(df))) {
        return(list(df = df, audit = list()))
    }
    chr <- as.character(df$chrom)
    chr <- str_remove(chr, regex("^chr", ignore_case = TRUE))
    chr <- str_remove(chr, regex("^ch", ignore_case = TRUE))
    chr <- str_to_upper(chr)
    chr[chr == "23"] <- "X"
    chr[chr == "24"] <- "Y"
    chr[chr == "M"] <- "MT"
    df$chrom <- chr
    audit <- list()
    if (dropNonstandardChr) {
        before <- nrow(df)
        standardChrs <- c(as.character(seq_len(22)), "X", "Y", "MT")
        df <- filter(df, is_in(chr, standardChrs))
        dropped <- before - nrow(df)
        if (dropped > 0L) audit$nonstandardChrDropped <- dropped
    }
    list(df = df, audit = audit)
}

# Drop rows missing any vital column (chrom/pos/A1/A2 + first signed stat).
.scDropMissData <- function(df, dropMissData) {
    if (!dropMissData || nrow(df) == 0L) {
        return(list(df = df, audit = list()))
    }
    vital <- intersect(c("chrom", "pos", "A1", "A2"), colnames(df))
    signedCol <- intersect(c("Z", "BETA"), colnames(df))[1L]
    if (!is.na(signedCol)) {
        vital <- c(vital, signedCol)
    }
    audit <- list()
    if (length(vital) > 0L) {
        before <- nrow(df)
        bad <- reduce(map(vital, .scColIsNa, df = df), `|`)
        if (any(bad)) {
            df <- filter(df, !bad)
        }
        dropped <- before - nrow(df)
        if (dropped > 0L) audit$missDataDropped <- dropped
    }
    list(df = df, audit = audit)
}

# Drop rows whose P is outside [0, 1].
.scDropPOutOfRange <- function(df, dropPOutOfRange) {
    if (!dropPOutOfRange || !is_in("P", colnames(df)) || nrow(df) == 0L) {
        return(list(df = df, audit = list()))
    }
    before <- nrow(df)
    p <- as.numeric(df$P)
    bad <- !is.na(p) & (p < 0 | p > 1)
    if (any(bad)) {
        df <- filter(df, !bad)
    }
    dropped <- before - nrow(df)
    audit <- list()
    if (dropped > 0L) {
        audit$pOutOfRangeDropped <- dropped
    }
    list(df = df, audit = audit)
}

# Clamp tiny P-values up to the floor.
.scClampSmallP <- function(df, clampSmallP, smallPFloor) {
    if (!clampSmallP || !is_in("P", colnames(df)) || nrow(df) == 0L) {
        return(list(df = df, audit = list()))
    }
    p <- as.numeric(df$P)
    smallMask <- !is.na(p) & p >= 0 & p < smallPFloor
    nClamped <- sum(smallMask)
    audit <- list()
    if (nClamped > 0L) {
        df$P[smallMask] <- smallPFloor
        audit$smallPClamped <- nClamped
    }
    list(df = df, audit = audit)
}

# Drop rows whose effect equals the null sentinel (0, or 1 for OR).
.scDropZeroEffect <- function(df, dropZeroEffect) {
    if (!dropZeroEffect || nrow(df) == 0L) {
        return(list(df = df, audit = list()))
    }
    effectCols <- intersect(
        c("BETA", "LOG_ODDS", "SIGNED_SUMSTAT", "OR"),
        colnames(df)
    )
    audit <- list()
    if (length(effectCols) > 0L) {
        before <- nrow(df)
        badMask <- rep(FALSE, nrow(df))
        for (col in effectCols) {
            vals <- as.numeric(df[[col]])
            sentinel <- if (col == "OR") 1 else 0
            badMask <- badMask | (!is.na(vals) & vals == sentinel)
        }
        if (any(badMask)) {
            df <- filter(df, !badMask)
        }
        dropped <- before - nrow(df)
        if (dropped > 0L) audit$zeroEffectDropped <- dropped
    }
    list(df = df, audit = audit)
}

# Drop rows with non-positive standard error.
.scDropNonpositiveSe <- function(df, dropNonpositiveSe) {
    if (!dropNonpositiveSe || !is_in("SE", colnames(df)) || nrow(df) == 0L) {
        return(list(df = df, audit = list()))
    }
    before <- nrow(df)
    se <- as.numeric(df$SE)
    bad <- !is.na(se) & se <= 0
    if (any(bad)) {
        df <- filter(df, !bad)
    }
    dropped <- before - nrow(df)
    audit <- list()
    if (dropped > 0L) {
        audit$nonpositiveSeDropped <- dropped
    }
    list(df = df, audit = audit)
}

# Per-row sanity checks: sequence the guarded checks, accumulating the audit.
.applySanityChecks <- function(
    df,
    coerceNumeric = TRUE,
    normalizeChr = TRUE,
    dropNonstandardChr = TRUE,
    dropMissData = TRUE,
    dropPOutOfRange = TRUE,
    clampSmallP = TRUE,
    smallPFloor = 5e-324,
    dropZeroEffect = TRUE,
    dropNonpositiveSe = TRUE
) {
    audit <- list()
    if (nrow(df) == 0L) {
        return(list(df = df, audit = audit))
    }
    r <- .scCoerceNumeric(df, coerceNumeric)
    df <- r$df
    audit <- modifyList(audit, r$audit)
    r <- .scNormalizeChr(df, normalizeChr, dropNonstandardChr)
    df <- r$df
    audit <- modifyList(audit, r$audit)
    r <- .scDropMissData(df, dropMissData)
    df <- r$df
    audit <- modifyList(audit, r$audit)
    r <- .scDropPOutOfRange(df, dropPOutOfRange)
    df <- r$df
    audit <- modifyList(audit, r$audit)
    r <- .scClampSmallP(df, clampSmallP, smallPFloor)
    df <- r$df
    audit <- modifyList(audit, r$audit)
    r <- .scDropZeroEffect(df, dropZeroEffect)
    df <- r$df
    audit <- modifyList(audit, r$audit)
    r <- .scDropNonpositiveSe(df, dropNonpositiveSe)
    df <- r$df
    audit <- modifyList(audit, r$audit)
    list(df = df, audit = audit)
}

# Apply ldMismatchQc (SLALOM/DENTIST) against the LD sketch. Returns the
# filtered df, outlier count, and the full per-variant diagnostics table
# (the data.frame returned by ldMismatchQc(), prepended with a
# variant_id column for downstream joins). Callers record `diagnostics`
# in the entry's qcInfo audit so the per-variant detail is available
# for plotting / postprocessing instead of just the outlier count.
.applyLdMismatchQcToEntry <- function(df, ldSketch, method) {
    variantIds <- df$SNP
    if (is.null(variantIds) || any(is.na(variantIds))) {
        abort("summaryStatsQc: ldMismatchQc requires SNP column on the entry.")
    }
    # Panel LD for the entry variants via the shared LD-from-sketch helper
    # (tuple match with chr-prefix tolerance, strand-ambiguous variants kept;
    # errors if any variant is absent from the panel).
    R <- .ldFromSketch(
        ldSketch,
        variantIds,
        label = "summaryStatsQc: zMismatchQc"
    )
    qc <- ldMismatchQc(
        zScore = df$Z,
        R = R,
        nSample = getNSamples(ldSketch),
        method = method
    )
    # slalom / dentist can leave NA in the outlier column when their
    # per-variant statistic is undefined (e.g. a degenerate dentist
    # chisq for variants effectively orthogonal to the lead). Treat NA as
    # "no evidence of being an outlier" (conservative: keep the variant)
    # so the downstream df / sum() / IRanges construction stay finite.
    outlierFlags <- qc$outlier
    outlierFlags[is.na(outlierFlags)] <- FALSE
    # Attach the variant_id column so the diagnostics data.frame stays
    # self-describing once it's separated from the input df.
    diagnostics <- if (is.data.frame(qc)) {
        cbind(
            variant_id = as.character(variantIds),
            qc,
            stringsAsFactors = FALSE
        )
    } else {
        NULL
    }
    list(
        df = filter(df, !outlierFlags),
        outliers = sum(outlierFlags),
        diagnostics = diagnostics
    )
}

# -----------------------------------------------------------------------------
# Signal screen: skip an entry/block with no strong signal by a chosen metric.
# -----------------------------------------------------------------------------
# The screen is driven by ONE metric at a time (enforced by
# .resolveScreenMetric): pip : max single-effect PIP (susie_ser $pip); cutoff<0
# => 3/nVar absZ : max |Z| (no model fit) logBf : max per-variant single-effect
# logBF (susie_ser $lbf_variable) bf : same evidence on the raw BF scale
# (compare maxlogBF > log(cutoff)) A "screen spec" flowing through the pipelines
# is EITHER a legacy PIP cutoff (numeric scalar, 0 = off -- the historical
# `pipCutoffToSkip`) OR a resolved screen object `list(metric, cutoff)`.
# .asScreen() canonicalizes either into `list(metric, cutoff)` or NULL (no
# screen); the pipeline channels stay untyped so only the screen primitives need
# to interpret the spec.

# Canonicalize a screen spec into list(metric, cutoff) or NULL (no screen).
.asScreen <- function(spec) {
    if (is.null(spec)) {
        return(NULL)
    }
    if (is.list(spec)) {
        # already a screen object
        if (
            is.null(spec$metric) ||
                is.null(spec$cutoff) ||
                is.na(spec$cutoff) ||
                spec$cutoff == 0
        ) {
            return(NULL)
        }
        return(spec)
    }
    # Legacy numeric: a PIP cutoff. Only a non-zero scalar activates it (a
    # non-scalar / NA / 0 means no screen), matching the historical behaviour.
    if (length(spec) != 1L || is.na(spec) || spec == 0) {
        return(NULL)
    }
    list(metric = "pip", cutoff = as.numeric(spec))
}

# Turn the four public cutoff arguments into a single screen object (or NULL).
# Enforces one-metric-at-a-time and rejects meaningless negative cutoffs. The
# pip metric keeps its `< 0 => 3 / nVariants` adaptive convention (resolved at
# screen time); absZ / bf must be > 0; logBf may be any non-zero value.
.resolveScreenMetric <- function(
    pipCutoffToSkip = 0,
    absZCutoffToSkip = 0,
    bfCutoffToSkip = 0,
    logBfCutoffToSkip = 0
) {
    scalar <- function(x) {
        if (is.null(x) || length(x) != 1L || is.na(x)) 0 else as.numeric(x)
    }
    cuts <- c(
        pip = scalar(pipCutoffToSkip),
        absZ = scalar(absZCutoffToSkip),
        bf = scalar(bfCutoffToSkip),
        logBf = scalar(logBfCutoffToSkip)
    )
    on <- cuts[cuts != 0]
    if (length(on) == 0L) {
        return(NULL)
    }
    if (length(on) > 1L) {
        nonZero <- str_flatten(sprintf("%s=%g", names(on), on), collapse = ", ")
        msg <- glue(
            "Only one signal screen may be enabled at a time, but these ",
            "are non-zero: {nonZero}. Set all but one of pipCutoffToSkip / ",
            "absZCutoffToSkip / bfCutoffToSkip / logBfCutoffToSkip to 0."
        )
        abort(msg)
    }
    metric <- names(on)
    cutoff <- unname(on[[1L]])
    if (metric == "absZ" && cutoff < 0) {
        abort("absZCutoffToSkip must be > 0 (it screens on max|Z|).")
    }
    if (metric == "bf" && cutoff < 0) {
        abort("bfCutoffToSkip must be > 0 (Bayes factors are positive).")
    }
    list(metric = metric, cutoff = cutoff)
}

# Decide whether the z-scores of one entry clear the chosen screen. Returns
# list(ok = logical, reason = character). susie_ser is fit at most once, and
# only for the metrics that need it (absZ stays model-free).
.entryScreenPass <- function(z, n, nVar, scr) {
    metric <- scr$metric
    cutoff <- scr$cutoff
    if (metric == "absZ") {
        m <- suppressWarnings(max(abs(as.numeric(z)), na.rm = TRUE))
        return(list(
            ok = is.finite(m) && m > cutoff,
            reason = sprintf(
                "no variant with |Z| above %g (max |Z| = %g)",
                cutoff,
                m
            )
        ))
    }
    ser <- susieR::susie_ser(z = z, n = n, coverage = NULL)
    if (metric == "pip") {
        eff <- if (cutoff < 0) 3 / nVar else cutoff
        return(list(
            ok = any(ser$pip > eff),
            reason = sprintf("no signals above PIP threshold %g", eff)
        ))
    }
    maxLbf <- suppressWarnings(max(as.numeric(ser$lbf_variable), na.rm = TRUE))
    if (metric == "logBf") {
        return(list(
            ok = is.finite(maxLbf) && maxLbf > cutoff,
            reason = sprintf(
                "no variant with logBF above %g (max logBF = %g)",
                cutoff,
                maxLbf
            )
        ))
    }
    # metric == "bf": compare in log space to avoid overflow of exp(maxLbf).
    list(
        ok = is.finite(maxLbf) && maxLbf > log(cutoff),
        reason = sprintf(
            "no variant with BF above %g (max BF = %g)",
            cutoff,
            exp(maxLbf)
        )
    )
}

# Per-entry signal screen. `screen` is a screen spec (see .asScreen). Skips
# (empties) the entry when the chosen metric shows no signal above its cutoff.
.applyEntryScreen <- function(df, n, screen) {
    scr <- .asScreen(screen)
    if (is.null(scr)) {
        return(list(df = df, skipped = FALSE))
    }
    res <- .entryScreenPass(df$Z, n = n, nVar = nrow(df), scr = scr)
    if (!res$ok) {
        return(list(
            df = slice(df, 0L),
            skipped = TRUE,
            reason = res$reason
        ))
    }
    list(df = df, skipped = FALSE)
}

# Back-compat thin wrapper for the original PIP-only screen (a bare numeric
# cutoff, 0 = off). Delegates to the generalized screen via .asScreen.
.applyPipScreen <- function(df, n, cutoff) {
    .applyEntryScreen(df, n = n, screen = cutoff)
}

# Prefix QC-track log lines with the entry label `lbl` (as `[lbl] ...`), or emit
# them bare when `lbl` is NA.
# @noRd
.qcEmit <- function(lbl, ...) {
    body <- str_c(...)
    if (is.na(lbl)) {
        inform(body)
    } else {
        msg <- glue("[{lbl}] {body}")
        inform(msg)
    }
}

# Internal: canonicalize the working per-variant `N`. The N source is resolved
# by a four-level priority: (1) per-variant N_CASE / N_CONTROL columns, (2)
# study -level nCase / nControl scalars, (3) a per-variant N column, (4) a
# study-level nSample scalar (total N). Levels 1-2 give the effective sample
# size (default) or the raw total (escape hatch). Returns list(df=, nSource=),
# where nSource is "effective" | "column" | "total" | "study-n" | NA_character_
# (no source). The entry label `lbl` prefixes the counts-win override log. A
# study-level scalar fills `df$N` even when the entry has no per-variant N
# column.
# --- .resolveEffectiveN helpers ---------------------------------------------

# Availability flags for the N-resolution decision tree.
.resolveNFlags <- function(df, opts) {
    hasScalar <- !is.null(opts$nCase) &&
        !is.null(opts$nControl) &&
        length(opts$nCase) == 1L &&
        length(opts$nControl) == 1L &&
        is.finite(opts$nCase) &&
        is.finite(opts$nControl) &&
        opts$nCase > 0 &&
        opts$nControl > 0
    hasNSample <- !is.null(opts$nSample) &&
        length(opts$nSample) == 1L &&
        is.finite(opts$nSample) &&
        opts$nSample > 0
    list(
        hasCols = all(is_in(c("N_CASE", "N_CONTROL"), colnames(df))),
        hasScalar = hasScalar,
        hasNSample = hasNSample,
        hasN = is_in("N", colnames(df)),
        nRow = nrow(df)
    )
}

# effectiveN off: prefer raw N -> raw total from counts -> study total N.
.resolveNRaw <- function(df, opts, f) {
    if (f$hasN) {
        return(list(df = df, nSource = "column"))
    }
    if (f$hasCols) {
        df$N <- as.numeric(df$N_CASE) + as.numeric(df$N_CONTROL)
        return(list(df = df, nSource = "total"))
    }
    if (f$hasScalar) {
        df$N <- rep(opts$nCase + opts$nControl, f$nRow)
        return(list(df = df, nSource = "total"))
    }
    if (f$hasNSample) {
        df$N <- rep(opts$nSample, f$nRow)
        return(list(df = df, nSource = "study-n"))
    }
    list(df = df, nSource = NA_character_)
}

# Default: per-variant c/c -> study c/c -> per-variant N -> study nSample.
.resolveNEffective <- function(df, opts, lbl, f) {
    if (f$hasCols) {
        if (f$hasN) {
            .qcEmit(
                lbl,
                "QC track: N overridden by effective N from per-variant ",
                "n_case/n_control."
            )
        }
        df$N <- effectiveN(df$N_CASE, df$N_CONTROL)
        return(list(df = df, nSource = "effective"))
    }
    if (f$hasScalar) {
        if (f$hasN) {
            .qcEmit(
                lbl,
                "QC track: N overridden by effective N from study ",
                "nCase/nControl."
            )
        }
        df$N <- rep(effectiveN(opts$nCase, opts$nControl), f$nRow)
        return(list(df = df, nSource = "effective"))
    }
    if (f$hasN) {
        return(list(df = df, nSource = "column"))
    }
    if (f$hasNSample) {
        df$N <- rep(opts$nSample, f$nRow)
        return(list(df = df, nSource = "study-n"))
    }
    list(df = df, nSource = NA_character_)
}

.resolveEffectiveN <- function(df, opts, lbl) {
    f <- .resolveNFlags(df, opts)
    if (!isTRUE(opts$effectiveN)) {
        return(.resolveNRaw(df, opts, f))
    }
    .resolveNEffective(df, opts, lbl, f)
}

# Internal: per-entry pipeline. Returns the cleaned GRanges and an audit list.
# --- .runEntrySummaryStatsQc: RAISS imputation step helpers ----------------

# Panel/dosage window indices for the entry, scoped per chromosome.
.qcRaissWindowIdx <- function(df, sketchSnpInfo, flank) {
    bounds <- tibble(
        chrom = str_remove(
            as.character(df$chrom),
            regex("^chr", ignore_case = TRUE)
        ),
        pos = as.integer(df$pos)
    ) |>
        group_by(.data$chrom) |>
        summarise(
            lo = min(.data$pos, na.rm = TRUE) - flank,
            hi = max(.data$pos, na.rm = TRUE) + flank,
            .groups = "drop"
        )
    # inner_join keeps only sketch SNPs whose chromosome appears in df (the
    # old is_in(skChrom, names(loByChr)) guard); filter keeps those inside the
    # [lo, hi] window. idx carries the original sketch row positions.
    tibble(
        chrom = str_remove(
            as.character(sketchSnpInfo$CHR),
            regex("^chr", ignore_case = TRUE)
        ),
        bp = as.integer(sketchSnpInfo$BP),
        idx = seq_along(sketchSnpInfo$CHR)
    ) |>
        inner_join(bounds, by = "chrom") |>
        filter(.data$bp >= .data$lo & .data$bp <= .data$hi) |>
        pull("idx")
}

# Assemble refPanel, knownZ table, and scaled dosage for the window.
.qcRaissBuildInputs <- function(df, ldSketch, windowIdx, sketchSnpInfo, opts) {
    refPanel <- .refVariantsFromSketch(ldSketch)[windowIdx, , drop = FALSE]
    refPanel$variant_id <- normalizeVariantId(refPanel$variant_id)
    refPanel <- refPanel[order(refPanel$pos), , drop = FALSE]
    knownVariantIds <- if (!is.null(df$SNP)) {
        as.character(df$SNP)
    } else {
        as.character(df$variant_id)
    }
    knownZ <- data.frame(
        chrom = as.character(df$chrom),
        pos = as.integer(df$pos),
        variant_id = knownVariantIds,
        A1 = as.character(df$A1),
        A2 = as.character(df$A2),
        z = as.numeric(df$Z),
        stringsAsFactors = FALSE
    )
    if (is_in("N", colnames(df))) {
        knownZ$n <- as.numeric(df$N)
    }
    if (is_in("BETA", colnames(df))) {
        knownZ$beta <- as.numeric(df$BETA)
    }
    if (is_in("SE", colnames(df))) {
        knownZ$se <- as.numeric(df$SE)
    }
    knownZ <- knownZ[order(knownZ$pos), , drop = FALSE]
    # meanImpute = FALSE so per-variant missingness is still visible; the
    # surviving columns are mean-imputed below, which is what meanImpute =
    # TRUE did.
    dosage <- .dosageMatrix(ldSketch, windowIdx, meanImpute = FALSE)
    colnames(dosage) <- normalizeVariantId(
        as.character(sketchSnpInfo$SNP[windowIdx])
    )
    dosage <- dosage[, refPanel$variant_id, drop = FALSE]
    keep <- .qcRaissTargetMask(refPanel, knownZ, dosage, opts)
    nDropped <- sum(!keep)
    refPanel <- refPanel[keep, , drop = FALSE]
    dosage <- dosage[, keep, drop = FALSE]
    scaledDosage <- scale(.qtlMeanImpute(dosage))
    scaledDosage[is.na(scaledDosage)] <- 0
    list(
        refPanel = refPanel,
        knownZ = knownZ,
        scaledDosage = scaledDosage,
        nDroppedTargets = nDropped
    )
}

# Which reference-panel variants to keep, given the MAF / MAC / missingness
# cutoffs in `imputeOpts`.
#
# The cutoffs govern what RAISS is willing to IMPUTE, not what it is willing to
# keep, so an observed variant survives whatever its frequency in the LD panel.
# That is not a convenience: `.raissSvdImpute` reaches `crossprod(V, zt)` with
# V from the panel's known columns and zt from `knownZscores$z`, so dropping an
# observed variant from the panel would leave those two out of step and pair
# z-scores with the wrong variants.
#
# Without this, every rare variant in the window of a large LD sketch becomes
# an imputation target, which is both slow and statistically pointless when the
# study is much smaller than the panel.
# @noRd
.qcRaissTargetMask <- function(refPanel, knownZ, dosage, opts) {
    mafCutoff <- opts$imputeOpts$mafCutoff %||% 0
    macCutoff <- opts$imputeOpts$macCutoff %||% 0
    imissCutoff <- opts$imputeOpts$imissCutoff %||% 1
    keep <- rep(TRUE, nrow(refPanel))
    if (mafCutoff <= 0 && macCutoff <= 0 && imissCutoff >= 1) {
        return(keep)
    }
    # Identified exactly as raissSingleMatrixFromX() does, so this mask cannot
    # disagree with the knowns/unknowns split it derives moments later.
    knownIds <- intersect(knownZ$variant_id, refPanel$variant_id)
    isTarget <- !is_in(refPanel$variant_id, knownIds)
    stats <- .qcRaissVariantStats(dosage)
    # A MAC cutoff is a MAF cutoff once expressed per panel sample; taking the
    # stricter of the two matches .qtlVariantFilters().
    nSamp <- nrow(dosage)
    effectiveMaf <- max(
        mafCutoff,
        if (nSamp > 0L) macCutoff / (2 * nSamp) else 0
    )
    fails <- (!is.na(stats$maf) & stats$maf < effectiveMaf) |
        stats$missRate > imissCutoff |
        is.na(stats$maf)
    keep & !(isTarget & fails)
}

# Per-variant MAF and missingness for the RAISS path. A thin alias over the
# shared `.panelVariantStats` (R/ld.R), so the analysis-time panel filter and
# the imputation-target filter cannot drift apart in how they measure a
# variant.
# @noRd
.qcRaissVariantStats <- function(dosage) {
    .panelVariantStats(dosage)
}

# Run RAISS with per-call option defaults.
.qcRaissRun <- function(inp, opts) {
    raiss(
        refPanel = inp$refPanel,
        knownZscores = inp$knownZ,
        genotypeMatrix = inp$scaledDosage,
        svdTol = opts$imputeOpts$svdTol %||% 1e-12,
        lamb = opts$imputeOpts$lamb %||% 0.01,
        r2Threshold = opts$imputeOpts$r2Threshold %||% 0.6,
        minimumLd = opts$imputeOpts$minimumLd %||% 5,
        verbose = FALSE
    )
}

# Convert the imputer output back into a QC tibble.
.qcRaissMerge <- function(imputed, knownZ, df) {
    if (is.null(imputed) || is.null(imputed$resultFilter)) {
        return(list(df = df, total = NA_integer_, imputed = 0L))
    }
    impDf <- imputed$resultFilter
    out <- tibble(
        chrom = impDf$chrom,
        pos = impDf$pos,
        SNP = impDf$variant_id,
        A1 = impDf$A1,
        A2 = impDf$A2,
        Z = impDf$z
    )
    if (is_in("n", colnames(impDf))) {
        out$N <- impDf$n
    }
    if (is_in("beta", colnames(impDf))) {
        out$BETA <- impDf$beta
    }
    if (is_in("se", colnames(impDf))) {
        out$SE <- impDf$se
    }
    if (is_in("N", colnames(out)) && any(is.na(out$N))) {
        out$N[is.na(out$N)] <- stats::median(out$N, na.rm = TRUE)
    }
    list(df = out, total = nrow(out), imputed = nrow(out) - nrow(knownZ))
}

# Emit the RAISS net-change QC track line.
.qcRaissReport <- function(imputeBefore, imputeAfter, lbl) {
    .qcEmit(
        lbl,
        "QC track: RAISS imputation ",
        imputeBefore,
        " -> ",
        imputeAfter,
        " variant(s) (net ",
        sprintf("%+d", imputeAfter - imputeBefore),
        ")."
    )
}

# Emit the harmonization QC track line (optional corrected/dropped detail).
.qcHarmonizeReport <- function(nOut, nHarmIn, counts, hasCounts, lbl) {
    detail <- if (hasCounts) {
        str_c(
            " (corrected: sign-flipped ",
            counts$harmCorrSign,
            ", strand-flipped ",
            counts$harmCorrStrand,
            "; dropped ",
            counts$harmDropped,
            ")"
        )
    } else {
        ""
    }
    .qcEmit(
        lbl,
        "QC track: harmonization kept ",
        nOut,
        " of ",
        nHarmIn,
        " variant(s)",
        detail,
        "."
    )
}

# Optional RAISS imputation against the ldSketch. Returns updated df + audit
# fields + before/after counts (imputation adds variants, so not monotonic).
.qcRaissImpute <- function(df, ldSketch, opts, lbl) {
    imputeBefore <- nrow(df)
    audit <- list()
    flank <- if (is.null(opts$imputeOpts$flank)) {
        0L
    } else {
        as.integer(opts$imputeOpts$flank)
    }
    sketchSnpInfo <- getSnpInfo(ldSketch)
    windowIdx <- .qcRaissWindowIdx(df, sketchSnpInfo, flank)
    if (length(windowIdx) == 0L) {
        .qcEmit(
            lbl,
            "QC track: RAISS imputation skipped ",
            "(no LD-panel variants in the ",
            "region window)."
        )
        audit$raissImputedVariants <- 0L
        return(list(
            df = df,
            audit = audit,
            imputeBefore = imputeBefore,
            imputeAfter = nrow(df)
        ))
    }
    inp <- .qcRaissBuildInputs(df, ldSketch, windowIdx, sketchSnpInfo, opts)
    if (inp$nDroppedTargets > 0L) {
        .qcEmit(
            lbl,
            "QC track: RAISS excluded ",
            as.character(inp$nDroppedTargets),
            " imputation target(s) below the MAF / missingness cutoffs."
        )
    }
    imputed <- .qcRaissRun(inp, opts)
    merged <- .qcRaissMerge(imputed, inp$knownZ, df)
    df <- merged$df
    if (!is.na(merged$total)) {
        audit$raissTotalVariants <- merged$total
    }
    audit$raissImputedVariants <- merged$imputed
    imputeAfter <- nrow(df)
    .qcRaissReport(imputeBefore, imputeAfter, lbl)
    list(
        df = df,
        audit = audit,
        imputeBefore = imputeBefore,
        imputeAfter = imputeAfter
    )
}

# --- .runEntrySummaryStatsQc: per-entry QC rollup helpers ------------------

# One "label N" removed-count segment, or NULL when nothing was removed.
.qcSeg <- function(val, label) {
    if (!is.null(val) && val > 0L) str_c(label, " ", val) else NULL
}

# Collect the per-step "removed" segments for the QC summary line.
.qcRemovedSegments <- function(entryAudit, qcCount, opts) {
    sc <- entryAudit$sanityChecks
    cf <- entryAudit$contentFilters
    segs <- c(
        .qcSeg(sc$nonstandardChrDropped, "nonstdChr"),
        .qcSeg(sc$missDataDropped, "missData"),
        .qcSeg(sc$pOutOfRangeDropped, "badP"),
        .qcSeg(sc$zeroEffectDropped, "zeroEffect"),
        .qcSeg(sc$nonpositiveSeDropped, "badSE"),
        .qcSeg(cf$mafDropped, "maf"),
        .qcSeg(cf$infoDropped, "info"),
        .qcSeg(cf$nDropped, "nCutoff"),
        .qcSeg(qcCount$harmDropped, "harmonization")
    )
    if (!identical(opts$zMismatchQc, "none")) {
        segs <- c(segs, str_c("mismatch ", qcCount$mismatchRemoved))
    }
    segs
}

# Emit the per-entry QC rollup: corrected (retained), removed, imputed.
.qcEmitRollup <- function(entryAudit, qcCount, opts, nStudyIn, nOut, lbl) {
    removedSegs <- .qcRemovedSegments(entryAudit, qcCount, opts)
    correctedSeg <- str_c(
        "sign-flip ",
        qcCount$harmCorrSign,
        ", strand-flip ",
        qcCount$harmCorrStrand
    )
    if (isTRUE(opts$alleleFlipKriging)) {
        correctedSeg <- str_c(
            correctedSeg,
            ", kriging-flip ",
            qcCount$krigingFlipped
        )
    }
    impSeg <- if (isTRUE(opts$impute) && !is.na(qcCount$imputeAfter)) {
        str_c(
            " | imputed ",
            sprintf("%+d", qcCount$imputeAfter - qcCount$imputeBefore)
        )
    } else {
        ""
    }
    .qcEmit(
        lbl,
        "QC summary: ",
        nStudyIn,
        " in -> ",
        nOut,
        " out",
        " | corrected: ",
        correctedSeg,
        if (length(removedSegs) > 0L) {
            str_c(" | removed: ", str_flatten(removedSegs, collapse = ", "))
        } else {
            ""
        },
        impSeg
    )
}

# --- .runEntrySummaryStatsQc: harmonization / kriging / mismatch steps ------

# Panel-vs-sumstats allele harmonization + counter bookkeeping.
.qcHarmonizeEntry <- function(df, ldSketch, opts, lbl) {
    nHarmIn <- nrow(df)
    df <- .matchAgainstSketch(
        df,
        ldSketch,
        matchMinProp = opts$matchMinProp,
        removeIndels = opts$removeIndels,
        removeStrandAmbiguous = opts$removeStrandAmbiguous,
        removeDups = TRUE
    )
    harmCounts <- attr(df, "qcCounts")
    attr(df, "qcCounts") <- NULL
    # Re-key SNP to the harmonized id: .matchAgainstSketch rewrites variant_id
    # to the panel orientation + sign-flips swapped variants but leaves SNP; a
    # stale SNP makes flipped variants miss the panel in later lookups.
    if (!is.null(df$variant_id)) {
        df$SNP <- df$variant_id
    }
    counts <- list(
        harmCorrSign = 0L,
        harmCorrStrand = 0L,
        harmDropped = nHarmIn - nrow(df)
    )
    if (!is.null(harmCounts)) {
        counts$harmCorrSign <- harmCounts$signFlip
        counts$harmCorrStrand <- harmCounts$strandFlip
    }
    .qcHarmonizeReport(nrow(df), nHarmIn, counts, !is.null(harmCounts), lbl)
    list(
        df = df,
        audit = list(matchedAgainstSketch = nrow(df)),
        counts = counts
    )
}

# Kriging allele-flip QC: sign-flip LD-inconsistent z-scores in place
# (susieR rule logLR > 2 & |z| > 2); variants are corrected and RETAINED.
.qcKrigingFlip <- function(df, ldSketch, opts, lbl) {
    if (!isTRUE(opts$alleleFlipKriging) || nrow(df) < 2L) {
        return(list(df = df, audit = list(), count = 0L))
    }
    nKrIn <- nrow(df)
    R <- .ldFromSketch(
        ldSketch,
        df$SNP,
        label = "summaryStatsQc: kriging prefilter"
    )
    nKrig <- if (!is.null(opts$nForPip) && is.finite(opts$nForPip)) {
        opts$nForPip
    } else {
        stats::median(as.numeric(df$N), na.rm = TRUE)
    }
    kr <- krigingOutlierQc(df$Z, R, n = nKrig, variantIds = df$SNP)
    nKr <- sum(kr$flip)
    if (nKr > 0L) {
        df$Z[kr$flip] <- -df$Z[kr$flip]
        if (is_in("BETA", colnames(df))) {
            df$BETA[kr$flip] <- -df$BETA[kr$flip]
        }
    }
    .qcEmit(
        lbl,
        "QC track: kriging sign-flipped ",
        nKr,
        " of ",
        nKrIn,
        " LD-inconsistent variant(s)."
    )
    list(
        df = df,
        count = nKr,
        audit = list(krigingFlipped = nKr, krigingDiagnostics = kr$diagnostics)
    )
}

# LD-mismatch QC (SLALOM / DENTIST). Retains full per-variant diagnostics.
.qcMismatchQc <- function(df, ldSketch, opts, lbl) {
    if (identical(opts$zMismatchQc, "none") || nrow(df) < 2L) {
        return(list(df = df, audit = list(), count = 0L))
    }
    nMmIn <- nrow(df)
    ldQc <- .applyLdMismatchQcToEntry(df, ldSketch, opts$zMismatchQc)
    df <- ldQc$df
    audit <- list(
        ldMismatchOutliersDropped = ldQc$outliers,
        ldMismatchMethod = opts$zMismatchQc
    )
    if (!is.null(ldQc$diagnostics)) {
        audit$ldMismatchDiagnostics <- ldQc$diagnostics
    }
    .qcEmit(
        lbl,
        "QC track: ",
        opts$zMismatchQc,
        " removed ",
        ldQc$outliers,
        " of ",
        nMmIn,
        " LD-mismatch outlier(s)."
    )
    list(df = df, audit = audit, count = ldQc$outliers)
}

# --- .runEntrySummaryStatsQc: setup + pre-harmonization step helpers --------

# Initialize per-entry QC state (df, audit, step counters, label).
.qcInitEntry <- function(gr, entryLabel) {
    df <- .entryGrangesToDf(gr)
    lbl <- if (!is.null(entryLabel) && isTRUE(str_length(entryLabel) > 0L)) {
        entryLabel
    } else {
        NA_character_
    }
    list(
        df = df,
        entryAudit = list(variantsIn = nrow(df)),
        qcCount = list(
            harmCorrSign = 0L,
            harmCorrStrand = 0L,
            harmDropped = 0L,
            krigingFlipped = 0L,
            mismatchRemoved = 0L,
            imputeBefore = NA_integer_,
            imputeAfter = NA_integer_
        ),
        lbl = lbl,
        nStudyIn = nrow(df)
    )
}

# Per-row sanity checks (bad P / zero effect / non-positive SE, small-P clamp,
# numeric coercion, CHR normalization, missing-data drop). Runs first.
.qcStepSanity <- function(df, opts, lbl) {
    nSanIn <- nrow(df)
    sanity <- .applySanityChecks(
        df,
        coerceNumeric = opts$coerceNumeric,
        normalizeChr = opts$normalizeChr,
        dropNonstandardChr = opts$dropNonstandardChr,
        dropMissData = opts$dropMissData,
        dropPOutOfRange = opts$dropPOutOfRange,
        clampSmallP = opts$clampSmallP,
        smallPFloor = opts$smallPFloor,
        dropZeroEffect = opts$dropZeroEffect,
        dropNonpositiveSe = opts$dropNonpositiveSe
    )
    df <- sanity$df
    audit <- list()
    if (length(sanity$audit) > 0L) {
        audit$sanityChecks <- sanity$audit
    }
    if (nSanIn > 0L && nrow(df) != nSanIn) {
        .qcEmit(
            lbl,
            "QC track: sanity checks kept ",
            nrow(df),
            " of ",
            nSanIn,
            " variant(s)."
        )
    }
    list(df = df, audit = audit)
}

# Canonicalize N to the effective sample size (case/control) before filtering,
# and keep nForPip consistent with the applied N. No-op for quantitative.
.qcStepEffectiveN <- function(df, opts, lbl) {
    nRes <- .resolveEffectiveN(df, opts, lbl)
    df <- nRes$df
    if (
        isTRUE(is_in(nRes$nSource, c("effective", "total"))) &&
            is_in("N", colnames(df))
    ) {
        opts$nForPip <- stats::median(as.numeric(df$N), na.rm = TRUE)
    }
    list(df = df, nSource = nRes$nSource, opts = opts)
}

# Variant-content filters (MAF / INFO / N).
.qcStepContentFilters <- function(df, opts, lbl) {
    nFiltIn <- nrow(df)
    cf <- .applyContentFilters(
        df,
        mafCutoff = opts$mafCutoff,
        infoCutoff = opts$infoCutoff,
        nCutoff = opts$nCutoff
    )
    df <- cf$df
    audit <- list()
    if (length(cf$audit) > 0L) {
        audit$contentFilters <- cf$audit
    }
    if (nFiltIn > 0L && nrow(df) != nFiltIn) {
        .qcEmit(
            lbl,
            "QC track: MAF/INFO/N filters kept ",
            nrow(df),
            " of ",
            nFiltIn,
            " variant(s)."
        )
    }
    list(df = df, audit = audit)
}

# Derive BETA/SE from signed Z, then P from Z (re-clamping tiny P).
.qcStepDerive <- function(df, opts, entryAudit) {
    derived <- .deriveBetaSeFromZ(df)
    df <- derived$df
    if (!is.null(derived$audit)) {
        entryAudit$betaSeFromZ <- derived$audit
    }
    if (is_in("Z", colnames(df)) && !is_in("P", colnames(df))) {
        df$P <- .zToPvalue(df$Z)
        entryAudit$pValueFromZ <- sum(!is.na(df$P))
        if (isTRUE(opts$clampSmallP) && nrow(df) > 0L) {
            smallMask <- !is.na(df$P) & df$P >= 0 & df$P < opts$smallPFloor
            nClamped <- sum(smallMask)
            if (nClamped > 0L) {
                df$P[smallMask] <- opts$smallPFloor
                prev <- entryAudit$sanityChecks$smallPClamped %||% 0L
                if (is.null(entryAudit$sanityChecks)) {
                    entryAudit$sanityChecks <- list()
                }
                entryAudit$sanityChecks$smallPClamped <- prev + nClamped
            }
        }
    }
    list(df = df, entryAudit = entryAudit)
}

# keepVariants subset + skipRegion drop.
.qcStepKeepSkip <- function(df, opts) {
    audit <- list()
    if (length(opts$keepVariants) > 0L) {
        before <- nrow(df)
        df <- filter(df, is_in(.data$SNP, opts$keepVariants))
        audit$keepVariantsDropped <- before - nrow(df)
    }
    if (!is.null(opts$skipRegion) && length(opts$skipRegion) > 0L) {
        before <- nrow(df)
        df <- .applySkipRegion(df, opts$skipRegion)
        audit$skipRegionDropped <- before - nrow(df)
    }
    list(df = df, audit = audit)
}

# Optional post-harmonization signal screen (PIP / |Z| / BF / logBF).
.qcStepScreen <- function(df, opts) {
    audit <- list()
    if (!is.null(opts$screen)) {
        scr <- .applyEntryScreen(df, n = opts$nForPip, screen = opts$screen)
        df <- scr$df
        audit$pipScreenSkipped <- isTRUE(scr$skipped)
        if (isTRUE(scr$skipped)) audit$pipScreenReason <- scr$reason
    }
    list(df = df, audit = audit)
}

# Optional RAISS imputation step wrapper (guarded; threads impute counters).
.qcRaissImputeStep <- function(df, ldSketch, opts, lbl, qcCount) {
    audit <- list()
    if (isTRUE(opts$impute) && nrow(df) >= 1L) {
        imp <- .qcRaissImpute(df, ldSketch, opts, lbl)
        df <- imp$df
        audit <- imp$audit
        qcCount$imputeBefore <- imp$imputeBefore
        qcCount$imputeAfter <- imp$imputeAfter
    }
    list(df = df, audit = audit, qcCount = qcCount)
}

# Early return when too few variants survive pre-harmonization QC.
.qcEarlyExit <- function(df, entryAudit, qcCount, opts, nIn, lbl) {
    entryAudit$earlyExit <-
        "fewer than two variants after pre-harmonization QC"
    # Still emit the rollup. An entry that QC empties is precisely the case a
    # user needs told about, and returning early used to make those drops
    # invisible in the log -- the audit recorded them, nothing said so.
    .qcEmitRollup(entryAudit, qcCount, opts, nIn, nrow(df), lbl)
    list(gr = .dfToEntryGranges(df), audit = entryAudit)
}

# Pre-harmonization phase: init + sanity + effective-N + content + derive +
# keep/skip. Returns the QC state carried into the harmonization phase.
.qcPreHarmonize <- function(gr, opts, entryLabel) {
    init <- .qcInitEntry(gr, entryLabel)
    df <- init$df
    entryAudit <- init$entryAudit
    lbl <- init$lbl
    san <- .qcStepSanity(df, opts, lbl)
    df <- san$df
    entryAudit <- modifyList(entryAudit, san$audit)
    eff <- .qcStepEffectiveN(df, opts, lbl)
    df <- eff$df
    entryAudit$nSource <- eff$nSource
    opts <- eff$opts
    cf <- .qcStepContentFilters(df, opts, lbl)
    df <- cf$df
    entryAudit <- modifyList(entryAudit, cf$audit)
    der <- .qcStepDerive(df, opts, entryAudit)
    df <- der$df
    entryAudit <- der$entryAudit
    ks <- .qcStepKeepSkip(df, opts)
    df <- ks$df
    entryAudit <- modifyList(entryAudit, ks$audit)
    list(
        df = df,
        entryAudit = entryAudit,
        opts = opts,
        qcCount = init$qcCount,
        lbl = lbl,
        nStudyIn = init$nStudyIn
    )
}

.runEntrySummaryStatsQc <- function(
    gr,
    ldSketch,
    refGenome,
    opts,
    entryLabel = NULL
) {
    pre <- .qcPreHarmonize(gr, opts, entryLabel)
    df <- pre$df
    entryAudit <- pre$entryAudit
    opts <- pre$opts
    qcCount <- pre$qcCount
    lbl <- pre$lbl
    if (nrow(df) < 2L) {
        return(.qcEarlyExit(df, entryAudit, qcCount, opts, length(gr), lbl))
    }
    harm <- .qcHarmonizeEntry(df, ldSketch, opts, lbl)
    df <- harm$df
    entryAudit <- modifyList(entryAudit, harm$audit)
    qcCount <- modifyList(qcCount, harm$counts)
    scr <- .qcStepScreen(df, opts)
    df <- scr$df
    entryAudit <- modifyList(entryAudit, scr$audit)
    kr <- .qcKrigingFlip(df, ldSketch, opts, lbl)
    df <- kr$df
    entryAudit <- modifyList(entryAudit, kr$audit)
    qcCount$krigingFlipped <- kr$count
    mm <- .qcMismatchQc(df, ldSketch, opts, lbl)
    df <- mm$df
    entryAudit <- modifyList(entryAudit, mm$audit)
    qcCount$mismatchRemoved <- mm$count
    imp <- .qcRaissImputeStep(df, ldSketch, opts, lbl, qcCount)
    df <- imp$df
    entryAudit <- modifyList(entryAudit, imp$audit)
    qcCount <- imp$qcCount
    .qcEmitRollup(entryAudit, qcCount, opts, pre$nStudyIn, nrow(df), lbl)
    entryAudit$variantsOut <- nrow(df)
    list(gr = .dfToEntryGranges(df), audit = entryAudit)
}

# Shrink an LD-sketch GenotypeHandle to the panel variants inside the summary
# statistics' per-chromosome position span. `entries` is a list/SimpleList of
# per-study (or per-tuple) GRanges. A genome-wide sketch otherwise carries a
# full-genome snpInfo; only variants inside [min,max] BP of each represented
# chromosome are reachable by harmonization or within-range imputation, so the
# rest is dropped at load time. NULL-safe; a no-op when the span already covers
# the panel. See [[.subsetGenotypeHandle]] for why this is read-safe.
# @noRd
.subsetSketchToRange <- function(ldSketch, entries) {
    if (is.null(ldSketch)) {
        return(NULL)
    }
    chrom <- unlist(map(entries, .entryChrom), use.names = FALSE)
    pos <- unlist(map(entries, .entryPos), use.names = FALSE)
    ok <- !is.na(chrom) & !is.na(pos)
    chrom <- chrom[ok]
    pos <- pos[ok]
    if (length(pos) == 0L) {
        return(ldSketch)
    }
    bounds <- tibble(chrom = chrom, pos = pos) |>
        group_by(chrom) |>
        summarise(lo = min(pos), hi = max(pos), .groups = "drop")
    si <- getSnpInfo(ldSketch)
    # left_join keeps every sketch SNP; one on a chromosome absent from the
    # entries gets lo/hi = NA -> inWindow FALSE (the old is_in guard). The
    # explicit is.na() guards force FALSE (never NA) for an absent chrom or an
    # NA bp, replacing the base keep[is.na(keep)] <- FALSE.
    keep <- tibble(
        chrom = canonChrom(as.character(si$CHR)),
        bp = as.integer(si$BP)
    ) |>
        left_join(bounds, by = "chrom") |>
        mutate(
            inWindow = !is.na(.data$lo) &
                !is.na(.data$bp) &
                .data$bp >= .data$lo &
                .data$bp <= .data$hi
        ) |>
        pull("inWindow")
    .subsetGenotypeHandle(ldSketch, keep)
}

# Shrink an LD-sketch GenotypeHandle to EXACTLY the variants present across the
# QC'd `entries` (imputation may have added variants; QC may have dropped some),
# matched by canonical variant id. Applied at the end of summaryStatsQc so the
# retained sketch mirrors the object's final variant set. NULL-safe.
# @noRd
.subsetSketchToIds <- function(ldSketch, entries) {
    if (is.null(ldSketch)) {
        return(NULL)
    }
    ids <- unlist(map(entries, .entrySnpIds), use.names = FALSE)
    if (length(ids) == 0L) {
        return(ldSketch)
    }
    si <- getSnpInfo(ldSketch)
    keep <- is_in(
        normalizeVariantId(as.character(si$SNP)),
        normalizeVariantId(unique(ids))
    )
    .subsetGenotypeHandle(ldSketch, keep)
}

# --- summaryStatsQc orchestration helpers ----------------------------------

# Validate input class + per-entry MAF/INFO column availability.
.ssqcCheckEntries <- function(sumstats, mafCutoff, infoCutoff) {
    if (
        !methods::is(sumstats, "QtlSumStats") &&
            !methods::is(sumstats, "GwasSumStats")
    ) {
        abort("summaryStatsQc requires a QtlSumStats or GwasSumStats input.")
    }
    for (i in seq_len(nrow(sumstats))) {
        cols <- colnames(S4Vectors::mcols(.collectionEntry(sumstats, i)))
        if (mafCutoff > 0 && !any(is_in(c("MAF", "FRQ"), cols))) {
            msg <- glue(
                "summaryStatsQc: mafCutoff > 0 requires every entry to ",
                "carry a MAF or FRQ column; entry {i} does not."
            )
            abort(msg)
        }
        if (infoCutoff > 0 && !is_in("INFO", cols)) {
            msg <- glue(
                "summaryStatsQc: infoCutoff > 0 requires every entry to ",
                "carry an INFO column; entry {i} does not."
            )
            abort(msg)
        }
    }
}

# Build the per-entry QC options list from the captured call parameters.
.ssqcBuildOpts <- function(p) {
    optNames <- c(
        "removeIndels",
        "removeStrandAmbiguous",
        "mafCutoff",
        "infoCutoff",
        "nCutoff",
        "skipRegion",
        "zMismatchQc",
        "alleleFlipKriging",
        "effectiveN",
        "impute",
        "imputeOpts",
        "matchMinProp",
        "coerceNumeric",
        "normalizeChr",
        "dropNonstandardChr",
        "dropMissData",
        "dropPOutOfRange",
        "clampSmallP",
        "smallPFloor",
        "dropZeroEffect",
        "dropNonpositiveSe"
    )
    opts <- p[optNames]
    opts$keepVariants <- as.character(p$keepVariants)
    opts$screen <- .resolveScreenMetric(
        p$pipCutoffToSkip,
        p$absZCutoffToSkip,
        p$bfCutoffToSkip,
        p$logBfCutoffToSkip
    )
    opts$nCase <- NULL
    opts$nControl <- NULL
    opts$nForPip <- NULL
    opts
}

# Per-entry sample-size options: median N for PIP, study case/control/total N.
.ssqcEntryOpts <- function(opts, sumstats, i) {
    mc <- S4Vectors::mcols(.collectionEntry(sumstats, i))
    opts$nForPip <- if (is_in("N", colnames(mc))) {
        stats::median(mc$N, na.rm = TRUE)
    } else {
        NULL
    }
    opts$nCase <- if (is_in("nCase", .tupleColumnNames(sumstats))) {
        as.numeric(sumstats$nCase)[[i]]
    } else {
        NULL
    }
    opts$nControl <- if (is_in("nControl", .tupleColumnNames(sumstats))) {
        as.numeric(sumstats$nControl)[[i]]
    } else {
        NULL
    }
    opts$nSample <- if (is_in("nSample", .tupleColumnNames(sumstats))) {
        as.numeric(sumstats$nSample)[[i]]
    } else {
        NULL
    }
    opts
}

# Per-entry log label: study/context/trait (QTL) or study (GWAS).
.ssqcEntryLabel <- function(sumstats, i, isQtl) {
    if (isQtl) {
        str_c(
            as.character(sumstats$study)[[i]],
            as.character(sumstats$context)[[i]],
            as.character(sumstats$trait)[[i]],
            sep = "/"
        )
    } else {
        as.character(sumstats$study)[[i]]
    }
}

# Run the per-entry QC pipeline across all entries.
.ssqcRunEntries <- function(sumstats, opts) {
    newEntries <- vector("list", nrow(sumstats))
    entryAudits <- vector("list", nrow(sumstats))
    isQtl <- methods::is(sumstats, "QtlSumStats")
    ldSketch <- getLdSketch(sumstats)
    refGenome <- getGenome(sumstats)
    for (i in seq_len(nrow(sumstats))) {
        opts <- .ssqcEntryOpts(opts, sumstats, i)
        result <- .runEntrySummaryStatsQc(
            gr = .collectionEntry(sumstats, i),
            ldSketch = ldSketch,
            refGenome = refGenome,
            opts = opts,
            entryLabel = .ssqcEntryLabel(sumstats, i, isQtl)
        )
        newEntries[[i]] <- result$gr
        entryAudits[[i]] <- result$audit
    }
    list(newEntries = newEntries, entryAudits = entryAudits)
}

# Assemble the qcInfo record (echoed options + per-entry audits).
.ssqcBuildQcInfo <- function(p, entryAudits) {
    optNames <- c(
        "removeIndels",
        "removeStrandAmbiguous",
        "mafCutoff",
        "infoCutoff",
        "nCutoff",
        "pipCutoffToSkip",
        "absZCutoffToSkip",
        "bfCutoffToSkip",
        "logBfCutoffToSkip",
        "zMismatchQc",
        "alleleFlipKriging",
        "effectiveN",
        "impute",
        "coerceNumeric",
        "normalizeChr",
        "dropNonstandardChr",
        "dropMissData",
        "dropPOutOfRange",
        "clampSmallP",
        "smallPFloor",
        "dropZeroEffect",
        "dropNonpositiveSe"
    )
    list(
        timestamp = NA_character_,
        options = p[optNames],
        entryAudit = entryAudits
    )
}

# Rebuild the SumStats object with QC'd entries, trimmed sketch, and qcInfo.
.ssqcRebuild <- function(sumstats, newEntries, newLdSketch, qcInfo) {
    has <- function(nm) is_in(nm, .tupleColumnNames(sumstats))
    nSample <- if (has("nSample")) as.numeric(sumstats$nSample) else NULL
    if (methods::is(sumstats, "GwasSumStats")) {
        GwasSumStats(
            study = as.character(sumstats$study),
            entry = newEntries,
            genome = getGenome(sumstats),
            ldSketch = newLdSketch,
            varY = as.numeric(sumstats$varY),
            nCase = if (has("nCase")) as.numeric(sumstats$nCase) else NULL,
            nControl = if (has("nControl")) {
                as.numeric(sumstats$nControl)
            } else {
                NULL
            },
            nSample = nSample,
            qcInfo = qcInfo,
            # Carry the existing block keys through: the entries are already
            # split, so without this the constructor re-derives them from
            # seqname and every block on one chromosome collapses to one key.
            blockId = if (has("blockId")) {
                as.character(sumstats$blockId)
            } else {
                NULL
            }
        )
    } else {
        QtlSumStats(
            study = as.character(sumstats$study),
            context = as.character(sumstats$context),
            trait = as.character(sumstats$trait),
            entry = newEntries,
            genome = getGenome(sumstats),
            ldSketch = newLdSketch,
            varY = as.numeric(sumstats$varY),
            nSample = nSample,
            qcInfo = qcInfo
        )
    }
}

#' Run QC on a SumStats Collection
#'
#' Applies a single QC pass to a \code{QtlSumStats} or \code{GwasSumStats}
#' collection: per-row sanity checks via \code{.applySanityChecks} (drop rows
#' with out-of-range / zero P, BETA == 0, SE <= 0, NA in vital columns; clamp
#' tiny P; normalize CHR; coerce signed columns to numeric), variant-content
#' filters (MAF / INFO / N) via \code{.applyContentFilters}, optional
#' \code{skipRegion} drop, optional PIP screen, panel-vs-sumstats allele
#' harmonization against the \code{ldSketch} via \code{harmonizeAlleles} (which
#' handles indels, strand-ambiguous variants, sign / strand flips, and duplicate
#' drops in a single sweep), optional SLALOM/DENTIST LD-mismatch QC, and
#' optional RAISS imputation. No Bioconductor genome / dbSNP packages required.
#'
#' The returned collection has its \code{qcInfo} slot populated with a per-entry
#' audit record (variant counts, drop counts at each step, which filters fired,
#' etc.). Fine-mapping and TWAS-weights pipelines reject SumStats inputs where
#' \code{length(getQcInfo(x)) == 0L}.
#'
#' Column-availability error contract: a non-zero \code{mafCutoff} requires
#' every entry to carry a \code{MAF} column; non-zero \code{infoCutoff} requires
#' \code{INFO}; non-zero \code{nCutoff} requires \code{N}. Missing column with a
#' non-zero cutoff is a hard error.
#'
#' @param sumstats A \code{QtlSumStats} or \code{GwasSumStats} collection.
#' @param removeIndels Logical (length 1). When \code{TRUE}, drop indels during
#'   panel harmonization. Default \code{FALSE}.
#' @param removeStrandAmbiguous Logical (length 1). When \code{TRUE}, drop A/T
#'   and C/G strand-ambiguous variants. Default \code{TRUE}.
#' @param mafCutoff Numeric (length 1). MAF threshold (variants with \code{MAF <
#'   mafCutoff} are dropped). Default 0. Requires \code{MAF} or \code{FRQ}
#'   column when non-zero.
#' @param infoCutoff Numeric (length 1). INFO score threshold. Default 0.
#'   Requires \code{INFO} column when non-zero.
#' @param nCutoff Numeric (length 1). Sample-size deviation threshold: drop
#'   variants whose \code{N} is more than \code{nCutoff}
#'   median-absolute-deviations from the median. Set to 0 to disable. Default 5.
#' @param keepVariants Optional character vector of variant IDs (SNP column) to
#'   retain prior to harmonization.
#' @param skipRegion Optional character vector of \code{"chr:start-end"}
#'   strings, or a \code{GRanges}, of regions to drop.
#' @param pipCutoffToSkip Numeric (length 1). When \code{!= 0}, run an
#'   LD-independent single-effect SER screen and skip the entry if no PIP
#'   exceeds the cutoff. \code{< 0} resolves to \code{3 / nVariants}. Default 0
#'   (no screen).
#' @param absZCutoffToSkip Numeric (length 1). Alternative signal screen: skip
#'   the entry when \code{max(abs(Z))} does not exceed the cutoff. No model fit.
#'   Default 0 (off).
#' @param bfCutoffToSkip Numeric (length 1). Alternative signal screen: skip the
#'   entry when the largest per-variant single-effect Bayes factor (from the
#'   same \code{susie_ser} fit as the PIP screen) does not exceed the cutoff.
#'   Compared in log space (\code{maxlogBF > log(cutoff)}). Must be \code{> 0}.
#'   Default 0 (off).
#' @param logBfCutoffToSkip Numeric (length 1). As \code{bfCutoffToSkip} but the
#'   cutoff is on the log Bayes factor scale (\code{maxlogBF > cutoff}). Default
#'   0 (off). Exactly one of \code{pipCutoffToSkip} / \code{absZCutoffToSkip} /
#'   \code{bfCutoffToSkip} / \code{logBfCutoffToSkip} may be non-zero (one
#'   screening metric at a time); enabling more than one is an error.
#' @param zMismatchQc One of \code{"none"} (default), \code{"slalom"},
#'   \code{"dentist"}.
#' @param alleleFlipKriging Logical (length 1). Opt-in kriging LD-consistency
#'   prefilter run before SLALOM/DENTIST. Default \code{FALSE}.
#' @param effectiveN Logical (length 1). When \code{TRUE} (default) and the
#'   input carries case/control counts --- per-variant \code{N_CASE} /
#'   \code{N_CONTROL} mcols, else the study-level \code{nCase} / \code{nControl}
#'   scalars --- the working per-variant \code{N} is set to the effective sample
#'   size \code{effectiveN(nCase, nControl)} BEFORE the N-cutoff filter, so the
#'   filter, kriging, and the downstream fit all use \code{N_eff}. When both
#'   counts and an \code{N} column are present the counts win: \code{N} is
#'   overridden and the override is logged. Inputs with no counts (quantitative
#'   traits) are unchanged. The escape hatch \code{effectiveN = FALSE} restores
#'   the raw \code{N} column (or, when there is no \code{N}, the raw total
#'   \code{nCase + nControl}) with no override. \code{qcInfo$options$effectiveN}
#'   records the setting and each entry's \code{nSource} is one of
#'   \code{"effective"}, \code{"column"}, \code{"total"}, or \code{NA}.
#' @param impute Logical (length 1). Run RAISS imputation against the
#'   \code{ldSketch}. Default \code{FALSE}. (Note: RAISS against the sketch is
#'   not yet fully wired for the new path; the option is accepted but currently
#'   emits a warning and is skipped.)
#' @param imputeOpts Named list of RAISS parameters. RAISS imputation scopes its
#'   reference panel to the analysis-region window (so a per-chromosome /
#'   genome-wide \code{ldSketch} does not materialize its full dosage);
#'   \code{flank} (default 0) widens that window by the given number of base
#'   pairs on each side to retain LD context for edge variants.
#'
#'   \code{mafCutoff} (default 0), \code{macCutoff} (default 0) and
#'   \code{imissCutoff} (default 1) bound which panel variants RAISS will
#'   attempt to impute. Without them every rare variant in the window of a
#'   large LD sketch becomes an imputation target, which is slow and of little
#'   value when the study is far smaller than the panel. The stricter of
#'   \code{mafCutoff} and \code{macCutoff / (2 * nSamples)} applies, matching
#'   \code{\link{QtlDataset}}.
#'
#'   These bound what is \strong{imputed}, not what is kept: a variant present
#'   in the sumstats survives whatever its frequency in the panel, both because
#'   it is observed data and because RAISS derives its LD basis from those same
#'   panel rows.
#' @param matchMinProp Minimum proportion of LD panel variants that must be
#'   matched by the sumstats; default 0.
#' @param coerceNumeric Logical. Coerce signed columns
#'   (Z/BETA/SE/OR/LOG_ODDS/SIGNED_SUMSTAT/P/MAF/FRQ/INFO/N) to numeric. Default
#'   \code{TRUE}.
#' @param normalizeChr Logical. Strip \code{"chr"} prefix, uppercase the
#'   chromosome label, and map 23->X, 24->Y, M->MT. Default \code{TRUE}.
#' @param dropNonstandardChr Logical. Drop variants whose CHR (after
#'   normalization) is outside 1..22, X, Y, MT. Default \code{TRUE}.
#' @param dropMissData Logical. Drop rows with NA in any vital column (chrom,
#'   pos, A1, A2, and at least one of Z / BETA). Default \code{TRUE}.
#' @param dropPOutOfRange Logical. Drop rows where \code{P < 0} or \code{P > 1}.
#'   Default \code{TRUE}.
#' @param clampSmallP Logical. Floor non-negative P values below
#'   \code{smallPFloor} to \code{smallPFloor} so \code{-log10(P)} stays finite.
#'   Applied to both input and Z-derived P values. Default \code{TRUE}.
#' @param smallPFloor Numeric (length 1). Floor for \code{clampSmallP}. Default
#'   \code{5e-324} (R's smallest positive double).
#' @param dropZeroEffect Logical. Drop rows where any effect column is exactly 0
#'   (\code{BETA}, \code{LOG_ODDS}, \code{SIGNED_SUMSTAT}) or \code{OR} is
#'   exactly 1. Default \code{TRUE}.
#' @param dropNonpositiveSe Logical. Drop rows where \code{SE <= 0}. Default
#'   \code{TRUE}.
#' @return A new \code{QtlSumStats} / \code{GwasSumStats} with cleaned entries
#'   and \code{qcInfo} populated.
#' @examples
#' data(gwasSumStatsS4Example)
#' summaryStatsQc(gwasSumStatsS4Example)
#' @export
summaryStatsQc <- function(
    sumstats,
    removeIndels = FALSE,
    removeStrandAmbiguous = TRUE,
    mafCutoff = 0,
    infoCutoff = 0,
    nCutoff = 5,
    keepVariants = NULL,
    skipRegion = NULL,
    pipCutoffToSkip = 0,
    absZCutoffToSkip = 0,
    bfCutoffToSkip = 0,
    logBfCutoffToSkip = 0,
    zMismatchQc = c("none", "slalom", "dentist"),
    alleleFlipKriging = FALSE,
    effectiveN = TRUE,
    impute = FALSE,
    imputeOpts = list(
        rcond = 0.01,
        r2Threshold = 0.6,
        minimumLd = 5,
        lamb = 0.01
    ),
    matchMinProp = 0,
    coerceNumeric = TRUE,
    normalizeChr = TRUE,
    dropNonstandardChr = TRUE,
    dropMissData = TRUE,
    dropPOutOfRange = TRUE,
    clampSmallP = TRUE,
    smallPFloor = 5e-324,
    dropZeroEffect = TRUE,
    dropNonpositiveSe = TRUE
) {
    zMismatchQc <- arg_match(zMismatchQc)
    .ssqcCheckEntries(sumstats, mafCutoff, infoCutoff)
    p <- as.list(environment())
    opts <- .ssqcBuildOpts(p)
    res <- .ssqcRunEntries(sumstats, opts)
    qcInfo <- .ssqcBuildQcInfo(p, res$entryAudits)
    newLdSketch <- .subsetSketchToIds(getLdSketch(sumstats), res$newEntries)
    .ssqcRebuild(sumstats, res$newEntries, newLdSketch, qcInfo)
}

# ---- slidingWindowLoop callbacks (distance mode) ------------------------
# `ctx` bundles the caller's segmentation state; see .segByDistRun.

# @noRd
.segDistMinBlock <- function(blockSize, ctx) {
    blockSize >= ctx$minBlockSize / 2 && (blockSize - ctx$minDim) >= 0
}

# @noRd
.segDistInitEnd <- function(startIdx, blockEnd, ctx) {
    min(.nthQuaterIdx(startIdx, 4, ctx$quaterIdx) + 1, blockEnd)
}

# Distance mode: fill is always q1 to q3 (inner 50% by distance); first/last
# corrections are handled by fix_block_fills in the loop.
# @noRd
.segDistFill <- function(
    startIdx,
    endIdx,
    notStartInterval,
    notLastInterval,
    ctx
) {
    list(
        start = .nthQuaterIdx(startIdx, 1, ctx$quaterIdx),
        end = .nthQuaterIdx(startIdx, 3, ctx$quaterIdx)
    )
}

# @noRd
.segDistStep <- function(startIdx, blockEnd, ctx) {
    .segByDistStep(startIdx, blockEnd, ctx$quaterIdx)
}

# If the last interval is small, go back one step.
# @noRd
.segDistAdjustLast <- function(startIdx, oldStartIdx, endIdx, blockEnd, ctx) {
    q1Old <- .nthQuaterIdx(oldStartIdx, 1, ctx$quaterIdx)
    small <- as.numeric(ctx$pos[min(endIdx - 1, ctx$n)]) -
        as.numeric(ctx$pos[q1Old]) <
        ctx$cutoff
    if (small) q1Old else startIdx
}

# ---- slidingWindowLoop callbacks (count mode) ---------------------------
# `ctx` bundles the caller's segmentation state; see segmentByCount.

# @noRd
.segCountMinBlock <- function(blockSize, ctx) {
    blockSize >= ctx$half
}

# @noRd
.segCountInitEnd <- function(startIdx, blockEnd, ctx) {
    if (blockEnd - ctx$half > startIdx + ctx$cutoff) {
        startIdx + ctx$cutoff
    } else {
        blockEnd
    }
}

# Count mode: fill based on index arithmetic (inner 50%).
# @noRd
.segCountFill <- function(
    startIdx,
    endIdx,
    notStartInterval,
    notLastInterval,
    ctx
) {
    list(
        start = if (notStartInterval) startIdx + ctx$quarter else startIdx,
        end = if (notLastInterval) endIdx - ctx$quarter else endIdx
    )
}

# @noRd
.segCountStep <- function(startIdx, blockEnd, ctx) {
    nextStart <- startIdx + ctx$half
    endIdx <- if (blockEnd - ctx$half > nextStart + ctx$cutoff) {
        nextStart + ctx$cutoff
    } else {
        blockEnd
    }
    list(startIdx = nextStart, endIdx = endIdx)
}

# ---- other map/apply helpers (lambda-free callbacks) --------------------

# One credible set's scalar summary row (top variant / PIP / z / p). The
# between-CS correlation columns are appended once by .csAppendCorrelationCols.
# @noRd
.extractCsInfoRow <- function(i, csNames, trimmed, variantNames, topLociTable) {
    csName <- csNames[i]
    indices <- trimmed$sets$cs[[csName]]
    csVariants <- variantNames[indices]
    csData <- filter(topLociTable, is_in(.data$variant_id, csVariants))
    topRow <- which.max(csData$pip)
    topVariant <- csData$variant_id[topRow]
    topZ <- csData$z[topRow]
    tibble(
        cs_name = csName,
        variants_per_cs = length(csVariants),
        top_variant = topVariant,
        top_variant_index = which(variantNames == topVariant),
        top_pip = csData$pip[topRow],
        top_z = topZ,
        p_value = .zToPvalue(topZ)
    )
}

# TRUE when credible set row `i` is a tagged (redundant) set.
# @noRd
.autoDecisionTagged <- function(i, df, highCorrCols) {
    if (df$top_cs[i]) {
        return(FALSE)
    }
    if (df$p_value[i] > 1e-4) {
        return(TRUE)
    }
    if (length(highCorrCols) == 0) {
        return(FALSE)
    }
    rowVals <- df |>
        filter(row_number() == i) |>
        select(all_of(highCorrCols))
    any(rowVals == 1)
}

# One block's variant-name data.frame from its LD matrix columns.
# @noRd
.ldVariantsDf <- function(ld) {
    # colnames(ld) is NULL for an empty / no-dimnames block. Make that explicit
    # as an empty character vector so the returned tibble ALWAYS carries a
    # `variants` column: .ldMergeVariants then reads $variants as character(0)
    # (length 0 -> the block is skipped) instead of hitting an absent column
    # (which a data.frame returned silently as NULL and a tibble warns on).
    variants <- colnames(ld)
    if (is.null(variants)) {
        variants <- character(0)
    }
    tibble(variants = variants)
}

# The is-NA mask of column `col` in `df` (for the vital-column drop reduce).
# @noRd
.scColIsNa <- function(col, df) {
    is.na(df[[col]])
}

# One GRanges entry's canonical chromosome vector.
# @noRd
.entryChrom <- function(gr) {
    canonChrom(as.character(GenomicRanges::seqnames(gr)))
}

# One GRanges entry's integer start positions.
# @noRd
.entryPos <- function(gr) {
    as.integer(GenomicRanges::start(gr))
}

# One GRanges entry's SNP ids (empty for NULL / no SNP mcol).
# @noRd
.entrySnpIds <- function(gr) {
    if (is.null(gr)) {
        return(character(0))
    }
    snp <- S4Vectors::mcols(gr)$SNP
    if (is.null(snp)) character(0) else as.character(snp)
}
