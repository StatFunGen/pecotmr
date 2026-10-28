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

#' @include MethodParam.R
NULL

#' @importFrom GenomicRanges seqnames
#' @importFrom S4Vectors mcols
NULL

# Allele harmonization (harmonizeAlleles) now lives in R/variantId.R, alongside
# the other variant-matching primitives (parseVariantId, matchVariants). It is
# package-internal and called from here (.matchAgainstLdRef), ctwasPipeline,
# and the pipeline join sites.

# Standardize a variant set (GRanges or data.frame) to a chrom/pos/alt/ref
# frame.
# @noRd
.variantsToDf <- function(x) {
    df <- if (is(x, "GRanges")) {
        as_tibble(mutate(
            as.data.frame(mcols(x)),
            chrom = as.character(seqnames(x)),
            pos = start(x)
        ))
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
#' @importFrom checkmate assertFlag
mergeVariantInfo <- function(variants1, variants2, all = TRUE) {
    assertFlag(all)
    df1 <- .variantsToDf(variants1)
    df2 <- .variantsToDf(variants2)

    key1 <- .canonicalAlleleKey(df1)
    key2 <- .canonicalAlleleKey(df2)

    # Detect flips: where df2's alt matches df1's ref at the same key
    matchIdx <- match(key2, key1)
    hasMatch <- !is.na(matchIdx)

    mi <- matchIdx[hasMatch]
    flip <- replace(
        rep(FALSE, nrow(df2)),
        hasMatch,
        df2$alt[hasMatch] == df1$ref[mi] & df2$ref[hasMatch] == df1$alt[mi]
    )

    # Apply flips to df2. Both replacements read the ORIGINAL `df2`, so the
    # pair swaps rather than each taking the other's already-swapped value.
    flipRows <- which(hasMatch)[flip[hasMatch]]
    flipped <- mutate(
        df2,
        alt = replace(.data$alt, flipRows, df2$ref[flipRows]),
        ref = replace(.data$ref, flipRows, df2$alt[flipRows])
    )

    if (all) {
        distinct(bind_rows(df1, flipped))
    } else {
        flipped
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
    .ldRequireOneSource(R, X)
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

# Require exactly one of R (LD matrix) or X (genotype matrix).
# @noRd
.ldRequireOneSource <- function(R, X) {
    if (is.null(R) && is.null(X)) {
        abort("Either R (LD matrix) or X (genotype matrix) must be provided.")
    }
    if (!is.null(R) && !is.null(X)) {
        abort("Provide either R or X, not both.")
    }
    invisible(NULL)
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
    withPos <- if (is_in("pos", lc)) {
        sumStat
    } else {
        `colnames<-`(
            sumStat,
            replace(colnames(sumStat), which(is_in(lc, "position")), "pos")
        )
    }
    named <- if (is_in("z", lc)) {
        withPos
    } else {
        `colnames<-`(
            withPos,
            replace(colnames(withPos), which(is_in(lc, "zscore")), "z")
        )
    }
    arrange(named, .data$pos)
}

# Run DENTIST on a single window, unpacking the shared tuning parameters.
# `ld` is .dentistLdInput()'s list(R) or list(F); a factor goes in as X, which
# dentistSingleWindow() reads without forming the LD.
.dentistCallSingle <- function(zScore, ld, nSample, dentistParam) {
    # dentistSingleWindow is what DentistParam() is built for -- ldMismatchQc
    # splices the same bundle into it -- so the tunables travel as one value
    # rather than being named again at every hop down from dentist().
    exec(
        dentistSingleWindow,
        zScore,
        R = ld$R,
        X = ld$F,
        nSample = nSample,
        !!!dentistParam
    )
}

# Segment into windows, run DENTIST per window, and merge the results.
# DENTIST on window `k`.
# @noRd
.dentistWindowAt <- function(
    k,
    sumStat,
    ld,
    nSample,
    windowDividedRes,
    dentistParam
) {
    # windowEndIdx is 1-based exclusive; convert to an inclusive range.
    idxRange <- windowDividedRes$windowStartIdx[
        k
    ]:(windowDividedRes$windowEndIdx[k] - 1L)
    .dentistCallSingle(
        sumStat$z[idxRange],
        .dentistLdSlice(ld, idxRange),
        nSample,
        dentistParam
    )
}

# One window's LD: the matrix's diagonal block, or the factor's columns.
# @noRd
.dentistLdSlice <- function(ld, idx) {
    if (is.null(ld$F)) {
        return(list(R = ld$R[idx, idx, drop = FALSE]))
    }
    list(F = ld$F[, idx, drop = FALSE])
}

# The LD a DENTIST run reads, with the sample size it needs. Genotypes under
# the sample estimator give their factor (.ldFactorFromDosage()), which
# DENTIST's factor kernel reads without forming the p x p LD; an LD matrix,
# or genotypes under another estimator, give the matrix (resolveLdInput()).
# Without `nSample`, genotypes count their own rows, as resolveLdInput()
# does.
# @noRd
.dentistLdInput <- function(R, X, nSample, ldMethod) {
    .ldRequireOneSource(R, X)
    if (is.null(R) && identical(ldMethod, "sample")) {
        X <- as.matrix(X)
        return(list(F = .ldFactorFromDosage(X), nSample = nSample %||% nrow(X)))
    }
    ld <- resolveLdInput(
        R = R,
        X = X,
        nSample = nSample,
        needNSample = TRUE,
        ldMethod = ldMethod
    )
    list(R = ld$R, nSample = ld$nSample)
}

.dentistWindows <- function(
    sumStat,
    ld,
    nSample,
    windowMode,
    windowSize,
    minDim,
    dentistParam
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
    dentistResultByWindow <- map(
        seq_len(nrow(windowDividedRes)),
        .dentistWindowAt,
        sumStat = sumStat,
        ld = ld,
        nSample = nSample,
        windowDividedRes = windowDividedRes,
        dentistParam = dentistParam
    )
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
#' @param X Genotype matrix (samples x SNPs). With \code{ldMethod = "sample"}
#'   DENTIST works on the standardized genotypes and never forms the LD
#'   matrix; any other \code{ldMethod} computes LD via \code{computeLd(X)}.
#'   \code{nSample} defaults to \code{nrow(X)}.
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
#' @param dentistParam The DENTIST algorithm settings, built with
#'   \code{\link{DentistParam}}: \code{pValueThreshold}, \code{propSVD},
#'   \code{gcControl}, \code{nIter}, \code{gPvalueThreshold},
#'   \code{duprThreshold}, \code{BPPARAM}, \code{correctChenEtAlBug}
#'   and \code{seed}. This is the same record
#'   \code{\link{ldMismatchQc}} takes as
#'   \code{method = DentistParam(...)}, so a setting means the same thing
#'   and carries the same default whichever entry point reaches it.
#' @param minDim In distance mode: minimum number of SNPs per block (default
#'   2000). In count mode: the number of variants per window (i.e., the window
#'   size).
#' @param ldMethod Character string specifying the LD computation method when
#'   \code{X} is provided: \code{"sample"} (default), which uses the
#'   genotypes directly, or \code{"population"} or \code{"gcta"}, which go
#'   through \code{computeLd}. Ignored when \code{R} is provided directly.
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
    minDim = 2000,
    ldMethod = "sample",
    dentistParam = DentistParam()
) {
    # The tunables arrive as the record rather than as nine formals of their
    # own: DentistParam() is already what ldMismatchQc(method =) takes and
    # what every helper below is spliced with, and duplicating its defaults
    # here left two places to change one number.
    .assertMethodParam(dentistParam, "DentistParam")
    # Genotypes become one factor of the LD here, and each window reads its
    # columns; the region-wide LD is never formed.
    ld <- .dentistLdInput(R, X, nSample, ldMethod)
    nSample <- ld$nSample
    sumStat <- .dentistResolveColumns(sumStat)
    windowMode <- arg_match(windowMode)
    if (nrow(sumStat) < minDim) {
        return(.dentistCallSingle(sumStat$z, ld, nSample, dentistParam))
    }
    .dentistWindows(
        sumStat,
        ld,
        nSample,
        windowMode,
        windowSize,
        minDim,
        dentistParam
    )
}

# --- dentistSingleWindow helpers -------------------------------------------

# Warn on small windows; validate the LD's shape against zScore: a square
# matrix over the variants, or a factor with one column per variant.
.dentistValidateInput <- function(zScore, ld) {
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
    if (!is.null(ld$F)) {
        if (ncol(ld$F) != length(zScore)) {
            abort(glue(
                "The genotypes must have one column per z-score ",
                "({ncol(ld$F)} columns, {length(zScore)} z-scores)."
            ))
        }
        return(invisible(NULL))
    }
    ldMat <- ld$R
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
.dentistDedup <- function(zScore, ld, duprThreshold) {
    if (duprThreshold >= 1.0) {
        return(list(zScore = zScore, ld = ld, dedupRes = NULL))
    }
    rThreshold <- round(sqrt(duprThreshold) * 1000) / 1000
    dedupRes <- if (is.null(ld$F)) {
        .findDuplicateVariants(zScore, ld$R, rThreshold)
    } else {
        .findDuplicateVariantsFactor(zScore, ld$F, rThreshold)
    }
    numDup <- sum(dedupRes$dupBearer != -1)
    if (numDup > 0) {
        nZ <- length(zScore)
        msg <- glue(
            "{numDup} duplicated variants out of a total of {nZ} ",
            "were found at r threshold of {rThreshold}"
        )
        inform(msg)
    }
    list(
        zScore = dedupRes$filteredZ,
        ld = if (is.null(ld$F)) {
            list(R = dedupRes$filteredLD)
        } else {
            list(F = dedupRes$filteredFactor)
        },
        dedupRes = dedupRes
    )
}

# Run the C++ iterative imputation and snake_case the output. The C++ returns
# any rsq values it capped at 1.0 in `rsqExceed`; summarize them into a single
# warning here (no warning handler / shared env needed).
.dentistRunImpute <- function(ld, nSample, zScore, p) {
    verboseIter <- getOption("pecotmr.dentist.verbose", FALSE)
    # The factor kernel imputes from the LD's factor without forming it.
    kernel <- if (is.null(ld$F)) {
        dentistIterativeImpute
    } else {
        dentistIterativeImputeFactor
    }
    raw <- kernel(
        # cpp11 requires exact integer types for int parameters
        ld$R %||% ld$F,
        as.integer(nSample),
        zScore,
        p$pValueThreshold,
        p$propSVD,
        p$gcControl,
        as.integer(p$nIter),
        p$gPvalueThreshold,
        p$threads,
        p$correctChenEtAlBug,
        verboseIter,
        if (is.null(p$seed)) NULL else as.integer(p$seed)
    )
    rsqExceed <- raw$rsqExceed
    res <- list_modify(raw, rsqExceed = zap())
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
#' @param X Genotype matrix (samples x SNPs). With \code{ldMethod = "sample"}
#'   DENTIST works on the standardized genotypes and never forms the LD
#'   matrix; any other \code{ldMethod} computes LD via \code{computeLd(X)}.
#'   \code{nSample} defaults to \code{nrow(X)}.
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
#' @param BPPARAM A \code{BiocParallelParam} whose workers number the
#'   threads DENTIST's OpenMP loops use when they run on this machine
#'   (\code{SerialParam}, \code{MulticoreParam}); any other backend gives one
#'   thread. Default \code{SerialParam()}.
#' @param correctChenEtAlBug Correct the original DENTIST operator! bug. Default
#'   is TRUE.
#' @param ldMethod Character string specifying the LD computation method when
#'   \code{X} is provided: \code{"sample"} (default), which uses the
#'   genotypes directly, or \code{"population"} or \code{"gcta"}, which go
#'   through \code{computeLd}. Ignored when \code{R} is provided directly.
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
    BPPARAM = SerialParam(),
    correctChenEtAlBug = TRUE,
    ldMethod = "sample",
    seed = NULL
) {
    .assertBpparam(BPPARAM)
    ld <- .dentistLdInput(R, X, nSample, ldMethod)
    nSample <- ld$nSample
    .dentistValidateInput(zScore, ld)
    p <- list(
        pValueThreshold = pValueThreshold,
        propSVD = propSVD,
        gcControl = gcControl,
        nIter = nIter,
        gPvalueThreshold = gPvalueThreshold,
        threads = .bpLocalThreads(BPPARAM, "dentistSingleWindow"),
        correctChenEtAlBug = correctChenEtAlBug,
        seed = seed
    )
    orgZscore <- zScore
    dedup <- .dentistDedup(zScore, ld, duprThreshold)
    imputed <- .dentistRunImpute(dedup$ld, nSample, dedup$zScore, p)
    res <- if (duprThreshold < 1.0) {
        addDupsBackDentist(orgZscore, imputed, dedup$dedupRes)
    } else {
        imputed
    }
    .dentistOutlierStat(res, pValueThreshold)
}

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
.dentistBuildAssignIdx <- function(dupBearer) {
    # Non-duplicates take the next free slot, which is just how many
    # non-duplicates have been seen so far; duplicates point at their bearer.
    isNew <- dupBearer == -1
    if_else(isNew, as.numeric(cumsum(isNew)), as.numeric(dupBearer))
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
    # Every row is independent of every other, so the whole table is built in
    # one shot. Duplicates sign-flip imputed_z and recompute z_diff from their
    # own z-score, so z_diff^2 matches the binary stat (DENTIST.h l706).
    isDup <- dupBearer != -1
    originalZ <- zScore[seq_len(nrowsDup)]
    rsqRow <- rsq[assignIdx]
    imputedRow <- imputedZ[assignIdx] * if_else(isDup, sign, 1)
    denom <- sqrt(pmax(1 - rsqRow, 1e-8))
    tibble(
        original_z = originalZ,
        imputed_z = imputedRow,
        iter_to_correct = iterToCorrect[assignIdx],
        rsq = rsqRow,
        z_diff = if_else(
            isDup,
            (originalZ - imputedRow) / denom,
            zDiff[assignIdx]
        ),
        is_duplicate = isDup
    )
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
addDupsBackDentist <- function(zScore, dentistOutput, findDupOutput) {
    dupBearer <- findDupOutput$dupBearer
    nrowsDup <- length(dupBearer)
    .dentistValidateDups(zScore, dentistOutput, dupBearer, nrowsDup)
    assignIdx <- .dentistBuildAssignIdx(dupBearer)
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
    allGaps <- c(1L, which(diffs > gapThreshold) + 1L, n + 1L)

    if (verbose && length(allGaps) - 2 > 0) {
        nGaps <- length(allGaps) - 2
        inform(glue("No. of gaps found: {nGaps}"))
        for (i in 2:(length(allGaps) - 1)) {
            gapNo <- i - 1
            startPos <- pos[allGaps[i] - 1]
            endPos <- pos[allGaps[i]]
            inform(glue("  Gap {gapNo}: {startPos} - {endPos}", .trim = FALSE))
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

# One field across a block's window records, typed as `c()` would have left
# it (empty when the block produced no window).
# @noRd
.swlField <- function(windows, field) {
    if (length(windows) == 0L) {
        return(integer(0))
    }
    list_c(map(windows, field))
}

# The loop itself stays: where the next window starts depends on the window
# just emitted, and the block ends when the stepper says so. What is gone is
# the four parallel `c()` accumulators -- each window is one record now, and
# the vectors are read off at the end.
.swlBlockWindows <- function(blockStart, blockEnd, fns, ctx) {
    startIdx <- blockStart
    endIdx <- fns$initEndFn(startIdx, blockEnd, ctx)
    oldStartIdx <- startIdx
    notStartInterval <- FALSE
    windows <- list()
    repeat {
        if (length(windows) >= 400) {
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
        windows[[length(windows) + 1L]] <- list(
            start = win$startIdx,
            end = min(endIdx, blockEnd),
            fillStart = win$fills$start,
            fillEnd = win$fills$end
        )
        if (!win$notLastInterval) {
            break
        }
        oldStartIdx <- win$startIdx
        stepped <- fns$stepFn(win$startIdx, blockEnd, ctx)
        startIdx <- stepped$startIdx
        endIdx <- stepped$endIdx
        notStartInterval <- TRUE
    }
    list(
        starts = .swlField(windows, "start"),
        ends = .swlField(windows, "end"),
        fillStarts = .swlField(windows, "fillStart"),
        fillEnds = .swlField(windows, "fillEnd")
    )
}

# Is block `k` long enough to window?
# @noRd
.swlBlockQualifies <- function(k, allGaps, minBlockFn, ctx) {
    minBlockFn(allGaps[k + 1] - allGaps[k], ctx)
}

# Block `k`'s windows, with the block's outer fill bounds snapped to the
# block's own edges. `replace()` returns a copy, so nothing is mutated.
# @noRd
.swlBlockAt <- function(k, allGaps, fns, ctx) {
    w <- .swlBlockWindows(allGaps[k], allGaps[k + 1], fns, ctx)
    # First window's fill starts at the window start; last window's fill ends
    # at the window end.
    list(
        starts = w$starts,
        ends = w$ends,
        fillStarts = replace(w$fillStarts, 1L, w$starts[1L]),
        fillEnds = replace(
            w$fillEnds,
            length(w$fillEnds),
            w$ends[length(w$ends)]
        )
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
    fns <- list(
        minBlockFn = minBlockFn,
        initEndFn = initEndFn,
        fillFn = fillFn,
        stepFn = stepFn,
        adjustLastFn = adjustLastFn
    )
    blocks <- keep(
        seq_len(length(allGaps) - 1),
        .swlBlockQualifies,
        allGaps = allGaps,
        minBlockFn = minBlockFn,
        ctx = ctx
    )
    windowed <- map(
        blocks,
        .swlBlockAt,
        allGaps = allGaps,
        fns = fns,
        ctx = ctx
    )
    buildSegmentResult(
        .swlField(windowed, "starts"),
        .swlField(windowed, "ends"),
        .swlField(windowed, "fillStarts"),
        .swlField(windowed, "fillEnds"),
        n,
        verbose
    )
}

# One quarter-index hop; the step index is unused, `reduce` just applies it
# `n` times.
# @noRd
.nthQuaterStep <- function(x, step, quaterIdx) {
    quaterIdx[x]
}

# Apply the quarter-distance index map `quaterIdx` `n` times to `x` (n = 1..4
# gives the 1st..4th quarter boundary from x). Used by segmentByDist.
# @noRd
.nthQuaterIdx <- function(x, n, quaterIdx) {
    reduce(seq_len(n), .nthQuaterStep, quaterIdx = quaterIdx, .init = x)
}

# --- segmentByDist helpers --------------------------------------------------

# For each SNP, the index of the last SNP within a quarter of `cutoff`,
# clamped to [1, n].
#
# The C++ original walks a second pointer forward, which is the same thing as
# counting how many positions fall strictly below each target -- that is
# `findInterval(..., left.open = TRUE)`, so no pointer has to be carried.
# Verified identical to the pointer walk over 600 random inputs, tied
# positions included.
.segByDistQuaterIdx <- function(pos, n, cutoff) {
    lastBelow <- findInterval(
        as.numeric(pos) + cutoff / 4,
        pos,
        left.open = TRUE
    )
    pmax(pmin(lastBelow, n), 1L)
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

# Window `k`'s rows, indexed globally and trimmed to that window's fill range.
# @noRd
.dentistMergeWindowAt <- function(k, dentistResultByWindow, windowDividedRes) {
    imputedK <- dentistResultByWindow[[k]]
    offset <- windowDividedRes$windowStartIdx[k] - 1
    imputedK |>
        mutate(
            index_within_window = seq_len(nrow(imputedK)),
            index_global = .data$index_within_window + offset
        ) |>
        filter(
            .data$index_global >= windowDividedRes$fillStartIdx[k] &
                .data$index_global < windowDividedRes$fillEndIdx[k]
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
    bind_rows(map(
        seq_len(nrow(windowDividedRes)),
        .dentistMergeWindowAt,
        dentistResultByWindow = dentistResultByWindow,
        windowDividedRes = windowDividedRes
    ))
}

# ## File-I/O functions (dentist_from_files, read_dentist_sumstat,
# parse_dentist_output)
### have been removed. Use the standard pipeline: load genotypes via
### loadGenotypeRegion(), compute LD via computeLd(), then call dentist()
### or ldMismatchQc() directly.

# =============================================================================
# SLALoM: Approximate Bayes Factor single-causal-variant outlier detection
# =============================================================================

# Approximate (Wakefield) Bayes factors from z-scores and standard errors:
# per-variant log-BF and normalized posterior probabilities.
# @noRd
.approxBayesFactor <- function(z, se, W = 0.04) {
    V <- se^2
    r <- W / (W + V)
    lbf <- 0.5 * (log(1 - r) + (r * z^2))
    # matrixStats is an Import and its logSumExp is the same algorithm in C;
    # pecotmr carried a hand-rolled R copy of it.
    denom <- matrixStats::logSumExp(lbf, na.rm = TRUE)
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

# Lead-variant LD column. From genotypes under the sample estimator it is read
# off X's LD factor -- crossprod(F, F[, lead]) / (n - 1), computeLd()'s column
# exactly -- without forming the p x p matrix. The population and GCTA
# estimators treat missing calls their own way, so they still build it
# through resolveLdInput(), as dentist does.
#' @importFrom checkmate assertMatrix
.slalomLeadR <- function(zScore, R, X, leadIdx, ldMethod) {
    if (is.null(R) && identical(ldMethod, "sample")) {
        factor <- .ldFactorFromDosage(as.matrix(X))
        assertMatrix(factor, ncols = length(zScore))
        lead <- crossprod(factor, factor[, leadIdx, drop = FALSE])
        # The lead's correlation with itself is 1 exactly, as on the
        # matrix's unit diagonal, not 1 give or take rounding.
        return(replace(drop(lead) / (nrow(factor) - 1), leadIdx, 1))
    }
    R <- resolveLdInput(R = R, X = X, ldMethod = ldMethod)$R
    assertMatrix(R, nrows = length(zScore), ncols = length(zScore))
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
    rawT <- (zScore - rLead * zScore[leadIdx])^2 / (1 - r2Lead)
    # A negative statistic means 1 - r2 went negative (|r| > 1 from a
    # mismatched panel); Inf routes it straight to the outlier branch.
    tDentistS <- replace(rawT, rawT < 0, Inf)
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
#' \code{X}. SLALOM needs only the lead variant's LD with the others, so with
#' \code{ldMethod = "sample"} it computes that one column from the genotypes
#' and never forms the LD matrix.
#'
#' @param zScore Numeric vector of z-scores corresponding to each variant.
#' @param R Square LD correlation matrix. Provide either \code{R} or \code{X}.
#' @param X Genotype matrix (samples x SNPs). With \code{ldMethod = "sample"}
#'   only the lead variant's LD column is computed from it; any other
#'   \code{ldMethod} computes the LD matrix via \code{computeLd(X)}.
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
#'   \code{X} is provided: \code{"sample"} (default), which uses the
#'   genotypes directly, or \code{"population"} or \code{"gcta"}, which go
#'   through \code{computeLd}. Ignored when \code{R} is provided directly.
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
    .ldRequireOneSource(R, X)
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
    rLead <- .slalomLeadR(zScore, R, X, leadIdx, ldMethod)
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
# Post-finemap credible-set update strategy
# =============================================================================

#' Process Credible Set Information and Determine Updating Strategy
#'
#' This function categorizes Credible Sets (CS) within a study block into
#' different updating strategies based on their statistical properties and
#' correlations.
#'
#' @param df Data frame. Contains information about Credible Sets for a specific
#'   study and block.
#' @param highCorrCols Character vector. Names of columns in df that represent
#'   high correlations; a non-top CS whose value in any of them is 1 is
#'   tagged. A CS's own column of \code{\link{extractCsInfo}}'s between-CS
#'   correlations (\code{cs_corr_<i>} for row \code{i}, its self-correlation)
#'   is never read, so passing every \code{cs_corr_*} column tags a CS only by
#'   its correlation with the other sets. Rows are taken in
#'   \code{extractCsInfo()}'s order.
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
#' @importFrom checkmate assertCharacter
#' @export
autoDecision <- function(df, highCorrCols) {
    assertCharacter(highCorrCols, any.missing = FALSE)
    # Identify top_cs
    topCsIndex <- which.max(abs(df$top_z))
    withTop <- mutate(df, top_cs = seq_len(nrow(df)) == topCsIndex)
    # Identify tagged_cs. `.autoDecisionTagged()` reads `top_cs`, so it has to
    # see the frame that already carries it.
    flagged <- mutate(
        withTop,
        tagged_cs = map_lgl(
            seq_len(nrow(withTop)),
            .autoDecisionTagged,
            df = withTop,
            highCorrCols = highCorrCols
        )
    )
    # Count total and remaining CS
    totalCs <- nrow(flagged)
    taggedCsCount <- sum(flagged$tagged_cs)
    remainingCs <- if (totalCs > 0) totalCs - 1 - taggedCsCount else 0
    # Determine method
    mutate(
        flagged,
        method = case_when(
            taggedCsCount == 0 & totalCs > 1 ~ "BVSR",
            (remainingCs == 0 & totalCs > 1) | (totalCs == 1) ~ "SER",
            remainingCs > 0 ~ "BCR",
            TRUE ~ NA_character_
        )
    )
}


# =============================================================================
# RAISS: regression-based sumstats imputation
# =============================================================================

# Panel rows that are unsafe to impute, over and above being already known.
#
# Both cases come from a reference panel that carries a variant and its own
# allele flip as two separate entries -- two distinct indels at one position
# that happen to be each other's flip is rare but biologically real.
#
#   * flip-of-known: the row is the allele flip of a variant the sumstats
#     already carry. It is not missing data; imputing it would invent a second
#     copy of a variant whose z-score is already in hand.
#   * flip-pair: the panel carries BOTH orientations at that position. Which
#     one an imputed value would refer to is undecidable -- the same ambiguity
#     that makes matchVariants() drop the sumstats variant there. Imputing
#     either would put the dropped variant straight back, which is what the
#     drop exists to prevent.
#
# Returns a logical over `refPanelIds`. Callers apply it to imputation
# CANDIDATES only, never to knowns: removing a known would break the
# correspondence between `knownZscores$z` and the LD rows indexed by `knowns`.
# @noRd
.raissUnsafeToImpute <- function(refPanelIds, knownIds) {
    .raissFlipPairMask(refPanelIds) |
        .raissFlipOfKnownMask(refPanelIds, knownIds)
}

# Rows whose own allele flip is also present in the panel.
#
# "Flip" is the relation matchVariants() uses, not a literal A1/A2 reversal:
# a panel that records the second orientation on the opposite strand (A:G and
# C:T) carries exactly the ambiguity one that records it plainly (A:G and G:A)
# does, and a string key sees only the second. Since the mask exists to stop
# imputation putting back a variant matchVariants() dropped, the two have to
# mean the same thing -- so the panel is matched against ITSELF and a row is
# ambiguous precisely when a sumstats variant sitting on it would not survive.
#
# Self-matching also subsumes the degenerate A1 == A2 row (its own flip) that
# the string key had to special-case: such a row matches itself exactly and
# never as a swap, so no sign conflict arises. Unparseable ids fall back to
# exact string identity, where every row matches itself and nothing is
# ambiguous.
#
# The SECOND of two identical entries comes back TRUE, which the string key
# missed. That is not ambiguity but it is still unsafe: imputing it would put
# one id into the sumstats twice, and matchVariants answers a repeated id
# once, so the next LD lookup would abort on the copy it could not place.
# @noRd
.raissFlipPairMask <- function(refPanelIds) {
    n <- length(refPanelIds)
    if (n == 0L) {
        return(logical(0))
    }
    m <- matchVariants(
        refPanelIds,
        refPanelIds,
        removeStrandAmbiguous = FALSE
    )
    !is_in(seq_len(n), m$idxA)
}

# The retained-component mask, truncated to at most `maxRank` components.
# @noRd
.svdCapRank <- function(keep, maxRank) {
    if (is.null(maxRank) || maxRank <= 0) {
        return(keep)
    }
    keepIdx <- which(keep)
    nKeep <- min(length(keepIdx), maxRank)
    if (length(keepIdx) <= nKeep) {
        return(keep)
    }
    replace(keep, keepIdx[(nKeep + 1):length(keepIdx)], FALSE)
}

# Rows that match a known sumstats variant in the opposite orientation.
# @noRd
.raissFlipOfKnownMask <- function(refPanelIds, knownIds) {
    mask <- rep(FALSE, length(refPanelIds))
    if (length(knownIds) == 0L) {
        return(mask)
    }
    m <- matchVariants(refPanelIds, knownIds, removeStrandAmbiguous = FALSE)
    replace(mask, m$idxA[m$sign < 0], TRUE)
}

# Positions the GWAS already typed.
#
# Imputation fills *un-observed* variants. A panel form at a site the GWAS
# measured -- a second alternate allele, or the opposite orientation -- is a
# redundant record at a measured position, not new information. On a large
# reference this dominates: 43 observed multi-allelic positions produced 2969
# imputed variants in one block before this guard.
#
# Keyed on (chrom, pos), which raiss() spells in UCSC style on both tables. A
# position with no GWAS observation still imputes every panel form it has,
# multi-allelic included.
# @noRd
.raissObservedPosition <- function(refPanel, knownZ) {
    n <- nrow(refPanel)
    if (is.null(refPanel$chrom) || is.null(knownZ$chrom) || n == 0L) {
        return(rep(FALSE, n))
    }
    is_in(
        .raissPositionKey(refPanel$chrom, refPanel$pos),
        .raissPositionKey(knownZ$chrom, knownZ$pos)
    )
}

# "chrom:pos" keys.
# @noRd
.raissPositionKey <- function(chrom, pos) {
    str_c(as.character(chrom), ":", as.character(pos))
}

# Imputation targets: panel rows the sumstats do not carry, minus the rows no
# imputed value could be trusted for and the rows sitting at a position the
# GWAS already typed.
# @noRd
.raissUnknownIdx <- function(refPanel, knownZ, knownsId) {
    candidate <- !is_in(refPanel$variant_id, knownsId)
    unsafe <- .raissUnsafeToImpute(refPanel$variant_id, knownZ$variant_id)
    observed <- .raissObservedPosition(refPanel, knownZ)
    which(candidate & !unsafe & !observed)
}

# Drop LD-matrix rows/cols for variants filtered out of the imputation result.
# @noRd
.raissFilterLd <- function(ldMatrix, refPanel, resultFilter) {
    filteredOutVariant <- setdiff(refPanel$variant_id, resultFilter$variant_id)
    if (length(filteredOutVariant) > 0) {
        filteredOutId <- match(filteredOutVariant, refPanel$variant_id)
        as.matrix(ldMatrix)[-filteredOutId, -filteredOutId]
    } else {
        as.matrix(ldMatrix)
    }
}

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
    unknowns <- .raissUnknownIdx(refPanel, knownZscores, knownsId)
    edge <- .raissEdgeCase(knowns, unknowns, knownZscores, verbose, ldMatrix)
    if (!is.null(edge)) {
        return(edge$value)
    }
    zt <- knownZscores$z
    sigT <- ldMatrix[knowns, knowns, drop = FALSE]
    sigIT <- ldMatrix[unknowns, knowns, drop = FALSE]
    results <- raissModel(zt, sigT, sigIT, lamb, rcond) |>
        formatRaissDf(refPanel, unknowns) |>
        filterRaissOutput(r2Threshold, minimumLd, verbose)
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
#
# The LD route pseudo-inverts the ridged LD of the known variants with
# MASS::ginv(tol = rcond), which drops every eigen-direction whose eigenvalue
# is below `rcond` times the largest. Here those eigenvalues are
# (d^2 + cReg) / (nSamples - 1), so the same directions are dropped and the
# two routes impute identically. Directions outside the genotypes' span carry
# no LD with the variants being imputed, so they change neither route. The LD
# score inverts nothing and keeps every direction.
.raissSvdImpute <- function(
    X,
    knowns,
    unknowns,
    zt,
    lamb,
    rcond,
    svdTol,
    nSamples
) {
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
    inverted <- d2PlusC > rcond * d2PlusC[[1L]]
    # w = U %*% diag(d / (d^2 + c)) %*% V' zt  (projection of zt through SVD),
    # over the directions the pseudo-inverse keeps.
    VtZt <- crossprod(V, zt)
    w <- U %*% (inverted * d / d2PlusC * VtZt)
    # Single dgemm X' %*% [w | U]: col 1 (rows unknowns) -> mu; rest -> A.
    XtWU <- crossprod(X, cbind(w, U))
    muRaw <- as.numeric(XtWU[unknowns, 1])
    A <- XtWU[unknowns, -1, drop = FALSE]
    rm(XtWU)
    # Variance and LD score in one pass over A^2.
    ASq <- A^2
    rm(A)
    dWeights <- cbind(inverted * d2 / d2PlusC, d2)
    scores <- ASq %*% dWeights
    rm(ASq)
    nm1 <- nSamples - 1
    varRaw <- (1 + lamb) - scores[, 1] / nm1
    raissLdScore <- scores[, 2] / nm1^2
    rm(scores)
    .raissResult(
        varRaw = varRaw,
        muRaw = muRaw,
        lamb = lamb,
        raissLdScore = raissLdScore,
        conditionNumber = rep(d[1] / d[length(d)], length(unknowns)),
        correctInversion = rep(TRUE, length(unknowns))
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
#' @param rcond Eigenvalue cutoff for the pseudo-inverse, as on the LD-based
#'   path: directions whose ridged eigenvalue is below \code{rcond} times the
#'   largest are dropped.
#' @param svdTol Relative tolerance below which singular values of Xt are
#'   treated as numerically zero.
#' @param r2Threshold R square threshold below which SNPs are filtered from the
#'   output.
#' @param minimumLd Minimum LD score threshold for SNP filtering.
#' @param verbose Logical indicating whether to print progress information.
#'
#' @return A list containing filtered and unfiltered results, and ldMat = NULL.
#' @importFrom dplyr arrange
#' @noRd
raissSingleMatrixFromX <- function(
    refPanel,
    knownZscores,
    X,
    lamb = 0.01,
    rcond = 0.01,
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
    unknowns <- .raissUnknownIdx(refPanel, knownZscores, knownsId)
    edge <- .raissEdgeCase(knowns, unknowns, knownZscores, verbose)
    if (!is.null(edge)) {
        return(edge$value)
    }
    imp <- .raissSvdImpute(
        X,
        knowns,
        unknowns,
        knownZscores$z,
        lamb = lamb,
        rcond = rcond,
        svdTol = svdTol,
        nSamples = nrow(X)
    )
    results <- formatRaissDf(imp, refPanel, unknowns) |>
        filterRaissOutput(r2Threshold, minimumLd, verbose)
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

# --- raiss: the genotype and LD-matrix routes --------------------------------

# The genotype route: SVD-based imputation from a standardized genotype
# matrix, which never forms the p x p LD matrix.
# @noRd
.raissGenotypePath <- function(refPanel, knownZscores, genotypeMatrix, p) {
    if (p$verbose) {
        inform("Processing genotype matrix via SVD-based imputation...")
    }
    raissSingleMatrixFromX(
        refPanel,
        knownZscores,
        genotypeMatrix,
        lamb = p$lamb,
        rcond = p$rcond,
        svdTol = p$svdTol,
        r2Threshold = p$r2Threshold,
        minimumLd = p$minimumLd,
        verbose = p$verbose
    )
}

# The LD-matrix route: imputation from the correlation matrix of refPanel's
# variants.
# @noRd
.raissLdPath <- function(refPanel, knownZscores, ldMatrix, p) {
    if (p$verbose) {
        inform("Processing LD matrix...")
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

# The LD reference a raiss() call names, if any: an LdMatrix passed as
# `ldMatrix`, or the genotype panel (or mixture of panels) passed as
# `genotypes`. Genotypes are taken only that way, never as a bare matrix:
# centring, scaling, column order and allele coding are all the caller's to get
# right in a bare matrix, and every one of them fails silently. A reference
# passed through the other argument is refused rather than read the wrong way.
# @noRd
.raissReferenceOf <- function(ldMatrix, genotypes) {
    if (.ldIsGenotypeRef(ldMatrix)) {
        abort("Pass genotypes as `genotypes`, not `ldMatrix`.")
    }
    if (!is.null(genotypes) && !.ldIsGenotypeRef(genotypes)) {
        abort(c(
            "`genotypes` must be a genotype panel or a GenotypeMixture.",
            i = "Read one with readGenotypes(), or pass LD as `ldMatrix`."
        ))
    }
    if (.ldIsMatrixRef(ldMatrix) || !is.null(genotypes)) {
        return(ldMatrix %||% genotypes)
    }
    if (!is.null(ldMatrix) && !is.matrix(ldMatrix)) {
        abort("`ldMatrix` must be a square LD matrix or an LdMatrix.")
    }
    NULL
}

# An LD reference brings its own variant table, which refPanel defaults to.
# Known and imputed rows are merged on chrom, so known z-scores on none of the
# reference's chromosomes have nothing to impute from; say so here rather
# than in an empty join.
# @noRd
.raissRefPanelOf <- function(ref, knownZscores) {
    panel <- .refVariantsFromLdRef(ref)
    if (!any(is_in(knownZscores[["chrom"]], panel$chrom))) {
        abort(c(
            "`knownZscores` shares no chromosome with the LD reference.",
            i = glue(
                "The reference covers ",
                "{str_flatten(unique(panel$chrom), ', ')}."
            )
        ))
    }
    panel
}

# knownZscores as raiss() reads it: a table as given, or a summary-statistics
# entry -- the GRanges sumStats() returns -- as the table summaryStatsQc()
# builds from one.
# @noRd
.raissKnownTable <- function(knownZscores) {
    if (!methods::is(knownZscores, "GRanges")) {
        return(knownZscores)
    }
    need <- c("SNP", "A1", "A2", "Z")
    absent <- setdiff(need, colnames(mcols(knownZscores)))
    if (length(absent) > 0L) {
        abort(c(
            glue(
                "`knownZscores` as a GRanges lacks column(s) ",
                "{str_flatten(absent, ', ')}."
            ),
            i = "Pass a summary-statistics entry, as sumStats(x) returns it."
        ))
    }
    .qcRaissKnownZ(.entryGrangesToDf(knownZscores))
}

# A table handed to raiss() with its `chrom` column in UCSC style, so its rows
# meet the reference's and each other's on equal labels. NULL stays NULL.
# @noRd
.raissUcscChrom <- function(df) {
    if (is.null(df) || is.null(df[["chrom"]])) {
        return(df)
    }
    `[[<-`(df, "chrom", value = withChrPrefix(df[["chrom"]]))
}

# The genotypes' LD factor for refPanel's variants (see .ldRefFactor()), the
# way summaryStatsQc() hands it to the SVD route, re-signed into refPanel's
# allele coding.
# @noRd
.raissScaledDosage <- function(genotypes, refPanel) {
    m <- .ldFromRefMatch(genotypes, refPanel$variant_id, "raiss", "error")
    sweep(.ldRefFactor(genotypes, m$idx), 2L, m$sign, `*`)
}

#' Impute Summary Statistics Using LD (RAISS)
#'
#' An R port of the RAISS summary-statistic imputation library, described in
#' Julienne H, Shi H, Pasaniuc B, Aschard H (2019), "RAISS: robust and accurate
#' imputation from summary statistics", Bioinformatics 35(22):4837-4839.
#' \doi{10.1093/bioinformatics/btz466}
#'
#' The imputation model it implements is the one introduced in Pasaniuc B,
#' Zaitlen N, et al. (2014), "Fast and accurate imputation of summary
#' statistics enhances evidence of functional enrichment", Bioinformatics
#' 30(20):2906-2914. \doi{10.1093/bioinformatics/btu416}
#'
#' One call imputes one region, from either of the two kinds of LD reference.
#' The LD can be given as a square correlation matrix or an
#' \code{\link{LdMatrix}}, and the imputation then inverts it. Or a genotype
#' panel from \code{\link{readGenotypes}} can be given instead, and the
#' imputation then works from the singular value decomposition of its
#' standardized dosages, never forming the p x p LD matrix (memory and
#' compute savings when n << p).
#'
#' @param refPanel A data frame containing 'chrom', 'pos', 'variant_id', 'A1',
#'   and 'A2'. Optional when the LD comes as an \code{LdMatrix} or a genotype
#'   panel, whose own variants are used. Required with a plain
#'   \code{ldMatrix}, whose rows and columns follow it.
#' @param knownZscores A data frame containing 'chrom', 'pos', 'variant_id',
#'   'A1', 'A2', and 'z' values, or a summary-statistics \code{GRanges} --
#'   one entry as \code{\link{sumStats}} returns it, with its \code{SNP},
#'   \code{A1}, \code{A2} and \code{Z} columns. Here and in \code{refPanel}
#'   a chromosome may be labelled in any style (\code{22}, \code{"22"},
#'   \code{"chr22"}); both are read in UCSC style (\code{"chr22"}), which
#'   the result uses.
#' @param ldMatrix A square LD matrix over \code{refPanel}'s variants, in its
#'   order, or an \code{LdMatrix} covering them. Provide either
#'   \code{ldMatrix} or \code{genotypes}, not both.
#' @param genotypes A genotype panel (see \code{\link{readGenotypes}}) or a
#'   \code{\link{GenotypeMixture}} covering the variants, as an alternative
#'   to \code{ldMatrix}. Its dosages are matched to \code{refPanel} by
#'   variant and allele, mean-imputed and standardized here.
#' @param lamb Regularization term added to the diagonal of the ldMatrix.
#' @param rcond Eigenvalue cutoff for the pseudo-inverse of the known
#'   variants' LD: directions whose ridged eigenvalue is below \code{rcond}
#'   times the largest are dropped. Applied on both routes, so a genotype
#'   panel and the LD matrix computed from it give the same imputation.
#' @param svdTol Relative tolerance below which singular values of the
#'   standardized genotypes are treated as numerically zero (genotype route
#'   only).
#' @param r2Threshold R square threshold below which SNPs are filtered from the
#'   output.
#' @param minimumLd Minimum LD score threshold for SNP filtering.
#' @param verbose Logical indicating whether to print progress information.
#'
#' @return A list containing filtered and unfiltered results, and filtered LD
#'   matrix (ldMat is NULL on the genotype route).
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
    refPanel = NULL,
    knownZscores,
    ldMatrix = NULL,
    genotypes = NULL,
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
    if (!is.null(genotypes) && !is.null(ldMatrix)) {
        abort("Provide either ldMatrix or genotypes, not both.")
    }
    if (is.null(genotypes) && is.null(ldMatrix)) {
        abort("Provide either ldMatrix or genotypes.")
    }
    knownZscores <- .raissUcscChrom(.raissKnownTable(knownZscores))
    refPanel <- .raissUcscChrom(refPanel)
    ref <- .raissReferenceOf(ldMatrix, genotypes)
    if (!is.null(ref)) {
        refPanel <- refPanel %||% .raissRefPanelOf(ref, knownZscores)
    } else if (is.null(refPanel)) {
        abort(str_c(
            "`refPanel` is required unless the LD comes as an LdMatrix or a ",
            "genotype panel."
        ))
    }
    if (!is.null(genotypes)) {
        X <- .raissScaledDosage(genotypes, refPanel)
        return(.raissGenotypePath(refPanel, knownZscores, X, p))
    }
    R <- if (is.null(ref)) {
        ldMatrix
    } else {
        .ldFromRef(ref, refPanel$variant_id, label = "raiss")
    }
    .raissLdPath(refPanel, knownZscores, R, p)
}

#' @rdname RaissParam
#' @aliases RaissParam-class
#' @exportClass RaissParam
setClass(
    "RaissParam",
    contains = "MethodParam",
    slots = c(
        lamb = "numeric",
        rcond = "numeric",
        svdTol = "numeric",
        r2Threshold = "numeric",
        minimumLd = "numeric",
        mafCutoff = "numeric",
        macCutoff = "numeric",
        imissCutoff = "numeric",
        flank = "numeric"
    )
)

#' @title Arguments For RAISS Imputation During Summary-Statistic QC
#' @description The settings \code{\link{summaryStatsQc}} uses when
#'   \code{impute = TRUE}. Five are forwarded to \code{\link{raiss}}; the
#'   other four are pecotmr's own, bounding which panel variants RAISS is
#'   asked to impute and how wide a window it reads.
#'
#'   \code{rcond} applies whatever the collection's LD reference is: QC hands
#'   RAISS a genotype panel's dosages, which route through the SVD solver, or
#'   an \code{\link{LdMatrix}}'s correlations, which are inverted directly,
#'   and both drop the same weak directions. \code{svdTol} only matters on
#'   the genotype route.
#' @param lamb Regularization term added to the LD diagonal. Default
#'   \code{0.01}.
#' @param rcond Eigenvalue cutoff for the pseudo-inverse of the known
#'   variants' LD: directions whose ridged eigenvalue is below \code{rcond}
#'   times the largest are dropped. Default \code{0.01}, as
#'   \code{\link{raiss}}'s own.
#' @param svdTol Relative tolerance below which singular values of the
#'   standardized genotypes are treated as numerically zero, used when the LD
#'   reference is a genotype panel. Default \code{1e-12}. Note this is
#'   tighter than \code{\link{raiss}}'s own \code{1e-8} default.
#' @param r2Threshold Imputation \eqn{R^2} below which a variant is dropped
#'   from the output. Default \code{0.6}.
#' @param minimumLd Minimum LD score for retaining an imputed variant.
#'   Default \code{5}.
#' @param mafCutoff Minor-allele-frequency floor for a variant to be an
#'   imputation TARGET. Default \code{0}.
#' @param macCutoff Minor-allele-count floor for a target; the stricter of
#'   this and \code{mafCutoff} applies. Default \code{0}.
#' @param imissCutoff Missingness ceiling for a target. Default \code{1}.
#' @param flank Base pairs by which to widen the analysis-region window on
#'   each side, retaining LD context for edge variants. The window never
#'   widens into another row of the same study on its chromosome -- the
#'   neighbouring LD block of a block-split GWAS -- since that row imputes
#'   its own variants, nor past a row's block range where that was given
#'   (\code{regionGiven}, see \code{\linkS4class{GwasSumStats}}). Default
#'   \code{0}.
#' @return A \code{RaissParam} object, a \code{\link{MethodParam}}.
#' @seealso \code{\link{summaryStatsQc}}, \code{\link{raiss}}
#' @examples
#' RaissParam(mafCutoff = 0.01, flank = 5e5)
#' @export
RaissParam <- function(
    lamb = 0.01,
    rcond = 0.01,
    svdTol = 1e-12,
    r2Threshold = 0.6,
    minimumLd = 5,
    mafCutoff = 0,
    macCutoff = 0,
    imissCutoff = 1,
    flank = 0
) {
    new(
        "RaissParam",
        lamb = lamb,
        rcond = rcond,
        svdTol = svdTol,
        r2Threshold = r2Threshold,
        minimumLd = minimumLd,
        mafCutoff = mafCutoff,
        macCutoff = macCutoff,
        imissCutoff = imissCutoff,
        flank = flank
    )
}

# The RAISS bundle as a plain list with every field present. A constructor
# result already carries all nine -- they are all formals with defaults --
# and the empty list that means "no options" becomes the defaults outright,
# so the read sites need no per-field fallback.
# @noRd
.ssqcResolveImputeArgs <- function(raissParam) {
    if (length(raissParam) == 0L) {
        return(as.list(RaissParam()))
    }
    as.list(raissParam)
}

#' @importFrom checkmate assertNumeric
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
    # Checked before sigT is used: the inversion retries on any error, so a
    # bad input would otherwise surface as base R's error from its third try.
    # sigT is the known x known LD and sigIT the unknown x known LD, so both
    # carry one column per known z-score.
    assertNumeric(zt, min.len = 1L)
    assertNumeric(sigT)
    assertNumeric(sigIT)
    k <- length(zt)
    assertMatrix(sigT, nrows = k, ncols = k)
    assertMatrix(sigIT, min.rows = 1L, ncols = k)
    sigTInv <- invertMatRecursive(sigT, lamb, rcond)
    checks <- .raissInversionChecks(
        sigT,
        sigTInv,
        lamb,
        nOut = if (batch) nrow(sigIT) else 1L,
        reportConditionNumber = reportConditionNumber
    )
    varRaissLdScore <- computeVar(sigIT, sigTInv, lamb, batch)
    .raissResult(
        varRaw = varRaissLdScore$var,
        muRaw = computeMu(sigIT, sigTInv, zt),
        lamb = lamb,
        raissLdScore = varRaissLdScore$raissLdScore,
        conditionNumber = checks$conditionNumber,
        correctInversion = checks$correctInversion
    )
}

# The inversion's diagnostics, repeated for each of the `nOut` variants a
# batch imputes: sigT's condition number when asked for (NA otherwise), and
# whether sigTInv inverts what was inverted. That is the ridged matrix
# (diagonal 1 + lamb), so it is what the inverse is checked against; against
# the raw sigT the check could only fail, off by O(lamb). The original RAiSS
# implementation checks the ridged matrix too: its fill_diagonal() writes the
# ridge into sig_t in place.
# @noRd
.raissInversionChecks <- function(
    sigT,
    sigTInv,
    lamb,
    nOut,
    reportConditionNumber
) {
    conditionNumber <- if (reportConditionNumber) {
        rep(kappa(sigT, exact = TRUE, norm = "2"), nOut)
    } else {
        NA
    }
    sigTRidged <- `diag<-`(sigT, 1 + lamb)
    list(
        conditionNumber = conditionNumber,
        correctInversion = rep(checkInversion(sigTRidged, sigTInv), nOut)
    )
}

# An imputation's result from its raw variances and means, as both RAiSS
# routes (raissModel(), .raissSvdImpute()) report it: the variances bounded
# by varInBoundaries(), and each mean scaled by sqrt(R2), with R2 = (1 + lamb)
# minus that variance.
# @noRd
.raissResult <- function(
    varRaw,
    muRaw,
    lamb,
    raissLdScore,
    conditionNumber,
    correctInversion
) {
    varNorm <- varInBoundaries(varRaw, lamb)
    R2 <- (1 + lamb) - varNorm
    list(
        var = varNorm,
        mu = muRaw / sqrt(R2),
        raissLdScore = raissLdScore,
        conditionNumber = conditionNumber,
        correctInversion = correctInversion
    )
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

    resolved <- mutate(
        mergedDf,
        # A known variant was not imputed, so it carries no imputation
        # quality: Var = -1 and an infinite LD score mark it as observed.
        Var = replace(.data$Var, fromKnown, -1),
        raissLdScore = replace(.data$raissLdScore, fromKnown, Inf),
        # Overlapping z columns resolved: z from knownZscores where available,
        # otherwise z from raissDf.
        z = if_else(fromKnown, .data$z.y, .data$z.x)
    ) |>
        # Remove the extra columns produced by the join (z.x, z.y).
        select(-all_of(c("z.x", "z.y"))) |>
        arrange(.data$pos)
    # Imputed variants' beta / se are NA to avoid confusion, since they are not
    # imputed. Both are optional (knownZscores may omit them), so guard on
    # column presence explicitly rather than relying on a data.frame silently
    # creating an all-NA column on `$col[mask] <- NA`.
    observed <- resolved$Var == -1
    mutate(
        resolved,
        !!!compact(list(
            beta = if (is_in("beta", colnames(resolved))) {
                replace(resolved$beta, observed, NA)
            },
            se = if (is_in("se", colnames(resolved))) {
                replace(resolved$se, observed, NA)
            }
        ))
    )
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
    selected <- select(
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
    zscoresNofilter <- mutate(selected, raissR2 = 1 - .data$Var)
    kept <- filter(
        zscoresNofilter,
        .data$raissR2 > r2Threshold & .data$raissLdScore >= minimumLd
    )
    if (verbose) {
        .filterRaissReport(zscoresNofilter, kept, r2Threshold, minimumLd)
    }
    list(zscoresNofilter = zscoresNofilter, zscores = kept)
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

# TRUE when sigTInv is a pseudo-inverse of sigT (sigT sigTInv sigT == sigT),
# else FALSE -- never all.equal()'s text describing the difference.
checkInversion <- function(sigT, sigTInv) {
    isTRUE(all.equal(sigT, sigT %*% (sigTInv %*% sigT), tolerance = 1e-5))
}

varInBoundaries <- function(var, lamb) {
    floored <- replace(var, var < 0, 0)
    replace(floored, floored > (0.99999 + lamb), 1)
}

#' @importFrom rlang try_fetch
invertMat <- function(mat, lamb, rcond) {
    try_fetch(
        {
            # Modify the diagonal elements of mat
            diag(mat) <- 1 + lamb
            # Compute the pseudo-inverse
            matInv <- ginv(mat, tol = rcond)
            return(matInv)
        },
        error = function(cnd) {
            # Second attempt with updated lamb and rcond in case of an error
            diag(mat) <- 1 + lamb * 1.1
            matInv <- ginv(mat, tol = rcond * 1.1)
            return(matInv)
        }
    )
}

invertMatRecursive <- function(mat, lamb, rcond) {
    try_fetch(
        {
            # Modify the diagonal elements of mat
            diag(mat) <- 1 + lamb
            # Compute the pseudo-inverse
            matInv <- ginv(mat, tol = rcond)
            return(matInv)
        },
        error = function(cnd) {
            # Recursive call with updated lamb and rcond in case of an error
            invertMat(mat, lamb * 1.1, rcond * 1.1)
        }
    )
}

# =============================================================================
# Top-level summaryStatsQc() pipeline + helpers
# =============================================================================

#' @rdname SlalomParam
#' @aliases SlalomParam-class
#' @exportClass SlalomParam
setClass(
    "SlalomParam",
    contains = "MethodParam",
    slots = c(
        standardError = "numeric_OR_NULL",
        abfPriorVariance = "numeric",
        nlog10pDentistSThreshold = "numeric",
        r2Threshold = "numeric",
        leadVariantChoice = "character"
    )
)

#' @title Arguments For The SLALOM LD-Mismatch Check
#' @description Options for \code{\link{slalom}}, pecotmr's SLALOM
#'   implementation. Every field is a formal of that function, so an unknown
#'   name is rejected by R as an unused argument; there is no \code{...}
#'   because the contract is complete. Pass the result as
#'   \code{ldMismatchQc(method = )} or
#'   \code{summaryStatsQc(ldMismatchQcMethod = )}.
#' @param standardError Per-variant standard errors. \code{NULL} (default)
#'   leaves \code{\link{slalom}}'s own data-dependent default in place.
#' @param abfPriorVariance Prior variance for the approximate Bayes factor.
#'   Default \code{0.04}.
#' @param nlog10pDentistSThreshold \code{-log10(p)} threshold on the
#'   DENTIST-S statistic above which a variant is an outlier. Default
#'   \code{4.0}.
#' @param r2Threshold Minimum r-squared with the lead variant for a variant to
#'   be tested. Default \code{0.6}.
#' @param leadVariantChoice How the lead variant is chosen. Default
#'   \code{"pvalue"}.
#' @return A \code{SlalomParam} object, a \code{\link{MethodParam}}.
#' @examples
#' SlalomParam(r2Threshold = 0.8)
#' @export
SlalomParam <- function(
    standardError = NULL,
    abfPriorVariance = 0.04,
    nlog10pDentistSThreshold = 4.0,
    r2Threshold = 0.6,
    leadVariantChoice = "pvalue"
) {
    new(
        "SlalomParam",
        standardError = standardError,
        abfPriorVariance = abfPriorVariance,
        nlog10pDentistSThreshold = nlog10pDentistSThreshold,
        r2Threshold = r2Threshold,
        leadVariantChoice = leadVariantChoice
    )
}

#' @rdname DentistParam
#' @aliases DentistParam-class
#' @exportClass DentistParam
setClass(
    "DentistParam",
    contains = "MethodParam",
    slots = c(
        pValueThreshold = "numeric",
        propSVD = "numeric",
        gcControl = "logical",
        nIter = "numeric",
        gPvalueThreshold = "numeric",
        duprThreshold = "numeric",
        BPPARAM = "BiocParallelParam",
        correctChenEtAlBug = "logical",
        seed = "numeric_OR_NULL"
    )
)

#' @title Arguments For The DENTIST LD-Mismatch Check
#' @description Options for \code{\link{dentistSingleWindow}}, pecotmr's
#'   DENTIST implementation. Every field is a formal of that function, so an
#'   unknown name is rejected by R as an unused argument; there is no
#'   \code{...} because the contract is complete. Pass the result as
#'   \code{ldMismatchQc(method = )} or
#'   \code{summaryStatsQc(ldMismatchQcMethod = )}.
#' @param pValueThreshold Significance threshold for the DENTIST statistic.
#'   Default \code{5e-8}.
#' @param propSVD Proportion of the singular value decomposition retained.
#'   Default \code{0.4}.
#' @param gcControl Logical; apply genomic control. Default \code{FALSE}.
#' @param nIter Number of DENTIST iterations. Default \code{10}.
#' @param gPvalueThreshold Grouping p-value threshold. Default \code{0.05}.
#' @param duprThreshold Duplicate r-squared threshold. Default \code{0.99}.
#' @param BPPARAM A \code{BiocParallelParam} whose workers number the
#'   threads DENTIST's OpenMP loops use when they run on this machine
#'   (\code{SerialParam}, \code{MulticoreParam}); any other backend gives one
#'   thread. Default \code{SerialParam()}.
#' @param correctChenEtAlBug Logical; apply the Chen et al. correction.
#'   Default \code{TRUE}.
#' @param seed Integer or \code{NULL}. Random seed for the iteration.
#' @return A \code{DentistParam} object, a \code{\link{MethodParam}}.
#' @examples
#' DentistParam(propSVD = 0.5)
#' @export
DentistParam <- function(
    pValueThreshold = 5e-8,
    propSVD = 0.4,
    gcControl = FALSE,
    nIter = 10,
    gPvalueThreshold = 0.05,
    duprThreshold = 0.99,
    BPPARAM = SerialParam(),
    correctChenEtAlBug = TRUE,
    seed = NULL
) {
    new(
        "DentistParam",
        pValueThreshold = pValueThreshold,
        propSVD = propSVD,
        gcControl = gcControl,
        nIter = nIter,
        gPvalueThreshold = gPvalueThreshold,
        duprThreshold = duprThreshold,
        BPPARAM = .assertBpparam(BPPARAM),
        correctChenEtAlBug = correctChenEtAlBug,
        seed = seed
    )
}

# The engines ldMismatchQc dispatches to, and the constructor for each. Keyed
# by the name `method` accepts, so a bare name resolves to an empty record.
# @noRd
.ldMismatchCtors <- function() {
    list(slalom = SlalomParam, dentist = DentistParam)
}

# Normalize summaryStatsQc's ldMismatchQcMethod: a name (including "none", for
# which there is no engine to configure) or the engine's constructor. Returns
# the value unchanged once validated, so it can be handed straight to
# ldMismatchQc(method = ).
# @noRd
.resolveLdMismatchChoice <- function(x) {
    allowed <- c("none", names(.ldMismatchCtors()))
    engine <- .engineOf(if (is.character(x)) x[[1L]] else x)
    if (is.null(engine) || !is_in(engine, allowed)) {
        .resolveEngineChoice(
            if (is.character(x)) x[[1L]] else x,
            allowed,
            "ldMismatchQcMethod"
        )
    }
    if (is.character(x)) x[[1L]] else x
}

#' Detect LD-Summary Statistic Mismatches
#'
#' Unified wrapper for detecting outlier variants due to LD-summary statistic
#' mismatches. Dispatches to either \code{\link{dentistSingleWindow}} or
#' \code{\link{slalom}} based on the \code{method} argument.
#'
#' @param zScore Numeric vector of z-scores.
#' @param R Square LD correlation matrix. Provide either \code{R} or \code{X}.
#' @param X Genotype matrix (samples x SNPs). With \code{ldMethod = "sample"}
#'   both methods work on the genotypes without forming the LD matrix; any
#'   other \code{ldMethod} computes LD via \code{computeLd}. \code{nSample}
#'   defaults to \code{nrow(X)}.
#' @param nSample Number of samples in the LD reference panel. Required when
#'   \code{R} is provided and \code{method = "dentist"}; inferred from \code{X}
#'   when \code{X} is provided.
#' @param method Which QC method to run: \code{"slalom"} (default) or
#'   \code{"dentist"}, or the matching constructor --
#'   \code{\link{SlalomParam}} / \code{\link{DentistParam}} -- to configure
#'   it at the same time. The constructor carries the choice, so there is no
#'   separate options argument that could disagree with it.
#' @param ldMethod Character string specifying the LD computation method when
#'   \code{X} is provided. One of \code{"sample"} (default),
#'   \code{"population"}, or \code{"gcta"}. Ignored when \code{R} is provided
#'   directly.
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
    ldMethod = "sample"
) {
    # `method` carries both the choice and its options: a bare name resolves
    # to that engine with no options set, a constructor result to that engine
    # configured. One value, so the two cannot disagree.
    chosen <- .resolveEngineChoice(
        if (is.character(method)) method[[1L]] else method,
        names(.ldMismatchCtors()),
        "method"
    )
    engineParam <- chosen$args
    if (chosen$engine == "dentist") {
        callArgs <- list_modify(
            list(
                zScore,
                R = R,
                X = X,
                nSample = nSample,
                ldMethod = ldMethod
            ),
            !!!engineParam
        )
        return(exec(dentistSingleWindow, !!!callArgs))
    } else {
        callArgs <- list_modify(
            list(zScore, R = R, X = X, ldMethod = ldMethod),
            !!!engineParam
        )
        qcResults <- exec(slalom, !!!callArgs)
        # Standardize output: slalom uses "outliers", rename to "outlier" for
        # consistency
        raw <- qcResults$data
        renameOutliers <- is_in("outliers", colnames(raw)) &&
            !is_in("outlier", colnames(raw))
        result <- if (renameOutliers) {
            rename(raw, outlier = "outliers")
        } else {
            raw
        }
        return(result)
    }
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
    replace(
        4 / (1 / nCase + 1 / nControl),
        is.na(nCase) | is.na(nControl) | nCase <= 0 | nControl <= 0,
        NA_real_
    )
}

# =============================================================================
# summaryStatsQc -- SumStats-input QC pipeline (replaces the previous
# data.frame- and QcResult-based summaryStatsQc and rssBasicQc).
# =============================================================================

# Convert one entry's GRanges into a flat tibble with the column shape
# harmonizeAlleles expects (lower-case chrom/pos plus the CapsCase mcols).
.entryGrangesToDf <- function(gr) {
    mc <- as.data.frame(S4Vectors::mcols(gr), stringsAsFactors = FALSE)
    out <- tibble(
        chrom = as.character(GenomicRanges::seqnames(gr)),
        pos = GenomicRanges::start(gr)
    )
    bind_cols(out, mc)
}

# Build a refVariants data.frame (chrom, pos, A1, A2, variant_id) from the
# panel's own variant ranges so harmonizeAlleles can join by (chrom, pos).
# The id comes from `.ldRefMatchIds()`, not the raw SNP label: a panel entry
# that spells a tag where an allele belongs still harmonizes off its A1/A2
# columns, and the repaired id is what the LD lookups downstream can find.
.refVariantsFromLdRef <- function(ldReference) {
    gr <- .ldRefRanges(ldReference)
    mc <- S4Vectors::mcols(gr)
    data.frame(
        chrom = .ldRefChrom(ldReference),
        pos = as.integer(GenomicRanges::start(gr)),
        A1 = as.character(mc$A1),
        A2 = as.character(mc$A2),
        variant_id = .ldRefMatchIds(ldReference),
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
    df <- mutate(
        df,
        !!!compact(list(
            SNP = if (
                is_in("variant_id", colnames(df)) &&
                    !is_in("SNP", colnames(df))
            ) {
                df$variant_id
            }
        ))
    )
    baseCols <- c("SNP", "A1", "A2", "Z", "N")
    # AF = directional effect-allele frequency (exported as af); MAF =
    # directionless QC frequency. Both carried when the loader resolved them.
    optCols <- c(
        "AF",
        "MAF",
        "INFO",
        "BETA",
        "SE",
        "P",
        "N_CASE",
        "N_CONTROL",
        # Lower case, unlike its neighbours, because it is not a statistic:
        # it is the per-variant flag saying whether `Z` was supplied by the
        # source or derived by pecotmr from BETA / SE, and GWAS-VCF's `EZ`
        # means only the former. Same name as the fine-mapping marginal
        # view's, so the vocabulary is one.
        "z_supplied"
    )
    use <- intersect(c(baseCols, optCols), colnames(df))
    S4Vectors::`mcols<-`(
        gr,
        value = S4Vectors::DataFrame(select(df, all_of(use)))
    )
}

# -----------------------------------------------------------------------------
# Shared entry-to-sumstat data.frame converter
# -----------------------------------------------------------------------------

# Base (variant_id/chrom/pos/A1/A2) frame for one sumstat-entry GRanges.
# @noRd
.entryDfBase <- function(gr, mc, chr) {
    n <- length(gr)
    tibble(
        variant_id = .entryDfCharColumn("SNP", mc, n),
        chrom = chr,
        pos = unname(as.integer(GenomicRanges::start(gr))),
        A1 = .entryDfCharColumn("A1", mc, n),
        A2 = .entryDfCharColumn("A2", mc, n)
    )
}

# One mcols column as character, or NA for each of the `n` variants when the
# entry does not carry it.
# @noRd
.entryDfCharColumn <- function(name, mc, n) {
    if (is_in(name, colnames(mc))) {
        as.character(mc[[name]])
    } else {
        rep(NA_character_, n)
    }
}

# One mcols column as numeric.
# @noRd
.entryStatColumn <- function(src, mc) {
    as.numeric(mc[[src]])
}

# Append the optional numeric stat columns present on the entry's mcols.
# @noRd
.entryDfAddStats <- function(df, mc) {
    statMap <- c(
        z = "Z",
        beta = "BETA",
        se = "SE",
        N = "N",
        maf = "MAF",
        af = "AF"
    )
    present <- statMap[is_in(statMap, colnames(mc))]
    # mutate() overwrites an existing column in place and appends a new one,
    # which is what the `df[[out]] <-` loop did.
    out <- mutate(
        df,
        !!!set_names(
            map(unname(present), .entryStatColumn, mc = mc),
            names(present)
        )
    )
    .entryDfDeriveMaf(out)
}

# AF is the DIRECTIONAL effect-allele frequency; MAF is its directionless
# form. A study that declared only `af:` therefore already carries the
# information, so derive the MAF rather than leaving the column absent --
# absent, the MR Wald ratio falls back to the z scale (see
# .cipGwasHasScale) even though a real frequency was supplied. Same
# pmin(af, 1 - af) the QC frequency filter uses.
#
# A declared MAF wins: it is what the study actually asserted, and AF is
# only a route to the same quantity.
# @noRd
.entryDfDeriveMaf <- function(df) {
    if (!is.null(df[["maf"]]) || is.null(df[["af"]])) {
        return(df)
    }
    mutate(df, maf = pmin(.data$af, 1 - .data$af))
}

# Convert one sumstat-entry GRanges into a flat data.frame with
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
# @noRd
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
    seqChr <- as.character(GenomicRanges::seqnames(gr))
    chr <- if (keepChrPrefix) {
        seqChr
    } else {
        str_remove(seqChr, regex("^chr", ignore_case = TRUE))
    }
    df <- .entryDfBase(gr, mc, chr) |> .entryDfAddStats(mc)
    deriveZ <- derive == "zFromBetaSe" &&
        is.null(df[["z"]]) &&
        !is.null(df[["beta"]]) &&
        !is.null(df[["se"]])
    if (!deriveZ) {
        return(df)
    }
    mutate(df, z = .data$beta / .data$se)
}

# One tuple of a GwasSumStats / QtlSumStats collection as the flat table
# .entryToSumstatDf() builds. `studyName` / `context` / `trait` select the
# tuple as sumStats() does; a GWAS takes `studyName` only.
# @noRd
.sumStatsDf <- function(
    x,
    studyName = NULL,
    context = NULL,
    trait = NULL,
    require = character(0),
    derive = c("none", "zFromBetaSe")
) {
    .entryToSumstatDf(
        sumStats(x, studyName = studyName, context = context, trait = trait),
        require = require,
        derive = derive,
        label = .sumStatsDfLabel(x, studyName, context, trait)
    )
}

# The tuple a .sumStatsDf() error names; "<auto>" marks a selector left to
# resolve on its own.
# @noRd
.sumStatsDfLabel <- function(x, studyName, context, trait) {
    parts <- c(
        studyName %||% "<auto>",
        context %||% "<auto>",
        trait %||% "<auto>"
    )
    if (is(x, "GwasSumStats")) {
        parts <- parts[[1L]]
    }
    glue("{class(x)[[1L]]}[{str_flatten(parts, '/')}]")
}

# Derive BETA and SE columns from signed Z when the entry has only Z.
# Formula (Zhu et al. 2016 / RAISS):
#   se   = 1 / sqrt(2 * maf * (1 - maf) * (N + z^2))
#   beta = z * se
# Requires Z, MAF, and N to all be present in `df`. No-op if BETA and SE
# are already there, or if any required column is missing. Returns:
#   list(df = <data.frame>, audit = NULL | list(nDerived = <int>))
# .zToPvalue (two-tailed normal p-value from a signed Z) is defined once in
# pvalCombineWrappers.R and shared package-wide.

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
    aboveTol <- if (tol > 0 && length(d) > 0) {
        out <- d / d[1] > tol
        if (!any(out)) {
            abort("All singular values are below the tolerance threshold.")
        }
        out
    } else {
        rep(TRUE, length(d))
    }
    keep <- .svdCapRank(aboveTol, maxRank)
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
    walk <- .dupWalk(z, rThreshold, list(R = ld))
    kept <- walk$dupBearer == -1
    c(
        list(
            filteredZ = z[kept],
            filteredLD = ld[kept, kept, drop = FALSE]
        ),
        walk
    )
}

# .findDuplicateVariants() from a factor of the LD: the same walk, its LD rows
# computed a stripe at a time from the factor, so the p x p matrix is never
# formed; the duplicates are dropped from the factor's columns.
# @noRd
.findDuplicateVariantsFactor <- function(z, factor, rThreshold) {
    walk <- .dupWalk(z, rThreshold, list(F = factor))
    kept <- walk$dupBearer == -1
    c(
        list(
            filteredZ = z[kept],
            filteredFactor = factor[, kept, drop = FALSE]
        ),
        walk
    )
}

# The duplicate walk: left to right, each variant not yet marked a duplicate
# anchors every later unmarked variant whose |r| with it exceeds rThreshold.
# Its LD rows are read from the LD source `ld` (.ldSourceRows()) a stripe of
# anchors at a time, so a factor never has to supply more than a stripe.
# @noRd
.dupWalk <- function(z, rThreshold, ld) {
    p <- length(z)
    dupBearer <- rep(-1, p)
    corABS <- rep(0, p)
    sign <- rep(1, p)
    count <- 1L
    minValue <- 1
    blockSize <- .ldStripeSize()
    block <- NULL
    blockStart <- 0L
    for (i in seq_len(p - 1L)) {
        if (dupBearer[i] != -1) {
            next
        }
        if (is.null(block) || i >= blockStart + nrow(block)) {
            blockStart <- i
            block <- .ldSourceRows(ld, i:min(i + blockSize - 1L, p - 1L))
        }
        row <- block[i - blockStart + 1L, ]
        idx <- (i + 1L):p
        corVec <- abs(row[idx])
        dupIdx <- which(dupBearer[idx] == -1 & corVec > rThreshold)
        if (length(dupIdx) > 0) {
            j <- idx[dupIdx]
            sign[j] <- if_else(row[j] < 0, -1, sign[j])
            corABS[j] <- corVec[dupIdx]
            dupBearer[j] <- count
        }
        minValue <- min(minValue, min(corVec))
        count <- count + 1L
    }
    list(
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
    derived <- mutate(
        df,
        !!!compact(list(
            BETA = if (!hasBeta) beta,
            SE = if (!hasSe) se
        ))
    )
    list(df = derived, audit = list(nDerived = sum(!is.na(se))))
}

# Drop variants whose (chrom, pos) overlaps any user-supplied skipRegion, a
# GRanges -- as every window given in R is.
# Normalise skipRegion to a chrom/start/end frame.
# @noRd
.parseSkipRegion <- function(skipRegion) {
    skipRegion <- .assertRegionGRanges(skipRegion, "`skipRegion`")
    tibble(
        chrom = as.character(GenomicRanges::seqnames(skipRegion)),
        start = GenomicRanges::start(skipRegion),
        end = GenomicRanges::end(skipRegion)
    )
}

# Rows falling inside skip region `i`.
# @noRd
.skipRegionMask <- function(i, parsed, dfChr, pos) {
    dfChr == parsed$chrom[i] & pos >= parsed$start[i] & pos <= parsed$end[i]
}

.applySkipRegion <- function(df, skipRegion) {
    if (is.null(skipRegion) || length(skipRegion) == 0L) {
        return(df)
    }
    parsed <- .parseSkipRegion(skipRegion)
    dfChr <- as.character(df$chrom)
    dropMask <- reduce(
        map(
            seq_len(nrow(parsed)),
            .skipRegionMask,
            parsed = parsed,
            dfChr = dfChr,
            pos = df$pos
        ),
        `|`,
        .init = rep(FALSE, nrow(df))
    )
    filter(df, !dropMask)
}

# Apply the panel-vs-sumstats allele harmonization using the slim
# harmonizeAlleles against the ldReference's variant info. Threads the
# variant-level filters (indels, strand-ambiguous, duplicates) through
# so the LD-panel-anchored pass handles them in a single sweep.
.matchAgainstLdRef <- function(
    df,
    ldReference,
    matchMinProp,
    removeIndels = FALSE,
    removeStrandAmbiguous = TRUE,
    removeDups = TRUE
) {
    refVariants <- .refVariantsFromLdRef(ldReference)
    flipCandidates <- c("Z", "BETA")
    colToFlip <- intersect(flipCandidates, colnames(df))
    if (length(colToFlip) == 0L) {
        msg <- glue(
            "summaryStatsQc: input entry must contain at least one of ",
            "Z or BETA before panel harmonization."
        )
        abort(msg)
    }
    # An allele swap sends af -> 1 - af, so the DIRECTIONAL AF is what must be
    # complemented. MAF is min(af, 1 - af), which is swap-invariant:
    # complementing it would corrupt it.
    colToComplement <- intersect("AF", colnames(df))
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
    raw <- res$harmonizedData
    out <- if (!is_in("chrom", colnames(raw)) && is_in("chr", colnames(raw))) {
        `colnames<-`(
            raw,
            replace(colnames(raw), colnames(raw) == "chr", "chrom")
        )
    } else {
        raw
    }
    `attr<-`(out, "qcCounts", attr(res, "qcCounts"))
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
    # Single source of truth: the directional AF when present, a directionless
    # MAF/FRQ only as fallback.
    freqCol <- intersect(c("AF", "MAF", "FRQ"), colnames(df))[1L]
    if (is.na(freqCol)) {
        # Skipped, not fatal: a frequency-less study can still run under a
        # default mafCutoff, and one warning says so.
        warn(str_c(
            "summaryStatsQc: mafCutoff > 0 but no af/MAF/FRQ frequency is ",
            "available; skipping the MAF filter."
        ))
        return(list(df = df, dropped = NULL))
    }
    before <- nrow(df)
    # Minor-allele frequency: min(af, 1 - af), direction-agnostic.
    mafVals <- pmin(
        as.numeric(df[[freqCol]]),
        1 - as.numeric(df[[freqCol]]),
        na.rm = FALSE
    )
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
    allN <- as.numeric(df$N)
    before <- nrow(df)
    # A variant with no N cannot be scored against the cohort median.
    hasN <- !is.na(allN)
    nVals <- allN[hasN]
    withN <- if (all(hasN)) df else filter(df, hasN)
    madN <- if (length(nVals) > 0L) stats::mad(nVals, constant = 1) else 0
    # A zero MAD means every retained N is identical, so no variant is an
    # outlier and the z-score would divide by zero.
    kept <- if (madN > 0) {
        filter(withN, abs(nVals - stats::median(nVals)) / madN <= nCutoff)
    } else {
        withN
    }
    list(df = kept, dropped = before - nrow(kept))
}

# An audit record with the unset entries omitted. An all-unset record is the
# empty list, not a zero-length named one, so callers can compare it directly.
# @noRd
.qcAudit <- function(...) {
    entries <- compact(list(...))
    if (length(entries) == 0L) {
        return(list())
    }
    entries
}

# The same, for counters that are reported only when something happened.
# @noRd
.qcAuditPositive <- function(...) {
    exec(.qcAudit, !!!discard(list(...), .qcCountIsZero))
}

# @noRd
.qcCountIsZero <- function(n) is.null(n) || n <= 0L

.applyContentFilters <- function(
    df,
    mafCutoff = 0,
    infoCutoff = 0,
    nCutoff = 5
) {
    maf <- .cfMaf(df, mafCutoff)
    info <- .cfInfo(maf$df, infoCutoff)
    n <- .cfN(info$df, nCutoff)
    list(
        df = n$df,
        audit = .qcAudit(
            mafDropped = maf$dropped,
            infoDropped = info$dropped,
            nDropped = n$dropped
        )
    )
}

# Per-row variant sanity / hygiene checks ported from MungeSumstats's
# check_*.R series but rewritten as pure data.frame operations with no
# genome / dbSNP dependency. Each step is gated by its own flag so a
# caller can disable any single check.
#
# Steps (in order; each contributes a count to audit):
#   - coerceNumeric: cast signed columns to numeric (catches stray "0.5"
#       strings). NA-introducing coercions are counted.
#   - normalizeChr:   spell chromosomes in UCSC style (withChrPrefix():
#       23->chrX, 24->chrY, M/MT->chrM). Optional dropNonstandardChr
#       removes rows whose chromosome is outside chr1..chr22, chrX, chrY,
#       chrM after normalization.
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
# @noRd
.scColumnIsNumeric <- function(col, df) {
    is.numeric(df[[col]])
}

# One column coerced to numeric, with the count of NAs that coercion
# introduced (values that were present but unparseable).
# @noRd
.scCoerceColumn <- function(col, df) {
    orig <- df[[col]]
    coerced <- suppressWarnings(as.numeric(orig))
    list(values = coerced, na = sum(is.na(coerced) & !is.na(orig)))
}

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
    toCoerce <- numericCols[
        !map_lgl(numericCols, .scColumnIsNumeric, df = df)
    ]
    coercions <- map(toCoerce, .scCoerceColumn, df = df)
    naIntroduced <- sum(map_int(coercions, "na"))
    list(
        df = mutate(df, !!!set_names(map(coercions, "values"), toCoerce)),
        audit = .qcAuditPositive(nonNumericCoerced = naIntroduced)
    )
}

# Normalize chromosome labels; optionally drop non-standard chromosomes.
.scNormalizeChr <- function(df, normalizeChr, dropNonstandardChr) {
    if (!normalizeChr || !is_in("chrom", colnames(df))) {
        return(list(df = df, audit = list()))
    }
    chr <- withChrPrefix(df$chrom)
    normalized <- mutate(df, chrom = chr)
    if (!dropNonstandardChr) {
        return(list(df = normalized, audit = list()))
    }
    standardChrs <- str_c("chr", c(as.character(seq_len(22)), "X", "Y", "M"))
    kept <- filter(normalized, is_in(chr, standardChrs))
    list(
        df = kept,
        audit = .qcAuditPositive(
            nonstandardChrDropped = nrow(normalized) - nrow(kept)
        )
    )
}

# Drop rows missing any vital column (chrom/pos/A1/A2 + first signed stat).
.scDropMissData <- function(df, dropMissData) {
    if (!dropMissData || nrow(df) == 0L) {
        return(list(df = df, audit = list()))
    }
    signedCol <- intersect(c("Z", "BETA"), colnames(df))[1L]
    vital <- c(
        intersect(c("chrom", "pos", "A1", "A2"), colnames(df)),
        if (!is.na(signedCol)) signedCol
    )
    if (length(vital) == 0L) {
        return(list(df = df, audit = list()))
    }
    bad <- reduce(map(vital, .scColIsNa, df = df), `|`)
    kept <- filter(df, !bad)
    list(
        df = kept,
        audit = .qcAuditPositive(missDataDropped = nrow(df) - nrow(kept))
    )
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
    list(
        df = df,
        audit = .qcAuditPositive(pOutOfRangeDropped = before - nrow(df))
    )
}

# Clamp tiny P-values up to the floor.
.scClampSmallP <- function(df, clampSmallP, smallPFloor) {
    if (!clampSmallP || !is_in("P", colnames(df)) || nrow(df) == 0L) {
        return(list(df = df, audit = list()))
    }
    p <- as.numeric(df$P)
    smallMask <- !is.na(p) & p >= 0 & p < smallPFloor
    nClamped <- sum(smallMask)
    if (nClamped == 0L) {
        return(list(df = df, audit = list()))
    }
    clamped <- mutate(df, P = replace(.data$P, smallMask, smallPFloor))
    list(df = clamped, audit = list(smallPClamped = nClamped))
}

# Drop rows whose effect equals the null sentinel (0, or 1 for OR).
# Rows whose effect column holds the no-effect sentinel (1 for OR, else 0).
# @noRd
.scZeroEffectMask <- function(col, df) {
    vals <- as.numeric(df[[col]])
    sentinel <- if (col == "OR") 1 else 0
    !is.na(vals) & vals == sentinel
}

.scDropZeroEffect <- function(df, dropZeroEffect) {
    if (!dropZeroEffect || nrow(df) == 0L) {
        return(list(df = df, audit = list()))
    }
    effectCols <- intersect(
        c("BETA", "LOG_ODDS", "SIGNED_SUMSTAT", "OR"),
        colnames(df)
    )
    if (length(effectCols) == 0L) {
        return(list(df = df, audit = list()))
    }
    badMask <- reduce(
        map(effectCols, .scZeroEffectMask, df = df),
        `|`,
        .init = rep(FALSE, nrow(df))
    )
    kept <- filter(df, !badMask)
    list(
        df = kept,
        audit = .qcAuditPositive(zeroEffectDropped = nrow(df) - nrow(kept))
    )
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
    list(
        df = df,
        audit = .qcAuditPositive(nonpositiveSeDropped = before - nrow(df))
    )
}

# Per-row sanity checks: sequence the guarded checks, accumulating the audit.
#' @importFrom purrr list_modify
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
    if (nrow(df) == 0L) {
        return(list(df = df, audit = list()))
    }
    # Order matters: numeric coercion first (later checks compare numbers),
    # then the label fixes, then the row drops, then the small-P clamp.
    steps <- list(
        list(fn = .scCoerceNumeric, args = list(coerceNumeric)),
        list(
            fn = .scNormalizeChr,
            args = list(normalizeChr, dropNonstandardChr)
        ),
        list(fn = .scDropMissData, args = list(dropMissData)),
        list(fn = .scDropPOutOfRange, args = list(dropPOutOfRange)),
        list(fn = .scClampSmallP, args = list(clampSmallP, smallPFloor)),
        list(fn = .scDropZeroEffect, args = list(dropZeroEffect)),
        list(fn = .scDropNonpositiveSe, args = list(dropNonpositiveSe))
    )
    reduce(steps, .scApplyStep, .init = list(df = df, audit = list()))
}

# Run one sanity-check step against the accumulated (df, audit) state. Every
# step takes the frame first and returns list(df, audit); a step that declines
# to run returns the frame unchanged and an empty audit.
# @noRd
.scApplyStep <- function(state, step) {
    r <- exec(step$fn, state$df, !!!step$args)
    list(df = r$df, audit = list_modify(state$audit, !!!r$audit))
}

# Apply ldMismatchQc (SLALOM/DENTIST) against the LD sketch. Returns the
# filtered df, outlier count, and the full per-variant diagnostics table
# (the data.frame returned by ldMismatchQc(), prepended with a
# variant_id column for downstream joins). Callers record `diagnostics`
# in the entry's QC audit so the per-variant detail is available
# for plotting / postprocessing instead of just the outlier count.
# Panel LD for an entry's variants, via the shared LD-from-sketch helper (tuple
# match with chr-prefix tolerance, strand-ambiguous variants kept), with the
# entry narrowed to the variants the panel actually carries.
#
# A variant can survive QC while its LD-panel partner is removed by the panel
# MAF/MAC/missingness filter (common in the study, rare in the panel --
# deletions especially). .panelVariantFilter passes such variants through and
# leaves the drop to onMissing here, so they are dropped here (keeping df
# aligned with the returned LD) rather than aborting the run.
#
# Returns list(R, X, df, dropped): a genotype reference gives X, the factor of
# its LD (.ldFactorFromRef()), which the mismatch QC reads directly; an
# LdMatrix gives the matrix R. Both are NULL when no variant is
# panel-supported, and `df` is empty then -- each caller shapes its own empty
# result. .qcPanelLdMatrix() turns either into the matrix for a step that
# needs one.
# @noRd
.qcPanelSupportedLd <- function(df, ldReference, label) {
    nIn <- nrow(df)
    isFactor <- .ldIsGenotypeRef(ldReference)
    ldFrom <- if (isFactor) .ldFactorFromRef else .ldFromRef
    raw <- ldFrom(ldReference, df$SNP, label = label, onMissing = "drop")
    if (is.null(raw)) {
        return(list(df = df[0L, , drop = FALSE], dropped = nIn))
    }
    keptIds <- attr(raw, "keptVariantIds")
    ld <- `attr<-`(raw, "keptVariantIds", NULL)
    supported <- if (is.null(keptIds) || length(keptIds) >= nIn) {
        list(df = df, dropped = 0L)
    } else {
        list(
            df = filter(df, is_in(.data$SNP, keptIds)),
            dropped = nIn - length(keptIds)
        )
    }
    c(if (isFactor) list(X = ld) else list(R = ld), supported)
}

# The LD matrix of a .qcPanelSupportedLd() result, formed from its factor
# when that is what the reference gave.
# @noRd
.qcPanelLdMatrix <- function(panel) {
    panel$R %||% (crossprod(panel$X) / (nrow(panel$X) - 1))
}

# DENTIST truncates its SVD at a rank bounded by the panel's real size, so it
# needs that size; SLALOM does not, so a reference that cannot say (an
# LdMatrix built without `nSamples`) passes NULL there.
# @noRd
.qcMismatchNSample <- function(ldReference, method) {
    if (identical(method, "dentist")) {
        return(.ldRefRequireNSamples(ldReference, "DENTIST"))
    }
    n <- .ldRefRealNSamples(ldReference)
    if (is.na(n)) NULL else n
}

.applyLdMismatchQcToEntry <- function(df, ldReference, method, panel = NULL) {
    if (is.null(df$SNP) || any(is.na(df$SNP))) {
        abort("summaryStatsQc: ldMismatchQc requires SNP column on the entry.")
    }
    panel <- panel %||%
        .qcPanelSupportedLd(
            df,
            ldReference,
            "summaryStatsQc: ldMismatchQcMethod"
        )
    if (is.null(panel$R) && is.null(panel$X)) {
        return(list(
            df = panel$df,
            outliers = 0L,
            diagnostics = NULL,
            panelUnsupportedDropped = panel$dropped
        ))
    }
    df <- panel$df
    variantIds <- df$SNP
    # A genotype reference's factor goes in as X: SLALOM reads only the lead
    # column of it.
    qc <- ldMismatchQc(
        zScore = df$Z,
        R = panel$R,
        X = panel$X,
        nSample = .qcMismatchNSample(ldReference, method),
        method = method
    )
    # slalom / dentist can leave NA in the outlier column when their
    # per-variant statistic is undefined (e.g. a degenerate dentist
    # chisq for variants effectively orthogonal to the lead). Treat NA as
    # "no evidence of being an outlier" (conservative: keep the variant)
    # so the downstream df / sum() / IRanges construction stay finite.
    outlierFlags <- replace(qc$outlier, is.na(qc$outlier), FALSE)
    list(
        df = filter(df, !outlierFlags),
        outliers = sum(outlierFlags),
        diagnostics = .ldMismatchDiagnostics(qc, variantIds),
        panelUnsupportedDropped = panel$dropped
    )
}

# The QC's per-variant statistics with each variant's id attached, so the
# table stays self-describing once it is separated from the entry; NULL when
# the method returned no table.
# @noRd
.ldMismatchDiagnostics <- function(qc, variantIds) {
    if (!is.data.frame(qc)) {
        return(NULL)
    }
    cbind(
        variant_id = as.character(variantIds),
        qc,
        stringsAsFactors = FALSE
    )
}

# -----------------------------------------------------------------------------
# Signal screen: skip an entry/block with no strong signal by a chosen metric.
# -----------------------------------------------------------------------------
# The screen is driven by ONE metric at a time (enforced by
# SignalScreenParam): pip : max single-effect PIP (susie_ser $pip); cutoff<0
# => 3/nVar absZ : max |Z| (no model fit) logBf : max per-variant single-effect
# logBF (susie_ser $lbf_variable) bf : same evidence on the raw BF scale
# (compare maxlogBF > log(cutoff)) A "screen spec" flowing through the pipelines
# is EITHER a legacy PIP cutoff (numeric scalar, 0 = off -- the historical
# `pipCutoffToSkip`) OR a resolved screen object `list(metric, cutoff)`.
# .asScreen() canonicalizes either into `list(metric, cutoff)` or NULL (no
# screen); the pipeline channels stay untyped so only the screen primitives need
# to interpret the spec. SignalScreenParam() is what the public entry points
# take; .screenResolve() turns one into the spec described above.

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
        return(list(
            df = mutate(
                df,
                N = as.numeric(.data$N_CASE) + as.numeric(.data$N_CONTROL)
            ),
            nSource = "total"
        ))
    }
    if (f$hasScalar) {
        return(list(
            df = mutate(df, N = rep(opts$nCase + opts$nControl, f$nRow)),
            nSource = "total"
        ))
    }
    if (f$hasNSample) {
        return(list(
            df = mutate(df, N = rep(opts$nSample, f$nRow)),
            nSource = "study-n"
        ))
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
        return(list(
            df = mutate(df, N = effectiveN(.data$N_CASE, .data$N_CONTROL)),
            nSource = "effective"
        ))
    }
    if (f$hasScalar) {
        if (f$hasN) {
            .qcEmit(
                lbl,
                "QC track: N overridden by effective N from study ",
                "nCase/nControl."
            )
        }
        return(list(
            df = mutate(
                df,
                N = rep(effectiveN(opts$nCase, opts$nControl), f$nRow)
            ),
            nSource = "effective"
        ))
    }
    if (f$hasN) {
        return(list(df = df, nSource = "column"))
    }
    if (f$hasNSample) {
        return(list(
            df = mutate(df, N = rep(opts$nSample, f$nRow)),
            nSource = "study-n"
        ))
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

# Panel/dosage window indices for the entry, scoped per chromosome. `clamp`
# (see .ssqcRaissBounds) is the stretch of a chromosome the entry may impute
# into when other rows of its study share that chromosome; the flank stops at
# it.
.qcRaissWindowIdx <- function(df, ldReference, flank, clamp = NULL) {
    bounds <- tibble(
        chrom = as.character(df$chrom),
        pos = as.integer(df$pos)
    ) |>
        group_by(.data$chrom) |>
        summarise(
            lo = min(.data$pos, na.rm = TRUE) - flank,
            hi = max(.data$pos, na.rm = TRUE) + flank,
            .groups = "drop"
        ) |>
        .qcRaissClampBounds(clamp)
    # inner_join keeps only sketch SNPs whose chromosome appears in df (the
    # old is_in(skChrom, names(loByChr)) guard); filter keeps those inside the
    # [lo, hi] window. idx carries the original sketch row positions.
    gr <- .ldRefRanges(ldReference)
    tibble(
        chrom = .ldRefChrom(ldReference),
        bp = as.integer(GenomicRanges::start(gr)),
        idx = seq_along(gr)
    ) |>
        inner_join(bounds, by = "chrom") |>
        filter(.data$bp >= .data$lo & .data$bp <= .data$hi) |>
        pull("idx")
}

# Narrow the flanked window to the entry's own stretch of the chromosome.
# @noRd
.qcRaissClampBounds <- function(bounds, clamp) {
    if (is.null(clamp)) {
        return(bounds)
    }
    on <- bounds$chrom == as.character(seqnames(clamp))
    mutate(
        bounds,
        lo = if_else(on, pmax(.data$lo, start(clamp)), .data$lo),
        hi = if_else(on, pmin(.data$hi, end(clamp)), .data$hi)
    )
}

# Assemble refPanel, knownZ table, and scaled dosage for the window.
# The sumstats side of the RAISS inputs: the required columns as a plain
# frame, plus whichever of N / BETA / SE the caller carries, sorted by position
# to match the reference panel.
# @noRd
.qcRaissNumericColumn <- function(nm, df) {
    as.numeric(df[[nm]])
}

.qcRaissKnownZ <- function(df) {
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
    # N / BETA / SE ride along when present: RAISS carries them through to the
    # merged output so the imputed rows sit in the same frame as the known
    # ones.
    optional <- c(N = "n", BETA = "beta", SE = "se")
    present <- optional[is_in(names(optional), colnames(df))]
    extras <- set_names(
        map(names(present), .qcRaissNumericColumn, df = df),
        unname(present)
    )
    knownZ |> mutate(!!!extras) |> arrange(.data$pos)
}

.qcRaissBuildInputs <- function(df, ldReference, windowIdx, opts) {
    # `refIdx` keeps each window variant's position in the reference through
    # the sort, so dosage, frequencies and LD are all read in panel order.
    windowPanel <- .refVariantsFromLdRef(ldReference)[
        windowIdx,
        ,
        drop = FALSE
    ] |>
        mutate(
            variant_id = normalizeVariantId(.data$variant_id),
            refIdx = windowIdx
        ) |>
        arrange(.data$pos)
    knownZ <- .qcRaissKnownZ(df)
    rp <- opts$raissParam
    windowed <- if (.ldIsPanel(ldReference)) {
        # meanImpute = FALSE so per-variant missingness is still visible; the
        # surviving columns are mean-imputed below.
        .ldRefDosage(ldReference, windowPanel$refIdx, meanImpute = FALSE)
    }
    keep <- .qcRaissTargetMask(
        windowPanel,
        knownZ,
        .qcRaissTargetStats(ldReference, windowPanel$refIdx, windowed, rp),
        opts,
        nSamples = .ldRefMacSamples(ldReference, rp$macCutoff)
    )
    c(
        list(
            refPanel = select(windowPanel[keep, , drop = FALSE], -"refIdx"),
            knownZ = knownZ,
            nDroppedTargets = sum(!keep)
        ),
        .qcRaissLdInput(ldReference, windowPanel$refIdx[keep], windowed, keep)
    )
}

# What RAISS imputes from: a panel's standardized dosages or a mixture's
# stacked factor, for its singular-value route, or an LdMatrix's block over
# the kept variants.
# @noRd
.qcRaissLdInput <- function(ldReference, refIdx, windowed, keep) {
    if (.ldIsMatrixRef(ldReference)) {
        return(list(ldMatrix = .ldRefLdAt(ldReference, refIdx)))
    }
    if (.ldIsMixtureRef(ldReference)) {
        return(list(genotypeMatrix = .ldRefFactor(ldReference, refIdx)))
    }
    scaled <- scale(.qtlMeanImpute(windowed[, keep, drop = FALSE]))
    list(genotypeMatrix = replace(scaled, is.na(scaled), 0))
}

# The per-variant statistics the target cutoffs are judged on, or NULL when
# no cutoff is set. A panel of real genotypes measures them on the window's
# dosages, already in hand; a sketch or an LdMatrix answers from what it
# records (see .ldRefVariantStats()).
# @noRd
.qcRaissTargetStats <- function(ldReference, refIdx, windowed, rp) {
    if (rp$mafCutoff <= 0 && rp$macCutoff <= 0 && rp$imissCutoff >= 1) {
        return(NULL)
    }
    if (!is.null(windowed) && !.ldRefIsProjection(ldReference)) {
        return(.qcRaissVariantStats(windowed))
    }
    .ldRefVariantStats(
        ldReference,
        refIdx,
        needMissing = rp$imissCutoff < 1,
        step = "RAISS target filtering"
    )
}

# Which reference-panel variants to keep, given the MAF / MAC / missingness
# cutoffs in `raissParam` and the per-variant `stats` they are judged on
# (NULL when no cutoff is set).
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
.qcRaissTargetMask <- function(refPanel, knownZ, stats, opts, nSamples) {
    keep <- rep(TRUE, nrow(refPanel))
    if (is.null(stats)) {
        return(keep)
    }
    rp <- opts$raissParam
    # Identified exactly as raissSingleMatrixFromX() does, so this mask cannot
    # disagree with the knowns/unknowns split it derives moments later.
    knownIds <- intersect(knownZ$variant_id, refPanel$variant_id)
    isTarget <- !is_in(refPanel$variant_id, knownIds)
    # Candidates only: a known must never be dropped here, or `knownZscores$z`
    # would stop lining up with the LD rows `knowns` indexes.
    unsafe <- isTarget &
        (.raissUnsafeToImpute(refPanel$variant_id, knownZ$variant_id) |
            .raissObservedPosition(refPanel, knownZ))
    effectiveMaf <- .panelEffectiveMaf(rp$mafCutoff, rp$macCutoff, nSamples)
    fails <- (!is.na(stats$maf) & stats$maf < effectiveMaf) | is.na(stats$maf)
    if (!is.null(stats$missRate)) {
        fails <- fails | stats$missRate > rp$imissCutoff
    }
    keep & !(isTarget & fails) & !unsafe
}

# Per-variant MAF and missingness for the RAISS path. A thin alias over the
# shared `.panelVariantStats` (R/ld.R), so the analysis-time panel filter and
# the imputation-target filter cannot drift apart in how they measure a
# variant.
# @noRd
.qcRaissVariantStats <- function(dosage) {
    .panelVariantStats(dosage)
}


# RaissParam() covers two audiences: five fields are raiss()'s own formals and
# four bound what RAISS is asked to impute. Splitting the bundle here lets the
# call site splice rather than name each field, and keeps the two groups from
# being confused for one another.
# @noRd
.raissForwardedNames <- function() {
    c("lamb", "rcond", "svdTol", "r2Threshold", "minimumLd")
}

# @noRd
.raissForwardedArgs <- function(raissParam) {
    as.list(raissParam)[.raissForwardedNames()]
}

# Run RAISS on the inputs .qcRaissBuildInputs() assembled: the window's
# standardized dosages on the genotype route, or its LD block on the matrix
# route. These were built here from the collection's own LD reference, so
# they go straight to the route rather than back through raiss()'s checks.
.qcRaissRun <- function(inp, opts) {
    p <- c(.raissForwardedArgs(opts$raissParam), list(verbose = FALSE))
    if (!is.null(inp$genotypeMatrix)) {
        return(.raissGenotypePath(
            inp$refPanel,
            inp$knownZ,
            inp$genotypeMatrix,
            p
        ))
    }
    .raissLdPath(inp$refPanel, inp$knownZ, inp$ldMatrix, p)
}

# Convert the imputer output back into a QC tibble.
.qcRaissMerge <- function(imputed, knownZ, df) {
    if (is.null(imputed) || is.null(imputed$resultFilter)) {
        return(list(df = df, total = NA_integer_, imputed = 0L))
    }
    impDf <- imputed$resultFilter
    rebuilt <- mutate(
        tibble(
            chrom = impDf$chrom,
            pos = impDf$pos,
            SNP = impDf$variant_id,
            A1 = impDf$A1,
            A2 = impDf$A2,
            Z = impDf$z
        ),
        !!!compact(list(
            N = if (is_in("n", colnames(impDf))) impDf$n,
            BETA = if (is_in("beta", colnames(impDf))) impDf$beta,
            SE = if (is_in("se", colnames(impDf))) impDf$se
        ))
    )
    # RAISS reconstructs the z-score only, so it has no frequency for an
    # imputed variant. But the OBSERVED variants came in with a (harmonized,
    # directional) AF, which the rebuilt frame above would otherwise discard --
    # leaving top_loci$af NA for the whole entry under --impute. Re-attach it by
    # SNP so observed variants keep their AF and imputed variants (absent from
    # `df`) get NA. No-op when the study declared no frequency.
    withAf <- if (is_in("AF", colnames(df))) {
        mutate(rebuilt, AF = as.numeric(df$AF)[match(.data$SNP, df$SNP)])
    } else {
        rebuilt
    }
    out <- .qcFillMissingN(withAf)
    list(df = out, total = nrow(out), imputed = nrow(out) - nrow(knownZ))
}

# Imputed variants carry no sample size, so they inherit the observed median.
# @noRd
.qcFillMissingN <- function(out) {
    if (!is_in("N", colnames(out)) || !any(is.na(out$N))) {
        return(out)
    }
    mutate(out, N = replace_na(.data$N, stats::median(.data$N, na.rm = TRUE)))
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

# Optional RAISS imputation against the ldReference. Returns updated df + audit
# fields + before/after counts (imputation adds variants, so not monotonic).
.qcRaissImpute <- function(df, ldReference, opts, lbl) {
    imputeBefore <- nrow(df)
    flank <- as.integer(opts$raissParam$flank)
    windowIdx <- .qcRaissWindowIdx(df, ldReference, flank, opts$raissBounds)
    if (length(windowIdx) == 0L) {
        .qcEmit(
            lbl,
            "QC track: RAISS imputation skipped ",
            "(no LD-panel variants in the ",
            "region window)."
        )
        return(list(
            df = df,
            audit = list(raissImputedVariants = 0L),
            imputeBefore = imputeBefore,
            imputeAfter = nrow(df)
        ))
    }
    inp <- .qcRaissBuildInputs(df, ldReference, windowIdx, opts)
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
    audit <- c(
        .qcAudit(
            raissTotalVariants = if (!is.na(merged$total)) merged$total
        ),
        list(raissImputedVariants = merged$imputed)
    )
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
    c(
        .qcSeg(sc$nonstandardChrDropped, "nonstdChr"),
        .qcSeg(sc$missDataDropped, "missData"),
        .qcSeg(sc$pOutOfRangeDropped, "badP"),
        .qcSeg(sc$zeroEffectDropped, "zeroEffect"),
        .qcSeg(sc$nonpositiveSeDropped, "badSE"),
        .qcSeg(cf$mafDropped, "maf"),
        .qcSeg(cf$infoDropped, "info"),
        .qcSeg(cf$nDropped, "nCutoff"),
        .qcSeg(qcCount$harmDropped, "harmonization"),
        if (!identical(.engineOf(opts$ldMismatchQcMethod), "none")) {
            str_c("mismatch ", qcCount$mismatchRemoved)
        }
    )
}

# Emit the per-entry QC rollup: corrected (retained), removed, imputed.
.qcEmitRollup <- function(entryAudit, qcCount, opts, nStudyIn, nOut, lbl) {
    removedSegs <- .qcRemovedSegments(entryAudit, qcCount, opts)
    correctedSeg <- str_c(
        "sign-flip ",
        qcCount$harmCorrSign,
        ", strand-flip ",
        qcCount$harmCorrStrand,
        if (isTRUE(opts$alleleFlipKriging)) {
            str_c(", kriging-flip ", qcCount$krigingFlipped)
        } else {
            ""
        }
    )
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
.qcHarmonizeEntry <- function(df, ldReference, opts, lbl) {
    nHarmIn <- nrow(df)
    matched <- .matchAgainstLdRef(
        df,
        ldReference,
        matchMinProp = opts$matchMinProp,
        removeIndels = opts$sumstatsFilterParam$removeIndels,
        removeStrandAmbiguous = opts$sumstatsFilterParam$removeStrandAmbiguous,
        removeDups = TRUE
    )
    harmCounts <- attr(matched, "qcCounts")
    # Re-key SNP to the harmonized id: .matchAgainstLdRef rewrites variant_id
    # to the panel orientation + sign-flips swapped variants but leaves SNP; a
    # stale SNP makes flipped variants miss the panel in later lookups.
    reKeyed <- mutate(
        `attr<-`(matched, "qcCounts", NULL),
        !!!compact(list(SNP = matched$variant_id))
    )
    counts <- list_assign(
        list(
            harmCorrSign = 0L,
            harmCorrStrand = 0L,
            harmDropped = nHarmIn - nrow(reKeyed)
        ),
        !!!compact(list(
            harmCorrSign = harmCounts$signFlip,
            harmCorrStrand = harmCounts$strandFlip
        ))
    )
    .qcHarmonizeReport(
        nrow(reKeyed),
        nHarmIn,
        counts,
        !is.null(harmCounts),
        lbl
    )
    list(
        df = reKeyed,
        audit = list(matchedAgainstLdReference = nrow(reKeyed)),
        counts = counts
    )
}

# Kriging allele-flip QC: sign-flip LD-inconsistent z-scores in place
# (susieR rule logLR > 2 & |z| > 2); variants are corrected and RETAINED.
# Report the variants dropped for having no LD-panel entry; silent when the
# panel supported every one of them.
# @noRd
.qcEmitPanelDrop <- function(lbl, dropped, nIn) {
    if (dropped == 0L) {
        return(invisible(NULL))
    }
    .qcEmit(
        lbl,
        "QC track: dropped ",
        dropped,
        " of ",
        nIn,
        " variant(s) with no LD-panel entry after panel filtering."
    )
}

# The sample size the kriging QC runs at: the caller's nForPip when it is
# usable, else the entry's median N.
# @noRd
.qcKrigingN <- function(df, opts) {
    if (!is.null(opts$nForPip) && is.finite(opts$nForPip)) {
        return(opts$nForPip)
    }
    stats::median(as.numeric(df$N), na.rm = TRUE)
}

# Negate Z (and BETA, where the entry carries one) for the variants kriging
# flagged as sign-flipped against the panel.
# @noRd
.qcApplyKrigingFlips <- function(df, flip) {
    if (!any(flip)) {
        return(df)
    }
    mutate(
        df,
        Z = replace(.data$Z, flip, -df$Z[flip]),
        !!!compact(list(
            BETA = if (is_in("BETA", colnames(df))) {
                replace(df$BETA, flip, -df$BETA[flip])
            }
        ))
    )
}

# @noRd
.qcKrigingFlip <- function(df, ldReference, opts, lbl, panel = NULL) {
    if (!.qcKrigingRuns(df, opts)) {
        return(list(df = df, audit = list(), count = 0L))
    }
    nKrIn <- nrow(df)
    panel <- panel %||%
        .qcPanelSupportedLd(
            df,
            ldReference,
            "summaryStatsQc: kriging prefilter"
        )
    if (is.null(panel$R) && is.null(panel$X)) {
        return(list(
            df = panel$df,
            count = 0L,
            audit = list(krigingFlipped = 0L, panelUnsupportedDropped = nKrIn)
        ))
    }
    supported <- panel$df
    nPanelDrop <- panel$dropped
    .qcEmitPanelDrop(lbl, nPanelDrop, nKrIn)
    # kriging_rss takes nothing but the matrix (dev/ld-factor-audit.md).
    kr <- krigingOutlierQc(
        supported$Z,
        .qcPanelLdMatrix(panel),
        n = .qcKrigingN(supported, opts),
        variantIds = supported$SNP
    )
    nKr <- sum(kr$flip)
    flipped <- .qcApplyKrigingFlips(supported, kr$flip)
    .qcEmit(
        lbl,
        "QC track: kriging sign-flipped ",
        nKr,
        " of ",
        nKrIn,
        " LD-inconsistent variant(s)."
    )
    list(
        df = flipped,
        count = nKr,
        audit = list(
            krigingFlipped = nKr,
            krigingDiagnostics = kr$diagnostics,
            panelUnsupportedDropped = nPanelDrop
        )
    )
}

# TRUE when the kriging allele-flip check runs on this entry.
# @noRd
.qcKrigingRuns <- function(df, opts) {
    isTRUE(opts$alleleFlipKriging) && nrow(df) >= 2L
}

# TRUE when an LD-mismatch check runs on this entry.
# @noRd
.qcMismatchRuns <- function(df, opts) {
    !identical(.engineOf(opts$ldMismatchQcMethod), "none") && nrow(df) >= 2L
}

# LD-mismatch QC (SLALOM / DENTIST). Retains full per-variant diagnostics.
# `panel` is the entry's panel-supported LD when the caller has already read
# it (see .qcPhaseLdChecks()).
.qcMismatchQc <- function(df, ldReference, opts, lbl, panel = NULL) {
    if (!.qcMismatchRuns(df, opts)) {
        return(list(df = df, audit = list(), count = 0L))
    }
    nMmIn <- nrow(df)
    ldQc <- .applyLdMismatchQcToEntry(
        df,
        ldReference,
        opts$ldMismatchQcMethod,
        panel = panel
    )
    df <- ldQc$df
    baseAudit <- list(
        ldMismatchOutliersDropped = ldQc$outliers,
        ldMismatchMethod = .engineOf(opts$ldMismatchQcMethod)
    )
    nPanelDrop <- ldQc$panelUnsupportedDropped %||% 0L
    audit <- c(
        baseAudit,
        .qcAudit(
            panelUnsupportedDropped = if (nPanelDrop > 0L) nPanelDrop,
            ldMismatchDiagnostics = ldQc$diagnostics
        )
    )
    if (nPanelDrop > 0L) {
        .qcEmit(
            lbl,
            "QC track: dropped ",
            nPanelDrop,
            " of ",
            nMmIn,
            " variant(s) with no LD-panel entry after panel filtering."
        )
    }
    .qcEmit(
        lbl,
        "QC track: ",
        .engineOf(opts$ldMismatchQcMethod),
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
    # The boundary: the bundle travels this far, and .applySanityChecks and
    # the per-check helpers below it keep plain scalars.
    cleaning <- opts$sumstatsCleaningParam
    sanity <- exec(.applySanityChecks, df, !!!cleaning)
    df <- sanity$df
    audit <- .qcAudit(
        sanityChecks = if (length(sanity$audit) > 0L) sanity$audit
    )
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
    applied <- isTRUE(is_in(nRes$nSource, c("effective", "total"))) &&
        is_in("N", colnames(df))
    list(
        df = df,
        nSource = nRes$nSource,
        opts = list_assign(
            opts,
            !!!compact(list(
                nForPip = if (applied) {
                    stats::median(as.numeric(df$N), na.rm = TRUE)
                }
            ))
        )
    )
}

# Variant-content filters (MAF / INFO / N).
.qcStepContentFilters <- function(df, opts, lbl) {
    nFiltIn <- nrow(df)
    # One mafCutoff drives two filters: the sumstats' own MAF column here,
    # and the LD panel in .ssqcPrunePanel. It has always been a single knob,
    # so it is read from panelFilter in both places rather than quietly
    # becoming two separate settings.
    cf <- .applyContentFilters(
        df,
        mafCutoff = opts$panelFilterParam$mafCutoff,
        infoCutoff = opts$sumstatsFilterParam$infoCutoff,
        nCutoff = opts$sumstatsFilterParam$nCutoff
    )
    df <- cf$df
    audit <- .qcAudit(
        contentFilters = if (length(cf$audit) > 0L) cf$audit
    )
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
    withBetaSe <- list_assign(
        entryAudit,
        !!!compact(list(betaSeFromZ = derived$audit))
    )
    if (!is_in("Z", colnames(derived$df)) || is_in("P", colnames(derived$df))) {
        return(list(df = derived$df, entryAudit = withBetaSe))
    }
    withP <- mutate(derived$df, P = .zToPvalue(.data$Z))
    clamped <- .qcClampDerivedP(withP, opts)
    list(
        df = clamped$df,
        entryAudit = list_assign(
            withBetaSe,
            pValueFromZ = sum(!is.na(withP$P)),
            !!!compact(list(
                sanityChecks = .qcAddClampCount(
                    entryAudit$sanityChecks,
                    clamped$nClamped
                )
            ))
        )
    )
}

# Re-clamp the P values just derived from Z, and report how many moved.
# @noRd
.qcClampDerivedP <- function(df, opts) {
    cleaning <- opts$sumstatsCleaningParam
    if (!isTRUE(cleaning$clampSmallP) || nrow(df) == 0L) {
        return(list(df = df, nClamped = 0L))
    }
    floorP <- cleaning$smallPFloor
    smallMask <- !is.na(df$P) & df$P >= 0 & df$P < floorP
    if (sum(smallMask) == 0L) {
        return(list(df = df, nClamped = 0L))
    }
    list(
        df = mutate(df, P = replace(.data$P, smallMask, floorP)),
        nClamped = sum(smallMask)
    )
}

# Fold this step's clamp count into the sanity-check record the earlier
# per-row checks may already have written. NULL when nothing was clamped, so
# the caller can drop the key rather than record an empty record.
# @noRd
.qcAddClampCount <- function(sanityChecks, nClamped) {
    if (nClamped == 0L) {
        return(NULL)
    }
    list_assign(
        sanityChecks %||% list(),
        smallPClamped = (sanityChecks$smallPClamped %||% 0L) + nClamped
    )
}

# keepVariants subset + skipRegion drop.
.qcStepKeepSkip <- function(df, opts) {
    nIn <- nrow(df)
    kept <- if (length(opts$keepVariants) > 0L) {
        filter(df, is_in(.data$SNP, opts$keepVariants))
    } else {
        df
    }
    nKept <- nrow(kept)
    skipped <- if (!is.null(opts$skipRegion) && length(opts$skipRegion) > 0L) {
        .applySkipRegion(kept, opts$skipRegion)
    } else {
        kept
    }
    list(
        df = skipped,
        audit = .qcAudit(
            keepVariantsDropped = if (length(opts$keepVariants) > 0L) {
                nIn - nKept
            },
            skipRegionDropped = if (
                !is.null(opts$skipRegion) && length(opts$skipRegion) > 0L
            ) {
                nKept - nrow(skipped)
            }
        )
    )
}

# Optional post-harmonization signal screen (PIP / |Z| / BF / logBF).
.qcStepScreen <- function(df, opts) {
    if (is.null(opts$screen)) {
        return(list(df = df, audit = list()))
    }
    scr <- .applyEntryScreen(df, n = opts$nForPip, screen = opts$screen)
    list(
        df = scr$df,
        audit = c(
            list(pipScreenSkipped = isTRUE(scr$skipped)),
            .qcAudit(pipScreenReason = if (isTRUE(scr$skipped)) scr$reason)
        )
    )
}

# Optional RAISS imputation step wrapper (guarded; threads impute counters).
.qcRaissImputeStep <- function(df, ldReference, opts, lbl, qcCount) {
    if (!isTRUE(opts$impute) || nrow(df) < 1L) {
        return(list(df = df, audit = list(), qcCount = qcCount))
    }
    imp <- .qcRaissImpute(df, ldReference, opts, lbl)
    list(
        df = imp$df,
        audit = imp$audit,
        qcCount = list_assign(
            qcCount,
            imputeBefore = imp$imputeBefore,
            imputeAfter = imp$imputeAfter
        )
    )
}

# Early return when too few variants survive pre-harmonization QC.
.qcEarlyExit <- function(df, entryAudit, qcCount, opts, nIn, lbl) {
    withExit <- list_assign(
        entryAudit,
        earlyExit = "fewer than two variants after pre-harmonization QC"
    )
    # Still emit the rollup. An entry that QC empties is precisely the case a
    # user needs told about, and returning early used to make those drops
    # invisible in the log -- the audit recorded them, nothing said so.
    .qcEmitRollup(withExit, qcCount, opts, nIn, nrow(df), lbl)
    list(gr = .dfToEntryGranges(df), audit = withExit)
}

# Each QC phase below is (state) -> state, so the sequence reads as a pipeline
# rather than a run of reassignments. `state` carries df + entryAudit + opts +
# qcCount + lbl (+ ldReference for the panel-aware phases); a phase touches only
# the parts it owns and passes the rest through.

# Per-row sanity checks.
# @noRd
.qcPhaseSanity <- function(state) {
    san <- .qcStepSanity(state$df, state$opts, state$lbl)
    list_assign(
        state,
        df = san$df,
        entryAudit = list_modify(state$entryAudit, !!!san$audit)
    )
}

# Effective-N canonicalization, which also rewrites `opts$nForPip`.
# @noRd
.qcPhaseEffectiveN <- function(state) {
    eff <- .qcStepEffectiveN(state$df, state$opts, state$lbl)
    list_assign(
        state,
        df = eff$df,
        entryAudit = list_assign(state$entryAudit, nSource = eff$nSource),
        opts = eff$opts
    )
}

# MAF / INFO / N content filters.
# @noRd
.qcPhaseContentFilters <- function(state) {
    cf <- .qcStepContentFilters(state$df, state$opts, state$lbl)
    list_assign(
        state,
        df = cf$df,
        entryAudit = list_modify(state$entryAudit, !!!cf$audit)
    )
}

# BETA/SE from Z, then P from Z. Owns the audit outright (it folds its own
# clamp count into the sanity record the earlier phase may have written).
# @noRd
.qcPhaseDerive <- function(state) {
    der <- .qcStepDerive(state$df, state$opts, state$entryAudit)
    list_assign(state, df = der$df, entryAudit = der$entryAudit)
}

# keepVariants subset + skipRegion drop.
# @noRd
.qcPhaseKeepSkip <- function(state) {
    ks <- .qcStepKeepSkip(state$df, state$opts)
    list_assign(
        state,
        df = ks$df,
        entryAudit = list_modify(state$entryAudit, !!!ks$audit)
    )
}

# Panel-vs-sumstats allele harmonization.
# @noRd
.qcPhaseHarmonize <- function(state) {
    harm <- .qcHarmonizeEntry(
        state$df,
        state$ldReference,
        state$opts,
        state$lbl
    )
    list_assign(
        state,
        df = harm$df,
        entryAudit = list_modify(state$entryAudit, !!!harm$audit),
        qcCount = list_modify(state$qcCount, !!!harm$counts)
    )
}

# Post-harmonization signal screen.
# @noRd
.qcPhaseScreen <- function(state) {
    scr <- .qcStepScreen(state$df, state$opts)
    list_assign(
        state,
        df = scr$df,
        entryAudit = list_modify(state$entryAudit, !!!scr$audit)
    )
}

# Kriging sign-flips and LD-mismatch QC. One phase, because the mismatch check
# reads the kriging-corrected frame and the two counts are reported together.
# @noRd
.qcPhaseLdChecks <- function(state) {
    # Both checks read the same panel-supported variants -- kriging only flips
    # signs -- so when both run the reference is read once, for the two.
    shared <- if (
        .qcKrigingRuns(state$df, state$opts) &&
            .qcMismatchRuns(state$df, state$opts)
    ) {
        .qcPanelSupportedLd(state$df, state$ldReference, "summaryStatsQc")
    }
    kr <- .qcKrigingFlip(
        state$df,
        state$ldReference,
        state$opts,
        state$lbl,
        panel = shared
    )
    mm <- .qcMismatchQc(
        kr$df,
        state$ldReference,
        state$opts,
        state$lbl,
        # Kriging has already recorded the panel's drops; its output carries
        # only supported variants.
        panel = if (!is.null(shared)) {
            list_assign(shared, df = kr$df, dropped = 0L)
        }
    )
    list_assign(
        state,
        df = mm$df,
        entryAudit = list_modify(
            list_modify(state$entryAudit, !!!kr$audit),
            !!!mm$audit
        ),
        qcCount = list_assign(
            state$qcCount,
            krigingFlipped = kr$count,
            mismatchRemoved = mm$count
        )
    )
}

# Optional RAISS imputation.
# @noRd
.qcPhaseImpute <- function(state) {
    imp <- .qcRaissImputeStep(
        state$df,
        state$ldReference,
        state$opts,
        state$lbl,
        state$qcCount
    )
    list_assign(
        state,
        df = imp$df,
        entryAudit = list_modify(state$entryAudit, !!!imp$audit),
        qcCount = imp$qcCount
    )
}

# Pre-harmonization phase: init + sanity + effective-N + content + derive +
# keep/skip. Returns the QC state carried into the harmonization phase.
.qcPreHarmonize <- function(gr, opts, entryLabel) {
    init <- .qcInitEntry(gr, entryLabel)
    final <- list(
        df = init$df,
        entryAudit = init$entryAudit,
        opts = opts,
        lbl = init$lbl
    ) |>
        .qcPhaseSanity() |>
        .qcPhaseEffectiveN() |>
        .qcPhaseContentFilters() |>
        .qcPhaseDerive() |>
        .qcPhaseKeepSkip()
    list(
        df = final$df,
        entryAudit = final$entryAudit,
        opts = final$opts,
        qcCount = init$qcCount,
        lbl = final$lbl,
        nStudyIn = init$nStudyIn
    )
}

.runEntrySummaryStatsQc <- function(
    gr,
    ldReference,
    opts,
    entryLabel = NULL
) {
    pre <- .qcPreHarmonize(gr, opts, entryLabel)
    if (nrow(pre$df) < 2L) {
        return(.qcEarlyExit(
            pre$df,
            pre$entryAudit,
            pre$qcCount,
            pre$opts,
            length(gr),
            pre$lbl
        ))
    }
    final <- list(
        df = pre$df,
        entryAudit = pre$entryAudit,
        qcCount = pre$qcCount,
        opts = pre$opts,
        lbl = pre$lbl,
        ldReference = ldReference
    ) |>
        .qcPhaseHarmonize() |>
        .qcPhaseScreen() |>
        .qcPhaseLdChecks() |>
        .qcPhaseImpute()
    .qcEmitRollup(
        final$entryAudit,
        final$qcCount,
        final$opts,
        pre$nStudyIn,
        nrow(final$df),
        final$lbl
    )
    list(
        gr = .dfToEntryGranges(final$df),
        audit = list_assign(final$entryAudit, variantsOut = nrow(final$df))
    )
}

# Number of ranges an entry contributes, NULL-safe.
# @noRd
.entryRangeCount <- function(gr) {
    if (is.null(gr)) 0L else length(gr)
}

# Shrink an LD-sketch GenotypeHandle to the panel variants inside the summary
# statistics' per-chromosome position span. `entries` is a list/SimpleList of
# per-study (or per-tuple) GRanges. A genome-wide sketch otherwise carries a
# full-genome snpInfo; only variants inside [min,max] BP of each represented
# chromosome are reachable by harmonization or within-range imputation, so the
# rest is dropped at load time. NULL-safe; a no-op when the span already covers
# the panel. See [[.subsetGenotypeHandle]] for why this is read-safe.
# @noRd
.subsetLdRefToRange <- function(ldReference, entries) {
    if (is.null(ldReference)) {
        return(NULL)
    }
    chromAll <- unname(list_c(map(entries, .entryChrom)))
    posAll <- unname(list_c(map(entries, .entryPos)))
    ok <- !is.na(chromAll) & !is.na(posAll)
    chrom <- chromAll[ok]
    pos <- posAll[ok]
    if (length(pos) == 0L) {
        # No span to keep: a zero-variant object references no LD, so the
        # retained panel is empty rather than the full genome-wide sketch.
        return(.emptyLdRef(ldReference))
    }
    bounds <- tibble(chrom = chrom, pos = pos) |>
        group_by(chrom) |>
        summarise(lo = min(pos), hi = max(pos), .groups = "drop")
    gr <- .ldRefRanges(ldReference)
    # left_join keeps every sketch SNP; one on a chromosome absent from the
    # entries gets lo/hi = NA -> inWindow FALSE (the old is_in guard). The
    # explicit is.na() guards force FALSE (never NA) for an absent chrom or an
    # NA bp, replacing the base keep[is.na(keep)] <- FALSE.
    keep <- tibble(
        chrom = .ldRefChrom(ldReference),
        bp = as.integer(GenomicRanges::start(gr))
    ) |>
        left_join(bounds, by = "chrom") |>
        mutate(
            inWindow = !is.na(.data$lo) &
                !is.na(.data$bp) &
                .data$bp >= .data$lo &
                .data$bp <= .data$hi
        ) |>
        pull("inWindow")
    .ldRefSubset(ldReference, keep)
}

# Shrink an LD-sketch GenotypeHandle to EXACTLY the variants present across the
# QC'd `entries` (imputation may have added variants; QC may have dropped some),
# matched by canonical variant id. Applied at the end of summaryStatsQc so the
# retained sketch mirrors the object's final variant set. NULL-safe.
# @noRd
.subsetLdRefToIds <- function(ldReference, entries) {
    if (is.null(ldReference)) {
        return(NULL)
    }
    # Genuine variant count (ranges) across entries, read independently of the
    # SNP mcol: zero ranges means the object references no LD, so the panel is
    # emptied. That is distinct from the pathological "ranges present but SNP
    # mcol absent" case below, which stays conservative and keeps the full
    # sketch rather than blanking a panel it simply could not key.
    nRanges <- sum(map_int(entries, .entryRangeCount))
    if (nRanges == 0L) {
        return(.emptyLdRef(ldReference))
    }
    ids <- unname(list_c(map(entries, .entrySnpIds)))
    if (length(ids) == 0L) {
        return(ldReference)
    }
    # Match on the (chrom, pos, allele) tuple rather than the raw id string.
    # The entry ids have just been canonicalized to chr:pos:A2:A1 by the QC
    # above, while panel ids are passed through verbatim by
    # `.repairVariantIds()` whenever their allele fields are already valid DNA
    # -- so a panel keyed chr:pos:A1:A2 (which is how a PLINK .bim writes ids;
    # a .pvar writes REF:ALT and is therefore already canonical) matches
    # nothing under string equality and the panel is silently emptied. This is
    # the matcher `.ldFromRefMatch()` and every other LD lookup already use.
    panelIds <- .ldRefMatchIds(ldReference)
    matched <- matchVariants(
        panelIds,
        unique(ids),
        removeStrandAmbiguous = FALSE
    )
    keep <- is_in(seq_along(panelIds), matched$idxA)
    .ldRefSubset(ldReference, keep)
}

# --- summaryStatsQc orchestration helpers ----------------------------------

# Validate input class + per-entry MAF/INFO column availability.
#' @importFrom checkmate assertMultiClass
.ssqcCheckEntries <- function(sumstats, infoCutoff) {
    assertMultiClass(sumstats, c("QtlSumStats", "GwasSumStats"))
    for (i in seq_len(nrow(sumstats))) {
        cols <- colnames(S4Vectors::mcols(.collectionEntry(sumstats, i)))
        # mafCutoff no longer pre-aborts on a missing frequency: .cfMaf
        # derives the MAF from the directional AF (or a fallback MAF/FRQ) and
        # skips with one warning when none is available.
        if (infoCutoff > 0 && !is_in("INFO", cols)) {
            msg <- glue(
                "summaryStatsQc: infoCutoff > 0 requires every entry to ",
                "carry an INFO column; entry {i} does not."
            )
            abort(msg)
        }
    }
}

# TRUE for a single finite non-negative cutoff value.
# @noRd
.ssqcIsCutoff <- function(x) {
    is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) && x >= 0
}

# Shape-check the panel-filter trio. Runs before the panel is opened, so a
# malformed cutoff fails on the call rather than after a dosage read.
# @noRd
.ssqcCheckPanelCutoffs <- function(panelFilterParam) {
    bad <- names(panelFilterParam)[
        !map_lgl(as.list(panelFilterParam), .ssqcIsCutoff)
    ]
    if (length(bad) > 0L) {
        msg <- glue(
            "summaryStatsQc: {str_flatten(bad, ', ')} must each be a single ",
            "finite number >= 0."
        )
        abort(msg)
    }
    invisible(NULL)
}

# Build the per-entry QC options list from the captured call parameters.
.ssqcBuildOpts <- function(
    sumstatsFilterParam,
    panelFilterParam,
    skipRegion,
    ldMismatchQcMethod,
    alleleFlipKriging,
    effectiveN,
    impute,
    raissParam,
    matchMinProp,
    sumstatsCleaningParam,
    keepVariants,
    signalScreenParam
) {
    # Built explicitly rather than by subsetting a captured environment with a
    # character vector: a renamed argument is now an error, not a NULL entry.
    opts <- list(
        sumstatsFilterParam = sumstatsFilterParam,
        panelFilterParam = panelFilterParam,
        skipRegion = skipRegion,
        ldMismatchQcMethod = ldMismatchQcMethod,
        alleleFlipKriging = alleleFlipKriging,
        effectiveN = effectiveN,
        impute = impute,
        raissParam = raissParam,
        matchMinProp = matchMinProp,
        sumstatsCleaningParam = sumstatsCleaningParam
    )
    # nCase / nControl / nForPip are deliberately absent here: they are
    # per-entry, and `.ssqcEntryOpts()` adds them for the entry being run.
    list_assign(
        opts,
        keepVariants = as.character(keepVariants),
        screen = .screenResolve(signalScreenParam)
    )
}

# Per-entry sample-size options: median N for PIP, study case/control/total N.
.ssqcEntryOpts <- function(opts, sumstats, i) {
    mc <- S4Vectors::mcols(.collectionEntry(sumstats, i))
    cols <- .tupleColumnNames(sumstats)
    list_assign(
        opts,
        !!!compact(list(
            nForPip = if (is_in("N", colnames(mc))) {
                stats::median(mc$N, na.rm = TRUE)
            },
            nCase = if (is_in("nCase", cols)) {
                as.numeric(sumstats$nCase)[[i]]
            },
            nControl = if (is_in("nControl", cols)) {
                as.numeric(sumstats$nControl)[[i]]
            },
            nSample = if (is_in("nSample", cols)) {
                as.numeric(sumstats$nSample)[[i]]
            }
        ))
    )
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

# Drop the reference-panel variants that fail the MAF / MAC / missingness
# cutoffs, returning the pruned sketch. NULL cutoffs (the default) are a no-op.
#
# Pruning the PANEL rather than masking a variant list is what makes one
# cutoff reach every reader: harmonization drops the study variants the
# retained panel no longer covers, and kriging, SLALOM/DENTIST and the RAISS
# known/unknown split then only ever see what survived.
# @noRd
.ssqcPrunePanel <- function(ldReference, cutoffs, label) {
    if (is.null(ldReference) || is.null(cutoffs)) {
        return(ldReference)
    }
    ids <- .ldRefMatchIds(ldReference)
    if (length(ids) == 0L) {
        return(ldReference)
    }
    keep <- .panelKeepMask(ids, ldReference, cutoffs, label)
    if (all(keep)) {
        return(ldReference)
    }
    if (.ldIsMatrixRef(ldReference) || .ldIsMixtureRef(ldReference)) {
        return(.ldRefSubset(ldReference, keep))
    }
    if (methods::is(ldReference, "GenotypeHandle")) {
        return(.subsetGenotypeHandle(ldReference, keep))
    }
    # A panel is cut where it can shed the dropped variants from its seed as
    # well as from its rows (.panelSubsetRows()). The mask is over the
    # PANEL's rows, so it is never applied to a handle whose rows differ.
    .panelSubsetRows(ldReference, keep)
}

# One entry's QC result: list(gr, audit). `opts` arrives as the shared base
# options; the per-entry fields are derived here rather than carried over from
# the previous entry.
# @noRd
.ssqcRunEntry <- function(i, sumstats, opts, refs, keys, isQtl, bounds) {
    .runEntrySummaryStatsQc(
        gr = .collectionEntry(sumstats, i),
        ldReference = .ldRefForKey(keys[[i]], refs),
        opts = list_assign(
            .ssqcEntryOpts(opts, sumstats, i),
            raissBounds = bounds[[i]]
        ),
        entryLabel = .ssqcEntryLabel(sumstats, i, isQtl)
    )
}

# The stretch of its chromosome each row may impute into, or NULL for a row
# with no window of its own that shares its chromosome with no other row of
# its tuple. A study split into LD blocks is several rows on one chromosome,
# and RAISS widening one row's window by `flank` would impute its neighbours'
# variants -- a second, imputed copy of a variant the next row observes. So
# each row keeps to its own window: the `region` it was given, which bounds it
# whether or not its neighbours were loaded, else the midpoints between its
# span and its neighbours'. Without a flank a row's window is its own
# variants' span, which needs no bound.
# @noRd
.ssqcRaissBounds <- function(sumstats, opts) {
    n <- nrow(sumstats)
    if (!isTRUE(opts$impute) || as.integer(opts$raissParam$flank) <= 0L) {
        return(vector("list", n))
    }
    chrom <- map_chr(as.list(sumstats), .rtlElementSeqname)
    map(
        seq_len(n),
        .ssqcRowBounds,
        group = str_c(.ssqcTupleKey(sumstats), chrom, sep = "\r"),
        spans = range(sumstats),
        blocks = .ssqcGivenBlocks(sumstats)
    )
}

# Each row's given window, or NULL where it has none: a computed window is
# only the span the midpoints are taken from.
# @noRd
.ssqcGivenBlocks <- function(sumstats) {
    regions <- .rowRegions(sumstats)
    map(seq_len(nrow(sumstats)), .ssqcGivenBlockAt, regions = regions)
}

# @noRd
.ssqcGivenBlockAt <- function(i, regions) {
    if (regions$regionGiven[[i]]) regions$region[i] else NULL
}

# The rows of one tuple: a GWAS study, or a QTL (study, context, trait).
# @noRd
.ssqcTupleKey <- function(sumstats) {
    cols <- intersect(
        c("study", "context", "trait"),
        .tupleColumnNames(sumstats)
    )
    exec(str_c, !!!map(cols, .ssqcKeyPart, sumstats = sumstats), sep = "\r")
}

# @noRd
.ssqcKeyPart <- function(col, sumstats) {
    as.character(.tupleColumn(sumstats, col))
}

# @noRd
.ssqcRowBounds <- function(i, group, spans, blocks) {
    if (!is.null(blocks[[i]])) {
        return(blocks[[i]])
    }
    if (is.na(group[[i]])) {
        return(NULL)
    }
    sibs <- setdiff(which(group == group[[i]]), i)
    if (length(sibs) == 0L) {
        return(NULL)
    }
    .ssqcMidpointBounds(spans[[i]], unlist(spans[sibs], use.names = FALSE))
}

# Halfway to the nearest neighbour on each side; a neighbour overlapping the
# row leaves it only its own span.
# @noRd
.ssqcMidpointBounds <- function(own, sib) {
    s <- start(own)
    e <- end(own)
    if (any(start(sib) <= e & end(sib) >= s)) {
        return(own)
    }
    below <- end(sib)[end(sib) < s]
    above <- start(sib)[start(sib) > e]
    GenomicRanges::GRanges(
        seqnames(own),
        IRanges::IRanges(
            if (length(below) == 0L) 1L else (max(below) + s) %/% 2L + 1L,
            if (length(above) == 0L) {
                .Machine$integer.max
            } else {
                (e + min(above)) %/% 2L
            }
        )
    )
}

# Run the per-entry QC pipeline across all entries, each against its own
# row's LD reference.
.ssqcRunEntries <- function(sumstats, opts) {
    isQtl <- methods::is(sumstats, "QtlSumStats")
    # Panel filter, once per distinct LD reference and BEFORE any entry is
    # harmonized, so a variant the panel cannot support is gone from the LD
    # and from the summary statistics alike.
    pruned <- map(
        ldReference(sumstats),
        .ssqcPrunePanel,
        cutoffs = .panelCutoffs(opts$panelFilterParam),
        label = "summaryStatsQc"
    )
    keys <- .ldRefKeys(sumstats)
    results <- map(
        seq_len(nrow(sumstats)),
        .ssqcRunEntry,
        sumstats = sumstats,
        opts = opts,
        refs = pruned,
        keys = keys,
        isQtl = isQtl,
        bounds = .ssqcRaissBounds(sumstats, opts)
    )
    list(
        newEntries = map(results, "gr"),
        entryAudits = map(results, "audit"),
        ldReferences = pruned,
        keys = keys
    )
}

# One screen metric's cutoff as qcOptions echoes it. The bundle leaves an unset
# metric absent; the echoed record has always spelled "off" as 0, and readers
# like qcDiagnostics() index these four names, so the output shape is kept
# even though the input is now one argument.
# @noRd
.ssqcScreenEcho <- function(signalScreenParam, metric) {
    signalScreenParam[[metric]] %||% 0
}

# Assemble the qcOptions record: the options QC ran with, echoed.
.ssqcBuildQcOptions <- function(
    sumstatsFilterParam,
    panelFilterParam,
    signalScreenParam,
    ldMismatchQcMethod,
    alleleFlipKriging,
    effectiveN,
    impute,
    sumstatsCleaningParam
) {
    list(
        timestamp = NA_character_,
        # Echoed FLAT, not nested: this record is part of the returned
        # object's shape, which readers like qcDiagnostics() index by
        # these names. The bundles are an argument convention, not an
        # output one.
        options = list(
            removeIndels = sumstatsFilterParam$removeIndels,
            removeStrandAmbiguous = sumstatsFilterParam$removeStrandAmbiguous,
            mafCutoff = panelFilterParam$mafCutoff,
            macCutoff = panelFilterParam$macCutoff,
            imissCutoff = panelFilterParam$imissCutoff,
            infoCutoff = sumstatsFilterParam$infoCutoff,
            nCutoff = sumstatsFilterParam$nCutoff,
            pipCutoffToSkip = .ssqcScreenEcho(signalScreenParam, "pip"),
            absZCutoffToSkip = .ssqcScreenEcho(signalScreenParam, "absZ"),
            bfCutoffToSkip = .ssqcScreenEcho(signalScreenParam, "bf"),
            logBfCutoffToSkip = .ssqcScreenEcho(signalScreenParam, "logBf"),
            ldMismatchQcMethod = ldMismatchQcMethod,
            alleleFlipKriging = alleleFlipKriging,
            effectiveN = effectiveN,
            impute = impute,
            # Echoed FLAT for the same reason the screen cutoffs are: these
            # nine names are part of the returned object's shape.
            coerceNumeric = sumstatsCleaningParam$coerceNumeric,
            normalizeChr = sumstatsCleaningParam$normalizeChr,
            dropNonstandardChr = sumstatsCleaningParam$dropNonstandardChr,
            dropMissData = sumstatsCleaningParam$dropMissData,
            dropPOutOfRange = sumstatsCleaningParam$dropPOutOfRange,
            clampSmallP = sumstatsCleaningParam$clampSmallP,
            smallPFloor = sumstatsCleaningParam$smallPFloor,
            dropZeroEffect = sumstatsCleaningParam$dropZeroEffect,
            dropNonpositiveSe = sumstatsCleaningParam$dropNonpositiveSe
        )
    )
}

# Rebuild the SumStats object with QC'd entries, trimmed sketches (one per
# row, as a list running parallel to `newEntries`), the qcOptions record and
# each row's audit (the qcAudit column). Every column the constructor does
# not rebuild itself is carried across as it was. So is each row's analysis
# window: QC changes the variants, not where the window is -- the one given,
# or the span the row was first built with -- and the span of what is left
# is range(x).
.ssqcRebuild <- function(
    sumstats,
    newEntries,
    newLdReference,
    qcOptions,
    audits
) {
    rebuilt <- if (methods::is(sumstats, "GwasSumStats")) {
        .ssqcRebuildGwas(sumstats, newEntries, newLdReference, qcOptions)
    } else {
        .ssqcRebuildQtl(sumstats, newEntries, newLdReference, qcOptions)
    }
    audited <- `$<-`(rebuilt, "qcAudit", S4Vectors::SimpleList(audits))
    .withRegionColumns(audited, .rowRegions(sumstats))
}

# @noRd
.ssqcRebuildQtl <- function(sumstats, newEntries, newLdReference, qcOptions) {
    QtlSumStats(
        studyName = as.character(sumstats$study),
        context = as.character(sumstats$context),
        trait = as.character(sumstats$trait),
        entry = newEntries,
        genome = unique(unname(GenomeInfoDb::genome(sumstats))),
        ldReference = newLdReference,
        varY = as.numeric(sumstats$varY),
        nSample = .ssqcNumericColumn(sumstats, "nSample"),
        qcOptions = qcOptions,
        traitPos = .ssqcColumnOrNull(sumstats, "traitPos"),
        extraCols = .ssqcExtraCols(sumstats, .SSQC_QTL_COLUMNS)
    )
}

# @noRd
.ssqcRebuildGwas <- function(sumstats, newEntries, newLdReference, qcOptions) {
    GwasSumStats(
        studyName = as.character(sumstats$study),
        entry = newEntries,
        genome = unique(unname(GenomeInfoDb::genome(sumstats))),
        ldReference = newLdReference,
        varY = as.numeric(sumstats$varY),
        nCase = .ssqcNumericColumn(sumstats, "nCase"),
        nControl = .ssqcNumericColumn(sumstats, "nControl"),
        nSample = .ssqcNumericColumn(sumstats, "nSample"),
        qcOptions = qcOptions,
        # Carry the existing block keys through: the entries are already
        # split, so without this the constructor re-derives them from
        # seqname and every block on one chromosome collapses to one key.
        blockId = .ssqcColumnOrNull(sumstats, "blockId"),
        extraCols = .ssqcExtraCols(sumstats, .SSQC_GWAS_COLUMNS)
    )
}

# The columns each constructor builds from its own arguments. Anything else a
# collection carries -- an `extraCols` column -- goes back in through
# `extraCols`.
# @noRd
.SSQC_SHARED_COLUMNS <- c(
    "study",
    "varY",
    "nSample",
    "ldReference",
    "region",
    "regionGiven",
    "qcAudit"
)

# @noRd
.SSQC_QTL_COLUMNS <- c(.SSQC_SHARED_COLUMNS, "context", "trait", "traitPos")

# @noRd
.SSQC_GWAS_COLUMNS <- c(.SSQC_SHARED_COLUMNS, "nCase", "nControl", "blockId")

# @noRd
.ssqcExtraCols <- function(sumstats, known) {
    extra <- setdiff(.tupleColumnNames(sumstats), known)
    set_names(map(extra, .tupleColumn, x = sumstats), extra)
}

# @noRd
.ssqcColumnOrNull <- function(sumstats, nm) {
    if (!is_in(nm, .tupleColumnNames(sumstats))) {
        return(NULL)
    }
    .tupleColumn(sumstats, nm)
}

# @noRd
.ssqcNumericColumn <- function(sumstats, nm) {
    if (!is_in(nm, .tupleColumnNames(sumstats))) {
        return(NULL)
    }
    as.numeric(.tupleColumn(sumstats, nm))
}

#' Run QC on a SumStats Collection
#'
#' Applies a single QC pass to a \code{QtlSumStats} or \code{GwasSumStats}
#' collection: per-row sanity checks via \code{.applySanityChecks} (drop rows
#' with out-of-range / zero P, BETA == 0, SE <= 0, NA in vital columns; clamp
#' tiny P; normalize CHR; coerce signed columns to numeric), variant-content
#' filters (MAF / INFO / N) via \code{.applyContentFilters}, optional
#' \code{skipRegion} drop, optional PIP screen, panel-vs-sumstats allele
#' harmonization against the \code{ldReference} via \code{harmonizeAlleles}
#' (which handles indels, strand-ambiguous variants, sign / strand flips, and
#' duplicate drops in a single sweep), optional SLALOM/DENTIST LD-mismatch QC,
#' and optional RAISS imputation. No Bioconductor genome / dbSNP packages
#' required.
#'
#' The returned collection records the options QC ran with in its
#' \code{qcOptions} slot, and what QC did to each row (variant counts, drop
#' counts at each step, which filters fired, etc.) in its \code{qcAudit}
#' column, read as a table with \code{\link{qcInfo}}. Fine-mapping and
#' TWAS-weights pipelines reject SumStats inputs where
#' \code{length(qcOptions(x)) == 0L}.
#'
#' Column-availability error contract: a non-zero \code{infoCutoff} requires
#' every entry to carry an \code{INFO} column, and a non-zero \code{nCutoff}
#' requires \code{N}; a missing column with a non-zero cutoff is a hard error.
#' \code{mafCutoff} is exempt --- a study with no frequency column is still
#' filtered against the reference panel (see below).
#'
#' @section Panel filters: \code{mafCutoff} / \code{macCutoff} /
#'   \code{imissCutoff} are measured against the \strong{LD reference panel},
#'   the same way \code{\link{fineMappingPipeline}} and
#'   \code{\link{twasWeightsPipeline}} measure them on the RSS path, so one
#'   number means the same thing wherever it is set. MAC is converted to a MAF
#'   equivalent and the stricter of the two applies, matching
#'   \code{\link{GenotypeFilterParam}}. Defaults (\code{0}, \code{0},
#'   \code{1}) filter nothing.
#'
#'   The panel is filtered once, before any entry is harmonized, so a variant
#'   the panel cannot support is absent from the LD used by kriging,
#'   SLALOM/DENTIST and RAISS, and is never an imputation target. Because
#'   harmonization drops summary-statistic variants the panel does not cover,
#'   this also \strong{discards observed variants} --- the intended behaviour:
#'   a variant whose LD cannot be estimated is not usable downstream.
#'
#'   \code{mafCutoff} additionally filters the summary statistics' own
#'   \code{AF} / \code{MAF} / \code{FRQ} column when one is present, so a
#'   variant common in the panel but rare in the study is dropped as well.
#'
#' @param sumstats A \code{QtlSumStats} or \code{GwasSumStats} collection.
#' @param sumstatsFilterParam Row filters applied to the summary statistics
#'   themselves, built with \code{\link{SumstatsFilterParam}}:
#'   \code{removeIndels} drops indels during panel harmonization,
#'   \code{removeStrandAmbiguous} drops A/T and C/G variants,
#'   \code{infoCutoff} is an INFO-score floor (requiring an \code{INFO}
#'   column when non-zero), and \code{nCutoff} drops variants whose \code{N}
#'   is more than that many median-absolute-deviations from the median (0
#'   disables it).
#' @param panelFilterParam LD-reference-panel filters, built with
#'   \code{\link{PanelFilterParam}}. \code{macCutoff} is converted to a MAF
#'   equivalent using \code{macCutoff / (2 * nSamples)} and the stricter of
#'   it and \code{mafCutoff} applies; \code{nSamples} counts the panel's real
#'   samples, which for a stochastic LD sketch are those its \code{.afreq}
#'   records rather than its pseudo-samples. \code{imissCutoff} is a
#'   per-variant missingness ceiling. See the panel-filters section below.
#'
#'   \code{mafCutoff} is the one setting that reaches past the panel: it is
#'   measured wherever it can be, against the summary statistics' own
#'   \code{AF} / \code{MAF} / \code{FRQ} column when one is present
#'   \emph{and} against the panel. It has always been a single knob, so it
#'   lives here rather than being split into two that could disagree.
#' @param keepVariants Optional character vector of variant IDs (SNP column) to
#'   retain prior to harmonization.
#' @param skipRegion Optional \code{GRanges} of regions to drop.
#' @param signalScreenParam Whether to skip an entry before fitting it,
#'   built with
#'   \code{\link{SignalScreenParam}} --- the same bundle
#'   \code{\link{fineMappingPipeline}} and \code{\link{colocboostPipeline}}
#'   take. \code{pip} runs an LD-independent single-effect SER screen and
#'   skips the entry if no PIP exceeds the cutoff (\code{< 0} resolves to
#'   \code{3 / nVariants}); \code{absZ} skips when \code{max(abs(Z))} does
#'   not exceed the cutoff, with no model fit; \code{bf} and \code{logBf}
#'   screen the largest per-variant single-effect Bayes factor from that same
#'   \code{susie_ser} fit, on the raw and log scales respectively. One metric
#'   at a time --- the constructor rejects a conflicting pair. Unset (the
#'   default) screens nothing.
#' @param ldMismatchQcMethod Which LD-mismatch check to run: \code{"none"}
#'   (default), \code{"slalom"} or \code{"dentist"}, or the matching
#'   constructor -- \code{\link{SlalomParam}} / \code{\link{DentistParam}} --
#'   to configure it at the same time. Exactly one engine is selected either
#'   way, since the constructor carries the choice.
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
#'   \code{nCase + nControl}) with no override.
#'   \code{qcOptions(x)$options$effectiveN} records the setting and each
#'   row's \code{nSource} (a \code{\link{qcInfo}} column) is one of
#'   \code{"effective"}, \code{"column"}, \code{"total"}, or \code{NA}.
#' @param impute Logical (length 1). Run RAISS imputation against the
#'   \code{ldReference}'s genotype panel as the last QC step, adding z-scores
#'   for panel variants the input lacks (see \code{raissParam}). The audit
#'   records how many variants were imputed. Default \code{FALSE}.
#' @param raissParam RAISS settings, built with \code{\link{RaissParam}}.
#'   A bare list is refused, since it cannot be checked.
#'
#'   RAISS imputation scopes its reference panel to the analysis-region window
#'   (so a per-chromosome / genome-wide \code{ldReference} does not materialize
#'   its full dosage); \code{flank} widens that window by the given number of
#'   base pairs on each side to retain LD context for edge variants, though
#'   never into the stretch of another row of the same study, nor past a
#'   block range that was given: a study split into LD blocks imputes each
#'   block once, in its own row.
#'
#'   \code{mafCutoff}, \code{macCutoff} and \code{imissCutoff} bound which
#'   panel variants RAISS will attempt to impute. Without them every rare
#'   variant in the window of a large LD sketch becomes an imputation target,
#'   which is slow and of little value when the study is far smaller than the
#'   panel. The stricter of \code{mafCutoff} and
#'   \code{macCutoff / (2 * nSamples)} applies, matching
#'   \code{\link{QtlDataset}}.
#'
#'   These are a \strong{further} tightening applied to imputation targets
#'   only, on top of the panel filter the top-level cutoffs already applied:
#'   the panel handed to RAISS holds nothing below those. Use them to impute
#'   only common variants from a panel deliberately kept wider than that ---
#'   e.g. \code{mafCutoff = 0.001} with
#'   \code{raissParam = RaissParam(mafCutoff = 0.01)}. A variant present in
#'   the sumstats is never dropped here: RAISS derives its LD basis from those
#'   same panel rows.
#'
#'   \code{lamb}, \code{svdTol}, \code{r2Threshold} and \code{minimumLd}
#'   are forwarded to \code{\link{raiss}} itself.
#' @param matchMinProp Minimum proportion of LD panel variants that must be
#'   matched by the sumstats; default 0.
#' @param sumstatsCleaningParam What makes a row a well-formed record, built
#'   with \code{\link{SumstatsCleaningParam}}: \code{coerceNumeric},
#'   \code{normalizeChr} / \code{dropNonstandardChr}, \code{dropMissData},
#'   \code{dropPOutOfRange}, \code{clampSmallP} / \code{smallPFloor},
#'   \code{dropZeroEffect} and \code{dropNonpositiveSe}. Applied before any
#'   filter or screen. Distinct from \code{sumstatsFilterParam}, which applies
#'   quality thresholds to rows that are already well formed.
#' @return A new \code{QtlSumStats} / \code{GwasSumStats} with cleaned entries,
#'   \code{qcOptions} filled and each row's QC audit attached.
#' @examples
#' data(gwasSumStatsS4Example)
#' summaryStatsQc(gwasSumStatsS4Example)
#' @export
summaryStatsQc <- function(
    sumstats,
    sumstatsFilterParam = SumstatsFilterParam(),
    panelFilterParam = PanelFilterParam(),
    keepVariants = NULL,
    skipRegion = NULL,
    signalScreenParam = SignalScreenParam(),
    ldMismatchQcMethod = c("none", "slalom", "dentist"),
    alleleFlipKriging = FALSE,
    effectiveN = TRUE,
    impute = FALSE,
    raissParam = RaissParam(),
    matchMinProp = 0,
    sumstatsCleaningParam = SumstatsCleaningParam()
) {
    r <- .ssqcResolveInputs(
        sumstats,
        sumstatsFilterParam,
        panelFilterParam,
        signalScreenParam,
        sumstatsCleaningParam,
        raissParam,
        ldMismatchQcMethod
    )
    opts <- .ssqcBuildOpts(
        sumstatsFilterParam = sumstatsFilterParam,
        panelFilterParam = panelFilterParam,
        skipRegion = skipRegion,
        ldMismatchQcMethod = r$ldMismatchQcMethod,
        alleleFlipKriging = alleleFlipKriging,
        effectiveN = effectiveN,
        impute = impute,
        raissParam = r$raissParam,
        matchMinProp = matchMinProp,
        sumstatsCleaningParam = r$cleaning,
        keepVariants = keepVariants,
        signalScreenParam = signalScreenParam
    )
    res <- .ssqcRunEntries(sumstats, opts)
    .ssqcAssemble(sumstats, res, opts, signalScreenParam)
}

# Validate every bundle and resolve the three that have a resolved form:
# the LD-mismatch choice (a name or the engine's own constructor -- the
# record carries both the choice and its options, so they cannot
# disagree), the cleaning defaults, and the imputation arguments.
# @noRd
.ssqcResolveInputs <- function(
    sumstats,
    sumstatsFilterParam,
    panelFilterParam,
    signalScreenParam,
    sumstatsCleaningParam,
    raissParam,
    ldMismatchQcMethod
) {
    .assertMethodParam(
        sumstatsFilterParam,
        "SumstatsFilterParam"
    )
    .assertMethodParam(panelFilterParam, "PanelFilterParam")
    .assertMethodParam(signalScreenParam, "SignalScreenParam")
    .assertMethodParam(raissParam, "RaissParam")
    .ssqcCheckEntries(sumstats, sumstatsFilterParam$infoCutoff)
    .ssqcCheckPanelCutoffs(panelFilterParam)
    list(
        ldMismatchQcMethod = .resolveLdMismatchChoice(ldMismatchQcMethod),
        cleaning = .sumstatsCleaningResolve(sumstatsCleaningParam),
        raissParam = .ssqcResolveImputeArgs(raissParam)
    )
}

# The QC record and the rebuilt object around the filtered entries.
# @noRd
.ssqcAssemble <- function(sumstats, res, opts, signalScreenParam) {
    qcOptions <- .ssqcBuildQcOptions(
        sumstatsFilterParam = opts$sumstatsFilterParam,
        panelFilterParam = opts$panelFilterParam,
        signalScreenParam = signalScreenParam,
        ldMismatchQcMethod = opts$ldMismatchQcMethod,
        alleleFlipKriging = opts$alleleFlipKriging,
        effectiveN = opts$effectiveN,
        impute = opts$impute,
        sumstatsCleaningParam = opts$sumstatsCleaningParam
    )
    # From the PRUNED panels, not the input's: narrowing an original again
    # would hand back a sketch the panel filter never reached. Each is cut to
    # the variants of the rows that use it.
    trimmed <- imap(
        res$ldReferences,
        .ssqcTrimToRows,
        entries = res$newEntries,
        keys = res$keys
    )
    perRow <- map(res$keys, .ldRefForKey, refs = trimmed)
    .ssqcRebuild(sumstats, res$newEntries, perRow, qcOptions, res$entryAudits)
}

# One LD reference cut to exactly the variants of the rows keyed to it.
# @noRd
.ssqcTrimToRows <- function(ldReference, key, entries, keys) {
    .subsetLdRefToIds(ldReference, entries[is_in(keys, key)])
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

# TRUE when credible set row `i` is a tagged (redundant) set.
# @noRd
.autoDecisionTagged <- function(i, df, highCorrCols) {
    if (df$top_cs[i]) {
        return(FALSE)
    }
    if (df$p_value[i] > 1e-4) {
        return(TRUE)
    }
    # A CS's own column -- cs_corr_<i> in extractCsInfo()'s layout -- holds its
    # self-correlation of 1, which says nothing about the other sets.
    otherCols <- setdiff(highCorrCols, str_c("cs_corr_", i))
    if (length(otherCols) == 0) {
        return(FALSE)
    }
    rowVals <- df |>
        filter(row_number() == i) |>
        select(all_of(otherCols))
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
    variants <- colnames(ld) %||% character(0)
    tibble(variants = variants)
}

# The is-NA mask of column `col` in `df` (for the vital-column drop reduce).
# @noRd
.scColIsNa <- function(col, df) {
    is.na(df[[col]])
}

# One GRanges entry's chromosome vector.
# @noRd
.entryChrom <- function(gr) {
    as.character(GenomicRanges::seqnames(gr))
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

#' @rdname SumstatsFilterParam
#' @aliases SumstatsFilterParam-class
#' @exportClass SumstatsFilterParam
setClass(
    "SumstatsFilterParam",
    contains = "MethodParam",
    slots = c(
        removeIndels = "logical",
        removeStrandAmbiguous = "logical",
        infoCutoff = "numeric",
        nCutoff = "numeric"
    )
)

#' @title Summary-Statistic Row Filtering Options
#' @description Which rows of a summary-statistics table survive QC, as one
#'   checked bundle. These act on the sumstats themselves, not on any genotype
#'   matrix, which is why they are separate from
#'   \code{\link{GenotypeFilterParam}}.
#' @param removeIndels Logical. Drop insertions and deletions. Default
#'   \code{FALSE}.
#' @param removeStrandAmbiguous Logical. Drop strand-ambiguous variants (A/T
#'   and C/G). Default \code{TRUE}.
#' @param infoCutoff Imputation-INFO floor. Default \code{0}.
#' @param nCutoff Per-variant sample-size floor. Default \code{5}.
#' @return A \code{SumstatsFilterParam} object, a \code{\link{MethodParam}}.
#' @seealso \code{\link{GenotypeFilterParam}}, \code{\link{PanelFilterParam}}
#' @examples
#' SumstatsFilterParam(removeIndels = TRUE, infoCutoff = 0.8)
#' @export
SumstatsFilterParam <- function(
    removeIndels = FALSE,
    removeStrandAmbiguous = TRUE,
    infoCutoff = 0,
    nCutoff = 5
) {
    new(
        "SumstatsFilterParam",
        removeIndels = removeIndels,
        removeStrandAmbiguous = removeStrandAmbiguous,
        infoCutoff = infoCutoff,
        nCutoff = nCutoff
    )
}

#' @rdname SumstatsCleaningParam
#' @aliases SumstatsCleaningParam-class
#' @exportClass SumstatsCleaningParam
setClass(
    "SumstatsCleaningParam",
    contains = "MethodParam",
    slots = c(
        coerceNumeric = "logical",
        normalizeChr = "logical",
        dropNonstandardChr = "logical",
        dropMissData = "logical",
        dropPOutOfRange = "logical",
        clampSmallP = "logical",
        smallPFloor = "numeric",
        dropZeroEffect = "logical",
        dropNonpositiveSe = "logical"
    )
)

#' @title Summary-Statistics Cleaning Settings
#' @description What makes a summary-statistics row a well-formed record:
#'   type coercion, chromosome-label normalization, dropping malformed or
#'   impossible rows, and flooring underflowed p-values.
#'
#'   Separate from \code{\link{SumstatsFilterParam}} on purpose. That one
#'   applies quality thresholds to \emph{valid} records (INFO, N, indels,
#'   strand-ambiguous variants); this one decides whether a row is a valid
#'   record at all. All of it is applied by \code{\link{summaryStatsQc}}
#'   before any filter or screen runs.
#' @param coerceNumeric Logical. Coerce the signed columns
#'   (Z/BETA/SE/OR/LOG_ODDS/SIGNED_SUMSTAT/P/MAF/FRQ/INFO/N) to numeric.
#'   Default \code{TRUE}.
#' @param normalizeChr Logical. Strip the \code{"chr"} prefix, uppercase the
#'   chromosome label, and map 23->X, 24->Y, M->MT. Default \code{TRUE}.
#' @param dropNonstandardChr Logical. Drop variants whose CHR (after
#'   normalization) is outside 1..22, X, Y, MT. Default \code{TRUE}.
#' @param dropMissData Logical. Drop rows with NA in any vital column (chrom,
#'   pos, A1, A2, and at least one of Z / BETA). Default \code{TRUE}.
#' @param dropPOutOfRange Logical. Drop rows where \code{P < 0} or
#'   \code{P > 1}. Default \code{TRUE}.
#' @param clampSmallP Logical. Floor non-negative P values below
#'   \code{smallPFloor} to \code{smallPFloor}, so \code{-log10(P)} stays
#'   finite. Applied to both input and Z-derived P values. Default
#'   \code{TRUE}.
#' @param smallPFloor Numeric (length 1). The floor \code{clampSmallP}
#'   applies. Default \code{5e-324} (R's smallest positive double).
#' @param dropZeroEffect Logical. Drop rows where any effect column is exactly
#'   0 (\code{BETA}, \code{LOG_ODDS}, \code{SIGNED_SUMSTAT}) or \code{OR}
#'   is exactly 1. Default \code{TRUE}.
#' @param dropNonpositiveSe Logical. Drop rows where \code{SE <= 0}. Default
#'   \code{TRUE}.
#' @return A \code{SumstatsCleaningParam} object, a \code{\link{MethodParam}}.
#' @examples
#' SumstatsCleaningParam(clampSmallP = FALSE, dropZeroEffect = FALSE)
#' @export
SumstatsCleaningParam <- function(
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
    new(
        "SumstatsCleaningParam",
        coerceNumeric = coerceNumeric,
        normalizeChr = normalizeChr,
        dropNonstandardChr = dropNonstandardChr,
        dropMissData = dropMissData,
        dropPOutOfRange = dropPOutOfRange,
        clampSmallP = clampSmallP,
        smallPFloor = smallPFloor,
        dropZeroEffect = dropZeroEffect,
        dropNonpositiveSe = dropNonpositiveSe
    )
}

# The cleaning settings as a plain list with every default applied.
# .assertMethodOptions accepts `list()` as "no options", and .applySanityChecks
# tests these with `if (!flag)`, so a missing field must not arrive as NULL.
# @noRd
.sumstatsCleaningResolve <- function(sumstatsCleaningParam) {
    .assertMethodParam(
        sumstatsCleaningParam,
        "SumstatsCleaningParam"
    )
    list(
        coerceNumeric = sumstatsCleaningParam$coerceNumeric %||% TRUE,
        normalizeChr = sumstatsCleaningParam$normalizeChr %||% TRUE,
        dropNonstandardChr = sumstatsCleaningParam$dropNonstandardChr %||% TRUE,
        dropMissData = sumstatsCleaningParam$dropMissData %||% TRUE,
        dropPOutOfRange = sumstatsCleaningParam$dropPOutOfRange %||% TRUE,
        clampSmallP = sumstatsCleaningParam$clampSmallP %||% TRUE,
        smallPFloor = sumstatsCleaningParam$smallPFloor %||% 5e-324,
        dropZeroEffect = sumstatsCleaningParam$dropZeroEffect %||% TRUE,
        dropNonpositiveSe = sumstatsCleaningParam$dropNonpositiveSe %||% TRUE
    )
}
