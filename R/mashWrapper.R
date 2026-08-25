# Filter rows of a z-score matrix by significance p-value cutoff.
# Returns integer indices of rows where any |z| exceeds the threshold.
# @noRd
filterBySignificance <- function(zMatrix, sigPCutoff) {
    zThreshold <- sqrt(stats::qchisq(sigPCutoff, df = 1, lower.tail = FALSE))
    which(apply(zMatrix, 1, .mashRowExceeds, zThreshold = zThreshold))
}

# Coerce to a numeric matrix and replace every NaN / Inf / NA cell with
# `replaceWith`. Accepts a data.frame or a matrix and always returns a matrix:
# mash operates on matrices, so the cleaning is a direct matrix op rather than a
# data.frame round-trip.
# @noRd
.mashReplaceValues <- function(x, replaceWith) {
    m <- as.matrix(x)
    storage.mode(m) <- "double"
    m[is.nan(m) | is.infinite(m) | is.na(m)] <- replaceWith
    m
}

# Coerce z-scores to a matrix (NaN/Inf/NA -> 0) and, when a missing-rate
# threshold is given, drop rows falling below it.
# @noRd
.mashProcessZ <- function(zData, filterByMissingRate) {
    zData <- .mashReplaceValues(zData, 0)

    if (!is.null(filterByMissingRate)) {
        proportionNonzero <- apply(zData, 1, .mashRowNonzeroRate)
        zData <- zData[proportionNonzero >= filterByMissingRate, , drop = FALSE]
    }

    return(zData)
}

#' Filter invalid summary statistics for mash input
#'
#' Assemble and clean a per-condition effect-size / z-score matrix for
#' \code{mashr}, dropping variants with invalid or insufficiently-observed
#' statistics.
#'
#' @param datList A named list of summary-statistic data frames / matrices.
#' @param bhat Optional name of the effect-size element in \code{datList}.
#' @param sbhat Optional name of the standard-error element in \code{datList}.
#' @param z Optional name of the z-score element in \code{datList}.
#' @param btoz Logical. If \code{TRUE}, derive z-scores from
#'   \code{bhat}/\code{sbhat}.
#' @param sigPCutoff Numeric. Significance p-value cutoff for selecting strong
#'   signals. Default \code{1e-6}.
#' @param filterByMissingRate Numeric in [0, 1]. Drop variants observed in fewer
#'   than this fraction of conditions. Default \code{0.2}.
#' @return A cleaned list of summary-statistic matrices suitable for mash.
#' @importFrom vroom vroom
#' @examples
#' datList <- list(strong = list(z = matrix(rnorm(9), 3, 3)))
#' filterInvalidSummaryStat(datList)
#' @export
filterInvalidSummaryStat <- function(
    datList,
    bhat = NULL,
    sbhat = NULL,
    z = NULL,
    btoz = FALSE,
    sigPCutoff = 1E-6,
    filterByMissingRate = 0.2
) {
    if (
        !is.null(bhat) &&
            !is.null(sbhat) &&
            all(is_in(c(bhat, sbhat), names(datList)))
    ) {
        datList <- .mashFilterBhatSbhat(
            datList,
            bhat,
            sbhat,
            filterByMissingRate
        )
    }
    if (btoz) {
        datList <- .mashFilterBtoz(datList, bhat, sbhat, sigPCutoff)
    }
    if (!is.null(z)) {
        datList <- .mashFilterZ(datList, filterByMissingRate, sigPCutoff)
    }
    datList
}

# Reset invalid bhat/sbhat cells (bhat -> 0, sbhat -> 1000) and, when a
# null/random partition is present, drop variants below `filterByMissingRate`
# non-missing.
# @noRd
.mashFilterBhatSbhat <- function(datList, bhat, sbhat, filterByMissingRate) {
    if (is.null(datList[[bhat]]) || is.null(datList[[sbhat]])) {
        return(datList)
    }
    datList[[bhat]] <- .mashReplaceValues(datList[[bhat]], 0)
    datList[[sbhat]] <- .mashReplaceValues(datList[[sbhat]], 1000)
    hasNullOrRandom <- is_in("null.b", names(datList)) ||
        is_in("random.b", names(datList))
    if (!hasNullOrRandom || is.null(filterByMissingRate)) {
        return(datList)
    }
    proportionNonzero <- apply(datList[[bhat]], 1, .mashRowNonzeroRate)
    keep <- proportionNonzero >= filterByMissingRate
    datList[[bhat]] <- datList[[bhat]][keep, ]
    datList[[sbhat]] <- datList[[sbhat]][keep, ]
    datList
}

# Derive z = bhat / sbhat (into a `<condition>.z` or `z` slot) and apply the
# significance cutoff to strong signals.
# @noRd
.mashFilterBtoz <- function(datList, bhat, sbhat, sigPCutoff) {
    if (any(str_detect(bhat, "\\.b$")) || any(str_detect(sbhat, "\\.s$"))) {
        zName <- str_c(str_remove(bhat, "\\.b$"), ".z")
        if (!is.null(datList[[bhat]]) && !is.null(datList[[sbhat]])) {
            datList[[zName]] <- as.matrix(datList[[bhat]] / datList[[sbhat]])
        } else {
            datList[zName] <- list(NULL)
        }
    } else if (!is.null(datList[[bhat]]) && !is.null(datList[[sbhat]])) {
        datList[["z"]] <- as.matrix(datList[[bhat]] / datList[[sbhat]])
    } else {
        datList["z"] <- list(NULL)
    }
    if (is_in("strong.z", names(datList)) && !is.null(sigPCutoff)) {
        keepIndex <- filterBySignificance(datList$strong.z, sigPCutoff)
        datList[["strong.z"]] <- datList$strong.z[keepIndex, ]
        datList[["strong.b"]] <- datList$strong.b[keepIndex, ]
        datList[["strong.s"]] <- datList$strong.s[keepIndex, ]
    }
    datList
}

# Process each partition's z-matrix (missing-rate filter) and apply the
# significance cutoff to strong z-scores.
# @noRd
.mashFilterZ <- function(datList, filterByMissingRate, sigPCutoff) {
    for (comp in c("strong", "random", "null")) {
        if (!is.null(datList[[comp]]) && !is.null(datList[[comp]]$z)) {
            datList[[comp]]$z <- .mashProcessZ(
                datList[[comp]]$z,
                filterByMissingRate
            )
        }
    }
    if (
        !is.null(datList$strong) &&
            !is.null(datList$strong$z) &&
            !is.null(sigPCutoff)
    ) {
        keepIndex <- filterBySignificance(datList$strong$z, sigPCutoff)
        datList$strong$z <- datList$strong$z[keepIndex, , drop = FALSE]
    }
    datList
}

#' Filter conditions from mash prior mixture components
#'
#' Drop the conditions not in \code{conditionsToKeep} from each prior covariance
#' matrix in \code{U}, optionally removing components whose weight is below
#' \code{wCutoff}.
#'
#' @param conditionsToKeep Character vector of condition names to retain.
#' @param U Named list of prior covariance matrices (one per mixture component).
#' @param w Optional numeric vector of mixture weights aligned to \code{U}.
#' @param wCutoff Numeric. Drop components with weight below this. Default
#'   \code{1e-4}.
#' @return A list with the filtered \code{U} (and \code{w} when supplied).
#' @importFrom purrr keep
#' @examples
#' conditionsToKeep <- c("cond1", "cond2")
#' cn <- c("cond1", "cond2", "cond3")
#' U <- list(shared = diag(3), corr = matrix(0.3, 3, 3) + diag(0.7, 3))
#' U <- lapply(U, function(m) {
#'   dimnames(m) <- list(cn, cn)
#'   m
#' })
#' filterMixtureComponents(conditionsToKeep = conditionsToKeep, U = U)
#' @export
filterMixtureComponents <- function(
    conditionsToKeep,
    U,
    w = NULL,
    wCutoff = 1e-04
) {
    conditionsToFilter <- setdiff(colnames(U[[1]]), conditionsToKeep)
    sumW <- sum(w)
    U <- .mashSubsetU(U, conditionsToKeep)
    # Drop all-zero matrices, then those below the weight cutoff.
    keepNames <- names(keep(U, .mashMatrixNonzero))
    if (!is.null(w)) {
        keepNames <- intersect(keepNames, names(w[w >= wCutoff]))
    }
    U <- U[keepNames]
    if (!is.null(w)) {
        w <- w[keepNames]
    }
    # Manually remove the U components driven by non-relevant contexts: the EM
    # can leave tiny non-zero diagonals, so all-zero removal alone won't drop
    # them, yet real diagonal signal must be kept.
    U[conditionsToFilter] <- NULL
    w <- w[!is_in(names(w), conditionsToFilter)]
    # Rescale the surviving weights back to the original total.
    w <- (w / sum(w)) * sumW
    msg <- glue("{length(U)} components of matrices remained after filtering.")
    inform(msg)
    list(U = U, w = w)
}

# Subset every U matrix to the kept conditions (erroring if a matrix lacks one).
# @noRd
.mashSubsetU <- function(U, conditionsToKeep) {
    map(U, .mashSubsetMatrix, conditionsToKeep = conditionsToKeep)
}


# Draw the random + null sub-samples used to estimate the null correlation.
# @noRd
.mashExtractOneData <- function(dat, nRandom, nNull) {
    if (is.null(dat)) {
        return(NULL)
    }
    if (is_in("z", names(dat))) {
        absZ <- abs(dat$z)
        zData <- dat$z
    } else {
        absZ <- abs(dat$bhat / dat$sbhat)
        zData <- NULL
    }
    random <- .mashSampleSubset(dat, zData, seq_len(nrow(absZ)), nRandom)
    null <- .mashSampleNull(dat, zData, absZ, nNull)
    list(random = random, null = null)
}

# Sample up to `n` rows from `poolIdx` and return them as a z (or bhat/sbhat)
# list, matching the source scale.
# @noRd
.mashSampleSubset <- function(dat, zData, poolIdx, n) {
    idx <- sample(poolIdx, min(n, length(poolIdx)), replace = FALSE)
    if (!is.null(zData)) {
        list(z = zData[idx, , drop = FALSE])
    } else {
        list(
            bhat = dat$bhat[idx, , drop = FALSE],
            sbhat = dat$sbhat[idx, , drop = FALSE]
        )
    }
}

# Null subset: variants with max|z| < 2. Empty (with a warning) when there are
# none, or too few to estimate the null correlation.
# @noRd
.mashSampleNull <- function(dat, zData, absZ, nNull) {
    nullId <- which(apply(absZ, 1, max) < 2)
    if (length(nullId) == 0) {
        msg <- glue(
            "no variants are included in the null dataset because absZ > 2 ",
            "for all variants in {dat$region %||% ''}"
        )
        warn(msg)
        return(list())
    }
    if (length(nullId) < ncol(absZ)) {
        msg <- glue(
            "not enough null data to estimate null correlation in ",
            "{dat$region %||% ''}"
        )
        warn(msg)
        return(list())
    }
    .mashSampleSubset(dat, zData, nullId, nNull)
}

#' Sample random and null variant subsets for mash
#'
#' Draw a random subset and a null (non-significant) subset of rows from a mash
#' data list, used to fit the mash prior and estimate the null correlation.
#'
#' @param dat A mash data list with \code{random} and \code{null} components.
#' @param nRandom Integer. Number of random rows to sample.
#' @param nNull Integer. Number of null rows to sample.
#' @param excludeCondition Optional character vector of conditions to exclude.
#' @param seed Optional integer random seed; \code{NULL} leaves the RNG
#'   unchanged.
#' @return A list with sampled \code{random} and \code{null} matrices.
#' @examples
#' cond <- c("brain", "blood", "muscle")
#' p <- 8
#' bhat <- matrix(rnorm(p * 3), p, 3,
#'   dimnames = list(sprintf("chr1:%d:A:G", 100L * (1:p)), cond))
#' sbhat <- matrix(abs(rnorm(p * 3)) + 0.1, p, 3,
#'   dimnames = list(sprintf("chr1:%d:A:G", 100L * (1:p)), cond))
#' dat <- list(bhat = bhat, sbhat = sbhat, Z = bhat / sbhat,
#'   snp = sprintf("chr1:%d:A:G", 100L * (1:p)))
#' mashRandNullSample(dat, nRandom = 2L, nNull = 2L,
#'   excludeCondition = character())
#' @export
mashRandNullSample <- function(
    dat,
    nRandom,
    nNull,
    excludeCondition,
    seed = NULL
) {
    if (!is.null(seed)) {
        withr::local_seed(seed)
    }

    if (length(excludeCondition) > 0) {
        colsToCheck <- if (is_in("z", names(dat))) "z" else "bhat"
        if (!all(is_in(excludeCondition, colnames(dat[[colsToCheck]])))) {
            msg <- glue(
                "Error: excludeCondition are not present in ",
                "{dat$region %||% ''}"
            )
            abort(msg)
        }
        for (key in intersect(names(dat), c("z", "bhat", "sbhat"))) {
            keep <- setdiff(colnames(dat[[key]]), excludeCondition)
            dat[[key]] <- dat[[key]][, keep, drop = FALSE]
        }
    }

    result <- .mashExtractOneData(dat, nRandom, nNull)
    return(result)
}

#' Merge two mash data lists
#'
#' Row-bind the components of two mash data lists, returning the non-empty one
#' when the other is empty.
#'
#' @param resData The accumulated mash data list (may be empty).
#' @param oneData The mash data list to merge in.
#' @return The merged mash data list.
#' @examples
#' # Each object's variants must be uniquely keyed (row names); the two
#' # objects share the same conditions (columns), which are aligned by name.
#' a <- list(strong = list(z = matrix(rnorm(9), 3, 3,
#'   dimnames = list(
#'     c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G"),
#'     c("t1", "t2", "t3")))))
#' b <- list(strong = list(z = matrix(rnorm(9), 3, 3,
#'   dimnames = list(
#'     c("chr1:400:A:G", "chr1:500:A:G", "chr1:600:A:G"),
#'     c("t1", "t2", "t3")))))
#' mergeMashData(a, b)
#' @export
mergeMashData <- function(resData, oneData) {
    if (length(resData) == 0 || is.null(resData)) {
        return(oneData)
    }
    if (length(oneData) == 0 || is.null(oneData)) {
        return(resData)
    }

    combinedData <- map(
        names(oneData),
        .mashCombineDatum,
        oneData = oneData,
        resData = resData
    )
    names(combinedData) <- names(oneData)
    return(combinedData)
}

# Build variants x conditions (Bhat, Shat) matrices for ONE object plus the
# row indices of its "strong" variants. Class dispatch:
#   QtlSumStats / GwasSumStats -> .mashSumStatsToMatrices; strong = the single
#       most significant variant (max|z|) per condition, unioned.
#   FineMappingResultBase      -> pivot getMarginalEffects() into a
#       variants x context (beta, se) pair; strong = the lead (max PIP) variant
#       of each credible set in each condition (getCs()), unioned. Conditions
#       with no credible set contribute no strong variant.
# For z-scale QtlSumStats the returned Shat is 1, so downstream code that forms
# z = b / s recovers the z-scores uniformly across both scales.
# @noRd
.mashObjectMatrices <- function(obj, inputScale, coverage) {
    if (methods::is(obj, "QtlSumStats") || methods::is(obj, "GwasSumStats")) {
        return(.mashSumStatsMatrices(obj, inputScale))
    }
    if (methods::is(obj, "FineMappingResultBase")) {
        return(.mashFmrMatrices(obj, coverage))
    }
    msg <- glue(
        "mashInput: each element of `objects` must be a QtlSumStats, ",
        "GwasSumStats, or FineMappingResult; got ",
        "{str_flatten(class(obj), '/')}."
    )
    abort(msg)
}

# (Bhat, Shat, strongRows) from a SumStats object; strong = the max|z| variant
# per condition column, unioned.
# @noRd
.mashSumStatsMatrices <- function(obj, inputScale) {
    mats <- .mashSumStatsToMatrices(obj, "mash input", inputScale = inputScale)
    z <- mats$b / mats$s
    strongRows <- sort(unique(apply(abs(z), 2L, which.max)))
    list(b = mats$b, s = mats$s, strongRows = strongRows)
}

# (Bhat, Shat, strongRows) from a FineMappingResult: pivot the marginal effects
# to a variants x contexts matrix pair; strong = each credible set's lead (max
# PIP) variant.
# @noRd
.mashFmrMatrices <- function(obj, coverage) {
    me <- getMarginalEffects(obj)
    cs <- getCs(obj, coverage = coverage)
    if (!all(is_in(c("variant_id", "context", "beta", "se"), names(me)))) {
        msg <- glue(
            "mashInput: getMarginalEffects() must return variant_id/context/",
            "beta/se columns; a FineMappingResult with >= 2 contexts is ",
            "required."
        )
        abort(msg)
    }
    pinned <- .mashFmrMethodPin(me, cs)
    me <- pinned$me
    cs <- pinned$cs
    contexts <- unique(me$context)
    variants <- unique(me$variant_id)
    b <- matrix(
        NA_real_,
        length(variants),
        length(contexts),
        dimnames = list(variants, contexts)
    )
    s <- b
    cell <- cbind(match(me$variant_id, variants), match(me$context, contexts))
    b[cell] <- me$beta
    s[cell] <- me$se
    list(b = b, s = s, strongRows = .mashFmrStrongRows(cs, variants))
}

# A multi-method FineMappingResult would duplicate (variant, context) cells;
# pin the first method so the pivot is unambiguous.
# @noRd
.mashFmrMethodPin <- function(me, cs) {
    if (!is_in("method", names(me)) || n_distinct(me$method) <= 1L) {
        return(list(me = me, cs = cs))
    }
    m1 <- me$method[[1L]]
    msg <- glue(
        "mashInput: FineMappingResult carries multiple methods; using ",
        "'{m1}'."
    )
    warn(msg)
    me <- filter(me, .data$method == m1)
    if (is_in("method", names(cs))) {
        cs <- filter(cs, .data$method == m1)
    }
    list(me = me, cs = cs)
}

# Row indices (into `variants`) of each credible set's lead (max PIP) variant.
# @noRd
.mashFmrStrongRows <- function(cs, variants) {
    strongVar <- character(0)
    csCol <- names(cs)[str_detect(names(cs), "^cs_")]
    if (nrow(cs) > 0L && length(csCol) > 0L && is_in("pip", names(cs))) {
        grp <- interaction(cs$context, cs[[csCol[[1L]]]], drop = TRUE)
        for (rows in split(seq_len(nrow(cs)), grp)) {
            strongVar <- c(
                strongVar,
                cs$variant_id[[rows[[which.max(cs$pip[rows])]]]]
            )
        }
    }
    strongRows <- sort(match(unique(strongVar), variants))
    strongRows[!is.na(strongRows)]
}

# Extract the strong / random / null partitions from ONE object as a flat
# list(strong.b, strong.s, random.b, random.s, null.b, null.s) of
# variants x conditions matrices. Random / null are drawn by the shared
# mashRandNullSample() over the object's (Bhat, Shat); strong is the
# deterministic class-specific selection from .mashObjectMatrices().
# @noRd
.mashObjectPartitions <- function(
    obj,
    nRandom,
    nNull,
    excludeCondition,
    coverage,
    inputScale,
    seed,
    independentVariants = NULL
) {
    mats <- .mashObjectMatrices(
        obj,
        inputScale = inputScale,
        coverage = coverage
    )
    keepCols <- setdiff(colnames(mats$b), excludeCondition)
    if (length(keepCols) < 2L) {
        msg <- glue(
            "mashInput: fewer than 2 conditions remain for an object (after ",
            "excludeCondition); mash operates across conditions and needs >= 2."
        )
        abort(msg)
    }
    pool <- .mashIndependentPool(mats, independentVariants)
    rn <- mashRandNullSample(
        list(bhat = pool$poolB, sbhat = pool$poolS),
        nRandom = nRandom,
        nNull = nNull,
        excludeCondition = excludeCondition,
        seed = seed
    )
    .mashPartitionOut(mats, keepCols, rn)
}

# Random / null candidate pool, optionally restricted to LD-independent
# variants (so the background carries no LD-correlated SNPs, which would bias
# Vhat / weights). Strong is always drawn from the full set elsewhere.
# @noRd
.mashIndependentPool <- function(mats, independentVariants) {
    if (is.null(independentVariants) || length(independentVariants) == 0L) {
        return(list(poolB = mats$b, poolS = mats$s))
    }
    # Rownames carry a "study::trait::" block prefix; strip to the bare variant
    # id before matching (proper chrom/pos/allele via matchVariants).
    rawIds <- str_remove(rownames(mats$b), ".*::")
    keepIdx <- matchVariants(
        rawIds,
        independentVariants,
        allowFlip = TRUE,
        removeStrandAmbiguous = FALSE
    )$idxA
    if (length(keepIdx) == 0L) {
        msg <- glue(
            "mashInput: no variants matched the independent-variant list; ",
            "the random/null background is empty for this object."
        )
        warn(msg)
    }
    list(
        poolB = mats$b[keepIdx, , drop = FALSE],
        poolS = mats$s[keepIdx, , drop = FALSE]
    )
}

# Assemble the flat strong/random/null (.b/.s) partition list, omitting any
# empty partition.
# @noRd
.mashPartitionOut <- function(mats, keepCols, rn) {
    out <- list()
    if (length(mats$strongRows) > 0L) {
        out[["strong.b"]] <- mats$b[mats$strongRows, keepCols, drop = FALSE]
        out[["strong.s"]] <- mats$s[mats$strongRows, keepCols, drop = FALSE]
    }
    if (!is.null(rn$random) && length(rn$random) > 0L) {
        out[["random.b"]] <- rn$random$bhat
        out[["random.s"]] <- rn$random$sbhat
    }
    if (!is.null(rn$null) && length(rn$null) > 0L) {
        out[["null.b"]] <- rn$null$bhat
        out[["null.s"]] <- rn$null$sbhat
    }
    out
}

#' Assemble MASH strong / random / null input from S4 objects
#'
#' Unified, S4-native replacement for the legacy
#' \code{load_multitrait_*_sumstat} + \code{mash_ran_null_sample} assembly.
#' Consumes a list of already-constructed objects (one per region) and returns
#' the flat \code{variants x conditions} matrix list consumed by the MASH
#' mixture-prior / fit / posterior steps.
#'
#' For EACH object three partitions are extracted:
#' \describe{
#'   \item{strong}{The high-signal variants (deterministic, class-specific).
#'     \code{QtlSumStats}: the single most significant variant (\eqn{\max|z|})
#'     per condition, unioned. \code{FineMappingResult}: the lead variant
#'     (\eqn{\max} PIP) of each credible set in each condition, unioned
#'     (conditions with no credible set contribute nothing).}
#'   \item{random}{\code{nRandom} variants sampled uniformly at random --
#'     represents the genome-wide mixture of effects and drives the mixture
#'     weights.}
#'   \item{null}{\code{nNull} variants sampled from those with \eqn{\max|z|<2}
#'     -- the noise floor used to estimate the residual correlation (Vhat).}
#' }
#' Random and null are selected identically for both classes
#' (\code{\link{mashRandNullSample}} over the object's \code{Bhat}/\code{Shat}).
#' Partitions are merged across objects (rownames disambiguated by region name),
#' cleaned + z-derived by \code{\link{filterInvalidSummaryStat}} (\code{btoz}),
#' and the strong \code{XtX} cross-product appended.
#'
#' @param objects A named \code{list} of \code{\link{QtlSumStats}} and/or
#'   \code{FineMappingResult} objects, one per region. Names disambiguate
#'   rownames across regions (defaults to \code{region1}, \code{region2}, ...).
#'   For \code{QtlSumStats} inputs \code{\link{summaryStatsQc}} must have been
#'   run (the matrix builder rejects un-QC'd SumStats).
#' @param nRandom,nNull Per-object random / null sample sizes (default 10 each).
#' @param excludeCondition Character vector of condition (column) names to drop.
#' @param coverage Credible-set coverage for \code{FineMappingResult} strong
#'   selection (default 0.95).
#' @param zOnly When \code{TRUE} the returned partitions carry only \code{.z}
#'   (the \code{.b}/\code{.s} matrices are dropped after z is derived).
#' @param sigPCutoff Significance cutoff applied to the strong partition
#'   (default 1e-6).
#' @param inputScale Matrix scale for \code{QtlSumStats} inputs
#'   (\code{"auto"}/\code{"beta"}/\code{"z"}); ignored for
#'   \code{FineMappingResult} (always effect-size scale).
#' @param independentVariants Optional character vector of variant ids (e.g. an
#'   LD-pruned independent SNP list). When supplied, the \emph{random} and
#'   \emph{null} background of every object is restricted to variants that match
#'   this set, so the background carries no LD-correlated SNPs (which would bias
#'   the residual correlation and the mixture weights). Matching is delegated to
#'   \code{matchVariants()} (chrom/pos/allele aware, ref/alt flips tolerated),
#'   \emph{not} a raw string compare, so a chr-prefix / separator / allele-order
#'   difference still matches. The \emph{strong} partition is never filtered.
#' @param seed RNG seed for the random / null sampling (default 999).
#'
#' @return A flat \code{list}: \code{strong.b}, \code{strong.s},
#'   \code{strong.z}, \code{random.*}, \code{null.*} (each a \code{variants x
#'   conditions} matrix) and \code{XtX} (a \code{conditions x conditions}
#'   matrix). The \code{.b} / \code{.s} matrices are omitted when \code{zOnly =
#'   TRUE}.
#' @seealso \code{\link{mashRandNullSample}}, \code{\link{mergeMashData}},
#'   \code{\link{filterInvalidSummaryStat}}
#' @examples
#' data(qtlSumStatsMulticontextExample)
#' ss <- qtlSumStatsMulticontextExample
#' mashInput(objects = list(strong = ss, random = ss))
#' @export
mashInput <- function(
    objects,
    nRandom = 10L,
    nNull = 10L,
    excludeCondition = character(0),
    coverage = 0.95,
    zOnly = FALSE,
    sigPCutoff = 1e-6,
    inputScale = c("auto", "beta", "z"),
    independentVariants = NULL,
    seed = 999L
) {
    inputScale <- arg_match(inputScale)
    if (!is.null(independentVariants)) {
        independentVariants <- as.character(independentVariants)
    }
    objects <- .mashPrepObjects(objects)
    cfg <- list(
        nRandom = nRandom,
        nNull = nNull,
        excludeCondition = excludeCondition,
        coverage = coverage,
        inputScale = inputScale,
        seed = seed,
        independentVariants = independentVariants
    )
    combined <- .mashCombinePartitions(objects, cfg)
    .mashFinalizeCombined(combined, sigPCutoff, zOnly)
}

# `objects` must be a non-empty list; unnamed lists get synthetic region names.
# @noRd
.mashPrepObjects <- function(objects) {
    if (!is.list(objects) || length(objects) == 0L) {
        msg <- glue(
            "mashInput: `objects` must be a non-empty list of QtlSumStats ",
            "and/or FineMappingResult objects."
        )
        abort(msg)
    }
    if (is.null(names(objects)) || any(str_length(names(objects)) == 0L)) {
        names(objects) <- str_c("region", seq_along(objects))
    }
    objects
}

# Extract + merge each object's strong/random/null partitions, disambiguating
# rownames by region before accumulating.
# @noRd
.mashCombinePartitions <- function(objects, cfg) {
    combined <- list()
    for (nm in names(objects)) {
        part <- .mashObjectPartitions(
            objects[[nm]],
            nRandom = cfg$nRandom,
            nNull = cfg$nNull,
            excludeCondition = cfg$excludeCondition,
            coverage = cfg$coverage,
            inputScale = cfg$inputScale,
            seed = cfg$seed,
            independentVariants = cfg$independentVariants
        )
        part <- map(part, .mashPrefixRownames, nm = nm)
        combined <- mergeMashData(combined, part)
    }
    combined
}

# Coerce to data.frame, clean each partition + derive z (random/null before
# strong so the strong significance filter runs once), restore the strong
# 1-row matrix shape, add the strong XtX, and optionally drop b/s slots.
# @noRd
.mashFinalizeCombined <- function(combined, sigPCutoff, zOnly) {
    combined <- map(combined, .mashAsDataFrameOrNull)
    for (cond in c("random", "null", "strong")) {
        bKey <- str_c(cond, ".b")
        sKey <- str_c(cond, ".s")
        if (!is.null(combined[[bKey]]) && !is.null(combined[[sKey]])) {
            combined <- filterInvalidSummaryStat(
                combined,
                bhat = bKey,
                sbhat = sKey,
                btoz = TRUE,
                sigPCutoff = sigPCutoff
            )
        }
    }
    combined <- .mashRestoreStrongShape(combined)
    combined <- .mashAddXtX(combined)
    if (zOnly) {
        combined[str_detect(names(combined), "\\.(b|s)$")] <- NULL
    }
    combined
}

# filterInvalidSummaryStat subsets strong without drop = FALSE, so a single
# surviving strong variant degrades to a vector; restore the 1-row matrix.
# @noRd
.mashRestoreStrongShape <- function(combined) {
    for (k in c("strong.b", "strong.s", "strong.z")) {
        v <- combined[[k]]
        if (!is.null(v) && is.null(dim(v))) {
            combined[[k]] <- matrix(
                v,
                nrow = 1L,
                dimnames = list(NULL, names(v))
            )
        }
    }
    combined
}

# Strong XtX cross-product (conditions x conditions), when strong.z is present.
# @noRd
.mashAddXtX <- function(combined) {
    if (
        !is.null(combined$strong.z) && nrow(as.matrix(combined$strong.z)) > 0L
    ) {
        sz <- as.matrix(combined$strong.z)
        combined$XtX <- crossprod(sz) / nrow(sz)
    }
    combined
}

#' @title Build a QtlSumStats from a Z-score matrix
#' @description Assemble a per-condition \code{\link{QtlSumStats}} from a
#'   \code{variants x conditions} Z-score matrix -- the input shape the mash
#'   pipeline uses when only Z is available. Each column is one condition, and
#'   conditions are distinguished by \code{context}, \code{trait}, or both: the
#'   columns may be different cell types / tissues (contexts), different
#'   molecular phenotypes (traits), or arbitrary context x trait pairs.
#'   Chromosome / position are decoded from the row (variant) identifiers via
#'   \code{\link{parseVariantId}} (with a synthetic-position fallback for ids
#'   that do not encode coordinates); \code{A1} / \code{A2} / \code{N} are
#'   placeholders because a Z-only input carries no alleles or sample sizes
#'   (mash reads only Z). A pass-through \code{qcInfo} record is set so the
#'   result clears the mash QC gate.
#' @param z Numeric matrix (variants x conditions). \code{rownames(z)} are
#'   variant ids (ideally \code{chr:pos:A2:A1}); \code{colnames(z)} label the
#'   conditions.
#' @param study Study identifier (recycled across conditions).
#' @param ldSketch A genotype panel (see \code{\link{readGenotypes}})
#'   embedded in the collection, or \code{NULL} (default) -- mash operates
#'   across
#'   conditions per variant and needs no LD reference.
#' @param context Condition context label(s): a single value recycled across
#'   every column, or a length-\code{ncol(z)} vector (one per condition).
#'   Defaults to \code{colnames(z)} -- one context per column.
#' @param trait Condition trait label(s): a single value recycled across every
#'   column, or a length-\code{ncol(z)} vector. Default \code{"mash"}. Pass
#'   \code{colnames(z)} here (with a constant \code{context}) when the columns
#'   are traits rather than contexts.
#' @param genome Genome build. Default \code{"GRCh38"}.
#' @param n Placeholder per-variant sample size. Default \code{1000}.
#' @param a1,a2 Placeholder alleles. Defaults \code{"A"} / \code{"G"}.
#' @param role Tag stored in the \code{qcInfo} record. Default \code{"mash"}.
#' @return A \code{\link{QtlSumStats}} with one entry per condition (column).
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @examples
#' panel <- readGenotypes(
#'   system.file("extdata", "toy_ref.bed", package = "pecotmr"))
#' z <- matrix(rnorm(6), 2, 3, dimnames = list(
#'   c("chr22:1:A:G", "chr22:2:A:G"), c("brain", "blood", "muscle")))
#' qtlSumStatsFromZMatrix(z = z, study = "s1", ldSketch = panel,
#'   context = colnames(z), trait = "g1", genome = "hg38", n = 100)
#' @export
qtlSumStatsFromZMatrix <- function(
    z,
    study,
    ldSketch = NULL,
    context = colnames(z),
    trait = "mash",
    genome = "GRCh38",
    n = 1000L,
    a1 = "A",
    a2 = "G",
    role = "mash"
) {
    if (!is.matrix(z) || !is.numeric(z)) {
        msg <- glue(
            "qtlSumStatsFromZMatrix: `z` must be a numeric variants x ",
            "conditions matrix."
        )
        abort(msg)
    }
    vids <- rownames(z)
    if (is.null(vids)) {
        vids <- str_c("var", seq_len(nrow(z)))
    }
    .qtlSumStatsFromMatrix(
        vids = vids,
        nCond = ncol(z),
        study = study,
        ldSketch = ldSketch,
        context = context,
        trait = trait,
        genome = genome,
        role = role,
        mcolFn = .mashZMcolFn,
        a1 = a1,
        a2 = a2,
        z = z,
        n = n
    )
}

#' @title Build a QtlSumStats from Bhat / Shat (effect-size) Matrices
#' @description Assemble a per-condition \code{\link{QtlSumStats}} from an
#'   aligned pair of \code{variants x conditions} effect-size (\code{Bhat}) and
#'   standard-error (\code{Shat}) matrices -- the beta-scale (EE) counterpart of
#'   \code{\link{qtlSumStatsFromZMatrix}}. Each entry carries \code{BETA},
#'   \code{SE}, and (derived) \code{Z = BETA / SE} mcols, so the result feeds
#'   \code{\link{mashPipeline}} / \code{\link{mashModelFit}} on either scale
#'   (\code{inputScale = "beta"} or \code{"z"}). Chromosome / position are
#'   decoded from the row (variant) ids exactly as in
#'   \code{\link{qtlSumStatsFromZMatrix}}.
#' @param bhat Numeric matrix (variants x conditions) of effect sizes.
#'   \code{rownames(bhat)} are variant ids; \code{colnames(bhat)} label the
#'   conditions.
#' @param shat Numeric matrix of standard errors, aligned with \code{bhat}
#'   (identical dimensions and row/column order).
#' @param study Study identifier (recycled across conditions).
#' @param ldSketch A genotype panel (see \code{\link{readGenotypes}})
#'   embedded in the collection, or \code{NULL} (default) -- mash operates
#'   across
#'   conditions per variant and needs no LD reference.
#' @param context,trait Condition labels; see
#'   \code{\link{qtlSumStatsFromZMatrix}}. Defaults \code{context =
#'   colnames(bhat)}, \code{trait = "mash"}.
#' @param genome Genome build. Default \code{"GRCh38"}.
#' @param n Placeholder per-variant sample size. Default \code{1000}.
#' @param a1,a2 Placeholder alleles. Defaults \code{"A"} / \code{"G"}.
#' @param role Tag stored in the \code{qcInfo} record. Default \code{"mash"}.
#' @return A \code{\link{QtlSumStats}} with one entry per condition (column).
#' @seealso \code{\link{qtlSumStatsFromZMatrix}}, \code{\link{mashModelFit}}
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @examples
#' panel <- readGenotypes(
#'   system.file("extdata", "toy_ref.bed", package = "pecotmr"))
#' bhat <- matrix(rnorm(6), 2, 3, dimnames = list(
#'   c("chr22:1:A:G", "chr22:2:A:G"), c("brain", "blood", "muscle")))
#' shat <- matrix(0.1, 2, 3, dimnames = dimnames(bhat))
#' qtlSumStatsFromBetaMatrix(bhat = bhat, shat = shat, study = "s1",
#'   ldSketch = panel, context = colnames(bhat), trait = "g1",
#'     genome = "hg38", n = 100)
#' @export
qtlSumStatsFromBetaMatrix <- function(
    bhat,
    shat,
    study,
    ldSketch = NULL,
    context = colnames(bhat),
    trait = "mash",
    genome = "GRCh38",
    n = 1000L,
    a1 = "A",
    a2 = "G",
    role = "mash"
) {
    .mashValidateBetaMatrix(bhat, shat)
    vids <- rownames(bhat)
    if (is.null(vids)) {
        vids <- str_c("var", seq_len(nrow(bhat)))
    }
    .qtlSumStatsFromMatrix(
        vids = vids,
        nCond = ncol(bhat),
        study = study,
        ldSketch = ldSketch,
        context = context,
        trait = trait,
        genome = genome,
        role = role,
        mcolFn = .mashBetaMcolFn,
        a1 = a1,
        a2 = a2,
        bhat = bhat,
        shat = shat,
        n = n
    )
}

# `bhat` / `shat` must be numeric variants x conditions matrices of identical
# dimension.
# @noRd
.mashValidateBetaMatrix <- function(bhat, shat) {
    if (!is.matrix(bhat) || !is.numeric(bhat)) {
        msg <- glue(
            "qtlSumStatsFromBetaMatrix: `bhat` must be a numeric ",
            "variants x conditions matrix."
        )
        abort(msg)
    }
    if (!is.matrix(shat) || !is.numeric(shat)) {
        msg <- glue(
            "qtlSumStatsFromBetaMatrix: `shat` must be a numeric ",
            "variants x conditions matrix."
        )
        abort(msg)
    }
    if (!identical(dim(bhat), dim(shat))) {
        msg <- glue(
            "qtlSumStatsFromBetaMatrix: ",
            "`bhat` ({str_flatten(dim(bhat), 'x')}) and ",
            "`shat` ({str_flatten(dim(shat), 'x')}) ",
            "must have identical dimensions."
        )
        abort(msg)
    }
}

# Internal: shared assembly for the z / beta matrix constructors. Recycles the
# context / trait labels, decodes chrom/pos from the variant ids (synthesising
# where they don't parse), builds one GRanges entry per condition with mcols
# from `mcolFn(j)`, and wraps the entries as a QtlSumStats.
# @noRd
.qtlSumStatsFromMatrix <- function(
    vids,
    nCond,
    study,
    ldSketch,
    context,
    trait,
    genome,
    role,
    mcolFn,
    ...
) {
    context <- .qszmRecycle(context, nCond, "context")
    trait <- .qszmRecycle(trait, nCond, "trait")
    # Decode chrom/pos from the variant ids; synthesise where they do not parse.
    parsed <- tryCatch(
        suppressWarnings(parseVariantId(vids)),
        error = function(e) NULL
    )
    chrom <- if (!is.null(parsed)) {
        as.character(parsed$chrom)
    } else {
        rep(NA_character_, length(vids))
    }
    pos <- if (!is.null(parsed)) {
        suppressWarnings(as.integer(parsed$pos))
    } else {
        rep(NA_integer_, length(vids))
    }
    chrom[is.na(chrom) | str_length(chrom) == 0L] <- "chr1"
    pos[is.na(pos)] <- seq_along(pos)[is.na(pos)]
    entries <- map(
        seq_len(nCond),
        .qszmEntry,
        chrom = chrom,
        pos = pos,
        vids = vids,
        mcolFn = mcolFn,
        mcolArgs = list(...)
    )
    QtlSumStats(
        study = rep(as.character(study), nCond),
        context = context,
        trait = trait,
        entry = entries,
        genome = genome,
        ldSketch = ldSketch,
        qcInfo = list(role = role, entryAudit = vector("list", nCond))
    )
}

# Internal: recycle a condition-label argument (context / trait) to one value
# per matrix column. Accepts a single value (recycled to every column) or a
# vector of length ncol; errors otherwise -- including NULL, which is how
# `context = colnames(x)` arrives when the matrix has no column names.
.qszmRecycle <- function(v, n, what) {
    if (is.null(v)) {
        msg <- glue(
            "qtlSumStats matrix constructor: `{what}` is NULL; pass a ",
            "length-1 or length-{n} value (or give the matrix column ",
            "names)."
        )
        abort(msg)
    }
    v <- as.character(v)
    if (length(v) == 1L) {
        return(rep(v, n))
    }
    if (length(v) == n) {
        return(v)
    }
    msg <- glue(
        "qtlSumStats matrix constructor: `{what}` must be length 1 or ",
        "ncol={n}, got {length(v)}."
    )
    abort(msg)
}

# Internal: convert a single SumStats object (post-QC) into a (Bhat, Shat)
# pair of matrices keyed by context. Each row of the matrix corresponds to
# one (variantId x (study, trait)) cell from the SumStats entries; each
# column corresponds to a context (from QtlSumStats $context; from
# GwasSumStats $study, which is the per-study mash column).
#
# For QtlSumStats:
#   * Pivots entries on (study, trait) so each (study, trait) becomes a
#     block of rows and each context becomes a column. Missing
#     (study, trait, context) cells are filled with NA.
# For GwasSumStats:
#   * Each row of the collection is one study; we treat each study as a
#     mash "context" (single block of rows per study, columns = studies).
#     This is rarely used on its own but lets a flat GwasSumStats pass
#     through alongside (or instead of) a QtlSumStats without special
#     casing further upstream.
#
# Variant alignment within a (study, trait) block uses the entry's
# variant order; missing variants in any one context are filled with NA.
# NA in Bhat is mapped to 0 and NA in Shat is mapped to a large value
# (1000) inside mashr::mash_set_data via its `zero_Bhat_Shat_reset`
# pathway, matching the prior pipeline's handling of incomplete cells.
# @noRd
.mashSumStatsToMatrices <- function(
    x,
    role,
    inputScale = c("auto", "beta", "z")
) {
    inputScale <- arg_match(inputScale)
    .mashValidateInput(x, role)
    setup <- .mashBlockSetup(x)
    resolvedScale <- .mashResolveScale(x, role, inputScale)
    blocks <- .mashBuildBlockMatrices(x, setup, resolvedScale)
    bhatBlocks <- blocks$bhat
    shatBlocks <- blocks$shat
    bhat <- exec(rbind, !!!bhatBlocks)
    shat <- exec(rbind, !!!shatBlocks)
    # bhat NA -> 0, shat NA / <= 0 -> 1000 (the mash_set_data
    # zero_Bhat_Shat_reset convention; missing-cell variants do not drive the
    # fit).
    bhat[is.na(bhat)] <- 0
    shat[is.na(shat) | shat <= 0] <- 1000
    list(b = bhat, s = shat)
}

# The SumStats input must be a QC'd, non-empty QtlSumStats / GwasSumStats.
# @noRd
.mashValidateInput <- function(x, role) {
    if (!methods::is(x, "QtlSumStats") && !methods::is(x, "GwasSumStats")) {
        msg <- glue(
            "mashPipeline: '{role}' input must be a QtlSumStats or ",
            "GwasSumStats; got {str_flatten(class(x), '/')}."
        )
        abort(msg)
    }
    if (length(getQcInfo(x)) == 0L) {
        msg <- glue(
            "mashPipeline: '{role}' SumStats has no QC info ",
            "(length(getQcInfo(x)) == 0L). ",
            "Run summaryStatsQc() on the SumStats before passing it to ",
            "mashPipeline()."
        )
        abort(msg)
    }
    if (nrow(x) == 0L) {
        msg <- glue(
            "mashPipeline: '{role}' SumStats has no entries (nrow == 0)."
        )
        abort(msg)
    }
}

# Block / column layout: QtlSumStats blocks by (study, trait) with context
# columns; GwasSumStats blocks by study with study columns.
# @noRd
.mashBlockSetup <- function(x) {
    isQtl <- methods::is(x, "QtlSumStats")
    studyCol <- as.character(x$study)
    if (isQtl) {
        traitCol <- as.character(x$trait)
        contextCol <- as.character(x$context)
        blockKeys <- str_c(studyCol, traitCol, sep = "::")
        columnLabels <- unique(contextCol)
    } else {
        traitCol <- NULL
        contextCol <- studyCol
        blockKeys <- studyCol
        columnLabels <- unique(studyCol)
    }
    list(
        isQtl = isQtl,
        studyCol = studyCol,
        traitCol = traitCol,
        contextCol = contextCol,
        blockKeys = blockKeys,
        columnLabels = columnLabels
    )
}

# Resolve which (Bhat, Shat) source to pull: "beta" (BETA/SE) or "z"
# (Z, Shat=1).
# "auto" picks beta when every entry has BETA+SE, else z; mixed inputs error.
# @noRd
.mashResolveScale <- function(x, role, inputScale) {
    entries <- .collectionEntries(x)
    caps <- map(entries, .mashEntryCaps)
    allHaveBetaSe <- all(map_lgl(caps, "hasBetaSe"))
    allHaveZ <- all(map_lgl(caps, "hasZ"))
    switch(
        inputScale,
        beta = .mashScaleBeta(allHaveBetaSe, role),
        z = .mashScaleZ(allHaveZ, role),
        auto = .mashScaleAuto(allHaveBetaSe, allHaveZ, role)
    )
}

# @noRd
.mashScaleBeta <- function(allHaveBetaSe, role) {
    if (!allHaveBetaSe) {
        msg <- glue(
            "mashPipeline: inputScale = 'beta' requires every '{role}' ",
            "entry to carry both BETA and SE mcols."
        )
        abort(msg)
    }
    "beta"
}

# @noRd
.mashScaleZ <- function(allHaveZ, role) {
    if (!allHaveZ) {
        msg <- glue(
            "mashPipeline: inputScale = 'z' requires every '{role}' entry ",
            "to carry a Z mcol."
        )
        abort(msg)
    }
    "z"
}

# @noRd
.mashScaleAuto <- function(allHaveBetaSe, allHaveZ, role) {
    if (allHaveBetaSe) {
        return("beta")
    }
    if (allHaveZ) {
        return("z")
    }
    msg <- glue(
        "mashPipeline: '{role}' SumStats has no usable scale - every ",
        "entry must carry (BETA, SE) or Z mcols."
    )
    abort(msg)
}

# Per (study, trait) block, a variant x context Bhat / Shat matrix pair.
# @noRd
.mashBuildBlockMatrices <- function(x, setup, resolvedScale) {
    blocks <- map(
        unique(setup$blockKeys),
        .mashBlockMatrix,
        x = x,
        setup = setup,
        resolvedScale = resolvedScale
    )
    list(bhat = map(blocks, "b"), shat = map(blocks, "s"))
}

# The per-context (Bhat, Shat) vectors for one block: the variant universe is
# the first-seen union of SNP ids across the block's contexts.
# @noRd
.mashBlockPerContext <- function(x, rowsInBlock, setup, resolvedScale) {
    requireCols <- if (resolvedScale == "beta") {
        c("SNP", "BETA", "SE")
    } else {
        c("SNP", "Z")
    }
    variantOrder <- character()
    perContextB <- list()
    perContextSe <- list()
    for (rIdx in rowsInBlock) {
        df <- .mashRowDf(x, rIdx, setup, requireCols)
        snps <- df$variant_id
        variantOrder <- c(variantOrder, setdiff(snps, variantOrder))
        ctx <- setup$contextCol[[rIdx]]
        if (resolvedScale == "beta") {
            perContextB[[ctx]] <- set_names(df$beta, snps)
            perContextSe[[ctx]] <- set_names(df$se, snps)
        } else {
            perContextB[[ctx]] <- set_names(df$z, snps)
            perContextSe[[ctx]] <- set_names(rep(1, length(snps)), snps)
        }
    }
    list(
        variantOrder = variantOrder,
        perContextB = perContextB,
        perContextSe = perContextSe
    )
}

# One row's sumstat data.frame (QtlSumStats keyed by study/context/trait;
# GwasSumStats by study).
# @noRd
.mashRowDf <- function(x, rIdx, setup, requireCols) {
    if (setup$isQtl) {
        getSumstatDf(
            x,
            study = setup$studyCol[[rIdx]],
            context = setup$contextCol[[rIdx]],
            trait = setup$traitCol[[rIdx]],
            require = requireCols
        )
    } else {
        getSumstatDf(x, study = setup$studyCol[[rIdx]], require = requireCols)
    }
}

# Assemble one block's (variant x context) Bhat / Shat matrices, disambiguating
# rownames by block key to avoid silent cross-block dedup.
# @noRd
.mashBlockMatrix <- function(bkey, x, setup, resolvedScale) {
    rowsInBlock <- which(setup$blockKeys == bkey)
    pc <- .mashBlockPerContext(x, rowsInBlock, setup, resolvedScale)
    dims <- list(pc$variantOrder, setup$columnLabels)
    nVar <- length(pc$variantOrder)
    nCol <- length(setup$columnLabels)
    bMat <- matrix(NA_real_, nrow = nVar, ncol = nCol, dimnames = dims)
    sMat <- matrix(NA_real_, nrow = nVar, ncol = nCol, dimnames = dims)
    for (ctx in names(pc$perContextB)) {
        bMat[names(pc$perContextB[[ctx]]), ctx] <- pc$perContextB[[ctx]]
        sMat[names(pc$perContextSe[[ctx]]), ctx] <- pc$perContextSe[[ctx]]
    }
    rn <- str_c(bkey, pc$variantOrder, sep = "::")
    rownames(bMat) <- rn
    rownames(sMat) <- rn
    list(b = bMat, s = sMat)
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# TRUE when any |z| in a matrix row reaches the significance threshold.
# @noRd
.mashRowExceeds <- function(row, zThreshold) {
    any(abs(row) >= zThreshold)
}

# The fraction of non-zero entries in a matrix row.
# @noRd
.mashRowNonzeroRate <- function(row) {
    mean(row != 0)
}

# TRUE when a covariance matrix is not identically zero.
# @noRd
.mashMatrixNonzero <- function(mat) {
    !all(mat == 0)
}

# Subset one U matrix to the kept conditions (erroring if any is absent).
# @noRd
.mashSubsetMatrix <- function(mat, conditionsToKeep) {
    missingConditions <- setdiff(conditionsToKeep, colnames(mat))
    if (length(missingConditions) > 0) {
        msg <- glue(
            "Condition(s) {str_flatten(missingConditions, ', ')} ",
            "not found in matrix"
        )
        abort(msg)
    }
    mat[conditionsToKeep, conditionsToKeep]
}

# Merge partition `d` of two mash data lists: a column-aligned row-bind (the two
# objects may measure different condition sets). bind_rows unions the columns,
# filling gaps with NA -> NaN. Returns a base data.frame because this backs the
# exported mergeMashData(), whose result is column-accessed (`$cond`); the
# mashInput pipeline then coerces these frames back to matrices. The variant-id
# rownames are load-bearing -- they survive (via as.matrix) as the output
# matrices' dimnames (tested, e.g. rownames(mashInput(...)$strong.z)), which a
# tibble (no rownames) would drop. Rows are APPENDED (each object's variants are
# distinct, disambiguated by the region prefix), so the row keys must be unique
# across the two sides -- a collision means the prefix invariant broke, and we
# error loudly rather than silently stack.
# @noRd
.mashCombineDatum <- function(d, oneData, resData) {
    od <- oneData[[d]]
    rd <- resData[[d]]
    if (length(od) == 0 || is.null(od)) {
        return(rd)
    }
    if (is.null(rd) || length(rd) == 0) {
        return(od)
    }
    rnRes <- rownames(as.data.frame(rd))
    rnOne <- rownames(as.data.frame(od))
    if (anyDuplicated(c(rnRes, rnOne)) > 0L) {
        abort(glue(
            "mergeMashData: duplicate variant ids across the merged ",
            "partitions -- each object's variants must be uniquely keyed ",
            "(the mashInput region prefix guarantees this). A collision ",
            "means two objects share a name or the prefix invariant broke."
        ))
    }
    combined <- bind_rows(as.data.frame(rd), as.data.frame(od))
    combined[is.na(combined)] <- NaN
    rownames(combined) <- c(rnRes, rnOne)
    combined
}

# Region-prefix one partition matrix's rownames (no-op for empty/NULL).
# @noRd
.mashPrefixRownames <- function(m, nm) {
    if (!is.null(m) && nrow(m) > 0L) {
        rownames(m) <- str_c(rownames(m), nm, sep = "_")
    }
    m
}

# Coerce one partition to a data.frame (NULL passes through). Kept as a base
# data.frame (not a tibble) so the variant-id rownames survive to the output
# matrices -- see .mashCombineDatum.
# @noRd
.mashAsDataFrameOrNull <- function(m) {
    if (is.null(m)) {
        return(NULL)
    } # nocov  (partitions are always matrices here, never NULL)
    as.data.frame(m)
}

# One condition's GRanges entry, mcols from `mcolFn(j, vids, <mcolArgs>)`.
# @noRd
.qszmEntry <- function(j, chrom, pos, vids, mcolFn, mcolArgs) {
    gr <- GenomicRanges::GRanges(
        seqnames = chrom,
        ranges = IRanges::IRanges(start = pos, width = 1L)
    )
    mcolCallArgs <- c(list(j, vids), mcolArgs)
    S4Vectors::mcols(gr) <- exec(mcolFn, !!!mcolCallArgs)
    gr
}

# mcols for condition `j` of a z-scale matrix (Z + placeholder N/alleles).
# @noRd
.mashZMcolFn <- function(j, vids, a1, a2, z, n) {
    S4Vectors::DataFrame(
        SNP = vids,
        A1 = rep(a1, length(vids)),
        A2 = rep(a2, length(vids)),
        Z = as.numeric(z[, j]),
        N = rep(as.integer(n), length(vids))
    )
}

# mcols for condition `j` of a beta-scale pair (BETA/SE + derived Z).
# @noRd
.mashBetaMcolFn <- function(j, vids, a1, a2, bhat, shat, n) {
    S4Vectors::DataFrame(
        SNP = vids,
        A1 = rep(a1, length(vids)),
        A2 = rep(a2, length(vids)),
        BETA = as.numeric(bhat[, j]),
        SE = as.numeric(shat[, j]),
        Z = as.numeric(bhat[, j] / shat[, j]),
        N = rep(as.integer(n), length(vids))
    )
}

# The (hasBetaSe, hasZ) scale capabilities of one sumstats entry.
# @noRd
.mashEntryCaps <- function(e) {
    mc <- S4Vectors::mcols(e)
    list(
        hasBetaSe = all(is_in(c("BETA", "SE"), colnames(mc))),
        hasZ = is_in("Z", colnames(mc))
    )
}
